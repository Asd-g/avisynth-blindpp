#include <algorithm>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <cstring>
#include <malloc.h>

#include "PostProcess.h"

// #define DEBUGMODE
// #define SELFCHECK
// #define PREFETCH

#ifdef PREFETCH
#define PREFETCH_AHEAD_V 8
#define PREFETCH_AHEAD_H 8
#define PREFETCH_ENABLE
#endif

#ifdef SELFCHECK
#define PP_SELF_CHECK
#define SELF_CHECK
#endif

#ifdef DEBUGMODE
#define SHOWDECISIONS_H
#define SHOW_DECISIONS
#define SHOWDECISIONS_V
#endif

void __stdcall fast_copy(
    const uint8_t* src, const int src_stride, uint8_t* dst, const int dst_stride, const int horizontal_size, int vertical_size) noexcept
{
    if (vertical_size == 0)
    {
        return;
    }
    else if (horizontal_size == src_stride && src_stride == dst_stride)
    {
        memcpy(dst, src, horizontal_size * vertical_size);
    }
    else
    {
        do
        {
            memcpy(dst, src, horizontal_size);
            dst += dst_stride;
            src += src_stride;
        } while (--vertical_size != 0);
    }
}

/* entry point for postprocessing */
void postprocess(unsigned char* src[], int src_stride, int UVsrc_stride, unsigned char* dst[], int dst_stride, int UVdst_stride,
    int horizontal_size, int vertical_size, QP_STORE_T* QP_store, int QP_stride, int mode, int moderate_h, int moderate_v, bool is422,
    bool iPP)
{

    uint8_t* puc_src;
    uint8_t* puc_dst;
    uint8_t* puc_flt;
    QP_STORE_T* QP_ptr;
    int y, i;

    if (iPP)
        vertical_size >>= 1; // field based post-processing

    /* this loop is (hopefully) going to improve performance */
    /* loop down the picture, copying and processing in vertical stripes, each four pixels high */
    for (y = 0; y < vertical_size; y += 4)
    {

        if (!(mode & PP_DONT_COPY))
        {
            if (!iPP)
            {
                puc_src = &((src[0])[y * src_stride]);
                puc_dst = &((dst[0])[y * dst_stride]);
                fast_copy(puc_src, src_stride, puc_dst, dst_stride, horizontal_size, 4);
            }
            else
            {
                puc_src = &((src[0])[y * 2 * src_stride]);
                puc_dst = &((dst[0])[y * 2 * dst_stride]);
                fast_copy(puc_src, src_stride, puc_dst, dst_stride, horizontal_size, 8);
            }
        }

        if (mode & PP_DEBLOCK_Y_H)
        {
            if (!iPP)
            {
                puc_flt = &((dst[0])[y * dst_stride]);
                QP_ptr = &(QP_store[(y >> 4) * QP_stride]);
                deblock_horiz(puc_flt, horizontal_size, dst_stride, QP_ptr, QP_stride, 0, moderate_h);
            }
            else
            {
                // top field
                puc_flt = &((dst[0])[y * 2 * dst_stride]);
                QP_ptr = &(QP_store[(y >> 4) * 2 * QP_stride]);
                deblock_horiz(puc_flt, horizontal_size, dst_stride * 2, QP_ptr, QP_stride * 2, 0, moderate_h);
                // bottom field
                puc_flt = &((dst[0])[(y * 2 + 1) * dst_stride]);
                QP_ptr = &(QP_store[((y >> 4) * 2 + 1) * QP_stride]);
                deblock_horiz(puc_flt, horizontal_size, dst_stride * 2, QP_ptr, QP_stride * 2, 0, moderate_h);
            }
        }

        if (mode & PP_DEBLOCK_Y_V)
        {
            if ((y % 8) && (y - 4) > 5)
            {
                if (!iPP)
                {
                    puc_flt = &((dst[0])[(y - 4) * dst_stride]);
                    QP_ptr = &(QP_store[(y >> 4) * QP_stride]);
                    deblock_vert(puc_flt, horizontal_size, dst_stride, QP_ptr, QP_stride, 0, moderate_v);
                }
                else
                {
                    // top field
                    puc_flt = &((dst[0])[(y - 4) * 2 * dst_stride]);
                    QP_ptr = &(QP_store[(y >> 4) * 2 * QP_stride]);
                    deblock_vert(puc_flt, horizontal_size, dst_stride * 2, QP_ptr, QP_stride * 2, 0, moderate_v);
                    // bottom field
                    puc_flt = &((dst[0])[((y - 4) * 2 + 1) * dst_stride]);
                    QP_ptr = &(QP_store[((y >> 4) * 2 + 1) * QP_stride]);
                    deblock_vert(puc_flt, horizontal_size, dst_stride * 2, QP_ptr, QP_stride * 2, 0, moderate_v);
                }
            }
        }
    }

    if (mode & PP_DERING_Y)
    {
        if (!iPP)
            dering(dst[0], horizontal_size, vertical_size, dst_stride, QP_store, QP_stride, 0);
        else
        {
            dering(dst[0], horizontal_size, vertical_size, dst_stride * 2, QP_store, QP_stride * 2, 0);
            dering(dst[0] + dst_stride, horizontal_size, vertical_size, dst_stride * 2, QP_store + QP_stride, QP_stride * 2, 0);
        }
    }

    // adjust for chroma
    if (!is422)
        vertical_size >>= 1;
    horizontal_size >>= 1;
    src_stride = UVsrc_stride;
    dst_stride = UVdst_stride;

    /* loop U then V */
    for (i = 1; i <= 2; i++)
    {
        for (y = 0; y < vertical_size; y += 4)
        {

            if (!(mode & PP_DONT_COPY))
            {
                if (!iPP)
                {
                    puc_src = &((src[i])[y * src_stride]);
                    puc_dst = &((dst[i])[y * dst_stride]);
                    fast_copy(puc_src, src_stride, puc_dst, dst_stride, horizontal_size, 4);
                }
                else
                {
                    puc_src = &((src[i])[y * 2 * src_stride]);
                    puc_dst = &((dst[i])[y * 2 * dst_stride]);
                    fast_copy(puc_src, src_stride, puc_dst, dst_stride, horizontal_size, 8);
                }
            }

            if (mode & PP_DEBLOCK_C_H)
            {
                if (!iPP)
                {
                    puc_flt = &((dst[i])[y * dst_stride]);
                    if (is422)
                        QP_ptr = &(QP_store[(y >> 4) * QP_stride]);
                    else
                        QP_ptr = &(QP_store[(y >> 3) * QP_stride]);
                    deblock_horiz(puc_flt, horizontal_size, dst_stride, QP_ptr, QP_stride, is422 ? 2 : 1, moderate_h);
                }
                else
                {
                    // top field
                    puc_flt = &((dst[i])[y * 2 * dst_stride]);
                    if (is422)
                        QP_ptr = &(QP_store[(y >> 4) * 2 * QP_stride]);
                    else
                        QP_ptr = &(QP_store[(y >> 3) * 2 * QP_stride]);
                    deblock_horiz(puc_flt, horizontal_size, dst_stride * 2, QP_ptr, QP_stride * 2, is422 ? 2 : 1, moderate_h);
                    // bottom field
                    puc_flt = &((dst[i])[(y * 2 + 1) * dst_stride]);
                    if (is422)
                        QP_ptr = &(QP_store[((y >> 4) * 2 + 1) * QP_stride]);
                    else
                        QP_ptr = &(QP_store[((y >> 3) * 2 + 1) * QP_stride]);
                    deblock_horiz(puc_flt, horizontal_size, dst_stride * 2, QP_ptr, QP_stride * 2, is422 ? 2 : 1, moderate_h);
                }
            }

            if (mode & PP_DEBLOCK_C_V)
            {
                if ((y % 8) && (y - 4) > 5)
                {
                    if (!iPP)
                    {
                        puc_flt = &((dst[i])[(y - 4) * dst_stride]);
                        if (is422)
                            QP_ptr = &(QP_store[(y >> 4) * QP_stride]);
                        else
                            QP_ptr = &(QP_store[(y >> 3) * QP_stride]);
                        deblock_vert(puc_flt, horizontal_size, dst_stride, QP_ptr, QP_stride, is422 ? 2 : 1, moderate_v);
                    }
                    else
                    {
                        // top field
                        puc_flt = &((dst[i])[(y - 4) * 2 * dst_stride]);
                        if (is422)
                            QP_ptr = &(QP_store[(y >> 4) * 2 * QP_stride]);
                        else
                            QP_ptr = &(QP_store[(y >> 3) * 2 * QP_stride]);
                        deblock_vert(puc_flt, horizontal_size, dst_stride * 2, QP_ptr, QP_stride * 2, is422 ? 2 : 1, moderate_v);
                        // bottom field
                        puc_flt = &((dst[i])[((y - 4) * 2 + 1) * dst_stride]);
                        if (is422)
                            QP_ptr = &(QP_store[((y >> 4) * 2 + 1) * QP_stride]);
                        else
                            QP_ptr = &(QP_store[((y >> 3) * 2 + 1) * QP_stride]);
                        deblock_vert(puc_flt, horizontal_size, dst_stride * 2, QP_ptr, QP_stride * 2, is422 ? 2 : 1, moderate_v);
                    }
                }
            }
        }
    }

    if (mode & PP_DERING_C)
    {
        if (!iPP)
        {
            dering(dst[1], horizontal_size, vertical_size, dst_stride, QP_store, QP_stride, is422 ? 2 : 1);
            dering(dst[2], horizontal_size, vertical_size, dst_stride, QP_store, QP_stride, is422 ? 2 : 1);
        }
        else
        {
            dering(dst[1], horizontal_size, vertical_size, dst_stride * 2, QP_store, QP_stride * 2, is422 ? 2 : 1);
            dering(dst[1] + dst_stride, horizontal_size, vertical_size, dst_stride * 2, QP_store + QP_stride, QP_stride * 2, is422 ? 2 : 1);
            dering(dst[2], horizontal_size, vertical_size, dst_stride * 2, QP_store, QP_stride * 2, is422 ? 2 : 1);
            dering(dst[2] + dst_stride, horizontal_size, vertical_size, dst_stride * 2, QP_store + QP_stride, QP_stride * 2, is422 ? 2 : 1);
        }
    }
}

/////////////////////////////////////////////////////////////////////////
// Post Processing Functions                                           //
//                                                                     //
/////////////////////////////////////////////////////////////////////////

/* this is a horizontal deblocking filter - i.e. it will smooth _vertical_ block edges */
void deblock_horiz(uint8_t* image, int width, int stride, QP_STORE_T* QP_store, int QP_stride, int chromaFlag, int moderate_h)
{
    /* loop over every block boundary in that row */
    for (int x = 8; x < width; x += 8)
    {
        int QP;
        /* extract QP from the decoder's array of QP values */
        if (chromaFlag == 0)
        {
            QP = QP_store[x / 16];
        }
        else
        {
            QP = QP_store[x / 8];
        }

        /* v points to pixel v0, in the left-hand block */
        uint8_t* v = image + x - 5;

        /* first decide whether to use default or DC offet mode */
        if (deblock_horiz_useDC(v, stride, moderate_h))
        { /* use DC offset mode */
            if (deblock_horiz_DC_on(v, stride, QP))
            {
                deblock_horiz_lpf9(v, stride, QP);
            }
        }
        else
        { /* use default mode */
            deblock_horiz_default_filter(v, stride, QP);
        }
    }
}

/* decide whether the DC filter should be turned on accoding to QP */
inline int deblock_horiz_DC_on(uint8_t* v, int stride, int QP)
{
    /* 99% of the time, this test turns out the same as the |max-min| strategy in the standard */
    QP *= 2;
    for (int i = 0; i < 4; ++i)
    {
        if (std::abs(v[0] - v[5]) >= QP)
            return false;
        if (std::abs(v[1] - v[8]) >= QP)
            return false;
        if (std::abs(v[1] - v[4]) >= QP)
            return false;
        if (std::abs(v[2] - v[7]) >= QP)
            return false;
        if (std::abs(v[3] - v[6]) >= QP)
            return false;
        v += stride;
    }
    return true;
}

/* horizontal deblocking filter used in default (non-DC) mode */
inline void deblock_horiz_default_filter(uint8_t* v, int stride, int QP)
{
    for (int y = 0; y < 4; ++y)
    {
        int q1 = v[4] - v[5];
        int q = q1 / 2;
        if (q != 0)
        {
            int a3_0 = 2 * (v[3] - v[6]) - 5 * q1;
            int aa3_0 = std::abs(a3_0);

            if (aa3_0 < 8 * QP)
            {
                int a3_1 = std::abs(5 * (v[3] - v[2]) + 2 * (v[1] - v[4]));
                int a3_2 = std::abs(5 * (v[7] - v[8]) + 2 * (v[5] - v[8]));

                int d = aa3_0 - std::min(a3_1, a3_2);

                if (d > 6)
                {
                    int delta = std::min(d, std::abs(q));

                    if ((q > 0 && a3_0 < 0) || (q < 0 && a3_0 > 0))
                    {
                        if (q > 0)
                        {
                            v[4] -= delta;
                            v[5] += delta;
                        }
                        else
                        {
                            v[4] += delta;
                            v[5] -= delta;
                        }
                    }
                }
            }
        }
        v += stride;
    }
}

/* The 9-tap low pass filter used in "DC" regions */
inline void deblock_horiz_lpf9(uint8_t* v, int stride, int QP)
{
    static const int C[8][10] = {
        // p1, v1, v2, v3, v4, v5, v6, v7, v8, p2
        {6, 4, 2, 2, 1, 1, 0, 0, 0, 0}, // v1' (sum=16)
        {4, 2, 4, 2, 2, 1, 1, 0, 0, 0}, // v2' (sum=16)
        {2, 2, 2, 4, 2, 2, 1, 1, 0, 0}, // v3' (sum=16)
        {1, 1, 2, 2, 4, 2, 2, 1, 1, 0}, // v4' (sum=16)
        {0, 1, 1, 2, 2, 4, 2, 2, 1, 1}, // v5' (sum=16)
        {0, 0, 1, 1, 2, 2, 4, 2, 2, 2}, // v6' (sum=16)
        {0, 0, 0, 1, 1, 2, 2, 4, 2, 4}, // v7' (sum=16)
        {0, 0, 0, 0, 1, 1, 2, 2, 4, 6}, // v8' (sum=16)
    };

    for (int y = 0; y < 4; ++y)
    {
        int p1 = (std::abs(v[0] - v[1]) < QP) ? v[0] : v[1];
        int p2 = (std::abs(v[8] - v[9]) < QP) ? v[9] : v[8];

        const int in[10] = {p1, v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8], p2};
        uint8_t out[8];

        for (int i = 0; i < 8; ++i)
        {
            int sum = 8; // Add 8 for rounding: (sum + 8) >> 4
            for (int j = 0; j < 10; ++j)
            {
                sum += in[j] * C[i][j];
            }
            out[i] = std::max(0, std::min(255, sum >> 4));
        }

        for (int i = 0; i < 8; ++i)
        {
            v[i + 1] = out[i];
        }

        v += stride;
    }
}

/* decide DC mode or default mode for the horizontal filter */
inline int deblock_horiz_useDC(uint8_t* v, int stride, int moderate_h)
{
    int eq_cnt = 0;
    for (int y = 0; y < 4; ++y)
    {
        const uint8_t* p = v + y * stride;
        for (int x = 1; x < 8; ++x)
        {
            if (std::abs((int)p[x] - (int)p[x + 1]) <= 1)
            {
                eq_cnt++;
            }
        }
    }
    return eq_cnt >= moderate_h;
}

/* this is a vertical deblocking filter - i.e. it will smooth _horizontal_ block edges */
void deblock_vert(uint8_t* image, int width, int stride, QP_STORE_T* QP_store, int QP_stride, int chromaFlag, int moderate_v)
{
    uint64_t v_local[20];
    uint64_t p1p2[4];
    int Bx, x, y;
    int QP;
    uint8_t* v;
    int useDC, DC_on;

    y = 0;

    /* loop over all blocks, left to right */
    for (Bx = 0; Bx < width; Bx += 8)
    {
        QP = chromaFlag == 1   ? QP_store[y / 8 * QP_stride + Bx / 8]
             : chromaFlag == 0 ? QP_store[y / 16 * QP_stride + Bx / 16]
                               : QP_store[y / 16 * QP_stride + Bx / 8];
        v = &(image[y * stride + Bx]) - 5 * stride;

        /* decide whether to use DC mode on a block-by-block basis */
        useDC = deblock_vert_useDC(v, stride, moderate_v);

        if (useDC)
        {
            /* we are in DC mode for this block.  But we only want to filter low-energy areas */
            DC_on = deblock_vert_DC_on(v, stride, QP);

            if (DC_on)
            {
                v = &(image[y * stride + Bx]) - 5 * stride;
                deblock_vert_copy_and_unpack(stride, &(v[stride]), &v_local[2], 8);
                deblock_vert_choose_p1p2(v, stride, p1p2, QP);
                deblock_vert_lpf9(v_local, p1p2, v, stride);
            }
        }
        else /* use the default filter */
        {
            x = Bx;
            v = &(image[y * stride + x]) - 5 * stride;
            deblock_vert_default_filter(v, stride, QP);
        }
    }
}

/* This function chooses the "endstops" for the vertial LPF9 filter: p1 and p2 */
/* We also convert these to 16-bit values here */
inline void deblock_vert_choose_p1p2(uint8_t* v, int stride, uint64_t* p1p2, int QP)
{
    uint16_t* p1_ptr = reinterpret_cast<uint16_t*>(p1p2);
    uint16_t* p2_ptr = reinterpret_cast<uint16_t*>(p1p2 + 2);

    const uint8_t* v0 = v;
    const uint8_t* v1 = v + stride;
    const uint8_t* v8 = v + 8 * stride;
    const uint8_t* v9 = v + 9 * stride;

    for (int i = 0; i < 8; ++i)
    {
        if (std::abs((int)v0[i] - (int)v1[i]) <= QP)
        {
            p1_ptr[i] = v0[i];
        }
        else
        {
            p1_ptr[i] = v1[i];
        }
        if (std::abs((int)v8[i] - (int)v9[i]) <= QP)
        {
            p2_ptr[i] = v9[i];
        }
        else
        {
            p2_ptr[i] = v8[i];
        }
    }
}

/* function using MMX to copy an 8-pixel wide column and unpack to 16-bit values */
/* n is the number of rows to copy - this must be even */
inline void deblock_vert_copy_and_unpack(int stride, uint8_t* source, uint64_t* dest, int n)
{
    auto* d = reinterpret_cast<uint16_t*>(dest);
    for (int i = 0; i < n; ++i)
    {
        const uint8_t* s_row = source + i * stride;
        uint16_t* d_row = d + i * 8;
        for (int j = 0; j < 8; ++j)
        {
            d_row[j] = s_row[j];
        }
    }
}

/* decide whether the DC filter should be turned on accoding to QP */
inline int deblock_vert_DC_on(uint8_t* v, int stride, int QP)
{
    QP *= 2;
    for (int i = 0; i < 8; ++i)
    {
        if (std::abs(v[i + 0 * stride] - v[i + 5 * stride]) >= QP)
            return false;
        if (std::abs(v[i + 1 * stride] - v[i + 4 * stride]) >= QP)
            return false;
        if (std::abs(v[i + 1 * stride] - v[i + 8 * stride]) >= QP)
            return false;
        if (std::abs(v[i + 2 * stride] - v[i + 7 * stride]) >= QP)
            return false;
        if (std::abs(v[i + 3 * stride] - v[i + 6 * stride]) >= QP)
            return false;
    }
    return true;
}

/* Vertical deblocking filter for use in non-flat picture regions */
inline void deblock_vert_default_filter(uint8_t* v, int stride, int QP)
{
    for (int i = 0; i < 8; ++i)
    {
        uint8_t* p = v + i;
        int v1 = p[1 * stride], v2 = p[2 * stride], v3 = p[3 * stride], v4 = p[4 * stride];
        int v5 = p[5 * stride], v6 = p[6 * stride], v7 = p[7 * stride], v8 = p[8 * stride];

        int q1 = v4 - v5;
        int q = q1 / 2;

        if (q != 0)
        {
            int a3_0 = 2 * (v3 - v6) - 5 * q1;
            int aa3_0 = std::abs(a3_0);

            if (aa3_0 < 8 * QP)
            {
                int a3_1 = std::abs(5 * (v3 - v2) + 2 * (v1 - v4));
                int a3_2 = std::abs(5 * (v7 - v6) + 2 * (v5 - v8));

                int d_abs = aa3_0 - std::min(a3_1, a3_2);

                if (d_abs > 0)
                {
                    int d = (5 * d_abs + 32) >> 6;
                    int final_delta = 0;

                    if ((q > 0 && a3_0 < 0) || (q < 0 && a3_0 >= 0))
                    {
                        int delta_mag = std::min(d, std::abs(q));
                        final_delta = (q > 0) ? delta_mag : -delta_mag;
                    }

                    int p4_new = (int)p[4 * stride] - final_delta;
                    int p5_new = (int)p[5 * stride] + final_delta;

                    p[4 * stride] = std::max(0, std::min(255, p4_new));
                    p[5 * stride] = std::max(0, std::min(255, p5_new));
                }
            }
        }
    }
}

/* Vertical 9-tap low-pass filter for use in "DC" regions of the picture */
inline void deblock_vert_lpf9(uint64_t* v_local, uint64_t* p1p2, uint8_t* v, int stride)
{
    auto* vv = reinterpret_cast<uint16_t*>(v_local) + 8;
    auto* p1p2_16 = reinterpret_cast<uint16_t*>(p1p2);

    for (int i = 0; i < 8; ++i)
    { // 8 columns
        int p1 = p1p2_16[i];
        int p2 = p1p2_16[i + 8];

        int v1 = vv[0 * 8 + i], v2 = vv[1 * 8 + i], v3 = vv[2 * 8 + i], v4 = vv[3 * 8 + i];
        int v5 = vv[4 * 8 + i], v6 = vv[5 * 8 + i], v7 = vv[6 * 8 + i], v8 = vv[7 * 8 + i];

        int psum = 4 + p1 * 3 + v1 + v2 + v3 + v4;
        int s;

        s = ((psum + v1) << 1) - (v4 - v5);
        v[(1 * stride) + i] = std::max(0, std::min(255, s >> 4));
        psum += v5 - p1;
        s = ((psum + v2) << 1) - (v5 - v6);
        v[(2 * stride) + i] = std::max(0, std::min(255, s >> 4));
        psum += v6 - p1;
        s = ((psum + v3) << 1) - (v6 - v7);
        v[(3 * stride) + i] = std::max(0, std::min(255, s >> 4));
        psum += v7 - p1;
        s = ((psum + v4) << 1) - (v7 - v8) + p1 - v1;
        v[(4 * stride) + i] = std::max(0, std::min(255, s >> 4));
        psum += v8 - v1;
        s = ((psum + v5) << 1) - (v8 - p2) + v1 - v2;
        v[(5 * stride) + i] = std::max(0, std::min(255, s >> 4));
        psum += p2 - v2;
        s = ((psum + v6) << 1) + v2 - v3;
        v[(6 * stride) + i] = std::max(0, std::min(255, s >> 4));
        psum += p2 - v3;
        s = ((psum + v7) << 1) + v3 - v4;
        v[(7 * stride) + i] = std::max(0, std::min(255, s >> 4));
        psum += p2 - v4;
        s = ((psum + v8) << 1) + v4 - v5;
        v[(8 * stride) + i] = std::max(0, std::min(255, s >> 4));
    }
}

/* decide DC mode or default mode in assembler */
inline int deblock_vert_useDC(uint8_t* v, int stride, int moderate_v)
{
    int eq_cnt = 0;
    for (int i = 1; i < 8; ++i)
    {
        const uint8_t* p1 = v + i * stride;
        const uint8_t* p2 = v + (i + 1) * stride;
        for (int j = 0; j < 8; ++j)
        {
            if (std::abs((int)p1[j] - (int)p2[j]) <= 1)
            {
                eq_cnt++;
            }
        }
    }
    return eq_cnt > moderate_v;
}

/* this is the deringing filter */
void dering(uint8_t* image, int width, int height, int stride, QP_STORE_T* QP_store, int QP_stride, int chroma)
{
    uint8_t b8x8filtered[64];

    /* loop over all the 8x8 blocks in the image... */
    /* don't process outer row of blocks for the time being. */
    for (int y = 8; y < height - 8; y += 8)
    {
        for (int x = 8; x < width - 8; x += 8)
        {
            /* QP for this block.. */
            QP_STORE_T QP = chroma == 1   ? QP_store[(y >> 3) * QP_stride + (x >> 3)]
                            : chroma == 0 ? QP_store[(y >> 4) * QP_stride + (x >> 4)]
                                          : QP_store[(y >> 4) * QP_stride + (x >> 3)];

            /* pointer to the top left pixel in 8x8   block */
            uint8_t* b8x8 = &(image[stride * y + x]);

            /* pointer to the top left pixel in 10x10 block */
            uint8_t* b10x10 = &(image[stride * (y - 1) + (x - 1)]);

            // Threshold determination - find min and max grey levels in the block
            uint8_t min_val = 255, max_val = 0;
            for (int r = 0; r < 8; ++r)
            {
                for (int c = 0; c < 8; ++c)
                {
                    uint8_t p = b8x8[r * stride + c];
                    min_val = std::min(min_val, p);
                    max_val = std::max(max_val, p);
                }
            }

            /* Threshold detirmination - compute threshold and dynamic range */
            uint8_t thr = (max_val + min_val + 1) >> 1;
            int max_diff = QP >> 1;

            // Filtering and clipping
            auto pavgb = [](int a, int b) { return (a + b + 1) >> 1; };
            for (int r = 0; r < 8; ++r)
            {
                for (int c = 0; c < 8; ++c)
                {
                    const uint8_t* p_row0 = &b10x10[(r + 0) * stride + c];
                    const uint8_t* p_row1 = &b10x10[(r + 1) * stride + c];
                    const uint8_t* p_row2 = &b10x10[(r + 2) * stride + c];

                    // Replicate the exact sequence of pavgb from assembly
                    int h_avg0 = pavgb(pavgb(p_row0[0], p_row0[2]), p_row0[1]);
                    int h_avg1 = pavgb(pavgb(p_row1[0], p_row1[2]), p_row1[1]);
                    int h_avg2 = pavgb(pavgb(p_row2[0], p_row2[2]), p_row2[1]);
                    int filtered_val = pavgb(pavgb(h_avg0, h_avg2), h_avg1);

                    uint8_t orig_val = b10x10[(r + 1) * stride + (c + 1)];
                    int lower_bound = std::max(0, (int)orig_val - max_diff);
                    int upper_bound = std::min(255, (int)orig_val + max_diff);

                    filtered_val = std::max(lower_bound, std::min(upper_bound, filtered_val));

                    b8x8filtered[r * 8 + c] = (uint8_t)filtered_val;
                }
            }

            // Decide whether to use filtered value
            for (int r = 0; r < 8; ++r)
            {
                for (int c = 0; c < 8; ++c)
                {
                    bool all_above = true;
                    bool all_below = true;
                    for (int dr = 0; dr < 3; ++dr)
                    {
                        for (int dc = 0; dc < 3; ++dc)
                        {
                            uint8_t p = b10x10[(r + dr) * stride + (c + dc)];
                            if (p < thr)
                                all_above = false;
                            if (p >= thr)
                                all_below = false;
                        }
                    }

                    if (all_above || all_below)
                    {
                        b8x8[r * stride + c] = b8x8filtered[r * 8 + c];
                    }
                }
            }
        }
    }
}
