#include <algorithm>
#include <cmath>
#include <cstring>

#include "PostProcess.h"

namespace
{
    void deblock_horiz(uint8_t* AVS_RESTRICT image, int width, int stride, const QP_STORE_T* AVS_RESTRICT QP_store, int QP_stride,
        int chroma_flag_for_qp, int moderate_h) noexcept;
    void deblock_vert(uint8_t* AVS_RESTRICT image, int width, int stride, const QP_STORE_T* AVS_RESTRICT QP_store, int QP_stride,
        int chroma_flag_for_qp, int moderate_v) noexcept;
    void dering(uint8_t* AVS_RESTRICT image, int width, int height, int stride, const QP_STORE_T* AVS_RESTRICT QP_store, int QP_stride,
        int chroma_flag_for_qp) noexcept;

    [[nodiscard]] inline int deblock_horiz_useDC(const uint8_t* AVS_RESTRICT v, int stride, int moderate_h) noexcept
    {
        int eq_cnt{};

        for (int y = 0; y < 4; ++y)
        {
            const uint8_t* p = v + y * stride;

            for (int x = 1; x < 8; ++x)
            {
                if (std::abs(static_cast<int>(p[x]) - static_cast<int>(p[x + 1])) <= 1)
                    eq_cnt++;
            }
        }

        return eq_cnt >= moderate_h;
    }

    [[nodiscard]] inline int deblock_horiz_DC_on(const uint8_t* AVS_RESTRICT v, int stride, int QP) noexcept
    {
        QP *= 2;

        for (int i = 0; i < 4; ++i)
        {
            if (std::abs(v[0] - v[5]) >= QP || std::abs(v[1] - v[8]) >= QP || std::abs(v[1] - v[4]) >= QP || std::abs(v[2] - v[7]) >= QP ||
                std::abs(v[3] - v[6]) >= QP)
                return false;

            v += stride;
        }

        return true;
    }

    inline void deblock_horiz_lpf9(uint8_t* AVS_RESTRICT v, int stride, int QP) noexcept
    {
        static constexpr int C[8][10]{
            {6, 4, 2, 2, 1, 1, 0, 0, 0, 0},
            {4, 2, 4, 2, 2, 1, 1, 0, 0, 0},
            {2, 2, 2, 4, 2, 2, 1, 1, 0, 0},
            {1, 1, 2, 2, 4, 2, 2, 1, 1, 0},
            {0, 1, 1, 2, 2, 4, 2, 2, 1, 1},
            {0, 0, 1, 1, 2, 2, 4, 2, 2, 2},
            {0, 0, 0, 1, 1, 2, 2, 4, 2, 4},
            {0, 0, 0, 0, 1, 1, 2, 2, 4, 6},
        };

        for (int y = 0; y < 4; ++y)
        {
            const int p1 = (std::abs(v[0] - v[1]) < QP) ? v[0] : v[1];
            const int p2 = (std::abs(v[8] - v[9]) < QP) ? v[9] : v[8];

            const int in[10]{p1, v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8], p2};
            uint8_t out[8];

            for (int i = 0; i < 8; ++i)
            {
                int sum{8};

                for (int j = 0; j < 10; ++j)
                    sum += in[j] * C[i][j];

                out[i] = static_cast<uint8_t>(std::clamp(sum >> 4, 0, 255));
            }

            for (int i = 0; i < 8; ++i)
                v[i + 1] = out[i];

            v += stride;
        }
    }

    inline void deblock_horiz_default_filter(uint8_t* AVS_RESTRICT v, int stride, int QP) noexcept
    {
        for (int y = 0; y < 4; ++y)
        {
            const int q1 = v[4] - v[5];
            const int q = q1 / 2;

            if (q != 0)
            {
                const int a3_0 = 2 * (v[3] - v[6]) - 5 * q1;
                const int aa3_0 = std::abs(a3_0);

                if (aa3_0 < 8 * QP)
                {
                    const int a3_1 = std::abs(5 * (v[3] - v[2]) + 2 * (v[1] - v[4]));
                    const int a3_2 = std::abs(5 * (v[7] - v[8]) + 2 * (v[5] - v[8]));
                    const int d = aa3_0 - std::min(a3_1, a3_2);

                    if (d > 6)
                    {
                        const int delta_mag = std::min(d, std::abs(q));
                        int delta = 0;

                        if ((q > 0 && a3_0 < 0) || (q < 0 && a3_0 > 0))
                            delta = (q > 0) ? -delta_mag : delta_mag;

                        if (delta != 0)
                        {
                            v[4] = static_cast<uint8_t>(std::clamp(static_cast<int>(v[4]) + delta, 0, 255));
                            v[5] = static_cast<uint8_t>(std::clamp(static_cast<int>(v[5]) - delta, 0, 255));
                        }
                    }
                }
            }

            v += stride;
        }
    }

    void deblock_horiz(uint8_t* AVS_RESTRICT image, int width, int stride, const QP_STORE_T* AVS_RESTRICT QP_store, int QP_stride,
        int chroma_flag_for_qp, int moderate_h) noexcept
    {
        for (int x = 8; x < width; x += 8)
        {
            const int QP = QP_store[(chroma_flag_for_qp == 0 || chroma_flag_for_qp == 3) ? (x / 16) : (x / 8)];
            uint8_t* v = image + x - 5;

            if (deblock_horiz_useDC(v, stride, moderate_h))
            {
                if (deblock_horiz_DC_on(v, stride, QP))
                    deblock_horiz_lpf9(v, stride, QP);
            }
            else
                deblock_horiz_default_filter(v, stride, QP);
        }
    }

    [[nodiscard]] inline int deblock_vert_useDC(const uint8_t* AVS_RESTRICT v, int stride, int moderate_v) noexcept
    {
        int eq_cnt{};

        for (int i = 1; i < 8; ++i)
        {
            const uint8_t* p1 = v + i * stride;
            const uint8_t* p2 = v + (i + 1) * stride;

            for (int j = 0; j < 8; ++j)
            {
                if (std::abs(static_cast<int>(p1[j]) - static_cast<int>(p2[j])) <= 1)
                    eq_cnt++;
            }
        }

        return eq_cnt > moderate_v;
    }

    [[nodiscard]] inline int deblock_vert_DC_on(const uint8_t* AVS_RESTRICT v, int stride, int QP) noexcept
    {
        QP *= 2;

        for (int i = 0; i < 8; ++i)
        {
            if (std::abs(v[i + 0 * stride] - v[i + 5 * stride]) >= QP || std::abs(v[i + 1 * stride] - v[i + 4 * stride]) >= QP ||
                std::abs(v[i + 1 * stride] - v[i + 8 * stride]) >= QP || std::abs(v[i + 2 * stride] - v[i + 7 * stride]) >= QP ||
                std::abs(v[i + 3 * stride] - v[i + 6 * stride]) >= QP)
                return false;
        }

        return true;
    }

    inline void deblock_vert_copy_and_unpack(int stride, const uint8_t* AVS_RESTRICT source, uint64_t* AVS_RESTRICT dest, int n) noexcept
    {
        auto* d = reinterpret_cast<uint16_t*>(dest);

        for (int i = 0; i < n; ++i)
        {
            const uint8_t* s_row = source + i * stride;
            uint16_t* d_row = d + i * 8;

            for (int j = 0; j < 8; ++j)
                d_row[j] = s_row[j];
        }
    }

    inline void deblock_vert_choose_p1p2(const uint8_t* AVS_RESTRICT v, int stride, uint64_t* AVS_RESTRICT p1p2, int QP) noexcept
    {
        auto* p1_ptr = reinterpret_cast<uint16_t*>(p1p2);
        auto* p2_ptr = reinterpret_cast<uint16_t*>(p1p2 + 2);

        const uint8_t* v0 = v;
        const uint8_t* v1 = v + stride;
        const uint8_t* v8 = v + 8 * stride;
        const uint8_t* v9 = v + 9 * stride;

        for (int i = 0; i < 8; ++i)
        {
            p1_ptr[i] = (std::abs(static_cast<int>(v0[i]) - static_cast<int>(v1[i])) <= QP) ? v0[i] : v1[i];
            p2_ptr[i] = (std::abs(static_cast<int>(v8[i]) - static_cast<int>(v9[i])) <= QP) ? v9[i] : v8[i];
        }
    }

    inline void deblock_vert_lpf9(
        const uint64_t* AVS_RESTRICT v_local, const uint64_t* AVS_RESTRICT p1p2, uint8_t* AVS_RESTRICT v, int stride) noexcept
    {
        const auto* vv = reinterpret_cast<const uint16_t*>(v_local) + 8;
        const auto* p1p2_16 = reinterpret_cast<const uint16_t*>(p1p2);

        for (int i = 0; i < 8; ++i)
        {
            const int p1 = p1p2_16[i];
            const int p2 = p1p2_16[i + 8];

            const int v1 = vv[0 * 8 + i], v2 = vv[1 * 8 + i], v3 = vv[2 * 8 + i], v4 = vv[3 * 8 + i];
            const int v5 = vv[4 * 8 + i], v6 = vv[5 * 8 + i], v7 = vv[6 * 8 + i], v8 = vv[7 * 8 + i];

            int psum = 4 + p1 * 3 + v1 + v2 + v3 + v4;

            int s = ((psum + v1) << 1) - (v4 - v5);
            v[(1 * stride) + i] = static_cast<uint8_t>(std::clamp(s >> 4, 0, 255));
            psum += v5 - p1;
            s = ((psum + v2) << 1) - (v5 - v6);
            v[(2 * stride) + i] = static_cast<uint8_t>(std::clamp(s >> 4, 0, 255));
            psum += v6 - p1;
            s = ((psum + v3) << 1) - (v6 - v7);
            v[(3 * stride) + i] = static_cast<uint8_t>(std::clamp(s >> 4, 0, 255));
            psum += v7 - p1;
            s = ((psum + v4) << 1) - (v7 - v8) + p1 - v1;
            v[(4 * stride) + i] = static_cast<uint8_t>(std::clamp(s >> 4, 0, 255));
            psum += v8 - v1;
            s = ((psum + v5) << 1) - (v8 - p2) + v1 - v2;
            v[(5 * stride) + i] = static_cast<uint8_t>(std::clamp(s >> 4, 0, 255));
            psum += p2 - v2;
            s = ((psum + v6) << 1) + v2 - v3;
            v[(6 * stride) + i] = static_cast<uint8_t>(std::clamp(s >> 4, 0, 255));
            psum += p2 - v3;
            s = ((psum + v7) << 1) + v3 - v4;
            v[(7 * stride) + i] = static_cast<uint8_t>(std::clamp(s >> 4, 0, 255));
            psum += p2 - v4;
            s = ((psum + v8) << 1) + v4 - v5;
            v[(8 * stride) + i] = static_cast<uint8_t>(std::clamp(s >> 4, 0, 255));
        }
    }

    inline void deblock_vert_default_filter(uint8_t* AVS_RESTRICT v, int stride, int QP) noexcept
    {
        for (int i = 0; i < 8; ++i)
        {
            uint8_t* p = v + i;
            const int v1 = p[1 * stride], v2 = p[2 * stride], v3 = p[3 * stride], v4 = p[4 * stride];
            const int v5 = p[5 * stride], v6 = p[6 * stride], v7 = p[7 * stride], v8 = p[8 * stride];

            const int q1 = v4 - v5;
            const int q = q1 / 2;

            if (q != 0)
            {
                const int a3_0 = 2 * (v3 - v6) - 5 * q1;
                const int aa3_0 = std::abs(a3_0);

                if (aa3_0 < 8 * QP)
                {
                    const int a3_1 = std::abs(5 * (v3 - v2) + 2 * (v1 - v4));
                    const int a3_2 = std::abs(5 * (v7 - v6) + 2 * (v5 - v8));
                    const int d_abs = aa3_0 - std::min(a3_1, a3_2);

                    if (d_abs > 0)
                    {
                        const int d = (5 * d_abs + 32) >> 6;
                        int final_delta{};

                        if ((q > 0 && a3_0 < 0) || (q < 0 && a3_0 >= 0))
                        {
                            const int delta_mag = std::min(d, std::abs(q));
                            final_delta = (q > 0) ? delta_mag : -delta_mag;
                        }

                        if (final_delta != 0)
                        {
                            p[4 * stride] = static_cast<uint8_t>(std::clamp(static_cast<int>(p[4 * stride]) - final_delta, 0, 255));
                            p[5 * stride] = static_cast<uint8_t>(std::clamp(static_cast<int>(p[5 * stride]) + final_delta, 0, 255));
                        }
                    }
                }
            }
        }
    }

    void deblock_vert(uint8_t* AVS_RESTRICT image, int width, int stride, const QP_STORE_T* AVS_RESTRICT QP_store, int QP_stride,
        int chroma_flag_for_qp, int moderate_v) noexcept
    {
        alignas(16) uint64_t v_local[20];
        alignas(16) uint64_t p1p2[4];

        for (int Bx = 0; Bx < width; Bx += 8)
        {
            const int QP = QP_store[(chroma_flag_for_qp == 0 || chroma_flag_for_qp == 3) ? (Bx / 16) : (Bx / 8)];
            uint8_t* v = image - 5 * stride + Bx;

            if (deblock_vert_useDC(v, stride, moderate_v))
            {
                if (deblock_vert_DC_on(v, stride, QP))
                {
                    deblock_vert_copy_and_unpack(stride, v + stride, &v_local[2], 8);
                    deblock_vert_choose_p1p2(v, stride, p1p2, QP);
                    deblock_vert_lpf9(v_local, p1p2, v, stride);
                }
            }
            else
                deblock_vert_default_filter(v, stride, QP);
        }
    }

    void dering(uint8_t* AVS_RESTRICT image, int width, int height, int stride, const QP_STORE_T* AVS_RESTRICT QP_store, int QP_stride,
        int chroma_flag_for_qp) noexcept
    {
        alignas(16) uint8_t b8x8filtered[64];

        for (int y = 8; y < height - 8; y += 8)
        {
            for (int x = 8; x < width - 8; x += 8)
            {
                int qp_x_idx;
                int qp_y_idx_shift;

                switch (chroma_flag_for_qp)
                {
                case 1: // YUV420
                    qp_x_idx = x >> 3;
                    qp_y_idx_shift = 3;
                    break;
                case 2: // YUV422
                    qp_x_idx = x >> 3;
                    qp_y_idx_shift = 4;
                    break;
                default: // Y8, YUV444
                    qp_x_idx = x >> 4;
                    qp_y_idx_shift = 4;
                    break;
                }

                const QP_STORE_T QP = QP_store[(y >> qp_y_idx_shift) * QP_stride + qp_x_idx];

                uint8_t* b8x8 = image + static_cast<ptrdiff_t>(stride) * y + x;
                const uint8_t* b10x10 = image + static_cast<ptrdiff_t>(stride) * (y - 1) + (x - 1);

                uint8_t min_val{255}, max_val{0};

                for (int r = 0; r < 8; ++r)
                {
                    for (int c = 0; c < 8; ++c)
                    {
                        const uint8_t p = b8x8[static_cast<ptrdiff_t>(r) * stride + c];
                        min_val = std::min(min_val, p);
                        max_val = std::max(max_val, p);
                    }
                }

                const uint8_t thr = (max_val + min_val + 1) >> 1;
                const int max_diff = QP >> 1;

                auto pavgb = [](int a, int b) { return (a + b + 1) >> 1; };

                for (int r = 0; r < 8; ++r)
                {
                    for (int c = 0; c < 8; ++c)
                    {
                        const uint8_t* p_row0 = &b10x10[static_cast<ptrdiff_t>(r + 0) * stride + c];
                        const uint8_t* p_row1 = &b10x10[static_cast<ptrdiff_t>(r + 1) * stride + c];
                        const uint8_t* p_row2 = &b10x10[static_cast<ptrdiff_t>(r + 2) * stride + c];

                        const int h_avg0 = pavgb(pavgb(p_row0[0], p_row0[2]), p_row0[1]);
                        const int h_avg1 = pavgb(pavgb(p_row1[0], p_row1[2]), p_row1[1]);
                        const int h_avg2 = pavgb(pavgb(p_row2[0], p_row2[2]), p_row2[1]);
                        int filtered_val = pavgb(pavgb(h_avg0, h_avg2), h_avg1);

                        const uint8_t orig_val = b10x10[static_cast<ptrdiff_t>(r + 1) * stride + (c + 1)];
                        filtered_val =
                            std::clamp(filtered_val, static_cast<int>(orig_val) - max_diff, static_cast<int>(orig_val) + max_diff);

                        b8x8filtered[static_cast<ptrdiff_t>(r) * 8 + c] = static_cast<uint8_t>(std::clamp(filtered_val, 0, 255));
                    }
                }

                for (int r = 0; r < 8; ++r)
                {
                    for (int c = 0; c < 8; ++c)
                    {
                        bool all_above{true};
                        bool all_below{true};

                        for (int dr = 0; dr < 3; ++dr)
                        {
                            for (int dc = 0; dc < 3; ++dc)
                            {
                                const uint8_t p = b10x10[static_cast<ptrdiff_t>(r + dr) * stride + (c + dc)];

                                if (p < thr)
                                    all_above = false;

                                if (p >= thr)
                                    all_below = false;
                            }
                        }

                        if (all_above || all_below)
                            b8x8[static_cast<ptrdiff_t>(r) * stride + c] = b8x8filtered[static_cast<ptrdiff_t>(r) * 8 + c];
                    }
                }
            }
        }
    }

} // anonymous namespace

void __stdcall fast_copy(const uint8_t* AVS_RESTRICT src, const int src_stride, uint8_t* AVS_RESTRICT dst, const int dst_stride,
    const int horizontal_size, int vertical_size) noexcept
{
    if (vertical_size <= 0)
        return;

    if (horizontal_size == src_stride && src_stride == dst_stride)
        memcpy(dst, src, static_cast<size_t>(horizontal_size) * vertical_size);
    else
    {
        const size_t bytes_to_copy = horizontal_size;

        for (int i = 0; i < vertical_size; ++i)
        {
            memcpy(dst, src, bytes_to_copy);
            dst += dst_stride;
            src += src_stride;
        }
    }
}

void postprocess(const uint8_t* const AVS_RESTRICT src[3], const int src_strides[3], uint8_t* const AVS_RESTRICT dst[3],
    const int dst_strides[3], const int width[2], int height[2], const QP_STORE_T* AVS_RESTRICT QP_store, int QP_stride, PPFlags mode,
    int moderate_h, int moderate_v, ChromaFormat format, bool iPP) noexcept
{
    const int num_planes = (format == ChromaFormat::Y8) ? 1 : 3;

    for (int plane_idx = 0; plane_idx < num_planes; ++plane_idx)
    {
        const bool is_luma = (plane_idx == 0);

        const uint8_t* p_src = src[plane_idx];
        uint8_t* p_dst = dst[plane_idx];
        const int src_stride = src_strides[plane_idx];
        const int dst_stride = dst_strides[plane_idx];
        const int plane_width = is_luma ? width[0] : width[1];
        const int plane_height = is_luma ? height[0] : height[1];

        const PPFlags deblock_h_flag = is_luma ? PPFlags::DeblockYH : PPFlags::DeblockCH;
        const PPFlags deblock_v_flag = is_luma ? PPFlags::DeblockYV : PPFlags::DeblockCV;
        const PPFlags dering_flag = is_luma ? PPFlags::DeringY : PPFlags::DeringC;

        int chroma_flag_for_qp;
        int qp_y_shift;

        if (is_luma)
        {
            chroma_flag_for_qp = 0;
            qp_y_shift = 4;
        }
        else
        {
            switch (format)
            {
            case ChromaFormat::YUV420:
                chroma_flag_for_qp = 1;
                qp_y_shift = 3;
                break;
            case ChromaFormat::YUV422:
                chroma_flag_for_qp = 2;
                qp_y_shift = 4;
                break;
            case ChromaFormat::YUV444:
                chroma_flag_for_qp = 3;
                qp_y_shift = 4;
                break;
            default:
                continue; // Should not happen with Y8 check above
            }
        }

        int v_size = plane_height;
        if (iPP)
        {
            v_size >>= 1;
        }

        // Deblocking loop
        for (int y = 0; y < v_size; y += 4)
        {
            if ((mode & PPFlags::DontCopy) != PPFlags::DontCopy)
            {
                if (!iPP)
                {
                    fast_copy(p_src + static_cast<ptrdiff_t>(y) * src_stride, src_stride, p_dst + static_cast<ptrdiff_t>(y) * dst_stride,
                        dst_stride, plane_width, 4);
                }
                else
                {
                    fast_copy(p_src + static_cast<ptrdiff_t>(y) * 2 * src_stride, src_stride,
                        p_dst + static_cast<ptrdiff_t>(y) * 2 * dst_stride, dst_stride, plane_width, 8);
                }
            }

            if ((mode & deblock_h_flag) == deblock_h_flag)
            {
                const auto* qp_row = &QP_store[(y >> qp_y_shift) * QP_stride];
                if (!iPP)
                {
                    deblock_horiz(p_dst + static_cast<ptrdiff_t>(y) * dst_stride, plane_width, dst_stride, qp_row, QP_stride,
                        chroma_flag_for_qp, moderate_h);
                }
                else
                {
                    deblock_horiz(p_dst + static_cast<ptrdiff_t>(y) * 2 * dst_stride, plane_width, dst_stride * 2,
                        &QP_store[(y >> qp_y_shift) * 2 * QP_stride], QP_stride * 2, chroma_flag_for_qp, moderate_h);
                    deblock_horiz(p_dst + (static_cast<ptrdiff_t>(y) * 2 + 1) * dst_stride, plane_width, dst_stride * 2,
                        &QP_store[((y >> qp_y_shift) * 2 + 1) * QP_stride], QP_stride * 2, chroma_flag_for_qp, moderate_h);
                }
            }

            if ((mode & deblock_v_flag) == deblock_v_flag)
            {
                if (y >= 12 && (y % 8) == 4)
                {
                    const int boundary_y = y - 4;
                    const auto* qp_row = &QP_store[(boundary_y >> qp_y_shift) * QP_stride];
                    if (!iPP)
                    {
                        deblock_vert(p_dst + static_cast<ptrdiff_t>(boundary_y) * dst_stride, plane_width, dst_stride, qp_row, QP_stride,
                            chroma_flag_for_qp, moderate_v);
                    }
                    else
                    {
                        deblock_vert(p_dst + static_cast<ptrdiff_t>(boundary_y) * 2 * dst_stride, plane_width, dst_stride * 2,
                            &QP_store[(boundary_y >> qp_y_shift) * 2 * QP_stride], QP_stride * 2, chroma_flag_for_qp, moderate_v);
                        deblock_vert(p_dst + (static_cast<ptrdiff_t>(boundary_y) * 2 + 1) * dst_stride, plane_width, dst_stride * 2,
                            &QP_store[((boundary_y >> qp_y_shift) * 2 + 1) * QP_stride], QP_stride * 2, chroma_flag_for_qp, moderate_v);
                    }
                }
            }
        }

        // Deringing
        if ((mode & dering_flag) == dering_flag)
        {
            if (!iPP)
            {
                dering(p_dst, plane_width, plane_height, dst_stride, QP_store, QP_stride, chroma_flag_for_qp);
            }
            else
            {
                dering(p_dst, plane_width, v_size, dst_stride * 2, QP_store, QP_stride * 2, chroma_flag_for_qp);
                dering(p_dst + dst_stride, plane_width, v_size, dst_stride * 2, QP_store + QP_stride, QP_stride * 2, chroma_flag_for_qp);
            }
        }
    }
}
