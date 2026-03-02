#pragma once

#include <tuple>

#include "PostProcess_common.h"
#include "vectorclass.h"

namespace SIMD_NAMESPACE
{
    template<typename T>
    [[nodiscard]] AVS_FORCEINLINE int deblock_horiz_useDC(const T* AVS_RESTRICT v, int stride, int moderate_h, int bit_depth) noexcept
    {
        const T* v_0{v};
        const T* v_1{v + stride};
        const T* v_2{v + 2 * stride};
        const T* v_3{v + 3 * stride};

        if constexpr (std::is_floating_point_v<T>)
        {
            const float thr{1.0f / 255.0f};

            Vec4f v1(v_0[1], v_1[1], v_2[1], v_3[1]);
            Vec4f v2(v_0[2], v_1[2], v_2[2], v_3[2]);
            Vec4f v3(v_0[3], v_1[3], v_2[3], v_3[3]);
            Vec4f v4(v_0[4], v_1[4], v_2[4], v_3[4]);
            Vec4f v5(v_0[5], v_1[5], v_2[5], v_3[5]);
            Vec4f v6(v_0[6], v_1[6], v_2[6], v_3[6]);
            Vec4f v7(v_0[7], v_1[7], v_2[7], v_3[7]);
            Vec4f v8(v_0[8], v_1[8], v_2[8], v_3[8]);

            Vec4f ones(1.0f);
            Vec4f zeros(0.0f);

            Vec4f sum{select(abs(v1 - v2) <= thr, ones, zeros)};
            sum += select(abs(v2 - v3) <= thr, ones, zeros);
            sum += select(abs(v3 - v4) <= thr, ones, zeros);
            sum += select(abs(v4 - v5) <= thr, ones, zeros);
            sum += select(abs(v5 - v6) <= thr, ones, zeros);
            sum += select(abs(v6 - v7) <= thr, ones, zeros);
            sum += select(abs(v7 - v8) <= thr, ones, zeros);

            return horizontal_add(sum) >= static_cast<float>(moderate_h);
        }
        else
        {
            const int thr{1 << (bit_depth - 8)};

            Vec4i v1(v_0[1], v_1[1], v_2[1], v_3[1]);
            Vec4i v2(v_0[2], v_1[2], v_2[2], v_3[2]);
            Vec4i v3(v_0[3], v_1[3], v_2[3], v_3[3]);
            Vec4i v4(v_0[4], v_1[4], v_2[4], v_3[4]);
            Vec4i v5(v_0[5], v_1[5], v_2[5], v_3[5]);
            Vec4i v6(v_0[6], v_1[6], v_2[6], v_3[6]);
            Vec4i v7(v_0[7], v_1[7], v_2[7], v_3[7]);
            Vec4i v8(v_0[8], v_1[8], v_2[8], v_3[8]);

            Vec4i ones(1);
            Vec4i zeros(0);

            Vec4i sum{select(abs(v1 - v2) <= thr, ones, zeros)};
            sum += select(abs(v2 - v3) <= thr, ones, zeros);
            sum += select(abs(v3 - v4) <= thr, ones, zeros);
            sum += select(abs(v4 - v5) <= thr, ones, zeros);
            sum += select(abs(v5 - v6) <= thr, ones, zeros);
            sum += select(abs(v6 - v7) <= thr, ones, zeros);
            sum += select(abs(v7 - v8) <= thr, ones, zeros);

            return horizontal_add(sum) >= moderate_h;
        }
    }

    template<typename T>
    AVS_FORCEINLINE void deblock_horiz_lpf9(T* AVS_RESTRICT v, int stride, int QP, int bit_depth, bool is_float_chroma) noexcept
    {
        const T* v_0{v};
        const T* v_1{v + stride};
        const T* v_2{v + 2 * stride};
        const T* v_3{v + 3 * stride};

        if constexpr (std::is_floating_point_v<T>)
        {
            const float qp_thresh{static_cast<float>(QP) / 255.0f};

            Vec4f v0(v_0[0], v_1[0], v_2[0], v_3[0]);
            Vec4f v1(v_0[1], v_1[1], v_2[1], v_3[1]);
            Vec4f v2(v_0[2], v_1[2], v_2[2], v_3[2]);
            Vec4f v3(v_0[3], v_1[3], v_2[3], v_3[3]);
            Vec4f v4(v_0[4], v_1[4], v_2[4], v_3[4]);
            Vec4f v5(v_0[5], v_1[5], v_2[5], v_3[5]);
            Vec4f v6(v_0[6], v_1[6], v_2[6], v_3[6]);
            Vec4f v7(v_0[7], v_1[7], v_2[7], v_3[7]);
            Vec4f v8(v_0[8], v_1[8], v_2[8], v_3[8]);
            Vec4f v9(v_0[9], v_1[9], v_2[9], v_3[9]);

            Vec4f p1{select(Vec4fb(abs(v0 - v1) < qp_thresh), v0, v1)};
            Vec4f p2{select(Vec4fb(abs(v8 - v9) < qp_thresh), v9, v8)};

            Vec4f out0{p1 * 6.0f + v1 * 4.0f + v2 * 2.0f + v3 * 2.0f + v4 + v5};
            Vec4f out1{p1 * 4.0f + v1 * 2.0f + v2 * 4.0f + v3 * 2.0f + v4 * 2.0f + v5 + v6};
            Vec4f out2{p1 * 2.0f + v1 * 2.0f + v2 * 2.0f + v3 * 4.0f + v4 * 2.0f + v5 * 2.0f + v6 + v7};
            Vec4f out3{p1 + v1 + v2 * 2.0f + v3 * 2.0f + v4 * 4.0f + v5 * 2.0f + v6 * 2.0f + v7 + v8};
            Vec4f out4{v1 + v2 + v3 * 2.0f + v4 * 2.0f + v5 * 4.0f + v6 * 2.0f + v7 * 2.0f + v8 + p2};
            Vec4f out5{v2 + v3 + v4 * 2.0f + v5 * 2.0f + v6 * 4.0f + v7 * 2.0f + v8 * 2.0f + p2 * 2.0f};
            Vec4f out6{v3 + v4 + v5 * 2.0f + v6 * 2.0f + v7 * 4.0f + v8 * 2.0f + p2 * 4.0f};
            Vec4f out7{v4 + v5 + v6 * 2.0f + v7 * 2.0f + v8 * 4.0f + p2 * 6.0f};

            const float inv16{1.0f / 16.0f};
            Vec4f vmin(is_float_chroma ? -0.5f : 0.0f);
            Vec4f vmax(is_float_chroma ? 0.5f : 1.0f);

            out0 = max(vmin, min(out0 * inv16, vmax));
            out1 = max(vmin, min(out1 * inv16, vmax));
            out2 = max(vmin, min(out2 * inv16, vmax));
            out3 = max(vmin, min(out3 * inv16, vmax));
            out4 = max(vmin, min(out4 * inv16, vmax));
            out5 = max(vmin, min(out5 * inv16, vmax));
            out6 = max(vmin, min(out6 * inv16, vmax));
            out7 = max(vmin, min(out7 * inv16, vmax));

            T* out_0{v};
            T* out_1{v + stride};
            T* out_2{v + 2 * stride};
            T* out_3{v + 3 * stride};

            out_0[1] = out0[0];
            out_1[1] = out0[1];
            out_2[1] = out0[2];
            out_3[1] = out0[3];
            out_0[2] = out1[0];
            out_1[2] = out1[1];
            out_2[2] = out1[2];
            out_3[2] = out1[3];
            out_0[3] = out2[0];
            out_1[3] = out2[1];
            out_2[3] = out2[2];
            out_3[3] = out2[3];
            out_0[4] = out3[0];
            out_1[4] = out3[1];
            out_2[4] = out3[2];
            out_3[4] = out3[3];
            out_0[5] = out4[0];
            out_1[5] = out4[1];
            out_2[5] = out4[2];
            out_3[5] = out4[3];
            out_0[6] = out5[0];
            out_1[6] = out5[1];
            out_2[6] = out5[2];
            out_3[6] = out5[3];
            out_0[7] = out6[0];
            out_1[7] = out6[1];
            out_2[7] = out6[2];
            out_3[7] = out6[3];
            out_0[8] = out7[0];
            out_1[8] = out7[1];
            out_2[8] = out7[2];
            out_3[8] = out7[3];
        }
        else
        {
            const int qp_scaled{QP * (1 << (bit_depth - 8))};

            Vec4i v0(v_0[0], v_1[0], v_2[0], v_3[0]);
            Vec4i v1(v_0[1], v_1[1], v_2[1], v_3[1]);
            Vec4i v2(v_0[2], v_1[2], v_2[2], v_3[2]);
            Vec4i v3(v_0[3], v_1[3], v_2[3], v_3[3]);
            Vec4i v4(v_0[4], v_1[4], v_2[4], v_3[4]);
            Vec4i v5(v_0[5], v_1[5], v_2[5], v_3[5]);
            Vec4i v6(v_0[6], v_1[6], v_2[6], v_3[6]);
            Vec4i v7(v_0[7], v_1[7], v_2[7], v_3[7]);
            Vec4i v8(v_0[8], v_1[8], v_2[8], v_3[8]);
            Vec4i v9(v_0[9], v_1[9], v_2[9], v_3[9]);

            Vec4i p1{select(Vec4ib(abs(v0 - v1) < qp_scaled), v0, v1)};
            Vec4i p2{select(Vec4ib(abs(v8 - v9) < qp_scaled), v9, v8)};

            Vec4i p1_2{p1 << 1};
            Vec4i p1_4{p1 << 2};
            Vec4i v1_2{v1 << 1};
            Vec4i v1_4{v1 << 2};
            Vec4i v2_2{v2 << 1};
            Vec4i v2_4{v2 << 2};
            Vec4i v3_2{v3 << 1};
            Vec4i v3_4{v3 << 2};
            Vec4i v4_2{v4 << 1};
            Vec4i v4_4{v4 << 2};
            Vec4i v5_2{v5 << 1};
            Vec4i v5_4{v5 << 2};
            Vec4i v6_2{v6 << 1};
            Vec4i v6_4{v6 << 2};
            Vec4i v7_2{v7 << 1};
            Vec4i v7_4{v7 << 2};
            Vec4i v8_2{v8 << 1};
            Vec4i v8_4{v8 << 2};
            Vec4i p2_2{p2 << 1};
            Vec4i p2_4{p2 << 2};

            Vec4i out0{p1_4 + p1_2 + v1_4 + v2_2 + v3_2 + v4 + v5 + 8};
            Vec4i out1{p1_4 + v1_2 + v2_4 + v3_2 + v4_2 + v5 + v6 + 8};
            Vec4i out2{p1_2 + v1_2 + v2_2 + v3_4 + v4_2 + v5_2 + v6 + v7 + 8};
            Vec4i out3{p1 + v1 + v2_2 + v3_2 + v4_4 + v5_2 + v6_2 + v7 + v8 + 8};
            Vec4i out4{v1 + v2 + v3_2 + v4_2 + v5_4 + v6_2 + v7_2 + v8 + p2 + 8};
            Vec4i out5{v2 + v3 + v4_2 + v5_2 + v6_4 + v7_2 + v8_2 + p2_2 + 8};
            Vec4i out6{v3 + v4 + v5_2 + v6_2 + v7_4 + v8_2 + p2_4 + 8};
            Vec4i out7{v4 + v5 + v6_2 + v7_2 + v8_4 + p2_4 + p2_2 + 8};

            Vec4i vmin(0);
            Vec4i vmax((1 << bit_depth) - 1);

            out0 = max(vmin, min(out0 >> 4, vmax));
            out1 = max(vmin, min(out1 >> 4, vmax));
            out2 = max(vmin, min(out2 >> 4, vmax));
            out3 = max(vmin, min(out3 >> 4, vmax));
            out4 = max(vmin, min(out4 >> 4, vmax));
            out5 = max(vmin, min(out5 >> 4, vmax));
            out6 = max(vmin, min(out6 >> 4, vmax));
            out7 = max(vmin, min(out7 >> 4, vmax));

            T* out_0{v};
            T* out_1{v + stride};
            T* out_2{v + 2 * stride};
            T* out_3{v + 3 * stride};

            out_0[1] = static_cast<T>(out0[0]);
            out_1[1] = static_cast<T>(out0[1]);
            out_2[1] = static_cast<T>(out0[2]);
            out_3[1] = static_cast<T>(out0[3]);
            out_0[2] = static_cast<T>(out1[0]);
            out_1[2] = static_cast<T>(out1[1]);
            out_2[2] = static_cast<T>(out1[2]);
            out_3[2] = static_cast<T>(out1[3]);
            out_0[3] = static_cast<T>(out2[0]);
            out_1[3] = static_cast<T>(out2[1]);
            out_2[3] = static_cast<T>(out2[2]);
            out_3[3] = static_cast<T>(out2[3]);
            out_0[4] = static_cast<T>(out3[0]);
            out_1[4] = static_cast<T>(out3[1]);
            out_2[4] = static_cast<T>(out3[2]);
            out_3[4] = static_cast<T>(out3[3]);
            out_0[5] = static_cast<T>(out4[0]);
            out_1[5] = static_cast<T>(out4[1]);
            out_2[5] = static_cast<T>(out4[2]);
            out_3[5] = static_cast<T>(out4[3]);
            out_0[6] = static_cast<T>(out5[0]);
            out_1[6] = static_cast<T>(out5[1]);
            out_2[6] = static_cast<T>(out5[2]);
            out_3[6] = static_cast<T>(out5[3]);
            out_0[7] = static_cast<T>(out6[0]);
            out_1[7] = static_cast<T>(out6[1]);
            out_2[7] = static_cast<T>(out6[2]);
            out_3[7] = static_cast<T>(out6[3]);
            out_0[8] = static_cast<T>(out7[0]);
            out_1[8] = static_cast<T>(out7[1]);
            out_2[8] = static_cast<T>(out7[2]);
            out_3[8] = static_cast<T>(out7[3]);
        }
    }

    template<typename T>
    AVS_FORCEINLINE void deblock_horiz_default_filter(T* AVS_RESTRICT v, int stride, int QP, int bit_depth, bool is_float_chroma) noexcept
    {
        const T* v_0{v};
        const T* v_1{v + stride};
        const T* v_2{v + 2 * stride};
        const T* v_3{v + 3 * stride};

        if constexpr (std::is_floating_point_v<T>)
        {
            const float qp_thresh{static_cast<float>(QP) / 255.0f};
            Vec4f vmin(is_float_chroma ? -0.5f : 0.0f);
            Vec4f vmax(is_float_chroma ? 0.5f : 1.0f);

            auto gather4f = [](const T* p0, const T* p1, const T* p2, const T* p3) { return Vec4f(*p0, *p1, *p2, *p3); };
            Vec4f v1{gather4f(v_0 + 1, v_1 + 1, v_2 + 1, v_3 + 1)};
            Vec4f v2{gather4f(v_0 + 2, v_1 + 2, v_2 + 2, v_3 + 2)};
            Vec4f v3{gather4f(v_0 + 3, v_1 + 3, v_2 + 3, v_3 + 3)};
            Vec4f v4{gather4f(v_0 + 4, v_1 + 4, v_2 + 4, v_3 + 4)};
            Vec4f v5{gather4f(v_0 + 5, v_1 + 5, v_2 + 5, v_3 + 5)};
            Vec4f v6{gather4f(v_0 + 6, v_1 + 6, v_2 + 6, v_3 + 6)};
            Vec4f v7{gather4f(v_0 + 7, v_1 + 7, v_2 + 7, v_3 + 7)};
            Vec4f v8{gather4f(v_0 + 8, v_1 + 8, v_2 + 8, v_3 + 8)};

            Vec4f q1{v4 - v5};
            Vec4f q{q1 * 0.5f};
            Vec4f a3_0{(v3 - v6) * 2.0f - q1 * 5.0f};
            Vec4f aa3_0{abs(a3_0)};

            Vec4fb cond1{aa3_0 < (8.0f * qp_thresh)};
            Vec4f a3_1{abs((v3 - v2) * 5.0f + (v1 - v4) * 2.0f)};
            Vec4f a3_2{abs((v7 - v8) * 5.0f + (v5 - v8) * 2.0f)};
            Vec4f d{aa3_0 - min(a3_1, a3_2)};

            Vec4fb cond2{d > (6.0f / 255.0f)};
            Vec4fb cond3{(q > 0.0f & a3_0 < 0.0f) | (q < 0.0f & a3_0 > 0.0f)};

            Vec4f delta_mag{min(d, abs(q))};
            Vec4f delta{select(cond3, select(q > 0.0f, -delta_mag, delta_mag), 0.0f)};

            delta = select((q != 0.0f) & cond1 & cond2, delta, 0.0f);

            Vec4f out4{max(vmin, min(v4 + delta, vmax))};
            Vec4f out5{max(vmin, min(v5 - delta, vmax))};

            auto scatter4f = [](Vec4f vec, T* p0, T* p1, T* p2, T* p3) {
                *p0 = vec[0];
                *p1 = vec[1];
                *p2 = vec[2];
                *p3 = vec[3];
            };
            scatter4f(out4, v + 4, v + stride + 4, v + 2 * stride + 4, v + 3 * stride + 4);
            scatter4f(out5, v + 5, v + stride + 5, v + 2 * stride + 5, v + 3 * stride + 5);
        }
        else
        {
            const int scale{1 << (bit_depth - 8)};
            const int qp_scaled{QP * scale};
            Vec4i vmin(0);
            Vec4i vmax((1 << bit_depth) - 1);

            auto gather4i = [](const T* p0, const T* p1, const T* p2, const T* p3) { return Vec4i(*p0, *p1, *p2, *p3); };
            Vec4i v1{gather4i(v_0 + 1, v_1 + 1, v_2 + 1, v_3 + 1)};
            Vec4i v2{gather4i(v_0 + 2, v_1 + 2, v_2 + 2, v_3 + 2)};
            Vec4i v3{gather4i(v_0 + 3, v_1 + 3, v_2 + 3, v_3 + 3)};
            Vec4i v4{gather4i(v_0 + 4, v_1 + 4, v_2 + 4, v_3 + 4)};
            Vec4i v5{gather4i(v_0 + 5, v_1 + 5, v_2 + 5, v_3 + 5)};
            Vec4i v6{gather4i(v_0 + 6, v_1 + 6, v_2 + 6, v_3 + 6)};
            Vec4i v7{gather4i(v_0 + 7, v_1 + 7, v_2 + 7, v_3 + 7)};
            Vec4i v8{gather4i(v_0 + 8, v_1 + 8, v_2 + 8, v_3 + 8)};

            Vec4i q1{v4 - v5};
            Vec4i q{q1 / 2};
            Vec4i a3_0{(v3 - v6) * 2 - q1 * 5};
            Vec4i aa3_0{abs(a3_0)};

            Vec4ib cond1{aa3_0 < (8 * qp_scaled)};
            Vec4i a3_1{abs((v3 - v2) * 5 + (v1 - v4) * 2)};
            Vec4i a3_2{abs((v7 - v8) * 5 + (v5 - v8) * 2)};
            Vec4i d{aa3_0 - min(a3_1, a3_2)};

            Vec4ib cond2{d > (6 * scale)};
            Vec4ib cond3{(q > 0 & a3_0 < 0) | (q < 0 & a3_0 > 0)};

            Vec4i delta_mag{min(d, abs(q))};
            Vec4i delta{select(cond3, select(q > 0, -delta_mag, delta_mag), 0)};

            delta = select((q != 0) & cond1 & cond2, delta, 0);

            Vec4i out4{max(vmin, min(v4 + delta, vmax))};
            Vec4i out5{max(vmin, min(v5 - delta, vmax))};

            auto scatter4i = [](Vec4i vec, T* p0, T* p1, T* p2, T* p3) {
                *p0 = static_cast<T>(vec[0]);
                *p1 = static_cast<T>(vec[1]);
                *p2 = static_cast<T>(vec[2]);
                *p3 = static_cast<T>(vec[3]);
            };
            scatter4i(out4, v + 4, v + stride + 4, v + 2 * stride + 4, v + 3 * stride + 4);
            scatter4i(out5, v + 5, v + stride + 5, v + 2 * stride + 5, v + 3 * stride + 5);
        }
    }

    template<typename T>
    void deblock_horiz(T* AVS_RESTRICT image, int width, int stride, const QP_STORE_T* AVS_RESTRICT QP_store, int QP_stride,
        int chroma_flag_for_qp, int moderate_h, int bit_depth, bool is_float_chroma) noexcept
    {
        for (int x{8}; x < width; x += 8)
        {
            const int QP{QP_store[(chroma_flag_for_qp == 0 || chroma_flag_for_qp == 3) ? (x / 16) : (x / 8)]};
            T* v{image + x - 5};

            if (deblock_horiz_useDC<T>(v, stride, moderate_h, bit_depth))
            {
                if (deblock_horiz_DC_on<T>(v, stride, QP, bit_depth))
                    deblock_horiz_lpf9<T>(v, stride, QP, bit_depth, is_float_chroma);
            }
            else
                deblock_horiz_default_filter<T>(v, stride, QP, bit_depth, is_float_chroma);
        }
    }

    template<typename T>
    [[nodiscard]] AVS_FORCEINLINE int deblock_vert_useDC(const T* AVS_RESTRICT v, int stride, int moderate_v, int bit_depth) noexcept
    {
        if constexpr (std::is_floating_point_v<T>)
        {
            const float thr{1.0f / 255.0f};
            Vec4f sum_lo(0.0f), sum_hi(0.0f);
            Vec4f ones(1.0f), zeros(0.0f);

            for (int i{1}; i < 8; ++i)
            {
                const T* p1{v + i * stride};
                const T* p2{v + (i + 1) * stride};

                Vec4f v1_lo;
                v1_lo.load(p1);
                Vec4f v1_hi;
                v1_hi.load(p1 + 4);
                Vec4f v2_lo;
                v2_lo.load(p2);
                Vec4f v2_hi;
                v2_hi.load(p2 + 4);

                sum_lo += select(abs(v1_lo - v2_lo) <= thr, ones, zeros);
                sum_hi += select(abs(v1_hi - v2_hi) <= thr, ones, zeros);
            }

            return horizontal_add(sum_lo + sum_hi) > static_cast<float>(moderate_v);
        }
        else
        {
            const int thr{1 << (bit_depth - 8)};
            Vec4i sum_lo(0), sum_hi(0);
            Vec4i ones(1), zeros(0);

            auto load4i = [](const T* p) { return Vec4i(p[0], p[1], p[2], p[3]); };

            for (int i{1}; i < 8; ++i)
            {
                const T* p1{v + i * stride};
                const T* p2{v + (i + 1) * stride};

                Vec4i v1_lo{load4i(p1)};
                Vec4i v1_hi{load4i(p1 + 4)};
                Vec4i v2_lo{load4i(p2)};
                Vec4i v2_hi{load4i(p2 + 4)};

                sum_lo += select(abs(v1_lo - v2_lo) <= thr, ones, zeros);
                sum_hi += select(abs(v1_hi - v2_hi) <= thr, ones, zeros);
            }

            return horizontal_add(sum_lo + sum_hi) > moderate_v;
        }
    }

    template<typename T, typename T_intermediate>
    AVS_FORCEINLINE void deblock_vert_lpf9(const T_intermediate* AVS_RESTRICT v_local, const T_intermediate* AVS_RESTRICT p1p2,
        T* AVS_RESTRICT v, int stride, int bit_depth, bool is_float_chroma) noexcept
    {
        if constexpr (std::is_floating_point_v<T>)
        {
            Vec4f vmin(is_float_chroma ? -0.5f : 0.0f);
            Vec4f vmax(is_float_chroma ? 0.5f : 1.0f);
            const float inv16{1.0f / 16.0f};

            for (int half{0}; half < 2; ++half)
            {
                int offset{half * 4};
                Vec4f p1;
                p1.load(p1p2 + offset);
                Vec4f p2;
                p2.load(p1p2 + 8 + offset);

                Vec4f v1;
                v1.load(v_local + 0 * 8 + offset);
                Vec4f v2;
                v2.load(v_local + 1 * 8 + offset);
                Vec4f v3;
                v3.load(v_local + 2 * 8 + offset);
                Vec4f v4;
                v4.load(v_local + 3 * 8 + offset);
                Vec4f v5;
                v5.load(v_local + 4 * 8 + offset);
                Vec4f v6;
                v6.load(v_local + 5 * 8 + offset);
                Vec4f v7;
                v7.load(v_local + 6 * 8 + offset);
                Vec4f v8;
                v8.load(v_local + 7 * 8 + offset);

                Vec4f psum{p1 * 3.0f + v1 + v2 + v3 + v4};
                Vec4f s;

                s = ((psum + v1) * 2.0f) - (v4 - v5);
                max(vmin, min(s * inv16, vmax)).store(v + 1 * stride + offset);

                psum += v5 - p1;
                s = ((psum + v2) * 2.0f) - (v5 - v6);
                max(vmin, min(s * inv16, vmax)).store(v + 2 * stride + offset);

                psum += v6 - p1;
                s = ((psum + v3) * 2.0f) - (v6 - v7);
                max(vmin, min(s * inv16, vmax)).store(v + 3 * stride + offset);

                psum += v7 - p1;
                s = ((psum + v4) * 2.0f) - (v7 - v8) + p1 - v1;
                max(vmin, min(s * inv16, vmax)).store(v + 4 * stride + offset);

                psum += v8 - v1;
                s = ((psum + v5) * 2.0f) - (v8 - p2) + v1 - v2;
                max(vmin, min(s * inv16, vmax)).store(v + 5 * stride + offset);

                psum += p2 - v2;
                s = ((psum + v6) * 2.0f) + v2 - v3;
                max(vmin, min(s * inv16, vmax)).store(v + 6 * stride + offset);

                psum += p2 - v3;
                s = ((psum + v7) * 2.0f) + v3 - v4;
                max(vmin, min(s * inv16, vmax)).store(v + 7 * stride + offset);

                psum += p2 - v4;
                s = ((psum + v8) * 2.0f) + v4 - v5;
                max(vmin, min(s * inv16, vmax)).store(v + 8 * stride + offset);
            }
        }
        else
        {
            Vec4i vmin(0);
            Vec4i vmax((1 << bit_depth) - 1);

            auto load4i = [](const T_intermediate* p) {
                if constexpr (sizeof(T_intermediate) == 2)
                {
                    return Vec4i().load_4us(p);
                }
                else
                {
                    Vec4i res;
                    res.load(p);
                    return res;
                }
            };

            auto store4i = [](Vec4i vec, T* p) {
                p[0] = static_cast<T>(vec[0]);
                p[1] = static_cast<T>(vec[1]);
                p[2] = static_cast<T>(vec[2]);
                p[3] = static_cast<T>(vec[3]);
            };

            for (int half{0}; half < 2; ++half)
            {
                int offset{half * 4};
                Vec4i p1{load4i(p1p2 + offset)};
                Vec4i p2{load4i(p1p2 + 8 + offset)};

                Vec4i v1{load4i(v_local + 0 * 8 + offset)};
                Vec4i v2{load4i(v_local + 1 * 8 + offset)};
                Vec4i v3{load4i(v_local + 2 * 8 + offset)};
                Vec4i v4{load4i(v_local + 3 * 8 + offset)};
                Vec4i v5{load4i(v_local + 4 * 8 + offset)};
                Vec4i v6{load4i(v_local + 5 * 8 + offset)};
                Vec4i v7{load4i(v_local + 6 * 8 + offset)};
                Vec4i v8{load4i(v_local + 7 * 8 + offset)};

                Vec4i psum{p1 * 3 + v1 + v2 + v3 + v4 + 4};
                Vec4i s;

                s = ((psum + v1) << 1) - (v4 - v5);
                store4i(max(vmin, min(s >> 4, vmax)), v + 1 * stride + offset);

                psum += v5 - p1;
                s = ((psum + v2) << 1) - (v5 - v6);
                store4i(max(vmin, min(s >> 4, vmax)), v + 2 * stride + offset);

                psum += v6 - p1;
                s = ((psum + v3) << 1) - (v6 - v7);
                store4i(max(vmin, min(s >> 4, vmax)), v + 3 * stride + offset);

                psum += v7 - p1;
                s = ((psum + v4) << 1) - (v7 - v8) + p1 - v1;
                store4i(max(vmin, min(s >> 4, vmax)), v + 4 * stride + offset);

                psum += v8 - v1;
                s = ((psum + v5) << 1) - (v8 - p2) + v1 - v2;
                store4i(max(vmin, min(s >> 4, vmax)), v + 5 * stride + offset);

                psum += p2 - v2;
                s = ((psum + v6) << 1) + v2 - v3;
                store4i(max(vmin, min(s >> 4, vmax)), v + 6 * stride + offset);

                psum += p2 - v3;
                s = ((psum + v7) << 1) + v3 - v4;
                store4i(max(vmin, min(s >> 4, vmax)), v + 7 * stride + offset);

                psum += p2 - v4;
                s = ((psum + v8) << 1) + v4 - v5;
                store4i(max(vmin, min(s >> 4, vmax)), v + 8 * stride + offset);
            }
        }
    }

    template<typename T>
    AVS_FORCEINLINE void deblock_vert_default_filter(T* AVS_RESTRICT v, int stride, int QP, int bit_depth, bool is_float_chroma) noexcept
    {
        if constexpr (std::is_floating_point_v<T>)
        {
            const float qp_thresh{static_cast<float>(QP) / 255.0f};
            Vec4f vmin(is_float_chroma ? -0.5f : 0.0f);
            Vec4f vmax(is_float_chroma ? 0.5f : 1.0f);

            for (int half{0}; half < 2; ++half)
            {
                const int offset{half * 4};
                Vec4f v1;
                v1.load(v + 1 * stride + offset);
                Vec4f v2;
                v2.load(v + 2 * stride + offset);
                Vec4f v3;
                v3.load(v + 3 * stride + offset);
                Vec4f v4;
                v4.load(v + 4 * stride + offset);
                Vec4f v5;
                v5.load(v + 5 * stride + offset);
                Vec4f v6;
                v6.load(v + 6 * stride + offset);
                Vec4f v7;
                v7.load(v + 7 * stride + offset);
                Vec4f v8;
                v8.load(v + 8 * stride + offset);

                Vec4f q1{v4 - v5};
                Vec4f q{q1 * 0.5f};
                Vec4f a3_0{(v3 - v6) * 2.0f - q1 * 5.0f};
                Vec4f aa3_0{abs(a3_0)};

                Vec4fb cond1{aa3_0 < (8.0f * qp_thresh)};
                Vec4f a3_1{abs((v3 - v2) * 5.0f + (v1 - v4) * 2.0f)};
                Vec4f a3_2{abs((v7 - v6) * 5.0f + (v5 - v8) * 2.0f)};
                Vec4f d_abs{aa3_0 - min(a3_1, a3_2)};

                Vec4fb cond2{d_abs > 0.0f};
                Vec4f d{(d_abs * 5.0f + (32.0f / 255.0f)) * (1.0f / 64.0f)};

                Vec4fb cond3{(q > 0.0f & a3_0 < 0.0f) | (q < 0.0f & a3_0 >= 0.0f)};
                Vec4f delta_mag{min(d, abs(q))};
                Vec4f final_delta{select(cond3, select(q > 0.0f, delta_mag, -delta_mag), 0.0f)};

                final_delta = select((q != 0.0f) & cond1 & cond2, final_delta, 0.0f);

                max(vmin, min(v4 - final_delta, vmax)).store(v + 4 * stride + offset);
                max(vmin, min(v5 + final_delta, vmax)).store(v + 5 * stride + offset);
            }
        }
        else
        {
            const int scale{1 << (bit_depth - 8)};
            const int qp_scaled{QP * scale};
            Vec4i vmin(0);
            Vec4i vmax((1 << bit_depth) - 1);

            auto load4i = [](const T* p) { return Vec4i(p[0], p[1], p[2], p[3]); };

            auto store4i = [](Vec4i vec, T* p) {
                p[0] = static_cast<T>(vec[0]);
                p[1] = static_cast<T>(vec[1]);
                p[2] = static_cast<T>(vec[2]);
                p[3] = static_cast<T>(vec[3]);
            };

            for (int half{0}; half < 2; ++half)
            {
                const int offset{half * 4};
                Vec4i v1{load4i(v + 1 * stride + offset)};
                Vec4i v2{load4i(v + 2 * stride + offset)};
                Vec4i v3{load4i(v + 3 * stride + offset)};
                Vec4i v4{load4i(v + 4 * stride + offset)};
                Vec4i v5{load4i(v + 5 * stride + offset)};
                Vec4i v6{load4i(v + 6 * stride + offset)};
                Vec4i v7{load4i(v + 7 * stride + offset)};
                Vec4i v8{load4i(v + 8 * stride + offset)};

                Vec4i q1{v4 - v5};
                Vec4i q{q1 / 2};
                Vec4i a3_0{(v3 - v6) * 2 - q1 * 5};
                Vec4i aa3_0{abs(a3_0)};

                Vec4ib cond1{aa3_0 < (8 * qp_scaled)};
                Vec4i a3_1{abs((v3 - v2) * 5 + (v1 - v4) * 2)};
                Vec4i a3_2{abs((v7 - v6) * 5 + (v5 - v8) * 2)};
                Vec4i d_abs{aa3_0 - min(a3_1, a3_2)};

                Vec4ib cond2{d_abs > 0};
                Vec4i d{(d_abs * 5 + (32 * scale)) >> 6};

                Vec4ib cond3{(q > 0 & a3_0 < 0) | (q < 0 & a3_0 >= 0)};
                Vec4i delta_mag{min(d, abs(q))};
                Vec4i final_delta{select(cond3, select(q > 0, delta_mag, -delta_mag), 0)};

                final_delta = select((q != 0) & cond1 & cond2, final_delta, 0);

                store4i(max(vmin, min(v4 - final_delta, vmax)), v + 4 * stride + offset);
                store4i(max(vmin, min(v5 + final_delta, vmax)), v + 5 * stride + offset);
            }
        }
    }

    template<typename T>
    void deblock_vert(T* AVS_RESTRICT image, int width, int stride, const QP_STORE_T* AVS_RESTRICT QP_store, int QP_stride,
        int chroma_flag_for_qp, int moderate_v, int bit_depth, bool is_float_chroma) noexcept
    {
        using intermediate_t = typename DeblockVertTypes<T>::intermediate;
        alignas(16) std::array<intermediate_t, 8 * 8> v_local;
        alignas(16) std::array<intermediate_t, 2 * 8> p1p2;

        for (int Bx{0}; Bx < width; Bx += 8)
        {
            const int QP{QP_store[(chroma_flag_for_qp == 0 || chroma_flag_for_qp == 3) ? (Bx / 16) : (Bx / 8)]};
            T* v{image - 5 * stride + Bx};

            if (deblock_vert_useDC<T>(v, stride, moderate_v, bit_depth))
            {
                if (deblock_vert_DC_on<T>(v, stride, QP, bit_depth))
                {
                    deblock_vert_copy_and_unpack<T, intermediate_t>(stride, v + stride, v_local.data(), 8);
                    deblock_vert_choose_p1p2<T, intermediate_t>(v, stride, p1p2.data(), QP, bit_depth);
                    deblock_vert_lpf9<T, intermediate_t>(v_local.data(), p1p2.data(), v, stride, bit_depth, is_float_chroma);
                }
            }
            else
                deblock_vert_default_filter<T>(v, stride, QP, bit_depth, is_float_chroma);
        }
    }

    template<typename T>
    void dering(T* AVS_RESTRICT image, int width, int height, int stride, const QP_STORE_T* AVS_RESTRICT QP_store, int QP_stride,
        int chroma_flag_for_qp, int bit_depth, bool is_float_chroma) noexcept
    {
        for (int y = 8; y < height - 8; y += 8)
        {
            for (int x = 8; x < width - 8; x += 8)
            {
                int qp_x_idx;
                int qp_y_idx_shift;

                switch (chroma_flag_for_qp)
                {
                case 1:
                    qp_x_idx = x >> 3;
                    qp_y_idx_shift = 3;
                    break;
                case 2:
                    qp_x_idx = x >> 3;
                    qp_y_idx_shift = 4;
                    break;
                default:
                    qp_x_idx = x >> 4;
                    qp_y_idx_shift = 4;
                    break;
                }

                const QP_STORE_T QP{QP_store[(y >> qp_y_idx_shift) * QP_stride + qp_x_idx]};
                T* b8x8{image + static_cast<ptrdiff_t>(stride) * y + x};
                const T* b10x10{image + static_cast<ptrdiff_t>(stride) * (y - 1) + (x - 1)};

                if constexpr (std::is_floating_point_v<T>)
                {
                    const float min_clamp{is_float_chroma ? -0.5f : 0.0f};
                    const float max_clamp{is_float_chroma ? 0.5f : 1.0f};
                    Vec4f vmin_clamp(min_clamp), vmax_clamp(max_clamp);

                    Vec4f min_lo;
                    min_lo.load(b8x8);
                    Vec4f min_hi;
                    min_hi.load(b8x8 + 4);
                    Vec4f max_lo{min_lo};
                    Vec4f max_hi{min_hi};

                    for (int r{1}; r < 8; ++r)
                    {
                        Vec4f p_lo;
                        p_lo.load(b8x8 + r * stride);
                        Vec4f p_hi;
                        p_hi.load(b8x8 + r * stride + 4);
                        min_lo = min(min_lo, p_lo);
                        min_hi = min(min_hi, p_hi);
                        max_lo = max(max_lo, p_lo);
                        max_hi = max(max_hi, p_hi);
                    }
                    float thr{(horizontal_max(max(max_lo, max_hi)) + horizontal_min(min(min_lo, min_hi))) * 0.5f};
                    float max_diff{(static_cast<float>(QP) / 255.0f) * 0.5f};
                    Vec4f vthr(thr), vmax_diff(max_diff);

                    auto h_avg_row = [](const T* p) {
                        Vec4f p0_lo;
                        p0_lo.load(p + 0);
                        Vec4f p0_hi;
                        p0_hi.load(p + 4);
                        Vec4f p1_lo;
                        p1_lo.load(p + 1);
                        Vec4f p1_hi;
                        p1_hi.load(p + 5);
                        Vec4f p2_lo;
                        p2_lo.load(p + 2);
                        Vec4f p2_hi;
                        p2_hi.load(p + 6);
                        return std::make_pair(((p0_lo + p2_lo) * 0.5f + p1_lo) * 0.5f, ((p0_hi + p2_hi) * 0.5f + p1_hi) * 0.5f);
                    };

                    auto check_row = [&](const T* p) {
                        Vec4f p0_lo;
                        p0_lo.load(p + 0);
                        Vec4f p0_hi;
                        p0_hi.load(p + 4);
                        Vec4f p1_lo;
                        p1_lo.load(p + 1);
                        Vec4f p1_hi;
                        p1_hi.load(p + 5);
                        Vec4f p2_lo;
                        p2_lo.load(p + 2);
                        Vec4f p2_hi;
                        p2_hi.load(p + 6);
                        Vec4fb below_lo{(p0_lo < vthr) & (p1_lo < vthr) & (p2_lo < vthr)};
                        Vec4fb above_lo{(p0_lo >= vthr) & (p1_lo >= vthr) & (p2_lo >= vthr)};
                        Vec4fb below_hi{(p0_hi < vthr) & (p1_hi < vthr) & (p2_hi < vthr)};
                        Vec4fb above_hi{(p0_hi >= vthr) & (p1_hi >= vthr) & (p2_hi >= vthr)};
                        return std::make_tuple(below_lo, above_lo, below_hi, above_hi);
                    };

                    Vec4f f_lo_arr[8];
                    Vec4f f_hi_arr[8];

                    for (int r{0}; r < 8; ++r)
                    {
                        auto [h0_lo, h0_hi] = h_avg_row(b10x10 + (r + 0) * stride);
                        auto [h1_lo, h1_hi] = h_avg_row(b10x10 + (r + 1) * stride);
                        auto [h2_lo, h2_hi] = h_avg_row(b10x10 + (r + 2) * stride);

                        Vec4f f_lo{((h0_lo + h2_lo) * 0.5f + h1_lo) * 0.5f};
                        Vec4f f_hi{((h0_hi + h2_hi) * 0.5f + h1_hi) * 0.5f};

                        Vec4f orig_lo;
                        orig_lo.load(b10x10 + (r + 1) * stride + 1);
                        Vec4f orig_hi;
                        orig_hi.load(b10x10 + (r + 1) * stride + 5);

                        f_lo_arr[r] = max(vmin_clamp, min(max(orig_lo - vmax_diff, min(f_lo, orig_lo + vmax_diff)), vmax_clamp));
                        f_hi_arr[r] = max(vmin_clamp, min(max(orig_hi - vmax_diff, min(f_hi, orig_hi + vmax_diff)), vmax_clamp));
                    }

                    for (int r{0}; r < 8; ++r)
                    {
                        auto [below0_lo, above0_lo, below0_hi, above0_hi] = check_row(b10x10 + (r + 0) * stride);
                        auto [below1_lo, above1_lo, below1_hi, above1_hi] = check_row(b10x10 + (r + 1) * stride);
                        auto [below2_lo, above2_lo, below2_hi, above2_hi] = check_row(b10x10 + (r + 2) * stride);

                        Vec4fb mask_lo{(below0_lo & below1_lo & below2_lo) | (above0_lo & above1_lo & above2_lo)};
                        Vec4fb mask_hi{(below0_hi & below1_hi & below2_hi) | (above0_hi & above1_hi & above2_hi)};

                        Vec4f dst_lo;
                        dst_lo.load(b8x8 + r * stride);
                        Vec4f dst_hi;
                        dst_hi.load(b8x8 + r * stride + 4);

                        select(mask_lo, f_lo_arr[r], dst_lo).store(b8x8 + r * stride);
                        select(mask_hi, f_hi_arr[r], dst_hi).store(b8x8 + r * stride + 4);
                    }
                }
                else if constexpr (sizeof(T) == 1)
                {
                    Vec16uc min_v;
                    min_v.loadl(b8x8);
                    Vec16uc max_v{min_v};
                    for (int r{1}; r < 8; ++r)
                    {
                        Vec16uc p;
                        p.loadl(b8x8 + r * stride);
                        min_v = min(min_v, p);
                        max_v = max(max_v, p);
                    }

                    uint8_t min_arr[16], max_arr[16];
                    min_v.store(min_arr);
                    max_v.store(max_arr);
                    uint8_t min_val = 255, max_val = 0;
                    for (int i{0}; i < 8; ++i)
                    {
                        min_val = std::min(min_val, min_arr[i]);
                        max_val = std::max(max_val, max_arr[i]);
                    }

                    int thr{(static_cast<int>(max_val) + static_cast<int>(min_val) + 1) >> 1};
                    int max_diff{QP >> 1};
                    Vec16uc vthr(thr), vmax_diff(max_diff);

                    auto h_avg_row = [&](const T* p) {
                        Vec16uc p0;
                        p0.loadl(p + 0);
                        Vec16uc p1;
                        p1.loadl(p + 1);
                        Vec16uc p2;
                        p2.loadl(p + 2);
                        return avg(avg(p0, p2), p1);
                    };

                    auto check_row = [&](const T* p) {
                        Vec16uc p0;
                        p0.loadl(p + 0);
                        Vec16uc p1;
                        p1.loadl(p + 1);
                        Vec16uc p2;
                        p2.loadl(p + 2);
                        Vec16cb below{(p0 < vthr) & (p1 < vthr) & (p2 < vthr)};
                        Vec16cb above = (p0 >= vthr) & (p1 >= vthr) & (p2 >= vthr);
                        return std::make_pair(below, above);
                    };

                    Vec16uc f_arr[8];

                    for (int r{0}; r < 8; ++r)
                    {
                        Vec16uc h0{h_avg_row(b10x10 + (r + 0) * stride)};
                        Vec16uc h1{h_avg_row(b10x10 + (r + 1) * stride)};
                        Vec16uc h2{h_avg_row(b10x10 + (r + 2) * stride)};

                        Vec16uc f{avg(avg(h0, h2), h1)};
                        Vec16uc orig;
                        orig.loadl(b10x10 + (r + 1) * stride + 1);

                        Vec16uc lower_bound{sub_saturated(orig, vmax_diff)};
                        Vec16uc upper_bound{add_saturated(orig, vmax_diff)};
                        f_arr[r] = max(lower_bound, min(f, upper_bound));
                    }

                    for (int r{0}; r < 8; ++r)
                    {
                        auto [below0, above0] = check_row(b10x10 + (r + 0) * stride);
                        auto [below1, above1] = check_row(b10x10 + (r + 1) * stride);
                        auto [below2, above2] = check_row(b10x10 + (r + 2) * stride);

                        Vec16cb mask{(below0 & below1 & below2) | (above0 & above1 & above2)};

                        Vec16uc dst;
                        dst.loadl(b8x8 + r * stride);
                        Vec16uc res{select(mask, f_arr[r], dst)};
                        Vec4i(res).storel(b8x8 + r * stride);
                    }
                }
                else if constexpr (sizeof(T) == 2)
                {
                    Vec8us min_v;
                    min_v.load(b8x8);
                    Vec8us max_v{min_v};
                    for (int r{1}; r < 8; ++r)
                    {
                        Vec8us p;
                        p.load(b8x8 + r * stride);
                        min_v = min(min_v, p);
                        max_v = max(max_v, p);
                    }

                    uint16_t min_arr[8], max_arr[8];
                    min_v.store(min_arr);
                    max_v.store(max_arr);
                    uint16_t min_val = 65535, max_val = 0;
                    for (int i{0}; i < 8; ++i)
                    {
                        min_val = std::min(min_val, min_arr[i]);
                        max_val = std::max(max_val, max_arr[i]);
                    }

                    int thr{(static_cast<int>(max_val) + static_cast<int>(min_val) + 1) >> 1};
                    int max_diff{(QP * (1 << (bit_depth - 8))) >> 1};
                    Vec8us vthr(thr), vmax_diff(max_diff);

                    auto h_avg_row = [&](const T* p) {
                        Vec8us p0;
                        p0.load(p + 0);
                        Vec8us p1;
                        p1.load(p + 1);
                        Vec8us p2;
                        p2.load(p + 2);
                        return avg(avg(p0, p2), p1);
                    };

                    auto check_row = [&](const T* p) {
                        Vec8us p0;
                        p0.load(p + 0);
                        Vec8us p1;
                        p1.load(p + 1);
                        Vec8us p2;
                        p2.load(p + 2);
                        Vec8sb below{(p0 < vthr) & (p1 < vthr) & (p2 < vthr)};
                        Vec8sb above = (p0 >= vthr) & (p1 >= vthr) & (p2 >= vthr);
                        return std::make_pair(below, above);
                    };

                    Vec8us f_arr[8];

                    for (int r{0}; r < 8; ++r)
                    {
                        Vec8us h0{h_avg_row(b10x10 + (r + 0) * stride)};
                        Vec8us h1{h_avg_row(b10x10 + (r + 1) * stride)};
                        Vec8us h2{h_avg_row(b10x10 + (r + 2) * stride)};

                        Vec8us f{avg(avg(h0, h2), h1)};
                        Vec8us orig;
                        orig.load(b10x10 + (r + 1) * stride + 1);

                        Vec8us lower_bound{sub_saturated(orig, vmax_diff)};
                        Vec8us upper_bound{add_saturated(orig, vmax_diff)};
                        f_arr[r] = max(lower_bound, min(f, upper_bound));
                    }

                    for (int r{0}; r < 8; ++r)
                    {
                        auto [below0, above0] = check_row(b10x10 + (r + 0) * stride);
                        auto [below1, above1] = check_row(b10x10 + (r + 1) * stride);
                        auto [below2, above2] = check_row(b10x10 + (r + 2) * stride);

                        Vec8sb mask{(below0 & below1 & below2) | (above0 & above1 & above2)};

                        Vec8us dst;
                        dst.load(b8x8 + r * stride);
                        select(mask, f_arr[r], dst).store(b8x8 + r * stride);
                    }
                }
            }
        }
    }
} // namespace SIMD_NAMESPACE
