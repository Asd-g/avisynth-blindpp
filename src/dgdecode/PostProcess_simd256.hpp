#include <tuple>

#include "PostProcess_common.h"
#include "vectorclass.h"

namespace SIMD_NAMESPACE
{
    template<typename T>
    [[nodiscard]] AVS_FORCEINLINE int deblock_vert_useDC(const T* AVS_RESTRICT v, int stride, int moderate_v, int bit_depth) noexcept
    {
        if constexpr (std::is_floating_point_v<T>)
        {
            const float thr{1.0f / 255.0f};
            Vec8f sum(0.0f);
            Vec8f ones(1.0f);
            Vec8f zeros(0.0f);

            for (int i{1}; i < 8; ++i)
            {
                Vec8f v1;
                v1.load(v + i * stride);
                Vec8f v2;
                v2.load(v + (i + 1) * stride);
                sum += select(abs(v1 - v2) <= thr, ones, zeros);
            }

            return horizontal_add(sum) > static_cast<float>(moderate_v);
        }
        else
        {
            const int thr{1 << (bit_depth - 8)};
            Vec8i sum(0);
            Vec8i ones(1);
            Vec8i zeros(0);

            auto load8i = [](const T* p) {
                Vec8i res;
                if constexpr (sizeof(T) == 1)
                    res.load_8uc(p);
                else
                    res.load_8us(p);
                return res;
            };

            for (int i{1}; i < 8; ++i)
            {
                Vec8i v1{load8i(v + i * stride)};
                Vec8i v2{load8i(v + (i + 1) * stride)};
                sum += select(abs(v1 - v2) <= thr, ones, zeros);
            }

            return horizontal_add(sum) > moderate_v;
        }
    }

    template<typename T, typename T_intermediate>
    AVS_FORCEINLINE void deblock_vert_lpf9(const T_intermediate* AVS_RESTRICT v_local, const T_intermediate* AVS_RESTRICT p1p2,
        T* AVS_RESTRICT v, int stride, int bit_depth, bool is_float_chroma) noexcept
    {
        if constexpr (std::is_floating_point_v<T>)
        {
            Vec8f vmin(is_float_chroma ? -0.5f : 0.0f);
            Vec8f vmax(is_float_chroma ? 0.5f : 1.0f);
            const float inv16{1.0f / 16.0f};

            Vec8f p1;
            p1.load(p1p2);
            Vec8f p2;
            p2.load(p1p2 + 8);

            Vec8f v1;
            v1.load(v_local + 0 * 8);
            Vec8f v2;
            v2.load(v_local + 1 * 8);
            Vec8f v3;
            v3.load(v_local + 2 * 8);
            Vec8f v4;
            v4.load(v_local + 3 * 8);
            Vec8f v5;
            v5.load(v_local + 4 * 8);
            Vec8f v6;
            v6.load(v_local + 5 * 8);
            Vec8f v7;
            v7.load(v_local + 6 * 8);
            Vec8f v8;
            v8.load(v_local + 7 * 8);

            Vec8f psum{p1 * 3.0f + v1 + v2 + v3 + v4};
            Vec8f s;

            s = ((psum + v1) * 2.0f) - (v4 - v5);
            max(vmin, min(s * inv16, vmax)).store(v + 1 * stride);

            psum += v5 - p1;
            s = ((psum + v2) * 2.0f) - (v5 - v6);
            max(vmin, min(s * inv16, vmax)).store(v + 2 * stride);

            psum += v6 - p1;
            s = ((psum + v3) * 2.0f) - (v6 - v7);
            max(vmin, min(s * inv16, vmax)).store(v + 3 * stride);

            psum += v7 - p1;
            s = ((psum + v4) * 2.0f) - (v7 - v8) + p1 - v1;
            max(vmin, min(s * inv16, vmax)).store(v + 4 * stride);

            psum += v8 - v1;
            s = ((psum + v5) * 2.0f) - (v8 - p2) + v1 - v2;
            max(vmin, min(s * inv16, vmax)).store(v + 5 * stride);

            psum += p2 - v2;
            s = ((psum + v6) * 2.0f) + v2 - v3;
            max(vmin, min(s * inv16, vmax)).store(v + 6 * stride);

            psum += p2 - v3;
            s = ((psum + v7) * 2.0f) + v3 - v4;
            max(vmin, min(s * inv16, vmax)).store(v + 7 * stride);

            psum += p2 - v4;
            s = ((psum + v8) * 2.0f) + v4 - v5;
            max(vmin, min(s * inv16, vmax)).store(v + 8 * stride);
        }
        else
        {
            Vec8i vmin(0);
            Vec8i vmax((1 << bit_depth) - 1);

            auto load8i = [](const T_intermediate* p) {
                if constexpr (sizeof(T_intermediate) == 2)
                {
                    Vec8i res;
                    res.load_8us(p);
                    return res;
                }
                else
                {
                    Vec8i res;
                    res.load(p);
                    return res;
                }
            };

            auto store8i = [](Vec8i vec, T* p) {
                if constexpr (sizeof(T) == 1)
                {
                    Vec16us v16{compress_saturated_s2u(vec, vec)};
                    Vec32uc v8{compress_saturated_s2u(Vec16s(v16), Vec16s(v16))};
                    Vec16uc(v8.get_low()).storel(p);
                }
                else if constexpr (sizeof(T) == 2)
                {
                    Vec16us v16{compress_saturated_s2u(vec, vec)};
                    v16.get_low().store(p);
                }
            };

            Vec8i p1{load8i(p1p2)};
            Vec8i p2{load8i(p1p2 + 8)};

            Vec8i v1{load8i(v_local + 0 * 8)};
            Vec8i v2{load8i(v_local + 1 * 8)};
            Vec8i v3{load8i(v_local + 2 * 8)};
            Vec8i v4{load8i(v_local + 3 * 8)};
            Vec8i v5{load8i(v_local + 4 * 8)};
            Vec8i v6{load8i(v_local + 5 * 8)};
            Vec8i v7{load8i(v_local + 6 * 8)};
            Vec8i v8{load8i(v_local + 7 * 8)};

            Vec8i psum{p1 * 3 + v1 + v2 + v3 + v4 + 4};
            Vec8i s;

            s = ((psum + v1) << 1) - (v4 - v5);
            store8i(max(vmin, min(s >> 4, vmax)), v + 1 * stride);

            psum += v5 - p1;
            s = ((psum + v2) << 1) - (v5 - v6);
            store8i(max(vmin, min(s >> 4, vmax)), v + 2 * stride);

            psum += v6 - p1;
            s = ((psum + v3) << 1) - (v6 - v7);
            store8i(max(vmin, min(s >> 4, vmax)), v + 3 * stride);

            psum += v7 - p1;
            s = ((psum + v4) << 1) - (v7 - v8) + p1 - v1;
            store8i(max(vmin, min(s >> 4, vmax)), v + 4 * stride);

            psum += v8 - v1;
            s = ((psum + v5) << 1) - (v8 - p2) + v1 - v2;
            store8i(max(vmin, min(s >> 4, vmax)), v + 5 * stride);

            psum += p2 - v2;
            s = ((psum + v6) << 1) + v2 - v3;
            store8i(max(vmin, min(s >> 4, vmax)), v + 6 * stride);

            psum += p2 - v3;
            s = ((psum + v7) << 1) + v3 - v4;
            store8i(max(vmin, min(s >> 4, vmax)), v + 7 * stride);

            psum += p2 - v4;
            s = ((psum + v8) << 1) + v4 - v5;
            store8i(max(vmin, min(s >> 4, vmax)), v + 8 * stride);
        }
    }

    template<typename T>
    AVS_FORCEINLINE void deblock_vert_default_filter(T* AVS_RESTRICT v, int stride, int QP, int bit_depth, bool is_float_chroma) noexcept
    {
        if constexpr (std::is_floating_point_v<T>)
        {
            const float qp_thresh{static_cast<float>(QP) / 255.0f};
            Vec8f vmin(is_float_chroma ? -0.5f : 0.0f);
            Vec8f vmax(is_float_chroma ? 0.5f : 1.0f);

            Vec8f v1;
            v1.load(v + 1 * stride);
            Vec8f v2;
            v2.load(v + 2 * stride);
            Vec8f v3;
            v3.load(v + 3 * stride);
            Vec8f v4;
            v4.load(v + 4 * stride);
            Vec8f v5;
            v5.load(v + 5 * stride);
            Vec8f v6;
            v6.load(v + 6 * stride);
            Vec8f v7;
            v7.load(v + 7 * stride);
            Vec8f v8;
            v8.load(v + 8 * stride);

            Vec8f q1{v4 - v5};
            Vec8f q{q1 * 0.5f};
            Vec8f a3_0{(v3 - v6) * 2.0f - q1 * 5.0f};
            Vec8f aa3_0{abs(a3_0)};

            Vec8fb cond1{aa3_0 < (8.0f * qp_thresh)};
            Vec8f a3_1{abs((v3 - v2) * 5.0f + (v1 - v4) * 2.0f)};
            Vec8f a3_2{abs((v7 - v6) * 5.0f + (v5 - v8) * 2.0f)};
            Vec8f d_abs{aa3_0 - min(a3_1, a3_2)};

            Vec8fb cond2{d_abs > 0.0f};
            Vec8f d{(d_abs * 5.0f + (32.0f / 255.0f)) * (1.0f / 64.0f)};

            Vec8fb cond3{(q > 0.0f & a3_0 < 0.0f) | (q < 0.0f & a3_0 >= 0.0f)};
            Vec8f delta_mag{min(d, abs(q))};
            Vec8f final_delta{select(cond3, select(q > 0.0f, delta_mag, -delta_mag), 0.0f)};

            final_delta = select((q != 0.0f) & cond1 & cond2, final_delta, 0.0f);

            max(vmin, min(v4 - final_delta, vmax)).store(v + 4 * stride);
            max(vmin, min(v5 + final_delta, vmax)).store(v + 5 * stride);
        }
        else
        {
            const int scale{1 << (bit_depth - 8)};
            const int qp_scaled{QP * scale};
            Vec8i vmin(0);
            Vec8i vmax((1 << bit_depth) - 1);

            auto load8i = [](const T* p) {
                Vec8i res;
                if constexpr (sizeof(T) == 1)
                    res.load_8uc(p);
                else
                    res.load_8us(p);
                return res;
            };

            auto store8i = [](Vec8i vec, T* p) {
                if constexpr (sizeof(T) == 1)
                {
                    Vec16us v16{compress_saturated_s2u(vec, vec)};
                    Vec32uc v8{compress_saturated_s2u(Vec16s(v16), Vec16s(v16))};
                    Vec16uc(v8.get_low()).storel(p);
                }
                else if constexpr (sizeof(T) == 2)
                {
                    Vec16us v16{compress_saturated_s2u(vec, vec)};
                    v16.get_low().store(p);
                }
            };

            Vec8i v1{load8i(v + 1 * stride)};
            Vec8i v2{load8i(v + 2 * stride)};
            Vec8i v3{load8i(v + 3 * stride)};
            Vec8i v4{load8i(v + 4 * stride)};
            Vec8i v5{load8i(v + 5 * stride)};
            Vec8i v6{load8i(v + 6 * stride)};
            Vec8i v7{load8i(v + 7 * stride)};
            Vec8i v8{load8i(v + 8 * stride)};

            Vec8i q1{v4 - v5};
            Vec8i q{q1 / 2};
            Vec8i a3_0{(v3 - v6) * 2 - q1 * 5};
            Vec8i aa3_0{abs(a3_0)};

            Vec8ib cond1{aa3_0 < (8 * qp_scaled)};
            Vec8i a3_1{abs((v3 - v2) * 5 + (v1 - v4) * 2)};
            Vec8i a3_2{abs((v7 - v6) * 5 + (v5 - v8) * 2)};
            Vec8i d_abs{aa3_0 - min(a3_1, a3_2)};

            Vec8ib cond2{d_abs > 0};
            Vec8i d{(d_abs * 5 + (32 * scale)) >> 6};

            Vec8ib cond3{(q > 0 & a3_0 < 0) | (q < 0 & a3_0 >= 0)};
            Vec8i delta_mag{min(d, abs(q))};
            Vec8i final_delta{select(cond3, select(q > 0, delta_mag, -delta_mag), 0)};

            final_delta = select((q != 0) & cond1 & cond2, final_delta, 0);

            store8i(max(vmin, min(v4 - final_delta, vmax)), v + 4 * stride);
            store8i(max(vmin, min(v5 + final_delta, vmax)), v + 5 * stride);
        }
    }

    template<typename T>
    void deblock_vert(T* AVS_RESTRICT image, int width, int stride, const QP_STORE_T* AVS_RESTRICT QP_store, int QP_stride,
        int chroma_flag_for_qp, int moderate_v, int bit_depth, bool is_float_chroma) noexcept
    {
        using intermediate_t = typename DeblockVertTypes<T>::intermediate;
        alignas(32) std::array<intermediate_t, 8 * 8> v_local;
        alignas(32) std::array<intermediate_t, 2 * 8> p1p2;

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
                    Vec8f vmin_clamp(is_float_chroma ? -0.5f : 0.0f);
                    Vec8f vmax_clamp(is_float_chroma ? 0.5f : 1.0f);

                    Vec8f min_v;
                    min_v.load(b8x8);
                    Vec8f max_v{min_v};
                    for (int r{1}; r < 8; ++r)
                    {
                        Vec8f p;
                        p.load(b8x8 + r * stride);
                        min_v = min(min_v, p);
                        max_v = max(max_v, p);
                    }
                    Vec8f vthr((horizontal_max(max_v) + horizontal_min(min_v)) * 0.5f);
                    Vec8f vmax_diff((static_cast<float>(QP) / 255.0f) * 0.5f);

                    auto h_avg_row = [](const T* p) {
                        Vec8f p0;
                        p0.load(p + 0);
                        Vec8f p1;
                        p1.load(p + 1);
                        Vec8f p2;
                        p2.load(p + 2);
                        return ((p0 + p2) * 0.5f + p1) * 0.5f;
                    };

                    auto check_row = [&](const T* p) {
                        Vec8f p0;
                        p0.load(p + 0);
                        Vec8f p1;
                        p1.load(p + 1);
                        Vec8f p2;
                        p2.load(p + 2);
                        Vec8fb below{(p0 < vthr) & (p1 < vthr) & (p2 < vthr)};
                        Vec8fb above{(p0 >= vthr) & (p1 >= vthr) & (p2 >= vthr)};
                        return std::make_pair(below, above);
                    };

                    Vec8f f_arr[8];

                    for (int r{0}; r < 8; ++r)
                    {
                        Vec8f h0{h_avg_row(b10x10 + (r + 0) * stride)};
                        Vec8f h1{h_avg_row(b10x10 + (r + 1) * stride)};
                        Vec8f h2{h_avg_row(b10x10 + (r + 2) * stride)};

                        Vec8f f{((h0 + h2) * 0.5f + h1) * 0.5f};
                        Vec8f orig;
                        orig.load(b10x10 + (r + 1) * stride + 1);

                        f_arr[r] = max(vmin_clamp, min(max(orig - vmax_diff, min(f, orig + vmax_diff)), vmax_clamp));
                    }

                    for (int r{0}; r < 8; ++r)
                    {
                        auto [below0, above0] = check_row(b10x10 + (r + 0) * stride);
                        auto [below1, above1] = check_row(b10x10 + (r + 1) * stride);
                        auto [below2, above2] = check_row(b10x10 + (r + 2) * stride);

                        Vec8fb mask{(below0 & below1 & below2) | (above0 & above1 & above2)};

                        Vec8f dst;
                        dst.load(b8x8 + r * stride);
                        select(mask, f_arr[r], dst).store(b8x8 + r * stride);
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

                    uint8_t min_arr[16];
                    uint8_t max_arr[16];
                    min_v.store(min_arr);
                    max_v.store(max_arr);
                    uint8_t min_val{255};
                    uint8_t max_val{0};
                    for (int i{0}; i < 8; ++i)
                    {
                        min_val = std::min(min_val, min_arr[i]);
                        max_val = std::max(max_val, max_arr[i]);
                    }

                    Vec16uc vthr((static_cast<int>(max_val) + static_cast<int>(min_val) + 1) >> 1);
                    Vec16uc vmax_diff(QP >> 1);

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
                        res.storel(b8x8 + r * stride);
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

                    uint16_t min_arr[8];
                    uint16_t max_arr[8];
                    min_v.store(min_arr);
                    max_v.store(max_arr);
                    uint16_t min_val{65535};
                    uint16_t max_val{0};
                    for (int i{0}; i < 8; ++i)
                    {
                        min_val = std::min(min_val, min_arr[i]);
                        max_val = std::max(max_val, max_arr[i]);
                    }

                    Vec8us vthr((static_cast<int>(max_val) + static_cast<int>(min_val) + 1) >> 1);
                    Vec8us vmax_diff((QP * (1 << (bit_depth - 8))) >> 1);

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
                        Vec8sb above{(p0 >= vthr) & (p1 >= vthr) & (p2 >= vthr)};
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
