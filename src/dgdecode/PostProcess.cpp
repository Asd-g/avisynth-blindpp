#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <type_traits>

#include "PostProcess.h"

namespace
{
    template<typename T>
    struct DeblockVertTypes;
    template<>
    struct DeblockVertTypes<uint8_t>
    {
        using intermediate = uint16_t;
    };
    template<>
    struct DeblockVertTypes<uint16_t>
    {
        using intermediate = uint32_t;
    };
    template<>
    struct DeblockVertTypes<float>
    {
        using intermediate = float;
    };

    template<typename T>
    [[nodiscard]] AVS_FORCEINLINE int deblock_horiz_useDC(const T* AVS_RESTRICT v, int stride, int moderate_h, int bit_depth) noexcept
    {
        int eq_cnt{};

        if constexpr (std::is_floating_point_v<T>)
        {
            const T thr{static_cast<T>(1.0 / 255.0)};

            for (int y{0}; y < 4; ++y)
            {
                const T* p{v + y * stride};

                for (int x{1}; x < 8; ++x)
                {
                    if (std::abs(p[x] - p[x + 1]) <= thr)
                        eq_cnt++;
                }
            }
        }
        else
        {
            const int thr{1 << (bit_depth - 8)};

            for (int y{0}; y < 4; ++y)
            {
                const T* p{v + y * stride};

                for (int x{1}; x < 8; ++x)
                {
                    if (std::abs(static_cast<int>(p[x]) - static_cast<int>(p[x + 1])) <= thr)
                        eq_cnt++;
                }
            }
        }

        return eq_cnt >= moderate_h;
    }

    template<typename T>
    [[nodiscard]] AVS_FORCEINLINE bool deblock_horiz_DC_on(const T* AVS_RESTRICT v, int stride, int QP, int bit_depth) noexcept
    {
        if constexpr (std::is_floating_point_v<T>)
        {
            const T qp_thresh{static_cast<T>(QP * 2) / 255.0f};

            for (int i{0}; i < 4; ++i)
            {
                if (std::abs(v[0] - v[5]) >= qp_thresh || std::abs(v[1] - v[8]) >= qp_thresh || std::abs(v[1] - v[4]) >= qp_thresh ||
                    std::abs(v[2] - v[7]) >= qp_thresh || std::abs(v[3] - v[6]) >= qp_thresh)
                    return false;

                v += stride;
            }
        }
        else
        {
            const int qp_thresh{QP * 2 * (1 << (bit_depth - 8))};

            for (int i{0}; i < 4; ++i)
            {
                if (std::abs(static_cast<int>(v[0]) - static_cast<int>(v[5])) >= qp_thresh ||
                    std::abs(static_cast<int>(v[1]) - static_cast<int>(v[8])) >= qp_thresh ||
                    std::abs(static_cast<int>(v[1]) - static_cast<int>(v[4])) >= qp_thresh ||
                    std::abs(static_cast<int>(v[2]) - static_cast<int>(v[7])) >= qp_thresh ||
                    std::abs(static_cast<int>(v[3]) - static_cast<int>(v[6])) >= qp_thresh)
                    return false;

                v += stride;
            }
        }

        return true;
    }

    template<typename T>
    AVS_FORCEINLINE void deblock_horiz_lpf9(T* AVS_RESTRICT v, int stride, int QP, int bit_depth, bool is_float_chroma) noexcept
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

        if constexpr (std::is_floating_point_v<T>)
        {
            const T qp_thresh{static_cast<T>(QP) / 255.0f};
            const T min_val{is_float_chroma ? static_cast<T>(-0.5) : static_cast<T>(0.0)};
            const T max_val{is_float_chroma ? static_cast<T>(0.5) : static_cast<T>(1.0)};

            for (int y{0}; y < 4; ++y)
            {
                const T p1{(std::abs(v[0] - v[1]) < qp_thresh) ? v[0] : v[1]};
                const T p2{(std::abs(v[8] - v[9]) < qp_thresh) ? v[9] : v[8]};

                const T in[10]{p1, v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8], p2};
                T out[8];

                for (int i{0}; i < 8; ++i)
                {
                    T sum{};

                    for (int j{0}; j < 10; ++j)
                        sum += in[j] * static_cast<T>(C[i][j]);

                    out[i] = std::clamp(sum / static_cast<T>(16.0), min_val, max_val);
                }

                for (int i{0}; i < 8; ++i)
                    v[i + 1] = out[i];

                v += stride;
            }
        }
        else
        {
            const int qp_scaled{QP * (1 << (bit_depth - 8))};
            const int max_pixel_val{(1 << bit_depth) - 1};
            using intermediate_t = int32_t;

            for (int y{0}; y < 4; ++y)
            {
                const intermediate_t p1{(std::abs(static_cast<int>(v[0]) - static_cast<int>(v[1])) < qp_scaled) ? v[0] : v[1]};
                const intermediate_t p2{(std::abs(static_cast<int>(v[8]) - static_cast<int>(v[9])) < qp_scaled) ? v[9] : v[8]};

                const intermediate_t in[10]{p1, v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8], p2};
                T out[8];

                for (int i{0}; i < 8; ++i)
                {
                    intermediate_t sum{8}; // for rounding

                    for (int j{0}; j < 10; ++j)
                        sum += in[j] * C[i][j];

                    out[i] = static_cast<T>(std::clamp(sum >> 4, 0, max_pixel_val));
                }

                for (int i{0}; i < 8; ++i)
                    v[i + 1] = out[i];

                v += stride;
            }
        }
    }

    template<typename T>
    AVS_FORCEINLINE void deblock_horiz_default_filter(T* AVS_RESTRICT v, int stride, int QP, int bit_depth, bool is_float_chroma) noexcept
    {
        if constexpr (std::is_floating_point_v<T>)
        {
            const T qp_thresh{static_cast<T>(QP) / 255.0f};
            const T min_val{is_float_chroma ? static_cast<T>(-0.5) : static_cast<T>(0.0)};
            const T max_val{is_float_chroma ? static_cast<T>(0.5) : static_cast<T>(1.0)};

            for (int y{0}; y < 4; ++y)
            {
                const T q1{v[4] - v[5]};
                const T q{q1 * static_cast<T>(0.5)};

                if (q != 0)
                {
                    const T a3_0{static_cast<T>(2.0) * (v[3] - v[6]) - static_cast<T>(5.0) * q1};
                    const T aa3_0{std::abs(a3_0)};

                    if (aa3_0 < static_cast<T>(8.0) * qp_thresh)
                    {
                        const T a3_1{std::abs(static_cast<T>(5.0) * (v[3] - v[2]) + static_cast<T>(2.0) * (v[1] - v[4]))};
                        const T a3_2{std::abs(static_cast<T>(5.0) * (v[7] - v[8]) + static_cast<T>(2.0) * (v[5] - v[8]))};
                        const T d{aa3_0 - std::min(a3_1, a3_2)};

                        if (d > static_cast<T>(6.0 / 255.0))
                        {
                            const T delta_mag{std::min(d, std::abs(q))};
                            T delta{0};

                            if ((q > 0 && a3_0 < 0) || (q < 0 && a3_0 > 0))
                                delta = (q > 0) ? -delta_mag : delta_mag;

                            if (delta != 0)
                            {
                                v[4] = std::clamp(v[4] + delta, min_val, max_val);
                                v[5] = std::clamp(v[5] - delta, min_val, max_val);
                            }
                        }
                    }
                }
                v += stride;
            }
        }
        else
        {
            const int scale{1 << (bit_depth - 8)};
            const int qp_scaled{QP * scale};
            const int max_pixel_val{(1 << bit_depth) - 1};
            for (int y{0}; y < 4; ++y)
            {
                const int q1{static_cast<int>(v[4]) - static_cast<int>(v[5])};
                const int q{q1 / 2};

                if (q != 0)
                {
                    const int a3_0{2 * (static_cast<int>(v[3]) - static_cast<int>(v[6])) - 5 * q1};
                    const int aa3_0{std::abs(a3_0)};

                    if (aa3_0 < 8 * qp_scaled)
                    {
                        const int a3_1 = std::abs(
                            5 * (static_cast<int>(v[3]) - static_cast<int>(v[2])) + 2 * (static_cast<int>(v[1]) - static_cast<int>(v[4])));
                        const int a3_2 = std::abs(
                            5 * (static_cast<int>(v[7]) - static_cast<int>(v[8])) + 2 * (static_cast<int>(v[5]) - static_cast<int>(v[8])));
                        const int d{aa3_0 - std::min(a3_1, a3_2)};

                        if (d > (6 * scale))
                        {
                            const int delta_mag{std::min(d, std::abs(q))};
                            int delta{0};

                            if ((q > 0 && a3_0 < 0) || (q < 0 && a3_0 > 0))
                                delta = (q > 0) ? -delta_mag : delta_mag;

                            if (delta != 0)
                            {
                                v[4] = static_cast<T>(std::clamp(static_cast<int>(v[4]) + delta, 0, max_pixel_val));
                                v[5] = static_cast<T>(std::clamp(static_cast<int>(v[5]) - delta, 0, max_pixel_val));
                            }
                        }
                    }
                }
                v += stride;
            }
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
        int eq_cnt{};

        if constexpr (std::is_floating_point_v<T>)
        {
            const T thr{static_cast<T>(1.0 / 255.0)};

            for (int i{1}; i < 8; ++i)
            {
                const T* p1{v + i * stride};
                const T* p2{v + (i + 1) * stride};

                for (int j{0}; j < 8; ++j)
                {
                    if (std::abs(p1[j] - p2[j]) <= thr)
                        eq_cnt++;
                }
            }
        }
        else
        {
            const int thr{1 << (bit_depth - 8)};

            for (int i{1}; i < 8; ++i)
            {
                const T* p1{v + i * stride};
                const T* p2{v + (i + 1) * stride};

                for (int j{0}; j < 8; ++j)
                {
                    if (std::abs(static_cast<int>(p1[j]) - static_cast<int>(p2[j])) <= thr)
                        eq_cnt++;
                }
            }
        }

        return eq_cnt > moderate_v;
    }

    template<typename T>
    [[nodiscard]] AVS_FORCEINLINE bool deblock_vert_DC_on(const T* AVS_RESTRICT v, int stride, int QP, int bit_depth) noexcept
    {
        if constexpr (std::is_floating_point_v<T>)
        {
            const T qp_thresh{static_cast<T>(QP * 2) / 255.0f};

            for (int i{0}; i < 8; ++i)
            {
                if (std::abs(v[i + 0 * stride] - v[i + 5 * stride]) >= qp_thresh ||
                    std::abs(v[i + 1 * stride] - v[i + 4 * stride]) >= qp_thresh ||
                    std::abs(v[i + 1 * stride] - v[i + 8 * stride]) >= qp_thresh ||
                    std::abs(v[i + 2 * stride] - v[i + 7 * stride]) >= qp_thresh ||
                    std::abs(v[i + 3 * stride] - v[i + 6 * stride]) >= qp_thresh)
                    return false;
            }
        }
        else
        {
            const int qp_thresh{QP * 2 * (1 << (bit_depth - 8))};

            for (int i{0}; i < 8; ++i)
            {
                if (std::abs(static_cast<int>(v[i + 0 * stride]) - static_cast<int>(v[i + 5 * stride])) >= qp_thresh ||
                    std::abs(static_cast<int>(v[i + 1 * stride]) - static_cast<int>(v[i + 4 * stride])) >= qp_thresh ||
                    std::abs(static_cast<int>(v[i + 1 * stride]) - static_cast<int>(v[i + 8 * stride])) >= qp_thresh ||
                    std::abs(static_cast<int>(v[i + 2 * stride]) - static_cast<int>(v[i + 7 * stride])) >= qp_thresh ||
                    std::abs(static_cast<int>(v[i + 3 * stride]) - static_cast<int>(v[i + 6 * stride])) >= qp_thresh)
                    return false;
            }
        }

        return true;
    }

    template<typename T, typename T_intermediate>
    AVS_FORCEINLINE void deblock_vert_copy_and_unpack(
        int stride, const T* AVS_RESTRICT source, T_intermediate* AVS_RESTRICT dest, int n) noexcept
    {
        for (int i{0}; i < n; ++i)
        {
            const T* s_row{source + i * stride};
            T_intermediate* d_row{dest + i * 8};

            for (int j{0}; j < 8; ++j)
                d_row[j] = s_row[j];
        }
    }

    template<typename T, typename T_intermediate>
    AVS_FORCEINLINE void deblock_vert_choose_p1p2(
        const T* AVS_RESTRICT v, int stride, T_intermediate* AVS_RESTRICT p1p2, int QP, int bit_depth) noexcept
    {
        T_intermediate* p1_ptr{p1p2};
        T_intermediate* p2_ptr{p1p2 + 8};

        const T* v0{v};
        const T* v1{v + stride};
        const T* v8{v + 8 * stride};
        const T* v9{v + 9 * stride};

        if constexpr (std::is_floating_point_v<T>)
        {
            const T qp_thresh{static_cast<T>(QP) / 255.0f};

            for (int i{0}; i < 8; ++i)
            {
                p1_ptr[i] = (std::abs(v0[i] - v1[i]) <= qp_thresh) ? v0[i] : v1[i];
                p2_ptr[i] = (std::abs(v8[i] - v9[i]) <= qp_thresh) ? v9[i] : v8[i];
            }
        }
        else
        {
            const int qp_scaled{QP * (1 << (bit_depth - 8))};

            for (int i{0}; i < 8; ++i)
            {
                p1_ptr[i] = (std::abs(static_cast<int>(v0[i]) - static_cast<int>(v1[i])) <= qp_scaled) ? v0[i] : v1[i];
                p2_ptr[i] = (std::abs(static_cast<int>(v8[i]) - static_cast<int>(v9[i])) <= qp_scaled) ? v9[i] : v8[i];
            }
        }
    }

    template<typename T, typename T_intermediate>
    AVS_FORCEINLINE void deblock_vert_lpf9(const T_intermediate* AVS_RESTRICT v_local, const T_intermediate* AVS_RESTRICT p1p2,
        T* AVS_RESTRICT v, int stride, int bit_depth, bool is_float_chroma) noexcept
    {
        const T_intermediate* vv{v_local};
        const T_intermediate* p1p2_16{p1p2};

        if constexpr (std::is_floating_point_v<T>)
        {
            const T min_val{is_float_chroma ? static_cast<T>(-0.5) : static_cast<T>(0.0)};
            const T max_val{is_float_chroma ? static_cast<T>(0.5) : static_cast<T>(1.0)};

            for (int i{0}; i < 8; ++i)
            {
                const T p1{p1p2_16[i]};
                const T p2{p1p2_16[i + 8]};

                const T v1{vv[0 * 8 + i]};
                const T v2{vv[1 * 8 + i]};
                const T v3{vv[2 * 8 + i]};
                const T v4{vv[3 * 8 + i]};
                const T v5{vv[4 * 8 + i]};
                const T v6{vv[5 * 8 + i]};
                const T v7{vv[6 * 8 + i]};
                const T v8{vv[7 * 8 + i]};

                T psum{p1 * 3 + v1 + v2 + v3 + v4};

                T s{((psum + v1) * 2) - (v4 - v5)};
                v[(1 * stride) + i] = std::clamp(s / 16, min_val, max_val);
                psum += v5 - p1;
                s = ((psum + v2) * 2) - (v5 - v6);
                v[(2 * stride) + i] = std::clamp(s / 16, min_val, max_val);
                psum += v6 - p1;
                s = ((psum + v3) * 2) - (v6 - v7);
                v[(3 * stride) + i] = std::clamp(s / 16, min_val, max_val);
                psum += v7 - p1;
                s = ((psum + v4) * 2) - (v7 - v8) + p1 - v1;
                v[(4 * stride) + i] = std::clamp(s / 16, min_val, max_val);
                psum += v8 - v1;
                s = ((psum + v5) * 2) - (v8 - p2) + v1 - v2;
                v[(5 * stride) + i] = std::clamp(s / 16, min_val, max_val);
                psum += p2 - v2;
                s = ((psum + v6) * 2) + v2 - v3;
                v[(6 * stride) + i] = std::clamp(s / 16, min_val, max_val);
                psum += p2 - v3;
                s = ((psum + v7) * 2) + v3 - v4;
                v[(7 * stride) + i] = std::clamp(s / 16, min_val, max_val);
                psum += p2 - v4;
                s = ((psum + v8) * 2) + v4 - v5;
                v[(8 * stride) + i] = std::clamp(s / 16, min_val, max_val);
            }
        }
        else
        {
            const int max_pixel_val{(1 << bit_depth) - 1};

            for (int i{0}; i < 8; ++i)
            {
                const int p1{static_cast<int>(p1p2_16[i])};
                const int p2{static_cast<int>(p1p2_16[i + 8])};

                const int v1{static_cast<int>(vv[0 * 8 + i])};
                const int v2{static_cast<int>(vv[1 * 8 + i])};
                const int v3{static_cast<int>(vv[2 * 8 + i])};
                const int v4{static_cast<int>(vv[3 * 8 + i])};
                const int v5{static_cast<int>(vv[4 * 8 + i])};
                const int v6{static_cast<int>(vv[5 * 8 + i])};
                const int v7{static_cast<int>(vv[6 * 8 + i])};
                const int v8{static_cast<int>(vv[7 * 8 + i])};

                int psum{4 + p1 * 3 + v1 + v2 + v3 + v4};

                int s{((psum + v1) << 1) - (v4 - v5)};
                v[(1 * stride) + i] = static_cast<T>(std::clamp(s >> 4, 0, max_pixel_val));
                psum += v5 - p1;
                s = ((psum + v2) << 1) - (v5 - v6);
                v[(2 * stride) + i] = static_cast<T>(std::clamp(s >> 4, 0, max_pixel_val));
                psum += v6 - p1;
                s = ((psum + v3) << 1) - (v6 - v7);
                v[(3 * stride) + i] = static_cast<T>(std::clamp(s >> 4, 0, max_pixel_val));
                psum += v7 - p1;
                s = ((psum + v4) << 1) - (v7 - v8) + p1 - v1;
                v[(4 * stride) + i] = static_cast<T>(std::clamp(s >> 4, 0, max_pixel_val));
                psum += v8 - v1;
                s = ((psum + v5) << 1) - (v8 - p2) + v1 - v2;
                v[(5 * stride) + i] = static_cast<T>(std::clamp(s >> 4, 0, max_pixel_val));
                psum += p2 - v2;
                s = ((psum + v6) << 1) + v2 - v3;
                v[(6 * stride) + i] = static_cast<T>(std::clamp(s >> 4, 0, max_pixel_val));
                psum += p2 - v3;
                s = ((psum + v7) << 1) + v3 - v4;
                v[(7 * stride) + i] = static_cast<T>(std::clamp(s >> 4, 0, max_pixel_val));
                psum += p2 - v4;
                s = ((psum + v8) << 1) + v4 - v5;
                v[(8 * stride) + i] = static_cast<T>(std::clamp(s >> 4, 0, max_pixel_val));
            }
        }
    }

    template<typename T>
    AVS_FORCEINLINE void deblock_vert_default_filter(T* AVS_RESTRICT v, int stride, int QP, int bit_depth, bool is_float_chroma) noexcept
    {
        if constexpr (std::is_floating_point_v<T>)
        {
            const T qp_thresh{static_cast<T>(QP) / 255.0f};
            const T min_val{is_float_chroma ? static_cast<T>(-0.5) : static_cast<T>(0.0)};
            const T max_val{is_float_chroma ? static_cast<T>(0.5) : static_cast<T>(1.0)};

            for (int i{0}; i < 8; ++i)
            {
                T* p{v + i};
                const T v1{p[1 * stride]};
                const T v2{p[2 * stride]};
                const T v3{p[3 * stride]};
                const T v4{p[4 * stride]};
                const T v5{p[5 * stride]};
                const T v6{p[6 * stride]};
                const T v7{p[7 * stride]};
                const T v8{p[8 * stride]};

                const T q1{v4 - v5};
                const T q{q1 * static_cast<T>(0.5)};

                if (q != 0)
                {
                    const T a3_0{static_cast<T>(2.0) * (v3 - v6) - static_cast<T>(5.0) * q1};
                    const T aa3_0{std::abs(a3_0)};

                    if (aa3_0 < static_cast<T>(8.0) * qp_thresh)
                    {
                        const T a3_1{std::abs(static_cast<T>(5.0) * (v3 - v2) + static_cast<T>(2.0) * (v1 - v4))};
                        const T a3_2{std::abs(static_cast<T>(5.0) * (v7 - v6) + static_cast<T>(2.0) * (v5 - v8))};
                        const T d_abs{aa3_0 - std::min(a3_1, a3_2)};

                        if (d_abs > 0)
                        {
                            const T d{(static_cast<T>(5.0) * d_abs + static_cast<T>(32.0 / 255.0)) / static_cast<T>(64.0)};
                            T final_delta{};

                            if ((q > 0 && a3_0 < 0) || (q < 0 && a3_0 >= 0))
                            {
                                const T delta_mag{std::min(d, std::abs(q))};
                                final_delta = (q > 0) ? delta_mag : -delta_mag;
                            }

                            if (final_delta != 0)
                            {
                                p[4 * stride] = std::clamp(p[4 * stride] - final_delta, min_val, max_val);
                                p[5 * stride] = std::clamp(p[5 * stride] + final_delta, min_val, max_val);
                            }
                        }
                    }
                }
            }
        }
        else
        {
            const int scale{1 << (bit_depth - 8)};
            const int qp_scaled{QP * scale};
            const int max_pixel_val{(1 << bit_depth) - 1};
            for (int i{0}; i < 8; ++i)
            {
                T* p{v + i};
                const int v1{p[1 * stride]};
                const int v2{p[2 * stride]};
                const int v3{p[3 * stride]};
                const int v4{p[4 * stride]};
                const int v5{p[5 * stride]};
                const int v6{p[6 * stride]};
                const int v7{p[7 * stride]};
                const int v8{p[8 * stride]};

                const int q1{v4 - v5};
                const int q{q1 / 2};

                if (q != 0)
                {
                    const int a3_0{2 * (v3 - v6) - 5 * q1};
                    const int aa3_0{std::abs(a3_0)};

                    if (aa3_0 < 8 * qp_scaled)
                    {
                        const int a3_1{std::abs(5 * (v3 - v2) + 2 * (v1 - v4))};
                        const int a3_2{std::abs(5 * (v7 - v6) + 2 * (v5 - v8))};
                        const int d_abs{aa3_0 - std::min(a3_1, a3_2)};

                        if (d_abs > 0)
                        {
                            const int d{(5 * d_abs + (32 * scale)) >> 6};
                            int final_delta{};

                            if ((q > 0 && a3_0 < 0) || (q < 0 && a3_0 >= 0))
                            {
                                const int delta_mag{std::min(d, std::abs(q))};
                                final_delta = (q > 0) ? delta_mag : -delta_mag;
                            }

                            if (final_delta != 0)
                            {
                                p[4 * stride] = static_cast<T>(std::clamp(static_cast<int>(p[4 * stride]) - final_delta, 0, max_pixel_val));
                                p[5 * stride] = static_cast<T>(std::clamp(static_cast<int>(p[5 * stride]) + final_delta, 0, max_pixel_val));
                            }
                        }
                    }
                }
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
        alignas(16) std::array<T, 64> b8x8filtered;

        for (int y{8}; y < height - 8; y += 8)
        {
            for (int x{8}; x < width - 8; x += 8)
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

                const QP_STORE_T QP{QP_store[(y >> qp_y_idx_shift) * QP_stride + qp_x_idx]};

                T* b8x8{image + static_cast<ptrdiff_t>(stride) * y + x};
                const T* b10x10{image + static_cast<ptrdiff_t>(stride) * (y - 1) + (x - 1)};

                T thr{};

                if constexpr (std::is_floating_point_v<T>)
                {
                    T min_val{1.0f}, max_val{is_float_chroma ? -1.0f : 0.0f};
                    const T min_clamp{is_float_chroma ? static_cast<T>(-0.5) : static_cast<T>(0.0)};
                    const T max_clamp{is_float_chroma ? static_cast<T>(0.5) : static_cast<T>(1.0)};

                    for (int r{0}; r < 8; ++r)
                    {
                        for (int c{0}; c < 8; ++c)
                        {
                            const T p{b8x8[static_cast<ptrdiff_t>(r) * stride + c]};
                            min_val = std::min(min_val, p);
                            max_val = std::max(max_val, p);
                        }
                    }

                    thr = (max_val + min_val) * static_cast<T>(0.5);
                    const T max_diff{static_cast<T>(QP) / 255.0f * static_cast<T>(0.5)};

                    auto pavgb{[](T a, T b) { return (a + b) * static_cast<T>(0.5); }};

                    for (int r{0}; r < 8; ++r)
                    {
                        for (int c{0}; c < 8; ++c)
                        {
                            const T* p_row0{&b10x10[static_cast<ptrdiff_t>(r + 0) * stride + c]};
                            const T* p_row1{&b10x10[static_cast<ptrdiff_t>(r + 1) * stride + c]};
                            const T* p_row2{&b10x10[static_cast<ptrdiff_t>(r + 2) * stride + c]};

                            const T h_avg0{pavgb(pavgb(p_row0[0], p_row0[2]), p_row0[1])};
                            const T h_avg1{pavgb(pavgb(p_row1[0], p_row1[2]), p_row1[1])};
                            const T h_avg2{pavgb(pavgb(p_row2[0], p_row2[2]), p_row2[1])};
                            T filtered_val{pavgb(pavgb(h_avg0, h_avg2), h_avg1)};

                            const T orig_val{b10x10[static_cast<ptrdiff_t>(r + 1) * stride + (c + 1)]};
                            filtered_val = std::clamp(filtered_val, orig_val - max_diff, orig_val + max_diff);

                            b8x8filtered[static_cast<ptrdiff_t>(r) * 8 + c] = std::clamp(filtered_val, min_clamp, max_clamp);
                        }
                    }
                }
                else
                {
                    T min_val{static_cast<T>((1 << bit_depth) - 1)}, max_val{0};
                    const int max_pixel_val{(1 << bit_depth) - 1};

                    for (int r{0}; r < 8; ++r)
                    {
                        for (int c{0}; c < 8; ++c)
                        {
                            const T p{b8x8[static_cast<ptrdiff_t>(r) * stride + c]};
                            min_val = std::min(min_val, p);
                            max_val = std::max(max_val, p);
                        }
                    }

                    thr = (static_cast<int>(max_val) + static_cast<int>(min_val) + 1) >> 1;
                    const int qp_scaled{QP * (1 << (bit_depth - 8))};
                    const int max_diff{qp_scaled >> 1};

                    auto pavgb{[](int a, int b) { return (a + b + 1) >> 1; }};

                    for (int r{0}; r < 8; ++r)
                    {
                        for (int c{0}; c < 8; ++c)
                        {
                            const T* p_row0{&b10x10[static_cast<ptrdiff_t>(r + 0) * stride + c]};
                            const T* p_row1{&b10x10[static_cast<ptrdiff_t>(r + 1) * stride + c]};
                            const T* p_row2{&b10x10[static_cast<ptrdiff_t>(r + 2) * stride + c]};

                            const int h_avg0{pavgb(pavgb(p_row0[0], p_row0[2]), p_row0[1])};
                            const int h_avg1{pavgb(pavgb(p_row1[0], p_row1[2]), p_row1[1])};
                            const int h_avg2{pavgb(pavgb(p_row2[0], p_row2[2]), p_row2[1])};
                            int filtered_val{pavgb(pavgb(h_avg0, h_avg2), h_avg1)};

                            const T orig_val{b10x10[static_cast<ptrdiff_t>(r + 1) * stride + (c + 1)]};
                            filtered_val =
                                std::clamp(filtered_val, static_cast<int>(orig_val) - max_diff, static_cast<int>(orig_val) + max_diff);

                            b8x8filtered[static_cast<ptrdiff_t>(r) * 8 + c] = static_cast<T>(std::clamp(filtered_val, 0, max_pixel_val));
                        }
                    }
                }

                for (int r{0}; r < 8; ++r)
                {
                    for (int c{0}; c < 8; ++c)
                    {
                        bool all_above{true};
                        bool all_below{true};

                        for (int dr{0}; dr < 3; ++dr)
                        {
                            for (int dc{0}; dc < 3; ++dc)
                            {
                                const T p{b10x10[static_cast<ptrdiff_t>(r + dr) * stride + (c + dc)]};

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

template<typename T>
void postprocess_impl(FrameData& frame, const PostProcessConfig& config) noexcept
{
    ChromaFormat format{config.format};
    const uint8_t** src{frame.src_planes.data()};
    uint8_t** dst{frame.dst_planes.data()};
    const int* src_strides{frame.src_strides.data()};
    const int* dst_strides{frame.dst_strides.data()};
    const int luma_width{frame.width[0]};
    const int chroma_width{frame.width[1]};
    const int luma_height{frame.height[0]};
    const int chroma_height{frame.height[1]};
    const bool iPP{config.iPP};
    const PPFlags mode{config.mode};
    const QP_STORE_T* AVS_RESTRICT QP_store{frame.qp_store};
    const int QP_stride{config.qp_stride};
    const int moderate_h{config.moderate_h};
    const int moderate_v{config.moderate_v};
    const int bit_depth{config.bit_depth};

    const int num_planes{(format == ChromaFormat::Y) ? 1 : 3};

    for (int plane_idx{0}; plane_idx < num_planes; ++plane_idx)
    {
        const bool is_luma{(plane_idx == 0)};
        const bool is_float_chroma{std::is_floating_point_v<T> && !is_luma};

        const uint8_t* AVS_RESTRICT p_src{src[plane_idx]};
        T* AVS_RESTRICT p_dst{reinterpret_cast<T*>(dst[plane_idx])};
        const int src_stride{static_cast<int>(src_strides[plane_idx] / sizeof(T))};
        const int dst_stride{static_cast<int>(dst_strides[plane_idx] / sizeof(T))};
        const int plane_width{static_cast<int>((is_luma ? luma_width : chroma_width) / sizeof(T))};
        const int plane_height{is_luma ? luma_height : chroma_height};
        [[maybe_unused]]
        const int src_stride_bytes{src_strides[plane_idx]};
        [[maybe_unused]]
        const int dst_stride_bytes{dst_strides[plane_idx]};
        [[maybe_unused]]
        const int width_bytes{(is_luma ? luma_width : chroma_width)};

        const PPFlags deblock_h_flag{is_luma ? PPFlags::DeblockYH : PPFlags::DeblockCH};
        const PPFlags deblock_v_flag{is_luma ? PPFlags::DeblockYV : PPFlags::DeblockCV};
        const PPFlags dering_flag{is_luma ? PPFlags::DeringY : PPFlags::DeringC};

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

        int v_size{plane_height};

        if (iPP)
            v_size >>= 1;

        // Deblocking loop
        for (int y{0}; y < v_size; y += 4)
        {
            if ((mode & PPFlags::DontCopy) != PPFlags::DontCopy)
            {
                auto* dst_bytes{reinterpret_cast<uint8_t*>(p_dst)};

                if (!iPP)
                {
                    fast_copy(p_src + static_cast<ptrdiff_t>(y) * src_stride_bytes, src_stride_bytes,
                        dst_bytes + static_cast<ptrdiff_t>(y) * dst_stride_bytes, dst_stride_bytes, width_bytes, 4);
                }
                else
                {
                    fast_copy(p_src + static_cast<ptrdiff_t>(y) * 2 * src_stride_bytes, src_stride_bytes,
                        dst_bytes + static_cast<ptrdiff_t>(y) * 2 * dst_stride_bytes, dst_stride_bytes, width_bytes, 8);
                }
            }

            if ((mode & deblock_h_flag) == deblock_h_flag)
            {
                const auto* qp_row{&QP_store[(y >> qp_y_shift) * QP_stride]};
                if (!iPP)
                {
                    deblock_horiz(p_dst + static_cast<ptrdiff_t>(y) * dst_stride, plane_width, dst_stride, qp_row, QP_stride,
                        chroma_flag_for_qp, moderate_h, bit_depth, is_float_chroma);
                }
                else
                {
                    deblock_horiz(p_dst + static_cast<ptrdiff_t>(y) * 2 * dst_stride, plane_width, dst_stride * 2,
                        &QP_store[(y >> qp_y_shift) * 2 * QP_stride], QP_stride * 2, chroma_flag_for_qp, moderate_h, bit_depth,
                        is_float_chroma);
                    deblock_horiz(p_dst + (static_cast<ptrdiff_t>(y) * 2 + 1) * dst_stride, plane_width, dst_stride * 2,
                        &QP_store[((y >> qp_y_shift) * 2 + 1) * QP_stride], QP_stride * 2, chroma_flag_for_qp, moderate_h, bit_depth,
                        is_float_chroma);
                }
            }

            if ((mode & deblock_v_flag) == deblock_v_flag)
            {
                if (y >= 12 && (y % 8) == 4)
                {
                    const int boundary_y{y - 4};
                    const auto* qp_row{&QP_store[(boundary_y >> qp_y_shift) * QP_stride]};

                    if (!iPP)
                    {
                        deblock_vert(p_dst + static_cast<ptrdiff_t>(boundary_y) * dst_stride, plane_width, dst_stride, qp_row, QP_stride,
                            chroma_flag_for_qp, moderate_v, bit_depth, is_float_chroma);
                    }
                    else
                    {
                        deblock_vert(p_dst + static_cast<ptrdiff_t>(boundary_y) * 2 * dst_stride, plane_width, dst_stride * 2,
                            &QP_store[(boundary_y >> qp_y_shift) * 2 * QP_stride], QP_stride * 2, chroma_flag_for_qp, moderate_v, bit_depth,
                            is_float_chroma);
                        deblock_vert(p_dst + (static_cast<ptrdiff_t>(boundary_y) * 2 + 1) * dst_stride, plane_width, dst_stride * 2,
                            &QP_store[((boundary_y >> qp_y_shift) * 2 + 1) * QP_stride], QP_stride * 2, chroma_flag_for_qp, moderate_v,
                            bit_depth, is_float_chroma);
                    }
                }
            }
        }

        // Deringing
        if ((mode & dering_flag) == dering_flag)
        {
            if (!iPP)
            {
                dering(p_dst, plane_width, plane_height, dst_stride, QP_store, QP_stride, chroma_flag_for_qp, bit_depth, is_float_chroma);
            }
            else
            {
                dering(p_dst, plane_width, v_size, dst_stride * 2, QP_store, QP_stride * 2, chroma_flag_for_qp, bit_depth, is_float_chroma);
                dering(p_dst + dst_stride, plane_width, v_size, dst_stride * 2, QP_store + QP_stride, QP_stride * 2, chroma_flag_for_qp,
                    bit_depth, is_float_chroma);
            }
        }
    }
}

// Explicit instantiations
template void postprocess_impl<uint8_t>(FrameData& frame, const PostProcessConfig& config) noexcept;
template void postprocess_impl<uint16_t>(FrameData& frame, const PostProcessConfig& config) noexcept;
template void postprocess_impl<float>(FrameData& frame, const PostProcessConfig& config) noexcept;

void __stdcall fast_copy(const uint8_t* AVS_RESTRICT src, int src_stride, uint8_t* AVS_RESTRICT dst, int dst_stride, int horizontal_size,
    int vertical_size) noexcept
{
    if (vertical_size <= 0)
        return;

    if (horizontal_size == src_stride && src_stride == dst_stride)
        memcpy(dst, src, static_cast<size_t>(horizontal_size) * vertical_size);
    else
    {
        const int bytes_to_copy{horizontal_size};

        for (int i{0}; i < vertical_size; ++i)
        {
            memcpy(dst, src, bytes_to_copy);
            dst += dst_stride;
            src += src_stride;
        }
    }
}
