#pragma once

#include <algorithm>

#include "PostProcess.h"

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

AVS_FORCEINLINE void fast_copy(const uint8_t* AVS_RESTRICT src, int src_stride, uint8_t* AVS_RESTRICT dst, int dst_stride,
    int horizontal_size, int vertical_size) noexcept
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
