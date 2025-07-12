#pragma once

#include <cstdint>

#include <avs/config.h>

enum class PPFlags : unsigned int
{
    None = 0,
    DeblockYH = 0x00000001, // Luma horizontal deblocking
    DeblockYV = 0x00000002, // Luma vertical deblocking
    DeblockCH = 0x00000004, // Chroma horizontal deblocking
    DeblockCV = 0x00000008, // Chroma vertical deblocking
    DeringY = 0x00000010,   // Luma deringing
    DeringC = 0x00000020,   // Chroma deringing
    DontCopy = 0x10000000,  // Postprocessor will not copy src -> dst
};

[[nodiscard]] inline constexpr PPFlags operator|(PPFlags a, PPFlags b) noexcept
{
    return static_cast<PPFlags>(static_cast<unsigned int>(a) | static_cast<unsigned int>(b));
}
[[nodiscard]] inline constexpr PPFlags operator&(PPFlags a, PPFlags b) noexcept
{
    return static_cast<PPFlags>(static_cast<unsigned int>(a) & static_cast<unsigned int>(b));
}
inline PPFlags& operator|=(PPFlags& a, PPFlags b) noexcept
{
    a = a | b;
    return a;
}

using QP_STORE_T = int;

enum class ChromaFormat
{
    Y8,
    YUV420,
    YUV422,
    YUV444
};

void postprocess(const uint8_t* const AVS_RESTRICT src[3], const int src_strides[3], uint8_t* const AVS_RESTRICT dst[3],
    const int dst_strides[3], const int width[2], int height[2], const QP_STORE_T* AVS_RESTRICT QP_store, int QP_stride, PPFlags mode,
    int moderate_h, int moderate_v, ChromaFormat format, bool iPP) noexcept;

void __stdcall fast_copy(const uint8_t* AVS_RESTRICT src, int src_stride, uint8_t* AVS_RESTRICT dst, int dst_stride, int horizontal_size,
    int vertical_size) noexcept;
