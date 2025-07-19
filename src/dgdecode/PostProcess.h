#pragma once

#include <array>
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

[[nodiscard]] AVS_FORCEINLINE constexpr PPFlags operator|(PPFlags a, PPFlags b) noexcept
{
    return static_cast<PPFlags>(static_cast<unsigned int>(a) | static_cast<unsigned int>(b));
}
[[nodiscard]] AVS_FORCEINLINE constexpr PPFlags operator&(PPFlags a, PPFlags b) noexcept
{
    return static_cast<PPFlags>(static_cast<unsigned int>(a) & static_cast<unsigned int>(b));
}
AVS_FORCEINLINE PPFlags& operator|=(PPFlags& a, PPFlags b) noexcept
{
    a = a | b;
    return a;
}

using QP_STORE_T = int;

enum class ChromaFormat
{
    Y,
    YUV420,
    YUV422,
    YUV444
};

// New struct for configuration parameters (constant per filter instance)
struct PostProcessConfig
{
    PPFlags mode;
    int moderate_h;
    int moderate_v;
    ChromaFormat format;
    bool iPP;
    int bit_depth;
    int qp_stride;
};

struct FrameData
{
    std::array<const uint8_t*, 3> src_planes;
    std::array<int, 3> src_strides;
    std::array<uint8_t*, 3> dst_planes;
    std::array<int, 3> dst_strides;
    std::array<int, 2> width;
    std::array<int, 2> height;
    const QP_STORE_T* qp_store;
};

using PostProcessFunction = void (*)(FrameData& frame, const PostProcessConfig& config) noexcept;

template<typename T>
void postprocess_impl(FrameData& frame, const PostProcessConfig& config) noexcept;

void __stdcall fast_copy(const uint8_t* AVS_RESTRICT src, int src_stride, uint8_t* AVS_RESTRICT dst, int dst_stride, int horizontal_size,
    int vertical_size) noexcept;
