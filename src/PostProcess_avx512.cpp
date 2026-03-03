#define SIMD_NAMESPACE avx512_opt
#include "PostProcess_simd256.hpp"

template void avx512_opt::deblock_vert<uint8_t>(uint8_t* AVS_RESTRICT image, int width, int stride, const QP_STORE_T* AVS_RESTRICT QP_store,
    int QP_stride, int chroma_flag_for_qp, int moderate_v, int bit_depth, bool is_float_chroma) noexcept;
template void avx512_opt::deblock_vert<uint16_t>(uint16_t* AVS_RESTRICT image, int width, int stride,
    const QP_STORE_T* AVS_RESTRICT QP_store, int QP_stride, int chroma_flag_for_qp, int moderate_v, int bit_depth,
    bool is_float_chroma) noexcept;
template void avx512_opt::deblock_vert<float>(float* AVS_RESTRICT image, int width, int stride, const QP_STORE_T* AVS_RESTRICT QP_store,
    int QP_stride, int chroma_flag_for_qp, int moderate_v, int bit_depth, bool is_float_chroma) noexcept;

template void avx512_opt::dering<uint8_t>(uint8_t* AVS_RESTRICT image, int width, int height, int stride,
    const QP_STORE_T* AVS_RESTRICT QP_store, int QP_stride, int chroma_flag_for_qp, int bit_depth, bool is_float_chroma) noexcept;
template void avx512_opt::dering<uint16_t>(uint16_t* AVS_RESTRICT image, int width, int height, int stride,
    const QP_STORE_T* AVS_RESTRICT QP_store, int QP_stride, int chroma_flag_for_qp, int bit_depth, bool is_float_chroma) noexcept;
template void avx512_opt::dering<float>(float* AVS_RESTRICT image, int width, int height, int stride,
    const QP_STORE_T* AVS_RESTRICT QP_store, int QP_stride, int chroma_flag_for_qp, int bit_depth, bool is_float_chroma) noexcept;
