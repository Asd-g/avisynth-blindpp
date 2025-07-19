#include <stdexcept>
#include <string>
#include <vector>

#include "BlindPP.h"

namespace
{
    ChromaFormat get_chroma_format(const VideoInfo& vi)
    {
        if (vi.IsY())
            return ChromaFormat::Y;
        if (vi.Is420())
            return ChromaFormat::YUV420;
        if (vi.Is422())
            return ChromaFormat::YUV422;
        if (vi.Is444())
            return ChromaFormat::YUV444;
        // Should be unreachable
        return ChromaFormat::YUV420;
    }

    PPFlags get_pp_flags(int cpu, std::string_view cpu2, IScriptEnvironment* env)
    {
        if (cpu2.empty())
        {
            switch (cpu)
            {
            case 0:
                return PPFlags::None;
            case 1:
                return PPFlags::DeblockYH;
            case 2:
                return PPFlags::DeblockYH | PPFlags::DeblockYV;
            case 3:
                return PPFlags::DeblockYH | PPFlags::DeblockYV | PPFlags::DeblockCH;
            case 4:
                return PPFlags::DeblockYH | PPFlags::DeblockYV | PPFlags::DeblockCH | PPFlags::DeblockCV;
            case 5:
                return PPFlags::DeblockYH | PPFlags::DeblockYV | PPFlags::DeblockCH | PPFlags::DeblockCV | PPFlags::DeringY;
            case 6:
                return PPFlags::DeblockYH | PPFlags::DeblockYV | PPFlags::DeblockCH | PPFlags::DeblockCV | PPFlags::DeringY |
                       PPFlags::DeringC;
            default:
                env->ThrowError("BlindPP: 'cpu' level must be between 0 and 6.");
            }
        }

        if (cpu2.length() != 6)
            env->ThrowError("BlindPP: 'cpu2' must be a string of 6 characters or empty.");

        PPFlags flags{PPFlags::None};
        if (cpu2[0] == 'x' || cpu2[0] == 'X')
            flags |= PPFlags::DeblockYH;
        if (cpu2[1] == 'x' || cpu2[1] == 'X')
            flags |= PPFlags::DeblockYV;
        if (cpu2[2] == 'x' || cpu2[2] == 'X')
            flags |= PPFlags::DeblockCH;
        if (cpu2[3] == 'x' || cpu2[3] == 'X')
            flags |= PPFlags::DeblockCV;
        if (cpu2[4] == 'x' || cpu2[4] == 'X')
            flags |= PPFlags::DeringY;
        if (cpu2[5] == 'x' || cpu2[5] == 'X')
            flags |= PPFlags::DeringC;
        return flags;
    }
} // namespace

BlindPP::BlindPP(PClip _child, int quant, int cpu, std::string_view cpu2, bool iPP, int moderate_h, int moderate_v, IScriptEnvironment* env)
    : GenericVideoFilter(_child),
      component_size_(vi.ComponentSize()),
      has_v8_(env->FunctionExists("propShow")),
      process_func_([&] {
          switch (component_size_)
          {
          case 1:
              return &postprocess_impl<uint8_t>;
          case 2:
              return &postprocess_impl<uint16_t>;
          case 4:
              return &postprocess_impl<float>;
          default:
              env->ThrowError("BlindPP: Unsupported pixel format component size.");
          }
      }()),
      pp_config_{.mode = get_pp_flags(cpu, cpu2, env),
          .moderate_h = moderate_h,
          .moderate_v = moderate_v,
          .format = get_chroma_format(vi),
          .iPP = iPP,
          .bit_depth = vi.BitsPerComponent(),
          .qp_stride = vi.width / 16}
{
    if (vi.width % 16 != 0)
        env->ThrowError("BlindPP: Need mod16 width.");

    if (vi.height % 16 != 0)
        env->ThrowError("BlindPP: Need mod16 height.");

    if (vi.IsRGB() || !vi.IsPlanar())
        env->ThrowError("BlindPP: Only planar Y/YUV formats are supported.");

    if (pp_config_.moderate_h < 0 || pp_config_.moderate_h > 255)
        env->ThrowError("BlindPP: 'moderate_h' must be between 0 and 255 (inclusive).");

    if (pp_config_.moderate_v < 0 || pp_config_.moderate_v > 255)
        env->ThrowError("BlindPP: 'moderate_v' must be between 0 and 255 (inclusive).");

    if (quant < 0 || quant > 63)
        env->ThrowError("BlindPP: 'quant' must be between 0 and 63 (inclusive).");

    const size_t n_blocks{static_cast<size_t>(vi.width) * vi.height / 256};
    QP_STORE_T* QP_raw{static_cast<QP_STORE_T*>(env->Allocate(n_blocks * sizeof(QP_STORE_T), 64, AVS_POOLED_ALLOC))};

    if (!QP_raw)
        env->ThrowError("BlindPP: Failed to allocate memory for QP store.");

    QP_ = std::unique_ptr<QP_STORE_T[], AvsDeleter>(QP_raw, AvsDeleter{env});
    std::fill_n(QP_.get(), n_blocks, quant);
}

PVideoFrame __stdcall BlindPP::GetFrame(int n, IScriptEnvironment* env)
{
    PVideoFrame src_frame{child->GetFrame(n, env)};
    PVideoFrame dst_frame{has_v8_ ? env->NewVideoFrameP(vi, &src_frame) : env->NewVideoFrame(vi)};

    FrameData frame_data{};
    frame_data.qp_store = QP_.get();

    const int num_planes{(pp_config_.format == ChromaFormat::Y) ? 1 : 3};
    constexpr int plane_indices[]{PLANAR_Y, PLANAR_U, PLANAR_V};

    for (int i{0}; i < num_planes; ++i)
    {
        const int plane{plane_indices[i]};
        frame_data.src_strides[i] = src_frame->GetPitch(plane);
        frame_data.dst_strides[i] = dst_frame->GetPitch(plane);
        frame_data.src_planes[i] = src_frame->GetReadPtr(plane);
        frame_data.dst_planes[i] = dst_frame->GetWritePtr(plane);

        if (i < 2)
        {
            frame_data.width[i] = src_frame->GetRowSize(plane);
            frame_data.height[i] = src_frame->GetHeight(plane);
        }
    }

    process_func_(frame_data, pp_config_);

    return dst_frame;
}

AVSValue __cdecl BlindPP_create(AVSValue args, void*, IScriptEnvironment* env)
{
    enum class Arg : int
    {
        Clip = 0,
        Quant,
        Cpu,
        Cpu2,
        IPP,
        ModerateH,
        ModerateV
    };

    return new BlindPP(args[static_cast<int>(Arg::Clip)].AsClip(), args[static_cast<int>(Arg::Quant)].AsInt(2),
        args[static_cast<int>(Arg::Cpu)].AsInt(6), args[static_cast<int>(Arg::Cpu2)].AsString(""),
        args[static_cast<int>(Arg::IPP)].AsBool(false), args[static_cast<int>(Arg::ModerateH)].AsInt(20),
        args[static_cast<int>(Arg::ModerateV)].AsInt(40), env);
}

const AVS_Linkage* AVS_linkage{};

extern "C" __declspec(dllexport) const char* __stdcall AvisynthPluginInit3(IScriptEnvironment* env, const AVS_Linkage* const vectors)
{
    AVS_linkage = vectors;
    env->AddFunction("BlindPP", "c[quant]i[cpu]i[cpu2]s[iPP]b[moderate_h]i[moderate_v]i", BlindPP_create, nullptr);
    return "BlindPP";
}
