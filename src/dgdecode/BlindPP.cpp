#include <string>
#include <vector>

#include "BlindPP.h"

namespace
{
    [[nodiscard]] ChromaFormat GetChromaFormat(const VideoInfo& vi) noexcept
    {
        if (vi.IsY8())
            return ChromaFormat::Y8;

        if (vi.IsYV12())
            return ChromaFormat::YUV420;

        if (vi.IsYV16())
            return ChromaFormat::YUV422;

        if (vi.IsYV24())
            return ChromaFormat::YUV444;

        return ChromaFormat::YUV420;
    }
} // namespace

BlindPP::BlindPP(PClip _child, int quant, int cpu, std::string_view cpu2, bool iPP, int moderate_h, int moderate_v, IScriptEnvironment* env)
    : GenericVideoFilter(_child),
      iPP_(iPP),
      pp_mode_(PPFlags::None),
      moderate_h_(moderate_h),
      moderate_v_(moderate_v),
      has_v8(env->FunctionExists("propShow")),
      format_(GetChromaFormat(vi))
{
    if (vi.width % 16 != 0)
        env->ThrowError("BlindPP: Need mod16 width.");

    if (vi.height % 16 != 0)
        env->ThrowError("BlindPP: Need mod16 height.");

    if (vi.IsRGB() || (!vi.IsPlanar() && !vi.IsRGB()))
        env->ThrowError("BlindPP: Only Y8, YV12, YV16, and YV24 are supported.");

    if (moderate_h_ < 0 || moderate_h_ > 255)
        env->ThrowError("BlindPP: 'moderate_h' must be between 0 and 255 (inclusive).");

    if (moderate_v_ < 0 || moderate_v_ > 255)
        env->ThrowError("BlindPP: 'moderate_v' must be between 0 and 255 (inclusive).");

    if (quant < 0 || quant > 31)
        env->ThrowError("BlindPP: 'quant' must be between 0 and 31 (inclusive).");

    const size_t n_blocks = static_cast<size_t>(vi.width) * vi.height / 256;
    QP_STORE_T* QP_raw = static_cast<QP_STORE_T*>(env->Allocate(n_blocks * sizeof(QP_STORE_T), 64, AVS_POOLED_ALLOC));

    if (!QP_raw)
        env->ThrowError("BlindPP: Failed to allocate memory for QP store.");

    QP_ = std::unique_ptr<QP_STORE_T[], AvsDeleter>(QP_raw, AvsDeleter{env});

    std::fill_n(QP_.get(), n_blocks, quant);

    switch (cpu)
    {
    case 0:
        pp_mode_ = PPFlags::None;
        break;
    case 1:
        pp_mode_ = PPFlags::DeblockYH;
        break;
    case 2:
        pp_mode_ = PPFlags::DeblockYH | PPFlags::DeblockYV;
        break;
    case 3:
        pp_mode_ = PPFlags::DeblockYH | PPFlags::DeblockYV | PPFlags::DeblockCH;
        break;
    case 4:
        pp_mode_ = PPFlags::DeblockYH | PPFlags::DeblockYV | PPFlags::DeblockCH | PPFlags::DeblockCV;
        break;
    case 5:
        pp_mode_ = PPFlags::DeblockYH | PPFlags::DeblockYV | PPFlags::DeblockCH | PPFlags::DeblockCV | PPFlags::DeringY;
        break;
    case 6:
        pp_mode_ = PPFlags::DeblockYH | PPFlags::DeblockYV | PPFlags::DeblockCH | PPFlags::DeblockCV | PPFlags::DeringY | PPFlags::DeringC;
        break;
    default:
        env->ThrowError("BlindPP: 'cpu' level must be between 0 and 6.");
    }

    if (cpu2.length() == 6)
    {
        pp_mode_ = PPFlags::None;
        if (cpu2[0] == 'x' || cpu2[0] == 'X')
        {
            pp_mode_ |= PPFlags::DeblockYH;
        }
        if (cpu2[1] == 'x' || cpu2[1] == 'X')
        {
            pp_mode_ |= PPFlags::DeblockYV;
        }
        if (cpu2[2] == 'x' || cpu2[2] == 'X')
        {
            pp_mode_ |= PPFlags::DeblockCH;
        }
        if (cpu2[3] == 'x' || cpu2[3] == 'X')
        {
            pp_mode_ |= PPFlags::DeblockCV;
        }
        if (cpu2[4] == 'x' || cpu2[4] == 'X')
        {
            pp_mode_ |= PPFlags::DeringY;
        }
        if (cpu2[5] == 'x' || cpu2[5] == 'X')
        {
            pp_mode_ |= PPFlags::DeringC;
        }
    }
}

PVideoFrame __stdcall BlindPP::GetFrame(int n, IScriptEnvironment* env)
{
    PVideoFrame src_frame = child->GetFrame(n, env);
    PVideoFrame dst_frame = has_v8 ? env->NewVideoFrameP(vi, &src_frame) : env->NewVideoFrameP(vi, &src_frame);

    const int num_planes = (format_ == ChromaFormat::Y8) ? 1 : 3;
    int src_strides[3]{};
    int dst_strides[3]{};
    int width[2]{};
    int height[2]{};
    const uint8_t* src_planes[3]{};
    uint8_t* dst_planes[3]{};

    constexpr int plane_indices[] = {PLANAR_Y, PLANAR_U, PLANAR_V};

    for (int i = 0; i < num_planes; ++i)
    {
        const int plane = plane_indices[i];
        src_strides[i] = src_frame->GetPitch(plane);
        dst_strides[i] = dst_frame->GetPitch(plane);
        src_planes[i] = src_frame->GetReadPtr(plane);
        dst_planes[i] = dst_frame->GetWritePtr(plane);

        if (i < 2)
        {
            width[i] = src_frame->GetRowSize(plane);
            height[i] = src_frame->GetHeight(plane);
        }
    }

    postprocess(src_planes, src_strides, dst_planes, dst_strides, width, height, QP_.get(), width[0] / 16, pp_mode_, moderate_h_,
        moderate_v_, format_, iPP_);

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
