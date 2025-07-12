#pragma once

#include <memory>
#include <string_view>

#include "PostProcess.h"
#include <avisynth.h>

struct AvsDeleter
{
    IScriptEnvironment* env;
    explicit AvsDeleter(IScriptEnvironment* e = nullptr) noexcept
        : env(e)
    {
    }
    void operator()(QP_STORE_T* ptr) const noexcept
    {
        if (ptr && env)
        {
            env->Free(ptr);
        }
    }
};

class BlindPP final : public GenericVideoFilter
{
    std::unique_ptr<QP_STORE_T[], AvsDeleter> QP_;
    const bool iPP_;
    PPFlags pp_mode_;
    const int moderate_h_;
    const int moderate_v_;
    const bool has_v8;
    const ChromaFormat format_;

public:
    BlindPP(PClip _child, int quant, int cpu, std::string_view cpu2, bool iPP, int moderate_h, int moderate_v, IScriptEnvironment* env);
    ~BlindPP() override = default;

    BlindPP(const BlindPP&) = delete;
    BlindPP& operator=(const BlindPP&) = delete;
    BlindPP(BlindPP&&) = delete;
    BlindPP& operator=(BlindPP&&) = delete;

    PVideoFrame __stdcall GetFrame(int n, IScriptEnvironment* env) override;
    int __stdcall SetCacheHints(int hints, int) noexcept override
    {
        return hints == CACHE_GET_MTMODE ? MT_MULTI_INSTANCE : 0;
    }
};
