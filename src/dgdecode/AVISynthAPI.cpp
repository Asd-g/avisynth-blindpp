/*
 *  Avisynth 2.5 API for MPEG2Dec3
 *
 *  Changes to fix frame dropping and random frame access are
 *  Copyright (C) 2003-2008 Donald A. Graft
 *
 *  Copyright (C) 2002-2003 Marc Fauconneau <marc.fd@liberysurf.fr>
 *
 *  based of the intial MPEG2Dec Avisytnh API Copyright (C) Mathias Born - May 2001
 *
 *  This file is part of DGMPGDec, a free MPEG-2 decoder
 *
 *  DGMPGDec is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2, or (at your option)
 *  any later version.
 *
 *  DGMPGDec is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with GNU Make; see the file COPYING.  If not, write to
 *  the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 */

#include <malloc.h>
#include "AvisynthAPI.h"

const AVS_Linkage* AVS_linkage = nullptr;


extern "C" __declspec(dllexport) const char* __stdcall
AvisynthPluginInit3(IScriptEnvironment* env, const AVS_Linkage* const vectors)
{
    AVS_linkage = vectors;

   env->AddFunction("BlindPP", "c[quant]i[cpu]i[cpu2]s[iPP]b[moderate_h]i[moderate_v]i", BlindPP::create, nullptr);

    return "BlindPP";
}
