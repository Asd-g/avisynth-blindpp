/*
 *  Avisynth 2.5 API for MPEG2Dec3
 *
 *  Copyright (C) 2002-2003 Marc Fauconneau <marc.fd@liberysurf.fr>
 *
 *  based of the intial MPEG2Dec Avisytnh API Copyright (C) Mathias Born - May 2001
 *
 *  This file is part of MPEG2Dec3, a free MPEG-2 decoder
 *
 *  MPEG2Dec3 is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2, or (at your option)
 *  any later version.
 *
 *  MPEG2Dec3 is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with GNU Make; see the file COPYING.  If not, write to
 *  the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 */

#ifndef MPEG2DECPLUS_AVS_API_H
#define MPEG2DECPLUS_AVS_API_H

#include <cstdint>
#include "avisynth.h"


class BlindPP : public GenericVideoFilter
{
    int* QP;
    bool iPP;
    int PP_MODE;
    int moderate_h, moderate_v;

public:
    BlindPP(AVSValue args, IScriptEnvironment* env);
    PVideoFrame __stdcall GetFrame(int n, IScriptEnvironment* env);
    ~BlindPP() {}
    int __stdcall SetCacheHints(int hints, int) { return hints == CACHE_GET_MTMODE ? MT_NICE_FILTER : 0; }
    static AVSValue __cdecl create(AVSValue args, void*, IScriptEnvironment* env);
};

#endif
