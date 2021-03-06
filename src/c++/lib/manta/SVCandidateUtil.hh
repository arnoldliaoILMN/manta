// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Manta
// Copyright (c) 2013 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/sequencing/licenses/>
//

///
/// \author Chris Saunders
///

#pragma once

#include "manta/SVCandidate.hh"


/// returns true if sv is below minimum size:
///
bool
isSVBelowMinSize(
    const SVCandidate& sv,
    const unsigned minSize);


namespace SV_TYPE
{
enum index_t
{
    UNKNOWN,
    INTERTRANSLOC,
    INVERSION,
    INDEL,
    TANDUP,
    COMPLEX
};

inline
const char*
label(const index_t idx)
{
    switch (idx)
    {
    case UNKNOWN:
        return "UNKNOWN";
    case INTERTRANSLOC:
        return "INTERTRANSLOC";
    case INVERSION:
        return "INVERSION";
    case INDEL:
        return "INDEL";
    case TANDUP:
        return "TANDUP";
    case COMPLEX:
        return "COMPLEX";
    default:
        return "UNKNOWN";
    }
}

}


SV_TYPE::index_t
getSVType(const SVCandidate& sv);
