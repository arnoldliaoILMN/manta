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

#include "blt_util/ReadKey.hh"

#include <iosfwd>
#include <string>
#include <set>


/// sample-specific evidence info
struct SVSampleInfo
{
    SVSampleInfo() :
        altAlleleBp1SpanReads(0),
        altAlleleBp2SpanReads(0),
        altAlleleSpanPairs(0),
        refAlleleBp1SpanPairs(0),
        refAlleleBp2SpanPairs(0)
    {}

    void
    clear()
    {
        altAlleleBp1SpanReads=0;
        altAlleleBp2SpanReads=0;
        altAlleleSpanPairs=0;
        refAlleleBp1SpanPairs=0;
        refAlleleBp2SpanPairs=0;
    }

    unsigned altAlleleBp1SpanReads;
    unsigned altAlleleBp2SpanReads;
    unsigned altAlleleSpanPairs;

    unsigned refAlleleBp1SpanPairs;
    unsigned refAlleleBp2SpanPairs;
};

std::ostream&
operator<<(std::ostream& os, const SVSampleInfo& si);


/// consolidate all somatic scoring results applied to an SV candidate
struct SomaticSVScoreInfo
{
    SomaticSVScoreInfo() :
        bp1MaxDepth(0),
        bp2MaxDepth(0),
        somaticScore(0)
    {}

    void
    clear()
    {
        tumor.clear();
        normal.clear();
        filters.clear();

        bp1MaxDepth=0;
        bp2MaxDepth=0;
        somaticScore=0;
    }

    SVSampleInfo tumor;
    SVSampleInfo normal;

    std::set<std::string> filters;

    unsigned bp1MaxDepth;
    unsigned bp2MaxDepth;
    unsigned somaticScore;
};

std::ostream&
operator<<(std::ostream& os, const SomaticSVScoreInfo& ssi);
