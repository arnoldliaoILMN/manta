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

#include "manta/SomaticSVScoreInfo.hh"

#include "boost/foreach.hpp"

#include <iostream>



std::ostream&
operator<<(
    std::ostream& os,
    const SVSampleInfo& si)
{
    static const char indent('\t');
    os << "SVSampleInfo:\n"
       << indent << "altAlleleBp1SpanReads: " << si.altAlleleBp1SpanReads << "\n"
       << indent << "altAlleleBp2SpanReads: " << si.altAlleleBp2SpanReads << "\n"
       << indent << " altAlleleSpanPairs: " << si.altAlleleSpanPairs << "\n"
       << indent << " refAlleleBp1SpanPairs: " << si.refAlleleBp1SpanPairs << "\n"
       << indent << " refAlleleBp2SpanPairs: " << si.refAlleleBp2SpanPairs << "\n"
    ;
    return os;
}



std::ostream&
operator<<(
    std::ostream& os,
    const SomaticSVScoreInfo& ssi)
{
    os << "SomaticSVScoreInfo bp1MaxDepth=" << ssi.bp1MaxDepth << " bp2MaxDepth=" << ssi.bp2MaxDepth << " somaticScore=" << ssi.somaticScore << "\n";
    os << "Tumor sample info " << ssi.tumor;
    os << "Normal sample info " << ssi.normal;
    BOOST_FOREACH(const std::string& filter, ssi.filters)
    {
        os << " " << filter;
    }
    return os;
}
