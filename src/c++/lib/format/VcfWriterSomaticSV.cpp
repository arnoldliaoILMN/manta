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

#include "format/VcfWriterSomaticSV.hh"

#include "boost/algorithm/string/join.hpp"

#include "manta/SomaticSVScoreInfo.hh"


void
VcfWriterSomaticSV::
addHeaderInfo() const
{
    _os << "##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description=\"Somatic mutation\">\n";
    _os << "##INFO=<ID=SOMATICSCORE,Number=1,Type=Integer,Description=\"Somatic variant Quality score\">\n";
}



void
VcfWriterSomaticSV::
addHeaderFilters() const
{
    if (_isMaxDepthFilter)
    {
        _os << "##FILTER=<ID=" << _somaticOpt.maxDepthFilterLabel << ",Description=\"Normal sample site depth is greater than " << _somaticOpt.maxDepthFactor << "x the mean chromosome depth near one or both variant breakends\">\n";
    }
}



void
VcfWriterSomaticSV::
modifyInfo(
    const bool isFirstOfPair,
    const SVCandidateSetData& /*svData*/,
    const SVCandidateAssemblyData& /*adata*/,
    std::vector<std::string>& infotags) const
{
    assert(_ssInfoPtr != NULL);
    const SomaticSVScoreInfo& ssInfo(*_ssInfoPtr);

    infotags.push_back("SOMATIC");
    infotags.push_back( str(boost::format("SOMATICSCORE=%i") % ssInfo.somaticScore) );
    infotags.push_back( str(boost::format("NORMAL_PAIR_SUPPORT=%i") % ssInfo.normal.altAlleleSpanPairs) );
    infotags.push_back( str(boost::format("TUMOR_PAIR_SUPPORT=%i") % ssInfo.tumor.altAlleleSpanPairs) );
    infotags.push_back( str(boost::format("NORMAL_BND_PAIR_SUPPORT=%i") %
                            (isFirstOfPair ? ssInfo.normal.altAlleleBp1SpanReads : ssInfo.normal.altAlleleBp2SpanReads) ) );
    infotags.push_back( str(boost::format("TUMOR_BND_PAIR_SUPPORT=%i") %
                            (isFirstOfPair ? ssInfo.tumor.altAlleleBp1SpanReads : ssInfo.tumor.altAlleleBp2SpanReads) ) );
    infotags.push_back( str(boost::format("BND_DEPTH=%i") %
                            (isFirstOfPair ? ssInfo.bp1MaxDepth : ssInfo.bp2MaxDepth) ) );
    infotags.push_back( str(boost::format("MATE_BND_DEPTH=%i") %
                            (isFirstOfPair ? ssInfo.bp2MaxDepth : ssInfo.bp1MaxDepth) ) );
}



std::string
VcfWriterSomaticSV::
getFilter() const
{
    assert(_ssInfoPtr != NULL);
    const SomaticSVScoreInfo& ssInfo(*_ssInfoPtr);

    if (ssInfo.filters.empty())
    {
        return "PASS";
    }
    else
    {
        return boost::algorithm::join(ssInfo.filters, ";");
    }
}


void
VcfWriterSomaticSV::
writeSV(
    const EdgeInfo& edge,
    const SVCandidateSetData& svData,
    const SVCandidateAssemblyData& adata,
    const SVCandidate& sv,
    const SomaticSVScoreInfo& ssInfo)
{
    //TODO: this is a lame way to customize subclass behavior:
    _ssInfoPtr=&ssInfo;
    writeSVCore(edge, svData, adata, sv);
    _ssInfoPtr=NULL;
}
