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
// <https://github.com/downloads/sequencing/licenses/>.
//

///
/// \author Chris Saunders
///

#include "manta/SVLocusScanner.hh"
#include "blt_util/align_path_bam_util.hh"
#include "blt_util/align_path_util.hh"
#include "blt_util/log.hh"
#include "common/Exceptions.hh"

#include "boost/foreach.hpp"



struct simpleAlignment
{
    simpleAlignment() :
        is_fwd_strand(true),
        pos(0)
    {}

    simpleAlignment(const bam_record& bamRead) :
        is_fwd_strand(bamRead.is_fwd_strand()),
        pos(bamRead.pos()-1)
    {
        bam_cigar_to_apath(bamRead.raw_cigar(),bamRead.n_cigar(),path);
    }

    bool is_fwd_strand;
    pos_t pos;
    ALIGNPATH::path_t path;
};



static
void
getSVCandidatesFromRead(
    const ReadScannerOptions& opt,
    const SVLocusScanner::CachedReadGroupStats& rstats,
    const simpleAlignment& align,
    std::vector<SVCandidate>& candidates)
{
    using namespace ALIGNPATH;

    const std::pair<unsigned,unsigned> ends(get_match_edge_segments(align.path));

    unsigned pathIndex(0);
    unsigned readOffset(0);
    pos_t refHeadPos(align.pos);

    const unsigned pathSize(align.path.size());
    while (pathIndex<pathSize)
    {
        const path_segment& ps(align.path[pathIndex]);
        const bool isBeginEdge(pathIndex<ends.first);
        const bool isEndEdge(pathIndex>ends.second);
        const bool isEdgeSegment(isBeginEdge || isEndEdge);

        const bool isSwapStart(is_segment_swap_start(align.path,pathIndex));

        assert(ps.type != SKIP);
        assert(! (isEdgeSegment && isSwapStart));


        unsigned nPathSegments(1); // number of path segments consumed
        if (isEdgeSegment)
        {
            // edge inserts are allowed for intron adjacent and grouper reads, edge deletions for intron adjacent only:
            if (ps.type == INSERT) {
            } else if (ps.type == DELETE) {
            }
        }
        else if (isSwapStart)
        {
            const swap_info sinfo(align.path,pathIndex);
            const unsigned swap_size(std::max(sinfo.insert_length,sinfo.delete_length));

            nPathSegments = sinfo.n_seg;
            nPathSegments = process_swap(max_indel_size,al.path,read_seq,
                                 sppr,obs,sample_no,
                                 pathIndex,readOffset,refHeadPos);

        }
        else if (is_segment_type_indel(align.path[pathIndex].type))
        {
            process_simple_indel(max_indel_size,al.path,read_seq,
                                 sppr,obs,sample_no,
                                 pathIndex,read_offset,ref_head_pos);

        }

        for (unsigned i(0); i<nPathSegments; ++i)
        {
            increment_path(align.path,pathIndex,readOffset,refHeadPos);
        }
    }
}



static
void
getReadBreakendsImpl(
    const ReadScannerOptions& opt,
    const SVLocusScanner::CachedReadGroupStats& rstats,
    const bam_record& localRead,
    const bam_record* remoteReadPtr,
    std::vector<SVCandidate>& candidates,
    known_pos_range2& localEvidenceRange)
{
    const simpleAlignment align(localRead);
    const unsigned readSize(apath_read_length(apath));


    // 1) process any large indels in the localRead:
    {
        getSVCandidatesFromRead(opt,rstats,localRead,candidates);
        BOOST_FOREACH(const ALIGNPATH::path_segment ps, apath)
        {
            ps.type
        }
    }


    const unsigned localRefLength(apath_ref_length(apath));

    unsigned thisReadNoninsertSize(0);
    if (localRead.is_fwd_strand())
    {
        thisReadNoninsertSize=(readSize-apath_read_trail_size(apath));
    }
    else
    {
        thisReadNoninsertSize=(readSize-apath_read_lead_size(apath));
    }

    candidates.resize(1);

    SVBreakend& localBreakend(candidates[0].bp1);
    SVBreakend& remoteBreakend(candidates[0].bp2);

    localBreakend.readCount = 1;

    // if remoteRead is not available, estimate mate localRead size to be same as local,
    // and assume no clipping on mate localRead:
    unsigned remoteReadNoninsertSize(readSize);
    unsigned remoteRefLength(localRefLength);

    if (NULL != remoteReadPtr)
    {
        // if remoteRead is available, we can more accurately determine the size:
        const bam_record& remoteRead(*remoteReadPtr);

        ALIGNPATH::path_t remoteApath;
        bam_cigar_to_apath(remoteRead.raw_cigar(),remoteRead.n_cigar(),remoteApath);

        const unsigned remoteReadSize(apath_read_length(remoteApath));
        remoteRefLength = (apath_ref_length(remoteApath));

        if (remoteRead.is_fwd_strand())
        {
            remoteReadNoninsertSize=(remoteReadSize-apath_read_trail_size(remoteApath));
        }
        else
        {
            remoteReadNoninsertSize=(remoteReadSize-apath_read_lead_size(remoteApath));
        }

        remoteBreakend.readCount = 1;

        localBreakend.pairCount = 1;
        remoteBreakend.pairCount = 1;
    }

    const pos_t totalNoninsertSize(thisReadNoninsertSize+remoteReadNoninsertSize);
    const pos_t breakendSize(std::max(
        static_cast<pos_t>(opt.minPairBreakendSize),
        static_cast<pos_t>(rstats.breakendRegion.max-totalNoninsertSize)));

    {
        localBreakend.interval.tid = (localRead.target_id());

        const pos_t startRefPos(localRead.pos()-1);
        const pos_t endRefPos(startRefPos+localRefLength);
        // expected breakpoint range is from the end of the localRead alignment to the (probabilistic) end of the fragment:
        if (localRead.is_fwd_strand())
        {
            localBreakend.state = SVBreakendState::RIGHT_OPEN;
            localBreakend.interval.range.set_begin_pos(endRefPos);
            localBreakend.interval.range.set_end_pos(endRefPos + breakendSize);
        }
        else
        {
            localBreakend.state = SVBreakendState::LEFT_OPEN;
            localBreakend.interval.range.set_end_pos(startRefPos);
            localBreakend.interval.range.set_begin_pos(startRefPos - breakendSize);
        }

        localEvidenceRange.set_range(startRefPos,endRefPos);
    }

    // get remote breakend estimate:
    {
        remoteBreakend.interval.tid = (localRead.mate_target_id());

        const pos_t startRefPos(localRead.mate_pos()-1);
        pos_t endRefPos(startRefPos+remoteRefLength);
        if (localRead.is_mate_fwd_strand())
        {
            remoteBreakend.state = SVBreakendState::RIGHT_OPEN;
            remoteBreakend.interval.range.set_begin_pos(endRefPos);
            remoteBreakend.interval.range.set_end_pos(endRefPos + breakendSize);
        }
        else
        {
            remoteBreakend.state = SVBreakendState::LEFT_OPEN;
            remoteBreakend.interval.range.set_end_pos(startRefPos);
            remoteBreakend.interval.range.set_begin_pos(startRefPos - breakendSize);
        }
    }
}



static
void
getSVLociImpl(
    const ReadScannerOptions& opt,
    const SVLocusScanner::CachedReadGroupStats& rstats,
    const bam_record& bamRead,
    std::vector<SVLocus>& loci)
{
    using namespace illumina::common;

    loci.clear();
    std::vector<SVCandidate> candidates;
    known_pos_range2 localEvidenceRange;

    getReadBreakendsImpl(opt, rstats, bamRead, NULL, candidates, localEvidenceRange);

    BOOST_FOREACH(const SVCandidate& cand, candidates)
    {
        const SVBreakend& localBreakend(cand.bp1);
        const SVBreakend& remoteBreakend(cand.bp2);
        if ((0==localBreakend.interval.range.size()) ||
            (0==remoteBreakend.interval.range.size()))
        {
            std::ostringstream oss;
            oss << "Empty Breakend proposed from bam record.\n"
                << "\tlocal_breakend: " << localBreakend << "\n"
                << "\tremote_breakend: " << remoteBreakend << "\n"
                << "\tbam_record: " << bamRead << "\n";
            BOOST_THROW_EXCEPTION(LogicException(oss.str()));
        }

        SVLocus locus;
        // set local breakend estimate:
        NodeIndexType localBreakendNode(0);
        {
            localBreakendNode = locus.addNode(localBreakend.interval);
            locus.setNodeEvidence(localBreakendNode,localEvidenceRange);
        }

        // set remote breakend estimate:
        {
            const NodeIndexType remoteBreakendNode(locus.addRemoteNode(remoteBreakend.interval));
            locus.linkNodes(localBreakendNode,remoteBreakendNode);
            locus.mergeSelfOverlap();
        }

        loci.push_back(locus);
    }
}



SVLocusScanner::
SVLocusScanner(
    const ReadScannerOptions& opt,
    const std::string& statsFilename,
    const std::vector<std::string>& alignmentFilename) :
    _opt(opt)
{
    // pull in insert stats:
    _rss.read(statsFilename.c_str());

    // cache the insert stats we'll be looking up most often:
    BOOST_FOREACH(const std::string& file, alignmentFilename)
    {
        const boost::optional<unsigned> index(_rss.getGroupIndex(file));
        assert(index);
        const ReadGroupStats rgs(_rss.getStats(*index));

        _stats.resize(_stats.size()+1);
        CachedReadGroupStats& stat(_stats.back());
        {
            Range& breakend(stat.breakendRegion);
            breakend.min=rgs.fragSize.quantile(_opt.breakendEdgeTrimProb);
            breakend.max=rgs.fragSize.quantile((1-_opt.breakendEdgeTrimProb));

            if (breakend.min<0.) breakend.min = 0;
            assert(breakend.max>0.);
        }
        {
            Range& ppair(stat.properPair);
            ppair.min=rgs.fragSize.quantile(_opt.properPairTrimProb);
            ppair.max=rgs.fragSize.quantile((1-_opt.properPairTrimProb));

            if (ppair.min<0.) ppair.min = 0;

            assert(ppair.max>0.);
        }
    }
}



bool
SVLocusScanner::
isReadFiltered(const bam_record& bamRead) const
{
    if (bamRead.is_filter()) return true;
    if (bamRead.is_dup()) return true;
    if (bamRead.is_secondary()) return true;
    if (bamRead.map_qual() < _opt.minMapq) return true;
    return false;
}



bool
SVLocusScanner::
isProperPair(
    const bam_record& bamRead,
    const unsigned defaultReadGroupIndex) const
{
    if (bamRead.is_unmapped() || bamRead.is_mate_unmapped()) return false;
    if (bamRead.target_id() != bamRead.mate_target_id()) return false;

    const Range& ppr(_stats[defaultReadGroupIndex].properPair);
    const int32_t fragmentSize(std::abs(bamRead.template_size()));
    if (fragmentSize > ppr.max || fragmentSize < ppr.min) return false;

    if     (bamRead.pos() < bamRead.mate_pos())
    {
        if (! bamRead.is_fwd_strand()) return false;
        if (  bamRead.is_mate_fwd_strand()) return false;
    }
    else if (bamRead.pos() > bamRead.mate_pos())
    {
        if (  bamRead.is_fwd_strand()) return false;
        if (! bamRead.is_mate_fwd_strand()) return false;
    }
    else
    {
        if (bamRead.is_fwd_strand() == bamRead.is_mate_fwd_strand()) return false;
    }

    return true;
}



void
SVLocusScanner::
getSVLoci(
    const bam_record& bamRead,
    const unsigned defaultReadGroupIndex,
    std::vector<SVLocus>& loci) const
{
    loci.clear();

    if (! bamRead.is_chimeric())
    {
        if (std::abs(bamRead.template_size())<2000) return;
    }

    const CachedReadGroupStats& rstats(_stats[defaultReadGroupIndex]);
    getSVLociImpl(_opt, rstats, bamRead, loci);
}



void
SVLocusScanner::
getBreakendPair(
    const bam_record& localRead,
    const bam_record* remoteReadPtr,
    const unsigned defaultReadGroupIndex,
    std::vector<SVCandidate>& candidates) const
{
    const CachedReadGroupStats& rstats(_stats[defaultReadGroupIndex]);

    // throw evidence range away in this case
    known_pos_range2 evidenceRange;
    getReadBreakendsImpl(_opt, rstats, localRead, remoteReadPtr, candidates, evidenceRange);
}
