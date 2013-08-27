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

#pragma once

#include "Alignment.hh"
#include "AlignmentScores.hh"

#include "blt_util/basic_matrix.hh"

#include "boost/foreach.hpp"

#include <cassert>


/// represents alignment of a query sequence which can switch over from reference1 to reference2
///
/// an empty alignment to one reference indicates that the entire alignment is to the other reference
///
template <typename ScoreType>
struct JumpAlignmentResult
{
    JumpAlignmentResult()
    {
        clear();
    }

    void
    clear()
    {
        score=0;
        jumpInsertSize=0;
        jumpRange=0;
        align1.clear();
        align2.clear();
    }


    ScoreType score;
    unsigned jumpInsertSize;
    unsigned jumpRange; ///< length of sequence over which jump would have the same score (left-most on align1 is reported)
    Alignment align1;
    Alignment align2;
};

template <typename ScoreType>
std::ostream& operator<<(std::ostream& os, JumpAlignmentResult<ScoreType>& alignment);



/// \brief a method to align a contig to two references
///
/// the alignment can make a single jump from reference1 to reference2
///
template <typename ScoreType>
struct GlobalJumpAligner
{
    GlobalJumpAligner(
        const AlignmentScores<ScoreType>& scores,
        const ScoreType jumpScore) :
        _scores(scores),
        _jumpScore(jumpScore)
    {}

    /// returns alignment path of query to reference
    template <typename SymIter>
    void
    align(
        const SymIter queryBegin, const SymIter queryEnd,
        const SymIter ref1Begin, const SymIter ref1End,
        const SymIter ref2Begin, const SymIter ref2End,
        JumpAlignmentResult<ScoreType>& result) const;


    /// read-only access to the aligner's scores:
    const AlignmentScores<ScoreType>&
    getScores() const
    {
        return _scores;
    }

    /// read-only access to the aligner's scores:
    const ScoreType&
    getJumpScore() const
    {
        return _jumpScore;
    }

    /// recorver a path alignment score without aligning, requires SEQ_MATCH style CIGAR
    ///
    ScoreType
    getPathScore(
        const ALIGNPATH::path_t& apath,
        const bool isScoreOffEdge = true) const
    {
        using namespace ALIGNPATH;

        ScoreType val(0);

        BOOST_FOREACH(const path_segment& ps, apath)
        {
            bool isIndel(false);
            switch (ps.type)
            {
            case MATCH:
                assert(0);
                break;
            case SEQ_MATCH:
                val += (_scores.match * ps.length);
                isIndel = false;
                break;
            case SEQ_MISMATCH:
                val += (_scores.mismatch * ps.length);
                isIndel = false;
                break;
            case INSERT:
            case DELETE:
                if (! isIndel) val += _scores.open;
                val += (_scores.extend * ps.length);
                isIndel = true;
                break;
            case SOFT_CLIP:
                if (isScoreOffEdge) val += (_scores.offEdge * ps.length);
                isIndel = false;
                break;
            default:
                break;
            }
        }
        return val;
    }

private:

    static
    uint8_t
    max3(
        ScoreType& max,
        const ScoreType v0,
        const ScoreType v1,
        const ScoreType v2)
    {
        max=v0;
        uint8_t ptr=0;
        if (v1>v0)
        {
            max=v1;
            ptr=1;
        }
        if (v2>max)
        {
            max=v2;
            ptr=2;
        }
        return ptr;
    }

    static
    uint8_t
    max4(
        ScoreType& max,
        const ScoreType v0,
        const ScoreType v1,
        const ScoreType v2,
        const ScoreType v3)
    {
        max=v0;
        uint8_t ptr=0;
        if (v1>v0)
        {
            max=v1;
            ptr=1;
        }
        if (v2>max)
        {
            max=v2;
            ptr=2;
        }
        if (v3>max)
        {
            max=v3;
            ptr=3;
        }
        return ptr;
    }

    const AlignmentScores<ScoreType> _scores;
    const ScoreType _jumpScore;

    // insert and delete are for seq1 wrt seq2
    struct ScoreVal
    {
        ScoreType match;
        ScoreType ins;
        ScoreType del;
        ScoreType jump;
    };

    struct PtrVal
    {
        uint8_t
        get(const AlignState::index_t i) const
        {
            switch (i)
            {
            case AlignState::MATCH:
                return match;
            case AlignState::INSERT:
                return ins;
            case AlignState::DELETE:
                return del;
            case AlignState::JUMP:
                return jump;
            default:
                assert(! "Unexpected Index Value");
                return 0;
            }
        }

        uint8_t match : 2;
        uint8_t ins : 2;
        uint8_t del : 2;
        uint8_t jump : 2;
    };

    // add the matrices here to reduce allocations over many alignment calls:
    typedef std::vector<ScoreVal> ScoreVec;
    mutable ScoreVec _score1;
    mutable ScoreVec _score2;

    typedef basic_matrix<PtrVal> PtrMat;
    mutable PtrMat _ptrMat1;
    mutable PtrMat _ptrMat2;
};


#include "alignment/GlobalJumpAlignerImpl.hh"

