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

#include "boost/test/unit_test.hpp"

#include "manta/SVLocusScanner.cpp"


BOOST_AUTO_TEST_SUITE( test_SVLocusScanner )


BOOST_AUTO_TEST_CASE( test_getSVCandidatesFromReadIndels )
{
    ReadScannerOptions opt;

    ALIGNPATH::path_t inputPath;
    cigar_to_apath("100M2000D100M",inputPath);

    bam_record bamRead;
    bam1_t* bamDataPtr(bamRead.get_data());
    edit_bam_cigar(inputPath,*bamDataPtr);

    SimpleAlignment align(bamRead);

    std::vector<SVCandidate> candidates;

    getSVCandidatesFromReadIndels(opt,align,0,candidates);


    BOOST_REQUIRE_EQUAL(candidates.size(),1u);
    BOOST_REQUIRE(candidates[0].bp1.interval.range.is_pos_intersect(100));
    BOOST_REQUIRE(candidates[0].bp2.interval.range.is_pos_intersect(2100));
}


BOOST_AUTO_TEST_CASE( test_getSVCandidatesFromSemiAligned )
{
    ReadScannerOptions opt;

    ALIGNPATH::path_t inputPath;
    cigar_to_apath("3=7X",inputPath);

    bam_record bamRead;
    bam1_t* bamDataPtr(bamRead.get_data());
    edit_bam_cigar(inputPath,*bamDataPtr);

    static const char testSeq[] = "ACGTTACGTT";

    // initialize test qual array to all Q30's:
    static const unsigned seqSize((sizeof(testSeq)-1)/sizeof(char));
    uint8_t qual[seqSize];
    for(unsigned i(0);i<seqSize;++i) { qual[i] = 30; }

    edit_bam_read_and_quality(testSeq,qual,*bamDataPtr);

    SimpleAlignment align(bamRead);
    align.pos = 500;

    std::vector<SVCandidate> candidates;

    getSVCandidatesFromSemiAligned(opt,bamRead,align,candidates);

    BOOST_REQUIRE_EQUAL(candidates.size(),1u);
    BOOST_REQUIRE(candidates[0].bp1.interval.range.is_pos_intersect(500));
    BOOST_REQUIRE_EQUAL(candidates[0].bp1.interval.range.begin_pos(),480);
    BOOST_REQUIRE_EQUAL(candidates[0].bp1.interval.range.end_pos(),520);
}


BOOST_AUTO_TEST_SUITE_END()

