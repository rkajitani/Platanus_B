/*
Copyright (C) 2018 Itoh Laboratory, Tokyo Institute of Technology

This file is part of Platanus_B.

Platanus_B is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

Platanus_B is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with Platanus_B; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

#include "solveDBG.h"
#include "seqlib.h"
#include "kmer.h"
#include <cmath>
#include <algorithm>
#include <sstream>
#include <climits>
#include <cfloat>
#include <iomanip>

using std::vector;
using std::string;
using std::unordered_map;
using std::cerr;
using std::endl;


//////////////////////////////////////////////////////////////////////////////////////
// const parameter define
//////////////////////////////////////////////////////////////////////////////////////
const double SolveDBG::MIN_LONG_READ_LENGTH_CUTOFF_FACTOR = 1;
const double SolveDBG::MAX_LONG_READ_LENGTH_CUTOFF_FACTOR = 8;
const double SolveDBG::LONG_READ_LENGTH_CUTOFF_UNIT_FACTOR = 0.1;
const long SolveDBG::LONG_READ_MIN_ALIGNMENT_LENGTH = 1000;
const double SolveDBG::LONG_READ_MIN_ALIGNMENT_COVERAGE = 0.8;
const double SolveDBG::LONG_READ_MIN_ALIGNMENT_IDENTITY = 0.8;

//////////////////////////////////////////////////////////////////////////////////////
// default constructor
//////////////////////////////////////////////////////////////////////////////////////
SolveDBG::SolveDBG()
: Scaffold()
{
    optionSingleArgs["-o"] = "out";
    optionSingleArgs["-e"] = "";
    optionSingleArgs["-v"] = "32";
    optionSingleArgs["-k"] = "1";
    optionSingleArgs["-l"] = "3";
    optionSingleArgs["-u"] = "0";
    optionSingleArgs["-L"] = "200000";
    optionSingleArgs["-t"] = "1";
    optionMultiArgs["-c"] = vector<string>();
    optionMultiArgs["-masked"] = vector<string>();
    optionMultiArgs["-b"] = vector<string>();
    optionMultiArgs["-p"] = vector<string>();
    optionMultiArgs["-ont"] = vector<string>();
    optionMultiArgs["-gc"] = vector<string>();
    optionMultiArgs["-x"] = vector<string>();
    optionMultiArgs["-X"] = vector<string>();

//    optionMultiArgs["-s"] = vector<string>(3);
//    optionMultiArgs["-s"][0] = "32";
//    optionMultiArgs["-s"][1] = "64";
//    optionMultiArgs["-s"][2] = "96";
    optionMultiArgs["-s"] = vector<string>(1);
    optionMultiArgs["-s"][0] = "32";

//    optionMultiArgs["-S"] = vector<string>(1);
//    optionMultiArgs["-S"][0] = "20";

	optionSingleArgs["-combine_l"] = "10000";
	optionSingleArgs["-combine_L"] = "100000";
	optionSingleArgs["-combine_t"] = "10000";
	optionSingleArgs["-combine_h"] = "10000";
	optionSingleArgs["-combine_s"] = "10";
	optionSingleArgs["-combine_g"] = "100000";
	optionSingleArgs["-combine_i"] = "0.9";

    optionBool["-no_scaffold"] = false;
    optionBool["-unphase"] = false;
	optionBool["-combine"] = false;
    optionBool["-reduce_redundancy"] = false;
    optionBool["-trim_overlap"] = false;
	optionBool["-divide_only"] = false;

    pairedEndSingleFileType.push_back("-ip");
    pairedEndSingleFileType.push_back("-op");
    pairedEndPairFileType.push_back("-IP");
    pairedEndPairFileType.push_back("-OP");
    optionSingleArgs["-tmp"] = ".";

    optionSingleArgs["-mapper"] = platanus::ConstParam::SUB_BIN_PATH + "/minimap2";
    optionBool["-minimap2_sensitive"] = false;
    optionBool["-minimap2"] = true;
    optionBool["-minialign"] = false;
}


//////////////////////////////////////////////////////////////////////////////////////
// usage
//////////////////////////////////////////////////////////////////////////////////////
void SolveDBG::usage(void) const
{
    std::cerr << "\nUsage: platanus_b solveDBG [Options]\n"
              << "Options:\n"
              << "    -o STR                             : prefix of output file (default " << optionSingleArgs.at("-o") << ", length <= " << platanus::ConstParam::MAX_FILE_LEN << ")\n"
              << "    -c FILE1 [FILE2 ...]               : contig_file (fasta format)\n"
              << "    -masked FILE1 [FILE2 ...]          : masked_contig_file for long-read mapping (fasta format)\n"
//              << "    -b FILE1 [FILE2 ...]               : bubble_seq_file (fasta format)\n"
              << "    -ip{INT} PAIR1 [PAIR2 ...]         : lib_id inward_pair_file (reads in 1 file, fasta or fastq)\n"
              << "    -IP{INT} FWD1 REV1 [FWD2 REV2 ...] : lib_id inward_pair_files (reads in 2 files, fasta or fastq)\n"
              << "    -op{INT} PAIR1 [PAIR2 ...]         : lib_id outward_pair_file (reads in 1 file, fasta or fastq)\n"
              << "    -OP{INT} FWD1 REV1 [FWD2 REV2 ...] : lib_id outward_pair_files (reads in 2 files, fasta or fastq)\n"
              << "    -p FILE1 [FILE2 ...]               : PacBio long-read file (fasta or fastq)\n"
              << "    -ont FILE1 [FILE2 ...]             : Oxford Nanopore long-read file (fasta or fastq)\n"
              << "    -gc FILE1 [FILE2 ...]              : Guiding contig file; i.e. other assemblies, synthetic long-reads or corrected reads (fasta or fastq)\n"
              << "    -x PAIR1 [PAIR2 ...]               : tagged_pair_files (10x Genomics) (reads in 1 file, fasta or fastq)\n"
              << "    -X FWD1 REV1 [FWD2 REV2 ...]       : tagged_pair_files (10x Genomics) (reads in 2 files, fasta or fastq)\n"
              << "    -n{INT} INT                        : lib_id minimum_insert_size\n"
              << "    -a{INT} INT                        : lib_id average_insert_size\n"
              << "    -d{INT} INT                        : lib_id SD_insert_size\n"
              << "    -e FLOAT                           : coverage depth of homozygous region (default auto)\n"
              << "    -L INT                             : maximum fragment length of tag (10x Genomics) (default " << optionSingleArgs.at("-L") << ")\n"
              << "    -s INT1 [INT2 ...]                 : mapping seed length for short reads (default " << optionMultiArgs.at("-s")[0] << ")\n"
//              << "    -S INT1 [INT2 ...]                 : mapping seed length for long reads (default " << optionMultiArgs.at("-S")[0] << ", only effective with -kmer_align option)\n"
//              << "    -k INT                             : minimum number of links to phase variants (default " << optionSingleArgs.at("-k") << ")\n"
              << "    -l INT                             : minimum number of links to scaffold (default " << optionSingleArgs.at("-l") << ")\n"
              << "    -u FLOAT                           : maximum difference for bubble crush (identity, default " << optionSingleArgs.at("-u") << ")\n"
              << "    -t INT                             : number of threads (<= " << optionSingleArgs.at("-t") << ", default 1)\n"
              << "    -unphase                           : not phase heterozygous regions and construct consensus scaffolds (default false)\n"
              << "    -combine                           : base-assembly (-c) is extended using other assemblies (-gc)  (default false)\n"
              << "    -mapper FILE                       : path of mapper executable file (default, minimap2; only effective with -p -ont -gc option)\n"
              << "    -minimap2_sensitive                : sensitive mode for minimap2 (default, off; only effective with -p -ont -gc option)\n"
              << "    -reduce_redundancy                 : reduce redundant sequences that exactly matche others (default, off)\n"
              << "    -trim_overlap                      : trim overlapping edges of scaffolds (default, off)\n"
              << "    -divide_only                       : only divide input sequences (default, off)\n"
//              << "    -minimap2                          : use minimap2 insterd of minimap (default) for alignment of long reads\n"
//              << "    -minialign                         : use minialign insterd of minimap (default) for alignment of long reads\n"
              << "    -tmp DIR                           : directory for temporary files (default " << optionSingleArgs.at("-tmp") << ")\n\n\n"

              << "Outputs:\n"
			  << "    PREFIX_*.fa\n"
             << std::endl;
}


//////////////////////////////////////////////////////////////////////////////////////
// initialize parameters
//////////////////////////////////////////////////////////////////////////////////////
void SolveDBG::initializeParameters(void)
{
	seedLength = platanus::ConstParam::SCAFFOLD_HASH_OVERLAP;
    multiSeedLengthForShortRead.clear();
	if (!(optionPairFile.empty())) {
		for (auto itr = optionMultiArgs["-s"].begin(); itr != optionMultiArgs["-s"].end(); ++itr) {
			unsigned length = std::stoi(*itr);
			multiSeedLengthForShortRead.push_back(length);
			if (seedLength > length)
				seedLength = length;
		}
	}

    keyLength = std::min(seedLength, platanus::ConstParam::SCAFFOLD_HASH_OVERLAP);
    bubbleThreshold = atof(optionSingleArgs["-u"].c_str());
    minLink = atoi(optionSingleArgs["-l"].c_str());
    minLinkToPhase = atoi(optionSingleArgs["-k"].c_str());
    minOverlapForScaffolding = atoi(optionSingleArgs["-v"].c_str());
    numThread = atoi(optionSingleArgs["-t"].c_str());
    pairedDBG.setSeedLength(seedLength);
    pairedDBG.setMinTolerenceFactor(MIN_TOL_FACTOR);
    pairedDBG.setMaxFragmentLengthOfTag(atoi(optionSingleArgs["-L"].c_str()));

    sort(optionPairFile.begin(), optionPairFile.end());
    numFilePerLibraryID.resize(optionPairFile.size());
    libraryIDList.resize(optionPairFile.size());
    unsigned numLibrary = 0;
    for (unsigned i = 0; i < optionPairFile.size(); ++i) {
        ++(numFilePerLibraryID[numLibrary]);
        libraryIDList[numLibrary] = optionPairFile[i].libraryID;
		if (i + 1 >= optionPairFile.size() || optionPairFile[i].libraryID != optionPairFile[i + 1].libraryID) {
            ++numLibrary;
        }
    }

    libraryMT.resize(numLibrary);
    omp_set_num_threads(numThread);

	platanus::setGlobalTmpFileDir(optionSingleArgs["-tmp"].c_str());
}


//////////////////////////////////////////////////////////////////////////////////////
// exec scaffold
//////////////////////////////////////////////////////////////////////////////////////
void SolveDBG::exec(void)
{
    initializeParameters();
    mapLibraryAndInitGraph(numThread);

	if (optionBool["-unphase"]) {
		pairedDBG.clearEdges();

		extendConsensus(false, true);

		if (!(libraryMT.empty()))
			pairedDBG.setTolerence(MAX_TOL_FACTOR * libraryMT.back()[0].getSDInsSize());
		else
			pairedDBG.setTolerence(this->contigMaxK);

		pairedDBG.loadResultSeq(this->contigMaxK, this->contigReadLength, optionBool["-trim_overlap"]);

		if (optionBool["-reduce_redundancy"])
			pairedDBG.markRedundantResultSeq(numThread);

		outputGraph("_consensusScaffold.fa");
		pairedDBG.outputResultSeqComponent(optionSingleArgs["-o"], "_consensusScaffoldComponent.bed");
		cerr << "solve_DBG completed!" << endl;
		return;
	}

	if (optionBool["-combine"]) {
		pairedDBG.clearEdges();

		combineAssembly();

		pairedDBG.loadResultSeq(this->contigMaxK, this->contigReadLength, optionBool["-trim_overlap"]);

		if (optionBool["-reduce_redundancy"])
			pairedDBG.markRedundantResultSeq(numThread);

		outputGraph("_combined.fa");
		pairedDBG.outputResultSeqComponent(optionSingleArgs["-o"], "_combinedComponent.bed");
		cerr << "solve_DBG completed!" << endl;
		return;
	}


	if (optionBool["-divide_only"]) {
		pairedDBG.clearEdges();

		extendConsensusToEstimateInsertSize();
		pairedDBG.resetGraph();

		pairedDBG.divideErroneousNodeBaseLevel(0, libraryMT.size(), this->contigReadLength, false, true, true, numThread);

		pairedDBG.loadDividedContigResultSeq(this->contigMaxK, this->contigReadLength);
		if (optionBool["-reduce_redundancy"])
			pairedDBG.markRedundantResultSeq(numThread);

		outputGraph("_divided.fa");
		pairedDBG.outputResultSeqComponent(optionSingleArgs["-o"], "_dividedComponent.bed");
		cerr << "solve_DBG completed!" << endl;
		return;
	}


	pairedDBG.extractDBGBubbleInformation();
	pairedDBG.clearEdges();

	extendConsensus(true, false);
	pairedDBG.resetGraph();


	pairedDBG.setMode(PairedDBG::OVERLAP_MODE);
	pairedDBG.setCutoffLength(0);
	pairedDBG.makeGraph(numThread);
	pairedDBG.extractDBGBubbleInformation();
	pairedDBG.setOppositeBubbleContigIDByEndMatch();
	pairedDBG.setOppositeBubbleContigIDOverlapped(numThread);
	pairedDBG.setBubbleJunctionContigIDOverlapped();
	pairedDBG.clearEdges();


	for (long outerIteration = 0; outerIteration < 4; ++outerIteration) {
		for (long iteration = 0; iteration < 2; ++iteration) {
			if (iteration == 0)
				pairedDBG.setMode(PairedDBG::OVERLAP_MODE | PairedDBG::PAIRED_END_LINK_MODE);
			else
				pairedDBG.setMode(PairedDBG::OVERLAP_MODE | PairedDBG::PAIRED_END_LINK_MODE | PairedDBG::NON_DBG_OVERLAP_MODE);
			pairedDBG.setCutoffLength(0);

			for (unsigned i = 0; i < libraryMT.size(); ++i) {
				cerr << "[LIBRARY " << i + 1 << "]\nAVE_INS = " << libraryMT[i][0].getAverageInsSize() << ", SD_INS = " << libraryMT[i][0].getSDInsSize() << endl;

				pairedDBG.setTargetLibraryIndex(i);
				unsigned tolerenceFactor = MAX_TOL_FACTOR;
				pairedDBG.setTolerence(tolerenceFactor * libraryMT[i][0].getSDInsSize());
				cerr << "TOLERENCE_OF_CONTIG_OVERLAP = " << pairedDBG.getTolerence() << endl;

				pairedDBG.setMinLink(minLinkToPhase);
				pairedDBG.joinUnambiguousNodePairIterative(numThread);
				pairedDBG.solveSimpleCrossStructureIterative(PairedDBG::CROSS_RESOLUTION_MODE::LINK, numThread);

				pairedDBG.setMinLink(minLink);
				pairedDBG.solveSimpleCrossStructureIterative(PairedDBG::CROSS_RESOLUTION_MODE::LINK, numThread);
			}
			if (longReadLibraryMT.size() > 0) {
				if (iteration == 0)
					pairedDBG.setMode(PairedDBG::OVERLAP_MODE | PairedDBG::LONG_READ_LINK_MODE);
				else
					pairedDBG.setMode(PairedDBG::OVERLAP_MODE | PairedDBG::LONG_READ_LINK_MODE | PairedDBG::NON_DBG_OVERLAP_MODE);
				cerr << "[LONG_READ_LIBRARY]\nAVE_READ_LENGTH = " << longReadLibraryMT[0].getAverageInsSize() << ", SD_READ_LENGTH = " << longReadLibraryMT[0].getSDInsSize() << endl;
				pairedDBG.setMinLink(minLinkToPhase);
				pairedDBG.joinUnambiguousNodePairIterative(numThread);
				pairedDBG.solveSimpleCrossStructureIterative(PairedDBG::CROSS_RESOLUTION_MODE::SCORE, numThread);

				pairedDBG.setMinLink(minLink);
				pairedDBG.solveSimpleCrossStructureIterative(PairedDBG::CROSS_RESOLUTION_MODE::TAG, numThread);
			}

			if (iteration == 0)
				pairedDBG.setMode(PairedDBG::OVERLAP_MODE | PairedDBG::PAIRED_END_LINK_MODE | PairedDBG::LONG_READ_LINK_MODE);
			else
				pairedDBG.setMode(PairedDBG::OVERLAP_MODE | PairedDBG::PAIRED_END_LINK_MODE | PairedDBG::LONG_READ_LINK_MODE | PairedDBG::NON_DBG_OVERLAP_MODE);
			for (unsigned i = libraryMT.size() - 1; i < libraryMT.size(); ++i) {
				pairedDBG.setTargetLibraryIndex(i);
				pairedDBG.setTolerenceFactor(MAX_TOL_FACTOR);
				pairedDBG.setMinLink(minLinkToPhase);
				pairedDBG.makeGraphAllLibraries(numThread);
				pairedDBG.joinUnambiguousNodePairIterative(numThread);
				pairedDBG.solveSimpleCrossStructureAllLibrariesIterative(PairedDBG::CROSS_RESOLUTION_MODE::LINK, numThread);

				pairedDBG.setMinLink(minLink);
				pairedDBG.solveSimpleCrossStructureAllLibrariesIterative(PairedDBG::CROSS_RESOLUTION_MODE::TAG, numThread);
			}
			if (optionBool["-no_scaffold"]) {
				pairedDBG.divideErroneousNode(minLink, PairedDBG::DIVISION_MODE::SWITCH, numThread);
				pairedDBG.divideNodeBasedOnBubblesIterative(true, numThread);
			}
		}

		if (optionBool["-no_scaffold"]) {
			pairedDBG.divideNodeBasedOnBubblesIterative(true, numThread);
			pairedDBG.divideErroneousNode(minLink, PairedDBG::DIVISION_MODE::SWITCH, numThread);
			pairedDBG.divideNodeBasedOnBubblesIterative(true, numThread);
			continue;
		}

		pairedDBG.setMode(PairedDBG::OVERLAP_MODE);
		pairedDBG.setCutoffLength(0);
		pairedDBG.setMinLink(minLink);
		pairedDBG.makeGraph(numThread);
		pairedDBG.clearEdges();


		for (long iteration = 0; iteration < 2; ++iteration) {
			for (unsigned i = 0; i < libraryMT.size(); ++i) {
				pairedDBG.setTargetLibraryIndex(i);
				if (iteration == 0)
					pairedDBG.setMinLink(minLinkToPhase);
				else
					pairedDBG.setMinLink(minLink);

				cerr << "[LIBRARY " << i + 1 << "]\nAVE_INS = " << libraryMT[i][0].getAverageInsSize() << ", SD_INS = " << libraryMT[i][0].getSDInsSize() << endl;

				pairedDBG.setMode(PairedDBG::PAIRED_END_LINK_MODE | PairedDBG::LENGTH_CUTOFF_MODE | PairedDBG::BUBBLE_AWARE_MODE);

				for (unsigned tolerenceFactor = MIN_TOL_FACTOR; tolerenceFactor <= MAX_TOL_FACTOR; ++tolerenceFactor) {
					pairedDBG.setCutoffLength(2.0 * tolerenceFactor * std::min(1.0 * libraryMT[i][0].getSDInsSize(), 0.1 * libraryMT[i][0].getAverageInsSize()));
					pairedDBG.setTolerence(pairedDBG.getCutoffLength() / 2);
					cerr << "TOLERENCE_LEVEL_OF_CONTIG_OVERLAP = " << pairedDBG.getTolerence() << endl;
					cerr << "LENGTH_CUTOFF = " << pairedDBG.getCutoffLength() << endl;

					pairedDBG.trimSparseEnd();
					if (iteration > 0)
						pairedDBG.joinUnambiguousNodePairGappedIterative(numThread);
					pairedDBG.solveSimpleGappedCrossStructureIterative(PairedDBG::CROSS_RESOLUTION_MODE::LINK, numThread);
					pairedDBG.solveSimpleGappedCrossStructureIterative(PairedDBG::CROSS_RESOLUTION_MODE::TAG, numThread);
				}
			}
			if (longReadLibraryMT.size() > 0) {
				pairedDBG.setMode(PairedDBG::LONG_READ_LINK_MODE);
				cerr << "[LONG_READ_LIBRARY]\nAVE_READ_LENGTH = " << longReadLibraryMT[0].getAverageLength() << endl;
				pairedDBG.setTolerence(2 * this->contigMaxK);
				if (iteration == 0)
					pairedDBG.setMinLink(minLinkToPhase);
				else
					pairedDBG.setMinLink(minLink);

				for (unsigned cutoffFactor = MIN_LONG_READ_LENGTH_CUTOFF_FACTOR; cutoffFactor <= MAX_LONG_READ_LENGTH_CUTOFF_FACTOR; ++cutoffFactor) {
					pairedDBG.setCutoffLength(cutoffFactor * LONG_READ_LENGTH_CUTOFF_UNIT_FACTOR * longReadLibraryMT[0].getAverageInsSize());
					pairedDBG.trimSparseEnd();
					pairedDBG.solveSimpleGappedCrossStructureIterative(PairedDBG::CROSS_RESOLUTION_MODE::SCORE, numThread);
					pairedDBG.solveSimpleGappedCrossStructureIterative(PairedDBG::CROSS_RESOLUTION_MODE::TAG, numThread);
				}
			}
			for (unsigned i = libraryMT.size() - 1; i < libraryMT.size(); ++i) {
				pairedDBG.setTargetLibraryIndex(i);
				pairedDBG.setMode(PairedDBG::PAIRED_END_LINK_MODE | PairedDBG::LONG_READ_LINK_MODE | PairedDBG::LENGTH_CUTOFF_MODE | PairedDBG::BUBBLE_AWARE_MODE);

				for (unsigned tolerenceFactor = MIN_TOL_FACTOR; tolerenceFactor <= MAX_TOL_FACTOR; ++tolerenceFactor) {
					pairedDBG.setTolerenceFactor(tolerenceFactor);
					pairedDBG.setCutoffLengthFactor(MIN_LONG_READ_LENGTH_CUTOFF_FACTOR * LONG_READ_LENGTH_CUTOFF_UNIT_FACTOR);
					pairedDBG.setCutoffLength(2.0 * tolerenceFactor * std::min(1.0 * libraryMT[i][0].getSDInsSize(), 0.1 * libraryMT[i][0].getAverageInsSize()));
					pairedDBG.setTolerence(pairedDBG.getCutoffLength() / 2);
					if (iteration == 0)
						pairedDBG.setMinLink(minLinkToPhase);
					else {
						pairedDBG.setMinLink(minLink);
						pairedDBG.joinUnambiguousNodePairGappedIterativeAllLibraries(numThread);
					}
					pairedDBG.trimSparseEnd();
					pairedDBG.solveSimpleGappedCrossStructureAllLibrariesIterative(PairedDBG::CROSS_RESOLUTION_MODE::LINK, numThread);
					pairedDBG.solveSimpleGappedCrossStructureAllLibrariesIterative(PairedDBG::CROSS_RESOLUTION_MODE::TAG, numThread);
				}
			}
		}


		pairedDBG.trimSparseEnd();
		pairedDBG.divideNodeUsingBubbleContigPair(numThread);

		pairedDBG.divideErroneousNode(minLink, PairedDBG::DIVISION_MODE::SWITCH, numThread);
		pairedDBG.divideNodeBasedOnBubblesIterative(true, numThread);
		pairedDBG.divideBubbleContigInNonHeteroNode();
		pairedDBG.divideBubbleJunctionNode(false);

		pairedDBG.setMinOverlap(minOverlapForScaffolding);
		for (long iteration = 0; iteration < 2; ++iteration) {
			pairedDBG.setMode(PairedDBG::PAIRED_END_LINK_MODE | PairedDBG::LENGTH_CUTOFF_MODE | PairedDBG::BUBBLE_AWARE_MODE | PairedDBG::PREVIOUS_DIVISION_AWARE_MODE);

			for (unsigned i = 0; i < libraryMT.size(); ++i) {
				pairedDBG.setTargetLibraryIndex(i);
				long linkThreshold;
				if (iteration % 2 == 0)
					linkThreshold = minLink;
				else
					linkThreshold = std::max(minLink, pairedDBG.estimateLink());

				cerr << "[LIBRARY " << i + 1 << "]\nAVE_INS = " << libraryMT[i][0].getAverageInsSize() << ", SD_INS = " << libraryMT[i][0].getSDInsSize() << endl;

				for (unsigned tolerenceFactor = MIN_TOL_FACTOR; tolerenceFactor <= MAX_TOL_FACTOR; ++tolerenceFactor) {
					pairedDBG.setCutoffLength(2.0 * tolerenceFactor * std::min(1.0 * libraryMT[i][0].getSDInsSize(), 0.1 * libraryMT[i][0].getAverageInsSize()));
					pairedDBG.setTolerence(pairedDBG.getCutoffLength() / 2);
					cerr << "TOLERENCE_LEVEL_OF_CONTIG_OVERLAP = " << pairedDBG.getTolerence() << endl;
					cerr << "LENGTH_CUTOFF = " << pairedDBG.getCutoffLength() << endl;

					pairedDBG.trimSparseEnd();

					pairedDBG.setMinLink(linkThreshold);
					pairedDBG.makeGraph(numThread);
					pairedDBG.deleteRepeatEdge();
					pairedDBG.detectRepeat(pairedDBG.getHeteroCoverage());
					pairedDBG.makeScaffold();

					pairedDBG.joinUnambiguousNodePairGappedIterative(numThread);
					pairedDBG.setMinLink(minLink);
					pairedDBG.solveSimpleGappedCrossStructureIterative(PairedDBG::CROSS_RESOLUTION_MODE::LINK, numThread);
					pairedDBG.solveSimpleGappedCrossStructureIterative(PairedDBG::CROSS_RESOLUTION_MODE::TAG, numThread);

					pairedDBG.setMinLink(linkThreshold);
					pairedDBG.makeGraph(numThread);
					pairedDBG.deleteErroneousEdgeNumTagRateIterative(numThread);
					pairedDBG.deleteDifferentBubbleEdgeIterative(numThread);
					pairedDBG.deleteRepeatEdge();
					pairedDBG.detectRepeat(pairedDBG.getHeteroCoverage());
					pairedDBG.makeScaffold();
				}
			}
			if (outerIteration < 2)
				pairedDBG.divideErroneousNode(minLink, PairedDBG::DIVISION_MODE::GAP, numThread, this->contigMaxK);
			else
				pairedDBG.divideErroneousNode(minLink, PairedDBG::DIVISION_MODE::SWITCH, numThread);
		}
		pairedDBG.setMinOverlap(this->contigMaxK - 1);


		pairedDBG.trimSparseEnd();
		pairedDBG.divideNodeUsingBubbleContigPair(numThread);
		pairedDBG.divideErroneousNode(minLink, PairedDBG::DIVISION_MODE::SWITCH, numThread);
		pairedDBG.divideNodeBasedOnBubblesIterative(false, numThread);


		pairedDBG.setMinOverlap(minOverlapForScaffolding);
		for (long iteration = 0; iteration < 2; ++iteration) {
			long linkThreashold = (iteration + 1) * minLink;
			pairedDBG.setMode(PairedDBG::PAIRED_END_LINK_MODE | PairedDBG::LONG_READ_LINK_MODE | PairedDBG::LENGTH_CUTOFF_MODE | PairedDBG::BUBBLE_AWARE_MODE | PairedDBG::PREVIOUS_DIVISION_AWARE_MODE);

			for (unsigned i = libraryMT.size() - 1; i < libraryMT.size(); ++i) {
				pairedDBG.setTargetLibraryIndex(i);
				for (unsigned tolerenceFactor = MIN_TOL_FACTOR; tolerenceFactor <= MAX_TOL_FACTOR; ++tolerenceFactor) {
					pairedDBG.setTolerenceFactor(tolerenceFactor);
					pairedDBG.setCutoffLengthFactor(MAX_LONG_READ_LENGTH_CUTOFF_FACTOR * LONG_READ_LENGTH_CUTOFF_UNIT_FACTOR);

					pairedDBG.trimSparseEnd();

					pairedDBG.setMinLink(linkThreashold);
					pairedDBG.setCutoffLength(2.0 * tolerenceFactor * std::min(1.0 * libraryMT[i][0].getSDInsSize(), 0.1 * libraryMT[i][0].getAverageInsSize()));
					pairedDBG.setTolerence(pairedDBG.getCutoffLength() / 2);
					pairedDBG.makeGraphAllLibraries(numThread);
					pairedDBG.deleteLongEdge(libraryMT[i][0].getAverageInsSize());
					pairedDBG.deleteRepeatEdge();
					pairedDBG.detectRepeat(pairedDBG.getHeteroCoverage());
					pairedDBG.makeScaffold();

					pairedDBG.setMinLink(minLink);
					pairedDBG.joinUnambiguousNodePairGappedIterativeAllLibraries(numThread);
					pairedDBG.deleteLongEdge(libraryMT[i][0].getAverageInsSize());
					pairedDBG.solveSimpleGappedCrossStructureAllLibrariesIterative(PairedDBG::CROSS_RESOLUTION_MODE::LINK, numThread);
					pairedDBG.solveSimpleGappedCrossStructureAllLibrariesIterative(PairedDBG::CROSS_RESOLUTION_MODE::TAG, numThread);

					pairedDBG.setMinLink(linkThreashold);
					pairedDBG.setCutoffLength(2.0 * tolerenceFactor * std::min(1.0 * libraryMT[i][0].getSDInsSize(), 0.1 * libraryMT[i][0].getAverageInsSize()));
					pairedDBG.setTolerence(pairedDBG.getCutoffLength() / 2);
					pairedDBG.makeGraphAllLibraries(numThread);
					pairedDBG.deleteLongEdge(libraryMT[i][0].getAverageInsSize());
					pairedDBG.deleteErroneousEdgeNumTagRateIterative(numThread);
					pairedDBG.deleteDifferentBubbleEdgeIterative(numThread);
					pairedDBG.deleteConflictingBubbleEdge(numThread);
					pairedDBG.deleteRepeatEdge();
					pairedDBG.detectRepeat(pairedDBG.getHeteroCoverage());
					pairedDBG.makeScaffold();
				}
			}
			if (outerIteration < 2)
				pairedDBG.divideErroneousNode(minLink, PairedDBG::DIVISION_MODE::GAP, numThread, this->contigMaxK);
			else
				pairedDBG.divideErroneousNode(minLink, PairedDBG::DIVISION_MODE::SWITCH, numThread);
		}
		pairedDBG.setMinOverlap(this->contigMaxK - 1);


		if (longReadLibraryMT.size() > 0) {
			pairedDBG.setMinOverlap(minOverlapForScaffolding);

			pairedDBG.divideNodeUsingBubbleContigPair(numThread);

			for (unsigned cutoffFactor = MIN_LONG_READ_LENGTH_CUTOFF_FACTOR; cutoffFactor <= MAX_LONG_READ_LENGTH_CUTOFF_FACTOR; ++cutoffFactor) {
				if (outerIteration == 0)
					pairedDBG.setMode(PairedDBG::LONG_READ_LINK_MODE | PairedDBG::LENGTH_CUTOFF_MODE | PairedDBG::PREVIOUS_DIVISION_AWARE_MODE);
				else
					pairedDBG.setMode(PairedDBG::LONG_READ_LINK_MODE | PairedDBG::LENGTH_CUTOFF_MODE | PairedDBG::BUBBLE_AWARE_MODE | PairedDBG::PREVIOUS_DIVISION_AWARE_MODE);

				pairedDBG.setTolerence(2 * this->contigMaxK);
				pairedDBG.setCutoffLength(cutoffFactor * LONG_READ_LENGTH_CUTOFF_UNIT_FACTOR *  longReadLibraryMT[0].getAverageInsSize());

				pairedDBG.trimSparseEnd();

				pairedDBG.setMinLink(minLink);
				pairedDBG.makeGraph(numThread);
				pairedDBG.deleteErroneousEdgeNumLinkRateIterative(numThread);
				pairedDBG.deleteRepeatEdge();
				pairedDBG.detectRepeat(pairedDBG.getHeteroCoverage());
				pairedDBG.makeScaffold();

				pairedDBG.setMinLink(minLink);
				pairedDBG.solveSimpleGappedCrossStructureIterative(PairedDBG::CROSS_RESOLUTION_MODE::LINK, numThread);
				pairedDBG.solveSimpleGappedCrossStructureIterative(PairedDBG::CROSS_RESOLUTION_MODE::TAG, numThread);

				pairedDBG.setMinLink(minLink);
				pairedDBG.makeGraph(numThread);
				pairedDBG.deleteErroneousEdgeNumTagRateIterative(numThread);
				pairedDBG.deleteDifferentBubbleEdgeIterative(numThread);
				pairedDBG.deleteErroneousEdgeNumLinkRateIterative(numThread);
				pairedDBG.deleteRepeatEdge();
				pairedDBG.detectRepeat(pairedDBG.getHeteroCoverage());
				pairedDBG.makeScaffold();

				pairedDBG.divideNodeUsingBubbleContigPair(numThread);
				pairedDBG.setMode(PairedDBG::OVERLAP_MODE | PairedDBG::LONG_READ_LINK_MODE);
				pairedDBG.setMinLink(minLinkToPhase);
				pairedDBG.setCutoffLength(0);
				pairedDBG.joinUnambiguousNodePairIterative(numThread);
				pairedDBG.solveSimpleCrossStructureIterative(PairedDBG::CROSS_RESOLUTION_MODE::LINK, numThread);
				pairedDBG.solveSimpleCrossStructureIterative(PairedDBG::CROSS_RESOLUTION_MODE::TAG, numThread);
			}

			pairedDBG.setMinOverlap(this->contigMaxK - 1);
		}


		pairedDBG.trimSparseEnd();
		pairedDBG.divideNodeBasedOnBubblesIterative(true, numThread);
		pairedDBG.divideBubbleJunctionNode(true);
		pairedDBG.divideErroneousNode(minLink, PairedDBG::DIVISION_MODE::SWITCH, numThread);
		pairedDBG.divideNodeBasedOnBubblesIterative(true, numThread);

		if (outerIteration < 2) {
			pairedDBG.joinUnambiguousNodePairIterative(numThread);
			pairedDBG.divideGappedNode((outerIteration + 1) * this->contigMaxK);
		}
	}

	pairedDBG.divideBubbleContigInNonHeteroNode();

	pairedDBG.setMode(0);
	pairedDBG.makeGraph(numThread);
	pairedDBG.adjustOppositeBubbleNodeIDDirection();
	pairedDBG.copyAllNodes(phasedGraph);

	pairedDBG.deleteSecondaryBubbleNodeAndEdge(numThread);
	pairedDBG.clearEdges();
	pairedDBG.makeScaffold();

	extendConsensus(false, true);

	if (!(libraryMT.empty()))
		pairedDBG.setTolerence(MAX_TOL_FACTOR * libraryMT.back()[0].getSDInsSize());
	else
		pairedDBG.setTolerence(this->contigMaxK);
    pairedDBG.loadResultSeq(this->contigMaxK, this->contigReadLength, false);
	outputGraph("_preliminaryConsensusScaffold.fa");

	pairedDBG.setMode(0);
	pairedDBG.remakeGraphRecoveringSecondaryBubble(phasedGraph);
	pairedDBG.makeGraph(numThread);
	pairedDBG.divideNodeBasedOnBubblesIterative(false, numThread);
	pairedDBG.makeGraph(numThread);
	pairedDBG.adjustOppositeBubbleNodeIDDirection();
	if (!(libraryMT.empty()))
		pairedDBG.setTolerence(MAX_TOL_FACTOR * libraryMT.back()[0].getSDInsSize());
	else
		pairedDBG.setTolerence(this->contigMaxK);
    pairedDBG.loadResultSeq(this->contigMaxK, this->contigReadLength, optionBool["-trim_overlap"]);

	if (optionBool["-reduce_redundancy"])
		pairedDBG.markRedundantResultSeq(numThread);

	if (!optionBool["-no_scaffold"])
		pairedDBG.outputResultSeqWithBubble(optionSingleArgs["-o"], "_primaryBubble.fa", "_secondaryBubble.fa", "_nonBubbleHetero.fa", "_nonBubbleOther.fa", "_bubbleRelation.tsv", this->contigMaxK, this->contigReadLength);
	else
		pairedDBG.outputResultSeqWithBubble(optionSingleArgs["-o"], "_primaryBubbleContig.fa", "_secondaryBubbleContig.fa", "_nonBubbleHeteroContig.fa", "_nonBubbleOtherContig.fa", "_bubbleContigRelation.tsv", this->contigMaxK, this->contigReadLength);

	pairedDBG.outputResultSeqComponent(optionSingleArgs["-o"], "_phasedScaffoldComponent.bed");
	
	cerr << "solve_DBG completed!" << endl;
	return;
}


//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
void SolveDBG::mapLibraryAndInitGraph(const int numThread)
{
    platanus::Contig contig;

    std::unique_ptr<HeteroMapper> mapper(new HeteroMapper(seedLength, keyLength));

    readLibrary(mapper, contig, numThread);
    cerr << "CONTIG_AVERAGE_COVERAGE = " << averageCoverage << endl;

	mapper->setMultiSeedLength(multiSeedLengthForShortRead);
    unsigned nowFileNumber = 0;
    for (unsigned i = 0; i < libraryMT.size(); ++i) {
		libraryMT[i][0].setAverageInsSize(0);
        int nowLibraryID = optionPairFile[nowFileNumber].libraryID;
        cerr << "[LIBRARY " << libraryIDList[i] << "]" << endl;
        // set average length and minimum insert size
        // estimate insert size
        libraryMT[i][0].setAverageLength((long)((double)(libraryMT[i][0].getTotalLength()) / (2 * libraryMT[i][0].getNumPair()) + 0.5));
        long minInsertion = optionMinIns.find(nowLibraryID) == optionMinIns.end() ? 0 : optionMinIns[nowLibraryID];

        mapper->contigMap.mapPairAndSaveReadLink(libraryMT[i], minInsertion, this->contigMaxK + 1, numThread);
//        mapper->contigMap.mapPairAndSaveReadLink(libraryMT[i], minInsertion, 0, numThread);

        libraryMT[i][0].setInsCutoffRate(optionInsCutoffRate.find(nowLibraryID) == optionInsCutoffRate.end() ? DEFAULT_INS_CUTOFF_RATE : optionInsCutoffRate[nowLibraryID]);
        if (optionAveIns.find(nowLibraryID) != optionAveIns.end() || optionSDIns.find(nowLibraryID) != optionSDIns.end()) {
            if (optionAveIns.find(nowLibraryID) != optionAveIns.end()) {
                libraryMT[i][0].setAverageInsSize(optionAveIns[nowLibraryID]);
                std::cerr << "Average insert size specified: AVE = " << libraryMT[i][0].getAverageInsSize() << std::endl;
            }
            if (optionSDIns.find(nowLibraryID) != optionSDIns.end()) {
                libraryMT[i][0].setSDInsSize(optionSDIns[nowLibraryID]);
            } else {
                libraryMT[i][0].setSDInsSize(static_cast<long>(static_cast<double>(libraryMT[i][0].getAverageInsSize()) / 10.0 + 0.5));
            }
        }
        nowFileNumber += numFilePerLibraryID[i];
    }

	if (longReadLibraryMT.size() > 0) {
		cerr << "[LONG_READ LIBRARY]" << endl;

		string alignerOutFilename(optionSingleArgs["-o"]);
		alignerOutFilename += "_longReadAlignment.tsv";

		string aligner = optionSingleArgs["-mapper"];
		if (optionSingleArgs["-mapper"].empty())
			aligner = "minimap2";

		vector<string> targetFilename;
		if (optionMultiArgs["-masked"].empty())
			targetFilename = optionMultiArgs["-c"];
		else
			targetFilename = optionMultiArgs["-masked"];

		unlink(alignerOutFilename.c_str());
			
		execMinimap2(targetFilename, optionMultiArgs["-b"], pacBioLongReadFilename, alignerOutFilename, numThread, aligner, "-p 1 -x map-pb");
		execMinimap2(targetFilename, optionMultiArgs["-b"], nanoporeLongReadFilename, alignerOutFilename, numThread, aligner, "-p 1 -x map-ont");
		execMinimap2(targetFilename, optionMultiArgs["-b"], guideContigLongReadFilename, alignerOutFilename, numThread, aligner, "-p 1 -x asm10");

		longReadLibraryMT[0].setAverageLength((long)((double)(longReadLibraryMT[0].getTotalLength()) / (2 * longReadLibraryMT[0].getNumPair()) + 0.5));

		if (!optionBool["-combine"])
			mapper->contigMap.readLongReadPAFfileAndSaveLink(alignerOutFilename, longReadLibraryMT, LONG_READ_MIN_ALIGNMENT_LENGTH, LONG_READ_MIN_ALIGNMENT_COVERAGE, LONG_READ_MIN_ALIGNMENT_IDENTITY, contigMaxK, -1, numThread);
		else
			mapper->contigMap.readLongReadPAFfileAndSaveLink(alignerOutFilename, longReadLibraryMT, atol(optionSingleArgs["-combine_l"].c_str()), LONG_READ_MIN_ALIGNMENT_COVERAGE, atof(optionSingleArgs["-combine_i"].c_str()), atol(optionSingleArgs["-combine_t"].c_str()), atol(optionSingleArgs["-combine_h"].c_str()), numThread);

//		std::ostringstream cmdSS;
//		cmdSS << "rm " << alignerOutFilename;
//		system(cmdSS.str().c_str());

		vector<long> insSizeDistribution;
		longReadLibraryMT[0].readInsertSizeFile(insSizeDistribution);
		longReadLibraryMT[0].estimateInsSize(insSizeDistribution, 0, platanus::ConstParam::LONG_READ_INS_SIZE_LOWER_BOUND_FACTOR, platanus::ConstParam::LONG_READ_INS_SIZE_UPPER_BOUND_FACTOR);

//		std::ostringstream outStream;
//		outStream << optionSingleArgs["-o"] << "_longReadLibrary" << "_readDistribution.tsv";
//		longReadLibraryMT[0].printInsertSizeFreq(outStream.str());
		cerr << "[LONG_READ_LIBRARY " << 1 << "]\nAVE_READ_LENGTH = " << longReadLibraryMT[0].getAverageInsSize()
			 << ", SD_READ_LENGTH = " << longReadLibraryMT[0].getSDInsSize() << endl;
	}	

	if (tagLibraryMT.size() > 0) {
		mapper->setMultiSeedLength(multiSeedLengthForShortRead);

		cerr << "[TAG LIBRARY]" << endl;
		mapper->contigMap.mapTagPairMT(tagLibraryMT, numThread);
	}	

	if (libraryMT.size() > 0) {
		pairedDBG.setAllLibraryMT(&libraryMT);
		pairedDBG.setTargetLibraryIndex(0);
		pairedDBG.initScaffolding(contig.coverage, mapper->contigMap, averageCoverage, bubbleThreshold);
		pairedDBG.setContigName(contig.name);

		vector<long> insSizeDistribution;
		libraryMT[0][0].readInsertSizeFile(insSizeDistribution);
		pairedDBG.insertSizeDistribution(libraryMT[0], insSizeDistribution, numThread);

		vector<long> seqLengths;
		pairedDBG.scaffoldLengthList(seqLengths);

		if (libraryMT[0][0].getAverageInsSize() == 0)
			libraryMT[0][0].estimateInsSize(insSizeDistribution, 0, platanus::ConstParam::SHORT_READ_INS_SIZE_LOWER_BOUND_FACTOR, platanus::ConstParam::SHORT_READ_INS_SIZE_UPPER_BOUND_FACTOR);

		std::ostringstream outStream;
		outStream << optionSingleArgs["-o"] << "_lib" << 1 << "_insFreq.tsv";
		libraryMT[0][0].printInsertSizeFreq(outStream.str());
		cerr << "[LIBRARY " << 1 << "]\nAVE_INS = " << libraryMT[0][0].getAverageInsSize()
			 << ", SD_INS = " << libraryMT[0][0].getSDInsSize() << endl;
	}
	else {
		pairedDBG.initScaffolding(contig.coverage, mapper->contigMap, averageCoverage, bubbleThreshold);
		pairedDBG.setContigName(contig.name);
	}

    pairedDBG.setContigMaxK(this->contigMaxK);
//    pairedDBG.setMinOverlap(this->contigMaxK - 1);
	pairedDBG.setMinOverlap(minOverlapForScaffolding);
    pairedDBG.saveOverlap(mapper->contigMap, this->contigMaxK - 1, this->contigMaxK, numThread);
    pairedDBG.classifyNode();

	if (longReadLibraryMT.size() > 0) {
		pairedDBG.setLongReadLibraryMT(&longReadLibraryMT);
	}

	if (tagLibraryMT.size() > 0) {
		pairedDBG.setTagLibraryMT(&tagLibraryMT);
		pairedDBG.countMappedTagForEachContig(numThread);
	}

    cerr << "destructing mapper objects..." << std::endl;
}


void SolveDBG::readLibrary(std::unique_ptr<HeteroMapper> &mapper, platanus::Contig &contig, const int numThread)
{
    omp_set_num_threads(numThread);
    # pragma omp parallel for schedule(static, 1)
    for (int i = -3; i < static_cast<int>(libraryMT.size()); ++i) {
        try {
            // load contig file
            if (i == -3) {
                for (unsigned i = 0; i < optionMultiArgs["-c"].size(); ++i)
                    contig.readFastaCoverage(optionMultiArgs["-c"][i]);

				long j = contig.numSeq;
                for (unsigned i = 0; i < optionMultiArgs["-b"].size(); ++i)
                    contig.readFastaCoverage(optionMultiArgs["-b"][i]);

				pairedDBG.setNumInputBubbleContig(contig.numSeq - j);
				contig.setNameIndex();


                this->contigMaxK = contig.getMaxKFromFastaHeader(optionMultiArgs["-c"][0]);
                this->contigReadLength = contig.getReadLengthFromFastaHeader(optionMultiArgs["-c"][0]);
                if (this->contigReadLength == 0)
                    this->contigReadLength = platanus::ConstParam::DEFAULT_CONTIG_READ_LEN;

				if (optionSingleArgs["-e"] == "")
					averageCoverage = contig.calculateAverageCoverageExcludingOutlier(contig.getSeqLengthMedian());
				else
					averageCoverage = atof(optionSingleArgs["-e"].c_str());

				pairedDBG.setAverageCoverage(averageCoverage);
				pairedDBG.setHeteroCoverage(averageCoverage/2.0);

                mapper->setContigMap(contig);
                mapper->makeKmerTableContigMap();
            }
			else if (i == -2) {
                if (optionMultiArgs["-p"].empty() && optionMultiArgs["-ont"].empty() && optionMultiArgs["-gc"].empty())
					continue;

				vector<string> filenames(optionMultiArgs["-p"]);
				pacBioLongReadFilename = optionMultiArgs["-p"];

				std::copy(optionMultiArgs["-ont"].begin(), optionMultiArgs["-ont"].end(), std::back_inserter(filenames));
				nanoporeLongReadFilename = optionMultiArgs["-ont"];

				std::copy(optionMultiArgs["-gc"].begin(), optionMultiArgs["-gc"].end(), std::back_inserter(filenames));
				guideContigLongReadFilename = optionMultiArgs["-gc"];

				longReadLibraryMT.resize(numThread);
				for (long j = 0; j < numThread; ++j)
					longReadLibraryMT[j].pairFP = platanus::makeTemporaryFile();

				for (unsigned j = 0; j < static_cast<long>(filenames.size()); ++j) {
					platanus::FILETYPE fileFormat = checkFileFormat(filenames[j]);
					bool isFastq;
					switch (fileFormat) {
					case platanus::UNKNOWN:
						throw platanus::FormatError();
						break;
					case platanus::FASTA:
						isFastq = false;
						break;
					default:
						isFastq = true;
						break;
					}
					ReadFastaSingleMT(longReadLibraryMT, filenames[j], numThread, false, isFastq, true);
				}
            }
			else if (i == -1) {
                if (optionMultiArgs["-x"].size() == 0 && optionMultiArgs["-X"].size() == 0)
					continue;

				vector<string> filenames(optionMultiArgs["-x"]);
				filenames.insert(filenames.end(), optionMultiArgs["-X"].begin(), optionMultiArgs["-X"].end());

				std::unordered_map<string, int> tagStringConverter;
				setTagStringConverter(filenames, tagStringConverter);

				tagLibraryMT.resize(numThread);
				for (long j = 0; j < numThread; ++j)
					tagLibraryMT[j].pairFP = platanus::makeTemporaryFile();

				for (unsigned j = 0; j < static_cast<long>(optionMultiArgs["-x"].size()); ++j) {
					platanus::FILETYPE fileFormat = checkFileFormat(optionMultiArgs["-x"][j]);
					bool isFastq;
					switch (fileFormat) {
					case platanus::UNKNOWN:
						throw platanus::FormatError();
						break;
					case platanus::FASTA:
						isFastq = false;
						break;
					default:
						isFastq = true;
						break;
					}
					ReadFastaSingleTaggedMT(tagLibraryMT, optionMultiArgs["-x"][j], numThread, false, isFastq, true, tagStringConverter);
				}

				for (unsigned j = 0; j < static_cast<long>(optionMultiArgs["-X"].size()); j += 2) {
					platanus::FILETYPE fileFormat1 = checkFileFormat(optionMultiArgs["-X"][j]);
					platanus::FILETYPE fileFormat2 = checkFileFormat(optionMultiArgs["-X"][j + 1]);
					if (fileFormat1 != fileFormat2) {
						throw platanus::FormatError("Different file type in paired-file (-X).");
					}
					bool isFastq;
					switch (fileFormat1) {
					case platanus::UNKNOWN:
						throw platanus::FormatError();
						break;
					case platanus::FASTA:
						isFastq = false;
						break;
					default:
						isFastq = true;
						break;
					}
					ReadFastaPairTaggedMT(tagLibraryMT, optionMultiArgs["-X"][j] , optionMultiArgs["-X"][j + 1], numThread, false, isFastq, tagStringConverter);
				}
            }
			else {
                unsigned nowFileNumber = 0;
                libraryMT[i].resize(numThread);
                for (int j = 0; j < i; ++j)
                    nowFileNumber += numFilePerLibraryID[j];

                for (int j = 0; j < numThread; ++j) {
                    libraryMT[i][j].makeTempPairFP();
                }
                for (int j = 0; j < numFilePerLibraryID[i]; ++j) {
                    bool isFastq, isMate;
                    isMate = optionPairFile[nowFileNumber+j].libraryType == "-op"
                           || optionPairFile[nowFileNumber+j].libraryType == "-OP" ? true : false;
                    if (optionPairFile[nowFileNumber+j].libraryType == "-ip" || optionPairFile[nowFileNumber+j].libraryType == "-op") {
                        platanus::FILETYPE fileFormat = checkFileFormat(optionPairFile[nowFileNumber+j].fileFirst);
                        switch (fileFormat) {
                        case platanus::UNKNOWN:
                            throw platanus::FormatError();
                            break;
                        case platanus::FASTA:
                            isFastq = false;
                            break;
                        default:
                            isFastq = true;
                            break;
                        }
                        ReadFastaSingleMT(libraryMT[i], optionPairFile[nowFileNumber+j].fileFirst, numThread, isMate, isFastq);
                    }
					else {
                        platanus::FILETYPE fileFormat1 = checkFileFormat(optionPairFile[nowFileNumber+j].fileFirst);
                        platanus::FILETYPE fileFormat2 = checkFileFormat(optionPairFile[nowFileNumber+j].fileSecond);
                        if (fileFormat1 != fileFormat2) {
                            throw platanus::FormatError("Different file type in paired-file.");
                        }
                        switch (fileFormat1) {
                        case platanus::UNKNOWN:
                            throw platanus::FormatError();
                            break;
                        case platanus::FASTA:
                            isFastq = false;
                            break;
                        default:
                            isFastq = true;
                            break;
                        }
                        ReadFastaPairMT(libraryMT[i], optionPairFile[nowFileNumber+j].fileFirst, optionPairFile[nowFileNumber+j].fileSecond, numThread, isMate, isFastq);
                    }
                }
            }
        } catch (platanus::ErrorBase &e) {
            e.showErrorMessage();
            exit(e.getID());
        }
    }
}


//////////////////////////////////////////////////////////////////////////////////////
// destroy graph
//////////////////////////////////////////////////////////////////////////////////////
void SolveDBG::outputAndAfterTreatment(void)
{
    cerr << "writing scaffold files..." << endl;
    string outFilename(optionSingleArgs["-o"]);
    string componentFilename(outFilename);
    outFilename += "_solvedContig.fa";
    componentFilename += "_solvedContigComponent.tsv";
    pairedDBG.cutAndPrintSeq(this->contigMaxK, this->contigReadLength, outFilename, componentFilename);
}

void SolveDBG::outputGraph(const char *suffix)
{
    cerr << "writing scaffold files..." << endl;
    string outFilename(optionSingleArgs["-o"]);
    outFilename += suffix;
    pairedDBG.printResultSeq(outFilename);
}

void SolveDBG::updateAndWriteInsertSize(const long libraryIndex)
{
	pairedDBG.updateInsertLengthFP(libraryMT[libraryIndex], numThread);
	vector<long> insSizeDistribution;
	libraryMT[libraryIndex][0].readInsertSizeFile(insSizeDistribution);
	pairedDBG.insertSizeDistribution(libraryMT[libraryIndex], insSizeDistribution, numThread);
	if (libraryIndex > 0)
		libraryMT[libraryIndex][0].estimateInsSize(insSizeDistribution, libraryMT[libraryIndex - 1][0].getAverageInsSize(), platanus::ConstParam::SHORT_READ_INS_SIZE_LOWER_BOUND_FACTOR, platanus::ConstParam::SHORT_READ_INS_SIZE_UPPER_BOUND_FACTOR);
	else
		libraryMT[libraryIndex][0].estimateInsSize(insSizeDistribution, 0, platanus::ConstParam::SHORT_READ_INS_SIZE_LOWER_BOUND_FACTOR, platanus::ConstParam::SHORT_READ_INS_SIZE_UPPER_BOUND_FACTOR);

	std::ostringstream outStream;
	outStream << optionSingleArgs["-o"] << "_lib" << (libraryIndex + 1) << "_insFreq.tsv";
	printInsertSizeFreq(outStream.str(), insSizeDistribution);
}

void SolveDBG::execMinialign(const string targetFilename, const string &readFilename, const string &outFilename, const long numThread, const string minialignExecutable)
{
	std::ostringstream oss;

	oss << minialignExecutable <<  " -x pacbio -m 0 -O paf" << " -t " << numThread << " " << targetFilename << " " << readFilename << " >" << outFilename;

	std::cerr << "Executing minialign ..." << std::endl;
	std::cerr << oss.str() << std::endl << std::endl;

	if (system(oss.str().c_str()) != 0) {
		throw platanus::AlignerError();
	}

	cerr << "minialign fineshed." << endl;
}

void SolveDBG::execMinimap(const vector<string> &contigFilenames, const vector<string> &bubbleFilenames, const vector<string> &readFilenames, const string &outFilename, const long minAlignmentLength, const long numThread, const string minimapExecutable)
{
	std::ostringstream oss;

	oss << "cat ";
	for (auto itr = contigFilenames.begin(); itr != contigFilenames.end(); ++itr)
		oss << " " << *itr;
	for (auto itr = bubbleFilenames.begin(); itr != bubbleFilenames.end(); ++itr)
		oss << " " << *itr;

	oss << " | " << minimapExecutable <<  " -t " << numThread <<  " -L " << minAlignmentLength << " - ";

	for (auto itr = readFilenames.begin(); itr != readFilenames.end(); ++itr)
		oss << " " << *itr;
	oss << " | cut -f1-11 ";
	oss << " >" << outFilename;

	std::cerr << "Executing minimap ..." << std::endl;
	std::cerr << oss.str() << std::endl << std::endl << std::endl;

	if (system(oss.str().c_str()) != 0) {
		throw platanus::AlignerError();
	}

	cerr << "minimap fineshed." << endl;
}

void SolveDBG::execMinimap2(const vector<string> &contigFilenames, const vector<string> &bubbleFilenames, const vector<string> &readFilenames, const string &outFilename, const long numThread, const string minimap2Executable, const string minimap2Option)
{
	if (readFilenames.empty())
		return;

	std::ostringstream oss;

	oss << "cat ";
	for (auto itr = contigFilenames.begin(); itr != contigFilenames.end(); ++itr)
		oss << " " << *itr;
	for (auto itr = bubbleFilenames.begin(); itr != bubbleFilenames.end(); ++itr)
		oss << " " << *itr;

	oss << " | " << minimap2Executable << " -c " << " -t " << numThread;
	if (optionBool["-minimap2_sensitive"])
		oss << " -p 0";
	
	oss << " " << minimap2Option << " - ";


	bool bzip2Flag = false;
	char bzip2TempFileName[platanus::globalTmpFileDir.size() + 7 + 1];
	platanus::FILECOMPRESSION format;

	for (auto itr = readFilenames.begin(); itr != readFilenames.end(); ++itr) {
		format = platanus::checkFileCompression(*itr);
		if (format == platanus::FILECOMPRESSION::BZIP2) {
			bzip2Flag = true;
			break;
		}
	}

	if (bzip2Flag) {
		strcpy(bzip2TempFileName, platanus::globalTmpFileDir.c_str());
		strcat(bzip2TempFileName, "/XXXXXX"); 

        int fd = mkstemp(bzip2TempFileName);
        if (fd == -1) {
            throw platanus::TMPError();
		}
        FILE *fp = fdopen(fd, "w+");
		fclose(fp);

		std::ostringstream bzip2Oss;
		bzip2Oss << "bzip2 -cd ";

		for (auto itr = readFilenames.begin(); itr != readFilenames.end(); ++itr) {
			format = platanus::checkFileCompression(*itr);
			if (format == platanus::FILECOMPRESSION::BZIP2)
				bzip2Oss << " " << *itr;
			else
				oss << " " << *itr;
		}

		bzip2Oss << " >"  << bzip2TempFileName;
		if (system(bzip2Oss.str().c_str()) != 0) {
			throw platanus::AlignerError();
		}

		oss << " " << bzip2TempFileName;
	}
	else {
		for (auto itr = readFilenames.begin(); itr != readFilenames.end(); ++itr)
			oss << " " << *itr;
	}

	oss << " | perl -pne \'s/cg:Z:\\S+//\' ";
	oss << " >>" << outFilename;

	std::cerr << "Executing minimap2 ..." << std::endl;
	std::cerr << oss.str() << std::endl << std::endl << std::endl;

	if (system(oss.str().c_str()) != 0) {
		throw platanus::AlignerError();
	}

	if (bzip2Flag) {
        unlink(bzip2TempFileName);
	}

	cerr << "minimap2 finished." << endl;
}


void SolveDBG::extendConsensusToEstimateInsertSize(void)
{
	pairedDBG.clearContigPreviousParentNodeID();

	pairedDBG.setMode(PairedDBG::OVERLAP_MODE);
	pairedDBG.setMinLink(minLink);
	pairedDBG.makeGraph(numThread);
	pairedDBG.joinUnambiguousNodePairIterative(numThread);

	pairedDBG.makeGraph(numThread);
	pairedDBG.setOppositeBubbleContigIDOverlapped(numThread);
	pairedDBG.deleteSecondaryBubbleNodeAndEdge(numThread);
	pairedDBG.clearEdges();
	pairedDBG.makeScaffold();
	pairedDBG.joinUnambiguousNodePairIterative(numThread);


	pairedDBG.setMinOverlap(minOverlapForScaffolding);

	pairedDBG.setMode(PairedDBG::PAIRED_END_LINK_MODE | PairedDBG::LENGTH_CUTOFF_MODE | PairedDBG::PREVIOUS_DIVISION_AWARE_MODE);

	for (long iteration = 0; iteration < 2; ++iteration) {
		for (unsigned libraryIndex = 0; libraryIndex < libraryMT.size(); ++libraryIndex) {
			pairedDBG.setTargetLibraryIndex(libraryIndex);
			updateAndWriteInsertSize(libraryIndex);

			if (iteration == 0)
				pairedDBG.setMinLink(std::max(minLink, pairedDBG.estimateLink()));
			else
				pairedDBG.setMinLink(minLink);

			cerr << "[LIBRARY " << libraryIndex + 1 << "]\nAVE_INS = " << libraryMT[libraryIndex][0].getAverageInsSize() << ", SD_INS = " << libraryMT[libraryIndex][0].getSDInsSize() << endl;
			for (unsigned tolerenceFactor = MIN_TOL_FACTOR; tolerenceFactor <= MAX_TOL_FACTOR; ++tolerenceFactor) {
				pairedDBG.setCutoffLength(2.0 * tolerenceFactor * libraryMT[libraryIndex][0].getSDInsSize());
				pairedDBG.setTolerence(pairedDBG.getCutoffLength() / 2);
				cerr << "TOLERENCE_LEVEL_OF_CONTIG_OVERLAP = " << pairedDBG.getTolerence() << endl;

				pairedDBG.makeGraph(numThread);
				pairedDBG.setOppositeBubbleContigIDGapped(numThread);
				pairedDBG.deleteSecondaryBubbleNodeAndEdge(numThread);

				pairedDBG.deleteErroneousEdgeNumTagRateIterative(numThread);
				if (iteration > 0)
					pairedDBG.deleteErroneousEdgeIterative(numThread);

				pairedDBG.deleteRepeatEdge();
				pairedDBG.detectRepeat(2 * pairedDBG.getHeteroCoverage());
				pairedDBG.makeScaffold();
			}
		}

//		pairedDBG.divideErroneousNode(minLink, PairedDBG::DIVISION_MODE::MIS, numThread);
	}

	pairedDBG.setMinOverlap(this->contigMaxK - 1);
}


void SolveDBG::extendConsensus(const bool bubbleRemovalFlag, const bool tagFlag)
{
	pairedDBG.clearContigPreviousParentNodeID();

/*
pairedDBG.setMode(PairedDBG::LONG_READ_LINK_MODE | PairedDBG::LENGTH_CUTOFF_MODE);
//pairedDBG.setCutoffLength(4 * 0.1 * longReadLibraryMT[0].getAverageInsSize());
pairedDBG.setCutoffLength(10000);
pairedDBG.setMinLink(minLink);
//pairedDBG.setTolerence(0.1 * longReadLibraryMT[0].getAverageInsSize());
pairedDBG.setTolerence(10000);
pairedDBG.makeGraph(numThread);
pairedDBG.detectRepeat(pairedDBG.getAverageCoverage());
pairedDBG.dumpAllEdges("test_graph.tsv");
exit(0);
*/

	pairedDBG.setMode(PairedDBG::OVERLAP_MODE);
	pairedDBG.setMinLink(minLink);
	pairedDBG.makeGraph(numThread);
	pairedDBG.joinUnambiguousNodePairIterative(numThread);
	if (bubbleRemovalFlag) {
		pairedDBG.makeGraph(numThread);
		pairedDBG.setOppositeBubbleContigIDOverlapped(numThread);
		pairedDBG.deleteSecondaryBubbleNodeAndEdge(numThread);
		pairedDBG.clearEdges();
		pairedDBG.makeScaffold();
		pairedDBG.joinUnambiguousNodePairIterative(numThread);
	}

	for (unsigned i = 0; i < libraryMT.size(); ++i)
		updateAndWriteInsertSize(i);

	for (long outerIteration = 0; outerIteration < 2; ++outerIteration) {
		for (long iteration = 0; iteration < 2; ++iteration) {
			for (unsigned i = 0; i < libraryMT.size(); ++i) {
				cerr << "[LIBRARY " << i + 1 << "]\nAVE_INS = " << libraryMT[i][0].getAverageInsSize() << ", SD_INS = " << libraryMT[i][0].getSDInsSize() << endl;
				pairedDBG.setTargetLibraryIndex(i);
				unsigned tolerenceFactor = MAX_TOL_FACTOR;
				pairedDBG.setTolerence(tolerenceFactor * libraryMT[i][0].getSDInsSize());
				cerr << "TOLERENCE_OF_CONTIG_OVERLAP = " << pairedDBG.getTolerence() << endl;

				pairedDBG.setMinLink(minLinkToPhase);
				pairedDBG.joinUnambiguousNodePairIterative(numThread);
				pairedDBG.solveSimpleCrossStructureIterative(PairedDBG::CROSS_RESOLUTION_MODE::LINK, numThread);

				pairedDBG.setMinLink(minLink);
				pairedDBG.solveSimpleCrossStructureIterative(PairedDBG::CROSS_RESOLUTION_MODE::LINK, numThread);
			}

			if (longReadLibraryMT.size() > 0) {
				if (iteration == 0)
					pairedDBG.setMode(PairedDBG::OVERLAP_MODE | PairedDBG::LONG_READ_LINK_MODE);
				else
					pairedDBG.setMode(PairedDBG::OVERLAP_MODE | PairedDBG::LONG_READ_LINK_MODE | PairedDBG::NON_DBG_OVERLAP_MODE);
				cerr << "[LONG_READ_LIBRARY]\nAVE_READ_LENGTH = " << longReadLibraryMT[0].getAverageInsSize() << ", SD_READ_LENGTH = " << longReadLibraryMT[0].getSDInsSize() << endl;
				pairedDBG.setMinLink(minLinkToPhase);
				pairedDBG.joinUnambiguousNodePairIterative(numThread);
				pairedDBG.solveSimpleCrossStructureIterative(PairedDBG::CROSS_RESOLUTION_MODE::SCORE, numThread);

				pairedDBG.setMinLink(minLink);
				pairedDBG.solveSimpleCrossStructureIterative(PairedDBG::CROSS_RESOLUTION_MODE::TAG, numThread);
			}

			if (iteration == 0)
				pairedDBG.setMode(PairedDBG::OVERLAP_MODE | PairedDBG::PAIRED_END_LINK_MODE | PairedDBG::LONG_READ_LINK_MODE);
			else
				pairedDBG.setMode(PairedDBG::OVERLAP_MODE | PairedDBG::PAIRED_END_LINK_MODE | PairedDBG::LONG_READ_LINK_MODE | PairedDBG::NON_DBG_OVERLAP_MODE);
			for (unsigned i = libraryMT.size() - 1; i < libraryMT.size(); ++i) {
				pairedDBG.setTargetLibraryIndex(i);
				pairedDBG.setTolerenceFactor(MAX_TOL_FACTOR);
				pairedDBG.setMinLink(minLinkToPhase);
				pairedDBG.makeGraphAllLibraries(numThread);
				pairedDBG.joinUnambiguousNodePairIterative(numThread);
				pairedDBG.solveSimpleCrossStructureAllLibrariesIterative(PairedDBG::CROSS_RESOLUTION_MODE::LINK, numThread);

				pairedDBG.setMinLink(minLink);
				pairedDBG.solveSimpleCrossStructureAllLibrariesIterative(PairedDBG::CROSS_RESOLUTION_MODE::TAG, numThread);
			}
		}

		pairedDBG.trimRepeatEnd();

		if (optionBool["-no_scaffold"])
			continue;


		for (long iteration = 0; iteration < 2; ++iteration) {
			for (unsigned i = 0; i < libraryMT.size(); ++i) {
				pairedDBG.setTargetLibraryIndex(i);
				if (iteration == 0)
					pairedDBG.setMinLink(minLinkToPhase);
				else
					pairedDBG.setMinLink(minLink);

				cerr << "[LIBRARY " << i + 1 << "]\nAVE_INS = " << libraryMT[i][0].getAverageInsSize() << ", SD_INS = " << libraryMT[i][0].getSDInsSize() << endl;

				pairedDBG.setMode(PairedDBG::PAIRED_END_LINK_MODE | PairedDBG::LENGTH_CUTOFF_MODE);

				for (unsigned tolerenceFactor = MIN_TOL_FACTOR; tolerenceFactor <= MAX_TOL_FACTOR; ++tolerenceFactor) {
					pairedDBG.setCutoffLength(tolerenceFactor * std::min(1.0 * libraryMT[i][0].getSDInsSize(), 0.1 * libraryMT[i][0].getAverageInsSize()));
					pairedDBG.setTolerence(pairedDBG.getCutoffLength());
					cerr << "TOLERENCE_LEVEL_OF_CONTIG_OVERLAP = " << pairedDBG.getTolerence() << endl;
					cerr << "LENGTH_CUTOFF = " << pairedDBG.getCutoffLength() << endl;

					pairedDBG.trimSparseEnd();
					pairedDBG.trimRepeatEnd();
					pairedDBG.solveSimpleGappedCrossStructureIterative(PairedDBG::CROSS_RESOLUTION_MODE::LINK, numThread);
					pairedDBG.solveSimpleGappedCrossStructureIterative(PairedDBG::CROSS_RESOLUTION_MODE::TAG, numThread);
				}
			}
			if (longReadLibraryMT.size() > 0) {
				pairedDBG.setMode(PairedDBG::LONG_READ_LINK_MODE);
				cerr << "[LONG_READ_LIBRARY]\nAVE_READ_LENGTH = " << longReadLibraryMT[0].getAverageLength() << endl;
				pairedDBG.setTolerence(2 * this->contigMaxK);
				if (iteration == 0)
					pairedDBG.setMinLink(minLinkToPhase);
				else
					pairedDBG.setMinLink(minLink);

				for (unsigned cutoffFactor = MIN_LONG_READ_LENGTH_CUTOFF_FACTOR; cutoffFactor <= MAX_LONG_READ_LENGTH_CUTOFF_FACTOR; ++cutoffFactor) {
					pairedDBG.setCutoffLength(cutoffFactor * LONG_READ_LENGTH_CUTOFF_UNIT_FACTOR * longReadLibraryMT[0].getAverageInsSize());
					pairedDBG.trimSparseEnd();
					pairedDBG.solveSimpleGappedCrossStructureIterative(PairedDBG::CROSS_RESOLUTION_MODE::SCORE, numThread);
					pairedDBG.solveSimpleGappedCrossStructureIterative(PairedDBG::CROSS_RESOLUTION_MODE::TAG, numThread);
				}
			}
			for (unsigned i = libraryMT.size() - 1; i < libraryMT.size(); ++i) {
				pairedDBG.setTargetLibraryIndex(i);
				pairedDBG.setMode(PairedDBG::PAIRED_END_LINK_MODE | PairedDBG::LONG_READ_LINK_MODE | PairedDBG::LENGTH_CUTOFF_MODE);

				for (unsigned tolerenceFactor = MIN_TOL_FACTOR; tolerenceFactor <= MAX_TOL_FACTOR; ++tolerenceFactor) {
					pairedDBG.setTolerenceFactor(tolerenceFactor);
					pairedDBG.setCutoffLengthFactor(MIN_LONG_READ_LENGTH_CUTOFF_FACTOR * LONG_READ_LENGTH_CUTOFF_UNIT_FACTOR);
					pairedDBG.setCutoffLength(2.0 * tolerenceFactor * std::min(1.0 * libraryMT[i][0].getSDInsSize(), 0.1 * libraryMT[i][0].getAverageInsSize()));
					pairedDBG.setTolerence(pairedDBG.getCutoffLength() / 2);
					if (iteration == 0)
						pairedDBG.setMinLink(minLinkToPhase);
					else
						pairedDBG.setMinLink(minLink);
					pairedDBG.trimSparseEnd();
					pairedDBG.trimRepeatEnd();
					pairedDBG.solveSimpleGappedCrossStructureAllLibrariesIterative(PairedDBG::CROSS_RESOLUTION_MODE::LINK, numThread);
					pairedDBG.solveSimpleGappedCrossStructureAllLibrariesIterative(PairedDBG::CROSS_RESOLUTION_MODE::TAG, numThread);
				}
			}
		}

		pairedDBG.trimSparseEnd();
		pairedDBG.trimRepeatEnd();
		pairedDBG.divideErroneousNodeBaseLevel(0, libraryMT.size(), this->contigReadLength, false, true, false, numThread);


		pairedDBG.setMinOverlap(minOverlapForScaffolding);

		for (long iteration = 0; iteration < 2; ++iteration) {
			pairedDBG.setMode(PairedDBG::PAIRED_END_LINK_MODE | PairedDBG::LENGTH_CUTOFF_MODE | PairedDBG::PREVIOUS_DIVISION_AWARE_MODE);

			for (unsigned libraryIndex = 0; libraryIndex < libraryMT.size(); ++libraryIndex) {
				updateAndWriteInsertSize(libraryIndex);

				pairedDBG.setTargetLibraryIndex(libraryIndex);
				pairedDBG.setMinLink(1);

				cerr << "[LIBRARY " << libraryIndex + 1 << "]\nAVE_INS = " << libraryMT[libraryIndex][0].getAverageInsSize() << ", SD_INS = " << libraryMT[libraryIndex][0].getSDInsSize() << endl;
				for (unsigned tolerenceFactor = MIN_TOL_FACTOR; tolerenceFactor <= MAX_TOL_FACTOR; ++tolerenceFactor) {
					pairedDBG.setCutoffLength(tolerenceFactor * libraryMT[libraryIndex][0].getSDInsSize());
					pairedDBG.setTolerence(pairedDBG.getCutoffLength());

					cerr << "TOLERENCE_LEVEL_OF_CONTIG_OVERLAP = " << pairedDBG.getTolerence() << endl;

					pairedDBG.makeGraph(numThread);
					if (iteration == 0)
						pairedDBG.deleteThinEdge(std::max(minLink, pairedDBG.estimateLink()));
					else
						pairedDBG.deleteThinEdge(minLink);

					if (bubbleRemovalFlag) {
						pairedDBG.setOppositeBubbleContigIDGapped(numThread);
						pairedDBG.deleteSecondaryBubbleNodeAndEdge(numThread);
					}

					if (tagFlag)
						pairedDBG.deleteErroneousEdgeNumTagRateIterative(numThread);
					if (iteration > 0)
						pairedDBG.deleteRepeatEdge();

					pairedDBG.detectRepeat(pairedDBG.getAverageCoverage());
					pairedDBG.makeScaffold();
				}

			}

			pairedDBG.divideErroneousNodeBaseLevel(0, libraryMT.size(), this->contigReadLength, false, true, false, numThread);
		}

		pairedDBG.trimSparseEnd();
		pairedDBG.trimRepeatEnd();



//		if (libraryMT.size() > 0) {
//			for (long iteration = 0; iteration < 2; ++iteration) {
//				pairedDBG.setTargetLibraryIndex(libraryMT.size() - 1);
//				long linkThreashold = (1 + iteration) * minLink;
//				pairedDBG.setMode(PairedDBG::LONG_READ_LINK_MODE | PairedDBG::PAIRED_END_LINK_MODE | PairedDBG::LENGTH_CUTOFF_MODE);
//
//				pairedDBG.setTolerenceFactor(MAX_TOL_FACTOR);
//				pairedDBG.setCutoffLengthFactor(MAX_LONG_READ_LENGTH_CUTOFF_FACTOR * LONG_READ_LENGTH_CUTOFF_UNIT_FACTOR);
//				pairedDBG.setMinLink(linkThreashold);
//				pairedDBG.setCutoffLength(MAX_TOL_FACTOR * libraryMT.back()[0].getSDInsSize());
//				pairedDBG.setTolerence(pairedDBG.getCutoffLength());
//				pairedDBG.makeGraphAllLibraries(numThread);
//				if (bubbleRemovalFlag) {
//					pairedDBG.setOppositeBubbleContigIDGapped(numThread);
//					pairedDBG.deleteSecondaryBubbleNodeAndEdge(numThread);
//				}
//				if (iteration > 0)
//					pairedDBG.deleteErroneousEdgeNumLinkRateIterative(numThread);
//				pairedDBG.deleteRepeatEdge();
//				pairedDBG.detectRepeat(pairedDBG.getAverageCoverage());
//				pairedDBG.makeScaffold();
//			}
//
//			pairedDBG.divideErroneousNodeBaseLevel(0, libraryMT.size(), this->contigReadLength, false, true, false, numThread);
//		}



		if (longReadLibraryMT.size() > 0) {
			pairedDBG.setMode(PairedDBG::LONG_READ_LINK_MODE | PairedDBG::LENGTH_CUTOFF_MODE);

			for (unsigned cutoffFactor = MIN_LONG_READ_LENGTH_CUTOFF_FACTOR; cutoffFactor <= MAX_LONG_READ_LENGTH_CUTOFF_FACTOR; ++cutoffFactor) {
				pairedDBG.setCutoffLength(cutoffFactor * LONG_READ_LENGTH_CUTOFF_UNIT_FACTOR *  longReadLibraryMT[0].getAverageInsSize());
				pairedDBG.setTolerence(std::min(MIN_LONG_READ_LENGTH_CUTOFF_FACTOR * longReadLibraryMT[0].getAverageInsSize(), 0.5 * pairedDBG.getCutoffLength()));
				pairedDBG.setMinLink(minLink);
				pairedDBG.makeGraph(numThread);
				if (bubbleRemovalFlag) {
					pairedDBG.setOppositeBubbleContigIDGapped(numThread);
					pairedDBG.deleteSecondaryBubbleNodeAndEdge(numThread);
				}
				pairedDBG.deleteErroneousEdgeScore(0.125, numThread);
				pairedDBG.deleteRepeatEdge();
				pairedDBG.detectRepeat(pairedDBG.getAverageCoverage());
				pairedDBG.makeScaffoldCombine();
			}

			pairedDBG.divideErroneousNodeBaseLevel(0, libraryMT.size(), this->contigReadLength, false, true, false, numThread);
		}

		pairedDBG.setMinOverlap(this->contigMaxK - 1);

		pairedDBG.trimSparseEnd();
		pairedDBG.trimRepeatEnd();
	}
}


void SolveDBG::combineAssembly()
{
	pairedDBG.clearContigPreviousParentNodeID();
	pairedDBG.setMinLink(1);
	pairedDBG.setMode(PairedDBG::LONG_READ_LINK_MODE | PairedDBG::LENGTH_CUTOFF_MODE);

/*
pairedDBG.setCutoffLength(atol(optionSingleArgs["-combine_l"].c_str()));
pairedDBG.setTolerence(atol(optionSingleArgs["-combine_t"].c_str()));
pairedDBG.makeGraph(numThread);
pairedDBG.detectRepeat(pairedDBG.getAverageCoverage());
pairedDBG.dumpAllEdges("test_graph.tsv");
*/

	long minLengthCutoff = atol(optionSingleArgs["-combine_l"].c_str());
	long maxLengthCutoff = atol(optionSingleArgs["-combine_L"].c_str());
	long numStep = atol(optionSingleArgs["-combine_s"].c_str());
	long lengthCutoffStep = (maxLengthCutoff - minLengthCutoff) / (numStep - 1);

	for (long lengthCutoff = minLengthCutoff; lengthCutoff <= maxLengthCutoff; lengthCutoff += lengthCutoffStep) {
		pairedDBG.setCutoffLength(lengthCutoff);
		pairedDBG.setTolerence(atol(optionSingleArgs["-combine_t"].c_str()));
		pairedDBG.makeGraph(numThread);
		pairedDBG.deleteErroneousEdgeScore(1.0, numThread);
		pairedDBG.makeScaffoldCombine();
		pairedDBG.divideGappedNode(atoi(optionSingleArgs["-combine_g"].c_str()));
	}
}
