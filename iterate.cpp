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

#include "iterate.h"
#include <sys/stat.h>
#include <sys/types.h>


const std::string IterateScaffold::CONTIG_FOOTER = "_contig.fa";
const std::string IterateScaffold::KMER_DIVIDE_FOOTER = "_kmerDivided.fa";
const std::string IterateScaffold::DIVIDE_FOOTER = "_divided.fa";
const std::string IterateScaffold::SCAF_FOOTER = "_consensusScaffold.fa";
const std::string IterateScaffold::POLISH_FOOTER = "_polished_consensusScaffold.fa";
const std::string IterateScaffold::GAP_FOOTER = "_gapClosed_polished_consensusScaffold.fa";
const std::string IterateScaffold::EX_FOOTER = "_extraContig.fa";
const std::string IterateScaffold::MERGE_FOOTER = "_merged.fa";
const std::string IterateScaffold::ITERATION_FOOTER = "_iterativeAssembly.fa";


//////////////////////////////////////////////////////////////////////////////////////
// constructor
//////////////////////////////////////////////////////////////////////////////////////
IterateScaffold::IterateScaffold()
{
    optionSingleArgs["-o"] = "out";
    optionSingleArgs["-i"] = "6";
    optionSingleArgs["-m"] = "1";
    optionSingleArgs["-l"] = "";
    optionSingleArgs["-t"] = "1";
    optionSingleArgs["-u"] = "0";
    optionSingleArgs["-sub_bin"] = platanus::ConstParam::SUB_BIN_PATH;
    optionMultiArgs["-c"] = std::vector<std::string>();
    optionMultiArgs["-p"] = std::vector<std::string>();
    optionMultiArgs["-ont"] = std::vector<std::string>();
    optionMultiArgs["-gc"] = std::vector<std::string>();
    pairedEndSingleFileType.emplace_back("-ip");
    pairedEndSingleFileType.emplace_back("-op");
    pairedEndPairFileType.emplace_back("-IP");
    pairedEndPairFileType.emplace_back("-OP");
    optionSingleArgs["-tmp"] = ".";
    optionBool["-trim_overlap"] = false;
    optionBool["-keep_file"] = false;
}




//////////////////////////////////////////////////////////////////////////////////////
// usage
//////////////////////////////////////////////////////////////////////////////////////
void IterateScaffold::usage(void) const
{

    std::cerr << "\nUsage: platanus_b iterate [Options]\n"
              << "Options:\n"
              << "    -o STR                             : prefix of output file and directory (do not use \"/\", default " << optionSingleArgs.at("-o") << ", length <= " << platanus::ConstParam::MAX_FILE_LEN << ")\n"
              << "    -c FILE1 [FILE2 ...]               : contig (or scaffold) file (fasta format)\n"
              << "    -i INT                             : number of iterations (default " << optionSingleArgs.at("-i") << ")\n"
              << "    -l INT                             : -l value of \"scaffold\" step\n"
              << "    -u FLOAT                           : maximum difference for bubble crush (identity, default " << optionSingleArgs.at("-u") << ")\n"
              << "    -ip{INT} PAIR1 [PAIR2 ...]         : lib_id inward_pair_file (reads in 1 file, fasta or fastq)\n"
              << "    -IP{INT} FWD1 REV1 [FWD2 REV2 ...] : lib_id inward_pair_files (reads in 2 files, fasta or fastq)\n"
              << "    -op{INT} PAIR1 [PAIR2 ...]         : lib_id outward_pair_file (reads in 1 file, fasta or fastq)\n"
              << "    -OP{INT} FWD1 REV1 [FWD2 REV2 ...] : lib_id outward_pair_files (reads in 2 files, fasta or fastq)\n"
              << "    -ont FILE1 [FILE2 ...]             : Oxford Nanopore long-read file (fasta or fastq)\n"
              << "    -p FILE1 [FILE2 ...]               : PacBio long-read file (fasta or fastq)\n"
              << "    -gc FILE1 [FILE2 ...]              : Guiding contig file; i.e. other assemblies, synthetic long-reads or corrected reads (fasta or fastq)\n"
              << "    -t INT                             : number of threads (default " << optionSingleArgs.at("-t") << ")\n"
              << "    -m INT                             : memory limit for making kmer distribution (GB, >=1, default " << optionSingleArgs.at("-m") << ")\n"
              << "    -tmp DIR                           : directory for temporary files (default " << optionSingleArgs.at("-tmp") << ")\n"
              << "    -sub_bin DIR                       : directory for binary files which platanus_b use internally (e.g. minimap2) (default " <<  optionSingleArgs.at("-sub_bin") << ")\n"
              << "    -keep_file                         : keep intermediate files (default, off)\n"
              << "    -trim_overlap                      : trim overlapping edges of scaffolds (default, off)\n"
              << "\n"

              << "Input format:\n"
              << "    Uncompressed and compressed (gzip or bzip2) files are accepted for -c, -ip, -IP, -op, -OP, -ont, -p and -gc options.\n"
              << "\n\n"

              << "Outputs:\n"
              << "    PREFIX_iterativeAssembly.fa\n"
              << std::endl;
}



//////////////////////////////////////////////////////////////////////////////////////
// main
//////////////////////////////////////////////////////////////////////////////////////
void IterateScaffold::exec(void)
{
	platanus::setGlobalTmpFileDir(optionSingleArgs["-tmp"].c_str());
	setSubBinPath(optionSingleArgs["-sub_bin"]);

    setIntermediateDirectoryName(optionSingleArgs["-o"] + "_iterateIntermediateResults");
    struct stat st;
	if (stat(intermediateDirectoryName.c_str(), &st) != 0) {
		int returnValue = mkdir(intermediateDirectoryName.c_str(), 0755);
		if (returnValue != 0) {
			throw platanus::CreateDirError(intermediateDirectoryName);
		}
	}

    const long iterateTimes = atoi(optionSingleArgs["-i"].c_str());

	countKmer();
    for (long times = 1; times <= iterateTimes; ++times) {
        setDirectoryName(intermediateDirectoryName, times);
        createDirectory();
        createContig(times);
        execKmerDivideMode();
        execScaffoldMode(times, iterateTimes);
        execPolishMode();
		execGapCloseMode(times, iterateTimes);
    }
	execFinalDivideMode();
	execFinalPolishMode();

    for (long times = iterateTimes/2; times <= iterateTimes - 1; ++times) {
		execCombineMode(times);
		execCombinatorialGapClosePerl(times);
		execRemoveRedundantSeqPerl(times);
	}

    std::ostringstream oss;
    oss << "mv " <<  directoryName << "/" << optionSingleArgs["-o"] << iterateTimes - 1 << "_closed.fa.rmred"
        << " " << optionSingleArgs["-o"] << ITERATION_FOOTER
		<< "; rm " << intermediateDirectoryName << "/" << optionSingleArgs["-o"] << "_kmer_occ.bin"
	;
    if (!optionBool["-keep_file"]) {
		oss << "; rm " << intermediateDirectoryName << "/*/*.fa"
			<< "; rm " << intermediateDirectoryName << "/*/*.bed"
			<< "; rm " << intermediateDirectoryName << "/*/*.tsv"
			<< "; rm " << intermediateDirectoryName << "/*.tsv"
			<< "; rm -rf " << intermediateDirectoryName << "/*/*_intermediates"
		;
	}

    if (system(oss.str().c_str()) != 0) {
        throw platanus::CreateLinkError();
    }
}


void IterateScaffold::setSubBinPath(std::string path)
{
	this->subBinPath = path;

	const std::vector<std::string> subBinList = {
		"close_gap.pl", \
		"combinatorial_gap_close.pl", \
		"fasta_around_gap.pl", \
		"fasta_grep.pl", \
		"get_flanked_region_info_outer.pl", \
		"get_flanked_region_info.pl", \
		"paf_contained_short_seq_list_bubble_aware.pl", \
		"paf_filter_flanking_pair.pl", \
		"paf_filter_qcov.pl", \
		"paf_match_short_seq_list_bubble_aware.pl", \
		"paf_max_match_unique.pl", \
		"reduce_filled_info.pl", \
		"remove_redundant_seq.pl", \
		"minimap2"
	};


	std::cerr << "checking sub-executables ..." << std::endl;
	for (auto it = subBinList.begin(); it != subBinList.end(); ++it) {
		std::ostringstream oss;
		oss << "which " << subBinPath << "/" << *it << ">&2";

		if (system(oss.str().c_str()) != 0) {
			throw platanus::SubBinError(*it);
		}
	}
	std::cerr << "all sub-executables found\n" << std::endl;
}


void IterateScaffold::createDirectory(void) const
{
    struct stat st;

    if (stat(directoryName.c_str(), &st) != 0) {
        int returnValue = mkdir(directoryName.c_str(), 0755);
        if (returnValue != 0) {
            throw platanus::CreateDirError(directoryName);
        }
    }
}


void IterateScaffold::createContig(const long times)
{
    if (times == 1) {
        std::ostringstream oss;
        oss << "cat ";
		for (auto argItr = optionMultiArgs["-c"].begin(); argItr != optionMultiArgs["-c"].end(); ++argItr) {
			oss << " " << *argItr;
		}
		oss << " > " << this->directoryName << "/" << optionSingleArgs["-o"] << CONTIG_FOOTER;

        if (system(oss.str().c_str()) != 0) {
            throw platanus::CreateLinkError();
        }
    } else {
        this->execMergeMode(times);
        std::ostringstream oss;
        oss << "cat " << this->directoryName << "/" << optionSingleArgs["-o"] << "_merged.fa"
            << " " << this->directoryName << "/" << optionSingleArgs["-o"] << "_mergedJunctionKmer.fa"
            << " >" << this->directoryName << "/" << optionSingleArgs["-o"] << CONTIG_FOOTER
            << "; rm " <<  this->directoryName << "/" << optionSingleArgs["-o"] << "_merged.fa"
            << " " << this->directoryName << "/" << optionSingleArgs["-o"] << "_mergedJunctionKmer.fa";

        if (system(oss.str().c_str()) != 0) {
            throw platanus::CreateLinkError();
        }
    }
}


void IterateScaffold::countKmer(void)
{
	long k = platanus::Contig().getMaxKFromFastaHeader(optionMultiArgs["-c"][0]);

	std::ostringstream oss;
    oss << this->platanusBinary << " assemble -kmer_occ_only -n 1"
        << " -t " << optionSingleArgs["-t"]
        << " -m " << optionSingleArgs["-m"]
        << " -tmp " << optionSingleArgs["-tmp"]
        << " -k " << k
        << " -o " << this->intermediateDirectoryName << "/" << optionSingleArgs["-o"]
        << " -f"
	;

    for (auto itr = optionPairFile.begin(), end = optionPairFile.end(); itr != end; ++itr) {
        oss << " " << itr->fileFirst;
        if (itr->fileSecond != "") {
            oss << " " << itr->fileSecond;
        }
    }

    oss  << " 2>" << this->intermediateDirectoryName << "/" << optionSingleArgs["-o"] << ".countKmerLog";

	if (system(oss.str().c_str()) != 0) {
		throw platanus::CreateLinkError();
	}
}


void IterateScaffold::execKmerDivideMode(void)
{
    std::ostringstream oss;
    oss << this->platanusBinary << " kmer_divide"
        << " -t " << optionSingleArgs["-t"]
        << " -tmp " << optionSingleArgs["-tmp"]
        << " -k " << this->intermediateDirectoryName << "/" << optionSingleArgs["-o"] << "_kmer_occ.bin"
        << " -f " << this->directoryName << "/" << optionSingleArgs["-o"] << CONTIG_FOOTER
        << " -o " << this->directoryName << "/" << optionSingleArgs["-o"]
		<< " 2>" << this->directoryName << "/" << optionSingleArgs["-o"] << ".kmerDivideLog"
	;

    if (system(oss.str().c_str()) != 0) {
        throw platanus::KmerDivideError();
    }
}


void IterateScaffold::execMergeMode(const long times)
{
    std::ostringstream oss;
    oss << this->platanusBinary << " merge"
        << " -m " << optionSingleArgs["-m"]
        << " -tmp " << optionSingleArgs["-tmp"]
        << " -f " << this->previousDirectoryName << "/" << optionSingleArgs["-o"] << GAP_FOOTER
        << " " << this->previousDirectoryName << "/" << optionSingleArgs["-o"] << EX_FOOTER
        << " -k " << 1.0 + 0.5 * ((times - 1) / 3)
        << " -l " << 1.0 + 0.5 * ((times - 1) / 3)
        << " -u " <<  optionSingleArgs["-u"]
        << " -o " << this->directoryName << "/" << optionSingleArgs["-o"]
        << " 2>" << this->directoryName << "/" << optionSingleArgs["-o"] << ".mergeLog";

    if (system(oss.str().c_str()) != 0) {
        throw platanus::MergeError();
    }
}


void IterateScaffold::execDivideMode(void)
{
    std::ostringstream oss;
    oss << this->platanusBinary << " divide" 
        << " -t " << optionSingleArgs["-t"]
        << " -tmp " << optionSingleArgs["-tmp"]
        << " -c " << this->directoryName << "/" << optionSingleArgs["-o"] << GAP_FOOTER
        << " -o " << this->directoryName << "/" << optionSingleArgs["-o"];

    for (auto itr = optionPairFile.begin(), end = optionPairFile.end(); itr != end; ++itr) {
        oss << " " << itr->libraryType << itr->libraryID << " " << itr->fileFirst;
        if (itr->fileSecond != "") {
            oss << " " << itr->fileSecond;
        }
    }

    oss << " 2>" << this->directoryName << "/" << optionSingleArgs["-o"] << ".divLog";

    if (system(oss.str().c_str()) != 0) {
        throw platanus::PolishError();
    }
}


void IterateScaffold::execScaffoldMode(const long times, const long numTimes)
{
    std::ostringstream oss;
    oss << this->platanusBinary << " solve_DBG -unphase"
        << " -t " << optionSingleArgs["-t"]
        << " -tmp " << optionSingleArgs["-tmp"]
        << " -c " << this->directoryName << "/" << optionSingleArgs["-o"] << KMER_DIVIDE_FOOTER
        << " -u " <<  optionSingleArgs["-u"]
        << " -o " << this->directoryName << "/" << optionSingleArgs["-o"]
		<< " -reduce_redundancy"
	;

    if (optionSingleArgs["-l"] != "") {
        oss << " -l " << optionSingleArgs["-l"];
    }
    if (optionSingleArgs["-r"] != "") {
        oss << " -r " << optionSingleArgs["-r"];
    }
	if (times == numTimes && optionBool["-trim_overlap"]) {
		 oss << " -trim_overlap";
	}

    for (auto itr = optionPairFile.begin(), end = optionPairFile.end(); itr != end; ++itr) {
        oss << " " << itr->libraryType << itr->libraryID << " " << itr->fileFirst;
        if (itr->fileSecond != "") {
            oss << " " << itr->fileSecond;
        }
    }

	const std::vector<std::string> longReadOption = {"-p", "-ont", "-gc"};
	for (auto it = longReadOption.begin(); it != longReadOption.end(); ++it) {
		if (!(optionMultiArgs[*it].empty()) && times >= numTimes/2) {
			oss << " " << *it;
			for (auto argItr = optionMultiArgs[*it].begin(); argItr != optionMultiArgs[*it].end(); ++argItr) {
				oss << " " << *argItr;
			}
		}
	}

    oss << " 2>" << this->directoryName << "/" << optionSingleArgs["-o"] << ".scafLog";

    if (system(oss.str().c_str()) != 0) {
        throw platanus::ScaffoldError();
    }
}


void IterateScaffold::execPolishMode(void)
{
    std::ostringstream oss;
    oss << this->platanusBinary << " polish"
        << " -t " << optionSingleArgs["-t"]
        << " -tmp " << optionSingleArgs["-tmp"]
        << " -c " << this->directoryName << "/" << optionSingleArgs["-o"] << SCAF_FOOTER
        << " -o " << this->directoryName << "/" << optionSingleArgs["-o"];

    for (auto itr = optionPairFile.begin(), end = optionPairFile.end(); itr != end; ++itr) {
        oss << " " << itr->libraryType << itr->libraryID << " " << itr->fileFirst;
        if (itr->fileSecond != "") {
            oss << " " << itr->fileSecond;
        }
    }

    oss << " 2>" << this->directoryName << "/" << optionSingleArgs["-o"] << ".polLog";

    if (system(oss.str().c_str()) != 0) {
        throw platanus::PolishError();
    }
}


void IterateScaffold::execFinalDivideMode(void)
{
    std::ostringstream oss;
    oss << this->platanusBinary << " divide"
        << " -t " << optionSingleArgs["-t"]
        << " -tmp " << optionSingleArgs["-tmp"]
        << " -c " << this->directoryName << "/" << optionSingleArgs["-o"] << GAP_FOOTER
        << " -o " << this->directoryName << "/" << optionSingleArgs["-o"] << "_final"
	;

    for (auto itr = optionPairFile.begin(), end = optionPairFile.end(); itr != end; ++itr) {
        oss << " " << itr->libraryType << itr->libraryID << " " << itr->fileFirst;
        if (itr->fileSecond != "") {
            oss << " " << itr->fileSecond;
        }
    }

	const std::vector<std::string> longReadOption = {"-p", "-ont", "-gc"};
	for (auto it = longReadOption.begin(); it != longReadOption.end(); ++it) {
		if (!(optionMultiArgs[*it].empty())) {
			oss << " " << *it;
			for (auto argItr = optionMultiArgs[*it].begin(); argItr != optionMultiArgs[*it].end(); ++argItr) {
				oss << " " << *argItr;
			}
		}
	}

    oss << " 2>" << this->intermediateDirectoryName << "/" << optionSingleArgs["-o"] << ".finalDivideLog";

    if (system(oss.str().c_str()) != 0) {
        throw platanus::DivideError();
    }
}


void IterateScaffold::execFinalPolishMode(void)
{
    std::ostringstream oss;
    oss << this->platanusBinary << " polish"
        << " -t " << optionSingleArgs["-t"]
        << " -tmp " << optionSingleArgs["-tmp"]
        << " -c " << this->directoryName << "/" << optionSingleArgs["-o"] << "_final" << DIVIDE_FOOTER
        << " -o " << this->directoryName << "/" << optionSingleArgs["-o"]
	;

    for (auto itr = optionPairFile.begin(), end = optionPairFile.end(); itr != end; ++itr) {
        oss << " " << itr->libraryType << itr->libraryID << " " << itr->fileFirst;
        if (itr->fileSecond != "") {
            oss << " " << itr->fileSecond;
        }
    }

    oss << " 2>" << this->intermediateDirectoryName << "/" << optionSingleArgs["-o"] << ".finalPolishLog";

    if (system(oss.str().c_str()) != 0) {
        throw platanus::PolishError();
    }
}


void IterateScaffold::execGapCloseMode(const long times, const long numTimes)
{
    std::ostringstream oss;
    oss << this->platanusBinary << " gap_close"
        << " -t " << optionSingleArgs["-t"]
        << " -tmp " << optionSingleArgs["-tmp"]
        << " -c " << this->directoryName << "/" << optionSingleArgs["-o"] << POLISH_FOOTER
        << " -o " << this->directoryName << "/" << optionSingleArgs["-o"]
		<< " -reduce_redundancy"
	;

	if (times < numTimes) {
		oss << " -extend";
	}

    for (auto itr = optionPairFile.begin(), end = optionPairFile.end(); itr != end; ++itr) {
        oss << " " << itr->libraryType << itr->libraryID << " " << itr->fileFirst;
        if (itr->fileSecond != "") {
            oss << " " << itr->fileSecond;
        }
    }
    oss << " 2>" << this->directoryName << "/" << optionSingleArgs["-o"] << ".gapLog";

    if (system(oss.str().c_str()) != 0) {
        throw platanus::GapError();
    }
}


void IterateScaffold::execCombineMode(const long times)
{
	std::ostringstream oss;
	oss << this->platanusBinary  << " solve_DBG -combine";

	std::ostringstream fileOss;
	fileOss << directoryName << "/" << optionSingleArgs["-o"] << times << "_combined.fa";
    struct stat st;
	if (stat(fileOss.str().c_str(), &st) != 0)
		oss <<  " -c " << directoryName << "/" << optionSingleArgs["-o"] << "_polished_final" << DIVIDE_FOOTER;
	else
		oss <<  " -c " << directoryName << "/" << optionSingleArgs["-o"] << times << "_combined.fa";

	oss << " -t " << optionSingleArgs["-t"] \
		<< " -tmp " << optionSingleArgs["-tmp"] \
		<< " -gc " << intermediateDirectoryName << "/" << optionSingleArgs["-o"] << times << "/" << optionSingleArgs["-o"] << "_gapClosed_polished_consensusScaffold.fa" \
		<< " -o " <<  directoryName << "/" << optionSingleArgs["-o"] << times \
		<< " > " <<  directoryName << "/" << optionSingleArgs["-o"] << times << ".combineLog 2>&1" \
	;

	std::cerr << oss.str() << std::endl;
    if (system(oss.str().c_str()) != 0) {
        throw platanus::ScaffoldError();
	}
}


void IterateScaffold::execCombinatorialGapClosePerl(const long times)
{
	std::ostringstream oss;
	oss << this->subBinPath  << "/combinatorial_gap_close.pl"
		<< " -t " << optionSingleArgs["-t"]
		<< " -b " << directoryName << "/" << optionSingleArgs["-o"] << times << "_combined.fa"
		<< " -c " << intermediateDirectoryName << optionSingleArgs["-o"] << times << "/" << optionSingleArgs["-o"] << "_gapClosed_polished_consensusScaffold.fa"
		<< " -o " << directoryName << "/" << optionSingleArgs["-o"] << times
		<< " > "  << directoryName << "/" << optionSingleArgs["-o"] << times << ".combinatorialGapCloseLog 2>&1"
	;

	std::cerr << oss.str() << std::endl;
    if (system(oss.str().c_str()) != 0) {
        throw platanus::SubBinError("combinatorial_gap_close.pl");
	}
}


void IterateScaffold::execRemoveRedundantSeqPerl(const long times)
{
	std::ostringstream oss;
	oss << this->subBinPath  << "/remove_redundant_seq.pl"
		<< " -t " << optionSingleArgs["-t"]
		<< " -o " << directoryName << "/" << optionSingleArgs["-o"] << times
		<< " " <<  directoryName << "/" << optionSingleArgs["-o"] << times << "_closed.fa"
		<< " > "  << directoryName << "/" << optionSingleArgs["-o"] << times << ".removeRedundantSeqLog 2>&1"
	;

	std::cerr << oss.str() << std::endl;
    if (system(oss.str().c_str()) != 0) {
        throw platanus::SubBinError("remove_redundant_seq.pl");
	}
}
