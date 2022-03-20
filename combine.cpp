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

#include "combine.h"
#include <sys/stat.h>
#include <sys/types.h>


//////////////////////////////////////////////////////////////////////////////////////
// constructor
//////////////////////////////////////////////////////////////////////////////////////
Combine::Combine()
{
    optionSingleArgs["-o"] = "out";
    optionSingleArgs["-t"] = "1";
    optionSingleArgs["-sub_bin"] = platanus::ConstParam::SUB_BIN_PATH;
    optionMultiArgs["-c"] = std::vector<std::string>();
    optionMultiArgs["-gc"] = std::vector<std::string>();
    optionSingleArgs["-tmp"] = ".";
    optionBool["-no_gap_close"] = false;
    optionBool["-keep_file"] = false;

	optionSingleArgs["-combine_l"] = "10000";
	optionSingleArgs["-combine_L"] = "100000";
	optionSingleArgs["-combine_t"] = "10000";
	optionSingleArgs["-combine_s"] = "10";
	optionSingleArgs["-combine_g"] = "100000";
	optionSingleArgs["-combine_i"] = "0.9";
}




//////////////////////////////////////////////////////////////////////////////////////
// usage
//////////////////////////////////////////////////////////////////////////////////////
void Combine::usage(void) const
{

    std::cerr << "\nUsage: platanus_b combine [Options]\n"
              << "Options:\n"
              << "    -o STR                             : prefix of output file and directory (do not use \"/\", default " << optionSingleArgs.at("-o") << ", length <= " << platanus::ConstParam::MAX_FILE_LEN << ")\n"
              << "    -c FILE1 [FILE2 ...]               : contig (or scaffold) file (fasta format)\n"
              << "    -gc FILE1 [FILE2 ...]              : Guiding contig file; i.e. other assemblies, synthetic long-reads or corrected reads (fasta or fastq)\n"
              << "    -p FILE1 [FILE2 ...]               : PacBio long-read file (fasta or fastq)\n"
              << "    -ont FILE1 [FILE2 ...]             : Oxford Nanopore long-read file (fasta or fastq)\n"
              << "    -t INT                             : number of threads (default " << optionSingleArgs.at("-t") << ")\n"
              << "    -tmp DIR                           : directory for temporary files (default " << optionSingleArgs.at("-tmp") << ")\n"
              << "    -sub_bin DIR                       : directory for binary files which platanus_b use internally (e.g. minimap2) (default " <<  optionSingleArgs.at("-sub_bin") << ")\n"
              << "    -no_gap_close                      : not close gaps by guiding contigs (default, false)\n"
              << "    -keep_file                         : keep intermediate files (default, off)\n"
              << "    -combine_g INT                     : maximug gap-size in scaffolding (default " << optionSingleArgs.at("-combine_g") << ")\n"
              << "    -combine_i INT                     : minimum identity in scaffolding (0-1, default " << optionSingleArgs.at("-combine_i") << ")\n"
              << "    -combine_l INT                     : minimum length-cutoff in scaffolding (default " << optionSingleArgs.at("-combine_l") << ")\n"
              << "    -combine_L INT                     : maximum length-cutoff in scaffolding (default " << optionSingleArgs.at("-combine_L") << ")\n"
              << "    -combine_t INT                     : length tolerance to detect conflicts (default " << optionSingleArgs.at("-combine_t") << ")\n"
              << "    -combine_s INT                     : number of steps in scaffolding (default " << optionSingleArgs.at("-combine_s") << ")\n"
              << "\n"

              << "Input format:\n"
              << "    Uncompressed and compressed (gzip or bzip2) files are accepted for -c, -gc, -ont and -p options.\n"
              << "\n\n"

              << "Outputs:\n"
              << "    PREFIX_combined.fa\n"
              << std::endl;
}



//////////////////////////////////////////////////////////////////////////////////////
// main
//////////////////////////////////////////////////////////////////////////////////////
void Combine::exec(void)
{
	platanus::setGlobalTmpFileDir(optionSingleArgs["-tmp"].c_str());
	setSubBinPath(optionSingleArgs["-sub_bin"]);

    setIntermediateDirectoryName(optionSingleArgs["-o"] + "_combineIntermediateResults");
    struct stat st;
	if (stat(intermediateDirectoryName.c_str(), &st) != 0) {
		int returnValue = mkdir(intermediateDirectoryName.c_str(), 0755);
		if (returnValue != 0) {
			throw platanus::CreateDirError(intermediateDirectoryName);
		}
	}

	execCombineMode();
    if (!optionBool["-no_gap_close"]) {
		execCombinatorialGapClosePerl();
		execRemoveRedundantSeqPerl();
	}

    std::ostringstream oss;
    if (!optionBool["-no_gap_close"]) {
		oss << "mv " <<  intermediateDirectoryName << "/" << optionSingleArgs["-o"] << "_closed.fa.rmred"
			<< " " << optionSingleArgs["-o"] << "_combined.fa" << ";"
		;
	}
	else {
		oss << "mv " <<  intermediateDirectoryName << "/" << optionSingleArgs["-o"] << "_combined.fa"
			<< " " << optionSingleArgs["-o"] << "_combined.fa" << ";"
		;
	}

    if (!optionBool["-keep_file"]) {
		oss << " rm " << intermediateDirectoryName << "/*.fa;"
			<< " rm " << intermediateDirectoryName << "/*.bed;"
			<< " rm " << intermediateDirectoryName << "/*.tsv;"
			<< " rm -rf " << intermediateDirectoryName << "/*_intermediates;"
		;
	}

	std::cerr << oss.str() << std::endl;
    if (system(oss.str().c_str()) != 0) {
        throw platanus::CreateLinkError();
    }
}


void Combine::setSubBinPath(std::string path)
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


void Combine::execCombineMode(void)
{
	std::vector<std::string> multiArgOptionList = {"-c", "-gc", "-p", "-ont"};

	std::ostringstream oss;
	oss << this->platanusBinary  << " solve_DBG -combine";

	for (auto optItr = multiArgOptionList.begin(); optItr != multiArgOptionList.end(); ++optItr) {
		if (!(optionMultiArgs[*optItr].empty())) {
			oss << " " << *optItr;
			for (auto argItr = optionMultiArgs[*optItr].begin(); argItr != optionMultiArgs[*optItr].end(); ++argItr) {
				oss << " " << *argItr;
			}
		}
	}

	oss << " -t " << optionSingleArgs["-t"]
		<< " -tmp " << optionSingleArgs["-tmp"]
		<< " -combine_l " << optionSingleArgs["-combine_l"]
		<< " -combine_L " << optionSingleArgs["-combine_L"]
		<< " -combine_t " << optionSingleArgs["-combine_t"]
		<< " -combine_s " << optionSingleArgs["-combine_s"]
		<< " -combine_g " << optionSingleArgs["-combine_g"]
		<< " -combine_i " << optionSingleArgs["-combine_i"]
		<< " -o " << intermediateDirectoryName << "/" << optionSingleArgs["-o"]
		<< " > "  << intermediateDirectoryName << "/" << optionSingleArgs["-o"] << ".scafLog 2>&1"
	;

	std::cerr << oss.str() << std::endl;
    if (system(oss.str().c_str()) != 0) {
        throw platanus::ScaffoldError();
	}
}


void Combine::execCombinatorialGapClosePerl(void)
{
	std::vector<std::string> multiArgOptionList = {"-gc", "-p", "-ont"};

	std::ostringstream oss;
	oss << this->subBinPath  << "/combinatorial_gap_close.pl"
		<< " -t " << optionSingleArgs["-t"]
		<< " -b " << intermediateDirectoryName << "/" << optionSingleArgs["-o"] << "_combined.fa"
		<< " -c "
	;
	for (auto optItr = multiArgOptionList.begin(); optItr != multiArgOptionList.end(); ++optItr) {
		for (auto argItr = optionMultiArgs[*optItr].begin(); argItr != optionMultiArgs[*optItr].end(); ++argItr) {
			oss << " " << *argItr;
		}
	}

	oss  << " -o " << intermediateDirectoryName << "/" << optionSingleArgs["-o"]
		<< " > "  << intermediateDirectoryName << "/" << optionSingleArgs["-o"] << ".combinatorialGapCloseLog 2>&1"
	;

	std::cerr << oss.str() << std::endl;
    if (system(oss.str().c_str()) != 0) {
        throw platanus::SubBinError("combinatorial_gap_close.pl");
	}
}


void Combine::execRemoveRedundantSeqPerl(void)
{
	std::ostringstream oss;
	oss << this->subBinPath  << "/remove_redundant_seq.pl"
		<< " -t " << optionSingleArgs["-t"]
		<< " -o " << intermediateDirectoryName << "/" << optionSingleArgs["-o"]
		<< " " << intermediateDirectoryName << "/" << optionSingleArgs["-o"] << "_closed.fa"
		<< " > "  << intermediateDirectoryName << "/" << optionSingleArgs["-o"] << ".removeRedundantSeqLog 2>&1"
	;

	std::cerr << oss.str() << std::endl;
    if (system(oss.str().c_str()) != 0) {
        throw platanus::SubBinError("remove_redundant_seq.pl");
	}
}
