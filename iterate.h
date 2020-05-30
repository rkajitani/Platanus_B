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

#ifndef ITERATE_H
#define ITERATE_H

#include "baseCommand.h"
#include "scaffold.h"
#include "gapClose.h"
#include "merge.h"
#include "polish.h"
#include <memory>
#include <sstream>
#include <vector>

class IterateScaffold : public BaseCommand
{
private:
    static const std::string CONTIG_FOOTER;
	static const std::string KMER_DIVIDE_FOOTER;
    static const std::string DIVIDE_FOOTER;
    static const std::string SCAF_FOOTER;
    static const std::string POLISH_FOOTER;
    static const std::string GAP_FOOTER;
    static const std::string EX_FOOTER;
    static const std::string DIV_FOOTER;
    static const std::string MERGE_FOOTER;
    static const std::string ITERATION_FOOTER;
    std::string directoryName;
    std::string previousDirectoryName;
    std::string intermediateDirectoryName;
    std::string platanusBinary;
    std::string subBinPath;


    virtual bool checkFileEnough(void)
    {
        if (optionMultiArgs["-c"].empty()) {
            std::cerr << "Error: No contig fasta files (-c) were specified!" << std::endl;
            return false;
        }
        return true;
    }

    virtual int checkOtherOption(char *argv) const
    {
        return 0;
    }


public:
    IterateScaffold();
    IterateScaffold(const IterateScaffold &) = delete;
    IterateScaffold &operator=(const IterateScaffold &) = delete;
    ~IterateScaffold() = default;

    void usage(void) const;
    void exec(void);

    void setIntermediateDirectoryName(const std::string name) { intermediateDirectoryName = name; }

    void setDirectoryName(const std::string intDirName, const long times)
    {
        std::ostringstream oss;
        oss <<  intDirName << "/" << optionSingleArgs["-o"] << times;
        previousDirectoryName = directoryName;
        directoryName = oss.str();
    }

    virtual bool parseArgs(int argc, char **argv)
    {
        if (BaseCommand::parseArgs(argc, argv)) {
            this->platanusBinary = argv[0];
            return true;
        }
        return false;
    }


    void createDirectory(void) const;
	void setSubBinPath(std::string path);
	void countKmer(void);
    void createContig(const long times);
    void execKmerDivideMode(void);
    void execDivideMode(void);
    void execMergeMode(const long times);
    void execScaffoldMode(const long times, const long numTimes);
    void execPolishMode(void);
	void execFinalDivideMode(void);
	void execFinalPolishMode(void);
    void execGapCloseMode(const long times, const long numTimes);
	void execClusterMode(void);
	void execPostClusteringScaffoldMode(void);
	void execCombineMode(const long times);
	void execCombinatorialGapClosePerl(const long times);
	void execRemoveRedundantSeqPerl(const long times);
};



#endif
