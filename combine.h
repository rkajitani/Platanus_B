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

#ifndef COMBINE_H
#define COMBINE_H

#include "baseCommand.h"
#include "scaffold.h"
#include "gapClose.h"
#include "merge.h"
#include "polish.h"
#include <memory>
#include <sstream>
#include <vector>

class Combine : public BaseCommand
{
private:
    std::string intermediateDirectoryName;
    std::string platanusBinary;
    std::string subBinPath;


    virtual bool checkFileEnough(void)
    {
        if (optionMultiArgs["-c"].empty() || (optionMultiArgs["-gc"].empty() && optionMultiArgs["-p"].empty() && optionMultiArgs["-ont"].empty())) {
            std::cerr << "Error: contig fasta files were specified; -c and (-gc or -p or -ont) are required." << std::endl;
            return false;
        }
        return true;
    }

    virtual int checkOtherOption(char *argv) const
    {
        return 0;
    }


public:
    Combine();
    Combine(const Combine &) = delete;
    Combine &operator=(const Combine &) = delete;
    ~Combine() = default;

    void usage(void) const;
    void exec(void);

    void setIntermediateDirectoryName(const std::string name) { intermediateDirectoryName = name; }

    virtual bool parseArgs(int argc, char **argv)
    {
        if (BaseCommand::parseArgs(argc, argv)) {
            this->platanusBinary = argv[0];
            return true;
        }
        return false;
    }

	void setSubBinPath(std::string path);
	void execCombineMode(void);
	void execCombinatorialGapClosePerl(void);
	void execRemoveRedundantSeqPerl(void);
};



#endif
