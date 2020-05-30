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

#include "kmer_divide.h"


//////////////////////////////////////////////////////////////////////////////////////
// default constructor
//////////////////////////////////////////////////////////////////////////////////////
ContigDivider::ContigDivider(): contig(), readLength(100), contigMaxK(100), medianCoverageAll(0), coverageRateThreshold(), maskCoverageRateThreshold(), occurencePointer(), breakPoint()
{
    optionSingleArgs["-o"] = "out";
    optionSingleArgs["-k"] = "";
    optionSingleArgs["-r"] = "0.1";
    optionSingleArgs["-R"] = "2.0";
    optionMultiArgs["-f"] = std::vector<std::string>();
    optionSingleArgs["-t"] = "1";
    optionSingleArgs["-tmp"] = ".";
    optionBool["-recalc_cov"] = false;
    optionBool["-mask"] = false;
}


//////////////////////////////////////////////////////////////////////////////////////
// usage
//////////////////////////////////////////////////////////////////////////////////////
void ContigDivider::usage(void) const
{
    std::cerr << "\nUsage: platanus_b kmer_divide [Options]\n"
              << "Options:\n"
              << "    -o STR               : prefix of output file (default " << optionSingleArgs.at("-o") << ", length <= " << platanus::ConstParam::MAX_FILE_LEN << ")\n"
              << "    -f FILE1 [FILE2 ...] : contig or scaffold_file (fasta format)\n"
              << "    -k STR               : k-mer-occurrence file (binary format)\n"
              << "    -r FLOAT             : coverage rate threshold to divide (divided if cov < FLOAT*median_cov; default " << optionSingleArgs.at("-r") << ")\n"
              << "    -R FLOAT             : coverage rate threshold to mask (masked if cov > FLOAT*overall_median_cov; effective with -mask; default " << optionSingleArgs.at("-R") << ")\n"
              << "    -recalc_cov          : re-falculate coverage depth, not divide sequences\n"
              << "    -mask                : mask sequences, not divide ones\n"
              << "    -t INT               : number of threads (<= " << platanus::ConstParam::MAX_THREAD << ", default " << optionSingleArgs.at("-t") << ")\n"
              << "    -tmp DIR             : directory for temporary files (default " << optionSingleArgs.at("-tmp") << ")\n\n\n"

              << "Output:\n"
              << "    PREFIX_kmerDivided.fa (when \"-recalc_cov\" not specified)\n"
              << "    PREFIX_kmerMasked.fa (when \"-mask\" not specified)\n"
              << "    PREFIX_recalc.fa (when \"-recalc_cov\" specified)\n"
              << std::endl;
}




//////////////////////////////////////////////////////////////////////////////////////
// exec divide
//////////////////////////////////////////////////////////////////////////////////////
void ContigDivider::exec(void)
{
    omp_set_num_threads(atoi(optionSingleArgs["-t"].c_str()));
	platanus::setGlobalTmpFileDir(optionSingleArgs["-tmp"].c_str());
    unsigned long long kmerLength = platanus::getKmerLengthFromBinary(optionSingleArgs["-k"]);

    this->readLength = platanus::Contig::getReadLengthFromFastaHeader(optionMultiArgs["-f"][0]);
	this->contigMaxK = contig.getMaxKFromFastaHeader(optionMultiArgs["-f"][0]);
	this->maskCoverageRateThreshold = atof(optionSingleArgs["-R"].c_str());
	this->coverageRateThreshold = atof(optionSingleArgs["-r"].c_str());
	if (!optionBool["-recalc_cov"]) {
		this->coverageRateThreshold = atof(optionSingleArgs["-r"].c_str());
	}
	else {
		this->coverageRateThreshold = DBL_MAX;
	}

    for (auto itr = optionMultiArgs["-f"].begin(), end = optionMultiArgs["-f"].end(); itr != end; ++itr) {
        this->contig.readFastaCoverage(*itr);
    }

    if (kmerLength <= 32) {
        Counter<Kmer31> counter(kmerLength);
        counter.readOccurrenceTableBinary(optionSingleArgs["-k"]);
        this->getOccurrenceArray(counter);
    } else if (kmerLength <= 64) {
        Counter<KmerN<Binstr63> > counter(kmerLength);
        counter.readOccurrenceTableBinary(optionSingleArgs["-k"]);
        this->getOccurrenceArray(counter);
    } else if (kmerLength <= 96) {
        Counter<KmerN<Binstr95> > counter(kmerLength);
        counter.readOccurrenceTableBinary(optionSingleArgs["-k"]);
        this->getOccurrenceArray(counter);
    } else if (kmerLength <= 128) {
        Counter<KmerN<Binstr127> > counter(kmerLength);
        counter.readOccurrenceTableBinary(optionSingleArgs["-k"]);
        this->getOccurrenceArray(counter);
    } else if (kmerLength <= 160) {
        Counter<KmerN<Binstr159> > counter(kmerLength);
        counter.readOccurrenceTableBinary(optionSingleArgs["-k"]);
        this->getOccurrenceArray(counter);
    } else {
        Counter<KmerN<binstr_t> > counter(kmerLength);
        counter.readOccurrenceTableBinary(optionSingleArgs["-k"]);
        this->getOccurrenceArray(counter);
    }

//	dumpKmerCoverage("out_kmer.depth");

	if (!optionBool["-mask"]) {
		this->setMedianCoverage();
		this->decideContigBreakPoint();

		std::string outputFilename = optionSingleArgs["-o"]; 
		if (!optionBool["-recalc_cov"]) {
			outputFilename += "_kmerDivided.fa";
		}
		else {
			outputFilename += "_recalc.fa";
		}
		this->divideAndPrintContig(outputFilename, kmerLength);
	}
	else {
		this->setMedianCoverageAll();
		this->maskHighCoverageKmer();

		std::string outputFilename = optionSingleArgs["-o"]; 
		outputFilename += "_kmerMasked.fa";
		this->printContig(outputFilename);
	}


	std::cerr << "divide completed" << std::endl;
}


//////////////////////////////////////////////////////////////////////////////////////
// exec divide
//////////////////////////////////////////////////////////////////////////////////////
template <typename KMER>
void ContigDivider::getOccurrenceArray(const Counter<KMER> &counter)
{
    unsigned long long kmerLength = counter.getKmerLength();
    const unsigned long long mask = kmerLength >= 32 ? ~static_cast<unsigned long long>(0) : ~(~static_cast<unsigned long long>(0) << (2 * kmerLength));
    KMER kmer(kmerLength);
    occurencePointer.reset(new OccurrenceArray(kmerLength, this->contig));

    for (unsigned long long i = 0; i < occurencePointer->getNumSeq(); ++i) {
        unsigned long long numKmer = occurencePointer->getNumKmer(i);
        if (numKmer == 0) continue;

        const platanus::SEQ &seq = contig.seq[i];
        unsigned long long start = 0;
        bool isInit = true;
        while (start < numKmer) {
            if (isInit) {
                unsigned long long j = 0;
                for (; j < kmerLength - 1; ++j) {
                    if (seq.base[start + j] == 4) break;

                    kmer.setForward(kmerLength - 2 - j, seq.base[start + j]);
                    kmer.setReverse(j + 1, 0x3 ^ seq.base[start + j]);
                }

                if (j == kmerLength - 1)
                    isInit = false;
                else {
                    start += j + 1;
                    continue;
                }
            }

            if (seq.base[start + kmerLength - 1] == 4) {
                start += kmerLength;
                isInit = true;
                continue;
            }
            kmer.forward <<= 2;
            kmer.maskForward(mask);
            kmer.reverse >>= 2;
            kmer.setForward(0, seq.base[start + kmerLength - 1]);
            kmer.setReverse(kmerLength - 1, 0x3 ^ seq.base[start + kmerLength - 1]);
            occurencePointer->setOccurrenceArray(i, start, std::max(counter.findValue(std::min(kmer.forward, kmer.reverse)), static_cast<unsigned short>(0)));
            ++start;
        }
    }
}


void ContigDivider::decideContigBreakPoint(void)
{
    const unsigned long long kmerLength = occurencePointer->getKmerLength();
    this->breakPoint.resize(this->contig.numSeq);

    # pragma omp parallel for schedule(dynamic)
    for (long i = 0; i < this->contig.numSeq; ++i) {
        std::vector<unsigned long long> occurrenceArray = occurencePointer->getOccurrenceArray(i);
        this->breakPoint[i].emplace_back(0);
        this->breakPoint[i].emplace_back(this->contig.seq[i].length - kmerLength + 1);

		const unsigned long long  breakCoverageCutoff = medianCoverage[i] * this->coverageRateThreshold;

		if (breakCoverageCutoff > 0) {
			for (unsigned long long j = 0; j < contig.seq[i].length - kmerLength + 1; ++j) {
				if (occurrenceArray[j] < breakCoverageCutoff) {
					this->breakPoint[i].emplace_back(j);
					this->breakPoint[i].emplace_back(j + 1);
				}
			}
		}
    }
}


void ContigDivider::maskHighCoverageKmer(void)
{
    const unsigned long long kmerLength = occurencePointer->getKmerLength();
    this->breakPoint.resize(this->contig.numSeq);
	const unsigned long long coverageThreshold = this->medianCoverageAll * this->maskCoverageRateThreshold;

	std::cerr << "masking high-coverage k-mers..." << std::endl;
	std::cerr << "OVERALL_MEDIAN_COVERAGE = " << this->medianCoverageAll << std::endl;
	std::cerr << "COVERAGE_THRESHOLD = " << coverageThreshold << std::endl;

    # pragma omp parallel for schedule(dynamic)
    for (long i = 0; i < this->contig.numSeq; ++i) {
        std::vector<unsigned long long> occurrenceArray = occurencePointer->getOccurrenceArray(i);
		for (unsigned long long j = 0; j < contig.seq[i].length - kmerLength + 1; ++j) {
			if (occurrenceArray[j] > coverageThreshold) {
				std::fill(this->contig.seq[i].base.begin() + j, this->contig.seq[i].base.begin() + j + kmerLength, 4);
			}
		}
    }
}


unsigned long long ContigDivider::findMedian(const std::vector<unsigned long long> &vec) const
{
    std::vector<unsigned long long> tmp = vec;
    std::partial_sort(tmp.begin(), tmp.begin() + tmp.size() / 2 + 1, tmp.end());
    return tmp[tmp.size() / 2];
}


void ContigDivider::setMedianCoverage(void)
{
    medianCoverage.assign(this->contig.numSeq, 0);
    for (long i = 0; i < this->contig.numSeq; ++i) {
		medianCoverage[i] = findMedian(occurencePointer->getOccurrenceArray(i));
	}
}


void ContigDivider::setMedianCoverageAll(void)
{
	std::vector<unsigned long long> buffer;
    for (long i = 0; i < this->contig.numSeq; ++i) {
		std::vector<unsigned long long> occurrenceArray = occurencePointer->getOccurrenceArray(i);
		std::copy(occurrenceArray.begin(), occurrenceArray.end(), std::back_inserter(buffer));
	}
	this->medianCoverageAll = findMedian(buffer);
}


void ContigDivider::divideAndPrintContig(const std::string &filename, const unsigned long long kmerLength)
{
    unsigned long long seqID = 1;
    std::ofstream ofs(filename.c_str());

    for (unsigned long long contigID = 0; contigID < this->contig.numSeq; ++contigID) {
		std::sort(breakPoint[contigID].begin(), breakPoint[contigID].end());
        for (auto itr = breakPoint[contigID].begin(), end = breakPoint[contigID].end(); itr != end; ++itr) {
            auto next = itr;
            ++next;
            if (next == end)
				break;

            if (*next - *itr <= 0)
				continue;

			const unsigned long long maxCoverageCutoff = std::max(medianCoverage[contigID] * this->coverageRateThreshold, 1.0);

			if (judgeMajorityGreaterOrEqualCoverage(contigID, *itr, *next, maxCoverageCutoff)) {
				unsigned long readCoverage = static_cast<unsigned long long>((this->calcAverageCoverage(contigID, *itr, *next) * this->readLength) / (this->readLength - kmerLength + 1) + 0.5);

				ofs << ">seq" << seqID << "_len" << *next - *itr + kmerLength - 1 << "_cov" << readCoverage << "_read" << this->readLength << "_maxK" << this->contigMaxK << "\n";

				unsigned long long j = 0;
				for (unsigned long long pos = *itr; pos < *next + kmerLength - 1; ++pos, ++j) {
					ofs.put(platanus::Bin2Char(this->contig.seq[contigID].base[pos]));
					if ((j + 1) % platanus::ConstParam::OUTPUT_LINE_LENGTH == 0)
						ofs.put('\n');
				}
				if (j % platanus::ConstParam::OUTPUT_LINE_LENGTH != 0)
					ofs.put('\n');
			}
        }
        ++seqID;
    }
    ofs.close();
}


void ContigDivider::printContig(const std::string &filename)
{
    std::ofstream ofs(filename.c_str());

    for (long i = 0; i < this->contig.numSeq; ++i) {
		ofs << '>' << this->contig.name[i] << '\n';

		std::string &seq = this->contig.seq[i].base;
		long length = this->contig.seq[i].base.size();
		long j = 0;
		for (; j < length; ++j) {
			ofs.put(platanus::Bin2Char(seq[j]));
			if ((j + 1) % platanus::ConstParam::OUTPUT_LINE_LENGTH == 0)
				ofs.put('\n');
		}
		if (j % platanus::ConstParam::OUTPUT_LINE_LENGTH != 0)
			ofs.put('\n');
    }
    ofs.close();
}


unsigned long long ContigDivider::calcAverageCoverage(const unsigned long long contigID, const unsigned long long start, const unsigned long long end)
{
    double coverage = 0;
    for (unsigned long long i = start; i < end; ++i)
        coverage += occurencePointer->getOccurrence(contigID, i);

    coverage /= (end - start);
    return static_cast<unsigned long>(coverage + 0.5);
}


unsigned long long ContigDivider::calcMaxCoverage(const unsigned long long contigID, const unsigned long long start, const unsigned long long end)
{
    unsigned long long max = 0;
    unsigned long long kmerLength = occurencePointer->getKmerLength();
    for (unsigned long long i = start; i < end - kmerLength + 1; ++i) {
        if (max < occurencePointer->getOccurrence(contigID, i))
			max = occurencePointer->getOccurrence(contigID, i);
    }

    return max;
}


bool ContigDivider::judgeMajorityGreaterOrEqualCoverage(const unsigned long long contigID, const unsigned long long start, const unsigned long long end, const unsigned long long threshold)
{
	unsigned long long n = 0;
    for (unsigned long long i = start; i < end; ++i) {
        if (threshold <= occurencePointer->getOccurrence(contigID, i)) {
			++n;
			if (n >= (end - start) / 2)
				return true;
		}
    }
	return false;
}


void ContigDivider::dumpKmerCoverage(const std::string outFile)
{
    const unsigned long long kmerLength = occurencePointer->getKmerLength();

    std::ofstream out(outFile);

    for (unsigned long long i = 0; i < this->contig.numSeq; ++i) {
		out << '>' << this->contig.name[i] << '\n';
        for (unsigned long long j = 0; j < contig.seq[i].length - kmerLength + 1; ++j) {
           out << occurencePointer->getOccurrence(i, j) << '\n';
        }
    }

	out.close();
}
