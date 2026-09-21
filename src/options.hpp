#pragma once

#include <iostream>
#include <getopt.h>
#include <vector>

#include "utils.hpp"

struct Options {
    Options() {}
    bool parse(int argc, const char* const* argv) {
        static struct option longOptions[] = {
            {"input", required_argument, 0, 'i'},
            {"left", required_argument, 0, 'l'},
            {"right", required_argument, 0, 'r'},
            {"min-count", required_argument, 0, 'c'},   // 新增：加载 k-mer 的最小频次
            {"max-count", required_argument, 0, 'C'},   // 新增：加载 k-mer 的最大频次
        //    {"kmer", required_argument, 0, 'k'},
            {"thread", required_argument, 0, 't'},
            {"output", required_argument, 0, 'o'},
            {"extend", required_argument, 0, 'e'},
            {"fluct_threshold", required_argument, 0, 'f'},
            {"debug", required_argument, 0, 'd'},
            {"help", no_argument, 0, 'h'},
            {0, 0, 0, 0}
        };

        int opt;
        int optionIndex = 0;
        while ((opt = getopt_long(argc, const_cast<char *const *>(argv), "i:l:r:c:C:t:o:e:f:d:h", longOptions, &optionIndex)) != -1) {
            switch (opt) {
                case 'i': inputFile = optarg; break;
                case 'l': lowCov = std::stoi(optarg); break;
                case 'r': highCov = std::stoi(optarg); break;
                case 'c': loadLowCov = std::stoi(optarg); break;
                case 'C': loadHighCov = std::stoi(optarg); break;
        //        case 'k': k = std::stoi(optarg); break;
                case 't': threads = std::stoi(optarg); break;
                case 'o': outputPath = SequenceUtils::ensureTrailingSlash(optarg); break;
                case 'e': extend = std::stoi(optarg); break;
                case 'f': deltaThreshold = std::stod(optarg); break;
                case 'd': debugKmer = optarg; break;
                case 'h':
                    std::cout << "Usage: fast2snp [options]\n";
                    return false;
            }
        }

        if (inputFile.empty() || lowCov == 0 || highCov == 0) {
            std::cerr << "Missing required arguments.\n";
            return false;
        }

        // other options
        //loadHighCov =  highCov;//10000;
        if (loadHighCov <= 0) loadHighCov = highCov;
        if (loadLowCov <= 0) loadLowCov = lowCov;
        std::cout << "Coverage: " << lowCov << " - " << highCov << " " << loadLowCov << "-" << loadHighCov << "\n";
        return true;
    }

    std::string getOutput() const { return outputPath + "bubbles.snp"; }
    std::string getDebug() const { return outputPath + "debug.txt"; }   
    std::string getCandSeedsFile() const { return outputPath + "cand_seeds.txt"; }   
    std::string getSeedsFile() const { return outputPath + "seeds.txt"; }   
    size_t getMaxSteps() const { return k + k + variantLen; }
    bool isHetCov(int cov) const { return cov >= lowCov && cov <= highCov; } 

    void setK(int k) const { 
        this->k = k; 
        if (extend <= 0) {
            extend = k - 1;
        }
    }
    bool isDebug() const { return !debugKmer.empty(); }
    std::string inputFile;
    
    int loadLowCov = 0;
    int loadHighCov = 0;
    int lowCov = 0;
    int highCov = 0;
    int threads = 1;
    std::string outputPath = "./";
    mutable int extend = -1;
    std::string debugKmer = "";
    const size_t MAX_WALKERS_SIZE = 20;
    int variantLen = 100;
    double minHetRatio = 0.5;

    mutable int k = -1;
    mutable double deltaThreshold = 0.0;
};
