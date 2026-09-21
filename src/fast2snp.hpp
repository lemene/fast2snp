#pragma once

#include <algorithm>
#include <bitset>
#include <cassert>
#include <chrono>
#include <cmath>
#include <ctime>
#include <fstream>
#include <iostream>
#include <iterator>
#include <stdlib.h>
#include <omp.h>
#include <sstream>
#include <string>
#include <thread>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <mutex>
#include <atomic>

#include "utils.hpp"
#include "options.hpp"
#include "graph_index.hpp"
#include "parallel_file_reader.hpp"




static const char ACGT[] = "ACGT";


class KmerProcessor {
public:
    // Bubble Walker Logic
    struct WalkState {
        enum State : uint8_t {
            Active = 0,
            DeadEnd,
            MERGED
        };
        std::string sequence;
        std::vector<uint16_t> kmerCounts;
        KmerInt currentK1mer;
        uint8_t k;
        uint8_t br;
        State state = State::Active;
        size_t totalDelta = 0;
        WalkState(std::string seq, std::vector<uint16_t> counts, KmerInt k1mer)
            : sequence(std::move(seq)), kmerCounts(std::move(counts)), currentK1mer(k1mer) { 
            k = sequence.size(); 
            br = SequenceUtils::kmerCharToInt(sequence.back());
        }

        WalkState() = default;
        WalkState(const WalkState& w) = default;

        bool isActive() const {
            return state == State::Active;
        }

        bool isMerged() const {
            return state == State::MERGED;
        }

        void pushKmerCount(uint16_t count) {
            if (kmerCounts.size() > 0) {
                size_t delta = std::abs((int)count - kmerCounts.back());

                totalDelta += std::min<size_t>(delta, 30);
            }
            kmerCounts.push_back(count);
        }


        void forward(uint8_t base, const KmerProcessor& graph, std::unordered_set<KmerInt, UllHash>& visited) {
            
            auto edgeCount = graph.getEdgeCount(currentK1mer, ACGT[base]);
            sequence += ACGT[base];
            pushKmerCount(edgeCount);
            currentK1mer = SequenceUtils::shiftRight(currentK1mer, k-1, ACGT[base]);

            assert(isActive());

            if (sequence.size() >= 2*(size_t)k-1) {
                visited.insert(currentK1mer);
            }
        }
        std::string countsAsString() const {
            return countsAsString(kmerCounts.size());
        }
        std::string countsAsString(size_t size) const {
            assert(size <= kmerCounts.size());

            std::stringstream ss;
            for (size_t i = 0; i < size; ++i) {
                ss << kmerCounts[i];
                if (i < size - 1) {
                    ss << ",";
                }
            }
            return ss.str();
        }


        double countsValidRate(size_t len, size_t low, size_t high) const {
            assert(len <= kmerCounts.size());
            size_t count = 0;
            for (size_t i = 0; i < len; ++i) {
                auto c = kmerCounts[i];
                if (c >= low && c <= high) {
                    count++;
                }
            }
            return count * 1.0 / len;
        }

        uint16_t kmerCount(size_t index) const {
            return kmerCounts[index];
        }

        bool isValidBranch(size_t len, const Options &opts) const {
            return countsValidRate(len, opts.lowCov, opts.highCov) >= opts.minHetRatio &&
               isBranchSmooth(len, opts) && isBranchHetSmooth(len, opts);
        }

        bool isBranchSmooth(size_t len, const Options &opts) const;


        bool isBranchSmooth_old(size_t len, const Options &opts) const {
            assert(len <= kmerCounts.size());
            
            if (!isBranchHetSmooth(len, opts)) return false;
            if (branchSteps(len, opts) > 3) return false;

            return true;
        }

        bool isBranchHetSmooth(size_t len, const Options &opts) const;
        size_t branchSteps(size_t len, const Options &opts) const;
    };

    struct Branch {
        std::vector<WalkState> walkers;
        std::unordered_set<KmerInt, UllHash> visited;
        void forward(const KmerProcessor& graph) {

            size_t originalSize = walkers.size();
            for (size_t i = 0; i < originalSize; ++i) {
                // Do not use aliases.
                if (!walkers[i].isActive()) continue;

                auto nexts = graph.getNextBases(walkers[i].currentK1mer);
                if (nexts.size() >= 1) {

                    size_t sub_origin_size = walkers.size();
                    for (size_t j = 1; j < nexts.size(); ++j) {
                        walkers.push_back(walkers[i]);   // copy of next.size() - 1 walker
                    }

                    walkers[i].forward(nexts[0], graph, visited);

                    for (size_t j = 1;  sub_origin_size + j - 1 < walkers.size(); j++) {
                        assert(j < walkers.size());
                        walkers[sub_origin_size + j - 1].forward(nexts[j], graph, visited);
                    }
                } else {
                    walkers[i].state = WalkState::DeadEnd;
                }
            }
        }
        bool isMerged(const Branch& other) {
            for (auto& w : walkers) {
                if (other.visited.find(w.currentK1mer) != other.visited.end()) {
                    return true;
                }
            }   
            return false;
        }

        void resize(const Options& opts) {
            if (walkers.size() > opts.MAX_WALKERS_SIZE) {
                std::nth_element(walkers.begin(), walkers.begin()+opts.MAX_WALKERS_SIZE, walkers.end(), 
                    [](const WalkState& a, const WalkState& b) { return a.totalDelta < b.totalDelta; });
                walkers.resize(opts.MAX_WALKERS_SIZE);
            }
        }
    };
    struct WalkerPair {
        WalkerPair() {}
        WalkerPair(size_t k, const WalkState& walkerA, const WalkState& walkerB) : k_(k), walkerA_(&walkerA), walkerB_(&walkerB) {

        }

        void findBestMatch(const Options& opts) { 
            
            std::string tailA = walkerA_->sequence.substr(k_-1);
            std::string tailB = walkerB_->sequence.substr(k_-1);

            for (size_t i = 0; i < tailA.length(); ++i) {
                for (size_t j = 0; j < tailB.length(); ++j) {
                    size_t currentLen = 0;
                    size_t idxA = i;
                    size_t idxB = j;
                    
                    while (idxA < tailA.length() && idxB < tailB.length()) {
                        if (tailA[idxA] == tailB[idxB]) {
                            currentLen++;
                            idxA++;
                            idxB++;
                        } else {
                            break; // Stop at first mismatch
                        }
                    }

                    if (currentLen > matchLen_) {
                        matchLen_ = currentLen;
                        offsetA_ = i;
                        offsetB_ = j;
                    }
                }
            }
            // Traceback
            for (; matchLen_ > 0; --matchLen_) { 
                if (walkerA_->isValidBranch(matchLen_ + offsetA_, opts) && walkerB_->isValidBranch(matchLen_ + offsetB_, opts)) { 
                    break;
                }
            }

        }


        int matchScore() const { 
            if (walkerA_ == nullptr || walkerB_ == nullptr) return 0;
            return (int)matchLen_ * 2 - (int)varSize();
        }
        size_t varSize() const {
            if (walkerA_ == nullptr || walkerB_ == nullptr) return INT32_MAX;
            return offsetA_ + offsetB_;
        }

        bool isValid(const Options& opts) {
            return matchLen_ >= opts.extend && 
                walkerA_->isValidBranch(offsetA_+matchLen_, opts) && 
                walkerB_->isValidBranch(offsetB_+matchLen_, opts);
        }

        std::string type() const {
            if (offsetA_ == 1 && offsetB_ == 1) {
                return "IsolatedSNP";
            } else if (offsetA_ != offsetB_) {
                return "INDEL";
            } else {
                return "NonIsolatedSNP";
            }
        }

        std::string toString(size_t id, const Options& opts) const {
            std::ostringstream oss;

            assert(matchLen_ >= opts.extend);
            size_t ext = std::min<size_t>(std::max<size_t>(opts.k - 1, opts.extend), opts.extend);

            oss << ">Bubble_" << id << "_" << type() << "_Path1_Ext" << ext
                    << " KmerCounts=" << walkerA_->countsAsString(offsetA_ + ext) << "\n";
            oss << walkerA_->sequence.substr(0, k_-1+offsetA_ + ext) << "\n";
            oss << ">Bubble_" << id << "_" << type() << "_Path2_Ext" << ext
                    << " KmerCounts=" << walkerB_->countsAsString(offsetB_ + ext) << "\n";
            oss << walkerB_->sequence.substr(0, k_-1+offsetB_ + ext) << "\n";
            return oss.str();
        }
    //protected:
        size_t k_ { 0 };
        const WalkState* walkerA_ { nullptr };
        const WalkState* walkerB_ { nullptr };
        size_t offsetA_ { 0 };
        size_t offsetB_ { 0 };
        size_t matchLen_ { 0 };
    };



private:
    const Options& opts_;
    GraphIndex graphIndex;

public:
    KmerProcessor(const Options& opts) : opts_(opts) { }
  
    void loadKmersFromTxt();

    // Optimized processSingleKmer that takes pre-converted KmerInt
    void processSingleKmerVal(KmerInt val, uint16_t count) ;
    
    // Helper to insert (k-1)-mer and its edge-count into the index
    // direction: 0 for Right, 1 for Left
    void insertToIndex(KmerInt kmerVal, char edgeBase, int direction, uint16_t count);

    uint16_t getEdgeCount(KmerInt k1mer, char base) const;
    std::vector<uint8_t> getNextBases(KmerInt k1mer) const;



    void computeKmerDeltaStats();

    std::vector<KmerInt> findBubbleSeed();
    void saveSeeds(const std::string& fname, const std::vector<KmerInt>& seeds);
    void findAndWalkBubbles();
    void branchForward(Branch& branch) const;
    WalkerPair detectBestMatch(const std::vector<Branch>& branchs) const ;
};
