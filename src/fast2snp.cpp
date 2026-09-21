#include "fast2snp.hpp"

#include <numeric> // accumulate

bool KmerProcessor::WalkState::isBranchHetSmooth(size_t len, const Options &opts) const { 
    double last = kmerCounts[0];
    double halfWin = (opts.highCov - opts.lowCov ) / 2;
    double lowCov = last - halfWin;
    double highCov = last + halfWin;
    assert(last >= opts.lowCov && last <= opts.highCov);
    for (size_t i = 1; i < len; ++i) {
        auto c = kmerCounts[i];
        if (c >= lowCov && c <= highCov) {
            auto d = std::abs(last - c) ;
            if (d > opts.deltaThreshold) {
                return false;
            }
            last = c;
        }
    }
    return true;
}

size_t KmerProcessor::WalkState::branchSteps(size_t len, const Options &opts) const {
    size_t steps = 0;
    double last = kmerCounts[0];
    for (size_t i = 1; i < len; ++i) {
        auto c = kmerCounts[i];
        double d = std::abs(last - c);
        double mid = (last + c) / 2;
        double threshold = mid <= opts.highCov ? opts.deltaThreshold : std::ceil(opts.deltaThreshold * double(mid / opts.highCov));
        if (d > threshold) {
            steps++;
        }
        last = c;
    }
    return steps;
}

bool KmerProcessor::WalkState::isBranchSmooth(size_t len, const Options &opts) const {
    size_t down_count = 0;
    size_t up_count = 0;
    bool down = false;
    double last = kmerCounts[0];
    for (size_t i = 1; i < len; ++i) {
        auto c = kmerCounts[i];
        
        double mid = (last + c) / 2;
        double threshold = mid <= opts.highCov ? opts.deltaThreshold : std::ceil(opts.deltaThreshold * double(mid / opts.highCov));
        
        double d = c - last;
        if (d < -threshold) {
            down_count++;
            down = true;
        } else if (d > threshold) {
            up_count++;
        }
        if (opts.isHetCov(c)) {
            down = false;
        } 
        if (opts.isDebug()) {
            std::cout << "i: " << i << " last: " << last << " c: " << c << " d: " << d << " down_count: " << down_count << " up_count: " << up_count << std::endl;
        }
        last = c;
    }
    if (opts.isDebug()) {
        std::cout << "len: " << len << " down: " << down << " down_count: " << down_count << " up_count: " << up_count << std::endl;
    }
    return !down;
}

void KmerProcessor::loadKmersFromTxt() {
    std::atomic<size_t> processedCount {0};

    auto getKmerLen = [this](const std::string& chunk) {
        std::istringstream iss(chunk);
        std::string line;
        while (std::getline(iss, line)) {
            std::istringstream iss(line);
            std::string kmer;
            int count;
            if (iss >> kmer >> count) {
                opts_.setK(kmer.length());
                break;
            }
        }
    };

    auto parseKmers = [this, &processedCount](const std::string& chunk) {
        std::istringstream iss(chunk);
        std::string line;
        while (std::getline(iss, line)) {
            std::istringstream iss(line);
            std::string kmer;
            int count;
            if (iss >> kmer >> count) {
                if (count >= opts_.loadLowCov && count <= opts_.loadHighCov) {
                    processSingleKmerVal(SequenceUtils::kmerToInt(kmer), count > UINT16_MAX ? UINT16_MAX : static_cast<uint16_t>(count));
                    processedCount.fetch_add(1);
                }
            }
        }
    };

    parallel_process(opts_.inputFile, opts_.threads, getKmerLen, parseKmers);
    
    std::cout << "Kmer Length " << opts_.k << std::endl;
    std::cout << "Processed " << processedCount.load() << " k-mers." << std::endl;
    std::cout << "Graph index size: " << graphIndex.size() << std::endl;
}


void KmerProcessor::processSingleKmerVal(KmerInt val, uint16_t count) {
    // CLEANUP: Mask the value to ensure no garbage in high bits
    KmerInt kmerMask = SequenceUtils::maskForBits(2 * opts_.k);
    val &= kmerMask;
    
    // 1. Split into Prefix (k-1) and Right Base
    KmerInt prefixVal = val >> 2;
    char rightBase = SequenceUtils::getBaseFromBitIndex(static_cast<int>(val & (KmerInt)0x3));
    insertToIndex(prefixVal, rightBase, 0, count); // 0 for Right

    // 2. Split into Suffix (k-1) and Left Base
    KmerInt suffixMask = SequenceUtils::maskForBits(2 * (opts_.k - 1));
    KmerInt suffixVal = val & suffixMask;
    char leftBase = SequenceUtils::getBaseFromBitIndex(static_cast<int>((val >> (2 * (opts_.k - 1))) & (KmerInt)0x3));
    insertToIndex(suffixVal, leftBase, 1, count); // 1 for Left
}

void KmerProcessor::insertToIndex(KmerInt kmerVal, char edgeBase, int direction, uint16_t count) {
    int km1 = opts_.k - 1;
    // Compute canonical form of the (k-1)-mer
    KmerInt canonicalVal = SequenceUtils::reverseInt(kmerVal, km1 * 2);
    int slot = 0;
    if (kmerVal < canonicalVal) {
        // Canonical form is smaller, so use it
        slot = baseToCountSlot(edgeBase);
        if (direction == 1) {
            slot += 4;
        }
    } else if (kmerVal > canonicalVal) {
        kmerVal = canonicalVal; 
        edgeBase = SequenceUtils::complement(edgeBase);
        slot = baseToCountSlot(edgeBase);
        if (direction == 0) {
            slot += 4;
        }
    } else if (kmerVal == canonicalVal) {
        if (direction == 1) {
            edgeBase = SequenceUtils::complement(edgeBase);
        }
        slot = baseToCountSlot(edgeBase);
    }

    PackedEdgeCounts initialValue = 0;
    setPackedCount(initialValue, slot, count);
    graphIndex.try_emplace_l(
        kmerVal,
        [&](auto& val) { setPackedCount(val.second, slot, count); },
        initialValue
    );
}

uint16_t KmerProcessor::getEdgeCount(KmerInt k1mer, char base) const {
    int k1 = opts_.k - 1;
    KmerInt canonicalVal = SequenceUtils::reverseInt(k1mer, k1 * 2);

    int slot = -1;
    if (k1mer <= canonicalVal) {
        slot = baseToCountSlot(base);
    } else {
        k1mer = canonicalVal;
        slot = baseToCountSlot(SequenceUtils::complement(base)) + 4;
    }

    auto it = graphIndex.find(k1mer);
    assert(it != graphIndex.end());
    assert(slot >= 0 && slot < 8);

    return getPackedCount(it->second, slot);
}


std::vector<uint8_t> KmerProcessor::getNextBases(KmerInt k1mer) const {
    std::vector<uint8_t> bases;
    int k1 = opts_.k - 1;
    KmerInt canonicalVal = SequenceUtils::reverseInt(k1mer, k1 * 2);
    bool flipped = false;
    
    if (k1mer <= canonicalVal) {
        canonicalVal = k1mer;
    } else {
        flipped = true;
        // k1mer <= canonicalVal 不需要处理 
    }
    
    auto it = graphIndex.find(canonicalVal);
    assert(it != graphIndex.end());
    
    PackedEdgeCounts packedCounts = it->second;
    if (!flipped) {
        if (getPackedCount(packedCounts, 0) > 0) bases.push_back(0);
        if (getPackedCount(packedCounts, 1) > 0) bases.push_back(1);
        if (getPackedCount(packedCounts, 2) > 0) bases.push_back(2);
        if (getPackedCount(packedCounts, 3) > 0) bases.push_back(3);
    } else {
        if (getPackedCount(packedCounts, 4) > 0) bases.push_back(3);
        if (getPackedCount(packedCounts, 5) > 0) bases.push_back(2);
        if (getPackedCount(packedCounts, 6) > 0) bases.push_back(1);
        if (getPackedCount(packedCounts, 7) > 0) bases.push_back(0);
    }
    
    return bases;
}
void KmerProcessor::computeKmerDeltaStats() { 

    std::vector<double> deltas;
    deltas.reserve(graphIndex.size());

    for (const auto& iter : graphIndex) {
        // KmerInt key = iter.first;
        PackedEdgeCounts packedCounts = iter.second;

        int maxLeftCount = getPackedCount(packedCounts, 0);
        uint16_t nonZeroLeft = 0;
        for (int slot = 0; slot < 4; ++slot) {
            auto count = getPackedCount(packedCounts, slot);
            nonZeroLeft += (count > 0);
            if (count > maxLeftCount) {
                maxLeftCount = count;
            }
        }

        int maxRightCount = getPackedCount(packedCounts, 4);
        uint16_t nonZeroRight = 0;
        for (int slot = 4; slot < 8; ++slot) {
            uint16_t count = getPackedCount(packedCounts, slot);
            nonZeroRight += (count > 0);
            if (count > maxRightCount) {
                maxRightCount = count;
            }   
        }

        if (nonZeroLeft == 1 && nonZeroRight == 1 && opts_.isHetCov(maxLeftCount) && opts_.isHetCov(maxRightCount)) {
            //deltas.push_back((double)std::abs(maxLeftCount - maxRightCount) * 2 / (double)(maxLeftCount + maxRightCount));
            deltas.push_back((double)std::abs(maxLeftCount - maxRightCount) );
        }
    }
    // compute mean variance of deltas
    double deltaMean = std::accumulate(deltas.begin(), deltas.end(), 0.0) / deltas.size();
    auto varianceDelta = std::accumulate(deltas.begin(), deltas.end(), 0.0, [deltaMean](double sum, uint16_t delta) { 
        return sum + std::pow(delta - deltaMean, 2); 
    }) / deltas.size();

    double deltaStdDev = std::sqrt(varianceDelta);
    opts_.deltaThreshold = std::ceil(deltaMean + deltaStdDev * 3);
    std::cout << "delta Mean: " << deltaMean << " StdDev: " << deltaStdDev  << " Threshold: " << opts_.deltaThreshold << std::endl;
}

std::vector<KmerInt> KmerProcessor::findBubbleSeed() {

    // Statistics for branch counts
    std::array<size_t, 3> branchCounts = {0, 0, 0}; // Index 0: 2-branches, Index 1: 3-branches, Index 2: 4+ branches
    
    std::vector<KmerInt> seeds;
    for (const auto& kv : graphIndex) {
        KmerInt key = kv.first;
        PackedEdgeCounts packedCounts = kv.second;

        KmerInt canonicalVal = SequenceUtils::reverseInt(key, (opts_.k-1)*2);
        assert(key <= canonicalVal);
        
        // Right branches (slots 0-3)， Left branches (slots 4-7)
        int rightBranches = countValidSlots(packedCounts, 0, 4, opts_.lowCov, opts_.highCov); // TODO
        int leftBranches = countValidSlots(packedCounts, 4, 8, opts_.lowCov, opts_.highCov);
        if (rightBranches >= 2 ) {
            assert(rightBranches <= 4);
            seeds.push_back(key);
            branchCounts[rightBranches-2]++;
        }
        
        if (leftBranches >= 2) {
            assert(leftBranches <= 4);
            KmerInt rcKey = SequenceUtils::reverseInt(key, (opts_.k-1)*2);
            seeds.push_back(rcKey);
            branchCounts[leftBranches-2]++;
        }
    }
    
    std::cout << "Found " << seeds.size() << " potential divergence points." << std::endl;
    size_t totalSeeds = branchCounts[0] + branchCounts[1] + branchCounts[2];
    if (totalSeeds > 0) {
        std::cout << "Branching Statistics:" << std::endl;
        std::cout << "  2-branches: " << branchCounts[0] << " (" << (branchCounts[0] * 100.0 / totalSeeds) << "%)" << std::endl;
        std::cout << "  3-branches: " << branchCounts[1] << " (" << (branchCounts[1] * 100.0 / totalSeeds) << "%)" << std::endl;
        std::cout << "  4-branches: " << branchCounts[2] << " (" << (branchCounts[2] * 100.0 / totalSeeds) << "%)" << std::endl;
    }
    return seeds;
}

void KmerProcessor::saveSeeds(const std::string& fname, const std::vector<KmerInt>& seeds) {
    std::ofstream seedOut(fname);
    
    if (!seedOut.is_open()) {
        std::cerr << "Error: Could not open file for writing: " << opts_.getSeedsFile() << std::endl;
        return;
    }

    for (const auto& key : seeds) {
        auto it = graphIndex.find(key);
        if (it == graphIndex.end()) {
            it = graphIndex.find(SequenceUtils::reverseInt(key, (opts_.k-1)*2));
        }
        assert(it != graphIndex.end());

        // 计算左右分支数
        PackedEdgeCounts packedCounts = it->second;
        int rightBranches = countValidSlots(packedCounts, 0, 4, opts_.lowCov, opts_.highCov);
        int leftBranches = countValidSlots(packedCounts, 4, 8, opts_.lowCov, opts_.highCov);
        
        // 获取 (k-1)-mer 的字符串表示用于可读性
        std::string seedSeq = SequenceUtils::kmerIntToString(key, opts_.k - 1);
        KmerInt rcKey = SequenceUtils::reverseInt(key, (opts_.k-1)*2);
        std::string reSeedSeq = SequenceUtils::kmerIntToString(rcKey, opts_.k - 1);
        
        seedOut << seedSeq << "\t" << reSeedSeq << "\t" << rightBranches << "\t" << leftBranches << std::endl;
    }
    
    seedOut.close();
    std::cout << "Saved " << seeds.size() << " seeds to " << opts_.getSeedsFile() << std::endl;
}

void KmerProcessor::findAndWalkBubbles() {
    std::ofstream snpOut(opts_.getOutput());
    std::ofstream debugOut(opts_.getDebug());

    std::atomic<size_t> bubbleCount(0);
    std::atomic<size_t> isolatedSnpCount(0);
    std::atomic<size_t> nonIsolatedSnpCount(0);
    std::atomic<size_t> indelCount(0);
    
    if (opts_.deltaThreshold <= 0)   computeKmerDeltaStats();

    std::vector<KmerInt> seeds = findBubbleSeed();
    saveSeeds(opts_.getSeedsFile(), seeds);

    #pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < seeds.size(); ++i) {
        
        KmerInt seedKey = seeds[i];
        std::string seedSeq = SequenceUtils::kmerIntToString(seedKey, opts_.k - 1);

        if (!opts_.debugKmer.empty() && seedSeq != opts_.debugKmer && seedSeq != SequenceUtils::reverseString(opts_.debugKmer)) {
            continue;
        }

        // Initial branches
        // Create a walker for each branch
        auto nextBases = getNextBases(seedKey);
        assert(nextBases.size() >= 2);
        std::vector<Branch> branches;
        
        for (auto base : nextBases) {
            std::string seq = seedSeq;
            seq += ACGT[base];
            auto firstEdgeCount = getEdgeCount(seedKey, ACGT[base]);
            if (firstEdgeCount >= opts_.lowCov && firstEdgeCount <= opts_.highCov) {
                KmerInt nextK1mer = SequenceUtils::shiftRight(seedKey, opts_.k-1, ACGT[base]);
                branches.push_back(Branch());
                branches.back().walkers.emplace_back(seq, std::vector<uint16_t>{firstEdgeCount}, nextK1mer);   
            }
        }
        assert(branches.size() >= 2);

        int maxSteps = opts_.getMaxSteps();
        for (int step = 0; step < maxSteps; ++step) {
            for (auto& branch : branches) {
                branch.forward(*this);
            }

            bool merged = false;
            for (size_t i = 0; i < branches.size(); ++i) {
                for (size_t j = 0; j < branches.size(); ++j) {
                    if (i == j) continue;
                    merged = branches[i].isMerged(branches[j]);
                    if (merged) break;
                }
                
                if (merged) break;
            }

            if (merged) break;

            for (auto& branch : branches) {
                branch.resize(opts_);
            }

        }

        if (opts_.isDebug()) {
            // #pragma omp critical
            for (auto& branch : branches) {
                for (const auto& w : branch.walkers) { 
                    std::cout <<"Bubble edge "<< (int)w.state << " " << w.sequence << " " << w.countsAsString() << std::endl;
                }
            }
        }

        // 检查walker之间有配对的吗
        WalkerPair bestWalkerPair = detectBestMatch(branches);

        if (bestWalkerPair.isValid(opts_)) { 
            if (bestWalkerPair.type() == "IsolatedSNP") {
                isolatedSnpCount++;
            } else if (bestWalkerPair.type() == "Indel") {
                indelCount++;
            } else if (bestWalkerPair.type() == "IsolatedIndel") {
                nonIsolatedSnpCount++;
            }
            
            std::string dumpInfo = bestWalkerPair.toString(bubbleCount.fetch_add(1), opts_);
            #pragma omp critical
            {
                snpOut << dumpInfo;
            }
        }
                    
    }
    
    snpOut.close();
    std::cout << "Finished walking bubbles. Total found: " << bubbleCount << std::endl;
    std::cout << "Classification Statistics:" << std::endl;
    std::cout << "  Isolated SNP: " << isolatedSnpCount << " (" << (isolatedSnpCount * 100.0 / bubbleCount) << "%)" << std::endl;
    std::cout << "  Non-Isolated SNP: " << nonIsolatedSnpCount << " (" << (nonIsolatedSnpCount * 100.0 / bubbleCount) << "%)" << std::endl;
    std::cout << "  Indel: " << indelCount << " (" << (indelCount * 100.0 / bubbleCount) << "%)" << std::endl;
}

auto KmerProcessor::detectBestMatch(const std::vector<Branch>& branches) const -> WalkerPair {

    WalkerPair bestWalkerPair;

    for (size_t i = 0; i < branches.size(); ++i) {
        for (size_t j = i+1; j < branches.size(); ++j) {
            for (const auto& w1 : branches[i].walkers) {
                for (const auto& w2 : branches[j].walkers) {
                    WalkerPair walkerPair(opts_.k, w1, w2);
                    walkerPair.findBestMatch(opts_);

                    if (walkerPair.isValid(opts_) && walkerPair.matchScore() > bestWalkerPair.matchScore()) {
                        bestWalkerPair = walkerPair;
                    }
                }
            }
        }
    }
    return bestWalkerPair;
}