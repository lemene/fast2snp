#pragma once

#include <chrono>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

// 128-bit k-mer integer (supports up to k=64 -> 128 bits)
using KmerInt = unsigned __int128;

class SequenceUtils {
public:
    static char getBaseFromBitIndex(int index) {
        static const char bases[] = {'A', 'C', 'G', 'T'};
        return bases[index];
    }

    static uint8_t kmerCharToInt(char c) { 
        switch (c) {
            case 'A': return 0;
            case 'C': return 1;
            case 'G': return 2;
            case 'T': return 3;
            default: return 4;
        }
    }

    static char complement(char c) {
        switch (c) {
            case 'A': return 'T';
            case 'C': return 'G';
            case 'G': return 'C';
            case 'T': return 'A';
            default: return 'N';
        }
    }
    
    static KmerInt maskForBits(int bits) {
        if (bits <= 0) return (KmerInt)0;
        if (bits >= 128) return ~((KmerInt)0);
        return (((KmerInt)1) << bits) - 1;
    }

    static KmerInt kmerToInt(const std::string& s) {
        KmerInt result = 0;
        for (char c : s) {
            result <<= 2;
            switch (c) {
                case 'A': result |= 0b00; break;
                case 'C': result |= 0b01; break;
                case 'G': result |= 0b10; break;
                case 'T': result |= 0b11; break;
            }
        }
        return result;
    }

    static std::string kmerIntToString(KmerInt val, int k) {
        std::string kmer;
        kmer.reserve(k);
        const char bases[] = {'A', 'C', 'G', 'T'};
        for (int i = k - 1; i >= 0; --i) {
            
            int twoBits = static_cast<int>((val >> (2 * i)) & (KmerInt)0x3);
            kmer.push_back(bases[twoBits]);
        }
        return kmer;
    }

    static std::vector<uint8_t> kmerIntToBytes(KmerInt val, int k) { 
        std::vector<uint8_t> kmerBytes;
        kmerBytes.reserve(k);
        for (int i = k - 1; i >= 0; --i) {
            int twoBits = static_cast<int>((val >> (2 * i)) & (KmerInt)0x3);
            kmerBytes.push_back(twoBits);
        }
        return kmerBytes;
    }
    
    static KmerInt reverseInt(KmerInt var, int k) {
        KmerInt mask = maskForBits(k);
        KmerInt result = var & mask;
        
        KmerInt reversed = 0;
        int groupCount = k / 2;
        for (int i = 0; i < groupCount; ++i) {
            KmerInt twoBits = (~result >> (2 * i)) & (KmerInt)0x3;   
            reversed |= twoBits  << (2 * (groupCount - 1 - i)); 
        }
        
        return reversed;
    }

    static std::string reverseString(const std::string& s) {
        std::string result;
        result.reserve(s.size());
        for (auto it = s.rbegin(); it != s.rend(); ++it) {
            result.push_back(complement(*it));
        }
        return result;
    }
    
    static KmerInt shiftRight(KmerInt var, int k, char base) {
        KmerInt nextK1mer = var;
        KmerInt mask = maskForBits(2*k);
        nextK1mer = (nextK1mer << 2) & mask;
        
        switch (base) {
            case 'A': nextK1mer |= 0; break;
            case 'C': nextK1mer |= 1; break;
            case 'G': nextK1mer |= 2; break;
            case 'T': nextK1mer |= 3; break;
        }
        return nextK1mer;
    }
    static std::string getFilename(const std::string& path) {
        size_t pos = path.find_last_of('/');
        return (pos == std::string::npos) ? path : path.substr(pos + 1);
    }
    
    static std::string ensureTrailingSlash(std::string path) {
        if (path.empty() || path.back() != '/') {
            path += '/';
        }
        return path;
    }
    
    static double getPeakRSSInGB() {
        std::ifstream file("/proc/self/status");
        std::string line;
        size_t peakRss = 0;
        while (std::getline(file, line)) {
            if (line.find("VmPeak:") == 0) {
                std::istringstream iss(line);
                std::string key;
                size_t value;
                std::string unit;
                iss >> key >> value >> unit;
                if (unit == "kB") {
                    peakRss = value;
                }
                break;
            }
        }
        return peakRss / (1024.0 * 1024.0);
    }

    // Helper to strip extension for KMC DB
    static std::string stripExtension(const std::string& filename) {
        size_t lastdot = filename.find_last_of(".");
        if (lastdot == std::string::npos) return filename;
        std::string ext = filename.substr(lastdot);
        if (ext == ".kmc_pre" || ext == ".kmc_suf") {
            return filename.substr(0, lastdot);
        }
        return filename;
    }
};

class AutoTimer {
public:
    AutoTimer(const char* name) : name_(name) {
        start_time_ = std::chrono::high_resolution_clock::now();
    }

    ~AutoTimer() {
        std::cout << "Peak Memory: " << SequenceUtils::getPeakRSSInGB() << " GB" << std::endl;

        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time_);
        std::cout << name_ << ": " << duration.count() << " ms" << std::endl;
    }

protected:
    const char* name_;
    std::chrono::time_point<std::chrono::high_resolution_clock> start_time_;

};