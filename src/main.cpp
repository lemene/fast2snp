#include <omp.h>

#include "options.hpp"
#include "fast2snp.hpp"

int main(int argc, const char* argv[]) {
    Options opts;
    
    if (!opts.parse(argc, argv)) {
        //opts.printUsage(std::cout);
        return 1;
    }

    omp_set_num_threads(opts.threads);
    AutoTimer timer("Total Execution Time");

    try {
        KmerProcessor processor(opts);
        processor.loadKmersFromTxt();
        std::string filename = SequenceUtils::getFilename(opts.inputFile);
        filename = SequenceUtils::stripExtension(filename);
        processor.findAndWalkBubbles();
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
