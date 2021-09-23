#include <iostream>
#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
#include <sys/stat.h>
#include <sastools/include/vector3.h>
#include <sastools/include/utils.h>
#include <Engine.h>

namespace po = boost::program_options;
namespace {
    const size_t ERROR_IN_COMMAND_LINE = 1;
    const size_t SUCCESS = 0;
    const size_t ERROR_UNHANDLED_EXCEPTION = 2;
}


bool fileExists(const std::string& file) {
    struct stat buf;
    return (stat(file.c_str(), &buf) == 0);
}

int main(int argc, char** argv) {

    std::string ref, tar, fileList;

    std::string descText =
            "\n  To align a set of PDB files : aligner --list files.inp\n";
    descText += "\n     To align against a model : aligner -r reference.pdb -t target.pdb ";
    descText += "\n                                                   ";
    descText += "\n     files.inp should be a list of pdb files to align";
    descText += "\n";
    descText += "\n   Program uses coherent point drift (CPD) algorithm for point set registration";
    descText += "\n   C++ library by Gadomski, P.J. (December 2016).";
    descText += "\n   CPD is default and can be switched to SUPCOMB like method with -c flag";
    descText += "\n   Mirror images are always checked by default, use --noMirror to switch off";
    descText += "\n   ";


    po::options_description desc(descText);

    desc.add_options()
            ("help,h", "Print help messages")
            ("reference,r", po::value<std::string>(&ref), "reference PDB file")
            ("target,t", po::value<std::string>(&tar), "target")
            ;

    try {

        po::variables_map vm;
        po::store(parse_command_line(argc, argv, desc), vm);
        // checking options list
        po::notify(vm);

        if (vm.count("help")) {
            std::cout << desc << "\n";
            return SUCCESS;
        }


        if (!fileExists(ref) ){ // check file exists
            SASTOOLS_UTILS_H::logger("CANNOT READ FILE", ref);
            throw std::invalid_argument("** ERROR => file not adequate " + ref);
        }

        if (!fileExists(tar) ){ // check file exists
            SASTOOLS_UTILS_H::logger("CANNOT READ FILE", tar);
            throw std::invalid_argument("** ERROR => file not adequate " + tar);
        }

        Engine motor = Engine(ref, tar);
        motor.search();
return 1;
    } catch (boost::program_options::required_option& e){
        std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
        return ERROR_IN_COMMAND_LINE;
    } catch (boost::program_options::error& e){
        std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
        return ERROR_IN_COMMAND_LINE;
    } catch (const std::invalid_argument& e){
        std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
        std::cerr<<"Type "<<typeid(e).name()<<std::endl;
    }

    return 0;
}