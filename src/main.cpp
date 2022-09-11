#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <filesystem>
// 3rd-Party Libraries
#include <boost/program_options.hpp>
// Original
#include "gmm.hpp"

// Define Namespace
namespace fs = std::filesystem;
namespace po = boost::program_options;

// Function Prototype
void Collect_Paths(const std::string root, const std::string sub, std::vector<std::string> &paths);
std::vector<std::string> Get_Paths(const std::string root);
std::vector<std::vector<double>> Get_Data(const std::vector<std::string> paths, const long long D);


// ----------------------
// 0. Argument Function
// ----------------------
po::options_description parse_arguments(){

    po::options_description args("Options", 200, 30);

    args.add_options()

        // (1) Define for General Parameter
        ("help", "produce help message")
        ("dataset", po::value<std::string>()->default_value("toy"), "dataset name")
        ("D", po::value<long long>()->default_value(2), "number of dimensions")
        ("K", po::value<long long>()->default_value(4), "number of normal distribution")
        ("verbose", po::value<bool>()->default_value(true), "verbose")

        // (2) Define for Training
        ("train_dir", po::value<std::string>()->default_value("train"), "training directory : ./datasets/<dataset>/<train_dir>/<all data>")

        // (3) Define for Test
        ("test_dir", po::value<std::string>()->default_value("test"), "test directory : ./datasets/<dataset>/<test_dir>/<all data>")

    ;

    return args;

}


// ------------------
// 1. Main Function
// ------------------
int main(int argc, const char *argv[]){

    // (1) Extract Arguments
    po::options_description args = parse_arguments();
    po::variables_map vm{};
    po::store(po::parse_command_line(argc, argv, args), vm);
    po::notify(vm);
    if (vm.empty() || vm.count("help")){
        std::cout << args << std::endl;
        return 1;
    }

    // (2.1) Get Training Data
    std::string train_dir;
    std::vector<std::string> train_paths;
    std::vector<std::vector<double>> train_data;
    /*****************************************************/
    train_dir = "datasets/" + vm["dataset"].as<std::string>() + "/" + vm["train_dir"].as<std::string>();
    train_paths = Get_Paths(train_dir);
    train_data = Get_Data(train_paths, vm["D"].as<long long>());

    // (2.2) Training for GMM
    GMM gmm(vm["verbose"].as<bool>());
    gmm.train(train_data, vm["D"].as<long long>(), vm["K"].as<long long>());

    // (3.1) Get Test Data
    std::string test_dir;
    std::vector<std::string> test_paths;
    std::vector<std::vector<double>> test_data;
    /*****************************************************/
    test_dir = "datasets/" + vm["dataset"].as<std::string>() + "/" + vm["test_dir"].as<std::string>();
    test_paths = Get_Paths(test_dir);
    test_data = Get_Data(test_paths, vm["D"].as<long long>());

    // (3.2) Test for GMM

    return 0;
}


// ------------------------------
// 2. Collecting Paths Function
// ------------------------------
void Collect_Paths(const std::string root, const std::string sub, std::vector<std::string> &paths){
    
    fs::path ROOT(root);
    
    for (auto &p : fs::directory_iterator(ROOT)){
        if (!fs::is_directory(p)){
            std::stringstream rpath;
            rpath << p.path().string();
            paths.push_back(rpath.str());
        }
        else{
            std::stringstream subsub;
            subsub << p.path().filename().string();
            Collect_Paths(root + '/' + subsub.str(), sub + subsub.str() + '/', paths);
        }
    }
    
    return;
}


// ---------------------------
// 3. Getting Paths Function
// ---------------------------
std::vector<std::string> Get_Paths(const std::string root){

    std::vector<std::string> paths;

    Collect_Paths(root, "", paths);
    std::sort(paths.begin(), paths.end());

    return paths;

}


// ---------------------------
// 4. Getting Data Function
// ---------------------------
std::vector<std::vector<double>> Get_Data(const std::vector<std::string> paths, const long long D){

    long long i;
    double element;
    std::ifstream ifs;
    std::vector<double> data_one;
    std::vector<std::vector<double>> data;

    for (std::string path : paths){
        
        ifs.open(path);

        data_one = std::vector<double>(D);
        for (i = 0; i < D; i++){
            ifs >> element;
            data_one[i] = element;
        }
        data.push_back(data_one);
        
        ifs.close();
    
    }

    return data;
}