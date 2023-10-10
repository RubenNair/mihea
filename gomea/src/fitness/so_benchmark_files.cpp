//
// Created by Damy Ha on 17-may-23.
//

#include "gomea/src/fitness/so_benchmarks.h"

//////////////////
/// File Paths ///
//////////////////
/**
 * Determines the path of the data based on the problem name and run index
 * @param problemName The (base) name of the problem
 * @param runIndex The run index
 * @return The path to the data
 */
string determinePathData(const string &problemName, int runIndex) {
    std::ostringstream result;
    result << "./data/" << problemName << "/" << problemName << "_run" << runIndex << ".txt";
    return result.str();
}

/**
 * Determines the path of the data info based on the problem name and run index
 * @param problemName The (base) name of the problem
 * @param runIndex The run index
 * @return The path to the info of the data
 */
string determinePathInfo(const string &problemName, int runIndex) {
    std::ostringstream result;
    result << "./data/" << problemName << "/" << problemName << "_run" << runIndex << "_info.txt";
    return result.str();
}

/**
 * Determines the path of the data info based on the problem name and run index
 * @param problemName The (base) name of the problem
 * @param runIndex The run index
 * @return The path to the info of the data
 */
string determinePathOptimalSolution(const string &problemName, int runIndex) {
    std::ostringstream result;
    result << "./data/" << problemName << "/" << problemName << "_run" << runIndex << "_optimalSolution.txt";
    return result.str();
}

///////////////////////
/// Copy data files ///
///////////////////////
void copyDataFilesToTargetDir(const string &targetDir, const string &problemName, int runIndex) {
    // Determine the file paths
    string pathData = determinePathData(problemName, runIndex);
    string pathDataInfo = determinePathInfo(problemName, runIndex);
    string pathOptimalSolution = determinePathOptimalSolution(problemName, runIndex);

    // Copy the files
    copyFile(pathData, targetDir + "data.txt");
    copyFile(pathDataInfo, targetDir + "data_info.txt");
    copyFile(pathOptimalSolution, targetDir + "optimal_solution.txt");
}

/**
 * Copies files from inputPath to outputPat
 * @param inputPath The input file
 * @param outputPath The output file
 * @return bool if the operation has succeeded
 */
bool copyFile(const string &inputPath, const string &outputPath) {
    ifstream inputFile(inputPath, ios::binary);
    ofstream outputFile(outputPath, ios::binary | std::ios::trunc); // Overwrite existing file

    // Check if the file exists
    if (!inputFile) {
        cout << "Input file does not exist." << endl;
        return false;
    }

    if (!outputFile) {
        cout << "Unable to create output file." << endl;
        return false;
    }

    // Copy the content
    outputFile << inputFile.rdbuf();

    // Close the files
    inputFile.close();
    outputFile.close();

    return true;
}

//////////////////
/// Read files ///
//////////////////
/**
 * Initializes data if necessary
 * @param initializeData Boolean if the data needs to be initialized
 * @param baseName The basename of the file where the data needs to be retrieved from.
 * @param runIndex The index of the run
 * @return The data or a null pointer
 */
shared_ptr<DataStructure<double>> initializeDataFromPath(bool initializeData, const string& baseName, int runIndex) {
    // Load the data if the data needs to be initialized
    shared_ptr<DataStructure<double>> result;
    if (initializeData) {
        result = make_shared<DataStructure<double>>(determinePathData(baseName, runIndex), determinePathInfo(baseName, runIndex));
    } else {
        result = nullptr;
    }

    return result;
}

/**
 * Retrieves the optimal solution and the boundaries.
 * Assumes that a file exists that contains at least the optimal network.
 * If the boundaries are unknown, these are expected to be set to "None".
 * @param pathBestSolution The path to the best solution file.
 * @param stringOriginalNetwork The original network used to create the data, retrieved from the file.
 * @param originalBoundaries The boundaries where the data is sampled from, retrieved from the file.
 */
void retrieveOptimalSolution(const string& pathBestSolution, string& stringOriginalNetwork, vector<vector<double>> &originalBoundaries) {

    // Check if path to data exists
    if (!path_exists(pathBestSolution)) {
        // Generate error message
        std::string error_message = std::string("Path to:") + pathBestSolution +
                std::string(" does not exist. Please run generate_data.py to generate data "
                            "(and put the file in the data folder).\n");
        throw std::runtime_error(error_message);
    }

    // Open the best solution file
    ifstream file(pathBestSolution, ios_base::binary);
    if (file.is_open()) {
        // Read the first line, which contains the optimal network
        getline(file, stringOriginalNetwork);

        // Read the second line, which contains the boundaries
        string line;
        getline(file, line);
        if (line != "None") {
            istringstream ss(line);
            string token;
            while (getline(ss, token, ';')) {
                // Initialize inner result
                vector<double> inner;

                // Check if the boundary has values
                token.erase(std::remove_if(token.begin(), token.end(), [](char c) { return isspace(c); }), token.end()); // Remove spaces
                if (!token.empty()) {
                    istringstream inner_ss(token);
                    string inner_token;
                    while (getline(inner_ss, inner_token, ',')) {
                        if (!inner_token.empty()) {
                            inner.push_back(std::stod(inner_token));
                        }
                    }
                }

                // Add the boundaries
                originalBoundaries.push_back(inner);
            }
        }

        // Remove the last (empty) element
        originalBoundaries.pop_back();

        // Close the file
        file.close();
    } else {
        cout << "Unable to open the file." << endl;
    }
}