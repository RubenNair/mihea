#include "gomea/src/mixed_integer/utils.hpp"

namespace gomea{
namespace mixedinteger{

void prepareFolder(string &folder)
{
    if (!filesystem::exists(folder))
    {
		filesystem::create_directories(folder);
    }
    // TODO Ruben: temporarily removing the creation of these folders, since I don't use them currently.
	// filesystem::create_directories(folder + "/fos");
	// filesystem::create_directories(folder + "/output");
}

void initStatisticsFile(string &folder, bool useBN)
{
    if (!filesystem::exists(folder))
    {
        prepareFolder(folder);
    }

    if(!filesystem::exists(folder + "/statistics.txt"))
    {
        ofstream outFile(folder + "/statistics.txt", ofstream::out);
        if (outFile.fail())
        {
            cerr << "Problems with opening file " << folder + "/elitists.txt!\n";
            exit(0);
        }
        if(useBN)
        {
            outFile << "#Evaluations " << "Time,sec. " << "Fitness " << "populationSize " << "#bins per continuous node " << "boundaries " << "VTR hit" << endl;
        } else 
        {
            outFile << "#Evaluations " << "Time,sec. " << "Fitness " << "populationSize " << "VTR hit" << endl;
        }
        outFile.close();
    }
}

void initElitistFile(string &folder)
{
    if (!filesystem::exists(folder))
    {
        prepareFolder(folder);
    }

    ofstream outFile(folder + "/elitists.txt", ofstream::out);
    if (outFile.fail())
    {
        cerr << "Problems with opening file " << folder + "/elitists.txt!\n";
        exit(0);
    }
    outFile << "#Evaluations " << "Time,sec. " << "Fitness " << "Solution" << endl;
    outFile.close();
}

void initLogFile(string &folder)
{
    if (!filesystem::exists(folder))
    {
        prepareFolder(folder);
    }

    ofstream outFile(folder + "/log.txt", ofstream::out);
    if (outFile.fail())
    {
        cerr << "Problems with opening file " << folder + "/log.txt!\n";
        exit(0);
    }
    // outFile << "#Evaluations " << "Time,sec. " << "Fitness " << endl;
    outFile.close();
}

void initBoundaryStatsFile(string &folder)
{
    if (!filesystem::exists(folder))
    {
        prepareFolder(folder);
    }

    ofstream outFile(folder + "/boundaryStats.txt", ofstream::out);
    if (outFile.fail())
    {
        cerr << "Problems with opening file " << folder + "/boundaryStats.txt!\n";
        exit(0);
    }
    outFile << "#bins per continuous node, " << "Fitness " << endl;
    outFile.close();
}

void writeStatisticsToFile(string &folder, long long numberOfEvaluations, long long time, solution_t<int> *solution, size_t populationSize, bool vtrHit)
{
    ofstream outFile(folder + "/statistics.txt", ofstream::app);
    if (outFile.fail())
    {
        cerr << "Problems with opening file " << folder + "/statistics.txt!\n";
        exit(0);
    }

    outFile << (int)numberOfEvaluations << " " << fixed << setprecision(3) << time/1000.0 << " " <<  setprecision(11) << solution->getObjectiveValue() << " " << populationSize << " " << vtrHit;
    outFile << endl;

    outFile.close();
}

void writeBNStatisticsToFile(string &folder, long long numberOfEvaluations, long long time, solution_BN *solution, size_t populationSize, bool vtrHit)
{
    ofstream outFile(folder + "/statistics.txt", ofstream::app);
    if (outFile.fail())
    {
        cerr << "Problems with opening file " << folder + "/statistics.txt!\n";
        exit(0);
    }

    string binsPerContinuousNode = "(";
    string boundaries = "";
    for(int i = 0; i < solution->getNumberOfNodesToDiscretize(); i++)
    {
        binsPerContinuousNode += to_string(solution->getBoundaries()[i].size()) + ",";
        for(int j = 0; j < solution->getBoundaries()[i].size(); j++)
        {
            boundaries += to_string(solution->getBoundaries()[i][j]) + ",";
        }
        boundaries += ";";
    }
    binsPerContinuousNode += ")";


    outFile << (int)numberOfEvaluations << " " << fixed << setprecision(3) << time/1000.0 << " " <<  setprecision(11) << solution->getObjectiveValue() << " " << populationSize << " " << binsPerContinuousNode << " " << boundaries << " " << vtrHit;
    outFile << endl;

    outFile.close();
}

void writeElitistSolutionToFile(string &folder, long long numberOfEvaluations, long long time, solution_mixed *solution)
{
    ofstream outFile(folder + "/elitists.txt", ofstream::app);
    if (outFile.fail())
    {
        cerr << "Problems with opening file " << folder + "/elitists.txt!\n";
        exit(0);
    }

    outFile << (int)numberOfEvaluations << " " << fixed << setprecision(3) << time/1000.0 << " " <<  setprecision(6) << solution->getObjectiveValue() << " ";
    for (int i = 0; i < solution->getNumberOfVariables(); ++i)
        outFile << +solution->variables[i];
    outFile << " ";
    for (int i = 0; i < solution->getNumberOfCVariables(); ++i)
        outFile << +solution->c_variables[i] << " ";
    outFile << endl;

    outFile.close();
}

void writePopulationToFile(string &folder, vec_t<solution_mixed*> population, string message, bool doLog)
{
    if(doLog)
    {
        ofstream outFile(folder + "/log.txt", ofstream::app);
        if (outFile.fail())
        {
            cerr << "Problems with opening file " << folder + "/log.txt!\n";
            exit(0);
        }
        outFile << message << endl;

        for (size_t i = 0; i < population.size(); ++i)
        {
            outFile << population[i]->getObjectiveValue() << "\t";
            for (int j = 0; j < population[i]->getNumberOfVariables(); ++j)
                outFile << +population[i]->variables[j] << " ";
            outFile << "\t";
            for(int j = 0; j < population[i]->getNumberOfCVariables(); ++j)
                outFile << +population[i]->c_variables[j] << " ";
            outFile << endl;
        }
        outFile << endl;

        outFile.close();
    }
}

void writePopulationBoundaryStatsToFile(string &folder, vec_t<solution_mixed*> population, string message)
{
    // ofstream outFile(folder + "/boundaryStats.txt", ofstream::app);
    // if (outFile.fail())
    // {
    //     cerr << "Problems with opening file " << folder + "/boundaryStats.txt!\n";
    //     exit(0);
    // }

    // outFile << message << endl;
    // for(size_t i = 0; i < population.size(); i++)
    // {
    //     solution_BN *solution = (solution_BN*)population[i];
    //     // Check if cast went well
    //     if(solution == NULL)
    //     {
    //         cerr << "Cast to solution_BN failed in writing boundary stats (utils.cpp)!" << endl;
    //         exit(0);
    //     }

    //     for(size_t j = 0; j < solution->getBoundaries().size(); j++)
    //     {
    //         outFile << solution->getBoundaries()[j].size() << " ";
    //     }
    //     outFile << ", ";
    //     outFile << solution->getObjectiveValue() << endl;
    // }
}

void writeBuildingBlocksToFile(string &folder, vec_t<solution_mixed*> population, string message, int k, bool doLog)
{
    if(doLog)
    {
        ofstream outFile(folder + "/log.txt", ofstream::app);
        streamsize prec = outFile.precision();
        if (outFile.fail())
        {
            cerr << "Problems with opening file " << folder + "/log.txt!\n";
            exit(0);
        }
        outFile << message << endl;

        for (size_t n = 0; n < population.size(); ++n)
        {            
            outFile.precision(1);
            outFile << '[' << setfill('0') << setw(3) << n << "] " << fixed << population[n]->getObjectiveValue() << "\t\t";
            outFile.precision(prec);
            for(int i = 0; i < population[n]->getNumberOfVariables() / k; i++)
            {
                bool is_opt_block = true;
                for(int j = i*k; j < (i+1)*k; j++)
                {
                    if(population[n]->variables[j] == '\000')
                    {
                        is_opt_block = false;
                        break;
                    }
                }
                outFile << (is_opt_block ? 1 : 0) << " ";
                
            }
            outFile << "\t";
            for(int j = 0; j < population[n]->getNumberOfCVariables(); ++j)
                outFile << +population[n]->c_variables[j] << " ";
            outFile << endl;
        }
        outFile << endl;

        outFile.close();
    }
}

bool allBuildingBlocksStillExist(vec_t<solution_mixed*> population, int k)
{
    set<int> blocks;
    for(int i = 0; i < population[0]->getNumberOfVariables(); i++)
    {
        blocks.insert(i);
    }

    for(size_t i = 0; i < population.size(); i++)
    {
        if(blocks.size() == 0)
        {
            return true;
        }
        vec_t<int> blocks_to_remove;
        for(int block : blocks)
        {
            if(population[i]->variables[block] == 1)
            {
                blocks_to_remove.push_back(block);
            }
        }
        for(int block : blocks_to_remove)
        {
            blocks.erase(block);
        }
    }
    if(blocks.size() != 0)
    {
        cout << "Indices that don't have a 1 anymore: ";
        for(int block : blocks)
        {
            cout << block << "(block " << block / k << "), ";
        }
        cout << endl;
        return false;
    }
    return true;
    // set<int> blocks;
    // for(int i = 0; i < population[0]->getNumberOfVariables() / k; i++)
    // {
    //     blocks.insert(i);
    // }

    // for(size_t i = 0; i < population.size(); i++)
    // {
    //     if(blocks.size() == 0)
    //     {
    //         return true;
    //     }
    //     vec_t<int> blocks_to_remove;
    //     for(int block : blocks)
    //     {
    //         bool is_opt_block = true;
    //         for(int j = block*k; j < (block+1)*k; j++)
    //         {
    //             if(population[i]->variables[j] == '\000')
    //             {
    //                 is_opt_block = false;
    //                 break;
    //             }
    //         }
    //         if(is_opt_block)
    //         {
    //             blocks_to_remove.push_back(block);
    //         }
    //     }
    //     for(int block : blocks_to_remove)
    //     {
    //         blocks.erase(block);
    //     }
    // }

    // return blocks.size() == 0;
}

string countBuildingBlocks(vec_t<solution_mixed*> population, int k)
{
    vec_t<int> block_counters(population[0]->getNumberOfVariables() / k, 0);
    vec_t<int> individual_counters(population[0]->getNumberOfVariables(), 0);

    for(size_t i = 0; i < population.size(); i++)
    {
        for(int j = 0; j < population[i]->getNumberOfVariables() / k; j++)
        {
            bool is_opt_block = true;
            for(int l = j*k; l < (j+1)*k; l++)
            {
                if(population[i]->variables[l] == '\000')
                {
                    is_opt_block = false;
                } else
                {
                    individual_counters[l]++;
                }
            }
            if(is_opt_block)
            {
                block_counters[j]++;
            }
        }
    }

    string result = "Block counters: ";
    for(size_t i = 0; i < block_counters.size(); i++)
    {
        result += to_string(block_counters[i]) + " ";
    }
    result += "\nIndividual counters: ";
    for(size_t i = 0; i < individual_counters.size(); i++)
    {
        result += to_string(individual_counters[i]) + " ";
    }
    return result;
}

void printPopulation(vec_t<solution_mixed *> &population)
{
    for (int i = 0; i < population.size(); ++i)
    {
        cout << "Solution " << i << ": ";
        for (int j = 0; j < population[i]->getNumberOfVariables(); ++j)
            cout << +population[i]->variables[j] << " ";
        cout << "\tfitness: " << population[i]->getObjectiveValue() << endl;
    }
}

void writeMatrixToFile(string &folder, double **matrix, int rows, int cols, string message, bool doLog)
{
    if(doLog)
    {
        ofstream outFile(folder + "/log.txt", ofstream::app);
        if (outFile.fail())
        {
            cerr << "Problems with opening file " << folder + "/log.txt!\n";
            exit(0);
        }
        outFile << message << endl;

        for (int i = 0; i < rows; ++i)
        {
            for (int j = 0; j < cols; ++j)
                outFile << matrix[i][j] << " ";
            outFile << endl;
        }
        outFile << endl;

        outFile.close();
    }
}

void writeVectorToFile(string &folder, double *vector, int length, string message, bool doLog)
{
    if(doLog)
    {
        ofstream outFile(folder + "/log.txt", ofstream::app);
        if (outFile.fail())
        {
            cerr << "Problems with opening file " << folder + "/log.txt!\n";
            exit(0);
        }
        outFile << message << endl;

        for (int i = 0; i < length; ++i)
        {
            outFile << vector[i] << " ";
        }
        outFile << endl;

        outFile.close();
    }
}

void writeMessageToLogFile(string &folder, string message, bool doLog)
{
    if(doLog) 
    {
        ofstream outFile(folder + "/log.txt", ofstream::app);
        if (outFile.fail())
        {
            cerr << "Problems with opening file " << folder + "/log.txt!\n";
            exit(0);
        }
        outFile << message << endl;

        outFile.close();
    }
}

void solutionsArchive::checkAlreadyEvaluated(vector<int> &genotype, archiveRecord *result)
{
    result->isFound = false;

    unordered_map<vector<int>, double, hashVector >::iterator it = archive.find(genotype);
    if (it != archive.end())
    {
        result->isFound = true;
        result->value = it->second;
    }
}

void solutionsArchive::insertSolution(vector<int> &genotype, double fitness)
{
    // #if DEBUG
    //  cout << "Inserting solution ";
    //  for (size_t i = 0; i < solution.size(); ++i)
    //      cout << solution[i];
    // #endif
    if (archive.size() >= maxArchiveSize)
        return;
    archive.insert(pair<vector<int>, double> (genotype, fitness));
}

/**
 * Calculates the number of links/archs based in the network
 * @param number_of_nodes The number of nodes (including Continuous and Discrete nodes)
 * @return The total number of links in the network
 */
size_t calculateNumberOfLinks(size_t number_of_nodes) {
    // Determine the number of links/archs
    size_t numberOfArchs = 0;
    for (int i = 0; i < number_of_nodes; ++i) { numberOfArchs += (number_of_nodes - i - 1); }
    return numberOfArchs;
}

/**
 * Find the maximum and minimum values in a 2D array (data).
 * @return A tuple containing the maximum and minimum values in the data.
 */
tuple<vec_t<double>, vec_t<double>> findMaxAndMinValuesInData(vec_t<vec_t<double>> &data) {
    // Initialize the maximum and minimum values
    vec_t<double> maxValues(data[0].size(), numeric_limits<double>::min());
    vec_t<double> minValues(data[0].size(), numeric_limits<double>::max());

    // Loop over the data
    for (const vec_t<double> &dataPoint : data) {
        // Loop over the data point
        for (size_t i = 0; i < dataPoint.size(); ++i) {
            // Update the maximum and minimum values
            maxValues[i] = max(maxValues[i], dataPoint[i]);
            minValues[i] = min(minValues[i], dataPoint[i]);
        }
    }

    return make_tuple(maxValues, minValues);
}

void writeRunCompletedFile(string &folder, const long long numberOfEvaluations, const clock_t startTime, bool doLog) {
    if(!doLog) 
    {
        return;
    }
    // Determine and create run completed file path
    string filepath = folder + "/completed.txt";
    ofstream completedFile;
    completedFile.open(filepath, ofstream::out | ofstream::trunc);

    clock_t currentTime = clock();

    // Write some closing statistics
    completedFile << "         Time / Evaluations / " << endl
        << setw(13) << scientific << setprecision(3) << (double(currentTime - startTime) / CLOCKS_PER_SEC) << " / "
        << setw(13) << fixed << numberOfEvaluations << " / ";

    completedFile.close();
}

///////////////////////
/// Copy data files ///
///////////////////////
/**
 * Copies the data files to the run folder for post-processing
 * @param pathToDataDir The path to the data dir
 * @param targetDir The target run dir
 * @param problemName The problem name
 * @param runIndex the run index
 */
void copyDataFilesToTargetDir(const string& pathToDataDir, const string &targetDir, const string &problemName, int runIndex) {
    // Determine the file paths
    string pathData = determinePathData(pathToDataDir, problemName, runIndex);
    string pathDataInfo = determinePathInfo(pathToDataDir, problemName, runIndex);
    string pathOptimalSolution = determinePathOptimalSolution(pathToDataDir, problemName, runIndex);

    // Copy the files
    copyFile(pathData, targetDir + "/data.txt");
    copyFile(pathDataInfo, targetDir + "/data_info.txt");
    copyFile(pathOptimalSolution, targetDir + "/optimal_solution.txt");
}

/**
 * Determines the path of the data based on the problem name and run index
 * @param pathToDataDir The path to the data folder
 * @param problemName The (base) name of the problem
 * @param runIndex The run index
 * @return The path to the data
 */
string determinePathData(const string& pathToDataDir, const string &problemName, int runIndex) {
    std::ostringstream result;
    result << pathToDataDir << "/" << problemName << "/" << problemName << "_run" << runIndex << ".txt";
    return result.str();
}

/**
 * Determines the path of the data info based on the problem name and run index
 * @param pathToDataDir The path to the data folder
 * @param problemName The (base) name of the problem
 * @param runIndex The run index
 * @return The path to the info of the data
 */
string determinePathInfo(const string& pathToDataDir, const string &problemName, int runIndex) {
    std::ostringstream result;
    result << pathToDataDir << "/" << problemName << "/" << problemName << "_run" << runIndex << "_info.txt";
    return result.str();
}

/**
 * Determines the path of the reference solution based on the problem name and run index
 * @param pathToDataDir The path to the data folder
 * @param problemName The (base) name of the problem
 * @param runIndex The run index
 * @return The path to the reference solution
 */
string determinePathOptimalSolution(const string& pathToDataDir, const string &problemName, int runIndex) {
    std::ostringstream result;
    result << pathToDataDir << "/" << problemName << "/" << problemName << "_run" << runIndex << "_optimalSolution.txt";
    return result.str();
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

/**
 * Writes the parameters of a run to a file
 * @param options The options given by the use
 * @param userInput The parameter input given by the user
 * @param fitnessFunction The fitness function that is retrieved
 */
void writeParametersFile(string &folder, Config *config, const Density *fitnessFunction, bool doLog) {
    if(!doLog)
    {
        return;
    }
    // Determine and create parameters file path
    string filepath = folder + "/parameters.txt";
    ofstream parameters_file;
    parameters_file.open(filepath, ofstream::out | ofstream::trunc);

    // Write parameters to file
    string text;
    text += "Algorithm:" + config->optimizerName + ", ";
    text += "Model:Linkage Tree - No Local Search, ";
    text += "DiscretizationPolicy:ParentPolicy, "; // TODO placeholder, maybe extract this from config
    text += "ProblemIndex:1004, "; // TODO if other problemindices are used, change this
    text += "ProblemName:BNStructureLearning, ";
    text += "FitnessFunction:Density, ";
    text += "PartialEvaluations:" + to_string(false) + ", ";
    text += "PopulationSize:" + to_string(-1) + ", "; // Todo: Change this if we make the population size a variable
    text += "Nodes:" + to_string(fitnessFunction->getNumberOfNodes()) + ", ";
    text += "NumberOfLinks:" + to_string(fitnessFunction->getNumberOfLinks()) + ", ";
    text += "NumberOfDiscretizations:" + to_string(fitnessFunction->getNumberOfNodesToDiscretize()) + ", ";
    text += "NumberOfParameters:" + to_string(config->numberOfdVariables + config->numberOfcVariables) + ", ";
    text += "SampleSize:" + to_string(config->data->getNumberOfDataRows()) + ", ";
    text += "MaximumNumberOfParents:" + to_string(fitnessFunction->getMaximumNumberOfParents()) + ", ";
    text += "MaximumNumberOfDiscretizations:" + to_string(fitnessFunction->getMaxNumberOfDiscretizations()) + ", ";
    text += "MaxEvaluations:" + to_string(config->maximumNumberOfEvaluations) + ", ";
    text += "MaxTime:" + to_string(config->maximumNumberOfSeconds) + ", ";
    text += "vtrUsed:" + to_string(config->fitness->use_vtr) + ", ";
    text += "vtr:" + to_string(config->vtr) + ", ";
    text += "Post-processingRequired:1, ";
    text += "PostRunDiscretizationPolicy:No Discretization, ";
    text += "LinkageModelIndex:-1, ";
    text += "LocalSearchIndex:-1, ";
    text += "FitnessFunctionIndex:-1, ";
    text += "DiscretizationPolicyIndex:0, ";
    text += "RunIndex:-1, ";
    text += "Seed:" + to_string(config->randomSeed) + ", ";
    text += "DataPath:" + fitnessFunction->getFitnessFunctionBaseName() + ", ";
    parameters_file << text << endl;

    // Close file
    parameters_file.close();
}

/**
 * Writes statistics of the current generation, given the sharedInformationPointer (contains the best solution found so far)
 * @param indexOptimizer The index of the optimizer
 * @param solutionToWrite The specific solution to write
 */
void write_multi_start_scheme_statistics(string &folder, solution_mixed *elitist, size_t optimizerIndex,
                                            string &optimizerName, size_t number_of_generations, clock_t startingRunTime,
                                            double avg_elitist_fitness) {
    // Write the solution
    if (elitist) {
        // Prepare variables
        // cast sharedInformationPointer->elitist to solution_BN type pointer
        solution_BN *casted_elitist_sol = (solution_BN *) elitist;

        // Write the statistics of the solution
        write_single_solution_to_multi_start_scheme_statistics(folder, casted_elitist_sol, optimizerName, optimizerIndex, number_of_generations, startingRunTime, avg_elitist_fitness);

        // Write the solution
        write_single_solution_to_multi_start_scheme_solutions(folder, casted_elitist_sol, optimizerName, optimizerIndex, number_of_generations, startingRunTime);
    }

}

/////////////////////////////////////
/// Multi-start scheme statistics ///
/////////////////////////////////////
/**
 * Initializes the generational statistics file
 */
void initialize_multi_start_scheme_statistics_file(string &folder) {
    // Create and open statistics file
    if(!filesystem::exists(folder + "/statistics_MSS.dat"))
    {
        ofstream outFile(folder + "/statistics_MSS.dat", ofstream::out);
        if (outFile.fail())
        {
            cerr << "Problems with opening file " << folder + "/statistics_MSS.dat!\n";
            exit(0);
        }

    // Add header to statistics file
    #ifdef DEBUG_STATISTICS
        outFile
                << "     Algorithm  Gen  OptInd  TotalEvaluations          Time          "
                "Best-objective       Average-objective       Accuracy    Sensitivity AvgHammingDistance       "
                "LogLikelihood       LLDifference    Dist.DistanceKL       ArcRatio      AvgNoDisc        AvgDiscDist        "
                "MaxDiscDist      Ham.Boundaries\n";
    #else
        outFile
                << "     Algorithm  Gen  OptInd  TotalEvaluations          Time          "
                "Best-objective       Average-objective\n";

    #endif
    outFile.close();
    }
}



/**
 * Writes a single solution to the multi start scheme statistics file.
 * @param solutionToWrite The solution to write to the file
 * @param optimizerOfSolution The optimizer that created the solution
 */
void write_single_solution_to_multi_start_scheme_statistics(string &folder, const solution_BN* solutionToWrite,
                                                                              string &optimizerName, size_t optimizerIndex, size_t number_of_generations, 
                                                                              clock_t startingRunTime, double avg_elitist_fitness) {
    ofstream outFile(folder + "/statistics_MSS.dat", ofstream::app);
    // Perform check
    if (outFile.fail())
    {
        cerr << "Problems with opening file " << folder + "/statistics_MSS.dat!\n";
        exit(0);
    }

    // Load variables
        double runTime = double(solutionToWrite->getTimeStamp() - startingRunTime) / CLOCKS_PER_SEC;

    // Write to statistics file
    outFile
            << setw(14) << optimizerName
            << " " << setw(4) << number_of_generations
            << " " << setw(7) << optimizerIndex
            << " " << setw(17) << solutionToWrite->getNumberOfFullEvaluations()
            << " " << setw(13) << scientific << setprecision(3) << runTime
            << " " << setw(23) << scientific << setprecision(16) << solutionToWrite->getObjectiveValue()
            << " " << setw(23) << scientific << setprecision(16) << avg_elitist_fitness
            << endl;
    
    outFile.close();
}


////////////////////////////////////
/// Multi-start scheme solutions ///
////////////////////////////////////
/**
 * Initializes the generational solution file
 */
void initialize_multi_start_scheme_solutions_file(string &folder) {
    // Create and open statistics file
    if(!filesystem::exists(folder + "/statistics_MSS_solutions.dat"))
    {
        ofstream outFile(folder + "/statistics_MSS_solutions.dat", ofstream::out);
        if (outFile.fail())
        {
            cerr << "Problems with opening file " << folder + "/statistics_MSS_solutions.dat!\n";
            exit(0);
        }

        // Add header to statistics file
    outFile
    << "     Algorithm  Gen  OptInd  TotalEvaluations          Time          "
       "Best-objective "
       "Network Discretizations Edges\n";
    
    outFile.close();
    }
}

/**
 * Writes a single solution to the multi start scheme solution file.
 * @param solutionToWrite The solution to write to the file
 * @param optimizerOfSolution The optimizer that created the solution
 */
void write_single_solution_to_multi_start_scheme_solutions(string &folder, const solution_BN* solutionToWrite,
                                                                             string &optimizerName, size_t optimizerIndex, size_t number_of_generations, 
                                                                             clock_t startingRunTime) {
    // Perform check
    ofstream outFile(folder + "/statistics_MSS_solutions.dat", ofstream::app);
    if (outFile.fail())
    {
        cerr << "Problems with opening file " << folder + "/statistics_MSS_solutions.dat!\n";
        exit(0);
    }

    // Load variables
    double runTime = double(solutionToWrite->getTimeStamp() - startingRunTime) / CLOCKS_PER_SEC;

    // Write to statistics file
    outFile
            << setw(14) << optimizerName
            << "/" << setw(4) << number_of_generations
            << "/" << setw(7) << optimizerIndex
            << "/" << setw(17) << solutionToWrite->getNumberOfFullEvaluations()
            << "/" << setw(13) << scientific << setprecision(3) << runTime
            << "/" << setw(23) << scientific << setprecision(16) << solutionToWrite->getObjectiveValue()
            << "/" << convertSolutionNetworkToString(solutionToWrite->getNetworkParameters())
            << "/" << convertInstantiationCountToString(getNumberOfInstantiations(solutionToWrite))
            << "/" << convertBoundariesToString(solutionToWrite->getBoundaries())
            << endl;

    outFile.close();
}


/////////////////////////////
/// Statistical Functions ///
/////////////////////////////
// /**
//  * Calculates the average best solution fitness
//  * @return The average fitness over all elitist solutions
//  */
// double getAverageElitistFitness() {
//     double sum_fitness = 0.0;
//     for (const auto &optimizer : this->optimizers) {
//         sum_fitness += optimizer->getElitistSolution()->getFitness();
//     }

//     return sum_fitness / (double) this->optimizers.size();
// }

///////////////////////
/// Other functions ///
///////////////////////
// /**
//  * Determines the current run time
//  * @return The run time
//  */
// double determineRunTime() {
//     // Determine time that has passed
//     clock_t currentTime = clock();
//     double result = double(currentTime - this->startingRunTime) / CLOCKS_PER_SEC;
//     return result;
// }

/**
 * Converts a network to a string
 * @param network The network to convert
 * @return The network as a string
 */
string convertSolutionNetworkToString(const vector<int> &network) {
    // Initialize string stream
    stringstream ss;
    ss << "[";

    // Add the elements of the network
    for (size_t i = 0; i < network.size(); ++i) {
        ss << network[i];
        // Skip the comma for the last item
        if (i != network.size() - 1) {
            ss << ",";
        }
    }

    ss << "]";
    return ss.str();
}

/**
 * Converts a list of instantiation counts to a string
 * @param instantiations The instantiations to write as a string
 * @return The instantiation count as a string
 */
string convertInstantiationCountToString(const vector<size_t> &instantiations) {
    // Initialize string stream
    stringstream ss;
    ss << "[";

    // Add the elements of the network
    for (size_t i = 0; i < instantiations.size(); ++i) {
        ss << instantiations[i];
        // Skip the comma for the last item
        if (i != instantiations.size() - 1) {
            ss << ",";
        }
    }

    ss << "]";
    return ss.str();

}

/**
 * Converts the boundaries to a string
 * @param boundaries The boundaries
 * @return The boundaries as a string
 */
string convertBoundariesToString(const vector<vector<double>> &boundaries) {
    // Initialize string stream
    stringstream ss;
    ss << "[";

    // Go over each node
    for (size_t i = 0; i < boundaries.size(); ++i) {
        ss << "[";

        // Add the values of the boundaries
        for (size_t j = 0; j < boundaries[i].size(); ++j) {
            ss << scientific << setprecision(16) << boundaries[i][j];

            // Skip the last item
            if (j != boundaries[i].size() - 1) { ss << ","; }
        }

        ss << "]";
        // Skip the last item
        if (i != boundaries.size() - 1) { ss << ";";}
    }

    ss << "]";

    return ss.str();
}

vector<size_t> getNumberOfInstantiations(const solution_BN *solution)
{
    vector<size_t> result;
    int continuousIndex = 0;
    for (size_t nodeIndex = 0; nodeIndex < solution->getNumberOfNodes(); ++nodeIndex) {
        if (solution->getNodeDataTypes()[nodeIndex] == Discrete) {
            result.push_back(solution->getNumberOfDiscretizationsperNode()[nodeIndex]);
        } else {
            result.push_back(solution->getBoundaries()[continuousIndex].size() + 1);
            continuousIndex++;
        }
    }

    return result;
}

}}