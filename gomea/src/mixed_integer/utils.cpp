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

void initStatisticsFile(string &folder)
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
        outFile << "#Evaluations " << "Time,sec. " << "Fitness " << "populationSize " << "VTR hit" << endl;
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


}}