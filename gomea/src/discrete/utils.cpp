#include "gomea/src/discrete/utils.hpp"

namespace gomea{
namespace discrete{

void prepareFolder(string &folder)
{
    if (!filesystem::exists(folder))
    {
		filesystem::create_directories(folder);
    }
	filesystem::create_directories(folder + "/fos");
	filesystem::create_directories(folder + "/output");
}

void initStatisticsFile(string &folder)
{
    ofstream outFile(folder + "/statistics.txt", ofstream::out);
    if (outFile.fail())
    {
        cerr << "Problems with opening file " << folder + "/elitists.txt!\n";
        exit(0);
    }
    outFile << "#Evaluations " << "Time,sec. " << "Fitness " << endl;
    outFile.close();
}

void initElitistFile(string &folder)
{
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

void writeStatisticsToFile(string &folder, long long numberOfEvaluations, long long time, solution_t<char> *solution)
{
    ofstream outFile(folder + "/statistics.txt", ofstream::app);
    if (outFile.fail())
    {
        cerr << "Problems with opening file " << folder + "/statistics.txt!\n";
        exit(0);
    }

    outFile << (int)numberOfEvaluations << " " << fixed << setprecision(3) << time/1000.0 << " " <<  setprecision(6) << solution->getObjectiveValue();
    outFile << endl;

    outFile.close();
}

void writeElitistSolutionToFile(string &folder, long long numberOfEvaluations, long long time, solution_t<char> *solution)
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
    outFile << endl;

    outFile.close();
}

void writePopulationToFile(string &folder, vec_t<solution_t<char>*> population, string message, bool doLog)
{
    if(filesystem::exists(folder + "/log.txt") && doLog)
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
            outFile << "\t" << endl;
        }
        outFile << endl;

        outFile.close();
    }
}

void writeBuildingBlocksToFile(string &folder, vec_t<solution_t<char>*> population, string message, int k, bool doLog)
{
    if(filesystem::exists(folder + "/log.txt") && doLog)
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
            outFile << "\t" << endl;
        }
        outFile << endl;

        outFile.close();
    }
}

string countBuildingBlocks(vec_t<solution_t<char>*> population, int k)
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

void printPopulation(vec_t<solution_t<char> *> &population)
{
    for (int i = 0; i < population.size(); ++i)
    {
        cout << "Solution " << i << ": ";
        for (int j = 0; j < population[i]->getNumberOfVariables(); ++j)
            cout << +population[i]->variables[j] << " ";
        cout << "\tfitness: " << population[i]->getObjectiveValue() << endl;
    }
}

void writeMessageToLogFile(string &folder, string message, bool doLog)
{
    if(filesystem::exists(folder + "/log.txt") && doLog) 
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

void solutionsArchive::checkAlreadyEvaluated(vector<char> &genotype, archiveRecord *result)
{
    result->isFound = false;

    unordered_map<vector<char>, double, hashVector >::iterator it = archive.find(genotype);
    if (it != archive.end())
    {
        result->isFound = true;
        result->value = it->second;
    }
}

void solutionsArchive::insertSolution(vector<char> &genotype, double fitness)
{
    // #if DEBUG
    //  cout << "Inserting solution ";
    //  for (size_t i = 0; i < solution.size(); ++i)
    //      cout << solution[i];
    // #endif
    if (archive.size() >= maxArchiveSize)
        return;
    archive.insert(pair<vector<char>, double> (genotype, fitness));
}


}}