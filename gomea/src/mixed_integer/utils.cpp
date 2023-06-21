#include "gomea/src/mixed_integer/utils.hpp"

namespace gomea{
namespace mixedinteger{

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
    if (!filesystem::exists(folder))
    {
        prepareFolder(folder);
    }

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

void writeStatisticsToFile(string &folder, long long numberOfEvaluations, long long time, solution_t<char> *solution)
{
    ofstream outFile(folder + "/statistics.txt", ofstream::app);
    if (outFile.fail())
    {
        cerr << "Problems with opening file " << folder + "/statistics.txt!\n";
        exit(0);
    }

    outFile << (int)numberOfEvaluations << " " << fixed << setprecision(3) << time/1000.0 << " " <<  setprecision(11) << solution->getObjectiveValue();
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

void writePopulationToFile(string &folder, vec_t<solution_mixed*> population, string message)
{
    if(filesystem::exists(folder + "/log.txt"))
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

void writeMatrixToFile(string &folder, double **matrix, int rows, int cols, string message)
{
    if(filesystem::exists(folder + "/log.txt"))
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

void writeVectorToFile(string &folder, double *vector, int length, string message)
{
    if(filesystem::exists(folder + "/log.txt"))
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

void writeMessageToLogFile(string &folder, string message)
{
    if(filesystem::exists(folder + "/log.txt")) 
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