#include <iostream>
#include <fstream>
using namespace std;

#include "gomea/src/mixed_integer/simpleGAMBIT.hpp"
#include "gomea/src/fitness/benchmarks-mixed.hpp"
#include "gomea/src/utils/time.hpp"

namespace gomea{
namespace mixedinteger{

simpleGAMBIT::simpleGAMBIT()
{
    // gomeaIMSInstance = new gomeaIMS();
    // iamalgamInstance = new iamalgam();
    return;
}

simpleGAMBIT::simpleGAMBIT(Config *config_): config(config_)
{
    maximumNumberOfGAMBITs   = config->maximumNumberOfGAMBITs;
    problemInstance         = config->fitness;
    basePopulationSize      = config->basePopulationSize;
    problemInstance->maximum_number_of_evaluations = config->maximumNumberOfEvaluations;
    problemInstance->maximum_number_of_seconds = config->maximumNumberOfSeconds;
    IMSsubgenerationFactor = 4;
    problemInstance->initializeRun();
    if( config->fix_seed )
    {
        utils::initializeRandomNumberGenerator(config->randomSeed);
    } else 
    {   
        utils::initializeRandomNumberGenerator();
    }
    // iamalgamInstance = new iamalgam(config);
}

simpleGAMBIT::~simpleGAMBIT()
{}

void simpleGAMBIT::initialize()
{
    cout << "[DEBUGGING] We are here now (simpleGAMBIT::initialize)" << endl;
    utils::initStartTime();
    clock_start_time = clock();
    utils::clearTimers();
    // output = output_statistics_t();

    // if( config->AnalyzeFOS )
    // {
    prepareFolder(config->folder);
    // }

    // gomeaIMSInstance = new gomeaIMS(config);
    // iamalgamInstance = new iamalgam(config);
    initElitistFile(config->folder);
    initStatisticsFile(config->folder, config->useBN);
    initialize_multi_start_scheme_statistics_file(config->folder);
    initialize_multi_start_scheme_solutions_file(config->folder);
    if(config->logDebugInformation)
    {
        initLogFile(config->folder);
    }
    if(config->useBN)
    {
        initBoundaryStatsFile(config->folder);
    }

    // gomeaIMSInstance->initialize();
    // iamalgamInstance->initialize();

    sharedInformationInstance = new sharedInformation(config->maxArchiveSize);

    isInitialized = true;
    gen = 0;
    prevGenElitistFitness = -1;
    prevGenElitistInitialized = false;

}

void simpleGAMBIT::ezilaitini()
{
    for (size_t i = 0; i < GAMBITs.size(); ++i)
    {
        delete GAMBITs[i];
    }
    GAMBITs.clear();
    if( isInitialized )
        delete sharedInformationInstance;
    isInitialized = false;
}

void simpleGAMBIT::run()
{
    initialize();

	try{
		while(!checkTermination())
		{
			if (numberOfGAMBITs < maximumNumberOfGAMBITs)
				initializeNewGAMBIT();

			generationalStepAllGAMBITs();

			numberOfGenerationsGAMBIT++;
		}
	}
	catch( utils::customException const& ){}
	hasTerminated = true;
	if(config->useBN)
    {
        solution_BN *elitist_BN = (solution_BN *) GAMBITs[currentGAMBITIndex]->sharedInformationPointer->elitist;
        writeBNStatisticsToFile(config->folder, GAMBITs[currentGAMBITIndex]->sharedInformationPointer->elitistSolutionHittingTimeEvaluations, GAMBITs[currentGAMBITIndex]->sharedInformationPointer->elitistSolutionHittingTimeMilliseconds, elitist_BN, GAMBITs[currentGAMBITIndex]->populationSize);
        gomea::fitness::BNStructureLearning *BN_fitness = (gomea::fitness::BNStructureLearning *) config->fitness;
        writeParametersFile(config->folder, config, BN_fitness->getDensity(), true);
        copyDataFilesToTargetDir("./data", config->folder, config->problemInstancePath, config->runIndex);
        writeRunCompletedFile(config->folder, problemInstance->full_number_of_evaluations, clock_start_time, true);

    } else
    {
        solution_t<int> *elitist_disc_solution = GAMBITs[currentGAMBITIndex]->sharedInformationPointer->elitist;
        writeStatisticsToFile(config->folder, GAMBITs[currentGAMBITIndex]->sharedInformationPointer->elitistSolutionHittingTimeEvaluations, GAMBITs[currentGAMBITIndex]->sharedInformationPointer->elitistSolutionHittingTimeMilliseconds, elitist_disc_solution, GAMBITs[currentGAMBITIndex]->populationSize);
    }

	ezilaitini();
}

/**
 * Copied from gomea/src/mixed_integer/gomeaIMS.cpp, but seems to be unused. 
 * Instead, the overload (runGeneration(int GAMBITIndex)) method seems to be used.
*/
void simpleGAMBIT::runGeneration()
{
	if( !isInitialized )
		initialize();

	if( checkTermination() )
		return;

	try{
		if (currentGAMBITIndex >= numberOfGAMBITs )
		{
			if( numberOfGAMBITs < maximumNumberOfGAMBITs)
				initializeNewGAMBIT();
			else
				currentGAMBITIndex = numberOfGAMBITs-1;
		}

		if(!GAMBITs[currentGAMBITIndex]->terminated)
			GAMBITs[currentGAMBITIndex]->terminated = checkTerminationGAMBIT(currentGAMBITIndex);
		
		if(!GAMBITs[currentGAMBITIndex]->terminated)
			runGeneration( currentGAMBITIndex );

		if( GAMBITs[currentGAMBITIndex]->numberOfGenerations % IMSsubgenerationFactor == 0 )
			currentGAMBITIndex++;
		else
			currentGAMBITIndex = minimumGAMBITIndex;
	}
	catch( utils::customException const& )
	{
		hasTerminated = true;
		// writeStatistics( currentGOMEAIndex );
		ezilaitini();
	}
}

// IMS version of running integrated algorithm: performs one generation for a GAMBIT instance.
// Most of old run() functionality is moved here.
void simpleGAMBIT::runGeneration(int GAMBITIndex)
{
    // TODO: move currently used run() functionality here (everything inside the while loop)
    Population *currGAMBIT = GAMBITs[GAMBITIndex];
    // writePopulationToFile(config->folder, currGAMBIT->population->solutions, "initial population --------", config->logDebugInformation);
    // writeMessageToLogFile(config->folder, countBuildingBlocks(currGAMBIT->population->solutions, 5), config->logDebugInformation);

    std::streamsize ss = std::cout.precision();
    cout << "[DEBUGGING] GEN: " << gen++ << "\t(probInst) Elitist Fitness: " << fixed << setprecision(7) << currGAMBIT->problemInstance->elitist_objective_value << setprecision(ss) << "\t(sharedInfoPointer) Elitist Fitness: " << currGAMBIT->sharedInformationPointer->elitistFitness  << "\tcurrGAMBIT popsize: " << currGAMBIT->populationSize << "\telitist discrete variables: ";
    for(int i = 0; i < currGAMBIT->problemInstance->number_of_variables; ++i)
    {
        cout << currGAMBIT->sharedInformationPointer->elitist->variables[i];
    }
    if(config->useBN)
    {
        solution_BN *elitist_BN = (solution_BN *) currGAMBIT->sharedInformationPointer->elitist;
        vec_t<vec_t<double>> el_boundaries = elitist_BN->getBoundaries();

        cout << "\t elitist boundaries (#boundaries: (" << fixed << setprecision(15);
        for(size_t i = 0; i < el_boundaries.size()-1; ++i)
        {
            cout << el_boundaries[i].size() << ", ";
        }
        cout << el_boundaries[el_boundaries.size()-1].size() << ")): ";
        for(size_t i = 0; i < el_boundaries.size(); ++i)
        {
            for(size_t j = 0; j < el_boundaries[i].size()-1; ++j)
            {
                cout << el_boundaries[i][j] << ",";
            }
            cout << el_boundaries[i][el_boundaries[i].size()-1] << ";";
        }
    }
    cout << setprecision(ss) << endl;
    // Learn discrete model
    if(config->numberOfVariables > 0)
    {
        currGAMBIT->learnDiscreteModel();

        // Loop over discrete FOS elements
        currGAMBIT->determineFOSOrder();
    }
    currGAMBIT->copyPopulationToOffspring();

    if(config->numberOfVariables > 0)
    {
        for(size_t i = 0; i < currGAMBIT->FOSInstance->size(); ++i)
        {
            // writeMessageToLogFile(config->folder, "####### [GEN " + to_string(gen) + "] FOS element " + to_string(i) + " of " + to_string(currGAMBIT->FOSInstance->size()) +  " -> FOS_index " + to_string(currGAMBIT->FOSInstance->FOSorder[i]) + " #######", config->logDebugInformation);
            
            // Learn continuous model
            //  Also generates new (continuous part of) population, evaluates all solutions and adapts distribution multipliers
            //  Basically does what "makePopulations()" does in iAMaLGaM C code.
            if(config->numberOfcVariables > 0) 
            {
                currGAMBIT->learnContinuousModel();
            }

            // In GAMBIT_K, update changes from continuous part (directly modified in population) to offspring, so discrete part with GOM uses the updated continuous part.
            if(config->dontUseOffspringPopulation)
                currGAMBIT->copyPopulationToOffspring();
            
            currGAMBIT->generateDiscretePopulation(i); // Generate discrete part of new population. Here, Single refers to a single FOS element, not single solution.
            
            // writePopulationToFile(config->folder, currGAMBIT->offspringPopulation->solutions, "OFFSPRING POPULATION After generating DISCRETE population ----------------------------------", config->logDebugInformation);
            // writeMessageToLogFile(config->folder, "\n" + countBuildingBlocks(currGAMBIT->offspringPopulation->solutions, 5) + "\n", config->logDebugInformation);
        }
    } else if(config->numberOfcVariables > 0) 
    {
        currGAMBIT->learnContinuousModel();
    }
    // update new population to previous offspringpopulation
    currGAMBIT->copyOffspringToPopulation();

    // calculate average fitness of population (to match steps from discrete GOMEA code)
    currGAMBIT->calculateAverageFitness();
    
    // check if improvement of elitist is found. If so, write elitist statistics to file.
    if(!prevGenElitistInitialized || currGAMBIT->sharedInformationPointer->elitistFitness < prevGenElitistFitness) // ASSUMING MINIMIZATION
    {
        prevGenElitistFitness = currGAMBIT->sharedInformationPointer->elitistFitness;
        write_multi_start_scheme_statistics(config->folder, sharedInformationInstance->elitist, sharedInformationInstance->optimizerIndex, config->optimizerName, numberOfGenerationsGAMBIT, clock_start_time, getAverageElitistFitness(), config->data->getColumnType());
        prevGenElitistInitialized = true;
    }

    currGAMBIT->numberOfGenerations++;
}

bool simpleGAMBIT::checkTermination()
{
    int i;

    if(checkEvaluationLimitTerminationCriterion() )
		hasTerminated = true;

	if( checkTimeLimitTerminationCriterion() )
		hasTerminated = true;

    if (numberOfGAMBITs == maximumNumberOfGAMBITs)
    {
        for (i = 0; i < maximumNumberOfGAMBITs; i++)
        {
            if (!GAMBITs[i]->terminated)
                return false;
        }

		hasTerminated = true;
    }
    
    return hasTerminated;
}

bool simpleGAMBIT::checkEvaluationLimitTerminationCriterion()
{
    if( !isInitialized )
		return( false );
	if( config->maximumNumberOfEvaluations > 0 && problemInstance->full_number_of_evaluations > config->maximumNumberOfEvaluations )
		hasTerminated = true;
	return hasTerminated; 
}

bool simpleGAMBIT::checkTimeLimitTerminationCriterion()
{
    if( !isInitialized )
        return( false );
    if( config->maximumNumberOfSeconds > 0 && utils::getElapsedTimeSinceStartSeconds() > config->maximumNumberOfSeconds )
        hasTerminated = true;
    return hasTerminated; 
}



void simpleGAMBIT::initializeNewGAMBIT()
{
    #ifdef DEBUG
        cout << "Current number Of GAMBITs is " << numberOfGOMEAs << " | Creating New GAMBIT!\n";
    #endif

    Population *newPopulation = NULL;

    if (numberOfGAMBITs == 0)
    {
        newPopulation = new Population(config, problemInstance, sharedInformationInstance, numberOfGAMBITs, basePopulationSize, numberOfGAMBITs, config->optimizerName);
        // cout << "Initial population:" << endl;
        // printPopulation(newPopulation->population);

    }
    else
        newPopulation = new Population(config, problemInstance, sharedInformationInstance, numberOfGAMBITs, 2 * GAMBITs[numberOfGAMBITs-1]->populationSize, numberOfGAMBITs, config->optimizerName, GAMBITs[0]->FOSInstance ); // Removed population increase, since that complicates things for iAMaLGaM. TODO check if it is easy/useful to include.
    
    GAMBITs.push_back(newPopulation);
    // // NOT IN GAMBIT: update iamalgam population with initial population
    // iamalgamInstance->updatePopulation(numberOfGAMBITs, newPopulation);
    
    numberOfGAMBITs++;
}

void simpleGAMBIT::generationalStepAllGAMBITs()
{
    int GAMBITIndexSmallest, GAMBITIndexBiggest;

    GAMBITIndexBiggest    = numberOfGAMBITs - 1;
    GAMBITIndexSmallest = 0;
    while(GAMBITIndexSmallest <= GAMBITIndexBiggest)
    {
        if (!GAMBITs[GAMBITIndexSmallest]->terminated)
            break;

        GAMBITIndexSmallest++;
    }

    GAMBITGenerationalStepAllGAMBITsRecursiveFold(GAMBITIndexSmallest, GAMBITIndexBiggest);
}

void simpleGAMBIT::GAMBITGenerationalStepAllGAMBITsRecursiveFold(int GAMBITIndexSmallest, int GAMBITIndexBiggest)
{
    int i, GAMBITIndex;

    for(i = 0; i < IMSsubgenerationFactor-1; i++)
    {
        for(GAMBITIndex = GAMBITIndexSmallest; GAMBITIndex <= GAMBITIndexBiggest; GAMBITIndex++)
        {
            if(!GAMBITs[GAMBITIndex]->terminated)
                GAMBITs[GAMBITIndex]->terminated = checkTerminationGAMBIT(GAMBITIndex);

            if((!GAMBITs[GAMBITIndex]->terminated) && (GAMBITIndex >= minimumGAMBITIndex))
			{
				runGeneration( GAMBITIndex );
			}
        }

        for(GAMBITIndex = GAMBITIndexSmallest; GAMBITIndex < GAMBITIndexBiggest; GAMBITIndex++)
            GAMBITGenerationalStepAllGAMBITsRecursiveFold(GAMBITIndexSmallest, GAMBITIndex);
    }
}

void simpleGAMBIT::FindCurrentGAMBIT()
{
    try{
		if (currentGAMBITIndex >= numberOfGAMBITs )
		{
			if( numberOfGAMBITs < maximumNumberOfGAMBITs)
				initializeNewGAMBIT();
			else
				currentGAMBITIndex = numberOfGAMBITs-1;
		}

		if(!GAMBITs[currentGAMBITIndex]->terminated)
			GAMBITs[currentGAMBITIndex]->terminated = checkTerminationGAMBIT(currentGAMBITIndex);
		
		// if(!GAMBITs[currentGAMBITIndex]->terminated)
			// runGeneration( currentGAMBITIndex );

		// if( GAMBITs[currentGAMBITIndex]->numberOfGenerations % iamalgamSubgenerationFactor == 0 )
		// 	currentGOMEAIndex++;
		else
        {
            // The next GAMBIT should be chosen (with larger population size), first make sure it exists.
            if(numberOfGAMBITs < maximumNumberOfGAMBITs)
                initializeNewGAMBIT();
            // Then update currentGAMBITIndex to point to next GAMBIT.
			currentGAMBITIndex = minimumGAMBITIndex;
        }   
	}
	catch( utils::customException const& ) 
    {
        cout << "VTR has been hit, terminating the program. (VTR exception has been caught)" << endl;
        hasTerminated = true;
        // writeStatistics( currentGAMBITIndex );
        ezilaitini();
    }
}

bool simpleGAMBIT::checkTerminationGAMBIT(int GAMBITIndex)
{
	if( checkTermination() )
		return true;

	if( config->maximumNumberOfGenerations > 0 && (int) GAMBITs[GAMBITIndex]->numberOfGenerations >= config->maximumNumberOfGenerations )
	{
        if( GAMBITIndex == minimumGAMBITIndex )
			minimumGAMBITIndex = GAMBITIndex+1;
		return true;
	}

	if( numberOfGAMBITs > 1 )
	{
		for (int i = GAMBITIndex+1; i < numberOfGAMBITs; i++)
		{        
			if (GAMBITs[i]->averageFitness < GAMBITs[GAMBITIndex]->averageFitness) // ASSUMING MINIMIZATION!
			{
				minimumGAMBITIndex = GAMBITIndex+1;
				return true;
			}
		}
	}

	if (!GAMBITs[GAMBITIndex]->allSolutionsAreEqual())
		return false;

	if( GAMBITIndex == minimumGAMBITIndex )
		minimumGAMBITIndex = GAMBITIndex+1;
    return true;
}


/**
 * Calculates the average best solution fitness
 * @return The average fitness over all elitist solutions
 */
double simpleGAMBIT::getAverageElitistFitness() {
    double sum_fitness = 0.0;
    for (size_t i = 0; i < numberOfGAMBITs; ++i) {
        sum_fitness += GAMBITs[i]->optimizerElitistFitness;
    }

    return sum_fitness / (double) numberOfGAMBITs;
}


}}