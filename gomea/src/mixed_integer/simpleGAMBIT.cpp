#include <iostream>
#include <fstream>
using namespace std;

#include "gomea/src/mixed_integer/simpleGAMBIT.hpp"

namespace gomea{
namespace mixedinteger{

simpleGAMBIT::simpleGAMBIT()
{
    // gomeaIMSInstance = new gomeaIMS();
    iamalgamInstance = new iamalgam();
    return;
}

simpleGAMBIT::simpleGAMBIT(Config *config_): config(config_)
{
    // gomeaIMSInstance = new gomeaIMS(config);
    iamalgamInstance = new iamalgam(config);
}

simpleGAMBIT::~simpleGAMBIT()
{}

void simpleGAMBIT::initialize()
{
    cout << "[DEBUGGING] We are here now (simpleGAMBIT::initialize)" << endl;
    utils::initStartTime();
    utils::clearTimers();
    output = output_statistics_t();
    maximumNumberOfGAMBITs   = config->maximumNumberOfGAMBITs;
    problemInstance         = config->fitness;
    basePopulationSize      = config->basePopulationSize;
    problemInstance->initializeRun();

    if( config->AnalyzeFOS )
    {
        prepareFolder(config->folder);
    }

    // gomeaIMSInstance = new gomeaIMS(config);
    // iamalgamInstance = new iamalgam(config);
    initElitistFile(config->folder);
    initStatisticsFile(config->folder);
    if(config->logDebugInformation)
    {
        initLogFile(config->folder);
    }

    // gomeaIMSInstance->initialize();
    // iamalgamInstance->initialize();

    sharedInformationInstance = new sharedInformation(config->maxArchiveSize);

    isInitialized = true;

}

void simpleGAMBIT::ezilaitini()
{
    for (size_t i = 0; i < GAMBITs.size(); ++i)
        delete GAMBITs[i];
    GAMBITs.clear();
    if( isInitialized )
        delete sharedInformationInstance;
    isInitialized = false;
}

// This method performs the integrated algorithm from the GAMBIT (2014) paper.
void simpleGAMBIT::run()
{

    initialize();
    
    int gen = 0;
    double prevGenElitistFitness = -1;
    int numGensNoChange = 0;
    while(!checkTermination() && (config->maximumNumberOfGenerations <= 0 || gen < config->maximumNumberOfGenerations ))
    {
        writeMessageToLogFile(config->folder, "\n\n\n############## Generation " + to_string(gen) + " ##############");
        // Moved everthing up to Learn discrete model inside the while loop, I think this is where it should be.
        // if (numberOfGAMBITs < maximumNumberOfGAMBITs)
        // 			initializeNewGAMBIT();
        FindCurrentGAMBIT(); // -> Finds the current GAMBIT to use for this iteration, and initializes a new GAMBIT with larger population size if necessary
        
        // Pass the GAMBITs (populations) to the GOMEA and iamalgam instances
        // _TODO: maybe it's better to figure out in here which population to use this run, then only pass that one? -> Done
        Population *currGAMBIT = GAMBITs[currentGAMBITIndex];
        writePopulationToFile(config->folder, currGAMBIT->population, "initial population --------");
        // gomeaIMSInstance->currGAMBIT = currGAMBIT;
        // gomeaIMSInstance->GAMBITs = GAMBITs;

        // _TODO: iamalgam's number_of_populations should be initialized to be maximumNumberOfGAMBITs -> Done I think
        // iamalgamInstance->currGAMBIT = currGAMBIT;
        // iamalgamInstance->GAMBITs = GAMBITs;

        // Learn discrete model
        // I assume this should be learning the linkage tree as happens in Population::makeOffspring() currently?
        // based on similarity matrix
        currGAMBIT->learnDiscreteModel();

        // Loop over discrete FOS elements
        currGAMBIT->determineFOSOrder();

        // assert((int) currGAMBIT->FOSInstance->size() == 2 * config->numberOfdVariables - 1);
        for(size_t i = 0; i < currGAMBIT->FOSInstance->size(); ++i)
        {
            writeMessageToLogFile(config->folder, "####### [GEN " + to_string(gen) + "] FOS element " + to_string(i) + " of " + to_string(currGAMBIT->FOSInstance->size()) + " #######");
            // cout << "[DEBUGGING] \tFOS element: " << i << endl;
            // // Perform truncation selection -> only learn continuous model based on top tau percent of solutions.
            // // Maybe easier to just do that in iamalgam.cpp instead of creating selection here?
            // double **selection = NULL; // TODO 

            // Learn continuous model
            //  Also generates new (continuous part of) population, evaluates all solutions and adapts distribution multipliers
            //  Basically does what "makePopulations()" does in iAMaLGaM C code.
            if(config->numberOfcVariables > 0) 
            {
                currGAMBIT->learnContinuousModel();
            }

            // Generate new population: replace whole continuous part of population with samples from (updated) continuous model,
            // use discrete model to find donors for each individual (and then replace only elements in current FOS element).
            // Idea: generate both parts first, then call some sort of evaluate all function from here. -> Might have to update some parts of iAMaLGaM after that still.
            currGAMBIT->generateDiscretePopulation(i); // Generate discrete part of new population. Here, Single refers to a single FOS element, not single solution.
            // cout << "[DEBUGGING] DONE GENERATING SINGLE OFFSPRING (discrete)" << endl;
            // currGAMBIT->generateNewContinuousPopulation();
            // cout << "[DEBUGGING] DONE GENERATING SINGLE OFFSPRING (continuous)" << endl;
            // evaluate all solutions in (new) population -> No longer here; instead in continuous and discrete subparts.
            // currGAMBIT->evaluateAllSolutionsInPopulation();

            // update new population to previous offspringpopulation
            currGAMBIT->copyOffspringToPopulation();

        }
        cout << "[DEBUGGING] GEN: " << gen++ << "\t(probInst) Elitist Fitness: " << currGAMBIT->problemInstance->elitist_objective_value << "\t(sharedInfoPointer) Elitist Fitness: " << currGAMBIT->sharedInformationPointer->elitistFitness  << "\tcurrGAMBIT popsize: " << currGAMBIT->populationSize << endl;
        if(abs(currGAMBIT->sharedInformationPointer->elitistFitness - prevGenElitistFitness) < 1e-12) 
        {
            numGensNoChange++;
            if(numGensNoChange >= 25)
            {
                cout << "[DEBUGGING] No change in elitist fitness for 25 generations, stopping." << endl;
                break;
            }
        }
        else {
            numGensNoChange = 0;
        }
        prevGenElitistFitness = currGAMBIT->sharedInformationPointer->elitistFitness;

    }


//     for (int i = 0; i < n; i++)
//     {
//         CreateRandomSolution();
//         EvaluateFitness(i);
//     }

//     while (!terminationCriterion)
//     {
//         LearnDiscreteModel();
//         for (int i = 0; i < 2 * l_d; i++)
//         {
//             selection = TruncationSelection(population[i], tao);
//             LearnContinuousModel(selection);
//             for(int j = 0; j < n; j++)
//             {
//                 newcPopulation[i] = GenerateContinuousPart(cPopulation[i]); // Not like this, keep populations in Population objects (GAMBITs)
//                 newdPopulation[i] = GenerateDiscretePart(j, newdPopulation[i], Population); // Not like this, keep populations in Population objects (GAMBITs)
//             }
//         }
//         cPopulation = newcPopulation; // Not like this, keep populations in Population objects (GAMBITs)
//         dPopulation = newdPopulation; // Not like this, keep populations in Population objects (GAMBITs)
//     }

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
	if( config->maximumNumberOfEvaluations > 0 && problemInstance->number_of_evaluations > config->maximumNumberOfEvaluations )
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
        newPopulation = new Population(config, problemInstance, sharedInformationInstance, numberOfGAMBITs, basePopulationSize);
    else
        newPopulation = new Population(config, problemInstance, sharedInformationInstance, numberOfGAMBITs, 2 * GAMBITs[numberOfGAMBITs-1]->populationSize, GAMBITs[0]->FOSInstance ); // Removed population increase, since that complicates things for iAMaLGaM. TODO check if it is easy/useful to include.
    
    GAMBITs.push_back(newPopulation);
    // // NOT IN GAMBIT: update iamalgam population with initial population
    // iamalgamInstance->updatePopulation(numberOfGAMBITs, newPopulation);
    
    numberOfGAMBITs++;
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
			if (GAMBITs[i]->averageFitness > GAMBITs[GAMBITIndex]->averageFitness) 
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



}}