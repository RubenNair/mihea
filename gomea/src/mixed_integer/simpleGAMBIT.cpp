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
    iamalgamInstance = new iamalgam();
    return;
}

simpleGAMBIT::simpleGAMBIT(Config *config_): config(config_)
{
    maximumNumberOfGAMBITs   = config->maximumNumberOfGAMBITs;
    problemInstance         = config->fitness;
    basePopulationSize      = config->basePopulationSize;
    problemInstance->maximum_number_of_evaluations = config->maximumNumberOfEvaluations;
    problemInstance->maximum_number_of_seconds = config->maximumNumberOfSeconds;
    problemInstance->initializeRun();
    if( config->fix_seed )
    {
        utils::initializeRandomNumberGenerator(config->randomSeed);
    } else 
    {   
        utils::initializeRandomNumberGenerator();
    }
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

    // if( config->AnalyzeFOS )
    // {
    prepareFolder(config->folder);
    // }

    // gomeaIMSInstance = new gomeaIMS(config);
    // iamalgamInstance = new iamalgam(config);
    initElitistFile(config->folder);
    initStatisticsFile(config->folder, config->useBN);
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
    clock_t clock_start_time = clock();
    initialize();
    
    int gen = 0;
    // double prevGenElitistFitness = -1;
    // int numGensNoChange = 0;
    try {
        while(!checkTermination() && (config->maximumNumberOfGenerations <= 0 || gen < config->maximumNumberOfGenerations ))
        {
            writeMessageToLogFile(config->folder, "\n\n\n############## Generation " + to_string(gen) + " ##############", config->logDebugInformation);
            // Moved everthing up to Learn discrete model inside the while loop, I think this is where it should be.
            // if (numberOfGAMBITs < maximumNumberOfGAMBITs)
            // 			initializeNewGAMBIT();
            FindCurrentGAMBIT(); // -> Finds the current GAMBIT to use for this iteration, and initializes a new GAMBIT with larger population size if necessary
            // Cover edge case for very small problems where optimum is found during initialization of population
            if(hasTerminated)
                break;

            // Pass the GAMBITs (populations) to the GOMEA and iamalgam instances
            // _TODO: maybe it's better to figure out in here which population to use this run, then only pass that one? -> Done
            Population *currGAMBIT = GAMBITs[currentGAMBITIndex];           
            writePopulationToFile(config->folder, currGAMBIT->population->solutions, "initial population --------", config->logDebugInformation);
            writeMessageToLogFile(config->folder, countBuildingBlocks(currGAMBIT->population->solutions, 5), config->logDebugInformation);

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
            // gomeaIMSInstance->currGAMBIT = currGAMBIT;
            // gomeaIMSInstance->GAMBITs = GAMBITs;

            // _TODO: iamalgam's number_of_populations should be initialized to be maximumNumberOfGAMBITs -> Done I think
            // iamalgamInstance->currGAMBIT = currGAMBIT;
            // iamalgamInstance->GAMBITs = GAMBITs;

            // Learn discrete model
            // I assume this should be learning the linkage tree as happens in Population::makeOffspring() currently?
            // based on similarity matrix
            if(config->numberOfVariables > 0)
            {
                currGAMBIT->learnDiscreteModel();

                // Loop over discrete FOS elements
                currGAMBIT->determineFOSOrder();
            }
            currGAMBIT->copyPopulationToOffspring();

            // assert((int) currGAMBIT->FOSInstance->size() == 2 * config->numberOfdVariables - 1);
            if(config->numberOfVariables > 0)
            {
                for(size_t i = 0; i < currGAMBIT->FOSInstance->size(); ++i)
                {
                    writeMessageToLogFile(config->folder, "####### [GEN " + to_string(gen) + "] FOS element " + to_string(i) + " of " + to_string(currGAMBIT->FOSInstance->size()) +  " -> FOS_index " + to_string(currGAMBIT->FOSInstance->FOSorder[i]) + " #######", config->logDebugInformation);
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

                    // In GAMBIT_K, update changes from continuous part (directly modified in population) to offspring, so discrete part with GOM uses the updated continuous part.
                    if(config->dontUseOffspringPopulation)
                        currGAMBIT->copyPopulationToOffspring();
                    // Generate new population: replace whole continuous part of population with samples from (updated) continuous model,
                    // use discrete model to find donors for each individual (and then replace only elements in current FOS element).
                    // Idea: generate both parts first, then call some sort of evaluate all function from here. -> Might have to update some parts of iAMaLGaM after that still.
                    currGAMBIT->generateDiscretePopulation(i); // Generate discrete part of new population. Here, Single refers to a single FOS element, not single solution.
                    // cout << "[DEBUGGING] DONE GENERATING SINGLE OFFSPRING (discrete)" << endl;
                    // currGAMBIT->generateNewContinuousPopulation();
                    // cout << "[DEBUGGING] DONE GENERATING SINGLE OFFSPRING (continuous)" << endl;
                    // evaluate all solutions in (new) population -> No longer here; instead in continuous and discrete subparts.
                    // currGAMBIT->evaluateAllSolutionsInPopulation();

                    writePopulationToFile(config->folder, currGAMBIT->offspringPopulation->solutions, "OFFSPRING POPULATION After generating DISCRETE population ----------------------------------", config->logDebugInformation);
                    writeMessageToLogFile(config->folder, "\n" + countBuildingBlocks(currGAMBIT->offspringPopulation->solutions, 5) + "\n", config->logDebugInformation);

                }
            } else if(config->numberOfcVariables > 0) 
            {
                currGAMBIT->learnContinuousModel();
            }
            // update new population to previous offspringpopulation
            currGAMBIT->copyOffspringToPopulation();

            // calculate average fitness of population (to match steps from discrete GOMEA code)
            currGAMBIT->calculateAverageFitness();
            
            currGAMBIT->numberOfGenerations++;

            if(currGAMBIT->iamalgamInstance->number_of_generations > 0 && currGAMBIT->iamalgamInstance->number_of_generations % 100 == 0)
            {
                writeElitistSolutionToFile(config->folder, currGAMBIT->sharedInformationPointer->elitistSolutionHittingTimeEvaluations, currGAMBIT->sharedInformationPointer->elitistSolutionHittingTimeMilliseconds, currGAMBIT->sharedInformationPointer->elitist);
            }

            // Check if discrete population has converged to a solution that is not the optimum. If so, terminate that GAMBIT.
            // NOTE: this assumes that for every problem, the optimum has all discrete variables set to 1.
            // That is not the case for BN problems, for example. Therefore, don't perform this check for BN problems.
            if(!config->useBN)
                currGAMBIT->checkDiscretePopulationConvergedNotOptimal();
            // if(gen >= 5)
            // {
            //     cout << "[DEBUGGING] SETTING WHOLE DISCRETE POPULATION TO 1!" << endl;
            //     // Set whole discrete population to all 1s, see if it can still find the optimum then.
            //     for(int i = 0; i < currGAMBIT->populationSize; ++i)
            //     {
            //         for(int j = 0; j < config->numberOfdVariables; j++)
            //         {
            //             currGAMBIT->population[i]->variables[j] = 1;
            //         }
            //     }
            //     currGAMBIT->evaluateAllSolutionsInPopulation();
            // }

            // if(abs(currGAMBIT->sharedInformationPointer->elitistFitness - prevGenElitistFitness) < 1e-12) 
            // {
            //     numGensNoChange++;
            //     if(numGensNoChange >= 25 && config->numberOfVariables > 0)
            //     {
            //         cout << "[DEBUGGING] No change in elitist fitness for 25 generations, stopping." << endl;
            //         break;
            //     }
            // }
            // else {
            //     numGensNoChange = 0;
            // }
            // prevGenElitistFitness = currGAMBIT->sharedInformationPointer->elitistFitness;

        }
        cout << "[EXIT] termination / max number of generations reached." << endl;
        if(config->useBN)
        {
            solution_BN *elitist_BN = (solution_BN *) GAMBITs[currentGAMBITIndex]->sharedInformationPointer->elitist;
            writeBNStatisticsToFile(config->folder, GAMBITs[currentGAMBITIndex]->sharedInformationPointer->elitistSolutionHittingTimeEvaluations, GAMBITs[currentGAMBITIndex]->sharedInformationPointer->elitistSolutionHittingTimeMilliseconds, elitist_BN, GAMBITs[currentGAMBITIndex]->populationSize);
            // gomea::fitness::BNStructureLearning *BN_fitness = (gomea::fitness::BNStructureLearning *) config->fitness;
            // writeParametersFile(config->folder, config, BN_fitness->getDensity());
            // copyDataFilesToTargetDir(determinePathData("./data", config->problemInstancePath, config->runIndex), config->folder, config->problemInstancePath, config->runIndex);
            // writeRunCompletedFile(config->folder, problemInstance->number_of_evaluations, clock_start_time);

        } else
        {
            solution_t<int> *elitist_disc_solution = GAMBITs[currentGAMBITIndex]->sharedInformationPointer->elitist;
            writeStatisticsToFile(config->folder, GAMBITs[currentGAMBITIndex]->sharedInformationPointer->elitistSolutionHittingTimeEvaluations, GAMBITs[currentGAMBITIndex]->sharedInformationPointer->elitistSolutionHittingTimeMilliseconds, elitist_disc_solution, GAMBITs[currentGAMBITIndex]->populationSize);
        }
    }
    catch( utils::customException const& )
    {
        cout << "VTR has been hit, terminating the program. (VTR exception has been caught)" << endl;
        // Temporary fix: for BN, no VTR is used. When timeout happens, this exception is thrown, but nothing written to file.
        if(config->useBN)
        {
            solution_BN *elitist_BN = (solution_BN *) GAMBITs[currentGAMBITIndex]->sharedInformationPointer->elitist;
            writeBNStatisticsToFile(config->folder, GAMBITs[currentGAMBITIndex]->sharedInformationPointer->elitistSolutionHittingTimeEvaluations, GAMBITs[currentGAMBITIndex]->sharedInformationPointer->elitistSolutionHittingTimeMilliseconds, elitist_BN, GAMBITs[currentGAMBITIndex]->populationSize);
            gomea::fitness::BNStructureLearning *BN_fitness = (gomea::fitness::BNStructureLearning *) config->fitness;
            writeParametersFile(config->folder, config, BN_fitness->getDensity());
            copyDataFilesToTargetDir(determinePathData("./data", config->problemInstancePath, config->runIndex), config->folder, config->problemInstancePath, config->runIndex);
            writeRunCompletedFile(config->folder, problemInstance->number_of_evaluations, clock_start_time);
        }
    }
    
	hasTerminated = true;
	ezilaitini();


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
    {
        newPopulation = new Population(config, problemInstance, sharedInformationInstance, numberOfGAMBITs, basePopulationSize);
        // cout << "Initial population:" << endl;
        // printPopulation(newPopulation->population);

    }
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