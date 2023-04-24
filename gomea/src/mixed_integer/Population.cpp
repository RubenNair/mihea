#include "gomea/src/mixed_integer/Population.hpp"

namespace gomea{
namespace mixedinteger{

Population::Population(Config *config_, fitness_t *problemInstance_, sharedInformation *sharedInformationPointer_, size_t GOMEAIndex_, size_t dPopulationSize_, linkage_model_pt FOSInstance_ ): 
        config(config_), 
        problemInstance(problemInstance_),
        sharedInformationPointer(sharedInformationPointer_),
        GOMEAIndex(GOMEAIndex_), 
        dPopulationSize(dPopulationSize_)
{
        terminated = false;
        numberOfGenerations = 0;
        averageFitness = 0.0;
        
        dPopulation.resize(dPopulationSize);
        offspringdPopulation.resize(dPopulationSize);
        noImprovementStretches.resize(dPopulationSize);
        
        vec_t<int> allGenes(problemInstance->number_of_variables);
        iota(allGenes.begin(), allGenes.end(), 0);

        for (size_t i = 0; i < dPopulationSize; ++i)
        {
            noImprovementStretches[i] = 0;

            dPopulation[i] = new solution_t<char>(problemInstance->number_of_variables, config->alphabetSize);
            dPopulation[i]->randomInit(&gomea::utils::rng);
            problemInstance->evaluate(dPopulation[i]);
            
            offspringdPopulation[i] = new solution_t<char>(problemInstance->number_of_variables, config->alphabetSize);
            *offspringdPopulation[i] = *dPopulation[i];
        }
			
		if( config->linkage_config != NULL )
		{
			FOSInstance = linkage_model_t::createFOSInstance( *config->linkage_config, problemInstance->number_of_variables );
            FOSInstance->initializeDependentSubfunctions( problemInstance->subfunction_dependency_map );
		}
		else if( FOSInstance_ == NULL )
		{
			FOSInstance = linkage_model_t::createLinkageTreeFOSInstance(config->FOSIndex, problemInstance->number_of_variables, config->linkage_config->lt_similarity_measure, config->linkage_config->lt_maximum_set_size);
		}
		else FOSInstance = FOSInstance_;
        
        #ifdef DEBUG
            cout << "New Population created! Population #" << GOMEAIndex << " dPopulationSize:" << dPopulationSize << endl;
            cout << this;
        #endif
}

Population::~Population()
{
    for (size_t i = 0; i < dPopulationSize; ++i)
    {
        delete dPopulation[i];
        delete offspringdPopulation[i];
    }
}

ostream & operator << (ostream &out, const Population &populationInstance)
{
    out << "Generation " << populationInstance.numberOfGenerations << ":" << endl;
    for (size_t i = 0; i < populationInstance.dPopulationSize; ++i)
        out << *populationInstance.dPopulation[i] << endl;
    out << endl;
    return out;
}

bool Population::allSolutionsAreEqual()
{
    for (size_t i = 1; i < dPopulationSize; i++)
    {
        for (int j = 0; j < problemInstance->number_of_variables; j++)
        {
            if (dPopulation[i]->variables[j] != dPopulation[0]->variables[j])
                return false;
        }
    }
    return true;
}

void Population::calculateAverageFitness()
{
    averageFitness = 0.0;
    for (size_t i = 0; i < dPopulationSize; ++i)
        averageFitness += dPopulation[i]->getObjectiveValue();
    averageFitness /= dPopulationSize;
}

double Population::getFitnessMean()
{
	double objective_avg = 0.0;
	for(int i = 0; i < dPopulationSize; i++ )
		objective_avg  += dPopulation[i]->getObjectiveValue();
	objective_avg = objective_avg / ((double) dPopulationSize);
	return( objective_avg );
}

double Population::getFitnessVariance()
{
	double objective_avg = getFitnessMean();
	double objective_var = 0.0;
	for(int i = 0; i < dPopulationSize; i++ )
		objective_var  += (dPopulation[i]->getObjectiveValue()-objective_avg)*(dPopulation[i]->getObjectiveValue()-objective_avg);
	objective_var = objective_var / ((double) dPopulationSize);

	if( objective_var <= 0.0 )
		objective_var = 0.0;
	return( objective_var );
}

double Population::getConstraintValueMean()
{
	double constraint_avg = 0.0;
	for(int i = 0; i < dPopulationSize; i++ )
		constraint_avg  += dPopulation[i]->getConstraintValue();
	constraint_avg = constraint_avg / ((double) dPopulationSize);

	return( constraint_avg );
}

double Population::getConstraintValueVariance()
{
	double constraint_avg = getConstraintValueMean();

	double constraint_var = 0.0;
	for(int i = 0; i < dPopulationSize; i++ )
		constraint_var  += (dPopulation[i]->getConstraintValue()-constraint_avg)*(dPopulation[i]->getConstraintValue()-constraint_avg);
	constraint_var = constraint_var / ((double) dPopulationSize);

	if( constraint_var <= 0.0 )
		constraint_var = 0.0;
	return( constraint_var );
}

solution_t<char> *Population::getBestSolution()
{
	int index_best = 0;
	for(int j = 1; j < dPopulationSize; j++ )
    {
        if( problemInstance->betterFitness( dPopulation[j]->getObjectiveValue(), dPopulation[j]->getConstraintValue(), dPopulation[index_best]->getObjectiveValue(), dPopulation[index_best]->getConstraintValue()) )
		{
			index_best = j;
        }
    }
	return( dPopulation[index_best] );
}

solution_t<char> *Population::getWorstSolution()
{
	int index_worst = 0;
	for(int j = 1; j < dPopulationSize; j++ )
    {
        if( problemInstance->betterFitness( dPopulation[index_worst]->getObjectiveValue(), dPopulation[index_worst]->getConstraintValue(), dPopulation[j]->getObjectiveValue(), dPopulation[j]->getConstraintValue()) )
		{
			index_worst = j;
        }
    }
	return( dPopulation[index_worst] );
}


void Population::copyOffspringToPopulation()
{
    for(size_t i = 0; i < dPopulationSize; i++)
    {
        *dPopulation[i] = *offspringdPopulation[i];
    }
}

void Population::makeOffspring()
{
    if( numberOfGenerations == 0 )
    {
        for (size_t i = 0; i < dPopulationSize; ++i)
            updateElitistAndCheckVTR(dPopulation[i]);
    }

    if( FOSInstance->type == linkage::LINKAGE_TREE )
    {
        if (FOSInstance->is_static)
        {
            if (FOSInstance->size() == 0)
            {
                FOSInstance->learnLinkageTreeFOS(problemInstance->getSimilarityMatrix(FOSInstance->getSimilarityMeasure()), false);
                FOSInstance->initializeDependentSubfunctions( problemInstance->subfunction_dependency_map );
            }
        }
        else
        {
            FOSInstance->learnLinkageTreeFOS(dPopulation, config->alphabetSize );
            FOSInstance->initializeDependentSubfunctions(problemInstance->subfunction_dependency_map);
        }
    }

    FOSInstance->setCountersToZero();
    if (config->AnalyzeFOS)
    {
        FOSInstance->writeToFileFOS(config->folder, GOMEAIndex, numberOfGenerations);
    }
 
    generateOffspring();

    if (config->AnalyzeFOS)
        FOSInstance->writeFOSStatistics(config->folder, GOMEAIndex, numberOfGenerations);

}

void Population::generateOffspring()
{
    vec_t<vec_t<int> > neighbors;
   
	assert( !config->useParallelFOSOrder || !config->fixFOSOrderForPopulation );
   	if( config->fixFOSOrderForPopulation )
    	FOSInstance->shuffleFOS(); 
	else if( config->useParallelFOSOrder )
    {
        assert( problemInstance->hasVariableInteractionGraph() );
		FOSInstance->determineParallelFOSOrder(problemInstance->variable_interaction_graph );
    }

    for (size_t i = 0; i < dPopulationSize; i++)
    {
        if (!config->useParallelFOSOrder && !config->fixFOSOrderForPopulation)
            FOSInstance->shuffleFOS();

        solution_t<char> backup = *dPopulation[i];

        bool solutionHasChanged;
        solutionHasChanged = GOM(i);

        /* Phase 2 (Forced Improvement): optimal mixing with elitist solution */
        if (config->useForcedImprovements)
        {
            if ((!solutionHasChanged) || (noImprovementStretches[i] > (1 + (log(dPopulationSize) / log(10)))))
                FI(i);
        }
        // RUBEN testing this setup. Previously, this seemed hardcoded for minimization.
        if((problemInstance->optimization_mode == opt_mode::MIN && !(offspringdPopulation[i]->getObjectiveValue() > dPopulation[i]->getObjectiveValue()))
        || (problemInstance->optimization_mode == opt_mode::MAX && !(offspringdPopulation[i]->getObjectiveValue() < dPopulation[i]->getObjectiveValue())))
            noImprovementStretches[i]++;
        else
            noImprovementStretches[i] = 0;
        
        // if (!(offspringdPopulation[i]->getObjectiveValue() > dPopulation[i]->getObjectiveValue())) // RUBEN hardcoded, assuming minimization?
        //     noImprovementStretches[i]++;
        // else
        //     noImprovementStretches[i] = 0;
    }
}

bool Population::GOM(size_t offspringIndex)
{
    size_t donorIndex;
    bool solutionHasChanged = false;
    bool thisIsTheElitistSolution = *offspringdPopulation[offspringIndex] == sharedInformationPointer->elitist;//(sharedInformationPointer->elitistSolutionGOMEAIndex == GOMEAIndex) && (sharedInformationPointer->elitistSolutionOffspringIndex == offspringIndex);
    
    *offspringdPopulation[offspringIndex] = *dPopulation[offspringIndex];
            
    vec_t<int> donorIndices(dPopulationSize);
    iota(donorIndices.begin(), donorIndices.end(), 0);

    for (size_t i = 0; i < FOSInstance->size(); i++)
    {
        int ind = FOSInstance->FOSorder[i];

        if (FOSInstance->elementSize(ind) == 0 || (int) FOSInstance->elementSize(ind) == problemInstance->number_of_variables)
            continue;

        bool donorEqualToOffspring = true;
        size_t indicesTried = 0;

        while (donorEqualToOffspring && indicesTried < donorIndices.size())
        {
            int j = gomea::utils::rng() % (donorIndices.size() - indicesTried);
            swap(donorIndices[indicesTried], donorIndices[indicesTried + j]);
            donorIndex = donorIndices[indicesTried];
            indicesTried++;

            if (donorIndex == offspringIndex)
                continue;

            vec_t<char> donorGenes;
            for(size_t j = 0; j < FOSInstance->elementSize(ind); j++)
            {
                int variableFromFOS = FOSInstance->FOSStructure[ind][j];
                //offspringPopulation[offspringIndex]->variables[variableFromFOS] = population[donorIndex]->variables[variableFromFOS];
                donorGenes.push_back(dPopulation[donorIndex]->variables[variableFromFOS]);
                if (donorGenes[j] != offspringdPopulation[offspringIndex]->variables[variableFromFOS]) // RUBEN: shouldn't donorGenes[j] be donorGenes[variableFromFOS]?
                    donorEqualToOffspring = false;
            }
            partial_solution_t<char> *partial_offspring = new partial_solution_t<char>(donorGenes, FOSInstance->FOSStructure[ind]);

            if (!donorEqualToOffspring)
            {
                //evaluateSolution(offspringPopulation[offspringIndex], backup, touchedGenes, backup->getObjectiveValue());
                //problemInstance->evaluatePartialSolution(offspringPopulation[offspringIndex], partial_offspring, FOSInstance->getDependentSubfunctions(ind) );
                problemInstance->evaluatePartialSolution(offspringdPopulation[offspringIndex], partial_offspring );

                // accept the change if this solution is not the elitist and the fitness is at least equally good (allows random walk in neutral fitness landscape)
                // however, if this is the elitist solution, only accept strict improvements, to avoid convergence problems
                if ((!thisIsTheElitistSolution && (partial_offspring->getObjectiveValue() >= offspringdPopulation[offspringIndex]->getObjectiveValue())) || 
                        (thisIsTheElitistSolution && (partial_offspring->getObjectiveValue() > offspringdPopulation[offspringIndex]->getObjectiveValue())))     
                {
                    offspringdPopulation[offspringIndex]->insertPartialSolution(partial_offspring);
                    // offspringPopulation[offspringIndex]->variables[variableFromFOS] = population[donorIndex]->variables[variableFromFOS];
                    //*backup = *offspringPopulation[offspringIndex];
                    
                    solutionHasChanged = true;
                    updateElitistAndCheckVTR(offspringdPopulation[offspringIndex]);

                    FOSInstance->improvementCounters[ind]++;
                }

                FOSInstance->usageCounters[ind]++;

            }
            delete partial_offspring;

            break;
        }
    }
    return solutionHasChanged;
}


bool Population::FI(size_t offspringIndex)
{
    if (!config->useParallelFOSOrder && !config->fixFOSOrderForPopulation)
        FOSInstance->shuffleFOS();

    bool solutionHasChanged = 0;

    for (size_t i = 0; i < FOSInstance->size(); i++)
    {
        int ind = FOSInstance->FOSorder[i];
        vec_t<char> touchedGenes = vec_t<char>(FOSInstance->elementSize(ind));
        bool donorEqualToOffspring = true;
        for (size_t j = 0; j < FOSInstance->elementSize(ind); j++)
        {
            int variableFromFOS = FOSInstance->FOSStructure[ind][j];
            touchedGenes[j] = sharedInformationPointer->elitist.variables[variableFromFOS];
            if (dPopulation[offspringIndex]->variables[variableFromFOS] != touchedGenes[j])
                donorEqualToOffspring = false;
        }
        gomea::partial_solution_t<char> *partial_offspring = new gomea::partial_solution_t<char>(touchedGenes, FOSInstance->FOSStructure[ind]);

        if (!donorEqualToOffspring)
        {
            problemInstance->evaluatePartialSolution(offspringdPopulation[offspringIndex], partial_offspring );

            if (partial_offspring->getObjectiveValue() > offspringdPopulation[offspringIndex]->getObjectiveValue() ) 
            {
                offspringdPopulation[offspringIndex]->insertPartialSolution(partial_offspring);
                updateElitistAndCheckVTR(offspringdPopulation[offspringIndex]);
                solutionHasChanged = true;
            }
        }
        delete partial_offspring;
        if (solutionHasChanged)
            break;
    }

    if (!solutionHasChanged)
    {
        *offspringdPopulation[offspringIndex] = sharedInformationPointer->elitist;
    }

    return solutionHasChanged;
}

/*void Population::evaluateSolution(solution_t<char> *parent, gomea::partial_solution_t<char> *solution ) 
{
    checkTimeLimit();

    //cout << "before eval" << solution -> fitness << endl;
    if (config->usePartialEvaluations && solution != NULL)
    {
        problemInstance->calculateFitnessPartialEvaluations(solution, solutionBefore, touchedGenes, fitnessBefore);
        sharedInformationPointer->numberOfEvaluations += (double)touchedGenes.size() / problemInstance->number_of_variables;
    }
    else
    {
        problemInstance->calculateFitness(solution);
        sharedInformationPointer->numberOfEvaluations += 1;
    }

    updateElitistAndCheckVTR(solution);
}*/

/*void Population::evaluateSolution(solution_t<char> *solution, solution_t<char> *solutionBefore, vec_t<int> &touchedGenes, double fitnessBefore)
{
    checkTimeLimit();

    // Do the actual evaluation
    archiveRecord searchResult;
    
    if (config->saveEvaluations)
        sharedInformationPointer->evaluatedSolutions->checkAlreadyEvaluated(solution->variables, &searchResult);
    
    if (searchResult.isFound)
        solution->setObjectiveValue(searchResult.value);
    else
    { 
        //cout << "before eval" << solution -> fitness << endl;
        if (config->usePartialEvaluations && solutionBefore != NULL)
        {
            assert(0);
            // TODO
            //problemInstance->evaluatePartialSolution(solution, solutionBefore, touchedGenes, fitnessBefore);
            sharedInformationPointer->numberOfEvaluations += (double)touchedGenes.size() / problemInstance->number_of_variables;
        }
        else
        {
            assert(0);
            // TODO
            //problemInstance->evaluate(solution);
            sharedInformationPointer->numberOfEvaluations += 1;
        }

        if (config->saveEvaluations)
            sharedInformationPointer->evaluatedSolutions->insertSolution(solution->variables, solution->getObjectiveValue());
    }

    updateElitistAndCheckVTR(solution);
}*/

void Population::checkTimeLimit()
{
    if ( config->maximumNumberOfSeconds > 0 && utils::getElapsedTimeSeconds(sharedInformationPointer->startTime) > config->maximumNumberOfSeconds)
    {
        terminated = true;
        throw utils::customException("time");
    }
}

void Population::updateElitistAndCheckVTR(solution_t<char> *solution)
{
    /* Update elitist solution */
    //if (sharedInformationPointer->firstEvaluationEver || (solution->getObjectiveValue() > sharedInformationPointer->elitist.getObjectiveValue()))
    if (sharedInformationPointer->firstEvaluationEver || problemInstance->betterFitness(solution,&sharedInformationPointer->elitist) )
    {
        sharedInformationPointer->elitistSolutionHittingTimeMilliseconds = utils::getElapsedTimeMilliseconds(sharedInformationPointer->startTime);
        sharedInformationPointer->elitistSolutionHittingTimeEvaluations = problemInstance->number_of_evaluations;

        sharedInformationPointer->elitist = *solution;
		sharedInformationPointer->elitistFitness = solution->getObjectiveValue();
		sharedInformationPointer->elitistConstraintValue = solution->getConstraintValue();
        
        /* Check the VTR */
        if (problemInstance->use_vtr && solution->getObjectiveValue() >= problemInstance->vtr)
        {
            writeStatisticsToFile(config->folder, sharedInformationPointer->elitistSolutionHittingTimeEvaluations, sharedInformationPointer->elitistSolutionHittingTimeMilliseconds, solution);
            writeElitistSolutionToFile(config->folder, sharedInformationPointer->elitistSolutionHittingTimeEvaluations, sharedInformationPointer->elitistSolutionHittingTimeMilliseconds, solution);
            cout << "VTR HIT!\n";
            terminated = true;
            throw utils::customException("vtr");
        }
    
        //writeStatisticsToFile(config->folder, sharedInformationPointer->elitistSolutionHittingTimeEvaluations, sharedInformationPointer->elitistSolutionHittingTimeMilliseconds, solution);
        //if( config->writeElitists )
			//writeElitistSolutionToFile(config->folder, sharedInformationPointer->elitistSolutionHittingTimeEvaluations, sharedInformationPointer->elitistSolutionHittingTimeMilliseconds, solution);
    }

    sharedInformationPointer->firstEvaluationEver = false;
}


}}