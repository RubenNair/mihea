#include "gomea/src/mixed_integer/Population.hpp"

namespace gomea{
namespace mixedinteger{


Population::Population(Config *config_, fitness_t *problemInstance_, sharedInformation *sharedInformationPointer_, size_t GOMEAIndex_, size_t populationSize_, size_t optimizerIndex, string optimizerName, linkage_model_pt FOSInstance_ ): 
        config(config_), 
        problemInstance(problemInstance_),
        sharedInformationPointer(sharedInformationPointer_),
        GOMEAIndex(GOMEAIndex_), 
        populationSize(populationSize_),
        optimizerIndex(optimizerIndex),
        optimizerName(optimizerName)
{
        terminated = false;
        numberOfGenerations = 0;
        averageFitness = 0.0;

        optimizerElitistFitnessInitialized = false;
        
        population = new Solutionset(populationSize);
        offspringPopulation = new Solutionset(populationSize);
        noImprovementStretches.resize(populationSize);


        for (size_t i = 0; i < populationSize; ++i)
        {
            noImprovementStretches[i] = 0;
            if(config->useBN) 
            {
                // If useBN flag is true, assume we're dealing with BNs. Static cast problemInstance to BNStructureLearning
                fitness::BNStructureLearning *BNproblemInstance = dynamic_cast<fitness::BNStructureLearning*>(problemInstance);
                
                if (BNproblemInstance == nullptr) {
                    throw std::runtime_error("Failed to cast problemInstance to BNStructureLearning");
                }

                vec_t<double> maxValuesData, minValuesData;
                vec_t<vec_t<double>> data = config->data->getDataMatrix().getRawMatrix();
                std::tie(maxValuesData, minValuesData) = findMaxAndMinValuesInData(data);

                population->solutions[i] = new solution_BN(config->numberOfVariables, config->alphabetSize, config->numberOfcVariables, 
                                                    config->data->getColumnType(), BNproblemInstance->getDensity()->getOriginalData()->getColumnNumberOfClasses(), 
                                                    config->discretization_policy_index, config->maxParents, config->maxInstantiations, 
                                                    BNproblemInstance, maxValuesData, minValuesData,
                                                    config->lower_user_range, config->upper_user_range, 
                                                    config->data, (double)i / populationSize, 
                                                    config->useNormalizedCVars, config->transformCVariables, config->useOptimalSolution, config->guaranteedInitSpread, config->extraCVarForNumberOfBins, config->forceNBoundariesUsed, config->problemInstancePath, config->runIndex);

                problemInstance->evaluate(population->solutions[i]);

                // Set offspringPopulation[i] to a new solution_BN, just to initialize the necessary memory.
                offspringPopulation->solutions[i] = new solution_BN(config->numberOfVariables, config->alphabetSize, config->numberOfcVariables, config->data->getColumnType(), BNproblemInstance->getDensity()->getOriginalData()->getColumnNumberOfClasses(), config->discretization_policy_index, config->maxParents, config->maxInstantiations, BNproblemInstance, maxValuesData, minValuesData, config->lower_user_range, config->upper_user_range, config->data);
                *offspringPopulation->solutions[i] = *population->solutions[i];
            } else
            {
                population->solutions[i] = new solution_mixed(config->numberOfVariables, config->alphabetSize, config->numberOfcVariables, problemInstance);
                population->solutions[i]->randomInit(&gomea::utils::rng);

                problemInstance->evaluate(population->solutions[i]);
                
                offspringPopulation->solutions[i] = new solution_mixed(config->numberOfVariables, config->alphabetSize, config->numberOfcVariables, problemInstance);
                *offspringPopulation->solutions[i] = *population->solutions[i];
            }
            

            // Update elitist based on these new solutions
            updateElitistAndCheckVTR(population->solutions[i]);
        }



		if( config->linkage_config != NULL && config->numberOfVariables > 0 )
		{
			FOSInstance = linkage_model_t::createFOSInstance( *config->linkage_config, problemInstance->number_of_variables );
            FOSInstance->initializeDependentSubfunctions( problemInstance->subfunction_dependency_map );
		}
		else if( FOSInstance_ == NULL )
		{
			FOSInstance = linkage_model_t::createLinkageTreeFOSInstance(config->FOSIndex, problemInstance->number_of_variables, config->linkage_config->lt_similarity_measure, config->linkage_config->lt_maximum_set_size);
		}
		else FOSInstance = FOSInstance_;
        
        if(config->numberOfcVariables > 0)
        {
            if(config->dontUseOffspringPopulation)
                iamalgamInstance = new iamalgam(config_, population);
            else
                iamalgamInstance = new iamalgam(config_, offspringPopulation);
            iamalgamInstance->initialize();
        }



        #ifdef DEBUG
            cout << "New Population created! Population #" << GOMEAIndex << " populationSize:" << populationSize << endl;
            cout << this;
        #endif
}

Population::~Population()
{
    for (size_t i = 0; i < populationSize; ++i)
    {
        delete population->solutions[i];
        delete offspringPopulation->solutions[i];
    }
    delete population;
    delete offspringPopulation;

    iamalgamInstance->ezilaitini();
    delete iamalgamInstance;
}

ostream & operator << (ostream &out, const Population &populationInstance)
{
    out << "Generation " << populationInstance.numberOfGenerations << ":" << endl;
    for (size_t i = 0; i < populationInstance.populationSize; ++i)
        out << *populationInstance.population->solutions[i] << endl;
    out << endl;
    return out;
}

bool Population::allSolutionsAreEqual()
{
    for (size_t i = 1; i < populationSize; i++)
    {
        if(*population->solutions[i] != *population->solutions[0])
            return false;
    }
    cout << "[DEBUGGING] All solutions are equal! Population #" << GOMEAIndex << endl;
    return true;
}

void Population::calculateAverageFitness()
{
    averageFitness = 0.0;
    for (size_t i = 0; i < populationSize; ++i)
        averageFitness += population->solutions[i]->getObjectiveValue();
    averageFitness /= populationSize;
}

double Population::getFitnessMean()
{
	double objective_avg = 0.0;
	for(size_t i = 0; i < populationSize; i++ )
		objective_avg  += population->solutions[i]->getObjectiveValue();
	objective_avg = objective_avg / ((double) populationSize);
	return( objective_avg );
}

double Population::getFitnessVariance()
{
	double objective_avg = getFitnessMean();
	double objective_var = 0.0;
	for(size_t i = 0; i < populationSize; i++ )
		objective_var  += (population->solutions[i]->getObjectiveValue()-objective_avg)*(population->solutions[i]->getObjectiveValue()-objective_avg);
	objective_var = objective_var / ((double) populationSize);

	if( objective_var <= 0.0 )
		objective_var = 0.0;
	return( objective_var );
}

double Population::getConstraintValueMean()
{
	double constraint_avg = 0.0;
	for(size_t i = 0; i < populationSize; i++ )
		constraint_avg  += population->solutions[i]->getConstraintValue();
	constraint_avg = constraint_avg / ((double) populationSize);

	return( constraint_avg );
}

double Population::getConstraintValueVariance()
{
	double constraint_avg = getConstraintValueMean();

	double constraint_var = 0.0;
	for(size_t i = 0; i < populationSize; i++ )
		constraint_var  += (population->solutions[i]->getConstraintValue()-constraint_avg)*(population->solutions[i]->getConstraintValue()-constraint_avg);
	constraint_var = constraint_var / ((double) populationSize);

	if( constraint_var <= 0.0 )
		constraint_var = 0.0;
	return( constraint_var );
}

solution_mixed *Population::getBestSolution()
{
	int index_best = 0;
	for(size_t j = 1; j < populationSize; j++ )
    {
        if( problemInstance->betterFitness( population->solutions[j]->getObjectiveValue(), population->solutions[j]->getConstraintValue(), population->solutions[index_best]->getObjectiveValue(), population->solutions[index_best]->getConstraintValue()) )
		{
			index_best = j;
        }
    }
	return( population->solutions[index_best] );
}

solution_mixed *Population::getWorstSolution()
{
	int index_worst = 0;
	for(size_t j = 1; j < populationSize; j++ )
    {
        if( problemInstance->betterFitness( population->solutions[index_worst]->getObjectiveValue(), population->solutions[index_worst]->getConstraintValue(), population->solutions[j]->getObjectiveValue(), population->solutions[j]->getConstraintValue()) )
		{
			index_worst = j;
        }
    }
	return( population->solutions[index_worst] );
}


void Population::copyOffspringToPopulation()
{
    for(size_t i = 0; i < populationSize; i++)
    {
        delete population->solutions[i];
        population->solutions[i] = offspringPopulation->solutions[i]->clone();
    }
}

void Population::copyPopulationToOffspring()
{
    for(size_t i = 0; i < populationSize; i++)
    {
        delete offspringPopulation->solutions[i];
        offspringPopulation->solutions[i] = population->solutions[i]->clone();
    }
}

void Population::makeOffspring()
{
    if( numberOfGenerations == 0 )
    {
        for (size_t i = 0; i < populationSize; ++i)
            updateElitistAndCheckVTR(population->solutions[i]);
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
            vec_t<solution_t<int>*> casted_population;
            for(solution_mixed* sol : population->solutions) 
            {
                casted_population.push_back(static_cast<solution_t<int>*>(sol));
            }

            FOSInstance->learnLinkageTreeFOS(casted_population, config->alphabetSize );
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

void Population::learnDiscreteModel()
{
    calculateAverageFitness();

    if( numberOfGenerations == 0 )
    {
        for (size_t i = 0; i < populationSize; ++i)
            updateElitistAndCheckVTR(population->solutions[i]);
    }

    if( FOSInstance->type == linkage::LINKAGE_TREE )
    {
        if (FOSInstance->is_static)
        {
            if (FOSInstance->size() == 0)
            {
                FOSInstance->learnLinkageTreeFOS(problemInstance->getSimilarityMatrix(FOSInstance->getSimilarityMeasure()), false);
            }
        }
        else
        {
            vec_t<solution_t<int>*> casted_population;
            for(solution_mixed* sol : population->solutions) 
            {
                casted_population.push_back(static_cast<solution_t<int>*>(sol));
            }

            FOSInstance->learnLinkageTreeFOS(casted_population, config->alphabetSize );
        }

        // For mixed-integer problems, the FOS set of all discrete elements is still valid (since it is not the whole solution (continuous part), and thus not a complete copy would be made)
        // So add the FOS set of all discrete elements to the FOSInstance.
        if(config->numberOfcVariables > 0)
        {
            vec_t<int> allDiscreteElements(problemInstance->number_of_variables);
            for (int i = 0; i < problemInstance->number_of_variables; ++i)
            {
                allDiscreteElements[i] = i;
            }
            FOSInstance->FOSStructure.push_back(allDiscreteElements);
        }
        FOSInstance->initializeDependentSubfunctions(problemInstance->subfunction_dependency_map);

    }

    FOSInstance->setCountersToZero();
}

void Population::determineFOSOrder()
{
    assert( !config->useParallelFOSOrder || !config->fixFOSOrderForPopulation );
   	if( config->fixFOSOrderForPopulation )
    	FOSInstance->shuffleFOS(); 
	else if( config->useParallelFOSOrder )
    {
        assert( problemInstance->hasVariableInteractionGraph() );
		FOSInstance->determineParallelFOSOrder(problemInstance->variable_interaction_graph );
    }

    if (!config->useParallelFOSOrder && !config->fixFOSOrderForPopulation)
            FOSInstance->shuffleFOS();
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
    int FIcount = 0;
    for (size_t i = 0; i < populationSize; i++)
    {
        if (!config->useParallelFOSOrder && !config->fixFOSOrderForPopulation)
            FOSInstance->shuffleFOS();

        solution_mixed backup = *population->solutions[i];

        bool solutionHasChanged;
        solutionHasChanged = GOM(i);

        /* Phase 2 (Forced Improvement): optimal mixing with elitist solution */
        if (config->useForcedImprovements)
        {
            if ((!solutionHasChanged) || (noImprovementStretches[i] > (1 + (log(populationSize) / log(10)))))
            {
                FI(i);
                FIcount++;
            }
        }
        
        if (!(offspringPopulation->solutions[i]->getObjectiveValue() < population->solutions[i]->getObjectiveValue())) // NOTE: assuming minimization
            noImprovementStretches[i]++;
        else
            noImprovementStretches[i] = 0;
    }
}

// Generate offspring for a single FOS element
void Population::generateSingleOffspring(int FOS_index)
{
    int FIcount = 0;
    for (size_t i = 0; i < populationSize; i++)
    {
        
        solution_mixed backup = *offspringPopulation->solutions[i];

        bool solutionHasChanged;
        solutionHasChanged = GOMSingleFOS(i, FOS_index);

        /* Phase 2 (Forced Improvement): optimal mixing with elitist solution */
        if (config->useForcedImprovements)
        {
            if ((!solutionHasChanged) || (noImprovementStretches[i] > (1 + (log(populationSize) / log(10)))))
            {
                FISingleFOS(i, FOS_index);
                FIcount++;
            }
        }
        
        if (!(offspringPopulation->solutions[i]->getObjectiveValue() < population->solutions[i]->getObjectiveValue())) // NOTE: assuming minimization
            noImprovementStretches[i]++;
        else
            noImprovementStretches[i] = 0;
    }
}

// NOTE: Confusing name perhaps, now only updates elitist and checks VTR.
void Population::evaluateAllSolutionsInPopulation()
{
    for (size_t i = 0; i < populationSize; i++)
    {
        updateElitistAndCheckVTR(offspringPopulation->solutions[i]);
    }
}

// Original GOM function, not used for MIHEA
bool Population::GOM(size_t offspringIndex)
{
    size_t donorIndex;
    bool solutionHasChanged = false;
    bool thisIsTheElitistSolution = *offspringPopulation->solutions[offspringIndex] == *sharedInformationPointer->elitist;
    
    *offspringPopulation->solutions[offspringIndex] = *population->solutions[offspringIndex];
            
    vec_t<int> donorIndices(populationSize);
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

            vec_t<int> donorGenes;
            for(size_t j = 0; j < FOSInstance->elementSize(ind); j++)
            {
                int variableFromFOS = FOSInstance->FOSStructure[ind][j];
                donorGenes.push_back(population->solutions[donorIndex]->variables[variableFromFOS]);
                if (donorGenes[variableFromFOS] != offspringPopulation->solutions[offspringIndex]->variables[variableFromFOS])
                    donorEqualToOffspring = false;
            }
            partial_solution_t<int> *partial_offspring = new partial_solution_t<int>(donorGenes, FOSInstance->FOSStructure[ind]);

            if (!donorEqualToOffspring)
            {
                problemInstance->evaluatePartialSolution(offspringPopulation->solutions[offspringIndex], partial_offspring );

                // accept the change if this solution is not the elitist and the fitness is at least equally good (allows random walk in neutral fitness landscape)
                // however, if this is the elitist solution, only accept strict improvements, to avoid convergence problems
                if ((!thisIsTheElitistSolution && (partial_offspring->getObjectiveValue() >= offspringPopulation->solutions[offspringIndex]->getObjectiveValue())) || 
                        (thisIsTheElitistSolution && (partial_offspring->getObjectiveValue() > offspringPopulation->solutions[offspringIndex]->getObjectiveValue())))     
                {
                    offspringPopulation->solutions[offspringIndex]->insertPartialSolution(partial_offspring);
                    
                    solutionHasChanged = true;
                    updateElitistAndCheckVTR(offspringPopulation->solutions[offspringIndex]);

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

bool Population::GOMSingleFOS(size_t offspringIndex, size_t FOSIndex)
{
    size_t donorIndex;
    bool solutionHasChanged = false;
    bool thisIsTheElitistSolution = *offspringPopulation->solutions[offspringIndex] == *sharedInformationPointer->elitist;
            
    vec_t<int> donorIndices(populationSize);
    iota(donorIndices.begin(), donorIndices.end(), 0);

    int ind = FOSInstance->FOSorder[FOSIndex];

    if (FOSInstance->elementSize(ind) == 0)
        return solutionHasChanged;

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

        vec_t<int> donorGenes = vec_t<int>(offspringPopulation->solutions[offspringIndex]->variables.size());
        // First copy all original genes to donorGenes
        for(size_t i = 0; i < offspringPopulation->solutions[offspringIndex]->variables.size(); i++)
            donorGenes[i] = offspringPopulation->solutions[offspringIndex]->variables[i];

        // Then copy the genes from the donor
        for(size_t j = 0; j < FOSInstance->elementSize(ind); j++)
        {
            int variableFromFOS = FOSInstance->FOSStructure[ind][j];
            donorGenes[variableFromFOS] = population->solutions[donorIndex]->variables[variableFromFOS];
            if (donorGenes[variableFromFOS] != offspringPopulation->solutions[offspringIndex]->variables[variableFromFOS])
                donorEqualToOffspring = false;
        }
        

        if (!donorEqualToOffspring)
        {   
            solution_mixed *offspring;
            if(config->useBN)
            {
                solution_BN *offspring_BN = dynamic_cast<solution_BN *>(offspringPopulation->solutions[offspringIndex])->clone();
                offspring_BN->reProcessParametersSolution(donorGenes);
                offspring = offspring_BN;
            } else {
                offspring = new solution_mixed(donorGenes, offspringPopulation->solutions[offspringIndex]->c_variables);
            }
            
            problemInstance->evaluate(offspring);

            // accept the change if this solution is not the elitist and the fitness is at least equally good (allows random walk in neutral fitness landscape)
            // however, if this is the elitist solution, only accept strict improvements, to avoid convergence problems
            // NOTE: assuming minimization
            if ((!thisIsTheElitistSolution && (offspring->getObjectiveValue() <= offspringPopulation->solutions[offspringIndex]->getObjectiveValue())) || 
                    (thisIsTheElitistSolution && (offspring->getObjectiveValue() < offspringPopulation->solutions[offspringIndex]->getObjectiveValue())))     
            {
                delete offspringPopulation->solutions[offspringIndex];
                offspringPopulation->solutions[offspringIndex] = offspring->clone();

                // For GAMBIT_K, also update population with better individual (so it can also be chosen as donor)
                if(config->dontUseOffspringPopulation) {
                    delete population->solutions[offspringIndex];
                    population->solutions[offspringIndex] = offspring->clone();
                }
                solutionHasChanged = true;
                updateElitistAndCheckVTR(offspringPopulation->solutions[offspringIndex]);

                FOSInstance->improvementCounters[ind]++;
            }
            
            delete offspring;

            FOSInstance->usageCounters[ind]++;

        }

        break;
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
        vec_t<int> touchedGenes = vec_t<int>(FOSInstance->elementSize(ind));
        bool donorEqualToOffspring = true;
        for (size_t j = 0; j < FOSInstance->elementSize(ind); j++)
        {
            int variableFromFOS = FOSInstance->FOSStructure[ind][j];
            touchedGenes[j] = sharedInformationPointer->elitist->variables[variableFromFOS];
            if (population->solutions[offspringIndex]->variables[variableFromFOS] != touchedGenes[j])
                donorEqualToOffspring = false;
        }
        gomea::partial_solution_t<int> *partial_offspring = new gomea::partial_solution_t<int>(touchedGenes, FOSInstance->FOSStructure[ind]);

        if (!donorEqualToOffspring)
        {
            problemInstance->evaluatePartialSolution(offspringPopulation->solutions[offspringIndex], partial_offspring );

            if (partial_offspring->getObjectiveValue() > offspringPopulation->solutions[offspringIndex]->getObjectiveValue() ) 
            {
                offspringPopulation->solutions[offspringIndex]->insertPartialSolution(partial_offspring);
                updateElitistAndCheckVTR(offspringPopulation->solutions[offspringIndex]);
                solutionHasChanged = true;
            }
        }
        delete partial_offspring;
        if (solutionHasChanged)
            break;
    }

    if (!solutionHasChanged)
    {
        *offspringPopulation->solutions[offspringIndex] = *sharedInformationPointer->elitist;
    }

    return solutionHasChanged;
}

bool Population::FISingleFOS(size_t offspringIndex, size_t FOSIndex)
{
    bool solutionHasChanged = 0;

    int ind = FOSInstance->FOSorder[FOSIndex];
    bool donorEqualToOffspring = true;

    vec_t<int> touchedGenes = vec_t<int>(offspringPopulation->solutions[offspringIndex]->variables.size());
    // First copy all original genes to touchedGenes
    for(size_t i = 0; i < offspringPopulation->solutions[offspringIndex]->variables.size(); i++)
        touchedGenes[i] = offspringPopulation->solutions[offspringIndex]->variables[i];

    for (size_t j = 0; j < FOSInstance->elementSize(ind); j++)
    {
        int variableFromFOS = FOSInstance->FOSStructure[ind][j];
        touchedGenes[variableFromFOS] = sharedInformationPointer->elitist->variables[variableFromFOS];
        if (offspringPopulation->solutions[offspringIndex]->variables[variableFromFOS] != touchedGenes[variableFromFOS]) 
            donorEqualToOffspring = false;
    }
    solution_mixed *offspring = new solution_mixed(touchedGenes, offspringPopulation->solutions[offspringIndex]->c_variables);

    if (!donorEqualToOffspring)
    {
        problemInstance->evaluate(offspring);

        if (offspring->getObjectiveValue() < offspringPopulation->solutions[offspringIndex]->getObjectiveValue() ) // NOTE: assuming minimization
        {
            offspringPopulation->solutions[offspringIndex] = offspring;
            updateElitistAndCheckVTR(offspringPopulation->solutions[offspringIndex]);
            solutionHasChanged = true;
        }
    }

    if (!solutionHasChanged)
    {
        *offspringPopulation->solutions[offspringIndex] = *sharedInformationPointer->elitist;
    }

    return solutionHasChanged;
}

void Population::checkTimeLimit()
{
    if ( config->maximumNumberOfSeconds > 0 && utils::getElapsedTimeSeconds(sharedInformationPointer->startTime) > config->maximumNumberOfSeconds)
    {
        terminated = true;
        throw utils::customException("time");
    }
}

void Population::updateElitistAndCheckVTR(solution_mixed *solution)
{
    // First update 'local' elitist (of current population) fitness
        if(!optimizerElitistFitnessInitialized || solution->getObjectiveValue() < optimizerElitistFitness) // NOTE: assuming minimization
        {
            optimizerElitistFitness = solution->getObjectiveValue();            
            optimizerElitistFitnessInitialized = true;
        }
    /* Update elitist solution */
    if (sharedInformationPointer->firstEvaluationEver || problemInstance->betterFitness(solution,sharedInformationPointer->elitist) )
    {
        if(config->printNewElitists || config->useOptimalSolution)
        {
            cout << "New elitist solution found: " << solution->getObjectiveValue() << ", c_variables: ";
            for(auto const& val : solution->c_variables)
                cout << val << " ";
            cout << "\t variables: ";
            for(int dval : solution->variables)
                cout << dval << " ";
            cout << endl;
        }
        sharedInformationPointer->elitistSolutionHittingTimeMilliseconds = utils::getElapsedTimeMilliseconds(sharedInformationPointer->startTime);
        sharedInformationPointer->elitistSolutionHittingTimeEvaluations = problemInstance->full_number_of_evaluations;

        // Replace sharedInformationPointer->elitist with (a new copy of) the solution
        delete sharedInformationPointer->elitist;
        sharedInformationPointer->elitist = solution->clone();
        
		sharedInformationPointer->elitistFitness = solution->getObjectiveValue();
		sharedInformationPointer->elitistConstraintValue = solution->getConstraintValue();
        sharedInformationPointer->optimizerIndex = this->optimizerIndex;
        
        /* Check the VTR */
        if (problemInstance->use_vtr && solution->getObjectiveValue() <= problemInstance->vtr + 1e-10) // NOTE: assuming minimization
        {
            writeStatisticsToFile(config->folder, sharedInformationPointer->elitistSolutionHittingTimeEvaluations, sharedInformationPointer->elitistSolutionHittingTimeMilliseconds, solution, populationSize, true);
            writeElitistSolutionToFile(config->folder, sharedInformationPointer->elitistSolutionHittingTimeEvaluations, sharedInformationPointer->elitistSolutionHittingTimeMilliseconds, solution);
            cout << "VTR HIT! (popsize: " << populationSize << ", stats -> (generation #evals): " << numberOfGenerations << " " << sharedInformationPointer->elitistSolutionHittingTimeEvaluations << ")\n";
            terminated = true;
            throw utils::customException("vtr");
        }
    }

    sharedInformationPointer->firstEvaluationEver = false;
}

void Population::generateDiscretePopulation(int FOS_index) 
{
    generateSingleOffspring(FOS_index);
}

// Learning continuous model part from 2014 paper MIHEA.
void Population::learnContinuousModel()
{
    iamalgamInstance->learnContinuousModel();

    iamalgamInstance->generateAndEvaluateNewSolutionsToFillPopulations();
    writePopulationToFile(config->folder, population->solutions, "POPULATION After generating CONTINUOUS population ----------------------------------", config->logDebugInformation);
    
    evaluateAllSolutionsInPopulation();

    iamalgamInstance->number_of_generations++;

    iamalgamInstance->adaptDistributionMultipliersForOnePopulation();
}

void Population::checkDiscretePopulationConvergedNotOptimal()
{
    bool allOnes = true;
    bool allEqual = true;
    for(size_t i = 1; i < populationSize; i++)
    {
        for(size_t j = 0; j < population->solutions[i]->variables.size(); j++)
        {
            if(population->solutions[i]->variables[j] != population->solutions[0]->variables[j])
            {
                allEqual = false;
                break;
            }
            if(allOnes && population->solutions[i]->variables[j] == 0)
            {
                allOnes = false;
            }
        }
        if(!allEqual)
            break;
    }

    if(allEqual && !allOnes)
    {
        cout << "Discrete parts of all solutions in population are equal, but not all 1. Optimum cannot be reached, terminating." << endl;
        terminated = true;
    }
}

bool Population::allContinuousEqual(solution_mixed *a, solution_mixed *b)
{
    for(size_t i = 0; i < a->c_variables.size(); i++)
    {
        if(a->c_variables[i] != b->c_variables[i])
            return false;
    }
    return true;
}

}}