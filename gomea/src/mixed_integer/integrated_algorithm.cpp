#include "gomea/src/mixed_integer/integrated_algorithm.hpp"

// main()
// {
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
//                 newcPopulation[i] = GenerateContinuousPart(cPopulation[i]);
//                 newdPopulation[i] = GenerateDiscretePart(j, newdPopulation[i], Population);
//             }
//         }
//         cPopulation = newcPopulation;
//         dPopulation = newdPopulation;
//     }
// }
