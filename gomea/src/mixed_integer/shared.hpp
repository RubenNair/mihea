#pragma once

#include <vector>
using namespace std;

#include "gomea/src/mixed_integer/utils.hpp"
#include "gomea/src/utils/time.hpp"
#include "gomea/src/common/solution_mixed.hpp"

namespace gomea{
namespace mixedinteger{

class sharedInformation
{
	public:
		time_t startTime;
		double elitistSolutionHittingTimeMilliseconds,
			   elitistSolutionHittingTimeEvaluations;

		solutionsArchive *evaluatedSolutions; // RUBEN: ignoring for now, seems to be unused currently anyway
		bool firstEvaluationEver;
		double elitistFitness;
		double elitistConstraintValue;
		solution_mixed *elitist = new solution_mixed(1,2,1);

    sharedInformation(int maxArchiveSize)
    {
        startTime = utils::getTimestamp();
        firstEvaluationEver = true;
        evaluatedSolutions = new solutionsArchive(maxArchiveSize);
        gomea::utils::clearTimers();
    }

    ~sharedInformation()
    {
        delete evaluatedSolutions;
    }
};

}}