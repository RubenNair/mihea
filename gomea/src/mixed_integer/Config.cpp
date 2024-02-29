#include "gomea/src/mixed_integer/Config.hpp"

#include "gomea/src/fitness/benchmarks-discrete.hpp"
#include "gomea/src/fitness/benchmarks-mixed.hpp"
#include "gomea/src/fitness/so_benchmarks.h"

namespace gomea{
namespace mixedinteger{

Config::Config(){}

Config::~Config()
{
    delete fitness;
    delete linkage_config;
}

/**
 * Parses the command line.
 * For options, see printUsage.
 */
void Config::splitString(const string &str, vector<string> &splitted, char delim)
{
    size_t current, previous = 0;
    current = str.find(delim);
    while (current != string::npos)
    {
        splitted.push_back(str.substr(previous, current - previous));
        previous = current + 1;
        current = str.find(delim, previous);
    }
    splitted.push_back(str.substr(previous, current - previous));
}

bool Config::isNumber(const string &str)
{
    return !str.empty() && all_of(str.begin(), str.end(), ::isdigit);
}

fitness_t *Config::getFitnessClassDiscrete( int problem_index, int number_of_variables)
{
	switch(problem_index) 
	{
        case 1:
        return new gomea::fitness::oneMaxSphere_t(number_of_variables, numberOfcVariables);
        case 2:
        return new gomea::fitness::oneMaxRotEllip_t(number_of_variables, numberOfcVariables);
        case 3:
        return new gomea::fitness::DT5Sphere_t(number_of_variables, numberOfcVariables);
        case 4:
        return new gomea::fitness::DT5RotEllip_t(number_of_variables, numberOfcVariables); //new gomea::fitness::maxCut_t(number_of_variables, problemInstancePath);
        case 5:
        return new gomea::fitness::DT5BlockRotEllip_t(number_of_variables, numberOfcVariables, a_value);
        case 50:
        return new gomea::fitness::DT5BlockRotEllipCentersZero_t(number_of_variables, numberOfcVariables, a_value);
        case 51:
        return new gomea::fitness::DT5BlockFullRotEllipWrongExponent_t(number_of_variables, numberOfcVariables, a_value);
        case 55:
        return new gomea::fitness::DT5BlockNOTRotEllip_t(number_of_variables, numberOfcVariables, a_value);
        case 551:
        return new gomea::fitness::DT5BlockNOTRotEllipWrongExponent_t(number_of_variables, numberOfcVariables, a_value);
        case 1000 ... 99999999:
        return new gomea::fitness::BNStructureLearning(numberOfdVariables, numberOfcVariables, problem_index, data, maxParents, maxDiscretizations, transformCVariables);
		default:
		return NULL;
	}
	return new gomea::fitness::oneMaxSphere_t(number_of_variables, numberOfcVariables);
}

void Config::setMethodInitParams(int settingIndex)
{
    // TODO set parameters for each method, easier than giving individual flags in command line.
    switch(settingIndex)
    {
        case 0:
        // Method 1, init 1
        // No flags needed, default settings
        break;
        case 1:
        // Method 1, init 2
        useNormalizedCVars = true;
        break;
        case 2:
        // Method 1, init 3
        guaranteedInitSpread = true;
        break;
        case 3:
        // Method 2, init 1
        transformCVariables = true;
        break;
        case 4:
        // Method 2, init 2
        transformCVariables = true;
        useNormalizedCVars = true;
        break;
        case 5:
        // Method 2, init 3
        transformCVariables = true;
        guaranteedInitSpread = true;
        break;
        case 6:
        // Method 3, init 1
        extraCVarForNumberOfBins = true;
        break;
        case 7:
        // Method 3, init 2
        extraCVarForNumberOfBins = true;
        useNormalizedCVars = true;
        break;
        case 8:
        // Method 3, init 3
        extraCVarForNumberOfBins = true;
        guaranteedInitSpread = true;
        break;
        default:
        return;
    }
}

void Config::setOptimizerName()
{
    if (useOptimalSolution) {
        optimizerName = "Ground_truth";
        return;
    }
    string name = "Method";

    if(transformCVariables) name += "2";
    else if(extraCVarForNumberOfBins) name += "3";
    else name += "1";

    if(forceNBoundariesUsed) name += "*";

    name += "_init";
    if(useNormalizedCVars) name += "2";
    else if(guaranteedInitSpread) name += "3";
    else name += "1";

    optimizerName = name;
}

bool Config::parseCommandLine(int argc, char **argv)
{
  const struct option longopts[] =
  {
    {"guaranteedInitSpread",        no_argument,        0, 'g'},
    {"maxDiscretizations",          required_argument,  0, 'X'},  
    {"extraCVarForNumberOfBins",    no_argument,        0, 'Y'},  
    {"popUpdatesDuringGOM",         no_argument,        0, 'Q'}, 
    {"analyzeFOS",                  no_argument,        0, 'w'},
    {"writeElitists",               no_argument,        0, 'e'},
    {"printNewElitists",            no_argument,        0, 'e'},
    {"dontUseOffspringPopulation",  no_argument,        0, 'k'},
    {"saveEvaluations",             no_argument,        0, 's'}, 
    {"printHelp",                   no_argument,        0, 'h'},
    {"lower_user_range",            required_argument,  0, 'a'},
    {"upper_user_range",            required_argument,  0, 'b'},
    {"forceNBoundariesUsed",        no_argument,        0, 'f'},
    {"runIndex",                    required_argument,  0, 'r'},
    {"basePopulationSize",          required_argument,  0, 'n'},
    {"maximumNumberOfGAMBITs",      required_argument,  0, 'N'},
    {"ProblemIndex",                required_argument,  0, 'P'},
    {"useParallelGOM",              required_argument,  0, 'p'},
    {"FOSIndex/linkage_config",     required_argument,  0, 'F'},
    {"methodInitParams",            required_argument,  0, 'm'},
    {"maximumNumberOfGenerations",  required_argument,  0, 'M'},
    {"useNormalizedCVars",          no_argument,        0, 'u'},
    {"logDebugInformation",         no_argument,        0, 'l'},
    {"numberOfVariables",           required_argument,  0, 'L'},
    {"useOptimalSolution",          no_argument,        0, 'o'},
    {"outputFolder",                      required_argument,  0, 'O'},
    {"transformCVariables",         no_argument,        0, 't'},
    {"maximumNumberOfSeconds",      required_argument,  0, 'T'},
    {"vtr",                         required_argument,  0, 'V'},
    {"seed",                        required_argument,  0, 'S'},
    {"instance",                    required_argument,  0, 'I'},
    {"LTsimilarityMeasure",           required_argument,  0, 'Z'},
    {"GPUIndex",                    required_argument,  0, 'G'},
    {0,                             0,                  0,  0 }
  };


  int c, index;
  while ((c = getopt_long(argc, argv, "h::k::n::p::X::Y::Q::g::w::e::s::f::r::P::F::m::u::l::L::o::O::t::T::S::V::I::B::Z::G::M::N::E::a::b::", longopts, &index)) != -1)
  {
    switch (c)
    {
        case 'g':
            guaranteedInitSpread = true;
            break;
		case 'X':
			maxDiscretizations = atoi(optarg);
			break;
		case 'Y':
			extraCVarForNumberOfBins = true;
			break;
		case 'Q':
			popUpdatesDuringGOM = 1;
			break;
        case 'w':
            AnalyzeFOS = 1;
            break;
        case 'e':
            writeElitists = 1;
            break;
        case 'E':
            printNewElitists = 1;
            break;
        case 'k':
            dontUseOffspringPopulation = 1;
            break;
        case 's':
            saveEvaluations = 1;
            break;
        case 'h':
            printHelp = 1;
            break;
        case 'a':
            lower_user_range = atof(optarg);
            break;
        case 'b':
            upper_user_range = atof(optarg);
            break;
        case 'f':
            forceNBoundariesUsed = true;
            break;
        case 'r':
            runIndex = atoi(optarg);
            break;
        case 'n':
            basePopulationSize = atoi(optarg);
            break;
        case 'N':
            maximumNumberOfGAMBITs = atoi(optarg);
            break;
        case 'P':
            {
                const string optarg_str = string(optarg);
                if (isNumber(optarg_str))   
                    problemIndex = atoi(optarg);
                else
                {
                    vector<string> tmp;
                    splitString(optarg_str, tmp, '_');
                    
                    if(tmp.size() == 3) 
                    {
                        problemIndex = atoi(tmp[0].c_str());
                        k = atoi(tmp[1].c_str());
                        s = atoi(tmp[2].c_str());
                    }
                    else if(tmp.size() == 2)
                    {
                        problemIndex = atoi(tmp[0].c_str());
                        a_value = atof(tmp[1].c_str());
                    }
                }
            }
            break;
        case 'p':
            useParallelGOM = atoi(optarg);
            break;
        case 'F':
        {
            const string optarg_str = string(optarg);
            if (isNumber(optarg_str)) 
            {   
                FOSIndex = atoi(optarg);
                // Assume that if it is just a number, the FOS is the default one (0, univariate)
                linkage_config = new linkage_config_t();
            } 
            else
            {
                vector<string> tmp;
                splitString(optarg_str, tmp, '_');
                FOSIndex = atoi(tmp[0].c_str());
                
                if(tmp.size() == 2 && FOSIndex == 1) { // MPM linkage
                    size_t block_size = atoi(tmp[1].c_str());
                    linkage_config = new linkage_config_t(block_size);
                } else if(tmp.size() >= 5 && FOSIndex == 2) { // Linkage Tree
                    int similarityMeasure = atoi(tmp[1].c_str());
                    bool filtered = atoi(tmp[2].c_str()) == 1;
                    int maximumSetSize = atoi(tmp[3].c_str());
                    bool is_static = atoi(tmp[4].c_str()) == 1;
                    linkage_config = new linkage_config_t(similarityMeasure, filtered, maximumSetSize, is_static); 
                }
            }
            break;
        }
        case 'm':
            setMethodInitParams(atoi(optarg));
            break;
        case 'M':
            maximumNumberOfGenerations = atoi(optarg);
            break;
        case 'u':
            useNormalizedCVars = true;
            break;
        case 'l':
            logDebugInformation = 1;
            break;
        case 'L':
            if(isNumber(string(optarg)))
            {
                numberOfVariables = atoi(optarg);
                numberOfdVariables = numberOfVariables;
                numberOfcVariables = numberOfVariables;
            } else 
            {
                vector<string> tmp;
                splitString(string(optarg), tmp, '_');
                if(tmp.size() == 2) {
                    numberOfdVariables = atoi(tmp[0].c_str());
                    numberOfcVariables = atoi(tmp[1].c_str());
                    // Default numberOfVariables to be equal to numberOfdVariables
                    numberOfVariables = numberOfdVariables;
                } else if(tmp.size() == 3) {
                    numberOfVariables = atoi(tmp[0].c_str());
                    numberOfdVariables = atoi(tmp[1].c_str());
                    numberOfcVariables = atoi(tmp[2].c_str());
                } 
            }
            break;
        case 'o':
            useOptimalSolution = true;
            break;
        case 'O':
            folder= string(optarg);
            break;
        case 't':
            transformCVariables = true;
            break;
        case 'T':
            maximumNumberOfSeconds = atof(optarg);
            break;
        case 'V':
            vtr = atof(optarg);
            break;
        case 'S':
			{
                fix_seed = true;
				randomSeed = atoll(optarg);
			}
            break;
        case 'I':
            this->problemInstancePath = string(optarg);
            cout << "problemInstancePath: " << this->problemInstancePath << endl;
            break;
        case 'Z':
        {
            if(linkage_config == NULL) {
                cout << "Please specify FOS (F flag) before similarity measure" << endl;
                exit(0);
            }
            linkage_config->lt_similarity_measure = atoi(optarg);
            break;
        }
        case 'G':
            GPUIndex = atoi(optarg);
            break;
        default:
            abort();
    }
  }

  if (printHelp)
  {
    printUsage();
    exit(0);
  }

   
    if(linkage_config == NULL) {
        linkage_config = new linkage_config_t();
    }

    // If problem instance path is passed (and problem index over 1000), we're dealing with a Bayesian Network problem. Update config parameters and parse input data.
    if(this->problemInstancePath != "" && problemIndex >= 1000) {
        this->useBN = true;
        this->alphabetSize = 3; // Discrete variables in solution can be 0, 1 or 2 (A<-/->B, A-->B, A<--B respectively)
        this->maxDiscretizations = maxDiscretizations == -1 ? 9 : maxDiscretizations; // Maximum number of discretizations for continuous variables. Default at 9, only set if not given as commmand line argument.

        // Parse input data (BN structure and data)
        this->data = initializeDataFromPath(true, this->problemInstancePath, runIndex);
        // determine number of d_variables / c_variables based on data + maxDiscretizations
        int numNodes = data->getNumberOfDataColumns();
        this->numberOfdVariables = ((numNodes-1)*numNodes) / 2;
        this->numberOfVariables = this->numberOfdVariables;
        if(this->extraCVarForNumberOfBins)
        {
            this->numberOfcVariables = data->getNumberOfContinuousVariables() * (maxDiscretizations + 1);
        } else {
            this->numberOfcVariables = data->getNumberOfContinuousVariables() * maxDiscretizations;
        }
    }

    fitness = getFitnessClassDiscrete(problemIndex, numberOfVariables);

    setOptimizerName();

    // Override vtr if it was given as command line parameter
    fitness->vtr = this->vtr != 1e+308 ? this->vtr : fitness->vtr; 

  return 1;
}

void Config::printUsage()
{
  cout << "  -h: Prints out this usage information.\n";
  cout << "  -P: Index of optimization problem to be solved (minimization). Default: 0 (oneMax)\n";

  cout << "  List of available problems:\n";
  cout << "    1: OneMax - Sphere (F1)\n";
  cout << "    2: OneMax - Rotated Ellipsoid (F2)\n";
  cout << "    3: Deceptive Trap 5 - Sphere (F3)\n";
  cout << "    4: Deceptive Trap 5 - Rotated Ellipsoid (F4)\n";
  cout << "    5: Blockwise Deceptive Trap 5 - Rotated Ellipsoid (F5) \n";
  cout << "    50: Blockwise Deceptive Trap 5 - Rotated Ellipsoid, but centers at 0 (essentially same as F4)\n";        
  cout << "    51: Blockwise Deceptive Trap 5 - Rotated Ellipsoid, different exponent\n";        
  cout << "    55: Blockwise Deceptive Trap 5 - UNROTATED Ellipsoid\n";        
  cout << "    551: Blockwise Deceptive Trap 5 - UNROTATED Ellipsoid, different exponent\n";        
  cout << "    1000 - 99999999: Bayesian Network structure learning\n";        
  cout << endl;
  
  cout << "  -g: Use guaranteed initialization spread (BNs)\n";
  cout << "  -X: Maximum number of discretizations for continuous variables (BNs). Default: 9\n";
  cout << "  -Y: Use an extra continuous variable per continuous node to indicate the number of bins (BNs)\n";
  cout << "  -Q: Use population updates during GOM\n";
  cout << "  -w: Write FOS to file\n";
  cout << "  -e: Write elitists to file\n";
  cout << "  -E: Print new elitists to console during optimization\n";
  cout << "  -k: Don't use offspring population: Krzysztofs version based on pseudocode in 2014 paper\n";
  cout << "  -s: Save all evaluations\n";
  cout << "  -a: Input lower user range for initialization of continuous variables. Default: 0\n";
  cout << "  -b: Input upper user range for initialization of continuous variables. Default: 1\n";
  cout << "  -f: If -Y is passed, force the amount of boundaries used to be exactly equal to value of extra parameter\n";
  cout << "  -r: Run index. Default: 0\n";
  cout << "  -n: Base population size. Default: 100\n";
  cout << "  -N: Maximum number of simultaneous optimizers for IMS. Default: 10\n";
  cout << "  -P: Problem index. Default: 0\n";
  cout << "  -p: Use parallel GOM. Default: 0\n";
  cout << "  -F: FOS type, 0 - Univariate, 1 - MPM, 2 - Linkage Tree. Default: 2\n";
  cout << "  -m: shortcut for setting method init params, 0 - Method 1, init 1, 1 - Method 1, init 2, 2 - Method 1, init 3, 3 - Method 2, init 1, 4 - Method 2, init 2, 5 - Method 2, init 3, 6 - Method 3, init 1, 7 - Method 3, init 2, 8 - Method 3, init 3. Default: 0\n";
  cout << "  -M: Maximum number of generations. Default: -1\n";
  cout << "  -u: Use normalized continuous variables\n";
  cout << "  -l: Log some debug information\n";
  cout << "  -L: Number of variables. Pass 1 value for equal discrete and continuous variables, 2 (separated by '_') for discrete and continuous different (in that order). Default: 10\n";
  cout << "  -o: Use optimal solution: all solutions will be set to optimal solution, and code will thus converge soon after evaluating. (BNs)\n";
  cout << "  -O: Output folder. Default: \"output\"\n";
  cout << "  -t: Transform continuous variables from [0, 1] to [1, inf) domain.\n";
  cout << "  -T: Maximum number of seconds to run. Default: -1\n";
  cout << "  -V: Value to reach. Default: not set\n";
  cout << "  -S: Random seed. Default: Random timestamp\n";
  cout << "  -I: Path to problem instance (where data for BNs can be found). Default: \"\"\n";
  cout << "  -Z: Similarity Measure to build Linkage Tree, 0 - Mutual Information, 1 - Normalized Mutual Information, 2 - Problem specific. Default: Mutual Information\n";      
  cout << "  -G: NOT USED: GPU index to use for parallel GOM. Default: 0\n";
}

void Config::printOverview()
{
  cout << "### Settings ######################################\n";
  cout << "#\n";
  cout << "# Use partial evaluations : " << (usePartialEvaluations ? "enabled" : "disabled")  << endl;
  cout << "# Write FOS to file : " << (AnalyzeFOS ? "enabled" : "disabled") << endl;
  cout << "# Save all evaluations : " << (saveEvaluations ? "enabled" : "disabled") << endl;
  if(linkage_config != NULL) 
  {
    cout << "# similarity measure : " << (linkage_config->lt_similarity_measure ? "normalized MI" : "MI") << endl;
  }
  
  
  cout << "#\n";
  cout << "###################################################\n";
  cout << "#\n";
  cout << "# Problem                      = " << fitness->name << endl;
  cout << "# Problem Instance Filename    = " << this->problemInstancePath << endl;
  cout << "# FOS                          = " << FOSName << endl;
  cout << "# Number of variables          = " << numberOfVariables << endl;
  cout << "# Time Limit (seconds)         = " << maximumNumberOfSeconds << endl;
  cout << "# VTR                          = " << ((vtr < 1e+308) ? to_string(vtr) : "not set") << endl;
  cout << "# Random seed                  = " << randomSeed << endl;
  cout << "# Folder                       = " << folder << endl;
  cout << "#\n";
  cout << "### Settings ######################################\n";
}

void Config::checkOptions()
{
}

}}