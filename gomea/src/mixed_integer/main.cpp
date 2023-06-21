#include "gomea/src/mixed_integer/Config.hpp"
#include "gomea/src/mixed_integer/gomea.hpp"
#include "gomea/src/mixed_integer/gomeaIMS.hpp"
#include "gomea/src/mixed_integer/simpleGAMBIT.hpp"

using namespace gomea;
int main(int argc, char **argv)
{
    mixedinteger::Config *config = new mixedinteger::Config();
    config->parseCommandLine(argc, argv);
    config->checkOptions();
    config->printOverview();
    // mixedinteger::GOMEA *gomeaInstance = new mixedinteger::gomeaIMS(config);
    mixedinteger::simpleGAMBIT *simpleGAMBITInstance = new mixedinteger::simpleGAMBIT(config);

    try
    {
        // gomeaInstance->run();
        cout << "[DEBUGGING] Starting the GAMBIT run!" << endl;
        simpleGAMBITInstance->run();
    }
    catch (utils::customException &ex)
    {
    }

    // delete gomeaInstance;
    delete simpleGAMBITInstance;
    delete config;
    
    return 0;
}

// main()
// {
    // mixedinteger::Config *config = new mixedinteger::Config();
    // config->parseCommandLine(argc, argv);
    // config->checkOptions();
    // config->printOverview();
    // mixedinteger::simpleGAMBIT *simpleGAMBITInstance = new mixedinteger::simpleGAMBIT(config);
// }
