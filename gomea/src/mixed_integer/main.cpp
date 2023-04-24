#include "gomea/src/mixed_integer/Config.hpp"
#include "gomea/src/mixed_integer/gomea.hpp"
#include "gomea/src/mixed_integer/gomeaIMS.hpp"

using namespace gomea;
int main(int argc, char **argv)
{
    mixedinteger::Config *config = new mixedinteger::Config();
    config->parseCommandLine(argc, argv);
    config->checkOptions();
    config->printOverview();
    mixedinteger::GOMEA *gomeaInstance = new mixedinteger::gomeaIMS(config);

    try
    {
        gomeaInstance->run();
    }
    catch (utils::customException &ex)
    {
    }

    delete gomeaInstance;
    delete config;
    
    return 0;
}
