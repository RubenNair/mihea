#include "gomea/src/discrete/Config.hpp"
#include "gomea/src/discrete/gomea.hpp"
#include "gomea/src/discrete/gomeaIMS.hpp"

// namespace gomea{
// namespace discrete{

int main(int argc, char **argv)
{
    gomea::discrete::Config *config = new gomea::discrete::Config();
    config->parseCommandLine(argc, argv);
    config->checkOptions();
    config->printOverview();

    gomea::discrete::GOMEA *gomeaInstance = new gomea::discrete::gomeaIMS(config);

    try
    {
        gomeaInstance->run();
    }
    catch (gomea::utils::customException &ex)
    {}

    delete gomeaInstance;
    delete config;

    return 0;
}

// }}