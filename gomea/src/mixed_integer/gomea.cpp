#include "gomea/src/mixed_integer/gomea.hpp"

namespace gomea{
namespace mixedinteger{

double GOMEA::readVTR(Config *config)
{
    string filename = config->folder + "/vtr.txt";
    ifstream inFile;
    inFile.open(filename);

    string vtr_str;
    inFile >> vtr_str;
    double vtr = stod(vtr_str);

    inFile.close();

    return vtr;
}

}}