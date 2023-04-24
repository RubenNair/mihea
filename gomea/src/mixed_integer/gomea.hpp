#pragma once

#include <string> 
#include <iostream>
using namespace std;

#include "gomea/src/mixed_integer/Config.hpp"

namespace gomea{
namespace mixedinteger{

class GOMEA
{
public:
    virtual void run() = 0;
    virtual ~GOMEA(){};

    double readVTR(Config *config);
};

}}