#pragma once

#include <string> 
#include <iostream>


#include "gomea/src/common/solution_mixed.hpp"

namespace gomea{
namespace mixedinteger{

class Solutionset
{
public:
    Solutionset(std::vector<solution_mixed*> solutions);
    Solutionset(size_t size) : solutions(size) {}
    virtual ~Solutionset() = default;

    solution_mixed* operator[](std::size_t idx)       { return solutions[idx]; }
    void resize(std::size_t size) { solutions.resize(size); }
    std::size_t size() const { return solutions.size(); }


    std::vector<solution_mixed*> solutions;

    
};

}}