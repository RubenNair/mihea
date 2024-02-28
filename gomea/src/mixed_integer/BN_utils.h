//
// Created by Damy Ha on 27-Sep-22.
//
#define _USE_MATH_DEFINES
#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iterator>
#include <memory>
#include <math.h>
#include <string>
#include <sys/stat.h>
#include <vector>


using namespace std;

// Global
bool path_exists(string path_to_file);                                  // Checks if a path exists

vector<string> split_string(string string_to_split, string delimiter);  // Splits a string
string strip(const string &input);                                      // Strips white spaces


// Math
double logfactorial(size_t number);         // The log of the factorial
double naiveLogFactorial(size_t number);    // The true log factorial
double log_factorial_ramanujan(size_t m);   // Ramanujan's approximation of the log factorial
double lutLogFactorial(size_t number);      // A look up table for log factorials

// Network
size_t calculateNumberOfLinks(size_t number_of_nodes);      // Calculates the number of archs in a network