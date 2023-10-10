//
// Created by Damy Ha on 01-Oct-22.
//

#ifndef IMPLEMENTATIONS_DATA_STRUCTURE_H
#define IMPLEMENTATIONS_DATA_STRUCTURE_H

#include <fstream>
#include <iostream>
#include <memory>
#include <sys/stat.h>

#include "gomea/src/utils/data_matrix.h"
#include "gomea/src/mixed_integer/BN_utils.h"


using namespace std;

enum ColumnDataType { Discrete, Continuous };

template<typename T>
class DataStructure {
public:
    DataStructure(const string &path_to_file,
                  const string &path_to_info);                          // Construct data from file
    DataStructure(const vector<vector<T>> &data,
                  const vector<string> &columns_names,
                  const vector<ColumnDataType> &column_type,
                  const vector<size_t> &number_of_classes);             // Construct data manually
    DataStructure(DataMatrix<T> data,
                  const vector<string> &columns_names,
                  const vector<ColumnDataType> &column_type,
                  const vector<size_t> &number_of_classes,
                  const vector<size_t> &number_of_unique_values);       // Construct data manually
    ~DataStructure();

    // Clone
    shared_ptr<DataStructure> clone();      // Clone this object

    vector<size_t> determineDiscreteNodeIndices() const;
    vector<size_t> determineContinuousNodeIndices() const;

    // Getters
    const string &getPathData() const;
    const string &getPathInfo() const;
    const DataMatrix<T> &getDataMatrix() const;
    const vector<string> &getColumnNames() const;
    const vector<ColumnDataType> &getColumnType() const;
    const vector<size_t> &getColumnNumberOfClasses() const;
    int getNumberOfDataRows();
    int getNumberOfDataColumns();
    size_t getNumberOfContinuousVariables();
    const vector<size_t> &getNumberOfUniqueValues() const;

private:
    // Data file
    string path_data;   // Path to the data
    string path_info;   // Path to information of the data

    // Basic information
    DataMatrix<T> data;                             // Contains the data
    vector<string> column_names;                    // Names of the columns
    vector<ColumnDataType> column_type;             // Discrete or Continuous data per column
    vector<size_t> column_number_of_classes;        // The number of classes a random variable can take

    // Derived information
    vector<size_t> number_of_unique_values;         // Used by the fitness functions (e.g. mdl). It is precalculated here over the fitness function to avoid recalculation when invoked.

    // Process and load functions
    void process_info_v1(vector<string> information,
                         string &column_name,
                         ColumnDataType &column_data_type,
                         size_t &column_number_of_instantiations,
                         const string& path_to_info);      // Processes a line of an information file

    void load_data_from_file_v1(const string& path_to_data);          // Loads data from a file into this->data
    void load_problem_data_from_file_v1(const string& path_to_info);  // Loads information about the problem

    // Preprocess functions
    void preProcessNumberOfUniqueValues();

    // Check functions
    void check_data_and_info_equal_columns();               // Performs a simple check to see if data.txt and info.txt are valid
    void check_discrete_nodes_have_number_of_variables();   // Performs a simple check if the discrete nodes have a valid number of discretizations

};

#endif //IMPLEMENTATIONS_DATA_STRUCTURE_H
