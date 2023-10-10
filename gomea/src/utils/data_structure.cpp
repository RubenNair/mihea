//
// Created by Damy Ha on 28-Sep-22.
//

#include <utility>

#include "gomea/src/utils/data_structure.h"

// Instantiation
template class DataStructure<double>;
template class DataStructure<size_t>;

/**
 * Create a data structure object
 * @param path_to_data The path to retrieve the data
 * @param path_to_info The path to retrieve information about the data
 */
template<typename T>
DataStructure<T>::DataStructure(const string& path_to_data, const string& path_to_info) {
    // Check if path to data exists
    if (!path_exists(path_to_data)) {
        // Generate error message
        std::string error_message = std::string("Path to:") + path_to_data +
                std::string(" does not exist. Please run generate_data.py to generate data "
                            "(and put the file in the data folder).\n");
        throw std::runtime_error(error_message);
    }

    // Check if path to info of data exists
    if (!path_exists(path_to_info)) {
        // Generate error message
        std::string error_message = std::string("Path to:") + path_to_info +
                std::string(" does not exist. Please run generate_data.py to generate info of the data "
                            "(and put the file in the data folder).\n");
        throw std::runtime_error(error_message);
    }


    // Retrieve the file path
    this->path_data = path_to_data;
    this->path_info = path_to_info;

    // Initialize variables
//    this->column_names.resize(0);
//    this->column_type.resize(0);
//    this->column_number_of_classes.resize(0);

    // Retrieve the data
    this->load_data_from_file_v1(path_to_data);
    this->load_problem_data_from_file_v1(path_to_info);

    // Process variables
    this->preProcessNumberOfUniqueValues();

    // Perform checks
    this->check_data_and_info_equal_columns();

}

/**
 * Create a data structure object manually. (Usually only used for debugging or testing)
 */
template<typename T>
DataStructure<T>::DataStructure(const vector<vector<T>> &data,
                             const vector<string> &columns_names,
                             const vector<ColumnDataType> &column_type,
                             const vector<size_t> &number_of_classes) {
    // Initialize variables
    this->path_data.clear();
    this->path_info.clear();

    // Copy over variables
    this->column_names = columns_names;
    this->column_type = column_type;
    this->column_number_of_classes = number_of_classes;

    // Initialize the data matrix
    this->data = DataMatrix<T>(data);

    // Process variables
    this->preProcessNumberOfUniqueValues();
}

/**
 * Create a data structure object manually. (Usually only used for debugging or testing)
 */
template<typename T>
DataStructure<T>::DataStructure(DataMatrix<T> data,
                                const vector<string> &columns_names,
                                const vector<ColumnDataType> &column_type,
                                const vector<size_t> &number_of_classes,
                                const vector<size_t> &number_of_unique_values) {
    // Initialize variables
    this->path_data.clear();
    this->path_info.clear();

    // Copy over variables
    this->column_names = columns_names;
    this->column_type = column_type;
    this->column_number_of_classes = number_of_classes;
    this->number_of_unique_values = number_of_unique_values;

    // Copy over the data
    this->data = std::move(data);

    // Perform checks
    this->check_discrete_nodes_have_number_of_variables();
}

template<typename T>
DataStructure<T>::~DataStructure() = default;

/**
 * Make a clone of the DataStructure
 * @return A clone of the data structure
 */
template<typename T>
shared_ptr<DataStructure<T>> DataStructure<T>::clone() {
    // Copy the data
    DataMatrix<T> copyOfData = this->data;

    // Create a new object
    shared_ptr<DataStructure> result = make_shared<DataStructure>(copyOfData,
                                                                  this->column_names,
                                                                  this->column_type,
                                                                  this->column_number_of_classes,
                                                                  this->number_of_unique_values);

    return result;
}

/**
 * Load data from the data file.
 * The data file should only contains columns of data (i.e. no columns names etc.), seperated by spaces.
 * All data elements should be discretezed starting from 0, 1, 2, .... At the end, the results should be translated back.
 * This is done, as 0, 1, 2, ...., makes the computation faster
 * Example:
 * 1 0 1 0
 * 3 0 2 1
 * 2 1 3 0
 * @param path_to_data The path to the data file
 */
template<typename T>
void DataStructure<T>::load_data_from_file_v1(const string& path_to_data) {
    // Load the data
    this->data.readFromFile(path_to_data);
}

/**
 * Loads information about the associated data file.
 * Information includes: column names, column data type and number of classes per column. Example:
 * Asia                 - Discrete - 2
 * Tuberculosis         - Discrete - 2
 * Smoker               - Discrete - 2
 * @param path_to_info The path to the info file of the data
 */
template<typename T>
void DataStructure<T>::load_problem_data_from_file_v1(const string& path_to_info){
    // Initialize variables
    string myText;

    // Read from the text file
    ifstream MyReadFile(path_to_info);
    while (getline (MyReadFile, myText)) {
        // Check if the line is empty
        if (myText.empty()) continue;

        // Split the line to retrieve data
        vector<string> information = split_string(myText, "-");

        // Process the data of a row
        string column_name;
        ColumnDataType data_type;
        size_t number_of_classes;
        process_info_v1(information,
                        column_name,
                        data_type,
                        number_of_classes,
                        path_to_info);

        // Assign the data
        this->column_names.push_back(column_name);
        this->column_type.push_back(data_type);
        this->column_number_of_classes.push_back(number_of_classes);
    }

    // Close the file
    MyReadFile.close();
}

/**
 * Process information of an info line.
 * Example of a line:
 * Asia                 - Discrete - 2
 * @param information A vector containing strings of information
 * @param column_name The column name
 * @param column_data_type enum of discrete or continuous data
 * @param column_number_of_instantiations Number of classes (read discretizations) of the data
 * @param path_to_info Path to the info file (used for warning messages)
 */
template<typename T>
void DataStructure<T>::process_info_v1(vector<string> information,
                                       string &column_name,
                                       ColumnDataType &column_data_type,
                                       size_t &column_number_of_instantiations,
                                       const string& path_to_info) {
    // Process the column name
    column_name = strip(information[0]);

    // Process the column data type
    string column_data_type_str = strip(information[1]); // Remove spaces
    if (column_data_type_str == string("Discrete")) {
        // Discrete
        column_data_type = Discrete;
    } else if (column_data_type_str == string("Continuous")) {
        // Continuous
        column_data_type = Continuous;
    } else {
        // Generate error message
        std::string error_message = string("In:") + path_to_info +
                std::string(", an unknown data type (Discrete/Continuous) has been detected:'" + information[1] + "'\n");
        throw std::runtime_error(error_message);
    }

    // Process the number of classes
    column_number_of_instantiations = stoi(information[2]);
}

/**
 * Calculates the number of unique values in the data (per node) and stores it in 'number_of_unique_values'
 */
template<typename T>
void DataStructure<T>::preProcessNumberOfUniqueValues() {

    // Calculate the number of unique values per node
    size_t numberOfNodes = this->getNumberOfDataColumns();
    vector<size_t> result(numberOfNodes);
    for (size_t node_index = 0; node_index < numberOfNodes; ++node_index) {
        // Retrieve the node type
        ColumnDataType node_type = this->column_type[node_index];
        if (node_type == Continuous) {
            // Get the unique values of the data
            vector<T> dataOfNode = this->data.getColumn(node_index);
            sort(dataOfNode.begin(), dataOfNode.end());
            dataOfNode.erase( unique( dataOfNode.begin(), dataOfNode.end() ), dataOfNode.end() );

            // Assign the result
            result[node_index] = dataOfNode.size();
        } else {
            // The number of unique values of discrete nodes is already given
            result[node_index] = this->column_number_of_classes[node_index];
        }
    }

    this->number_of_unique_values = result;

}

/**
 * Checks that the data and info file have an equal number of columns
 */
template<typename T>
void DataStructure<T>::check_data_and_info_equal_columns() {
    // Determine variables
    size_t number_of_columns_data = this->data.getNumberOfColumns();
    assert(number_of_columns_data == this->column_names.size());
    assert(number_of_columns_data == this->column_type.size());
    assert(number_of_columns_data == this->column_number_of_classes.size());
}

/**
 * Checks that the discrete classes have been parsed the number of classes it can take
 */
template<typename T>
void DataStructure<T>::check_discrete_nodes_have_number_of_variables() {
    // Perform general checks
    this->check_data_and_info_equal_columns();

    // Go over each column
    for (size_t node_index = 0; node_index < this->getNumberOfDataColumns(); ++node_index) {
        // Check for discrete nodes only
        ColumnDataType nodeType = this->column_type[node_index];
        if (nodeType == Discrete) {
            // Retrieve the number of discretizations
            size_t number_of_discretizations = this->getColumnNumberOfClasses()[node_index];
            assert(number_of_discretizations != 0);
        }
    }
}

template<typename T>
vector<size_t> DataStructure<T>::determineDiscreteNodeIndices() const {
    // Initialize variables
    vector<size_t> result;
    result.reserve(this->column_type.size());

    // Add the nodes that need to be discretized
    for (size_t node_index = 0; node_index < this->column_type.size(); ++node_index) {
        if (this->column_type[node_index] == Discrete) {
            result.push_back(node_index);
        }
    }

    return result;
}

template<typename T>
vector<size_t> DataStructure<T>::determineContinuousNodeIndices() const {
    // Initialize variables
    vector<size_t> result;
    result.reserve(this->column_type.size());

    // Add the nodes that need to be discretized
    for (size_t node_index = 0; node_index < this->column_type.size(); ++node_index) {
        if (this->column_type[node_index] == Continuous) {
            result.push_back(node_index);
        }
    }

    return result;

}

// Getters
template<typename T> const string &DataStructure<T>::getPathData() const { return path_data; }
template<typename T> const string &DataStructure<T>::getPathInfo() const { return path_info; }
template<typename T> const DataMatrix<T> &DataStructure<T>::getDataMatrix() const { return data; }
template<typename T> const vector<string> &DataStructure<T>::getColumnNames() const { return column_names; }
template<typename T> const vector<ColumnDataType> &DataStructure<T>::getColumnType() const { return column_type; }
template<typename T> const vector<size_t> &DataStructure<T>::getColumnNumberOfClasses() const { return column_number_of_classes; }
template<typename T> int DataStructure<T>::getNumberOfDataRows() { return this->data.getNumberOfRows(); }
template<typename T> int DataStructure<T>::getNumberOfDataColumns() { return this->data.getNumberOfColumns(); }
template<typename T> const vector<size_t> &DataStructure<T>::getNumberOfUniqueValues() const { return number_of_unique_values; }
template<typename T> size_t DataStructure<T>::getNumberOfContinuousVariables() { return std::count(this->column_type.begin(), this->column_type.end(),Continuous); }



