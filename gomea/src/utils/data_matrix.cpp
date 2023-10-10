//
// Created by Damy Ha on 21-Dec-22.
//

#include "gomea/src/utils/data_matrix.h"

// Instantiation
template class DataMatrix<double>;
template class DataMatrix<size_t>;

/**
 * Constructor
 */
template<typename T>
DataMatrix<T>::DataMatrix() = default;

/**
 * Constructor
 * @param matrix A predefined matrix
 */
template<typename T>
DataMatrix<T>::DataMatrix(const vector<vector<T>> &matrix) : matrix(matrix) {
    // Check that the given matrix is valid
    this->checkMatrixHasEqualColumns();
}

/**
 * Constructor that reserves space
 * @param numberOfRows The number of rows to reserve
 * @param numberOfColumns The number of columns to reserve
 */
template<typename T>
DataMatrix<T>::DataMatrix(size_t numberOfRows, size_t numberOfColumns) {
    vector<vector<T>> result(numberOfRows, vector<T>(numberOfColumns));
    this->matrix = result;
}

/**
 * Constructor that reserves space and initializes all cells to a value
 * @param numberOfRows The number of rows to reserve
 * @param numberOfColumns The number of columns to reserve
 * @param value The value to set each cell to
 */
template<typename T>
DataMatrix<T>::DataMatrix(size_t numberOfRows, size_t numberOfColumns, T value) {
    vector<vector<T>> result(numberOfRows, vector<T>(numberOfColumns, value));
    this->matrix = result;
}

/**
 * Reads in a file of numbers and transform it into a 2D matrix
 * @param filename The name of the file.
 */
template<typename T>
void DataMatrix<T>::readFromFile(const string &filename) {
    // Clear the matrix
    this->matrix.clear();

    // Open the input file
    ifstream input(filename);

    // Read the file line by line
    string line;
    while (getline(input, line)) {
        // Split the line by space and store the numbers in a temporary vector
        istringstream iss(line);
        vector<T> rowOfValues;
        T number;
        while (iss >> number) {
            rowOfValues.push_back(number);
        }

        // Add the row to the results
        this->matrix.push_back(rowOfValues);
    }

    // Perform check that matrix is equal sized
    this->checkMatrixHasEqualColumns();

    // Close the input file
    input.close();
}

// Getters
template<typename T>
vector<vector<T>> DataMatrix<T>::getRawMatrix() const {
    return this->matrix;
}

template<typename T>
vector<T> DataMatrix<T>::getRow(size_t index) const {
    // Check if the index is within the bounds
    assert( index < this->matrix.size());

    return this->matrix[index];
}

template<typename T>
vector<T> DataMatrix<T>::getColumn(size_t index) const {
    // Check if the index is within the bounds of the row
    size_t numberOfColumns = this->matrix[0].size();
    size_t numberOfSamples = this->matrix.size();
    assert( index < numberOfColumns);

    // Create a vector to store the column
    vector<T> column(numberOfSamples);

    // Iterate over the rows in the 2D matrix
    for (size_t i = 0; i < numberOfSamples; ++i) {
        column[i] = this->matrix[i][index];
    }

    return column;
}

template<typename T>
DataMatrix<T> DataMatrix<T>::getColumns(vector<size_t> indices) const {
    // Retrieve the data
    vector<vector<T>> data(this->matrix.size(), vector<T>(indices.size()));
    for(size_t localNodeIndex = 0; localNodeIndex < indices.size(); ++localNodeIndex) {
        // Check that the indices to be retrieved are smaller than the max number of columns
        size_t nodeIndex = indices[localNodeIndex];
        assert(nodeIndex < matrix[0].size());

        // Set the relevant columns and rows
        for(size_t sampleIndex = 0; sampleIndex < this->matrix.size(); ++sampleIndex) {
            data[sampleIndex][localNodeIndex] = this->matrix[sampleIndex][nodeIndex];
        }
    }

    // Create another matrix
    DataMatrix<T> result(data);
    return result;
}

template<typename T>
T DataMatrix<T>::getElement(size_t i, size_t j) const {
    // Perform checks
    assert( i < this->matrix.size());
    assert( j < this->matrix[i].size());

    // Set the value
    return this->matrix[i][j];
}

template<typename T>
size_t DataMatrix<T>::getNumberOfRows() const { return matrix.size(); }

template<typename T>
size_t DataMatrix<T>::getNumberOfColumns() const { return this->matrix[0].size(); }

// Setters
template<typename T>
void DataMatrix<T>::setRow(size_t index, const vector<T> &row) {
    // Perform check
    assert (index < this->matrix.size());

    // Set the row
    matrix[index] = row;
}

template<typename T>
void DataMatrix<T>::setColumn(size_t index, const vector<T> &column) {
    // Perform check
    assert (index < this->matrix[0].size());

    // Iterate over the rows in the 2D matrix
    for (size_t i = 0; i < matrix.size(); i++) {
        matrix[i][index] = column[i];
    }
}

template<typename T>
void DataMatrix<T>::setElement(size_t i, size_t j, T value) {
    // Perform checks
    assert( i < this->matrix.size());
    assert( j < this->matrix[i].size());

    // Set the value
    this->matrix[i][j] = value;
}

/**
 * Overloaded operator that adds to matrices
 * @param other The other matrix
 * @return
 */
template<typename T>
DataMatrix<T> DataMatrix<T>::operator+(const DataMatrix<T> &other) const {
    // Determine variables
    size_t numberOfRows = matrix.size();
    size_t numberOfColumns = matrix[0].size();

    // Check if the matrices have the same dimensions
    assert(matrix.size() == other.matrix.size());
    assert(matrix[0].size() == other.matrix[0].size());

    // Iterate over the elements of the matrices and add them
    vector<vector<T>> data = this->matrix;
    for (size_t i = 0; i < numberOfRows; i++) {
        for (size_t j = 0; j < numberOfColumns; j++) {
            data[i][j] += other.matrix[i][j];
        }
    }

    // Create another matrix
    DataMatrix<T> result(data);
    return result;
}

/**
 * Overloaded operator that adds the values from the other matrix
 * @param other The other matrix to add
 * @return A pointer to this object
 */
template<typename T>
DataMatrix<T> DataMatrix<T>::operator+=(const DataMatrix<T> &other) {
    // Determine variables
    size_t numberOfRows = matrix.size();
    size_t numberOfColumns = matrix[0].size();

    // Check if the matrices have the same dimensions
    assert(matrix.size() == other.matrix.size());
    assert(matrix[0].size() == other.matrix[0].size());

    // Iterate over the elements of the matrices and add them
    for (size_t i = 0; i < numberOfRows; i++) {
        for (size_t j = 0; j < numberOfColumns; j++) {
            this->matrix[i][j] += other.matrix[i][j];
        }
    }

    return *this;
}

//////////////
/// Checks ///
//////////////
/**
 * Checks if the matrix has an equal number of elements in each row
 */
template<typename T>
void DataMatrix<T>::checkMatrixHasEqualColumns() {
    size_t rowSize = this->matrix[0].size();
    for (const auto &row : this->matrix) {
        assert(row.size() == rowSize && "Each row must have the same number of columns");
    }
}

//////////////////
/// Conversion ///
//////////////////
/**
 * Converts a data matrix of type T to type U
 * @tparam T The original data type
 * @tparam U The target data type
 * @param dataMatrix The data matrix to convert
 * @return The converted data matrix
 */
template<typename T, typename U>
DataMatrix<U> convertToType(const DataMatrix<T>& dataMatrix) {
    // Initialize variables
    size_t numRows = dataMatrix.getNumberOfDataRows();
    size_t numCols = dataMatrix.getNumberOfDataColumns();

    // Create a new DataStructure<U> with the same dimensions
    DataMatrix<U> result(numRows, numCols);

    // Iterate over each element and convert from type T to type U
    for (size_t i = 0; i < numRows; ++i) {
        for (size_t j = 0; j < numCols; ++j) {
            // Convert and assign each element
            U convertedValue = static_cast<U>(dataMatrix.getElement(i, j));
            result.setElement(i, j, convertedValue);
        }
    }

    return result;
}

template<typename T, typename U>
vector<U> convertToType(const vector<T>& dataColumn) {
    // Initialize variables
    size_t numberOfSamples = dataColumn.size();
    vector<U> result(numberOfSamples);

    // Iterate over each element and convert from type T to type U
    for (size_t i = 0; i < numberOfSamples; ++i) {
        result[i] = static_cast<U>(dataColumn[i]);
    }

    return result;
}

// Explicitly instantiate the functions
template vector<size_t> convertToType<int, size_t>(const vector<int>& dataColumn);
template vector<size_t> convertToType<double, size_t>(const vector<double>& dataColumn);