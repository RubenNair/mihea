//
// Created by Damy Ha on 21-Dec-22.
// Temporary class to circumvent Armadillo
//

#ifndef IMPLEMENTATIONS_DATA_MATRIX_H
#define IMPLEMENTATIONS_DATA_MATRIX_H

#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <cassert>

using namespace std;

template<typename T>
class DataMatrix {
public:
    // Constructor
    DataMatrix();
    DataMatrix(const vector<vector<T>> &matrix);
    DataMatrix(size_t numberOfRows, size_t numberOfColumns);
    DataMatrix(size_t numberOfRows, size_t numberOfColumns, T value);

    // Data Processing method
    void readFromFile(const string &filename);

    DataMatrix operator+(const DataMatrix &other) const;
    DataMatrix operator+=(const DataMatrix &other);

    // Counting

    // Getters
    vector<vector<T>> getRawMatrix() const;
    vector<T> getRow(size_t index) const;
    vector<T> getColumn(size_t index) const;
    DataMatrix getColumns(vector<size_t> indices) const;
    T getElement(size_t i, size_t j) const;
    size_t getNumberOfRows() const;
    size_t getNumberOfColumns() const;

    // Setters
    void setRow(size_t index, const vector<T> &row);
    void setColumn(size_t index, const vector<T> &column);
    void setElement(size_t i, size_t j, T value);

private:
    // Variables
    vector<vector<T>> matrix;

    // Check
    void checkMatrixHasEqualColumns();  // Checks if the matrix has an equal number of elements in each row
};

/// Conversion
template<typename T, typename U> DataMatrix<U> convertToType(const DataMatrix<T>& dataMatrix);
template<typename T, typename U> vector<U> convertToType(const vector<T>& dataColumn);

#endif //IMPLEMENTATIONS_DATA_MATRIX_H
