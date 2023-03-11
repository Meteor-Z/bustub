//===----------------------------------------------------------------------===//
//
//                         BusTub
//
// p0_starter.h
//
// Identification: src/include/primer/p0_starter.h
//
// Copyright (c) 2015-2020, Carnegie Mellon University Database Group
//
//===----------------------------------------------------------------------===//

#pragma once

#include <memory>
#include <stdexcept>
#include <vector>

#include "common/exception.h"

namespace bustub
{

/**
 * The Matrix type defines a common
 * interface for matrix operations.
 */
    template<typename T>
    class Matrix
    {
    protected:
        /**
         * TODO(P0): Add implementation ok
         * 增加一个构造的实例
         * Construct a new Matrix instance.
         * @param rows The number of rows
         * @param cols The number of columns
         *
         */
        Matrix(int rows, int cols) : rows_(rows), cols_(cols)
        {
            linear_ = new T[rows * cols];
        }

        /** The number of rows in the matrix */
        int rows_;
        /** The number of columns in the matrix */
        int cols_;

        /**
         * 在构造函数中创建数组，在析构函数中释放数组
         * TODO(P0): Allocate the array in the constructor.  ok
         * TODO(P0): Deallocate the array in the destructor. ok
         * A flattened array containing the elements of the matrix.
         */
        T* linear_;

    public:
        /** @return The number of rows in the matrix */
        virtual int GetRowCount() const = 0;

        /** @return The number of columns in the matrix */
        virtual int GetColumnCount() const = 0;

        /**
         * Get the (i,j)th matrix element.
         *
         * Throw OUT_OF_RANGE if either index is out of range.
         *
         * @param i The row index
         * @param j The column index
         * @return The (i,j)th matrix element
         * @throws OUT_OF_RANGE if either index is out of range
         */
        virtual T GetElement(int i, int j) const = 0;

        /**
         * Set the (i,j)th matrix element.
         *
         * Throw OUT_OF_RANGE if either index is out of range.
         *
         * @param i The row index
         * @param j The column index
         * @param val The value to insert
         * @throws OUT_OF_RANGE if either index is out of range
         */
        virtual void SetElement(int i, int j, T val) = 0;

        /**
         * Fill the elements of the matrix from `source`.
         *
         * Throw OUT_OF_RANGE in the event that `source`
         * does not contain the required number of elements.
         *
         * @param source The source container
         * @throws OUT_OF_RANGE if `source` is incorrect size
         */
        virtual void FillFrom(const std::vector <T>& source) = 0;

        /**
         * Destroy a matrix instance.
         * TODO(P0): Add implementation
         */
        virtual ~Matrix()
        {
            delete[] linear_;
        }
    };

/**
 * The RowMatrix type is a concrete matrix implementation.
 * It implements the interface defined by the Matrix type.
 */
    template<typename T>
    class RowMatrix : public Matrix<T>
    {
    public:
        /**
         * TODO(P0): Add implementation
         *
         * Construct a new RowMatrix instance.
         * @param rows The number of rows
         * @param cols The number of columns
         */
        RowMatrix(int rows, int cols) : Matrix<T>(rows, cols)
        {
            this->data_ = new T* [rows]; // 二级指针指向一级指针，然后一级指针是指向就是 T
            for (int i = 0; i < rows; i++)
            {
                this->data_[i] = new T[cols];
            }
        }

        /**
         * TODO(P0): Add implementation  ok
         * @return The number of rows in the matrix
         */
        int GetRowCount() const override
        {
            return this->rows_;
        }

        /**
         * TODO(P0): Add implementation  ok
         * @return The number of columns in the matrix
         */
        int GetColumnCount() const override
        {
            return this->cols_;
        }

        /**
         * TODO(P0): Add implementation  ok
         * 这个是返回i,j 这个值
         * 如果出界了，那么就返回一个异常
         * Get the (i,j)th matrix element.
         *
         * Throw OUT_OF_RANGE if either index is out of range.
         *
         * @param i The row index
         * @param j The column index
         * @return The (i,j)th matrix element
         * @throws OUT_OF_RANGE if either index is out of range
         */
        T GetElement(int i, int j) const override
        {
            if (0 <= i && i < GetRowCount() && 0 <= j && j < GetColumnCount())
            {
                return this->data_[i][j];
            } else
            {
                throw NotImplementedException{"RowMatrix::GetElement() not implemented."};
            }
        }

        /**
         * Set the (i,j)th matrix element. ok
         * 如果出界了，那么就抛出异常
         * Throw OUT_OF_RANGE if either index is out of range.
         *
         * @param i The row index
         * @param j The column index
         * @param val The value to insert
         * @throws OUT_OF_RANGE if either index is out of range
         */
        void SetElement(int i, int j, T val) override
        {
            if (0 <= i && i < GetRowCount() && 0 <= j && j <= GetColumnCount())
            {
                this->data_[i][j] = val;
            } else
            {
                throw Exception(ExceptionType::OUT_OF_RANGE, "出界了");
            }
        }

        /**
         * TODO(P0): Add implementation  ok
         * 将这个矩阵进行填充，如果不能填充，要返回异常
         * Fill the elements of the matrix from `source`.
         *
         * Throw OUT_OF_RANGE in the event that `source`
         * does not contain the required number of elements.
         *
         * @param source The source container
         * @throws OUT_OF_RANGE if `source` is incorrect size
         */
        void FillFrom(const std::vector <T>& source) override
        {
            if (GetColumnCount() * GetRowCount() != source.size())
            {
                throw NotImplementedException{"RowMatrix::FillFrom() not implemented."};
            }
            // 通过循环赋值
            int cnt = 0;
            for (int i = 0; i < GetColumnCount(); i++)
            {
                for (int j = 0; j < GetRowCount(); j++)
                {
                    this->data_[i][j] = source[cnt++];
                }
            }
        }

        /**
         *
         * TODO(P0): Add implementation  ok
         * 将元素进行释放
         * Destroy a RowMatrix instance.
         */
        ~RowMatrix() override
        {
            delete[] data_;
        }

    private:
        /**
         * A 2D array containing the elements of the matrix in row-major format.
         *
         * TODO(P0):
         * - Allocate the array of row pointers in the constructor.
         * - Use these pointers to point to corresponding elements of the `linear` array.
         * - Don't forget to deallocate the array in the destructor.
         */
        T** data_;
    };

/**
 * The RowMatrixOperations class defines operations
 * that may be performed on instances of `RowMatrix`.
 */
    template<typename T>
    class RowMatrixOperations
    {
    public:
        /**
         * Compute (`matrixA` + `matrixB`) and return the result.
         * Return `nullptr` if dimensions mismatch for input matrices.
         * @param matrixA Input matrix
         * @param matrixB Input matrix
         * @return The result of matrix addition
         */
        static std::unique_ptr <RowMatrix<T>> Add(const RowMatrix<T>* matrixA, const RowMatrix<T>* matrixB)
        {

            // TODO(P0): Add implementation
            // 计算

            return std::unique_ptr<RowMatrix<T>>(nullptr);
        }

        /**
         * Compute the matrix multiplication (`matrixA` * `matrixB` and return the result.
         * Return `nullptr` if dimensions mismatch for input matrices.
         * @param matrixA Input matrix
         * @param matrixB Input matrix
         * @return The result of matrix multiplication
         */
        static std::unique_ptr <RowMatrix<T>> Multiply(const RowMatrix<T>* matrixA, const RowMatrix<T>* matrixB)
        {
            // TODO(P0): Add implementation
            return std::unique_ptr<RowMatrix<T>>(nullptr);
        }

        /**
         * Simplified General Matrix Multiply operation. Compute (`matrixA` * `matrixB` + `matrixC`).
         * Return `nullptr` if dimensions mismatch for input matrices.
         * @param matrixA Input matrix
         * @param matrixB Input matrix
         * @param matrixC Input matrix
         * @return The result of general matrix multiply
         */
        static std::unique_ptr <RowMatrix<T>> GEMM(const RowMatrix<T>* matrixA, const RowMatrix<T>* matrixB,
                                                   const RowMatrix<T>* matrixC)
        {
            // TODO(P0): Add implementation
            return std::unique_ptr<RowMatrix<T>>(nullptr);
        }
    };
}  // namespace bustub
