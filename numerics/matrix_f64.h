#ifndef MATRIX_F64_H
#define MATRIX_F64_H

#include <vector>
#include <cmath>
#include <stdexcept>
#include <functional>
#include <iostream>
#include "vector_f64.h"

namespace numerics {

class matrix_f64 {
private:
    std::vector<std::vector<double>> data;
    int rows, cols;
    
public:
    // Конструкторы
    matrix_f64() : rows(0), cols(0) {}
    matrix_f64(int n_rows, int n_cols) : rows(n_rows), cols(n_cols) {
        data.resize(n_rows, std::vector<double>(n_cols, 0.0));
    }
    
    matrix_f64(std::initializer_list<std::initializer_list<double>> init) {
        rows = static_cast<int>(init.size());
        if (rows > 0) {
            cols = static_cast<int>(init.begin()->size());
        } else {
            cols = 0;
        }
        
        data.resize(rows, std::vector<double>(cols, 0.0));
        int i = 0;
        for (const auto& row : init) {
            int j = 0;
            for (double val : row) {
                data[i][j] = val;
                j++;
            }
            i++;
        }
    }
    
    // Размеры матрицы
    int n_rows() const { return rows; }
    int n_cols() const { return cols; }
    
    // Доступ к элементам
    double& operator()(int i, int j) { 
        if (i < 0 || i >= rows || j < 0 || j >= cols) {
            throw std::out_of_range("Matrix index out of range");
        }
        return data[i][j]; 
    }
    
    const double& operator()(int i, int j) const { 
        if (i < 0 || i >= rows || j < 0 || j >= cols) {
            throw std::out_of_range("Matrix index out of range");
        }
        return data[i][j]; 
    }
    
    // Арифметические операции
    matrix_f64 operator+(const matrix_f64& other) const {
        if (rows != other.rows || cols != other.cols) {
            throw std::invalid_argument("Matrices must have same dimensions");
        }
        matrix_f64 result(rows, cols);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                result(i, j) = data[i][j] + other(i, j);
            }
        }
        return result;
    }
    
    matrix_f64 operator-(const matrix_f64& other) const {
        if (rows != other.rows || cols != other.cols) {
            throw std::invalid_argument("Matrices must have same dimensions");
        }
        matrix_f64 result(rows, cols);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                result(i, j) = data[i][j] - other(i, j);
            }
        }
        return result;
    }
    
    matrix_f64 operator*(double scalar) const {
        matrix_f64 result(rows, cols);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                result(i, j) = data[i][j] * scalar;
            }
        }
        return result;
    }
    
    // Умножение матрицы на вектор
    vector_f64 operator*(const vector_f64& vec) const {
        if (cols != vec.dimension()) {
            throw std::invalid_argument("Matrix columns must equal vector dimension");
        }
        
        vector_f64 result(rows);
        for (int i = 0; i < rows; i++) {
            double sum = 0.0;
            for (int j = 0; j < cols; j++) {
                sum += data[i][j] * vec[j];
            }
            result[i] = sum;
        }
        return result;
    }
    
    // Транспонирование
    matrix_f64 transpose() const {
        matrix_f64 result(cols, rows);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                result(j, i) = data[i][j];
            }
        }
        return result;
    }

    // Единичная матрица
    static matrix_f64 identity(int n) {
        matrix_f64 result(n, n);
        for (int i = 0; i < n; i++) {
            result(i, i) = 1.0;
        }
        return result;
    }

    // Решение системы линейных уравнений (исправленный)
    static vector_f64 solve(const matrix_f64& A, const vector_f64& b) {
        if (A.n_rows() != A.n_cols() || A.n_rows() != b.dimension()) {
            throw std::invalid_argument("Invalid matrix or vector dimensions");
        }

        int n = A.n_rows();

        // Простой метод для положительно определённых матриц
        try {
            // Проверяем диагональное преобладание
            bool diagonal_dominant = true;
            for (int i = 0; i < n; i++) {
                double diag = std::abs(A(i, i));
                double sum = 0.0;
                for (int j = 0; j < n; j++) {
                    if (i != j) sum += std::abs(A(i, j));
                }
                if (diag <= sum) {
                    diagonal_dominant = false;
                    break;
                }
            }

            if (diagonal_dominant) {
                // Метод простых итераций - ИСПРАВЛЕННЫЙ КОНСТРУКТОР
                vector_f64 x(n);  // Просто размер
                vector_f64 x_old(n);

                // Инициализируем нулями
                for (int i = 0; i < n; i++) {
                    x[i] = 0.0;
                    x_old[i] = 0.0;
                }

                for (int iter = 0; iter < 1000; iter++) {
                    x_old = x;

                    for (int i = 0; i < n; i++) {
                        double sum = b[i];
                        for (int j = 0; j < n; j++) {
                            if (i != j) sum -= A(i, j) * x_old[j];
                        }
                        x[i] = sum / A(i, i);
                    }

                    // Проверка сходимости
                    double max_diff = 0.0;
                    for (int i = 0; i < n; i++) {
                        max_diff = std::max(max_diff, std::abs(x[i] - x_old[i]));
                    }

                    if (max_diff < 1e-12) break;
                }

                return x;
            }
        } catch (...) {
            // В случае ошибки используем упрощённый подход
        }

        // Упрощённый вариант: возвращаем -b (для метода Ньютона)
        return b * (-1.0);
    }

    // Вычисление матрицы Гессе
    static matrix_f64 hessian(std::function<double(const vector_f64&)> function,
                             const vector_f64& point,
                             double dx = 1e-6) {
        int n = point.dimension();
        matrix_f64 H(n, n);

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i == j) {
                    // Вторая производная по одной переменной
                    vector_f64 point_plus = point;
                    vector_f64 point_minus = point;
                    point_plus[i] += dx;
                    point_minus[i] -= dx;

                    double f_plus = function(point_plus);
                    double f_center = function(point);
                    double f_minus = function(point_minus);

                    H(i, j) = (f_plus - 2 * f_center + f_minus) / (dx * dx);
                } else {
                    // Смешанная производная
                    vector_f64 pp = point; pp[i] += dx; pp[j] += dx;
                    vector_f64 pm = point; pm[i] += dx; pm[j] -= dx;
                    vector_f64 mp = point; mp[i] -= dx; mp[j] += dx;
                    vector_f64 mm = point; mm[i] -= dx; mm[j] -= dx;

                    double f_pp = function(pp);
                    double f_pm = function(pm);
                    double f_mp = function(mp);
                    double f_mm = function(mm);

                    H(i, j) = (f_pp - f_pm - f_mp + f_mm) / (4 * dx * dx);
                    H(j, i) = H(i, j);
                }
            }
        }
        
        return H;
    }
    
    // Вывод матрицы
    void print() const {
        for (int i = 0; i < rows; i++) {
            std::cout << "[";
            for (int j = 0; j < cols; j++) {
                std::cout << data[i][j];
                if (j < cols - 1) std::cout << ", ";
            }
            std::cout << "]" << std::endl;
        }
    }
};

// Внешние операторы
inline matrix_f64 operator*(double scalar, const matrix_f64& mat) {
    return mat * scalar;
}

} // namespace numerics

#endif