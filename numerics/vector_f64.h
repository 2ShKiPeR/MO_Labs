#ifndef VECTOR_F64_H
#define VECTOR_F64_H

#include <vector>
#include <cmath>
#include <stdexcept>

namespace numerics {

class vector_f64 {
private:
    std::vector<double> data;

public:
    // Конструкторы
    vector_f64() = default;
    vector_f64(int dimension) : data(dimension, 0.0) {}
    vector_f64(std::initializer_list<double> init) : data(init) {}

    // Размерность вектора
    int dimension() const { return static_cast<int>(data.size()); }

    // Доступ к элементам
    double& operator[](int index) {
        if (index < 0 || index >= dimension()) {
            throw std::out_of_range("Vector index out of range");
        }
        return data[index];
    }

    const double& operator[](int index) const {
        if (index < 0 || index >= dimension()) {
            throw std::out_of_range("Vector index out of range");
        }
        return data[index];
    }

    // Арифметические операции
    vector_f64 operator+(const vector_f64& other) const {
        if (dimension() != other.dimension()) {
            throw std::invalid_argument("Vectors must have same dimension");
        }
        vector_f64 result(dimension());
        for (int i = 0; i < dimension(); i++) {
            result[i] = data[i] + other[i];
        }
        return result;
    }

    vector_f64 operator-(const vector_f64& other) const {
        if (dimension() != other.dimension()) {
            throw std::invalid_argument("Vectors must have same dimension");
        }
        vector_f64 result(dimension());
        for (int i = 0; i < dimension(); i++) {
            result[i] = data[i] - other[i];
        }
        return result;
    }

    vector_f64 operator*(double scalar) const {
        vector_f64 result(dimension());
        for (int i = 0; i < dimension(); i++) {
            result[i] = data[i] * scalar;
        }
        return result;
    }

    vector_f64 operator/(double scalar) const {
        if (scalar == 0.0) {
            throw std::invalid_argument("Division by zero");
        }
        vector_f64 result(dimension());
        for (int i = 0; i < dimension(); i++) {
            result[i] = data[i] / scalar;
        }
        return result;
    }

    // Статические методы
    static double distance(const vector_f64& a, const vector_f64& b) {
        if (a.dimension() != b.dimension()) {
            throw std::invalid_argument("Vectors must have same dimension");
        }
        double sum = 0.0;
        for (int i = 0; i < a.dimension(); i++) {
            double diff = a[i] - b[i];
            sum += diff * diff;
        }
        return std::sqrt(sum);
    }

    static vector_f64 direction(const vector_f64& from, const vector_f64& to) {
        if (from.dimension() != to.dimension()) {
            throw std::invalid_argument("Vectors must have same dimension");
        }
        vector_f64 result = to - from;
        double length = distance(from, to);
        if (length == 0.0) {
            throw std::invalid_argument("Vectors are identical");
        }
        return result / length;
    }

    static vector_f64 basis(int dimension, int coord) {
        if (coord < 0 || coord >= dimension) {
            throw std::out_of_range("Coordinate index out of range");
        }
        vector_f64 result(dimension);
        result[coord] = 1.0;
        return result;
    }

    // Установка конкретной координаты
    vector_f64 coord(int coordinate, double value) const {
        if (coordinate < 0 || coordinate >= dimension()) {
            throw std::out_of_range("Coordinate index out of range");
        }
        vector_f64 result = *this;
        result[coordinate] = value;
        return result;
    }
};

// Внешний оператор умножения скаляра на вектор
inline vector_f64 operator*(double scalar, const vector_f64& vec) {
    return vec * scalar;
}

} // namespace numerics

#endif