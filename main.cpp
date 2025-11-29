#include <iostream>
#include <cmath>
#include "lab2/lab2.h"

// Многомерные тестовые функции
double sphere_function(const numerics::vector_f64& x) {
    double sum = 0;
    for (int i = 0; i < x.dimension(); i++) {
        sum += x[i] * x[i];
    }
    return sum;
}

double rosenbrock_function(const numerics::vector_f64& x) {
    return 100 * std::pow(x[1] - x[0] * x[0], 2) + std::pow(1 - x[0], 2);
}

int main() {
    std::cout << "=== Lab 2: All 4 Multidimensional Search Methods ===" << std::endl;

    multidim_search_result result;
    numerics::vector_f64 lhs({-2.0, -2.0});
    numerics::vector_f64 rhs({2.0, 2.0});
    numerics::vector_f64 start_point({-1.5, 1.5});

    // 1. Многомерная дихотомия
    std::cout << "\n1. Multidimensional Bisection Method (Sphere):" << std::endl;
    lab2::multidim_bisect(result, sphere_function, lhs, rhs, 1e-6);
    std::cout << "Minimum point: (" << result.result[0] << ", " << result.result[1] << ")" << std::endl;
    std::cout << "Accuracy: " << result.accuracy << std::endl;
    std::cout << "Iterations: " << result.iterations << std::endl;
    std::cout << "Function probes: " << result.function_probes << std::endl;

    // 2. Многомерное золотое сечение
    std::cout << "\n2. Multidimensional Golden Ratio Method (Sphere):" << std::endl;
    lab2::multidim_golden_ratio(result, sphere_function, lhs, rhs, 1e-6);
    std::cout << "Minimum point: (" << result.result[0] << ", " << result.result[1] << ")" << std::endl;
    std::cout << "Accuracy: " << result.accuracy << std::endl;
    std::cout << "Iterations: " << result.iterations << std::endl;
    std::cout << "Function probes: " << result.function_probes << std::endl;

    // 3. Многомерный метод Фибоначчи
    std::cout << "\n3. Multidimensional Fibonacci Method (Sphere):" << std::endl;
    lab2::multidim_fibonacci(result, sphere_function, lhs, rhs, 1e-6);
    std::cout << "Minimum point: (" << result.result[0] << ", " << result.result[1] << ")" << std::endl;
    std::cout << "Accuracy: " << result.accuracy << std::endl;
    std::cout << "Iterations: " << result.iterations << std::endl;
    std::cout << "Function probes: " << result.function_probes << std::endl;

    // 4. Покоординатный спуск
    std::cout << "\n4. Coordinate Descent Method (Rosenbrock):" << std::endl;
    lab2::coordinate_descent(result, rosenbrock_function, start_point, 1e-6);
    std::cout << "Minimum point: (" << result.result[0] << ", " << result.result[1] << ")" << std::endl;
    std::cout << "Accuracy: " << result.accuracy << std::endl;
    std::cout << "Iterations: " << result.iterations << std::endl;
    std::cout << "Function probes: " << result.function_probes << std::endl;

    return 0;
}