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
    std::cout << "=== Lab 2: Multidimensional Search Methods (Golden Ratio in Coordinate Descent) ===" << std::endl;

    multidim_search_result result;

    // Тест 1: Многомерная дихотомия
    std::cout << "\n1. Multidimensional Bisection Method (Sphere):" << std::endl;
    numerics::vector_f64 lhs({-2.0, -2.0});
    numerics::vector_f64 rhs({2.0, 2.0});

    lab2::multidim_bisect(result, sphere_function, lhs, rhs, 1e-6);
    std::cout << "Minimum point: (" << result.result[0] << ", " << result.result[1] << ")" << std::endl;
    std::cout << "Accuracy: " << result.accuracy << std::endl;
    std::cout << "Iterations: " << result.iterations << std::endl;
    std::cout << "Function probes: " << result.function_probes << std::endl;

    // Тест 2: Многомерное золотое сечение
    std::cout << "\n2. Multidimensional Golden Ratio Method (Sphere):" << std::endl;
    lab2::multidim_golden_ratio(result, sphere_function, lhs, rhs, 1e-6);
    std::cout << "Minimum point: (" << result.result[0] << ", " << result.result[1] << ")" << std::endl;
    std::cout << "Accuracy: " << result.accuracy << std::endl;
    std::cout << "Iterations: " << result.iterations << std::endl;
    std::cout << "Function probes: " << result.function_probes << std::endl;

    // Тест 3: Многомерный метод Фибоначчи
    std::cout << "\n3. Multidimensional Fibonacci Method (Sphere):" << std::endl;
    lab2::multidim_fibonacci(result, sphere_function, lhs, rhs, 1e-6);
    std::cout << "Minimum point: (" << result.result[0] << ", " << result.result[1] << ")" << std::endl;
    std::cout << "Accuracy: " << result.accuracy << std::endl;
    std::cout << "Iterations: " << result.iterations << std::endl;
    std::cout << "Function probes: " << result.function_probes << std::endl;

    // Тест 4: Покоординатный спуск С ЗОЛОТЫМ СЕЧЕНИЕМ (основной тест)
    std::cout << "\n4. Coordinate Descent with GOLDEN RATIO Method (Rosenbrock):" << std::endl;
    multidim_search_result result2;
    numerics::vector_f64 start_point({-1.5, 1.5});

    // Оптимальные параметры для золотого сечения
    lab2::coordinate_descent(result2, rosenbrock_function, start_point,
                             5.0,   // λ = 5.0 (оптимально)
                             1e-4,  // eps = 0.0001 (для скорости)
                             200);  // max_iterations = 200

    std::cout << "Minimum point: (" << result2.result[0] << ", " << result2.result[1] << ")" << std::endl;
    std::cout << "Accuracy: " << result2.accuracy << std::endl;
    std::cout << "Iterations: " << result2.iterations << std::endl;
    std::cout << "Function probes: " << result2.function_probes << std::endl;

    return 0;
}