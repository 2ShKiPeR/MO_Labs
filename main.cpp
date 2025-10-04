#include <iostream>
#include <cmath>
#include "lab1/lab1.h"

F64 quadratic_function(F64 x) {
    return x * x - 4 * x + 4;
}

int main() {
    std::cout << "=== Lab 1: One-Dimensional Search Methods ===" << std::endl;

    search_result result;

    // Тестирование метода дихотомии
    std::cout << "\n1. Bisection Method:" << std::endl;
    lab1::bisect(result, quadratic_function, -10, 10, 1e-6);
    std::cout << "Minimum point: " << result.result << std::endl;
    std::cout << "Accuracy: " << result.accuracy << std::endl;
    std::cout << "Iterations: " << result.iterations << std::endl;
    std::cout << "Function probes: " << result.function_probes << std::endl;

    // Тестирование метода золотого сечения
    std::cout << "\n2. Golden Ratio Method:" << std::endl;
    lab1::golden_ratio(result, quadratic_function, -10, 10, 1e-6);
    std::cout << "Minimum point: " << result.result << std::endl;
    std::cout << "Accuracy: " << result.accuracy << std::endl;
    std::cout << "Iterations: " << result.iterations << std::endl;
    std::cout << "Function probes: " << result.function_probes << std::endl;

    // Тестирование метода Фибоначчи
    std::cout << "\n3. Fibonacci Method:" << std::endl;
    lab1::fibonacci(result, quadratic_function, -10, 10, 1e-6);
    std::cout << "Minimum point: " << result.result << std::endl;
    std::cout << "Accuracy: " << result.accuracy << std::endl;
    std::cout << "Iterations: " << result.iterations << std::endl;
    std::cout << "Function probes: " << result.function_probes << std::endl;

    return 0;
}