#include <iostream>
#include <cmath>
#include <vector>
#include "lab3/lab3.h"

// Тестовые функции
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

// Функции для штрафов
double constrained_function(const numerics::vector_f64& x) {
    return x[0] * x[0] + x[1] * x[1];
}

double constrained_function_ext(const numerics::vector_f64& x) {
    return std::pow(x[0], 2) + std::pow(x[1] - 3.0, 2);
}

// Ограничения
double inequality_constraint1(const numerics::vector_f64& x) {
    return x[0] + x[1] - 1.0;
}

double inequality_constraint2(const numerics::vector_f64& x) {
    return x[0] - 0.5;
}

double inequality_constraint_ext(const numerics::vector_f64& x) {
    return x[0] + x[1] - 2.0;
}

double equality_constraint_ext(const numerics::vector_f64& x) {
    return x[0] - x[1];
}

int main() {
    std::cout << "=== Lab 3: FINAL VERSION ===" << std::endl;

    higher_order_search_result result;

    // Тест 1: Все методы на СФЕРЕ (гарантированно работают)
    std::cout << "\n=== TESTING ON SPHERE FUNCTION ===" << std::endl;
    numerics::vector_f64 start_point({3.0, 4.0});

    std::cout << "\n1. Gradient Descent:" << std::endl;
    lab3::gradient_descent(result, sphere_function, start_point, 0.1, 1e-8);
    std::cout << "Point: (" << result.result[0] << ", " << result.result[1] << ")" << std::endl;
    std::cout << "Value: " << result.function_value << " (ideal: 0)" << std::endl;

    std::cout << "\n2. Conjugate Gradient:" << std::endl;
    lab3::conjugate_gradient(result, sphere_function, start_point, 1e-8);
    std::cout << "Point: (" << result.result[0] << ", " << result.result[1] << ")" << std::endl;
    std::cout << "Value: " << result.function_value << " (ideal: 0)" << std::endl;

    std::cout << "\n3. Newton-Raphson:" << std::endl;
    lab3::newton_raphson(result, sphere_function, start_point, 1e-12);
    std::cout << "Point: (" << result.result[0] << ", " << result.result[1] << ")" << std::endl;
    std::cout << "Value: " << result.function_value << " (ideal: 0)" << std::endl;

    // Тест 2: Методы штрафа
    std::cout << "\n=== TESTING PENALTY METHODS ===" << std::endl;

    std::cout << "\n4. Interior Penalty:" << std::endl;
    std::cout << "Problem: min(x²+y²) with x+y≥1, x≥0.5" << std::endl;
    std::vector<constraint_func> ineq_constraints = {inequality_constraint1, inequality_constraint2};
    numerics::vector_f64 start_point2({2.0, 2.0});

    lab3::interior_penalty(result, constrained_function, ineq_constraints, start_point2, 1.0, 1e-6);
    std::cout << "Point: (" << result.result[0] << ", " << result.result[1] << ")" << std::endl;
    std::cout << "Value: " << result.function_value << " (ideal: 0.5 at (0.5,0.5))" << std::endl;
    std::cout << "Constraints: " << inequality_constraint1(result.result) << ", "
              << inequality_constraint2(result.result) << std::endl;

    std::cout << "\n5. Exterior Penalty:" << std::endl;
    std::cout << "Problem: min((x-3)²+(y-2)²) with x+y≤2, x=y" << std::endl;
    std::vector<constraint_func> ext_ineq = {inequality_constraint_ext};
    std::vector<constraint_func> ext_eq = {equality_constraint_ext};
    numerics::vector_f64 start_point3({0.0, 0.0});

    lab3::exterior_penalty(result, constrained_function_ext, ext_ineq, ext_eq, start_point3, 1.0, 1e-6);
    std::cout << "Point: (" << result.result[0] << ", " << result.result[1] << ")" << std::endl;
    std::cout << "Value: " << result.function_value << " (on boundary)" << std::endl;
    std::cout << "Constraints: " << inequality_constraint_ext(result.result) << ", "
              << equality_constraint_ext(result.result) << std::endl;

    return 0;
}