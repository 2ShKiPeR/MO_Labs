#include "lab2.h"
#include <cmath>
#include <iostream>

namespace lab2 {

// Константа золотого сечения
static const double PSI = (std::sqrt(5.0) - 1.0) / 2.0;

void multidim_bisect(
    multidim_search_result & result,
    function_nd function,
    const numerics::vector_f64& lhs,
    const numerics::vector_f64& rhs,
    double eps,
    int max_iterations)
{
    result.clear();

    if (lhs.dimension() != rhs.dimension()) {
        std::cerr << "Error: Vectors have different dimensions!" << std::endl;
        return;
    }

    numerics::vector_f64 current_lhs = lhs;
    numerics::vector_f64 current_rhs = rhs;

    for (; result.iterations < max_iterations; result.iterations++) {
        double segment_length = numerics::vector_f64::distance(current_lhs, current_rhs);
        result.accuracy = segment_length;

        if (segment_length < 2 * eps) {
            break;
        }

        // Центральная точка
        numerics::vector_f64 xc = (current_lhs + current_rhs) * 0.5;

        // Вектор направления
        numerics::vector_f64 dir = numerics::vector_f64::direction(current_lhs, current_rhs);
        dir = dir * (0.1 * eps);

        // Точки xl и xr
        numerics::vector_f64 xl = xc - dir;
        numerics::vector_f64 xr = xc + dir;

        double f_xl = function(xl);
        double f_xr = function(xr);
        result.function_probes += 2;

        if (f_xl > f_xr) {
            current_lhs = xl;
        } else {
            current_rhs = xr;
        }
    }

    result.result = (current_lhs + current_rhs) * 0.5;
    result.accuracy = numerics::vector_f64::distance(current_lhs, current_rhs) * 0.5;
}

void multidim_golden_ratio(
    multidim_search_result & result,
    function_nd function,
    const numerics::vector_f64& lhs,
    const numerics::vector_f64& rhs,
    double eps,
    int max_iterations)
{
    result.clear();

    if (lhs.dimension() != rhs.dimension()) {
        std::cerr << "Error: Vectors have different dimensions!" << std::endl;
        return;
    }

    numerics::vector_f64 current_lhs = lhs;
    numerics::vector_f64 current_rhs = rhs;

    // Начальные точки по правилу золотого сечения
    numerics::vector_f64 xl = current_rhs - (current_rhs - current_lhs) * PSI;
    numerics::vector_f64 xr = current_lhs + (current_rhs - current_lhs) * PSI;

    double f_xl = function(xl);
    double f_xr = function(xr);
    result.function_probes += 2;

    for (; result.iterations < max_iterations; result.iterations++) {
        double segment_length = numerics::vector_f64::distance(current_lhs, current_rhs);
        result.accuracy = segment_length;

        if (segment_length < 2 * eps) {
            break;
        }

        if (f_xl > f_xr) {
            current_lhs = xl;
            xl = xr;
            f_xl = f_xr;
            xr = current_lhs + (current_rhs - current_lhs) * PSI;
            f_xr = function(xr);
            result.function_probes++;
        } else {
            current_rhs = xr;
            xr = xl;
            f_xr = f_xl;
            xl = current_rhs - (current_rhs - current_lhs) * PSI;
            f_xl = function(xl);
            result.function_probes++;
        }
    }

    result.result = (current_lhs + current_rhs) * 0.5;
    result.accuracy = numerics::vector_f64::distance(current_lhs, current_rhs) * 0.5;
}

void multidim_fibonacci(
    multidim_search_result & result,
    function_nd function,
    const numerics::vector_f64& lhs,
    const numerics::vector_f64& rhs,
    double eps)
{
    result.clear();

    if (lhs.dimension() != rhs.dimension()) {
        std::cerr << "Error: Vectors have different dimensions!" << std::endl;
        return;
    }

    double condition = numerics::vector_f64::distance(lhs, rhs) / eps;
    double fib_t = 0.0, fib_1 = 1.0, fib_2 = 1.0;

    // Определяем количество итераций
    while (fib_2 < condition) {
        result.iterations++;
        fib_t = fib_1;
        fib_1 = fib_2;
        fib_2 += fib_t;
    }
    result.function_probes = result.iterations + 2;

    numerics::vector_f64 current_lhs = lhs;
    numerics::vector_f64 current_rhs = rhs;

    // Начальные точки по правилу Фибоначчи
    numerics::vector_f64 xl = current_lhs + (current_rhs - current_lhs) * ((fib_2 - fib_1) / fib_2);
    numerics::vector_f64 xr = current_lhs + (current_rhs - current_lhs) * (fib_1 / fib_2);

    double f_xl = function(xl);
    double f_xr = function(xr);

    for (int index = result.iterations; index; index--) {
        fib_t = fib_2 - fib_1;
        fib_2 = fib_1;
        fib_1 = fib_t;

        if (f_xl > f_xr) {
            current_lhs = xl;
            f_xl = f_xr;
            xl = xr;
            xr = current_lhs + (current_rhs - current_lhs) * (fib_1 / fib_2);
            f_xr = function(xr);
            result.function_probes++;
        } else {
            current_rhs = xr;
            xr = xl;
            f_xr = f_xl;
            xl = current_lhs + (current_rhs - current_lhs) * ((fib_2 - fib_1) / fib_2);
            f_xl = function(xl);
            result.function_probes++;
        }
    }

    result.result = (current_lhs + current_rhs) * 0.5;
    result.accuracy = numerics::vector_f64::distance(current_lhs, current_rhs);
}

void coordinate_descent(
    multidim_search_result & result,
    function_nd function,
    const numerics::vector_f64& start_point,
    double eps,
    int max_iterations)
{
    result.clear();
    result.result = start_point;
    int n = start_point.dimension();

    numerics::vector_f64 current_point = start_point;
    double prev_value = function(current_point);
    result.function_probes++;

    for (; result.iterations < max_iterations; result.iterations++) {
        numerics::vector_f64 old_point = current_point;

        // Проходим по всем координатам
        for (int coord = 0; coord < n; coord++) {
            // Базисный вектор для текущей координаты
            numerics::vector_f64 basis = numerics::vector_f64::basis(n, coord);

            // Используем многомерную дихотомию вдоль этого направления
            multidim_search_result line_search_result;
            multidim_bisect(
                line_search_result,
                function,
                current_point - basis * 10.0,  // Левая граница
                current_point + basis * 10.0,  // Правая граница
                eps,
                100
            );

            // Обновляем текущую точку
            current_point = line_search_result.result;
            result.function_probes += line_search_result.function_probes;
        }

        // Проверяем условие сходимости
        double movement = numerics::vector_f64::distance(old_point, current_point);
        double current_value = function(current_point);
        result.function_probes++;

        if (movement < eps && std::abs(current_value - prev_value) < eps) {
            break;
        }

        prev_value = current_value;
    }

    result.result = current_point;
    result.accuracy = numerics::vector_f64::distance(start_point, current_point);
}

} // namespace lab2