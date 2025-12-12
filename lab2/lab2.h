#ifndef LAB2_H
#define LAB2_H

#include <functional>
#include <vector>
#include "../numerics/vector_f64.h"

// Для многомерных функций
typedef double(*function_nd)(const numerics::vector_f64&);

// Структура для результатов многомерного поиска
struct multidim_search_result {
    numerics::vector_f64 result;    // Найденная точка
    double accuracy;                // Достигнутая точность
    int iterations;                 // Количество итераций
    int function_probes;            // Количество измерений функции

    void clear() {
        result = numerics::vector_f64();
        accuracy = 0;
        iterations = 0;
        function_probes = 0;
    }
};

namespace lab2 {
    // Константа золотого сечения
    static const double PSI = (std::sqrt(5.0) - 1.0) / 2.0;

    // 1. Многомерная дихотомия
    void multidim_bisect(
            multidim_search_result & result,
            function_nd function,
            const numerics::vector_f64& lhs,
            const numerics::vector_f64& rhs,
            double eps = 1e-6,
            int max_iterations = 1000
    );

    // 2. Многомерное золотое сечение
    void multidim_golden_ratio(
            multidim_search_result & result,
            function_nd function,
            const numerics::vector_f64& lhs,
            const numerics::vector_f64& rhs,
            double eps = 1e-6,
            int max_iterations = 1000
    );

    // 3. Многомерный метод Фибоначчи
    void multidim_fibonacci(
            multidim_search_result & result,
            function_nd function,
            const numerics::vector_f64& lhs,
            const numerics::vector_f64& rhs,
            double eps = 1e-6
    );

    // 4. Метод покоординатного спуска
    void coordinate_descent(
            multidim_search_result & result,
            function_nd function,
            const numerics::vector_f64& start_point,
            double lambda = 1.0,
            double eps = 1e-6,
            int max_iterations = 1000
    );

    // Вспомогательная функция для вычисления чисел Фибоначчи
    std::vector<double> generate_fibonacci_sequence(int n);
}

#endif