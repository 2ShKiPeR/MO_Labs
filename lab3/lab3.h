#ifndef LAB3_H
#define LAB3_H

#include <functional>
#include <vector>
#include "../numerics/vector_f64.h"
#include "../numerics/matrix_f64.h"

// Типы функций для ограничений
typedef double(*constraint_func)(const numerics::vector_f64&);

// Структура для результатов поиска
struct higher_order_search_result {
    numerics::vector_f64 result;    // Найденная точка
    double accuracy;                // Достигнутая точность
    int iterations;                 // Количество итераций
    int function_probes;            // Количество измерений функции
    double function_value;          // Значение функции в найденной точке

    void clear() {
        result = numerics::vector_f64();
        accuracy = 0;
        iterations = 0;
        function_probes = 0;
        function_value = 0;
    }
};

namespace lab3 {
    // 1. Градиентный спуск (улучшенный)
    void gradient_descent(
        higher_order_search_result & result,
        std::function<double(const numerics::vector_f64&)> function,
        const numerics::vector_f64& start_point,
        double lambda = 0.1,        // Начальная длина шага
        double eps = 1e-8,
        int max_iterations = 1000
    );

    // 2. Метод сопряжённых градиентов (улучшенный)
    void conjugate_gradient(
        higher_order_search_result & result,
        std::function<double(const numerics::vector_f64&)> function,
        const numerics::vector_f64& start_point,
        double eps = 1e-8,
        int max_iterations = 2000
    );

    // 3. Метод Ньютона-Рафсона
    void newton_raphson(
        higher_order_search_result & result,
        std::function<double(const numerics::vector_f64&)> function,
        const numerics::vector_f64& start_point,
        double eps = 1e-12,
        int max_iterations = 100
    );

    // 4. Функция внутреннего штрафа (улучшенный)
    void interior_penalty(
        higher_order_search_result & result,
        std::function<double(const numerics::vector_f64&)> function,
        const std::vector<constraint_func>& inequality_constraints,
        const numerics::vector_f64& start_point,
        double initial_penalty = 10.0,
        double eps = 1e-6,
        int max_iterations = 500
    );

    // 5. Функция внешнего штрафа
    void exterior_penalty(
        higher_order_search_result & result,
        std::function<double(const numerics::vector_f64&)> function,
        const std::vector<constraint_func>& inequality_constraints,
        const std::vector<constraint_func>& equality_constraints,
        const numerics::vector_f64& start_point,
        double penalty_param = 1.0,
        double eps = 1e-6,
        int max_iterations = 1000
    );

    // Вспомогательная функция для поиска λ
    double line_search_golden_ratio(
        std::function<double(const numerics::vector_f64&)> function,
        const numerics::vector_f64& start_point,
        const numerics::vector_f64& direction,
        double eps = 1e-6
    );
}

#endif