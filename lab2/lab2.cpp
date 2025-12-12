#include "lab2.h"
#include <cmath>
#include <iostream>
#include <algorithm>

namespace lab2 {

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
            // Вычисляем длину отрезка между векторами
            double segment_length = numerics::vector_f64::distance(current_lhs, current_rhs);
            result.accuracy = segment_length;

            // Проверяем условие остановки
            if (segment_length < 2 * eps) {
                break;
            }

            // Центральная точка
            numerics::vector_f64 xc = (current_lhs + current_rhs) * 0.5;

            // Вектор направления (нормированный) * 0.1 * eps
            numerics::vector_f64 dir = numerics::vector_f64::direction(current_lhs, current_rhs);
            dir = dir * (0.1 * eps);

            // Точки xl и xr
            numerics::vector_f64 xl = xc - dir;
            numerics::vector_f64 xr = xc + dir;

            // Вычисляем значения функции
            double f_xl = function(xl);
            double f_xr = function(xr);
            result.function_probes += 2;

            // Преобразуем промежуток неопределённости
            if (f_xl > f_xr) {
                current_lhs = xl;
            } else {
                current_rhs = xr;
            }
        }

        // Финальный результат - середина отрезка
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

    std::vector<double> generate_fibonacci_sequence(int n) {
        std::vector<double> fib(n + 2);
        fib[0] = 1.0;
        fib[1] = 1.0;
        for (int i = 2; i < n + 2; i++) {
            fib[i] = fib[i - 1] + fib[i - 2];
        }
        return fib;
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

        double distance_lr = numerics::vector_f64::distance(lhs, rhs);
        double condition = distance_lr / eps;

        // Генерируем последовательность Фибоначчи
        int n = 2;
        std::vector<double> fib = generate_fibonacci_sequence(n);
        while (fib[n] <= condition && n < 100) {  // Ограничиваем n для избежания переполнения
            n++;
            fib = generate_fibonacci_sequence(n);
        }

        // Если n слишком большое, ограничиваем
        if (n >= 100) {
            n = 100;
            fib = generate_fibonacci_sequence(n);
        }

        numerics::vector_f64 current_lhs = lhs;
        numerics::vector_f64 current_rhs = rhs;

        // Начальные точки по правилу Фибоначчи
        double Ln = fib[n];
        double Ln_1 = fib[n - 1];

        numerics::vector_f64 xl = current_lhs + (current_rhs - current_lhs) * ((Ln - Ln_1) / Ln);
        numerics::vector_f64 xr = current_lhs + (current_rhs - current_lhs) * (Ln_1 / Ln);

        double f_xl = function(xl);
        double f_xr = function(xr);
        result.function_probes += 2;

        // Мера близости точек Δ
        double delta = distance_lr / 100.0;
        numerics::vector_f64 dir_vector = (rhs - lhs) / distance_lr;

        for (int index = n - 1; index >= 1; index--) {
            result.iterations++;

            // Сдвиг пары чисел Фибоначчи
            double Ln_2 = fib[index - 1];
            Ln = Ln_1;
            Ln_1 = Ln_2;

            // Проверяем близость точек
            double distance_xl_xr = numerics::vector_f64::distance(xl, xr);

            if (f_xl > f_xr) {
                current_lhs = xl;
                xl = xr;
                f_xl = f_xr;
                xr = current_lhs + (current_rhs - current_lhs) * (Ln_1 / Ln);

                // Если точки слишком близки
                if (distance_xl_xr < delta) {
                    xr = xr + dir_vector * delta;
                }

                f_xr = function(xr);
                result.function_probes++;
            } else {
                current_rhs = xr;
                xr = xl;
                f_xr = f_xl;
                xl = current_lhs + (current_rhs - current_lhs) * ((Ln - Ln_1) / Ln);

                // Если точки слишком близки
                if (distance_xl_xr < delta) {
                    xl = xl - dir_vector * delta;
                }

                f_xl = function(xl);
                result.function_probes++;
            }
        }

        result.result = (current_lhs + current_rhs) * 0.5;
        result.accuracy = numerics::vector_f64::distance(current_lhs, current_rhs) * 0.5;
    }

    void coordinate_descent(
            multidim_search_result & result,
            function_nd function,
            const numerics::vector_f64& start_point,
            double lambda,
            double eps,
            int max_iterations)
    {
        result.clear();
        result.result = start_point;

        int n = start_point.dimension();
        numerics::vector_f64 current_point = start_point;

        for (; result.iterations < max_iterations; result.iterations++) {
            numerics::vector_f64 old_point = current_point;

            // Проходим по всем координатам
            for (int coord = 0; coord < n; coord++) {
                // Базисный вектор для текущей координаты
                numerics::vector_f64 basis = numerics::vector_f64::basis(n, coord);

                // Исследуем окрестность для определения направления
                numerics::vector_f64 xl = current_point - basis * eps;
                numerics::vector_f64 xr = current_point + basis * eps;

                double f_xl = function(xl);
                double f_xr = function(xr);
                result.function_probes += 2;

                // Определяем направление поиска
                numerics::vector_f64 search_end_point;
                if (f_xl > f_xr) {
                    search_end_point = current_point + basis * lambda;
                } else {
                    search_end_point = current_point - basis * lambda;
                }

                // Используем многомерную дихотомию для поиска минимума на отрезке
                multidim_search_result line_search_result;
                multidim_bisect(
                        line_search_result,
                        function,
                        current_point,
                        search_end_point,
                        eps,
                        100
                );

                // Обновляем текущую точку
                current_point = line_search_result.result;
                result.function_probes += line_search_result.function_probes;
            }

            // Проверяем условие остановки
            double movement = numerics::vector_f64::distance(old_point, current_point);
            result.accuracy = movement;

            if (movement < 2 * eps) {
                break;
            }
        }

        result.result = current_point;
    }

} // namespace lab2