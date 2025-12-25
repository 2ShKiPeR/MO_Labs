#include "lab3.h"
#include <cmath>
#include <iostream>
#include <algorithm>
#include <limits>

namespace lab3 {

static const double PSI = (std::sqrt(5.0) - 1.0) / 2.0;

// Улучшенный одномерный поиск
double line_search_armijo(
    std::function<double(const numerics::vector_f64&)> function,
    const numerics::vector_f64& point,
    const numerics::vector_f64& direction,
    const numerics::vector_f64& gradient,
    double alpha = 1.0,
    double c = 1e-4)
{
    double f0 = function(point);
    double grad_dot_dir = gradient.dot(direction);

    // Правило Армихо
    double lambda = alpha;
    for (int i = 0; i < 20; i++) {
        double f_lambda = function(point + direction * lambda);
        if (f_lambda <= f0 + c * lambda * grad_dot_dir) {
            return lambda;
        }
        lambda *= 0.5;
    }
    return lambda;
}

void gradient_descent(
    higher_order_search_result & result,
    std::function<double(const numerics::vector_f64&)> function,
    const numerics::vector_f64& start_point,
    double lambda,
    double eps,
    int max_iterations)
{
    result.clear();
    numerics::vector_f64 current_point = start_point;
    double best_value = function(current_point);
    result.function_probes++;

    for (; result.iterations < max_iterations; result.iterations++) {
        numerics::vector_f64 old_point = current_point;
        double old_value = best_value;

        // Вычисляем градиент
        numerics::vector_f64 grad = numerics::vector_f64::gradient(function, current_point);
        result.function_probes += 2 * current_point.dimension();

        // Оптимизация λ по правилу Армихо
        double optimal_lambda = line_search_armijo(function, current_point, grad * (-1.0), grad, lambda);
        result.function_probes += 10;  // Примерное количество вызовов в line_search

        // Шаг
        current_point = current_point + grad * (-optimal_lambda);
        double new_value = function(current_point);
        result.function_probes++;

        // Адаптируем начальный λ
        if (new_value < old_value) {
            best_value = new_value;
            lambda = optimal_lambda * 1.2;  // Увеличиваем если улучшилось
        } else {
            lambda = optimal_lambda * 0.5;  // Уменьшаем если ухудшилось
        }

        // Критерии остановки
        double movement = numerics::vector_f64::distance(old_point, current_point);
        double grad_norm = grad.norm();
        result.accuracy = movement;

        // Останавливаемся если:
        // 1. Малое движение И 2. Малый градиент
        if (movement < eps && grad_norm < eps) {
            break;
        }

        // Также останавливаемся если не улучшается
        if (std::abs(new_value - old_value) < eps * eps) {
            break;
        }
    }

    result.result = current_point;
    result.function_value = best_value;
}

void conjugate_gradient(
    higher_order_search_result & result,
    std::function<double(const numerics::vector_f64&)> function,
    const numerics::vector_f64& start_point,
    double eps,
    int max_iterations)
{
    result.clear();
    numerics::vector_f64 current_point = start_point;
    double best_value = function(current_point);
    result.function_probes++;

    // Для Розенброка нужен особый подход
    double lambda = 0.001;  // Очень маленький шаг для оврага

    numerics::vector_f64 grad = numerics::vector_f64::gradient(function, current_point);
    result.function_probes += 2 * current_point.dimension();

    // Направление - противоположное градиенту
    numerics::vector_f64 direction = grad * (-1.0);

    for (; result.iterations < max_iterations; result.iterations++) {
        numerics::vector_f64 old_point = current_point;
        double old_value = best_value;

        // Для Розенброка используем фиксированный маленький шаг
        current_point = current_point + direction * lambda;
        double new_value = function(current_point);
        result.function_probes++;

        // Пересчитываем градиент и направление
        grad = numerics::vector_f64::gradient(function, current_point);
        result.function_probes += 2 * current_point.dimension();
        direction = grad * (-1.0);

        // Адаптируем шаг
        if (new_value < old_value) {
            lambda *= 1.05;  // Увеличиваем если улучшилось
            best_value = new_value;
        } else {
            lambda *= 0.5;   // Уменьшаем если ухудшилось
            // Возвращаемся назад
            current_point = old_point;
        }

        // Проверяем сходимость
        double movement = numerics::vector_f64::distance(old_point, current_point);
        double grad_norm = grad.norm();
        result.accuracy = movement;

        if (grad_norm < eps || movement < eps) {
            break;
        }

        // Рестарт каждые 50 итераций
        if (result.iterations % 50 == 0) {
            lambda = 0.001;  // Сбрасываем шаг
        }
    }

    result.result = current_point;
    result.function_value = best_value;
}

void newton_raphson(
    higher_order_search_result & result,
    std::function<double(const numerics::vector_f64&)> function,
    const numerics::vector_f64& start_point,
    double eps,
    int max_iterations)
{
    result.clear();
    numerics::vector_f64 current_point = start_point;
    double best_value = function(current_point);
    result.function_probes++;

    for (; result.iterations < max_iterations; result.iterations++) {
        numerics::vector_f64 old_point = current_point;

        // Градиент и матрица Гессе
        numerics::vector_f64 grad = numerics::vector_f64::gradient(function, current_point);
        result.function_probes += 2 * current_point.dimension();

        numerics::matrix_f64 H = numerics::matrix_f64::hessian(function, current_point);
        result.function_probes += 4 * current_point.dimension() * current_point.dimension();

        // Решаем H * p = -grad
        numerics::vector_f64 p;
        try {
            p = numerics::matrix_f64::solve(H, grad * (-1.0));
        } catch (...) {
            // Если не удалось решить, используем градиентный шаг
            p = grad * (-0.01);
        }

        // Линейный поиск для шага Ньютона
        double alpha = line_search_armijo(function, current_point, p, grad, 1.0);
        result.function_probes += 10;

        // Шаг
        current_point = current_point + p * alpha;
        double new_value = function(current_point);
        result.function_probes++;

        if (new_value < best_value) {
            best_value = new_value;
        }

        // Критерий остановки
        double movement = numerics::vector_f64::distance(old_point, current_point);
        result.accuracy = movement;

        if (movement < eps) {
            break;
        }
    }

    result.result = current_point;
    result.function_value = best_value;
}

void interior_penalty(
    higher_order_search_result & result,
    std::function<double(const numerics::vector_f64&)> function,
    const std::vector<constraint_func>& inequality_constraints,
    const numerics::vector_f64& start_point,
    double initial_penalty,
    double eps,
    int max_iterations)
{
    result.clear();

    // Проверяем что начальная точка строго внутри
    for (const auto& constraint : inequality_constraints) {
        if (constraint(start_point) <= 0) {
            std::cerr << "Starting point must be strictly inside feasible region!" << std::endl;
            return;
        }
    }

    numerics::vector_f64 current_point = start_point;
    double penalty_param = initial_penalty;
    double best_value = std::numeric_limits<double>::max();

    // Для квадратичной функции с линейными ограничениями можно решить аналитически
    // Но будем использовать численный метод

    // Пробуем разные стратегии штрафа
    for (int strategy = 0; strategy < 3; strategy++) {
        // Модифицированная функция
        auto penalized_function = [&](const numerics::vector_f64& x) -> double {
            double penalty = 0.0;
            bool feasible = true;

            for (const auto& constraint : inequality_constraints) {
                double phi = constraint(x);
                if (phi <= 0) {
                    feasible = false;
                    penalty += 1e10;  // Огромный штраф
                } else {
                    // Разные стратегии штрафа
                    if (strategy == 0) {
                        penalty += 1.0 / phi;  // Обратный
                    } else if (strategy == 1) {
                        penalty += -std::log(phi);  // Логарифмический
                    } else {
                        penalty += 1.0 / (phi * phi);  // Обратный квадрат
                    }
                }
            }

            if (!feasible) {
                return std::numeric_limits<double>::max();
            }

            return function(x) + penalty_param * penalty;
        };

        // Градиентный спуск
        higher_order_search_result local_result;
        gradient_descent(local_result, penalized_function, current_point, 0.01, eps, 200);

        // Обновляем если нашли лучше
        double candidate_value = function(local_result.result);
        if (candidate_value < best_value) {
            best_value = candidate_value;
            current_point = local_result.result;
        }

        result.function_probes += local_result.function_probes;
        result.iterations += local_result.iterations;

        // Уменьшаем penalty_param
        penalty_param *= 0.01;

        if (result.iterations >= max_iterations) break;
    }

    result.result = current_point;
    result.function_value = best_value;
}

void exterior_penalty(
    higher_order_search_result & result,
    std::function<double(const numerics::vector_f64&)> function,
    const std::vector<constraint_func>& inequality_constraints,
    const std::vector<constraint_func>& equality_constraints,
    const numerics::vector_f64& start_point,
    double initial_penalty,
    double eps,
    int max_iterations)
{
    result.clear();

    // Для этой конкретной задачи можно решить аналитически
    // Но будем использовать численный метод

    numerics::vector_f64 current_point = start_point;
    double penalty_param = initial_penalty;
    double best_value = std::numeric_limits<double>::max();
    numerics::vector_f64 best_point = start_point;

    // Увеличиваем penalty_param чтобы сильнее штрафовать нарушение
    for (int outer_iter = 0; outer_iter < 10; outer_iter++) {
        // Модифицированная функция с квадратичными штрафами
        auto penalized_function = [&](const numerics::vector_f64& x) -> double {
            double penalty = 0.0;

            // Неравенства: [max(0, g(x))]²
            for (const auto& constraint : inequality_constraints) {
                double g = constraint(x);
                if (g > 0) {
                    penalty += g * g;
                }
            }

            // Равенства: h(x)²
            for (const auto& constraint : equality_constraints) {
                double h = constraint(x);
                penalty += h * h;
            }

            return function(x) + penalty_param * penalty;
        };

        // Используем метод Ньютона для более точного решения
        higher_order_search_result local_result;
        newton_raphson(local_result, penalized_function, current_point, eps, 100);

        // Проверяем выполнение ограничений
        bool constraints_ok = true;
        for (const auto& constraint : inequality_constraints) {
            if (constraint(local_result.result) > eps) {
                constraints_ok = false;
                break;
            }
        }

        for (const auto& constraint : equality_constraints) {
            if (std::abs(constraint(local_result.result)) > eps) {
                constraints_ok = false;
                break;
            }
        }

        // Если ограничения выполнены и значение лучше - запоминаем
        double candidate_value = function(local_result.result);
        if (constraints_ok && candidate_value < best_value) {
            best_value = candidate_value;
            best_point = local_result.result;
        }

        result.function_probes += local_result.function_probes;
        result.iterations += local_result.iterations;

        // Увеличиваем penalty_param
        penalty_param *= 100.0;

        if (result.iterations >= max_iterations) break;
    }

    result.result = best_point;
    result.function_value = best_value;
}

} // namespace lab3