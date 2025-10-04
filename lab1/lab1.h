#ifndef LAB1_H
#define LAB1_H

#include <cmath>
#include <functional>

typedef double F64;
typedef int I32;

enum class search_method_type {
    bisection,
    golden_ratio,
    fibonacchi
};

struct search_result {
    F64 result;
    F64 accuracy;
    I32 iterations;
    I32 function_probes;
    search_method_type type;

    void clear() {
        result = 0;
        accuracy = 0;
        iterations = 0;
        function_probes = 0;
    }
};

typedef F64(*function_1d)(F64);

namespace lab1 {
    // Константа золотого сечения
    static const F64 PSI = (std::sqrt(5.0) - 1.0) / 2.0;

    // Метод дихотомии
    void bisect(
        search_result & result,
        function_1d function,
        F64 left,
        F64 right,
        F64 eps = 1e-6,
        I32 max_iterations = 1000
    );

    // Метод золотого сечения
    void golden_ratio(
        search_result & result,
        function_1d function,
        F64 left,
        F64 right,
        F64 eps = 1e-6,
        I32 max_iterations = 1000
    );

    // Метод Фибоначчи
    void fibonacci(
        search_result & result,
        function_1d function,
        F64 left,
        F64 right,
        F64 eps = 1e-6
    );
}

#endif