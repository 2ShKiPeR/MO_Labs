#include "lab1.h"
#include <cmath>
#include <algorithm>

namespace lab1 {

void bisect(
    search_result & result,
    function_1d function,
    F64 left,
    F64 right,
    F64 eps,
    I32 max_iterations)
{
    result.clear();
    result.type = search_method_type::bisection;

    if (left > right)
        std::swap(left, right);

    for (; result.iterations != max_iterations && (result.accuracy = std::abs(right - left)) > 2 * eps;
           result.iterations++, result.function_probes += 2)
    {
        result.result = (right + left) * 0.5;

        if (function(result.result - eps * 1e-1) > function(result.result + eps * 1e-1))
            left = result.result;
        else
            right = result.result;
    }
}

void golden_ratio(
    search_result & result,
    function_1d function,
    F64 left,
    F64 right,
    F64 eps,
    I32 max_iterations)
{
    result.clear();
    result.type = search_method_type::golden_ratio;

    if (left > right)
        std::swap(left, right);

    I32 iteration = 0;
    F64 x_l = right - (right - left) * PSI;
    F64 x_r = left + (right - left) * PSI;
    F64 f_l = function(x_l);
    F64 f_r = function(x_r);

    for (; result.iterations != max_iterations && (result.accuracy = std::abs(right - left)) > 2 * eps; result.iterations++)
    {
        if (f_l > f_r)
        {
            left = x_l;
            x_l = x_r;
            f_l = f_r;
            x_r = left + (right - left) * PSI;
            f_r = function(x_r);
        }
        else
        {
            right = x_r;
            x_r = x_l;
            f_r = f_l;
            x_l = right - (right - left) * PSI;
            f_l = function(x_l);
        }
    }
    result.result = (right + left) * 0.5;
    result.function_probes = result.iterations + 2;
}

void fibonacci(
    search_result & result,
    function_1d function,
    F64 left,
    F64 right,
    F64 eps)
{
    result.clear();
    result.type = search_method_type::fibonacchi;

    if (left > right)
        std::swap(left, right);

    F64 condition = (right - left) / eps;
    F64 fib_t = 0.0, fib_1 = 1.0, fib_2 = 1.0;

    while (fib_2 < condition)
    {
        result.iterations++;
        fib_t = fib_1;
        fib_1 = fib_2;
        fib_2 += fib_t;
    }
    result.function_probes = result.iterations + 2;

    F64 x_l = left + (right - left) * ((fib_2 - fib_1) / fib_2);
    F64 x_r = left + (right - left) * (fib_1 / fib_2);

    F64 f_l = function(x_l);
    F64 f_r = function(x_r);

    for (I32 index = result.iterations; index; index--)
    {
        fib_t = fib_2 - fib_1;
        fib_2 = fib_1;
        fib_1 = fib_t;

        if (f_l > f_r)
        {
            left = x_l;
            f_l = f_r;
            x_l = x_r;
            x_r = left + (right - left) * (fib_1 / fib_2);
            f_r = function(x_r);
        }
        else
        {
            right = x_r;
            x_r = x_l;
            f_r = f_l;
            x_l = left + (right - left) * ((fib_2 - fib_1) / fib_2);
            f_l = function(x_l);
        }
    }
    result.result = (right + left) * 0.5;
    result.accuracy = (right - left);
}

} // namespace lab1