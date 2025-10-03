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

}