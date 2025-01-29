#include <iostream>

#include "adcg_ascalar.hpp"
#include "adcg_binary_ops.hpp"
#include "adcg_math_ops.hpp"

int main() {
    
    ascalar<real> x(1.0, 0);
    ascalar<real> y(2.0, 1);
    ascalar<real> z(3.0, 2);

    ascalar<real> res = x + z;
    ascalar<real> res2 = res + y - x;
    ascalar<real> res3 = tan(sin(cos(x * res2 / y)));

    std::cout << "X: " << x.getValue() << std::endl;
    std::cout << "Y: " << y.getValue() << std::endl;
    std::cout << "Z: " << z.getValue() << std::endl;

    std::cout << res3.getValue() << std::endl;

    std::map<std::size_t, real> J = res3.getDerivative();

    for (auto i = J.begin(); i != J.end(); i++)
        std::cout << i->first << " : " << i->second << '\n';

    return 0;

}