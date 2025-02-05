#include <iostream>

#include "adcg_ascalar.hpp"
#include "adcg_binary_ops.hpp"
#include "adcg_math_ops.hpp"

typedef ascalar<real> scalar;

int main() {
    
    scalar x(1.0, 0);
    scalar y(2.0, 1);
    scalar z(3.0, 2);

    scalar res = x + z;
    scalar res2 = res + y - x;
    scalar res3 = tan(sin(cos(x * res2 / y)));
    scalar res4 = pow(2, 2 - 5/log(y * exp(pow(atan(asin(acos(res3) / z) / z), 3.5))) - 0.1 );

    std::cout << "X: " << x.getValue() << std::endl;
    std::cout << "Y: " << y.getValue() << std::endl;
    std::cout << "Z: " << z.getValue() << std::endl;

    std::cout << res4.getValue() << std::endl;

    std::map<std::size_t, real> J = res4.getDerivative();

    for (auto i = J.begin(); i != J.end(); i++)
        std::cout << i->first << " : " << i->second << '\n';

    return 0;

}