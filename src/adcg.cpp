#include <iostream>

#include "adcg_ascalar.hpp"

int main() {
    
    ascalar<real> x(1.0, 0);
    ascalar<real> y(2.0, 1);
    ascalar<real> z(3.0, 2);

    ascalar<real> res(0.0);
    ascalar<real> res2(0.0);
    ascalar<real> res3(0.0);

    res = x + z;
    res2 = res + y + x;
    res3 = x * y * res2;

    std::cout << "X: " << x.getValue() << std::endl;
    std::cout << "Y: " << y.getValue() << std::endl;
    std::cout << "Z: " << z.getValue() << std::endl;

    std::cout << res3.getValue() << std::endl;

    std::map<std::size_t, real> J = res3.getDerivative();

    for (auto i = J.begin(); i != J.end(); i++)
        std::cout << i->first << " : " << i->second << '\n';

    return 0;

}