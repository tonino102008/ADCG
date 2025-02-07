#include <iostream>

#include "adcg_ascalar.hpp"
#include "adcg_binary_ops.hpp"
#include "adcg_math_ops.hpp"

typedef ascalar<ascalar<real>> scalar;

int main() {

    scalar x(1.0, 0);
    scalar y(2.0, 1);

    x.getValue().getDerivative()[0] = 1;
    y.getValue().getDerivative()[1] = 1;

    scalar res = y * x * x + y * y / x;

    std::cout << "X: " << x.getValue().getValue() << std::endl;
    std::cout << "Y: " << y.getValue().getValue() << std::endl;
    std::cout << std::endl;

    for (auto i = x.getDerivative().begin(); i != x.getDerivative().end(); i++) {
        std::cout << i->first << " : " << i->second.getValue() << '\n';
        for (auto j = i->second.getDerivative().begin(); j != i->second.getDerivative().end(); j++)
            std::cout << j->first << " : " << j->second << '\n';
        std::cout << std::endl;
    }
    std::cout << std::endl;

    for (auto i = y.getDerivative().begin(); i != y.getDerivative().end(); i++) {
        std::cout << i->first << " : " << i->second.getValue() << '\n';
        for (auto j = i->second.getDerivative().begin(); j != i->second.getDerivative().end(); j++)
            std::cout << j->first << " : " << j->second << '\n';
        std::cout << std::endl;
    }
    std::cout << std::endl;

    std::cout << "FUNCTION" << std::endl;
    std::cout << res.getValue().getValue() << std::endl;
    std::cout << std::endl;

    std::map<std::size_t, ascalar<real>> J = res.getDerivative();

    std::cout << "JACOBIAN" << std::endl;
    for (auto i = J.begin(); i != J.end(); i++) {
        std::cout << i->first << " : " << i->second.getValue() << '\n';
        for (auto j = i->second.getDerivative().begin(); j != i->second.getDerivative().end(); j++)
            std::cout << j->first << " : " << j->second << '\n';
        std::cout << std::endl;
    }
    std::cout << std::endl;

    return 0;

}