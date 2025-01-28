#ifndef ADCG_MATH_OPS_H_
#define ADCG_MATH_OPS_H_

#include "math.h"

#include "adcg_types.hpp"
#include "adcg_ascalar.hpp"

/**
 * @brief 
 * 
 * @tparam T 
 * @tparam Op
 */

template<typename T, typename Op>
class ascalarSin {

private:

    const Op& op_;

public:

    ascalarSin(const Op& op) : op_(op) {};
    
    T getValue() const { return sin(op_.getValue()); };
    std::map<std::size_t, T> getDerivative() const { 
        std::map<std::size_t, T> der(op_.getDerivative());
        for (auto i = der.begin(); i != der.end(); i++)
            der[i->first] *= cos(op_.getValue());
        return der;
    };

};

/**
 * @brief 
 * 
 * @tparam T 
 * @tparam R1 
 * @param op
 * @return ascalar<T, ascalarSin<T, R>>
 */

template<typename T, typename R>
ascalar<T, ascalarSin<T, R>> sin(const ascalar<T, R>& op) {
    return ascalar<T, ascalarSin<T, R>>(ascalarSin<T, R>(op.getData()));
};

/**
 * @brief 
 * 
 * @tparam T 
 * @tparam Op
 */

template<typename T, typename Op>
class ascalarCos {

private:

    const Op& op_;

public:

    ascalarCos(const Op& op) : op_(op) {};
    
    T getValue() const { return cos(op_.getValue()); };
    std::map<std::size_t, T> getDerivative() const { 
        std::map<std::size_t, T> der(op_.getDerivative());
        for (auto i = der.begin(); i != der.end(); i++)
            der[i->first] *= -sin(op_.getValue());
        return der;
    };

};

/**
 * @brief 
 * 
 * @tparam T 
 * @tparam R
 * @param op
 * @return ascalar<T, ascalarCos<T, R>>
 */

template<typename T, typename R>
ascalar<T, ascalarCos<T, R>> cos(const ascalar<T, R>& op) {
    return ascalar<T, ascalarCos<T, R>>(ascalarCos<T, R>(op.getData()));
};

/**
 * @brief 
 * 
 * @tparam T 
 * @tparam Op
 */

template<typename T, typename Op>
class ascalarTan {

private:

    const Op& op_;

public:

    ascalarTan(const Op& op) : op_(op) {};
    
    T getValue() const { return tan(op_.getValue()); };
    std::map<std::size_t, T> getDerivative() const { 
        std::map<std::size_t, T> der(op_.getDerivative());
        for (auto i = der.begin(); i != der.end(); i++)
            der[i->first] *= (1 + tan(op_.getValue()) * tan(op_.getValue()));
        return der;
    };

};

/**
 * @brief 
 * 
 * @tparam T 
 * @tparam R
 * @param op
 * @return ascalar<T, ascalarTan<T, R>>
 */

template<typename T, typename R>
ascalar<T, ascalarTan<T, R>> tan(const ascalar<T, R>& op) {
    return ascalar<T, ascalarTan<T, R>>(ascalarTan<T, R>(op.getData()));
};

#endif // ADCG_MATH_OPS_H_