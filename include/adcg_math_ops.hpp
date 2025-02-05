#ifndef ADCG_MATH_OPS_H_
#define ADCG_MATH_OPS_H_

#include "math.h"
#include "cassert"

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

/**
 * @brief 
 * 
 * @tparam T 
 * @tparam Op
 */

template<typename T, typename Op>
class ascalarAsin {

private:

    const Op& op_;

public:

    ascalarAsin(const Op& op) : op_(op) {};
    
    T getValue() const { return asin(op_.getValue()); };
    std::map<std::size_t, T> getDerivative() const { 
        assert((op_.getValue() > -1) && (op_.getValue() < 1));
        std::map<std::size_t, T> der(op_.getDerivative());
        for (auto i = der.begin(); i != der.end(); i++)
            der[i->first] /= sqrt(1 - pow(op_.getValue(), 2));
        return der;
    };

};

/**
 * @brief 
 * 
 * @tparam T 
 * @tparam R1 
 * @param op
 * @return ascalar<T, ascalarAsin<T, R>>
 */

template<typename T, typename R>
ascalar<T, ascalarAsin<T, R>> asin(const ascalar<T, R>& op) {
    return ascalar<T, ascalarAsin<T, R>>(ascalarAsin<T, R>(op.getData()));
};

/**
 * @brief 
 * 
 * @tparam T 
 * @tparam Op
 */

template<typename T, typename Op>
class ascalarAcos {

private:

    const Op& op_;

public:

    ascalarAcos(const Op& op) : op_(op) {};
    
    T getValue() const { return acos(op_.getValue()); };
    std::map<std::size_t, T> getDerivative() const { 
        assert((op_.getValue() > -1) && (op_.getValue() < 1));
        std::map<std::size_t, T> der(op_.getDerivative());
        for (auto i = der.begin(); i != der.end(); i++)
            der[i->first] /= -sqrt(1 - pow(op_.getValue(), 2));
        return der;
    };

};

/**
 * @brief 
 * 
 * @tparam T 
 * @tparam R1 
 * @param op
 * @return ascalar<T, ascalarAcos<T, R>>
 */

template<typename T, typename R>
ascalar<T, ascalarAcos<T, R>> acos(const ascalar<T, R>& op) {
    return ascalar<T, ascalarAcos<T, R>>(ascalarAcos<T, R>(op.getData()));
};

/**
 * @brief 
 * 
 * @tparam T 
 * @tparam Op
 */

template<typename T, typename Op>
class ascalarAtan {

private:

    const Op& op_;

public:

    ascalarAtan(const Op& op) : op_(op) {};
    
    T getValue() const { return atan(op_.getValue()); };
    std::map<std::size_t, T> getDerivative() const { 
        std::map<std::size_t, T> der(op_.getDerivative());
        for (auto i = der.begin(); i != der.end(); i++)
            der[i->first] /= (1 + pow(op_.getValue(), 2));
        return der;
    };

};

/**
 * @brief 
 * 
 * @tparam T 
 * @tparam R1 
 * @param op
 * @return ascalar<T, ascalarAtan<T, R>>
 */

template<typename T, typename R>
ascalar<T, ascalarAtan<T, R>> atan(const ascalar<T, R>& op) {
    return ascalar<T, ascalarAtan<T, R>>(ascalarAtan<T, R>(op.getData()));
};

/**
 * @brief 
 * 
 * @tparam T 
 * @tparam Op 
 */

template<typename T, typename Op>
class ascalarPow {

private:

    const Op& op_;
    const double& n_;

public:

    ascalarPow(const Op& op, const double& n) : op_(op), n_(n) {};
    
    T getValue() const { return pow(op_.getValue(), n_); };
    std::map<std::size_t, T> getDerivative() const { 
        std::map<std::size_t, T> der(op_.getDerivative());
        for (auto i = der.begin(); i != der.end(); i++)
            der[i->first] *= n_ * pow(op_.getValue(), n_ - 1);
        return der;
    };

};

/**
 * @brief 
 * 
 * @tparam T 
 * @tparam R 
 * @param op 
 * @param n 
 * @return ascalar<T, ascalarPow<T, R>> 
 */

template<typename T, typename R>
ascalar<T, ascalarPow<T, R>> pow(const ascalar<T, R>& op, const double& n) {
    return ascalar<T, ascalarPow<T, R>>(ascalarPow<T, R>(op.getData(), n));
};

/**
 * @brief 
 * 
 * @tparam T 
 * @tparam Op
 */

template<typename T, typename Op>
class ascalarExp {

private:

    const Op& op_;

public:

    ascalarExp(const Op& op) : op_(op) {};
    
    T getValue() const { return exp(op_.getValue()); };
    std::map<std::size_t, T> getDerivative() const { 
        std::map<std::size_t, T> der(op_.getDerivative());
        for (auto i = der.begin(); i != der.end(); i++)
            der[i->first] *= exp(op_.getValue());
        return der;
    };

};

/**
 * @brief 
 * 
 * @tparam T 
 * @tparam R1 
 * @param op
 * @return ascalar<T, ascalarExp<T, R>>
 */

template<typename T, typename R>
ascalar<T, ascalarExp<T, R>> exp(const ascalar<T, R>& op) {
    return ascalar<T, ascalarExp<T, R>>(ascalarExp<T, R>(op.getData()));
};

/**
 * @brief 
 * 
 * @tparam T 
 * @tparam Op
 */

template<typename T, typename Op>
class ascalarLog {

private:

    const Op& op_;

public:

    ascalarLog(const Op& op) : op_(op) {};
    
    T getValue() const { return log(op_.getValue()); };
    std::map<std::size_t, T> getDerivative() const { 
        assert(op_.getValue() > 0);
        std::map<std::size_t, T> der(op_.getDerivative());
        for (auto i = der.begin(); i != der.end(); i++)
            der[i->first] /= op_.getValue();
        return der;
    };

};

/**
 * @brief 
 * 
 * @tparam T 
 * @tparam R1 
 * @param op
 * @return ascalar<T, ascalarLog<T, R>>
 */

template<typename T, typename R>
ascalar<T, ascalarLog<T, R>> log(const ascalar<T, R>& op) {
    return ascalar<T, ascalarLog<T, R>>(ascalarLog<T, R>(op.getData()));
};

/**
 * @brief 
 * 
 * @tparam T 
 * @tparam Op
 */

template<typename T, typename Op>
class ascalarExpG {

private:

    const double& n_;
    const Op& op_;

public:

    ascalarExpG(const double& n, const Op& op) : n_(n), op_(op) {};
    
    T getValue() const { return pow(n_, op_.getValue()); };
    std::map<std::size_t, T> getDerivative() const { 
        assert(n_ > 0);
        std::map<std::size_t, T> der(op_.getDerivative());
        for (auto i = der.begin(); i != der.end(); i++)
            der[i->first] *= pow(n_, op_.getValue()) * log(n_);
        return der;
    };

};

/**
 * @brief 
 * 
 * @tparam T 
 * @tparam R 
 * @param n 
 * @param op 
 * @return ascalar<T, ascalarExpG<T, R>> 
 */

template<typename T, typename R>
ascalar<T, ascalarExpG<T, R>> pow(const double& n, const ascalar<T, R>& op) {
    return ascalar<T, ascalarExpG<T, R>>(ascalarExpG<T, R>(n, op.getData()));
};

#endif // ADCG_MATH_OPS_H_