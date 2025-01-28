#ifndef ADCG_BINARY_OPS_H_
#define ADCG_BINARY_OPS_H_

#include "adcg_types.hpp"
#include "adcg_ascalar.hpp"

/**
 * @brief 
 * 
 * @tparam T 
 * @tparam Op1 
 * @tparam Op2 
 */

template<typename T, typename Op1, typename Op2>
class ascalarAdd {

private:

    const Op1& op1_;
    const Op2& op2_;

public:

    ascalarAdd(const Op1& op1, const Op2& op2) : op1_(op1), op2_(op2) {};
    
    T getValue() const { return op1_.getValue() + op2_.getValue(); };
    std::map<std::size_t, T> getDerivative() const { 
        std::map<std::size_t, T> der(op1_.getDerivative());
        std::map<std::size_t, T> op2 = op2_.getDerivative();
        for (auto i = op2.begin(); i != op2.end(); i++)
            der[i->first] += i->second;
        return der; 
    };

};

/**
 * @brief 
 * 
 * @tparam T 
 * @tparam R1 
 * @tparam R2 
 * @param op1 
 * @param op2 
 * @return ascalar<T, ascalarAdd<T, R1, R2>>
 */

template<typename T, typename R1, typename R2>
ascalar<T, ascalarAdd<T, R1, R2>> operator+(const ascalar<T, R1>& op1, const ascalar<T, R2>& op2) {
    return ascalar<T, ascalarAdd<T, R1, R2>>(ascalarAdd<T, R1, R2>(op1.getData(), op2.getData()));
};

/**
 * @brief 
 * 
 * @tparam T 
 * @tparam Op1 
 * @tparam Op2 
 */

template<typename T, typename Op1, typename Op2>
class ascalarSub {

private:

    const Op1& op1_;
    const Op2& op2_;

public:

    ascalarSub(const Op1& op1, const Op2& op2) : op1_(op1), op2_(op2) {};
    
    T getValue() const { return op1_.getValue() - op2_.getValue(); };
    std::map<std::size_t, T> getDerivative() const { 
        std::map<std::size_t, T> der(op1_.getDerivative());
        std::map<std::size_t, T> op2 = op2_.getDerivative();
        for (auto i = op2.begin(); i != op2.end(); i++)
            der[i->first] -= i->second;
        return der; 
    };

};

/**
 * @brief 
 * 
 * @tparam T 
 * @tparam R1 
 * @tparam R2 
 * @param op1 
 * @param op2 
 * @return ascalar<T, ascalarSub<T, R1, R2>>
 */

template<typename T, typename R1, typename R2>
ascalar<T, ascalarSub<T, R1, R2>> operator-(const ascalar<T, R1>& op1, const ascalar<T, R2>& op2) {
    return ascalar<T, ascalarSub<T, R1, R2>>(ascalarSub<T, R1, R2>(op1.getData(), op2.getData()));
};

/**
 * @brief 
 * 
 * @tparam T 
 * @tparam Op1 
 * @tparam Op2 
 */

template<typename T, typename Op1, typename Op2>
class ascalarMul {

private:

    const Op1& op1_;
    const Op2& op2_;

public:

    ascalarMul(const Op1& op1, const Op2& op2) : op1_(op1), op2_(op2) {};

    T getValue() const { return op1_.getValue() * op2_.getValue(); };
    std::map<std::size_t, T> getDerivative() const { 
        std::map<std::size_t, T> der;
        std::map<std::size_t, T> op1 = op1_.getDerivative();
        std::map<std::size_t, T> op2 = op2_.getDerivative();
        for (auto i = op1.begin(); i != op1.end(); i++)
            der[i->first] = op1_.getValue() * op2[i->first] + op2_.getValue() * i->second;
        for (auto i = op2.begin(); i != op2.end(); i++)
            der[i->first] = op1_.getValue() * i->second + op2_.getValue() * op1[i->first];
        return der;
    };

};

/**
 * @brief 
 * 
 * @tparam T 
 * @tparam R1 
 * @tparam R2 
 * @param op1 
 * @param op2 
 * @return ascalar<T, ascalarMul<T, R1, R2>> 
 */

template<typename T, typename R1, typename R2>
ascalar<T, ascalarMul<T, R1, R2>> operator*(const ascalar<T, R1>& op1, const ascalar<T, R2>& op2) {
    return ascalar<T, ascalarMul<T, R1, R2>>(ascalarMul<T, R1, R2>(op1.getData(), op2.getData()));
};

/**
 * @brief 
 * 
 * @tparam T 
 * @tparam Op1 
 * @tparam Op2 
 */

template<typename T, typename Op1, typename Op2>
class ascalarDiv {

private:

    const Op1& op1_;
    const Op2& op2_;

public:

    ascalarDiv(const Op1& op1, const Op2& op2) : op1_(op1), op2_(op2) {};

    T getValue() const { return op1_.getValue() / op2_.getValue(); };
    std::map<std::size_t, T> getDerivative() const { 
        std::map<std::size_t, T> der;
        std::map<std::size_t, T> op1 = op1_.getDerivative();
        std::map<std::size_t, T> op2 = op2_.getDerivative();
        for (auto i = op1.begin(); i != op1.end(); i++)
            der[i->first] = (i->second * op2_.getValue() - op1_.getValue() * op2[i->first]) / (op2_.getValue() * op2_.getValue());
        for (auto i = op2.begin(); i != op2.end(); i++)
            der[i->first] = (op1[i->first] * op2_.getValue() - op1_.getValue() * i->second) / (op2_.getValue() * op2_.getValue());
        return der;
    };

};

/**
 * @brief 
 * 
 * @tparam T 
 * @tparam R1 
 * @tparam R2 
 * @param op1 
 * @param op2 
 * @return ascalar<T, ascalarDiv<T, R1, R2>> 
 */

template<typename T, typename R1, typename R2>
ascalar<T, ascalarDiv<T, R1, R2>> operator/(const ascalar<T, R1>& op1, const ascalar<T, R2>& op2) {
    return ascalar<T, ascalarDiv<T, R1, R2>>(ascalarDiv<T, R1, R2>(op1.getData(), op2.getData()));
};

#endif // ADCG_BINARY_OPS_H_