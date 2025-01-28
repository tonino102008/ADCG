#ifndef ADCG_ASCALAR_H_
#define ADCG_ASCALAR_H_

#include <map>

#include "adcg_types.hpp"

/**
 * @brief 
 * 
 * @tparam T 
 */

template<typename T>
class adouble {

private:
    
    T val_;
    std::map<std::size_t, T> der_;
    const std::size_t n_;

public:

    adouble() : val_(0), der_(std::map<std::size_t, T>()), n_(0) { der_[n_] = 0; };
    adouble(const T& val) : val_(val), der_(std::map<std::size_t, T>()), n_(0) { der_[n_] = 0; };
    adouble(const T& val, const std::size_t& n) : val_(val), der_(std::map<std::size_t, T>()), n_(n) { der_[n_] = 1; };
    adouble(const adouble<T>& adouble) : val_(adouble.getValue()), der_(adouble.getDerivative()), n_(adouble.getNo()) {};

    const std::size_t& getNo() const { return n_; };
    T getValue() const { return val_; };
    T& getValue() { return val_; };
    std::map<std::size_t, T> getDerivative() const { return der_; };
    std::map<std::size_t, T>& getDerivative() { return der_; };

};

/**
 * @brief 
 * 
 * @tparam T 
 * @tparam C 
 */

template<typename T, typename C = class adouble<T>>
class ascalar {

private:
    
    C adouble_;

public:

    ascalar() : adouble_(C()) {};
    ascalar(const T& val) : adouble_(C(val)) {};
    ascalar(const T& val, const std::size_t& n) : adouble_(C(val, n)) {};
    ascalar(const C& adouble) : adouble_(adouble) {};
    template<typename T2, typename C2>
    ascalar(const ascalar<T2, C2>& ascalar) : adouble_(ascalar.getData()) {};

    const C& getData() const { return adouble_; };

    T getValue() const { return adouble_.getValue(); };
    T& getValue() { return adouble_.getValue(); };
    std::map<std::size_t, T> getDerivative() const { return adouble_.getDerivative(); };
    std::map<std::size_t, T>& getDerivative() { return adouble_.getDerivative(); };

    template<typename T2, typename C2>
    ascalar& operator= (const ascalar<T2, C2>& op) {
        adouble_.getValue() = op.getValue();
        adouble_.getDerivative() = op.getDerivative();
        return *this;
    };
    
};

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

#endif // ADCG_ASCALAR_H_