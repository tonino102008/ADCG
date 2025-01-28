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

#endif // ADCG_ASCALAR_H_