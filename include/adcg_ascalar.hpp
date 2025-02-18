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

public:

    adouble() : val_(0), der_(std::map<std::size_t, T>()) {};
    adouble(const T& val) : val_(val), der_(std::map<std::size_t, T>()) {};
    adouble(const T& val, const std::map<std::size_t, T>& der) : val_(val), der_(der) {};
    adouble(const T& val, const std::size_t& n) : val_(val), der_(std::map<std::size_t, T>()) { der_[n] = 1; };
    adouble(const adouble<T>& adouble) : val_(adouble.getValue()), der_(adouble.getDerivative()) {};

    T getValue() const { return val_; };
    T& getValue() { return val_; };
    std::map<std::size_t, T> getDerivative() const { return der_; };
    std::map<std::size_t, T>& getDerivative() { return der_; };

    template<typename T2>
    adouble& operator+= (const adouble<T2>& op) {
        val_ += op.getValue();
        std::map<std::size_t, T> der2 = op.getDerivative();
        for (auto i = der2.begin(); i != der2.end(); i++)
            der_[i->first] += i->second;
        return *this;
    };

    template<typename T2>
    adouble& operator-= (const adouble<T2>& op) {
        val_ -= op.getValue();
        std::map<std::size_t, T> der2 = op.getDerivative();
        for (auto i = der2.begin(); i != der2.end(); i++)
            der_[i->first] -= i->second;
        return *this;
    };

    template<typename T2>
    adouble& operator*= (const T2& op) {
        std::map<std::size_t, T> der = der_;
        std::map<std::size_t, T> der2 = op.getDerivative();
        for (auto i = der_.begin(); i != der_.end(); i++)
            der_[i->first] = val_ * der2[i->first] + op.getValue() * i->second;
        for (auto i = der2.begin(); i != der2.end(); i++)
            der_[i->first] = val_ * i->second + op.getValue() * der[i->first];
        val_ *= op.getValue();
        return *this;
    };

    template<typename T2>
    adouble& operator/= (const T2& op) {
        std::map<std::size_t, T> der = der_;
        std::map<std::size_t, T> der2 = op.getDerivative();
        for (auto i = der_.begin(); i != der_.end(); i++)
            der_[i->first] = (i->second * op.getValue() - val_ * der2[i->first]) / (op.getValue() * op.getValue());
        for (auto i = der2.begin(); i != der2.end(); i++)
            der_[i->first] = (der[i->first] * op.getValue() - val_ * i->second) / (op.getValue() * op.getValue());
        val_ /= op.getValue();
        return *this;
    };

    adouble& operator*= (const double& n) {
        val_ *= n;
        for (auto i = der_.begin(); i != der_.end(); i++)
            der_[i->first] *= n;
        return *this;
    };

    adouble& operator/= (const double& n) {
        val_ /= n;
        for (auto i = der_.begin(); i != der_.end(); i++)
            der_[i->first] /= n;
        return *this;
    };

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
    ascalar(const ascalar<T2, C2>& ascalar) : adouble_(C(ascalar.getValue(), ascalar.getDerivative())) {};

    const C& getData() const { return adouble_; };

    T getValue() const { return adouble_.getValue(); };
    T& getValue() { return adouble_.getValue(); };
    std::map<std::size_t, T> getDerivative() const { return adouble_.getDerivative(); };
    std::map<std::size_t, T>& getDerivative() { return adouble_.getDerivative(); };

    template<typename T2, typename C2>
    ascalar& operator+= (const ascalar<T2, C2>& op) {
        adouble_ += op.getData();
        return *this;
    };

    template<typename T2, typename C2>
    ascalar& operator-= (const ascalar<T2, C2>& op) {
        adouble_ -= op.getData();
        return *this;
    };

    template<typename T2, typename C2>
    ascalar& operator*= (const ascalar<T2, C2>& op) {
        adouble_ *= op.getData();
        return *this;
    };

    template<typename T2, typename C2>
    ascalar& operator/= (const ascalar<T2, C2>& op) {
        adouble_ /= op.getData();
        return *this;
    };

    ascalar& operator*= (const double& n) {
        adouble_ *= n;
        return *this;
    };

    ascalar& operator/= (const double& n) {
        adouble_ /= n;
        return *this;
    };

    template<typename T2, typename C2>
    ascalar& operator= (const ascalar<T2, C2>& op) {
        adouble_.getValue() = op.getValue();
        adouble_.getDerivative() = op.getDerivative();
        return *this;
    };
    
};

template<std::size_t n, std::size_t m>
void densifyDerivative(const std::array<std::map<std::size_t, real>, n>& spJ, std::array<real, m>& J) {
    std::size_t k = m / n;
    for (std::size_t i = 0; i < n; i++)
        for (auto j = spJ[i].begin(); j != spJ[i].end(); j++)
            J[j->first + i * k] = j->second;
}

#endif // ADCG_ASCALAR_H_