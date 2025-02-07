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

/**
 * @brief 
 * 
 * @tparam T 
 * @tparam Op1 
 */

template<typename T, typename Op1>
class ascalarAddC {

private:

    const Op1& op1_;
    const double& n_;

public:

    ascalarAddC(const Op1& op1, const double& n) : op1_(op1), n_(n) {};
    
    T getValue() const { return op1_.getValue() + n_; };
    std::map<std::size_t, T> getDerivative() const { 
        std::map<std::size_t, T> der(op1_.getDerivative());
        return der; 
    };

};

/**
 * @brief 
 * 
 * @tparam T 
 * @tparam R1 
 * @param op1 
 * @return ascalar<T, ascalarAddC<T, R1>>
 */

template<typename T, typename R1>
ascalar<T, ascalarAddC<T, R1>> operator+(const ascalar<T, R1>& op1, const double& n) {
    return ascalar<T, ascalarAddC<T, R1>>(ascalarAddC<T, R1>(op1.getData(), n));
};

/**
 * @brief 
 * 
 * @tparam T 
 * @tparam R1 
 * @param n 
 * @param op1 
 * @return ascalar<T, ascalarAddC<T, R1>> 
 */

template<typename T, typename R1>
ascalar<T, ascalarAddC<T, R1>> operator+(const double& n, const ascalar<T, R1>& op1) {
    return ascalar<T, ascalarAddC<T, R1>>(ascalarAddC<T, R1>(op1.getData(), n));
};

/**
 * @brief 
 * 
 * @tparam T 
 * @tparam Op1 
 */

template<typename T, typename Op1>
class ascalarSubC {

private:

    const Op1& op1_;
    const double& n_;

public:

    ascalarSubC(const Op1& op1, const double& n) : op1_(op1), n_(n) {};
    
    T getValue() const { return n_ - op1_.getValue(); };
    std::map<std::size_t, T> getDerivative() const { 
        std::map<std::size_t, T> der(op1_.getDerivative());
        for (auto i = der.begin(); i != der.end(); i++)
            der[i->first] *= -1;
        return der; 
    };

};

/**
 * @brief 
 * 
 * @tparam T 
 * @tparam R1 
 * @param n 
 * @param op1 
 * @return ascalar<T, ascalarSubC<T, R1>> 
 */

template<typename T, typename R1>
ascalar<T, ascalarSubC<T, R1>> operator-(const double& n, const ascalar<T, R1>& op1) {
    return ascalar<T, ascalarSubC<T, R1>>(ascalarSubC<T, R1>(op1.getData(), n));
};

/**
 * @brief 
 * 
 * @tparam T 
 * @tparam Op1 
 */

template<typename T, typename Op1>
class ascalarSubC1 {

private:

    const Op1& op1_;
    const double& n_;

public:

    ascalarSubC1(const Op1& op1, const double& n) : op1_(op1), n_(n) {};
    
    T getValue() const { return op1_.getValue() - n_; };
    std::map<std::size_t, T> getDerivative() const { 
        std::map<std::size_t, T> der(op1_.getDerivative());
        return der; 
    };

};

/**
 * @brief 
 * 
 * @tparam T 
 * @tparam R1 
 * @param n 
 * @param op1 
 * @return ascalar<T, ascalarSubC1<T, R1>> 
 */

template<typename T, typename R1>
ascalar<T, ascalarSubC1<T, R1>> operator-(const ascalar<T, R1>& op1, const double& n) {
    return ascalar<T, ascalarSubC1<T, R1>>(ascalarSubC1<T, R1>(op1.getData(), n));
};

template<typename T, typename Op1>
class ascalarSubC2 {

private:

    const Op1& op1_;

public:

    ascalarSubC2(const Op1& op1) : op1_(op1) {};
    
    T getValue() const { return -op1_.getValue(); };
    std::map<std::size_t, T> getDerivative() const { 
        std::map<std::size_t, T> der(op1_.getDerivative());
        for (auto i = der.begin(); i != der.end(); i++)
            der[i->first] *= -1;
        return der; 
    };

};

/**
 * @brief 
 * 
 * @tparam T 
 * @tparam R1 
 * @param n 
 * @param op1 
 * @return ascalar<T, ascalarSubC2<T, R1>> 
 */

template<typename T, typename R1>
ascalar<T, ascalarSubC2<T, R1>> operator-(const ascalar<T, R1>& op1) {
    return ascalar<T, ascalarSubC2<T, R1>>(ascalarSubC2<T, R1>(op1.getData()));
};

/**
 * @brief 
 * 
 * @tparam T 
 * @tparam Op1 
 */

template<typename T, typename Op1>
class ascalarMulC {

private:

    const Op1& op1_;
    const double& n_;

public:

    ascalarMulC(const Op1& op1, const double& n) : op1_(op1), n_(n) {};
    
    T getValue() const { return op1_.getValue() * n_; };
    std::map<std::size_t, T> getDerivative() const { 
        std::map<std::size_t, T> der(op1_.getDerivative());
        for (auto i = der.begin(); i != der.end(); i++)
            der[i->first] *= n_;
        return der; 
    };

};

/**
 * @brief 
 * 
 * @tparam T 
 * @tparam R1 
 * @param op1 
 * @return ascalar<T, ascalarMulC<T, R1>>
 */

template<typename T, typename R1>
ascalar<T, ascalarMulC<T, R1>> operator*(const ascalar<T, R1>& op1, const double& n) {
    return ascalar<T, ascalarMulC<T, R1>>(ascalarMulC<T, R1>(op1.getData(), n));
};

/**
 * @brief 
 * 
 * @tparam T 
 * @tparam R1 
 * @param n 
 * @param op1 
 * @return ascalar<T, ascalarMulC<T, R1>> 
 */

template<typename T, typename R1>
ascalar<T, ascalarMulC<T, R1>> operator*(const double& n, const ascalar<T, R1>& op1) {
    return ascalar<T, ascalarMulC<T, R1>>(ascalarMulC<T, R1>(op1.getData(), n));
};

/**
 * @brief 
 * 
 * @tparam T 
 * @tparam Op1 
 */

template<typename T, typename Op1>
class ascalarDivC {

private:

    const Op1& op1_;
    const double& n_;

public:

    ascalarDivC(const Op1& op1, const double& n) : op1_(op1), n_(n) {};
    
    T getValue() const { return op1_.getValue() / n_; };
    std::map<std::size_t, T> getDerivative() const { 
        std::map<std::size_t, T> der(op1_.getDerivative());
        for (auto i = der.begin(); i != der.end(); i++)
            der[i->first] /= n_;
        return der; 
    };

};

/**
 * @brief 
 * 
 * @tparam T 
 * @tparam R1 
 * @param op1 
 * @return ascalar<T, ascalarDivC<T, R1>>
 */

template<typename T, typename R1>
ascalar<T, ascalarDivC<T, R1>> operator/(const ascalar<T, R1>& op1, const double& n) {
    return ascalar<T, ascalarDivC<T, R1>>(ascalarDivC<T, R1>(op1.getData(), n));
};

template<typename T, typename Op1>
class ascalarDivC1 {

private:

    const Op1& op1_;
    const double& n_;

public:

    ascalarDivC1(const Op1& op1, const double& n) : op1_(op1), n_(n) {};
    
    T getValue() const { return n_ / op1_.getValue(); };
    std::map<std::size_t, T> getDerivative() const { 
        std::map<std::size_t, T> der(op1_.getDerivative());
        for (auto i = der.begin(); i != der.end(); i++)
            der[i->first] *= -n_ / (op1_.getValue() * op1_.getValue());
        return der; 
    };

};

/**
 * @brief 
 * 
 * @tparam T 
 * @tparam R1 
 * @param n 
 * @param op1 
 * @return ascalar<T, ascalarDivC1<T, R1>> 
 */

template<typename T, typename R1>
ascalar<T, ascalarDivC1<T, R1>> operator/(const double& n, const ascalar<T, R1>& op1) {
    return ascalar<T, ascalarDivC1<T, R1>>(ascalarDivC1<T, R1>(op1.getData(), n));
};

#endif // ADCG_BINARY_OPS_H_