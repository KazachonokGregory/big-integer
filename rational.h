#pragma once

#include "biginteger.h"

/*
 * A class that implements arithmetic for arbitrary rational numbers.
 */
class Rational {
public:
    // constructor by default: sets to 0
    Rational();
    // constructor from a 64-bit integer
    Rational(int64_t);
    // constructor from a BigInteger
    Rational(const BigInteger&);
    // constructor from numerator and denominator
    Rational(const BigInteger&, const BigInteger&);
    Rational(const Rational&) = default;
    ~Rational() = default;

    BigInteger::Sign& sign();
    BigInteger::Sign sign() const;

    void change_sign();
    bool operator ==(const Rational&) const = default;
    std::strong_ordering operator <=>(const Rational&) const;
    Rational& operator =(const Rational&) = default;
    string toString() const;
    // display as a decimal number with given precision
    string asDecimal(size_t) const;
    Rational& operator +=(const Rational&);
    Rational& operator *=(const Rational&);
    Rational& operator -=(const Rational&);
    Rational& operator /=(const Rational&);
    explicit operator double() const;
private:
    BigInteger m_numerator, m_denominator;
    // divide by the greatest common factor
    void simplify();
};

void Rational::simplify() {
    if (m_numerator.sign() == BigInteger::Sign::ZERO) {
        m_denominator = 1;
    }
    if (m_denominator.sign() == BigInteger::Sign::NEGATIVE) {
        m_denominator.change_sign();
        m_numerator.change_sign();
    }
    BigInteger factor = gcd(abs(m_numerator), m_denominator);
    m_numerator /= factor;
    m_denominator /= factor;
}

Rational::Rational():
    m_numerator(0), m_denominator(1)
{}

Rational::Rational(int64_t number):
    m_numerator(number), m_denominator(1)
{}

Rational::Rational(const BigInteger& number):
    m_numerator(number), m_denominator(1)
{}

Rational::Rational(const BigInteger& first, const BigInteger& second):
    m_numerator(first), m_denominator(second)
{
    assert(m_denominator.sign() != BigInteger::Sign::ZERO && "Error: can't divide by 0!");
    simplify();
}

BigInteger::Sign& Rational::sign() {
    return m_numerator.sign();
}

BigInteger::Sign Rational::sign() const {
    return m_numerator.sign();
}

string Rational::toString() const {
    if (m_denominator == 1_bi) {
        return m_numerator.toString();
    }
    return m_numerator.toString() + "/" + m_denominator.toString();
}

string Rational::asDecimal(size_t precision = 0) const {
    BigInteger number = m_numerator;
    for (size_t i = 0; i < precision; ++i) {
        number *= BigInteger::S_TEN;
    }
    number /= m_denominator;
    string result = number.toString();
    size_t first_digit = result[0] == '-' ? 1 : 0;
    size_t leading_zeros = 0;
    if (precision + first_digit + 1 >= result.size()) {
        leading_zeros = precision + first_digit + 1 - result.size();
    }
    result.insert(first_digit, leading_zeros, '0');
    if (precision) {
        result.insert(result.size() - precision, 1, '.');
    }

    return result;
}

Rational& Rational::operator +=(const Rational& other) {
    m_numerator = m_numerator * other.m_denominator + other.m_numerator * m_denominator;
    m_denominator *= other.m_denominator;
    simplify();
    return *this;
}

Rational& Rational::operator *=(const Rational& other) {
    m_numerator *= other.m_numerator;
    m_denominator *= other.m_denominator;
    simplify();
    return *this;
}

void Rational::change_sign() {
    m_numerator.change_sign();
}

Rational operator -(const Rational& number) {
    Rational copy = number;
    copy *= -1;
    return copy;
}

Rational& Rational::operator -=(const Rational& other) {
    *this += -other;
    return *this;
}

Rational& Rational::operator /=(const Rational& other) {
    m_numerator *= other.m_denominator;
    m_denominator *= other.m_numerator;
    simplify();
    return *this;
}

Rational operator +(const Rational& lhs, const Rational& rhs) {
    Rational copy = lhs;
    copy += rhs;
    return copy;
}

Rational operator -(Rational lhs, const Rational& rhs) {
    Rational copy = lhs;
    copy -= rhs;
    return copy;
}

Rational operator *(Rational lhs, const Rational& rhs) {
    Rational copy = lhs;
    copy *= rhs;
    return copy;
}

Rational operator /(Rational lhs, const Rational& rhs) {
    Rational copy = lhs;
    copy /= rhs;
    return copy;
}

// convert rational to double
Rational::operator double() const {
    const size_t double_precision = 30;
    string s = asDecimal(double_precision);
    return stod(s);
}

std::strong_ordering Rational::operator <=>(const Rational& other) const {
    if (sign() != other.sign()) {
        return sign() <=> other.sign();
    }
    return m_numerator * other.m_denominator <=> other.m_numerator * m_denominator;
}

