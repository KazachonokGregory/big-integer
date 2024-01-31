#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
#include <compare>
#include <cassert>
#include <cmath>
#include <complex>
#include <numbers>

using std::vector;
using std::string;

/*
 * A class that implements arithmetic for arbitrarily large integer numbers.
 */
class BigInteger {
public:
    enum class Sign;

    static const int S_TEN = 10;

    // constructor by default: sets to 0
    BigInteger();
    // constructor from a 64-bit integer
    BigInteger(int64_t);
    // constructor from a string (if the string isn't a valid number, behavior is undefined)
    explicit BigInteger(const string&);
    BigInteger(const BigInteger&) = default;
    ~BigInteger() = default;

    bool operator ==(const BigInteger&) const;
    BigInteger& operator =(const BigInteger&) = default;
    BigInteger::Sign sign() const;
    BigInteger::Sign& sign();
    // converts to bool like regular int
    explicit operator bool() const;
    string toString() const;
    BigInteger& operator +=(const BigInteger&);
    BigInteger& operator -=(const BigInteger&);
    // switches the sign to the opposite
    void change_sign();
    BigInteger operator *(const BigInteger&) const;
    BigInteger& operator *=(const BigInteger&);
    BigInteger& operator *=(int64_t);
    BigInteger& operator /=(const BigInteger&);
    BigInteger& operator /=(int64_t);
    BigInteger& operator %=(const BigInteger&);
    BigInteger& operator %=(int64_t);
    int64_t operator %(int64_t) const;

    // compares the absolute value with the other
    std::strong_ordering abs_compare(const BigInteger&) const;

    // the number of digits
    size_t size() const;

private:
    typedef std::complex<long double> cld;
    // for numbers with more digits than S_MAGIC_CONSTANT, FFT is used for multiplication
    static const int64_t S_MAGIC_CONSTANT = 400; 
    static constexpr long double S_PI = std::numbers::pi_v<long double>;

    vector<int64_t> m_digits;
    BigInteger::Sign m_sign;

    // constructor from a vector of digits between 0 and S_BASE - 1
    BigInteger(const vector<int64_t>&, BigInteger::Sign);

    static BigInteger::Sign sign(int64_t);
    // divides and rounds up
    static size_t ceil_division(size_t, size_t);
    static size_t reverse_bits(size_t, size_t);
    // the Fast Fourier Transform (https://en.wikipedia.org/wiki/Fast_Fourier_transform)
    static void fft(size_t, cld*, cld);
    static vector<int64_t> multiply_polynomials(const vector<int64_t>&, const vector<int64_t>&);

    // use FFT for multiplication
    BigInteger fft_multiply(const BigInteger&) const;
    // use the simple algorithm for multiplication
    BigInteger quadratic_multiply(const BigInteger&) const;
    // divide and return both the quotient and the remainder
    std::pair<BigInteger, BigInteger> div(const BigInteger&) const;
    // BigInteger is stored in base S_BASE
    static const int64_t S_BASE = 10000;
    // S_BASE is 10 to the power S_BASE_LENGTH
    static const size_t S_BASE_LENGTH = 4; 
};

enum class BigInteger::Sign {
    POSITIVE = 1,
    ZERO = 0,
    NEGATIVE = -1
};

std::strong_ordering BigInteger::abs_compare(const BigInteger& other) const {
    if (size() != other.size()) {
        return size() <=> other.size();
    }
    for (size_t i = size(); i > 0; --i) {
        if (m_digits[i - 1] != other.m_digits[i - 1]) {
            return m_digits[i - 1] <=> other.m_digits[i - 1];
        }
    }
    return std::strong_ordering::equivalent;
}

BigInteger::Sign BigInteger::sign(int64_t x) {
    if(x > 0) {
        return BigInteger::Sign::POSITIVE;
    }
    if(x < 0) {
        return BigInteger::Sign::NEGATIVE;
    }
    return BigInteger::Sign::ZERO;
}

size_t BigInteger::ceil_division(size_t x, size_t y) {
    return (x + y - 1) / y;
}

BigInteger::BigInteger():
    m_digits(1, 0), m_sign(BigInteger::Sign::ZERO)
{}

BigInteger::BigInteger(int64_t number):
    m_sign(sign(number)) {
        number = number < 0 ? -number : number;
        if (number == 0) {
            m_digits = {0};
            return;
        }
        while (number > 0) {
            m_digits.push_back(number % S_BASE);
            number /= S_BASE;
        }
    }

BigInteger::BigInteger(const string& s) {
    if (s == "0" || s == "-0") {
        m_sign = BigInteger::Sign::ZERO;
        m_digits = {0};
        return;
    }

    size_t start = 0;
    m_sign = BigInteger::Sign::POSITIVE;
    if (s[0] == '-') {
        m_sign = BigInteger::Sign::NEGATIVE;
        ++start;
    }
    while (s[start] == '0' && s.size() - start > 1) {
        ++start;
    }
    size_t new_size = ceil_division(s.size() - start, S_BASE_LENGTH);
    m_digits.resize(new_size);
    for (size_t i = s.size(); i > start; --i) {
        m_digits[(i - start - 1) / S_BASE_LENGTH] = m_digits[(i - start - 1) / S_BASE_LENGTH] 
            * BigInteger::S_TEN + s[s.size() + start - i] - '0';
    }
}

BigInteger::BigInteger(const vector<int64_t>& m_digits, BigInteger::Sign sign):
    m_digits(m_digits), m_sign(sign)
{}

size_t BigInteger::size() const {
    return m_digits.size();
}

BigInteger::Sign BigInteger::sign() const {
    return m_sign;
}

BigInteger::Sign& BigInteger::sign() {
    return m_sign;
}

bool BigInteger::operator ==(const BigInteger& other) const {
    if (sign() != other.sign()) {
        return false;
    }
    return m_digits == other.m_digits;
}

std::strong_ordering operator <=>(const BigInteger& lhs, const BigInteger& rhs) {
    if (lhs.sign() != rhs.sign()) {
        return lhs.sign() <=> rhs.sign();
    }
    std::strong_ordering order = lhs.abs_compare(rhs);
    if (lhs.sign() == BigInteger::Sign::NEGATIVE) {
        return 0 <=> order;
    }    
    return order;
}

std::istream& operator >>(std::istream& in, BigInteger& number) {
    string s;
    in >> s;
    number = static_cast<BigInteger>(s);
    return in;
}

std::ostream& operator <<(std::ostream& out, BigInteger number) {
    out << number.toString();
    return out;
}

BigInteger::operator bool() const {
    return m_sign != BigInteger::Sign::ZERO;
}

string BigInteger::toString() const {
    std::stringstream ss;
    if (sign() == BigInteger::Sign::NEGATIVE) {
        ss << '-';
    }
    ss << m_digits.back();
    for (size_t i = m_digits.size() - 1; i > 0; --i) {
        ss << std::setw(BigInteger::S_BASE_LENGTH) << std::setfill('0') << m_digits[i - 1];
    }    
    string result;
    ss >> result;
    return result;
}

void BigInteger::change_sign() {
    if (m_sign == BigInteger::Sign::POSITIVE) {
        m_sign = BigInteger::Sign::NEGATIVE;
    }
    else if (m_sign == BigInteger::Sign::NEGATIVE) {
        m_sign = BigInteger::Sign::POSITIVE;
    }
}

BigInteger operator -(const BigInteger& number) {
    BigInteger copy = number;
    copy.change_sign();
    return copy;
}

BigInteger& BigInteger::operator +=(const BigInteger& other) {
    if (other.m_sign == BigInteger::Sign::ZERO) {
        return *this;
    }
    if (m_sign == BigInteger::Sign::ZERO) {
        *this = other;
        return *this;
    }
    int different_signs = m_sign != other.m_sign ? -1 : 1;
    int need_change_sign = different_signs == -1 && abs_compare(other) == std::strong_ordering::less ? -1 : 1;
    bool carry = 0;
    size_t max_size = std::max(size(), other.size());
    for (size_t i = 0; i < max_size || carry; ++i) {
        if (i >= other.size() && carry == 0) {
            break;
        }
        if (i >= size()) {
            m_digits.push_back(0);
        }
        int64_t add = (i < other.size() ? other.m_digits[i] : 0);
        m_digits[i] += add * different_signs;
        m_digits[i] = m_digits[i] * need_change_sign;
        m_digits[i] += carry * different_signs;
        carry = (m_digits[i] >= S_BASE || m_digits[i] < 0);
        if (carry) {
            m_digits[i] -= S_BASE * different_signs;
        }
    }
    while(m_digits.size() > 1 && m_digits.back() == 0) {
        m_digits.pop_back();
    }
    if (need_change_sign == -1) {
        change_sign();
    }
    if (m_digits.size() == 1 && m_digits.back() == 0) {
        m_sign = BigInteger::Sign::ZERO;
    }
    return *this;
}

BigInteger abs(const BigInteger& number) {
    BigInteger copy = number;
    if(copy.sign() == BigInteger::Sign::NEGATIVE) {
        copy.change_sign();
    }
    return copy;
}

BigInteger& BigInteger::operator -=(const BigInteger& other) {
    change_sign();
    *this += other;
    change_sign();
    return *this;
}

BigInteger operator +(const BigInteger& lhs, const BigInteger& rhs) {
    BigInteger copy = lhs;
    copy += rhs;
    return copy;
}

BigInteger operator -(const BigInteger& lhs, const BigInteger& rhs) {
    BigInteger copy = lhs;
    copy -= rhs;
    return copy;
}

size_t BigInteger::reverse_bits(size_t x, size_t logn) {
    size_t answer = 0;
    for (size_t i = 0; i < logn; ++i) {
        if (x & (1 << i)) {
            answer |= (1 << (logn - 1 - i));
        }
    }
    return answer;
}

void BigInteger::fft(size_t n, cld* a, cld q) {
    size_t logn = 0;
    while (static_cast<size_t>(1 << logn) < n) {
        ++logn;
    }

    vector<cld> q_powers(logn);
    q_powers[0] = q;
    for (size_t i = 1; i < logn; ++i) {
        q_powers[i] = q_powers[i-1] * q_powers[i-1];
    }

    for (size_t i = 0; i < n; ++i) {
        size_t reversed_i = reverse_bits(i, logn);
        if (i < reversed_i) {
            std::swap(a[i], a[reversed_i]);
        }
    }

    for (size_t i = 1; i <= logn; ++i) {
        size_t block_size = (1 << (i - 1));
        cld q_current = q_powers[logn - i];
        for (size_t j = 0; j < n; j += block_size * 2) {
            cld q_degree = 1;
            cld* b = a + j;
            cld* c = a + j + block_size;
            for (size_t k = 0; k < block_size; ++k) {
                cld first = b[k], second = c[k];
                b[k] = first + q_degree * second;
                c[k] = first - q_degree * second;
                q_degree *= q_current;
            }
        }
    }
}

vector<int64_t> BigInteger::multiply_polynomials(const vector<int64_t>& A, const vector<int64_t>& B) {
    size_t n = 1;
    while (n < A.size() + B.size()) {
        n *= 2;
    }
    cld* a = new cld[n];
    cld* b = new cld[n];
    std::copy(A.begin(), A.end(), a);
    std::copy(B.begin(), B.end(), b);

    cld q = {cosl(S_PI * 2 / n), sinl(S_PI * 2 / n)};
    fft(n, a, q);
    fft(n, b, q);

    for (size_t i = 0; i < n; ++i) {
        a[i] *= b[i];
    }
    fft(n, a, conj(q));
    vector<int64_t> C;
    C.reserve(n);
    for (size_t i = 0; i < n; ++i) {
        a[i] /= n;
        C.push_back(static_cast<int64_t>(roundl(a[i].real())));
    }
    while(C.back() == 0 && C.size() > 1) {
        C.pop_back();
    }
    delete [] a;
    delete [] b;
    return C;
}

BigInteger BigInteger::fft_multiply(const BigInteger& other) const {
    BigInteger result(multiply_polynomials(m_digits, other.m_digits), m_sign);
    if (other.m_sign == BigInteger::Sign::NEGATIVE) {
        result.change_sign();
    }
    int64_t carry = 0;
    for (size_t i = 0; i < result.size() || carry; ++i) {
        if (i >= result.size()) {
            result.m_digits.push_back(0);
        }
        result.m_digits[i] += carry;
        carry = result.m_digits[i] / S_BASE;
        result.m_digits[i] -= carry * S_BASE;
    }

    return result;
}

BigInteger BigInteger::quadratic_multiply(const BigInteger& other) const {
    BigInteger result(vector<int64_t>(size() + other.size()), m_sign);
    if (other.m_sign == BigInteger::Sign::NEGATIVE) {
        result.change_sign();
    }
    for (size_t i = 0; i < size(); ++i) {
        int64_t carry = 0;
        for (size_t j = 0; j < other.size() || carry; ++j) {
            int64_t current = result.m_digits[i+j] + (j < other.size() ? m_digits[i] * other.m_digits[j] : 0ll) + carry;
            carry = current / S_BASE;
            result.m_digits[i + j] = current - carry * S_BASE;
        }
    }

    while (result.size() > 1 && result.m_digits.back() == 0) {
        result.m_digits.pop_back();
    }

    return result;
}

BigInteger BigInteger::operator *(const BigInteger& other) const {
    if (m_sign == BigInteger::Sign::ZERO || other.m_sign == BigInteger::Sign::ZERO) {
        return 0;
    }
    uint64_t operations_expected = size() * other.size();
    uint64_t max_size = static_cast<uint64_t>(std::max(size(), other.size()));
    if (operations_expected > S_MAGIC_CONSTANT * S_MAGIC_CONSTANT && max_size > S_MAGIC_CONSTANT) {
        return fft_multiply(other);
    }
    return quadratic_multiply(other);
}

BigInteger& BigInteger::operator *=(const BigInteger& other) {
    *this = *this * other;
    return *this;
}

BigInteger& BigInteger::operator *=(int64_t number) {
    if (number == 0) {
        *this = 0;
        return *this;
    }

    if (number < 0) {
        change_sign();
        number = -number;
    }

    int64_t carry = 0;
    for (size_t i = 0; i < m_digits.size() || carry; ++i) {
        if (i >= m_digits.size()) {
            m_digits.push_back(0);
        }
        int64_t cur = carry + m_digits[i] * number;
        carry = cur / S_BASE;
        m_digits[i] = cur - carry * S_BASE;
    }
    while (m_digits.size() > 1 && m_digits.back() == 0) {
        m_digits.pop_back();
    }
    return *this;
}

BigInteger operator *(const BigInteger& lhs, int64_t rhs) {
    BigInteger copy = lhs;
    copy *= rhs;
    return copy;
}

// integer literal that indicates BigInteger
BigInteger operator ""_bi(unsigned long long x) {
    return BigInteger(static_cast<int64_t>(x));
}

// string literal that indicates BigInteger
BigInteger operator ""_bi(const char* x, size_t) {
    return BigInteger(x);
}

std::pair<BigInteger, BigInteger> BigInteger::div(const BigInteger& other) const {
    assert(other != 0_bi);
    BigInteger result;
    size_t n = size(), m = other.size();
    if (n < m) {
        return {0, *this};
    }
    BigInteger remainder = 0;
    BigInteger quotient = 0;
    BigInteger other_abs = abs(other);
    for (size_t i = 0; i < n; ++i) {
        remainder = remainder * S_BASE + m_digits[n - 1 - i];
        if (i >= m - 1) {
            int64_t left = 0;
            int64_t right = S_BASE - 1;
            while (left < right) {
                int64_t middle = (left + right + 1) / 2;
                if (other_abs * middle <= remainder) {
                    left = middle;
                }
                else {
                    right = middle - 1;
                }
            }
            remainder -= other_abs * left;
            quotient = quotient * S_BASE + left;
        }
    }

    if (m_sign != other.m_sign) {
        quotient.change_sign();
    }
    if (m_sign == BigInteger::Sign::NEGATIVE) {
        remainder.change_sign();
    }

    return {quotient, remainder};
}

BigInteger operator /(const BigInteger& lhs, const BigInteger& rhs) {
    BigInteger copy = lhs;
    copy /= rhs;
    return copy;
}

BigInteger operator %(const BigInteger lhs, const BigInteger& rhs) {
    BigInteger copy = lhs;
    copy %= rhs;
    return copy;
}

BigInteger& BigInteger::operator /=(const BigInteger& other) {
    *this = this->div(other).first;
    return *this;
}

BigInteger& BigInteger::operator %=(const BigInteger& other) {
    *this = this->div(other).second;
    return *this;
}

BigInteger& operator ++(BigInteger& number) {
    return number += 1;
}

BigInteger operator ++(BigInteger& number, int) {
    BigInteger previous_value = number;
    number += 1;
    return previous_value;
}

BigInteger& operator --(BigInteger& number) {
    return number -= 1;
}

BigInteger operator --(BigInteger& number, int) {
    BigInteger previous_value = number;
    number -= 1;
    return previous_value;
}

BigInteger& BigInteger::operator /=(int64_t number) {
    assert(number != 0);

    int64_t carry = 0;
    for (size_t i = m_digits.size(); i > 0; --i) {
        long long current = m_digits[i - 1] + carry * S_BASE;
        m_digits[i - 1] = current / number;
        carry = current - m_digits[i - 1] * number;
    }
    while (m_digits.size() > 1 && m_digits.back() == 0) {
        m_digits.pop_back();
    }
    if (number < 0) {
        change_sign();
    }
    return *this;
}

BigInteger operator /(const BigInteger& lhs, int64_t rhs) {
    BigInteger copy = lhs;
    copy /= rhs;
    return copy;
}

int64_t BigInteger::operator %(int64_t rhs) const {
    int64_t result = 0;
    for (size_t i = size(); i > 0; --i) {
        result = (result * S_BASE + m_digits[i - 1]) % rhs;
    }
    if (sign() == BigInteger::Sign::NEGATIVE) {
        result = -result;
    }
    return result;
}

BigInteger& BigInteger::operator %=(int64_t number) {
    *this = *this % number;
    return *this;
}

// compute the greatest common divisor using binary Euclidean algorithm
BigInteger gcd(BigInteger A, BigInteger B) {
    if (A == 0_bi) {
        return B;
    }
    if (B == 0_bi) {
        return A;
    }

    BigInteger power_2 = 1;
    while (A % 2 == 0 && B % 2 == 0) {
        A /= 2;
        B /= 2;
        power_2 *= 2;
    }
    while (A % 2 == 0) {
        A /= 2;
    }
    while (B % 2 == 0) {
        B /= 2;
    }

    while (true) {
        if (B < A) {
            std::swap(A, B);
        }
        B -= A;
        if (B == 0_bi) {
            return A * power_2;
        }
        while (B % 2 == 0) {
            B /= 2;
        }
    }
}

//
