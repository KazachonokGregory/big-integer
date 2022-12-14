#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
#include <compare>
#include <cassert>
#include <cmath>
#include <complex>

using std::vector;
using std::string;

class BigInteger {
public:
    enum class Sign;

    static const int64_t s_base = 10000;
    static const int s_ten = 10;
    static const size_t s_base_length = 4; // the base is a power of 10

    BigInteger();
    BigInteger(int64_t);
    BigInteger(string);
    BigInteger(const vector<int64_t>&, BigInteger::Sign);
    BigInteger(const BigInteger&);
    ~BigInteger();

    void swap(BigInteger&);
    BigInteger& operator =(BigInteger);
    size_t size() const;
    const vector<int64_t>& digits() const; // returns the digits vector
    vector<int64_t>& digits();
    int64_t& operator[](size_t);
    int64_t operator[](size_t) const;
    BigInteger::Sign sign() const;
    BigInteger::Sign& sign();
    explicit operator bool() const;
    string toString() const;
    BigInteger& operator +=(const BigInteger&);
    BigInteger& operator -=(const BigInteger&);
    void change_sign();
    BigInteger operator *(const BigInteger&) const;
    BigInteger& operator *=(const BigInteger&);
    BigInteger& operator *=(int64_t);
    BigInteger& operator /=(const BigInteger&);
    BigInteger& operator /=(int64_t);
    BigInteger& operator %=(const BigInteger&);
    BigInteger& operator %=(int64_t);

    bool abs_less(const BigInteger&) const;

private:
    typedef std::complex<long double> cld;
    static const int64_t s_magic_constant = 400; // regulates for what sizes of data we will use FFT for multiplication
    static constexpr long double s_pi_doubled = 3.14159265358979323846L * 2;

    vector<int64_t> m_digits;
    BigInteger::Sign m_sign;

    static BigInteger::Sign sign(int64_t);
    static size_t ceil_division(size_t, size_t);
    static size_t reverse_bits(size_t, size_t);
    static void fft(size_t, cld*, cld);
    static vector<int64_t> multiply_polynomials(const vector<int64_t>&, const vector<int64_t>&);

    BigInteger fft_multiply(const BigInteger&) const;
    BigInteger quadratic_multiply(const BigInteger&) const;
    std::pair<BigInteger, BigInteger> div(const BigInteger&) const;
};

enum class BigInteger::Sign {
    POSITIVE = 1,
    ZERO = 0,
    NEGATIVE = -1
};

bool BigInteger::abs_less(const BigInteger& other) const {
    if (size() != other.size()) {
        return size() < other.size();
    }
    for (size_t i = size(); i > 0; --i) {
        if (m_digits[i - 1] != other.m_digits[i - 1]) {
            return m_digits[i - 1] < other.m_digits[i - 1];
        }
    }
    return 0;
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
        m_digits.push_back(number % s_base);
        number /= s_base;
    }
}

BigInteger::BigInteger(string s) {
    if (s == "0" || s == "-0") {
        m_sign = BigInteger::Sign::ZERO;
        m_digits = {0};
        return;
    }

    if (s[0] == '-') {
        m_sign = BigInteger::Sign::NEGATIVE;
        s.erase(0, 1);
    }
    else {
        m_sign = BigInteger::Sign::POSITIVE;
    }
    while (s[0] == '0' && s.size() > 2) {
	s.erase(0, 1);
    }
    size_t new_size = ceil_division(s.size(), s_base_length);
    m_digits.resize(new_size);
    for (size_t i = s.size(); i > 0; --i) {
        m_digits[(i - 1) / s_base_length] = m_digits[(i - 1) / s_base_length] * BigInteger::s_ten + s[s.size() - i] - '0';
    }
}

BigInteger::BigInteger(const vector<int64_t>& m_digits, BigInteger::Sign sign):
m_digits(m_digits), m_sign(sign)
{}

BigInteger::BigInteger(const BigInteger& other):
m_digits(other.m_digits), m_sign(other.m_sign)
{}

BigInteger::~BigInteger() {}

void BigInteger::swap(BigInteger& other) {
    std::swap(m_digits, other.m_digits);
    std::swap(m_sign, other.m_sign);
}

BigInteger& BigInteger::operator =(BigInteger other) {
    swap(other);
    return *this;
}

size_t BigInteger::size() const {
    return m_digits.size();
}

const vector<int64_t>& BigInteger::digits() const {
    return m_digits;
}

vector<int64_t>& BigInteger::digits() {
    return m_digits;
}

int64_t& BigInteger::operator [](size_t index) {
    return m_digits[index];
}

int64_t BigInteger::operator [](size_t index) const {
    return m_digits[index];
}

BigInteger::Sign BigInteger::sign() const {
    return m_sign;
}

BigInteger::Sign& BigInteger::sign() {
    return m_sign;
}

bool operator ==(const BigInteger& lhs, const BigInteger& rhs) {
    if (lhs.sign() != rhs.sign()) {
        return false;
    }
    return lhs.digits() == rhs.digits();
}

std::strong_ordering operator <=>(const BigInteger& lhs, const BigInteger& rhs) {
    if (lhs.sign() != rhs.sign()) {
        return lhs.sign() <=> rhs.sign();
    }
    size_t n = lhs.size(), m = rhs.size();
    if (n != m) {
        if (lhs.sign() == BigInteger::Sign::POSITIVE) {
            return n <=> m;
        }
        else {
            return m <=> n;
        }
    }
    for (size_t i = n; i > 0; --i) {
        if (lhs[i - 1] != rhs[i - 1]) {
            if (lhs.sign() == BigInteger::Sign::POSITIVE) {
                return lhs[i - 1] <=> rhs[i - 1];
            }
            else {
                return rhs[i - 1] <=> lhs[i - 1];
            }
        }
    }
    return std::strong_ordering::equivalent;
}

std::istream& operator >>(std::istream& in, BigInteger& number) {
    string s;
    in >> s;
    number = s;
    return in;
}
std::ostream& operator <<(std::ostream& out, BigInteger number) {
    if (number.sign() == BigInteger::Sign::NEGATIVE) {
        out << '-';
    }
    //assert(isdigit(number.digits().back()));
    out << number.digits().back();
    if (number.size() == 1) {
        return out;
    }
    for (size_t i = number.size() - 1; i > 0; --i) {
        out << std::setw(BigInteger::s_base_length) << std::setfill('0') << number[i - 1];
    }
    return out;
}

BigInteger::operator bool() const {
    return m_sign != BigInteger::Sign::ZERO;
}

string BigInteger::toString() const {
    string result;
    std::stringstream ss;
    ss << *this;
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

BigInteger operator -(BigInteger number) {
    number.change_sign();
    return number;
}

BigInteger& BigInteger::operator +=(const BigInteger& other) {
    if (other.m_sign == BigInteger::Sign::ZERO) {
        return *this;
    }
    if (m_sign == BigInteger::Sign::ZERO) {
        *this = other;
        return *this;
    }
    bool different_signs = (m_sign != other.m_sign);
    bool need_change_sign = different_signs & abs_less(other);
    bool carry = 0;
    size_t max_size = std::max(size(), other.size());
    for (size_t i = 0; i < max_size || carry; ++i) {
        if (i >= other.size() && carry == 0) {
            break;
        }
        if (i >= size()) {
            m_digits.push_back(0);
        }
        int64_t add = (i < other.size() ? other.m_digits[i]:0);
        m_digits[i] += (different_signs ? -add : add);
        m_digits[i] = (need_change_sign ? -m_digits[i] : m_digits[i]);
        m_digits[i] += (different_signs ? -carry : carry);
        carry = (m_digits[i] >= s_base || m_digits[i] < 0);
        if (carry) {
            m_digits[i] += (different_signs ? s_base : -s_base);
        }
    }
    while(m_digits.size() > 1 && m_digits.back() == 0) {
        m_digits.pop_back();
    }
    if (need_change_sign) {
        change_sign();
    }
    if (m_digits.size() == 1 && m_digits.back() == 0) {
        m_sign = BigInteger::Sign::ZERO;
    }
    return *this;
}

BigInteger abs(BigInteger number) {
    if(number.sign() == BigInteger::Sign::NEGATIVE) {
        number.change_sign();
    }
    return number;
}

BigInteger& BigInteger::operator -=(const BigInteger& other) {
    *this += -other;
    return *this;
}

BigInteger operator +(BigInteger lhs, const BigInteger& rhs) {
    lhs += rhs;
    return lhs;
}

BigInteger operator -(BigInteger lhs, const BigInteger& rhs) {
    lhs -= rhs;
    return lhs;
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

void BigInteger::fft(size_t n, cld* m_digits, cld q) {
    size_t logn = 0;
    while (static_cast<size_t>(1 << logn) < n) {
        ++logn;
    }

    cld* q_powers = new cld[logn];
    q_powers[0] = q;
    for (size_t i = 1; i < logn; ++i) {
        q_powers[i] = q_powers[i-1] * q_powers[i-1];
    }

    for (size_t i = 0; i < n; ++i) {
        size_t reversed_i = reverse_bits(i, logn);
        if (i < reversed_i) {
            std::swap(m_digits[i], m_digits[reversed_i]);
        }
    }

    for (size_t i = 1; i <= logn; ++i) {
        size_t block_size = (1 << (i - 1));
        cld q_current = q_powers[logn - i];
        for (size_t j = 0; j < n; j += block_size * 2) {
            cld q_degree = 1;
            cld* b = m_digits + j;
            cld* c = m_digits + j + block_size;
            for (size_t k = 0; k < block_size; ++k) {
                cld first = b[k], second = c[k];
                b[k] = first + q_degree * second;
                c[k] = first - q_degree * second;
                q_degree *= q_current;
            }
        }
    }

    delete[] q_powers;
}

vector<int64_t> BigInteger::multiply_polynomials(const vector<int64_t>& A, const vector<int64_t>& B) {
    size_t n = 1;
    while (n < A.size() + B.size()) {
        n *= 2;
    }
    cld* a = new cld[n];
    cld* b = new cld[n];
    for (size_t i = 0; i < n; ++i) {
        a[i] = (i < A.size() ? A[i] : 0);
        b[i] = (i < B.size() ? B[i] : 0);
    }
    cld q = {cosl(s_pi_doubled / n), sinl(s_pi_doubled / n)};
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
        result[i] += carry;
        carry = result[i] / s_base;
        result[i] -= carry * s_base;
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
            int64_t current = result[i+j] + (j < other.size() ? m_digits[i] * other[j] : 0ll) + carry;
            carry = current / s_base;
            result[i+j] = current - carry * s_base;
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
    if (operations_expected > s_magic_constant * s_magic_constant && max_size > s_magic_constant) {
        return fft_multiply(other);
    }
    else {
        return quadratic_multiply(other);
    }
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
        carry = cur / s_base;
        m_digits[i] = cur - carry * s_base;
    }
    while (m_digits.size() > 1 && m_digits.back() == 0) {
        m_digits.pop_back();
    }
    return *this;
}

BigInteger operator *(BigInteger lhs, int64_t rhs) {
    lhs *= rhs;
    return lhs;
}

BigInteger operator ""_bi(const char* x) {
    return BigInteger(x);
}

std::pair<BigInteger, BigInteger> BigInteger::div(const BigInteger& other) const {
    // if *this = p * other + r, this function returns {p, r}
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
        remainder = remainder * s_base + m_digits[n - 1 - i];
        if (i >= m - 1) {
            int64_t left = 0;
            int64_t right = s_base - 1;
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
            quotient = quotient * s_base + left;
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

BigInteger operator /(BigInteger lhs, const BigInteger& rhs) {
    lhs /= rhs;
    return lhs;
}

BigInteger operator %(BigInteger lhs, const BigInteger& rhs) {
    lhs %= rhs;
    return lhs;
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
        long long current = m_digits[i - 1] + carry * s_base;
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
    BigInteger result = lhs;
    result /= rhs;
    return result;
}

int64_t operator %(const BigInteger& lhs, int64_t rhs) {
    int64_t result = 0;
    for (size_t i = lhs.size(); i > 0; --i) {
        result = (result * lhs.s_base + lhs[i - 1]) % rhs;
    }
    if (lhs.sign() == BigInteger::Sign::NEGATIVE) {
        result = -result;
    }
    return result;
}

BigInteger& BigInteger::operator %=(int64_t number) {
    *this = *this % number;
    return *this;
}

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
            A.swap(B);
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

class Rational {
public:
    Rational();
    Rational(int64_t);
    Rational(const BigInteger&);
    Rational(const BigInteger&, const BigInteger&);
    Rational(const Rational&);
    ~Rational();

    BigInteger& numerator();
    const BigInteger& numerator() const;
    BigInteger& denominator();
    const BigInteger& denominator() const;
    BigInteger::Sign& sign();
    BigInteger::Sign sign() const;

    void swap(Rational& other);
    Rational& operator =(Rational);
    string toString() const;
    string asDecimal(size_t) const;
    Rational& operator +=(const Rational&);
    Rational& operator *=(const Rational&);
    Rational& operator -=(const Rational&);
    Rational& operator /=(const Rational&);
    explicit operator double() const;
private:
    BigInteger m_numerator, m_denominator;
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
    assert(m_denominator.sign() != BigInteger::Sign::ZERO);
    simplify();
}

Rational::Rational(const Rational& other):
m_numerator(other.m_numerator), m_denominator(other.m_denominator)
{}

Rational::~Rational() {}

void Rational::swap(Rational& other) {
    m_numerator.swap(other.m_numerator);
    m_denominator.swap(other.m_denominator);
}

Rational& Rational::operator =(Rational other) {
    swap(other);
    return *this;
}

BigInteger& Rational::numerator() {
    return m_numerator;
}

const BigInteger& Rational::numerator() const {
    return m_numerator;
}

BigInteger& Rational::denominator() {
    return m_denominator;
}

const BigInteger& Rational::denominator() const {
    return m_denominator;
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
        number *= BigInteger::s_ten;
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

Rational operator -(Rational number) {
    number.denominator().change_sign();
    return number;
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

Rational operator +(Rational lhs, const Rational& rhs) {
    lhs += rhs;
    return lhs;
}

Rational operator -(Rational lhs, const Rational& rhs) {
    lhs -= rhs;
    return lhs;
}

Rational operator *(Rational lhs, const Rational& rhs) {
    lhs *= rhs;
    return lhs;
}

Rational operator /(Rational lhs, const Rational& rhs) {
    lhs /= rhs;
    return lhs;
}

Rational::operator double() const {
    const size_t double_precision = 30;
    string s = asDecimal(double_precision);
    return stod(s);
}

std::strong_ordering operator <=>(const Rational& lhs, const Rational& rhs) {
    if (lhs.sign() != rhs.sign()) {
        return lhs.sign() <=> rhs.sign();
    }
    return lhs.numerator() * rhs.denominator() <=> rhs.numerator() * lhs.denominator();
}

bool operator ==(const Rational& lhs, const Rational& rhs) {
    return lhs.numerator() * rhs.denominator() == rhs.numerator() * lhs.denominator();
}

//
