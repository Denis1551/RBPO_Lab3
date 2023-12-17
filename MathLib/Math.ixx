#include <cmath>
#include <iostream>
#include <cstdlib>
#define _USE_MATH_DEFINES
#include <math.h>
export module Math;

export namespace Math {

    int FindGreatestCommonDivisor(int a, int b) {
        while (b != 0) {
            int temp = b;
            b = a % b;
            a = temp;
        }
        return std::abs(a);
    }

    int FindLeastCommonMultiple(int x, int y) {
        int gcd = FindGreatestCommonDivisor(x, y);
        return std::abs(x * y) / gcd;
    }

	class Complex {
	public:
		Complex() 
            : r(0), argm(0) {}

		Complex(double re, double im)
			: r(std::sqrt(re * re + im * im)), argm(std::atan2(im, re)) {}

        Complex(double r)
            : r(r), argm(0) {}

		static Complex FromExponentialForm(double r, double argm) {
			return Complex(r * cos(argm), r * sin(argm));
		}

		static Complex FromAlgebraicForm(double re, double im) {
			return Complex(re, im);
		}

		double Re() const {
			return r * cos(argm);
		}

		double Im() const {
			return r * sin(argm);
		}

		double Mod() const {
			return r;
		}

		double Arg() const {
			return argm;
		}

        explicit operator double() const {
            return r * cos(argm);
        }

        Complex operator-() const {
            return Complex(-Re(), -Im());
        }

        Complex& operator=(const Complex& other) {
            r = other.r;
            argm = other.argm;
            return *this;
        }

        Complex& operator++() {
            *this = Complex(r * cos(argm) + 1, r * sin(argm));
            return *this;
        }

        Complex operator++(int) {
            Complex temp = *this;
            ++(*this); 
            return temp;
        }

        Complex& operator--() {
            *this = Complex(r * cos(argm) - 1, r * sin(argm)); 
            return *this;
        }

        Complex operator--(int) {
            Complex temp = *this;
            --(*this); 
            return temp;
        }

        Complex& operator+=(const Complex& other) {
            double re1 = r * cos(argm), im1 = r * sin(argm);
            double re2 = other.r * cos(other.argm), im2 = other.r * sin(other.argm);
            *this = Complex(re1 + re2, im1 + im2);
            return *this;
        }

        Complex& operator-=(const Complex& other) {
            double re1 = r * cos(argm), im1 = r * sin(argm);
            double re2 = other.r * cos(other.argm), im2 = other.r * sin(other.argm);
            *this = Complex(re1 - re2, im1 - im2); 
            return *this;
        }

        Complex& operator*=(const Complex& other) {
            r *= other.r;
            argm += other.argm;
            return *this;
        }

        Complex& operator/=(const Complex& other) {
            r /= other.r;
            argm -= other.argm;
            return *this;
        }

        friend Complex operator+(const Complex& lhs, const Complex& rhs);
        friend Complex operator-(const Complex& lhs, const Complex& rhs);
        friend Complex operator*(const Complex& lhs, const Complex& rhs);
        friend Complex operator/(const Complex& lhs, const Complex& rhs);
        friend Complex operator"" i(long double im);
        friend Complex operator"" i(unsigned long long im);
        friend std::ostream& operator<<(std::ostream& out, const Complex& c);
	private:
		double r, argm;
	};

    Complex operator+(const Complex& lhs, const Complex& rhs) {
        double re = lhs.Re() + rhs.Re();
        double im = lhs.Im() + rhs.Im();
        if (std::abs(re) < 1e-6) re = 0.0;
        if (std::abs(im) < 1e-6) im = 0.0;
        return Complex(re, im);
    }

    Complex operator-(const Complex& lhs, const Complex& rhs) {
        return Complex(lhs.Re() - rhs.Re(), lhs.Im() - rhs.Im());
    }

    Complex operator*(const Complex& lhs, const Complex& rhs) {
        double r = lhs.Mod() * rhs.Mod();
        double argm = lhs.Arg() + rhs.Arg();
        return Complex::FromExponentialForm(r, argm);
    }

    Complex operator/(const Complex& lhs, const Complex& rhs) {
        double r = lhs.Mod() / rhs.Mod();
        double argm = lhs.Arg() - rhs.Arg();
        return Complex::FromExponentialForm(r, argm);
    }

    Complex operator"" i(long double im) {
        return Complex(0, im);
    }

    Complex operator"" i(unsigned long long im) {
        return Complex(0, static_cast<long double>(im));
    }

    std::ostream& operator<<(std::ostream& out, const Complex& c) {
        out << c.Re(); 
        if (c.Im() == 0) {
            return out; 
        }
        out << (c.Im() > 0 ? " + " : " - ") << std::abs(c.Im()) << "i";
        return out;
    }


	class Rational {
    public: 
        Rational() : nominator(0), denominator(1) {}

        Rational(int nominator, int denominator) : nominator(nominator), denominator(denominator) {
            if (denominator == 0) {
                throw std::invalid_argument("Знаменатель не может быть равен 0.");
            }
            normalize();
        }

        Rational(int number) : nominator(number), denominator(1) {}

        int Nominator() const {
            return nominator;
        }

        int Denominator() const {
            return denominator;
        }

        explicit operator double() const {
            return static_cast<double>(nominator) / denominator;
        }

        Rational operator-() const {
            return Rational(-nominator, denominator);
        }

        Rational& operator++() {
            nominator += denominator;
            return *this;
        }

        Rational operator++(int) {
            Rational temp(*this);
            ++(*this);
            return temp;
        }

        Rational& operator--() {
            nominator -= denominator;
            return *this;
        }

        Rational operator--(int) {
            Rational temp(*this);
            --(*this);
            return temp;
        }

        Rational& operator+=(const Rational& rhs) {
            nominator = nominator * rhs.denominator + rhs.nominator * denominator;
            denominator = denominator * rhs.denominator;
            normalize();
            return *this;
        }

        Rational& operator-=(const Rational& rhs) {
            nominator = nominator * rhs.denominator - rhs.nominator * denominator;
            denominator = denominator * rhs.denominator;
            normalize();
            return *this;
        }

        Rational& operator*=(const Rational& rhs) {
            nominator *= rhs.nominator;
            denominator *= rhs.denominator;
            normalize();
            return *this;
        }

        Rational& operator/=(const Rational& rhs) {
            nominator *= rhs.denominator;
            denominator *= rhs.nominator;
            normalize();
            return *this;
        }

        friend Rational operator+(const Rational& lhs, const Rational& rhs);
        friend Rational operator-(const Rational& lhs, const Rational& rhs);
        friend Rational operator*(const Rational& lhs, const Rational& rhs);
        friend Rational operator/(const Rational& lhs, const Rational& rhs);
        friend bool operator==(const Rational& lhs, const Rational& rhs);
        friend bool operator<(const Rational& lhs, const Rational& rhs);
        friend bool operator<=(const Rational& lhs, const Rational& rhs);
        friend bool operator>(const Rational& lhs, const Rational& rhs);
        friend bool operator>=(const Rational& lhs, const Rational& rhs);
        friend std::ostream& operator<<(std::ostream& os, const Rational& r);

        private:
            int nominator, denominator;

            void normalize() {
                int gcd = FindGreatestCommonDivisor(nominator, denominator);
                nominator /= gcd;
                denominator /= gcd;
                if (denominator < 0) {
                    nominator = -nominator;
                    denominator = -denominator;
                }
            }
	};

    Rational operator+(const Rational& lhs, const Rational& rhs) {
        int newNominator = lhs.nominator * rhs.denominator + rhs.nominator * lhs.denominator;
        int newDenominator = lhs.denominator * rhs.denominator;
        return Rational(newNominator, newDenominator);
    }

    Rational operator-(const Rational& lhs, const Rational& rhs) {
        int newNominator = lhs.nominator * rhs.denominator - rhs.nominator * lhs.denominator;
        int newDenominator = lhs.denominator * rhs.denominator;
        return Rational(newNominator, newDenominator);
    }

    Rational operator*(const Rational& lhs, const Rational& rhs) {
        int newNominator = lhs.nominator * rhs.nominator;
        int newDenominator = lhs.denominator * rhs.denominator;
        return Rational(newNominator, newDenominator);
    }

    Rational operator/(const Rational& lhs, const Rational& rhs) {
        int newNominator = lhs.nominator * rhs.denominator;
        int newDenominator = lhs.denominator * rhs.nominator;
        return Rational(newNominator, newDenominator);
    }


    bool operator==(const Rational& lhs, const Rational& rhs) {
        return lhs.nominator * rhs.denominator == lhs.denominator * rhs.nominator;
    }


    bool operator<(const Rational& lhs, const Rational& rhs) {
        return lhs.nominator * rhs.denominator < lhs.denominator * rhs.nominator;
    }

    bool operator<=(const Rational& lhs, const Rational& rhs) {
        return lhs.nominator * rhs.denominator <= lhs.denominator * rhs.nominator;
    }

    bool operator>(const Rational& lhs, const Rational& rhs) {
        return rhs < lhs;
    }

    bool operator>=(const Rational& lhs, const Rational& rhs) {
        return rhs <= lhs;
    }

    std::ostream& operator<<(std::ostream& os, const Rational& r) {
        os << r.nominator;
        if (r.denominator != 1) {
            os << '/' << r.denominator;
        }
        return os;
    }
}