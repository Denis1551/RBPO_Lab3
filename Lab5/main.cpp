import Math;
#include <iostream>
#include <cmath>
#include <stdexcept>
using namespace Math;

Complex sin(Complex a) {
    return Complex(sin(a.Re()) * cosh(a.Im()), cos(a.Re()) * sinh(a.Im()));
}

Rational sin(Rational drob) {
    double sin_val = sin(static_cast<double>(drob));
    int int_part = std::round(sin_val / 1e-6);
    return Rational(int_part, 1000000);
}

Complex f(const Complex& z) {
    Complex a(0.5);
    return 2 * z + sin(z - a);
}
Rational f(const Rational& r) {
    const Rational a(1, 2);
    return r * 2 + sin(r - a);
}

double f(double x) {
    return 2 * x + sin(x - 0.5);
}

int main() {
    setlocale(LC_ALL, "ru");
    double re, im, real;
    int num, den;
    char continue_input;

    do {
        std::cout << "Введите действительную и мнимую части комплексного числа (через пробел): ";
        std::cin >> re >> im;
        Complex z(re, im);
        std::cout << "Результат для комплексного числа: " << f(z) << "\n\n";

        do {
            std::cout << "Введите числитель и знаменатель рационального числа (через пробел): ";
            std::cin >> num >> den;

            try {
                Rational r(num, den);
                Rational res = f(r);
                std::cout << "Результат для рационального числа: " << static_cast<double>(res) << "\n\n";
                break; 
            }
            catch (const std::invalid_argument& e) {
                std::cout << "Ошибка: " << e.what() << std::endl;
                std::cin.clear(); 
                std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            }
        } while (true);

        std::cout << "Введите вещественное число: ";
        std::cin >> real;
        std::cout << "Результат для вещественного числа: " << f(real) << "\n\n";

        std::cout << "Хотите ввести другие значения? (y/n): ";
        std::cin >> continue_input;
    } while (continue_input == 'y' || continue_input == 'Y');

    system("pause");
    return 0;
}