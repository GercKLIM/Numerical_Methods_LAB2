// ### Лабораторная работа №2 ###
// "Численное решение краевых задач для одномерного уравнения теплопроводности"
// Authors:
// @GercKLIM - Климов Олег ФН2-61б
// @Ship-Vano - Шаманов Иван ФН2-61б
//

#include <iostream>
#include "algebra.h"    // Алгебра векторов и матриц
#include "TESTS.cpp"    // Класс условий тестов
#include "SolvePDE.cpp" // Методы решения PDE

int main() {

    // Проверка работы алгебры
//    std::vector<std::vector<double>> a = create_identity_matrix<double>(2); // Создание единичной матрицы
//    std::vector<std::vector<double>> b = {{2., 2.}, {2., 2.}};
//    std::cout << a * b - a << std::endl; // Операции над матрицами и вывод в консоль
//    std::cout << "Complete!" << std::endl;
//    return 0;

    // Создание первого теста
    PDE_data test1;
    test1.c = 1;
    test1.rho = 1;
    test1.L = 10;
    test1.T = 10;
    test1.set_K([](double x) { return x; });
    test1.set_G_left([](double x) { return x; });
    test1.set_G_right([](double x) { return x; });
    //test1.show(); // Вывод информации о тесте

    // Тест: Вариант 5
    PDE_data test5;
    test5.c = 2;
    test5.rho = 0.25;
    test5.L = 1;
    test5.t0 = 0.5;
    test5.T = 10;
    test5.u0 = 0.2;
    test5.set_K([](double x) {

        double x1 = 0.5, x2 = 2./3.;
        double k1 = 2. , k2 = 0.5;
        double L = 1;

        if (x <= x1) {
            return k1;

        } else if (x < x2) {
            return (k1 * ((x - x2) / (x1 - x2)) + k2 * ((x - x1) / (x2 - x1)));

        } else if (x <= L) {
            return k2;

        } else {
            return 0.;

        }
    });

    test5.set_G_left([](double x) { return 0.2; });
    test5.set_G_right([](double x) { return 10; });

    /* Случай 1. */
    //ExplicitScheme(2, 1, 1, test1);
    SolvePDE_1(test5, 0.01, 0.01, 0.1, "ExpScheme_test5");




    std::cout << std::endl << "Complete!" << std::endl;
}
