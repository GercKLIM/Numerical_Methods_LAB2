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
#include "Config.h"     // Константы для тестов

int main() {
    // Тест 1: Медь, фиксированная температура на концах
    PDE_data test1;
    test1.c = COPPER_C;
    test1.rho = COPPER_RHO;
    test1.L = 50.;
    test1.T = 1000.;
    test1.u0 = 300.;
    test1.set_K([](double x) { return COPPER_K; });
    test1.K_type = false; //решаем явным методом
    //test1.set_K([](double x){return 416.31 - 0.05904*x + 7.0872*1e7/(x*x*x);});
    //test1.K_type = true; //Решаем неявным методом
    //test1.set_K([](double x) { return 401*(1+0.003861*(x-293)); });
    //test1.K_type = true; //решаем неявным методом
    test1.set_G_left([&](double x) { return test1.u0; });
    test1.G_left_type = false;
    test1.set_G_right([&](double x) { return test1.u0; });
    test1.G_right_type = false;
    test1.set_init_func([&](double x){ return test1.u0 + x*(test1.L-x); });
    ExplicitScheme(2., 1., 0., test1, "test1");
    test1.show(); // Вывод информации о тесте

    // Тест 2: потоки на концах
    PDE_data test2;
    test2.c = 10.;
    test2.rho = 1000.;
    test2.L = 50.;
    test2.T = 1000.;
    test2.set_K([](double x) { return 1000.; });
    test2.set_G_left([](double x) { return 15.; });
    test2.G_left_type = true;
    test2.set_G_right([](double x) { return 15.; });
    test2.G_right_type = true;
    test2.set_init_func([](double x){ return 1.; });
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
    test5.G_left_type = false;
    test5.set_G_right([](double x) { return 10.; });
    test5.G_right_type = false;
    // начальная температура по всему стержню
    test5.set_init_func([](double x){ return 20.; });
    //test5 end

    /* Случай test 5. */
    //ExplicitScheme(0.01, 0.1, 1., test5);
    //SolvePDE_1(test5, 0.01, 0.1, 0.5, "ExpScheme_test5");




    std::cout << std::endl << "Complete!" << std::endl;
}
