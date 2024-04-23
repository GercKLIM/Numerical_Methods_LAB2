// ### Лабораторная работа №2 ###
// "Численное решение краевых задач для одномерного уравнения теплопроводности"
// Authors:
// @GercKLIM - Климов Олег ФН2-61б
// @Ship-Vano - Шаманов Иван ФН2-61б
//

#include <iostream>
#include "algebra.h"    // Алгебра векторов и матриц
#include "SolvePDE.cpp" // Методы решения PDE
#include "Config.h"     // Константы для тестов



int main() {

    /* Тест 1: Алюминий, фиксированная температура на концах */
    PDE_data test1;
    test1.c = ALUMINUM_C;
    test1.rho = ALUMINUM_RHO;
    test1.h = 0.05;
    test1.L = 1.;
    test1.tau = 0.009;
    //test1.tau = 14;
    test1.T = 100.;
    test1.u0 = 800.;
    test1.set_G_left([&](double x) { return test1.u0; });
    test1.G_left_type = false;
    //test1.set_G_left([&](double t) {return t*t*12;});
    //test1.G_left_type = true;
    test1.set_G_right([&](double x) { return test1.u0; });
    test1.G_right_type = false;
    test1.set_init_func([&](double x){ return test1.u0-500 - x*(test1.L-x); });
    test1.set_K([&](double x, double u){return 237*(1+0.0034*(u-293));});
    //test1.set_K([&](double x, double u){return 13.7+0.0017*u+0.000003*u*u;});
    test1.K_type = true; //Решаем итерационным методом
    if(!test1.K_type) {
        FiniteScheme(test1.tau, test1.h,1.,test1,"test1");
    } else {
        infoIterationScheme(test1.tau, test1.h, 1., test1, "test1_iterational",1e-4);
    }
    test1.set_K([&](double x, double u) { return ALUMINUM_K; });
    test1.K_type = false;
    if(!test1.K_type) {
        FiniteScheme(test1.tau, test1.h,0.5,test1,"test1");
    } else {
        IterationScheme(test1.tau, test1.h, 1., test1, "test1_iterational");
    }
    test1.show(); // Вывод информации о тесте



    /* Тест 2: Алюминий, постоянная температура на левом конце + нулевой поток на правом (теплоизоляция) */
    PDE_data test2;
    test2.c = ALUMINUM_C;
    test2.rho = ALUMINUM_RHO;
    test2.L = 1.;
    test2.T = 100.;
    test2.h = 0.005;
    test2.tau = 0.1;
    test2.u0 = 800.;
    test2.set_G_left([&](double x) { return test2.u0; });
    test2.G_left_type = false;
    test2.set_G_right([&](double x) { return 0.; });
    test2.G_right_type = true;
    test2.set_init_func([&](double x){ return test2.u0-500 - x*(test2.L-x); });
    test2.set_K([&](double x, double u) { return ALUMINUM_K; });
    test2.K_type = false;
    if(!test2.K_type) {
        FiniteScheme(test2.tau, test2.h,1.,test2,"test2");
    } else {
        IterationScheme(test2.tau, test2.h, 1., test2, "test2_iterational");
    }

    /* Тест квазилинейного уравнения (методичка, п.3)*/
    PDE_data test3;
    test3.c = 1.;
    test3.rho = 1.;
    test3.L = 10.;
    test3.T = 1.;
    test3.h = 0.2;
    test3.tau = 0.0002;
    test3.u0 = sqrt(2*5*5/0.5);
    test3.set_G_left([&](double t) { return test3.u0*sqrt(t); });
    test3.G_left_type = false;
    test3.set_G_right([&](double x) { return 0.; });
    test3.G_right_type = true;
    test3.set_init_func([&](double x){ return 0.; });
    test3.set_K([&](double x, double u) { return 0.5*u*u; });
    test3.K_type = true;
    if(!test3.K_type) {
        FiniteScheme(test3.tau, test3.h,1.,test3,"test3");
    } else {
        infoIterationScheme(test3.tau, test3.h, 1., test3, "test3_iterational", 1e-10);
    }

    /* Тест 4: Алюминий, теплоизоляция на концах */
    PDE_data test4;
    test4.c = ALUMINUM_C;
    test4.rho = ALUMINUM_RHO;
    test4.h = 0.005;
    test4.L = 1.;
    test4.tau = 0.0002;
    //test1.tau = 14;
    test4.T = 1.;
    test4.u0 = 800.;
    test4.set_G_left([&](double t) {return 0.;});
    test4.G_left_type = true;
    test4.set_G_right([&](double t) { return 0.; });
    test4.G_right_type = true;
    test4.set_init_func([&](double x){ return test4.u0-500 - x*(test4.L-x); });
    test4.set_K([&](double x, double u) { return ALUMINUM_K; });
    test4.K_type = false;
    if(!test4.K_type) {
        FiniteScheme(test4.tau, test4.h,0.5,test4,"test4-sigma05");
    } else {
        IterationScheme(test4.tau, test4.h, 1., test4, "test4_iterational");
    }
    test4.show(); // Вывод информации о тесте


    /* Тест: Вариант 5 */
    PDE_data test5;
    test5.c = 2.;
    test5.t0 = 0.5;
    test5.rho = 0.25;
    test5.L = 1.;
    test5.T = 10.;
    test3.h = 0.2;
    test3.tau = 0.002;
    test5.set_K([](double x, double u) {
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
    test5.K_type = false;
    test5.set_G_left([](double x) { return 0.2; });
    test5.G_left_type = false;
    test5.set_G_right([](double x) { return 10.; });
    test5.G_right_type = false;
    test5.set_init_func([](double x){ return 20.; });  // Начальная температура по всему стержню
    //FiniteScheme(test5.tau, test5.h, 1., test5, "test5iters");
    /* Случай test 5. */
    //ExplicitScheme(0.01, 0.1, 1., test5);
    //SolvePDE_1(test5, 0.01, 0.1, 0.5, "ExpScheme_test5");



    std::cout << std::endl << "Complete!" << std::endl;
}
