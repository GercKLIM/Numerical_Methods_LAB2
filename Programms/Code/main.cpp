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


void test0() {

    /* Тест 0: C анал решением */
    PDE_data test0;
    test0.c = 1;
    test0.rho = 1;
    test0.L = 1.;
    test0.T = 1.;
    //test0.u0 = 800.;

    test0.set_G_left([&](double x) { return 0; });
    test0.G_left_type = false;
    test0.set_G_right([&](double x) { return 0; });
    test0.G_right_type = false;
    test0.set_init_func([&](double x) { return sin(3.145657 * x); });
    test0.set_K([&](double x, double u) { return 1; });
    test0.K_type = true; //Решаем итерационным методом

    FiniteScheme(test0.tau, test0.h, 1., test0, "test1");

    //test1.show(); // Вывод информации о тесте
}


void test1() {

    /* Тест 1: Алюминий, фиксированная температура на концах */
    PDE_data test1;
    test1.c = ALUMINUM_C;
    test1.rho = ALUMINUM_RHO;
    test1.h = 0.005;
    test1.L = 1.;
    test1.tau = 0.009;
    test1.tau = 1.;
    test1.T = 100.;
    test1.u0 = 800.;
    test1.set_G_left([&](double x) { return test1.u0; });
    test1.G_left_type = false;
    test1.set_G_right([&](double x) { return test1.u0; });
    test1.G_right_type = false;
    test1.set_init_func([&](double x) { return test1.u0 - 500 - x * (test1.L - x); });
    test1.set_K([&](double x, double u) { return 237 * (1 + 0.0034 * (u - 293)); });
    test1.K_type = true; //Решаем итерационным методом
    if (!test1.K_type) {
        FiniteScheme(test1.tau, test1.h, 1., test1, "test1");
    } else {
        IterationScheme(test1.tau, test1.h, 0., test1, "test1_iterational");
    }
    test1.set_K([&](double x, double u) { return ALUMINUM_K; });
    test1.K_type = false;
    if (!test1.K_type) {
        FiniteScheme(test1.tau, test1.h, 1., test1, "test1");
    } else {
        IterationScheme(test1.tau, test1.h, 0., test1, "test1_iterational");
    }
    //test1.show(); // Вывод информации о тесте
}

void test2() {

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
    test2.set_G_right([&](double x) { return 0; });
    test2.G_right_type = true;
    test2.set_init_func([&](double x) { return test2.u0 - 500 - x * (test2.L - x); });
    test2.set_K([&](double x, double u) { return ALUMINUM_K; });
    test2.K_type = false;
    if (!test2.K_type) {
        FiniteScheme(test2.tau, test2.h, 1., test2, "test2");
    } else {
        IterationScheme(test2.tau, test2.h, 0., test2, "test2_iterational");
    }
}

void test3() {
    /* Тест: Вариант 5 */
    PDE_data test5;
    test5.c = 2;
    test5.rho = 0.25;
    test5.L = 1;
    test5.t0 = 0.5;
    test5.T = 10;
    test5.u0 = 0.2;
    test5.set_K([](double x, double u) {
        double x1 = 0.5, x2 = 2. / 3.;
        double k1 = 2., k2 = 0.5;
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
    test5.set_init_func([](double x) { return 20.; });  // Начальная температура по всему стержню

    /* Случай test 5. */
    //ExplicitScheme(0.01, 0.1, 1., test5);
    //SolvePDE_1(test5, 0.01, 0.1, 0.5, "ExpScheme_test5");
//    if (!test5.K_type) {
//        FiniteScheme(test5.tau, test5.h, 1., test5, "test2");
//    } else {
//        IterationScheme(test5.tau, test5.h, 0., test5, "test2_iterational");
//    }
}

bool MyTest(std::string filepath){

    // Открытие файлы для чтения
    std::ifstream file(filepath);

    if (file.is_open()) {
        file.close();

        vector<double> param= ImportData(filepath);

        PDE_data test;
        test.c = param[0];
        test.rho = param[1];
        test.h = param[2];
        test.L = param[3];
        test.tau = param[4];
        test.T = param[5];
        test.u0 = param[6];
        test.G_left_type = param[7];
        test.G_right_type = param[8];
        test.K_type = param[9]; // Выбор метода

        test.set_G_left([&](double x) { return test.u0; });
        test.set_G_right([&](double x) { return test.u0; });

        int num_K, num_init_func;
        std::cout << "Choice K:";
        std::cin >> num_K;
        std::cout << "Choice init function:";
        std::cin >> num_init_func;


        switch (num_K) {
            case 1:
                test.set_init_func([&](double x) { return test.u0 - 500 - x * (test.L - x); });
                break;
            case 2:
                test.set_init_func([&](double x) { return 1; });
                break;
            case 3:
                test.set_init_func([&](double x) { return x; });
                break;
            default:
                test.set_init_func([&](double x) { return 1; });
                break;
        }

        switch (num_K) {
            case 1:
                test.set_K([&](double x, double u) { return 237 * (1 + 0.0034 * (u - 293)); });
                break;
            case 2:
                test.set_K([&](double x, double u) { return 1; });
                break;
            case 3:
                test.set_K([&](double x, double u) { return x; });
                break;
            default:
                test.set_K([&](double x, double u) { return 1; });
                break;
        }


        // TODO: Написать вызов метода и выбор функций

        if (!test.K_type) {
            FiniteScheme(test.tau, test.h, 1., test, "testNULL");
        } else {
            IterationScheme(test.tau, test.h, 0., test, "testNULL_iterational");
        }




        return true;
    } else {
        return false;
    }
}



/* Функция, которая делает наш код программой с потоком данных */
void programm(){

    // Приветствие
    std::cout << std::endl;
    std::cout << "WELKOME TO MY RPOGRAMM! :) "<< std::endl;
    std::cout << std::endl;

    // Инструкция
    std::cout << std::endl;

    // Бесконечный цикл
    while (true){
        int num_test;
        std::cout << "Enter test number:  ";
        std::cin >> num_test;
        std::cout << std::endl;


        if (num_test == 1) {
            test1();
            std::cout << "Complete! What's next, Boss?" << std::endl << std::endl;

        } else if (num_test == 2){
            test2();
            std::cout << "Complete! What's next, Boss?" << std::endl << std::endl;

        } else if (num_test == 3){
            test3();
            std::cout << "Complete! What's next, Boss?" << std::endl << std::endl;

        } else if (num_test == 4){

            // Снова цикл, пока не введем "0"
            while (true) {
                std::cout << "Enter filepath your parametres: ";
                string test_filepath;
                std::cin >> test_filepath;


                if (test_filepath == "0") {
                    break;
                }

                if (!MyTest(test_filepath)) {
                    std::cout << "Error: file not found, please try again :( " << std::endl << std::endl;
                } else {
                    std::cout << "Complete! What's next, Boss?" << std::endl << std::endl;
                }

            }

        } else if (num_test == 0){
            break;
        }

    }

    std::cout << "Thanks for used! See u again ;)" << std::endl << std::endl;

}



/* Функция для подготовки данных для создания таблиц */
void make_data_for_tables() {

    /* Тест 0: наипростейшая задача теплопроводности с существующим аналитическим решением */
    PDE_data test0;
    test0.c = 1;
    test0.rho = 1;
    test0.L = 0.8;
    test0.T = 0.8;

    // Граничные условия
    test0.set_G_left([&](double x) { return 0; });
    test0.G_left_type = false;
    test0.set_G_right([&](double x) { return 0; });
    test0.G_right_type = false;

    // Начальное условие
    test0.set_init_func([&](double x) { return sin(M_PI * x); });

    // K = 1
    test0.set_K([&](double x, double u) { return 1; });


    // При sigma = 0.0 -> O(h^2 + tau)
    //     sigma = 0.5 -> O(h^2 + tau^2)
    //     sigma = 1.0 -> O(h^2 + tau)

    std::vector<double> sigmas = {0.0, 0.5, 1.0};



    for (int n = 0; n < sigmas.size(); n++) {

        double tau = 0.02;
        double h = 0.2;

        for (int i = 0; i < 6; i++) {

            double sigma = sigmas[n];
            test0.tau = tau;
            test0.h = h;
            FiniteScheme(test0.tau, test0.h, sigma, test0, "data_for_tables/test0/sigma_" + to_string(n) + "/test0_" + to_string(i) + ".txt");

            // Обновляем шаги
            if (n == 1) {
                tau = tau / 4.;
                h = h / 4.;
            } else {
                tau = tau / 4.;
                h = h / 2.;
            }

        }

    }

//    /* Запись аналитического решения test0 */
//
//    auto anal_sol = [&](double x, double t) -> double {
//        return std::exp(-M_PI * M_PI * t) * std::sin(M_PI * x);
//    };



}


int main(){
    programm();
    //make_data_for_tables();

    return 0;
}