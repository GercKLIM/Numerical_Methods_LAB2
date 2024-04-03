// ### Лабораторная работа №2 ###
// "Численное решение краевых задач для одномерного уравнения теплопроводности"
// Authors:
// @GercKLIM - Климов Олег ФН2-61б
// @Ship-Vano - Шаманов Иван ФН2-61б
//

#include <iostream>
#include "algebra.h"  // Алгебра векторов и матриц
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
    test1.show(); // Вывод информации о тесте
    //std::cout << test1.info();
    //std::cout << "K(2) = " << test1.K(2) << std::endl;

    std::cout << std::endl << "Complete!" << std::endl;

    ExplicitScheme(2, 1, 1, test1);
}

