// ### Лабораторная работа №2 ###
// "Численное решение краевых задач для одномерного уравнения теплопроводности"
// Authors:
// @GercKLIM - Климов Олег ФН2-61б
// @Ship-Vano - Шаманов Иван ФН2-61б
//

#include <iostream>
#include "algebra.cpp"

int main() {

    // Проверка работы алгебры
    std::vector<std::vector<double>> a = create_identity_matrix<double>(2); // Создание единичной матрицы
    std::vector<std::vector<double>> b = {{2., 2.}, {2., 2.}};
    std::cout << a * b - a << std::endl; // Операции над матрицами и вывод в консоль
    return 0;
}

