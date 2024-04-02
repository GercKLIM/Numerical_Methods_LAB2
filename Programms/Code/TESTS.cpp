/* ### Тесты ### */



#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <functional>


/* Класс условий задачи */
class PDE_data {

public:
    double c = 0;    // Удельная теплоемкость    (мб в точке)
    double rho = 0;  // Линейная плотность массы (мб в точке)
    double L = 0;    // Длина стержня
    double T = 0;    // Конец временного интервала

    /* Получение функции теплопроводности K(x) */

    // Функция, для присваивания лямбда-функции к функции double K(double)
    void set_K(std::function<double(double)> func) {
        myFunctionK = func;
        K_is_set = true;
    }

    // Функция - коэффициента теплопроводности
    double K(double x) {
        if (myFunctionK) {
            return myFunctionK(x);
        } else {
            return 0;
        }
    }


    /* Получение Граничных условий G_left, G_right*/

    // Функция, для задания функции левой границы
    void set_G_left(std::function<double(double)> func) {
        myFunction_G_left = func;
        G_left_is_set = true;
    }

    // Функция, для задания фунции правой границы
    void set_G_right(std::function<double(double)> func) {
        myFunction_G_right = func;
        G_right_is_set = true;
    }

    // Функция - Левое граничное условие
    double G_left(double x) {
        if (myFunction_G_left) {
            return myFunction_G_left(x);
        } else {
            return 0;
        }
    }

    // Функция - Правое граничное условие
    double G_right(double x) {
        if (myFunction_G_right) {
            return myFunction_G_right(x);
        } else {
            return 0;
        }
    }



    // Вывод информации об объекте
    void show() {
        std::cout << "PDE_data object info:" << std::endl;
        std::cout << "c   = " << c << std::endl;
        std::cout << "rho = " << c << std::endl;
        std::cout << "L   = " << c << std::endl;
        std::cout << "T   = " << c << std::endl;
        std::cout << "K       is " <<  ((K_is_set) ? "set" : "NOT set") << std::endl;
        std::cout << "K       is " <<  ((K_is_set) ? "set" : "NOT set") << std::endl;
        std::cout << "G_left  is " <<  ((G_left_is_set) ? "set" : "NOT set") << std::endl;
        std::cout << "G_right is " <<  ((G_right_is_set) ? "set" : "NOT set") << std::endl;
    }

    std::vector<bool> info(){
        return {K_is_set, G_left_is_set, G_right_is_set};
    }
private:
    // Хранение лямбда-функции как std::function
    std::function<double(double)> myFunctionK;
    std::function<double(double)> myFunction_G_left;
    std::function<double(double)> myFunction_G_right;

    // Состояние заданности условий
    bool K_is_set = false;
    bool G_left_is_set = false;
    bool G_right_is_set = false;
};
