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

    double tau;        // Шаг по времени
    double h;          // Шаг по пространству
    double c = 0;      // Удельная теплоемкость    (мб в точке)
    double rho = 0;    // Линейная плотность массы (мб в точке)
    double x0 = 0;
    double L = 0;      // Длина стержня
    double t0 = 0;     // Начальная температура
    double T = 0;      // Конец временного интервала
    double u0 = 0;     // Начальная температура
    bool G_left_type = false;  // Тип граничных условий слева (0 - первого рода, 1 - второго рода)
    bool G_right_type = false; // Тип граничных условий справа (0 - первого рода, 1 - второго рода)
    bool K_type = false;       // Коэффициент теплопроводности зависит от температуры (0 - K = const, 1 - K = K(T))

    /* Получение функции теплопроводности K(x) */

    // Функция, для присваивания лямбда-функции к функции double K(double)
    void set_K(std::function<double(double, double)> func) {
        myFunctionK = func;
        K_is_set = true;
    }

    // Функция - коэффициента теплопроводности
    double K(double x, double u) {
        if (myFunctionK) {
            return myFunctionK(x, u);
        } else {
            return 0;
        }
    }

    std::function<double(double, double)> K_ptr = [&] (double x, double u) {return K(x, u);};

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

    // Функция для задания начального состояния системы
    void set_init_func(std::function<double(double)> func){
        initFunction = func;
        init_is_set = true;
    }

    // Функция - Левое граничное условие
    double G_left(double t) {
        if (myFunction_G_left) {
            return myFunction_G_left(t);
        } else {
            return 0;
        }
    }

    // Функция - Правое граничное условие
    double G_right(double t) {
        if (myFunction_G_right) {
            return myFunction_G_right(t);
        } else {
            return 0;
        }
    }



    // Вывод информации об объекте
    void show() {
        std::cout << "PDE_data object info:" << std::endl;
        std::cout << "c   = " << c << std::endl;
        std::cout << "rho = " << rho << std::endl;
        std::cout << "L   = " << L << std::endl;
        std::cout << "T   = " << T << std::endl;
        std::cout << "K       is " <<  ((K_is_set) ? "set" : "NOT set") << std::endl;
        std::cout << "G_left  is " <<  ((G_left_is_set) ? "set" : "NOT set") << std::endl;
        std::cout << "G_right is " <<  ((G_right_is_set) ? "set" : "NOT set") << std::endl;
    }

    // Вывод вектора значений класса: заданы ли функции
    std::vector<bool> info(){
        return {G_left_type, G_left_is_set, G_right_type, G_right_is_set, K_is_set};
    }


public:
    std::function<double(double)> initFunction;
private:
    // Хранение лямбда-функции как std::function
    std::function<double(double, double)> myFunctionK;
    std::function<double(double)> myFunction_G_left;
    std::function<double(double)> myFunction_G_right;

    // Состояние заданности условий
    bool K_is_set = false;
    bool G_left_is_set = false;
    bool G_right_is_set = false;
    bool init_is_set = false;
};
