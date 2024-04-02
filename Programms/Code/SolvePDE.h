//
// Объявление функций для решения Уравнений в Частных Производных(УЧП)
//


#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <functional>

// TODO: функция численного интегрирования (реализовать в algebra.h/cpp)

// TODO: функция прогонки (реализовать в algebra.h/cpp)

// TODO: функция для получения данных для теста из файла (реализовать в algebra.h/cpp) *хз нужно ли, удобней пока без нее

// TODO: допилить класс PDE_data, чтобы он работал в двух случаях постановки задачи (реализовать в TESTS.cpp)


/* Функция для решения PDE в случае 1 (по методичке) */
void SolvePDE_1(PDE_data test, double h, double tau, double sigma, std::string filename);