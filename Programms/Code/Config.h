//
// Created by Иван on 4/14/2024.
//

#ifndef CODE_CONFIG_H
#define CODE_CONFIG_H

// ВСЕ ЗНАЧЕНИЯ ДАНЫ В СИ
// [rho] = [кг/м^3] плотность
// [c] = [Дж/(кг*К)] теплоёмкость
// [K] = [Вт/(м*К)] теплопроводность

//COPPER: (ГОСТ 859-78)
double COPPER_RHO = 8500;
double COPPER_C = 4200;
double COPPER_K = 407;

//ALUMINUM:  (ГОСТ 22233-83)
double ALUMINUM_RHO = 2600;
double ALUMINUM_C = 840;
double ALUMINUM_K = 221;

//STEEL: (Сталь стержневая арматурная (ГОСТ 10884-81))
double STEEL_RHO = 7850;
double STEEL_C = 482;
double STEEL_K = 58;

//
#endif //CODE_CONFIG_H
