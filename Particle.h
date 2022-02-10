//
//  Particle.h
//  Calculation of partition functions and  heats capacity
//
//  Created by Дмитрий on 26.12.2021.
//

#ifndef Particle_h
#define Particle_h

#include <cmath>
#include <vector>
#include <string>
using namespace std;

// Класс "Частица"
class Particle{
    
protected:
    
    //Для перевода табличных значений к нужному порядку
    double cm_to_Joule = 1.986e-23;
    double cm_to_m = 100;
    
    double k = 1.3806e-23; // Постоянная Больцмана, Дж / K
    double h = 6.6261e-34; // Постоянная Планка, Дж * c
        
public:
    virtual vector<double> Z_tr(int temperature, string type){
        vector<double> Z_tr_array(temperature);
        return Z_tr_array;
    }// Поступательная статистическая сумма
    virtual vector<double> Z_int(int temperature, string type){
        vector<double> Z_int_array(temperature);
        return Z_int_array;
    }// Внутренняя статистическая сумма
    virtual vector<double> Z(int temperature, string type){
        vector<double> Z_array(temperature);
        return Z_array;
    } //  Полная статистическая сумма
    
    virtual vector<double> E_tr(int temperature, string type){
        vector<double> E_tr_array(temperature);
        return E_tr_array;
    } // Поступательная удельная энергия
    virtual vector<double> E_int(int temperature, string type){
        vector<double> E_int_array(temperature);
        return E_int_array;
    } //  Внутренняя удельная энергия
    virtual vector<double> E(int temperature, string type){
        vector<double> E_array(temperature);
        return E_array;
    } // Полная удельная энергия
    
    virtual vector<double> c_V_tr(int temperature, string type){
        vector<double> c_V_tr_array(temperature);
        return c_V_tr_array;
    } // Поступательная удельная теплоёмкость при постоянном объёме
    virtual vector<double> c_V_int(int temperature, string type){
        vector<double> c_V_int_array(temperature);
        return c_V_int_array;
    } // Внутренняя удельная теплоёмкость при постоянном объёме
    virtual vector<double> c_V(int temperature, string type){
        vector<double> c_V_array(temperature);
        return c_V_array;
    } // Полная удельная теплоёмкость при постоянном объёме
};
#endif /* Particle_h */
