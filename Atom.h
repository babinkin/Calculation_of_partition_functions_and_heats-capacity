//
//  Atom.h
//  Calculation of partition functions and  heats capacity
//
//  Created by Дмитрий on 26.12.2021.
//

#ifndef Atom_h
#define Atom_h

#include <cmath>
#include "Particle.h"

// Класс "Атом", наследуемый от "Частица"
class Atom : public Particle{
    
public:
    
    //Поля класса
    string name; // Название атома
    double e_el[431]; // Энергия электрона молекулы, Дж
    double g_el[431]; // Статистический вес
    double m; //  Молекулярная масса, кг
    
    //Сеттеры:
    void Set_e_el(double e_el_input[]){
        for(int i = 0; i < 431; i++){
            e_el[i] = e_el_input[i] * 1.986e-23; 
        }
    }
    void Set_g_el(double g_el_input[]){
        for(int i = 0; i < 431; i++){
            g_el[i] = g_el_input[i];
        }
    }
    void Set_m(double m_input){
        
        m = m_input;
    }
    
    //Методы класса:
    
    // Поступательная стат. сумма
    vector<double> Z_tr(int temperature, string type) override{
        int temp;
        if(type == "Multiply"){
            temp = temperature;
        }
        else if(type == "Single"){
            temp = 1;
        }
        else temp = 0;
        vector<double> Z_tr_array(temp);
        if(type == "Multiply"){
            for(int i = 0; i < temperature; i++){
                Z_tr_array[i] = pow(2 * 3.14 * m * k * i / pow(h, 2), 3 / 2);
            }
        }
        else if(type == "Single"){
            
            Z_tr_array[0] = pow(2 * 3.14 * m * k * temperature / pow(h, 2), 3 / 2);
            return Z_tr_array;
            
        }
        return Z_tr_array;
    }
    // Внутренняя стат. сумма
    vector<double> Z_int(int temperature, string type) override{
        int temp;
        if(type == "Multiply"){
            temp = temperature;
        }
        else if(type == "Single"){
            temp = 1;
        }
        else temp = 0;
        vector<double> Z_int_array(temp);
        if(type == "Multiply"){
            for(int r = 0; r <  temperature; r++){
                Z_int_array[r] = 0;
            }
            for(int t = 0; t < temperature; t++){
                for(int n = 0; n < 20; n ++){
                    Z_int_array[t] = Z_int_array[t] + g_el[n]  * exp((- e_el[n]) / (k * t));
                           
            }
        }
        }
        else if(type == "Single"){
            Z_int_array[0] = 0;
            for(int n = 0; n < 20; n ++){
                        Z_int_array[0] = Z_int_array[0] + g_el[n] * exp((- e_el[n]) / (k * temperature));
            }
            return Z_int_array;
        }
        return Z_int_array;
    
        
    }
    // Полная стат. сумма
    vector<double> Z(int temperature, string type) override{
        int temp;
        if(type == "Multiply"){
            temp = temperature;
        }
        else if(type == "Single"){
            temp = 1;
        }
        else temp = 0;
        vector<double> Z_array(temp);
        if(type == "Multiply"){
            for(int t = 0; t < temperature; t++){
                vector<double> a = Z_tr(t, "Single");
                vector<double> b = Z_int(t, "Single");
                Z_array[t] = a[0] * b[0];
            }
    }
        else if(type == "Single"){
            vector<double> a = Z_tr(temperature, "Single");
            vector<double> b = Z_int(temperature, "Single");
            Z_array[0] = a[0] * b[0];
            return Z_array;
        }
        return Z_array;
    }
    
    //Поступательная удельная энергия
    vector<double> E_tr(int temperature, string type) override{
        int temp;
        if(type == "Multiply"){
            temp = temperature;
        }
        else if(type == "Single"){
            temp = 1;
        }
        else temp = 0;
        vector<double> E_tr_array(temp);
        if(type == "Multiply"){
            for(int r = 0; r <  temperature; r++){
                E_tr_array[r] = 0;
            }
            for(int t = 0; t < temperature; t++){
                E_tr_array[t] = 3 * k * t / (2 * m);
            }
            
    }
        else if(type == "Single"){
            E_tr_array[0] = 3 * k * temperature / (2 * m);
            return E_tr_array;
        }
        return E_tr_array;
    }
    
    //Внутренняя удельная энергия
    vector<double> E_int(int temperature, string type) override{
        int temp;
        if(type == "Multiply"){
            temp = temperature;
        }
        else if(type == "Single"){
            temp = 1;
        }
        else temp = 0;
        vector<double> E_int_array(temp);
        if(type == "Multiply"){
            for(int r = 0; r <  temperature; r++){
                E_int_array[r] = 0;
            }
            for(int t = 0; t < temperature; t++){
                for(int n = 0; n < 20; n ++){
                    E_int_array[t] = E_int_array[t] + e_el[n] + g_el[n] * exp((- e_el[n]) / (k * t)) / (m * Z_int(t, "Single")[0]);
                }
            }
            
    }
        else if(type == "Single"){
            E_int_array[0] = 0;
            for(int n = 0; n < 20; n ++){
            E_int_array[0] = E_int_array[0] + e_el[n] + g_el[n] * exp((- e_el[n]) / (k * temperature)) / (m * Z_int(temperature, "Single")[0]);
            }
            return E_int_array;
        }
        return E_int_array;
    }
    
    // Полная удельная энергия
    vector<double> E(int temperature, string type) override{
        int temp;
        if(type == "Multiply"){
            temp = temperature;
        }
        else if(type == "Single"){
            temp = 1;
        }
        else temp = 0;
        vector<double> E_array(temp);
        if(type == "Multiply"){
            for(int r = 0; r <  temperature; r++){
                E_array[r] = 0;
            }
            for(int t = 0; t < temperature; t++){
                E_array[t] = E_tr(t, "Single")[0] * E_int(t, "Single")[0];
            }
            
    }
        else if(type == "Single"){
            E_array[0] = E_tr(temperature, "Single")[0] * E_int(temperature, "Single")[0];
            return E_array;
        }
        return E_array;
    }
    
    // Поступательная удельная теплоёмкость при постоянном объёме
    vector<double> c_V_tr(int temperature, string type) override{
        int temp;
        if(type == "Multiply"){
            temp = temperature;
        }
        else if(type == "Single"){
            temp = 1;
        }
        else temp = 0;
        vector<double> c_V_tr_array(temp);
        if(type == "Multiply"){
            for(int r = 0; r <  temperature; r++){
                c_V_tr_array[r] = 0;
            }
            for(int t = 0; t < temperature; t++){
                c_V_tr_array[t] = 3 * k  / (2 * m);
            }
            
    }
        else if(type == "Single"){
            c_V_tr_array[0] = 3 * k  / (2 * m);
            return c_V_tr_array;
        }
        return c_V_tr_array;
    }
    // Внутренняя удельная теплоёмкость при постоянном объёме
    vector<double> c_V_int(int temperature, string type) override{
        int temp;
        if(type == "Multiply"){
            temp = temperature;
        }
        else if(type == "Single"){
            temp = 1;
        }
        else temp = 0;
        vector<double> c_V_int_array(temp);
        vector<double> sum1(temp);
        vector<double> sum2(temp);
        if(type == "Multiply"){
            for(int r = 0; r <  temperature; r++){
                c_V_int_array[r] = 0;
                sum1[r] = 0;
                sum2[r] = 0;
            }
            for(int t = 0; t < temperature; t++){
                for(int n = 0; n < 20; n ++){
                    sum1[t] = pow(e_el[n], 2) * g_el[n] * exp((- e_el[n]) / (k * t)) / (pow(k, 2) * pow(t, 2));
                    sum2[t] = e_el[n] * g_el[n] * exp((- e_el[n]) / (k * t)) / (k * t);
                }
                c_V_int_array[t] = k / m * (sum1[t] / Z_int(t, "Single")[0] - pow(sum2[t] / Z_int(t, "Single")[0], 2));
            }
            
    }
        else if(type == "Single"){
            vector<double> sum1(1);
            vector<double> sum2(1);
            for(int n = 0; n < 20; n ++){
                sum1[0] = pow(e_el[n], 2) * g_el[n] * exp((- e_el[n]) / (k * temperature)) / (pow(k, 2) * pow(temperature, 2));
                sum2[0] = e_el[n] * g_el[n] * exp((- e_el[n]) / (k * temperature)) / (k * temperature);
            }
            c_V_int_array[0] = k / m * (sum1[0] / Z_int(temperature, "Single")[0] - pow(sum2[0] / Z_int(temperature, "Single")[0], 2));
            return c_V_int_array;
        }
        return c_V_int_array;
    }
    // Полная удельная теплоёмкость при постоянном объёме
    vector<double> c_V(int temperature, string type) override{
        int temp;
        if(type == "Multiply"){
            temp = temperature;
        }
        else if(type == "Single"){
            temp = 1;
        }
        else temp = 0;
        vector<double> c_V_array(temp);
        if(type == "Multiply"){
            for(int r = 0; r <  temperature; r++){
                c_V_array[r] = 0;
            }
            for(int t = 0; t < temperature; t++){
                c_V_array[t] = c_V_tr(t, "Single")[0] + c_V_int(t, "Single")[0];
            }
            
    }
        else if(type == "Single"){
            c_V_array[0] = c_V_tr(temperature, "Single")[0] + c_V_int(temperature, "Single")[0];
            return c_V_array;
        }
        return c_V_array;
    }
};

#endif /* Atom_h */
