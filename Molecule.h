//
//  Molecule.h
//  Calculation of partition functions and  heats capacity
//
//  Created by Дмитрий on 26.12.2021.
//

#ifndef Molecule_h
#define Molecule_h

#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include "Particle.h"
#include <string>
using namespace std;


// Класс "Молекула", наследуемый от "Частица"
class Molecule : public Particle{
    
public:
    
    //Поля класса
    string name; // Название молекулы
    double a_e[11];
    double e_el[11]; // Энергия электрона молекулы, Дж
    double g_el[11]; // Статистический вес
    
    // Спектроскопические константы, Дж
    double w_e[11];
    double wx_e[11];
    double wy_e[11];
    double B_e[11];
    double D_e[11];
    double b_e[11];
    
    double e_diss[11]; // Энергия диссоциации, Дж
    double m; //  Молекулярная масса, кг
    
    //Сеттеры:
    void Set_a_e(double a_e_input[]){
        for(int i = 0; i < 11; i++){
            a_e[i] = a_e_input[i] * cm_to_Joule / 100;
        }
    }
    void Set_e_el(double e_el_input[]){
        for(int i = 0; i < 11; i++){
            e_el[i] = e_el_input[i] * cm_to_Joule;
        }
    }
    void Set_g_el(double g_el_input[]){
        for(int i = 0; i < 11; i++){
            g_el[i] = g_el_input[i];
        }
    }
    void Set_w_e(double w_e_input[]){
        for(int i = 0; i < 11; i++){
            w_e[i] = w_e_input[i] * cm_to_Joule;
        }
    }
    void Set_wx_e(double wx_e_input[]){
        for(int i = 0; i < 11; i++){
            wx_e[i] = wx_e_input[i] * cm_to_Joule;
        }
    }
    void Set_wy_e(double wy_e_input[]){
        for(int i = 0; i < 11; i++){
            wy_e[i] = wy_e_input[i] * cm_to_Joule;
        }
    }
    void Set_B_e(double B_e_input[]){
        for(int i = 0; i < 11; i++){
            B_e[i] = B_e_input[i] * cm_to_Joule;
        }
    }
    void Set_D_e(double D_e_input[]){
        for(int i = 0; i < 11; i++){
            D_e[i] = D_e_input[i] * cm_to_Joule / 1000000;
        }
    }
    void Set_b_e(double b_e_input[]){
        for(int i = 0; i < 11; i++){
            b_e[i] = b_e_input[i] * cm_to_Joule;
        }
    }
    void Set_e_diss(double e_diss_input[]){
        for(int i = 0; i < 11; i++){
            e_diss[i] = abs(e_diss_input[i]) * cm_to_Joule + e_el[i];
            //e_diss[i] = e_diss_input[i] * cm_to_Joule;
        }
    }
    void Set_m(double m_input){
        
        m = m_input;
    }
    
    //Методы класса:
    
    // Колебательная энергия
    double e_vibr(int n, int i){
        return w_e[n] * (i + 0.5) - wx_e[n] * pow(i + 0.5, 2) + wy_e[n] * pow(i + 0.5, 3);
    }
    
    // Вращательная энергия
    double e_rot(int n, int i, int j){
        double B_i = B_e[n] - a_e[n] * (i + 0.5);
        double D_i = D_e[n] - b_e[n] * (i + 0.5);
        return B_i * j * (j + 1) - D_i * pow(j, 2) * pow(j + 1, 2);
    }
    
    // Внутренняя энергия
    double e_int(int n, int i, int j){
        return e_el[n] + e_vibr(n, i) + e_rot(n, i, j);
    }
    
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
            int i = 0;
            int j = 0;
            for(int t = 0; t < temperature; t++){
                for(int n = 0; n < 5; n ++){
                    while(e_int(n, i, j) < e_diss[n]){
                        while(e_int(n, i, j) < e_diss[n]){
                            Z_int_array[t] = Z_int_array[t] + g_el[n] * 1 * (2 * j + 1) * exp((- e_int(n, i, j)) / (k * t));
                            j = j + 1;
                        }
                        i = i + 1;
                        j = 0;
                    }
                    i = 0;
            }
        }
        }
        else if(type == "Single"){
            Z_int_array[0] = 0;
            int i = 0;
            int j = 0;
            for(int n = 0; n < 5; n ++){
                while(e_int(n, i, j) < e_diss[n]){
                    while(e_int(n, i, j) < e_diss[n]){
                        Z_int_array[0] = Z_int_array[0] + g_el[n] * 1 * (2 * j + 1) * exp((- e_int(n, i, j)) / (k * temperature));
                        j = j + 1;
                    }
                    i = i + 1;
                    j = 0;
                }
                i = 0;
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
            int i = 0;
            int j = 0;
            double i_sum = 0;
            double j_sum = 0;
            double n_sum = 0;
            
            for(int t = 0; t < temperature; t++){
                for(int n = 0; n < 5; n ++){
                    double e_in = e_int(n, i, j);
                    while(e_in < e_diss[n]){
                        while(e_in < e_diss[n]){
                            j = j + 1;
                            j_sum = j_sum + e_rot(n, i, j) * (2 * j + 1) * exp((- e_rot(n, i, j)) / (k * t));
                            e_in = e_int(n, i, j);
                        }
                        i_sum = i_sum + e_vibr(n, i) * 1 * exp((- e_vibr(n, i)) / (k * t)) * j_sum;
                        j_sum = 0;
                        i = i + 1;
                        j = 0;
                        e_in = e_int(n, i, j);
                    }
                    n_sum = n_sum + e_el[n] * g_el[n] * exp((- e_el[n]) / (k * t)) * i_sum;
                    i_sum = 0;
                    i = 0;
                    j = 0;
                }
                E_int_array[t] = n_sum / (m * Z_int(t, "Single")[0]);
                n_sum = 0;
            }
    }
        else if(type == "Single"){
            int i = 0;
            int j = 0;
            double i_sum = 0;
            double j_sum = 0;
            double n_sum = 0;
            for(int n = 0; n < 5; n ++){
                double e_in = e_int(n, i, j);
                while(e_in < e_diss[n]){
                    while(e_in < e_diss[n]){
                        j = j + 1;
                        j_sum = j_sum + e_rot(n, i, j) * (2 * j + 1) * exp((- e_rot(n, i, j)) / (k * temperature));
                        e_in = e_int(n, i, j);
                    }
                    i_sum = i_sum + e_vibr(n, i) * 1 * exp((- e_vibr(n, i)) / (k * temperature)) * j_sum;
                    j_sum = 0;
                    i = i + 1;
                    j = 0;
                    e_in = e_int(n, i, j);
                }
                n_sum = n_sum + e_el[n] * g_el[n] * exp((- e_el[n]) / (k * temperature)) * i_sum;
                i_sum = 0;
                i = 0;
                j = 0;
            }
            E_int_array[0] = n_sum / (m * Z_int(temperature, "Single")[0]);
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
        if(type == "Multiply"){
            for(int r = 0; r <  temperature; r++){
                c_V_int_array[r] = 0;
            }
            int i = 0;
            int j = 0;
            double i_sum = 0;
            double j_sum = 0;
            double n_sum = 0;
            
            for(int t = 0; t < temperature; t++){
                vector<double> sum1(temp);
                vector<double> sum2(temp);
                for(int n = 0; n < 5; n ++){
                    double e_in = e_int(n, i, j);
                    while(e_in < e_diss[n]){
                        while(e_in < e_diss[n]){
                            j = j + 1;
                            j_sum = j_sum + pow(e_rot(n, i, j), 2) * (2 * j + 1) / (pow(k, 2) * pow(t, 2)) * exp((- e_rot(n, i, j)) / (k * t));
                            e_in = e_int(n, i, j);
                        }
                        i_sum = i_sum + pow(e_vibr(n, i), 2) * 1 / (pow(k, 2) * pow(t, 2))  * exp((- e_vibr(n, i)) / (k * t)) * j_sum;
                        j_sum = 0;
                        i = i + 1;
                        j = 0;
                        e_in = e_int(n, i, j);
                    }
                    n_sum = n_sum + pow(e_el[n], 2) * g_el[n] / (pow(k, 2) * pow(t, 2)) * exp((- e_el[n]) / (k * t)) * i_sum;
                    i_sum = 0;
                    i = 0;
                    j = 0;
                }
                sum1[t] = n_sum / Z_int(t, "Single")[0];
                n_sum = 0;
                for(int n = 0; n < 5; n ++){
                    double e_in = e_int(n, i, j);
                    e_in = e_int(n, i, j);
                    while(e_in < e_diss[n]){
                        while(e_in < e_diss[n]){
                            j = j + 1;
                            j_sum = j_sum + e_rot(n, i, j) * (2 * j + 1) / (k * t) * exp(- e_rot(n, i, j) / (k * t));
                            e_in = e_int(n, i, j);
                        }
                        i_sum = i_sum + e_vibr(n, i) * 1 / (k * t) * exp(- e_vibr(n, i) / (k * t)) * j_sum;
                        j_sum = 0;
                        i = i + 1;
                        j = 0;
                        e_in = e_int(n, i, j);
                    }
                    n_sum = n_sum + e_el[n] * g_el[n] / (k * t) * exp(- e_el[n] / (k * t)) * i_sum;
                    i_sum = 0;
                    i = 0;
                    j = 0;
                }
                sum2[t] = pow((n_sum / Z_int(t, "Single")[0]), 2);
                c_V_int_array[t] = k / m * (sum1[t] - sum2[t]);
            }
    }
        else if(type == "Single"){
            int i = 0;
            int j = 0;
            double i_sum = 0;
            double j_sum = 0;
            double n_sum = 0;
                vector<double> sum1(1);
                vector<double> sum2(1);
                for(int n = 0; n < 5; n ++){
                    double e_in = e_int(n, i, j);
                    while(e_in < e_diss[n]){
                        while(e_in < e_diss[n]){
                            j = j + 1;
                            j_sum = j_sum + pow(e_rot(n, i, j), 2) * (2 * j + 1) / (pow(k, 2) * pow(temperature, 2)) * exp((- e_rot(n, i, j)) / (k * temperature));
                            e_in = e_int(n, i, j);
                        }
                        i_sum = i_sum + pow(e_vibr(n, i), 2) * 1 / (pow(k, 2) * pow(temperature, 2))  * exp((- e_vibr(n, i)) / (k * temperature)) * j_sum;
                        j_sum = 0;
                        i = i + 1;
                        j = 0;
                        e_in = e_int(n, i, j);
                    }
                    n_sum = n_sum + pow(e_el[n], 2) * g_el[n] / (pow(k, 2) * pow(temperature, 2)) * exp((- e_el[n]) / (k * temperature)) * i_sum;
                    i_sum = 0;
                    i = 0;
                    j = 0;
                }
                sum1[0] = n_sum / Z_int(temperature, "Single")[0];
                n_sum = 0;
                for(int n = 0; n < 5; n ++){
                    double e_in = e_int(n, i, j);
                    e_in = e_int(n, i, j);
                    while(e_in < e_diss[n]){
                        while(e_in < e_diss[n]){
                            j = j + 1;
                            j_sum = j_sum + e_rot(n, i, j) * (2 * j + 1) / (k * temperature) * exp(- e_rot(n, i, j) / (k * temperature));
                            e_in = e_int(n, i, j);
                        }
                        i_sum = i_sum + e_vibr(n, i) * 1 / (k * temperature) * exp(- e_vibr(n, i) / (k * temperature)) * j_sum;
                        j_sum = 0;
                        i = i + 1;
                        j = 0;
                        e_in = e_int(n, i, j);
                    }
                    n_sum = n_sum + e_el[n] * g_el[n] / (k * temperature) * exp(- e_el[n] / (k * temperature)) * i_sum;
                    i_sum = 0;
                    i = 0;
                    j = 0;
                }
                sum2[0] = pow((n_sum / Z_int(temperature, "Single")[0]), 2);
                c_V_int_array[0] = k / m * (sum1[0] - sum2[0]);
            
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

#endif /* Molecule_h */

