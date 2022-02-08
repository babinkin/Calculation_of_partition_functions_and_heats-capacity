//
//  main.cpp
//  Calculation of partition functions and  heats capacity
//
//  Created by Дмитрий on 26.12.2021.
//

#include <iostream>
#include <cmath>
#include "Particle.h"
#include "Molecule.h"
#include "Atom.h"
#include <fstream>
#include <vector>

using namespace std;

int main() {
    
    setlocale(LC_ALL, "ru");
    
    //Экземпляр класса молекулы N2 (Динитроген)
    Molecule Dinitrogen;
    
    // Заполняем массив значений спектроскопических постоянных данными из файлов в папке data for N2
    
    //alpha_e
   double a_e_input[11];
    ifstream a_e_data("data for N2/a_e.txt");
    if(a_e_data)
    {
        int n;
        double t;
        a_e_data >> n;
        for (int i = 0; i < n; i++) {
            a_e_data >> t;
            a_e_input[i] = t;
            Dinitrogen.Set_a_e(a_e_input);
        }
        a_e_data.close();
    }else
        cout<<"Файл a_e.txt не открылся" << endl;
    
    //e_el
    double e_el_input[11];
     ifstream e_el_data("data for N2/e_el.txt");
     if(e_el_data)
     {
         int n;
         double t;
         e_el_data >> n;
         for (int i = 0; i < n; i++) {
             e_el_data >> t;
             e_el_input[i] = t;
             Dinitrogen.Set_e_el(e_el_input);
         }
         e_el_data.close();
     }else
         cout<<"Файл e_el.txt не открылся" << endl;
     
    //g_el
    double g_el_input[11];
     ifstream g_el_data("data for N2/g_el.txt");
     if(g_el_data)
     {
         int n;
         double t;
         g_el_data >> n;
         for (int i = 0; i < n; i++) {
             g_el_data >> t;
             g_el_input[i] = t;
             Dinitrogen.Set_g_el(g_el_input);
         }
         g_el_data.close();
     }else
         cout<<"Файл g_el.txt не открылся" << endl;
    
    //w_e
    double w_e_input[11];
     ifstream w_e_data("data for N2/w_e.txt");
     if(w_e_data)
     {
         int n;
         double t;
         w_e_data >> n;
         for (int i = 0; i < n; i++) {
             w_e_data >> t;
             w_e_input[i] = t;
             Dinitrogen.Set_w_e(w_e_input);
         }
         w_e_data.close();
     }else
         cout<<"Файл w_e.txt не открылся" << endl;
     
    //wx_e
    double wx_e_input[11];
     ifstream wx_e_data("data for N2/wx_e.txt");
     if(wx_e_data)
     {
         int n;
         double t;
         wx_e_data >> n;
         for (int i = 0; i < n; i++) {
             wx_e_data >> t;
             wx_e_input[i] = t;
             Dinitrogen.Set_wx_e(wx_e_input);
         }
         wx_e_data.close();
     }else
         cout<<"Файл wx_e.txt не открылся" << endl;
     
    //wy_e
    double wy_e_input[11];
     ifstream wy_e_data("data for N2/wy_e.txt");
     if(wy_e_data)
     {
         int n;
         double t;
         wy_e_data >> n;
         for (int i = 0; i < n; i++) {
             wy_e_data >> t;
             wy_e_input[i] = t;
             Dinitrogen.Set_wy_e(wy_e_input);
         }
         wy_e_data.close();
     }else
         cout<<"Файл wy_e.txt не открылся" << endl;
     
    //B_e
    double B_e_input[11];
     ifstream B_e_data("data for N2/B_e.txt");
     if(B_e_data)
     {
         int n;
         double t;
         B_e_data >> n;
         for (int i = 0; i < n; i++) {
             B_e_data >> t;
             B_e_input[i] = t;
             Dinitrogen.Set_B_e(B_e_input);
         }
         B_e_data.close();
     }else
         cout<<"Файл B_e.txt не открылся" << endl;
     
    //D_e
    double D_e_input[11];
     ifstream D_e_data("data for N2/D_e.txt");
     if(D_e_data)
     {
         int n;
         double t;
         D_e_data >> n;
         for (int i = 0; i < n; i++) {
             D_e_data >> t;
             D_e_input[i] = t;
             Dinitrogen.Set_D_e(D_e_input);
         }
         D_e_data.close();
     }else
         cout<<"Файл D_e.txt не открылся" << endl;
    
    //b_e
    double b_e_input[11];
     ifstream b_e_data("data for N2/beta_e.txt");
     if(b_e_data)
     {
         int n;
         double t;
         b_e_data >> n;
         for (int i = 0; i < n; i++) {
             b_e_data >> t;
             b_e_input[i] = t;
             Dinitrogen.Set_b_e(b_e_input);
         }
         b_e_data.close();
     }else
         cout<<"Файл beta_e.txt не открылся" << endl;
     
    //e_diss
    double e_diss_input[11];
     ifstream e_diss_data("data for N2/e_diss.txt");
     if(e_diss_data)
     {
         int n;
         double t;
         e_diss_data >> n;
         for (int i = 0; i < n; i++) {
             e_diss_data >> t;
             e_diss_input[i] = t;
             Dinitrogen.Set_e_diss(e_diss_input);
         }
         e_diss_data.close();
     }else
         cout<<"Файл e_diss.txt не открылся" << endl;
    
    //m
    double m_input;
     ifstream m_data("data for N2/m.txt");
     if(m_data)
     {
         int n;
         double t;
         m_data >> n;
         for (int i = 0; i < n; i++) {
             m_data >> t;
             m_input = t;
             Dinitrogen.Set_m(m_input);
         }
         m_data.close();
     }else
         cout<<"Файл m.txt не открылся" << endl;
    
    
    //Экземпляр класса атома углерода C
    Atom Carbon;
    
    //e_el
    double e_el_Atom_input[431];
     ifstream e_el_Atom_data("data for C/e_el.txt");
     if(e_el_Atom_data)
     {
         int n;
         double t;
         e_el_Atom_data >> n;
         for (int i = 0; i < n; i++) {
             e_el_Atom_data >> t;
             e_el_Atom_input[i] = t;
             Carbon.Set_e_el(e_el_Atom_input);
         }
         e_el_Atom_data.close();
     }else
         cout<<"Файл e_el.txt не открылся" << endl;

    
    //g_el
    double g_el_Atom_input[431];
     ifstream g_el_Atom_data("data for C/g_el.txt");
     if(g_el_Atom_data)
     {
         int n;
         double t;
         g_el_Atom_data >> n;
         for (int i = 0; i < n; i++) {
             g_el_Atom_data >> t;
             g_el_Atom_input[i] = t;
             Carbon.Set_g_el(g_el_Atom_input);
         }
         g_el_Atom_data.close();
     }else
         cout<<"Файл g_el.txt не открылся" << endl;
    
    //m
    double m_Atom_input;
     ifstream m_Atom_data("data for C/m.txt");
     if(m_Atom_data)
     {
         int n;
         double t;
         m_Atom_data >> n;
         for (int i = 0; i < n; i++) {
             m_Atom_data >> t;
             m_Atom_input = t;
             Carbon.Set_m(m_Atom_input);
         }
         m_Atom_data.close();
     }else
         cout<<"Файл m.txt не открылся" << endl;
    
    double T;
    cout << "Введите T: "<< endl;
    cin >> T;
    cout <<"N2:" <<endl;
    cout <<"Поступательная стат. сумма для N2: " << endl;
    cout << Dinitrogen.Z_tr(T, "Single")[0] << endl;
    cout <<"Внутренняя стат. сумма для N2: " << endl;
    cout << Dinitrogen.Z_int(T, "Single")[0] << endl;
    cout <<"Полная стат. сумма для N2: " << endl;
    cout << Dinitrogen.Z(T, "Single")[0] << endl;
    cout <<"Поступательная удельная энергия для N2: " << endl;
    cout << Dinitrogen.E_tr(T, "Single")[0] << endl;
    cout <<"Внутренняя удельная энергия для N2: " << endl;
    cout << Dinitrogen.E_int(T, "Single")[0] << endl;
    cout <<"Полная удельная энергия для N2: " << endl;
    cout << Dinitrogen.E(T, "Single")[0] << endl;
    cout <<"Поступательная удельная теплоёмкость при постоянном объёме для N2: " << endl;
    cout << Dinitrogen.c_V_tr(T, "Single")[0] << endl;
    cout <<"Внутренняя удельная теплоёмкость при постоянном объёме для N2: " << endl;
    cout << Dinitrogen.c_V_int(T, "Single")[0] << endl;
    cout <<"Полная удельная теплоёмкость при постоянном объёме для N2: " << endl;
    cout << Dinitrogen.c_V(T, "Single")[0] << endl;
    cout <<"C:" <<endl;
    cout <<"Поступательная стат. сумма для C: " << endl;
    cout << Carbon.Z_tr(T, "Single")[0] << endl;
    cout <<"Внутренняя стат. сумма для C: " << endl;
    cout << Carbon.Z_int(T, "Single")[0] << endl;
    cout <<"Полная стат. сумма для C: " << endl;
    cout << Carbon.Z(T, "Single")[0] << endl;
    cout <<"Поступательная удельная энергия для C: " << endl;
    cout << Carbon.E_tr(T, "Single")[0] << endl;
    cout <<"Внутренняя удельная энергия C: " << endl;
    cout << Carbon.E_int(T, "Single")[0] << endl;
    cout <<"Полная удельная энергия C: " << endl;
    cout << Carbon.E(T, "Single")[0] << endl;
    cout <<"Поступательная удельная теплоёмкость при постоянном объёме для C: " << endl;
    cout << Carbon.c_V_tr(T, "Single")[0] << endl;
    cout <<"Внутренняя удельная теплоёмкость при постоянном объёме для C: " << endl;
    cout << Carbon.c_V_int(T, "Single")[0] << endl;
    cout <<"Полная удельная теплоёмкость при постоянном объёме для C: " << endl;
    cout << Carbon.c_V(T, "Single")[0] << endl;
    return 0;
}
