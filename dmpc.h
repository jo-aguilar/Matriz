#ifndef DMPC_H
#define DMPC_H
#include "matriz.h"

double malha_fechada(int Nc, int Np, Matriz* F, Matriz* Phi, Matriz* R_barra, 
		   Matriz* Rs_barra, double r_ki, double rw, Matriz* x_ki);
void func_deltau(Matriz* Phi, Matriz* R_barra, Matriz* Rs_barra, Matriz* F, double r_ki, Matriz* x_ki, Matriz* ret); //pronto
void func_Phi(Matriz* A, Matriz* B, Matriz* C, int Np, int Nc, Matriz* dest); //pronto
void func_F(Matriz* A, Matriz* C, Matriz* dest, int Np); //pronto
void ret_A(Matriz* Am, Matriz* Cm, Matriz* A); //pronto
void ret_B(Matriz* Bm, Matriz* Cm, Matriz* B); //pronto
void ret_C(Matriz* Cm, Matriz* C); //pronto 
void discret(double, Matriz*, Matriz*, Matriz*, Matriz*); //pronto

#endif
