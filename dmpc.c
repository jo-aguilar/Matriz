#include "dmpc.h"

double malha_fechada (int Nc, int Np, Matriz* F, Matriz* Phi, Matriz* R_barra, 
		      Matriz* Rs_barra, double r_ki, double rw, Matriz* x_ki){
//Malha fechada para execução real, não manipulando entrada por meio de ganhos,
//enviando para fora da função somente o sinal de controle
//[!!!] DEVE SER TESTADO POSTERIORMENTE	
	static double u_anterior = NAN;
	static Matriz* x_anterior;
	if(isnan(u_anterior)){
		x_anterior = x_ki;
		u_anterior = 0;
	}
	Matriz Delta_U = MATRIZ_NULA(); 
	func_deltau(Phi, R_barra, Rs_barra, F, r_ki, x_ki, &Delta_U);
	double du_calc = u_anterior + Delta_U.ret(&Delta_U, 0, 0);
	
	if(du_calc >0.5)
		du_calc = 0.5;
	else if(du_calc < 0)
		du_calc = 0;
	
	u_anterior = du_calc;
	//não se está fazendo update de x_ki, aparentemente
	return du_calc;
}


void func_deltau(Matriz* Phi, Matriz* R_barra, Matriz* Rs_barra,
	       Matriz* F, double r_ki, Matriz* x_ki, Matriz* ret){
//Obtenção do incremento do sinal de controle Delta_U
	Matriz Phi_T     = MATRIZ_NULA(); transposta(Phi, &Phi_T);
	Matriz Phi_T_Phi = MATRIZ_NULA(); mat_mult(&Phi_T, Phi, &Phi_T_Phi);
	Matriz M0        = MATRIZ_NULA(); sm(&Phi_T_Phi, R_barra, &M0);
	Matriz M1        = MATRIZ_NULA(); inv(&M0, &M1);
	Matriz M2        = MATRIZ_NULA(); mat_mult(&M1, &Phi_T, &M2);
	Matriz Rs_barra_r_ki = MATRIZ_NULA(); mult(Rs_barra, &Rs_barra_r_ki, r_ki);
	Matriz Fx_ki     = MATRIZ_NULA(); mat_mult(F, x_ki, &Fx_ki);
	Matriz M3        = MATRIZ_NULA(); sub(&Rs_barra_r_ki, &Fx_ki, &M3);
	mat_mult(&M2, &M3, ret);
	//inv(Phi_T*Phi + R_barra)*Phi^T*(Rs_barra*r_ki - F*x_ki)
}

void func_Phi(Matriz* A, Matriz* B, Matriz* C, int Np, int Nc, Matriz* dest){
//Matriz deslizada de Toeplitz formada a partir de um vetor de fatores de
//potência a serem deslizados

	Matriz Phi = matriz((double[BUFF]){0}, 1, Np);
	Matriz Phi_primario = matriz((double[BUFF]){}, 1, Np);
	Matriz Temp1 = MATRIZ_NULA(),
	       Temp2 = MATRIZ_NULA(),
	       Temp3 = MATRIZ_NULA();
	int indice = Np - 1;

	//Crição do vetor (1xNp) a ser deslizado para a criação de Toeplitz
	for(int i = 1; i <= Np; i++){
		if((Np-i) != 0) {
			int pot = Np - i;
			mat_pot(A, &Temp1, pot);     //Temp1.print(&Temp1);
			mat_mult(C, &Temp1, &Temp2); //Temp2.print(&Temp2);
			mat_mult(&Temp2, B, &Temp3); //Temp3.print(&Temp3);
			  
		} else {
			mat_mult(C, B, &Temp3); 
			}
		Phi_primario.matriz[indice--] = Temp3.ret(&Temp3, 0, 0);
	}
	
	dest->linhas = Np;
	dest->colunas = Nc;
	int indice2 = 0;
	for(int i = 0; i < Np; i++){
		int indice3 = i;
		for(int j = 0; j < Nc; j++){
			if(indice3>=0){
				dest->matriz[indice2++] = 
					Phi_primario.ret(&Phi_primario, 0, indice3--); 
			} else {
				dest->matriz[indice2++] = 0;
			}
		}
		
	}
}

void func_F(Matriz* A, Matriz* C, Matriz* dest, int Np){
//Entrega a matriz compacta F (representativa das operações inerentes
//ao sistema
	Matriz M0 = MATRIZ_NULA();
	Matriz M1 = MATRIZ_NULA();
	Matriz M2 = MATRIZ_NULA();
	Matriz F  = MATRIZ_NULA();
	dest->linhas = Np;
	dest->colunas = 2;

	int indice = 0;
	for(int i = 1; i <= Np; i++){
		mat_pot(A, &M0, i);    
		mat_mult(C, &M0, &M1); 
		transposta(&M1, &M2);  
		dest->matriz[indice++] = M2.ret(&M2, 0, 0);
		dest->matriz[indice++] = M2.ret(&M2, 0, 1); 
	}
}

void ret_A(Matriz* Am, Matriz* Cm, Matriz* A){
//Entrega da versão estendida de A para CPDM
	Matriz M0 = zeros(Am->linhas, 1);
	Matriz M1   = zeros(Am->linhas, 1);                   //M1.print(&M1);       
	Matriz M2   = MATRIZ_NULA(); extend(Am, &M1, &M2);    //M2.print(&M2);              
	Matriz CmAm = MATRIZ_NULA(); mat_mult(Cm, Am, &CmAm); //CmAm.print(&CmAm);           
	Matriz Intg = MATRIZ_NULA(); ident(&Intg, 1);
	Matriz M3   = MATRIZ_NULA(); extend(&CmAm, &Intg, &M3); //M3.print(&M3);
	Matriz M4 = MATRIZ_NULA(); transposta(&M2, &M4); //M4.print(&M4); 
	Matriz M5 = MATRIZ_NULA(); transposta(&M3, &M5); //M5.print(&M5);
	Matriz M6 = MATRIZ_NULA(); extend(&M4, &M5, &M6); //M6.print(&M6);
	transposta(&M6, A);
}


void ret_B(Matriz* Bm, Matriz* Cm, Matriz* B){
//Entrega da versão estendida de B para CPDM
	Matriz CmBm = MATRIZ_NULA(); mat_mult(Cm, Bm, &CmBm);
	Matriz M0   = MATRIZ_NULA(); transposta(Bm, &M0);
	Matriz M1   = MATRIZ_NULA(); extend(&M0, &CmBm, &M1);
	transposta(&M1, B);
}


void ret_C(Matriz* Cm, Matriz* C){
//Entrega da versão estendida de C para CPDM
	Matriz I = MATRIZ_NULA(); ident(&I, 1);
	Matriz M0 = zeros(1, Cm->colunas);
	extend(&M0, &I, C);
}


void discret(double T, Matriz* Am, Matriz* Bm, Matriz* Ad, Matriz* Bd){
//Modifica os valores internos das matrizes Bd, e Ad a partir da discretização
//de Tustin feita sobre as Matrizes Am, Bm e sobre o valor de tempo amostrado T
	Matriz I = MATRIZ_NULA(); ident(&I, Am->colunas);
	Matriz A_modif = MATRIZ_NULA(); mult(Am, &A_modif, T/2);
	Matriz M0 = MATRIZ_NULA(); sub(&I, &A_modif, &M0);
	Matriz M1 = MATRIZ_NULA(); inv(&M0, &M1);          
	Matriz M2 = MATRIZ_NULA(); sm(&I, &A_modif, &M2);  
	mat_mult(&M1, &M2, Ad);

	Matriz M3 = MATRIZ_NULA(); mult(Bm, &M3, T);
	mat_mult(&M1, &M3, Bd);

}

