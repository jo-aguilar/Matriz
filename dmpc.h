#ifndef DMPC_H
#define DMPC_H

#include "matriz.h"

Matriz ret_A(Matriz* Am, Matriz* Cm);
Matriz ret_B(Matriz* Bm, Matriz* Cm);
Matriz ret_C(Matriz* Cm);
void discret(double, Matriz*, Matriz*, Matriz*, Matriz*);

#endif
