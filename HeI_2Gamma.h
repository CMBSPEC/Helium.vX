//-------------------------------------------------------------------
// Author: Jens Chluba 
// Date: June 2007
//-------------------------------------------------------------------
#ifndef HEI_2GAMMA_H
#define HEI_2GAMMA_H

#include <vector>
#include "physical_consts.h"
#include "Cosmos.h"
#include "File.h"

using namespace std;

double phi_HeI_2g(double nu);
double phi_HeI_2s_decay_spectrum(double y)
{
	if(y<0.0000001 || y> 0.9999999) return 0.0;
    double w=y*(1.0-y);
    
    return 19.6039602*pow(w, 3)*(1.742-7.2*w+12.8*w*w)/pow(w+0.03, 2);
}

double phi_HeI_2g_ind(double nu, double Tg);

#endif

