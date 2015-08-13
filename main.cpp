//===========================================================================================================
// Simple code to test helium atom model
//
// Author: Jens Chluba
// Date             : August 2015
// last modification: August 2015
//===========================================================================================================

//===========================================================================================================
// several libs
//===========================================================================================================
#include "HeI_Atom.h"
#include "He4_Quantum_Defects.h"

#include "physical_consts.h"
#include "routines.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

using namespace std;

//===========================================================================================================
int main(int narg, char *args[])
{
	if(!(narg==2))
	{
		cout << "error: usage $ ./Helium nShells" << endl;
		exit(0);
	}

	string str=args[narg-1];
    int nShells=atof(str.c_str());

    int njres=10;
    int nQ=10;
    int nTS=10;
    
    Gas_of_HeI_Atoms HeIA(nShells, njres, nQ, nTS, -2);
    
    for(int n=2; n<=nShells; n++)
        for(int l=0; l<=min(6, n-1); l++)
            for(int s=0; s<=1; s++)
            {
                double E1, E2;
                if(s==0)
                {
                    E1=compute_DEc_QD(n, l, s, l);
                    E2=HeIA.Get_nu_ion(HeIA.Get_Level_index(n, l, s, l))*const_h/const_e;
                    
                    cout << n << " " << l << " " << s << " " << l << " " << E1/E2-1.0
                              << " " << E1*const_e/const_h/1.0e+6 << " " << E2*const_e/const_h/1.0e+6 << endl;
               }
                else
                {
                    int jmin=l-1;
                    if(l==0) jmin=1;
                    
                    for(int j=jmin; j<=l+1; j++)
                    {
                        E1=compute_DEc_QD(n, l, s, j);
                        E2=HeIA.Get_nu_ion(HeIA.Get_Level_index(n, l, s, j))*const_h/const_e;
                        
                        
                        cout << n << " " << l << " " << s << " " << j << " " << E1/E2-1.0
                                  << " " << E1*const_e/const_h/1.0e+6 << " " << E2*const_e/const_h/1.0e+6 << endl;
                    }
                }
            }
    
    cout << " number of levels: " << HeIA.Get_total_number_of_Levels() << endl;
    cout << " number of transitions: " << HeIA.Get_n_down() << endl;
    wait_f_r();
    
    cout << " specific transition data " << endl;
    HeIA.Sing.Level(8, 0).display_all_downward_transitions(); wait_f_r();
    HeIA.Sing.Level(8, 5).display_all_downward_transitions(); wait_f_r();
    HeIA.Sing.Level(10, 5).display_all_downward_transitions(); wait_f_r();
    HeIA.Sing.Level(10, 8).display_all_downward_transitions(); wait_f_r();
    
    HeIA.Sing.Level(14, 8).display_all_downward_transitions(); wait_f_r();
    HeIA.Sing.Level(15, 2).display_all_downward_transitions(); wait_f_r();
    
    HeIA.Trip.Level(8, 0, 1).display_all_downward_transitions(); wait_f_r();
    HeIA.Trip.Level(8, 5, 4).display_all_downward_transitions(); wait_f_r();
    HeIA.Trip.Level(10, 5, 4).display_all_downward_transitions(); wait_f_r();
    HeIA.Trip.Level(10, 8, 7).display_all_downward_transitions(); wait_f_r();
    
    HeIA.Trip_no_j.Level(14, 7).display_all_downward_transitions(); wait_f_r();
    HeIA.Trip_no_j.Level(14, 8).display_all_downward_transitions(); wait_f_r();
    HeIA.Trip_no_j.Level(11, 7).display_all_downward_transitions(); wait_f_r();

    return 0;
}
//===========================================================================================================
//===========================================================================================================
