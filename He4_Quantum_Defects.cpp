//===========================================================================================================
// Author: Jens Chluba
// first implementation: June   2007
// last modification   : August 2015
//===========================================================================================================
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

#include "physical_consts.h"
#include "He4_Quantum_Defects.h"

using namespace std;

//===========================================================================================================
// Helium constants
//===========================================================================================================
//const double const_RM=3289391007.44e+6;   // Rydberg constant for Helium in [Hz] (from Drake book)
//const double const_mu_M=1.37074562e-4;    // ration mu/M for 4^He                (from Drake book)
const double const_RM=3289391006.715e+6;    // Rydberg constant for Helium in [Hz] (from Drake book)
const double const_mu_M=1.370745641e-4;     // ratio mu/M for 4^He                 (from Drake book)

//===========================================================================================================
double Drake_DW(int n)
{
    return const_RM*(-0.75*pow(const_alpha/n/n, 2)
                     +pow(const_mu_M/n, 2)*(1.0+5.0/6.0*pow(const_alpha*2.0, 2)));
}

double QD_d(double n, double d, double d0, double d2, double d4, double d6, double d8)
{
    double nmd=n-d;
    return (d0+d2/pow(nmd, 2)+d4/pow(nmd, 4)+d6/pow(nmd, 6)+d8/pow(nmd, 8));
}

//===========================================================================================================
// returns the energy of the level in eV
//===========================================================================================================
double compute_DEc_QD(int n, int l, int s, int j)
{
    if(l>6)
    {
        cout << " this level energy cannot be obtained using quantum defects: "
        << n << " " << l << " " << s << " " << j << endl;
        return 0.0;
    }
    
    // load the quantum defect coefficients
    double d0, d2, d4, d6, d8;
    int index;
    if(s==0)
    {
        index=l*5;
        d0=DrakeQD_SLL[index+0].QD;
        d2=DrakeQD_SLL[index+1].QD;
        d4=DrakeQD_SLL[index+2].QD;
        d6=DrakeQD_SLL[index+3].QD;
        d8=DrakeQD_SLL[index+4].QD;
    }
    else if(s==1 && j==l-1 && l>0)
    {
        index=(l-1)*5;
        d0=DrakeQD_TLLm1[index+0].QD;
        d2=DrakeQD_TLLm1[index+1].QD;
        d4=DrakeQD_TLLm1[index+2].QD;
        d6=DrakeQD_TLLm1[index+3].QD;
        d8=DrakeQD_TLLm1[index+4].QD;
    }
    else if(s==1 && j==l && l>0)
    {
        index=(l-1)*5;
        d0=DrakeQD_TLL[index+0].QD;
        d2=DrakeQD_TLL[index+1].QD;
        d4=DrakeQD_TLL[index+2].QD;
        d6=DrakeQD_TLL[index+3].QD;
        d8=DrakeQD_TLL[index+4].QD;
    }
    else if(s==1 && j==l+1)
    {
        index=l*5;
        d0=DrakeQD_TLLp1[index+0].QD;
        d2=DrakeQD_TLLp1[index+1].QD;
        d4=DrakeQD_TLLp1[index+2].QD;
        d6=DrakeQD_TLLp1[index+3].QD;
        d8=DrakeQD_TLLp1[index+4].QD;
    }
    else{ cerr << " compute_DEc_QD:: oops " << endl; return 0.0; }
    
    double ns=n, nso, Dn=1.0e+300;
    double eps=1.0e-15;
    
    while(fabs(Dn/ns)>eps)
    {
        nso=ns;
        ns=n-QD_d(n, n-ns, d0, d2, d4, d6, d8);
        Dn=ns-nso;
    }
    
    return (const_RM/ns/ns+Drake_DW(n))*const_h/const_e;
}

//===========================================================================================================
void write_DEc_QD(int nmax, string path)
{
    string name=path+"./DEc_QD_Helium.dat";
    ofstream ofile(name.c_str());
    ofile.precision(16);
    double d;
    
    cout.precision(16);
    
    // Singlet-states
    for(int n=2; n<=nmax; n++)
        for(int l=0; l<(int) min(n, 7); l++)
        {
            d=compute_DEc_QD(n, l, 0, l);
            cout << n << " " << l << " " << 0 << " " << l << "    " << scientific << d << endl;
            ofile << n << " " << l << " " << 0 << " " << l << "    " << scientific << d << endl;
        }
    
    // Triplet-states
    for(int n=2; n<=nmax; n++)
        for(int l=0; l<(int) min(n, 7); l++)
            for(int j=l-1; j<=l+1; j++)
            {
                if(l==0) j=0+1;
                d=compute_DEc_QD(n, l, 1, j);
                cout << n << " " << l << " " << 1 << " " << j << "    " << scientific << d << endl; 
                ofile << n << " " << l << " " << 1 << " " << j << "    " << scientific << d << endl; 
            }
    
    ofile.close();
    return;
}
//===========================================================================================================
//===========================================================================================================
