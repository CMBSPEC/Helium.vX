//===========================================================================================================
// Author: Jens Chluba
// first implementation: June   2007
// last modification   : August 2015
//===========================================================================================================
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "physical_consts.h"
#include "Definitions.h"
#include "routines.h"
#include "File.h"
#include "Patterson.h"
#include "HeI_Ric_Topbase.h"

using namespace std;

int mess_Topbase=0;
double nuc_fac_TopBase=1.0e+10;   // old calculations used == 2.0 and ==10.0
int npol=2;

//===========================================================================================================
// convert the crossections into Ric
//===========================================================================================================
// some data
//===========================================================================================================
struct TopBase_Level
{
  int n;
  int l;
  double Ec;
  string fname;
};

TopBase_Level TopBase_Levels_S[25]=
{
    // S-states
    {1, 0, 1.78678E+00, "sig_c.HeIS.1s.dat"},
    {2, 0, 2.89696E-01, "sig_c.HeIS.2s.dat"},
    {3, 0, 1.21938E-01, "sig_c.HeIS.3s.dat"},
    {4, 0, 6.69194E-02, "sig_c.HeIS.4s.dat"},
    {5, 0, 4.22279E-02, "sig_c.HeIS.5s.dat"},
    {6, 0, 2.90552E-02, "sig_c.HeIS.6s.dat"},
    {7, 0, 2.12076E-02, "sig_c.HeIS.7s.dat"},
    {8, 0, 1.61581E-02, "sig_c.HeIS.8s.dat"},
    {9, 0, 1.27188E-02, "sig_c.HeIS.9s.dat"},
    {10, 0, 1.02713E-02, "sig_c.HeIS.10s.dat"},
    // P-states
    {2, 1, 2.47481E-01, "sig_c.HeIS.2p.dat"},
    {3, 1, 1.10225E-01, "sig_c.HeIS.3p.dat"},
    {4, 1, 6.21098E-02, "sig_c.HeIS.4p.dat"},
    {5, 1, 3.97964E-02, "sig_c.HeIS.5p.dat"},
    {6, 1, 2.76588E-02, "sig_c.HeIS.6p.dat"},
    {7, 1, 2.03328E-02, "sig_c.HeIS.7p.dat"},
    {8, 1, 1.55743E-02, "sig_c.HeIS.8p.dat"},
    {9, 1, 1.23100E-02, "sig_c.HeIS.9p.dat"},
    //{10, 1, 9.97393E-03, "sig_c.HeIS.10p.dat"},
    // D-states
    {3, 2, 1.11240E-01, "sig_c.HeIS.3d.dat"},
    {4, 2, 6.25585E-02, "sig_c.HeIS.4d.dat"},
    {5, 2, 4.00308E-02, "sig_c.HeIS.5d.dat"},
    {6, 2, 2.77960E-02, "sig_c.HeIS.6d.dat"},
    {7, 2, 2.04198E-02, "sig_c.HeIS.7d.dat"},
    {8, 2, 1.56328E-02, "sig_c.HeIS.8d.dat"},
    {9, 2, 1.23512E-02, "sig_c.HeIS.9d.dat"}
    //{10, 2, 1.00041E-02, "sig_c.HeIS.10d.dat"}
};

TopBase_Level TopBase_Levels_T[18]=
{
    // S-states
    {2, 0, 3.50319E-01, "sig_c.HeIT.2s.dat"},
    {3, 0, 1.37361E-01, "sig_c.HeIT.3s.dat"},
    {4, 0, 7.30168E-02, "sig_c.HeIT.4s.dat"},
    {5, 0, 4.52344E-02, "sig_c.HeIT.5s.dat"},
    {6, 0, 3.07531E-02, "sig_c.HeIT.6s.dat"},
    {7, 0, 2.22588E-02, "sig_c.HeIT.7s.dat"},
    {8, 0, 1.68535E-02, "sig_c.HeIT.8s.dat"},
    {9, 0, 1.32025E-02, "sig_c.HeIT.9s.dat"},
    {10, 0, 1.06212E-02, "sig_c.HeIT.10s.dat"},
    // P-states
    {2, 1, 2.66159E-01, "sig_c.HeIT.2p.dat"},
    {3, 1, 1.16111E-01, "sig_c.HeIT.3p.dat"},
    //{4, 1, 6.46272E-02, "sig_c.HeIT.4p.dat"},
    //{5, 1, 4.10912E-02, "sig_c.HeIT.5p.dat"},
    //{6, 1, 2.84094E-02, "sig_c.HeIT.6p.dat"},
    //{7, 1, 2.08058E-02, "sig_c.HeIT.7p.dat"},
    //{8, 1, 1.58913E-02, "sig_c.HeIT.8p.dat"},
    //{9, 1, 1.25326E-02, "sig_c.HeIT.9p.dat"},
    //{10, 1, 1.01362E-02, "sig_c.HeIT.10p.dat"},
    // D-states
    {3, 2, 1.11272E-01, "sig_c.HeIT.3d.dat"},
    {4, 2, 6.25767E-02, "sig_c.HeIT.4d.dat"},
    {5, 2, 4.00413E-02, "sig_c.HeIT.5d.dat"},
    {6, 2, 2.78024E-02, "sig_c.HeIT.6d.dat"},
    {7, 2, 2.04240E-02, "sig_c.HeIT.7d.dat"},
    {8, 2, 1.56357E-02, "sig_c.HeIT.8d.dat"},
    {9, 2, 1.23533E-02, "sig_c.HeIT.9d.dat"}
    //{10, 2, 1.00056E-02, "sig_c.HeIT.10d.dat"}
};

struct TopBase_Level_data
{
    vector<double> lgnu_Hz;
    vector<double> lgsig_ic;
};

vector<TopBase_Level_data> TopBase_Level_data_S;
vector<TopBase_Level_data> TopBase_Level_data_T;

//===========================================================================================================
//
// load Topbase data
//
//===========================================================================================================
void load_Topbase_data(string fname, double Ec, vector<TopBase_Level_data> &TD)
{
    if(mess_Topbase>=1)
    {
        cout << " *****************************************************************" << endl;
        cout << " reading from file: " << fname << endl;
        
        Ec*=const_Ry_inf_icm*const_cl;                         // conversion E/Ryd --> Hz
        cout << " Ec: " << Ec << endl;
    }
    
    ifstream file(fname.c_str());
    string str;
    TopBase_Level_data dum;
    
    char *pEnd;
    double d1, d2;
    
    // read the data
    while(!file.eof())
    {
        getline(file, str);
        
        d1= strtod(str.c_str(), &pEnd);
        d2= strtod(pEnd,NULL);
        
        dum.lgnu_Hz.push_back(log(d1*const_Ry_inf_icm*const_cl));  // conversion E/Ryd --> Hz
        dum.lgsig_ic.push_back(log(d2*1.0e-18));                   // conversion Mbarn --> cm^2
    }
    
    TD.push_back(dum);
    file.close();
    
    return;
}

void load_all_Topbase_data(string path)
{
    TopBase_Level_data_S.clear();
    TopBase_Level_data_T.clear();
    
    // Singlet-states
    for(int i=0; i<25; i++)
        load_Topbase_data(path+TopBase_Levels_S[i].fname, TopBase_Levels_S[i].Ec, TopBase_Level_data_S);

    // Triplet-states
    for(int i=0; i<18; i++)
        load_Topbase_data(path+TopBase_Levels_T[i].fname, TopBase_Levels_T[i].Ec, TopBase_Level_data_T);

    return;
}

//===========================================================================================================
//
// compute Ric
//
//===========================================================================================================
struct TopBase_integration_data
{
    TopBase_Level_data *TDp;
    int np;
    double Tg;
};

double N_nu_pl_Topbase(double nu, double Tg)
{ return 2.0*pow(nu/const_cl, 2)/( exp( min(700.0, const_h_kb*nu/Tg) )-1.0); }

double dRic_lin(double lgnu, void *p)
{
    TopBase_integration_data *d=((TopBase_integration_data *) p);
    
    double nu=exp(lgnu), y=0, dy=0;
    polint_JC(&(d->TDp->lgnu_Hz[0]), &(d->TDp->lgsig_ic[0]), d->np, lgnu, npol, &y, &dy);

    return nu*exp(y)*N_nu_pl_Topbase(nu, d->Tg);
}

double integrate_Ric(TopBase_Level_data &TD, double Ec, double Tg)
{
    TopBase_integration_data d;
    d.Tg=Tg;
    d.TDp=&TD;
    d.np=TD.lgnu_Hz.size();
    void *p=&d;
    
    double r=0.0;
    double epsrel=1.0e-5, epsabs=1.0e-20;
    double nuc=Ec*const_Ry_inf_icm*const_cl, num=min(nuc*nuc_fac_TopBase, exp(TD.lgnu_Hz.back()));
    double ab=log(nuc), bb=log(num);
    
    r=Integrate_using_Patterson_adaptive(ab, bb, epsrel, epsabs, dRic_lin, p);
    
    return FOURPI*r;    
}

//===========================================================================================================
//
// photoionization rate in 1/sec
//
//===========================================================================================================
double Ric_Topbase(int n, int l, int s, double Tg)
{
    int i;
    if(s==0)
    {
        for(i=0; i<(int)TopBase_Level_data_S.size(); i++)
            if(TopBase_Levels_S[i].n==n && TopBase_Levels_S[i].l==l) break;
        
        return integrate_Ric(TopBase_Level_data_S[i], TopBase_Levels_S[i].Ec, Tg);
    }
    else if(s==1)
    {
        for(i=0; i<(int)TopBase_Level_data_T.size(); i++)
            if(TopBase_Levels_T[i].n==n && TopBase_Levels_T[i].l==l) break;
        
        return integrate_Ric(TopBase_Level_data_T[i], TopBase_Levels_T[i].Ec, Tg);
    }
    else{ cerr << " no Topbase data found for (n, l, s) == ("
               << n << ", " << l << ", " << s << ")! Exiting... " << endl; exit(0); }
    
    return 0.0;
}

//===========================================================================================================
//
// photoionization cross section
//
//===========================================================================================================
double sig_ic_Topbase(int n, int l, int s, double nu)
{
    int i;
    double y=0, dy=0;
    
    if(s==0)
    {
        for(i=0; i<(int)TopBase_Level_data_S.size(); i++)
            if(TopBase_Levels_S[i].n==n && TopBase_Levels_S[i].l==l) break;
        
        if(TopBase_Levels_S[i].Ec*const_Ry_inf_icm*const_cl>nu) return 0;
        
        polint_JC(&(TopBase_Level_data_S[i].lgnu_Hz[0]),
                  &(TopBase_Level_data_S[i].lgsig_ic[0]),
                  TopBase_Level_data_S[i].lgnu_Hz.size(),
                  log(nu), npol, &y, &dy);
        
        return exp(y);
    }
    else if(s==1)
    {
        for(i=0; i<(int)TopBase_Level_data_T.size(); i++)
            if(TopBase_Levels_T[i].n==n && TopBase_Levels_T[i].l==l) break;
        
        if(TopBase_Levels_T[i].Ec*const_Ry_inf_icm*const_cl>nu) return 0;

        polint_JC(&(TopBase_Level_data_T[i].lgnu_Hz[0]),
                  &(TopBase_Level_data_T[i].lgsig_ic[0]),
                  TopBase_Level_data_T[i].lgnu_Hz.size(),
                  log(nu), npol, &y, &dy);
        
        return exp(y);
    }
    else{ cerr << " no Topbase data found for (n, l, s) == ("
               << n << ", " << l << ", " << s << ")! Exiting... " << endl; exit(0); }
    
    return 0.0;
}

//===========================================================================================================
//===========================================================================================================
