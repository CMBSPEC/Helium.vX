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
const double nuc_fac_TopBase=1.0e+10;   // old calculations used == 2.0 and ==10.0

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

const int nTopBase_Levels_S=sizeof(TopBase_Levels_S)/sizeof(TopBase_Levels_S[0]);
const int nTopBase_Levels_T=sizeof(TopBase_Levels_T)/sizeof(TopBase_Levels_T[0]);

static void check_Topbase_data_loaded(const vector<TopBase_Level_data> &TD, int expected, int s)
{
    if((int)TD.size()!=expected)
    {
        cerr << " Topbase data for s=" << s << " has not been loaded correctly. "
             << "Expected " << expected << " tables but found " << TD.size()
             << ". Exiting... " << endl;
        exit(0);
    }

    return;
}

static int find_Topbase_level_index(int n, int l, int s)
{
    const TopBase_Level *levels=NULL;
    int nlevels=0;

    if(s==0)
    {
        check_Topbase_data_loaded(TopBase_Level_data_S, nTopBase_Levels_S, s);
        levels=TopBase_Levels_S;
        nlevels=nTopBase_Levels_S;
    }
    else if(s==1)
    {
        check_Topbase_data_loaded(TopBase_Level_data_T, nTopBase_Levels_T, s);
        levels=TopBase_Levels_T;
        nlevels=nTopBase_Levels_T;
    }
    else
    {
        cerr << " no Topbase data found for (n, l, s) == ("
             << n << ", " << l << ", " << s << ")! Exiting... " << endl;
        exit(0);
    }

    for(int i=0; i<nlevels; i++)
        if(levels[i].n==n && levels[i].l==l) return i;

    cerr << " no Topbase data found for (n, l, s) == ("
         << n << ", " << l << ", " << s << ")! Exiting... " << endl;
    exit(0);

    return 0;
}

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
    if(!file.is_open())
    {
        cerr << " could not open Topbase data file: " << fname << ". Exiting... " << endl;
        exit(0);
    }

    string str;
    TopBase_Level_data dum;

    char *pEnd, *pEnd2;
    double d1, d2;

    // read the data
    while(getline(file, str))
    {
        d1= strtod(str.c_str(), &pEnd);
        if(pEnd==str.c_str()) continue;

        d2= strtod(pEnd, &pEnd2);
        if(pEnd2==pEnd) continue;

        if(d1<=0.0 || d2<0.0)
        {
            cerr << " invalid Topbase data in file: " << fname << ". Exiting... " << endl;
            exit(0);
        }

        dum.lgnu_Hz.push_back(log(d1*const_Ry_inf_icm*const_cl));  // conversion E/Ryd --> Hz
        // some of the table values are ==0.0 so make sure to not have that...
        dum.lgsig_ic.push_back(log(d2*1.0e-18+1.0e-300));          // conversion Mbarn --> cm^2
    }

    if(dum.lgnu_Hz.size()<2)
    {
        cerr << " insufficient Topbase data in file: " << fname << ". Exiting... " << endl;
        exit(0);
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
    for(int i=0; i<nTopBase_Levels_S; i++)
        load_Topbase_data(path+TopBase_Levels_S[i].fname, TopBase_Levels_S[i].Ec, TopBase_Level_data_S);

    // Triplet-states
    for(int i=0; i<nTopBase_Levels_T; i++)
        load_Topbase_data(path+TopBase_Levels_T[i].fname, TopBase_Levels_T[i].Ec, TopBase_Level_data_T);

    return;
}

//===========================================================================================================
double interpolate_Topbase(double lgnu,
                               const vector<double> &lgnu_Hz,
                               const vector<double> &lgsig,
                               long unsigned int np)
{
    long unsigned int j=0;
    if(np<2) throw_error("interpolate_Topbase", "insufficient data", 1);

    // find index corresponding to x (start at i=istart)
    locate_JC(&lgnu_Hz[0], np, lgnu, &j);
    if(j==np-1) throw_error("interpolate_Topbase", "out of bound", 1);

    double dlnsig_dln_nu=(lgsig[j+1]-lgsig[j])/(lgnu_Hz[j+1]-lgnu_Hz[j]);

    return exp(lgsig[j]+dlnsig_dln_nu*(lgnu-lgnu_Hz[j]));
}

//===========================================================================================================
//
// photoionization cross section
//
//===========================================================================================================
double sig_ic_Topbase(int n, int l, int s, double nu)
{
    int i=find_Topbase_level_index(n, l, s);

    if(s==0)
    {
        if(TopBase_Levels_S[i].Ec*const_Ry_inf_icm*const_cl>nu) return 0;

        return interpolate_Topbase(log(nu),
                                TopBase_Level_data_S[i].lgnu_Hz,
                                TopBase_Level_data_S[i].lgsig_ic,
                                TopBase_Level_data_S[i].lgnu_Hz.size());
    }
    else
    {
        if(TopBase_Levels_T[i].Ec*const_Ry_inf_icm*const_cl>nu) return 0;

        return interpolate_Topbase(log(nu),
                                       TopBase_Level_data_T[i].lgnu_Hz,
                                       TopBase_Level_data_T[i].lgsig_ic,
                                       TopBase_Level_data_T[i].lgnu_Hz.size());
    }
    return 0.0;
}

//===========================================================================================================
//
// compute Ric
//
//===========================================================================================================
struct TopBase_integration_data
{
    const TopBase_Level_data *TDp;
    int np;
    double Tg, rho, nuc;
};

//===========================================================================================================
double dRic_lin(double lgnu, void *p)
{
    TopBase_integration_data *d=((TopBase_integration_data *) p);
    
    double nu=exp(lgnu);
    double sig_nu=interpolate_Topbase(lgnu, d->TDp->lgnu_Hz, d->TDp->lgsig_ic, d->np);
    double Tg=d->Tg, x=const_h_kb*nu/Tg, xc=const_h_kb*d->nuc/Tg;

    return nu*pow(nu/const_cl, 2)*sig_nu*exp(xc-x)/one_minus_exp_mx(x);
}

//===========================================================================================================
double integrate_Ric_Patt(const TopBase_Level_data &TD, double Ec, double nucHe, double Tg)
{
    TopBase_integration_data d;
    d.nuc=Ec*const_Ry_inf_icm*const_cl;
    d.TDp=&TD;
    d.np=TD.lgnu_Hz.size();

    double T_scale=d.nuc/nucHe;  // to fix the small mismatch in the level energies use Tg
    Tg*=T_scale;
    d.Tg=Tg;

    double r=0.0;
    double epsrel=1.0e-5, epsabs=1.0e-30;
    double nuc=d.nuc, num=min(nuc*nuc_fac_TopBase, exp(TD.lgnu_Hz.back()));
    double ab=nuc, bb=num;
    
    //-------------------------------------------------------------------
    // recombination cross section exponentially cuts of for x-xi >> 1
    //-------------------------------------------------------------------
    double nu_x=const_h_kb/Tg, xc=nu_x*d.nuc;
    if(xc>=3.0) bb=min(bb, (xc -log(1.0e-30))/nu_x);
    else bb=min(bb, (3.0-log(1.0e-30))/nu_x);

    r=Integrate_using_Patterson_adaptive(log(ab), log(bb), epsrel, epsabs, dRic_lin, &d);

    return 2.0*FOURPI*r*pow(T_scale, -3)*exp(-xc);
}

//===========================================================================================================
double integrate_Ric(const TopBase_Level_data &TD, double Ec, double nucHe, double Tg)
{ return integrate_Ric_Patt(TD, Ec, nucHe, Tg); }

//===========================================================================================================
//
// photoionization rate in 1/sec
//
//===========================================================================================================
double Ric_Topbase(int n, int l, int s, double nucHe, double Tg)
{
    int i=find_Topbase_level_index(n, l, s);
    if(s==0)
        return integrate_Ric(TopBase_Level_data_S[i], TopBase_Levels_S[i].Ec, nucHe, Tg);
    else
        return integrate_Ric(TopBase_Level_data_T[i], TopBase_Levels_T[i].Ec, nucHe, Tg);
    
    return 0.0;
}

//===========================================================================================================
//
// compute Rci
//
//===========================================================================================================
double dIci_lin(double lgnu, void *p)
{
    TopBase_integration_data *d=((TopBase_integration_data *) p);

    double nu=exp(lgnu);
    double sig_nu=interpolate_Topbase(lgnu, d->TDp->lgnu_Hz, d->TDp->lgsig_ic, d->np);
    double Tg=d->Tg, x=const_h_kb*nu/Tg, xc=const_h_kb*d->nuc/Tg;

    return nu*pow(nu/const_cl, 2)*sig_nu*exp((xc-x)/d->rho)*nbbp1_func(x);
}

//===========================================================================================================
double integrate_Ici_Patt(const TopBase_Level_data &TD, double Ec, double nucHe, double Tg, double rho)
{
    TopBase_integration_data d;
    d.nuc=Ec*const_Ry_inf_icm*const_cl;
    d.rho=rho; // == Te/Tg
    d.TDp=&TD;
    d.np=TD.lgnu_Hz.size();

    double T_scale=d.nuc/nucHe;  // to fix the small mismatch in the level energies use Tg
    Tg*=T_scale;
    d.Tg=Tg;

    double r=0.0;
    double epsrel=1.0e-5, epsabs=1.0e-30;
    double nuc=d.nuc, num=min(nuc*nuc_fac_TopBase, exp(TD.lgnu_Hz.back()));
    double ab=nuc, bb=num;

    //-------------------------------------------------------------------
    // recombination cross section exponentially cuts of for x-xi >> rho
    //-------------------------------------------------------------------
    double nu_x=const_h_kb/Tg, xc=nu_x*d.nuc;
    if(xc/rho>=3.0) bb=min(bb, (xc -rho*log(1.0e-30))/nu_x);
    else bb=min(bb, (3.0-rho*log(1.0e-30))/nu_x);

    r=Integrate_using_Patterson_adaptive(log(ab), log(bb), epsrel, epsabs, dIci_lin, &d);

    return 2.0*FOURPI*r*pow(T_scale, -3);
}

//===========================================================================================================
//
// recombination rate integral 8pi * int (nu/c)^2 sig_nu exp(xce-xe) (1+ng) d nu
//
//===========================================================================================================
double integrate_Ici(const TopBase_Level_data &TD, double Ec, double nucHe, double Tg, double rho)
{ return integrate_Ici_Patt(TD, Ec, nucHe, Tg, rho); }

//===========================================================================================================
//
// recombination rate integral int ...
//
//===========================================================================================================
double Ici_Topbase(int n, int l, int s, double nucHe, double Tg, double rho)
{
    int i=find_Topbase_level_index(n, l, s);
    if(s==0)
        return integrate_Ici(TopBase_Level_data_S[i], TopBase_Levels_S[i].Ec, nucHe, Tg, rho);
    else
        return integrate_Ici(TopBase_Level_data_T[i], TopBase_Levels_T[i].Ec, nucHe, Tg, rho);

    return 0.0;
}

// gT=f(T) exp(-xce)
double Rci_Topbase(int n, int l, int s, double gT, double nucHe, double Tg, double rho)
{ return gT*Ici_Topbase(n, l, s, nucHe, Tg, rho); }

//===========================================================================================================
//===========================================================================================================
