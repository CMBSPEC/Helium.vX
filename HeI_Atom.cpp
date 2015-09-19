//========================================================================================
// Author: Jens Chluba
// Date: June 2007
// last modification: August 2015
//========================================================================================
// 18.08.2015: Added photo-ionization cross section routines
// 13.08.2015: Added Ric setup. Smith rates work as before with gw/(2s+1)/(2l+1) factor.
//             Hydrogenic rates are more accurate now, since not only 10Ec is used for
//             integration. Topbase rates also up to Em instead of 2Ec.
// 12.08.2015: Added all the transition rate setups. Data from Drake & Morton is used for
//             levels with n<=10 (aside from the gaps in that data). Checked that the
//             number of transitions is the same as well as some explicit transition
//             values. This made the helium setup significantly faster + reduced the
//             data that is required.
// 01.08.2015: Got rid of quantum defect data files. These where also not as accurate.
//             Now things match the book of Drake very well. Checked that hydrogenic
//             energies are already pretty accurate at n>=5.
// 31.07.2015: - fixed issue with setting up several HeI atoms. This was related to global
//               variable 'HeI_Atom_njresolved'
//             - loading additional quadrupole and Triplet-Singlet transitions is now
//               controlled without global functions (more safe...)
// 30.07.2015: added reduced mass
//========================================================================================

#include <iostream>
#include <string>
#include <cmath>

#include "HeI_Atom.h"
#include "He4_Quantum_Defects.h"
#include "Oscillator_strength.h"
#include "Smits_He_recomb_data.h"
#include "HeI_Ric_Topbase.h"

#include "routines.h"
#include "physical_consts.h"
#include "File.h"
#include "Voigtprofiles.h"

using namespace std;

//===================================================================================
const bool check_Transition_Data_Tables=false;

//===================================================================================
// He 1s2 ionization potential in eV
//===================================================================================
const double He1s2_ion=198310.6690*const_cl*const_h/const_e;  

//===================================================================================
// location of the atomic model for helium
//===================================================================================
#ifdef HEIDATADIR
#define HEIDATA ((string)HEIDATADIR)
#else
#define HEIDATA ((string)"")
#endif

const string path=HEIDATA+"./Helium.Data/";

//===================================================================================
// Rydberg for hydrogenic levels in neutral helium
//===================================================================================
const double const_EHeI=const_EH_inf/(1.0+const_me_malp);

//===================================================================================
// for default data is loaded
//===================================================================================
static int HeI_Atom_read_transition_data=1;
const bool switch_Bauman=0;
const bool switch_Hydrogenic=0;

//===================================================================================
// files for transition rates and energies from Drake & Morton
//===================================================================================
const string nameA_SS=path+"He4_SS.dat";
const string nameA_TS=path+"He4_TS.dat";
const string nameA_TT=path+"He4_TT.dat";
const string nameA_ST=path+"He4_ST.dat";
const string nameA_add=path+"He4_add.dat";

//===================================================================================
// 25th May 2009: added other n^3 P_1 - 1^1 S_0 intecombination lines for 3<=n<=10
//===================================================================================
const string nameA_add_Int=path+"He4_add.intercombination_lines.dat";

//===================================================================================
// 30th May 2009: added n^1D_2-1^1S_0 quadrupole lines for 3<=n<=10 
//===================================================================================
const string nameA_add_Quad=path+"He4_Quadrupole_lines.dat";

//===================================================================================
// do not change these global variables!
//===================================================================================
int _HeI_add_Intercomb=0;
int _HeI_add_Quad=0;
Atom_HeI_Triplet *Trip_glob=NULL;

//===================================================================================
// Read the tables for transition-rates and energies
//===================================================================================
void skip_header(inputFile &f, int nn)
{
    string str;
    for(int i=0; i<nn; i++) str=f.get_next_line();
    
    return;
}

//===================================================================================
void read_transition_inf(string fname, int n, int l, int s, int j, 
                         vector<Transition_Data_HeI_A> &v)
{
    inputFile f(fname);
    skip_header(f, 7);
    
    int nv, sv, lv, jv;
    double dum;
    Transition_Data_HeI_A dd;
    for(int k=0; k<5; k++) dd.xxx[k]=0.0;
    
    while(!f.iseof())
    {
        f.get_next_int();
        // lower level
        dd.np=f.get_next_int();
        dd.sp=f.get_next_int();
        dd.lp=f.get_next_int();
        dd.jp=f.get_next_int();
        
        // upper level
        nv=f.get_next_int();
        sv=f.get_next_int();
        lv=f.get_next_int();
        jv=f.get_next_int();
        
        if(nv==n && sv==s && lv==l && jv==j) 
        {
            f.get_next();
            dum=f.get_next();
            dd.DE=f.get_next()-dum;
            dd.A21=f.get_next();
            
            // set derived values
            dd.Dnu=dd.DE*const_e/const_h;
            dd.gwp=(2*dd.jp+1);
            dd.lambda21=const_cl/dd.Dnu; // better use my natural constants
            
            v.push_back(dd);
        }
        
        if(!f.iseof()) f.get_next_line();
    }
    
    f.close();   
    return;
}

//===================================================================================
void read_level_energy_DM(string fname, int n, int l, int s, int j, double &E1)
{
    // open the file with transition information 
    inputFile f(fname);
    
    skip_header(f, 7);
    
    int nv, sv, lv, jv;
    while(!f.iseof())
    {
        f.get_next_int();
        // lower level
        nv=f.get_next_int();
        sv=f.get_next_int();
        lv=f.get_next_int();
        jv=f.get_next_int();
        
        if(nv==n && lv==l && sv==s && jv==j) 
        {
            // skip next entries
            f.get_next_int();
            f.get_next_int();
            f.get_next_int();
            f.get_next_int();
            
            f.get_next();
            E1=f.get_next();
            f.close();   
            return;
        }
        else 
        {
            // upper level
            nv=f.get_next_int();
            sv=f.get_next_int();
            lv=f.get_next_int();
            jv=f.get_next_int();
            
            if(nv==n && lv==l && sv==s && jv==j) 
            {
                // skip next entries
                f.get_next();
                f.get_next();
                E1=f.get_next();
                f.close();   
                return;
            }
        }
        
        if(!f.iseof()) f.get_next_line();
    }
    
    cout << " read_level_energy_DM:\n no data on level: (" 
         << n << ", " << l << ", " << s << ", " << j << ") found in file: " 
         << fname << endl;
    
    E1=-6000.0;
    
    f.close();   
    return;
}

//===================================================================================
void read_transition_add_Quadrupole(string fname, int n, int l, int s, int j, 
                                    vector<Transition_Data_HeI_A> &v, int mflag)
{
    if(mflag>=1) cout << " read_transition_add_Quadrupole::"
                      << " loading n^1 D_2 - 1^1 S_0 quadrupole transitions. " << endl;
    inputFile f(fname);
    skip_header(f, 7);
    
    int nv, sv, lv, jv;
    double dum;
    Transition_Data_HeI_A dd;
    for(int k=0; k<5; k++) dd.xxx[k]=0.0;
    
    while(!f.iseof())
    {
        f.get_next_int();
        // lower level
        dd.np=f.get_next_int();
        dd.sp=f.get_next_int();
        dd.lp=f.get_next_int();
        dd.jp=f.get_next_int();
        
        // upper level
        nv=f.get_next_int();
        sv=f.get_next_int();
        lv=f.get_next_int();
        jv=f.get_next_int();
        
        if(nv==n && sv==s && lv==l && jv==j) 
        {
            f.get_next();
            dum=f.get_next();
            dd.DE=f.get_next()-dum;
            dd.A21=f.get_next();
            double fv=f.get_next();
            
            // set derived values
            dd.Dnu=dd.DE*const_e/const_h;
            dd.gwp=(2*dd.jp+1);
            dd.lambda21=const_cl/dd.Dnu; // better use my natural constants

            dd.A21=FOURPI*2.0*const_PIe2_mec*fv/dd.lambda21/dd.lambda21*1.0/5.0;        
            
            if(mflag>=1) cout << " (" << nv << " " << lv << " " << sv << " " << jv << ") --> ("
                              << dd.np << " " << dd.lp << " " << dd.sp << " " << dd.jp << ") A="
                              << dd.A21 << " " << 1.0/(FOURPI*2.0*const_PIe2_mec*1.0e+16) << endl;

            v.push_back(dd);
        }
        
        if(!f.iseof()) f.get_next_line();
    }
    
    f.close();   
    return;
}

//===================================================================================
void read_transition_add_TS(string fname, int n, int l, int s, int j, 
                            vector<Transition_Data_HeI_A> &v, int mflag)
{
    if(mflag>=1) cout << " read_transition_add_TS::"
                      << " loading additional n^3 P_1 - 1^1 S_0 transitions. " << endl; 

    inputFile f(fname);
    skip_header(f, 7);
    
    int nv, sv, lv, jv;
    double dum;
    Transition_Data_HeI_A dd;
    for(int k=0; k<5; k++) dd.xxx[k]=0.0;
    
    while(!f.iseof())
    {
        f.get_next_int();
        // lower level
        dd.np=f.get_next_int();
        dd.sp=f.get_next_int();
        dd.lp=f.get_next_int();
        dd.jp=f.get_next_int();
        
        // upper level
        nv=f.get_next_int();
        sv=f.get_next_int();
        lv=f.get_next_int();
        jv=f.get_next_int();
        
        if(nv==n && sv==s && lv==l && jv==j) 
        {
            f.get_next();
            dum=f.get_next();
            dd.DE=f.get_next()-dum;
            dd.A21=f.get_next();

            // set derived values
            dd.Dnu=dd.DE*const_e/const_h;
            dd.gwp=(2*dd.jp+1);
            dd.lambda21=const_cl/dd.Dnu; // better use my natural constants
            
            if(mflag>=1) cout << " (" << nv << " " << lv << " " << sv << " " 
                              << jv << ") --> (" << dd.np << " " << dd.lp << " " 
                              << dd.sp << " " << dd.jp 
                              << ") A=" << dd.A21 << endl;

            v.push_back(dd);
        }
        
        if(!f.iseof()) f.get_next_line();
    }
    
    f.close();   
    return;
}

//===================================================================================
//
// Added separate setup routines for hydrogenic levels (Aug 2015). These avoid loading
// data which makes startup slower + dependent on the maximum number of shells that
// were included by the tables.
//
//===================================================================================
double DnuH_func_Helium(int nup, int nlow)
{ return const_EH_inf_Hz/(1.0+const_me_malp)*(1.0/nlow/nlow-1.0/nup/nup); }

void check_Hydrogenic_transition(int n, int l, int s, int j, int np, int lp, int sp, int jp,
                                 double A21, vector<Transition_Data_HeI_A> &v)
{
    cout << n << " " << l << " " << s << " " << j << " --> " << np << " " << lp << " " << sp << " " << jp;

    bool was_found=0;
    for(int m=0; m<(int)v.size(); m++)
    if(v[m].np==np && v[m].lp==lp && v[m].jp==jp && v[m].sp==sp)
    {
        //cout << " found: " << v[m].A21 << " " << A21 << " " << v[m].A21/A21-1.0;
        if(fabs(v[m].A21/A21-1.0)>=0.005)
        { cout << " found: " << v[m].A21 << " " << A21 << " " << v[m].A21/A21-1.0; wait_f_r(); }
        was_found=1;
        break;
    }
    
    cout << endl;
    
    if(!was_found){ cout << " transition not found!" << endl; wait_f_r(); }
    
    return;
}

//===================================================================================
void compute_Hydrogenic_transition_inf_S(int n, int l,
                                         Atom_HeI_Singlet &Sing,
                                         vector<Transition_Data_HeI_A> &v)
{
    // simple sanity checks
    if(l>=n || l<0 || n<1) return;

    Transition_Data_HeI_A dd;
    
    for(int k=0; k<5; k++) dd.xxx[k]=0.0;
    dd.sp=0;
    
    double nun=Sing.Level(n, l).Get_nu_ion();

    for(int np=1; np<n; np++)
        for(int lp=l-1; lp<=l+1; lp+=2)
            if(lp>=0 && lp<np)
            {
                dd.np=np;
                dd.lp=lp;
                dd.jp=dd.lp;
                dd.gwp=2.0*dd.jp+1.0;
                dd.Dnu=Sing.Level(dd.np, dd.lp).Get_nu_ion()-nun;
                dd.DE=dd.Dnu*const_h/const_e;
                dd.lambda21=const_cl/dd.Dnu;
                
                double A21H=A_SH(1, const_malpha_mp, n, l, dd.np, dd.lp);
                double DnuH=DnuH_func_Helium(n, np);
                dd.A21=A21H*pow(dd.Dnu/DnuH, 3);
                v.push_back(dd);
                
                //check_Hydrogenic_transition(n, l, 0, l, np, lp, 0, lp, dd.A21, v);
            }
    
    //wait_f_r(" done transition singlet ");

    return;
}

//===================================================================================
// approximate transtion line ratios using Wigner 6J symbols
//===================================================================================
double Line_Ratio(int l, int j, int lp, int jp)
{
    int Dl=lp-l;
    int Dj=jp-j;
    double r;
    
    if(Dl==-1)
    {
        r=(2.0+l+j)*(1.0+l+j)*(-1.0+l+j)*(-2.0+l+j)/(4.0*j*(2.0*j+1.0)*l*(2.0*l-1.0));
        if(Dj>=0) r*=(2.0*j+1.0)*(-1.0-l+j)*(-2.0+l-j)/((1.0+j)*(1.0+l+j)*(-2.0+l+j));
        if(Dj==1) r*=(-3.0+l-j)*(j-l)*j/((2.0*j+1.0)*(-1.0+l+j)*(2.0+l+j));
    }
    else if(Dl==1)
    {
        r=(4.0+l+j)*(3.0+l+j)*(1.0+l+j)*(l+j)/(4.0*(j+1.0)*(2.0*j+1.0)*(l+1.0)*(2.0*l+3.0));
        if(Dj<=0) r*=(2.0*j+1.0)*(-1.0+l-j)*(-2.0-l+j)/(j*(1.0+l+j)*(4.0+l+j));
        if(Dj==-1) r*=(-3.0-l+j)*(l-j)*(1.0+j)/((2.0*j+1.0)*(l+j)*(3.0+l+j));
    }
    else{ cerr << " Line_Ratio::This should not happen..." << endl; exit(0); }
    
//    cout << l << " " << j << " || " << lp << " " << jp << endl;
//    cout << Dl << " " << Dj << " " << r << " " << (2.0*jp+1.0)/(2.0*lp+1.0)/3.0 << endl;
//    wait_f_r();
    
    return r;
}

//===================================================================================
void compute_Hydrogenic_transition_inf_T(int n, int l, int j, int np, int lp,
                                         Atom_HeI_Triplet &Trip,
                                         vector<Transition_Data_HeI_A> &v)
{
    // simple sanity checks
    if(np>=n) return;
    if(l==0 && j!=1) return;
    if(l>=n || l<0 || n<2 || j<l-1 || j>l+1) return;
    if(lp>=np || lp<0 || np<2) return;
    
    Transition_Data_HeI_A dd;
    
    for(int k=0; k<5; k++) dd.xxx[k]=0.0;
    dd.sp=1;
    dd.np=np;
    dd.lp=lp;
    
    double nun=Trip.Level(n, l, j).Get_nu_ion();
    
    for(int jp=lp-1; jp<=lp+1; jp++)
    {
        if(lp==0) jp=0+1;
        
        if(lp>=0 && lp<np && (j-jp==-1 || j-jp==0 || j-jp==1)) // additional check of Dj!
        {
            dd.jp=jp;
            dd.gwp=2.0*dd.jp+1.0;
            dd.Dnu=Trip.Level(dd.np, dd.lp, dd.jp).Get_nu_ion()-nun;
            dd.DE=dd.Dnu*const_h/const_e;
            dd.lambda21=const_cl/dd.Dnu;
            
            double A21H=A_SH(1, const_malpha_mp, n, l, dd.np, dd.lp);
            double DnuH=DnuH_func_Helium(n, np);
            // compute triplet transion with Wigner-expression
            if(switch_Hydrogenic) dd.A21=A21H*pow(dd.Dnu/DnuH, 3)*Line_Ratio(l, j, lp, jp);
            else dd.A21=A21H*pow(dd.Dnu/DnuH, 3)*(2.0*jp+1.0)/(2.0*lp+1.0)/3.0;
            v.push_back(dd);
        }
    }
    
    return;
}

//===================================================================================
void compute_Hydrogenic_transition_inf_T(int n, int l, int j,
                                         Atom_HeI_Triplet &Trip,
                                         vector<Transition_Data_HeI_A> &v)
{
    // simple sanity checks
    if(l==0 && j!=1) return;
    if(l>=n || l<0 || n<2 || j<l-1 || j>l+1) return;

    Transition_Data_HeI_A dd;
    
    for(int k=0; k<5; k++) dd.xxx[k]=0.0;
    dd.sp=1;
    
    double nun=Trip.Level(n, l, j).Get_nu_ion();
    
    for(int np=2; np<n; np++) // no np=1 for triplet case
        for(int lp=l-1; lp<=l+1; lp+=2)
            for(int jp=lp-1; jp<=lp+1; jp++)
            {
                if(lp==0) jp=0+1;
                
                if(lp>=0 && lp<np && (j-jp==-1 || j-jp==0 || j-jp==1)) // additional check of Dj!
                {
                    dd.np=np;
                    dd.lp=lp;
                    dd.jp=jp;
                    dd.gwp=2.0*dd.jp+1.0;
                    dd.Dnu=Trip.Level(dd.np, dd.lp, dd.jp).Get_nu_ion()-nun;
                    dd.DE=dd.Dnu*const_h/const_e;
                    dd.lambda21=const_cl/dd.Dnu;
                    
                    double A21H=A_SH(1, const_malpha_mp, n, l, dd.np, dd.lp);
                    double DnuH=DnuH_func_Helium(n, np);
                    if(switch_Hydrogenic) dd.A21=A21H*pow(dd.Dnu/DnuH, 3)*Line_Ratio(l, j, lp, jp);
                    else dd.A21=A21H*pow(dd.Dnu/DnuH, 3)*(2.0*jp+1.0)/(2.0*lp+1.0)/3.0;
                    v.push_back(dd);
                    
                    //check_Hydrogenic_transition(n, l, 1, j, np, lp, 1, jp, dd.A21, v);
                }
            }

    //wait_f_r(" done transition triplet ");
    
    return;
}

//===================================================================================
void compute_Hydrogenic_transition_inf_T_no_j(int n, int l,
                                              Atom_HeI_Triplet_no_j &Trip_no_j,
                                              vector<Transition_Data_HeI_A> &v)
{
    // simple sanity checks
    if(l>=n || l<0 || n<=Trip_no_j.Get_njres()) return;
    
    Transition_Data_HeI_A dd;
    
    for(int k=0; k<5; k++) dd.xxx[k]=0.0;
    dd.sp=1;
    
    double nun=Trip_no_j.Level(n, l).Get_nu_ion();
    
    // non-j-resolved --> j-resolved
    if(Trip_glob!=NULL)
    {
        for(int np=2; np<=Trip_no_j.Get_njres(); np++)
            for(int lp=l-1; lp<=l+1; lp+=2)
                for(int jp=lp-1; jp<=lp+1; jp++)
                {
                    if(lp==0) jp=0+1;
                    
                    if(lp>=0 && lp<np)
                    {
                        dd.np=np;
                        dd.lp=lp;
                        dd.jp=jp;
                        dd.gwp=2.0*dd.jp+1.0;
                        dd.Dnu=Trip_glob->Level(dd.np, dd.lp, dd.jp).Get_nu_ion()-nun;
                        dd.DE=dd.Dnu*const_h/const_e;
                        dd.lambda21=const_cl/dd.Dnu;
                        
                        double A21H=A_SH(1, const_malpha_mp, n, l, dd.np, dd.lp);
                        double DnuH=DnuH_func_Helium(n, np);
                        dd.A21=A21H*pow(dd.Dnu/DnuH, 3)*(2.0*jp+1.0)/(2.0*lp+1.0)/3.0;
                        v.push_back(dd);
                        
                        //check_Hydrogenic_transition(n, l, 1, -10, np, lp, 1, jp, dd.A21, v);
                    }
                }
    }
    
    // non-j-resolved --> non-j-resolved
    for(int np=Trip_no_j.Get_njres()+1; np<n; np++)
        for(int lp=l-1; lp<=l+1; lp+=2)
            if(lp>=0 && lp<np)
            {
                dd.np=np;
                dd.lp=lp;
                dd.jp=-10;
                dd.gwp=3.0*(2*dd.lp+1);
                dd.Dnu=Trip_no_j.Level(dd.np, dd.lp).Get_nu_ion()-nun;
                dd.DE=dd.Dnu*const_h/const_e;
                dd.lambda21=const_cl/dd.Dnu;
                
                double A21H=A_SH(1, const_malpha_mp, n, l, dd.np, dd.lp);
                double DnuH=DnuH_func_Helium(n, np);
                dd.A21=A21H*pow(dd.Dnu/DnuH, 3);
                v.push_back(dd);
                
                //check_Hydrogenic_transition(n, l, 1, -10, np, lp, 1, -10, dd.A21, v);
            }

    //wait_f_r(" done transition non-j ");
    
    return;
}

//===================================================================================
//===================================================================================


//###################################################################################
//
// class: Electron_Level_HeI_Singlet
//
//###################################################################################

//===================================================================================
// Konstructors and Destructors for Electron_Level_HeI_Singlet
//===================================================================================
void Electron_Level_HeI_Singlet::init(int n, int l, int mflag)
{
    nn=n; ll=l; 
    mess_flag=mflag;
    gw=g_l();
    A_values.clear();
    
    //===============================================================================
    // energies
    //===============================================================================
    // quantum defect values
    if(ll<=6 && nn>10 && nn<=30) DE=He1s2_ion-compute_DEc_QD(nn, ll, 0, ll);
    // hydrogenic energies for all levels l>=7 & n>=11
    else if((ll>6 && nn>10) || (ll<=6 && nn>30)) DE=He1s2_ion-const_EHeI/nn/nn;
    else if((nn==9 && ll==8) || (nn==10 && ll==7) || (nn==10 && ll==8) || (nn==10 && ll==9)) 
        DE=He1s2_ion-const_EHeI/nn/nn;
    // the rest is in Drake & Morton
    else read_level_energy_DM(nameA_SS, nn, ll, 0, ll, DE);
    //===============================================================================
    
    if(DE==-6000.0){ Dnu=Eion=Eion_ergs=nuion=0.0; cout << " BAAAADDD " << endl; }
    else 
    {
        Dnu=DE*const_e/const_h;
        Eion=He1s2_ion-DE;
        Eion_ergs=Eion*const_e;
        nuion=Eion_ergs/const_h;
        
        //===========================================================================
        // all transitions from that level
        //===========================================================================
        if(HeI_Atom_read_transition_data==1)
        {
            // set Drake and Morton values
            read_transition_inf(nameA_SS, nn, ll, 0, ll, A_values);
            read_transition_inf(nameA_TS, nn, ll, 0, ll, A_values);
            read_transition_inf(nameA_add, nn, ll, 0, ll, A_values);

            if(_HeI_add_Quad>=nn && ll==2)
                read_transition_add_Quadrupole(nameA_add_Quad, nn, ll, 0, ll, A_values, mflag);
        }
    }
    
    create_Transition_lookup_table();
    
    if(mess_flag>0)
    {
        cout << "\n %##############################################################%" << endl;
        
        display_level_information();
        cout << endl;
        display_all_downward_transitions();
        
        cout << " %##############################################################%" << endl;
    }
    
    return;
}

//===================================================================================
Electron_Level_HeI_Singlet::Electron_Level_HeI_Singlet(int n, int l, int mflag)
{ init(n, l, mflag); }

//===================================================================================
Electron_Level_HeI_Singlet::~Electron_Level_HeI_Singlet()
{ A_values.clear(); }

//===================================================================================
void Electron_Level_HeI_Singlet::Set_hydrogenic_transitions(Atom_HeI_Singlet &Sing)
{
    compute_Hydrogenic_transition_inf_S(nn, ll, Sing, A_values);
    create_Transition_lookup_table();
    return;
}

//===================================================================================
void Electron_Level_HeI_Singlet::create_Transition_lookup_table()
{
    ZERO_Data.np=ZERO_Data.sp=ZERO_Data.lp=ZERO_Data.jp=ZERO_Data.gwp=0;
    ZERO_Data.A21=ZERO_Data.lambda21=ZERO_Data.Dnu=ZERO_Data.DE=0.0;
    
    Transition_Data_HeI_A_lookup_table.clear();
    vector<int> dum1;
    vector<int> dum2;
    vector<vector<int> > dum22;

    //=============================================================
    // get all transitions to s=0
    //=============================================================
    dum1.clear();
    for(int m=0; m<Get_n_down(); m++)
        if(A_values[m].sp==0) dum1.push_back(m);
    
    //=============================================================
    // sort according to np
    //=============================================================
    int nmin=100000, nmax=0;
    for(int m=0; m<(int)dum1.size(); m++)
    {   
        if(A_values[dum1[m]].np<nmin) nmin=A_values[dum1[m]].np;
        if(A_values[dum1[m]].np>nmax) nmax=A_values[dum1[m]].np;
    }
    if(dum1.size()==0) nmax=nmin=0;

    dum2.clear();
    dum2.push_back(-1); // n=0..nmin transtion
    dum22.clear();
    for(int np=0; np<nmin; np++) dum22.push_back(dum2);

    if(nmin>=1) for(int np=nmin; np<=nmax; np++)
    {
        dum2.clear();
        
        for(int m=0; m<(int)dum1.size(); m++)
            if(A_values[dum1[m]].np==np) dum2.push_back(dum1[m]);
        
        if(dum2.size()==0) dum2.push_back(-1); 
        
        dum22.push_back(dum2);
    }       
    Transition_Data_HeI_A_lookup_table.push_back(dum22);
        
    //=============================================================
    // get all transitions to s=1
    //=============================================================
    dum1.clear();
    for(int m=0; m<Get_n_down(); m++)
        if(A_values[m].sp==1) dum1.push_back(m);
    
    //=============================================================
    // sort according to np
    //=============================================================
    nmin=100000; nmax=0;
    for(int m=0; m<(int)dum1.size(); m++)
    {   
        if(A_values[dum1[m]].np<nmin) nmin=A_values[dum1[m]].np;
        if(A_values[dum1[m]].np>nmax) nmax=A_values[dum1[m]].np;
    }
    if(dum1.size()==0) nmax=nmin=0;
    
    dum2.clear();
    dum2.push_back(-1); // n=0..nmin transtion
    dum22.clear();
    for(int np=0; np<nmin; np++) dum22.push_back(dum2);
    
    if(nmin>=1) for(int np=nmin; np<=nmax; np++)
    {
        dum2.clear();
        
        for(int m=0; m<(int)dum1.size(); m++)
            if(A_values[dum1[m]].np==np) dum2.push_back(dum1[m]);
        
        if(dum2.size()==0) dum2.push_back(-1); 
        
        dum22.push_back(dum2);
    }       
    Transition_Data_HeI_A_lookup_table.push_back(dum22);
    
    // check the data
    if(check_Transition_Data_Tables)
    {
        Transition_Data_HeI_A T;
        for(int m=0; m<Get_n_down(); m++)
        {
            T=Get_Trans_Data(A_values[m].np, A_values[m].lp, A_values[m].sp, A_values[m].jp);
//          cout << m << " " << T.A21 << " " << A_values[m].A21 << endl;
            if(T.np==0 || T.A21-A_values[m].A21!=0.0) wait_f_r("create_Transition_lookup_table::NOOOO");        
            if(T.np>nn){ cout << "create_Transition_lookup_table::Singlet:"
                              << " there is some transtions with np>nn. Here np= " 
                              << T.np << " and nn= " << nn << endl; wait_f_r(); }
        }
    }
    
    return;
}

//===================================================================================
void Electron_Level_HeI_Singlet::display_Trans_Data(const Transition_Data_HeI_A &v)
{
  cout.precision(16);

  cout << " %==============================================================%" << endl;
  cout << " % Electron_Level_HeI_Singlet::display_Trans_Data:\n % (" 
       << nn << ", " << ll << ", 0, " << ll << ") --> (" 
       << v.np << ", " << v.lp << ", " << v.sp << ", " << v.jp << "), gw:" << gw << endl;
  cout << " %==============================================================%" << endl;
  cout << scientific << " Dnu: " << v.Dnu*1.0e-9 << " GHz, DE: " << v.DE << " eV, lambda21: " 
       << v.lambda21 << " cm"<< endl;
  //cout << " A21: " << v.A21  << " sec^-1, f: " << v.f << endl << endl; 
  cout << scientific << " A21: " << v.A21  << " sec^-1, gw: " << v.gwp << endl << endl; 
  
  return;
}

void Electron_Level_HeI_Singlet::display_level_information()
{
  cout.precision(16);

  cout << " %==============================================================%" << endl;
  cout << " % Electron_Level_HeI_Singlet::display_level_information:\n % (n, l, s, j) == (" 
       << nn << ", " << ll << ", 0, " << ll << "), gw:" << gw  
       << endl;
  cout << " %==============================================================%" << endl;
  
  cout << "\n gw: " << gw << endl;
  cout << scientific << " Dnu_1s: " << Dnu*1.0e-9 << " GHz, nuion: " << nuion*1.0e-9 << " GHz" << endl;
  cout << scientific << " DE_1s: " << DE << " eV, Eion: " << Eion << " eV " << endl;
  cout << " Number of transitions: " << Get_n_down() << endl << endl;
  cout << " %==============================================================%" << endl;

  return;
}

//===================================================================================
void Electron_Level_HeI_Singlet::display_downward_transition(int i)
{
    if(i>=(int)A_values.size())
    { 
    cout << " Electron_Level::display_downward_transition::no transition tabulated " << i << endl;
    return;
    }
    
    display_Trans_Data(A_values[i]);

    return;
}

//===================================================================================
void Electron_Level_HeI_Singlet::display_all_downward_transitions()
{
    for(int l=0; l<(int)A_values.size(); l++) display_Trans_Data(A_values[l]);

    return;
}

//======================================================================
const Transition_Data_HeI_A& Electron_Level_HeI_Singlet::Get_Trans_Data(int np, int lp, 
                                                                        int sp, int jp) const
{
    if(sp==0 && (int)Transition_Data_HeI_A_lookup_table[0].size()>np)
    {
        int ntrans=Transition_Data_HeI_A_lookup_table[0][np].size();
        
        for(int m=0; m<ntrans; m++) 
        {
            if(Transition_Data_HeI_A_lookup_table[0][np][m]!=-1)
            {
                if(np==A_values[Transition_Data_HeI_A_lookup_table[0][np][m]].np && 
                   sp==A_values[Transition_Data_HeI_A_lookup_table[0][np][m]].sp && 
                   lp==A_values[Transition_Data_HeI_A_lookup_table[0][np][m]].lp && 
                   jp==A_values[Transition_Data_HeI_A_lookup_table[0][np][m]].jp)
                    return A_values[Transition_Data_HeI_A_lookup_table[0][np][m]];      
            }
        }
    }
    if(sp==1 && (int)Transition_Data_HeI_A_lookup_table[1].size()>np)
    {
        int ntrans=Transition_Data_HeI_A_lookup_table[1][np].size();
        
        for(int m=0; m<ntrans; m++) 
        {
            if(Transition_Data_HeI_A_lookup_table[1][np][m]!=-1)
            {
                if(np==A_values[Transition_Data_HeI_A_lookup_table[1][np][m]].np && 
                   sp==A_values[Transition_Data_HeI_A_lookup_table[1][np][m]].sp && 
                   lp==A_values[Transition_Data_HeI_A_lookup_table[1][np][m]].lp && 
                   jp==A_values[Transition_Data_HeI_A_lookup_table[1][np][m]].jp)
                    return A_values[Transition_Data_HeI_A_lookup_table[1][np][m]];      
            }
        }
    }
    
    return ZERO_Data;
}

//===================================================================================
double Electron_Level_HeI_Singlet::Get_A21(int np, int lp, int sp, int jp) const
{
  for(int i=0; i<(int)A_values.size(); i++) 
    {
      if(A_values[i].np==np && A_values[i].lp==lp && A_values[i].sp==sp && A_values[i].jp==jp) 
          return A_values[i].A21;
      //else return 0.0;
    }

  return 0.0;
}

//===================================================================================
double Electron_Level_HeI_Singlet::Get_nu21(int np, int lp, int sp, int jp) const
{
    for(int i=0; i<(int)A_values.size(); i++) 
    {
        if(A_values[i].np==np && A_values[i].lp==lp && A_values[i].sp==sp && A_values[i].jp==jp) 
            return A_values[i].Dnu;
        //else return 0.0;
    }
    
    return 0.0;
}

//===================================================================================
double Electron_Level_HeI_Singlet::Get_lambda21(int np, int lp, int sp, int jp) const
{
    for(int i=0; i<(int)A_values.size(); i++) 
    {
//      cout << i << " " 
//           << A_values[i].np << " " << A_values[i].lp << " " 
//           << A_values[i].sp << " " << A_values[i].jp << endl;
        
        if(A_values[i].np==np && A_values[i].lp==lp && A_values[i].sp==sp && A_values[i].jp==jp) 
            return A_values[i].lambda21;
        //else return 0.0;
    }
    
    return 0.0;
}

//===================================================================================
void Electron_Level_HeI_Singlet::Set_xxx_Transition_Data(int m, int d, double v)
{ 
  if(d<5)
    {
      A_values[m].xxx[d]=v; 
    }
  else cout << " Electron_Level_HeI_Singlet::Set_xxx_Transition_Data: check index " << endl;

  return; 
}


//==============================================================================
// Saha-relations with continuum
//==============================================================================
// f(T) == (Ni/[Ne Nc])_LTE
double Electron_Level_HeI_Singlet::Ni_NeNc_LTE(double TM) const
{
    // HeI --> gc=2 --> ge*gc=4
    // 30.07.2015: added reduced mass factor
    return gw/4.0*pow(const_lambdac, 3)*pow(2.0*PI*const_kb_mec2*Get_mu_red()*TM, -1.5)
                 *exp(Eion_ergs/const_kB/TM );
}

double Electron_Level_HeI_Singlet::Xi_Saha(double Xe, double Xc, double NH, double TM) const
{ return Xe*Xc*NH*Ni_NeNc_LTE(TM); }

double Electron_Level_HeI_Singlet::Ni_Saha(double Ne, double Nc, double TM) const
{ return Ne*Nc*Ni_NeNc_LTE(TM); }


//###################################################################################
//
// class: Electron_Level_HeI_Triplet
//
//###################################################################################

//===================================================================================
// Konstructors and Destructors for Electron_Level_HeI_Triplet
//===================================================================================
void Electron_Level_HeI_Triplet::init(int n, int l, int j, int mflag)
{
    nn=n; ll=l; jj=j; 
    mess_flag=mflag;
    gw=g_l();
    A_values.clear();
    
    //===============================================================================
    // energies
    //===============================================================================
    // quantum defect values
    if(ll<=6 && nn>10 && nn<=30) DE=He1s2_ion-compute_DEc_QD(nn, ll, 1, jj);
    // hydrogenic energies
    else if((ll>6 && nn>10) || (ll<=6 && nn>30)) DE=He1s2_ion-const_EHeI/nn/nn;
    else if((nn==8 && ll==7 && jj==6)) read_level_energy_DM(nameA_TS, nn, ll, 1, jj, DE);
    else if((nn==8 && ll==7) || (nn==9 && ll==8) || 
            (nn==10 && ll==7)|| (nn==10 && ll==8)|| (nn==10 && ll==9)) DE=He1s2_ion-const_EHeI/nn/nn;
    //
    // the rest is in Drake & Morton
    else read_level_energy_DM(nameA_TT, nn, ll, 1, jj, DE);
    //===============================================================================
    
    if(DE==-6000.0){ Dnu=Eion=Eion_ergs=nuion=0.0; cout << " BAAAADDD " << endl; }
    else
    {
        Dnu=DE*const_e/const_h;
        Eion=He1s2_ion-DE;
        Eion_ergs=Eion*const_e;
        nuion=Eion_ergs/const_h;
        
        //===========================================================================
        // all transitions from that level
        //===========================================================================
        if(HeI_Atom_read_transition_data==1)
        {
            // set Drake and Morton values (j-resolved)
            read_transition_inf(nameA_TT, nn, ll, 1, jj, A_values);
            read_transition_inf(nameA_ST, nn, ll, 1, jj, A_values);
            read_transition_inf(nameA_add, nn, ll, 1, jj, A_values);
            if(_HeI_add_Intercomb>=nn && nn>2 && ll==1 && jj==1)
                read_transition_add_TS(nameA_add_Int, nn, ll, 1, jj, A_values, mflag); 
        }
    }
    
    create_Transition_lookup_table();

    if(mess_flag>0)
    {
        cout << "\n %##############################################################%" << endl;
        
        display_level_information();
        cout << endl;
        display_all_downward_transitions();
        
        cout << " %##############################################################%" << endl;
    }
    
    return;
}

//===================================================================================
Electron_Level_HeI_Triplet::Electron_Level_HeI_Triplet(int n, int l, int j, int mflag)
{ init(n, l, j, mflag); }

Electron_Level_HeI_Triplet::~Electron_Level_HeI_Triplet()
{ A_values.clear(); }

//===================================================================================
void Electron_Level_HeI_Triplet::Set_hydrogenic_transitions(Atom_HeI_Triplet &Trip)
{
    compute_Hydrogenic_transition_inf_T(nn, ll, jj, Trip, A_values);
    create_Transition_lookup_table();
    
    return;
}

void Electron_Level_HeI_Triplet::Set_hydrogenic_transitions(int np, int jp, Atom_HeI_Triplet &Trip)
{
    compute_Hydrogenic_transition_inf_T(nn, ll, jj, np, jp, Trip, A_values);
    create_Transition_lookup_table();
    
    return;
}

//===================================================================================
void Electron_Level_HeI_Triplet::create_Transition_lookup_table()
{
    ZERO_Data.np=ZERO_Data.sp=ZERO_Data.lp=ZERO_Data.jp=ZERO_Data.gwp=0;
    ZERO_Data.A21=ZERO_Data.lambda21=ZERO_Data.Dnu=ZERO_Data.DE=0.0;
    
    Transition_Data_HeI_A_lookup_table.clear();
    vector<int> dum1;
    vector<int> dum2;
    vector<vector<int> > dum22;
    
    //=============================================================
    // get all transitions to s=0
    //=============================================================
    dum1.clear();
    for(int m=0; m<Get_n_down(); m++)
        if(A_values[m].sp==0) dum1.push_back(m);
    
    //=============================================================
    // sort according to np
    //=============================================================
    int nmin=100000, nmax=0;
    for(int m=0; m<(int)dum1.size(); m++)
    {   
        if(A_values[dum1[m]].np<nmin) nmin=A_values[dum1[m]].np;
        if(A_values[dum1[m]].np>nmax) nmax=A_values[dum1[m]].np;
    }
    if(dum1.size()==0) nmax=nmin=0;
    
    dum2.clear();
    dum2.push_back(-1); // n=0..nmin transtion
    dum22.clear();
    for(int np=0; np<nmin; np++) dum22.push_back(dum2);
    
    if(nmin>=1) for(int np=nmin; np<=nmax; np++)
    {
        dum2.clear();
        
        for(int m=0; m<(int)dum1.size(); m++)
            if(A_values[dum1[m]].np==np) dum2.push_back(dum1[m]);
        
        if(dum2.size()==0) dum2.push_back(-1); 
        
        dum22.push_back(dum2);
    }       
    Transition_Data_HeI_A_lookup_table.push_back(dum22);
    
    //=============================================================
    // get all transitions to s=1
    //=============================================================
    dum1.clear();
    for(int m=0; m<Get_n_down(); m++)
        if(A_values[m].sp==1) dum1.push_back(m);
    
    //=============================================================
    // sort according to np
    //=============================================================
    nmin=100000; nmax=0;
    for(int m=0; m<(int)dum1.size(); m++)
    {   
        if(A_values[dum1[m]].np<nmin) nmin=A_values[dum1[m]].np;
        if(A_values[dum1[m]].np>nmax) nmax=A_values[dum1[m]].np;
    }
    if(dum1.size()==0) nmax=nmin=0;
    
    dum2.clear();
    dum2.push_back(-1); // n=0..nmin transtion
    dum22.clear();
    for(int np=0; np<nmin; np++) dum22.push_back(dum2);
    
    if(nmin>=1) for(int np=nmin; np<=nmax; np++)
    {
        dum2.clear();
        
        for(int m=0; m<(int)dum1.size(); m++)
            if(A_values[dum1[m]].np==np) dum2.push_back(dum1[m]);
        
        if(dum2.size()==0) dum2.push_back(-1); 
        
        dum22.push_back(dum2);
    }       
    Transition_Data_HeI_A_lookup_table.push_back(dum22);
    
    // check the data
    if(check_Transition_Data_Tables)
    {
        Transition_Data_HeI_A T;
        for(int m=0; m<Get_n_down(); m++)
        {
            T=Get_Trans_Data(A_values[m].np, A_values[m].lp, A_values[m].sp, A_values[m].jp);
//          cout << m << " " << T.A21 << " " << A_values[m].A21 << " " << T.np << " " 
//               << T.lp << " " << T.sp << " " << T.jp << " <-- " 
//               << nn << " " << ll << " " << 1 << " " << jj << endl;
            
            if(T.np==0 || T.A21-A_values[m].A21!=0.0) 
                wait_f_r(" create_Transition_lookup_table::NOOOO");
            if(T.np>nn)
            { 
                cout << " create_Transition_lookup_table::Triplet:"
                     << " there is some transtions with np>nn. Here np= " << T.np 
                     << " and nn= " << nn << endl; 
                wait_f_r(); 
            }
        }
    }
    
    return;
}

//===================================================================================
void Electron_Level_HeI_Triplet::display_Trans_Data(const Transition_Data_HeI_A &v)
{
  cout.precision(16);

  cout << " %==============================================================%" << endl;
  cout << " % Electron_Level_HeI_Triplet::display_Trans_Data:\n % (" 
       << nn << ", " << ll << ", 1, " << jj << ") --> (" 
       << v.np << ", " << v.lp << ", " << v.sp << ", " << v.jp << "), gw:" << gw << endl;
  cout << " %==============================================================%" << endl;
  cout << scientific << " Dnu: " << v.Dnu*1.0e-9 << " GHz, DE: " << v.DE << " eV, lambda21: " 
       << v.lambda21 << " cm"<< endl;
  //cout << " A21: " << v.A21  << " sec^-1, f: " << v.f << endl << endl; 
  cout << scientific << " A21: " << v.A21  << " sec^-1, gw: " << v.gwp << endl << endl; 
  
  return;
}

//===================================================================================
void Electron_Level_HeI_Triplet::display_level_information()
{
  cout.precision(16);

  cout << " %==============================================================%" << endl;
  cout << " % Electron_Level_HeI_Triplet::display_level_information:\n % (n, l, s, j) == (" 
       << nn << ", " << ll << ", 1, " << jj << "), gw:" << gw  
       << endl;
  cout << " %==============================================================%" << endl;
  
  cout << "\n gw: " << gw << endl;
  cout << scientific << " Dnu_1s: " << Dnu*1.0e-9 << " GHz, nuion: " 
       << nuion*1.0e-9 << " GHz" << endl;
  cout << scientific << " DE_1s: " << DE << " eV, Eion: " << Eion << " eV " << endl;
  cout << " Number of transitions: " << Get_n_down() << endl << endl;
  cout << " %==============================================================%" << endl;

  return;
}

void Electron_Level_HeI_Triplet::display_downward_transition(int i)
{
    if(i>=(int)A_values.size())
    { 
        cout << " Electron_Level_HeI_Triplet::display_downward_transition::"
             << " no transition tabulated " << i << endl;
        return;
    }
    
    display_Trans_Data(A_values[i]);

    return;
}

void Electron_Level_HeI_Triplet::display_all_downward_transitions()
{
    for(int l=0; l<(int)A_values.size(); l++) display_Trans_Data(A_values[l]);

    return;
}

const Transition_Data_HeI_A& Electron_Level_HeI_Triplet::Get_Trans_Data(int np, int lp, 
                                                                        int sp, int jp) const
{
    if(sp==0 && (int)Transition_Data_HeI_A_lookup_table[0].size()>np)
    {
        int ntrans=Transition_Data_HeI_A_lookup_table[0][np].size();
        
        for(int m=0; m<ntrans; m++) 
        {           
            if(Transition_Data_HeI_A_lookup_table[0][np][m]!=-1)
            {
                if(np==A_values[Transition_Data_HeI_A_lookup_table[0][np][m]].np && 
                   sp==A_values[Transition_Data_HeI_A_lookup_table[0][np][m]].sp && 
                   lp==A_values[Transition_Data_HeI_A_lookup_table[0][np][m]].lp && 
                   jp==A_values[Transition_Data_HeI_A_lookup_table[0][np][m]].jp)
                    return A_values[Transition_Data_HeI_A_lookup_table[0][np][m]];      
            }
        }
    }
    if(sp==1 && (int)Transition_Data_HeI_A_lookup_table[1].size()>np)
    {
        int ntrans=Transition_Data_HeI_A_lookup_table[1][np].size();
        
        for(int m=0; m<ntrans; m++) 
        {
            if(Transition_Data_HeI_A_lookup_table[1][np][m]!=-1)
            {
                if(np==A_values[Transition_Data_HeI_A_lookup_table[1][np][m]].np && 
                   sp==A_values[Transition_Data_HeI_A_lookup_table[1][np][m]].sp && 
                   lp==A_values[Transition_Data_HeI_A_lookup_table[1][np][m]].lp && 
                   jp==A_values[Transition_Data_HeI_A_lookup_table[1][np][m]].jp)
                    return A_values[Transition_Data_HeI_A_lookup_table[1][np][m]];      
            }
        }
    }
    
    return ZERO_Data;
}

//======================================================================
double Electron_Level_HeI_Triplet::Get_A21(int np, int lp, int sp, int jp) const
{
  for(int i=0; i<(int)A_values.size(); i++) 
    {
      if(A_values[i].np==np && A_values[i].lp==lp && A_values[i].sp==sp && A_values[i].jp==jp) 
          return A_values[i].A21;
      //else return 0.0;
    }

  return 0.0;
}

double Electron_Level_HeI_Triplet::Get_nu21(int np, int lp, int sp, int jp) const
{
    for(int i=0; i<(int)A_values.size(); i++) 
    {
        if(A_values[i].np==np && A_values[i].lp==lp && A_values[i].sp==sp && A_values[i].jp==jp) 
            return A_values[i].Dnu;
        //else return 0.0;
    }
    
    return 0.0;
}

double Electron_Level_HeI_Triplet::Get_lambda21(int np, int lp, int sp, int jp) const
{
    for(int i=0; i<(int)A_values.size(); i++) 
    {
        if(A_values[i].np==np && A_values[i].lp==lp && A_values[i].sp==sp && A_values[i].jp==jp) 
            return A_values[i].lambda21;
        //else return 0.0;
    }
    
    return 0.0;
}

void Electron_Level_HeI_Triplet::Set_xxx_Transition_Data(int m, int d, double v) 
{ 
  if(d<5)
    {
      A_values[m].xxx[d]=v; 
    }
  else cout << " Electron_Level_HeI_Triplet::Set_xxx_Transition_Data: check index " << endl;

  return; 
}

//==============================================================================
// Saha-relations with continuum
//==============================================================================
// f(T) == (Ni/[Ne Nc])_LTE
double Electron_Level_HeI_Triplet::Ni_NeNc_LTE(double TM) const
{ 
    // HeI --> gc=2 --> ge*gc=4
    // 30.07.2015: added reduced mass factor
    return gw/4.0*pow(const_lambdac, 3)*pow(2.0*PI*const_kb_mec2*Get_mu_red()*TM, -1.5)
                 *exp(Eion_ergs/const_kB/TM );
}

double Electron_Level_HeI_Triplet::Xi_Saha(double Xe, double Xc, double NH, double TM) const
{ return Xe*Xc*NH*Ni_NeNc_LTE(TM); }

double Electron_Level_HeI_Triplet::Ni_Saha(double Ne, double Nc, double TM) const
{ return Ne*Nc*Ni_NeNc_LTE(TM); }



//###################################################################################
//
// class: Electron_Level_HeI_Triplet_no_j
//
//###################################################################################

//===================================================================================
// Konstructors and Destructors for Electron_Level_HeI_Triplet_no_j
//===================================================================================
void Electron_Level_HeI_Triplet_no_j::init(int n, int l, int njres, int mflag)
{
    nn=n; ll=l;
    mess_flag=mflag;
    gw=g_l();
    A_values.clear();
    
    //===============================================================================
    // above njresolved the triplet states will not be treated j-resolved.
    // this should not be smaller than 10
    //===============================================================================
    njresolved=(int)min(max(njres, 10), n);

    //===============================================================================
    // energies
    //===============================================================================
    // quantum defect values
    if(ll<=6 && nn>10 && nn<=30)
    { 
        double dum;
        dum=He1s2_ion-compute_DEc_QD(nn, ll, 1, ll+1);
        DE=dum*(2*(ll+1)+1)/(2*ll+1);
        if(l>0)
        {
            dum=He1s2_ion-compute_DEc_QD(nn, ll, 1, ll);
            DE+=dum;
            dum=He1s2_ion-compute_DEc_QD(nn, ll, 1, ll-1);
            DE+=dum*(2*(ll-1)+1)/(2*ll+1);
        }
        DE/=3.0;
    }
    
    // hydrogenic energies
    else if((ll>6 && nn>10) || (ll<=6 && nn>30)) DE=He1s2_ion-const_EHeI/nn/nn;
    else{ cout << " n<=10 is not allowed !" << nn << endl; exit(0); }
    //===============================================================================
    
    if(DE==-6000.0){ Dnu=Eion=Eion_ergs=nuion=0.0; cout << " BAAAADDD " << endl; }
    else
    {
        Dnu=DE*const_e/const_h;
        Eion=He1s2_ion-DE;
        Eion_ergs=Eion*const_e;
        nuion=Eion_ergs/const_h;
    }
    
    create_Transition_lookup_table();

    if(mess_flag>0)
    {
        cout << "\n %##############################################################%" << endl;
        
        display_level_information();
        cout << endl;
        display_all_downward_transitions();
        
        cout << " %##############################################################%" << endl;
    }
    
    return;
}

Electron_Level_HeI_Triplet_no_j::Electron_Level_HeI_Triplet_no_j(int n, int l, int njres, int mflag)
{ init(n, l, njres, mflag); }

Electron_Level_HeI_Triplet_no_j::~Electron_Level_HeI_Triplet_no_j()
{ A_values.clear(); }

//===================================================================================
void Electron_Level_HeI_Triplet_no_j::Set_hydrogenic_transitions(Atom_HeI_Triplet_no_j &Trip_no_j)
{
    compute_Hydrogenic_transition_inf_T_no_j(nn, ll, Trip_no_j, A_values);
    create_Transition_lookup_table();
    
    return;
}

//===================================================================================
void Electron_Level_HeI_Triplet_no_j::create_Transition_lookup_table()
{
    ZERO_Data.np=ZERO_Data.sp=ZERO_Data.lp=ZERO_Data.jp=ZERO_Data.gwp=0;
    ZERO_Data.A21=ZERO_Data.lambda21=ZERO_Data.Dnu=ZERO_Data.DE=0.0;
    
    Transition_Data_HeI_A_lookup_table.clear();
    vector<int> dum1;
    vector<int> dum2;
    vector<vector<int> > dum22;
    
    //=============================================================
    // get all transitions to s=0
    //=============================================================
    dum1.clear();
    for(int m=0; m<Get_n_down(); m++)
        if(A_values[m].sp==0) dum1.push_back(m);
    
    //=============================================================
    // sort according to np
    //=============================================================
    int nmin=100000, nmax=0;
    for(int m=0; m<(int)dum1.size(); m++)
    {   
        if(A_values[dum1[m]].np<nmin) nmin=A_values[dum1[m]].np;
        if(A_values[dum1[m]].np>nmax) nmax=A_values[dum1[m]].np;
    }
    if(dum1.size()==0) nmax=nmin=0;
    
    dum2.clear();
    dum2.push_back(-1); // n=0..nmin transtion
    dum22.clear();
    for(int np=0; np<nmin; np++) dum22.push_back(dum2);
    
    if(nmin>=1) for(int np=nmin; np<=nmax; np++)
    {
        dum2.clear();
        
        for(int m=0; m<(int)dum1.size(); m++)
            if(A_values[dum1[m]].np==np) dum2.push_back(dum1[m]);
        
        if(dum2.size()==0) dum2.push_back(-1); 
        
        dum22.push_back(dum2);
    }       
    Transition_Data_HeI_A_lookup_table.push_back(dum22);
    
    //=============================================================
    // get all transitions to s=1
    //=============================================================
    dum1.clear();
    for(int m=0; m<Get_n_down(); m++)
        if(A_values[m].sp==1) dum1.push_back(m);
    
    //=============================================================
    // sort according to np
    //=============================================================
    nmin=100000; nmax=0;
    for(int m=0; m<(int)dum1.size(); m++)
    {   
        if(A_values[dum1[m]].np<nmin) nmin=A_values[dum1[m]].np;
        if(A_values[dum1[m]].np>nmax) nmax=A_values[dum1[m]].np;
    }
    if(dum1.size()==0) nmax=nmin=0;
    
    dum2.clear();
    dum2.push_back(-1); // n=0..nmin transtion
    dum22.clear();
    for(int np=0; np<nmin; np++) dum22.push_back(dum2);
    
    if(nmin>=1) for(int np=nmin; np<=nmax; np++)
    {
        dum2.clear();
        
        for(int m=0; m<(int)dum1.size(); m++)
            if(A_values[dum1[m]].np==np) dum2.push_back(dum1[m]);
        
        if(dum2.size()==0) dum2.push_back(-1); 
        
        dum22.push_back(dum2);
    }       
    Transition_Data_HeI_A_lookup_table.push_back(dum22);
    
    // check the data
    if(check_Transition_Data_Tables)
    {
        Transition_Data_HeI_A T;
        for(int m=0; m<Get_n_down(); m++)
        {
            T=Get_Trans_Data(A_values[m].np, A_values[m].lp, A_values[m].sp, A_values[m].jp);
//          cout << m << " " << T.A21 << " " << A_values[m].A21 << endl;

            if(T.np==0 || T.A21-A_values[m].A21!=0.0) 
                wait_f_r(" create_Transition_lookup_table::NOOOO");
            if(T.np>nn)
            { 
                cout << " create_Transition_lookup_table::Triplet_no_j:"
                     << " there is some transtions with np>nn. Here np= " << T.np 
                     << " and nn= " << nn << endl;
                wait_f_r(); 
            }
        }
    }
    
    return;
}

void Electron_Level_HeI_Triplet_no_j::display_Trans_Data(const Transition_Data_HeI_A &v)
{
  cout.precision(16);

  cout << " %==============================================================%" << endl;
  cout << " % Electron_Level_HeI_Triplet_no_j::display_Trans_Data:\n % (" 
       << nn << ", " << ll << ", 1, -10) --> (";
  // is the final Level j-resolved?
  if(v.jp!=-10) cout << v.np << ", " << v.lp << ", " << v.sp << ", " << v.jp << ") " << endl;
  else cout << v.np << ", " << v.lp << ", " << v.sp << "), gw:" << gw << endl;
  cout << " %==============================================================%" << endl;
  cout << scientific << " Dnu: " << v.Dnu*1.0e-9 << " GHz, DE: " << v.DE << " eV, lambda21: " 
       << v.lambda21 << " cm"<< endl;
  //cout << " A21: " << v.A21  << " sec^-1, f: " << v.f << endl << endl; 
  cout << scientific << " A21: " << v.A21  << " sec^-1, gw: " << v.gwp << endl << endl; 
  
  return;
}

void Electron_Level_HeI_Triplet_no_j::display_level_information()
{
  cout.precision(16);

  cout << " %==============================================================%" << endl;
  cout << " % Electron_Level_HeI_Triplet_no_j::display_level_information:\n % (n, l, s) == (" 
       << nn << ", " << ll << ", 1), gw:" << gw  
       << endl;
  cout << " %==============================================================%" << endl;
  
  cout << "\n gw: " << gw << endl;
  cout << scientific << " Dnu_1s: " << Dnu*1.0e-9 << " GHz, nuion: " 
       << nuion*1.0e-9 << " GHz" << endl;
  cout << scientific << " DE_1s: " << DE << " eV, Eion: " << Eion << " eV " << endl;
  cout << " Number of transitions: " << Get_n_down() << endl << endl;
  cout << " %==============================================================%" << endl;

  return;
}

void Electron_Level_HeI_Triplet_no_j::display_downward_transition(int i)
{
    if(i>=(int)A_values.size())
    { 
        cout << " Electron_Level_HeI_Triplet_no_j::display_downward_transition::"
             << " no transition tabulated " << i << endl;
        return;
    }
    
    display_Trans_Data(A_values[i]);

    return;
}

void Electron_Level_HeI_Triplet_no_j::display_all_downward_transitions()
{
    for(int l=0; l<(int)A_values.size(); l++) display_Trans_Data(A_values[l]);

    return;
}

//======================================================================
const Transition_Data_HeI_A& Electron_Level_HeI_Triplet_no_j::Get_Trans_Data(int np, int lp, 
                                                                             int sp, int jp) const
{
    if(sp==0 && (int)Transition_Data_HeI_A_lookup_table[0].size()>np)
    {
        int ntrans=Transition_Data_HeI_A_lookup_table[0][np].size();
        
        for(int m=0; m<ntrans; m++) 
        {
            if(Transition_Data_HeI_A_lookup_table[0][np][m]!=-1)
            {
                if(np==A_values[Transition_Data_HeI_A_lookup_table[0][np][m]].np && 
                   sp==A_values[Transition_Data_HeI_A_lookup_table[0][np][m]].sp && 
                   lp==A_values[Transition_Data_HeI_A_lookup_table[0][np][m]].lp && 
                   jp==A_values[Transition_Data_HeI_A_lookup_table[0][np][m]].jp)
                    return A_values[Transition_Data_HeI_A_lookup_table[0][np][m]];      
            }
        }
    }
    if(sp==1 && (int)Transition_Data_HeI_A_lookup_table[1].size()>np)
    {
        int ntrans=Transition_Data_HeI_A_lookup_table[1][np].size();
        
        for(int m=0; m<ntrans; m++) 
        {
            if(Transition_Data_HeI_A_lookup_table[1][np][m]!=-1)
            {
                if(np==A_values[Transition_Data_HeI_A_lookup_table[1][np][m]].np && 
                   sp==A_values[Transition_Data_HeI_A_lookup_table[1][np][m]].sp && 
                   lp==A_values[Transition_Data_HeI_A_lookup_table[1][np][m]].lp && 
                   jp==A_values[Transition_Data_HeI_A_lookup_table[1][np][m]].jp)
                    return A_values[Transition_Data_HeI_A_lookup_table[1][np][m]];      
            }
        }
    }
    
    return ZERO_Data;
}

//======================================================================
double Electron_Level_HeI_Triplet_no_j::Get_A21(int np, int lp, int sp, int jp) const
{
    if(np<=(int)njresolved)
        for(int i=0; i<(int)A_values.size(); i++) 
        {
            if(A_values[i].np==np && A_values[i].lp==lp && A_values[i].sp==sp && A_values[i].jp==jp) 
                return A_values[i].A21;
            //else return 0.0;
        }
    else 
        for(int i=0; i<(int)A_values.size(); i++) 
        {
            if(A_values[i].np==np && A_values[i].lp==lp && A_values[i].sp==sp && A_values[i].jp==-10) 
                return A_values[i].A21;
            //else return 0.0;
        }
    
  return 0.0;
}

double Electron_Level_HeI_Triplet_no_j::Get_nu21(int np, int lp, int sp, int jp) const
{
    for(int i=0; i<(int)A_values.size(); i++) 
    {
        if(A_values[i].np==np && A_values[i].lp==lp && A_values[i].sp==sp && A_values[i].jp==jp) 
            return A_values[i].Dnu;
        //else return 0.0;
    }
    
    return 0.0;
}

double Electron_Level_HeI_Triplet_no_j::Get_lambda21(int np, int lp, int sp, int jp) const
{
    for(int i=0; i<(int)A_values.size(); i++) 
    {
        if(A_values[i].np==np && A_values[i].lp==lp && A_values[i].sp==sp && A_values[i].jp==jp) 
            return A_values[i].lambda21;
        //else return 0.0;
    }
    
    return 0.0;
}

void Electron_Level_HeI_Triplet_no_j::Set_xxx_Transition_Data(int m, int d, double v) 
{ 
  if(d<5)
    {
      A_values[m].xxx[d]=v; 
    }
  else cout << " Electron_Level_HeI_Triplet_no_j::Set_xxx_Transition_Data: check index " << endl;

  return; 
}

//==============================================================================
// Saha-relations with continuum
//==============================================================================
// f(T) == (Ni/[Ne Nc])_LTE
double Electron_Level_HeI_Triplet_no_j::Ni_NeNc_LTE(double TM) const
{ 
    // HeI --> gc=2 --> ge*gc=4
    // 30.07.2015: added reduced mass factor
    return gw/4.0*pow(const_lambdac, 3)*pow(2.0*PI*const_kb_mec2*Get_mu_red()*TM, -1.5)
                 *exp(Eion_ergs/const_kB/TM );
}

double Electron_Level_HeI_Triplet_no_j::Xi_Saha(double Xe, double Xc, double NH, double TM) const
{ return Xe*Xc*NH*Ni_NeNc_LTE(TM); }

double Electron_Level_HeI_Triplet_no_j::Ni_Saha(double Ne, double Nc, double TM) const
{ return Ne*Nc*Ni_NeNc_LTE(TM); }

//###################################################################################
//
// class: Atomic_Shell_HeI_Singlet
//
//###################################################################################

//===================================================================================
//Konstructors and Destructors for class: Atomic_Shell_HeI_Singlet
//===================================================================================
void Atomic_Shell_HeI_Singlet::init(int n, int mflag)
{
    nn=n;
    mess_flag=mflag;
    create_Electron_Levels();
}

Atomic_Shell_HeI_Singlet::Atomic_Shell_HeI_Singlet(int n, int mflag)
{ init(n, mflag); }
    
Atomic_Shell_HeI_Singlet::~Atomic_Shell_HeI_Singlet()
{ Angular_Momentum_Level.clear(); }

//-----------------------------------------------------------------------------------
void Atomic_Shell_HeI_Singlet::create_Electron_Levels()
{
    if(mess_flag>-1)
    {
        cout << "\n %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%\n %" << endl;
        cout << " % Atomic_Shell_HeI_Singlet::create_Electron_Levels: filling shell: " 
             << nn <<"\n % " << endl;
        cout << " %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%\n" << endl;
    }
    
    // free memory (if necessary)
    Angular_Momentum_Level.clear();
    
    Electron_Level_HeI_Singlet v;
    // fill with empty electron-states
    for(int l=0; l<nn; l++) Angular_Momentum_Level.push_back(v);
    // now initialize each state
    for(int l=0; l<nn; l++) Angular_Momentum_Level[l].init(nn, l, mess_flag);
    
    return;
}

void Atomic_Shell_HeI_Singlet::display_general_data_of_level(int i)
{
    if(i>=(int)(nn) || Angular_Momentum_Level.size()==0)
    {
        cout << " This level does not exist inside shell " << nn << endl;
        return;
    }

    cout << " %==============================================================%" << endl;
    cout << " % Atomic_Shell_HeI_Singlet::display_general_data_of_level:\n %\n" << endl; 
  
    Angular_Momentum_Level[i].display_level_information();

    return;
}

//###################################################################################
//
//class: Atomic_Shell_HeI_Triplet
//
//###################################################################################
//===================================================================================
// Konstructors and Destructors for class: Atomic_Shell_HeI_Triplet
//===================================================================================
void Atomic_Shell_HeI_Triplet::init(int n, int mflag)
{
    nn=n;
    mess_flag=mflag;
    create_Electron_Levels();
}

Atomic_Shell_HeI_Triplet::Atomic_Shell_HeI_Triplet(int n, int mflag)
{ init(n, mflag); }
    
Atomic_Shell_HeI_Triplet::~Atomic_Shell_HeI_Triplet()
{ Angular_Momentum_Level.clear(); }

//-----------------------------------------------------------------------------------
void Atomic_Shell_HeI_Triplet::create_Electron_Levels()
{
    if(mess_flag>-1)
    {
        cout << "\n %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%\n %" << endl;
        cout << " % Atomic_Shell_HeI_Triplet::create_Electron_Levels: filling shell: " 
             << nn <<"\n % " << endl;
        cout << " %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%\n" << endl;
    }
    
    // free memory (if necessary)
    Angular_Momentum_Level.clear();

    Electron_Level_HeI_Triplet v;
    vector<Electron_Level_HeI_Triplet> v3;
    v3.push_back(v);
    v3.push_back(v);
    v3.push_back(v);

    // fill with empty electron-states
    for(int l=0; l<nn; l++) Angular_Momentum_Level.push_back(v3);
    // now initialize each state
    // n3S -state
    Angular_Momentum_Level[0][2].init(nn, 0, 1, mess_flag);
    // n3X -states
    for(int l=1; l<nn; l++)
      {
          Angular_Momentum_Level[l][0].init(nn, l, l-1, mess_flag);
          Angular_Momentum_Level[l][1].init(nn, l, l, mess_flag);
          Angular_Momentum_Level[l][2].init(nn, l, l+1, mess_flag);
      }

    return;
}

void Atomic_Shell_HeI_Triplet::display_general_data_of_level(int i)
{
    if(i>=(int)(nn) || Angular_Momentum_Level.size()==0)
    {
        cout << " This level does not exist inside shell " << nn << endl;
        return;
    }
    
    cout << " %==============================================================%" << endl;
    cout << " % Atomic_Shell_HeI_Triplet::display_general_data_of_level:\n %\n" << endl; 
    
    Angular_Momentum_Level[i][2].display_level_information();
    if(i!=0)  
    {
        Angular_Momentum_Level[i][1].display_level_information();
        Angular_Momentum_Level[i][0].display_level_information();
    }
    return;
}

//###################################################################################
//
// class: Atomic_Shell_HeI_Triplet_no_j
//
//###################################################################################

//===================================================================================
// Konstructors and Destructors for class: Atomic_Shell_HeI_Triplet_no_j
//===================================================================================
void Atomic_Shell_HeI_Triplet_no_j::init(int n, int njres, int mflag)
{
    nn=n;
    mess_flag=mflag;

    //===============================================================================
    // above njresolved the triplet states will not be treated j-resolved.
    // this should not be smaller than 10
    //===============================================================================
    njresolved=(int)min(max(njres, 10), n);

    create_Electron_Levels();
}

Atomic_Shell_HeI_Triplet_no_j::Atomic_Shell_HeI_Triplet_no_j(int n, int njres, int mflag)
{ init(n, njres, mflag); }
    
Atomic_Shell_HeI_Triplet_no_j::~Atomic_Shell_HeI_Triplet_no_j()
{ Angular_Momentum_Level.clear(); }

//-----------------------------------------------------------------------------------
void Atomic_Shell_HeI_Triplet_no_j::create_Electron_Levels()
{
    if(mess_flag>-1)
    {
        cout << "\n %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%\n %" << endl;
        cout << " % Atomic_Shell_HeI_Triplet_no_j::create_Electron_Levels: filling shell: " 
             << nn <<"\n % " << endl;
        cout << " % The levels will not be j-resolved \n % " << endl;
        cout << " %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%\n" << endl;
    }
    
    // free memory (if necessary)
    Angular_Momentum_Level.clear();

    Electron_Level_HeI_Triplet_no_j v;
    // fill with empty electron-states
    for(int l=0; l<nn; l++) Angular_Momentum_Level.push_back(v);
    // now initialize each state
    for(int l=0; l<nn; l++) Angular_Momentum_Level[l].init(nn, l, njresolved, mess_flag);

    return;
}

void Atomic_Shell_HeI_Triplet_no_j::display_general_data_of_level(int i)
{
    if(i>=(int)(nn) || Angular_Momentum_Level.size()==0)
    {
        cout << " This level does not exist inside shell " << nn << endl;
        return;
    }
    
    cout << " %==============================================================%" << endl;
    cout << " % Atomic_Shell_HeI_Triplet_no_j::display_general_data_of_level:\n %\n" << endl; 
    
    Angular_Momentum_Level[i].display_level_information();
    
    return;
}

//###################################################################################
//
// class: Atom_HeI_Singlet
//
//###################################################################################

//===================================================================================
// Konstructors and Destructors for Atom_HeI_Singlet
//===================================================================================
void Atom_HeI_Singlet::init(int nS, int mflag)
{
    nShells=nS;
    mess_flag=mflag;

    fill_Level_Map();

    create_Shells();
    if(mess_flag>2) display_Level_Map();
}

Atom_HeI_Singlet:: Atom_HeI_Singlet(){}

Atom_HeI_Singlet:: Atom_HeI_Singlet(int nS, int mflag){ init(nS, mflag); }

Atom_HeI_Singlet::~Atom_HeI_Singlet()
{ 
  Level_Map.clear(); 
  Shell.clear();
}

//-----------------------------------------------------------------------------------
void Atom_HeI_Singlet::create_Shells()
{
    if(mess_flag>-1)
    {
        cout << "\n %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%\n %" << endl;
        cout << " % Atom_HeI_Singlet::create_Shells:creating all the shells up to n=" 
             << nShells << endl;
        cout << " %\n %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%\n" << endl;
    }
    
    int m=mess_flag-1;
    if(mess_flag>=1) m++;
    if(mess_flag>=2) m++;
    
    // free memory (if necessary)
    Shell.clear();

    Atomic_Shell_HeI_Singlet v;
    // fill with empty shells
    for(int n=0; n<=nShells; n++) Shell.push_back(v);
    // create each shells
    for(int n=1; n<=nShells; n++) Shell[n].init(n, m);
    
    // add hydrogenic transitions for DM singlet gap
    if(nShells>=9)
    {
        Level(9, 8).Set_hydrogenic_transitions(*this);
    
        if(nShells>=10)
        {
            Level(10, 7).Set_hydrogenic_transitions(*this);
            Level(10, 8).Set_hydrogenic_transitions(*this);
            Level(10, 9).Set_hydrogenic_transitions(*this);
        }
    }
    
    // add all hydrogenic transitions for levels with n>10
    for(int n=11; n<=nShells; n++)
        for(int l=0; l<n; l++) Level(n, l).Set_hydrogenic_transitions(*this);
    
    return;
}

//------------------------------------------------------------------------------------------------------
void Atom_HeI_Singlet::fill_Level_Map()
{
    // free memory (if necessary)
    Level_Map.clear(); 
    
    // this function creates the map of indicies i-> (n,l)
    Level_nl v;
    
    for(int n=1; n<=nShells; n++) 
        for(int l=0; l<n; l++)
        {
            v.n=n; v.l=l;         
            Level_Map.push_back(v);
        }
    
    return;
}

void Atom_HeI_Singlet::display_Level_Map()
{
    cout << "\n %==============================================================%" << endl;
    cout << " % Atom_HeI_Singlet::display_Level_Map:" << endl; 
    cout << " % Format: i [check of inversion] --> (n, l, s, j)" << endl; 
    cout << " %==============================================================%" << endl;
    
    for(int i=0; i<(int)Level_Map.size(); i++)
    {
        if(Get_Level_index(Level_Map[i].n, 0)==(int)i) cout << " Shell " << Level_Map[i].n << endl;
        
        cout << " " << i << " [" << Get_Level_index(Level_Map[i].n, Level_Map[i].l) << "] --> (" 
             << Level_Map[i].n << ", " << Level_Map[i].l << ", 0, " 
             << Level_Map[i].l <<")" << endl; 
    }
    cout << endl;
    
    return;
}

//================================================================================
// to check level access
//================================================================================
void Atom_HeI_Singlet::check_Level(int n, int l, string mess) const
{
    if(n<1 || n>nShells || l>=n || l<0)
    {
        cout << " Atom_HeI_Singlet::" << mess << ": you are trying to access a non-existing level: " << n 
             << ", " << l << ", 0, " << l << endl;
        exit(0);
    }
    return;
}

void Atom_HeI_Singlet::check_Level(int i, string mess) const
{
    if(i> Get_total_number_of_Levels())
    {
        cout << " Atom_HeI_Singlet::" << mess << ": you are trying to access a non-existing level: " << i 
             << " total number of levels: " << Get_total_number_of_Levels() << endl;
        exit(0);
    }
    return;
}

const Electron_Level_HeI_Singlet& Atom_HeI_Singlet::Level(int i) const
{ 
    check_Level(i, "Level");

    //cout << " Atom_HeI_Singlet: " << i << " " << Level_Map[i].n << " " << Level_Map[i].l << endl;

    return Shell[Level_Map[i].n].Angular_Momentum_Level[Level_Map[i].l]; 
}
 
const Electron_Level_HeI_Singlet& Atom_HeI_Singlet::Level(int n, int l) const
{ 
    check_Level(n, l, "Level");

    return Shell[n].Angular_Momentum_Level[l]; 
}
 
Electron_Level_HeI_Singlet& Atom_HeI_Singlet::Level(int i) 
{ 
   check_Level(i, "Level");

   //cout << " Atom_HeI_Singlet: " << i << " " << Level_Map[i].n << " " << Level_Map[i].l << endl;

   return Shell[Level_Map[i].n].Angular_Momentum_Level[Level_Map[i].l]; 
}
 
Electron_Level_HeI_Singlet& Atom_HeI_Singlet::Level(int n, int l)
{ 
    check_Level(n, l, "Level");

    return Shell[n].Angular_Momentum_Level[l]; 
}
 
//###################################################################################
//
// class: Atom_HeI_Triplet
//
//###################################################################################

//===================================================================================
// Konstructors and Destructors for Atom_HeI_Triplet
//===================================================================================
void Atom_HeI_Triplet::init(int nS, int mflag)
{
    nShells=nS;
    mess_flag=mflag;

    fill_Level_Map();
    create_Shells();

    if(mess_flag>2) display_Level_Map();
}

Atom_HeI_Triplet:: Atom_HeI_Triplet(){}

Atom_HeI_Triplet:: Atom_HeI_Triplet(int nS, int mflag){ init(nS, mflag); }

Atom_HeI_Triplet::~Atom_HeI_Triplet()
{ 
  Level_Map.clear(); 
  Shell.clear();
}

//-----------------------------------------------------------------------------------
void Atom_HeI_Triplet::create_Shells()
{
    if(mess_flag>-1)
    {
        cout << "\n %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%\n %" << endl;
        cout << " % Atom_HeI_Triplet::create_Shells:creating all the shells up to n=" 
             << nShells << endl;
        cout << " %\n %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%\n" << endl;
    }
    
    int m=mess_flag-1;
    if(mess_flag>=1) m++;
    if(mess_flag>=2) m++;
    
    // free memory (if necessary)
    Shell.clear();
    
    Atomic_Shell_HeI_Triplet v;
    // fill with empty shells
    for(int n=0; n<=nShells; n++) Shell.push_back(v);
    // create each shells
    for(int n=2; n<=nShells; n++) Shell[n].init(n, m); 
    
    // add hydrogenic transitions for DM triplet gap
    if(nShells>=8)
    {
        for(int j=7-1; j<=7+1; j++) Level(8, 7, j).Set_hydrogenic_transitions(*this);

        if(nShells>=9)
        {
            for(int j=8-1; j<=8+1; j++) Level(9, 8, j).Set_hydrogenic_transitions(*this);
            for(int j=6-1; j<=6+1; j++) Level(9, 6, j).Set_hydrogenic_transitions(8, 7, *this);

            if(nShells>=10)
            {
                for(int j=7-1; j<=7+1; j++) Level(10, 7, j).Set_hydrogenic_transitions(*this);
                for(int j=8-1; j<=8+1; j++) Level(10, 8, j).Set_hydrogenic_transitions(*this);
                for(int j=9-1; j<=9+1; j++) Level(10, 9, j).Set_hydrogenic_transitions(*this);
                for(int j=6-1; j<=6+1; j++) Level(10, 6, j).Set_hydrogenic_transitions(8, 7, *this);
            }
        }
    }
    
    // add all hydrogenic transitions for levels with n>10
    for(int n=11; n<=nShells; n++)
        for(int l=0; l<n; l++)
            for(int j=l-1; j<=l+1; j++)
                Level(n, l, j).Set_hydrogenic_transitions(*this);

    return;
}

int Atom_HeI_Triplet::Get_Level_index(int n, int l, int j) const 
{ 
    if(n<2 || n>nShells || l>=n || l<0 || l+1<j || l>j+1 || (l==0 && j!=1)) 
    {
        cerr << " Atom_HeI_Triplet::Get_Level_index: No level with (n, l, s, j)= ("
             << n << ", " << l << ", " << 1 << ", " << j << ")!!! " << endl;
        exit(0);
    }

  return Get_number_of_Levels_until(n-1)+3*(l+1)-2-(l+1-j)-1;
}

//------------------------------------------------------------------------------------------------------
void Atom_HeI_Triplet::fill_Level_Map()
{
    // free memory (if necessary)
    Level_Map.clear(); 
    
    // this function creates the map of indicies i-> (n,l)
    Level_nl v;
    
    for(int n=2; n<=nShells; n++) 
    {
        v.n=n; v.l=0; v.j=1;      
        Level_Map.push_back(v);
        
        for(int l=1; l<n; l++)
        {
            v.n=n; v.l=l; v.j=l-1;        
            Level_Map.push_back(v);
            v.n=n; v.l=l; v.j=l;      
            Level_Map.push_back(v);
            v.n=n; v.l=l; v.j=l+1;        
            Level_Map.push_back(v);
        }
    }
    
    return;
}

void Atom_HeI_Triplet::display_Level_Map()
{
    cout << "\n %==============================================================%" << endl;
    cout << " % Atom_HeI_Triplet::display_Level_Map:" << endl; 
    cout << " % Format: i --> (n, l, s, j)" << endl; 
    cout << " %==============================================================%" << endl;
    
    for(int i=0; i<(int)Level_Map.size(); i++)
    {
        if(Get_Level_index(Level_Map[i].n, 0, 1)==(int)i) cout << " Shell " << Level_Map[i].n << endl;
        
        cout << " " << i << " [" << Get_Level_index(Level_Map[i].n, Level_Map[i].l, Level_Map[i].j) << "] --> (" 
             << Level_Map[i].n << ", " << Level_Map[i].l << ", 1, " 
             << Level_Map[i].j <<")" << endl; 
    }
    
    cout << endl;
    
    return;
}

//================================================================================
// to check level access
//================================================================================
void Atom_HeI_Triplet::check_Level(int n, int l, int j, string mess) const
{
    if(n<2 || n>nShells || l>=n || l<0 || l+1<j || l>j+1 || (l==0 && j!=1))
    {
        cout << " Atom_HeI_Triplet::" << mess 
             << ": you are trying to access a non-existing level: " << n 
             << ", " << l << ", 1, " << j << endl;
        exit(0);
    }
    return;
}

void Atom_HeI_Triplet::check_Level(int i, string mess) const
{
    if(i> Get_total_number_of_Levels())
    {
        cout << " Atom_HeI_Triplet::" << mess 
             << ": you are trying to access a non-existing level: " << i 
             << " total number of levels: " << Get_total_number_of_Levels() << endl;
        exit(0);
    }
    return;
}

const Electron_Level_HeI_Triplet& Atom_HeI_Triplet::Level(int i) const
{ 
    check_Level(i, "Level");

    // the access of data is such that for l==0 --> j-(l-1)=2
    return Shell[Level_Map[i].n].Angular_Momentum_Level[Level_Map[i].l][Level_Map[i].j-(Level_Map[i].l-1)]; 
}
 
const Electron_Level_HeI_Triplet& Atom_HeI_Triplet::Level(int n, int l, int j) const
{ 
    check_Level(n, l, j, "Level");

    // the access of data is such that for l==0 --> j-(l-1)=2
    return Shell[n].Angular_Momentum_Level[l][j-(l-1)]; 
}
 
Electron_Level_HeI_Triplet& Atom_HeI_Triplet::Level(int i) 
{ 
    check_Level(i, "Level");

    // the access of data is such that for l==0 --> j-(l-1)=2
    return Shell[Level_Map[i].n].Angular_Momentum_Level[Level_Map[i].l][Level_Map[i].j-(Level_Map[i].l-1)]; 
}
 
Electron_Level_HeI_Triplet& Atom_HeI_Triplet::Level(int n, int l, int j)
{ 
    check_Level(n, l, j, "Level");

    // the access of data is such that for l==0 --> j-(l-1)=2
    return Shell[n].Angular_Momentum_Level[l][j-(l-1)]; 
}

//###################################################################################
//
// class: Atom_HeI_Triplet_no_j
//
//###################################################################################

//===================================================================================
// Konstructors and Destructors for Atom_HeI_Triplet_no_j
//===================================================================================
void Atom_HeI_Triplet_no_j::init(int nS, int njres, int mflag)
{
    nShells=nS;
    mess_flag=mflag;

    //===============================================================================
    // above njresolved the triplet states will not be treated j-resolved.
    // this should not be smaller than 10
    //===============================================================================
    njresolved=(int)min(max(njres, 10), nS);

    fill_Level_Map();
    create_Shells();

    if(mess_flag>2) display_Level_Map();
}

Atom_HeI_Triplet_no_j:: Atom_HeI_Triplet_no_j(){}

Atom_HeI_Triplet_no_j:: Atom_HeI_Triplet_no_j(int nS, int njres, int mflag){ init(nS, njres, mflag); }

Atom_HeI_Triplet_no_j::~Atom_HeI_Triplet_no_j()
{ 
  Level_Map.clear(); 
  Shell.clear();
}

//-----------------------------------------------------------------------------------
void Atom_HeI_Triplet_no_j::create_Shells()
{
    if(mess_flag>-1)
    {
        cout << "\n %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%\n %" << endl;
        cout << " % Atom_HeI_Triplet_no_j::create_Shells:creating all the shells from " 
             << njresolved+1 << " up to n=" << nShells << endl;
        cout << " %\n %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%\n" << endl;
    }
    
    int m=mess_flag-1;
    if(mess_flag>=1) m++;
    if(mess_flag>=2) m++;
    
    // free memory (if necessary)
    Shell.clear();
    
    Atomic_Shell_HeI_Triplet_no_j v;
    // fill with empty shells
    for(int n=0; n<=nShells; n++) Shell.push_back(v);
    // create each shells
    for(int n=njresolved+1; n<=nShells; n++) Shell[n].init(n, njresolved, m);
    
    // add hydrogenic transitions for non-j-resolved --> j-resolved & non-j-resolved
    for(int n=njresolved+1; n<=nShells; n++)
        for(int l=0; l<n; l++) Level(n, l).Set_hydrogenic_transitions(*this);
    
    return;
}

int Atom_HeI_Triplet_no_j::Get_number_of_Levels_until(int nmax) const
{ return nmax*(nmax+1)/2 - njresolved*(njresolved+1)/2; }

int Atom_HeI_Triplet_no_j::Get_Level_index(int n, int l) const 
{ 
    if(n>nShells || n<=njresolved || l>=n || l<0)
    {
        cerr << " Atom_HeI_Triplet_no_j::Get_Level_index: No level with (n, l, s, j)= ("
             << n << ", " << l << ", " << 1 << ", " << -10 << ")!!! " << endl;
        exit(0);
    }
    
    return Get_number_of_Levels_until(n-1)+l; 
}

//------------------------------------------------------------------------------------------------------
void Atom_HeI_Triplet_no_j::fill_Level_Map()
{
    // free memory (if necessary)
    Level_Map.clear(); 

    // this function creates the map of indicies i-> (n,l)
    Level_nl v;
    
    for(int n=njresolved+1; n<=nShells; n++)
      for(int l=0; l<n; l++)
    {
      v.n=n; v.l=l;     
      Level_Map.push_back(v);
    }
   
    return;
}

void Atom_HeI_Triplet_no_j::display_Level_Map()
{
    cout << "\n %==============================================================%" << endl;
    cout << " % Atom_HeI_Triplet_no_j::display_Level_Map:" << endl; 
    cout << " % Format: i --> (n, l, s, -10)" << endl; 
    cout << " %==============================================================%" << endl;
    
    for(int i=0; i<(int)Level_Map.size(); i++)
    {
        if(Get_Level_index(Level_Map[i].n, 0)==(int)i) cout << " Shell " << Level_Map[i].n << endl;
        
        cout << " " << i << " [" << Get_Level_index(Level_Map[i].n, Level_Map[i].l) << "] --> (" 
             << Level_Map[i].n << ", " << Level_Map[i].l << ", 1, -10)" << endl; 
    }
    
    cout << endl;
    
    return;
}

//================================================================================
// to check level access
//================================================================================
void Atom_HeI_Triplet_no_j::check_Level(int n, int l, string mess) const
{
    if(n>nShells || n<=njresolved || l>=n || l<0)
    {
        cout << " Atom_HeI_Triplet_no_j::" << mess 
             << ": you are trying to access a non-existing level: " << n 
             << ", " << l << ", 1 " << endl;
        exit(0);
    }
    return;
}

void Atom_HeI_Triplet_no_j::check_Level(int i, string mess) const
{
    if(i> Get_total_number_of_Levels())
    {
        cout << " Atom_HeI_Triplet_no_j::" << mess 
             << ": you are trying to access a non-existing level: " << i 
             << " total number of levels: " << Get_total_number_of_Levels() << endl;
        exit(0);
    }
    return;
}

const Electron_Level_HeI_Triplet_no_j& Atom_HeI_Triplet_no_j::Level(int i) const
{ 
    check_Level(i, "Level");

    return Shell[Level_Map[i].n].Angular_Momentum_Level[Level_Map[i].l]; 
}
 
const Electron_Level_HeI_Triplet_no_j& Atom_HeI_Triplet_no_j::Level(int n, int l) const
{ 
    check_Level(n, l, "Level");

    return Shell[n].Angular_Momentum_Level[l]; 
}
 
Electron_Level_HeI_Triplet_no_j& Atom_HeI_Triplet_no_j::Level(int i) 
{ 
    check_Level(i, "Level");

    return Shell[Level_Map[i].n].Angular_Momentum_Level[Level_Map[i].l];
}
 
Electron_Level_HeI_Triplet_no_j& Atom_HeI_Triplet_no_j::Level(int n, int l)
{ 
    check_Level(n, l, "Level");

    return Shell[n].Angular_Momentum_Level[l]; 
}

//###################################################################################
//
// class: Gas_of_HeI_Atoms
//
//###################################################################################

//===================================================================================
// Konstructors and Destructors for class: Gas_of_HeI_Atoms
//===================================================================================
void Gas_of_HeI_Atoms::init(int nS, int njres, int nQ, int nTS, int mflag)
{ 
    //===============================================================================
    // do not read transition data but only energies!
    //===============================================================================
    if(nS<0){ HeI_Atom_read_transition_data=0; nS=-nS; }
    
    //===============================================================================
    // above njresolved the triplet states will not be treated j-resolved.
    // this should not be smaller than 10
    //===============================================================================
    njresolved=(int)min(max(njres, 10), nS);
    
    //===============================================================================
    // set information flags
    //===============================================================================
    add_Q_lines_loaded=add_TS_lines_loaded=0;
    n_HeI_add_Intercomb=max(0, min(nTS, 10));
    n_HeI_add_Quad=max(0, min(nQ, 10));
    
    if(n_HeI_add_Quad>0) add_Q_lines_loaded=1;
    if(n_HeI_add_Intercomb>0) add_TS_lines_loaded=1;
    
    _HeI_add_Intercomb=n_HeI_add_Intercomb;
    _HeI_add_Quad=n_HeI_add_Quad;
    
    //===============================================================================
    Sing.init(nS, mflag);
    Trip.init((int)min(nS, njresolved), mflag);
    Trip_glob=&Trip;
    if(nS>(int)njresolved) Trip_no_j.init(nS, njresolved, mflag);
    Trip_glob=NULL;
    
    indexT=Sing.Get_total_number_of_Levels();
    nl=indexT_no_j=indexT+Trip.Get_total_number_of_Levels();
    if(nS>(int)njresolved) nl+=Trip_no_j.Get_total_number_of_Levels();
    
    if(HeI_Atom_read_transition_data!=0) 
    {
        //============================================================
        // check the consistency of transition data (17.05.2009)
        //============================================================
        check_transition_data();        
        
        //===========================================================
        // create the Voigt profiles
        //===========================================================
        double A21, f=0.0, nu21, lam21, Gamma;
       
        //===========================================================
        // Transition-Data for n^1 P_1 - 1^1 S_0
        //===========================================================
        for(int n=2; n<=(int)min(nS, 10); n++)
        {
            A21=Sing.Level(n, 1).Get_A21(1, 0, 0, 0);
            nu21=Sing.Level(n, 1).Get_nu21(1, 0, 0, 0);
            lam21=Sing.Level(n, 1).Get_lambda21(1, 0, 0, 0);
            
            Gamma=0.0;
            for(int m=0; m<(int)Sing.Level(n, 1).Get_n_down(); m++) 
                Gamma+=Sing.Level(n, 1).Get_Trans_Data(m).A21;
            
            if(mflag>=1) cout << " Initializing profile for " << n 
                              << "^1 P_1 - 1^1 S_0 transition:" 
                              << nu21 << " " << lam21 << " " << A21 << " " 
                              << Gamma << " " << f << endl;

            phi_HeI_nP_S[n].Set_atomic_data_Gamma(nu21, lam21, A21, f, Gamma, 4);
        }
        
        //===========================================================
        // Transition-Data for n^3 P_1 - 1^1 S_0
        //===========================================================
        for(int n=2; n<=(int)min(nS, n_HeI_add_Intercomb); n++)
        {
            A21=Trip.Level(n, 1, 1).Get_A21(1, 0, 0, 0);
            nu21=Trip.Level(n, 1, 1).Get_nu21(1, 0, 0, 0);
            lam21=Trip.Level(n, 1, 1).Get_lambda21(1, 0, 0, 0);
            
            Gamma=0.0;
            for(int m=0; m<(int)Trip.Level(n, 1, 1).Get_n_down(); m++) 
                Gamma+=Trip.Level(n, 1, 1).Get_Trans_Data(m).A21;
            
            if(mflag>=1) cout << " Initializing profile for " << n 
                              << "^3 P_1 - 1^1 S_0 transition:" 
                              << nu21 << " " << lam21 << " " << A21 << " " 
                              << Gamma << " " << f << endl;

            phi_HeI_nP_T[n].Set_atomic_data_Gamma(nu21, lam21, A21, f, Gamma, 4);
        }

        //===========================================================
        // Transition-Data for n^1 D_1 - 1^1 S_0
        //===========================================================
        for(int n=3; n<=(int)min(nS, n_HeI_add_Quad); n++)
        {
            A21=Sing.Level(n, 2).Get_A21(1, 0, 0, 0);
            nu21=Sing.Level(n, 2).Get_nu21(1, 0, 0, 0);
            lam21=Sing.Level(n, 2).Get_lambda21(1, 0, 0, 0);
            
            Gamma=0.0;
            for(int m=0; m<(int)Sing.Level(n, 2).Get_n_down(); m++) 
                Gamma+=Sing.Level(n, 2).Get_Trans_Data(m).A21;
            
            if(mflag>=1) cout << " Initializing profile for " << n 
                              << "^1 D_1 - 1^1 S_0 transition:" 
                              << nu21 << " " << lam21 << " " << A21 << " " 
                              << Gamma << " " << f << endl;
            
            phi_HeI_nD_S[n].Set_atomic_data_Gamma(nu21, lam21, A21, f, Gamma, 4);
        }
    }
    
    //========================================
    // create HI Ly-c cross section
    //========================================
    HILyc.init(1, 0, 1, 1, 0);

    //========================================
    // restore default state for data read
    //========================================
    HeI_Atom_read_transition_data=1;
    
    _HeI_add_Intercomb=_HeI_add_Quad=0;
    
    return;
}

Gas_of_HeI_Atoms::Gas_of_HeI_Atoms(){}

Gas_of_HeI_Atoms::Gas_of_HeI_Atoms(int nS, int njres, int mflag)
{ init(nS, njres, 0, 0, mflag); }

void Gas_of_HeI_Atoms::init(int nS, int njres, int mflag)
{
    init(nS, njres, 0, 0, mflag);
    return;
}

Gas_of_HeI_Atoms::Gas_of_HeI_Atoms(int nS, int njres, int nQ, int nTS, int mflag)
{ init(nS, njres, nQ, nTS, mflag); }


    
Gas_of_HeI_Atoms::~Gas_of_HeI_Atoms(){}

//===================================================================================
// Access to Voigt-profiles
//===================================================================================
const Voigtprofile_Dawson& Gas_of_HeI_Atoms::nP_S_profile(int n) const 
{ 
    if(n>1 && n<=(int)min(Get_nShells(), 10)) return phi_HeI_nP_S[n]; 
    
    return phi_HeI_nP_S[0]; 
}

Voigtprofile_Dawson& Gas_of_HeI_Atoms::nP_S_profile(int n)
{ 
    if(n>1 && n<=(int)min(Get_nShells(), 10)) return phi_HeI_nP_S[n]; 
    
    return phi_HeI_nP_S[0]; 
}

//===================================================================================
const Voigtprofile_Dawson& Gas_of_HeI_Atoms::nP_T_profile(int n) const 
{ 
    if(n>1 && n<=(int)min(Get_nShells(), n_HeI_add_Intercomb)) return phi_HeI_nP_T[n];
    
    return phi_HeI_nP_T[0]; 
}

Voigtprofile_Dawson& Gas_of_HeI_Atoms::nP_T_profile(int n)
{ 
    if(n>1 && n<=(int)min(Get_nShells(), n_HeI_add_Intercomb)) return phi_HeI_nP_T[n];
    
    return phi_HeI_nP_T[0]; 
}

//===================================================================================
const Voigtprofile_Dawson& Gas_of_HeI_Atoms::nD_S_profile(int n) const 
{ 
    if(n>1 && n<=(int)min(Get_nShells(), n_HeI_add_Quad)) return phi_HeI_nD_S[n];
    
    return phi_HeI_nD_S[0]; 
}

Voigtprofile_Dawson& Gas_of_HeI_Atoms::nD_S_profile(int n)
{ 
    if(n>1 && n<=(int)min(Get_nShells(), n_HeI_add_Quad)) return phi_HeI_nD_S[n];
    
    return phi_HeI_nD_S[0]; 
}

//===================================================================================
void Gas_of_HeI_Atoms:: check_transition_data()
{
    double nucu, nucl, Dnu, val;
    int ip, count=0;
    
    for(int k=0; k<(int)nl; k++)
    {
        nucu=Get_nu_ion(k);
        for(int m=0; m<(int)Get_n_down(k); m++)
        {
            ip=Get_Level_index(Get_Trans_Data(k, m).np, Get_Trans_Data(k, m).lp, Get_Trans_Data(k, m).sp, Get_Trans_Data(k, m).jp);
            
            nucl=Get_nu_ion(ip);
            Dnu=nucl-nucu;
            
            val=Get_Trans_Data(k, m).Dnu;
            
            if(fabs(Dnu/val-1.0)>=1.0e-8) 
            {
                count++;
                cout << "\n Level " << k << " == (" << Get_n(k) << ", " <<  Get_l(k) << ", "
                                                    << Get_S(k) << ", " <<  Get_J(k) << ")" << endl;
                cout << " There is an inconsistency of " << fabs(Dnu/val-1.0)
                     << " in the transition frequency " << val << " Hz to level "
                     << ip << " == (" << Get_Trans_Data(k, m).np << ", " << Get_Trans_Data(k, m).lp << ", "
                                      << Get_Trans_Data(k, m).sp << ", " << Get_Trans_Data(k, m).jp << ")" << endl;
            }
        }
    }
    
    if(count!=0)
    {
        cout << " There were " << count
             << " inconsistencies in the transition frequencies. Maybe you want to recompute your Helium model... "
             << endl;
        wait_f_r();
    }
    
    return;
}

//===================================================================================
const Transition_Data_HeI_A& Gas_of_HeI_Atoms::Get_Trans_Data(int i, int np, int lp,
                                                              int sp, int jp) const
{
    int n=Get_n(i);
    int l=Get_l(i);
    int j=Get_J(i);
    int s=Get_S(i);
    
    if(s==0) return Sing.Level(n, l).Get_Trans_Data(np, lp, sp, jp);
    else if(s==1 && n<=(int)njresolved)
        return Trip.Level(n, l, j).Get_Trans_Data(np, lp, sp, jp);
    else if(s==1 && n>(int)njresolved)
        return Trip_no_j.Level(n, l).Get_Trans_Data(np, lp, sp, jp);
    
    cout << " Gas_of_HeI_Atoms:: fail " << endl;
    
    return Sing.Level(i).Get_Trans_Data(0, lp, sp, jp);
}

double Gas_of_HeI_Atoms:: Get_A(int n, int l, int s, int j, int np, int lp, int sp, int jp) const
{
  if(s==0) return Sing.Level(n, l).Get_A21(np, lp, sp, jp);
  else if(s==1 && n<=njresolved) return Trip.Level(n, l, j).Get_A21(np, lp, sp, jp);
  else if(s==1 && n>njresolved) return Trip_no_j.Level(n, l).Get_A21(np, lp, sp, jp);

  cout << " Gas_of_HeI_Atoms:: fail " << endl;
    
  return 0.0;
}

double Gas_of_HeI_Atoms:: Get_nu21(int n, int l, int s, int j,
                                   int np, int lp, int sp, int jp) const
{
    if(s==0) return Sing.Level(n, l).Get_nu21(np, lp, sp, jp);
    else if(s==1 && n<=njresolved) return Trip.Level(n, l, j).Get_nu21(np, lp, sp, jp);
    else if(s==1 && n>njresolved) return Trip_no_j.Level(n, l).Get_nu21(np, lp, sp, jp);
    
    cout << " Gas_of_HeI_Atoms:: fail " << endl;
    return 0.0;
}

double Gas_of_HeI_Atoms:: Get_lambda21(int n, int l, int s, int j, 
                                       int np, int lp, int sp, int jp) const
{
    if(s==0) return Sing.Level(n, l).Get_lambda21(np, lp, sp, jp);
    else if(s==1 && n<=njresolved) return Trip.Level(n, l, j).Get_lambda21(np, lp, sp, jp);
    else if(s==1 && n>njresolved) return Trip_no_j.Level(n, l).Get_lambda21(np, lp, sp, jp);
    
    cout << " Gas_of_HeI_Atoms:: fail " << endl;
    return 0.0;
}

double Gas_of_HeI_Atoms::Xi(int n, int l, int s, int j) const
{
  if(s==0) return Sing.Level(n, l).Get_Xi();
  else if(s==1 && n<=njresolved) return Trip.Level(n, l, j).Get_Xi();
  else if(s==1 && n>njresolved) return Trip_no_j.Level(n, l).Get_Xi();

  cout << " Gas_of_HeI_Atoms:: fail " << endl;
  return 0.0;
}

double Gas_of_HeI_Atoms::Ric(int n, int l, int s, int j) const
{
  if(s==0) return Sing.Level(n, l).Get_Ric();
  else if(s==1 && n<=njresolved) return Trip.Level(n, l, j).Get_Ric();
  else if(s==1 && n>njresolved) return Trip_no_j.Level(n, l).Get_Ric();

  cout << " Gas_of_HeI_Atoms:: fail " << endl;
  return 0.0;
}

void Gas_of_HeI_Atoms::Set_Xi(int n, int l, int s, int j, double X)
{
  if(s==0) Sing.Level(n, l).Set_Xi(X);
  else if(s==1 && n<=njresolved) Trip.Level(n, l, j).Set_Xi(X);
  else if(s==1 && n>njresolved) Trip_no_j.Level(n, l).Set_Xi(X);
  return;
}

void Gas_of_HeI_Atoms::Set_Ric(int n, int l, int s, int j, double X)
{
  if(s==0) Sing.Level(n, l).Set_Ric(X);
  else if(s==1 && n<=njresolved) Trip.Level(n, l, j).Set_Ric(X);
  else if(s==1 && n>njresolved) Trip_no_j.Level(n, l).Set_Ric(X);
  return;
}

//=================================================================================================
// show level info
//=================================================================================================
void Gas_of_HeI_Atoms::Show_NLSJ(int i) const
{
    cout << " Gas_of_HeI_Atoms::Show_NLSJ: (N, L, S, J) == (" 
         << Get_n(i) << ", " << Get_l(i) << ", " 
         << Get_S(i) << ", " << Get_J(i) << ") " 
         << endl;
    
    return;
}

//=================================================================================================
// versions that follow Singlet & Triplet sub-sequently
//=================================================================================================
int Gas_of_HeI_Atoms::Get_Level_index(int n, int l, int s, int j) const
{
    if(s==0 && j==l) return Sing.Get_Level_index(n, l);
    else if(s==1 && n<=njresolved) return Trip.Get_Level_index(n, l, j)+indexT;
    else if(s==1 && n>njresolved) return Trip_no_j.Get_Level_index(n, l)+indexT_no_j;
    
    cerr << " Get_as_of_HeI_Atoms::Get_Level_index: No level with (n, l, s, j)= ("
         << n << ", " << l << ", " << s << ", " << j << ")!!! " << endl;
    exit(0);
    
    return 0; 
}

//===================================================================================
int Gas_of_HeI_Atoms::Get_n(int i) const
{
  if(i<indexT) return Sing.Level(i).Get_n();
  else if(i<indexT_no_j) return Trip.Level(i-indexT).Get_n();
  else if(i<nl) return Trip_no_j.Level(i-indexT_no_j).Get_n();

  cout << " Gas_of_HeI_Atoms:: fail " << endl;
  return 1;
}

int Gas_of_HeI_Atoms::Get_l(int i) const
{
  if(i<indexT) return Sing.Level(i).Get_l();
  else if(i<indexT_no_j) return Trip.Level(i-indexT).Get_l();
  else if(i<nl) return Trip_no_j.Level(i-indexT_no_j).Get_l();

  cout << " Gas_of_HeI_Atoms:: fail " << endl;
  return 1;
}

int Gas_of_HeI_Atoms::Get_S(int i) const
{
  if(i<indexT) return Sing.Level(i).Get_S();
  else if(i<indexT_no_j) return Trip.Level(i-indexT).Get_S();
  else if(i<nl) return Trip_no_j.Level(i-indexT_no_j).Get_S();

  cout << " Gas_of_HeI_Atoms:: fail " << endl;
  return 1;
}

int Gas_of_HeI_Atoms::Get_J(int i) const
{
  if(i<indexT) return Sing.Level(i).Get_J();
  else if(i<indexT_no_j) return Trip.Level(i-indexT).Get_J();
  else if(i<nl) return -10;

  cout << " Gas_of_HeI_Atoms:: fail " << endl;
  return 1;
}

double Gas_of_HeI_Atoms::Get_gw(int i) const
{
  if(i<indexT) return Sing.Level(i).Get_gw();
  else if(i<indexT_no_j) return Trip.Level(i-indexT).Get_gw();
  else if(i<nl) return Trip_no_j.Level(i-indexT_no_j).Get_gw();

  cout << " Gas_of_HeI_Atoms:: fail " << endl;
  return 0.0;
}

double Gas_of_HeI_Atoms::Get_nu_ion(int i) const
{
    if(i<indexT) return Sing.Level(i).Get_nu_ion();
    else if(i<indexT_no_j) return Trip.Level(i-indexT).Get_nu_ion();
    else if(i<nl) return Trip_no_j.Level(i-indexT_no_j).Get_nu_ion();
    
    cout << " Gas_of_HeI_Atoms:: Get_nu_ion: fail " << endl;
    return 0.0;
}
    
//===================================================================================
int Gas_of_HeI_Atoms::Get_n_down() const
{
  int s=0;
  for(int i=0; i<indexT; i++) s+=Sing.Level(i).Get_n_down();
  for(int i=indexT; i<indexT_no_j; i++) s+=Trip.Level(i-indexT).Get_n_down();
  for(int i=indexT_no_j; i<nl; i++) s+=Trip_no_j.Level(i-indexT_no_j).Get_n_down();

  return s;
}

int Gas_of_HeI_Atoms::Get_n_down(int i) const
{
  if(i<indexT) return Sing.Level(i).Get_n_down();
  else if(i<indexT_no_j) return Trip.Level(i-indexT).Get_n_down();
  else if(i<nl) return Trip_no_j.Level(i-indexT_no_j).Get_n_down();

  cout << " Gas_of_HeI_Atoms:: fail " << endl;
  return 1;
}

const Transition_Data_HeI_A& Gas_of_HeI_Atoms::Get_Trans_Data(int i, int m) const
{
  if(i<indexT) return Sing.Level(i).Get_Trans_Data(m);
  else if(i<indexT_no_j) return Trip.Level(i-indexT).Get_Trans_Data(m);
  else if(i<nl) return Trip_no_j.Level(i-indexT_no_j).Get_Trans_Data(m);

  cout << " Gas_of_HeI_Atoms:: fail " << endl;
  return Sing.Level(0).Get_Trans_Data(0);
}

//===================================================================================
double Gas_of_HeI_Atoms::Xi(int i) const
{
  if(i<indexT) return Sing.Level(i).Get_Xi();
  else if(i<indexT_no_j) return Trip.Level(i-indexT).Get_Xi();
  else if(i<nl) return Trip_no_j.Level(i-indexT_no_j).Get_Xi();

  cout << " Gas_of_HeI_Atoms:: fail " << endl;
  return 0.0;
}

double Gas_of_HeI_Atoms::Ric(int i) const
{
  if(i<indexT) return Sing.Level(i).Get_Ric();
  else if(i<indexT_no_j) return Trip.Level(i-indexT).Get_Ric();
  else if(i<nl) return Trip_no_j.Level(i-indexT_no_j).Get_Ric();

  cout << " Gas_of_HeI_Atoms:: fail " << endl;
  return 0.0;
}

void Gas_of_HeI_Atoms::Set_Xi(int i, double X)
{
  if(i<indexT) Sing.Level(i).Set_Xi(X);
  else if(i<indexT_no_j) Trip.Level(i-indexT).Set_Xi(X);
  else if(i<nl) Trip_no_j.Level(i-indexT_no_j).Set_Xi(X);

  return;
}

void Gas_of_HeI_Atoms::Set_Ric(int i, double X)
{
  if(i<indexT) Sing.Level(i).Set_Ric(X);
  else if(i<indexT_no_j) Trip.Level(i-indexT).Set_Ric(X);
  else if(i<nl) Trip_no_j.Level(i-indexT_no_j).Set_Ric(X);

  return;
}

double Gas_of_HeI_Atoms::X_tot() const
{
    double r=0.0;   
    // smallest terms first!
    for(int i=nl-1; i>=indexT_no_j; i--) r+=Trip_no_j.Level(i-indexT_no_j).Get_Xi();
    for(int i=indexT_no_j-1; i>=indexT; i--) r+=Trip.Level(i-indexT).Get_Xi();
    // carefull here! i=-1 is non-sense; need signed int here!
    for(int i=indexT-1; i>=0; i--) r+=Sing.Level(i).Get_Xi();
    return r; 
}

void Gas_of_HeI_Atoms::Set_xxx_Transition_Data(int i, int m, int d, double v)
{
  if(i<indexT) Sing.Level(i).Set_xxx_Transition_Data(m, d, v);
  else if(i<indexT_no_j) Trip.Level(i-indexT).Set_xxx_Transition_Data(m, d, v);
  else if(i<nl) Trip_no_j.Level(i-indexT_no_j).Set_xxx_Transition_Data(m, d, v);
  return;
}

//===================================================================================
// Saha-relations with continuum
//===================================================================================
double Gas_of_HeI_Atoms::Ni_NeNc_LTE(int i, double TM) const
{
  if(i<indexT) return Sing.Level(i).Ni_NeNc_LTE(TM);
  else if(i<indexT_no_j) return Trip.Level(i-indexT).Ni_NeNc_LTE(TM);
  else if(i<nl) return Trip_no_j.Level(i-indexT_no_j).Ni_NeNc_LTE(TM);

  cout << " Gas_of_HeI_Atoms:: fail " << endl;
  return 0.0;
}

double Gas_of_HeI_Atoms::Xi_Saha(int i, double Xe, double Xc, double NH, double TM)  const
{
  if(i<indexT) return Sing.Level(i).Xi_Saha(Xe, Xc, NH, TM);
  else if(i<indexT_no_j) return Trip.Level(i-indexT).Xi_Saha(Xe, Xc, NH, TM);
  else if(i<nl) return Trip_no_j.Level(i-indexT_no_j).Xi_Saha(Xe, Xc, NH, TM);

  cout << " Gas_of_HeI_Atoms:: fail " << endl;
  return 0.0;
}

double Gas_of_HeI_Atoms::Ni_Saha(int i, double Ne, double Nc, double TM)  const
{
  if(i<indexT) return Sing.Level(i).Ni_Saha(Ne, Nc, TM);
  else if(i<indexT_no_j) return Trip.Level(i-indexT).Ni_Saha(Ne, Nc, TM);
  else if(i<nl) return Trip_no_j.Level(i-indexT_no_j).Ni_Saha(Ne, Nc, TM);

  cout << " Gas_of_HeI_Atoms:: fail " << endl;
  return 0.0;
}
 
//===================================================================================
//
// photoionization rate setup
//
//===================================================================================
// This is the structure for the rates (T==TopBase, S==Smits, H==Hydrogenic)
//
// HeI-Singlet:
//
//    "T",
//    "T", "T",
//    "T", "T", "T",
//    "T", "T", "T", "S",
//    "T", "T", "T", "S", "S", // n=5
//    "T", "T", "T", "H", "H", "H",
//    "T", "T", "T", "H", "H", "H", "H",
//    "T", "T", "T", "H", "H", "H", "H", "H",
//    "T", "T", "T", "H", "H", "H", "H", "H", "H",
//    "T", "H", "H", "H", "H", "H", "H", "H", "H", "H" // n=10
//
// HeI-Triplet:
//
//    "T", "T",
//    "T", "T", "T",
//    "T", "S", "T", "S",
//    "T", "S", "T", "S", "S", // n=5
//    "T", "H", "T", "H", "H", "H",
//    "T", "H", "T", "H", "H", "H", "H",
//    "T", "H", "T", "H", "H", "H", "H", "H",
//    "T", "H", "T", "H", "H", "H", "H", "H", "H",
//    "T", "H", "H", "H", "H", "H", "H", "H", "H", "H" // n=10
//
//===================================================================================
void Gas_of_HeI_Atoms::clear_Interaction_w_photons()
{
    for(int k=0; k<(int)Interaction_with_Photons_SH_QSP.size(); k++)
        for(int l=0; l<(int)Interaction_with_Photons_SH_QSP[k].size(); l++)
            Interaction_with_Photons_SH_QSP[k][l].clear();
    
    Interaction_with_Photons_SH_QSP.clear();
    
    return;
}

//===================================================================================
// initialize photoionization rates
//===================================================================================
void Gas_of_HeI_Atoms::init_photoionization_rates(int mflag)
{
    if(mflag>0) cout << " Gas_of_HeI_Atoms::setting up photoionization rates " << endl;

    vector<double> dum(8, 0.0);
    Ric_norm=Tg_norm_ref=dum;

    // setup memory
    clear_Interaction_w_photons();
    Interaction_with_Photons_SH_QSP.resize(Get_nShells());
    for(int n=1; n<=Get_nShells(); n++)
    {
        Interaction_with_Photons_SH_QSP[n-1].resize(n);
        
        for(int l=0; l<n; l++)
            Interaction_with_Photons_SH_QSP[n-1][l].init(n, l, 1.0, const_malpha_mp, mflag);
    }
    
    if(mflag>0) cout << " Gas_of_HeI_Atoms::init_photoionization_rates: done with hydrogenic " << endl;
    
    // topbase data
    load_all_Topbase_data(path+"TopBase_data/");
    
    if(mflag>0) cout << " Gas_of_HeI_Atoms::init_photoionization_rates: done with Topbase " << endl;
    
    return;
}

//===================================================================================
// use hydrogenic value Ric = (nuHe/nuH)^3 * RicH(T*nucH/nucHe)
//===================================================================================
double Gas_of_HeI_Atoms::R_ic(int i, double T_g)
{
    int n=Get_n(i), l=Get_l(i), s=Get_S(i);
    bool include_stim=1;
    double facBauman=(Get_J(i)!=-10 && switch_Bauman ? (2.0*Get_J(i)+1.0)/(2.0*l+1.0)/(2.0*s+1.0) : 1.0);
    
    if(n<=10)
    {
        //===========================================================================
        // Topbase
        //===========================================================================
        if(l==0) return facBauman*Ric_Topbase(n, l, s, T_g);
        else if(l==1)
        {
            if(s==0 && n<=9) return facBauman*Ric_Topbase(n, l, s, T_g);
            if(s==1 && n<=3) return facBauman*Ric_Topbase(n, l, s, T_g);
        }
        else if(l==2 && n<=9) return facBauman*Ric_Topbase(n, l, s, T_g);
        
        //===========================================================================
        // Smith recombination coefficients without stimulated recombination
        //===========================================================================
        if(s==0 && ((n==4 && l==3) ||
                    (n==5 && l==3) || (n==5 && l==4)) )
        {
            double Ric=Smits_Rec_Rate(n, l, s, T_g)/Ni_NeNc_LTE(i, T_g);
            if(include_stim) Ric*=1.0+1.0/(exp( min(700.0, const_h_kb*Get_nu_ion(i)/T_g))-1.0);
            
            return facBauman*Ric;
        }
        
        if(s==1 && ((n==4 && l==1) || (n==4 && l==3) ||
                    (n==5 && l==1) || (n==5 && l==3) || (n==5 && l==4)) )
        {
            double Ric=Smits_Rec_Rate(n, l, s, T_g)/Ni_NeNc_LTE(i, T_g);
            if(include_stim) Ric*=1.0+1.0/(exp( min(700.0, const_h_kb*Get_nu_ion(i)/T_g))-1.0);
            Ric*=Get_gw(i)/3.0/(2.0*l+1.0);  // factor from j average
            
            // setting like for Rubino-Martion et al paper
            //double Ric=Smits_Rec_Rate(n, l, s, T_g)/Trip.Level(n, l, l+1).Ni_NeNc_LTE(T_g);
            //if(include_stim) Ric*=1.0+1.0/(exp( min(700.0, const_h_kb*Trip.Level(n, l, l+1).Get_nu_ion()/T_g))-1.0);
            //Ric*=(2.0*(l+1.0)+1.0)/3.0/(2.0*l+1.0);
            
            return facBauman*Ric;
        }
    }
    
    //===============================================================================
    // Hydrogenic coefficients
    //===============================================================================
    double nucHe=Get_nu_ion(i), nucH=Interaction_with_Photons_SH_QSP[n-1][l].Get_nu_ionization();
    double chi=nucHe/nucH;

    return facBauman*pow(chi, 3)*Interaction_with_Photons_SH_QSP[n-1][l].R_nl_c_Int(T_g/chi);
}

double Gas_of_HeI_Atoms::R_ic(int n, int l, int s, int j, double T_g)
{ return R_ic(Get_Level_index(n, l, s, j), T_g); }

//===================================================================================
// use Saha relation
//===================================================================================
double Gas_of_HeI_Atoms::R_ci(int i, double T_g)
{ return R_ic(i, T_g)*Ni_NeNc_LTE(i, T_g); }

double Gas_of_HeI_Atoms::R_ci(int n, int l, int s, int j, double T_g)
{ return R_ci(Get_Level_index(n, l, s, j), T_g); }

//===================================================================================
// cross sections
//===================================================================================
double Gas_of_HeI_Atoms::sig_ic_Hyd(int i, double nu)
{
    int n=Get_n(i), l=Get_l(i);
    double facBauman=(Get_J(i)!=-10 && switch_Bauman ? (2.0*Get_J(i)+1.0)/(2.0*l+1.0)/(2.0*Get_S(i)+1.0) : 1.0);
    
    double nucHe=Get_nu_ion(i), nucH=Interaction_with_Photons_SH_QSP[n-1][l].Get_nu_ionization();
    double chi=nucHe/nucH;
    return facBauman*Interaction_with_Photons_SH_QSP[n-1][l].sig_phot_ion_lim(nu/chi);
}

double Gas_of_HeI_Atoms::sig_ic(int i, double nu, double Tg)
{
    //return sig_ic_Hyd(i, nu);
    
    int n=Get_n(i), l=Get_l(i), s=Get_S(i);
    double facBauman=(Get_J(i)!=-10 && switch_Bauman ? (2.0*Get_J(i)+1.0)/(2.0*l+1.0)/(2.0*s+1.0) : 1.0);
    
    if(n<=10)
    {
        //===========================================================================
        // Topbase
        //===========================================================================
        if(l==0) return facBauman*sig_ic_Topbase(n, l, s, nu);
        else if(l==1)
        {
            if(s==0 && n<=9) return facBauman*sig_ic_Topbase(n, l, s, nu);
            if(s==1 && n<=3) return facBauman*sig_ic_Topbase(n, l, s, nu);
        }
        else if(l==2 && n<=9) return facBauman*sig_ic_Topbase(n, l, s, nu);
        
        //===========================================================================
        // for Smith use hydrogenic shape but renormalized
        //===========================================================================
        if(s==0 && ((n==4 && l==3) ||
                    (n==5 && l==3) || (n==5 && l==4)) )
        {
            double nucHe=Get_nu_ion(i), nucH=Interaction_with_Photons_SH_QSP[n-1][l].Get_nu_ionization();
            if(nu<nucHe) return 0.0;
            double chi=nucHe/nucH;
            
            int index_norm;
            if(n==4 && l==3) index_norm=0;
            if(n==5 && l==3) index_norm=1;
            if(n==5 && l==4) index_norm=2;
            
            if(Tg!=Tg_norm_ref[index_norm])
            {
                double RicS=R_ic(i, Tg);
                double RicH=pow(chi, 3)*Interaction_with_Photons_SH_QSP[n-1][l].R_nl_c_Int(Tg/chi);

                Ric_norm[index_norm]=RicS/RicH;
                Tg_norm_ref[index_norm]=Tg;
            }

            return facBauman*Ric_norm[index_norm]*Interaction_with_Photons_SH_QSP[n-1][l].sig_phot_ion_lim(nu/chi);
        }
        
        if(s==1 && ((n==4 && l==1) || (n==4 && l==3) ||
                    (n==5 && l==1) || (n==5 && l==3) || (n==5 && l==4)) )
        {
            double nucHe=Get_nu_ion(i), nucH=Interaction_with_Photons_SH_QSP[n-1][l].Get_nu_ionization();
            if(nu<nucHe) return 0.0;
            double chi=nucHe/nucH;
            
            int index_norm;
            if(n==4 && l==1) index_norm=3;
            if(n==4 && l==3) index_norm=4;
            if(n==5 && l==1) index_norm=5;
            if(n==5 && l==3) index_norm=6;
            if(n==5 && l==4) index_norm=7;
            
            if(Tg!=Tg_norm_ref[index_norm])
            {
                double RicS=R_ic(i, Tg);
                double RicH=pow(chi, 3)*Interaction_with_Photons_SH_QSP[n-1][l].R_nl_c_Int(Tg/chi);
                
                Ric_norm[index_norm]=RicS/RicH;
                Tg_norm_ref[index_norm]=Tg;
            }
            
            return facBauman*Ric_norm[index_norm]*Interaction_with_Photons_SH_QSP[n-1][l].sig_phot_ion_lim(nu/chi);
        }
    }
    
    //===============================================================================
    // Hydrogenic value
    //===============================================================================
    double nucHe=Get_nu_ion(i), nucH=Interaction_with_Photons_SH_QSP[n-1][l].Get_nu_ionization();
    double chi=nucHe/nucH;
    
    return facBauman*Interaction_with_Photons_SH_QSP[n-1][l].sig_phot_ion_lim(nu/chi);
}

//===================================================================================
//===================================================================================
