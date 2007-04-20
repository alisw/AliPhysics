#ifndef ROOT_HydCommon
#define ROOT_HydCommon
//****************************************************************************//
//      -------------------------------------------------------------         //
//      HYDJET, fast MC code to simulate flow effects, jet production         // 
//      and jet quenching in heavy ion AA collisions at the LHC               //
//      -------------------------------------------------------------         //
//      This code is merging HYDRO (flow effects), PYTHIA6.4 (hard jet        // 
//      production) and PYQUEN (jet quenching)                                //
//      --------------------------------------------------------------        //
//                                                                            //
//      Igor Lokhtin, SINP MSU, Moscow, RU                                    //
//        e-mail: Igor.Lokhtin@cern.ch                                        // 
//                                                                            //
//      Reference for HYDJET:                                                 // 
//      I.P. Lokhtin, A.M. Snigirev,                                          // 
//      Eur. Phys. J. C 46 (2006) 211.                                        //
//                                                                            //
//      References for HYDRO:                                                 // 
//      N.A.Kruglov, I.P.Lokhtin, L.I.Sarycheva, A.M.Snigirev,                // 
//      Z. Phys. C 76 (1997) 99;                                              //  
//      I.P.Lokhtin, L.I.Sarycheva, A.M.Snigirev,                             // 
//      Phys. Lett. B 537 (2002) 261;                                         //   
//    I.P.Lokhtin, A.M.Snigirev, Preprint SINP MSU 2004-14/753,hep-ph/0312204.//
//                                                                            //
//      References for PYQUEN:                                                // 
//      I.P.Lokhtin, A.M.Snigirev, Eur.Phys.J. C16 (2000) 527;                //   
//    I.P.Lokhtin, A.M.Snigirev, Preprint SINP MSU 2004-13/752, hep-ph/0406038.//
//                                                                             //
//      References for PYTHIA:                                                 //
//      T.Sjostrand et al., Comput.Phys.Commun. 135 (2001) 238;                // 
//      T.Sjostrand, S. Mrena and P. Skands, hep-ph/0603175.                   //
//                                                                             //
//      Reference for JETSET event format:                                     //
//      T.Sjostrand, Comput.Phys.Commun. 82 (1994) 74.                         //
//                                                                             // 
//      --------------------------------------------------------------         //
//      Web-page:                                                              //
//      http://cern.ch/lokhtin/hydro                                           //
//      --------------------------------------------------------------         //
//                                                                             //
//**************************************************************************** //
#ifndef __CFORTRAN_LOADED
//*KEEP,cfortran.
#include "cfortran.h"
//*KEND.
#endif

extern "C" {


/*=========================================================*/
/* COMMON/HYFLOW/YTFL,YLFL,FPART */
/*---------------------------------------------------------*/
typedef struct {
   float    ytfl;
   float    ylfl;
   float    fpart;
} HyflowCommon;

#define HYFLOW COMMON_BLOCK(HYFLOW,hyflow)
COMMON_BLOCK_DEF(HyflowCommon,HYFLOW);

/**************************************************************************/
/*           D E S C R I P T I O N :                                      */
/*------------------------------------------------------------------------*/
/*COMMON /hyflow/ ytfl,ylfl,fpart
ytfl - maximum transverse collective rapidity, controls slope of low-pt spectra
(allowed range is 0.01<ytfl<3.0, default value is ytfl=1.).
ylfl - maximum longitudinal collective rapidity, controls width of eta-spectra
(allowed range is 0.01<ylfl<7.0, default value is ylfl=5.).
fpart - fraction of soft multiplicity proportional to the number of nucleon
participants; then (1.-fpart) will be the fraction of soft multiplicity
proportional to the number of nucleon-nucleon binary sub-collisions
(allowed range is 0.01<fpart<1.0, default value is fpart=1.).
========================================================================*/

/*========================================================================*/
/* COMMON/HYJPAR/NHSEL,PTMIN,NJET                                         */
/*------------------------------------------------------------------------*/
typedef struct {
   Int_t      nhsel;
   float      ptmin;
   Int_t      njet;
} HyjparCommon;

#define HYJPAR COMMON_BLOCK(HYJPAR,hyjpar)
COMMON_BLOCK_DEF(HyjparCommon,HYJPAR);
/*************************************************************************/
/*           D E S C R I P T I O N :                                     */
/*-----------------------------------------------------------------------*/
/*COMMON /hyjpar/ nhsel,ptmin,njet
Input Parameters:
nhsel - flag to include jet production in hydro event:
nhsel=0 - jet production off (pure HYDRO event);
nhsel=1 - jet production on, jet quenching off (HYDRO+njet*PYTHIA events);
nhsel=2 - jet production & jet quenching on (HYDRO+njet*PYQUEN events);
nhsel=3 - jet production on, jet quenching off, HYDRO off (njet*PYTHIA events);
nhsel=4 - jet production & jet quenching on, HYDRO off (njet*PYQUEN events);
(default value is nhsel = 0).
ptmin - minimal pt of parton-parton scattering in PYTHIA event, parameter
ckin(3) in PYTHIA common block /pysubs/ should be set equal to ptmin
(allowed range is 5 GeV < ptmin < 500 GeV, default value ptmin=10 GeV).

Output Parameters:
njet - number of hard parton-parton scatterings with pt>ptmin in event.
*/
/*=======================================================================*/



/*========================================================*/
/* COMMON/HYFPAR/ BGEN,NBCOL,NPART.NPYT,NHYD              */
/*--------------------------------------------------------*/
typedef struct {
   float  bgen;
   int  nbcol;
   int  npart;
   int  npyt;
   int  nhyd;
} HyfparCommon;

#define HYFPAR COMMON_BLOCK(HYFPAR,hyfpar)
COMMON_BLOCK_DEF(HyfparCommon,HYFPAR);
/*************************************************************************/
/*           D E S C R I P T I O N :                                     */
/*-----------------------------------------------------------------------*/
/*common /hyfpar/ bgen,nbcol,npart,npyt,nhyd
bgen - generated value of impact parameter in units of nucleus radius RA
(i.e the value in [fm] will be bgen*RA).
nbcol - mean number of nucleon-nucleon binary sub-collisions at given 'bgen'.
npart - mean number of nucleon participants at given 'bgen'.
npyt - multiplicity of hard PYTHIA/PYQUEN-induced particles in event
       (including full parton story).
nhyd - multiplicity of soft HYDRO-induced hadrons in event.
									 */
/*=======================================================================*/

/*=======================================================================*/
// COMMON/LUJETS/ N,K(150000,5),P(150000,5),V(150000,5)
/*-----------------------------------------------------------------------*/
typedef struct {
   Int_t    n;
   Int_t    k[5][150000];
   Float_t p[5][150000];
   Float_t v[5][150000];
} LujetsCommon;

#define LUJETS COMMON_BLOCK(LUJETS,lujets)
COMMON_BLOCK_DEF(LujetsCommon,LUJETS);
/*************************************************************************/
/*	     D E S C R I P T I O N :			                 */
/*-----------------------------------------------------------------------*/
/*COMMON /lujets/ n,k(150000,5),p(150000,5),v(150000,5)
n - total event multiplicity
k(i,1-5) - particle codes
p(i,1-5) - particle four-momentum and mass
v(i,1-5) - particle vertex, production time and lifetime

NOTE! First 'npyt' lines in event list correspond to PYTHIA/PYQUEN-induced
      particles, last 'nhyd' lines -- HYDRO-induced hadrons.
									 */
/*=======================================================================*/

/*=======================================================================*/
// COMMON/HYJETS/ NL,KL(150000,5),PL(150000,5),VL(150000,5)
/*-----------------------------------------------------------------------*/
typedef struct {
   Int_t     nl;
   Int_t     kl[5][150000];
   Float_t  pl[5][150000];
   Float_t  vl[5][150000];
} HyjetsCommon;

#define HYJETS COMMON_BLOCK(HYJETS,hyjets)
COMMON_BLOCK_DEF(HyjetsCommon,HYJETS);
/*************************************************************************/
/*           D E S C R I P T I O N :                                     */
/*-----------------------------------------------------------------------*/
/*COMMON /hyjets/ nl,kl(150000,5),pl(150000,5),vl(150000,5)
contains list of parton history of event in the same format as /lujets/ */
/*=======================================================================*/

/* COMMON from Pythia */

/*=======================================================================*/
/* COMMON/PYDAT1/ MSTU(200),PARU(200),MSTJ(200),PARJ(200)                */
/*-----------------------------------------------------------------------*/
typedef struct {
  Int_t     mstu[200];
  Double_t  paru[200];
  Int_t     mstj[200];
  Double_t  parj[200];
} Pydat1Common;

#define PYDAT1 COMMON_BLOCK(PYDAT1,pydat1)
COMMON_BLOCK_DEF(Pydat1Common,PYDAT1);


/*=======================================================================*/
// COMMON/PYSUBS/ MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
/*-----------------------------------------------------------------------*/
typedef struct {
  Int_t    msel;
  Int_t    mselpd;
  Int_t    msub[500];
  Int_t    kfin[2][81];
  Double_t ckin[200];
} PysubsCommon;

#define PYSUBS COMMON_BLOCK(PYSUBS,pysubs)
COMMON_BLOCK_DEF(PysubsCommon,PYSUBS);

/*=======================================================================*/
// COMMON/PYPARS/ MSTP(200),PARP(200),MSTI(200),PARI(200)
/*-----------------------------------------------------------------------*/
typedef struct {
  Int_t     mstp[200];
  Double_t  parp[200];
  Int_t     msti[200];
  Double_t  pari[200];
} PyparsCommon;

#define PYPARS COMMON_BLOCK(PYPARS,pypars)
COMMON_BLOCK_DEF(PyparsCommon,PYPARS);

}

#endif
