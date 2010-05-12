// This header file declares structure objects for
// accessing common blocks from the fortran side of HYDJET++ (PYQUEN and PYTHIA).
// The input parameters and output information are transmitted 
// through these structures.

#ifndef HYJETCOMMON
#define HYJETCOMMON

extern "C" {

#define f2cFortran
#include "cfortran.h"

  //----------------------------------------------------------------
  // common /hyipar/ bminh,bmaxh,AW,RA,npar0,nbco0,Apb,Rpb,np,init,ipr        
  typedef struct // for HYIPAR common block
  {
    Double_t bminh;  // minimum impact parameter in units of nucleus radius RA   
    Double_t bmaxh;  // maximum impact parameter in units of nucleus radius RA
    Double_t AW;     // beam and target nucleus atomic weight
    Double_t RA;     // nucleus radius
    Double_t npar0;  // number of participants  at "reference point" (Pb,b=0) 
    Double_t nbco0;  // number of binary sub-collisions at "reference point" (Pb,b=0)
    Double_t Apb;    // Apb=207
    Double_t Rpb;    //Rpb=1.15*Apb^1/3
    Double_t np;     //internal variable of HYDJET
    Int_t init;      //internal flag to start HYDJET initialization
    Int_t ipr;      //internal flag to start HYDJET initialization
  } HYIPARCommon;
  
#define HYIPAR COMMON_BLOCK(HYIPAR,hyipar)
  COMMON_BLOCK_DEF(HYIPARCommon, HYIPAR);
  //----------------------------------------------------------------

  //      common/service/iseed_fromC,iPythDecay,parPYTH(100)
  typedef struct // for SERVICE common block
  {
    Int_t iseed_fromC;      // random number generator seed
    Int_t iPythDecay; //internal parameter for swith on/off PYTHIA decayer not used under ALIROOT 
    Double_t parPYTH[100]; //array for PYTHIA internal parameters
  } SERVICECommon;
 
#define SERVICE COMMON_BLOCK(SERVICE,service)
  COMMON_BLOCK_DEF(SERVICECommon, SERVICE);
  //----------------------------------------------------------------

  //  common/SERVICEEV/ipdg,delta
 // variables for program MYDELTA which allows to transmit delta -limits for BW from PYTHIA to c++ part
  typedef struct // for SERVICEEV common block
  {
    Float_t delta;    // limits for the resonance Breit-Wigner mass (from PYTHIA)
    Int_t KC;  //PYTHIA internal compressed code of particle        
    Int_t ipdg; //particle pdg code transmitted from from C++
  }SERVICEEVCommon;
 
#define SERVICEEV COMMON_BLOCK(SERVICEEV,serviceev)
  COMMON_BLOCK_DEF(SERVICEEVCommon, SERVICEEV);

  //----------------------------------------------------------------
  // common /hyjpar/ ptmin,sigin,sigjet,nhsel,ishad,njet 
  typedef struct // for HYJPAR common block
  {
    Double_t ptmin;   // minimum pt for N-N collision in PYTHIA 
    Double_t sigin;   // total inelastic NN cross section at given c.m.s. energy (in mb)
    Double_t sigjet;  // hard scattering NN cross section at given ptmin and energy (in mb).
    Int_t nhsel;      // nhsel parameter for the event type (0-hydro only; 1-hydro+unquenched jets; 
                      //                 2-hydro+quenched jets; 3-unquenched jets only; 4-quenched jets only)
    Int_t ishad;      // flag for turning on/off nuclear shadowing
    Int_t njet;       // number of hard parton-parton scatterings with pt>ptmin in event
  } HYJPARCommon;
 
#define HYJPAR COMMON_BLOCK(HYJPAR,hyjpar)
  COMMON_BLOCK_DEF(HYJPARCommon, HYJPAR);
  //----------------------------------------------------------------


  //      common /hypyin/ ene,rzta,rnta,bfix,ifb,nh
  typedef struct //for HYPYIN common block
  {
    Double_t ene; //CMS energy per nucleon [GeV]
    Double_t rzta;//fraction of protons in nucleus
    Double_t rnta;//fraction of neutrons in nucleus
    Double_t bfix; //fix impact parameter
    Int_t ifb; // flag of type of centrality generation 
               // =0 impact parameter is fixed (bfix)  
                // >0 or <0 impact parameter is generated with standard Glauber geometry 
                //              between minimum (bmin) and maximum (bmax) values; 
                              
    Int_t nh; //mean soft mult. in central PbPb
  }HYPYINCommon;
  
#define HYPYIN COMMON_BLOCK(HYPYIN,hypyin)
  COMMON_BLOCK_DEF(HYPYINCommon, HYPYIN);
  
  
  //----------------------------------------------------------------
  //  common /hyfpar/ bgen,nbcol,npart,npyt,nhyd,npart0        
  typedef struct // for HYFPAR common block
  {
    Double_t bgen;        // impact parameter 
    Int_t nbcol;          // number of binary collisions
    Double_t npart;       // n participants
    Double_t npart0;      // n participants
    Int_t npyt;           // number of pythia particles
    Int_t nhyd;           // number of hydjet particles
  }HYFPARCommon;
  
#define HYFPAR COMMON_BLOCK(HYFPAR,hyfpar)
  COMMON_BLOCK_DEF(HYFPARCommon, HYFPAR);
  
  //----------------------------------------------------------------
  typedef struct //for HYPART common block
  {
    Double_t ppart[50000][10];   // particle information
    Double_t bmin;               // min impact parameter
    Double_t bmax;               // max impact parameter
    Int_t njp;                   // 
  }HYPARTCommon;
  
#define HYPART COMMON_BLOCK(HYPART,hypart)
  COMMON_BLOCK_DEF(HYPARTCommon, HYPART);
  //----------------------------------------------------------------

  //      common /pyqpar/ T0,tau0,nf,ienglu,ianglu 
  typedef struct // for PYQPAR common block
  {
    Double_t T0;     // initial temperature of QGP for central Pb+Pb collisions  
    Double_t tau0;   // proper time of quark-gluon plasma formation 
    Int_t nf;        // number of flavours
    Int_t ienglu;    // type of energy loss
    Int_t ianglu;    // type of angular distribution for emitted gluons
  }PYQPARCommon;
  
#define PYQPAR COMMON_BLOCK(PYQPAR,pyqpar)
  COMMON_BLOCK_DEF(PYQPARCommon, PYQPAR);
  
} 
#endif                     
