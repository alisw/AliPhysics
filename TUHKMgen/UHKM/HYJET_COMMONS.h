
#ifndef HYJETCOMMON
#define HYJETCOMMON

extern "C" {

#define f2cFortran
#include "cfortran.h"


//----------------------------------------------------------------
// common /hyipar/ bminh,bmaxh,AW,RA,npar0,nbco0,Apb,Rpb,np,init,ipr        
 typedef struct //HYIPAR
   {
 Double_t bminh;
 Double_t bmaxh; 
 Double_t AW;
 Double_t RA;
 Double_t npar0;
 Double_t nbco0;
 Double_t Apb;
 Double_t Rpb;
 Double_t np;
 Int_t init;
 Int_t ipr;
  }HYIPARCommon;
 
#define HYIPAR COMMON_BLOCK(HYIPAR,hyipar)
COMMON_BLOCK_DEF(HYIPARCommon, HYIPAR);
//----------------------------------------------------------------

//      common/service/iseed_fromC,iPythDecay,parPYTH(100)
 typedef struct //SERVICE
   {
 Int_t iseed_fromC; 
 Int_t iPythDecay;
 Double_t parPYTH[100];
  }SERVICECommon;
 
#define SERVICE COMMON_BLOCK(SERVICE,service)
COMMON_BLOCK_DEF(SERVICECommon, SERVICE);
//----------------------------------------------------------------

//  common/SERVICEEV/ipdg,delta

 typedef struct //SERVICEEV
   {
 Float_t delta;
 Int_t KC;
 Int_t ipdg;
  }SERVICEEVCommon;
 
#define SERVICEEV COMMON_BLOCK(SERVICEEV,serviceev)
COMMON_BLOCK_DEF(SERVICEEVCommon, SERVICEEV);

//----------------------------------------------------------------

  // common /hyjpar/ ptmin,sigin,sigjet,nhsel,ishad,njet 
 typedef struct //HYJPAR
   {
 Double_t ptmin;
 Double_t sigin;
 Double_t sigjet;
 Int_t nhsel;
 Int_t ishad;
 Int_t njet;
  }HYJPARCommon;
 
#define HYJPAR COMMON_BLOCK(HYJPAR,hyjpar)
COMMON_BLOCK_DEF(HYJPARCommon, HYJPAR);
//----------------------------------------------------------------


//      common /hypyin/ ene,rzta,rnta,bfix,ifb,nh
 typedef struct //HYPYIN
   {
 Double_t ene;
 Double_t rzta;
 Double_t rnta;
 Double_t bfix; 
 Int_t ifb;
 Int_t nh;
  }HYPYINCommon;
 
#define HYPYIN COMMON_BLOCK(HYPYIN,hypyin)
COMMON_BLOCK_DEF(HYPYINCommon, HYPYIN);


//----------------------------------------------------------------
 //  common /hyfpar/ bgen,nbcol,npart,npyt,nhyd,npart0        
 typedef struct //HYFPAR
   {
 Double_t bgen;
 Double_t nbcol;
 Double_t npart;
 Double_t npart0;
 Int_t npyt;
 Int_t nhyd;
  }HYFPARCommon;
 
#define HYFPAR COMMON_BLOCK(HYFPAR,hyfpar)
COMMON_BLOCK_DEF(HYFPARCommon, HYFPAR);

//----------------------------------------------------------------
 typedef struct //HYPART
   {
 Double_t ppart[50000][10];
 Double_t bmin;
 Double_t bmax;
 Int_t njp;
  }HYPARTCommon;
 
#define HYPART COMMON_BLOCK(HYPART,hypart)
COMMON_BLOCK_DEF(HYPARTCommon, HYPART);
//----------------------------------------------------------------

//      common /pyqpar/ T0,tau0,nf,ienglu,ianglu 

 typedef struct //PYQPAR
   {
 Double_t T0;
 Double_t tau0;
 Int_t nf;
 Int_t ienglu; 
 Int_t ianglu;
  }PYQPARCommon;
 
#define PYQPAR COMMON_BLOCK(PYQPAR,pyqpar)
COMMON_BLOCK_DEF(PYQPARCommon, PYQPAR);

} 
#endif                     
