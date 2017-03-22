#ifndef ROOT_PyqCommon
#define ROOT_PyqCommon

#ifndef __CFORTRAN_LOADED
//*KEEP,cfortran.
#include "cfortran.h"
//*KEND.
#endif

extern "C" {

//  common /pyqpar/ T0u,tau0u,nfu,ienglu,ianglu
    
typedef struct {
    Double_t t0;
    Double_t tau0;
    Int_t    nf;
    Int_t    iengl;
    Int_t    iangl;
} PyqparCommon;

#define PYQPAR COMMON_BLOCK(PYQPAR,pyqpar)
 COMMON_BLOCK_DEF(PyqparCommon,PYQPAR);
}


typedef struct {
    Double_t b1;
    Double_t psib1;
    Double_t rb1;
    Double_t rb2;
    Int_t    noquen;
} parimpCommon;
#define PARIMP COMMON_BLOCK(PARIMP,parimp)
COMMON_BLOCK_DEF(parimpCommon,PARIMP);

/*
Parameters in COMMON BLOCK PYQPAR can be varied by user: 


COMMON /pyqpar/ T0,tau0,nf,ienglu,ianglu

T0 - initial temparature of quark-gluon plasma  
(allowed range is 0.2 GeV < T0 < 2 GeV, default value is T0=1 GeV);

tau0 - proper time of quark-gluon plasma formation
(allowed range is 0.01 < tau0 < 10 fm/c, default value is tau0=0.1 fm/c)

nf - number of active quark flavours in quark-gluon plasma
(nf=0, 1, 2 or 3, default value is nf=0);

ienglu - flag to fix type of medium-induced partonic energy loss
(ienglu=0 - radiative and collisional loss,
ienglu=1 - radiative loss only, ienglu=2 - collisional loss only,
default value is ienglu=0);

ianglu - flag to fix type of angular distribution of emitted gluons
(ianglu=0 - small-angular, ianglu=1 - wide-angular, ianglu=2 - collinear,
default value is ianglu-0).  

NOTE! If specified by user value of such parameter is out of allowed range, 
the default value is used in PYQUEN run. 

NOTE! Default parameters of quark-gluon plasma (T0, tau0, nf) were selected as 
an estimation for LHC heavy ion beam energies. 
*/


#endif
