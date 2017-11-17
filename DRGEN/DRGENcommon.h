#ifndef ROOT_DRGENcommon
#define ROOT_DRGENcommon
//------------------------------------------------------------------------
// DRGENcommon is an interface COMMON blocks of the fortran event generator
// of two-pomeron processes in pp collisions
//%
// Sergey.Evdokimov@cern.ch
// 09 Oct 2013
//------------------------------------------------------------------------

#ifndef __CFORTRAN_LOADED
#include "cfortran.h"
#endif

//COMMON /pp2init/ Iproc1, Iproc2, sqrtS, aMmin, aMmax, Weight, ProcXsec(20) 
extern "C" {
  typedef struct {
    Int_t   Iproc1;
    Int_t   Iproc2;
    Double_t sqrtS;
    Double_t aMmin;
    Double_t aMmax;
    Double_t Weight;
    Double_t ProcXsec[20];
    Double_t F2Polarization[8];
  } pp2initCommon;
#ifdef IN_TDRGEN_CXX
#define PP2INIT  COMMON_BLOCK(PP2INIT,pp2init)
  COMMON_BLOCK_DEF(pp2initCommon,PP2INIT);
#endif

  //      COMMON /pp2evnt/ ievent,Iproc,Npart,iParticleNumber,pParticle(10,5),iParticleCode(10),iParticleStatus(10) 
  typedef struct {
    Int_t   ievent;
    Int_t   Iproc;
    Int_t   Npart;
    Int_t   iParticleNumber;
    Double_t pParticle[10][5];
    Int_t iParticleCode[10];
    Int_t iParticleStatus[10];
  } pp2evntCommon;
#define PP2EVNT COMMON_BLOCK(PP2EVNT,pp2evnt)
#ifdef IN_TDRGEN_CXX
  COMMON_BLOCK_DEF(pp2evntCommon,PP2EVNT);
#endif

}

#endif //ROOT_DRGENcommon
