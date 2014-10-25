#ifndef ROOT_DCommon
#define ROOT_DCommon

#ifndef __CFORTRAN_LOADED
//*KEEP,cfortran.
#include "cfortran.h"
//*KEND.
#endif

/*
ccccc hepevt output
      integer nmxhep,kk
      parameter (nmxhep=4000)
      integer nevhep,nhep,isthep,idhep,jmohep,jdahep
      double precision phep,vhep
      common /hepevt/ nevhep,nhep,isthep(nmxhep),idhep(nmxhep),
     &jmohep(2,nmxhep),jdahep(2,nmxhep),phep(5,nmxhep),vhep(4,nmxhep)
*/
  const Int_t nmxhep = 4000;  
  typedef struct {
    Int_t nevhep;
    Int_t nhep;
    Int_t isthep[nmxhep];
    Int_t idhep[nmxhep];
    Int_t jmohep[nmxhep][2];
    Int_t jdahep[nmxhep][2];
    Double_t phep[nmxhep][5];
    Double_t vhep[nmxhep][4];
  } hepevtCommon;

#define HEPEVT COMMON_BLOCK(HEPEVT, hepevt)
COMMON_BLOCK_DEF(hepevtCommon, HEPEVT);

/*
ccccc Les Houches  Event Common Block
      INTEGER MAXNUP
      PARAMETER (MAXNUP=500)
      INTEGER NUP,IDPRUP,IDUP,ISTUP,MOTHUP,ICOLUP
      DOUBLE PRECISION XWGTUP,SCALUP,AQEDUP,AQCDUP,PUP,VTIMUP,SPINUP
      COMMON/HEPEUP/NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP,
     &              IDUP(MAXNUP),ISTUP(MAXNUP),MOTHUP(2,MAXNUP),
     &              ICOLUP(2,MAXNUP),PUP(5,MAXNUP),VTIMUP(MAXNUP),
     &              SPINUP(MAXNUP)
  */
 extern "C" {
  const Int_t MAXNUP = 500;
  typedef struct {
    Int_t NUP;
    Int_t IDPRUP;
    Double_t XWGTUP;
    Double_t SCALUP;
    Double_t AQEDUP;
    Double_t AQCDUP;
    Int_t IDUP[MAXNUP];
    Int_t ISTUP[MAXNUP];
    Int_t MOTHUP[MAXNUP][2];
    Int_t ICOLUP[MAXNUP][2];
    Double_t PUP[MAXNUP][5];
    Double_t VTIMUP[MAXNUP];
    Double_t SPINUP[MAXNUP];
  } hepeupCommon;

#define HEPEUP COMMON_BLOCK(HEPEUP, hepeup)
COMMON_BLOCK_DEF(hepeupCommon, HEPEUP);

  /*
    common/vars/s,rts,mmes,yx
  */
    typedef struct {
      Double_t s;
      Double_t rts;
      Double_t mmes;
      Double_t yx;
      Int_t    iin;
    } varsCommon;

#define VARS COMMON_BLOCK(VARS, vars)
  COMMON_BLOCK_DEF(varsCommon, VARS);


  //ccc 
  //      common/cuts/etaelmax,etaelmin,ptelmin,ptphmin,ecut,rmax,rmin,mcut
  typedef struct {
    Double_t etaelmax;
    Double_t etaelmin;
    Double_t ptelmin;
    Double_t ptphmin;
    Double_t ecut;
    Double_t rmax;
    Double_t rmin;
    Double_t mcut;
  } cutsCommon;
#define CUTS COMMON_BLOCK(CUTS, cuts)
 COMMON_BLOCK_DEF(cutsCommon, CUTS);

  /*
 character prefix*50,fsp*10,order*10,pflag*10,fsi*10,formf*10
     &,ppbar*10,output*10,mregge*10,cuts*10,unw*10
      common/flags/ iin, pflag, fsi, ppbar, output, cuts, unw
      common/ff/formf
  */
   typedef struct {
     char pflag[10];
     char fsi[10];
     char ppbar[10];
     char output[10];
     char cuts[10];
     char unw[10];
    } flagsCommon;

#define FLAGS COMMON_BLOCK(FLAGS, flags)
  COMMON_BLOCK_DEF(flagsCommon, FLAGS);

  typedef struct {
    char formf[10];
  } ffCommon;
#define FF COMMON_BLOCK(FF, ff)
  COMMON_BLOCK_DEF(ffCommon, FF);


#endif
}

