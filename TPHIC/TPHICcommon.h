#ifndef ROOT_TPHICcommon
#define ROOT_TPHICcommon
//------------------------------------------------------------------------
// TPHICcommon is an interface COMMON blocks of the fortran event generator
// of two-photon processes in ultraperipheral ion collisions
//%
// Yuri.Kharlov@cern.ch
// 15 April 2003
//------------------------------------------------------------------------

#ifndef __CFORTRAN_LOADED
#include "cfortran.h"
#endif

//       COMMON /ggini/ iproc,nevent,ilumf,lumfil,ebmn,eb,iz,ia,amas,
//      &               amin,amax,ymin,ymax,nmas,ny, kferm,
//      &               kf_onium,xmres,xgtres,xggres, xlumint, moddcy,
//      &               thetamin, costhv1, kv1,kv2,gvpar(4)
//       CHARACTER lumfil*80
extern "C" {
  typedef struct {
    Int_t   iproc;
    Int_t   nevent;
    Int_t   ilumf;
    char    lumfil[80];
    Float_t ebmn;
    Float_t eb;
    Int_t   iz;
    Int_t   ia;
    Float_t amas;
    Float_t amin;
    Float_t amax;
    Float_t ymin;
    Float_t ymax;
    Int_t   nmas;
    Int_t   ny;
    Int_t   kferm;
    Int_t   kfonium;
    Float_t xmres;
    Float_t xgtres;
    Float_t xggres;
    Float_t xlumint;
    Int_t   moddcy;
    Float_t thetamin;
    Float_t costhv1;
    Int_t   kv1;
    Int_t   kv2;
    Float_t qvpar[4];
  } GGiniCommon;
#ifdef IN_TPHICGEN_CXX
#define GGINI  COMMON_BLOCK(GGINI,ggini)
  COMMON_BLOCK_DEF(GGiniCommon,GGINI);
#endif

//       COMMON /ggevnt/ nrun,ievent,wsq,ygg,xmg1,xmg2, p2g(5),
//      &                ptag1(4),ptag2(4), ngg, kgg(10),pgg(20,5)
  typedef struct {
    Int_t   nrun;
    Int_t   ievent;
    Float_t wsq;
    Float_t ygg;
    Float_t xmg1;
    Float_t xmg2;
    Float_t p2g[5];
    Float_t ptag1[4];
    Float_t ptag2[4];
    Int_t   ngg;
    Int_t   kgg[10];
    Float_t pgg[5][20];
  } GGevntCommon;
#define GGEVNT COMMON_BLOCK(GGEVNT,ggevnt)
#ifdef IN_TPHICGEN_CXX
  COMMON_BLOCK_DEF(GGevntCommon,GGEVNT);
#endif

//       COMMON /ggxs/ xsmax0, xscur0, xscur, xsbra, xssum, ntry, xstot,
//      &  xstote, ssbr(10)
  typedef struct {
    Float_t xsmax0;
    Float_t xscur0;
    Float_t xscur;
    Float_t xsbra;
    Float_t xssum;
    Int_t   ntry;
    Float_t xstot;
    Float_t xstote;
    Float_t ssbr[10];
  } GGxsCommon;
#ifdef IN_TPHICGEN_CXX
#define GGXS   COMMON_BLOCK(GGXS,ggxs)
  COMMON_BLOCK_DEF(GGxsCommon,GGXS);
#endif

}

#endif
