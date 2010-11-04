#ifndef ROOT_ACommon
#define ROOT_ACommon

#ifndef __CFORTRAN_LOADED
//*KEEP,cfortran.
#include "cfortran.h"
//*KEND.
#endif

#define _MAXNPARTICLE_ 150001

extern "C" {

typedef struct {
   Float_t    hipr1[100];
   Int_t      ihpr2[50];
   Float_t    hint1[100];
   Int_t      ihnt2[50];
} HparntCommon;
#define HPARNT COMMON_BLOCK(HPARNT,hparnt)
COMMON_BLOCK_DEF(HparntCommon,HPARNT);

typedef struct {
  Int_t      natt;
  Float_t    eatt;
  Int_t      jatt;
  Int_t      nt;
  Int_t      np;
  Int_t      n0;
  Int_t      n01;
  Int_t      n10;
  Int_t      n11;
} Hmain1Common;
#define HMAIN1 COMMON_BLOCK(HMAIN1,hmain1)
COMMON_BLOCK_DEF(Hmain1Common,HMAIN1);

typedef struct {
  Float_t    bb;   
} ExtraCentCommon;
#define EXTRACENT COMMON_BLOCK(EXTRACENT,extracent)
COMMON_BLOCK_DEF(ExtraCentCommon,EXTRACENT);

typedef struct {
   Int_t    katt[4][_MAXNPARTICLE_];
   Float_t  patt[4][_MAXNPARTICLE_];
   Float_t  vatt[4][_MAXNPARTICLE_];    
} Hmain2Common;
#define HMAIN2 COMMON_BLOCK(HMAIN2,hmain2)
COMMON_BLOCK_DEF(Hmain2Common,HMAIN2);

typedef struct {
   Int_t    npj[300];
   Int_t    kfpj[500][300];
   Float_t  pjpx[500][300];
   Float_t  pjpy[500][300];
   Float_t  pjpz[500][300];
   Float_t  pjpe[500][300];
   Float_t  pjpm[500][300];
   Int_t    ntj[300];
   Int_t    kftj[500][300];
   Float_t  pjtx[500][300];
   Float_t  pjty[500][300];
   Float_t  pjtz[500][300];
   Float_t  pjte[500][300];
   Float_t  pjtm[500][300];
} Hjjet1Common;
#define HJJET1 COMMON_BLOCK(HJJET1,hjjet1)
COMMON_BLOCK_DEF(Hjjet1Common,HJJET1);

typedef struct {
   Int_t      nsg;
   Int_t      njsg[_MAXNPARTICLE_];
   Int_t      iasg[3][_MAXNPARTICLE_];
   Int_t      k1sg[100][_MAXNPARTICLE_];
   Int_t      k2sg[100][_MAXNPARTICLE_];
   Float_t    pxsg[100][_MAXNPARTICLE_];
   Float_t    pysg[100][_MAXNPARTICLE_];
   Float_t    pzsg[100][_MAXNPARTICLE_];
   Float_t    pesg[100][_MAXNPARTICLE_];
   Float_t    pmsg[100][_MAXNPARTICLE_];
} Hjjet2Common;
#define HJJET2 COMMON_BLOCK(HJJET2,hjjet2)
COMMON_BLOCK_DEF(Hjjet2Common,HJJET2);

typedef struct {
   Int_t    nfp[15][300];
   Float_t  pp[15][300];
   Int_t    nft[15][300];
   Float_t  pt[15][300];
} HstrngCommon;
#define HSTRNG COMMON_BLOCK(HSTRNG,hstrng)
COMMON_BLOCK_DEF(HstrngCommon,HSTRNG);

typedef struct {
   Float_t  yp[300][3];
   Float_t  yt[300][3];
} HjcrdnCommon;
#define HJCRDN COMMON_BLOCK(HJCRDN,hjcrdn)
COMMON_BLOCK_DEF(HjcrdnCommon,HJCRDN);

typedef struct {
  Int_t    mstu[200];
  Float_t  paru[200];
  Int_t    mstj[200];
  Float_t  parj[200];
} Ludat1Common;
#define LUDAT1 COMMON_BLOCK(LUDAT1A,ludat1a)
COMMON_BLOCK_DEF(Ludat1Common,LUDAT1);

typedef struct {
  Int_t   lblast[_MAXNPARTICLE_];
  Float_t xlast[_MAXNPARTICLE_][4];
  Float_t plast[_MAXNPARTICLE_][4];
  Int_t   nlast;
} HBTCommon;
#define HBT COMMON_BLOCK(HBT,hbt)
COMMON_BLOCK_DEF(HBTCommon,HBT);

typedef struct {
  Int_t nevent;
  Int_t isoft;
  Int_t isflag;
  Int_t izpc;
} AnimCommon;
#define ANIM COMMON_BLOCK(ANIM,anim)
COMMON_BLOCK_DEF(AnimCommon,ANIM);

typedef struct {
  Float_t xmp; 
  Float_t xmu; 
  Float_t alpha;
  Float_t rscut2;
  Float_t cutof2;
} Para2Common;
#define PARA2 COMMON_BLOCK(PARA2,para2)
COMMON_BLOCK_DEF(Para2Common,PARA2);

typedef struct {
  Float_t masspr;
  Float_t massta;
  Int_t iseed;
  Int_t iavoid;
  Float_t dt;
} Input1Common;
#define INPUT1 COMMON_BLOCK(INPUT1,input1)
COMMON_BLOCK_DEF(Input1Common,INPUT1);

typedef struct {
  Int_t ilab;
  Int_t manyb;
  Int_t ntmax;
  Int_t icoll;
  Int_t insys;
  Int_t ipot;
  Int_t mode; 
  Int_t imomen;
  Int_t nfreq;
  Int_t icflow;
  Int_t icrho;
  Int_t icou;
  Int_t kpoten;
  Int_t kmul; 
} Input2Common;
#define INPUT2 COMMON_BLOCK(INPUT2,input2)
COMMON_BLOCK_DEF(Input2Common,INPUT2);

typedef struct {
  Int_t ipop;
} PopcornCommon;
#define POPCORN COMMON_BLOCK(POPCORN,popcorn)
COMMON_BLOCK_DEF(PopcornCommon,POPCORN);

#undef _MAXNPARTICLE_ 

} /*extern*/
#endif

