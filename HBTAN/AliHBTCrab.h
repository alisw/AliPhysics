/* $Id$ */

// This class introduces the weight's calculation
// according to the Lednicky's algorithm.
// The detailed description of the algorithm can be found
// in comments to fortran code:
// fsiw.f, fsiini.f

#ifndef ALIHBTCrab_H
#define ALIHBTCrab_H

#include "AliHBTWeights.h"

//class Complex;
//typedef Complex double_complex;
#ifdef __DECCXX
#include <complex.h>
#endif
#include <math.h>

class AliHBTPair;

class AliHBTCrab: public AliHBTWeights
 {
   public:

     virtual ~AliHBTCrab(){fgCrab =0x0;}
     static AliHBTCrab* Instance();
     void Set();

     Double_t GetWeight(const AliHBTPair* partpair);
     void Init(Int_t pid1,Int_t pid2); //put the initial values in fortran commons fsiini, led_bldata
   private:
     AliHBTCrab();
     AliHBTCrab(const AliHBTCrab &/*source*/);
     AliHBTCrab & operator=(const AliHBTCrab& /*source*/);

     void get_com_quantities(const AliHBTPair* pair, double *qred,double *r,double *qdotr,double *mom, int *test);
     double  corrcalc(double trueqred,double trueqdotr,double truer);
     
     Bool_t fBreitWigner;
     Bool_t fReducedMom;
     Float_t fMaxMomentum;
     
     Bool_t  SetConfig(const AliHBTPair* pair);
     
     Int_t fPid1;
     Int_t fPid2;
     
     Double_t MASS1;
     Double_t MASS2;
 
     Float_t INTERACTION_WSYM;/* fractions of symmetric and antisym weights of the various spin channels */
     Float_t INTERACTION_WANTI;
     Float_t INTERACTION_WNOSYM;
     
     Float_t INTERACTION_DELK;
     Int_t INTERACTION_NKMAX;/* number of momentum points in mesh for strong/coul. interaction */
     
#ifdef __DECCXX
     static const complex ci;
#else
     static const double_complex ci;
#endif
     static const Double_t fgkROOT2;//! some const
     static const Double_t fgkWcons; //constant for fm->GeV conversion 1/0.1973
     
#ifdef __DECCXX
     complex cgamma(complex c);
#else
     double_complex cgamma(double_complex c);
#endif
     
     static AliHBTCrab* fgCrab;
     ClassDef(AliHBTCrab,1)
 };
 
#endif
