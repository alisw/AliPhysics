/* $Id$ */

//__________________________________________________________________________
////////////////////////////////////////////////////////////////////////////
//
// class AliHBTCrab
//
// This class introduces the weight's calculation
// according to the Lednicky's algorithm.
// The detailed description of the algorithm can be found
// in comments to fortran code:
// fsiw.f, fsiini.f
//
// Piotr.Skowronski@cern.ch
////////////////////////////////////////////////////////////////////////////

#ifndef ALIHBTCrab_H
#define ALIHBTCrab_H

#include "AliHBTWeights.h"

#ifdef __DECCXX
 #include <complex.h>
#else
 class Complex;
 typedef Complex doublecomplex;
#endif
 
//#include <math.h>
 
class AliHBTPair;

class AliHBTCrab: public AliHBTWeights
 {
   public:

     AliHBTCrab();
     virtual ~AliHBTCrab(){fgCrab =0x0;}
     static AliHBTCrab* Instance();
     void Set();

     Double_t GetWeight(AliHBTPair* partpair);
     void Init(Int_t pid1,Int_t pid2); //put the initial values in fortran commons fsiini, led_bldata
     
   private:
     AliHBTCrab(const AliHBTCrab &/*source*/);
     AliHBTCrab & operator=(const AliHBTCrab& /*source*/);

     void GetComQuantities(const AliHBTPair* pair, double *qred,double *r,double *qdotr,double *mom, int *test);
     double  CorrCalc(double trueqred,double trueqdotr,double truer);
     
     Bool_t fBreitWigner;//switch if to calculated BW
     Bool_t fReducedMom;//switch if
     Float_t fMaxMomentum;//switch if
     
     Bool_t  SetConfig(const AliHBTPair* pair);
     
     Int_t fPid1;//PID of the first particle
     Int_t fPid2;//PID of the second particle
     
     Double_t fMass1;//mass of the first particle
     Double_t fMass2;//mass of the second particle
 
     Float_t fInteractionWsym;// fractions of symmetric and antisym weights of the various spin channels 
     Float_t fInteractionWanti;//comment
     Float_t fInteractionWnosym;//comment
     
     Float_t fInteractionDelk;//comment
     Int_t fInteractionNkmax;// number of momentum points in mesh for strong/coul. interaction 
     
#ifdef __DECCXX
     static const complex fgkCI;//complex (1,0)
#else
     static const doublecomplex fgkCI;//complex (1,0)
#endif
     static const Double_t fgkROOT2;//! some const
     static const Double_t fgkWcons; //constant for fm->GeV conversion 1/0.1973
     
#ifdef __DECCXX
     complex CGamma(complex c);
#else
     doublecomplex CGamma(doublecomplex c);
#endif
     
     static AliHBTCrab* fgCrab; //pointer to instance of this class - singleton
     ClassDef(AliHBTCrab,1)
 };
 
#endif
