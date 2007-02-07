#ifndef AliAODPair_H
#define AliAODPair_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//_________________________________________________________________________
///////////////////////////////////////////////////////////////////////////
//
// class AliAODPair
//
// class implements pair of particles and taking care of caluclation (almost)
// all of pair properties (Qinv, InvMass,...)
// more info: http://aliweb.cern.ch/people/skowron/analyzer/index.html
//
////////////////////////////////////////////////////////////////////////////

#include <TMath.h>
#include <TObject.h>

#include "AliVAODParticle.h"

class AliAODPair: public TObject
{
 public:
   AliAODPair(Bool_t rev = kFALSE); //contructor
   AliAODPair(AliVAODParticle* part1, AliVAODParticle* part2, Bool_t rev = kFALSE); //contructor
   AliAODPair(const AliAODPair& in);
   
   virtual ~AliAODPair(){}
   
   AliAODPair& operator=(const AliAODPair& in);
   
   void SetParticles(AliVAODParticle* p1,AliVAODParticle* p2); //sets particles in the pair
   AliAODPair* GetSwappedPair() {return fSwappedPair;} //returns pair with swapped particles
   
   AliVAODParticle* Particle1() const {return fPart1;} //returns pointer to first particle
   AliVAODParticle* Particle2() const {return fPart2;} //returns pointer to decond particle
   
   virtual void     Changed();
   //Center Mass System - Longitudinally Comoving
   
   virtual Double_t GetInvMass(); //returns invariant mass of the pair
   virtual Double_t GetMt();
   virtual Double_t GetQInv(); //returns Q invariant
   virtual Double_t GetQSideLCMS(); //returns Q Side CMS longitudionally co-moving
   virtual Double_t GetQOutLCMS(); //returns Q out CMS longitudionally co-moving
   virtual Double_t GetQLongLCMS(); //returns Q Long CMS longitudionally co-moving
   virtual Double_t GetQtLCMS(); //returns Q transverse CMS longitudionally co-moving
   
   virtual Double_t GetQt(); //returns Q transverse to Kt
   
   
   virtual Double_t GetKt();  //returns K transverse
   virtual Double_t GetKStar();
   virtual Double_t GetKStarOut();  //z.ch.
   virtual Double_t GetKStarSide(); //z.ch.
   virtual Double_t GetKStarLong(); //z.ch.
   

   virtual Double_t GetAvarageDistance();//returns avarage distnace between two tracks
   
   virtual Double_t GetDeltaE(); //return difference of Energies
   virtual Double_t GetDeltaP(); //return difference of momenta (scalar difference)
   virtual Double_t GetDeltaPvector(); //return legth of difference vector of momenta
   virtual Double_t GetDeltaPt();
   virtual Double_t GetDeltaPx();
   virtual Double_t GetDeltaPy();
   virtual Double_t GetDeltaPz();
   
   virtual Double_t GetDeltaTheta();
   virtual Double_t GetDeltaPhi();
   
   virtual Double_t GetGammaToLCMS();
   virtual Double_t GetGammaToTransverse();
   virtual Double_t GetPIDProb() const {return fPart1->GetPidProb()*fPart2->GetPidProb();}
   
   virtual Double_t GetRStar() ;
   virtual Double_t GetR() ;//returns distance between particle production points   
   
   void   MirrorSecond();
   void   DeleteSecond();
   
   void   Print(const Option_t* option ) const {TObject::Print(option);}
   void   Print() ;
   
 protected:
   AliVAODParticle* fPart1;  //pointer to first particle
   AliVAODParticle* fPart2;  //pointer to second particle
  
   AliAODPair* fSwappedPair; //pointer to swapped pair
   
/************************************************************/
/************CMS (LC) Q's   *********************************/
/************************************************************/
   //Center Mass System - Longitudinally Comoving
   
   Double_t fQSideLCMS;  //value of Q side CMS longitudially co-moving
   Bool_t   fQSideLCMSNotCalc; //flag indicating if fQSideLCMS is already calculated for this pair
   
   Double_t fQOutLCMS; //value of Q out CMS longitudially co-moving
   Bool_t   fQOutLCMSNotCalc;//flag indicating if fQOutLCMS is already calculated for this pair
   
   Double_t fQLongLCMS; //value of Q long CMS longitudially co-moving
   Bool_t   fQLongLCMSNotCalc;//flag indicating if fQLongLCMS is already calculated for this pair
   
   Double_t fQtLCMS; //value of Qt CMS longitudially co-moving (hypot(qsidelcms,qoutlcms))
   Bool_t   fQtLCMSNotCalc;//flag indicating if fQLongLCMS is already calculated for this pair

   Double_t fQt; //value of Qt, projection of 3-mom diff to Kt
   Bool_t   fQtNotCalc;//flag indicating if fQt is already calculated for this pair
   
/************************************************************/
/************************************************************/
   Double_t fQInv;  //half of differnece of 4-momenta
   Bool_t   fQInvNotCalc;//flag indicating if fQInv is already calculated for this pair
   
   Double_t fInvMass;  //invariant mass
   Bool_t   fInvMassNotCalc;//flag indicating if fInvMass is already calculated for this pair
   
   Double_t fKt; //K == sum vector of particle's momenta. Kt transverse component
   Bool_t   fKtNotCalc;//flag indicating if fKt is already calculated for this pair
   
   Double_t fKStar; // KStar
   Bool_t   fKStarNotCalc;// flag indicating if fKStar is calculated
   Double_t fKStarOut; // KStarOut   z.ch.
   Double_t fKStarSide;// KStarSide  z.ch.
   Double_t fKStarLong;// KStarLong  z.ch.

   Bool_t   fKStarCompNotCalc; // flag indicating if CalcuteKStarComp() is calculated  z.ch.

   Double_t fPInv;  //invariant momentum
   
   Double_t fQSide; //Q Side
   Double_t fOut;//Q Out
   Double_t fQLong;//Q Long

   Double_t fMt;//Transverse coordinate of Inv. Mass
   Bool_t   fMtNotCalc;//flag indicating if Mt is calculated for current pair
      
   Double_t fInvMassSqr;//squre of invariant mass
   Bool_t   fMassSqrNotCalc; //flag indicating if fInvMassSqr for this pair
   void     CalculateInvMassSqr();
   
   Double_t fQInvL; //Qinv in longitudional direction
   Bool_t   fQInvLNotCalc;//flag indicating if fQInvL is calculated for current pair
   void     CalculateQInvL();

   Double_t fAvarageDistance;//value of the avarage distance calculated out of track points
   Bool_t   fAvarageDistanceNotCalc;//flag indicating if the avarage distance is calculated
   
   Double_t fPxSum;// Sum of Px momenta
   Double_t fPySum;// Sum of Py momenta
   Double_t fPzSum;// Sum of Pz momenta
   Double_t fESum;// Sum of energies
   Bool_t   fSumsNotCalc;//flag indicating if fPxSum,fPxSum,fPxSum and fESum is calculated for current pair
   void     CalculateSums();
   void     CalculateKStarComp();
   
   Double_t fPxDiff;// Difference of Px momenta
   Double_t fPyDiff;// Difference of Px momenta
   Double_t fPzDiff;// Difference of Px momenta
   Double_t fEDiff;// Difference of Px momenta
   Bool_t   fDiffsNotCalc;//flag indicating if fPxDiff,fPxDiff,fPxDiff and fEDiff is calculated for current pair
   void     CalculateDiffs();
   
   Double_t fGammaLCMS;//gamma of boost in LCMS
   Bool_t   fGammaLCMSNotCalc;//flag indicating if fGammaLCMS is calculated for current pair
   /***************************************************/
   Bool_t   fChanged;//flag indicating if object has been changed

   void     CalculateBase();
   Double_t AvDistance();
   
   
 private:
  ClassDef(AliAODPair,1)
};
/****************************************************************/
inline
void AliAODPair::SetParticles(AliVAODParticle* p1,AliVAODParticle* p2)
{
 //sets the particle to the pair
 
 fPart1 = p1; 
 fPart2 = p2;
 if (fSwappedPair) //if we have Swapped (so we are not)
   fSwappedPair->SetParticles(fPart2,p1); //set particles for him too
 Changed();
 //and do nothing until will be asked for
} 
/****************************************************************/

inline
void AliAODPair::Changed()
{
 // Resel all calculations (flags)
 fChanged           = kTRUE;
 fSumsNotCalc       = kTRUE;
 fDiffsNotCalc      = kTRUE;
 fMassSqrNotCalc    = kTRUE;
 fInvMassNotCalc    = kTRUE;
 fQInvNotCalc       = kTRUE;
 fMtNotCalc         = kTRUE;
 fQSideLCMSNotCalc = kTRUE;
 fQOutLCMSNotCalc  = kTRUE;
 fQLongLCMSNotCalc = kTRUE;
 fQtLCMSNotCalc    = kTRUE;
 fQtNotCalc        = kTRUE;
 fKtNotCalc         = kTRUE;
 fKStarNotCalc      = kTRUE;
 fKStarCompNotCalc  = kTRUE;
 fQInvLNotCalc      = kTRUE;
 fGammaLCMSNotCalc = kTRUE;
 fAvarageDistanceNotCalc = kTRUE;
}
/****************************************************************/
inline 
void AliAODPair::CalculateInvMassSqr()
 {
  //calculates square of qinv
  if (fMassSqrNotCalc)
   {
     CalculateSums();
 
     Double_t fPart12s= (fPxSum*fPxSum) + (fPySum*fPySum) + (fPzSum*fPzSum);
 
     fInvMassSqr=fESum*fESum-fPart12s;

     fMassSqrNotCalc = kFALSE;
   }
 }
/****************************************************************/
inline 
void AliAODPair::CalculateQInvL()
 {
 //Calculates square root of Qinv
  if (fQInvLNotCalc)
  {
   CalculateDiffs();
   fQInvL = fEDiff*fEDiff - ( fPxDiff*fPxDiff + fPyDiff*fPyDiff + fPzDiff*fPzDiff );
   fQInvLNotCalc = kFALSE;
  }
 }
/****************************************************************/ 
inline 
void AliAODPair::CalculateSums()
 {
   //calculates momenta and energy sums
   if(fSumsNotCalc)
    {
     fPxSum = fPart1->Px()+fPart2->Px();
     fPySum = fPart1->Py()+fPart2->Py();
     fPzSum = fPart1->Pz()+fPart2->Pz();
     fESum  = fPart1->E() + fPart2->E();
     fSumsNotCalc = kFALSE;
    }
 }
/****************************************************************/
inline
void AliAODPair::CalculateKStarComp()
{
  
  if (fKStarCompNotCalc)
    {
      CalculateSums();

      Double_t ptrans = fPxSum*fPxSum + fPySum*fPySum;
      Double_t mtrans = fESum*fESum - fPzSum*fPzSum;
      Double_t pinv  =  TMath::Sqrt(mtrans - ptrans);
      ptrans         =  TMath::Sqrt(ptrans);
      mtrans         =  TMath::Sqrt(mtrans);
      
      Double_t px1   = fPart1->Px();
      Double_t py1   = fPart1->Py();
      Double_t pz1   = fPart1->Pz();
      Double_t pE1   = fPart1->E();

      // boost to LCMS
      Double_t beta  = fPzSum / fESum;
      Double_t gamma = fESum / mtrans;

      fKStarLong     = gamma * (pz1 - beta * pE1);
      double   temp  = gamma * (pE1 - beta * pz1);

      // rotate in transverse plane
      fKStarSide = (-px1*fPySum + py1*fPxSum)/ptrans;
      fKStarOut  = ( px1*fPxSum + py1*fPySum)/ptrans;
 
      // go from LCMS to CMS
      gamma = mtrans/pinv;
      beta  = ptrans/mtrans;
      fKStarOut  = gamma * (fKStarOut - beta * temp);

      fKStarCompNotCalc = kFALSE;
    }
}

/****************************************************************/
inline 
void AliAODPair::CalculateDiffs()
 {
   //calculates momenta and energy differences 
   if(fDiffsNotCalc)
    {
     fPxDiff = fPart1->Px()-fPart2->Px();
     fPyDiff = fPart1->Py()-fPart2->Py();
     fPzDiff = fPart1->Pz()-fPart2->Pz();
     fEDiff  = fPart1->E() - fPart2->E();
     fDiffsNotCalc = kFALSE;
    }
 }

/****************************************************************/

inline 
Double_t AliAODPair::GetDeltaE() 
{
 //returns difference of energies
  return fPart1->E() - fPart2->E();
}
/****************************************************************/

inline 
Double_t AliAODPair::GetDeltaP() 
{
 //returns difference of momenta (scalars)
 return fPart1->P() - fPart2->P();
}
/****************************************************************/

inline 
Double_t AliAODPair::GetDeltaPvector() //return difference of momenta
{
 //returns legth of the momenta difference vector 
 CalculateDiffs();
 return TMath::Sqrt(fPxDiff*fPxDiff + fPyDiff*fPyDiff + fPzDiff*fPzDiff);
}
/****************************************************************/

inline 
Double_t AliAODPair::GetDeltaPt()
 {
   //returns difference of Pz
   return fPart1->Pt()-fPart2->Pt();
 }
/****************************************************************/

inline 
Double_t AliAODPair::GetDeltaPx()
 {
   //returns difference of Pz
   CalculateDiffs();
   return fPxDiff;
 }
/****************************************************************/
inline 
Double_t AliAODPair::GetDeltaPy()
 {
   //returns difference of Py
   CalculateDiffs();
   return fPyDiff;
 }

/****************************************************************/
inline 
Double_t AliAODPair::GetDeltaPz()
 {
   //returns difference of Pz
   CalculateDiffs();
   return fPzDiff;
 }
/****************************************************************/

inline 
Double_t AliAODPair::GetDeltaPhi()
 {
   //returns difference of Phi
   Double_t phi1 = fPart1->Phi();
   Double_t phi2 = fPart2->Phi();
   Double_t diff = phi1-phi2;
   if (TMath::Abs(diff) > TMath::Pi())
    {
      if (phi1 > TMath::Pi())
       {
         phi1-=TMath::TwoPi();
       }
      else
       {
         phi2-=TMath::TwoPi();
       }
      diff = phi1-phi2; 
    }
   return diff;
 }
/****************************************************************/

inline 
Double_t AliAODPair::GetDeltaTheta()
 {
   //returns difference of Theta
   return fPart1->Theta()-fPart2->Theta();
 }
/****************************************************************/


#endif
