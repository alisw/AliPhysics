/* $Id$ */

//_____________________________________________________________________________
///////////////////////////////////////////////////////////////////////////////
//
// class AliHBTWeightTheorQInvFctn
//
// This function allows to obtain Q_inv correlation function with weights
// calculated by Lednicky's alghorithm.
// Numerator is filled with weighted events. Weights are attributed to simulated particles.
// Weights are calculated with corresponding simulated particles momenta.
// Denominator is filled with mixing unweighted simulated particles.
// One needs only simulated pairs, so 
// this function is of class AliHBTOnePairFctn1D.
//-----------------------------------------------------------
// This class introduces the weights calculated according 
// with functions of efficiency of identification (TPC+TOF) 
// (calculated by B.V. Batyunia).
// Author: Ludmila Malinina, JINR (malinina@sunhe.jinr.ru)
//-----------------------------------------------------------

#include "AliHBTWeightTheorFctn.h"

//--for test--AliHBTWeightQInvFctn* yyy= new AliHBTWeightQInvFctn();
/*************************************************************/
/*************************************************************/
/*************************************************************/

ClassImp(AliHBTWeightTheorQInvFctn)  
/*************************************************************/

AliHBTWeightTheorQInvFctn::AliHBTWeightTheorQInvFctn(Int_t nbins, Double_t maxXval, Double_t minXval):
  AliHBTOnePairFctn1D(nbins,maxXval,minXval)
{
  //ctor
 fWriteNumAndDen = kTRUE;//change default behaviour
 Rename("wqinvtheorcf","Q_{inv} Weight Theoretical Correlation Function");
}
/**************************************************************/

void  AliHBTWeightTheorQInvFctn::ProcessSameEventParticles(AliHBTPair* partpair)
{
  //Processes Particles and tracks Same different event
  partpair  = CheckPair(partpair);
  if (partpair == 0x0) return;
  Double_t weight = partpair->GetWeight();
  if(TMath::Abs(weight)<=10.) fNumerator->Fill(partpair->GetQInv(),weight);
} 
/**************************************************************/

TH1* AliHBTWeightTheorQInvFctn::GetResult() 
{
 //returns the scaled ratio
 delete fRatio;
 fRatio = GetRatio(Scale());
 return fRatio;
}                    
                                                              
/*************************************************************/
/*************************************************************/
/*************************************************************/

ClassImp(AliHBTWeightTheorQOutFctn)  
/*************************************************************/

AliHBTWeightTheorQOutFctn::AliHBTWeightTheorQOutFctn(Int_t nbins, Double_t maxXval, Double_t minXval):
  AliHBTOnePairFctn1D(nbins,maxXval,minXval)
{
  //ctor
 fWriteNumAndDen = kTRUE;//change default behaviour
 Rename("wqouttheorcf","Q_{out} Weight Theoretical Correlation Function");
}
/****************************************************************/

void  AliHBTWeightTheorQOutFctn::ProcessSameEventParticles(AliHBTPair* partpair)
{
  //Processes Particles and tracks Same different even
  partpair  = CheckPair(partpair);
  if (partpair == 0x0) return;
  Double_t weight = partpair->GetWeight();
  if(TMath::Abs(weight)<=10.) fNumerator->Fill(partpair->GetQOutLCMS(),weight);
} 

/**************************************************************/
TH1* AliHBTWeightTheorQOutFctn::GetResult() 
{
 //returns the scaled ratio
 delete fRatio;
 fRatio = GetRatio(Scale());
 return fRatio;
}                    

/*************************************************************/
/*************************************************************/
/*************************************************************/

ClassImp(AliHBTWeightTheorQSideFctn)  
/*************************************************************/

AliHBTWeightTheorQSideFctn::AliHBTWeightTheorQSideFctn(Int_t nbins, Double_t maxXval, Double_t minXval):
  AliHBTOnePairFctn1D(nbins,maxXval,minXval)
{
  //ctor
 fWriteNumAndDen = kTRUE;//change default behaviour
 Rename("wqsidetheorcf","Q_{side} Weight Theoretical Correlation Function");
}
/****************************************************************/

void  AliHBTWeightTheorQSideFctn::ProcessSameEventParticles(AliHBTPair* partpair)
{
  //Processes Particles and tracks Same different even
  partpair  = CheckPair(partpair);
  if (partpair == 0x0) return;
  Double_t weight = partpair->GetWeight();
  if(TMath::Abs(weight)<=10.) fNumerator->Fill(partpair->GetQSideLCMS(),weight);
} 
/**************************************************************/

TH1* AliHBTWeightTheorQSideFctn::GetResult() 
{
 //returns the scaled ratio
 delete fRatio;
 fRatio = GetRatio(Scale());
 return fRatio;
}                    

/*************************************************************/
/*************************************************************/
/*************************************************************/

ClassImp(AliHBTWeightTheorQLongFctn)  
/*************************************************************/

AliHBTWeightTheorQLongFctn::AliHBTWeightTheorQLongFctn(Int_t nbins, Double_t maxXval, Double_t minXval):
  AliHBTOnePairFctn1D(nbins,maxXval,minXval)
{
  //ctor
 fWriteNumAndDen = kTRUE;//change default behaviour
 Rename("wqlongtheorcf","Q_{long} Weight Theoretical Correlation Function");
}
/****************************************************************/

void  AliHBTWeightTheorQLongFctn::ProcessSameEventParticles(AliHBTPair* partpair)
{
  //Processes Particles and tracks Same different even
  partpair  = CheckPair(partpair);
  if (partpair == 0x0) return;
  Double_t weight = partpair->GetWeight();
  if(TMath::Abs(weight)<=10.) fNumerator->Fill(partpair->GetQLongLCMS(),weight);
} 
/**************************************************************/

TH1* AliHBTWeightTheorQLongFctn::GetResult() 
{
 //returns the scaled ratio
 delete fRatio;
 fRatio = GetRatio(Scale());
 return fRatio;
}                    

/*************************************************************/
/*************************************************************/
/*************************************************************/

ClassImp(AliHBTWeightTheorOSLFctn)

AliHBTWeightTheorOSLFctn::AliHBTWeightTheorOSLFctn(Int_t nXbins, Double_t maxXval, Double_t minXval,
                                                   Int_t nYbins, Double_t maxYval, Double_t minYval,
                                                   Int_t nZbins, Double_t maxZval, Double_t minZval):
 AliHBTOnePairFctn3D(nXbins,maxXval,minXval,nYbins,maxYval,minYval,nZbins,maxZval,minZval)
{
  //ctor
  fWriteNumAndDen = kTRUE;//change default behaviour
  Rename("wqosltheorcf","Q_{out}-Q_{side}-Q_{long} Weight Theoretical Correlation Fctn");
}

/*************************************************************/

void AliHBTWeightTheorOSLFctn::ProcessSameEventParticles(AliHBTPair* partpair)
{
//Fills numerator
  partpair  = CheckPair(partpair);
  if (partpair == 0x0) return;
  Double_t weight = partpair->GetWeight();
  Double_t out = TMath::Abs(partpair->GetQOutLCMS());
  Double_t side = TMath::Abs(partpair->GetQSideLCMS());
  Double_t lon = TMath::Abs(partpair->GetQLongLCMS());

/*  
  if (out < 0.01)
    if (side < 0.01)
      if (lon < 0.01)
       {
         Info("TheorOSL Num","================================================================");
         Info("TheorOSL Num","o:%f, s:%f, l:%f, w%f",out,side,lon,weight);
         Info("TheorOSL Num","First");
         partpair->Particle1()->Print();
         Info("TheorOSL Num","Second");
         partpair->Particle2()->Print();
         fflush(0);
       }
*/ 
  fNumerator->Fill(out,side,lon,weight);
}
/*************************************************************/

TH1* AliHBTWeightTheorOSLFctn::GetResult()
{
 //returns the scaled ratio
 delete fRatio;
 fRatio = GetRatio(Scale());
 return fRatio;
}                    
