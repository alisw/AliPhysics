#include "AliHBTWeightFctn.h"
/* $Id$ */
//_________________________________________________________________________
//
// class AliHBTWeightQInvFctn
//
// This class allows to obtain Q_inv correlation function with weights
// calculated by Lednicky's alghorithm.
// Numerator is filled with weighted events. Weights are attributed to reconstructed tracks.
// Weights are calculated with corresponding simulated particles momenta.
// Denominator is filled with mixing unweighted reconstructed tracks.
// One needs both pairs 
// (simulated and recontructed), thus function is of class AliHBTTwoPairFctn1D.
// Author: Ludmila Malinina, JINR (malinina@sunhe.jinr.ru)
//
////////////////////////////////////////////////////////////////////////////////

ClassImp( AliHBTWeightQInvFctn )

  
/****************************************************************/
AliHBTWeightQInvFctn::AliHBTWeightQInvFctn(Int_t nbins, Double_t maxXval, Double_t minXval):
 AliHBTTwoPairFctn1D(nbins,maxXval,minXval)
{
//ctor
 fWriteNumAndDen = kTRUE;//change default behaviour
 Rename("wqinvcf","Q_{inv} Weight Correlation Function");
}
/****************************************************************/

void  AliHBTWeightQInvFctn::ProcessSameEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair)
{
  //process particles from same events (fills numerator)
  trackpair = CheckPair(trackpair);
//  partpair  = CheckPair(partpair);
  if ( trackpair && partpair)     
  {
//    Double_t weightPID=1.;
    Double_t weight = 1.0;
    if ( ( trackpair->Particle1()->GetPdgCode() == partpair->Particle1()->GetPdgCode()) &&
         ( trackpair->Particle2()->GetPdgCode() == partpair->Particle2()->GetPdgCode())    )
      {   
         weight=partpair->GetWeight();
      }   
//    Double_t weight=weightHBT*weightPID;
    fNumerator->Fill(trackpair->GetQInv(),weight);
  }
} 
/****************************************************************/

void  AliHBTWeightQInvFctn::ProcessDiffEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair)
{
  // Fills the denominator using mixed pairs
  trackpair = CheckPair(trackpair);
//  partpair  = CheckPair(partpair);
  if ( trackpair && partpair)  
  {
     fDenominator->Fill(trackpair->GetQInv());
  }
}
/**************************************************************/

TH1* AliHBTWeightQInvFctn::GetResult()
{ 
//returns ratio of numerator and denominator                                    
 //returns the scaled ratio
 delete fRatio;
 fRatio = GetRatio(Scale());
 return fRatio;
}                    
                                                              
/**************************************************************************************/
/**************************************************************************************/
/**************************************************************************************/
/**************************************************************************************/

ClassImp(AliHBTWeightQOutFctn)
    
AliHBTWeightQOutFctn::AliHBTWeightQOutFctn(Int_t nbins, Double_t maxXval, Double_t minXval):
 AliHBTTwoPairFctn1D(nbins,maxXval,minXval)
{
//ctor
 fWriteNumAndDen = kTRUE;//change default behaviour
 Rename("wqoutcf","Q_{out} Weight Correlation Function");
}
/****************************************************************/

void AliHBTWeightQOutFctn::ProcessSameEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair)
{
  //process particles from same events (fills numerator)  
  trackpair = CheckPair(trackpair);
//  partpair  = CheckPair(partpair);
  if ( trackpair && partpair)     
  {
//    Double_t weightPID=1.;
    Double_t weight = 1.0;
    if ( ( trackpair->Particle1()->GetPdgCode() == partpair->Particle1()->GetPdgCode()) &&
         ( trackpair->Particle2()->GetPdgCode() == partpair->Particle2()->GetPdgCode())    )
      {   
         weight=partpair->GetWeight();
      }   
//    Double_t weight=weightHBT*weightPID;
    fNumerator->Fill(trackpair->GetQOutLCMS(),weight);
  }
} 
/****************************************************************/

void AliHBTWeightQOutFctn::ProcessDiffEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair)
{
  //process particles from diff events (fills denominator)
  trackpair = CheckPair(trackpair);
//  partpair  = CheckPair(partpair);
  if ( trackpair && partpair)  
  {
     fDenominator->Fill(trackpair->GetQOutLCMS());
  }
}
/**************************************************************/

TH1* AliHBTWeightQOutFctn::GetResult() 
{ 
 //returns the scaled ratio
 delete fRatio;
 fRatio = GetRatio(Scale());
 return fRatio;
}                    
                                                              
/*************************************************************************************/ 
/*************************************************************************************/ 
/*************************************************************************************/ 

ClassImp(AliHBTWeightQLongFctn)
AliHBTWeightQLongFctn::AliHBTWeightQLongFctn(Int_t nbins, Double_t maxXval, Double_t minXval):
 AliHBTTwoPairFctn1D(nbins,maxXval,minXval)
{
//ctor
 fWriteNumAndDen = kTRUE;//change default behaviour
 Rename("wqlongcf","Q_{long} Weight Correlation Function");
}
/****************************************************************/

void AliHBTWeightQLongFctn::ProcessSameEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair)
{
 //process particles from same events (fills numerator)
  trackpair = CheckPair(trackpair);
//  partpair  = CheckPair(partpair);
  if ( trackpair && partpair)     
  {
//    Double_t weightPID=1.;
    Double_t weight = 1.0;
    if ( ( trackpair->Particle1()->GetPdgCode() == partpair->Particle1()->GetPdgCode()) &&
         ( trackpair->Particle2()->GetPdgCode() == partpair->Particle2()->GetPdgCode())    )
      {   
         weight=partpair->GetWeight();
      }   
//    Double_t weight=weightHBT*weightPID;
    fNumerator->Fill(trackpair->GetQLongLCMS(),weight);
  }
} 
/****************************************************************/

void AliHBTWeightQLongFctn::ProcessDiffEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair)
{
  //process particles from diff events (fills denominator)
  trackpair = CheckPair(trackpair);
//  partpair  = CheckPair(partpair);
  if ( trackpair && partpair)  
  {
     fDenominator->Fill(trackpair->GetQLongLCMS());
  }
}
/**************************************************************/
TH1* AliHBTWeightQLongFctn::GetResult()
{ 
 //returns the scaled ratio
 delete fRatio;
 fRatio = GetRatio(Scale());
 return fRatio;
}                    
                                                              
/*************************************************************************************/ 
/*************************************************************************************/ 
/*************************************************************************************/ 

ClassImp(AliHBTWeightQSideFctn)
/*************************************************************************************/ 

AliHBTWeightQSideFctn::AliHBTWeightQSideFctn(Int_t nbins, Double_t maxXval, Double_t minXval):
 AliHBTTwoPairFctn1D(nbins,maxXval,minXval)
{
//ctor
 fWriteNumAndDen = kTRUE;//change default behaviour
 Rename("wqsidecf","Q_{side} Weight Correlation Function");
}
/****************************************************************/

void AliHBTWeightQSideFctn::ProcessSameEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair)
{
  //process particles from same events (fills numerator)
  trackpair = CheckPair(trackpair);
//  partpair  = CheckPair(partpair);
  if ( trackpair && partpair)     
  {
//    Double_t weightPID=1.;
    Double_t weight = 1.0;
    if ( ( trackpair->Particle1()->GetPdgCode() == partpair->Particle1()->GetPdgCode()) &&
         ( trackpair->Particle2()->GetPdgCode() == partpair->Particle2()->GetPdgCode())    )
      {   
         weight=partpair->GetWeight();
      }   
//    Double_t weight=weightHBT*weightPID;
    fNumerator->Fill(trackpair->GetQSideLCMS(),weight);
  }
} 
/****************************************************************/

void  AliHBTWeightQSideFctn::ProcessDiffEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair)
{
  //process particles from diff events (fills denominator)
  trackpair = CheckPair(trackpair);
//  partpair  = CheckPair(partpair);
  if ( trackpair && partpair)  
  {
     fDenominator->Fill(trackpair->GetQSideLCMS());
  }
}
/**************************************************************/

TH1* AliHBTWeightQSideFctn::GetResult() 
{ 
 //returns the scaled ratio
 delete fRatio;
 fRatio = GetRatio(Scale());
 return fRatio;
}                    
                                                              
/*************************************************************************************/ 
/*************************************************************************************/ 
/*************************************************************************************/ 

ClassImp(AliHBTWeightTwoKStarFctn)
/*************************************************************************************/ 
AliHBTWeightTwoKStarFctn::AliHBTWeightTwoKStarFctn(Int_t nbins, Double_t maxXval, Double_t minXval):
 AliHBTTwoPairFctn1D(nbins,maxXval,minXval)
{
//ctor
 fWriteNumAndDen = kTRUE;//change default behaviour
 Rename("wtwokstarcf","2*K^{*} Weight Correlation Function");
}
/****************************************************************/

void AliHBTWeightTwoKStarFctn::ProcessSameEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair)
{
  //process particles from same events (fills numerator)
  trackpair = CheckPair(trackpair);
//  partpair  = CheckPair(partpair);
  if ( trackpair && partpair)     
  {
//    Double_t weightPID=1.;
    Double_t weight = 1.0;
    if ( ( trackpair->Particle1()->GetPdgCode() == partpair->Particle1()->GetPdgCode()) &&
         ( trackpair->Particle2()->GetPdgCode() == partpair->Particle2()->GetPdgCode())    )
      {   
         weight=partpair->GetWeight();
      }   
//    Double_t weight=weightHBT*weightPID;
    fNumerator->Fill(2.0*(trackpair->GetKStar()),weight);
  }
} 
/****************************************************************/

void  AliHBTWeightTwoKStarFctn::ProcessDiffEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair)
{
  //process particles from diff events (fills denominator)
  trackpair = CheckPair(trackpair);
//  partpair  = CheckPair(partpair);
  if ( trackpair && partpair)  
  {
     fDenominator->Fill(2.0*(trackpair->GetKStar()));
  }
}
/**************************************************************/
TH1* AliHBTWeightTwoKStarFctn::GetResult() 
                                                                               
{ 
 //returns the scaled ratio
 delete fRatio;
 fRatio = GetRatio(Scale());
 return fRatio;
}                    
                                                              
/*************************************************************************************/ 
/*************************************************************************************/ 
/*************************************************************************************/ 

ClassImp(AliHBTWeightQOutQSideFctn)
/*************************************************************************************/ 
    
AliHBTWeightQOutQSideFctn::AliHBTWeightQOutQSideFctn(Int_t nxbins, Double_t maxXval, Double_t minXval,
                                                         Int_t nybins, Double_t maxYval, Double_t minYval):
 AliHBTTwoPairFctn2D(nxbins,maxXval,minXval,nybins,maxYval,minYval)
{
  //ctor
 fWriteNumAndDen = kTRUE;//change default behaviour
 Rename("wqoutqsidecf","Q_{out} Q_{side} Weight Correlation Function 2D");
}    
/*************************************************************************************/ 

void AliHBTWeightQOutQSideFctn::ProcessSameEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair)
{
  //process particles from same events (fills numerator)
  trackpair = CheckPair(trackpair);
//  partpair  = CheckPair(partpair);
  if ( trackpair && partpair)     
  {
//    Double_t weightPID=1.;
    Double_t weight = 1.0;
    if ( ( trackpair->Particle1()->GetPdgCode() == partpair->Particle1()->GetPdgCode()) &&
         ( trackpair->Particle2()->GetPdgCode() == partpair->Particle2()->GetPdgCode())    )
      {   
         weight=partpair->GetWeight();
      }   
//    Double_t weight=weightHBT*weightPID;
    fNumerator->Fill(trackpair->GetQOutLCMS(),trackpair->GetQSideLCMS(),weight);
  }
} 
/****************************************************************/

void AliHBTWeightQOutQSideFctn::ProcessDiffEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair)
{
  //process particles from diff events (fills denominator)
  trackpair = CheckPair(trackpair);
//  partpair  = CheckPair(partpair);
  if ( trackpair && partpair)  
  {
     fDenominator->Fill(trackpair->GetQOutLCMS(),trackpair->GetQSideLCMS());
  }
}
/**************************************************************/

TH1* AliHBTWeightQOutQSideFctn::GetResult()
{
 //returns the scaled ratio
 delete fRatio;
 fRatio = GetRatio(Scale());
 return fRatio;
}

/*************************************************************************************/ 
/*************************************************************************************/ 
/*************************************************************************************/ 

ClassImp(AliHBTWeightQOutQLongFctn)
/*************************************************************************************/ 
    
AliHBTWeightQOutQLongFctn::AliHBTWeightQOutQLongFctn(Int_t nxbins, Double_t maxXval, Double_t minXval,
                                                         Int_t nybins, Double_t maxYval, Double_t minYval):
 AliHBTTwoPairFctn2D(nxbins,maxXval,minXval,nybins,maxYval,minYval)
{
  //ctor
 fWriteNumAndDen = kTRUE;//change default behaviour
 Rename("wqoutqlongcf","Q_{out} Q_{long} Weight Correlation Function 2D");
}    
/*************************************************************************************/ 

void AliHBTWeightQOutQLongFctn::ProcessSameEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair)
{
  //process particles from same events (fills numerator)
  trackpair = CheckPair(trackpair);
//  partpair  = CheckPair(partpair);
  if ( trackpair && partpair)     
  {
//    Double_t weightPID=1.;
    Double_t weight = 1.0;
    if ( ( trackpair->Particle1()->GetPdgCode() == partpair->Particle1()->GetPdgCode()) &&
         ( trackpair->Particle2()->GetPdgCode() == partpair->Particle2()->GetPdgCode())    )
      {   
         weight=partpair->GetWeight();
      }   
//    Double_t weight=weightHBT*weightPID;
    fNumerator->Fill(trackpair->GetQOutLCMS(),trackpair->GetQLongLCMS(),weight);
  }
} 
/****************************************************************/

void AliHBTWeightQOutQLongFctn::ProcessDiffEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair)
{
  //process particles from diff events (fills denominator)
  trackpair = CheckPair(trackpair);
//  partpair  = CheckPair(partpair);
  if ( trackpair && partpair)  
  {
     fDenominator->Fill(trackpair->GetQOutLCMS(),trackpair->GetQLongLCMS());
  }
}
/**************************************************************/

TH1* AliHBTWeightQOutQLongFctn::GetResult()
{
 //returns the scaled ratio
 delete fRatio;
 fRatio = GetRatio(Scale());
 return fRatio;
}

/*************************************************************************************/ 
/*************************************************************************************/ 
/*************************************************************************************/ 

ClassImp(AliHBTWeightQSideQLongFctn)
/*************************************************************************************/ 
    
AliHBTWeightQSideQLongFctn::AliHBTWeightQSideQLongFctn(Int_t nxbins, Double_t maxXval, Double_t minXval,
                                                         Int_t nybins, Double_t maxYval, Double_t minYval):
 AliHBTTwoPairFctn2D(nxbins,maxXval,minXval,nybins,maxYval,minYval)
{
  //ctor
 fWriteNumAndDen = kTRUE;//change default behaviour
 Rename("wqsideqlongcf","Q_{side} Q_{long} Weight Correlation Function 2D");
}    
/*************************************************************************************/ 

void AliHBTWeightQSideQLongFctn::ProcessSameEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair)
{
  //process particles from same events (fills numerator)
  trackpair = CheckPair(trackpair);
//  partpair  = CheckPair(partpair);
  if ( trackpair && partpair)     
  {
//    Double_t weightPID=1.;
    Double_t weight = 1.0;
    if ( ( trackpair->Particle1()->GetPdgCode() == partpair->Particle1()->GetPdgCode()) &&
         ( trackpair->Particle2()->GetPdgCode() == partpair->Particle2()->GetPdgCode())    )
      {   
         weight=partpair->GetWeight();
      }   
//    Double_t weight=weightHBT*weightPID;
    fNumerator->Fill(trackpair->GetQSideLCMS(),trackpair->GetQLongLCMS(),weight);
  }
} 
/****************************************************************/

void AliHBTWeightQSideQLongFctn::ProcessDiffEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair)
{
  //process particles from diff events (fills denominator)
  trackpair = CheckPair(trackpair);
//  partpair  = CheckPair(partpair);
  if ( trackpair && partpair)  
  {
     fDenominator->Fill(trackpair->GetQSideLCMS(),trackpair->GetQLongLCMS());
  }
}
/**************************************************************/

TH1* AliHBTWeightQSideQLongFctn::GetResult()
{
 //returns the scaled ratio
 delete fRatio;
 fRatio = GetRatio(Scale());
 return fRatio;
}
/*************************************************************/
/*************************************************************/
/*************************************************************/ 

ClassImp(AliHBTWeightQOutSQideQLongFctn)

AliHBTWeightQOutSQideQLongFctn::AliHBTWeightQOutSQideQLongFctn(Int_t nXbins, Double_t maxXval, Double_t minXval,
                                                   Int_t nYbins, Double_t maxYval, Double_t minYval,
                                                   Int_t nZbins, Double_t maxZval, Double_t minZval):
 AliHBTTwoPairFctn3D(nXbins,maxXval,minXval,nYbins,maxYval,minYval,nZbins,maxZval,minZval)
{
  //ctor
  fWriteNumAndDen = kTRUE;//change default behaviour
  Rename("wqoslcf","Q_{out}-Q_{side}-Q_{long} Weight Correlation Fctn");
}
/*************************************************************/

void AliHBTWeightQOutSQideQLongFctn::ProcessSameEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair)
{
//Fills numerator
  //process particles from same events (fills numerator)
  trackpair = CheckPair(trackpair);
//  partpair  = CheckPair(partpair);
  if ( trackpair && partpair)     
  {
//    Double_t weightPID=1.;
    Double_t weight = 1.0;
    if ( ( trackpair->Particle1()->GetPdgCode() == partpair->Particle1()->GetPdgCode()) &&
         ( trackpair->Particle2()->GetPdgCode() == partpair->Particle2()->GetPdgCode())    )
      {   
         weight=partpair->GetWeight();//here we take weight from the particle pair
      }   
//    Double_t weight=weightHBT*weightPID;
    Double_t out = TMath::Abs(trackpair->GetQOutLCMS());
    Double_t side = TMath::Abs(trackpair->GetQSideLCMS());
    Double_t lon = TMath::Abs(trackpair->GetQLongLCMS());
    
    fNumerator->Fill(out,side,lon,weight);//here we fill in q's corresponding to track pair 
                                          //weight calculated for the simulated one
  }
}
/*************************************************************/

void AliHBTWeightQOutSQideQLongFctn::ProcessDiffEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair)
{
  //process particles from diff events (fills denominator)
  trackpair = CheckPair(trackpair);
//  partpair  = CheckPair(partpair);
  if ( trackpair && partpair)  
   {
     Double_t out = TMath::Abs(trackpair->GetQOutLCMS());
     Double_t side = TMath::Abs(trackpair->GetQSideLCMS());
     Double_t lon = TMath::Abs(trackpair->GetQLongLCMS());
   
     fDenominator->Fill(out,side,lon);
   }
}
/*************************************************************/

TH1* AliHBTWeightQOutSQideQLongFctn::GetResult()
{
 //returns the scaled ratio
 delete fRatio;
 fRatio = GetRatio(Scale());
 return fRatio;
}                    
