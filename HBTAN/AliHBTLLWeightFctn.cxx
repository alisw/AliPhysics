#include "AliHBTLLWeightFctn.h"
/* $Id$ */
//_________________________________________________________________________
//
// class 
//
//This class allows to obtain Q_inv correlation function with weights
//calculated by Lednicky's alghorithm.
//Numerator is filled with weighted events. Weights are attributed to reconstructed tracks.
//Weights are calculated with corresponding simulated particles momenta.
//Denominator is filled with mixing unweighted reconstructed tracks.
//One needs both pairs 
//(simulated and recontructed), thus function is of class AliHBTTwoPairFctn1D.
//Author: Ludmila Malinina, JINR (malinina@sunhe.jinr.ru)
#include "AliHBTLLWeights.h"
#include "AliHBTLLWeightsPID.h"

//--for test--AliHBTLLWeightQInvFctn* yyy= new AliHBTLLWeightQInvFctn();

ClassImp( AliHBTLLWeightQInvFctn )

  
/****************************************************************/
AliHBTLLWeightQInvFctn::AliHBTLLWeightQInvFctn(Int_t nbins, Double_t maxXval, Double_t minXval):
 AliHBTTwoPairFctn1D(nbins,maxXval,minXval)
{
//ctor
 fWriteNumAndDen = kTRUE;//change default behaviour
 Rename("wqinvcf","Q_{inv} Weight Correlation Function");
}
/****************************************************************/
void  AliHBTLLWeightQInvFctn::ProcessSameEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair)
{
  //process particles from same events (fills numerator)
  trackpair = CheckPair(trackpair);
  partpair  = CheckPair(partpair);
  if ( trackpair && partpair)     
  {
    Double_t weightPID=1.;
    Double_t weightHBT=AliHBTLLWeights::Instance()->GetWeight(partpair);
    Double_t weight=weightHBT*weightPID;
    if(TMath::Abs(weight)<=10.)fNumerator->Fill(trackpair->GetQInv(),weight);
  }
} 
/****************************************************************/

void  AliHBTLLWeightQInvFctn::ProcessDiffEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair)
{
  // Fills the denominator using mixed pairs
  trackpair = CheckPair(trackpair);
  partpair  = CheckPair(partpair);
  if ( trackpair && partpair)  
  {
     fDenominator->Fill(trackpair->GetQInv());
  }
}
/**************************************************************/
TH1* AliHBTLLWeightQInvFctn::GetResult()
{ 
//returns ratio of numerator and denominator                                    
 return GetRatio(Scale());                                                  
}                    
                                                              
/**************************************************************************************/
/**************************************************************************************/
/**************************************************************************************/
/**************************************************************************************/

ClassImp(AliHBTLLWeightQOutFctn)
    
AliHBTLLWeightQOutFctn::AliHBTLLWeightQOutFctn(Int_t nbins, Double_t maxXval, Double_t minXval):
 AliHBTTwoPairFctn1D(nbins,maxXval,minXval)
{
//ctor
 fWriteNumAndDen = kTRUE;//change default behaviour
 Rename("wqoutcf","Q_{out} Weight Correlation Function");
}
/****************************************************************/
void AliHBTLLWeightQOutFctn::ProcessSameEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair)
{
  //process particles from same events (fills numerator)  
  trackpair = CheckPair(trackpair);
  partpair  = CheckPair(partpair);
  if ( trackpair && partpair)     
  {
    Double_t weightPID=1.;
    Double_t weightHBT=AliHBTLLWeights::Instance()->GetWeight(partpair);
    Double_t weight=weightHBT*weightPID;
    if(TMath::Abs(weight)<=10.) fNumerator->Fill(trackpair->GetQOutCMSLC(),weight);
  }
} 
/****************************************************************/

void AliHBTLLWeightQOutFctn::ProcessDiffEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair)
{
  //process particles from diff events (fills denominator)
  trackpair = CheckPair(trackpair);
  partpair  = CheckPair(partpair);
  if ( trackpair && partpair)  
  {
     fDenominator->Fill(trackpair->GetQOutCMSLC());
  }
}
/**************************************************************/
TH1* AliHBTLLWeightQOutFctn::GetResult() 
                                                                               
{ 
//returns ratio of numerator and denominator                                    
 return GetRatio(Scale());                                                  
}                    
                                                              
/*************************************************************************************/ 
/*************************************************************************************/ 
/*************************************************************************************/ 

ClassImp(AliHBTLLWeightQLongFctn)
AliHBTLLWeightQLongFctn::AliHBTLLWeightQLongFctn(Int_t nbins, Double_t maxXval, Double_t minXval):
 AliHBTTwoPairFctn1D(nbins,maxXval,minXval)
{
//ctor
 fWriteNumAndDen = kTRUE;//change default behaviour
 Rename("wqlongcf","Q_{long} Weight Correlation Function");
}
/****************************************************************/
void AliHBTLLWeightQLongFctn::ProcessSameEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair)
{
  //process particles from same events (fills numerator)
  trackpair = CheckPair(trackpair);
  partpair  = CheckPair(partpair);
  if ( trackpair && partpair)     
  {
    Double_t weightPID=1.;
    Double_t weightHBT=AliHBTLLWeights::Instance()->GetWeight(partpair);
    Double_t weight=weightHBT*weightPID;
    if(TMath::Abs(weight)<=10.) fNumerator->Fill(trackpair->GetQLongCMSLC(),weight);
  }
} 
/****************************************************************/

void AliHBTLLWeightQLongFctn::ProcessDiffEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair)
{
  //process particles from diff events (fills denominator)
  trackpair = CheckPair(trackpair);
  partpair  = CheckPair(partpair);
  if ( trackpair && partpair)  
  {
     fDenominator->Fill(trackpair->GetQLongCMSLC());
  }
}
/**************************************************************/
TH1* AliHBTLLWeightQLongFctn::GetResult() 
                                                                               
{ 
//returns ratio of numerator and denominator                                    
 return GetRatio(Scale());                                                  
}                    
                                                              
/*************************************************************************************/ 
/*************************************************************************************/ 
/*************************************************************************************/ 

ClassImp(AliHBTLLWeightQSideFctn)
/*************************************************************************************/ 

AliHBTLLWeightQSideFctn::AliHBTLLWeightQSideFctn(Int_t nbins, Double_t maxXval, Double_t minXval):
 AliHBTTwoPairFctn1D(nbins,maxXval,minXval)
{
//ctor
 fWriteNumAndDen = kTRUE;//change default behaviour
 Rename("wqsidecf","Q_{side} Weight Correlation Function");
}
/****************************************************************/
void AliHBTLLWeightQSideFctn::ProcessSameEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair)
{
  //process particles from same events (fills numerator)
  trackpair = CheckPair(trackpair);
  partpair  = CheckPair(partpair);
  if ( trackpair && partpair)     
  {
    Double_t weightPID=1.;
    Double_t weightHBT=AliHBTLLWeights::Instance()->GetWeight(partpair);
    Double_t weight=weightHBT*weightPID;
    if(TMath::Abs(weight)<=10.) fNumerator->Fill(trackpair->GetQSideCMSLC(),weight);
  }
} 
/****************************************************************/

void  AliHBTLLWeightQSideFctn::ProcessDiffEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair)
{
  //process particles from diff events (fills denominator)
  trackpair = CheckPair(trackpair);
  partpair  = CheckPair(partpair);
  if ( trackpair && partpair)  
  {
     fDenominator->Fill(trackpair->GetQSideCMSLC());
  }
}
/**************************************************************/
TH1* AliHBTLLWeightQSideFctn::GetResult() 
                                                                               
{ 
//returns ratio of numerator and denominator                                    
 return GetRatio(Scale());                                                  
}                    
                                                              
/*************************************************************************************/ 
/*************************************************************************************/ 
/*************************************************************************************/ 

ClassImp(AliHBTLLWeightTwoKStarFctn)
/*************************************************************************************/ 
AliHBTLLWeightTwoKStarFctn::AliHBTLLWeightTwoKStarFctn(Int_t nbins, Double_t maxXval, Double_t minXval):
 AliHBTTwoPairFctn1D(nbins,maxXval,minXval)
{
//ctor
 fWriteNumAndDen = kTRUE;//change default behaviour
 Rename("wtwokstarcf","2*K^{*} Weight Correlation Function");
}
/****************************************************************/
void AliHBTLLWeightTwoKStarFctn::ProcessSameEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair)
{
  //process particles from same events (fills numerator)
  trackpair = CheckPair(trackpair);
  partpair  = CheckPair(partpair);
  if ( trackpair && partpair)     
  {
    Double_t weightPID=1.;
    Double_t weightHBT=AliHBTLLWeights::Instance()->GetWeight(partpair);
    Double_t weight=weightHBT*weightPID;
    if(TMath::Abs(weight)<=10.) fNumerator->Fill(2.0*(trackpair->GetKStar()),weight);
  }
} 
/****************************************************************/

void  AliHBTLLWeightTwoKStarFctn::ProcessDiffEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair)
{
  //process particles from diff events (fills denominator)
  trackpair = CheckPair(trackpair);
  partpair  = CheckPair(partpair);
  if ( trackpair && partpair)  
  {
     fDenominator->Fill(2.0*(trackpair->GetKStar()));
  }
}
/**************************************************************/
TH1* AliHBTLLWeightTwoKStarFctn::GetResult() 
                                                                               
{ 
//returns ratio of numerator and denominator                                    
 return GetRatio(Scale());                                                  
}                    
                                                              
/*************************************************************************************/ 
/*************************************************************************************/ 
/*************************************************************************************/ 

ClassImp(AliHBTLLWeightQOutQSideFctn)
/*************************************************************************************/ 
    
AliHBTLLWeightQOutQSideFctn::AliHBTLLWeightQOutQSideFctn(Int_t nxbins, Double_t maxXval, Double_t minXval,
                                                         Int_t nybins, Double_t maxYval, Double_t minYval):
 AliHBTTwoPairFctn2D(nxbins,maxXval,minXval,nybins,maxYval,minYval)
{
  //ctor
 fWriteNumAndDen = kTRUE;//change default behaviour
 Rename("wqoutqsidecf","Q_{out} Q_{side} Weight Correlation Function 2D");
}    
/*************************************************************************************/ 
void AliHBTLLWeightQOutQSideFctn::ProcessSameEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair)
{
  //process particles from same events (fills numerator)
  trackpair = CheckPair(trackpair);
  partpair  = CheckPair(partpair);
  if ( trackpair && partpair)     
  {
    Double_t weightPID=1.;
    Double_t weightHBT=AliHBTLLWeights::Instance()->GetWeight(partpair);
    Double_t weight=weightHBT*weightPID;
    if(TMath::Abs(weight)<=10.) 
      fNumerator->Fill(trackpair->GetQOutCMSLC(),trackpair->GetQSideCMSLC(),weight);
  }
} 
/****************************************************************/

void AliHBTLLWeightQOutQSideFctn::ProcessDiffEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair)
{
  //process particles from diff events (fills denominator)
  trackpair = CheckPair(trackpair);
  partpair  = CheckPair(partpair);
  if ( trackpair && partpair)  
  {
     fDenominator->Fill(trackpair->GetQOutCMSLC(),trackpair->GetQSideCMSLC());
  }
}
/**************************************************************/
TH1* AliHBTLLWeightQOutQSideFctn::GetResult()
{
  //returns result
  return GetRatio(Scale());
}

/*************************************************************************************/ 
/*************************************************************************************/ 
/*************************************************************************************/ 

ClassImp(AliHBTLLWeightQOutQLongFctn)
/*************************************************************************************/ 
    
AliHBTLLWeightQOutQLongFctn::AliHBTLLWeightQOutQLongFctn(Int_t nxbins, Double_t maxXval, Double_t minXval,
                                                         Int_t nybins, Double_t maxYval, Double_t minYval):
 AliHBTTwoPairFctn2D(nxbins,maxXval,minXval,nybins,maxYval,minYval)
{
  //ctor
 fWriteNumAndDen = kTRUE;//change default behaviour
 Rename("wqoutqlongcf","Q_{out} Q_{long} Weight Correlation Function 2D");
}    
/*************************************************************************************/ 
void AliHBTLLWeightQOutQLongFctn::ProcessSameEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair)
{
  //process particles from same events (fills numerator)
  trackpair = CheckPair(trackpair);
  partpair  = CheckPair(partpair);
  if ( trackpair && partpair)     
  {
    Double_t weightPID=1.;
    Double_t weightHBT=AliHBTLLWeights::Instance()->GetWeight(partpair);
    Double_t weight=weightHBT*weightPID;
    if(TMath::Abs(weight)<=10.) 
      fNumerator->Fill(trackpair->GetQOutCMSLC(),trackpair->GetQLongCMSLC(),weight);
  }
} 
/****************************************************************/

void AliHBTLLWeightQOutQLongFctn::ProcessDiffEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair)
{
  //process particles from diff events (fills denominator)
  trackpair = CheckPair(trackpair);
  partpair  = CheckPair(partpair);
  if ( trackpair && partpair)  
  {
     fDenominator->Fill(trackpair->GetQOutCMSLC(),trackpair->GetQLongCMSLC());
  }
}
/**************************************************************/

TH1* AliHBTLLWeightQOutQLongFctn::GetResult()
{
  //returns result
  return GetRatio(Scale());
}

/*************************************************************************************/ 
/*************************************************************************************/ 
/*************************************************************************************/ 

ClassImp(AliHBTLLWeightQSideQLongFctn)
/*************************************************************************************/ 
    
AliHBTLLWeightQSideQLongFctn::AliHBTLLWeightQSideQLongFctn(Int_t nxbins, Double_t maxXval, Double_t minXval,
                                                         Int_t nybins, Double_t maxYval, Double_t minYval):
 AliHBTTwoPairFctn2D(nxbins,maxXval,minXval,nybins,maxYval,minYval)
{
  //ctor
 fWriteNumAndDen = kTRUE;//change default behaviour
 Rename("wqsideqlongcf","Q_{side} Q_{long} Weight Correlation Function 2D");
}    
/*************************************************************************************/ 
void AliHBTLLWeightQSideQLongFctn::ProcessSameEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair)
{
  //process particles from same events (fills numerator)
  trackpair = CheckPair(trackpair);
  partpair  = CheckPair(partpair);
  if ( trackpair && partpair)     
  {
    Double_t weightPID=1.;
    Double_t weightHBT=AliHBTLLWeights::Instance()->GetWeight(partpair);
    Double_t weight=weightHBT*weightPID;
    if(TMath::Abs(weight)<=10.) 
      fNumerator->Fill(trackpair->GetQSideCMSLC(),trackpair->GetQLongCMSLC(),weight);
  }
} 
/****************************************************************/

void AliHBTLLWeightQSideQLongFctn::ProcessDiffEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair)
{
  //process particles from diff events (fills denominator)
  trackpair = CheckPair(trackpair);
  partpair  = CheckPair(partpair);
  if ( trackpair && partpair)  
  {
     fDenominator->Fill(trackpair->GetQSideCMSLC(),trackpair->GetQLongCMSLC());
  }
}
/**************************************************************/
TH1* AliHBTLLWeightQSideQLongFctn::GetResult()
{
  //returns result
  return GetRatio(Scale());
}
