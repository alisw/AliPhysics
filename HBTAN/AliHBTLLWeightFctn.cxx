/* $Id$ */

//This class allows to obtain Q_inv correlation function with weights
//calculated by Lednicky's alghorithm.
//Numerator is filled with weighted events. Weights are attributed to reconstructed tracks.
//Weights are calculated with corresponding simulated particles momenta.
//Denominator is filled with mixing unweighted reconstructed tracks.
//One needs both pairs 
//(simulated and recontructed), thus function is of class AliHBTTwoPairFctn1D.
//Author: Ludmila Malinina, JINR (malinina@sunhe.jinr.ru)
#include "AliHBTLLWeightFctn.h"
#include "AliHBTLLWeights.h"

//--for test--AliHBTLLWeightQInvFctn* yyy= new AliHBTLLWeightQInvFctn();

ClassImp( AliHBTLLWeightQInvFctn )  
/****************************************************************/
AliHBTLLWeightQInvFctn::AliHBTLLWeightQInvFctn(Int_t nbins, Double_t maxXval, Double_t minXval):
           AliHBTTwoPairFctn1D(nbins,maxXval,minXval)
{
//ctor
}
/****************************************************************/
void  AliHBTLLWeightQInvFctn::ProcessSameEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair)
{
//Processes Particles and tracks Same different even
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
/****************************************************************/
void AliHBTLLWeightQOutFctn::ProcessSameEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair)
{
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
/****************************************************************/
void  AliHBTLLWeightQLongFctn::ProcessSameEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair)
{
//Processes Particles and tracks Same different even
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

void  AliHBTLLWeightQLongFctn::ProcessDiffEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair)
{
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
/****************************************************************/
void  AliHBTLLWeightQSideFctn::ProcessSameEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair)
{
//Processes Particles and tracks Same different even
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

ClassImp(AliHBTLLWeightTwoKStarFctn)
/****************************************************************/
void  AliHBTLLWeightTwoKStarFctn::ProcessSameEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair)
{
//Processes Particles and tracks Same different even
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

