/* $Id$ */

//-----------------------------------------------------------
// This class introduces the weights calculated according 
// with functions of efficiency of identification (TPC+TOF) 
// (calculated by B.V. Batyunia).
// Author: Ludmila Malinina, JINR (malinina@sunhe.jinr.ru)
//-----------------------------------------------------------

#include "AliHBTLLWeightTheorFctn.h"
#include "AliHBTLLWeights.h"

//--for test--AliHBTLLWeightQInvFctn* yyy= new AliHBTLLWeightQInvFctn();
/*************************************************************/
/*************************************************************/
/*************************************************************/

ClassImp(AliHBTLLWeightTheorQInvFctn)  
/*************************************************************/

AliHBTLLWeightTheorQInvFctn::
AliHBTLLWeightTheorQInvFctn(Int_t nbins, Double_t maxXval, Double_t minXval):
  AliHBTOnePairFctn1D(nbins,maxXval,minXval)
{
  //ctor
 fWriteNumAndDen = kTRUE;//change default behaviour
 Rename("wqinvtheorcf","Q_{inv} Weight Theoretical Correlation Function");
}
/****************************************************************/
void  AliHBTLLWeightTheorQInvFctn::ProcessSameEventParticles(AliHBTPair* partpair)
{
  //Processes Particles and tracks Same different even
  partpair  = CheckPair(partpair);
  Double_t weight = AliHBTLLWeights::Instance()->GetWeight(partpair);
  if(TMath::Abs(weight)<=10.) fNumerator->Fill(partpair->GetQInv(),weight);
} 

/**************************************************************/
TH1* AliHBTLLWeightTheorQInvFctn::GetResult() 
{
  //returns ratio of numerator and denominator
  return GetRatio(Scale());
}                    
                                                              
/*************************************************************/
/*************************************************************/
/*************************************************************/

ClassImp(AliHBTLLWeightTheorQOutFctn)  
/*************************************************************/

AliHBTLLWeightTheorQOutFctn::
AliHBTLLWeightTheorQOutFctn(Int_t nbins, Double_t maxXval, Double_t minXval):
  AliHBTOnePairFctn1D(nbins,maxXval,minXval)
{
  //ctor
 fWriteNumAndDen = kTRUE;//change default behaviour
 Rename("wqouttheorcf","Q_{out} Weight Theoretical Correlation Function");
}
/****************************************************************/
void  AliHBTLLWeightTheorQOutFctn::ProcessSameEventParticles(AliHBTPair* partpair)
{
  //Processes Particles and tracks Same different even
  partpair  = CheckPair(partpair);
  Double_t weight = AliHBTLLWeights::Instance()->GetWeight(partpair);
  if(TMath::Abs(weight)<=10.) fNumerator->Fill(partpair->GetQOutCMSLC(),weight);
} 

/**************************************************************/
TH1* AliHBTLLWeightTheorQOutFctn::GetResult() 
{
  //returns ratio of numerator and denominator
  return GetRatio(Scale());
}                    

/*************************************************************/
/*************************************************************/
/*************************************************************/

ClassImp(AliHBTLLWeightTheorQSideFctn)  
/*************************************************************/

AliHBTLLWeightTheorQSideFctn::
AliHBTLLWeightTheorQSideFctn(Int_t nbins, Double_t maxXval, Double_t minXval):
  AliHBTOnePairFctn1D(nbins,maxXval,minXval)
{
  //ctor
 fWriteNumAndDen = kTRUE;//change default behaviour
 Rename("wqsidetheorcf","Q_{side} Weight Theoretical Correlation Function");
}
/****************************************************************/
void  AliHBTLLWeightTheorQSideFctn::ProcessSameEventParticles(AliHBTPair* partpair)
{
  //Processes Particles and tracks Same different even
  partpair  = CheckPair(partpair);
  Double_t weight = AliHBTLLWeights::Instance()->GetWeight(partpair);
  if(TMath::Abs(weight)<=10.) fNumerator->Fill(partpair->GetQSideCMSLC(),weight);
} 

/**************************************************************/
TH1* AliHBTLLWeightTheorQSideFctn::GetResult() 
{
  //returns ratio of numerator and denominator
  return GetRatio(Scale());
}                    

/*************************************************************/
/*************************************************************/
/*************************************************************/

ClassImp(AliHBTLLWeightTheorQLongFctn)  
/*************************************************************/

AliHBTLLWeightTheorQLongFctn::
AliHBTLLWeightTheorQLongFctn(Int_t nbins, Double_t maxXval, Double_t minXval):
  AliHBTOnePairFctn1D(nbins,maxXval,minXval)
{
  //ctor
 fWriteNumAndDen = kTRUE;//change default behaviour
 Rename("wqlongtheorcf","Q_{long} Weight Theoretical Correlation Function");
}
/****************************************************************/
void  AliHBTLLWeightTheorQLongFctn::ProcessSameEventParticles(AliHBTPair* partpair)
{
  //Processes Particles and tracks Same different even
  partpair  = CheckPair(partpair);
  Double_t weight = AliHBTLLWeights::Instance()->GetWeight(partpair);
  if(TMath::Abs(weight)<=10.) fNumerator->Fill(partpair->GetQLongCMSLC(),weight);
} 

/**************************************************************/
TH1* AliHBTLLWeightTheorQLongFctn::GetResult() 
{
  //returns ratio of numerator and denominator
  return GetRatio(Scale());
}                    
