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

ClassImp(AliHBTLLWeightTheorQInvFctn)  
/*************************************************************/

AliHBTLLWeightTheorQInvFctn::
AliHBTLLWeightTheorQInvFctn(Int_t nbins, Double_t maxXval, Double_t minXval):
  AliHBTOnePairFctn1D(nbins,maxXval,minXval)
{
  //ctor
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
                                                              
