#include <Riostream.h>
#include "AliHBTLLWeightFctn.h"
#include "AliHBTLLWeights.h"

//--for test--AliHBTLLWeightQInvFctn* yyy= new AliHBTLLWeightQInvFctn();

ClassImp( AliHBTLLWeightQInvFctn )  
/****************************************************************/
AliHBTLLWeightQInvFctn::AliHBTLLWeightQInvFctn(Int_t nbins, Double_t maxXval, Double_t minXval):
           AliHBTTwoPairFctn1D(nbins,maxXval,minXval)
{
  Rename("Correlation function, method of weights(Lednicky's algorithm)");
}
/****************************************************************/
void  AliHBTLLWeightQInvFctn::ProcessSameEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair)
{

  trackpair = CheckPair(trackpair);
  partpair  = CheckPair(partpair);
  if ( trackpair && partpair)     
  {
     fNumerator->Fill(trackpair->GetQInv(),
          AliHBTLLWeights::Instance()->GetWeight(partpair));
  }
  
} 
                

/****************************************************************/
 void  AliHBTLLWeightQInvFctn::ProcessDiffEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair)
{
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
 TH1* res = GetRatio(Scale());                                                  
  
 if(res)
   {                                                                             
    res->GetXaxis()->SetTitle("Qinv [GeV/c]");                                       
    res->GetYaxis()->SetTitle("C(Qinv)");                                          
    res->SetTitle(GetTitle());
   }                                                                             
  return res;                                                                     
}                    
                                                              
