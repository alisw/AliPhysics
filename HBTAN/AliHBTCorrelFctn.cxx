#include "AliHBTCorrelFctn.h"



ClassImp(AliHBTQInvCorrelFctn)

//Corroleation function is created from dividing two histograms of QInvariant:
//  of particles from the same evnt
//by 
//  of particles from different events

TH1* AliHBTQInvCorrelFctn::GetResult()
{
 return GetRatio(GetDenominator()->GetMaximum()/GetNumerator()->GetMaximum());
}


ClassImp(AliHBTInvMassCorrelFctn)

AliHBTInvMassCorrelFctn::
AliHBTInvMassCorrelFctn(Int_t nbins, Double_t maxXval, Double_t minXval):
                        AliHBTTwoPartFctn1D(nbins,maxXval,minXval)
{
  Rename("InvMass CF","Invariant Mass Correlation Function");
}

TH1* AliHBTInvMassCorrelFctn::GetResult()
{
 return GetRatio(GetDenominator()->GetMaximum()/GetNumerator()->GetMaximum());
}
