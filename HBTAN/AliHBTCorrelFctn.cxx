#include "AliHBTCorrelFctn.h"
//_____________________________________________________________
///////////////////////////////////////////////////////////////
//
//  Set of Correlation fuctions
//  AliHBTQInvCorrelFctn - Q Invariant correlatyion function
//
//  Corroleation function is created from dividing two histograms of QInvariant:
//    of particles from the same evnt
//  by 
//    of particles from different events
//
///////////////////////////////////////////////////////////////

ClassImp(AliHBTQInvCorrelFctn)


AliHBTQInvCorrelFctn::
AliHBTQInvCorrelFctn(Int_t nbins, Double_t maxXval, Double_t minXval):
                     AliHBTOnePairFctn1D(nbins,maxXval,minXval)
{
 Rename("qinvcf","Q_{inv} Correlation Function");
}


TH1* AliHBTQInvCorrelFctn::GetResult()
{  
 return GetRatio(Scale());
}
/*************************************************************************************/ 
/*************************************************************************************/ 
/*************************************************************************************/ 

ClassImp(AliHBTQOutCMSLCCorrelFctn)
TH1* AliHBTQOutCMSLCCorrelFctn::GetResult()
{
 return GetRatio(Scale());
}
/*************************************************************************************/ 
/*************************************************************************************/ 
/*************************************************************************************/ 

ClassImp(AliHBTQLongCMSLCCorrelFctn)
TH1* AliHBTQLongCMSLCCorrelFctn::GetResult()
{
 return GetRatio(Scale());
}
/*************************************************************************************/ 
/*************************************************************************************/ 
/*************************************************************************************/ 

ClassImp(AliHBTQSideCMSLCCorrelFctn)
TH1* AliHBTQSideCMSLCCorrelFctn::GetResult()
{
 return GetRatio(Scale());
}


/*************************************************************************************/ 
/*************************************************************************************/ 
/*************************************************************************************/ 

ClassImp(AliHBTInvMassCorrelFctn)

AliHBTInvMassCorrelFctn::
AliHBTInvMassCorrelFctn(Int_t nbins, Double_t maxXval, Double_t minXval):
                        AliHBTOnePairFctn1D(nbins,maxXval,minXval)
{
  Rename("InvMass CF","Invariant Mass Correlation Function");
}

TH1* AliHBTInvMassCorrelFctn::GetResult()
{
 TString name = fName + " Result";
 return (TH1*)GetNumerator()->Clone(name.Data());
}
/*************************************************************************************/ 
/*************************************************************************************/ 
/*************************************************************************************/ 

ClassImp(AliHBTTwoKStarCorrelFctn)

AliHBTTwoKStarCorrelFctn::
AliHBTTwoKStarCorrelFctn(Int_t nbins, Double_t maxXval, Double_t minXval):
                     AliHBTOnePairFctn1D(nbins,maxXval,minXval)
{
 Rename("twokstarcf","2K^{*} Correlation Function");
}

/*************************************************************************************/ 

TH1* AliHBTTwoKStarCorrelFctn::GetResult()
{  
 return GetRatio(Scale());
}

/*************************************************************************************/ 
