#include "AliHBTCorrelFctn.h"
//____________________________________________________________________________
//////////////////////////////////////////////////////////////////////////////
//
// class AliHBTQInvCorrelFctn
// class AliHBTQOutCMSLCCorrelFctn
// class AliHBTQLongCMSLCCorrelFctn
// class AliHBTQSideCMSLCCorrelFctn
// class AliHBTInvMassCorrelFctn
// class AliHBTTwoKStarCorrelFctn
//
// Set of functions:
//   Q Invaraint Correlation Function
//   Invariant Mass Function
//
// more info: http://aliweb.cern.ch/people/skowron/analyzer/index.html
// Piotr.Skowronski@cern.ch
//
//////////////////////////////////////////////////////////////////////////////

ClassImp(AliHBTQInvCorrelFctn)

//Corroleation function is created from dividing two histograms of QInvariant:
//  of particles from the same evnt
//by 
//  of particles from different events

AliHBTQInvCorrelFctn::AliHBTQInvCorrelFctn(Int_t nbins, Double_t maxXval, Double_t minXval):
 AliHBTOnePairFctn1D(nbins,maxXval,minXval)
{
 fWriteNumAndDen = kTRUE;//change default behaviour
 Rename("qinvcf","Q_{inv} Correlation Function");
}
/*************************************************************************************/ 

TH1* AliHBTQInvCorrelFctn::GetResult()
{  
 return GetRatio(Scale());
}
/*************************************************************************************/ 
/*************************************************************************************/ 
/*************************************************************************************/ 

ClassImp(AliHBTQOutCMSLCCorrelFctn)
    
AliHBTQOutCMSLCCorrelFctn::AliHBTQOutCMSLCCorrelFctn(Int_t nbins, Double_t maxXval, Double_t minXval):
 AliHBTOnePairFctn1D(nbins,maxXval,minXval)
{
  //ctor
 fWriteNumAndDen = kTRUE;//change default behaviour
 Rename("qoutcf","Q_{out} Correlation Function");
}
/*************************************************************************************/ 
    
TH1* AliHBTQOutCMSLCCorrelFctn::GetResult()
{
 //returns result of the function
 return GetRatio(Scale());
}
/*************************************************************************************/ 
/*************************************************************************************/ 
/*************************************************************************************/ 

ClassImp(AliHBTQLongCMSLCCorrelFctn)
    
AliHBTQLongCMSLCCorrelFctn::AliHBTQLongCMSLCCorrelFctn(Int_t nbins, Double_t maxXval, Double_t minXval):
 AliHBTOnePairFctn1D(nbins,maxXval,minXval)
{
  //ctor
 fWriteNumAndDen = kTRUE;//change default behaviour
 Rename("qlongcf","Q_{long} Correlation Function");
}
/*************************************************************************************/ 
    
TH1* AliHBTQLongCMSLCCorrelFctn::GetResult()
{
 //returns result of the function
 return GetRatio(Scale());
}
/*************************************************************************************/ 
/*************************************************************************************/ 
/*************************************************************************************/ 

ClassImp(AliHBTQSideCMSLCCorrelFctn)
    
AliHBTQSideCMSLCCorrelFctn::AliHBTQSideCMSLCCorrelFctn(Int_t nbins, Double_t maxXval, Double_t minXval):
 AliHBTOnePairFctn1D(nbins,maxXval,minXval)
{
 //ctor
 fWriteNumAndDen = kTRUE;//change default behaviour
 Rename("qsidecf","Q_{side} Correlation Function");
}
/*************************************************************************************/ 
    
TH1* AliHBTQSideCMSLCCorrelFctn::GetResult()
{
 //returns result
 return GetRatio(Scale());
}


/*************************************************************************************/ 
/*************************************************************************************/ 
/*************************************************************************************/ 

ClassImp(AliHBTInvMassCorrelFctn)

AliHBTInvMassCorrelFctn::AliHBTInvMassCorrelFctn(Int_t nbins, Double_t maxXval, Double_t minXval):
 AliHBTOnePairFctn1D(nbins,maxXval,minXval)
{
 //ctor 
 fWriteNumAndDen = kTRUE;//change default behaviour
 Rename("InvMass CF","Invariant Mass Correlation Function");
}

TH1* AliHBTInvMassCorrelFctn::GetResult()
{
 //returns result
 TString name = fName + " Result";
 return (TH1*)GetNumerator()->Clone(name.Data());
}
/*************************************************************************************/ 
/*************************************************************************************/ 
/*************************************************************************************/ 

ClassImp(AliHBTTwoKStarCorrelFctn)

AliHBTTwoKStarCorrelFctn::AliHBTTwoKStarCorrelFctn(Int_t nbins, Double_t maxXval, Double_t minXval):
 AliHBTOnePairFctn1D(nbins,maxXval,minXval)
{
 //ctor 
 fWriteNumAndDen = kTRUE;//change default behaviour
 Rename("twokstarcf","2K^{*} Correlation Function");
}

/*************************************************************************************/ 

TH1* AliHBTTwoKStarCorrelFctn::GetResult()
{  
 //returns result
 return GetRatio(Scale());
}

/*************************************************************************************/ 
