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
 //returns the scaled ratio
 delete fRatio;
 fRatio = GetRatio(Scale());
 return fRatio;
}
/*************************************************************************************/ 
/*************************************************************************************/ 
/*************************************************************************************/ 

ClassImp(AliHBTOutSideLongFctn)

AliHBTOutSideLongFctn::AliHBTOutSideLongFctn(Int_t nXbins, Double_t maxXval, Double_t minXval,
                                                   Int_t nYbins, Double_t maxYval, Double_t minYval,
                                                   Int_t nZbins, Double_t maxZval, Double_t minZval):
 AliHBTOnePairFctn3D(nXbins,maxXval,minXval,nYbins,maxYval,minYval,nZbins,maxZval,minZval),
 fAbs(kTRUE)
{
//ctor
  fWriteNumAndDen = kTRUE;//change default behaviour
  Rename("qoslcf","Q_{out}-Q_{side}-Q_{long} Correlation Fctn");
}
/*************************************************************************************/ 

TH1* AliHBTOutSideLongFctn::GetResult()
{
 //returns the scaled ratio
 delete fRatio;
 fRatio = GetRatio(Scale());
 return fRatio;
}

void AliHBTOutSideLongFctn::GetValues(AliHBTPair* pair, Double_t& x, Double_t& y, Double_t& z) const
{ 
  //calculates values of that function
  //qout qside and qlong
  
  x=pair->GetQOutCMSLC(); 
  y=pair->GetQSideCMSLC(); 
  z=pair->GetQLongCMSLC();
  if (fAbs)
   {
     x = TMath::Abs(x);
     y = TMath::Abs(y);
     z = TMath::Abs(z);
   }
} 

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
 //returns the scaled ratio
 delete fRatio;
 fRatio = GetRatio(Scale());
 return fRatio;
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
 //returns the scaled ratio
 delete fRatio;
 fRatio = GetRatio(Scale());
 return fRatio;
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
 //returns the scaled ratio
 delete fRatio;
 fRatio = GetRatio(Scale());
 return fRatio;
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
 return GetNumerator();
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
 //returns the scaled ratio
 delete fRatio;
 fRatio = GetRatio(Scale());
 return fRatio;
}

/*************************************************************************************/ 
/*************************************************************************************/ 
/*************************************************************************************/ 
ClassImp(AliHBTAvSeparCorrelFctn)

AliHBTAvSeparCorrelFctn::AliHBTAvSeparCorrelFctn(Int_t nbins, Double_t maxXval, Double_t minXval):
 AliHBTOnePairFctn1D(nbins,maxXval,minXval)
{
 //ctor 
 fWriteNumAndDen = kTRUE;//change default behaviour
 Rename("avsepcf","Avarage separation Correlation Function");
}

/*************************************************************************************/ 

TH1* AliHBTAvSeparCorrelFctn::GetResult()
{  
 //returns the scaled ratio
 delete fRatio;
 fRatio = GetRatio(Scale());
 return fRatio;
}

/*************************************************************************************/ 
