#include "AliHBTCorrelFctn.h"
//____________________________________________________________________________
//////////////////////////////////////////////////////////////////////////////
//
// class AliHBTQInvCorrelFctn
// class AliHBTQOutLCMSCorrelFctn
// class AliHBTQLongLCMSCorrelFctn
// class AliHBTQSideLCMSCorrelFctn
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
 AliHBTOnePairFctn3D(nXbins,maxXval,minXval,nYbins,maxYval,minYval,nZbins,maxZval,minZval)
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
  
  x=pair->GetQOutLCMS(); 
  y=pair->GetQSideLCMS(); 
  z=pair->GetQLongLCMS();
  if (fAbs)
   {
     x = TMath::Abs(x);
     y = TMath::Abs(y);
     z = TMath::Abs(z);
   }
} 

/*************************************************************************************/ 

ClassImp(AliHBTQOutLCMSCorrelFctn)
    
AliHBTQOutLCMSCorrelFctn::AliHBTQOutLCMSCorrelFctn(Int_t nbins, Double_t maxXval, Double_t minXval):
 AliHBTOnePairFctn1D(nbins,maxXval,minXval)
{
  //ctor
 fWriteNumAndDen = kTRUE;//change default behaviour
 Rename("qoutcf","Q_{out} Correlation Function");
}
/*************************************************************************************/ 
    
TH1* AliHBTQOutLCMSCorrelFctn::GetResult()
{
 //returns the scaled ratio
 delete fRatio;
 fRatio = GetRatio(Scale());
 return fRatio;
}
/*************************************************************************************/ 
/*************************************************************************************/ 
/*************************************************************************************/ 

ClassImp(AliHBTQLongLCMSCorrelFctn)
    
AliHBTQLongLCMSCorrelFctn::AliHBTQLongLCMSCorrelFctn(Int_t nbins, Double_t maxXval, Double_t minXval):
 AliHBTOnePairFctn1D(nbins,maxXval,minXval)
{
  //ctor
 fWriteNumAndDen = kTRUE;//change default behaviour
 Rename("qlongcf","Q_{long} Correlation Function");
}
/*************************************************************************************/ 
    
TH1* AliHBTQLongLCMSCorrelFctn::GetResult()
{
 //returns the scaled ratio
 delete fRatio;
 fRatio = GetRatio(Scale());
 return fRatio;
}
/*************************************************************************************/ 
/*************************************************************************************/ 
/*************************************************************************************/ 

ClassImp(AliHBTQSideLCMSCorrelFctn)
    
AliHBTQSideLCMSCorrelFctn::AliHBTQSideLCMSCorrelFctn(Int_t nbins, Double_t maxXval, Double_t minXval):
 AliHBTOnePairFctn1D(nbins,maxXval,minXval)
{
 //ctor
 fWriteNumAndDen = kTRUE;//change default behaviour
 Rename("qsidecf","Q_{side} Correlation Function");
}
/*************************************************************************************/ 
    
TH1* AliHBTQSideLCMSCorrelFctn::GetResult()
{
 //returns the scaled ratio
 delete fRatio;
 fRatio = GetRatio(Scale());
 return fRatio;
}


/*************************************************************************************/ 
/*************************************************************************************/ 
/*************************************************************************************/ 
ClassImp(AliHBTQtLCMSCorrelFctn)
    
AliHBTQtLCMSCorrelFctn::AliHBTQtLCMSCorrelFctn(Int_t nbins, Double_t maxXval, Double_t minXval):
 AliHBTOnePairFctn1D(nbins,maxXval,minXval)
{
  //ctor
 fWriteNumAndDen = kTRUE;//change default behaviour
 Rename("Qtcf","Q_{t}(LCMS) Correlation Function");
}
/*************************************************************************************/ 
    
TH1* AliHBTQtLCMSCorrelFctn::GetResult()
{
 //returns the scaled ratio
 delete fRatio;
 fRatio = GetRatio(Scale());
 return fRatio;
}
/*************************************************************************************/ 
/*************************************************************************************/ 
/*************************************************************************************/ 
ClassImp(AliHBTQtCorrelFctn)
    
AliHBTQtCorrelFctn::AliHBTQtCorrelFctn(Int_t nbins, Double_t maxXval, Double_t minXval):
 AliHBTOnePairFctn1D(nbins,maxXval,minXval)
{
  //ctor
 fWriteNumAndDen = kTRUE;//change default behaviour
 Rename("qtcf","Q_{t} Correlation Function");
}
/*************************************************************************************/ 
    
TH1* AliHBTQtCorrelFctn::GetResult()
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

ClassImp(AliHBTAvSeparVsQInvCorrelFctn)

AliHBTAvSeparVsQInvCorrelFctn::AliHBTAvSeparVsQInvCorrelFctn(Int_t nXbins, Double_t maxXval, Double_t minXval,
                                                             Int_t nYbins, Double_t maxYval, Double_t minYval):
 AliHBTOnePairFctn2D(nXbins,maxXval,minXval,nYbins,maxYval,minYval)
{
 //ctor 
 fWriteNumAndDen = kTRUE;//change default behaviour
 Rename("avsepvsqinv","Avarage Separation VS Q_{inv} Correlation Function");
}


TH1* AliHBTAvSeparVsQInvCorrelFctn::GetResult()
{  
 //returns the scaled ratio
 delete fRatio;
 fRatio = GetRatio(Scale());
 return fRatio;
}

/**************************************************************/
/**************************************************************/
/**************************************************************/


ClassImp(AliHBTQOutQSideFctn)


AliHBTQOutQSideFctn::AliHBTQOutQSideFctn(Int_t nxbins, Double_t maxXval, Double_t minXval,
                                         Int_t nybins, Double_t maxYval, Double_t minYval):
 AliHBTOnePairFctn2D(nxbins,maxXval,minXval,nybins,maxYval,minYval)
{
  //ctor
 fWriteNumAndDen = kTRUE;//change default behaviour
 Rename("qoutqsidecf","Q_{out} Q_{side} Correlation Function 2D");
}    
/**************************************************************/

TH1* AliHBTQOutQSideFctn::GetResult()
{
 //returns the scaled ratio
 delete fRatio;
 fRatio = GetRatio(Scale());
 return fRatio;
}
/**************************************************************/
/**************************************************************/

ClassImp(AliHBTQOutQLongFctn)

AliHBTQOutQLongFctn::AliHBTQOutQLongFctn(Int_t nxbins, Double_t maxXval, Double_t minXval,
                                                         Int_t nybins, Double_t maxYval, Double_t minYval):
 AliHBTOnePairFctn2D(nxbins,maxXval,minXval,nybins,maxYval,minYval)
{
  //ctor
 fWriteNumAndDen = kTRUE;//change default behaviour
 Rename("qoutqlongcf","Q_{out} Q_{long} Correlation Function 2D");
}    


/**************************************************************/

TH1* AliHBTQOutQLongFctn::GetResult()
{
 //returns the scaled ratio
 delete fRatio;
 fRatio = GetRatio(Scale());
 return fRatio;
}
/**************************************************************/
/**************************************************************/
/**************************************************************/

ClassImp(AliHBTQSideQLongFctn)

/**************************************************************/
AliHBTQSideQLongFctn::AliHBTQSideQLongFctn(Int_t nxbins, Double_t maxXval, Double_t minXval,
                                                         Int_t nybins, Double_t maxYval, Double_t minYval):
 AliHBTOnePairFctn2D(nxbins,maxXval,minXval,nybins,maxYval,minYval)
{
  //ctor
 fWriteNumAndDen = kTRUE;//change default behaviour
 Rename("qsideqlongcf","Q_{side} Q_{long} Correlation Function 2D");
}    

TH1* AliHBTQSideQLongFctn::GetResult()
{
 //returns the scaled ratio
 delete fRatio;
 fRatio = GetRatio(Scale());
 return fRatio;
}
