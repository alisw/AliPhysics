#include "AliHBTTwoTrackEffFctn.h"
//____________________________________________________________________
//////////////////////////////////////////////////////////////////////
//                                                                  //
//  class AliHBTTwoTrackEffFctn                                     //
//                                                                  //
//  classes for calculating two track efficiency of the tracking    //
//  binning is done using value of simulated pair montum difference // 
//  pair must be recontructed, that is why we need both pairs       //
//  (simulated and recontructed), thus functions are "two pair"     //
//  Piotr.Skowronski@cern.ch                                        //
//                                                                  //
//////////////////////////////////////////////////////////////////////


ClassImp(AliHBTTwoTrackEffFctn)
/******************************************************************/

AliHBTTwoTrackEffFctn::AliHBTTwoTrackEffFctn()
{
  //def ctor
}
/******************************************************************/

AliHBTTwoTrackEffFctn::AliHBTTwoTrackEffFctn(Int_t nbins, Double_t maxval, Double_t minval):
     AliHBTOnePairFctn1D("TwoTrackEff","Two Track Efficiency",nbins,maxval,minval)
{
//contructor
//nbins - numner of bins of the function
//maxval - max X of the fctn
//minval - min X of the fctn
 GetNumerator()->GetXaxis()->SetTitle("dP [GeV]");
 GetDenominator()->GetXaxis()->SetTitle("dP [GeV]");
}
/******************************************************************/

TH1* AliHBTTwoTrackEffFctn::GetResult()
{
//returns ratio of numerator and denominator
 delete fRatio;
 fRatio = GetRatio(Scale());
 if(fRatio)
  {
   fRatio->GetXaxis()->SetTitle("dP [GeV]");
   fRatio->GetYaxis()->SetTitle("C(dP)");
   fRatio->SetTitle("Double Track Resolution: dP Correlation Fctn.");
  }
 return fRatio;
}
/******************************************************************/
/******************************************************************/
/******************************************************************/
ClassImp(AliHBTTwoTrackEffFctnPxPyPz)

AliHBTTwoTrackEffFctnPxPyPz::AliHBTTwoTrackEffFctnPxPyPz(Int_t nXbins, Double_t maxXval, Double_t minXval,
                                                   Int_t nYbins, Double_t maxYval, Double_t minYval,
                                                   Int_t nZbins, Double_t maxZval, Double_t minZval):
 AliHBTOnePairFctn3D(nXbins,maxXval,minXval,nYbins,maxYval,minYval,nZbins,maxZval,minZval)
{
//ctor
//Set Axis Title
 fWriteNumAndDen = kTRUE;
 Rename("tteffpxpypz","P_{x} P_{y} P_{z} Two Track Efficiency Function");
 if(fNumerator)
  {
   fNumerator->GetXaxis()->SetTitle("\\Delta P_{x} [GeV]");
   fNumerator->GetYaxis()->SetTitle("\\Delta P_{y} [GeV]");
   fNumerator->GetZaxis()->SetTitle("\\Delta P_{z} [GeV]");
  }

 if(fDenominator)
  {
   fDenominator->GetXaxis()->SetTitle("\\Delta P_{x} [GeV]");
   fDenominator->GetYaxis()->SetTitle("\\Delta P_{y} [GeV]");
   fDenominator->GetZaxis()->SetTitle("\\Delta P_{z} [GeV]");
  }

}
/******************************************************************/

void AliHBTTwoTrackEffFctnPxPyPz::GetValues(AliHBTPair* pair, Double_t& x, Double_t&y ,Double_t& z) const
{
//Returns values to be histogrammed
//it does not 
 x = pair->GetDeltaPx();
 y = pair->GetDeltaPy();
 z = pair->GetDeltaPz();
}
/******************************************************************/

TH1* AliHBTTwoTrackEffFctnPxPyPz::GetResult()
{
//returns ratio of numerator and denominator
 delete fRatio;
 fRatio = GetRatio(Scale());
 return fRatio;
}

/******************************************************************/
/******************************************************************/
/******************************************************************/
ClassImp(AliHBTTwoTrackEffFctnPtThetaPhi)

AliHBTTwoTrackEffFctnPtThetaPhi::AliHBTTwoTrackEffFctnPtThetaPhi(Int_t nXbins, Double_t maxXval, Double_t minXval,
                                                   Int_t nYbins, Double_t maxYval, Double_t minYval,
                                                   Int_t nZbins, Double_t maxZval, Double_t minZval):
 AliHBTOnePairFctn3D(nXbins,maxXval,minXval,nYbins,maxYval,minYval,nZbins,maxZval,minZval)
{
//ctor
//Set Axis Title
 fWriteNumAndDen = kTRUE;
 Rename("tteffptthetaphi","P_{t} \\theta \\phi Two Track Efficiency Function");
 if(fNumerator)
  {
   fNumerator->GetXaxis()->SetTitle("\\Delta P_{t} [GeV]");
   fNumerator->GetYaxis()->SetTitle("\\Delta \\theta [rad]");
   fNumerator->GetZaxis()->SetTitle("\\Delta \\phi [rad]");
  }

 if(fDenominator)
  {
   fDenominator->GetXaxis()->SetTitle("\\Delta P_{t} [GeV]");
   fDenominator->GetYaxis()->SetTitle("\\Delta \\theta [rad]");
   fDenominator->GetZaxis()->SetTitle("\\Delta \\phi [rad]");
  }
}
/******************************************************************/

void AliHBTTwoTrackEffFctnPtThetaPhi::GetValues(AliHBTPair* pair, Double_t& x, Double_t&y ,Double_t& z) const
{
//Returns values to be histogrammed
//it does not 
 x = pair->GetDeltaPt();
 y = pair->GetDeltaTheta();
 z = pair->GetDeltaPhi();
}
/******************************************************************/

TH1* AliHBTTwoTrackEffFctnPtThetaPhi::GetResult()
{
//returns ratio of numerator and denominator
 delete fRatio;
 fRatio = GetRatio(Scale());
 if(fRatio)
  {
   fRatio->GetXaxis()->SetTitle("\\Delta P_{t} [GeV]");
   fRatio->GetYaxis()->SetTitle("\\Delta \\theta [rad]");
   fRatio->GetZaxis()->SetTitle("\\Delta \\phi [rad]");
  }
 return fRatio;
}
