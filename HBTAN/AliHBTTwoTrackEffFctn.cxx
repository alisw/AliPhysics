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
ClassImp(AliHBTTwoTrackEffFctn3D)

AliHBTTwoTrackEffFctn3D::AliHBTTwoTrackEffFctn3D()
{
//Set Axis Title
}
/******************************************************************/

void AliHBTTwoTrackEffFctn3D::GetValues(AliHBTPair* pair, Double_t& x, Double_t&y ,Double_t& z)
{
//Returns values to be histogrammed
//it does not 
 x = pair->GetDeltaPx();
 y = pair->GetDeltaPy();
 z = pair->GetDeltaPz();
}
