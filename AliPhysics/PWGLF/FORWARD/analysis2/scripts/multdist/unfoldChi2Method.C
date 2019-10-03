#ifndef __CINT__
#include <iostream>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TList.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TPad.h>
#include <THistPainter.h>
#include <TObject.h>
#include <TMath.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TLatex.h>
#include <TLine.h>
#include <TMarker.h>
#include <TStyle.h>
#include <TVirtualFitter.h>
#include "/home/caz/ALICE/AliRoot/ANALYSIS/AliUnfolding.h"
#include "/home/caz/ALICE/AliRoot/ANALYSIS/AliUnfolding.cxx"
#include <AliPWG0Helper.h>
// #include <AliMultiplicityCorrection.h>
//#include <TMinuit.h>
#include <TBox.h>
#include <TGaxis.h>
#else
class TF1;
class TH1D;
class TH1F;
class TGraphErrors;
class TGaxis;
// class AliUnfolding;
#endif
using namespace std;


void UnfoldChi2Min(TH1F*   data, 
		   TH1F*&  unfolded, 
		   TH2D*   response=0 ,
		   TH1F*   eff =0,
		   Int_t   regularization=1,
		   Float_t   regularizationWeight=1000,
		   Float_t minInitialValue=0.0001, 
		   Bool_t  skipBin0InChi2= kFALSE,
		   Int_t   skipNbins=0 ){
#ifdef __CINT__
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libPWG0base");
  gSystem->Load("libPWG0dep");
  gSystem->Load("libPWG0selectors");
#endif

  cout << data << "  " << unfolded << endl;
  
  AliUnfolding::RegularizationType type = AliUnfolding::kNone;
  switch (regularization) {
  case 0: type = AliUnfolding::kNone; break; 
  case 1: type = AliUnfolding::kPol0; break; 
  case 2: type = AliUnfolding::kPol1; break; 
  case 3: type = AliUnfolding::kLog; break; 
  case 4: type = AliUnfolding::kEntropy; break; 
  case 5: type = AliUnfolding::kCurvature; break; 
  case 6: type = AliUnfolding::kRatio; break; 
  };

  AliUnfolding::SetUnfoldingMethod(AliUnfolding::kChi2Minimization);
  AliUnfolding::SetNbins(data->GetNbinsX(), unfolded->GetNbinsX()); 
  AliUnfolding::SetMinimumInitialValue(kTRUE, minInitialValue);
  AliUnfolding::SetSkip0BinInChi2(skipBin0InChi2);
  AliUnfolding::SetSkipBinsBegin(skipNbins);
  AliUnfolding::SetChi2Regularization(type, regularizationWeight);
  //AliUnfolding::SetCreateOverflowBin(1e-6);
    
  AliUnfolding::Unfold(response, eff, data, data, unfolded, kFALSE);
 
 
}
