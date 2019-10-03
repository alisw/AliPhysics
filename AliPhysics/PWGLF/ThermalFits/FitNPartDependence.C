#if !defined (__CINT__) || (defined(__MAKECINT__))
#include <iostream>
#include <map>
#include <fstream>
#include "TString.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "AliParticleYield.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TROOT.h"
#include "TVirtualFitter.h"
#include "TMath.h"
#include "TFitResult.h"
#include "TMinuit.h"
#include "TLatex.h"
#include "TDatabasePDG.h"
#include "TPad.h"
#include "TCanvas.h"

#endif




std::map<TString,Float_t> npartPbPb;
std::map<TString,Float_t> npartPbPbErr;

void ReadCentralityFromFile() ;
Double_t FitShiftedGraphAndExtrapolate(TGraphErrors * gr, Int_t shift, TF1 * f1, const char * centr, Color_t color,  Bool_t draw = kFALSE) ;
Double_t ErrorFunction (Double_t *xx, Double_t *p ) ;
Double_t FuncPlusErr (Double_t *xx, Double_t *p) ;

enum {kNoShift, kShiftUp, kShiftDown, kShiftHarder, kShiftSofter};
enum {kCorrelated, kAntiCorrelated};

TF1 * fForErrors = 0;
const int npar = 3;
Double_t matrix[npar][npar];
TString centrFile;
TString systemAndEnergy;
Double_t maxy = 0;
TString label = "#frac{d#it{N}}{d#it{y}} = #it{a} + #it{b} #times (#it{N}_{part})^{#it{c}}";
//TString label = "#frac{d#it{N}}{d#it{y}} = #it{a} + #it{b} #times (d#it{N}/d#it{#eta})^{#it{c}}";
Int_t collSystem  = 2;
Float_t energy = 2760;

void FitNPartDependence() {

  // This macro can be used to fit the centrality dependence of
  // different particles in order to extrapolate or interpolate the
  // yield as a function of centrality.
  // The value of npart or dndeta is read from a text file, which must
  // be in the format "centrTag value error"

  // The particle yields are taken from the machine readable files.
  // If you want to fit a different particle, you have to uncomment
  // the corresponding section below.  if you want to try a different
  // function, scroll down to the definition of f1

  // For the estimate of the systematic uncertainties, the syst graph
  // is shifted up and down within its uncertainties. Few strategies
  // are possible (see the call to FitShiftedGraphAndExtrapolate below
  // to change the behavior):
  // 1. All points are coherently shifted up and down by 1 syst bar
  // (kShiftUp, kShiftDown)
  // 2. A worst case scenario is implemented, where the first point is
  // shifted by +- 1 sigma, the last one by -+ 1 sigma and all the
  // other points are shifted smoothly in between (kShiftHarder,
  // kShiftSofter)

  //__________________________________________________________________//

  // WARNING: check isSum

  // KStar
  // centrFile = "npart_PbPb.txt";
  // //  centrFile = "dndeta_PbPb.txt";
  // maxy = 50;
  // systemAndEnergy = "Pb-Pb #sqrt{#it{s}}_{NN} = 2.76 TeV";
  // const char * centralityToPlot[] = {   "V0M0020" ,  "V0M2040" ,  "V0M4060" ,  "V0M6080" , 0};
  // const char * centrToExtrapolate = "V0M0010";
  //    Int_t errorType = kAntiCorrelated;

  // Int_t pdg = 313;
  // TClonesArray * arr = AliParticleYield::ReadFromASCIIFile("./PbPb_2760_Kstar892.txt");  
  // Deuteron
  // gROOT->ProcessLine(".x figTemplate.C(0,0.00001,400,0.2)");  
  // const char * centralityToPlot[] = {   "V0M0010" ,  "V0M1020" ,  "V0M2040" ,  "V0M4060" ,  "V0M6080",0 };
  // const char * centrToExtrapolate = "V0M0010";
  // Int_t pdg = 1000010020;
  // TClonesArray * arr = AliParticleYield::ReadFromASCIIFile("PbPb_2760_DeuHelium3.txt");
  //  Deuteron pPb
  // centrFile = "dndeta_pPb.txt";
  // const char * centralityToPlot[] = {   "V0A0010", "V0A1020", "V0A2040", "V0A4060", "V0A6000" ,0};
  // const char * centrToExtrapolate = "V0A0005";
  // //const char * centrToExtrapolate = "V0A6080";
  // Int_t pdg = -1000010020;
  // TClonesArray * arr = AliParticleYield::ReadFromASCIIFile("pPb_5020_deuteron.txt");
  // maxy = 0.01;
  // systemAndEnergy = "p-Pb #sqrt{#it{s}}_{NN} = 5.02 TeV";
  // energy = 5020;
  // collSystem = 1;
  // K* pPb
  // centrFile =  "dndeta_pPb.txt";
  // //  centrFile = "dndeta_PbPb.txt";
  // maxy = 1;
  // systemAndEnergy = "p-Pb #sqrt{#it{s}}_{NN} = 5.02 TeV";
  // const char * centralityToPlot[] = {   "V0A0020" ,  "V0A2040" ,  "V0A4060" ,  "V0A6080", "V0A8000" , 0};
  // const char * centrToExtrapolate = "V0A0005";
  // Int_t pdg = 313;
  // TClonesArray * arr = AliParticleYield::ReadFromASCIIFile("./pPb_5020_Kstar.txt");  
  // energy = 5020;
  // collSystem = 1;


  // Helium3
  centrFile = "npart_PbPb.txt";
  // //  centrFile = "dndeta_PbPb.txt";

  gROOT->ProcessLine(".x figTemplate.C(0,0.00001,400,0.2)");  
  const char * centralityToPlot[] = {   "V0M0020" ,  "V0M2080" ,0};
  const char * centrToExtrapolate = "V0M1020"; // "V0M0010";
  Int_t pdg = 1000020030;
  TClonesArray * arr = AliParticleYield::ReadFromASCIIFile("PbPb_2760_DeuHelium3.txt");
  maxy=1e-3;
  Int_t errorType = kCorrelated;
  

  //  gPad->GetCanvas()->SetTitle(Form("c%s", TDatabasePDG::Instance()->GetParticle(pdg)->GetName()));

  ReadCentralityFromFile();
  TGraphErrors * grSyst = new TGraphErrors; 
  TGraphErrors * grStat = new TGraphErrors; 
  Double_t maxx = 1.1*npartPbPb[centrToExtrapolate];
  //  Double_t maxx = 1.1*npartPbPb["V0A0005"];
  // Function
  // TF1 * f1 = new TF1 ("f1", "[0] + [1]*x", 0, maxx);
  // f1->SetParameters(1,1);
  TF1 * f1 = new TF1 ("f1", "[0] + [1]*x^[2]", 0, maxx);
  f1->SetParameters(0, 6.3266e-04, 8.99883e-01);
  f1->FixParameter(0,0);
  maxx=400;
  std::cout << "Centrality "<< maxx << " " << maxy << std::endl;
  
  gROOT->ProcessLine(Form(".x figTemplate.C(0,0.0000001,%f,%f)", maxx, maxy));
  gPad->GetCanvas()->SetTitle("cKstar");

  //  f1->FixParameter(0, 0);
  Int_t icentr = 0;
  AliParticleYield * part = 0;
  while (centralityToPlot[icentr]) {
    part =  AliParticleYield::FindParticle (arr, pdg, collSystem, energy, centralityToPlot[icentr]);
    if(part) {
      grSyst->SetPoint     (icentr , npartPbPb[centralityToPlot[icentr]]    , part->GetYield());
      grSyst->SetPointError(icentr , npartPbPbErr[centralityToPlot[icentr]] , part->GetSystError());
      grStat->SetPoint     (icentr , npartPbPb[centralityToPlot[icentr]]    , part->GetYield());
      grStat->SetPointError(icentr , npartPbPbErr[centralityToPlot[icentr]] , part->GetStatError());
      
      //      part->Print();
    }
    icentr++;
  }
  grStat->Draw("P");
  grStat->SetMarkerStyle(24);
  grSyst->Draw("PE2");
  grSyst->SetFillStyle(0);
  Double_t yield          = FitShiftedGraphAndExtrapolate(grStat, kNoShift  , f1, centrToExtrapolate, kRed  );

  // We need to cache the covariance matrix and the function here, otherwise things get messed up
  fForErrors = (TF1*) f1->Clone();
  fForErrors->SetLineColor(kRed);
  fForErrors->SetLineStyle(0);
  fForErrors->Draw("same");
  gMinuit->mnemat(&matrix[0][0],npar);
  // The statistical error is computed propagating the uncertainty on the parameters to the function via the covariance matrix
  TF1 * fShiftedPlus   = new TF1("fShiftedPlus"  , FuncPlusErr, 0, maxx,1); fShiftedPlus->SetParameter(0,+1);
  TF1 * fShiftedMinus  = new TF1("fShiftedMinus" , FuncPlusErr, 0, maxx,1); fShiftedMinus->SetParameter(0,-1);
  fShiftedMinus->SetLineStyle(kDotted);
  fShiftedPlus->SetLineStyle (kDotted);
  fShiftedPlus ->Draw("same");
  fShiftedMinus->Draw("same");
  TF1 * fError = new TF1("fError", ErrorFunction, 0,maxx, 0);

  // The uncertainty on the systematics is computed shifting the graph up and down + refitting
  Double_t errorSystPlus  = 0;
  Double_t errorSystMinus = 0;
  if (errorType == kCorrelated) {
    errorSystPlus  = TMath::Abs(FitShiftedGraphAndExtrapolate(grSyst, kShiftUp  , f1, centrToExtrapolate, kRed)-yield);
    errorSystMinus = TMath::Abs(FitShiftedGraphAndExtrapolate(grSyst, kShiftDown, f1, centrToExtrapolate, kRed)-yield);
  }
  else if (errorType == kAntiCorrelated) {
    errorSystPlus  = TMath::Abs(FitShiftedGraphAndExtrapolate(grSyst, kShiftHarder, f1, centrToExtrapolate, kRed)-yield);
    errorSystMinus = TMath::Abs(FitShiftedGraphAndExtrapolate(grSyst, kShiftSofter, f1, centrToExtrapolate, kRed)-yield);
  }
  else {
    std::cout << "Wrong error type" << std::endl;
    return;
  }

  Double_t errorStat = fError->Eval(npartPbPb[centrToExtrapolate]);


  std::cout << "Yield from fit: ("<<centrToExtrapolate<< "="<<npartPbPb[centrToExtrapolate]<<")" 
            << yield 
            << "   +" << errorStat << " " 
            << "   +" << errorSystPlus << " " << errorSystMinus//<< std::endl;
            << "   +" << fError->Eval(npartPbPb[centrToExtrapolate])  << std::endl;
  std::cout << yield << "    " 
            << errorStat << "     " 
            << (errorSystPlus+errorSystMinus)/2 << std::endl;
  

  TGraphErrors * gExtrap = new TGraphErrors();
  gExtrap->SetMarkerStyle(20);
  gExtrap->SetPoint(0, npartPbPb[centrToExtrapolate], yield);
  gExtrap->SetPointError(0, 0, (errorSystPlus+errorSystMinus)/2);
  gExtrap->Draw("P");

  TLatex * text = new TLatex (0.2,0.81,systemAndEnergy);
  text->SetNDC();
  text->Draw();
  TLatex * text2 = new TLatex (0.2,0.72,part->GetLatexName());
  text2->SetNDC();
  text2->Draw(); 

  TLatex * func = new TLatex(0.53,0.23,label);
  func->SetNDC();
  func->Draw();


}

Double_t ErrorFunction (Double_t *xx, Double_t *p ) {



  //  Double_t x = xx[0];
  Double_t func[npar];
  // func[0] = 1;
  // func[1] = TMath::Power(x, fForErrors->GetParameter(2));
  // func[2] = fForErrors->GetParameter(1)*TMath::Power(x, fForErrors->GetParameter(2))*TMath::Log(x);
  //  In general, one can compute the derivative numerically using root.
  for(Int_t ipar = 0; ipar < npar; ipar++){
        func[ipar] =  fForErrors->GradientPar(ipar, xx);
  }

  

  Double_t variance = 0;
  
  for(Int_t ipar = 0; ipar < npar; ipar++){
    for(Int_t jpar = 0; jpar < npar; jpar++){
      variance = variance + func[ipar]*func[jpar] * matrix[ipar][jpar];
    }
    
  }
  Double_t error = TMath::Sqrt(variance);
  return error;
}

Double_t FuncPlusErr (Double_t *xx, Double_t *p) {

  Double_t value = fForErrors->Eval(xx[0]);
  Double_t error = ErrorFunction(xx, p);

  if(p[0] > 0) {
    return value + error;
  }
  return value -error;
}


Double_t FitShiftedGraphAndExtrapolate(TGraphErrors * gr, Int_t shift, TF1 * f1, const char * centr, Color_t color, Bool_t draw) {
  // Shift graph
  // 1 = up
  // 2 = down

  TGraphErrors * grLocal = (TGraphErrors*) gr->Clone();
  Int_t npoint = grLocal->GetN() ;
  if(shift != kNoShift) {
    for(Int_t ipoint = 0; ipoint < npoint; ipoint++){
      Double_t relError = grLocal->GetEY()[ipoint]/grLocal->GetY()[ipoint];
      Double_t value    = grLocal->GetY() [ipoint];
      if(shift == kShiftUp  ) value += grLocal->GetEY()[ipoint];        
      if(shift == kShiftDown) value -= grLocal->GetEY()[ipoint];   
      if(shift == kShiftHarder) {        
        Double_t valueShift = -1. + 2.*Float_t(ipoint)/(npoint-1);
        value += valueShift*grLocal->GetEY()[ipoint];
        std::cout << "valueShift " << valueShift << std::endl;
      }
      if(shift == kShiftSofter) {        
        Double_t valueShift = +1. - 2.*Float_t(ipoint)/(npoint-1);
        value += valueShift*grLocal->GetEY()[ipoint];
        std::cout << "valueShift " << valueShift << std::endl;
      }
      grLocal->SetPoint     (ipoint, grLocal->GetX() [ipoint], value);
      grLocal->SetPointError(ipoint, grLocal->GetEX()[ipoint], relError*value);
    }
  }
  grLocal->Fit(f1, "", "Q");


  Double_t yield = f1->Eval(npartPbPb[centr]);
  TF1 * clone = (TF1*) f1->Clone();
  clone->SetLineWidth(1);
  clone->SetLineStyle(kDashed);
  clone->SetLineColor(color);
  clone->Draw("Same");
  if(draw) {
    TVirtualPad * oldPad = gPad;
    new TCanvas(Form("shifted_%d", shift), Form("shifted_%d", shift));
    grLocal->Draw("AP");
    clone->Draw("same");
    oldPad->cd();
  } else {
    delete grLocal;
  }
  return yield;

}

void ReadCentralityFromFile() {
  ifstream filein (centrFile);
  TString line;
  while (line.ReadLine((filein))) {    
    TObjArray * tokens = line.Tokenize(" ");
    if(tokens->GetEntries() != 3) {
      delete tokens;
      continue;
    }
    TString str        = ((TObjString*)tokens->At(0))->String();
    Float_t npart     = ((TObjString*)tokens->At(1))->String().Atof();
    Float_t npartErr  = ((TObjString*)tokens->At(2))->String().Atof();
    //    std::cout << "["<<str.Data() <<"]" << npart << " " << npartErr << std::endl;
    if(npart) {
      npartPbPb   [str] = npart;
      npartPbPbErr[str]  = npartErr;
    }
    delete tokens;
  }

}
