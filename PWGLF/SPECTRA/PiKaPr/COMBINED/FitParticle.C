#if !defined(__CINT__) || defined(__MAKECINT__)
#include <iostream>

#include "AliPWGHistoTools.h"
#include "AliPWGFunc.h"
#include "AliLatexTable.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1.h"
#include "TROOT.h"
#include "TDatabasePDG.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TSystem.h"
#include "THashList.h"
#include "TMinuit.h"
#include "TLatex.h"

using namespace std;
#endif

enum {kFitExpPt, kFitLevi, fFitExpMt, kFitBoltzmann, kFitBlastWave, kFitBoseEinstein, kFitFermiDirac};
Bool_t skipMean = 0;

void FitParticle(TH1 * h, const char * partName, Float_t min = 0, Float_t max =3, Float_t scaleHisto = -1., Int_t fitFunc = kFitLevi, Int_t vartype = AliPWGFunc::kdNdpt, const char * fileOut = 0, Bool_t wait = 0, Float_t meanMin = 0., Float_t meanMax = 100) ;
void FitParticle(const char * file, const char * histo, const char * partName,  const char * listname=0, Float_t min = 0, Float_t max =3, Float_t scaleHisto = -1., Int_t fitFunc = kFitLevi, Int_t vartype = AliPWGFunc::kdNdpt, const char * fileOut = 0, Bool_t wait = 0);

void FitParticle(const char * file, const char * histo, const char * partName,  const char * listname, Float_t min, Float_t max, Float_t scaleHisto, Int_t fitFunc, Int_t vartype, const char * fileOut, Bool_t wait) {

  // Generic Macro to fit any particle using the PWG2/SPECTRA/Fit macros
  //
  // You have to provide:
  // - file: file name. If not set, uses gFile
  // - histo: histogram name name (assumend to be in the main folder of the file)
  // - listname: if different from 0, histo is looked in this tlist
  // - partName: it is used to get the mass of the particle from
  //   TDatabasePDG. If you are not sure, just use a part of the name: a
  //   list of names matching it is printed on screen
  // 
  // You can optionally provide:
  // - min, max: the fit range
  // - scaleHisto: a scaling factor for the histo. If negative, it is
  //   ignored. If the histo is scaled, the bin width is also
  //   divided).
  // - fitFunc: id of the function, levi is the default
  //   valide options: kFitExpPt, kFitLevi, fFitExpMt, kFitBoltzmann, kFitBlastWave
  // - varType: the variable used in the pt spectrum (see AliPWGFunc.h)

  // Author: Michele Floris, CERN

  // load stuff
  gSystem->Load("libTree.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libMinuit.so");
  gSystem->Load("libSTEERBase.so");
  gSystem->Load("libESD.so");
  gSystem->Load("libAOD.so");
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  gSystem->Load("libCORRFW.so");
  gSystem->Load("libPWG2spectra.so");


  // open file
  TFile * f = file ? new TFile(file) : gFile;  
  TH1 * h = 0;
  if(listname){
    TList * l = (TList*) gDirectory->Get(listname);
    h = (TH1*) l->FindObject(histo);
  }
  else{
    h = (TH1*) gDirectory->Get(histo); // 
  }

  FitParticle(h, partName, min, max, scaleHisto, fitFunc, vartype, fileOut, wait);
  


}

void FitParticle(TH1 * h, const char * partName, Float_t min , Float_t max, Float_t scaleHisto, Int_t fitFunc, Int_t vartype, const char * fileOut, Bool_t wait, Float_t meanMin, Float_t meanMax) { 


  // get histo and draw
  AliPWGFunc * fm = new AliPWGFunc;
  fm->SetVarType(AliPWGFunc::VarType_t(vartype));
  //  fm->SetVarType(AliPWGFunc::VarType_t(0));//FIXME
  //  cout << "Warning: hacked vartype" << endl;
  
  if (!TDatabasePDG::Instance()->GetParticle(partName)) {
    cout << "Wrong particle name " << partName << endl;

    const THashList * l = TDatabasePDG::Instance()->ParticleList();
    Int_t npart = l->GetSize();
    for(Int_t ipart = 0; ipart < npart; ipart++){
      TString name = l->At(ipart)->GetName();
      if(name.Contains(partName, TString::kIgnoreCase))
	cout << " - Did you mean [" << name.Data() << "] ?" << endl;		 
    }
    
    cout << "Try [ TDatabasePDG::Instance()->GetParticle(partName) ] by hand "<< endl;
    return;
    
  }
  Double_t mass = TDatabasePDG::Instance()->GetParticle(partName)->Mass();
  cout << "Fitting ["<<partName<<"], mass: " << mass << endl;
  

  AliLatexTable * table;

  TF1* func = 0;
  Int_t normPar = -1;
  Int_t slopePar=0;
  if (fitFunc == kFitLevi)   {
    func = fm->GetLevi (mass, 0.4, 750,3);
    func->SetParLimits(1,0.0001,20000);
    table = new AliLatexTable(10,"c|ccccccccc");
    table->InsertCustomRow("Part & Func & Yield & FD Slope & $\\Chi^2$/NDF & Min X & Fit Min & Fit Max & Frac Below  & Mean ");
    normPar = 0;
    slopePar = 2;
  }
  if (fitFunc == kFitExpPt)  {
    func = fm->GetPTExp(0.2, 20);
    table = new AliLatexTable(10,"c|ccccccccc");
    table->InsertCustomRow("Part & Func & Yield & pT Slope & $\\Chi^2$/NDF & Min X & Fit Min & Fit Max & Frac Below & Mean ");
    normPar = 0;
    slopePar = 1;
  }
  if (fitFunc == fFitExpMt)  {
    func = fm->GetMTExp(mass,0.2, 20);
    table = new AliLatexTable(10,"c|ccccccccc");
    table->InsertCustomRow("Part & Func & Yield & mT Slope & $\\Chi^2$/NDF & Min X & Fit Min & Fit Max & Frac Below & Mean  ");
    normPar = 0;
    slopePar = 1;
  }
  if (fitFunc == kFitBoltzmann)  {
    func = fm->GetBoltzmann(mass, 0.2, 20);
    table = new AliLatexTable(10,"c|ccccccccc");
    table->InsertCustomRow("Part & Func & Yield & Bz Slope & $\\Chi^2$/NDF & Min X & Fit Min & Fit Max & Frac Below & Mean  ");
    normPar = 0;
    slopePar = 1;

  }
  if (fitFunc == kFitBlastWave)  {
    func = fm->GetBGBW(mass,0.6,0.3, 1, 1e5);// beta, T, n, norm 
    table = new AliLatexTable(12,"cc|cccccccccccc");
    table->InsertCustomRow("Part & Func & Yield & T & beta & n & $\\Chi^2$/NDF & Min X & Fit Min & Fit Max & Frac Below  & Mean ");
    func->SetParLimits(1, 0.1, 0.99);
    func->SetParLimits(2, 0.01, 1);
    func->SetParLimits(3, 0.01, 2);
    normPar = 4;
    slopePar = 2;
  }
  if (fitFunc == kFitBoseEinstein)  {
    func = fm->GetBoseEinstein(mass,0.3,20);
    table = new AliLatexTable(10,"c|ccccccccc");
    table->InsertCustomRow("Part & Func & Yield & BE Slope & $\\Chi^2$/NDF & Min X & Fit Min & Fit Max & Frac Below & Mean ");
    normPar = 0;
    slopePar = 1;
  }
  if (fitFunc == kFitFermiDirac)  {
    func = fm->GetFermiDirac(mass,0.3,20);
    table = new AliLatexTable(10,"c|ccccccccc");
    table->InsertCustomRow("Part & Func & Yield & FD Slope & $\\Chi^2$/NDF & Min X & Fit Min & Fit Max & Frac Below  & Mean  ");
    normPar = 0;
    slopePar = 1;
  }

    h->Sumw2();
  if (scaleHisto > 0) h->Scale(scaleHisto, "width");
//   TH1 * h = AliPWGHistoTools::GetdNdPtFromOneOverPt((TH1*) gDirectory->Get(histo)); // FIXME
//   cout << "WARNING SCALING2PI" << endl;
//   h->Scale(2*TMath::Pi());//Fixme

   h->Fit(func); // FIXME
   h->Fit(func); // FIXME
   h->Fit(func); // FIXME
   h->Fit(func,"IME0","",min,max);      
   h->Fit(func,"IME0","",min,max);      

   //   gMinuit->Command("SET STRATEGY 2"); // FIXME

  // if (!AliPWGHistoTools::Fit(h,func,min,max)) {
  //   cout << "Fitting error!" << endl;
  //   //    return;    
  // } FIXME
  cout << "Drawing" << endl;

  TCanvas * c = new TCanvas();
  cout << "1" << endl;
  c->SetLogy();
  cout << "2" << endl;
  h->Draw();
  cout << "3" << endl;
  func->Draw("same");
  cout << "4" << endl;

  // Print results nicely
  // populate table

  Double_t yield=0,yieldE=0;
  cout << "Y" << endl;
  AliPWGHistoTools::GetYield(h,func,yield,yieldE);  
  cout << "YE" << endl;
  TLatex * l = new TLatex(2,(h->GetMaximum()+h->GetMinimum())/2,Form("%3.3f #pm %3.3f (%s)", yield,yieldE, func->GetName()));
  l->Draw();

//   Float_t yield  = func->Integral(0.45,1.05);
//   Float_t yieldE = func->IntegralError(0.45,1.05);
  cout << "Slope Par: "<< slopePar << endl;

  Double_t tslope   = func->GetParameter(slopePar);
  Double_t tslopeE = func->GetParError(slopePar);	

  table->SetNextCol(Form("%s (%s)", partName, h->GetName()));
  table->SetNextCol(func->GetName());
  table->SetNextCol(yield,yieldE,-4);
  table->SetNextCol(tslope,tslopeE,-4);
  if(fitFunc == kFitBlastWave) {
    table->SetNextCol(func->GetParameter(1),func->GetParError(1),-4);
    table->SetNextCol(func->GetParameter(3),func->GetParError(3),-4);
  }
  //  table->SetNextCol(func->GetParameter(1),func->GetParError(1),-4);
  table->SetNextCol(Form("%2.2f/%d",func->GetChisquare(),func->GetNDF()));
  Float_t lowestPoint = TMath::Max(AliPWGHistoTools::GetLowestNotEmptyBinEdge(h),min);
  //  Float_t yieldAbove  = func->Integral(lowestPoint,100);
  Float_t yieldBelow  = func->Integral(0,lowestPoint);
  table->SetNextCol(lowestPoint,-3);
  table->SetNextCol(min,-2);
  table->SetNextCol(max,-2);
  table->SetNextCol(yieldBelow/yield,-3);

  Double_t mean=0, meane=0;
  // Float_t mean2=0, mean2e=0;
  //  AliPWGHistoTools::GetMean      (func, mean,  meane , 0.,100., normPar);
  AliPWGHistoTools::GetMeanDataAndExtrapolation      (h, func, mean,  meane , meanMin, meanMax);
  //  AliPWGHistoTools::GetMeanDataAndExtrapolation      (h, func, mean,  meane , 0.,4.5);
  // AliPWGHistoTools::GetMeanSquare(func, mean2, mean2e, 0.,100., normPar);
  if(skipMean) table->SetNextCol("N/A");
  else table->SetNextCol(mean,  meane ,-4);
  // table->SetNextCol(mean2, mean2e,-4);
  //			 fMean2->IntegralError(0,100)/func->Integral(0,100),-7);
  table->InsertRow();
  table->PrintTable("ASCII");
  
  if(fileOut) {
    TFile * fout = new TFile(fileOut, "update");    
    c->SetName(Form("c_%s_%s_fit_%2.2f_%2.2f", h->GetName(), func->GetName(), min, max));
    c->Write();
    TH1F * hFit = new TH1F ("hFitResults", "hFitResults", 2,0,1);
    hFit->SetBinContent(1, yield); hFit->SetBinError(1, yieldE);
    hFit->SetBinContent(2, mean ); hFit->SetBinError(2, meane);
    hFit->GetXaxis()->SetBinLabel(1, "yield");
    hFit->GetXaxis()->SetBinLabel(2, "<pt>");
    hFit->Write(Form("hFit_%s_%s_fit_%2.2f_%2.2f", h->GetName(), func->GetName(), min, max));
    fout->Purge(); // remove old canvases
    fout->Close();
  }


  if(wait)  c->WaitPrimitive();


}

// Double_t MicheleFunction(Double_t *x, Double_t *p) {

//   Double_t xloc = x[0];
//   Double_t y = 0;
//   if(xloc < p[0]) {
//     y = xloc*TMath::Exp(-xloc/p[2]);
//   } else if (xloc > p[0] && xloc < p[1]) {
    
//   }


// }
