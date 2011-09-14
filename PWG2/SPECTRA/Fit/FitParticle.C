#if !defined(__CINT__) || defined(__MAKECINT__)
#include <iostream>

#include "AliBWTools.h"
#include "AliBWFunc.h"
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

using namespace std;
#endif

enum {kFitExpPt, kFitLevi, fFitExpMt, kFitBoltzmann, kFitBlastWave};

void FitParticle(const char * file, const char * histo, const char * partName,  const char * listname=0, Float_t min = 0, Float_t max =3, Float_t scaleHisto = -1., Int_t fitFunc = kFitLevi, Int_t vartype = AliBWFunc::kdNdpt) {

  // Generic Macro to fit any particle using the PWG2/SPECTRA/Fit macros
  //
  // You have to provide:
  // - file: file name 
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
  // - varType: the variable used in the pt spectrum (see AliBWFunc.h)

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


  // get histo and draw
  AliBWFunc * fm = new AliBWFunc;
  fm->SetVarType(AliBWFunc::VarType_t(vartype));
  //  fm->SetVarType(AliBWFunc::VarType_t(0));//FIXME
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
  
  TF1* func = 0;
  Int_t normPar = -1;
  if (fitFunc == kFitLevi)   {
    func = fm->GetLevi (mass, 0.4, 20,3);
    func->SetParLimits(1,1,100);
    normPar = 0;
  }
  if (fitFunc == kFitExpPt)  {
    func = fm->GetPTExp(0.2, 20);
  }
  if (fitFunc == fFitExpMt)  {
    func = fm->GetMTExp(mass,0.2, 20);
  }
  if (fitFunc == kFitBoltzmann)  {
    func = fm->GetBoltzmann(mass, 0.2, 20);
  }
  if (fitFunc == kFitBlastWave)  {
    func = fm->GetBGBW(mass,0.6,0.3, 20);
  }


  TFile * f = new TFile(file);  
  TH1 * h = 0;
  if(listname){
    TList * l = (TList*) gDirectory->Get(listname);
    h = (TH1*) l->FindObject(histo);
  }
  else{
    h = (TH1*) gDirectory->Get(histo); // 
  }
  h->Sumw2();
  if (scaleHisto > 0) h->Scale(scaleHisto, "width");
//   TH1 * h = AliBWTools::GetdNdPtFromOneOverPt((TH1*) gDirectory->Get(histo)); // FIXME
//   cout << "WARNING SCALING2PI" << endl;
//   h->Scale(2*TMath::Pi());//Fixme

   h->Fit(func); // FIXME
   gMinuit->Command("SET STRATEGY 2"); // FIXME

  if (!AliBWTools::Fit(h,func,min,max)) {
    cout << "Fitting error!" << endl;
    //    return;    
  }
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
  AliLatexTable table(9,"c|ccccccc");
  table.InsertCustomRow("Part & Integral & T Slope & n & $\\Chi^2$/NDF & Min X & Frac Above & \\langle p_{t} \\rangle  & \\langle p_{t}^{2} \\rangle");
  // populate table
  // Float_t yield  = func->Integral(0,100);
  // Float_t yieldE = func->IntegralError(0,100);
  Double_t yield=0,yieldE=0;
  cout << "Y" << endl;
  AliBWTools::GetYield(h,func,yield,yieldE);  
  cout << "YE" << endl;
//   Float_t yield  = func->Integral(0.45,1.05);
//   Float_t yieldE = func->IntegralError(0.45,1.05);
  Double_t tslope   = func->GetParameter(2);
  Double_t tslopeE = func->GetParError(2);	

  table.SetNextCol(partName);
  table.SetNextCol(yield,yieldE,-4);
  table.SetNextCol(tslope,tslopeE,-4);
  table.SetNextCol(func->GetParameter(1),func->GetParError(1),-4);
  table.SetNextCol(Form("%2.2f/%d",func->GetChisquare(),func->GetNDF()));
  cout << "5" << endl;
  Float_t lowestPoint = TMath::Max(AliBWTools::GetLowestNotEmptyBinEdge(h),min);
  cout << "6" << endl;
  //Float_t yieldAbove = 0;
  Float_t yieldAbove  = func->Integral(lowestPoint,100);
  table.SetNextCol(lowestPoint,-3);
  table.SetNextCol(yieldAbove/yield,-3);
  Float_t mean=0, meane=0;
  Float_t mean2=0, mean2e=0;
  cout << "6" << endl;
  AliBWTools::GetMean      (func, mean,  meane , 0.,100., normPar);
  AliBWTools::GetMeanSquare(func, mean2, mean2e, 0.,100., normPar);
  cout << "8" << endl;
  table.SetNextCol(mean,  meane ,-4);
  table.SetNextCol(mean2, mean2e,-4);
  //			 fMean2->IntegralError(0,100)/func->Integral(0,100),-7);
  table.InsertRow();
  table.PrintTable("ASCII");
  


}
