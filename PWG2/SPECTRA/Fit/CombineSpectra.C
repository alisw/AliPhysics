#include <iostream>
#include "TH1F.h"
#include "TGraphErrors.h"
#include "AliBWFunc.h"
#include "AliBWTools.h"
#include "TF1.h"
#include "TFile.h"
#include "TDatabasePDG.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TFolder.h"
#include "TStyle.h"
#include "AliLatexTable.h"
#include "TLegend.h"
#include "TVirtualFitter.h"
#include "TMath.h"
#include "TH2F.h"
#include "TSystem.h"
#include "TLine.h"
#include "TLatex.h"
#include "TMath.h"

#include "TASImage.h"
#include "TPaveText.h"
#include "ALICEWorkInProgress.C"
#include "StarPPSpectra.C"
#include "GetE735Ratios.C"

using namespace std;


// A bunch of useful enums and constants
enum {kPion=0,kKaon,kProton,kNPart};
enum {kTPC=0,kTOF,kITS,kITSTPC,kK0,kKinks,kCombTOFTPC,kNHist};// "k0" listed here as a kind of PID method...
const Int_t kNDet = kITS+2;
enum {kPos=0,kNeg,kNCharge};
enum {kPhojet=0,kPyTuneAtlasCSC, kPyTuneCMS6D6T, kPyTunePerugia0, kNTunes} ;
enum {kFitLevi=0, kFitUA1, kFitPowerLaw,
      kFitPhojet, kFitAtlasCSC, kFitCMS6D6T, kFitPerugia0,
      kNFit};


// flags, labels and names
const char * partFlag[] = {"Pion", "Kaon", "Proton"};
const char * detFlag[]  = {"TPC", "TOF", "ITS", "ITS Global", "K0", "Kinks", "Combined TOF + TPC"};
const char * chargeFlag[]  = {"Pos", "Neg"};
const char * chargeLabel[]  = {"Positive", "Negative"};
const char * partLabel[kNPart][kNCharge] = {{"#pi^{+}", "#pi^{-}"}, 
					    //					    {"K^{+} (#times 2)", "K^{-} (#times 2)"}, 
					    {"K^{+}", "K^{-}"}, 
					    {"p" ,  "#bar{p}"}};
const char * mcTuneName[] = {"Phojet", 
			     "Pythia - CSC 306", 
			     "Pythia - D6T 109", 
			     "Pythia - Perugia0 - 320", };
const char * funcName[] = { "Levi", "UA1", "PowerLaw", "Phojet", "AtlasCSC", "CMS6D6T", "Perugia0"};

// Style
//const Int_t marker[] = {25,24,28,20,21}; // set for highlithining marek
const Int_t marker[] = {20,24,25,28,21}; // standard set
const Int_t color [] = {1,2,4};

const Int_t mcLineColor[] = {kGreen,kRed,kBlue,kBlack};
const Int_t mcLineStyle[] = {1,2,3,4};


// Template needed to combine different detectors
const Float_t templBins[] = {0.05,0.1,0.15,0.20,0.25,0.30,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,2,2.2,2.4,2.6};
Int_t nbinsTempl=31;

TH1F * htemplate = new TH1F("htemplate", "htemplate",nbinsTempl, templBins );

//  Globals
TH1F * hSpectra[kNHist][kNPart][kNCharge];
TH1F * hSpectraMC[kNTunes][kNPart][kNCharge];
Double_t mass[kNPart];

// Functions:
// Loading
void LoadSpectra() ;
void LoadMC() ;

// Additional tasks (may be uncommented below)
void DrawStar(Int_t icharge);
void GetITSResiduals();
void DrawWithModels() ;
void DrawAllAndKaons();
void DrawWithJacek();
void DrawRatioToStar();
void DrawRatios();

// External stuff
//void ALICEWorkInProgress(TCanvas *c,TString today, TString label);

// Used to tag plots
TDatime dt;
TString today = "";



// Switches
Bool_t convertToMT = 0;
Bool_t doPrint = 0;
Int_t  fitFuncID = kFitLevi;
Bool_t scaleKaons =  kFALSE;
Bool_t correctSecondaries  = 1;
Bool_t correctGeantFlukaXS = 1;

void CombineSpectra() {

  // This macro is used to combine the 900 GeV spectra from 2009

  // Load stuff
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

  // Set date
  today = today +  long(dt.GetDay()) +"/" + long(dt.GetMonth()) +"/"+ long(dt.GetYear());


  // Set Masses
  mass[kPion]   = TDatabasePDG::Instance()->GetParticle("pi+")->Mass();
  mass[kKaon]   = TDatabasePDG::Instance()->GetParticle("K+")->Mass();
  mass[kProton] = TDatabasePDG::Instance()->GetParticle("proton")->Mass();

  // Load histos
  LoadSpectra();
  LoadMC();

  // Additional tasks
  //DrawStar(icharge);
  //  GetITSResiduals();
  //DrawAllAndKaons();  
  //  DrawWithModels() ;
  //DrawWithJacek();
  //DrawRatioToStar();
  //DrawRatios();
  //return;


  // Draw combined & Fit
  AliBWFunc * fm = new AliBWFunc;
  fm->SetVarType(AliBWFunc::kdNdpt);
  if (convertToMT) fm->SetVarType(AliBWFunc::kOneOverMtdNdmt);

  // table to print results
  AliLatexTable table(10,"c|cccccccc");
  if (fitFuncID == kFitLevi) {
    table.InsertCustomRow("Part & Yield & Yield (FIT) &  T Slope & n & $\\Chi^2$/NDF & Min X & Frac Above & \\langle p_{t} \\rangle  & \\langle p_{t}^{2} \\rangle");
  }  else if (fitFuncID == kFitPowerLaw) {
    table.InsertCustomRow("Part & Yield & Norm &  n & pt0 & $\\Chi^2$/NDF & Min X & Frac Above & \\langle p_{t} \\rangle  & \\langle p_{t}^{2} \\rangle");    
  } else {
    table.InsertCustomRow("Part & Yield & Par0 & Par2  & Par1 & $\\Chi^2$/NDF & Min X & Frac Above & \\langle p_{t} \\rangle  & \\langle p_{t}^{2} \\rangle");

  }

  AliLatexTable tempTable(4,"c|ccc");
  tempTable.InsertCustomRow("Part & Yield & Yield Below & Frac Above\\\\");


  //  Fit all  
  for(Int_t icharge = 0; icharge < kNCharge; icharge++){
    
    TCanvas * c2 = new TCanvas(TString("cCombined")+chargeFlag[icharge]+"_"+funcName[fitFuncID], TString("cCombined")+chargeFlag[icharge]);
    TH2F * hempty = new TH2F(TString("hempty")+long(icharge),"hempty",100,0.,2.9, 100, 0.001,5);
    hempty->SetXTitle("p_{t} (GeV/c)");
    hempty->SetYTitle("1/N_{ev} d^{2}N / dydp_{t} (GeV/c)^{-1}");
    hempty->Draw();
    //    DrawStar(icharge);
    c2->SetLogy();
    TLegend * l = new TLegend(0.516779, 0.729021, 0.89094 ,0.916084, chargeLabel[icharge]);
    l->SetFillColor(kWhite);

    for(Int_t ipart = 0; ipart < kNPart; ipart++){
      Float_t fitmin = 0;
      Float_t fitmax = 3;

      // Get functions
      TF1 * func = 0;
      if(fitFuncID == kFitLevi)          {
	if (ipart == kPion)
	  func = fm->GetLevi(mass[ipart], 0.12, 7, 1.5);
	if (ipart == kKaon)
	  func = fm->GetLevi(mass[ipart], 0.17, 7, 0.17);
	if (ipart == kProton)
	  func = fm->GetLevi(mass[ipart], 0.15, 8, 0.08);
      }
      else if(fitFuncID == kFitUA1)      func = fm->GetUA1(mass[ipart],0.2,1.25,1.25,0.2,1.5);
      else if(fitFuncID == kFitPowerLaw) {
	if (ipart == kPion)
	  func = fm->GetPowerLaw(1.0, 8.6, 7);
	if (ipart == kKaon)
	  func = fm->GetPowerLaw(3.0, 12, 2.6);
	if (ipart == kProton)
	  func = fm->GetPowerLaw(24, 72, 0.8);
      }
      else if(fitFuncID == kFitPhojet)   func = fm->GetHistoFunc(hSpectraMC[kPhojet]        [ipart][icharge]);
      else if(fitFuncID == kFitAtlasCSC) func = fm->GetHistoFunc(hSpectraMC[kPyTuneAtlasCSC][ipart][icharge]);
      else if(fitFuncID == kFitPerugia0) func = fm->GetHistoFunc(hSpectraMC[kPyTunePerugia0][ipart][icharge]);
      else if(fitFuncID == kFitCMS6D6T)  func = fm->GetHistoFunc(hSpectraMC[kPyTuneCMS6D6T] [ipart][icharge]);
      else {
	cout << "Unknown fit ID " << fitFuncID << endl;
	return;	
      }

      if(fitFuncID >= kFitPhojet){
	fitmin = 0.0;
	fitmax = 1.0;
      }

      if(!AliBWTools::Fit(hSpectra[kCombTOFTPC][ipart][icharge],func,fitmin,fitmax)) {
	cout << " FIT ERROR " << endl;
	return;      
      }
      hSpectra[kCombTOFTPC][ipart][icharge]->Draw("same");    
      hSpectra[kCombTOFTPC][ipart][icharge]->GetListOfFunctions()->At(0)->Draw("same");
      ((TF1*)hSpectra[kCombTOFTPC][ipart][icharge]->GetListOfFunctions()->At(0))->SetRange(0,4);
      l->AddEntry(hSpectra[kCombTOFTPC][ipart][icharge], 
		  scaleKaons && ipart == kKaon ? 
		  (TString(partLabel[ipart][icharge])+" #times 2").Data() 
		  : partLabel[ipart][icharge]);
//       TF1 * fClone = (TF1*) hSpectra[kCombTOFTPC][ipart][icharge]->GetListOfFunctions()->At(0)->Clone();
//       hSpectra[kCombTOFTPC][ipart][icharge]->GetListOfFunctions()->Add(fClone);
//       fClone->SetLineStyle(kDashed);
//       fClone->SetRange(0,100);
//       fClone->Draw("same");
      
      // populate table
      //     Float_t yield  = func->Integral(0.45,1.05);
      //     Float_t yieldE = func->IntegralError(0.45,1.05);
      
      Float_t yield  = func->Integral(0.,100);
      //Float_t yieldE = func->IntegralError(0.,100);

      Double_t yieldTools, yieldETools;
      Double_t partialYields[3],partialYieldsErrors[3]; 
      AliBWTools::GetYield(hSpectra[kCombTOFTPC][ipart][icharge], func, yieldTools, yieldETools, 
			   0, 100, partialYields,partialYieldsErrors);
      Double_t tslope   = func->GetParameter(2);
      Double_t tslopeE  = func->GetParError(2);	
      
      table.SetNextCol(partLabel[ipart][icharge]);
      //table.SetNextCol(yield,yieldE,-4);
      table.SetNextCol(yieldTools, yieldETools,-4);
      table.SetNextCol(func->GetParameter(0));
      table.SetNextCol(tslope,tslopeE,-4);
      table.SetNextCol(func->GetParameter(1)); 
      table.SetNextCol(Form("%2.2f/%d",func->GetChisquare(),func->GetNDF())); 
      Float_t lowestPoint = AliBWTools::GetLowestNotEmptyBinEdge(hSpectra[kCombTOFTPC][ipart][icharge]);
      //Float_t lowestPoint = AliBWTools::GetLowestNotEmptyBinEdge(hSpectra[kITS][ipart][icharge]);
      Float_t yieldAbove  = func->Integral(lowestPoint,100);
      table.SetNextCol(lowestPoint,-2);
      table.SetNextCol(yieldAbove/yield,-2);
      Float_t mean, meane;
      Float_t mean2, mean2e;
      AliBWTools::GetMean      (func, mean,  meane );
      AliBWTools::GetMeanSquare(func, mean2, mean2e);
      table.SetNextCol(mean,  meane ,-4);
      table.SetNextCol(mean2, mean2e,-4);
      
      //			 fMean2->IntegralError(0,100)/func->Integral(0,100),-7);
      table.InsertRow();
      

      /// TEMP TABLE
      tempTable.SetNextCol(partLabel[ipart][icharge]);
      tempTable.SetNextCol(yieldTools, yieldETools, -4);
      tempTable.SetNextCol(partialYields[1], partialYieldsErrors[1], -4);
      tempTable.SetNextCol(yieldAbove/yield,-2);
      tempTable.InsertRow();
    }
    l->Draw();
    if (doPrint) {
      c2->Update();
      gSystem->ProcessEvents();
      c2->Print(TString(c2->GetName()) + ".eps");
      ALICEWorkInProgress(c2,"","#splitline{ALICE Preliminary}{Statistical Error Only}");
      c2->Print(TString(c2->GetName()) + ".png");
    }
    
  }

  
  table.PrintTable("ASCII");
  
  cout << "" << endl;
  tempTable.PrintTable("ASCII");

}

void LoadSpectra() {

  TFile * f=0;

  // Load


  // TOF
  //  f = new TFile("./Files/spectra-pos-y.root");
  f = new TFile("./Files/spectra-pos-y_20100615.root");
  hSpectra[kTOF][kPion]  [kPos]= (TH1F*) f->Get("hpi");
  hSpectra[kTOF][kProton][kPos]= (TH1F*) f->Get("hpr");
  hSpectra[kTOF][kKaon]  [kPos]= (TH1F*) f->Get("hka");
  hSpectra[kTOF][kPion]  [kPos]->SetName("hpiPos");
  hSpectra[kTOF][kProton][kPos]->SetName("hprPos");
  hSpectra[kTOF][kKaon]  [kPos]->SetName("hkaPos");
  //f = new TFile("./Files/spectra-neg-y.root");
  f = new TFile("./Files/spectra-neg-y_20100615.root");
  hSpectra[kTOF][kPion]  [kNeg]= (TH1F*) f->Get("hpi");
  hSpectra[kTOF][kProton][kNeg]= (TH1F*) f->Get("hpr");
  hSpectra[kTOF][kKaon]  [kNeg]= (TH1F*) f->Get("hka");
  hSpectra[kTOF][kPion]  [kNeg]->SetName("hpiNeg");
  hSpectra[kTOF][kProton][kNeg]->SetName("hprNeg");
  hSpectra[kTOF][kKaon]  [kNeg]->SetName("hkaNeg");


  // Clean UP TPC spectra, removing unwanted points
  cout << "Cleaning Up TOF spectra" << endl;
  Int_t nbin =  hSpectra[kTOF][kKaon][kPos]->GetNbinsX();
  for(Int_t ibin = 1; ibin <= nbin; ibin++){
    Float_t pt =  hSpectra[kTOF][kKaon][kPos]->GetBinCenter(ibin);
    for(Int_t icharge = 0; icharge < kNCharge; icharge++){
      if (pt > 2.35) {
	hSpectra[kTOF][kKaon][icharge]->SetBinContent(ibin,0);
	hSpectra[kTOF][kKaon][icharge]->SetBinError  (ibin,0);	
	hSpectra[kTOF][kProton][icharge]->SetBinContent(ibin,0);
	hSpectra[kTOF][kProton][icharge]->SetBinError  (ibin,0);	
      }            
    }
  }
//   cout << "Scaling TOF to TPC" << endl;  
//   // Scale protons to TPC
//   hSpectra[kTOF][kProton][kPos]->Scale(1./1.05);
//   // Scale pbar so that pbar/p is compatible with Panos
//   hSpectra[kTOF][kProton][kNeg]->Scale(1./1.1);
  
  
  


  // ITS SA (Emanuele)
  //  f = new TFile("./Files/ITSsaSPECTRA_3clusters20100619.root");
  f = new TFile("./Files/ITSsaSPECTRAperMICHELE_20100703.root");
  hSpectra[kITS][kPion]  [kPos]= (TH1F*) f->Get("hSpectraPos0");
  hSpectra[kITS][kKaon]  [kPos]= (TH1F*) f->Get("hSpectraPos1");
  hSpectra[kITS][kProton][kPos]= (TH1F*) f->Get("hSpectraPos2");
  hSpectra[kITS][kPion]  [kNeg]= (TH1F*) f->Get("hSpectraNeg0");
  hSpectra[kITS][kKaon]  [kNeg]= (TH1F*) f->Get("hSpectraNeg1");
  hSpectra[kITS][kProton][kNeg]= (TH1F*) f->Get("hSpectraNeg2");

  for(Int_t ipart = 0; ipart < kNPart; ipart++){
    for(Int_t icharge = 0; icharge < kNCharge; icharge++){
      Int_t nbin = hSpectra[kITS][ipart][icharge]->GetNbinsX();
      //      hSpectra[kITS][ipart][icharge]->Scale(276004.);// Emanule divided his spectra...
      for(Int_t ibin = 1; ibin <= nbin; ibin++){
	if(hSpectra[kITS][ipart][icharge]->GetBinContent(ibin) < 0 ) {
	  hSpectra[kITS][ipart][icharge]->SetBinContent(ibin,0);
	  hSpectra[kITS][ipart][icharge]->SetBinError  (ibin,0);
	}
// 	if ((ipart == kKaon && ibin >= 12) || (ipart == kProton && ibin >= 20)) {
// 	  hSpectra[kITS][ipart][icharge]->SetBinContent(ibin,0);
// 	  hSpectra[kITS][ipart][icharge]->SetBinError  (ibin,0);
// 	}
      }
      
    }
  }


  // ITS + TPC (Marek)
  f = TFile::Open("./Files/SpectraCorrectedITSBeforeProtons_20100618.root");
  TList * list = (TList*) gDirectory->Get("output");
  hSpectra[kITSTPC][kPion]  [kPos]= (TH1F*) list->FindObject("Pions");
  hSpectra[kITSTPC][kKaon]  [kPos]= (TH1F*) list->FindObject("Kaons");
  hSpectra[kITSTPC][kProton][kPos]= (TH1F*) list->FindObject("Protons");
  hSpectra[kITSTPC][kPion]  [kNeg]= (TH1F*) list->FindObject("AntiPions");
  hSpectra[kITSTPC][kKaon]  [kNeg]= (TH1F*) list->FindObject("AntiKaons");
  hSpectra[kITSTPC][kProton][kNeg]= (TH1F*) list->FindObject("AntiProtons");

  // TPC
  //  htemplate =  hSpectra[kITS][kProton][kNeg]; //FIXME
  f = new TFile("./Files/protonSpectra_20100615.root");
  hSpectra[kTPC][kProton][kPos]= AliBWTools::GetHistoFromGraph((TGraphErrors*)f->Get("protonPosClassic"),htemplate);
  hSpectra[kTPC][kProton][kNeg]= AliBWTools::GetHistoFromGraph((TGraphErrors*)f->Get("protonNegClassic"),htemplate);
  f = new TFile("./Files/pionSpectra_20100615.root");
  hSpectra[kTPC][kPion][kPos]= AliBWTools::GetHistoFromGraph((TGraphErrors*)f->Get("pionPosClassic"),htemplate);
  hSpectra[kTPC][kPion][kNeg]= AliBWTools::GetHistoFromGraph((TGraphErrors*)f->Get("pionNegClassic"),htemplate);
  f = new TFile("./Files/kaonSpectra_20100615.root");
  hSpectra[kTPC][kKaon][kPos]= AliBWTools::GetHistoFromGraph((TGraphErrors*)f->Get("kaonPosClassic"),htemplate);
  hSpectra[kTPC][kKaon][kNeg]= AliBWTools::GetHistoFromGraph((TGraphErrors*)f->Get("kaonNegClassic"),htemplate);

  // Clean UP TPC spectra, removing unwanted points
  cout << "Cleaning Up TPC spectra" << endl;
  nbin =  hSpectra[kTPC][kKaon][kPos]->GetNbinsX();
  for(Int_t ibin = 0; ibin < nbin; ibin++){
    Float_t pt =  hSpectra[kTPC][kKaon][kPos]->GetBinCenter(ibin);
    if (pt > 0.45) {
      for(Int_t icharge = 0; icharge < kNCharge; icharge++){
	hSpectra[kTPC][kKaon][icharge]->SetBinContent(ibin,0);
	hSpectra[kTPC][kKaon][icharge]->SetBinError  (ibin,0);	
      }      
    }
    if (pt < 0.45) {
      for(Int_t icharge = 0; icharge < kNCharge; icharge++){
	hSpectra[kTPC][kProton][icharge]->SetBinContent(ibin,0);
	hSpectra[kTPC][kProton][icharge]->SetBinError  (ibin,0);	
      }      
    }
  }
  


 
  
  // K0s
  f = new TFile ("./Files/PtSpectraCorrectedK0sOff-dNdy-Jun02.root");
  //  hSpectra[kK0][kKaon][kPos] = (TH1F*) AliBWTools::GetdNdPtFromOneOverPt((TH1*) gDirectory->Get("hSpectraOff")); 
  hSpectra[kK0][kKaon][kPos] = (TH1F*) gDirectory->Get("hSpectraOff"); 
  //  hSpectra[kK0][kKaon][kPos]->Scale(2*TMath::Pi());
  //  hSpectra[kK0][kKaon][kPos]->Scale(1./272463);
  hSpectra[kK0][kKaon][kNeg] = hSpectra[kK0][kKaon][kPos];

  // Kinks: TO BE FIXED WITH POSITIVES AND NEGATIVES
  //  f = new TFile ("./Files/PtAllKaonKinkRap6Apr24.root");
  f = new TFile ("./Files/PtKaonKinkJune13AllPN_20100615.root");
  hSpectra[kKinks][kKaon][kPos] = (TH1F*)gDirectory->Get("fptallKPA");
  hSpectra[kKinks][kKaon][kNeg] = (TH1F*)gDirectory->Get("fptallKNA");


  // Apply correction factrs
  // Secondaries for protons
  f = new TFile ("./Files/corrFactorProtons_20100615.root");
  TH1F * hCorrSecondaries = (TH1F*) gDirectory->Get("corrFactorProtons");
  if(correctSecondaries) {
    cout << "CORRECTING SECONDARIES" << endl;
    
    for(Int_t idet = 0; idet <= kTOF; idet++){ // TPC and TOF only
      for(Int_t icharge = 0; icharge < kNCharge; icharge++){
	Int_t ipart = kProton;
	TH1* h = hSpectra[idet][ipart][icharge]; // lighter notation below
	if (h){
	  Int_t nbins = h->GetNbinsX();
	  for(Int_t ibin = 0; ibin < nbins; ibin++){
	    //	    Float_t pt = convertToMT ? TMath::Sqrt(h->GetBinCenter(ibin)*h->GetBinCenter(ibin)-mass[kProton]*mass[kProton]) : h->GetBinCenter(ibin);
	    Float_t pt = h->GetBinCenter(ibin);
	    if (icharge == kNeg) pt = -pt;
	    Int_t binCorrection = hCorrSecondaries->FindBin(pt);
	    Float_t correction    = hCorrSecondaries->GetBinContent(binCorrection);
	    //	    cout << pt << " " << correction << endl;
	    
	    if (correction != 0) {// If the bin is empty this is a  0
	      h->SetBinContent(ibin,h->GetBinContent(ibin)/correction);
	      h->SetBinError  (ibin,h->GetBinError  (ibin)/correction);
	    } else if (h->GetBinContent(ibin) > 0) { // If we are skipping a non-empty bin, we notify the user
	      cout << "Not correcting bin "<<ibin << " for protons secondaries, " << detFlag[idet] << " " << chargeFlag[icharge] << endl;
	      cout << " Bin content: " << h->GetBinContent(ibin)  << endl;
	      
	    }
	  }	
	}
      }
    }
  }
  // geant/fluka absorption
  if(correctGeantFlukaXS) {
    cout << "CORRECTING GEANT3/FLUKA" << endl;
    
    f = new TFile ("./Files/correctionForCrossSection_20100615.root");
    TH2D * hCorrFluka[kNCharge];
    hCorrFluka[kPos] = (TH2D*)gDirectory->Get("gHistCorrectionForCrossSectionProtons");
    hCorrFluka[kNeg] = (TH2D*)gDirectory->Get("gHistCorrectionForCrossSectionAntiProtons");
    for(Int_t idet = 0; idet < kNDet; idet++){
      if (idet == kITS) continue; // skip ITS sa
      for(Int_t icharge = 0; icharge < kNCharge; icharge++){
	Int_t ipart = kProton;
	TH1 * h = hSpectra[idet][ipart][icharge]; // lighter notation below
	if (h){
	  Int_t nbins = h->GetNbinsX();
	  for(Int_t ibin = 0; ibin < nbins; ibin++){
// 	    Float_t pt = convertToMT ? 
// 	      TMath::Sqrt(h->GetBinCenter(ibin)*h->GetBinCenter(ibin)-mass[kProton]*mass[kProton]) :
// 	      h->GetBinCenter(ibin);
	    Float_t pt = h->GetBinCenter(ibin);
	    Float_t minPtCorrection = hCorrFluka[icharge]->GetYaxis()->GetBinLowEdge(1);
	    Float_t maxPtCorrection = hCorrFluka[icharge]->GetYaxis()->GetBinLowEdge(h->GetNbinsY()+1);
	    if (pt < minPtCorrection) pt = minPtCorrection+0.0001;
	    if (pt > maxPtCorrection) pt = maxPtCorrection;
	    Float_t correction = hCorrFluka[icharge]->GetBinContent(1,hCorrFluka[icharge]->GetYaxis()->FindBin(pt));
	    if (idet == kTOF && icharge == kNeg) {

	      // Apply parametrized correction computed by francesco
	      // Fitted panos correction and using momentum at the outer radius of the TPC

	      Float_t ptav = pt; // Just to use the same name francesco uses...

	      // from pT constrained at P.V. (ptav) to pT TPC outer (ptTPCout)	      
	      Float_t ptTPCout=ptav*(1-6.81059e-01*TMath::Exp(-ptav*4.20094));
	      
	      // traking correction (fit to Panos)
	      Float_t antiprotonEC = 1 - 0.129758 *TMath::Exp(-ptav*0.679612);
	      
	      // TOF matching efficiency correction (derived from Panos one scaled for M.B.(TOF)/M.B.(TPC)).
	      Float_t antiprotonEC2 = TMath::Power(1 -  0.129758*TMath::Exp(-ptTPCout*0.679612),0.07162/0.03471);
	      correction = antiprotonEC * antiprotonEC2;
	    }
	    //	    cout << icharge<< " " << h->GetBinCenter(ibin) << " " << pt << " " << correction << endl;
	    if (correction != 0) {// If the bin is empty this is a  0
	      h->SetBinContent(ibin,h->GetBinContent(ibin)*correction);
	      h->SetBinError  (ibin,h->GetBinError  (ibin)*correction);
// 	      if (idet == kTOF) {
// 		cout << "CORRECTING TOF TWICE" << endl;
// 		h->SetBinContent(ibin,h->GetBinContent(ibin)*correction);
// 		h->SetBinError  (ibin,h->GetBinError  (ibin)*correction);		
// 	      }
	    } else if (h->GetBinContent(ibin) > 0) { // If we are skipping a non-empty bin, we notify the user
	      cout << "Fluka/GEANT: Not correcting bin "<<ibin << " for protons secondaries, " << detFlag[idet] << " " << chargeFlag[icharge] << endl;
	      cout << " Bin content: " << h->GetBinContent(ibin)  << endl;
	    }
	    
	  }
	  
	}
      }
    }    
  }


  // Set style and scale
  for(Int_t idet = 0; idet < kNDet; idet++){
    for(Int_t ipart = 0; ipart < kNPart; ipart++){
      for(Int_t icharge = 0; icharge < kNCharge; icharge++){
	if (hSpectra[idet][ipart][icharge]){
	  hSpectra[idet][ipart][icharge]->SetStats(0); // disable stats
	  hSpectra[idet][ipart][icharge]->SetMarkerColor (color[ipart] );
	  hSpectra[idet][ipart][icharge]->SetLineColor   (color[ipart] );
	  hSpectra[idet][ipart][icharge]->SetMarkerStyle (marker[ipart]);
	  hSpectra[idet][ipart][icharge]->SetXTitle("p_{t} (GeV/c)");
	  hSpectra[idet][ipart][icharge]->SetYTitle("1/N_{ev} dN/dp_{t} (GeV/c)^{-1}");
	  if (convertToMT) {
	    TH1F * htmp = (TH1F*) AliBWTools::GetOneOverPtdNdPt(hSpectra[idet][ipart][icharge]);
	    hSpectra[idet][ipart][icharge] = (TH1F*)AliBWTools::GetdNdmtFromdNdpt(htmp,mass[ipart]);
	    hSpectra[idet][ipart][icharge]->SetXTitle("m_{t} (GeV)");
	    hSpectra[idet][ipart][icharge]->SetYTitle("1/N_{ev} 1/m_{t} dN/dp_{t} (GeV)^{-1}");
	  }
// 	  if (idet == kTOF || idet == kTPC) {
// 	      hSpectra[idet][ipart][icharge]->Scale(1./236763);	      
// 	  } 
// 	  if (idet == kITS){
// 	    hSpectra[idet][ipart][icharge]->Scale(202945./236763);	      	      
// 	  }
	  if (scaleKaons && ipart == kKaon) {
	    hSpectra[idet][ipart][icharge]->Scale(2.);	      	    
	  }
	} else {
	  printf("Cannot load %s,%s,%s\n",detFlag[idet],partFlag[ipart],chargeFlag[icharge]);
	}
      }
    }
  }

  


  // Combine detectors
  for(Int_t ipart = 0; ipart < kNPart; ipart++){
    for(Int_t icharge = 0; icharge < kNCharge; icharge++){
      TH1F * htemplLocal = htemplate; // If we are converting to 1/mt dNdmt we need to convert the template as well...
      
      if (convertToMT) {
	TH1F * htmp = (TH1F*) AliBWTools::GetOneOverPtdNdPt(htemplate);
	htemplLocal = (TH1F*)AliBWTools::GetdNdmtFromdNdpt(htmp,mass[ipart]);

      }
      hSpectra[kCombTOFTPC][ipart][icharge] = AliBWTools::CombineHistos(hSpectra[kTOF][ipart][icharge],
									hSpectra[kTPC][ipart][icharge],
									htemplLocal,1.);;
//       if (convertToMT) {
// 	TH1F * htmp = (TH1F*) AliBWTools::GetOneOverPtdNdPt(hSpectra[kCombTOFTPC][ipart][icharge]);
// 	hSpectra[kCombTOFTPC][ipart][icharge] = (TH1F*)AliBWTools::GetdNdmtFromdNdpt(htmp,mass[ipart]);
// 	hSpectra[kCombTOFTPC][ipart][icharge]->SetXTitle("m_{t} (GeV)");
// 	hSpectra[kCombTOFTPC][ipart][icharge]->SetYTitle("1/N_{ev} 1/m_{t} dN/dp_{t} (GeV)^{-1}");
//       }
    }
  }


  // Scale for the number of inelastic collisions and correct for
  // efficiency losses due to physics selection:

  Double_t effPhysSel[kNPart];
  effPhysSel[kPion]   = 1.012;
  effPhysSel[kKaon]   = 1.013;
  effPhysSel[kProton] = 1.014;


  for(Int_t idet = 0; idet < kNHist; idet++){
    for(Int_t ipart = 0; ipart < kNPart; ipart++){
      for(Int_t icharge = 0; icharge < kNCharge; icharge++){
	if(hSpectra[idet][ipart][icharge]) {
	  //	  cout << "Scaling!" << endl;
	  hSpectra[idet][ipart][icharge]->Scale(1.*effPhysSel[ipart]/278366.15); // Scale PhysSel tutti? // FIXME
	}
      }
    }
  }


}

void LoadMC() {

  TFile * f = 0;
  const char * evClass= "INEL";
  const char * files[] = {"./Files/dndeta_Phojet_10M_900GeV.root", 
			  "./Files/dndeta_AtlasCSC306_10M_900GeV.root", 
			  "./Files/dndeta_CMSD6T109_10M_900GeV.root", 
			  "./Files/dndeta_Perugia0320_10M_900GeV.root", };
  
  // Phojet
  for(Int_t itune = 0; itune < kNTunes; itune++){
    f = new TFile(files[itune]);
    for(Int_t ipart = 0; ipart < kNPart; ipart++){
      for(Int_t icharge = 0; icharge < kNCharge; icharge++){
	hSpectraMC[itune][ipart][icharge] = (TH1F*) f->Get(Form("fHistPtID_%s_%s%s",evClass,partFlag[ipart],icharge==0 ? "Pos" : "Neg"));
      }
    }
  }
  

  // Set style
  for(Int_t itune = 0; itune < kNTunes; itune++){
    for(Int_t ipart = 0; ipart < kNPart; ipart++){
      for(Int_t icharge = 0; icharge < kNCharge; icharge++){
	if (hSpectraMC[itune][ipart][icharge]){
	  hSpectraMC[itune][ipart][icharge]->SetName(TString(hSpectraMC[itune][ipart][icharge]->GetName())+mcTuneName[itune]);
	  hSpectraMC[itune][ipart][icharge]->SetMarkerColor (mcLineColor[itune] );
	  hSpectraMC[itune][ipart][icharge]->SetLineColor   (mcLineColor[itune] );
	  hSpectraMC[itune][ipart][icharge]->SetLineStyle   (mcLineStyle[itune] );
	  hSpectraMC[itune][ipart][icharge]->SetMarkerStyle (1);
	  if (convertToMT) {
	    TH1F * htmp = (TH1F*)AliBWTools::GetOneOverPtdNdPt(hSpectraMC[itune][ipart][icharge]);
	    hSpectraMC[itune][ipart][icharge] = (TH1F*)AliBWTools::GetdNdmtFromdNdpt(htmp,mass[ipart]);
	    hSpectraMC[itune][ipart][icharge]->SetXTitle("m_{t} (GeV)");
	    hSpectraMC[itune][ipart][icharge]->SetYTitle("1/N_{ev} 1/m_{t} dN/dp_{t} (GeV)^{-1}");
	  }

	} else {
	  printf("Cannot load MC # %d,%s,%s\n",itune,partFlag[ipart],chargeFlag[icharge]);
	}
      }
    }
  }
  
}

void DrawStar(Int_t icharge) {
  //  cout << icharge << endl;
  
  //  gROOT->LoadMacro("StarPPSpectra.C");
  TGraphErrors ** gStar = StarPPSpectra();
  
  for(Int_t istar = 0; istar < 6; istar++){
    gStar[istar]->SetMarkerStyle(kOpenStar);
    if      (icharge==kPos &&  (istar%2) ) gStar[istar]->Draw("P");
    else if (icharge==kNeg && !(istar%2) ) gStar[istar]->Draw("P");
    else cout << "Skipping Star " << istar << endl;    
  }
}

void GetITSResiduals() {

 
  for(Int_t ipart = 0; ipart < kNPart; ipart++){
    for(Int_t icharge = 0; icharge < kNCharge; icharge++){
      //      cout << "1 " << ipart << " " << icharge << endl;
      
//       gSystem->ProcessEvents();
//       gSystem->Sleep(1000);
      TF1* f = (TF1*)   hSpectra[kCombTOFTPC][ipart][icharge]->GetListOfFunctions()->At(0);
      TH1F * hres1, *hres2;
      AliBWTools::GetResiduals(hSpectra[kITS][ipart][icharge], f, &hres1, &hres2);

      TCanvas * c1 = new TCanvas(TString(partFlag[ipart])+"_"+chargeFlag[icharge],TString(partFlag[ipart])+"_"+chargeFlag[icharge]);
      c1->SetLogy();
      hSpectra[kCombTOFTPC][ipart][icharge]->Draw();
      hSpectra[kITS][ipart][icharge]->SetMarkerStyle(24);
      hSpectra[kCombTOFTPC][ipart][icharge]->SetMarkerStyle(20);
      hSpectra[kITS][ipart][icharge]->Draw("same");
      hSpectra[kCombTOFTPC][ipart][icharge]->GetListOfFunctions()->At(0)->Draw("same");
      TLegend* l = new TLegend(    0.182886,    0.192308,    0.505034,0.384615, TString(partLabel[ipart][icharge])+" "+chargeFlag[icharge]);
      l->AddEntry(hSpectra[kCombTOFTPC][ipart][icharge], "TOF+TPC");
      l->AddEntry(hSpectra[kITS][ipart][icharge],        "ITS");
      l->AddEntry(f,        "Fit to TOF+TPC");
      l->Draw();

      TCanvas * c2 = new TCanvas(TString(partFlag[ipart])+"_"+chargeFlag[icharge]+"_res",
				 TString(partFlag[ipart])+"_"+chargeFlag[icharge]+"_res");  
      c2->SetGridy();
      hres2->SetMinimum(-1);
      hres2->SetMaximum(1);
      hres2->Draw();
      hres2->GetYaxis()->SetTitleOffset(1.2);
      Float_t x = AliBWTools::GetLowestNotEmptyBinEdge(hSpectra[kCombTOFTPC][ipart][icharge]);
      TLine * line = new TLine(x,-1,x,1);
      line->SetLineStyle(kDashed);
      line->Draw("same");
      
      if (doPrint) {
	c1->Update();
	c2->Update();
	gSystem->ProcessEvents();
	c1->Print(TString(c1->GetName()) + ".png");
	c2->Print(TString(c2->GetName()) + ".png");
      }
    }
  }
}

void DrawWithModels() {

  for(Int_t icharge = 0; icharge < kNCharge; icharge++){
 
    // Draw with models
    for(Int_t ipart = 0; ipart < kNPart; ipart++){
      // Multipad canvas
      TCanvas * c1 = new TCanvas(TString("cSpectra")+partFlag[ipart]+chargeFlag[icharge],TString("cSpectra")+partFlag[ipart]+chargeFlag[icharge] );
      TPad *p1 = new TPad(TString("p1")+partFlag[ipart]+chargeFlag[icharge], "p1", 0.0, 0.3, 1.0,  0.95, 0, 0, 0);
      p1->SetBottomMargin(0);
      p1->Draw();
      
      TPad *p2 = new TPad(TString("p2")+partFlag[ipart]+chargeFlag[icharge], "p2", 0.0, 0.05, 1.0,  0.3, 0, 0, 0);
      p2->SetTopMargin(0);
      p2->SetBottomMargin(0.4);
      p2->Draw();

      Float_t scaleFonts = (0.95-0.3)/(0.3-0.05);

      // Draw spectra
      p1->cd();
      p1->SetLogy();
      TH2F * hempty = new TH2F(TString("hempty")+long(ipart),"hempty",100,0.,4, 100, 0.0015,5);
      hempty->SetXTitle("p_{t} (GeV/c)");
      hempty->SetYTitle("1/N_{ev} d^{2}N / dydp_{t} (GeV/c)^{-1}");
      hempty->Draw();
      c1->SetLogy();


      TLegend * l =new TLegend(       0.543624,  0.431818,  0.892617,0.629371);
      l->SetFillColor(kWhite);
      hSpectra[kCombTOFTPC][ipart][icharge]->Draw("same");
      l->AddEntry(hSpectra[kTOF][ipart][icharge],TString ("Data: ")+partLabel[ipart][icharge]);
      for(Int_t itune = 0; itune < kNTunes; itune++){      
	l->AddEntry(hSpectraMC[itune][ipart][icharge],mcTuneName[itune]);
	hSpectraMC[itune][ipart][icharge]->SetLineWidth(2);    
	hSpectraMC[itune][ipart][icharge]->Draw("same,chist");    
      }
      l->Draw("same");

      // Draw ratio
      p2->cd();
      TH2F * hemptyr = new TH2F(TString("hemptyratio")+long(ipart),"hempty",100,0.,4, 100, 0.01,2.99);
      hemptyr->SetXTitle("p_{t} (GeV/c)");
      hemptyr->SetYTitle("Data/MC");
      hemptyr->GetXaxis()->SetLabelSize(0.04*scaleFonts);
      hemptyr->GetYaxis()->SetLabelSize(0.04*scaleFonts);
      hemptyr->GetYaxis()->SetTitleSize(0.05*scaleFonts);
      hemptyr->GetYaxis()->SetTitleOffset(1.4/scaleFonts);
      hemptyr->GetXaxis()->SetTitleSize(0.05*scaleFonts);
      hemptyr->SetTickLength(0.03*scaleFonts, "X");
      hemptyr->SetTickLength(0.02*scaleFonts, "Y");
      //      hemptyr->GetXaxis()->SetTitleOffset(1.4/scaleFonts);
      hemptyr->GetYaxis()->SetNdivisions(5);
      hemptyr->Draw("");

      AliBWFunc fm;
      for(Int_t itune = 0; itune < kNTunes; itune++){      
	TF1* f = fm.GetHistoFunc(hSpectraMC[itune][ipart][icharge], TString("f")+mcTuneName[itune]);

	//	l->AddEntry(hSpectraMC[itune][ipart][icharge],mcTuneName[itune]);
	TH1F* hRatio = AliBWTools::DivideHistoByFunc(hSpectra[kCombTOFTPC][ipart][icharge],f);
	hRatio->SetLineStyle(hSpectraMC[itune][ipart][icharge]->GetLineStyle());
	hRatio->SetLineColor(hSpectraMC[itune][ipart][icharge]->GetLineColor());
	hRatio->SetLineWidth(hSpectraMC[itune][ipart][icharge]->GetLineWidth());
	hRatio->Draw("lhist,same");
      }


      // print
      if(doPrint) {
	c1->Update();
	gSystem->ProcessEvents();
	c1->Print(TString(c1->GetName())+".eps");
	ALICEWorkInProgress(c1,"","#splitline{ALICE Preliminary}{Statistical Error Only}");
	c1->Print(TString(c1->GetName())+".png");
	
      }
    }
  }



}

void DrawAllAndKaons() {


  //  gROOT->LoadMacro("ALICEWorkInProgress.C");

  gStyle->SetOptFit(0);

  TH1F * hKaonsAllTPCTOF = (TH1F*) hSpectra[kCombTOFTPC][kKaon][kPos]->Clone();
  hKaonsAllTPCTOF->Add(hSpectra[kCombTOFTPC][kKaon][kNeg]);
  
  TH1F * hK0Scaled    = (TH1F*) hSpectra[kK0][kKaon][kPos]->Clone();
  hK0Scaled->Add(hSpectra[kK0][kKaon][kPos]);

  hSpectra[kKinks][kKaon][kPos]->SetMarkerStyle(25);
  hSpectra[kKinks][kKaon][kPos]->SetStats(0);
  TH1F * hKinksAll = (TH1F*) hSpectra[kKinks][kKaon][kPos]->Clone();
  hKinksAll->Add(hSpectra[kKinks][kKaon][kNeg]);
  
  TCanvas * c1 = new TCanvas("cKaons","cKaons");
  c1->SetLogy();
  TH2F * hempty = new TH2F("hempty_allkaons","hempty",100,0.,3, 100, 1e-4,10);
  hempty->SetXTitle("p_{t} (GeV/c)");
  hempty->SetYTitle("dN / dp_{t} (A.U.)");
  hempty->Draw();
  hKaonsAllTPCTOF->Draw("same");
  hK0Scaled->Draw("same");
  hKinksAll->Draw("same");

  TLegend * leg = new TLegend(0.2013423,0.2255245,0.5503356,0.4335664,NULL,"brNDC");
//    leg->SetBorderSize(0);
//    leg->SetLineColor(1);
//    leg->SetLineStyle(1);
//    leg->SetLineWidth(1);
//    leg->SetFillColor(19);
//    leg->SetFillStyle(1001);
   TLegendEntry *entry=leg->AddEntry("hkaPos_h_tpcKaonPos","K^{+} + K^{-}, TPC+TOF ","lpf");
   entry=leg->AddEntry("h1PtSpectraOff_inv","K^{0} #times 2","lpf");
   entry=leg->AddEntry("fptallK","K^{+} + K ^{-}, Kinks","lpf");
   leg->Draw();

   ALICEWorkInProgress(c1,today.Data(),"#splitline{ALICE Performance}{Not fully corrected}");
   TLatex * tex = new TLatex(0.2120805,0.01288336,"Statistical error only");
   tex->SetTextColor(2);
   tex->SetTextFont(42);
   tex->SetTextSize(0.03496503);
   tex->Draw();

   c1->Update();
   if(doPrint) c1->Print(TString(c1->GetName())+".png");

  // Draw all "stable" hadrons
  for(Int_t icharge = 0; icharge < kNCharge; icharge++){
    TCanvas * c1 = new TCanvas(TString("cAll_")+chargeFlag[icharge],TString("cAll_")+chargeFlag[icharge]);
    c1->SetLogy();
    TH2F * hempty = new TH2F(TString("hempty")+long(icharge),"hempty",100,0.,4, 100, 1e-4,10);
    hempty->SetXTitle("p_{t} (GeV/c)");
    hempty->SetYTitle("dN / dp_{t} (A.U.)");
    hempty->Draw();
    leg = new TLegend(  0.645973,  0.325175,  0.892617,0.636364, NULL,"brNDC");
    for(Int_t ipart = 0; ipart < kNPart; ipart++) {
      for(Int_t idet = 0; idet <= kITSTPC; idet++){
	//	if (idet == kITS) continue;
	// 	if (idet == kITSTPC) hSpectra[idet][ipart][icharge]->SetMarkerColor(kGreen);
 	hSpectra[idet][ipart][icharge]->SetMarkerStyle(marker[idet]);
	hSpectra[idet][ipart][icharge]->Draw("same");
	leg->AddEntry(hSpectra[idet][ipart][icharge],TString(partLabel[ipart][icharge])+" (" + detFlag[idet]  + ")","lpf");
      }
      //      leg->AddLine();
    }    
    leg->Draw();
    ALICEWorkInProgress(c1,today.Data(),"#splitline{ALICE Performance}{Not Fully Corrected}");
    c1->Update();
    if(doPrint) c1->Print(TString(c1->GetName())+".png");
  }
 

  //  Draw ratios (tmp)

//   new TCanvas;
//   hSpectra[kTPC][kKaon][kNeg]->Divide(hSpectra[kTPC][kKaon][kPos]);
//   hSpectra[kTPC][kKaon][kNeg]->Draw();
//   new TCanvas;
//   hSpectra[kITS][kKaon][kNeg]->Divide(hSpectra[kITS][kKaon][kPos]);
//   hSpectra[kITS][kKaon][kNeg]->Draw("");
  //    hSpectra[kTOF][kProton][kPos]->Draw("same");
  
//   for(Int_t icharge = 0; icharge < kNCharge; icharge++){
//     TCanvas * c1 = new TCanvas(TString("cAllRatio_")+chargeFlag[icharge],TString("cAllRatio_")+chargeFlag[icharge]);
//     c1->SetGridy();
//     TH2F * hempty = new TH2F(TString("hempty")+long(icharge),"hempty",100,0.2,1, 100, 0.0,2.0);
//     hempty->SetXTitle("p_{t} (GeV/c)");
//     hempty->SetYTitle("ITSsa / TPC");
//     hempty->Draw();
    

//     for(Int_t ipart = 0; ipart < kNPart; ipart++) {
//       hSpectra[kITS][ipart][icharge]->Divide(hSpectra[kTPC][ipart][icharge]);
//       hSpectra[kITS][ipart][icharge]->Draw("same");
//     }
//     if(doPrint) c1->Print(TString(c1->GetName())+".png");
//   }
    
// end of tmp

}

void DrawWithJacek() {

  //1. Convert spectra to dNdeta and sum
  TH1F * hsum = (TH1F*) htemplate->Clone();
  hsum->Reset();
  Int_t idet= kCombTOFTPC;
  for(Int_t icharge = 0; icharge < kNCharge; icharge++){
    for(Int_t ipart = 0; ipart < kNPart; ipart++){
      TH1 * h = hSpectra[idet][ipart][icharge];
      Int_t nbin = h->GetNbinsX();
      for(Int_t ibin = 1; ibin <= nbin; ibin++){
	Double_t pt = h->GetBinCenter(ibin);
	Double_t mt = TMath::Sqrt(pt*pt + mass[ipart]*mass[ipart]);
	Double_t jacobian = pt/mt;
	h->SetBinContent(ibin,h->GetBinContent(ibin)*jacobian);
	h->SetBinError  (ibin,h->GetBinError  (ibin)*jacobian);
	Int_t ibinSum = hsum->FindBin(pt);
	Double_t epsilon = 0.0001;
	if ( h->GetBinContent(ibin) > 0 && 
	     (TMath::Abs(h->GetBinLowEdge(ibin)   - hsum->GetBinLowEdge(ibinSum)) > epsilon || 
	      TMath::Abs(h->GetBinLowEdge(ibin+1) - hsum->GetBinLowEdge(ibinSum+1)) )
	     ) {
	  cout << "DISCREPANCY IN BIN RANGES" << endl;
	  cout << pt << " " << ibinSum << " " << ibin  << "; " << h->GetBinContent(ibin) << endl
	       << h->GetBinLowEdge(ibin) << "-"  << h->GetBinLowEdge(ibin+1) << endl
	       << hsum->GetBinLowEdge(ibinSum) << "-"  << hsum->GetBinLowEdge(ibinSum+1) << endl;
	  cout << "" << endl;	    
	}
	

	hsum->SetBinContent(ibinSum,hsum->GetBinContent(ibinSum)+h->GetBinContent(ibin)); // EROOR FIXME
	hsum->SetBinError  (ibinSum,0);
      }
      //	hsum->Add(h);
    }      
  }    


  // Load Jacek and Draw both:  
  new TFile ("./Files/dNdPt_Data_Points_ALICE_900GeV.root");
  TGraphErrors * gJacek = (TGraphErrors*) gDirectory->Get("inel");
  gJacek->Draw("AP");
  hsum->Draw("same");

  TGraphErrors * gRatio = AliBWTools::DivideGraphByHisto(gJacek,hsum);
  
 
  new TCanvas();
  gRatio->Draw("AP");

  

}


void DrawRatioToStar() {

  // Star data
  //  gROOT->LoadMacro("StarPPSpectra.C");
  TGraphErrors ** gStar = StarPPSpectra();

  // ALICE, INEL -> NSD
  Double_t scaleYield = 3.58/3.02; // from paper 2
  for(Int_t ipart = 0; ipart < kNPart; ipart++){
    for(Int_t icharge = 0; icharge < kNCharge; icharge++){
      hSpectra[kCombTOFTPC][ipart][icharge]->Scale(scaleYield);
    }    
  }
  
    
  TCanvas * c1 = new TCanvas("cRatioToStarNeg","cRatioToStarNeg");
  TH2F * hempty = new TH2F(TString("hemptyNeg"),"hemptyNeg",100,0.,1.5, 100, 0.001,1.8);
  hempty->SetXTitle("p_{t} (GeV/c)");
  hempty->SetYTitle("ALICE/STAR (NSD)");
  
  TLegend *leg = new TLegend(0.6510067,0.1853147,0.8892617,0.4178322,"Negative","brNDC");
  leg->SetBorderSize(0);
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(0);
  leg->SetFillStyle(1001);


  hempty->Draw();
  TGraphErrors * g ;
  g = AliBWTools::DivideGraphByHisto(gStar[0],hSpectra[kCombTOFTPC][kPion][kNeg],1);
  g->SetMarkerStyle(kFullCircle);
  g->SetMarkerColor(kBlack);
  g->Draw("p");
  leg->AddEntry(g,"#pi^{-}","lp");
  g = AliBWTools::DivideGraphByHisto(gStar[2],hSpectra[kCombTOFTPC][kKaon][kNeg],1);
  g->SetMarkerStyle(kOpenCircle);
  g->SetMarkerColor(kRed);
  g->Draw("p");
  leg->AddEntry(g,"K^{-}","lp");
  g = AliBWTools::DivideGraphByHisto(gStar[4],hSpectra[kCombTOFTPC][kProton][kNeg],1);
  g->SetMarkerStyle(kOpenSquare);
  g->SetMarkerColor(kBlue);
  g->Draw("p");
  leg->AddEntry(g,"#bar{p}","lp");  
  leg->Draw();




  TCanvas * c2 = new TCanvas("cRatioToStarPos","cRatioToStarPos");
  hempty->Draw(); 
  leg = new TLegend(0.6510067,0.1853147,0.8892617,0.4178322,"Positive","brNDC");
  leg->SetBorderSize(0);
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(0);
  leg->SetFillStyle(1001);
 //  TGraphErrors * g ;
  g = AliBWTools::DivideGraphByHisto(gStar[1],hSpectra[kCombTOFTPC][kPion][kPos],1);
  g->SetMarkerStyle(kFullCircle);
  g->SetMarkerColor(kBlack);
  g->Draw("p");
  leg->AddEntry(g,"#pi^{+}","lp");
  g = AliBWTools::DivideGraphByHisto(gStar[3],hSpectra[kCombTOFTPC][kKaon][kPos],1);
  g->SetMarkerStyle(kOpenCircle);
  g->SetMarkerColor(kRed);
  g->Draw("p");
  leg->AddEntry(g,"K^{+}","lp");
  g = AliBWTools::DivideGraphByHisto(gStar[5],hSpectra[kCombTOFTPC][kProton][kPos],1);
  g->SetMarkerStyle(kOpenSquare);
  g->SetMarkerColor(kBlue);
  g->Draw("p");
  leg->AddEntry(g,"p","lp");


  leg->Draw();

  c1->Update();
  c2->Update();
  gSystem->ProcessEvents();
  c1->Print(TString(c1->GetName()) + ".eps");
  c2->Print(TString(c2->GetName()) + ".eps");
  ALICEWorkInProgress(c1,today.Data(),"#splitline{ALICE Preliminary}{Statistical Error Only}");
  ALICEWorkInProgress(c2,today.Data(),"#splitline{ALICE Preliminary}{Statistical Error Only}");
  c1->Update();
  c2->Update();
  c1->Print(TString(c1->GetName()) + ".png");
  c2->Print(TString(c2->GetName()) + ".png");




}



void DrawRatios() {

  // Draws ratios of combined spectra

  // Compute ratios
  TH1F * hPosNegRatio[kNPart];
  
  for(Int_t ipart = 0; ipart < kNPart; ipart++){
    hPosNegRatio[ipart] = (TH1F*) hSpectra[kCombTOFTPC][ipart][kPos]->Clone();
    hPosNegRatio[ipart]->Divide(hSpectra[kCombTOFTPC][ipart][kNeg]);
    hPosNegRatio[ipart]->SetYTitle(TString(partLabel[ipart][kPos])+"/"+partLabel[ipart][kNeg]);    
    hPosNegRatio[ipart]->SetMinimum(0.5);
    hPosNegRatio[ipart]->SetMaximum(1.5);
  }
  
  TH1F * hKPiRatio = (TH1F*) hSpectra[kCombTOFTPC][kKaon][kPos]->Clone();
  hKPiRatio->Add(hSpectra[kCombTOFTPC][kKaon][kNeg]);
  TH1F * htmp = (TH1F*) hSpectra[kCombTOFTPC][kPion][kPos]->Clone();
  htmp->Add(hSpectra[kCombTOFTPC][kPion][kNeg]);
  hKPiRatio->Divide(htmp);
  hKPiRatio->SetYTitle("(K^{+}+K^{-})/(#pi^{+}+#pi^{-})");

  TH1F * hKPiRatioMC[kNTunes];
  for(Int_t itune = 0; itune < kNTunes; itune++){
    hKPiRatioMC[itune] = (TH1F*) hSpectraMC[itune][kKaon][kPos]->Clone();
    hKPiRatioMC[itune]->Add(hSpectraMC[itune][kKaon][kNeg]);
    TH1F * htmp = (TH1F*) hSpectraMC[itune][kPion][kPos]->Clone();
    htmp->Add(hSpectraMC[itune][kPion][kNeg]);
    hKPiRatioMC[itune]->Divide(htmp);
    hKPiRatioMC[itune]->SetYTitle("(K^{+}+K^{-})/(#pi^{+}+#pi^{-})");    
  }
  


  TH1F * hPPiRatio = (TH1F*) hSpectra[kCombTOFTPC][kProton][kPos]->Clone();
  hPPiRatio->Add(hSpectra[kCombTOFTPC][kProton][kNeg]);
  htmp = (TH1F*) hSpectra[kCombTOFTPC][kPion][kPos]->Clone();
  htmp->Add(hSpectra[kCombTOFTPC][kPion][kNeg]);
  hPPiRatio->Divide(htmp);
  hPPiRatio->SetYTitle("(p+#bar{p})/(#pi^{+}+#pi^{-})");

  TH1F * hPPiRatioMC[kNTunes];
  for(Int_t itune = 0; itune < kNTunes; itune++){
    hPPiRatioMC[itune] = (TH1F*) hSpectraMC[itune][kProton][kPos]->Clone();
    hPPiRatioMC[itune]->Add(hSpectraMC[itune][kProton][kNeg]);
    TH1F * htmp = (TH1F*) hSpectraMC[itune][kPion][kPos]->Clone();
    htmp->Add(hSpectraMC[itune][kPion][kNeg]);
    hPPiRatioMC[itune]->Divide(htmp);
    hPPiRatioMC[itune]->SetYTitle("(p+#bar{p})/(#pi^{+}+#pi^{-})");    
  }
 
  // Draw
//   TH2F * hempty = new TH2F(TString("hempty"),"hempty",100,0.,1.5, 100, 0.001,1.8);
//   hempty->SetXTitle("p_{t} (GeV/c)");
  //  hempty->SetYTitle("");

  // tmp: overlay levi fits
  AliBWFunc * fm2 = new AliBWFunc;
  fm2->SetVarType(AliBWFunc::kdNdpt);
  TF1 * fLevi[kNPart]  [kNCharge];
  fLevi[kPion]  [kPos] = fm2->GetLevi (mass[0], 0.1243, 7.614785, 1.524167, "fLeviPiPlus");
  fLevi[kKaon]  [kPos] = fm2->GetLevi (mass[1], 0.1625, 5.869318, 0.186361, "fLeviKPlus");
  fLevi[kProton][kPos] = fm2->GetLevi (mass[2], 0.1773, 6.918065, 0.086389, "fLeviPPlus");
  fLevi[kPion]  [kNeg] = fm2->GetLevi (mass[0], 0.1267, 7.979582, 1.515908, "fLeviPiNeg");
  fLevi[kKaon]  [kNeg] = fm2->GetLevi (mass[1], 0.1721, 6.927956, 0.191140, "fLeviKNeg");
  fLevi[kProton][kNeg] = fm2->GetLevi (mass[2], 0.1782, 8.160362, 0.082091, "fLeviPNeg");
  for(Int_t ipart = 0; ipart < kNPart; ipart++){
    for(Int_t icharge = 0; icharge < kNCharge; icharge++){
      fLevi[ipart][icharge]->SetRange(0,4);
    }    
  }
  

  for(Int_t ipart = 0; ipart < kNPart; ipart++){
    TCanvas * c1 = new TCanvas(TString("cRatio_")+partFlag[ipart], TString("cRatio_")+partFlag[ipart]);
    hPosNegRatio[ipart]->Draw();
    TF1 * fRatio = new TF1 (TString("fRatio")+partFlag[ipart], TString(fLevi[ipart][kPos]->GetName())+"/"+fLevi[ipart][kNeg]->GetName());
    //    fRatio->Draw("same");
    fRatio->SetRange(0,5);
    if (doPrint) {
      c1->Update();
      gSystem->ProcessEvents();
      c1->Print(TString(c1->GetName()) + ".eps");
    }
    
  }


  TCanvas * c2 = new TCanvas(TString("cRatio_KPi"),TString("cRatio_KPi"));  
  hKPiRatio->Draw();
  TLegend * lMC = new TLegend(0.526846, 0.18007, 0.887584,0.407343);
  lMC->SetFillColor(kWhite);

  //  gROOT->LoadMacro("GetE735Ratios.C");
  GetE735Ratios(0,0)->Draw("EX0,same");
  GetE735Ratios(0,1)->Draw("EX0,same");
  GetE735Ratios(0,2)->Draw("EX0,same");
  GetE735Ratios(0,3)->Draw("EX0,same");
  hKPiRatio->SetMarkerStyle(20);
  hKPiRatio->Draw("same");
  
  for(Int_t itune = 0; itune < kNTunes; itune++){
    lMC->AddEntry(hKPiRatioMC[itune],mcTuneName[itune]);
    hKPiRatioMC[itune]->SetLineWidth(2);    
    hKPiRatioMC[itune]->Draw("same,chist");    	
  }
  
  lMC->Draw();


  TLegend * l = new TLegend(  0.1879,  0.68,  0.54,0.92);
  l->SetFillColor(kWhite);
  l->AddEntry(hKPiRatio, "ALICE, #sqrt{s} = 900 GeV","lpf");
  l->AddEntry(GetE735Ratios(0,0), "E735, #sqrt{s} = 300 GeV","lpf");
  l->AddEntry(GetE735Ratios(0,1), "E735, #sqrt{s} = 540 GeV","lpf");
  l->AddEntry(GetE735Ratios(0,2), "E735, #sqrt{s} = 1000 GeV","lpf");
  l->AddEntry(GetE735Ratios(0,3), "E735, #sqrt{s} = 1800 GeV","lpf");
  l->Draw();


  TCanvas * c3 = new TCanvas(TString("cRatio_PPi"),TString("cRatio_PPi"));  
  hPPiRatio->Draw();
  hPPiRatio->SetMaximum(0.39);
  lMC = new TLegend(0.526846, 0.18007, 0.887584,0.407343);
  lMC->SetFillColor(kWhite);
  for(Int_t itune = 0; itune < kNTunes; itune++){
    lMC->AddEntry(hKPiRatioMC[itune],mcTuneName[itune]);
    hKPiRatioMC[itune]->SetLineWidth(2);    
    hKPiRatioMC[itune]->Draw("same,chist");    	
  }
  
  lMC->Draw();

  if (doPrint) {
    c2->Update();
    gSystem->ProcessEvents();
    c2->Print(TString(c2->GetName()) + ".eps");
    c3->Update();
    gSystem->ProcessEvents();
    c3->Print(TString(c3->GetName()) + ".eps");
  }


}

