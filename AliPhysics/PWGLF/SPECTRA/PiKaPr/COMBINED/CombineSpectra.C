#if !defined(__CINT__) || defined(__MAKECINT__)
#include <iostream>
#include "TH1F.h"
#include "TGraphErrors.h"
#include "AliPWGFunc.h"
#include "AliPWGHistoTools.h"
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
#include "StarPPSpectra.C"
#include "GetE735Ratios.C"
#include "TString.h"
#include "TObjString.h"
#include "TGraphAsymmErrors.h"

#endif


using namespace std;


// A bunch of useful enums and constants
enum {kPion=0,kKaon,kProton,kNPart};
enum {kTPC=0,kTOF,kITS,kITSTPC,kK0,kKinks,kCombTOFTPC,kCombAll,kNHist};// "k0" listed here as a kind of PID method...
const Int_t kNDet = kITS+2;
enum {kPos=0,kNeg,kNCharge};
enum {kPhojet=0,kPyTuneAtlasCSC, kPyTuneCMS6D6T, kPyTunePerugia0, kNTunes} ;
enum {kFitLevi=0, kFitUA1, kFitPowerLaw,
      kFitPhojet, kFitAtlasCSC, kFitCMS6D6T, kFitPerugia0,
      kNFit};
enum {kDoFits=0, kDoRatios, kDoSuperposition, kDoDrawWithModels, kDoCompareToStar, kDoDrawSyst, kDoCompareToAllCharged, kDoHelp};
enum {kStatErrors = 0, kSystErrors, kStatSystErrors}; // which errors do we put in the histo that we fit? stat,syst or stat+syst?

// flags, labels and names
const char * partFlag[] = {"Pion", "Kaon", "Proton"};
const char * detFlag[]  = {"TPC", "TOF", "ITS standalone", "ITS-TPC", "K0", "Kinks", "Combined TOF + TPC", "Combined TOF + TPC + ITS"};
const char * chargeFlag[]  = {"Pos", "Neg"};
const char * chargeLabel[]  = {"Positive", "Negative"};
const char * partLabel[kNPart][kNCharge] = {{"#pi^{+}", "#pi^{-}"}, 
					    //					    {"K^{+} (#times 2)", "K^{-} (#times 2)"}, 
					    {"K^{+}", "K^{-}"}, 
					    {"p" ,  "#bar{p}"}};
const char * partLatex[kNPart][kNCharge] = {{"$\\pi^{+}$", "$\\pi^{-}$"}, 
					    //					    {"K^{+} (#times 2)", "K^{-} (#times 2)"}, 
					    {"$K^{+}$", "$K^{-}$"}, 
					    {"$p$" ,  "$\\bar{p}$"}};
const char * mcTuneName[] = {"Phojet", 
			     "Pythia - CSC 306", 
			     "Pythia - D6T 109", 
			     "Pythia - Perugia0 - 320", };
const char * funcName[] = { "Levi", "UA1", "PowerLaw", "Phojet", "AtlasCSC", "CMS6D6T", "Perugia0"};

// Style
//const Int_t marker[] = {25,24,28,20,21}; // set for highlithining marek
//const Int_t marker[] = {20,24,25,28,21}; // standard set
const Int_t marker[] = {24,25,28,27,21}; // No full symbols
const Int_t color [] = {1,2,4};

const Int_t mcLineColor[] = {kGreen,kRed,kBlue,kBlack};
const Int_t mcLineStyle[] = {1,2,3,4};


// Template needed to combine different detectors
const Float_t templBins[] = {0.05,0.1,0.12,0.14,0.16,0.18,0.20,0.25,0.30,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,2,2.2,2.4,2.6};
Int_t nbinsTempl=34;

TH1F * htemplate = new TH1F("htemplate", "htemplate",nbinsTempl, templBins );

//  Globals
TH1F * hSpectra[kNHist][kNPart][kNCharge];
TH1F * hSpectraMC[kNTunes][kNPart][kNCharge];
TH1F * hSystError[kNHist][kNPart][kNCharge];
Double_t mass[kNPart];

// Functions:
// Loading
void LoadSpectra() ;
void LoadMC() ;

// Additional tasks (may be uncommented below)
void LoadLibs();
void DrawStar(Int_t icharge);
void GetITSResiduals();
void DrawWithModels() ;
void DrawAllAndKaons();
void DrawWithJacek();
void DrawRatioToStar();
void DrawRatios();
void DrawWithSyst();
void FitCombined();
void PrintCanvas(TCanvas* c,const TString formats) ;
void RebinHistosForRatios() ;
void Help();

// External stuff
void ALICEWorkInProgress(TCanvas *c,TString today="11/05/2010", TString label = "ALICE performance"){

  TPad *myPadLogo = new TPad("myPadLogo", "Pad for ALICE Logo",0.72,0.72,0.89,0.89);
  //TPad *myPadLogo = new TPad("myPadLogo", "Pad for ALICE Logo",0.72,0.62,0.85,0.75);
  myPadLogo->SetFillColor(0); 
  myPadLogo->SetBorderMode(0);
  myPadLogo->SetBorderSize(2);
  myPadLogo->SetFrameBorderMode(0);
  myPadLogo->SetLeftMargin(0.0);
  myPadLogo->SetTopMargin(0.0);
  myPadLogo->SetBottomMargin(0.0);
  myPadLogo->SetRightMargin(0.0);
  myPadLogo->SetFillStyle(0);
  myPadLogo->Draw();
  myPadLogo->cd();
  TASImage *myAliceLogo = new TASImage("alice_logo.png");
  myAliceLogo->Draw();
  c->cd();  
  TPaveText* t1=new TPaveText(0.418103, 0.837798, 0.656609, 0.888393,"NDC");
  t1->SetFillStyle(0);
  t1->SetBorderSize(0);
  t1->AddText(0.,0.,label);
  t1->SetTextColor(kRed);
  t1->SetTextFont(42);
  t1->SetTextSize(0.04);
  t1->Draw();
  TPaveText* t2=new TPaveText(0.418103, 0.80, 0.656609, 0.84,"NDC");
  t2->SetFillStyle(0);
  t2->SetBorderSize(0);
  t2->AddText(0.,0.,"pp at #sqrt{s} = 900 GeV (2009 data)");
  t2->SetTextColor(kRed);
  t2->SetTextFont(42);
  t2->SetTextSize(0.027);
  t2->Draw();
  TPaveText* t3=new TPaveText(0.418103, 0.76, 0.656609, 0.80,"NDC");
  t3->SetFillStyle(0);
  t3->SetBorderSize(0);
  t3->AddText(0.,0.,"Statistical and systematic errors");
  t3->SetTextColor(kRed);
  t3->SetTextFont(42);
  t3->SetTextSize(0.027);
  t3->Draw(); 
  // TPaveText* t2=new TPaveText(0.65,0.65,0.89,0.7,"NDC");
  // t2->SetFillStyle(0);
  // t2->SetBorderSize(0);
  // t2->SetTextColor(kRed);
  // t2->SetTextFont(52);
  // t2->AddText(0.,0.,today.Data());
  // t2->Draw();
}

// Used to tag plots
TDatime dt;
TString today = "";



// Switches
Bool_t convertToMT = 0;
Bool_t sumCharge = 0;
Int_t whatToFit = kStatErrors;
Bool_t scaleKaons =  kFALSE;
Bool_t drawStar =  kFALSE; // Overlay star when doing fits
Bool_t correctSecondaries  = 1;
Bool_t correctGeantFlukaXS = 1;
Int_t iCombInStudy = kCombAll; //kCombTOFTPC
Int_t fitFuncID = kFitLevi;
Bool_t showMC=kTRUE;
Bool_t showE735=kTRUE;
Bool_t useSecCorrFromDCA=kTRUE;
Bool_t flagPreliminary=kFALSE; // Add "preliminary" flag and logo to the plots
Bool_t doPrint = 1;
Bool_t applyEtaCutCorrection=kTRUE; // Fixes a bug in the ITSsa analysis
TString printFormats = "C,eps,root"; // format in which canvases will be printed, if PrintCanvas is called (not all prints are based on printcanvas at the moment). This is a comma separated list.


void CombineSpectra(Int_t analysisType=kDoHelp, Int_t  locfitFuncID = kFitLevi) {  //kDoSuperposition;//kDoDrawWithModels;// kDoFits; //kDoRatios;  

  // This macro is used to combine the 900 GeV spectra from 2009

  fitFuncID=locfitFuncID;

  // Load stuff
  LoadLibs();

  // Print Help and quit
  if (analysisType == kDoHelp) {
    Help();
    return;
  }


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
  if(analysisType==kDoSuperposition) DrawAllAndKaons();  
  else if(analysisType==kDoDrawWithModels)  DrawWithModels() ;
  else if(analysisType==kDoCompareToAllCharged) DrawWithJacek();
  else if(analysisType==kDoCompareToStar) DrawRatioToStar();
  else if(analysisType==kDoRatios) DrawRatios();
  else if(analysisType==kDoDrawSyst) DrawWithSyst();
  else if(analysisType==kDoFits) FitCombined();
  return;
}


void FitCombined() {
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

  // Draw combined & Fit
  AliPWGFunc * fm = new AliPWGFunc;
  fm->SetVarType(AliPWGFunc::kdNdpt);
  if (convertToMT) fm->SetVarType(AliPWGFunc::kOneOverMtdNdmt);

  // table to print results
  AliLatexTable table(10,"c|ccccccccc");
  if (fitFuncID == kFitLevi) {
    table.InsertCustomRow("Part & Yield & Yield (FIT) &  T Slope & n & $\\chi^2$/NDF & Min X & Extrap (%) & $\\langle p_{t} \\rangle (Fit)$  & $\\langle p_{t} \\rangle$ (FIT+DATA)\\\\");
  }  else if (fitFuncID == kFitPowerLaw) {
    table.InsertCustomRow("Part & Yield & Norm &  n & pt0 & $\\chi^2$/NDF & Min X & Frac Above  & $\\langle p_{t} \\rangle (Fit)$  & $\\langle p_{t} \\rangle$ (FIT+DATA)  \\\\");    
  } else {
    table.InsertCustomRow("Part & Yield & Par0 & Par2  & Par1 & $\\chi^2$/NDF & Min X & Frac Above  & $\\langle p_{t} \\rangle (Fit)$  & $\\langle p_{t} \\rangle$ (FIT+DATA)  \\\\");

  }
  table.InsertHline();
  AliLatexTable tempTable(3,"c|cc");
  tempTable.InsertCustomRow("Part & Yield &  $\\langle p_{t} \\rangle$ \\\\");
  tempTable.InsertHline();

  TH1F* hRatiosToFit[kNPart][kNCharge]; // Ratio data/fit function, stat error by default
  TH1F* hRatiosToFitSyst[kNPart][kNCharge]; // Ratio data/fit, stat + syst
  //  Fit all 
  Int_t chargeLoop = sumCharge ? 1 : 2; 
  Bool_t divideRatiosCanvas=kTRUE;
  for(Int_t icharge = 0; icharge < chargeLoop; icharge++){

    TCanvas * c2 = new TCanvas(TString("cCombined")+chargeFlag[icharge]+"_"+funcName[fitFuncID], TString("cCombined")+chargeFlag[icharge],700,700);
    c2->SetTickx();
    c2->SetTicky();
    c2->SetLeftMargin(0.14);
    TCanvas * c2r = 0;
    if(divideRatiosCanvas) c2r = new TCanvas(TString("cCombinedRatio")+chargeFlag[icharge]+"_"+funcName[fitFuncID], TString("cCombinedRatio")+chargeFlag[icharge],1200,500);
    else                   c2r = new TCanvas(TString("cCombinedRatio")+chargeFlag[icharge]+"_"+funcName[fitFuncID], TString("cCombinedRatio")+chargeFlag[icharge],700,700);
    c2->cd();
    gPad->SetLogy();
    TH2F * hempty = new TH2F(TString("hempty")+long(icharge),"hempty",100,0.,2.9, 100, 0.0005,5);
    hempty->SetXTitle("p_{t} (GeV/c)");
    hempty->SetYTitle("1/N_{ev} d^{2}N / dydp_{t} (GeV/c)^{-1}");
    hempty->GetYaxis()->SetTitleOffset(1.35);
    hempty->GetXaxis()->SetTitleOffset(1.1);
    hempty->Draw();
    c2r->cd();
    if(divideRatiosCanvas) {
      c2r->Divide(3,1);
      c2r->cd(1);
    }
    gPad->SetGridy();
    TH2F * hemptyR = new TH2F(TString("hemptyR")+long(icharge),"hemptyR",100,0.,2.9, 100, 0.5,1.5);
    hemptyR->SetXTitle("p_{t} (GeV/c)");
    hemptyR->SetYTitle("Data/Fit");
    hemptyR->Draw();
    

    TLegend * l = new TLegend(0.176724, 0.181548, 0.477011, 0.462798,chargeLabel[icharge]);
    l->SetFillColor(kWhite);
    l->SetTextSize(0.035);
    

    
    TPaveText* tf=new TPaveText(0.321839, 0.175595, 0.685345, 0.299107,"NDC");
    if(fitFuncID == kFitLevi){
      tf->AddText("#frac{dN}{dp_{t}} #propto p_{t} #left(1+#frac{#sqrt{m^{2}+p_{t}^{2}} -m}{nT} #right)^{-n}");
      //      tf->SetNDC();
      tf->SetTextFont(12);
      tf->SetTextSize(0.032);
    }

    for(Int_t ipart = 0; ipart < kNPart; ipart++){
      printf(" ----- Fit %s %s ------\n",partFlag[ipart],chargeFlag[icharge]);
      Float_t fitmin = 0;
      Float_t fitmax = 3;

      // Get functions
      TF1 * func = 0;
      Int_t normPar = -1;
      if(fitFuncID == kFitLevi)          {
	normPar = 0; // The levi is normalized by parameter 0
	if (ipart == kPion)
	  func = fm->GetLevi(mass[ipart], 0.12, 7, 1.5);
	if (ipart == kKaon)
	  func = fm->GetLevi(mass[ipart], 0.17, 7, 0.17);
	if (ipart == kProton)
	  func = fm->GetLevi(mass[ipart], 0.15, 8.5, 0.09);
      }      
      else if(fitFuncID == kFitUA1)      func = fm->GetUA1(mass[ipart],0.2,1.25,1.25,0.2,1.5);
      else if(fitFuncID == kFitPowerLaw) {
	// fitmin=0;
	// fitmax=1;
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

      // Temp: fit histo with sist errors 
      TH1F * hsyst     = new TH1F(*htemplate);
      TH1F * hsyststat = 0;
      hsyst->SetFillColor(kGray);
      AliPWGHistoTools::GetValueAndError(hsyst,hSpectra[iCombInStudy][ipart][icharge],hSystError[iCombInStudy][ipart][icharge],kTRUE);
      hsyststat= new TH1F(*hsyst);
      AliPWGHistoTools::GetHistoCombinedErrors(hsyststat,hSpectra[iCombInStudy][ipart][icharge]); // combine syst and stat
      hsyststat->SetFillColor(kGray);

      TH1F * hToFit = 0;
      if (whatToFit == kStatErrors) hToFit = hSpectra[iCombInStudy][ipart][icharge]; // Shorthand
      if (whatToFit == kStatSystErrors) {
	hToFit = hsyststat;// Shorthand
      }
      if (whatToFit == kSystErrors) hToFit = hsyst;
      if (sumCharge) {
	if(whatToFit == kStatErrors)
	  hToFit->Add(hSpectra[iCombInStudy][ipart][1]);
	else if (whatToFit == kStatSystErrors) {
	  TH1F * hsyst2     = new TH1F(*htemplate);
	  TH1F * hsyststat2 = 0;
	  hsyst2->SetFillColor(kGray);
	  AliPWGHistoTools::GetValueAndError(hsyst2,hSpectra[iCombInStudy][ipart][1],hSystError[iCombInStudy][ipart][1],kTRUE);
	  hsyststat2= new TH1F(*hsyst2);
	  AliPWGHistoTools::GetHistoCombinedErrors(hsyststat2,hSpectra[iCombInStudy][ipart][1]); // combine syst and stat
	  hToFit->Add(hsyststat2);
	}
      }

      

      if(!AliPWGHistoTools::Fit(hToFit,func,fitmin,fitmax)) {
	cout << " FIT ERROR " << endl;
	return;      
      }
      cout << "DRAWING" << endl;
      c2->cd();      
      hsyststat->Draw("same,e2");    
      // TH1F * hsyststat2 = new TH1F(*hsyststat);
      //      hsyststat2->SetFillColor(kWhite);
      //      hsyststat2->SetLineColor(kBlack);
      //      hsyststat2->Draw("same,e2");    
      hToFit->Draw("same");    
      TF1* fitfunc=(TF1*)hToFit->GetListOfFunctions()->At(0);
      fitfunc->Draw("same");
      fitfunc->SetRange(0,4);
      fitfunc->SetLineColor(hSpectra[iCombInStudy][ipart][icharge]->GetLineColor());
      if(drawStar)    DrawStar(icharge);
      hRatiosToFit[ipart][icharge]=(TH1F*)hToFit->Clone(Form("hRatio%s%s",chargeFlag[icharge],partFlag[icharge])); // Ratio data/fit
      hRatiosToFitSyst[ipart][icharge]=(TH1F*)hsyststat->Clone(Form("hRatioSyst%s%s",chargeFlag[icharge],partFlag[icharge])); // Ratio data/fit
      // Syst error on ratio with no fill
      hRatiosToFitSyst[ipart][icharge]->SetLineColor(hToFit->GetLineColor());
      hRatiosToFitSyst[ipart][icharge]->SetFillStyle(0);
      

      // Compute ratio data/function integrating the function in the bin width
      for(Int_t iBin=1; iBin<hToFit->GetNbinsX(); iBin++){
	Double_t lowLim=hToFit->GetBinLowEdge(iBin);
	Double_t highLim=hToFit->GetBinLowEdge(iBin+1);
	Double_t contFunc=fitfunc->Integral(lowLim,highLim)/(highLim-lowLim);
	Double_t ratio=hToFit->GetBinContent(iBin)/contFunc;
	Double_t eratio=hToFit->GetBinError(iBin)/contFunc;
	Double_t eratioSyst=hsyststat->GetBinError(iBin)/contFunc;
	hRatiosToFitSyst[ipart][icharge]->SetBinContent(iBin,ratio);
	hRatiosToFitSyst[ipart][icharge]->SetBinError(iBin,eratioSyst);
	hRatiosToFit[ipart][icharge]->SetBinContent(iBin,ratio);
	hRatiosToFit[ipart][icharge]->SetBinError(iBin,eratio);
      }
      c2->Update();
      l->AddEntry(hToFit, 
		  scaleKaons && ipart == kKaon ? 
		  (TString(partLabel[ipart][icharge])+" #times 2").Data() 
		  : partLabel[ipart][icharge]);
      
      // populate table
      //     Float_t yield  = func->Integral(0.45,1.05);
      //     Float_t yieldE = func->IntegralError(0.45,1.05);
      
      Float_t yield  = func->Integral(0.,100);
      //Float_t yieldE = func->IntegralError(0.,100);

      Double_t yieldTools, yieldETools;
      Double_t partialYields[3],partialYieldsErrors[3]; 
      AliPWGHistoTools::GetYield(hToFit, func, yieldTools, yieldETools, 
			   0, 100, partialYields,partialYieldsErrors);
      Double_t tslope   = func->GetParameter(2);
      Double_t tslopeE  = func->GetParError(2);	
      
      table.SetNextCol(partLatex[ipart][icharge]);
      //table.SetNextCol(yield,yieldE,-4);
      table.SetNextCol(yieldTools, yieldETools,-4);
      table.SetNextCol(func->GetParameter(0));
      table.SetNextCol(tslope,tslopeE,-4);
      table.SetNextCol(func->GetParameter(1),func->GetParError(1),-2); 
      table.SetNextCol(Form("%2.2f/%d",func->GetChisquare(),func->GetNDF())); 
      Float_t lowestPoint = AliPWGHistoTools::GetLowestNotEmptyBinEdge(hToFit);
      //Float_t lowestPoint = AliPWGHistoTools::GetLowestNotEmptyBinEdge(hSpectra[kITS][ipart][icharge]);
      Float_t yieldAbove  = func->Integral(lowestPoint,100);
      table.SetNextCol(lowestPoint,-2);
      table.SetNextCol((1 - yieldAbove/yield)*100,-1);
      Float_t mean, meane;  // Mean using only fit function
      Double_t meanDF, meanDFe; // Mean from data + fit 
      Float_t mean2, mean2e;
      if (fitFuncID < kFitPhojet) {
	AliPWGHistoTools::GetMean      (func, mean,  meane , 0.,100., normPar);
	AliPWGHistoTools::GetMeanSquare(func, mean2, mean2e, 0.,100., normPar);
      }
      //      AliPWGHistoTools::GetMeanDataAndExtrapolation(hToFit, func, meanDF, meanDFe, 0.,100.);
      AliPWGHistoTools::GetMeanDataAndExtrapolation(hToFit, func, meanDF, meanDFe, 0.,100.);
      table.SetNextCol(mean,   meane  ,-4);
      table.SetNextCol(meanDF, meanDFe,-4);
      
      //			 fMean2->IntegralError(0,100)/func->Integral(0,100),-7);
      table.InsertRow();
      

      /// TEMP TABLE
      tempTable.SetNextCol(partLatex[ipart][icharge]);
      tempTable.SetNextCol(yieldTools, yieldETools, -4);
      tempTable.SetNextCol(meanDF,  meanDFe ,-4);
      //      tempTable.SetNextCol(partialYields[1], partialYieldsErrors[1], -4);
      //      tempTable.SetNextCol(yieldAbove/yield,-2);
      tempTable.InsertRow();
      c2r->cd();
      if(divideRatiosCanvas) {
	c2r->cd(ipart+1);
	gPad->SetGridy();
	TH2F * heloc = new TH2F(*hemptyR);
	heloc->Draw();
      }
      hRatiosToFitSyst[ipart][icharge]->Draw("e2same");
      hRatiosToFit[ipart][icharge]->Draw("esame");

    }
    c2->cd();
    if (flagPreliminary) ALICEWorkInProgress(c2,"","ALICE Preliminary");
    l->Draw();
    if (flagPreliminary) tf->Draw();
    c2r->cd();
    if(flagPreliminary) ALICEWorkInProgress(c2r,"","ALICE Preliminary");
    l->Draw();
    PrintCanvas(c2,printFormats);
    PrintCanvas(c2r,printFormats);
    
    
  }

  
  table.PrintTable("ASCII");
  cout << "" << endl;
  table.PrintTable("");
  
  cout << "" << endl;
  tempTable.PrintTable("ASCII");



}

void LoadSpectra() {

  // Loads spectra for different detectors and corresponding systematic errors.

  TFile * f=0;

  // Systematic errors: initialize histos
  gROOT->cd();
  for(Int_t ipart = 0; ipart < kNPart; ipart++){
    for(Int_t idet = 0; idet < kNHist; idet++){
      for(Int_t icharge = 0; icharge < kNCharge; icharge++){
	hSystError[idet][ipart][icharge] = new TH1F (TString("hSyst_")+detFlag[idet]+partFlag[ipart]+chargeFlag[icharge],
						     TString("hSyst_")+detFlag[idet]+partFlag[ipart]+chargeFlag[icharge],
						     nbinsTempl,templBins);
	hSystError[idet][ipart][icharge]->SetMarkerColor (color[ipart] );
	hSystError[idet][ipart][icharge]->SetLineColor   (color[ipart] );

      }
      
    }
    
  }

  // Load


  // TOF
  // Load Efficiencies
  f =  new TFile("./Files/effTOFhistosSmooth.root");
  TFile * ftrk =  new TFile("./Files/effhistos.root");
  TH1D * hEffTrackTOF[kNPart][kNCharge];
  TH1D * hEffMatchTOF[kNPart][kNCharge];
  hEffTrackTOF[kPion]  [kPos] = (TH1D*) ftrk->Get("hpitrk_pos");
  hEffTrackTOF[kKaon]  [kPos] = (TH1D*) ftrk->Get("hkatrk_pos");
  hEffTrackTOF[kProton][kPos] = (TH1D*) ftrk->Get("hprtrk_pos");
  hEffMatchTOF[kPion]  [kPos] = (TH1D*) f->Get("hpieff_pos");
  hEffMatchTOF[kKaon]  [kPos] = (TH1D*) f->Get("hkaeff_pos");
  hEffMatchTOF[kProton][kPos] = (TH1D*) f->Get("hpreff_pos");
  hEffTrackTOF[kPion]  [kNeg] = (TH1D*) ftrk->Get("hpitrk_neg");
  hEffTrackTOF[kKaon]  [kNeg] = (TH1D*) ftrk->Get("hkatrk_neg");
  hEffTrackTOF[kProton][kNeg] = (TH1D*) ftrk->Get("hprtrk_neg");
  hEffMatchTOF[kPion]  [kNeg] = (TH1D*) f->Get("hpieff_neg");
  hEffMatchTOF[kKaon]  [kNeg] = (TH1D*) f->Get("hkaeff_neg");
  hEffMatchTOF[kProton][kNeg] = (TH1D*) f->Get("hpreff_neg");

  //  f = new TFile("./Files/spectra-pos-y.root");
  f = new TFile("./Files/spectraRaw-pos-y.root");
  hSpectra[kTOF][kPion]  [kPos]= (TH1F*) f->Get("hpi");
  hSpectra[kTOF][kProton][kPos]= (TH1F*) f->Get("hpr");
  hSpectra[kTOF][kKaon]  [kPos]= (TH1F*) f->Get("hka");
  hSpectra[kTOF][kPion]  [kPos]->SetName("hpiPos");
  hSpectra[kTOF][kProton][kPos]->SetName("hprPos");
  hSpectra[kTOF][kKaon]  [kPos]->SetName("hkaPos");
  //f = new TFile("./Files/spectra-neg-y.root");
  f = new TFile("./Files/spectraRaw-neg-y.root");
  hSpectra[kTOF][kPion]  [kNeg]= (TH1F*) f->Get("hpi");
  hSpectra[kTOF][kProton][kNeg]= (TH1F*) f->Get("hpr");
  hSpectra[kTOF][kKaon]  [kNeg]= (TH1F*) f->Get("hka");
  hSpectra[kTOF][kPion]  [kNeg]->SetName("hpiNeg");
  hSpectra[kTOF][kProton][kNeg]->SetName("hprNeg");
  hSpectra[kTOF][kKaon]  [kNeg]->SetName("hkaNeg");

  // Divide for efficiency
  hSpectra[kTOF][kPion]  [kPos]->Divide(hEffTrackTOF[kPion]  [kPos]);
  hSpectra[kTOF][kKaon]  [kPos]->Divide(hEffTrackTOF[kKaon]  [kPos]);
  hSpectra[kTOF][kProton][kPos]->Divide(hEffTrackTOF[kProton][kPos]);
  hSpectra[kTOF][kPion]  [kPos]->Divide(hEffMatchTOF[kPion]  [kPos]);
  hSpectra[kTOF][kKaon]  [kPos]->Divide(hEffMatchTOF[kKaon]  [kPos]);
  hSpectra[kTOF][kProton][kPos]->Divide(hEffMatchTOF[kProton][kPos]);
  hSpectra[kTOF][kPion]  [kNeg]->Divide(hEffTrackTOF[kPion]  [kNeg]);
  hSpectra[kTOF][kKaon]  [kNeg]->Divide(hEffTrackTOF[kKaon]  [kNeg]);
  hSpectra[kTOF][kProton][kNeg]->Divide(hEffTrackTOF[kProton][kNeg]);
  hSpectra[kTOF][kPion]  [kNeg]->Divide(hEffMatchTOF[kPion]  [kNeg]);
  hSpectra[kTOF][kKaon]  [kNeg]->Divide(hEffMatchTOF[kKaon]  [kNeg]);
  hSpectra[kTOF][kProton][kNeg]->Divide(hEffMatchTOF[kProton][kNeg]);


  // Clean UP TOF spectra, removing unwanted points
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
      if (pt < 0.7) {
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
  
  // TOF: systematics
  // Load TOF systematic errors:
  f = new TFile ("./Files/systMatchingPos.root");
  AliPWGHistoTools::AddHisto(hSystError[kTOF][kPion][kPos]  ,(TH1*)gDirectory->Get("hErrPionMatch")  );
  AliPWGHistoTools::AddHisto(hSystError[kTOF][kProton][kPos],(TH1*)gDirectory->Get("hErrProtonMatch"));
  AliPWGHistoTools::AddHisto(hSystError[kTOF][kKaon][kPos]  ,(TH1*)gDirectory->Get("hErrKaonMatch")  );
  f = new TFile ("./Files/systMatchingNeg.root");
  AliPWGHistoTools::AddHisto(hSystError[kTOF][kPion]  [kNeg],(TH1*)gDirectory->Get("hErrPionMatch"));
  AliPWGHistoTools::AddHisto(hSystError[kTOF][kProton][kNeg],(TH1*)gDirectory->Get("hErrProtonMatch"));
  AliPWGHistoTools::AddHisto(hSystError[kTOF][kKaon]  [kNeg],(TH1*)gDirectory->Get("hErrKaonMatch"));  
  f = new TFile ("./Files/systPIDall.root");
  AliPWGHistoTools::AddHisto(hSystError[kTOF][kPion]  [kPos],(TH1*)gDirectory->Get("hpiCorr"));
  AliPWGHistoTools::AddHisto(hSystError[kTOF][kProton][kPos],(TH1*)gDirectory->Get("hpCorr"));
  AliPWGHistoTools::AddHisto(hSystError[kTOF][kKaon]  [kPos],(TH1*)gDirectory->Get("hkCorr"));  
  AliPWGHistoTools::AddHisto(hSystError[kTOF][kPion]  [kNeg],(TH1*)gDirectory->Get("hpiCorr"));
  AliPWGHistoTools::AddHisto(hSystError[kTOF][kProton][kNeg],(TH1*)gDirectory->Get("hpCorr"));
  AliPWGHistoTools::AddHisto(hSystError[kTOF][kKaon]  [kNeg],(TH1*)gDirectory->Get("hkCorr"));  
  
  
  // ITS SA 
  f = new TFile("./Files/ITSsaSpectraCorr-28Sep10.root");
  hSpectra[kITS][kPion]  [kPos]= (TH1F*) f->Get("hSpectraPos0");
  hSpectra[kITS][kKaon]  [kPos]= (TH1F*) f->Get("hSpectraPos1");
  hSpectra[kITS][kProton][kPos]= (TH1F*) f->Get("hSpectraPos2");
  hSpectra[kITS][kPion]  [kNeg]= (TH1F*) f->Get("hSpectraNeg0");
  hSpectra[kITS][kKaon]  [kNeg]= (TH1F*) f->Get("hSpectraNeg1");
  hSpectra[kITS][kProton][kNeg]= (TH1F*) f->Get("hSpectraNeg2");

  // Remove unwanted bins
  for(Int_t ipart = 0; ipart < kNPart; ipart++){
    for(Int_t icharge = 0; icharge < kNCharge; icharge++){
      Int_t nbinITS = hSpectra[kITS][ipart][icharge]->GetNbinsX();
      for(Int_t ibin = 1; ibin <= nbinITS; ibin++){
	if(hSpectra[kITS][ipart][icharge]->GetBinContent(ibin) < 0 ) {
	  hSpectra[kITS][ipart][icharge]->SetBinContent(ibin,0);
	  hSpectra[kITS][ipart][icharge]->SetBinError  (ibin,0);
	}
	//	if(ipart == kProton && (ibin==9 || ibin==10)){ // FIXME
	if(ipart == kProton && (ibin==9)){
	  printf("Kill bin %d (%f - %f GeV/c)for ITS protons\n",ibin,hSpectra[kITS][ipart][icharge]->GetBinLowEdge(ibin),hSpectra[kITS][ipart][icharge]->GetBinLowEdge(ibin)+hSpectra[kITS][ipart][icharge]->GetBinWidth(ibin));
	  hSpectra[kITS][ipart][icharge]->SetBinContent(ibin,0);
	  hSpectra[kITS][ipart][icharge]->SetBinError  (ibin,0);
	}
	// FIXME: uncomment?
	// if(ipart == kKaon && ibin==7){
	//   printf("Kill bin %d (%f - %f GeV/c)for ITS kaons\n",ibin,hSpectra[kITS][ipart][icharge]->GetBinLowEdge(ibin),hSpectra[kITS][ipart][icharge]->GetBinLowEdge(ibin)+hSpectra[kITS][ipart][icharge]->GetBinWidth(ibin));
	//   hSpectra[kITS][ipart][icharge]->SetBinContent(ibin,0);
	//   hSpectra[kITS][ipart][icharge]->SetBinError  (ibin,0);
	// }
	
// 	if ((ipart == kKaon && ibin >= 12) || (ipart == kProton && ibin >= 20)) {
// 	  hSpectra[kITS][ipart][icharge]->SetBinContent(ibin,0);
// 	  hSpectra[kITS][ipart][icharge]->SetBinError  (ibin,0);
// 	}
      }
      
    }
  }

  if(useSecCorrFromDCA){
    TFile* fseccorr = new TFile("./Files/CorrFac-SecProtons3Meth-ITSsa.root");
    TH1F* hcorrp=(TH1F*)fseccorr->Get("hSecPCorrFromDCAexpo");
    TH1F* hcorrpbar=(TH1F*)fseccorr->Get("hSecPbarCorrFromDCAexpo");
    hSpectra[kITS][kProton][kPos]->Multiply(hcorrp);
    hSpectra[kITS][kProton][kNeg]->Multiply(hcorrpbar);
    fseccorr->Close();
  }

  // Correct for the eta bug (see mails around 30/01/2011)
  if(applyEtaCutCorrection) {
    TH1F * hCorrectionITSsaEtaCut[kNPart][kNCharge] = {0};
    f= TFile::Open("./Files/ITSsaEtaCutCorr.root");
    hCorrectionITSsaEtaCut[kPion  ][kPos] = (TH1F*) f->Get("ITSsaPos0");
    hCorrectionITSsaEtaCut[kKaon  ][kPos] = (TH1F*) f->Get("ITSsaPos1");
    hCorrectionITSsaEtaCut[kProton][kPos] = (TH1F*) f->Get("ITSsaPos2");
    hCorrectionITSsaEtaCut[kPion  ][kNeg] = (TH1F*) f->Get("ITSsaNeg0");
    hCorrectionITSsaEtaCut[kKaon  ][kNeg] = (TH1F*) f->Get("ITSsaNeg1");
    hCorrectionITSsaEtaCut[kProton][kNeg] = (TH1F*) f->Get("ITSsaNeg2");
    for(Int_t ipart = 0; ipart < kNPart; ipart++){
      for(Int_t icharge = 0; icharge < kNCharge; icharge++){
	Int_t nbin = hSpectra[kITS][ipart][icharge]->GetNbinsX();
	for(Int_t ibin = 0; ibin < nbin; ibin++){
	  Int_t content = hSpectra[kITS][ipart][icharge]->GetBinContent(ibin) * hCorrectionITSsaEtaCut[ipart][icharge]->GetBinContent(ibin);
	  Int_t error   = hSpectra[kITS][ipart][icharge]->GetBinError(ibin) * hCorrectionITSsaEtaCut[ipart][icharge]->GetBinContent(ibin);
	  hSpectra[kITS][ipart][icharge]->SetBinContent(ibin, content);
	  hSpectra[kITS][ipart][icharge]->SetBinError  (ibin, error);
	}
      }    
    }
  }
    
  // Load ITS sa systematics, only pt dependent ones
  //  f = TFile::Open("./Files/ITSsa-systematics_20100930.root");
  f = TFile::Open("./Files/ITSsa-systematics-20101014.root");
  for(Int_t ipart = 0; ipart < kNPart; ipart++){
      for(Int_t icharge = 0; icharge < kNCharge; icharge++){
	AliPWGHistoTools::AddHisto(hSystError[kITS][ipart][icharge], (TH1*) gDirectory->Get(Form("hSystTot%s%s",chargeFlag[icharge],partFlag[ipart]))); // Using total error
      }
    }

  // ITS + TPC (Marek)
  //f = TFile::Open("./Files/SpectraCorrectedITSBeforeProtons20100720.root");
  //  f = TFile::Open("./Files/SpectraCorrectedITSBeforeProtons14092010new.root");
  f = TFile::Open("./Files/SpectraCorrectedITSBeforeProtonsstat29102010.root");
  TList * list = (TList*) gDirectory->Get("output");
  hSpectra[kITSTPC][kPion]  [kPos]= (TH1F*) list->FindObject("Pions");
  hSpectra[kITSTPC][kKaon]  [kPos]= (TH1F*) list->FindObject("Kaons");
  hSpectra[kITSTPC][kProton][kPos]= (TH1F*) list->FindObject("Protons");
  hSpectra[kITSTPC][kPion]  [kNeg]= (TH1F*) list->FindObject("AntiPions");
  hSpectra[kITSTPC][kKaon]  [kNeg]= (TH1F*) list->FindObject("AntiKaons");
  hSpectra[kITSTPC][kProton][kNeg]= (TH1F*) list->FindObject("AntiProtons");
  // Systematics
  f = TFile::Open("./Files/SpectraCorrectedITSBeforeProtonssys29102010.root");
  list = (TList*) gDirectory->Get("output");
  hSystError[kITSTPC][kPion]  [kPos]= (TH1F*) AliPWGHistoTools::GetRelativeError((TH1*) list->FindObject("Pions"));
  hSystError[kITSTPC][kKaon]  [kPos]= (TH1F*) AliPWGHistoTools::GetRelativeError((TH1*) list->FindObject("Kaons"));
  hSystError[kITSTPC][kProton][kPos]= (TH1F*) AliPWGHistoTools::GetRelativeError((TH1*) list->FindObject("Protons"));
  hSystError[kITSTPC][kPion]  [kNeg]= (TH1F*) AliPWGHistoTools::GetRelativeError((TH1*) list->FindObject("AntiPions"));
  hSystError[kITSTPC][kKaon]  [kNeg]= (TH1F*) AliPWGHistoTools::GetRelativeError((TH1*) list->FindObject("AntiKaons"));
  hSystError[kITSTPC][kProton][kNeg]= (TH1F*) AliPWGHistoTools::GetRelativeError((TH1*) list->FindObject("AntiProtons"));



  // TPC
  f = new TFile("./Files/protonSpectra_20100615.root");
  hSpectra[kTPC][kProton][kPos]= AliPWGHistoTools::GetHistoFromGraph((TGraphErrors*)f->Get("protonPosClassic"),htemplate);
  hSpectra[kTPC][kProton][kNeg]= AliPWGHistoTools::GetHistoFromGraph((TGraphErrors*)f->Get("protonNegClassic"),htemplate);
  f = new TFile("./Files/pionSpectra_20100615.root");
  hSpectra[kTPC][kPion][kPos]= AliPWGHistoTools::GetHistoFromGraph((TGraphErrors*)f->Get("pionPosClassic"),htemplate);
  hSpectra[kTPC][kPion][kNeg]= AliPWGHistoTools::GetHistoFromGraph((TGraphErrors*)f->Get("pionNegClassic"),htemplate);
  //  f = new TFile("./Files/kaonSpectra_20100615.root");
  f = new TFile("./Files/kaonsTPCcorr_2010_08_31.root");
  hSpectra[kTPC][kKaon][kPos]= AliPWGHistoTools::GetHistoFromGraph((TGraphErrors*)f->Get("kaonPosClassic"),htemplate);
  hSpectra[kTPC][kKaon][kNeg]= AliPWGHistoTools::GetHistoFromGraph((TGraphErrors*)f->Get("kaonNegClassic"),htemplate);

  // Clean UP TPC spectra, removing unwanted points
  cout << "Cleaning Up TPC spectra" << endl;
  nbin =  hSpectra[kTPC][kKaon][kPos]->GetNbinsX();
  for(Int_t ibin = 0; ibin < nbin; ibin++){
    Float_t pt =  hSpectra[kTPC][kKaon][kPos]->GetBinCenter(ibin);
    if (pt > 0.45 || pt<0.25) {
      for(Int_t icharge = 0; icharge < kNCharge; icharge++){
      	hSpectra[kTPC][kKaon][icharge]->SetBinContent(ibin,0);
      	hSpectra[kTPC][kKaon][icharge]->SetBinError  (ibin,0);	
      }      
    }
    // if (pt < 0.25) {
    //   hSpectra[kTPC][kKaon][kNeg]->SetBinContent(ibin,0);
    //   hSpectra[kTPC][kKaon][kNeg]->SetBinError  (ibin,0);	
    // }
    if (pt < 0.45) {
      for(Int_t icharge = 0; icharge < kNCharge; icharge++){
	hSpectra[kTPC][kProton][icharge]->SetBinContent(ibin,0);
	hSpectra[kTPC][kProton][icharge]->SetBinError  (ibin,0);	
      }      
    }
  }
  
  // Load TPC systematics
  cout << "WARNING: TPC SYST FOR NEGATIVES TO BE CORRECTED" << endl;
  f = TFile ::Open("./Files/pionsSystSum.root");
  AliPWGHistoTools::AddHisto(hSystError[kTPC][kPion][kPos],(TH1*) gDirectory->Get("pionsSyst"));
  AliPWGHistoTools::AddHisto(hSystError[kTPC][kPion][kNeg],(TH1*) gDirectory->Get("pionsSyst"));
  f = TFile ::Open("./Files/kaonsSystSum.root");
  AliPWGHistoTools::AddHisto(hSystError[kTPC][kKaon][kPos],(TH1*) gDirectory->Get("kaonsSyst"));
  AliPWGHistoTools::AddHisto(hSystError[kTPC][kKaon][kNeg],(TH1*) gDirectory->Get("kaonsSyst"));
  f = TFile ::Open("./Files/ProtonsSystSum.root");
  AliPWGHistoTools::AddHisto(hSystError[kTPC][kProton][kPos],(TH1*) gDirectory->Get("ProtonsSyst"));
  AliPWGHistoTools::AddHisto(hSystError[kTPC][kProton][kNeg],(TH1*) gDirectory->Get("ProtonsSyst"));
    
  // K0s
  f = new TFile ("./Files/PtSpectraCorrectedK0sOff_20100803.root");
  //  hSpectra[kK0][kKaon][kPos] = (TH1F*) AliPWGHistoTools::GetdNdPtFromOneOverPt((TH1*) gDirectory->Get("hSpectraOff")); 
  hSpectra[kK0][kKaon][kPos] = (TH1F*) gDirectory->Get("hSpectraOff"); 
  //  hSpectra[kK0][kKaon][kPos]->Scale(2*TMath::Pi());
  //  hSpectra[kK0][kKaon][kPos]->Scale(1./272463);
  hSpectra[kK0][kKaon][kNeg] = hSpectra[kK0][kKaon][kPos];

  // Kinks
  //  f = new TFile ("./Files/PtAllKaonKinkRap6Apr24.root");
  //  f = new TFile ("./Files/PtKaonKinkJune13AllPN_20100615.root");
  //  f = new TFile ("./Files/KaonKinkJun16PhySel2N.root");
  f = new TFile ("./Files/KaonKinkJun16PhySel2NRebPass6.root");
  hSpectra[kKinks][kKaon][kPos] = (TH1F*)gDirectory->Get("fptallKPA");
  hSpectra[kKinks][kKaon][kNeg] = (TH1F*)gDirectory->Get("fptallKNA");
  // hSpectra[kKinks][kKaon][kPos]->Scale(0.5/0.7); // different rapidity range for kinks
  // hSpectra[kKinks][kKaon][kNeg]->Scale(0.5/0.7); // different rapidity range for kinks
  // hSpectra[kKinks][kKaon][kPos]->Scale(276004./263345.); // different N of events
  // hSpectra[kKinks][kKaon][kNeg]->Scale(276004./263345.); // different N of events
  // This normalization if for KaonKinkJun16PhySel2N
  // hSpectra[kKinks][kKaon][kPos]->Scale(1./303214); // different N of events
  // hSpectra[kKinks][kKaon][kNeg]->Scale(1./303214); // different N of events
  // This normalization is for KaonKinkJun16PhySel2NRebPass6
  hSpectra[kKinks][kKaon][kPos]->Scale(1./260432.26); // different N of events
  hSpectra[kKinks][kKaon][kNeg]->Scale(1./260432.26); // different N of events

  // Apply correction factors
  // Secondaries for protons

  //f = new TFile ("./Files/corrFactorProtons_20100615.root");
  f = new TFile ("./Files/corrFactorProtons_2010_09_24.root");
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
	    
	    if (correction) {// If the bin is empty this is a  0
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
    // common file for Kaons
    TFile *fFlukakaons = TFile::Open("./Files/correctionForCrossSection.321.root");
    TH1D * hCorrFlukakaon[kNCharge];
    hCorrFlukakaon[kPos] = (TH1D*)fFlukakaons->Get("gHistCorrectionForCrossSectionParticles");
    hCorrFlukakaon[kNeg] = (TH1D*)fFlukakaons->Get("gHistCorrectionForCrossSectionAntiParticles");

    for(Int_t idet = 0; idet < kNDet; idet++){
      if( idet != kITS) continue; // comment to use fluka for kaons for all dets
      if (idet == kTOF) continue; // TOF already corrected
      for(Int_t icharge = 0; icharge < kNCharge; icharge++){
	Int_t ipart = kKaon;
	TH1 * h = hSpectra[idet][ipart][icharge]; // only ITS sa
	if (h){
	  Int_t nbins = h->GetNbinsX();
	  Int_t nbinsy=hCorrFlukakaon[icharge]->GetNbinsY();
	  for(Int_t ibin = 0; ibin < nbins; ibin++){
	    Float_t pt = h->GetBinCenter(ibin);
	    Float_t minPtCorrection = hCorrFlukakaon[icharge]->GetXaxis()->GetBinLowEdge(1);
	    Float_t maxPtCorrection = hCorrFlukakaon[icharge]->GetXaxis()->GetBinLowEdge(nbinsy+1);
	    if (pt < minPtCorrection) pt = minPtCorrection+0.0001;
	    if (pt > maxPtCorrection) pt = maxPtCorrection;
	    Float_t correction = hCorrFlukakaon[icharge]->GetBinContent(hCorrFlukakaon[icharge]->GetXaxis()->FindBin(pt));
	    if (correction != 0) {// If the bin is empty this is a  0
	      h->SetBinContent(ibin,h->GetBinContent(ibin)*correction);
	      h->SetBinError  (ibin,h->GetBinError  (ibin)*correction);
	    } else if (h->GetBinContent(ibin) > 0) { // If we are skipping a non-empty bin, we notify the user
	      cout << "Fluka/GEANT: Not correcting bin "<<ibin << " for protons secondaries, ITS, " << chargeFlag[icharge] << endl;
	      cout << " Bin content: " << h->GetBinContent(ibin)  << endl;
	    }
	  }
	}
      }
    }
 
    // PROTONS

    // ITS specific file for protons/antiprotons
    TFile* fITS = new TFile ("./Files/correctionForCrossSectionITS_20100719.root");
    TH2D * hCorrFlukaITS[kNCharge];
    hCorrFlukaITS[kPos] = (TH2D*)fITS->Get("gHistCorrectionForCrossSectionProtons");
    hCorrFlukaITS[kNeg] = (TH2D*)fITS->Get("gHistCorrectionForCrossSectionAntiProtons");
    
    for(Int_t icharge = 0; icharge < kNCharge; icharge++){
      Int_t ipart = kProton;
      TH1 * h = hSpectra[kITS][ipart][icharge]; // only ITS sa
      if (h){
	Int_t nbins = h->GetNbinsX();
	Int_t nbinsy=hCorrFlukaITS[icharge]->GetNbinsY();
	for(Int_t ibin = 0; ibin < nbins; ibin++){
	  Float_t pt = h->GetBinCenter(ibin);
	  Float_t minPtCorrection = hCorrFlukaITS[icharge]->GetYaxis()->GetBinLowEdge(1);
	  Float_t maxPtCorrection = hCorrFlukaITS[icharge]->GetYaxis()->GetBinLowEdge(nbinsy+1);
	  if (pt < minPtCorrection) pt = minPtCorrection+0.0001;
	  if (pt > maxPtCorrection) pt = maxPtCorrection;
	  Float_t correction = hCorrFlukaITS[icharge]->GetBinContent(1,hCorrFlukaITS[icharge]->GetYaxis()->FindBin(pt));
	  if (correction != 0) {// If the bin is empty this is a  0
	    h->SetBinContent(ibin,h->GetBinContent(ibin)*correction);
	    h->SetBinError  (ibin,h->GetBinError  (ibin)*correction);
	  } else if (h->GetBinContent(ibin) > 0) { // If we are skipping a non-empty bin, we notify the user
	    cout << "Fluka/GEANT: Not correcting bin "<<ibin << " for protons secondaries, ITS, " << chargeFlag[icharge] << endl;
	    cout << " Bin content: " << h->GetBinContent(ibin)  << endl;
	  }
	}
      }
    }
      
    

    //f = new TFile ("./Files/correctionForCrossSection_20100615.root");
    f = new TFile ("./Files/correctionForCrossSection_20100920.root");
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
	    //	    hCorrFluka[icharge]->Draw();
	    TH1D * htmpFluka=hCorrFluka[icharge]->ProjectionY();
	    Float_t maxPtCorrection = AliPWGHistoTools::GetHighestNotEmptyBinEdge(htmpFluka);//->GetYaxis()->GetBinLowEdge(hCorrFluka[icharge]->GetNbinsY()+1);
	    if (pt < minPtCorrection) pt = minPtCorrection+0.0001;
	    if (pt > maxPtCorrection) pt = maxPtCorrection-0.0001;
	    Float_t correction = hCorrFluka[icharge]->GetBinContent(1,hCorrFluka[icharge]->GetYaxis()->FindBin(pt));
	    
	    // already in the efficiency correction (F. Noferini)
	    if (idet == kTOF) {
	      if (icharge == kNeg)  correction = 1; // antiprotons already corrected in efficiency
	      // Scale correction for the different material budget. Recipe by Francesco Noferini
	      else {
		correction = TMath::Power(correction,0.07162/0.03471);
	      }
	    }	    
	    if (correction != 0) {// If the bin is empty this is a  0
	      
	      h->SetBinContent(ibin,h->GetBinContent(ibin)*correction);
	      h->SetBinError  (ibin,h->GetBinError  (ibin)*correction);
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
	    TH1F * htmp = (TH1F*) AliPWGHistoTools::GetOneOverPtdNdPt(hSpectra[idet][ipart][icharge]);
	    hSpectra[idet][ipart][icharge] = (TH1F*)AliPWGHistoTools::GetdNdmtFromdNdpt(htmp,mass[ipart]);
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


  // Create fake weights for the mean; To be update once we have syst errors
  // Using syste from table in paper. It would be better to have this as a function of pt.
  TH1F * hWeights[3][kNPart];
  const Double_t kWeights[3][kNPart] =  
    // {{4,  3,  10.2},  // TPC
    //  {4.1,8.8,7.0 },  //TOF
    //  {4.5,5.6,7.0 }}; // ITS
    {{0.1,0.1,0.1},  // TPC
     {0.1,0.1,0.1},  //TOF
     {0.1,0.1,0.1}}; // ITS
  for(Int_t ipart = 0; ipart <= kNPart ; ipart++){
    for(Int_t idet = 0; idet <= kITS ; idet++){
      hWeights[idet][ipart] = (TH1F*) hSpectra[idet][ipart][kPos]->Clone();
      Int_t nbinW = hWeights[idet][ipart]->GetNbinsX();
      for(Int_t ibin = 1; ibin <= nbinW; ibin++){
	hWeights[idet][ipart]->SetBinError(ibin, kWeights[idet][ipart]);
      }    
    }
  }
  const Double_t scaleDet[] = {1.00,1.00,1.00}; // During the weekly meeting 19/08/2010 we decided not to apply any scaling.
  //  const Double_t scaleDet[] = {0.98,1,1.00}; // Scaling factor for the different detectors. Estimated from ratios, it has an estimated uncertainty of ~2% 
  //  const Double_t scaleDet[] = {0.88,1,0.88}; // Scaling factor for the different detectors. Estimated from ratios, it has an estimated uncertainty of ~2% 


  // Combine detectors
  for(Int_t ipart = 0; ipart < kNPart; ipart++){
    for(Int_t icharge = 0; icharge < kNCharge; icharge++){
      TH1F * htemplLocal = htemplate; // If we are converting to 1/mt dNdmt we need to convert the template as well...
      
      if (convertToMT) {
	TH1F * htmp = (TH1F*) AliPWGHistoTools::GetOneOverPtdNdPt(htemplate);
	htemplLocal = (TH1F*)AliPWGHistoTools::GetdNdmtFromdNdpt(htmp,mass[ipart]);

      }
      hSpectra[kCombTOFTPC][ipart][icharge] = AliPWGHistoTools::CombineHistos(hSpectra[kTOF][ipart][icharge],
									hSpectra[kTPC][ipart][icharge],
									htemplLocal,1.);;

      hSpectra[kCombAll][ipart][icharge]   = 
	AliPWGHistoTools::Combine3HistosWithErrors(hSpectra[kITS][ipart][icharge],  // Histos to combine
					     hSpectra[kTPC][ipart][icharge],
					     hSpectra[kTOF][ipart][icharge],
					     hSystError[kITS][ipart][icharge], // Errors (weights) used for the average
					     hSystError[kTPC][ipart][icharge],
					     hSystError[kTOF][ipart][icharge],
					     // hWeights[kITS][ipart],
					     // hWeights[kTPC][ipart],
					     // hWeights[kTOF][ipart],
					     htemplLocal,1,   // template, take statistical error from TPC in overlap
					     scaleDet[kITS],  // relative scaling
					     scaleDet[kTPC],
					     scaleDet[kTOF],
					     (TH1**)&hSystError[kCombAll][ipart][icharge], // put combined syst error here
					     1 // weights histos contain error in bin content
					     );
//       if (convertToMT) {
// 	TH1F * htmp = (TH1F*) AliPWGHistoTools::GetOneOverPtdNdPt(hSpectra[kCombTOFTPC][ipart][icharge]);
// 	hSpectra[kCombTOFTPC][ipart][icharge] = (TH1F*)AliPWGHistoTools::GetdNdmtFromdNdpt(htmp,mass[ipart]);
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
	  if(idet != kKinks && idet != kK0){ // Kinks use a different run list, k0s normalized by Helene
	    if(idet!=kITSTPC) hSpectra[idet][ipart][icharge]->Scale(1.*effPhysSel[ipart]/278366.15); // Scale PhysSel efficiency (computed by Alexander)
	    else hSpectra[idet][ipart][icharge]->Scale(1./278366.15); // Scale PhysSel efficiency (computed by Alexander)
	  }
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
	    TH1F * htmp = (TH1F*)AliPWGHistoTools::GetOneOverPtdNdPt(hSpectraMC[itune][ipart][icharge]);
	    hSpectraMC[itune][ipart][icharge] = (TH1F*)AliPWGHistoTools::GetdNdmtFromdNdpt(htmp,mass[ipart]);
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
      AliPWGHistoTools::GetResiduals(hSpectra[kITS][ipart][icharge], f, &hres1, &hres2);

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
      Float_t x = AliPWGHistoTools::GetLowestNotEmptyBinEdge(hSpectra[kCombTOFTPC][ipart][icharge]);
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

  Int_t chargeLoop = sumCharge ? 1 : 2; 
  for(Int_t icharge = 0; icharge < chargeLoop; icharge++){
 
    // Draw with models
    for(Int_t ipart = 0; ipart < kNPart; ipart++){
      // Multipad canvas
      TString chargeFlagLocal = chargeFlag[icharge];
      if(sumCharge) chargeFlagLocal += chargeFlag[1];
      TCanvas * c1 = new TCanvas(TString("cSpectra")+partFlag[ipart]+chargeFlagLocal,TString("cSpectra")+partFlag[ipart]+chargeFlagLocal,700,700);
      TPad *p1 = new TPad(TString("p1")+partFlag[ipart]+chargeFlagLocal, "p1", 0.0, 0.35, 1.0,  0.95, 0, 0, 0);
      p1->SetBottomMargin(0);
      p1->Draw();
      
      TPad *p2 = new TPad(TString("p2")+partFlag[ipart]+chargeFlagLocal, "p2", 0.0, 0.05, 1.0,  0.35, 0, 0, 0);
      p2->SetTopMargin(0);
      p2->SetBottomMargin(0.3);
      p2->Draw();

      Float_t scaleFonts = (0.95-0.3)/(0.3-0.05);

      // Draw spectra
      p1->cd();
      p1->SetLogy();
      Float_t maxy = sumCharge ? 10 : 5;
      TH2F * hempty = new TH2F(TString("hempty")+long(icharge)+long(ipart),"hempty",100,0.,3, 100, 0.0015,maxy);
      hempty->SetXTitle("p_{t} (GeV/c)");
      hempty->SetYTitle("1/N_{ev} d^{2}N / dydp_{t} (GeV/c)^{-1}");
      hempty->Draw();
      c1->SetLogy();
      


      TLegend * l =new TLegend(0.502874, 0.493056, 0.892241, 0.904762);
      l->SetFillColor(kWhite);      
      if(sumCharge) 	{
	hSpectra[iCombInStudy][ipart][kPos]->Add(hSpectra[iCombInStudy][ipart][kNeg]); // Draw pos+neg
      }
      hSpectra[iCombInStudy][ipart][icharge]->Draw("same");
      l->AddEntry(hSpectra[kTOF][ipart][icharge],TString ("Data"));
      for(Int_t itune = 0; itune < kNTunes; itune++){      
	if(sumCharge)	hSpectraMC[itune][ipart][icharge]->Add(hSpectraMC[itune][ipart][kNeg]); // Draw pos+neg;    
	l->AddEntry(hSpectraMC[itune][ipart][icharge],mcTuneName[itune]);
	hSpectraMC[itune][ipart][icharge]->SetLineWidth(2);    
	AliPWGHistoTools::GetGraphFromHisto(hSpectraMC[itune][ipart][icharge])->Draw("CX");
      }
      l->Draw("same");
     
      TLatex * tex = new TLatex(0.6712643,2.353486,sumCharge ? Form("%s+%s",partLabel[ipart][icharge],partLabel[ipart][1]) : partLabel[ipart][icharge]);
      tex->SetTextFont(42);
      tex->SetTextSize(0.07936508);
      tex->SetLineWidth(2);
      tex->Draw();

      // Draw ratio
      p2->cd();
      Float_t maxyr = sumCharge ? 3.99 : 2.99;
      TH2F * hemptyr = new TH2F(TString("hemptyratio")+long(icharge)+long(ipart),"hempty",100,0.,3, 100, 0.01,maxyr);
      hemptyr->SetXTitle("p_{t} (GeV/c)");
      hemptyr->SetYTitle("Data/MC");
      hemptyr->GetXaxis()->SetLabelSize(0.04*scaleFonts);
      hemptyr->GetYaxis()->SetLabelSize(0.04*scaleFonts);
      hemptyr->GetYaxis()->SetTitleSize(0.05*scaleFonts);
      hemptyr->GetYaxis()->SetTitleOffset(1.4/scaleFonts);
      hemptyr->GetXaxis()->SetTitleSize(0.05*scaleFonts);
      hemptyr->GetXaxis()->SetTitleOffset(1.05);
      hemptyr->SetTickLength(0.03*scaleFonts, "X");
      hemptyr->SetTickLength(0.02*scaleFonts, "Y");
      //      hemptyr->GetXaxis()->SetTitleOffset(1.4/scaleFonts);
      hemptyr->GetYaxis()->SetNdivisions(5);
      hemptyr->Draw("");

      AliPWGFunc fm;
      for(Int_t itune = 0; itune < kNTunes; itune++){      
	TF1* f = fm.GetHistoFunc(hSpectraMC[itune][ipart][icharge], TString("f")+mcTuneName[itune]);

	//	l->AddEntry(hSpectraMC[itune][ipart][icharge],mcTuneName[itune]);
	TH1F* hRatio = AliPWGHistoTools::DivideHistoByFunc(hSpectra[iCombInStudy][ipart][icharge],f);
	hRatio->SetLineStyle(hSpectraMC[itune][ipart][icharge]->GetLineStyle());
	hRatio->SetLineColor(hSpectraMC[itune][ipart][icharge]->GetLineColor());
	hRatio->SetLineWidth(hSpectraMC[itune][ipart][icharge]->GetLineWidth());
	hRatio->SetMarkerStyle(0);
	hRatio->Draw("same");
	//	AliPWGHistoTools::GetGraphFromHisto(hRatio)->Draw("CX");
      }


      // print
      PrintCanvas(c1,printFormats);
    }
  }



}

void DrawAllAndKaons() {

  

  //  gROOT->LoadMacro("ALICEWorkInProgress.C");

  //  gStyle->SetOptFit(0);
  gStyle->SetStatX(0.9);
  gStyle->SetStatY(0.88);
  gStyle->SetStatW(0.4);
  gStyle->SetStatH(0.1);

  TH1F * hKaonsAllTPCTOF = (TH1F*) hSpectra[iCombInStudy][kKaon][kPos]->Clone();
  hKaonsAllTPCTOF->Add(hSpectra[iCombInStudy][kKaon][kNeg]);
  
  TH1F * hK0Scaled    = (TH1F*) hSpectra[kK0][kKaon][kPos]->Clone();
  hK0Scaled->Add(hSpectra[kK0][kKaon][kPos]);

  hSpectra[kKinks][kKaon][kPos]->SetMarkerStyle(25);
  hSpectra[kKinks][kKaon][kPos]->SetLineColor(4);
  hSpectra[kKinks][kKaon][kPos]->SetStats(0);
  TH1F * hKinksAll = (TH1F*) hSpectra[kKinks][kKaon][kPos]->Clone();
  hKinksAll->Add(hSpectra[kKinks][kKaon][kNeg]);
  
  TCanvas * c1 = new TCanvas("cKaons","cKaons",700,700);
  c1->SetLogy();
  TH2F * hempty = new TH2F("hempty_allkaons","hempty",100,0.,3, 100, 1e-3,6);
  hempty->SetXTitle("p_{t} (GeV/c)");
  //  hempty->SetYTitle("dN / dp_{t} (A.U.)");
  hempty->SetYTitle("1/N_{ev} d^{2}N / dydp_{t} (GeV/c)^{-1}");
  hempty->Draw();
  hK0Scaled->Draw("same");
  hKaonsAllTPCTOF->Draw("same");
  hKinksAll->Draw("same");
  

  TLegend * leg = new TLegend(0.456897, 0.71875, 0.892241, 0.936012,NULL,"brNDC");
  //    leg->SetBorderSize(0);
//    leg->SetLineColor(1);
//    leg->SetLineStyle(1);
//    leg->SetLineWidth(1);
//    leg->SetFillColor(19);
    leg->SetFillColor(0);
   TLegendEntry *entry=leg->AddEntry(hKaonsAllTPCTOF,"K^{+} + K^{-}, ITS+TPC+TOF ","lpf");
   entry=leg->AddEntry(hK0Scaled,"K^{0} #times 2","lpf");
   entry=leg->AddEntry(hKinksAll,"K^{+} + K ^{-}, Kinks","lpf");
   leg->Draw();

   //   ALICEWorkInProgress(c1,today.Data(),"#splitline{ALICE Prelimiary}{Statistical Error Only}");
   // TLatex * tex = new TLatex(0.2120805,0.01288336,"Statistical error only");
   // tex->SetTextColor(2);
   // tex->SetTextFont(42);
   // tex->SetTextSize(0.03496503);
   // tex->Draw();

   c1->Update();
   PrintCanvas(c1, printFormats);

  // Draw all "stable" hadrons
   Bool_t divideCanvas=kFALSE; 
   TH1F * hSpectraSystError[kNHist][kNPart][kNCharge];// used to store spectra with syst error only
       
  for(Int_t icharge = 0; icharge < kNCharge; icharge++){
    TCanvas * c1h = 0;
    if(divideCanvas) c1h = new TCanvas(TString("cAll_")+chargeFlag[icharge],TString("cAll_")+chargeFlag[icharge],1200,500);
    else {
      c1h = new TCanvas(TString("cAll_")+chargeFlag[icharge],TString("cAll_")+chargeFlag[icharge],700,700);
      c1h->SetLogy();
    }
    c1h->SetLeftMargin(0.14);
    if(divideCanvas) {
      c1h->Divide(3,1);
      c1h->cd(1);
    }
    TH2F * hemptyLoc = new TH2F(TString("hempty")+long(icharge),"hempty",100,0.,4, 100, 1e-2,5);
    hemptyLoc->SetXTitle("p_{t} (GeV/c)");
    hemptyLoc->SetYTitle("1/N_{ev} d^{2}N / dydp_{t} (GeV/c)^{-1}");
    hemptyLoc->GetYaxis()->SetTitleOffset(1.35);
    hemptyLoc->GetXaxis()->SetTitleOffset(1.1);
    hemptyLoc->GetXaxis()->SetRangeUser(0,1);
    hemptyLoc->Draw();    

    for(Int_t ipart = 0; ipart < kNPart; ipart++) {                  
      if(ipart == kPion)         leg = new TLegend(0.558908, 0.78125, 0.889368, 0.924107, NULL,"brNDC");
      else if (ipart == kKaon)   leg = new TLegend(0.16523, 0.177083, 0.50431, 0.31994, NULL,"brNDC");
      else if (ipart == kProton) leg = new TLegend(0.51, 0.177083, 0.84, 0.31994, NULL,"brNDC");
      leg->SetFillColor(0);
      if(divideCanvas) {
	cout << "CD!!" << ipart+1 << endl;	
	c1h->cd(ipart+1);
	TH2F * hClone = new TH2F(*hemptyLoc);
	hClone->GetXaxis()->SetRangeUser(0,1);
	if(ipart == kPion)   hClone->GetYaxis()->SetRangeUser(0,3.5);
	if(ipart == kKaon)   hClone->GetYaxis()->SetRangeUser(0,0.28);
	if(ipart == kProton) hClone->GetYaxis()->SetRangeUser(0,0.1);
	hClone->Draw();
      }
      for(Int_t idet = 0; idet <= kITSTPC; idet++){
	//	if (idet == kITS) continue;
	// 	if (idet == kITSTPC) hSpectra[idet][ipart][icharge]->SetMarkerColor(kGreen);
 	hSpectra[idet][ipart][icharge]->SetMarkerStyle(marker[idet]);
	TH1F * hsyst = new TH1F(*hSpectra[idet][ipart][icharge]);
	hsyst->Reset();
	AliPWGHistoTools::GetValueAndError(hsyst, hSpectra[idet][ipart][icharge],hSystError[idet][ipart][icharge],kTRUE);
	//	hSpectra[idet][ipart][icharge]->Draw("same");
	hsyst->Draw("same");
	hSpectraSystError[idet][ipart][icharge] = hsyst;
	leg->AddEntry(hSpectra[idet][ipart][icharge],TString(partLabel[ipart][icharge])+" (" + detFlag[idet]  + ")","lpf");
      }
      //      leg->AddLine();
      leg->Draw();
    } 
    if(divideCanvas) c1h->cd(1);
    //    leg->Draw();
    //    ALICEWorkInProgress(c1h,today.Data(),"#splitline{ALICE Preliminary}{Statistical Error Only}");
    PrintCanvas(c1h,printFormats);
  }
 

  //  Draw ratios 

  // K-/K+ in the different detectors
  TCanvas * cpm=new TCanvas("cpm","Kminus/Kplus",700,700);
  cpm->Divide(2,2);
  cpm->cd(1);
  TH1F* hRatioKPKM_TPC=new TH1F(*(hSpectra[kTPC][kKaon][kNeg]));
  hRatioKPKM_TPC->SetMinimum(0.5);
  hRatioKPKM_TPC->SetMaximum(1.5);
  hRatioKPKM_TPC->Divide(hSpectra[kTPC][kKaon][kPos]);
  hRatioKPKM_TPC->GetYaxis()->SetTitle("K-/K+ (TPC)");
  hRatioKPKM_TPC->Draw();
  cpm->cd(2);
  TH1F* hRatioKPKM_ITS=new TH1F(*(hSpectra[kITS][kKaon][kNeg]));
  hRatioKPKM_ITS->Divide(hSpectra[kITS][kKaon][kPos]);
  hRatioKPKM_ITS->SetMinimum(0.5);
  hRatioKPKM_ITS->SetMaximum(1.5);
  hRatioKPKM_ITS->GetYaxis()->SetTitle("K-/K+ (ITSsa)");
  hRatioKPKM_ITS->Draw("");
  cpm->cd(3);
  TH1F* hRatioKPKM_TOF=new TH1F(*(hSpectra[kTOF][kKaon][kNeg]));
  hRatioKPKM_TOF->Divide(hSpectra[kTOF][kKaon][kPos]);
  hRatioKPKM_TOF->SetMinimum(0.5);
  hRatioKPKM_TOF->SetMaximum(1.5);
  hRatioKPKM_TOF->GetYaxis()->SetTitle("K-/K+ (TOF)");
  hRatioKPKM_TOF->Draw("");
  cpm->cd(4);
  TH1F* hRatioKPKM_ITSTPC=new TH1F(*(hSpectra[kITSTPC][kKaon][kNeg]));
  hRatioKPKM_ITSTPC->Divide(hSpectra[kITSTPC][kKaon][kPos]);
  hRatioKPKM_ITSTPC->SetMinimum(0.5);
  hRatioKPKM_ITSTPC->SetMaximum(1.5);
  hRatioKPKM_ITSTPC->GetYaxis()->SetTitle("K-/K+ (ITS Global)");
  hRatioKPKM_ITSTPC->Draw("");
  

  // ITS/TPC
  TH1F * hRatioITSTPC[kNPart][kNCharge];
  for(Int_t icharge = 0; icharge < kNCharge; icharge++){ // loop over charges
    // Create canvas
    TCanvas * c1h = new TCanvas(TString("cITSTPCRatio_")+chargeFlag[icharge],TString("cITSTPCRatio_")+chargeFlag[icharge],1200,500);
    c1h->Divide(3,1);
    c1h->SetGridy();
    TH2F * hemptyLoc = new TH2F(TString("hemptyR")+long(icharge),"ITSsa/TPC ",100,0.,1., 100, 0.5,1.5);
    hemptyLoc->SetXTitle("p_{t} (GeV/c)");
    hemptyLoc->SetYTitle("ITSsa / TPC");
    // Loop over particles
    for(Int_t ipart = 0; ipart < kNPart; ipart++) {
      // Clone histo
      hRatioITSTPC[ipart][icharge]=new TH1F(*hSpectraSystError[kITS][ipart][icharge]);
      Int_t nBinsITS=hSpectraSystError[kITS][ipart][icharge]->GetNbinsX();
      Int_t nBinsTPC=hSpectraSystError[kTPC][ipart][icharge]->GetNbinsX();
      // Loop over ITS bins, 
      for(Int_t iBin=1; iBin<=nBinsITS; iBin++){
	hRatioITSTPC[ipart][icharge]->SetBinContent(iBin,0.);
	hRatioITSTPC[ipart][icharge]->SetBinError(iBin,0.);
	Float_t lowPtITS=hSpectraSystError[kITS][ipart][icharge]->GetBinLowEdge(iBin);
	Float_t binWidITS=hSpectraSystError[kITS][ipart][icharge]->GetBinWidth(iBin);
	// Loop over TPC bins and find overlapping bins to ITS
	for(Int_t jBin=1; jBin<=nBinsTPC; jBin++){
	  Float_t lowPtTPC=hSpectraSystError[kTPC][ipart][icharge]->GetBinLowEdge(jBin);
	  Float_t binWidTPC=hSpectraSystError[kTPC][ipart][icharge]->GetBinWidth(jBin);
	  if(TMath::Abs(lowPtITS-lowPtTPC)<0.001 && TMath::Abs(binWidTPC-binWidITS)<0.001){
	    Float_t numer=hSpectraSystError[kITS][ipart][icharge]->GetBinContent(iBin);
	    Float_t denom=hSpectraSystError[kTPC][ipart][icharge]->GetBinContent(jBin);
	    Float_t enumer=hSpectraSystError[kITS][ipart][icharge]->GetBinError(iBin);
	    Float_t edenom=hSpectraSystError[kTPC][ipart][icharge]->GetBinError(jBin);
	    Double_t ratio=0.;
	    Double_t eratio=0.;
	    if(numer>0. && denom>0.){
	      ratio=numer/denom;
	      eratio=ratio*TMath::Sqrt((enumer/numer)*(enumer/numer)+(edenom/denom)*(edenom/denom));
	    }
	    hRatioITSTPC[ipart][icharge]->SetBinContent(iBin,ratio);
	    hRatioITSTPC[ipart][icharge]->SetBinError(iBin,eratio);
	    break;
	  }
	}
      }
      c1h->cd(ipart+1);
      // hemptyLoc->SetStats(1);
      // hemptyLoc->Draw(); 
      hRatioITSTPC[ipart][icharge]->SetYTitle("ITSsa/TPC");
      hRatioITSTPC[ipart][icharge]->SetStats(1);
      hRatioITSTPC[ipart][icharge]->GetYaxis()->SetRangeUser(0.5,1.5);
      hRatioITSTPC[ipart][icharge]->Draw("");
      TF1 * f = new TF1("fpol0","[0]");
      f->SetParameter(0,1);
      hRatioITSTPC[ipart][icharge]->Fit(f,"","same");
       
    }
    PrintCanvas(c1h,printFormats);
  }

  // ITS/TPC
  TH1F * hRatioITSGlobalTPC[kNPart][kNCharge];
  for(Int_t icharge = 0; icharge < kNCharge; icharge++){ // loop over charges
    // Create canvas
    TCanvas * c1h = new TCanvas(TString("cITSGlobalTPCRatio_")+chargeFlag[icharge],TString("cITSGlobalTPCRatio_")+chargeFlag[icharge],1200,500);
    c1h->Divide(3,1);
    c1h->SetGridy();
    TH2F * hemptyLoc = new TH2F(TString("hemptyR")+long(icharge),"ITS Global/TPC ",100,0.,1., 100, 0.5,1.5);
    hemptyLoc->SetXTitle("p_{t} (GeV/c)");
    hemptyLoc->SetYTitle("ITS Global / TPC");
    // Loop over particles
    for(Int_t ipart = 0; ipart < kNPart; ipart++) {
      // Clone histo
      hRatioITSGlobalTPC[ipart][icharge]=new TH1F(*hSpectraSystError[kITSTPC][ipart][icharge]);
      Int_t nBinsITS=hSpectraSystError[kITSTPC][ipart][icharge]->GetNbinsX();
      Int_t nBinsTPC=hSpectraSystError[kTPC][ipart][icharge]->GetNbinsX();
      // Loop over ITS bins, 
      for(Int_t iBin=1; iBin<=nBinsITS; iBin++){
	hRatioITSGlobalTPC[ipart][icharge]->SetBinContent(iBin,0.);
	hRatioITSGlobalTPC[ipart][icharge]->SetBinError(iBin,0.);
	Float_t lowPtITS=hSpectraSystError[kITSTPC][ipart][icharge]->GetBinLowEdge(iBin);
	Float_t binWidITS=hSpectraSystError[kITSTPC][ipart][icharge]->GetBinWidth(iBin);
	// Loop over TPC bins and find overlapping bins to ITS
	for(Int_t jBin=1; jBin<=nBinsTPC; jBin++){
	  Float_t lowPtTPC=hSpectraSystError[kTPC][ipart][icharge]->GetBinLowEdge(jBin);
	  Float_t binWidTPC=hSpectraSystError[kTPC][ipart][icharge]->GetBinWidth(jBin);
	  if(TMath::Abs(lowPtITS-lowPtTPC)<0.001 && TMath::Abs(binWidTPC-binWidITS)<0.001){
	    Float_t numer=hSpectraSystError[kITSTPC][ipart][icharge]->GetBinContent(iBin);
	    Float_t denom=hSpectraSystError[kTPC][ipart][icharge]->GetBinContent(jBin);
	    Float_t enumer=hSpectraSystError[kITSTPC][ipart][icharge]->GetBinError(iBin);
	    Float_t edenom=hSpectraSystError[kTPC][ipart][icharge]->GetBinError(jBin);
	    Double_t ratio=0.;
	    Double_t eratio=0.;
	    if(numer>0. && denom>0.){
	      ratio=numer/denom;
	      eratio=ratio*TMath::Sqrt((enumer/numer)*(enumer/numer)+(edenom/denom)*(edenom/denom));
	    }
	    hRatioITSGlobalTPC[ipart][icharge]->SetBinContent(iBin,ratio);
	    hRatioITSGlobalTPC[ipart][icharge]->SetBinError(iBin,eratio);
	    break;
	  }
	}
      }
      c1h->cd(ipart+1);
      // hemptyLoc->SetStats(1);
      // hemptyLoc->Draw(); 
      hRatioITSGlobalTPC[ipart][icharge]->SetYTitle("ITS Global/TPC");
      hRatioITSGlobalTPC[ipart][icharge]->SetStats(1);
      hRatioITSGlobalTPC[ipart][icharge]->GetYaxis()->SetRangeUser(0.5,1.5);
      hRatioITSGlobalTPC[ipart][icharge]->Draw("");
      TF1 * f = new TF1("fpol0","[0]");
      f->SetParameter(0,1);
      hRatioITSGlobalTPC[ipart][icharge]->Fit(f,"","same");
       
    }
    PrintCanvas(c1h,printFormats);
  }

  // TOF/TPC
  TH1F * hRatioTOFTPC[kNPart][kNCharge];
  for(Int_t icharge = 0; icharge < kNCharge; icharge++){ // loop over charges
    // create canvas
    TCanvas * c1t = new TCanvas(TString("cTOFTPCRatio_")+chargeFlag[icharge],TString("cTOFTPCRatio_")+chargeFlag[icharge],1200,500);
    c1t->SetGridy();
    c1t->Divide(3,1);
    TH2F * hemptyt = new TH2F(TString("hemptyRt")+long(icharge),"TOF/TPC ",100,0.,1., 100, 0.5,1.5);
    hemptyt->SetXTitle("p_{t} (GeV/c)");
    hemptyt->SetYTitle("TOF / TPC");
    //    hemptyt->Draw();    
    for(Int_t ipart = 0; ipart < kNPart; ipart++) { // loop over particles
      // Clone histo
      hRatioTOFTPC[ipart][icharge]=new TH1F(*hSpectraSystError[kTOF][ipart][icharge]);
      Int_t nBinsTOF=hSpectraSystError[kTOF][ipart][icharge]->GetNbinsX();
      Int_t nBinsTPC=hSpectraSystError[kTPC][ipart][icharge]->GetNbinsX();
      // Loop over TOF bins
      for(Int_t iBin=1; iBin<=nBinsTOF; iBin++){
	hRatioTOFTPC[ipart][icharge]->SetBinContent(iBin,0.);
	hRatioTOFTPC[ipart][icharge]->SetBinError(iBin,0.);
	Float_t lowPtTOF=hSpectraSystError[kTOF][ipart][icharge]->GetBinLowEdge(iBin);
	Float_t binWidTOF=hSpectraSystError[kTOF][ipart][icharge]->GetBinWidth(iBin);
	// Loop over TPC bins and find overlapping bins to ITS
	for(Int_t jBin=1; jBin<=nBinsTPC; jBin++){
	  Float_t lowPtTPC=hSpectraSystError[kTPC][ipart][icharge]->GetBinLowEdge(jBin);
	  Float_t binWidTPC=hSpectraSystError[kTPC][ipart][icharge]->GetBinWidth(jBin);
	  if(TMath::Abs(lowPtTOF-lowPtTPC)<0.001 && TMath::Abs(binWidTPC-binWidTOF)<0.001){
	    Float_t numer=hSpectraSystError[kTOF][ipart][icharge]->GetBinContent(iBin);
	    Float_t denom=hSpectraSystError[kTPC][ipart][icharge]->GetBinContent(jBin);
	    Float_t enumer=hSpectraSystError[kTOF][ipart][icharge]->GetBinError(iBin);
	    Float_t edenom=hSpectraSystError[kTPC][ipart][icharge]->GetBinError(jBin);
	    Double_t ratio=0.;
	    Double_t eratio=0.;
	    if(numer>0. && denom>0.){
	      ratio=numer/denom;
	      eratio=ratio*TMath::Sqrt((enumer/numer)*(enumer/numer)+(edenom/denom)*(edenom/denom));
	    }
	    hRatioTOFTPC[ipart][icharge]->SetBinContent(iBin,ratio);
	    hRatioTOFTPC[ipart][icharge]->SetBinError(iBin,eratio);
	    break;
	  }
	}
      }
      c1t->cd(ipart+1);
      hRatioTOFTPC[ipart][icharge]->SetYTitle("TOF/TPC");
      hRatioTOFTPC[ipart][icharge]->SetStats(1);
      hRatioTOFTPC[ipart][icharge]->GetYaxis()->SetRangeUser(0.5,1.5);
      hRatioTOFTPC[ipart][icharge]->Draw("");
      hRatioTOFTPC[ipart][icharge]->Print();
      TF1 * f = new TF1("fpol0","[0]");
      f->SetParameter(0,1);
      hRatioTOFTPC[ipart][icharge]->Fit(f,"","same");
    }
    PrintCanvas(c1t,printFormats);
  }

}

void DrawWithJacek() {

  //1. Convert spectra to dNdeta and sum
  TH1F * hsum = (TH1F*) htemplate->Clone();
  hsum->Reset();
  Int_t idet= iCombInStudy;
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
  new TFile ("./Files/pt_allcharged_published.root");
  TGraphErrors * gJacek = (TGraphErrors*) gDirectory->Get("dndpt");
  gJacek->SetMarkerStyle(20);
  gJacek->SetMarkerColor(kGray);
  gJacek->Draw("AP");
  gJacek->GetHistogram()->SetXTitle("p_{t} (GeV)");
  gJacek->GetHistogram()->SetYTitle("d^{2}N/dp_{t}d#eta (GeV^{-1})");
  hsum->Draw("same");

  TGraphErrors * gRatio = AliPWGHistoTools::DivideGraphByHisto(gJacek,hsum);  
  new TCanvas();
  gRatio->Draw("AP");

  

}


void DrawRatioToStar() {

  // Star data
  //  gROOT->LoadMacro("StarPPSpectra.C");
  TGraphErrors ** gStar = StarPPSpectra();
  gStar[0]->SetMarkerStyle(kOpenStar);
  gStar[1]->SetMarkerStyle(kFullStar);
  gStar[2]->SetMarkerStyle(kOpenStar);
  gStar[3]->SetMarkerStyle(kFullStar);
  gStar[4]->SetMarkerStyle(kOpenStar);
  gStar[5]->SetMarkerStyle(kFullStar);

  // ALICE, INEL -> NSD
  Double_t scaleYield = 3.58/3.02; // from paper 2
  for(Int_t ipart = 0; ipart < kNPart; ipart++){
    for(Int_t icharge = 0; icharge < kNCharge; icharge++){
      hSpectra[iCombInStudy][ipart][icharge]->Scale(scaleYield);
    }    
  }
  
    
  TCanvas * c1 = new TCanvas("cRatioToStarNeg","cRatioToStarNeg");
  TH2F * hempty = new TH2F(TString("hemptyNeg"),"hemptyNeg",100,0.,1.5, 100, 0.001,1.8);
  hempty->SetXTitle("p_{t} (GeV/c)");
  hempty->SetYTitle("ALICE/STAR (NSD)");
  hempty->Draw();

  TCanvas * c1Comp = new TCanvas("cCompToStarNeg","cCompToStarNeg");
  c1Comp->SetLogy();
  TH2F * hempty2 = new TH2F(TString("hemptyCompNeg"),"hemptyCompNeg",100,0.,1.5, 100, 0.001,10);
  hempty2->SetXTitle("p_{t} (GeV/c)");
  hempty2->SetYTitle("1/N_{ev} d^{2}N / dydp_{t} (GeV/c)^{-1} [NSD]");
  hempty2->Draw();
  
  TLegend *leg = new TLegend(0.6510067,0.1853147,0.8892617,0.4178322,"Negative","brNDC");
  leg->SetBorderSize(0);
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(0);
  leg->SetFillStyle(1001);



  c1->cd();
  TGraphErrors * g ;
  g = AliPWGHistoTools::DivideGraphByHisto(gStar[0],hSpectra[iCombInStudy][kPion][kNeg],1);
  g->SetMarkerStyle(kFullCircle);
  g->SetMarkerColor(kBlack);
  g->Draw("p");
  leg->AddEntry(g,"#pi^{-}","lp");
  g = AliPWGHistoTools::DivideGraphByHisto(gStar[2],hSpectra[iCombInStudy][kKaon][kNeg],1);
  g->SetMarkerStyle(kOpenCircle);
  g->SetMarkerColor(kRed);
  g->Draw("p");
  leg->AddEntry(g,"K^{-}","lp");
  g = AliPWGHistoTools::DivideGraphByHisto(gStar[4],hSpectra[iCombInStudy][kProton][kNeg],1);
  g->SetMarkerStyle(kOpenSquare);
  g->SetMarkerColor(kBlue);
  g->Draw("p");
  leg->AddEntry(g,"#bar{p}","lp");  
  leg->Draw();

  c1Comp->cd();
  gStar[0]->Draw("p");
  hSpectra[iCombInStudy][kPion][kNeg]->Draw("same");
  gStar[2]->Draw("p");
  hSpectra[iCombInStudy][kKaon][kNeg]->Draw("same");
  gStar[4]->Draw("p");
  hSpectra[iCombInStudy][kProton][kNeg]->Draw("same");



  TCanvas * c2 = new TCanvas("cRatioToStarPos","cRatioToStarPos");
  hempty->Draw(); 
  leg = new TLegend(0.6510067,0.1853147,0.8892617,0.4178322,"Positive","brNDC");
  leg->SetBorderSize(0);
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(0);
  leg->SetFillStyle(1001);

  TCanvas * c2Comp = new TCanvas("cCompToStarPos","cCompToStarPos");
  c2Comp->SetLogy();
  hempty2->Draw();

  c2->cd();
 //  TGraphErrors * g ;
  g = AliPWGHistoTools::DivideGraphByHisto(gStar[1],hSpectra[iCombInStudy][kPion][kPos],1);
  g->SetMarkerStyle(kFullCircle);
  g->SetMarkerColor(kBlack);
  g->Draw("p");
  leg->AddEntry(g,"#pi^{+}","lp");
  g = AliPWGHistoTools::DivideGraphByHisto(gStar[3],hSpectra[iCombInStudy][kKaon][kPos],1);
  g->SetMarkerStyle(kOpenCircle);
  g->SetMarkerColor(kRed);
  g->Draw("p");
  leg->AddEntry(g,"K^{+}","lp");
  g = AliPWGHistoTools::DivideGraphByHisto(gStar[5],hSpectra[iCombInStudy][kProton][kPos],1);
  g->SetMarkerStyle(kOpenSquare);
  g->SetMarkerColor(kBlue);
  g->Draw("p");
  leg->AddEntry(g,"p","lp");
  leg->Draw();


  c2Comp->cd();
  gStar[1]->Draw("p");
  hSpectra[iCombInStudy][kPion][kPos]->Draw("same");
  gStar[3]->Draw("p");
  hSpectra[iCombInStudy][kKaon][kPos]->Draw("same");
  gStar[5]->Draw("p");
  hSpectra[iCombInStudy][kProton][kPos]->Draw("same");


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

  // Rebin histograms (hi pt was fluctuating too much in the TOF)
  RebinHistosForRatios();

  // Compute ratios
  TH1F * hPosNegRatio[kNPart];
  
  for(Int_t ipart = 0; ipart < kNPart; ipart++){
    hPosNegRatio[ipart] = (TH1F*) hSpectra[iCombInStudy][ipart][kPos]->Clone();
    hPosNegRatio[ipart]->Divide(hSpectra[iCombInStudy][ipart][kNeg]);
    hPosNegRatio[ipart]->SetYTitle(TString(partLabel[ipart][kPos])+"/"+partLabel[ipart][kNeg]);    
    hPosNegRatio[ipart]->SetMinimum(0.5);
    hPosNegRatio[ipart]->SetMaximum(1.5);
  }
  
  
  TH1F * hKPiRatio = (TH1F*) hSpectra[iCombInStudy][kKaon][kPos]->Clone();
  hKPiRatio->Add(hSpectra[iCombInStudy][kKaon][kNeg]);
  TH1F * htmp = (TH1F*) hSpectra[iCombInStudy][kPion][kPos]->Clone();
  htmp->Add(hSpectra[iCombInStudy][kPion][kNeg]);
  hKPiRatio->Divide(htmp);
  hKPiRatio->SetYTitle("(K^{+}+K^{-})/(#pi^{+}+#pi^{-})");
  // Get ratio of syst errors
  // FIXME: commented combined error
  TH1F * hKPiRatioSyst = new TH1F(*htemplate);;
  AliPWGHistoTools::GetValueAndError(hKPiRatioSyst,hSpectra[iCombInStudy][kKaon][kPos],hSystError[iCombInStudy][kKaon][kPos],kTRUE);
  //  AliPWGHistoTools::GetHistoCombinedErrors(hKPiRatioSyst,hSpectra[iCombInStudy][kKaon][kPos]);
  AliPWGHistoTools::GetValueAndError(htmp,hSpectra[iCombInStudy][kKaon][kNeg],hSystError[iCombInStudy][kKaon][kNeg],kTRUE);
  // //  AliPWGHistoTools::GetHistoCombinedErrors(htmp,hSpectra[iCombInStudy][kKaon][kNeg]);
  hKPiRatioSyst->Add(htmp);
  AliPWGHistoTools::GetValueAndError(htmp,hSpectra[iCombInStudy][kPion][kNeg],hSystError[iCombInStudy][kPion][kNeg],kTRUE);
  //  AliPWGHistoTools::GetHistoCombinedErrors(htmp,hSpectra[iCombInStudy][kPion][kNeg]);
  TH1F * htmp2 = new TH1F(*htemplate);
  AliPWGHistoTools::GetValueAndError(htmp2,hSpectra[iCombInStudy][kPion][kPos],hSystError[iCombInStudy][kPion][kPos],kTRUE);
  //  AliPWGHistoTools::GetHistoCombinedErrors(htmp2,hSpectra[iCombInStudy][kPion][kPos]);
  htmp->Add(htmp2);
  hKPiRatioSyst->Divide(htmp);
  hKPiRatioSyst->SetYTitle("(K^{+}+K^{-})/(#pi^{+}+#pi^{-})");
  // Multiply error for sqrt(2): the syst is not reduced by combining positive and negative charges
  Int_t nbin = hKPiRatioSyst->GetNbinsX();
  for(Int_t ibin = 1; ibin <= nbin; ibin++){
    hKPiRatioSyst->SetBinError(ibin, hKPiRatioSyst->GetBinError(ibin)*TMath::Sqrt(2));
  }
  



  TH1F * hKPiRatioMC[kNTunes];
  if(showMC){
    for(Int_t itune = 0; itune < kNTunes; itune++){
      hKPiRatioMC[itune] = (TH1F*) hSpectraMC[itune][kKaon][kPos]->Clone();
      hKPiRatioMC[itune]->Add(hSpectraMC[itune][kKaon][kNeg]);
      TH1F * htmpLoc = (TH1F*) hSpectraMC[itune][kPion][kPos]->Clone();
      htmpLoc->Add(hSpectraMC[itune][kPion][kNeg]);
      hKPiRatioMC[itune]->Divide(htmpLoc);
      hKPiRatioMC[itune]->SetYTitle("(K^{+}+K^{-})/(#pi^{+}+#pi^{-})");    
    }
  }
  


  TH1F * hPPiRatio = (TH1F*) hSpectra[iCombInStudy][kProton][kPos]->Clone();
  hPPiRatio->Add(hSpectra[iCombInStudy][kProton][kNeg]);
  htmp = (TH1F*) hSpectra[iCombInStudy][kPion][kPos]->Clone();
  htmp->Add(hSpectra[iCombInStudy][kPion][kNeg]);
  hPPiRatio->Divide(htmp);
  hPPiRatio->SetYTitle("(p+#bar{p})/(#pi^{+}+#pi^{-})");

  TH1F * hPPiRatioMC[kNTunes];
  if(showMC){
    for(Int_t itune = 0; itune < kNTunes; itune++){
      hPPiRatioMC[itune] = (TH1F*) hSpectraMC[itune][kProton][kPos]->Clone();
      hPPiRatioMC[itune]->Add(hSpectraMC[itune][kProton][kNeg]);
      TH1F * htmpLoc = (TH1F*) hSpectraMC[itune][kPion][kPos]->Clone();
      htmpLoc->Add(hSpectraMC[itune][kPion][kNeg]);
      hPPiRatioMC[itune]->Divide(htmpLoc);
      hPPiRatioMC[itune]->SetYTitle("(p+#bar{p})/(#pi^{+}+#pi^{-})");    
    }
  }
 
  // Phenix ratio
  TFile * fPhenix = new TFile("./Files/phenixdata.root");
  TH1F * hPPiRatioPhenix = (TH1F*) fPhenix->Get("ProtonPHNX");
  hPPiRatioPhenix->Add((TH1F*) fPhenix->Get("PbarPHNX"));
  htmp = (TH1F*) fPhenix->Get("PiPlusPHNX");
  htmp->Add((TH1F*) fPhenix->Get("PiMinusPHNX"));
  hPPiRatioPhenix = AliPWGHistoTools::DivideHistosDifferentBins(hPPiRatioPhenix,htmp);
  hPPiRatioPhenix->SetMarkerStyle(24);
  hPPiRatioPhenix->SetBinContent(14,0);
  hPPiRatioPhenix->SetBinError(14,0);
  // Draw
//   TH2F * hempty = new TH2F(TString("hempty"),"hempty",100,0.,1.5, 100, 0.001,1.8);
//   hempty->SetXTitle("p_{t} (GeV/c)");
  //  hempty->SetYTitle("");

  // tmp: overlay levi fits
  AliPWGFunc * fm2 = new AliPWGFunc;
  fm2->SetVarType(AliPWGFunc::kdNdpt);
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
    TString detName = detFlag[iCombInStudy];
    detName.ReplaceAll(" ", "_");
    detName.ReplaceAll("+", "");

    TCanvas * c1 = new TCanvas(TString("cRatio_")+detName+TString("_")+partFlag[ipart], TString("cRatio_")+detName+partFlag[ipart]);
    c1->SetGridy();
    hPosNegRatio[ipart]->Draw();
    TF1 * fRatio = new TF1 (TString("fRatio")+partFlag[ipart], TString(fLevi[ipart][kPos]->GetName())+"/"+fLevi[ipart][kNeg]->GetName());
    //    fRatio->Draw("same");
    fRatio->SetRange(0,5);
    PrintCanvas(c1,printFormats);
    
  }


  TCanvas * c2 = new TCanvas(TString("cRatio_KPi"),TString("cRatio_KPi"));  
  //  c2->SetGridy();

  hKPiRatioSyst->SetMaximum(0.62);  
  hKPiRatioSyst->SetFillColor(kGray);
  hKPiRatioSyst->SetMarkerStyle(24);
  //  hKPiRatioSyst->Draw("e2");
  // TMP FIXME: ERROR FROM FP:
  TFile * f = TFile::Open("./Files/ITSsa-systematics-20101014.root");
  TH1F * hsystRatio = (TH1F*) f->Get("hSystRatioKpi");
  TH1F * hsyst = new TH1F(*htemplate);
  AliPWGHistoTools::GetValueAndError(hsyst,hKPiRatio,hsystRatio,1);
  //  AliPWGHistoTools::GetHistoCombinedErrors(hsyst,hKPiRatio);
  //  hsyst->Draw("samee1");
  // END TMP
  hKPiRatio->SetMaximum(0.62);  
  hKPiRatio->Draw("");
  TLegend * lMC = new TLegend(0.526846, 0.18007, 0.887584,0.407343);
  lMC->SetFillColor(kWhite);
  TGraphErrors ** gStar = StarPPSpectra();
  cout << gStar[6] << endl;
  
  gStar[6]->SetMarkerStyle(kOpenStar);
  gStar[6]->Draw("P");
  
  if(showE735){
    gROOT->LoadMacro("GetE735Ratios.C");
    // GetE735Ratios(0,0)->Draw("EX0,same");
    // GetE735Ratios(0,1)->Draw("EX0,same");
    // GetE735Ratios(0,2)->Draw("EX0,same");
    GetE735Ratios(0,3)->Draw("EX0,same");
  }
  hKPiRatio->SetMarkerStyle(20);
  hKPiRatio->SetMarkerColor(kRed);
  hKPiRatio->SetLineColor(kRed);
  //hKPiRatio->Draw("same");
  
  if(showMC){
    for(Int_t itune = 0; itune < kNTunes; itune++){
      lMC->AddEntry(hKPiRatioMC[itune],mcTuneName[itune]);
      hKPiRatioMC[itune]->SetLineWidth(2);    
      hKPiRatioMC[itune]->Draw("same,chist");    	
    }
  
    lMC->Draw();
  }
  cout << "12" << endl;

  if(showE735){
    TLegend * l = new TLegend(  0.1879,  0.68,  0.54,0.92);
    l->SetFillColor(kWhite);
    l->AddEntry(hKPiRatio, "ALICE, #sqrt{s} = 900 GeV","lpf");
    l->AddEntry(gStar[6],  "STAR, #sqrt{s} = 200 GeV", "lp");
    // l->AddEntry(GetE735Ratios(0,0), "E735, #sqrt{s} = 300 GeV","lpf");
    // l->AddEntry(GetE735Ratios(0,1), "E735, #sqrt{s} = 540 GeV","lpf");
    // l->AddEntry(GetE735Ratios(0,2), "E735, #sqrt{s} = 1000 GeV","lpf");
    l->AddEntry(GetE735Ratios(0,3), "E735, #sqrt{s} = 1800 GeV","lpf");
    l->Draw();
  }


  TCanvas * c3 = new TCanvas(TString("cRatio_PPi"),TString("cRatio_PPi"));  
  cout << "DRAWING" << endl;
  //  c3->SetGridy();
  hPPiRatio->SetMarkerStyle(20);
  hPPiRatio->SetMarkerColor(kRed);
  hPPiRatio->SetLineColor(kRed);
  hPPiRatio->Draw();
  hPPiRatio->SetMaximum(0.39);
  if(showMC){
    lMC = new TLegend(0.526846, 0.18007, 0.887584,0.407343);
    lMC->SetFillColor(kWhite);

    for(Int_t itune = 0; itune < kNTunes; itune++){
      lMC->AddEntry(hKPiRatioMC[itune],mcTuneName[itune]);
      hPPiRatioMC[itune]->SetLineWidth(2);    
      hPPiRatioMC[itune]->Draw("same,chist");    	
    }
  
    lMC->Draw();
  }
  hPPiRatioPhenix->Draw("same");

  
  TLegend * l = new TLegend(0.186242, 0.781469, 0.538591, 0.921329);
  l->SetFillColor(kWhite);
  l->AddEntry(hPPiRatio, "ALICE, #sqrt{s} = 900 GeV","lpf");
  l->AddEntry(hPPiRatioPhenix, "PHENIX, #sqrt{s} = 200 GeV","lpf");
  //  TLegend * l = new TLegend(0.206376, 0.77972, 0.600671, 0.909091);
  l->Draw();

  PrintCanvas(c2,printFormats);
  PrintCanvas(c3,printFormats);

}

void DrawWithSyst() {

  // Draws detector and combined with syst errors. 

  for(Int_t idet = 0; idet < kNHist; idet++){
    if(idet > kITSTPC && idet != iCombInStudy) continue;
    for(Int_t icharge = 0; icharge < kNCharge; icharge++){
      TCanvas * c = new TCanvas(TString("cWithSyst")+detFlag[idet]+chargeFlag[icharge],TString("cWithSyst")+detFlag[idet]+chargeFlag[icharge]);
      TH2F * hempty = new TH2F(TString("hempty")+long(icharge),"hempty",100,0.,2.9, 100, 0.0005,5);
      hempty->SetXTitle("p_{t} (GeV/c)");
      hempty->SetYTitle("1/N_{ev} d^{2}N / dydp_{t} (GeV/c)^{-1}");
      hempty->GetYaxis()->SetTitleOffset(1.35);
      hempty->GetXaxis()->SetTitleOffset(1.1);
      hempty->Draw();      
      c->SetLogy();
      for(Int_t ipart = 0; ipart < kNPart; ipart++){
	//	cout << detFlag[idet] << " " << chargeFlag[icharge] << " " << partFlag[ipart] << endl;
	
	TString opt = ipart ? "" : "same";
	TH1F * hsyst = new TH1F(*htemplate);
	AliPWGHistoTools::GetValueAndError(hsyst,hSpectra[idet][ipart][icharge],hSystError[idet][ipart][icharge],kTRUE);
	hsyst->SetFillColor(kGray);
	hsyst->Draw("e5,same");
	hSpectra[idet][ipart][icharge]->Draw("same");
	// Do not draw syst error outside of spectra range
	Float_t lowEdge = AliPWGHistoTools::GetLowestNotEmptyBinEdge(hSpectra[idet][ipart][icharge]);
	Float_t hiEdge  = AliPWGHistoTools::GetHighestNotEmptyBinEdge(hSpectra[idet][ipart][icharge]);
	Int_t lowBin = hSystError[idet][ipart][icharge]->FindBin(lowEdge);
	Int_t hiBin = hSystError[idet][ipart][icharge]->FindBin(hiEdge);
	Int_t nbin = 	hSystError[idet][ipart][icharge]->GetNbinsX();
	for(Int_t ibin = 1; ibin < lowBin; ibin++){
	  hSystError[idet][ipart][icharge]->SetBinContent(ibin,0);	  
	  hSystError[idet][ipart][icharge]->SetBinError(ibin,0);	  
	}
	for(Int_t ibin = hiBin; ibin <= nbin; ibin++){
	  hSystError[idet][ipart][icharge]->SetBinContent(ibin,0);	  
	  hSystError[idet][ipart][icharge]->SetBinError(ibin,0);	  
	}
	// Draw relative syst error, convert to graph to avoid vertical bars at the edges (empty bins in histo)
	AliPWGHistoTools::GetGraphFromHisto(hSystError[idet][ipart][icharge])->Draw("LX");

				     

      }
      PrintCanvas(c,printFormats);
    }
  }
}

void Help() {

  cout << "Macro: void CombineSpectra(Int_t analysisType=kDoFits, Int_t  fitFuncID = kFitLevi) " << endl;
  cout << "" << endl;
  cout << "WARNING: some function modify the histos upon execution. You are advised to quit root" << endl;
  cout << "         before running a different function" << endl;
  cout << "" << endl;
  cout << "Possible Arguments" << endl;
  cout << "- analysisType:" << endl;  
  cout << "    kDoFits:           Fit Combined Spectra " << endl;
  cout << "    kDoRatios:         Particle ratios K/pi and p/pi" << endl;
  cout << "    kDoSuperposition:  Compare different detectors (superimpose and ratios)" << endl;
  cout << "    kDoCompareStar:    Compare combined spectra to star results" << endl;
  cout << "    kDoDrawWithModels: Compare combined spectra and models" << endl;
  cout << "    kDoDrawSyst:       Draws spectra from individual detectors with their systematic error" << endl;
  cout << "    kDoCompareToAllCharged:       convert to dn/deta/dpt, sum up and compare to all charged" << endl;
  cout << "    kDoHelp:           This help" << endl;
  cout << "- fitFuncID, function used to extrapolate and compute yields" << endl;
  cout << "    An analitic fit function [kFitLevi, kFitUA1, kFitPowerLaw]" << endl;
  cout << "    Or a shape from a MC moder [kFitPhojet, kFitAtlasCSC, kFitCMS6D6T, kFitPerugia0]" << endl;
  cout << "    Which is fitted to the data at low pt and used to extrapolate at low pt" << endl;


}

void PrintCanvas(TCanvas* c,const TString formats) {
  // print a canvas in every of the given comma-separated formats
  // ensure the canvas is updated
  if(!doPrint) return;
  c->Update();
  gSystem->ProcessEvents();
  TObjArray * arr = formats.Tokenize(",");
  TIterator * iter = arr->MakeIterator();
  TObjString * element = 0;
  TString name  =c ->GetName();
  name.ReplaceAll(" ","_");
  name.ReplaceAll("+","Plus");
  name.ReplaceAll("-","");
  while ((element = (TObjString*) iter->Next())) {
    c->Print(name+ "."+element->GetString());
  }
}


void LoadLibs(){

  gSystem->Load("libTree");
  gSystem->Load("libVMC");
  gSystem->Load("libMinuit");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libCORRFW");
  gSystem->Load("libPWG2spectra");

}

void RebinHistosForRatios() {
  // rebins all spectra histograms 

  cout << "WARNING: REBINNING HISTOS BEFORE RATIOS" << endl;
  

  const Float_t rebinBins[] = {0.05,0.1,0.12,0.14,0.16,0.18,0.20,0.25,0.30,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,1.1,1.2,1.3,1.4,1.6,1.8,2.4};
  Int_t nrebinBins=29;



  Int_t idet = iCombInStudy;
  for(Int_t ipart = 0; ipart < kNPart; ipart++){
    for(Int_t icharge = 0; icharge < kNCharge; icharge++){  

      TH1F * h =  hSpectra[idet][ipart][icharge];
      TString hname = h->GetName();
      hname += "_rebin";
     
      TH1F *hrebin   = new TH1F(hname, hname,nrebinBins,rebinBins);
      hrebin->Sumw2();

      for(Int_t i=1;i <= h->GetNbinsX();i++){
	Float_t x = h->GetBinCenter(i);
	Int_t j=0;
	while(x > rebinBins[j+1]) j++;
	j++;
	hrebin->AddBinContent(j,h->GetBinContent(i));

	Float_t err = hrebin->GetBinError(j) + h->GetBinError(i)*h->GetBinError(i);
	hrebin->SetBinError(j,err);
      }

      for(Int_t i=1;i <= hrebin->GetNbinsX();i++){
	hrebin->SetBinError(i,TMath::Sqrt(hrebin->GetBinError(i)));
      }
      hSpectra[idet][ipart][icharge] = hrebin;
    }
  }

}
