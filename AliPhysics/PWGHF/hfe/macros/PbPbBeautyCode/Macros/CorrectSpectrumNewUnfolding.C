// This Script takes the output from the fit calculation,
// applies efficiency correction and unfolds, creates
// pT-spectrum and uses a pp input spectrum for the RAA

#include <MartinsStyle/mystyle.h>
#include <comparison.h>
#include <AliHFEUnfolding.h>
#include <TMath.h>

int gAliFirstBin;
int gAliLastBin;
int gAliNBins;
int gAliNPara;

TH2D * gAliR;
TH1D * gAliUnfoldingSpectrum;



void SetStyle(Bool_t graypalette=kFALSE);
void myPadSetUp(TPad *currentPad, float currentLeft=0.11, float currentTop=0.04, float currentRight=0.04, float currentBottom=0.15);
void DrawLogo (Int_t logo=0, Double_t xmin =  0.28, Double_t ymin= 0.68) ;
void FakeHistosOnlyForExample(TH1*&hstat, TH1*&hsyst, TH1*&hsystCorr);
void LoadLibs();


// Preferred colors and markers
const Int_t fillColors[] = {kGray+1,  kRed-10, kBlue-9, kGreen-8, kMagenta-9, kOrange-9,kCyan-8,kYellow-7}; // for syst bands
const Int_t colors[]     = {kBlack, kRed+1 , kBlue+1, kGreen+3, kMagenta+1, kOrange-1,kCyan+2,kYellow+2};
const Int_t markers[]    = {kFullCircle, kFullSquare,kOpenCircle,kOpenSquare,kOpenDiamond,kOpenCross,kFullCross,kFullDiamond,kFullStar,kOpenStar};


void figTemplate() {
  // Load necessary libraries
  LoadLibs();
  // Set the default style
  SetStyle();

  // Prepare Figure, please stick to the default canvas size(s) unless absolutely necessary in your case
  // Rectangular
  TCanvas *cfig = new TCanvas("cfig", "Alice Figure Template", 800, 600); 
  // Square
  //TCanvas *cfig = new TCanvas("cfig", "Alice Figure Template", 800, 800); 
  cfig->SetLogy();
  // Set Titles etc..
  TH1 * h = cfig->DrawFrame(0,10,10,1000);
  // === Commonly used x/ titles: ===
  // pt invariant yields
  const char *  texPtX="#it{p}_{T} (GeV/#it{c})";
  const char *  texPtY="1/#it{N}_{ev} 1/(2#pi#it{p}_{T}) d#it{N}/(d#it{p}_{T}d#it{y}) ((GeV/#it{c})^{-2})";
  // mt invariant yields
  const char *  texMtX="#it{m}_{T} (GeV/#it{c^2})";
  const char *  texMtY="1/#it{N}_{ev} 1/(2#pi#it{m}_{T}) d#it{N}/(d#it{m}_{T}d#it{y}) ((GeV/#it{c^2})^{-2}) "; 
  // Invariant mass with decay products K and pi
  const char *  texMassX="#it{M}_{K#pi} (GeV/#it{c}^2)";
  const char *  texMassY="d#it{N}/(d#it{M}_{K#pi}";
  // <pt>, npart
  const char * texMeanPt    = "#LT#it{p}_{T}#GT (GeV/#it{c})";
  const char * texMeanNpart = "#LT#it{N}_{part}#GT";


  // Set titles
  h->SetXTitle(texPtX);
  // Please be consistent on the y label
  h->SetYTitle(texPtY);

  // Draw your histos here:
  TH1* hstat, *hsyst, *hsystCorr;
  Int_t icolor=0;
  FakeHistosOnlyForExample(hstat,hsyst,hsystCorr);
  hsystCorr->SetFillColor(fillColors[icolor]);
  hsyst    ->SetFillColor(fillColors[icolor]);
  hsystCorr->Draw("E3,same");
  hsyst->SetFillStyle(0); // To draw empty boxes
  hsyst->SetLineColor(colors[icolor]); // To draw empty boxes
  hsystCorr->SetLineColor(colors[icolor]); // To draw empty boxes
  hsyst->Draw("E2,same");
  hstat->Draw("E,same");
  hstat->SetMarkerStyle(markers[0]);
  // use the same color for markers and lines
  hstat->SetMarkerColor(colors [icolor]);
  hstat->SetLineColor  (colors [icolor]);

  // Draw the logo   
  //  0: Just "ALICE" (for final data), to be added only if ALICE does not appear otherwise (e.g. in the legend)
  //  >0: ALICE Preliminary
  DrawLogo(1, 0.59, 0.81);

  // You should always specify the colliding system
  // NOTATION: pp, p-Pb, Pb-Pb. 
  // Don't forget to use #sqrt{s_{NN}} for p-Pb and Pb-Pb
  // You can change the position of this with
  TLatex * text = new TLatex (0.55,87,"p-Pb #sqrt{#it{s}}_{NN} = 5.02 TeV");
  text->Draw();
  TLatex * text2 = new TLatex (0.55,55,"V0A Multiplicity Classes (Pb-Side)");
  text2->SetTextSizePixels(24);
  text2->Draw();
 
  //Legend, if needed
  TLegend * leg = new TLegend(  0.19,  0.19,  0.57, 0.42);
  leg->AddEntry(hstat,     "0-5\%, stat errors",   "LPE");
  leg->AddEntry(hsyst,     "syst error (Uncorrelated)",  "F");
  leg->AddEntry(hsystCorr, "syst error (Correlated)",    "F" );
  leg->SetFillColor(0);
  leg->SetTextSize(gStyle->GetTextSize()*0.8);
  leg->Draw();
  // Save to HEP data

  AliHEPDataParser * hepParser = new AliHEPDataParser(hstat, hsyst);
  hepParser->SetTitle("pt distribution of pi+-, arXiv:XXXX.YYYY");
  hepParser->SetName("1/Nev 1/p_T 1/2pi d^2N/(dp_Tdy) (GeV/c)^{-1}"); 
  hepParser->SetXaxisName("PT IN GEV/c");
  hepParser->SetReaction("RE: P PB --> PI + X");
  hepParser->SetEnergy("SQRT(SNN) : 5020.0 GeV");
  hepParser->SetRapidityRange("YRAP : -0.5 - +0.5");
  hepParser->SaveHEPDataFile("figTemplateHEPData.txt");    // it must be specified explicity if graphs are to be used






}

//________________________________
void LoadLibs() {
  gSystem->Load("libCore.so");  
  gSystem->Load("libGeom.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libVMC");
  gSystem->Load("libTree");
  gSystem->Load("libMinuit");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libCORRFW");
  gSystem->Load("libPWGTools");
}

void DrawLogo (Int_t logo, Double_t xmin, Double_t ymin) {

  // Logo is not needed anymore, now we only write alice preliminary
  // Logo:
  // 0: Justr writes "ALICE" (for final data)
  // Anything eles: writes "ALICE Preliminary"

  TLatex *   tex = new TLatex(xmin,ymin, logo ? "ALICE Preliminary (requested)" : "ALICE");
  tex->SetNDC();
  tex->SetTextFont(42);
  tex->Draw();

  // OLD logo
  //  TPad * currentPad = gPad;
  // Double_t AliLogo_LowX =xmin;
  // Double_t AliLogo_LowY = ymin;
  // Double_t AliLogo_Height = size;
  // //ALICE logo is a  file that is 821x798 pixels->should be wider than a square
  // Double_t AliLogo_Width  = (821./798.) * AliLogo_Height * gPad->GetWh() / gPad->GetWw();
  
  // TPad *myPadLogo = new TPad("myPadLogo", "Pad for ALICE Logo",AliLogo_LowX,AliLogo_LowY,AliLogo_LowX+AliLogo_Width,AliLogo_LowY+AliLogo_Height);
  // myPadSetUp(myPadLogo,0,0,0,0);
  // //myPadLogo->SetFixedAspectRatio(1);
  // myPadLogo->Draw();
  // myPadLogo->cd();
  // if (logo == 0) {
  //   myPadLogo->SetFillColor(2); // color to first figure out where is the pad then comment !
  // } else if (logo == 1){
  //   TASImage *myAliceLogo = new TASImage(performanceLogoPath);
  //   myAliceLogo->Draw();
  // } else if (logo == 2) {
  //   TASImage *myAliceLogo = new TASImage(preliminaryLogoPath);
  //   myAliceLogo->Draw();
  // }
  // // go back to the old pad
  // currentPad->cd();

}

void myPadSetUp(TPad *currentPad, float currentLeft, float currentTop, float currentRight, float currentBottom){
  currentPad->SetLeftMargin(currentLeft);
  currentPad->SetTopMargin(currentTop);
  currentPad->SetRightMargin(currentRight);
  currentPad->SetBottomMargin(currentBottom);
  return;
}


void SetStyle(Bool_t graypalette) {
  cout << "Setting style!" << endl;
  
  gStyle->Reset("Plain");
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  if(graypalette) gStyle->SetPalette(8,0);
  else gStyle->SetPalette(1);
  gStyle->SetCanvasColor(10);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetFrameLineWidth(1);
  gStyle->SetFrameFillColor(kWhite);
  gStyle->SetPadColor(10);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetHistLineWidth(1);
  gStyle->SetHistLineColor(kRed);
  gStyle->SetFuncWidth(2);
  gStyle->SetFuncColor(kGreen);
  gStyle->SetLineWidth(2);
  gStyle->SetLabelSize(0.045,"xyz");
  gStyle->SetLabelOffset(0.01,"y");
  gStyle->SetLabelOffset(0.01,"x");
  gStyle->SetLabelColor(kBlack,"xyz");
  gStyle->SetTitleSize(0.05,"xyz");
  gStyle->SetTitleOffset(1.25,"y");
  gStyle->SetTitleOffset(1.2,"x");
  gStyle->SetTitleFillColor(kWhite);
  gStyle->SetTextSizePixels(26);
  gStyle->SetTextFont(42);
  //  gStyle->SetTickLength(0.04,"X");  gStyle->SetTickLength(0.04,"Y"); 

  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(kWhite);
  //  gStyle->SetFillColor(kWhite);
  gStyle->SetLegendFont(42);


}






void UnfoldingLikelihood(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{ // one bin added at beginning and end in gAliUnfoldingSpectrum
  //npar should be lastbin-firstbin+3
  double estimatedBinCounts, difference, error;
  f=0;
  for(int i=gAliFirstBin;i<=gAliLastBin;i++) // only measured contribute to Chi2
  {
    estimatedBinCounts=0;
    for(int j=gAliFirstBin-1;j<=gAliLastBin+1;j++)
    {
      estimatedBinCounts+=R->GetBinContent(j,i)*par[j-gAliFirstBin+1];  // par: "true" distribution!
    }
    difference=estimatedBinCounts-gAliUnfoldingSpectrum->GetBinContent(i);
    //cout << "bin: " << i << " estimated: " << estimatedBinCounts << " measured: " << gAliUnfoldingSpectrum->GetBinContent(i) << endl;
    error = gAliUnfoldingSpectrum->GetBinError(i);
    f+= (difference * difference)/(error*error);   // Chi2 calculation
  }
  f+= (par[1]*2.-par[0])*(par[1]*2.-par[0])/(par[1]*par[1]);
  f+= (par[gAliNPara-2]/2.-par[gAliNPara-1])*(par[gAliNPara-2]/2.-par[gAliNPara-1])/(par[gAliNPara-2]*par[gAliNPara-2]/4.);
}

void NormalizeResponseMatrix(void)
{
  double sum;
  for(int i=1;i<=gAliNBins;i++)
  {
    sum=0;
    for(int j=1;j<=gAliNBins;j++)
      sum +=gAliR->GetBinContent(i,j);
    if(sum>0)
    for(int j=1;j<=gAliNBins;j++)
     gAliR->SetBinContent(i,j,gAliR->GetBinContent(i,j)/sum);
  }
}

void CorrectSpectrumNewUnfolding(void)
{
  mystyle();
  TString * OutputFileName = new TString("BSpecCorr.root");
  
  TFile * inFileSpectrum = new TFile("FitOutput/Mar19/SpectraFitPbPb-05HalfRAAlastHadron.root");
  TFile * inFileCorrection = new TFile("correctionPbPb/correctionDiagramsPbPbenh.root");
  TH2D * Cuts = (TH2D*) inFileCorrection->Get("individualCuts");
  TH1D * efficiency = (TH1D*)Cuts->ProjectionX("eff", 8,8); // 8,8 with TOF!
  efficiency->Divide(Cuts->ProjectionX("tot",1,1));
  efficiency->Scale(0.73);
 gAliR = (TH2D*) inFileCorrection->Get("pTMCdatabins");
  
  double nEvents = ((TH1D*)inFileSpectrum->Get("evtsb0"))->Integral();
  cout << "Number of Events: " << nEvents << endl;
  TH1D * spectrum = inFileSpectrum->Get("BSpectrum");
  TH1D * spectrumsys = inFileSpectrum->Get("BSpectrumSys");
  gAliUnfoldingSpectrum = (TH1D*) spectrum->Clone(); gAliUnfoldingSpectrum->SetName("gAliUnfoldingSpectrum");
  gAliNBins = spectrum->GetNbinsX();
  gAliFirstBin=1;
  gAliLastBin=gAliNBins;
  
  

  
  
  
  // Here  do the unfolding calculation!
  int gAliNParameters = gAliLastBin- gAliFirstBin + 1;
  while(spectrum->GetBinContent(gAliFirstBin)<10.1) gAliFirstBin++;
  while(spectrum->GetBinContent(gAliLastBin)<10.1) gAliLastBin--;
  gAliNPara=gAliLastBin-gAliFirstBin+3;
  // Sum Errors for Chi2 method!
  //for(int i=gAliFirstBin;i<=gAliLastBin;i++)
    //gAliUnfoldingSpectrum->SetBinError(i, TMath::Sqrt(gAliUnfoldingSpectrum->GetBinError(i)*gAliUnfoldingSpectrum->GetBinError(i)+ spectrumsys->GetBinError(i)*spectrumsys->GetBinError(i)));
  // end of error summation
  cout << "First bin: " << gAliFirstBin << endl;
  cout << "Last bin: " << gAliLastBin << endl;
  
  // Diagrams for systematic uncertainties
  TF1 * FunctionMotherpT = new TF1("FunctionMotherpT", "(.12)", 1., 6.);
  TF1 * FunctionFakeConversions = new TF1("FunctionFakeConversions", "(.1/x)", 1., 6.);
  TF1 * FunctionContamination = new TF1("FunctionContamination", "( TMath::Exp(-x+1)*0.3)", 1., 6.);
  TF1 * FunctionResolution = new TF1("FunctionResolution", ".06", 1., 6.);
  TF1 * FunctionStrangeness = new TF1("FunctionStrangeness", "(.05/x)", 1., 6.);
  double sysPID = 0.12;
  //double sysHadrCont = 0.025;
  double sysEnergyScaling = 0.1;
  double additionalErrorInFirstBin = 0.03;
  double additionalErrorInLastBin = 0.08;
  //TF1 * shiftError = new TF1("shiftError", "0.4/x", 0.5, 10.0);
  TF1 * shiftError = new TF1("shiftError", "0.13", 0.5, 10.0);
  
  TH1D * UncertPID = (TH1D*) spectrum->Clone(); UncertPID->SetName("UncertPID"); UncertPID->Scale(0.);
  setOPT_hists((TH1F*)UncertPID, "p_{T} (GeV/c)", "sys. uncertainty", 510, 20, 1.3, kGray+2);
  UncertPID->SetLineStyle(2);
  TH1D * UncertErrorBSpec = (TH1D*) spectrum->Clone(); UncertErrorBSpec->SetName("UncertErrorBSpec"); UncertErrorBSpec->Scale(0.);
  setOPT_hists((TH1F*)UncertErrorBSpec, "p_{T} (GeV/c)", "sys. uncertainty", 510, 20, 1.3, colors[7]);
  UncertErrorBSpec->SetLineStyle(4);
  UncertErrorBSpec->SetLineWidth(3);
  TH1D * UncertHadrCont = (TH1D*) spectrum->Clone(); UncertHadrCont->SetName("UncertHadrCont"); UncertHadrCont->Scale(0.);
  setOPT_hists((TH1F*)UncertHadrCont, "p_{T} (GeV/c)", "sys. uncertainty", 510, 20, 1.3, colors[6]);
  UncertHadrCont->SetLineStyle(5);
  TH1D * UncertIPRes = (TH1D*) spectrum->Clone(); UncertIPRes->SetName("UncertIPRes"); UncertIPRes->Scale(0.);
  setOPT_hists((TH1F*)UncertIPRes, "p_{T} (GeV/c)", "sys. uncertainty", 510, 20, 1.3, colors[7]);
  UncertIPRes->SetLineStyle(1);
  TH1D * UncertStrangeness = (TH1D*) spectrum->Clone(); UncertStrangeness->SetName("UncertStrangeness"); UncertStrangeness->Scale(0.);
  setOPT_hists((TH1F*)UncertStrangeness, "p_{T} (GeV/c)", "sys. uncertainty", 510, 20, 1.3, colors[8]);
  UncertStrangeness->SetLineStyle(1);
  TH1D * UncertEnergyScaling = (TH1D*) spectrum->Clone(); UncertEnergyScaling->SetName("UncertEnergyScaling"); UncertEnergyScaling->Scale(0.);
  setOPT_hists((TH1F*)UncertEnergyScaling, "p_{T} (GeV/c)", "sys. uncertainty", 510, 20, 1.3, colors[5]);
  UncertEnergyScaling->SetLineStyle(6);
  UncertEnergyScaling->SetLineWidth(3);
  TH1D * UncertUnfolding = (TH1D*) spectrum->Clone(); UncertUnfolding->SetName("UncertUnfolding"); UncertUnfolding->Scale(0.);
  setOPT_hists((TH1F*)UncertUnfolding, "p_{T} (GeV/c)", "sys. uncertainty", 510, 20, 1.3, colors[4]);
  UncertUnfolding->SetBinContent(gAliFirstBin, additionalErrorInFirstBin); UncertUnfolding->SetBinContent(gAliLastBin, additionalErrorInLastBin);
  UncertUnfolding->SetLineStyle(7);
  TH1D * UncertV0Shift = (TH1D*) spectrum->Clone(); UncertV0Shift->SetName("UncertV0Shift"); UncertV0Shift->Scale(0.);
  setOPT_hists((TH1F*)UncertV0Shift, "p_{T} (GeV/c)", "sys. uncertainty", 510, 20, 1.3, colors[3]);
  UncertV0Shift->SetLineWidth(3);
  UncertV0Shift->SetLineStyle(3);
  TH1D * UncertFit = (TH1D*) spectrum->Clone(); UncertFit->SetName("UncertFit"); UncertFit->Scale(0.);
  setOPT_hists((TH1F*)UncertFit, "p_{T} (GeV/c)", "sys. uncertainty", 510, 20, 1.3, colors[2]);
  TH1D * Uncertpp = (TH1D*) spectrum->Clone(); Uncertpp->SetName("Uncertpp"); Uncertpp->Scale(0.);
  setOPT_hists((TH1F*)Uncertpp, "p_{T} (GeV/c)", "sys. uncertainty", 510, 20, 1.3, colors[1]);
  Uncertpp->SetLineStyle(8);
  TH1D * UncertSum = (TH1D*) spectrum->Clone(); UncertSum->SetName("UncertSum"); UncertSum->Scale(0.);
  setOPT_hists((TH1F*)UncertSum, "p_{T} (GeV/c)", "sys. uncertainty", 510, 20, 1.3, colors[0]);
  Uncertpp->SetLineStyle(9);
  
  
  double ctr=0;
  
  for(int i=gAliFirstBin;i<=gAliLastBin;i++)
  {
    ctr=UncertV0Shift->GetBinCenter(i);
    UncertPID->SetBinContent(i, sysPID);
    UncertIPRes->SetBinContent(i, FunctionResolution->Eval(ctr));
    //UncertErrorBSpec->SetBinContent(i, ErrorBSpec);
    UncertErrorBSpec->SetBinContent(i, FunctionMotherpT->Eval(ctr));
    UncertStrangeness->SetBinContent(i, FunctionStrangeness->Eval(ctr));
    //UncertHadrCont->SetBinContent(i, sysHadrCont);
    UncertHadrCont->SetBinContent(i, FunctionContamination->Eval(ctr));
    UncertEnergyScaling->SetBinContent(i, sysEnergyScaling);
    //UncertV0Shift->SetBinContent(i, shiftError->Eval(ctr));
    UncertV0Shift->SetBinContent(i, FunctionFakeConversions->Eval(ctr));
    UncertFit->SetBinContent(i, spectrumsys->GetBinError(i)/spectrumsys->GetBinContent(i));
  }
  // end of diagram definitions
  
  
  if(!true){  // This is for comparison of unfolding!
  for(int i=1;i<gAliFirstBin-1;i++)  for(int j=1; j<=19;j++)  {R->SetBinContent(j,i,0.);}
  for(int i=gAliLastBin+2;i<=19-1;i++)  for(int j=1; j<=19;j++)  {R->SetBinContent(j,i,0.);}}
  NormalizeResponseMatrix();
  gAliUnfoldingSpectrum->SetBinContent(gAliFirstBin-1,gAliUnfoldingSpectrum->GetBinContent(gAliFirstBin)*2); // This sets non-measured bins to reasonable values
  gAliUnfoldingSpectrum->SetBinError(gAliFirstBin-1,gAliUnfoldingSpectrum->GetBinError(gAliFirstBin)*2);
  gAliUnfoldingSpectrum->SetBinContent(gAliLastBin+1,gAliUnfoldingSpectrum->GetBinContent(gAliLastBin)/2);
  gAliUnfoldingSpectrum->SetBinError(gAliLastBin+1,gAliUnfoldingSpectrum->GetBinError(gAliLastBin)/2);
  
  // Prepare Unfolding
  AliHFEUnfolding * UnfoldingObject = new AliHFEUnfolding(); 
  UnfoldingObject->SetData(gAliUnfoldingSpectrum, gAliFirstBin-1,gAliLastBin+1);   //(-1, +1)
  UnfoldingObject->SetResponseMatrix(gAliR);
  //TFile * CheckOut = new TFile("checkHists.root", "RECREATE");
  //gAliUnfoldingSpectrum->Write();
  //R->Write();
  //CheckOut->Close();
  // Unfold
  TH1D * CorrectedSpectrum = UnfoldingObject->Unfold();
  
  // Here comes the check with TUnfold!
  //Make Smaller Matrix
  cout << "first bin: " << gAliFirstBin << endl;
  TH2D * SmallResponse = new TH2D("SmallResponse", "SmallResponse", gAliLastBin-gAliFirstBin+3, 0., 1., gAliLastBin-gAliFirstBin+3, 0., 1.);
  TH1D * SmallData = new TH1D("SmallData", "SmallData", gAliLastBin-gAliFirstBin+3, 0., 1.);
  for(int i = 1; i<= gAliLastBin-gAliFirstBin+3;i++)
  {SmallData->SetBinContent(i, gAliUnfoldingSpectrum->GetBinContent(i+gAliFirstBin-2)); SmallData->SetBinError(i, gAliUnfoldingSpectrum->GetBinError(i+gAliFirstBin-2));}
  for(int i = 1; i<= gAliLastBin-gAliFirstBin+3;i++)
    for(int j = 1; j<= gAliLastBin-gAliFirstBin+3;j++)
  {SmallResponse->SetBinContent(i, j,gAliR->GetBinContent(i+gAliFirstBin-2, j+gAliFirstBin-2));}
  //TH1D * iniMC = SmallResponse->ProjectionX(); iniMC->SetName("iniMC");
  TH1D * UnfOut = SmallResponse->ProjectionX(); UnfOut->SetName("UnfOut");
  
  TUnfold TUnfoldingObject(SmallResponse, TUnfold::kHistMapOutputHoriz, TUnfold::kRegModeNone, TUnfold::kEConstraintNone);
  TUnfoldingObject.SetInput(SmallData);
  TUnfoldingObject.DoUnfold(0.);
  TUnfoldingObject.GetOutput(UnfOut);
  TH1D * UnfComp = CorrectedSpectrum->Clone(); UnfComp->SetName("UnfComp");
  TH1D * UnfErrComp = CorrectedSpectrum->Clone(); UnfErrComp->SetName("UnfErrComp"); UnfErrComp->SetLineColor(kGreen);
  //Divide Matrix unfolding by Chi2
  for(int i=gAliFirstBin;i<=gAliLastBin+1;i++){
    UnfErrComp->SetBinContent(i, UnfComp->GetBinError(i)/ UnfOut->GetBinError(i-gAliFirstBin+2));
    UnfErrComp->SetBinError(i, 0.);    
    UnfComp->SetBinContent(i, UnfComp->GetBinContent(i)/ UnfOut->GetBinContent(i-gAliFirstBin+2));
    UnfComp->SetBinError(i, UnfComp->GetBinError(i)/ UnfOut->GetBinContent(i-gAliFirstBin+2));
  }
  UnfComp->SetAxisRange(0.8, 1.2, "Y");
  UnfComp->SetAxisRange(0., 7.9, "X");
  // End of comarison check
  

  CorrectedSpectrum->SetLineColor(kGreen);
  TH1D * reconstructedSpectrum = (TH1D*) CorrectedSpectrum->Clone(); reconstructedSpectrum->SetName("reconstructedSpectrum");reconstructedSpectrum->SetLineColor(kOrange);
  reconstructedSpectrum->Scale(2.);
  double sum;
  for(int i=gAliFirstBin;i<=gAliLastBin;i++)
  {
    sum=0;
    for(int j=gAliFirstBin;j<=gAliLastBin;j++)
    {
      sum+=gAliR->GetBinContent(j,i)*CorrectedSpectrum->GetBinContent(j);
    }
    sum+= gAliR->GetBinContent(gAliLastBin+1,i)*CorrectedSpectrum->GetBinContent(gAliLastBin)/2.;
    sum+= gAliR->GetBinContent(gAliFirstBin-1,i)*CorrectedSpectrum->GetBinContent(gAliFirstBin)*2.;
    reconstructedSpectrum->SetBinContent(i, sum);
  }
  TH1D *CorrectedSpectrumSys = (TH1D*) CorrectedSpectrum->Clone();CorrectedSpectrumSys->SetName("CorrectedSpectrumSys");
  // Unfolding has been done. Now use the relative errors from the original distribution
  
  for(int i=gAliFirstBin;i<=gAliLastBin;i++)
  {
    //CorrectedSpectrum->SetBinError(i, CorrectedSpectrum->GetBinContent(i)/spectrum->GetBinContent(i)*spectrum->GetBinError(i));
    CorrectedSpectrumSys->SetBinError(i, CorrectedSpectrumSys->GetBinContent(i)/spectrumsys->GetBinContent(i)*spectrumsys->GetBinError(i));
  }
  gAliUnfoldingSpectrum->Divide(reconstructedSpectrum);
  
  
  // Here the spectrum has been unfolded, now efficiency correction!
  // Check later that binning is the same!
  double sysAdd = sysPID; // From inclusive publication, other systematic errors
  //double ErrorHadronTest = sysHadrCont;
  
  //double EnergyScaling = sysEnergyScaling;
  //sysAdd = TMath::Sqrt(sysAdd*sysAdd+ErrorHadronTest*ErrorHadronTest+EnergyScaling*EnergyScaling+ErrorBSpec*ErrorBSpec);
  
  
  double TotalErr=0;
  
  for(int i=gAliFirstBin;i<=gAliLastBin;i++)
  {
    TotalErr = UncertPID->GetBinContent(i)*UncertPID->GetBinContent(i); // PID
    if(i==gAliFirstBin) TotalErr += additionalErrorInFirstBin*additionalErrorInFirstBin; // unfolding
    if(i==gAliLastBin) TotalErr +=additionalErrorInLastBin*additionalErrorInLastBin;
    TotalErr += UncertIPRes->GetBinContent(i)*UncertIPRes->GetBinContent(i);
    TotalErr += UncertV0Shift->GetBinContent(i)*UncertV0Shift->GetBinContent(i);
    TotalErr += UncertHadrCont->GetBinContent(i)*UncertHadrCont->GetBinContent(i);
    TotalErr += UncertErrorBSpec->GetBinContent(i)*UncertErrorBSpec->GetBinContent(i);
    TotalErr += UncertStrangeness->GetBinContent(i)*UncertStrangeness->GetBinContent(i);
    TotalErr=TMath::Sqrt(TotalErr);
    
    spectrum->SetBinContent(i, spectrum->GetBinContent(i)/efficiency->GetBinContent(i));
    spectrum->SetBinError(i, spectrum->GetBinError(i)/efficiency->GetBinContent(i));
    spectrumsys->SetBinContent(i, spectrumsys->GetBinContent(i)/efficiency->GetBinContent(i));
    spectrumsys->SetBinError(i, TMath::Sqrt(spectrumsys->GetBinError(i)*spectrumsys->GetBinError(i) +
                                                      spectrumsys->GetBinContent(i)*spectrumsys->GetBinContent(i))/efficiency->GetBinContent(i));
    CorrectedSpectrum->SetBinContent(i, CorrectedSpectrum->GetBinContent(i)/efficiency->GetBinContent(i));
    CorrectedSpectrum->SetBinError(i, CorrectedSpectrum->GetBinError(i)/efficiency->GetBinContent(i));
    reconstructedSpectrum->SetBinContent(i, reconstructedSpectrum->GetBinContent(i)/efficiency->GetBinContent(i));
    reconstructedSpectrum->SetBinError(i, reconstructedSpectrum->GetBinError(i)/efficiency->GetBinContent(i));
    CorrectedSpectrumSys->SetBinContent(i, CorrectedSpectrumSys->GetBinContent(i)/efficiency->GetBinContent(i));
    CorrectedSpectrumSys->SetBinError(i, CorrectedSpectrumSys->GetBinError(i)/efficiency->GetBinContent(i));
    CorrectedSpectrumSys->SetBinError(i, TMath::Sqrt(CorrectedSpectrumSys->GetBinError(i)*CorrectedSpectrumSys->GetBinError(i) + 
                                                      CorrectedSpectrumSys->GetBinContent(i)*CorrectedSpectrumSys->GetBinContent(i)*TotalErr*TotalErr));
  }
  // Here unfolding and efficiency correction are applied
  // now do bin width correction!

  for(int i=gAliFirstBin;i<=gAliLastBin;i++)
  {
    spectrum->SetBinContent(i, spectrum->GetBinContent(i)/spectrum->GetBinWidth(i));
    spectrum->SetBinError(i, spectrum->GetBinError(i)/spectrum->GetBinWidth(i));
    spectrumsys->SetBinContent(i, spectrumsys->GetBinContent(i)/spectrum->GetBinWidth(i));
    spectrumsys->SetBinError(i, spectrumsys->GetBinError(i)/spectrum->GetBinWidth(i));
    CorrectedSpectrum->SetBinContent(i, CorrectedSpectrum->GetBinContent(i)/CorrectedSpectrum->GetBinWidth(i));
    CorrectedSpectrum->SetBinError(i, CorrectedSpectrum->GetBinError(i)/CorrectedSpectrum->GetBinWidth(i));
    reconstructedSpectrum->SetBinContent(i, reconstructedSpectrum->GetBinContent(i)/reconstructedSpectrum->GetBinWidth(i));
    reconstructedSpectrum->SetBinError(i, reconstructedSpectrum->GetBinError(i)/reconstructedSpectrum->GetBinWidth(i));
    CorrectedSpectrumSys->SetBinContent(i, CorrectedSpectrumSys->GetBinContent(i)/CorrectedSpectrumSys->GetBinWidth(i));
    CorrectedSpectrumSys->SetBinError(i, CorrectedSpectrumSys->GetBinError(i)/CorrectedSpectrumSys->GetBinWidth(i));
  }
  // Unfolding, efficiency correction and bin width correction are applied.
  // Now add the constant factors for cross section
  
  bool divideBypT = true;
  double twopi =  3.14159 * 2.;
  //double crossSection = 62.2; // mb
  double TAA = 18.915; // mb
  double correctionFactor = 1./(nEvents*1.6*twopi)/2.;  // e+ / e- also
  
  if(divideBypT)
  for(int i=gAliFirstBin;i<=gAliLastBin;i++)
  {
    spectrum->SetBinContent(i, spectrum->GetBinContent(i)/spectrum->GetBinCenter(i)*correctionFactor);
    spectrum->SetBinError(i, spectrum->GetBinError(i)/spectrum->GetBinCenter(i)*correctionFactor);
    spectrumsys->SetBinContent(i, spectrumsys->GetBinContent(i)/spectrum->GetBinCenter(i)*correctionFactor);
    spectrumsys->SetBinError(i, spectrumsys->GetBinError(i)/spectrum->GetBinCenter(i)*correctionFactor);
    CorrectedSpectrum->SetBinContent(i, CorrectedSpectrum->GetBinContent(i)/spectrum->GetBinCenter(i)*correctionFactor);
    CorrectedSpectrum->SetBinError(i, CorrectedSpectrum->GetBinError(i)/spectrum->GetBinCenter(i)*correctionFactor);
    reconstructedSpectrum->SetBinContent(i, reconstructedSpectrum->GetBinContent(i)/spectrum->GetBinCenter(i)*correctionFactor);
    reconstructedSpectrum->SetBinError(i, reconstructedSpectrum->GetBinError(i)/spectrum->GetBinCenter(i)*correctionFactor);
    CorrectedSpectrumSys->SetBinContent(i, CorrectedSpectrumSys->GetBinContent(i)/spectrum->GetBinCenter(i)*correctionFactor);
    CorrectedSpectrumSys->SetBinError(i, CorrectedSpectrumSys->GetBinError(i)/spectrum->GetBinCenter(i)*correctionFactor);
  }
  
  // now add error from shift analysis!
  
  double pT=0;
  for(int i=gAliFirstBin;i<=gAliLastBin;i++)
  {
    pT = spectrumsys->GetXaxis()->GetBinCenter(i);
    
    spectrumsys->SetBinError(i, TMath::Sqrt(spectrumsys->GetBinError(i)*spectrumsys->GetBinError(i) + spectrumsys->GetBinContent(i)*shiftError->Eval(pT)*spectrumsys->GetBinContent(i)*shiftError->Eval(pT)));
    
    //CorrectedSpectrumSys->SetBinError(i, CorrectedSpectrumSys->GetBinError(i)/spectrum->GetBinCenter(i)*correctionFactor);
    //cout << "old Error: " << CorrectedSpectrumSys->GetBinError(i)/CorrectedSpectrumSys->GetBinContent(i) << endl;
    //cout << "change: " << shiftError->Eval(pT) << endl;
    CorrectedSpectrumSys->SetBinError(i, TMath::Sqrt(CorrectedSpectrumSys->GetBinError(i)*CorrectedSpectrumSys->GetBinError(i) + CorrectedSpectrumSys->GetBinContent(i)*shiftError->Eval(pT)*CorrectedSpectrumSys->GetBinContent(i)*shiftError->Eval(pT)));
    //cout << "new Error: " << CorrectedSpectrumSys->GetBinError(i)/CorrectedSpectrumSys->GetBinContent(i) << endl;
  }
  // end of error from shift analysis
  
  
  //spectrum->Divide(efficiency);
  //spectrumsys->Divide(efficiency);
  
  /*spectrum->Draw("e");
  //spectrumsys->Draw("SAME");
  //gAliUnfoldingSpectrum->Draw("SAME");
  CorrectedSpectrum->Draw("SAME");
  reconstructedSpectrum->Draw("SAME");
  gPad->SetLogy();*/
  
  
  
  
  setOPT_hists((TH1F*)CorrectedSpectrum, "p_{T} (GeV/c)", "efficiency", 510, 20, 1.3, kRed);
  setOPT_hists((TH1F*)CorrectedSpectrumSys, "p_{T} (GeV/c)", "efficiency", 510, 20, 1.3, kRed);
  // now make TGraphs from the histograms!
  TGraphErrors * BSpectrumG = new TGraphErrors(gAliLastBin-gAliFirstBin+1);
  TGraphErrors * BSpectrumGSys = new TGraphErrors(gAliLastBin-gAliFirstBin+1);
  for(int i=gAliFirstBin;i<=gAliLastBin;i++)
  {
    BSpectrumG->SetPoint(i-gAliFirstBin+1, CorrectedSpectrum->GetBinCenter(i), CorrectedSpectrum->GetBinContent(i));
    BSpectrumG->SetPointError(i-gAliFirstBin+1, CorrectedSpectrum->GetBinWidth(i)/2., CorrectedSpectrum->GetBinError(i));
    BSpectrumGSys->SetPoint(i-gAliFirstBin+1, CorrectedSpectrumSys->GetBinCenter(i), CorrectedSpectrumSys->GetBinContent(i));
    BSpectrumGSys->SetPointError(i-gAliFirstBin+1, CorrectedSpectrumSys->GetBinWidth(i)/2., CorrectedSpectrumSys->GetBinError(i));
  }
  setOPT_graph(BSpectrumG, "p_{T}", "dN/dp_{T}");
  setOPT_graph(BSpectrumGSys, "p_{T}", "dN/dp_{T}");
  BSpectrumG->SetMarkerColor(kRed);
  BSpectrumG->SetMarkerSize(1.0);
  BSpectrumGSys->SetMarkerSize(0);
  BSpectrumGSys->SetFillColor(kRed);
  BSpectrumGSys->SetFillStyle(3002);
  
  TH1D * emptyhist = new TH1D("emptyhist", "", 100, 0., 8.);
  //setOPT_hists((TH1F*)emptyhist, "p_{T} (GeV/c)", "#frac{1}{2 #pi p_{T}} #frac{d^{2}#sigma}{dp_{T}d#eta}", 510, 20, 1.3, kBlack);
  setOPT_hists((TH1F*)emptyhist, "p_{T} (GeV/c)", "#frac{1}{2 #pi p_{T} N_{coll}} #frac{d^{2}N}{dp_{T}d#eta}", 510, 20, 1.3, kBlack);
  emptyhist->SetAxisRange(0., 8., "X");
  emptyhist->SetAxisRange(CorrectedSpectrum->GetBinContent(gAliLastBin)/2., CorrectedSpectrum->GetMaximum()*2., "Y");
  new TCanvas;
  emptyhist->Draw();
  BSpectrumGSys->Draw("SAMEP2");
  BSpectrumG->Draw("SAMEP");
  //CorrectedSpectrum->Draw("SAME");
  gPad->SetLogy();
  
  
  
  TFile * MJfile = new TFile("MJoutput/MJb2eptspectra.root");
  TH1D * MJSpectrum = (TH1D*) MJfile->Get("hstat");
  TGraphAsymmErrors *  MJSyst= (TGraphAsymmErrors*) MJfile->Get("gsys");
  setOPT_graph(MJSyst, "p_{T}", "dN/dp_{T}");
  MJSyst->SetMarkerSize(0);
  MJSyst->SetFillColor(kGreen);
  MJSyst->SetFillStyle(3002);
  TH1D * MJSystH = (TH1D*) MJSpectrum->Clone(); MJSystH->SetName("MJSystH"); MJSystH->Scale(0);
  double tx,ty;
  for(int i=0;i<MJSyst->GetN();i++)
  {
    MJSyst->GetPoint(i,tx,ty);
    MJSystH->SetBinContent(MJSpectrum->FindBin(tx), ty);
    MJSystH->SetBinError(MJSpectrum->FindBin(tx), MJSyst->GetErrorY(i));
  }
  setOPT_hists((TH1F*)MJSpectrum, "p_{T} (GeV/c)", "spectrum", 510, 20, 1.0, kGreen);
  
  TH1D * rebinnedMJoutput = MakeComparisonDiagram(MJSpectrum, CorrectedSpectrum);
  rebinnedMJoutput->SetBinContent(6,0);
  TH1D * rebinnedMJoutputSys = MakeComparisonDiagramCorrelated(MJSystH, CorrectedSpectrum);
  rebinnedMJoutputSys->SetBinContent(6,0);
  //setOPT_hists((TH1F*)rebinnedMJoutputSys, "p_{T} (GeV/c)", "efficiency", 510, 20, 0.0, kBlue);
  rebinnedMJoutputSys->SetMarkerSize(0.1);rebinnedMJoutputSys->SetMarkerColor(kBlue);rebinnedMJoutputSys->SetLineColor(kBlue);
  setOPT_hists((TH1F*)rebinnedMJoutput, "p_{T} (GeV/c)", "spectrum", 510, 20, 1.0, kGreen+2);
  
  for(int i=gAliFirstBin;i<=gAliLastBin;i++)
  {
   Uncertpp->SetBinContent(i, rebinnedMJoutputSys->GetBinError(i)/rebinnedMJoutputSys->GetBinContent(i));
  }
  
  // Correct Spectrum for Energy difference!
  TF1 *scalefct = new TF1("scalefct","[0]/TMath::Power(([1]+x/[2]),[3])+[4]*TMath::Exp(-x/[5])",0,30);
  scalefct->SetParameter(0,0.6379);
  scalefct->SetParameter(1,1.44);
  scalefct->SetParameter(2,6.503);
  scalefct->SetParameter(3,0.7975);
  scalefct->SetParameter(4,0.2135);
  scalefct->SetParameter(5,1.37);
  rebinnedMJoutput->Multiply(scalefct);
  rebinnedMJoutputSys->Multiply(scalefct);  
  // end of correction
  
  
  //rebinnedMJoutputSys->Draw("SAMEE1");
  //rebinnedMJoutput->Draw("SAME");
  
   TFile * inShingoRAA = new TFile("ShingoRAA/Raa010Rebin.root");
  TGraphErrors * ShingoRAA = (TGraphErrors *) inShingoRAA->Get("RaaComb");
  TGraphAsymmErrors * sysShingoRAA = (TGraphAsymmErrors *) inShingoRAA->Get("RaaCombSys");
  setOPT_graph(ShingoRAA, "p_{T}", "dN/dp_{T}");
  setOPT_graph(sysShingoRAA, "p_{T}", "dN/dp_{T}");
  ShingoRAA->SetMarkerColor(kBlue);
  ShingoRAA->SetMarkerSize(1.0);
  sysShingoRAA->SetMarkerSize(0);
  sysShingoRAA->SetFillColor(kBlue);
  sysShingoRAA->SetFillStyle(3003);
  
  
 new TCanvas;
  TH1D * RAA = (TH1D*) CorrectedSpectrum->Clone(); RAA->SetName("RAA");
  TH1D * RAAsys = (TH1D*) CorrectedSpectrumSys->Clone(); RAAsys->SetName("RAAsys");
  RAA->Scale(1./TAA);
  RAAsys->Scale(1./TAA);
  RAA->Divide(rebinnedMJoutput);
  RAAsys->Divide(rebinnedMJoutputSys);
  
  TGraphErrors * RAAG = new TGraphErrors(gAliLastBin-gAliFirstBin+1);
  TGraphErrors * RAAGSys = new TGraphErrors(gAliLastBin-gAliFirstBin+1);
  for(int i=gAliFirstBin;i<=gAliLastBin;i++)
  {
    RAAG->SetPoint(i-gAliFirstBin+1, RAA->GetBinCenter(i), RAA->GetBinContent(i));
    RAAG->SetPointError(i-gAliFirstBin+1, RAA->GetBinWidth(i)/2., RAA->GetBinError(i));
    RAAGSys->SetPoint(i-gAliFirstBin+1, RAAsys->GetBinCenter(i), RAAsys->GetBinContent(i));
    RAAGSys->SetPointError(i-gAliFirstBin+1, RAAsys->GetBinWidth(i)/2., RAAsys->GetBinError(i));
  }
  
  TH1D * emptyhist2 = new TH1D("emptyhist2", "", 100, 0., 8.);
  setOPT_hists((TH1F*)emptyhist2, "p_{T} (GeV/c)", "R_{AA}", 510, 20, 1.3, kBlack);
  emptyhist2->SetAxisRange(0., 8., "X");
  emptyhist2->SetAxisRange(0., 2.0, "Y");
  setOPT_graph(RAAG, "p_{T}", "dN/dp_{T}");
  setOPT_graph(RAAGSys, "p_{T}", "dN/dp_{T}");
  RAAG->SetMarkerColor(kRed);
  RAAG->SetMarkerSize(1.0);
  RAAGSys->SetMarkerSize(0);
  RAAGSys->SetFillColor(kRed);
  RAAGSys->SetFillStyle(3002);
  
  emptyhist2->Draw();
  //sysShingoRAA->Draw("SAMEP2");
  //ShingoRAA->Draw("SAMEP");
  RAAGSys->Draw("SAMEP2");
  RAAG->Draw("SAMEP");
  
  
  
  TLegend * legRatio=plotLegend("right_top"," ",0.7,0.5,0.0,+0.05, Form("Beauty R_{AA}"),1);
  legRatio->AddEntry(RAA, "RAA Fit Method 0-20%");
  //legRatio->AddEntry(ShingoRAA, "RAA TPC-EMCal (D+B) 0-10%");
  legRatio->Draw("SAME");
  
  TFile * outFile = new TFile(OutputFileName->Data(), "RECREATE");
  CorrectedSpectrum->Write();
  CorrectedSpectrumSys->Write();
  BSpectrumG->Write();
  BSpectrumGSys->Write();
  RAA->Write();
  RAAsys->Write();
  
  new TCanvas;
 gAliR->SetAxisRange(0.1,6.,"X");
 gAliR->SetAxisRange(0.1,6.,"Y");
  setOPT_hists((TH1F*)gAliR, "p_{T, true} (GeV/c)", "p_{T, reconstructed} (GeV/c)", 510, 20, 1.0, kGreen+2);
 gAliR->Draw("COLZ");
  
  // Load necessary libraries
  LoadLibs();
  // Set the default style
  SetStyle();
  
  
  TCanvas *cfig = new TCanvas("systematics", "Systematic Uncertainties", 800, 600); 
  //UncertSum  
  double totsyserr=0;
  for(int i=gAliFirstBin;i<=gAliLastBin;i++)
  {
    totsyserr=0;
    totsyserr += UncertPID->GetBinContent(i)*UncertPID->GetBinContent(i);
    totsyserr += UncertIPRes->GetBinContent(i)*UncertIPRes->GetBinContent(i);
    totsyserr += UncertV0Shift->GetBinContent(i)*UncertV0Shift->GetBinContent(i);
    totsyserr += UncertHadrCont->GetBinContent(i)*UncertHadrCont->GetBinContent(i);
    totsyserr += UncertErrorBSpec->GetBinContent(i)*UncertErrorBSpec->GetBinContent(i);
    totsyserr += UncertStrangeness->GetBinContent(i)*UncertStrangeness->GetBinContent(i);
    totsyserr += UncertUnfolding->GetBinContent(i)*UncertUnfolding->GetBinContent(i);
    totsyserr += UncertFit->GetBinContent(i)*UncertFit->GetBinContent(i);
    totsyserr += UncertEnergyScaling->GetBinContent(i)*UncertEnergyScaling->GetBinContent(i);
    totsyserr += Uncertpp->GetBinContent(i)*Uncertpp->GetBinContent(i);
    UncertSum->SetBinContent(i, TMath::Sqrt(totsyserr));
    
  }
  
  TH1 * h = cfig->DrawFrame(0,0.0001,11,0.6);
  // === Commonly used x/ titles: ===
  // pt invariant yields
  const char *  texPtX="#it{p}_{T} (GeV/#it{c})";
  const char *  texPtY="rel. sys. uncertainty";
  h->SetXTitle(texPtX);
  h->SetYTitle(texPtY);
  
  TH1D * emptyhistsys = new TH1D("emptyhist", "", 100, 0., 8.);
  setOPT_hists((TH1F*)emptyhistsys, "p_{T} (GeV/c)", "rel. sys. uncertainty", 510, 20, 1.3, kBlack);
  emptyhistsys->SetAxisRange(0., 11., "X");
  emptyhistsys->SetAxisRange(0., 0.5, "Y");
  emptyhistsys->Draw("SAME");
  UncertPID->Draw("SAME");
  UncertIPRes->Draw("SAME");
  UncertV0Shift->Draw("SAME");
  UncertHadrCont->Draw("SAME");
  UncertErrorBSpec->Draw("SAME");
  UncertStrangeness->Draw("SAME");
  UncertUnfolding->Draw("SAME");
  UncertEnergyScaling->Draw("SAME");
  UncertFit->Draw("SAME");
  Uncertpp->Draw("SAME");
  UncertSum->Draw("SAME");
  TLegend * legsys=plotLegend("right_bottom"," ",1.0,1.0,0.04,+0.14, Form(""),1);
  legsys->AddEntry(UncertFit, "Fit", "L");
  legsys->AddEntry(UncertPID, "PID, Tracking", "L");
  legsys->AddEntry(UncertIPRes, "IP Resolution", "L");
  legsys->AddEntry(UncertV0Shift, "Conversion Electrons", "L");
  legsys->AddEntry(UncertHadrCont, "Hadron Contamination", "L");
  legsys->AddEntry(UncertErrorBSpec, "HF Mother spectra", "L");
  legsys->AddEntry(UncertStrangeness, "Strangeness", "L");
  legsys->AddEntry(UncertUnfolding, "Unfolding", "L");
  legsys->AddEntry(UncertEnergyScaling, "Energy Scaling", "L");
  legsys->AddEntry(Uncertpp, "pp Reference", "L");
  legsys->AddEntry(UncertSum, "total", "L");
  legsys->Draw("SAME");
  
  TLatex *   tex = new TLatex(0.2, 0.81, 1. ? "ALICE Preliminary" : "ALICE");
  tex->SetNDC();
  tex->SetTextFont(42);
  //tex->Draw("SAME");
  TLatex *   texSystem1 = new TLatex(0.58, 0.72, "Pb-Pb, 0-20% centrality");
  texSystem1->SetNDC();
  texSystem1->SetTextFont(42);
  texSystem1->SetTextSize(0.04);
  texSystem1->Draw("SAME");
  TLatex *   texSystem2 = new TLatex(0.61, 0.67, "#sqrt{#it{s}}_{NN} = 2.76 TeV");
  texSystem2->SetNDC();
  texSystem2->SetTextFont(42);
  texSystem2->SetTextSize(0.04);
  texSystem2->Draw("SAME");
  
  
  new TCanvas;
  TLine * LowpTLine = new TLine(1.3, 0.8, 1.3, 1.2); LowpTLine->SetLineColor(kRed);
  TLine * HighpTLine = new TLine(6.0, 0.8, 6.0, 1.2); HighpTLine->SetLineColor(kRed);
  TLegend * legUnfolding=plotLegend("left_top"," ",1.5,0.5,0.15,+0.05, Form(""),1);
  setOPT_hists((TH1F*)UnfComp, "p_{T} (GeV/c)", "ratio", 510, 20, 1.3, kBlack);
  legUnfolding->AddEntry(UnfComp, "AliHFEUnfolding/TUnfold (full Response)");
  UnfComp->Draw();
  LowpTLine->Draw("SAME");
  HighpTLine->Draw("SAME");
  legUnfolding->Draw("SAME");
  
}