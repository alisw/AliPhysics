#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"

#include "AliCFContainer.h"
#include "AliCFEffGrid.h"
#include "AliCFDataGrid.h"
#include "AliPID.h"

#include <iostream>

#include "THnSparseDefinitions.h"

enum type { kTrackPt = 0, kZ = 1, kXi = 2, kDistance = 3, kJT = 4, kNtypes = 5 };
enum typeMCSysErrors { kNoErrors = 0, kErrorsWithoutMultDep = 1, kErrorsIncludingMultDep = 2, kErrorsOnlyMultDep = 3, kErrorsForMerging = 4,
                       kErrorsForMergingOnlyMultDep = 5 };

const Int_t numParamsMult = 3;

const TString obsString[kNtypes] = {"Pt", "z", "#xi", "R", "Jt"};
const TString obsStringMCsysError[kNtypes] = {"Pt", "Z", "Xi", "R", "Jt"};

Int_t iPt     = 0;
Int_t iMCid   = 0;
Int_t iEta    = 0;
Int_t iCharge = 0;
Int_t iMult   = 0;
Int_t iJetPt  = 0;
Int_t iZ = 0;
Int_t iXi = 0;
Int_t iDistance = 0;
Int_t iJt = 0;

Int_t iObsAxis = 0;

// Values from Xianguo:
Double_t PsNLHC10eDATA = 131136086;
Double_t PsVrecNLHC10eDATA = 118545498;
Double_t PsVrecVzNLHC10eDATA = 104749839;
Double_t PsNLHC10dDATA = 158428731;
Double_t PsVrecNLHC10dDATA = 143011819;
Double_t PsVrecVzNLHC10dDATA = 127089415;


//const Double_t corrNevVertexcut = (PsVrecNLHC10eDATA + PsVrecNLHC10dDATA) / (PsNLHC10eDATA + PsNLHC10dDATA); 

const Double_t corrNevVertexcut = 0.900896856759954789; // My value (agrees on permille level with the one of Xianguo)
//const Double_t corrNevVertexcut = 0.900790018966733386; // My value for 10d only (agrees on permille level with the one of Xianguo and on sub-permille level with my value for 10d+10e)


// Systmatic uncertainties for (residual) multiplicity of correction factors:
const Double_t sysErrMultDepSecCorr_pr = 0.01;
const Double_t sysErrMultDepEffCorr = 0.01;
const Double_t sysErrMultDepEffCorrLowestMult = 0.03;




// Values for the efficiency correction 10d<->10e:
// (values obtained for V0M analysis, but ratios 10d/10e ~the same for ref. mult.)
const Double_t statFor10d = 2.085004e8;
const Double_t statFor10e = 1.676471e8;
const Double_t statFor10de = statFor10d + statFor10e;
const Double_t ratioStat10d_10de = statFor10d / statFor10de;
const Double_t ratioStat10e_10de = statFor10e / statFor10de;

const Double_t effRatio10eOver10d = 1. - 0.023;
const Double_t effCorrFactor10d10e = 1. * ratioStat10d_10de + effRatio10eOver10d * ratioStat10e_10de;


// Values for the relative systematic uncertainties of the MC secondary correction
Double_t secCorrRelSysError[kNtypes][AliPID::kSPECIES]; // Initialised in InitialiseRelSysErrOfSecCorr()!

// New method with division by kaons, fitting etc and using only the low-pT part. z is just taken from pT and assumed to be the same.
//**** Inclusive with track bit 4096
const Double_t secCorrRelSysErrPtInclusive[AliPID::kSPECIES] = { 0.045, 0.045, 0.381, 0.381, 0.157 };
//**** Inclusive with track bit 528
//const Double_t secCorrRelSysErrPtInclusive[AliPID::kSPECIES] = { 0.035, 0.035, 0.282, 0.282, 0.101 };

//**** Jets with track bit 528 (corr with 14b6)
const Double_t secCorrRelSysErrPtJets[AliPID::kSPECIES]      = { 0.006, 0.006, 0.269, 0.269, 0.163 };
const Double_t secCorrRelSysErrZ[AliPID::kSPECIES]           = { 0.006, 0.006, 0.269, 0.269, 0.163 };
//**** Jets with track bit 528 (corr with 10f6a)
//const Double_t secCorrRelSysErrPtJets[AliPID::kSPECIES]      = { 0.004, 0.004, 0.234, 0.234, 0.140 };
//const Double_t secCorrRelSysErrZ[AliPID::kSPECIES]           = { 0.004, 0.004, 0.234, 0.234, 0.140 };
//**** Jets with track bit 272
//const Double_t secCorrRelSysErrPtJets[AliPID::kSPECIES]      = { 0.004, 0.004, 0.161, 0.161, 0.090 };
//const Double_t secCorrRelSysErrZ[AliPID::kSPECIES]           = { 0.004, 0.004, 0.161, 0.161, 0.090 };


const Double_t secCorrRelSysErrXi[AliPID::kSPECIES] = { 0.22, 0.22, 0.22, 0.22, 0.22};// No fits from my side, just stick to Oliver's values

//TODO FULLJETS determine secondary correction rel errors. At the moment (since no serious data set used), just take dummy z values
const Double_t secCorrRelSysErrDistance[AliPID::kSPECIES]           = { 0.006, 0.006, 0.269, 0.269, 0.163 };
const Double_t secCorrRelSysErrJt[AliPID::kSPECIES]                 = { 0.006, 0.006, 0.269, 0.269, 0.163 };

/*OLD
const Double_t secCorrRelSysErrPtInclusive[AliPID::kSPECIES] = { 0.031, 0.031, 0.667, 0.667, 0.172 };
const Double_t secCorrRelSysErrPtJets[AliPID::kSPECIES]      = { 0.008, 0.008, 0.624, 0.624, 0.224 };
const Double_t secCorrRelSysErrZ[AliPID::kSPECIES]  = { 0.010, 0.010, 0.542, 0.542, 0.240 };
const Double_t secCorrRelSysErrXi[AliPID::kSPECIES] = { 0.22, 0.22, 0.22, 0.22, 0.22};// No fits from my side, just stick to Oliver's values
*/
// Olivers values: pT, z, xi = { 0.2, 0.17, 0.22 };

// Values for the relative systematic uncertainty of the MC eff and res correction for inclusive
const Double_t relMCSysErrEffInclusive = 0.05;
const Double_t relMCSysErrResInclusive = 0.02;


// Thresholds for corrFac(Eff,acc,res) and total rel. sys. error. If thresholds are exceeded, data points are rejected
const Double_t thresholdCF = 9./4.;
const Double_t thresholdCFerror = 1./3.;

const Int_t nPtBinsType2 = 44;
const Double_t binsPtType2[nPtBinsType2+1] = {0., 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45,
          0.5, 0.55, 0.6, 0.65, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4,
          1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.4, 3.8,
          4.5, 5.5, 6.5, 8.0, 10.0, 12.0, 14.0, 16.0, 20.0, 24.0,
          28.0, 32.0, 36.0, 45.0, 50.0 };

//___________________________________________________________________
void setupHist(TH1* h, TString histName, TString histTitle, TString xAxisTitle, TString yAxisTitle, Int_t color)
{
  if (histName != "")
    h->SetName(histName.Data());
  h->SetTitle(histTitle.Data());
  
  if (xAxisTitle != "")
    h->GetXaxis()->SetTitle(xAxisTitle.Data());
  if (yAxisTitle != "")
    h->GetYaxis()->SetTitle(yAxisTitle.Data());
  
  h->SetMarkerStyle(24);
  h->SetLineColor(color);
  h->SetMarkerColor(color);
  
  h->SetStats(kFALSE);
}



//___________________________________________________________________
void setupHistMCrelSysErr(TH1* h, TH1* hStyle)
{
  if (!h)
    return;
  
  if (hStyle) {
    const Int_t lineColour = hStyle->GetLineColor();
    const Int_t fillColour = hStyle->GetFillColor();
    const Int_t markerColour = hStyle->GetMarkerColor();
    const Int_t lineWidth = hStyle->GetLineWidth();
    
    h->SetLineColor(lineColour);
    h->SetFillColor(fillColour);
    h->SetMarkerColor(markerColour);
    h->SetLineWidth(lineWidth);
  
    TAxis* xAxis = hStyle->GetXaxis();
    h->GetXaxis()->SetTitle(xAxis->GetTitle());
    h->GetXaxis()->SetMoreLogLabels(xAxis->GetMoreLogLabels());
    h->GetXaxis()->SetNoExponent(xAxis->GetNoExponent());
  }
  
  h->GetXaxis()->SetLabelSize(0.05);
  h->GetYaxis()->SetLabelSize(0.05);
  h->GetXaxis()->SetTitleSize(0.06);
  h->GetYaxis()->SetTitleSize(0.06);
  h->GetXaxis()->SetTitleOffset(1.0);
  h->GetYaxis()->SetTitleOffset(0.77);
  
  h->GetYaxis()->SetTitle("Rel. Sys. Error");
  
  
  for (Int_t i = 1; i <= h->GetNbinsX(); i++)
    h->SetBinError(i, 0);
  
  h->SetMinimum();
  h->SetMaximum();
  h->GetYaxis()->SetRange(0, -1);
}


//___________________________________________________________________
Bool_t HasNonVanishingContent(TH1D* h)
{
  if (h) {
    for (Int_t i = 1; i <= h->GetNbinsX(); i++) {
      if (h->GetBinContent(i) > 0)
        return kTRUE;
    }
  }
  
  return kFALSE;
}


//___________________________________________________________________
void InitialiseRelSysErrOfSecCorr(Bool_t inclusive)
{
  // Initialise the rel. systematic errors for the secondary correction factors

  printf("Initialising rel. sys. errors of sec. corr. factors for %s...\n", inclusive ? "inclusive" : "jets");
  
  for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
    secCorrRelSysError[kTrackPt][species] = inclusive ? secCorrRelSysErrPtInclusive[species] : secCorrRelSysErrPtJets[species];
    secCorrRelSysError[kZ][species] = secCorrRelSysErrZ[species];
    secCorrRelSysError[kXi][species] = secCorrRelSysErrXi[species];
    secCorrRelSysError[kDistance][species] = secCorrRelSysErrDistance[species];
    secCorrRelSysError[kJT][species] = secCorrRelSysErrJt[species];
  }
}


//___________________________________________________________________
//------ For plotting the sys errors stacked -------
TCanvas* plotSysErrorsStacked(Int_t selectSp, TH1D* hSysErrEff[], TH1D* hSysErrRes[], TH1D* hSysErrBbB[], TH1D* hSysErrSec[],
                              TH1D* hSysErrSS[], TH1D* hSysErrGFL[], TH1D* hSysErrTot[], TH1D* hSysErrEff10d10e[], TH1D* hSysErrMu,
                              TH1D* hSysErrSecMultDep[], TH1D* hSysErrEffMultDep[],
                              const TString canvName, TString header, const TString saveName, Double_t jetPtUp, Int_t iObs,
                              typeMCSysErrors sysErrorTypeMC)
{
  // Plot stacked errors.
  
  // plot only relative errors 
  // stacked represenation adding errros linear, show also quadratic sum
  
  Bool_t setLogX = iObs == kTrackPt;
  TH1D* hSysErrAccEff[AliPID::kSPECIES] = { 0x0, };
  TH1D* hSysErrAccEff10d10e[AliPID::kSPECIES] = { 0x0, };
  TH1D* hSysErrAccRes[AliPID::kSPECIES] = { 0x0, };
  TH1D* hSysErrAccBbB[AliPID::kSPECIES] = { 0x0, };
  TH1D* hSysErrAccSec[AliPID::kSPECIES] = { 0x0, };
  TH1D* hSysErrAccSS[AliPID::kSPECIES] = { 0x0, };
  TH1D* hSysErrAccGFL[AliPID::kSPECIES] = { 0x0, };
  
  TH1D* hSysErrAccSecMultDep[AliPID::kSPECIES] = { 0x0, };
  TH1D* hSysErrAccEffMultDep[AliPID::kSPECIES] = { 0x0, };

  TH1D* hSysErrAccMu = 0x0; 
  
  TH1D* hSysErrAccTot[AliPID::kSPECIES] = { 0x0, };
  
  
  // Special treatment if only the sys errors due to multiplicity depend are to be plotted
  const Bool_t onlyMult = (sysErrorTypeMC == kErrorsOnlyMultDep || sysErrorTypeMC == kErrorsForMergingOnlyMultDep);
  
  TH1D* hLast = 0x0;
  Int_t numHists = 0;
  
  const Int_t lineWidthStacked = 1;
  const Int_t sp = selectSp;
  
  Bool_t hasNonVanishingSysErrEff = HasNonVanishingContent(hSysErrEff[sp]);
  Bool_t hasNonVanishingSysErrRes = HasNonVanishingContent(hSysErrRes[sp]);
  Bool_t hasNonVanishingSysErrGFL = HasNonVanishingContent(hSysErrGFL[sp]);
  
  if (!onlyMult) {
    if (hasNonVanishingSysErrEff) {
      hSysErrEff[sp]->SetStats(0);
      hSysErrAccEff[sp] = (TH1D*) hSysErrEff[sp]->Clone(Form("hSysErrAccEff_%d",sp));
      numHists++;
      hLast = hSysErrAccEff[sp];
    }
    
    if (hSysErrEff10d10e && hSysErrEff10d10e[sp]) {
      if (!hLast) {
        hSysErrEff10d10e[sp]->SetStats(0);
        hSysErrAccEff10d10e[sp] = (TH1D*) hSysErrEff10d10e[sp]->Clone(Form("hSysErrAccEff10d10e_%d",sp));
      }
      else {
        hSysErrAccEff10d10e[sp] = (TH1D*) hLast->Clone(Form("hSysErrAccEff10d10e_%d",sp));
        hSysErrAccEff10d10e[sp]->Add(hSysErrEff10d10e[sp]);
      }
      numHists++;
      hLast = hSysErrAccEff10d10e[sp];
    }
  }
  
  if (hSysErrEffMultDep && hSysErrEffMultDep[sp]) {
    if (!hLast) {
      hSysErrEffMultDep[sp]->SetStats(0);
      hSysErrAccEffMultDep[sp] = (TH1D*) hSysErrEffMultDep[sp]->Clone(Form("hSysErrAccEffMultDep_%d",sp));
    }
    else {
      hSysErrAccEffMultDep[sp] = (TH1D*) hLast->Clone(Form("hSysErrAccEffMultDep_%d",sp));
      hSysErrAccEffMultDep[sp]->Add(hSysErrEffMultDep[sp]);
    }
    numHists++;
    hLast = hSysErrAccEffMultDep[sp];
  }
  
  if (!onlyMult) {
    if (hasNonVanishingSysErrRes) {
      if (!hLast) {
        hSysErrRes[sp]->SetStats(0);
        hSysErrAccRes[sp] = (TH1D*) hSysErrRes[sp]->Clone(Form("hSysErrAccRes_%d",sp));
      }
      else {
        hSysErrAccRes[sp] = (TH1D*) hLast->Clone(Form("hSysErrAccRes_%d",sp));
        hSysErrAccRes[sp]->Add(hSysErrRes[sp]);
      }
      hLast = hSysErrAccRes[sp];
      numHists++;
    }

    // For inclusive, there is no uncertainty from shape
    if (jetPtUp > 0) {
      if (!hLast) {
        hSysErrBbB[sp]->SetStats(0);
        hSysErrAccBbB[sp] = (TH1D*) hSysErrBbB[sp]->Clone(Form("hSysErrAccBbB_%d",sp));
      }
      else {
        hSysErrAccBbB[sp] = (TH1D*) hLast->Clone(Form("hSysErrAccBbB_%d",sp));
        hSysErrAccBbB[sp]->Add(hSysErrBbB[sp]);
      }
      hLast = hSysErrAccBbB[sp];
      numHists++;
    }
    
    if (!hLast) {
      hSysErrSec[sp]->SetStats(0);
      hSysErrAccSec[sp] = (TH1D*) hSysErrSec[sp]->Clone(Form("hSysErrAccSec_%d",sp));
    }
    else {
      hSysErrAccSec[sp] = (TH1D*) hLast->Clone(Form("hSysErrAccSec_%d",sp));
      hSysErrAccSec[sp]->Add(hSysErrSec[sp]);
    }
    hLast = hSysErrAccSec[sp];
    numHists++;
  }
  
  if (hSysErrSecMultDep && hSysErrSecMultDep[sp]) {
    if (!hLast) {
      hSysErrSecMultDep[sp]->SetStats(0);
      hSysErrAccSecMultDep[sp] = (TH1D*) hSysErrSecMultDep[sp]->Clone(Form("hSysErrAccSecMultDep_%d",sp));
    }
    else {
      hSysErrAccSecMultDep[sp] = (TH1D*) hLast->Clone(Form("hSysErrAccSecMultDep_%d",sp));
      hSysErrAccSecMultDep[sp]->Add(hSysErrSecMultDep[sp]);
    }
    numHists++;
    hLast = hSysErrAccSecMultDep[sp];
  }
  
  if (!onlyMult) {
    if (hSysErrSS[sp]) {
      hSysErrAccSS[sp] = (TH1D*) hLast->Clone(Form("hSysErrAccSS_%d",sp));
      hSysErrAccSS[sp]->Add(hSysErrSS[sp]);
      hLast = hSysErrAccSS[sp];
      numHists++;
    }
    
    if (hSysErrGFL[sp] && hasNonVanishingSysErrGFL) {
      hSysErrAccGFL[sp] = (TH1D*) hLast->Clone(Form("hSysErrAccGFL_%d",sp));
      hSysErrAccGFL[sp]->Add(hSysErrGFL[sp]);
      hLast = hSysErrAccGFL[sp];
      numHists++;
    }
    
    /*if (hSysErrSS[sp] && hSysErrGFL[sp]) {
        hSysErrAccSS[sp] = (TH1D*) hSysErrAccSec[sp]->Clone(Form("hSysErrAccSS_%d",sp));
        hSysErrAccSS[sp]->Add(hSysErrSS[sp]);

        hSysErrAccGFL[sp] = (TH1D*) hSysErrAccSS[sp]->Clone(Form("hSysErrAccGFL_%d",sp));
        hSysErrAccGFL[sp]->Add(hSysErrGFL[sp]);
    }*/
    
    if (hSysErrMu) {
      hSysErrAccMu = (TH1D*) hLast->Clone(Form("%s_%d",hSysErrMu->GetName(), sp));
      hSysErrAccMu->Add(hSysErrMu);
      hLast = hSysErrAccMu;
      numHists++;
    }
  }
  
  
  hSysErrAccTot[sp] = (TH1D*) hSysErrTot[sp]->Clone(Form("hSysErrAccTot_%d",sp));
  numHists++;
  
  if (!hLast)
    hLast = hSysErrTot[sp];
    
  // ---- plots ---

  //gStyle->SetOptStat(0);
  
  TCanvas *c1 = new TCanvas(canvName.Data(),"",760,420);
  c1->SetGrid(0,0);
  c1->SetTopMargin(0.02);
  c1->SetRightMargin(0.018);
  c1->SetLeftMargin(0.09);
  c1->SetBottomMargin(0.14);
  
  c1->SetLogx(setLogX);
  
  const Int_t numColumnsLeg = 2;
  const Int_t numHistRowsLeg = ((Int_t)(numHists / numColumnsLeg)) + ((numHists % numColumnsLeg) > 0);
  
  const Double_t legYLowDefault = 0.6;//0.5;
  const Double_t legYHighDefault = 0.95;
  // 1 row for the header
  Double_t legHeight = (legYHighDefault - legYLowDefault) / (4 + 1) * (numHistRowsLeg + 1);
  TLegend* leg = new TLegend(0.35, legYHighDefault - legHeight, 0.95, legYHighDefault);
  //TLegend* leg = new TLegend(0.28, legYHighDefault - legHeight, 0.95, legYHighDefault);
  leg->SetNColumns(numColumnsLeg);
  
  leg->SetTextSize(0.05);//0.06);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetMargin(0.15);
  
  leg->SetHeader(header.Data());

  TH1D* histForAxis = hLast;
  
  setupHistMCrelSysErr(hLast, 0x0);
  
  histForAxis->GetYaxis()->SetRangeUser(0,0.59);
  if (jetPtUp > 0 && iObs == kTrackPt)
    histForAxis->GetXaxis()->SetRangeUser(0.15, TMath::Max(30, (Int_t)jetPtUp));
  
  if (iObs == kZ)
    histForAxis->GetXaxis()->SetRangeUser(0., 1.);

  TString drawStr = "";
  
  if (hSysErrAccMu) {
    hSysErrAccMu->SetLineWidth(lineWidthStacked);
    hSysErrAccMu->SetFillColor(kGray);
    hSysErrAccMu->SetLineColor(kBlack);
    hSysErrAccMu->Draw(drawStr.Data());
    drawStr = "same";
  }
  
  if (hSysErrAccGFL[selectSp]) {
    hSysErrAccGFL[selectSp]->SetLineWidth(lineWidthStacked);
    hSysErrAccGFL[selectSp]->SetFillColor(9);
    hSysErrAccGFL[selectSp]->SetLineColor(kBlack);
    hSysErrAccGFL[selectSp]->Draw(drawStr.Data());
    drawStr = "same";
  }


  if (hSysErrAccSS[selectSp]) {
    hSysErrAccSS[selectSp]->SetLineWidth(lineWidthStacked);
    hSysErrAccSS[selectSp]->SetFillColor(3);
    hSysErrAccSS[selectSp]->SetLineColor(kBlack);
    hSysErrAccSS[selectSp]->Draw(drawStr.Data());
    drawStr = "same";
  }
  
  if (hSysErrAccSecMultDep[selectSp]) {
    hSysErrAccSecMultDep[selectSp]->SetLineWidth(lineWidthStacked);
    hSysErrAccSecMultDep[selectSp]->SetFillColor(kOrange - 3);
    hSysErrAccSecMultDep[selectSp]->SetLineColor(kBlack);
    hSysErrAccSecMultDep[selectSp]->Draw(drawStr.Data());
    drawStr = "same";
  }
  
  if (hSysErrAccSec[selectSp]) {
    hSysErrAccSec[selectSp]->SetLineWidth(lineWidthStacked);
    hSysErrAccSec[selectSp]->SetFillColor(5);
    hSysErrAccSec[selectSp]->SetLineColor(kBlack);
    hSysErrAccSec[selectSp]->Draw(drawStr.Data());
    drawStr = "same";
  }
  
  if (hSysErrAccBbB[selectSp]) {
    hSysErrAccBbB[selectSp]->SetLineWidth(lineWidthStacked);
    hSysErrAccBbB[selectSp]->SetFillColor(7);
    hSysErrAccBbB[selectSp]->SetLineColor(kBlack);
    hSysErrAccBbB[selectSp]->Draw(drawStr.Data());
    drawStr = "same";
  }
  
  if (hSysErrAccRes[selectSp]) {
    hSysErrAccRes[selectSp]->SetLineWidth(lineWidthStacked);
    hSysErrAccRes[selectSp]->SetFillColor(4);
    hSysErrAccRes[selectSp]->SetLineColor(kBlack);
    hSysErrAccRes[selectSp]->Draw(drawStr.Data());
    drawStr = "same";
  }
  
  if (hSysErrAccEffMultDep[selectSp]) {
    hSysErrAccEffMultDep[selectSp]->SetLineWidth(lineWidthStacked);
    hSysErrAccEffMultDep[selectSp]->SetFillColor(kMagenta);
    hSysErrAccEffMultDep[selectSp]->SetLineColor(kBlack);
    hSysErrAccEffMultDep[selectSp]->Draw(drawStr.Data());
    drawStr = "same";
  }
  
  if (hSysErrAccEff10d10e[selectSp]) {
    hSysErrAccEff10d10e[selectSp]->SetLineWidth(lineWidthStacked);
    hSysErrAccEff10d10e[selectSp]->SetFillColor(kCyan);
    hSysErrAccEff10d10e[selectSp]->SetLineColor(kBlack);
    hSysErrAccEff10d10e[selectSp]->Draw(drawStr.Data());
    drawStr = "same";
  }
  
  if (hSysErrAccEff[selectSp]) {
    hSysErrAccEff[selectSp]->SetLineWidth(lineWidthStacked);
    hSysErrAccEff[selectSp]->SetFillColor(2);
    hSysErrAccEff[selectSp]->SetLineColor(kBlack);
    hSysErrAccEff[selectSp]->Draw(drawStr.Data());
    drawStr = "same";
  }
  
  for (Int_t i = 1; i <= histForAxis->GetNbinsX(); i++) {
    Double_t s = 0;
    
    if (hSysErrEff && hSysErrEff[selectSp])
      s += hSysErrEff[selectSp]->GetBinContent(i)*hSysErrEff[selectSp]->GetBinContent(i);
    
    if (hSysErrRes && hSysErrRes[selectSp])
      s += hSysErrRes[selectSp]->GetBinContent(i)*hSysErrRes[selectSp]->GetBinContent(i);
    if (jetPtUp > 0 && hSysErrBbB && hSysErrBbB[selectSp]) // No shape uncertainties for inclusive case
      s += hSysErrBbB[selectSp]->GetBinContent(i)*hSysErrBbB[selectSp]->GetBinContent(i);
    if (hSysErrSec && hSysErrSec[selectSp])
      s += hSysErrSec[selectSp]->GetBinContent(i)*hSysErrSec[selectSp]->GetBinContent(i);
    if (hSysErrSS && hSysErrSS[selectSp])
      s += hSysErrSS[selectSp]->GetBinContent(i)*hSysErrSS[selectSp]->GetBinContent(i);
    if (hSysErrGFL && hSysErrGFL[selectSp])
      s += hSysErrGFL[selectSp]->GetBinContent(i)*hSysErrGFL[selectSp]->GetBinContent(i);
    if (hSysErrMu)
      s += hSysErrMu->GetBinContent(i)*hSysErrMu->GetBinContent(i);
    if (hSysErrEff10d10e && hSysErrEff10d10e[selectSp])
      s += hSysErrEff10d10e[selectSp]->GetBinContent(i)*hSysErrEff10d10e[selectSp]->GetBinContent(i);
    if (hSysErrEffMultDep && hSysErrEffMultDep[selectSp])
      s += hSysErrEffMultDep[selectSp]->GetBinContent(i)*hSysErrEffMultDep[selectSp]->GetBinContent(i);
    if (hSysErrSecMultDep && hSysErrSecMultDep[selectSp])
      s += hSysErrSecMultDep[selectSp]->GetBinContent(i)*hSysErrSecMultDep[selectSp]->GetBinContent(i);
    
    
    Double_t sOther = hSysErrAccTot[selectSp]->GetBinContent(i) * hSysErrAccTot[selectSp]->GetBinContent(i);
    
    if (TMath::Abs(s - sOther) > 1e-10) {
      if (TMath::Abs(sOther) > 1e-10) {
        printf("Mismatch %s, bin centre %f: sqrt(s) %e <-> sqrt(sOther) %e. Should only happen due to zero correction factor.\n\n", 
               AliPID::ParticleShortName(selectSp), hSysErrAccEff[selectSp]->GetXaxis()->GetBinCenter(i), TMath::Sqrt(s), TMath::Sqrt(sOther));
      }
      
      // In any case, just set the error to the visible sum, such that people don't get confused by the plot!
      hSysErrAccTot[selectSp]->SetBinContent(i, TMath::Sqrt(s));
      
      
      /* NOTE: Won't help for corr equals zero, since after multiplying the error with zero, it is zero again....
      // And also the error is relative to its individual correction factor. So, things get confused if it is just applied relative
      // to the total correction factor
      hSysErrAccTot[selectSp]->SetBinContent(i, TMath::Sqrt(s));
      hSysErrTot[selectSp]->SetBinContent(i, TMath::Sqrt(s));
      
      Double_t oldVal = hSysErrTotAbs[selectSp]->GetBinError(i);
      hSysErrTotAbs[selectSp]->SetBinError(i, hSysErrTotAbs[selectSp]->GetBinContent(i) * TMath::Sqrt(s));
      printf("Before: %e -> After: %e\n\n", oldVal, hSysErrTotAbs[selectSp]->GetBinError(i));
      */
    }
  }
  
  
  hSysErrAccTot[selectSp]->SetLineColor(1);
  hSysErrAccTot[selectSp]->SetLineWidth(3);  
  hSysErrAccTot[selectSp]->SetLineStyle(2);  
  hSysErrAccTot[selectSp]->Draw("same");
  
  
  TString prefix = "";
  
  if (!onlyMult && hasNonVanishingSysErrEff) {
    leg->AddEntry(hSysErrAccEff[selectSp],Form("%sEfficiency", prefix.Data()),"F");
    prefix = "+ ";
  }
  if (hSysErrAccEff10d10e[selectSp]) {
    leg->AddEntry(hSysErrAccEff10d10e[selectSp],Form("%sEfficiency (data set)", prefix.Data()),"F");
    prefix = "+ ";
  }
  if (hSysErrAccEffMultDep[selectSp]) {
    leg->AddEntry(hSysErrAccEffMultDep[selectSp],Form("%sEfficiency (mult.)", prefix.Data()),"F");
    prefix = "+ ";
  }
  if (hSysErrAccRes[selectSp]) {
    leg->AddEntry(hSysErrAccRes[selectSp],Form("%sResolution", prefix.Data()),"F");
    prefix = "+ ";
  }
  if (hSysErrAccBbB[selectSp]) {
    leg->AddEntry(hSysErrAccBbB[selectSp],Form("%sShape", prefix.Data()),"F");
  //leg->AddEntry(hSysErrAccEff[selectSp],Form("%s", !hSysErrAccGFL[selectSp] ? "#splitline{efficiency,resolution,shape}{+ GEANT/FLUKA}" 
  //                                                                          : "efficiency,resolution,shape"),"F");
  //leg->AddEntry(hSysErrAccSec[selectSp],Form("%s", !hSysErrAccSS[selectSp] ? "#splitline{+ TODOsecondaries}{+ TODOstrangeness scaling}"
  //                                                                         : "+ secondaries"),"F");
    prefix = "+ ";
  }
  if (hSysErrAccSec[selectSp]) {
    leg->AddEntry(hSysErrAccSec[selectSp],Form("%sSecondaries", prefix.Data()),"F");
    prefix = "+ ";
  }
  if (hSysErrAccSecMultDep[selectSp]) {
    leg->AddEntry(hSysErrAccSecMultDep[selectSp],Form("%sSecondaries (mult.)", prefix.Data()),"F");
    prefix = "+ ";
  }
  if (hSysErrAccSS[selectSp]) {
    leg->AddEntry(hSysErrAccSS[selectSp],Form("%sStrangeness scaling", prefix.Data()),"F");
    prefix = "+ ";
  }
  if (hSysErrAccGFL[selectSp]) {
    leg->AddEntry(hSysErrAccGFL[selectSp],Form("%sGEANT/FLUKA", prefix.Data()),"F");
    prefix = "+ ";
  }
  if (hSysErrAccMu) {
    leg->AddEntry(hSysErrAccMu,Form("%sPion/muon", prefix.Data()),"F");
  }
  leg->AddEntry(hSysErrAccTot[selectSp],"Total (quadr. sum)","l");
  leg->Draw();
  
  ClearTitleFromHistoInCanvas(c1);
  
  // Axis ticks and/or grid may be hidden by filled histogram - draw them on top again
  c1->RedrawAxis(""); // To get the ticks visible again
  c1->RedrawAxis("G"); // To get the grid visible again
  
  if (saveName != "")
    c1->SaveAs(saveName.Data());
  
  return c1;
}

//___________________________________________________________________
void setupYaxisEff(TH1* h, Bool_t isJet)
{
  if (!h)
    return;
  
  h->GetYaxis()->SetTitle("Eff. x Acc. x #it{p}_{T} Res.");//"Efficiency x Acceptance x #it{p}_{T} Resolution");
  //h->GetYaxis()->SetTitleSize(0.04);
  //h->GetYaxis()->SetLabelSize(0.04);
  //h->GetYaxis()->SetTitleOffset(1.2);
  
  if (isJet)
    h->GetYaxis()->SetRangeUser(0., 1.5);
  else
    h->GetYaxis()->SetRangeUser(0., 1.01);
}


//___________________________________________________________________
void setupYaxisEffBbB(TH1* h, Bool_t isJet)
{
  if (!h)
    return;
  
  h->GetYaxis()->SetTitle("Corr. Factor (Eff. x Acc. x #it{p}_{T} Res.)");
  //h->GetYaxis()->SetTitleSize(0.04);
  //h->GetYaxis()->SetLabelSize(0.04);
  //h->GetYaxis()->SetTitleOffset(1.2);
  
  if (isJet)
    h->GetYaxis()->SetRangeUser(0.6, 2.);
  else
    h->GetYaxis()->SetRangeUser(0.99, 2.);
}


//___________________________________________________________________
TCanvas* createCanvasForCorrFactor(const TString name, const TString title, Int_t wtopx = 0, Int_t wtopy = 300, Int_t ww = 900, Int_t wh = 900)
{
  TCanvas* c = new TCanvas(name.Data(), title.Data(), wtopx, wtopy, ww, wh);
  c->SetTopMargin(0.02);
  c->SetLeftMargin(0.14);
  c->SetRightMargin(0.03);
  c->SetBottomMargin(0.14);
  
  return c;
}


//___________________________________________________________________
void setupHistCorrFactor(TH1D* h)
{
  if (!h)
    return;
  
  h->SetStats(kFALSE);
  h->SetLineWidth(2.);
  h->GetXaxis()->SetMoreLogLabels(kTRUE);
  h->GetXaxis()->SetNoExponent(kTRUE);
  
  h->GetXaxis()->SetLabelSize(0.05);
  h->GetXaxis()->SetTitleSize(0.06);
  h->GetXaxis()->SetTitleOffset(1.0);
  
  h->GetYaxis()->SetLabelSize(0.05);
  h->GetYaxis()->SetTitleSize(0.06);
  h->GetYaxis()->SetTitleOffset(1.1);
  
  setAxisTitlesItalic(h);
}


//___________________________________________________________________
TLegend* createLegend(TCanvas* c, TString opt = "l", Double_t x1 = 0.4, Double_t y1 = 0.8, Double_t x2 = 0.6, Double_t y2 = 0.94)
{
  if (!c)
    return 0x0;
  
  TLegend* leg = new TLegend(x1, y1, x2, y2);
  
  for (Int_t i = 0; i < c->GetListOfPrimitives()->GetSize(); i++) {
    if (c->GetListOfPrimitives()->At(i)->InheritsFrom(TH1::Class())) {
      TObject* obj = c->GetListOfPrimitives()->At(i);
      leg->AddEntry(obj, obj->GetTitle(), opt.Data());
    }
  }
  
  leg->SetNColumns(2);
  leg->SetTextSize(0.06);
  leg->SetFillStyle(1001);
  leg->SetFillColor(kWhite);
  leg->SetBorderSize(0);
  leg->SetMargin(0.3);
  
  c->cd();
  leg->Draw();
  
  return leg;
}


//___________________________________________________________________
void cleanUpHistEntriesForJets(TH1D* h, Bool_t isPt, Int_t upperJetPt)
{
  // If tracks are selected via ideal cone instead of track refs, it can happen that there are entries with pT > pTjet.
  // Just remove these entries in the pT histos (are anyway not used, since this bin is not filled in data (similar cut there).
  
  if (!h || !isPt || upperJetPt < 0)
    return;
  
  for (Int_t i = 1; i <= h->GetNbinsX(); i++) {
    if (h->GetXaxis()->GetBinLowEdge(i) >= upperJetPt) {
      h->SetBinContent(i, 0);
      h->SetBinError(i, 0);
    }
  }
}


//___________________________________________________________________
TH1D* convertYieldHistFromEtaToY(TH1D* hYieldEta, Int_t species, Double_t etaAbsCut = 0.9)
{
  if (!hYieldEta)
    return 0x0;
  
  TH1D* hYieldRapidity = new TH1D(*hYieldEta);
  hYieldRapidity->SetName(Form("%s_rapidityConverted", hYieldEta->GetName()));
  TString yAxisTitle = hYieldRapidity->GetYaxis()->GetTitle();
  yAxisTitle = yAxisTitle.ReplaceAll("#eta", "y");
  hYieldRapidity->GetYaxis()->SetTitle(yAxisTitle.Data());
  
  // For conversion from dN/deta to dN/dy, divide by the value of this function:
  // par[0] = particle mass
  // par[1] = eta interval (e.g.. 0.9 => |eta| < 0.9)
  // x = pT;
  
  // Use formula to calculate y(eta). Calculate delta_y / delta_eta = (y(eta_2) - y(eta_1)) / (eta_2 - eta_1).
  // Here: eta_1 = -eta_2 and eta_2 = etaAbsCut > 0.
  
  TF1 fY("fY", "0.5 * TMath::Log((TMath::Sqrt([0]*[0]+x*x*TMath::CosH([1])*TMath::CosH([1]))+x*TMath::SinH([1])) / (TMath::Sqrt([0]*[0]+x*x*TMath::CosH([1])*TMath::CosH([1]))-x*TMath::SinH([1])))",
         0.15, 50.);
  fY.SetParameter(0, AliPID::ParticleMass(species));
  
  for (Int_t i = 1; i <= hYieldRapidity->GetNbinsX(); i++) {
    const Double_t pT = hYieldRapidity->GetXaxis()->GetBinCenter(i);
    
    fY.SetParameter(1, TMath::Abs(etaAbsCut));
    const Double_t y2 = fY.Eval(pT);
    
    fY.SetParameter(1, -TMath::Abs(etaAbsCut));
    const Double_t y1 = fY.Eval(pT);
    
    const Double_t deltaY = y2 - y1;
    const Double_t deltaEta = 2. * TMath::Abs(etaAbsCut);
    const Double_t corrFactor = deltaEta / deltaY;
    
    hYieldRapidity->SetBinContent(i, hYieldRapidity->GetBinContent(i) * corrFactor);
    hYieldRapidity->SetBinError(i, hYieldRapidity->GetBinError(i) * corrFactor);
  }
  
  return hYieldRapidity;
}



//___________________________________________________________________
TH1D* convertToPionRatioHistFromEtaToY(TH1D* hRatioEta, Int_t species, Double_t etaAbsCut = 0.9)
{
  if (!hRatioEta)
    return 0x0;
  
  TH1D* hRatioRapidity = new TH1D(*hRatioEta);
  hRatioRapidity->SetName(Form("%s_rapidityConverted", hRatioEta->GetName()));
  TString yAxisTitle = hRatioRapidity->GetYaxis()->GetTitle();
  yAxisTitle = yAxisTitle.ReplaceAll("#eta", "y");
  hRatioRapidity->GetYaxis()->SetTitle(yAxisTitle.Data());
  
  // For conversion from dN/deta to dN/dy, divide by the value of this function for the species and multiply by that for pions:
  // par[0] = particle mass
  // par[1] = eta interval (e.g.. 0.9 => |eta| < 0.9)
  // x = pT;
  
  // Use formula to calculate y(eta). Calculate delta_y / delta_eta = (y(eta_2) - y(eta_1)) / (eta_2 - eta_1).
  // Here: eta_1 = -eta_2 and eta_2 = etaAbsCut > 0.
  // Note: In the ratio, delta_eta drops out!
  
  TF1 fY("fY", "0.5 * TMath::Log((TMath::Sqrt([0]*[0]+x*x*TMath::CosH([1])*TMath::CosH([1]))+x*TMath::SinH([1])) / (TMath::Sqrt([0]*[0]+x*x*TMath::CosH([1])*TMath::CosH([1]))-x*TMath::SinH([1])))",
         0.15, 50.);
  fY.SetParameter(0, AliPID::ParticleMass(species));
  
  TF1 fYpi("fYpi", "0.5 * TMath::Log((TMath::Sqrt([0]*[0]+x*x*TMath::CosH([1])*TMath::CosH([1]))+x*TMath::SinH([1])) / (TMath::Sqrt([0]*[0]+x*x*TMath::CosH([1])*TMath::CosH([1]))-x*TMath::SinH([1])))",
         0.15, 50.);
  fYpi.SetParameter(0, AliPID::ParticleMass(AliPID::kPion));
  
  for (Int_t i = 1; i <= hRatioRapidity->GetNbinsX(); i++) {
    const Double_t pT = hRatioRapidity->GetXaxis()->GetBinCenter(i);
    
    fY.SetParameter(1, TMath::Abs(etaAbsCut));
    const Double_t y2 = fY.Eval(pT);
    
    fY.SetParameter(1, -TMath::Abs(etaAbsCut));
    const Double_t y1 = fY.Eval(pT);
    
    const Double_t deltaY = y2 - y1;
    
    
    fYpi.SetParameter(1, TMath::Abs(etaAbsCut));
    const Double_t y2pi = fYpi.Eval(pT);
    
    fYpi.SetParameter(1, -TMath::Abs(etaAbsCut));
    const Double_t y1pi = fYpi.Eval(pT);
    
    const Double_t deltaYpi = y2pi - y1pi;
    
    
    const Double_t corrFactor = deltaYpi / deltaY;
    
    hRatioRapidity->SetBinContent(i, hRatioRapidity->GetBinContent(i) * corrFactor);
    hRatioRapidity->SetBinError(i, hRatioRapidity->GetBinError(i) * corrFactor);
  }
  
  return hRatioRapidity;
}


//___________________________________________________________________
const Double_t* getBins(Int_t type, Int_t& nBins)
{
  if (type == 2) {
    nBins = nPtBinsType2;
    
    return binsPtType2;
  }
  
  return 0x0;
}


//___________________________________________________________________
Double_t trackingPtGeantFlukaCorrectionPrMinus(Double_t pTmc)
{
  return (1. - 0.129758 * TMath::Exp(-pTmc * 0.679612));
}


//___________________________________________________________________
Double_t trackingPtGeantFlukaCorrectionKaMinus(Double_t pTmc)
{
  return TMath::Min((0.972865 + 0.0117093 * pTmc), 1.);
}


//___________________________________________________________________
Bool_t geantFlukaCorrection(AliCFContainer* data, Int_t genStepToDownscale)
{
  // Issue: GEANT/FLUKA correction factor is for MC_pT. 
  // Finally the effeciency should be DIVIDED by this correction factor.
  // To include resolution effects, it is therefore the best way to just
  // multiply the generated step with the correction factor.
  
  if (!data) {
    printf("No CFContainer for GEANT/FLUKA correction!\n");
    return kFALSE;
  }
  
  if (iPt < 0 || iMCid < 0 || iCharge < 0) {
    printf("Data axis for GEANT/FLUKA correction not found!\n");
    return kFALSE;
  }
  
  if (!data->GetGrid(genStepToDownscale)) {
    printf("Step for downscaling (GEANT/FLUKA) not found!\n");
    return kFALSE;
  }
  
  const Int_t nDim = data->GetNVar();
  Int_t coord[nDim];
  Double_t binCenterCoord[nDim];
  
  Long64_t nBinsGrid = data->GetGrid(genStepToDownscale)->GetGrid()->GetNbins();
  
  for (Long64_t iBin = 0; iBin < nBinsGrid; iBin++) {
    Double_t binContent = data->GetGrid(genStepToDownscale)->GetGrid()->GetBinContent(iBin, coord);
    Double_t binError  = data->GetGrid(genStepToDownscale)->GetGrid()->GetBinError(iBin);
    
    for (Int_t iDim = 0; iDim < nDim; iDim++) 
      binCenterCoord[iDim] = data->GetBinCenter(iDim, coord[iDim]);

    if (binCenterCoord[iCharge] < 0) {
      Double_t corrFactor = 1.;
      
      if (binCenterCoord[iMCid] - 0.5 == AliPID::kProton) 
        corrFactor = trackingPtGeantFlukaCorrectionPrMinus(binCenterCoord[iPt]);
      else if (binCenterCoord[iMCid] - 0.5 == AliPID::kKaon)
        corrFactor = trackingPtGeantFlukaCorrectionKaMinus(binCenterCoord[iPt]);
      else
        continue;
      
      data->GetGrid(genStepToDownscale)->GetGrid()->SetBinContent(iBin, binContent * corrFactor);
      data->GetGrid(genStepToDownscale)->GetGrid()->SetBinError(iBin, binError * corrFactor);
    }
  }
  
  return kTRUE;
}


//___________________________________________________________________
void applyMultCorrection(TH1D* hist[AliPID::kSPECIES], TH1D* hMultParaSec[AliPID::kSPECIES][numParamsMult], Double_t mult)
{
  TF1 multPara("multPara", "[0]+[1]*exp([2]*x)", 0., 100.);
  for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
    if (!hist[species] || !hMultParaSec[species][0])
      continue;
    
    // Check if all histos are available
    for (Int_t iPara = 0; iPara < numParamsMult; iPara++) {
      if (!hMultParaSec[species][iPara]) {
        printf("Fatal: Histogram for mult correction missing!\n");
        exit(-1);
      }
    }
    
    // Each histogram has the parameters for a fixed pT bin. The multiplicity is fixed from the data sample.
    // Loop over pT bins...
    for (Int_t iX = 1; iX <= hist[species]->GetNbinsX(); iX++) {
      for (Int_t iPara = 0; iPara < numParamsMult; iPara++)
        multPara.SetParameter(iPara, hMultParaSec[species][iPara]->GetBinContent(iX));
      
      // ... and apply the correction
      const Double_t corr = multPara.Eval(mult);
      hist[species]->SetBinContent(iX, corr * hist[species]->GetBinContent(iX));
      hist[species]->SetBinError(iX, corr * hist[species]->GetBinError(iX));
      
      //printf("%f-%f: %e\n", hist[species]->GetXaxis()->GetBinLowEdge(iX), hist[species]->GetXaxis()->GetBinUpEdge(iX), corr);
    }
  }
}


//____________________________________________________________________________________________________________________
void undoJetNormalisationHist2D(TH2* hData, TH2* hNumJets, const Int_t lowerCentralityBinLimit, const Int_t upperCentralityBinLimit)
{
  // Undo normalisation to 1/numJets. NOTE: jetPt binning of hData and hNumJets assumed to be the same!
  
  if (!hData || !hNumJets)
    return;
  
  for (Int_t binJetPt = 0; binJetPt <= hData->GetNbinsY() + 1; binJetPt++) {
    const Double_t numJets = hNumJets->Integral(lowerCentralityBinLimit, upperCentralityBinLimit, binJetPt, binJetPt);
    Bool_t noJets = numJets < 1e-13;
    
    for (Int_t binObs = 0; binObs <= hData->GetNbinsX() + 1; binObs++) {
      if (noJets) 
        continue;

      hData->SetBinContent(binObs, binJetPt, hData->GetBinContent(binObs, binJetPt) * numJets);
      hData->SetBinError(binObs, binJetPt, hData->GetBinError(binObs, binJetPt) * numJets);
    }
  }
}


//___________________________________________________________________
void normaliseHist(TH1* h, Double_t scale = 1.0)
{
  if (h->GetSumw2N() <= 0)
    h->Sumw2();
  
  h->Scale(scale);
  
  for (Int_t i = 1; i <= h->GetNbinsX(); i++) {
    Double_t normFactor = h->GetBinWidth(i);
    h->SetBinContent(i, h->GetBinContent(i) / normFactor);
    h->SetBinError(i, h->GetBinError(i) / normFactor);
  }
}


//___________________________________________________________________
void multiplyHistsDifferentBinning(TH1* h1, const TH1* h2, Double_t c1, Double_t c2, Bool_t ignoreErrorOfSecondHist = kFALSE)
{
   Int_t nbinsx = h1->GetNbinsX();
   Int_t nbinsy = h1->GetNbinsY();
   Int_t nbinsz = h1->GetNbinsZ();
   
   if (h1->GetDimension() < 2)
     nbinsy -= 1;
   if (h1->GetDimension() < 3)
     nbinsz -= 1;
   
   if (h1->GetSumw2N() == 0)
     h1->Sumw2();
   
   Int_t bin, bin2, binx, biny, binz;
   Double_t b1,b2,w,d1,d2;
   d1 = c1*c1;
   d2 = c2*c2;
   
   // Loop over bins (including underflows/overflows)
   for (binz = 0; binz <= nbinsz + 1; binz++) {
      for (biny = 0; biny <= nbinsy + 1; biny++) {
         for (binx = 0; binx <= nbinsx + 1;binx++) {
            bin = binx + (nbinsx + 2) * (biny + (nbinsy + 2) * binz);
            bin2 = h2->FindFixBin(h1->GetXaxis()->GetBinCenter(binx), h1->GetYaxis()->GetBinCenter(biny),
                                  h1->GetZaxis()->GetBinCenter(binz));
            b1  = h1->GetBinContent(bin);
            b2  = h2->GetBinContent(bin2);
            w   = (c1*b1)*(c2*b2);
            h1->SetBinContent(bin, w);
            Double_t e1 = h1->GetBinError(bin);
            Double_t e2 = ignoreErrorOfSecondHist ? 0. : h2->GetBinError(bin2);
            h1->SetBinError(bin, TMath::Sqrt(d1*d2*(e1*e1*b2*b2 + e2*e2*b1*b1)));
         }
      }
   }
}


//___________________________________________________________________
void divideHistsDifferentBinning(TH1* h1, const TH1* h2, Double_t c1, Double_t c2, Bool_t ignoreErrorOfSecondHist = kFALSE)
{
   Int_t nbinsx = h1->GetNbinsX();
   Int_t nbinsy = h1->GetNbinsY();
   Int_t nbinsz = h1->GetNbinsZ();
   
   if (h1->GetDimension() < 2)
     nbinsy -= 1;
   if (h1->GetDimension() < 3)
     nbinsz -= 1;
   
   if (h1->GetSumw2N() == 0)
     h1->Sumw2();
   
   Int_t bin, bin2, binx, biny, binz;
   Double_t b1,b2,w,d1,d2;
   d1 = c1*c1;
   d2 = c2*c2;
   
   // Loop over bins (including underflows/overflows)
   for (binz = 0; binz <= nbinsz + 1; binz++) {
      for (biny = 0; biny <= nbinsy + 1; biny++) {
         for (binx = 0; binx <= nbinsx + 1;binx++) {
            bin = binx + (nbinsx + 2) * (biny + (nbinsy + 2) * binz);
            bin2 = h2->FindFixBin(h1->GetXaxis()->GetBinCenter(binx), h1->GetYaxis()->GetBinCenter(biny),
                                  h1->GetZaxis()->GetBinCenter(binz));
            b1  = h1->GetBinContent(bin);
            b2  = h2->GetBinContent(bin2);
            if (b2)
              w = (c1*b1)/(c2*b2);
            else 
              w = 0;
            h1->SetBinContent(bin, w);
            
            if (!b2) {
              h1->SetBinError(bin, 0);
              continue;
            }
            
            Double_t b22 = b2*b2*d2;
            Double_t e1 = h1->GetBinError(bin);
            Double_t e2 = ignoreErrorOfSecondHist ? 0. : h2->GetBinError(bin2);
            h1->SetBinError(bin, TMath::Sqrt(d1*d2*(e1*e1*b2*b2 + e2*e2*b1*b1)/(b22*b22)));
         }
      }
   }
}


//___________________________________________________________________
TH1D* getRelErrorHist(TH1D* h, TString histName, TString yTitle)
{
  if (!h)
    return 0x0;
  
  TH1D* hRel = new TH1D(*h);
  hRel->Reset();
  hRel->SetMinimum();
  hRel->SetMaximum();
  hRel->GetYaxis()->SetRange(0, -1);
  hRel->SetName(histName.Data());
  hRel->GetYaxis()->SetTitle(yTitle.Data());
      
  for (Int_t i = 1; i <= hRel->GetNbinsX(); i++) {
    const Double_t val = h->GetBinContent(i);
    const Double_t err = h->GetBinError(i);
    hRel->SetBinContent(i, val > 0 ? err / val : 0);
    hRel->SetBinError(i, 0);
  }
  
  return hRel;
}


//___________________________________________________________________
TH1D* getSysErrorHisto(const TH1D* h)
{
  // Take original histogram with statistical errors and clone. Set errors to zero.

  if (!h)
    return 0x0;
  
  TH1D* hSys = new TH1D(*h);
  hSys->SetName(Form("%sSys", h->GetName()));
  for (Int_t i = 0; i <= hSys->GetNbinsX(); i++)
    hSys->SetBinError(i, 0);
  
  return hSys;
}


//___________________________________________________________________
void getToPionRatioMCrelSysError(TH1D* hMCRelSysErrorToPiRatio[], TH1D* hMCRelSysError[], const TString histNameSuffix)
{
  // Calculates the relative sys. errors of the MC corrections for the to-pion ratios.
  // Assume 100% correlation of the correction factors (CFs), i.e. if CF A of pions goes up by dA, so the CF B for another species will with
  // dB, i.e. the signs of dA and dB are expected to be the same.
  // Looking at the ratio R = A/B, it is changed to R' = A/B * ( (1+dA/A) / (1+dB/B) ). As before, take 50% of the maximum deviations
  // of all these correction factors as systematic error:
  // 0.5 * A/B * |max_variations{(1+dA/A)/(1+dB/B)} - min_variations{(1+dA/A)/(1+dB/B)}|.
  // For the relative sys error, R=A/B drops out completely
  
  for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
    if (!hMCRelSysError[species] || species == AliPID::kPion)
      continue;
    
    hMCRelSysErrorToPiRatio[species] = new TH1D(*hMCRelSysError[species]);
    hMCRelSysErrorToPiRatio[species]->Reset();
    hMCRelSysErrorToPiRatio[species]->SetMaximum();
    hMCRelSysErrorToPiRatio[species]->SetMinimum();
    hMCRelSysErrorToPiRatio[species]->SetName(Form("hMCRelSysErrorToPiRatio%s_%s", histNameSuffix.Data(),
                                                   AliPID::ParticleShortName(species)));
    for (Int_t i = 1; i <= hMCRelSysErrorToPiRatio[species]->GetNbinsX(); i++) {
      // R is not needed, since it drops out for the relative error
      const Double_t relErrCorrFactorSpecies = hMCRelSysError[species]->GetBinContent(i); // = dA/A
      const Double_t relErrCorrFactorPion    = hMCRelSysError[AliPID::kPion]->GetBinContent(i); // = dB/B
      
      // Consider all 3 variations: dA = dB = 0; dA, dB > 0; dA, dB < 0
      const Int_t arrSize = 3;
      Double_t arr[arrSize] = { 1., 
                                (1. + relErrCorrFactorPion) > 0 ? (1. + relErrCorrFactorSpecies) / (1. + relErrCorrFactorPion)
                                                                : 1.,
                                (1. - relErrCorrFactorPion) > 0 ? (1. - relErrCorrFactorSpecies) / (1. - relErrCorrFactorPion)
                                                                : 1.};
      Double_t maxEl = TMath::MaxElement(arrSize, arr);
      Double_t minEl = TMath::MinElement(arrSize, arr);
      
      const Double_t relErrCorrFactorToPiRatio = 0.5 * TMath::Abs(maxEl - minEl);
      hMCRelSysErrorToPiRatio[species]->SetBinContent(i, relErrCorrFactorToPiRatio);
      hMCRelSysErrorToPiRatio[species]->SetBinError(i, 0);
    }
    
    hMCRelSysErrorToPiRatio[species]->GetYaxis()->SetRange(-1, 0);
  }
}


//___________________________________________________________________
Bool_t extractMuonCorrectionFactor(AliCFContainer* data, TString obs, TH1D** hPiFracInMuPi, TH1D** hPiFracInMuPiStrangeScale)
{
  // Extract correction factor for muon contamination in pi+mu sample (primaries AND secondaries!)
 
  if (!hPiFracInMuPi || !hPiFracInMuPiStrangeScale)
    return kFALSE;
  
  // Take all = secondaries and primaries
  TH2D* hMCobs = (TH2D*)data->Project(kStepRecWithRecCutsMeasuredObs, iObsAxis, iMCid);
  TH1D* hMuAll = hMCobs->ProjectionX("hMuAll", AliPID::kMuon + 1, AliPID::kMuon + 1, "e");
  TH1D* hPiAll = hMCobs->ProjectionX("hPiAll", AliPID::kPion + 1, AliPID::kPion + 1, "e");
  
  TH1D* hMuPiAll = new TH1D(*hPiAll);
  hMuPiAll->SetName("hMuPiAll");
  hMuPiAll->Add(hMuAll);
  
  *hPiFracInMuPi = new TH1D(*hPiAll);
  // Binomial errors since nom real subset of den
  (*hPiFracInMuPi)->Divide(hPiAll, hMuPiAll, 1., 1., "B"); 
  setupHist(*hPiFracInMuPi, Form("hPiFracInMuPi_%s", obs.Data()), "#pi_{all} / (#pi_{all} + #mu_{all})",
            data->GetVarTitle(iObsAxis), "Correction factor", kRed);
  
  delete hMCobs;
  delete hMuAll;
  delete hPiAll;
  delete hMuPiAll;
  
  
  // Same with strangeness scaling
  TH2D* hMCobsStrangeScale = (TH2D*)data->Project(kStepRecWithRecCutsMeasuredObsStrangenessScaled, iObsAxis, iMCid);
  TH1D* hMuAllStrangeScale = hMCobsStrangeScale->ProjectionX("hMuAllStrangeScale", AliPID::kMuon + 1, AliPID::kMuon + 1, "e");
  TH1D* hPiAllStrangeScale = hMCobsStrangeScale->ProjectionX("hPiAllStrangeScale", AliPID::kPion + 1, AliPID::kPion + 1, "e");
  
  TH1D* hMuPiAllStrangeScale = new TH1D(*hPiAllStrangeScale);
  hMuPiAllStrangeScale->SetName("hMuPiAllStrangeScale");
  hMuPiAllStrangeScale->Add(hMuAllStrangeScale);
  
  *hPiFracInMuPiStrangeScale = new TH1D(*hPiAllStrangeScale);
  // Binomial errors since nom real subset of den
  (*hPiFracInMuPiStrangeScale)->Divide(hPiAllStrangeScale, hMuPiAllStrangeScale, 1., 1., "B"); 
  setupHist(*hPiFracInMuPiStrangeScale, Form("hPiFracInMuPiStrangeScale_%s", obs.Data()), 
            "#pi_{all} / (#pi_{all} + #mu_{all})", data->GetVarTitle(iObsAxis), "Correction factor", kRed);
  
  delete hMCobsStrangeScale;
  delete hMuAllStrangeScale;
  delete hPiAllStrangeScale;
  delete hMuPiAllStrangeScale;
  
  return kTRUE;
}


//___________________________________________________________________
void convertEffToBinByBinCorr(TH1D* hEfficiency[AliPID::kSPECIES], TH1D* hEffBinByBinCorr[AliPID::kSPECIES], Bool_t isJet)
{
  for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
    if (!hEfficiency[species])
      continue;
    
    hEffBinByBinCorr[species] = new TH1D(*hEfficiency[species]);
    TString histName = hEfficiency[species]->GetName();
    histName = histName.ReplaceAll("Efficiency", "EffBinByBinCorr");
    hEffBinByBinCorr[species]->SetName(histName.Data());
    
    for (Int_t i = 1; i <= hEfficiency[species]->GetNbinsX(); i++) {
      const Double_t eff = hEfficiency[species]->GetBinContent(i);
      const Double_t effErr = hEfficiency[species]->GetBinError(i);
      
      const Double_t effBbB = eff > 0. ? 1. / eff : 0.;
      const Double_t effBbBErr = eff > 0. ? effErr / (eff * eff) : 0.;
      
      hEffBinByBinCorr[species]->SetBinContent(i, effBbB);
      hEffBinByBinCorr[species]->SetBinError(i, effBbBErr);
    }
    
    setupYaxisEffBbB(hEffBinByBinCorr[species], isJet);
  }
}


//___________________________________________________________________
void convertMuCorrForToPiRatio(TH1D* hMuCorr, TH1D** hMuCorrToPiRatio)
{
  // For to-pion ratio, the correction factor needs to be inverted
  if (!hMuCorr)
    return;
  
  *hMuCorrToPiRatio = new TH1D(*hMuCorr);
  TString histName = hMuCorr->GetName();
  histName = Form("%s_toPiRatio", histName.Data());
  (*hMuCorrToPiRatio)->SetName(histName.Data());
  
  for (Int_t i = 1; i <= (*hMuCorrToPiRatio)->GetNbinsX(); i++) {
    const Double_t corr = hMuCorr->GetBinContent(i);
    const Double_t corrErr = hMuCorr->GetBinError(i);
    
    const Double_t corrToPiRatio = corr > 0. ? 1. / corr : 0.;
    const Double_t corrTOPiRatioErr = corr > 0. ? corrErr / (corr * corr) : 0.;
    
    (*hMuCorrToPiRatio)->SetBinContent(i, corrToPiRatio);
    (*hMuCorrToPiRatio)->SetBinError(i, corrTOPiRatioErr);
  }
}


//___________________________________________________________________
Bool_t extractEfficiencies(AliCFContainer* data, TFile* saveFile, TString suffixGF, Int_t iObs, Int_t genStepEff,
                           Bool_t restrictJetPtAxis, Double_t actualUpperJetPt, Double_t nJetsGen, Double_t nJetsRec,
                           TH1D* hYield[AliPID::kSPECIES],
                           TH1D* hEffBinByBinCorr[AliPID::kSPECIES], TH1D* hEffBinByBinCorrToPiRatio[AliPID::kSPECIES],
                           TH1D** hEffAll)
{
  // NOTE: Takes and returns the bin-by-bin correction factors for the efficiency, i.e. the INVERTED efficiency!
  // NOTE: Except for the total efficiency!
  
  if (!saveFile) {
    printf("No save file!\n");
    return kFALSE;
  }
  
  TH1D* hEfficiency[AliPID::kSPECIES] = { 0x0, };
  TH1D* hEfficiencyToPiRatio[AliPID::kSPECIES] = { 0x0, };
  
  // Get single track efficiencies (x acceptance), i.e. without the resolution correction (pT (, eta)), all generated observables
  // -> Needed to estimate sys. error of MC correction.
  // Only relevant in case of inclusive!
  
  TCanvas* cEffSingleTrack = 0x0;
  TH1D* hSingleTrackEfficiency[AliPID::kSPECIES] = { 0x0, };
  
  //if (!restrictJetPtAxis)
  {
    // Construct the efficiency grid from the data container 
    // Clone data since the scaling for some reason changes data!
    AliCFContainer* dataForSingleTrackEff = new AliCFContainer(*data);
    dataForSingleTrackEff->SetName(Form("dataForSingleTrackEff%s", suffixGF.Data()));
    AliCFEffGrid* effSingleTrack = new AliCFEffGrid(Form("effSingleTrack%s", suffixGF.Data()), "Efficiency x Acceptance",
                                                    *dataForSingleTrackEff);
    
    effSingleTrack->CalculateEfficiency(kStepRecWithRecCutsPrimaries, genStepEff);
    
    // If the jet axis is restricted (i.e. jet input), scale with the corresponding number of jets.
    // Note: Since this is supposed to be a real scaling, set the error of the scale factor to zero
    // (second element of factor array)
    if (restrictJetPtAxis) {
      Double_t factor_Numerator[2] = { nJetsRec > 0 ? 1. / nJetsRec : 0., 0.  };
      Double_t factor_Denominator[2] = { nJetsGen > 0 ? 1. / nJetsGen : 0., 0.  };
      effSingleTrack->GetNum()->Scale(factor_Numerator);
      effSingleTrack->GetDen()->Scale(factor_Denominator);
    }
    

    // Get the efficiencies vs. pT for each species
    TH2D* hEffSingleTrackID2Pt = (TH2D*)effSingleTrack->Project(iPt, iMCid);
    hEffSingleTrackID2Pt->SetName(Form("hEffSingleTrackID2Pt%s", suffixGF.Data()));
    for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
      hSingleTrackEfficiency[species] = hEffSingleTrackID2Pt->ProjectionX(Form("hSingleTrackEfficiency%s_%s", suffixGF.Data(), 
                                                                              AliPID::ParticleShortName(species)), species + 1,
                                                                              species + 1, "e");
      hSingleTrackEfficiency[species]->SetTitle(Form("%s", AliPID::ParticleLatexName(species)));
      hSingleTrackEfficiency[species]->SetLineColor(hYield[species]->GetLineColor());
      hSingleTrackEfficiency[species]->SetMarkerColor(hYield[species]->GetLineColor());
      hSingleTrackEfficiency[species]->GetXaxis()->SetRangeUser(0.15, 50);
      hSingleTrackEfficiency[species]->GetYaxis()->SetRangeUser(0., 1.01);
      hSingleTrackEfficiency[species]->GetYaxis()->SetTitle("Efficiency x Acceptance");
      setupHistCorrFactor(hSingleTrackEfficiency[species]);
      cleanUpHistEntriesForJets(hSingleTrackEfficiency[species], kTRUE, -1);
    }
    
    delete hEffSingleTrackID2Pt;
    hEffSingleTrackID2Pt = 0x0;
    
    
    cEffSingleTrack = createCanvasForCorrFactor(Form("cEffSingleTrack%s", suffixGF.Data()), "Efficiency x Acceptance for different species");
    cEffSingleTrack->SetGridx(0);
    cEffSingleTrack->SetGridy(1);
    cEffSingleTrack->SetLogx(1);
    
    hSingleTrackEfficiency[0]->Draw("E1");
    
    for (Int_t i = 1; i < AliPID::kSPECIES; i++) {
      if (i == AliPID::kMuon)
        continue;
      hSingleTrackEfficiency[i]->Draw("E1 same");
    }
    createLegend(cEffSingleTrack, "l", 0.4, 0.25, 0.6, 0.39);
    
    ClearTitleFromHistoInCanvas(cEffSingleTrack);
    
    delete dataForSingleTrackEff;
    dataForSingleTrackEff = 0x0;
    
    delete effSingleTrack;
    effSingleTrack = 0x0;
  }
  
  
  // Construct the efficiency grid from the data container 
  // Clone data since the scaling for some reason changes data!
  AliCFContainer* dataForEff = new AliCFContainer(*data);
  dataForEff->SetName(Form("dataForEff%s", suffixGF.Data()));
  AliCFEffGrid* eff = new AliCFEffGrid(Form("eff%s", suffixGF.Data()), "Efficiency x Acceptance x pT Resolution", *dataForEff);
  
  // Either one can take kStepRecWithGenCutsMeasuredObs or, what I prefer, one can take
  // kStepRecWithRecCutsMeasuredObsPrimaries => The difference is only the eta cut, which is on the rec level
  // in the latter case, i.e. one corrects for eta resolution (although the effect should be very small)
  eff->CalculateEfficiency(kStepRecWithRecCutsMeasuredObsPrimaries, genStepEff);
  
  // If the jet axis is restricted (i.e. jet input), scale with the corresponding number of jets.
  // Note: Since this is supposed to be a real scaling, set the error of the scale factor to zero
  // (second element of factor array)
  if (restrictJetPtAxis) {
    Double_t factor_Numerator[2] = { nJetsRec > 0 ? 1. / nJetsRec : 0., 0.  };
    Double_t factor_Denominator[2] = { nJetsGen > 0 ? 1. / nJetsGen : 0., 0.  };
    eff->GetNum()->Scale(factor_Numerator);
    eff->GetDen()->Scale(factor_Denominator);
  }
  
  //The efficiency along obs and eta, and 2-D projection
  TCanvas* cEff = createCanvasForCorrFactor(Form("cEff%s", suffixGF.Data()), "Efficiency x Acceptance x pT Resolution");
  cEff->Divide(2, 1);
  if (iObs == kTrackPt)
    cEff->GetPad(1)->SetLogx(kTRUE);
  cEff->cd(1);
  *hEffAll = (TH1D*)eff->Project(iObsAxis); //the efficiency vs obs
  (*hEffAll)->SetName(Form("hEfficiency%s_all", suffixGF.Data()));
  (*hEffAll)->SetDirectory(0);
  if (iObs == kTrackPt)
    (*hEffAll)->GetXaxis()->SetRangeUser(0.15, 50.);
  setupYaxisEff((*hEffAll), restrictJetPtAxis);
  setupHistCorrFactor(*hEffAll);
  (*hEffAll)->Draw("E1");
  cEff->cd(2);
  TH1D* hEffEta = (TH1D*)eff->Project(iEta); //the efficiency vs eta
  hEffEta->SetName(Form("hEfficiencyEta%s_all", suffixGF.Data()));
  setupYaxisEff(hEffEta, restrictJetPtAxis);
  setupHistCorrFactor(hEffEta);
  hEffEta->Draw("E1");
  TH2D* hEffID2 = (TH2D*)eff->Project(iObsAxis, iMCid);

  // Get the efficiencies vs. obs for each species
  for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
    hEfficiency[species] = hEffID2->ProjectionX(Form("hEfficiency%s_%s", suffixGF.Data(), AliPID::ParticleShortName(species)),
                                                species + 1, species + 1, "e");
    hEfficiency[species]->SetTitle(Form("%s", AliPID::ParticleLatexName(species)));
    hEfficiency[species]->SetLineColor(hYield[species]->GetLineColor());
    hEfficiency[species]->SetMarkerColor(hYield[species]->GetLineColor());
    if (iObs == kTrackPt)
      hEfficiency[species]->GetXaxis()->SetRangeUser(0.15, actualUpperJetPt > 0 ? TMath::Max(30., actualUpperJetPt) : 50.);
    setupHistCorrFactor(hEfficiency[species]);
    setupYaxisEff(hEfficiency[species], restrictJetPtAxis);
    cleanUpHistEntriesForJets(hEfficiency[species], iObs == kTrackPt, actualUpperJetPt);
  }
  
  delete hEffID2;
  hEffID2 = 0x0;
  
  TCanvas* cEff2 = createCanvasForCorrFactor(Form("cEff2%s", suffixGF.Data()), "Efficiency x Acceptance x pT Resolution for different species");
  cEff2->SetGridx(0);
  cEff2->SetGridy(1);
  if (iObs == kTrackPt)
    cEff2->SetLogx(1);
  
  hEfficiency[0]->Draw("E1");
  
  for (Int_t i = 1; i < AliPID::kSPECIES; i++) {
    if (i == AliPID::kMuon)
        continue;
    hEfficiency[i]->Draw("E1 same");
  }
  createLegend(cEff2);
  
  ClearTitleFromHistoInCanvas(cEff2);
  
  delete dataForEff;
  dataForEff = 0x0;
  
  delete eff;
  eff = 0x0;
  
  
  
  // Efficiency correction for to-pi-ratios
  for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
    hEfficiencyToPiRatio[species] = 0x0;
    
    if (species == AliPID::kPion)
      continue; // Do not consider pion-to-pion ratio
      
    hEfficiencyToPiRatio[species] = new TH1D(*hEfficiency[species]);
    hEfficiencyToPiRatio[species]->Reset();
    hEfficiencyToPiRatio[species]->SetName(Form("hEfficiencyToPionRatio%s_%s", suffixGF.Data(), AliPID::ParticleShortName(species)));
    hEfficiencyToPiRatio[species]->SetTitle(Form("%s/#pi", AliPID::ParticleLatexName(species)));
    
    // Samples for different species are independent, so just divide correction factors
    hEfficiencyToPiRatio[species]->Divide(hEfficiency[species], hEfficiency[AliPID::kPion], 1., 1., ""); 
    hEfficiencyToPiRatio[species]->GetYaxis()->SetRangeUser(0., 2.0);
  }
  
  TCanvas* cEffToPiRatio = createCanvasForCorrFactor(Form("cEffToPiRatio%s", suffixGF.Data()),
                                                          "Efficiency x Acceptance x pT Resolution of to-#pi-ratio for different species");
  cEffToPiRatio->SetGridx(0);
  cEffToPiRatio->SetGridy(1);
  if (iObs == kTrackPt)
    cEffToPiRatio->SetLogx(1);
  
  for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
    if (i == AliPID::kPion)
      continue;
    
    if (i == AliPID::kMuon)
        continue;
    
    hEfficiencyToPiRatio[i]->Draw(Form("E1%s", i == 0 ? "" : " same"));
  }
  
  createLegend(cEffToPiRatio, "l", 0.4, 0.8, 0.64, 0.94);
  
  ClearTitleFromHistoInCanvas(cEffToPiRatio);
  
  
  
  // Convert to bin-by-bin correction factors (i.e. invert)
  convertEffToBinByBinCorr(hEfficiency, hEffBinByBinCorr, restrictJetPtAxis);
  
  TCanvas* cEffBbB = createCanvasForCorrFactor(Form("cEffBbB%s", suffixGF.Data()),
                                               "Bin-by-Bin corection for Efficiency x Acceptance x pT Resolution for different species");
  cEffBbB->SetGridx(0);
  cEffBbB->SetGridy(1);
  if (iObs == kTrackPt)
    cEffBbB->SetLogx(1);
  
  hEffBinByBinCorr[0]->Draw("E1");
  
  for (Int_t i = 1; i < AliPID::kSPECIES; i++) {
    if (i == AliPID::kMuon)
        continue;
    
    hEffBinByBinCorr[i]->Draw("E1 same");
  }
  createLegend(cEffBbB, "l", 0.6, 0.8, 0.8, 0.94);
  
  ClearTitleFromHistoInCanvas(cEffBbB);
  
  
  
  convertEffToBinByBinCorr(hEfficiencyToPiRatio, hEffBinByBinCorrToPiRatio, restrictJetPtAxis);

  TCanvas* cEffBbBToPiRatio = createCanvasForCorrFactor(Form("cEffBbBToPiRatio%s", suffixGF.Data()),
                                                        "Bin-by-Bin correction for Efficiency x Acceptance x pT Resolution of to-#pi-ratio for different species");
  cEffBbBToPiRatio->SetGridx(0);
  cEffBbBToPiRatio->SetGridy(1);
  if (iObs == kTrackPt)
    cEffBbBToPiRatio->SetLogx(1);
  
  for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
    if (i == AliPID::kPion)
      continue;
    
    if (i == AliPID::kMuon)
        continue;
    
    hEffBinByBinCorrToPiRatio[i]->Draw(Form("E1%s", i == 0 ? "" : " same"));
  }
  
  createLegend(cEffBbBToPiRatio, "l", 0.6, 0.8, 0.84, 0.94);
  
  ClearTitleFromHistoInCanvas(cEffBbBToPiRatio);
  
  
  // Save to file
  saveFile->cd();
  
  if (cEffSingleTrack)
    cEffSingleTrack->Write();
  
  for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
    if (hSingleTrackEfficiency[i])
      hSingleTrackEfficiency[i]->Write();
  }
  
  if (cEff)
    cEff->Write();
  
  if (cEff2)
    cEff2->Write();
  
  if (cEffToPiRatio)
    cEffToPiRatio->Write();
  
  if (*hEffAll)
    (*hEffAll)->Write();
  
  if (hEffEta)
    hEffEta->Write();
  
  for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
    if (hEfficiency[i])
      hEfficiency[i]->Write();
  }
  
  for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
    if (hEfficiencyToPiRatio[i])
      hEfficiencyToPiRatio[i]->Write();
  }
  
  if (cEffBbB)
    cEffBbB->Write();
  
  if (cEffBbBToPiRatio)
    cEffBbBToPiRatio->Write();
  
  for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
    if (hEffBinByBinCorr[i])
      hEffBinByBinCorr[i]->Write();
  }
  
  for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
    if (hEffBinByBinCorrToPiRatio[i])
      hEffBinByBinCorrToPiRatio[i]->Write();
  }
  
  return kTRUE;
}


//___________________________________________________________________
// Efficiency for inclusive spectra vs. pT (or also jets, but without using the generation)
// E.g. a 'calcEfficiency.C+("finalCuts/MC_pp/7TeV/LHC10f6a/corrected/finalisedSplines/analytical/Jets/nclCut/noMCidForGen/bhess_PID_Jets_efficiency.root", "finalCuts/pp/7TeV/10d_10e_merged.pass2/finalisedSplines/finalMapsAndTail/Jets/nclCut/outputSystematicsTotal_SummedSystematicErrors__2014_02_12.root", "OliversMacros/sysErr/files", kTRUE, kTRUE, kTRUE, 0, -2, -2, -2, -2, -1, -1, 0, -100, 1, 0.9, 0.852, kTRUE)' -b -q
Int_t calcEfficiency(TString pathNameEfficiency, TString pathNameData, TString pathMCsysErrors, 
                     Bool_t correctGeantFluka, Bool_t scaleStrangeness,
                     Bool_t applyMuonCorrection,
                     Int_t chargeMode /*kNegCharge = -1, kAllCharged = 0, kPosCharge = 1*/,
                     Double_t lowerCentralityData /*= -2*/, Double_t upperCentralityData /*= -2*/,
                     Double_t lowerCentrality /*= -2*/, Double_t upperCentrality /*= -2*/,
                     Double_t lowerJetPt /*= -1*/ , Double_t upperJetPt/* = -1*/,
                     Int_t iObs,
                     Double_t constantCorrectionAbovePtThreshold,
                     Int_t rebinEfficiencyObs,
                     Double_t etaAbsCut /*= 0.9*/,
                     Double_t eps_trigger /*= 0.852*/,
                     typeMCSysErrors sysErrorTypeMC,
                     Bool_t normaliseToNInel = kTRUE,
                     Bool_t correctMCID = kFALSE, // Load results from MC ID instead of fit and correct+compare them. NOTE: Do not use the systematics and other histograms in this case. Only the corrected results and the generated MC truth
                     Bool_t individualMultCorr = kFALSE,
                     TString multCorrPathName = "",
                     Bool_t correctEff10d10e = kFALSE,
                     Bool_t isLowestMultBin = kFALSE)
{
  PrintSettingsAxisRangeForMultiplicityAxisForMB();
  
  Bool_t addMCsysErrors = sysErrorTypeMC != kNoErrors;
  
  if (iObs < 0 || iObs >= kNtypes) {
    printf("Unknown observable: %d\n", iObs);
    return -1;
  }
  
  printf("observable: %s\n", obsString[iObs].Data());
  TString pathData = pathNameData;
  pathData.Replace(pathData.Last('/'), pathData.Length(), "");
  
  TString subDir = "";
  if (sysErrorTypeMC == kErrorsOnlyMultDep)
    subDir = "/results_divided_by_MB";
  else if (sysErrorTypeMC == kErrorsForMerging)
    subDir = "/results_for_merging";
  else if (sysErrorTypeMC == kErrorsForMergingOnlyMultDep)
    subDir = "/results_for_merging/results_divided_by_MB";
  
  TString pathSaveData = Form("%s%s", pathData.Data(), subDir.Data());
  
  TFile* fileEff = TFile::Open(pathNameEfficiency.Data());
  if (!fileEff) {
    printf("Failed to open efficiency file \"%s\"\n", pathNameEfficiency.Data());
    return -1;
  }
  
  AliCFContainer* data = (AliCFContainer*)(fileEff->Get("containerEff"));
  if (!data) {
    printf("Failed to load efficiency container!\n");
    return -1;
  }
  
  // For backward compatibility:
  // Check whether "P_{T}" or "p_{T}" is used
  TString momentumString = "p";
  for (Int_t i = 0; i < data->GetNVar(); i++) {
    TString temp = data->GetVarTitle(i);
    if (temp.Contains("P_{")) {
      momentumString = "P";
      break;
    }
    else if (temp.Contains("p_{")) {
      momentumString = "p";
      break;
    }
  }
  
  iPt     = data->GetVar(Form("%s_{T} (GeV/c)", momentumString.Data()));
  iMCid   = data->GetVar("MC ID");
  iEta    = data->GetVar("#eta");
  iCharge = data->GetVar("Charge (e_{0})");
  iMult   = data->GetVar("Centrality Percentile");
  // Will be set later, if jet pT is restricted
  iJetPt  = 0;
  iZ = 0;
  iXi = 0;
  iDistance = 0;
  iJt = 0;
  
  iObsAxis = iPt; // To be set to other values later;
  
  if (lowerJetPt >= 0 && upperJetPt >= 0) {
    iJetPt    = data->GetVar(Form("%s_{T}^{jet} (GeV/c)", momentumString.Data()));
    iZ        = data->GetVar(Form("z = %s_{T}^{track} / %s_{T}^{jet}", momentumString.Data(), momentumString.Data()));
    iXi       = data->GetVar(Form("#xi = ln(%s_{T}^{jet} / %s_{T}^{track})", momentumString.Data(), momentumString.Data()));
    iDistance = data->GetVar("R");
    iJt       = data->GetVar("j_{T} (GeV/c)");
    
    if (iObs == kZ)
      iObsAxis = iZ;
    else if (iObs == kXi)
      iObsAxis = iXi;
    else if (iObs == kDistance)
      iObsAxis = iDistance;
    else if (iObs == kJT)
      iObsAxis = iJt;
  }
  
  TFile* fileData = TFile::Open(pathNameData.Data());
  if (!fileData) {
    printf("Failed to open data file \"%s\"\n", pathNameData.Data());
    return -1;
  }
  
  TH1D* hYield[AliPID::kSPECIES] = { 0x0, };
  TH1D* hYieldSysError[AliPID::kSPECIES] = { 0x0, };
  TGraphAsymmErrors* gYieldSysError[AliPID::kSPECIES] = { 0x0, };
  TH1D* hMCgenPrimYield[AliPID::kSPECIES] = { 0x0, };
  Int_t numMCgenPrimYieldHistsFound = 0;
  
  for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
    TString speciesName = AliPID::ParticleName(species);
    TString firstLetter = speciesName(0);
    firstLetter.ToUpper();
    speciesName.Replace(0, 1, firstLetter.Data());
    TString histName = Form("hYield%ss%s", speciesName.Data(), correctMCID ? "MC" : "");
    hYield[species] = (TH1D*)fileData->Get(histName.Data());
    if (!hYield[species]) {
      printf("Failed to load hist \"%s\"\n", histName.Data());
      return -1;
    }
    hYield[species]->SetFillStyle(0);
    
    TString graphName = Form("systematicErrorYields_%s", AliPID::ParticleName(species));
    gYieldSysError[species] = (TGraphAsymmErrors*)fileData->Get(graphName.Data());
    
    // In case of MC also retrieve the MC truth generated yields
    TString histNameMCgenYields = Form("hMCgenYieldsPrimSpecies_%s", AliPID::ParticleShortName(species));
    hMCgenPrimYield[species] = (TH1D*)fileData->Get(histNameMCgenYields.Data());
    if (hMCgenPrimYield[species])
      numMCgenPrimYieldHistsFound++;
  }
  
  if (numMCgenPrimYieldHistsFound > 0 && numMCgenPrimYieldHistsFound != AliPID::kSPECIES) {
    printf("Error: Unable to retrieve all MC generated prim yield histos! Got %d.\n", numMCgenPrimYieldHistsFound);
    return -1;
  }
  
  TH1D* hRatioToPi[AliPID::kSPECIES] = { 0x0, };
  TH1D* hRatioToPiSysError[AliPID::kSPECIES] = { 0x0, };
  TGraphAsymmErrors* gRatioToPiSysError[AliPID::kSPECIES] = { 0x0, };
  
  for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
    if (species == AliPID::kPion)
      continue;
    TString speciesName = AliPID::ParticleName(species);
    TString firstLetter = speciesName(0);
    firstLetter.ToUpper();
    speciesName.Replace(0, 1, firstLetter.Data());
    TString histName = Form("hRatioToPi%ss%s", speciesName.Data(), correctMCID ? "MC" : "");
    hRatioToPi[species] = (TH1D*)fileData->Get(histName.Data());
    if (!hRatioToPi[species]) {
      printf("Failed to load hist \"%s\"\n", histName.Data());
    }
    else
      hRatioToPi[species]->SetFillStyle(0);
    
    TString graphName = Form("systematicErrorToPiRatio_%s", AliPID::ParticleName(species));
    gRatioToPiSysError[species] = (TGraphAsymmErrors*)fileData->Get(graphName.Data());
  }
  
  
  TCanvas* cFractions = (TCanvas*)fileData->Get("cFractionsWithTotalSystematicError");
  if (!cFractions)
    cFractions = (TCanvas*)fileData->Get("cFractions");
  
  
  // Convert graphs with systematic errors into histograms (assume symmetric error, which is currently true)
  for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
    if (gYieldSysError[species]) {
      hYieldSysError[species] = new TH1D(*hYield[species]);
      hYieldSysError[species]->SetName(Form("%s_sysError", hYieldSysError[species]->GetName()));
      hYieldSysError[species]->SetFillStyle(0);
      
      for (Int_t binX = 1; binX <= hYieldSysError[species]->GetNbinsX(); binX++) {
        hYieldSysError[species]->SetBinContent(binX, gYieldSysError[species]->GetY()[binX - 1]);
        hYieldSysError[species]->SetBinError(binX, gYieldSysError[species]->GetErrorY(binX - 1));
      }
    }
    
    if (gRatioToPiSysError[species]) {
      if (hRatioToPi[species]) {
        hRatioToPiSysError[species] = new TH1D(*hRatioToPi[species]);
        hRatioToPiSysError[species]->SetName(Form("%s_sysError", hRatioToPiSysError[species]->GetName()));
        hRatioToPiSysError[species]->SetFillStyle(0);
        
        for (Int_t binX = 1; binX <= hRatioToPiSysError[species]->GetNbinsX(); binX++) {
          hRatioToPiSysError[species]->SetBinContent(binX, gRatioToPiSysError[species]->GetY()[binX - 1]);
          hRatioToPiSysError[species]->SetBinError(binX, gRatioToPiSysError[species]->GetErrorY(binX - 1));
        }
      }
    }
  }
  
  
  
  if (applyMuonCorrection) {
    // If muons are corrected, they are removed and should also disappear from the generated MC truth
    if (hMCgenPrimYield[AliPID::kMuon]) {
      for (Int_t i = 1; i <= hMCgenPrimYield[AliPID::kMuon]->GetNbinsX(); i++) {
        hMCgenPrimYield[AliPID::kMuon]->SetBinContent(i, 0.);
        hMCgenPrimYield[AliPID::kMuon]->SetBinError(i, 0.);
      }
    }
    
    if (correctMCID) {
      // Muons are not included in the pions in case of taking MC-ID as "data". Add them and set original yield to zero.
      hYield[AliPID::kPion]->Add(hYield[AliPID::kMuon]);

      for (Int_t i = 1; i <= hMCgenPrimYield[AliPID::kMuon]->GetNbinsX(); i++) {
        hYield[AliPID::kMuon]->SetBinContent(i, 0.);
        hYield[AliPID::kMuon]->SetBinError(i, 0.);
      }
      
      // Sys. error meaningless in case of MC-ID. Set everything to zero.
      for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
        if (hYieldSysError[species]) {
          for (Int_t i = 1; i <= hYieldSysError[species]->GetNbinsX(); i++) {
            hYieldSysError[species]->SetBinContent(i, 0.);
            hYieldSysError[species]->SetBinError(i, 0.);
          }
        }
      }
      
      // A little bit more complicated for the to-pion ratios. Here take the new pion+muon yield and calculate the ratios freshly.
      // Just normal division, since counts statistically independent.
      for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
        if (species == AliPID::kPion || !hRatioToPi[species])
          continue;
        
        
        hRatioToPi[species]->Reset();
        hRatioToPi[species]->Divide(hRatioToPi[species] , hRatioToPi[AliPID::kPion]);
      }
      
      
      // Sys. error meaningless in case of MC-ID. Set everything to zero.
      for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
        if (hRatioToPiSysError[species]) {
          for (Int_t i = 1; i <= hRatioToPiSysError[species]->GetNbinsX(); i++) {
            hRatioToPiSysError[species]->SetBinContent(i, 0.);
            hRatioToPiSysError[species]->SetBinError(i, 0.);
          }
        }
      }
    }
  }
  
  // Take binning from pion yield (binning for all species the same) and create a new AliCFContainer with this new binning
  
  // NOTE: Only the binning of the considered observable matters!
  TAxis* axis = 0x0;
  TArrayD* binsNew = 0x0;
  
  if (iObs == kTrackPt) {
    TH1D hDummyPt(*hYield[AliPID::kPion]);
    hDummyPt.SetName("hDummyPt");
    
    if (rebinEfficiencyObs > 1) {
      Int_t nBinsNew = 0;
      const Double_t* newBins = getBins(rebinEfficiencyObs, nBinsNew);
      
      hDummyPt.SetBins(nBinsNew, newBins);
      
      axis = hDummyPt.GetXaxis();
    
      //axis = hDummyPt.Rebin(rebinEfficiencyObs, "", 0)->GetXaxis();
    }
    else
      axis = hDummyPt.GetXaxis();
    
    const TArrayD* binsPtCurrent = axis->GetXbins();
    binsNew = new TArrayD(*binsPtCurrent);
  }
  else {
    TH1D hDummyObs(*hYield[AliPID::kPion]);
    hDummyObs.SetName("hDummyObs");
    
    if (rebinEfficiencyObs > 1)
      axis = hDummyObs.Rebin(rebinEfficiencyObs, "", 0)->GetXaxis();
    else
      axis = hDummyObs.GetXaxis();

    const TArrayD* binsObsCurrent = axis->GetXbins();
    binsNew = new TArrayD(*binsObsCurrent);
  }
  
  
  const Int_t nEffDims = data->GetNVar();
  Int_t nEffBins[nEffDims];
  
  for (Int_t iDim = 0; iDim < nEffDims; iDim++) {
    if ((iDim == iPt && iObs == kTrackPt) ||
        (iDim == iZ && iObs == kZ) ||
        (iDim == iXi && iObs == kXi) ||
        (iDim == iDistance && iObs == kDistance) ||
        (iDim == iJt && iObs == kJT))
      nEffBins[iDim] = axis->GetNbins();
    else 
      nEffBins[iDim] = data->GetNBins(iDim);
  }
  
 
  // Just make one large pT bin above some threshold, if desired
  if (iObs == kTrackPt && binsNew->fN != 0 && constantCorrectionAbovePtThreshold > 0) {
    for (Int_t iBin = 0; iBin <= nEffBins[iPt] + 1; iBin++) {
      // Find the first bin edged really larger than the threshold.
      // If the bin edge before equals the threshold, just set the
      // current bin edge to the right end of the spectrum -> Done.
      // If the bin edge before is different, set the bin edge to the
      // threshold
      if (binsNew->fArray[iBin] > constantCorrectionAbovePtThreshold) {
        if (binsNew->fArray[iBin - 1] == constantCorrectionAbovePtThreshold) {
          binsNew->fArray[iBin] = binsNew->fArray[nEffBins[iPt]];
          nEffBins[iPt] = iBin;
          break;
        }
        else {
          binsNew->fArray[iBin] = constantCorrectionAbovePtThreshold;
        }
      }
    }
  }
  
  //TODO large xi bin BELOW some threshold
  /*
  // Just make one large Xi bin above some threshold, if desired
  if (iObs == kXi && binsNew->fN != 0 && constantCorrectionAboveXiThreshold > 0) {
    for (Int_t iBin = 0; iBin <= nEffBins[iXi] + 1; iBin++) {
      // Find the first bin edged really larger than the threshold.
      // If the bin edge before equals the threshold, just set the
      // current bin edge to the right end of the spectrum -> Done.
      // If the bin edge before is different, set the bin edge to the
      // threshold
      if (binsNew->fArray[iBin] > constantCorrectionAboveXiThreshold) {
        if (binsNew->fArray[iBin - 1] == constantCorrectionAboveXiThreshold) {
          binsNew->fArray[iBin] = binsNew->fArray[nEffBins[iXi]];
          nEffBins[iXi] = iBin;
          break;
        }
        else {
          binsNew->fArray[iBin] = constantCorrectionAboveXiThreshold;
        }
      }
    }
  }*/
  
  AliCFContainer *dataRebinned = new AliCFContainer(Form("%s_rebinned", data->GetName()), Form("%s (rebinned)", data->GetTitle()),
                                                    data->GetNStep(), nEffDims, nEffBins);
  
  for (Int_t iDim = 0; iDim < nEffDims; iDim++) {
    dataRebinned->SetVarTitle(iDim, data->GetVarTitle(iDim));
    
    if ((iDim == iPt && iObs == kTrackPt) ||
        (iDim == iZ && iObs == kZ) ||
        (iDim == iXi && iObs == kXi) ||
        (iDim == iDistance && iObs == kDistance) ||
        (iDim == iJt && iObs == kJT)) {
      if (binsNew->fN == 0)
        dataRebinned->SetBinLimits(iDim, axis->GetXmin(), axis->GetXmax());
      else
        dataRebinned->SetBinLimits(iDim, binsNew->fArray);
    }
    else {
      dataRebinned->SetBinLimits(iDim, data->GetBinLimits(iDim));
    }
  }
  
  for (Int_t iStep = 0; iStep < data->GetNStep(); iStep++)
    dataRebinned->SetStepTitle(iStep, data->GetStepTitle(iStep));
  
  Int_t coord[nEffDims];
  Double_t binCenterCoord[nEffDims];
  
  // Fill content from old grid into the new grid with proper binning
  for (Int_t iStep = 0; iStep < data->GetNStep(); iStep++) {
    Long64_t nBinsGrid = data->GetGrid(iStep)->GetGrid()->GetNbins();
    
    for (Long64_t iBin = 0; iBin <= nBinsGrid + 1; iBin++) {
      Double_t binContent = data->GetGrid(iStep)->GetGrid()->GetBinContent(iBin, coord);
      Double_t binError2  = data->GetGrid(iStep)->GetGrid()->GetBinError2(iBin);
      
      for (Int_t iDim = 0; iDim < nEffDims; iDim++) {
        binCenterCoord[iDim] = data->GetBinCenter(iDim, coord[iDim]);
      }

      Long64_t iBinRebinned = dataRebinned->GetGrid(iStep)->GetGrid()->GetBin(binCenterCoord);
      dataRebinned->GetGrid(iStep)->GetGrid()->AddBinContent(iBinRebinned, binContent);
      dataRebinned->GetGrid(iStep)->GetGrid()->AddBinError2(iBinRebinned, binError2);
    }
  }
  
  // If desired, restrict centrality axis
  Int_t lowerCentralityBinLimit = -1;
  Int_t upperCentralityBinLimit = -2; // Integral(lowerCentBinLimit, uppCentBinLimit) will not be restricted if these values are kept. In particular, under- and overflow bin will be used!
  Bool_t restrictCentralityAxis = kFALSE;
  Double_t actualLowerCentrality = -1.;
  Double_t actualUpperCentrality = -1.;
  
  if (lowerCentrality >= -1 && upperCentrality >= -1) {
    // Add subtract a very small number to avoid problems with values right on the border between two bins
    lowerCentralityBinLimit = dataRebinned->GetAxis(iMult, 0)->FindFixBin(lowerCentrality + 0.001);
    upperCentralityBinLimit = dataRebinned->GetAxis(iMult, 0)->FindFixBin(upperCentrality - 0.001);
    
    // Check if the values look reasonable
    if (lowerCentralityBinLimit <= upperCentralityBinLimit && lowerCentralityBinLimit >= 0
        && upperCentralityBinLimit <= dataRebinned->GetAxis(iMult, 0)->GetNbins() + 1) {
      restrictCentralityAxis = kTRUE;
    }
    else {
      std::cout << std::endl;
      std::cout << "Requested centrality range out of limits or upper and lower limit are switched!" << std::endl;
      return -1;
    }
  }
  
  if (!restrictCentralityAxis) 
    GetAxisRangeForMultiplicityAxisForMB(dataRebinned->GetAxis(iMult, 0), lowerCentralityBinLimit, upperCentralityBinLimit);
  
  dataRebinned->SetRangeUser(iMult, lowerCentralityBinLimit, upperCentralityBinLimit, kTRUE);
  actualLowerCentrality = dataRebinned->GetAxis(iMult, 0)->GetBinLowEdge(dataRebinned->GetAxis(iMult, 0)->GetFirst());
  actualUpperCentrality = dataRebinned->GetAxis(iMult, 0)->GetBinUpEdge(dataRebinned->GetAxis(iMult, 0)->GetLast());
  
  std::cout << "centrality: ";
  if (restrictCentralityAxis)
    std::cout << actualLowerCentrality << " - " << actualUpperCentrality << std::endl;
  else
    std::cout << "MB (" << actualLowerCentrality << " - " << actualUpperCentrality << ")" << std::endl;
  
  const TString centralityString = !restrictCentralityAxis ? "" : Form("_cent_%.4f_%.4f_", actualLowerCentrality, actualUpperCentrality);
  
  
  // If desired, restrict jetPt axis
  Int_t lowerJetPtBinLimit = -1;
  Int_t upperJetPtBinLimit = -2;
  Bool_t restrictJetPtAxis = kFALSE;
  Double_t actualLowerJetPt = -1.;
  Double_t actualUpperJetPt = -1.;
  
  if (lowerJetPt >= 0 && upperJetPt >= 0) {
    // Add subtract a very small number to avoid problems with values right on the border between to bins
    lowerJetPtBinLimit = dataRebinned->GetAxis(iJetPt, 0)->FindFixBin(lowerJetPt + 0.001);
    upperJetPtBinLimit = dataRebinned->GetAxis(iJetPt, 0)->FindFixBin(upperJetPt - 0.001);
    
    // Check if the values look reasonable
    if (lowerJetPtBinLimit <= upperJetPtBinLimit && lowerJetPtBinLimit >= 0 &&
        upperJetPtBinLimit <= dataRebinned->GetAxis(iJetPt, 0)->GetNbins() + 1) {
      actualLowerJetPt = dataRebinned->GetAxis(iJetPt, 0)->GetBinLowEdge(lowerJetPtBinLimit);
      actualUpperJetPt = dataRebinned->GetAxis(iJetPt, 0)->GetBinUpEdge(upperJetPtBinLimit);

      restrictJetPtAxis = kTRUE;
    }
    else {
      std::cout << std::endl;
      std::cout << "Requested jet pT range out of limits or upper and lower limit are switched!" << std::endl;
      return -1;
    }
  }
  
  std::cout << "jet pT: ";
  if (restrictJetPtAxis) {
    std::cout << actualLowerJetPt << " - " << actualUpperJetPt << std::endl;
    dataRebinned->SetRangeUser(iJetPt, lowerJetPtBinLimit, upperJetPtBinLimit, kTRUE);
  }
  else {
    std::cout << "All" << std::endl;
  }
  
  InitialiseRelSysErrOfSecCorr(!restrictJetPtAxis);
  
  const TString jetString = !restrictJetPtAxis ? "_inclusive" : Form("_jetPt_%.0f_%.0f", actualLowerJetPt, actualUpperJetPt);
  
  // If desired, restrict charge axis
  std::cout << "Charge selection (efficiency): ";
  if (chargeMode == kAllCharged)
    std::cout << "All charged particles" << std::endl;
  else if (chargeMode == kNegCharge)
    std::cout << "Negative particles only" << std::endl;
  else if (chargeMode == kPosCharge)
    std::cout << "Positive particles only" << std::endl;
  else {
    std::cout << "Unknown -> ERROR" << std::endl;
    return -1;
  }
  
  const Bool_t restrictCharge = (chargeMode != kAllCharged);
  
  TString chargeString = "";
  if (chargeMode == kPosCharge)
    chargeString = "_posCharge";
  else if (chargeMode == kNegCharge)
    chargeString = "_negCharge";
  
  
  Int_t lowerChargeBinLimit = -1;
  Int_t upperChargeBinLimit = -2;
  Double_t actualLowerCharge = -999;
  Double_t actualUpperCharge = -999;
  
  if (restrictCharge) {
    // Add subtract a very small number to avoid problems with values right on the border between to bins
    if (chargeMode == kNegCharge) {
      lowerChargeBinLimit = dataRebinned->GetAxis(iCharge, 0)->FindFixBin(-1. + 0.001);
      upperChargeBinLimit = dataRebinned->GetAxis(iCharge, 0)->FindFixBin(0. - 0.001);
    }
    else if (chargeMode == kPosCharge) {
      lowerChargeBinLimit = dataRebinned->GetAxis(iCharge, 0)->FindFixBin(0. + 0.001);
      upperChargeBinLimit = dataRebinned->GetAxis(iCharge, 0)->FindFixBin(1. - 0.001);
    }
    
    // Check if the values look reasonable
    if (lowerChargeBinLimit <= upperChargeBinLimit && lowerChargeBinLimit >= 0
        && upperChargeBinLimit <= dataRebinned->GetAxis(iCharge, 0)->GetNbins() + 1) {
      actualLowerCharge = dataRebinned->GetAxis(iCharge, 0)->GetBinLowEdge(lowerChargeBinLimit);
      actualUpperCharge = dataRebinned->GetAxis(iCharge, 0)->GetBinUpEdge(upperChargeBinLimit);
      
      std::cout << "Charge range (efficiency): " << actualLowerCharge << " - " << actualUpperCharge << std::endl;
    }
    else {
      std::cout << std::endl;
      std::cout << "Requested charge range (efficiency) out of limits or upper and lower limit are switched!" << std::endl;
      return -1;
    }
    
    dataRebinned->SetRangeUser(iCharge, lowerChargeBinLimit, upperChargeBinLimit, kTRUE);
    data->SetRangeUser(iCharge, lowerChargeBinLimit, upperChargeBinLimit, kTRUE);
  }
  
  // If desired, lowed files with multiplicity dependence of correction factors (only pT at the moment!)
  TH1D* hMultParaSec[AliPID::kSPECIES][numParamsMult];
  for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
    for (Int_t iPara = 0; iPara < numParamsMult; iPara++)
      hMultParaSec[species][iPara] = 0x0;
  }
  
  Bool_t doMultCorrectionSec = kFALSE;
  Double_t averageMultOfSample = -1;
  
  if (multCorrPathName != "" && iObs == kTrackPt) {
    // Add the corresponding charge string to the file name
    if (chargeMode == kPosCharge)
      multCorrPathName.ReplaceAll(".root", "_posCharge.root");
    else if (chargeMode == kNegCharge)
      multCorrPathName.ReplaceAll(".root", "_negCharge.root");
    
    TFile* fMultCorr = TFile::Open(multCorrPathName.Data(), "READ");
    if (!fMultCorr) {
      std::cout << "Failed to load file with mult corrections: " << multCorrPathName.Data() << "!" << std::endl;
      return -1;
    }
    
    std::cout << "Loading mult corrections from file: " << multCorrPathName.Data() << std::endl;
    
    for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
      Int_t numLoaded = 0;
      
      for (Int_t iPara = 0; iPara < numParamsMult; iPara++) {
        hMultParaSec[species][iPara] = (TH1D*)fMultCorr->Get(Form("h_%s_p%d", AliPID::ParticleShortName(species), iPara));
        if (!hMultParaSec[species][iPara])
          break;
        numLoaded++;
      }
      
      if (numLoaded != 0 && numLoaded != numParamsMult) {
        std::cout << "Error. Found " << numLoaded << " histos for mult correction, but expected " << numParamsMult << "!" << std::endl;
        return -1;
      }
      else if (numLoaded == numParamsMult) {
        std::cout << "Loaded params for " << AliPID::ParticleShortName(species) << std::endl;
        doMultCorrectionSec = kTRUE;
      }
    }
    
    if (doMultCorrectionSec) {
      averageMultOfSample = 0.5 * (lowerCentralityData + upperCentralityData);
      std::cout << "Applying mult correction assuming average mult " << averageMultOfSample << std::endl;
    }
  }
  else
    std::cout << "Mult corrections disabled" << std::endl;
  
  std::cout << "Correct efficiency 10d<->10e is " << (correctEff10d10e ? "enabled" : "disabled") << std::endl;
  
  std::cout << std::endl;
  
  
  // If jet axis is restricted (i.e. jet input), also need num jet histos for proper normalisation.
  // Load numJet histos from bhess_PID*.root file related to efficiency file.
  TH2D* hNjetsGen = 0x0;
  TH2D* hNjetsRec = 0x0;

  if (restrictJetPtAxis) {
    TString pathNameDataMC = pathNameEfficiency;
    pathNameDataMC.ReplaceAll("_efficiency", "");
    
    TFile* fDataMC = TFile::Open(pathNameDataMC.Data());
    if (!fDataMC)  {
      std::cout << std::endl;
      std::cout << "Failed to open file \"" << pathNameDataMC.Data() << "\" to obtain num of rec/gen jets!" << std::endl;
      
      return -1;
    }
    
    TString listName = pathNameDataMC;
    listName.Replace(0, listName.Last('/') + 1, "");
    listName.ReplaceAll(".root", "");
      
    TObjArray* histList = (TObjArray*)(fDataMC->Get(listName.Data()));
    if (!histList) {
      std::cout << std::endl;
      std::cout << "Failed to load list \"" << listName.Data() << "\" to obtain num of rec/gen jets!" << std::endl;
      return -1;
    }
    
    hNjetsGen = (TH2D*)histList->FindObject("fh2FFJetPtGen");
    hNjetsRec = (TH2D*)histList->FindObject("fh2FFJetPtRec");
    
    if (!hNjetsRec || !hNjetsGen) {
      std::cout << "Failed to load number of jets histos!" << std::endl;
      
      
      // For backward compatibility (TODO REMOVE IN FUTURE): Load info from fixed AnalysisResults file (might be wrong, if other
      // period is considered; also: No multiplicity information)
      TString pathEfficiency = pathNameEfficiency;
      pathEfficiency.Replace(pathEfficiency.Last('/'), pathEfficiency.Length(), "");
      TString pathBackward = Form("%s/AnalysisResults.root", pathEfficiency.Data());
      TFile* fBackward = TFile::Open(pathBackward.Data());
      
      TString dirDataInFile = "";
      TDirectory* dirData = fBackward ? (TDirectory*)fBackward->Get(fBackward->GetListOfKeys()->At(0)->GetName()) : 0x0;
    
      TList* list = dirData ? (TList*)dirData->Get(dirData->GetListOfKeys()->At(0)->GetName()) : 0x0;

      TH1D* hFFJetPtRec = list ? (TH1D*)list->FindObject("fh1FFJetPtRecCuts") : 0x0;
      TH1D* hFFJetPtGen = list ? (TH1D*)list->FindObject("fh1FFJetPtGen") : 0x0;
      
      if (hFFJetPtRec && hFFJetPtGen) {
        printf("***WARNING: For backward compatibility, using file \"%s\" to get number of jets. BUT: Might be wrong period and has no mult info!***\n",
          pathBackward.Data());
        
        hNjetsRec = new TH2D("fh2FFJetPtRec", "", 1, -1, 1, dataRebinned->GetAxis(iJetPt, 0)->GetNbins(),
                            dataRebinned->GetAxis(iJetPt, 0)->GetXbins()->GetArray());
        
        for (Int_t iJet = 1; iJet <= hNjetsRec->GetNbinsY(); iJet++) {
          Int_t lowerBin = hFFJetPtRec->FindFixBin(hNjetsRec->GetYaxis()->GetBinLowEdge(iJet) + 1e-3);
          Int_t upperBin = hFFJetPtRec->FindFixBin(hNjetsRec->GetYaxis()->GetBinUpEdge(iJet) - 1e-3);
          hNjetsRec->SetBinContent(1, iJet, hFFJetPtRec->Integral(lowerBin, upperBin));
        }
        
        hNjetsGen = new TH2D("fh2FFJetPtGen", "", 1, -1, 1,  dataRebinned->GetAxis(iJetPt, 0)->GetNbins(),
                            dataRebinned->GetAxis(iJetPt, 0)->GetXbins()->GetArray());
        
        for (Int_t iJet = 1; iJet <= hNjetsGen->GetNbinsY(); iJet++) {
          Int_t lowerBin = hFFJetPtGen->FindFixBin(hNjetsGen->GetYaxis()->GetBinLowEdge(iJet) + 1e-3);
          Int_t upperBin = hFFJetPtGen->FindFixBin(hNjetsGen->GetYaxis()->GetBinUpEdge(iJet) - 1e-3);
          hNjetsGen->SetBinContent(1, iJet, hFFJetPtGen->Integral(lowerBin, upperBin));
        }
      }
      
      if (!hNjetsRec || ! hNjetsGen)
        return -1;
    }
  }

  // For normalisation to number of jets
  // NOTE: These numbers are for the efficiency only! The data will be normalised to its own number!!!
  const Double_t nJetsGen = hNjetsGen ? hNjetsGen->Integral(lowerCentralityBinLimit, upperCentralityBinLimit, lowerJetPtBinLimit,
                                                            upperJetPtBinLimit) : 1.;
  const Double_t nJetsRec = hNjetsRec ? hNjetsRec->Integral(lowerCentralityBinLimit, upperCentralityBinLimit, lowerJetPtBinLimit,
                                                            upperJetPtBinLimit) : 1.;
  
  
  // Save results to file
  TString saveFileName = pathNameData;
  saveFileName.ReplaceAll(Form("%s/", pathData.Data()), "");
  saveFileName.Prepend("output_EfficiencyCorrection_");
  if (individualMultCorr)
    saveFileName.ReplaceAll("__20", "_individualCorrFactor__20");
  saveFileName.ReplaceAll(".root", Form("%s.root", chargeString.Data()));
  
  if (correctMCID) {
    saveFileName.ReplaceAll(".root", "_MCIDcorrected.root");
  }
  
  TString saveFilePathName = Form("%s/%s", pathSaveData.Data(), saveFileName.Data());
  TFile* saveFile = TFile::Open(saveFilePathName.Data(), "RECREATE");
  
  
  
  
  
  // Correction of muon contamination
  TH1D* hPiFracInMuPi = 0x0;
  TH1D* hPiFracInMuPiStrangeScale = 0x0;
  extractMuonCorrectionFactor(dataRebinned, obsString[iObs], &hPiFracInMuPi, &hPiFracInMuPiStrangeScale);
  
  // Secondary correction
  AliCFEffGrid* sec = new AliCFEffGrid("sec", "Secondary Contamination", *dataRebinned);
  AliCFEffGrid* secStrangeScale = new AliCFEffGrid("secStrangeScale", "Secondary Contamination with Strangeness Scaling",
                                                   *dataRebinned);
  
  // Either one can take kStepRecWithGenCutsMeasuredObs or, what I prefer, one can take
  // kStepRecWithRecCutsMeasuredObsPrimaries => The difference is only the eta cut, which is on the rec level
  // in the latter case, i.e. one corrects for eta resolution (although the effect should be very small)
  // => TESTED: There is NO difference (only 1 bin vs. pt shows a deviation by one entry), so this cut has,
  // in principle, more or less no effect
  
  // For testing with data set w/o strangeness secStrangeScale->CalculateEfficiency(kStepRecWithRecCutsMeasuredObsPrimaries, kStepRecWithRecCutsMeasuredObs);
  secStrangeScale->CalculateEfficiency(kStepRecWithRecCutsMeasuredObsPrimaries, kStepRecWithRecCutsMeasuredObsStrangenessScaled);
  sec->CalculateEfficiency(kStepRecWithRecCutsMeasuredObsPrimaries, kStepRecWithRecCutsMeasuredObs);
  
  // QA plots for secondaries
  TCanvas* cSec = createCanvasForCorrFactor("cSec", "Secondary Contamination");
  cSec->Divide(2, 1);
  if (iObs == kTrackPt)
    cSec->GetPad(1)->SetLogx(kTRUE);
  cSec->cd(1);
  TH1D* hSecAll = (TH1D*)sec->Project(iObsAxis); 
  hSecAll->SetName("hSecAll");
  if (iObs == kTrackPt)
    hSecAll->GetXaxis()->SetRangeUser(0.15, 50.);
  hSecAll->GetYaxis()->SetTitle("Primary Fraction");
  setupHistCorrFactor(hSecAll);
  hSecAll->Draw("E1");
  cSec->cd(2);
  TH1D* hEtaSec = (TH1D*)sec->Project(iEta);
  hEtaSec->SetName("hEtaSec");
  hEtaSec->GetYaxis()->SetTitle("Primary Fraction");
  setupHistCorrFactor(hEtaSec);
  hEtaSec->Draw("E1");
  TH2D* hSecID2 = (TH2D*)sec->Project(iObsAxis, iMCid);
  hSecID2->SetName("hSecID2");
  
  // Get the secondary contamination vs. pT for each species
  TH1D* hSec[AliPID::kSPECIES];
  for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
    hSec[species] = hSecID2->ProjectionX(Form("hSec_%s", AliPID::ParticleShortName(species)), species + 1, species + 1, "e");
    hSec[species]->SetTitle(Form("%s", AliPID::ParticleLatexName(species)));
    hSec[species]->SetLineColor(hYield[species]->GetLineColor());
    hSec[species]->SetMarkerColor(hYield[species]->GetLineColor());
    if (iObs == kTrackPt)
      hSec[species]->GetXaxis()->SetRangeUser(0.15, actualUpperJetPt > 0 ? TMath::Max(30., actualUpperJetPt) : 50);
    hSec[species]->GetYaxis()->SetRangeUser(0., 1.01);
    setupHistCorrFactor(hSec[species]);
    cleanUpHistEntriesForJets(hSec[species], iObs == kTrackPt, upperJetPt);
    hSec[species]->GetYaxis()->SetTitle("Primary Fraction");
  }
  
  TCanvas* cSec2 = createCanvasForCorrFactor("cSec2", "Primary fraction for different species");
  cSec2->SetGridx(0);
  cSec2->SetGridy(1);
  if (iObs == kTrackPt)
    cSec2->SetLogx(1);
  
  hSec[0]->DrawCopy("E1");
  
  for (Int_t i = 1; i < AliPID::kSPECIES; i++) {
    if (i == AliPID::kMuon)
        continue;
    
    hSec[i]->DrawCopy("E1 same");
  }
  
  createLegend(cSec2, "l", 0.69, 0.21, 0.89, 0.35);
  
  ClearTitleFromHistoInCanvas(cSec2);
  
  TCanvas* cSec2MultCorr = 0x0;
  if (doMultCorrectionSec) {
    applyMultCorrection(hSec, hMultParaSec, averageMultOfSample);
    cSec2MultCorr = createCanvasForCorrFactor("cSec2_multCorrected", "Primary fraction for different species");
    cSec2MultCorr->SetGridx(0);
    cSec2MultCorr->SetGridy(1);
    if (iObs == kTrackPt)
      cSec2MultCorr->SetLogx(1);
    
    hSec[0]->Draw("E1");
    
    for (Int_t i = 1; i < AliPID::kSPECIES; i++) {
      if (i == AliPID::kMuon)
          continue;
      
      hSec[i]->Draw("E1 same");
    }
    
    createLegend(cSec2MultCorr, "l", 0.69, 0.21, 0.89, 0.35);
  }
  
  // QA plots for secondaries with strangeness scaling
  TCanvas* cSecSS = createCanvasForCorrFactor("cSecSS", "Secondary Contamination (strangeness scaled)");
  cSecSS->Divide(2, 1);
  if (iObs == kTrackPt)
    cSecSS->GetPad(1)->SetLogx(kTRUE);
  cSecSS->cd(1);
  TH1D* hSecSSall = (TH1D*)secStrangeScale->Project(iObsAxis); 
  hSecSSall->SetName("hSecSSall");
  if (iObs == kTrackPt)
    hSecSSall->GetXaxis()->SetRangeUser(0.15, 50.);
  hSecSSall->GetYaxis()->SetTitle("Primary Fraction");
  setupHistCorrFactor(hSecSSall);
  hSecSSall->Draw("E1");
  cSecSS->cd(2);
  TH1D* hEtaSecSS = (TH1D*)secStrangeScale->Project(iEta);
  hEtaSecSS->SetName("hEtaSecSS");
  hEtaSecSS->GetYaxis()->SetTitle("Primary Fraction");
  setupHistCorrFactor(hEtaSecSS);
  hEtaSecSS->Draw("E1");
  TH2D* hSecSSID2 = (TH2D*)secStrangeScale->Project(iObsAxis, iMCid);
  hSecSSID2->SetName("hSecSSID2");
  
  // Get the secondary contamination vs. pT for each species
  TH1D* hSecSS[AliPID::kSPECIES];
  for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
    hSecSS[species] = hSecSSID2->ProjectionX(Form("hSecSS_%s", AliPID::ParticleShortName(species)), species + 1, species + 1, "e");
    hSecSS[species]->SetTitle(Form("%s", AliPID::ParticleLatexName(species)));
    hSecSS[species]->SetLineColor(hYield[species]->GetLineColor());
    hSecSS[species]->SetMarkerColor(hYield[species]->GetLineColor());
    if (iObs == kTrackPt)
      hSecSS[species]->GetXaxis()->SetRangeUser(0.15, actualUpperJetPt > 0 ? TMath::Max(30., actualUpperJetPt) : 50);
    hSecSS[species]->GetYaxis()->SetRangeUser(0., 1.01);
    setupHistCorrFactor(hSecSS[species]);
    cleanUpHistEntriesForJets(hSecSS[species], iObs == kTrackPt, upperJetPt);
    hSecSS[species]->GetYaxis()->SetTitle("Primary Fraction");
  }
  
  TCanvas* cSecSS2 = createCanvasForCorrFactor("cSecSS2", "Primary fraction for different species (strangeness scaled)");
  cSecSS2->SetGridx(0);
  cSecSS2->SetGridy(1);
  if (iObs == kTrackPt)
    cSecSS2->SetLogx(1);
  
  hSecSS[0]->DrawCopy("E1");
  
  for (Int_t i = 1; i < AliPID::kSPECIES; i++) {
    if (i == AliPID::kMuon)
        continue;
    
    hSecSS[i]->DrawCopy("E1 same");
  }
  
  createLegend(cSecSS2, "l", 0.69, 0.21, 0.89, 0.35);
  
  ClearTitleFromHistoInCanvas(cSecSS2);
  
  TCanvas* cSecSS2MultCorr = 0x0;
  if (doMultCorrectionSec) {
    applyMultCorrection(hSecSS, hMultParaSec, averageMultOfSample);
    cSecSS2MultCorr = createCanvasForCorrFactor("cSecSS2_multCorrected", "Primary fraction for different species (strangeness scaled)");
    cSecSS2MultCorr->SetGridx(0);
    cSecSS2MultCorr->SetGridy(1);
    if (iObs == kTrackPt)
      cSecSS2MultCorr->SetLogx(1);
    
    hSecSS[0]->Draw("E1");
    
    for (Int_t i = 1; i < AliPID::kSPECIES; i++) {
      if (i == AliPID::kMuon)
          continue;
      
      hSecSS[i]->Draw("E1 same");
    }
    
    createLegend(cSecSS2MultCorr, "l", 0.69, 0.21, 0.89, 0.35);
  }
  
  
  // Secondary correction for to-pi-ratios
  TH1D* hSecToPiRatio[AliPID::kSPECIES] = { 0x0, };
  TH1D* hSecSSToPiRatio[AliPID::kSPECIES] = { 0x0, };
  
  for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
    if (species == AliPID::kPion)
      continue; // Do not consider pion-to-pion ratio
      
    hSecToPiRatio[species] = new TH1D(*hSec[species]);
    hSecToPiRatio[species]->Reset();
    hSecToPiRatio[species]->SetName(Form("hSecToPionRatio_%s", AliPID::ParticleShortName(species)));
    hSecToPiRatio[species]->SetTitle(Form("%s/#pi", AliPID::ParticleLatexName(species)));
    hSecToPiRatio[species]->SetLineColor(getLineColorAliPID(species));
    hSecToPiRatio[species]->SetMarkerColor(getLineColorAliPID(species));
    hSecToPiRatio[species]->SetMarkerStyle(24);
    hSecToPiRatio[species]->GetYaxis()->SetTitle("Primary Fraction of Ratio");
    
    hSecSSToPiRatio[species] = new TH1D(*hSecToPiRatio[species]);
    hSecSSToPiRatio[species]->SetName(Form("hSecSSToPionRatio_%s", AliPID::ParticleShortName(species)));
    
    // Samples for different species are independent, so just divide correction factors
    hSecToPiRatio[species]->Divide(hSec[species], hSec[AliPID::kPion], 1., 1., ""); 
    hSecToPiRatio[species]->GetYaxis()->SetRangeUser(0., 1.1);
    
    hSecSSToPiRatio[species]->Divide(hSecSS[species], hSecSS[AliPID::kPion], 1., 1., ""); 
    hSecSSToPiRatio[species]->GetYaxis()->SetRangeUser(0., 1.1);
  }
  
  TCanvas* cSecToPiRatio = createCanvasForCorrFactor("cSecToPiRatio", "Primary fraction of to-#pi-ratio for different species");
  cSecToPiRatio->SetGridx(0);
  cSecToPiRatio->SetGridy(1);
  if (iObs == kTrackPt)
    cSecToPiRatio->SetLogx(1);
  
  for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
    if (i == AliPID::kPion)
      continue;
    
    if (i == AliPID::kMuon)
        continue;
    
    hSecToPiRatio[i]->DrawCopy(Form("E1%s", i == 0 ? "" : " same"));
  }
  
  createLegend(cSecToPiRatio, "l", 0.68, 0.21, 0.96, 0.35);
  
  ClearTitleFromHistoInCanvas(cSecToPiRatio);
  
  TCanvas* cSecToPiRatioMultCorr = 0x0;
  if (doMultCorrectionSec) {
    applyMultCorrection(hSecToPiRatio, hMultParaSec, averageMultOfSample);
    cSecToPiRatioMultCorr = createCanvasForCorrFactor("cSecToPiRatioMultCorr", "Primary fraction of to-#pi-ratio for different species");
    cSecToPiRatioMultCorr->SetGridx(0);
    cSecToPiRatioMultCorr->SetGridy(1);
    if (iObs == kTrackPt)
      cSecToPiRatioMultCorr->SetLogx(1);
    
    for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
      if (i == AliPID::kPion)
        continue;
      
      if (i == AliPID::kMuon)
          continue;
      
      hSecToPiRatio[i]->DrawCopy(Form("E1%s", i == 0 ? "" : " same"));
    }
    
    createLegend(cSecToPiRatioMultCorr, "l", 0.68, 0.21, 0.96, 0.35);
    
    ClearTitleFromHistoInCanvas(cSecToPiRatioMultCorr);
  }
  
  TCanvas* cSecSSToPiRatio = createCanvasForCorrFactor("cSecSSToPiRatio",
                                                       "Primary fraction of to-#pi-ratio for different species (strangeness scaled)");
  cSecSSToPiRatio->SetGridx(0);
  cSecSSToPiRatio->SetGridy(1);
  if (iObs == kTrackPt)
    cSecSSToPiRatio->SetLogx(1);
  
  for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
    if (i == AliPID::kPion)
      continue;
    
    if (i == AliPID::kMuon)
        continue;
    
    hSecSSToPiRatio[i]->DrawCopy(Form("E1%s", i == 0 ? "" : " same"));
  }
  
  createLegend(cSecSSToPiRatio, "l", 0.68, 0.21, 0.96, 0.35);
  
  ClearTitleFromHistoInCanvas(cSecSSToPiRatio);
  
  
  TCanvas* cSecSSToPiRatioMultCorr = 0x0;
  if (doMultCorrectionSec) {
    applyMultCorrection(hSecSSToPiRatio, hMultParaSec, averageMultOfSample);
    cSecSSToPiRatioMultCorr = createCanvasForCorrFactor("cSecSSToPiRatioMultCorr",
                                                        "Primary fraction of to-#pi-ratio for different species (strangeness scaled)");
    cSecSSToPiRatioMultCorr->SetGridx(0);
    cSecSSToPiRatioMultCorr->SetGridy(1);
    if (iObs == kTrackPt)
      cSecSSToPiRatioMultCorr->SetLogx(1);
    
    for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
      if (i == AliPID::kPion)
        continue;
      
      if (i == AliPID::kMuon)
          continue;
      
      hSecSSToPiRatio[i]->Draw(Form("E1%s", i == 0 ? "" : " same"));
    }
    
    createLegend(cSecSSToPiRatioMultCorr, "l", 0.68, 0.21, 0.96, 0.35);
    
    ClearTitleFromHistoInCanvas(cSecSSToPiRatioMultCorr);
  }
  
  
  // Efficiency
  const Int_t genStepEff = kStepGenWithGenCuts;
  
  TH1D* hEfficiency[AliPID::kSPECIES] = { 0x0, };
  TH1D* hEfficiencyToPiRatio[AliPID::kSPECIES] = { 0x0, };
  
  TH1D* hEfficiencyNoGF[AliPID::kSPECIES] = { 0x0, };
  TH1D* hEfficiencyToPiRatioNoGF[AliPID::kSPECIES] = { 0x0, };
  
  TH1D* hEfficiencyWithGF[AliPID::kSPECIES] = { 0x0, };
  TH1D* hEfficiencyToPiRatioWithGF[AliPID::kSPECIES] = { 0x0, };
  
  TH1D* hEfficiencyAllNoPIDNoGF = 0x0;
  TH1D* hEfficiencyAllNoPIDWithGF = 0x0;

  // Extract efficiency without GEANT-FLUKA correction
  extractEfficiencies(dataRebinned, saveFile, "_noGF", iObs, genStepEff, restrictJetPtAxis, actualUpperJetPt, nJetsGen, nJetsRec, hYield, 
                      hEfficiencyNoGF, hEfficiencyToPiRatioNoGF, &hEfficiencyAllNoPIDNoGF);

  
  // If desired, apply the GEANT/FLUKA correction first
  // NOTE: This will change dataRebinned! If anything else should also be calculated, it must be done before.
  // Otherwise, the GEANT/FLUKA correction factor for the efficiency will also affect it!
  if (correctGeantFluka) {
    printf("Applying GEANT/FLUKA correction...\n");
    if (!geantFlukaCorrection(dataRebinned, genStepEff)) {
      printf("GEANT/FLUKA correction could not be applied!\n");
      return kFALSE;
    }
    
    // Extract efficiency with GEANT-FLUKA correction
    extractEfficiencies(dataRebinned, saveFile, "_GF", iObs, genStepEff, restrictJetPtAxis, actualUpperJetPt, nJetsGen, nJetsRec, hYield, 
                        hEfficiencyWithGF, hEfficiencyToPiRatioWithGF, &hEfficiencyAllNoPIDWithGF);
    
    for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
      hEfficiency[species] = hEfficiencyWithGF[species];
      hEfficiencyToPiRatio[species] = hEfficiencyToPiRatioWithGF[species];
    }
  }
  else {
    for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
        hEfficiency[species] = hEfficiencyNoGF[species];
        hEfficiencyToPiRatio[species] = hEfficiencyToPiRatioNoGF[species];
      }
  }
  
  if (correctEff10d10e) {
    printf("\n\nCorrecting efficiency 10d/10e. Scaling efficiency with (10d: %f, 10e: %f, 10d/10ed: %f, 10e/10de: %f, eff(10e)-eff(10d): %f) factor %f...\n\n",
           statFor10d, statFor10e, ratioStat10d_10de, ratioStat10e_10de, effRatio10eOver10d, effCorrFactor10d10e);
    
    // Note: the ratio of the efficiencies between 10d and 10e is the same for all species. Therefore, the to-pion ratios are not affected.
    // Also note: the ratio is basically independent of pT and multiplicity.
    
    for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
      if (hEfficiencyWithGF[species])
        hEfficiencyWithGF[species]->Scale(effCorrFactor10d10e);
      if (hEfficiencyNoGF[species])
        hEfficiencyNoGF[species]->Scale(effCorrFactor10d10e);
    }
  }
  
  
  
  
  // Extract the sys errors on the MC corrections and add them accordingly to the total systematic error
  
  //NOTE: For muon correction, strangeness scaling and GEANT-FLUKA (GF) correction, the difference between the results with and without
  // correction are used. This is mathematical equivalent to taking the (1-f_mu_corr), eps_w/_SS - eps_w/o_S and prim_w/_GF - prim_w/o_GF
  // as errors and using error propagation, AS LONG as these factors are multiplied and NOT divided (i.e. one has to take the real bin-by-
  // bin correction factor difference and not those of the efficiencies).
  // To be on the save site, take the spread of the results and assign as asymmetric error.
  
  
  TH1D* hSecSys[AliPID::kSPECIES];
  TH1D* hEfficiencySys[AliPID::kSPECIES];
  
  for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
    hSecSys[species] = getSysErrorHisto(scaleStrangeness ? hSecSS[species] : hSec[species]);
    hEfficiencySys[species] = getSysErrorHisto(correctGeantFluka ? hEfficiency[species] : hEfficiencyNoGF[species]);
  }
  
  TFile* fMCsys_eff = 0x0;
  TFile* fMCsys_res = 0x0;
  TFile* fMCsys_shape = 0x0;

  TH1D* hMCsysErrorMuCorr = 0x0;
  TH1D* hMCRelSysErrorMuCorr = 0x0;
  //TH1D* hMCRelSysErrorEff[AliPID::kSPECIES] = { 0x0, };
  TH1D* hMCRelSysErrorEff_eff[AliPID::kSPECIES] = { 0x0, };
  TH1D* hMCRelSysErrorEff_res[AliPID::kSPECIES] = { 0x0, };
  TH1D* hMCRelSysErrorEff_shape[AliPID::kSPECIES] = { 0x0, };
  TH1D* hMCRelSysErrorEff10d10e[AliPID::kSPECIES] = { 0x0, };
  TH1D* hMCRelSysErrorSec[AliPID::kSPECIES] = { 0x0, };
  TH1D* hMCRelSysErrorSS[AliPID::kSPECIES] = { 0x0, };
  TH1D* hMCRelSysErrorGF[AliPID::kSPECIES] = { 0x0, };
  
  TH1D* hMCRelSysErrorSecMultDep[AliPID::kSPECIES] = { 0x0, };
  TH1D* hMCRelSysErrorEffMultDep[AliPID::kSPECIES] = { 0x0, };
  
  
  // Correction factors and systematics for to-pion ratios
  
  // NOTE: It is equivalent to take the difference of the to-pion CFs w/ and w/o some correction to propagating the errors with Oliver's
  // formula (and also the assumption 100% positive correlation is the same in both cases)
  TH1D* hSecToPiRatioSys[AliPID::kSPECIES];
  TH1D* hEfficiencyToPiRatioSys[AliPID::kSPECIES];
  
  for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
    hSecToPiRatioSys[species] = getSysErrorHisto(scaleStrangeness ? hSecSSToPiRatio[species] : hSecToPiRatio[species]);
    hEfficiencyToPiRatioSys[species] = getSysErrorHisto(hEfficiencyToPiRatio[species]);
  }
  
  
  TH1D* hPiFracInMuPiToPiRatioStrangeScale = 0x0;
  TH1D* hPiFracInMuPiToPiRatio = 0x0;
  
  convertMuCorrForToPiRatio(hPiFracInMuPiStrangeScale, &hPiFracInMuPiToPiRatioStrangeScale);
  convertMuCorrForToPiRatio(hPiFracInMuPi, &hPiFracInMuPiToPiRatio);
  
  TH1D* hMCsysErrorToPiRatioMuCorr = 0x0;
  
  TH1D* hMCRelSysErrorToPiRatioMuCorr = 0x0;
  //TH1D* hMCRelSysErrorToPiRatioEff[AliPID::kSPECIES] = { 0x0, };
  TH1D* hMCRelSysErrorToPiRatioEff_eff[AliPID::kSPECIES] = { 0x0, };
  TH1D* hMCRelSysErrorToPiRatioEff_res[AliPID::kSPECIES] = { 0x0, };
  TH1D* hMCRelSysErrorToPiRatioEff_shape[AliPID::kSPECIES] = { 0x0, };
  TH1D* hMCRelSysErrorToPiRatioSec[AliPID::kSPECIES] = { 0x0, };
  TH1D* hMCRelSysErrorToPiRatioSS[AliPID::kSPECIES] = { 0x0, };
  TH1D* hMCRelSysErrorToPiRatioGF[AliPID::kSPECIES] = { 0x0, };
  
  //Error drops out: TH1D* hMCRelSysErrorToPiRatioEffMultDep[AliPID::kSPECIES] = { 0x0, };
  TH1D* hMCRelSysErrorToPiRatioSecMultDep[AliPID::kSPECIES] = { 0x0, };
  
  if (addMCsysErrors) {
    // Efficiency
    TH1F* hMCsysErrorEffTmp[AliPID::kSPECIES] = { 0x0, };
    TH1F* hMCsysErrorResTmp[AliPID::kSPECIES] = { 0x0, };
    TH1F* hMCsysErrorShapeTmp[AliPID::kSPECIES] = { 0x0, };
    
    if (restrictJetPtAxis) { // In case of jets, load rel. sys. errors from corresponding file
      // Note no error for "noGF" required, since errors not used anyway in that case
      TString pathNameMCeff = Form("%s/outSysErr_eff.root", pathMCsysErrors.Data());
      fMCsys_eff = TFile::Open(pathNameMCeff.Data(), "READ");
      if (!fMCsys_eff) {
        printf("Failed to load file with MC sys errors: %s\n!", pathNameMCeff.Data());
        return -1;
      }
      
      TString pathNameMCres = Form("%s/outSysErr_res.root", pathMCsysErrors.Data());
      fMCsys_res = TFile::Open(pathNameMCres.Data(), "READ");
      if (!fMCsys_res) {
        printf("Failed to load file with MC sys errors: %s\n!", pathNameMCres.Data());
        return -1;
      }
      
      TString pathNameMCshape = Form("%s/outSysErr_BbB.root", pathMCsysErrors.Data());
      fMCsys_shape = TFile::Open(pathNameMCshape.Data(), "READ");
      if (!fMCsys_shape) {
        printf("Failed to load file with MC sys errors: %s\n!", pathNameMCshape.Data());
        return -1;
      }
      
      
      for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
        hMCsysErrorEffTmp[species] = (TH1F*)fMCsys_eff->Get(Form("hSysErrEff%s_%02.0f_%2.0f_%s", obsStringMCsysError[iObs].Data(), lowerJetPt,
                                                                upperJetPt, AliPID::ParticleShortName(species)));
        hMCsysErrorResTmp[species] = (TH1F*)fMCsys_res->Get(Form("hSysErrRes%s_%02.0f_%2.0f_%s", obsStringMCsysError[iObs].Data(), lowerJetPt,
                                                                upperJetPt, AliPID::ParticleShortName(species)));
        hMCsysErrorShapeTmp[species] = (TH1F*)fMCsys_shape->Get(Form("hSysErrBbB%s_%02.0f_%2.0f_%s", obsStringMCsysError[iObs].Data(),
                                                                    lowerJetPt, upperJetPt, AliPID::ParticleShortName(species)));
        if (!hMCsysErrorEffTmp[species] || !hMCsysErrorResTmp[species] || !hMCsysErrorShapeTmp[species]) {
          printf("ERROR: No MC sys error for %s found!\n", AliPID::ParticleShortName(species));
          for (Int_t i = 1; i <= hEfficiencySys[species]->GetNbinsX(); i++)
            hEfficiencySys[species]->SetBinError(i, 0.);
        }
        else {
          // Check, whether binning is consistent
          for (Int_t i = 1; i <= hMCsysErrorEffTmp[species]->GetNbinsX(); i++) {
            if (TMath::Abs(hMCsysErrorEffTmp[species]->GetXaxis()->GetBinCenter(i) -
                          hEfficiencySys[species]->GetXaxis()->GetBinCenter(i)) > 1e-5 ||
              TMath::Abs(hMCsysErrorResTmp[species]->GetXaxis()->GetBinCenter(i) -
                          hEfficiencySys[species]->GetXaxis()->GetBinCenter(i)) > 1e-5 ||
              TMath::Abs(hMCsysErrorShapeTmp[species]->GetXaxis()->GetBinCenter(i) -
                          hEfficiencySys[species]->GetXaxis()->GetBinCenter(i)) > 1e-5) {
              printf("Error: Inconsistent binning MC corr <-> MC sys error!\n");
              return kFALSE;
            }
            
            // Add rel. errors in quadrature
            Double_t relSysErrEff = hMCsysErrorEffTmp[species]->GetBinContent(i);
            Double_t relSysErrRes = hMCsysErrorResTmp[species]->GetBinContent(i);
            Double_t relSysErrShape = hMCsysErrorShapeTmp[species]->GetBinContent(i);
            Double_t relSysErr = TMath::Sqrt(relSysErrEff*relSysErrEff + relSysErrRes*relSysErrRes + relSysErrShape*relSysErrShape);
            
            // Set error to zero in case of mult errors only! Easiest way of doing it consistently in the following...
            // For the merging, only keep the error from the resolution
            if (sysErrorTypeMC == kErrorsOnlyMultDep || sysErrorTypeMC == kErrorsForMergingOnlyMultDep)
              relSysErr = 0.;
            else if (sysErrorTypeMC == kErrorsForMerging)
              relSysErr = relSysErrRes;
            
            
            // Multiply with correction factor to obtain absolute sys error and set it for the efficiency histo
            hEfficiencySys[species]->SetBinError(i, hEfficiencySys[species]->GetBinContent(i) * relSysErr);
            
            hMCRelSysErrorEff_eff[species] = (TH1D*)hMCsysErrorEffTmp[species]->Clone(Form("hMCRelSysErrorEff_eff_%s",
                                                                                          AliPID::ParticleShortName(species)));
            setupHistMCrelSysErr(hMCRelSysErrorEff_eff[species], hEfficiencySys[species]);
            
            hMCRelSysErrorEff_res[species] = (TH1D*)hMCsysErrorResTmp[species]->Clone(Form("hMCRelSysErrorEff_res_%s",
                                                                                          AliPID::ParticleShortName(species)));
            setupHistMCrelSysErr(hMCRelSysErrorEff_res[species], hEfficiencySys[species]);
            
            hMCRelSysErrorEff_shape[species] = (TH1D*)hMCsysErrorShapeTmp[species]->Clone(Form("hMCRelSysErrorEff_shape_%s",
                                                                                              AliPID::ParticleShortName(species)));
            setupHistMCrelSysErr(hMCRelSysErrorEff_shape[species], hEfficiencySys[species]);
          }
          
          //hMCRelSysErrorEff[species] = getRelErrorHist(hEfficiencySys[species], Form("hMCRelSysErrorEff_%s", AliPID::ParticleShortName(species)),
          //                                             "Rel. sys. error (MC corr, eff)");
        }
      }
    }
    else { // In case of inclusive, no error for shape and pT-independent errors for eff and res
      for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
        if (!hEfficiencySys[species])
          continue;
        
        hMCRelSysErrorEff_eff[species] = (TH1D*)hEfficiencySys[species]->Clone(Form("hMCRelSysErrorEff_eff_%s",
                                                                                    AliPID::ParticleShortName(species)));
        hMCRelSysErrorEff_res[species] = (TH1D*)hEfficiencySys[species]->Clone(Form("hMCRelSysErrorEff_res_%s",
                                                                                    AliPID::ParticleShortName(species)));
        
        setupHistMCrelSysErr(hMCRelSysErrorEff_eff[species], hEfficiencySys[species]);
        setupHistMCrelSysErr(hMCRelSysErrorEff_res[species], hEfficiencySys[species]);
        
        for (Int_t i = 1; i <= hEfficiencySys[species]->GetNbinsX(); i++) {
          // Set rel. errors for eff and res (none of these errors for onlyMultDep, only efficiency for the merging)
          hMCRelSysErrorEff_eff[species]->SetBinContent(i, (sysErrorTypeMC == kErrorsOnlyMultDep || sysErrorTypeMC == kErrorsForMerging ||
                                                            sysErrorTypeMC == kErrorsForMergingOnlyMultDep) 
                                                           ? 0. : relMCSysErrEffInclusive);
          hMCRelSysErrorEff_res[species]->SetBinContent(i, (sysErrorTypeMC == kErrorsOnlyMultDep ||
                                                            sysErrorTypeMC == kErrorsForMergingOnlyMultDep)
                                                           ? 0. : relMCSysErrResInclusive);
          
          // Add rel. errors in quadrature
          Double_t relSysErrEff = hMCRelSysErrorEff_eff[species]->GetBinContent(i);
          Double_t relSysErrRes = hMCRelSysErrorEff_res[species]->GetBinContent(i);
          Double_t relSysErr = TMath::Sqrt(relSysErrEff*relSysErrEff + relSysErrRes*relSysErrRes);
          
          // Multiply with correction factor to obtain absolute sys error and set it for the efficiency histo
          hEfficiencySys[species]->SetBinError(i, hEfficiencySys[species]->GetBinContent(i) * relSysErr);
        }
      }
    }
    
    
    // Compute the secondary correction errors. The relative errors are stored hard-coded in an array. Relative means relative error
    // of (1-hSec*), i.e. relative w.r.t. the deviation from unity
    for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
      if (hSecSys[species]) {
        for (Int_t i = 1; i <= hSecSys[species]->GetNbinsX(); i++) {
          const Double_t value = hSecSys[species]->GetBinContent(i);
          Double_t corr = TMath::Abs(1. - value);
          
          // Set error to zero in case of mult errors only! Easiest way of doing it consistently in the following...
          if (sysErrorTypeMC == kErrorsOnlyMultDep || sysErrorTypeMC == kErrorsForMergingOnlyMultDep)
            corr = 0.;
          
          hSecSys[species]->SetBinError(i, corr * secCorrRelSysError[iObs][species] * value);
        }
        
        hMCRelSysErrorSec[species] = getRelErrorHist(hSecSys[species], Form("hMCRelSysErrorSec_%s", AliPID::ParticleShortName(species)),
                                                     "Rel. sys. error (MC corr, sec)");
      }
    }
    
    // Sys error of muon correction
    if (applyMuonCorrection) {
      hMCsysErrorMuCorr = new TH1D(scaleStrangeness ? *hPiFracInMuPiStrangeScale : *hPiFracInMuPi);
      hMCsysErrorMuCorr->SetName("hMCsysErrorMuCorr");
      hMCsysErrorMuCorr->SetTitle("#mu correction");
      
      // The real correction factor is (1-hist) and for this correction, 100% error is assumed. This means
      // the total error is just 1-hist (is the same as comparing the results obtained w/ and w/o this correction,
      // but with all other corrections, since in the relative error, all other factors drop out and then one multiplies
      // again with the correction factor and propagates....)
      
      for (Int_t i = 1; i <= hMCsysErrorMuCorr->GetNbinsX(); i++) {
        // Effect of strangeness scaling on muon correction is negligible. To be on the save side, take the most extreme
        // corr (w/ or w/o scaling) as systematic error.
        const Double_t value = hPiFracInMuPiStrangeScale->GetBinContent(i);
        Double_t corr = TMath::Abs(1. - value);
        
        const Double_t value2 = hPiFracInMuPi->GetBinContent(i);
        Double_t corr2 = TMath::Abs(1. - value2);
        
        // Set error to zero in case of mult errors only! Easiest way of doing it consistently in the following...
        if (sysErrorTypeMC == kErrorsOnlyMultDep || sysErrorTypeMC == kErrorsForMergingOnlyMultDep) {
          corr = 0.;
          corr2 = 0.;
        }
          
        hMCsysErrorMuCorr->SetBinError(i, TMath::Max(corr, corr2));
      }
      
      hMCRelSysErrorMuCorr = getRelErrorHist(hMCsysErrorMuCorr, "hMCRelSysErrorMuCorr", "Rel. sys. error (MC corr, #mu)");
    }
    
    
    // -------------------------------------------------------------------------------------------------------
    // Before continuing, calculate the errors for the to-pion ratios. This is necessary because GF and SS will
    // be added in the following and they should be decoupled first!
    
    // Extract rel. sys. errors separately for eff, res and shape
    getToPionRatioMCrelSysError(hMCRelSysErrorToPiRatioEff_eff, hMCRelSysErrorEff_eff, "Eff_eff");
    getToPionRatioMCrelSysError(hMCRelSysErrorToPiRatioEff_res, hMCRelSysErrorEff_res, "Eff_res");
    getToPionRatioMCrelSysError(hMCRelSysErrorToPiRatioEff_shape, hMCRelSysErrorEff_shape, "Eff_shape");
    
    // Sum them up in quadrature and obtain the total systematic error from these sources
    for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
      if (!hEfficiencyToPiRatioSys[species])
        continue;
      
      // NOTE: The efficiency histograms contain now the bin-by-bin correction factor, i.e. are already inverted!
      for (Int_t i = 1; i <= hEfficiencyToPiRatioSys[species]->GetNbinsX(); i++) {
        const Double_t corrFactorToPiRatio = hEfficiencyToPiRatioSys[species]->GetBinContent(i);
        
        const Double_t relErrEff   = hMCRelSysErrorToPiRatioEff_eff[species]->GetBinContent(i);
        const Double_t relErrRes   = hMCRelSysErrorToPiRatioEff_res[species]->GetBinContent(i);
        const Double_t relErrShape = hMCRelSysErrorToPiRatioEff_shape[species] ? hMCRelSysErrorToPiRatioEff_shape[species]->GetBinContent(i)
                                                                               : 0.;
        
        const Double_t relErrTot = TMath::Sqrt(relErrEff*relErrEff + relErrRes*relErrRes + relErrShape*relErrShape);
        Double_t errTot = corrFactorToPiRatio * relErrTot;
        
        // Set error to zero in case of mult errors only! Easiest way of doing it consistently in the following...
        // For the merging, only take the error from the resolution
        if (sysErrorTypeMC == kErrorsOnlyMultDep || sysErrorTypeMC == kErrorsForMergingOnlyMultDep)
          errTot = 0.;
        else if (sysErrorTypeMC == kErrorsForMerging)
          errTot = corrFactorToPiRatio * relErrRes;
         
          
        hEfficiencyToPiRatioSys[species]->SetBinError(i, errTot);
      }
      
      //hMCRelSysErrorToPiRatioEff[species] = getRelErrorHist(hEfficiencyToPiRatioSys[species],
      //                                                      Form("hMCRelSysErrorToPiRatioEff_%s", AliPID::ParticleShortName(species)),
      //                                                      "Rel. sys. error (MC corr, eff)");
    }
    
    
    /*OLD, not separate eff,res,shape
    // Assume 100% correlation of the correction factors (CFs), i.e. if CF A of pions goes up by dA, so the CF B for another species will with
    // dB, i.e. the signs of dA and dB are expected to be the same.
    // Looking at the ratio R = A/B, it is changed to R' = A/B * ( (1+dA/A) / (1+dB/B) ). As before, take 50% of the maximum deviations
    // of all these correction factors as systematic error:
    // 0.5 * A/B * |max_variations{(1+dA/A)/(1+dB/B)} - min_variations{(1+dA/A)/(1+dB/B)}|
    
    for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
      if (!hEfficiencyToPiRatioSys[species])
        continue;
      
      // NOTE: The efficiency histograms contain now the bin-by-bin correction factor, i.e. are already inverted!
      // Anyway, this doesn't change the error formula....
      for (Int_t i = 1; i <= hEfficiencyToPiRatioSys[species]->GetNbinsX(); i++) {
        const Double_t corrFactorToPiRatio = hEfficiencyToPiRatioSys[species]->GetBinContent(i); // = R = A/B
        
        const Double_t corrFactorSpecies = hEfficiencySys[species]->GetBinContent(i); // = A
        const Double_t relErrCorrFactorSpecies = corrFactorSpecies > 0 ? hEfficiencySys[species]->GetBinError(i) / corrFactorSpecies : 0; // = dA/A
        
        const Double_t corrFactorPion = hEfficiencySys[AliPID::kPion]->GetBinContent(i); // = B
        const Double_t relErrCorrFactorPion = corrFactorPion > 0 ? hEfficiencySys[AliPID::kPion]->GetBinError(i) / corrFactorPion : 0; // = dB/B
        
        
        // Consider all 3 variations: dA = dB = 0; dA, dB > 0; dA, dB < 0
        const Int_t arrSize = 3;
        Double_t arr[arrSize] = { 1., 
                                  (1. + relErrCorrFactorPion) > 0 ? (1. + relErrCorrFactorSpecies) / (1. + relErrCorrFactorPion)
                                                                  : 1.,
                                  (1. - relErrCorrFactorPion) > 0 ? (1. - relErrCorrFactorSpecies) / (1. - relErrCorrFactorPion)
                                                                  : 1.};
        Double_t maxEl = TMath::MaxElement(arrSize, arr);
        Double_t minEl = TMath::MinElement(arrSize, arr);
        
        const Double_t errCorrFactorToPiRatio = 0.5 * corrFactorToPiRatio * TMath::Abs(maxEl - minEl);
        hEfficiencyToPiRatioSys[species]->SetBinError(i, errCorrFactorToPiRatio);
      }
      
      hMCRelSysErrorToPiRatioEff[species] = getRelErrorHist(hEfficiencyToPiRatioSys[species],
                                                            Form("hMCRelSysErrorToPiRatioEff_%s", AliPID::ParticleShortName(species)),
                                                            "Rel. sys. error (MC corr, eff)");
    }*/
    
    
    
    // Systematics of the secondary correction:
    // Assume independent sys errors and just propagate the individual sys errors of the yields to the to-pion ratios.
    for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
      if (hSecToPiRatioSys[species]) {
        for (Int_t i = 1; i <= hSecToPiRatioSys[species]->GetNbinsX(); i++) {
          const Double_t corrFactorToPiRatio = hSecToPiRatioSys[species]->GetBinContent(i);
        
          const Double_t corrFactorSpecies = hSecSys[species]->GetBinContent(i);
          const Double_t relErrCorrFactorSpecies = corrFactorSpecies > 0 ? hSecSys[species]->GetBinError(i) / corrFactorSpecies : 0;
          
          const Double_t corrFactorPion = hSecSys[AliPID::kPion]->GetBinContent(i);
          const Double_t relErrCorrFactorPion = corrFactorPion > 0 ? hSecSys[AliPID::kPion]->GetBinError(i) / corrFactorPion : 0;
          
          Double_t errCorrFactorToPiRatio = corrFactorToPiRatio * TMath::Sqrt(relErrCorrFactorSpecies * relErrCorrFactorSpecies +
                                                                              relErrCorrFactorPion * relErrCorrFactorPion);
          
          // Set error to zero in case of mult errors only! Easiest way of doing it consistently in the following...
          if (sysErrorTypeMC == kErrorsOnlyMultDep || sysErrorTypeMC == kErrorsForMergingOnlyMultDep)
            errCorrFactorToPiRatio = 0.;
          
          hSecToPiRatioSys[species]->SetBinError(i, errCorrFactorToPiRatio);
        }
        
        hMCRelSysErrorToPiRatioSec[species] = getRelErrorHist(hSecToPiRatioSys[species],
                                                              Form("hMCRelSysErrorToPiRatioSec_%s", AliPID::ParticleShortName(species)),
                                                              "Rel. sys. error (MC corr, sec)");
      }
    }
    
    /*OLD: Like Oliver is doing it....
    // NOTE: Think about this:
    // Assume 100% correlation of the correction factors (CFs), i.e. if CF A of pions goes up by dA, so the CF B for another species will with
    // dB, i.e. the signs of dA and dB are expected to be the same.
    // Looking at the ratio R = A/B, it is changed to R' = A/B * ( (1+dA/A) / (1+dB/B) ).
    // The error is defined as the difference R' - R, i.e. the relative error becomes ( (1+dA/A) / (1+dB/B) - 1 ).
    
    for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
      if (!hEfficiencyToPiRatioSys[species])
        continue;
      
      // NOTE: The efficiency histograms contain now the bin-by-bin correction factor, i.e. are already inverted!
      // Anyway, this doesn't change the error formula....
      for (Int_t i = 1; i <= hEfficiencyToPiRatioSys[species]->GetNbinsX(); i++) {
        const Double_t corrFactorToPiRatio = hEfficiencyToPiRatioSys[species]->GetBinContent(i); // = R = A/B
        
        const Double_t corrFactorSpecies = hEfficiencySys[species]->GetBinContent(i); // = A
        const Double_t relErrCorrFactorSpecies = corrFactorSpecies > 0 ? hEfficiencySys[species]->GetBinError(i) / corrFactorSpecies : 0; // = dA/A
        
        const Double_t corrFactorPion = hEfficiencySys[AliPID::kPion]->GetBinContent(i); // = B
        const Double_t relErrCorrFactorPion = corrFactorPion > 0 ? hEfficiencySys[AliPID::kPion]->GetBinError(i) / corrFactorPion : 0; // = dB/B
        
        const Double_t relErrCorrFactorToPiRatio =  (1. + relErrCorrFactorPion) > 0
                                                      ? (1. + relErrCorrFactorSpecies) / (1. + relErrCorrFactorPion) - 1;
                                                      : 0,
        
        const Double_t errCorrFactorToPiRatio = corrFactorToPiRatio * relErrCorrFactorToPiRatio;
        hEfficiencyToPiRatioSys[species]->SetBinError(i, errCorrFactorToPiRatio);
      }
      
      hMCRelSysErrorToPiRatioEff[species] = getRelErrorHist(hEfficiencyToPiRatioSys[species],
                                                            Form("hMCRelSysErrorToPiRatioEff_%s", AliPID::ParticleShortName(species)),
                                                            "Rel. sys. error (MC corr, eff)");
    }
    */
    /*Even OLDER
    // No dedicated errors for to-pion ratios from Oliver's analysis. Not so clear, how much one double-counts, if independent errors
    // on the yields are assumed. But, on the other hand, varying both, the pion and species efficiency/resolution etc. in the same direction,
    // could mean that one misses part of the error (e.g. resolution is dEdx dependent, maybe resolution underestimated for pions,
    // but overestimated for protons, etc.). 
    // Therefore, it is most convenient and conservative to assume independent sys errors and
    // just propagate the individual sys errors of the yields to the to-pion ratios
    
    for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
      if (!hEfficiencyToPiRatioSys[species])
        continue;
      
      // NOTE: The efficiency histograms contain now the bin-by-bin correction factor, i.e. are already inverted!
      // Anyway, this doesn't change the error formula....
      for (Int_t i = 1; i <= hEfficiencyToPiRatioSys[species]->GetNbinsX(); i++) {
        const Double_t corrFactorToPiRatio = hEfficiencyToPiRatioSys[species]->GetBinContent(i);
        
        const Double_t corrFactorSpecies = hEfficiencySys[species]->GetBinContent(i);
        const Double_t relErrCorrFactorSpecies = corrFactorSpecies > 0 ? hEfficiencySys[species]->GetBinError(i) / corrFactorSpecies : 0;
        
        const Double_t corrFactorPion = hEfficiencySys[AliPID::kPion]->GetBinContent(i);
        const Double_t relErrCorrFactorPion = corrFactorPion > 0 ? hEfficiencySys[AliPID::kPion]->GetBinError(i) / corrFactorPion : 0;
        
        const Double_t errCorrFactorToPiRatio = corrFactorToPiRatio * TMath::Sqrt(relErrCorrFactorSpecies * relErrCorrFactorSpecies +
                                                                                  relErrCorrFactorPion * relErrCorrFactorPion);
        hEfficiencyToPiRatioSys[species]->SetBinError(i, errCorrFactorToPiRatio);
      }
    }*/
    
    // Sys error of muon correction
    if (applyMuonCorrection) {
      hMCsysErrorToPiRatioMuCorr = new TH1D(scaleStrangeness ? *hPiFracInMuPiToPiRatioStrangeScale : *hPiFracInMuPiToPiRatio);
      hMCsysErrorToPiRatioMuCorr->SetName("hMCsysErrorToPiRatioMuCorr");
      hMCsysErrorToPiRatioMuCorr->SetTitle("#mu correction");
      
      // The real correction factor is (1-hist) and for this correction, 100% error is assumed. This means
      // the total error is just 1-hist (is the same as comparing the results obtained w/ and w/o this correction,
      // but with all other corrections, since in the relative error, all other factors drop out and then one multiplies
      // again with the correction factor and propagates....)
      
      // Again: Equivalent to propagating with as Oliver is doing it....
      
      for (Int_t i = 1; i <= hMCsysErrorToPiRatioMuCorr->GetNbinsX(); i++) {
        // Effect of strangeness scaling on muon correction is negligible. To be on the save side, take the most extreme
        // corr (w/ or w/o scaling) as systematic error.
        const Double_t value = hPiFracInMuPiToPiRatioStrangeScale->GetBinContent(i);
        Double_t corr = TMath::Abs(1. - value);
        
        const Double_t value2 = hPiFracInMuPiToPiRatio->GetBinContent(i);
        Double_t corr2 = TMath::Abs(1. - value2);
        
        // Set error to zero in case of mult errors only! Easiest way of doing it consistently in the following...
        if (sysErrorTypeMC == kErrorsOnlyMultDep || sysErrorTypeMC == kErrorsForMergingOnlyMultDep) {
          corr = 0.;
          corr2 = 0.;
        }
      
        hMCsysErrorToPiRatioMuCorr->SetBinError(i, TMath::Max(corr, corr2));
      }
      
      hMCRelSysErrorToPiRatioMuCorr = getRelErrorHist(hMCsysErrorToPiRatioMuCorr, "hMCRelSysErrorToPiRatioMuCorr",
                                                      "Rel. sys. error (MC corr, #mu)");
    }
    
    
    // -------------------------------------------------------------------------------------------------------
    // Continue with the yield systematic error estimation
    
    
    // Sys error of strangeness scaling.
    // Since it is an multiplicative factor, taking the difference w/ and w/o that correction is the same as directly taking
    // the difference of the correction factors for the error.
    // If strangeness scaling is not used, don't add corresponding error
    if (scaleStrangeness &&
        sysErrorTypeMC != kErrorsOnlyMultDep && sysErrorTypeMC != kErrorsForMergingOnlyMultDep) {
      for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
        if (!hSecSys[species])
          continue;
        
        TH1D* hDummy = new TH1D(*hSecSys[species]);
        hDummy->SetDirectory(0);
        
        for (Int_t i = 1; i <= hSecSys[species]->GetNbinsX(); i++) {
          // Add the error in quadrature to the sys error
          const Double_t oldSysError = hSecSys[species]->GetBinError(i);
          const Double_t newSysError = TMath::Abs(hSecSS[species]->GetBinContent(i) - hSec[species]->GetBinContent(i));
          
          hSecSys[species]->SetBinError(i, TMath::Sqrt(oldSysError*oldSysError + newSysError*newSysError));
          hDummy->SetBinError(i, newSysError);
        }
        
        hMCRelSysErrorSS[species] = getRelErrorHist(hDummy, Form("hMCRelSysErrorSS_%s", AliPID::ParticleShortName(species)),
                                                    "Rel. sys. error (MC corr, SS)");
        delete hDummy;
      }
    }
    
    
    // Sys error of GEANT-FLUKA correction.
    // Similarly as for strangeness scaling....
    
    if (correctGeantFluka &&
        sysErrorTypeMC != kErrorsOnlyMultDep && sysErrorTypeMC != kErrorsForMerging && sysErrorTypeMC != kErrorsForMergingOnlyMultDep) {
      // Although only K and p are affected, just loop over all species...
      for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
        if (!hEfficiencySys[species])
          continue;
        
        TH1D* hDummy = new TH1D(*hEfficiencySys[species]);
        hDummy->SetDirectory(0);
        
        for (Int_t i = 1; i <= hEfficiencySys[species]->GetNbinsX(); i++) {
          // Add the error in quadrature to the sys error
          const Double_t oldSysError = hEfficiencySys[species]->GetBinError(i);
          const Double_t newSysError = TMath::Abs(hEfficiency[species]->GetBinContent(i) - hEfficiencyNoGF[species]->GetBinContent(i));
          
          hEfficiencySys[species]->SetBinError(i, TMath::Sqrt(oldSysError*oldSysError + newSysError*newSysError));
          hDummy->SetBinError(i, newSysError);
        }
        
        hMCRelSysErrorGF[species] = getRelErrorHist(hDummy, Form("hMCRelSysErrorGF%s", AliPID::ParticleShortName(species)),
                                                    "Rel. sys. error (MC corr, GF)");
        delete hDummy;
      }
    }
    
    // Sys. error of the 10d/10e efficiency correction (not necessary for the to-pion ratios)
    if (correctEff10d10e &&
        sysErrorTypeMC != kErrorsOnlyMultDep && sysErrorTypeMC != kErrorsForMergingOnlyMultDep) {
      for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
        if (!hEfficiencySys[species])
          continue;
        
        TH1D* hDummy = new TH1D(*hEfficiencySys[species]);
        
        for (Int_t i = 1; i <= hEfficiencySys[species]->GetNbinsX(); i++) {
            // Add the error in quadrature to the sys error
            const Double_t oldSysError = hEfficiencySys[species]->GetBinError(i);
            const Double_t newSysError = hDummy->GetBinContent(i) * TMath::Abs(1. - effCorrFactor10d10e);
            
            hEfficiencySys[species]->SetBinError(i, TMath::Sqrt(oldSysError*oldSysError + newSysError*newSysError));
            hDummy->SetBinError(i, newSysError);
        }
        
        hMCRelSysErrorEff10d10e[species] = getRelErrorHist(hDummy,
                                                          Form("hMCRelSysErrorEff10d10e_%s", AliPID::ParticleShortName(species)),
                                                          "Rel. sys. error (Eff. 10d<->10e)");
        delete hDummy;
      }
    }
    
    
    // Multiplicity dependencies, efficiency
    if (sysErrorTypeMC != kErrorsWithoutMultDep && sysErrorTypeMC != kErrorsForMerging && sysErrorTypeMC != kErrorsForMergingOnlyMultDep) {
      for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
        if (!hEfficiencySys[species])
          continue;
        
        TH1D* hDummy = new TH1D(*hEfficiencySys[species]);
        
        for (Int_t i = 1; i <= hEfficiencySys[species]->GetNbinsX(); i++) {
            // Add the error in quadrature to the sys error
            const Double_t oldSysError = hEfficiencySys[species]->GetBinError(i);
            const Double_t newSysError = hDummy->GetBinContent(i) * (isLowestMultBin ? sysErrMultDepEffCorrLowestMult : sysErrMultDepEffCorr);
            
            hEfficiencySys[species]->SetBinError(i, TMath::Sqrt(oldSysError*oldSysError + newSysError*newSysError));
            hDummy->SetBinError(i, newSysError);
        }
        
        hMCRelSysErrorEffMultDep[species] = getRelErrorHist(hDummy,
                                                            Form("hMCRelSysErrorEffMultDep_%s", AliPID::ParticleShortName(species)),
                                                            "Rel. sys. error (mult. dep. of eff. corr.)");
        delete hDummy;
      }
    }
      
    
    // Sys. error of mult. dep., secondaries, protons only!
    if (sysErrorTypeMC != kErrorsWithoutMultDep /* && sysErrorTypeMC != kErrorsForMerging && sysErrorTypeMC != kErrorsForMergingOnlyMultDep*/) {
      Int_t species = AliPID::kProton;
      
      if (hSecSys[species]) {
        TH1D* hDummy = new TH1D(*hSecSys[species]);
        
        for (Int_t i = 1; i <= hSecSys[species]->GetNbinsX(); i++) {
          // Add the error in quadrature to the sys error
          const Double_t oldSysError = hSecSys[species]->GetBinError(i);
          const Double_t newSysError = hDummy->GetBinContent(i) * sysErrMultDepSecCorr_pr;
          
          hSecSys[species]->SetBinError(i, TMath::Sqrt(oldSysError*oldSysError + newSysError*newSysError));
          hDummy->SetBinError(i, newSysError);
        }
        
        hMCRelSysErrorSecMultDep[species] = getRelErrorHist(hDummy,
                                                            Form("hMCRelSysErrorSecMultDep_%s", AliPID::ParticleShortName(species)),
                                                            "Rel. sys. error (mult. dep. of secondary corr.)");
        delete hDummy;
      }
    }
    
    
    // -------------------------------------------------------------------------------------------------------
    // Same for the to-pion ratios
    
    // Sys error of strangeness scaling.
    // Same comments as for the yields apply.
    if (scaleStrangeness &&
        sysErrorTypeMC != kErrorsOnlyMultDep && sysErrorTypeMC != kErrorsForMergingOnlyMultDep) {
      for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
        if (!hSecToPiRatioSys[species])
          continue;
        
        TH1D* hDummy = new TH1D(*hSecToPiRatioSys[species]);
        hDummy->SetDirectory(0);
        
        for (Int_t i = 1; i <= hSecToPiRatioSys[species]->GetNbinsX(); i++) {
          // Add the error in quadrature to the sys error
          const Double_t oldSysError = hSecToPiRatioSys[species]->GetBinError(i);
          const Double_t newSysError = TMath::Abs(hSecSSToPiRatio[species]->GetBinContent(i) - hSecToPiRatio[species]->GetBinContent(i));
          
          hSecToPiRatioSys[species]->SetBinError(i, TMath::Sqrt(oldSysError*oldSysError + newSysError*newSysError));
          hDummy->SetBinError(i, newSysError);
        }
        
        hMCRelSysErrorToPiRatioSS[species] = getRelErrorHist(hDummy, Form("hMCRelSysErrorToPiRatioSS_%s", AliPID::ParticleShortName(species)),
                                                             "Rel. sys. error (MC corr, SS)");
        delete hDummy;
      }
    }
    
    // Sys error of GEANT-FLUKA correction.
    // Similarly as for strangeness scaling....
    if (correctGeantFluka &&
        sysErrorTypeMC != kErrorsOnlyMultDep && sysErrorTypeMC != kErrorsForMerging && sysErrorTypeMC != kErrorsForMergingOnlyMultDep) {
      // Although only K and p are affected, just loop over all species...
      for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
        if (!hEfficiencyToPiRatioSys[species])
          continue;
        
        TH1D* hDummy = new TH1D(*hEfficiencyToPiRatioSys[species]);
        hDummy->SetDirectory(0);
        
        for (Int_t i = 1; i <= hEfficiencyToPiRatioSys[species]->GetNbinsX(); i++) {
          // Add the error in quadrature to the sys error
          const Double_t oldSysError = hEfficiencyToPiRatioSys[species]->GetBinError(i);
          const Double_t newSysError = TMath::Abs(hEfficiencyToPiRatioWithGF[species]->GetBinContent(i) -
                                                  hEfficiencyToPiRatioNoGF[species]->GetBinContent(i));
          
          hEfficiencyToPiRatioSys[species]->SetBinError(i, TMath::Sqrt(oldSysError*oldSysError + newSysError*newSysError));
          hDummy->SetBinError(i, newSysError);
        }
        
        hMCRelSysErrorToPiRatioGF[species] = getRelErrorHist(hDummy, Form("hMCRelSysErrorToPiRatioGF%s", AliPID::ParticleShortName(species)),
                                                             "Rel. sys. error (MC corr, GF)");
        delete hDummy;
      }
    }
  }
  
  /* Sys. error drops out because efficiency mult dep. ~species-independent
  // Multiplicity dependencies, efficiency
  if (sysErrorTypeMC != kErrorsWithoutMultDep) {
    for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
      if (!hEfficiencyToPiRatioSys[species])
        continue;
      
      TH1D* hDummy = new TH1D(*hEfficiencyToPiRatioSys[species]);
      
      for (Int_t i = 1; i <= hEfficiencyToPiRatioSys[species]->GetNbinsX(); i++) {
          // Add the error in quadrature to the sys error
          const Double_t oldSysError = hEfficiencyToPiRatioSys[species]->GetBinError(i);
          const Double_t newSysError = hDummy->GetBinContent(i) * (isLowestMultBin ? sysErrMultDepEffCorrLowestMult : sysErrMultDepEffCorr);
          
          hEfficiencyToPiRatioSys[species]->SetBinError(i, TMath::Sqrt(oldSysError*oldSysError + newSysError*newSysError));
          hDummy->SetBinError(i, newSysError);
      }
      
      hMCRelSysErrorToPiRatioEffMultDep[species] = getRelErrorHist(hDummy,
                                                                    Form("hMCRelSysErrorToPiRatioEffMultDep_%s", 
                                                                        AliPID::ParticleShortName(species)),
                                                                    "Rel. sys. error (mult. dep. of eff. corr.)");
      delete hDummy;
    }
  }
  */
  
  // Sys. error of mult. dep., secondaries, protons only (identical to yield, since pions do not really change)!
  if (sysErrorTypeMC != kErrorsWithoutMultDep /* && sysErrorTypeMC != kErrorsForMerging && sysErrorTypeMC != kErrorsForMergingOnlyMultDep*/) {
    Int_t species = AliPID::kProton;
    
    if (hSecToPiRatioSys[species]) {
      TH1D* hDummy = new TH1D(*hSecToPiRatioSys[species]);
      
      for (Int_t i = 1; i <= hSecToPiRatioSys[species]->GetNbinsX(); i++) {
        // Add the error in quadrature to the sys error
        const Double_t oldSysError = hSecToPiRatioSys[species]->GetBinError(i);
        const Double_t newSysError = hDummy->GetBinContent(i) * sysErrMultDepSecCorr_pr;
        
        hSecToPiRatioSys[species]->SetBinError(i, TMath::Sqrt(oldSysError*oldSysError + newSysError*newSysError));
        hDummy->SetBinError(i, newSysError);
      }
      
      hMCRelSysErrorToPiRatioSecMultDep[species] = getRelErrorHist(hDummy,
                                                                   Form("hMCRelSysErrorToPiRatioSecMultDep_%s", 
                                                                        AliPID::ParticleShortName(species)),
                                                                   "Rel. sys. error (mult. dep. of secondary corr.)");
      delete hDummy;
    }
  }
  
  
  // Obtain the total correction factors and errors
  TH1D* hMCcorrTotal[AliPID::kSPECIES] = { 0x0, };
  TH1D* hMCcorrTotalSysError[AliPID::kSPECIES] = { 0x0, };
  TH1D* hMCcorrTotalRelSysError[AliPID::kSPECIES] = { 0x0, };
  
  for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
    // Correction factor histos can have different binning than data histos (constant factor above some threshold)
    // -> Need special functions to multiply and divide such histos
    
    // Bin-by-bin correction factors are multiplied, errors assumed to be independent
    
    // Eff
    hMCcorrTotal[species] = new TH1D(*hEfficiency[species]);
    hMCcorrTotal[species]->SetName(Form("hMCcorrTotal_%s", AliPID::ParticleShortName(species)));
    
    // Sec
    multiplyHistsDifferentBinning(hMCcorrTotal[species], scaleStrangeness ? hSecSS[species] : hSec[species], 1., 1.);
    
    // Mu corr
    if (applyMuonCorrection && species == AliPID::kPion)
      multiplyHistsDifferentBinning(hMCcorrTotal[species], scaleStrangeness ? hPiFracInMuPiStrangeScale : hPiFracInMuPi, 1., 1.);
    
    
    if (addMCsysErrors) {
      // Eff (+ GF + 10d10e + mult. dep., if used)
      hMCcorrTotalSysError[species] = new TH1D(*hEfficiencySys[species]);
      hMCcorrTotalSysError[species]->SetName(Form("hMCcorrTotalSysError_%s", AliPID::ParticleShortName(species)));
      
      // Sec (+ SS + mult. dep., if used)
      multiplyHistsDifferentBinning(hMCcorrTotalSysError[species], hSecSys[species], 1., 1.);
      
      // Mu corr (partially + SS, if used)
      if (applyMuonCorrection && species == AliPID::kPion)
        multiplyHistsDifferentBinning(hMCcorrTotalSysError[species], hMCsysErrorMuCorr, 1., 1.);
      
      hMCcorrTotalRelSysError[species] = getRelErrorHist(hMCcorrTotalSysError[species],
                                                         Form("hMCcorrTotalRelSysError_%s", AliPID::ParticleShortName(species)),
                                                         "Rel. sys. error (MC corr, tot)");
    }
    
    hMCcorrTotal[species]->GetYaxis()->SetRangeUser(0.5, 2.0);
    hMCcorrTotal[species]->GetYaxis()->SetTitle("Total correction factor");
    //hMCcorrTotal[species]->GetYaxis()->SetTitleSize(0.05);
    //hMCcorrTotal[species]->GetYaxis()->SetTitleOffset(0.8);
    if (iObs == kTrackPt)
      hMCcorrTotal[species]->GetXaxis()->SetRangeUser(0.15, actualUpperJetPt > 0 ? TMath::Max(30., actualUpperJetPt) : 50);
  }
  
  if (!restrictCharge) {
    TCanvas* cMCTotCorrYields = createCanvasForCorrFactor("cMCTotCorrYields",""); //ww=760, wh=420
    cMCTotCorrYields->SetLogx(iObs == kTrackPt);
    cMCTotCorrYields->SetGrid(0,0);
    hMCcorrTotal[AliPID::kPion]->Draw();
    hMCcorrTotal[AliPID::kKaon]->Draw("same");
    hMCcorrTotal[AliPID::kProton]->Draw("same");
    Bool_t isNotInclusive = actualUpperJetPt > 0;
    TLegend* legYield = createLegend(cMCTotCorrYields, "l", isNotInclusive ? 0.4 : 0.5, 0.7, isNotInclusive ? 0.8 : 0.9, 0.94);
    legYield->SetHeader(isNotInclusive ? Form("#it{p}_{T}^{jet, ch} = %.0f-%.0f GeV/#it{c}", actualLowerJetPt, actualUpperJetPt) 
                                       : "Inclusive");
    legYield->Draw();
    
    ClearTitleFromHistoInCanvas(cMCTotCorrYields);
    
    TString saveNameYield = Form("%s/MCCorrTot_%s_yield%s%s%s.pdf", pathSaveData.Data(), obsStringMCsysError[iObs].Data(),
                                     centralityString.Data(), jetString.Data(), chargeString.Data());
    cMCTotCorrYields->SaveAs(saveNameYield.Data());
  }
  
  TCanvas* cMCSysErrorsStackedYields[AliPID::kSPECIES] = { 0x0, };
  
  // Plot stacked errors. ALSO properly adjusts error for zero correction factor!
  for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
    if (species == AliPID::kMuon || !addMCsysErrors) continue;
    
    TString saveName = Form("%s/MCCorrRelSysError_%s_yield%s%s_%s%s.pdf", pathSaveData.Data(), obsStringMCsysError[iObs].Data(),
                            centralityString.Data(), jetString.Data(), AliPID::ParticleShortName(species), chargeString.Data());
    cMCSysErrorsStackedYields[species] = plotSysErrorsStacked(species, hMCRelSysErrorEff_eff, hMCRelSysErrorEff_res, hMCRelSysErrorEff_shape, 
                                                              hMCRelSysErrorSec, hMCRelSysErrorSS,
                                                              hMCRelSysErrorGF, hMCcorrTotalRelSysError,
                                                              hMCRelSysErrorEff10d10e,
                                                              species == AliPID::kPion ? hMCRelSysErrorMuCorr : 0x0,
                                                              hMCRelSysErrorSecMultDep, hMCRelSysErrorEffMultDep,
                                                              Form("cMCSysErrorsStackedYields_%s", AliPID::ParticleShortName(species)),
                                                              Form("%s in pp #sqrt{s}=7 TeV, %s",
                                                                   AliPID::ParticleLatexName(species),
                                                                   actualLowerJetPt < 0 && actualUpperJetPt < 0
                                                                     ? "inclusive"
                                                                     : Form("#it{p}_{T}^{jet, ch} = %.0f-%.0f GeV/#it{c}", actualLowerJetPt, 
                                                                            actualUpperJetPt)),
                                                              saveName.Data(), actualUpperJetPt, iObs, sysErrorTypeMC);
  }
  
  
  // Set corr. factors to negative value if corrFac(Eff,acc,res) larger than threshold
  // (or if factor is zero (i.e. there is no efficiency corr factor)) or (total) rel. sys error larger than threshold.
  // Do not do this in case of correction of MC-ID (otherwise fractions etc. will be wrong in these bins).
  
  if (!correctMCID) {
    printf("\n**********************\nChecking validity of bins:\n");
    
    for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
      if (!hMCcorrTotalRelSysError[species] || !hMCcorrTotalSysError[species]) {
        printf("Skipping species %s due to lacking histograms for MC!\n", AliPID::ParticleShortName(species));
        continue;
      }
      for (Int_t bin = 1; bin <= hMCcorrTotal[species]->GetNbinsX(); bin++) {
        const Bool_t corrFacTooLarge = hEfficiency[species]->GetBinContent(bin) > thresholdCF || hEfficiency[species]->GetBinContent(bin) <= 0.;
        const Bool_t corrFacErrTooLarge = hMCcorrTotalRelSysError[species]->GetBinContent(bin) > thresholdCFerror;
        if (corrFacTooLarge || corrFacErrTooLarge) {
          hMCcorrTotal[species]->SetBinContent(bin, -999);
          hMCcorrTotal[species]->SetBinError(bin, 0);
          
          hMCcorrTotalSysError[species]->SetBinContent(bin, -999);
          hMCcorrTotalSysError[species]->SetBinError(bin, 0);
          
          if (iObs != kTrackPt ||
              (iObs == kTrackPt && hMCcorrTotal[species]->GetXaxis()->GetBinLowEdge(bin) >= 0.15 &&
              ((actualUpperJetPt < 0 && actualLowerJetPt < 0) ||
                (actualUpperJetPt > 0 && hMCcorrTotal[species]->GetXaxis()->GetBinLowEdge(bin) < actualUpperJetPt))))
          printf("species %s, bin %f-%f not used: corrFacTooLarge %d, corrFacErrTooLarge %d\n", AliPID::ParticleShortName(species),
                hMCcorrTotal[species]->GetXaxis()->GetBinLowEdge(bin), hMCcorrTotal[species]->GetXaxis()->GetBinUpEdge(bin),
                corrFacTooLarge, corrFacErrTooLarge);
        }
      }
    }
    
    printf("\n**********************\n");
  }

  
  // Correct the yields with the efficiencies and primary fractions
  
  // Just for an estimate (without proper errors), sum up uncorrected yields (gives total raw yield by construction!) and correct
  // with correction factor for all particles!
  TH1D* hYieldTotalNoPID = new TH1D(*hYield[AliPID::kPion]);
  hYieldTotalNoPID->SetName("hYieldTotalNoPID");
  hYieldTotalNoPID->SetLineColor(kBlack);
  hYieldTotalNoPID->SetMarkerColor(kBlack);
  hYieldTotalNoPID->Add(hYield[AliPID::kElectron]);
  hYieldTotalNoPID->Add(hYield[AliPID::kMuon]);
  hYieldTotalNoPID->Add(hYield[AliPID::kKaon]);
  hYieldTotalNoPID->Add(hYield[AliPID::kProton]);
  
  TH1D* hYieldTotalNoPIDCorrected = new TH1D(*hYieldTotalNoPID);
  hYieldTotalNoPIDCorrected->SetName("hYieldTotalNoPIDCorrected");
  
  divideHistsDifferentBinning(hYieldTotalNoPIDCorrected, correctGeantFluka ? hEfficiencyAllNoPIDWithGF : hEfficiencyAllNoPIDNoGF, 1., 1.);
  multiplyHistsDifferentBinning(hYieldTotalNoPIDCorrected, scaleStrangeness ? hSecSSall : hSecAll, 1., 1.);
  
  
  TH1D* hYieldCorrected[AliPID::kSPECIES] = { 0x0, };
  TH1D* hYieldCorrectedSysError[AliPID::kSPECIES] = { 0x0, };
  for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
    hYieldCorrected[species] = new TH1D(*hYield[species]);
    hYieldCorrected[species]->SetName(Form("%s_corrected", hYield[species]->GetName()));
    hYieldCorrected[species]->SetTitle(Form("%s", hYield[species]->GetTitle()));
    //hYieldCorrected[species]->SetTitle(Form("%s, secondary and efficiency x acceptance x pT resolution corrected",
    //                                        hYield[species]->GetTitle()));
    
    // Correction factor histos can have different binning than data histos (constant factor above some threshold)
    // -> Need special functions to multiply and divide such histos
    
    multiplyHistsDifferentBinning(hYieldCorrected[species], hMCcorrTotal[species], 1., 1.);
    
    if (hYieldSysError[species]) {
      hYieldCorrectedSysError[species] = new TH1D(*hYieldSysError[species]);
      hYieldCorrectedSysError[species]->SetName(Form("%s_corrected", hYieldSysError[species]->GetName()));
      hYieldCorrectedSysError[species]->SetTitle(Form("%s", hYieldSysError[species]->GetTitle()));
      
      // Ignore the sys errors, if addMCsysErrors is kFALSE and just use the histos with the correction factors + stat errors
      multiplyHistsDifferentBinning(hYieldCorrectedSysError[species],
                                    addMCsysErrors ? hMCcorrTotalSysError[species] : hMCcorrTotal[species], 1., 1., !addMCsysErrors);
    }
  }
  
  
  // Calculate the total corrected yield. The original total yield had no error (just the number of detected tracks in a obs bin),
  // but due to the correction there is some error for the total yield. Also the error of the fractions introduces uncertainties
  // for the yields of individual species
  // Note that due to correlations (not taken into account) the error of the total yield is not correct
  TH1D* hYieldCorrectedTotal = new TH1D(*hYieldCorrected[0]);
  hYieldCorrectedTotal->SetLineColor(kBlack);
  hYieldCorrectedTotal->SetMarkerColor(kBlack);
  hYieldCorrectedTotal->SetName("hYieldCorrectedTotal");
  hYieldCorrectedTotal->SetTitle("Total");
  //hYieldCorrectedTotal->SetTitle("Total yield, secondary and efficiency x acceptance x pT resolution corrected");
  
  for (Int_t bin = 1; bin <= hYieldCorrectedTotal->GetNbinsX(); bin++) {
    // Some bins (e.g. those that are rejected) have value -999 -> This should not be summed up, but set to zero!
    if (hYieldCorrectedTotal->GetBinContent(bin) <= 0) {
      hYieldCorrectedTotal->SetBinContent(bin, 0.);
      hYieldCorrectedTotal->SetBinError(bin, 0.);
    }
    
    // Assume uncorrelated errors (actually, this is not true)
    for (Int_t i = 1; i < AliPID::kSPECIES; i++) {
      const Double_t contSpec = hYieldCorrected[i]->GetBinContent(bin);
      
      if (contSpec > 0) {
        const Double_t contTot = hYieldCorrectedTotal->GetBinContent(bin);
        const Double_t errTot = hYieldCorrectedTotal->GetBinError(bin);
        
        const Double_t errSpec = hYieldCorrected[i]->GetBinError(bin);
        
        hYieldCorrectedTotal->SetBinContent(bin, contTot + contSpec);
        hYieldCorrectedTotal->SetBinError(bin, TMath::Sqrt(errTot*errTot + errSpec*errSpec));
      }
    }
  }
  
  // Calculate the corrected fractions
  TH1D* hFractionCorrected[AliPID::kSPECIES] = { 0x0, };
  TH1D* hFractionCorrectedSysError[AliPID::kSPECIES] = { 0x0, };
  for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
    hFractionCorrected[species] = new TH1D(*hYield[species]);
    TString oldName = hYield[species]->GetName();
    TString newName = oldName.ReplaceAll("Yield", "Fraction");
    TString oldTitle = hYield[species]->GetTitle();
    TString newTitle = oldTitle.ReplaceAll("Yield", "Fraction");
    hFractionCorrected[species]->SetName(Form("%s_corrected", newName.Data()));
    hFractionCorrected[species]->SetTitle(Form("%s", AliPID::ParticleLatexName(species)));
    hFractionCorrected[species]->GetYaxis()->SetTitle("Corrected Fraction");
    hFractionCorrected[species]->GetXaxis()->SetMoreLogLabels(kTRUE);
    hFractionCorrected[species]->GetXaxis()->SetNoExponent(kTRUE);

    // Binomial error as for efficiencies (numerator and denominator are statistically not independent) for correct error calculation
    // (numerator is a "selected" subset of the denominator). It doesn't matter that the error of the histos is not just "sqrt(content)"
    // because the error formula also works for weighted histograms (which means that the error can be more or less anything)
    hFractionCorrected[species]->Divide(hYieldCorrected[species], hYieldCorrectedTotal, 1., 1., "B"); 
    

    //  The systematic errors just need to be scaled in the same way as the fractions.
    // So, just take the ratio to the uncorrected fraction and scale the sys. error accordingly
    // or, in this case, just divide by the same total yield as for yield -> fractions
    if (hYieldCorrectedSysError[species]) {
      hFractionCorrectedSysError[species] = new TH1D(*hFractionCorrected[species]);
      hFractionCorrectedSysError[species]->SetName(Form("%s_sysError", hFractionCorrected[species]->GetName()));
      hFractionCorrectedSysError[species]->SetTitle(Form("%s (sys. error)", hFractionCorrected[species]->GetTitle()));
      
      for (Int_t binX = 1; binX <= hFractionCorrectedSysError[species]->GetNbinsX(); binX++) {
        const Double_t corrTotalYield = hYieldCorrectedTotal->GetBinContent(binX);
        const Double_t scaleFactor = corrTotalYield > 0 ? 1.0 / corrTotalYield : 1.;
        hFractionCorrectedSysError[species]->SetBinError(binX,   hYieldCorrectedSysError[species]->GetBinError(binX) * scaleFactor);
      }
    }
  }
  
  // It may happen that a rejected bin has zero content before the correction. In such a case, the point is still visible with huge error bar.
  // Therefore, manually set points for these cases.
  for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
    if (!hMCcorrTotal[species])
      continue;
    for (Int_t bin = 1; bin <= hMCcorrTotal[species]->GetNbinsX(); bin++) {
      if (hMCcorrTotal[species]->GetBinContent(bin) < 0) {
        if (hYieldCorrected[species]) {
          hYieldCorrected[species]->SetBinContent(bin, -999);
          hYieldCorrected[species]->SetBinError(bin, 0);
        }
        if (hYieldCorrectedSysError[species]) {
          hYieldCorrectedSysError[species]->SetBinContent(bin, -999);
          hYieldCorrectedSysError[species]->SetBinError(bin, 0);
        }
        if (hFractionCorrected[species]) {
          hFractionCorrected[species]->SetBinContent(bin, -999);
          hFractionCorrected[species]->SetBinError(bin, 0);
        }
        if (hFractionCorrectedSysError[species]) {
          hFractionCorrectedSysError[species]->SetBinContent(bin, -999);
          hFractionCorrectedSysError[species]->SetBinError(bin, 0);
        }
      }
    }
  }
  
  
  
  // To-pion ratios: Total uncertainties and corrections
  
  // Obtain the total correction factors and errors
  TH1D* hMCcorrTotalToPiRatio[AliPID::kSPECIES] = { 0x0, };
  TH1D* hMCcorrTotalToPiRatioSysError[AliPID::kSPECIES] = { 0x0, };
  TH1D* hMCcorrTotalToPiRatioRelSysError[AliPID::kSPECIES] = { 0x0, };
  
  
  for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
    if (species == AliPID::kPion)
      continue;
    // Correction factor histos can have different binning than data histos (constant factor above some threshold)
    // -> Need special functions to multiply and divide such histos
    
    // Bin-by-bin correction factors are multiplied, errors assumed to be independent
    
    // Eff
    hMCcorrTotalToPiRatio[species] = new TH1D(*hEfficiencyToPiRatio[species]);
    hMCcorrTotalToPiRatio[species]->SetName(Form("hMCcorrTotalToPiRatio_%s", AliPID::ParticleShortName(species)));
    
    // Sec
    multiplyHistsDifferentBinning(hMCcorrTotalToPiRatio[species], scaleStrangeness ? hSecSSToPiRatio[species] : hSecToPiRatio[species],
                                  1., 1.);
    
    // Mu corr
    if (applyMuonCorrection)
      multiplyHistsDifferentBinning(hMCcorrTotalToPiRatio[species],
                                    scaleStrangeness ? hPiFracInMuPiToPiRatioStrangeScale : hPiFracInMuPiToPiRatio, 1., 1.);
    
    
    if (addMCsysErrors) {
      // Eff (+ GF, if used)
      hMCcorrTotalToPiRatioSysError[species] = new TH1D(*hEfficiencyToPiRatioSys[species]);
      hMCcorrTotalToPiRatioSysError[species]->SetName(Form("hMCcorrTotalToPiRatioSysError_%s", AliPID::ParticleShortName(species)));
      
      // Sec (+ SS, if used)
      multiplyHistsDifferentBinning(hMCcorrTotalToPiRatioSysError[species], hSecToPiRatioSys[species], 1., 1.);
      
      // Mu corr (partially + SS, if used)
      if (applyMuonCorrection)
        multiplyHistsDifferentBinning(hMCcorrTotalToPiRatioSysError[species], hMCsysErrorToPiRatioMuCorr, 1., 1.);
      
      hMCcorrTotalToPiRatioRelSysError[species] = getRelErrorHist(hMCcorrTotalToPiRatioSysError[species],
                                                         Form("hMCcorrTotalToPiRatioRelSysError_%s", AliPID::ParticleShortName(species)),
                                                         "Rel. sys. error (MC corr, tot)");
    }
    
    hMCcorrTotalToPiRatio[species]->GetYaxis()->SetRangeUser(0.5, 2.0);
    hMCcorrTotalToPiRatio[species]->GetYaxis()->SetTitle("Total correction factor");
    //hMCcorrTotalToPiRatio[species]->GetYaxis()->SetTitleSize(0.05);
    //hMCcorrTotalToPiRatio[species]->GetYaxis()->SetTitleOffset(0.8);
    if (iObs == kTrackPt)
      hMCcorrTotalToPiRatio[species]->GetXaxis()->SetRangeUser(0.15, actualUpperJetPt > 0 ? TMath::Max(30., actualUpperJetPt) : 50);
  }
  
  if (!restrictCharge) {
    TCanvas* cMCTotCorrToPiRatios = createCanvasForCorrFactor("cMCTotCorrToPiRatios",""); //ww=760, wh=420
    cMCTotCorrToPiRatios->SetLogx(iObs == kTrackPt);
    cMCTotCorrToPiRatios->SetGrid(0,0);
    hMCcorrTotalToPiRatio[AliPID::kKaon]->Draw("");
    hMCcorrTotalToPiRatio[AliPID::kProton]->Draw("same");
    Bool_t isNotInclusive = actualUpperJetPt > 0;
    TLegend* legToPiRatio = createLegend(cMCTotCorrToPiRatios, "l", isNotInclusive ? 0.4 : 0.5, 0.7, isNotInclusive ? 0.8 : 0.9, 0.94);
    legToPiRatio->SetHeader(isNotInclusive ? Form("#it{p}_{T}^{jet, ch} = %.0f-%.0f GeV/#it{c}", actualLowerJetPt, actualUpperJetPt) 
                                           : "Inclusive");
    
    ClearTitleFromHistoInCanvas(cMCTotCorrToPiRatios);
    
    TString saveNameToPiRatio = Form("%s/MCCorrTot_%s_toPiRatio%s%s%s.pdf", pathSaveData.Data(), obsStringMCsysError[iObs].Data(),
                                     centralityString.Data(), jetString.Data(), chargeString.Data());
    cMCTotCorrToPiRatios->SaveAs(saveNameToPiRatio.Data());
  }
  
  
  
  
  TCanvas* cMCSysErrorsStackedToPiRatios[AliPID::kSPECIES] = { 0x0, };
  
  // Plot stacked errors. ALSO properly adjusts error for zero correction factor!
  for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
    if (species == AliPID::kMuon || species == AliPID::kPion || !addMCsysErrors) continue;
    // Always forward mu corr histo for to-pion ratios!
    TString saveName = Form("%s/MCCorrRelSysError_%s_toPiRatio%s%s_%s%s.pdf", pathSaveData.Data(), obsStringMCsysError[iObs].Data(),
                            centralityString.Data(), jetString.Data(), AliPID::ParticleShortName(species), chargeString.Data());
    cMCSysErrorsStackedToPiRatios[species] = plotSysErrorsStacked(species, hMCRelSysErrorToPiRatioEff_eff, hMCRelSysErrorToPiRatioEff_res,
                                                                  hMCRelSysErrorToPiRatioEff_shape, hMCRelSysErrorToPiRatioSec,
                                                                  hMCRelSysErrorToPiRatioSS, hMCRelSysErrorToPiRatioGF,
                                                                  hMCcorrTotalToPiRatioRelSysError,
                                                                  0x0,
                                                                  hMCRelSysErrorToPiRatioMuCorr,
                                                                  hMCRelSysErrorToPiRatioSecMultDep, 0x0/*hMCRelSysErrorToPiRatioEffMultDep*/,
                                                                  Form("cMCSysErrorsStackedToPiRatios_%s", AliPID::ParticleShortName(species)),
                                                                  Form("%s/#pi in pp #sqrt{s}=7 TeV, %s",
                                                                   AliPID::ParticleLatexName(species),
                                                                   actualLowerJetPt < 0 && actualUpperJetPt < 0
                                                                     ? "inclusive"
                                                                     : Form("#it{p}_{T}^{jet, ch} = %.0f-%.0f GeV/#it{c}", actualLowerJetPt, 
                                                                            actualUpperJetPt)),
                                                                  saveName.Data(), actualUpperJetPt, iObs, sysErrorTypeMC);
  }
  
  
  // Take the same cut offs also for to-pion ratios (although errors are smaller, but we want to have the same ranges)
  for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
    if (!hMCcorrTotalToPiRatio[species])
      continue;
    for (Int_t bin = 1; bin <= hMCcorrTotalToPiRatio[species]->GetNbinsX(); bin++) {
      if (hMCcorrTotal[species]->GetBinContent(bin) < 0) {
        hMCcorrTotalToPiRatio[species]->SetBinContent(bin, -999);
        hMCcorrTotalToPiRatio[species]->SetBinError(bin, 0);
        
        hMCcorrTotalToPiRatioSysError[species]->SetBinContent(bin, -999);
        hMCcorrTotalToPiRatioSysError[species]->SetBinError(bin, 0);
      }
    }
  }
  
  // Correct the to-pi-ratios with the efficiencies and primary fractions
 
  TH1D* hRatioToPiCorrected[AliPID::kSPECIES] = { 0x0, };
  TH1D* hRatioToPiCorrectedSysError[AliPID::kSPECIES] = { 0x0, };
  for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
    if (!hRatioToPi[species])
      continue;
    
    hRatioToPiCorrected[species] = new TH1D(*hRatioToPi[species]);
    hRatioToPiCorrected[species]->SetName(Form("%s_corrected", hRatioToPi[species]->GetName()));
    hRatioToPiCorrected[species]->SetTitle(Form("%s", hRatioToPi[species]->GetTitle()));
    //hRatioToPiCorrected[species]->SetTitle(Form("%s, secondary and efficiency x acceptance x pT resolution corrected",
    //                                            hRatioToPi[species]->GetTitle()));
    
    // Correction factor histos can have different binning than data histos (constant factor above some threshold)
    // -> Need special functions to multiply and divide such histos
    
    
    
    multiplyHistsDifferentBinning(hRatioToPiCorrected[species], hMCcorrTotalToPiRatio[species], 1., 1.);
    
    if (hRatioToPiSysError[species]) {
      hRatioToPiCorrectedSysError[species] = new TH1D(*hRatioToPiSysError[species]);
      hRatioToPiCorrectedSysError[species]->SetName(Form("%s_corrected", hRatioToPiSysError[species]->GetName()));
      hRatioToPiCorrectedSysError[species]->SetTitle(Form("%s", hRatioToPiSysError[species]->GetTitle()));
      
      // Ignore the sys errors, if addMCsysErrors is kFALSE and just use the histos with the correction factors + stat errors
      multiplyHistsDifferentBinning(hRatioToPiCorrectedSysError[species],
                                    addMCsysErrors ? hMCcorrTotalToPiRatioSysError[species] : hMCcorrTotalToPiRatio[species], 1., 1.,
                                    !addMCsysErrors);
    }
  }
  
  
  // It may happen that a rejected bin has zero content before the correction. In such a case, the point is still visible with huge error bar.
  // Therefore, manually set points for these cases.
  for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
    if (!hMCcorrTotalToPiRatio[species])
      continue;
    for (Int_t bin = 1; bin <= hMCcorrTotalToPiRatio[species]->GetNbinsX(); bin++) {
      if (hMCcorrTotal[species]->GetBinContent(bin) < 0) {
        if (hRatioToPiCorrected[species]) {
          hRatioToPiCorrected[species]->SetBinContent(bin, -999);
          hRatioToPiCorrected[species]->SetBinError(bin, 0);
        }
        if (hRatioToPiCorrectedSysError[species]) {
          hRatioToPiCorrectedSysError[species]->SetBinContent(bin, -999);
          hRatioToPiCorrectedSysError[species]->SetBinError(bin, 0);
        }
      }
    }
  }
  
  
  
  if (addMCsysErrors) {
    saveFile->mkdir("MCsystematics");
    saveFile->cd("MCsystematics");
    
    if (hMCsysErrorMuCorr)
      hMCsysErrorMuCorr->Write();
    
    if (hMCsysErrorToPiRatioMuCorr)
      hMCsysErrorToPiRatioMuCorr->Write();
    
    for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
      if (hMCcorrTotal[species])
        hMCcorrTotal[species]->Write();
      
      
      if (hMCRelSysErrorEff_eff[species])
        hMCRelSysErrorEff_eff[species]->Write();
      
      if (hMCRelSysErrorEff_res[species])
        hMCRelSysErrorEff_res[species]->Write();
      
      if (hMCRelSysErrorEff_shape[species])
        hMCRelSysErrorEff_shape[species]->Write();
      
      if (hMCRelSysErrorEff10d10e[species])
        hMCRelSysErrorEff10d10e[species]->Write();
      
      if (hMCRelSysErrorSec[species])
        hMCRelSysErrorSec[species]->Write();
      
      if (hMCRelSysErrorSS[species])
        hMCRelSysErrorSS[species]->Write();
      
      if (hMCRelSysErrorGF[species])
        hMCRelSysErrorGF[species]->Write();
      
      if (hMCRelSysErrorMuCorr && species == AliPID::kPion)
        hMCRelSysErrorMuCorr->Write();
      
      if (hMCRelSysErrorSecMultDep[species])
        hMCRelSysErrorSecMultDep[species]->Write();
      
      if (hMCRelSysErrorEffMultDep[species])
        hMCRelSysErrorEffMultDep[species]->Write();
      
      if (hMCcorrTotalSysError[species])
        hMCcorrTotalSysError[species]->Write();
      
      if (hMCcorrTotalRelSysError[species])
        hMCcorrTotalRelSysError[species]->Write();
      
      
      if (hMCcorrTotalToPiRatio[species])
        hMCcorrTotalToPiRatio[species]->Write();
      
      if (hMCRelSysErrorToPiRatioEff_eff[species])
        hMCRelSysErrorToPiRatioEff_eff[species]->Write();
      
      if (hMCRelSysErrorToPiRatioEff_res[species])
        hMCRelSysErrorToPiRatioEff_res[species]->Write();
      
      if (hMCRelSysErrorToPiRatioEff_shape[species])
        hMCRelSysErrorToPiRatioEff_shape[species]->Write();
      
      if (hMCRelSysErrorToPiRatioSec[species])
        hMCRelSysErrorToPiRatioSec[species]->Write();
      
      if (hMCRelSysErrorToPiRatioSS[species])
        hMCRelSysErrorToPiRatioSS[species]->Write();
      
      if (hMCRelSysErrorToPiRatioGF[species])
        hMCRelSysErrorToPiRatioGF[species]->Write();
      
      if (hMCRelSysErrorToPiRatioMuCorr && species == AliPID::kPion)
        hMCRelSysErrorToPiRatioMuCorr->Write();
      
      //if (hMCRelSysErrorToPiRatioEffMultDep[species])
      //  hMCRelSysErrorToPiRatioEffMultDep[species]->Write();
      
      if (hMCRelSysErrorToPiRatioSecMultDep[species])
        hMCRelSysErrorToPiRatioSecMultDep[species]->Write();
        
      if (hMCcorrTotalToPiRatioSysError[species])
        hMCcorrTotalToPiRatioSysError[species]->Write();
      
      if (hMCcorrTotalToPiRatioRelSysError[species])
        hMCcorrTotalToPiRatioRelSysError[species]->Write();
      
      
       if (cMCSysErrorsStackedYields[species])
         cMCSysErrorsStackedYields[species]->Write();
       
       if (cMCSysErrorsStackedToPiRatios[species])
         cMCSysErrorsStackedToPiRatios[species]->Write();
    }
    
    saveFile->cd("");
  }
  
  
  // If MC is available, calculate the generated fractions
  TH1D* hMCgenPrimYieldTotal = 0x0;
  TH1D* hMCgenPrimFraction[AliPID::kSPECIES];
  for (Int_t i = 0; i < AliPID::kSPECIES; i++)
    hMCgenPrimFraction[i] = 0x0;
  
  if (numMCgenPrimYieldHistsFound > 0) {
    hMCgenPrimYieldTotal = new TH1D(*hMCgenPrimYield[0]);
    hMCgenPrimYieldTotal->SetLineColor(kBlack);
    hMCgenPrimYieldTotal->SetMarkerColor(kBlack);
    hMCgenPrimYieldTotal->SetName("hMCgenPrimYieldTotal");
    hMCgenPrimYieldTotal->SetTitle("Total (MC truth)");
    //hMCgenPrimYieldTotal->SetTitle("Total generated primary yield (MC truth)");
    
    for (Int_t i = 1; i < AliPID::kSPECIES; i++)
      hMCgenPrimYieldTotal->Add(hMCgenPrimYield[i], 1.);
    
    // Calculate the MC fractions
    for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
      hMCgenPrimYield[species]->SetTitle(Form("%s (MC truth)", AliPID::ParticleLatexName(species)));
      
      hMCgenPrimFraction[species] = new TH1D(*hMCgenPrimYield[species]);
      TString oldName = hMCgenPrimYield[species]->GetName();
      TString newName = oldName.ReplaceAll("Yield", "Fraction");
      TString oldTitle = hMCgenPrimYield[species]->GetTitle();
      TString newTitle = oldTitle.ReplaceAll("yield", "fraction");
      hMCgenPrimFraction[species]->SetName(newName.Data());
      hMCgenPrimFraction[species]->SetTitle(newTitle.Data());

      // Binomial error as for efficiencies (numerator and denominator are statistically not independent) for correct error calculation
      // (numerator is a "selected" subset of the denominator).
      hMCgenPrimFraction[species]->Divide(hMCgenPrimFraction[species], hMCgenPrimYieldTotal, 1., 1., "B"); 
    }
  }
  
  TCanvas* cCorrData = new TCanvas("cCorrData", "Corrected data", 0, 300, 900, 900);
  cCorrData->Divide(2, 1);//, 0., 0.01);
  if (iObs == kTrackPt)
    cCorrData->GetPad(1)->SetLogx(1);
  cCorrData->GetPad(1)->SetLogy(1);
  if (iObs == kTrackPt)
    cCorrData->GetPad(2)->SetLogx(1);
  cCorrData->GetPad(2)->SetLogy(1);
  
  cCorrData->GetPad(1)->SetRightMargin(0.001);
  cCorrData->GetPad(2)->SetRightMargin(0.001);
  
  cCorrData->GetPad(1)->SetLeftMargin(0.2);
  cCorrData->GetPad(2)->SetLeftMargin(0.2);
  
  cCorrData->cd(1); // uncorrected
  hYield[AliPID::kPion]->GetYaxis()->SetTitleOffset(1.4);
  hYield[AliPID::kPion]->Draw("E1");
  
  for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
    if (hYieldSysError[i])
      hYieldSysError[i]->Draw("E2 same");
    
    if (i == AliPID::kPion) continue;
    hYield[i]->Draw("E1 same");
  }
  
  ClearTitleFromHistoInCanvas(cCorrData, 1);
  
  cCorrData->cd(2); // corrected
  hYieldCorrectedTotal->GetYaxis()->SetTitleOffset(1.4);
  hYieldCorrectedTotal->GetXaxis()->SetMoreLogLabels(kTRUE);
  hYieldCorrectedTotal->GetXaxis()->SetNoExponent(kTRUE);
  hYieldCorrectedTotal->GetYaxis()->SetRangeUser(hYieldCorrected[AliPID::kMuon]->GetBinContent(hYieldCorrected[AliPID::kMuon]->FindLastBinAbove(0.)) * 0.1,
                                                 hYieldCorrectedTotal->GetBinContent(hYieldCorrectedTotal->GetMaximumBin()) * 10.);
  
  if (hMCgenPrimYieldTotal) {
    hMCgenPrimYieldTotal->GetYaxis()->SetTitleOffset(1.4);
    hMCgenPrimYieldTotal->GetXaxis()->SetMoreLogLabels(kTRUE);
    hMCgenPrimYieldTotal->GetXaxis()->SetNoExponent(kTRUE);
    hMCgenPrimYieldTotal->GetXaxis()->SetTitle(hYieldCorrectedTotal->GetXaxis()->GetTitle());
    hMCgenPrimYieldTotal->GetYaxis()->SetRangeUser(hYieldCorrected[AliPID::kMuon]->GetBinContent(hYieldCorrected[AliPID::kMuon]->FindLastBinAbove(0.)) * 0.1,
                                                   hYieldCorrectedTotal->GetBinContent(hYieldCorrectedTotal->GetMaximumBin()) * 10.);
  }
  
  if (hMCgenPrimYieldTotal) {
    hMCgenPrimYieldTotal->Draw("E1");
    hYieldCorrectedTotal->Draw("E1 same");
  }
  else
    hYieldCorrectedTotal->Draw("E1");
  
  for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
    if (hMCgenPrimYield[i])
      hMCgenPrimYield[i]->Draw("E1 same");
    hYieldCorrected[i]->Draw("E1 same");
  }
  
  TLegend* legTemp = cCorrData->cd(2)->BuildLegend(0.25, 0.16, 0.65, 0.51);
  legTemp->SetNColumns(2);
  legTemp->SetFillColor(kWhite);
  
  // Do not include in legend
  for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
    if (hYieldCorrectedSysError[i])
      hYieldCorrectedSysError[i]->Draw("E2 same");
  }
  
  ClearTitleFromHistoInCanvas(cCorrData, 2);
  
  TCanvas* cCorrYieldsRatio = 0x0;
  
  TH1D* hYieldCorrectedRatioToMC[AliPID::kSPECIES];
  for (Int_t i = 0; i < AliPID::kSPECIES; i++) 
    hYieldCorrectedRatioToMC[i] = 0x0;

  TH1D* hYieldCorrectedTotalRatioToMC = 0x0;
  
  if (numMCgenPrimYieldHistsFound > 0) {
    // Compare with MC truth
    for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
      hYieldCorrectedRatioToMC[species] = new TH1D(*hYieldCorrected[species]);
      hYieldCorrectedRatioToMC[species]->SetName(Form("hYieldCorrectedRatioToMC_%s", AliPID::ParticleShortName(species)));
      hYieldCorrectedRatioToMC[species]->SetTitle(Form("%s", AliPID::ParticleLatexName(species)));
      hYieldCorrectedRatioToMC[species]->GetYaxis()->SetTitle("Corrected Yield: Fit / MC Truth");
      hYieldCorrectedRatioToMC[species]->Divide(hMCgenPrimYield[species]);
      hYieldCorrectedRatioToMC[species]->SetMarkerStyle(20);
    }
    
    hYieldCorrectedTotalRatioToMC = new TH1D(*hYieldCorrectedTotal);
    hYieldCorrectedTotalRatioToMC->SetName("hYieldCorrectedTotalRatioToMC");
    hYieldCorrectedTotalRatioToMC->SetTitle("Total");
    hYieldCorrectedTotalRatioToMC->GetYaxis()->SetTitle("Corrected Yield: Fit / MC Truth");
    hYieldCorrectedTotalRatioToMC->Divide(hMCgenPrimYieldTotal);
    hYieldCorrectedTotalRatioToMC->SetMarkerStyle(20);
    
    cCorrYieldsRatio = new TCanvas("cCorrYieldsRatio", "Corrected Yields Comparison to MC", 0, 300, 900, 900);
    cCorrYieldsRatio->SetGridx(0);
    cCorrYieldsRatio->SetGridy(1);
    if (iObs == kTrackPt)
    cCorrYieldsRatio->SetLogx(1);
    
    hYieldCorrectedTotalRatioToMC->GetYaxis()->SetRangeUser(0.8, 1.2);
    hYieldCorrectedTotalRatioToMC->GetYaxis()->SetTitleOffset(0.85);
    hYieldCorrectedTotalRatioToMC->Draw("E1");
    
    for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
      if (applyMuonCorrection && species == AliPID::kMuon)
        continue;
      hYieldCorrectedRatioToMC[species]->Draw("E1 same");
    }
    
    cCorrYieldsRatio->BuildLegend()->SetFillColor(kWhite);
    
    ClearTitleFromHistoInCanvas(cCorrYieldsRatio);
  }
  
  
  if (cFractions)
    cFractions->Draw();
  
  TCanvas* cCorrFractions = new TCanvas("cCorrFractions", "Corrected particleFractions", 0, 300, 900, 900);
  if (iObs == kTrackPt)
    cCorrFractions->SetLogx(1);
  hFractionCorrected[0]->GetYaxis()->SetRangeUser(0., 1.);
  hFractionCorrected[0]->Draw("E1");
  if (hMCgenPrimFraction[0])
    hMCgenPrimFraction[0]->Draw("E1 same");
  
  for (Int_t i = 1; i < AliPID::kSPECIES; i++) {
    hFractionCorrected[i]->Draw("E1 same");
    if (hMCgenPrimFraction[i])
      hMCgenPrimFraction[i]->Draw("E1 same");
  }
  
  cCorrFractions->BuildLegend()->SetFillColor(kWhite);
  
  // Do not include in legend!!
  for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
    if (hFractionCorrectedSysError[i])
      hFractionCorrectedSysError[i]->Draw("E2 same");
  }
  
  
  ClearTitleFromHistoInCanvas(cCorrFractions);
  
  
  
  
  
  
  TCanvas* cCorrDataToPiRatio = 0x0;
  if (hRatioToPi[AliPID::kKaon]) {
    cCorrDataToPiRatio = new TCanvas("cCorrDataToPiRatio", "Corrected data", 0, 300, 900, 900);
    cCorrDataToPiRatio->Divide(2, 1);//, 0., 0.01);
    if (iObs == kTrackPt)
      cCorrDataToPiRatio->GetPad(1)->SetLogx(1);
    if (iObs == kTrackPt)
      cCorrDataToPiRatio->GetPad(2)->SetLogx(1);
    
    cCorrDataToPiRatio->GetPad(1)->SetRightMargin(0.001);
    cCorrDataToPiRatio->GetPad(2)->SetRightMargin(0.001);
    
    cCorrDataToPiRatio->GetPad(1)->SetLeftMargin(0.2);
    cCorrDataToPiRatio->GetPad(2)->SetLeftMargin(0.2);
    
    cCorrDataToPiRatio->cd(1); // uncorrected
    hRatioToPi[AliPID::kKaon]->GetYaxis()->SetTitleOffset(1.4);
    hRatioToPi[AliPID::kKaon]->GetYaxis()->SetRangeUser(0, 0.7);
    hRatioToPi[AliPID::kKaon]->Draw("E1");
    
    for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
      if (hRatioToPiSysError[i]) {
        hRatioToPiSysError[i]->GetYaxis()->SetRangeUser(0, 0.7);
        hRatioToPiSysError[i]->Draw("E2 same");
      }
      
      if (i == AliPID::kKaon || !hRatioToPi[i]) continue;
      
      hRatioToPi[i]->GetYaxis()->SetRangeUser(0, 0.7);
      hRatioToPi[i]->Draw("E1 same");
    }
    
    ClearTitleFromHistoInCanvas(cCorrDataToPiRatio, 1);
    
    cCorrDataToPiRatio->cd(2); // corrected
    
    for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
      if (hRatioToPiCorrected[i]) {
        hRatioToPiCorrected[i]->GetYaxis()->SetRangeUser(0, 0.7);
        hRatioToPiCorrected[i]->GetYaxis()->SetTitleOffset(1.4);
        hRatioToPiCorrected[i]->Draw(i == 0 ? "E1" : "E1 same");
      }
    }
    
    //legTemp = cCorrDataToPiRatio->cd(2)->BuildLegend(0.25, 0.16, 0.65, 0.51);
    legTemp = cCorrDataToPiRatio->cd(2)->BuildLegend(0.25, 0.75, 0.75, 0.88);
    legTemp->SetNColumns(2);
    legTemp->SetFillColor(kWhite);
    
    // Do not include in legend
    for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
      if (hRatioToPiCorrectedSysError[i]) {
        hRatioToPiCorrectedSysError[i]->GetYaxis()->SetRangeUser(0, 0.7);
        hRatioToPiCorrectedSysError[i]->Draw("E2 same");
      }
    }
    
    ClearTitleFromHistoInCanvas(cCorrDataToPiRatio, 2);
  }
  
  
  
  TH1D* hYieldCorrectedRapidity[AliPID::kSPECIES] = { 0x0, };
  TH1D* hYieldCorrectedSysErrorRapidity[AliPID::kSPECIES] = { 0x0, };
  TH1D* hYieldCorrectedTotalRapidity = 0x0;
  TH1D* hPileUpFraction = 0x0;
  TH1D* hYieldCorrectedSysErrorPileUp[AliPID::kSPECIES] = { 0x0, };
  TH1D* hYieldCorrectedSysErrorPileUpRapidity[AliPID::kSPECIES] = { 0x0, };
  
  TH1D* hRatioToPiCorrectedRapidity[AliPID::kSPECIES] = { 0x0, };
  TH1D* hRatioToPiCorrectedSysErrorRapidity[AliPID::kSPECIES] = { 0x0, };
  
  // Only do the following for inclusive
  if (!restrictJetPtAxis) {
    // Normalise to 1/N_inel
    // The yields are already normalised to 1/N_ev_trigger&vtx&zvtx. Just need to apply correction factor N_ev_trigger&vtx / N_ev_trigger
    // and eps_trigger
    
    // Try to get it from the histos directly. If the information is missing (old version), use global value defined at the top
    // of this file
    Double_t nEventCorrFactor = 0.;
    
    TH1* hNumEventsTriggerSel = dynamic_cast<TH1*>(fileData->Get("fhEventsTriggerSel"));
    TH1* hNumEventsTriggerSelVtxCut = dynamic_cast<TH1*>(fileData->Get("fhEventsTriggerSelVtxCut"));
    TH1* hNumEventsTriggerSelVtxCutZ = dynamic_cast<TH1*>(fileData->Get("fhEventsProcessedNoPileUpRejection"));
    TH1* hNumEventsTriggerSelVtxCutZPileUpRejected = dynamic_cast<TH1*>(fileData->Get("fhEventsProcessed"));
    
    Bool_t extractPileUpSysError = (hNumEventsTriggerSelVtxCutZ && hNumEventsTriggerSelVtxCutZPileUpRejected);
    Double_t pileUpFraction = 0.;
    
    if (hNumEventsTriggerSel && hNumEventsTriggerSelVtxCut) {
      printf("Loading correction factors for number of events from histograms....\n");
      
      // If desired, restrict centrality axis
      Int_t lowerCentralityBinLimitData = -1;
      Int_t upperCentralityBinLimitData = -2; // Integral(lowerCentBinLimit, uppCentBinLimit) will not be restricted if these values are kept. In particular, under- and overflow bin will be used!
      Bool_t restrictCentralityAxisData = kFALSE;
      Double_t actualLowerCentralityData = -1.;
      Double_t actualUpperCentralityData = -1.;
      
      if (lowerCentralityData >= -1 && upperCentralityData >= -1) {
        // Add subtract a very small number to avoid problems with values right on the border between to bins
        lowerCentralityBinLimitData = hNumEventsTriggerSel->GetXaxis()->FindFixBin(lowerCentralityData + 0.001);
        upperCentralityBinLimitData = hNumEventsTriggerSel->GetXaxis()->FindFixBin(upperCentralityData - 0.001);
        
        // Check if the values look reasonable
        if (lowerCentralityBinLimitData <= upperCentralityBinLimitData && lowerCentralityBinLimitData >= 0
            && upperCentralityBinLimitData <= hNumEventsTriggerSel->GetXaxis()->GetNbins() + 1) {
          actualLowerCentralityData = hNumEventsTriggerSel->GetXaxis()->GetBinLowEdge(lowerCentralityBinLimitData);
          actualUpperCentralityData = hNumEventsTriggerSel->GetXaxis()->GetBinUpEdge(upperCentralityBinLimitData);

          restrictCentralityAxisData = kTRUE;
        }
        else {
          std::cout << std::endl;
          std::cout << "Requested centrality range (data) out of limits or upper and lower limit are switched!" << std::endl;
          return kFALSE;
        }
      }
      
      if (!restrictCentralityAxisData) {
        GetAxisRangeForMultiplicityAxisForMB(hNumEventsTriggerSel->GetXaxis(), lowerCentralityBinLimitData, upperCentralityBinLimitData);
        actualLowerCentralityData = hNumEventsTriggerSel->GetXaxis()->GetBinLowEdge(TMath::Max(0, lowerCentralityBinLimitData));
        actualUpperCentralityData = hNumEventsTriggerSel->GetXaxis()->GetBinUpEdge(upperCentralityBinLimitData);
      }
      else if (doMultCorrectionSec) {
        // Consistency check
        Double_t actualCentralMultValue = 0.5 * (actualUpperCentralityData + actualLowerCentralityData);
        if (TMath::Abs(actualCentralMultValue - averageMultOfSample) > 1e-4) {
          std::cout << "Mismatch between actual central multiplicity value (" << actualCentralMultValue << ") and average multiplicity ("
                    << averageMultOfSample << ")!" << std::endl;
          exit(-1);
        }
      }
  
      std::cout << std::endl;
      if (restrictCentralityAxisData)
        std::cout << "Requested centrality range (data): " << actualLowerCentralityData << " - " << actualUpperCentralityData << std::endl; 
      else {
        std::cout << "Requested centrality range (data): MB (" << actualLowerCentralityData << " - " << actualUpperCentralityData<< ")"
                  << std::endl;
      }
      
      
      Double_t numEventsTriggerSel = hNumEventsTriggerSel->Integral(lowerCentralityBinLimitData, upperCentralityBinLimitData);
      Double_t numEventsTriggerSelVtxCut = hNumEventsTriggerSelVtxCut->Integral(lowerCentralityBinLimitData, upperCentralityBinLimitData);
      
      printf("Found events: trigger sel %f, trigger sel + vtx cut %f\n", numEventsTriggerSel, numEventsTriggerSelVtxCut);
      
      if (numEventsTriggerSel <= 0 || numEventsTriggerSelVtxCut <= 0) {
        printf("Error: Something is wrong with these event numbers!\n");
        return kFALSE;
      }
      else
        nEventCorrFactor = numEventsTriggerSelVtxCut / numEventsTriggerSel;
      
      if (extractPileUpSysError) {
        hPileUpFraction = new TH1D(*((TH1D*)hNumEventsTriggerSelVtxCutZPileUpRejected));
        hPileUpFraction->SetName("hPileUpFraction");
        hPileUpFraction->Reset();
        // Fraction of events surviving the pile-up rejection
        hPileUpFraction->Divide(hNumEventsTriggerSelVtxCutZPileUpRejected, hNumEventsTriggerSelVtxCutZ, 1., 1., "B");
        
        // Fraction of events rejected due to pile-up detection -> Also take into account under- and overflow!
        // -> Take TMath::Abs to be on the save side, also it should never happen that survival propability is larger unity
        for (Int_t i = 0; i <= hPileUpFraction->GetNbinsX() + 1; i++)
          hPileUpFraction->SetBinContent(i, TMath::Abs(1. - hPileUpFraction->GetBinContent(i)));
        
        // Cannot take integral of pile-up fraction, but must take calculate integrals for each class first and then form the ratio
        const Double_t numEvtsTot = hNumEventsTriggerSelVtxCutZ->Integral(lowerCentralityBinLimitData, upperCentralityBinLimitData);
        const Double_t numEvtsAfterRej = hNumEventsTriggerSelVtxCutZPileUpRejected->Integral(lowerCentralityBinLimitData,
                                                                                             upperCentralityBinLimitData);
        pileUpFraction = (numEvtsTot > 0) ? TMath::Abs(1. - numEvtsAfterRej / numEvtsTot) : 0.;
      }
    }
    else {
      printf("No histos found for correction factors for number of events from histograms....\n");
      nEventCorrFactor = corrNevVertexcut;
    }
    
    if (normaliseToNInel) {
      printf("Applying Nevt correction factor: %f\nAnd trigger efficiency: %f\n\n", nEventCorrFactor, eps_trigger);
      
      
      hYieldCorrectedTotal->Scale(nEventCorrFactor * eps_trigger);
      if (hMCgenPrimYieldTotal)
        hMCgenPrimYieldTotal->Scale(nEventCorrFactor * eps_trigger);
          
      for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
        hYieldCorrected[species]->Scale(nEventCorrFactor * eps_trigger);
        if (hYieldCorrectedSysError[species])
          hYieldCorrectedSysError[species]->Scale(nEventCorrFactor * eps_trigger);
        
        if (hMCgenPrimYield[species])
          hMCgenPrimYield[species]->Scale(nEventCorrFactor * eps_trigger);
      }
    }
    else 
      printf("NOT normalising to N_inel, but just to N_evt.....\n");
    
    if (extractPileUpSysError) {
      for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
        if (!hYieldCorrected[species])
          continue;
        
        hYieldCorrectedSysErrorPileUp[species] = new TH1D(*hYieldCorrected[species]);
        hYieldCorrectedSysErrorPileUp[species]->SetName(Form("%s_sysErrorPileUp", hYieldCorrected[species]->GetName()));
        hYieldCorrectedSysErrorPileUp[species]->Reset();
        
        for (Int_t i = 1; i <= hYieldCorrected[species]->GetNbinsX(); i++) {
          const Double_t yieldSpecies = hYieldCorrected[species]->GetBinContent(i);
          hYieldCorrectedSysErrorPileUp[species]->SetBinContent(i, yieldSpecies);
          // Take 100% of pile-up fraction as relative sys error from that source:
          // Delta yield = yield * Delta N / N =: yield * pileUpFraction
          hYieldCorrectedSysErrorPileUp[species]->SetBinError(i, yieldSpecies * pileUpFraction); 
        }
      }
      
      printf("Sys. error from pile-up has been estimated (pile-up fraction is %f)!\n", pileUpFraction);
    }
    
    // Set proper y-axis title
    TString yAxisTitle = Form("1/%s 1/(2#pip_{T}) d^{2}N/d#etadp_{T} (GeV/c)^{-2}", normaliseToNInel ? "N_{INEL}" : "N_{evt}");
    hYieldCorrectedTotal->GetYaxis()->SetTitle(yAxisTitle.Data());
    if (hMCgenPrimYieldTotal)
      hMCgenPrimYieldTotal->GetYaxis()->SetTitle(yAxisTitle.Data());
    
    for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
      hYieldCorrected[species]->GetYaxis()->SetTitle(yAxisTitle.Data());
      
      if (hYieldCorrectedSysError[species])
        hYieldCorrectedSysError[species]->GetYaxis()->SetTitle(yAxisTitle.Data());
      
      if (hMCgenPrimYield[species])
        hMCgenPrimYield[species]->GetYaxis()->SetTitle(yAxisTitle.Data());
    }
    
    // Convert every identified yield from dN/deta to dN/dy. For the total yield, take the new sum.
    // NOTE: The error of the total yield will be wrong because correlations are not taken into account!
    hYieldCorrectedTotalRapidity = new TH1D(*hYieldCorrectedTotal);
    hYieldCorrectedTotalRapidity->SetName("hYieldCorrectedTotalRapidity");
    hYieldCorrectedTotalRapidity->Reset();
    
    for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
      hYieldCorrectedRapidity[i] = convertYieldHistFromEtaToY(hYieldCorrected[i], i, etaAbsCut);
      hYieldCorrectedSysErrorRapidity[i] = convertYieldHistFromEtaToY(hYieldCorrectedSysError[i], i, etaAbsCut);
      hYieldCorrectedTotalRapidity->Add(hYieldCorrectedRapidity[i]);
      
      hYieldCorrectedSysErrorPileUpRapidity[i] = convertYieldHistFromEtaToY(hYieldCorrectedSysErrorPileUp[i], i, etaAbsCut);
    }
    
    // Similarly for the to-pion ratios
    for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
      hRatioToPiCorrectedRapidity[i] = convertToPionRatioHistFromEtaToY(hRatioToPiCorrected[i], i, etaAbsCut);
      hRatioToPiCorrectedSysErrorRapidity[i] = convertToPionRatioHistFromEtaToY(hRatioToPiCorrectedSysError[i], i, etaAbsCut);
    }
  }
  
  // Save results to file
  saveFile->cd();
  
  if (cSec)
    cSec->Write();
  
  if (cSecSS)
    cSecSS->Write();
  
  if (cSec2)
    cSec2->Write();
  
  if (cSec2MultCorr)
    cSec2MultCorr->Write();
  
  if (cSecSS2)
    cSecSS2->Write();
  
  if (cSecSS2MultCorr)
    cSecSS2MultCorr->Write();
  
  if (cSecToPiRatio)
    cSecToPiRatio->Write();
  
  if (cSecToPiRatioMultCorr)
    cSecToPiRatioMultCorr->Write();
  
  if (cSecSSToPiRatio)
    cSecSSToPiRatio->Write();
  
  if (cSecSSToPiRatioMultCorr)
    cSecSSToPiRatioMultCorr->Write();
  
  if (cCorrData)
    cCorrData->Write();
  
  if (cCorrYieldsRatio)
    cCorrYieldsRatio->Write();
  
  if (cFractions)
    cFractions->Write();
  
  if (cCorrFractions)
    cCorrFractions->Write();
  
  if (cCorrDataToPiRatio)
    cCorrDataToPiRatio->Write();
  
  if (hSecAll)
    hSecAll->Write();
  
  if (hSecSSall)
    hSecSSall->Write();
  
  for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
    if (hSec[i])
      hSec[i]->Write();
    
    if (hSecSS[i])
      hSecSS[i]->Write();
  }
  
  for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
    if (hSecToPiRatio[i])
      hSecToPiRatio[i]->Write();
    
    if (hSecSSToPiRatio[i])
      hSecSSToPiRatio[i]->Write();
  }
  
  if (hPiFracInMuPi)
    hPiFracInMuPi->Write();
  
  if (hPiFracInMuPiStrangeScale)
    hPiFracInMuPiStrangeScale->Write();
  
  if (hPiFracInMuPiToPiRatio)
    hPiFracInMuPiToPiRatio->Write();
  
  if (hPiFracInMuPiToPiRatioStrangeScale)
    hPiFracInMuPiToPiRatioStrangeScale->Write();
  
  if (hPileUpFraction)
    hPileUpFraction->Write();
  
  for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
    if (hYield[i])
      hYield[i]->Write();
    
    if (hYieldSysError[i])
      hYieldSysError[i]->Write();
  }
  
  for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
    if (hRatioToPi[i])
      hRatioToPi[i]->Write();
    
    if (hRatioToPiSysError[i])
      hRatioToPiSysError[i]->Write();
  }
  
  for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
    if (hRatioToPiCorrected[i])
      hRatioToPiCorrected[i]->Write();
    
    if (hRatioToPiCorrectedSysError[i])
      hRatioToPiCorrectedSysError[i]->Write();
  }
  
  if (hYieldTotalNoPID)
    hYieldTotalNoPID->Write();
  
  if (hYieldTotalNoPIDCorrected)
    hYieldTotalNoPIDCorrected->Write();
  
  if (hYieldCorrectedTotal)
    hYieldCorrectedTotal->Write();
  
  for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
    if (hYieldCorrected[i])
      hYieldCorrected[i]->Write();
    
    if (hYieldCorrectedSysError[i])
      hYieldCorrectedSysError[i]->Write();
    
    if (hYieldCorrectedSysErrorPileUp[i])
      hYieldCorrectedSysErrorPileUp[i]->Write();
  }
  
  for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
    if (hMCgenPrimYield[i])
      hMCgenPrimYield[i]->Write();
  }
  
  if (hMCgenPrimYieldTotal)
      hMCgenPrimYieldTotal->Write();
  
  for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
    if (hYieldCorrectedRatioToMC[i])
      hYieldCorrectedRatioToMC[i]->Write();
  }
  
  if (hYieldCorrectedTotalRatioToMC)
      hYieldCorrectedTotalRatioToMC->Write();
  
  for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
    if (hFractionCorrected[i])
      hFractionCorrected[i]->Write();
    
    if (hFractionCorrectedSysError[i])
      hFractionCorrectedSysError[i]->Write();
  }
  
  for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
    if (hMCgenPrimFraction[i])
      hMCgenPrimFraction[i]->Write();
  }
  
  if (hYieldCorrectedTotalRapidity)
    hYieldCorrectedTotalRapidity->Write();
  
  for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
    if (hYieldCorrectedRapidity[i])
      hYieldCorrectedRapidity[i]->Write();
    
    if (hYieldCorrectedSysErrorRapidity[i])
      hYieldCorrectedSysErrorRapidity[i]->Write();
    
    if (hYieldCorrectedSysErrorPileUpRapidity[i])
      hYieldCorrectedSysErrorPileUpRapidity[i]->Write();
    
    if (hRatioToPiCorrectedRapidity[i])
      hRatioToPiCorrectedRapidity[i]->Write();
    
    if (hRatioToPiCorrectedSysErrorRapidity[i])
      hRatioToPiCorrectedSysErrorRapidity[i]->Write();
    
  }
  
  TNamed* settings = new TNamed(
      Form("Settings: Efficiency file \"%s\", Data file \"%s\", MC sys error path \"%s\", correctGeantFluka %d, scaleStrangeness %d, applyMuonCorrection %d, lowerCentralityData %.3f, upperCentralityData %.3f, lowerCentrality %.3f, upperCentrality %.3f, lowerJetPt %.1f, upperJetPt %.1f, constantCorrectionAbovePtThreshold %.3f, rebinEfficiencyObs %d, iObs %d, eps_trigger %f, sysErrorTypeMC %d, normaliseToNInel %d, correctMCID %d, individualMultCorr %d, Mult dependence of sec corr file \"%s\", correctEff10d10e %d, isLowestMultBin %d\n",
           pathNameEfficiency.Data(), pathNameData.Data(), pathMCsysErrors.Data(), correctGeantFluka, scaleStrangeness, 
           applyMuonCorrection, lowerCentralityData, upperCentralityData, lowerCentrality, upperCentrality, lowerJetPt, upperJetPt,
           constantCorrectionAbovePtThreshold, rebinEfficiencyObs, iObs, eps_trigger, sysErrorTypeMC, normaliseToNInel, correctMCID,
           individualMultCorr, multCorrPathName.Data(), correctEff10d10e, isLowestMultBin), "");
  settings->Write();
  
  saveFile->Close();
  
  printf("\n");
  PrintSettingsAxisRangeForMultiplicityAxisForMB();
  printf("\n\n####IMPORTANT: Do NOT forget to adjust the sys errors according to filterbit, jet inclusive <-> normal inclusive!\n");
  return 0;
}