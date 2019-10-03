#ifndef PrepareRawSpectra_H
#define PrepareRawSpectra_H

#include "AliUtilTOFParams.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1F.h"
#include "TList.h"
#include "TOFMacroUtils.h"
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"
#include "UtilFiles.h"
#include "UtilHisto.h"
#include "UtilMessages.h"
#include "Utilities/UtilPlots.h"
#include "iostream"
#include "libs/TOFsignal.C"

#include "AliTOFTemplateFitter.h"
using std::cout;
using std::endl;

////////////////////////////////////////////////////////////////////////////
///                                                                       //
///                                                                       //
/// Include macro to compute the raw yields of pions, kaons and protons   //
/// Using the TOF signal                                                  //
///                                                                       //
/// Authors:                                                              //
/// N. Jacazio,  nicolo.jacazio[AROBASe]bo.infn.it                        //
////////////////////////////////////////////////////////////////////////////

void SetColor(TH1* h, Int_t c, Int_t f);

void PrepareRawSpectra(TString dirname = latestDir[0], TString dataname = "TOFDataHisto.root", TString fitprefix = "", UInt_t iCharge = 0, UInt_t iSpecies = 0, UInt_t iMult = 11, UInt_t fitmode = 0, TString prefix = "", const Bool_t fitsigma = kTRUE, const Bool_t save = kFALSE)
{
  Infomsg("PrepareRawSpectra", Form("TString dirname %s, TString dataname %s, TString fitprefix %s, const UInt_t iCharge %i, const UInt_t iSpecies %i, const UInt_t iMult %i, UInt_t fitmode %i, TString prefix %s, const Bool_t fitsigma %i", dirname.Data(), dataname.Data(), fitprefix.Data(), iCharge, iSpecies, iMult, fitmode, prefix.Data(), fitsigma));

  if (iCharge >= kCharges)
    Fatalmsg("PrepareRawSpectra", Form("Wrong charge: %i", iCharge));
  //
  if (iSpecies >= kSpecies)
    Fatalmsg("PrepareRawSpectra", Form("Wrong particle: %i", iSpecies));
  //
  if (fitmode >= kFModes)
    Fatalmsg("PrepareRawSpectra", Form("Wrong fitmode: %i", fitmode));
  //
  if (iMult >= nMultBin)
    Fatalmsg("PrepareRawSpectra", Form("Multiplicity bin %i out of bounds", iMult));
  //
  const Bool_t drawresult = kTRUE;                                                                 //To draw the magenta histogram of the total fit result
  const Bool_t samemismatch = iMult == 11 || iMult < 2 ? kFALSE : kFALSE;                          //To use the mismatch template from the Centrality integrated case
  const Int_t mismatchindex = 0;                                                                   //Which mismatch template to use in case we want to use the same one as for other multiplicity bins
  const Int_t useReduced = 0;                                                                      //1: reduced 2: replaced To use the mismatch template reduced in the signal contribution
  const Bool_t geteventinfo = kFALSE;                                                              //Set to take the event wide information, this shall not be mandatory as it is more for checks on normalization
  const Bool_t drawlines = kTRUE;                                                                  //Set to draw lines representing fitrange in the fitted distribution
  const Int_t filltemplates = 3001;                                                                //Set to fill the templates with shaded colors if negative nothing is done
  const Bool_t normfromlist = kFALSE;                                                              //Flag to normalize from the list
  const Bool_t showchi2 = kFALSE;                                                                  //Flag to write the chi2 of the fit on each plot
  const Bool_t showratiomean = kFALSE;                                                             //Flag to write the mean and other parameters of the ratio between the fit and the data
  const Bool_t limitrange = kTRUE;                                                                 //Flag to limit the plotting range of the th2 to the requested range
  const Bool_t showtemplateusage = kTRUE;                                                          //Flag to show the template usage for each pt bin
  const Bool_t T0fillOnly = dataname.Contains("0T0") || dataname.Contains("_T0") ? kTRUE : kFALSE; //Flag to check if the file has T0 fill only, in this case it should work as for peripheral Pb--Pb

  gROOT->SetStyle("Plain");
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

  //
  const TString speclab = pC[iCharge] + pS[iSpecies];
  const TString speclabel = pCharge[iCharge] + pSpecies[iSpecies];
  const TString spectitle = pCharge[iCharge] + " " + pSpecies[iSpecies];
  //
  TString fNameDataMismatch = "";
  if (dataname.Contains(" ")) {
    fNameDataMismatch = dataname;
    fNameDataMismatch.Replace(0, fNameDataMismatch.Index(" ") + 1, "");
    fNameDataMismatch = dirname + Form("Processed/%s", fNameDataMismatch.Data());
    while (dataname.Contains(" ") && dataname.Sizeof() > 2)
      dataname.Chop();
  }
  TString fNameData = dirname + Form("Processed/%s", dataname.Data());
  TString fNameDataEvt = dirname + "List/TListTOF.root";

  TFile* finData = GetFile(fNameData, "READ");
  TList* linData = 0x0;
  if (outputintolists)
    GetListFromFile(finData, Form("TOFList_%s", MultBinString[iMult].Data()), linData);

  TFile* finDataMismatch = 0x0;
  if (!fNameDataMismatch.EqualTo(""))
    finDataMismatch = GetFile(fNameDataMismatch, "READ");

  //Common path for macro output
  const TString outpath = Form("./Spectra/%s/Yields/%s/", systemString[optpp].Data(), speclabel.Data());
  TString fname = "Yield";
  TString identifier = pCharge[iCharge];
  identifier.Append(pSpecies[iSpecies]);
  identifier.Append("_");
  identifier.Append(fitmodes[fitmode]);
  identifier.Append("_");
  identifier.Append(MultBinString[iMult]);
  if (!prefix.IsNull() && !prefix.BeginsWith('_'))
    identifier.Append("_");
  //
  identifier.Append(prefix);
  //
  fname.Append(identifier);
  fname.Append(".root");
  fname.Prepend(outpath);
  //
  TFile* fout = nullptr;
  if (save) {
    fout = GetFile(fname, "RECREATE");
    Infomsg("PrepareRawSpectra", Form("Saving output to file %s", fout->GetName()));
  } else
    Infomsg("PrepareRawSpectra", Form("NOT Saving output. Would like to save it in %s", fname.Data()));
  //
  TFile* foutfits = nullptr;
  if (save && !fitprefix.EqualTo("")) {
    foutfits = GetFile(Form("%sFits/Fits%s_%s.root", outpath.Data(), identifier.Data(), fitprefix.Data()), "RECREATE");
    Infomsg("PrepareRawSpectra", Form("Saving fit results to file %s", foutfits->GetName()));
  }
  //
  if (fout)
    fout->cd();
  //
  TList* lHistograms = new TList();
  lHistograms->SetOwner();

  TList* lCanvas = new TList();
  lCanvas->SetOwner();

  //
  //Define all needed histograms
  //
  TH1D* hNEvt = nullptr;
  TH1D* hEvtMult = nullptr;
  TH1D* hEvtMultAftEvSel = nullptr;

  TH1F* fEntries = nullptr;                  //Event information
  TH1F* hTOF[kPtBins] = { nullptr };         //TOF distribution (time or n-sigma)
  TH1F* hTOFMismatch[kPtBins] = { nullptr }; //TOF mismatch template
  TH1F* hTOFExpected[kPtBins][kExpSpecies];  //TOF signal/background templates

  //Fitted results
  TH1F* hFitted[kPtBins] = { nullptr };         //Total fitted result
  TH1F* hMismatchFitted[kPtBins] = { nullptr }; //Mismatch fitted
  TH1F* hExpectedFitted[kPtBins][kExpSpecies];  //signal/background fitted

  //Ratio to fitted
  TH1F* hRatioToFitted[kPtBins] = { nullptr };

  //Total yield
  TH1F* hYield = nullptr;

  //Total yield residual in the +- 3 sigma range
  TH1F* hYieldResidual = nullptr;

  //Total background
  TH1F* hBackground[kExpSpecies + 1] = { nullptr };
  TH1F* hBackgroundOverlap[kExpSpecies + 1] = { nullptr };

  TH2F* hTOFPt = nullptr;
  TH2F* hTOFPtMismatch = nullptr;
  TH2F* hTOFPtSignal[kExpSpecies] = { nullptr };
  TH2F* hTOFPtFittedMismatch = nullptr;
  TH2F* hTOFPtFittedSignal[kExpSpecies] = { nullptr };
  TH2F* hTOFPtFitted = nullptr;

  //Functions
  TF1* fTOFsignal[kExpSpecies] = { nullptr };
  TF1* fTOFsignalSum = nullptr;
  TF1* fTOFbackground = nullptr;
  //
  //Track information
  //

  //Get the histograms
  TString hname = "";
//Helper macro to get
#define GetMyH(H, L, F)                               \
  {                                                   \
    if (outputintolists && L || !F)                   \
      GetHistogram(static_cast<TList*>(L), hname, H); \
    else                                              \
      GetHistogram(static_cast<TFile*>(F), hname, H); \
    H->SetDirectory(0);                               \
    lHistograms->Add(H);                              \
  }
  hname = Form("fEntries_%s", MultBinString[iMult].Data());
  GetMyH(fEntries, linData, finData);

  for (Int_t ptbin = 0; ptbin < kPtBins; ptbin++) { //Pt loop
    hname = "hTOF";
    if (fitsigma)
      hname += "Sigma";
    //
    hname += Form("%s_%i_%s", speclab.Data(), ptbin, MultBinString[iMult].Data());
    GetMyH(hTOF[ptbin], linData, finData);
    //
    hname = "hTOFMismatch";
    if (fitsigma)
      hname += "Sigma";
    //
    hname += Form("%s_%i_%s", speclab.Data(), ptbin, MultBinString[samemismatch ? mismatchindex : iMult].Data());
    if (fitsigma && (useReduced == 1 || useReduced == 2)) { //Take the coocked background
      TFile auxifile(Form("./TOFMismatch/SignalRemoved/%s/%s%i%s.root", pSpecies[iSpecies].Data(), pCharge[iCharge].Data(), ptbin, useReduced == 2 ? "Bkg" : "Sig"), "READ");
      if (auxifile.IsOpen()) {
        GetHistogram(&auxifile, useReduced == 2 ? "FitReplaced" : "FitReduced", hTOFMismatch[ptbin]);
        hTOFMismatch[ptbin]->SetName(hname);
        hTOFMismatch[ptbin]->SetDirectory(0);
        auxifile.Close();
      }
    } else {
      if (finDataMismatch) { //Get histogram for mismatch from other file
        GetMyH(hTOFMismatch[ptbin], 0x0, finDataMismatch);
      } else { //Get the mismacth from the same file
        GetMyH(hTOFMismatch[ptbin], linData, finData);
      }
    }
    //Templates for each species
    for (Int_t species = 0; species < kExpSpecies; species++) { //Species loop
      hname = "hTOF";
      if (fitsigma)
        hname += "Sigma";
      //
      hname += Form("Expected%s_%s_%i_%s", speclab.Data(), pS_all[species].Data(), ptbin, MultBinString[iMult].Data());
      GetMyH(hTOFExpected[ptbin][species], linData, finData);
    }
  }

  //
  //Event wide information
  //

  //Get the histograms for the whole event
  if (geteventinfo) {
    TFile finDataEvt(fNameDataEvt, "READ");
    TList* lin = 0x0;
    GetListFromFile(&finDataEvt, "cOutputList", lin);
    //
    hname = Form("hNEvt");
    GetMyH(hNEvt, lin, 0x0);
    //

    hname = Form("hEvtMult");
    GetMyH(hEvtMult, lin, 0x0);
    //

    hname = Form("hEvtMultAftEvSel");
    GetMyH(hEvtMultAftEvSel, lin, 0x0);
    //
    delete lin;
    finDataEvt.Close();
  }
#undef GetMyH

  //Closing input files
  if (finDataMismatch) {
    finDataMismatch->Close();
    delete finDataMismatch;
  }
  finData->Close();
  delete finData;

  //
  //Define histogram to show which template is included in the fit
  //
  TH2I* husage = 0x0;
  if (showtemplateusage) {
    husage = new TH2I("husage", Form("TemplateUsage;Template;%s", ptstring.Data()), kExpSpecies + 1, -0.5, -0.5 + kExpSpecies + 1, kPtBins, -0.5, -0.5 + kPtBins);
    for (Int_t i = 0; i < kExpSpecies; i++)
      husage->GetXaxis()->SetBinLabel(i + 1, Form("Temp. %s", speciesRoot_all[3 * i + iCharge].Data()));
    husage->GetXaxis()->SetBinLabel(kExpSpecies + 1, "Temp. Mismatch");

    for (Int_t i = 0; i < kPtBins; i++)
      husage->GetYaxis()->SetBinLabel(i + 1, Form("[%.2f,%.2f]", fBinPt[i], fBinPt[i + 1]));
  }

  //
  //Fitting section
  //
  Int_t FI = TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, kExpSpecies + 1);

  //Define the histograms
  hYield = new TH1F("hYield" + speclab, Form("Yield for %s;%s;Counts", spectitle.Data(), ptstring.Data()), kPtBins, fBinPt);
  lHistograms->AddFirst(hYield);

  hYieldResidual = new TH1F("hYieldResidual" + speclab, Form("Residual Yield for %s;%s;Counts", spectitle.Data(), ptstring.Data()), kPtBins, fBinPt);
  hYieldResidual->SetLineColor(kRed);
  hYieldResidual->SetMarkerColor(kRed);
  hYieldResidual->Sumw2();
  lHistograms->AddFirst(hYieldResidual);

  for (Int_t i = 0; i < kExpSpecies + 1; i++) {
    hBackground[i] = new TH1F(Form("hBackground%s%i", speclab.Data(), i), Form("Background for %s %i;%s;%s Yield", spectitle.Data(), i, ptstring.Data(), i == 0 ? "Mismatch" : pSpecies_all[i - 1].Data()), kPtBins, fBinPt);
    lHistograms->Add(hBackground[i]);

    hBackgroundOverlap[i] = new TH1F(Form("hBackgroundOverlap%s%i", speclab.Data(), i), Form("Background Overlap for %s %i;%s;%s Overlap Yield", spectitle.Data(), i, ptstring.Data(), i == 0 ? "Mismatch" : pSpecies_all[i - 1].Data()), kPtBins, fBinPt);
    lHistograms->Add(hBackgroundOverlap[i]);
  }

  Int_t rows = 7, columns = 7;
  TCanvas* cTOF = new TCanvas("cTOF" + speclab, "Canvas with TOF distributions for " + spectitle, screendim[0], screendim[1]);
  cTOF->Divide(rows, columns);
  lCanvas->Add(cTOF);

  TCanvas* cFit = new TCanvas("cFit" + speclab, "Fitted distribution for " + spectitle, screendim[0], screendim[1]);
  lCanvas->Add(cFit);
  cFit->Divide(rows, columns);

  TCanvas* cFitOnly = new TCanvas("cFitOnly" + speclab, "Only Fitted distribution for " + spectitle, screendim[0], screendim[1]);
  lCanvas->Add(cFitOnly);
  cFitOnly->Divide(rows, columns);

  TCanvas* cFitRatio = new TCanvas("cFitRatio" + speclab, "Fitted distribution Ratio for " + spectitle, screendim[0], screendim[1]);
  lCanvas->Add(cFitRatio);
  cFitRatio->Divide(rows, columns);

  TCanvas* cTOFPt = new TCanvas("cTOFPt" + speclab, "Canvas with Pt TOF distributions for " + spectitle, screendim[0], screendim[1]);
  lCanvas->Add(cTOFPt);

  TCanvas* cYield = new TCanvas("cYield" + speclab, "Canvas with Yield for " + spectitle);
  cYield->Divide(3);
  lCanvas->Add(cYield);

  TCanvas* cYieldResidual = new TCanvas("cYieldResidual" + speclab, "Canvas with YieldResidual for " + spectitle);
  cYieldResidual->Divide(2, 2);
  lCanvas->Add(cYieldResidual);

  TCanvas* cBackground = new TCanvas("cBackground" + speclab, "Canvas with Background for " + spectitle);
  cBackground->Divide(4, 2);
  lCanvas->Add(cBackground);

  TCanvas* cBackgroundOverlap = new TCanvas("cBackgroundOverlap" + speclab, "Canvas with BackgroundOverlap for " + spectitle);
  cBackgroundOverlap->Divide(4, 2);
  lCanvas->Add(cBackgroundOverlap);

  TCanvas* cEvent = new TCanvas("cEvent", "Canvas with Events");
  cEvent->Divide(2, 2);
  lCanvas->Add(cEvent);

  Int_t padcounter = 1;

  //Pt ranges where to fit
  //RooFit pt ranges
  const Int_t MinPtRoofit_PbPb[nSpecies] = { 10, 10, 12, 26 };
  const Int_t MaxPtRoofit_PbPb[nSpecies] = { 43, 43, 43, 59 };
  const Int_t MinPtRoofit_pp[nSpecies] = { 10, 10, 16, 26 };
  const Int_t MaxPtRoofit_pp[nSpecies] = { 43, 43, 43, 59 };
  //Functions pt ranges
  const Int_t MinPtFunctions[nSpecies] = { 10, 12, 21, 26 };
  const Int_t MaxPtFunctions[nSpecies] = { 37, 40, 40, 59 };
  //
  Int_t MinPtFit = optpp ? MinPtRoofit_pp[iSpecies] : MinPtRoofit_PbPb[iSpecies];
  Int_t MaxPtFit = optpp ? MaxPtRoofit_pp[iSpecies] : MaxPtRoofit_PbPb[iSpecies];
  for (Int_t ptbin = 0; ptbin < kPtBins; ptbin++) { //Pt loop
    Bool_t status = ptbin >= MinPtFit && ptbin <= MaxPtFit;
    Printf("###############################");
    Printf("Pt bin %i/%i: [%.2f, %.2f]", ptbin + 1, kPtBins, fBinPt[ptbin], fBinPt[ptbin + 1]);
    Printf("###############################");

    //Drawing Delta T
    Int_t colorcounter = 0;
    cTOF->cd(padcounter);
    gPad->SetLogy();

    SetColor(hTOF[ptbin], kBlack, 0);
    hTOF[ptbin]->DrawCopy();

    SetColor(hTOFMismatch[ptbin], FI + colorcounter, filltemplates);
    colorcounter++;
    hTOFMismatch[ptbin]->DrawCopy("same");

    for (Int_t species = 0; species < kExpSpecies; species++) { //Species loop
      SetColor(hTOFExpected[ptbin][species], FI + colorcounter, filltemplates);
      colorcounter++;
      hTOFExpected[ptbin][species]->DrawCopy("same");
    }

    DrawLabel(Form("pt %i [%.2f,%.2f]", ptbin, fBinPt[ptbin], fBinPt[ptbin + 1]));

    //
    //Define data and template histograms
    //
    TH1F* datahisto = hTOF[ptbin];
    TH1F* mismatchhisto = hTOFMismatch[ptbin];
    TH1F* templatehhisto[kExpSpecies] = { 0x0 };
    for (Int_t species = 0; species < kExpSpecies; species++)
      templatehhisto[species] = hTOFExpected[ptbin][species];

    //
    //Fit Parameters
    //
    const Double_t showrange[2] = { fitsigma ? -40. : -3000., fitsigma ? 40. : 3000. };            //Range where to show the actual fit in the saved figure
    const Bool_t setstart = optpp ? kTRUE : kTRUE;                                                 //Gets the first bin above zero as a fitrange
    const Bool_t setstop = optpp ? kFALSE : (iMult < 3 || iMult == nMultBin - 1) ? kFALSE : kTRUE; //Gets the last bin above zero as a fitrange

    Double_t range[2] = { fitsigma ? -300. : -10000., fitsigma ? 300. : 10000. };  //Integration range to extract yields only if using TFF
    Double_t fitrange[2] = { fitsigma ? -300. : -3000., fitsigma ? 300. : 6000. }; //Fit range
    const Int_t rebin = (fitmode == 3) ? -2 : -2;
    const Bool_t setzeroerrors = kFALSE;
    const Bool_t useMismatch = optpp ? kTRUE : kTRUE;
    const Bool_t fixMismatch = kFALSE;                       //If true Mismatch template is scaled before the fit and then fixed
    const Bool_t useOptions = kTRUE;                         //Using the pt/Centrality dependet options
    const Bool_t setTemplateWithRanges = kTRUE;              //Using the fitranges to set if or not to use the templates
    Bool_t useit[kExpSpecies] = { 0, 0, 1, 1, 1, 0 };        //Choose if use or not the template for that particular particle species
    const Bool_t usefun[kExpSpecies] = { 0, 0, 0, 0, 0, 0 }; //Choose if use or not a full function instead of the template for fitting
    Int_t ntemplates = 0;
    Int_t indexoffset = 0;
    if (useOptions && status) { //Custom pt and mult dependent parameters
      if (fitsigma) {
        if (ptbin < 25)
          useit[0] = kTRUE; //Low pt electrons
        if (ptbin < 20)
          useit[1] = kTRUE; //Low pt muons
        if (ptbin < 5)
          useit[4] = kFALSE; //low pt protons
        if (ptbin > 30)
          useit[5] = kTRUE; //high pt Deuteron
        if (fitmode == 1) { //RooFit
          if (iSpecies == 0) {
            // 	    rebin = 4;
            fitrange[0] = ptbin > 10 ? -10 : -25;
            if (ptbin >= 27)
              fitrange[0] = -5;
            if (ptbin > 40)
              fitrange[0] = -3;

            if (ptbin < 34)
              fitrange[1] = 200;
            else if (ptbin < 37)
              fitrange[1] = 100;
            else
              fitrange[1] = 40;

            if (0) {
              if (ptbin < 10)
                useit[4] = kFALSE;
            } else {

              fitrange[0] = -30;
              if (ptbin >= 29)
                fitrange[0] = -5;
              fitrange[1] = 50;

              if (ptbin < 12)
                useit[3] = kFALSE; // low pt Kaons
              if (ptbin < 24)
                useit[4] = kFALSE; // low pt Protons
              if (ptbin < 32)
                useit[5] = kFALSE; // low pt Deuterons
              if (ptbin > 39) {
                fitrange[1] = 10;
                useit[5] = kFALSE; // low pt Deuterons
              }
            }
            if (iMult >= 9 && iMult != nMultBin - 1 && ptbin >= 32)
              fitrange[1] = 20;
            if (iMult >= 9 && iMult != nMultBin - 1 && ptbin > 33)
              fitrange[1] = 15;
            if (iMult >= 9 && iMult != nMultBin - 1 && ptbin >= 35)
              fitrange[1] = 10;
            if (iMult >= 9 && iMult != nMultBin - 1 && ptbin >= 38)
              fitrange[1] = 8;
            if (iMult >= 8 && iMult != nMultBin - 1 && ptbin > 39)
              fitrange[1] = 5;
            if (iMult >= 8 && iMult != nMultBin - 1 && ptbin < 10)
              useit[3] = kFALSE; //peripheral Kaons
            if (iMult >= 9 && iMult != nMultBin - 1 && ptbin < 43)
              useit[5] = kFALSE; //peripheral Deuteron

          } else if (iSpecies == 1) {
            // 	    rebin = 5;
            fitrange[0] = ptbin > 37 ? -5 : ptbin > 33 ? -8 : ptbin > 20 ? -10 : ptbin > 10 ? -25 : -40;
            if (iCharge == 1 && ptbin > 30)
              fitrange[0] = -7;
            if (iCharge == 1 && ptbin >= 37)
              fitrange[0] = -4;
            fitrange[1] = ptbin > 20 ? 100 : 200;
            if (ptbin >= 40)
              fitrange[1] = 40;
            if (iMult >= 7 && iMult != nMultBin - 1 && ptbin > 15)
              fitrange[1] = 10;
            if (iMult >= 7 && iMult != nMultBin - 1 && ptbin > 19)
              fitrange[1] = 5;
            if (optpp) {
              fitrange[0] = -70;
              fitrange[1] = 50;
            }
            useit[0] = kFALSE; //Low pt electrons
            useit[1] = kFALSE; //Low pt muons
            // 	if(ptbin < 15) useit[2] = kFALSE;//low pt pions
            if (ptbin < 5)
              useit[4] = kFALSE; //low pt protons
            if (ptbin <= 34)
              useit[5] = kFALSE; //high pt Deuteron
            if (iCharge == 1 && ptbin < 37)
              useit[5] = kFALSE; //high pt Deuteron
            if (iMult >= 8 && iMult != nMultBin - 1 && ptbin < 26)
              useit[4] = kFALSE; //peripheral Protons
            if (iMult >= 7 && iMult != nMultBin - 1 && ptbin > 15 && ptbin < 35)
              useit[4] = kFALSE; //peripheral Protons
            if (iMult >= 7 && iMult != nMultBin - 1 && ptbin < 43)
              useit[5] = kFALSE; //peripheral Deuteron
            if (iMult == 9 && ptbin < 15) {
              // useit[2] = kFALSE;//peripheral Deuteron
              fitrange[1] = 5; //peripheral Deuteron
            }
            if (optpp)
              useit[5] = kFALSE; //high pt Deuteron
          } else if (iSpecies == 2) {
            // 	    rebin = 5;

            //
            //Setting fit range
            //
            if (0) {
              fitrange[0] = -50;
              if (ptbin >= 15)
                fitrange[0] = -40;
              if (ptbin >= 20)
                fitrange[0] = -30;
              if (ptbin >= 25)
                fitrange[0] = -25;
              if (ptbin >= 28)
                fitrange[0] = -10;
              if (ptbin >= 40)
                fitrange[0] = -5;
            } else {
              fitrange[0] = -30;
              fitrange[0] = templatehhisto[0]->GetXaxis()->GetBinCenter(templatehhisto[0]->GetMaximumBin());
            }
            fitrange[1] = ptbin > 35 ? 30 : 150;
            // fitrange[0] = -15;
            // fitrange[1] = 15;

            fitrange[1] = 40;
            fitrange[0] = -20;
            fitrange[1] = 50;
            if (iMult >= 8 && iMult != nMultBin - 1 && ptbin >= 35)
              fitrange[1] = 6;
            if (iMult >= 9 && iMult != nMultBin - 1 && ptbin > 40)
              fitrange[1] = 4;

            //
            //Setting templates
            //
            //Unused templates
            useit[0] = kFALSE; //Low pt electrons
            useit[1] = kFALSE; //Low pt muons

            //Concurrent templates
            // // 	if(ptbin < 25) useit[3] = kFALSE;
            // if(ptbin < 30) useit[2] = kFALSE;//high pt Deuteron
            // if(ptbin < 30) useit[3] = kFALSE;//high pt Deuteron
            // if(ptbin > 40) useit[5] = kTRUE;//high pt Deuteron
            // else useit[5] = kFALSE;//high pt Deuteron
            // if(/*iCharge == 0 && */iMult >= 6 && iMult != nMultBin -1 && ptbin < 43) useit[5] = kFALSE;//peripheral Deuteron

          } else if (iSpecies == 3) {
            // 	    rebin = 5;
            fitrange[0] = -50;
            if (ptbin >= 15)
              fitrange[0] = -40;
            if (ptbin >= 20)
              fitrange[0] = -30;
            if (ptbin >= 25)
              fitrange[0] = -25;
            if (ptbin >= 28)
              fitrange[0] = -10;
            if (ptbin >= 40)
              fitrange[0] = -5;
            fitrange[1] = ptbin > 35 ? 30 : 150;
            useit[5] = kTRUE; //Include Deuterons!!
          }
        } else if (fitmode == 3) { //Functions
          for (Int_t species = 0; species < kExpSpecies; species++)
            useit[species] = (species == kpi || species == kK || species == kp) ? kTRUE : kFALSE;
          if (ptbin < 20)
            useit[kp] = kFALSE;
        }

      } else {
        if (fitmode == 1) { //RooFit
          useit[kpi] = kTRUE;
          switch (iSpecies) {
          case 0: //Pions
          {
            if (ptbin < 12)
              useit[ke] = kFALSE;
            useit[kpi] = kTRUE;
            if (ptbin < 12)
              useit[kK] = kFALSE;
            if (ptbin < 15)
              useit[kp] = kFALSE;
            else
              useit[kp] = kTRUE;
            useit[4] = kFALSE;
            useit[kp] = kTRUE;
            break;
          }
          case 1: //Kaons
          {
            useit[kK] = kTRUE;

            useit[ke] = kFALSE;
            useit[kmu] = kFALSE;
            if (ptbin < 15)
              useit[kpi] = kFALSE;
            else
              useit[kpi] = kTRUE;

            if (ptbin < 20)
              useit[kp] = kFALSE;
            else
              useit[kp] = kTRUE;

            if (ptbin < 15) {
              fitrange[0] = -1000;
              fitrange[1] = 1000;
            }
            break;
          }
          case 2: //Protons
          {
            if (ptbin < 20) {
              useit[2] = kFALSE;
              useit[3] = kFALSE;
            }
            break;
          }
          }

          if (ptbin > 10)
            useit[3] = kTRUE;
          // if(ptbin > 12) useit[4] = kTRUE;
        } else if (fitmode == 3) { //Functions
          for (Int_t species = 0; species < kExpSpecies; species++)
            useit[species] = (species == kpi || species == kK || species == kp) ? kTRUE : kFALSE;
          if (ptbin < 16)
            useit[kK] = kFALSE;
          if (ptbin < 20)
            useit[kp] = kFALSE;
        }
      }
    }
    if (setTemplateWithRanges) {
      for (Int_t i = 0; i < kExpSpecies; i++) {
        if (useit[i] && !IsHistogramInRange(templatehhisto[i], fitrange[0], fitrange[1], /*threshold*/ 1., /*verbose*/ kFALSE)) {
          useit[i] = kFALSE;
          cout << "Template " << templatehhisto[i]->GetName() << " is not in range" << endl;
        } else if (useit[i])
          cout << "Template " << templatehhisto[i]->GetName() << " is in range" << endl;
      }
    }

    //
    //Signal definition
    //
    if (fitmode == 3 && status) { //Fit with functions
      for (Int_t species = 0; species < kExpSpecies; species++) {
        fTOFsignal[species] = new TF1(Form("expTOFsignal%s%s", pCharge[iCharge].Data(), pSpecies_all[species].Data()), TOFsignal, range[0], range[1], 4);
        fTOFsignal[species]->SetLineColor(templatehhisto[species]->GetLineColor());
        //Normalization
        fTOFsignal[species]->SetParameter(0, templatehhisto[species]->GetMaximum());
        fTOFsignal[species]->SetParLimits(0, 0, (Double_t)datahisto->GetEntries());
        //Mean
        fTOFsignal[species]->SetParameter(1, templatehhisto[species]->GetMean());
        fTOFsignal[species]->SetParLimits(1, templatehhisto[species]->GetMean() - 2. * templatehhisto[species]->GetRMS(), templatehhisto[species]->GetMean() + 2. * templatehhisto[species]->GetRMS());
        //Sigma
        fTOFsignal[species]->SetParameter(2, fitsigma ? 1. : tofReso);
        fTOFsignal[species]->SetParLimits(2, fitsigma ? 0.5 : 40., fitsigma ? 2.0 : 300.);
        //Tail
        fTOFsignal[species]->SetParameter(3, fitsigma ? 1. : tofTail);
        fTOFsignal[species]->SetParLimits(3, fitsigma ? 0.5 : 40., fitsigma ? 2.0 : 200.);
        fTOFsignal[species]->SetParNames("Norm", "Mean", "Sigma", "Tail");
      }

      Int_t usable = 0;
      for (Int_t species = 0; species < kExpSpecies; species++) {
        cout << "Species " << species << " : " << useit[species] << endl;
        if (useit[species])
          usable++;
      }

      if (useMismatch) {
        if (usable == 1)
          fTOFsignalSum = new TF1(Form("expTOFsignalSum"), TOFsignalBkg, range[0], range[1], 4 + 1);
        else if (usable == 2)
          fTOFsignalSum = new TF1(Form("expTOFsignalSum"), TOFsignalBkg_double, range[0], range[1], 8 + 1);
        else if (usable == 3)
          fTOFsignalSum = new TF1(Form("expTOFsignalSum"), TOFsignalBkg_triple, range[0], range[1], 12 + 1);
        else
          Fatalmsg("PrepareRawSpectra", Form("Too many (%i) requested functions!", usable));

        fTOFbackground = new TF1(Form("expTOFbackground"), "pol0", range[0], range[1]);
        fTOFbackground->SetParameter(0, 0);
        fTOFbackground->SetParLimits(0, 0, datahisto->GetEffectiveEntries());
        fTOFbackground->SetParName(0, "Bkg");

      } else {
        if (usable == 1)
          fTOFsignalSum = new TF1(Form("expTOFsignalSum"), TOFsignal, range[0], range[1], 4);
        else if (usable == 2)
          fTOFsignalSum = new TF1(Form("expTOFsignalSum"), TOFsignal_double, range[0], range[1], 8);
        else if (usable == 3)
          fTOFsignalSum = new TF1(Form("expTOFsignalSum"), TOFsignal_triple, range[0], range[1], 12);
        else
          Fatalmsg("PrepareRawSpectra", Form("Too many (%i) requested functions!", usable));
      }

      Int_t counter = 0;
      for (Int_t species = 0; species < kExpSpecies; species++) {
        if (!useit[species])
          continue;
        for (Int_t i = 0; i < fTOFsignal[species]->GetNpar(); i++) {
          cout << "For parameter " << counter << " using " << species << endl;

          Double_t l[2];
          fTOFsignal[species]->GetParLimits(i, l[0], l[1]);
          fTOFsignalSum->SetParameter(counter, fTOFsignal[species]->GetParameter(i));
          fTOFsignalSum->SetParLimits(counter, l[0], l[1]);
          fTOFsignalSum->SetParName(counter, fTOFsignal[species]->GetParName(i));

          counter++;
        }
      }
      if (useMismatch) {
        for (Int_t i = 0; i < fTOFbackground->GetNpar(); i++) {
          Double_t l[2];
          fTOFbackground->GetParLimits(i, l[0], l[1]);
          fTOFsignalSum->SetParameter(counter, fTOFbackground->GetParameter(i));
          fTOFsignalSum->SetParLimits(counter, l[0], l[1]);
          fTOFsignalSum->SetParName(counter, fTOFbackground->GetParName(i));

          counter++;
        }
      }

      if (counter != fTOFsignalSum->GetNpar())
        Fatalmsg("PrepareRawSpectra", "Problem in the initial conditions for conv. fit funct.");
      fTOFsignalSum->SetLineColor(kMagenta);
    }

    //Checking the templates before running the fit
    const Bool_t checkrange = kFALSE;
    Bool_t inrange[kExpSpecies] = { 0 };

    //To choose the which templates are in range it is important that they all have the same bin limits, the selection is made with overflow and underflow fraction
    if (checkrange) {
      Double_t templateentries[kExpSpecies] = { -999 };
      Double_t templateoverflow[kExpSpecies] = { -999 };
      Double_t templateunderflow[kExpSpecies] = { -999 };
      Double_t templatemean[kExpSpecies] = { -999 };
      Double_t templatesigma[kExpSpecies] = { -999 };

      for (Int_t species = 0; species < kExpSpecies; species++) { //Species loop
        templateentries[species] = templatehhisto[species]->GetEntries();
        templateoverflow[species] = templatehhisto[species]->GetBinContent(templatehhisto[species]->GetNbinsX() + 1);
        templateunderflow[species] = templatehhisto[species]->GetBinContent(0);
        templatemean[species] = templatehhisto[species]->GetMean();
        templatesigma[species] = templatehhisto[species]->GetRMS();

        if (templateentries[species] == 0)
          continue;
        Double_t overflowfraction = templateoverflow[species] / templateentries[species];
        Double_t underflowfraction = templateunderflow[species] / templateentries[species];
        Double_t diffplus = templatemean[species] + templatesigma[species];
        Double_t diffminus = templatemean[species] - templatesigma[species];
        Printf("AliPID %i mean %f rms %f difference %f", species, templatemean[species], templatesigma[species], diffplus);
        Printf("\t\t entries %f overflow %f (frac. %f) underflow %f (frac. %f)", templateentries[species], templateoverflow[species], overflowfraction, templateunderflow[species], underflowfraction);

        if (overflowfraction < 0.99 && underflowfraction < 0.99)
          inrange[species] = kTRUE;
        else
          inrange[species] = kFALSE;

        if (useit[species] && !inrange[species])
          useit[species] = kFALSE;
      }
    }

    //
    //Count the templates used
    //
    if (useMismatch) { //Mismatch
      if (husage)
        husage->Fill(kExpSpecies, ptbin);
      indexoffset++;
      ntemplates++;
      Printf("Using Mismatch template: %i/%i", indexoffset, ntemplates);
    }
    for (UInt_t species = 0; species < kExpSpecies; species++) { //Signal for particles
      if (useit[species]) {
        if (husage)
          husage->Fill(species, ptbin);
        ntemplates++;
        if (species < iSpecies + kpi)
          indexoffset++;
        Printf("Using template %i %s: %i/%i", species, pSpecies_all[species].Data(), indexoffset, ntemplates);
      }
    }

    if (rebin > 0)
      datahisto->Rebin(rebin);

    TObjArray* templates = new TObjArray(ntemplates);
    TObjArray* prediction = new TObjArray(ntemplates + 1); //Has also the total fitted template
    if (fitmode == 3) {                                    //Fit with functions only
      for (Int_t species = 0; species < kExpSpecies; species++) {
        if (useit[species])
          templates->Add(fTOFsignal[species]);
      }
      if (useMismatch)
        templates->Add(fTOFbackground);
    } else { //Fit with templates
      if (useMismatch) {
        if (rebin > 0)
          mismatchhisto->Rebin(rebin);
        if (setzeroerrors)
          for (Int_t i = 1; i <= mismatchhisto->GetNbinsX(); i++)
            mismatchhisto->SetBinError(i, 0);
        templates->Add(mismatchhisto);
      }
      for (Int_t species = 0; species < kExpSpecies; species++) {
        if (rebin > 0)
          templatehhisto[species]->Rebin(rebin);
        if (!useit[species])
          continue;
        if (usefun[species]) {

          // TF1 *funct = new TF1(Form("signal%s%s", pC[iCharge].Data(), pS_all[species].Data()), TOFsignalNorm, range[0], range[1], 3);
          // funct->SetTitle(Form("%s", speciesRoot_all[3*species + iCharge].Data()));
          const Double_t m = templatehhisto[species]->GetXaxis()->GetBinCenter(templatehhisto[species]->GetMaximumBin()); //templatehhisto[species]->GetMean() > 0 ? templatehhisto[species]->GetMean() : 0;
          const Double_t mLow = fitsigma ? m - 3 : m - 1;
          const Double_t mUp = fitsigma ? m + 3 : m + 1;
          const Double_t s = fitsigma ? 1 : 56;
          const Double_t sLow = fitsigma ? s - .5 : s - 1;
          const Double_t sUp = fitsigma ? s + .5 : s + 1;
          const Double_t t = 1.25;
          const Double_t tLow = t - 1;
          const Double_t tUp = t + .01;

          //
          RooRealVar* x = new RooRealVar("x", "x" + speciesRoot_all[3 * species + iCharge],
              datahisto->GetXaxis()->GetBinLowEdge(1),
              datahisto->GetXaxis()->GetBinUpEdge(datahisto->GetNbinsX()));
          RooRealVar* mean = new RooRealVar(
              "mean" + speciesRoot_all[3 * species + iCharge],
              "mean" + speciesRoot_all[3 * species + iCharge],
              m);
          // m, mLow, mUp);
          RooRealVar* sigma = new RooRealVar(
              "sigma" + speciesRoot_all[3 * species + iCharge],
              "sigma" + speciesRoot_all[3 * species + iCharge],
              // s);
              s, s * .95, 400);
          // s, sLow, sUp);
          cout << "RMS IS " << templatehhisto[species]->GetRMS() << endl;
          RooRealVar* qL = new RooRealVar(
              "qL" + speciesRoot_all[3 * species + iCharge],
              "qL" + speciesRoot_all[3 * species + iCharge],
              1.0);
          // 1.0, 0.9999, 1.0001);
          RooRealVar* qR = new RooRealVar(
              "qR" + speciesRoot_all[3 * species + iCharge],
              "qR" + speciesRoot_all[3 * species + iCharge],
              1.667);
          // 1.667, 1.5, 1.7);
          RooqGaussian* funct = new RooqGaussian("qGaussian" + speciesRoot_all[3 * species + iCharge], speciesRoot_all[3 * species + iCharge], *x, *mean, *sigma, *qR, *qL);
          // funct->SetParameters(m, s, t);
          // funct->SetParLimits(0, mLow, mUp);
          // funct->SetParLimits(1, sLow, sUp);
          // funct->SetParLimits(2, tLow, tUp);
          // funct->SetParNames("Mean", "Sigma", "Tail");

          templates->Add(funct);

          // new TCanvas();
          // RooPlot *plot = x.frame();
          // funct->plotOn(plot);
          // plot->Draw();
          // lHistograms->Add(funct);
          // return;
        } else
          templates->Add(templatehhisto[species]);
      }
    }

    //
    //Set fit ranges
    //
    if (useOptions && status) { //Auto-Set start/stop for fit
      //Sets the first bin above a certain value as the start for the fit
      if (setstart) {
        const Double_t mincounts = 10;
        const Int_t shift = 0;
        const Int_t databin = GetHistoLowRange(datahisto, mincounts);
        if (fitmode == 3) { //With functions only data is important
          fitrange[0] = datahisto->GetXaxis()->GetBinLowEdge(databin);
          Printf("Autoset start for function fit to bin %i [%f, %f]", databin, datahisto->GetXaxis()->GetBinLowEdge(databin), datahisto->GetXaxis()->GetBinUpEdge(databin));
        } else {
          const Int_t firstindex = GetFirstHistogram(templates, mincounts); //First template from the left with non zero values
          TH1F* firsthisto = dynamic_cast<TH1F*>(templates->At(firstindex));
          const Int_t templatebin = GetHistoLowRange(firsthisto, mincounts);
          const Int_t templatebinAfter = GetHistoLowRangeAfter(firsthisto, mincounts, templatebin);
          const Int_t templatebinNoHoles = GetHistoNoHolesAfter(firsthisto, templatebin);
          if (templatebin > templatebinAfter)
            Fatalmsg("PrepareRawSpectra", "Wrongly defined bins");
          if (!SameBinning(datahisto, firsthisto))
            Fatalmsg("PrepareRawSpectra", Form("Histograms have different binning : datahisto %s (%i) and the %s (%i)", datahisto->GetName(), datahisto->GetNbinsX(), firsthisto->GetName(), firsthisto->GetNbinsX()));

          const Int_t binlimits[2] = { databin, templatebinNoHoles };
          const TString bins[2] = {
            Form("databin %i [%f, %f] -> %f", binlimits[0], datahisto->GetXaxis()->GetBinLowEdge(binlimits[0]), datahisto->GetXaxis()->GetBinUpEdge(binlimits[0]), datahisto->GetBinContent(binlimits[0])),
            Form("templates %i [%f, %f] -> %f", binlimits[1], datahisto->GetXaxis()->GetBinLowEdge(binlimits[1]), datahisto->GetXaxis()->GetBinUpEdge(binlimits[1]), datahisto->GetBinContent(binlimits[1]))
          };
          TString type = "";

          cout << "First histogram with required counts is " << firstindex << " " << templates->At(firstindex)->GetName() << endl;
          cout << "Setting start to first bin above " << mincounts << " " << bins[0] << ", " << bins[1] << endl;
          Double_t x = -999;
          if (binlimits[0] > 0 && binlimits[1] < binlimits[0]) {
            x = datahisto->GetXaxis()->GetBinLowEdge(binlimits[0] + shift);
            type = "from templates";
          } else if (binlimits[0] > 0) {
            x = datahisto->GetXaxis()->GetBinLowEdge(binlimits[1] + shift);
            type = "from data";
          }

          Bool_t inrange = kTRUE;
          const Double_t min = iSpecies == 0 ? 2.2 : 2.;
          if (fitsigma) {
            if (x != -999 && x > -min) {
              inrange = kFALSE;
              fitrange[0] = -min;
            }
          } else {
            if (x != -999 && x > -min * 80) {
              inrange = kFALSE;
              fitrange[0] = -min * 80;
            }
          }

          if (inrange && x != -999 && fitrange[0] < x)
            fitrange[0] = x;

          // if(x != -999 && fitrange[0] < x) fitrange[0] = x;

          const Int_t bin = datahisto->GetXaxis()->FindBin(fitrange[0]);
          Infomsg("PrepareRawSpectra", Form("Required to autoset the start point %s for the fit from %f (%f) bin %i content %f", type.Data(), fitrange[0], x, bin, datahisto->GetBinContent(bin)));
        }
      }

      //Sets the last bin above a certain value as the stop for the fit
      if (setstop) {
        const Double_t mincounts = 5;
        const Int_t databin = datahisto->FindLastBinAbove(mincounts);
        Double_t x = datahisto->GetXaxis()->GetBinUpEdge(databin);
        x = GetHistoNoHolesBefore(datahisto, x);

        cout << "Setting stop to last bin above " << mincounts << endl;
        Bool_t inrange = kTRUE;
        if (fitsigma) {
          if (x > 0 && x < 5.) {
            inrange = kFALSE;
            fitrange[1] = 5.;
          }
        } else {
          if (x > 0 && x < 5. * 80) {
            inrange = kFALSE;
            fitrange[1] = 30. * 80;
          }
        }

        if (inrange && databin > 0 && fitrange[1] > x)
          fitrange[1] = x;

        const Int_t bin = datahisto->GetXaxis()->FindBin(fitrange[1]);
        Infomsg("PrepareRawSpectra", Form("Required to autoset the stop point for the fit to %f bin %i content %f", fitrange[1], bin, datahisto->GetBinContent(bin)));
      }

      //Check the definition of the fit ranges
      if (fitrange[1] < fitrange[0]) {
        TCanvas* RangeAttempt = new TCanvas("RangeAttempt", "Attempt");
        datahisto->DrawCopy();
        TLine* l[2];
        for (Int_t counter = 0; counter < 2; counter++) {
          l[counter] = new TLine(fitrange[counter], 0, fitrange[counter], datahisto->GetMaximum());
          l[counter]->Draw();
        }
        DrawLabel(Form("Wrong range defined [%f,%f]", fitrange[0], fitrange[1]));
        RangeAttempt->SaveAs("./Images/RangeAttempt.pdf");
        Fatalmsg("PrepareRawSpectra", Form("Fit range auto defined but wrongly.. %f < %f and should be the opposite", fitrange[1], fitrange[0]));
      }
    }

    TArrayD fraction(ntemplates);
    TArrayD fractionErr(ntemplates);
    for (Int_t i = 0; i < ntemplates; i++) {
      fraction[i] = 0;
      fractionErr[i] = 0;
      if (fitmode == 1) { //fraction are used to set the initial conditions
        if (useMismatch && i == 0) {
          fraction[i] = datahisto->GetEffectiveEntries() * 0.01;
          if (fixMismatch) {
            Int_t binlow = datahisto->GetXaxis()->FindBin(150);
            Int_t binup = datahisto->GetXaxis()->FindBin(210);
            Double_t tempintegral = mismatchhisto->Integral(fitrange[0], fitrange[1]);
            fraction[i] = -tempintegral * datahisto->Integral(binlow, binup) / (datahisto->GetXaxis()->GetBinUpEdge(210) - datahisto->GetXaxis()->GetBinLowEdge(150));
          }

        } else
          fraction[i] = datahisto->GetEffectiveEntries();
      }
    }

    cFit->cd(padcounter);

    //
    //Fitting
    //
    Double_t chi2 = showchi2 ? indexoffset : -10;
    // Double_t chi2 = -1;

    if (status == kTRUE) {
      switch (fitmode) {
      case 0:
        status = PerformFitWithTFF(datahisto, templates, range, fitrange, fraction, fractionErr, prediction);
        break;
      case 1:
        status = PerformFitWithRooFit(datahisto, templates, range, fitrange, fraction, fractionErr, prediction, chi2);
        break;
      case 2:
        status = PerformFitWithCD(datahisto, templates, range, fitrange, fraction, fractionErr, prediction);
        break;
      case 3:
        status = PerformFitWithFunctions(datahisto, templates, fTOFsignalSum, range, fitrange, fraction, fractionErr, prediction);
        break;
      default:
        Fatalmsg("PrepareRawSpectra", Form("The fitmode selected (%i) is not available", fitmode));
        break;
      }
    } else
      Infomsg("PrepareRawSpectra", "Required not to fit!");
    gPad->SetLogy();

    //Drawing labels on the pad
    TString label = Form("%s %i [%.2f,%.2f]", ptstring.Data(), ptbin, fBinPt[ptbin], fBinPt[ptbin + 1]);
    if (!status)
      label += "Fit not performed!";
    DrawLabel(label);
    if (showchi2)
      DrawLabel(Form("Chi2 = %f", chi2), 0.14, 0.82, 0.33, 0.92);

    //Saving pad to file
    SavePad(Form("%sFits/%i_%s%s.pdf", outpath.Data(), ptbin, MultBinString[iMult].Data(), fitmodes[fitmode].Data()), showrange[0], showrange[1], 1, -1, kTRUE);

    //
    //Get the yield from the fit
    //
    const Double_t ptmean = (fBinPt[ptbin + 1] + fBinPt[ptbin]) / 2.;
    const Int_t ptmeanindex = hYield->GetXaxis()->FindBin(ptmean);
    const Double_t Yield = ntemplates > 0 ? fraction[indexoffset] : 0;
    const Double_t YieldError = ntemplates > 0 ? fractionErr[indexoffset] : 0;
    const Double_t YieldResidual = status ? GetResidualYield(datahisto, (TH1*)prediction->At(0), -3, 3) : 0;

    //Yield histograms
    if (ptmeanindex != ptbin + 1)
      Fatalmsg("PrepareRawSpectra", Form("Wrong index %i compare to %i", ptmeanindex, ptbin + 1));
    hYield->SetBinContent(ptmeanindex, status ? Yield : -10);
    hYield->SetBinError(ptmeanindex, status ? YieldError : 0);
    hYieldResidual->SetBinContent(ptmeanindex, YieldResidual);

    Int_t fractionindex = 0;
    for (Int_t index = 0; index < kExpSpecies + 1; index++) {
      hBackground[index]->SetBinContent(ptmeanindex, -10);
      hBackground[index]->SetBinError(ptmeanindex, 0);
      if (index == 0 && !useMismatch)
        continue;
      else if (index >= 1 && !useit[index - 1])
        continue;
      hBackground[index]->SetBinContent(ptmeanindex, fraction[fractionindex]);
      hBackground[index]->SetBinError(ptmeanindex, fractionErr[fractionindex]);
      fractionindex++;
    }

    //Retrieve fit histograms
    if (status) {
      Infomsg("PrepareRawSpectra", Form("Fit went OK with %i templates", ntemplates));
      cout << "In output:" << endl;
      prediction->ls();
      cout << "In input:" << endl;
      templates->ls();
      Int_t counter = 0;
      hFitted[ptbin] = (TH1F*)prediction->At(counter++);
      if (useMismatch)
        hMismatchFitted[ptbin] = (TH1F*)prediction->At(counter++);
      else {
        hMismatchFitted[ptbin] = (TH1F*)(mismatchhisto->Clone(Form("unfittedprediction%i", counter)));
        hMismatchFitted[ptbin]->Reset();
      }

      for (Int_t species = 0; species < kExpSpecies; species++) {
        if (useit[species])
          hExpectedFitted[ptbin][species] = (TH1F*)prediction->At(counter++);
        else {
          hExpectedFitted[ptbin][species] = (TH1F*)(templatehhisto[species]->Clone(Form("unfittedprediction%i", counter)));
          hExpectedFitted[ptbin][species]->Reset();
        }
      }

    } else {
      Int_t counter = 0;
      hFitted[ptbin] = (TH1F*)datahisto->Clone("results");
      hFitted[ptbin]->Reset();
      hMismatchFitted[ptbin] = (TH1F*)(mismatchhisto->Clone(Form("prediction%i", counter++)));
      hMismatchFitted[ptbin]->Reset();
      for (Int_t species = 0; species < kExpSpecies; species++) {
        hExpectedFitted[ptbin][species] = (TH1F*)(templatehhisto[species]->Clone(Form("prediction%i", counter++)));
        hExpectedFitted[ptbin][species]->Reset();
      }
    }

    TH1F* hlines[2] = { 0x0 };
    if (drawlines) { //Histograms with the fit range
      hlines[0] = (TH1F*)datahisto->Clone(Form("fitranges%i_%i_%i_%i", iCharge, iSpecies, iMult, ptbin));
      hlines[0]->SetDirectory(0);
      hlines[0]->Reset();
      hlines[0]->SetTitle("Fit range");
      hlines[0]->SetLineColor(kCyan);
      hlines[0]->Fill(fitrange[0], datahisto->GetMaximum());
      hlines[0]->Fill(fitrange[1], datahisto->GetMaximum());

      hlines[1] = (TH1F*)hlines[0]->Clone(Form("ranges%i_%i_%i_%i", iCharge, iSpecies, iMult, ptbin));
      hlines[1]->SetDirectory(0);
      hlines[1]->Reset();
      hlines[1]->SetTitle("Integration range");
      hlines[1]->SetLineColor(kOrange);
      hlines[1]->Fill(range[0], datahisto->GetMaximum());
      hlines[1]->Fill(range[1], datahisto->GetMaximum());
    }

    //
    //Finally drawing the result of the fit
    //
    cFitOnly->cd(padcounter);
    gPad->SetLogy();
    datahisto->SetTitle("Data");
    datahisto->DrawCopy()->GetXaxis()->SetRangeUser(showrange[0], showrange[1]);
    hFitted[ptbin]->SetLineColor(kMagenta);
    hFitted[ptbin]->SetMarkerColor(kMagenta);
    if (filltemplates > 0) {
      hFitted[ptbin]->SetFillStyle(filltemplates);
      hFitted[ptbin]->SetFillColor(kMagenta);
    }
    hFitted[ptbin]->DrawCopy("same");
    hMismatchFitted[ptbin]->DrawCopy("same");
    Double_t bkgerror = 0;
    hBackgroundOverlap[0]->SetBinContent(ptmeanindex, GetOverlapFraction(hExpectedFitted[ptbin][iSpecies + 2], hMismatchFitted[ptbin], bkgerror, ptbin == 20 ? kFALSE : kFALSE));
    hBackgroundOverlap[0]->SetBinError(ptmeanindex, bkgerror);

    for (Int_t species = 0; species < kExpSpecies; species++) {
      hBackgroundOverlap[species + 1]->SetBinContent(ptmeanindex, GetOverlapFraction(hExpectedFitted[ptbin][iSpecies + 2], hExpectedFitted[ptbin][species], bkgerror));
      hBackgroundOverlap[species + 1]->SetBinError(ptmeanindex, bkgerror);
      if (useit[species])
        hExpectedFitted[ptbin][species]->DrawCopy("same");
    }
    if (drawlines) { //Histograms with the fit range
      hlines[0]->DrawCopy("same");
      hlines[1]->DrawCopy("same");
    }
    DrawLabel(label);
    if (showchi2)
      DrawLabel(Form("Chi2 = %f", chi2), 0.14, 0.82, 0.33, 0.92);

    //Saving pads to figures
    SavePad(Form("%sFits/FITONLY_%i_%s%s.pdf", outpath.Data(), ptbin, MultBinString[iMult].Data(), fitmodes[fitmode].Data()), showrange[0], showrange[1], 1, -1, kTRUE, Form("Fit with %s", fitmodes[fitmode].Data()), kTRUE);
    gPad->SetLogy(kFALSE);
    SavePad(Form("%sFits/FITONLY_linear_%i_%s%s.pdf", outpath.Data(), ptbin, MultBinString[iMult].Data(), fitmodes[fitmode].Data()), showrange[0], showrange[1], 1, -1, kTRUE, Form("Fit with %s", fitmodes[fitmode].Data()), kTRUE);
    gPad->SetLogy();

    //Drawing ratio between the fit and the data
    cFitRatio->cd(padcounter);
    hRatioToFitted[ptbin] = (TH1F*)hFitted[ptbin]->Clone(Form("RatioToFitted%s%s_%i", pS[iSpecies].Data(), pS[iSpecies].Data(), ptbin));
    hRatioToFitted[ptbin]->SetTitle(Form("Ratio to Fit;%s;Fitted/Data", fitsigma ? nsigmastringSpecies[iSpecies + kpi].Data() : tofsignalstringSpecies[iSpecies + kpi].Data()));
    hRatioToFitted[ptbin]->Sumw2();
    hRatioToFitted[ptbin]->Divide(hFitted[ptbin], datahisto, 1, 1);
    hRatioToFitted[ptbin]->DrawCopy()->GetXaxis()->SetRangeUser(fitrange[0], fitrange[1]);
    ;
    DrawLabel(label);
    if (showratiomean) {
      TF1* meanvalue = new TF1("fmean", "pol0", -100, 100);
      meanvalue->SetParameter(0, 0);
      meanvalue->SetParLimits(0, 0, 100);
      hRatioToFitted[ptbin]->Fit(meanvalue, "0 N R M S");
      meanvalue->DrawCopy("same");
      DrawLabel(Form("Mean = %f", meanvalue->GetParameter(0)), 0.14, 0.82, 0.33, 0.92);
    }

    padcounter++;
  }

  cTOFPt->cd();
  //2D histograms
  TH1F* auxiarray[kPtBins];

  for (Int_t ptbin = 0; ptbin < kPtBins; ptbin++)
    auxiarray[ptbin] = hTOF[ptbin];
  hTOFPt = Form2DHisto(auxiarray, kPtBins, fBinPt, Form("hTOFPt%s%s", pC[iCharge].Data(), pS[iSpecies].Data()), Form("hTOFPt%s%s;%s;%s", pC[iCharge].Data(), pS[iSpecies].Data(), fitsigma ? nsigmastringSpecies[iSpecies + kpi].Data() : tofsignalstringSpecies[iSpecies + kpi].Data(), ptstring.Data()));
  CopyStyle(hTOFPt, auxiarray[0]);

  for (Int_t ptbin = 0; ptbin < kPtBins; ptbin++)
    auxiarray[ptbin] = hMismatchFitted[ptbin];
  hTOFPtFittedMismatch = Form2DHisto(auxiarray, kPtBins, fBinPt, Form("hTOFPtFittedMismatch%s%s", pC[iCharge].Data(), pS[iSpecies].Data()), Form("hTOFPtFittedMismatch%s%s", pC[iCharge].Data(), pS[iSpecies].Data()));
  CopyStyle(hTOFPtFittedMismatch, auxiarray[0]);

  for (Int_t species = 0; species < kExpSpecies; species++) {
    for (Int_t ptbin = 0; ptbin < kPtBins; ptbin++)
      auxiarray[ptbin] = hExpectedFitted[ptbin][species];
    hTOFPtFittedSignal[species] = Form2DHisto(auxiarray, kPtBins, fBinPt, Form("hTOFPtFittedSignal%s%s%s", pC[iCharge].Data(), pS[iSpecies].Data(), pS_all[species].Data()), Form("hTOFPtFittedSignal%s%s%s", pC[iCharge].Data(), pS[iSpecies].Data(), pS_all[species].Data()));
    CopyStyle(hTOFPtFittedSignal[species], auxiarray[0]);
  }

  for (Int_t ptbin = 0; ptbin < kPtBins; ptbin++)
    auxiarray[ptbin] = hFitted[ptbin];
  hTOFPtFitted = Form2DHisto(auxiarray, kPtBins, fBinPt, Form("hTOFPtFitted%s%s", pC[iCharge].Data(), pS[iSpecies].Data()), Form("hTOFPtFitted%s%s", pC[iCharge].Data(), pS[iSpecies].Data()));
  CopyStyle(hTOFPtFitted, auxiarray[0]);

  if (limitrange) {
    hTOFPt->GetXaxis()->SetRangeUser(fitsigma ? -10 : -500, fitsigma ? 10 : 500);
    hTOFPt->GetYaxis()->SetRangeUser(PtRange[optpp][2 * iSpecies], PtRange[optpp][2 * iSpecies + 1]);
  }
  hTOFPt->DrawCopy();
  hTOFPtFittedMismatch->DrawCopy("same");
  for (Int_t species = 0; species < kExpSpecies; species++)
    hTOFPtFittedSignal[species]->DrawCopy("same");
  if (drawresult)
    hTOFPtFitted->DrawCopy("same");
  cTOFPt->AddExec("dynamic", "DynamicExec()");

  //
  //Computing normalization factor
  //
  Double_t evtnorm = 1;
  if (geteventinfo) {            //Compute the normalization from the list of events
    if (iMult != nMultBin - 1) { //MB case
      Double_t temp = hEvtMultAftEvSel->GetEffectiveEntries();
      hEvtMultAftEvSel->GetXaxis()->SetRangeUser(MultBin[iMult] + 0.0001, MultBin[iMult + 1] - 0.0001);
      hEvtMult->GetXaxis()->SetRangeUser(MultBin[iMult] + 0.0001, MultBin[iMult + 1] - 0.0001);
      if (temp == hEvtMultAftEvSel->GetEffectiveEntries())
        Fatalmsg("PrepareRawSpectra", "Normalization");
    } else { //Centrality bin
      hEvtMultAftEvSel->GetXaxis()->SetRangeUser(-1, hEvtMultAftEvSel->GetNbinsX() + 1);
      hEvtMult->GetXaxis()->SetRangeUser(-1, hEvtMult->GetNbinsX() + 1);
      if (hEvtMultAftEvSel->GetSumOfWeights() != hEvtMultAftEvSel->GetEffectiveEntries())
        Fatalmsg("PrepareRawSpectra", "Normalization");
    }
    evtnorm = 1. / hEvtMultAftEvSel->GetEffectiveEntries();
    //   evtnorm = 1./hEvtMultAftEvSel->GetSumOfWeights()
    //   evtnorm = 1./hEvtMult->GetEffectiveEntries()
  }

  //Compute the normalization from the list of analysed tracks
  const Double_t norm = 1. / fEntries->GetBinContent(1);

  if (geteventinfo && norm != evtnorm) {
    Warningmsg("PrepareRawSpectra", Form("The two Normalizations differ! (evt) %f (processed) %f evt/processed ratio %f%%\n\t The number of events in the multiplicity class is: (evt) %f (processed) %f", evtnorm, norm, norm > 0 ? evtnorm / norm : -1, hEvtMultAftEvSel->GetEffectiveEntries(), fEntries->GetBinContent(1)));
  }

  //
  //Draw the yield
  //

  cYield->cd(1);
  gPad->SetLogy();
  hYield->SetTitle("Raw Yield Unbinned");
  hYield->GetXaxis()->SetRangeUser(0, 5);
  hYield->DrawCopy()->GetYaxis()->SetTitle("Unscaled unbinned counts");
  lHistograms->AddFirst(((TH1F*)hYield->Clone(Form("hUnscaledYieldUnbinned%s%s", pC[iCharge].Data(), pS[iSpecies].Data()))));
  DrawLabel("Raw Yield Unbinned");
  cYield->cd(2);
  gPad->SetLogy();
  hYield->Scale(1, "width"); //Scale to the bin width
  hYield->SetTitle("Raw Yield");
  hYield->DrawCopy()->GetYaxis()->SetTitle("Unscaled counts");
  lHistograms->AddFirst(((TH1F*)hYield->Clone(Form("hUnscaledYield%s%s", pC[iCharge].Data(), pS[iSpecies].Data()))));
  DrawLabel("Unscaled counts");
  cYield->cd(3);
  gPad->SetLogy();
  if (normfromlist) {
    hYield->Scale(evtnorm);
    hYieldResidual->Scale(evtnorm, "width");
    hYield->SetTitle("Yield List Evt. Norm.");
  } else {
    hYield->Scale(norm);
    hYieldResidual->Scale(norm, "width");
    hYield->SetTitle("Yield Evt. Norm.");
  }
  hYield->DrawCopy()->GetYaxis()->SetTitle("Evt. Scaled counts");
  DrawLabel("Evt. Scaled counts");
  cYield->SaveAs(Form("%s%s%s_%s.pdf", outpath.Data(), pCharge[iCharge].Data(), pSpecies[iSpecies].Data(), MultBinString[iMult].Data()));

  //
  //Draw the residual of the yield
  //
  TH1* H = 0x0;
  cYieldResidual->cd(1);
  H = hYieldResidual->DrawCopy("");
  H->GetXaxis()->SetRangeUser(0, 5);
  H->SetLineColor(kRed);
  H->GetYaxis()->SetTitle("Residual");
  DrawLabel("Residual yield");
  cYieldResidual->cd(2);
  H->DrawCopy("")->Divide(hYield);
  DrawLabel("Ratio to yield");
  cYieldResidual->cd(3);
  H = hYieldResidual->DrawCopy("");
  H->SetLineColor(kRed);
  H->GetXaxis()->SetRangeUser(0, 5);
  hYield->DrawCopy("same");
  H->Add(hYield);
  DrawLabel("Yield and Yield + Residual");
  cYieldResidual->cd(4);
  H = H->DrawCopy("");
  H->Divide(hYield);
  H->GetYaxis()->SetRangeUser(.9, 1.1);
  DrawLabel("Yield + Residual Ratio to Yield");

  //
  //Draw the background
  //
  for (Int_t i = 0; i < kExpSpecies + 1; i++) {
    cBackground->cd(i + 1);
    hBackground[i]->SetMinimum(0);
    hBackground[i]->GetXaxis()->SetRangeUser(0, 5);
    hBackground[i]->DrawCopy()->Scale(1, "width");
    DrawLabel(i == 0 ? "Mismatch" : pSpecies_all[i - 1].Data());

    cBackgroundOverlap->cd(i + 1);
    hBackgroundOverlap[i]->SetMinimum(0);
    hBackgroundOverlap[i]->GetXaxis()->SetRangeUser(0, 5);
    hBackgroundOverlap[i]->DrawCopy()->Scale(1, "width");
    DrawLabel(i == 0 ? "Mismatch" : pSpecies_all[i - 1].Data());
  }

  //
  //Draw event info
  //
  cEvent->cd(1);
  if (geteventinfo)
    hNEvt->DrawCopy();
  else
    DrawLabel("Required with no event info");

  cEvent->cd(2);
  if (geteventinfo) {
    hEvtMult->GetXaxis()->SetRangeUser(0, 210);
    hEvtMult->DrawCopy();
  } else
    DrawLabel("Required with no event info");

  cEvent->cd(3);
  if (geteventinfo) {
    hEvtMultAftEvSel->GetXaxis()->SetRangeUser(0, 102);
    hEvtMultAftEvSel->DrawCopy();
  } else
    DrawLabel("Required with no event info");

  cEvent->cd(4);
  fEntries->DrawCopy();

  //Draw used templates
  if (husage) {
    TCanvas* cUsage = new TCanvas("cUsage", "Templates");
    husage->DrawCopy("COLZ");
    if (productionmode)
      delete cUsage;
  }

  if (!fitprefix.EqualTo("")) { //Save fitted histos to file for fits
    if (foutfits)
      foutfits->cd();
    for (Int_t ptbin = 0; ptbin < kPtBins; ptbin++) {
      TH1F* auxihisto = hTOF[ptbin];
      auxihisto->SetName(Form("TOFData%s%s%s_pt%i", pC[iCharge].Data(), pS[iSpecies].Data(), MultBinString[iMult].Data(), ptbin));
      auxihisto->SetTitle(Form("TOF for %s %s in %s ptbin %i [%.2f,%.2f]", pCharge[iCharge].Data(), pSpecies[iSpecies].Data(), MultBinString[iMult].Data(), ptbin, fBinPt[ptbin], fBinPt[ptbin + 1]));
      if (foutfits)
        auxihisto->Write(auxihisto->GetName());

      hFitted[ptbin]->SetName(Form("TOFFitted%s%s%s_pt%i", pC[iCharge].Data(), pS[iSpecies].Data(), MultBinString[iMult].Data(), ptbin));
      hFitted[ptbin]->SetTitle(Form("Fitted TOF for %s %s in %s ptbin %i [%.2f,%.2f]", pCharge[iCharge].Data(), pSpecies[iSpecies].Data(), MultBinString[iMult].Data(), ptbin, fBinPt[ptbin], fBinPt[ptbin + 1]));
      if (foutfits)
        hFitted[ptbin]->Write(hFitted[ptbin]->GetName());

      hMismatchFitted[ptbin]->SetName(Form("TOFFitted_Mismatch_%s%s%s_pt%i", pC[iCharge].Data(), pS[iSpecies].Data(), MultBinString[iMult].Data(), ptbin));
      hMismatchFitted[ptbin]->SetTitle(Form("Fitted TOF Mismatch for %s %s in %s ptbin %i [%.2f,%.2f]", pCharge[iCharge].Data(), pSpecies[iSpecies].Data(), MultBinString[iMult].Data(), ptbin, fBinPt[ptbin], fBinPt[ptbin + 1]));
      if (hMismatchFitted[ptbin]->GetEffectiveEntries() > 1 && foutfits)
        hMismatchFitted[ptbin]->Write(hMismatchFitted[ptbin]->GetName());

      for (Int_t species = 0; species < kExpSpecies; species++) {
        hExpectedFitted[ptbin][species]->SetName(Form("TOFFitted_%s_%s%s%s_pt%i", pS_all[species].Data(), pC[iCharge].Data(), pS[iSpecies].Data(), MultBinString[iMult].Data(), ptbin));
        hExpectedFitted[ptbin][species]->SetTitle(Form("Fitted TOF %s for %s %s in %s ptbin %i [%.2f,%.2f]", pSpecies[species].Data(), pCharge[iCharge].Data(), pSpecies[iSpecies].Data(), MultBinString[iMult].Data(), ptbin, fBinPt[ptbin], fBinPt[ptbin + 1]));
        if (hExpectedFitted[ptbin][species]->GetEffectiveEntries() <= 1)
          continue;
        if (foutfits)
          hExpectedFitted[ptbin][species]->Write(hExpectedFitted[ptbin][species]->GetName());
      }

      hRatioToFitted[ptbin]->SetName(Form("TOFFitted_Ratio_%s%s%s_pt%i", pC[iCharge].Data(), pS[iSpecies].Data(), MultBinString[iMult].Data(), ptbin));
      hRatioToFitted[ptbin]->SetTitle(Form("Fitted TOF Mismatch for %s %s in %s ptbin %i [%.2f,%.2f]", pCharge[iCharge].Data(), pSpecies[iSpecies].Data(), MultBinString[iMult].Data(), ptbin, fBinPt[ptbin], fBinPt[ptbin + 1]));
      if (foutfits)
        hRatioToFitted[ptbin]->Write(hRatioToFitted[ptbin]->GetName());
    }

    if (foutfits) {
      foutfits->Close();
      delete foutfits;
    }
  }

  if (fout) {
    fout->cd();
    lCanvas->Write(Form("lCanvasRawSpectra%s%s", pCharge[iCharge].Data(), pSpecies[iSpecies].Data()), 1);
    lHistograms->Write(Form("lHistogramsRawSpectra%s%s", pCharge[iCharge].Data(), pSpecies[iSpecies].Data()), 1);
  }

  if (productionmode) {
    if (lHistograms)
      delete lHistograms;
    if (lCanvas)
      delete lCanvas;
  }

  if (linData)
    delete linData;
  if (fout) {
    fout->Close();
    delete fout;
  }
}

void SetColor(TH1* h, Int_t c, Int_t f)
{
  h->SetMarkerColor(c);
  h->SetLineColor(c);
  if (f > 0) {
    h->SetFillStyle(f);
    h->SetFillColor(c);
  }
}

#endif
