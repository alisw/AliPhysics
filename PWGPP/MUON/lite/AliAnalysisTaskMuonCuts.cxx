/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id: AliAnalysisTaskMuonCuts.cxx 47782 2011-02-24 18:37:31Z martinez $ */

/// Task for tuning cuts for single muons in the spectrometer.
/// The output is a list of histograms.
/// Task need to be run in steps:
/// Step 1: find <DCAx> and <DCAy> in MC
/// Step 2: find signal sigma_pDCA in MC
/// Step 3: find <DCAx> and <DCAy> in data
/// Step 4: apply sigma_pDCA cut on simulations and check the
///         resulting kinematic variables
///
/// \author Diego Stocco

#define AliAnalysisTaskMuonCuts_cxx

#include "AliAnalysisTaskMuonCuts.h"

// ROOT includes
#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMath.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TF1.h"
#include "TStyle.h"
//#include "TMCProcess.h"

// STEER includes
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAODMCParticle.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
//#include "AliStack.h"
#include "AliESDEvent.h"
#include "AliESDMuonTrack.h"
#include "AliInputEventHandler.h"

#include "AliAnalysisManager.h"

// PWG includes
#include "AliMergeableCollection.h"
#include "AliMuonTrackCuts.h"


/// \cond CLASSIMP
ClassImp(AliAnalysisTaskMuonCuts) // Class implementation in ROOT context
/// \endcond


//________________________________________________________________________
AliAnalysisTaskMuonCuts::AliAnalysisTaskMuonCuts() :
  AliVAnalysisMuon(),
  fHistoTypeKeys(0x0),
  fThetaAbsKeys(0x0),
  fSigmaCuts(TArrayD())
{
  /// Default ctor.
}

//________________________________________________________________________
AliAnalysisTaskMuonCuts::AliAnalysisTaskMuonCuts(const char *name, const AliMuonTrackCuts& cuts) :
  AliVAnalysisMuon(name, cuts),
  fHistoTypeKeys(0x0),
  fThetaAbsKeys(0x0),
  fSigmaCuts(TArrayD())
{
  //
  /// Constructor.
  //
  TString histoTypeKeys = "hDCAxVsP hDCAyVsP hPdcaVsP hPDCAVsPCheck hDCAVsPCheck hChiProbVsP hSigmaVsPt hSigmaVsEta";
  fHistoTypeKeys = histoTypeKeys.Tokenize(" ");

  TString thetaAbsKeys = "ThetaAbs23 ThetaAbs310";
  fThetaAbsKeys = thetaAbsKeys.Tokenize(" ");

  SetSigmaCuts();
}


//________________________________________________________________________
AliAnalysisTaskMuonCuts::~AliAnalysisTaskMuonCuts()
{
  //
  /// Destructor
  //

  delete fHistoTypeKeys;
  delete fThetaAbsKeys;
}

//___________________________________________________________________________
void AliAnalysisTaskMuonCuts::MyUserCreateOutputObjects()
{
  TH1* histo = 0x0;
  TString histoName = "";
  for ( Int_t itheta=0; itheta<kNthetaAbs; ++itheta ) {
    for ( Int_t isrc=0; isrc<kNtrackSources; ++isrc ) {
      histoName = GetHistoName(kDCAxVsP, itheta, isrc);
      histo = new TH2F(histoName.Data(), histoName.Data(), 50, 0., 200., 100, -50., 50.);
      histo->SetXTitle("p (GeV/c)");
      histo->SetYTitle("DCA_{x} (cm)");
      AddObjectToCollection(histo);

      histoName = GetHistoName(kDCAyVsP, itheta, isrc);
      histo = new TH2F(histoName.Data(), histoName.Data(), 50, 0., 200., 100, -50., 50.);
      histo->SetXTitle("p (GeV/c)");
      histo->SetYTitle("DCA_{y} (cm)");
      AddObjectToCollection(histo);

      histoName = GetHistoName(kPdcaVsP, itheta, isrc);
      histo = new TH2F(histoName.Data(), histoName.Data(), 100, 0., 400., 50, 0., 500.);
      histo->SetXTitle("p (GeV/c)");
      histo->SetYTitle("p#timesDCA (cm #times GeV/c)");
      AddObjectToCollection(histo);

      histoName = GetHistoName(kPDCAVsPCheck, itheta, isrc);
      histo = new TH2F(histoName.Data(), histoName.Data(), 100, 0., 800., 100, 0., 5000.);
      histo->SetXTitle("p (GeV/c)");
      histo->SetYTitle("p#timesDCA (cm #times GeV/c)");
      AddObjectToCollection(histo);

      histoName = GetHistoName(kDCAVsPCheck, itheta, isrc);
      histo = new TH2F(histoName.Data(), histoName.Data(), 100, 0., 800., 100, 0., 200.);
      histo->SetXTitle("p (GeV/c)");
      histo->SetYTitle("DCA (cm)");
      AddObjectToCollection(histo);

      histoName = GetHistoName(kChiProbVsP, itheta, isrc);
      histo = new TH2F(histoName.Data(), histoName.Data(), 50, 0., 200., 50, 0., 1.);
      histo->SetXTitle("p (GeV/c)");
      histo->SetYTitle("Chisquare prob.");
      AddObjectToCollection(histo);

      histoName = GetHistoName(kSigmaVsPt, itheta, isrc);
      histo = new TH2F(histoName.Data(), histoName.Data(), 100, 0., 100., fSigmaCuts.GetSize(), 0.5, 0.5+(Double_t)fSigmaCuts.GetSize());
      histo->SetXTitle("p_{t} (GeV/c)");
      histo->SetYTitle("#sigma_{p#timesDCA}");
      for ( Int_t ibin=0; ibin<fSigmaCuts.GetSize(); ++ibin ) {
        histo->GetYaxis()->SetBinLabel(ibin+1,Form("%g", fSigmaCuts[ibin]));
      }
      AddObjectToCollection(histo);

      histoName = GetHistoName(kSigmaVsEta, itheta, isrc);
      histo = new TH2F(histoName.Data(), histoName.Data(), 25, -4.5, -2., fSigmaCuts.GetSize(), 0.5, 0.5+(Double_t)fSigmaCuts.GetSize());
      histo->SetXTitle("#eta");
      histo->SetYTitle("#sigma_{p#timesDCA}");
      for ( Int_t ibin=0; ibin<fSigmaCuts.GetSize(); ++ibin ) {
        histo->GetYaxis()->SetBinLabel(ibin+1,Form("%g", fSigmaCuts[ibin]));
      }
      AddObjectToCollection(histo);
    } // loop on track sources
  } // loop on theta abs

  fMuonTrackCuts->Print();

}

//________________________________________________________________________
TString AliAnalysisTaskMuonCuts::GetHistoName(Int_t histoTypeIndex, Int_t thetaAbsIndex, Int_t srcIndex)
{
  /// Get local histogram name
  TString histoName = Form("%s%s%s", fHistoTypeKeys->At(histoTypeIndex)->GetName(), fThetaAbsKeys->At(thetaAbsIndex)->GetName(), fSrcKeys->At(srcIndex)->GetName());

  return histoName;
}

//________________________________________________________________________
void AliAnalysisTaskMuonCuts::ProcessEvent(TString physSel, const TObjArray& selectTrigClasses, TString centrality)
{
  //
  /// Fill histogram
  //

  if ( GetVertexSPD()->GetNContributors() < fMinNvtxContirbutors ) return;

  TString histoName = "";
  AliVParticle* track = 0x0;
  Int_t nTracks = GetNTracks();
  for (Int_t itrack = 0; itrack < nTracks; itrack++) {
    track = GetTrack(itrack);
    fMuonTrackCuts->SetNSigmaPdca(1.e10);
    if ( ! fMuonTrackCuts->IsSelected(track) ) continue;

    Double_t rAbsEnd =  ( fAODEvent ) ? ((AliAODTrack*)track)->GetRAtAbsorberEnd(): ((AliESDMuonTrack*)track)->GetRAtAbsorberEnd();
    Double_t thetaAbsEndDeg = TMath::ATan( rAbsEnd / 505. ) * TMath::RadToDeg();
    Int_t thetaAbsBin = ( thetaAbsEndDeg < 3. ) ? kThetaAbs23 : kThetaAbs310;

    Int_t trackSrc = GetParticleType(track);

    TVector3 dcaAtVz = fMuonTrackCuts->GetCorrectedDCA(track);
    Double_t pTotMean = fMuonTrackCuts->GetAverageMomentum(track);
    Double_t dca = dcaAtVz.Mag();
    Double_t pDca = pTotMean * dca;

    Double_t chi2 = pDca / fMuonTrackCuts->GetSigmaPdca(rAbsEnd) ;
    chi2 *= chi2;
    Double_t chiProb = TMath::Prob(chi2, 2);

    Double_t pTot = track->P();
    Double_t pt = track->Pt();
    Double_t eta = track->Eta();

    for ( Int_t itrig=0; itrig<selectTrigClasses.GetEntries(); ++itrig ) {
      TString trigClassName = ((TObjString*)selectTrigClasses.At(itrig))->GetString();
      histoName = GetHistoName(kDCAxVsP, thetaAbsBin, trackSrc);
      ((TH2*)GetMergeableObject(physSel, trigClassName, centrality, histoName))->Fill(pTot, dcaAtVz.X());

      histoName = GetHistoName(kDCAyVsP, thetaAbsBin, trackSrc);
      ((TH2*)GetMergeableObject(physSel, trigClassName, centrality, histoName))->Fill(pTot, dcaAtVz.Y());
  
      histoName = GetHistoName(kPdcaVsP, thetaAbsBin, trackSrc);
      ((TH2*)GetMergeableObject(physSel, trigClassName, centrality, histoName))->Fill(pTot, pDca);

      histoName = GetHistoName(kPDCAVsPCheck, thetaAbsBin, trackSrc);
      ((TH2*)GetMergeableObject(physSel, trigClassName, centrality, histoName))->Fill(pTot, pDca);

      histoName = GetHistoName(kDCAVsPCheck, thetaAbsBin, trackSrc);
      ((TH2*)GetMergeableObject(physSel, trigClassName, centrality, histoName))->Fill(pTot, dca);

      histoName = GetHistoName(kChiProbVsP, thetaAbsBin, trackSrc);
      ((TH2*)GetMergeableObject(physSel, trigClassName, centrality, histoName))->Fill(pTot, chiProb);
    } // loop on selected trigger classes

    for ( Int_t isigma=0; isigma<fSigmaCuts.GetSize(); ++isigma) {
      fMuonTrackCuts->SetNSigmaPdca(fSigmaCuts[isigma]);
      if ( ! fMuonTrackCuts->IsSelected(track) ) continue;
      for ( Int_t itrig=0; itrig<selectTrigClasses.GetEntries(); ++itrig ) {
        TString trigClassName = ((TObjString*)selectTrigClasses.At(itrig))->GetString();
        histoName = GetHistoName(kSigmaVsPt, thetaAbsBin, trackSrc);
        ((TH2*)GetMergeableObject(physSel, trigClassName, centrality, histoName))->Fill(pt, isigma+1);
        histoName = GetHistoName(kSigmaVsEta, thetaAbsBin, trackSrc);
        ((TH2*)GetMergeableObject(physSel, trigClassName, centrality, histoName))->Fill(eta, isigma+1);
      } // loop on selected trigger classes 
    } // loop on sigmas
  }
}

//________________________________________________________________________
void AliAnalysisTaskMuonCuts::SetSigmaCuts(Int_t nSigmaCuts, Double_t* sigmaCuts)
{
  /// Set number of sigmas
  if ( ! sigmaCuts || nSigmaCuts < 0 ) {
//     if ( defaultChiSquare ) {
//       Double_t cuts[] = {0.3, 0.1, 0.05, 0.03, 0.01, 0.005, 0.003, 0.001, 0.0005, 0.0003, 0.0001, 0.};
//       Int_t nCuts = sizeof(cuts)/sizeof(cuts[0]);
//       fSigmaCuts.Set(nCuts, cuts);
//     }
//     else {
      Double_t cuts[] = {2., 3., 4., 5., 6., 7., 10., 15., 20., 25., 30., 1.e10};
      Int_t nCuts = sizeof(cuts)/sizeof(cuts[0]);
      fSigmaCuts.Set(nCuts, cuts);
      //    }
  }
  else {
    fSigmaCuts.Set(nSigmaCuts, sigmaCuts);
  }
}

//________________________________________________________________________
void AliAnalysisTaskMuonCuts::Terminate(Option_t *) {
  //
  /// Draw some histogram at the end.
  //

  AliVAnalysisMuon::Terminate("");

  if ( ! fMergeableCollection ) return;

  TString physSel = fTerminateOptions->At(0)->GetName();
  TString trigClassName = fTerminateOptions->At(1)->GetName();
  TString centralityRange = fTerminateOptions->At(2)->GetName();
  TString furtherOpt = fTerminateOptions->At(3)->GetName();
  furtherOpt.ToUpper();

  furtherOpt.ReplaceAll("  ", " ");
  furtherOpt.ReplaceAll(" =", "=");
  furtherOpt.ReplaceAll("= ", "=");
  TObjArray* optArray = furtherOpt.Tokenize(" ");
  Double_t refSigmaCut = 6.;
  for ( Int_t iopt=0; iopt<optArray->GetEntries(); ++iopt ) {
    TString currOpt = optArray->At(iopt)->GetName();
    if ( currOpt.Contains("REFSIGMA") ) {
      currOpt.Remove(0,currOpt.Index("=")+1);
      refSigmaCut = currOpt.Atof();
    }
  }
  delete optArray;

  Int_t srcColors[kNtrackSources] = {kBlack, kRed, kGreen, kBlue, kViolet, 7, kOrange};

  TCanvas* can = 0x0;
  Int_t xshift = 100;
  Int_t yshift = 100;
  Int_t igroup1 = -1;
  Int_t igroup2 = 0;

  //////////////
  // Reco DCA //
  //////////////
  igroup1++;
  igroup2 = 0;
  TString histoName = "", currName = "", histoPattern = "", drawOpt = "";
  currName = Form("%s_recoDCA", GetName());
  can = new TCanvas(currName.Data(),"Reco DCA",igroup1*xshift,igroup2*yshift,600,600);
  can->Divide(2,2);
  igroup2++;
  Int_t recoDcaHisto[2] = {kDCAxVsP, kDCAyVsP};
  TString dcaName[2] = {"DCAx", "DCAy"};
  Double_t averageDca[4] = {0., 0.};
  printf("\nAverage reconstructed DCA:\n");
  TF1* fitFuncMeanDca = new TF1("fitFuncMeanDca","gausn",-20.,20.);
  fitFuncMeanDca->SetParNames("Norm", "Mean", "Sigma");
  for ( Int_t itheta=0; itheta<kNthetaAbs; ++itheta ) {
    for ( Int_t ihisto=0; ihisto<2; ++ihisto ) {
      TH2* histo = 0x0;
      histoPattern = "";
      histoPattern = Form("%s & %s", fHistoTypeKeys->At(recoDcaHisto[ihisto])->GetName(), fThetaAbsKeys->At(itheta)->GetName());
      histo = (TH2*)GetSum(physSel, trigClassName, centralityRange, histoPattern);
      if ( ! histo ) continue;

      TH1* meanDcaVsP = histo->ProjectionX(Form("mean%sVsP_%s", dcaName[ihisto].Data(), fThetaAbsKeys->At(itheta)->GetName()));
      meanDcaVsP->Reset();
      meanDcaVsP->SetTitle(Form("Mean %s vs. p %s", dcaName[ihisto].Data(), fThetaAbsKeys->At(itheta)->GetName()));
      meanDcaVsP->SetYTitle(Form("<%s> (cm)", dcaName[ihisto].Data()));
      meanDcaVsP->SetStats(kFALSE);

      Int_t nPbins = histo->GetXaxis()->GetNbins();
      //Int_t nPadX = (Int_t)TMath::Sqrt(nPbins);
      //Int_t nPadY = nPadX;
      //if ( nPadX * nPadY < nPbins ) nPadX++;
      TCanvas* meanDcaFitCan = 0x0;
      Int_t nPadX = 5;
      Int_t nPadY = 5;
      Int_t ipad = 0;
      Int_t ican = 0;

      for ( Int_t ibin=2; ibin<=nPbins; ++ibin ) {
        currName = Form("hMean%s_%s_%s", dcaName[ihisto].Data(), physSel.Data(), trigClassName.Data());
        Int_t minBin = ( ibin == 0 ) ? 1 : ibin;
        Int_t maxBin = ( ibin == 0 ) ? nPbins : ibin;
        if ( ibin > 0 ) currName += Form("_pBin%i", ibin);
        TH1* projHisto = histo->ProjectionY(currName.Data(), minBin, maxBin, "e");
        projHisto->SetTitle(Form("%s %s %g < p < %g (GeV/c)", dcaName[ihisto].Data(), fThetaAbsKeys->At(itheta)->GetName(), meanDcaVsP->GetXaxis()->GetBinLowEdge(minBin), meanDcaVsP->GetXaxis()->GetBinUpEdge(maxBin)));

        if ( projHisto->GetEntries() == 0 ) continue;
        if ( ipad % (nPadX*nPadY) == 0 ) {
          currName = histo->GetName();
          currName += Form("Fit_can_%i", ican++);
          meanDcaFitCan = new TCanvas(currName.Data(), currName.Data(),igroup1*xshift,igroup2*yshift,600,600);
          meanDcaFitCan->Divide(nPadX,nPadY);
          ipad = 0;
        }
        meanDcaFitCan->cd(++ipad);
        gPad->SetLogy();
        if ( projHisto->Integral() > 50 ) {
          fitFuncMeanDca->SetParameter(0, projHisto->Integral());
          fitFuncMeanDca->SetParameter(1, projHisto->GetMean());
          fitFuncMeanDca->SetParameter(2, projHisto->GetRMS());
          Double_t fitDcaLim = ( ibin == 0 ) ?  40. :  40./((Double_t)ibin);
          fitDcaLim = TMath::Max(5., fitDcaLim);
          projHisto->Fit(fitFuncMeanDca, "RQ", "e", -fitDcaLim, fitDcaLim);
          Double_t chi2 = fitFuncMeanDca->GetChisquare();
          Double_t ndf = fitFuncMeanDca->GetNDF();
          if ( ndf <= 0.) continue;
          if ( chi2 / ndf > 100. ) continue;
          Double_t meanDca = fitFuncMeanDca->GetParameter(1);
          Double_t meanDcaErr = fitFuncMeanDca->GetParError(1);
          meanDcaVsP->SetBinContent(ibin, meanDca);
          meanDcaVsP->SetBinError(ibin, meanDcaErr);
          //if ( ibin == 0 ) printf("%s  meanDca = %g +- %g\n", fThetaAbsKeys->At(itheta)->GetName(),  meanDca, meanDcaErr);
        }
        else projHisto->Draw("e");
      } // loop on momentum bins
      can->cd(2*itheta+ihisto+1);
      //meanDcaVsP->SetLineColor(srcColors[isrc]);
      //meanDcaVsP->SetMarkerColor(srcColors[isrc]);
      //meanDcaVsP->SetMarkerStyle(20+isrc);
      meanDcaVsP->Fit("pol0","Q");
      TF1* trendFit = (TF1*)meanDcaVsP->GetListOfFunctions()->FindObject("pol0");
      if ( trendFit ) {
        averageDca[2*itheta+ihisto] = trendFit->GetParameter(0);
        printf("  %s: mean %s = %g +- %g  (chi2/%i = %g)\n",  fThetaAbsKeys->At(itheta)->GetName(), dcaName[ihisto].Data(), trendFit->GetParameter(0), trendFit->GetParError(0), trendFit->GetNDF(), trendFit->GetChisquare()/((Double_t)trendFit->GetNDF()));
      }
        //drawOpt = "esame";
      //leg->AddEntry(meanDcaVsP, fSrcKeys->At(isrc)->GetName(), "lp");
      //} // loop on src
    //can->cd(itheta+1);
    //leg->Draw("same");

      //can->cd(ipad++);
      //histo->Draw();
      //meanDca[ihisto] = histo->GetMean();
    } // loop on histo type
  } // loop on theta abs

  for ( Int_t itheta=0; itheta<kNthetaAbs; ++itheta ) {  
    printf("muonCuts->SetMeanDCA(%g, %g, 0.); // %s\n", averageDca[2*itheta], averageDca[2*itheta+1], fThetaAbsKeys->At(itheta)->GetName());
  }

  //////////////
  // Fit pDCA //
  //////////////
  Double_t nSigmaCut = fMuonTrackCuts->GetNSigmaPdca(); //6.;
  Double_t sigmaMeasCut[2] = { fMuonTrackCuts->GetSigmaPdca(AliMuonTrackCuts::kThetaAbs23), fMuonTrackCuts->GetSigmaPdca(AliMuonTrackCuts::kThetaAbs23)}; //{99., 54.}; //{120., 63.};
  Double_t relPResolution = fMuonTrackCuts->GetRelPResolution(); //4.5e-4; //6.e-3;//8.e-4;
  Double_t angleResolution = 535.*fMuonTrackCuts->GetSlopeResolution(); //6.e-4;
  Double_t pMinCut = 0.1;
  Double_t pMaxCut =  800.;
  const Int_t kNCutFuncs = 2;
  Int_t nShowFuncs = 1;
  TString cutFormula = "[1]*TMath::Sqrt( ( [0] / ( 1. - [1]*[2]*x / ( 1.+[1]*[2]*x ) ) ) * ( [0] / ( 1. - [1]*[2]*x / ( 1.+[1]*[2]*x ) ) ) + [3]*[3]*x*x)";
  Double_t cutParam[kNCutFuncs][4] = {{sigmaMeasCut[0], nSigmaCut, relPResolution, angleResolution}, {sigmaMeasCut[0], nSigmaCut, 0., 0.32}};
  Int_t cutColor[kNCutFuncs] = {kBlack, kRed};
  igroup1++;
  igroup2 = 0;
  can = new TCanvas("pdcaSigmaFit","Sigma vs P fit",igroup1*xshift,igroup2*yshift,600,600);
  can->Divide(2,1);
  igroup2++;
  TF1* fitFunc = new TF1("fitFunc", "x*gausn", 0., 400.);
  fitFunc->SetParNames("Norm", "Mean", "Sigma");
  gStyle->SetOptFit(1111);
  Double_t xMinFit[2] = {0., 0.};
  Double_t xMaxFit[2] = {320., 150.}; // {360., 180.};
  printf("\nSigma p x DCA:\n");
  Double_t averageSigmaPdca[kNtrackSources*kNthetaAbs] = {0.};
  for ( Int_t itheta=0; itheta<kNthetaAbs; ++itheta ) {
    TLegend* leg = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg->SetBorderSize(1);
    for ( Int_t isrc=0; isrc<kNtrackSources; ++isrc ) {
      histoPattern = GetHistoName(kPdcaVsP, itheta, isrc);
      TH2* histo = (TH2*)GetSum(physSel, trigClassName, centralityRange, histoPattern);
      if ( ! histo ) continue;
      if ( histo->Integral() < 200 ) {
        delete histo;
        continue;
      }

      TH1* sigmaVsP = histo->ProjectionX(Form("sigma%s_%s_%s", fHistoTypeKeys->At(kPdcaVsP)->GetName(), fThetaAbsKeys->At(itheta)->GetName(), fSrcKeys->At(isrc)->GetName()));
      sigmaVsP->Reset();
      sigmaVsP->SetTitle(Form("#sigma_{p#timesDCA} vs. p %s %s", fThetaAbsKeys->At(itheta)->GetName(), fSrcKeys->At(isrc)->GetName()));
      sigmaVsP->SetYTitle("#sigma_{p#timesDCA} (cm #times GeV/c)");
      sigmaVsP->SetStats(kFALSE);

      Int_t nPbins = histo->GetXaxis()->GetNbins();
      //Int_t nPadX = (Int_t)TMath::Sqrt(nPbins);
      //Int_t nPadY = nPadX;
      //if ( nPadX * nPadY < nPbins ) nPadX++;
      TCanvas* pdcaFitCan = 0x0;
      Int_t nPadX = 5;
      Int_t nPadY = 5;
      Int_t ipad = 0;
      Int_t ican = 0;

      for ( Int_t ibin=2; ibin<=nPbins; ++ibin ) {
        currName = Form("hPDCA_%s_%s_%s", fSrcKeys->At(isrc)->GetName(), physSel.Data(), trigClassName.Data());
        Int_t minBin = ( ibin == 0 ) ? 1 : ibin;
        Int_t maxBin = ( ibin == 0 ) ? nPbins : ibin;
        if ( ibin > 0 ) currName += Form("_pBin%i", ibin);
        TH1* projHisto = histo->ProjectionY(currName.Data(), minBin, maxBin, "e");
        if ( projHisto->GetEntries() == 0 ) continue;
        projHisto->SetTitle(Form("P DCA %s %s %g < p < %g (GeV/c)", fThetaAbsKeys->At(itheta)->GetName(), fSrcKeys->At(isrc)->GetName(), sigmaVsP->GetXaxis()->GetBinLowEdge(minBin), sigmaVsP->GetXaxis()->GetBinUpEdge(maxBin)));
        if ( ipad % (nPadX*nPadY) == 0 ) {
          currName = histo->GetName();
          currName += Form("Fit_can_%i", ican++);
          pdcaFitCan = new TCanvas(currName.Data(), currName.Data(),igroup1*xshift,igroup2*yshift,600,600);
          pdcaFitCan->Divide(nPadX,nPadY);
          ipad = 0;
        }
        pdcaFitCan->cd(++ipad);
        if ( projHisto->Integral() > 0.) projHisto->Scale(1./projHisto->Integral());
        gPad->SetLogy();
        if ( projHisto->GetEntries() > 50 ) {
          fitFunc->SetParameter(0, projHisto->Integral());
          fitFunc->SetParameter(1, projHisto->GetMean());
          fitFunc->SetParameter(2, projHisto->GetRMS());
          projHisto->Fit(fitFunc, "RQ", "e", xMinFit[itheta], xMaxFit[itheta]);
          Double_t chi2 = fitFunc->GetChisquare();
          Double_t ndf = fitFunc->GetNDF();
          if ( ndf <= 0.) continue;
          if ( chi2 / ndf > 100. ) continue;
          Double_t sigma = TMath::Abs(fitFunc->GetParameter(2));
          Double_t sigmaErr = fitFunc->GetParError(2);
          sigmaVsP->SetBinContent(ibin, sigma);
          sigmaVsP->SetBinError(ibin, sigmaErr);
          //if ( ibin == 0 ) printf("%s %s  sigma = %g +- %g\n", fThetaAbsKeys->At(itheta)->GetName(), fSrcKeys->At(isrc)->GetName(), sigma, sigmaErr);
        }
        else projHisto->Draw("e");
      } // loop on momentum bins
      can->cd(itheta+1);
      sigmaVsP->SetLineColor(srcColors[isrc]);
      sigmaVsP->SetMarkerColor(srcColors[isrc]);
      sigmaVsP->SetMarkerStyle(20+isrc);
      drawOpt = "e";
      if ( gPad->GetListOfPrimitives()->GetEntries() != 0 ) drawOpt += "same";
      sigmaVsP->Draw(drawOpt.Data());
      sigmaVsP->Fit("pol0","Q",drawOpt.Data());
      TF1* trendFit = (TF1*)sigmaVsP->GetListOfFunctions()->FindObject("pol0");
      if ( trendFit ) {
        averageSigmaPdca[kNtrackSources*itheta+isrc] = trendFit->GetParameter(0);
        printf("  %s %s: Sigma = %g +- %g  (chi2/%i = %g)\n",  fThetaAbsKeys->At(itheta)->GetName(), fSrcKeys->At(isrc)->GetName(), trendFit->GetParameter(0), trendFit->GetParError(0), trendFit->GetNDF(), trendFit->GetChisquare()/((Double_t)trendFit->GetNDF()));
      }
      leg->AddEntry(sigmaVsP, fSrcKeys->At(isrc)->GetName(), "lp");

      // Plot 2D function for check!
      histoPattern = GetHistoName(kPDCAVsPCheck, itheta, isrc);
      TH2* histoCheck = (TH2*)GetSum(physSel, trigClassName, centralityRange, histoPattern);
      if ( ! histoCheck ) histoCheck = histo; // needed for old data
      currName = histoCheck->GetName();
      currName.Append("_plotCut");
      TCanvas* pdcaVsPcan = new TCanvas(currName.Data(), currName.Data(), igroup1*xshift,(igroup2+1)*yshift,600,600);
      pdcaVsPcan->SetLogz();
      pdcaVsPcan->SetRightMargin(0.12);
      histoCheck->Draw("COLZ");

      for ( Int_t icut=0; icut<nShowFuncs; ++icut ) {
        currName = Form("%s_cutFunc%i", histo->GetName(), icut);
        TF1* cutFunction = new TF1(currName.Data(),cutFormula.Data(), pMinCut, pMaxCut);
        cutParam[icut][0] = sigmaMeasCut[itheta];
        cutParam[icut][1] = nSigmaCut;
        cutFunction->SetParameters(cutParam[icut]);
        cutFunction->SetLineWidth(2);
        cutFunction->SetLineColor(cutColor[icut]);
        cutFunction->Draw("same");
      } // loop on cut func
    } // loop on src
    can->cd(itheta+1);
    leg->Draw("same");

    for ( Int_t icut=0; icut<nShowFuncs; ++icut ) {
      currName = Form("sigmaCut_%s_%i", fThetaAbsKeys->At(itheta)->GetName(), icut);
      TF1* cutFunction = new TF1(currName.Data(), cutFormula.Data(), pMinCut, pMaxCut);
      cutParam[icut][0] = sigmaMeasCut[itheta];
      cutParam[icut][1] = 1.;
      cutFunction->SetParameters(cutParam[icut]);
      cutFunction->SetLineColor(cutColor[icut]);
      cutFunction->SetLineWidth(2);
      cutFunction->Draw("same");
    }
  } // loop on theta abs

  for ( Int_t isrc=0; isrc<kNtrackSources; ++isrc ) {
    printf("muonCuts->SetSigmaPdca(%g, %g); // %s\n", averageSigmaPdca[isrc], averageSigmaPdca[kNtrackSources+isrc], fSrcKeys->At(isrc)->GetName());
  }
  printf("\n");

  igroup2++;
  for ( Int_t itheta=0; itheta<kNthetaAbs; ++itheta ) {
    for ( Int_t isrc=0; isrc<kNtrackSources; ++isrc ) {
      histoPattern = GetHistoName(kDCAVsPCheck, itheta, isrc);
      TH2* histoCheck = (TH2*)GetSum(physSel, trigClassName, centralityRange, histoPattern);
      if ( ! histoCheck ) continue;
      currName = histoCheck->GetName();
      currName.Append("_plotCut");
      TCanvas* pdcaVsPcan = new TCanvas(currName.Data(), currName.Data(), igroup1*xshift,(igroup2+1)*yshift,600,600);
      pdcaVsPcan->SetRightMargin(0.12);
      pdcaVsPcan->SetLogz();
      histoCheck->Draw("COLZ");

      for ( Int_t icut=0; icut<nShowFuncs; ++icut ) {
        currName = histoCheck->GetName();
        currName.Append(Form("_cutFunc%i",icut));
        TString currFormula = cutFormula;
        currFormula.Append("/x");
        TF1* cutFunction = new TF1(currName.Data(),currFormula.Data(), pMinCut, pMaxCut);
        cutParam[icut][0] = sigmaMeasCut[itheta];
        cutParam[icut][1] = nSigmaCut;
        cutFunction->SetParameters(cutParam[icut]);
        cutFunction->SetLineWidth(2);
        cutFunction->SetLineColor(cutColor[icut]);
        cutFunction->Draw("same");
      } // loop on cut functions
    } // loop on src
  } //loop on theta


  ///////////////////////////
  // Check Chi square prob //
  ///////////////////////////
  igroup1++;
  igroup2 = 0;
  for ( Int_t itheta=0; itheta<kNthetaAbs; ++itheta ) {
    for ( Int_t isrc=0; isrc<kNtrackSources; ++isrc ) {
      histoPattern = GetHistoName(kChiProbVsP, itheta, isrc);
      TH2* histo = (TH2*)GetSum(physSel, trigClassName, centralityRange, histoPattern);
      if ( ! histo ) continue;

      Int_t nPbins = histo->GetXaxis()->GetNbins();
      //Int_t nPadX = (Int_t)TMath::Sqrt(nPbins);
      //Int_t nPadY = nPadX;
      //if ( nPadX * nPadY < nPbins ) nPadX++;
      TCanvas* chi2probCan = 0x0;
      Int_t nPadX = 5;
      Int_t nPadY = 5;
      Int_t ipad = 0;
      Int_t ican = 0;

      for ( Int_t ibin=0; ibin<=nPbins; ++ibin ) {
        currName = Form("hChiProb_%s_%s_%s", fSrcKeys->At(isrc)->GetName(), physSel.Data(), trigClassName.Data());
        Int_t minBin = ( ibin == 0 ) ? 1 : ibin;
        Int_t maxBin = ( ibin == 0 ) ? nPbins : ibin;
        if ( ibin > 0 ) currName += Form("_pBin%i", ibin);
        TH1* projHisto = histo->ProjectionY(currName.Data(), minBin, maxBin);
        projHisto->SetTitle(Form("Chisquare prob %s %s %g < p < %g (GeV/c)", fThetaAbsKeys->At(itheta)->GetName(), fSrcKeys->At(isrc)->GetName(), histo->GetXaxis()->GetBinLowEdge(minBin), histo->GetXaxis()->GetBinUpEdge(maxBin)));
        if ( projHisto->GetEntries() == 0 ) continue;
        if ( ipad % (nPadX*nPadY) == 0 ) {
          currName = histo->GetName();
          currName += Form("Fit_can_%i", ican++);
          chi2probCan = new TCanvas(currName.Data(), currName.Data(),igroup1*xshift,igroup2*yshift,600,600);
          chi2probCan->Divide(nPadX,nPadY);
          ipad = 0;
        }
        chi2probCan->cd(++ipad);
        gPad->SetLogy();
        projHisto->Draw("e");
      } // loop on momentum bins
    } // loop on src
  } // loop on theta abs


  //////////////////////
  // Check sigma cuts //
  //////////////////////
  printf("\nReference sigma cut %g\n", refSigmaCut);

  Float_t fracOfHeight = 0.35;
  Float_t rightMargin = 0.03;
  Int_t cutColors[14] = {kBlack, kRed, kBlue, kGreen, kCyan, kMagenta, kOrange, kViolet, kSpring, kGray, kSpring, kAzure, kPink, kYellow};
  Int_t* orderCuts = 0x0;
  Int_t nSigmaCuts = 0;
  Int_t checkHistos[2] = {kSigmaVsPt, kSigmaVsEta};
  Bool_t useCustomSigma = furtherOpt.Contains("CUSTOMSIGMA");

  igroup1++;
  igroup2 = 0;

  for ( Int_t ihisto=0; ihisto<2; ++ihisto ) {
    for ( Int_t isrc=0; isrc<kNtrackSources; ++isrc ) {
      can = 0x0;
      for ( Int_t itheta=0; itheta<kNthetaAbs; ++itheta ) {
        histoPattern = GetHistoName(checkHistos[ihisto], itheta, isrc);
        TH2* histo = (TH2*)GetSum(physSel, trigClassName, centralityRange, histoPattern);
        if ( ! histo ) continue;
        if ( histo->Integral() == 0. ) {
          delete histo;
          continue;
        }
        if ( ! orderCuts ) {
          // Re-order axis
          TAxis* axis = histo->GetYaxis();
          nSigmaCuts = ( useCustomSigma ) ? fSigmaCuts.GetSize() : axis->GetNbins();
          orderCuts = new Int_t[nSigmaCuts];
          Int_t countBin = 0;
          for ( Int_t isigma=0; isigma<axis->GetNbins(); ++isigma ) {
            TString currLabel = axis->GetBinLabel(isigma+1);
            Double_t currSigma = currLabel.Atof();
            if ( useCustomSigma ) {
              Bool_t foundMatch = kFALSE;
              for ( Int_t jsigma=0; jsigma<fSigmaCuts.GetSize(); ++jsigma ) {
                if ( TMath::Abs(currSigma - fSigmaCuts[jsigma]) < 1e-4 ) {
                  foundMatch = kTRUE;
                  break;
                }
              }
              if ( ! foundMatch ) continue; 
            }
            Int_t currBin = ( TMath::Abs(currSigma - refSigmaCut) < 1e-4) ? 0 : ++countBin;
            orderCuts[currBin] = isigma+1;
          }
        }
        currName = histo->GetName();
        currName.Append("_can");
        can = new TCanvas(currName.Data(), currName.Data(), igroup1*xshift,igroup2*yshift, 600, 600);
        TLegend* leg = new TLegend(0.6,0.6,0.9,0.9);
        leg->SetBorderSize(1);
        can->Divide(1,2,0,0);
        can->cd(1);
        gPad->SetPad(0., fracOfHeight, 0.99, 0.99);
        //gPad->SetTopMargin(0.08);
        gPad->SetTopMargin(0.03);
        gPad->SetRightMargin(rightMargin);
        gPad->SetLogy();

        can->cd(2);
        gPad->SetPad(0., 0., 0.99, fracOfHeight);
        gPad->SetRightMargin(rightMargin);
        gPad->SetBottomMargin(0.08/fracOfHeight);

        TH1* refCutHisto = 0x0;
        TString legTitle = "";
        for ( Int_t isigma=0; isigma<nSigmaCuts; ++isigma ) {
          currName = histo->GetName();
          currName.Append(Form("_sigma%i", isigma));
          Int_t currBin = orderCuts[isigma];
          TH1* projectHisto = histo->ProjectionX(currName.Data(), currBin, currBin);
          projectHisto->SetLineColor(cutColors[isigma]);
          projectHisto->SetMarkerColor(cutColors[isigma]);
          projectHisto->SetMarkerStyle(20+isigma);
          projectHisto->SetStats(kFALSE);
          projectHisto->SetTitle("");
          can->cd(1);
          drawOpt = "e";
          if ( gPad->GetListOfPrimitives()->GetEntries() != 0 ) drawOpt += "same";
          projectHisto->Draw(drawOpt.Data());
          TString currLabel = histo->GetYaxis()->GetBinLabel(currBin);
          TString cutName = ( TMath::Abs(currLabel.Atof() ) > 1000. ) ? "Total" :  Form("%s #sigma cut", currLabel.Data());
          leg->AddEntry(projectHisto, cutName.Data(), "lp");
          if ( ! refCutHisto ) {
            refCutHisto = projectHisto;
            legTitle = Form("(pass n #sigma) / (pass %s #sigma)", currLabel.Data());
            continue;
          }
          currName.Append("_ratio");
          TH1* histoRatio = (TH1*)projectHisto->Clone(currName.Data());
          histoRatio->Sumw2();
          histoRatio->Divide(refCutHisto);
          histoRatio->SetTitle("");
          histoRatio->GetYaxis()->SetTitle(legTitle.Data());
          histoRatio->GetXaxis()->SetLabelSize(0.04/fracOfHeight);
          histoRatio->GetXaxis()->SetTitleSize(0.035/fracOfHeight);
          histoRatio->GetYaxis()->SetLabelSize(0.03/fracOfHeight);
          histoRatio->GetYaxis()->SetTitleSize(0.03/fracOfHeight);
          histoRatio->GetXaxis()->SetTitleOffset(1.);
          histoRatio->GetYaxis()->SetTitleOffset(0.6);
          histoRatio->GetYaxis()->SetRangeUser(0.5,1.5);

          can->cd(2);
          drawOpt = "e";
          if ( gPad->GetListOfPrimitives()->GetEntries() != 0 ) drawOpt += "same";
          histoRatio->Draw(drawOpt.Data());
        }// loop on sigma cuts
        can->cd(1);
        leg->Draw("same");
      } // loop on theta abs
    } // loop on src
  } // loop on histo type

  igroup1++;
  igroup2=0;
  Double_t ptMin[] = {0., 2., 4., 15.};
  Int_t nPtMins = sizeof(ptMin)/sizeof(ptMin[0]);
  for ( Int_t iptmin=0; iptmin<nPtMins; ++iptmin) {
    for ( Int_t isrc=0; isrc<kNtrackSources; ++isrc ) {
      can = 0x0;
      for ( Int_t itheta=0; itheta<kNthetaAbs; ++itheta ) {
        TLegend* leg = 0x0;
        histoName = GetHistoName(kSigmaVsPt, itheta, isrc);
        for ( Int_t icent=1; icent<=fCentralityClasses->GetNbins(); ++icent ) {
          TH2* histo = (TH2*)GetSum(physSel, trigClassName, fCentralityClasses->GetBinLabel(icent), histoName);
          if ( ! histo ) continue;
          Int_t ptMinBin = histo->GetXaxis()->FindBin(ptMin[iptmin]);
          Int_t ptMaxBin = histo->GetXaxis()->GetNbins()+1;
          currName = histo->GetName();
          currName += Form("_%s_ptMin%i", fCentralityClasses->GetBinLabel(icent), TMath::Nint(ptMin[iptmin]));
          TH1* projectHisto = histo->ProjectionY(currName.Data(), ptMinBin, ptMaxBin);
          if ( projectHisto->GetMaximum() < 50. ) {
            delete projectHisto;
            continue;
          }
          if ( ! can ) {
            currName = Form("CutEffect_%s_ptMin%i", fSrcKeys->At(isrc)->GetName(), TMath::Nint(ptMin[iptmin]));
            can = new TCanvas(currName.Data(), currName.Data(), igroup1*xshift,igroup2*yshift, 600, 600);
            can->Divide(2,1);
          }
          if ( ! leg ) {
            leg = new TLegend(0.6,0.6,0.9,0.9);
            leg->SetBorderSize(1);
          }
          can->cd(itheta+1);
          projectHisto->SetTitle(Form("Cut effect %s %s %g < p_{t} < %g (GeV/c)", fSrcKeys->At(isrc)->GetName(), fThetaAbsKeys->At(itheta)->GetName(), histo->GetXaxis()->GetBinLowEdge(ptMinBin), histo->GetXaxis()->GetBinUpEdge(ptMaxBin)));
          projectHisto->SetLineColor(cutColors[icent-1]);
          projectHisto->SetMarkerColor(cutColors[icent-1]);
          projectHisto->SetMarkerStyle(19+icent);
          projectHisto->SetStats(0);
          drawOpt = "e";
          if ( gPad->GetListOfPrimitives()->GetEntries() != 0 ) drawOpt += "same";
          projectHisto->Draw(drawOpt.Data());
          leg->AddEntry(projectHisto, fCentralityClasses->GetBinLabel(icent), "lp");
          Double_t keptEvents = projectHisto->GetBinContent(orderCuts[0]);
          Double_t totalEvents = projectHisto->GetBinContent(orderCuts[nSigmaCuts-1]);
          Double_t accepted = ( totalEvents > 0. ) ? keptEvents / totalEvents : 1.;
          Double_t acceptedErr = ( totalEvents > 0. ) ? TMath::Sqrt(accepted*(1.-accepted)/totalEvents) : 1.;
          printf("%12s %11s %6s (pt>%g) rejected evts : %6.2f +- %6.3f %%\n", fSrcKeys->At(isrc)->GetName(), fThetaAbsKeys->At(itheta)->GetName(), fCentralityClasses->GetBinLabel(icent), ptMin[iptmin], (1.-accepted)*100., acceptedErr*100.);
          //printf("  rejected %g  total %g   (%s vs %s)\n",totalEvents-keptEvents,totalEvents,projectHisto->GetXaxis()->GetBinLabel(orderCuts[0]),projectHisto->GetXaxis()->GetBinLabel(orderCuts[nSigmaCuts-1]));
        } // loop on centrality
        if ( leg ) leg->Draw("same");
      } // loop on theta abs
    } // loop on sources
  } // loop on pt min
  delete [] orderCuts;
}
