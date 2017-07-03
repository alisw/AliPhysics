/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
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

///////////////////////////////////////////////////////////////////////////
// AliAnalysisTaskSE for the extraction of the various histograms to
// study the pt spectra of identified hadrons:
// - log(dEdx)-log(dEdxBB) distributions for pions, kaons and protons in pt bins
// - Pt distributions of pions, kaons and protons with nSigma PID
// Authors:
// E. Biolcati, biolcati@to.infn.it
// L. Milano, milano@to.infn.it
// F. Prino, prino@to.infn.it
// Y. Corrales, corrales@to.infn.it
// N. Jacazio, jacazio@to.infn.it
///////////////////////////////////////////////////////////////////////////

#include <Rtypes.h>
#include <TChain.h>
#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <THnSparse.h>
#include <TNtuple.h>
#include <TParticle.h>
#include <TRandom3.h>


#include "AliAnalysisManager.h"
#include "AliAnalysisUtils.h"
#include "AliCentrality.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliESDUtils.h"
#include "AliITSPIDResponse.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliMultiplicity.h"
#include "AliMultSelection.h"
#include "AliPID.h"
//#include "AliStack.h"

#include "AliAnalysisTaskSEITSsaSpectra.h"

ClassImp(AliAnalysisTaskSEITSsaSpectra)
/* $Id$ */

//
//
//________________________________________________________________________
AliAnalysisTaskSEITSsaSpectra::AliAnalysisTaskSEITSsaSpectra():
  AliAnalysisTaskSE("TaskITSsaSpectra"),
  fESD(NULL),
  fITSPidParams(NULL),
  fITSPIDResponse(NULL),
  fEventCuts(0),
  fOutput(NULL),
  fListCuts(NULL),
  fListTree(NULL),
  fListPriors(NULL),
  fDCAxyCutFunc(NULL),
  fDCAzCutFunc(NULL),
  fNtupleData(NULL),
  fNtupleMC(NULL),
  fHistNEvents(NULL),
  fHistMCEvents(NULL),
  fHistMultBefEvtSel(NULL),
  fHistMultAftEvtSel(NULL),
  fHistVtxZ(NULL),
  fHistDEDX(NULL),
  fHistDEDXdouble(NULL),
  fTriggerSel(AliVEvent::kMB),
  fMaxVtxZCut(10.),
  fChkIsSDDIn(kTRUE),
  fRejIncDAQ(kTRUE),
  fDoSPDCvsTCut(kTRUE),
  fChkVtxSPDRes(kTRUE),
  fChkVtxZSep(kFALSE),
  fReqBothVtx(kFALSE),
  fExtEventCuts(kFALSE),
  fMinSPDPts(1),
  fMinNdEdxSamples(3),
  fAbsEtaCut(.8),
  fMinRapCut(-.5),
  fMaxRapCut( .5),
  fCMSRapFct(.0),
  fMindEdx(0.),
  fMinNSigma(1.5),
  fMaxChi2Clu(2.5),
  fNSigmaDCAxy(7.),
  fNSigmaDCAz(7.),
  fMultMethod(0),
  fMultEstimator("V0M"),
  fLowMult(-1.),
  fUpMult(-1.),
  fEvtMult(-999),
  fYear(2010),
  fPidMethod(kMeanCut),
  fUseDefaultPriors(kTRUE),
  fIsMC(kFALSE),
  fFillNtuple(kFALSE),
  fFillIntDistHist(kFALSE),
  fPlpType(kNoPileup),
  fMinPlpContribMV(5),
  fMaxPlpChi2MV(5.),
  fMinWDistMV(15.),
  fCheckPlpFromDifferentBCMV(kFALSE),
  fMinPlpContribSPD(5),
  fMinPlpZdistSPD(.8),
  fRandGener(0x0),
  fSmearMC(kFALSE),
  fSmearP(0.),
  fSmeardEdx(0.)
{
  //Constructor
  fRandGener = new TRandom3(0);

  Double_t xbins[kNbins + 1] = {0.08, 0.10, 0.12, 0.14, 0.16, 0.18, 0.20, 0.25,
                                0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65,
                                0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.0
                               };
  for (Int_t iBin = 0; iBin < (kNbins + 1); ++iBin)
    fPtBinLimits[iBin] = xbins[iBin];

  for (int iChg = 0; iChg < kNchg; ++iChg)
    fHistNTracks[iChg] = NULL;

  for (Int_t iSpc = 0; iSpc < kNspc; ++iSpc) {
    for (Int_t iChg = 0; iChg < kNchg; ++iChg) {
      int index = iSpc * kNchg + iChg;

      fHistNSigmaSep[index]  = NULL;

      fHistPrimMCGenVtxZall[index] = NULL;
      fHistPrimMCGenVtxZcut[index] = NULL;

      fHistReco      [index] = NULL;
      fHistMCReco    [index] = NULL;
      fHistPrimMCReco[index] = NULL;

      fHistPrimMCReco[index] = NULL;
      fHistSstrMCReco[index] = NULL;
      fHistSmatMCReco[index] = NULL;

      for (Int_t iptbin = 0; iptbin < kNbins; ++iptbin) {
        fHistRecoDCA[index][iptbin] = NULL; //! histo with DCA distibution

        fHistPrimDCA[index][iptbin] = NULL; //! histo with DCA distibution and dedx PID
        fHistSstrDCA[index][iptbin] = NULL; //! histo with DCA distibution and dedx PID
        fHistSmatDCA[index][iptbin] = NULL; //! histo with DCA distibution and dedx PID

        fHistMCtruthPrimDCA[index][iptbin] = NULL; //! histo with DCA distibution, MC truth
        fHistMCtruthSstrDCA[index][iptbin] = NULL; //! histo with DCA distibution, MC truth
        fHistMCtruthSmatDCA[index][iptbin] = NULL; //! histo with DCA distibution, MC truth
      }
    }
  }
  for (Int_t j = 0; j < 4; j++) fHistCharge[j] = NULL;

  for (Int_t iptbin = 0; iptbin < kNbins; ++iptbin) {
    //dEdx distributions
    fHistPosHypPi[iptbin] = NULL;
    fHistPosHypKa[iptbin] = NULL;
    fHistPosHypPr[iptbin] = NULL;
    fHistNegHypPi[iptbin] = NULL;
    fHistNegHypKa[iptbin] = NULL;
    fHistNegHypPr[iptbin] = NULL;

    //dEdx distributions for MC
    fHistMCPosOtherHypPion[iptbin] = NULL;
    fHistMCPosOtherHypKaon[iptbin] = NULL;
    fHistMCPosOtherHypProt[iptbin] = NULL;
    fHistMCPosElHypPion[iptbin]    = NULL;
    fHistMCPosElHypKaon[iptbin]    = NULL;
    fHistMCPosElHypProt[iptbin]    = NULL;
    fHistMCPosPiHypPion[iptbin]    = NULL;
    fHistMCPosPiHypKaon[iptbin]    = NULL;
    fHistMCPosPiHypProt[iptbin]    = NULL;
    fHistMCPosKaHypPion[iptbin]    = NULL;
    fHistMCPosKaHypKaon[iptbin]    = NULL;
    fHistMCPosKaHypProt[iptbin]    = NULL;
    fHistMCPosPrHypPion[iptbin]    = NULL;
    fHistMCPosPrHypKaon[iptbin]    = NULL;
    fHistMCPosPrHypProt[iptbin]    = NULL;

    fHistMCNegOtherHypPion[iptbin] = NULL;
    fHistMCNegOtherHypKaon[iptbin] = NULL;
    fHistMCNegOtherHypProt[iptbin] = NULL;
    fHistMCNegElHypPion[iptbin]    = NULL;
    fHistMCNegElHypKaon[iptbin]    = NULL;
    fHistMCNegElHypProt[iptbin]    = NULL;
    fHistMCNegPiHypPion[iptbin]    = NULL;
    fHistMCNegPiHypKaon[iptbin]    = NULL;
    fHistMCNegPiHypProt[iptbin]    = NULL;
    fHistMCNegKaHypPion[iptbin]    = NULL;
    fHistMCNegKaHypKaon[iptbin]    = NULL;
    fHistMCNegKaHypProt[iptbin]    = NULL;
    fHistMCNegPrHypPion[iptbin]    = NULL;
    fHistMCNegPrHypKaon[iptbin]    = NULL;
    fHistMCNegPrHypProt[iptbin]    = NULL;
  }

  DefineInput(0, TChain::Class());
  if (!fUseDefaultPriors) DefineInput(1, TList::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
  if (fFillNtuple) DefineOutput(3, TList::Class());

  AliInfo("End of AliAnalysisTaskSEITSsaSpectra");
}

//
//
//___________________________________________________________________________
AliAnalysisTaskSEITSsaSpectra::~AliAnalysisTaskSEITSsaSpectra()
{
  // Destructor in case not running on proof
  if (!AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
    if (fOutput) {
      delete fOutput;
      fOutput = NULL;
    }
    if (fListCuts) {
      delete fListCuts;
      fListCuts = NULL;
    }
    if (fListTree) {
      delete fListTree;
      fListTree = NULL;
    }
    if (fRandGener) {
      delete fRandGener;
      fRandGener = NULL;
    }
    if (fITSPIDResponse) {
      delete fITSPIDResponse;
      fITSPIDResponse = NULL;
    }
  }

  AliInfo("End of AliAnalysisTaskSEITSsaSpectra destructor");
}

//
//
//________________________________________________________________________
void AliAnalysisTaskSEITSsaSpectra::UserCreateOutputObjects()
{
  // Create a TList with histograms and a TNtuple
  // Called once
  if (fFillNtuple) {
    fListTree = new TList();
    fListTree->SetOwner();

    fNtupleData = new TNtuple("fNtupleData", "fNtupleData", "p:pt:s0:s1:s2:s3:dEdx:sign:eta:dcaXY:dcaZ:clumap");
    fListTree->Add(fNtupleData);
    fNtupleMC   = new TNtuple("fNtupleMC", "fNtupleMC", "mcPt:pdgcode:sign:mcEta:mcRap:isph:run");
    fListTree->Add(fNtupleMC);
  }

  if (!fUseDefaultPriors) {
    fListPriors = dynamic_cast<TList*>(GetInputData(1));  //FIXME
    if (!fListPriors) {
      AliDebug(3, "Errors use priors set but no priors list found");
      PostAllData();
      return;
    }
  }

  fOutput = new TList();
  fOutput->SetOwner();
  fOutput->SetName("Spiderman");

  if (fExtEventCuts) { //configure AliEventCuts
    //fEventCuts.SetManualMode();
    //Histograms for event Selection
    fEventCuts.AddQAplotsToList(fOutput);
  }

  TString plpName[3] = {"Sel", "SPD", "MV"};
  fHistNEvents = new TH1I("fHistNEvents", "Number of processed events;Ev. Sel. Step;Counts", kNEvtCuts, .5, (kNEvtCuts + .5));
  fHistNEvents->Sumw2();
  fHistNEvents->SetMinimum(0);
  fHistNEvents->GetXaxis()->SetBinLabel(kIsReadable,  "Readable");
  fHistNEvents->GetXaxis()->SetBinLabel(kIsNotIncDAQ, "PassIncDAQ");
  fHistNEvents->GetXaxis()->SetBinLabel(kPassTrig, "PassPhysSelTrig");
  fHistNEvents->GetXaxis()->SetBinLabel(kIsNotPileup, Form("IsNotPileup_%s", plpName[fPlpType].Data()));
  fHistNEvents->GetXaxis()->SetBinLabel(kCorrelation, "Correlation");
  fHistNEvents->GetXaxis()->SetBinLabel(kHasRecVtx,   "HasVertex");
  fHistNEvents->GetXaxis()->SetBinLabel(kHasGoodVtxZ, "HasGoodVertex");
  fHistNEvents->GetXaxis()->SetBinLabel(kIsSDDIn,     "HasSDDIn");
  fHistNEvents->GetXaxis()->SetBinLabel(kPassSPDclsVsTCut, "PassClsVsTrackletBG");
  fHistNEvents->GetXaxis()->SetBinLabel(kPassMultSel, "PassMultSel");
  fHistNEvents->GetXaxis()->SetBinLabel(kNEvtCuts,    "IsSelected");
  fOutput->Add(fHistNEvents);

  fHistMCEvents = new TH1I("fHistMCEvents", "Number of processed events;Ev. Sel. Step;Counts", kNEvtCuts, .5, (kNEvtCuts + .5));
  fHistMCEvents->Sumw2();
  fHistMCEvents->SetMinimum(0);
  fHistMCEvents->GetXaxis()->SetBinLabel(kIsReadable, "Readable");
  fHistMCEvents->GetXaxis()->SetBinLabel(kPassTrig, "PassPhysSelTrig");
  fHistMCEvents->GetXaxis()->SetBinLabel(kIsNotIncDAQ, "PassIncDAQ");
  fHistMCEvents->GetXaxis()->SetBinLabel(kIsNotPileup, Form("IsNotPileup_%s", plpName[fPlpType].Data()));
  fHistMCEvents->GetXaxis()->SetBinLabel(kCorrelation, "Correlation");
  fHistMCEvents->GetXaxis()->SetBinLabel(kHasRecVtx,   "HasVertex");
  fHistMCEvents->GetXaxis()->SetBinLabel(kHasGoodVtxZ, "HasGoodVertex");
  fHistMCEvents->GetXaxis()->SetBinLabel(kIsSDDIn,     "HasSDDIn");
  fHistMCEvents->GetXaxis()->SetBinLabel(kPassSPDclsVsTCut, "PassClsVsTrackletBG");
  fHistMCEvents->GetXaxis()->SetBinLabel(kPassMultSel, "PassMultSel");
  fHistMCEvents->GetXaxis()->SetBinLabel(kNEvtCuts,    "IsSelected");
  fOutput->Add(fHistMCEvents);

  Int_t kNMultBin = 115;
  Float_t lMultBinLimit[kNMultBin];
  for (Int_t ibin = 0; ibin < kNMultBin; ++ibin)
    lMultBinLimit[ibin] = (ibin < 100) ? (0. + (1.) * ibin) : (100. + 100 * (ibin - 100));

  fHistMultBefEvtSel = new TH1F("fHistMultBefEvtSel", "Event Multiplicity before event selection", kNMultBin - 1, lMultBinLimit);
  fHistMultBefEvtSel->Sumw2();
  fHistMultBefEvtSel->SetMinimum(0);
  fOutput->Add(fHistMultBefEvtSel);

  fHistMultAftEvtSel = new TH1F("fHistMultAftEvtSel", "Event Multiplicity after event selection", kNMultBin - 1, lMultBinLimit);
  fHistMultAftEvtSel->Sumw2();
  fHistMultAftEvtSel->SetMinimum(0);
  fOutput->Add(fHistMultAftEvtSel);

  fHistVtxZ = new TH1F("fHistVtxZ", "Vtx Z distribution", 400, -20, 20);
  fHistVtxZ->Sumw2();
  fHistVtxZ->SetMinimum(0);
  fOutput->Add(fHistVtxZ);

  const char spcName[kNspc][3] = {"Pi", "Ka", "Pr"};
  const char chgName[kNchg][4] = {"Pos", "Neg"};

  char* histName;
  //Histo with track cuts
  for (Int_t iChg = 0; iChg < kNchg; ++iChg) {
    histName = Form("fHistNTracks%s", chgName[iChg]);
    fHistNTracks[iChg] = new TH2F(histName, "Number of ITSsa tracks; Trk Sel; pt [GeV/c]",
                                  20, .5, 20.5, kNbins, fPtBinLimits);
    fHistNTracks[iChg]->Sumw2();

    TString label("no selection");//1
    fHistNTracks[iChg]->GetXaxis()->SetBinLabel(kHasNoSelection, label.Data());
    label = "ITSsa";//2
    fHistNTracks[iChg]->GetXaxis()->SetBinLabel(kIsITSsa, label.Data());
    label = "ITSrefit";//3
    fHistNTracks[iChg]->GetXaxis()->SetBinLabel(kIsITSrefit, label.Data());
    label = "neutral particle";//4
    fHistNTracks[iChg]->GetXaxis()->SetBinLabel(kIsNotNeutralParticle, label.Data());
    label = "SPDcls";//7
    fHistNTracks[iChg]->GetXaxis()->SetBinLabel(kPassSPD, label.Data());
    label = "SDD+SSD cls";//8
    fHistNTracks[iChg]->GetXaxis()->SetBinLabel(kPassPIDcls, label.Data());
    label = "chi2/ncls";//9
    fHistNTracks[iChg]->GetXaxis()->SetBinLabel(kPassChi2Ncls, label.Data());
    label = "eta";//11
    fHistNTracks[iChg]->GetXaxis()->SetBinLabel(kIsInEta, label.Data());
    label = "dE/dx < 0";//12
    fHistNTracks[iChg]->GetXaxis()->SetBinLabel(kPassdEdx, label.Data());
    label = "Pt cut";
    fHistNTracks[iChg]->GetXaxis()->SetBinLabel(kPassPtCut, label.Data());
    label = "DCAz";//13
    fHistNTracks[iChg]->GetXaxis()->SetBinLabel(kPassDCAzcut, label.Data());
    label = "DCAxy";//14
    fHistNTracks[iChg]->GetXaxis()->SetBinLabel(kPassDCAxycut, label.Data());
    fOutput->Add(fHistNTracks[iChg]);
  }

  //binning for the histogram
  const Int_t hnbins = 400;
  Double_t hxmin = 0.01;
  Double_t hxmax = 10;
  Double_t hlogxmin = TMath::Log10(hxmin);
  Double_t hlogxmax = TMath::Log10(hxmax);
  Double_t hbinwidth = (hlogxmax - hlogxmin) / hnbins;
  Double_t hxbins[hnbins + 1];
  hxbins[0] = 0.01;
  for (Int_t i = 1; i <= hnbins; i++) {
    hxbins[i] = hxmin + TMath::Power(10, hlogxmin + i * hbinwidth);
  }

  fHistDEDX = new TH2F("fHistDEDX", "", hnbins, hxbins, 900, 0, 1000);
  fOutput->Add(fHistDEDX);

  fHistDEDXdouble = new TH2F("fHistDEDXdouble", "", 500, -5, 5, 900, 0, 1000);
  fOutput->Add(fHistDEDXdouble);

  for (Int_t ispc = 0; ispc < kNspc; ++ispc) {
    for (Int_t ichg = 0; ichg < kNchg; ++ichg) {
      Int_t index = ispc * kNchg + ichg;

      histName = Form("fHistNSigmaSep%s%s", spcName[ispc], chgName[ichg]);
      fHistNSigmaSep[index] = new TH2F(histName, histName, hnbins, hxbins, 1000, -10, 10);
      fOutput->Add(fHistNSigmaSep[index]);

      //Reconstructed
      histName = Form("fHistReco%s%s", spcName[ispc], chgName[ichg]);
      fHistReco[index] = new TH1F(histName, histName, kNbins, fPtBinLimits);
      fOutput->Add(fHistReco[index]);

      if (fIsMC) {
        //
        //Histograms MC part Gen bef and afte all selection Good Vertex Gen.
        histName = Form("fHistPrimMCGenVtxZall%s%s", spcName[ispc], chgName[ichg]);
        fHistPrimMCGenVtxZall[index] = new TH2F(histName, histName, kNEvtCuts, .5, (kNEvtCuts + .5), kNbins, fPtBinLimits);

        histName = Form("fHistPrimMCGenVtxZcut%s%s", spcName[ispc], chgName[ichg]);
        fHistPrimMCGenVtxZcut[index] = new TH2F(histName, histName, kNEvtCuts, .5, (kNEvtCuts + .5), kNbins, fPtBinLimits);

        fOutput->Add(fHistPrimMCGenVtxZall[index]);
        fOutput->Add(fHistPrimMCGenVtxZcut[index]);

        histName = Form("fHistMCReco%s%s", spcName[ispc], chgName[ichg]);
        fHistMCReco[index] = new TH2F(histName, histName, kNbins, fPtBinLimits, 4, -1.5, 2.5);
        fOutput->Add(fHistMCReco[index]);

        histName = Form("fHistMCPrimReco%s%s", spcName[ispc], chgName[ichg]);
        fHistMCPrimReco[index] = new TH2F(histName, histName, kNbins, fPtBinLimits, 4, -1.5, 2.5);
        fOutput->Add(fHistMCPrimReco[index]);

        //
        //Histograms MC part Rec.
        histName = Form("fHistPrimMCReco%s%s", spcName[ispc], chgName[ichg]);
        fHistPrimMCReco[index] = new TH1F(histName, histName, kNbins, fPtBinLimits);

        histName = Form("fHistSstrMCReco%s%s", spcName[ispc], chgName[ichg]);
        fHistSstrMCReco[index] = new TH1F(histName, histName, kNbins, fPtBinLimits);

        histName = Form("fHistSmatMCReco%s%s", spcName[ispc], chgName[ichg]);
        fHistSmatMCReco[index] = new TH1F(histName, histName, kNbins, fPtBinLimits);

        fOutput->Add(fHistPrimMCReco[index]);
        fOutput->Add(fHistSstrMCReco[index]);
        fOutput->Add(fHistSmatMCReco[index]);
      }//end IsMC

      for (Int_t iptbin = 0; iptbin < kNbins; iptbin++) {
        histName = Form("fHistRecoDCA%s%s_%d", spcName[ispc], chgName[ichg], iptbin);
        fHistRecoDCA[index][iptbin] = new TH1F(histName, histName, 2000, -2, 2);
        fOutput->Add(fHistRecoDCA[index][iptbin]);

        if (fIsMC) {
          histName = Form("fHistPrimDCA%s%s_%d", spcName[ispc], chgName[ichg], iptbin);
          fHistPrimDCA[index][iptbin] = new TH1F(histName, histName, 2000, -2, 2);

          histName = Form("fHistSstrDCA%s%s_%d", spcName[ispc], chgName[ichg], iptbin);
          fHistSstrDCA[index][iptbin] = new TH1F(histName, histName, 2000, -2, 2);

          histName = Form("fHistSmatDCA%s%s_%d", spcName[ispc], chgName[ichg], iptbin);
          fHistSmatDCA[index][iptbin] = new TH1F(histName, histName, 2000, -2, 2);

          histName = Form("fHistMCtruthPrimDCA%s%s_%d", spcName[ispc], chgName[ichg], iptbin);
          fHistMCtruthPrimDCA[index][iptbin] = new TH1F(histName, histName, 2000, -2, 2);

          histName = Form("fHistMCtruthSstrDCA%s%s_%d", spcName[ispc], chgName[ichg], iptbin);
          fHistMCtruthSstrDCA[index][iptbin] = new TH1F(histName, histName, 2000, -2, 2);

          histName = Form("fHistMCtruthSmatDCA%s%s_%d", spcName[ispc], chgName[ichg], iptbin);
          fHistMCtruthSmatDCA[index][iptbin] = new TH1F(histName, histName, 2000, -2, 2);

          fOutput->Add(fHistPrimDCA[index][iptbin]);
          fOutput->Add(fHistSstrDCA[index][iptbin]);
          fOutput->Add(fHistSmatDCA[index][iptbin]);
          fOutput->Add(fHistMCtruthPrimDCA[index][iptbin]);
          fOutput->Add(fHistMCtruthSstrDCA[index][iptbin]);
          fOutput->Add(fHistMCtruthSmatDCA[index][iptbin]);
        }//end isMC
      }
    }
  }

  for (Int_t i = 0; i < 4; i++) {
    fHistCharge[i] = new TH1F(Form("fHistChargeLay%d", i), Form("fHistChargeLay%d", i), 100, 0, 300);
    fOutput->Add(fHistCharge[i]);
  }

  if (fFillIntDistHist) {
    for (Int_t i = 0; i < kNbins; i++) {
      //dEdx distributions
      fHistPosHypPi[i] = new TH1F(Form("fHistPosHypPi%d", i), Form("fHistPosHypPi%d", i), 2000, -1, 1); //DCA distr. with NSigma PID
      fHistPosHypKa[i] = new TH1F(Form("fHistPosHypKa%d", i), Form("fHistPosHypKa%d", i), 2000, -1, 1);
      fHistPosHypPr[i] = new TH1F(Form("fHistPosHypPr%d", i), Form("fHistPosHypPr%d", i), 2000, -1, 1);
      fHistNegHypPi[i] = new TH1F(Form("fHistNegHypPi%d", i), Form("fHistNegHypPi%d", i), 2000, -1, 1);
      fHistNegHypKa[i] = new TH1F(Form("fHistNegHypKa%d", i), Form("fHistNegHypKa%d", i), 2000, -1, 1);
      fHistNegHypPr[i] = new TH1F(Form("fHistNegHypPr%d", i), Form("fHistNegHypPr%d", i), 2000, -1, 1);

      fOutput->Add(fHistPosHypPi[i]); //DCA distr
      fOutput->Add(fHistPosHypKa[i]);
      fOutput->Add(fHistPosHypPr[i]);
      fOutput->Add(fHistNegHypPi[i]);
      fOutput->Add(fHistNegHypKa[i]);
      fOutput->Add(fHistNegHypPr[i]);

      if (fIsMC) {
        fHistMCPosOtherHypPion[i] = new TH1F(Form("fHistMCPosOtherHypPion%d", i), Form("fHistMCPosOtherHypPion%d", i), 175, -3.5, 3.5); //MC truth
        fHistMCPosOtherHypKaon[i] = new TH1F(Form("fHistMCPosOtherHypKaon%d", i), Form("fHistMCPosOtherHypKaon%d", i), 175, -3.5, 3.5);
        fHistMCPosOtherHypProt[i] = new TH1F(Form("fHistMCPosOtherHypProt%d", i), Form("fHistMCPosOtherHypProt%d", i), 175, -3.5, 3.5);
        fHistMCPosElHypPion[i]    = new TH1F(Form("fHistMCPosElHypPion%d", i), Form("fHistMCPosElHypPion%d", i), 175, -3.5, 3.5);
        fHistMCPosElHypKaon[i]    = new TH1F(Form("fHistMCPosElHypKaon%d", i), Form("fHistMCPosElHypKaon%d", i), 175, -3.5, 3.5);
        fHistMCPosElHypProt[i]    = new TH1F(Form("fHistMCPosElHypProt%d", i), Form("fHistMCPosElHypProt%d", i), 175, -3.5, 3.5);
        fHistMCPosPiHypPion[i]    = new TH1F(Form("fHistMCPosPiHypPion%d", i), Form("fHistMCPosPiHypPion%d", i), 175, -3.5, 3.5);
        fHistMCPosPiHypKaon[i]    = new TH1F(Form("fHistMCPosPiHypKaon%d", i), Form("fHistMCPosPiHypKaon%d", i), 175, -3.5, 3.5);
        fHistMCPosPiHypProt[i]    = new TH1F(Form("fHistMCPosPiHypProton%d", i), Form("fHistMCPosPiHypProton%d", i), 175, -3.5, 3.5);
        fHistMCPosKaHypPion[i]    = new TH1F(Form("fHistMCPosKaHypPion%d", i), Form("fHistMCPosKaHypPion%d", i), 175, -3.5, 3.5);
        fHistMCPosKaHypKaon[i]    = new TH1F(Form("fHistMCPosKaHypKaon%d", i), Form("fHistMCPosKaHypKaon%d", i), 175, -3.5, 3.5);
        fHistMCPosKaHypProt[i]    = new TH1F(Form("fHistMCPosKaHypProt%d", i), Form("fHistMCPosKaHypProt%d", i), 175, -3.5, 3.5);
        fHistMCPosPrHypPion[i]    = new TH1F(Form("fHistMCPosPrHypPion%d", i), Form("fHistMCPosPrHypPion%d", i), 175, -3.5, 3.5);
        fHistMCPosPrHypKaon[i]    = new TH1F(Form("fHistMCPosPrHypKaon%d", i), Form("fHistMCPosPrHypKaon%d", i), 175, -3.5, 3.5);
        fHistMCPosPrHypProt[i]    = new TH1F(Form("fHistMCPosPrHypProt%d", i), Form("fHistMCPosPrHypProt%d", i), 175, -3.5, 3.5);

        fHistMCNegOtherHypPion[i] = new TH1F(Form("fHistMCNegOtherHypPion%d", i), Form("fHistMCNegOtherHypPion%d", i), 175, -3.5, 3.5); //MC truth
        fHistMCNegOtherHypKaon[i] = new TH1F(Form("fHistMCNegOtherHypKaon%d", i), Form("fHistMCNegOtherHypKaon%d", i), 175, -3.5, 3.5);
        fHistMCNegOtherHypProt[i] = new TH1F(Form("fHistMCNegOtherHypProt%d", i), Form("fHistMCNegOtherHypProt%d", i), 175, -3.5, 3.5);
        fHistMCNegElHypPion[i]    = new TH1F(Form("fHistMCNegElHypPion%d", i), Form("fHistMCNegElHypPion%d", i), 175, -3.5, 3.5);
        fHistMCNegElHypKaon[i]    = new TH1F(Form("fHistMCNegElHypKaon%d", i), Form("fHistMCNegElHypKaon%d", i), 175, -3.5, 3.5);
        fHistMCNegElHypProt[i]    = new TH1F(Form("fHistMCNegElHypProt%d", i), Form("fHistMCNegElHypProt%d", i), 175, -3.5, 3.5);
        fHistMCNegPiHypPion[i]    = new TH1F(Form("fHistMCNegPiHypPion%d", i), Form("fHistMCNegPiHypPion%d", i), 175, -3.5, 3.5);
        fHistMCNegPiHypKaon[i]    = new TH1F(Form("fHistMCNegPiHypKaon%d", i), Form("fHistMCNegPiHypKaon%d", i), 175, -3.5, 3.5);
        fHistMCNegPiHypProt[i]    = new TH1F(Form("fHistMCNegPiHypProt%d", i), Form("fHistMCNegPiHypProt%d", i), 175, -3.5, 3.5);
        fHistMCNegKaHypPion[i]    = new TH1F(Form("fHistMCNegKaHypPion%d", i), Form("fHistMCNegKaHypPion%d", i), 175, -3.5, 3.5);
        fHistMCNegKaHypKaon[i]    = new TH1F(Form("fHistMCNegKaHypKaon%d", i), Form("fHistMCNegKaHypKaon%d", i), 175, -3.5, 3.5);
        fHistMCNegKaHypProt[i]    = new TH1F(Form("fHistMCNegKaHypProt%d", i), Form("fHistMCNegKaHypProt%d", i), 175, -3.5, 3.5);
        fHistMCNegPrHypPion[i]    = new TH1F(Form("fHistMCNegPrHypPion%d", i), Form("fHistMCNegPrHypPion%d", i), 175, -3.5, 3.5);
        fHistMCNegPrHypKaon[i]    = new TH1F(Form("fHistMCNegPrHypKaon%d", i), Form("fHistMCNegPrHypKaon%d", i), 175, -3.5, 3.5);
        fHistMCNegPrHypProt[i]    = new TH1F(Form("fHistMCNegPrHypProt%d", i), Form("fHistMCNegPrHypProt%d", i), 175, -3.5, 3.5);

        fOutput->Add(fHistMCPosOtherHypPion[i]);//MC truth
        fOutput->Add(fHistMCPosOtherHypKaon[i]);
        fOutput->Add(fHistMCPosOtherHypProt[i]);
        fOutput->Add(fHistMCPosElHypPion[i]);
        fOutput->Add(fHistMCPosElHypKaon[i]);
        fOutput->Add(fHistMCPosElHypProt[i]);
        fOutput->Add(fHistMCPosPiHypPion[i]);
        fOutput->Add(fHistMCPosPiHypKaon[i]);
        fOutput->Add(fHistMCPosPiHypProt[i]);
        fOutput->Add(fHistMCPosKaHypPion[i]);
        fOutput->Add(fHistMCPosKaHypKaon[i]);
        fOutput->Add(fHistMCPosKaHypProt[i]);
        fOutput->Add(fHistMCPosPrHypPion[i]);
        fOutput->Add(fHistMCPosPrHypKaon[i]);
        fOutput->Add(fHistMCPosPrHypProt[i]);

        fOutput->Add(fHistMCNegOtherHypPion[i]);//MC truth
        fOutput->Add(fHistMCNegOtherHypKaon[i]);
        fOutput->Add(fHistMCNegOtherHypProt[i]);
        fOutput->Add(fHistMCNegElHypPion[i]);
        fOutput->Add(fHistMCNegElHypKaon[i]);
        fOutput->Add(fHistMCNegElHypProt[i]);
        fOutput->Add(fHistMCNegPiHypPion[i]);
        fOutput->Add(fHistMCNegPiHypKaon[i]);
        fOutput->Add(fHistMCNegPiHypProt[i]);
        fOutput->Add(fHistMCNegKaHypPion[i]);
        fOutput->Add(fHistMCNegKaHypKaon[i]);
        fOutput->Add(fHistMCNegKaHypProt[i]);
        fOutput->Add(fHistMCNegPrHypPion[i]);
        fOutput->Add(fHistMCNegPrHypKaon[i]);
        fOutput->Add(fHistMCNegPrHypProt[i]);
      }//end IsMC
    }//end FillIntDistHist
  }

  // Post output data container
  PostData(1, fOutput);
  if (fFillNtuple) PostData(3, fListTree);
  AliInfo("End of CreateOutputObjects");
}

//
//
//___________________________________________________________________________
void AliAnalysisTaskSEITSsaSpectra::CreateDCAcutFunctions()
{
  fListCuts = new TList();
  fListCuts->SetOwner();
  Double_t xyP[3];
  Double_t zP[3];
  if (fYear == 2009) {
    if (fIsMC) {
      xyP[0] = 88.63; //MC LHC10a12
      xyP[1] = 19.57;
      xyP[2] = 1.65;
      zP[0] = 140.98;
      zP[1] = 62.33;
      zP[2] = 1.15;
    } else {
      xyP[0] = 85.28; //DATA 900 GeV pass6
      xyP[1] = 25.78;
      xyP[2] = 1.55;
      zP[0] = 146.80;
      zP[1] = 70.07;
      zP[2] = 1.11;
    }
  } else if (fYear == 2010) {
    if (fIsMC) {
      xyP[0] = 36.; //MC LHC10d1
      xyP[1] = 43.9;
      xyP[2] = 1.3;
      zP[0] = 111.9;
      zP[1] = 59.8;
      zP[2] = 1.2;
    } else {
      xyP[0] = 32.7;//DATA 7 TeV pass2
      xyP[1] = 44.8;
      xyP[2] = 1.3;
      zP[0] = 117.3;
      zP[1] = 66.8;
      zP[2] = 1.2;
    }
  }
  fDCAxyCutFunc = new TF1("fDCAxyCutFunc", "[3]*([0]+[1]/TMath::Power(TMath::Abs(x),[2]))", .05, 10.);
  for (Int_t ipar = 0; ipar < 3; ipar++) fDCAxyCutFunc->SetParameter(ipar, xyP[ipar]);
  fDCAxyCutFunc->SetParameter(3, fNSigmaDCAxy);
  fDCAxyCutFunc->SetParName(3, "Sigmas");

  fDCAzCutFunc = new TF1("fDCAzCutFunc", "[3]*([0]+[1]/TMath::Power(TMath::Abs(x),[2]))", .05, 10.);
  for (Int_t ipar = 0; ipar < 3; ipar++) fDCAzCutFunc->SetParameter(ipar, zP[ipar]);
  fDCAzCutFunc->SetParameter(3, fNSigmaDCAz);
  fDCAzCutFunc->SetParName(3, "Sigmas");

  fListCuts->Add(fDCAxyCutFunc);
  fListCuts->Add(fDCAzCutFunc);
}

//
//
//________________________________________________________________________
void AliAnalysisTaskSEITSsaSpectra::Init()
{
  // Initialization
  AliInfo("Tracks selections");
  AliInfoF(" y = yLab + %.3f,  Ymin %.1f, Ymax %.1f, Eabs %.1f, DCAxyCut %.1f, DCAzCut %.1f, Chi2 %.1f,   nSPD %d,   nPID %d",
           fCMSRapFct, fMinRapCut, fMaxRapCut, fAbsEtaCut, fNSigmaDCAxy, fNSigmaDCAz, fMaxChi2Clu, fMinSPDPts, fMinNdEdxSamples);

  if (fMultMethod)
    AliInfoF("Cent. %.f %.f %s", fLowMult, fUpMult, fMultEstimator.Data());

  CreateDCAcutFunctions(); //Creating kParamContainer data
  // Post parameter data container
  PostData(2, fListCuts);
  return;
}

//
//
//________________________________________________________________________
void AliAnalysisTaskSEITSsaSpectra::UserExec(Option_t*)
{
  // Main loop
  // Called for each event
  //
  fESD = dynamic_cast<AliESDEvent*>(InputEvent());
  if (!fESD) {
    AliDebug(3, "fESD not available");
    PostAllData();
    return;
  }

  //Fill some histograms before event selection
  //Check Monte Carlo information and other access first:
  //AliStack*   lMCstack = NULL;
  AliMCEvent* lMCevent = NULL;

  Bool_t lHasGoodVtxGen = kFALSE;
  if (fIsMC) {
    AliMCEventHandler* lMCevtHandler = dynamic_cast<AliMCEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
    if (!lMCevtHandler) {
      AliDebug(3, "Could not retrieve MC event handler");
      PostAllData();
      return;
    }

    lMCevent = lMCevtHandler->MCEvent();
    if (!lMCevent) {
      AliDebug(3, "Could not retrieve MC event");
      PostAllData();
      return;
    }

    //lMCstack = lMCevent->Stack();
    //if (!lMCstack) {
    //  AliDebug(3, "MC stack not available");
    //  PostAllData();
    //  return;
    //}

    const AliVVertex* lGenVtx = lMCevent->GetPrimaryVertex();
    if (!lGenVtx) {
      AliDebug(3, "MC Vtx not available");
      PostAllData();
      return;
    }

    //Selection of the MC sample
    if (!(TMath::Abs(lGenVtx->GetZ()) > fMaxVtxZCut)) {
      lHasGoodVtxGen = kTRUE;
    }
  }

  //Event selection
  EEvtCut_Type lastEvtCutPassed = kIsReadable;
  Bool_t lIsEventSelected = IsEventAccepted(lastEvtCutPassed);
  if (fIsMC) AnalyseMCParticles(lMCevent, lastEvtCutPassed, lHasGoodVtxGen);
  for (int istep = (int)kIsReadable; istep < (int)lastEvtCutPassed; ++istep) {
    fHistNEvents->Fill(istep);
    if (lHasGoodVtxGen) fHistMCEvents->Fill(istep);
  }

  if (!lIsEventSelected)  {
    AliDebugF(3, "Event rejected by event selection after %d step", (int)lastEvtCutPassed);
    PostAllData();
    return;
  }
	fHistNEvents->Fill(kNEvtCuts);
  if (lHasGoodVtxGen) fHistMCEvents->Fill(kNEvtCuts);

  if (fMultMethod) //Fill fHistMultAftEvtSel after the event Selection
    fHistMultAftEvtSel->Fill(fEvtMult);

  if (!fITSPIDResponse)
    fITSPIDResponse = new AliITSPIDResponse(fIsMC);

  TAxis* lAxis = new TAxis(kNbins, fPtBinLimits); //utils

  //loop on tracks
  for (Int_t iTrk = 0; iTrk < fESD->GetNumberOfTracks(); ++iTrk) {
    AliESDtrack* track = dynamic_cast<AliESDtrack*>(fESD->GetTrack(iTrk));
    if (!track) continue;

    if (fIsMC) {
      Int_t trkLabel = TMath::Abs(track->GetLabel());

      TParticle* mcTrk = ((AliMCParticle*)lMCevent->GetTrack(trkLabel))->Particle();
      Int_t pdg = mcTrk->GetPdgCode();
      if (TMath::Abs(pdg) > 1E10) //protection to remove High ionization part
        continue;
    }

    //track selection
    Double_t trkPt  = track->Pt();
    Double_t dEdxLay[4]; track->GetITSdEdxSamples(dEdxLay);
    Double_t dEdx = track->GetITSsignal();

    Int_t   iChg   = (track->GetSign() > 0) ? 0 : 1; //0: Pos; 1:Neg
    ULong_t status = track->GetStatus();
    UInt_t  clumap = track->GetITSClusterMap();

    Int_t nSPD = 0;
    Int_t nPtsForPid = 0;
    for (Int_t il = 0; il < 6; il++) {
      if (TESTBIT(clumap, il)) {
        if (il < 2) nSPD++;
        else nPtsForPid++;
      }
    }

    ETrkCut_Type trkSel = kHasNoSelection;

    //"no selection"
    fHistNTracks[iChg]->Fill(trkSel, trkPt);

    //"ITSsa"
    if (!(status & AliESDtrack::kITSpureSA)) continue;
    trkSel = kIsITSsa;
    fHistNTracks[iChg]->Fill(trkSel, trkPt);

    //"ITSrefit"
    if (!(status & AliESDtrack::kITSrefit)) continue;
    trkSel = kIsITSrefit;
    fHistNTracks[iChg]->Fill(trkSel, trkPt);

    //"neutral particle"
    if (TMath::Abs(track->GetSign()) < 0.0001) continue;
    trkSel = kIsNotNeutralParticle;
    fHistNTracks[iChg]->Fill(trkSel, trkPt);

    //"SPDcls"
    if (nSPD < fMinSPDPts) continue;//At least one point in the SPD
    trkSel = kPassSPD;
    fHistNTracks[iChg]->Fill(trkSel, trkPt);

    //"SDD+SSD cls" at least 3 points on SSD/SDD
    if (nPtsForPid < fMinNdEdxSamples) continue;
    trkSel = kPassPIDcls;
    fHistNTracks[iChg]->Fill(trkSel, trkPt);

    //"chi2/ncls"->chisquare/nclusters
    Int_t nclu = nSPD + nPtsForPid;
    if (track->GetITSchi2() / nclu > fMaxChi2Clu) continue;
    trkSel = kPassChi2Ncls;
    fHistNTracks[iChg]->Fill(trkSel, trkPt);

    //"eta"->pseudorapidity
    if (TMath::Abs(track->Eta()) > fAbsEtaCut) continue;
    trkSel = kIsInEta;
    fHistNTracks[iChg]->Fill(trkSel, trkPt);

    //"dE/dx < 0"->truncated mean
    if (dEdx < 0) continue;
    trkSel = kPassdEdx;
    fHistNTracks[iChg]->Fill(trkSel, trkPt);

    //fill propaganda plot with dedx before pt cut
    fHistDEDX->Fill(track->GetP(), dEdx);
    fHistDEDXdouble->Fill(track->GetP()*track->GetSign(), dEdx);

    //"ptCut"
    Int_t trkPtBin = lAxis->FindFixBin(trkPt);
    if (!trkPtBin || trkPtBin > lAxis->GetNbins()) continue;
    trkPtBin--;
    trkSel = kPassPtCut;
    fHistNTracks[iChg]->Fill(trkSel, trkPt);

    /////////////////////////////////////////////////////////////////////////////
    // Because we fit the DCA distributions, we need them after the DCAz cut,  //
    // otherwise the templates and distributions are modified by this cut      //
    // (if it is done with the one of the DCAxy)                               //
    /////////////////////////////////////////////////////////////////////////////
    Float_t impactXY, impactZ;
    track->GetImpactParameters(impactXY, impactZ);
    //"DCAz"
    if (!DCAcutZ(impactZ, trkPt)) continue;
    trkSel = kPassDCAzcut;
    fHistNTracks[iChg]->Fill(trkSel, trkPt);

    if (fFillNtuple) {
      Float_t xnt[12];
      Int_t index = 0;
      /*1 */ xnt[index++] = (float)track->GetP();
      /*2 */ xnt[index++] = (float)track->Pt();
      /*3 */ xnt[index++] = (float)dEdxLay[0];
      /*4 */ xnt[index++] = (float)dEdxLay[1];
      /*5 */ xnt[index++] = (float)dEdxLay[2];
      /*6 */ xnt[index++] = (float)dEdxLay[3];
      /*7 */ xnt[index++] = (float)dEdx;
      /*8*/ xnt[index++]  = (float)track->GetSign();
      /*9*/ xnt[index++]  = (float)track->Eta();
      /*10*/ xnt[index++] = (float)impactXY;
      /*11*/ xnt[index++] = (float)impactZ;
      /*12*/ xnt[index++] = (float)clumap;

      fNtupleData->Fill(xnt);
    } else {
      //track PID aproach
      Double_t logdiff[4];
      Int_t fPid = GetTrackPid(track, logdiff);
      //End PID approach

      //Compute y
      Double_t y[AliPID::kSPECIES];
      //loop per specie x4
      for (Int_t ispc = 0; ispc < AliPID::kSPECIES; ++ispc) {
        Float_t mass = AliPID::ParticleMass(AliPID::EParticleType(ispc));
        y[ispc] = Eta2y(trkPt, mass, track->Eta());
        y[ispc] += fCMSRapFct;
      }
      Bool_t lIsGoodTrack = ((fPid > AliPID::kMuon) && IsRapIn(y[fPid]));
      Int_t  lPidIndex = lIsGoodTrack ? ((fPid - 2) * kNchg + iChg) : -1;

      Int_t  lMCtrk = -999;
      Int_t  lMCpdg = -999;
      Int_t  lMCspc = AliPID::kElectron;
      if (fIsMC) {
        lMCtrk = TMath::Abs(track->GetLabel());
        TParticle* trkMC = ((AliMCParticle*)lMCevent->GetTrack(lMCtrk))->Particle();
        lMCpdg = trkMC->GetPdgCode();

        //        if (TMath::Abs(lMCpdg) ==   11 && fPid == AliPID::kPion) lMCspc = AliPID::kPion;
        //        if (TMath::Abs(lMCpdg) ==   13 && fPid == AliPID::kPion) lMCspc = AliPID::kPion;
        if (TMath::Abs(lMCpdg) ==  211) lMCspc = AliPID::kPion; // select Pi+/Pi- only
        if (TMath::Abs(lMCpdg) ==  321) lMCspc = AliPID::kKaon; // select K+/K- only
        if (TMath::Abs(lMCpdg) == 2212) lMCspc = AliPID::kProton; // select p+/p- only
      }
      Bool_t lIsGoodPart = ((lMCspc > AliPID::kMuon) && IsRapIn(y[lMCspc]));
      Int_t lMCtIndex = lIsGoodPart ? ((lMCspc - 2) *  kNchg + iChg) : -1;

      if (lIsGoodTrack) {
        //DCA distributions, before the DCA cuts, based on PID approach
        fHistRecoDCA[lPidIndex][trkPtBin]->Fill(impactXY);
        //DCA distributions, before the DCAxy cuts from the MC kinematics
        //Filling DCA distribution with MC truth Physics values
        if (fIsMC) {
          if (lMCevent->IsPhysicalPrimary(lMCtrk))        fHistPrimDCA[lPidIndex][trkPtBin]->Fill(impactXY);
          if (lMCevent->IsSecondaryFromWeakDecay(lMCtrk)) fHistSstrDCA[lPidIndex][trkPtBin]->Fill(impactXY);
          if (lMCevent->IsSecondaryFromMaterial(lMCtrk))  fHistSmatDCA[lPidIndex][trkPtBin]->Fill(impactXY);
        }
      }// end lIsGoodTrack
      if (fIsMC && lIsGoodPart) {
        //DCA distributions, before the DCAxy cuts from the MC kinematics
        //Filling DCA distribution with MC truth Physics values
        if (lMCevent->IsPhysicalPrimary(lMCtrk))        fHistMCtruthPrimDCA[lMCtIndex][trkPtBin]->Fill(impactXY);
        if (lMCevent->IsSecondaryFromWeakDecay(lMCtrk)) fHistMCtruthSstrDCA[lMCtIndex][trkPtBin]->Fill(impactXY);
        if (lMCevent->IsSecondaryFromMaterial(lMCtrk))  fHistMCtruthSmatDCA[lMCtIndex][trkPtBin]->Fill(impactXY);
      }// end lIsGoodPart
      //"DCAxy"
      if (!DCAcutXY(impactXY, trkPt)) continue;
      trkSel = kPassDCAxycut;
      fHistNTracks[iChg]->Fill(trkSel, trkPt);

      if (lIsGoodTrack) fHistReco[lPidIndex]->Fill(trkPt);
      //Filling Histos for Reco Efficiency
      //information from the MC kinematics
      if (fIsMC && lIsGoodPart) {
        if (lMCevent->IsPhysicalPrimary(lMCtrk))         fHistPrimMCReco[lMCtIndex]->Fill(trkPt);
        if (lMCevent->IsSecondaryFromWeakDecay(lMCtrk))  fHistSstrMCReco[lMCtIndex]->Fill(trkPt);
        if (lMCevent->IsSecondaryFromMaterial(lMCtrk))   fHistSmatMCReco[lMCtIndex]->Fill(trkPt);
      }

      Int_t binPart = (lMCspc > AliPID::kMuon) ? (lMCspc - 2) : -1;
      if (fIsMC && lIsGoodTrack) {
        fHistMCReco[lPidIndex]->Fill(trkPt, binPart);
        if (lMCevent->IsPhysicalPrimary(lMCtrk))
          fHistMCPrimReco[lPidIndex]->Fill(trkPt, binPart);
      }//end y

      if (lIsGoodTrack && fFillIntDistHist) {
        //
        //integral approach histograms
        //
        Int_t lAbsPdgCode = TMath::Abs(lMCpdg);
        if (track->GetSign() > 0) {
          if (IsRapIn(y[AliPID::kPion]))  fHistPosHypPi[trkPtBin]->Fill(logdiff[0]);
          if (IsRapIn(y[AliPID::kKaon]))  fHistPosHypKa[trkPtBin]->Fill(logdiff[1]);
          if (IsRapIn(y[AliPID::kProton]))fHistPosHypPr[trkPtBin]->Fill(logdiff[2]);
          if (fIsMC) {
            if (IsRapIn(y[AliPID::kPion])) {
              if ((lAbsPdgCode != 11) && (lAbsPdgCode != 211) && (lAbsPdgCode != 321) && (lAbsPdgCode != 2212))
                fHistMCPosOtherHypPion[trkPtBin]->Fill(logdiff[0]);
              if (lAbsPdgCode == 11)  fHistMCPosElHypPion[trkPtBin]->Fill(logdiff[0]);
              if (lAbsPdgCode == 211) fHistMCPosPiHypPion[trkPtBin]->Fill(logdiff[0]);
              if (lAbsPdgCode == 321) fHistMCPosKaHypPion[trkPtBin]->Fill(logdiff[0]);
              if (lAbsPdgCode == 2212)fHistMCPosPrHypPion[trkPtBin]->Fill(logdiff[0]);
            }
            if (IsRapIn(y[AliPID::kKaon])) {
              if ((lAbsPdgCode != 11) && (lAbsPdgCode != 211) && (lAbsPdgCode != 321) && (lAbsPdgCode != 2212))
                fHistMCPosOtherHypKaon[trkPtBin]->Fill(logdiff[1]);
              if (lAbsPdgCode == 11)  fHistMCPosElHypKaon[trkPtBin]->Fill(logdiff[1]);
              if (lAbsPdgCode == 211) fHistMCPosPiHypKaon[trkPtBin]->Fill(logdiff[1]);
              if (lAbsPdgCode == 321) fHistMCPosKaHypKaon[trkPtBin]->Fill(logdiff[1]);
              if (lAbsPdgCode == 2212)fHistMCPosPrHypKaon[trkPtBin]->Fill(logdiff[1]);
            }
            if (IsRapIn(y[AliPID::kProton])) {
              if ((lAbsPdgCode != 11) && (lAbsPdgCode != 211) && (lAbsPdgCode != 321) && (lAbsPdgCode != 2212))
                fHistMCPosOtherHypProt[trkPtBin]->Fill(logdiff[2]);
              if (lAbsPdgCode == 11)  fHistMCPosElHypProt[trkPtBin]->Fill(logdiff[2]);
              if (lAbsPdgCode == 211) fHistMCPosPiHypProt[trkPtBin]->Fill(logdiff[2]);
              if (lAbsPdgCode == 321) fHistMCPosKaHypProt[trkPtBin]->Fill(logdiff[2]);
              if (lAbsPdgCode == 2212)fHistMCPosPrHypProt[trkPtBin]->Fill(logdiff[2]);
            }
          }
        } else {
          if (IsRapIn(y[AliPID::kPion]))  fHistNegHypPi[trkPtBin]->Fill(logdiff[0]);
          if (IsRapIn(y[AliPID::kKaon]))  fHistNegHypKa[trkPtBin]->Fill(logdiff[1]);
          if (IsRapIn(y[AliPID::kProton]))fHistNegHypPr[trkPtBin]->Fill(logdiff[2]);
          if (fIsMC) {
            if (IsRapIn(y[AliPID::kPion])) {
              if ((lAbsPdgCode != 11) && (lAbsPdgCode != 211) && (lAbsPdgCode != 321) && (lAbsPdgCode != 2212))
                fHistMCNegOtherHypPion[trkPtBin]->Fill(logdiff[0]);
              if (lAbsPdgCode == 11)  fHistMCNegElHypPion[trkPtBin]->Fill(logdiff[0]);
              if (lAbsPdgCode == 211) fHistMCNegPiHypPion[trkPtBin]->Fill(logdiff[0]);
              if (lAbsPdgCode == 321) fHistMCNegKaHypPion[trkPtBin]->Fill(logdiff[0]);
              if (lAbsPdgCode == 2212)fHistMCNegPrHypPion[trkPtBin]->Fill(logdiff[0]);
            }
            if (IsRapIn(y[AliPID::kKaon])) {
              if ((lAbsPdgCode != 11) && (lAbsPdgCode != 211) && (lAbsPdgCode != 321) && (lAbsPdgCode != 2212))
                fHistMCNegOtherHypKaon[trkPtBin]->Fill(logdiff[1]);
              if (lAbsPdgCode == 11)  fHistMCNegElHypKaon[trkPtBin]->Fill(logdiff[1]);
              if (lAbsPdgCode == 211) fHistMCNegPiHypKaon[trkPtBin]->Fill(logdiff[1]);
              if (lAbsPdgCode == 321) fHistMCNegKaHypKaon[trkPtBin]->Fill(logdiff[1]);
              if (lAbsPdgCode == 2212)fHistMCNegPrHypKaon[trkPtBin]->Fill(logdiff[1]);
            }
            if (IsRapIn(y[AliPID::kProton])) {
              if ((lAbsPdgCode != 11) && (lAbsPdgCode != 211) && (lAbsPdgCode != 321) && (lAbsPdgCode != 2212))
                fHistMCNegOtherHypProt[trkPtBin]->Fill(logdiff[2]);
              if (lAbsPdgCode == 11)  fHistMCNegElHypProt[trkPtBin]->Fill(logdiff[2]);
              if (lAbsPdgCode == 211) fHistMCNegPiHypProt[trkPtBin]->Fill(logdiff[2]);
              if (lAbsPdgCode == 321) fHistMCNegKaHypProt[trkPtBin]->Fill(logdiff[2]);
              if (lAbsPdgCode == 2212)fHistMCNegPrHypProt[trkPtBin]->Fill(logdiff[2]);
            }
          }
        }
      }//end lIsGooTrack
    }//end else fill Ntuple
  }//end track loop
  delete lAxis;

  // Post output data.
  PostAllData();
  AliDebug(4, "............. end of Exec");
}

//
//
//________________________________________________________________________
void AliAnalysisTaskSEITSsaSpectra::Terminate(Option_t*)
{
  // Merge output
  // Called once at the end of the query

  AliInfo("end of Terminate");
  return;
}

//
//
//________________________________________________________________________
void AliAnalysisTaskSEITSsaSpectra::SetSPDPileupSelection(Int_t cont, Float_t distance)
{
  fPlpType          = kPileupSPD;
  fMinPlpContribSPD = cont;
  fMinPlpZdistSPD   = distance;
};

//________________________________________________________________________
void AliAnalysisTaskSEITSsaSpectra::SetMVPileUpSelection(Int_t cont, Float_t chi2, Float_t wgthZdiff, Bool_t chkDiffBC)
{
  fPlpType                   = kPileupMV;
  fMinPlpContribMV           = cont;
  fMaxPlpChi2MV              = chi2;
  fMinWDistMV                = wgthZdiff;
  fCheckPlpFromDifferentBCMV = chkDiffBC;
}


//
//
//________________________________________________________________________
Bool_t AliAnalysisTaskSEITSsaSpectra::CheckExtraEvtSelStep(EEvtCut_Type& evtSel)
{
	//Check if has SDD info (if requiered)
	if (fChkIsSDDIn && !fIsMC) {
		TString firedTriggerClasses(fESD->GetFiredTriggerClasses());
    if (!(firedTriggerClasses.Contains("ALL") || firedTriggerClasses.Contains("CENT"))) {
			AliDebug(3, "Event dont accepted by AliEventCuts");
      AliDebug(3, "Event Rejected: SDD out trigger cluster");
      PostAllData();
			evtSel = kIsSDDIn;
      return kFALSE;
    }
	}

  if (fDoSPDCvsTCut) {
    AliAnalysisUtils utils;
    if (utils.IsSPDClusterVsTrackletBG(fESD)) {
    	AliDebug(3, "Event with incompatible SPD clusters and tracklet");
      PostAllData();
    	evtSel = kPassSPDclsVsTCut;
      return kFALSE;
    }
  }

	if (fMultEstimator >= 0 && !IsMultSelected()){
		AliDebug(3, "Event doesn't pass multiplicity selection");
    PostAllData();
    evtSel = kPassMultSel;
    return kFALSE;
	}
	else {
		fHistMultBefEvtSel->Fill(fEvtMult);
	}

	return kTRUE;
}

//
//
//________________________________________________________________________
Bool_t AliAnalysisTaskSEITSsaSpectra::IsEventAccepted(EEvtCut_Type& evtSel)
{

  if (fExtEventCuts) {
    if (fEventCuts.AcceptEvent(fESD)) {
			if (!CheckExtraEvtSelStep(evtSel))
				return kFALSE;
			evtSel = kNEvtCuts;
			return kTRUE;
    }
    AliDebug(3, "Event dont accepted by AliEventCuts");

    if (!fEventCuts.PassedCut(AliEventCuts::kDAQincomplete)) {
      AliDebug(3, "Event with incomplete DAQ");
      PostAllData();
    	evtSel = kIsNotIncDAQ;
      return kFALSE;
    }

		if (!fEventCuts.PassedCut(AliEventCuts::kTrigger)) {
    	AliDebug(3, "Event doesn't pass physics evt. sel. for trigger");
    	PostAllData();
  		evtSel = kPassTrig;
    	return kFALSE;
  	}

    if (!fEventCuts.PassedCut(AliEventCuts::kPileUp)) {
      AliDebug(3, "Event with PileUp");
      PostAllData();
	    evtSel = kIsNotPileup;
      return kFALSE;
    }

		if (!fEventCuts.PassedCut(AliEventCuts::kCorrelations)) {
      AliDebug(3, "Event with PileUp");
      PostAllData();
    	evtSel = kCorrelation;
      return kFALSE;
    }

    if (!fEventCuts.PassedCut(AliEventCuts::kVertex) || !fEventCuts.PassedCut(AliEventCuts::kVertexQuality)) {
      AliDebug(3, "Event doesn't pass has good vtx sel");
      PostAllData();
      evtSel = kHasRecVtx;
      return kFALSE;
    }

    if (!fEventCuts.PassedCut(AliEventCuts::kVertexPosition)) {
      AliDebugF(3, "Vertex with Z>%f cm", fMaxVtxZCut);
      PostAllData();
    	evtSel = kHasGoodVtxZ;
      return kFALSE;
    }
		evtSel = kNEvtCuts;
		return kFALSE;
  } else {
		if (fRejIncDAQ && fESD->IsIncompleteDAQ()) {
      AliDebug(3, "Event with incomplete DAQ");
      PostAllData();
    	evtSel = kIsNotIncDAQ;
      return kFALSE;
    }

		UInt_t maskPhysSel = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  maskPhysSel &= fTriggerSel;
  	if (!maskPhysSel) {
    	AliDebugF(3, "Event doesn't pass physics evt. sel. for trigger %d", fTriggerSel);
    	PostAllData();
	  	evtSel = kPassTrig;
    	return kFALSE;
  	}

    if (fDoSPDCvsTCut) {
      AliAnalysisUtils utils;
      if (utils.IsSPDClusterVsTrackletBG(fESD)) {
        AliDebug(3, "Event with incompatible SPD clusters and tracklet");
        PostAllData();
    		evtSel = kPassSPDclsVsTCut;
        return kFALSE;
      }
    }

    if (IsPileUp()) {
      AliDebug(3, "Event with PileUp");
      PostAllData();
    	evtSel = kIsNotPileup;
      return kFALSE;
    }

    if (!IsVtxReconstructed()) {
      AliDebug(3, "Event doesn't pass vtx quality sel");
      PostAllData();
    	evtSel = kHasRecVtx;
      return kFALSE;
    }

    if (!IsGoodVtxZ()) {
      AliDebugF(3, "Vertex with Z>%f cm", fMaxVtxZCut);
      PostAllData();
    	evtSel = kHasGoodVtxZ;
      return kFALSE;
    }

		if (!CheckExtraEvtSelStep(evtSel))
			return kFALSE;
	}

	evtSel = kNEvtCuts;
  return kTRUE;
}

//
//
//________________________________________________________________________
void AliAnalysisTaskSEITSsaSpectra::SetupStandardEventCutsForRun1()
{
  fChkIsSDDIn   = kTRUE;
  fMultMethod   = 0;
  fExtEventCuts = kFALSE;
  fRejIncDAQ    = kFALSE;
  fDoSPDCvsTCut = kFALSE;
  fPlpType      = kNoPileup;
  fTriggerSel   = AliVEvent::kMB;
  fMaxVtxZCut   = 10.;
  fReqBothVtx   = kFALSE;
  fChkVtxSPDRes = kTRUE;
  fChkVtxZSep   = kFALSE;
}

//
//
//________________________________________________________________________
void AliAnalysisTaskSEITSsaSpectra::SetupEventCutsForRun1pPb()
{
  fChkIsSDDIn   = kTRUE;
  fMultMethod   = 2;
  fExtEventCuts = kFALSE;
  fRejIncDAQ    = kFALSE;
  fDoSPDCvsTCut = kFALSE;
  fPlpType      = kNoPileup;
  fTriggerSel   = AliVEvent::kINT7;
  fMaxVtxZCut   = 10.;
  fReqBothVtx   = kFALSE;
  fChkVtxSPDRes = kTRUE;
  fChkVtxZSep   = kFALSE;
}

//
//
//________________________________________________________________________
void AliAnalysisTaskSEITSsaSpectra::SetupStandardEventCutsForRun2()
{
  fChkIsSDDIn   = kTRUE;
  fMultMethod   = 0;
  fExtEventCuts = kFALSE;
  fRejIncDAQ    = kTRUE;
  fDoSPDCvsTCut = kTRUE;
  fPlpType      = kPileupSPD;
  fMinPlpContribSPD = 3;
  fMinPlpZdistSPD   = .8;
  fTriggerSel   = AliVEvent::kINT7;
  fMaxVtxZCut   = 10.;
  fReqBothVtx   = kFALSE;
  fChkVtxSPDRes = kTRUE;
  fChkVtxZSep   = kTRUE;
}

//
//
//________________________________________________________________________
Bool_t AliAnalysisTaskSEITSsaSpectra::IsMultSelected()
{
	if (fMultMethod > 5) {
		AliWarning(". Skipping multiplicity selection");
		fMultMethod = 0;
	}
  if (!fMultMethod) return kTRUE; 		// skip multiplicity check

  if (fMultMethod == 1) { //New multiplicity/centrality class framework
  	AliMultSelection* fMultSel = (AliMultSelection*) fESD->FindListObject("MultSelection");
    if (!fMultSel) {
      //If you get this warning (and lPercentiles 300) please check that the AliMultSelectionTask actually ran (before your task)
      AliWarning("AliMultSelection object not found!");
      return kFALSE;
    } else {
      //Event selection is embedded in the Multiplicity estimator so that the Multiplicity percentiles are well defined and refer to the same sample
      fEvtMult = fMultSel->GetMultiplicityPercentile(fMultEstimator.Data(), kFALSE);
    }
	} else if (fMultMethod == 2) { //OLD multiplicity/centrality class framework
    AliCentrality* centrality = fESD->GetCentrality();
    fEvtMult = centrality->GetCentralityPercentile(fMultEstimator.Data());
 	} else if (fMultMethod == 3){ //selection on the event multiplicity based on global tracks
		// tracks+tracklets
    fEvtMult = (float)AliESDtrackCuts::GetReferenceMultiplicity(fESD, AliESDtrackCuts::kTrackletsITSTPC, 0.8);
  } else if (fMultMethod == 4) {
    // tracklets
    fEvtMult = (float)AliESDtrackCuts::GetReferenceMultiplicity(fESD, AliESDtrackCuts::kTracklets, 0.8);
  } else if (fMultMethod == 5) {
    // clusters in SPD1
    const AliMultiplicity* mult = fESD->GetMultiplicity();
    Float_t nClu1 = (Float_t)mult->GetNumberOfITSClusters(1);
    fEvtMult = AliESDUtils::GetCorrSPD2(nClu1, fESD->GetPrimaryVertexSPD()->GetZ()) + 0.5;
  }

	if (fEvtMult < fLowMult || fEvtMult >= fUpMult)
		return kFALSE;

  return kTRUE;
}

//
//
//________________________________________________________________________
Bool_t AliAnalysisTaskSEITSsaSpectra::IsPileUp()
{
  AliAnalysisUtils utils;

  switch (fPlpType) {

    case kPileupSPD :
      utils.SetMinPlpContribSPD(fMinPlpContribSPD);
      utils.SetMinPlpZdistSPD(fMinPlpZdistSPD);
      utils.SetUseMVPlpSelection(kFALSE);
      break;

    case kPileupMV :
      utils.SetMinPlpContribMV(fMinPlpContribMV);
      utils.SetMaxPlpChi2MV(fMaxPlpChi2MV);
      utils.SetMinWDistMV(fMinWDistMV);
      utils.SetCheckPlpFromDifferentBCMV(fCheckPlpFromDifferentBCMV);
      utils.SetUseMVPlpSelection(kTRUE);
      break;

    default :
      return kFALSE;
      break;
  }

  return utils.IsPileUpEvent(fESD);
}

//
//
//________________________________________________________________________
Bool_t AliAnalysisTaskSEITSsaSpectra::IsGoodVtxZ()
{
  const AliESDVertex* vtx = fESD->GetPrimaryVertex();
  fHistVtxZ->Fill(vtx->GetZ());
  if (TMath::Abs(vtx->GetZ()) > fMaxVtxZCut) return kFALSE;

  return kTRUE;
}

//
//
//________________________________________________________________________
Bool_t AliAnalysisTaskSEITSsaSpectra::IsGoodSPDvtxRes(const AliESDVertex* spdVtx)
{
  if (!spdVtx)
    return kFALSE;

  Double_t cov[6] = {0};
  spdVtx->GetCovarianceMatrix(cov);
  Double_t zRes = TMath::Sqrt(cov[5]);
  if (spdVtx->IsFromVertexerZ() && (spdVtx->GetDispersion() > 0.04 || zRes > 0.25))
    return kFALSE;

  return kTRUE;
}

//
//
//________________________________________________________________________
Bool_t AliAnalysisTaskSEITSsaSpectra::IsVtxReconstructed()
{
  const AliESDVertex* trkVertex = fESD->GetPrimaryVertexTracks();
  const AliESDVertex* spdVertex = fESD->GetPrimaryVertexSPD();
  Bool_t hasTrk = trkVertex->GetStatus();
  Bool_t hasSPD = spdVertex->GetStatus();

  //Note that AliVertex::GetStatus checks that N_contributors is > 0
  //reject events if both are explicitly requested and none is available
  if (fReqBothVtx && !(hasSPD && hasTrk))
    return kFALSE;

  //reject events if none between the SPD or track verteces are available
  //if no trk vertex, try to fall back to SPD vertex;
  if (!hasTrk) {
    if (!hasSPD)
      return kFALSE;
    //on demand check the spd vertex resolution and reject if not satisfied
    if (fChkVtxSPDRes && !IsGoodSPDvtxRes(spdVertex))
      return kFALSE;
  } else {
    if (hasSPD) {
      //if enabled check the spd vertex resolution and reject if not satisfied
      //if enabled, check the proximity between the spd vertex and trak vertex, and reject if not satisfied
      if (fChkVtxSPDRes && !IsGoodSPDvtxRes(spdVertex))
        return kFALSE;
      if ((fChkVtxZSep && TMath::Abs(spdVertex->GetZ() - trkVertex->GetZ()) > 0.5))
        return kFALSE;
    }
  }

  return kTRUE;
}

//
//
//________________________________________________________________________
Bool_t AliAnalysisTaskSEITSsaSpectra::DCAcut(Double_t impactXY, Double_t impactZ, Double_t pt) const
{
  return (DCAcutXY(impactXY, pt) && DCAcutZ(impactZ, pt));
}

//
//
//________________________________________________________________________
Bool_t AliAnalysisTaskSEITSsaSpectra::DCAcutXY(Double_t impactXY, Double_t pt) const
{
  // cut on transverse impact parameter updated on 20-5-2010
  // from the study of L. Milano, F. Prino on the ITS standalone tracks
  // using the common binning of the TPC tracks
  Double_t xyMax = fDCAxyCutFunc->Eval(pt); //in micron
  AliDebugF(3, "Max value for the DCAxy Cut is:%f Measured value for DCAxy is:%f cut to %.0f sigmas\n",
            xyMax,
            TMath::Abs(impactXY) * 10000,
            fDCAxyCutFunc->GetParameter(3));
  if ((TMath::Abs(impactXY) * 10000) > xyMax) return kFALSE;

  return kTRUE;
}

//
//
//________________________________________________________________________
Bool_t AliAnalysisTaskSEITSsaSpectra::DCAcutZ(Double_t impactZ, Double_t pt) const
{
  // cut on transverse impact parameter updated on 20-5-2010
  // from the study of L. Milano, F. Prino on the ITS standalone tracks
  // using the common binning of the TPC tracks

  Double_t zMax = fDCAzCutFunc->Eval(pt); //in micron
  AliDebugF(3, "Max value for the DCAz Cut is:%f Measured value for DCAz is:%f cut to %.0f sigmas\n",
            zMax,
            TMath::Abs(impactZ) * 10000,
            fDCAzCutFunc->GetParameter(3));
  if ((TMath::Abs(impactZ) * 10000) > zMax) return kFALSE;

  return kTRUE;
}

//
//
//________________________________________________________________________
Double_t AliAnalysisTaskSEITSsaSpectra::CookdEdx(Double_t* s) const
{
  // truncated mean for the dEdx
  Int_t nc = 0;
  Double_t dedx[4] = {0., 0., 0., 0.};
  for (Int_t il = 0; il < 4; il++) { // count good (>0) dE/dx values
    if (s[il] > fMindEdx) {
      dedx[nc] = s[il];
      nc++;
    }
  }
  if (nc < fMinNdEdxSamples) return -1.;

  Double_t tmp;
  Int_t swap; // sort in ascending order
  do {
    swap = 0;
    for (Int_t i = 0; i < nc - 1; i++) {
      if (dedx[i] <= dedx[i + 1]) continue;
      tmp = dedx[i];
      dedx[i] = dedx[i + 1];
      dedx[i + 1] = tmp;
      swap++;
    }
  } while (swap);

  Double_t sumamp = 0, sumweight = 0;
  Double_t weight[4] = {1., 1., 0., 0.};
  if (nc == 3) weight[1] = 0.5;
  else if (nc < 3) weight[1] = 0.;
  for (Int_t i = 0; i < nc; i++) {
    sumamp += dedx[i] * weight[i];
    sumweight += weight[i];
  }
  return sumamp / sumweight;
}

//
//
//________________________________________________________________________
Double_t AliAnalysisTaskSEITSsaSpectra::Eta2y(Double_t pt, Double_t m, Double_t eta) const
{
  // convert eta to y
  Double_t mt = TMath::Sqrt(m * m + pt * pt);
  return TMath::ASinH(pt / mt * TMath::SinH(eta));
}

//
//
//________________________________________________________________________
void AliAnalysisTaskSEITSsaSpectra::AnalyseMCParticles(AliMCEvent* lMCevent,
                                                       EEvtCut_Type lastEvtCutPassed,
                                                       bool lHasGoodVtxGen)
{
  //first loop on stack, and get primary pi,k,p spectra with truth MC values
  Int_t nTrackMC = (lMCevent) ? lMCevent->GetNumberOfTracks() : 0;
  for (Int_t i_mcTrk = 0; i_mcTrk < nTrackMC; ++i_mcTrk) {

    TParticle* mcTrk = ((AliMCParticle*)lMCevent->GetTrack(i_mcTrk))->Particle();
    Int_t pdg  = mcTrk->GetPdgCode();

    Int_t iChg = pdg >= 0 ? 0 : 1; // only works for charged pi,K,p (0 Pos, 1 Neg)

    Int_t iPart = -1;
    if      (TMath::Abs(pdg) == 211)  iPart = 0; // select Pi+/Pi- only
    else if (TMath::Abs(pdg) == 321)  iPart = 1; // select K+/K- only
    else if (TMath::Abs(pdg) == 2212) iPart = 2; // select p+/p- only
    else  continue;

    Double_t mcPt  = mcTrk->Pt();
    if (mcPt < 0.08 || mcPt > 1.0) continue; // pt cut

    Double_t mcEta = mcTrk->Eta();
    Double_t mcRap = Eta2y(mcPt, mcTrk->GetMass(), mcEta) + fCMSRapFct;

    Bool_t lIsPhysPrimary = lMCevent->IsPhysicalPrimary(i_mcTrk);
    if (fFillNtuple) {
      //filling MC ntuple
      Float_t xntMC[8];
      Int_t indexMC = 0;
      xntMC[indexMC++] = (Float_t)mcPt;
      xntMC[indexMC++] = (Float_t)pdg;
      xntMC[indexMC++] = (Float_t)iChg;
      xntMC[indexMC++] = (Float_t)mcEta;
      xntMC[indexMC++] = (Float_t)mcRap;
      xntMC[indexMC++] = (Float_t)lIsPhysPrimary;
      xntMC[indexMC++] = (Float_t)lastEvtCutPassed;
      xntMC[indexMC++] = (Float_t)fESD->GetRunNumber();

      fNtupleMC->Fill(xntMC);
    }

    if (!IsRapIn(mcRap)) continue;  //rapidity cut
    if (!lIsPhysPrimary) continue;  //primary cut

    int index = iPart * kNchg + iChg;
    for (int istep = (int)kIsReadable; istep <= (int)lastEvtCutPassed; ++istep) {
      fHistPrimMCGenVtxZall[index]->Fill(istep, TMath::Abs(mcPt));
      if (lHasGoodVtxGen)
        fHistPrimMCGenVtxZcut[index]->Fill(istep, TMath::Abs(mcPt));
    }
  }//end MC stack loop
}

//
//
//________________________________________________________________________
void AliAnalysisTaskSEITSsaSpectra::PostAllData()
{
  PostData(1, fOutput);
  PostData(2, fListCuts);
  if (fFillNtuple) PostData(3, fListTree);

  return;
}

//
//
//________________________________________________________________________
Int_t AliAnalysisTaskSEITSsaSpectra::GetTrackPid(AliESDtrack* track, Double_t* logdiff) const
{
  AliPID::EParticleType iType[4] = {AliPID::kPion, AliPID::kKaon, AliPID::kProton, AliPID::kDeuteron};

  Int_t pid = -1;

  Double_t dEdxLay[4];
  track->GetITSdEdxSamples(dEdxLay);
  Double_t dedx = track->GetITSsignal();
  Float_t     p = track->GetP();
  if (fIsMC && fSmearMC) {
    dedx = fRandGener->Gaus(dedx, fSmeardEdx * dedx);
    p = fRandGener->Gaus(p, fSmearP * p);
  }

  Double_t bbtheo[4];
  for (Int_t i = 0; i < 4; i++) {
    Float_t mass = AliPID::ParticleMass(iType[i]);
    bbtheo[i]  = fITSPIDResponse->BetheITSsaHybrid(p, mass);
    logdiff[i] = TMath::Log(dedx) - TMath::Log(bbtheo[i]);
  }

  UInt_t clumap = track->GetITSClusterMap();
  Int_t nPtsForPid = 0;
  for (Int_t j = 2; j < 6; j++)
    if (TESTBIT(clumap, j)) nPtsForPid++;

  Float_t resodedx = fITSPIDResponse->GetResolution(1, nPtsForPid, kTRUE);
  switch (fPidMethod) {
    case kNSigCut : {
      Double_t nSigmaMin = 99999;
      for (Int_t ispc = 0; ispc < kNspc; ispc++) {
        Double_t bb = bbtheo[ispc];
        Double_t nSigma = TMath::Abs((dedx - bb) / (resodedx * bb));
        if (nSigma < nSigmaMin) {
          nSigmaMin = nSigma;
          pid = iType[ispc];
        }
      }
      pid = (nSigmaMin < fMinNSigma) ? pid : -1;
      break;
    }
    case kMeanCut : {
      for (Int_t ispc = 0; ispc < kNspc; ispc++) {
        if (dedx < bbtheo[0]) {
          Double_t nsigma = TMath::Abs((dedx - bbtheo[0]) / (resodedx * bbtheo[0]));
          pid = (nsigma < 2.) ? AliPID::kPion : -1;
          break;
        }
        Double_t bb = bbtheo[ispc];
        Int_t edge = (dedx < bb) ? ispc - 1 : ispc + 1;
        Double_t bbdistance = TMath::Abs((bbtheo[ispc] - bbtheo[edge]) / 2);

        if (bbdistance != 0.f) {
          Double_t nsigma = TMath::Abs((dedx - bb) / bbdistance);
          if (nsigma < 1.) pid = iType[ispc];
        } else {
          AliInfo("Zero distance");
          pid = AliPID::kPion;
        }
      }
      break;
    }
    case kLanGaus : {
      Double_t prior[AliPID::kSPECIES];
      Double_t probITS[AliPID::kSPECIES];
      GetPriors(track, prior);
      if (!fITSPidParams)
        fITSPIDResponse->GetITSProbabilities(track->GetP(), dEdxLay, probITS, fIsMC);
      else
        fITSPIDResponse->GetITSProbabilities(track->GetP(), dEdxLay, probITS, fITSPidParams);
      pid = GetMostProbable(probITS, prior);
      break;
    }
    default:
      AliInfo("Pid method not known");
      return -1;
      break;
  }

  //Sigma Separation
  for (Int_t ispc = 0; ispc < kNspc; ++ispc) {
    int ichg  = (track->GetSign() > 0) ? 0 : 1;
    int index = ispc * kNchg + ichg;
    Double_t bb = bbtheo[ispc];
    fHistNSigmaSep[index]->Fill(track->Pt(), ((dedx - bb) / (resodedx * bb)));
  }
  //fill charge distribution histo to check the calibration
  for (Int_t j = 0; j < 4; j++) {
    if (dEdxLay[j] < fMindEdx) continue;
    fHistCharge[j]->Fill(dEdxLay[j]);
  }

  return (pid == -1) ? 0 : pid;
}

//
//
//________________________________________________________________________
Int_t AliAnalysisTaskSEITSsaSpectra::GetMostProbable(const Double_t* pDens, const Double_t* prior) const
{
  // get the most probable particle id hypothesis
  // assuming the a priori probabilities "prior"
  Double_t max = 0.;
  Int_t id = AliPID::kPion;
  for (Int_t i = 0; i < AliPID::kSPECIES; ++i) {
    Double_t prob = pDens[i] * prior[i];
    if (prob > max) {
      max = prob;
      id = i;
    }
  }
  if (max == 0)
    AliError("Invalid probability densities or priors");

  return id;
}

//________________________________________________________________________
void AliAnalysisTaskSEITSsaSpectra::GetPriors(const AliVTrack* track, Double_t* prior) const
{
  const char* prtName[AliPID::kSPECIES] = {"El", "Mu", "Pi", "Ka", "Pr"};
  Double_t pt = track->Pt();

  if      (pt < 0.08) pt = 0.0801;
  else if (pt > 0.9999) pt = 0.9999;

  Float_t usedCent = fEvtMult;//+ 5*(cent>0);
  if (usedCent < -.99) usedCent = -0.9;
  if (usedCent > 99.9) usedCent = 99.9;

  for (Int_t i = 0; i < AliPID::kSPECIES; ++i)
    prior[i] = 1. / AliPID::kSPECIES;

  if (!fUseDefaultPriors) {
    if (!fListPriors) {
      AliInfo("Errors use priors set but no priors list found.");
      AliInfo("Default(equal) priors used.");
      return;
    } else {
      for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
        TH2F* hPriors = dynamic_cast<TH2F*>(fListPriors->FindObject(Form("hDefaultITSsa%sPriors", prtName[i])));
        if (!hPriors) {
          AliInfo("Errors use priors set but no histogram with priors found.");
          AliInfo("Default(equal) priors used.");
          return;
        }
        prior[i] = hPriors->Interpolate(usedCent, pt);
      }
    }
  } else {
    if (pt > .160) prior[0] = 0;
    prior[1] = 0.;
  }

  Double_t sumP = 0;
  for (Int_t i = 0; i < AliPID::kSPECIES; i++) sumP += prior[i];
  for (Int_t i = 0; i < AliPID::kSPECIES; i++) prior[i] /= sumP;

  return;
}

//
//
//________________________________________________________________________
void AliAnalysisTaskSEITSsaSpectra::ComputeBayesProbabilities(Double_t* probs, const Double_t* pDens, const Double_t* prior)
{
  //Calculate bayesian probabilities
  Double_t sum = 0.;
  for ( Int_t i = 0; i < AliPID::kSPECIES; i++) sum += pDens[i] * prior[i];

  if (sum <= 0) {
    AliError("Invalid probability densities or priors");
    for (Int_t i = 0; i < AliPID::kSPECIES; i++)
      probs[i] = -1;
    return;
  }

  for (Int_t i = 0; i < AliPID::kSPECIES; i++)
    probs[i] = pDens[i] * prior[i] / sum;
}
