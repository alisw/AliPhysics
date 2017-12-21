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
#include <TH3F.h>
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
  fCentBins(),
  fDCABins(),
  fPtBins(),
  fTriggerSel(AliVEvent::kMB),
  fMaxVtxZCut(10.),
	fChkIsEventINELgtZERO(kFALSE),
  fChkIsSDDIn(kTRUE),
  fRejIncDAQ(kTRUE),
  fDoSPDCvsTCut(kTRUE),
	fUseSelectVertex2015pp(kTRUE),
  fChkVtxSPDRes(kTRUE),
  fChkVtxZSep(kFALSE),
  fReqBothVtx(kFALSE),
  fExtEventCuts(kFALSE),
  fMultMethod(0),
  fMultEstimator("V0M"),
	fMultEvSel(kFALSE),
  fLowMult(-5.),
  fUpMult(500),
  fEvtMult(-1),
  fPlpType(BIT(kNoPileup)),
  fMinPlpContribSPD(5),
  fMinPlpZdistSPD(.8),
  fMinPlpContribMV(5),
  fMaxPlpChi2MV(5.),
  fMinWDistMV(15.),
  fCheckPlpFromDifferentBCMV(kFALSE),
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
  fYear(2010),
  fPidMethod(kMeanCut),
  fUseDefaultPriors(kTRUE),
  fIsMC(kFALSE),
  fFillNtuple(kFALSE),
  fFillIntDistHist(kFALSE),
  fRandGener(0x0),
  fSmearMC(kFALSE),
  fSmearP(0.),
  fSmeardEdx(0.)
{
  //Constructor
  fRandGener = new TRandom3(0);

  // Inizialize bins array
  fCentBins.Set(3);
  fCentBins[0]=-5.f;
  fCentBins[1]=0.f;
  fCentBins[2]=100.f;

  float ptBins[kNbins + 1] = {
    0.08, 0.10, 0.12, 0.14, 0.16, 0.18, 0.20, 0.25,
    0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65,
    0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.0
  };
  fPtBins.Set(kNbins+1, ptBins);

  SetDCABins(2000, -2, 2);

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

/*      for (Int_t iptbin = 0; iptbin < kNbins; ++iptbin) {
        fHistRecoDCA[index][iptbin] = NULL; //! histo with DCA distibution

        fHistPrimDCA[index][iptbin] = NULL; //! histo with DCA distibution and dedx PID
        fHistSstrDCA[index][iptbin] = NULL; //! histo with DCA distibution and dedx PID
        fHistSmatDCA[index][iptbin] = NULL; //! histo with DCA distibution and dedx PID

        fHistMCtruthPrimDCA[index][iptbin] = NULL; //! histo with DCA distibution, MC truth
        fHistMCtruthSstrDCA[index][iptbin] = NULL; //! histo with DCA distibution, MC truth
        fHistMCtruthSmatDCA[index][iptbin] = NULL; //! histo with DCA distibution, MC truth
      } */
    }
  }
  for (Int_t iL = 0; iL < 4; ++iL) fHistCharge[iL] = NULL;

  //dEdx distributions
  fHistPosHypPi = NULL;
  fHistPosHypKa = NULL;
  fHistPosHypPr = NULL;
  fHistNegHypPi = NULL;
  fHistNegHypKa = NULL;
  fHistNegHypPr = NULL;

  //dEdx distributions for MC
  fHistMCPosOtherHypPion = NULL;
  fHistMCPosOtherHypKaon = NULL;
  fHistMCPosOtherHypProt = NULL;
  fHistMCPosElHypPion    = NULL;
  fHistMCPosElHypKaon    = NULL;
  fHistMCPosElHypProt    = NULL;
  fHistMCPosPiHypPion    = NULL;
  fHistMCPosPiHypKaon    = NULL;
  fHistMCPosPiHypProt    = NULL;
  fHistMCPosKaHypPion    = NULL;
  fHistMCPosKaHypKaon    = NULL;
  fHistMCPosKaHypProt    = NULL;
  fHistMCPosPrHypPion    = NULL;
  fHistMCPosPrHypKaon    = NULL;
  fHistMCPosPrHypProt    = NULL;

  fHistMCNegOtherHypPion = NULL;
  fHistMCNegOtherHypKaon = NULL;
  fHistMCNegOtherHypProt = NULL;
  fHistMCNegElHypPion    = NULL;
  fHistMCNegElHypKaon    = NULL;
  fHistMCNegElHypProt    = NULL;
  fHistMCNegPiHypPion    = NULL;
  fHistMCNegPiHypKaon    = NULL;
  fHistMCNegPiHypProt    = NULL;
  fHistMCNegKaHypPion    = NULL;
  fHistMCNegKaHypKaon    = NULL;
  fHistMCNegKaHypProt    = NULL;
  fHistMCNegPrHypPion    = NULL;
  fHistMCNegPrHypKaon    = NULL;
  fHistMCNegPrHypProt    = NULL;

  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
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
  
  const int nPtBins = fPtBins.GetSize() - 1;
  const int nCentBins = fCentBins.GetSize() - 1;
  const int nDCABins = fDCABins.GetSize() - 1;
  const float* ptBins = fPtBins.GetArray();
  const float* centBins = fCentBins.GetArray();
  const float* dcaBins = fDCABins.GetArray();

  float evBins[kNEvtCuts + 1];
  SetBins(kNEvtCuts,.5f,kNEvtCuts+.5f,evBins);

	const char* notApp= "_notApplied";
  fHistNEvents = new TH2I("fHistNEvents", "Number of processed events;Centrality (%);",
      nCentBins,centBins,kNEvtCuts,evBins);
  fHistNEvents->Sumw2();
  fHistNEvents->SetMinimum(0);
  fHistNEvents->GetYaxis()->SetBinLabel(kIsReadable,  "Readable");
  fHistNEvents->GetYaxis()->SetBinLabel(kPassMultSel, Form("PassMultSel%s",(fMultMethod ? "" : notApp)));
  fHistNEvents->GetYaxis()->SetBinLabel(kIsSDDIn,     Form("HasSDDIn%s",   (fChkIsSDDIn ? "" : notApp)));
  fHistNEvents->GetYaxis()->SetBinLabel(kIsNotIncDAQ, Form("PassIncDAQ%s", (fRejIncDAQ  ? "" : notApp)));
  fHistNEvents->GetYaxis()->SetBinLabel(kPassTrig, "PassPhysSelTrig");
  fHistNEvents->GetYaxis()->SetBinLabel(kPassINELgtZERO, Form("PassINELgtZERO%s", (fChkIsEventINELgtZERO ? "" : notApp)));
  fHistNEvents->GetYaxis()->SetBinLabel(kCorrelations, Form("Correlations%s", (fExtEventCuts ? "" : notApp)));
  fHistNEvents->GetYaxis()->SetBinLabel(kPassSPDclsVsTCut, Form("PassClsVsTrackletBG%s", (fDoSPDCvsTCut ? "" : notApp)));
  fHistNEvents->GetYaxis()->SetBinLabel(kIsPileupSPD,           Form("PassIsPileupSPD%s", ((fPlpType & BIT(kPileupSPD) ? "" : notApp))));
  fHistNEvents->GetYaxis()->SetBinLabel(kIsPileupSPDinMultBins, Form("PassIsPileupSPDinMultBins%s", ((fPlpType & BIT(kPileupInMultBins) ? "" : notApp))));
  fHistNEvents->GetYaxis()->SetBinLabel(kIsPileupMV,            Form("PassIsPileupMV%s", ((fPlpType & BIT(kPileupMV) ? "" : notApp))));
  fHistNEvents->GetYaxis()->SetBinLabel(kHasRecVtx,   "HasVertex");
  fHistNEvents->GetYaxis()->SetBinLabel(kHasGoodVtxZ, "HasGoodVertex");
  fHistNEvents->GetYaxis()->SetBinLabel(kNEvtCuts,    "IsSelected");
  fOutput->Add(fHistNEvents);

  fHistMCEvents = new TH2I("fHistMCEvents", "Number of processed events;Centrality (%);",
      nCentBins,centBins,kNEvtCuts,evBins);
  fHistMCEvents->Sumw2();
  fHistMCEvents->SetMinimum(0);
	fHistMCEvents->GetYaxis()->SetBinLabel(kIsReadable,  "Readable");
  fHistMCEvents->GetYaxis()->SetBinLabel(kPassMultSel, Form("PassMultSel%s",(fMultMethod ? "" : notApp)));
  fHistMCEvents->GetYaxis()->SetBinLabel(kIsSDDIn,     Form("HasSDDIn%s",   (fChkIsSDDIn ? "" : notApp)));
  fHistMCEvents->GetYaxis()->SetBinLabel(kIsNotIncDAQ, Form("PassIncDAQ%s", (fRejIncDAQ  ? "" : notApp)));
  fHistMCEvents->GetYaxis()->SetBinLabel(kPassTrig, "PassPhysSelTrig");
  fHistMCEvents->GetYaxis()->SetBinLabel(kPassINELgtZERO, Form("PassINELgtZERO%s", (fChkIsEventINELgtZERO ? "" : notApp)));
  fHistMCEvents->GetYaxis()->SetBinLabel(kCorrelations, Form("Correlations%s", (fExtEventCuts ? "" : notApp)));
  fHistMCEvents->GetYaxis()->SetBinLabel(kPassSPDclsVsTCut, Form("PassClsVsTrackletBG%s", (fDoSPDCvsTCut ? "" : notApp)));
  fHistMCEvents->GetYaxis()->SetBinLabel(kIsPileupSPD,           Form("PassIsPileupSPD%s", ((fPlpType & BIT(kPileupSPD) ? "" : notApp))));
  fHistMCEvents->GetYaxis()->SetBinLabel(kIsPileupSPDinMultBins, Form("PassIsPileupSPDinMultBins%s", ((fPlpType & BIT(kPileupInMultBins) ? "" : notApp))));
  fHistMCEvents->GetYaxis()->SetBinLabel(kIsPileupMV,            Form("PassIsPileupMV%s", ((fPlpType & BIT(kPileupMV) ? "" : notApp))));
  fHistMCEvents->GetYaxis()->SetBinLabel(kHasRecVtx,   "HasVertex");
  fHistMCEvents->GetYaxis()->SetBinLabel(kHasGoodVtxZ, "HasGoodVertex");
  fHistMCEvents->GetYaxis()->SetBinLabel(kNEvtCuts,    "IsSelected");
  fOutput->Add(fHistMCEvents);

  Int_t nMultBin = 100;
  Float_t lMultBinLimit[nMultBin+1];
  SetBins(nMultBin,0,100,lMultBinLimit);

  fHistMultBefEvtSel = new TH1F("fHistMultBefEvtSel", "Event Multiplicity before event selection;Centrality (%)",
      nMultBin, lMultBinLimit);
  fHistMultBefEvtSel->Sumw2();
  fHistMultBefEvtSel->SetMinimum(0);
  fOutput->Add(fHistMultBefEvtSel);

  fHistMultAftEvtSel = new TH1F("fHistMultAftEvtSel", "Event Multiplicity after event selection;Centrality (%)",
      nMultBin, lMultBinLimit);
  fHistMultAftEvtSel->Sumw2();
  fHistMultAftEvtSel->SetMinimum(0);
  fOutput->Add(fHistMultAftEvtSel);
  const int nVtxBins = 400;
  float vtxBins[nVtxBins + 1];
  SetBins(nVtxBins,-20,20,vtxBins);
  fHistVtxZ = new TH2F("fHistVtxZ", "Vtx Z distribution;Centrality (%);Z_vtx",
      nCentBins,centBins,nVtxBins,vtxBins);
  fHistVtxZ->Sumw2();
  fHistVtxZ->SetMinimum(0);
  fOutput->Add(fHistVtxZ);

  const char spcName[kNspc][3] = {"Pi", "Ka", "Pr"};
  const char chgName[kNchg][4] = {"Pos", "Neg"};

  char* histName;
  const int nTrkBins = 20;
  float trkBins[nTrkBins + 1];
  SetBins(nTrkBins,.5f,nTrkBins+.5f,trkBins);
  //Histo with track cuts
  for (Int_t iChg = 0; iChg < kNchg; ++iChg) {
    histName = Form("fHistNTracks%s", chgName[iChg]);
    fHistNTracks[iChg] = new TH3F(histName, "Number of ITSsa tracks;Centrality (%);#it{p}_{T} (GeV/#it{c}); Trk Selection",
        nCentBins,centBins,nPtBins,ptBins,nTrkBins,trkBins);
    fHistNTracks[iChg]->Sumw2();

    TString label("no selection");//1
    fHistNTracks[iChg]->GetZaxis()->SetBinLabel(kHasNoSelection, label.Data());
    label = "ITSsa";//2
    fHistNTracks[iChg]->GetZaxis()->SetBinLabel(kIsITSsa, label.Data());
    label = "ITSrefit";//3
    fHistNTracks[iChg]->GetZaxis()->SetBinLabel(kIsITSrefit, label.Data());
    label = "neutral particle";//4
    fHistNTracks[iChg]->GetZaxis()->SetBinLabel(kIsNotNeutralParticle, label.Data());
    label = "SPDcls";//7
    fHistNTracks[iChg]->GetZaxis()->SetBinLabel(kPassSPD, label.Data());
    label = "SDD+SSD cls";//8
    fHistNTracks[iChg]->GetZaxis()->SetBinLabel(kPassPIDcls, label.Data());
    label = "chi2/ncls";//9
    fHistNTracks[iChg]->GetZaxis()->SetBinLabel(kPassChi2Ncls, label.Data());
    label = "eta";//11
    fHistNTracks[iChg]->GetZaxis()->SetBinLabel(kIsInEta, label.Data());
    label = "dE/dx < 0";//12
    fHistNTracks[iChg]->GetZaxis()->SetBinLabel(kPassdEdx, label.Data());
    label = "Pt cut";
    fHistNTracks[iChg]->GetZaxis()->SetBinLabel(kPassPtCut, label.Data());
    label = "DCAz";//13
    fHistNTracks[iChg]->GetZaxis()->SetBinLabel(kPassDCAzcut, label.Data());
    label = "DCAxy";//14
    fHistNTracks[iChg]->GetZaxis()->SetBinLabel(kPassDCAxycut, label.Data());
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
      fHistReco[index] = new TH2F(histName,";Centrality (%);#it{p}_{T} (GeV/#it{c});",
          nCentBins,centBins,nPtBins,ptBins);
      fOutput->Add(fHistReco[index]);

      histName = Form("fHistRecoDCA%s%s", spcName[ispc], chgName[ichg]);
      fHistRecoDCA[index] = new TH3F(histName,";Centrality (%);#it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)",
          nCentBins,centBins,nPtBins,ptBins,nDCABins,dcaBins);
      fOutput->Add(fHistRecoDCA[index]);

      if (fIsMC) {
        //
        //Histograms MC part Gen bef and afte all selection Good Vertex Gen.
        histName = Form("fHistPrimMCGenVtxZall%s%s", spcName[ispc], chgName[ichg]);
        fHistPrimMCGenVtxZall[index] = new TH3F(histName,";Centrality (%);#it{p}_{T} (GeV/#it{c});",
          nCentBins,centBins,nPtBins,ptBins,kNEvtCuts,evBins);

        histName = Form("fHistPrimMCGenVtxZcut%s%s", spcName[ispc], chgName[ichg]);
        fHistPrimMCGenVtxZcut[index] = new TH3F(histName,";Centrality (%);#it{p}_{T} (GeV/#it{c});",
          nCentBins,centBins,nPtBins,ptBins,kNEvtCuts,evBins);

        fOutput->Add(fHistPrimMCGenVtxZall[index]);
        fOutput->Add(fHistPrimMCGenVtxZcut[index]);

        histName = Form("fHistMCReco%s%s", spcName[ispc], chgName[ichg]);
        fHistMCReco[index] = new TH3F(histName,";Centrality (%);#it{p}_{T} (GeV/#it{c});",
          nCentBins,centBins,nPtBins,ptBins,kNEvtCuts,evBins);
        fOutput->Add(fHistMCReco[index]);

        histName = Form("fHistMCPrimReco%s%s", spcName[ispc], chgName[ichg]);
        fHistMCPrimReco[index] = new TH3F(histName,";Centrality (%);#it{p}_{T} (GeV/#it{c});",
          nCentBins,centBins,nPtBins,ptBins,kNEvtCuts,evBins);
        fOutput->Add(fHistMCPrimReco[index]);

        //Histograms MC part Rec.
        histName = Form("fHistPrimMCReco%s%s", spcName[ispc], chgName[ichg]);
        fHistPrimMCReco[index] = new TH2F(histName,";Centrality (%);#it{p}_{T} (GeV/#it{c})",
          nCentBins,centBins,nPtBins,ptBins);
        histName = Form("fHistSstrMCReco%s%s", spcName[ispc], chgName[ichg]);
        fHistSstrMCReco[index] = new TH2F(histName,";Centrality (%);#it{p}_{T} (GeV/#it{c})",
          nCentBins,centBins,nPtBins,ptBins);
        histName = Form("fHistSmatMCReco%s%s", spcName[ispc], chgName[ichg]);
        fHistSmatMCReco[index] = new TH2F(histName,";Centrality (%);#it{p}_{T} (GeV/#it{c})",
          nCentBins,centBins,nPtBins,ptBins);

        fOutput->Add(fHistPrimMCReco[index]);
        fOutput->Add(fHistSstrMCReco[index]);
        fOutput->Add(fHistSmatMCReco[index]);

        //Histograms MC DCAxy
        histName = Form("fHistPrimDCA%s%s", spcName[ispc], chgName[ichg]);
        fHistPrimDCA[index] = new TH3F(histName,";Centrality (%);#it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)",
          nCentBins,centBins,nPtBins,ptBins,nDCABins,dcaBins);

        histName = Form("fHistSstrDCA%s%s", spcName[ispc], chgName[ichg]);
        fHistSstrDCA[index] = new TH3F(histName,";Centrality (%);#it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)",
          nCentBins,centBins,nPtBins,ptBins,nDCABins,dcaBins);

        histName = Form("fHistSmatDCA%s%s", spcName[ispc], chgName[ichg]);
        fHistSmatDCA[index] = new TH3F(histName,";Centrality (%);#it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)",
          nCentBins,centBins,nPtBins,ptBins,nDCABins,dcaBins);

        histName = Form("fHistMCtruthPrimDCA%s%s", spcName[ispc], chgName[ichg]);
        fHistMCtruthPrimDCA[index] = new TH3F(histName,";Centrality (%);#it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)",
          nCentBins,centBins,nPtBins,ptBins,nDCABins,dcaBins);

        histName = Form("fHistMCtruthSstrDCA%s%s", spcName[ispc], chgName[ichg]);
        fHistMCtruthSstrDCA[index] = new TH3F(histName,";Centrality (%);#it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)",
          nCentBins,centBins,nPtBins,ptBins,nDCABins,dcaBins);

        histName = Form("fHistMCtruthSmatDCA%s%s", spcName[ispc], chgName[ichg]);
        fHistMCtruthSmatDCA[index] = new TH3F(histName,";Centrality (%);#it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)",
          nCentBins,centBins,nPtBins,ptBins,nDCABins,dcaBins);

        fOutput->Add(fHistPrimDCA[index]);
        fOutput->Add(fHistSstrDCA[index]);
        fOutput->Add(fHistSmatDCA[index]);
        fOutput->Add(fHistMCtruthPrimDCA[index]);
        fOutput->Add(fHistMCtruthSstrDCA[index]);
        fOutput->Add(fHistMCtruthSmatDCA[index]);

      }//end IsMC
    }
  }

  for (Int_t i = 0; i < 4; i++) {
    fHistCharge[i] = new TH1F(Form("fHistChargeLay%d", i), Form("fHistChargeLay%d", i), 100, 0, 300);
    fOutput->Add(fHistCharge[i]);
  }

  if (fFillIntDistHist) {
 
    //dEdx distributions

    fHistPosHypPi = new TH2F("fHistPosHypPi", "fHistPosHypPi", nPtBins,ptBins,nDCABins,dcaBins); //DCA distr. with NSigma PID
    fHistPosHypKa = new TH2F("fHistPosHypKa", "fHistPosHypKa", nPtBins,ptBins,nDCABins,dcaBins);
    fHistPosHypPr = new TH2F("fHistPosHypPr", "fHistPosHypPr", nPtBins,ptBins,nDCABins,dcaBins);
    fHistNegHypPi = new TH2F("fHistNegHypPi", "fHistNegHypPi", nPtBins,ptBins,nDCABins,dcaBins);
    fHistNegHypKa = new TH2F("fHistNegHypKa", "fHistNegHypKa", nPtBins,ptBins,nDCABins,dcaBins);
    fHistNegHypPr = new TH2F("fHistNegHypPr", "fHistNegHypPr", nPtBins,ptBins,nDCABins,dcaBins);

    fOutput->Add(fHistPosHypPi); //DCA distr
    fOutput->Add(fHistPosHypKa);
    fOutput->Add(fHistPosHypPr);
    fOutput->Add(fHistNegHypPi);
    fOutput->Add(fHistNegHypKa);
    fOutput->Add(fHistNegHypPr);

    if (fIsMC) {
      const int nBins = 175;
      float dEdxBins[nBins + 1];
      SetBins(nBins,-3.5f,3.5f,dEdxBins);
      fHistMCPosOtherHypPion = new TH2F("fHistMCPosOtherHypPion", "fHistMCPosOtherHypPion",nPtBins,ptBins,nBins,dEdxBins); //MC truth
      fHistMCPosOtherHypKaon = new TH2F("fHistMCPosOtherHypKaon", "fHistMCPosOtherHypKaon",nPtBins,ptBins,nBins,dEdxBins);
      fHistMCPosOtherHypProt = new TH2F("fHistMCPosOtherHypProt", "fHistMCPosOtherHypProt",nPtBins,ptBins,nBins,dEdxBins);
      fHistMCPosElHypPion    = new TH2F("fHistMCPosElHypPion", "fHistMCPosElHypPion",nPtBins,ptBins,nBins,dEdxBins);
      fHistMCPosElHypKaon    = new TH2F("fHistMCPosElHypKaon", "fHistMCPosElHypKaon",nPtBins,ptBins,nBins,dEdxBins);
      fHistMCPosElHypProt    = new TH2F("fHistMCPosElHypProt", "fHistMCPosElHypProt",nPtBins,ptBins,nBins,dEdxBins);
      fHistMCPosPiHypPion    = new TH2F("fHistMCPosPiHypPion", "fHistMCPosPiHypPion",nPtBins,ptBins,nBins,dEdxBins);
      fHistMCPosPiHypKaon    = new TH2F("fHistMCPosPiHypKaon", "fHistMCPosPiHypKaon",nPtBins,ptBins,nBins,dEdxBins);
      fHistMCPosPiHypProt    = new TH2F("fHistMCPosPiHypProton", "fHistMCPosPiHypProton",nPtBins,ptBins,nBins,dEdxBins);
      fHistMCPosKaHypPion    = new TH2F("fHistMCPosKaHypPion", "fHistMCPosKaHypPion",nPtBins,ptBins,nBins,dEdxBins);
      fHistMCPosKaHypKaon    = new TH2F("fHistMCPosKaHypKaon", "fHistMCPosKaHypKaon",nPtBins,ptBins,nBins,dEdxBins);
      fHistMCPosKaHypProt    = new TH2F("fHistMCPosKaHypProt", "fHistMCPosKaHypProt",nPtBins,ptBins,nBins,dEdxBins);
      fHistMCPosPrHypPion    = new TH2F("fHistMCPosPrHypPion", "fHistMCPosPrHypPion",nPtBins,ptBins,nBins,dEdxBins);
      fHistMCPosPrHypKaon    = new TH2F("fHistMCPosPrHypKaon", "fHistMCPosPrHypKaon",nPtBins,ptBins,nBins,dEdxBins);
      fHistMCPosPrHypProt    = new TH2F("fHistMCPosPrHypProt", "fHistMCPosPrHypProt",nPtBins,ptBins,nBins,dEdxBins);

      fHistMCNegOtherHypPion = new TH2F("fHistMCNegOtherHypPion","fHistMCNegOtherHypPion",nPtBins,ptBins,nBins,dEdxBins); //MC truth
      fHistMCNegOtherHypKaon = new TH2F("fHistMCNegOtherHypKaon","fHistMCNegOtherHypKaon",nPtBins,ptBins,nBins,dEdxBins);
      fHistMCNegOtherHypProt = new TH2F("fHistMCNegOtherHypProt","fHistMCNegOtherHypProt",nPtBins,ptBins,nBins,dEdxBins);
      fHistMCNegElHypPion    = new TH2F("fHistMCNegElHypPion", "fHistMCNegElHypPion",nPtBins,ptBins,nBins,dEdxBins);
      fHistMCNegElHypKaon    = new TH2F("fHistMCNegElHypKaon", "fHistMCNegElHypKaon",nPtBins,ptBins,nBins,dEdxBins);
      fHistMCNegElHypProt    = new TH2F("fHistMCNegElHypProt", "fHistMCNegElHypProt",nPtBins,ptBins,nBins,dEdxBins);
      fHistMCNegPiHypPion    = new TH2F("fHistMCNegPiHypPion", "fHistMCNegPiHypPion",nPtBins,ptBins,nBins,dEdxBins);
      fHistMCNegPiHypKaon    = new TH2F("fHistMCNegPiHypKaon", "fHistMCNegPiHypKaon",nPtBins,ptBins,nBins,dEdxBins);
      fHistMCNegPiHypProt    = new TH2F("fHistMCNegPiHypProt", "fHistMCNegPiHypProt",nPtBins,ptBins,nBins,dEdxBins);
      fHistMCNegKaHypPion    = new TH2F("fHistMCNegKaHypPion", "fHistMCNegKaHypPion",nPtBins,ptBins,nBins,dEdxBins);
      fHistMCNegKaHypKaon    = new TH2F("fHistMCNegKaHypKaon", "fHistMCNegKaHypKaon",nPtBins,ptBins,nBins,dEdxBins);
      fHistMCNegKaHypProt    = new TH2F("fHistMCNegKaHypProt", "fHistMCNegKaHypProt",nPtBins,ptBins,nBins,dEdxBins);
      fHistMCNegPrHypPion    = new TH2F("fHistMCNegPrHypPion", "fHistMCNegPrHypPion",nPtBins,ptBins,nBins,dEdxBins);
      fHistMCNegPrHypKaon    = new TH2F("fHistMCNegPrHypKaon", "fHistMCNegPrHypKaon",nPtBins,ptBins,nBins,dEdxBins);
      fHistMCNegPrHypProt    = new TH2F("fHistMCNegPrHypProt", "fHistMCNegPrHypProt",nPtBins,ptBins,nBins,dEdxBins);

      fOutput->Add(fHistMCPosOtherHypPion);//MC truth
      fOutput->Add(fHistMCPosOtherHypKaon);
      fOutput->Add(fHistMCPosOtherHypProt);
      fOutput->Add(fHistMCPosElHypPion);
      fOutput->Add(fHistMCPosElHypKaon);
      fOutput->Add(fHistMCPosElHypProt);
      fOutput->Add(fHistMCPosPiHypPion);
      fOutput->Add(fHistMCPosPiHypKaon);
      fOutput->Add(fHistMCPosPiHypProt);
      fOutput->Add(fHistMCPosKaHypPion);
      fOutput->Add(fHistMCPosKaHypKaon);
      fOutput->Add(fHistMCPosKaHypProt);
      fOutput->Add(fHistMCPosPrHypPion);
      fOutput->Add(fHistMCPosPrHypKaon);
      fOutput->Add(fHistMCPosPrHypProt);

      fOutput->Add(fHistMCNegOtherHypPion);//MC truth
      fOutput->Add(fHistMCNegOtherHypKaon);
      fOutput->Add(fHistMCNegOtherHypProt);
      fOutput->Add(fHistMCNegElHypPion);
      fOutput->Add(fHistMCNegElHypKaon);
      fOutput->Add(fHistMCNegElHypProt);
      fOutput->Add(fHistMCNegPiHypPion);
      fOutput->Add(fHistMCNegPiHypKaon);
      fOutput->Add(fHistMCNegPiHypProt);
      fOutput->Add(fHistMCNegKaHypPion);
      fOutput->Add(fHistMCNegKaHypKaon);
      fOutput->Add(fHistMCNegKaHypProt);
      fOutput->Add(fHistMCNegPrHypPion);
      fOutput->Add(fHistMCNegPrHypKaon);
      fOutput->Add(fHistMCNegPrHypProt);
    }//end IsMC
  }//end FillIntDistHist

  PostData(1, fOutput);

  CreateDCAcutFunctions(); //Creating kParamContainer data
  PostData(2, fListCuts);

  // Post output data container
  if (fFillNtuple) {
    fListTree = new TList();
    fListTree->SetOwner();

    fNtupleData = new TNtuple("fNtupleData", "fNtupleData", "mult:p:pt:s0:s1:s2:s3:dEdx:sign:eta:dcaXY:dcaZ:clumap");
    fListTree->Add(fNtupleData);
    fNtupleMC   = new TNtuple("fNtupleMC", "fNtupleMC", "mult:mcPt:pdgcode:sign:mcEta:mcRap:isph:run");
    fListTree->Add(fNtupleMC);
  }
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
  Printf("Inizializing Task, be sure to run after all configuration have been set...");
  if (!fUseDefaultPriors) DefineInput(1, TList::Class());
  if (fFillNtuple) DefineOutput(3, TList::Class());
  AliInfo("Tracks selections");
  AliInfoF(" y = yLab + %.3f,  Ymin %.1f, Ymax %.1f, Eabs %.1f, DCAxyCut %.1f, DCAzCut %.1f, Chi2 %.1f,   nSPD %d,   nPID %d",
           fCMSRapFct, fMinRapCut, fMaxRapCut, fAbsEtaCut, fNSigmaDCAxy, fNSigmaDCAz, fMaxChi2Clu, fMinSPDPts, fMinNdEdxSamples);

  if (fMultMethod)
    AliInfoF("Cent. %.f %.f %s", fLowMult, fUpMult, fMultEstimator.Data());

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
    fHistNEvents->Fill(fEvtMult,istep);
    if (lHasGoodVtxGen) fHistMCEvents->Fill(fEvtMult,istep);
  }

  if (!lIsEventSelected)  {
    AliDebugF(3, "Event rejected by event selection after %d step", (int)lastEvtCutPassed);
    PostAllData();
    return;
  }
	fHistNEvents->Fill(fEvtMult,kNEvtCuts);
  if (lHasGoodVtxGen) fHistMCEvents->Fill(fEvtMult,kNEvtCuts);

  if (fMultMethod) //Fill fHistMultAftEvtSel after the event Selection
    fHistMultAftEvtSel->Fill(fEvtMult);

  if (!fITSPIDResponse)
    fITSPIDResponse = new AliITSPIDResponse(fIsMC);

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
    fHistNTracks[iChg]->Fill(fEvtMult,trkPt,trkSel);

    //"ITSsa"
    if (!(status & AliESDtrack::kITSpureSA)) continue;
    trkSel = kIsITSsa;
    fHistNTracks[iChg]->Fill(fEvtMult,trkPt,trkSel);

    //"ITSrefit"
    if (!(status & AliESDtrack::kITSrefit)) continue;
    trkSel = kIsITSrefit;
    fHistNTracks[iChg]->Fill(fEvtMult,trkPt,trkSel);

    //"neutral particle"
    if (TMath::Abs(track->GetSign()) < 0.0001) continue;
    trkSel = kIsNotNeutralParticle;
    fHistNTracks[iChg]->Fill(fEvtMult,trkPt,trkSel);

    //"SPDcls"
    if (nSPD < fMinSPDPts) continue;//At least one point in the SPD
    trkSel = kPassSPD;
    fHistNTracks[iChg]->Fill(fEvtMult,trkPt,trkSel);

    //"SDD+SSD cls" at least 3 points on SSD/SDD
    if (nPtsForPid < fMinNdEdxSamples) continue;
    trkSel = kPassPIDcls;
    fHistNTracks[iChg]->Fill(fEvtMult,trkPt,trkSel);

    //"chi2/ncls"->chisquare/nclusters
    Int_t nclu = nSPD + nPtsForPid;
    if (track->GetITSchi2() / nclu > fMaxChi2Clu) continue;
    trkSel = kPassChi2Ncls;
    fHistNTracks[iChg]->Fill(fEvtMult,trkPt,trkSel);

    //"eta"->pseudorapidity
    if (TMath::Abs(track->Eta()) > fAbsEtaCut) continue;
    trkSel = kIsInEta;
    fHistNTracks[iChg]->Fill(fEvtMult,trkPt,trkSel);

    //"dE/dx < 0"->truncated mean
    if (dEdx < 0) continue;
    trkSel = kPassdEdx;
    fHistNTracks[iChg]->Fill(fEvtMult,trkPt,trkSel);

    //fill propaganda plot with dedx before pt cut
    fHistDEDX->Fill(track->GetP(), dEdx);
    fHistDEDXdouble->Fill(track->GetP()*track->GetSign(), dEdx);

    //"ptCut"
    if ((trkPt < fPtBins[0]) || (trkPt >= fPtBins[fPtBins.GetSize()-1])) continue;
    trkSel = kPassPtCut;
    fHistNTracks[iChg]->Fill(fEvtMult,trkPt,trkSel);

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
    fHistNTracks[iChg]->Fill(fEvtMult,trkPt,trkSel);

    if (fFillNtuple) {
      Float_t xnt[12];
      Int_t index = 0;
      /*1 */ xnt[index++] = (float)fEvtMult;
      /*2 */ xnt[index++] = (float)track->GetP();
      /*3 */ xnt[index++] = (float)track->Pt();
      /*4 */ xnt[index++] = (float)dEdxLay[0];
      /*5 */ xnt[index++] = (float)dEdxLay[1];
      /*6 */ xnt[index++] = (float)dEdxLay[2];
      /*7 */ xnt[index++] = (float)dEdxLay[3];
      /*8 */ xnt[index++] = (float)dEdx;
      /*9*/ xnt[index++]  = (float)track->GetSign();
      /*10*/ xnt[index++]  = (float)track->Eta();
      /*11*/ xnt[index++] = (float)impactXY;
      /*12*/ xnt[index++] = (float)impactZ;
      /*13*/ xnt[index++] = (float)clumap;

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
        fHistRecoDCA[lPidIndex]->Fill(fEvtMult,trkPt,impactXY);

        //DCA distributions, before the DCAxy cuts from the MC kinematics
        //Filling DCA distribution with MC truth Physics values
        if (fIsMC) {
          if (lMCevent->IsPhysicalPrimary(lMCtrk))        fHistPrimDCA[lPidIndex]->Fill(fEvtMult,trkPt,impactXY);
          if (lMCevent->IsSecondaryFromWeakDecay(lMCtrk)) fHistSstrDCA[lPidIndex]->Fill(fEvtMult,trkPt,impactXY);
          if (lMCevent->IsSecondaryFromMaterial(lMCtrk))  fHistSmatDCA[lPidIndex]->Fill(fEvtMult,trkPt,impactXY);
        }
      }// end lIsGoodTrack
      if (fIsMC && lIsGoodPart) {
        //DCA distributions, before the DCAxy cuts from the MC kinematics
        //Filling DCA distribution with MC truth Physics values
        if (lMCevent->IsPhysicalPrimary(lMCtrk))        fHistMCtruthPrimDCA[lMCtIndex]->Fill(fEvtMult,trkPt,impactXY);
        if (lMCevent->IsSecondaryFromWeakDecay(lMCtrk)) fHistMCtruthSstrDCA[lMCtIndex]->Fill(fEvtMult,trkPt,impactXY);
        if (lMCevent->IsSecondaryFromMaterial(lMCtrk))  fHistMCtruthSmatDCA[lMCtIndex]->Fill(fEvtMult,trkPt,impactXY);
      }// end lIsGoodPart
      //"DCAxy"
      if (!DCAcutXY(impactXY, trkPt)) continue;
      trkSel = kPassDCAxycut;
      fHistNTracks[iChg]->Fill(fEvtMult,trkPt,trkSel);

      if (lIsGoodTrack) fHistReco[lPidIndex]->Fill(fEvtMult,trkPt);
      //Filling Histos for Reco Efficiency
      //information from the MC kinematics
      if (fIsMC && lIsGoodPart) {
        if (lMCevent->IsPhysicalPrimary(lMCtrk))         fHistPrimMCReco[lMCtIndex]->Fill(fEvtMult,trkPt);
        if (lMCevent->IsSecondaryFromWeakDecay(lMCtrk))  fHistSstrMCReco[lMCtIndex]->Fill(fEvtMult,trkPt);
        if (lMCevent->IsSecondaryFromMaterial(lMCtrk))   fHistSmatMCReco[lMCtIndex]->Fill(fEvtMult,trkPt);
      }

      Int_t binPart = (lMCspc > AliPID::kMuon) ? (lMCspc - 2) : -1;
      if (fIsMC && lIsGoodTrack) {
        fHistMCReco[lPidIndex]->Fill(fEvtMult,fEvtMult,trkPt, binPart);
        if (lMCevent->IsPhysicalPrimary(lMCtrk))
          fHistMCPrimReco[lPidIndex]->Fill(fEvtMult,trkPt, binPart);
      }//end y

      if (lIsGoodTrack && fFillIntDistHist) {
        //
        //integral approach histograms
        //
        Int_t lAbsPdgCode = TMath::Abs(lMCpdg);
        if (track->GetSign() > 0) {
          if (IsRapIn(y[AliPID::kPion]))  fHistPosHypPi->Fill(trkPt,logdiff[0]);
          if (IsRapIn(y[AliPID::kKaon]))  fHistPosHypKa->Fill(trkPt,logdiff[1]);
          if (IsRapIn(y[AliPID::kProton]))fHistPosHypPr->Fill(trkPt,logdiff[2]);
          if (fIsMC) {
            if (IsRapIn(y[AliPID::kPion])) {
              if ((lAbsPdgCode != 11) && (lAbsPdgCode != 211) && (lAbsPdgCode != 321) && (lAbsPdgCode != 2212))
                fHistMCPosOtherHypPion->Fill(trkPt,logdiff[0]);
              if (lAbsPdgCode == 11)  fHistMCPosElHypPion->Fill(trkPt,logdiff[0]);
              if (lAbsPdgCode == 211) fHistMCPosPiHypPion->Fill(trkPt,logdiff[0]);
              if (lAbsPdgCode == 321) fHistMCPosKaHypPion->Fill(trkPt,logdiff[0]);
              if (lAbsPdgCode == 2212)fHistMCPosPrHypPion->Fill(trkPt,logdiff[0]);
            }
            if (IsRapIn(y[AliPID::kKaon])) {
              if ((lAbsPdgCode != 11) && (lAbsPdgCode != 211) && (lAbsPdgCode != 321) && (lAbsPdgCode != 2212))
                fHistMCPosOtherHypKaon->Fill(trkPt,logdiff[1]);
              if (lAbsPdgCode == 11)  fHistMCPosElHypKaon->Fill(trkPt,logdiff[1]);
              if (lAbsPdgCode == 211) fHistMCPosPiHypKaon->Fill(trkPt,logdiff[1]);
              if (lAbsPdgCode == 321) fHistMCPosKaHypKaon->Fill(trkPt,logdiff[1]);
              if (lAbsPdgCode == 2212)fHistMCPosPrHypKaon->Fill(trkPt,logdiff[1]);
            }
            if (IsRapIn(y[AliPID::kProton])) {
              if ((lAbsPdgCode != 11) && (lAbsPdgCode != 211) && (lAbsPdgCode != 321) && (lAbsPdgCode != 2212))
                fHistMCPosOtherHypProt->Fill(trkPt,logdiff[2]);
              if (lAbsPdgCode == 11)  fHistMCPosElHypProt->Fill(trkPt,logdiff[2]);
              if (lAbsPdgCode == 211) fHistMCPosPiHypProt->Fill(trkPt,logdiff[2]);
              if (lAbsPdgCode == 321) fHistMCPosKaHypProt->Fill(trkPt,logdiff[2]);
              if (lAbsPdgCode == 2212)fHistMCPosPrHypProt->Fill(trkPt,logdiff[2]);
            }
          }
        } else {
          if (IsRapIn(y[AliPID::kPion]))  fHistNegHypPi->Fill(trkPt,logdiff[0]);
          if (IsRapIn(y[AliPID::kKaon]))  fHistNegHypKa->Fill(trkPt,logdiff[1]);
          if (IsRapIn(y[AliPID::kProton]))fHistNegHypPr->Fill(trkPt,logdiff[2]);
          if (fIsMC) {
            if (IsRapIn(y[AliPID::kPion])) {
              if ((lAbsPdgCode != 11) && (lAbsPdgCode != 211) && (lAbsPdgCode != 321) && (lAbsPdgCode != 2212))
                fHistMCNegOtherHypPion->Fill(trkPt,logdiff[0]);
              if (lAbsPdgCode == 11)  fHistMCNegElHypPion->Fill(trkPt,logdiff[0]);
              if (lAbsPdgCode == 211) fHistMCNegPiHypPion->Fill(trkPt,logdiff[0]);
              if (lAbsPdgCode == 321) fHistMCNegKaHypPion->Fill(trkPt,logdiff[0]);
              if (lAbsPdgCode == 2212)fHistMCNegPrHypPion->Fill(trkPt,logdiff[0]);
            }
            if (IsRapIn(y[AliPID::kKaon])) {
              if ((lAbsPdgCode != 11) && (lAbsPdgCode != 211) && (lAbsPdgCode != 321) && (lAbsPdgCode != 2212))
                fHistMCNegOtherHypKaon->Fill(trkPt,logdiff[1]);
              if (lAbsPdgCode == 11)  fHistMCNegElHypKaon->Fill(trkPt,logdiff[1]);
              if (lAbsPdgCode == 211) fHistMCNegPiHypKaon->Fill(trkPt,logdiff[1]);
              if (lAbsPdgCode == 321) fHistMCNegKaHypKaon->Fill(trkPt,logdiff[1]);
              if (lAbsPdgCode == 2212)fHistMCNegPrHypKaon->Fill(trkPt,logdiff[1]);
            }
            if (IsRapIn(y[AliPID::kProton])) {
              if ((lAbsPdgCode != 11) && (lAbsPdgCode != 211) && (lAbsPdgCode != 321) && (lAbsPdgCode != 2212))
                fHistMCNegOtherHypProt->Fill(trkPt,logdiff[2]);
              if (lAbsPdgCode == 11)  fHistMCNegElHypProt->Fill(trkPt,logdiff[2]);
              if (lAbsPdgCode == 211) fHistMCNegPiHypProt->Fill(trkPt,logdiff[2]);
              if (lAbsPdgCode == 321) fHistMCNegKaHypProt->Fill(trkPt,logdiff[2]);
              if (lAbsPdgCode == 2212)fHistMCNegPrHypProt->Fill(trkPt,logdiff[2]);
            }
          }
        }
      }//end lIsGooTrack
    }//end else fill Ntuple
  }//end track loop

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
Bool_t AliAnalysisTaskSEITSsaSpectra::IsEventAccepted(EEvtCut_Type& evtSel)
{
	//Check multiplicity selection first
	if (!IsMultSelected()){
		AliDebug(3, "Event doesn't pass multiplicity selection");
    PostAllData();
    evtSel = kPassMultSel;
    return kFALSE;
	} else {
		fHistMultBefEvtSel->Fill(fEvtMult);
	}

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

  if (fExtEventCuts) {
    if (fEventCuts.AcceptEvent(fESD)) {
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
	    evtSel = kIsPileupMV;
      return kFALSE;
    }

		if (!fEventCuts.PassedCut(AliEventCuts::kCorrelations)) {
      AliDebug(3, "Event with PileUp");
      PostAllData();
    	evtSel = kCorrelations;
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

  } else {  //If not Eventcut used

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

		if(fChkIsEventINELgtZERO && !(AliESDtrackCuts::GetReferenceMultiplicity(fESD, AliESDtrackCuts::kTracklets, 1.0) >= 1)){
    	AliDebug(3, "Event doesn't pass IsEventINELgtZERO selection");
    	PostAllData();
	  	evtSel = kPassINELgtZERO;
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

    if (! (fPlpType & BIT(kNoPileup)) ) {
			UInt_t lFlag = IsPileup();
			if (lFlag & BIT(kPileupSPD)) {
      	AliDebug(3, "Pileup event from IsPileupFromSPD");
     	 	PostAllData();
    		evtSel = kIsPileupSPD;
      	return kFALSE;
			}
			if (lFlag & BIT(kPileupInMultBins)) {
      	AliDebug(3, "Pileup event from IsPileupSPDinMultBins");
     	 	PostAllData();
    		evtSel = kIsPileupSPDinMultBins;
      	return kFALSE;
			}
			if (lFlag & BIT(kPileupMV)) {
      	AliDebug(3, "Pileup event from IsPileupMV");
     	 	PostAllData();
    		evtSel = kIsPileupSPD;
      	return kFALSE;
			}
    }

    if (fUseSelectVertex2015pp && !SelectVertex2015pp()) {
      AliDebug(3, "Event doesn't pass vtx 2015 pp selection sel");
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
  fPlpType      = BIT(kNoPileup);
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
  fPlpType      = BIT(kNoPileup);
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
      fEvtMult = fMultSel->GetMultiplicityPercentile(fMultEstimator.Data(),fMultEvSel);
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
UInt_t AliAnalysisTaskSEITSsaSpectra::IsPileup()
{
	UInt_t lReturn = 0u;

	if      (fPlpType & BIT(kNoPileup)  )
		return 1u; //skip pileup check;
	else if (fPlpType & BIT(kPileupSPD) ) {
		if (fESD->IsPileupFromSPD(fMinPlpContribSPD,fMinPlpZdistSPD))
			lReturn |= BIT(kPileupSPD);
	}
	else if (fPlpType & BIT(kPileupInMultBins)) {
		if (fESD->IsPileupFromSPDInMultBins())
			lReturn |= BIT(kPileupInMultBins);;
	}
	else if (fPlpType & BIT(kPileupMV)  ) {
		AliAnalysisUtils utils;
    utils.SetMinPlpContribMV(fMinPlpContribMV);
    utils.SetMaxPlpChi2MV(fMaxPlpChi2MV);
    utils.SetMinWDistMV(fMinWDistMV);
    utils.SetCheckPlpFromDifferentBCMV(fCheckPlpFromDifferentBCMV);
    utils.SetUseMVPlpSelection(kTRUE);

		if (utils.IsPileUpEvent(fESD))
			lReturn |= BIT(kPileupMV);
	}

  return lReturn;
}

//
//
//________________________________________________________________________
Bool_t AliAnalysisTaskSEITSsaSpectra::IsGoodVtxZ()
{
  const AliESDVertex* vtx = fESD->GetPrimaryVertex();
  fHistVtxZ->Fill(fEvtMult,vtx->GetZ());
  if (TMath::Abs(vtx->GetZ()) > fMaxVtxZCut) return kFALSE;

  return kTRUE;
}

//
//
//________________________________________________________________________
Bool_t AliAnalysisTaskSEITSsaSpectra::IsGoodSPDvtxRes(const AliESDVertex* spdVtx)
{
  if (!spdVtx) return kFALSE;
  if (spdVtx->IsFromVertexerZ() && !(spdVtx->GetDispersion() < 0.04 || spdVtx->GetZRes() < 0.25)) return kFALSE;
  return kTRUE;
}

//
//
//________________________________________________________________________
Bool_t AliAnalysisTaskSEITSsaSpectra::SelectVertex2015pp()
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
    if (mcPt > 1.0) continue; // pt cut

    Double_t mcEta = mcTrk->Eta();
    Double_t mcRap = Eta2y(mcPt, mcTrk->GetMass(), mcEta) + fCMSRapFct;

    Bool_t lIsPhysPrimary = lMCevent->IsPhysicalPrimary(i_mcTrk);
    if (fFillNtuple) {
      //filling MC ntuple
      Float_t xntMC[8];
      Int_t indexMC = 0;
      xntMC[indexMC++] = (Float_t)fEvtMult;
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
      fHistPrimMCGenVtxZall[index]->Fill(fEvtMult,TMath::Abs(mcPt), istep);
      if (lHasGoodVtxGen)
        fHistPrimMCGenVtxZcut[index]->Fill(fEvtMult,TMath::Abs(mcPt), istep);
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

//
//
//________________________________________________________________________
void AliAnalysisTaskSEITSsaSpectra::SetBins(int nbins, float min, float max, float* bins) {
  const float delta = (max - min) / nbins;
  for (int iB = 0; iB < nbins; ++iB) {
    bins[iB] = min + iB * delta;
  }
  bins[nbins] = max;
}

//
//
//________________________________________________________________________
void AliAnalysisTaskSEITSsaSpectra::SetCentBins(int nbins, float *bins) {
  fCentBins.Set(nbins + 1, bins);
}


//
//
//________________________________________________________________________
void AliAnalysisTaskSEITSsaSpectra::SetDCABins(int nbins, float min, float max) {
  const float delta = (max - min) / nbins;
  fDCABins.Set(nbins + 1);
  for (int iB = 0; iB < nbins; ++iB) {
    fDCABins[iB] = min + iB * delta;
  }
  fDCABins[nbins] = max;
}

//
//
//________________________________________________________________________
void AliAnalysisTaskSEITSsaSpectra::SetDCABins(int nbins, float *bins) {
  fDCABins.Set(nbins + 1, bins);
}

//
//
//________________________________________________________________________
void AliAnalysisTaskSEITSsaSpectra::SetPtBins(int nbins, float *bins) {
  fPtBins.Set(nbins + 1, bins);
}
