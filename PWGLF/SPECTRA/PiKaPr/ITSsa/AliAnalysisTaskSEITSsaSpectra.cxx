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

#include <string>

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
  AliAnalysisTaskSEITSsaSpectra::AliAnalysisTaskSEITSsaSpectra(bool __def_prior, bool __fill_ntuple)
  : AliAnalysisTaskSE("TaskITSsaSpectra"),
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
    fHistDEDXGen(NULL),
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
    fMaxRapCut(.5),
    fCMSRapFct(.0),
    fMindEdx(0.),
    fMinNSigma(1.5),
    fMaxChi2Clu(2.5),
    fNSigmaDCAxy(7.),
    fNSigmaDCAz(7.),
    fYear(2010),
    fPidMethod(kMeanCut),
    fUseDefaultPriors(__def_prior),
    fFillNtuple(__fill_ntuple),
    fIsMC(kFALSE),
    fIsDCAUnfoldHistoEnabled(kFALSE),
    fIsNominalBfield(kTRUE),
    fFillIntDistHist(kFALSE),
    fRandGener(0x0),
    fSmearMC(kFALSE),
    fSmearP(0.),
    fSmeardEdx(0.)
{
  // Constructor
  fRandGener = new TRandom3(0);

  // Inizialize bins array
  fCentBins.Set(3);
  fCentBins[0] = -5.f;
  fCentBins[1] = 0.f;
  fCentBins[2] = 100.f;

  double ptBins[kNbins + 1] = { 0.08, 0.10, 0.12, 0.14, 0.16, 0.18, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45,
                               0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.0 };
  fPtBins.Set(kNbins + 1, ptBins);
  const int nDCAbins = 2000;
  double dcaBins[nDCAbins];
  SetBins(nDCAbins, -2, 2, dcaBins);
  SetDCABins(nDCAbins, dcaBins);

  for (int i_chg = 0; i_chg < kNchg; ++i_chg)
    fHistNTracks[i_chg] = NULL;

  for (int i_spc = 0; i_spc < kNspc; ++i_spc) {
    for (int i_chg = 0; i_chg < kNchg; ++i_chg) {
      int index = i_spc * kNchg + i_chg;
      fHistNSigmaSep[index] = NULL;
      fHistSepPowerReco[index] = NULL;
      fHistSepPowerTrue[index] = NULL;

      fHistMCPart[index] = NULL;
      fHistMCPartGoodGenVtxZ[index] = NULL;

      fHistReco[index] = NULL;
      fHistRecoMC[index] = NULL;
      fHistRecoTrueMC[index] = NULL;
      fHistMCDCA[index] = NULL;

      fHistTruePIDMCReco[index] = NULL;
      fHistTruePIDMCGen[index] = NULL;
    }
  }
  for (int iL = 0; iL < 4; ++iL)
    fHistCharge[iL] = NULL;

  fHistMCGenCharged = NULL;
  fHistRecoChargedMC = NULL;

  // dEdx distributions
  fHistPosHypPi = NULL;
  fHistPosHypKa = NULL;
  fHistPosHypPr = NULL;
  fHistNegHypPi = NULL;
  fHistNegHypKa = NULL;
  fHistNegHypPr = NULL;

  // dEdx distributions for MC
  fHistMCPosOtherHypPion = NULL;
  fHistMCPosOtherHypKaon = NULL;
  fHistMCPosOtherHypProt = NULL;
  fHistMCPosElHypPion = NULL;
  fHistMCPosElHypKaon = NULL;
  fHistMCPosElHypProt = NULL;
  fHistMCPosPiHypPion = NULL;
  fHistMCPosPiHypKaon = NULL;
  fHistMCPosPiHypProt = NULL;
  fHistMCPosKaHypPion = NULL;
  fHistMCPosKaHypKaon = NULL;
  fHistMCPosKaHypProt = NULL;
  fHistMCPosPrHypPion = NULL;
  fHistMCPosPrHypKaon = NULL;
  fHistMCPosPrHypProt = NULL;

  fHistMCNegOtherHypPion = NULL;
  fHistMCNegOtherHypKaon = NULL;
  fHistMCNegOtherHypProt = NULL;
  fHistMCNegElHypPion = NULL;
  fHistMCNegElHypKaon = NULL;
  fHistMCNegElHypProt = NULL;
  fHistMCNegPiHypPion = NULL;
  fHistMCNegPiHypKaon = NULL;
  fHistMCNegPiHypProt = NULL;
  fHistMCNegKaHypPion = NULL;
  fHistMCNegKaHypKaon = NULL;
  fHistMCNegKaHypProt = NULL;
  fHistMCNegPrHypPion = NULL;
  fHistMCNegPrHypKaon = NULL;
  fHistMCNegPrHypProt = NULL;

  //Define input
  DefineInput(0, TChain::Class());
  if (!fUseDefaultPriors) DefineInput(1, TList::Class());

  //Define output
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
  if (!fUseDefaultPriors) {
    fListPriors = dynamic_cast<TList *>(GetInputData(1)); // FIXME
    if (!fListPriors) {
      AliDebug(3, "Errors use priors set but no priors list found");
      PostAllData();
      return;
    }
  }

  fOutput = new TList();
  fOutput->SetOwner();
  fOutput->SetName("Spiderman");

  if (fExtEventCuts) { // configure AliEventCuts
    // fEventCuts.SetManualMode();
    // Histograms for event Selection
    fEventCuts.AddQAplotsToList(fOutput);
  }

  const int nPtBins = fPtBins.GetSize() - 1;
  const int nCentBins = fCentBins.GetSize() - 1;
  const int nDCABins = fDCABins.GetSize() - 1;
  const double *ptBins = fPtBins.GetArray();
  const double *centBins = fCentBins.GetArray();
  const double *dcaBins = fDCABins.GetArray();

  double evBins[kNEvtCuts + 1];
  SetBins(kNEvtCuts, .5, kNEvtCuts + .5, evBins);

  const char *notApp = "_notApplied";
  fHistNEvents =
    new TH2I("fHistNEvents", "Number of processed events;Centrality (%);", nCentBins, centBins, kNEvtCuts, evBins);
  fHistNEvents->Sumw2();
  fHistNEvents->SetMinimum(0);
  fHistNEvents->GetYaxis()->SetBinLabel(kIsReadable, "Readable");
  fHistNEvents->GetYaxis()->SetBinLabel(kPassMultSel, Form("PassMultSel%s", (fMultMethod ? "" : notApp)));
  fHistNEvents->GetYaxis()->SetBinLabel(kIsSDDIn, Form("HasSDDIn%s", (fChkIsSDDIn ? "" : notApp)));
  fHistNEvents->GetYaxis()->SetBinLabel(kIsNotIncDAQ, Form("PassIncDAQ%s", (fRejIncDAQ ? "" : notApp)));
  fHistNEvents->GetYaxis()->SetBinLabel(kPassTrig, "PassPhysSelTrig");
  fHistNEvents->GetYaxis()->SetBinLabel(kPassINELgtZERO,
                                        Form("PassINELgtZERO%s", (fChkIsEventINELgtZERO ? "" : notApp)));
  fHistNEvents->GetYaxis()->SetBinLabel(kCorrelations, Form("Correlations%s", (fExtEventCuts ? "" : notApp)));
  fHistNEvents->GetYaxis()->SetBinLabel(kPassSPDclsVsTCut,
                                        Form("PassClsVsTrackletBG%s", (fDoSPDCvsTCut ? "" : notApp)));
  fHistNEvents->GetYaxis()->SetBinLabel(kIsNotPileupSPD,
                                        Form("PassIsPileupSPD%s", ((fPlpType & BIT(kPileupSPD) ? "" : notApp))));
  fHistNEvents->GetYaxis()->SetBinLabel(
    kIsNotPileupSPDinMultBins, Form("PassIsPileupSPDinMultBins%s", ((fPlpType & BIT(kPileupInMultBins) ? "" : notApp))));
  fHistNEvents->GetYaxis()->SetBinLabel(kIsNotPileupMV,
                                        Form("PassIsPileupMV%s", ((fPlpType & BIT(kPileupMV) ? "" : notApp))));
  fHistNEvents->GetYaxis()->SetBinLabel(kIsNotPileup, "PassIsNoPileup");
  fHistNEvents->GetYaxis()->SetBinLabel(kHasRecVtx, "HasVertex");
  fHistNEvents->GetYaxis()->SetBinLabel(kHasGoodVtxZ, "HasGoodVertex");
  fHistNEvents->GetYaxis()->SetBinLabel(kNEvtCuts, "IsSelected");
  fOutput->Add(fHistNEvents);

  fHistMCEvents =
    new TH2I("fHistMCEvents", "Number of processed events;Centrality (%);", nCentBins, centBins, kNEvtCuts, evBins);
  fHistMCEvents->Sumw2();
  fHistMCEvents->SetMinimum(0);
  fHistMCEvents->GetYaxis()->SetBinLabel(kIsReadable, "Readable");
  fHistMCEvents->GetYaxis()->SetBinLabel(kPassMultSel, Form("PassMultSel%s", (fMultMethod ? "" : notApp)));
  fHistMCEvents->GetYaxis()->SetBinLabel(kIsSDDIn, Form("HasSDDIn%s", (fChkIsSDDIn ? "" : notApp)));
  fHistMCEvents->GetYaxis()->SetBinLabel(kIsNotIncDAQ, Form("PassIncDAQ%s", (fRejIncDAQ ? "" : notApp)));
  fHistMCEvents->GetYaxis()->SetBinLabel(kPassTrig, "PassPhysSelTrig");
  fHistMCEvents->GetYaxis()->SetBinLabel(kPassINELgtZERO,
                                         Form("PassINELgtZERO%s", (fChkIsEventINELgtZERO ? "" : notApp)));
  fHistMCEvents->GetYaxis()->SetBinLabel(kCorrelations, Form("Correlations%s", (fExtEventCuts ? "" : notApp)));
  fHistMCEvents->GetYaxis()->SetBinLabel(kPassSPDclsVsTCut,
                                         Form("PassClsVsTrackletBG%s", (fDoSPDCvsTCut ? "" : notApp)));
  fHistMCEvents->GetYaxis()->SetBinLabel(kIsNotPileupSPD,
                                         Form("PassIsPileupSPD%s", ((fPlpType & BIT(kPileupSPD) ? "" : notApp))));
  fHistMCEvents->GetYaxis()->SetBinLabel(
    kIsNotPileupSPDinMultBins, Form("PassIsPileupSPDinMultBins%s", ((fPlpType & BIT(kPileupInMultBins) ? "" : notApp))));
  fHistMCEvents->GetYaxis()->SetBinLabel(kIsNotPileupMV,
                                         Form("PassIsPileupMV%s", ((fPlpType & BIT(kPileupMV) ? "" : notApp))));
  fHistMCEvents->GetYaxis()->SetBinLabel(kIsNotPileup, "PassIsNoPileup");
  fHistMCEvents->GetYaxis()->SetBinLabel(kHasRecVtx, "HasVertex");
  fHistMCEvents->GetYaxis()->SetBinLabel(kHasGoodVtxZ, "HasGoodVertex");
  fHistMCEvents->GetYaxis()->SetBinLabel(kNEvtCuts, "IsSelected");
  fOutput->Add(fHistMCEvents);

  fHistMultBefEvtSel =
    new TH1F("fHistMultBefEvtSel", "Event Multiplicity before event selection;Centrality (%)", nCentBins, centBins);
  fHistMultBefEvtSel->Sumw2();
  fHistMultBefEvtSel->SetMinimum(0);
  fOutput->Add(fHistMultBefEvtSel);

  fHistMultAftEvtSel =
    new TH1F("fHistMultAftEvtSel", "Event Multiplicity after event selection;Centrality (%)", nCentBins, centBins);
  fHistMultAftEvtSel->Sumw2();
  fHistMultAftEvtSel->SetMinimum(0);
  fOutput->Add(fHistMultAftEvtSel);

  const int nVtxBins = 400;
  double vtxBins[nVtxBins + 1];
  SetBins(nVtxBins, -20, 20, vtxBins);

  fHistVtxZ = new TH2F("fHistVtxZ", "Vtx Z distribution;Centrality (%);Z_vtx", nCentBins, centBins, nVtxBins, vtxBins);
  fHistVtxZ->Sumw2();
  fHistVtxZ->SetMinimum(0);
  fOutput->Add(fHistVtxZ);

  std::string spc_name[kNspc] = { "Pi", "Ka", "Pr" };
  std::string chg_name[kNchg] = { "Pos", "Neg" };

  std::string hist_name;

  const int nTrkBins = 20;
  double trkBins[nTrkBins + 1];
  SetBins(nTrkBins, .5, nTrkBins + .5, trkBins);

  // Histo with track cuts
  for (int i_chg = 0; i_chg < kNchg; ++i_chg) {
    hist_name = Form("fHistNTracks%s", chg_name[i_chg].data());
    fHistNTracks[i_chg] =
      new TH3F(hist_name.data(), "Number of ITSsa tracks;Centrality (%);#it{p}_{T} (GeV/#it{c}); Trk Selection",
               nCentBins, centBins, nPtBins, ptBins, nTrkBins, trkBins);
    fHistNTracks[i_chg]->Sumw2();

    TString label("no selection"); // 1
    fHistNTracks[i_chg]->GetZaxis()->SetBinLabel(kHasNoSelection, label.Data());
    label = "ITSsa"; // 2
    fHistNTracks[i_chg]->GetZaxis()->SetBinLabel(kIsITSsa, label.Data());
    label = "ITSrefit"; // 3
    fHistNTracks[i_chg]->GetZaxis()->SetBinLabel(kIsITSrefit, label.Data());
    label = "neutral particle"; // 4
    fHistNTracks[i_chg]->GetZaxis()->SetBinLabel(kIsNotNeutralParticle, label.Data());
    label = "SPDcls"; // 7
    fHistNTracks[i_chg]->GetZaxis()->SetBinLabel(kPassSPD, label.Data());
    label = "SDD+SSD cls"; // 8
    fHistNTracks[i_chg]->GetZaxis()->SetBinLabel(kPassPIDcls, label.Data());
    label = "chi2/ncls"; // 9
    fHistNTracks[i_chg]->GetZaxis()->SetBinLabel(kPassChi2Ncls, label.Data());
    label = "eta"; // 11
    fHistNTracks[i_chg]->GetZaxis()->SetBinLabel(kIsInEta, label.Data());
    label = "dE/dx < 0"; // 12
    fHistNTracks[i_chg]->GetZaxis()->SetBinLabel(kPassdEdx, label.Data());
    label = "Pt cut";
    fHistNTracks[i_chg]->GetZaxis()->SetBinLabel(kPassPtCut, label.Data());
    label = "DCAz"; // 13
    fHistNTracks[i_chg]->GetZaxis()->SetBinLabel(kPassDCAzcut, label.Data());
    label = "DCAxy"; // 14
    fHistNTracks[i_chg]->GetZaxis()->SetBinLabel(kPassDCAxycut, label.Data());
    fOutput->Add(fHistNTracks[i_chg]);
  }

  // binning for the histogram
  const int hnbins = 400;
  double hxmin = 0.01;
  double hxmax = 10;
  double hlogxmin = TMath::Log10(hxmin);
  double hlogxmax = TMath::Log10(hxmax);
  double hbinwidth = (hlogxmax - hlogxmin) / hnbins;
  double hxbins[hnbins + 1];
  hxbins[0] = 0.01;
  for (int i = 1; i <= hnbins; i++) {
    hxbins[i] = hxmin + TMath::Power(10, hlogxmin + i * hbinwidth);
  }

  fHistDEDXGen = new TH2F("fHistDEDXGen", "", hnbins, hxbins, 900, 0, 1000);
  fOutput->Add(fHistDEDXGen);

  fHistDEDX = new TH2F("fHistDEDX", "", hnbins, hxbins, 900, 0, 1000);
  fOutput->Add(fHistDEDX);

  fHistDEDXdouble = new TH2F("fHistDEDXdouble", "", 500, -5, 5, 900, 0, 1000);
  fOutput->Add(fHistDEDXdouble);

  if (fIsMC) { //for correlation between momenta (MC)
    const UInt_t nDimsP = 5;                                         // cent, recP, genP, IsPrim/Sec
    int nBinsP[nDimsP] = { nCentBins, hnbins, hnbins, 4, 900}; //
    double minBinP[nDimsP] = { 0., 0.01, 0.01, -.5, 0.};         // Dummy limits for cent, recP, genP
    double maxBinP[nDimsP] = { 1., 10., 10., 3.5, 1000.};           // Dummy limits for cent, recP, genP
    fHistRecoChargedMC =
      new THnSparseF("fHistRecoChargedMC", ";Centrality (%);#it{p} (GeV/#it{c});#it{p} (GeV/#it{c});", nDimsP,
                     nBinsP, minBinP, maxBinP);
    fHistRecoChargedMC->GetAxis(0)->Set(nCentBins, centBins); // Real limits for cent
    fHistRecoChargedMC->GetAxis(1)->Set(hnbins, hxbins);     // Real limits for rec p
    fHistRecoChargedMC->GetAxis(2)->Set(hnbins, hxbins);     // Real limits for gen p
    fOutput->Add(fHistRecoChargedMC);

    //for efficiency calculation
    fHistMCGenCharged = new TH3F("fHistMCGenCharged", ";Centrality (%);#it{p} (GeV/#it{c});", nCentBins, centBins,
                                  hnbins, hxbins, kNEvtCuts, evBins);
    fOutput->Add(fHistMCGenCharged);
  }

  for (int i_spc = 0; i_spc < kNspc; ++i_spc) {
    for (int i_chg = 0; i_chg < kNchg; ++i_chg) {
      int index = i_spc * kNchg + i_chg;

      const int nDEDXbins = 1000;
      double dedxBins[nDEDXbins + 1];
      SetBins(nDEDXbins, 0., 1000., dedxBins);
      hist_name = Form("fHistNSigmaSep%s%s", spc_name[i_spc].data(), chg_name[i_chg].data());
      fHistNSigmaSep[index] = new TH2F(hist_name.data(), hist_name.data(), hnbins, hxbins, 1000, -10., 10.);
      fOutput->Add(fHistNSigmaSep[index]);
      hist_name = Form("fHistSepPowerReco%s%s", spc_name[i_spc].data(), chg_name[i_chg].data());
      fHistSepPowerReco[index] = new TH2F(hist_name.data(), hist_name.data(), nPtBins, ptBins, nDEDXbins, dedxBins);
      fOutput->Add(fHistSepPowerReco[index]);
      hist_name = Form("fHistSepPowerTrue%s%s", spc_name[i_spc].data(), chg_name[i_chg].data());
      fHistSepPowerTrue[index] = new TH2F(hist_name.data(), hist_name.data(), nPtBins, ptBins, nDEDXbins, dedxBins);
      fOutput->Add(fHistSepPowerTrue[index]);
      // Reconstructed
      hist_name = Form("fHistReco%s%s", spc_name[i_spc].data(), chg_name[i_chg].data());
      fHistReco[index] =
        new TH2F(hist_name.data(), ";Centrality (%);#it{p}_{T} (GeV/#it{c});", nCentBins, centBins, nPtBins, ptBins);
      fOutput->Add(fHistReco[index]);

      hist_name = Form("fHistDCAReco%s%s", spc_name[i_spc].data(), chg_name[i_chg].data());
      fHistDCAReco[index] = new TH3F(hist_name.data(), ";Centrality (%);#it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)",
                                     nCentBins, centBins, nPtBins, ptBins, nDCABins, dcaBins);
      fOutput->Add(fHistDCAReco[index]);

      if (fIsMC) {
        // Histograms MC part Gen bef and afte all selection Good Vertex Gen.
        hist_name = Form("fHistMCPart%s%s", spc_name[i_spc].data(), chg_name[i_chg].data());
        fHistMCPart[index] = new TH3F(hist_name.data(), ";Centrality (%);#it{p}_{T} (GeV/#it{c});", nCentBins, centBins,
                                      nPtBins, ptBins, kNEvtCuts, evBins);
        hist_name = Form("fHistMCPartGoodGenVtxZ%s%s", spc_name[i_spc].data(), chg_name[i_chg].data());
        fHistMCPartGoodGenVtxZ[index] = new TH3F(hist_name.data(), ";Centrality (%);#it{p}_{T} (GeV/#it{c});",
                                                 nCentBins, centBins, nPtBins, ptBins, kNEvtCuts, evBins);
        fOutput->Add(fHistMCPart[index]);
        fOutput->Add(fHistMCPartGoodGenVtxZ[index]);

        hist_name = Form("fHistRecoMC%s%s", spc_name[i_spc].data(), chg_name[i_chg].data());
        const UInt_t nDims = 6;                                         // cent, recPt, genPt, truePID, IsPrim
        int nBins[nDims] = { nCentBins, nPtBins, nPtBins, 1000, 4, 2 }; //
        double minBin[nDims] = { 0., 0., 0., -200, -1.5, -.5 };         // Dummy limits for cent, recPt, genPt
        double maxBin[nDims] = { 1., 1., 1., 200, 2.5, 1.5 };           // Dummy limits for cent, recPt, genPt
        fHistRecoMC[index] =
          new THnSparseF(hist_name.data(), ";Centrality (%);#it{p}_{T} (GeV/#it{c});#it{p}_{T} (GeV/#it{c});", nDims,
                         nBins, minBin, maxBin);
        fHistRecoMC[index]->GetAxis(0)->Set(nCentBins, centBins); // Real limits for cent
        fHistRecoMC[index]->GetAxis(1)->Set(nPtBins, ptBins);     // Real limits for rec pt
        fHistRecoMC[index]->GetAxis(2)->Set(nPtBins, ptBins);     // Real limits for gen pt
        fOutput->Add(fHistRecoMC[index]);

        hist_name = Form("fHistRecoTrueMC%s%s", spc_name[i_spc].data(), chg_name[i_chg].data());
        const UInt_t nDims2 = 4;                                         // cent, recPt, genPt, IsPrim/Sec
        int nBins2[nDims2] = { nCentBins, nPtBins, nPtBins, 4 }; //
        double minBin2[nDims2] = { 0., 0., 0., -.5 };         // Dummy limits for cent, recPt, genPt
        double maxBin2[nDims2] = { 1., 1., 1., 3.5 };           // Dummy limits for cent, recPt, genPt
        fHistRecoTrueMC[index] =
          new THnSparseF(hist_name.data(), ";Centrality (%);#it{p}_{T} (GeV/#it{c});#it{p}_{T} (GeV/#it{c});", nDims2,
                         nBins2, minBin2, maxBin2);
        fHistRecoTrueMC[index]->GetAxis(0)->Set(nCentBins, centBins); // Real limits for cent
        fHistRecoTrueMC[index]->GetAxis(1)->Set(nPtBins, ptBins);     // Real limits for rec pt
        fHistRecoTrueMC[index]->GetAxis(2)->Set(nPtBins, ptBins);     // Real limits for gen pt
        fOutput->Add(fHistRecoTrueMC[index]);

        //For DCAxy
        hist_name = Form("fHistMCDCA%s%s", spc_name[i_spc].data(), chg_name[i_chg].data());
        int nBinsDCA[nDims] = { nCentBins, nPtBins, nPtBins, 4, 4, nDCABins }; //
        double minBinDCA[nDims] = { 0., 0., 0., -1.5, -0.5, -2. };         // Dummy limits for cent, recPt, genPt
        double maxBinDCA[nDims] = { 1., 1., 1., 2.5, 3.5, 2. };           // Dummy limits for cent, recPt, genPt
        fHistMCDCA[index] =
          new THnSparseF(hist_name.data(), ";Centrality (%);#it{p}_{T} (GeV/#it{c});#it{p}_{T} (GeV/#it{c});;;DCAxy", nDims,
                         nBinsDCA, minBinDCA, maxBinDCA);
        fHistMCDCA[index]->GetAxis(0)->Set(nCentBins, centBins); // Real limits for cent
        fHistMCDCA[index]->GetAxis(1)->Set(nPtBins, ptBins);     // Real limits for rec pt
        fHistMCDCA[index]->GetAxis(2)->Set(nPtBins, ptBins);     // Real limits for gen pt
        fOutput->Add(fHistMCDCA[index]);

        //        // Histograms MC part Rec.
        const int nPhysBins = 4;
        double physBins[nPhysBins + 1] = {-0.5, 0.5, 1.5, 2.5, 3.5};
        hist_name = Form("fHistTruePIDMCReco%s%s", spc_name[i_spc].data(), chg_name[i_chg].data());
        fHistTruePIDMCReco[index] = new TH3F(hist_name.data(), ";Centrality (%);#it{p}_{T} (GeV/#it{c});", nCentBins,
                                             centBins, nPtBins, ptBins, nPhysBins, physBins);
        fOutput->Add(fHistTruePIDMCReco[index]);

        hist_name = Form("fHistTruePIDMCGen%s%s", spc_name[i_spc].data(), chg_name[i_chg].data());
        fHistTruePIDMCGen[index] = new TH3F(hist_name.data(), ";Centrality (%);#it{p}_{T} (GeV/#it{c});", nCentBins,
                                             centBins, nPtBins, ptBins, nPhysBins, physBins);
        fOutput->Add(fHistTruePIDMCGen[index]);

        // Histograms MC DCAxy
        hist_name = Form("fHistDCARecoPID_prim%s%s", spc_name[i_spc].data(), chg_name[i_chg].data());
        fHistDCARecoPID_prim[index] =
          new TH3F(hist_name.data(), ";Centrality (%);#it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", nCentBins, centBins,
                   nPtBins, ptBins, nDCABins, dcaBins);

        hist_name = Form("fHistDCARecoPID_sstr%s%s", spc_name[i_spc].data(), chg_name[i_chg].data());
        fHistDCARecoPID_sstr[index] =
          new TH3F(hist_name.data(), ";Centrality (%);#it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", nCentBins, centBins,
                   nPtBins, ptBins, nDCABins, dcaBins);

        hist_name = Form("fHistDCARecoPID_smat%s%s", spc_name[i_spc].data(), chg_name[i_chg].data());
        fHistDCARecoPID_smat[index] =
          new TH3F(hist_name.data(), ";Centrality (%);#it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", nCentBins, centBins,
                   nPtBins, ptBins, nDCABins, dcaBins);

        hist_name = Form("fHistDCATruePID_prim%s%s", spc_name[i_spc].data(), chg_name[i_chg].data());
        fHistDCATruePID_prim[index] =
          new TH3F(hist_name.data(), ";Centrality (%);#it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", nCentBins, centBins,
                   nPtBins, ptBins, nDCABins, dcaBins);

        hist_name = Form("fHistDCATruePID_sstr%s%s", spc_name[i_spc].data(), chg_name[i_chg].data());
        fHistDCATruePID_sstr[index] =
          new TH3F(hist_name.data(), ";Centrality (%);#it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", nCentBins, centBins,
                   nPtBins, ptBins, nDCABins, dcaBins);

        hist_name = Form("fHistDCATruePID_smat%s%s", spc_name[i_spc].data(), chg_name[i_chg].data());
        fHistDCATruePID_smat[index] =
          new TH3F(hist_name.data(), ";Centrality (%);#it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", nCentBins, centBins,
                   nPtBins, ptBins, nDCABins, dcaBins);

        fOutput->Add(fHistDCARecoPID_prim[index]);
        fOutput->Add(fHistDCARecoPID_sstr[index]);
        fOutput->Add(fHistDCARecoPID_smat[index]);
        fOutput->Add(fHistDCATruePID_prim[index]);
        fOutput->Add(fHistDCATruePID_sstr[index]);
        fOutput->Add(fHistDCATruePID_smat[index]);
      } // end IsMC
    }
  }

  for (int i = 0; i < 4; i++) {
    fHistCharge[i] = new TH1F(Form("fHistChargeLay%d", i), Form("fHistChargeLay%d", i), 100, 0, 300);
    fOutput->Add(fHistCharge[i]);
  }

  if (fFillIntDistHist) {

    // dEdx distributions

    fHistPosHypPi =
      new TH2F("fHistPosHypPi", "fHistPosHypPi", nPtBins, ptBins, nDCABins, dcaBins); // DCA distr. with NSigma PID
    fHistPosHypKa = new TH2F("fHistPosHypKa", "fHistPosHypKa", nPtBins, ptBins, nDCABins, dcaBins);
    fHistPosHypPr = new TH2F("fHistPosHypPr", "fHistPosHypPr", nPtBins, ptBins, nDCABins, dcaBins);
    fHistNegHypPi = new TH2F("fHistNegHypPi", "fHistNegHypPi", nPtBins, ptBins, nDCABins, dcaBins);
    fHistNegHypKa = new TH2F("fHistNegHypKa", "fHistNegHypKa", nPtBins, ptBins, nDCABins, dcaBins);
    fHistNegHypPr = new TH2F("fHistNegHypPr", "fHistNegHypPr", nPtBins, ptBins, nDCABins, dcaBins);

    fOutput->Add(fHistPosHypPi); // DCA distr
    fOutput->Add(fHistPosHypKa);
    fOutput->Add(fHistPosHypPr);
    fOutput->Add(fHistNegHypPi);
    fOutput->Add(fHistNegHypKa);
    fOutput->Add(fHistNegHypPr);

    if (fIsMC) {
      const int nBins = 175;
      double dEdxBins[nBins + 1];
      SetBins(nBins, -3.5, 3.5, dEdxBins);
      fHistMCPosOtherHypPion =
        new TH2F("fHistMCPosOtherHypPion", "fHistMCPosOtherHypPion", nPtBins, ptBins, nBins, dEdxBins); // MC truth
      fHistMCPosOtherHypKaon =
        new TH2F("fHistMCPosOtherHypKaon", "fHistMCPosOtherHypKaon", nPtBins, ptBins, nBins, dEdxBins);
      fHistMCPosOtherHypProt =
        new TH2F("fHistMCPosOtherHypProt", "fHistMCPosOtherHypProt", nPtBins, ptBins, nBins, dEdxBins);
      fHistMCPosElHypPion = new TH2F("fHistMCPosElHypPion", "fHistMCPosElHypPion", nPtBins, ptBins, nBins, dEdxBins);
      fHistMCPosElHypKaon = new TH2F("fHistMCPosElHypKaon", "fHistMCPosElHypKaon", nPtBins, ptBins, nBins, dEdxBins);
      fHistMCPosElHypProt = new TH2F("fHistMCPosElHypProt", "fHistMCPosElHypProt", nPtBins, ptBins, nBins, dEdxBins);
      fHistMCPosPiHypPion = new TH2F("fHistMCPosPiHypPion", "fHistMCPosPiHypPion", nPtBins, ptBins, nBins, dEdxBins);
      fHistMCPosPiHypKaon = new TH2F("fHistMCPosPiHypKaon", "fHistMCPosPiHypKaon", nPtBins, ptBins, nBins, dEdxBins);
      fHistMCPosPiHypProt =
        new TH2F("fHistMCPosPiHypProton", "fHistMCPosPiHypProton", nPtBins, ptBins, nBins, dEdxBins);
      fHistMCPosKaHypPion = new TH2F("fHistMCPosKaHypPion", "fHistMCPosKaHypPion", nPtBins, ptBins, nBins, dEdxBins);
      fHistMCPosKaHypKaon = new TH2F("fHistMCPosKaHypKaon", "fHistMCPosKaHypKaon", nPtBins, ptBins, nBins, dEdxBins);
      fHistMCPosKaHypProt = new TH2F("fHistMCPosKaHypProt", "fHistMCPosKaHypProt", nPtBins, ptBins, nBins, dEdxBins);
      fHistMCPosPrHypPion = new TH2F("fHistMCPosPrHypPion", "fHistMCPosPrHypPion", nPtBins, ptBins, nBins, dEdxBins);
      fHistMCPosPrHypKaon = new TH2F("fHistMCPosPrHypKaon", "fHistMCPosPrHypKaon", nPtBins, ptBins, nBins, dEdxBins);
      fHistMCPosPrHypProt = new TH2F("fHistMCPosPrHypProt", "fHistMCPosPrHypProt", nPtBins, ptBins, nBins, dEdxBins);

      fHistMCNegOtherHypPion =
        new TH2F("fHistMCNegOtherHypPion", "fHistMCNegOtherHypPion", nPtBins, ptBins, nBins, dEdxBins); // MC truth
      fHistMCNegOtherHypKaon =
        new TH2F("fHistMCNegOtherHypKaon", "fHistMCNegOtherHypKaon", nPtBins, ptBins, nBins, dEdxBins);
      fHistMCNegOtherHypProt =
        new TH2F("fHistMCNegOtherHypProt", "fHistMCNegOtherHypProt", nPtBins, ptBins, nBins, dEdxBins);
      fHistMCNegElHypPion = new TH2F("fHistMCNegElHypPion", "fHistMCNegElHypPion", nPtBins, ptBins, nBins, dEdxBins);
      fHistMCNegElHypKaon = new TH2F("fHistMCNegElHypKaon", "fHistMCNegElHypKaon", nPtBins, ptBins, nBins, dEdxBins);
      fHistMCNegElHypProt = new TH2F("fHistMCNegElHypProt", "fHistMCNegElHypProt", nPtBins, ptBins, nBins, dEdxBins);
      fHistMCNegPiHypPion = new TH2F("fHistMCNegPiHypPion", "fHistMCNegPiHypPion", nPtBins, ptBins, nBins, dEdxBins);
      fHistMCNegPiHypKaon = new TH2F("fHistMCNegPiHypKaon", "fHistMCNegPiHypKaon", nPtBins, ptBins, nBins, dEdxBins);
      fHistMCNegPiHypProt = new TH2F("fHistMCNegPiHypProt", "fHistMCNegPiHypProt", nPtBins, ptBins, nBins, dEdxBins);
      fHistMCNegKaHypPion = new TH2F("fHistMCNegKaHypPion", "fHistMCNegKaHypPion", nPtBins, ptBins, nBins, dEdxBins);
      fHistMCNegKaHypKaon = new TH2F("fHistMCNegKaHypKaon", "fHistMCNegKaHypKaon", nPtBins, ptBins, nBins, dEdxBins);
      fHistMCNegKaHypProt = new TH2F("fHistMCNegKaHypProt", "fHistMCNegKaHypProt", nPtBins, ptBins, nBins, dEdxBins);
      fHistMCNegPrHypPion = new TH2F("fHistMCNegPrHypPion", "fHistMCNegPrHypPion", nPtBins, ptBins, nBins, dEdxBins);
      fHistMCNegPrHypKaon = new TH2F("fHistMCNegPrHypKaon", "fHistMCNegPrHypKaon", nPtBins, ptBins, nBins, dEdxBins);
      fHistMCNegPrHypProt = new TH2F("fHistMCNegPrHypProt", "fHistMCNegPrHypProt", nPtBins, ptBins, nBins, dEdxBins);

      fOutput->Add(fHistMCPosOtherHypPion); // MC truth
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

      fOutput->Add(fHistMCNegOtherHypPion); // MC truth
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
    } // end IsMC
  }   // end FillIntDistHist

  PostData(1, fOutput);

  CreateDCAcutFunctions(); // Creating kParamContainer data
  PostData(2, fListCuts);

  // Post output data container
  if (fFillNtuple) {
    fListTree = new TList();
    fListTree->SetOwner();

    fNtupleData = new TNtuple("fNtupleData", "fNtupleData", "mult:p:pt:s0:s1:s2:s3:dEdx:sign:eta:dcaXY:dcaZ:clumap:MCisph:MCdpg:MCpt");
    fListTree->Add(fNtupleData);
    fNtupleMC = new TNtuple("fNtupleMC", "fNtupleMC", "mult:mcPt:pdgcode:sign:mcEta:mcRap:isph:run");
    fListTree->Add(fNtupleMC);
    PostData(3, fListTree);
  }

  AliInfo("End of CreateOutputObjects");
}

//
//
//___________________________________________________________________________
void AliAnalysisTaskSEITSsaSpectra::CreateDCAcutFunctions()
{
  fListCuts = new TList();
  fListCuts->SetOwner();
  double xyP[3];
  double zP[3];
  if (fYear == 2009) {
    if (fIsMC) {
      xyP[0] = 88.63; // MC LHC10a12
      xyP[1] = 19.57;
      xyP[2] = 1.65;
      zP[0] = 140.98;
      zP[1] = 62.33;
      zP[2] = 1.15;
    } else {
      xyP[0] = 85.28; // DATA 900 GeV pass6
      xyP[1] = 25.78;
      xyP[2] = 1.55;
      zP[0] = 146.80;
      zP[1] = 70.07;
      zP[2] = 1.11;
    }
  } else if (fYear == 2010) {
    if (fIsMC) {
      xyP[0] = 36.; // MC LHC10d1
      xyP[1] = 43.9;
      xyP[2] = 1.3;
      zP[0] = 111.9;
      zP[1] = 59.8;
      zP[2] = 1.2;
    } else {
      xyP[0] = 32.7; // DATA 7 TeV pass2
      xyP[1] = 44.8;
      xyP[2] = 1.3;
      zP[0] = 117.3;
      zP[1] = 66.8;
      zP[2] = 1.2;
    }
  }
  fDCAxyCutFunc = new TF1("fDCAxyCutFunc", "[3]*([0]+[1]/TMath::Power(TMath::Abs(x),[2]))", .05, 10.);
  for (int ipar = 0; ipar < 3; ipar++)
    fDCAxyCutFunc->SetParameter(ipar, xyP[ipar]);
  fDCAxyCutFunc->SetParameter(3, fNSigmaDCAxy);
  fDCAxyCutFunc->SetParName(3, "Sigmas");

  fDCAzCutFunc = new TF1("fDCAzCutFunc", "[3]*([0]+[1]/TMath::Power(TMath::Abs(x),[2]))", .05, 10.);
  for (int ipar = 0; ipar < 3; ipar++)
    fDCAzCutFunc->SetParameter(ipar, zP[ipar]);
  fDCAzCutFunc->SetParameter(3, fNSigmaDCAz);
  fDCAzCutFunc->SetParName(3, "Sigmas");

  fListCuts->Add(fDCAxyCutFunc);
  fListCuts->Add(fDCAzCutFunc);
}

//
//
//________________________________________________________________________
void AliAnalysisTaskSEITSsaSpectra::Initialization()
{
  // Initialization
  Printf("Inizializing Task, be sure to run after all configuration have been set...");
  AliInfo("Tracks selections");
  AliInfoF(
    " y = yLab + %.3f,  Ymin %.1f, Ymax %.1f, Eabs %.1f, DCAxyCut %.1f, DCAzCut %.1f, Chi2 %.1f,   nSPD %d,   nPID %d",
    fCMSRapFct, fMinRapCut, fMaxRapCut, fAbsEtaCut, fNSigmaDCAxy, fNSigmaDCAz, fMaxChi2Clu, fMinSPDPts,
    fMinNdEdxSamples);

  if (fMultMethod)
    AliInfoF("Cent. %.f %.f %s", fLowMult, fUpMult, fMultEstimator.Data());

  return;
}

//
//
//________________________________________________________________________
void AliAnalysisTaskSEITSsaSpectra::UserExec(Option_t *)
{
  // Main loop
  // Called for each event
  //
  fESD = dynamic_cast<AliESDEvent *>(InputEvent());
  if (!fESD) {
    AliDebug(3, "fESD not available");
    PostAllData();
    return;
  }

  // Fill some histograms before event selection
  // Check Monte Carlo information and other access first:
  // AliStack*   lMCstack = NULL;
  AliMCEvent *lMCevent = NULL;
  bool lHasGoodVtxGen = kFALSE;
  if (fIsMC) {
    AliMCEventHandler *lMCevtHandler =
      dynamic_cast<AliMCEventHandler *>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
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
    if (!lMCevent->GetPrimaryVertex()) {
      AliDebug(3, "MC Vtx not available");
      PostAllData();
      return;
    }
    // Selection of the MC sample
    if (!(TMath::Abs(lMCevent->GetPrimaryVertex()->GetZ()) > fMaxVtxZCut)) {
      lHasGoodVtxGen = kTRUE;
    }
  }

  // Event selection
  EEvtCut_Type lastEvtCutPassed = kIsReadable;
  bool lIsEventSelected = IsEventAccepted(lastEvtCutPassed);
  for (int istep = (int)kIsReadable; istep <= (int)lastEvtCutPassed; ++istep) {
    fHistNEvents->Fill(fEvtMult, istep);
    if (lHasGoodVtxGen)
      fHistMCEvents->Fill(fEvtMult, istep);
  }
  if (fIsMC)
    AnalyseMCParticles(lMCevent, lastEvtCutPassed, lHasGoodVtxGen);

  if (!lIsEventSelected) {
    AliDebugF(3, "Event rejected by event selection after %d step", (int)lastEvtCutPassed);
    PostAllData();
    return;
  }

  fHistNEvents->Fill(fEvtMult, kNEvtCuts);
  if (lHasGoodVtxGen)
    fHistMCEvents->Fill(fEvtMult, kNEvtCuts);

  if (fMultMethod) // Fill fHistMultAftEvtSel after the event Selection
    fHistMultAftEvtSel->Fill(fEvtMult);

  if (!fITSPIDResponse)
    fITSPIDResponse = new AliITSPIDResponse(fIsMC);

  // loop on tracks
  for (int i_trk = 0; i_trk < fESD->GetNumberOfTracks(); ++i_trk) {
    AliESDtrack *track = dynamic_cast<AliESDtrack *>(fESD->GetTrack(i_trk));
    if (!track)
      continue;

    if (fIsMC) {
      int trkLabel = TMath::Abs(track->GetLabel());

      AliMCParticle *mcTrk = ((AliMCParticle *)lMCevent->GetTrack(trkLabel));
      int pdg = mcTrk->PdgCode();
      if (TMath::Abs(pdg) > 1E10) // protection to remove High ionization part
        continue;
    }

    // track selection
    double trkPt = track->Pt();
    double dEdxLay[4];
    track->GetITSdEdxSamples(dEdxLay);
    double dEdx = track->GetITSsignal();

    int i_chg = (track->GetSign() > 0) ? 0 : 1; // 0: Pos; 1:Neg
    ULong_t status = track->GetStatus();
    UInt_t clumap = track->GetITSClusterMap();

    int nSPD = 0;
    int nPtsForPid = 0;
    for (int il = 0; il < 6; il++) {
      if (TESTBIT(clumap, il)) {
        if (il < 2)
          nSPD++;
        else
          nPtsForPid++;
      }
    }

    ETrkCut_Type trkSel = kHasNoSelection;

    //"no selection"
    fHistDEDXGen->Fill(track->GetP(), dEdx);
    fHistNTracks[i_chg]->Fill(fEvtMult, trkPt, trkSel);

    //"ITSsa"
    if (!(status & AliESDtrack::kITSpureSA))
      continue;
    trkSel = kIsITSsa;
    fHistNTracks[i_chg]->Fill(fEvtMult, trkPt, trkSel);

    //"ITSrefit"
    if (!(status & AliESDtrack::kITSrefit))
      continue;
    trkSel = kIsITSrefit;
    fHistNTracks[i_chg]->Fill(fEvtMult, trkPt, trkSel);

    //"neutral particle"
    if (TMath::Abs(track->GetSign()) < 0.0001)
      continue;
    trkSel = kIsNotNeutralParticle;
    fHistNTracks[i_chg]->Fill(fEvtMult, trkPt, trkSel);

    //"SPDcls"
    if (nSPD < fMinSPDPts)
      continue; // At least one point in the SPD
    trkSel = kPassSPD;
    fHistNTracks[i_chg]->Fill(fEvtMult, trkPt, trkSel);

    //"SDD+SSD cls" at least 3 points on SSD/SDD
    if (nPtsForPid < fMinNdEdxSamples)
      continue;
    trkSel = kPassPIDcls;
    fHistNTracks[i_chg]->Fill(fEvtMult, trkPt, trkSel);

    //"chi2/ncls"->chisquare/nclusters
    int nclu = nSPD + nPtsForPid;
    if (track->GetITSchi2() / nclu > fMaxChi2Clu)
      continue;
    trkSel = kPassChi2Ncls;
    fHistNTracks[i_chg]->Fill(fEvtMult, trkPt, trkSel);

    //"eta"->pseudorapidity
    if (TMath::Abs(track->Eta()) > fAbsEtaCut)
      continue;
    trkSel = kIsInEta;
    fHistNTracks[i_chg]->Fill(fEvtMult, trkPt, trkSel);

    //"dE/dx < 0"->truncated mean
    if (dEdx < 0)
      continue;
    trkSel = kPassdEdx;
    fHistNTracks[i_chg]->Fill(fEvtMult, trkPt, trkSel);

    // fill propaganda plot with dedx before pt cut
    fHistDEDX->Fill(track->GetP(), dEdx);
    fHistDEDXdouble->Fill(track->GetP() * track->GetSign(), dEdx);

    if(fIsMC){//correlation between momenta (measured and true ones) --> before pt cut!
      int lMCtrk = TMath::Abs(track->GetLabel());
      AliMCParticle *trkMC = (AliMCParticle *)lMCevent->GetTrack(lMCtrk);
      float pMC   = trkMC->P();
      int ptype = 0;
      if (lMCevent->IsPhysicalPrimary(lMCtrk)){
        ptype = 0;
      }
      else if (lMCevent->IsSecondaryFromWeakDecay(lMCtrk)){
        ptype = 1;
      }
      else if (lMCevent->IsSecondaryFromMaterial(lMCtrk)){
        ptype = 2;
      }
      else {
        ptype = 3;
      }
      double tmp_vect[5] = {fEvtMult, track->GetP(), pMC, static_cast<double>(ptype), dEdx};
      fHistRecoChargedMC->Fill(tmp_vect);
    }

    //"ptCut"
    if ((trkPt < fPtBins[0]) || (trkPt >= fPtBins[fPtBins.GetSize() - 1]))
      continue;
    trkSel = kPassPtCut;
    fHistNTracks[i_chg]->Fill(fEvtMult, trkPt, trkSel);

    /////////////////////////////////////////////////////////////////////////////
    // Because we fit the DCA distributions, we need them after the DCAz cut,  //
    // otherwise the templates and distributions are modified by this cut      //
    // (if it is done with the one of the DCAxy)                               //
    /////////////////////////////////////////////////////////////////////////////
    float impactXY, impactZ;
    track->GetImpactParameters(impactXY, impactZ);
    //"DCAz"
    if (!DCAcutZ(impactZ, trkPt))
      continue;
    trkSel = kPassDCAzcut;
    fHistNTracks[i_chg]->Fill(fEvtMult, trkPt, trkSel);

    if (fFillNtuple) {
      float xnt[16];
      int index = 0;
      /*1 */ xnt[index++] = (float)fEvtMult;
      /*2 */ xnt[index++] = (float)track->GetP();
      /*3 */ xnt[index++] = (float)track->Pt();
      /*4 */ xnt[index++] = (float)dEdxLay[0];
      /*5 */ xnt[index++] = (float)dEdxLay[1];
      /*6 */ xnt[index++] = (float)dEdxLay[2];
      /*7 */ xnt[index++] = (float)dEdxLay[3];
      /*8 */ xnt[index++] = (float)dEdx;
      /*9*/ xnt[index++] = (float)track->GetSign();
      /*10*/ xnt[index++] = (float)track->Eta();
      /*11*/ xnt[index++] = (float)impactXY;
      /*12*/ xnt[index++] = (float)impactZ;
      /*13*/ xnt[index++] = (float)clumap;

      float lMCpt = -999;
      int lMCpdg = -999;
      int lMCisph = -999;
      if (fIsMC) {
        int lMCtrk = TMath::Abs(track->GetLabel());
        lMCisph=lMCevent->IsPhysicalPrimary(lMCtrk);

        AliMCParticle *trkMC = (AliMCParticle *)lMCevent->GetTrack(lMCtrk);
        lMCpdg = trkMC->PdgCode();
        lMCpt =  trkMC->Pt();
      }
      /*14*/ xnt[index++] = (float)lMCisph;
      /*15*/ xnt[index++] = (float)lMCpdg;
      /*16*/ xnt[index++] = (float)lMCpt;
      fNtupleData->Fill(xnt);
    } else {
      // track PID aproach
      double logdiff[4];
      int fPid = GetTrackPid(track, logdiff);
      // End PID approach

      // Compute y
      double y[AliPID::kSPECIES];
      // loop per specie x4
      for (int i_spc = 0; i_spc < AliPID::kSPECIES; ++i_spc) {
        float mass = AliPID::ParticleMass(AliPID::EParticleType(i_spc));
        y[i_spc] = Eta2y(trkPt, mass, track->Eta());
        y[i_spc] += fCMSRapFct;
      }
      bool lIsGoodTrack = ((fPid > AliPID::kMuon) && IsRapIn(y[fPid]));
      int lPidIndex = lIsGoodTrack ? ((fPid - 2) * kNchg + i_chg) : -1;

      int lMCtrk = -999;
      int lMCpdg = -999;
      int lMCspc = AliPID::kElectron;
      float lMCpt = -999;
      float lMCp = -999;
      if (fIsMC) {
        lMCtrk = TMath::Abs(track->GetLabel());
        AliMCParticle *trkMC = (AliMCParticle *)lMCevent->GetTrack(lMCtrk);
        lMCpdg = trkMC->PdgCode();
        lMCpt  = trkMC->Pt();
        lMCp   = trkMC->P();

        //        if (TMath::Abs(lMCpdg) ==   11 && fPid == AliPID::kPion) lMCspc = AliPID::kPion;
        //        if (TMath::Abs(lMCpdg) ==   13 && fPid == AliPID::kPion) lMCspc = AliPID::kPion;
        if (TMath::Abs(lMCpdg) == 211)
          lMCspc = AliPID::kPion; // select Pi+/Pi- only
        if (TMath::Abs(lMCpdg) == 321)
          lMCspc = AliPID::kKaon; // select K+/K- only
        if (TMath::Abs(lMCpdg) == 2212)
          lMCspc = AliPID::kProton; // select p+/p- only
        // otherwise considered as electron and skkipped  later
      }
      bool lIsGoodPart = ((lMCspc > AliPID::kMuon) && IsRapIn(y[lMCspc]));
      int lMCtIndex = lIsGoodPart ? ((lMCspc - 2) * kNchg + i_chg) : -1;

      if (lIsGoodTrack) {
        // DCA distributions, before the DCA cuts, based on PID approach
        fHistDCAReco[lPidIndex]->Fill(fEvtMult, trkPt, impactXY);

        // DCA distributions, before the DCAxy cuts from the MC kinematics
        // Filling DCA distribution with MC truth Physics values
        if (fIsMC && fIsDCAUnfoldHistoEnabled) {
          int ptype = 0;
          if (lMCevent->IsPhysicalPrimary(lMCtrk)){
            fHistDCARecoPID_prim[lPidIndex]->Fill(fEvtMult, trkPt, impactXY);
            ptype = 0;
          }
          else if (lMCevent->IsSecondaryFromWeakDecay(lMCtrk)){
            fHistDCARecoPID_sstr[lPidIndex]->Fill(fEvtMult, trkPt, impactXY);
            ptype = 1;
          }
          else if (lMCevent->IsSecondaryFromMaterial(lMCtrk)){
            fHistDCARecoPID_smat[lPidIndex]->Fill(fEvtMult, trkPt, impactXY);
            ptype = 2;
          }
          else {
            ptype = 3;
          }

          int binPart = (lMCspc > AliPID::kMuon) ? (lMCspc - 2) : -1;
          double tmp_vect[6] = { fEvtMult,
                                 trkPt,
                                 lMCpt,
                                 static_cast<double>(binPart),
                                 static_cast<double>(ptype),
                                 static_cast<double>(impactXY)
                               };
          fHistMCDCA[lPidIndex]->Fill(tmp_vect);
        }
      } // end lIsGoodTrack

      if (fIsMC && lIsGoodPart) {
        // DCA distributions, before the DCAxy cuts from the MC kinematics
        // Filling DCA distribution with MC truth Physics values
        if (lMCevent->IsPhysicalPrimary(lMCtrk))
          fHistDCATruePID_prim[lMCtIndex]->Fill(fEvtMult, trkPt, impactXY);
        if (lMCevent->IsSecondaryFromWeakDecay(lMCtrk))
          fHistDCATruePID_sstr[lMCtIndex]->Fill(fEvtMult, trkPt, impactXY);
        if (lMCevent->IsSecondaryFromMaterial(lMCtrk))
          fHistDCATruePID_smat[lMCtIndex]->Fill(fEvtMult, trkPt, impactXY);
      } // end lIsGoodPart
      //"DCAxy"
      if (!DCAcutXY(impactXY, trkPt))
        continue;
      trkSel = kPassDCAxycut;
      fHistNTracks[i_chg]->Fill(fEvtMult, trkPt, trkSel);

      if (lIsGoodTrack)
        fHistReco[lPidIndex]->Fill(fEvtMult, trkPt);
      fHistSepPowerReco[lPidIndex]->Fill(trkPt, dEdx);
      // Filling Histos for Reco Efficiency
      // information from the MC kinematics (truth PID)
      int ptype = 0;
      if (fIsMC && lIsGoodPart) {
        fHistSepPowerTrue[lMCtIndex]->Fill(lMCpt, dEdx);
        if (lMCevent->IsPhysicalPrimary(lMCtrk)){
          fHistTruePIDMCReco[lMCtIndex]->Fill(fEvtMult, trkPt, 0);
          fHistTruePIDMCGen[lMCtIndex]->Fill(fEvtMult, lMCpt, 0);
          ptype = 0;
        }
        else if (lMCevent->IsSecondaryFromWeakDecay(lMCtrk)){
          fHistTruePIDMCReco[lMCtIndex]->Fill(fEvtMult, trkPt, 1);
          fHistTruePIDMCGen[lMCtIndex]->Fill(fEvtMult, lMCpt, 1);
          ptype = 1;
        }
        else if (lMCevent->IsSecondaryFromMaterial(lMCtrk)){
          fHistTruePIDMCReco[lMCtIndex]->Fill(fEvtMult, trkPt, 2);
          fHistTruePIDMCGen[lMCtIndex]->Fill(fEvtMult, lMCpt, 2);
          ptype = 2;
        }
        else {
          fHistTruePIDMCReco[lMCtIndex]->Fill(fEvtMult, trkPt, 3);
          fHistTruePIDMCGen[lMCtIndex]->Fill(fEvtMult, lMCpt, 3);
          ptype = 3;
          AliWarning("Weird particle physics");
        }

        double tmp_vect2[4] = {fEvtMult, trkPt, lMCpt, static_cast<double>(ptype)};
        fHistRecoTrueMC[lMCtIndex]->Fill(tmp_vect2);
      }

      int binPart = (lMCspc > AliPID::kMuon) ? (lMCspc - 2) : -1;
      if (fIsMC && lIsGoodTrack) {
        double invPtDifference = 1. / trkPt - 1. / lMCpt;
        double tmp_vect[6] = { fEvtMult,
                               trkPt,
                               lMCpt,
                               invPtDifference,
                               static_cast<double>(binPart),
                               (lMCevent->IsPhysicalPrimary(lMCtrk)) ? 0. : 1. };
        fHistRecoMC[lPidIndex]->Fill(tmp_vect);
      } // end y

      if (lIsGoodTrack && fFillIntDistHist) {
        //
        // integral approach histograms
        //
        int lAbsPdgCode = TMath::Abs(lMCpdg);
        if (track->GetSign() > 0) {
          if (IsRapIn(y[AliPID::kPion]))
            fHistPosHypPi->Fill(trkPt, logdiff[0]);
          if (IsRapIn(y[AliPID::kKaon]))
            fHistPosHypKa->Fill(trkPt, logdiff[1]);
          if (IsRapIn(y[AliPID::kProton]))
            fHistPosHypPr->Fill(trkPt, logdiff[2]);
          if (fIsMC) {
            if (IsRapIn(y[AliPID::kPion])) {
              if ((lAbsPdgCode != 11) && (lAbsPdgCode != 211) && (lAbsPdgCode != 321) && (lAbsPdgCode != 2212))
                fHistMCPosOtherHypPion->Fill(trkPt, logdiff[0]);
              if (lAbsPdgCode == 11)
                fHistMCPosElHypPion->Fill(trkPt, logdiff[0]);
              if (lAbsPdgCode == 211)
                fHistMCPosPiHypPion->Fill(trkPt, logdiff[0]);
              if (lAbsPdgCode == 321)
                fHistMCPosKaHypPion->Fill(trkPt, logdiff[0]);
              if (lAbsPdgCode == 2212)
                fHistMCPosPrHypPion->Fill(trkPt, logdiff[0]);
            }
            if (IsRapIn(y[AliPID::kKaon])) {
              if ((lAbsPdgCode != 11) && (lAbsPdgCode != 211) && (lAbsPdgCode != 321) && (lAbsPdgCode != 2212))
                fHistMCPosOtherHypKaon->Fill(trkPt, logdiff[1]);
              if (lAbsPdgCode == 11)
                fHistMCPosElHypKaon->Fill(trkPt, logdiff[1]);
              if (lAbsPdgCode == 211)
                fHistMCPosPiHypKaon->Fill(trkPt, logdiff[1]);
              if (lAbsPdgCode == 321)
                fHistMCPosKaHypKaon->Fill(trkPt, logdiff[1]);
              if (lAbsPdgCode == 2212)
                fHistMCPosPrHypKaon->Fill(trkPt, logdiff[1]);
            }
            if (IsRapIn(y[AliPID::kProton])) {
              if ((lAbsPdgCode != 11) && (lAbsPdgCode != 211) && (lAbsPdgCode != 321) && (lAbsPdgCode != 2212))
                fHistMCPosOtherHypProt->Fill(trkPt, logdiff[2]);
              if (lAbsPdgCode == 11)
                fHistMCPosElHypProt->Fill(trkPt, logdiff[2]);
              if (lAbsPdgCode == 211)
                fHistMCPosPiHypProt->Fill(trkPt, logdiff[2]);
              if (lAbsPdgCode == 321)
                fHistMCPosKaHypProt->Fill(trkPt, logdiff[2]);
              if (lAbsPdgCode == 2212)
                fHistMCPosPrHypProt->Fill(trkPt, logdiff[2]);
            }
          }
        } else {
          if (IsRapIn(y[AliPID::kPion]))
            fHistNegHypPi->Fill(trkPt, logdiff[0]);
          if (IsRapIn(y[AliPID::kKaon]))
            fHistNegHypKa->Fill(trkPt, logdiff[1]);
          if (IsRapIn(y[AliPID::kProton]))
            fHistNegHypPr->Fill(trkPt, logdiff[2]);
          if (fIsMC) {
            if (IsRapIn(y[AliPID::kPion])) {
              if ((lAbsPdgCode != 11) && (lAbsPdgCode != 211) && (lAbsPdgCode != 321) && (lAbsPdgCode != 2212))
                fHistMCNegOtherHypPion->Fill(trkPt, logdiff[0]);
              if (lAbsPdgCode == 11)
                fHistMCNegElHypPion->Fill(trkPt, logdiff[0]);
              if (lAbsPdgCode == 211)
                fHistMCNegPiHypPion->Fill(trkPt, logdiff[0]);
              if (lAbsPdgCode == 321)
                fHistMCNegKaHypPion->Fill(trkPt, logdiff[0]);
              if (lAbsPdgCode == 2212)
                fHistMCNegPrHypPion->Fill(trkPt, logdiff[0]);
            }
            if (IsRapIn(y[AliPID::kKaon])) {
              if ((lAbsPdgCode != 11) && (lAbsPdgCode != 211) && (lAbsPdgCode != 321) && (lAbsPdgCode != 2212))
                fHistMCNegOtherHypKaon->Fill(trkPt, logdiff[1]);
              if (lAbsPdgCode == 11)
                fHistMCNegElHypKaon->Fill(trkPt, logdiff[1]);
              if (lAbsPdgCode == 211)
                fHistMCNegPiHypKaon->Fill(trkPt, logdiff[1]);
              if (lAbsPdgCode == 321)
                fHistMCNegKaHypKaon->Fill(trkPt, logdiff[1]);
              if (lAbsPdgCode == 2212)
                fHistMCNegPrHypKaon->Fill(trkPt, logdiff[1]);
            }
            if (IsRapIn(y[AliPID::kProton])) {
              if ((lAbsPdgCode != 11) && (lAbsPdgCode != 211) && (lAbsPdgCode != 321) && (lAbsPdgCode != 2212))
                fHistMCNegOtherHypProt->Fill(trkPt, logdiff[2]);
              if (lAbsPdgCode == 11)
                fHistMCNegElHypProt->Fill(trkPt, logdiff[2]);
              if (lAbsPdgCode == 211)
                fHistMCNegPiHypProt->Fill(trkPt, logdiff[2]);
              if (lAbsPdgCode == 321)
                fHistMCNegKaHypProt->Fill(trkPt, logdiff[2]);
              if (lAbsPdgCode == 2212)
                fHistMCNegPrHypProt->Fill(trkPt, logdiff[2]);
            }
          }
        }
      } // end lIsGooTrack
    }   // end else fill Ntuple
  }     // end track loop

  // Post output data.
  PostAllData();
  AliDebug(4, "............. end of Exec");
}

//
//
//________________________________________________________________________
void AliAnalysisTaskSEITSsaSpectra::Terminate(Option_t *)
{
  // Merge output
  // Called once at the end of the query

  AliInfo("end of Terminate");
  return;
}

//
//
//________________________________________________________________________
bool AliAnalysisTaskSEITSsaSpectra::IsEventAccepted(EEvtCut_Type &evtSel)
{
  evtSel = kIsReadable;
  // Check multiplicity selection first
  if (!IsMultSelected()) {
    AliDebug(3, "Event doesn't pass multiplicity selection");
    return kFALSE;
  } else {
    fHistMultBefEvtSel->Fill(fEvtMult);
  }
  evtSel = kPassMultSel;

  // Check if has SDD info (if requiered)
  if (fChkIsSDDIn && !fIsMC) {
    TString firedTriggerClasses(fESD->GetFiredTriggerClasses());
    if (!(firedTriggerClasses.Contains("ALL") || firedTriggerClasses.Contains("CENT"))) {
      AliDebug(3, "Event dont accepted by AliEventCuts");
      AliDebug(3, "Event Rejected: SDD out trigger cluster");
      return kFALSE;
    }
  }
  evtSel = kIsSDDIn;

  if (fExtEventCuts) {
    if (fEventCuts.AcceptEvent(fESD)) {
      evtSel = static_cast<EEvtCut_Type>(kNEvtCuts - 1);
      return kTRUE;
    }
    AliDebug(3, "Event dont accepted by AliEventCuts");

    if (!fEventCuts.PassedCut(AliEventCuts::kDAQincomplete)) {
      AliDebug(3, "Event with incomplete DAQ");
      return kFALSE;
    }
    evtSel = kIsNotIncDAQ;

    if (!fEventCuts.PassedCut(AliEventCuts::kTrigger)) {
      AliDebug(3, "Event doesn't pass physics evt. sel. for trigger");
      return kFALSE;
    }
    evtSel = kPassTrig;

    if (!fEventCuts.PassedCut(AliEventCuts::kINELgt0)) {
      AliDebug(3, "Event doesn't pass INEL>0 cut");
      return kFALSE;
    }
    evtSel = kPassINELgtZERO;

    if (!fEventCuts.PassedCut(AliEventCuts::kCorrelations)) {
      AliDebug(3, "Event with PileUp");
      return kFALSE;
    }
    evtSel = kCorrelations;

    if (!fEventCuts.PassedCut(AliEventCuts::kPileUp)) {
      AliDebug(3, "Event with PileUp");
      return kFALSE;
    }
    evtSel = kIsNotPileup;


    if (!fEventCuts.PassedCut(AliEventCuts::kVertex) || !fEventCuts.PassedCut(AliEventCuts::kVertexQuality)) {
      AliDebug(3, "Event doesn't pass has good vtx sel");
      return kFALSE;
    }
    evtSel = kHasRecVtx;

    if (!fEventCuts.PassedCut(AliEventCuts::kVertexPosition)) {
      AliDebugF(3, "Vertex with Z>%f cm", fMaxVtxZCut);
      return kFALSE;
    }
    evtSel = kHasGoodVtxZ;

    return kFALSE;
  } else { // If not Eventcut used

    if (fRejIncDAQ && fESD->IsIncompleteDAQ()) {
      AliDebug(3, "Event with incomplete DAQ");
      return kFALSE;
    }
    evtSel = kIsNotIncDAQ;

    UInt_t maskPhysSel =
      ((AliInputEventHandler *)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
    maskPhysSel &= fTriggerSel;
    if (!maskPhysSel) {
      AliDebugF(3, "Event doesn't pass physics evt. sel. for trigger %d", fTriggerSel);
      return kFALSE;
    }
    evtSel = kPassTrig;

    if (fChkIsEventINELgtZERO &&
        !(AliESDtrackCuts::GetReferenceMultiplicity(fESD, AliESDtrackCuts::kTracklets, 1.0) >= 1)) {
      AliDebug(3, "Event doesn't pass IsEventINELgtZERO selection");
      return kFALSE;
    }
    evtSel = kPassINELgtZERO;

    if (fDoSPDCvsTCut) {
      AliAnalysisUtils utils;
      if (utils.IsSPDClusterVsTrackletBG(fESD)) {
        AliDebug(3, "Event with incompatible SPD clusters and tracklet");
        return kFALSE;
      }
    }
    evtSel = kPassSPDclsVsTCut;

    if (!(fPlpType & BIT(kNoPileup))) {
      UInt_t lFlag = IsPileup();
      if (lFlag & BIT(kPileupSPD)) {
        AliDebug(3, "Pileup event from IsPileupFromSPD");
        return kFALSE;
      }
      evtSel = kIsNotPileupSPD;
      if (lFlag & BIT(kPileupInMultBins)) {
        AliDebug(3, "Pileup event from IsPileupSPDinMultBins");
        return kFALSE;
      }
      evtSel = kIsNotPileupSPDinMultBins;
      if (lFlag & BIT(kPileupMV)) {
        AliDebug(3, "Pileup event from IsPileupMV");
        return kFALSE;
      }
      evtSel = kIsNotPileupMV;
    }
    evtSel = kIsNotPileup;

    if (fUseSelectVertex2015pp && !SelectVertex2015pp()) {
      AliDebug(3, "Event doesn't pass vtx 2015 pp selection sel");
      return kFALSE;
    }
    evtSel = kHasRecVtx;

    if (!IsGoodVtxZ()) {
      AliDebugF(3, "Vertex with Z>%f cm", fMaxVtxZCut);
      return kFALSE;
    }
    evtSel = kHasGoodVtxZ;

    return kTRUE;
  }
  // Must ever be here
  AliDebug(3, "Event is rejected by unknown reason");
  return kFALSE;
}

//
//
//________________________________________________________________________
void AliAnalysisTaskSEITSsaSpectra::SetupStandardEventCutsForRun1()
{
  fChkIsSDDIn = kTRUE;
  fMultMethod = 0;
  fExtEventCuts = kFALSE;
  fRejIncDAQ = kFALSE;
  fDoSPDCvsTCut = kFALSE;
  fPlpType = BIT(kNoPileup);
  fTriggerSel = AliVEvent::kMB;
  fMaxVtxZCut = 10.;
  fReqBothVtx = kFALSE;
  fChkVtxSPDRes = kTRUE;
  fChkVtxZSep = kFALSE;
}

//
//
//________________________________________________________________________
void AliAnalysisTaskSEITSsaSpectra::SetupEventCutsForRun1pPb()
{
  fChkIsSDDIn = kTRUE;
  fMultMethod = 2;
  fExtEventCuts = kFALSE;
  fRejIncDAQ = kFALSE;
  fDoSPDCvsTCut = kFALSE;
  fPlpType = BIT(kNoPileup);
  fTriggerSel = AliVEvent::kINT7;
  fMaxVtxZCut = 10.;
  fReqBothVtx = kFALSE;
  fChkVtxSPDRes = kTRUE;
  fChkVtxZSep = kFALSE;
}

//
//
//________________________________________________________________________
void AliAnalysisTaskSEITSsaSpectra::SetupStandardEventCutsForRun2()
{
  fChkIsSDDIn = kTRUE;
  fMultMethod = 0;
  fExtEventCuts = kFALSE;
  fRejIncDAQ = kTRUE;
  fDoSPDCvsTCut = kTRUE;
  fPlpType = kPileupSPD;
  fMinPlpContribSPD = 3;
  fMinPlpZdistSPD = .8;
  fTriggerSel = AliVEvent::kINT7;
  fMaxVtxZCut = 10.;
  fReqBothVtx = kFALSE;
  fChkVtxSPDRes = kTRUE;
  fChkVtxZSep = kTRUE;
}

//
//
//________________________________________________________________________
bool AliAnalysisTaskSEITSsaSpectra::IsMultSelected()
{
  if (fMultMethod > 5) {
    AliWarning(". Skipping multiplicity selection");
    fMultMethod = 0;
  }
  if (!fMultMethod)
    return kTRUE;              // skip multiplicity check
  else if (fMultMethod == 1) { // New multiplicity/centrality class framework
    AliMultSelection *fMultSel = (AliMultSelection *)fESD->FindListObject("MultSelection");
    if (!fMultSel) {
      // If you get this warning (and lPercentiles 300) please check that the AliMultSelectionTask actually ran (before
      // your task)
      AliWarning("AliMultSelection object not found!");
      return kFALSE;
    } else {
      // Event selection is embedded in the Multiplicity estimator so that the Multiplicity percentiles are well defined
      // and refer to the same sample
      fEvtMult = fMultSel->GetMultiplicityPercentile(fMultEstimator.Data(), fMultEvSel);
    }
  } else if (fMultMethod == 2) { // OLD multiplicity/centrality class framework
    AliCentrality *centrality = fESD->GetCentrality();
    fEvtMult = centrality->GetCentralityPercentile(fMultEstimator.Data());
  } else if (fMultMethod == 3) { // selection on the event multiplicity based on global tracks
    // tracks+tracklets
    fEvtMult = (float)AliESDtrackCuts::GetReferenceMultiplicity(fESD, AliESDtrackCuts::kTrackletsITSTPC, 0.8);
  } else if (fMultMethod == 4) {
    // tracklets
    fEvtMult = (float)AliESDtrackCuts::GetReferenceMultiplicity(fESD, AliESDtrackCuts::kTracklets, 0.8);
  } else if (fMultMethod == 5) {
    // clusters in SPD1
    const AliMultiplicity *mult = fESD->GetMultiplicity();
    float nClu1 = (float)mult->GetNumberOfITSClusters(1);
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

  if (fPlpType & BIT(kNoPileup))
    return 1u; // skip pileup check;
  else if (fPlpType & BIT(kPileupSPD)) {
    if (fESD->IsPileupFromSPD(fMinPlpContribSPD, fMinPlpZdistSPD))
      lReturn |= BIT(kPileupSPD);
  } else if (fPlpType & BIT(kPileupInMultBins)) {
    if (fESD->IsPileupFromSPDInMultBins())
      lReturn |= BIT(kPileupInMultBins);
    ;
  } else if (fPlpType & BIT(kPileupMV)) {
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
bool AliAnalysisTaskSEITSsaSpectra::IsGoodVtxZ()
{
  const AliESDVertex *vtx = fESD->GetPrimaryVertex();
  fHistVtxZ->Fill(fEvtMult, vtx->GetZ());
  if (TMath::Abs(vtx->GetZ()) > fMaxVtxZCut)
    return kFALSE;

  return kTRUE;
}

//
//
//________________________________________________________________________
bool AliAnalysisTaskSEITSsaSpectra::IsGoodSPDvtxRes(const AliESDVertex *spdVtx)
{
  if (!spdVtx)
    return kFALSE;
  if (spdVtx->IsFromVertexerZ() && !(spdVtx->GetDispersion() < 0.04 && spdVtx->GetZRes() < 0.25))
    return kFALSE;
  return kTRUE;
}

//
//
//________________________________________________________________________
bool AliAnalysisTaskSEITSsaSpectra::SelectVertex2015pp()
{
  const AliESDVertex *trkVertex = fESD->GetPrimaryVertexTracks();
  const AliESDVertex *spdVertex = fESD->GetPrimaryVertexSPD();
  bool hasTrk = trkVertex->GetStatus();
  bool hasSPD = spdVertex->GetStatus();

  // Note that AliVertex::GetStatus checks that N_contributors is > 0
  // reject events if both are explicitly requested and none is available
  if (fReqBothVtx && !(hasSPD && hasTrk))
    return kFALSE;

  // reject events if none between the SPD or track verteces are available
  // if no trk vertex, try to fall back to SPD vertex;
  if (!hasTrk) {
    if (!hasSPD)
      return kFALSE;
    // on demand check the spd vertex resolution and reject if not satisfied
    if (fChkVtxSPDRes && !IsGoodSPDvtxRes(spdVertex))
      return kFALSE;
  } else {
    if (hasSPD) {
      // if enabled check the spd vertex resolution and reject if not satisfied
      // if enabled, check the proximity between the spd vertex and trak vertex, and reject if not satisfied
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
bool AliAnalysisTaskSEITSsaSpectra::DCAcut(double impactXY, double impactZ, double pt) const
{
  return (DCAcutXY(impactXY, pt) && DCAcutZ(impactZ, pt));
}

//
//
//________________________________________________________________________
bool AliAnalysisTaskSEITSsaSpectra::DCAcutXY(double impactXY, double pt) const
{
  // cut on transverse impact parameter updated on 20-5-2010
  // from the study of L. Milano, F. Prino on the ITS standalone tracks
  // using the common binning of the TPC tracks
  double xyMax = fDCAxyCutFunc->Eval(pt); // in micron
  AliDebugF(3, "Max value for the DCAxy Cut is:%f Measured value for DCAxy is:%f cut to %.0f sigmas\n", xyMax,
            TMath::Abs(impactXY) * 10000, fDCAxyCutFunc->GetParameter(3));
  if ((TMath::Abs(impactXY) * 10000) > xyMax)
    return kFALSE;

  return kTRUE;
}

//
//
//________________________________________________________________________
bool AliAnalysisTaskSEITSsaSpectra::DCAcutZ(double impactZ, double pt) const
{
  // cut on transverse impact parameter updated on 20-5-2010
  // from the study of L. Milano, F. Prino on the ITS standalone tracks
  // using the common binning of the TPC tracks

  double zMax = fDCAzCutFunc->Eval(pt); // in micron
  AliDebugF(3, "Max value for the DCAz Cut is:%f Measured value for DCAz is:%f cut to %.0f sigmas\n", zMax,
            TMath::Abs(impactZ) * 10000, fDCAzCutFunc->GetParameter(3));
  if ((TMath::Abs(impactZ) * 10000) > zMax)
    return kFALSE;

  return kTRUE;
}

//
//
//________________________________________________________________________
double AliAnalysisTaskSEITSsaSpectra::CookdEdx(double *s) const
{
  // truncated mean for the dEdx
  int nc = 0;
  double dedx[4] = { 0., 0., 0., 0. };
  for (int il = 0; il < 4; il++) { // count good (>0) dE/dx values
    if (s[il] > fMindEdx) {
      dedx[nc] = s[il];
      nc++;
    }
  }
  if (nc < fMinNdEdxSamples)
    return -1.;

  double tmp;
  int swap; // sort in ascending order
  do {
    swap = 0;
    for (int i = 0; i < nc - 1; i++) {
      if (dedx[i] <= dedx[i + 1])
        continue;
      tmp = dedx[i];
      dedx[i] = dedx[i + 1];
      dedx[i + 1] = tmp;
      swap++;
    }
  } while (swap);

  double sumamp = 0, sumweight = 0;
  double weight[4] = { 1., 1., 0., 0. };
  if (nc == 3)
    weight[1] = 0.5;
  else if (nc < 3)
    weight[1] = 0.;
  for (int i = 0; i < nc; i++) {
    sumamp += dedx[i] * weight[i];
    sumweight += weight[i];
  }
  return sumamp / sumweight;
}

//
//
//________________________________________________________________________
double AliAnalysisTaskSEITSsaSpectra::Eta2y(double pt, double m, double eta) const
{
  // convert eta to y
  double mt = TMath::Sqrt(m * m + pt * pt);
  return TMath::ASinH(pt / mt * TMath::SinH(eta));
}

//
//
//________________________________________________________________________
void AliAnalysisTaskSEITSsaSpectra::AnalyseMCParticles(AliMCEvent *lMCevent, EEvtCut_Type lastEvtCutPassed,
                                                       bool lHasGoodVtxGen)
{
  // first loop on stack, and get primary pi,k,p spectra with truth MC values
  int nTrackMC = (lMCevent) ? lMCevent->GetNumberOfTracks() : 0;
  for (int i_mcTrk = 0; i_mcTrk < nTrackMC; ++i_mcTrk) {
    AliMCParticle *mcTrk = ((AliMCParticle *)lMCevent->GetTrack(i_mcTrk));
    int pdg = mcTrk->PdgCode();
    int i_chg = pdg >= 0 ? 0 : 1; // only works for charged pi,K,p (0 Pos, 1 Neg)
    int iPart = -1;
    if (TMath::Abs(pdg) == 211)
      iPart = 0; // select Pi+/Pi- only
    else if (TMath::Abs(pdg) == 321)
      iPart = 1; // select K+/K- only
    else if (TMath::Abs(pdg) == 2212)
      iPart = 2; // select p+/p- only
    else
      continue;

    double mcPt = mcTrk->Pt();
    double mcP  = mcTrk->P();
    bool lIsPhysPrimary = lMCevent->IsPhysicalPrimary(i_mcTrk);
    for (int istep = (int)kIsReadable; istep <= (int)lastEvtCutPassed; ++istep)//before the cut on pt and rap is applied
      if(lIsPhysPrimary) fHistMCGenCharged->Fill(fEvtMult, TMath::Abs(mcP), istep);

    if (mcPt > 1.0)
      continue; // pt cut
    double mcEta = mcTrk->Eta();
    double mcRap = mcTrk->Y() + fCMSRapFct;
    if (fFillNtuple) {
      // filling MC ntuple
      float xntMC[8];
      int indexMC = 0;
      xntMC[indexMC++] = (float)fEvtMult;
      xntMC[indexMC++] = (float)mcPt;
      xntMC[indexMC++] = (float)pdg;
      xntMC[indexMC++] = (float)i_chg;
      xntMC[indexMC++] = (float)mcEta;
      xntMC[indexMC++] = (float)mcRap;
      xntMC[indexMC++] = (float)lIsPhysPrimary;
      //xntMC[indexMC++] = (float)lastEvtCutPassed;
      xntMC[indexMC++] = (float)fESD->GetRunNumber();

      fNtupleMC->Fill(xntMC);
    }

    if (!IsRapIn(mcRap))
      continue; // rapidity cut
    if (!lIsPhysPrimary)
      continue; // primary cut

    int index = iPart * kNchg + i_chg;
    for (int istep = (int)kIsReadable; istep <= (int)lastEvtCutPassed; ++istep) {
      fHistMCPart[index]->Fill(fEvtMult, TMath::Abs(mcPt), istep);
      if (lHasGoodVtxGen)
        fHistMCPartGoodGenVtxZ[index]->Fill(fEvtMult, TMath::Abs(mcPt), istep);
    }
  } // end MC stack loop
}

//
//
//________________________________________________________________________
void AliAnalysisTaskSEITSsaSpectra::PostAllData()
{
  PostData(1, fOutput);
  PostData(2, fListCuts);
  if (fFillNtuple)
    PostData(3, fListTree);

  return;
}

//
//
//________________________________________________________________________
double AliAnalysisTaskSEITSsaSpectra::BetheITSsaHybrid(double p, double mass) const
{
  //
  // returns AliExternalTrackParam::BetheBloch normalized to
  // fgMIP at the minimum. The PHOBOS parameterization is used for beta*gamma>0.76.
  // For beta*gamma<0.76 a polinomial function is used

  Double_t fBBsaHybrid[10];
  Double_t fBBsaElectron[6];

  if(!fIsMC) {//DATA
    if(fIsNominalBfield){//nominal magnetic field (0.5T)
      fBBsaHybrid[0]=1.43505E7;  //PHOBOS+Polinomial parameterization
      fBBsaHybrid[1]=49.3402;
      fBBsaHybrid[2]=1.77741E-7;
      fBBsaHybrid[3]=1.77741E-7;
      fBBsaHybrid[4]=1.01311E-7;
      fBBsaHybrid[5] = -2.; //exponent of beta in the PHOBOS
      fBBsaHybrid[6]=77.2777;
      fBBsaHybrid[7]=33.4099;
      fBBsaHybrid[8]=46.0089;
      fBBsaHybrid[9]=-2.26583;
      fBBsaElectron[0]=4.05799E6;  //electrons in the ITS
      fBBsaElectron[1]=38.5713;
      fBBsaElectron[2]=1.46462E-7;
      fBBsaElectron[3]=1.46462E-7;
      fBBsaElectron[4]=4.40284E-7;
      fBBsaElectron[5]=-2.;
    }
    else{//DATA: lowB field (0.2T)
      //PIONS
      if(mass>0.13 && mass<0.14){
        fBBsaHybrid[0]=2.0014e7;//E0  //PHOBOS+Polinomial parameterization
        fBBsaHybrid[1]=23.446696;//b
        fBBsaHybrid[2]=2.4831e-7;//a
        fBBsaHybrid[3]=2.4831e-7;//c
        fBBsaHybrid[4]=1.4825e-7;//d
        fBBsaHybrid[5]=-2.16;
        fBBsaHybrid[6]=4.913234;//p0
        fBBsaHybrid[7]=242.095214;//p1
        fBBsaHybrid[8]=-227.007773;//p2
        fBBsaHybrid[9]=131.060250;//p3
      }
      //KAONS
      else if(mass>0.4 && mass<0.5){
        fBBsaHybrid[0]=1.1603e7;//E0  //PHOBOS+Polinomial parameterization
        fBBsaHybrid[1]=26.908492;//b
        fBBsaHybrid[2]=2.1257e-4;//a
        fBBsaHybrid[3]=-1.6542e-4;//c
        fBBsaHybrid[4]=8.3836e-8;//d
        fBBsaHybrid[5]=-1.98;
        fBBsaHybrid[6]=-85.739134;//p0
        fBBsaHybrid[7]=379.078774;//p1
        fBBsaHybrid[8]=-241.702618;//p2
        fBBsaHybrid[9]=91.297155;//p3
      }
      else{//PROTONS and DEUTERONS
        fBBsaHybrid[0]=1.1644e7;//E0  //PHOBOS+Polinomial parameterization
        fBBsaHybrid[1]=42.148199;//b
        fBBsaHybrid[2]=1.6044e-4;//a
        fBBsaHybrid[3]=-9.9356e-5;//c
        fBBsaHybrid[4]=8.3214e-8;//d
        fBBsaHybrid[5]=-1.82;
        fBBsaHybrid[6]=-37.982015;//p0
        fBBsaHybrid[7]=245.348026;//p1
        fBBsaHybrid[8]=-116.456305;//p2
        fBBsaHybrid[9]=50.786350;//p3
      }

      fBBsaElectron[0]=66.524096;//E0 //electrons in the ITS
      fBBsaElectron[1]=106.793254;//b
      fBBsaElectron[2]=1.0939e-1;//a
      fBBsaElectron[3]=-3.8033e-3;//c
      fBBsaElectron[4]=-1.6387e-3;//d
      fBBsaElectron[5]=-1.60;//exp
    }

  } else {//MC

    if(fIsNominalBfield){//nominal magnetic field
      fBBsaHybrid[0]=1.05381E7; //PHOBOS+Polinomial parameterization
      fBBsaHybrid[1]=89.3933;
      fBBsaHybrid[2]=2.4831E-7;
      fBBsaHybrid[3]=2.4831E-7;
      fBBsaHybrid[4]=7.80591E-8;
      fBBsaHybrid[5]=-2.;
      fBBsaHybrid[6]=62.9214;
      fBBsaHybrid[7]=32.347;
      fBBsaHybrid[8]=58.7661;
      fBBsaHybrid[9]=-3.39869;
      fBBsaElectron[0]=2.26807E6; //electrons in the ITS
      fBBsaElectron[1]=99.985;
      fBBsaElectron[2]=0.000714841;
      fBBsaElectron[3]=0.000259585;
      fBBsaElectron[4]=1.39412E-7;
      fBBsaElectron[5]=-2.;
    }
    else{//MC low B field

      //PIONS
      if(mass>0.13 && mass<0.14){

        fBBsaHybrid[0]=1.9993e7;//E0  //PHOBOS+Polinomial parameterization
        fBBsaHybrid[1]=23.503097;//b
        fBBsaHybrid[2]=2.4831e-7;//a
        fBBsaHybrid[3]=2.4831e-7;//c
        fBBsaHybrid[4]=1.4809e-7;//d
        fBBsaHybrid[5]=-2.16;
        fBBsaHybrid[6]=4.930021;//p0
        fBBsaHybrid[7]=241.625126;//p1
        fBBsaHybrid[8]=-225.955937;//p2
        fBBsaHybrid[9]=130.683475;//p3
      }
      //KAONS
      else if(mass>0.4 && mass<0.5){

        fBBsaHybrid[0]=1.1916e7;//E0  //PHOBOS+Polinomial parameterization
        fBBsaHybrid[1]=18.519683;//b
        fBBsaHybrid[2]=2.5525e-4;//a
        fBBsaHybrid[3]=-2.1776e-4;//c
        fBBsaHybrid[4]=8.4678e-8;//d
        fBBsaHybrid[5]=-2.12;
        fBBsaHybrid[6]=-87.655919;//p0
        fBBsaHybrid[7]=383.014648;//p1
        fBBsaHybrid[8]=-244.355671;//p2
        fBBsaHybrid[9]=92.010403;//p3
      }
      else{//PROTONS and DEUTERONS
        fBBsaHybrid[0]=1.1729e7;//E0  //PHOBOS+Polinomial parameterization
        fBBsaHybrid[1]=42.145678;//b
        fBBsaHybrid[2]=1.5571e-4;//a
        fBBsaHybrid[3]=-9.4656e-5;//c
        fBBsaHybrid[4]=8.6941e-8;//d
        fBBsaHybrid[5]=-1.82;
        fBBsaHybrid[6]=-35.460903;//p0
        fBBsaHybrid[7]=239.560138;//p1
        fBBsaHybrid[8]=-112.309480;//p2
        fBBsaHybrid[9]=49.977767;//p3
      }

      fBBsaElectron[0]=67.425232;//E0 //electrons in the ITS
      fBBsaElectron[1]=106.744327;//b
      fBBsaElectron[2]=1.0971e-1;//a
      fBBsaElectron[3]=-4.1163e-3;//c
      fBBsaElectron[4]=-1.9255e-3;//d
      fBBsaElectron[5]=-1.60;//exp of beta
    }

}

  Double_t bg=p/mass;
  Double_t beta = bg/TMath::Sqrt(1.+ bg*bg);
  Double_t gamma=bg/beta;
  Double_t bb=1.;

  Double_t betagcut = 0.76;

  if(!fIsNominalBfield){
    if(mass>0.13 && mass<0.14) betagcut = 1.50;
    else if(mass>0.4 && mass<0.5) betagcut = 1.04;
    else betagcut = 1.0;
  }

  Double_t par[10];
  //parameters for pi, K, p
  for(Int_t ip=0; ip<10;ip++) par[ip]=fBBsaHybrid[ip];
    //if it is an electron the PHOBOS part of the parameterization is tuned for e
    //in the range used for identification beta*gamma is >0.76 for electrons
    //To be used only between 100 and 160 MeV/c
    if(mass>0.0005 && mass<0.00052) for(Int_t ip=0; ip<6; ip++) par[ip]=fBBsaElectron[ip];

      if(gamma>=0. && beta>0. && bg>0.1){
        if(bg>betagcut){//PHOBOS
          Double_t eff=1.0;
          if(bg<par[2])
           eff=(bg-par[3])*(bg-par[3])+par[4];
          else
           eff=(par[2]-par[3])*(par[2]-par[3])+par[4];
          bb=(par[1]+2.0*TMath::Log(gamma)-beta*beta)*(par[0]/(TMath::Power(beta, -par[5])))*eff;
        }else{//Polinomial
          bb=par[6] + par[7]/bg + par[8]/(bg*bg) + par[9]/(bg*bg*bg);
        }
      }
  return bb;
}

//
//
//________________________________________________________________________
int AliAnalysisTaskSEITSsaSpectra::GetTrackPid(AliESDtrack *track, double *logdiff) const
{
  AliPID::EParticleType iType[4] = { AliPID::kPion, AliPID::kKaon, AliPID::kProton, AliPID::kDeuteron };

  int pid = -1;

  double dEdxLay[4];
  track->GetITSdEdxSamples(dEdxLay);
  double dedx = track->GetITSsignal();
  float p = track->GetP();
  if (fIsMC && fSmearMC) {
    dedx = fRandGener->Gaus(dedx, fSmeardEdx * dedx);
    p = fRandGener->Gaus(p, fSmearP * p);
  }

  double bbtheo[4];
  for (int i = 0; i < 4; i++) {
    float mass = AliPID::ParticleMass(iType[i]);
    //bbtheo[i] = fITSPIDResponse->BetheITSsaHybrid(p, mass);
    bbtheo[i] = BetheITSsaHybrid(p, mass);
    logdiff[i] = TMath::Log(dedx) - TMath::Log(bbtheo[i]);
  }

  UInt_t clumap = track->GetITSClusterMap();
  int nPtsForPid = 0;
  for (int j = 2; j < 6; j++)
    if (TESTBIT(clumap, j))
      nPtsForPid++;

  float resodedx = fITSPIDResponse->GetResolution(1, nPtsForPid, kTRUE);
  switch (fPidMethod) {
    case kNSigCut: {
      double nSigmaMin = 99999;
      for (int i_spc = 0; i_spc < kNspc; i_spc++) {
        double bb = bbtheo[i_spc];
        double nSigma = TMath::Abs((dedx - bb) / (resodedx * bb));
        if (nSigma < nSigmaMin) {
          nSigmaMin = nSigma;
          pid = iType[i_spc];
        }
      }
      pid = (nSigmaMin < fMinNSigma) ? pid : -1;
      break;
    }
    case kMeanCut: {
      for (int i_spc = 0; i_spc < kNspc; i_spc++) {
        if (dedx < bbtheo[0]) {
          double nsigma = TMath::Abs((dedx - bbtheo[0]) / (resodedx * bbtheo[0]));
          pid = (nsigma < 2.) ? AliPID::kPion : -1;
          break;
        }
        double bb = bbtheo[i_spc];
        int edge = (dedx < bb) ? i_spc - 1 : i_spc + 1;
        double bbdistance = TMath::Abs((bbtheo[i_spc] - bbtheo[edge]) / 2);

        if (bbdistance != 0.f) {
          double nsigma = TMath::Abs((dedx - bb) / bbdistance);
          if (nsigma < 1.)
            pid = iType[i_spc];
        } else {
          AliInfo("Zero distance");
          pid = AliPID::kPion;
        }
      }
      break;
    }
    case kLanGaus: {
      double prior[AliPID::kSPECIES];
      double probITS[AliPID::kSPECIES];
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

  // Sigma Separation
  for (int i_spc = 0; i_spc < kNspc; ++i_spc) {
    int i_chg = (track->GetSign() > 0) ? 0 : 1;
    int index = i_spc * kNchg + i_chg;
    double bb = bbtheo[i_spc];
    fHistNSigmaSep[index]->Fill(track->Pt(), ((dedx - bb) / (resodedx * bb)));
  }

  // fill charge distribution histo to check the calibration
  for (int j = 0; j < 4; j++) {
    if (dEdxLay[j] < fMindEdx)
      continue;
    fHistCharge[j]->Fill(dEdxLay[j]);
  }

  return (pid == -1) ? 0 : pid;
}

//
//
//________________________________________________________________________
int AliAnalysisTaskSEITSsaSpectra::GetMostProbable(const double *pDens, const double *prior) const
{
  // get the most probable particle id hypothesis
  // assuming the a priori probabilities "prior"
  double max = 0.;
  int id = AliPID::kPion;
  for (int i = 0; i < AliPID::kSPECIES; ++i) {
    double prob = pDens[i] * prior[i];
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
void AliAnalysisTaskSEITSsaSpectra::GetPriors(const AliVTrack *track, double *prior) const
{
  const char *prtName[AliPID::kSPECIES] = { "El", "Mu", "Pi", "Ka", "Pr" };
  double pt = track->Pt();

  if (pt < 0.08)
    pt = 0.0801;
  else if (pt > 0.9999)
    pt = 0.9999;

  float usedCent = fEvtMult; //+ 5*(cent>0);
  if (usedCent < -.99)
    usedCent = -0.9;
  if (usedCent > 99.9)
    usedCent = 99.9;

  for (int i = 0; i < AliPID::kSPECIES; ++i)
    prior[i] = 1. / AliPID::kSPECIES;

  if (!fUseDefaultPriors) {
    if (!fListPriors) {
      AliInfo("Errors use priors set but no priors list found.");
      AliInfo("Default(equal) priors used.");
      return;
    } else {
      for (int i = 0; i < AliPID::kSPECIES; i++) {
        TH2F *hPriors = dynamic_cast<TH2F *>(fListPriors->FindObject(Form("hDefaultITSsa%sPriors", prtName[i])));
        if (!hPriors) {
          AliInfo("Errors use priors set but no histogram with priors found.");
          AliInfo("Default(equal) priors used.");
          return;
        }
        prior[i] = hPriors->Interpolate(usedCent, pt);
      }
    }
  } else {
    if (pt > .160)
      prior[0] = 0;
    prior[1] = 0.;
  }

  double sumP = 0;
  for (int i = 0; i < AliPID::kSPECIES; i++)
    sumP += prior[i];
  for (int i = 0; i < AliPID::kSPECIES; i++)
    prior[i] /= sumP;

  return;
}

//
//
//________________________________________________________________________
void AliAnalysisTaskSEITSsaSpectra::ComputeBayesProbabilities(double *probs, const double *pDens, const double *prior)
{
  // Calculate bayesian probabilities
  double sum = 0.;
  for (int i = 0; i < AliPID::kSPECIES; i++)
    sum += pDens[i] * prior[i];

  if (sum <= 0) {
    AliError("Invalid probability densities or priors");
    for (int i = 0; i < AliPID::kSPECIES; i++)
      probs[i] = -1;
    return;
  }

  for (int i = 0; i < AliPID::kSPECIES; i++)
    probs[i] = pDens[i] * prior[i] / sum;
}

//
//
//________________________________________________________________________
void AliAnalysisTaskSEITSsaSpectra::SetBins(const int nbins, double min, double max, double *bins)
{
  const double delta = (max - min) / nbins;
  for (int iB = 0; iB < nbins; ++iB) {
    bins[iB] = min + iB * delta;
  }
  bins[nbins] = max;
}

//
//
//________________________________________________________________________
void AliAnalysisTaskSEITSsaSpectra::SetCentBins(int nbins, double *bins) { fCentBins.Set(nbins + 1, bins); }

//
//
//________________________________________________________________________
void AliAnalysisTaskSEITSsaSpectra::SetDCABins(int nbins, double *bins) { fDCABins.Set(nbins + 1, bins); }

//
//
//________________________________________________________________________
void AliAnalysisTaskSEITSsaSpectra::SetPtBins(int nbins, double *bins) { fPtBins.Set(nbins + 1, bins); }
