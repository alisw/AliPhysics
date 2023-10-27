/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

//______________________________________________________________________________
// Analysis task for high pt particle correlations
// author: Jasper Parkila, D.J. Kim
// ALICE Group University of Jyvaskyla
// Finland
//////////////////////////////////////////////////////////////////////////////

#include "TRandom.h"
#include "TList.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAODVertex.h"
#include "TLorentzVector.h"
#include <TFile.h>
#include <TRandom3.h>
#include <TH3D.h>

#include "AliAnalysisUtils.h"
#include "AliAODMCParticle.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliESDtrackCuts.h"

#include "AliAODMCHeader.h"
#include "AliJetContainer.h"
#include "TMath.h"
#include "AliJHFTagging.h"
#include "AliExternalTrackParam.h"
#include "AliVertexerTracks.h"
#include "AliHFJetsTagging.h"
#include "AliHFJetsTaggingVertex.h"
#include "AliRDHFJetsCutsVertex.h"
#include "AliFJWrapper.h"
#include <vector>
#include <algorithm>

ClassImp(AliJHFTagging);

//______________________________________________________________________________
AliJHFTagging::AliJHFTagging() : AliAnalysisTaskEmcalJet(),
                                 fInputList(0),
                                 fInputListb(0x0),
                                 fInputListc(0x0),
                                 fInputListlf(0x0),
                                 fFillControlHists(false),
                                 fThresholdIP(2.5),
                                 fMaxDespersion(0.03),
                                 fThresholdLxy(7.),
                                 fEmbeddPerpendicular(false),
                                 fEvtNum(0),
                                 fDiamond(0x0),
                                 fVertexer(0x0),
                                 fRandom(new TRandom3(0)),
                                 fFastJetWrapper(0x0),
                                 fTrackGenerator(0x0),
                                 fMCArray(0x0),
                                 // Bjet cuts
                                 fTCMinTrackPt(0.5),
                                 fTCMinClusTPC(80),
                                 fTCMinHitsITS(2),
                                 fTCMaxChi2pNDF(5.),
                                 fTCMaxIPxy(1.),
                                 fTCMaxIPz(5.),
                                 fTCMaxDecayLength(5),
                                 fTCMaxDCATrackJet(0.07),
                                 fhistInclusiveJetCuts(0x0),
                                 fhistbJetCuts(0x0),
                                 fhistcJetCuts(0x0),
                                 fhistlfJetCuts(0x0),
                                 fh1dJetRecPtAcceptedunCorr(0x0),
                                 f2histRhoVsDeltaPt(0x0),
                                 f2histRhoVsDeltaPtFirst(0x0),
                                 f2histRhoVsDeltaPtSecond(0x0),
                                 f2histRhoVsDeltaPtThird(0x0),
                                 fh2DeltaPtEmbeddCorrelationPerpendicular(0x0),
                                 fh1dJetGenPt(0x0),
                                 fh1dJetGenPtUnidentified(0x0),
                                 fh1dJetGenPtudsg(0x0),
                                 fh1dJetGenPtc(0x0),
                                 fh1dJetGenPtb(0x0),
                                 fh1dJetRecPt(0x0),
                                 fh1dJetRecPtAccepted(0x0),
                                 fh1dJetRecEtaPhiAccepted(0x0),
                                 fh1dJetRecPtUnidentified(0x0),
                                 fh1dJetRecPtudsg(0x0),
                                 fh1dJetRecPtc(0x0),
                                 fh1dJetRecPtb(0x0),
                                 fh1dJetRecPtUnidentifiedAccepted(0x0),
                                 fh1dJetRecPtudsgAccepted(0x0),
                                 fh1dJetRecPtcAccepted(0x0),
                                 fh1dJetRecPtbAccepted(0x0),
                                 fh2dJetGenPtVsJetRecPt(0x0),
                                 fh2dJetGenPtVsJetRecPtb(0x0),
                                 fh2dJetGenPtVsJetRecPtc(0x0),
                                 fh2dJetGenPtVsJetRecPtudsg(0x0),
                                 fh2dJetSignedImpParXY(0x0),
                                 fh2dJetSignedImpParXYUnidentified(0x0),
                                 fh2dJetSignedImpParXYudsg(0x0),
                                 fh2dJetSignedImpParXYb(0x0),
                                 fh2dJetSignedImpParXYc(0x0),
                                 fh2dJetSignedImpParXYSignificance(0x0),
                                 fh2dJetSignedImpParXYSignificanceUnidentified(0x0),
                                 fh2dJetSignedImpParXYSignificanceudsg(0x0),
                                 fh2dJetSignedImpParXYSignificanceb(0x0),
                                 fh2dJetSignedImpParXYSignificancec(0x0),
                                 fh2dJetSignedImpParXYZ(0x0),
                                 fh2dJetSignedImpParXYZUnidentified(0x0),
                                 fh2dJetSignedImpParXYZudsg(0x0),
                                 fh2dJetSignedImpParXYZb(0x0),
                                 fh2dJetSignedImpParXYZc(0x0),
                                 fh2dJetSignedImpParXYZSignificance(0x0),
                                 fh2dJetSignedImpParXYZSignificanceUnidentified(0x0),
                                 fh2dJetSignedImpParXYZSignificanceudsg(0x0),
                                 fh2dJetSignedImpParXYZSignificanceb(0x0),
                                 fh2dJetSignedImpParXYZSignificancec(0x0),
                                 fh2dJetSignedImpParXYSignificanceFirst(0x0),
                                 fh2dJetSignedImpParXYSignificanceUnidentifiedFirst(0x0),
                                 fh2dJetSignedImpParXYSignificanceudsgFirst(0x0),
                                 fh2dJetSignedImpParXYSignificancebFirst(0x0),
                                 fh2dJetSignedImpParXYSignificancecFirst(0x0),
                                 // Second
                                 fh2dJetSignedImpParXYSignificanceSecond(0x0),
                                 fh2dJetSignedImpParXYSignificanceUnidentifiedSecond(0x0),
                                 fh2dJetSignedImpParXYSignificanceudsgSecond(0x0),
                                 fh2dJetSignedImpParXYSignificancebSecond(0x0),
                                 fh2dJetSignedImpParXYSignificancecSecond(0x0),
                                 // Third
                                 fh2dJetSignedImpParXYSignificanceThird(0x0),
                                 fh2dJetSignedImpParXYSignificanceUnidentifiedThird(0x0),
                                 fh2dJetSignedImpParXYSignificanceudsgThird(0x0),
                                 fh2dJetSignedImpParXYSignificancebThird(0x0),
                                 fh2dJetSignedImpParXYSignificancecThird(0x0),
                                 fJetContainerMC(0x0),
                                 fJetContainerData(0x0),
                                 fAODIn(0x0),
                                 fPrimaryVertex(0x0),
                                 // Secondary Vertex
                                 fDoSVAnalysis(false),
                                 fDoTrackCountingAnalysis(true),
                                 fUsePartonDef(true),
                                 fVtxTagger3Prong(0x0),
                                 fVtxTagger2Prong(0x0),
                                 fEsdTrackCuts(0x0),
                                 fHistSV2Prong(0x0),
                                 fHistSV2ProngUnidentified(0x0),
                                 fHistSV2Prongb(0x0),
                                 fHistSV2Prongc(0x0),
                                 fHistSV2Pronglf(0x0),
                                 fHistDispersion2Prong(0x0),
                                 fHistDispersion2ProngUnidentified(0x0),
                                 fHistDispersion2Prongb(0x0),
                                 fHistDispersion2Prongc(0x0),
                                 fHistDispersion2Pronglf(0x0),
                                 fHistSV3Prong(0x0),
                                 fHistSV3ProngUnidentified(0x0),
                                 fHistSV3Prongb(0x0),
                                 fHistSV3Prongc(0x0),
                                 fHistSV3Pronglf(0x0),
                                 fHistDispersion3Prong(0x0),
                                 fHistDispersion3ProngUnidentified(0x0),
                                 fHistDispersion3Prongb(0x0),
                                 fHistDispersion3Prongc(0x0),
                                 fHistDispersion3Pronglf(0x0)
{
  DefineOutput(1, TList::Class());
  SetMakeGeneralHistograms(true);
}

//______________________________________________________________________________
AliJHFTagging::AliJHFTagging(const char *name) : AliAnalysisTaskEmcalJet(name, true),
                                                 fInputList(0),
                                                 fInputListb(0x0),
                                                 fInputListc(0x0),
                                                 fInputListlf(0x0),
                                                 fFillControlHists(false),
                                                 fThresholdIP(2.5),
                                                 fMaxDespersion(0.03),
                                                 fThresholdLxy(7.),
                                                 fEmbeddPerpendicular(false),
                                                 fEvtNum(0),
                                                 fDiamond(0x0),
                                                 fVertexer(0x0),
                                                 fRandom(new TRandom3(0)),
                                                 fFastJetWrapper(0x0),
                                                 fTrackGenerator(0x0),
                                                 fMCArray(0x0),
                                                 // Bjet cuts
                                                 fTCMinTrackPt(0.5),
                                                 fTCMinClusTPC(80),
                                                 fTCMinHitsITS(2),
                                                 fTCMaxChi2pNDF(5.),
                                                 fTCMaxIPxy(1.),
                                                 fTCMaxIPz(5.),
                                                 fTCMaxDecayLength(5),
                                                 fTCMaxDCATrackJet(0.07),
                                                 fhistInclusiveJetCuts(0x0),
                                                 fhistbJetCuts(0x0),
                                                 fhistcJetCuts(0x0),
                                                 fhistlfJetCuts(0x0),
                                                 fh1dJetRecPtAcceptedunCorr(0x0),
                                                 f2histRhoVsDeltaPt(0x0),
                                                 f2histRhoVsDeltaPtFirst(0x0),
                                                 f2histRhoVsDeltaPtSecond(0x0),
                                                 f2histRhoVsDeltaPtThird(0x0),
                                                 fh2DeltaPtEmbeddCorrelationPerpendicular(0x0),
                                                 fh1dJetGenPt(0x0),
                                                 fh1dJetGenPtUnidentified(0x0),
                                                 fh1dJetGenPtudsg(0x0),
                                                 fh1dJetGenPtc(0x0),
                                                 fh1dJetGenPtb(0x0),
                                                 fh1dJetRecPt(0x0),
                                                 fh1dJetRecPtAccepted(0x0),
                                                 fh1dJetRecEtaPhiAccepted(0x0),
                                                 fh1dJetRecPtUnidentified(0x0),
                                                 fh1dJetRecPtudsg(0x0),
                                                 fh1dJetRecPtc(0x0),
                                                 fh1dJetRecPtb(0x0),
                                                 fh1dJetRecPtUnidentifiedAccepted(0x0),
                                                 fh1dJetRecPtudsgAccepted(0x0),
                                                 fh1dJetRecPtcAccepted(0x0),
                                                 fh1dJetRecPtbAccepted(0x0),
                                                 fh2dJetGenPtVsJetRecPt(0x0),
                                                 fh2dJetGenPtVsJetRecPtb(0x0),
                                                 fh2dJetGenPtVsJetRecPtc(0x0),
                                                 fh2dJetGenPtVsJetRecPtudsg(0x0),
                                                 fh2dJetSignedImpParXY(0x0),
                                                 fh2dJetSignedImpParXYUnidentified(0x0),
                                                 fh2dJetSignedImpParXYudsg(0x0),
                                                 fh2dJetSignedImpParXYb(0x0),
                                                 fh2dJetSignedImpParXYc(0x0),
                                                 fh2dJetSignedImpParXYSignificance(0x0),
                                                 fh2dJetSignedImpParXYSignificanceUnidentified(0x0),
                                                 fh2dJetSignedImpParXYSignificanceudsg(0x0),
                                                 fh2dJetSignedImpParXYSignificanceb(0x0),
                                                 fh2dJetSignedImpParXYSignificancec(0x0),
                                                 fh2dJetSignedImpParXYZ(0x0),
                                                 fh2dJetSignedImpParXYZUnidentified(0x0),
                                                 fh2dJetSignedImpParXYZudsg(0x0),
                                                 fh2dJetSignedImpParXYZb(0x0),
                                                 fh2dJetSignedImpParXYZc(0x0),
                                                 fh2dJetSignedImpParXYZSignificance(0x0),
                                                 fh2dJetSignedImpParXYZSignificanceUnidentified(0x0),
                                                 fh2dJetSignedImpParXYZSignificanceudsg(0x0),
                                                 fh2dJetSignedImpParXYZSignificanceb(0x0),
                                                 fh2dJetSignedImpParXYZSignificancec(0x0),
                                                 fh2dJetSignedImpParXYSignificanceFirst(0x0),
                                                 fh2dJetSignedImpParXYSignificanceUnidentifiedFirst(0x0),
                                                 fh2dJetSignedImpParXYSignificanceudsgFirst(0x0),
                                                 fh2dJetSignedImpParXYSignificancebFirst(0x0),
                                                 fh2dJetSignedImpParXYSignificancecFirst(0x0),
                                                 // Second
                                                 fh2dJetSignedImpParXYSignificanceSecond(0x0),
                                                 fh2dJetSignedImpParXYSignificanceUnidentifiedSecond(0x0),
                                                 fh2dJetSignedImpParXYSignificanceudsgSecond(0x0),
                                                 fh2dJetSignedImpParXYSignificancebSecond(0x0),
                                                 fh2dJetSignedImpParXYSignificancecSecond(0x0),
                                                 // Third
                                                 fh2dJetSignedImpParXYSignificanceThird(0x0),
                                                 fh2dJetSignedImpParXYSignificanceUnidentifiedThird(0x0),
                                                 fh2dJetSignedImpParXYSignificanceudsgThird(0x0),
                                                 fh2dJetSignedImpParXYSignificancebThird(0x0),
                                                 fh2dJetSignedImpParXYSignificancecThird(0x0),
                                                 fJetContainerMC(0x0),
                                                 fJetContainerData(0x0),
                                                 fAODIn(0x0),
                                                 fPrimaryVertex(0x0),
                                                 // Secondary Vertex
                                                 fDoSVAnalysis(false),
                                                 fDoTrackCountingAnalysis(true),
                                                 fUsePartonDef(true),
                                                 fVtxTagger3Prong(0x0),
                                                 fVtxTagger2Prong(0x0),
                                                 fEsdTrackCuts(0x0),
                                                 fHistSV2Prong(0x0),
                                                 fHistSV2ProngUnidentified(0x0),
                                                 fHistSV2Prongb(0x0),
                                                 fHistSV2Prongc(0x0),
                                                 fHistSV2Pronglf(0x0),
                                                 fHistDispersion2Prong(0x0),
                                                 fHistDispersion2ProngUnidentified(0x0),
                                                 fHistDispersion2Prongb(0x0),
                                                 fHistDispersion2Prongc(0x0),
                                                 fHistDispersion2Pronglf(0x0),
                                                 fHistSV3Prong(0x0),
                                                 fHistSV3ProngUnidentified(0x0),
                                                 fHistSV3Prongb(0x0),
                                                 fHistSV3Prongc(0x0),
                                                 fHistSV3Pronglf(0x0),
                                                 fHistDispersion3Prong(0x0),
                                                 fHistDispersion3ProngUnidentified(0x0),
                                                 fHistDispersion3Prongb(0x0),
                                                 fHistDispersion3Prongc(0x0),
                                                 fHistDispersion3Pronglf(0x0)
{
  DefineOutput(1, TList::Class());
  SetMakeGeneralHistograms(true);
}

//____________________________________________________________________________
AliJHFTagging::AliJHFTagging(const AliJHFTagging &ap) : AliAnalysisTaskEmcalJet(ap.GetName()),
                                                        fInputList(ap.fInputList),
                                                        fInputListb(ap.fInputListb),
                                                        fInputListc(ap.fInputListc),
                                                        fInputListlf(ap.fInputListlf),
                                                        fHFJetUtils(ap.fHFJetUtils.get()),
                                                        fFillControlHists(ap.fFillControlHists),
                                                        fThresholdIP(ap.fThresholdIP),
                                                        fMaxDespersion(ap.fMaxDespersion),
                                                        fThresholdLxy(ap.fThresholdLxy),
                                                        fEmbeddPerpendicular(ap.fEmbeddPerpendicular),
                                                        fDiamond(ap.fDiamond),
                                                        fVertexer(ap.fVertexer),
                                                        fRandom(ap.fRandom),
                                                        fFastJetWrapper(ap.fFastJetWrapper),
                                                        fTrackGenerator(ap.fTrackGenerator),
                                                        fMCArray(ap.fMCArray),
                                                        // Bjet cuts
                                                        fTCMinTrackPt(ap.fTCMinTrackPt),
                                                        fTCMinClusTPC(ap.fTCMinClusTPC),
                                                        fTCMinHitsITS(ap.fTCMinHitsITS),
                                                        fTCMaxChi2pNDF(ap.fTCMaxChi2pNDF),
                                                        fTCMaxIPxy(ap.fTCMaxIPxy),
                                                        fTCMaxIPz(ap.fTCMaxIPz),
                                                        fTCMaxDecayLength(ap.fTCMaxDecayLength),
                                                        fTCMaxDCATrackJet(ap.fTCMaxDCATrackJet),
                                                        fJetContainerMC(ap.fJetContainerMC),
                                                        fJetContainerData(ap.fJetContainerData),
                                                        fAODIn(ap.fAODIn),
                                                        fPrimaryVertex(ap.fPrimaryVertex),
                                                        // Secondary Vertex
                                                        fDoSVAnalysis(ap.fDoSVAnalysis),
                                                        fDoTrackCountingAnalysis(ap.fDoTrackCountingAnalysis),
                                                        fUsePartonDef(ap.fUsePartonDef),
                                                        fVtxTagger3Prong(ap.fVtxTagger3Prong),
                                                        fVtxTagger2Prong(ap.fVtxTagger2Prong),
                                                        fjetCuts3Prong(ap.fjetCuts3Prong.get()),
                                                        fjetCuts2Prong(ap.fjetCuts2Prong.get()),
                                                        fEsdTrackCuts(ap.fEsdTrackCuts)
{
  DefineOutput(1, TList::Class());
  SetMakeGeneralHistograms(true);
  AliInfo("----DEBUG AliJHFTagging COPY ----");
}

//_____________________________________________________________________________
AliJHFTagging &AliJHFTagging::operator=(const AliJHFTagging &ap)
{
  // assignment operator
  AliInfo("----DEBUG AliJHFTagging operator= ----");
  this->~AliJHFTagging();
  new (this) AliJHFTagging(ap);
  return *this;
}

//______________________________________________________________________________
void AliJHFTagging::Init()
{
  AliInfo("Doing initialization");
}

//______________________________________________________________________________
void AliJHFTagging::Terminate(Option_t *)
{
  AliInfo("AliJHFTagging Analysis DONE !!");
}

//______________________________________________________________________________
AliJHFTagging::~AliJHFTagging()
{
  delete fInputList;
  delete fInputListb;
  delete fInputListc;
  delete fInputListlf;

  delete fOutput;
  delete fJetContainerMC;
  delete fJetContainerData;
  delete fAODIn;
  delete fPrimaryVertex;
  delete fVtxTagger3Prong;
  delete fVtxTagger2Prong;
  delete fEsdTrackCuts;
  delete fRandom;
  delete fFastJetWrapper;
  delete fTrackGenerator;
  delete fVertexer;
  delete fDiamond;
}

// ########################################################################################
void AliJHFTagging::UserCreateOutputObjects()
{

  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  fInputList = new TClonesArray("AliEmcalJet", 2500);
  fInputList->SetOwner(true);

  if (fIsPythia)
  {
    fInputListb = new TClonesArray("AliEmcalJet", 2500);
    fInputListb->SetOwner(true);

    fInputListc = new TClonesArray("AliEmcalJet", 2500);
    fInputListc->SetOwner(true);

    fInputListlf = new TClonesArray("AliEmcalJet", 2500);
    fInputListlf->SetOwner(true);
  }

  if (fDoSVAnalysis)
  {
    fEsdTrackCuts = new AliESDtrackCuts("AliESDtrackCuts", "default");

    fEsdTrackCuts->SetRequireSigmaToVertex(false);
    fEsdTrackCuts->SetMinNClustersTPC(70);
    fEsdTrackCuts->SetMaxChi2PerClusterTPC(4);
    fEsdTrackCuts->SetRequireTPCRefit(true);
    fEsdTrackCuts->SetRequireITSRefit(true);
    fEsdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
    fEsdTrackCuts->SetMinDCAToVertexXY(fDoSVAnalysis ? 0. : 0.008);
    fEsdTrackCuts->SetEtaRange(-0.9, 0.9);
    fEsdTrackCuts->SetPtRange(1., 1.e10);

    fEsdTrackCuts->SetEtaRange(-0.8, 0.8);
    fEsdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);

    fjetCuts3Prong = std::make_unique<AliRDHFJetsCutsVertex>("jetCuts3Prong");
    fjetCuts3Prong->AddTrackCuts(fEsdTrackCuts);
    fjetCuts3Prong->SetNprongs(3);
    fjetCuts3Prong->SetMinPtHardestTrack(1.);
    fjetCuts3Prong->SetIsElec(false);       // true to select e in jet vertex
    fjetCuts3Prong->SetSecVtxWithKF(false); // default with StrLinMinDist

    fjetCuts2Prong = std::make_unique<AliRDHFJetsCutsVertex>("jetCuts2Prong");
    fjetCuts2Prong->AddTrackCuts(fEsdTrackCuts);
    fjetCuts2Prong->SetNprongs(2);
    fjetCuts2Prong->SetMinPtHardestTrack(1.);

    fVtxTagger3Prong = new AliHFJetsTaggingVertex();
    fVtxTagger3Prong->SetCuts(fjetCuts3Prong.get());

    fVtxTagger2Prong = new AliHFJetsTaggingVertex();
    fVtxTagger2Prong->SetCuts(fjetCuts2Prong.get());
  }

  if (fEmbeddPerpendicular)
  {
    fFastJetWrapper = new AliFJWrapper("FJWrapper", "FJWrapper");
    fFastJetWrapper->SetAreaType(fastjet::active_area);
    fFastJetWrapper->SetGhostArea(0.005);
    fFastJetWrapper->SetR(0.4);
    fFastJetWrapper->SetAlgorithm(fastjet::antikt_algorithm);
    fFastJetWrapper->SetRecombScheme(fastjet::pt_scheme);
    fTrackGenerator = new TRandom(0);
  }

  if (fIsPythia)
    fHFJetUtils = std::make_unique<AliHFJetsTagging>("fHFJetUtils");

  if (fFillControlHists)
  {
    MakeControlHistograms();
  }
}

// ################################################## Init histograms
void AliJHFTagging::MakeControlHistograms()
{

  TString BjetCuts[15] = {
      "all" /*0*/,
      "FilterBit 4" /*1*/,
      "p_{T} cut" /*2*/,
      "#eta cut" /*3*/,
      "TPC refit" /*4*/,
      "ITS refit" /*5*/,
      "SPD hits" /*6*/,
      "TPC clusters" /*7*/,
      "ITS clusters" /*8*/,
      "Chi2/NDF" /*9*/,
      "IP_{xy}" /*10*/,
      "IP_{z}" /*11*/,
      "DecayLength" /*12*/,
      "DCA Track-Jet" /*13*/,
      "+IP_{xy}" /*14*/
  };

  const int nBins2dSignificance = 400;
  const int nBins3dSignificance = 250;
  const int nBins2d = 500;
  const int nBins3d = 250;

  if (!fOutput)
    fOutput = new AliEmcalList();
  fOutput->SetName("fJHFTaggingOutput");
  fOutput->SetOwner(true);

  fhistInclusiveJetCuts = new TH1D("fhistInclusiveJetCuts", "Number of Inclusive jets after cuts", 15, 0, 15);
  if (fIsPythia)
  {
    fhistbJetCuts = new TH1D("fhistbJetCuts", "Number of b jets after cuts", 15, 0, 15);
    fhistcJetCuts = new TH1D("fhistcJetCuts", "Number of c jets after cuts", 15, 0, 15);
    fhistlfJetCuts = new TH1D("fhistlfJetCuts", "Number of lf jets after cuts", 15, 0, 15);
  }

  for (int j = 0; j < 15; j++)
  {
    fhistInclusiveJetCuts->GetXaxis()->SetBinLabel(j + 1, BjetCuts[j].Data());
    if (fIsPythia)
    {
      fhistbJetCuts->GetXaxis()->SetBinLabel(j + 1, BjetCuts[j].Data());
      fhistcJetCuts->GetXaxis()->SetBinLabel(j + 1, BjetCuts[j].Data());
      fhistlfJetCuts->GetXaxis()->SetBinLabel(j + 1, BjetCuts[j].Data());
    }
  }

  fh1dJetRecEtaPhiAccepted = new TH2D("fh1dJetRecEtaPhiAccepted", "detector level jet;#eta;phi", 200, -1.0, 1.0, 200, 0., TMath::TwoPi());

  fh1dJetRecPtAcceptedunCorr = new TH1D("fh1dJetRecPtAcceptedunCorr", "Rec Jet Pt uncorrected;#it{p}_{T,jet} (GeV/#it{c})", 250, 0, 250);

  f2histRhoVsDeltaPt = new TH2D("f2histRhoVsDeltaPt", "Rho Vs Delta Pt;#delta P_{T}^{RC} (Gev/c);#rho (Gev/c)", 170, -20, 150, 30, 0, 30);
  f2histRhoVsDeltaPtFirst = new TH2D("f2histRhoVsDeltaPtFirst", "Rho Vs Delta Pt N=1 tagged Events;#delta P_{T}^{RC} (Gev/c);#rho (Gev/c)", 170, -20, 150, 30, 0, 30);
  f2histRhoVsDeltaPtSecond = new TH2D("f2histRhoVsDeltaPtSecond", "Rho Vs Delta Pt N=2 tagged Events;#delta P_{T}^{RC} (Gev/c);#rho (Gev/c)", 170, -20, 150, 30, 0, 30);
  f2histRhoVsDeltaPtThird = new TH2D("f2histRhoVsDeltaPtThird", "Rho Vs Delta Pt N=3 tagged Events;#delta P_{T}^{RC} (Gev/c);#rho (Gev/c)", 170, -20, 150, 30, 0, 30);

  if (fEmbeddPerpendicular)
  {
    fh2DeltaPtEmbeddCorrelationPerpendicular = new TH2D("fh2DeltaPtEmbeddCorrelationPerpendicular", "Rho Vs Delta embedding;#delta P_{T}^{RC} (Gev/c);#rho (Gev/c)", 170, -20, 150, 30, 0, 30);
  }

  if (fIsPythia)
  {
    fh1dJetGenPt = new TH1D("fh1dJetGenPt", "generator level jets;#it{p}_{T} (GeV/#it{c}); count", 500, 0, 250);
    fh1dJetGenPtUnidentified = new TH1D("fh1dJetGenPtUnidentified", "generator level jets (no flavour assigned);#it{p}_{T} (GeV/#it{c}); count", 500, 0, 250);
    fh1dJetGenPtudsg = new TH1D("fh1dJetGenPtudsg", "generator level udsg jets;#it{p}_{T} (GeV/#it{c}); count", 500, 0, 250);
    fh1dJetGenPtc = new TH1D("fh1dJetGenPtc", "generator level c jets;#it{p}_{T} (GeV/#it{c}); count", 500, 0, 250);
    fh1dJetGenPtb = new TH1D("fh1dJetGenPtb", "generator level b jets;#it{p}_{T} (GeV/#it{c}); count", 500, 0, 250);

    fh2dJetGenPtVsJetRecPt = new TH2D("fh2dJetGenPtVsJetRecPt", "detector momentum response;rec pt;gen pt", 500, 0, 250, 500, 0, 250);
    fh2dJetGenPtVsJetRecPtb = new TH2D("fh2dJetGenPtVsJetRecPtb", "detector momentum response;rec pt;gen pt", 500, 0, 250, 500, 0, 250);
    fh2dJetGenPtVsJetRecPtc = new TH2D("fh2dJetGenPtVsJetRecPtc", "detector momentum response;rec pt;gen pt", 500, 0, 250, 500, 0, 250);
    fh2dJetGenPtVsJetRecPtudsg = new TH2D("fh2dJetGenPtVsJetRecPtudsg", "detector momentum response;rec pt;gen pt", 500, 0, 250, 500, 0, 250);

    // Track Counting Analysis
    if (fDoTrackCountingAnalysis)
    {
      fh2dJetSignedImpParXYUnidentified = new TH2D("fh2dJetSignedImpParXYUnidentified", "fh2dJetSignedImpParXYUnidentified;#it{p}_{T,jet} (GeV/#it{c}); 2D Impact Paramter (cm);a.u.", 500, 0., 250, nBins2d, -1, 1);
      fh2dJetSignedImpParXYZUnidentified = new TH2D("fh2dJetSignedImpParXYZUnidentified", "fh2dJetSignedImpParXYZUnidentified;#it{p}_{T,jet} (GeV/#it{c}); 3d imp. parameter (cm);a.u.", 500, 0., 250, nBins3d, -1, 1);
      fh2dJetSignedImpParXYSignificanceUnidentified = new TH2D("fh2dJetSignedImpParXYSignificanceUnidentified", "fh2dJetSignedImpParXYSignificanceUnidentified;#it{p}_{T,jet} (GeV/#it{c}); 2D Impact Paramter significance;a.u.", 500, 0., 250, nBins2dSignificance, -100, 100);
      fh2dJetSignedImpParXYZSignificanceUnidentified = new TH2D("fh2dJetSignedImpParXYZSignificanceUnidentified", "fh2dJetSignedImpParXYZSignificanceUnidentified;#it{p}_{T,jet} (GeV/#it{c}); 3d imp. parameter significance;a.u.", 500, 0., 250, nBins3dSignificance, -100, 100);

      fh2dJetSignedImpParXYudsg = new TH2D("fh2dJetSignedImpParXYudsg", "fh2dJetSignedImpParXYudsg;#it{p}_{T,jet} (GeV/#it{c}); 2D Impact Paramter (cm);a.u.", 500, 0., 250, nBins2d, -1, 1);
      fh2dJetSignedImpParXYZudsg = new TH2D("fh2dJetSignedImpParXYZudsg", "fh2dJetSignedImpParXYZudsg;#it{p}_{T,jet} (GeV/#it{c}); 3d imp. parameter (cm);a.u.", 500, 0., 250, nBins3d, -1, 1);
      fh2dJetSignedImpParXYSignificanceudsg = new TH2D("fh2dJetSignedImpParXYSignificanceudsg", "fh2dJetSignedImpParXYSignificanceudsg;#it{p}_{T,jet} (GeV/#it{c}); 2D Impact Paramter significance;a.u.", 500, 0., 250, nBins2dSignificance, -100, 100);
      fh2dJetSignedImpParXYZSignificanceudsg = new TH2D("fh2dJetSignedImpParXYZSignificanceudsg", "fh2dJetSignedImpParXYZSignificanceudsg;#it{p}_{T,jet} (GeV/#it{c}); 3d imp. parameter significance;a.u.", 500, 0., 250, nBins3dSignificance, -100, 100);

      fh2dJetSignedImpParXYc = new TH2D("fh2dJetSignedImpParXYc", "fh2dJetSignedImpParXYc;#it{p}_{T,jet} (GeV/#it{c}); 2D Impact Paramter (cm);a.u.", 500, 0., 250, nBins2d, -1, 1);
      fh2dJetSignedImpParXYZc = new TH2D("fh2dJetSignedImpParXYZc", "fh2dJetSignedImpParXYZc;#it{p}_{T,jet} (GeV/#it{c}); 3d imp. parameter (cm);a.u.", 500, 0., 250, nBins3d, -1, 1);
      fh2dJetSignedImpParXYSignificancec = new TH2D("fh2dJetSignedImpParXYSignificancec", "fh2dJetSignedImpParXYSignificancec;#it{p}_{T,jet} (GeV/#it{c}); 2D Impact Paramter significance;a.u.", 500, 0., 250, nBins2dSignificance, -100, 100);
      fh2dJetSignedImpParXYZSignificancec = new TH2D("fh2dJetSignedImpParXYZSignificancec", "fh2dJetSignedImpParXYZSignificancec;#it{p}_{T,jet} (GeV/#it{c}); 3d imp. parameter significance;a.u.", 500, 0., 250, nBins3dSignificance, -100, 100);

      fh2dJetSignedImpParXYb = new TH2D("fh2dJetSignedImpParXYb", "fh2dJetSignedImpParXYb;#it{p}_{T,jet} (GeV/#it{c}); 2D Impact Paramter (cm);a.u.", 500, 0., 250, nBins2d, -1, 1);
      fh2dJetSignedImpParXYZb = new TH2D("fh2dJetSignedImpParXYZb", "fh2dJetSignedImpParXYZb;#it{p}_{T,jet} (GeV/#it{c}); 3d imp. parameter (cm);a.u.", 500, 0., 250, nBins3d, -1, 1);
      fh2dJetSignedImpParXYSignificanceb = new TH2D("fh2dJetSignedImpParXYSignificanceb", "fh2dJetSignedImpParXYSignificanceb;#it{p}_{T,jet} (GeV/#it{c}); 2D Impact Paramter significance;a.u.", 500, 0., 250, nBins2dSignificance, -100, 100);
      fh2dJetSignedImpParXYZSignificanceb = new TH2D("fh2dJetSignedImpParXYZSignificanceb", "fh2dJetSignedImpParXYZSignificanceb;#it{p}_{T,jet} (GeV/#it{c}); 3d imp. parameter significance;a.u.", 500, 0., 250, nBins3dSignificance, -100, 100);
      // N=1
      fh2dJetSignedImpParXYSignificanceUnidentifiedFirst = new TH2D("fh2dJetSignedImpParXYSignificanceUnidentifiedFirst", "fh2dJetSignedImpParXYSignificanceUnidentifiedFirst;#it{p}_{T,jet} (GeV/#it{c}); 2D Impact Paramter significance;a.u.", 500, 0., 250, nBins2dSignificance, -100, 100);
      fh2dJetSignedImpParXYSignificanceudsgFirst = new TH2D("fh2dJetSignedImpParXYSignificanceudsgFirst", "fh2dJetSignedImpParXYSignificanceudsgFirst;#it{p}_{T,jet} (GeV/#it{c}); 2D Impact Paramter significance;a.u.", 500, 0., 250, nBins2dSignificance, -100, 100);
      fh2dJetSignedImpParXYSignificancecFirst = new TH2D("fh2dJetSignedImpParXYSignificancecFirst", "fh2dJetSignedImpParXYSignificancecFirst;#it{p}_{T,jet} (GeV/#it{c}); 2D Impact Paramter significance;a.u.", 500, 0., 250, nBins2dSignificance, -100, 100);
      fh2dJetSignedImpParXYSignificancebFirst = new TH2D("fh2dJetSignedImpParXYSignificancebFirst", "fh2dJetSignedImpParXYSignificancebFirst;#it{p}_{T,jet} (GeV/#it{c}); 2D Impact Paramter significance;a.u.", 500, 0., 250, nBins2dSignificance, -100, 100);

      // N=2
      fh2dJetSignedImpParXYSignificanceUnidentifiedSecond = new TH2D("fh2dJetSignedImpParXYSignificanceUnidentifiedSecond", "fh2dJetSignedImpParXYSignificanceUnidentifiedSecond;#it{p}_{T,jet} (GeV/#it{c}); 2D Impact Paramter significance;a.u.", 500, 0., 250, nBins2dSignificance, -100, 100);
      fh2dJetSignedImpParXYSignificanceudsgSecond = new TH2D("fh2dJetSignedImpParXYSignificanceudsgSecond", "fh2dJetSignedImpParXYSignificanceudsgSecond;#it{p}_{T,jet} (GeV/#it{c}); 2D Impact Paramter significance;a.u.", 500, 0., 250, nBins2dSignificance, -100, 100);
      fh2dJetSignedImpParXYSignificancecSecond = new TH2D("fh2dJetSignedImpParXYSignificancecSecond", "fh2dJetSignedImpParXYSignificancecSecond;#it{p}_{T,jet} (GeV/#it{c}); 2D Impact Paramter significance;a.u.", 500, 0., 250, nBins2dSignificance, -100, 100);
      fh2dJetSignedImpParXYSignificancebSecond = new TH2D("fh2dJetSignedImpParXYSignificancebSecond", "fh2dJetSignedImpParXYSignificancebSecond;#it{p}_{T,jet} (GeV/#it{c}); 2D Impact Paramter significance;a.u.", 500, 0., 250, nBins2dSignificance, -100, 100);
      // N=3
      fh2dJetSignedImpParXYSignificanceUnidentifiedThird = new TH2D("fh2dJetSignedImpParXYSignificanceUnidentifiedThird", "fh2dJetSignedImpParXYSignificanceUnidentifiedThird;#it{p}_{T,jet} (GeV/#it{c}); 2D Impact Paramter significance;a.u.", 500, 0., 250, nBins2dSignificance, -100, 100);
      fh2dJetSignedImpParXYSignificanceudsgThird = new TH2D("fh2dJetSignedImpParXYSignificanceudsgThird", "fh2dJetSignedImpParXYSignificanceudsgThird;#it{p}_{T,jet} (GeV/#it{c}); 2D Impact Paramter significance;a.u.", 500, 0., 250, nBins2dSignificance, -100, 100);
      fh2dJetSignedImpParXYSignificancecThird = new TH2D("fh2dJetSignedImpParXYSignificancecThird", "fh2dJetSignedImpParXYSignificancecThird;#it{p}_{T,jet} (GeV/#it{c}); 2D Impact Paramter significance;a.u.", 500, 0., 250, nBins2dSignificance, -100, 100);
      fh2dJetSignedImpParXYSignificancebThird = new TH2D("fh2dJetSignedImpParXYSignificancebThird", "fh2dJetSignedImpParXYSignificancebThird;#it{p}_{T,jet} (GeV/#it{c}); 2D Impact Paramter significance;a.u.", 500, 0., 250, nBins2dSignificance, -100, 100);
    }
  }
  // Jet histograms
  fh1dJetRecPt = new TH1D("fh1dJetRecPt", "detector level jets;#it{p}_{T} (GeV/#it{c}); count", 500, 0, 250);
  fh1dJetRecPtAccepted = new TH1D("fh1dJetRecPtAccepted", "accepted detector level jets;#it{p}_{T} (GeV/#it{c}); count", 500, 0, 250);

  if (fIsPythia)
  {
    fh1dJetRecPtUnidentified = new TH1D("fh1dJetRecPtUnidentified", "detector level jets;#it{p}_{T} (GeV/#it{c}); count", 500, 0, 250);
    fh1dJetRecPtudsg = new TH1D("fh1dJetRecPtudsg", "detector level jets;#it{p}_{T} (GeV/#it{c}); count", 500, 0, 250);
    fh1dJetRecPtc = new TH1D("fh1dJetRecPtc", "detector level jets;#it{p}_{T} (GeV/#it{c}); count", 500, 0, 250);
    fh1dJetRecPtb = new TH1D("fh1dJetRecPtb", "detector level jets;#it{p}_{T} (GeV/#it{c}); count", 500, 0, 250);
    fh1dJetRecPtUnidentifiedAccepted = new TH1D("fh1dJetRecPtUnidentifiedAccepted", "detector level jets;#it{p}_{T} (GeV/#it{c}); count", 500, 0, 250);
    fh1dJetRecPtudsgAccepted = new TH1D("fh1dJetRecPtudsgAccepted", "detector level jets;#it{p}_{T} (GeV/#it{c}); count", 500, 0, 250);
    fh1dJetRecPtcAccepted = new TH1D("fh1dJetRecPtcAccepted", "detector level jets;#it{p}_{T} (GeV/#it{c}); count", 500, 0, 250);
    fh1dJetRecPtbAccepted = new TH1D("fh1dJetRecPtbAccepted", "detector level jets;#it{p}_{T} (GeV/#it{c}); count", 500, 0, 250);
  }

  if (fDoTrackCountingAnalysis)
  {
    fh2dJetSignedImpParXY = new TH2D("fh2dJetSignedImpParXY", "fh2dJetSignedImpParXY;#it{p}_{T,jet} (GeV/#it{c}); 2D Impact Paramter (cm);a.u.", 500, 0., 250, nBins2d, -1, 1);
    fh2dJetSignedImpParXYZ = new TH2D("fh2dJetSignedImpParXYZ", "fh2dJetSignedImpParXYZ;#it{p}_{T,jet} (GeV/#it{c}); 3d imp. parameter (cm);a.u.", 500, 0., 250, nBins3d, -1, 1);
    fh2dJetSignedImpParXYSignificance = new TH2D("fh2dJetSignedImpParXYSignificance", "fh2dJetSignedImpParXYSignificance;#it{p}_{T,jet} (GeV/#it{c}); 2D Impact Paramter significance;a.u.", 500, 0., 250, nBins2dSignificance, -100, 100);
    fh2dJetSignedImpParXYZSignificance = new TH2D("fh2dJetSignedImpParXYZSignificance", "fh2dJetSignedImpParXYZSignificance;#it{p}_{T,jet} (GeV/#it{c}); 3d imp. parameter significance;a.u.", 500, 0., 250, nBins3dSignificance, -100, 100);

    // N=1
    fh2dJetSignedImpParXYSignificanceFirst = new TH2D("fh2dJetSignedImpParXYSignificanceFirst", "fh2dJetSignedImpParXYSignificanceFirst;#it{p}_{T,jet} (GeV/#it{c}); 2D Impact Paramter significance;a.u.", 500, 0., 250, nBins2dSignificance, -100, 100);
    // N=2
    fh2dJetSignedImpParXYSignificanceSecond = new TH2D("fh2dJetSignedImpParXYSignificanceSecond", "fh2dJetSignedImpParXYSignificanceSecond;#it{p}_{T,jet} (GeV/#it{c}); 2D Impact Paramter significance;a.u.", 500, 0., 250, nBins2dSignificance, -100, 100);
    // N=3
    fh2dJetSignedImpParXYSignificanceThird = new TH2D("fh2dJetSignedImpParXYSignificanceThird", "fh2dJetSignedImpParXYSignificanceThird;#it{p}_{T,jet} (GeV/#it{c}); 2D Impact Paramter significance;a.u.", 500, 0., 250, nBins2dSignificance, -100, 100);
  }

  if (fDoSVAnalysis)
  {
    fHistSV2Prong = new TH3D("fHistSV2Prong", "Secondary vertex 2Prong;#it{p}_{T,jet} (GeV/#it{c});L_{xy}/#sigma;M_{Vtx}(GeV/c^2)", 500, 0., 250., 160, 0, 80, 100, 0, 10);
    fHistSV3Prong = new TH3D("fHistSV3Prong", "Secondary vertex 3Prong;#it{p}_{T,jet} (GeV/#it{c});L_{xy}/#sigma;M_{Vtx}(GeV/c^2)", 500, 0., 250., 160, 0, 80, 100, 0, 10);

    fHistDispersion2Prong = new TH2D("fHistDispersion2Prong", "Secondary vertex Dispersion 2Prong;#it{p}_{T,jet} (GeV/#it{c});Vtx Dispersion", 80, 0., 250., 100, 0, 0.5);
    fHistDispersion3Prong = new TH2D("fHistDispersion3Prong", "Secondary vertex Dispersion 3Prong;#it{p}_{T,jet} (GeV/#it{c});Vtx Dispersion", 80, 0., 250., 100, 0, 0.5);

    if (fIsPythia)
    {

      fHistSV2ProngUnidentified = new TH3D("fHistSV2ProngUnidentified", "Secondary vertex 2Prong Unidentified;#it{p}_{T,jet} (GeV/#it{c});L_{xy}/#sigma;M_{Vtx}(GeV/c^2)", 500, 0., 250., 400, 0, 80, 100, 0, 10);
      fHistSV3ProngUnidentified = new TH3D("fHistSV3ProngUnidentified", "Secondary vertex 3Prong Unidentified;#it{p}_{T,jet} (GeV/#it{c});L_{xy}/#sigma;M_{Vtx}(GeV/c^2)", 500, 0., 250., 400, 0, 80, 100, 0, 10);

      fHistSV2Prongb = new TH3D("fHistSV2Prongb", "Secondary vertex 2Prong b-jet;#it{p}_{T,jet} (GeV/#it{c});L_{xy}/#sigma;M_{Vtx}(GeV/c^2)", 500, 0., 250., 160, 0, 80, 100, 0, 10);
      fHistSV3Prongb = new TH3D("fHistSV3Prongb", "Secondary vertex 3Prong b-jet;#it{p}_{T,jet} (GeV/#it{c});L_{xy}/#sigma;M_{Vtx}(GeV/c^2)", 500, 0., 250., 160, 0, 80, 100, 0, 10);

      fHistSV2Prongc = new TH3D("fHistSV2Prongc", "Secondary vertex 2Prong c-jet;#it{p}_{T,jet} (GeV/#it{c});L_{xy}/#sigma;M_{Vtx}(GeV/c^2)", 500, 0., 250., 160, 0, 80, 100, 0, 10);
      fHistSV3Prongc = new TH3D("fHistSV3Prongc", "Secondary vertex 3Prong c-jet;#it{p}_{T,jet} (GeV/#it{c});L_{xy}/#sigma;M_{Vtx}(GeV/c^2)", 500, 0., 250., 160, 0, 80, 100, 0, 10);

      fHistSV2Pronglf = new TH3D("fHistSV2Pronglf", "Secondary vertex 2Prong lf-jet;#it{p}_{T,jet} (GeV/#it{c});L_{xy}/#sigma;M_{Vtx}(GeV/c^2)", 500, 0., 250., 160, 0, 80, 100, 0, 10);
      fHistSV3Pronglf = new TH3D("fHistSV3Pronglf", "Secondary vertex 3Prong lf-jet;#it{p}_{T,jet} (GeV/#it{c});L_{xy}/#sigma;M_{Vtx}(GeV/c^2)", 500, 0., 250., 160, 0, 80, 100, 0, 10);

      fHistDispersion2ProngUnidentified = new TH2D("fHistDispersion2ProngUnidentified", "Secondary vertex Dispersion 2Prong Unidentified;#it{p}_{T,jet} (GeV/#it{c});Vtx Dispersion", 500, 0., 250., 100, 0, 0.5);
      fHistDispersion3ProngUnidentified = new TH2D("fHistDispersion3ProngUnidentified", "Secondary vertex Dispersion 3Prong Unidentified;#it{p}_{T,jet} (GeV/#it{c});Vtx Dispersion", 500, 0., 250., 100, 0, 0.5);

      fHistDispersion2Prongb = new TH2D("fHistDispersion2Prongb", "Secondary vertex Dispersion 2Prong b-jet;#it{p}_{T,jet} (GeV/#it{c});Vtx Dispersion", 500, 0., 250., 100, 0, 0.5);
      fHistDispersion3Prongb = new TH2D("fHistDispersion3Prongb", "Secondary vertex Dispersion 3Prong b-jet;#it{p}_{T,jet} (GeV/#it{c});Vtx Dispersion", 500, 0., 250., 100, 0, 0.5);

      fHistDispersion2Prongc = new TH2D("fHistDispersion2Prongc", "Secondary vertex Dispersion 2Prong c-jet;#it{p}_{T,jet} (GeV/#it{c});Vtx Dispersion", 500, 0., 250., 100, 0, 0.5);
      fHistDispersion3Prongc = new TH2D("fHistDispersion3Prongc", "Secondary vertex Dispersion 3Prong c-jet;#it{p}_{T,jet} (GeV/#it{c});Vtx Dispersion", 500, 0., 250., 100, 0, 0.5);

      fHistDispersion2Pronglf = new TH2D("fHistDispersion2Pronglf", "Secondary vertex Dispersion 2Prong lf-jet;#it{p}_{T,jet} (GeV/#it{c});Vtx Dispersion", 500, 0., 250., 100, 0, 0.5);
      fHistDispersion3Pronglf = new TH2D("fHistDispersion3Pronglf", "Secondary vertex Dispersion 3Prong lf-jet;#it{p}_{T,jet} (GeV/#it{c});Vtx Dispersion", 500, 0., 250., 100, 0, 0.5);
    }
  }

  // Add to output list
  fOutput->Add(fhistInclusiveJetCuts);
  if (fIsPythia)
  {
    fOutput->Add(fhistbJetCuts);
    fOutput->Add(fhistcJetCuts);
    fOutput->Add(fhistlfJetCuts);
  }
  //____

  fOutput->Add(fh1dJetRecPtAcceptedunCorr);

  fOutput->Add(f2histRhoVsDeltaPt);
  fOutput->Add(f2histRhoVsDeltaPtFirst);
  fOutput->Add(f2histRhoVsDeltaPtSecond);
  fOutput->Add(f2histRhoVsDeltaPtThird);

  if (fEmbeddPerpendicular)
  {
    fOutput->Add(fh2DeltaPtEmbeddCorrelationPerpendicular);
  }

  if (fIsPythia)
  {
    fOutput->Add(fh1dJetGenPt);
    fOutput->Add(fh1dJetGenPtUnidentified);
    fOutput->Add(fh1dJetGenPtudsg);
    fOutput->Add(fh1dJetGenPtc);
    fOutput->Add(fh1dJetGenPtb);
    fOutput->Add(fh2dJetGenPtVsJetRecPt);
    fOutput->Add(fh2dJetGenPtVsJetRecPtb);
    fOutput->Add(fh2dJetGenPtVsJetRecPtc);
    fOutput->Add(fh2dJetGenPtVsJetRecPtudsg);

    if (fDoTrackCountingAnalysis)
    {
      fOutput->Add(fh2dJetSignedImpParXYUnidentified);
      fOutput->Add(fh2dJetSignedImpParXYZUnidentified);
      fOutput->Add(fh2dJetSignedImpParXYSignificanceUnidentified);
      fOutput->Add(fh2dJetSignedImpParXYZSignificanceUnidentified);
      fOutput->Add(fh2dJetSignedImpParXYudsg);
      fOutput->Add(fh2dJetSignedImpParXYZudsg);
      fOutput->Add(fh2dJetSignedImpParXYSignificanceudsg);
      fOutput->Add(fh2dJetSignedImpParXYZSignificanceudsg);
      fOutput->Add(fh2dJetSignedImpParXYc);
      fOutput->Add(fh2dJetSignedImpParXYZc);
      fOutput->Add(fh2dJetSignedImpParXYSignificancec);
      fOutput->Add(fh2dJetSignedImpParXYZSignificancec);
      fOutput->Add(fh2dJetSignedImpParXYb);
      fOutput->Add(fh2dJetSignedImpParXYZb);
      fOutput->Add(fh2dJetSignedImpParXYSignificanceb);
      fOutput->Add(fh2dJetSignedImpParXYZSignificanceb);

      fOutput->Add(fh2dJetSignedImpParXYSignificanceUnidentifiedFirst);
      fOutput->Add(fh2dJetSignedImpParXYSignificanceudsgFirst);
      fOutput->Add(fh2dJetSignedImpParXYSignificancecFirst);
      fOutput->Add(fh2dJetSignedImpParXYSignificancebFirst);

      fOutput->Add(fh2dJetSignedImpParXYSignificanceUnidentifiedSecond);
      fOutput->Add(fh2dJetSignedImpParXYSignificanceudsgSecond);
      fOutput->Add(fh2dJetSignedImpParXYSignificancecSecond);
      fOutput->Add(fh2dJetSignedImpParXYSignificancebSecond);

      fOutput->Add(fh2dJetSignedImpParXYSignificanceUnidentifiedThird);
      fOutput->Add(fh2dJetSignedImpParXYSignificanceudsgThird);
      fOutput->Add(fh2dJetSignedImpParXYSignificancecThird);
      fOutput->Add(fh2dJetSignedImpParXYSignificancebThird);
    }
  }
  fOutput->Add(fh1dJetRecPt);
  fOutput->Add(fh1dJetRecPtAccepted);
  fOutput->Add(fh1dJetRecEtaPhiAccepted);

  if (fIsPythia)
  {

    fOutput->Add(fh1dJetRecPtUnidentified);
    fOutput->Add(fh1dJetRecPtudsg);
    fOutput->Add(fh1dJetRecPtc);
    fOutput->Add(fh1dJetRecPtb);
    fOutput->Add(fh1dJetRecPtUnidentifiedAccepted);
    fOutput->Add(fh1dJetRecPtudsgAccepted);
    fOutput->Add(fh1dJetRecPtcAccepted);
    fOutput->Add(fh1dJetRecPtbAccepted);
  }

  if (fDoTrackCountingAnalysis)
  {
    fOutput->Add(fh2dJetSignedImpParXY);
    fOutput->Add(fh2dJetSignedImpParXYZ);
    fOutput->Add(fh2dJetSignedImpParXYSignificance);
    fOutput->Add(fh2dJetSignedImpParXYZSignificance);

    fOutput->Add(fh2dJetSignedImpParXYSignificanceFirst);
    fOutput->Add(fh2dJetSignedImpParXYSignificanceSecond);
    fOutput->Add(fh2dJetSignedImpParXYSignificanceThird);
  }

  if (fDoSVAnalysis)
  {
    fOutput->Add(fHistDispersion2Prong);
    fOutput->Add(fHistSV2Prong);

    if (fIsPythia)
    {
      fOutput->Add(fHistDispersion2ProngUnidentified);
      fOutput->Add(fHistDispersion2Prongb);
      fOutput->Add(fHistDispersion2Prongc);
      fOutput->Add(fHistDispersion2Pronglf);

      fOutput->Add(fHistSV2ProngUnidentified);
      fOutput->Add(fHistSV2Prongb);
      fOutput->Add(fHistSV2Prongc);
      fOutput->Add(fHistSV2Pronglf);
    }

    fOutput->Add(fHistDispersion3Prong);
    fOutput->Add(fHistSV3Prong);

    if (fIsPythia)
    {
      fOutput->Add(fHistDispersion3ProngUnidentified);
      fOutput->Add(fHistDispersion3Prongb);
      fOutput->Add(fHistDispersion3Prongc);
      fOutput->Add(fHistDispersion3Pronglf);

      fOutput->Add(fHistSV3ProngUnidentified);
      fOutput->Add(fHistSV3Prongb);
      fOutput->Add(fHistSV3Prongc);
      fOutput->Add(fHistSV3Pronglf);
    }
  }

  TIter next(fOutput);
  while (TObject *obj = next.Next())
  {
    if (obj->IsA() == TH1D::Class() || obj->IsA() == TH2D::Class())
    {
      ((TH1 *)obj)->Sumw2();
    }
  }
  PostData(1, fOutput); // Post data for ALL output slots > 0 here.
}

//______________________________________________________________________________
Bool_t AliJHFTagging::Run()
{
  if (fDoTrackCountingAnalysis && fDoSVAnalysis)
  {
    AliFatal("Cannot use two tagging methods. Please choose only one");
  }

  // Processing of one event
  if (!((Entry() - 1) % 100))
    AliInfo(Form(" Processing event # %lld", Entry()));

  fAODIn = NULL;

  fAODIn = dynamic_cast<AliAODEvent *>(InputEvent());
  fPrimaryVertex = dynamic_cast<AliAODVertex *>(fAODIn->GetPrimaryVertex());

  // initializing variables from last event
  fInputList->Clear();
  if (fIsPythia)
  {
    fInputListb->Clear();
    fInputListc->Clear();
    fInputListlf->Clear();
  }

  fEvtNum++;
  if (fEvtNum % 1000 == 0)
    std::cout << "evt : " << fEvtNum << ", fEntry: " << fEntry << std::endl;

  fMCArray = NULL;

  if (fIsPythia)
  {
    fJetContainerMC = static_cast<AliJetContainer *>(fJetCollArray.At(1));
    fMCArray = dynamic_cast<TClonesArray *>(fAODIn->FindListObject(AliAODMCParticle::StdBranchName()));
  }

  fVertexer = new AliVertexerTracks(fAODIn->GetMagneticField());
  fVertexer->SetITSMode();
  fVertexer->SetMinClusters(3);
  fVertexer->SetConstraintOn();

  const AliVVertex *trkVtx = dynamic_cast<const AliVVertex *>(fAODIn->GetPrimaryVertex());
  TString vtxTtl = trkVtx->GetTitle();
  if (vtxTtl.Contains("WithConstraint"))
  {
    Float_t diamondcovxy[3];
    fAODIn->GetDiamondCovXY(diamondcovxy);
    double pos[3] = {fAODIn->GetDiamondX(), fAODIn->GetDiamondY(), 0.};
    double cov[6] = {diamondcovxy[0], diamondcovxy[1], diamondcovxy[2], 0., 0., 10. * 10.};
    fDiamond = new AliESDVertex(pos, cov, 1., 1);
    fVertexer->SetVtxStart(fDiamond);
  }

  // Main part jet analysis
  // preparation
  fJetContainerData = static_cast<AliJetContainer *>(fJetCollArray.At(0));

  double randomConePt(0.);
  if (fFillControlHists)
  {
    randomConePt = GetDeltaPtRandomCone();
    f2histRhoVsDeltaPt->Fill(randomConePt, fJetContainerData->GetRhoVal());
  }

  if (fIsPythia)
  {
    // SetContainer
    AliEmcalJet *jetgen = 0x0;
    AliAODMCParticle *partonAOD = NULL;

    if (!MatchJetsGeometricDefault())
      std::cout << "Error running jet matching!" << std::endl;

    fJetContainerMC->ResetCurrentID();

    // Fill gen. level jet histograms
    while ((jetgen = fJetContainerMC->GetNextAcceptJet()))
    {
      if (!jetgen)
        continue;
      int MCJetflavour = 0;
      int partonpdg = 0;

      if (fUsePartonDef)
      {
        partonAOD = fHFJetUtils->IsMCJetParton(fMCArray, jetgen, 0.4);
        if (!(partonAOD))
          MCJetflavour = 0;
        else
        {
          partonpdg = abs(partonAOD->PdgCode());
          if (partonpdg == 1 || partonpdg == 2 || partonpdg == 3 || partonpdg == 21)
            MCJetflavour = 1;
          else if (partonpdg == 4)
            MCJetflavour = 2;
          else if (partonpdg == 5)
            MCJetflavour = 3;
        }
      }
      else
      {
        partonAOD = fHFJetUtils->IsMCJetMeson(fMCArray, jetgen, 0.4);
        if (!(partonAOD))
          MCJetflavour = 0;
        else
        {
          partonpdg = abs(partonAOD->PdgCode());
          if (fHFJetUtils->IsBMeson(partonpdg))
            MCJetflavour = 3;
          else if (fHFJetUtils->IsDMeson(partonpdg))
            MCJetflavour = 2;
          else
            MCJetflavour = 1;
        }
      }

      double genpt = jetgen->Pt();
      if (!(fJetContainerMC->GetRhoParameter() == 0x0))
      {
        genpt = genpt - fJetContainerMC->GetRhoVal() * jetgen->Area();
      }

      if (fFillControlHists)
      {
        fh1dJetGenPt->Fill(genpt);
        if (MCJetflavour == 0)
          fh1dJetGenPtUnidentified->Fill(genpt);
        else if (MCJetflavour == 1)
          fh1dJetGenPtudsg->Fill(genpt);
        else if (MCJetflavour == 2)
          fh1dJetGenPtc->Fill(genpt);
        else if (MCJetflavour == 3)
          fh1dJetGenPtb->Fill(genpt);
      }
    }
  }

  // loop rec level jets
  AliEmcalJet *jetrec = 0x0;
  AliEmcalJet *jetmatched = 0x0;
  fJetContainerData->ResetCurrentID();
  double jetPt = 0;

  AliAODMCParticle *partonAOD = NULL;

  std::vector<double> sImpParXY, sImpParXYZ, sImpParXYSig, sImpParXYZSig;

  bool TaggedFirst(0), TaggedSecond(0), TaggedThird(0); // IPs Tagging
  bool taggedSV2Prong(0), taggedSV3Prong(0);

  bool firstJetFound = false;

  while ((jetrec = fJetContainerData->GetNextJet()))
  {
    if (!jetrec)
      continue;

    jetPt = jetrec->Pt();

    if (!(fJetContainerData->GetRhoParameter() == 0x0))
    {
      jetPt = jetPt - fJetContainerData->GetRhoVal() * jetrec->Area();
    }

    // make inclusive signed imp. parameter constituent histograms
    int ntracks = (int)jetrec->GetNumberOfTracks();

    double dca[2] = {-99999, -99999};
    double cov[3] = {-99999, -99999, -99999};
    double sign = 0;

    short jetFlavor = 0;
    int partonpdg = 0;
    if (fIsPythia)
    {
      jetmatched = 0x0;
      jetmatched = jetrec->MatchedJet();

      if (jetmatched)
      {
        if (fUsePartonDef)
        {
          partonAOD = fHFJetUtils->IsMCJetParton(fMCArray, jetmatched, 0.4);

          if ((!partonAOD))
            jetFlavor = 0;
          else
          {
            partonpdg = abs(partonAOD->PdgCode());
            if (partonpdg == 1 || partonpdg == 2 || partonpdg == 3 || partonpdg == 21)
              jetFlavor = 1;
            else if (partonpdg == 4)
              jetFlavor = 2;
            else if (partonpdg == 5)
              jetFlavor = 3;
          }
        }
        else
        {
          partonAOD = fHFJetUtils->IsMCJetMeson(fMCArray, jetmatched, 0.4);
          if (!(partonAOD))
            jetFlavor = 0;
          else
          {
            partonpdg = abs(partonAOD->PdgCode());
            if (fHFJetUtils->IsBMeson(partonpdg))
              jetFlavor = 3;
            else if (fHFJetUtils->IsDMeson(partonpdg))
              jetFlavor = 2;
            else
              jetFlavor = 1;
          }
        }
      }
    }
    //	Printf("%s:%i",__FUNCTION__,__LINE__);

    if (fFillControlHists)
    {
      fh1dJetRecPt->Fill(jetrec->Pt());
      if (fIsPythia)
      {
        if (jetFlavor == 0)
          fh1dJetRecPtUnidentified->Fill(jetPt);
        else if (jetFlavor == 1)
          fh1dJetRecPtudsg->Fill(jetPt);
        else if (jetFlavor == 2)
          fh1dJetRecPtc->Fill(jetPt);
        else if (jetFlavor == 3)
          fh1dJetRecPtb->Fill(jetPt);
      }
    }

    UInt_t rejectionReason = 0;
    if (!(fJetContainerData->AcceptJet(jetrec, rejectionReason)))
      continue;

    if (fFillControlHists)
    {
      fh1dJetRecEtaPhiAccepted->Fill(jetrec->Eta(), jetrec->Phi());
      fh1dJetRecPtAcceptedunCorr->Fill(jetrec->Pt());
      fh1dJetRecPtAccepted->Fill(jetPt);
      if (fIsPythia)
      {
        if (jetrec->MatchedJet())
        {
          double genpt = jetrec->MatchedJet()->Pt();
          if (!(fJetContainerMC->GetRhoParameter() == 0x0))
          {
            genpt = genpt - fJetContainerMC->GetRhoVal() * jetrec->MatchedJet()->Area();
          }
          fh2dJetGenPtVsJetRecPt->Fill(jetPt, genpt);

          if (jetFlavor == 0)
          {
            fh1dJetRecPtUnidentifiedAccepted->Fill(jetPt);
          }
          else if (jetFlavor == 1)
          {
            fh1dJetRecPtudsgAccepted->Fill(jetPt);
            fh2dJetGenPtVsJetRecPtudsg->Fill(jetPt, genpt);
          }
          else if (jetFlavor == 2)
          {
            fh1dJetRecPtcAccepted->Fill(jetPt);
            fh2dJetGenPtVsJetRecPtc->Fill(jetPt, genpt);
          }
          else if (jetFlavor == 3)
          {
            fh1dJetRecPtbAccepted->Fill(jetPt);
            fh2dJetGenPtVsJetRecPtb->Fill(jetPt, genpt);
          }
        }
      }
    }

    if (fDoTrackCountingAnalysis)
    {
      for (int itrack = 0; itrack < ntracks; ++itrack)
      {
        double dcatrackjet = 999;
        double lineardecaylenth = 999;

        AliAODTrack *trackAOD = (AliAODTrack *)((fJetContainerData->GetParticleContainer())->GetParticle(jetrec->TrackAt(itrack)));

        if (!trackAOD)
          continue;

        if (!IsTrackAcceptedBJetCuts(trackAOD, jetFlavor))
          continue;

        if (!CalculateJetSignedTrackImpactParameter(trackAOD, jetrec, dca, cov, sign, dcatrackjet, lineardecaylenth))
          continue;

        // Select only tracks with dca rphi <1 cm and dca z < 5 cm
        //  linear decay length < 5 cm
        //  dca track to jet < 0.07 cm

        if (abs(dca[0]) > fTCMaxIPxy)
          continue;
        FillCandidateJet(10, jetFlavor);
        if (abs(dca[1]) > fTCMaxIPz)
          continue;
        FillCandidateJet(11, jetFlavor);
        if (lineardecaylenth > fTCMaxDecayLength)
          continue;
        FillCandidateJet(12, jetFlavor);
        if (dcatrackjet > fTCMaxDCATrackJet)
          continue;
        FillCandidateJet(13, jetFlavor);

        double cursImParXY = TMath::Abs(GetValImpactParameter(kXY, dca, cov)) * sign;
        double cursImParXYZ = TMath::Abs(GetValImpactParameter(kXYZ, dca, cov)) * sign;
        double cursImParXYSig = TMath::Abs(GetValImpactParameter(kXYSig, dca, cov)) * sign;
        double cursImParXYZSig = TMath::Abs(GetValImpactParameter(kXYZSig, dca, cov)) * sign;

        if (cursImParXYSig > fThresholdIP)
          FillCandidateJet(14, jetFlavor);

        if (!fDoTrackCountingAnalysis)
          continue;

        if (fFillControlHists)
        {
          fh2dJetSignedImpParXY->Fill(jetPt, cursImParXY);
          fh2dJetSignedImpParXYZ->Fill(jetPt, cursImParXYZ);
          fh2dJetSignedImpParXYSignificance->Fill(jetPt, cursImParXYSig);
          fh2dJetSignedImpParXYZSignificance->Fill(jetPt, cursImParXYZSig);

          if (fIsPythia)
          {
            if (jetFlavor == 0)
            {
              fh2dJetSignedImpParXYUnidentified->Fill(jetPt, cursImParXY);
              fh2dJetSignedImpParXYZUnidentified->Fill(jetPt, cursImParXYZ);
              fh2dJetSignedImpParXYSignificanceUnidentified->Fill(jetPt, cursImParXYSig);
              fh2dJetSignedImpParXYZSignificanceUnidentified->Fill(jetPt, cursImParXYZSig);
            }
            else if (jetFlavor == 1)
            {
              fh2dJetSignedImpParXYudsg->Fill(jetPt, cursImParXY);
              fh2dJetSignedImpParXYZudsg->Fill(jetPt, cursImParXYZ);
              fh2dJetSignedImpParXYSignificanceudsg->Fill(jetPt, cursImParXYSig);
              fh2dJetSignedImpParXYZSignificanceudsg->Fill(jetPt, cursImParXYZSig);
            }
            else if (jetFlavor == 2)
            {
              fh2dJetSignedImpParXYc->Fill(jetPt, cursImParXY);
              fh2dJetSignedImpParXYZc->Fill(jetPt, cursImParXYZ);
              fh2dJetSignedImpParXYSignificancec->Fill(jetPt, cursImParXYSig);
              fh2dJetSignedImpParXYZSignificancec->Fill(jetPt, cursImParXYZSig);
            }
            else if (jetFlavor == 3)
            {
              fh2dJetSignedImpParXYb->Fill(jetPt, cursImParXY);
              fh2dJetSignedImpParXYZb->Fill(jetPt, cursImParXYZ);
              fh2dJetSignedImpParXYSignificanceb->Fill(jetPt, cursImParXYSig);
              fh2dJetSignedImpParXYZSignificanceb->Fill(jetPt, cursImParXYZSig);
            }
          }
        }

        sImpParXY.push_back(cursImParXY);
        sImpParXYZ.push_back(cursImParXYZ);
        sImpParXYSig.push_back(cursImParXYSig);
        sImpParXYZSig.push_back(cursImParXYZSig);
      } // end of track loop

      // Ordered n=1,2,3 sip
      std::sort(sImpParXY.begin(), sImpParXY.end(), std::greater<double>());
      std::sort(sImpParXYZ.begin(), sImpParXYZ.end(), std::greater<double>());
      std::sort(sImpParXYSig.begin(), sImpParXYSig.end(), std::greater<double>());
      std::sort(sImpParXYZSig.begin(), sImpParXYZSig.end(), std::greater<double>());

      std::vector<double> DefaultDiscriminator = sImpParXYSig;

      // First largest
      if (DefaultDiscriminator.size() > 0)
      {
        if (DefaultDiscriminator.at(0) > fThresholdIP)
          TaggedFirst = true;

        if (fFillControlHists)
        {
          fh2dJetSignedImpParXYSignificanceFirst->Fill(jetPt, sImpParXYSig.at(0));

          if (fIsPythia)
          {
            if (jetFlavor == 0)
            {
              fh2dJetSignedImpParXYSignificanceUnidentifiedFirst->Fill(jetPt, sImpParXYSig.at(0));
            }
            else if (jetFlavor == 1)
            {
              fh2dJetSignedImpParXYSignificanceudsgFirst->Fill(jetPt, sImpParXYSig.at(0));
            }
            else if (jetFlavor == 2)
            {
              fh2dJetSignedImpParXYSignificancecFirst->Fill(jetPt, sImpParXYSig.at(0));
            }
            else if (jetFlavor == 3)
            {
              fh2dJetSignedImpParXYSignificancebFirst->Fill(jetPt, sImpParXYSig.at(0));
            }
          }
        }
      } // N=1

      // Second largest
      if (DefaultDiscriminator.size() > 1)
      {
        if (DefaultDiscriminator.at(1) > fThresholdIP)
        {
          TaggedSecond = true;
          if (!firstJetFound && fEmbeddPerpendicular && fFillControlHists)
          {
            firstJetFound = true;
            double EmbeddingPt = GetDeltaPtPerpEmbedding(jetrec->Eta(), jetrec->Phi());
            fh2DeltaPtEmbeddCorrelationPerpendicular->Fill(EmbeddingPt, fJetContainerData->GetRhoVal());
          }
        }

        if (fFillControlHists)
        {
          fh2dJetSignedImpParXYSignificanceSecond->Fill(jetPt, sImpParXYSig.at(1));

          if (fIsPythia)
          {
            if (jetFlavor == 0)
            {
              fh2dJetSignedImpParXYSignificanceUnidentifiedSecond->Fill(jetPt, sImpParXYSig.at(1));
            }
            else if (jetFlavor == 1)
            {
              fh2dJetSignedImpParXYSignificanceudsgSecond->Fill(jetPt, sImpParXYSig.at(1));
            }
            else if (jetFlavor == 2)
            {
              fh2dJetSignedImpParXYSignificancecSecond->Fill(jetPt, sImpParXYSig.at(1));
            }
            else if (jetFlavor == 3)
            {
              fh2dJetSignedImpParXYSignificancebSecond->Fill(jetPt, sImpParXYSig.at(1));
            }
          }
        }
      } // N=2

      // Third largest
      if (DefaultDiscriminator.size() > 2)
      {
        if (DefaultDiscriminator.at(2) > fThresholdIP)
          TaggedThird = true;

        if (fFillControlHists)
        {
          fh2dJetSignedImpParXYSignificanceThird->Fill(jetPt, sImpParXYSig.at(2));

          if (fIsPythia)
          {
            if (jetFlavor == 0)
            {
              fh2dJetSignedImpParXYSignificanceUnidentifiedThird->Fill(jetPt, sImpParXYSig.at(2));
            }
            else if (jetFlavor == 1)
            {
              fh2dJetSignedImpParXYSignificanceudsgThird->Fill(jetPt, sImpParXYSig.at(2));
            }
            else if (jetFlavor == 2)
            {
              fh2dJetSignedImpParXYSignificancecThird->Fill(jetPt, sImpParXYSig.at(2));
            }
            else if (jetFlavor == 3)
            {
              fh2dJetSignedImpParXYSignificancebThird->Fill(jetPt, sImpParXYSig.at(2));
            }
          }
        }
      } // N=3

      sImpParXY.clear();
      sImpParXYZ.clear();
      sImpParXYSig.clear();
      sImpParXYZSig.clear();
      DefaultDiscriminator.clear();
    } // end track counting

    // Secondary Vertex Tagger, for testing the DATA driven approach
    if (fDoSVAnalysis)
    {
      double vtxPos[3] = {fPrimaryVertex->GetX(), fPrimaryVertex->GetY(), fPrimaryVertex->GetZ()};
      double covMatrix[6] = {0};
      fPrimaryVertex->GetCovarianceMatrix(covMatrix);
      AliESDVertex *esdVtx = new AliESDVertex(vtxPos, covMatrix, fPrimaryVertex->GetChi2(), fPrimaryVertex->GetNContributors());

      // 3 Prong Vertex
      TClonesArray *secVertexArr3Prong = 0;
      std::vector<pair<double, int>> arrDispersion3Prong;
      arrDispersion3Prong.reserve(5);

      secVertexArr3Prong = new TClonesArray("AliAODVertex");
      int nDauRejCount3Prong = 0;
      int nVtx3Prong = fVtxTagger3Prong->FindVertices(jetrec,
                                                      static_cast<AliParticleContainer *>(fParticleCollArray.At(0))->GetArray(),
                                                      fAODIn,
                                                      esdVtx,
                                                      fAODIn->GetMagneticField(),
                                                      secVertexArr3Prong,
                                                      0,
                                                      arrDispersion3Prong,
                                                      nDauRejCount3Prong);

      if (nVtx3Prong > 0)
      {
        double *decLenXY = new double[nVtx3Prong];
        double *errdecLenXY = new double[nVtx3Prong];
        double *sigdecLenXY = new double[nVtx3Prong];
        double *invMasses = new double[nVtx3Prong];
        double *sigmavertex = new double[nVtx3Prong];

        int *idxLxy = new int[nVtx3Prong];

        for (int iv = 0; iv < secVertexArr3Prong->GetEntriesFast(); iv++)
        {
          AliAODVertex *secVtx = (AliAODVertex *)(secVertexArr3Prong->UncheckedAt(iv));

          // Calculate vtx distance
          double effX = secVtx->GetX() - esdVtx->GetX();
          double effY = secVtx->GetY() - esdVtx->GetY();
          // double effZ = secVtx->GetZ() - esdVtx->GetZ();

          // ##### Vertex properties
          // vertex dispersion
          sigmavertex[iv] = arrDispersion3Prong[iv].first;

          // invariant mass
          invMasses[iv] = fVtxTagger3Prong->GetVertexInvariantMass(secVtx);

          // signed length
          decLenXY[iv] = fPrimaryVertex->DistanceXYToVertex(secVtx);
          double jetP[3];
          jetrec->PxPyPz(jetP);
          double signLxy = effX * jetP[0] + effY * jetP[1];
          if (signLxy < 0.)
            decLenXY[iv] *= -1.;

          errdecLenXY[iv] = fPrimaryVertex->ErrorDistanceXYToVertex(secVtx);

          sigdecLenXY[iv] = decLenXY[iv] / errdecLenXY[iv];
        }

        TMath::Sort(nVtx3Prong, decLenXY, idxLxy);

        double invariantMass = invMasses[idxLxy[0]];
        double dispersion = sigmavertex[idxLxy[0]];
        double LxySign = sigdecLenXY[idxLxy[0]];

        if (LxySign > fThresholdLxy && dispersion < fMaxDespersion)
        {
          taggedSV3Prong = true;
          if (!firstJetFound && fEmbeddPerpendicular && fFillControlHists)
          {
            firstJetFound = true;
            double EmbeddingPt = GetDeltaPtPerpEmbedding(jetrec->Eta(), jetrec->Phi());
            fh2DeltaPtEmbeddCorrelationPerpendicular->Fill(EmbeddingPt, fJetContainerData->GetRhoVal());
          }
        }

        if (fFillControlHists)
        {
          fHistDispersion3Prong->Fill(jetPt, dispersion);
          if (fIsPythia)
          {
            if (jetFlavor == 0)
              fHistDispersion3ProngUnidentified->Fill(jetPt, dispersion);
            if (jetFlavor == 1)
              fHistDispersion3Pronglf->Fill(jetPt, dispersion);
            if (jetFlavor == 2)
              fHistDispersion3Prongc->Fill(jetPt, dispersion);
            if (jetFlavor == 3)
              fHistDispersion3Prongb->Fill(jetPt, dispersion);
          }
          if (dispersion < fMaxDespersion)
          {
            fHistSV3Prong->Fill(jetPt, LxySign, invariantMass);

            if (fIsPythia)
            {
              if (jetFlavor == 0)
                fHistSV3ProngUnidentified->Fill(jetPt, LxySign, invariantMass);
              if (jetFlavor == 1)
                fHistSV3Pronglf->Fill(jetPt, LxySign, invariantMass);
              if (jetFlavor == 2)
                fHistSV3Prongc->Fill(jetPt, LxySign, invariantMass);
              if (jetFlavor == 3)
                fHistSV3Prongb->Fill(jetPt, LxySign, invariantMass);
            }
          }
        }

        secVertexArr3Prong->Clear();
        delete secVertexArr3Prong;
        delete[] idxLxy;
        delete[] decLenXY;
        delete[] invMasses;
        delete[] errdecLenXY;
        delete[] sigdecLenXY;
        delete[] sigmavertex;
      }
      else
      {
        secVertexArr3Prong->Clear();
        delete secVertexArr3Prong;
      }

      // 2 Prong Vertex
      TClonesArray *secVertexArr2Prong = 0;
      std::vector<pair<double, int>> arrDispersion2Prong;
      arrDispersion2Prong.reserve(5);

      secVertexArr2Prong = new TClonesArray("AliAODVertex");
      int nDauRejCount2Prong = 0;
      int nVtx2Prong = fVtxTagger2Prong->FindVertices(jetrec,
                                                      static_cast<AliParticleContainer *>(fParticleCollArray.At(0))->GetArray(),
                                                      fAODIn,
                                                      esdVtx,
                                                      fAODIn->GetMagneticField(),
                                                      secVertexArr2Prong,
                                                      0,
                                                      arrDispersion2Prong,
                                                      nDauRejCount2Prong);

      if (nVtx2Prong > 0)
      {
        double *decLenXY = new double[nVtx2Prong];
        double *errdecLenXY = new double[nVtx2Prong];
        double *sigdecLenXY = new double[nVtx2Prong];
        double *invMasses = new double[nVtx2Prong];
        double *sigmavertex = new double[nVtx2Prong];

        int *idxLxy = new int[nVtx2Prong];

        for (int iv = 0; iv < secVertexArr2Prong->GetEntriesFast(); iv++)
        {
          AliAODVertex *secVtx = (AliAODVertex *)(secVertexArr2Prong->UncheckedAt(iv));

          // Calculate vtx distance
          double effX = secVtx->GetX() - esdVtx->GetX();
          double effY = secVtx->GetY() - esdVtx->GetY();
          // double effZ = secVtx->GetZ() - esdVtx->GetZ();

          // ##### Vertex properties
          // vertex dispersion
          sigmavertex[iv] = arrDispersion2Prong[iv].first;

          // invariant mass
          invMasses[iv] = fVtxTagger2Prong->GetVertexInvariantMass(secVtx);

          // signed length
          decLenXY[iv] = fPrimaryVertex->DistanceXYToVertex(secVtx);
          double jetP[3];
          jetrec->PxPyPz(jetP);
          double signLxy = effX * jetP[0] + effY * jetP[1];
          if (signLxy < 0.)
            decLenXY[iv] *= -1.;

          errdecLenXY[iv] = fPrimaryVertex->ErrorDistanceXYToVertex(secVtx);

          sigdecLenXY[iv] = decLenXY[iv] / errdecLenXY[iv];
        }

        TMath::Sort(nVtx2Prong, decLenXY, idxLxy);

        double invariantMass = invMasses[idxLxy[0]];
        double dispersion = sigmavertex[idxLxy[0]];
        double LxySign = sigdecLenXY[idxLxy[0]];

        if (LxySign > fThresholdLxy && dispersion < fMaxDespersion)
        {
          taggedSV2Prong = true;
        }

        if (fFillControlHists)
        {
          fHistDispersion2Prong->Fill(jetPt, dispersion);
          if (fIsPythia)
          {
            if (jetFlavor == 0)
              fHistDispersion2ProngUnidentified->Fill(jetPt, dispersion);
            if (jetFlavor == 1)
              fHistDispersion2Pronglf->Fill(jetPt, dispersion);
            if (jetFlavor == 2)
              fHistDispersion2Prongc->Fill(jetPt, dispersion);
            if (jetFlavor == 3)
              fHistDispersion2Prongb->Fill(jetPt, dispersion);
          }
          if (dispersion < fMaxDespersion)
          {
            fHistSV2Prong->Fill(jetPt, LxySign, invariantMass);

            if (fIsPythia)
            {
              if (jetFlavor == 0)
                fHistSV2ProngUnidentified->Fill(jetPt, LxySign, invariantMass);
              if (jetFlavor == 1)
                fHistSV2Pronglf->Fill(jetPt, LxySign, invariantMass);
              if (jetFlavor == 2)
                fHistSV2Prongc->Fill(jetPt, LxySign, invariantMass);
              if (jetFlavor == 3)
                fHistSV2Prongb->Fill(jetPt, LxySign, invariantMass);
            }
          }
        }

        secVertexArr2Prong->Clear();
        delete secVertexArr2Prong;
        delete[] idxLxy;
        delete[] decLenXY;
        delete[] invMasses;
        delete[] errdecLenXY;
        delete[] sigdecLenXY;
        delete[] sigmavertex;
      }
      else
      {
        secVertexArr2Prong->Clear();
        delete secVertexArr2Prong;
      }

      delete esdVtx;

    } // end SV

    if (TaggedSecond)
    {
      FillCandidateHFJet(jetrec, fInputList);
      if (fIsPythia)
      {
        FillCandidateHFJetMC(jetrec, jetFlavor);
      }
    }

    if (taggedSV3Prong)
    {
      FillCandidateHFJet(jetrec, fInputList);
      if (fIsPythia)
      {
        FillCandidateHFJetMC(jetrec, jetFlavor);
      }
    }

  } // End jet loop

  if (fFillControlHists)
  {
    if (TaggedFirst)
    {
      f2histRhoVsDeltaPtFirst->Fill(randomConePt, fJetContainerData->GetRhoVal());
    }
    if (TaggedSecond)
    {
      f2histRhoVsDeltaPtSecond->Fill(randomConePt, fJetContainerData->GetRhoVal());
    }
    if (TaggedThird)
    {
      f2histRhoVsDeltaPtThird->Fill(randomConePt, fJetContainerData->GetRhoVal());
    }
  }

  return true;
}

// ######################################################################################## Jet matching 1/4
bool AliJHFTagging::MatchJetsGeometricDefault()
{
  double matchingpar1 = 0.25;
  double matchingpar2 = 0.25;
  if (!fJetContainerData || !fJetContainerData->GetArray() || !fJetContainerMC || !fJetContainerMC->GetArray())
  {
    return false;
  }
  DoJetLoop();
  AliEmcalJet *jet1 = 0;
  AliEmcalJet *jet2 = 0;

  fJetContainerData->ResetCurrentID();
  while ((jet1 = fJetContainerData->GetNextJet()))
  {
    jet2 = jet1->ClosestJet();
    if (!jet2)
      continue;
    if (jet2->ClosestJet() != jet1)
      continue;
    if (jet1->ClosestJetDistance() > matchingpar1 || jet2->ClosestJetDistance() > matchingpar2)
      continue;
    // Matched jet found
    jet1->SetMatchedToClosest(1);
    jet2->SetMatchedToClosest(1);
  }

  return true;
}
// ######################################################################################## Jet matching 2/4
void AliJHFTagging::DoJetLoop()
{
  // Do the jet loop.
  double minjetpt = 1.;

  AliEmcalJet *jet1 = 0;
  AliEmcalJet *jet2 = 0;
  fJetContainerMC->ResetCurrentID();
  while ((jet2 = fJetContainerMC->GetNextJet()))
    jet2->ResetMatching();
  fJetContainerData->ResetCurrentID();
  while ((jet1 = fJetContainerData->GetNextJet()))
  {
    jet1->ResetMatching();
    if (jet1->MCPt() < minjetpt)
      continue;
    fJetContainerMC->ResetCurrentID();
    while ((jet2 = fJetContainerMC->GetNextJet()))
    {
      SetMatchingLevel(jet1, jet2, 1);
    } // jet2 loop
  }   // jet1 loop
}
// ######################################################################################## Jet matching 3/4
void AliJHFTagging::SetMatchingLevel(AliEmcalJet *jet1, AliEmcalJet *jet2, int matching)
{
  double d1 = -1;
  double d2 = -1;

  switch (matching)
  {
  case 1:
    GetGeometricalMatchingLevel(jet1, jet2, d1);
    d2 = d1;
    break;
  default:
    break;
  }
  if (d1 >= 0)
  {

    if (d1 < jet1->ClosestJetDistance())
    {
      jet1->SetSecondClosestJet(jet1->ClosestJet(), jet1->ClosestJetDistance());
      jet1->SetClosestJet(jet2, d1);
    }
    else if (d1 < jet1->SecondClosestJetDistance())
    {
      jet1->SetSecondClosestJet(jet2, d1);
    }
  }
  if (d2 >= 0)
  {
    if (d2 < jet2->ClosestJetDistance())
    {
      jet2->SetSecondClosestJet(jet2->ClosestJet(), jet2->ClosestJetDistance());
      jet2->SetClosestJet(jet1, d2);
    }
    else if (d2 < jet2->SecondClosestJetDistance())
    {
      jet2->SetSecondClosestJet(jet1, d2);
    }
  }
}

// ######################################################################################## Jet matching 4/4
void AliJHFTagging::GetGeometricalMatchingLevel(AliEmcalJet *jet1, AliEmcalJet *jet2, double &d) const
{
  double deta = jet2->Eta() - jet1->Eta();
  double dphi = jet2->Phi() - jet1->Phi();
  dphi = TVector2::Phi_mpi_pi(dphi);
  d = TMath::Sqrt(deta * deta + dphi * dphi);
}

////////////////////////////////////////////////////////////////////////////////
double AliJHFTagging::GetDeltaPtRandomCone()
{

  double deltaPt = -1000.;
  AliParticleContainer *partcont = 0x0;
  partcont = static_cast<AliParticleContainer *>(fParticleCollArray.At(0));

  double jetradius = fJetContainerData->GetJetRadius();
  double minEta = -0.5;
  double maxEta = 0.5;
  double tmpRandConeEta = -999;
  double tmpRandConePhi = -999;
  double tmpConePt = -1.;

  AliEmcalJet *LeadingJet = NULL;
  AliEmcalJet *SubLeadingJet = NULL;
  double LJeta = 999;
  double LJphi = 999;

  double SLJeta = 999;
  double SLJphi = 999;

  if (fJetContainerData)
  {

    Float_t maxJetPts[] = {0, 0};

    fJetContainerData->ResetCurrentID();

    AliEmcalJet *jet = 0x0;

    while ((jet = fJetContainerData->GetNextAcceptJet()))
    {

      if (!jet)
        continue;

      if (jet->Pt() > maxJetPts[0])
      {
        maxJetPts[1] = maxJetPts[0];
        SubLeadingJet = LeadingJet;
        maxJetPts[0] = jet->Pt();
        LeadingJet = jet;
      }
      else if (jet->Pt() > maxJetPts[1])
      {
        maxJetPts[1] = jet->Pt();
        SubLeadingJet = jet;
      }
    }

    if (LeadingJet)
    {
      LJeta = LeadingJet->Eta();
      LJphi = LeadingJet->Phi();
    }
    if (SubLeadingJet)
    {
      SLJeta = SubLeadingJet->Eta();
      SLJphi = SubLeadingJet->Phi();
    }
  }

  double dLJ = 0;
  double dSLJ = 0;
  int repeats = 0;

  do
  {
    tmpRandConeEta = minEta + fRandom->Rndm() * (maxEta - minEta);
    tmpRandConePhi = fRandom->Rndm() * TMath::TwoPi();
    dLJ = TMath::Sqrt((LJeta - tmpRandConeEta) * (LJeta - tmpRandConeEta) + (LJphi - tmpRandConePhi) * (LJphi - tmpRandConePhi));
    dSLJ = TMath::Sqrt((SLJeta - tmpRandConeEta) * (SLJeta - tmpRandConeEta) + (SLJphi - tmpRandConePhi) * (SLJphi - tmpRandConePhi));
    repeats++;

  } while (dLJ < 0.45 || dSLJ < 0.45);

  partcont->ResetCurrentID();
  while (AliAODTrack *trackAOD = static_cast<AliAODTrack *>(partcont->GetNextAcceptParticle()))
  {
    if (!((trackAOD)->TestFilterBit(1 << 4)) && !((trackAOD)->TestFilterBit(1 << 9)))
      continue;

    if (fabs(trackAOD->Eta()) > 0.9)
      continue;

    if (trackAOD->Pt() < 0.15)
      continue;

    if (sqrt((trackAOD->Eta() - tmpRandConeEta) * (trackAOD->Eta() - tmpRandConeEta) +
             TVector2::Phi_mpi_pi((trackAOD->Phi() - tmpRandConePhi)) *
                 TVector2::Phi_mpi_pi((trackAOD->Phi() - tmpRandConePhi))) < jetradius)
    {
      tmpConePt += trackAOD->Pt();
    }
  }

  if (tmpConePt > 0)
  {
    deltaPt = tmpConePt - jetradius * jetradius * TMath::Pi() * fJetContainerData->GetRhoVal();
    return deltaPt;
  }
  return deltaPt;
}

//________________________________________________________________________
double AliJHFTagging::GetDeltaPtPerpEmbedding(double signalEta, double signalPhi)
{

  fFastJetWrapper->Clear();
  //----------------Generating NEW perpendicular track

  AliParticleContainer *partcont = 0x0;
  partcont = static_cast<AliParticleContainer *>(fParticleCollArray.At(0));

  double gen_pt = fTrackGenerator->Uniform(0, 100); // this will be pT of the track that you embedd
  TLorentzVector lVec;
  lVec.SetPtEtaPhiM(gen_pt, signalEta, signalPhi + TMath::Pi() / 2, 0);               // here ignalEta,signalPhi  are some directions that you choose
  fFastJetWrapper->AddInputVector(lVec.Px(), lVec.Py(), lVec.Pz(), lVec.E(), -99999); // fill embedded track to the array of proto-jets

  //-----Filling   fFastJetWrapper with tracks from track container

  partcont->ResetCurrentID();
  while (AliAODTrack *trk = static_cast<AliAODTrack *>(partcont->GetNextAcceptParticle()))
  {
    if (!((trk)->TestFilterBit(1 << 4)) && !((trk)->TestFilterBit(1 << 9)))
      continue;

    if (fabs(trk->Eta()) > 0.9)
      continue;

    if (trk->Pt() < 0.15)
      continue;

    fFastJetWrapper->AddInputVector(trk->Px(), trk->Py(), trk->Pz(), trk->P(), 1); // fill reconstructed tracks to the array of proto-jets
  }

  fFastJetWrapper->Run(); // this creates jets including the embedded track

  //--------------- DelPt analysis of the new container
  std::vector<fastjet::PseudoJet> jets_incl = fFastJetWrapper->GetInclusiveJets(); // this is the new jet array which included also the jet with embedded track
  double deltaPtEmb = -1000.;
  double sumTrkEmbeddedPt = 0;

  for (UInt_t ijet = 0; ijet < jets_incl.size(); ++ijet)
  {                                                                                          // loop over these jets and search for jet with embedded track.
    std::vector<fastjet::PseudoJet> constituents(fFastJetWrapper->GetJetConstituents(ijet)); // get list of constituents of each jet
    sumTrkEmbeddedPt = 0;
    for (UInt_t ic = 0; ic < constituents.size(); ++ic)
    { // loop over constituents
      if (constituents[ic].user_index() == -99999)
      {                                            // this is the constituent which was embedded
        sumTrkEmbeddedPt += constituents[ic].pt(); // get pT of this embedded track
        break;                                     // we have embedded just 1 track so we can break, there will be not other in this event
      }
    }

    if (sumTrkEmbeddedPt > 0)
    {                                                                                                                       // the jet with embedded track was found
      deltaPtEmb = jets_incl.at(ijet).pt() - jets_incl.at(ijet).area() * fJetContainerData->GetRhoVal() - sumTrkEmbeddedPt; // calculate delta pT and subtract pT of embedded track
      break;
    }
  }

  return deltaPtEmb;
}

// ########################################################################################Track Selection
bool AliJHFTagging::IsTrackAcceptedBJetCuts(AliAODTrack *track, int jetFlavour)
{

  int iCutIndex = 0; // indicator of current selection step

  FillCandidateJet(iCutIndex, jetFlavour);
  iCutIndex++;

  if (!(((AliAODTrack *)track)->TestFilterBit(9) || ((AliAODTrack *)track)->TestFilterBit(4)))
    return false;
  FillCandidateJet(iCutIndex, jetFlavour);
  iCutIndex++;

  if (track->Pt() < fTCMinTrackPt)
    return false; // 0.5
  FillCandidateJet(iCutIndex, jetFlavour);
  iCutIndex++;

  if (TMath::Abs(track->Eta()) > 0.9)
    return false;
  FillCandidateJet(iCutIndex, jetFlavour);
  iCutIndex++;

  ULong_t status = track->GetStatus();
  if (!(status & AliAODTrack::kTPCrefit))
    return false;
  FillCandidateJet(iCutIndex, jetFlavour);
  iCutIndex++;

  if (!(status & AliAODTrack::kITSrefit))
    return false;
  FillCandidateJet(iCutIndex, jetFlavour);
  iCutIndex++;

  if (!track->HasPointOnITSLayer(0) && !track->HasPointOnITSLayer(1))
    return false;
  FillCandidateJet(iCutIndex, jetFlavour);
  iCutIndex++;

  if (track->GetNcls(1) < fTCMinClusTPC)
    return false; // 80
  FillCandidateJet(iCutIndex, jetFlavour);
  iCutIndex++;

  int SPDSSDHits = track->HasPointOnITSLayer(0) + track->HasPointOnITSLayer(1) + track->HasPointOnITSLayer(4) + track->HasPointOnITSLayer(5);
  if (SPDSSDHits < fTCMinHitsITS)
    return false; // 2
  FillCandidateJet(iCutIndex, jetFlavour);
  iCutIndex++;

  if (track->Chi2perNDF() >= fTCMaxChi2pNDF)
    return false; // 5
  FillCandidateJet(iCutIndex, jetFlavour);
  iCutIndex++;

  return true;
}
//____________________________________________________
void AliJHFTagging::FillCandidateJet(int CutIndex, int JetFlavor)
{

  if (!fFillControlHists)
  {
    return;
  }

  fhistInclusiveJetCuts->Fill(CutIndex);
  if (fIsPythia)
  {
    if (JetFlavor == 3)
      fhistbJetCuts->Fill(CutIndex);
    if (JetFlavor == 2)
      fhistcJetCuts->Fill(CutIndex);
    if (JetFlavor == 1)
      fhistlfJetCuts->Fill(CutIndex);
  }
}

// ######################################################################################## Post-process ImpPar
double AliJHFTagging::GetValImpactParameter(TTypeImpPar type, double *impar, double *cov)
{
  double result = -999;
  double dFdx = 0;
  double dFdy = 0;

  switch (type)
  {
  case kXY:
    result = impar[0];
    break;
  case kXYSig:
    result = impar[0] / TMath::Sqrt(cov[0]);
    break;
  case kXYZ:
    result = TMath::Sqrt(impar[0] * impar[0] + impar[1] * impar[1]);
    break;
  case kXYZSig:
    result = TMath::Sqrt(impar[0] * impar[0] + impar[1] * impar[1]);
    dFdx = 2 * impar[0] / result;
    dFdy = 2 * impar[1] / result;
    result /= TMath::Sqrt(cov[0] * dFdx * dFdx + cov[2] * dFdy * dFdy + 2 * cov[1] * dFdx * dFdy);
    break;
  case kXYZSigmaOnly:
    dFdx = 2 * impar[0] / result;
    dFdy = 2 * impar[1] / result;
    result = TMath::Sqrt(cov[0] * dFdx * dFdx + cov[2] * dFdy * dFdy + 2 * cov[1] * dFdx * dFdy);
    break;

  default:
    break;
  }
  return result;
}

// ######################################################################################## Calculate signed  impact parameters
bool AliJHFTagging::CalculateJetSignedTrackImpactParameter(AliAODTrack *track, AliEmcalJet *jet, double *impar, double *cov, double &sign, double &dcajetrack, double &lineardecaylength)
{

  double vtxPos[3];
  AliESDVertex *vtxESDNew = 0x0;

  AliExternalTrackParam etp;
  etp.CopyFromVTrack(track);

  int skipped[1] = {-1};
  int id = (int)track->GetID();
  if (id < 0)
    return false;
  skipped[0] = id;
  fVertexer->SetSkipTracks(1, skipped);
  vtxESDNew = fVertexer->FindPrimaryVertex(fAODIn);
  if (!vtxESDNew)
    return false;
  if (vtxESDNew->GetNContributors() <= 0)
  {
    vtxESDNew = NULL;
    delete vtxESDNew;
    return false;
  }
  // convert to AliAODVertex
  double cova[6];
  vtxESDNew->GetXYZ(vtxPos);     // position
  vtxESDNew->GetCovMatrix(cova); // covariance matrix

  // Calculate Impact Parameters
  if (!etp.PropagateToDCA(vtxESDNew, fAODIn->GetMagneticField(), 3., impar, cov))
    return false;

  delete vtxESDNew;

  // Calculate Sign
  double posdcatrack[3] = {0., 0., 0.};
  etp.GetXYZ(posdcatrack);
  double ipvector3[3] = {posdcatrack[0] - vtxPos[0], posdcatrack[1] - vtxPos[1], posdcatrack[2] - vtxPos[2]};
  sign = TMath::Sign(1., ipvector3[0] * jet->Px() + ipvector3[1] * jet->Py() + ipvector3[2] * jet->Pz());

  // Calculate decay legnth and track jet DCA against new vertex
  double bpxpypz[3] = {jet->Px(), jet->Py(), jet->Pz()};
  double bcv[21] = {0};
  AliExternalTrackParam bjetparam(vtxPos, bpxpypz, bcv, (Short_t)0);
  double xa = 0., xb = 0.;
  bjetparam.GetDCA(&etp, fAODIn->GetMagneticField(), xa, xb);
  double xyz[3] = {0., 0., 0.};
  double xyzb[3] = {0., 0., 0.};
  bjetparam.GetXYZAt(xa, fAODIn->GetMagneticField(), xyz);
  etp.GetXYZAt(xb, fAODIn->GetMagneticField(), xyzb);
  double bdecaylength =
      TMath::Sqrt(
          (vtxPos[0] - xyzb[0]) * (vtxPos[0] - xyzb[0]) +
          (vtxPos[1] - xyzb[1]) * (vtxPos[1] - xyzb[1]) +
          (vtxPos[2] - xyzb[2]) * (vtxPos[2] - xyzb[2]));
  dcajetrack =
      TMath::Sqrt(
          (xyzb[0] - xyz[0]) * (xyzb[0] - xyz[0]) +
          (xyzb[1] - xyz[1]) * (xyzb[1] - xyz[1]) +
          (xyzb[2] - xyz[2]) * (xyzb[2] - xyz[2]));
  if (bdecaylength > 0)
    lineardecaylength = bdecaylength;

  return true;
}

void AliJHFTagging::FillCandidateHFJet(AliEmcalJet *jet, TClonesArray *inputList)
{
  new ((*inputList)[inputList->GetEntriesFast()])(AliEmcalJet *)(jet);
}

void AliJHFTagging::FillCandidateHFJetMC(AliEmcalJet *jet, short jetFlavor)
{
  if (jetFlavor == 3)
    new ((*fInputListb)[fInputListc->GetEntriesFast()])(AliEmcalJet *)(jet);
  if (jetFlavor == 2)
    new ((*fInputListc)[fInputListc->GetEntriesFast()])(AliEmcalJet *)(jet);
  if (jetFlavor == 0 || jetFlavor == 1)
    new ((*fInputListlf)[fInputListlf->GetEntriesFast()])(AliEmcalJet *)(jet);
}