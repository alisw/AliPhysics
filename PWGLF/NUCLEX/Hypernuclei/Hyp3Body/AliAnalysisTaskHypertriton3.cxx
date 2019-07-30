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
//
// Author:
// P. Fecchio, pfecchio@cern.ch
///////////////////////////////////////////////////////////////////////////

#include <climits>
#include <numeric>
#include <utility>
#include <vector>

/// ROOT includes
#include <Riostream.h>
#include <TChain.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TList.h>
#include <TObjArray.h>
#include <TVector3.h>

/// AliRoot icludes
#include "AliAnalysisTaskSE.h"
#include "AliCentrality.h"
#include "AliESD.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDVertex.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliExternalTrackParam.h"
#include "AliInputEventHandler.h"
#include "AliMCEvent.h"
#include "AliMultEstimator.h"
#include "AliMultInput.h"
#include "AliMultSelection.h"
#include "AliMultVariable.h"
#include "AliMultiplicity.h"
#include "AliOADBContainer.h"
#include "AliOADBMultSelection.h"
#include "AliPDG.h"
#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliPhysicsSelection.h"
#include "TLorentzVector.h"

#include "AliAnalysisTaskHypertriton3.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskHypertriton3);

namespace {
struct TrackMC {
  AliESDtrack *track;
  AliVParticle *mother;
  AliVParticle *particle;
  int motherId;
  int partId;
};

struct CandidateMC {
  AliESDtrack *track_deu;
  AliESDtrack *track_p;
  AliESDtrack *track_pi;
  AliVParticle *part1;
  AliVParticle *part2;
  AliVParticle *part3;
  AliVParticle *mother;
  int motherId;
  unsigned char status;
};

bool IsHyperTriton3Daughter(AliMCEvent *mcEvent, const AliVParticle *vPart) {

  int nDaughters = 0;

  int lLabelMother = vPart->GetMother();
  if (lLabelMother < 0 || !mcEvent->IsPhysicalPrimary(lLabelMother)) return false;

  AliVParticle *vMotherPart = mcEvent->GetTrack(lLabelMother);
  int lMotherPDG            = vMotherPart->PdgCode();
  if (std::abs(lMotherPDG) != 1010010030) return false;

  for (int iD = vMotherPart->GetDaughterFirst(); iD <= vMotherPart->GetDaughterLast(); iD++) {
    AliVParticle *dPart = mcEvent->GetTrack(iD);
    int dPartPDG        = dPart->PdgCode();
    if (std::abs(dPartPDG) != 11) nDaughters++;
  }
  if (nDaughters != 3) return false;

  return true;
}

bool IsFakeCandidate(AliMCEvent *mcEvent, int mId, AliVParticle *p1, AliVParticle *p2, AliVParticle *p3) {

  AliVParticle *vMother = mcEvent->GetTrack(mId);
  bool fake             = false;

  for (int iD = vMother->GetDaughterFirst(); iD <= vMother->GetDaughterLast(); iD++) {

    AliVParticle *dPart = mcEvent->GetTrack(iD);
    int dPartPDG        = dPart->PdgCode();
    if (vMother->PdgCode() == 1010010030) {
      if (dPartPDG == 1000010020) {
        if (dPart->GetLabel() != p1->GetLabel()) {
          fake = true;
          break;
        }
      }
      if (dPartPDG == 2212) {
        if (dPart->GetLabel() != p2->GetLabel()) {
          fake = true;
          break;
        }
      }
      if (dPartPDG == -211) {
        if (dPart->GetLabel() != p3->GetLabel()) {
          fake = true;
          break;
        }
      }
    }

    if (vMother->PdgCode() == -1010010030) {
      if (dPartPDG == -1000010020) {
        if (dPart->GetLabel() != p1->GetLabel()) {
          fake = true;
          break;
        }
      }
      if (dPartPDG == -2212) {
        if (dPart->GetLabel() != p2->GetLabel()) {
          fake = true;
          break;
        }
      }
      if (dPartPDG == 211) {
        if (dPart->GetLabel() != p3->GetLabel()) {
          fake = true;
          break;
        }
      }
    }
  }
  return fake;
}

bool IsReflectionCandidate(AliMCEvent *mcEvent, int mId, AliVParticle *p1, AliVParticle *p2, AliVParticle *p3) {

  int n_deu = 0;
  int n_p   = 0;
  int n_pi  = 0;

  AliVParticle *vMother = mcEvent->GetTrack(mId);
  int mPartPDG          = vMother->PdgCode();

  int vPDG[3] = {p1->PdgCode(), p2->PdgCode(), p3->PdgCode()};

  for (int iPDG = 0; iPDG < 3; iPDG++) {

    if (mPartPDG == 1010010030) {
      if (vPDG[iPDG] == 1000010020) n_deu++;
      if (vPDG[iPDG] == 2212) n_p++;
      if (vPDG[iPDG] == -211) n_pi++;
    }

    if (mPartPDG == -1010010030) {
      if (vPDG[iPDG] == -1000010020) n_deu++;
      if (vPDG[iPDG] == -2212) n_p++;
      if (vPDG[iPDG] == 211) n_pi++;
    }
  }
  return (n_deu == 1 && n_p == 1 && n_pi == 1);
}

bool HasTOF(AliVTrack *track) {
  const bool hasTOFout  = track->GetStatus() & AliVTrack::kTOFout;
  const bool hasTOFtime = track->GetStatus() & AliVTrack::kTIME;
  const float len       = track->GetIntegratedLength();
  return hasTOFout && hasTOFtime && (len > 350.);
}

} // namespace


//________________________________________________________________________
AliAnalysisTaskHypertriton3::AliAnalysisTaskHypertriton3(TString taskname)
    : AliAnalysisTaskSE(taskname.Data()),
      // support objects
      fEventCuts{},            //
      fPIDResponse{nullptr},   //
      fESDtrackCuts{nullptr},  //
      fPrimaryVertex{nullptr}, //
      // setting parameters
      fCosPoiningAngleLimit{0}, //
      // output objects
      fOutputList{nullptr},                            //
      fTree{nullptr},                                  //
      fTreeHyp3BodyVarTracks{nullptr},                 //
      fTreeHyp3BodyVarNclsTPC{0},                      //
      fTreeHyp3BodyVarNclsITS{0},                      //
      fTreeHyp3BodyVarGlobalChi2{0},                   //
      fTreeHyp3BodyVarNsigmaTPC{0},                    //
      fTreeHyp3BodyVarNsigmaTOF{0},                    //
      fTreeHyp3BodyVarFlags{0},                        //
      fTreeHyp3BodyVarPDGcodes{0},                     //
      fTreeHyp3BodyVarEventId{0},                      //
      fTreeHyp3BodyVarMotherId{0},                     //
      fTreeHyp3BodyVarCandStat{0},                     //
      fTreeHyp3BodyVarTruePx{0},                       //
      fTreeHyp3BodyVarTruePy{0},                       //
      fTreeHyp3BodyVarTruePz{0},                       //
      fTreeHyp3BodyVarDecayVx{0},                      //
      fTreeHyp3BodyVarDecayVy{0},                      //
      fTreeHyp3BodyVarDecayVz{0},                      //
      fTreeHyp3BodyVarDecayT{0},                       //
      fTreeHyp3BodyVarPVx{0},                          //
      fTreeHyp3BodyVarPVy{0},                          //
      fTreeHyp3BodyVarPVz{0},                          //
      fTreeHyp3BodyVarPVt{0},                          //
      fTreeHyp3BodyVarMagneticField{0},                //
      fHistEventCounter{nullptr},                      //
      fHistCentrality{nullptr},                        //
      fHistGeneratedPtVsYVsCentralityHypTrit{nullptr}, //
      fHistGeneratedPtVsYVsCentralityAntiHypTrit{nullptr} {

  // Standard Output
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
}

//________________________________________________________________________
AliAnalysisTaskHypertriton3::~AliAnalysisTaskHypertriton3() {
  // destructor
  if (fOutputList) {
    delete fOutputList;
    fOutputList = nullptr;
  }
  if (fTree) {
    delete fTree;
    fTree = nullptr;
  }
}

//________________________________________________________________________
void AliAnalysisTaskHypertriton3::UserCreateOutputObjects() {

  AliAnalysisManager *fMgr = AliAnalysisManager::GetAnalysisManager();
  if (!fMgr) AliFatal("Could not find analysis manager.");
  AliInputEventHandler *fHandl = (AliInputEventHandler *)fMgr->GetInputEventHandler();
  if (!fHandl) AliFatal("No input event handler.");
  fPIDResponse = fHandl->GetPIDResponse();
  fHandl->SetNeedField();

  // Multiplicity
  if (!fESDtrackCuts) {
    fESDtrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE, kFALSE);
    fESDtrackCuts->SetPtRange(0.15); // adding pt cut
    fESDtrackCuts->SetEtaRange(-1.0, 1.0);
  }

  // Create a TList with Histograms
  fOutputList = new TList();
  fOutputList->SetOwner();
  // fEventCuts.AddQAplotsToList(fOutputList);

  // Histogram Output: Event-by-Event
  fHistCentrality =
      new TH1D("fHistCentrality", "WARNING: no pileup rejection applied!;Centrality;Event Count", 100, 0, 100);
  fOutputList->Add(fHistCentrality);

  // Histogram Output: Efficiency Denominator
  fHistGeneratedPtVsYVsCentralityHypTrit =
      new TH3D("fHistGeneratedPtVsYVsCentralityHypTrit", ";pT;y;centrality", 500, 0, 25, 40, -1.0, 1.0, 100, 0, 100);
  fOutputList->Add(fHistGeneratedPtVsYVsCentralityHypTrit);
  fHistGeneratedPtVsYVsCentralityAntiHypTrit = new TH3D("fHistGeneratedPtVsYVsCentralityAntiHypTrit",
                                                        ";pT;y;centrality", 500, 0, 25, 40, -1.0, 1.0, 100, 0, 100);
  fOutputList->Add(fHistGeneratedPtVsYVsCentralityAntiHypTrit);

  // Histogram Output: Event-by-Event
  fHistEventCounter = new TH1D("fHistEventCounter", ";Evt. Sel. Step;Count", 2, 0, 2);
  fHistEventCounter->GetXaxis()->SetBinLabel(1, "Processed");
  fHistEventCounter->GetXaxis()->SetBinLabel(2, "Selected");
  fOutputList->Add(fHistEventCounter);

  fTree = new TTree("fTreeHyperTriton3Body", "HyperTriton3BodyCandidates");

  fTree->Branch("fTreeHyp3BodyVarTrack0", &fTreeHyp3BodyVarTracks[0], 16000, 99);
  fTree->Branch("fTreeHyp3BodyVarTrack1", &fTreeHyp3BodyVarTracks[1], 16000, 99);
  fTree->Branch("fTreeHyp3BodyVarTrack2", &fTreeHyp3BodyVarTracks[2], 16000, 99);

  fTree->Branch("fTreeHyp3BodyVarNclsTPC0", &fTreeHyp3BodyVarNclsTPC[0], "fTreeHyp3BodyVarNclsTPC0/I");
  fTree->Branch("fTreeHyp3BodyVarNclsTPC1", &fTreeHyp3BodyVarNclsTPC[1], "fTreeHyp3BodyVarNclsTPC1/I");
  fTree->Branch("fTreeHyp3BodyVarNclsTPC2", &fTreeHyp3BodyVarNclsTPC[2], "fTreeHyp3BodyVarNclsTPC2/I");

  fTree->Branch("fTreeHyp3BodyVarNclsITS0", &fTreeHyp3BodyVarNclsITS[0], "fTreeHyp3BodyVarNclsITS0/I");
  fTree->Branch("fTreeHyp3BodyVarNclsITS1", &fTreeHyp3BodyVarNclsITS[1], "fTreeHyp3BodyVarNclsITS1/I");
  fTree->Branch("fTreeHyp3BodyVarNclsITS2", &fTreeHyp3BodyVarNclsITS[2], "fTreeHyp3BodyVarNclsITS2/I");

  fTree->Branch("fTreeHyp3BodyVarGlobalChi20", &fTreeHyp3BodyVarGlobalChi2[0], "fTreeHyp3BodyVarGlobalChi20/F");
  fTree->Branch("fTreeHyp3BodyVarGlobalChi21", &fTreeHyp3BodyVarGlobalChi2[1], "fTreeHyp3BodyVarGlobalChi21/F");
  fTree->Branch("fTreeHyp3BodyVarGlobalChi22", &fTreeHyp3BodyVarGlobalChi2[2], "fTreeHyp3BodyVarGlobalChi22/F");

  fTree->Branch("fTreeHyp3BodyVarNsigmaTPC0", &fTreeHyp3BodyVarNsigmaTPC[0], "fTreeHyp3BodyVarNsigmaTPC0/F");
  fTree->Branch("fTreeHyp3BodyVarNsigmaTPC1", &fTreeHyp3BodyVarNsigmaTPC[1], "fTreeHyp3BodyVarNsigmaTPC1/F");
  fTree->Branch("fTreeHyp3BodyVarNsigmaTPC2", &fTreeHyp3BodyVarNsigmaTPC[2], "fTreeHyp3BodyVarNsigmaTPC2/F");

  fTree->Branch("fTreeHyp3BodyVarNsigmaTOF0", &fTreeHyp3BodyVarNsigmaTOF[0], "fTreeHyp3BodyVarNsigmaTOF0/F");
  fTree->Branch("fTreeHyp3BodyVarNsigmaTOF1", &fTreeHyp3BodyVarNsigmaTOF[1], "fTreeHyp3BodyVarNsigmaTOF1/F");
  fTree->Branch("fTreeHyp3BodyVarNsigmaTOF2", &fTreeHyp3BodyVarNsigmaTOF[2], "fTreeHyp3BodyVarNsigmaTOF2/F");

  fTree->Branch("fTreeHyp3BodyVarFlags0", &fTreeHyp3BodyVarFlags[0], "fTreeHyp3BodyVarFlags0/l");
  fTree->Branch("fTreeHyp3BodyVarFlags1", &fTreeHyp3BodyVarFlags[1], "fTreeHyp3BodyVarFlags1/l");
  fTree->Branch("fTreeHyp3BodyVarFlags2", &fTreeHyp3BodyVarFlags[2], "fTreeHyp3BodyVarFlags2/l");

  fTree->Branch("fTreeHyp3BodyVarPDGcode0", &fTreeHyp3BodyVarPDGcodes[0], "fTreeHyp3BodyVarPDGcode0/I");
  fTree->Branch("fTreeHyp3BodyVarPDGcode1", &fTreeHyp3BodyVarPDGcodes[1], "fTreeHyp3BodyVarPDGcode1/I");
  fTree->Branch("fTreeHyp3BodyVarPDGcode2", &fTreeHyp3BodyVarPDGcodes[2], "fTreeHyp3BodyVarPDGcode2/I");

  fTree->Branch("fTreeHyp3BodyVarEventId", &fTreeHyp3BodyVarEventId, "fTreeHyp3BodyVarEventId/l");
  fTree->Branch("fTreeHyp3BodyVarMotherId", &fTreeHyp3BodyVarMotherId, "fTreeHyp3BodyVarMotherId/I");

  fTree->Branch("fTreeHyp3BodyVarCandStat", &fTreeHyp3BodyVarCandStat, "fTreeHyp3BodyVarCandStat/b");

  fTree->Branch("fTreeHyp3BodyVarTruePx", &fTreeHyp3BodyVarTruePx, "fTreeHyp3BodyVarTruePx/F");
  fTree->Branch("fTreeHyp3BodyVarTruePy", &fTreeHyp3BodyVarTruePy, "fTreeHyp3BodyVarTruePy/F");
  fTree->Branch("fTreeHyp3BodyVarTruePz", &fTreeHyp3BodyVarTruePz, "fTreeHyp3BodyVarTruePz/F");

  fTree->Branch("fTreeHyp3BodyVarDecayVx", &fTreeHyp3BodyVarDecayVx, "fTreeHyp3BodyVarDecayVx/F");
  fTree->Branch("fTreeHyp3BodyVarDecayVy", &fTreeHyp3BodyVarDecayVy, "fTreeHyp3BodyVarDecayVy/F");
  fTree->Branch("fTreeHyp3BodyVarDecayVz", &fTreeHyp3BodyVarDecayVz, "fTreeHyp3BodyVarDecayVz/F");
  fTree->Branch("fTreeHyp3BodyVarDecayT", &fTreeHyp3BodyVarDecayT, "fTreeHyp3BodyVarDecayT/F");

  fTree->Branch("fTreeHyp3BodyVarPVx", &fTreeHyp3BodyVarPVx, "fTreeHyp3BodyVarPVx/F");
  fTree->Branch("fTreeHyp3BodyVarPVy", &fTreeHyp3BodyVarPVy, "fTreeHyp3BodyVarPVy/F");
  fTree->Branch("fTreeHyp3BodyVarPVz", &fTreeHyp3BodyVarPVz, "fTreeHyp3BodyVarPVz/F");
  fTree->Branch("fTreeHyp3BodyVarPVt", &fTreeHyp3BodyVarPVt, "fTreeHyp3BodyVarPVt/F");

  fTree->Branch("fTreeHyp3BodyVarMagneticField", &fTreeHyp3BodyVarMagneticField, "fTreeHyp3BodyVarMagneticField/F");

  PostData(1, fOutputList);
  PostData(2, fTree);
}

//________________________________________________________________________
void AliAnalysisTaskHypertriton3::UserExec(Option_t *) {
  // main loop called for each analized event

  AliESDEvent *esdEvent = dynamic_cast<AliESDEvent *>(InputEvent());
  if (!esdEvent) {
    ::Fatal("AliAnalysisTaskStrangenessLifetimes::UserExec", "AliESDEvent not found.");
    return;
  }

  AliMCEvent *mcEvent = MCEvent();
  if (!mcEvent) {
    ::Fatal("AliAnalysisTaskStrangenessLifetimes::UserExec", "Could not retrieve MC event");
    return;
  }

  auto ComputeRapidity = [](double rE, double rPz) {
    double rValue = -100;
    if ((rE - rPz + 1.e-13) != 0 && (rE + rPz) != 0) {
      rValue = 0.5 * TMath::Log((rE + rPz) / (rE - rPz + 1.e-13));
    }
    return rValue;
  };

  long vNTracks        = esdEvent->GetNumberOfTracks();
  float vMagneticField = esdEvent->GetMagneticField();

  // total number of analyzed events
  fHistEventCounter->Fill(0.5);

  //------------------------------------------------
  // Multiplicity Information Acquistion
  //------------------------------------------------
  float vCentrality               = 500;
  int lEvSelCode                  = 100;
  AliMultSelection *MultSelection = (AliMultSelection *)esdEvent->FindListObject("MultSelection");
  if (!MultSelection) {
    AliCentrality *centrality = 0x0;
    centrality                = esdEvent->GetCentrality();
    if (centrality) {
      vCentrality = centrality->GetCentralityPercentileUnchecked("V0M");
      lEvSelCode  = 0;
      if (centrality->GetQuality() > 1) {
        // Not good!
        lEvSelCode = 999;
      }
    }
  } else {
    // V0M Multiplicity Percentile
    vCentrality = MultSelection->GetMultiplicityPercentile("V0M");
    // Event Selection Code
    lEvSelCode = MultSelection->GetEvSelCode();
  }

  fHistCentrality->Fill(vCentrality);

  if (lEvSelCode != 0) {
    PostData(1, fOutputList);
    PostData(2, fTree);
    return;
  }

  // number of selected events
  fHistEventCounter->Fill(1.5); // selected events

  //--------------------------------------------------------------------------------
  // Part 1: fill the vector of the MC hypertritons
  //--------------------------------------------------------------------------------
  for (Int_t iPart = 0; iPart < mcEvent->GetNumberOfTracks(); iPart++) {
    AliVParticle *vPart = mcEvent->GetTrack(iPart);
    if (!vPart) {
      ::Warning("AliAnalysisTaskHyperTriton2He3piML::UserExec",
                "Generated loop %i - MC TParticle pointer to current stack particle = 0x0 ! Skipping.", iPart);
      continue;
    }
    if (mcEvent->IsPhysicalPrimary(iPart) != true) continue;

    // fill the histos of generated particles for efficiency denominator
    int vPartPDG    = vPart->PdgCode();
    double vPartPt  = vPart->Pt();
    double vPartRap = ComputeRapidity(vPart->E(), vPart->Pz());
    if (vPartPDG == 1010010030) fHistGeneratedPtVsYVsCentralityHypTrit->Fill(vPartPt, vPartRap, vCentrality);
    if (vPartPDG == -1010010030) fHistGeneratedPtVsYVsCentralityAntiHypTrit->Fill(vPartPt, vPartRap, vCentrality);
  }

  //--------------------------------------------------------------------------------
  // Part 2: establish list of tracks coming from hypertriton in the 3 body channel
  //--------------------------------------------------------------------------------
  std::vector<TrackMC> vDeuteronP;
  std::vector<TrackMC> vDeuteronM;
  std::vector<TrackMC> vProtonP;
  std::vector<TrackMC> vProtonM;
  std::vector<TrackMC> vPionP;
  std::vector<TrackMC> vPionM;
  vDeuteronP.reserve(vNTracks);
  vDeuteronM.reserve(vNTracks);
  vProtonP.reserve(vNTracks);
  vProtonM.reserve(vNTracks);
  vPionP.reserve(vNTracks);
  vPionM.reserve(vNTracks);

  for (Int_t iTrack = 0; iTrack < vNTracks; iTrack++) {
    AliESDtrack *esdTrack = esdEvent->GetTrack(iTrack);
    if (!esdTrack) continue;

    /// quality track selections
    if (((esdTrack->GetStatus() & AliVTrack::kTPCrefit) == 0) || esdTrack->GetKinkIndex(0) > 0 ||
        esdTrack->GetTPCNcls() < 70 || esdTrack->GetTPCchi2() > 4 * esdTrack->GetTPCNcls() ||
        std::fabs(esdTrack->Eta()) > 0.9)
      continue;

    int lLabel          = (int)TMath::Abs(esdTrack->GetLabel());
    AliVParticle *vPart = mcEvent->GetTrack(lLabel);

    if (IsHyperTriton3Daughter(mcEvent, vPart)) {

      int lLabelMother          = vPart->GetMother();
      AliVParticle *vMotherPart = mcEvent->GetTrack(lLabelMother);

      /// PID for hypertriton daughter tracks
      float nSigmaDeu = fPIDResponse->NumberOfSigmasTPC(esdTrack, AliPID::kDeuteron);
      float nSigmaP   = fPIDResponse->NumberOfSigmasTPC(esdTrack, AliPID::kProton);
      float nSigmaPi  = fPIDResponse->NumberOfSigmasTPC(esdTrack, AliPID::kPion);

      TrackMC tmc{esdTrack, vMotherPart, vPart, lLabelMother, iTrack};

      float trackSign = esdTrack->GetSign();

      if (trackSign > 0) {
        if (nSigmaDeu < 5.) vDeuteronP.push_back(tmc);
        if (nSigmaP < 5.) vProtonP.push_back(tmc);
        if (nSigmaPi < 5.) vPionP.push_back(tmc);
      }

      if (trackSign < 0) {
        if (nSigmaDeu < 5.) vDeuteronM.push_back(tmc);
        if (nSigmaP < 5.) vProtonM.push_back(tmc);
        if (nSigmaPi < 5.) vPionM.push_back(tmc);
      }
    }
  }

  //--------------------------------------------------------------------------------
  // Part 3: find the triplets of reconstructed daughters and fill the tree
  //--------------------------------------------------------------------------------
  std::vector<CandidateMC> candidates;
  candidates.reserve(vNTracks);

  fTreeHyp3BodyVarMagneticField = vMagneticField;
  fTreeHyp3BodyVarEventId++;
  fTreeHyp3BodyVarPVt = vDeuteronP.back().mother->Tv();
  fTreeHyp3BodyVarPVx = vDeuteronP.back().mother->Xv();
  fTreeHyp3BodyVarPVy = vDeuteronP.back().mother->Yv();
  fTreeHyp3BodyVarPVz = vDeuteronP.back().mother->Zv();

  for (size_t iTrack = 0; iTrack < vDeuteronP.size(); iTrack++) {
    for (size_t jTrack = 0; jTrack < vProtonP.size(); jTrack++) {
      /// consider tracks from the same mother only: part 1
      if (vDeuteronP[iTrack].motherId != vProtonP[jTrack].motherId) continue;
      /// reject candidates with the same track: part 1
      if (vDeuteronP[iTrack].partId == vProtonP[jTrack].partId) continue;
      for (size_t zTrack = 0; zTrack < vPionM.size(); zTrack++) {
        /// consider tracks from the same mother only: part 2
        if (vDeuteronP[iTrack].motherId != vPionM[zTrack].motherId) continue;
        /// reject candidates with the same track: part 2
        if (vDeuteronP[iTrack].partId == vPionM[zTrack].partId || vProtonP[jTrack].partId == vPionM[zTrack].partId)
          continue;

        /// define the status of the reconstructed candidate
        unsigned char stat  = 0u;
        bool goodCand       = false;
        bool reflectionCand = false;

        goodCand = !IsFakeCandidate(mcEvent, vDeuteronP[iTrack].motherId, vDeuteronP[iTrack].particle,
                                    vProtonP[jTrack].particle, vPionM[zTrack].particle);
        if (!goodCand)
          reflectionCand = IsReflectionCandidate(mcEvent, vDeuteronP[iTrack].motherId, vDeuteronP[iTrack].particle,
                                                 vProtonP[jTrack].particle, vPionM[zTrack].particle);

        if (goodCand) stat |= g;
        if (reflectionCand) stat |= r;

        CandidateMC c;
        c.track_deu = vDeuteronP[iTrack].track;
        c.track_p   = vProtonP[jTrack].track;
        c.track_pi  = vPionM[zTrack].track;
        c.part1     = vDeuteronP[iTrack].particle;
        c.part2     = vProtonP[jTrack].particle;
        c.part3     = vPionM[zTrack].particle;
        c.mother    = vDeuteronP[iTrack].mother;
        c.motherId  = vDeuteronP[iTrack].motherId;
        c.status    = stat;
        candidates.push_back(c);
      }
    }
  }

  for (size_t iTrack = 0; iTrack < vDeuteronM.size(); iTrack++) {
    for (size_t jTrack = 0; jTrack < vProtonM.size(); jTrack++) {
      /// consider tracks from the same mother only: part 1
      if (vDeuteronM[iTrack].motherId != vProtonM[jTrack].motherId) continue;
      /// reject candidates with the same track: part 1
      if (vDeuteronM[iTrack].partId == vProtonM[jTrack].partId) continue;
      for (size_t zTrack = 0; zTrack < vPionP.size(); zTrack++) {
        /// consider tracks from the same mother only: part 2
        if (vDeuteronM[iTrack].motherId != vPionP[zTrack].motherId) continue;
        /// reject candidates with the same track: part 2
        if (vDeuteronM[iTrack].partId == vPionP[zTrack].partId || vProtonM[jTrack].partId == vPionP[zTrack].partId)
          continue;

        /// define the status of the reconstructed candidate
        unsigned char stat  = 0u;
        bool goodCand       = false;
        bool reflectionCand = false;

        goodCand = !IsFakeCandidate(mcEvent, vDeuteronM[iTrack].motherId, vDeuteronM[iTrack].particle,
                                    vProtonM[jTrack].particle, vPionP[zTrack].particle);
        if (!goodCand)
          reflectionCand = IsReflectionCandidate(mcEvent, vDeuteronM[iTrack].motherId, vDeuteronM[iTrack].particle,
                                                 vProtonM[jTrack].particle, vPionP[zTrack].particle);

        if (goodCand) stat |= g;
        if (reflectionCand) stat |= r;

        CandidateMC c;
        c.track_deu = vDeuteronM[iTrack].track;
        c.track_p   = vProtonM[jTrack].track;
        c.track_pi  = vPionP[zTrack].track;
        c.part1     = vDeuteronM[iTrack].particle;
        c.part2     = vProtonM[jTrack].particle;
        c.part3     = vPionP[zTrack].particle;
        c.mother    = vDeuteronM[iTrack].mother;
        c.motherId  = vDeuteronM[iTrack].motherId;
        c.status    = stat;
        candidates.push_back(c);
      }
    }
  }

  /// sorting hypertriton candidates respect the motherId
  std::sort(candidates.begin(), candidates.end(),
            [](const CandidateMC &a, const CandidateMC &b) { return a.motherId < b.motherId; });

  for (auto cand : candidates) {
    /// fill the tree of findable
    fTreeHyp3BodyVarTracks[0] = static_cast<AliExternalTrackParam *>(cand.track_deu);
    fTreeHyp3BodyVarTracks[1] = static_cast<AliExternalTrackParam *>(cand.track_p);
    fTreeHyp3BodyVarTracks[2] = static_cast<AliExternalTrackParam *>(cand.track_pi);

    fTreeHyp3BodyVarNclsTPC[0] = (Int_t)cand.track_deu->GetTPCNcls();
    fTreeHyp3BodyVarNclsTPC[1] = (Int_t)cand.track_p->GetTPCNcls();
    fTreeHyp3BodyVarNclsTPC[2] = (Int_t)cand.track_pi->GetTPCNcls();

    fTreeHyp3BodyVarNclsITS[0] = (Int_t)cand.track_deu->GetITSNcls();
    fTreeHyp3BodyVarNclsITS[1] = (Int_t)cand.track_p->GetITSNcls();
    fTreeHyp3BodyVarNclsITS[2] = (Int_t)cand.track_pi->GetITSNcls();

    fTreeHyp3BodyVarGlobalChi2[0] = (Float_t)cand.track_deu->GetGlobalChi2();
    fTreeHyp3BodyVarGlobalChi2[1] = (Float_t)cand.track_p->GetGlobalChi2();
    fTreeHyp3BodyVarGlobalChi2[2] = (Float_t)cand.track_pi->GetGlobalChi2();

    fTreeHyp3BodyVarNsigmaTPC[0] = std::abs(fPIDResponse->NumberOfSigmasTPC(cand.track_deu, AliPID::kDeuteron));
    fTreeHyp3BodyVarNsigmaTPC[1] = std::abs(fPIDResponse->NumberOfSigmasTPC(cand.track_p, AliPID::kProton));
    fTreeHyp3BodyVarNsigmaTPC[2] = std::abs(fPIDResponse->NumberOfSigmasTPC(cand.track_pi, AliPID::kPion));

    HasTOF(cand.track_deu)
        ? fTreeHyp3BodyVarNsigmaTOF[0] = fPIDResponse->NumberOfSigmasTOF(cand.track_deu, AliPID::kDeuteron)
        : fTreeHyp3BodyVarNsigmaTOF[0] = -1.0;
    HasTOF(cand.track_p)
        ? fTreeHyp3BodyVarNsigmaTOF[1] = fPIDResponse->NumberOfSigmasTOF(cand.track_p, AliPID::kProton)
        : fTreeHyp3BodyVarNsigmaTOF[1] = -1.0;
    HasTOF(cand.track_pi)
        ? fTreeHyp3BodyVarNsigmaTOF[2] = fPIDResponse->NumberOfSigmasTOF(cand.track_pi, AliPID::kPion)
        : fTreeHyp3BodyVarNsigmaTOF[2] = -1.0;

    fTreeHyp3BodyVarFlags[0] = (ULong64_t)cand.track_deu->GetStatus();
    fTreeHyp3BodyVarFlags[1] = (ULong64_t)cand.track_p->GetStatus();
    fTreeHyp3BodyVarFlags[2] = (ULong64_t)cand.track_pi->GetStatus();

    fTreeHyp3BodyVarPDGcodes[0] = cand.part1->PdgCode();
    fTreeHyp3BodyVarPDGcodes[1] = cand.part2->PdgCode();
    fTreeHyp3BodyVarPDGcodes[2] = cand.part3->PdgCode();

    AliVParticle *vHyperTriton = cand.mother;
    fTreeHyp3BodyVarTruePx     = vHyperTriton->Px();
    fTreeHyp3BodyVarTruePy     = vHyperTriton->Py();
    fTreeHyp3BodyVarTruePz     = vHyperTriton->Pz();

    AliVParticle *vProng    = cand.part1;
    fTreeHyp3BodyVarDecayVx = vProng->Xv();
    fTreeHyp3BodyVarDecayVy = vProng->Yv();
    fTreeHyp3BodyVarDecayVz = vProng->Zv();
    fTreeHyp3BodyVarDecayT  = vProng->Tv();

    fTreeHyp3BodyVarMotherId = cand.motherId;

    fTreeHyp3BodyVarCandStat = cand.status;
    fTree->Fill();
  }

  PostData(1, fOutputList);
  PostData(2, fTree);
  return;
}

//________________________________________________________________________
void AliAnalysisTaskHypertriton3::Terminate(Option_t *) {
  // Merge output
  // Called once at the end of the query

  fOutputList = dynamic_cast<TList *>(GetOutputData(1));
  if (!fOutputList) {
    printf("ERROR: fOutputList not available\n");
    return;
  }

  printf("end of Terminate");
  return;

} // end of Terminate