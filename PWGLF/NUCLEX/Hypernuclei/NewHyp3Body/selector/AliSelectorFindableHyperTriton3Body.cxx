#define AliSelectorFindableHyperTriton3Body_cxx

#include "AliSelectorFindableHyperTriton3Body.h"
#include "AliESDVertex.h"
#include "AliESDtrack.h"
#include "AliExternalTrackParam.h"
#include "AliPID.h"
#include "AliVTrack.h"
#include <TH1D.h>
#include <TH2D.h>
#include <TList.h>
#include <TLorentzVector.h>
#include <TObjArray.h>
#include <TStyle.h>
#include <TVector3.h>
#include <cmath>

AliSelectorFindableHyperTriton3Body::AliSelectorFindableHyperTriton3Body(TString outputName, TString outputPath,
                                                                         TTree *)
    : fOutputFileName{outputName}, fOutputFilePath{outputPath}, fHypertritonVertexer(), fHypertritonVertexerHard() {
  fHypertritonVertexerHard.SetToleranceGuessCompatibility(2);
  fHypertritonVertexerHard.SetMaxDinstanceInit(1.5);
}

//______________________________________________________________________________
void AliSelectorFindableHyperTriton3Body::Begin(TTree * /*tree*/) {
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).

  TString option = GetOption();
}

//______________________________________________________________________________
void AliSelectorFindableHyperTriton3Body::SlaveBegin(TTree * /*tree*/) {
  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).

  TString option = GetOption();

  // create and initialize PIDResponse object
  fPIDResponse = new AliPIDResponse(kTRUE);
  fPIDResponse->SetOADBPath("~/alice/AliPhysics/OADB/COMMON/PID/MC/PassNameMap.root");
  fPIDResponse->SetMCperiod("LHC16h7a");
  fPIDResponse->SetRecoPass(1);

  /// check histograms
  fHistTrackCheck[0] = new TH1D("fHistAllTrackCheck", "", 2, -1, 1);
  GetOutputList()->Add(fHistTrackCheck[0]);
  fHistTrackCheck[1] = new TH1D("fHistDeuTrackCheck", "", 2, -1, 1);
  GetOutputList()->Add(fHistTrackCheck[1]);
  fHistTrackCheck[2] = new TH1D("fHistPTrackCheck", "", 2, -1, 1);
  GetOutputList()->Add(fHistTrackCheck[2]);
  fHistTrackCheck[3] = new TH1D("fHistPiTrackCheck", "", 2, -1, 1);
  GetOutputList()->Add(fHistTrackCheck[3]);

  fHistPDGCheck[0] = new TH1D("fHistPDG_all", "", 2, -1, 1);
  GetOutputList()->Add(fHistPDGCheck[0]);
  fHistPDGCheck[1] = new TH1D("fHistPDG_deu", "", 2, -1, 1);
  GetOutputList()->Add(fHistPDGCheck[1]);
  fHistPDGCheck[2] = new TH1D("fHistPDG_p", "", 2, -1, 1);
  GetOutputList()->Add(fHistPDGCheck[2]);
  fHistPDGCheck[3] = new TH1D("fHistPDG_pi", "", 2, -1, 1);
  GetOutputList()->Add(fHistPDGCheck[3]);
  fHistPDG = new TH2D("fHistDaughterPDG", "", 8, 0, 8, 3, 0, 3);
  GetOutputList()->Add(fHistPDG);

  fHistChargeCheck = new TH1D("fHistChargeCheck", "", 2, -1, 1);
  GetOutputList()->Add(fHistChargeCheck);

  fHistPassCheck = new TH1D("fHistCheckPass", "", 6, -1, 5);
  GetOutputList()->Add(fHistPassCheck);

  fHistClonesCheck = new TH1D("fHistCheckCVSizeTrack", "", 6, -1, 5);
  GetOutputList()->Add(fHistClonesCheck);

  fHistSameTrackCheck[0] = new TH1D("fHistSameTrackDeu", "", 2, -1, 1);
  GetOutputList()->Add(fHistSameTrackCheck[0]);
  fHistSameTrackCheck[1] = new TH1D("fHistSameTrackP", "", 2, -1, 1);
  GetOutputList()->Add(fHistSameTrackCheck[1]);
  fHistSameTrackCheck[2] = new TH1D("fHistSameTrackPi", "", 2, -1, 1);
  GetOutputList()->Add(fHistSameTrackCheck[2]);
  fHistSameTrackCheck[3] = new TH1D("fHistSameTrackAll", "", 2, -1, 1);
  GetOutputList()->Add(fHistSameTrackCheck[3]);

  fHistVertexChi2 = new TH1D("fHistVertexChi2", ";Vertex #chi^{2}/NDF; Counts", 500, 0, 100);
  GetOutputList()->Add(fHistVertexChi2);

  const char lAM[3]{"AM"};
  const char *lSpecies[3]{"d", "p", "pi"};
  const char *lSpeciesPair[3]{"dp", "ppi", "pid"};
  const char *lCuts[4]{"ref", "std", "cuts", "hard"};
  const char lCoords[4]{"xyz"};

  for (int iMatter = 0; iMatter < 2; iMatter++) {
    for (int iCuts = 0; iCuts < 4; iCuts++) {
      fHistInvMass[iMatter][iCuts] = new TH2D(Form("fHistInvMass_%c_%s", lAM[iMatter], lCuts[iCuts]),
                                              ";m (dp#pi) [GeV/#it{c}^{2}];Counts", 200, 2.95, 3.35, 20, 0, 10);
      fHistPt[iMatter][iCuts] =
          new TH1D(Form("fHistPt_%c_%s", lAM[iMatter], lCuts[iCuts]), ";#it{p}_{T} [GeV/#it{c}];Counts", 100, 0, 10);
      GetOutputList()->Add(fHistInvMass[iMatter][iCuts]);
      GetOutputList()->Add(fHistPt[iMatter][iCuts]);
    }
  }
  for (int iCoord = 0; iCoord < 3; iCoord++) {
    fHistResDecayVtx[iCoord] =
        new TH2D(Form("fHistResDecayVtx%c", lCoords[iCoord]),
                 Form(";#Delta%c [mm];#it{p}_{T} [GeV/#it{c}]", lCoords[iCoord]), 500, -10, 10, 10, 0, 10);
    GetOutputList()->Add(fHistResDecayVtx[iCoord]);
  }

  for (int iCoord = 0; iCoord < 3; iCoord++) {
    fHistTrackDistance[iCoord] =
        new TH1D(Form("fTrackDistance_%s-%s", lSpeciesPair[iCoord], lSpeciesPair[(iCoord + 1) % 3]), ";distance [mm]",
                 500, 0, 2000);
    GetOutputList()->Add(fHistTrackDistance[iCoord]);
  }

  for (int iSpecies = 0; iSpecies < 3; iSpecies++) {
    fHistDCA2pV[iSpecies] =
        new TH2D(Form("fDCA2Primaryvtx_%s", lSpecies[iSpecies]), ";mm ;#it{p}_{T} [GeV/#it{c}]", 600, 0, 30, 10, 0, 10);
    fHistDCA2dV[iSpecies] =
        new TH2D(Form("fDCA2Decayvtx_%s", lSpecies[iSpecies]), ";mm ;#it{p}_{T} [GeV/#it{c}]", 600, 0, 30, 10, 0, 10);
    GetOutputList()->Add(fHistDCA2pV[iSpecies]);
    GetOutputList()->Add(fHistDCA2dV[iSpecies]);
  }

  for (int iSpecies = 0; iSpecies < 3; iSpecies++) {
    fHistNSigma[iSpecies]      = new TH1D(Form("fHistNSigma_%s", lSpecies[iSpecies]), "n#sigma", 10, 0, 10);
    fHistNSigmaCheck[iSpecies] = new TH1D(Form("fHistNSigmaCheck_%s", lSpecies[iSpecies]), "", 2, -1, 1);
    GetOutputList()->Add(fHistNSigma[iSpecies]);
    GetOutputList()->Add(fHistNSigmaCheck[iSpecies]);
  }
  fHistNSigmaCheck[3] = new TH1D("fHistNSigmaCheckAll", "", 2, -1, 1);
  GetOutputList()->Add(fHistNSigmaCheck[3]);

  for (int iSpecies = 0; iSpecies < 3; iSpecies++) {
    fHistNSigmaFakeCheck[iSpecies] = new TH1D(Form("fHistNSigmaFakeCheck_%s", lSpecies[iSpecies]), "", 2, -1, 1);
    GetOutputList()->Add(fHistNSigmaFakeCheck[iSpecies]);
  }
  fHistNSigmaFakeCheck[3] = new TH1D("fHistNSigmaFakeCheckAll", "", 2, -1, 1);
  GetOutputList()->Add(fHistNSigmaFakeCheck[3]);

  fHistCosPAngle = new TH1D("fCosPointingAngle", ";#it{cos#theta_{pointing}} ;Counts", 5000, 0.5, 1.);
  GetOutputList()->Add(fHistCosPAngle);
}

//______________________________________________________________________________
Bool_t AliSelectorFindableHyperTriton3Body::Process(Long64_t entry) {
  // The Process() function is called for each entry in the tree (or possibly
  // keyed object in the case of PROOF) to be processed. The entry argument
  // specifies which entry in the currently loaded tree is to be processed.
  // When processing keyed objects with PROOF, the object is already loaded
  // and is available via the fObject pointer.
  //
  // This function should contain the \"body\" of the analysis. It can contain
  // simple or elaborate selection criteria, run algorithms on the data
  // of the event and typically fill histograms.
  //
  // The processing can be stopped by calling Abort().
  //
  // Use fStatus to set the return value of TTree::Process().
  //
  // The return value is currently not used.

  fReader.SetEntry(entry);

  // define support stuff
  bool lIsGoodTrack[4]      = {false};
  bool lIsGoodCharge        = false;
  bool lIsGoodPDG[4]        = {false};
  bool lIsGoodPID[4]        = {false};
  bool lIsGoodVertex        = false;
  bool lIsGoodPAngle        = false;
  bool lIsGoodCandidate     = false;
  bool lIsGoodVertexHard    = false;
  bool lIsGoodPAngleHard    = false;
  bool lIsGoodCandidateHard = false;

  float nSigma[3] = {0., 0., 0.};

  double fPrimaryVertexPos[3]   = {0., 0., 0.};
  double fPrimaryVertexCov[6]   = {0.};
  double fDecayVertexPos[3]     = {0., 0., 0.};
  double fDecayVertexPosHard[3] = {0., 0., 0.};
  double fDecayLenght[3]        = {0., 0., 0.};
  double fDecayLenghtHard[3]    = {0., 0., 0.};
  double fDCA2dv[3][3]          = {{0.}};
  double fDCA2pv[3][3]          = {{0.}};
  double fDCA2dvcov[3][6]       = {{0.}};
  double fDCA2pvcov[3][6]       = {{0.}};
  TLorentzVector fHypertriton   = {0., 0., 0., 0.};
  TLorentzVector fLVDaughter[3];

  AliExternalTrackParam fParamTrack[3] = {*fTreeHyp3BodyVarTracks[0], *fTreeHyp3BodyVarTracks[1],
                                          *fTreeHyp3BodyVarTracks[2]};
  const float lMasses[3]{AliPID::ParticleMass(AliPID::kDeuteron), AliPID::ParticleMass(AliPID::kProton),
                         AliPID::ParticleMass(AliPID::kPion)};
  const int pdgCodes[6]{AliPID::ParticleCode(AliPID::kElectron), AliPID::ParticleCode(AliPID::kMuon),
                        AliPID::ParticleCode(AliPID::kPion),     AliPID::ParticleCode(AliPID::kKaon),
                        AliPID::ParticleCode(AliPID::kProton),   AliPID::ParticleCode(AliPID::kDeuteron)};
  AliPID::EParticleType lSpecNum[3]{AliPID::kDeuteron, AliPID::kProton, AliPID::kPion};

  AliESDtrack *t1 = new AliESDtrack(*fTreeHyp3BodyVarTracks[0]);
  AliESDtrack *t2 = new AliESDtrack(*fTreeHyp3BodyVarTracks[1]);
  AliESDtrack *t3 = new AliESDtrack(*fTreeHyp3BodyVarTracks[2]);
  AliESDtrack *t[3]{t1, t2, t3};

  //______________________________________________
  /// clones check
  bool lIsPossibleClone = (fLastMother >= 0) && (fCurrentEventId == *fTreeHyp3BodyVarEventId) && (fLastMother == *fTreeHyp3BodyVarMotherId);
  fCurrentEventId = *fTreeHyp3BodyVarEventId;
  fLastMother     = *fTreeHyp3BodyVarMotherId;

  // create candidate object
  HyperTritonCandidate fCand(t1, t2, t3);
  fCand.fPDG[0]          = *fTreeHyp3BodyVarPDGcodes[0];
  fCand.fPDG[1]          = *fTreeHyp3BodyVarPDGcodes[1];
  fCand.fPDG[2]          = *fTreeHyp3BodyVarPDGcodes[2];
  fCand.fPrimaryVtx[0]   = *fTreeHyp3BodyVarPVtx[0];
  fCand.fPrimaryVtx[1]   = *fTreeHyp3BodyVarPVtx[1];
  fCand.fPrimaryVtx[2]   = *fTreeHyp3BodyVarPVtx[2];
  fCand.fPrimaryVtx[3]   = *fTreeHyp3BodyVarPVtx[3];
  fCand.fTrueDecayVtx[0] = *fTreeHyp3BodyVarDecayVtx[0];
  fCand.fTrueDecayVtx[1] = *fTreeHyp3BodyVarDecayVtx[1];
  fCand.fTrueDecayVtx[2] = *fTreeHyp3BodyVarDecayVtx[2];
  fCand.fTrueDecayVtx[3] = *fTreeHyp3BodyVarDecayVtx[3];
  fCand.fTrueP[0]        = *fTreeHyp3BodyVarTrueP[0];
  fCand.fTrueP[1]        = *fTreeHyp3BodyVarTrueP[1];
  fCand.fTrueP[2]        = *fTreeHyp3BodyVarTrueP[2];
  fCand.fMagField        = *fTreeHyp3BodyVarMagneticField;
  fCand.fIsClone         = lIsPossibleClone;

  /// pdg check
  int pdg_sum = 0;
  for (int iS = 0; iS < 3; iS++) {
    lIsGoodPDG[iS] = (std::abs(*fTreeHyp3BodyVarPDGcodes[iS]) == AliPID::ParticleCode(lSpecNum[iS]));
    fHistPDGCheck[iS+1]->Fill(lIsGoodPDG[iS] - 0.5);
    pdg_sum += lIsGoodPDG[iS];
  }
  bool lIsHyp = (*fTreeHyp3BodyVarPDGcodes[0] == 1000010020) && (*fTreeHyp3BodyVarPDGcodes[1] == 2212) &&
                (*fTreeHyp3BodyVarPDGcodes[2] == -211);
  bool lIsAntiHyp = (*fTreeHyp3BodyVarPDGcodes[0] == -1000010020) && (*fTreeHyp3BodyVarPDGcodes[1] == -2212) &&
                    (*fTreeHyp3BodyVarPDGcodes[2] == 211);
  lIsGoodPDG[3] = (lIsHyp || lIsAntiHyp);
  fHistPDGCheck[0]->Fill(lIsGoodPDG[3] - 0.5);

  // findalble pdg check
  for (int iD = 0; iD < 3; iD++) {
    for (int iP = 0; iP < 6; iP++) {
      if (std::abs(*fTreeHyp3BodyVarPDGcodes[iD]) == pdgCodes[iP]) fHistPDG->Fill(iP, iD);
    }
  }

  /// track cuts check
  int trk_sum = 0;
  for (int iS = 0; iS < 3; iS++) {
    lIsGoodTrack[iS] = IsTrackCutsGood(t[iS]);
    fHistTrackCheck[iS + 1]->Fill(lIsGoodTrack[iS] - 0.5);
    trk_sum += lIsGoodPDG[iS];
  }
  lIsGoodTrack[3] = (trk_sum == 3);
  fHistTrackCheck[0]->Fill(lIsGoodTrack[3] - 0.5);

  // PID selection on candidates
  if (lIsGoodPDG[3] && lIsGoodTrack[3]) {
    int sum = 0;
    for (int iS = 0; iS < 3; iS++) {
      nSigma[iS]     = std::abs(fPIDResponse->NumberOfSigmasTPC(t[iS], lSpecNum[iS]));
      lIsGoodPID[iS] = nSigma[iS] < 3.;
      fHistNSigmaCheck[iS]->Fill(lIsGoodPID[iS] - 0.5);
      fHistNSigma[iS]->Fill(nSigma[iS]);
      sum += lIsGoodPID[iS];
    }
    lIsGoodPID[3] = (sum == 3);
    fHistNSigmaCheck[3]->Fill(lIsGoodPID[3] - 0.5);
  }
  // PID selection on fake candidates
  if (!lIsGoodPDG && lIsGoodTrack) {
    int sum = 0;
    for (int iS = 0; iS < 3; iS++) {
      if (std::abs(*fTreeHyp3BodyVarPDGcodes[iS]) != AliPID::ParticleCode(lSpecNum[iS])) {
        nSigma[iS]     = std::abs(fPIDResponse->NumberOfSigmasTPC(t[iS], lSpecNum[iS]));
        lIsGoodPID[iS] = nSigma[iS] < 3.;
        fHistNSigmaFakeCheck[iS]->Fill(lIsGoodPID[iS] - 0.5);
        sum += lIsGoodPID[iS];
      }
    }
    bool rejectedCandidate = (sum < 3);
    fHistNSigmaFakeCheck[3]->Fill(rejectedCandidate - 0.5);
  }

  /// charge check
  if ((t1->Charge() == 1) && (t2->Charge() == 1) && (t3->Charge() == -1)) lIsGoodCharge = true;
  if ((t1->Charge() == -1) && (t2->Charge() == -1) && (t3->Charge() == 1)) lIsGoodCharge = true;
  fHistChargeCheck->Fill(lIsGoodCharge - 0.5);

  // create the TLorentzVector of the hyper-triton candidate
  for (int iTrack = 0; iTrack < 3; iTrack++) {
    fLVDaughter[iTrack].SetXYZM(fParamTrack[iTrack].Px(), fParamTrack[iTrack].Py(), fParamTrack[iTrack].Pz(),
                                lMasses[iTrack]);
    fHypertriton += fLVDaughter[iTrack];
  }
  double hypMass = fHypertriton.M();
  double hypPt   = fHypertriton.Pt();

  // total number of findable candidates
  fHistPassCheck->Fill(-0.5);

  if (lIsGoodTrack[3] && lIsGoodPDG[3] && lIsGoodPID[3] && lIsGoodCharge) {
    fHistPassCheck->Fill(0.5);
    // reconstruct the decay vertex with the dedicated vertexer
    bool recoVertex     = fHypertritonVertexer.FindDecayVertex(&fCand.fDaughterTrack[0], &fCand.fDaughterTrack[1],
                                                           &fCand.fDaughterTrack[2], fCand.fMagField);
    bool recoVertexHard = fHypertritonVertexerHard.FindDecayVertex(&fCand.fDaughterTrack[0], &fCand.fDaughterTrack[1],
                                                                   &fCand.fDaughterTrack[2], fCand.fMagField);
    AliESDVertex *fDecayVertex     = dynamic_cast<AliESDVertex *>(fHypertritonVertexer.GetCurrentVertex());
    AliESDVertex *fDecayVertexHard = dynamic_cast<AliESDVertex *>(fHypertritonVertexerHard.GetCurrentVertex());
    if (recoVertex) {
      double vertexChi2NDF = fDecayVertex->GetChi2perNDF() * 1.;
      fHistVertexChi2->Fill(vertexChi2NDF);
    }
    // possibility to constraint the goodness of the reconstructed decay vertex
    // lIsGoodVertex = recoVertex && (fDecayVertex->GetChi2perNDF() < 20.);
    lIsGoodVertex     = recoVertex;
    lIsGoodVertexHard = recoVertexHard;

    // compute decay vertex position and compute cos(pointingangle)
    if (lIsGoodVertex) {
      fHistPassCheck->Fill(1.5);
      fDecayVertex->GetXYZ(fDecayVertexPos);
      for (int iTrack = 0; iTrack < 3; iTrack++) {
        fHistResDecayVtx[iTrack]->Fill((fDecayVertexPos[iTrack] - *fTreeHyp3BodyVarDecayVtx[iTrack]) * 10., hypPt);
        fDecayLenght[iTrack] = fDecayVertexPos[iTrack] - fCand.fPrimaryVtx[iTrack];
      }
      TVector3 v(fDecayLenght[0], fDecayLenght[1], fDecayLenght[2]);
      float pointAngle = fHypertriton.Angle(v);
      float cospa      = std::cos(pointAngle);
      fHistCosPAngle->Fill(cospa);
      lIsGoodPAngle    = (cospa > 0.80);
      lIsGoodCandidate = lIsGoodPAngle;
      if (lIsGoodPAngle) fHistPassCheck->Fill(2.5);

      // compute the DCA of the 3 tracks from the primary and decay vertex
      fPrimaryVertexPos[0] = *fTreeHyp3BodyVarPVtx[0] * 1.;
      fPrimaryVertexPos[1] = *fTreeHyp3BodyVarPVtx[1] * 1.;
      fPrimaryVertexPos[2] = *fTreeHyp3BodyVarPVtx[2] * 1.;
      AliESDVertex fPV(fPrimaryVertexPos, fPrimaryVertexCov, 1., 1000);
      for (int iTrack = 0; iTrack < 3; iTrack++) {
        fParamTrack[iTrack].PropagateToDCA(fDecayVertex, fCand.fMagField, 1000., fDCA2dv[iTrack], fDCA2dvcov[iTrack]);
        fParamTrack[iTrack].PropagateToDCA(&fPV, fCand.fMagField, 1000., fDCA2pv[iTrack], fDCA2pvcov[iTrack]);
        float ddv = std::sqrt(fDCA2dv[iTrack][0] * fDCA2dv[iTrack][0] + fDCA2dv[iTrack][1] * fDCA2dv[iTrack][1] +
                              fDCA2dv[iTrack][2] * fDCA2dv[iTrack][2]);
        float dpv = std::sqrt(fDCA2pv[iTrack][0] * fDCA2pv[iTrack][0] + fDCA2pv[iTrack][1] * fDCA2pv[iTrack][1] +
                              fDCA2pv[iTrack][2] * fDCA2pv[iTrack][2]);
        fHistDCA2dV[iTrack]->Fill(ddv * 10., hypPt);
        fHistDCA2pV[iTrack]->Fill(dpv * 10., hypPt);
      }

      // compute the track2track distance used in the vertexer
      float pPM[3][3];
      for (int iTrack = 0; iTrack < 3; iTrack++) {
        fHypertritonVertexer.Find2ProngClosestPoint(&fParamTrack[iTrack], &fParamTrack[(iTrack + 1) % 3],
                                                    fCand.fMagField, pPM[iTrack]);
      }

      for (int iPerm = 0; iPerm < 3; iPerm++) {
        float d = 0.;
        for (int iDim = 0; iDim < 3; ++iDim) {
          d += (pPM[iPerm][iDim] - pPM[(iPerm + 1) % 3][iDim]) * (pPM[iPerm][iDim] - pPM[(iPerm + 1) % 3][iDim]);
        }
        fHistTrackDistance[iPerm]->Fill(std::sqrt(d) * 10.);
      }
    }
    // again with harder selections
    if (lIsGoodVertexHard) {
      fDecayVertexHard->GetXYZ(fDecayVertexPosHard);
      for (int iTrack = 0; iTrack < 3; iTrack++) {
        fDecayLenghtHard[iTrack] = fDecayVertexPosHard[iTrack] - fCand.fPrimaryVtx[iTrack];
      }
      TVector3 v(fDecayLenghtHard[0], fDecayLenghtHard[1], fDecayLenghtHard[2]);
      float pointAngle     = fHypertriton.Angle(v);
      float cospa          = std::cos(pointAngle);
      lIsGoodPAngleHard    = (cospa > 0.98);
      lIsGoodCandidateHard = lIsGoodPAngleHard;
    }
  }

  // add PID selection
  lIsGoodCandidate     = lIsGoodCandidate && lIsGoodPID[3];
  lIsGoodCandidateHard = lIsGoodCandidateHard && lIsGoodPID[3];

  if (lIsGoodPDG[3] && lIsGoodCharge && !lIsPossibleClone) {
    // invariant mass distribution for efficiency
    if (lIsAntiHyp) fHistInvMass[0][0]->Fill(hypMass, hypPt);
    if (lIsHyp) fHistInvMass[1][0]->Fill(hypMass, hypPt);
    if (lIsGoodTrack[3] && lIsGoodPID[3] && lIsAntiHyp) fHistInvMass[0][1]->Fill(hypMass, hypPt);
    if (lIsGoodCandidate && lIsAntiHyp) fHistInvMass[0][2]->Fill(hypMass, hypPt);
    if (lIsGoodCandidateHard && lIsAntiHyp) fHistInvMass[0][3]->Fill(hypMass, hypPt);
    if (lIsGoodTrack[3] && lIsGoodPID[3] && lIsHyp) fHistInvMass[1][1]->Fill(hypMass, hypPt);
    if (lIsGoodCandidate && lIsHyp) fHistInvMass[1][2]->Fill(hypMass, hypPt);
    if (lIsGoodCandidateHard && lIsHyp) fHistInvMass[1][3]->Fill(hypMass, hypPt);

    // transverse momentum distribution
    fHistPt[lIsHyp][0]->Fill(hypPt);
    if (lIsGoodTrack[3 && lIsGoodPID[3]]) fHistPt[lIsHyp][1]->Fill(hypPt);
    if (lIsGoodCandidate) fHistPt[lIsHyp][2]->Fill(hypPt);
    if (lIsGoodCandidateHard) fHistPt[lIsHyp][3]->Fill(hypPt);
  }

  // clones check
  if (fClonesVector.size() == 0) {
    if (lIsGoodCandidate && lIsGoodPDG[3] && lIsGoodTrack[3] && lIsGoodCharge) fClonesVector.push_back(fCand);
  } else {
    if (lIsPossibleClone) {
      if (lIsGoodCandidate && lIsGoodPDG[3] && lIsGoodTrack[3] && lIsGoodCharge) fClonesVector.push_back(fCand);
    } else {
      fHistClonesCheck->Fill(fClonesVector.size() - 0.5);
      // loop sul vettore qui
      std::vector<AliESDtrack> fDeuTrack;
      std::vector<AliESDtrack> fProTrack;
      std::vector<AliESDtrack> fPioTrack;
      for (const auto cand : fClonesVector) {
        AliESDtrack d(cand.fDaughterTrack[0]);
        AliESDtrack p(cand.fDaughterTrack[1]);
        AliESDtrack pi(cand.fDaughterTrack[2]);
        fDeuTrack.push_back(d);
        fProTrack.push_back(p);
        fPioTrack.push_back(pi);
      }
      if (fClonesVector.size() == 2) {
        bool sameAll = IsSameTrack(&fDeuTrack[0], &fDeuTrack[1]) && IsSameTrack(&fProTrack[0], &fProTrack[1]) &&
                       IsSameTrack(&fPioTrack[0], &fPioTrack[1]);
        bool sameDeu = IsSameTrack(&fDeuTrack[0], &fDeuTrack[1]);
        bool sameP   = IsSameTrack(&fProTrack[0], &fProTrack[1]);
        bool samePi  = IsSameTrack(&fPioTrack[0], &fPioTrack[1]);
        fHistSameTrackCheck[0]->Fill(sameDeu - 0.5);
        fHistSameTrackCheck[1]->Fill(sameP - 0.5);
        fHistSameTrackCheck[2]->Fill(samePi - 0.5);
        fHistSameTrackCheck[3]->Fill(sameAll - 0.5);
      }
      fClonesVector.clear();
      if (lIsGoodCandidate && lIsGoodPDG[3] && lIsGoodTrack[3] && lIsGoodCharge) fClonesVector.push_back(fCand);
    }
  }

  return kTRUE;
}

//______________________________________________________________________________
void AliSelectorFindableHyperTriton3Body::SlaveTerminate() {
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.
}

//______________________________________________________________________________
void AliSelectorFindableHyperTriton3Body::Terminate() {
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.
  TFile output(Form("%s/%s", fOutputFilePath.Data(), fOutputFileName.Data()), "RECREATE");
  GetOutputList()->Write();
  output.Close();
}

//______________________________________________________________________________
bool AliSelectorFindableHyperTriton3Body::IsTrackCutsGood(AliESDtrack *lTrack) {
  bool isCandidateGood = true;
  if (((lTrack->GetStatus() & AliVTrack::kTPCrefit) == 0) || (lTrack->GetTPCNcls() < 70) ||
      (lTrack->GetTPCchi2() > 4 * lTrack->GetTPCNcls()) || (std::fabs(lTrack->Eta()) > 0.9) ||
      (lTrack->GetKinkIndex(0) > 0)) {
    isCandidateGood = false;
  }
  return isCandidateGood;
}

//______________________________________________________________________________
bool AliSelectorFindableHyperTriton3Body::IsSameTrack(AliESDtrack *track1, AliESDtrack *track2) {
  return ((std::fabs(track1->GetTPCchi2() - track2->GetTPCchi2()) < 1.e-10) &&
          (track1->GetTPCNcls() == track2->GetTPCNcls()) && (track1->GetITSNcls() == track2->GetITSNcls()) &&
          (track1->GetStatus() == track2->GetStatus()));
}

//______________________________________________________________________________
HyperTritonCandidate::HyperTritonCandidate(AliESDtrack *track1, AliESDtrack *track2, AliESDtrack *track3)
    : fDaughterTrack{*track1, *track2, *track3}, fPDG(), fMId(), fTrueP(), fTrueDecayVtx(), fPrimaryVtx(), fMagField() {
}
