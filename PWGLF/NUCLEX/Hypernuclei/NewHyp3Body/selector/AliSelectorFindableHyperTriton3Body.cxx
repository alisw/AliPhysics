#define AliSelectorFindableHyperTriton3Body_cxx

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

#include "AliAnalysisTaskHypertriton3New.h"
#include "AliSelectorFindableHyperTriton3Body.h"

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

  const char lAM[3]{"AM"};
  const char *lSpecies[3]{"d", "p", "pi"};
  const char *lSpeciesPair[3]{"dp", "ppi", "pid"};
  const char *lCuts[4]{"ref", "std", "cuts", "hard"};
  const char *lStat[3]{"", "_fake", "_reflection'"};
  const char lCoords[4]{"xyz"};

  /// invariant mass distributions
  for (int iMatter = 0; iMatter < 2; iMatter++) {
    for (int iCuts = 0; iCuts < 4; iCuts++) {
      for (int iStat = 0; iStat < 3; iStat++) {
        fHistInvMass[iMatter][iCuts][iStat] =
            new TH2D(Form("fHistInvMass_%c_%s%s", lAM[iMatter], lCuts[iCuts], lStat[iStat]),
                     ";m (dp#pi) [GeV/#it{c}^{2}];Counts", 200, 2.95, 3.35, 20, 0, 10);
        GetOutputList()->Add(fHistInvMass[iMatter][iCuts][iStat]);
      }
    }
  }

  /// transverse momentum distributions
  // for (int iMatter = 0; iMatter < 2; iMatter++) {
  //   for (int iCuts = 0; iCuts < 4; iCuts++) {
  //     for (int iFake = 0; iFake < 2; iFake++) {
  //       fHistPt[iMatter][iCuts][iFake] = new TH1D(Form("fHistPt_%c_%s%s", lAM[iMatter], lCuts[iCuts], lFake[iFake]),
  //                                                 ";#it{p}_{T} [GeV/#it{c}];Counts", 100, 0, 10);
  //       GetOutputList()->Add(fHistPt[iMatter][iCuts][iFake]);
  //     }
  //   }
  // }

  /// daughters transverse momentum distributions
  // for (int iSpecies = 0; iSpecies < 3; iSpecies++) {
  //   for (int iCuts = 0; iCuts < 4; iCuts++) {
  //     for (int iFake = 0; iFake < 2; iFake++) {
  //       fHistDaughterPt[iSpecies][iCuts][iFake] =
  //           new TH1D(Form("fHistPt_%s_%s%s", lSpecies[iSpecies], lCuts[iCuts], lFake[iFake]),
  //                    ";#it{p}_{T} [GeV/#it{c}];Counts", 100, 0, 10);
  //       GetOutputList()->Add(fHistDaughterPt[iSpecies][iCuts][iFake]);
  //     }
  //   }
  // }

  /// daughter N Cluster TPC distribution
  // for (int iSpecies = 0; iSpecies < 3; iSpecies++) {
  //   for (int iCuts = 0; iCuts < 4; iCuts++) {
  //     for (int iFake = 0; iFake < 2; iFake++) {
  //       fHistNClsTPC[iSpecies][iCuts][iFake] =
  //           new TH1D(Form("fHistNClsTPC_%s_%s%s", lSpecies[iSpecies], lCuts[iCuts], lFake[iFake]), ":N_{cluster
  //           TPC}",
  //                    100, 0, 200);
  //       GetOutputList()->Add(fHistNClsTPC[iSpecies][iCuts][iFake]);
  //     }
  //   }
  // }

  /// track check histograms
  // fHistTrackCheck[0] = new TH1D("fHistAllTrackCheck", "", 2, -1, 1);
  // GetOutputList()->Add(fHistTrackCheck[0]);
  // fHistTrackCheck[1] = new TH1D("fHistDeuTrackCheck", "", 2, -1, 1);
  // GetOutputList()->Add(fHistTrackCheck[1]);
  // fHistTrackCheck[2] = new TH1D("fHistPTrackCheck", "", 2, -1, 1);
  // GetOutputList()->Add(fHistTrackCheck[2]);
  // fHistTrackCheck[3] = new TH1D("fHistPiTrackCheck", "", 2, -1, 1);
  // GetOutputList()->Add(fHistTrackCheck[3]);

  /// PID check histograms
  // for (int iSpecies = 0; iSpecies < 3; iSpecies++) {
  //   for (int iFake = 0; iFake < 2; iFake++) {
  //     fHistNSigma[iSpecies][iFake] =
  //         new TH1D(Form("fHistNSigma_%s%s", lSpecies[iSpecies], lFake[iFake]), "n#sigma", 10, 0, 10);
  //     fHistCheckPID[iSpecies][iFake] =
  //         new TH1D(Form("fHistCheckPID_%s%s", lSpecies[iSpecies], lFake[iFake]), "", 2, -1, 1);
  //     GetOutputList()->Add(fHistNSigma[iSpecies][iFake]);
  //     GetOutputList()->Add(fHistCheckPID[iSpecies][iFake]);
  //   }
  // }
  // fHistCheckPID[3][0] = new TH1D("fHistNSigmaCheckAll", "", 2, -1, 1);
  // fHistCheckPID[3][1] = new TH1D("fHistNSigmaCheckAll_Fake", "", 2, -1, 1);
  // GetOutputList()->Add(fHistCheckPID[3][0]);
  // GetOutputList()->Add(fHistCheckPID[3][1]);

  // PDG of the reco tracks histogram
  // fHistPDG = new TH2D("fHistDaughterPDG", "", 8, 0, 8, 3, 0, 3);
  // GetOutputList()->Add(fHistPDG);

  /// check on the reco tracks charge
  fHistChargeCheck = new TH1D("fHistChargeCheck", "", 2, -1, 1);
  GetOutputList()->Add(fHistChargeCheck);

  // pass selections histogram
  // fHistPassCheck[0] = new TH1D("fHistCheckPass", "", 6, 0, 6);
  // fHistPassCheck[0]->GetXaxis()->SetBinLabel(1, "no selection");
  // fHistPassCheck[0]->GetXaxis()->SetBinLabel(2, "track and charge");
  // fHistPassCheck[0]->GetXaxis()->SetBinLabel(3, "PID");
  // fHistPassCheck[0]->GetXaxis()->SetBinLabel(4, "good decay vertex");
  // fHistPassCheck[0]->GetXaxis()->SetBinLabel(5, "#it{cos#theta_{pointing}} selection");
  // fHistPassCheck[0]->GetXaxis()->LabelsOption("v");
  // GetOutputList()->Add(fHistPassCheck[0]);

  // fHistPassCheck[1] = new TH1D("fHistCheckFakePass", "", 6, 0, 6);
  // fHistPassCheck[1]->GetXaxis()->SetBinLabel(1, "no selection");
  // fHistPassCheck[1]->GetXaxis()->SetBinLabel(2, "track and charge");
  // fHistPassCheck[1]->GetXaxis()->SetBinLabel(3, "PID");
  // fHistPassCheck[1]->GetXaxis()->SetBinLabel(4, "good decay vertex");
  // fHistPassCheck[1]->GetXaxis()->SetBinLabel(5, "#it{cos#theta_{pointing}} selection");
  // fHistPassCheck[1]->GetXaxis()->LabelsOption("v");
  // GetOutputList()->Add(fHistPassCheck[1]);

  /// Chi2 vertex and cos(theta_pointing) distributions
  // for (int iFake = 0; iFake < 2; iFake++) {
  //   fHistVertexChi2[iFake] = new TH1D(Form("fHistVertexChi2%s", lFake[iFake]), "", 100, 0, 200);
  //   fHistCosPAngle[iFake] =
  //       new TH1D(Form("fCosPointingAngle%s", lFake[iFake]), ";#it{cos#theta_{pointing}} ;Counts", 5000, 0.5, 1.);
  //   GetOutputList()->Add(fHistVertexChi2[iFake]);
  //   GetOutputList()->Add(fHistCosPAngle[iFake]);
  // }

  /// decay vertex resolution histo
  // for (int iCoord = 0; iCoord < 3; iCoord++) {
  //   fHistResDecayVtx[iCoord] =
  //       new TH2D(Form("fHistResDecayVtx%c", lCoords[iCoord]),
  //                Form(";#Delta%c [mm];#it{p}_{T} [GeV/#it{c}]", lCoords[iCoord]), 500, -10, 10, 10, 0, 10);
  //   GetOutputList()->Add(fHistResDecayVtx[iCoord]);
  // }

  /// DCA to primary vertex
  // for (int iSpecies = 0; iSpecies < 3; iSpecies++) {
  //   for (int iFake = 0; iFake < 2; iFake++) {
  //     fHistDCA2pV[iSpecies][iFake] = new TH2D(Form("fDCA2Primaryvtx_%s%s", lSpecies[iSpecies], lFake[iFake]),
  //                                             ";mm ;#it{p}_{T} [GeV/#it{c}]", 600, 0, 30, 10, 0, 10);
  //     GetOutputList()->Add(fHistDCA2pV[iSpecies][iFake]);
  //   }
  // }

  /// DCA to secondary vertex
  // for (int iSpecies = 0; iSpecies < 3; iSpecies++) {
  //   for (int iFake = 0; iFake < 2; iFake++) {
  //     fHistDCA2dV[iSpecies][iFake] = new TH2D(Form("fDCA2Decayvtx_%s%s", lSpecies[iSpecies], lFake[iFake]),
  //                                             ";mm ;#it{p}_{T} [GeV/#it{c}]", 600, 0, 30, 10, 0, 10);
  //     GetOutputList()->Add(fHistDCA2dV[iSpecies][iFake]);
  //   }
  // }

  /// track distance to decay vertex histo
  // for (int iCoord = 0; iCoord < 3; iCoord++) {
  //   for (int iFake = 0; iFake < 2; iFake++) {
  //     fHistTrackDistance[iCoord][iFake] =
  //         new TH1D(Form("fTrackDistance%s_%s-%s", lSpeciesPair[iCoord], lSpeciesPair[(iCoord + 1) % 3],
  //         lFake[iFake]),
  //                  ";distance [mm]", 500, 0, 2000);
  //     GetOutputList()->Add(fHistTrackDistance[iCoord][iFake]);
  //   }
  // }

  /// decay length distribution
  // for (int iFake = 0; iFake < 2; iFake++) {
  //   fHistDecayLenght[iFake] = new TH1D(Form("fHistDecayLenght%s", lFake[iFake]), "", 100, 0, 200);
  //   GetOutputList()->Add(fHistDecayLenght[iFake]);
  // }
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
  bool lIsGoodTrack[4] = {false};
  bool lIsGoodCharge   = false;
  bool lIsGoodPID[4]   = {false};

  bool lIsGoodDCApi      = false;
  bool lIsGoodVertex     = false;
  bool lIsGoodPAngle     = false;
  bool lIsGoodVertexHard = false;
  bool lIsGoodPAngleHard = false;

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

  TLorentzVector fHypertriton = {0., 0., 0., 0.};
  TLorentzVector fLVDaughter[3];

  const float lMasses[3]{AliPID::ParticleMass(AliPID::kDeuteron), AliPID::ParticleMass(AliPID::kProton),
                         AliPID::ParticleMass(AliPID::kPion)};

  const int pdgCodes[7]{AliPID::ParticleCode(AliPID::kElectron), AliPID::ParticleCode(AliPID::kMuon),
                        AliPID::ParticleCode(AliPID::kPion),     AliPID::ParticleCode(AliPID::kKaon),
                        AliPID::ParticleCode(AliPID::kProton),   AliPID::ParticleCode(AliPID::kDeuteron),
                        AliPID::ParticleCode(AliPID::kHe3)};

  AliPID::EParticleType lSpecNum[3]{AliPID::kDeuteron, AliPID::kProton, AliPID::kPion};

  //------------------------------------------------------------
  /// Get stuff from tree
  //------------------------------------------------------------
  float lMagField      = *fTreeHyp3BodyVarMagneticField;
  unsigned char status = *fTreeHyp3BodyVarCandStat;

  bool lIsGood       = (status & g);
  bool lIsReflection = (status & r);
  bool lIsFake       = (!lIsGood && !lIsReflection);

  AliExternalTrackParam fParamTrack[3] = {*fTreeHyp3BodyVarTracks[0], *fTreeHyp3BodyVarTracks[1],
                                          *fTreeHyp3BodyVarTracks[2]};
  AliESDtrack t[3]{*fTreeHyp3BodyVarTracks[0], *fTreeHyp3BodyVarTracks[1], *fTreeHyp3BodyVarTracks[2]};

  // bool lIsHyp = (*fTreeHyp3BodyVarPDGcodes[0] == 1000010020) && (*fTreeHyp3BodyVarPDGcodes[1] == 2212) &&
  //               (*fTreeHyp3BodyVarPDGcodes[2] == -211);
  // bool lIsAntiHyp = (*fTreeHyp3BodyVarPDGcodes[0] == -1000010020) && (*fTreeHyp3BodyVarPDGcodes[1] == -2212) &&
  //                   (*fTreeHyp3BodyVarPDGcodes[2] == 211);
  // lIsGoodPDG = (lIsHyp || lIsAntiHyp);

  /// charge of the candidate
  bool lCharge = t[0].GetSign() < 0;

  /// findalble fake pdg check
  // for (int iD = 0; iD < 3; iD++) {
  //   for (int iP = 0; iP < 7; iP++) {
  //     if (std::abs(*fTreeHyp3BodyVarPDGcodes[iD]) == pdgCodes[iP]) fHistPDG->Fill(iP, iD);
  //   }
  // }

  /// track cuts check
  // int trk_sum = 0;
  // for (int iS = 0; iS < 3; iS++) {
  //   lIsGoodTrack[iS] = IsTrackCutsGood(&t[iS]);
  //   fHistTrackCheck[iS + 1]->Fill(lIsGoodTrack[iS] - 0.5);
  //   trk_sum += lIsGoodTrack[iS];
  // }
  // lIsGoodTrack[3] = (trk_sum == 3);
  // fHistTrackCheck[0]->Fill(lIsGoodTrack[3] - 0.5);

  /// PID selection on candidates
  // if (lIsGoodTrack[3]) {
  //   int pid_sum = 0;
  //   for (int iS = 0; iS < 3; iS++) {
  //     nSigma[iS]     = std::abs(fPIDResponse->NumberOfSigmasTPC(&t[iS], lSpecNum[iS]));
  //     lIsGoodPID[iS] = nSigma[iS] < 3.;
  //     fHistNSigma[iS][lIsFake]->Fill(nSigma[iS]);
  //     fHistCheckPID[iS][lIsFake]->Fill(lIsGoodPID[iS] - 0.5);
  //     pid_sum += lIsGoodPID[iS];
  //   }
  //   lIsGoodPID[3] = (pid_sum == 3);
  //   fHistCheckPID[3][lIsFake]->Fill(lIsGoodPID[3] - 0.5);
  // }

  /// charge check TODO: strano che non vengano già selezionati tutti dal task
  if (!(t[0].GetSign() != t[1].GetSign() || t[1].GetSign() == t[2].GetSign())) lIsGoodCharge = true;
  fHistChargeCheck->Fill(lIsGoodCharge - 0.5);

  // create the TLorentzVector of the hyper-triton candidate
  for (int iTrack = 0; iTrack < 3; iTrack++) {
    fLVDaughter[iTrack].SetXYZM(fParamTrack[iTrack].Px(), fParamTrack[iTrack].Py(), fParamTrack[iTrack].Pz(),
                                lMasses[iTrack]);
    fHypertriton += fLVDaughter[iTrack];
  }
  double hypMass = fHypertriton.M();
  double hypPt   = fHypertriton.Pt();

  // if (lIsGoodTrack[3] && lIsGoodCharge && lIsGoodPID[3]) {
  if (lIsGoodCharge) {
    // reconstruct the decay vertex with the dedicated vertexer
    bool recoVertex     = fHypertritonVertexer.FindDecayVertex(&t[0], &t[1], &t[2], *fTreeHyp3BodyVarMagneticField);
    bool recoVertexHard = fHypertritonVertexerHard.FindDecayVertex(&t[0], &t[1], &t[2], *fTreeHyp3BodyVarMagneticField);
    AliESDVertex *fDecayVertex     = dynamic_cast<AliESDVertex *>(fHypertritonVertexer.GetCurrentVertex());
    AliESDVertex *fDecayVertexHard = dynamic_cast<AliESDVertex *>(fHypertritonVertexerHard.GetCurrentVertex());
    if (recoVertex) {
      double vertexChi2NDF = fDecayVertex->GetChi2perNDF() * 1.;
      // fHistVertexChi2[lIsFake]->Fill(vertexChi2NDF);
    }
    /// possibility to constraint the goodness of the reconstructed decay vertex
    /// lIsGoodVertex = recoVertex && (fDecayVertex->GetChi2perNDF() < 20.);
    lIsGoodVertex     = recoVertex;
    lIsGoodVertexHard = recoVertexHard;

    /// compute decay vertex position and compute cos(pointingangle)
    if (lIsGoodVertex) {
      fDecayVertex->GetXYZ(fDecayVertexPos);
      for (int iTrack = 0; iTrack < 3; iTrack++) {
        // fHistResDecayVtx[iTrack]->Fill((fDecayVertexPos[iTrack] - *fTreeHyp3BodyVarDecayVtx[iTrack]) * 10., hypPt);
        fDecayLenght[iTrack] = fDecayVertexPos[iTrack] - *fTreeHyp3BodyVarPVtx[iTrack];
      }
      TVector3 v(fDecayLenght[0], fDecayLenght[1], fDecayLenght[2]);
      float pointAngle = fHypertriton.Angle(v);
      float cospa      = std::cos(pointAngle);
      // fHistCosPAngle[lIsFake]->Fill(cospa);
      lIsGoodPAngle = (cospa > 0.80);

      // fHistDecayLenght[lIsFake]->Fill(v.Mag());

      /// compute the DCA of the 3 tracks from the primary and decay vertex
      fPrimaryVertexPos[0] = *fTreeHyp3BodyVarPVtx[0] * 1.;
      fPrimaryVertexPos[1] = *fTreeHyp3BodyVarPVtx[1] * 1.;
      fPrimaryVertexPos[2] = *fTreeHyp3BodyVarPVtx[2] * 1.;
      AliESDVertex fPV(fPrimaryVertexPos, fPrimaryVertexCov, 1., 1000);
      for (int iTrack = 0; iTrack < 3; iTrack++) {
        fParamTrack[iTrack].PropagateToDCA(fDecayVertex, lMagField, 1000., fDCA2dv[iTrack], fDCA2dvcov[iTrack]);
        fParamTrack[iTrack].PropagateToDCA(&fPV, lMagField, 1000., fDCA2pv[iTrack], fDCA2pvcov[iTrack]);
        float ddv = std::sqrt(fDCA2dv[iTrack][0] * fDCA2dv[iTrack][0] + fDCA2dv[iTrack][1] * fDCA2dv[iTrack][1] +
                              fDCA2dv[iTrack][2] * fDCA2dv[iTrack][2]);
        float dpv = std::sqrt(fDCA2pv[iTrack][0] * fDCA2pv[iTrack][0] + fDCA2pv[iTrack][1] * fDCA2pv[iTrack][1] +
                              fDCA2pv[iTrack][2] * fDCA2pv[iTrack][2]);
        // fHistDCA2dV[iTrack][lIsFake]->Fill(ddv * 10., hypPt);
        // fHistDCA2pV[iTrack][lIsFake]->Fill(dpv * 10., hypPt);
        if (iTrack == 2) lIsGoodDCApi = dpv > 0.1;
      }

      /// compute the track2track distance used in the vertexer
      float pPM[3][3];
      for (int iTrack = 0; iTrack < 3; iTrack++) {
        fHypertritonVertexer.Find2ProngClosestPoint(&fParamTrack[iTrack], &fParamTrack[(iTrack + 1) % 3], lMagField,
                                                    pPM[iTrack]);
      }

      // for (int iPerm = 0; iPerm < 3; iPerm++) {
      //   float d = 0.;
      //   for (int iDim = 0; iDim < 3; ++iDim) {
      //     d += (pPM[iPerm][iDim] - pPM[(iPerm + 1) % 3][iDim]) * (pPM[iPerm][iDim] - pPM[(iPerm + 1) % 3][iDim]);
      //   }
      //   fHistTrackDistance[iPerm][lIsFake]->Fill(std::sqrt(d) * 10.);
      // }

      /// again with harder selections
      if (lIsGoodVertexHard) {
        fDecayVertexHard->GetXYZ(fDecayVertexPosHard);
        for (int iTrack = 0; iTrack < 3; iTrack++) {
          fDecayLenghtHard[iTrack] = fDecayVertexPosHard[iTrack] - *fTreeHyp3BodyVarPVtx[iTrack];
        }
        TVector3 v(fDecayLenghtHard[0], fDecayLenghtHard[1], fDecayLenghtHard[2]);
        float pointAngle  = fHypertriton.Angle(v);
        float cospa       = std::cos(pointAngle);
        lIsGoodPAngleHard = (cospa > 0.997);
      }
    }
  }

  /// invariant mass distribution for efficiency
  if (lIsGood) fHistInvMass[lCharge][1][0]->Fill(hypMass, hypPt);
  if (lIsFake) fHistInvMass[lCharge][1][1]->Fill(hypMass, hypPt);
  if (lIsReflection) fHistInvMass[lCharge][1][2]->Fill(hypMass, hypPt);
  /// track and PID selections
  // if (lIsGoodTrack[3] && lIsGoodCharge && lIsGoodPID[3]) fHistInvMass[lCharge][1][lIsFake]->Fill(hypMass, hypPt);
  /// vertex and cos(PA) selections
  if (lIsGoodDCApi && lIsGoodVertex && lIsGoodPAngle) {
    if (lIsGood) fHistInvMass[lCharge][2][0]->Fill(hypMass, hypPt);
    if (lIsFake) fHistInvMass[lCharge][2][1]->Fill(hypMass, hypPt);
    if (lIsReflection) fHistInvMass[lCharge][2][2]->Fill(hypMass, hypPt);
  }
  /// hard vertex and cos(PA) selections
  // if (lIsGoodTrack[3] && lIsGoodCharge && lIsGoodPID[3] && lIsGoodDCApi && lIsGoodVertex && lIsGoodPAngleHard)
  //   fHistInvMass[lCharge][3][lIsFake]->Fill(hypMass, hypPt);

  // transverse momentum distribution
  // fHistPt[lCharge][0][lIsFake]->Fill(hypPt);
  // if (lIsGoodTrack[3] && lIsGoodCharge && lIsGoodPID[3]) fHistPt[lCharge][1][lIsFake]->Fill(hypPt);
  // if (lIsGoodTrack[3] && lIsGoodCharge && lIsGoodPID[3] && lIsGoodDCApi && lIsGoodVertex && lIsGoodPAngle)
  //   fHistPt[lCharge][2][lIsFake]->Fill(hypPt);
  // if (lIsGoodTrack[3] && lIsGoodCharge && lIsGoodPID[3] && lIsGoodDCApi && lIsGoodVertex && lIsGoodPAngleHard)
  //   fHistPt[lCharge][3][lIsFake]->Fill(hypPt);

  /// daughter transverse momentum distribution
  // for (int iSpecies = 0; iSpecies < 3; iSpecies++) {
  //   double pt = t[iSpecies].Pt() * 1.;
  //   fHistDaughterPt[iSpecies][0][lIsFake]->Fill(pt);
  //   if (lIsGoodTrack[3] && lIsGoodCharge && lIsGoodPID[3]) fHistDaughterPt[iSpecies][1][lIsFake]->Fill(pt);
  //   if (lIsGoodTrack[3] && lIsGoodCharge && lIsGoodPID[3] && lIsGoodDCApi && lIsGoodVertex && lIsGoodPAngle)
  //     fHistDaughterPt[iSpecies][2][lIsFake]->Fill(pt);
  //   if (lIsGoodTrack[3] && lIsGoodCharge && lIsGoodPID[3] && lIsGoodDCApi && lIsGoodVertex && lIsGoodPAngleHard)
  //     fHistDaughterPt[iSpecies][3][lIsFake]->Fill(pt);
  // }

  /// daughter N Cluster TPC distribution
  // for (int iSpecies = 0; iSpecies < 3; iSpecies++) {
  //   int ncls = (int)t[iSpecies].GetTPCNcls();
  //   fHistNClsTPC[iSpecies][0][lIsFake]->Fill(ncls);
  //   if (lIsGoodTrack[3] && lIsGoodCharge && lIsGoodPID[3]) fHistNClsTPC[iSpecies][1][lIsFake]->Fill(ncls);
  //   if (lIsGoodTrack[3] && lIsGoodCharge && lIsGoodPID[3] && lIsGoodDCApi && lIsGoodVertex && lIsGoodPAngle)
  //     fHistNClsTPC[iSpecies][2][lIsFake]->Fill(ncls);
  //   if (lIsGoodTrack[3] && lIsGoodCharge && lIsGoodPID[3] && lIsGoodDCApi && lIsGoodVertex && lIsGoodPAngleHard)
  //     fHistNClsTPC[iSpecies][3][lIsFake]->Fill(ncls);
  // }

  // total number of findable candidates and clones
  // fHistPassCheck[lIsFake]->Fill(0.5);
  // if (lIsGoodTrack[3] && lIsGoodCharge) fHistPassCheck[lIsFake]->Fill(1.5);
  // if (lIsGoodTrack[3] && lIsGoodCharge && lIsGoodPID[3]) fHistPassCheck[lIsFake]->Fill(2.5);
  // if (lIsGoodTrack[3] && lIsGoodCharge && lIsGoodPID[3] && lIsGoodDCApi && lIsGoodVertex)
  //   fHistPassCheck[lIsFake]->Fill(3.5);
  // if (lIsGoodTrack[3] && lIsGoodCharge && lIsGoodPID[3] && lIsGoodDCApi && lIsGoodVertex && lIsGoodPAngle)
  //   fHistPassCheck[lIsFake]->Fill(4.5);
  // if (lIsGoodTrack[3] && lIsGoodCharge && lIsGoodPID[3] && lIsGoodDCApi && lIsGoodVertex && lIsGoodPAngleHard)
  //   fHistPassCheck[lIsFake]->Fill(5.5);

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