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

#include "AliSelectorFindableHyperTriton3Body.h"

const unsigned char g = 0x1; // on if is a good candidate
const unsigned char r = 0x2; // on if is reflection candidate

/// usefull functions
template <typename T> double Pot2(T a) { return a * a; }

template <typename T> double Distance(T pX, T pY, T pZ, T dX, T dY, T dZ) {
  return std::sqrt(Pot2(pX - dX) + Pot2(pY - dY) + Pot2(pZ - dZ));
}

AliSelectorFindableHyperTriton3Body::AliSelectorFindableHyperTriton3Body(TString outputName, TString outputPath,
                                                                         TTree *)
    : fOutputFileName{outputName}, //
      fOutputFilePath{outputPath}, //
      fHypertritonVertexer(),      //
      fHypertritonVertexerHard() {
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

  const char lAM[3]{"AM"};
  const char *lSpecies[3]{"d", "p", "pi"};
  const char *lSpeciesPair[3]{"dp", "ppi", "pid"};
  const char *lCuts[3]{"std", "cuts", "hard"};
  const char *lStat[3]{"", "_reflection", "_fake"};
  const char *lFake[2]{"", "_fake"};
  const char lCoords[4]{"xyz"};

  /// invariant mass distributions
  for (int iMatter = 0; iMatter < 2; iMatter++) {
    for (int iCuts = 0; iCuts < 3; iCuts++) {
      for (int iStat = 0; iStat < 3; iStat++) {
        fHistInvMass[iMatter][iCuts][iStat] =
            new TH2D(Form("fHistInvMass_%c_%s%s", lAM[iMatter], lCuts[iCuts], lStat[iStat]),
                     ";m (dp#pi) [GeV/#it{c}^{2}];Counts", 200, 2.95, 3.35, 20, 0, 10);
        GetOutputList()->Add(fHistInvMass[iMatter][iCuts][iStat]);
      }
    }
  }

  /// transverse momentum distributions
  for (int iMatter = 0; iMatter < 2; iMatter++) {
    for (int iCuts = 0; iCuts < 3; iCuts++) {
      for (int iStat = 0; iStat < 3; iStat++) {
        fHistPt[iMatter][iCuts][iStat] = new TH1D(Form("fHistPt_%c_%s%s", lAM[iMatter], lCuts[iCuts], lStat[iStat]),
                                                  ";#it{p}_{T} [GeV/#it{c}];Counts", 100, 0, 10);
        GetOutputList()->Add(fHistPt[iMatter][iCuts][iStat]);
      }
    }
  }

  /// daughters transverse momentum distributions
  for (int iSpecies = 0; iSpecies < 3; iSpecies++) {
    for (int iCuts = 0; iCuts < 3; iCuts++) {
      for (int iStat = 0; iStat < 3; iStat++) {
        fHistDaughterPt[iSpecies][iCuts][iStat] =
            new TH1D(Form("fHistPt_%s_%s%s", lSpecies[iSpecies], lCuts[iCuts], lStat[iStat]),
                     ";#it{p}_{T} [GeV/#it{c}];Counts", 100, 0, 10);
        GetOutputList()->Add(fHistDaughterPt[iSpecies][iCuts][iStat]);
      }
    }
  }

  /// daughter N Cluster ITS distribution
  for (int iSpecies = 0; iSpecies < 3; iSpecies++) {
    for (int iFake = 0; iFake < 2; iFake++) {
      fHistNclsITS[iSpecies][iFake] =
          new TH1D(Form("fHistNclsITS_%s%s", lSpecies[iSpecies], lFake[iFake]), ":N_{cluster ITS}", 9, -1, 8);
      GetOutputList()->Add(fHistNclsITS[iSpecies][iFake]);
    }
  }

  /// daughter N Cluster TPC distribution
  for (int iSpecies = 0; iSpecies < 3; iSpecies++) {
    for (int iFake = 0; iFake < 2; iFake++) {
      fHistNclsTPC[iSpecies][iFake] =
          new TH1D(Form("fHistNclsTPC_%s%s", lSpecies[iSpecies], lFake[iFake]), ":N_{cluster TPC}", 100, 0, 200);
      GetOutputList()->Add(fHistNclsTPC[iSpecies][iFake]);
    }
  }

  /// tracks global Chi2 distributions
  for (int iSpecies = 0; iSpecies < 3; iSpecies++) {
    for (int iFake = 0; iFake < 2; iFake++) {
      fHistGlobalTrackChi2[iSpecies][iFake] = new TH1D(
          Form("fHistGlobalTrackChi2_%s%s", lSpecies[iSpecies], lFake[iFake]), ":N_{cluster TPC}", 100, 0, 200);
      GetOutputList()->Add(fHistGlobalTrackChi2[iSpecies][iFake]);
    }
  }

  /// PID check histograms
  for (int iSpecies = 0; iSpecies < 3; iSpecies++) {
    for (int iStat = 0; iStat < 2; iStat++) {
      fHistNSigmaTPC[iSpecies][iStat] =
          new TH1D(Form("fHistNSigmaTPC_%s%s", lSpecies[iSpecies], lStat[iStat]), "n#sigma", 10, 0, 10);
      fHistNSigmaTOF[iSpecies][iStat] =
          new TH1D(Form("fHistNSigmaTOF_%s%s", lSpecies[iSpecies], lStat[iStat]), "n#sigma", 10, 0, 10);
      fHistCheckPID[iSpecies][iStat] =
          new TH1D(Form("fHistCheckPID_%s%s", lSpecies[iSpecies], lStat[iStat]), "", 2, -1, 1);
      GetOutputList()->Add(fHistNSigmaTPC[iSpecies][iStat]);
      GetOutputList()->Add(fHistNSigmaTOF[iSpecies][iStat]);
      GetOutputList()->Add(fHistCheckPID[iSpecies][iStat]);
    }
  }
  fHistCheckPID[3][0] = new TH1D("fHistNSigmaCheckAll", "", 2, -1, 1);
  fHistCheckPID[3][1] = new TH1D("fHistNSigmaCheckAll_Fake", "", 2, -1, 1);
  GetOutputList()->Add(fHistCheckPID[3][0]);
  GetOutputList()->Add(fHistCheckPID[3][1]);

  /// Chi2 vertex and cos(theta_pointing) distributions
  for (int iFake = 0; iFake < 2; iFake++) {
    fHistVertexChi2[iFake] = new TH1D(Form("fHistVertexChi2%s", lFake[iFake]), "", 100, 0, 200);
    fHistCosPAngle[iFake] =
        new TH1D(Form("fCosPointingAngle%s", lFake[iFake]), ";#it{cos#theta_{pointing}} ;Counts", 5000, 0.5, 1.);
    GetOutputList()->Add(fHistVertexChi2[iFake]);
    GetOutputList()->Add(fHistCosPAngle[iFake]);
  }

  /// decay vertex resolution histo
  for (int iCoord = 0; iCoord < 3; iCoord++) {
    fHistResDecayVtx[iCoord] =
        new TH1D(Form("fHistResDecayVtx%c", lCoords[iCoord]), Form(";#Delta%c [mm]", lCoords[iCoord]), 500, -10, 10);
    GetOutputList()->Add(fHistResDecayVtx[iCoord]);
  }

  /// DCA to primary vertex
  for (int iSpecies = 0; iSpecies < 3; iSpecies++) {
    for (int iFake = 0; iFake < 2; iFake++) {
      fHistDCA2pV[iSpecies][iFake] =
          new TH1D(Form("fDCA2Primaryvtx_%s%s", lSpecies[iSpecies], lFake[iFake]), ";mm ", 600, 0, 30);
      GetOutputList()->Add(fHistDCA2pV[iSpecies][iFake]);
    }
  }

  /// DCA to secondary vertex
  for (int iSpecies = 0; iSpecies < 3; iSpecies++) {
    for (int iFake = 0; iFake < 2; iFake++) {
      fHistDCA2dV[iSpecies][iFake] =
          new TH1D(Form("fDCA2Decayvtx_%s%s", lSpecies[iSpecies], lFake[iFake]), ";mm ", 600, 0, 30);
      GetOutputList()->Add(fHistDCA2dV[iSpecies][iFake]);
    }
  }

  /// track distance to decay vertex histo
  for (int iCoord = 0; iCoord < 3; iCoord++) {
    for (int iFake = 0; iFake < 2; iFake++) {
      fHistTrackDistance[iCoord][iFake] =
          new TH1D(Form("lTrackDistance%s_%s-%s", lSpeciesPair[iCoord], lSpeciesPair[(iCoord + 1) % 3], lFake[iFake]),
                   ";distance [mm]", 500, 0, 2000);
      GetOutputList()->Add(fHistTrackDistance[iCoord][iFake]);
    }
  }

  /// ct distribution
  for (int iFake = 0; iFake < 2; iFake++) {
    fHistCT[iFake] = new TH1D(Form("fHistCT%s", lFake[iFake]), "", 50, 0, 100);
    GetOutputList()->Add(fHistCT[iFake]);
  }
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
  bool lIsGoodPID[4] = {false};

  bool lIsGoodDCApi      = false;
  bool lIsGoodVertex     = false;
  bool lIsGoodPAngle     = false;
  bool lIsGoodVertexHard = false;
  bool lIsGoodPAngleHard = false;

  double lPrimaryVertexCov[6]   = {0.};
  double lDecayVertexPos[3]     = {0., 0., 0.};
  double lDecayVertexPosHard[3] = {0., 0., 0.};
  double lDecayLenght[3]        = {0., 0., 0.};
  double lDecayLenghtHard[3]    = {0., 0., 0.};
  double lDCA2dv[3][3]          = {{0.}};
  double lDCA2pv[3][3]          = {{0.}};
  double lDCA2dvcov[3][6]       = {{0.}};
  double lDCA2pvcov[3][6]       = {{0.}};

  TLorentzVector lHypertriton = {0., 0., 0., 0.};
  TLorentzVector lLVDaughter[3];

  const float lMasses[3]{AliPID::ParticleMass(AliPID::kDeuteron), AliPID::ParticleMass(AliPID::kProton),
                         AliPID::ParticleMass(AliPID::kPion)};

  //------------------------------------------------------------
  // Get stuff from tree
  //------------------------------------------------------------
  float lMagField = *fTreeHyp3BodyVarMagneticField;
  unsigned char s = *fTreeHyp3BodyVarCandStat;

  bool lIsGood       = (s & g);
  bool lIsReflection = (s & r);
  // bool lIsFake       = (!lIsGood && !lIsReflection);

  int cStatus = 0;
  if (lIsReflection) cStatus = 1;
  if (!lIsGood && !lIsReflection) cStatus = 2;

  AliExternalTrackParam *lTrackDeu = new AliExternalTrackParam(*fTreeHyp3BodyVarTracks[0]);
  AliExternalTrackParam *lTrackP   = new AliExternalTrackParam(*fTreeHyp3BodyVarTracks[1]);
  AliExternalTrackParam *lTrackPi  = new AliExternalTrackParam(*fTreeHyp3BodyVarTracks[2]);

  AliExternalTrackParam *lTrack[3]{lTrackDeu, lTrackP, lTrackPi};

  int lNclsITS[3]{*fTreeHyp3BodyVarNclsITS[0], *fTreeHyp3BodyVarNclsITS[1], *fTreeHyp3BodyVarNclsITS[2]};
  int lNclsTPC[3]{*fTreeHyp3BodyVarNclsTPC[0], *fTreeHyp3BodyVarNclsTPC[1], *fTreeHyp3BodyVarNclsTPC[2]};

  double lTrackGlobalChi2[3]{*fTreeHyp3BodyVarGlobalChi2[0] * 1., *fTreeHyp3BodyVarGlobalChi2[1] * 1.,
                             *fTreeHyp3BodyVarGlobalChi2[2] * 1.};

  double lNsigmaTPC[3]{*fTreeHyp3BodyVarNsigmaTPC[0] * 1., *fTreeHyp3BodyVarNsigmaTPC[1] * 1.,
                       *fTreeHyp3BodyVarNsigmaTPC[2] * 1.};
  double lNsigmaTOF[3]{*fTreeHyp3BodyVarNsigmaTOF[0] * 1., *fTreeHyp3BodyVarNsigmaTOF[1] * 1.,
                       *fTreeHyp3BodyVarNsigmaTOF[2] * 1.};

  // unsigned long lTrackFlag[3]{*fTreeHyp3BodyVarFlags[0], *fTreeHyp3BodyVarFlags[1], *fTreeHyp3BodyVarFlags[2]};

  double lTrueDecayVtx[3]{*fTreeHyp3BodyVarDecayVtx[0] * 1., *fTreeHyp3BodyVarDecayVtx[1] * 1.,
                          *fTreeHyp3BodyVarDecayVtx[2] * 1.};
  double lTruePrimaryVtx[3]{*fTreeHyp3BodyVarPVtx[0] * 1., *fTreeHyp3BodyVarPVtx[1] * 1.,
                            *fTreeHyp3BodyVarPVtx[2] * 1.};

  /// charge of the candidate
  bool lCharge = lTrack[0]->GetSign() < 0;

  //------------------------------------------------------------
  // N cluster in ITS, TPC and Track Chi2 distributions
  //------------------------------------------------------------
  for (int iTrack = 0; iTrack < 3; iTrack++) {
    fHistNclsITS[iTrack][lIsGood]->Fill(lNclsITS[iTrack]);
    fHistNclsTPC[iTrack][lIsGood]->Fill(lNclsTPC[iTrack]);
    fHistGlobalTrackChi2[iTrack][lIsGood]->Fill(lTrackGlobalChi2[iTrack]);
  }

  //------------------------------------------------------------
  // PID selection of candidates
  //------------------------------------------------------------
  int pid_sum = 0;
  for (int iS = 0; iS < 3; iS++) {
    fHistNSigmaTPC[iS][lIsGood]->Fill(lNsigmaTPC[iS]);
    fHistNSigmaTOF[iS][lIsGood]->Fill(lNsigmaTOF[iS]);
    lIsGoodPID[iS] = lNsigmaTPC[iS] < 3.;
    fHistCheckPID[iS][lIsGood]->Fill(lIsGoodPID[iS] - 0.5);
    pid_sum += lIsGoodPID[iS];
  }
  lIsGoodPID[3] = (pid_sum == 3);
  fHistCheckPID[3][lIsGood]->Fill(lIsGoodPID[3] - 0.5);

  //------------------------------------------------------------
  // Secondary vertex reconstruction
  //------------------------------------------------------------

  /// create the TLorentzVector of the hyper-triton candidate
  for (int iTrack = 0; iTrack < 3; iTrack++) {
    lLVDaughter[iTrack].SetXYZM(lTrack[iTrack]->Px(), lTrack[iTrack]->Py(), lTrack[iTrack]->Pz(), lMasses[iTrack]);
    lHypertriton += lLVDaughter[iTrack];
  }
  double hypMass = lHypertriton.M();
  double hypPt   = lHypertriton.Pt();
  // double hypP    = lHypertriton.P();

  if (lIsGoodPID[3]) {
    // reconstruct the decay vertex with the dedicated vertexer
    bool recoVertex            = fHypertritonVertexer.FindDecayVertex(lTrack[0], lTrack[1], lTrack[2], lMagField);
    bool recoVertexHard        = fHypertritonVertexerHard.FindDecayVertex(lTrack[0], lTrack[1], lTrack[2], lMagField);
    AliESDVertex *lDecayVertex = static_cast<AliESDVertex *>(fHypertritonVertexer.GetCurrentVertex());
    AliESDVertex *lDecayVertexHard = static_cast<AliESDVertex *>(fHypertritonVertexerHard.GetCurrentVertex());
    if (recoVertex) {
      double vertexChi2NDF = lDecayVertex->GetChi2perNDF() * 1.;
      fHistVertexChi2[lIsGood]->Fill(vertexChi2NDF);
    }
    /// possibility to constraint the goodness of the reconstructed decay vertex
    /// lIsGoodVertex = recoVertex && (fDecayVertex->GetChi2perNDF() < 20.);
    lIsGoodVertex     = recoVertex;
    lIsGoodVertexHard = recoVertexHard;

    /// compute decay vertex position and compute cos(pointingangle)
    if (lIsGoodVertex) {
      lDecayVertex->GetXYZ(lDecayVertexPos);
      for (int iTrack = 0; iTrack < 3; iTrack++) {
        fHistResDecayVtx[iTrack]->Fill((lDecayVertexPos[iTrack] - lTrueDecayVtx[iTrack]) * 10.);
        lDecayLenght[iTrack] = lDecayVertexPos[iTrack] - lTruePrimaryVtx[iTrack];
      }

      TVector3 v(lDecayLenght[0], lDecayLenght[1], lDecayLenght[2]);
      fHistCT[lIsGood]->Fill(v.Mag());

      float pointAngle = lHypertriton.Angle(v);
      float cospa      = std::cos(pointAngle);

      fHistCosPAngle[lIsGood]->Fill(cospa);
      lIsGoodPAngle = (cospa > 0.90);

      /// compute the DCA of the 3 tracks from the primary and decay vertex
      AliESDVertex lPV(lTruePrimaryVtx, lPrimaryVertexCov, 1., 1000);

      for (int iTrack = 0; iTrack < 3; iTrack++) {
        lTrack[iTrack]->PropagateToDCA(lDecayVertex, lMagField, 1000., lDCA2dv[iTrack], lDCA2dvcov[iTrack]);
        lTrack[iTrack]->PropagateToDCA(&lPV, lMagField, 1000., lDCA2pv[iTrack], lDCA2pvcov[iTrack]);

        float ddv = std::sqrt(Pot2(lDCA2dv[iTrack][0]) + Pot2(lDCA2dv[iTrack][1]) + Pot2(lDCA2dv[iTrack][2]));
        float dpv = std::sqrt(Pot2(lDCA2pv[iTrack][0]) + Pot2(lDCA2pv[iTrack][1]) + Pot2(lDCA2pv[iTrack][2]));
        fHistDCA2dV[iTrack][lIsGood]->Fill(ddv * 10.);
        fHistDCA2pV[iTrack][lIsGood]->Fill(dpv * 10.);
        // TODO: provare a vedere cosa succede all'efficienza facendo questo taglio su tutte e 3 le specie
        if (iTrack == 2) lIsGoodDCApi = dpv > 0.1;
      }
    }

    /// compute the track2track distance used in the vertexer
    // float pPM[3][3];
    // for (int iTrack = 0; iTrack < 3; iTrack++) {
    //   fHypertritonVertexer.Find2ProngClosestPoint(&lTrack[iTrack], &lTrack[(iTrack + 1) % 3], lMagField,
    //   pPM[iTrack]);
    // }

    // for (int iPerm = 0; iPerm < 3; iPerm++) {
    //   float d = 0.;
    //   for (int iDim = 0; iDim < 3; ++iDim) {
    //     d += (pPM[iPerm][iDim] - pPM[(iPerm + 1) % 3][iDim]) * (pPM[iPerm][iDim] - pPM[(iPerm + 1) % 3][iDim]);
    //   }
    //   fHistTrackDistance[iPerm][lIsGood]->Fill(std::sqrt(d) * 10.);
    // }

    /// again with harder selections
    if (lIsGoodVertexHard) {
      lDecayVertexHard->GetXYZ(lDecayVertexPosHard);
      for (int iTrack = 0; iTrack < 3; iTrack++) {
        lDecayLenghtHard[iTrack] = lDecayVertexPosHard[iTrack] - lTruePrimaryVtx[iTrack];
      }
      TVector3 v(lDecayLenghtHard[0], lDecayLenghtHard[1], lDecayLenghtHard[2]);

      float pointAngle  = lHypertriton.Angle(v);
      float cospa       = std::cos(pointAngle);
      lIsGoodPAngleHard = (cospa > 0.998);
    }
  }

  /// invariant mass distribution for efficiency
  fHistInvMass[lCharge][0][cStatus]->Fill(hypMass, hypPt);
  if (lIsGoodDCApi && lIsGoodVertex && lIsGoodPAngle) fHistInvMass[lCharge][1][cStatus]->Fill(hypMass, hypPt);
  if (lIsGoodDCApi && lIsGoodVertexHard && lIsGoodPAngleHard) fHistInvMass[lCharge][2][cStatus]->Fill(hypMass, hypPt);

  /// transverse momentum distribution
  fHistPt[lCharge][0][cStatus]->Fill(hypPt);
  if (lIsGoodDCApi && lIsGoodVertex && lIsGoodPAngle) fHistPt[lCharge][1][cStatus]->Fill(hypPt);
  if (lIsGoodDCApi && lIsGoodVertexHard && lIsGoodPAngleHard) fHistPt[lCharge][2][cStatus]->Fill(hypPt);

  /// daughter transverse momentum distribution
  for (int iSpecies = 0; iSpecies < 3; iSpecies++) {
    double pt = lTrack[iSpecies]->Pt() * 1.;
    fHistDaughterPt[iSpecies][0][cStatus]->Fill(pt);
    if (lIsGoodDCApi && lIsGoodVertex && lIsGoodPAngle) fHistDaughterPt[iSpecies][1][cStatus]->Fill(pt);
    if (lIsGoodDCApi && lIsGoodVertexHard && lIsGoodPAngleHard) fHistDaughterPt[iSpecies][2][cStatus]->Fill(pt);
  }

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