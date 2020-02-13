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

/// usefull functions

namespace {


template <typename T> double Sq(T a) { return a * a; }

template <typename T> double Point2PointDistance(T *p0, T *p1) {
  double d2 = 0.;
  for (int iDim = 0; iDim < 3; ++iDim) {
    d2 += Sq(p0[iDim] - p1[iDim]);
  }
  return std::sqrt(d2);
}

template <typename T> double Norm(T x, T y) { return std::sqrt(Sq(x) + Sq(y)); }

template <typename T> double Norm(T x, T y, T z) { return std::sqrt(Sq(x) + Sq(y) + Sq(z)); }
} // namespace

AliSelectorFindableHyperTriton3Body::AliSelectorFindableHyperTriton3Body(TString outputName, TString outputPath,TTree *):
fOutputFileName{outputName}, //
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
  const char *lCuts[3]{"findable", "cuts", "hard"};
  const char lCoords[3]{'x','y','z'};

  /// invariant mass distributions
  for (int iMatter = 0; iMatter < 2; iMatter++) {
    for (int iCuts = 0; iCuts < 3; iCuts++) {
      fHistInvMass[iMatter][iCuts] = new TH2D(Form("fHistInvMass_%c_%s", lAM[iMatter], lCuts[iCuts]), ";#it{m} (dp#pi) (GeV/#it{c}^{2});#it{p}_{T} (GeV/#it{c})", 200, 2.95, 3.35, 20, 0, 10);
      GetOutputList()->Add(fHistInvMass[iMatter][iCuts]);
    }
  }

  /// transverse momentum distributions
  for (int iMatter = 0; iMatter < 2; iMatter++) {
    for (int iCuts = 0; iCuts < 3; iCuts++) {
      fHistPt[iMatter][iCuts] = new TH1D(Form("fHistPt_%c_%s", lAM[iMatter], lCuts[iCuts]), ";#it{p}_{T} (GeV/#it{c});Counts", 100, 0, 10);
      GetOutputList()->Add(fHistPt[iMatter][iCuts]);
    }
  }

  /// daughters transverse momentum distributions
  for (int iSpecies = 0; iSpecies < 3; iSpecies++) {
    for (int iCuts = 0; iCuts < 3; iCuts++) {
      fHistDaughterPt[iSpecies][iCuts] = new TH1D(Form("fHistPt_%s_%s", lSpecies[iSpecies], lCuts[iCuts]), ";#it{p}_{T} (GeV/#it{c});Counts", 100, 0, 10);
      GetOutputList()->Add(fHistDaughterPt[iSpecies][iCuts]);
    }
  }

  /// daughter N Cluster ITS distribution
  for (int iSpecies = 0; iSpecies < 3; iSpecies++) {
    fHistNclsITS[iSpecies] = new TH1D(Form("fHistNclsITS_%s", lSpecies[iSpecies]), ";N_{cluster ITS};Counts", 9, -1, 8);
    GetOutputList()->Add(fHistNclsITS[iSpecies]);
  }

  /// check on Ncls ITS for all the three tracks
  fHistNclsITSCheck = new TH1D("fHistNclsITSCheck", "", 2, -1, 1);
  GetOutputList()->Add(fHistNclsITSCheck);

  /// daughter N Cluster TPC distribution
  for (int iSpecies = 0; iSpecies < 3; iSpecies++) {
    fHistNclsTPC[iSpecies] = new TH1D(Form("fHistNclsTPC_%s", lSpecies[iSpecies]), ";N_{cluster TPC};Counts", 100, 0, 200);
    GetOutputList()->Add(fHistNclsTPC[iSpecies]);
  }

  /// PID check histograms
  for (int iSpecies = 0; iSpecies < 3; iSpecies++) {
    fHistNSigmaTPC[iSpecies] = new TH1D(Form("fHistNSigmaTPC_%s", lSpecies[iSpecies]), "n#sigma_{TPC}", 10, 0, 10);
    fHistNSigmaTOF[iSpecies] = new TH1D(Form("fHistNSigmaTOF_%s", lSpecies[iSpecies]), "n#sigma_{TOF}", 10, 0, 10);
    fHistCheckPID = new TH1D(Form("fHistCheckPID_%s", lSpecies[iSpecies]), "", 2, -1, 1);
    GetOutputList()->Add(fHistNSigmaTPC[iSpecies]);
    GetOutputList()->Add(fHistNSigmaTOF[iSpecies]);
    GetOutputList()->Add(fHistCheckPID);
  }
  fHistCheckPID = new TH1D("fHistNSigmaCheckAll", "", 2, -1, 1);
  GetOutputList()->Add(fHistCheckPID);

  /// Chi2 vertex and cos(theta_pointing) distributions
  fHistVertexChi2 = new TH1D("fHistVertexChi2", "", 100, 0, 200);
  fHistCosPAngle  = new TH1D("fCosPointingAngle", ";#it{cos#theta_{pointing}};Counts", 5000, 0.5, 1.);
  GetOutputList()->Add(fHistVertexChi2);
  GetOutputList()->Add(fHistCosPAngle);

  /// decay vertex resolution histo
  for (int iCoord = 0; iCoord < 3; iCoord++) {
    fHistResDecayVtx[iCoord] =
        new TH1D(Form("fHistResDecayVtx%c", lCoords[iCoord]), Form(";#Delta%c (mm)", lCoords[iCoord]), 500, -10, 10);
    GetOutputList()->Add(fHistResDecayVtx[iCoord]);
  }

  /// DCA to primary vertex
  for (int iSpecies = 0; iSpecies < 3; iSpecies++) {
    fHistDCA2pvXY[iSpecies] = new TH1D(Form("fDCA2PrimaryvtxXY_%s", lSpecies[iSpecies]), ";DCA_{xy} (mm);Counts", 600, 0, 30);
    fHistDCA2pvZ[iSpecies]  = new TH1D(Form("fDCA2PrimaryvtxZ_%s", lSpecies[iSpecies]), ";DCA_{z} (mm);Counts", 600, 0, 30);
    fHistDCA2pv[iSpecies]   = new TH1D(Form("fDCA2Primaryvtx_%s", lSpecies[iSpecies]), ";DCA (mm);Counts", 600, 0, 30);
    GetOutputList()->Add(fHistDCA2pvXY[iSpecies]);
    GetOutputList()->Add(fHistDCA2pvZ[iSpecies]);
    GetOutputList()->Add(fHistDCA2pv[iSpecies]);
  }

  /// DCA to secondary vertex
  for (int iSpecies = 0; iSpecies < 3; iSpecies++) {
    fHistDCA2dvXY[iSpecies] = new TH1D(Form("fDCA2DecayvtxXY_%s", lSpecies[iSpecies]), ";DCA_{xy} (mm);Counts", 600, 0, 30);
    fHistDCA2dvZ[iSpecies]  = new TH1D(Form("fDCA2DecayvtxZ_%s", lSpecies[iSpecies]), ";DCA_{z} (mm);Counts", 600, 0, 30);
    fHistDCA2dv[iSpecies]   = new TH1D(Form("fDCA2Decayvtx_%s", lSpecies[iSpecies]), ";DCA (mm);Counts", 600, 0, 30);
    GetOutputList()->Add(fHistDCA2dvXY[iSpecies]);
    GetOutputList()->Add(fHistDCA2dvZ[iSpecies]);
    GetOutputList()->Add(fHistDCA2dv[iSpecies]);
  }

  /// track distance to decay vertex histo
  for (int iCoord = 0; iCoord < 3; iCoord++) {
    fHistTrackDistance[iCoord] = new TH1D(Form("lTrackDistance_%s-%s", lSpeciesPair[iCoord], lSpeciesPair[(iCoord + 1) % 3]), ";distance [mm]", 200, 0, 200);
    GetOutputList()->Add(fHistTrackDistance[iCoord]);
  }

  /// histograms for efficiencies
  for(int iMatter=0; iMatter<2; iMatter++){
    fHistGen[iMatter] = new TH2D(Form("fHistGen_%c",lAM[iMatter]),";#it{p}_{T} (GeV/#it{c});#it{ct} (cm);",100,0.,10.,100,0.,100.);
    fHistRecSingle[iMatter] = new TH2D(Form("fHistRecSingle_%c",lAM[iMatter]),";#it{p}_{T} (GeV/#it{c});#it{ct} (cm);",100,0.,10.,100,0.,100.);
    fHistRecFake[iMatter] = new TH2D(Form("fHistRecFake_%c",lAM[iMatter]),";#it{p}_{T} (GeV/#it{c});#it{ct} (cm);",100,0.,10.,100,0.,100.);
    fHistRecClones[iMatter] = new TH2D(Form("fHistRecClones_%c",lAM[iMatter]),";#it{p}_{T} (GeV/#it{c});#it{ct} (cm);",100,0.,10.,100,0.,100.);

    GetOutputList()->Add(fHistGen[iMatter]);
    GetOutputList()->Add(fHistRecSingle[iMatter]);
    GetOutputList()->Add(fHistRecFake[iMatter]);
    GetOutputList()->Add(fHistRecClones[iMatter]);
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
  bool lIsGoodITS = true;
  bool lIsGoodTPCpid = true;
  bool lIsGoodTOFpid = true;
  bool lIsGoodDCA = true;

  bool lIsGoodVertex     = false;
  bool lIsGoodPAngle     = false;
  bool lIsGoodVertexHard = false;
  bool lIsGoodPAngleHard = false;

  double lPrimaryVertexCov[6]   = {0.};
  double lRecDecayVtx[3]        = {0.};
  double lRecDecayVtxHard[3]    = {0.};
  double lRecDecayLenght[3]     = {0.};
  double lRecDecayLenghtHard[3] = {0.};

  float lDCAmax[3] = {0.05, 0.05, 0.1};

  TLorentzVector lLVhyp = {0., 0., 0., 0.};
  TLorentzVector lLVdaughter[3];

  const float lMasses[3]{AliPID::ParticleMass(AliPID::kDeuteron), AliPID::ParticleMass(AliPID::kProton), AliPID::ParticleMass(AliPID::kPion)};

  //------------------------------------------------------------
  // Get stuff from tree
  //------------------------------------------------------------
  float lMagField = *fTreeHyp3BodyVarMagneticField;

  double lRecPrimaryVtx[3] = {0.};
  fPrimaryVertex->GetXYZ(lRecPrimaryVtx);
  double lTruePrimaryVtx[3]{*fTreeHyp3BodyVarPVtx[0],*fTreeHyp3BodyVarPVtx[1],*fTreeHyp3BodyVarPVtx[2]};

  AliESDtrack *lTrack[3]{&*fTreeHyp3BodyVarTracks[0],&*fTreeHyp3BodyVarTracks[1],&*fTreeHyp3BodyVarTracks[2]};

  int lNclsITS[3]{lTrack[0]->GetITSNcls(), lTrack[1]->GetITSNcls(), lTrack[2]->GetITSNcls()};
  int lNclsTPC[3]{lTrack[0]->GetTPCNcls(), lTrack[1]->GetTPCNcls(), lTrack[2]->GetTPCNcls()};

  float lNsigmaTPC[3] = {*fTreeHyp3BodyVarNsigmaTPC[0],*fTreeHyp3BodyVarNsigmaTPC[1],*fTreeHyp3BodyVarNsigmaTPC[2]};
  float lNsigmaTOF[3] = {*fTreeHyp3BodyVarNsigmaTOF[0],*fTreeHyp3BodyVarNsigmaTOF[1],*fTreeHyp3BodyVarNsigmaTOF[2]};

  double lTrueDecayVtx[3]{*fTreeHyp3BodyVarDecayVtx[0],*fTreeHyp3BodyVarDecayVtx[1],*fTreeHyp3BodyVarDecayVtx[2]};

  fCurrentEventId = *fTreeHyp3BodyVarEventId;
  fCurrentMotherId = *fTreeHyp3BodyVarMotherId;
  fFakeCand = *fTreeHyp3BodyVarIsFakeCand;
  

  /// charge of the candidate
  bool lCharge = lTrack[0]->GetSign() < 0;

  /// create the TLorentzVector of the hyper-triton candidate
  for (int iTrack = 0; iTrack < 3; iTrack++) {
    lLVdaughter[iTrack].SetXYZM(lTrack[iTrack]->Px(), lTrack[iTrack]->Py(), lTrack[iTrack]->Pz(), lMasses[iTrack]);
    lLVhyp += lLVdaughter[iTrack];
  }

  double lHypMassRec = lLVhyp.M();
  double lHypPtRec   = lLVhyp.Pt();

  //------------------------------------------------------------
  //                  Generated hypertritons                   
  //------------------------------------------------------------

  double lHypPtGen = std::hypot(*fTreeHyp3BodyVarTrueP[0],*fTreeHyp3BodyVarTrueP[1]);
  double lHypCtGen = Norm(lTrueDecayVtx[0]-lTruePrimaryVtx[0],lTrueDecayVtx[1]-lTruePrimaryVtx[1],lTrueDecayVtx[2]-lTruePrimaryVtx[2]);

  if(fCurrentMotherId != fLastMotherId){
    fLastMotherId = fCurrentMotherId;
    fHistGen[lCharge]->Fill(lHypPtGen,lHypCtGen);
    if(fCurrentEventId != fLastEventId) fLastEventId = fCurrentEventId;
    fNclones = 0;
  } else {
    if(fCurrentEventId != fLastEventId){
      fLastEventId = fCurrentEventId;
      fHistGen[lCharge]->Fill(lHypPtGen,lHypCtGen);
      fNclones = 0;
    }
  }

  //------------------------------------------------------------
  // N cluster in ITS, TPC and Track Chi2 distributions
  // and ITS cluster and PID selectios
  //------------------------------------------------------------

  for (int iTrack = 0; iTrack < 3; iTrack++) {
    fHistNclsITS[iTrack]->Fill(lNclsITS[iTrack] - 0.5);
    fHistNclsTPC[iTrack]->Fill(lNclsTPC[iTrack] - 0.5);
    lIsGoodITS &= lNclsITS[iTrack] > 0;

    fHistNSigmaTPC[iTrack]->Fill(lNsigmaTPC[iTrack]);
    fHistNSigmaTOF[iTrack]->Fill(lNsigmaTOF[iTrack]);
    lIsGoodTPCpid &= std::abs(lNsigmaTPC[iTrack]) < 3;
    lIsGoodTOFpid &= std::abs(lNsigmaTOF[iTrack]) < 3;
  }

  fHistNclsITSCheck->Fill((int)lIsGoodITS - 0.5);

  fHistCheckPID->Fill((int)lIsGoodTPCpid - 0.5);

  //------------------------------------------------------------
  // Secondary vertex reconstruction
  //------------------------------------------------------------
  // reconstruct the decay vertex with the dedicated vertexer
  lIsGoodVertex     = fHypertritonVertexer.FindDecayVertex(lTrack[0], lTrack[1], lTrack[2], lMagField);
  lIsGoodVertexHard = fHypertritonVertexerHard.FindDecayVertex(lTrack[0], lTrack[1], lTrack[2], lMagField);

  AliESDVertex *lDecayVertex     = static_cast<AliESDVertex *>(fHypertritonVertexer.GetCurrentVertex());
  AliESDVertex *lDecayVertexHard = static_cast<AliESDVertex *>(fHypertritonVertexerHard.GetCurrentVertex());    

  /// compute decay vertex position and compute cos(pointingangle)
  if (lIsGoodVertex) {
    
    lDecayVertex->GetXYZ(lRecDecayVtx);

    for (int iCoord = 0; iCoord < 3; iCoord++) {
      fHistResDecayVtx[iCoord]->Fill((lRecDecayVtx[iCoord] - lTrueDecayVtx[iCoord]) * 10.);
      lRecDecayLenght[iCoord] = lRecDecayVtx[iCoord] - lRecPrimaryVtx[iCoord];
    }

    TVector3 vRecDecayLenght(lRecDecayLenght[0], lRecDecayLenght[1], lRecDecayLenght[2]);
    float lHypCtRec = vRecDecayLenght.Mag();

    // Fake candidates
    if(fFakeCand){
      fHistRecFake[lCharge]->Fill(lHypPtRec,lHypCtRec);
    } else {
      if(fNclones==0){
        fHistRecSingle[lCharge]->Fill(lHypPtRec,lHypCtRec);
        fNclones++;
      } else {
        fHistRecClones[lCharge]->Fill(lHypPtRec,lHypCtRec);
      }
    }

    float pointAngle = lLVhyp.Angle(vRecDecayLenght);
    float cospa      = std::cos(pointAngle);

    fHistCosPAngle->Fill(cospa);
    lIsGoodPAngle = (cospa > 0.98);

    /// compute the DCA of the 3 tracks from the primary and decay vertex
    AliESDVertex lPV(lTruePrimaryVtx, lPrimaryVertexCov, 1., 1000);

    for (int iTrack = 0; iTrack < 3; iTrack++) {
      double dca2dv[2]    = {0.};
      double dca2pv[2]    = {0.};
      double dca2dvcov[3] = {0.};
      double dca2pvcov[3] = {0.};

      lTrack[iTrack]->PropagateToDCA(lDecayVertex, lMagField, 1000., dca2dv, dca2dvcov);
      lTrack[iTrack]->PropagateToDCA(&lPV, lMagField, 1000., dca2pv, dca2pvcov);

      float dcaXYdv = std::abs(dca2dv[0]) * 10.;    // in mm
      float dcaZdv  = std::abs(dca2dv[1]) * 10.;    // in mm
      float dcadv   = Norm(dcaXYdv, dcaZdv) * 10.;  // in mm

      fHistDCA2dvXY[iTrack]->Fill(dcaXYdv);
      fHistDCA2dvZ[iTrack]->Fill(dcaZdv);
      fHistDCA2dv[iTrack]->Fill(dcadv);

      float dcaXYpv = std::abs(dca2pv[0]) * 10.;    // in mm
      float dcaZpv  = std::abs(dca2pv[1]) * 10.;    // in mm
      float dcapv   = Norm(dcaXYpv, dcaZpv) * 10.;  // in mm

      lIsGoodDCA &= dcapv > lDCAmax[iTrack];

      fHistDCA2pvXY[iTrack]->Fill(dcaXYpv);
      fHistDCA2pvZ[iTrack]->Fill(dcaZpv);
      fHistDCA2pv[iTrack]->Fill(dcapv);
    }

    /// compute the track2track distance used in the vertexer
    float pPM[3][3];

    for (int iPerm = 0; iPerm < 3; iPerm++) {
      fHypertritonVertexer.Find2ProngClosestPoint(lTrack[iPerm], lTrack[(iPerm + 1) % 3], lMagField, pPM[iPerm]);
    }

    for (int iPerm = 0; iPerm < 3; iPerm++) {
      float distance = Point2PointDistance(pPM[iPerm], pPM[(iPerm + 1) % 3]);
      fHistTrackDistance[iPerm]->Fill(distance * 10.);
    }

    double vertexChi2NDF = lDecayVertex->GetChi2perNDF();
    fHistVertexChi2->Fill(vertexChi2NDF);
  }

  /// again with harder selections
  if (lIsGoodVertexHard) {
    lDecayVertexHard->GetXYZ(lRecDecayVtxHard);

    for (int iTrack = 0; iTrack < 3; iTrack++) {
      lRecDecayLenghtHard[iTrack] = lRecDecayVtxHard[iTrack] - lTruePrimaryVtx[iTrack];
    }

    TVector3 vRecDecayLenght(lRecDecayLenghtHard[0], lRecDecayLenghtHard[1], lRecDecayLenghtHard[2]);

    float pointAngle = lLVhyp.Angle(vRecDecayLenght);
    float cospa      = std::cos(pointAngle);

    lIsGoodPAngleHard = (cospa > 0.998);
  }

  bool lIsGoodCuts     = lIsGoodDCA && lIsGoodVertex && lIsGoodPAngle;
  bool lIsGoodHardCuts = lIsGoodDCA && lIsGoodVertexHard && lIsGoodPAngleHard;

  /// invariant mass distribution for efficiency
  fHistInvMass[lCharge][0]->Fill(lHypMassRec, lHypPtRec);
  if (lIsGoodCuts) fHistInvMass[lCharge][1]->Fill(lHypMassRec, lHypPtRec);
  if (lIsGoodHardCuts) fHistInvMass[lCharge][2]->Fill(lHypMassRec, lHypPtRec);

  /// transverse momentum distribution
  fHistPt[lCharge][0]->Fill(lHypPtRec);
  if (lIsGoodCuts) fHistPt[lCharge][1]->Fill(lHypPtRec);
  if (lIsGoodHardCuts) fHistPt[lCharge][2]->Fill(lHypPtRec);

  /// daughter transverse momentum distribution
  for (int iSpecies = 0; iSpecies < 3; iSpecies++) {
    double pt = lTrack[iSpecies]->Pt() * 1.;
    fHistDaughterPt[iSpecies][0]->Fill(pt);
    if (lIsGoodCuts) fHistDaughterPt[iSpecies][1]->Fill(pt);
    if (lIsGoodHardCuts) fHistDaughterPt[iSpecies][2]->Fill(pt);
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