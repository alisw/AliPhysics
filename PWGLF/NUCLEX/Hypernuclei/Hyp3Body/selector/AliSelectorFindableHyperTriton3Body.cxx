#define AliSelectorFindableHyperTriton3Body_cxx

#include "Hyp3FindConfig.h"
#include "AliESDVertex.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
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

AliSelectorFindableHyperTriton3Body::AliSelectorFindableHyperTriton3Body(TString outputName, TString outputPath,TTree *):
fOutputFileName{outputName},
fOutputFilePath{outputPath}, 
fHypertritonVertexer(){
  fESDtrackCuts = AliESDtrackCuts::GetStandardV0DaughterCuts();
}

//______________________________________________________________________________
void AliSelectorFindableHyperTriton3Body::Begin(TTree * /*tree*/){
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).

  TString option = GetOption();
}

//______________________________________________________________________________
void AliSelectorFindableHyperTriton3Body::SlaveBegin(TTree * /*tree*/){
  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).

  TString option = GetOption();

  const char lAM[3]{"AM"};
  const char *lSpecies[3]{"d", "p", "pi"};
  const char *lCuts[3]{"tree", "selection", "vertexer"};
  const char lCoords[3]{'x','y','z'};

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

  /// Histograms for selection

  for(int iCuts = 0; iCuts < 3; iCuts++){
    for(int iMatter = 0; iMatter < 2; iMatter++){
      fHistInvMassPt[iMatter][iCuts] = new TH2D(Form("fHistInvMassPt_%c_%s", lAM[iMatter], lCuts[iCuts]), ";#it{m} (dp#pi) (GeV/#it{c}^{2});#it{p}_{T} (GeV/#it{c})", 200, 2.95, 3.35, 20, 0, 10);
      GetOutputList()->Add(fHistInvMassPt[iMatter][iCuts]);
    }
    for(int iSpecies = 0; iSpecies < 3; iSpecies++){
      fHistDaughterPt[iSpecies][iCuts] = new TH1D(Form("fHistPt_%s_%s", lSpecies[iSpecies], lCuts[iCuts]), ";#it{p}_{T} (GeV/#it{c});Counts", 100, 0, 10);
      fHistDaughterTPCchi2[iSpecies][iCuts] = new TH1D(Form("fHistTPCchi2_%s_%s", lSpecies[iSpecies], lCuts[iCuts]), ";#chi^{2};Counts", 100, 0, 10);
      fHistDaughterITSchi2[iSpecies][iCuts] = new TH1D(Form("fHistITSchi2_%s_%s", lSpecies[iSpecies], lCuts[iCuts]), ";#chi^{2};Counts", 100, 0, 10);
      fHistNclsITS[iSpecies][iCuts] = new TH1D(Form("fHistNclsITS_%s_%s", lSpecies[iSpecies], lCuts[iCuts]), ";N_{cluster ITS};Counts", 7, -0.5, 6.5);
      fHistNclsTPC[iSpecies][iCuts] = new TH1D(Form("fHistNclsTPC_%s_%s", lSpecies[iSpecies], lCuts[iCuts]), ";N_{cluster TPC};Counts", 101, -0.5, 200.5);
      fHistNSigmaTPC[iSpecies][iCuts] = new TH1D(Form("fHistNSigmaTPC_%s_%s", lSpecies[iSpecies], lCuts[iCuts]), "n#sigma_{TPC}", 10, 0, 10);
      fHistNSigmaTOF[iSpecies][iCuts] = new TH1D(Form("fHistNSigmaTOF_%s_%s", lSpecies[iSpecies], lCuts[iCuts]), "n#sigma_{TOF}", 10, 0, 10);
      GetOutputList()->Add(fHistDaughterPt[iSpecies][iCuts]);
      GetOutputList()->Add(fHistDaughterTPCchi2[iSpecies][iCuts]);
      GetOutputList()->Add(fHistDaughterITSchi2[iSpecies][iCuts]);
      GetOutputList()->Add(fHistNclsITS[iSpecies][iCuts]);
      GetOutputList()->Add(fHistNclsTPC[iSpecies][iCuts]);
      GetOutputList()->Add(fHistNSigmaTPC[iSpecies][iCuts]);
      GetOutputList()->Add(fHistNSigmaTOF[iSpecies][iCuts]);
    } 
  }

  /// Histograms after vertexer

  fHistVertexChi2 = new TH1D("fHistVertexChi2", "", 100, 0, 200);
  fHistCosPAngle  = new TH1D("fCosPointingAngle", ";#it{cos#theta_{pointing}};Counts", 5000, 0.5, 1.);
  GetOutputList()->Add(fHistVertexChi2);
  GetOutputList()->Add(fHistCosPAngle);

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
    fHistDCA2dvXY[iSpecies] = new TH1D(Form("fDCA2DecayvtxXY_%s", lSpecies[iSpecies]), ";DCA_{xy} (mm);Counts", 600, 0, 30);
    fHistDCA2dvZ[iSpecies]  = new TH1D(Form("fDCA2DecayvtxZ_%s", lSpecies[iSpecies]), ";DCA_{z} (mm);Counts", 600, 0, 30);
    fHistDCA2dv[iSpecies]   = new TH1D(Form("fDCA2Decayvtx_%s", lSpecies[iSpecies]), ";DCA (mm);Counts", 600, 0, 30);
    fHistTrackDistance[iSpecies] = new TH1D(Form("lTrackDistance_%s-%s",lSpecies[iSpecies], lSpecies[(iSpecies + 1) % 3]), ";distance (mm);counts", 200, 0, 200);
    GetOutputList()->Add(fHistDCA2pvXY[iSpecies]);
    GetOutputList()->Add(fHistDCA2pvZ[iSpecies]);
    GetOutputList()->Add(fHistDCA2pv[iSpecies]);
    GetOutputList()->Add(fHistDCA2dvXY[iSpecies]);
    GetOutputList()->Add(fHistDCA2dvZ[iSpecies]);
    GetOutputList()->Add(fHistDCA2dv[iSpecies]);
    GetOutputList()->Add(fHistTrackDistance[iSpecies]);
  }
}

//______________________________________________________________________________
Bool_t AliSelectorFindableHyperTriton3Body::Process(Long64_t entry){
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

  // Support variables

  double lPrimaryVtxCov[6] = {0.};
  double lRecPrimaryVtx[3] = {0.};
  double lRecDecayVtx[3] = {0.};
  double lRecDecayLenght[3] = {0.};

  TLorentzVector lLVhyp = {0., 0., 0., 0.};
  TLorentzVector lLVdaughter[3];

  const float lMasses[3]{AliPID::ParticleMass(AliPID::kDeuteron), AliPID::ParticleMass(AliPID::kProton), AliPID::ParticleMass(AliPID::kPion)};

  //------------------------------------------------------------
  // Get main observables from tree
  //------------------------------------------------------------
  float lMagField = *fTreeHyp3BodyVarMagneticField;
  fPrimaryVertex->GetXYZ(lRecPrimaryVtx);
  double lTruePrimaryVtx[3]{*fTreeHyp3BodyVarPVtx[0],*fTreeHyp3BodyVarPVtx[1],*fTreeHyp3BodyVarPVtx[2]};

  AliESDtrack *lTrack[3]{&*fTreeHyp3BodyVarTracks[0],&*fTreeHyp3BodyVarTracks[1],&*fTreeHyp3BodyVarTracks[2]};
  int lDaughterNclsITS[3]{lTrack[0]->GetITSNcls(), lTrack[1]->GetITSNcls(), lTrack[2]->GetITSNcls()};
  int lDaughterNclsTPC[3]{lTrack[0]->GetTPCNcls(), lTrack[1]->GetTPCNcls(), lTrack[2]->GetTPCNcls()};
  double lDaughterChi2TPC[3]{lTrack[0]->GetTPCchi2(), lTrack[1]->GetTPCchi2(), lTrack[2]->GetTPCchi2()};
  double lDaughterChi2ITS[3]{lTrack[0]->GetITSchi2(), lTrack[1]->GetITSchi2(), lTrack[2]->GetITSchi2()};
  float lDaughterNsigmaTPC[3] = {*fTreeHyp3BodyVarNsigmaTPC[0],*fTreeHyp3BodyVarNsigmaTPC[1],*fTreeHyp3BodyVarNsigmaTPC[2]};
  float lDaughterNsigmaTOF[3] = {*fTreeHyp3BodyVarNsigmaTOF[0],*fTreeHyp3BodyVarNsigmaTOF[1],*fTreeHyp3BodyVarNsigmaTOF[2]};

  double lTrueDecayVtx[3]{*fTreeHyp3BodyVarDecayVtx[0],*fTreeHyp3BodyVarDecayVtx[1],*fTreeHyp3BodyVarDecayVtx[2]};
  fCurrentEventId = *fTreeHyp3BodyVarEventId;
  fCurrentMotherId = *fTreeHyp3BodyVarMotherId;
  fFakeCand = *fTreeHyp3BodyVarIsFakeCand;
  bool lCharge = lTrack[0]->GetSign() < 0;

  float lDaughterPt[3] = {0.};
  for (int iTrack = 0; iTrack < 3; iTrack++) {
    lLVdaughter[iTrack].SetXYZM(lTrack[iTrack]->Px(), lTrack[iTrack]->Py(), lTrack[iTrack]->Pz(), lMasses[iTrack]);
    lDaughterPt[iTrack] = lLVdaughter[iTrack].Pt();
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
  //                    Before selection
  //------------------------------------------------------------

  fHistInvMassPt[lCharge][0]->Fill(lHypMassRec,lHypPtRec);
  for(int iTrack = 0; iTrack < 3; iTrack++){
    fHistDaughterPt[iTrack][0]->Fill(lDaughterPt[iTrack]);
    fHistNclsITS[iTrack][0]->Fill(lDaughterNclsITS[iTrack]);
    fHistNclsTPC[iTrack][0]->Fill(lDaughterNclsTPC[iTrack]);
    fHistNSigmaTPC[iTrack][0]->Fill(lDaughterNsigmaTPC[iTrack]);
    fHistNSigmaTOF[iTrack][0]->Fill(lDaughterNsigmaTOF[iTrack]);
    fHistDaughterTPCchi2[iTrack][0]->Fill(lDaughterChi2TPC[iTrack]);
    fHistDaughterITSchi2[iTrack][0]->Fill(lDaughterChi2ITS[iTrack]);
  }

  //------------------------------------------------------------
  //                      After selection
  //------------------------------------------------------------

  if(AcceptCandidate()){
    fHistInvMassPt[lCharge][1]->Fill(lHypMassRec,lHypPtRec);
    for(int iTrack = 0; iTrack < 3; iTrack++){
      fHistDaughterPt[iTrack][1]->Fill(lDaughterPt[iTrack]);
      fHistNclsITS[iTrack][1]->Fill(lDaughterNclsITS[iTrack]);
      fHistNclsTPC[iTrack][1]->Fill(lDaughterNclsTPC[iTrack]);
      fHistNSigmaTPC[iTrack][1]->Fill(lDaughterNsigmaTPC[iTrack]);
      fHistNSigmaTOF[iTrack][1]->Fill(lDaughterNsigmaTOF[iTrack]);
      fHistDaughterTPCchi2[iTrack][1]->Fill(lDaughterChi2TPC[iTrack]);
      fHistDaughterITSchi2[iTrack][1]->Fill(lDaughterChi2ITS[iTrack]);
    }
  } else return true;

  //------------------------------------------------------------
  // Secondary vertex reconstruction
  //------------------------------------------------------------
  // reconstruct the decay vertex with the dedicated vertexer
  if(!fHypertritonVertexer.FindDecayVertex(lTrack[0], lTrack[1], lTrack[2], lMagField)) return true;

  fHistInvMassPt[lCharge][2]->Fill(lHypMassRec,lHypPtRec);
  for(int iTrack = 0; iTrack < 3; iTrack++){
    fHistDaughterPt[iTrack][2]->Fill(lDaughterPt[iTrack]);
    fHistNclsITS[iTrack][2]->Fill(lDaughterNclsITS[iTrack]);
    fHistNclsTPC[iTrack][2]->Fill(lDaughterNclsTPC[iTrack]);
    fHistNSigmaTPC[iTrack][2]->Fill(lDaughterNsigmaTPC[iTrack]);
    fHistNSigmaTOF[iTrack][2]->Fill(lDaughterNsigmaTOF[iTrack]);
    fHistDaughterTPCchi2[iTrack][2]->Fill(lDaughterChi2TPC[iTrack]);
    fHistDaughterITSchi2[iTrack][2]->Fill(lDaughterChi2ITS[iTrack]);
  }
  
  AliESDVertex *lDecayVertex = static_cast<AliESDVertex *>(fHypertritonVertexer.GetCurrentVertex());   
  lDecayVertex->GetXYZ(lRecDecayVtx);
  for (int iCoord = 0; iCoord < 3; iCoord++) {
    fHistResDecayVtx[iCoord]->Fill((lRecDecayVtx[iCoord] - lTrueDecayVtx[iCoord]) * 10.);
    lRecDecayLenght[iCoord] = lRecDecayVtx[iCoord] - lRecPrimaryVtx[iCoord];
  }
  TVector3 vRecDecayLenght(lRecDecayLenght[0], lRecDecayLenght[1], lRecDecayLenght[2]);
  float lHypCtRec = vRecDecayLenght.Mag();

  /// Efficiency histograms
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

  /// compute the DCA of the 3 tracks from the primary and decay vertex
  AliESDVertex lPV(lTruePrimaryVtx, lPrimaryVtxCov, 1., 1000);

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

    fHistDCA2pvXY[iTrack]->Fill(dcaXYpv);
    fHistDCA2pvZ[iTrack]->Fill(dcaZpv);
    fHistDCA2pv[iTrack]->Fill(dcapv);
  }
  /// compute the track2track distance used in the vertexer
  float pPM[3][3];
  for (int iPerm = 0; iPerm < 3; iPerm++) {
    fHypertritonVertexer.Find2ProngClosestPoint(lTrack[iPerm], lTrack[(iPerm + 1) % 3], lMagField, pPM[iPerm]);
    fHistTrackDistance[iPerm]->Fill(Point2PointDistance(pPM[iPerm], pPM[(iPerm + 1) % 3]) * 10.);
  }
  double vertexChi2NDF = lDecayVertex->GetChi2perNDF();
  fHistVertexChi2->Fill(vertexChi2NDF);

  return true;
}

//______________________________________________________________________________
void AliSelectorFindableHyperTriton3Body::SlaveTerminate(){
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.
}

//______________________________________________________________________________
void AliSelectorFindableHyperTriton3Body::Terminate(){
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.
  TFile output(Form("%s/%s", fOutputFilePath.Data(), fOutputFileName.Data()), "RECREATE");
  GetOutputList()->Write();
  output.Close();
}

//______________________________________________________________________________
bool AliSelectorFindableHyperTriton3Body::AcceptCandidate(){
  for(int iTrack = 0; iTrack < 3; iTrack++){
    if(!fESDtrackCuts->AcceptTrack(&*fTreeHyp3BodyVarTracks[iTrack])) return false;
    if(*fTreeHyp3BodyVarNsigmaTPC[iTrack] > kNsigmaTPC[iTrack]) return false;
    if((*fTreeHyp3BodyVarTracks)->GetITSNcls() < kNclusITS[iTrack]) return false;
    if((*fTreeHyp3BodyVarTracks)->GetTPCNcls() < kNclusTPC[iTrack]) return false;
  }
  return true;
}