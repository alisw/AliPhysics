#define AliSelectorFindableHyperTriton3Body_cxx

#include "AliSelectorFindableHyperTriton3Body.h"
#include "AliPID.h"
#include "AliVTrack.h"
#include <cmath>
#include <TH1D.h>
#include <TList.h>
#include <TLorentzVector.h>
#include <TObjArray.h>
#include <TStyle.h>

void AliSelectorFindableHyperTriton3Body::Begin(TTree * /*tree*/)
{
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).

  TString option = GetOption();
}

void AliSelectorFindableHyperTriton3Body::SlaveBegin(TTree * /*tree*/)
{
  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).

  TString option = GetOption();
  const char* lClone[2]{"","Clone"};
  for (int iClone = 0; iClone < 2; ++iClone) {
    fHistInvMass[iClone] = new TH1D(Form("fHistInvMass%s",lClone[iClone]),";m (dp#pi) [GeV/#it{c}^{2}];Counts",200,2.95,3.35);
    fHistPt[iClone] = new TH1D(Form("fHistPt%s",lClone[iClone]),";#it{p}_{T} [GeV/#it{c}];Counts",100,0,10);
    fHistVertexChi2[iClone] = new TH1D(Form("fHistVertexChi2%s",lClone[iClone]),";Vertex #chi^{2}/NDF; Counts",500,0,100);

    GetOutputList()->Add(fHistInvMass[iClone]);
    GetOutputList()->Add(fHistPt[iClone]);
    GetOutputList()->Add(fHistVertexChi2[iClone]);

    const char lCoords[4]{"xyz"};
    for (int iCoord = 0; iCoord < 3; ++iCoord) {
      fHistResDecayVtx[iCoord][iClone] = new TH1D(Form("fHistResDecayVtx%c%s",lCoords[iCoord],lClone[iClone]),Form(";#Delta%c [mm]",lCoords[iCoord]),500,-10,10);
      GetOutputList()->Add(fHistResDecayVtx[iCoord][iClone]);
    }

  }
}

Bool_t AliSelectorFindableHyperTriton3Body::Process(Long64_t entry)
{
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
  fVertexer.SetFieldkG(*fTreeHyp3BodyVarMagneticField);

  for (int iTrack = 0; iTrack < 3; ++iTrack) {
    if ((fTreeHyp3BodyVarTracks[iTrack]->GetStatus() & AliVTrack::kTPCrefit) == 0 ||
        fTreeHyp3BodyVarTracks[iTrack]->GetTPCNcls() < 50 ||
        fTreeHyp3BodyVarTracks[iTrack]->GetTPCchi2() > 4 * fTreeHyp3BodyVarTracks[iTrack]->GetTPCNcls() ||
        std::abs(fTreeHyp3BodyVarTracks[iTrack]->Eta()) > 0.9 ||
        fTreeHyp3BodyVarTracks[iTrack]->GetKinkIndex(0) > 0)
      return true;
  }

  bool lIsClone = (fLastMother >= 0) && (fCurrentEventId == *fTreeHyp3BodyVarEventId) && (fLastMother == *fTreeHyp3BodyVarMotherId);
  fCurrentEventId = *fTreeHyp3BodyVarEventId;
  fLastMother = *fTreeHyp3BodyVarMotherId;

  float lTruePt = std::hypot(*fTreeHyp3BodyVarTrueP[0], *fTreeHyp3BodyVarTrueP[1]);
  float lMasses[3]{AliPID::ParticleMass(AliPID::kDeuteron),AliPID::ParticleMass(AliPID::kProton),AliPID::ParticleMass(AliPID::kPion)};
  TObjArray lTrackArray(3);
  fHistPt[lIsClone]->Fill(lTruePt);
  TLorentzVector lVectors[2];
  for (int iVec = 0; iVec < 3; ++iVec) {
    lVectors[1].SetPtEtaPhiM(fTreeHyp3BodyVarTracks[iVec]->Pt(), fTreeHyp3BodyVarTracks[iVec]->Eta(), fTreeHyp3BodyVarTracks[iVec]->Phi(), lMasses[iVec]);
    lVectors[0] += lVectors[1];
    lTrackArray.Add(&*fTreeHyp3BodyVarTracks[iVec]);
  }
  fHistInvMass[lIsClone]->Fill(lVectors[0].M());
  AliESDVertex* lDecayVtx = fVertexer.VertexForSelectedESDTracks(&lTrackArray);
  fHistVertexChi2[lIsClone]->Fill(lDecayVtx->GetChi2perNDF());

  double lXYZ[3]{0.};
  lDecayVtx->GetXYZ(lXYZ);
  for (int iCoord{0}; iCoord < 3; ++iCoord) {
    fHistResDecayVtx[iCoord][lIsClone]->Fill((lXYZ[iCoord] - *fTreeHyp3BodyVarDecayVtx[iCoord]) * 10.);
  }

  return kTRUE;
}

void AliSelectorFindableHyperTriton3Body::SlaveTerminate()
{
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.

}

void AliSelectorFindableHyperTriton3Body::Terminate()
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.
  TFile output("output.root", "recreate");
  GetOutputList()->Write();
  output.Close();
}
