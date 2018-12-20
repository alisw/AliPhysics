#define AliSelectorFindableHyperTriton3Body_cxx

#include "AliSelectorFindableHyperTriton3Body.h"
#include "AliPID.h"
#include <cmath>
#include <TH1D.h>
#include <TList.h>
#include <TLorentzVector.h>
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
   fHistInvMass = new TH1D("fHistInvMass",";m (dp#pi) [GeV/#it{c}^{2}];Counts",200,2.6,3.2);
   fHistClones = new TH1D("fHistClones",";#it{p}_{T} [GeV/#it{c}];Counts",100,0,10);

   GetOutputList()->Add(fHistInvMass);
   GetOutputList()->Add(fHistClones);
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

   bool lIsClone = (fLastMother >= 0) && (fCurrentEventId == *fTreeHyp3BodyVarEventId) && (fLastMother == *fTreeHyp3BodyVarMotherId);
   fCurrentEventId = *fTreeHyp3BodyVarEventId;
   fLastMother = *fTreeHyp3BodyVarMotherId;

   float lTruePt = std::hypot(*fTreeHyp3BodyVarTruePx, *fTreeHyp3BodyVarTruePy);
   float lMasses[3]{AliPID::ParticleMass(AliPID::kDeuteron),AliPID::ParticleMass(AliPID::kProton),AliPID::ParticleMass(AliPID::kPion)};
   if (lIsClone)
      fHistClones->Fill(lTruePt);
   else {
      TLorentzVector lVectors[2];
      for (int iVec = 0; iVec < 3; ++iVec) {
         lVectors[1].SetPtEtaPhiM(fTreeHyp3BodyVarTracks[iVec]->Pt(), fTreeHyp3BodyVarTracks[iVec]->Eta(), fTreeHyp3BodyVarTracks[iVec]->Phi(), lMasses[iVec]);
         lVectors[0] += lVectors[1];
      }
      fHistInvMass->Fill(lVectors[0].M());
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