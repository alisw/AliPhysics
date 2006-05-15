#include "AlidNdEtaEffSelector.h"

#include <TStyle.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TParticle.h>
#include <TParticlePDG.h>

#include <AliLog.h>
#include <../PYTHIA6/AliGenPythiaEventHeader.h>
#include <../TPC/AliTPCtrack.h>
#include <AliTracker.h>

#include <iostream>
using namespace std;

ClassImp(AlidNdEtaEffSelector)

AlidNdEtaEffSelector::AlidNdEtaEffSelector(TTree *) :
  AliSelector(),
  fEsdTrackCuts(0),
  fdNdEtaCorrection(0)
{
  // Constructor. Initialization of pointers
}

AlidNdEtaEffSelector::~AlidNdEtaEffSelector()
{
  // Remove all pointers

  // histograms are in the output list and deleted when the output
  // list is deleted by the TSelector dtor
}

void AlidNdEtaEffSelector::Begin(TTree * tree)
{
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).

  AliSelector::Begin(tree);

  fdNdEtaCorrection = new dNdEtaCorrection();
}

void AlidNdEtaEffSelector::SlaveBegin(TTree * tree)
{
  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).

  AliSelector::SlaveBegin(tree);

  fEsdTrackCuts = new ESDtrackQualityCuts();
  fEsdTrackCuts->DefineHistograms(1);

  fEsdTrackCuts->SetMinNClustersTPC(50);
  fEsdTrackCuts->SetMaxChi2PerClusterTPC(3.5);
  fEsdTrackCuts->SetMaxCovDiagonalElements(2,2,0.5,0.5,2);
  fEsdTrackCuts->SetRequireTPCRefit(kTRUE);

  fEsdTrackCuts->SetMinNsigmaToVertex(3);
  fEsdTrackCuts->SetAcceptKingDaughters(kFALSE);

  AliLog::SetClassDebugLevel("ESDtrackQualityCuts",1);
}

Bool_t AlidNdEtaEffSelector::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. Typically here the branch pointers
  // will be retrieved. It is normaly not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed.

  if (AliSelector::Notify() == kFALSE)
    return kFALSE;

  // ########################################################
  // Magnetic field
  AliTracker::SetFieldMap(GetAliRun()->Field(), kTRUE); // kTRUE means uniform magnetic field

  return kTRUE;
}

Bool_t AlidNdEtaEffSelector::IsPrimary(const TParticle* aParticle, Int_t aTotalPrimaries)
{
  // if the particle has a daughter primary, we do not want to count it
  if (aParticle->GetFirstDaughter() != -1 && aParticle->GetFirstDaughter() < aTotalPrimaries)
    return kFALSE;

  Int_t pdgCode = TMath::Abs(aParticle->GetPdgCode());

  // skip quarks and gluon
  if (pdgCode > 10 && pdgCode != 21)
    return kTRUE;

  return kFALSE;
}

Bool_t AlidNdEtaEffSelector::Process(Long64_t entry)
{
  // The Process() function is called for each entry in the tree (or possibly
  // keyed object in the case of PROOF) to be processed. The entry argument
  // specifies which entry in the currently loaded tree is to be processed.
  // It can be passed to either TTree::GetEntry() or TBranch::GetEntry()
  // to read either all or the required parts of the data. When processing
  // keyed objects with PROOF, the object is already loaded and is available
  // via the fObject pointer.
  //
  // This function should contain the "body" of the analysis. It can contain
  // simple or elaborate selection criteria, run algorithms on the data
  // of the event and typically fill histograms.

  // WARNING when a selector is used with a TChain, you must use
  //  the pointer to the current TTree to call GetEntry(entry).
  //  The entry is always the local entry number in the current tree.
  //  Assuming that fChain is the pointer to the TChain being processed,
  //  use fChain->GetTree()->GetEntry(entry).

  if (AliSelector::Process(entry) == kFALSE)
    return kFALSE;

  if (!fESD || !fHeader)
    return kFALSE;

  // ########################################################
  // get the EDS vertex
  const AliESDVertex* vtxESD = fESD->GetVertex();

  Double_t vtx[3];
  Double_t vtx_res[3];
  vtxESD->GetXYZ(vtx);

  vtx_res[0] = vtxESD->GetXRes();
  vtx_res[1] = vtxESD->GetYRes();
  vtx_res[2] = vtxESD->GetZRes();

  // the vertex should be reconstructed
  if (strcmp(vtxESD->GetName(),"default")==0)
    return kTRUE;

  // the resolution should be reasonable???
  if (vtx_res[2]==0 || vtx_res[2]>0.1)
    return kTRUE;

  // ########################################################
  // get the MC vertex
  AliGenPythiaEventHeader* genHeader = (AliGenPythiaEventHeader*) fHeader->GenEventHeader();

  TArrayF vtxMC(3);
  genHeader->PrimaryVertex(vtxMC);
  vtx[0] = vtxMC[0];
  vtx[1] = vtxMC[1];
  vtx[2] = vtxMC[2];

  // ########################################################
  // loop over mc particles
  TTree* particleTree = GetKinematics();
  TParticle* particle = 0;
  particleTree->SetBranchAddress("Particles", &particle);

  Int_t nPrim  = fHeader->GetNprimary();
  Int_t nTotal = fHeader->GetNtrack();

  for (Int_t i_mc = nTotal - nPrim; i_mc < nTotal; ++i_mc)
  {
    particleTree->GetEntry(i_mc);

    if (!particle)
      continue;

    if (strcmp(particle->GetName(),"XXX") == 0)
    {
      printf("WARNING: There is a particle named XXX (%d).\n", i_mc);
      continue;
    }

    TParticlePDG* pdgPart = particle->GetPDG();

    if (strcmp(pdgPart->ParticleClass(),"Unknown") == 0)
    {
      printf("WARNING: There is a particle with an unknown particle class (%d pdg code %d).\n", i_mc, particle->GetPdgCode());
      continue;
    }

    if (IsPrimary(particle, nPrim) == kFALSE)
      continue;

    if (pdgPart->Charge() == 0)
      continue;

    fdNdEtaCorrection->FillGene(vtx[2], particle->Eta());

  }// end of mc particle

  // ########################################################
  // loop over esd tracks
  Int_t nTracks = fESD->GetNumberOfTracks();
  for (Int_t t=0; t<nTracks; t++)
  {
    AliESDtrack* esdTrack = fESD->GetTrack(t);

    // cut the esd track?
    if (!fEsdTrackCuts->AcceptTrack(esdTrack))
      continue;

    AliTPCtrack* tpcTrack = new AliTPCtrack(*esdTrack);
    if (tpcTrack->GetAlpha()==0)
    {
      cout << " Warning esd track alpha = 0" << endl;
      continue;
    }

    //Float_t theta = tpcTrack->Theta();
    //Float_t eta   = -TMath::Log(TMath::Tan(theta/2.));

    // using the eta of the mc particle
    Int_t label = TMath::Abs(esdTrack->GetLabel());
    if (label<0) 
    {
      cout << " cannot find corresponding mc part !!! " << label << endl;
      continue;
    }
    particleTree->GetEntry(nTotal - nPrim + label);

    fdNdEtaCorrection->FillMeas(vtx[2], particle->Eta());

  } // end of track loop

  return kTRUE;
}

void AlidNdEtaEffSelector::SlaveTerminate()
{
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.

  AliSelector::SlaveTerminate();

  // Add the histograms to the output on each slave server
  if (!fOutput)
  {
    printf("ERROR: Output list not initialized\n");
    return;
  }
}

void AlidNdEtaEffSelector::Terminate()
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

  AliSelector::Terminate();

  fdNdEtaCorrection->Finish();

  TFile* fout = new TFile("correction_map.root","RECREATE");

  fEsdTrackCuts->SaveHistograms("esd_track_cuts");
  fdNdEtaCorrection->SaveHistograms();

  fout->Write();
  fout->Close();

  fdNdEtaCorrection->DrawHistograms();
}
