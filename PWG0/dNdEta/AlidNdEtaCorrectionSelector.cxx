/* $Id$ */

#include "AlidNdEtaCorrectionSelector.h"

#include <TStyle.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TParticle.h>
#include <TParticlePDG.h>
#include <TH1F.h>

#include <TChain.h>
#include <TSelector.h>
#include <TFile.h>

#include <AliLog.h>
#include <AliTracker.h>
#include <AliESDVertex.h>
#include <AliESD.h>
#include <AliESDtrack.h>
#include <AliRunLoader.h>
#include <AliStack.h>

#include <AliHeader.h>
#include <AliGenEventHeader.h>
#include <AliGenPythiaEventHeader.h>
#include <AliGenCocktailEventHeader.h>

#include "esdTrackCuts/AliESDtrackCuts.h"
#include "dNdEta/AlidNdEtaCorrection.h"
#include "AliPWG0Helper.h"

ClassImp(AlidNdEtaCorrectionSelector)

AlidNdEtaCorrectionSelector::AlidNdEtaCorrectionSelector() :
  AliSelectorRL(),
  fEsdTrackCuts(0),
  fdNdEtaCorrection(0),
  fPIDParticles(0),
  fPIDTracks(0),
  fClustersITSPos(0),
  fClustersTPCPos(0),
  fClustersITSNeg(0),
  fClustersTPCNeg(0),
  fSignMode(0)
{
  //
  // Constructor. Initialization of pointers
  //
}

AlidNdEtaCorrectionSelector::~AlidNdEtaCorrectionSelector()
{
  //
  // Destructor
  //

  // histograms are in the output list and deleted when the output
  // list is deleted by the TSelector dtor
}

Bool_t AlidNdEtaCorrectionSelector::SignOK(TParticlePDG* particle)
{
  // returns if a particle with this sign should be counted
  // this is determined by the value of fSignMode, which should have the same sign
  // as the charge
  // if fSignMode is 0 all particles are counted

  if (fSignMode == 0)
    return kTRUE;

  if (!particle)
  {
    printf("WARNING: not counting a particle that does not have a pdg particle\n");
    return kFALSE;
  }

  Double_t charge = particle->Charge();

  if (fSignMode > 0)
    if (charge < 0)
      return kFALSE;

  if (fSignMode < 0)
    if (charge > 0)
      return kFALSE;

  return kTRUE;
}

void AlidNdEtaCorrectionSelector::ReadUserObjects(TTree* tree)
{
  // read the user objects, called from slavebegin and begin

  if (!fEsdTrackCuts && fInput)
    fEsdTrackCuts = dynamic_cast<AliESDtrackCuts*> (fInput->FindObject("AliESDtrackCuts"));

  if (!fEsdTrackCuts && tree)
    fEsdTrackCuts = dynamic_cast<AliESDtrackCuts*> (tree->GetUserInfo()->FindObject("AliESDtrackCuts"));

  if (!fEsdTrackCuts)
     AliDebug(AliLog::kError, "ERROR: Could not read EsdTrackCuts from input list.");
}

void AlidNdEtaCorrectionSelector::Begin(TTree * tree)
{
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).

  AliSelectorRL::Begin(tree);

  ReadUserObjects(tree);

  TString option = GetOption();
  AliInfo(Form("Called with option %s.", option.Data()));

  if (option.Contains("only-positive"))
  {
    AliInfo("Processing only positive particles.");
    fSignMode = 1;
  }
  else if (option.Contains("only-negative"))
  {
    AliInfo("Processing only negative particles.");
    fSignMode = -1;
  }
}

void AlidNdEtaCorrectionSelector::SlaveBegin(TTree * tree)
{
  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).

  AliSelectorRL::SlaveBegin(tree);

  ReadUserObjects(tree);

  fdNdEtaCorrection = new AlidNdEtaCorrection("dndeta_correction", "dndeta_correction");

  fPIDParticles = new TH1F("pid_particles", "PID of generated primary particles", 10001, -5000.5, 5000.5);
  fPIDTracks = new TH1F("pid_tracks", "MC PID of reconstructed tracks", 10001, -5000.5, 5000.5);

  fClustersITSPos = new TH1F("clusters_its_pos", "clusters_its_pos", 7, -0.5, 6.5);
  fClustersTPCPos = new TH1F("clusters_tpc_pos", "clusters_tpc_pos", 160, -0.5, 159.5);

  fClustersITSNeg = new TH1F("clusters_its_neg", "clusters_its_neg", 7, -0.5, 6.5);
  fClustersTPCNeg = new TH1F("clusters_tpc_neg", "clusters_tpc_neg", 160, -0.5, 159.5);
}

Bool_t AlidNdEtaCorrectionSelector::Process(Long64_t entry)
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
  //  Assuming that fTree is the pointer to the TChain being processed,
  //  use fTree->GetTree()->GetEntry(entry).

  AliDebug(AliLog::kDebug+1,"Processing event ...\n");


  if (AliSelectorRL::Process(entry) == kFALSE)
    return kFALSE;

  // check prerequesites
  if (!fESD)
  {
    AliDebug(AliLog::kError, "ESD branch not available");
    return kFALSE;
  }

  AliHeader* header = GetHeader();
  if (!header)
  {
    AliDebug(AliLog::kError, "Header not available");
    return kFALSE;
  }

  AliStack* stack = GetStack();
  if (!stack)
  {
    AliDebug(AliLog::kError, "Stack not available");
    return kFALSE;
  }

  if (!fEsdTrackCuts)
  {
    AliDebug(AliLog::kError, "fESDTrackCuts not available");
    return kFALSE;
  }

  Bool_t vertexReconstructed = AliPWG0Helper::IsVertexReconstructed(fESD);

  Bool_t eventTriggered = AliPWG0Helper::IsEventTriggered(fESD);

  // get the MC vertex
  AliGenEventHeader* genHeader = header->GenEventHeader();

  TArrayF vtxMC(3);
  genHeader->PrimaryVertex(vtxMC);

  // loop over mc particles
  Int_t nPrim  = stack->GetNprimary();

  for (Int_t iMc = 0; iMc < nPrim; ++iMc)
  {
    AliDebug(AliLog::kDebug+1, Form("MC Loop: Processing particle %d.", iMc));

    TParticle* particle = stack->Particle(iMc);

    if (!particle)
    {
      AliDebug(AliLog::kError, Form("UNEXPECTED: particle with label %d not found in stack (mc loop).", iMc));
      continue;
    }

    if (AliPWG0Helper::IsPrimaryCharged(particle, nPrim) == kFALSE)
      continue;

    if (SignOK(particle->GetPDG()) == kFALSE)
      continue;

    Float_t eta = particle->Eta();
    Float_t pt = particle->Pt();

    if (vertexReconstructed) {
      fdNdEtaCorrection->FillParticle(vtxMC[2], eta, pt);

      if (pt > 0.1 && pt < 0.2)
	fPIDParticles->Fill(particle->GetPdgCode());
    }
  }// end of mc particle

  // ########################################################
  // loop over esd tracks
  Int_t nTracks = fESD->GetNumberOfTracks();

  // count the number of "good" tracks as parameter for vertex reconstruction efficiency
  Int_t nGoodTracks = 0;
  for (Int_t t=0; t<nTracks; t++)
  {
    AliDebug(AliLog::kDebug+1, Form("ESD Loop: Processing track %d.", t));

    AliESDtrack* esdTrack = fESD->GetTrack(t);

    // cut the esd track?
    if (!fEsdTrackCuts->AcceptTrack(esdTrack))
      continue;

    nGoodTracks++;

    // using the properties of the mc particle
    Int_t label = TMath::Abs(esdTrack->GetLabel());
    if (label == 0)
    {
      AliDebug(AliLog::kWarning, Form("WARNING: cannot find corresponding mc part for track %d.", t));
      continue;
    }

    TParticle* particle = stack->Particle(label);
    if (!particle)
    {
      AliDebug(AliLog::kError, Form("UNEXPECTED: particle with label %d not found in stack (track loop).", label));
      continue;
    }

    if (SignOK(particle->GetPDG()) == kFALSE)
        continue;

    if (vertexReconstructed)
    {
      fdNdEtaCorrection->FillParticleWhenMeasuredTrack(vtxMC[2], particle->Eta(), particle->Pt());
      if (particle->Pt() > 0.1 && particle->Pt() < 0.2)
	{
	  fPIDTracks->Fill(particle->GetPdgCode());
        if (particle->GetPDG()->Charge() > 0)
        {
          fClustersITSPos->Fill(esdTrack->GetITSclusters(0));
          fClustersTPCPos->Fill(esdTrack->GetTPCclusters(0));
        }
        else
        {
          fClustersITSNeg->Fill(esdTrack->GetITSclusters(0));
          fClustersTPCNeg->Fill(esdTrack->GetTPCclusters(0));
        }
      }
    }
  } // end of track loop

  // stuff regarding the vertex reco correction and trigger bias correction
  if (eventTriggered) {
    fdNdEtaCorrection->FillEventWithTrigger(vtxMC[2], nGoodTracks);
    if (vertexReconstructed)
      fdNdEtaCorrection->FillEventWithTriggerWithReconstructedVertex(vtxMC[2], nGoodTracks);
  }

  // getting process information
  Int_t processtype = AliPWG0Helper::GetPythiaEventProcessType(header);
  AliDebug(AliLog::kDebug+1,Form(" Found pythia procces type %d", processtype));

  if (processtype<0)
    AliDebug(AliLog::kError, Form("Unkown Pythia process type %d.", processtype));
  
  fdNdEtaCorrection->FillEventAll(vtxMC[2], nGoodTracks, "INEL");
  
  if (processtype!=92 && processtype!=93)
    fdNdEtaCorrection->FillEventAll(vtxMC[2], nGoodTracks, "NSD");
  
  if (processtype!=92 && processtype!=93 && processtype!=94)
    fdNdEtaCorrection->FillEventAll(vtxMC[2], nGoodTracks, "ND");

  return kTRUE;
}

void AlidNdEtaCorrectionSelector::SlaveTerminate()
{
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.

  AliSelectorRL::SlaveTerminate();

  // Add the histograms to the output on each slave server
  if (!fOutput)
  {
    AliDebug(AliLog::kError, "ERROR: Output list not initialized");
    return;
  }

  fOutput->Add(fdNdEtaCorrection);
}

void AlidNdEtaCorrectionSelector::Terminate()
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

  AliSelectorRL::Terminate();

  fdNdEtaCorrection = dynamic_cast<AlidNdEtaCorrection*> (fOutput->FindObject("dndeta_correction"));
  if (!fdNdEtaCorrection)
  {
    AliDebug(AliLog::kError, "Could not read object from output list");
    return;
  }

  fdNdEtaCorrection->Finish();

  TFile* fout = new TFile(Form("correction_map%s.root", GetOption()), "RECREATE");

  if (fEsdTrackCuts)
    fEsdTrackCuts->SaveHistograms("esd_track_cuts");
  fdNdEtaCorrection->SaveHistograms();

  fout->Write();
  fout->Close();

  fdNdEtaCorrection->DrawHistograms();

  if (fPIDParticles && fPIDTracks)
  {
    new TCanvas("pidcanvas", "pidcanvas", 500, 500);

    fPIDParticles->Draw();
    fPIDTracks->SetLineColor(2);
    fPIDTracks->Draw("SAME");

    TDatabasePDG* pdgDB = new TDatabasePDG;

    for (Int_t i=0; i <= fPIDParticles->GetNbinsX()+1; ++i)
      if (fPIDParticles->GetBinContent(i) > 0)
        printf("PDG = %d (%s): generated: %d, reconstructed: %d, ratio: %f\n", (Int_t) fPIDParticles->GetBinCenter(i), pdgDB->GetParticle((Int_t) fPIDParticles->GetBinCenter(i))->GetName(), (Int_t) fPIDParticles->GetBinContent(i), (Int_t) fPIDTracks->GetBinContent(i), ((fPIDTracks->GetBinContent(i) > 0) ? fPIDParticles->GetBinContent(i) / fPIDTracks->GetBinContent(i) : -1));

    delete pdgDB;
    pdgDB = 0;
  }

  if (fClustersITSPos && fClustersITSNeg && fClustersTPCPos && fClustersTPCNeg)
  {
    TCanvas* canvas = new TCanvas("clusters", "clusters", 1000, 500);
    canvas->Divide(2, 1);

    canvas->cd(1);
    fClustersITSPos->Draw();
    fClustersITSNeg->SetLineColor(kRed);
    fClustersITSNeg->Draw("SAME");

    canvas->cd(2);
    fClustersTPCPos->Draw();
    fClustersTPCNeg->SetLineColor(kRed);
    fClustersTPCNeg->Draw("SAME");
  }
}
