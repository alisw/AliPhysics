/* $Id$ */

#include "AlidNdEtaSystematicsSelector.h"

#include <TStyle.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TVector3.h>
#include <TChain.h>
#include <TFile.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TParticle.h>
#include <TParticlePDG.h>

#include <AliLog.h>
#include <AliESD.h>
#include <AliStack.h>
#include <AliHeader.h>
#include <AliGenEventHeader.h>
#include <../STEER/AliGenPythiaEventHeader.h>
#include <../STEER/AliGenCocktailEventHeader.h>

#include "AliESDtrackCuts.h"
#include "AliPWG0Helper.h"
#include "dNdEta/AlidNdEtaCorrection.h"

ClassImp(AlidNdEtaSystematicsSelector)

AlidNdEtaSystematicsSelector::AlidNdEtaSystematicsSelector() :
  AliSelectorRL(),
  fSecondaries(0),
  fSigmaVertex(0),
  fEsdTrackCuts(0),
  fPIDParticles(0),
  fPIDTracks(0),
  fSignMode(0),
  fMultiplicityMode(0)
{
  //
  // Constructor. Initialization of pointers
  //

  for (Int_t i=0; i<4; ++i)
    fdNdEtaCorrectionSpecies[i] = 0;

  for (Int_t i=0; i<3; ++i) {
    fdNdEtaCorrectionVertexReco[i] = 0;
    fdNdEtaCorrectionTriggerBias[i] = 0;
  }
}

AlidNdEtaSystematicsSelector::~AlidNdEtaSystematicsSelector()
{
  //
  // Destructor
  //
}

void AlidNdEtaSystematicsSelector::Begin(TTree* tree)
{
  // Begin function

  AliSelectorRL::Begin(tree);

  ReadUserObjects(tree);
}

void AlidNdEtaSystematicsSelector::ReadUserObjects(TTree* tree)
{
  // read the user objects, called from slavebegin and begin

  if (!fEsdTrackCuts && fInput)
    fEsdTrackCuts = dynamic_cast<AliESDtrackCuts*> (fInput->FindObject("AliESDtrackCuts"));

  if (!fEsdTrackCuts && tree)
    fEsdTrackCuts = dynamic_cast<AliESDtrackCuts*> (tree->GetUserInfo()->FindObject("AliESDtrackCuts"));

  if (!fEsdTrackCuts)
     AliDebug(AliLog::kError, "ERROR: Could not read EsdTrackCuts from input list.");
}

void AlidNdEtaSystematicsSelector::SlaveBegin(TTree* tree)
{
  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).

  AliSelector::SlaveBegin(tree);

  ReadUserObjects(tree);

  TString option(GetOption());

  printf("Running AlidNdEtaSystematicsSelector with options %s\n", option.Data());

  // Options: secondaries particle-composition sigma-vertex vertexreco triggerbias

  if (option.Contains("secondaries"))
  {
    fSecondaries = new TH2F("fSecondaries", "fSecondaries;Case;N", 8, 0.5, 8.5, 2001, -0.5, 2000.5);
    fSecondaries->GetXaxis()->SetBinLabel(1, "All Primaries");
    fSecondaries->GetXaxis()->SetBinLabel(2, "All Secondaries");
    fSecondaries->GetXaxis()->SetBinLabel(3, "Primaries pT > 0.3 GeV/c");
    fSecondaries->GetXaxis()->SetBinLabel(4, "Secondaries pT > 0.3 GeV/c");
    fSecondaries->GetXaxis()->SetBinLabel(5, "Tracks from Primaries");
    fSecondaries->GetXaxis()->SetBinLabel(6, "Tracks from Secondaries");
    fSecondaries->GetXaxis()->SetBinLabel(7, "Accepted Tracks from Primaries");
    fSecondaries->GetXaxis()->SetBinLabel(8, "Accepted Tracks from Secondaries");
  }

  if (option.Contains("particle-composition"))
  {
    for (Int_t i=0; i<4; ++i)
    {
      TString name;
      name.Form("correction_%d", i);
      fdNdEtaCorrectionSpecies[i] = new AlidNdEtaCorrection(name, name);
    }
  }

  if (option.Contains("sigma-vertex"))
  {
    fSigmaVertex = new TH1F("fSigmaVertex", "fSigmaVertex;Nsigma2vertex;NacceptedTracks", 50, 0.05, 5.05);
    printf("WARNING: sigma-vertex analysis enabled. This will produce weird results in the AliESDtrackCuts histograms\n");
  }

  if (option.Contains("vertexreco")) {
    fdNdEtaCorrectionVertexReco[0] = new AlidNdEtaCorrection("vertexRecoND", "vertexRecoND");
    fdNdEtaCorrectionVertexReco[1] = new AlidNdEtaCorrection("vertexRecoSD", "vertexRecoSD");
    fdNdEtaCorrectionVertexReco[2] = new AlidNdEtaCorrection("vertexRecoDD", "vertexRecoDD");
  }

  if (option.Contains("triggerbias")) {
    fdNdEtaCorrectionTriggerBias[0] = new AlidNdEtaCorrection("triggerBiasND", "triggerBiasND");
    fdNdEtaCorrectionTriggerBias[1] = new AlidNdEtaCorrection("triggerBiasSD", "triggerBiasSD");
    fdNdEtaCorrectionTriggerBias[2] = new AlidNdEtaCorrection("triggerBiasDD", "triggerBiasDD");
  }

  fPIDParticles = new TH1F("pid_particles", "PID of generated primary particles", 10001, -5000.5, 5000.5);
  fPIDTracks = new TH1F("pid_tracks", "MC PID of reconstructed tracks", 10001, -5000.5, 5000.5);

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

  if (option.Contains("low-multiplicity"))
  {
    AliInfo("Processing only events with low multiplicity.");
    fMultiplicityMode = 1;
  }
  else if (option.Contains("high-multiplicity"))
  {
    AliInfo("Processing only events with high multiplicity.");
    fMultiplicityMode = 2;
  }
}

Bool_t AlidNdEtaSystematicsSelector::Process(Long64_t entry)
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

  if (AliSelectorRL::Process(entry) == kFALSE)
    return kFALSE;

  // Check prerequisites
  if (!fESD)
  {
    AliDebug(AliLog::kError, "ESD branch not available");
    return kFALSE;
  }

  if (!fEsdTrackCuts)
  {
    AliDebug(AliLog::kError, "fESDTrackCuts not available");
    return kFALSE;
  }

  AliStack* stack = GetStack();
  if (!stack)
  {
    AliDebug(AliLog::kError, "Stack not available");
    return kFALSE;
  }

  // --------------------------------------------------------------
  // related to events:
  //
  // stuff for vertex reconstruction correction systematics
  Bool_t vertexRecoStudy  = kFALSE;
  Bool_t triggerBiasStudy = kFALSE;
  if (fdNdEtaCorrectionVertexReco[0])  vertexRecoStudy = kTRUE;
  if (fdNdEtaCorrectionTriggerBias[0]) triggerBiasStudy = kTRUE;

  AliHeader* header = GetHeader();
  if (!header) {
    AliDebug(AliLog::kError, "Header not available");
    return kFALSE;
  }
  
  // getting process information
  Int_t processtype = AliPWG0Helper::GetEventProcessType(header);

  // can only read pythia headers, either directly or from cocktalil header
  AliGenEventHeader* genHeader = (AliGenEventHeader*)(header->GenEventHeader());
  
  if (!genHeader) {
    AliDebug(AliLog::kError, "Gen header not available");
    return kFALSE;
  }
  
  // get the MC vertex
  TArrayF vtxMC(3);
  genHeader->PrimaryVertex(vtxMC);
  
  TObjArray* listOfTracks    = fEsdTrackCuts->GetAcceptedTracks(fESD);
  Int_t nGoodTracks          = listOfTracks->GetEntries();
  
  Bool_t eventTriggered      = AliPWG0Helper::IsEventTriggered(fESD);
  Bool_t vertexReconstructed = AliPWG0Helper::IsVertexReconstructed(fESD);
  
  // non diffractive
  if (processtype == AliPWG0Helper::kND) {
    // NB: passing the wrong process type here (1), since the process type is defined by the index in the array (here non-diffractive) // ???? CKB ???
    if (triggerBiasStudy) fdNdEtaCorrectionTriggerBias[0]->FillEvent(vtxMC[2], nGoodTracks, eventTriggered, vertexReconstructed, 1);      
    if (vertexRecoStudy)  fdNdEtaCorrectionVertexReco[0] ->FillEvent(vtxMC[2], nGoodTracks, eventTriggered, vertexReconstructed, 1);
  }
  
  // single diffractive
  if (processtype == AliPWG0Helper::kSD) { 
    if (triggerBiasStudy) fdNdEtaCorrectionTriggerBias[1]->FillEvent(vtxMC[2], nGoodTracks, eventTriggered, vertexReconstructed, 1);     
    if (vertexRecoStudy)  fdNdEtaCorrectionVertexReco[1] ->FillEvent(vtxMC[2], nGoodTracks, eventTriggered, vertexReconstructed, 1);
  }
  
  // double diffractive
  if (processtype==AliPWG0Helper::kDD) { 
    if (triggerBiasStudy) fdNdEtaCorrectionTriggerBias[2]->FillEvent(vtxMC[2], nGoodTracks, eventTriggered, vertexReconstructed, 1);     
    if (vertexRecoStudy)  fdNdEtaCorrectionVertexReco[2] ->FillEvent(vtxMC[2], nGoodTracks, eventTriggered, vertexReconstructed, 1);
  }

  for (Int_t i=0; i<4; ++i) {
    if (fdNdEtaCorrectionSpecies[i])
      fdNdEtaCorrectionSpecies[i]->FillEvent(vtxMC[2], nGoodTracks, eventTriggered, vertexReconstructed, 1);
  }

  // --------------------------------------------------------------
  // MC particle loop
  //

  Int_t nPrim  = stack->GetNprimary();

  for (Int_t iMc = 0; iMc < nPrim; ++iMc)
  {
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
    
    Int_t id = -1;
    switch (TMath::Abs(particle->GetPdgCode()))
      {
      case 211: id = 0; break;
      case 321: id = 1; break;
      case 2212: id = 2; break;
      case 11:
        {
          /*if (pt < 0.1 && particle->GetMother(0) > -1)
          {
            TParticle* mother = stack->Particle(particle->GetMother(0));
            printf("Mother pdg code is %d\n", mother->GetPdgCode());
          } */
          //particle->Dump();
          //if (particle->GetUniqueID() == 1)

          //printf("Found illegal particle mc id = %d file = %s event = %d\n", iMc,  fTree->GetCurrentFile()->GetName(), fTree->GetTree()->GetReadEntry());
        }
      default: id = 3; break;
    }

    if (eventTriggered & vertexReconstructed) {
      if (fdNdEtaCorrectionSpecies[id])
	fdNdEtaCorrectionSpecies[id]->FillMCParticle(vtxMC[2], eta, pt, eventTriggered, vertexReconstructed, 1);
      //if (pt < 0.1)
      if (fPIDParticles)
	fPIDParticles->Fill(particle->GetPdgCode());
    }
    
    // non diffractive
    if (processtype==AliPWG0Helper::kND) { 
      if (triggerBiasStudy) fdNdEtaCorrectionTriggerBias[0]->FillMCParticle(vtxMC[2], eta, pt, eventTriggered, vertexReconstructed, 1);
      if (vertexRecoStudy)  fdNdEtaCorrectionVertexReco[0] ->FillMCParticle(vtxMC[2], eta, pt, eventTriggered, vertexReconstructed, 1);
    }
    // single diffractive
    if (processtype==AliPWG0Helper::kSD) { 
      if (triggerBiasStudy) fdNdEtaCorrectionTriggerBias[1]->FillMCParticle(vtxMC[2], eta, pt, eventTriggered, vertexReconstructed, 1);
      if (vertexRecoStudy)  fdNdEtaCorrectionVertexReco[1] ->FillMCParticle(vtxMC[2], eta, pt, eventTriggered, vertexReconstructed, 1);
    }
    // double diffractive
    if (processtype==AliPWG0Helper::kDD) { 
      if (triggerBiasStudy) fdNdEtaCorrectionTriggerBias[2]->FillMCParticle(vtxMC[2], eta, pt, eventTriggered, vertexReconstructed, 1);
      if (vertexRecoStudy)  fdNdEtaCorrectionVertexReco[2] ->FillMCParticle(vtxMC[2], eta, pt, eventTriggered, vertexReconstructed, 1);
    }

  }// end of mc particle


  // --------------------------------------------------------------
  // ESD track loop
  //

  if (fMultiplicityMode == 1 && listOfTracks->GetEntries() > 20 ||
      fMultiplicityMode == 2 && listOfTracks->GetEntries() < 40)
  {
    delete listOfTracks;
    listOfTracks = 0;
    return kTRUE;
  }

  // loop over esd tracks
  TIterator* iter = listOfTracks->MakeIterator();
  TObject* obj = 0;
  while ((obj = iter->Next()) != 0)
  {
    AliESDtrack* esdTrack = dynamic_cast<AliESDtrack*> (obj);
    if (!esdTrack)
      continue;    

    // using the properties of the mc particle
    Int_t label = TMath::Abs(esdTrack->GetLabel());
    TParticle* particle = stack->Particle(label);
    if (!particle)
    {
      AliDebug(AliLog::kError, Form("UNEXPECTED: particle with label %d not found in stack (track loop).", label));
      continue;
    }

    Float_t eta = particle->Eta();
    Float_t pt  = particle->Pt();
    
    // non diffractive
    if (processtype==AliPWG0Helper::kND) { 
      if (triggerBiasStudy) fdNdEtaCorrectionTriggerBias[0]->FillTrackedParticle(vtxMC[2], eta, pt);
      if (vertexRecoStudy)  fdNdEtaCorrectionVertexReco[0] ->FillTrackedParticle(vtxMC[2], eta, pt);
    }
    
    // single diffractive
    if (processtype==AliPWG0Helper::kSD) { 
      if (triggerBiasStudy) fdNdEtaCorrectionTriggerBias[1]->FillTrackedParticle(vtxMC[2], eta, pt);
      if (vertexRecoStudy)  fdNdEtaCorrectionVertexReco[1] ->FillTrackedParticle(vtxMC[2], eta, pt);
    }

    // double diffractive
    if (processtype==AliPWG0Helper::kDD) { 
      if (triggerBiasStudy) fdNdEtaCorrectionTriggerBias[2]->FillTrackedParticle(vtxMC[2], eta, pt);
      if (vertexRecoStudy)  fdNdEtaCorrectionVertexReco[2] ->FillTrackedParticle(vtxMC[2], eta, pt);
    }

    // find primary particle that created this particle
    TParticle* mother = AliPWG0Helper::FindPrimaryMother(stack, label);
    if (!mother)
      continue;

    if (SignOK(mother->GetPDG()) == kFALSE)
      continue;

    //printf("We continue with particle %d (pdg %d)\n", label, mother->GetPdgCode());

    Int_t id = -1;
    switch (TMath::Abs(mother->GetPdgCode()))
    {
      case 211:  id = 0; break;
      case 321:  id = 1; break;
      case 2212: id = 2; break;
      default:   id = 3; break;
    }

    if (vertexReconstructed && eventTriggered) {
      if (fdNdEtaCorrectionSpecies[id])
	fdNdEtaCorrectionSpecies[id]->FillTrackedParticle(vtxMC[2], eta, pt);
      //if (pt < 0.1)
      if (fPIDTracks)
	fPIDTracks->Fill(particle->GetPdgCode());
    }
  } // end of track loop

  
  delete iter;
  iter = 0;

  if (fSecondaries)
    FillSecondaries();

  if (fSigmaVertex)
    FillSigmaVertex();


  delete listOfTracks;
  listOfTracks = 0;

  return kTRUE;
}

void AlidNdEtaSystematicsSelector::FillSecondaries()
{
  // fills the secondary histograms
  
  AliStack* stack = GetStack();

  Int_t particlesPrimaries = 0;
  Int_t particlesSecondaries = 0;
  Int_t particlesPrimariesPtCut = 0;
  Int_t particlesSecondariesPtCut = 0;

  for (Int_t iParticle = 0; iParticle < stack->GetNtrack(); iParticle++)
  {
    TParticle* particle = stack->Particle(iParticle);
    if (!particle)
    {
      AliDebug(AliLog::kError, Form("UNEXPECTED: particle with label %d not found in stack (particle loop).", iParticle));
      continue;
    }

    if (TMath::Abs(particle->Eta()) > 0.9)
      continue;

    Bool_t isPrimary = kFALSE;
    if (iParticle < stack->GetNprimary())
    {
      if (AliPWG0Helper::IsPrimaryCharged(particle, stack->GetNprimary()) == kFALSE)
        continue;

      isPrimary = kTRUE;
    }

    if (isPrimary)
      particlesPrimaries++;
    else
      particlesSecondaries++;

    if (particle->Pt() < 0.3)
      continue;

    if (isPrimary)
      particlesPrimariesPtCut++;
    else
      particlesSecondariesPtCut++;
  }

  fSecondaries->Fill(1, particlesPrimaries);
  fSecondaries->Fill(2, particlesSecondaries);
  fSecondaries->Fill(3, particlesPrimariesPtCut);
  fSecondaries->Fill(4, particlesSecondariesPtCut);

  Int_t allTracksPrimaries = 0;
  Int_t allTracksSecondaries = 0;
  Int_t acceptedTracksPrimaries = 0;
  Int_t acceptedTracksSecondaries = 0;

  for (Int_t iTrack = 0; iTrack < fESD->GetNumberOfTracks(); iTrack++)
  {
    AliESDtrack* esdTrack = fESD->GetTrack(iTrack);
    if (!esdTrack)
    {
      AliDebug(AliLog::kError, Form("UNEXPECTED: Could not get track %d.", iTrack));
      continue;
    }

    Int_t label = TMath::Abs(esdTrack->GetLabel());
    TParticle* particle = stack->Particle(label);
    if (!particle)
    {
      AliDebug(AliLog::kError, Form("UNEXPECTED: particle with label %d not found in stack (track loop).", label));
      continue;
    }

   if (label < stack->GetNprimary())
      allTracksPrimaries++;
    else
      allTracksSecondaries++;

    if (!fEsdTrackCuts->AcceptTrack(esdTrack))
      continue;

    if (label < stack->GetNprimary())
      acceptedTracksPrimaries++;
    else
      acceptedTracksSecondaries++;
  }

  fSecondaries->Fill(5, allTracksPrimaries);
  fSecondaries->Fill(6, allTracksSecondaries);
  fSecondaries->Fill(7, acceptedTracksPrimaries);
  fSecondaries->Fill(8, acceptedTracksSecondaries);

  //printf("P = %d, S = %d, P_pt = %d, S_pt = %d, T_P = %d, T_S = %d, T_P_acc = %d, T_S_acc = %d\n", particlesPrimaries, particlesSecondaries, particlesPrimariesPtCut, particlesSecondariesPtCut, allTracksPrimaries, allTracksSecondaries, acceptedTracksPrimaries, acceptedTracksSecondaries);
}

void AlidNdEtaSystematicsSelector::FillSigmaVertex()
{
  // fills the fSigmaVertex histogram

  // save the old value
  Float_t oldSigmaVertex = fEsdTrackCuts->GetMinNsigmaToVertex();

  // set to maximum
  fEsdTrackCuts->SetMinNsigmaToVertex(5);

  TObjArray* list = fEsdTrackCuts->GetAcceptedTracks(fESD);

  TIterator* iter = list->MakeIterator();
  TObject* obj = 0;
  while ((obj = iter->Next()) != 0)
  {
    AliESDtrack* esdTrack = dynamic_cast<AliESDtrack*> (obj);
    if (!esdTrack)
      continue;

    Float_t sigma2Vertex = fEsdTrackCuts->GetSigmaToVertex(esdTrack);

    for (Double_t nSigma = 0.1; nSigma < 5.05; nSigma += 0.1)
    {
      if (sigma2Vertex < nSigma)
        fSigmaVertex->Fill(nSigma);
    }
  }

  delete iter;
  iter = 0;

  delete list;
  list = 0;

  // set back the old value
  fEsdTrackCuts->SetMinNsigmaToVertex(oldSigmaVertex);
}

void AlidNdEtaSystematicsSelector::SlaveTerminate()
{
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.

  AliSelector::SlaveTerminate();

  // Add the histograms to the output on each slave server
  if (!fOutput)
  {
    AliDebug(AliLog::kError, Form("ERROR: Output list not initialized."));
    return;
  }

  if (fSecondaries)
    fOutput->Add(fSecondaries);

  for (Int_t i=0; i<4; ++i)
    if (fdNdEtaCorrectionSpecies[i])
      fOutput->Add(fdNdEtaCorrectionSpecies[i]);

  if (fSigmaVertex)
    fOutput->Add(fSigmaVertex);

  for (Int_t i=0; i<3; ++i) {
    if (fdNdEtaCorrectionVertexReco[i])
      fOutput->Add(fdNdEtaCorrectionVertexReco[i]);
    
    if (fdNdEtaCorrectionTriggerBias[i])
      fOutput->Add(fdNdEtaCorrectionTriggerBias[i]);
  }
  
}

void AlidNdEtaSystematicsSelector::Terminate()
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

  AliSelector::Terminate();

  fSecondaries = dynamic_cast<TH2F*> (fOutput->FindObject("fSecondaries"));
  for (Int_t i=0; i<4; ++i)
    fdNdEtaCorrectionSpecies[i] = dynamic_cast<AlidNdEtaCorrection*> (fOutput->FindObject(Form("correction_%d", i)));
  fSigmaVertex = dynamic_cast<TH1F*> (fOutput->FindObject("fSigmaVertex"));

  fdNdEtaCorrectionVertexReco[0] = dynamic_cast<AlidNdEtaCorrection*> (fOutput->FindObject("vertexRecoND"));
  fdNdEtaCorrectionVertexReco[1] = dynamic_cast<AlidNdEtaCorrection*> (fOutput->FindObject("vertexRecoSD"));
  fdNdEtaCorrectionVertexReco[2] = dynamic_cast<AlidNdEtaCorrection*> (fOutput->FindObject("vertexRecoDD"));

  fdNdEtaCorrectionTriggerBias[0] = dynamic_cast<AlidNdEtaCorrection*> (fOutput->FindObject("triggerBiasND"));
  fdNdEtaCorrectionTriggerBias[1] = dynamic_cast<AlidNdEtaCorrection*> (fOutput->FindObject("triggerBiasSD"));
  fdNdEtaCorrectionTriggerBias[2] = dynamic_cast<AlidNdEtaCorrection*> (fOutput->FindObject("triggerBiasDD"));

  if (fPIDParticles)
  {
    TDatabasePDG* pdgDB = new TDatabasePDG;

    for (Int_t i=0; i <= fPIDParticles->GetNbinsX()+1; ++i)
      if (fPIDParticles->GetBinContent(i) > 0)
        printf("PDG = %d (%s): generated: %d, reconstructed: %d, ratio: %f\n", 
	       (Int_t) fPIDParticles->GetBinCenter(i), 
	       pdgDB->GetParticle((Int_t) fPIDParticles->GetBinCenter(i))->GetName(), 
	       (Int_t) fPIDParticles->GetBinContent(i), 
	       (Int_t) fPIDTracks->GetBinContent(i), 
	       ((fPIDTracks->GetBinContent(i) > 0) ? fPIDParticles->GetBinContent(i) / fPIDTracks->GetBinContent(i) : -1));

    delete pdgDB;
    pdgDB = 0;
  }

  for (Int_t i=0; i<3; ++i) {
    if (fdNdEtaCorrectionVertexReco[i])
      fdNdEtaCorrectionVertexReco[i]->Finish();
    
    if (fdNdEtaCorrectionTriggerBias[i])
      fdNdEtaCorrectionTriggerBias[i]->Finish();
  }

  TFile* fout = TFile::Open("systematics.root", "RECREATE");

  if (fEsdTrackCuts)
    fEsdTrackCuts->SaveHistograms("esd_track_cuts");

  if (fSecondaries)
    fSecondaries->Write();

  if (fSigmaVertex)
    fSigmaVertex->Write();
  
  for (Int_t i=0; i<4; ++i)
    if (fdNdEtaCorrectionSpecies[i])
      fdNdEtaCorrectionSpecies[i]->SaveHistograms();

  for (Int_t i=0; i<3; ++i) {
    if (fdNdEtaCorrectionVertexReco[i])
      fdNdEtaCorrectionVertexReco[i]->SaveHistograms();

    if (fdNdEtaCorrectionTriggerBias[i])
      fdNdEtaCorrectionTriggerBias[i]->SaveHistograms();
  }

  fout->Write();
  fout->Close();

  if (fSecondaries)
  {
    new TCanvas;
    fSecondaries->Draw("COLZ");
  }

  if (fSigmaVertex)
  {
    new TCanvas;
    fSigmaVertex->Draw();
  }
}

Bool_t AlidNdEtaSystematicsSelector::SignOK(TParticlePDG* particle)
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
