#include "AliAnaTaskV0EffDecomposition.h"

// ROOT includes
#include <TList.h>
#include <TH1.h>


// AliRoot includes
#include <AliMCEventHandler.h>
#include <AliMCEvent.h>
#include <AliAnalysisManager.h>


// STL includes
#include <iostream>
using namespace std;


ClassImp(AliAnaTaskV0EffDecomposition)

//_____________________________________________________________________________
AliAnaTaskV0EffDecomposition::AliAnaTaskV0EffDecomposition():
  AliAnalysisTaskSE(),
  fAOD(0x0),
  fMCArray(0x0),
  fTrackFilterBit(1), // default is TPC filter bit
  fTrigBit(AliVEvent::kMB), // default is MB
  fVtxCut(10),        // default is 10 cm
  fEtaCut(0.8),       // default is 0.8
  fMinCent(0),        // default is 0%
  fMaxCent(10),       // default is 10%
  fPdgV0(3122),       // default is Lambda
  fPdgPos(2212),      // default is p
  fPdgNeg(-211),      // default is pi-
  fListOfObjects(0x0),
  hV0Gen(0x0),
  hV0Rec(0x0),
  hDaughterRec(0x0)
{
  // Default constructor (should not be used)
}

//______________________________________________________________________________
AliAnaTaskV0EffDecomposition::AliAnaTaskV0EffDecomposition(const char *name):
  AliAnalysisTaskSE(name),
  fAOD(0x0),
  fMCArray(0x0),
  fTrackFilterBit(1), // default is TPC filter bit
  fTrigBit(AliVEvent::kMB), // default is MB
  fVtxCut(10),        // default is 10 cm
  fEtaCut(0.8),       // default is 0.8
  fMinCent(0),        // default is 0%
  fMaxCent(10),       // default is 10%
  fPdgV0(3122),       // default is Lambda
  fPdgPos(2212),      // default is p
  fPdgNeg(-211),      // default is pi-
  fListOfObjects(0x0),
  hV0Gen(0x0),
  hV0Rec(0x0),
  hDaughterRec(0x0)
{
  // Output slot #1 writes into a TList
  DefineOutput(1, TList::Class());
}

//_____________________________________________________________________________
AliAnaTaskV0EffDecomposition::~AliAnaTaskV0EffDecomposition()
{
  // Destructor
  // histograms are in the output list and deleted when the output
  // list is deleted by the TSelector dtor
}

//______________________________________________________________________________
void AliAnaTaskV0EffDecomposition::UserCreateOutputObjects()
{ 
  // This method is called once per worker node
  // Here we define the output: histograms and debug tree if requested 

  OpenFile(1);
  fListOfObjects = new TList();
  fListOfObjects->SetOwner();
  
  //
  // Histograms
  //  
  hV0Gen = new TH1D("hV0Gen", "Number of generated V0s; p_{T} [GeV/c]; Counts", 
		    40, 0, 20);
  fListOfObjects->Add(hV0Gen);

  hV0Rec = new TH1D("hV0Rec", "Number of reconstructed V0s; p_{T} [GeV/c]; Counts", 
		    40, 0, 20);
  fListOfObjects->Add(hV0Rec);

  hDaughterRec = new TH1D("hDaughterRec", "Number of reconstructed track pairs; V0 p_{T} [GeV/c]; Counts", 
		       40, 0, 20);
  fListOfObjects->Add(hDaughterRec);

  hV0Ghost = new TH1D("hV0Ghost", "Number of reconstructed V0 ghosts; V0 p_{T} [GeV/c]; Counts", 
		      40, 0, 20);
  fListOfObjects->Add(hV0Ghost);

  hTrackGhost = new TH1D("hTrackGhost", "Number of reconstructed track pair ghosts; V0 p_{T} [GeV/c]; Counts", 
			 40, 0, 20);
  fListOfObjects->Add(hTrackGhost);

  hV0ButNoTracks = new TH1D("hV0ButNoTracks", "Strange V0s: V0 rec but no daughter tracks; V0 p_{T} [GeV/c]; Counts", 
			 40, 0, 20);
  fListOfObjects->Add(hV0ButNoTracks);
  
  // One could add much more differential histograms if needed

  // Post output data.
  PostData(1, fListOfObjects);
}

//______________________________________________________________________________
void AliAnaTaskV0EffDecomposition::UserExec(Option_t *) 
{
  // Main loop
  
  //
  // Comment: This method matches completely the same method for the high pT
  // tracks
  //


  //
  // First we make sure that we have valid input(s)!
  //
  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  if(!fAOD){
    Printf("%s:%d AODEvent not found in Input Manager",(char*)__FILE__,__LINE__);
    this->Dump();
    return;
  }

  // fMC = dynamic_cast<AliMCEvent*>(MCEvent());
  // if(fMC)
  //   fMC->Dump();
  
  fMCArray = (TClonesArray*)fAOD->FindListObject("mcparticles");
  if(!fMCArray){
    Printf("%s:%d AOD MC array not found in Input Manager",(char*)__FILE__,__LINE__);
    this->Dump();
    return;
  }    

  // Here we check: 
  // 1) that the event was triggered
  if(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))
     ->IsEventSelected() & fTrigBit ) {

    // 2) if the reconstrcuted vertex was in the correct range
    Float_t zvtx = -999;
    const AliVVertex* primaryVertex = fAOD->GetPrimaryVertex(); 
    if(primaryVertex->GetNContributors()>0)
      zvtx = primaryVertex->GetZ();
    
    if (TMath::Abs(zvtx) < fVtxCut) {	
      
      
      // 3) if the centrality is in the correct range
      Float_t centrality = -10;
      AliCentrality *centObject = fAOD->GetCentrality();
      if(centObject)
	centrality = centObject->GetCentralityPercentile("V0M"); 
      if((centrality < fMaxCent) && (centrality>=fMinCent)) {
	
	// as we focus here on Pb-Pb data I decide that the simplest would be
	// to demand that all 3 are ok before going on
	
	// Fill MC gen histograms
	ProcessMCTruthAOD();
	
	// Fill V0 histograms (loop over V0s)
	AnalyzeV0AOD();
	
	// Fill track pair histograms (loop over tracks)
	AnalyzeDaughtersAOD();	

	// In this final step we try to see for each generated V0 how many
	// times it was reconstructed as a V0 and/or a track pair
	AnalyzeRecMothersAOD();
      }
    }
  }

  // Post output data.
  PostData(1, fListOfObjects);
}

//_____________________________________________________________________________
void AliAnaTaskV0EffDecomposition::ProcessMCTruthAOD() 
{
  // Fill the special MC histogram with the MC truth info

  const Int_t nTracksMC = fMCArray->GetEntriesFast();
  
  for (Int_t iTracks = 0; iTracks < nTracksMC; iTracks++) {
    
    AliAODMCParticle* trackMC = dynamic_cast<AliAODMCParticle*>(fMCArray->At(iTracks));
    // We will use the Generator Index as a counter to see what we reconstruct
    trackMC->SetGeneratorIndex(0);


    Int_t pdgCode = trackMC->PdgCode();
    if(pdgCode != fPdgV0)
      continue;
    
    // Select only primaries
    if(!(trackMC->IsPhysicalPrimary()))
      continue;
    
    if (TMath::Abs(trackMC->Eta()) > fEtaCut )
      continue;
    
    Float_t ptMC      = trackMC->Pt();
    // Float_t pMC       = trackMC->P();
    // Float_t etaMC     = trackMC->Eta();
    // Float_t phiMC     = trackMC->Phi();
    
    hV0Gen->Fill(ptMC);
  }
}




//____________________________________________________________________
AliAODMCParticle* AliAnaTaskV0EffDecomposition::FindPrimaryMotherAOD(AliAODMCParticle* startParticle, Int_t& nSteps)
{
  //
  // Finds the first mother among the primary particles of the particle
  // identified by <label>, i.e. the primary that "caused" this particle
  //
  // Taken from AliPWG0Helper class
  //

  nSteps = 0;

  AliAODMCParticle* mcPart = startParticle;

  while (mcPart) {
    
    if(mcPart->IsPrimary())
      return mcPart;
    
    Int_t mother = mcPart->GetMother();
    
    mcPart = dynamic_cast<AliAODMCParticle*>(fMCArray->At(mother));
    nSteps++; // 1 level down
  }
  
  return 0;
}

//_____________________________________________________________________
void AliAnaTaskV0EffDecomposition::AnalyzeV0AOD() {
  
  Int_t nV0s = fAOD->GetNumberOfV0s();
  if(nV0s < 1)
    return;
  
  // Check that a primary vertex was reconstructed
  AliAODVertex *myBestPrimaryVertex = fAOD->GetPrimaryVertex();
  if (!myBestPrimaryVertex) return;
  
  
  for (Int_t iV0 = 0; iV0 < nV0s; iV0++) {
    
    AliAODv0 *aodV0 = fAOD->GetV0(iV0);
    if (!aodV0) continue;
    
    // common part
    
    // AliAODTrack (V0 Daughters)
    AliAODVertex* vertex = aodV0->GetSecondaryVtx();
    if (!vertex) {
      Printf("ERROR: Could not retrieve vertex");
      continue;
    }
    
    AliAODTrack *p_track = (AliAODTrack*)vertex->GetDaughter(0);
    AliAODTrack *n_track = (AliAODTrack*)vertex->GetDaughter(1);
    if (!p_track || !n_track) {
      Printf("ERROR: Could not retrieve one of the daughter track");
      continue;
    }
    
    // Remove like-sign
    if (p_track->Charge() == n_track->Charge()) {
      continue;
    } 
    
    // Make sure charge ordering is ok
    if (p_track->Charge() < 0) {
      AliAODTrack* helpTrack = p_track;
      p_track = n_track;
      n_track = helpTrack;
    } 
    
    // // Eta cut on decay products ?
    // if(TMath::Abs(p_track->Eta()) > fEtaCut || TMath::Abs(n_track->Eta()) > fEtaCut)
    //   continue;
    
    AliAODMCParticle* p_mother = ValidateTrack(p_track, fPdgPos);
    if(!p_mother)
      continue;
    AliAODMCParticle* n_mother = ValidateTrack(n_track, fPdgNeg);
    if(!n_mother)
      continue;
    
    // check that mother is the same
    if(p_mother != n_mother)
      continue;

    // check that mother has good eta
    if (TMath::Abs(p_mother->Eta()) > fEtaCut )
      continue;
    
    // One could also consider to fill the pT of the reconstructed mother
    hV0Rec->Fill(p_mother->Pt());
    p_mother->SetGeneratorIndex(p_mother->GetGeneratorIndex() + 1);
  }//end loop over v0's
}


//________________________________________________________________________
AliAODMCParticle* AliAnaTaskV0EffDecomposition::ValidateTrack(AliAODTrack* track, 
							   Int_t pdgDaughter)
{
  // Validate V0 daughter track
  
  // Apply track quality cuts
  if(!track->TestFilterBit(fTrackFilterBit)) {
    return 0;
  }
  
  // Do we want to check the invarinat mass here? e.g. a la
  // Double_t deltaInvMassK0s   = aodV0->MassK0Short()-0.498;
  
  const Int_t label = TMath::Abs(track->GetLabel());
  AliAODMCParticle* mcTrack = dynamic_cast<AliAODMCParticle*>(fMCArray->At(label));
  if (!mcTrack)
    return 0;
  if(mcTrack->IsPhysicalPrimary())
    return 0;
  Int_t pdgCode = mcTrack->GetPdgCode();
  if(pdgCode != pdgDaughter)
    return 0;
  
  // mother_steps is the number of steps we have to go backe to find the
  // primary mother. Here we only accept primary V0s so it has to be 1
  Int_t mother_steps = 0;
  AliAODMCParticle* mother = FindPrimaryMotherAOD(mcTrack, mother_steps);
  if(!mother) 
    return 0;
  if(mother_steps != 1)
    return 0;
  Int_t pdgMother = mother->GetPdgCode();
  if(pdgMother != fPdgV0)
    return 0;
  
  return mother;
}



//________________________________________________________________________
void AliAnaTaskV0EffDecomposition::AnalyzeDaughtersAOD()
{
  // check to see if two daughter tracks of the same mother was reconstructed
  // (here there is no V0 requirement)
  const Int_t nAODTracks = fAOD->GetNumberOfTracks();
  // Bool_t wasUsed[nAODTracks];
  // for(Int_t iT = 0; iT < nAODTracks; iT++)
  //   wasUsed[iT] = kFALSE;
  
  for(Int_t iT = 0; iT < nAODTracks; iT++) {
    
    // if(wasUsed[iT] == kTRUE)
    //   continue; // already part of a matched pair
    AliAODTrack* track1 = (AliAODTrack*)fAOD->GetTrack(iT);
    
    Int_t charge  = track1->Charge();

    Int_t pdgDaughter1 = fPdgPos;
    Int_t pdgDaughter2 = fPdgNeg;

    if(charge < 0) {

       pdgDaughter1 = fPdgNeg;
       pdgDaughter2 = fPdgPos;
    }
    
    AliAODMCParticle* mother1 = ValidateTrack(track1, pdgDaughter1);
    if(!mother1)
      continue;

    // loop over the remaining tracks and see of the other daughter track was
    // also reconstructed
    for(Int_t jT = iT+1; jT < nAODTracks; jT++) {
      
      // if(wasUsed[jT] == kTRUE)
      // 	continue; // already part of a matched pair
      AliAODTrack* track2 = (AliAODTrack*)fAOD->GetTrack(jT);
    
      AliAODMCParticle* mother2 = ValidateTrack(track2, pdgDaughter2);
      if(!mother2)
	continue;
      
      // check that mother is the same
      if(mother1 != mother2)
	continue;
      
      //      wasUsed[jT] = kTRUE;
      
      // check that mother has good eta
      if (TMath::Abs(mother1->Eta()) < fEtaCut ) {
	
	hDaughterRec->Fill(mother1->Pt());
	mother1->SetGeneratorIndex(mother1->GetGeneratorIndex() + 100);

      }
    }    
  }
}
  
//_____________________________________________________________________________
void AliAnaTaskV0EffDecomposition::AnalyzeRecMothersAOD() 
{
  // See how many mothers were reconstruted and how many times

  const Int_t nTracksMC = fMCArray->GetEntriesFast();
  
  for (Int_t iTracks = 0; iTracks < nTracksMC; iTracks++) {
    
    AliAODMCParticle* trackMC = dynamic_cast<AliAODMCParticle*>(fMCArray->At(iTracks));

    Int_t pdgCode = trackMC->PdgCode();
    if(pdgCode != fPdgV0)
      continue;
    
    // Select only primaries
    if(!(trackMC->IsPhysicalPrimary()))
      continue;
    
    if (TMath::Abs(trackMC->Eta()) > fEtaCut )
      continue;
    
    Float_t ptMC      = trackMC->Pt();
    // Float_t pMC       = trackMC->P();
    // Float_t etaMC     = trackMC->Eta();
    // Float_t phiMC     = trackMC->Phi();
    
    Int_t nV0rec = trackMC->GetGeneratorIndex()%100;
    Int_t nTrackrec = Int_t(trackMC->GetGeneratorIndex()/100);

    if(nV0rec > 1)
      hV0Ghost->Fill(ptMC);
    if(nTrackrec > 1)
      hTrackGhost->Fill(ptMC);
    if(nV0rec==1 && nTrackrec==0)
      hV0ButNoTracks->Fill(ptMC);      
  }
}
