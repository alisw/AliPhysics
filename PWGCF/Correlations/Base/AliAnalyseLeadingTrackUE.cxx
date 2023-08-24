/*************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: A.Abrahantes, E.Lopez, S.Vallero                               *
 * Version 1.0                                                            *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include <TList.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TObject.h>
#include <TVector3.h>

#include "AliAnalyseLeadingTrackUE.h"

#include "AliAODEvent.h"
#include "AliAODMCParticle.h"
#include "AliAODTrack.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliInputEventHandler.h"
#include "AliMCEvent.h"
#include "AliVParticle.h"
#include "AliAODMCHeader.h"
#include "TFormula.h"

#include "AliAnalysisManager.h"
#include "AliMCEventHandler.h"
#include "AliStack.h"
#include "AliPIDResponse.h"
#include "AliHelperPID.h"


////////////////////////////////////////////////
//--------------------------------------------- 
// Class for transverse regions analysis
//---------------------------------------------
////////////////////////////////////////////////


using namespace std;

ClassImp(AliAnalyseLeadingTrackUE)

//-------------------------------------------------------------------
AliAnalyseLeadingTrackUE::AliAnalyseLeadingTrackUE() :
  TObject(),
  fDebug(0),
  fFilterBit(16),
  fTrackStatus(0),
  fOnlyHadrons(kFALSE),
  fCheckMotherPDG(kTRUE),
  fTrackEtaCut(0.8),
  fTrackEtaCutMin(-1.),
  fTrackPhiCutEvPlMin(0.),
  fTrackPhiCutEvPlMax(0.),
  fTrackPtMin(0),
  fEventSelection(AliVEvent::kMB|AliVEvent::kUserDefined),
  fDCAXYCut(0),
  fDCAZCut(0),
  fSharedClusterCut(-1),
  fCrossedRowsCut(-1),
  fFoundFractionCut(-1),
  fEsdTrackCuts(0x0), 
  fEsdTrackCutsExtra1(0x0), 
  fEsdTrackCutsExtra2(0x0), 
  fHelperPID(0x0),
  fEventCounter(0)
{
  // constructor
}


//-------------------------------------------------------------------
AliAnalyseLeadingTrackUE & AliAnalyseLeadingTrackUE::operator = (const AliAnalyseLeadingTrackUE & /*source*/)
{
  // assignment operator
  return *this;
}


//-------------------------------------------------------------------
AliAnalyseLeadingTrackUE::~AliAnalyseLeadingTrackUE()
{

  //clear memory
  
}


//____________________________________________________________________
Bool_t AliAnalyseLeadingTrackUE::ApplyCuts(TObject* track)
{
  
  // select track according to set of cuts
  if (!fEsdTrackCuts->IsSelected(track) )return kFALSE;
  if (fEsdTrackCutsExtra1 && fEsdTrackCutsExtra2 && !fEsdTrackCutsExtra1->IsSelected(track) && !fEsdTrackCutsExtra2->IsSelected(track)) return kFALSE;

  return kTRUE;
}


//____________________________________________________________________
void AliAnalyseLeadingTrackUE::DefineESDCuts(Int_t filterbit) {
  
  // Reproduces the cuts of the corresponding bit in the ESD->AOD filtering
  // (see $ALICE_ROOT/ANALYSIS/macros/AddTaskESDFilter.C)
  
  if (filterbit == -1)
    filterbit = fFilterBit;

  if (filterbit == 128)
  {
    fEsdTrackCuts = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
    fEsdTrackCuts->SetMinNClustersTPC(70);
  }
  else if (filterbit == 256)
  {
    // syst study
    fEsdTrackCuts = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
    fEsdTrackCuts->SetMinNClustersTPC(80);
    fEsdTrackCuts->SetMaxChi2PerClusterTPC(3);
    fEsdTrackCuts->SetMaxDCAToVertexZ(2.7);
    fEsdTrackCuts->SetMaxDCAToVertexXY(1.9);
  }
  else if (filterbit == 512)
  {
    // syst study
    fEsdTrackCuts = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
    fEsdTrackCuts->SetMinNClustersTPC(60);
    fEsdTrackCuts->SetMaxChi2PerClusterTPC(5);
    fEsdTrackCuts->SetMaxDCAToVertexZ(3.7);
    fEsdTrackCuts->SetMaxDCAToVertexXY(2.9);
  }
  else if (filterbit == 1024)
  {
    fEsdTrackCuts = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
    fEsdTrackCuts->SetMinNClustersTPC(-1);
    fEsdTrackCuts->SetMinNCrossedRowsTPC(70);
    fEsdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
  }
  else if (filterbit == 2048) // mimic hybrid tracks
  {
    // correspond to esdTrackCutsHTG, but WITHOUT spd constraint. this is checked with the next object
    fEsdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE);
    fEsdTrackCuts->SetName("Global Hybrid tracks, loose DCA");
    fEsdTrackCuts->SetMaxDCAToVertexXY(2.4);
    fEsdTrackCuts->SetMaxDCAToVertexZ(3.2);
    fEsdTrackCuts->SetDCAToVertex2D(kTRUE);
    fEsdTrackCuts->SetMaxChi2TPCConstrainedGlobal(36);
    fEsdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kOff);
    fEsdTrackCuts->SetMaxFractionSharedTPCClusters(0.4);

    // Add SPD requirement 
    fEsdTrackCutsExtra1 = new AliESDtrackCuts("SPD", "Require 1 cluster in SPD");
    fEsdTrackCutsExtra1->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
    // A track passing fEsdTrackCuts and fEsdTrackCutsExtra1 corresponds to esdTrackCutsHTG

    fEsdTrackCutsExtra2 = new AliESDtrackCuts("No_SPD", "Reject tracks with cluster in SPD");
    fEsdTrackCutsExtra2->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kNone);
    // A track passing fEsdTrackCuts and fEsdTrackCutsExtra2 corresponds to esdTrackCutsHTGC and needs to be constrained
  }
  else if (filterbit == 2049) // mimic hybrid tracks for LHC11a p4 AOD113
  {
    // from http://svn.cern.ch/guest/AliRoot/tags/v5-02-Rev-20/ANALYSIS/macros/AddTaskESDFilter.C
    
    fEsdTrackCuts = new AliESDtrackCuts("AliESDtrackCuts"); 

    TFormula *f1NClustersTPCLinearPtDep = new TFormula("f1NClustersTPCLinearPtDep","70.+30./20.*x");
    fEsdTrackCuts->SetMinNClustersTPCPtDep(f1NClustersTPCLinearPtDep,20.);
    fEsdTrackCuts->SetMinNClustersTPC(70);
    fEsdTrackCuts->SetMaxChi2PerClusterTPC(4);
    fEsdTrackCuts->SetRequireTPCStandAlone(kTRUE); //cut on NClustersTPC and chi2TPC Iter1
    fEsdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    fEsdTrackCuts->SetRequireTPCRefit(kTRUE);
    fEsdTrackCuts->SetMaxFractionSharedTPCClusters(0.4);
    // ITS
    // moved down! fEsdTrackCuts->SetRequireITSRefit(kTRUE);
    fEsdTrackCuts->SetRequireITSRefit(kFALSE);
    //accept secondaries
    fEsdTrackCuts->SetMaxDCAToVertexXY(2.4);
    fEsdTrackCuts->SetMaxDCAToVertexZ(3.2);
    fEsdTrackCuts->SetDCAToVertex2D(kTRUE);
    //reject fakes
    fEsdTrackCuts->SetMaxChi2PerClusterITS(36);
    fEsdTrackCuts->SetMaxChi2TPCConstrainedGlobal(36);

    fEsdTrackCuts->SetRequireSigmaToVertex(kFALSE);

    // Add SPD requirement 
    fEsdTrackCutsExtra1 = new AliESDtrackCuts("SPD", "Require 1 cluster in SPD");
    fEsdTrackCutsExtra1->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
    fEsdTrackCutsExtra1->SetRequireITSRefit(kTRUE);
    // A track passing fEsdTrackCuts and fEsdTrackCutsExtra1 corresponds to esdTrackCutsHG0 (bit 16)

    fEsdTrackCutsExtra2 = new AliESDtrackCuts("No_ITS", "No ITS refit");
    // already implicit above: fEsdTrackCutsExtra2->SetRequireITSRefit(kFALSE);
    // A track passing fEsdTrackCuts and fEsdTrackCutsExtra2 corresponds to esdTrackCutsGCOnly and needs to be constrained (bit 256)
  }  
  else if (filterbit == 4096) // require TOF matching
  {
    fEsdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kTRUE);
    // FIXME: TOF REQUIREMENTS ARE IN GetParticleSpecies FOR THE MOMENT
  }
  else
  {
    fEsdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011();
    fEsdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kOff);

    // Add SPD requirement 
    fEsdTrackCutsExtra1 = new AliESDtrackCuts("SPD", "Require 1 cluster in SPD");
    fEsdTrackCutsExtra1->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);

    // Add SDD requirement 
    fEsdTrackCutsExtra2 = new AliESDtrackCuts("SDD", "Require 1 cluster in first layer SDD");
    fEsdTrackCutsExtra2->SetClusterRequirementITS(AliESDtrackCuts::kSDD,AliESDtrackCuts::kFirst);
  }
}

//____________________________________________________________________
TObjArray*  AliAnalyseLeadingTrackUE::FindLeadingObjects(TObject *obj)
{

  // Returns an array of charged particles (or jets) ordered according to their pT.

  Int_t nTracks = NParticles(obj);


  if( !nTracks ) return 0;
 
  // Define array of AliVParticle objects
  TObjArray* tracks = new TObjArray(nTracks);

  // Loop over tracks or jets
  for (Int_t ipart=0; ipart<nTracks; ++ipart) {
  	AliVParticle* part = ParticleWithCuts( obj, ipart );
        if (!part) continue;
	// Accept leading-tracks in a limited pseudo-rapidity range	
  	if( TMath::Abs(part->Eta()) > fTrackEtaCut || TMath::Abs(part->Eta()) < fTrackEtaCutMin ) continue;
  	tracks->AddLast( part );
  	}
  // Order tracks by pT	
  QSortTracks( *tracks, 0, tracks->GetEntriesFast() );

  nTracks = tracks->GetEntriesFast();
  if( !nTracks ) return 0;

  return tracks;
  }


void AliAnalyseLeadingTrackUE::RemoveInjectedSignals(TObjArray* tracks, TObject* mcObj, Int_t maxLabel)
{
  // remove injected signals (primaries above <maxLabel>)
  // <tracks> can be the following cases:
  // a. tracks: in this case the label is taken and then case b.
  // b. particles: the first stable mother is searched and checked if it is <= <maxLabel>
  // <mcObj> can be AOD (TClonesArray) or ESD (AliMCEvent)
  
  TClonesArray* arrayMC = 0;
  AliMCEvent* mcEvent = 0;
  if (mcObj->InheritsFrom("AliMCEvent"))
    mcEvent = static_cast<AliMCEvent*>(mcObj);
  else if (mcObj->InheritsFrom("TClonesArray"))
    arrayMC = static_cast<TClonesArray*>(mcObj);
  else
  {
    mcObj->Dump();
    AliFatal("Invalid object passed");
  }
  
  Int_t before = tracks->GetEntriesFast();

  for (Int_t i=0; i<before; ++i) 
  {
    AliVParticle* part = (AliVParticle*) tracks->At(i);
    
    if (part->InheritsFrom("AliESDtrack") || part->InheritsFrom("AliAODTrack"))
      part = ((mcEvent) ? mcEvent->GetTrack(TMath::Abs(part->GetLabel())) : (AliVParticle*)arrayMC->At(TMath::Abs(part->GetLabel())));
      
    AliVParticle* mother = part;
    if (mcEvent)
    {
      while (!mcEvent->IsPhysicalPrimary(mother->GetLabel()))
      {
	if (((AliMCParticle*)mother)->GetMother() < 0)
	{
	  mother = 0;
	  break;
	}

	mother = (AliMCParticle*) mcEvent->GetTrack(((AliMCParticle*)mother)->GetMother());
	if (!mother)
	  break;
      }
    }
    else
    {
      // find the primary mother
      while (!((AliAODMCParticle*)mother)->IsPhysicalPrimary())
      {
	if (((AliAODMCParticle*)mother)->GetMother() < 0)
	{
	  mother = 0;
	  break;
	}
	  
	mother = (AliVParticle*) arrayMC->At(((AliAODMCParticle*)mother)->GetMother());
	if (!mother)
	  break;
      }
    }
    
    if (!mother)
      AliWarning(Form("WARNING: No mother found for particle %d. Removing it from list.", part->GetLabel()));

//     Printf("%d %d %d", i, part->GetLabel(), mother->GetLabel());
    if (!mother || mother->GetLabel() >= maxLabel)
    {
//       Printf("Removing %d with label %d", i, part->GetLabel()); ((AliMCParticle*)part)->Particle()->Print(); ((AliMCParticle*)mother)->Particle()->Print();
      TObject* object = tracks->RemoveAt(i);
      if (tracks->IsOwner())
	delete object;
    }
  }
 
  tracks->Compress();
  
  AliInfo(Form("Reduced from %d to %d", before, tracks->GetEntriesFast())); 
}

//-------------------------------------------------------------------
void AliAnalyseLeadingTrackUE::RemoveWeakDecays(TObjArray* tracks, TObject* mcObj)
{
  // remove particles from weak decays
  // <tracks> can be the following cases:
  // a. tracks: in this case the label is taken and then case b.
  // b. particles: it is checked if IsSecondaryFromWeakDecay is true
  // <mcObj> can be AOD (TClonesArray) or ESD (AliMCEvent)
  
  TClonesArray* arrayMC = 0;
  AliMCEvent* mcEvent = 0;
  if (mcObj->InheritsFrom("AliMCEvent"))
    mcEvent = static_cast<AliMCEvent*>(mcObj);
  else if (mcObj->InheritsFrom("TClonesArray"))
    arrayMC = static_cast<TClonesArray*>(mcObj);
  else
  {
    mcObj->Dump();
    AliFatal("Invalid object passed");
  }
  
  Int_t before = tracks->GetEntriesFast();

  for (Int_t i=0; i<before; ++i) 
  {
    AliVParticle* part = (AliVParticle*) tracks->At(i);
    
    if (part->InheritsFrom("AliESDtrack") || part->InheritsFrom("AliAODTrack"))
      part = ((mcEvent) ? mcEvent->GetTrack(TMath::Abs(part->GetLabel())) : (AliVParticle*)arrayMC->At(TMath::Abs(part->GetLabel())));
    
    if (part->InheritsFrom("AliAODMCParticle"))
    {
      if (!((AliAODMCParticle*) part)->IsSecondaryFromWeakDecay())
	continue;
    }
    else if (part->InheritsFrom("AliMCParticle") && mcEvent)
    {
      if (!(mcEvent->Stack()->IsSecondaryFromWeakDecay(part->GetLabel())))
	continue;
    }
    else
    {
      part->Dump();
      AliFatal("Unknown particle");
    }
    
//     Printf("Removing %d with label %d", i, part->GetLabel()); part->Dump();
    TObject* object = tracks->RemoveAt(i);
    if (tracks->IsOwner())
      delete object;
  }
 
  tracks->Compress();
  
  if (before > tracks->GetEntriesFast())
    AliInfo(Form("Reduced from %d to %d", before, tracks->GetEntriesFast())); 
}

//-------------------------------------------------------------------
TObjArray* AliAnalyseLeadingTrackUE::GetAcceptedParticles(TObject* obj, TObject* arrayMC, Bool_t onlyprimaries, Int_t particleSpecies, Bool_t useEtaPtCuts, Bool_t speciesOnTracks, Double_t evtPlane, Bool_t onlyCharged, ULong64_t generatorIndexMask)
{
  // Returns an array of particles that pass the cuts, if arrayMC is given each reconstructed particle is replaced by its corresponding MC particles, depending on the parameter onlyprimaries only for primaries 
  // particleSpecies: -1 all particles are returned
  //                  0 (pions) 1 (kaons) 2 (protons) 3 (others) particles
  // speciesOnTracks if kFALSE, particleSpecies is only applied on the matched MC particle (not on the track itself)
  // Passing down the Double_t* evtPlane (range [-pi/2,pi/2]) will apply a phi cut with respect to the eventplane between fTrackPhiCutEvPlMin and fTrackPhiCutEvPlMax. For values outside [-pi/2,pi/2], this will be ignored.
  
  Int_t nTracks = NParticles(obj);
  TObjArray* tracks = new TObjArray;
  
  // for TPC only tracks
  Bool_t hasOwnership = kFALSE;
  if ((fFilterBit == 128 || fFilterBit == 256 || fFilterBit == 512 || fFilterBit == 1024 || fFilterBit == 2048) && obj->InheritsFrom("AliESDEvent"))
    hasOwnership = kTRUE;
  
  if (!arrayMC)
    tracks->SetOwner(hasOwnership);
 
  // Loop over tracks or jets
  for (Int_t ipart=0; ipart<nTracks; ++ipart) {
    AliVParticle* part = ParticleWithCuts( obj, ipart, onlyprimaries, (speciesOnTracks) ? particleSpecies : -1, onlyCharged);
    if (!part) continue;
    
    if (TMath::Abs(evtPlane)<=TMath::Pi()/2) { //evtPlane range: (-pi/2,pi/2)
      Double_t phiPart = part->Phi(); //range: [0,2*pi)
      if(phiPart>TMath::Pi()) phiPart-=2*TMath::Pi();

      Double_t dPhi = 0; //range: [0,pi/2], i.e. the difference over the shortest angle.
      Double_t diff = TMath::Abs(phiPart-evtPlane);
      if(diff<=TMath::Pi()/2) dPhi = diff;
      else if(diff<=TMath::Pi()) dPhi = TMath::Pi()-diff;
      else dPhi = diff-TMath::Pi();
      
      if(dPhi<fTrackPhiCutEvPlMin || dPhi>fTrackPhiCutEvPlMax) {
        if (hasOwnership)
          delete part;
        continue;
      }
    }
    
    if (useEtaPtCuts)
      if (TMath::Abs(part->Eta()) > fTrackEtaCut || TMath::Abs(part->Eta()) < fTrackEtaCutMin || part->Pt() < fTrackPtMin)
      {
	if (hasOwnership)
	  delete part;
        continue;
      }
      
//     Printf("%p %p %d Accepted %d %f %f %f", obj, arrayMC, particleSpecies, ipart, part->Eta(), part->Phi(), part->Pt());
    
    if (arrayMC) {
      Int_t label = part->GetLabel();
      if (hasOwnership)
	delete part;
      // re-define part as the matched MC particle
      part = ParticleWithCuts(arrayMC, TMath::Abs(label),onlyprimaries, particleSpecies, onlyCharged);
      if (!part)continue;
    }
    
    ULong64_t genIndex = (ULong64_t)part->GetGeneratorIndex();
    if (generatorIndexMask != 0ull && (genIndex >= 64ull || (1ull<<genIndex & generatorIndexMask) == 0ull)){
      if (hasOwnership)
        delete part;
      continue;
    }
    
    tracks->Add(part);
  }

  return tracks;
}

//-------------------------------------------------------------------
TObjArray* AliAnalyseLeadingTrackUE::GetFakeParticles(TObject* obj, TObject* arrayMC, Bool_t onlyprimaries, Int_t particleSpecies, Bool_t useEtaPtCuts)
{
  // particleSpecies: -1 all particles are returned
  //                  0 (pions) 1 (kaons) 2 (protons) 3 (others) particles

  Int_t nTracks = NParticles(obj);
  TObjArray* tracksReconstructed = new TObjArray;
  TObjArray* tracksOriginal = new TObjArray;
  TObjArray* tracksFake = new TObjArray;

  // for TPC only tracks
  Bool_t hasOwnership = kFALSE;
  if ((fFilterBit == 128 || fFilterBit == 256 || fFilterBit == 512 || fFilterBit == 1024 || fFilterBit == 2048) && obj->InheritsFrom("AliESDEvent"))
    hasOwnership = kTRUE;

  tracksReconstructed->SetOwner(hasOwnership);
  tracksFake->SetOwner(hasOwnership);

  // Loop over tracks or jets
  for (Int_t ipart=0; ipart<nTracks; ++ipart) {
    AliVParticle* partReconstructed = ParticleWithCuts( obj, ipart, onlyprimaries, particleSpecies );
    if (!partReconstructed) continue;

    if (useEtaPtCuts)
      if (TMath::Abs(partReconstructed->Eta()) > fTrackEtaCut || TMath::Abs(partReconstructed->Eta()) < fTrackEtaCutMin || partReconstructed->Pt() < fTrackPtMin)
      {
        if (hasOwnership)
          delete partReconstructed;
        continue;
      }

    Int_t label = partReconstructed->GetLabel();
    if (label == 0)
    {
      /*
      Printf(">>> TPC only track:");
      partReconstructed->Print();
      partReconstructed->Dump();
      Printf(">>> Global track:");
      ((AliESDEvent*) obj)->GetTrack(ipart)->Dump();
      Printf("Fake (TPC only): eta = %f, phi = %f, pT = %f, ncl = %d, dedx = %f", partReconstructed->Eta(), partReconstructed->Phi(), partReconstructed->Pt(), ((AliESDtrack*) partReconstructed)->GetTPCclusters(0), ((AliESDtrack*) partReconstructed)->GetTPCsignal());
      Printf("Fake (global  ): eta = %f, phi = %f, pT = %f, ncl = %d, dedx = %f", ((AliESDEvent*) obj)->GetTrack(ipart)->Eta(), ((AliESDEvent*) obj)->GetTrack(ipart)->Phi(), ((AliESDEvent*) obj)->GetTrack(ipart)->Pt(), ((AliESDEvent*) obj)->GetTrack(ipart)->GetTPCclusters(0), ((AliESDEvent*) obj)->GetTrack(ipart)->GetTPCsignal());
      */
      tracksFake->AddLast(partReconstructed);
      continue;
    }

    AliVParticle* partOriginal = ParticleWithCuts(arrayMC, TMath::Abs(label),onlyprimaries, particleSpecies);
    if (!partOriginal)
    {
      if (hasOwnership)
	delete partReconstructed;
      continue;
    }

    tracksReconstructed->AddLast(partReconstructed);
    tracksOriginal->AddLast(partOriginal);
  }
  TObjArray* pairs = new TObjArray;
  pairs->SetOwner(kTRUE);
  pairs->Add(tracksReconstructed);
  pairs->Add(tracksOriginal);
  pairs->Add(tracksFake);
  return pairs;
}

//-------------------------------------------------------------------
TObjArray* AliAnalyseLeadingTrackUE::GetMinMaxRegion(TList *transv1, TList *transv2)
{
  
 // Returns two lists of particles, one for MIN and one for MAX region
  Double_t sumpT1 = 0.;
  Double_t sumpT2 = 0.;

  Int_t particles1 = transv1->GetEntries();
  Int_t particles2 = transv2->GetEntries();
  
  // Loop on transverse region 1
  for (Int_t i=0; i<particles1; i++){
        AliVParticle *part = (AliVParticle*)transv1->At(i);
	sumpT1 +=  part->Pt();
  	}

  // Loop on transverse region 2
  for (Int_t i=0; i<particles2; i++){
        AliVParticle *part = (AliVParticle*)transv2->At(i);
	sumpT2 +=  part->Pt();
  	}

  TObjArray *regionParticles = new TObjArray;
  if ( sumpT2 >= sumpT1 ){
  	regionParticles->AddLast(transv1); // MIN
  	regionParticles->AddLast(transv2); // MAX 
  }else {
  	regionParticles->AddLast(transv2); // MIN
  	regionParticles->AddLast(transv1); // MAX
  	}

 return regionParticles;
}

//-------------------------------------------------------------------
Int_t  AliAnalyseLeadingTrackUE::NParticles(TObject* obj)
{
 
  //Returns the number of particles in AliAODMCParticle array  or AliAODTracks or AliESDTracks 

  Int_t nTracks;
  
  if (obj->InheritsFrom("TClonesArray")){ // MC particles
  	TClonesArray *arrayMC = static_cast<TClonesArray*>(obj);
        nTracks = arrayMC->GetEntriesFast();
  }else if (obj->InheritsFrom("TObjArray")){ // list of AliVParticle
  	TObjArray *array = static_cast<TObjArray*>(obj);
        nTracks = array->GetEntriesFast();
  }else if (obj->InheritsFrom("AliAODEvent")){  // RECO AOD tracks
  	AliAODEvent *aodEvent = static_cast<AliAODEvent*>(obj);
        nTracks = aodEvent->GetNumberOfTracks();
  }else if (obj->InheritsFrom("AliESDEvent")){  // RECO ESD tracks
  	AliESDEvent *esdEvent = static_cast<AliESDEvent*>(obj);
        nTracks = esdEvent->GetNumberOfTracks();
  }else if (obj->InheritsFrom("AliMCEvent")){  // RECO ESD tracks
  	AliMCEvent *mcEvent = static_cast<AliMCEvent*>(obj);
        nTracks = mcEvent->GetNumberOfTracks();
  }else {
  	if (fDebug > 1) AliFatal(" Analysis type not defined !!! ");
	return 0;
	}
  
  return nTracks;
}

//-------------------------------------------------------------------
AliVParticle*  AliAnalyseLeadingTrackUE::ParticleWithCuts(TObject* obj, Int_t ipart, Bool_t onlyprimaries, Int_t particleSpecies, Bool_t onlyCharged)
{
  // Returns track or MC particle at position "ipart" if passes selection criteria
  // particleSpecies: -1 all particles are returned
  //                  0 (pions) 1 (kaons) 2 (protons) 3 (others) particles
  AliVParticle *part=0;
  
  if (obj->InheritsFrom("TClonesArray")){ // AOD-MC PARTICLE
  	TClonesArray *arrayMC = static_cast<TClonesArray*>(obj);
        part = (AliVParticle*)arrayMC->At( ipart );
	if (!part)return 0;
	// eventually only primaries
	if (onlyprimaries && !( ((AliAODMCParticle*)part)->IsPhysicalPrimary()) )return 0;
        // eventually only hadrons
        if (fOnlyHadrons){
	        Int_t pdgCode = ((AliAODMCParticle*)part)->GetPdgCode();
  		Bool_t isHadron = TMath::Abs(pdgCode)==211 ||  // Pion
  				  TMath::Abs(pdgCode)==2212 || // Proton
  	                          TMath::Abs(pdgCode)==321;    // Kaon
                if (!isHadron) return 0;				  
		}
        if (particleSpecies != -1) {
                // find the primary mother
                AliVParticle* mother = part;
                if(fCheckMotherPDG) {
		  while (!((AliAODMCParticle*)mother)->IsPhysicalPrimary())
		  {
                  if (((AliAODMCParticle*)mother)->GetMother() < 0)
                  {
                    mother = 0;
                    break;
                  }
                    
                  mother = (AliVParticle*) arrayMC->At(((AliAODMCParticle*)mother)->GetMother());
                  if (!mother)
                    break;
                 }
                }
                if (mother)
                {
                  Int_t pdgCode = ((AliAODMCParticle*)mother)->GetPdgCode();
                  if (particleSpecies == 0 && TMath::Abs(pdgCode)!=211)
                          return 0;
                  if (particleSpecies == 1 && TMath::Abs(pdgCode)!=321)
                          return 0;
                  if (particleSpecies == 2 && TMath::Abs(pdgCode)!=2212)
                          return 0;
                  if (particleSpecies == 3 && (TMath::Abs(pdgCode)==211 || TMath::Abs(pdgCode)==321 || TMath::Abs(pdgCode)==2212))
                          return 0;
                }
                else
                {
                  // if mother not found, accept particle only in case of particleSpecies == 3. To include it in all or no sample is no solution
                  AliWarning(Form("WARNING: No mother found for particle %d:", part->GetLabel()));
                  //part->Print();
  
                  /*
                  // this code prints the details of the mother that is missing in the AOD
                  AliMCEventHandler* fMcHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
  
                  AliMCEvent* fMcEvent = fMcHandler->MCEvent();
                  
                  fMcEvent->Stack()->Particle(fMcEvent->Stack()->Particle(part->GetLabel())->GetMother(0))->Print();
                  fMcEvent->Stack()->Particle(fMcEvent->Stack()->Particle(fMcEvent->Stack()->Particle(part->GetLabel())->GetMother(0))->GetMother(0))->Print();
                  Printf("eta = %f", fMcEvent->Stack()->Particle(fMcEvent->Stack()->Particle(fMcEvent->Stack()->Particle(part->GetLabel())->GetMother(0))->GetMother(0))->Eta());
                  */
                  
                  if (particleSpecies != 3)
                    return 0;
                }
        }
  
  }else if (obj->InheritsFrom("TObjArray")){ // list of AliVParticle
  	TObjArray *array = static_cast<TObjArray*>(obj);
        part = (AliVParticle*)array->At( ipart );
	if (!part)return 0;
  }else if (obj->InheritsFrom("AliMCEvent")){ // MC PARTICLE
        AliMCEvent* mcEvent =  static_cast<AliMCEvent*>(obj);
	part = mcEvent->GetTrack( ipart );
	if (!part) return 0;
	// eventually only primaries
	if (onlyprimaries && !( mcEvent->IsPhysicalPrimary(ipart)) )return 0;
        // eventually only hadrons
        //TO-DO
	/*if (fOnlyHadrons){
	        Int_t pdgCode =  part->GetPdgCode();
  		Bool_t isHadron = TMath::Abs(pdgCode)==211 ||  // Pion
  				  TMath::Abs(pdgCode)==2212 || // Proton
  	                          TMath::Abs(pdgCode)==321;    // Kaon
                if (!isHadron) return 0;				  
		}
	*/
       if (particleSpecies != -1) {
                // find the primary mother
                AliMCParticle* mother = (AliMCParticle*) part;
// 		Printf("");
                if(fCheckMotherPDG) {
		  while (!mcEvent->IsPhysicalPrimary(mother->GetLabel()))
                  {
// 		  Printf("pdg = %d; mother = %d", mother->PdgCode(), mother->GetMother());
                  if (mother->GetMother() < 0)
                  {
                    mother = 0;
                    break;
                  }
                    
                  mother = (AliMCParticle*) mcEvent->GetTrack(mother->GetMother());
                  if (!mother)
                    break;
                  }
                }
                if (mother)
                {
                  Int_t pdgCode = mother->PdgCode();
                  if (particleSpecies == 0 && TMath::Abs(pdgCode)!=211)
                          return 0;
                  if (particleSpecies == 1 && TMath::Abs(pdgCode)!=321)
                          return 0;
                  if (particleSpecies == 2 && TMath::Abs(pdgCode)!=2212)
                          return 0;
                  if (particleSpecies == 3 && (TMath::Abs(pdgCode)==211 || TMath::Abs(pdgCode)==321 || TMath::Abs(pdgCode)==2212))
                          return 0;
                }
                else
                {
                  // if mother not found, accept particle only in case of particleSpecies == 3. To include it in all or no sample is no solution
                  Printf("WARNING: No mother found for particle %d:", part->GetLabel());
                  //part->Dump();
                  //part->Print();
  
                  /*
                  // this code prints the details of the mother that is missing in the AOD
                  AliMCEventHandler* fMcHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
  
                  AliMCEvent* fMcEvent = fMcHandler->MCEvent();
                  
                  fMcEvent->Stack()->Particle(fMcEvent->Stack()->Particle(part->GetLabel())->GetMother(0))->Print();
                  fMcEvent->Stack()->Particle(fMcEvent->Stack()->Particle(fMcEvent->Stack()->Particle(part->GetLabel())->GetMother(0))->GetMother(0))->Print();
                  Printf("eta = %f", fMcEvent->Stack()->Particle(fMcEvent->Stack()->Particle(fMcEvent->Stack()->Particle(part->GetLabel())->GetMother(0))->GetMother(0))->Eta());
                  */
                  
                  if (particleSpecies != 3)
                    return 0;
                }
        }
  }else if (obj->InheritsFrom("AliAODEvent")){ // RECO AOD TRACKS
  	AliAODEvent *aodEvent = static_cast<AliAODEvent*>(obj);
        part = aodEvent->GetTrack(ipart);
        
        AliAODTrack* aodTrack = dynamic_cast<AliAODTrack*> (part);
	
	// track selection cuts
	if (fFilterBit != 0 && !(aodTrack->TestFilterBit(fFilterBit))) return 0; 
	if (fTrackStatus != 0 && !CheckTrack(part)) return 0;
	
	// DCA XY
	if (fDCAXYCut || fDCAZCut)
	{
	  const AliVVertex* vertex = aodEvent->GetPrimaryVertex();
	  if (!vertex)
	    return 0;
	  
	  Double_t pos[2];
	  Double_t covar[3];
	  AliAODTrack* clone = dynamic_cast<AliAODTrack*> (part->Clone());
	  Bool_t success = clone->PropagateToDCA(vertex, ((AliVAODHeader*)aodEvent->GetHeader())->GetMagneticField(), 3, pos, covar);
	  delete clone;
	  if (!success)
	    return 0;

// 	  Printf("%f", ((AliAODTrack*)part)->DCA());
// 	  Printf("%f", pos[0]);
	  if ((fDCAXYCut && TMath::Abs(pos[0]) > fDCAXYCut->Eval(part->Pt()))
		|| (fDCAZCut && TMath::Abs(pos[1]) > fDCAZCut->Eval(part->Pt())))
	    return 0;
	}
	
	if (fSharedClusterCut >= 0)
	{
	  Double_t frac = Double_t(aodTrack->GetTPCnclsS()) / Double_t(aodTrack->GetTPCncls());
	  if (frac > fSharedClusterCut)
	    return 0;
	}
	
	if (fCrossedRowsCut >= 0)
	{
	  if (aodTrack->GetTPCNCrossedRows() < fCrossedRowsCut)
	    return 0;
	}
	
	if (fFoundFractionCut >= 0)
	{
	  UInt_t findableClusters = aodTrack->GetTPCNclsF();
	  if (findableClusters == 0)
	    return 0;
	  if (((Double_t) aodTrack->GetTPCNCrossedRows() / findableClusters) < fFoundFractionCut)
	    return 0;
	}

	// eventually only hadrons
	if (fOnlyHadrons){
		Bool_t isHadron = aodTrack->GetMostProbablePID()==AliAODTrack::kPion ||
	                          aodTrack->GetMostProbablePID()==AliAODTrack::kKaon ||
		                  aodTrack->GetMostProbablePID()==AliAODTrack::kProton;
        	if (!isHadron) return 0;				  
		}
		
	if (particleSpecies != -1 && fHelperPID->GetParticleSpecies((AliVTrack*) part,kTRUE) != particleSpecies) return 0;
  
  }else if (obj->InheritsFrom("AliESDEvent")){ // RECO ESD TRACKS
  	AliESDEvent *esdEvent = static_cast<AliESDEvent*>(obj);
        part = esdEvent->GetTrack(ipart);
	if (!part)return 0;

	// track selection cuts
	if (!( ApplyCuts(part)) )
	  return 0; 
	
	if (fTrackStatus != 0 && !CheckTrack(part)) return 0;

	if (fFilterBit == 128 || fFilterBit == 256 || fFilterBit == 512 || fFilterBit == 1024)
	{
	  // create TPC only tracks constrained to the SPD vertex

	  const AliESDVertex *vtxSPD = esdEvent->GetPrimaryVertexSPD();

	  AliESDtrack* track = AliESDtrackCuts::GetTPCOnlyTrack(esdEvent, ipart);
	  if(!track) return 0;
    
// 	  Printf(">%f %f %f", track->Eta(), track->Phi(), track->Pt());
	  
	  if(track->Pt()>0.){
	    // only constrain tracks above threshold
	    AliExternalTrackParam exParam;
	    // take the B-feild from the ESD, no 3D fieldMap available at this point
	    Bool_t relate = kFALSE;
	    relate = track->RelateToVertexTPC(vtxSPD,esdEvent->GetMagneticField(),kVeryBig,&exParam);
	    if(!relate)
	    {
//                 Printf("relating failed");
	      delete track;
	      return 0;
	    }
	    track->Set(exParam.GetX(),exParam.GetAlpha(),exParam.GetParameter(),exParam.GetCovariance());
	  }
	  
// 	  Printf(">%f %f %f\n", track->Eta(), track->Phi(), track->Pt());

	  part = track;
	}
	else if (fFilterBit == 2048)
	{
	  // hybrid tracks
	  
	  // clone
	  AliESDtrack* esdTrack = new AliESDtrack(*((AliESDtrack*) part));
// 	  Printf("%d %d %d %d %d", fEsdTrackCuts->IsSelected(esdTrack), fEsdTrackCutsExtra1->IsSelected(esdTrack), fEsdTrackCutsExtra2->IsSelected(esdTrack), esdTrack->HasPointOnITSLayer(0), esdTrack->HasPointOnITSLayer(1));
	  
	  if (fEsdTrackCutsExtra2->IsSelected(esdTrack))
	  {
// 	    Float_t ptBefore = esdTrack->Pt();
	    // set constrained pT as default pT
	    if (!esdTrack->GetConstrainedParam())
	      return 0;
	    esdTrack->CopyFromVTrack(esdTrack->GetConstrainedParam());
// 	    Printf("%f %f", ptBefore, esdTrack->Pt());
	  }
	  part = esdTrack;
	}
	
	// eventually only hadrons
	//TO-DO
	/*if (fOnlyHadrons){
		Bool_t isHadron = ((AliAODTrack*)part)->GetMostProbablePID()==AliAODTrack::kPion ||
	                          ((AliAODTrack*)part)->GetMostProbablePID()==AliAODTrack::kKaon ||
		                  ((AliAODTrack*)part)->GetMostProbablePID()==AliAODTrack::kProton;
        	if (!isHadron) return 0;				  
		}
	*/	
	if (particleSpecies != -1 && fHelperPID->GetParticleSpecies((AliVTrack*) part,kTRUE) != particleSpecies) return 0; // If it is -1 you take all the particles

  }else {
  	if (fDebug > 1) AliFatal(" Analysis type not defined !!! ");
	return 0;
	}
  
  // only charged
  if (onlyCharged && !part->Charge())return 0;
  
  part->SetUniqueID((fEventCounter * 50000 + ipart) * 10);
  return part;
}


//-------------------------------------------------------------------
void  AliAnalyseLeadingTrackUE::QSortTracks(TObjArray &a, Int_t first, Int_t last)
{
  // Sort array of TObjArray of tracks by Pt using a quicksort algorithm.
  
  static TObject *tmp;
  static int i;           // "static" to save stack space
  int j;
  
  while (last - first > 1) {
    i = first;
    j = last;
    for (;;) {
      while (++i < last && ((AliVParticle*)a[i])->Pt() > ((AliVParticle*)a[first])->Pt() )
        ;
      while (--j > first && ((AliVParticle*)a[j])->Pt() < ((AliVParticle*)a[first])->Pt() )
        ;
      if (i >= j)
        break;
      
      tmp  = a[i];
      a[i] = a[j];
      a[j] = tmp;
    }
    if (j == first) {
      ++first;
      continue;
    }
    tmp = a[first];
    a[first] = a[j];
    a[j] = tmp;
    if (j - first < last - (j + 1)) {
      QSortTracks(a, first, j);
      first = j + 1;   // QSortTracks(j + 1, last);
    } else {
      QSortTracks(a, j + 1, last);
      last = j;        // QSortTracks(first, j);
    }
  }
}

//____________________________________________________________________
TObjArray*  AliAnalyseLeadingTrackUE::SortRegions(const AliVParticle* leading, TObject* obj, TObject* arrayMC, Bool_t onlyprimaries)
{

  // Assign particles to towards, away or transverse regions.
  // Returns a lists of particles for each region.

  static const Double_t k60rad  = 60.*TMath::Pi()/180.;
  static const Double_t k120rad = 120.*TMath::Pi()/180.;

  // Define output lists of particles
  TList *toward = new TList();
  TList *away = new TList();
  // Two transverse regions, for the moment those are not yet MIN and MAX!!! 
  // MIN and MAX can be sorted in GetMinMaxRegion function
  TList *transverse1 = new TList();
  TList *transverse2 = new TList();
  
  TObjArray *regionParticles = new TObjArray;
  regionParticles->SetOwner(kTRUE);
  
  regionParticles->AddLast(toward);
  regionParticles->AddLast(away);
  regionParticles->AddLast(transverse1);
  regionParticles->AddLast(transverse2);
  
  if (!leading)
    return regionParticles;
 
  // Switch to vector for leading particle
  TVector3 leadVect(leading->Px(),leading->Py(),leading->Pz());
  
  Int_t nTracks = NParticles(obj);
  if( !nTracks ) return 0;
  // Loop over tracks 
  for (Int_t ipart=0; ipart<nTracks; ++ipart) {
  	AliVParticle* part = ParticleWithCuts(obj, ipart);
	if (!part)continue;
	//Switch to vectors for particles 
  	TVector3 partVect(part->Px(), part->Py(), part->Pz());
 
  	Int_t region = 0;
  	if( TMath::Abs(partVect.Eta()) > fTrackEtaCut || TMath::Abs(partVect.Eta()) < fTrackEtaCutMin) continue;
  	// transverse regions
  	if (leadVect.DeltaPhi(partVect) < -k60rad && leadVect.DeltaPhi(partVect) > -k120rad )region = -1; //left
  	if (leadVect.DeltaPhi(partVect) > k60rad && leadVect.DeltaPhi(partVect) < k120rad ) region = 1;   //right

  	if (TMath::Abs(leadVect.DeltaPhi(partVect)) < k60rad ) region = 2;    //forward
  	if (TMath::Abs(leadVect.DeltaPhi(partVect)) > k120rad ) region = -2;  //backward
  	
  	// skip leading particle 
        if (leading == part)
  	  continue;
  	
	if (!region)continue;
        if (arrayMC && arrayMC->InheritsFrom("TClonesArray") && obj->InheritsFrom("AliAODEvent")){
		Int_t label = ((AliAODTrack*)part)->GetLabel();
		// re-define part as the matched MC particle
		part = (AliAODMCParticle*)ParticleWithCuts(arrayMC, TMath::Abs(label),onlyprimaries);
		if (!part)continue;
                // skip leading particle 
                if (leading == part)
                  continue;
		}
        if (arrayMC && arrayMC->InheritsFrom("AliMCEvent") && obj->InheritsFrom("AliESDEvent")){
		Int_t label = ((AliESDtrack*)part)->GetLabel();
		// look for the matched MC particle (but do not re-define part)
		if (!ParticleWithCuts(arrayMC, TMath::Abs(label),onlyprimaries)) continue;
		}

	if ( region == 1 ) transverse1->Add(part);
	if ( region == -1 ) transverse2->Add(part);
	if ( region == 2 ) toward->Add(part);
	if ( region == -2 ) away->Add(part);

	}//end loop on tracks
  
  return regionParticles;
  
}


//____________________________________________________________________
Bool_t  AliAnalyseLeadingTrackUE::TriggerSelection(const TObject* obj)
{
  if (!obj) // MC
    return kFALSE;

  // Use AliPhysicsSelection to select good events, works for ESD and AOD
  if (!(((AliInputEventHandler*)obj)->IsEventSelected()&(fEventSelection)))
    return kFALSE;

  return kTRUE;
}

//____________________________________________________________________
Bool_t  AliAnalyseLeadingTrackUE::VertexSelection(const TObject* obj, Int_t ntracks, Double_t zed)
{

  //Require 1 vertex (no TPC stand-alone) with a minimum number of tracks and z-coordinate in a limited range
  
  if (obj->InheritsFrom("AliAODEvent")){ 
  	Int_t nVertex = ((AliAODEvent*)obj)->GetNumberOfVertices();
  	if( nVertex > 0 ) { 
  		AliAODVertex* vertex = (AliAODVertex*)((AliAODEvent*)obj)->GetPrimaryVertex();
		Int_t nTracksPrim = vertex->GetNContributors();
  		Double_t zVertex = vertex->GetZ();
  		if (fDebug > 1)AliInfo(Form(" Vertex in = %f with %d particles by  %s data ...",zVertex,nTracksPrim,vertex->GetName()));
  		// Reject TPC only vertex
		TString name(vertex->GetName());
		if (name.CompareTo("PrimaryVertex") && name.CompareTo("SPDVertex"))return kFALSE;

		// Select a quality vertex by number of tracks?
  		if( nTracksPrim < ntracks || TMath::Abs(zVertex) >= zed ) {
  			if (fDebug > 1) AliInfo(" Primary-vertex Selection: event REJECTED ...");
  			return kFALSE;
  			}
  		// TODO remove vertexer Z events with dispersion > 0.02: Doesn't work for AOD at present
                //if (strcmp(vertex->GetTitle(), "AliVertexerZ") == 0 && vertex->GetDispersion() > 0.02)
                //  return kFALSE;
  		if (fDebug > 1) AliInfo(" Primary-vertex Selection: event ACCEPTED...");
		} else {
  			if (fDebug > 1) AliInfo(" Primary-vertex Selection: event REJECTED ...");
  			return kFALSE;
  			}
	}

  if (obj->InheritsFrom("AliMCEvent"))
  { 
    if (TMath::Abs(((AliMCEvent*) obj)->GetPrimaryVertex()->GetZ()) >= zed)
    {
      if (fDebug > 1) AliInfo(" Primary-vertex Selection: event (based on MC) REJECTED ...");
      return kFALSE;
    }
  }

  if (obj->InheritsFrom("AliAODMCHeader"))
  { 
    if (TMath::Abs(((AliAODMCHeader*) obj)->GetVtxZ()) >= zed)
    {
      if (fDebug > 1) AliInfo(" Primary-vertex Selection: event (based on MC) REJECTED ...");
      return kFALSE;
    }
  }

  // ESD case for DCA studies
  if (obj->InheritsFrom("AliESDEvent")){
       AliESDVertex* vertex = (AliESDVertex*)((AliESDEvent*)obj)->GetPrimaryVertex();
       if ( vertex){
               Int_t nTracksPrim = vertex->GetNContributors();
               Double_t zVertex = vertex->GetZ();
               if (fDebug > 1)AliInfo(Form(" Vertex in = %f with %d particles by  %s data ...",zVertex,nTracksPrim,vertex->GetName()));
               // Reject SPD or TPC only vertex
               TString name(vertex->GetName());
               if (name.CompareTo("PrimaryVertex") && name.CompareTo("SPDVertex"))return kFALSE;

               // Select a quality vertex by number of tracks?
               if( nTracksPrim < ntracks || TMath::Abs(zVertex) >= zed ) {
                       if (fDebug > 1) AliInfo(" Primary-vertex Selection: event REJECTED ...");
                       return kFALSE;
                       }
               // TODO remove vertexer Z events with dispersion > 0.02: Doesn't work for AOD at present
                //if (strcmp(vertex->GetTitle(), "AliVertexerZ") == 0 && vertex->GetDispersion() > 0.02)
                //  return kFALSE;
               if (fDebug > 1) AliInfo(" Primary-vertex Selection: event ACCEPTED...");
               } else {
                       if (fDebug > 1) AliInfo(" Primary-vertex Selection: event REJECTED ...");
                       return kFALSE;
                       }
       }
	
  return kTRUE;
}

//____________________________________________________________________

Bool_t AliAnalyseLeadingTrackUE::CheckTrack(AliVParticle * part)
{
  // check if the track status flags are set
  
  UInt_t status=((AliVTrack*)part)->GetStatus();
  if ((status & fTrackStatus) == fTrackStatus)
    return kTRUE;
  return kFALSE;
}
