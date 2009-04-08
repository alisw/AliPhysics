/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

//
//  Analysis task for Kinematic filtering
//  Fill AODtrackMC tracks from Kinematic stack
//
 
#include <TChain.h>
#include <TFile.h>
#include "TParticle.h"
#include "TString.h"

#include "AliAnalysisTaskMCParticleFilter.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisFilter.h"
#include "AliHeader.h"
#include "AliStack.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliAODEvent.h"
#include "AliAODHeader.h"
#include "AliAODMCHeader.h"
#include "AliAODHandler.h"
#include "AliAODVertex.h"
#include "AliAODMCParticle.h"


#include "AliLog.h"

ClassImp(AliAnalysisTaskMCParticleFilter)

////////////////////////////////////////////////////////////////////////

//____________________________________________________________________
AliAnalysisTaskMCParticleFilter::AliAnalysisTaskMCParticleFilter():
AliAnalysisTaskSE(),
  fTrackFilterMother(0x0)
{
  // Default constructor
}

//____________________________________________________________________
AliAnalysisTaskMCParticleFilter::AliAnalysisTaskMCParticleFilter(const char* name):
    AliAnalysisTaskSE(name),
    fTrackFilterMother(0x0)
{
  // Default constructor
}


//____________________________________________________________________
AliAnalysisTaskMCParticleFilter::AliAnalysisTaskMCParticleFilter(const AliAnalysisTaskMCParticleFilter& obj):
    AliAnalysisTaskSE(obj),
    fTrackFilterMother(obj.fTrackFilterMother)
{
  // Copy constructor
}

//____________________________________________________________________
AliAnalysisTaskMCParticleFilter::~AliAnalysisTaskMCParticleFilter()
{
  //  if( fTrackFilterMother ) delete fTrackFilterMother;
}


//____________________________________________________________________
AliAnalysisTaskMCParticleFilter& AliAnalysisTaskMCParticleFilter::operator=(const AliAnalysisTaskMCParticleFilter& other)
{
// Assignment
  if(this!=&other) {
    AliAnalysisTaskSE::operator=(other);
    fTrackFilterMother = other.fTrackFilterMother;
  }
  return *this;
}

//____________________________________________________________________
void AliAnalysisTaskMCParticleFilter::UserCreateOutputObjects()
{
  // Create the output container
    if (OutputTree()&&fTrackFilterMother) 
	OutputTree()->GetUserInfo()->Add(fTrackFilterMother);

    // this part is mainly needed to set the MCEventHandler
    // to the AODHandler, this will not be needed when
    // AODHandler goes to ANALYSISalice
    // setting in the steering macro will not work on proof :(
    // so we have to do it in a task

    // the branch booking can also go into the AODHandler then


    // mcparticles
    TClonesArray *tca = new TClonesArray("AliAODMCParticle", 0);
    tca->SetName(AliAODMCParticle::StdBranchName());
    AddAODBranch("TClonesArray",&tca);

    // MC header...
    AliAODMCHeader *mcHeader = new AliAODMCHeader();
    Printf("AODMCHeader %p",mcHeader);
    Printf("AODMCHeader ** %p",&mcHeader);
    mcHeader->SetName(AliAODMCHeader::StdBranchName());
    AddAODBranch("AliAODMCHeader",&mcHeader);    

    

    AliMCEventHandler *mcH = (AliMCEventHandler*) ((AliAnalysisManager::GetAnalysisManager())->GetMCtruthEventHandler()); 
    AliAODHandler *aodH = dynamic_cast<AliAODHandler*> ((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());
    if(!aodH){
      Printf("%s:&d Could not get AODHandler",(char*)__FILE__,__LINE__);
      return;
    }
    aodH->SetMCEventHandler(mcH);

    
}

//____________________________________________________________________
void AliAnalysisTaskMCParticleFilter::UserExec(Option_t */*option*/)
{
// Execute analysis for current event
//

// Fill AOD tracks from Kinematic stack
    
  // get AliAOD Event 
  AliAODEvent* aod = AODEvent();
  if (!aod) {
      AliWarning("No Output Handler connected, doing nothing !") ;
      return;
  }

  AliMCEventHandler *mcH = (AliMCEventHandler*) ((AliAnalysisManager::GetAnalysisManager())->GetMCtruthEventHandler()); 
  if(!mcH){ 
    AliWarning("No MC handler Found");
    return;
  }


  // fetch the output 
  // Fill control histos

  // tmp array for holding the mctracks
  // need to set mother and duagther __before__
  // filling in the tree

  AliMCEvent *mcE = MCEvent();

  if(!mcE){
    AliWarning("No MC event Found");
    return;
  }

  AliStack* stack = mcE->Stack();
  Int_t np    = mcE->GetNumberOfTracks();
  Int_t nprim = stack->GetNprimary();
  // TODO ADD MC VERTEX
  
  // We take all real primaries


  Int_t j=0;
  for (Int_t ip = 0; ip < np; ip++){
    AliMCParticle* mcpart =  mcE->GetTrack(ip);
    TParticle* part = mcpart->Particle();
    Float_t xv = part->Vx();
    Float_t yv = part->Vy();
    Float_t zv = part->Vz();
    Float_t rv = TMath::Sqrt(xv * xv + yv * yv);
      
    Bool_t write = kFALSE;
    
    if (ip < nprim) {
      // Select the primary event
      write = kTRUE;
    } else if (part->GetUniqueID() == 4) {
      // Particles from decay
      // Check that the decay chain ends at a primary particle
      TParticle* mother = part;
      Int_t imo = part->GetFirstMother();
      while((imo >= nprim) && (mother->GetUniqueID() == 4)) {
	mother =  mcE->Stack()->Particle(imo);
	imo =  mother->GetFirstMother();
      }
      // Select according to pseudorapidity and production point of primary ancestor
      if (imo < nprim && Select(mcE->Stack()->Particle(imo), rv, zv))write = kTRUE;         
    } else if (part->GetUniqueID() == 5) {
      // Now look for pair production
      Int_t imo = part->GetFirstMother();
      if (imo < nprim) {
	// Select, if the gamma is a primary
	write = kTRUE;
      } else {
	// Check if the gamma comes from the decay chain of a primary particle
	TParticle* mother = mcE->Stack()->Particle(imo);
	imo = mother->GetFirstMother();
	while((imo >= nprim) && (mother->GetUniqueID() == 4)) {
	  mother =  mcE->Stack()->Particle(imo);
	  imo =  mother->GetFirstMother();
	}
	// Select according to pseudorapidity and production point 
	if (imo < nprim && Select(mother, rv, zv)) 
	  write = kTRUE;
      }
    }
    /*
    else if (part->GetUniqueID() == 13){
      // Evaporation
      // Check that we end at a primary particle
      TParticle* mother = part;
      Int_t imo = part->GetFirstMother();
      while((imo >= nprim) && ((mother->GetUniqueID() == 4 ) || ( mother->GetUniqueID() == 13))) {
	mother =  mcE->Stack()->Particle(imo);
	imo =  mother->GetFirstMother();
      }
      // Select according to pseudorapidity and production point 
	if (imo < nprim && Select(mother, rv, zv)) 
	  write = kTRUE;
    }
    */    
    if (write) {
      if(mcH)mcH->SelectParticle(ip);
      j++;
    }
  }


    // check for charm daughters
  static int  iTaken = 0;
  static int  iAll = 0;
  static int  iCharm = 0;
  for (Int_t ip = 0; ip < np; ip++){
    AliMCParticle* mcpart =  mcE->GetTrack(ip);
    TParticle* part = mcpart->Particle();

    //    if((TMath::Abs(part->GetPdgCode())/400)==1){
    if((TMath::Abs(part->GetPdgCode()))==411){
      // cases 
      iCharm++;
      Printf("Decay Mother %s",part->GetPDG()->GetName());
      Int_t d0 =  part->GetFirstDaughter();
      Int_t d1 =  part->GetLastDaughter();
      if(d0>0&&d1>0){
	for(int id = d0;id <= d1;id++){
	  TParticle* daughter =  mcE->Stack()->Particle(id);
	  Printf("Decay Daughter %s",daughter->GetPDG()->GetName());
	  iAll++;
	  if(mcH->IsParticleSelected(id))iTaken++;
	}
      }
    }
  }
  Printf("Taken daughters %d/%d of %d charm",iTaken,iAll,iCharm);


  return;
}

Bool_t AliAnalysisTaskMCParticleFilter::Select(TParticle* part, Float_t rv, Float_t zv)
{
  // Selection accoring to eta of the mother and production point
  // This has to be refined !!!!!!

  // Esp if we don't have collisison in the central barrel but e.g. beam gas
  //  return kTRUE;

  Float_t eta = part->Eta();

  // central barrel consider everything in the ITS...
  // large eta window for smeared vertex and SPD acceptance (2 at z = 0)
  // larger for V0s in the TPC
  //  if (TMath::Abs(eta) < 2.5 && rv < 250. && TMath::Abs(zv)<255)return kTRUE;

  if (TMath::Abs(eta) < 2.5 && rv < 170)return kTRUE;   

  // Andreas' Cuts
  //  if (TMath::Abs(eta) < 1. && rv < 170)return kTRUE;   



  // Muon arm
  if(eta > -4.2 && eta < -2.3 && zv > -500.)return kTRUE; // Muon arms

  // PMD acceptance 2.3 <= eta < = 3.5
  //  if(eta>2.0&&eta<3.7)return kTRUE; 

  return kFALSE;
 
}

