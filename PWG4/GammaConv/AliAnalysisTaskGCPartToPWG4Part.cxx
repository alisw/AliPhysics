#include <iostream>
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TString.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskGCPartToPWG4Part.h"

#include "AliESDEvent.h"
#include "AliESDCaloCluster.h"
#include "AliESDInputHandler.h"

#include "AliAODPWG4ParticleCorrelation.h"
#include "AliAODEvent.h"
#include "AliAODHandler.h"
#include "AliAODCaloCluster.h"
#include "AliGammaConversionAODObject.h"
#include "AliAODConversionParticle.h"
#include "AliAODJet.h"

#include "AliAODInputHandler.h"

#include "AliAODMCParticle.h"

#include "AliMCAnalysisUtils.h"

// Gamma - jet correlation analysis task
// Authors: Svein Lindal


using namespace std;

ClassImp(AliAnalysisTaskGCPartToPWG4Part)

//________________________________________________________________________
AliAnalysisTaskGCPartToPWG4Part::AliAnalysisTaskGCPartToPWG4Part() 
: AliAnalysisTaskSE(), 
  fDeltaAODFileName(""),
  fAODBranchName("GammaConv_gamma"),
  fAODPWG4Particles(NULL),
  fAnaUtils(NULL)
{
  // Dummy Constructor
}

//________________________________________________________________________________
AliAnalysisTaskGCPartToPWG4Part::~AliAnalysisTaskGCPartToPWG4Part() {

  if(fAODPWG4Particles)
    fAODPWG4Particles = NULL;
  delete fAODPWG4Particles;

  if(fAnaUtils)
    delete fAnaUtils;
  fAnaUtils = NULL;

}



//________________________________________________________________________
AliAnalysisTaskGCPartToPWG4Part::AliAnalysisTaskGCPartToPWG4Part(const char *name) : 
  AliAnalysisTaskSE(name), 
  fDeltaAODFileName(""),
  fAODBranchName("GammaConv_gamma"),
  fAODPWG4Particles(NULL),
  fAnaUtils(NULL)
{
  // Constructor
  DefineInput(0, TChain::Class());
  // Output slot #0 id reserved by the base class for AOD

  // Output slot #1 writes into a TH1 container
  DefineOutput(1, TList::Class());

  fAnaUtils = new AliMCAnalysisUtils();
}



//________________________________________________________________________
void AliAnalysisTaskGCPartToPWG4Part::UserCreateOutputObjects() {
  fAODPWG4Particles = new TClonesArray("AliAODPWG4ParticleCorrelation", 0);
  fAODPWG4Particles->SetName("ConversionGamma");
  AddAODBranch("TClonesArray", &fAODPWG4Particles);

}

//________________________________________________________________________
void AliAnalysisTaskGCPartToPWG4Part::UserExec(Option_t *) 
{
  
  //AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler()->SetFillAOD(kTRUE);

  //Clear stuff for new event
  CleanUp();

  ///Get AOD event
  AliAODEvent * aodEvent = GetAODEvent();
  if(!aodEvent) {
    AliError("No AOD event!!");
    return;
  }
  

  ProcessConvGamma(aodEvent);
    

  //PostData(1, fOutputList);
        
}



//___________________________________________________________________________________________
void AliAnalysisTaskGCPartToPWG4Part::ProcessConvGamma( const AliAODEvent * const aodEvent ) {
  
  TClonesArray * tracks = aodEvent->GetTracks();
  if(!tracks) {
    cout << "No tracks!!!"<<endl;
    return;
  }


  TClonesArray * convGamma = GetConversionGammas(aodEvent);
  if(!convGamma) {
    AliError(Form("No branch by name %s found in file %s", fAODBranchName.Data(), fDeltaAODFileName.Data()));
    return;
  }
  
  TClonesArray * arrayMC = dynamic_cast<TClonesArray*>(aodEvent->GetList()->FindObject(AliAODMCParticle::StdBranchName()));
  
  for (Int_t iPhot = 0; iPhot < convGamma->GetEntriesFast(); iPhot++) {
    //AliAODPWG4Particle * photon = dynamic_cast<AliAODPWG4Particle*>(convGamma->At(iPhot));
    
    
    //if(!photon) {
      AliGammaConversionAODObject * aodO = dynamic_cast<AliGammaConversionAODObject*>(convGamma->At(iPhot));
      if (!aodO) {

	AliError(Form("ERROR: Could not receive ga %d\n", iPhot));
	continue;
      }
      
      AliAODPWG4ParticleCorrelation * photon = AddToAOD(aodO, fAODPWG4Particles, "ConvGamma");
      Int_t tag = CheckTag(photon, tracks, arrayMC);
      if(tag > 0 ) photon->SetTag(tag);
     
  }
}


//////_________________________________________________________________________________________
Int_t AliAnalysisTaskGCPartToPWG4Part::CheckTag(AliAODPWG4ParticleCorrelation * particle, TClonesArray * tracks, TClonesArray * arrayMC) {

  Int_t tag = 0;

  Int_t l1 = particle->GetTrackLabel(0);
  Int_t l2 = particle->GetTrackLabel(1);
  
  AliAODTrack * track1 = NULL;
  AliAODTrack * track2 = NULL;


  for(int i = 0; i < tracks->GetEntriesFast(); i++) {
    
    AliAODTrack * track = (AliAODTrack*)tracks->At(i);
    if (track->GetID() == l1) {
      track1 = track;
    } else if (track->GetID() == l2) {
      track2 = track; 
    }
    
    if(track1 && track2) break;
  }


   
  if(track1->GetLabel() < 0 || track2->GetLabel() < 0) {
    //cout << "error balla"<< endl; 
  
  } else { 

    AliAODMCParticle * mcPart1 = dynamic_cast<AliAODMCParticle*>(arrayMC->At(track1->GetLabel()));
    AliAODMCParticle * mcPart2 = dynamic_cast<AliAODMCParticle*>(arrayMC->At(track2->GetLabel()));
    
    if (mcPart1 && mcPart2) {
      
      //if(mcPart1)
      // cout << mcPart1->GetMother()<< " " <<mcPart2->GetMother() << endl;
      //cout << mcPart1->GetPdgCode()<< " " << mcPart2->GetPdgCode() << endl;
      
      if(mcPart1->GetMother() == mcPart2->GetMother()) {
	//cout << mcPart1->GetMother()->GetLabel();
	
	Int_t motherIndex = mcPart1->GetMother();
	AliAODMCParticle * mother = dynamic_cast<AliAODMCParticle*>(arrayMC->At(motherIndex));
	// cout << "yeay  " <<  mother->GetPdgCode() << endl;
	//  cout << mother->IsPrimary();
	
	
	tag = fAnaUtils->CheckOriginInAOD(&motherIndex, 1, arrayMC);
	
	if(fAnaUtils->CheckTagBit(tag, AliMCAnalysisUtils::kMCPrompt)) {
	  cout << tag << endl;
	  cout << "prim : " << mother->IsPrimary();
	  cout << " " << mother->GetLabel() << endl;
	} 

	if(! mother->IsPrimary()) {
	  //AliAODMCParticle * gp = dynamic_cast<AliAODMCParticle*>(arrayMC->At(mother->GetMother()));
	  //cout << gp->GetPdgCode()  << endl;
	} else if (motherIndex < 10) {
	  cout << "MI: " << motherIndex << endl;
	}

      }
    }
  }
  //cout << "REturn tag " << tag << endl;
  return tag;
}


///__________________________________________________________________________________
AliAODPWG4ParticleCorrelation * AliAnalysisTaskGCPartToPWG4Part::AddToAOD(AliGammaConversionAODObject * aodO, TClonesArray * branch, TString detector) {
  new((*branch)[branch->GetEntriesFast()]) AliAODPWG4ParticleCorrelation(aodO->Px(), aodO->Py(), aodO->Pz(), aodO->E());
  AliAODPWG4ParticleCorrelation * photon = dynamic_cast<AliAODPWG4ParticleCorrelation*>(branch->Last());
  photon->SetTagged(aodO->IsTagged());
  photon->SetTrackLabel(aodO->GetLabel1(), aodO->GetLabel2());
  photon->SetDetector(detector);
  return photon;
}

///__________________________________________________________________________________
AliAODPWG4ParticleCorrelation * AliAnalysisTaskGCPartToPWG4Part::AddToAOD(AliAODConversionParticle * aodO, TClonesArray * branch, TString detector) {
  new((*branch)[branch->GetEntriesFast()]) AliAODPWG4ParticleCorrelation(aodO->Px(), aodO->Py(), aodO->Pz(), aodO->E());
  AliAODPWG4ParticleCorrelation * photon = dynamic_cast<AliAODPWG4ParticleCorrelation*>(branch->Last());
  photon->SetTrackLabel(aodO->GetLabel1(), aodO->GetLabel2());
  photon->SetDetector(detector);
  return photon;
}




//_____________________________________________________________________
void AliAnalysisTaskGCPartToPWG4Part::Terminate(Option_t *) {
  // Draw result to the screen
  // Called once at the end of the query
}

//_____________________________________________________________________
AliAODEvent * AliAnalysisTaskGCPartToPWG4Part::GetAODEvent() {
  //Get the AOD event from whereever it might be
  AliAODEvent * aodEvent = dynamic_cast<AliAODEvent*>(InputEvent());
  if(!aodEvent) {
    aodEvent = AODEvent();
  }
  
  return aodEvent;

}

//_____________________________________________________________________
TClonesArray * AliAnalysisTaskGCPartToPWG4Part::GetConversionGammas(const AliAODEvent * aodEvent) {

  //Get Conversion gamma branch of AOD. First try standard AOD
  TClonesArray * convGamma = dynamic_cast<TClonesArray*>(aodEvent->FindListObject(fAODBranchName.Data()));
  

  //If it's there, send it back
  if(convGamma)  return convGamma;


  //If AOD not in standard file have to locate it in delta AOD
  if( !(fDeltaAODFileName.Length() > 0)  ) return NULL;
  
  AliAODHandler * aodHandler = dynamic_cast<AliAODHandler*>(AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler()); 
  if(aodHandler) {
    AliAODExtension * gExt = dynamic_cast<AliAODExtension*>(aodHandler->GetExtensions()->FindObject(fDeltaAODFileName));
    if(gExt) {
      AliAODEvent * gcEvent = gExt->GetAOD();
      return dynamic_cast<TClonesArray*>(gcEvent->FindListObject("GammaConv_gamma"));
    }
  }  
  return NULL;
}

//_________________________________________________________________________
void AliAnalysisTaskGCPartToPWG4Part::CleanUp() {
  fAODPWG4Particles->Delete();
}


