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

// Gamma - jet correlation analysis task
// Authors: Svein Lindal


using namespace std;

ClassImp(AliAnalysisTaskGCPartToPWG4Part)

//________________________________________________________________________
AliAnalysisTaskGCPartToPWG4Part::AliAnalysisTaskGCPartToPWG4Part() 
: AliAnalysisTaskSE(), 
  fDeltaAODFileName("AliAODConversionGamma.root"),
  fAODBranchName("ConvGamma_gamma"),
  fAODPWG4Particles(NULL)
{
  // Dummy Constructor
}

//________________________________________________________________________________
AliAnalysisTaskGCPartToPWG4Part::~AliAnalysisTaskGCPartToPWG4Part() {

  if(fAODPWG4Particles)
    fAODPWG4Particles = NULL;
  delete fAODPWG4Particles;

}



//________________________________________________________________________
AliAnalysisTaskGCPartToPWG4Part::AliAnalysisTaskGCPartToPWG4Part(const char *name) : 
  AliAnalysisTaskSE(name), 
  fDeltaAODFileName("AliAODConversionGamma.root"),
  fAODBranchName("ConvGamma_gamma"),
  fAODPWG4Particles(NULL)
{
  // Constructor
  DefineInput(0, TChain::Class());
  // Output slot #0 id reserved by the base class for AOD

  // Output slot #1 writes into a TH1 container
  DefineOutput(1, TList::Class());
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
  


  for (Int_t iPhot = 0; iPhot < convGamma->GetEntriesFast(); iPhot++) {
    AliAODPWG4Particle * photon = dynamic_cast<AliAODPWG4Particle*>(convGamma->At(iPhot));
    
    
    if(!photon) {
      AliGammaConversionAODObject * aodO = dynamic_cast<AliGammaConversionAODObject*>(convGamma->At(iPhot));
      if (!aodO) {
	AliError(Form("ERROR: Could not receive ga %d\n", iPhot));
	continue;
      }
      
      //if(aodO->Pt() < GetMinPt()) continue;
      AddToAOD(aodO, fAODPWG4Particles, "ConvGamma");

    }
  }
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


