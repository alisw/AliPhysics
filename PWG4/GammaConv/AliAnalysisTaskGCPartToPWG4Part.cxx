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
#include "AliAODConversionPhoton.h"
#include "AliAODJet.h"

#include "AliAODInputHandler.h"

#include "AliAODMCParticle.h"


#include "AliAODMCHeader.h"
// Gamma - jet correlation analysis task
// Authors: Svein Lindal


using namespace std;

ClassImp(AliAnalysisTaskGCPartToPWG4Part)

//________________________________________________________________________
AliAnalysisTaskGCPartToPWG4Part::AliAnalysisTaskGCPartToPWG4Part() 
: AliAnalysisTaskSE(), 
  fDeltaAODFileName(""),
  fGammaCutString("GammaConv"),
  fPionCutString("GammaConv"),
  fAODBranchName("GammaConv_gamma"),
  fAODPWG4Photons(NULL),
  fAODPWG4Pi0(NULL),
   fDebugLevel(0)
{
  // Dummy Constructor
}

//________________________________________________________________________________
AliAnalysisTaskGCPartToPWG4Part::~AliAnalysisTaskGCPartToPWG4Part() {

  if(fAODPWG4Photons)
    delete fAODPWG4Photons;
  fAODPWG4Photons = NULL;

  if(fAODPWG4Pi0)
    delete fAODPWG4Pi0;
  fAODPWG4Pi0 = NULL;
 

}



//________________________________________________________________________
AliAnalysisTaskGCPartToPWG4Part::AliAnalysisTaskGCPartToPWG4Part(const char *name) : 
  AliAnalysisTaskSE(name), 
  fDeltaAODFileName(""),
  fGammaCutString("GammaConv"),
  fPionCutString("GammaConv"),
  fAODBranchName("GammaConv_gamma"),
  fAODPWG4Photons(NULL),
  fAODPWG4Pi0(NULL),
  fDebugLevel(0)
{
  // Constructor
  DefineInput(0, TChain::Class());
  // Output slot #0 id reserved by the base class for AOD

  // Output slot #1 writes into a TH1 container
  DefineOutput(1, TList::Class());

}



//________________________________________________________________________
void AliAnalysisTaskGCPartToPWG4Part::UserCreateOutputObjects() {
  fAODPWG4Photons = new TClonesArray("AliAODPWG4ParticleCorrelation", 0);
  fAODPWG4Photons->SetName("PhotonsCTS");
  AddAODBranch("TClonesArray", &fAODPWG4Photons);

  fAODPWG4Pi0 = new TClonesArray("AliAODPWG4ParticleCorrelation", 0);
  fAODPWG4Pi0->SetName("Pi0sCTS");
  AddAODBranch("TClonesArray", &fAODPWG4Pi0);

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


  //TClonesArray * arrayMC = dynamic_cast<TClonesArray*>(aodEvent->GetList()->FindObject(AliAODMCParticle::StdBranchName()));  
  for (Int_t iPhot = 0; iPhot < convGamma->GetEntriesFast(); iPhot++) {

    AliAODPWG4ParticleCorrelation * photon = NULL;
    AliAODConversionPhoton * convParticle = dynamic_cast<AliAODConversionPhoton*>(convGamma->At(iPhot));
    if (convParticle && BothTracksPresent(convParticle, tracks)) {
      photon = AddToAOD(convParticle, fAODPWG4Photons, "ConvGamma");
      
    } else {
      continue;
    }
    
    if(photon && fDebugLevel > 2) {
      printf("Added conversion photon number %d, pt: %f \n", iPhot, photon->Pt());
    }

  }


  TClonesArray * pions = GetPions(aodEvent);
  if(!pions) {
    AliError(Form("No branch by name %s found in file %s", fAODBranchName.Data(), fDeltaAODFileName.Data()));
    return;
  }

  for (Int_t iPhot = 0; iPhot < pions->GetEntriesFast(); iPhot++) {
    AliAODPWG4ParticleCorrelation * pion = NULL;
    AliAODConversionPhoton * convParticle = dynamic_cast<AliAODConversionPhoton*>(pions->At(iPhot));
    if (convParticle && BothGammaPresent(convParticle, convGamma, tracks)) {
      pion = AddPionToAOD(convParticle, fAODPWG4Pi0, "ConvGamma", convGamma);
      
    } else {
      continue;
    }
    
    if(pion && fDebugLevel > 2) {
      printf("Added conversion pion number %d, pt: %f \n", iPhot, pion->Pt());
    }
  }

  

}




///__________________________________________________________________________________
AliAODPWG4ParticleCorrelation * AliAnalysisTaskGCPartToPWG4Part::AddToAOD(AliGammaConversionAODObject * aodO, TClonesArray * branch, TString detector) {
  new((*branch)[branch->GetEntriesFast()]) AliAODPWG4ParticleCorrelation(aodO->Px(), aodO->Py(), aodO->Pz(), aodO->E());
  AliAODPWG4ParticleCorrelation * photon = dynamic_cast<AliAODPWG4ParticleCorrelation*>(branch->Last());
  if(photon) {
    photon->SetTagged(aodO->IsTagged());
    photon->SetTrackLabel(aodO->GetLabel1(), aodO->GetLabel2());
    photon->SetDetector(detector);
    return photon;
  } else {
    return NULL;
  }
}

///__________________________________________________________________________________
AliAODPWG4ParticleCorrelation * AliAnalysisTaskGCPartToPWG4Part::AddToAOD(AliAODConversionPhoton * aodO, TClonesArray * branch, TString detector) {
  new((*branch)[branch->GetEntriesFast()]) AliAODPWG4ParticleCorrelation(aodO->Px(), aodO->Py(), aodO->Pz(), aodO->E());
  AliAODPWG4ParticleCorrelation * photon = dynamic_cast<AliAODPWG4ParticleCorrelation*>(branch->Last());
  if(photon) {
    photon->SetTrackLabel(aodO->GetLabel1(), aodO->GetLabel2());
    photon->SetDetector(detector);
    return photon;
  } else {
    return NULL;
  }
  
}


///__________________________________________________________________________________
AliAODPWG4ParticleCorrelation * AliAnalysisTaskGCPartToPWG4Part::AddPionToAOD(AliAODConversionPhoton * pion, TClonesArray * branch, TString detector, TClonesArray * photons) {
  new((*branch)[branch->GetEntriesFast()]) AliAODPWG4ParticleCorrelation(pion->Px(), pion->Py(), pion->Pz(), pion->E());
  AliAODPWG4ParticleCorrelation * pwg4Pion = dynamic_cast<AliAODPWG4ParticleCorrelation*>(branch->Last());
  if(pwg4Pion) {
    Int_t tl[4] = {-1, -1, -1, -1};
    //pion->GetGrandChildren(photons, tl);
    pwg4Pion->SetTrackLabel(tl[0], tl[1], tl[2], tl[3]);
    pwg4Pion->SetDetector(detector);
    for(Int_t i = 0; i < 4; i++) {
      cout << tl[i] << " ";
    }
    cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
    return pwg4Pion;
  } else {
    return NULL;
  }
  
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
TClonesArray * AliAnalysisTaskGCPartToPWG4Part::GetAODBranch(const AliAODEvent * aodEvent, TString branchName) const {

  //Get Conversion gamma branch of AOD. First try standard AOD
  TClonesArray * branch = dynamic_cast<TClonesArray*>(aodEvent->FindListObject(branchName.Data()));
  

  //If it's there, send it back
  if(branch)  return branch;
  

  //If AOD not in standard file have to locate it in delta AOD
  if( !(fDeltaAODFileName.Length() > 0)  ) return NULL;
  
  AliAODHandler * aodHandler = dynamic_cast<AliAODHandler*>(AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler()); 
  if(aodHandler) {
    AliAODExtension * gExt = dynamic_cast<AliAODExtension*>(aodHandler->GetExtensions()->FindObject(fDeltaAODFileName));
    if(gExt) {
      AliAODEvent * gcEvent = gExt->GetAOD();
      return dynamic_cast<TClonesArray*>(gcEvent->FindListObject(branchName.Data()));
    }
  }  
  return NULL;

}
//_____________________________________________________________________
TClonesArray * AliAnalysisTaskGCPartToPWG4Part::GetConversionGammas(const AliAODEvent * aodEvent) const {
  return GetAODBranch(aodEvent, Form("%s_gamma", fGammaCutString.Data()));
 
}

//_____________________________________________________________________
TClonesArray * AliAnalysisTaskGCPartToPWG4Part::GetPions(const AliAODEvent * aodEvent) const {
    return GetAODBranch(aodEvent, Form("%s_Pi0", fPionCutString.Data()));
 
}


//_________________________________________________________________________
void AliAnalysisTaskGCPartToPWG4Part::CleanUp() {
  fAODPWG4Photons->Delete();
  fAODPWG4Pi0->Delete();
}


//______________________________________________________________________________________________
Bool_t AliAnalysisTaskGCPartToPWG4Part::BothTracksPresent(const AliAODConversionPhoton * const photon, const TClonesArray * const tracks)  const {

  AliAODTrack * track1 = NULL;
  AliAODTrack * track2 = NULL;
  for(Int_t i = 0; i < tracks->GetEntriesFast(); i++) {
    AliAODTrack * track = dynamic_cast<AliAODTrack*>(tracks->At(i));
    if(track) {
      if(track->GetID() == photon->GetLabel1()) track1 = track;
      else if (track->GetID() == photon->GetLabel2()) track2 = track;
      if(track1 && track2) break;
    }
  }
  
  if(track1 && track2) {
    return kTRUE;
  }
  cout << "Could not get both tracks!!! labels "  << photon->GetLabel1() << " " << photon->GetLabel2()  <<endl;
  return kFALSE;
  

}

//______________________________________________________________________________________________
Bool_t AliAnalysisTaskGCPartToPWG4Part::BothGammaPresent(const AliAODConversionPhoton * const pion, const TClonesArray * const photons, const TClonesArray * const tracks)  const {

  AliAODConversionPhoton * photon1 = dynamic_cast<AliAODConversionPhoton*>(photons->At(pion->GetLabel1()));
  AliAODConversionPhoton * photon2 = dynamic_cast<AliAODConversionPhoton*>(photons->At(pion->GetLabel2()));

  if(photon1 && photon2) {
    if( BothTracksPresent(photon1, tracks) &&  BothTracksPresent(photon1, tracks)) {
      return kTRUE;
    }
  } else {
    cout << "can't find both photons "<< endl;
  }
  
  return kFALSE;
}
