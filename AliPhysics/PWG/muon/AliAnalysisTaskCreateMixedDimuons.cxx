#include "TChain.h"
#include "TTree.h"

#include "AliAnalysisTaskME.h"
#include "AliAnalysisManager.h"

#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAODVertex.h"
#include "AliMultiEventInputHandler.h"
#include "AliAODHandler.h"

#include "AliAnalysisTaskCreateMixedDimuons.h"
#include "AliEventPoolMuon.h"
#include "TDatabasePDG.h"
#include "TRandom.h"

// Example of an analysis task creating aod events filled with mixed muon pairs
//
// Authors Alessandro De Falco and Antonio Uras, INFN Cagliari
// alessandro.de.falco@ca.infn.it  antonio.uras@ca.infn.it

#define AliAnalysisTaskCreateMixedDimuons_CXX

ClassImp(AliAnalysisTaskCreateMixedDimuons)

//=================================================================================

AliAnalysisTaskCreateMixedDimuons::AliAnalysisTaskCreateMixedDimuons(const char *name) 
: AliAnalysisTaskME(name),
  fBufferSize(0),
  fOutputUserHandler(0x0),
  fOutputUserAOD(0X0),
  fOutputUserAODTree(0X0),
  fPoolMuon(0X0),
  fDebug(0X0)
{

  // Constructor

  // Default input and output containers
  DefineInput (0, TChain::Class());

  // User-defined input and output containers
  DefineOutput(1, TTree::Class());

  // ---------------------------------------------------------------

  fDebug = kFALSE;
  fBufferSize = 0;

  RequireFreshBuffer();

  for (Int_t i=0; i<100; i++) fInputAOD[i] = 0;

  fOutputUserHandler = 0;
  fOutputUserAOD     = 0;
  fOutputUserAODTree = 0;  

}

//=================================================================================

void AliAnalysisTaskCreateMixedDimuons::ConnectInputData(Option_t *) {

  // Connect ESD or AOD here
  // Called once

  printf("-> AliAnalysisTaskCreateMixedDimuons::ConnectInputData\n");

  TTree* tree = (TTree*) GetInputData(0);

  if (!tree) {
    Printf("ERROR: Could not read chain from input slot 0");
  } 
  else {

    fInputHandler = (AliMultiEventInputHandler*) AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler();
    fPoolMuon     = (AliEventPoolMuon*)          AliAnalysisManager::GetAnalysisManager()->GetEventPool();

    if (!fInputHandler) {
      Printf("ERROR: Could not get AliMultiAODInputHandler");
    }
    else {  
      fBufferSize = fInputHandler->GetBufferSize();
      if (fBufferSize>100) {
        printf("\n*** WARNING AliAnalysisTaskCreateMixedDimuons::ConnectInputData -> Trying to set fBufferSize>100, forcing fBufferSize=100 ***\n\n");
        fBufferSize = 100;
      }
      for (Int_t i=0; i<fBufferSize; i++) fInputAOD[i] = (AliAODEvent*) fInputHandler->GetEvent(i);
    }
  }

  printf("<- AliAnalysisTaskCreateMixedDimuons::ConnectInputData\n");

}

//=================================================================================

void AliAnalysisTaskCreateMixedDimuons::UserCreateOutputObjects() {

  // Here the user-defined output containers should be created!!!
  // Called once

  fOutputUserHandler = new AliAODHandler();

  fOutputUserHandler -> Init("");

  fOutputUserAOD = fOutputUserHandler -> GetAOD();

  fOutputUserAODTree = fOutputUserHandler -> GetTree();

}

//=================================================================================

void AliAnalysisTaskCreateMixedDimuons::UserExec(Option_t *) {

  if (!fOutputUserAOD) {
    Printf("ERROR: fOutputUserAOD not available\n");
    return;
  }

  printf("Calling USER EXEC\n\n");

  for (Int_t iEv=0; iEv<fBufferSize; iEv++) {
    for (Int_t jEv=0; jEv<iEv; jEv++) {

      Int_t nTracksEv[2]  = {0};
      Int_t nFWMuonsEv[2] = {0};
      
      nTracksEv[0] = fInputAOD[iEv]->GetNumberOfTracks();
      nTracksEv[1] = fInputAOD[jEv]->GetNumberOfTracks();
      // Check if both events contain std aod tracks
      if (!dynamic_cast<AliAODTrack*>(fInputAOD[iEv]->GetTrack(0)) || dynamic_cast<AliAODTrack*>(fInputAOD[jEv]->GetTrack(0))){
        AliFatal("Not a standard AOD");
      }

      for (Int_t i=0; i<nTracksEv[0]; i++) if(((AliAODTrack*)fInputAOD[iEv]->GetTrack(i))->IsMuonTrack()) nFWMuonsEv[0]++;
      for (Int_t i=0; i<nTracksEv[1]; i++) if(((AliAODTrack*)fInputAOD[jEv]->GetTrack(i))->IsMuonTrack()) nFWMuonsEv[1]++;
      
      // Muon track mixing to fill a mass spectrum
      
      if (nFWMuonsEv[0] && nFWMuonsEv[1]) {
	
	Int_t rndMuonTrack[2] = {0};
	rndMuonTrack[0] = gRandom->Integer(nFWMuonsEv[0]);
	rndMuonTrack[1] = gRandom->Integer(nFWMuonsEv[1]);
	
	Int_t nFWMUonsAdded = 0;
	Int_t nPosTracksAdded = 0;
	Int_t nNegTracksAdded = 0;
	
	AliAODVertex *vertex = new AliAODVertex();
	vertex -> SetX(0.0);
	vertex -> SetY(0.0);
	vertex -> SetZ(fPoolMuon->GetMeanPrimaryVertexZ());
	
	Int_t muonCounter[2] = {0};
	
	// adding tracks and vertex to the output event...
	
	for (Int_t i=0; i<nTracksEv[0]; i++) {
	  if(((AliAODTrack*)fInputAOD[iEv]->GetTrack(i))->IsMuonTrack()) {
	    if (fDebug) printf("fInputAOD[%d]->GetTrack(%d) = %p    pt = %f     uniqueID = %d\n",
			       iEv,i,fInputAOD[iEv]->GetTrack(i),fInputAOD[iEv]->GetTrack(i)->Pt(),
			       fInputAOD[iEv]->GetTrack(i)->GetUniqueID());
	    if (muonCounter[0]==rndMuonTrack[0]) {
	      fOutputUserAOD->AddTrack((AliAODTrack*)fInputAOD[iEv]->GetTrack(i));
	      nFWMUonsAdded++;
	      if (fInputAOD[iEv]->GetTrack(i)->Charge()>0) nPosTracksAdded++;
	      else nNegTracksAdded++;
	    }
	    muonCounter[0]++;
	  }
	}
	
	for (Int_t i=0; i<nTracksEv[1]; i++) {
	  if(((AliAODTrack*)fInputAOD[jEv]->GetTrack(i))->IsMuonTrack()) {
	    if (fDebug) printf("fInputAOD[%d]->GetTrack(%d) = %p    pt = %f     uniqueID = %d\n",
			       jEv,i,fInputAOD[jEv]->GetTrack(i),fInputAOD[jEv]->GetTrack(i)->Pt(),
			       fInputAOD[jEv]->GetTrack(i)->GetUniqueID());
	    if (muonCounter[1]==rndMuonTrack[1]) {
	      fOutputUserAOD->AddTrack((AliAODTrack*)fInputAOD[jEv]->GetTrack(i));
	      nFWMUonsAdded++;
	      if (fInputAOD[jEv]->GetTrack(i)->Charge()>0) nPosTracksAdded++;
	      else nNegTracksAdded++;
	    }
	    muonCounter[1]++;	
	  }
	}
	
	fOutputUserAOD->AddVertex(vertex);
	
	// ... done!
	
	if (fDebug) {
	  for (Int_t i=0; i<nFWMUonsAdded; i++) {
	    AliAODTrack *tr = dynamic_cast<AliAODTrack*>( fOutputUserAOD->GetTrack(i));
	    if(!tr) AliFatal("Not a standard AOD");
	    printf("fOutputUserAOD->GetTrack(%d) = %p    pt = %f\n",i,tr,tr->Pt());
	  }
	}
	

        AliAODHeader * header = dynamic_cast<AliAODHeader*>(fOutputUserAOD->GetHeader());
        if(!header) AliFatal("Not a standard AOD");
	header->SetRefMultiplicity(nFWMUonsAdded); 
	header->SetRefMultiplicityPos(nPosTracksAdded);
	header->SetRefMultiplicityNeg(nNegTracksAdded);
	
	fOutputUserHandler -> FinishEvent();
	
      }

      PostData(1, fOutputUserAODTree);
      
    }
    
  }
  
}      

//===================================================================================

void AliAnalysisTaskCreateMixedDimuons::Terminate(Option_t *) {

  // Called once at the end of the query

  printf("\n\nCalling TERMINATE \n\n\n");

  fOutputUserHandler -> Terminate();

}

//===================================================================================
