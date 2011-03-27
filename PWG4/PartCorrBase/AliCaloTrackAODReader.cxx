
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
/* $Id: $ */

//_________________________________________________________________________
// Class for reading data (AODs) in order to do prompt gamma
//  or other particle identification and correlations.
// This part is commented: Mixing analysis can be done, input AOD with events
// is opened in the AliCaloTrackReader::Init()
// 
//
//*-- Author: Gustavo Conesa (LNF-INFN) 
//////////////////////////////////////////////////////////////////////////////


// --- ROOT system ---
//#include "Riostream.h"

//---- ANALYSIS system ----
#include "AliCaloTrackAODReader.h" 
#include "AliAODInputHandler.h"
#include "AliMultiEventInputHandler.h"
#include "AliAnalysisManager.h"
#include "AliMixedEvent.h"

ClassImp(AliCaloTrackAODReader)

//____________________________________________________________________________
AliCaloTrackAODReader::AliCaloTrackAODReader() : 
  AliCaloTrackReader(), fOrgInputEvent(0x0)
{
  //Default Ctor
  
  //Initialize parameters
  fDataType=kAOD;
  fReadStack          = kTRUE;
  fReadAODMCParticles = kFALSE;
 
}

//____________________________________________________________________________
//void AliCaloTrackAODReader::GetSecondInputAODVertex(Double_t  v[3]) const {
//	//Return vertex position of second AOD input
//	
//	fSecondInputAODEvent->GetPrimaryVertex()->GetXYZ(v);
//
//}

//____________________________________________________________________________
AliCentrality* AliCaloTrackAODReader::GetCentrality() const {
  // recover centrality object.
  AliAODEvent* event    = dynamic_cast<AliAODEvent*> (fInputEvent);
  AliAODEvent* orgevent = dynamic_cast<AliAODEvent*> (fOrgInputEvent);

  if(event && !fSelectEmbeddedClusters) {
    //Normal AOD event
    return event->GetHeader()->GetCentralityP() ;
  }
  else if(fSelectEmbeddedClusters && fOrgInputEvent) {
    // centrality in AOD from input, not in embedded event
    // temporary fix until this object is copied to the output event in embedding analysis
    return orgevent->GetHeader()->GetCentralityP();
  }
  else {
    return 0x0 ; 
  }
}


//____________________________________________________________________________
void AliCaloTrackAODReader::SetInputOutputMCEvent(AliVEvent* input, AliAODEvent* aod, AliMCEvent* mc) {
  // Connect the data pointers
  // If input is AOD, do analysis with input, if not, do analysis with the output aod.

  //printf("AODInputHandler %p, MergeEvents %d \n",aodIH, aodIH->GetMergeEvents());

  Bool_t tesd = kFALSE ; 
  Bool_t taod = kTRUE ; 
  if ( strcmp(input->GetName(), "AliMixedEvent") == 0 ) {
    AliMultiEventInputHandler* multiEH = dynamic_cast<AliMultiEventInputHandler*>((AliAnalysisManager::GetAnalysisManager())->GetInputEventHandler());
    if(multiEH){
      if (multiEH->GetFormat() == 0 ) {
        tesd = kTRUE ; 
      } else if (multiEH->GetFormat() == 1) {
        taod = kTRUE ; 
      }
    }
    else{
      printf("AliCaloTrackAODReader::SetInputOutputMCEvent() - MultiEventHandler is NULL");
      abort();
    }
  }
  if (strcmp(input->GetName(),"AliESDEvent") == 0) {
    tesd = kTRUE ; 
  } else if (strcmp(input->GetName(),"AliAODEvent") == 0) {
    taod = kTRUE ; 
  }
  

  if(tesd)   {
    SetInputEvent(aod);
    SetOutputEvent(aod);
    fOrgInputEvent = input;
  }
  else if(taod){
    AliAODInputHandler* aodIH = dynamic_cast<AliAODInputHandler*>((AliAnalysisManager::GetAnalysisManager())->GetInputEventHandler());
	  if (aodIH && aodIH->GetMergeEvents()) {
		  //Merged events, use output AOD.
		  SetInputEvent(aod);
		  SetOutputEvent(aod);
      fOrgInputEvent = input;
	  }
	  else{
		  SetInputEvent(input);
		  SetOutputEvent(aod);
	  }
  }
  else{ 
    AliFatal(Form("AliCaloTrackAODReader::SetInputOutputMCEvent() - STOP : Wrong data format: %s\n",input->GetName()));
  }
  
  SetMC(mc);
  
}

