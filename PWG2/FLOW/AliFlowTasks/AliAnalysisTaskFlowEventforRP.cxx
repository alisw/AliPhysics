/*************************************************************************
* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
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

////////////////////////////////////////////////////
// AliAnalysisTaskFlowEventforRP:
//
// analysis task for filling the flow event
// from MCEvent, ESD
// and put it in an output stream so the calculated
// Reaction Plane can be stored in the AODHeader
// when the AOD is made from the ESD 
// for cuts the correction framework is used
// which also outputs QA histograms to view
// the effects of the cuts
////////////////////////////////////////////////////

#include "Riostream.h" //needed as include
#include "TChain.h"
#include "TTree.h"
#include "TFile.h" //needed as include
#include "TList.h"
#include "TRandom3.h"
#include "TTimeStamp.h"

// ALICE Analysis Framework
class AliAnalysisTaskSE;
#include "AliAnalysisManager.h"

// ESD interface
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"

// AOD interface
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"

// Monte Carlo Eventp
#include "AliAODHandler.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"

// ALICE Correction Framework
#include "AliCFManager.h"

// Interface to Event generators to get Reaction Plane Angle
#include "AliGenCocktailEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliGenGeVSimEventHeader.h"

// Interface to make the Flow Event Simple used in the flow analysis methods
#include "AliFlowEventSimpleMaker.h"
#include "AliFlowVector.h"
#include "AliAnalysisTaskFlowEventforRP.h"


ClassImp(AliAnalysisTaskFlowEventforRP)
  
//________________________________________________________________________
AliAnalysisTaskFlowEventforRP::AliAnalysisTaskFlowEventforRP(const char *name) : 
  AliAnalysisTaskSE(name), 
//  fOutputFile(NULL),
  //fESD(NULL),
  //fAOD(NULL),
  fEventMaker(NULL),
  fAnalysisType("ESD"),
  fCFManager1(NULL),
  fCFManager2(NULL),
  fMinMult(0),
  fMaxMult(10000000),
  fMCReactionPlaneAngle(0.)
  
  
{
  // Constructor
  cout<<"AliAnalysisTaskFlowEventforRP::AliAnalysisTaskFlowEventforRP(const char *name)"<<endl;

  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Define here the flow event output
  DefineOutput(0, AliFlowEventSimple::Class());  
  
}

//________________________________________________________________________
AliAnalysisTaskFlowEventforRP::AliAnalysisTaskFlowEventforRP() : 
  //  fOutputFile(NULL),
  //fESD(NULL),
  //fAOD(NULL),
  fEventMaker(NULL),
  fAnalysisType("ESD"),
  fCFManager1(NULL),
  fCFManager2(NULL),
  fMinMult(0),
  fMaxMult(10000000),
  fMCReactionPlaneAngle(0.)
{
  // Constructor
  cout<<"AliAnalysisTaskFlowEventforRP::AliAnalysisTaskFlowEventforRP()"<<endl;
}


//________________________________________________________________________
AliAnalysisTaskFlowEventforRP::~AliAnalysisTaskFlowEventforRP()
{
  //
  // Destructor
  //
  
}


//________________________________________________________________________
void AliAnalysisTaskFlowEventforRP::UserCreateOutputObjects() 
{
  // Called at every worker node to initialize
  cout<<"AliAnalysisTaskFlowEventforRP::UserCreateOutputObjects()"<<endl;

  if (!(fAnalysisType == "ESD")) {
    cout<<"WRONG ANALYSIS TYPE! only ESD for this method."<<endl;
    exit(1);
  }

  // Flow Event maker
  fEventMaker = new AliFlowEventSimpleMaker();
  
}

//________________________________________________________________________
void AliAnalysisTaskFlowEventforRP::UserExec(Option_t *) 
{
  // Main loop

  AliESDEvent* esd = dynamic_cast<AliESDEvent*>(InputEvent());
  //AliESD*      old = esd->GetAliESDOld();

  // Called for each event
  AliFlowEventSimple* fEvent = NULL;
  Double_t fRP = 0.; // the monte carlo reaction plane angle
  
  AliMCEvent* mcEvent = NULL;
  // See if we can get Monte Carlo Information and if so get the reaction plane
  AliMCEventHandler* eventHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
  if (eventHandler) {
  mcEvent = eventHandler->MCEvent();
  if (mcEvent) {
    //COCKTAIL with HIJING
    if (!strcmp(mcEvent-> GenEventHeader()->GetName(),"Cocktail Header")) { //returns 0 if matches
      AliGenCocktailEventHeader *headerC = dynamic_cast<AliGenCocktailEventHeader *> (mcEvent-> GenEventHeader()); 
      if (headerC) {
	TList *lhd = headerC->GetHeaders();
	if (lhd) {
	  AliGenHijingEventHeader *hdh = dynamic_cast<AliGenHijingEventHeader *> (lhd->At(0)); 
	  if (hdh) {
	    fRP = hdh->ReactionPlaneAngle();
	    //cout<<"The reactionPlane from Hijing (Cocktail) is: "<< fRP <<endl;
	  }
	}
      }
      //else { cout<<"headerC is NULL"<<endl; }
    }
    //GEVSIM
    else if (!strcmp(mcEvent-> GenEventHeader()->GetName(),"GeVSim header")) { //returns 0 if matches
      AliGenGeVSimEventHeader* headerG = (AliGenGeVSimEventHeader*)(mcEvent->GenEventHeader());
      if (headerG) {
	fRP = headerG->GetEventPlane();
	//cout<<"The reactionPlane from GeVSim is: "<< fRP <<endl;
      }
      //else { cout<<"headerG is NULL"<<endl; }
    }
    //HIJING
    else if (!strcmp(mcEvent-> GenEventHeader()->GetName(),"Hijing")) { //returns 0 if matches
      AliGenHijingEventHeader* headerH = (AliGenHijingEventHeader*)(mcEvent->GenEventHeader());
      if (headerH) {
	fRP = headerH->ReactionPlaneAngle();
	//cout<<"The reactionPlane from Hijing is: "<< fRP <<endl;
      }
      //else { cout<<"headerH is NULL"<<endl; }
    }
  }
  else {cout<<"No MC event!"<<endl; }
  }
  else {cout<<"No eventHandler!"<<endl; }

  fEventMaker->SetMCReactionPlaneAngle(fRP);
  
  //setting event cuts
  fEventMaker->SetMinMult(fMinMult);
  fEventMaker->SetMaxMult(fMaxMult);
  
  // Fill the FlowEventSimple for ESD input  
  //else if (fAnalysisType == "ESD") {
  if (fAnalysisType == "ESD") {
    if (!fCFManager1) {cout << "ERROR: No pointer to correction framework cuts! " << endl; return; }
    if (!fCFManager2) {cout << "ERROR: No pointer to correction framework cuts! " << endl; return; }
    
    if (!esd) { Printf("ERROR: esd not available"); return;}
    Printf("There are %d tracks in this event", esd->GetNumberOfTracks());
    
    // analysis
    fEvent = fEventMaker->FillTracks(esd,fCFManager1,fCFManager2);
        
    AliFlowVector vQ = fEvent->GetQ();                      
    Double_t dRP[1] = {0.0};                      // Phi is een Double_t, maar SetQTheta heeft een Double_t* nodig, dus een double in array vorm. 
    dRP[0] = vQ.Phi()/2; 
      
    cout<<"The reaction plane from MC is "<<fRP<<endl;
    cout<<"The calculated reaction plane is "<<dRP[0]<<endl;

    
    // Update the header
    AliAODHeader* header = AODEvent()->GetHeader();
    header->SetRunNumber(esd->GetRunNumber());
    header->SetQTheta(dRP,1);
        
  }
    
  PostData(0,fEvent);
  
} 

//________________________________________________________________________
void AliAnalysisTaskFlowEventforRP::Terminate(Option_t *) 
{
  // Called once at the end of the query -- do not call in case of CAF
  
}


