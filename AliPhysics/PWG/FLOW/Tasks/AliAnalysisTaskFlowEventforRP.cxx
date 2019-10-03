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
#include "AliAnalysisTaskSE.h"
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
#include "AliGenEposEventHeader.h"

// Interface to make the Flow Event Simple used in the flow analysis methods
#include "AliFlowEvent.h"
#include "AliFlowVector.h"
#include "AliAnalysisTaskFlowEventforRP.h"

using std::cout;
using std::endl;
ClassImp(AliAnalysisTaskFlowEventforRP)
  
//________________________________________________________________________
AliAnalysisTaskFlowEventforRP::AliAnalysisTaskFlowEventforRP(const char *name) : 
  AliAnalysisTaskSE(name), 
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
    
}

//________________________________________________________________________
void AliAnalysisTaskFlowEventforRP::UserExec(Option_t *) 
{
  // Main loop

  AliESDEvent* esd = dynamic_cast<AliESDEvent*>(InputEvent());
  AliFlowEvent* fEvent = NULL;
  AliMCEvent* mcEvent  = MCEvent();
  Double_t fRP = 0.; // the monte carlo reaction plane angle
  
  // Fill the FlowEventSimple for ESD input  
  if (fAnalysisType == "ESD") {
    if (!(fCFManager1&&fCFManager2))
      {
	cout << "ERROR: No pointer to correction framework cuts! " << endl; 
	return; 
      }
    if (!esd)
      {
	AliError("ERROR: ESD not available");
	return;
      }
    
    //check the offline trigger (check if the event has the correct trigger)
    //AliInfo(Form("ESD has %d tracks", fInputEvent->GetNumberOfTracks()));
    
    //check multiplicity
    if (!fCFManager1->CheckEventCuts(AliCFManager::kEvtRecCuts,esd))
      {
	cout << "Event does not pass multiplicity cuts" << endl;
	return;
      }
    
    // make the flowevent
    fEvent = new AliFlowEvent(esd,fCFManager1,fCFManager2);
    
    if (mcEvent && mcEvent->GenEventHeader()) 
      {
	fEvent->SetMCReactionPlaneAngle(mcEvent);
	fRP = fEvent->GetMCReactionPlaneAngle();
      }
    
    //check final event cuts
    Int_t mult = fEvent->NumberOfTracks();
    cout << "FlowEvent has "<<mult<<" tracks"<<endl;
    if (mult<fMinMult || mult>fMaxMult)
      {
	cout << "FlowEvent cut on multiplicity" << endl;
	return;
      }

  
    // get the flow vector     
    AliFlowVector vQ = fEvent->GetQ();                      
    Double_t dRP[1] = {0.0};   
    // Phi is a Double_t, but SetQTheta() needs as input Double_t*, 
    // an array of doubles. 
    dRP[0] = vQ.Phi()/2; 
      
    cout<<"The reaction plane from MC is "<<fRP<<endl;
    cout<<"The calculated reaction plane is "<<dRP[0]<<endl;

    
    // Update the header
    AliAODHeader* header = dynamic_cast<AliAODHeader*>(AODEvent()->GetHeader());
    if(!header) AliFatal("Not a standard AOD");
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


