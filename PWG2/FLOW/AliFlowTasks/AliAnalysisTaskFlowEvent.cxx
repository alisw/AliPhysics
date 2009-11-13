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
// AliAnalysisTaskFlowEvent:
//
// analysis task for filling the flow event
// from MCEvent, ESD, AOD ....
// and put it in an output stream so it can 
// be used by the various flow analysis methods 
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
class AliAnalysisTask;
#include "AliAnalysisManager.h"

// ESD interface
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"

// AOD interface
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"

// Monte Carlo Event
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
#include "AliFlowEventSimpleMaker.h"

#include "AliAnalysisTaskFlowEvent.h"


ClassImp(AliAnalysisTaskFlowEvent)
  
//________________________________________________________________________
AliAnalysisTaskFlowEvent::AliAnalysisTaskFlowEvent(const char *name, Bool_t on) : 
  AliAnalysisTask(name, ""), 
//  fOutputFile(NULL),
  fESD(NULL),
  fAOD(NULL),
  fEventMaker(NULL),
  fAnalysisType("ESD"),
  fCFManager1(NULL),
  fCFManager2(NULL),
  fQAInt(NULL),
  fQADiff(NULL),
  fMinMult(0),
  fMaxMult(10000000),
  fQA(on),
  fMCReactionPlaneAngle(0.),
  fCount(0),
  fNoOfLoops(1),
  fEllipticFlowValue(0.),
  fSigmaEllipticFlowValue(0.),
  fMultiplicityOfEvent(1000000000),
  fSigmaMultiplicityOfEvent(0),
  fMyTRandom3(NULL)
{
  // Constructor
  cout<<"AliAnalysisTaskFlowEvent::AliAnalysisTaskFlowEvent(const char *name)"<<endl;

  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Define here the flow event output
  DefineOutput(0, AliFlowEventSimple::Class());  
  if(on) {
    DefineOutput(1, TList::Class());
    DefineOutput(2, TList::Class()); }  
  // and for testing open an output file
  //  fOutputFile = new TFile("FlowEvents.root","RECREATE");

}

//________________________________________________________________________
AliAnalysisTaskFlowEvent::AliAnalysisTaskFlowEvent() : 
  //  fOutputFile(NULL),
  fESD(NULL),
  fAOD(NULL),
  fEventMaker(NULL),
  fAnalysisType("ESD"),
  fCFManager1(NULL),
  fCFManager2(NULL),
  fQAInt(NULL),
  fQADiff(NULL),
  fMinMult(0),
  fMaxMult(10000000),
  fQA(kFALSE),
  fMCReactionPlaneAngle(0.),
  fCount(0),
  fNoOfLoops(1),
  fEllipticFlowValue(0.),
  fSigmaEllipticFlowValue(0.),
  fMultiplicityOfEvent(1000000000),
  fSigmaMultiplicityOfEvent(0),
  fMyTRandom3(NULL)
{
  // Constructor
  cout<<"AliAnalysisTaskFlowEvent::AliAnalysisTaskFlowEvent()"<<endl;
}

//________________________________________________________________________
AliAnalysisTaskFlowEvent::AliAnalysisTaskFlowEvent(const char *name, Bool_t on, UInt_t iseed) : 
  AliAnalysisTask(name, ""), 
//  fOutputFile(NULL),
  fESD(NULL),
  fAOD(NULL),
  fEventMaker(NULL),
  fAnalysisType("ESD"),
  fCFManager1(NULL),
  fCFManager2(NULL),
  fQAInt(NULL),
  fQADiff(NULL),
  fMinMult(0),
  fMaxMult(10000000),
  fQA(on),
  fMCReactionPlaneAngle(0.),
  fCount(0),
  fNoOfLoops(1),
  fEllipticFlowValue(0.),
  fSigmaEllipticFlowValue(0.),
  fMultiplicityOfEvent(1000000000),
  fSigmaMultiplicityOfEvent(0),
  fMyTRandom3(NULL)
{
  // Constructor
  cout<<"AliAnalysisTaskFlowEvent::AliAnalysisTaskFlowEvent(const char *name, Bool_t on, UInt_t iseed)"<<endl;

  fMyTRandom3 = new TRandom3(iseed);   
  gRandom->SetSeed(fMyTRandom3->Integer(65539));
  
  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Define here the flow event output
  DefineOutput(0, AliFlowEventSimple::Class());  
  if(on) {
    DefineOutput(1, TList::Class());
    DefineOutput(2, TList::Class()); }  
  // and for testing open an output file
  //  fOutputFile = new TFile("FlowEvents.root","RECREATE");

}

//________________________________________________________________________
AliAnalysisTaskFlowEvent::~AliAnalysisTaskFlowEvent()
{
  //
  // Destructor
  //
  if (fMyTRandom3) delete fMyTRandom3;
  // objects in the output list are deleted 
  // by the TSelector dtor (I hope)

}

//________________________________________________________________________
void AliAnalysisTaskFlowEvent::ConnectInputData(Option_t *) 
{
  // Connect ESD or AOD here
  // Called once
  cout<<"AliAnalysisTaskFlowEvent::ConnectInputData(Option_t *)"<<endl;

  TTree* tree = dynamic_cast<TTree*> (GetInputData(0));
  if (!tree) {
    Printf("ERROR: Could not read chain from input slot 0");
  } 
  else {
    if (fAnalysisType == "ESD" || fAnalysisType == "ESDMC0" || fAnalysisType == "ESDMC1" || fAnalysisType == "MC") {
      AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
      if (!esdH) {
	Printf("ERROR: Could not get ESDInputHandler");
      } 
      else {fESD = esdH->GetEvent();}
    }
    else if (fAnalysisType == "AOD") {
      AliAODInputHandler *aodH = dynamic_cast<AliAODInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
      if (!aodH) {
	Printf("ERROR: Could not get AODInputHandler");
      }
      else { fAOD = aodH->GetEvent();}
    }
    else {
      Printf("!!!!!Wrong analysis type: Only ESD, ESDMC0, ESDMC1, AOD and MC types are allowed!");
      exit(1);
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskFlowEvent::CreateOutputObjects() 
{
  // Called at every worker node to initialize
  cout<<"AliAnalysisTaskFlowEvent::CreateOutputObjects()"<<endl;

  if (!(fAnalysisType == "AOD" || fAnalysisType == "ESD" || fAnalysisType == "ESDMC0"  || fAnalysisType == "ESDMC1" || fAnalysisType == "MC")) {
    cout<<"WRONG ANALYSIS TYPE! only ESD, ESDMC0, ESDMC1, AOD and MC are allowed."<<endl;
    exit(1);
  }

  // Flow Event maker
  fEventMaker = new AliFlowEventSimpleMaker();
}

//________________________________________________________________________
void AliAnalysisTaskFlowEvent::Exec(Option_t *) 
{
  // Main loop
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
	      //cout<<"The reactionPlane from Hijing is: "<< fRP <<endl;
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

      //EPOS
      else if (!strcmp(mcEvent->GenEventHeader()->GetName(),"EPOS")) {
	AliGenEposEventHeader* headerE = (AliGenEposEventHeader*)(mcEvent->GenEventHeader());
	if (headerE) {
	  fRP = headerE->ReactionPlaneAngle();
	  //cout<<"The reactionPlane from EPOS is: "<< fR <<endl;
	}
	//else { cout<<"headerE is NULL"<<endl; }
      }

    }
    //else {cout<<"No MC event!"<<endl; }
    
  }
  //else {cout<<"No eventHandler!"<<endl; }


  fEventMaker->SetMCReactionPlaneAngle(fRP);
  //setting event cuts
  fEventMaker->SetMinMult(fMinMult);
  fEventMaker->SetMaxMult(fMaxMult);

  if (fEllipticFlowValue != 0.) {  
    // set the value of the monte carlo event plane for the flow event
    cout << "settings for afterburner in TaskFlowEvent.cxx:" << endl;
    cout << "fCount = " << fCount << endl;
    cout << "fNoOfLoops = " << fNoOfLoops << endl;
    cout << "fEllipticFlowValue = " << fEllipticFlowValue << endl;
    cout << "fSigmaEllipticFlowValue = " << fSigmaEllipticFlowValue << endl;
    cout << "fMultiplicityOfEvent = " << fMultiplicityOfEvent << endl;
    cout << "fSigmaMultiplicityOfEvent = " << fSigmaMultiplicityOfEvent << endl;

    Double_t xRPangle=0.;
    Double_t xNewFlowValue = 0.;
    Int_t nNewMultOfEvent = 100000000;

    if (fMyTRandom3) {  
      xRPangle = TMath::TwoPi()*(fMyTRandom3->Rndm());
      xNewFlowValue = fMyTRandom3->Gaus(fEllipticFlowValue,fSigmaEllipticFlowValue);
      nNewMultOfEvent =  TMath::Nint(fMyTRandom3->Gaus(fMultiplicityOfEvent,fSigmaMultiplicityOfEvent));
    }
    else {
      cout << "no random generator pointer initialized " << endl;
    }
    cout << "xRPangle = " << xRPangle << endl;
    cout << "xNewFlowValue = " << xNewFlowValue << endl;
    cout << "nNewMultOfEvent = " << nNewMultOfEvent << endl;
    cout << "settings for after burner" << endl;  

    fEventMaker->SetMCReactionPlaneAngle(xRPangle);
    fEventMaker->SetNoOfLoops(fNoOfLoops);
    fEventMaker->SetEllipticFlowValue(xNewFlowValue);
    fEventMaker->SetMultiplicityOfEvent(nNewMultOfEvent);  
    //end settings afterburner
  }  

  
  // Fill the FlowEventSimple for MC input          
  if (fAnalysisType == "MC") {
    if (!fCFManager1) {cout << "ERROR: No pointer to correction framework cuts! " << endl; return; }
    if (!fCFManager2) {cout << "ERROR: No pointer to correction framework cuts! " << endl; return; }

    // Process MC truth, therefore we receive the AliAnalysisManager and ask it for the AliMCEventHandler
    // This handler can return the current MC event
    if (!mcEvent) { Printf("ERROR: Could not retrieve MC event"); return;}

    fCFManager1->SetMCEventInfo(mcEvent);
    fCFManager2->SetMCEventInfo(mcEvent);

    // analysis 
    Printf("Number of MC particles: %d", mcEvent->GetNumberOfTracks());
    fEvent = fEventMaker->FillTracks(mcEvent,fCFManager1,fCFManager2);
    // here we have the fEvent and want to make it available as an output stream
    // so no delete fEvent;
  }
  // Fill the FlowEventSimple for ESD input  
  else if (fAnalysisType == "ESD") {
    if (!fCFManager1) {cout << "ERROR: No pointer to correction framework cuts! " << endl; return; }
    if (!fCFManager2) {cout << "ERROR: No pointer to correction framework cuts! " << endl; return; }

    if (!fESD) { Printf("ERROR: fESD not available"); return;}
    Printf("There are %d tracks in this event", fESD->GetNumberOfTracks());
    
    // analysis 
    fEvent = fEventMaker->FillTracks(fESD,fCFManager1,fCFManager2);
  }
  // Fill the FlowEventSimple for ESD input combined with MC info  
  else if (fAnalysisType == "ESDMC0" || fAnalysisType == "ESDMC1" ) {
    if (!fCFManager1) {cout << "ERROR: No pointer to correction framework cuts! " << endl; return; }
    if (!fCFManager2) {cout << "ERROR: No pointer to correction framework cuts! " << endl; return; }
    if (!fESD) { Printf("ERROR: fESD not available"); return;}
    Printf("There are %d tracks in this event", fESD->GetNumberOfTracks());
    
    if (!mcEvent) {Printf("ERROR: Could not retrieve MC event"); return;}

    fCFManager1->SetMCEventInfo(mcEvent);
    fCFManager2->SetMCEventInfo(mcEvent);


    if (fAnalysisType == "ESDMC0") { 
      fEvent = fEventMaker->FillTracks(fESD, mcEvent, fCFManager1, fCFManager2, 0); //0 = kine from ESD, 1 = kine from MC
    } else if (fAnalysisType == "ESDMC1") {
      fEvent = fEventMaker->FillTracks(fESD, mcEvent, fCFManager1, fCFManager2, 1); //0 = kine from ESD, 1 = kine from MC
    }
  }
  // Fill the FlowEventSimple for AOD input  
  else if (fAnalysisType == "AOD") {
    if (!fAOD) {Printf("ERROR: fAOD not available"); return;}
    Printf("There are %d tracks in this event", fAOD->GetNumberOfTracks());

    // analysis 
    //For the moment don't use CF //AliFlowEventSimple* fEvent = fEventMaker->FillTracks(fAOD,fCFManager1,fCFManager2);
    fEvent = fEventMaker->FillTracks(fAOD);
  }

  //fListHistos->Print();
  //  fOutputFile->WriteObject(fEvent,"myFlowEventSimple");	
  PostData(0,fEvent);
  if (fQA) {
    PostData(1,fQAInt);
    PostData(2,fQADiff); }
} 

//________________________________________________________________________
void AliAnalysisTaskFlowEvent::Terminate(Option_t *) 
{
  // Called once at the end of the query -- do not call in case of CAF

}


