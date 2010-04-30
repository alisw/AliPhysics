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
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"

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
#include "AliFlowEvent.h"

#include "AliAnalysisTaskFlowEvent.h"

ClassImp(AliAnalysisTaskFlowEvent)

//________________________________________________________________________
AliAnalysisTaskFlowEvent::AliAnalysisTaskFlowEvent() :
  AliAnalysisTaskSE(),
  //  fOutputFile(NULL),
  fAnalysisType("ESD"),
  fCFManager1(NULL),
  fCFManager2(NULL),
  fQAInt(NULL),
  fQADiff(NULL),
  fMinMult(0),
  fMaxMult(10000000),
  fMinA(-1.0),
  fMaxA(-0.01),
  fMinB(0.01),
  fMaxB(1.0),
  fQA(kFALSE),
  fMCReactionPlaneAngle(0.),
  fCount(0),
  fNoOfLoops(1),
  fEllipticFlowValue(0.),
  fSigmaEllipticFlowValue(0.),
  fMultiplicityOfEvent(1000000000),
  fSigmaMultiplicityOfEvent(0),
  fMyTRandom3(NULL),
  fbAfterburnerOn(kFALSE)
{
  // Constructor
  cout<<"AliAnalysisTaskFlowEvent::AliAnalysisTaskFlowEvent()"<<endl;
}

//________________________________________________________________________
AliAnalysisTaskFlowEvent::AliAnalysisTaskFlowEvent(const char *name, Bool_t on, UInt_t iseed) :
  AliAnalysisTaskSE(name),
  //  fOutputFile(NULL),
  fAnalysisType("ESD"),
  fCFManager1(NULL),
  fCFManager2(NULL),
  fQAInt(NULL),
  fQADiff(NULL),
  fMinMult(0),
  fMaxMult(10000000),
  fMinA(-1.0),
  fMaxA(-0.01),
  fMinB(0.01),
  fMaxB(1.0),
  fQA(on),
  fMCReactionPlaneAngle(0.),
  fCount(0),
  fNoOfLoops(1),
  fEllipticFlowValue(0.),
  fSigmaEllipticFlowValue(0.),
  fMultiplicityOfEvent(1000000000),
  fSigmaMultiplicityOfEvent(0),
  fMyTRandom3(NULL),
  fbAfterburnerOn(kFALSE)
{
  // Constructor
  cout<<"AliAnalysisTaskFlowEvent::AliAnalysisTaskFlowEvent(const char *name, Bool_t on, UInt_t iseed)"<<endl;
  fMyTRandom3 = new TRandom3(iseed);
  gRandom->SetSeed(fMyTRandom3->Integer(65539));


  // Define output slots here
  // Define here the flow event output
  DefineOutput(1, AliFlowEventSimple::Class());
  if(on)
  {
    DefineOutput(2, TList::Class());
    DefineOutput(3, TList::Class());
  }
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
void AliAnalysisTaskFlowEvent::UserCreateOutputObjects()
{
  // Called at every worker node to initialize
  cout<<"AliAnalysisTaskFlowEvent::CreateOutputObjects()"<<endl;

  if (!(fAnalysisType == "AOD" || fAnalysisType == "ESD" || fAnalysisType == "ESDMCkineESD"  || fAnalysisType == "ESDMCkineMC" || fAnalysisType == "MC"))
  {
    cout<<"WRONG ANALYSIS TYPE! only ESD, ESDMCkineESD, ESDMCkineMC, AOD and MC are allowed."<<endl;
    exit(1);
  }
}

//________________________________________________________________________
void AliAnalysisTaskFlowEvent::UserExec(Option_t *)
{
  // Main loop
  // Called for each event
  AliFlowEvent* flowEvent = NULL;
  AliMCEvent*  mcEvent = MCEvent();                              // from TaskSE
  AliESDEvent* myESD = dynamic_cast<AliESDEvent*>(InputEvent()); // from TaskSE
  AliAODEvent* myAOD = dynamic_cast<AliAODEvent*>(InputEvent()); // from TaskSE

  // Make the FlowEvent for MC input
  if (fAnalysisType == "MC")
  {
    // Process MC truth, therefore we receive the AliAnalysisManager and ask it for the AliMCEventHandler
    // This handler can return the current MC event
    if (!(fCFManager1&&fCFManager2))
    {
      cout << "ERROR: No pointer to correction framework cuts! " << endl;
      return;
    }
    if (!mcEvent)
    {
      Printf("ERROR: Could not retrieve MC event");
      return;
    }

    fCFManager1->SetMCEventInfo(mcEvent);
    fCFManager2->SetMCEventInfo(mcEvent);
    
    Printf("Number of MC particles: %d", mcEvent->GetNumberOfTracks());
    // analysis
    flowEvent = new AliFlowEvent(mcEvent,fCFManager1,fCFManager2);
  }
  // Make the FlowEvent for ESD input
  else if (fAnalysisType == "ESD")
  {
    if (!(fCFManager1&&fCFManager2))
    {
      cout << "ERROR: No pointer to correction framework cuts! " << endl;
      return;
    }
    if (!myESD)
    {
      Printf("ERROR: ESD not available");
      return;
    }
    //check the offline trigger (check if the event has the correct trigger)
    Printf("There are %d tracks in this event", fInputEvent->GetNumberOfTracks());
    // analysis
    flowEvent = new AliFlowEvent(myESD,fCFManager1,fCFManager2);
  }
  // Make the FlowEvent for ESD input combined with MC info
  else if (fAnalysisType == "ESDMCkineESD" || fAnalysisType == "ESDMCkineMC" )
  {
    if (!(fCFManager1&&fCFManager2))
    {
      cout << "ERROR: No pointer to correction framework cuts! " << endl;
      return;
    }
    if (!myESD)
    {
      Printf("ERROR: ESD not available");
      return;
    }
    Printf("There are %d tracks in this event", fInputEvent->GetNumberOfTracks());

    if (!mcEvent)
    {
      Printf("ERROR: Could not retrieve MC event");
      return;
    }

    fCFManager1->SetMCEventInfo(mcEvent);
    fCFManager2->SetMCEventInfo(mcEvent);

    if (fAnalysisType == "ESDMCkineESD")
    {
      flowEvent = new AliFlowEvent(myESD, mcEvent, AliFlowEvent::kESDkine, fCFManager1, fCFManager2 );
    }
    else if (fAnalysisType == "ESDMCkineMC")
    {
      flowEvent = new AliFlowEvent(myESD, mcEvent, AliFlowEvent::kMCkine, fCFManager1, fCFManager2 );
    }
  }
  // Make the FlowEventSimple for AOD input
  else if (fAnalysisType == "AOD")
  {
    if (!myAOD)
    {
      Printf("ERROR: AOD not available");
      return;
    }
    Printf("There are %d tracks in this event", myAOD->GetNumberOfTracks());
    flowEvent = new AliFlowEvent(myAOD);
  }

  //check event cuts
  Int_t mult = flowEvent->NumberOfTracks();
  if (mult<fMinMult && mult>fMaxMult) return;

  //tag subevents
  flowEvent->TagSubeventsInEta(fMinA,fMaxA,fMinB,fMaxB);

  ////TODO afterburner
  //if (fbAfterburnerOn && fMyTRandom3) {
  //  // set the new value of the values using a after burner
  //  cout << "settings for afterburner in TaskFlowEvent.cxx:" << endl;
  //  cout << "fCount = " << fCount << endl;
  //  cout << "fNoOfLoops = " << fNoOfLoops << endl;
  //  cout << "fEllipticFlowValue = " << fEllipticFlowValue << endl;
  //  cout << "fSigmaEllipticFlowValue = " << fSigmaEllipticFlowValue << endl;
  //  cout << "fMultiplicityOflowEvent = " << fMultiplicityOflowEvent << endl;
  //  cout << "fSigmaMultiplicityOflowEvent = " << fSigmaMultiplicityOflowEvent << endl;
  //  Double_t xRPangle;
  //  if (!mcEvent) { xRPangle = TMath::TwoPi()*(fMyTRandom3->Rndm()); }
  //  else { xRPangle = fRP; }
  //  Double_t xNewFlowValue = fMyTRandom3->Gaus(fEllipticFlowValue,fSigmaEllipticFlowValue);
  //  Int_t nNewMultOflowEvent =  TMath::Nint(fMyTRandom3->Gaus(fMultiplicityOflowEvent,fSigmaMultiplicityOflowEvent));
  //  cout << "xRPangle = " << xRPangle << endl;
  //  cout << "xNewFlowValue = " << xNewFlowValue << endl;
  //  cout << "nNewMultOflowEvent = " << nNewMultOflowEvent << endl;
  //  cout << "settings for after burner" << endl;
  //  flowEventMaker->SetMCReactionPlaneAngle(xRPangle);
  //  flowEventMaker->SetNoOfLoops(fNoOfLoops);
  //  flowEventMaker->SetEllipticFlowValue(xNewFlowValue);
  //  flowEventMaker->SetMultiplicityOflowEvent(nNewMultOflowEvent);
  //  //end settings afterburner
  //}

  // if monte carlo event get reaction plane from monte carlo (depends on generator)
  if (mcEvent && mcEvent->GenEventHeader()) flowEvent->SetMCReactionPlaneAngle(mcEvent);

  //fListHistos->Print();
  //  fOutputFile->WriteObject(flowEvent,"myFlowEventSimple");
  PostData(1,flowEvent);
  if (fQA)
  {
    PostData(2,fQAInt);
    PostData(3,fQADiff);
  }
}

//________________________________________________________________________
void AliAnalysisTaskFlowEvent::Terminate(Option_t *)
{
  // Called once at the end of the query -- do not call in case of CAF
}


