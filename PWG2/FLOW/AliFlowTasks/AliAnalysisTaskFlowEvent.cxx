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
#include "TH2F.h"
#include "TRandom3.h"
#include "TTimeStamp.h"

// ALICE Analysis Framework
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliESDtrackCuts.h"

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

// Interface to Load short life particles
#include "TObjArray.h"
#include "AliFlowCandidateTrack.h"

// Interface to make the Flow Event Simple used in the flow analysis methods
#include "AliFlowEvent.h"
#include "AliFlowTrackCuts.h"
#include "AliFlowEventCuts.h"
#include "AliFlowCommonConstants.h"
#include "AliAnalysisTaskFlowEvent.h"

#include "AliLog.h"

ClassImp(AliAnalysisTaskFlowEvent)

//________________________________________________________________________
AliAnalysisTaskFlowEvent::AliAnalysisTaskFlowEvent() :
  AliAnalysisTaskSE(),
  //  fOutputFile(NULL),
  fAnalysisType("AUTOMATIC"),
  fRPType("Global"),
  fCFManager1(NULL),
  fCFManager2(NULL),
  fCutsEvent(NULL),
  fCutsRP(NULL),
  fCutsPOI(NULL),
  fQAInt(NULL),
  fQADiff(NULL),
  fMinMult(0),
  fMaxMult(10000000),
  fMinA(-1.0),
  fMaxA(-0.01),
  fMinB(0.01),
  fMaxB(1.0),
  fQA(kFALSE),
  fLoadCandidates(kFALSE),
  fNbinsMult(10000),
  fNbinsPt(100),   
  fNbinsPhi(100),
  fNbinsEta(200),
  fNbinsQ(500),
  fMultMin(0.),            
  fMultMax(10000.),
  fPtMin(0.),	     
  fPtMax(10.),
  fPhiMin(0.),	     
  fPhiMax(TMath::TwoPi()),
  fEtaMin(-5.),	     
  fEtaMax(5.),	     
  fQMin(0.),	     
  fQMax(3.),
  fExcludedEtaMin(0.), 
  fExcludedEtaMax(0.), 
  fExcludedPhiMin(0.),
  fExcludedPhiMax(0.),
  fAfterburnerOn(kFALSE),
  fNonFlowNumberOfTrackClones(0),
  fV1(0.),
  fV2(0.),
  fV3(0.),
  fV4(0.),
  fMyTRandom3(NULL)
{
  // Constructor
  cout<<"AliAnalysisTaskFlowEvent::AliAnalysisTaskFlowEvent()"<<endl;
}

//________________________________________________________________________
AliAnalysisTaskFlowEvent::AliAnalysisTaskFlowEvent(const char *name, TString RPtype, Bool_t on, UInt_t iseed, Bool_t bCandidates) :
  AliAnalysisTaskSE(name),
  //  fOutputFile(NULL),
  fAnalysisType("AUTOMATIC"),
  fRPType(RPtype),
  fCFManager1(NULL),
  fCFManager2(NULL),
  fCutsEvent(NULL),
  fCutsRP(NULL),
  fCutsPOI(NULL),
  fQAInt(NULL),
  fQADiff(NULL),
  fMinMult(0),
  fMaxMult(10000000),
  fMinA(-1.0),
  fMaxA(-0.01),
  fMinB(0.01),
  fMaxB(1.0),
  fQA(on),
  fLoadCandidates(bCandidates),
  fNbinsMult(10000),
  fNbinsPt(100),   
  fNbinsPhi(100),
  fNbinsEta(200),
  fNbinsQ(500),
  fMultMin(0.),            
  fMultMax(10000.),
  fPtMin(0.),	     
  fPtMax(10.),
  fPhiMin(0.),	     
  fPhiMax(TMath::TwoPi()),
  fEtaMin(-5.),	     
  fEtaMax(5.),	     
  fQMin(0.),	     
  fQMax(3.),
  fExcludedEtaMin(0.), 
  fExcludedEtaMax(0.), 
  fExcludedPhiMin(0.),
  fExcludedPhiMax(0.),
  fAfterburnerOn(kFALSE),
  fNonFlowNumberOfTrackClones(0),
  fV1(0.),
  fV2(0.),
  fV3(0.),
  fV4(0.),
  fMyTRandom3(NULL)
{
  // Constructor
  cout<<"AliAnalysisTaskFlowEvent::AliAnalysisTaskFlowEvent(const char *name, Bool_t on, UInt_t iseed)"<<endl;
  fMyTRandom3 = new TRandom3(iseed);
  gRandom->SetSeed(fMyTRandom3->Integer(65539));

  int availableINslot=1;
  //FMD input slot
  if (strcmp(RPtype,"FMD")==0) {
    DefineInput(availableINslot++, TList::Class());
  }
  //Candidates input slot
  if( fLoadCandidates )
    DefineInput(availableINslot, TObjArray::Class());

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
  delete fMyTRandom3;
  // objects in the output list are deleted
  // by the TSelector dtor (I hope)

}

//________________________________________________________________________
void AliAnalysisTaskFlowEvent::NotifyRun()
{
  //at the beginning of each new run
	AliESDEvent* fESD = dynamic_cast<AliESDEvent*> (InputEvent());
  if (!fESD) return;

  Int_t run = fESD->GetRunNumber();  
  AliInfo(Form("Stariting run #%i",run));
}

//________________________________________________________________________
void AliAnalysisTaskFlowEvent::UserCreateOutputObjects()
{
  // Called at every worker node to initialize
  cout<<"AliAnalysisTaskFlowEvent::CreateOutputObjects()"<<endl;

  if (!(fAnalysisType == "AOD" || fAnalysisType == "ESD" || fAnalysisType == "ESDMCkineESD"  || fAnalysisType == "ESDMCkineMC" || fAnalysisType == "MC" || fAnalysisType == "AUTOMATIC"))
  {
    AliError("WRONG ANALYSIS TYPE! only ESD, ESDMCkineESD, ESDMCkineMC, AOD, MC and AUTOMATIC are allowed.");
    exit(1);
  }

  //set the common constants
  AliFlowCommonConstants* cc = AliFlowCommonConstants::GetMaster();
  cc->SetNbinsMult(fNbinsMult);
  cc->SetNbinsPt(fNbinsPt);
  cc->SetNbinsPhi(fNbinsPhi); 
  cc->SetNbinsEta(fNbinsEta);
  cc->SetNbinsQ(fNbinsQ);
  cc->SetMultMin(fMultMin);
  cc->SetMultMax(fMultMax);
  cc->SetPtMin(fPtMin);
  cc->SetPtMax(fPtMax);
  cc->SetPhiMin(fPhiMin);
  cc->SetPhiMax(fPhiMax);
  cc->SetEtaMin(fEtaMin);
  cc->SetEtaMax(fEtaMax);
  cc->SetQMin(fQMin);
  cc->SetQMax(fQMax);

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
  AliMultiplicity* myTracklets = NULL;
  AliESDPmdTrack* pmdtracks = NULL;//pmd      

  int availableINslot=1;

  if (!(fCutsRP&&fCutsPOI&&fCutsEvent))
  {
    AliError("cuts not set");
    return;
  }

  //use the new and temporarily incomplete way of doing things
  if (fAnalysisType == "AUTOMATIC")
  {
    //check event cuts
    if (!fCutsEvent->IsSelected(InputEvent())) return;

    //first attach all possible information to the cuts
    fCutsRP->SetEvent( InputEvent(), MCEvent() );  //attach event
    fCutsPOI->SetEvent( InputEvent(), MCEvent() );

    //then make the event
    flowEvent = new AliFlowEvent( fCutsRP, fCutsPOI );
    if (myESD)
      flowEvent->SetReferenceMultiplicity(fCutsEvent->GetReferenceMultiplicity(InputEvent()));
    if (mcEvent && mcEvent->GenEventHeader()) flowEvent->SetMCReactionPlaneAngle(mcEvent);
  }

  // Make the FlowEvent for MC input
  else if (fAnalysisType == "MC")
  {
    // Process MC truth, therefore we receive the AliAnalysisManager and ask it for the AliMCEventHandler
    // This handler can return the current MC event
    if (!(fCFManager1&&fCFManager2))
    {
      AliError("ERROR: No pointer to correction framework cuts! ");
      return;
    }
    if (!mcEvent)
    {
      AliError("ERROR: Could not retrieve MC event");
      return;
    }

    fCFManager1->SetMCEventInfo(mcEvent);
    fCFManager2->SetMCEventInfo(mcEvent);
    
    AliInfo(Form("Number of MC particles: %d", mcEvent->GetNumberOfTracks()));

    //check multiplicity 
    if (!fCFManager1->CheckEventCuts(AliCFManager::kEvtGenCuts,mcEvent))
    {
      AliWarning("Event does not pass multiplicity cuts"); return;
    }
    //make event
    flowEvent = new AliFlowEvent(mcEvent,fCFManager1,fCFManager2);
  }

  // Make the FlowEvent for ESD input
  else if (fAnalysisType == "ESD")
  {
    if (!(fCFManager1&&fCFManager2))
    {
      AliError("ERROR: No pointer to correction framework cuts!");
      return;
    }
    if (!myESD)
    {
      AliError("ERROR: ESD not available");
      return;
    }

    //check the offline trigger (check if the event has the correct trigger)
    AliInfo(Form("ESD has %d tracks", fInputEvent->GetNumberOfTracks()));

    //check multiplicity
    if (!fCFManager1->CheckEventCuts(AliCFManager::kEvtRecCuts,myESD))
    {
      AliWarning("Event does not pass multiplicity cuts"); return;
    }

    //make event
    if (fRPType == "Global") {
      flowEvent = new AliFlowEvent(myESD,fCFManager1,fCFManager2);
    }
    else if (fRPType == "TPCOnly") {
      flowEvent = new AliFlowEvent(myESD,fCFManager2,kFALSE);
    }
    else if (fRPType == "TPCHybrid") {
      flowEvent = new AliFlowEvent(myESD,fCFManager2,kTRUE);
    }
    else if (fRPType == "Tracklet"){
      flowEvent = new AliFlowEvent(myESD,myTracklets,fCFManager2);
    }
    else if (fRPType == "FMD"){
      TList* FMDdata = dynamic_cast<TList*>(GetInputData(availableINslot++));
      if(!FMDdata) {
        cout<<" No FMDdata "<<endl;
        return;
      }
      TH2F* histFMD = dynamic_cast<TH2F*>(FMDdata->FindObject("dNdetadphiHistogramTrVtx"));
      if (!histFMD) {
        cout<< "No histFMD"<<endl;
        return;
      }
      flowEvent = new AliFlowEvent(myESD,histFMD,fCFManager2);
    }
    else if (fRPType == "PMD"){
      flowEvent = new AliFlowEvent(myESD,pmdtracks,fCFManager2);
    }
    else return;
    
    // if monte carlo event get reaction plane from monte carlo (depends on generator)
    if (mcEvent && mcEvent->GenEventHeader()) flowEvent->SetMCReactionPlaneAngle(mcEvent);
    //set reference multiplicity, TODO: maybe move it to the constructor?
    flowEvent->SetReferenceMultiplicity(AliESDtrackCuts::GetReferenceMultiplicity(myESD,kTRUE));
  }

  // Make the FlowEvent for ESD input combined with MC info
  else if (fAnalysisType == "ESDMCkineESD" || fAnalysisType == "ESDMCkineMC" )
  {
    if (!(fCFManager1&&fCFManager2))
    {
      AliError("ERROR: No pointer to correction framework cuts! ");
      return;
    }
    if (!myESD)
    {
      AliError("ERROR: ESD not available");
      return;
    }
    AliInfo(Form("There are %d tracks in this event", fInputEvent->GetNumberOfTracks()));

    if (!mcEvent)
    {
      AliError("ERROR: Could not retrieve MC event");
      return;
    }

    fCFManager1->SetMCEventInfo(mcEvent);
    fCFManager2->SetMCEventInfo(mcEvent);

    //check multiplicity
    if (!fCFManager1->CheckEventCuts(AliCFManager::kEvtRecCuts,myESD))
    {
      AliWarning("Event does not pass multiplicity cuts"); return;
    }

    //make event
    if (fAnalysisType == "ESDMCkineESD")
    {
      flowEvent = new AliFlowEvent(myESD, mcEvent, AliFlowEvent::kESDkine, fCFManager1, fCFManager2 );
    }
    else if (fAnalysisType == "ESDMCkineMC")
    {
      flowEvent = new AliFlowEvent(myESD, mcEvent, AliFlowEvent::kMCkine, fCFManager1, fCFManager2 );
    }
    if (!flowEvent) return;
    // if monte carlo event get reaction plane from monte carlo (depends on generator)
    if (mcEvent && mcEvent->GenEventHeader()) flowEvent->SetMCReactionPlaneAngle(mcEvent);
    //set reference multiplicity, TODO: maybe move it to the constructor?
    flowEvent->SetReferenceMultiplicity(AliESDtrackCuts::GetReferenceMultiplicity(myESD,kTRUE));
  }
  // Make the FlowEventSimple for AOD input
  else if (fAnalysisType == "AOD")
  {
    if (!myAOD)
    {
      AliError("ERROR: AOD not available");
      return;
    }
    AliInfo(Form("AOD has %d tracks", myAOD->GetNumberOfTracks()));
    flowEvent = new AliFlowEvent(myAOD);
  }

  //inject candidates
  if(fLoadCandidates) {
    TObjArray* Candidates = dynamic_cast<TObjArray*>(GetInputData(availableINslot++));
    //if(Candidates->GetEntriesFast()) 
      //printf("I received %d candidates\n",Candidates->GetEntriesFast());
    if (Candidates)
    {
      for(int iCand=0; iCand!=Candidates->GetEntriesFast(); ++iCand ) {
        AliFlowCandidateTrack *cand = dynamic_cast<AliFlowCandidateTrack*>(Candidates->At(iCand));
        if (!cand) continue;
        cand->SetForPOISelection(kTRUE);
        cand->SetForRPSelection(kFALSE);
        //printf(" Ⱶ Checking at candidate %d with %d daughters\n",iCand,cand->GetNDaughters());
        for(int iDau=0; iDau!=cand->GetNDaughters(); ++iDau) {
          //printf("    Ⱶ Daughter %d with fID %d", iDau, cand->GetIDDaughter(iDau) );
          for(int iRPs=0; iRPs!=flowEvent->NumberOfTracks(); ++iRPs ) {
            AliFlowTrack *iRP = dynamic_cast<AliFlowTrack*>(flowEvent->GetTrack( iRPs ));
            if (!iRP) continue;
            if( !iRP->InRPSelection() )
              continue;
            if( cand->GetIDDaughter(iDau) == iRP->GetID() ) {
              //printf(" was in RP set");
              cand->SetDaughter( iDau, iRP );
              //temporarily untagging all daugters
              iRP->SetForRPSelection(kFALSE);
            }
          }
          //printf("\n");
        }
        flowEvent->AddTrack(cand);
      }
    }
  }

  if (!flowEvent) return; //shuts up coverity

  //check final event cuts
  Int_t mult = flowEvent->NumberOfTracks();
  AliInfo(Form("FlowEvent has %i tracks",mult));
  if (mult<fMinMult || mult>fMaxMult)
  {
    AliWarning("FlowEvent cut on multiplicity"); return;
  }

  //define dead zone
  flowEvent->DefineDeadZone(fExcludedEtaMin, fExcludedEtaMax, fExcludedPhiMin, fExcludedPhiMax );


  //////////////////////////////////////////////////////////////////////////////
  ///////////////////////////AFTERBURNER
  if (fAfterburnerOn)
  {
    //if reaction plane not set from elsewhere randomize it before adding flow
    if (!flowEvent->IsSetMCReactionPlaneAngle())
      flowEvent->SetMCReactionPlaneAngle(gRandom->Uniform(0.0,TMath::TwoPi()));

    flowEvent->AddFlow(fV1,fV2,fV3,fV4);     //add flow
    flowEvent->CloneTracks(fNonFlowNumberOfTrackClones); //add nonflow by cloning tracks
  }
  //////////////////////////////////////////////////////////////////////////////

  //tag subEvents
  flowEvent->TagSubeventsInEta(fMinA,fMaxA,fMinB,fMaxB);

  //fListHistos->Print();
  //fOutputFile->WriteObject(flowEvent,"myFlowEventSimple");
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

