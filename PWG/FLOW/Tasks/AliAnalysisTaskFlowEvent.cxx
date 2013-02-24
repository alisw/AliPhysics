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
#include "TF1.h"
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

using std::cout;
using std::endl;
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
  fCutContainer(NULL),
  fQAList(NULL),
  fMinMult(0),
  fMaxMult(10000000),
  fMinA(-1.0),
  fMaxA(-0.01),
  fMinB(0.01),
  fMaxB(1.0),
  fQAon(kFALSE),
  fLoadCandidates(kFALSE),
  fNbinsMult(10000),
  fNbinsPt(100),   
  fNbinsPhi(100),
  fNbinsEta(200),
  fNbinsQ(500),
  fNbinsMass(1),
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
  fMassMin(-1.),	     
  fMassMax(0.),
  fHistWeightvsPhiMin(0.),
  fHistWeightvsPhiMax(3.),
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
  fV5(0.),
  fDifferentialV2(0),
  fFlowEvent(NULL),
  fShuffleTracks(kFALSE),
  fMyTRandom3(NULL)
{
  // Constructor
  AliDebug(2,"AliAnalysisTaskFlowEvent::AliAnalysisTaskFlowEvent()");
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
  fCutContainer(new TList()),
  fQAList(NULL),
  fMinMult(0),
  fMaxMult(10000000),
  fMinA(-1.0),
  fMaxA(-0.01),
  fMinB(0.01),
  fMaxB(1.0),
  fQAon(on),
  fLoadCandidates(bCandidates),
  fNbinsMult(10000),
  fNbinsPt(100),   
  fNbinsPhi(100),
  fNbinsEta(200),
  fNbinsQ(500),
  fNbinsMass(1),
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
  fMassMin(-1.),	     
  fMassMax(0.),
  fHistWeightvsPhiMin(0.),
  fHistWeightvsPhiMax(3.),
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
  fV5(0.),
  fDifferentialV2(0),
  fFlowEvent(NULL),
  fShuffleTracks(kFALSE),
  fMyTRandom3(NULL)
{
  // Constructor
  AliDebug(2,"AliAnalysisTaskFlowEvent::AliAnalysisTaskFlowEvent(const char *name, Bool_t on, UInt_t iseed)");
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
  DefineOutput(2, TList::Class());

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
  delete fFlowEvent;
  delete fCutsEvent;
  delete fQAList;
  if (fCutContainer) fCutContainer->Delete(); delete fCutContainer;
  // objects in the output list are deleted
  // by the TSelector dtor (I hope)

}

//________________________________________________________________________
void AliAnalysisTaskFlowEvent::NotifyRun()
{
  //at the beginning of each new run
  if (fCutsRP)   fCutsRP->SetRunsMuon(fInputHandler);  // XZhang 20120604
  if (fCutsPOI) fCutsPOI->SetRunsMuon(fInputHandler);  // XZhang 20120604
	AliESDEvent* fESD = dynamic_cast<AliESDEvent*> (InputEvent());
  if (!fESD) return;

  Int_t run = fESD->GetRunNumber();  
  AliInfo(Form("Stariting run #%i",run));
}

//________________________________________________________________________
void AliAnalysisTaskFlowEvent::UserCreateOutputObjects()
{
  // Called at every worker node to initialize
  AliDebug(2,"AliAnalysisTaskFlowEvent::CreateOutputObjects()");

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
  cc->SetNbinsMass(fNbinsMass);
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
  cc->SetMassMin(fMassMin);
  cc->SetMassMax(fMassMax);
  cc->SetHistWeightvsPhiMax(fHistWeightvsPhiMax);
  cc->SetHistWeightvsPhiMin(fHistWeightvsPhiMin);

  fFlowEvent = new AliFlowEvent(10000);

  if (fQAon)
  {
    fQAList=new TList();
    fQAList->SetOwner();
    fQAList->SetName(Form("%s QA",GetName()));
    if (fCutsEvent->GetQA()) fQAList->Add(fCutsEvent->GetQA()); //0
    if (fCutsRP->GetQA()) fQAList->Add(fCutsRP->GetQA());  //1
    if (fCutsPOI->GetQA())fQAList->Add(fCutsPOI->GetQA()); //2
    fQAList->Add(new TH1F("event plane angle","event plane angle;angle [rad];",100,0.,TMath::TwoPi())); //3
    PostData(2,fQAList);
  }
}

//________________________________________________________________________
void AliAnalysisTaskFlowEvent::UserExec(Option_t *)
{
  // Main loop
  // Called for each event
  //delete fFlowEvent;
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

  //DEFAULT - automatically takes care of everything
  if (fAnalysisType == "AUTOMATIC")
  {
    //check event cuts
    if (InputEvent() && !fCutsEvent->IsSelected(InputEvent())) return;

    //first attach all possible information to the cuts
    fCutsRP->SetEvent( InputEvent(), MCEvent() );  //attach event
    fCutsPOI->SetEvent( InputEvent(), MCEvent() );

    //then make the event
    fFlowEvent->Fill( fCutsRP, fCutsPOI );
    //fFlowEvent = new AliFlowEvent( fCutsRP, fCutsPOI );

    //    if (myESD)
      fFlowEvent->SetReferenceMultiplicity(fCutsEvent->GetReferenceMultiplicity(InputEvent()));
    if (mcEvent && mcEvent->GenEventHeader()) fFlowEvent->SetMCReactionPlaneAngle(mcEvent);
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
    fFlowEvent = new AliFlowEvent(mcEvent,fCFManager1,fCFManager2);
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
      fFlowEvent = new AliFlowEvent(myESD,fCFManager1,fCFManager2);
    }
    else if (fRPType == "TPCOnly") {
      fFlowEvent = new AliFlowEvent(myESD,fCFManager2,kFALSE);
    }
    else if (fRPType == "TPCHybrid") {
      fFlowEvent = new AliFlowEvent(myESD,fCFManager2,kTRUE);
    }
    else if (fRPType == "Tracklet"){
      fFlowEvent = new AliFlowEvent(myESD,myTracklets,fCFManager2);
    }
    else if (fRPType == "FMD"){
      TList* dataFMD = dynamic_cast<TList*>(GetInputData(availableINslot++));
      if(!dataFMD) {
        cout<<" No dataFMD "<<endl;
        return;
      }
      TH2F* histFMD = dynamic_cast<TH2F*>(dataFMD->FindObject("dNdetadphiHistogramTrVtx"));
      if (!histFMD) {
        cout<< "No histFMD"<<endl;
        return;
      }
      fFlowEvent = new AliFlowEvent(myESD,histFMD,fCFManager2);
    }
    else if (fRPType == "PMD"){
      fFlowEvent = new AliFlowEvent(myESD,pmdtracks,fCFManager2);
    }
    else return;
    
    // if monte carlo event get reaction plane from monte carlo (depends on generator)
    if (mcEvent && mcEvent->GenEventHeader()) fFlowEvent->SetMCReactionPlaneAngle(mcEvent);
    //set reference multiplicity, TODO: maybe move it to the constructor?
    fFlowEvent->SetReferenceMultiplicity(AliESDtrackCuts::GetReferenceMultiplicity(myESD,kTRUE));
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
      fFlowEvent = new AliFlowEvent(myESD, mcEvent, AliFlowEvent::kESDkine, fCFManager1, fCFManager2 );
    }
    else if (fAnalysisType == "ESDMCkineMC")
    {
      fFlowEvent = new AliFlowEvent(myESD, mcEvent, AliFlowEvent::kMCkine, fCFManager1, fCFManager2 );
    }
    if (!fFlowEvent) return;
    // if monte carlo event get reaction plane from monte carlo (depends on generator)
    if (mcEvent && mcEvent->GenEventHeader()) fFlowEvent->SetMCReactionPlaneAngle(mcEvent);
    //set reference multiplicity, TODO: maybe move it to the constructor?
    fFlowEvent->SetReferenceMultiplicity(AliESDtrackCuts::GetReferenceMultiplicity(myESD,kTRUE));
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
    fFlowEvent = new AliFlowEvent(myAOD);
  }

  //inject candidates
  if(fLoadCandidates) {
    TObjArray* candidates = dynamic_cast<TObjArray*>(GetInputData(availableINslot++));
    //if(candidates->GetEntriesFast()) 
    //  printf("I received %d candidates\n",candidates->GetEntriesFast());
    if (candidates)
    {
      for(int iCand=0; iCand!=candidates->GetEntriesFast(); ++iCand ) {
        AliFlowCandidateTrack *cand = dynamic_cast<AliFlowCandidateTrack*>(candidates->At(iCand));
        if (!cand) continue;
        //printf(" - Checking at candidate %d with %d daughters: mass %f\n",iCand,cand->GetNDaughters(),cand->Mass());
        for(int iDau=0; iDau!=cand->GetNDaughters(); ++iDau) {
          //printf("    - Daughter %d with fID %d", iDau, cand->GetIDDaughter(iDau) );
          for(int iRPs=0; iRPs!=fFlowEvent->NumberOfTracks(); ++iRPs ) {
            AliFlowTrack *iRP = dynamic_cast<AliFlowTrack*>(fFlowEvent->GetTrack( iRPs ));
            if (!iRP) continue;
            if( !iRP->InRPSelection() )
              continue;
            if( cand->GetIDDaughter(iDau) == iRP->GetID() ) {
              //printf(" was in RP set");
              //cand->SetDaughter( iDau, iRP );
              //temporarily untagging all daugters
              iRP->SetForRPSelection(kFALSE);
	      fFlowEvent->SetNumberOfRPs( fFlowEvent->GetNumberOfRPs() -1 );
            }
          }
          //printf("\n");
        }
	cand->SetForPOISelection(kTRUE);
	fFlowEvent->InsertTrack( ((AliFlowTrack*) cand) );
      }
    }
  }

  if (!fFlowEvent) return; //shuts up coverity

  //check final event cuts
  Int_t mult = fFlowEvent->NumberOfTracks();
  //  AliInfo(Form("FlowEvent has %i tracks",mult));
  if (mult<fMinMult || mult>fMaxMult)
  {
    AliWarning("FlowEvent cut on multiplicity"); return;
  }

  //define dead zone
  fFlowEvent->DefineDeadZone(fExcludedEtaMin, fExcludedEtaMax, fExcludedPhiMin, fExcludedPhiMax );


  //////////////////////////////////////////////////////////////////////////////
  ///////////////////////////AFTERBURNER
  if (fAfterburnerOn)
  {
    //if reaction plane not set from elsewhere randomize it before adding flow
    if (!fFlowEvent->IsSetMCReactionPlaneAngle())
      fFlowEvent->SetMCReactionPlaneAngle(gRandom->Uniform(0.0,TMath::TwoPi()));

    if(fDifferentialV2)
      fFlowEvent->AddV2(fDifferentialV2);
    else 
      fFlowEvent->AddFlow(fV1,fV2,fV3,fV4,fV5);     //add flow
    fFlowEvent->CloneTracks(fNonFlowNumberOfTrackClones); //add nonflow by cloning tracks
  }
  //////////////////////////////////////////////////////////////////////////////

  //tag subEvents
  fFlowEvent->TagSubeventsInEta(fMinA,fMaxA,fMinB,fMaxB);

  //QA
  if (fQAon)
  {
    TH1* h1 = static_cast<TH1*>(fQAList->FindObject("event plane angle"));
    h1->Fill(fFlowEvent->GetMCReactionPlaneAngle());
  }

  //do we want to serve shullfed tracks to everybody?
  fFlowEvent->SetShuffleTracks(fShuffleTracks);

  //fListHistos->Print();
  //fOutputFile->WriteObject(fFlowEvent,"myFlowEventSimple");
  PostData(1,fFlowEvent);
}

//________________________________________________________________________
void AliAnalysisTaskFlowEvent::Terminate(Option_t *)
{
  // Called once at the end of the query -- do not call in case of CAF
}

