//_________________________________________________________________________
//  Utility Class for transverse energy studies; charged hadrons
//  Task for analysis
//  - reconstruction and MC output
// implementation file
//
//Created by Christine Nattrass, Rebecca Scott, Irakli Martashvili
//University of Tennessee at Knoxville
//_________________________________________________________________________
//Necessary to read config macros
#include <TROOT.h>
#include <TSystem.h>
#include <TInterpreter.h>

#include "TChain.h"
#include "TList.h"
#include "TH2F.h"

#include "AliESDEvent.h"
#include "AliMCEvent.h"
#include "AliESDtrackCuts.h"

#include "AliAnalysisTaskHadEt.h"
#include "AliAnalysisHadEtReconstructed.h"
#include "AliAnalysisHadEtMonteCarlo.h"
#include "AliPWG0Helper.h"
#include "AliTriggerAnalysis.h"
#include "AliAnalysisManager.h"
#include "AliPIDResponse.h"
#include "AliTPCPIDResponse.h" 
#include "AliInputEventHandler.h"

#include <iostream>
#include "AliLog.h"

using namespace std;

ClassImp(AliAnalysisTaskHadEt)



//________________________________________________________________________
  AliAnalysisTaskHadEt::AliAnalysisTaskHadEt(const char *name, Bool_t isMc, TString recoConfigFile, TString mcConfigFile) :
        AliAnalysisTaskTransverseEnergy(name, isMc)
	,fPIDResponse(0)
	,fRecAnalysis(0)
	,fMCAnalysis(0)
	,fIsSim(isMc)
	,kIsOfflineV0AND(0)
	,kIsOfflineMB(0)
{
    // Constructor
  //input hander
  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  if (!man) {
    AliFatal("Analysis manager needed");
    return;
  }

  AliInputEventHandler *inputHandler=dynamic_cast<AliInputEventHandler*>(man->GetInputEventHandler());
  if (!inputHandler) {
    AliFatal("Input handler needed");
    return;
  }

  //pid response object
  fPIDResponse=inputHandler->GetPIDResponse();
  if (!fPIDResponse) AliError("PIDResponse object was not created");
  else{cout<<"PIDResponse was created!"<<endl;}


  fMCConfigFile = mcConfigFile;
  fRecoConfigFile = recoConfigFile;

  if(fMCAnalysis) delete fMCAnalysis;
  if(fRecAnalysis) delete fRecAnalysis;

  if (fRecoConfigFile.Length()) {
    cout<<"Rereading AliAnalysisHadEtReconstructed configuration file "<<fRecoConfigFile<<endl;
    gROOT->LoadMacro(fRecoConfigFile);
    fRecAnalysis = (AliAnalysisHadEtReconstructed *) gInterpreter->ProcessLine("ConfigHadEtReconstructed()");
  }

  if (fMCConfigFile.Length()) {
    cout<<"Rereading AliAnalysisHadEtMonteCarlo configuration file "<<fMCConfigFile<<endl;
    gROOT->LoadMacro(fMCConfigFile);
    fMCAnalysis = (AliAnalysisHadEtMonteCarlo *) gInterpreter->ProcessLine("ConfigHadEtMonteCarlo()");
    fMCAnalysis->SetHadEtReconstructed(fRecAnalysis);
  }

    // Define input and output slots here
    // Input slot #0 works with a TChain
    DefineInput(0, TChain::Class());
    // Output slot #1 writes into a TH1 container

    DefineOutput(1, TList::Class());
}
AliAnalysisTaskHadEt::~AliAnalysisTaskHadEt(){//Destructor
  delete fRecAnalysis;
  delete fMCAnalysis;
  delete fPIDResponse;
}


//________________________________________________________________________
void AliAnalysisTaskHadEt::UserCreateOutputObjects()
{
    // Create histograms

    // Called once


  //input hander
  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  if (!man) {
    AliFatal("Analysis manager needed");
    return;
  }

  AliInputEventHandler *inputHandler=dynamic_cast<AliInputEventHandler*>(man->GetInputEventHandler());
  if (!inputHandler) {
    AliFatal("Input handler needed");
    return;
  }

  //pid response object
  fPIDResponse=inputHandler->GetPIDResponse();
  if (!fPIDResponse) AliError("PIDResponse object was not created");
  else{cout<<"PIDResponse was created!"<<endl;}


  fOutputList = new TList;
  fOutputList->SetOwner();
  fMCAnalysis->SetHistoList(fOutputList);
  fRecAnalysis->SetHistoList(fOutputList);
  if(fIsSim) fMCAnalysis->CreateHistograms();
  fRecAnalysis->CreateHistograms();


  if(fRecAnalysis->DataSet() != fMCAnalysis->DataSet()){
    cout<<"Warning: Reconstruction data set and Monte Carlo data set are not the same!  Setting data set to "<<fRecAnalysis->DataSet()<<endl;
  }

  Bool_t selectPrimaries=kTRUE;
  if(fEsdtrackCutsITSTPC) delete fEsdtrackCutsITSTPC;
  if(fEsdtrackCutsITS) delete fEsdtrackCutsITS;
  if(fEsdtrackCutsTPC) delete fEsdtrackCutsTPC;
  //We do not use these because we are using the 2010 900 GeV data
//   if(fRecAnalysis->DataSet()==2009){
//     cout<<"Setting track cuts for the 2009 p+p collisions at 900 GeV"<<endl;
//     fEsdtrackCutsITSTPC = AliESDtrackCuts::GetStandardITSTPCTrackCuts2009(selectPrimaries);
//     fEsdtrackCutsITSTPC->SetName("fEsdTrackCuts");
//     fEsdtrackCutsTPC = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
//     fEsdtrackCutsTPC->SetName("fEsdTrackCutsTPCOnly");
//     //ITS stand alone cuts - similar to 2009 cuts but with only ITS hits required
//     fEsdtrackCutsITS =  AliESDtrackCuts::GetStandardITSPureSATrackCuts2009(kTRUE,kFALSE);//we do want primaries but we do not want to require PID info
//     fEsdtrackCutsITS->SetName("fEsdTrackCutsITS");
//   }
  if(fRecAnalysis->DataSet()==2010 || fRecAnalysis->DataSet()==20111||fRecAnalysis->DataSet()==2009 || fRecAnalysis->DataSet()==2012 || fRecAnalysis->DataSet()==2013){
    // AliAnalysisTaskSE::	SelectCollisionCandidates(AliVEvent::kINT7 ) ;
    if(fRecAnalysis->DataSet()==2010)cout<<"Setting track cuts for the 2010 p+p collisions at 7 TeV"<<endl;
    else{
      if(fRecAnalysis->DataSet()==2012)cout<<"Setting track cuts for the 2012 p+p collisions at 8 TeV"<<endl;
      else{
	if(fRecAnalysis->DataSet()==2013)cout<<"Setting track cuts for the 2013 p+Pb collisions at 5 TeV"<<endl;
	else{
	  if(fRecAnalysis->DataSet()==2009){cout<<"Setting track cuts for the 2010 p+p collisions at 900 GeV"<<endl;}
	  else{cout<<"Setting track cuts for the 2011 p+p collisions at 2.76 TeV"<<endl;}
	}
      }
    }
    //cout<<"Warning:  Have not set 2010 track cuts yet!!"<<endl;
    fEsdtrackCutsITSTPC = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(selectPrimaries);
    fEsdtrackCutsITSTPC->SetName("fEsdTrackCuts");
    fEsdtrackCutsTPC = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
    fEsdtrackCutsTPC->SetName("fEsdTrackCutsTPCOnly");
    //ITS stand alone cuts - similar to 2009 cuts but with only ITS hits required
    fEsdtrackCutsITS =  AliESDtrackCuts::GetStandardITSSATrackCuts2010(kTRUE,kFALSE);//we do want primaries but we do not want to require PID info
    fEsdtrackCutsITS->SetName("fEsdTrackCutsITS");
  }
  if(fRecAnalysis->DataSet()==20100){
    cout<<"Setting track cuts for the 2010 Pb+Pb collisions at 2.76 TeV"<<endl;
    //cout<<"Warning:  Have not set 2010 track cuts yet!!"<<endl;
    fEsdtrackCutsITSTPC = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(selectPrimaries);
    fEsdtrackCutsITSTPC->SetName("fEsdTrackCuts");
    fEsdtrackCutsTPC = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
    fEsdtrackCutsTPC->SetName("fEsdTrackCutsTPCOnly");
    //ITS stand alone cuts - similar to 2009 cuts but with only ITS hits required
    fEsdtrackCutsITS =  AliESDtrackCuts::GetStandardITSSATrackCutsPbPb2010(kTRUE,kFALSE);//we do want primaries but we do not want to require PID info
    // fEsdtrackCutsITS =  AliESDtrackCuts::GetStandardITSPureSATrackCuts2010(kTRUE,kFALSE);//we do want primaries but we do not want to require PID info
   fEsdtrackCutsITS->SetName("fEsdTrackCutsITS");
  }

  fOutputList->Add(fEsdtrackCutsITSTPC);
  fOutputList->Add(fEsdtrackCutsTPC);
  fOutputList->Add(fEsdtrackCutsITS);
  if(fEsdtrackCutsITSTPC && fEsdtrackCutsTPC){
    fRecAnalysis->SetITSTrackCuts( GetITSTrackCuts());
    fMCAnalysis->SetITSTrackCuts( GetITSTrackCuts());
    fRecAnalysis->SetTPCITSTrackCuts( GetTPCITSTrackCuts());
    fMCAnalysis->SetTPCITSTrackCuts( GetTPCITSTrackCuts());
    fRecAnalysis->SetTPCOnlyTrackCuts( GetTPCOnlyTrackCuts());
    fMCAnalysis->SetTPCOnlyTrackCuts( GetTPCOnlyTrackCuts());
    //add ITS stuff!
  }
  else{
    Printf("Error: no track cuts!");
  }



 PostData(1, fOutputList);
}

//________________________________________________________________________
void AliAnalysisTaskHadEt::UserExec(Option_t *)
{ // execute method
  fESDEvent = dynamic_cast<AliESDEvent*>(InputEvent());
if (!fESDEvent) {
  Printf("ERROR: Could not retrieve event");
  return;
 }
//cout<<"AliAnalysisTaskHadEt 165"<<endl;

//Int_t res = CheckPhysicsSelection(fESDEvent->GetRunNumber()); // Check if the physics selection is valid for this run

//AliCentrality *cent = GetCentralityObject();

//if(res == 0 && cent){
//if(cent){
 AliTriggerAnalysis *fTriggerAnalysis = new AliTriggerAnalysis();

  kIsOfflineV0AND = fTriggerAnalysis->IsOfflineTriggerFired(fESDEvent, AliTriggerAnalysis::kV0AND);  
  kIsOfflineMB = fTriggerAnalysis->IsOfflineTriggerFired(fESDEvent, AliTriggerAnalysis::kMB1);  
  fRecAnalysis->SetIsOfflineV0AND(kIsOfflineV0AND);
  fMCAnalysis->SetIsOfflineV0AND(kIsOfflineV0AND);
  fMCAnalysis->SetIsOfflineMB(kIsOfflineMB);

  Int_t eventtype = 	AliPWG0Helper::kInvalidProcess;
  if(fIsSim && fRecAnalysis->DataSet()!=20100) eventtype = (Int_t) AliPWG0Helper::GetEventProcessType(MCEvent()->Header());
  fRecAnalysis->AnalyseEvent(fESDEvent,eventtype);
  if(fIsSim){
    AliMCEvent* mcEvent = MCEvent();
    if(!mcEvent){  
      AliFatal("ERROR: MC Event does not exist");
      return;
    }
    if (fESDEvent){
      ((AliAnalysisHadEtMonteCarlo*)fMCAnalysis)->AnalyseEvent((AliVEvent*)mcEvent,(AliVEvent*)fESDEvent);

    if( fRecAnalysis->DataSet()==20100 || AliPWG0Helper::GetEventProcessType(mcEvent->Header()) == AliPWG0Helper::kND){//either non-diffractive or Pb+Pb
      if(fMCAnalysis->Full()){
	fMCAnalysis->FillSimTotEtMinusRecoTotEtFullAcceptanceTPC( fRecAnalysis->GetCorrectedTotEtFullAcceptanceTPC() );
	fMCAnalysis->FillSimTotEtMinusRecoTotEtFullAcceptanceITS( fRecAnalysis->GetCorrectedTotEtFullAcceptanceITS() );
	fMCAnalysis->FillSimTotEtMinusRecoTotEtFullAcceptanceTPCNoPID( fRecAnalysis->GetCorrectedTotEtFullAcceptanceTPCNoPID() );
	fMCAnalysis->FillSimTotEtMinusRecoTotEtFullAcceptanceITSNoPID( fRecAnalysis->GetCorrectedTotEtFullAcceptanceITSNoPID() );
	fMCAnalysis->FillSimHadEtMinusRecoHadEtFullAcceptanceTPC( fRecAnalysis->GetCorrectedHadEtFullAcceptanceTPC() );
	fMCAnalysis->FillSimHadEtMinusRecoHadEtFullAcceptanceITS( fRecAnalysis->GetCorrectedHadEtFullAcceptanceITS() );
	fMCAnalysis->FillSimHadEtMinusRecoHadEtFullAcceptanceTPCNoPID( fRecAnalysis->GetCorrectedHadEtFullAcceptanceTPCNoPID() );
	fMCAnalysis->FillSimHadEtMinusRecoHadEtFullAcceptanceITSNoPID( fRecAnalysis->GetCorrectedHadEtFullAcceptanceITSNoPID() );

	fMCAnalysis->FillSimTotEtVsRecoTotEtFullAcceptanceTPC( fRecAnalysis->GetCorrectedTotEtFullAcceptanceTPC() );
	fMCAnalysis->FillSimTotEtVsRecoTotEtFullAcceptanceITS( fRecAnalysis->GetCorrectedTotEtFullAcceptanceITS() );
	fMCAnalysis->FillSimTotEtVsRecoTotEtFullAcceptanceTPCNoPID( fRecAnalysis->GetCorrectedTotEtFullAcceptanceTPCNoPID() );
	fMCAnalysis->FillSimTotEtVsRecoTotEtFullAcceptanceITSNoPID( fRecAnalysis->GetCorrectedTotEtFullAcceptanceITSNoPID() );
	fMCAnalysis->FillSimHadEtVsRecoHadEtFullAcceptanceTPC( fRecAnalysis->GetCorrectedHadEtFullAcceptanceTPC() );
	fMCAnalysis->FillSimHadEtVsRecoHadEtFullAcceptanceITS( fRecAnalysis->GetCorrectedHadEtFullAcceptanceITS() );
	fMCAnalysis->FillSimHadEtVsRecoHadEtFullAcceptanceTPCNoPID( fRecAnalysis->GetCorrectedHadEtFullAcceptanceTPCNoPID() );
	fMCAnalysis->FillSimHadEtVsRecoHadEtFullAcceptanceITSNoPID( fRecAnalysis->GetCorrectedHadEtFullAcceptanceITSNoPID() );

	fMCAnalysis->FillSimPiKPEtVsRecoPiKPEtFullAcceptanceTPC( fRecAnalysis->GetCorrectedPiKPEtFullAcceptanceTPC() );
	fMCAnalysis->FillSimPiKPEtVsRecoPiKPEtFullAcceptanceITS( fRecAnalysis->GetCorrectedPiKPEtFullAcceptanceITS() );
	fMCAnalysis->FillSimPiKPEtVsRecoPiKPEtFullAcceptanceTPCNoPID( fRecAnalysis->GetCorrectedPiKPEtFullAcceptanceTPCNoPID() );
	fMCAnalysis->FillSimPiKPEtVsRecoPiKPEtFullAcceptanceITSNoPID( fRecAnalysis->GetCorrectedPiKPEtFullAcceptanceITSNoPID() );//Had
      }
      if(fMCAnalysis->EMCAL()){
	fMCAnalysis->FillSimTotEtMinusRecoTotEtEMCALAcceptanceTPC( fRecAnalysis->GetCorrectedTotEtEMCALAcceptanceTPC() );
	fMCAnalysis->FillSimTotEtMinusRecoTotEtEMCALAcceptanceITS( fRecAnalysis->GetCorrectedTotEtEMCALAcceptanceITS() );
	fMCAnalysis->FillSimTotEtMinusRecoTotEtEMCALAcceptanceTPCNoPID( fRecAnalysis->GetCorrectedTotEtEMCALAcceptanceTPCNoPID() );
	fMCAnalysis->FillSimTotEtMinusRecoTotEtEMCALAcceptanceITSNoPID( fRecAnalysis->GetCorrectedTotEtEMCALAcceptanceITSNoPID() );
	fMCAnalysis->FillSimHadEtMinusRecoHadEtEMCALAcceptanceTPC( fRecAnalysis->GetCorrectedHadEtEMCALAcceptanceTPC() );
	fMCAnalysis->FillSimHadEtMinusRecoHadEtEMCALAcceptanceITS( fRecAnalysis->GetCorrectedHadEtEMCALAcceptanceITS() );
	fMCAnalysis->FillSimHadEtMinusRecoHadEtEMCALAcceptanceTPCNoPID( fRecAnalysis->GetCorrectedHadEtEMCALAcceptanceTPCNoPID() );
	fMCAnalysis->FillSimHadEtMinusRecoHadEtEMCALAcceptanceITSNoPID( fRecAnalysis->GetCorrectedHadEtEMCALAcceptanceITSNoPID() );
      }
      if(fMCAnalysis->PHOS()){
	fMCAnalysis->FillSimTotEtMinusRecoTotEtPHOSAcceptanceTPC( fRecAnalysis->GetCorrectedTotEtPHOSAcceptanceTPC() );
	fMCAnalysis->FillSimTotEtMinusRecoTotEtPHOSAcceptanceITS( fRecAnalysis->GetCorrectedTotEtPHOSAcceptanceITS() );
	fMCAnalysis->FillSimTotEtMinusRecoTotEtPHOSAcceptanceTPCNoPID( fRecAnalysis->GetCorrectedTotEtPHOSAcceptanceTPCNoPID() );
	fMCAnalysis->FillSimTotEtMinusRecoTotEtPHOSAcceptanceITSNoPID( fRecAnalysis->GetCorrectedTotEtPHOSAcceptanceITSNoPID() );
	fMCAnalysis->FillSimHadEtMinusRecoHadEtPHOSAcceptanceTPC( fRecAnalysis->GetCorrectedHadEtPHOSAcceptanceTPC() );
	fMCAnalysis->FillSimHadEtMinusRecoHadEtPHOSAcceptanceITS( fRecAnalysis->GetCorrectedHadEtPHOSAcceptanceITS() );
	fMCAnalysis->FillSimHadEtMinusRecoHadEtPHOSAcceptanceTPCNoPID( fRecAnalysis->GetCorrectedHadEtPHOSAcceptanceTPCNoPID() );
	fMCAnalysis->FillSimHadEtMinusRecoHadEtPHOSAcceptanceITSNoPID( fRecAnalysis->GetCorrectedHadEtPHOSAcceptanceITSNoPID() );
      }
      if(fMCAnalysis->PiKP() && fMCAnalysis->Full()){
	fMCAnalysis->FillSimPiKPMinusRecoPiKPFullAcceptanceTPC(fRecAnalysis->GetCorrectedPiKPEtFullAcceptanceTPC());
	fMCAnalysis->FillSimPiKPMinusRecoPiKPFullAcceptanceITS(fRecAnalysis->GetCorrectedPiKPEtFullAcceptanceITS());
	fMCAnalysis->FillSimPiKPMinusRecoPiKPFullAcceptanceTPCNoPID(fRecAnalysis->GetCorrectedPiKPEtFullAcceptanceTPCNoPID());
	fMCAnalysis->FillSimPiKPMinusRecoPiKPFullAcceptanceITSNoPID(fRecAnalysis->GetCorrectedPiKPEtFullAcceptanceITSNoPID());
      }
    }
  }
  }
  delete fTriggerAnalysis;
  //}
//cout<<"End Event"<<endl<<endl;
// Post output data.
 PostData(1, fOutputList);
}

//________________________________________________________________________
void AliAnalysisTaskHadEt::Terminate(Option_t *)
{
    // Draw result to the screen
    // Called once at the end of the query

    fOutputList = dynamic_cast<TList*> (GetOutputData(1));
    if (!fOutputList) {
        printf("ERROR: Output list not available\n");
        return;
    }
}


