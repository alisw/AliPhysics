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

#include <iostream>

using namespace std;

ClassImp(AliAnalysisTaskHadEt)



//________________________________________________________________________
  AliAnalysisTaskHadEt::AliAnalysisTaskHadEt(const char *name, Bool_t isMc, TString recoConfigFile, TString mcConfigFile) :
        AliAnalysisTaskTransverseEnergy(name, isMc)
	,fRecAnalysis(0)
	,fMCAnalysis(0)
	,fIsSim(isMc)
{
    // Constructor
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
}


//________________________________________________________________________
void AliAnalysisTaskHadEt::UserCreateOutputObjects()
{
    // Create histograms

    // Called once



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
  if(fRecAnalysis->DataSet()==2010 || fRecAnalysis->DataSet()==20111||fRecAnalysis->DataSet()==2009){
    AliAnalysisTaskSE::	SelectCollisionCandidates(AliVEvent::kINT7 ) ;
    if(fRecAnalysis->DataSet()==2010)cout<<"Setting track cuts for the 2010 p+p collisions at 7 TeV"<<endl;
    else{
      if(fRecAnalysis->DataSet()==2009){cout<<"Setting track cuts for the 2010 p+p collisions at 900 GeV"<<endl;}
      else{cout<<"Setting track cuts for the 2011 p+p collisions at 2.76 TeV"<<endl;}
    }
    //cout<<"Warning:  Have not set 2010 track cuts yet!!"<<endl;
    fEsdtrackCutsITSTPC = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(selectPrimaries);
    fEsdtrackCutsITSTPC->SetName("fEsdTrackCuts");
    fEsdtrackCutsTPC = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
    fEsdtrackCutsTPC->SetName("fEsdTrackCutsTPCOnly");
    //ITS stand alone cuts - similar to 2009 cuts but with only ITS hits required
    fEsdtrackCutsITS =  AliESDtrackCuts::GetStandardITSPureSATrackCuts2010(kTRUE,kFALSE);//we do want primaries but we do not want to require PID info
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

Int_t res = CheckPhysicsSelection(fESDEvent->GetRunNumber()); // Check if the physics selection is valid for this run

AliCentrality *cent = GetCentralityObject();

if(res == 0 && cent){
  
  //cout<<"New Event"<<endl;  

  AliMCEvent* mcEvent = MCEvent();
  Int_t eventtype = (Int_t) AliPWG0Helper::GetEventProcessType(mcEvent->Header());
  fRecAnalysis->AnalyseEvent(fESDEvent,eventtype);

// if (!mcEvent) {
//   Printf("ERROR: Could not retrieve MC event");
//  }
  if (mcEvent && fESDEvent && fIsSim){
      ((AliAnalysisHadEtMonteCarlo*)fMCAnalysis)->AnalyseEvent((AliVEvent*)mcEvent,(AliVEvent*)fESDEvent);
      if(fMCAnalysis->Full()){
	fMCAnalysis->FillSimTotEtMinusRecoTotEtFullAcceptanceTPC( fRecAnalysis->GetCorrectedTotEtFullAcceptanceTPC() );
	fMCAnalysis->FillSimTotEtMinusRecoTotEtFullAcceptanceITS( fRecAnalysis->GetCorrectedTotEtFullAcceptanceITS() );
	fMCAnalysis->FillSimTotEtMinusRecoTotEtFullAcceptanceTPCNoPID( fRecAnalysis->GetCorrectedTotEtFullAcceptanceTPCNoPID() );
	fMCAnalysis->FillSimTotEtMinusRecoTotEtFullAcceptanceITSNoPID( fRecAnalysis->GetCorrectedTotEtFullAcceptanceITSNoPID() );
	fMCAnalysis->FillSimHadEtMinusRecoHadEtFullAcceptanceTPC( fRecAnalysis->GetCorrectedHadEtFullAcceptanceTPC() );
	fMCAnalysis->FillSimHadEtMinusRecoHadEtFullAcceptanceITS( fRecAnalysis->GetCorrectedHadEtFullAcceptanceITS() );
	fMCAnalysis->FillSimHadEtMinusRecoHadEtFullAcceptanceTPCNoPID( fRecAnalysis->GetCorrectedHadEtFullAcceptanceTPCNoPID() );
	fMCAnalysis->FillSimHadEtMinusRecoHadEtFullAcceptanceITSNoPID( fRecAnalysis->GetCorrectedHadEtFullAcceptanceITSNoPID() );
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


