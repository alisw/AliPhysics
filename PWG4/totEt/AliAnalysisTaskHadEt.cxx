//_________________________________________________________________________
//  Utility Class for transverse energy studies; charged hadrons
//  Task for analysis
//  - reconstruction and MC output
// implementation file
//
//Created by Christine Nattrass, Rebecca Scott, Irakli Martashvili
//University of Tennessee at Knoxville
//_________________________________________________________________________

#include "TChain.h"
#include "TList.h"
#include "TH2F.h"

#include "AliESDEvent.h"
#include "AliMCEvent.h"
#include "AliESDtrackCuts.h"

#include "AliAnalysisTaskHadEt.h"
#include "AliAnalysisHadEtReconstructed.h"
#include "AliAnalysisHadEtMonteCarlo.h"

#include <iostream>

using namespace std;

ClassImp(AliAnalysisTaskHadEt)



//________________________________________________________________________
AliAnalysisTaskHadEt::AliAnalysisTaskHadEt(const char *name) :
        AliAnalysisTaskSE(name)
        ,fOutputList(0)
        ,fRecAnalysis(0)
        ,fMCAnalysis(0)
        ,fHistEtRecvsEtMC(0)
        ,fEsdtrackCutsITSTPC(0)
        ,fEsdtrackCutsTPC(0)
        ,fEsdtrackCutsITS(0)
{
    // Constructor

    fMCAnalysis = new AliAnalysisHadEtMonteCarlo();
    fMCAnalysis->Init();

    fRecAnalysis = new AliAnalysisHadEtReconstructed();
    fRecAnalysis->Init();

    // Define input and output slots here
    // Input slot #0 works with a TChain
    DefineInput(0, TChain::Class());
    // Output slot #1 writes into a TH1 container

    DefineOutput(1, TList::Class());
}
AliAnalysisTaskHadEt::~AliAnalysisTaskHadEt(){//Destructor
  fOutputList->Clear();
  delete fOutputList;
  delete fRecAnalysis;
  delete fMCAnalysis;
  delete fEsdtrackCutsITSTPC;
  delete fEsdtrackCutsTPC;
  delete fEsdtrackCutsITS;
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
  fMCAnalysis->CreateHistograms();
  fRecAnalysis->CreateHistograms();
  fRecAnalysis->FillOutputList();
  fMCAnalysis->FillOutputList();


    Bool_t selectPrimaries=kTRUE;
    fEsdtrackCutsITSTPC = AliESDtrackCuts::GetStandardITSTPCTrackCuts2009(selectPrimaries);
    fEsdtrackCutsITSTPC->SetName("fEsdTrackCuts");
    fEsdtrackCutsTPC = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
    fEsdtrackCutsTPC->SetName("fEsdTrackCutsTPCOnly");
    //ITS stand alone cuts - similar to 2009 cuts but with only ITS hits required
    fEsdtrackCutsITS =  AliESDtrackCuts::GetStandardITSPureSATrackCuts2009(kTRUE,kFALSE);//we do want primaries but we do not want to require PID info
    fEsdtrackCutsITS->SetName("fEsdTrackCutsITS");

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
  AliESDEvent *event = dynamic_cast<AliESDEvent*>(InputEvent());
if (!event) {
  Printf("ERROR: Could not retrieve event");
  return;
 }
fRecAnalysis->AnalyseEvent(event);

AliMCEvent* mcEvent = MCEvent();
if (mcEvent)
  {
    ((AliAnalysisHadEtMonteCarlo*)fMCAnalysis)->AnalyseEvent((AliVEvent*)mcEvent,(AliVEvent*)event);
    //cout<<"Simulated Hadronic Et "<<fMCAnalysis->GetSimulatedHadronicEt()<<" Reconstructed Hadronic Et "<<fRecAnalysis->GetCorrectedHadEtFullAcceptanceITS()<<endl;
    //cout<<"Simulated Total Et "<<fMCAnalysis->GetSimulatedTotalEt()<<" Reconstructed Total Et "<<fRecAnalysis->GetCorrectedTotEtFullAcceptanceITS()<<endl;
    fMCAnalysis->FillSimTotEtVsRecoTotEtFullAcceptanceTPC( fRecAnalysis->GetCorrectedTotEtFullAcceptanceTPC() );
    fMCAnalysis->FillSimTotEtVsRecoTotEtFullAcceptanceITS( fRecAnalysis->GetCorrectedTotEtFullAcceptanceITS() );
//     fMCAnalysis->FillSimTotEtVsRecoTotEtEMCALAcceptanceTPC( fRecAnalysis->GetCorrectedTotEtEMCALAcceptanceTPC() );
//     fMCAnalysis->FillSimTotEtVsRecoTotEtEMCALAcceptanceITS( fRecAnalysis->GetCorrectedTotEtEMCALAcceptanceITS() );
//     fMCAnalysis->FillSimTotEtVsRecoTotEtPHOSAcceptanceTPC( fRecAnalysis->GetCorrectedTotEtPHOSAcceptanceTPC() );
//     fMCAnalysis->FillSimTotEtVsRecoTotEtPHOSAcceptanceITS( fRecAnalysis->GetCorrectedTotEtPHOSAcceptanceITS() );
    fMCAnalysis->FillSimTotEtVsRecoTotEtFullAcceptanceTPCNoPID( fRecAnalysis->GetCorrectedTotEtFullAcceptanceTPCNoPID() );
    fMCAnalysis->FillSimTotEtVsRecoTotEtFullAcceptanceITSNoPID( fRecAnalysis->GetCorrectedTotEtFullAcceptanceITSNoPID() );
//     fMCAnalysis->FillSimTotEtVsRecoTotEtEMCALAcceptanceTPCNoPID( fRecAnalysis->GetCorrectedTotEtEMCALAcceptanceTPCNoPID() );
//     fMCAnalysis->FillSimTotEtVsRecoTotEtEMCALAcceptanceITSNoPID( fRecAnalysis->GetCorrectedTotEtEMCALAcceptanceITSNoPID() );
//     fMCAnalysis->FillSimTotEtVsRecoTotEtPHOSAcceptanceTPCNoPID( fRecAnalysis->GetCorrectedTotEtPHOSAcceptanceTPCNoPID() );
//     fMCAnalysis->FillSimTotEtVsRecoTotEtPHOSAcceptanceITSNoPID( fRecAnalysis->GetCorrectedTotEtPHOSAcceptanceITSNoPID() );
    fMCAnalysis->FillSimHadEtVsRecoHadEtFullAcceptanceTPC( fRecAnalysis->GetCorrectedHadEtFullAcceptanceTPC() );
    fMCAnalysis->FillSimHadEtVsRecoHadEtFullAcceptanceITS( fRecAnalysis->GetCorrectedHadEtFullAcceptanceITS() );
//     fMCAnalysis->FillSimHadEtVsRecoHadEtEMCALAcceptanceTPC( fRecAnalysis->GetCorrectedHadEtEMCALAcceptanceTPC() );
//     fMCAnalysis->FillSimHadEtVsRecoHadEtEMCALAcceptanceITS( fRecAnalysis->GetCorrectedHadEtEMCALAcceptanceITS() );
//     fMCAnalysis->FillSimHadEtVsRecoHadEtPHOSAcceptanceTPC( fRecAnalysis->GetCorrectedHadEtPHOSAcceptanceTPC() );
//     fMCAnalysis->FillSimHadEtVsRecoHadEtPHOSAcceptanceITS( fRecAnalysis->GetCorrectedHadEtPHOSAcceptanceITS() );
    fMCAnalysis->FillSimHadEtVsRecoHadEtFullAcceptanceTPCNoPID( fRecAnalysis->GetCorrectedHadEtFullAcceptanceTPCNoPID() );
    fMCAnalysis->FillSimHadEtVsRecoHadEtFullAcceptanceITSNoPID( fRecAnalysis->GetCorrectedHadEtFullAcceptanceITSNoPID() );
//     fMCAnalysis->FillSimHadEtVsRecoHadEtEMCALAcceptanceTPCNoPID( fRecAnalysis->GetCorrectedHadEtEMCALAcceptanceTPCNoPID() );
//     fMCAnalysis->FillSimHadEtVsRecoHadEtEMCALAcceptanceITSNoPID( fRecAnalysis->GetCorrectedHadEtEMCALAcceptanceITSNoPID() );
//     fMCAnalysis->FillSimHadEtVsRecoHadEtPHOSAcceptanceTPCNoPID( fRecAnalysis->GetCorrectedHadEtPHOSAcceptanceTPCNoPID() );
//     fMCAnalysis->FillSimHadEtVsRecoHadEtPHOSAcceptanceITSNoPID( fRecAnalysis->GetCorrectedHadEtPHOSAcceptanceITSNoPID() );


    fMCAnalysis->FillSimTotEtMinusRecoTotEtFullAcceptanceTPC( fRecAnalysis->GetCorrectedTotEtFullAcceptanceTPC() );
    fMCAnalysis->FillSimTotEtMinusRecoTotEtFullAcceptanceITS( fRecAnalysis->GetCorrectedTotEtFullAcceptanceITS() );
//     fMCAnalysis->FillSimTotEtMinusRecoTotEtEMCALAcceptanceTPC( fRecAnalysis->GetCorrectedTotEtEMCALAcceptanceTPC() );
//     fMCAnalysis->FillSimTotEtMinusRecoTotEtEMCALAcceptanceITS( fRecAnalysis->GetCorrectedTotEtEMCALAcceptanceITS() );
//     fMCAnalysis->FillSimTotEtMinusRecoTotEtPHOSAcceptanceTPC( fRecAnalysis->GetCorrectedTotEtPHOSAcceptanceTPC() );
//     fMCAnalysis->FillSimTotEtMinusRecoTotEtPHOSAcceptanceITS( fRecAnalysis->GetCorrectedTotEtPHOSAcceptanceITS() );
    fMCAnalysis->FillSimTotEtMinusRecoTotEtFullAcceptanceTPCNoPID( fRecAnalysis->GetCorrectedTotEtFullAcceptanceTPCNoPID() );
    fMCAnalysis->FillSimTotEtMinusRecoTotEtFullAcceptanceITSNoPID( fRecAnalysis->GetCorrectedTotEtFullAcceptanceITSNoPID() );
//     fMCAnalysis->FillSimTotEtMinusRecoTotEtEMCALAcceptanceTPCNoPID( fRecAnalysis->GetCorrectedTotEtEMCALAcceptanceTPCNoPID() );
//     fMCAnalysis->FillSimTotEtMinusRecoTotEtEMCALAcceptanceITSNoPID( fRecAnalysis->GetCorrectedTotEtEMCALAcceptanceITSNoPID() );
//     fMCAnalysis->FillSimTotEtMinusRecoTotEtPHOSAcceptanceTPCNoPID( fRecAnalysis->GetCorrectedTotEtPHOSAcceptanceTPCNoPID() );
//     fMCAnalysis->FillSimTotEtMinusRecoTotEtPHOSAcceptanceITSNoPID( fRecAnalysis->GetCorrectedTotEtPHOSAcceptanceITSNoPID() );
    fMCAnalysis->FillSimHadEtMinusRecoHadEtFullAcceptanceTPC( fRecAnalysis->GetCorrectedHadEtFullAcceptanceTPC() );
    fMCAnalysis->FillSimHadEtMinusRecoHadEtFullAcceptanceITS( fRecAnalysis->GetCorrectedHadEtFullAcceptanceITS() );
//     fMCAnalysis->FillSimHadEtMinusRecoHadEtEMCALAcceptanceTPC( fRecAnalysis->GetCorrectedHadEtEMCALAcceptanceTPC() );
//     fMCAnalysis->FillSimHadEtMinusRecoHadEtEMCALAcceptanceITS( fRecAnalysis->GetCorrectedHadEtEMCALAcceptanceITS() );
//     fMCAnalysis->FillSimHadEtMinusRecoHadEtPHOSAcceptanceTPC( fRecAnalysis->GetCorrectedHadEtPHOSAcceptanceTPC() );
//     fMCAnalysis->FillSimHadEtMinusRecoHadEtPHOSAcceptanceITS( fRecAnalysis->GetCorrectedHadEtPHOSAcceptanceITS() );
    fMCAnalysis->FillSimHadEtMinusRecoHadEtFullAcceptanceTPCNoPID( fRecAnalysis->GetCorrectedHadEtFullAcceptanceTPCNoPID() );
    fMCAnalysis->FillSimHadEtMinusRecoHadEtFullAcceptanceITSNoPID( fRecAnalysis->GetCorrectedHadEtFullAcceptanceITSNoPID() );
//     fMCAnalysis->FillSimHadEtMinusRecoHadEtEMCALAcceptanceTPCNoPID( fRecAnalysis->GetCorrectedHadEtEMCALAcceptanceTPCNoPID() );
//     fMCAnalysis->FillSimHadEtMinusRecoHadEtEMCALAcceptanceITSNoPID( fRecAnalysis->GetCorrectedHadEtEMCALAcceptanceITSNoPID() );
//     fMCAnalysis->FillSimHadEtMinusRecoHadEtPHOSAcceptanceTPCNoPID( fRecAnalysis->GetCorrectedHadEtPHOSAcceptanceTPCNoPID() );
//     fMCAnalysis->FillSimHadEtMinusRecoHadEtPHOSAcceptanceITSNoPID( fRecAnalysis->GetCorrectedHadEtPHOSAcceptanceITSNoPID() );


    fMCAnalysis->FillSimTotEtMinusRawEtFullAcceptanceTPC( fRecAnalysis->GetRawEtFullAcceptanceTPC() );
    fMCAnalysis->FillSimTotEtMinusRawEtFullAcceptanceITS( fRecAnalysis->GetRawEtFullAcceptanceITS() );
//     fMCAnalysis->FillSimTotEtMinusRawTotEtEMCALAcceptanceTPC( fRecAnalysis->GetRawEtEMCALAcceptanceTPC() );
//     fMCAnalysis->FillSimTotEtMinusRawTotEtEMCALAcceptanceITS( fRecAnalysis->GetRawEtEMCALAcceptanceITS() );
//     fMCAnalysis->FillSimTotEtMinusRawTotEtPHOSAcceptanceTPC( fRecAnalysis->GetRawEtPHOSAcceptanceTPC() );
//     fMCAnalysis->FillSimTotEtMinusRawTotEtPHOSAcceptanceITS( fRecAnalysis->GetRawEtPHOSAcceptanceITS() );

    fMCAnalysis->FillSimHadEtMinusRawEtFullAcceptanceTPC( fRecAnalysis->GetRawEtFullAcceptanceTPC() );
    fMCAnalysis->FillSimHadEtMinusRawEtFullAcceptanceITS( fRecAnalysis->GetRawEtFullAcceptanceITS() );
//     fMCAnalysis->FillSimHadEtMinusRawHadEtEMCALAcceptanceTPC( fRecAnalysis->GetRawEtEMCALAcceptanceTPC() );
//     fMCAnalysis->FillSimHadEtMinusRawHadEtEMCALAcceptanceITS( fRecAnalysis->GetRawEtEMCALAcceptanceITS() );
//     fMCAnalysis->FillSimHadEtMinusRawHadEtPHOSAcceptanceTPC( fRecAnalysis->GetRawEtPHOSAcceptanceTPC() );
//     fMCAnalysis->FillSimHadEtMinusRawHadEtPHOSAcceptanceITS( fRecAnalysis->GetRawEtPHOSAcceptanceITS() );

    fMCAnalysis->FillSimPiKPMinusRecoPiKPFullAcceptanceTPC(fRecAnalysis->GetCorrectedPiKPEtFullAcceptanceTPC());
    fMCAnalysis->FillSimPiKPMinusRecoPiKPFullAcceptanceITS(fRecAnalysis->GetCorrectedPiKPEtFullAcceptanceITS());
    fMCAnalysis->FillSimPiKPMinusRecoPiKPFullAcceptanceTPCNoPID(fRecAnalysis->GetCorrectedPiKPEtFullAcceptanceTPCNoPID());
    fMCAnalysis->FillSimPiKPMinusRecoPiKPFullAcceptanceITSNoPID(fRecAnalysis->GetCorrectedPiKPEtFullAcceptanceITSNoPID());

  }

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


