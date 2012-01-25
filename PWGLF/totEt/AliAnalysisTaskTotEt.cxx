//_________________________________________________________________________
//  Utility Class for transverse energy studies
//  Task for analysis
//  - reconstruction and MC output
// implementation file
//
//*-- Authors: Oystein Djuvsland (Bergen), David Silvermyr (ORNL)
	       //_________________________________________________________________________
	       //Necessary to read config macros
#include <TROOT.h>
#include <TSystem.h>
#include <TInterpreter.h>

#include "TChain.h"
#include "TList.h"
#include "TFile.h"
#include "TH2F.h"
#include "THnSparse.h"

#include "AliESDEvent.h"
#include "AliMCEvent.h"
#include "AliESDtrackCuts.h"

#include "AliAnalysisTaskTotEt.h"
#include "AliAnalysisEtReconstructedPhos.h"
#include "AliAnalysisEtReconstructedEmcal.h"
#include "AliAnalysisEtMonteCarloPhos.h"
#include "AliAnalysisEtMonteCarloEmcal.h"
#include "AliAnalysisEmEtMonteCarlo.h"
#include "AliAnalysisEmEtReconstructed.h"

#include <iostream>
#include <AliCentrality.h>

  using namespace std;

ClassImp(AliAnalysisTaskTotEt)

//________________________________________________________________________
  AliAnalysisTaskTotEt::AliAnalysisTaskTotEt(const char *name, Bool_t isMc) :
    AliAnalysisTaskTransverseEnergy(name, isMc)
    ,fRecAnalysis(0)
    ,fMCAnalysis(0)
					    // ,fSparseHistRecVsMc(0)
					    //,fSparseRecVsMc(0)
{
  // Constructor
  // select if we should use EMCal or PHOS class
  // PHOS by default, EMCal if name string contains EMC
  TString t(name);
  //t.ToUpper();
  if (t.Contains("EMC")) {
    if (t.Contains("Detail")) {

      cout<<"Rereading AliAnalysisEtMonteCarlo configuration file..."<<endl;
      gROOT->LoadMacro(fMCConfigFile);
      fMCAnalysis = (AliAnalysisEmEtMonteCarlo *) gInterpreter->ProcessLine("ConfigEtMonteCarlo(true,true)");
			
      cout << "Instantiating AliAnalysisEmEtMonteCarlo class..."<< endl;
    }
    else if (fMCConfigFile.Length()) {
      cout<<"Rereading AliAnalysisEtMonteCarloEmcal configuration file..."<<endl;
      gROOT->LoadMacro(fMCConfigFile);
      fMCAnalysis = (AliAnalysisEtMonteCarloEmcal *) gInterpreter->ProcessLine("ConfigEtMonteCarlo()");
    }
		
    if (t.Contains("Detail")) {
      cout<<"Rereading AliAnalysisEmEtReconstructed configuration file..."<<endl;
      gROOT->LoadMacro(fRecoConfigFile);
      fRecAnalysis = (AliAnalysisEmEtReconstructed *) gInterpreter->ProcessLine("ConfigEtReconstructed(true,true)");

    }
    else if (fRecoConfigFile.Length()) {
      cout<<"Rereading AliAnalysisEtReconstructedEmcal configuration file..."<<endl;
      gROOT->LoadMacro(fRecoConfigFile);
      fRecAnalysis = (AliAnalysisEtReconstructedEmcal *) gInterpreter->ProcessLine("ConfigEtReconstructed()");
    }
  }
  else {
    if (fMCConfigFile.Length()) {
      cout<<"Rereading AliAnalysisEtMonteCarloPhos configuration file..."<<endl;
      gROOT->LoadMacro(fMCConfigFile);
			
      fMCAnalysis = (AliAnalysisEtMonteCarloPhos *) gInterpreter->ProcessLine("ConfigEtMonteCarlo(false)");
      cout << fMCAnalysis << endl;
    }
		
    if (fRecoConfigFile.Length()) {
      cout<<"Rereading AliAnalysisEtReconstructedPhos configuration file..."<<endl;
      gROOT->LoadMacro(fRecoConfigFile);
      fRecAnalysis = (AliAnalysisEtReconstructedPhos *) gInterpreter->ProcessLine("ConfigEtReconstructed(false)");
    }
  }
  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #1 writes into a TH1 container
	
  DefineOutput(1, TList::Class());
	
}
AliAnalysisTaskTotEt::~AliAnalysisTaskTotEt() {//Destructor
  //    fOutputList->Clear();
  delete fRecAnalysis;
  delete fMCAnalysis;
  //delete fSparseHistRecVsMc;
  //delete fSparseRecVsMc;
}

//________________________________________________________________________
void AliAnalysisTaskTotEt::UserCreateOutputObjects()
{
  // Create histograms
  // Called once
  if (fMCAnalysis)
    fMCAnalysis->CreateHistograms();
  fRecAnalysis->CreateHistograms();
  fOutputList = new TList;
  fOutputList->SetOwner();
  fRecAnalysis->FillOutputList(fOutputList);
  if (fMCAnalysis)
    fMCAnalysis->FillOutputList(fOutputList);
  fHistEtRecvsEtMC = new TH2F("fHistEtRecvsEtMC", "Reconstructed E_{T} vs MC E_{T}", 1000, 0.000, 100, 1000, 0.0001, 100);
  fHistEtRecOverEtMC = new TH2F("fHistEtRecOverEtMC", "Reconstructed E_{T} over MC E_{T} vs centrality", 1000, 0.00, 2.0, 11, -0.5, 10.5); 
  fHistDiffEtRecEtMCOverEtMC = new TH2F("fHistDiffEtRecEtMCOverEtMC", "fHistDiffEtRecEtMCOverEtMC", 10000, 0.0, 1000, 1000, -5, 5); 
  fOutputList->Add(fHistEtRecvsEtMC);
  fOutputList->Add(fHistEtRecOverEtMC);
  fOutputList->Add(fHistDiffEtRecEtMCOverEtMC);

  Bool_t selectPrimaries=kTRUE;
  if(fRecAnalysis->DataSet()==2010 || fRecAnalysis->DataSet()==20111||fRecAnalysis->DataSet()==2009){
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
  if (fEsdtrackCutsITSTPC && fEsdtrackCutsTPC) {
    fRecAnalysis->SetITSTrackCuts( GetITSTrackCuts());
    if (fMCAnalysis)
      fMCAnalysis->SetITSTrackCuts( GetITSTrackCuts());
    fRecAnalysis->SetTPCITSTrackCuts( GetTPCITSTrackCuts());
    if (fMCAnalysis)
      fMCAnalysis->SetTPCITSTrackCuts( GetTPCITSTrackCuts());
    fRecAnalysis->SetTPCOnlyTrackCuts( GetTPCOnlyTrackCuts());
    if (fMCAnalysis)
      fMCAnalysis->SetTPCOnlyTrackCuts( GetTPCOnlyTrackCuts());
    //add ITS stuff!
  }
  else {
    Printf("Error: no track cuts!");
  }

}

//________________________________________________________________________
void AliAnalysisTaskTotEt::UserExec(Option_t *)
{ // execute method
	
  fESDEvent = dynamic_cast<AliESDEvent*>(InputEvent());
  if (!fESDEvent)
    {
      Printf("ERROR: Could not retrieve event");
      return;
    }
	
  //Int_t res = CheckPhysicsSelection(fESDEvent->GetRunNumber());
	
  AliCentrality *cent = GetCentralityObject();
	
  //if (res == 0 && cent)
  //{
  if (IsPhysicsSelected())
    {
      fRecAnalysis->SetCentralityObject(cent);
      fRecAnalysis->AnalyseEvent(fESDEvent);
		
      AliMCEvent* mcEvent = MCEvent();
      if (mcEvent)
	{
	  fMCAnalysis->SetCentralityObject(cent);
	  fMCAnalysis->AnalyseEvent(mcEvent, fESDEvent);
	  //fMCAnalysis->AnalyseEvent(mcEvent);
	}
      if(fMCAnalysis)
	{
	  fHistEtRecvsEtMC->Fill(fRecAnalysis->GetTotNeutralEtAcc(), fMCAnalysis->GetTotNeutralEtAcc());
	  if(fMCAnalysis->GetTotNeutralEtAcc()) fHistEtRecOverEtMC->Fill(fRecAnalysis->GetTotNeutralEt()/fMCAnalysis->GetTotNeutralEtAcc(), cent->GetCentralityClass10("V0M"));
	  if(fMCAnalysis->GetTotNeutralEtAcc()) fHistDiffEtRecEtMCOverEtMC->Fill(fMCAnalysis->GetTotNeutralEt(), (fRecAnalysis->GetTotNeutralEt()-fMCAnalysis->GetTotNeutralEt())/fMCAnalysis->GetTotNeutralEt());
	}
    }
  //}
  // Post output data.
  PostData(1, fOutputList);
}

//________________________________________________________________________
void AliAnalysisTaskTotEt::Terminate(Option_t *)
{
  // Draw result to the screen
  // Called once at the end of the query
	
  fOutputList = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutputList) {
    printf("ERROR: Output list not available\n");
    return;
  }
}



