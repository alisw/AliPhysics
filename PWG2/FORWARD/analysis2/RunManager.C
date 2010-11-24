/** 
 * Script to set-up a train 
 * 
 * @param esd           ESD file 
 * @param mc            Whether to do MC or not
 * @param nEvents       Number of events
 * @param nCutBins      Bins to cut away 
 * @param correctionCut 
 *
 * @ingroup pwg2_forward_analysis_scripts
 */
void RunManager(const char* esd, Bool_t mc=kFALSE, Int_t nEvents=1000,
		Int_t nCutBins=1, Float_t correctionCut=0.1, 
		const char* mode="local")
{
  gSystem->Load("libVMC");
  gSystem->Load("libTree");
  
  gSystem->Load("libSTEERBase");
  
  gSystem->Load("libESD") ;
  gSystem->Load("libAOD") ;
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  
  gSystem->Load("libPhysics");
  gSystem->Load("libPWG0base");
  gSystem->Load("libPWG0dep");
  gSystem->Load("libPWG2forward");
  gSystem->Load("libPWG2forward2");
  
  //You can expand this chain if you have more data :-)
  TChain* chain = new TChain("esdTree");
  chain->Add(esd);
  
  //Creating the manager and handlers
  AliAnalysisManager *mgr  = new AliAnalysisManager("Analysis Train", 
						    "FMD analysis train");
  AliESDInputHandler *esdHandler = new AliESDInputHandler();
  esdHandler->SetInactiveBranches("AliESDACORDE "
				  "AliRawDataErrorLogs "
				  "CaloClusters "
				  "Cascades "
				  "EMCALCells "
				  "EMCALTrigger "
				  "Kinks "
				  "Cascades "
				  "MuonTracks "
				  "TrdTracks "
				  "CaloClusters");
  mgr->SetInputEventHandler(esdHandler);      
       
       
  // Monte Carlo handler
  // AliMCEventHandler* mcHandler = new AliMCEventHandler();
  // mgr->SetMCtruthEventHandler(mcHandler);
  // mcHandler->SetReadTR(readTR);    
  
  // AOD output handler
  AliAODHandler* aodHandler   = new AliAODHandler();
  mgr->SetOutputEventHandler(aodHandler);
  aodHandler->SetFillAOD(kTRUE);
  aodHandler->SetFillAODforRun(kTRUE);
  aodHandler->SetOutputFileName("AliAODs.root");
    
  gROOT->LoadMacro("$ALICE_ROOT/PWG2/FORWARD/analysis2/AddTaskFMD.C");
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
  AddTaskFMD(nCutBins, correctionCut);
  AddTaskPhysicsSelection(mc, kTRUE, kTRUE);
  
  // Run the analysis
    
  TStopwatch t;
  if (!mgr->InitAnalysis()) {
    Error("RunManager", "Failed to initialize analysis train!");
    return;
  }
  // Some informative output 
  mgr->PrintStatus();
  mgr->SetUseProgressBar(kTRUE);

  // Write train to file - a test 
  TDirectory* savDir = gDirectory;
  TFile* file = TFile::Open("analysisTrain.root", "RECREATE");
  mgr->Write();
  file->Close();
  savDir->cd();
  
  // Run the train 
  t.Start();
  mgr->StartAnalysis(mode, chain, nEvents);
  t.Stop();
  t.Print();
}
//
// EOF
//
