AliAnalysisTaskTrackingSysPropagation *AddTaskTrackingSysPropagation(Int_t system=0 /*0=pp,1=PbPb*/,
                                                                     AliAnalysisTaskTrackingSysPropagation::DecChannel ch=AliAnalysisTaskTrackingSysPropagation::kD0toKpi,
                                                                     TString filename = "",
                                                                     TString filenameHistME = "",
                                                                     TString filenameHistTrEff = "",
                                                                     TString postname = "",
                                                                     Double_t maxPt=60.,
                                                                     TString filenameHistMEPr = "")
{
    
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTask", "No analysis manager to connect to.");
  }
    
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskImpParDistrib", "This task requires an input event handler");
    return NULL;
  }
    
    
  Bool_t stdcuts=kFALSE;
  TFile* filecuts;
  if( filename.EqualTo("") ) {
    stdcuts=kTRUE;
  } else {
    filecuts=TFile::Open(filename.Data());
    if(!filecuts ||(filecuts&& !filecuts->IsOpen())){
      Printf("FATAL: Cut object not found: analysis will not start!\n"); return NULL;
    }
    else printf("Cut object correctly found\n");
  }
    
  //Analysis Task
    
  AliRDHFCuts* analysiscuts;
    
  if(ch==AliAnalysisTaskTrackingSysPropagation::kDstartoKpipi) analysiscuts = (AliRDHFCuts*)filecuts->Get("DStartoKpipiCuts");
  else if(ch==AliAnalysisTaskTrackingSysPropagation::kD0toKpi) analysiscuts = (AliRDHFCuts*)filecuts->Get("D0toKpiCuts");
  else if(ch==AliAnalysisTaskTrackingSysPropagation::kLctopK0s) analysiscuts = (AliRDHFCuts*)filecuts->Get("LctoK0spCuts");
  else if(ch==AliAnalysisTaskTrackingSysPropagation::kLctopKpi) analysiscuts = (AliRDHFCuts*)filecuts->Get("LctopKpiAnalysisCuts");
  else if(ch==AliAnalysisTaskTrackingSysPropagation::kLctopKpiFromSc) analysiscuts = (AliRDHFCuts*)filecuts->Get("XictopKpiProdCuts");
   
  else analysiscuts = (AliRDHFCuts*)filecuts->Get("AnalysisCuts");
    
  TFile* fileMESys;
  if( filenameHistME.EqualTo("") ) {
    Printf("FATAL: Histo with ME syst. not found: analysis will not start!\n"); return NULL;
  } else {
    fileMESys=TFile::Open(filenameHistME.Data());
    if(!fileMESys ||(fileMESys&& !fileMESys->IsOpen())){
      Printf("FATAL: Histo with ME syst. not found: analysis will not start!\n"); return NULL;
    }
  }
  TH1F *histoME = (TH1F*)fileMESys->Get("h");
  if(!histoME){ Printf("FATAL: Histo with ME syst. not found: analysis will not start!\n"); return NULL;}
    
  TH1F *histoMEPr = 0x0;
  TFile* fileMESysPr;
  if(!filenameHistMEPr.EqualTo("")) {
    fileMESysPr=TFile::Open(filenameHistMEPr.Data());
    histoMEPr = (TH1F*)fileMESysPr->Get("hPr");
    if(!histoMEPr){ Printf("FATAL: Histo with ME syst. for proton not found: analysis will not start!\n"); return NULL;}
  }
  
  TFile* fileTrEffSys;
  if( filenameHistTrEff.EqualTo("") ) {
    Printf("FATAL: Histo with TrEff syst. not found: analysis will not start!\n"); return NULL;
  } else {
    fileTrEffSys=TFile::Open(filenameHistTrEff.Data());
    if(!fileTrEffSys ||(fileTrEffSys&& !fileTrEffSys->IsOpen())){
      Printf("FATAL: Histo with TrEff syst. not found: analysis will not start!\n"); return NULL;
    }
  }
  TH1F *histoTrEff = (TH1F*)fileTrEffSys->Get("hTrEff");
  if(!histoTrEff) {Printf("FATAL: Histo TrEff not found: analysis will not start!\n"); return NULL;}
    
    
  AliAnalysisTaskTrackingSysPropagation *Task = new AliAnalysisTaskTrackingSysPropagation(ch, analysiscuts, histoME, histoTrEff, histoMEPr);
    
  Task->SetMaximumPt(maxPt);
  Task->SetDebugLevel(1);
  mgr->AddTask(Task);
    
  // Create containers for input/output
  TString outputFileName = AliAnalysisManager::GetCommonFileName();
  outputFileName += Form(":MatchEff_%s",postname.Data());
    
  AliAnalysisDataContainer *coutput =0x0;
  coutput = mgr->CreateContainer(Form("systUnc_%s",postname.Data()),
				 TList::Class(),
				 AliAnalysisManager::kOutputContainer,
				 outputFileName );
    
    
  mgr->ConnectInput(Task,0,mgr->GetCommonInputContainer());
  mgr->ConnectOutput(Task,1,coutput);
    
  return Task;
}
