//Based on the file $ALICE_ROOT/PWGLF/STRANGENESS/Cascades/macros/runProofCascadePbPb()
void runProofLambdaOverK0sJets(TString  proofCluster  = "xsanchez@skaf.saske.sk",
			       TString  alirootVer    = "VO_ALICE@AliRoot::v5-03-70-AN",
			       TString  rootVer       = "VO_ALICE@ROOT::v5-34-02", 
			       TString  path          = "/alice/data/LHC10h_000138624_p2_AOD049",
			       TString  name          = "LambdaOverK0sRatio", 
			       Int_t    data          = 2010,
			       Float_t  minCen        = 0.,
			       Float_t  maxCen        = 90.,
			       Float_t  ptMinTrig     = 5.,
			       Float_t  ptMaxTrig     = 10.,
			       Float_t  etaMaxTrig    = 0.75,
			       Float_t  rapMaxV0      = 0.75,
			       Bool_t   sepInjec      = kTRUE,
			       Bool_t   isMC          = kFALSE,
			       Bool_t   usePID        = kTRUE,
			       Bool_t   doQA          = kFALSE){
  
  Printf("   \nThe parameters of the programm are : \n ");
  Printf(" \t Analysis mode:\t %s\n \t Centrality:\t %.1lf - %.1lf\n \t Use MC Data?:\t %s\n \t Use PID?:\t %s\n",
	 "Proof", minCen,maxCen,
	 (isMC) ? "Yes" : "No",
	 (usePID) ? "Yes" : "No");
  
  // _____________________________________________________ // 

  gEnv->SetValue("XSec.GSI.DelegProxy", "2");
    
  TProof::Mgr(proofCluster.Data())->SetROOTVersion(rootVer.Data());

  TString alirootMode = ""; 
  TString extraLibs;
  TList *list = new TList();
  alirootMode="ALIROOT";
  extraLibs+= "ANALYSIS:OADB:ANALYSISalice:CORRFW";  
  list->Add(new TNamed("ALIROOT_MODE", alirootMode.Data()));
  list->Add(new TNamed("ALIROOT_EXTRA_LIBS", extraLibs.Data()));
   
  TProof::Reset(proofCluster.Data());
  TProof::Open(proofCluster.Data()); 
  //TProof::Open(proofCluster.Data(),"workers=1"); 
  gProof->ClearPackages();
  gProof->EnablePackage(alirootVer.Data(),list);

  // _____________________________________________________ //

  gSystem->AddIncludePath("-I$ALICE_ROOT/include");
  
  AliAnalysisManager *mgr = new AliAnalysisManager("Manager");
  
  AliAODInputHandler* aodH = new AliAODInputHandler;
  mgr->SetInputEventHandler(aodH);
   
  //PID
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
  AliAnalysisTask *pidTask = AddTaskPIDResponse(isMC);
  //AliAnalysisTask *pidTask = AddTaskPIDResponse(isMC,kTRUE);
  if(!pidTask) { printf("no PIDtask\n"); return; }

  Float_t checkIDTrig= kTRUE;
    
  // My task
  gROOT->LoadMacro("AliAnalysisTaskLambdaOverK0sJets.cxx+g"); 
  gROOT->LoadMacro("AddTaskLambdaOverK0sJets.C");
  AliAnalysisTaskLambdaOverK0sJets *task = AddTaskLambdaOverK0sJets(name,data,minCen,maxCen,ptMinTrig,ptMaxTrig,etaMaxTrig,checkIDTrig,rapMaxV0,sepInjec,isMC,usePID,doQA);
  
  // _____________________________________________________ //
  
  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  mgr->StartAnalysis("proof",path); 
  //mgr->StartAnalysis("proof",path,1,1); 
  
}

