void AddTasksFlavourJet(const Int_t iCandType = 1 /*0 = D0, 1=Dstar...*/,
   const TString sCutFile = "cutsHF/D0toKpiCutsppRecVtxNoPileupRejNoEMCAL.root",
   const Double_t dJetPtCut   = 1.,
   const Double_t dJetAreaCut = 0.,
   const char *acctype = "TPC",
   const TString sRunPeriod = "LHC10b",
   const Int_t    uBeamType = 0,
   const UInt_t uTriggerMask = AliVEvent::kMB, /*for jets; the D mesons trigger is defined in the cut object*/
   const Bool_t bIsMC = kFALSE,
   const Bool_t bIsReco = kFALSE,
   const Bool_t bIsMap = kFALSE,
   TString sText=""/*completes the name of the candidate task lists*/
   )
{
   const TString sInputTrkMC  = "MCParticlesSelected";
   const TString sInputTrkRec  = "tracks";
   const TString sUsedTrks  = "PicoTracks";
   const TString sUsedClus  = "";
   TString sInputTrk = bIsReco ? sInputTrkRec : sInputTrkMC;
   const Int_t iJetAlgo = 1;
   const Int_t iJetType = 1;
   /*
   const Int_t    nRadius = 3;
   const Double_t aRadius[] = {  0.2,   0.4,   0.6  };
   const TString  sRadius[] = { "R02", "R04", "R06" };
   */
   const Int_t    nRadius = 1;
   const Double_t aRadius[] = {  0.4  };
   const TString  sRadius[] = { "R04" };
   
   //=============================================================================
   
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   
   if (!mgr) {
      ::Error("AddTasksFlavourJet.C::AddTasksFlavourJet", "No analysis manager to connect to.");
      return;
   }
   
   TString type = mgr->GetInputEventHandler()->GetDataType();
   if (!type.Contains("ESD") && !type.Contains("AOD")) {
      ::Error("AddTasksFlavourJet.C::AddTasksFlavourJet", "Task manager to have an ESD or AOD input handler.");
      return;
   }
   
   if (!mgr->GetInputEventHandler()) {
      ::Error("AddTasksFlavourJet.C::AddTasksFlavourJet", "This task requires an input event handler");
      return;
   }
   //=============================================================================
   
   UInt_t uAnaType = (((iJetType==0) ||     (iJetType==2)) ? 1 : 0);
   Int_t  iLeading =  ((iJetType==0) ? 3 : ((iJetType==1)  ? 0 : 1));
   Int_t leadHadType=0; /* 0=charged, 1=neutral, 2=both*/
   //D mesons -- PID
   gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
   AliAnalysisTaskSE *taskRespPID = AddTaskPIDResponse(bIsMC);
   
   // -- D meson selection
  
   gROOT->LoadMacro("AddTaskDFilterAndCorrelationsExch.C");

   if(bIsMap) {
      AliAnalysisTaskSEDmesonsFilterCJ *taskMCDmesonsFilter = AddTaskSEDmesonsFilterCJ(iCandType,sCutFile,bIsMC,kFALSE,sText);
   }
   // EMCal framework
   // -- Physics selection task
   gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskEmcalPhysicsSelection.C");
   AliPhysicsSelectionTask *physSelTask = AddTaskEmcalPhysicsSelection(kTRUE, kTRUE, uTriggerMask, 5, 5, 10, kTRUE, -1, -1, -1, -1);
   
   if (!physSelTask) {
      cout << "no physSelTask"; 
      return; 
   }
   
   // -- 
   gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskEmcalSetup.C");
   AliEmcalSetupTask *taskSetupEMCal = AddTaskEmcalSetup();
   taskSetupEMCal->SetGeoPath("$ALICE_ROOT/OADB/EMCAL");
   taskSetupEMCal->SelectCollisionCandidates(uTriggerMask);
   
   // Jet preparation
   //gROOT->LoadMacro("/data/Work/jets/testEMCalJetFramework/code/v4/AddTaskJetPreparation.C");
   gROOT->LoadMacro("$ALICE_ROOT/PWGJE/EMCALJetTasks/macros/AddTaskJetPreparation.C");
   AddTaskJetPreparation(sRunPeriod,"PicoTracks",bIsMC ? "MCParticlesSelected" : "",/*next 7 emcal default settings*/"","",2.,0.,0.03,0.015,0.15,uTriggerMask, kFALSE /*track cluster*/,kFALSE /*do histos*/,kTRUE /*make pico tracks*/,kFALSE /*make emcal trigger*/,kFALSE /*is emcal train*/);
   
   gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskEmcalPicoTrackMaker.C");
   AliEmcalPicoTrackMaker *taskPicoTrack = AddTaskEmcalPicoTrackMaker(sUsedTrks.Data(),sInputTrk.Data());
   taskPicoTrack->SelectCollisionCandidates(uTriggerMask);
   
   gROOT->LoadMacro("$ALICE_ROOT/PWGJE/EMCALJetTasks/macros/AddTaskEmcalJet.C");
   //gROOT->LoadMacro("$ALICE_ROOT/PWGJE/EMCALJetTasks/macros/AddTaskEmcalJetSample.C");
   
   gROOT->LoadMacro("$ALICE_ROOT/PWGJE/EMCALJetTasks/macros/AddTaskJetResponseMaker.C");
   
   for (Int_t i=0; i<nRadius; i++) {
      //jet reconstruction for correlation
      //AliEmcalJetTask *taskFJ = AddTaskEmcalJet(bIsMC&&!bIsReco ? sInputTrkMC.Data() : sUsedTrks.Data(),sUsedClus.Data(),iJetAlgo,aRadius[i],iJetType);
      //I think like follows it's fine: sUsedTrks contains the picotracks of MC or reco tracks according to IsReco
      AliEmcalJetTask *taskFJ = AddTaskEmcalJet(sUsedTrks.Data(),sUsedClus.Data(),iJetAlgo,aRadius[i],iJetType);
      
      taskFJ->SelectCollisionCandidates(uTriggerMask);
      
            
     // AliAnalysisTaskEmcalJetSample *taskjetsample=AddTaskEmcalJetSample(sUsedTrks.Data(),"",taskFJ->GetName(),
      	// "Rho",aRadius[i],dJetPtCut,dJetAreaCut,AliAnalysisTaskEmcalJet::kTPC,iLeading,
      	// "AliAnalysisTaskEmcalJetSample");
      //taskjetsample->SelectCollisionCandidates(uTriggerMask);
      
      //Filter and correlation with D meson
      
      AddTaskDFilterAndCorrelationsExch(
	 iCandType,
      	 sCutFile,
      	 bIsMC,
      	 bIsReco,
	 "",
      	 taskFJ->GetName(),
      	 //Form("JetR%s",sRadius[i].Data()),
      	 iLeading,
      	 leadHadType,
      	 aRadius[i],
      	 dJetPtCut,
      	 acctype
      	//percjetareacut=1.
	  );
      
      AliEmcalJetTask *taskMCJ;
      //jet reconstruction for correction map
      if(bIsMap){ 
      	 taskMCJ = AddTaskEmcalJet(sInputTrkMC.Data(),sUsedClus.Data(),iJetAlgo,aRadius[i],
      	    iJetType);
      	 
      	 AliAnalysisTaskFlavourJetCorrelations *taskMCDmesonCJ = AddTaskFlavourJetCorrelations(
      	    iCandType,
      	    sCutFile,
      	    kTRUE,
      	    kFALSE,
      	    taskMCJ->GetName(),
      	    Form("JetR%s",sRadius[i].Data()),
      	    iLeading,
      	    leadHadType,
      	    aRadius[i],
      	    dJetPtCut,
      	    acctype
      	    /*percjetareacut=1.*/);
      	 
      	 taskMCDmesonCJ->SetName(Form("AliAnalysisTaskSEEmcalJetMCDmesonsCJ_%s",sRadius[i].Data()));
      	 taskMCDmesonCJ->SetForceBeamType(uBeamType);
      	 taskMCDmesonCJ->SetAnaType(uAnaType);
      	 taskMCDmesonCJ->SetLeadingHadronType(iLeading);
      	 //  taskDmesonCJ->SelectCollisionCandidates(uTriggerMask);
      	 
      	 // definition of correction map
      	 Int_t tag=0;
      	 if(iCandType == 0) tag=AliEmcalJet::kD0;
      	 if(iCandType == 1) tag=AliEmcalJet::kDStar;
      	 Printf("************** tag = %d", tag&0x1);
      	 AliJetResponseMaker* taskResp=AddTaskJetResponseMaker(
      	    sUsedTrks.Data(),sUsedClus.Data(),taskFJ->GetName(),"",aRadius[i],
      	    sInputTrkMC.Data(),"",taskMCJ->GetName(),"",aRadius[i],dJetPtCut,dJetAreaCut,5,0,AliJetResponseMaker::kGeometrical, 0.25,0.25,"TPC",-999,-999,-999,"AliJetResponseMaker", kFALSE, 0, -10,10, tag );
      	 taskResp->SetMinJetMCPt(0); //added to bypass a return not needed (feature of PrepareJetTask)
      	 //taskResp->SetHistoType(1);
      	 
      }
   }
   
   return;
}
