/**
 * @file   QA.C
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Wed Sep 24 15:04:38 2014
 * 
 * @brief  Master script for QA train 
 * 
 * @note Do not modify this script. 
 *
 *
 * This script reads in 4 other scripts 
 *
 * - GRP.C to load the global run parameters for the selected run,
 *   such as collision system, energy, etc.
 * 
 * - AODConfig.C which defines a number of functions that return
 *   either true or false.  The tasks added depends on these functions
 *
 * - BaseConfig.C which defines some base classes 
 * 
 * - DetConfig.C which defines which detectors are active and on. 
 * 
 * Users can customize QAConfig.C and DetConfig.C according to their
 * needs
 */
// Trigger mask.
UInt_t kTriggerInt        = AliVEvent::kAnyINT;
UInt_t kTriggerMuonAll    = (AliVEvent::kMUL7       | 
			     AliVEvent::kMUSH7      | 
			     AliVEvent::kMUU7       |
			     AliVEvent::kMUS7       |
			     AliVEvent::kMUSPB      |
			     AliVEvent::kMUSHPB     |
			     AliVEvent::kMuonLikePB | 
			     AliVEvent::kMuonUnlikePB);
UInt_t kTriggerMuonBarell = AliVEvent::kMUU7;
UInt_t kTriggerEMC        = (AliVEvent::kEMC7   | 
			     AliVEvent::kEMC8   | 
			     AliVEvent::kEMCEJE | 
			     AliVEvent::kEMCEGA);
UInt_t kTriggerHM         = AliVEvent::kHighMult;
UInt_t kTriggerMask       = kTriggerInt;

/**
 * Interface (pure virtual) that all configuration classes must
 * implement.
 */
struct VirtualQACfg
{
  /** @return  */
  virtual Bool_t DoCDBconnect()  const = 0;
  /** @return  */
  virtual Bool_t DoEventStat()   const = 0;
  /** @return  */
  virtual Bool_t DoCentrality()  const = 0;
  /** @return  */
  virtual Bool_t DoQAsym()       const = 0;
  /** @return  there is a 2nd file */
  virtual Bool_t DoVZERO()       const = 0;
  /** @return  */
  virtual Bool_t DoVZEROPbPb()   const = 0;
  /** @return  */
  virtual Bool_t DoVertex()      const = 0;
  /** @return  needs RP    */
  virtual Bool_t DoSPD()         const = 0;
  /** @return  */
  virtual Bool_t DoTPC()         const = 0;
  /** @return  */
  virtual Bool_t DoHLT()         const = 0;
  /** @return  needs RP */
  virtual Bool_t DoSDD()         const = 0;
  /** @return  */
  virtual Bool_t DoSSDdEdx()     const = 0;
  /** @return  */
  virtual Bool_t DoTRD()         const = 0;
  /** @return  */
  virtual Bool_t DoITS()         const = 0;
  /** @return  */
  virtual Bool_t DoITSsaTracks() const = 0;
  /** @return  */
  virtual Bool_t DoITSalign()    const = 0;
  /** @return  */
  virtual Bool_t DoCALO()        const = 0;
  /** @return  */
  virtual Bool_t DoMUONTrig()    const = 0;
  /** @return  */
  virtual Bool_t DoImpParRes()   const = 0;
  /** @return  */
  virtual Bool_t DoMUON()        const = 0;
  /** @return  */
  virtual Bool_t DoTOF()         const = 0;
  /** @return  */
  virtual Bool_t DoHMPID()       const = 0;
  /** @return  */
  virtual Bool_t DoT0()          const = 0;
  /** @return  */
  virtual Bool_t DoZDC()         const = 0;
  /** @return  */
  virtual Bool_t DoPIDResponse() const = 0;
  /** @return  */
  virtual Bool_t DoPIDqa()       const = 0;
  /** @return  */
  virtual Bool_t DoFWD()         const = 0;
  /** @return  */
  virtual Bool_t DoPHOS()        const = 0;
  /** @return  */
  virtual Bool_t DoPHOSTrig()    const = 0;
  /** @return  */
  virtual Bool_t DoEMCAL()       const = 0;
  /** @return  */
  virtual Bool_t DoFBFqa()       const = 0;
  /** @return  NEEDS geometry */
  virtual Bool_t DoMUONEff()     const = 0;
  /** @return  NEEDS MCtruth  */
  virtual Bool_t DoV0()          const = 0;
  /** @return Get Debug level */
  virtual Int_t DebugLevel() const = 0;

  virtual void PrintOne(const char* title, Bool_t use) const 
  {
    Printf("%-30s : %3s", title, use ? "yes" : "no");
  }
  virtual void Print() const 
  {
    PrintOne("CDBconnect ",  		DoCDBconnect());
    PrintOne("EventStat ",   		DoEventStat());
    PrintOne("Centrality ",  		DoCentrality());
    PrintOne("QAsym ",       		DoQAsym());
    PrintOne("VZERO",       		DoVZERO());
    PrintOne("VZEROPbPb ",   		DoVZEROPbPb());
    PrintOne("Vertex ",      		DoVertex());
    PrintOne("SPD  needs RP   ",        DoSPD());
    PrintOne("TPC ",         		DoTPC());
    PrintOne("HLT ",         		DoHLT());
    PrintOne("SDD  needs RP",         	DoSDD());
    PrintOne("SSDdEdx ",     		DoSSDdEdx());
    PrintOne("TRD ",         		DoTRD());
    PrintOne("ITS ",         		DoITS());
    PrintOne("ITSsaTracks ", 		DoITSsaTracks());
    PrintOne("ITSalign ",    		DoITSalign());
    PrintOne("CALO ",        		DoCALO());
    PrintOne("MUONTrig ",    		DoMUONTrig());
    PrintOne("ImpParRes ",   		DoImpParRes());
    PrintOne("MUON ",        		DoMUON());
    PrintOne("TOF ",         		DoTOF());
    PrintOne("HMPID ",       		DoHMPID());
    PrintOne("T0 ",          		DoT0());
    PrintOne("ZDC ",         		DoZDC());
    PrintOne("PIDResponse ", 		DoPIDResponse());
    PrintOne("PIDqa ",       		DoPIDqa());
    PrintOne("FWD ",         		DoFWD());
    PrintOne("PHOS ",        		DoPHOS());
    PrintOne("PHOSTrig ",    		DoPHOSTrig());
    PrintOne("EMCAL ",       		DoEMCAL());
    PrintOne("FBFqa ",       		DoFBFqa());
    PrintOne("MUONEff NEEDS geometry",	DoMUONEff());
    PrintOne("V0  NEEDS MCtruth ",      DoV0());
  }
};
VirtualQACfg* qaCfg = 0;

//====================================================================
/** 
 * Load the needed libraries 
 * 
 */
void LoadLibraries()
{
  Bool_t is10h = grp->period.EqualTo("LHC10h",TString::kIgnoreCase);
  gSystem->SetIncludePath("-I. " 
			"-I$ROOTSYS/include " 
			"-I$ALICE_ROOT/include " 
			"-I$ALICE_ROOT " 
			"-I$ALICE_ROOT/ITS " 
			"-I$ALICE_ROOT/TRD " 
			"-I$ALICE_PHYSICS/PWGPP " 
			"-I$ALICE_PHYSICS/PWGPP/TRD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libOADB");
  gSystem->Load("libESDfilter");
  gSystem->Load("libCORRFW");
  gSystem->Load("libTender");
  gSystem->Load("libPWGPP");
  gSystem->Load("libAliHLTTrigger");

  if ((qaCfg->DoEMCAL() && detCfg->UseEMCAL()) || 
      ((qaCfg->DoPHOS() || qaCfg->DoPHOSTrig())  && detCfg->UsePHOS()) || 
      (qaCfg->DoCALO() && !is10h)) {
    gSystem->Load("libEMCALUtils");
    gSystem->Load("libPHOSUtils");
    gSystem->Load("libPWGCaloTrackCorrBase");
    gSystem->Load("libPWGGACaloTrackCorrelations");
    gSystem->Load("libPWGGACaloTasks");
    gSystem->Load("libPWGGAPHOSTasks");
    gSystem->Load("libPWGTools");
    gSystem->Load("libPWGEMCAL");
    gSystem->Load("libPWGGAEMCALTasks");
  }  
  if((qaCfg->DoMUON() || qaCfg->DoMUONTrig()) && detCfg->UseMUON()) {
    gSystem->Load("libPWGmuon");
    gSystem->Load("libPWGPPMUONlite");
    gSystem->Load("libPWGmuondep");
  }
  if (qaCfg->DoFWD() && detCfg->UseFMD()) {
    gSystem->Load("libPWGLFforward2");
  }      
}
 
//====================================================================
/** 
 * Add the analysis tasks 
 * 
 * @param cdb_location 
 */
void AddAnalysisTasks(const char *cdb_location)
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  mgr->SetCommonFileName("QAresults.root");

  Bool_t is10h = grp->period.EqualTo("LHC10h",TString::kIgnoreCase);
  // --- Some short-hands --------------------------------------------
  TString ali   = "$ALICE_PHYSICS";
  TString ana   = "$ALICE_ROOT/ANALYSIS";
  TString oadb  = ali + "/OADB";
  TString pwghf = ali + "/PWGHF";
  TString pwglf = ali + "/PWGLF";
  TString pwgje = ali + "/PWGJE";
  TString pwgdq = ali + "/PWGDQ";
  TString pwgpp = ali + "/PWGPP";
  TString pwgga = ali + "/PWGGA";
  
  // --- Statistics task ---------------------------------------------
  mgr->AddStatisticsTask(kTriggerMask);

  // --- CDB connection ----------------------------------------------
  if (qaCfg->DoCDBconnect()) {
    gROOT->LoadMacro(pwgpp+"/PilotTrain/AddTaskCDBconnect.C");
    AliTaskCDBconnect *taskCDB = AddTaskCDBconnect(cdb_location, grp->run);
    if (!taskCDB) return;
  }    
  
  // --- Event Statistics (Jan Fiete) --------------------------------
  if (qaCfg->DoEventStat()) {
    gROOT->LoadMacro(oadb+"/macros/AddTaskPhysicsSelection.C");
    AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection(kTRUE/*MC*/);
    // Hack by Alexander for LHC10h
    // gROOT->LoadMacro("LHC10hPS.C");
    // AliOADBPhysicsSelection* ops = LHC10hPS(grp->period, grp->run);
    // if (ops) 
    //   physSelTask->GetPhysicsSelection()->SetCustomOADBObjects(ops,0);
  }
  // --- PIDResponse(JENS) -------------------------------------------
  if (qaCfg->DoPIDResponse() && !is10h) {
    gROOT->LoadMacro(ana+"/macros/AddTaskPIDResponse.C"); 
    AliAnalysisTaskPIDResponse *PIDResponse = AddTaskPIDResponse(kTRUE);
    PIDResponse->SelectCollisionCandidates(kTriggerMask);
  }    
  // --- Centrality (A. Toia) ----------------------------------------
  if (qaCfg->DoCentrality()) {
    gROOT->LoadMacro(oadb+"/macros/AddTaskCentrality.C");
    AliCentralitySelectionTask *taskCentrality = AddTaskCentrality();
    taskCentrality->SetMCInput();        
  }   
  
  // --- Vertexing (A. Dainese) --------------------------------------
  if (qaCfg->DoVertex()) {
    gROOT->LoadMacro(pwgpp+"/macros/AddTaskVertexESD.C");
    // Specific setting for MC
    AliAnalysisTaskVertexESD* taskvertexesd =  
      AddTaskVertexESD(kTRUE, kTriggerMask);
    taskvertexesd->SelectCollisionCandidates(kTriggerMask);
  }  

  // --- TPC QA (E. Sicking) -----------------------------------------
  if (qaCfg->DoQAsym()) {
  // offline trigger in AddTask
    gROOT->LoadMacro(pwgpp+"/PilotTrain/AddTaskQAsym.C");
    AliAnalysisTaskSE * taskqasim = AddTaskQAsym(0, 
						 kTriggerMask, 
						 kTriggerHM, 
						 kTriggerEMC, 
						 kTriggerMuonBarell);
  }  
  // --- VZERO QA  (C. Cheshkov) -------------------------------------
  if (qaCfg->DoVZERO() && detCfg->UseVZERO()) {
    gROOT->LoadMacro(pwgpp+"/PilotTrain/AddTaskVZEROQA.C");
    AliAnalysisTaskSE * taskv0qa = AddTaskVZEROQA(0);
  }
  if (qaCfg->DoVZEROPbPb() && detCfg->UseVZERO() && grp->IsAA()) {
    gROOT->LoadMacro(pwgpp+"/VZERO/AddTaskVZEROPbPb.C");
    AliAnaVZEROPbPb* taskV0PbPb = 
      (AliAnaVZEROPbPb*)AddTaskVZEROPbPb(Int_t(grp->run));
  }
  // --- TPC (Jacek Otwinowski & Michael Knichel) --------------------
  //
  //
  // - Optionally MC information can be used by setting the 1st
  //   argument to true
  // - Optionally friends information can be switched off by setting
  //   the 2st argument to false
  // - Optionally highMult axis can be used by setting the 3st
  //   argument to true (for PbPb)
  if (qaCfg->DoTPC() && detCfg->UseTPC()) {
    gROOT->LoadMacro(pwgpp+"/TPC/macros/AddTaskPerformanceTPCdEdxQA.C");
    AliPerformanceTask *tpcQA = 0;
    if (grp->IsAA()) {
       // High multiplicity Pb-Pb
       tpcQA = AddTaskPerformanceTPCdEdxQA(kTRUE, kTRUE, kTRUE);
    } else {
      // Low multiplicity (pp)
       tpcQA = AddTaskPerformanceTPCdEdxQA(kTRUE, kTRUE, kFALSE);
    }
    tpcQA->SelectCollisionCandidates(kTriggerMask);
    AliPerformanceRes::SetMergeEntriesCut(5000000); 
  }  

  // --- HLT (Alberica Toia) -----------------------------------------
  if (qaCfg->DoHLT() && detCfg->UseTPC()) {
    gROOT->LoadMacro(pwgpp+"/TPC/macros/AddTaskPerformanceTPCdEdxQA.C");
    AliPerformanceTask *hltQA = AddTaskPerformanceTPCdEdxQA(kTRUE, kTRUE, 
							    kFALSE,0,kTRUE);
    hltQA->SelectCollisionCandidates(kTriggerMask);
  }  
  // --- SPD (A. Mastroserio) ----------------------------------------
  if (qaCfg->DoSPD() && detCfg->UseITS()) {
    gROOT->LoadMacro(pwgpp+"/PilotTrain/AddTaskSPDQA.C");
    AliAnalysisTaskSPD* taskspdqa = (AliAnalysisTaskSPD*)AddTaskSPDQA();
    // Request from Annalisa
    if (grp->IsAA()) taskspdqa->SetHeavyIonMode();
    taskspdqa->SelectCollisionCandidates(kTriggerMask);
    taskspdqa->SetOCDBInfo(grp->run, "raw://");
  }  
  // --- SDD (F. Prino) ----------------------------------------------
  if (qaCfg->DoSDD() && detCfg->UseITS()) {
    gROOT->LoadMacro(pwgpp+"/PilotTrain/AddSDDPoints.C");
    AliAnalysisTaskSE* tasksdd = AddSDDPoints();
    tasksdd->SelectCollisionCandidates(kTriggerMask);
  }
  // --- SSD dEdx (Marek Chojnacki) ----------------------------------
  if (qaCfg->DoSSDdEdx() && detCfg->UseITS()) {
    gROOT->LoadMacro(pwgpp+"/PilotTrain/AddTaskdEdxSSDQA.C");
    AliAnalysisTaskSE* taskssddedx = AddTaskdEdxSSDQA();
    taskssddedx->SelectCollisionCandidates(kTriggerMask);
  }

  // --- ITS ---------------------------------------------------------
  if (qaCfg->DoITS() && detCfg->UseITS()) {
    // hardcoded non-zero trigger mask
    gROOT->LoadMacro(pwgpp+"/macros/AddTaskPerformanceITS.C");
    AliAnalysisTaskITSTrackingCheck *itsQA = 0;
    AliAnalysisTaskITSTrackingCheck *itsQACent0010 = 0;
    AliAnalysisTaskITSTrackingCheck *itsQACent3050 = 0;
    AliAnalysisTaskITSTrackingCheck *itsQACent6080 = 0;
    if(grp->IsPP()) {
      itsQA = AddTaskPerformanceITS(kTRUE);
    } else {
      itsQA = AddTaskPerformanceITS(kTRUE);
      itsQACent0010 = AddTaskPerformanceITS(kTRUE,kFALSE,kFALSE,3500,10000);
      itsQACent3050 = AddTaskPerformanceITS(kTRUE,kFALSE,kFALSE,590,1570);
      itsQACent6080 = AddTaskPerformanceITS(kTRUE,kFALSE,kFALSE,70,310);
    }
  }
  // --- ITS saTracks, align (F.Prino) -------------------------------
  if (qaCfg->DoITSsaTracks() && detCfg->UseITS()) {
    // offline trigger in AddTask
    gROOT->LoadMacro(pwgpp+"/macros/AddTaskITSsaTracks.C");
    AliAnalysisTaskITSsaTracks *itssaTracks = AddTaskITSsaTracks(kTRUE,kFALSE);
    itssaTracks->SelectCollisionCandidates(kTriggerMask);
  }   
  if (qaCfg->DoITSalign() && detCfg->UseITS()) {
    // no offline trigger selection
     gROOT->LoadMacro(pwgpp+"/macros/AddTaskITSAlign.C");
     AliAnalysisTaskITSAlignQA *itsAlign = AddTaskITSAlign(0,2011);
  }   

  // --- TRD (Alex Bercuci, M. Fasel) --------------------------------
  if(qaCfg->DoTRD() && detCfg->UseTRD()) {
    // no offline trigger selection
    gROOT->LoadMacro(pwgpp+"/macros/AddTrainPerformanceTRD.C");
    // steer individual TRD tasks
    Bool_t 
      doCheckESD(kTRUE),  // AliTRDcheckESD
      doCheckDET(kTRUE),  // AliTRDcheckDET
      doEffic(kTRUE),     // AliTRDefficiency
      doResolution(kTRUE),// AliTRDresolution
      doCheckPID(kTRUE),  // AliTRDcheckPID
      doV0Monitor(kFALSE);// AliTRDv0Monitor
    AddTrainPerformanceTRD(Translate(doCheckESD, doCheckDET, doEffic, 
				     doResolution, doCheckPID, doV0Monitor));
  }

  // --- ZDC (Chiara Oppedisano) -------------------------------------
  if(qaCfg->DoZDC() && detCfg->UseZDC()) {
    // hardcoded kMB trigger mask
     gROOT->LoadMacro(pwgpp+"/ZDC/AddTaskZDCQA.C");
     AliAnalysisTaskSE *taskZDC = AddTaskZDCQA();
     taskZDC->SelectCollisionCandidates(kTriggerMask);
  }   

  // --- Calorimetry (Gustavo Conesa) --------------------------------
  if(qaCfg->DoCALO() && !is10h) {
    gROOT->LoadMacro(pwgga+
		     "/CaloTrackCorrelations/macros/QA/AddTaskCalorimeterQA.C");
    AliAnalysisTaskCaloTrackCorrelation *taskCaloQA = 
      AddTaskCalorimeterQA("default");
    taskCaloQA->SetDebugLevel(0);
    // offline mask set in AddTask to kMB
    taskCaloQA->SelectCollisionCandidates(kTriggerMask);
    // Add a new calo task with EMC1 trigger only
    if (!is10h) {
      taskCaloQA = AddTaskCalorimeterQA("trigEMC");
      taskCaloQA->SetDebugLevel(0);
      taskCaloQA->SelectCollisionCandidates(kTriggerEMC);
    }
  }

  // --- Muon Trigger ------------------------------------------------
  if(qaCfg->DoMUONTrig() && detCfg->UseMUON()) {
    // no offline trigger selection
    gROOT->LoadMacro(pwgpp+"/macros/AddTaskMTRchamberEfficiency.C");
    AliAnalysisTaskTrigChEff *taskMuonTrig = AddTaskMTRchamberEfficiency();
  }

  // --- Muon Efficiency (not used) ----------------------------------
  if(qaCfg->DoMUONEff() && detCfg->UseMUON()) {
      gROOT->LoadMacro(ali+"/PWG3/muondep/AddTaskMUONTrackingEfficiency.C");
      AliAnalysisTaskMuonTrackingEff *taskMuonTrackEff = 
	AddTaskMUONTrackingEfficiency();
  }
  
  // --- V0-Decay Reconstruction (Ana Marin) (not used) --------------
  if (qaCfg->DoV0()) {
    gROOT->LoadMacro(pwgpp+"/macros/AddTaskV0QA.C");
    AliAnalysisTaskV0QA *taskv0QA = AddTaskV0QA(kTRUE);
  }
  
  // -- Impact parameter resolution ----------------------------------
  // (xianbao.yuan@pd.infn.it, andrea.dainese@pd.infn.it)
  if (qaCfg->DoImpParRes()) {
    gROOT->LoadMacro(pwgpp+"/macros/AddTaskImpParRes.C");
    AliAnalysisTaskSE* taskimpparres=0;
    // Specific setting for MC
    if(grp->IsPP()) {
      taskimpparres= AddTaskImpParRes(kTRUE);
    } else {
      taskimpparres= AddTaskImpParRes(kTRUE,-1,kTRUE,kFALSE);
    }
    taskimpparres->SelectCollisionCandidates(kTriggerMask);
  }  

  // --- MUON QA (Philippe Pillot) -----------------------------------
  if (qaCfg->DoMUON() && detCfg->UseMUON()) {
    // trigger analysis internal
    gROOT->LoadMacro(pwgpp+"/PilotTrain/AddTaskMuonQA.C");
    AliAnalysisTaskSE* taskmuonqa= AddTaskMuonQA();
  }  

  // --- TOF (Francesca Bellini) -------------------------------------
  if (qaCfg->DoTOF() && detCfg->UseTOF()) {
    gROOT->LoadMacro(pwgpp+"/TOF/AddTaskTOFQA.C");
    AliAnalysisTaskTOFqa *tofQA = AddTaskTOFQA(kFALSE);
    tofQA->SelectCollisionCandidates(kTriggerMask);
  } 

  // --- PIDqa(JENS) -------------------------------------------------
  if (qaCfg->DoPIDqa() && !is10h) {
    gROOT->LoadMacro(ana+"/macros/AddTaskPIDqa.C");
    AliAnalysisTaskPIDqa *PIDQA = AddTaskPIDqa();
    PIDQA->SelectCollisionCandidates(kTriggerMask);
  }  
 
  // --- HMPID QA (Giacomo Volpe) ------------------------------------
  //
  if (qaCfg->DoHMPID() && detCfg->UseHMPID()) {
    gROOT->LoadMacro(pwgpp+"/HMPID/AddTaskHmpidQA.C");
    AliAnalysisTaskSE* taskhmpidqa= AddTaskHmpidQA(kTRUE);
    // offline mask set in AddTask to kMB
    taskhmpidqa->SelectCollisionCandidates(kTriggerMask);
  }      

  // --- T0 QA (Alla Mayevskaya) -------------------------------------
  if (qaCfg->DoT0() && detCfg->UseT0()) {
    // no offline trigger selection
    gROOT->LoadMacro(pwgpp+"/T0/AddTaskT0QA.C");
    AliT0AnalysisTaskQA* taskt0qa= AddTaskT0QA();
    taskt0qa->SelectCollisionCandidates(kTriggerMask);
  }      

  // ---- FMD QA (Christian Holm Christiansen) -----------------------
  if (qaCfg->DoFWD() && detCfg->UseFMD()) {
    gROOT->LoadMacro(pwglf+"/FORWARD/analysis2/AddTaskForwardQA.C");
    // Parameters: usemc, usecentrality
    // AliAnalysisTaskSE *forwardQA = (AliAnalysisTaskSE *)
    AddTaskForwardQA(kTRUE, qaCfg->DoCentrality());
    // HACK: to read corrections from current directory
    const char* hack="AliForwardCorrectionManager::Instance().SetPrefix(\".\")";
    gROOT->ProcessLine(hack);
    const char* hack2="AliForwardCorrectionManager::Instance().Print(\"R\")";
    gROOT->ProcessLine(hack2);
    // No offline trigger config. needed (see #84077)
  }
  
  // --- PHOS QA (Boris Polishchuk) ----------------------------------
  if (qaCfg->DoPHOS()&& detCfg->UsePHOS()) {
    gROOT->LoadMacro(pwgga+"/PHOSTasks/CaloCellQA/macros/AddTaskCaloCellsQA.C");
    AliAnalysisTaskCaloCellsQA *taskPHOSCellQA1 = 
      AddTaskCaloCellsQA(4, 1, NULL,"PHOSCellsQA_AnyInt"); 
    taskPHOSCellQA1->SelectCollisionCandidates(kTriggerMask);
    taskPHOSCellQA1->GetCaloCellsQA()->SetClusterEnergyCuts(0.3,0.3,1.0);

    AliAnalysisTaskCaloCellsQA *taskPHOSCellQA2 = 
      AddTaskCaloCellsQA(4, 1, NULL,"PHOSCellsQA_PHI7"); 
    taskPHOSCellQA2->SelectCollisionCandidates(AliVEvent::kPHI7);
    taskPHOSCellQA2->GetCaloCellsQA()->SetClusterEnergyCuts(0.3,0.3,1.0);

    // Pi0 QA fo PbPb
    if (grp->IsAA()) {
      gROOT->LoadMacro(pwgga+"/PHOSTasks/PHOS_PbPbQA/macros/AddTaskPHOSPbPb.C");
      AliAnalysisTaskPHOSPbPbQA* phosPbPb = AddTaskPHOSPbPbQA(0);
    }
  }
  if (qaCfg->DoPHOSTrig() && detCfg->UsePHOS()) {
    gROOT->LoadMacro(pwgga+
		     "/PHOSTasks/PHOS_TriggerQA/macros/AddTaskPHOSTriggerQA.C");
    AliAnalysisTaskPHOSTriggerQA *taskPHOSTrig = AddTaskPHOSTriggerQA(0,0);
  }   

  // --- EMCAL QA (Gustavo Conesa) -----------------------------------
  if (qaCfg->DoEMCAL() && detCfg->UseEMCAL()) {
     gROOT->LoadMacro(pwgga+"/EMCALTasks/macros/AddTaskEMCALTriggerQA.C");
     AliAnalysisTaskEMCALTriggerQA *emctrig = AddTaskEMCALTriggerQA();
  }   

  // --- FLOW and BF QA (C.Perez && A.Rodriguez) ---------------------
  if (qaCfg->DoFBFqa()) {
    gROOT->LoadMacro(pwgpp+"/macros/AddTaskFBFqa.C");
    AliAnalysisTaskSE *qaFBFMB = (AliAnalysisTaskSE*)AddTaskFBFqa("qaFBFmb",
								  kFALSE);
    qaFBFMB->SelectCollisionCandidates(AliVEvent::kMB);
    AliAnalysisTaskSE *qaFBFSC = (AliAnalysisTaskSE*)AddTaskFBFqa("qaFBFsc",
								  kFALSE);
    qaFBFSC->SelectCollisionCandidates(AliVEvent::kSemiCentral);
    AliAnalysisTaskSE *qaFBFCE = (AliAnalysisTaskSE*)AddTaskFBFqa("qaFBFce",
								  kFALSE);
    qaFBFCE->SelectCollisionCandidates(AliVEvent::kCentral);
  }
}
/** 
 * Helper function to make @c outputs_valid file 
 * 
 */
void ValidateOutput()
{
  std::ofstream out;
  out.open("outputs_valid", ios::out);
  out.close();    
}  

//====================================================================
/** 
 * Run QA merging 
 * 
 * @param dir      directory 
 * @param stage    stage 
 */
void QAMerge(const char *dir, Int_t stage)
{
// Merging method
  TStopwatch  timer;     timer.Start();
  TString     outputDir     = dir;
  TObjArray   outputFiles;
  outputFiles.Add(new TObjString("QAresults.root"));
  outputFiles.Add(new TObjString("EventStat_temp.root"));

  TString     mergeExcludes = "";
  TIter       iter(&outputFiles);
  TObjString* str           = 0;
  Bool_t      merged        = kTRUE;
  while((str = static_cast<TObjString*>(iter()))) {
    TString& outputFile = str->GetString();
    // Skip already merged outputs
    if (!gSystem->AccessPathName(outputFile)) {
      ::Warning("Merge","Output file <%s> found. Not merging again.",
		outputFile.Data());
      continue;
    }
    if (mergeExcludes.Contains(outputFile.Data())) continue;
    merged = AliAnalysisAlien::MergeOutput(outputFile, 
					   outputDir, 
					   10, 
					   stage);
    if (!merged) {
      ::Error("Merge", "Cannot merge %s\n", outputFile.Data());
      continue;
    }
  }
  TString infolog = "fileinfo.log";
  AliAnalysisAlien::MergeInfo(infolog, dir); 

  if (!outputDir.Contains("Stage")) {
    ValidateOutput();
    timer.Print();
    return;
  }
  // --- Set up to run terminate -------------------------------------
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  mgr->SetSkipTerminate(kFALSE);
  if (!mgr->InitAnalysis()) return;

  mgr->PrintStatus();
  mgr->StartAnalysis("gridterminate", (TTree*)0);
  ValidateOutput();
  timer.Print();
}

//====================================================================
/** 
 * Run QA trains 
 * 
 * @param run         Run number 
 * @param xmlfile     Collection file 
 * @param stage       Stage 
 * @param cdb         CDB location 
 */
void QA(UInt_t      run, 
	const char* xmlfile   = "wn.xml",
	Int_t       stage     = 0, /*0 = QA train, 1...n - merging stage*/
	const char* cdb       = "raw://")
{
  // -----------------------------------------------------------------
  // 
  // Get GRP parameters.  Defines global "grp" as a pointer to GRPData
  //
  gROOT->Macro(Form("GRP.C(%d)", run));
  gROOT->Macro("BaseConfig.C");
  gROOT->Macro("QAConfig.C");
  gROOT->Macro("DetConfig.C");
  qaCfg->Print();
  Int_t   debug_level = qaCfg->DebugLevel();   // Debugging
  TString cdbString(cdb);
  if (cdbString.Contains("raw://")) {
    TGrid::Connect("alien://");
    if (!gGrid || !gGrid->IsConnected()) {
      ::Error("QAtrain", "No grid connection");
      return;
    }  
  }  

  // --- Some settings -----------------------------------------------
  // Set temporary merging directory to current one
  gSystem->Setenv("TMPDIR", gSystem->pwd());
  // Set temporary compilation directory to current one
  gSystem->SetBuildDir(gSystem->pwd(), kTRUE);
  // Load common libraries and set include path
  LoadLibraries();
  printf("Include path: %s\n", gSystem->GetIncludePath());

  // === Make the analysis manager and connect event handlers ========
  // 
  // --- Analysis manager and load libraries -------------------------
  AliAnalysisManager *mgr  = new AliAnalysisManager("QA", "Production train");
  mgr->SetRunFromPath(grp->run);

  // --- Create ESD input handler ------------------------------------
  AliESDInputHandlerRP *esdHandler = new AliESDInputHandlerRP();
  esdHandler->SetReadFriends(kTRUE);
  esdHandler->SetActiveBranches("ESDfriend");
  mgr->SetInputEventHandler(esdHandler);
  
  // --- Monte Carlo handler -----------------------------------------
  if (true) {
    AliMCEventHandler* mcHandler = new AliMCEventHandler();
    mgr->SetMCtruthEventHandler(mcHandler);
    mcHandler->SetPreReadMode(1);
    mcHandler->SetReadTR(true);
  }

  // === Set up tasks ================================================
  //
  // --- Create tasks ------------------------------------------------
  AddAnalysisTasks(cdb);

  // --- Debugging if needed -----------------------------------------
  if (debug_level > 0) mgr->SetDebugLevel(debug_level);

  // --- If merging, do so here and exit -----------------------------
  if (stage>0) {
    QAMerge(xmlfile, stage);
    return;
  }   

  // === Run the analysis ============================================
  //
  // --- Make our chain ----------------------------------------------
  TChain *chain = new TChain("esdTree");
  chain->Add("AliESDs.root");

  // --- Run the thing -----------------------------------------------
  TStopwatch timer;
  timer.Start();
  if (!mgr->InitAnalysis()) return;

  mgr->PrintStatus(); 
  mgr->SetSkipTerminate(kTRUE);
  mgr->SetNSysInfo(1);
  mgr->StartAnalysis("local", chain);
  timer.Print();
}

// 
// EOF
// 

