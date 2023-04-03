/**************************************************************************/
// Dukhishyam Mallick
// Macro to configure the KStarPlusMinus analysis task
 // It calls all configs desired by the user, by means
 // of the boolean switches defined in the first lines.
 // ---
 // Inputs:
 //  1) flag to know if running on MC or data
 //  2) collision system, whether pp, pPb or PbPb
 // --
 // Returns:
 //  kTRUE  --> initialization successful
 //  kFALSE --> initialization failed (some config gave errors)
 //


/****************************************************************************/

//enum ERsnCollType_t { kPP=0,      kPPb,      kPbPb};
#ifdef __CLING__
R__ADD_INCLUDE_PATH($ALICE_PHYSICS)
#include <PWGLF/RESONANCES/macros/mini/ConfigKsKs.C>
//#include "ConfigKsKs.C"
#endif

enum pairYCutSet { kPairDefault,    // USED ONLY FOR pA
    kNegative,       // USED ONLY FOR pA
    kCentral         // USED ONLY FOR pA
};

enum eventCutSet { kEvtDefault=0,
    kNoPileUpCut, //=1
    kDefaultVtx12,//=2
    kDefaultVtx8, //=3
    kDefaultVtx5, //=4
    kMCEvtDefault, //=5
    kSpecial1, //=6
    kSpecial2, //=7
    kNoEvtSel, //=8
    kSpecial3, //=9
    kSpecial4, //=10
    kSpecial5 //=11
};

enum eventMixConfig { kDisabled = -1,
    kMixDefault,     //=0 //10 events, Dvz = 1cm, DC = 10
    k5Evts,          //=1 //5 events, Dvz = 1cm, DC = 10
    k5Cent,          //=2 //10 events, Dvz = 1cm, DC = 5
};
AliRsnMiniAnalysisTask *AddTaskKsKs
(
 TString     outNameSuffix = "KStarPlusMinus_V0Mass_Pt",
 Bool_t      isMC =0,
 Int_t       isAOD = 1,
 Bool_t      isPP=1,
 AliPIDResponse::EBeamType collSys = AliPIDResponse::kPP,
 Int_t       evtCutSetID = 0,
 Int_t       pairCutSetID = 0,
 Bool_t      enaMultSel = kTRUE,    //enable multiplicity axis
 Int_t       aodFilterBit = 0,
 Bool_t      enableMonitor=kTRUE,
 TString     monitorOpt="pp",
 Float_t     piPIDCut = 3.0,
 Int_t       customQualityCutsID=1,
 AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutPiCandidate = AliRsnCutSetDaughterParticle::kTPCTOFpidphikstarpPb2016,
 Int_t       polaxis=0,
 Float_t     pi_k0s_PIDCut = 5.0,
 Float_t     massTol = 0.03,
 Float_t     massTolVeto = 0.0043,
 Int_t       tol_switch = 1,
 Double_t    tol_sigma = 6,
 Float_t     pLife = 20,
 Float_t     radiuslow = 0.5,
 Bool_t      Switch = kTRUE,
 Float_t     k0sDCA = 1.0,
 Float_t     k0sCosPoinAn = 0.97,
 Float_t     k0sDaughDCA = 1.0,
 Float_t     DCAxy = 0.06,
 Bool_t      enableSys = kFALSE,
 Int_t       Sys= 0,
 UInt_t      triggerMask=AliVEvent::kINT7, 
 Float_t     masslow = 0.5,//inv mass axis low edge
 Float_t     massup = 2.5, //inv mass axis upper edge
 Int_t       nbins = 200, //inv mass axis n bin
 Float_t     ptlow = 0.0, //pT axis low edge
 Float_t     ptup = 20.0,  //pT axis upper edge
 Int_t       nbinspt = 200, //pT axis n bins
 Float_t     coslow = 0.0,  //costheta axis low edge
 Float_t     cosup = 1.0,  //costheta axis upper edge
 Int_t       nbinscos = 10  // ncosbins
 )
{

  Float_t     cutV = 10.0;
  Bool_t      isGT = 0;
  Int_t       NTPCcluster = 70;
  Float_t     crossedRows = 70;
  Float_t     rowsbycluster = 0.8;
  Float_t      v0rapidity= 0.5;
    
  //-------------------------------------------
    // event cuts
    //-------------------------------------------
   Double_t vtxZcut            = 10.0; //cm, default cut on vtx z
  Bool_t rejectPileUp         = kTRUE; //rejects pileup from SPD
  Bool_t useMVPileUpSelection = kFALSE; //alternative pile-up rejection, default is SPD
  Int_t  MinPlpContribSPD     = 5; //default value if used
  Int_t  MinPlpContribMV      = 5; //default value if used
  // Bool_t selectDPMJETevtNSDpA = kFALSE; //cut to select true NSD events in DPMJET

  if (evtCutSetID==eventCutSet::kDefaultVtx12) vtxZcut = 12.0; //cm
  if (evtCutSetID==eventCutSet::kDefaultVtx5) vtxZcut = 5.0; //cm
  if (evtCutSetID==eventCutSet::kNoPileUpCut) rejectPileUp = kFALSE;

  /*
  if (evtCutSetID>=eventCutSet::kMCEvt) {
    vtxZcut = 1.0e3; //cm
  }
  
  if (evtCutSetID==eventCutSet::kMCEvtDefault) {
    vtxZcut = 10.0; //cm
  }
  */
  //-------------------------------------------
  //pair rapidity cut
  //-------------------------------------------
  Double_t minYlab = -0.5;
  Double_t maxYlab =  0.5;
  /*
  if (pairCutSetID==pairYCutSet::kCentralTight) { //|y_cm|<0.3
    minYlab = -0.3;    maxYlab = 0.3;
  }

  if (pairCutSetID==pairYCutSet::kpADefault) { //-0.5<y_cm<0.0
    minYlab = -0.465;    maxYlab = 0.035;
  }
  
  if (pairCutSetID==pairYCutSet::kpACentral) { //|y_cm|<0.3
    minYlab = -0.765;    maxYlab = -0.165;
  }
  */
  //-------------------------------------------
  //mixing settings
  //-------------------------------------------
  Int_t       nmix = 5;
  Float_t     maxDiffVzMix = 1.0;
  Float_t     maxDiffMultMix = 5.0;
  if (collSys==AliPIDResponse::kPBPB) maxDiffMultMix = 10.0;
  
  // -- INITIALIZATION ----------------------------------------------------------------------------
  //
  // retrieve analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskF0", "No analysis manager to connect to.");
    return NULL;
  } 

  // create the task and configure 
  TString taskName = Form("RsnTaskF0");
  AliRsnMiniAnalysisTask *task = new AliRsnMiniAnalysisTask(taskName.Data(), isMC);

  //trigger 
  if (!isAOD) task->UseESDTriggerMask(triggerMask); //ESD
  else task->SelectCollisionCandidates(triggerMask); //AOD
  
  //-----------------------------------------------------------------------------------------------
  // -- MULTIPLICITY/CENTRALITY -------------------------------------------------------------------
  //-----------------------------------------------------------------------------------------------
  if (collSys==AliPIDResponse::kPP) task->UseMultiplicity("AliMultSelection_V0M");
  if (collSys==AliPIDResponse::kPPB) task->UseMultiplicity("AliMultSelection_V0A");
  if (collSys==AliPIDResponse::kPBPB) task->UseMultiplicity("AliMultSelection_V0M");
  
  //-----------------------------------------------------------------------------------------------
  // -- EVENT MIXING CONFIG -----------------------------------------------------------------------
  //-----------------------------------------------------------------------------------------------
  task->UseContinuousMix();
  task->SetNMix(nmix);
  task->SetMaxDiffVz(maxDiffVzMix);
  task->SetMaxDiffMult(maxDiffMultMix);
  //::Info("AddTaskF0", Form("Event mixing configuration: \n events to mix = %i \n max diff. vtxZ = cm %5.3f \n max diff multi = %5.3f", nmix, maxDiffVzMix, maxDiffMultMix));
   
  //-----------------------------------------------------------------------------------------------
  // -- EVENT SELECTION ---------------------------------------------------------------------------
  //-----------------------------------------------------------------------------------------------
  AliRsnCutPrimaryVertex *cutVertex = new AliRsnCutPrimaryVertex("cutVertex", 10.0, 0, kFALSE);
  if ((collSys==AliPIDResponse::kPP) && (!isMC)) cutVertex->SetCheckPileUp(rejectPileUp);   // set the check for pileup
  cutVertex->SetCheckZResolutionSPD(); 
  Printf("AddTaskF0 - CheckZResolutionSPD:              ON");
  cutVertex->SetCheckDispersionSPD(); 
  Printf("AddTaskF0 - CheckDispersionSPD:               ON");
  cutVertex->SetCheckZDifferenceSPDTrack(); 
  Printf("AddTaskF0 - CheckZDifferenceSPDTrack:         ON");

  //set check for pileup in 2013
  AliRsnCutEventUtils *cutEventUtils = new AliRsnCutEventUtils("cutEventUtils", kFALSE, rejectPileUp);
  cutEventUtils->SetCheckIncompleteDAQ(kTRUE);
  cutEventUtils->SetCheckSPDClusterVsTrackletBG();
  Printf("AddTaskF0 - CheckIncompleteDAQ:                  ON");
  Printf("AddTaskF0 - SetCheckSPDClusterVsTrackletBG:      ON");
  
  if (useMVPileUpSelection){
    cutEventUtils->SetUseMVPlpSelection(useMVPileUpSelection);
    cutEventUtils->SetMinPlpContribMV(MinPlpContribMV);
    cutEventUtils->SetMinPlpContribSPD(MinPlpContribSPD);
    //::Info("AddTaskF0", Form("Multiple-vtx Pile-up rejection:      ON \nSettings: MinPlpContribMV = %i, MinPlpContribSPD = %i", MinPlpContribMV, MinPlpContribSPD));
  } else {
    cutEventUtils->SetMinPlpContribSPD(MinPlpContribSPD);
    Printf("AddTaskF0 - SPD Pile-up rejection:     ON \nSettings: MinPlpContribSPD = %i", MinPlpContribSPD);
  }
  Printf("AddTaskF0 - Pile-up rejection mode:     %i", rejectPileUp);   
  
  AliRsnCutSet *eventCuts = new AliRsnCutSet("eventCuts", AliRsnTarget::kEvent);
  eventCuts->AddCut(cutEventUtils);
  eventCuts->AddCut(cutVertex);
  eventCuts->SetCutScheme(Form("%s&%s", cutEventUtils->GetName(), cutVertex->GetName()));
  task->SetEventCuts(eventCuts);
   
  //connect task
  mgr->AddTask(task);
    
    // -- EVENT-ONLY COMPUTATIONS -------------------------------------------------------------------
    //vertex
    Int_t vtxID=task->CreateValue(AliRsnMiniValue::kVz,kFALSE);
    AliRsnMiniOutput* outVtx=task->CreateOutput("eventVtx","HIST","EVENT");
    outVtx->AddAxis(vtxID,240,-12.0,12.0);
    
    //multiplicity
    Int_t multID=task->CreateValue(AliRsnMiniValue::kMult,kFALSE);
    AliRsnMiniOutput* outMult=task->CreateOutput("eventMult","HIST","EVENT");
    if(isPP)
      outMult->AddAxis(multID,400,0.5,400.5);
    else outMult->AddAxis(multID,100,0.,100.);
    
    TH2F* hvz=new TH2F("hVzVsCent","",100,0.,100., 240,-12.0,12.0);
    task->SetEventQAHist("vz",hvz);//plugs this histogram into the fHAEventVz data member
    
    TH2F* hmc=new TH2F("MultiVsCent","", 100,0.,100., 400,0.5,400.5);
    hmc->GetYaxis()->SetTitle("QUALITY");
    task->SetEventQAHist("multicent",hmc);//plugs this histogram into the fHAEventMultiCent data member
    
    
    // -- PAIR CUTS (common to all resonances) ------------------------------------------------------
    
    AliRsnCutMiniPair* cutY=new AliRsnCutMiniPair("cutRapidity", AliRsnCutMiniPair::kRapidityRange);
    if (collSys==AliPIDResponse::kPP) cutY->SetRangeD(minYlab,maxYlab);
    if (collSys==AliPIDResponse::kPBPB) cutY->SetRangeD(minYlab,maxYlab);
    if (collSys==AliPIDResponse::kPPB) cutY->SetRangeD(-0.465,0.035);
    //cutY->SetRangeD(-0.465, 0.035);// 0 < y_cm < 0.5; y_cm = y_lab + 0.465
    
    AliRsnCutMiniPair* cutV0=new AliRsnCutMiniPair("cutV0", AliRsnCutMiniPair::kContainsV0Daughter);
    
    AliRsnCutSet* PairCutsSame=new AliRsnCutSet("PairCutsSame",AliRsnTarget::kMother);
    PairCutsSame->AddCut(cutY);
    PairCutsSame->AddCut(cutV0);
    PairCutsSame->SetCutScheme(TString::Format("%s&(!%s)",cutY->GetName(),cutV0->GetName()).Data());
    //note the use of the ! operator in this cut scheme
    
    AliRsnCutSet* PairCutsMix=new AliRsnCutSet("PairCutsMix",AliRsnTarget::kMother);
    PairCutsMix->AddCut(cutY);
    PairCutsMix->SetCutScheme(cutY->GetName());
    //
    // -- CONFIG ANALYSIS --------------------------------------------------------------------------
    if (isMC) {
        Printf("========================== MC analysis - PID cuts not used");
    } else
        Printf("========================== DATA analysis - PID cuts used");

    // #ifdef __CINT__
    //  gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/ConfigKsKs.C");
      // gROOT->LoadMacro("ConfigKsKs.C");
    // #endif 
          if (!ConfigKsKs(task, isPP, isMC,customQualityCutsID, aodFilterBit,enableMonitor,enaMultSel, monitorOpt.Data(),cutPiCandidate,polaxis, piPIDCut, pi_k0s_PIDCut, massTol, massTolVeto, tol_switch, tol_sigma, pLife, radiuslow, Switch, k0sDCA, k0sCosPoinAn, k0sDaughDCA, "", PairCutsSame, PairCutsMix, DCAxy, enableSys, Sys,masslow,massup,nbins,ptlow,ptup,nbinspt,coslow,cosup,nbinscos)) return 0x0;
    // -- CONTAINERS --------------------------------------------------------------------------------  //
   

    TString outputFileName = AliAnalysisManager::GetCommonFileName();
    //  outputFileName += ":Rsn";
    Printf("AddTaskKStarPlusMinus - Set OutputFileName : \n %s\n", outputFileName.Data() );
    AliAnalysisDataContainer *output = mgr->CreateContainer(Form("RsnOut_%s_%.1f_%.1f_%.1f_%.2f_%.3f_%.f_%.f_%.f_%.1f_%.2f_%.1f_%.3f_%.1f", outNameSuffix.Data(),piPIDCut,customQualityCutsID,pi_k0s_PIDCut,massTol,massTolVeto,pLife,radiuslow, k0sDCA,k0sCosPoinAn,k0sDaughDCA, DCAxy, Sys), TList::Class(), AliAnalysisManager::kOutputContainer, outputFileName);
    
    mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(task, 1, output);
    
    return task;
}
