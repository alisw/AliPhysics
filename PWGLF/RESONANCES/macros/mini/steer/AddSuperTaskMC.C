/***************************************************************************
              fbellini@cern.ch - last modified on 28/11/2013
//
//Lauches KStar analysis with rsn mini package
//Allows basic configuration of pile-up check and event cuts
//
****************************************************************************/
enum pairYCutSet { kPairDefault,
		   kpADefault,
		   kpANegative,
		   kpACentral,
		   kCentralTight
                 };

enum eventCutSet { kOld = -1, 
		   kEvtDefault,
		   kDefaultVtx8,
		   kDefaultVtx5,
		   kMCEvt,
		   kMCEvtDefault, //
		   kpAEvtDefault,		
		   kpANoPileUpCut,
		   kpAPileUpMV,
		   kpAPileUpSPD3,		      
		   kMCEvtpA2013,
		   kMCpAEvtDefault,
		   kMCEvtNSD,
		   kMCEvtNSDpA2013, 
		   kMCEvtNSDdefault
};

enum ERsnSpecie_t{kRho=0,
		  kKstar,
		  kAKstar,
		  kPhi,
		  kDelta0,
		  kDeltaPP,
		  kADeltaMM,
		  kLstar,
		  kALstar,
		  kF0,
		  kNrsn};

//-------------------------------------------------------------------
Bool_t ConfigSuperOutput(AliRsnMiniAnalysisTask *task, 
			 Bool_t                 isPP,
			 AliRsnCutSet           *cutsPair,
			 TString                partname="kStar",
			 Int_t                  pdgCode=313,
			 Float_t                mass = 0.89445,
			 Float_t                masslow = 0.7,
			 Float_t                massup = 1.1,
			 Int_t                  nbins = 400,
			 Char_t                 charge1 = '+',
			 Char_t                 charge2 = '-',
			 RSNPID                 d1 = AliRsnDaughter::kKaon,
			 RSNPID                 d2 = AliRsnDaughter::kPion,
			 Int_t                  icut1 = 0,
			 Int_t                  icut2 = 0);

AliRsnMiniAnalysisTask * AddSuperTaskMC(Bool_t      isPP = kFALSE,
					TString     outNameSuffix = "q2011",
					Int_t       evtCutSetID = 0,
					Int_t       pairCutSetID = 0,
					Int_t       aodFilterBit = 5,
					Bool_t      enableMonitor = kFALSE,
					TString     monitorOpt = "NoSIGN")
{  

  //-------------------------------------------
  // event cuts
  //-------------------------------------------
  UInt_t      triggerMask = AliVEvent::kMB;
  Bool_t      rmFirstEvtChunk = kFALSE; //needed for pA 2013
  Bool_t      rejectPileUp = kFALSE; //best if used, for pA 2013
  Int_t       MinPlpContribSPD = 5; //default value if used
  Bool_t      useMVPileUpSelection = kFALSE; //
  Int_t       MinPlpContribMV = 5; //default value if used
  Bool_t      useVtxCut2013pA = kFALSE; //default use recommended 2013 pA vtx selection
  Double_t    vtxZcut = 10.0; //cm, default cut on vtx z
  Bool_t      selectDPMJETevtNSDpA = kFALSE; //cut to select in DPMJET true NSD events
  
  if (evtCutSetID==eventCutSet::kOld) {
    triggerMask = AliVEvent::kAnyINT;
  }
    
  if (evtCutSetID==eventCutSet::kDefaultVtx8){
    vtxZcut = 8.0; //cm
  } 
  
  if (evtCutSetID==eventCutSet::kDefaultVtx5){
    vtxZcut = 5.0; //cm
  }

 if (evtCutSetID==eventCutSet::kMCEvt) {
    rmFirstEvtChunk = kFALSE;
    rejectPileUp = kFALSE;
    useVtxCut2013pA = kFALSE;
    vtxZcut = 1.0e3; //cm
    selectDPMJETevtNSDpA = kFALSE;
  }

 if (evtCutSetID==eventCutSet::kMCEvtDefault) {
   rmFirstEvtChunk = kFALSE;
   rejectPileUp = kFALSE;
   useVtxCut2013pA = kFALSE;
   vtxZcut = 10.0; //cm
   selectDPMJETevtNSDpA = kFALSE;
 }
 
 if (evtCutSetID>=eventCutSet::kpAEvtDefault) {
   triggerMask = AliVEvent::kINT7;
   rmFirstEvtChunk = kTRUE;
   rejectPileUp = kTRUE;
   useVtxCut2013pA = kTRUE;
 }
 
 if (evtCutSetID==eventCutSet::kpANoPileUpCut) {
   rmFirstEvtChunk = kTRUE;
   rejectPileUp = kFALSE;
 }
  
 if (evtCutSetID==eventCutSet::kpAPileUpMV) {
   useMVPileUpSelection = kTRUE;
   MinPlpContribSPD = 3;
   //MinPlpContribMV = 5; //already set as default
 }
  
 if (evtCutSetID==eventCutSet::kpAPileUpSPD3) {
   MinPlpContribSPD = 3;
 }
  
 if (evtCutSetID==eventCutSet::kMCEvtpA2013) {
   rmFirstEvtChunk = kFALSE;
   rejectPileUp = kFALSE;
   useVtxCut2013pA = kTRUE;
   vtxZcut = 1.0e3; //cm
   selectDPMJETevtNSDpA = kFALSE;
 }
  
 if (evtCutSetID==eventCutSet::kMCpAEvtDefault) {
   rmFirstEvtChunk = kFALSE;
   rejectPileUp = kFALSE;
   useVtxCut2013pA = kTRUE;
   vtxZcut = 10.0; //cm
   selectDPMJETevtNSDpA = kFALSE;
 }
 
 if (evtCutSetID>=eventCutSet::kMCEvtNSD) {
   rmFirstEvtChunk = kFALSE;
   rejectPileUp = kFALSE;
   useVtxCut2013pA = kFALSE;
   vtxZcut = 1.0e3; //cm
   selectDPMJETevtNSDpA = kTRUE;
 }
  
 if (evtCutSetID==eventCutSet::kMCEvtNSDpA2013) {
   rmFirstEvtChunk = kFALSE;
   rejectPileUp = kFALSE;
   vtxZcut = 1.0e3; //cm
   selectDPMJETevtNSDpA = kTRUE;
   useVtxCut2013pA = kTRUE;
 }
  
 if (evtCutSetID==eventCutSet::kMCEvtNSDdefault) {
   rmFirstEvtChunk = kFALSE;
   rejectPileUp = kFALSE;
   selectDPMJETevtNSDpA = kTRUE;
   useVtxCut2013pA = kTRUE;
   vtxZcut = 10.0; //cm
 }
  
  //-------------------------------------------
  //pair cuts
  //-------------------------------------------
  Double_t    minYlab =  -0.5;
  Double_t    maxYlab =  0.5;

  if (pairCutSetID==pairYCutSet::kpADefault) { //-0.5<y_cm<0.0
    minYlab = -0.465;    maxYlab = 0.035;
  }
  
  if (pairCutSetID==pairYCutSet::kpANegative) { //-0.5<y_cm<0.0
    minYlab = -0.965;    maxYlab = -0.465;
  }
  
  if (pairCutSetID==pairYCutSet::kpACentral) { //|y_cm|<0.3
    minYlab = -0.765;    maxYlab = -0.165;
  }
  
  if (pairCutSetID==pairYCutSet::kCentralTight) { //|y_cm|<0.3
    minYlab = -0.3;    maxYlab = 0.3;
  }

  //-------------------------------------------
  //mixing settings
  //-------------------------------------------
  Int_t       nmix = 5;
  Float_t     maxDiffVzMix = 1.0;
  Float_t     maxDiffMultMix = 5.0;

  //
  // -- INITIALIZATION ----------------------------------------------------------------------------
  // retrieve analysis manager
  //
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddAnalysisTaskTOFKStar", "No analysis manager to connect to.");
      return NULL;
   } 

   // create the task and configure 
   TString taskName = Form("taskRsnMC");
   AliRsnMiniAnalysisTask *task = new AliRsnMiniAnalysisTask(taskName.Data(), isMC);
   task->UseESDTriggerMask(triggerMask); //ESD
   // task->SelectCollisionCandidates(triggerMask); //AOD
   
   // set event mixing options
   task->UseContinuousMix();
   task->SetNMix(nmix);
   task->SetMaxDiffVz(maxDiffVzMix);
   task->SetMaxDiffMult(maxDiffMultMix);
   ::Info("AddTaskMC", Form("Event mixing configuration: \n events to mix = %i \n max diff. vtxZ = cm %5.3f \n max diff multi = %5.3f", nmix, maxDiffVzMix, maxDiffMultMix));
   mgr->AddTask(task);
   
   //
   // -- EVENT CUTS (same for all configs) ---------------------------------------------------------
   //
   AliRsnCutPrimaryVertex *cutVertex = new AliRsnCutPrimaryVertex("cutVertex", vtxZcut, 0, kFALSE);
   AliRsnCutEventUtils *cutEventUtils = new AliRsnCutEventUtils("cutEventUtils", rmFirstEvtChunk, rejectPileUp);
   cutEventUtils->SetUseVertexSelection2013pA(useVtxCut2013pA, vtxZcut);
   ::Info("AddTaskKStarPPB", Form(":::::::::::::::::: Vertex cut as pA 2013 (max Vz = %4.2f cm): %s", vtxZcut, (useVtxCut2013pA?"ON":"OFF")));  
   if (isMC) {
     cutEventUtils->SetFilterNSDeventsDPMJETpA2013(selectDPMJETevtNSDpA);
     ::Info("AddTaskKStarPPB", Form(":::::::::::::::::: NSD selection in DPMJET pA: %s", (selectDPMJETevtNSDpA?"ON":"OFF")));  
   }
 
   if (useMVPileUpSelection){
     cutEventUtils->SetUseMVPlpSelection(useMVPileUpSelection);
     cutEventUtils->SetMinPlpContribMV(MinPlpContribMV);
     cutEventUtils->SetMinPlpContribSPD(MinPlpContribSPD);
     ::Info("AddTaskKStarPPB", Form("Multiple-vtx Pile-up rejection settings: MinPlpContribMV = %i, MinPlpContribSPD = %i", MinPlpContribMV, MinPlpContribSPD));
   } else {
     cutEventUtils->SetMinPlpContribSPD(MinPlpContribSPD);
     ::Info("AddTaskKStarPPB", Form("SPD Pile-up rejection settings: MinPlpContribSPD = %i", MinPlpContribSPD));
   }
   ::Info("AddTaskKStarPPB", Form(":::::::::::::::::: Pile-up rejection mode: %s", (rejectPileUp?"ON":"OFF")));   
   ::Info("AddTaskKStarPPB", Form("::::::::::::: Remove first event in chunk: %s", (rmFirstEvtChunk?"ON":"OFF")));   
   
   AliRsnCutSet *eventCuts = new AliRsnCutSet("eventCuts", AliRsnTarget::kEvent);
   eventCuts->AddCut(cutEventUtils);
   eventCuts->AddCut(cutVertex);
   eventCuts->SetCutScheme(Form("%s&%s", cutEventUtils->GetName(), cutVertex->GetName()));
   task->SetEventCuts(eventCuts);
   
   //
   // -- EVENT-ONLY COMPUTATIONS -------------------------------------------------------------------
   //   
   //vertex
   Int_t vtxID = task->CreateValue(AliRsnMiniValue::kVz, kFALSE);
   AliRsnMiniOutput *outVtx = task->CreateOutput("eventVtx", "HIST", "EVENT");
   outVtx->AddAxis(vtxID, 500, -50.0, 50.0);
   
   // -- PAIR CUTS (common to all resonances) ------------------------------------------------------
   //
   AliRsnCutMiniPair *cutY = new AliRsnCutMiniPair("cutRapidity", AliRsnCutMiniPair::kRapidityRange);
   cutY->SetRangeD(minYlab, maxYlab);
   
   AliRsnCutSet *cutsPair = new AliRsnCutSet("pairCuts", AliRsnTarget::kMother);
   cutsPair->AddCut(cutY);
   cutsPair->SetCutScheme(cutY->GetName());
   
   //
   // -- CONFIG ANALYSIS --------------------------------------------------------------------------
   //
   // set daughter cuts
   //Use default quality cuts std 2010 with crossed rows TPC
   //Default use loose or no PID selection 
   Bool_t useCrossedRows = 1; 
   Float_t nsigma = 3.0;  
   AliRsnCutSetDaughterParticle * cutSetPi = new AliRsnCutSetDaughterParticle("cutPi", AliRsnCutSetDaughterParticle::kQualityStd2010, AliPID::kPion, nsigma, aodFilterBit, useCrossedRows);
   AliRsnCutSetDaughterParticle * cutSetKa = new AliRsnCutSetDaughterParticle("cutKa", AliRsnCutSetDaughterParticle::kQualityStd2010, AliPID::kKaon, nsigma, aodFilterBit, useCrossedRows);
   AliRsnCutSetDaughterParticle * cutSetPro = new AliRsnCutSetDaughterParticle("cutPro", AliRsnCutSetDaughterParticle::kQualityStd2010, AliPID::kProton, nsigma, aodFilterBit, useCrossedRows);
   Int_t icutPi = task->AddTrackCuts(cutSetPi);
   Int_t icutKa = task->AddTrackCuts(cutSetKa);
   Int_t icutPro = task->AddTrackCuts(cutSetPro);
   
   if (enableMonitor){
    Printf("======== Cut monitoring enabled");
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/AddMonitorOutput.C");
    AddMonitorOutput(isMC, cutSetPi->GetMonitorOutput(), monitorOpt.Data());
    AddMonitorOutput(isMC, cutSetKa->GetMonitorOutput()), monitorOpt.Data();
    AddMonitorOutput(isMC, cutSetPro->GetMonitorOutput()), monitorOpt.Data();
   }  

   //configure output for each rsn
   TString partname [10] = {"Rho","Kstar","aKstar", "Phi", "Delta0","Deltapp","Deltamm", "Lstar", "aLstar", "f0"};
   Int_t    pdgCode [10] = {113, 313, -313, 333, 2114, 2224, -2224, 3124, -3124, 9010221};
   Float_t mass     [10] = {0.770, 0.89495, 0.89445, 1.019461, 1.232, 1.2349, 1.2349, 1.51954, 1.51954, 0.990};
   Float_t masslow  [10] = {0.5, 0.7, 0.7, 0.9, 0.8, 0.8, 0.8, 1.3, 1.3, 0.8};
   Float_t massup   [10] = {1.0, 1.0, 1.0, 1.1, 1.4, 1.4, 1.4, 1.7, 1.7, 1.2};
   Int_t   nbins    [10] = {250, 300, 300,  200, 300, 300, 300, 200, 200, 200};
   Char_t  charge1  [10] = {'+', '+','-','+','+','+','-','+','-','+'};
   Char_t  charge2  [10] = {'-', '-','+','-','-','+','-','-','+','-'};
   RSNPID d1[10] = { AliRsnDaughter::kPion,  AliRsnDaughter::kKaon,  AliRsnDaughter::kKaon, 
		     AliRsnDaughter::kKaon,  AliRsnDaughter::kProton,AliRsnDaughter::kProton,
		     AliRsnDaughter::kProton, AliRsnDaughter::kProton, AliRsnDaughter::kProton, 
		     AliRsnDaughter::kPion};
   RSNPID d2[10] = { AliRsnDaughter::kPion,  AliRsnDaughter::kPion,  AliRsnDaughter::kPion,
		     AliRsnDaughter::kKaon, AliRsnDaughter::kPion, AliRsnDaughter::kPion, 
		     AliRsnDaughter::kPion, AliRsnDaughter::kKaon, AliRsnDaughter::kKaon, 
		     AliRsnDaughter::kPion};
   
   for (Int_t i=0;i<ERsnSpecie_t::kNrsn;i++){
     // set daughter cuts according to decay
     Int_t icut1, icut2;
     if (d1[i]==AliRsnDaughter::kProton) icut2 = icutPro;
     else if (d1[i]==AliRsnDaughter::kKaon) icut2 = icutKa;
     else icut1 = icutPi;
     
     if (d2[i]==AliRsnDaughter::kProton) icut2 = icutPro;
     else if (d2[i]==AliRsnDaughter::kKaon) icut2 = icutKa;
     else icut2 = icutPi;
     
     //configure output
     ConfigSuperOutput(task, isPP, cutsPair, partname[i].Data(), pdgCode[i], mass[i], masslow[i], massup[i], nbins[i], charge1[i], charge2[i], d1[i],d2[i], icut1, icut2);
   }
   
   // -- CONTAINERS --------------------------------------------------------------------------------
   TString outputFileName = AliAnalysisManager::GetCommonFileName();
   Printf("AddSuperTaskMC - Set OutputFileName : \n %s\n", outputFileName.Data() );
   AliAnalysisDataContainer *output = mgr->CreateContainer(Form("RsnOut_%s",outNameSuffix.Data()), 
							   TList::Class(), 
							   AliAnalysisManager::kOutputContainer, 
							   outputFileName);
   mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(task, 1, output);
   
   return task;
}

//-------------------------------------------------------------------
Bool_t ConfigSuperOutput(AliRsnMiniAnalysisTask *task, 
			 Bool_t                 isPP,
			 AliRsnCutSet           *cutsPair,
			 TString                partname,
			 Int_t                  pdgCode,
			 Float_t                mass,
			 Float_t                masslow,
			 Float_t                massup,
			 Int_t                  nbins,
			 Char_t                 charge1,
			 Char_t                 charge2,
			 RSNPID                 d1,
			 RSNPID                 d2,
			 Int_t                  icut1,
			 Int_t                  icut2)
{
  // -- Values ------------------------------------------------------------------------------------
  AliRsnMiniValue::EType yaxisVar = AliRsnMiniValue::kPt;
  /* invariant mass   */ Int_t imID   = task->CreateValue(AliRsnMiniValue::kInvMass, kFALSE);
  /* IM resolution    */ Int_t resID  = task->CreateValue(AliRsnMiniValue::kInvMassRes, kTRUE);
  /* transv. momentum */ Int_t ptID   = task->CreateValue(AliRsnMiniValue::kPt, kFALSE);
  /* pseudorapidity   */ Int_t etaID  = task->CreateValue(AliRsnMiniValue::kEta, kFALSE);
  /* rapidity         */ Int_t yID    = task->CreateValue(AliRsnMiniValue::kY, kTRUE);
  TString output = "HIST";

  //TRUE RECO PAIRS - MASS vs PT
  AliRsnMiniOutput * outtrue = task->CreateOutput(Form("true_pt_%s", partname.Data()), output.Data(),"TRUE");
  outtrue->SetCutID(0, icut1);
  outtrue->SetCutID(1, icut2);
  outtrue->SetCharge(0, charge1);
  outtrue->SetCharge(1, charge2);
  outtrue->SetDaughter(0, d1);
  outtrue->SetDaughter(1, d2);
  outtrue->SetMotherPDG(pdgCode);
  outtrue->SetMotherMass(mass);
  outtrue->SetPairCuts(cutsPair);
  outtrue->AddAxis(ptID, 200, 0.0, 20.0); //default use mother pt
 
  //TRUE RECO PAIRS - MASS vs eta
  AliRsnMiniOutput * outtrueEta = task->CreateOutput(Form("true_eta_%s", partname.Data()), output.Data(),"TRUE");
  outtrueEta->SetCutID(0, icut1);
  outtrueEta->SetCutID(1, icut2);
  outtrueEta->SetCharge(0, charge1);
  outtrueEta->SetCharge(1, charge2);
  outtrueEta->SetDaughter(0, d1);
  outtrueEta->SetDaughter(1, d2);
  outtrueEta->SetMotherPDG(pdgCode);
  outtrueEta->SetMotherMass(mass);
  outtrueEta->SetPairCuts(cutsPair);
  outtrueEta->AddAxis(etaID, 40, -1.0, 1.0); 
  
  //TRUE RECO PAIRS - MASS RESOLUTION
  AliRsnMiniOutput * outres = task->CreateOutput(Form("res_pt_%s", partname.Data()), output.Data(),"TRUE");
  outres->SetCutID(0, icut1);
  outres->SetCutID(1, icut2);
  outres->SetCharge(0, charge1);
  outres->SetCharge(1, charge2);
  outres->SetDaughter(0, d1);
  outres->SetDaughter(1, d2);
  outres->SetMotherPDG(pdgCode);
  outres->SetMotherMass(mass);
  outres->SetPairCuts(cutsPair);
  outres->AddAxis(resID, 200, -0.01, 0.01);
  outres->AddAxis(ptID, 200, 0.0, 20.0);
     
  //GENERATED PAIRS
  AliRsnMiniOutput * outm = task->CreateOutput(Form("mother_pt_%s", partname.Data()), output.Data(),"MOTHER");
  outm->SetDaughter(0, d1);
  outm->SetDaughter(1, d2);
  outm->SetMotherPDG(pdgCode);
  outm->SetMotherMass(mass);
  outm->SetPairCuts(cutsPair);
  outm->AddAxis(ptID, 200, 0.0, 20.0);

  AliRsnMiniOutput * outmEta = task->CreateOutput(Form("mother_eta_%s", partname.Data()), output.Data(),"MOTHER");
  outmEta->SetDaughter(0, d1);
  outmEta->SetDaughter(1, d2);
  outmEta->SetMotherPDG(pdgCode);
  outmEta->SetMotherMass(mass);
  outmEta->SetPairCuts(cutsPair);
  outmEta->AddAxis(etaID, 40, -1., 1.);
  
  AliRsnMiniOutput * outmY = task->CreateOutput(Form("mother_y_%s", partname.Data()), output.Data(),"MOTHER");
  outmY->SetDaughter(0, d1);
  outmY->SetDaughter(1, d2);
  outmY->SetMotherPDG(pdgCode);
  outmY->SetMotherMass(mass);
  outmY->SetPairCuts(cutsPair);
  outmY->AddAxis(yID, 40, -1.0, 1.0); 

  return kTRUE;
}
