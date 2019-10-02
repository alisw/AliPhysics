/***************************************************************************
fbellini@cern.ch - last modified on 04/09/2018
****************************************************************************/

AliRsnMiniAnalysisTask * AddSuperTaskMC(Int_t   collSystem    = AliPIDResponse::kPBPB,
					UInt_t  triggerMask   = AliVEvent::kINT7,
					Int_t   aodFilterBit  = 5,
					Bool_t  enableMonitor = kFALSE,
					TString outNameSuffix = "")
{  

  // -------------------------------------

  
  // ADD task for MC QA - simple settings
  // --------------------------------------
  
  // event cuts
  Double_t    vtxZcut = 10.0; //cm, default cut on vtx z
  TString     monitorOpt = "NoSIGN";
  
  //pair cuts
  Double_t    minYlab = -0.5;
  Double_t    maxYlab =  0.5;

  if (collSystem == AliPIDResponse::kPPB) {
    //-0.5<y_cm<0.0
    minYlab = -0.465;
    maxYlab = 0.035;
  }
  
  //
  // -- INITIALIZATION ----------------------------------------------------------------------------
  // retrieve analysis manager
  //
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddSuperTaskMC", "No analysis manager to connect to.");
    return NULL;
  } 

  // create the task and configure
  Bool_t isMC = kTRUE;
  TString taskName = Form("RsnTaskMCqa");
  AliRsnMiniAnalysisTask *task = new AliRsnMiniAnalysisTask(taskName.Data(), isMC);
  task->UseMultiplicity("AliMultSelection_V0M");
  //task->UseESDTriggerMask(triggerMask); //ESD
  //task->SelectCollisionCandidates(triggerMask); //AOD
  mgr->AddTask(task);

  //
  // -- EVENT CUTS (same for all configs) ---------------------------------------------------------
  //
  AliRsnCutEventUtils* cutEventUtils = new AliRsnCutEventUtils("cutEventUtils", kFALSE, kFALSE);
  cutEventUtils->SetCheckAcceptedMultSelection(); 
  AliRsnCutPrimaryVertex *cutVertex = new AliRsnCutPrimaryVertex("cutVertex", vtxZcut, 0, kFALSE);
  AliRsnCutSet *eventCuts = new AliRsnCutSet("eventCuts", AliRsnTarget::kEvent);
  eventCuts->AddCut(cutVertex);
  eventCuts->AddCut(cutEventUtils);
  eventCuts->SetCutScheme(Form("%s&%s", cutVertex->GetName(), cutEventUtils->GetName()));
  task->SetEventCuts(eventCuts);
   
  //
  // -- PAIR CUTS (common to all resonances) ------------------------------------------------------
  //
  AliRsnCutMiniPair *cutY = new AliRsnCutMiniPair("cutRapidity", AliRsnCutMiniPair::kRapidityRange);
  cutY->SetRangeD(minYlab, maxYlab);
   
  AliRsnCutSet *cutsPair = new AliRsnCutSet("pairCuts", AliRsnTarget::kMother);
  //cutsPair->AddCut(cutY);
  //cutsPair->SetCutScheme(cutY->GetName());
  
  //
  // SET EVENT QA MONITORING
  //
  TH1F* hEventsVsMulti = new TH1F("hAEventsVsMulti", "; multiplicity (%); events", 101, 0, 101.0);
  task->SetEventQAHist("EventsVsMulti", hEventsVsMulti);//custom binning for fHAEventsVsMulti
  
  TH2F* hAEventVzCent = new TH2F("hAEventVzCent", "Accepted events Vz vs multiplicity ; multiplicity (%); Vz (cm)", 101, 0, 101.0, 400, -20., 20.);
  task->SetEventQAHist("vz", hAEventVzCent);//custom binning for fHAEventsVsMulti
  
  //
  // -- CONFIG ANALYSIS --------------------------------------------------------------------------
  //
  // set daughter cuts
  // Use default quality cuts std 2011 with crossed rows TPC
  // Default use loose or no PID selection
  
  Bool_t useCrossedRows = 1; 
  Float_t nsigma = 3.0;  
  AliRsnCutSetDaughterParticle * cutSetPi = new AliRsnCutSetDaughterParticle("cutPi", AliRsnCutSetDaughterParticle::kQualityStd2011, AliPID::kPion, nsigma, aodFilterBit, useCrossedRows);
  AliRsnCutSetDaughterParticle * cutSetKa = new AliRsnCutSetDaughterParticle("cutKa", AliRsnCutSetDaughterParticle::kQualityStd2011, AliPID::kKaon, nsigma, aodFilterBit, useCrossedRows);
  AliRsnCutSetDaughterParticle * cutSetPro = new AliRsnCutSetDaughterParticle("cutPro", AliRsnCutSetDaughterParticle::kQualityStd2011, AliPID::kProton, nsigma, aodFilterBit, useCrossedRows);
  Int_t icutPi = task->AddTrackCuts(cutSetPi);
  Int_t icutKa = task->AddTrackCuts(cutSetKa);
  Int_t icutPro = task->AddTrackCuts(cutSetPro);
  
  // if (enableMonitor){
  //   Printf("======== Cut monitoring enabled");
  //   gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/AddMonitorOutput.C");
  //   AddMonitorOutput(isMC, cutSetPi->GetMonitorOutput(), monitorOpt.Data());
  //   AddMonitorOutput(isMC, cutSetKa->GetMonitorOutput()), monitorOpt.Data();
  //   AddMonitorOutput(isMC, cutSetPro->GetMonitorOutput()), monitorOpt.Data();
  // }  
 
  //configure output for each rsn
  const Int_t Np = 10;
  TString partname [Np] = {"Rho","Kstar","aKstar", "Phi", "Lstar", "aLstar", "f0", "f2", "K2star0","aK2star0"};
  Int_t    pdgCode [Np] = {113, 313, -313, 333, 3124, -3124, 9010221, 225, 315, -315 };
  Float_t mass     [Np] = {0.770, 0.89495, 0.89445, 1.019461, 1.51954, 1.51954, 0.990, 1.270, 1.432, 1.432};
  Float_t masslow  [Np] = {0.5, 0.7, 0.7, 0.9, 1.3, 1.3, 0.8, 0.8,1.0,1.0};
  Float_t massup   [Np] = {1.0, 1.0, 1.0, 1.1, 1.7, 1.7, 1.2, 1.6,2.0,2.0};
  Char_t  charge1  [Np] = {'+', '+','-','+','+','-','+','+','+','-'};
  Char_t  charge2  [Np] = {'-', '-','+','-','-','+','-','-','-','+'};
  
  RSNPID d1[Np] = { AliRsnDaughter::kPion, AliRsnDaughter::kKaon,  AliRsnDaughter::kKaon, 
		    AliRsnDaughter::kKaon, AliRsnDaughter::kProton, AliRsnDaughter::kProton, 
		    AliRsnDaughter::kPion, AliRsnDaughter::kPion, AliRsnDaughter::kKaon, AliRsnDaughter::kPion};
  RSNPID d2[Np] = { AliRsnDaughter::kPion, AliRsnDaughter::kPion,  AliRsnDaughter::kPion,
		    AliRsnDaughter::kKaon, AliRsnDaughter::kKaon, AliRsnDaughter::kKaon, 
		    AliRsnDaughter::kPion, AliRsnDaughter::kPion, AliRsnDaughter::kPion, AliRsnDaughter::kKaon};

  // invariant mass
  Int_t imID   = task->CreateValue(AliRsnMiniValue::kInvMass, kFALSE);
  // IM resolution
  Int_t resID  = task->CreateValue(AliRsnMiniValue::kInvMassRes, kTRUE);
  // transv. momentum
  Int_t ptID   = task->CreateValue(AliRsnMiniValue::kPt, kFALSE);
  // pseudorapidity
  Int_t etaID  = task->CreateValue(AliRsnMiniValue::kEta, kFALSE);
  // rapidity
  Int_t yID    = task->CreateValue(AliRsnMiniValue::kY, kTRUE);

  AliRsnMiniOutput * out[Np][3];
  AliRsnMiniOutput * outG[Np][2];
  TString axes[3] = {"invm_pt", "pt_rap", "res_pt"};
  // TString axesGen[3] = {"invm_pt", "pt_rap", "res_pt"};
  
  for (Int_t i = 0; i<10; i++){
    Int_t icut1, icut2;
    if (d1[i]==AliRsnDaughter::kProton) icut2 = icutPro;
    else if (d1[i]==AliRsnDaughter::kKaon) icut2 = icutKa;
    else icut1 = icutPi;
    
    if (d2[i]==AliRsnDaughter::kProton) icut2 = icutPro;
    else if (d2[i]==AliRsnDaughter::kKaon) icut2 = icutKa;
    else icut2 = icutPi;

    Int_t nbins = (massup[i]-masslow[i])*1000;
   
    //configure output
    for (Int_t j = 0; j<3; j++){
      out[i][j] = task->CreateOutput(Form("%s_true_%s", partname[i].Data(), axes[j].Data()), "HIST","TRUE");
      out[i][j]->SetCutID(0, icut1);
      out[i][j]->SetCutID(1, icut2);
      out[i][j]->SetCharge(0, charge1[i]);
      out[i][j]->SetCharge(1, charge2[i]);
      out[i][j]->SetDaughter(0, d1[i]);
      out[i][j]->SetDaughter(1, d2[i]);
      out[i][j]->SetMotherPDG(pdgCode[i]);
      out[i][j]->SetMotherMass(mass[i]);
      out[i][j]->SetPairCuts(cutsPair);
      if (axes[j].Contains("invm")) out[i][j]->AddAxis(imID, nbins, masslow[i], massup[i]);
      if (axes[j].Contains("res")) out[i][j]->AddAxis(resID, 400, -0.02, 0.02);
      if (axes[j].Contains("pt")) out[i][j]->AddAxis(ptID, 150, 0.0, 15.0);
      if (axes[j].Contains("rap")) out[i][j]->AddAxis(yID, 20, -1.0, 1.0);
      if (axes[j].Contains("eta")) out[i][j]->AddAxis(etaID, 40, -2.0, 2.0);
      
      if (j > 1) continue;

      outG[i][j] = task->CreateOutput(Form("%s_gen_%s", partname[i].Data(), axes[j].Data()), "HIST","MOTHER");
      outG[i][j]->SetDaughter(0, d1[i]);
      outG[i][j]->SetDaughter(1, d2[i]);
      outG[i][j]->SetCharge(0, charge1[i]);
      outG[i][j]->SetCharge(1, charge2[i]);
      outG[i][j]->SetMotherPDG(pdgCode[i]);
      outG[i][j]->SetMotherMass(mass[i]);
      outG[i][j]->SetPairCuts(cutsPair);
      if (axes[j].Contains("invm")) outG[i][j]->AddAxis(imID, nbins, masslow[i], massup[i]);
      if (axes[j].Contains("pt")) outG[i][j]->AddAxis(ptID, 150, 0.0, 15.0);
      if (axes[j].Contains("rap")) outG[i][j]->AddAxis(yID, 20, -1.0, 1.0);
      if (axes[j].Contains("eta")) outG[i][j]->AddAxis(etaID, 40, -2.0, 2.0);
    }
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
