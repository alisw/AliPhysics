/***************************************************************************
              Anders Knospe - last modified on 26 March 2016

//Lauches phi analysis with rsn mini package
//Allows basic configuration of pile-up check and event cuts
****************************************************************************/

enum pairYCutSet { kPairDefault=0,
		   kCentral //=1
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
		      kMixDefault,//=0 //10 events, Dvz = 1cm, DC = 10
		      k5Evts, //=1 //5 events, Dvz = 1cm, DC = 10
		      k5Cent,  //=2 //10 events, Dvz = 1cm, DC = 5
		      k5Evts5Cent
                    };

AliRsnMiniAnalysisTask * AddTaskPhiPP5TeV_PID
(
 Bool_t      isMC=kFALSE,
 Bool_t      isPP=kTRUE,
 TString     outNameSuffix="tpc2stof3sveto",
 Int_t       evtCutSetID=0,
 Int_t       pairCutSetID=0,
 Int_t       mixingConfigID=0,
 Int_t       aodFilterBit=5,
 Int_t       customQualityCutsID=1,
 AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutKaCandidate=AliRsnCutSetDaughterParticle::kTPCpidTOFveto3s,
 Float_t     nsigmaKa=2.,
 Bool_t      enableMonitor=kTRUE,
 Bool_t      IsMcTrueOnly=kFALSE,
 TString     monitorOpt="NoSIGN",
 Bool_t      useMixLS=0,
 Bool_t      checkReflex=0,
 AliRsnMiniValue::EType yaxisvar=AliRsnMiniValue::kPt,
 TString     polarizationOpt="" /* J - Jackson,T - Transversity */
)
{  
  //-------------------------------------------
  // event cuts
  //-------------------------------------------
  UInt_t      triggerMask=AliVEvent::kINT7;
  Bool_t      rejectPileUp=kTRUE;
  Double_t    vtxZcut=10.0;//cm, default cut on vtx z
  Int_t       MultBins=aodFilterBit/100;

  if(evtCutSetID==eventCutSet::kDefaultVtx12) vtxZcut=12.0; //cm
  if(evtCutSetID==eventCutSet::kDefaultVtx8) vtxZcut=8.0; //cm
  if(evtCutSetID==eventCutSet::kDefaultVtx5) vtxZcut=5.0; //cm
  if(evtCutSetID==eventCutSet::kNoPileUpCut) rejectPileUp=kFALSE;
  if(evtCutSetID==eventCutSet::kSpecial2) vtxZcut=1.e6;//off

  if(!isPP || isMC || MultBins) rejectPileUp=kFALSE;

  //-------------------------------------------
  //pair cuts
  //-------------------------------------------
  Double_t    minYlab=-0.5;
  Double_t    maxYlab= 0.5;

  if(pairCutSetID==pairYCutSet::kCentral){//|y_cm|<0.3
    minYlab=-0.3; maxYlab=0.3;
  }

  //-------------------------------------------
  //mixing settings
  //-------------------------------------------
  Int_t       nmix=0;
  Float_t     maxDiffVzMix=1.;
  Float_t     maxDiffMultMix=10.;

  if(mixingConfigID==eventMixConfig::kMixDefault) nmix=10;
  if(mixingConfigID==eventMixConfig::k5Evts) nmix=5;
  if(mixingConfigID==eventMixConfig::k5Cent) maxDiffMultMix=5;
  if(mixingConfigID==eventMixConfig::k5Evts5Cent){nmix=5; maxDiffMultMix=5;}

  // -- INITIALIZATION ----------------------------------------------------------------------------
  // retrieve analysis manager

  AliAnalysisManager* mgr=AliAnalysisManager::GetAnalysisManager();
  if(!mgr){
    ::Error("AddTaskPhiPP5TeV_PID", "No analysis manager to connect to.");
    return NULL;
  }

  // create the task and configure 
  TString taskName=Form("phi%s%s_%i%i",(isPP? "pp" : "PbPb"),(isMC ? "MC" : "Data"),(Int_t)cutKaCandidate);
  AliRsnMiniAnalysisTask* task=new AliRsnMiniAnalysisTask(taskName.Data(),isMC);
  if(evtCutSetID==eventCutSet::kSpecial4 || evtCutSetID==eventCutSet::kSpecial5) task->UseESDTriggerMask(triggerMask); //ESD ****** check this *****
  if(evtCutSetID!=eventCutSet::kNoEvtSel && evtCutSetID!=eventCutSet::kSpecial3 && evtCutSetID!=eventCutSet::kSpecial4) task->SelectCollisionCandidates(triggerMask); //AOD

  if(isPP){
    if(MultBins==1) task->UseMultiplicity("AliMultSelection_V0M");
    else if(MultBins==2) task->UseMultiplicity("AliMultSelection_RefMult08");
    else task->UseMultiplicity("QUALITY");
  }else task->UseCentrality("V0M");

  // set event mixing options
  task->UseContinuousMix();
  //task->UseBinnedMix();
  task->SetNMix(nmix);
  task->SetMaxDiffVz(maxDiffVzMix);
  task->SetMaxDiffMult(maxDiffMultMix);
  ::Info("AddTaskPhiPP5TeV_PID", Form("Event mixing configuration: \n events to mix = %i \n max diff. vtxZ = cm %5.3f \n max diff multi = %5.3f", nmix, maxDiffVzMix, maxDiffMultMix));

  mgr->AddTask(task);

  // cut on primary vertex:
  // - 2nd argument --> |Vz| range
  // - 3rd argument --> minimum required number of contributors to vtx
  // - 4th argument --> tells if TPC stand-alone vertexes must be accepted

  AliRsnCutPrimaryVertex* cutVertex=0;
  if(evtCutSetID!=eventCutSet::kSpecial1 && evtCutSetID!=eventCutSet::kNoEvtSel &&  (!MultBins || fabs(vtxZcut-10.)>1.e-10)){
    cutVertex=new AliRsnCutPrimaryVertex("cutVertex",vtxZcut,0,kFALSE);
    if(!MultBins && evtCutSetID!=eventCutSet::kSpecial3){
    cutVertex->SetCheckZResolutionSPD();
    cutVertex->SetCheckDispersionSPD();
    cutVertex->SetCheckZDifferenceSPDTrack();
    }
   if (evtCutSetID==eventCutSet::kSpecial3) cutVertex->SetCheckGeneratedVertexZ();
  }

  AliRsnCutEventUtils* cutEventUtils=0;
  if(evtCutSetID!=eventCutSet::kNoEvtSel && evtCutSetID!=eventCutSet::kSpecial3){
    cutEventUtils=new AliRsnCutEventUtils("cutEventUtils",kTRUE,rejectPileUp);
    if(!MultBins){
      cutEventUtils->SetCheckIncompleteDAQ();
      cutEventUtils->SetCheckSPDClusterVsTrackletBG();
    }else{
      //cutEventUtils->SetCheckInelGt0SPDtracklets();
      cutEventUtils->SetRemovePileUppA2013(kFALSE);
      cutEventUtils->SetCheckAcceptedMultSelection();
    }
  }

  if(isPP && (!isMC) && cutVertex){ 
    cutVertex->SetCheckPileUp(rejectPileUp);// set the check for pileup  
    ::Info("AddTaskPhiPP5TeV_PID", Form(":::::::::::::::::: Pile-up rejection mode: %s", (rejectPileUp)?"ON":"OFF"));
  }

  // define and fill cut set for event cut
  AliRsnCutSet* eventCuts=0;
  if(cutEventUtils || cutVertex){
    eventCuts=new AliRsnCutSet("eventCuts",AliRsnTarget::kEvent);

    if(cutEventUtils && cutVertex){
      eventCuts->AddCut(cutEventUtils);
      eventCuts->AddCut(cutVertex);
      eventCuts->SetCutScheme(Form("%s&%s",cutEventUtils->GetName(),cutVertex->GetName()));
    }else if(cutEventUtils && !cutVertex){
      eventCuts->AddCut(cutEventUtils);
      eventCuts->SetCutScheme(Form("%s",cutEventUtils->GetName()));
    }else if(!cutEventUtils && cutVertex){
      eventCuts->AddCut(cutVertex);
      eventCuts->SetCutScheme(Form("%s",cutVertex->GetName()));
    }

    task->SetEventCuts(eventCuts);
  }

  //-- EVENT-ONLY COMPUTATIONS -------------------------------------------------------------------

  //vertex
  Int_t vtxID=task->CreateValue(AliRsnMiniValue::kVz,kFALSE);
  AliRsnMiniOutput* outVtx=task->CreateOutput("eventVtx","HIST","EVENT");
  outVtx->AddAxis(vtxID,240,-12.0,12.0);

  //multiplicity or centrality
  Int_t multID=task->CreateValue(AliRsnMiniValue::kMult,kFALSE);
  AliRsnMiniOutput* outMult=task->CreateOutput("eventMult","HIST","EVENT");
  if(isPP && !MultBins) outMult->AddAxis(multID,400,0.5,400.5);
  else outMult->AddAxis(multID,110,0.,110.);

  TH2F* hvz=new TH2F("hVzVsCent","",110,0.,110., 240,-12.0,12.0);
  task->SetEventQAHist("vz",hvz);//plugs this histogram into the fHAEventVz data member

  TH2F* hmc=new TH2F("MultiVsCent","", 110,0.,110., 400,0.5,400.5);
  hmc->GetYaxis()->SetTitle("QUALITY");
  task->SetEventQAHist("multicent",hmc);//plugs this histogram into the fHAEventMultiCent data member

  // -- PAIR CUTS (common to all resonances) ------------------------------------------------------

  AliRsnCutMiniPair* cutY=new AliRsnCutMiniPair("cutRapidity",AliRsnCutMiniPair::kRapidityRange);
  cutY->SetRangeD(minYlab,maxYlab);

  AliRsnCutSet* cutsPair=new AliRsnCutSet("pairCuts",AliRsnTarget::kMother);
  cutsPair->AddCut(cutY);
  cutsPair->SetCutScheme(cutY->GetName());

  // -- CONFIG ANALYSIS --------------------------------------------------------------------------

  gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/ConfigPhiPP5TeV_PID.C");
  if(!ConfigPhiPP5TeV_PID(task,isMC,isPP,"",cutsPair,aodFilterBit,customQualityCutsID,cutKaCandidate,nsigmaKa,enableMonitor,isMC&IsMcTrueOnly,monitorOpt.Data(),useMixLS,isMC&checkReflex,yaxisvar,polarizationOpt)) return 0x0;

  // -- CONTAINERS --------------------------------------------------------------------------------

  TString outputFileName = AliAnalysisManager::GetCommonFileName();
  //  outputFileName += ":Rsn";
  Printf("AddAnalysisTaskPhiPP5TeV_PID - Set OutputFileName : \n %s\n",outputFileName.Data());

  AliAnalysisDataContainer* output=mgr->CreateContainer(Form("RsnOut_%s",outNameSuffix.Data()),
							TList::Class(),
							AliAnalysisManager::kOutputContainer,
							outputFileName);
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, output);

  return task;
}

