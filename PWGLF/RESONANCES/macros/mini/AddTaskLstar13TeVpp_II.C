/***************************************************************************
    Modified by himani.bhatt@cern.ch  - last modified on 30/04/2018 

   priyanka.sett@cern.ch - last modified on 27/10/2016
  
 

  for L* in pp 13 TeV analysis


// General macro to configure the RSN analysis task.
// It calls all configs desired by the user, by means
// of the boolean switches defined in the first lines.
// ---
// Inputs:
//  1) flag to know if running on MC or data
//  2) path where all configs are stored
// ---
// Returns:
//  kTRUE  --> initialization successful
//  kFALSE --> initialization failed (some config gave errors)
//
****************************************************************************/
enum pairYCutSet { kPairDefault=0,
		   kCentral //=1
                 };
enum eventCutSet { kDefaultVtx = 0,
		   kDefaultVtx9, //=1
		   kDefaultVtx8, //=2
		   kDefaultVtx12, //=3
		   kDefaultVtx11,  //=4
		   kNoPileUpCut, //=5                 
		   kNoEvtSel, //=6   //No event selection, only INEL events  
		   kSpecial2,//=7   // No Vz cut on vtx, f_vtx = kSpecial2/default_data
		   kSpecial3 //= 8  // No event selection, INEL events with Vz cut on vtx  kSpecial3/default_MC = f_SL
};

enum eventMixConfig { kDisabled = -1,
		      k5Evts5Cent, //=0 //5 events, Dvz = 1cm, DC = 5
		      kMixDefault, //=1 //10 events, Dvz = 1cm, DC = 10
		      k5Evts, //=2 //5 events, Dvz = 1cm, DC = 10
		      k5Cent  //=3 //10 events, Dvz = 1cm, DC = 5
		     
};

AliRsnMiniAnalysisTask * AddTaskLstar13TeVpp_II
(
 Bool_t      isMC = kFALSE,
 Bool_t      isPP = kTRUE,
 Int_t       aodFilterBit = 5,
 Int_t       pairCutSetID=0,
 Int_t       evtCutSetID = 0,
 Int_t       mixingConfigID = 0,
 Int_t       MultBins = 0,// default for V0_M and MultBins = 2 for RefMult0_8 
 Int_t       customQualityCutsID=1, // for default
 Float_t     nsigmaPr = 2.0,
 Float_t     nsigmaKa = 2.0,
 Bool_t      enableMonitor = kTRUE,
 Bool_t      IsMcTrueOnly = kFALSE,
 TString     outNameSuffix = "Default"
 )

{

  AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutPrCandidate = AliRsnCutSetDaughterParticle::kTPCTOFpidLstar13ppTeV;
  AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutKaCandidate = AliRsnCutSetDaughterParticle::kTPCTOFpidLstar13ppTeV;
  
  Int_t       signedPdg = 3124;
  TString     monitorOpt = "NoSIGN";  //Flag for AddMonitorOutput.C e.g."NoSIGN"
  Bool_t      useCrossedRows = kTRUE;
  const char *yaxisVar = "";  //yaxisVar = "PtDaughter_PDaughter_cent"
  Bool_t      useMixLS = 0;
  //-------------------------------------------
  // event cuts
  //-------------------------------------------
  
  UInt_t      triggerMask = AliVEvent::kINT7;
  if(evtCutSetID>=100){
    triggerMask=AliVEvent::kHighMultV0;
    evtCutSetID=evtCutSetID%100;
  }

  Bool_t      rejectPileUp=kTRUE;
  if(!isPP || isMC || MultBins) rejectPileUp=kFALSE;

  Double_t  vtxZcut = 10.0;//default cut on vtx z
  
  if(evtCutSetID==eventCutSet::kDefaultVtx12) vtxZcut=12.0; //cm
  if(evtCutSetID==eventCutSet::kDefaultVtx11) vtxZcut=11.0; //cm
  if(evtCutSetID==eventCutSet::kDefaultVtx9) vtxZcut=9.0; //cm
  if(evtCutSetID==eventCutSet::kDefaultVtx8) vtxZcut=8.0; //cm
  if(evtCutSetID==eventCutSet::kNoPileUpCut) rejectPileUp=kFALSE;
  if(evtCutSetID==eventCutSet::kSpecial2) vtxZcut=1.e6;//off


  //-------------------------------------------
  //mixing settings
  //-------------------------------------------
  Int_t       nmix=10;
  Float_t     maxDiffVzMix=1.;
  Float_t     maxDiffMultMix=10.;
  
  if(mixingConfigID==eventMixConfig::kMixDefault) nmix=10;
  if(mixingConfigID==eventMixConfig::k5Evts) nmix=5;
  if(mixingConfigID==eventMixConfig::k5Cent) maxDiffMultMix=5;
  if(mixingConfigID==eventMixConfig::k5Evts5Cent) { maxDiffMultMix=5; nmix = 5;}
 
//-------------------------------------------
  //pair cuts
  //-------------------------------------------
  Double_t    minYlab=-0.5;
  Double_t    maxYlab= 0.5;

  if(pairCutSetID==pairYCutSet::kCentral){//|y_cm|<0.3
    minYlab=-0.3; maxYlab=0.3;
  } 
  
  //
  // -- INITIALIZATION ----------------------------------------------------------------------------
  // retrieve analysis manager
  //
  
  
  TString outputFileName = Form("%s", AliAnalysisManager::GetCommonFileName());
  // Objects name
  
  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddAnalysisTaskTPCLStarSyst", "No analysis manager to connect to.");
    return NULL;
  } 
  
  
  // create the task and configure 
  TString taskName = outNameSuffix.Data();
  AliRsnMiniAnalysisTask *task = new AliRsnMiniAnalysisTask(taskName.Data(), isMC);
  
  
  task->UseESDTriggerMask(triggerMask);  // for ESD 
  // task->SelectCollisionCandidates(triggerMask);  // Priyanka used for AODs
  
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
  //if (!isPP) task->SetMaxDiffAngle(maxDiffAngleMixDeg*TMath::DegToRad()); //set angle diff in rad
  ::Info("AddAnalysisTasLStar", Form("Event mixing configuration: \n events to mix = %i \n max diff. vtxZ = cm %5.3f \n max diff multi = %5.3f \n ", nmix, maxDiffVzMix, maxDiffMultMix));
  
  mgr->AddTask(task);
  
  //
  // -- EVENT CUTS (same for all configs) ---------------------------------------------------------
  //  
  // cut on primary vertex:
  // - 2nd argument --> |Vz| range
  // - 3rd argument --> minimum required number of contributors
  // - 4th argument --> tells if TPC stand-alone vertexes must be accepted
  
  
  // vertex cuts
  
  Bool_t rejectPileUp=kTRUE;
  if (!isPP && !isMC) rejectPileUp=kFALSE;
  AliRsnCutPrimaryVertex* cutVertex=0;
  if(evtCutSetID!=eventCutSet::kNoEvtSel &&  (MultBins == 0 || fabs(vtxZcut-10.)>1.e-10))
    { // passes the vertex quality cut
      cutVertex=new AliRsnCutPrimaryVertex("cutVertex",vtxZcut,0,kFALSE);
      if(MultBins == 0 && evtCutSetID!=eventCutSet::kSpecial3){
	cutVertex->SetCheckZResolutionSPD();
	cutVertex->SetCheckDispersionSPD();
	cutVertex->SetCheckZDifferenceSPDTrack();
      }
      if(evtCutSetID==eventCutSet::kSpecial3) cutVertex->SetCheckGeneratedVertexZ();  // only Vz cut on vtx for INEL events
    }
  
  
  
  // other event selection cuts
  AliRsnCutEventUtils* cutEventUtils=0;
  if(evtCutSetID!=eventCutSet::kNoEvtSel && evtCutSetID!=eventCutSet::kSpecial3)
    {
      cutEventUtils=new AliRsnCutEventUtils("cutEventUtils",kTRUE,rejectPileUp);
      if(MultBins == 0){
	cutEventUtils->SetCheckIncompleteDAQ();
	cutEventUtils->SetCheckSPDClusterVsTrackletBG();
      }else{
	cutEventUtils->SetRemovePileUppA2013(kFALSE);
	cutEventUtils->SetCheckAcceptedMultSelection();
      }
    }
  
  
  // set the check for pileup 
  if(isPP && (!isMC) && cutVertex){ 
    cutVertex->SetCheckPileUp(rejectPileUp); 
    ::Info("AddAnalysisTaskLStar", Form(":::::::::::::::::: Pile-up rejection mode: %s", (rejectPileUp)?"ON":"OFF"));
  }
  
  // define and fill cut set for event cuts
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
  
  
  
  //
  // -- EVENT-ONLY COMPUTATIONS -------------------------------------------------------------------
  //   
  //vertex
  Int_t vtxID = task->CreateValue(AliRsnMiniValue::kVz, kFALSE);
  AliRsnMiniOutput *outVtx = task->CreateOutput("eventVtx", "HIST", "EVENT");
  outVtx->AddAxis(vtxID, 400, -20.0, 20.0);
  
  //multiplicity or centrality
  Int_t multID = task->CreateValue(AliRsnMiniValue::kMult, kFALSE);
  AliRsnMiniOutput *outMult = task->CreateOutput("eventMult", "HIST", "EVENT");
  if (isPP && !MultBins)   outMult->AddAxis(multID, 400, 0.0, 400.0);
  else  outMult->AddAxis(multID, 110, 0.0, 110.0);
  
  Double_t multbins[200];
  int j,nmult=0;
  for(j=0;j<10;j++){multbins[nmult]=0.0001*j; nmult++;}
  for(j=1;j<10;j++){multbins[nmult]=0.001*j; nmult++;}
  for(j=1;j<50;j++){multbins[nmult]=0.01*j; nmult++;}
  for(j=5;j<10;j++){multbins[nmult]=0.1*j; nmult++;}
  for(j=1;j<=100;j++){multbins[nmult]=j; nmult++;}
  nmult--;
  
  TH1F* hEventsVsMulti=new TH1F("hAEventsVsMulti","",nmult,multbins);
  // task->SetEventQAHist("EventsVsMulti",hEventsVsMulti);//custom binning for fHAEventsVsMulti

  TH2F* hvz=new TH2F("hVzVsCent","",110,0.,110., 240,-12.0,12.0);
  task->SetEventQAHist("vz",hvz);//plugs this histogram into the fHAEventVz data member

  double ybins[500];
  for(j=0;j<=401;j++) ybins[j]=j-0.5;

  TH2F* hmc=new TH2F("MultiVsCent","", nmult,multbins, 401,ybins);
  hmc->GetYaxis()->SetTitle("QUALITY");
  task->SetEventQAHist("multicent",hmc);//plugs this histogram into the fHAEventMultiCent data member

  //  <<==============
  
  
  //
  // -- PAIR CUTS  -------------------------------------------------------------------------------
  //
  AliRsnCutMiniPair* cutY=new AliRsnCutMiniPair("cutRapidity",AliRsnCutMiniPair::kRapidityRange);
  cutY->SetRangeD(minYlab,maxYlab);

  AliRsnCutSet* cutsPair=new AliRsnCutSet("pairCuts",AliRsnTarget::kMother);
  cutsPair->AddCut(cutY);
  cutsPair->SetCutScheme(cutY->GetName());

  //
  // -- CONFIG ANALYSIS --------------------------------------------------------------------------
  
  
  //for systematic checks
  {
    
    // gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/ConfigureLstar13TeVpp_II.C");
    gROOT->LoadMacro("ConfigureLstar13TeVpp_II.C");
    if (!ConfigureLstar13TeVpp_II(task, isMC, isPP, "", cutsPair, aodFilterBit, customQualityCutsID, cutPrCandidate, cutKaCandidate, nsigmaPr, nsigmaKa,  enableMonitor, isMC&IsMcTrueOnly, signedPdg, monitorOpt, useCrossedRows, yaxisVar ,useMixLS)) 
      return 0x0;  
  }
  
  
  //
  // -- CONTAINERS --------------------------------------------------------------------------------
  //
  TString outputFileName = AliAnalysisManager::GetCommonFileName();
  Printf("AddAnalysisTaskLstar13TeVpp - Set OutputFileName : \n %s\n", outputFileName.Data() );
  
  AliAnalysisDataContainer *output = mgr->CreateContainer(Form("RsnOut_%s",outNameSuffix.Data()),TList::Class(),AliAnalysisManager::kOutputContainer, outputFileName);
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, output);
  cout<<" taskname  =  "<<taskName.Data()<<endl;
  return task;
}
