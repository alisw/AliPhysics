#if !defined (__CINT__) || defined (__CLING__)
R__ADD_INCLUDE_PATH($ALICE_PHYSICS)
#include <PWGLF/RESONANCES/macros/mini/ConfigKStarMult13TeVpp_AOD.C>
#endif
/***************************************************************************
Modified on: 09/05/2023
By: Sonali Padhan

 Modified By: Pragati Sahoo on: 04/08/2019
 

//Launches Kstar analysis with rsn mini package
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
           kTriggered, //=6
           kNoVzCut, //=7
           kNoEvtSel, //=8
           kINEL10, //=9
           kIGZ10, //=10
           kIGZ //=11
                 };
enum eventMixConfig { kDisabled = -1,
kMixDefault,//=0 //10 events, Dvz = 1cm, DC = 10
k5Evts, //=1 //5 events, Dvz = 1cm, DC = 10
k5Cent,  //=2 //10 events, Dvz = 1cm, DC = 5
k5Evts5Cent
      };

AliRsnMiniAnalysisTask * AddTaskKStarMult13TeVpp_AOD
(
 Bool_t      isMC=kFALSE,
 Bool_t      isPP=kTRUE,
 TString     outNameSuffix="tofveto3stpc2s",
 Int_t       evtCutSetID=0,
 Int_t       pairCutSetID=0,
 Int_t       mixingConfigID=0,
 Int_t       aodFilterBit=5,
 Int_t       customQualityCutsID=1,
AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutKaCandidate=AliRsnCutSetDaughterParticle::kTPCpidTOFveto3s,
 Float_t     nsigmaPi       = 2.0,
 Float_t     nsigmaK        = 2.0,
 Float_t     nsigmaTOF       = 3.0,
 Bool_t      enableMonitor  = kTRUE
 )
{
    //-------------------------------------------
    // event cuts
    //-------------------------------------------
    bool isAOD=(evtCutSetID>=1000);
    UInt_t      triggerMask=AliVEvent::kINT7;
  if(evtCutSetID%1000>=100){
       triggerMask=AliVEvent::kHighMultV0;
     }
     evtCutSetID=evtCutSetID%100;
  Bool_t      rejectPileUp = kTRUE;
  Double_t    vtxZcut      = 10.0;//cm, default cut on vtx z
  Int_t       MultBins= aodFilterBit/100;

  if(evtCutSetID==eventCutSet::kDefaultVtx12) vtxZcut=12.0; //cm
  if(evtCutSetID==eventCutSet::kDefaultVtx8) vtxZcut=8.0; //cm
  if(evtCutSetID==eventCutSet::kDefaultVtx5) vtxZcut=5.0; //cm
  if(evtCutSetID==eventCutSet::kNoPileUpCut) rejectPileUp=kFALSE;
  if(evtCutSetID==eventCutSet::kNoVzCut) vtxZcut=1.e6;//off
  
  
  if(!isPP || isMC || MultBins) rejectPileUp=kFALSE;


  //-------------------------------------------
  //pair cuts
  //-------------------------------------------
  Double_t    minYlab = -0.5;
  Double_t    maxYlab = 0.5;
  
  if(pairCutSetID==pairYCutSet::kCentral){//|y_cm|<0.3
    minYlab=-0.3; maxYlab=0.3;
  }
    Bool_t CheckDecay=true;
     if(customQualityCutsID==99){customQualityCutsID=1; CheckDecay=false;}
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
  //
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddTaskKStarMult13TeVpp_AOD", "No analysis manager to connect to.");
      return NULL;
   }

   // create the task and configure
   
   TString taskName = Form("KStar%s%s", (isPP? "pp" : "PbPb"), (isMC ? "MC" : "Data"));
   
   AliRsnMiniAnalysisTask *task = new AliRsnMiniAnalysisTask(taskName.Data(), isMC);
   if(evtCutSetID!=eventCutSet::kNoEvtSel && evtCutSetID!=eventCutSet::kINEL10 && evtCutSetID!=eventCutSet::kIGZ10 && evtCutSetID!=eventCutSet::kIGZ){
     if(!isAOD)task->UseESDTriggerMask(triggerMask); //ESD
     task->SelectCollisionCandidates(triggerMask); //AOD
   }
   
   if(isPP){
     if(MultBins==1) task->UseMultiplicity("AliMultSelection_V0M");
     else if(MultBins==2) task->UseMultiplicity("AliMultSelection_SPDTracklets08");// C here
     else if(MultBins==3) task->UseMultiplicity("AliMultSelection_SPDTracklets08to15");
     else task->UseMultiplicity("QUALITY");
   }else task->UseCentrality("V0M");
    
   // set event mixing options
   task->UseContinuousMix();
   //task->UseBinnedMix();
   task->SetNMix(nmix);
   task->SetMaxDiffVz(maxDiffVzMix);
   task->SetMaxDiffMult(maxDiffMultMix);
   task->UseMC(isMC);
   ::Info("AddTaskKStarMult13TeVpp_AOD", Form("Event mixing configuration: \n events to mix = %i \n max diff. vtxZ = cm %5.3f \n max diff multi = %5.3f \n", nmix, maxDiffVzMix, maxDiffMultMix));
   
   mgr->AddTask(task);
   
   //
   // -- EVENT CUTS (same for all configs) ---------------------------------------------------------
   //
   // cut on primary vertex:
   // - 2nd argument --> |Vz| range
   // - 3rd argument --> minimum required number of contributors
   // - 4th argument --> tells if TPC stand-alone vertexes must be accepted
    if(!isPP || isMC) rejectPileUp=kFALSE;

   AliRsnCutPrimaryVertex* cutVertex=0;
   if(evtCutSetID!=eventCutSet::kTriggered && evtCutSetID!=eventCutSet::kNoEvtSel && evtCutSetID!=eventCutSet::kIGZ){
     if(evtCutSetID==eventCutSet::kINEL10 || evtCutSetID==eventCutSet::kIGZ10){
       cutVertex=new AliRsnCutPrimaryVertex("cutVertex",vtxZcut,0,kFALSE);
       cutVertex->SetCheckGeneratedVertexZ();
       
     }else if(!MultBins || fabs(vtxZcut-10.)>1.e-10){
       cutVertex=new AliRsnCutPrimaryVertex("cutVertex",vtxZcut,0,kFALSE);
       if(!MultBins){
     cutVertex->SetCheckZResolutionSPD();
     cutVertex->SetCheckDispersionSPD();
     cutVertex->SetCheckZDifferenceSPDTrack();
       }
     }
   }

  AliRsnCutEventUtils* cutEventUtils=0;
  if(evtCutSetID!=eventCutSet::kNoEvtSel && evtCutSetID!=eventCutSet::kINEL10){
    cutEventUtils=new AliRsnCutEventUtils("cutEventUtils",kTRUE,rejectPileUp);
    if(evtCutSetID==eventCutSet::kIGZ10 || evtCutSetID==eventCutSet::kIGZ) cutEventUtils->SetCheckInelGt0MC();
    else if(!MultBins){
      cutEventUtils->SetCheckIncompleteDAQ();
      cutEventUtils->SetCheckSPDClusterVsTrackletBG();
    }else{
      //cutEventUtils->SetCheckInelGt0SPDtracklets();
      cutEventUtils->SetRemovePileUppA2013(kFALSE);
      if(evtCutSetID!=eventCutSet::kTriggered) cutEventUtils->SetCheckAcceptedMultSelection();
    }
  }

  if(isPP && (!isMC) && cutVertex){
    cutVertex->SetCheckPileUp(rejectPileUp);// set the check for pileup
    ::Info("AddTaskKStarMult13TeVpp_AOD", "%s", Form(":::::::::::::::::: Pile-up rejection mode: %s", (rejectPileUp)?"ON":"OFF"));
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
  //
  // -- EVENT-ONLY COMPUTATIONS -------------------------------------------------------------------
  //
  //vertex
  Int_t vtxID = task->CreateValue(AliRsnMiniValue::kVz, kFALSE);
  AliRsnMiniOutput *outVtx = task->CreateOutput("eventVtx", "HIST", "EVENT");
  outVtx->AddAxis(vtxID, 400, -20.0, 20.0);
  
  //multiplicity or centrality
  Int_t multID=task->CreateValue(AliRsnMiniValue::kMult,kFALSE);
  AliRsnMiniOutput* outMult=task->CreateOutput("eventMult","HIST","EVENT");
  if(isPP && !MultBins) outMult->AddAxis(multID,400,0.5,400.5);
  else outMult->AddAxis(multID,110,0.,110.);

  Double_t multbins[200];
  int j,nmult=0;
  for(j=0;j<10;j++){multbins[nmult]=0.0001*j; nmult++;}
  for(j=1;j<10;j++){multbins[nmult]=0.001*j; nmult++;}
  for(j=1;j<50;j++){multbins[nmult]=0.01*j; nmult++;}
  for(j=5;j<10;j++){multbins[nmult]=0.1*j; nmult++;}
  for(j=1;j<=100;j++){multbins[nmult]=j; nmult++;}
  nmult--;
    TH1F* hEventsVsMulti=new TH1F("hAEventsVsMulti","",nmult,multbins);
    task->SetEventQAHist("EventsVsMulti",hEventsVsMulti);//custom binning for fHAEventsVsMulti
    TH2F* hvz=new TH2F("hVzVsCent","",110,0.,110., 400,-20.0,20.0);
     task->SetEventQAHist("vz",hvz);//plugs this histogram into the fHAEventVz data member
    
    double ybins[500];
    for(j=0;j<=401;j++) ybins[j]=j-0.5;
    //Here
    TH2F* hmc=new TH2F("MultiVsCent","", nmult,multbins, 401,ybins);
    hmc->GetYaxis()->SetTitle("QUALITY");
    task->SetEventQAHist("multicent",hmc);//plugs this histogram into the fHAEventMultiCent data member
    
  
  //
  // -- PAIR CUTS (common to all resonances) ------------------------------------------------------
  //
  
   AliRsnCutMiniPair* cutY=new AliRsnCutMiniPair("cutRapidity",AliRsnCutMiniPair::kRapidityRange);
   cutY->SetRangeD(minYlab,maxYlab);

   AliRsnCutSet* cutsPair=new AliRsnCutSet("pairCuts",AliRsnTarget::kMother);
   cutsPair->AddCut(cutY);
   cutsPair->SetCutScheme(cutY->GetName());

   task->SetCheckDecay(CheckDecay);

   //   gROOT->LoadMacro("ConfigKStarMult13TeVpp_AOD.C"); 
   // -- CONFIG ANALYSIS --------------------------------------------------------------------------
   //
   
  
       #if !defined (__CINT__) || defined (__CLING__)
   if (!ConfigKStarMult13TeVpp_AOD(task, isMC, isPP, "", cutsPair, aodFilterBit, customQualityCutsID, cutKaCandidate, nsigmaPi,nsigmaK,nsigmaTOF, enableMonitor,triggerMask)) return 0x0;
    #else
      //gROOT->LoadMacro("ConfigKStarMult13TeVpp_AOD.C");
gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/ConfigKStarMult13TeVpp_AOD.C");
      if (!ConfigKStarMult13TeVpp_AOD(task, isMC, isPP, "", cutsPair, aodFilterBit, customQualityCutsID, cutKaCandidate, nsigmaPi,nsigmaK,nsigmaTOF, enableMonitor,triggerMask)) return 0x0;
       #endif
       //
   // -- CONTAINERS --------------------------------------------------------------------------------
   //
   TString outputFileName = AliAnalysisManager::GetCommonFileName();
   //  outputFileName += ":Rsn";
   Printf("AddTaskKStarMult13TeVpp_AOD - Set OutputFileName : \n %s\n", outputFileName.Data() );
   
   AliAnalysisDataContainer *output = mgr->CreateContainer(Form("RsnOut_%s",outNameSuffix.Data()),
                               TList::Class(),
                               AliAnalysisManager::kOutputContainer,
                               outputFileName);
   mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(task, 1, output);
   
   return task;
}
