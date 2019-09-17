/***************************************************************************
Modified on: 04/08/2019
By: Pragati Sahoo

Modified by himani.bhatt@cern.ch  - last modified on 30/04/2018 
priyanka.sett@cern.ch - last modified on 27/10/2016

//Launches LStar analysis with rsn mini package
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

AliRsnMiniAnalysisTask * AddTaskLstar13TeVpp_II
(
 Bool_t      isMC=kFALSE,
 Bool_t      isPP=kTRUE,
 TString     outNameSuffix="tpc2stof3veto",
 Int_t       evtCutSetID=0,
 Int_t       pairCutSetID=0,
 Int_t       mixingConfigID=0,
 Int_t       aodFilterBit=5,
 Int_t       customQualityCutsID=1,
 Float_t     nsigmaPr = 2.0,
 Float_t     nsigmaKa=2.,
 Bool_t      enableMonitor=kTRUE,
 Bool_t      IsMcTrueOnly=kFALSE
 )
{  

AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutPrCandidate = AliRsnCutSetDaughterParticle::kTPCTOFpidLstar13ppTeV;
  AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutKaCandidate = AliRsnCutSetDaughterParticle::kTPCTOFpidLstar13ppTeV;
  
  //-------------------------------------------
  // event cuts
  //-------------------------------------------
  UInt_t      triggerMask=AliVEvent::kINT7;
  if(evtCutSetID>=100){
    triggerMask=AliVEvent::kHighMultV0;
    evtCutSetID=evtCutSetID%100;
  }

  Bool_t      rejectPileUp=kTRUE;
  Double_t    vtxZcut=10.0;//cm, default cut on vtx z
  Int_t       MultBins=aodFilterBit/100;

  if(evtCutSetID==eventCutSet::kDefaultVtx12) vtxZcut=12.0; //cm
  if(evtCutSetID==eventCutSet::kDefaultVtx8) vtxZcut=8.0; //cm
  if(evtCutSetID==eventCutSet::kDefaultVtx5) vtxZcut=5.0; //cm
  if(evtCutSetID==eventCutSet::kNoPileUpCut) rejectPileUp=kFALSE;
  if(evtCutSetID==eventCutSet::kNoVzCut) vtxZcut=1.e6;//off

  if(!isPP || isMC || MultBins) rejectPileUp=kFALSE;

  //-------------------------------------------
  //pair cuts
  //-------------------------------------------
  Double_t    minYlab=-0.5;
  Double_t    maxYlab= 0.5;

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

  AliAnalysisManager* mgr=AliAnalysisManager::GetAnalysisManager();
  if(!mgr){
    ::Error("AddTaskLstar13TeVpp_II", "No analysis manager to connect to.");
    return NULL;
  }

  // create the task and configure
  TString taskName = Form("LStar%s%s_%i",(isPP? "pp" : "PbPb"),(isMC ? "MC" : "Data"));
    //Form("%s", AliAnalysisManager::GetCommonFileName());
    //Form("LStar%s%s_%i",(isPP? "pp" : "PbPb"),(isMC ? "MC" : "Data"),(Int_t)cutKaCandidate);
 
  // Objects name
  AliRsnMiniAnalysisTask* task=new AliRsnMiniAnalysisTask(taskName.Data(),isMC);
  
  if(evtCutSetID!=eventCutSet::kNoEvtSel && evtCutSetID!=eventCutSet::kINEL10 && evtCutSetID!=eventCutSet::kIGZ10 && evtCutSetID!=eventCutSet::kIGZ){
    task->UseESDTriggerMask(triggerMask); //ESD
    //task->SelectCollisionCandidates(triggerMask); //AOD
  }

  if(isPP){
    if(MultBins==1) task->UseMultiplicity("AliMultSelection_V0M");
    else if(MultBins==2) task->UseMultiplicity("AliMultSelection_RefMult08");
    else if(MultBins==3) task->UseMultiplicity("AliMultSelection_SPDTracklets08to15");
    else task->UseMultiplicity("QUALITY");
  }else task->UseCentrality("V0M");

  // set event mixing options
  task->UseContinuousMix();
  //task->UseBinnedMix();
  task->SetNMix(nmix);
  task->SetMaxDiffVz(maxDiffVzMix);
  task->SetMaxDiffMult(maxDiffMultMix);
  ::Info("AddTaskLstar13TeVpp_II", Form("Event mixing configuration: \n events to mix = %i \n max diff. vtxZ = cm %5.3f \n max diff multi = %5.3f", nmix, maxDiffVzMix, maxDiffMultMix));
  //task->SaveRsnTreeInFile(kTRUE);

  mgr->AddTask(task);

  // cut on primary vertex:
  // - 2nd argument --> |Vz| range
  // - 3rd argument --> minimum required number of contributors to vtx
  // - 4th argument --> tells if TPC stand-alone vertexes must be accepted

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
    ::Info("AddTaskLstar13TeVpp_II", Form(":::::::::::::::::: Pile-up rejection mode: %s", (rejectPileUp)?"ON":"OFF"));
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

  // -- EVENT-ONLY COMPUTATIONS -------------------------------------------------------------------

  //vertex
  Int_t vtxID=task->CreateValue(AliRsnMiniValue::kVz,kFALSE);
  AliRsnMiniOutput* outVtx=task->CreateOutput("eventVtx","HIST","EVENT");
  outVtx->AddAxis(vtxID,240,-12.0,12.0);

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

  TH2F* hvz=new TH2F("hVzVsCent","",110,0.,110., 240,-12.0,12.0);
  task->SetEventQAHist("vz",hvz);//plugs this histogram into the fHAEventVz data member

  double ybins[500];
  for(j=0;j<=401;j++) ybins[j]=j-0.5;

  TH2F* hmc=new TH2F("MultiVsCent","", nmult,multbins, 401,ybins);
  hmc->GetYaxis()->SetTitle("QUALITY");
  task->SetEventQAHist("multicent",hmc);//plugs this histogram into the fHAEventMultiCent data member

  // -- PAIR CUTS (common to all resonances) ------------------------------------------------------

  AliRsnCutMiniPair* cutY=new AliRsnCutMiniPair("cutRapidity",AliRsnCutMiniPair::kRapidityRange);
  cutY->SetRangeD(minYlab,maxYlab);

  AliRsnCutSet* cutsPair=new AliRsnCutSet("pairCuts",AliRsnTarget::kMother);
  cutsPair->AddCut(cutY);
  cutsPair->SetCutScheme(cutY->GetName());

  task->SetCheckDecay(CheckDecay);

  // -- CONFIG ANALYSIS --------------------------------------------------------------------------

  gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/ConfigureLstar13TeVpp_II.C");
  if (!ConfigureLstar13TeVpp_II(task, isMC, isPP, "", cutsPair, aodFilterBit, customQualityCutsID, cutPrCandidate, cutKaCandidate, nsigmaPr, nsigmaKa,  enableMonitor, isMC&IsMcTrueOnly)) return 0x0;
 

  // -- CONTAINERS --------------------------------------------------------------------------------

  TString outputFileName = AliAnalysisManager::GetCommonFileName();
  //  outputFileName += ":Rsn";
   Printf("AddAnalysisTaskLstar13TeVpp - Set OutputFileName : \n %s\n", outputFileName.Data() );
 
  AliAnalysisDataContainer* output=mgr->CreateContainer(Form("RsnOut_%s",outNameSuffix.Data()),TList::Class(),AliAnalysisManager::kOutputContainer, outputFileName);
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, output);

  return task;
}
