/******************
Modify the macro for ROOT6 enviroment   
****************************************************************************/




#ifdef __CLING__
R__ADD_INCLUDE_PATH($ALICE_PHYSICS)
#include <PWGLF/RESONANCES/macros/mini/AddMonitorOutput.C>
#include <PWGLF/RESONANCES/macros/mini/qa/AddMonitorOutputV0.C>
#include <PWGLF/RESONANCES/macros/mini/qa/AddMonitorOutputCascade.C>
#endif

#include "TDatabasePDG.h"
#include "TParticlePDG.h"



Bool_t ConfigRsnQA_ROOT6(AliRsnMiniAnalysisTask *task, 
			    Bool_t                 isMC, 
			    Bool_t                 isPP,
			    AliRsnCutSet           *cutsPair,
			    Int_t                  Strcut = 2011,
			    Int_t                  customQualityCutsID = AliRsnCutSetDaughterParticle::kDisableCustom,
			    AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutKaCandidate=AliRsnCutSetDaughterParticle::kTPCTOFpidphikstarpPb2016,
			    Float_t                nsigmaPi = 2.0,
			    Float_t                nsigmaK  = 2.0,
			    Float_t                nsigmaTOF= 3.0,
			    Bool_t                 enableMonitor = kTRUE
			    );
Bool_t SetCustomQualityCut(AliRsnCutTrackQuality * trkQualityCut, 
                           Int_t customQualityCutsID = 1,
                           Int_t trCut = 2011);


AliRsnMiniAnalysisTask * AddTaskRsnQA_ROOT6(
						Bool_t      isMC                = kFALSE,
						Bool_t      isPP                = kFALSE,
						Int_t       Strcut              = 2011,
						Int_t       customQualityCutsID = AliRsnCutSetDaughterParticle::kDisableCustom,		     
            AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutKaCandidate=AliRsnCutSetDaughterParticle::kTPCTOFpidphikstarpPb2016,

						Float_t     nsigmaPi            = 2.0,
						Float_t     nsigmaK             = 2.0,
						Float_t     nsigmaTOF           = 3.0,
						Bool_t      enableMonitor       = kTRUE,
						Int_t       nmix                = 5,
						Float_t     maxDiffVzMix        = 1.0,
						Float_t     maxDiffMultMix      = 5.0,
						TString     outNameSuffix       = "pPb"
						)
{  
  Bool_t      rejectPileUp = kTRUE;




  //
  // -- INITIALIZATION ----------------------------------------------------------------------------
  // retrieve analysis manager
  //
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddTaskPhipp", "No analysis manager to connect to.");
      return NULL;
   } 


   // create the task and configure 
   TString taskName = Form("Phi%s%s", (isPP? "pp" : "PPb"), (isMC ? "MC" : "Data"));
   
   AliRsnMiniAnalysisTask *task = new AliRsnMiniAnalysisTask(taskName.Data(), isMC);
   task->UseESDTriggerMask(AliVEvent::kINT7);
   if (isPP) 
   task->UseMultiplicity("QUALITY");
   else
   task->UseMultiplicity("AliMultSelection_V0M");//Only for RunII
   // set event mixing options
   task->UseContinuousMix();
   //task->UseBinnedMix();
   task->SetNMix(nmix);
   task->SetMaxDiffVz(maxDiffVzMix);
   task->SetMaxDiffMult(maxDiffMultMix);
   task->UseMC(isMC);
   ::Info("AddTaskPhiRsnQA", "%s", Form("Event mixing configuration: \n events to mix = %i \n max diff. vtxZ = cm %5.3f \n max diff multi = %5.3f \n", nmix, maxDiffVzMix, maxDiffMultMix));
   
   mgr->AddTask(task);
   
   AliRsnCutEventUtils* cutEventUtils=new AliRsnCutEventUtils("cutEventUtils",kTRUE,rejectPileUp);
   cutEventUtils->SetCheckAcceptedMultSelection();
   AliRsnCutSet *eventCuts = new AliRsnCutSet("eventCuts", AliRsnTarget::kEvent);
   eventCuts->AddCut(cutEventUtils);
   eventCuts->SetCutScheme(Form("%s",cutEventUtils->GetName()));
   task->SetEventCuts(eventCuts);
   
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
   if (isPP) 
     outMult->AddAxis(multID, 400, 0.0, 400.0);
   else
     outMult->AddAxis(multID, 100, 0.0, 100.0);
   
   //
   // -- PAIR CUTS (common to all resonances) ------------------------------------------------------
   //
   AliRsnCutMiniPair *cutY = new AliRsnCutMiniPair("cutRapidity", AliRsnCutMiniPair::kRapidityRange);
   cutY->SetRangeD(-0.5, 0.5);
   AliRsnCutSet *cutsPair = new AliRsnCutSet("pairCuts", AliRsnTarget::kMother);
   cutsPair->AddCut(cutY);
   cutsPair->SetCutScheme(cutY->GetName());
   

   //
   // -- CONFIG ANALYSIS --------------------------------------------------------------------------
   //
   
     if (!ConfigRsnQA_ROOT6(task, isMC, isPP, cutsPair,Strcut,customQualityCutsID,cutKaCandidate,nsigmaPi,nsigmaK,nsigmaTOF,enableMonitor))
       {
	 printf(" returning 0 from Config");
	 return 0x0;
   
       }   
   // -- CONTAINERS --------------------------------------------------------------------------------
   //

   TString outputFileName = AliAnalysisManager::GetCommonFileName();
   //outputFileName += ":Rsn";
   Printf("AddTaskPhiRsnQA - Set OutputFileName : \n %s\n", outputFileName.Data() );
   
   AliAnalysisDataContainer *output = mgr->CreateContainer(Form("RsnOut_%s",outNameSuffix.Data()), 
							   TList::Class(), 
							   AliAnalysisManager::kOutputContainer, 
							   outputFileName);
   mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(task, 1, output);
   
   return task;
}


Bool_t ConfigRsnQA_ROOT6(AliRsnMiniAnalysisTask *task, Bool_t isMC, Bool_t isPP, AliRsnCutSet *cutsPair, Int_t Strcut, Int_t customQualityCutsID, AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutKaCandidate, Float_t nsigmaPi, Float_t nsigmaK, Float_t nsigmaTOF, Bool_t enableMonitor)
{
  // retrieve mass from PDG database
  Int_t         pdg  = 333;
  TDatabasePDG *db   = NULL;
  db=TDatabasePDG::Instance();
  TParticlePDG *part = NULL;
  part=db->GetParticle(pdg);
  Double_t mass      = part->Mass();

  // set daughter cuts
  AliRsnCutSetDaughterParticle* cutSetPi=NULL;
  AliRsnCutSetDaughterParticle* cutSetK=NULL;
  AliRsnCutTrackQuality* trkQualityCut= new AliRsnCutTrackQuality("myQualityCut");
  if(!trkQualityCut) return kFALSE;


     // if(SetCustomQualityCut(trkQualityCut,customQualityCutsID,Strcut)){
    cutSetPi=new AliRsnCutSetDaughterParticle(Form("cutPi%i_%2.1fsigmaTPC_%2.1fsigmaTOF",cutKaCandidate,nsigmaPi,nsigmaTOF),trkQualityCut,cutKaCandidate,AliPID::kPion,nsigmaPi,nsigmaTOF);
    cutSetK=new AliRsnCutSetDaughterParticle(Form("cutK%i_%2.1fsigma_%2.1fsigmaTOF",cutKaCandidate, nsigmaK,nsigmaTOF),trkQualityCut,cutKaCandidate,AliPID::kKaon,nsigmaK,nsigmaTOF);
    //}
    /*else{
    printf("Daughter Track cuts has been selected =================\n");
    return kFALSE;
    }*/

  Int_t iCutPi = task->AddTrackCuts(cutSetPi);
  Int_t iCutK  = task->AddTrackCuts(cutSetK);
  if(enableMonitor){
    Printf("======== Monitoring cut AliRsnCutSetDaughterParticle enabled");

    AddMonitorOutput(isMC, cutSetPi->GetMonitorOutput());
    AddMonitorOutput(isMC, cutSetK->GetMonitorOutput());
  }

  // -- Values ------------------------------------------------------------------------------------
  

/* invariant mass   */ Int_t imID   = task->CreateValue(AliRsnMiniValue::kInvMass,    kFALSE);
/* IM resolution    */ Int_t resID  = task->CreateValue(AliRsnMiniValue::kInvMassRes, kTRUE);
/* transv. momentum */  Int_t ptID   = task->CreateValue(AliRsnMiniValue::kPt,         kFALSE);
/* centrality       */  Int_t centID = task->CreateValue(AliRsnMiniValue::kMult,       kFALSE);
/* pseudorapidity   */  Int_t etaID  = task->CreateValue(AliRsnMiniValue::kEta,        kFALSE);
/* rapidity         */  Int_t yID    = task->CreateValue(AliRsnMiniValue::kY,          kFALSE);






  Bool_t  use     [12] = {1         ,1         ,1         ,1         ,1       ,1       ,isMC     ,isMC     ,isMC     ,isMC     ,isMC    ,isMC    };
  Bool_t  useIM   [12] = {1         ,1         ,1         ,1         ,1       ,1       ,1        ,1        ,1        ,1        ,0       ,0       };
  TString name    [12] = {"UnlikePM","UnlikeMP","MixingPM","MixingMP","LikePP","LikeMM","MCGenPM","MCGenMP","TruesPM","TruesMP","ResPM" ,"ResMP" };
  TString comp    [12] = {"PAIR"    ,"PAIR"    ,"MIX"     ,"MIX"     ,"PAIR"  ,"PAIR"  ,"MOTHER" ,"MOTHER" ,"TRUE"   ,"TRUE"   ,"TRUE"  ,"TRUE"  };
  TString output  [12] = {"SPARSE"  ,"SPARSE"  ,"SPARSE"  ,"SPARSE"  ,"SPARSE","SPARSE","SPARSE" ,"SPARSE" ,"SPARSE" ,"SPARSE" ,"SPARSE","SPARSE"};
  Char_t  charge1 [12] = {'+'       ,'-'       ,'+'       ,'-'       ,'+'     ,'-'     ,'+'      ,'-'      ,'+'      ,'-'      ,'+'     ,'-'     };
  Char_t  charge2 [12] = {'-'       ,'+'       ,'-'       ,'+'       ,'+'     ,'-'     ,'-'      ,'+'      ,'_'      ,'+'      ,'-'     ,'+'     };
  Int_t   cutIDK  [12] = {iCutK     ,iCutK     ,iCutK     ,iCutK     ,iCutK   ,iCutK   ,iCutK    ,iCutK    ,iCutK    ,iCutK    ,iCutK   ,iCutK   };
  Int_t   cutIDPi [12] = {iCutPi    ,iCutPi    ,iCutPi    ,iCutPi    ,iCutPi  ,iCutPi  ,iCutPi   ,iCutPi   ,iCutPi   ,iCutPi   ,iCutPi  ,iCutPi  };
  Int_t   PDGCode [12] = {333       ,-333      ,333       ,333       ,333     ,333     ,333      ,-333     ,333      ,-333      ,333     ,-333   };

  for (Int_t i = 0; i < 12; i++) {
    if (!use[i]) continue;
    AliRsnMiniOutput *out = task->CreateOutput(Form("CustomId%d_%s", customQualityCutsID, name[i].Data()), output[i].Data(), comp[i].Data());
    out->SetDaughter(0, AliRsnDaughter::kKaon);
    out->SetDaughter(1, AliRsnDaughter::kKaon);
    out->SetCutID(0, cutIDK[i]);
    out->SetCutID(1, cutIDK[i]);
    out->SetCharge(0, charge1[i]);
    out->SetCharge(1, charge2[i]);
    out->SetMotherPDG(PDGCode[i]);
    out->SetMotherMass(mass);
    out->SetPairCuts(cutsPair);

    // axis X: invmass (or resolution)
    if (useIM[i]) 
      out->AddAxis(imID, 180, 0.9, 1.8);
     else
     out->AddAxis(resID, 200, -0.02, 0.02);
    // axis Y: transverse momentum
    out->AddAxis(ptID, 500, 0.0, 50.0);
    // axis Z: centrality-multiplicity
    //if (!isPP)
    out->AddAxis(centID, 100, 0.0, 100.0);
    //else 
    //out->AddAxis(centID, 400, 0.0, 400.0);
    //axis W: pseudorapidity
    //out->AddAxis(etaID, 20, -1.0, 1.0);
    //axis J: rapidity
    //out->AddAxis(yID, 90, -4.5, 4.5);
  }
   return kTRUE;
  //  return kFALSE;
}
Bool_t SetCustomQualityCut(AliRsnCutTrackQuality * trkQualityCut, Int_t customQualityCutsID,Int_t trCut)
{
  //Sets configuration for track quality object different from std quality cuts.
  //Returns kTRUE if track quality cut object is successfully defined,
  //returns kFALSE if an invalid set of cuts (customQualityCutsID) is chosen or if the
  //object to be configured does not exist.

  if ((!trkQualityCut)){
    Printf("::::: SetCustomQualityCut:: use default quality cuts specified in task configuration.");
    return kFALSE;
  }

  if(customQualityCutsID>=1 && customQualityCutsID<100 && customQualityCutsID!=2){
    if(trCut == 2011){
      trkQualityCut->SetDefaults2011(kTRUE,kTRUE);
      Printf("%s", Form("::::: SetCustomQualityCut:: using standard 2011 track quality cuts"));
    }
    else if(trCut == 2015){
      trkQualityCut->SetDefaults2011(kTRUE,kTRUE);
      trkQualityCut->GetESDtrackCuts()->SetCutGeoNcrNcl(3., 130., 1.5, 0.85, 0.7);
      Printf("%s", Form("::::: SetCustomQualityCut:: using standard 2015 track quality cuts"));
    }
    if(customQualityCutsID==3){trkQualityCut->GetESDtrackCuts()->SetMaxDCAToVertexXYPtDep("0.015+0.05/pt^1.1");}//10Sig // D = 7*(0.0015+0.0050/pt^1.1)
    else if(customQualityCutsID==4){trkQualityCut->GetESDtrackCuts()->SetMaxDCAToVertexXYPtDep("0.006+0.02/pt^1.1");}//4Sig
    else if(customQualityCutsID==5){trkQualityCut->GetESDtrackCuts()->SetMaxDCAToVertexZ(3.);}// D = 2.
    else if(customQualityCutsID==6){trkQualityCut->GetESDtrackCuts()->SetMaxDCAToVertexZ(1.);} 
    else if(customQualityCutsID==7){trkQualityCut->GetESDtrackCuts()->SetMaxDCAToVertexZ(0.2);}
    else if(customQualityCutsID==8){trkQualityCut->GetESDtrackCuts()->SetMaxChi2PerClusterTPC(5.);}// D = 4
    else if(customQualityCutsID==9){trkQualityCut->GetESDtrackCuts()->SetMaxChi2PerClusterTPC(3);}
    else if(customQualityCutsID==10){trkQualityCut->GetESDtrackCuts()->SetMinNCrossedRowsTPC(60);}// D = 70
    else if(customQualityCutsID==11){trkQualityCut->GetESDtrackCuts()->SetMinNCrossedRowsTPC(80);}
    else if(customQualityCutsID==12){trkQualityCut->GetESDtrackCuts()->SetMinNCrossedRowsTPC(100);}
    else if(customQualityCutsID==13){trkQualityCut->GetESDtrackCuts()->SetMinRatioCrossedRowsOverFindableClustersTPC(0.7);}// D = 8
    else if(customQualityCutsID==14){trkQualityCut->GetESDtrackCuts()->SetMinRatioCrossedRowsOverFindableClustersTPC(0.9);}
    else if(customQualityCutsID==15){trkQualityCut->GetESDtrackCuts()->SetMaxChi2PerClusterITS(49.);}// D = 36
    else if(customQualityCutsID==16){trkQualityCut->GetESDtrackCuts()->SetMaxChi2PerClusterITS(25.);}
    else if(customQualityCutsID==17){trkQualityCut->GetESDtrackCuts()->SetMaxChi2TPCConstrainedGlobal(49.);}// D = 36
    else if(customQualityCutsID==18){trkQualityCut->GetESDtrackCuts()->SetMaxChi2TPCConstrainedGlobal(25.);}
    else if(customQualityCutsID==19){trkQualityCut->GetESDtrackCuts()->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kOff);}
    
    trkQualityCut->Print();
    return kTRUE;
  }else if(customQualityCutsID==2 || (customQualityCutsID>=100 && customQualityCutsID<200)){
    trkQualityCut->SetDefaultsTPCOnly(kTRUE);
    Printf("%s", Form("::::: SetCustomQualityCut:: using TPC-only track quality cuts"));
    if(customQualityCutsID==103){trkQualityCut->GetESDtrackCuts()->SetMaxDCAToVertexXY(3.);}
    else if(customQualityCutsID==104){trkQualityCut->GetESDtrackCuts()->SetMaxDCAToVertexXY(1.);}
    else if(customQualityCutsID==105){trkQualityCut->GetESDtrackCuts()->SetMaxDCAToVertexZ(4.);}
    else if(customQualityCutsID==106){trkQualityCut->GetESDtrackCuts()->SetMaxDCAToVertexZ(1.);}
    else if(customQualityCutsID==107){trkQualityCut->GetESDtrackCuts()->SetMaxChi2PerClusterTPC(7.);}
    else if(customQualityCutsID==108){trkQualityCut->GetESDtrackCuts()->SetMaxChi2PerClusterTPC(2.5);}
    else if(customQualityCutsID==109){trkQualityCut->GetESDtrackCuts()->SetMinNClustersTPC(30);}
    else if(customQualityCutsID==110){trkQualityCut->GetESDtrackCuts()->SetMinNClustersTPC(85);}

    trkQualityCut->Print();
    return kTRUE;
  }else{
    Printf("::::: SetCustomQualityCut:: use default quality cuts specified in task configuration.");
    return kFALSE;
  }
  trkQualityCut->SetPtRange(0.15, 100000.0);
  trkQualityCut->SetEtaRange(-0.8, 0.8);
  
  Printf("%s", Form("::::: SetCustomQualityCut:: using custom track quality cuts #%i",customQualityCutsID));
  trkQualityCut->Print();
  return kTRUE;
}
