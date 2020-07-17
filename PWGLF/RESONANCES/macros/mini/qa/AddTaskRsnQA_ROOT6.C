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
#include "AliRsnCutCascade.h"





void AddMonitorOutput_P(TString s="",TObjArray* m=0,TString o="",AliRsnLoopDaughter* l=0);
void AddMonitorOutput_Pt(TString s="",TObjArray* m=0,TString o="",AliRsnLoopDaughter* l=0);
void AddMonitorOutput_Eta(TString s="",TObjArray* m=0,TString o="",AliRsnLoopDaughter* l=0);
void AddMonitorOutput_DCAxy(TString s="",TObjArray* m=0,TString o="",AliRsnLoopDaughter* l=0);
void AddMonitorOutput_DCAz(TString s="",TObjArray* m=0,TString o="",AliRsnLoopDaughter* l=0);
void AddMonitorOutput_chi2TPC(TString s="",TObjArray* m=0,TString o="",AliRsnLoopDaughter* l=0);
void AddMonitorOutput_NclTPC(TString s="",TObjArray* m=0,TString o="",AliRsnLoopDaughter* l=0);
void AddMonitorOutput_V0Mass(TString s="",TObjArray* m=0,TString o="",AliRsnLoopDaughter* l=0);
void AddMonitorOutput_V0DCA(TString n="",TObjArray* m=0,TString o="",AliRsnLoopDaughter* l=0);
void AddMonitorOutput_V0Radius(TString s="",TObjArray* m=0,TString o="",AliRsnLoopDaughter* l=0);
void AddMonitorOutput_V0Lifetime(TString s="",TObjArray* m=0,TString o="",AliRsnLoopDaughter* l=0);
void AddMonitorOutput_V0DaughterDCA(TString s="",TObjArray* m=0,TString o="",AliRsnLoopDaughter* l=0);
void AddMonitorOutput_V0DCA2TPV(TString s="",TObjArray* m=0,TString o="",AliRsnLoopDaughter* l=0);










Bool_t ConfigRsnXiQA_ROOT6(AliRsnMiniAnalysisTask *task, 
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




Bool_t ConfigRsnChrgKStarQA_ROOT6(AliRsnMiniAnalysisTask *task, 
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





Bool_t ConfigRsnPhiQA_ROOT6(AliRsnMiniAnalysisTask *task, 
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
						TString     outNameSuffix       = "pp",
						Int_t       pidofdaugh1         =310,
						Int_t       pidofdaugh2         =211

					    )
{  //start of add task
  Bool_t      rejectPileUp = kTRUE;




  //
  // -- INITIALIZATION ----------------------------------------------------------------------------
  // retrieve analysis manager
  //
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddTaskRsnQA", "No analysis manager to connect to.");
      return NULL;
   } 


   // create the task and configure 
   TString taskName = Form("PhiKSTARXi_%s%s", (isPP? "pp" : "PPb"), (isMC ? "MC" : "Data"));
   
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
   ::Info("AddTaskRsnQA", "%s", Form("Event mixing configuration: \n events to mix = %i \n max diff. vtxZ = cm %5.3f \n max diff multi = %5.3f \n", nmix, maxDiffVzMix, maxDiffMultMix));
   
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
   if(pidofdaugh1==321 && pidofdaugh2==321)   //daughters of phi, daugh1=Kaon, daugh2=Kaon
	{
	  printf("------------------i am in phi configuration-------------------");
     if (!ConfigRsnPhiQA_ROOT6(task, isMC, isPP, cutsPair,Strcut,customQualityCutsID,cutKaCandidate,nsigmaPi,nsigmaK,nsigmaTOF,enableMonitor))
       {
	 printf("-------------------- returning 0 from phi Config-------------------");
	 return 0x0;
   
       }   
	}


      else if(pidofdaugh1==310 && pidofdaugh2==211)  //daughters of K*+/-, daugh1=KShort, daugh2=pion
	{
	  printf("------------------i am in KStar configuration-------------------");
   if (!ConfigRsnChrgKStarQA_ROOT6(task, isMC, isPP, cutsPair,Strcut,customQualityCutsID,cutKaCandidate,nsigmaPi,nsigmaK,nsigmaTOF,enableMonitor))
       {
	 printf("------------------returning 0 from chargekstar Config-------------------");
	 return 0x0;
   
       }   
	}


    else if(pidofdaugh1==3312 && pidofdaugh2==211)  //daughters of Xi*0, daugh1=Xi-, daugh2=pion+
	{
	  printf("------------------i am in Xi configuration-------------------");
   if (!ConfigRsnXiQA_ROOT6(task, isMC, isPP, cutsPair,Strcut,customQualityCutsID,cutKaCandidate,nsigmaPi,nsigmaK,nsigmaTOF,enableMonitor))
       {
	 printf("------------------returning 0 from chargekstar Config-------------------");
	 return 0x0;
   
       }   
	}









   // -- CONTAINERS --------------------------------------------------------------------------------
   //

   TString outputFileName = AliAnalysisManager::GetCommonFileName();
   //outputFileName += ":Rsn";
   Printf("AddTaskPhiPbPbRunTwo - Set OutputFileName : \n %s\n", outputFileName.Data() );
   
   AliAnalysisDataContainer *output = mgr->CreateContainer(Form("RsnOut_%s",outNameSuffix.Data()), 
							   TList::Class(), 
							   AliAnalysisManager::kOutputContainer, 
							   outputFileName);
   mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(task, 1, output);
   
   return task;
}




Bool_t ConfigRsnPhiQA_ROOT6(AliRsnMiniAnalysisTask *task, Bool_t isMC, Bool_t isPP, AliRsnCutSet *cutsPair, Int_t Strcut, Int_t customQualityCutsID, AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutKaCandidate, Float_t nsigmaPi, Float_t nsigmaK, Float_t nsigmaTOF, Bool_t enableMonitor)
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
  }
   return kTRUE;
  //  return kFALSE;
}




Bool_t ConfigRsnChrgKStarQA_ROOT6(AliRsnMiniAnalysisTask *task, Bool_t isMC, Bool_t isPP, AliRsnCutSet *cutsPair, Int_t Strcut, Int_t customQualityCutsID, AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutKaCandidate, Float_t nsigmaPi, Float_t nsigmaK, Float_t nsigmaTOF, Bool_t enableMonitor)

{


  AliRsnCutSetDaughterParticle* cutSetPi=NULL;
  AliRsnCutSetDaughterParticle* cutSetK=NULL;
  AliRsnCutTrackQuality* trkQualityCut= new AliRsnCutTrackQuality("myQualityCut");
  if(!trkQualityCut) return kFALSE;


  // if(SetCustomQualityCut(trkQualityCut,customQualityCutsID,Strcut)){                                                                     
  cutSetPi=new AliRsnCutSetDaughterParticle(Form("cutPi%i_%2.1fsigmaTPC_%2.1fsigmaTOF",cutKaCandidate,nsigmaPi,nsigmaTOF),trkQualityCut,cutKaCandidate,AliPID::kPion,nsigmaPi,nsigmaTOF);
  cutSetK=new AliRsnCutSetDaughterParticle(Form("cutK%i_%2.1fsigma_%2.1fsigmaTOF",cutKaCandidate, nsigmaK,nsigmaTOF),trkQualityCut,cutKaCandidate,AliPID::kKaon,nsigmaK,nsigmaTOF);
  //}

  Float_t aodFilterBit=5.0;

  // cutSetK=new AliRsnCutSetDaughterParticle(Form("cutK_bit%2.1fsigma",aodFilterBit),trkQualityCut,AliRsnCutSetDaughterParticle::kQualityStd2011,AliPID::kPion,-1.);
  //cutSetPi=new AliRsnCutSetDaughterParticle(Form("cutPi%i_%2.1fsigma_%2.1fsigma",cutKaCandidate, nsigmaPi,nsigmaTOF),trkQualityCut,cutKaCandidate,AliPID::kPion,nsigmaPi,nsigmaTOF);

                                                                                                                                        
  
  Int_t iCutPi = task->AddTrackCuts(cutSetPi);
  Int_t iCutK  = task->AddTrackCuts(cutSetK);



  ////////////////////////////////////////////////////////////                                                                               
  // selections for K0s and for the daughters of K0s                                                                                         
  /////////////////////////////////////////////////////////////                                                                              
  // 
  Float_t crossedRows=70;
  Float_t rowsbycluster=0.8;
  Float_t DCAxy=0.06;

                                                                                                                                        
  // selections for pion daugthers of K0s                                                                                                    
  AliESDtrackCuts *esdTrackCuts = new AliESDtrackCuts("qualityDaughterK0s");
  esdTrackCuts->SetEtaRange(-0.8,0.8);
  esdTrackCuts->SetRequireTPCRefit();
  esdTrackCuts->SetAcceptKinkDaughters(0); //                                                                                                
  esdTrackCuts->SetMinNCrossedRowsTPC(crossedRows);
  esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(rowsbycluster);
  esdTrackCuts->SetMaxChi2PerClusterTPC(100);
  esdTrackCuts->SetMinDCAToVertexXY(DCAxy); //Use one of the two - pt dependent or fixed value cut.                                          

  //                                                                                                                                         
  /////////////////////////////////////////////////                                                                                          
  // selections for K0s                      
  Float_t pi_k0s_PIDCut=5.0;
  Float_t vorapidity=0.5;
  Float_t k0sDaughDCA=0.3;
  Float_t k0sDCA=0.3;
  Float_t k0sCosPoinAn=0.98;
  Float_t massTol=0.03;
  Float_t massTolVeto=0.0043;
  Bool_t Switch=kTRUE;
  Float_t pLife=20;
  Float_t radiuslow=0.5;
  Float_t radiushigh=200;
  Float_t v0rapidity=0.5;
  Int_t   tol_switch = 1;
  Double_t tol_sigma = 6.0;



                                                                                                
  AliRsnCutV0 *cutK0s = new AliRsnCutV0("cutK0s", kK0Short, AliPID::kPion, AliPID::kPion);
  cutK0s->SetPIDCutPion(pi_k0s_PIDCut);        // PID for the pion daughter of K0s                                                           
  cutK0s->SetESDtrackCuts(esdTrackCuts);
  cutK0s->SetMaxDaughtersDCA(k0sDaughDCA);
  cutK0s->SetMaxDCAVertex(k0sDCA);
  cutK0s->SetMinCosPointingAngle(k0sCosPoinAn);
  cutK0s->SetTolerance(massTol);
  cutK0s->SetToleranceVeto(massTolVeto);   //Rejection range for Competing V0 Rejection                                                      
  cutK0s->SetSwitch(Switch);
  cutK0s->SetfLife(pLife);
  cutK0s->SetfLowRadius(radiuslow);
  cutK0s->SetfHighRadius(radiushigh);
  cutK0s->SetMaxRapidity(v0rapidity);
  cutK0s->SetpT_Tolerance(tol_switch);
  cutK0s->SetMassTolSigma(tol_sigma);

 
  AliRsnCutSet *cutSetK0s = new AliRsnCutSet("setK0s", AliRsnTarget::kDaughter);
  cutSetK0s->AddCut(cutK0s);
  cutSetK0s->SetCutScheme(cutK0s->GetName());
  Int_t iCutK0s = task->AddTrackCuts(cutSetK0s);


  TString pname="k0s";
  if(enableMonitor){
    Printf("======== Monitoring cut AliRsnCutSetDaughterParticle enabled");

    AddMonitorOutput(isMC, cutSetPi->GetMonitorOutput());
    AddMonitorOutput(isMC, cutSetK->GetMonitorOutput());
  

    AddMonitorOutput_P(pname,cutSetK0s->GetMonitorOutput());
    AddMonitorOutput_Pt(pname,cutSetK0s->GetMonitorOutput());
    AddMonitorOutput_V0Mass(pname,cutSetK0s->GetMonitorOutput());
    AddMonitorOutput_V0DCA(pname,cutSetK0s->GetMonitorOutput());
    AddMonitorOutput_V0Radius(pname,cutSetK0s->GetMonitorOutput());
    AddMonitorOutput_V0Lifetime(pname,cutSetK0s->GetMonitorOutput());
    AddMonitorOutput_V0DaughterDCA(pname,cutSetK0s->GetMonitorOutput());
    AddMonitorOutput_V0DCA2TPV(pname,cutSetK0s->GetMonitorOutput());


}



  AliRsnCutMiniPair* cutYKS=new AliRsnCutMiniPair("cutRapidity", AliRsnCutMiniPair::kRapidityRange);
  cutYKS->SetRangeD(-0.5,0.5);
  
    
  AliRsnCutMiniPair* cutV0=new AliRsnCutMiniPair("cutV0",AliRsnCutMiniPair::kContainsV0Daughter);
    
  AliRsnCutSet* cutsPairSameKS=new AliRsnCutSet("pairCutsSameKS",AliRsnTarget::kMother);
  cutsPairSameKS->AddCut(cutYKS);
  cutsPairSameKS->AddCut(cutV0);
  cutsPairSameKS->SetCutScheme(TString::Format("%s&(!%s)",cutYKS->GetName(),cutV0->GetName()).Data());
    
  AliRsnCutSet* cutsPairMixKS=new AliRsnCutSet("pairCutsMixKS", AliRsnTarget::kMother);
  cutsPairMixKS->AddCut(cutYKS);
  cutsPairMixKS->SetCutScheme(cutYKS->GetName());





/* invariant mass   */ Int_t imID   = task->CreateValue(AliRsnMiniValue::kInvMass,    kFALSE);
/* IM resolution    */ Int_t resID  = task->CreateValue(AliRsnMiniValue::kInvMassRes, kTRUE);
/* transv. momentum */  Int_t ptID   = task->CreateValue(AliRsnMiniValue::kPt,         kFALSE);
/* centrality       */  Int_t centID = task->CreateValue(AliRsnMiniValue::kMult,       kFALSE);
/* pseudorapidity   */  Int_t etaID  = task->CreateValue(AliRsnMiniValue::kEta,        kFALSE);
/* rapidity         */  Int_t yID    = task->CreateValue(AliRsnMiniValue::kY,          kFALSE);







  Bool_t  use     [10] = {1         ,1         ,1         ,1       ,isMC     ,isMC     ,isMC     ,isMC     ,isMC    ,isMC    };
  Bool_t  useIM   [10] = {1         ,1         ,1         ,1       ,1        ,1        ,1        ,1        ,0       ,0       };
  TString name    [10] = {"KStar","AKSTar","KStarmix","AKStarmix","KStarMCGen","AKStarMCGen","KStarTrues","AKStarTrues","ResPM" ,"ResMP" };
  TString comp    [10] = {"PAIR"    ,"PAIR"    ,"MIX"     ,"MIX" ,"MOTHER" ,"MOTHER" ,"TRUE"   ,"TRUE"   ,"TRUE"  ,"TRUE"  };
  TString output  [10] = {"SPARSE"  ,"SPARSE"  ,"SPARSE"  ,"SPARSE"  ,"SPARSE","SPARSE","SPARSE" ,"SPARSE" ,"SPARSE" ,"SPARSE"};
  Char_t  charge1 [10] = {'0'       ,'0'       ,'0'       ,'0'       ,'0'     ,'0'     ,'0'      ,'0'      ,'0'      ,'0'     };
  Char_t  charge2 [10] = {'+'       ,'-'       ,'+'       ,'-'       ,'+'     ,'-'     ,'+'      ,'-'      ,'+'      ,'-'     };  
  Int_t   cutIDK  [10] = {iCutK0s     ,iCutK0s     ,iCutK0s     ,iCutK0s     ,iCutK0s   ,iCutK0s   ,iCutK0s    ,iCutK0s    ,iCutK0s    ,iCutK0s };
  Int_t   cutIDPi [10] = {iCutPi    ,iCutPi    ,iCutPi    ,iCutPi    ,iCutPi  ,iCutPi  ,iCutPi   ,iCutPi   ,iCutPi   ,iCutPi  };
  Int_t   PDGCode [10] = {323       ,-323      ,323       ,-323       ,323     ,-323     ,323      ,-323     ,323      ,-323    };
  AliRsnCutSet* paircuts[10] = {cutsPairSameKS,  cutsPairSameKS,   cutsPairMixKS,   cutsPairMixKS,    cutsPairSameKS,   cutsPairSameKS, cutsPairSameKS,   cutsPairSameKS, cutsPairSameKS,   cutsPairSameKS  };






 for (Int_t i = 0; i < 10; i++) {
    if (!use[i]) continue;
    AliRsnMiniOutput *out = task->CreateOutput(Form("CustomId%d_%s", customQualityCutsID, name[i].Data()), output[i].Data(), comp[i].Data());
    out->SetDaughter(0, AliRsnDaughter::kKaon0);
    out->SetDaughter(1, AliRsnDaughter::kPion);
    out->SetCutID(0, cutIDK[i]);
    out->SetCutID(1, cutIDPi[i]);
    out->SetCharge(0, charge1[i]);
    out->SetCharge(1, charge2[i]);
    out->SetMotherPDG(PDGCode[i]);
    out->SetMotherMass(0.891);
    out->SetPairCuts(paircuts[i]);

    // axis X: invmass (or resolution)
    if (useIM[i]) 
      out->AddAxis(imID, 180, 0.6, 1.5);
     else
     out->AddAxis(resID, 200, -0.02, 0.02);
    // axis Y: transverse momentum
    out->AddAxis(ptID, 500, 0.0, 50.0);
    // axis Z: centrality-multiplicity
    //if (!isPP)
    out->AddAxis(centID, 100, 0.0, 100.0);
  }
   return kTRUE;
  //  return kFALSE;
}






Bool_t ConfigRsnXiQA_ROOT6(AliRsnMiniAnalysisTask *task, Bool_t isMC, Bool_t isPP, AliRsnCutSet *cutsPair, Int_t Strcut, Int_t customQualityCutsID, AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutKaCandidate, Float_t nsigmaPi, Float_t nsigmaK, Float_t nsigmaTOF, Bool_t enableMonitor)

{








  AliRsnCutSetDaughterParticle* cutSetPi=NULL;
  AliRsnCutTrackQuality* trkQualityCut= new AliRsnCutTrackQuality("myQualityCut");
  if(!trkQualityCut) return kFALSE;


  // if(SetCustomQualityCut(trkQualityCut,customQualityCutsID,Strcut)){                                                                       

 cutSetPi=new AliRsnCutSetDaughterParticle(Form("cutPi%i_%2.1fsigmaTPC_%2.1fsigmaTOF",cutKaCandidate,nsigmaPi,nsigmaTOF),trkQualityCut,cutKaCandidate,AliPID::kPion,nsigmaPi,nsigmaTOF);
  //}                                                                                                                                       
  
Int_t iCutPi = task->AddTrackCuts(cutSetPi);



  ////////////////////////////////////////////////////////////                                                                               
  // selections for K0s and for the daughters of K0s                                                                                         
  /////////////////////////////////////////////////////////////                                                                              
  // 
  Float_t crossedRows=70;
  Float_t rowsbycluster=0.8;
  Float_t DCAxy=0.06;

                                                                                                                                        
  // selections for pion daugthers of K0s                                                                                                    
  AliESDtrackCuts *esdTrackCuts = new AliESDtrackCuts("qualityDaughterK0s");
  esdTrackCuts->SetEtaRange(-0.8,0.8);
  esdTrackCuts->SetRequireTPCRefit();
  esdTrackCuts->SetAcceptKinkDaughters(0); //                                                                                                
  esdTrackCuts->SetMinNCrossedRowsTPC(crossedRows);
  esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(rowsbycluster);
  esdTrackCuts->SetMaxChi2PerClusterTPC(100);
  esdTrackCuts->SetMinDCAToVertexXY(DCAxy); //Use one of the two - pt dependent or fixed value cut.                                          








  // selections for Xi
  Float_t XiPIDcut=3.;
  Float_t V0dDCA=1.6;
  Float_t XidDCA=1.6;
  Float_t XiMinDCA=0.07;
  Float_t Xi_massTol=0.007;
  Float_t Xi_massTolVeto=0.007;
  Float_t V0CosPoinAn=0.97;
  Float_t XiCosPoinAn=0.97;
    
  AliRsnCutCascade* cutXi=new AliRsnCutCascade("cutXi",kXiMinus);
  cutXi->SetPIDCutV0Proton(XiPIDcut);
  cutXi->SetPIDCutV0Pion(XiPIDcut);
  cutXi->SetPIDCutBachelor(XiPIDcut);
  cutXi->SetESDtrackCuts(esdTrackCuts);
  cutXi->SetV0MaxDaughtersDCA(V0dDCA);
  cutXi->SetCascadeMaxDaughtersDCA(XidDCA);
  cutXi->SetV0MaxDCAVertex(1e5); // not using
  cutXi->SetV0MinDCAVertex(XiMinDCA);
  cutXi->SetCascadeMaxDCAVertex(1e5); // not using
  cutXi->SetCascadeMinDCAVertex(-1e5); // not using
  cutXi->SetV0LowRadius(0); // not using
  cutXi->SetV0HighRadius(1e5); // not using
  cutXi->SetCascadeLowRadius(0); // not using
  cutXi->SetCascadeHighRadius(1e5); // not using
  cutXi->SetMassTolerance(Xi_massTol);
  cutXi->SetMassToleranceVeto(Xi_massTolVeto);//Rejection range for Competing Xi Rejection
  cutXi->SetSwitch(kFALSE); // not using
  cutXi->SetV0MinCosPointingAngle(V0CosPoinAn);
  cutXi->SetCascadeMinCosPointingAngle(XiCosPoinAn);
  //cutXi->SetMaxPseudorapidity(0.8);
  cutXi->SetMaxRapidity(0.8);
  cutXi->SetMinTPCcluster(-1);
    
  AliRsnCutCascade* cutXibar=new AliRsnCutCascade("cutXibar",kXiPlusBar);

  cutXibar->SetPIDCutV0Proton(XiPIDcut);
  cutXibar->SetPIDCutV0Pion(XiPIDcut);
  cutXibar->SetPIDCutBachelor(XiPIDcut);
  cutXibar->SetESDtrackCuts(esdTrackCuts);
  cutXibar->SetV0MaxDaughtersDCA(V0dDCA);
  cutXibar->SetCascadeMaxDaughtersDCA(XidDCA);
  cutXibar->SetV0MaxDCAVertex(1e5); // not using
  cutXibar->SetV0MinDCAVertex(XiMinDCA);
  cutXibar->SetCascadeMaxDCAVertex(1e5); // not using
  cutXibar->SetCascadeMinDCAVertex(-1e5); // not using
  cutXibar->SetV0LowRadius(0); // not using
  cutXibar->SetV0HighRadius(1e5); // not using
  cutXibar->SetCascadeLowRadius(0); // not using
  cutXibar->SetCascadeHighRadius(1e5); // not using
  cutXibar->SetMassTolerance(Xi_massTol);
  cutXibar->SetMassToleranceVeto(Xi_massTolVeto);//Rejection range for Competing Xi Rejection
  cutXibar->SetSwitch(kFALSE); // not using
  cutXibar->SetV0MinCosPointingAngle(V0CosPoinAn);
  cutXibar->SetCascadeMinCosPointingAngle(XiCosPoinAn);
  //cutXibar->SetMaxPseudorapidity(0.8);
  cutXibar->SetMaxRapidity(0.8);
  cutXibar->SetMinTPCcluster(-1);


  AliRsnCutSet* cutSetXi=new AliRsnCutSet("setXi",AliRsnTarget::kDaughter);
  cutSetXi->AddCut(cutXi);
  cutSetXi->SetCutScheme(cutXi->GetName());
  Int_t iCutXi=task->AddTrackCuts(cutSetXi);
    
  AliRsnCutSet* cutSetXibar=new AliRsnCutSet("setXibar",AliRsnTarget::kDaughter);
  cutSetXibar->AddCut(cutXibar);
  cutSetXibar->SetCutScheme(cutXibar->GetName());
  Int_t iCutXibar=task->AddTrackCuts(cutSetXibar);



  if(enableMonitor){
    Printf("======== Monitoring cut AliRsnCutSetDaughterParticle enabled");
        
    AddMonitorOutput(isMC,cutSetPi->GetMonitorOutput());
    AddMonitorOutputV0(isMC, cutSetXi->GetMonitorOutput(), "Lambda", "nokine");
    AddMonitorOutputV0(isMC, cutSetXibar->GetMonitorOutput(), "AntiLambda", "nokine");
    AddMonitorOutputCascade(isMC, cutSetXi->GetMonitorOutput(), "Xi");
    AddMonitorOutputCascade(isMC, cutSetXibar->GetMonitorOutput(), "AntiXi");
  }






  // pair cuts
  AliRsnCutMiniPair* cutY2=new AliRsnCutMiniPair("cutRapidity", AliRsnCutMiniPair::kRapidityRange);
  cutY2->SetRangeD(-0.5,0.5);
  
    
  AliRsnCutMiniPair* cutV0=new AliRsnCutMiniPair("cutV0",AliRsnCutMiniPair::kContainsV0Daughter);
    
  AliRsnCutSet* cutsPairSame=new AliRsnCutSet("pairCutsSame",AliRsnTarget::kMother);
  cutsPairSame->AddCut(cutY2);
  cutsPairSame->AddCut(cutV0);
  cutsPairSame->SetCutScheme(TString::Format("%s&(!%s)",cutY2->GetName(),cutV0->GetName()).Data());
    
  AliRsnCutSet* cutsPairMix=new AliRsnCutSet("pairCutsMix", AliRsnTarget::kMother);
  cutsPairMix->AddCut(cutY2);
  cutsPairMix->SetCutScheme(cutY2->GetName());





/* invariant mass   */  Int_t imID   = task->CreateValue(AliRsnMiniValue::kInvMass,    kFALSE);
/* IM resolution    */  Int_t resID  = task->CreateValue(AliRsnMiniValue::kInvMassRes, kTRUE);
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
  Int_t   cutIDK  [12] = {iCutXibar     ,iCutXi     ,iCutXibar     ,iCutXi     ,iCutXibar   ,iCutXi   ,iCutXibar    ,iCutXi    ,iCutXibar    ,iCutXi    ,iCutXibar   ,iCutXi   };
  Int_t   cutIDPi [12] = {iCutPi    ,iCutPi    ,iCutPi    ,iCutPi    ,iCutPi  ,iCutPi  ,iCutPi   ,iCutPi   ,iCutPi   ,iCutPi   ,iCutPi  ,iCutPi  };
  Int_t   PDGCode [12] = {-3324       ,3324      ,-3324       ,+3324       ,3324     ,3324     ,-3324      ,3324     ,-3324      ,3324      ,-3324     ,3324   };

  AliRsnCutSet* paircuts[12] = {cutsPairSame,  cutsPairSame,   cutsPairMix,   cutsPairMix,    cutsPairSame,   cutsPairSame, cutsPairSame,   cutsPairSame, cutsPairSame,   cutsPairSame, cutsPairSame,   cutsPairSame   };






  for (Int_t i = 0; i < 12; i++) {
    if (!use[i]) continue;
    AliRsnMiniOutput *out = task->CreateOutput(Form("CustomId%d_%s", customQualityCutsID, name[i].Data()), output[i].Data(), comp[i].Data());
    out->SetDaughter(0, AliRsnDaughter::kXi);
    out->SetDaughter(1, AliRsnDaughter::kPion);
    out->SetCutID(0, cutIDK[i]);
    out->SetCutID(1, cutIDPi[i]);
    out->SetCharge(0, charge1[i]);
    out->SetCharge(1, charge2[i]);
    out->SetMotherPDG(PDGCode[i]);
    out->SetMotherMass(1.53);
    out->SetPairCuts(paircuts[i]);

    // axis X: invmass (or resolution)
    if (useIM[i]) 
      out->AddAxis(imID, 180, 1.4, 2.3);
     else
     out->AddAxis(resID, 200, -0.02, 0.02);
    // axis Y: transverse momentum
    out->AddAxis(ptID, 500, 0.0, 50.0);
    // axis Z: centrality-multiplicity
    //if (!isPP)
    out->AddAxis(centID, 100, 0.0, 100.0);
  }
   return kTRUE;
  //  return kFALSE;

}
















void AddMonitorOutput_P(TString name,TObjArray *mon,TString opt,AliRsnLoopDaughter *loop){
  AliRsnValueDaughter* a=new AliRsnValueDaughter(Form("%s_mom",name.Data()),AliRsnValueDaughter::kP);
  a->SetBins(0.,10.0,0.05);
  AliRsnListOutput* o=new AliRsnListOutput(Form("out_%s",a->GetName()),AliRsnListOutput::kHistoDefault);
  o->AddValue(a);
  if (mon) mon->Add(o);
  if (loop) loop->AddOutput(o);
}


void AddMonitorOutput_Pt(TString name,TObjArray *mon,TString opt,AliRsnLoopDaughter *loop){
  AliRsnValueDaughter* a=new AliRsnValueDaughter(Form("%s_pt",name.Data()),AliRsnValueDaughter::kPt);
  a->SetBins(0.,10.0,0.05);
  AliRsnListOutput* o=new AliRsnListOutput(Form("out_%s",a->GetName()),AliRsnListOutput::kHistoDefault);
  o->AddValue(a);
  if (mon) mon->Add(o);
  if (loop) loop->AddOutput(o);
}

void AddMonitorOutput_Eta(TString name,TObjArray *mon,TString opt,AliRsnLoopDaughter *loop){
  AliRsnValueDaughter* a=new AliRsnValueDaughter(Form("%s_eta",name.Data()),AliRsnValueDaughter::kEta);
  a->SetBins(-2.,2.,0.01);
  AliRsnListOutput* o=new AliRsnListOutput(Form("out_%s",a->GetName()),AliRsnListOutput::kHistoDefault);
  o->AddValue(a);
  if (mon) mon->Add(o);
  if (loop) loop->AddOutput(o);
}

void AddMonitorOutput_DCAxy(TString name,TObjArray *mon,TString opt,AliRsnLoopDaughter *loop){
  AliRsnValueDaughter* a=new AliRsnValueDaughter(Form("%s_dcaxy",name.Data()),AliRsnValueDaughter::kDCAXY);
  a->SetBins(-0.5,0.5,0.001);
  AliRsnListOutput* o=new AliRsnListOutput(Form("out_%s",a->GetName()),AliRsnListOutput::kHistoDefault);
  o->AddValue(a);
  if (mon) mon->Add(o);
  if (loop) loop->AddOutput(o);
}

void AddMonitorOutput_DCAz(TString name,TObjArray *mon,TString opt,AliRsnLoopDaughter *loop){
  AliRsnValueDaughter* a=new AliRsnValueDaughter(Form("%s_dcaz",name.Data()),AliRsnValueDaughter::kDCAZ);
  a->SetBins(-2.5,2.5,0.005);
  AliRsnListOutput* o=new AliRsnListOutput(Form("out_%s",a->GetName()),AliRsnListOutput::kHistoDefault);
  o->AddValue(a);
  if (mon) mon->Add(o);
  if (loop) loop->AddOutput(o);
}

void AddMonitorOutput_NclTPC(TString name,TObjArray *mon,TString opt,AliRsnLoopDaughter *loop){
  AliRsnValueDaughter* a=new AliRsnValueDaughter(Form("%s_NclTPC",name.Data()),AliRsnValueDaughter::kNTPCclusters);
  a->SetBins(-0.5,199.5,1);
  AliRsnListOutput* o=new AliRsnListOutput(Form("out_%s",a->GetName()),AliRsnListOutput::kHistoDefault);
  o->AddValue(a);
  if (mon) mon->Add(o);
  if (loop) loop->AddOutput(o);
}

void AddMonitorOutput_chi2TPC(TString name,TObjArray *mon,TString opt,AliRsnLoopDaughter *loop){
  AliRsnValueDaughter* a=new AliRsnValueDaughter(Form("%s_chi2TPC",name.Data()),AliRsnValueDaughter::kTPCchi2);
  a->SetBins(0.0,6,.1);
  AliRsnListOutput* o=new AliRsnListOutput(Form("out_%s",a->GetName()),AliRsnListOutput::kHistoDefault);
  o->AddValue(a);
  if (mon) mon->Add(o);
  if (loop) loop->AddOutput(o);
}



void AddMonitorOutput_V0Mass(TString name,TObjArray *mon,TString opt,AliRsnLoopDaughter *loop){
  AliRsnValueDaughter* a=new AliRsnValueDaughter(Form("%s_v0mass",name.Data()),AliRsnValueDaughter::kV0Mass);
  name.ToLower();
  if(name.Contains("k0")) a->SetBins(0.4,0.6,0.001);
  else if(name.Contains("lambda")) a->SetBins(1.08,1.16,0.001);
  else a->SetBins(0.,3.,0.01);
  AliRsnListOutput* o=new AliRsnListOutput(Form("out_%s",a->GetName()),AliRsnListOutput::kHistoDefault);
  o->AddValue(a);
  if (mon) mon->Add(o);
  if (loop) loop->AddOutput(o);
}

void AddMonitorOutput_V0DCA(TString name,TObjArray *mon,TString opt,AliRsnLoopDaughter *loop){
  AliRsnValueDaughter* a=new AliRsnValueDaughter(Form("%s_v0dca",name.Data()),AliRsnValueDaughter::kV0DCA);
  a->SetBins(0.0,0.4,0.001);
  AliRsnListOutput* o=new AliRsnListOutput(Form("out_%s",a->GetName()),AliRsnListOutput::kHistoDefault);
  o->AddValue(a);
  if (mon) mon->Add(o);
  if (loop) loop->AddOutput(o);
}

void AddMonitorOutput_V0Radius(TString name,TObjArray *mon,TString opt,AliRsnLoopDaughter *loop){
  AliRsnValueDaughter* a=new AliRsnValueDaughter(Form("%s_v0radius",name.Data()),AliRsnValueDaughter::kV0Radius);
  a->SetBins(0.0,200,0.2);
  AliRsnListOutput* o=new AliRsnListOutput(Form("out_%s",a->GetName()),AliRsnListOutput::kHistoDefault);
  o->AddValue(a);
  if (mon) mon->Add(o);
  if (loop) loop->AddOutput(o);
}

void AddMonitorOutput_V0Lifetime(TString name,TObjArray *mon,TString opt,AliRsnLoopDaughter *loop){
  AliRsnValueDaughter* a=new AliRsnValueDaughter(Form("%s_v0lifetime",name.Data()),AliRsnValueDaughter::kV0Lifetime);
  a->SetBins(0.0,200,0.1);
  AliRsnListOutput* o=new AliRsnListOutput(Form("out_%s",a->GetName()),AliRsnListOutput::kHistoDefault);
  o->AddValue(a);
  if (mon) mon->Add(o);
  if (loop) loop->AddOutput(o);
}

void AddMonitorOutput_V0DaughterDCA(TString name,TObjArray *mon,TString opt,AliRsnLoopDaughter *loop){
  AliRsnValueDaughter* a=new AliRsnValueDaughter(Form("%s_v0ddca",name.Data()),AliRsnValueDaughter::kDaughterDCA);
  a->SetBins(0.0,2,0.001);
  AliRsnListOutput* o=new AliRsnListOutput(Form("out_%s",a->GetName()),AliRsnListOutput::kHistoDefault);
  o->AddValue(a);
  if (mon) mon->Add(o);
  if (loop) loop->AddOutput(o);
}

void AddMonitorOutput_V0DCA2TPV(TString name,TObjArray *mon,TString opt,AliRsnLoopDaughter *loop){
  //DCA of secondary tracks to primary vertex
  AliRsnValueDaughter* a=new AliRsnValueDaughter(Form("%s_v0dca2tpv",name.Data()),AliRsnValueDaughter::kV0DCAXY);
  a->SetBins(-10.,10.,0.01);
  AliRsnListOutput* o=new AliRsnListOutput(Form("out_%s",a->GetName()),AliRsnListOutput::kHistoDefault);
  o->AddValue(a);
  if (mon) mon->Add(o);
  if (loop) loop->AddOutput(o);
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
