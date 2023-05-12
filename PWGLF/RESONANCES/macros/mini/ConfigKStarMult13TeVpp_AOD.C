//#ifdef __CLING__
#if !defined (__CINT__) || defined (__CLING__)
R__ADD_INCLUDE_PATH($ALICE_PHYSICS)
#include <PWGLF/RESONANCES/macros/mini/AddMonitorOutput.C>
#endif
/***************************************************************************
  priyanka.sett@cern.ch - last modified on 27/10/2016
  himani.bhatt@cern.ch  - last modified on 30/04/2018
 Modified by Sonali Padhan on 09/05/2023
// *** Configuration script for K* analysis for 13 TeV pp data  ***
//
// A configuration script for RSN package needs to define the followings:
//
// (1) decay tree of each resonance to be studied, which is needed to select
//     true pairs and to assign the right mass to all candidate daughters
// (2) cuts at all levels: single daughters, tracks, events
// (3) output objects: histograms or trees
****************************************************************************/
Bool_t SetCustomQualityCut(AliRsnCutTrackQuality* trkQualityCut,
Int_t customQualityCutsID = 0,
Int_t customFilterBit = 0);
Bool_t ConfigKStarMult13TeVpp_AOD(
            AliRsnMiniAnalysisTask *task,
            Bool_t                 isMC,
            Bool_t                 isPP,
            const char             *suffix,
            AliRsnCutSet           *cutsPair,
            Int_t                  aodFilterBit = 5,
            Int_t                  customQualityCutsID=AliRsnCutSetDaughterParticle::kDisableCustom,
            AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutKaCandidate=AliRsnCutSetDaughterParticle::kTPCpidTOFveto3s,
            Float_t                nsigmaPi  = 2.0,
            Float_t                nsigmaK   = 2.0,
            Float_t                nsigmaTOF = 3.0,
            Bool_t                 enableMonitor = kTRUE,
            UInt_t triggerMask = AliVEvent::kINT7
)
        
{

  Int_t       signedPdg = 313;
  TString     monitorOpt = "NoSIGN";  //Flag for AddMonitorOutput.C e.g."NoSIGN"
  Bool_t      useCrossedRows = kTRUE;
  const char *yaxisVar = "";  //yaxisVar = "PtDaughter_PDaughter_cent"
  Bool_t      useMixLS = 0;
  
//  // retrieve mass from PDG database
  Int_t         pdg  = 313;
TDatabasePDG *db   = TDatabasePDG::Instance();
    TParticlePDG *part = db->GetParticle(pdg);
 Double_t mass      = part->Mass();
    TString opt(suffix);
    if (strlen(suffix) > 0) opt = Form("_%s", suffix);
    Bool_t isPtDaughter  = opt.Contains("PTDAUGHTER") ;
    Bool_t isPDaughter   = opt.Contains("PDAUGHTER") ;
    Bool_t iscent        = opt.Contains("CENT") ;
    Bool_t iseta         = opt.Contains("ETA") ;
    Bool_t israpidity    = opt.Contains("RAPIDITY") ;
    // manage suffix
    if (strlen(suffix) > 0) suffix = Form("_%s", suffix);
    Int_t MultBins = aodFilterBit / 100;
    aodFilterBit = aodFilterBit % 100;

    // set daughter cuts
  AliRsnCutSetDaughterParticle* cutSetPi;
  AliRsnCutSetDaughterParticle* cutSetK;
  AliRsnCutTrackQuality* trkQualityCut= new AliRsnCutTrackQuality("myQualityCut");
  if(!trkQualityCut) return kFALSE;

  if(SetCustomQualityCut(trkQualityCut,customQualityCutsID,aodFilterBit)){
    cutSetPi=new AliRsnCutSetDaughterParticle(Form("cutPi%i_%2.1fsigma",cutKaCandidate,nsigmaPi),trkQualityCut,cutKaCandidate,AliPID::kPion,nsigmaPi, nsigmaTOF);
    cutSetK=new AliRsnCutSetDaughterParticle(Form("cutK%i_%2.1fsigma",cutKaCandidate, nsigmaK),trkQualityCut,cutKaCandidate,AliPID::kKaon,nsigmaK, nsigmaTOF);
  }
  else{
     printf("Doughter Track cuts has been selected =================\n");
       return kFALSE;
  }
  Int_t iCutPi = task->AddTrackCuts(cutSetPi);
  Int_t iCutK  = task->AddTrackCuts(cutSetK);
  if(enableMonitor){
    Printf("======== Monitoring cut AliRsnCutSetDaughterParticle enabled");
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/AddMonitorOutput.C");
    AddMonitorOutput(isMC, cutSetPi->GetMonitorOutput(),monitorOpt.Data());
    AddMonitorOutput(isMC, cutSetK->GetMonitorOutput(),monitorOpt.Data());
  }

  // -- Values ------------------------------------------------------------------------------------
  /* invariant mass   */ Int_t imID   = task->CreateValue(AliRsnMiniValue::kInvMass,    kFALSE);
  /* IM resolution    */ Int_t resID  = task->CreateValue(AliRsnMiniValue::kInvMassRes, kTRUE);
  /* transv. momentum */ Int_t ptID   = task->CreateValue(AliRsnMiniValue::kPt,         kFALSE);
  /* centrality       */ Int_t centID = task->CreateValue(AliRsnMiniValue::kMult,       kFALSE);
  /* pseudorapidity   */ Int_t etaID  = task->CreateValue(AliRsnMiniValue::kEta,        kFALSE);
  /* rapidity         */ Int_t yID    = task->CreateValue(AliRsnMiniValue::kY,          kFALSE);
    /* 1st daughter pt  */ Int_t fdpt   = task->CreateValue(AliRsnMiniValue::kFirstDaughterPt, kFALSE);
    /* 2nd daughter pt  */ Int_t sdpt   = task->CreateValue(AliRsnMiniValue::kSecondDaughterPt, kFALSE);
    /* 1st daughter p   */ Int_t fdp    = task->CreateValue(AliRsnMiniValue::kFirstDaughterP, kFALSE);
    /* 2nd daughter p   */ Int_t sdp    = task->CreateValue(AliRsnMiniValue::kSecondDaughterP, kFALSE);

    /* Dip Angle        */ Int_t OpAn   = task->CreateValue(AliRsnMiniValue::kDipAngle, kFALSE);
    /* kPtRatio         */ Int_t PtRat  = task->CreateValue(AliRsnMiniValue::kPtRatio, kFALSE);
    /* kptresolution    */ Int_t respT  = task->CreateValue(AliRsnMiniValue::kPairPtRes, kFALSE);
    /* kmassresolution  */ Int_t resmass  = task->CreateValue(AliRsnMiniValue::kInvMassRes, kFALSE);
    /* IM diff          */ Int_t IMdif   = task->CreateValue(AliRsnMiniValue::kInvMassDiff, kTRUE);
    
    Double_t multbins[200];
    int j,nmult=0;
    if(isMC){
      for(j=0;j<10;j++){multbins[nmult]=0.001*j; nmult++;}
      for(j=1;j<50;j++){multbins[nmult]=0.01*j; nmult++;}
      for(j=5;j<10;j++){multbins[nmult]=0.1*j; nmult++;}
      for(j=1;j<10;j++){multbins[nmult]=j; nmult++;}
      for(j=2;j<=20;j++){multbins[nmult]=5.*j; nmult++;}
    }
    else if(triggerMask==AliVEvent::kHighMultV0){
      for(j=0;j<10;j++){multbins[nmult]=0.001*j; nmult++;}
      for(j=1;j<50;j++){multbins[nmult]=0.01*j; nmult++;}
      for(j=5;j<=10;j++){multbins[nmult]=0.1*j; nmult++;}
    }
    else{
      for(j=0;j<10;j++){multbins[nmult]=0.1*j; nmult++;}
      for(j=1;j<10;j++){multbins[nmult]=j; nmult++;}
      for(j=2;j<=20;j++){multbins[nmult]=5.*j; nmult++;}
    }
    
  // -- Create all needed outputs -----------------------------------------------------------------
  // use an array for more compact writing, which are different on mixing and charges
  // [0] = unlike
  // [1] = mixing
  // [2] = like ++
  // [3] = like --

  Bool_t  use     [12] = {1         ,1         ,1         ,1         ,1       ,1       ,isMC     ,isMC     ,isMC     ,isMC     ,useMixLS    ,useMixLS};
  Bool_t  useIM   [12] = {1         ,1         ,1         ,1         ,1       ,1       ,1        ,1        ,1        ,1        ,1       ,1       };
  TString name    [12] = {"UnlikePM","UnlikeMP","MixingPM","MixingMP","LikePP","LikeMM","MCGenPM","MCGenMP","TruesPM","TruesMP","ResPM" ,"ResMP" };
  TString comp    [12] = {"PAIR"    ,"PAIR"    ,"MIX"     ,"MIX"     ,"PAIR"  ,"PAIR"  ,"MOTHER" ,"MOTHER" ,"TRUE"   ,"TRUE"   ,"TRUE"  ,"TRUE"  };
  TString output  [12] = {"SPARSE"  ,"SPARSE"  ,"SPARSE"  ,"SPARSE"  ,"SPARSE","SPARSE","SPARSE" ,"SPARSE" ,"SPARSE" ,"SPARSE" ,"SPARSE","SPARSE"};
  Char_t  charge1 [12] = {'+'       ,'-'       ,'+'       ,'-'       ,'+'     ,'-'     ,'+'      ,'-'      ,'+'      ,'-'      ,'+'     ,'-'     };
  Char_t  charge2 [12] = {'-'       ,'+'       ,'-'       ,'+'       ,'+'     ,'-'     ,'-'      ,'+'      ,'_'      ,'+'      ,'-'     ,'+'     };
  Int_t   cutIDPi [12] = {iCutPi    ,iCutPi    ,iCutPi    ,iCutPi    ,iCutPi  ,iCutPi  ,iCutPi   ,iCutPi   ,iCutPi   ,iCutPi   ,iCutPi  ,iCutPi  };
  Int_t   cutIDK  [12] = {iCutK     ,iCutK     ,iCutK     ,iCutK     ,iCutK   ,iCutK   ,iCutK    ,iCutK    ,iCutK    ,iCutK    ,iCutK   ,iCutK   };
  Int_t   PDGCode [12] = {313       ,313       ,313       ,313       ,313     ,313     ,313      ,-313     ,313      ,-313      ,313     ,-313    };

  
  for (Int_t i = 0; i < 12; i++) {
    if (!use[i]) continue;
    AliRsnMiniOutput *out = task->CreateOutput(Form("KSTAR_%s%s", name[i].Data(), suffix), output[i].Data(), comp[i].Data());
    out->SetDaughter(0, AliRsnDaughter::kKaon);
    out->SetDaughter(1, AliRsnDaughter::kPion);
    out->SetCutID(0, cutIDK[i]);
    out->SetCutID(1, cutIDPi[i]);
    out->SetCharge(0, charge1[i]);
    out->SetCharge(1, charge2[i]);
    out->SetMotherPDG(PDGCode[i]);
    out->SetMotherMass(0.89594);
    out->SetPairCuts(cutsPair);

    // axis X: invmass (or resolution)
    if (useIM[i])
      out->AddAxis(imID, 90, 0.6, 1.5);
    else
      out->AddAxis(resID, 200, -0.02, 0.02);
    
    // axis Y: transverse momentum
    out->AddAxis(ptID, 500, 0.0, 50.0);
    
    // axis Z: centrality-multiplicity
    if(!isPP || MultBins) out->AddAxis(centID,nmult,multbins);
       else out->AddAxis(centID,nmult,multbins);
     // if (!isPP)out->AddAxis(centID, 100, 0.0, 100.0);
    //else out->AddAxis(centID, 400, 0.0, 400.0);
      
    // axis W: pseudorapidity
    // out->AddAxis(etaID, 20, -1.0, 1.0);
    // axis J: rapidity
    //out->AddAxis(yID, 90, -4.5, 4.5);
    
  }
  return kTRUE;
}

Bool_t SetCustomQualityCut(AliRsnCutTrackQuality * trkQualityCut, Int_t customQualityCutsID = 0,Int_t trCut = 2011)
{
  //Sets configuration for track quality object different from std quality cuts.
  //Returns kTRUE if track quality cut object is successfully defined,
  //returns kFALSE if an invalid set of cuts (customQualityCutsID) is chosen or if the
  //object to be configured does not exist.

  if ((!trkQualityCut)){
    Printf("::::: SetCustomQualityCut:: use default quality cuts specified in task configuration.");
    return kFALSE;
  }

  if(customQualityCutsID>=1 && customQualityCutsID<100 && customQualityCutsID!=12){
    trkQualityCut->SetDefaults2011(kTRUE,kTRUE);
    Printf("::::: SetCustomQualityCut:: using standard 2011 track quality cuts.............");
      {
    
  if(customQualityCutsID==11){trkQualityCut->GetESDtrackCuts()->SetMaxDCAToVertexXYPtDep("0.015+0.05/pt^1.1");}//10Sig // D = 7*(0.0015+0.0050/pt^1.1)
    else if(customQualityCutsID==2){trkQualityCut->GetESDtrackCuts()->SetMaxDCAToVertexXYPtDep("0.006+0.02/pt^1.1");}//4Sig
    else if(customQualityCutsID==3){trkQualityCut->GetESDtrackCuts()->SetMaxDCAToVertexZ(0.2);}//Default=2cm
    else if(customQualityCutsID==4){trkQualityCut->GetESDtrackCuts()->SetMaxDCAToVertexZ(1.0);}
    else if(customQualityCutsID==5){trkQualityCut->GetESDtrackCuts()->SetMaxDCAToVertexZ(3.0);}
    else if(customQualityCutsID==6){trkQualityCut->GetESDtrackCuts()->SetMaxChi2PerClusterTPC(5.);}// D = 4
    else if(customQualityCutsID==7){trkQualityCut->GetESDtrackCuts()->SetMaxChi2PerClusterTPC(3);}
    else if(customQualityCutsID==8){trkQualityCut->GetESDtrackCuts()->SetMinNCrossedRowsTPC(100);}// D = 70
    else if(customQualityCutsID==9){trkQualityCut->GetESDtrackCuts()->SetMinNCrossedRowsTPC(80);}
    else if(customQualityCutsID==10){trkQualityCut->GetESDtrackCuts()->SetMinNCrossedRowsTPC(60);}
    else
    if(customQualityCutsID==13){trkQualityCut->GetESDtrackCuts()->SetMaxChi2PerClusterITS(49.);}// D = 36
        else
    if(customQualityCutsID==14){trkQualityCut->GetESDtrackCuts()->SetMaxChi2PerClusterITS(25.);}
    if(customQualityCutsID==15){trkQualityCut->GetESDtrackCuts()->SetMinRatioCrossedRowsOverFindableClustersTPC(0.7);}// D = 8
    else if(customQualityCutsID==16){trkQualityCut->GetESDtrackCuts()->SetMinRatioCrossedRowsOverFindableClustersTPC(0.9);}
     
    else if(customQualityCutsID==17){trkQualityCut->GetESDtrackCuts()->SetMaxChi2TPCConstrainedGlobal(49.);}// D = 36
    else if(customQualityCutsID==18){trkQualityCut->GetESDtrackCuts()->SetMaxChi2TPCConstrainedGlobal(25.);}
    else if(customQualityCutsID==19){trkQualityCut->GetESDtrackCuts()->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kOff);}
      }
    trkQualityCut->Print();
    return kTRUE;
  }else if(customQualityCutsID==12 || (customQualityCutsID>=100 && customQualityCutsID<200)){
    trkQualityCut->SetDefaultsTPCOnly(kTRUE);
    Printf(Form("::::: SetCustomQualityCut:: using TPC-only track quality cuts"));
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
      
      
      trkQualityCut->SetPtRange(0.15, 20.0);
      trkQualityCut->SetEtaRange(-0.8, 0.8);
      
      Printf("%s",Form("::::: SetCustomQualityCut:: using custom track quality cuts #%i",customQualityCutsID));
      trkQualityCut->Print();
      return kTRUE;
    }

