
/***************************************************************************
nasir.mehdi.malik@cern.ch - last modified on 
  priyanka.sett@cern.ch - last modified on 27/10/2016
  himani.bhatt@cern.ch  - last modified on 30/04/2018
// *** Configuration script for L*, anti-L*, syst. analysis for 13 TeV pp data  ***
// 
// A configuration script for RSN package needs to define the followings:
//
// (1) decay tree of each resonance to be studied, which is needed to select
//     true pairs and to assign the right mass to all candidate daughters
// (2) cuts at all levels: single daughters, tracks, events
// (3) output objects: histograms or trees
****************************************************************************/

#if !defined (__CINT__) || defined (__CLING__)
R__ADD_INCLUDE_PATH($ALICE_PHYSICS)
#include <PWGLF/RESONANCES/macros/mini/AddMonitorOutput.C>

Bool_t SetCustomQualityCut(AliRsnCutTrackQuality* trkQualityCut,
                           Int_t customQualityCutsID = 0,
                           Int_t customFilterBit = 0);

#endif
Bool_t ConfigureLstarpp_sp
(  
    AliRsnMiniAnalysisTask *task, 
    Bool_t                 isMC, 
    Bool_t                 isPP,
    const char             *suffix,
    AliRsnCutSet           *cutsPair,
    Int_t                  aodFilterBit = 5,
     Int_t                  _spBin=500,
    Int_t                  customQualityCutsID=AliRsnCutSetDaughterParticle::kDisableCustom,
    AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutPrCandidate = AliRsnCutSetDaughterParticle::kTPCTOFpidLstar13ppTeV,
    AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutKaCandidate = AliRsnCutSetDaughterParticle::kTPCTOFpidLstar13ppTeV,
    Float_t                nsigmaPr = 3.0,
    Float_t                nsigmaKa = 3.0,
    Float_t                nsigmaTOFPr = 3.0,
    Float_t                nsigmaTOFKa = 3.0,
    Bool_t                 enableMonitor = kTRUE,
    UInt_t triggerMask = AliVEvent::kINT7
 )
{

 
   // retrieve mass from PDG database
  // Int_t         pdg  = 3124;
  TDatabasePDG *db   = TDatabasePDG::Instance();
  TParticlePDG *part = db->GetParticle(3124);
  Double_t mass      = part->Mass();
  
  TString     monitorOpt = "NoSIGN";  //Flag for AddMonitorOutput.C e.g."NoSIGN"
  Bool_t      useCrossedRows = kTRUE;
  const char *yaxisVar = "";  //yaxisVar = "PtDaughter_PDaughter_cent"
  Bool_t      useMixLS = 0;
 
  // TString opt(yaxisVar);
  //opt.ToUpper();

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
  AliRsnCutSetDaughterParticle * cutSetP;
  AliRsnCutSetDaughterParticle * cutSetK;


  AliRsnCutTrackQuality * trkQualityCut =  new AliRsnCutTrackQuality("myQualityCut");
  if(SetCustomQualityCut(trkQualityCut,customQualityCutsID,aodFilterBit)){
    //Set custom quality cuts for systematic checks
    cutSetK=new AliRsnCutSetDaughterParticle(Form("cutKaon_%i_%2.1fsigma",cutKaCandidate, nsigmaKa),trkQualityCut,cutKaCandidate,AliPID::kKaon,nsigmaKa,nsigmaTOFKa);
    cutSetP=new AliRsnCutSetDaughterParticle(Form("cutProton_%i_%2.1fsigma",cutPrCandidate, nsigmaPr),trkQualityCut,cutPrCandidate,AliPID::kProton,nsigmaPr,nsigmaTOFPr);
    
  }
  
  else{
    
    Bool_t useCrossedRows = 1;
    cutSetK=new AliRsnCutSetDaughterParticle(Form("cutK%i_%2.1fsigma",cutKaCandidate,nsigmaKa),cutKaCandidate,AliPID::kKaon,nsigmaKa,aodFilterBit,useCrossedRows);
    cutSetP=new AliRsnCutSetDaughterParticle(Form("cutPr%i_%2.1fsigma",cutPrCandidate,nsigmaPr),cutKaCandidate,AliPID::kProton,nsigmaPr,aodFilterBit,useCrossedRows);
  }

  Int_t iCutP = task->AddTrackCuts(cutSetP);
  Int_t iCutK = task->AddTrackCuts(cutSetK);
  
  if (enableMonitor){
    Printf("======== Cut monitoring enabled");
#ifdef __CINT__
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/AddMonitorOutput.C");
#endif
    AddMonitorOutput(isMC, cutSetP->GetMonitorOutput(), monitorOpt.Data());
    AddMonitorOutput(isMC, cutSetK->GetMonitorOutput(), monitorOpt.Data());
  }  
  
  // -- Values ------------------------------------------------------------------------------------
  /* invariant mass   */ Int_t imID   = task->CreateValue(AliRsnMiniValue::kInvMass, kFALSE);
  /* IM resolution    */ Int_t resID  = task->CreateValue(AliRsnMiniValue::kInvMassRes, kTRUE);
  /* transv. momentum */ Int_t ptID   = task->CreateValue(AliRsnMiniValue::kPt, kFALSE);
  /* centrality       */ Int_t centID = task->CreateValue(AliRsnMiniValue::kMult, kFALSE);
  /* pseudorapidity   */ Int_t etaID  = task->CreateValue(AliRsnMiniValue::kEta, kFALSE);
  /* rapidity         */ Int_t yID    = task->CreateValue(AliRsnMiniValue::kY, kFALSE);
   /* spherocity      */ Int_t SpherocityID = task->CreateValue(AliRsnMiniValue::kSpherocity,kFALSE);
  /* 1st daughter pt  */ Int_t fdpt   = task->CreateValue(AliRsnMiniValue::kFirstDaughterPt, kFALSE);
  /* 2nd daughter pt  */ Int_t sdpt   = task->CreateValue(AliRsnMiniValue::kSecondDaughterPt, kFALSE);
  /* 1st daughter p   */ Int_t fdp    = task->CreateValue(AliRsnMiniValue::kFirstDaughterP, kFALSE);
  /* 2nd daughter p   */ Int_t sdp    = task->CreateValue(AliRsnMiniValue::kSecondDaughterP, kFALSE);

  /* Dip Angle        */ Int_t OpAn   = task->CreateValue(AliRsnMiniValue::kDipAngle, kFALSE);
  /* kPtRatio         */ Int_t PtRat  = task->CreateValue(AliRsnMiniValue::kPtRatio, kFALSE);


  /* kptresolution    */ Int_t respT  = task->CreateValue(AliRsnMiniValue::kPairPtRes, kFALSE);
  /* kmassresolution  */ Int_t resmass  = task->CreateValue(AliRsnMiniValue::kInvMassRes, kFALSE);
  /* IM diff          */ Int_t IMdif   = task->CreateValue(AliRsnMiniValue::kInvMassDiff, kTRUE);


  // -- Create all needed outputs -----------------------------------------------------------------
  // use an array for more compact writing, which are different on mixing and charges
  // [0] = unlike
  // [1] = mixing
  // [2] = like ++
  // [3] = like --

  Bool_t  use     [12] = {1         ,1         ,1         ,1         ,1       ,1       ,isMC     ,isMC     ,isMC     ,isMC     ,isMC    ,isMC};
  Bool_t  useIM   [12] = {1         ,1         ,1         ,1         ,1       ,1       ,1        ,1        ,1        ,1        ,0       ,0       };
  TString name    [12] = {"UnlikePM","UnlikeMP","MixingPM","MixingMP","LikePP","LikeMM","MCGenPM","MCGenMP","TruesPM","TruesMP","ResPM" ,"ResMP" };
  TString comp    [12] = {"PAIR"    ,"PAIR"    ,"MIX"     ,"MIX"     ,"PAIR"  ,"PAIR"  ,"MOTHER" ,"MOTHER" ,"TRUE"   ,"TRUE"   ,"TRUE"  ,"TRUE"  };
  TString output  [12] = {"SPARSE"  ,"SPARSE"  ,"SPARSE"  ,"SPARSE"  ,"SPARSE","SPARSE","SPARSE" ,"SPARSE" ,"SPARSE" ,"SPARSE" ,"SPARSE","SPARSE"};
  Char_t  charge1 [12] = {'+'       ,'-'       ,'+'       ,'-'       ,'+'     ,'-'     ,'+'      ,'-'      ,'+'      ,'-'      ,'+'     ,'-'     };
  Char_t  charge2 [12] = {'-'       ,'+'       ,'-'       ,'+'       ,'+'     ,'-'     ,'-'      ,'+'      ,'_'      ,'+'      ,'-'     ,'+'     };
  Int_t   cutID1  [12] = {iCutP     ,iCutP     ,iCutP    ,iCutP      ,iCutP   ,iCutP   ,iCutP    ,iCutP    ,iCutP    ,iCutP    ,iCutP   ,iCutP   };
  Int_t   cutID2  [12] = {iCutK     ,iCutK     ,iCutK    ,iCutK      ,iCutK   ,iCutK   ,iCutK    ,iCutK    ,iCutK    ,iCutK    ,iCutK   ,iCutK   };
  Int_t   PDGCode [12] = {3124      ,3124      ,3124     ,3124       ,3124    ,3124    ,3124     ,-3124    ,3124     ,-3124    ,3124    ,-3124   };

   Double_t multbins[200];
  int j,nmult=0;
  if(triggerMask==AliVEvent::kHighMultV0){
    for(j=0;j<10;j++){multbins[nmult]=0.001*j; nmult++;}
    for(j=1;j<10;j++){multbins[nmult]=0.01*j; nmult++;}
    for(j=1;j<=10;j++){multbins[nmult]=0.1*j; nmult++;}
  }else{
    for(j=0;j<10;j++){multbins[nmult]=0.1*j; nmult++;}
    for(j=1;j<=10;j++){multbins[nmult]=j; nmult++;}
  }

  for (Int_t i = 0; i < 12; i++) {
    if (!use[i]) continue;
    AliRsnMiniOutput *out = task->CreateOutput(Form("Lstar_%s%s", name[i].Data(), suffix), output[i].Data(), comp[i].Data());
    out->SetCutID(0, cutID1[i]);
    out->SetCutID(1, cutID2[i]);
    out->SetDaughter(0, AliRsnDaughter::kProton);
    out->SetDaughter(1, AliRsnDaughter::kKaon);
    out->SetCharge(0, charge1[i]);
    out->SetCharge(1, charge2[i]);
    out->SetMotherPDG(PDGCode[i]);
    out->SetMotherMass(mass);
    out->SetPairCuts(cutsPair);

    // axis X: invmass (or resolution)
    if (useIM[i]) 
      out->AddAxis(imID, 120, 1.4, 2.0);
    else
      out->AddAxis(resID, 200, -0.02, 0.02);
    
    //axis Y: transverse momentum of pair as default - else chosen value
    out->AddAxis(ptID,200,0.,10.);//default use mother pt

    // axis Z: centrality-multiplicity
    // axis Z: centrality-multiplicity

     
     if(isPP && !MultBins)  out->AddAxis(centID, 400, 0., 400.);
     else out->AddAxis(centID,nmult,multbins);
    
   	  
    // axis W: pseudorapidity
    // out->AddAxis(etaID, 20, -1.0, 1.0);
    // axis J: rapidity
    //out->AddAxis(yID, 12, -0.6, 0.6);

     out->AddAxis(SpherocityID, _spBin,0.,1.);
    
  }
  return kTRUE;
}

Bool_t SetCustomQualityCut(AliRsnCutTrackQuality * trkQualityCut, Int_t customQualityCutsID = 0, Int_t customFilterBit = 0)
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
    
   if(!customFilterBit)
    {//ESD
      if(customQualityCutsID==11){trkQualityCut->GetESDtrackCuts()->SetMaxDCAToVertexXYPtDep("0.0150+0.0500/pt^1.1");}
      else if(customQualityCutsID==2){trkQualityCut->GetESDtrackCuts()->SetMaxDCAToVertexXYPtDep("0.006+0.0200/pt^1.1");}

      else if(customQualityCutsID==3){trkQualityCut->GetESDtrackCuts()->SetMaxDCAToVertexZ(0.2);}
      else if(customQualityCutsID==4){trkQualityCut->GetESDtrackCuts()->SetMaxDCAToVertexZ(1.0);}

      else if(customQualityCutsID==5){trkQualityCut->GetESDtrackCuts()->SetMaxChi2PerClusterTPC(2.3);}
      else if(customQualityCutsID==6){trkQualityCut->GetESDtrackCuts()->SetMaxChi2PerClusterTPC(3.5);}
      else if(customQualityCutsID==7){trkQualityCut->GetESDtrackCuts()->SetMaxChi2PerClusterTPC(3.);}

      else if(customQualityCutsID==8){trkQualityCut->GetESDtrackCuts()->SetMinNCrossedRowsTPC(100);}
      else if(customQualityCutsID==9){trkQualityCut->GetESDtrackCuts()->SetMinNCrossedRowsTPC(80);}
      else if(customQualityCutsID==10){trkQualityCut->GetESDtrackCuts()->SetMinNCrossedRowsTPC(60);}
     
      else if(customQualityCutsID==13){trkQualityCut->GetESDtrackCuts()->SetMaxChi2PerClusterITS(49.);}
      else if(customQualityCutsID==14){trkQualityCut->GetESDtrackCuts()->SetMaxChi2PerClusterITS(4.);}
      else if(customQualityCutsID==15){trkQualityCut->GetESDtrackCuts()->SetMaxChi2PerClusterITS(25.);}
      
      else if(customQualityCutsID==16){trkQualityCut->GetESDtrackCuts()->SetMaxChi2TPCConstrainedGlobal(49.);}
      else if(customQualityCutsID==17){trkQualityCut->GetESDtrackCuts()->SetMaxChi2TPCConstrainedGlobal(25.);}
      //else if(customQualityCutsID==18){trkQualityCut->GetESDtrackCuts()->SetSPDminNClusters(0);}
      
      else if(customQualityCutsID==19){trkQualityCut->GetESDtrackCuts()->SetMinRatioCrossedRowsOverFindableClustersTPC(0.9);}
      else if(customQualityCutsID==20){trkQualityCut->GetESDtrackCuts()->SetMinRatioCrossedRowsOverFindableClustersTPC(0.7);}
    }else{//AOD
       Printf("::::: SetCustomQualityCut:: using standard 2011 track quality cuts + AOD custom track cuts......");
            trkQualityCut->SetCheckOnlyFilterBit(kFALSE);
	     if(customQualityCutsID==11){trkQualityCut->SetDCARPtFormula("0.0150+0.0500/pt^1.1");}
	     else if(customQualityCutsID==2){trkQualityCut->SetDCARPtFormula("0.006+0.0200/pt^1.1");}
             else if(customQualityCutsID==3){trkQualityCut->SetDCAZmax(0.2);}
	     else if(customQualityCutsID==4){trkQualityCut->SetDCAZmax(1.0);}
             else if(customQualityCutsID==5){trkQualityCut->SetTPCmaxChi2(2.3);}
	     else if(customQualityCutsID==6){trkQualityCut->SetTPCmaxChi2(3.5);}
             else if(customQualityCutsID==7){trkQualityCut->SetTPCmaxChi2(3.);}   //
             else if(customQualityCutsID==8){trkQualityCut->SetMinNCrossedRowsTPC(100,kTRUE);}
	     else if(customQualityCutsID==9){trkQualityCut->SetMinNCrossedRowsTPC(80,kTRUE);}
	     else if(customQualityCutsID==10){trkQualityCut->SetMinNCrossedRowsTPC(60,kTRUE);}
	      else if(customQualityCutsID==13){trkQualityCut->SetITSmaxChi2(49.);}
	      else if(customQualityCutsID==14){trkQualityCut->SetITSmaxChi2(4.);}
	      else if(customQualityCutsID==15){trkQualityCut->SetITSmaxChi2(25.);}
	      else if(customQualityCutsID==16){trkQualityCut->SetMaxChi2TPCConstrainedGlobal(49.);}
	      else if(customQualityCutsID==17){trkQualityCut->SetMaxChi2TPCConstrainedGlobal(25.);}
              else if(customQualityCutsID==19){trkQualityCut->SetMinNCrossedRowsOverFindableClsTPC(0.9,kTRUE);}
	      else if(customQualityCutsID==20){trkQualityCut->SetMinNCrossedRowsOverFindableClsTPC(0.7,kTRUE);}
	      else if(customQualityCutsID==21){trkQualityCut->SetTrackMaxChi2(2.3);}
	      else if(customQualityCutsID==22){trkQualityCut->SetTrackMaxChi2(3.5);}
         
       
        }
    

    
    trkQualityCut->Print();
    return kTRUE;
  }else if(customQualityCutsID==12 || (customQualityCutsID>=100 && customQualityCutsID<200)){
    trkQualityCut->SetDefaultsTPCOnly(kTRUE);
    Printf("::::: SetCustomQualityCut:: using TPC-only track quality cuts");
    
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

