/***************************************************************************
              adash@cern.ch - last modified on 03/04/2013

// *** Configuration script for Phi-Meson analysis with 2013 PPb runs ***
// 
// A configuration script for RSN package needs to define the followings:
//
// (1) decay tree of each resonance to be studied, which is needed to select
//     true pairs and to assign the right mass to all candidate daughters
// (2) cuts at all levels: single daughters, tracks, events
// (3) output objects: histograms or trees
****************************************************************************/

Bool_t Configkstarpp13mult(  
			AliRsnMiniAnalysisTask *task, 
			Bool_t                 isMC, 
			Bool_t                 isPP,
			const char             *suffix,
			AliRsnCutSet           *cutsPair,
			Float_t                nsigmaPi = 3.0,
			Float_t                nsigmaK  = 3.0,
			Bool_t                 enableMonitor = kTRUE,
			TString                optSyt="DefaultITSTPC2011",
			UInt_t                 triggerMask=AliVEvent::kINT7
			  )
{
  
  //These are the Default values for 2011 ESD track cuts
  
  AliPID::EParticleType  type1   = AliPID::kPion;
  AliPID::EParticleType  type2   = AliPID::kKaon;

  // manage suffix
  if (strlen(suffix) > 0) suffix = Form("_%s", suffix);
  
  // retrieve mass from PDG database
  Int_t         pdg  = 313;
  TDatabasePDG *db   = TDatabasePDG::Instance();
  TParticlePDG *part = db->GetParticle(pdg);
  Double_t mass      = part->Mass();
  
  Float_t nsigmaKaTOF=3.0;
Float_t nsigmaPiTOF=3.0;

  // set daughter cuts
  AliRsnCutSetDaughterParticle* cutSetPi;
  AliRsnCutSetDaughterParticle* cutSetK;

  AliRsnCutTrackQuality* trkQualityCut= new AliRsnCutTrackQuality("myQualityCut");
  if(optSyt.Contains("DefaultITSTPC2011")){ 
    trkQualityCut->SetDefaults2011(kTRUE,kTRUE);
    Printf(Form("::::: SetCustomQualityCut:: using standard 2011 track quality cuts"));
  }
  else if(optSyt.Contains("DefaultTPCOnly")){
    trkQualityCut->SetDefaultsTPCOnly(kTRUE);
    Printf(Form("::::: SetCustomQualityCut:: using TPC-only track quality cuts"));
  }
  trkQualityCut->Print();
  
  //kTPCpidTOFveto4s
  //kTPCpidTOFveto3s
  //kFastTPCpidNsigma
    cutSetPi=new AliRsnCutSetDaughterParticle(Form("cutPi%i_%2.1fsigma",AliRsnCutSetDaughterParticle::kTPCTOFpidphipp2015,nsigmaPi),trkQualityCut,AliRsnCutSetDaughterParticle::kTPCTOFpidphipp2015,AliPID::kPion,nsigmaPi,nsigmaPiTOF);
  cutSetK=new AliRsnCutSetDaughterParticle(Form("cutK%i_%2.1fsigma",AliRsnCutSetDaughterParticle::kTPCTOFpidphipp2015, nsigmaK),trkQualityCut,AliRsnCutSetDaughterParticle::kTPCTOFpidphipp2015,AliPID::kKaon,nsigmaK,nsigmaKaTOF);
  //cutSetPi=new AliRsnCutSetDaughterParticle(Form("cutPi%i_%2.1fsigma",AliRsnCutSetDaughterParticle::kTPCpidTOFveto4s,nsigmaPi),trkQualityCut,AliRsnCutSetDaughterParticle::kTPCpidTOFveto4s,AliPID::kPion,nsigmaPi);
  //cutSetK=new AliRsnCutSetDaughterParticle(Form("cutK%i_%2.1fsigma",AliRsnCutSetDaughterParticle::kTPCpidTOFveto4s, nsigmaK),trkQualityCut,AliRsnCutSetDaughterParticle::kTPCpidTOFveto4s,AliPID::kKaon,nsigmaK);

  Int_t iCutPi = task->AddTrackCuts(cutSetPi);
  Int_t iCutK  = task->AddTrackCuts(cutSetK);
  
  if(enableMonitor){
    Printf("======== Monitoring cut AliRsnCutSetDaughterParticle enabled");
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/AddMonitorOutput.C");
    AddMonitorOutput(isMC, cutSetPi->GetMonitorOutput());
    AddMonitorOutput(isMC, cutSetK->GetMonitorOutput());
  }  
  // -- Values ------------------------------------------------------------------------------------
  /* invariant mass   */ Int_t imID   = task->CreateValue(AliRsnMiniValue::kInvMass,    kFALSE);
  /* IM resolution    */ Int_t resID  = task->CreateValue(AliRsnMiniValue::kInvMassRes, kTRUE);
  /* transv. momentum */ Int_t ptID   = task->CreateValue(AliRsnMiniValue::kPt,         kFALSE);
  /* centrality       */ Int_t centID = task->CreateValue(AliRsnMiniValue::kMult,       kFALSE);
  /* pseudorapidity   */ Int_t etaID  = task->CreateValue(AliRsnMiniValue::kEta,        kFALSE);
  /* rapidity         */ Int_t yID    = task->CreateValue(AliRsnMiniValue::kY,          kFALSE);
  // -- Create all needed outputs -----------------------------------------------------------------
  // use an array for more compact writing, which are different on mixing and charges
  // [0] = unlike
  // [1] = mixing
  // [2] = like ++
  // [3] = like --

  Bool_t  use     [12] = {1         ,1         ,1         ,1         ,1       ,1       ,isMC     ,isMC     ,isMC     ,isMC     ,isMC    ,isMC    };
  Bool_t  useIM   [12] = {1         ,1         ,1         ,1         ,1       ,1       ,1        ,1        ,1        ,1        ,0       ,0       };
  TString name    [12] = {"UnlikePM","UnlikeMP","MixingPM","MixingMP","LikePP","LikeMM","MCGenPM","MCGenMP","TruesPM","TruesMP","ResPM" ,"ResMP" };
  TString comp    [12] = {"PAIR"    ,"PAIR"    ,"MIX"     ,"MIX"     ,"PAIR"  ,"PAIR"  ,"MOTHER" ,"MOTHER" ,"TRUE"   ,"TRUE"   ,"TRUE"  ,"TRUE"  };
  TString output  [12] = {"SPARSE"  ,"SPARSE"  ,"SPARSE"  ,"SPARSE"  ,"SPARSE","SPARSE","SPARSE" ,"SPARSE" ,"SPARSE" ,"SPARSE" ,"SPARSE","SPARSE"};
  Char_t  charge1 [12] = {'+'       ,'-'       ,'+'       ,'-'       ,'+'     ,'-'     ,'+'      ,'-'      ,'+'      ,'-'      ,'+'     ,'-'     };
  Char_t  charge2 [12] = {'-'       ,'+'       ,'-'       ,'+'       ,'+'     ,'-'     ,'-'      ,'+'      ,'_'      ,'+'      ,'-'     ,'+'     };
  Int_t   cutIDPi [12] = {iCutPi    ,iCutPi    ,iCutPi    ,iCutPi    ,iCutPi  ,iCutPi  ,iCutPi   ,iCutPi   ,iCutPi   ,iCutPi   ,iCutPi  ,iCutPi  };
  Int_t   cutIDK  [12] = {iCutK     ,iCutK     ,iCutK     ,iCutK     ,iCutK   ,iCutK   ,iCutK    ,iCutK    ,iCutK    ,iCutK    ,iCutK   ,iCutK   };
  Int_t   PDGCode [12] = {313       ,313       ,313       ,313       ,313     ,313     ,313      ,-313     ,313      ,-313      ,313     ,-313    };
  
  Double_t multbins[200];
  int j,nmult=0;
  if(triggerMask==AliVEvent::kHighMultV0){
    for(j=0;j<10;j++){multbins[nmult]=0.001*j; nmult++;}
    for(j=1;j<10;j++){multbins[nmult]=0.01*j; nmult++;}
    for(j=1;j<=10;j++){multbins[nmult]=0.1*j; nmult++;}
  }else{
    for(j=0;j<10;j++){multbins[nmult]=0.1*j; nmult++;}
    for(j=1;j<=100;j++){multbins[nmult]=j; nmult++;}
  }
  
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
    out->SetMotherMass(mass);
    out->SetPairCuts(cutsPair);

    // axis X: invmass (or resolution)
    if (useIM[i]) 
      out->AddAxis(imID, 240, 0.6, 3.0);
    else
      out->AddAxis(resID, 200, -0.02, 0.02);
    
    // axis Y: transverse momentum
    out->AddAxis(ptID, 300, 0.0, 30.0);
    
    // axis Z: centrality-multiplicity
    if (!isPP)
      out->AddAxis(centID,nmult,multbins);
    else 
      out->AddAxis(centID,nmult,multbins);
    // axis W: pseudorapidity
    // out->AddAxis(etaID, 20, -1.0, 1.0);
    // axis J: rapidity
    //out->AddAxis(yID, 90, -4.5, 4.5);
    
  }
  return kTRUE;
}
