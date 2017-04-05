/***************************************************************************
// Created on 22.03.2017 by
// Francesca Bellini f(bellini@cern.ch)
// Sourav Kundu (sourav.kundu@cern.ch)
//
// Configures Phi analysis with rsn mini package for QA
// Allows basic configuration of pile-up check and event cuts
//
****************************************************************************/

Bool_t ConfigRsnQA(AliRsnMiniAnalysisTask *task, 
		   Bool_t                 isMC, 
		   Bool_t                 isPP,
		   AliRsnCutSet           *cutsPair,
		   Int_t                  aodFilterBit = 5,
		   Bool_t                 enableMonitor = kTRUE,
		   TString                monitorOpt = "NoSIGN")
{
  //Configure phi
  TString    partname = "kPhi";
  Int_t      pdgCode  = 333;
  Float_t    mass     = 1.019461;
  Float_t    masslow  = 0.980;
  Float_t    massup   = 1.200;
  Int_t      nbins    = 220;
  RSNPID     d1 = AliRsnDaughter::kKaon;
  RSNPID     d2 = AliRsnDaughter::kKaon;

  // -- Values ------------------------------------------------------------------------------------
  AliRsnMiniValue::EType yaxisVar = AliRsnMiniValue::kPt;
  /* invariant mass   */ Int_t imID   = task->CreateValue(AliRsnMiniValue::kInvMass, kFALSE);
  /* IM resolution    */ Int_t resID  = task->CreateValue(AliRsnMiniValue::kInvMassRes, kTRUE);
  /* transv. momentum */ Int_t ptID   = task->CreateValue(AliRsnMiniValue::kPt, kFALSE);
  /* centrality       */ Int_t centID = task->CreateValue(AliRsnMiniValue::kMult, kFALSE);
  /* pseudorapidity   */ Int_t etaID  = task->CreateValue(AliRsnMiniValue::kEta, kFALSE);
  /* rapidity         */ Int_t yID    = task->CreateValue(AliRsnMiniValue::kY, kFALSE);

  // set daughter cuts
  //use default quality cuts std 2011 with crossed rows TPC  
  Bool_t useCrossedRows = 1; Float_t nsigma = 3.0;
  AliRsnCutSetDaughterParticle * cutSetKa = new AliRsnCutSetDaughterParticle("cutKa", AliRsnCutSetDaughterParticle::kFastTPCpidNsigma, AliPID::kKaon, nsigma, aodFilterBit, useCrossedRows);
  cutSetKa->SetUse2011StdQualityCuts(kTRUE);
  
  if (enableMonitor){
    Printf("======== Cut monitoring enabled");
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/AddMonitorOutput.C");
    AddMonitorOutput(isMC, cutSetKa->GetMonitorOutput()), monitorOpt.Data();
  }
  
  Int_t iCutKa = task->AddTrackCuts(cutSetKa);

 // -- Create all needed outputs -----------------------------------------------------------------
  // use an array for more compact writing, which are different on mixing and charges
  // [0] = unlike
  // [1] = mixing

  Bool_t  use    [6] = { 1      , 1      , 1      , 1      , isMC   , isMC};
  Bool_t  useIM  [6] = { 1      , 1      , 1      , 1      , 1      , 0};
  TString name   [6] = {"Unlike","Mixing","LikePP","LikeMM","Trues" ,"Res"};
  TString comp   [6] = {"PAIR"  , "MIX"  ,"PAIR"  ,"PAIR"  ,"TRUE"  ,"TRUE"};
  Char_t  charge1[6] = {'+'     , '+'    ,'+'     ,'-'     , '+'    , '+'};
  Char_t  charge2[6] = {'-'     , '-'    ,'+'     ,'-'     , '-'    , '-'};
  TString output = "SPARSE";

  for(Int_t i=0;i<6;i++){
    if(!use[i]) continue;

    AliRsnMiniOutput* out = task->CreateOutput(Form("phi_%s",name[i].Data()),output.Data(),comp[i].Data());
    out->SetCutID(0,iCutKa);
    out->SetCutID(1,iCutKa);
    out->SetDaughter(0, d1);
    out->SetDaughter(1, d2);
    out->SetCharge(0,charge1[i]);
    out->SetCharge(1,charge2[i]);
    out->SetMotherPDG(pdgCode);
    out->SetMotherMass(mass);
    out->SetPairCuts(cutsPair);
    //axis X: invmass (or resolution)
    if(useIM[i]) out->AddAxis(imID,nbins, masslow, massup);
    else out->AddAxis(resID, 200, -0.02, 0.02);
    //axis Y: transverse momentum of pair as default 
    out->AddAxis(ptID, 200, 0., 20.);
    // axis Z: centrality-multiplicity
    out->AddAxis(centID, 100, 0., 100.);
  }

  //GENERATED PAIRS
  AliRsnMiniOutput * outm = task->CreateOutput("phi_mother", output.Data(),"MOTHER");
  outm->SetDaughter(0, d1);
  outm->SetDaughter(1, d2);
  outm->SetMotherPDG(pdgCode);
  outm->SetMotherMass(mass);
  outm->SetPairCuts(cutsPair);
  outm->AddAxis(imID, nbins, masslow, massup);
  outm->AddAxis(ptID, 200, 0.0, 20.0);
  outm->AddAxis(centID, 100, 0.0, 100.0);

  return kTRUE;

}

  // if (enableMonitor){
  //   Printf("======== Cut monitoring enabled");
  //   gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/AddMonitorOutput.C");
  //   AddMonitorOutput(isMC, cutSetKa->GetMonitorOutput()), monitorOpt.Data();
  //   AddMonitorOutput(isMC, cutSetPi->GetMonitorOutput(), monitorOpt.Data());
  //   AddMonitorOutput(isMC, cutSetPro->GetMonitorOutput()), monitorOpt.Data();
  //}  

   /* Aditional cuts for pion and proton candidates
     AliRsnCutSetDaughterParticle * cutSetPi = new AliRsnCutSetDaughterParticle("cutPi", AliRsnCutSetDaughterParticle::kQualityStd2011, AliPID::kPion, nsigma, aodFilterBit, useCrossedRows);
     Int_t icutPi = task->AddTrackCuts(cutSetPi);
    
     AliRsnCutSetDaughterParticle * cutSetPro = new AliRsnCutSetDaughterParticle("cutPro", AliRsnCutSetDaughterParticle::kQualityStd2011, AliPID::kProton, nsigma, aodFilterBit, useCrossedRows);
     Int_t icutPro = task->AddTrackCuts(cutSetPro);

     if (d1==AliRsnDaughter::kProton) icut2 = icutPro;
     else if (d1==AliRsnDaughter::kKaon) icut2 = icutKa;
     else icut1 = icutPi;

     if (d2==AliRsnDaughter::kProton) icut2 = icutPro;
     else if (d2==AliRsnDaughter::kKaon) icut2 = icutKa;
     else icut2 = icutPi;
   */
