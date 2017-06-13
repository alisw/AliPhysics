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
		   AliRsnCutSet           *cutsPair,
		   Int_t                  aodFilterBit = 5,
		   Bool_t                 enableMonitor = kTRUE,
		   TString                monitorOpt = "NoSIGN",
		   Bool_t                 useGeoCutsPbPb2015 = kFALSE)
{
  //Configure phi
  TString    partname = "kPhi";
  Int_t      pdgCode  = 333;
  Float_t    mass     = 1.019461;
  Float_t    masslow  = 0.980;
  Float_t    massup   = 1.100;
  Int_t      nbins    = 120;
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

  AliRsnCutTrackQuality* trkQualityCut = new AliRsnCutTrackQuality("myQualityCut");
  trkQualityCut->SetDefaults2011(kTRUE,kTRUE);
  if (useGeoCutsPbPb2015) trkQualityCut->GetESDtrackCuts()->SetCutGeoNcrNcl(3., 130., 1.5, 0.85, 0.7);

  AliRsnCutSetDaughterParticle*  cutSetKa = new AliRsnCutSetDaughterParticle(Form("cutK%i_%2.1fsigma",AliRsnCutSetDaughterParticle::kFastTPCpidNsigma, nsigma),trkQualityCut,AliRsnCutSetDaughterParticle::kFastTPCpidNsigma,AliPID::kKaon,nsigma);
  Int_t iCutKa = task->AddTrackCuts(cutSetKa);

  if (enableMonitor){
    Printf("======== Cut monitoring enabled");
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/AddMonitorOutput.C");
    AddMonitorOutput(isMC, cutSetKa->GetMonitorOutput()), monitorOpt.Data();
  }

  // -- Create all needed outputs -----------------------------------------------------------------
  // use an array for more compact writing, which are different on mixing and charges
  // [0] = unlike
  // [1] = mixing

  Bool_t  use    [8] = { 1      , 1      , 1      ,isMC    , isMC   , isMC, isMC , isMC};
  TString xvar   [8] = { "M",    "M",    "M",     "M",     "M",     "Res",  "Eta", "Eta"};
  TString yvar   [8] = { "Pt",   "Pt",   "Pt",    "Pt",    "Pt",    "Pt",    "Pt", "Pt"};
  Bool_t  useIM  [8] = { 1      , 1      , 1      ,1       , 1      , 0   ,  1    , 1   };
  TString name   [8] = {"Unlike","LikePP","LikeMM","Gen"   ,"Trues" , "Res" , "Gen", "Trues"};
  TString comp   [8] = {"PAIR"  ,"PAIR"  ,"PAIR"  ,"MOTHER","TRUE"  , "TRUE", "MOTHER", "TRUE"};
  Char_t  charge1[8] = {'+'     ,'+'     ,'-'     , '+'    , '+'    , '+'   , '+'     , '+'};
  Char_t  charge2[8] = {'-'     ,'+'     ,'-'     , '-'    , '-'    , '-'   , '-'     , '-'};
  TString output = "HIST";

  for(Int_t i=0;i<8;i++){
    if(!use[i]) continue;
    TString histname = Form("%s_%s%s", name[i].Data(), xvar[i].Data(), yvar[i].Data());
    AliRsnMiniOutput* outpt = task->CreateOutput(histname.Prepend("phi_"), output.Data(), comp[i].Data());
    outpt->SetCutID(0,iCutKa);
    outpt->SetCutID(1,iCutKa);
    outpt->SetDaughter(0, d1);
    outpt->SetDaughter(1, d2);
    outpt->SetCharge(0,charge1[i]);
    outpt->SetCharge(1,charge2[i]);
    outpt->SetMotherPDG(pdgCode);
    outpt->SetMotherMass(mass);
    outpt->SetPairCuts(cutsPair);
    if (i<6) {
      if(useIM[i]) outpt->AddAxis(imID,nbins, masslow, massup);
      else outpt->AddAxis(resID, 200, -0.02, 0.02);
    } else {
      outpt->AddAxis(etaID, 20,-1.,1.);
    }
    outpt->AddAxis(ptID, 100, 0., 20.);
  }

  // if (isMC) {
  //   AliRsnMiniOutput * outgenrap = task->CreateOutput("phirap_gen", output.Data(),"MOTHER");
  //   outgenrap->SetDaughter(0, d1);
  //   outgenrap->SetDaughter(1, d2);
  //   outgenrap->SetMotherPDG(pdgCode);
  //   outgenrap->SetMotherMass(mass);
  //   outgenrap->SetPairCuts(cutsPair);
  //   outgenpt->AddAxis(ptID, 100, 0., 10.);
  //   outgenrap->AddAxis(etaID, 16,-0.8,0.8);

  //   AliRsnMiniOutput * outrecrap = task->CreateOutput("phirap_rec", output.Data(),"Trues");
  //   outrecrap->SetDaughter(0, d1);
  //   outrecrap->SetDaughter(1, d2);
  //   outrecrap->SetMotherPDG(pdgCode);
  //   outrecrap->SetMotherMass(mass);
  //   outrecrap->SetPairCuts(cutsPair);
  //   outrecpt->AddAxis(ptID, 100, 0., 10.);
  //   outrecrap->AddAxis(etaID, 16,-0.8,0.8);
  // }

  return kTRUE;
}
