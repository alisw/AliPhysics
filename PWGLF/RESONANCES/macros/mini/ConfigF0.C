/***************************************************************************
fbellini@cern.ch - created on 25/08/2017
Configuration script for f0(980) analysis
****************************************************************************/
#ifdef __CLING__
R__ADD_INCLUDE_PATH($ALICE_PHYSICS)
#include <PWGLF/RESONANCES/macros/mini/AddMonitorOutput.C>
#endif

Bool_t ConfigF0(AliRsnMiniAnalysisTask *task = 0x0, 
		Bool_t                 isMC = kFALSE, 
		AliPIDResponse::EBeamType collSys = AliPIDResponse::kPP, //=0, kPPB=1, kPBPB=2
		AliRsnCutSet           *cutsPair = 0x0,             //cuts on the pair
		Bool_t                 enaMultSel = kTRUE,    //enable multiplicity axis
      		Float_t                masslow = 0.3,         //inv mass axis low edge 
		Float_t                massup = 1.3,          //inv mass axis upper edge 
		Int_t                  nbins = 1000,           //inv mass axis n bins
		Int_t                  aodFilterBit = 5,      //filter bit for AOD analysis
		AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutPid = AliRsnCutSetDaughterParticle::kTPCpidTOFveto3s, // pid cut set
		Float_t                nsigma = 3.0,          //nsigma of TPC PID cut
		Bool_t                 enableMonitor = kTRUE) //enable single track QA plots
{
  //-----------------------
  //General 
  //-----------------------
  TString partname="f0";
  Int_t   pdgCode=9010221;
  Float_t mass = 0.990;// 2014 PDG: M = 0.990 ± 20 GeV
  RSNPID  d1 = AliRsnDaughter::kPion;
  RSNPID  d2 = AliRsnDaughter::kPion;

  //Additional options for monitoring plots
  TString monitorOpt = "NoSIGN";
  
  //-----------------------
  // CUTS
  //-----------------------
  //use default quality cuts std 2010 with crossed rows TPC
  Bool_t useCrossedRows = 1; 
  AliRsnCutSetDaughterParticle * cutSetPi = new AliRsnCutSetDaughterParticle("cutPi", cutPid, AliPID::kPion, nsigma, aodFilterBit, useCrossedRows);
  cutSetPi->SetUse2011StdQualityCuts(kTRUE);
  Int_t icutPi = task->AddTrackCuts(cutSetPi);

  //set daughter cuts
  Int_t icut1 = icutPi;
  Int_t icut2 = icutPi;

  //monitor single-track selection based on track quality cuts only
  AliRsnCutSetDaughterParticle * cutSetQuality = new AliRsnCutSetDaughterParticle("cutQuality", AliRsnCutSetDaughterParticle::kQualityStd2011, AliPID::kPion, 10.0, aodFilterBit, useCrossedRows);
  Int_t icutQuality = task->AddTrackCuts(cutSetQuality);
 
  //QA plots 
  if (enableMonitor){
    Printf("======== Cut monitoring enabled");
#ifdef __CINT__
   gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/AddMonitorOutput.C");
#endif
    AddMonitorOutput(isMC, cutSetQuality->GetMonitorOutput(), monitorOpt.Data(), 0);
    AddMonitorOutput(isMC, cutSetPi->GetMonitorOutput(), monitorOpt.Data(), 0);
  }  

  //-----------------------
  // OUTPUT
  //-----------------------
  /* invariant mass   */ Int_t imID   = task->CreateValue(AliRsnMiniValue::kInvMass, kFALSE);
  /* IM resolution    */ Int_t resID  = task->CreateValue(AliRsnMiniValue::kInvMassRes, kTRUE);
  /* transv. momentum */ Int_t ptID   = task->CreateValue(AliRsnMiniValue::kPt, kFALSE);
  /* centrality       */ Int_t multID = task->CreateValue(AliRsnMiniValue::kMult, kFALSE);
  /* pseudorapidity   */ Int_t etaID  = task->CreateValue(AliRsnMiniValue::kEta, kFALSE);
  /* rapidity         */ Int_t yID    = task->CreateValue(AliRsnMiniValue::kY, kFALSE);
  /* 1st daughter pt  */ Int_t fdpt   = task->CreateValue(AliRsnMiniValue::kFirstDaughterPt, kFALSE);
  /* 2nd daughter pt  */ Int_t sdpt   = task->CreateValue(AliRsnMiniValue::kSecondDaughterPt, kFALSE);
  /* 1st daughter p   */ Int_t fdp    = task->CreateValue(AliRsnMiniValue::kFirstDaughterP, kFALSE);
  /* 2nd daughter p   */ Int_t sdp    = task->CreateValue(AliRsnMiniValue::kSecondDaughterP, kFALSE);
  TString output = "HIST"; // or "SPARSE"
  TString name[6] = {"UnlikePM", "MixingPM", "LikePP", "LikeMM", "MixingPP", "MixingMM"};
  TString comp[6] = {"PAIR"    , "MIX",      "PAIR"  , "PAIR"  , "MIX"     , "MIX"     };
  Char_t charge1[6] = {'+', '+', '+', '-', '+', '-'};
  Char_t charge2[6] = {'-', '-', '+', '-', '+', '-'};

  //DATA 
  for (Int_t i = 0; i < 6; i++) {
    AliRsnMiniOutput *out = task->CreateOutput(Form("f0_%s", name[i].Data()), output.Data(), comp[i].Data());
    out->SetCutID(0, icut1);
    out->SetCutID(1, icut2);
    out->SetDaughter(0, d1);
    out->SetDaughter(1, d2);
    out->SetCharge(0, charge1[i]);
    out->SetCharge(1, charge2[i]);
    out->SetMotherPDG(pdgCode);
    out->SetMotherMass(mass);
    out->SetPairCuts(cutsPair);
    // axis X: invmass 
    out->AddAxis(imID, nbins, masslow, massup);
    //axis Y: mother pt
    out->AddAxis(ptID, 200, 0.0, 20.0); //default use mother pt
    //axis Z: multiplicity
    if (enaMultSel) out->AddAxis(multID, 100, 0.0, 100.0);
  }
  

  //Template for BG
  TString bgTemplate[7]  = {"rho", "omega", "Kstar", "antiKstar", "K0s", "phi","f2"};
  Char_t bgTemplateC1[7] = {'+', '+', '+', '+', '+', '+', '+'};
  Char_t bgTemplateC2[7] = {'-', '-', '-', '-', '-', '-', '-'};
  Int_t bgTemplatePDG[7] = {113, 223, 313, -313, 310, 333, 225};
  Double_t bgTemplateM[7]   = {775.26, 8.49, 891.66, 891.66, 497.611, 1019.461, 1275.5};
  RSNPID bgID1[7] = {AliRsnDaughter::kPion, AliRsnDaughter::kPion,AliRsnDaughter::kKaon, AliRsnDaughter::kPion, AliRsnDaughter::kPion, AliRsnDaughter::kKaon, AliRsnDaughter::kPion};
  RSNPID bgID2[7] = {AliRsnDaughter::kPion, AliRsnDaughter::kPion,AliRsnDaughter::kPion, AliRsnDaughter::kKaon, AliRsnDaughter::kPion, AliRsnDaughter::kKaon, AliRsnDaughter::kPion};

  if (isMC) {
    //TRUE RECO PAIRS - TEMPLATE FOR BG
    for (Int_t ibg = 0; ibg <7; ibg++) {
      AliRsnMiniOutput * outtempl = task->CreateOutput(Form("bg_%s", bgTemplate[ibg].Data()), output.Data(),"TRUE");
      outtempl->SetCutID(0, icut1);
      outtempl->SetCutID(1, icut2);
      outtempl->SetCharge(0, bgTemplateC1[0]);
      outtempl->SetCharge(1, bgTemplateC2[0]);
      outtempl->SetDaughter(0, bgID1[ibg]);
      outtempl->SetDaughter(1, bgID2[ibg]);
      outtempl->SetMotherPDG(bgTemplatePDG[ibg]);
      outtempl->SetMotherMass(bgTemplateM[ibg]);
      outtempl->SetPairCuts(cutsPair);
      // axis X: invmass 
      outtempl->AddAxis(imID, nbins, masslow, massup);
      //axis Y: mother pt
      outtempl->AddAxis(ptID, 200, 0.0, 20.0); //default use mother pt
      // axis Z: multrality-multiplicity
      if (enaMultSel) outtempl->AddAxis(multID, 100, 0.0, 100.0);
    }
    //TRUE RECO PAIRS - MASS
    AliRsnMiniOutput * outtrue = task->CreateOutput(Form("truef0_%s", partname.Data()), output.Data(),"TRUE");
    outtrue->SetCutID(0, icut1);
    outtrue->SetCutID(1, icut2);
    outtrue->SetCharge(0, charge1[0]);
    outtrue->SetCharge(1, charge2[0]);
    outtrue->SetDaughter(0, d1);
    outtrue->SetDaughter(1, d2);
    outtrue->SetMotherPDG(pdgCode);
    outtrue->SetMotherMass(mass);
    outtrue->SetPairCuts(cutsPair);
    // axis X: invmass 
    outtrue->AddAxis(imID, nbins, masslow, massup);
    //axis Y: mother pt
    outtrue->AddAxis(ptID, 200, 0.0, 20.0); //default use mother pt
    // axis Z: multiplicity
    if (enaMultSel) outtrue->AddAxis(multID, 100, 0.0, 100.0);

    
    //TRUE RECO PAIRS - MASS RESOLUTION
    AliRsnMiniOutput * outres = task->CreateOutput(Form("Mres_%s", partname.Data()), output.Data(),"TRUE");
    outres->SetCutID(0, icut1);
    outres->SetCutID(1, icut2);
    outres->SetCharge(0, charge1[0]);
    outres->SetCharge(1, charge2[0]);
    outres->SetDaughter(0, d1);
    outres->SetDaughter(1, d2);
    outres->SetMotherPDG(pdgCode);
    outres->SetMotherMass(mass);
    outres->SetPairCuts(cutsPair);
    // axis X: invmass resolution
    outres->AddAxis(resID, 200, -0.01, 0.01);
    //axis Y: mother pt
    outres->AddAxis(ptID, 200, 0.0, 20.0);
    // axis Z: multiplicity
    if (enaMultSel) outres->AddAxis(multID, 100, 0.0, 100.0);

    //TRUE RECO PAIRS - rapidity
    AliRsnMiniOutput * outrap = task->CreateOutput(Form("trueRap_%s", partname.Data()), output.Data(),"TRUE");
    outrap->SetCutID(0, icut1);
    outrap->SetCutID(1, icut2);
    outrap->SetCharge(0, charge1[0]);
    outrap->SetCharge(1, charge2[0]);
    outrap->SetDaughter(0, d1);
    outrap->SetDaughter(1, d2);
    outrap->SetMotherPDG(pdgCode);
    outrap->SetMotherMass(mass);
    outrap->SetPairCuts(cutsPair);
    outrap->AddAxis(ptID, 160, 0.0, 16.0);
    outrap->AddAxis(yID,  120, -0.6, 0.6);
    outrap->AddAxis(etaID, 200, -1., 1.);

    
    //GENERATED PAIRS
    AliRsnMiniOutput * outm = task->CreateOutput(Form("motherf0_%s", partname.Data()), output.Data(),"MOTHER");
    outm->SetDaughter(0, d1);
    outm->SetDaughter(1, d2);
    outm->SetMotherPDG(pdgCode);
    outm->SetMotherMass(mass);
    outm->SetPairCuts(cutsPair);
    outm->AddAxis(imID, nbins, masslow, massup);
    outm->AddAxis(ptID, 200, 0.0, 20.0);
    if (enaMultSel) outm->AddAxis(multID, 100, 0.0, 100.0);

    //GENERATED PAIRS
    AliRsnMiniOutput * outmy = task->CreateOutput(Form("motherRap_%s", partname.Data()), output.Data(),"MOTHER");
    outmy->SetDaughter(0, d1);
    outmy->SetDaughter(1, d2);
    outmy->SetMotherPDG(pdgCode);
    outmy->SetMotherMass(mass);
    outmy->SetPairCuts(cutsPair);
    outmy->AddAxis(ptID, 160, 0.0, 16.0);
    outmy->AddAxis(yID,  120, -0.6, 0.6);
    outmy->AddAxis(etaID, 200, -1., 1.);

    //f2 GENERATED PAIRS
    AliRsnMiniOutput * outm2 = task->CreateOutput("motherf2", output.Data(),"MOTHER");
    outm2->SetDaughter(0, d1);
    outm2->SetDaughter(1, d2);
    outm2->SetMotherPDG(bgTemplatePDG[6]);
    outm2->SetMotherMass(bgTemplateM[6]);
    outm2->SetPairCuts(cutsPair);
    outm2->AddAxis(imID, nbins, masslow, massup);
    outm2->AddAxis(ptID, 200, 0.0, 20.0);
    if (enaMultSel) outm2->AddAxis(multID, 100, 0.0, 100.0);

    //f2 GENERATED PAIRS
    AliRsnMiniOutput * outmy2 = task->CreateOutput("motherf2Rap", output.Data(),"MOTHER");
    outmy2->SetDaughter(0, d1);
    outmy2->SetDaughter(1, d2);
    outmy2->SetMotherPDG(bgTemplatePDG[6]);
    outmy2->SetMotherMass(bgTemplateM[6]);
    outmy2->SetPairCuts(cutsPair);
    outmy2->AddAxis(ptID, 160, 0.0, 16.0);
    outmy2->AddAxis(yID,  120, -0.6, 0.6);
    outmy2->AddAxis(etaID, 200, -1., 1.);
  }

  return kTRUE;
}
