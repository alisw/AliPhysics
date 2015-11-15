/***************************************************************************
              fbellini@cern.ch - last modified on 02/07/2014

	      *** Configuration script for K*, anti-K* analysis of 2010 pp 7TeV datasets ***
This analysis task is used to extend the pT reach of the K* spectra published in Eur.
Phys. J. C72(2012)2183. 
****************************************************************************/

Bool_t ConfigRsnMC(AliRsnMiniAnalysisTask *task, 
		   Bool_t                 isMC, 
		   Bool_t                 isPP,
		   AliRsnCutSet           *cutsPair,
		   TString                partname="kStar",
		   Int_t                  pdgCode=313,
		   Float_t                mass = 0.89445,
		   Float_t                masslow = 0.7,
		   Float_t                massup = 1.1,
		   Int_t                  nbins = 400,
		   Char_t                 charge1 = '+',
		   Char_t                 charge2 = '-',
		   RSNPID                 d1 = AliRsnDaughter::kKaon,
		   RSNPID                 d2 = AliRsnDaughter::kPion,
		   Int_t                  aodFilterBit = 5,
		   Bool_t                 enableMonitor = kTRUE,
		   TString                monitorOpt = "NoSIGN")
{
  // -- Values ------------------------------------------------------------------------------------
  AliRsnMiniValue::EType yaxisVar = AliRsnMiniValue::kPt;
  /* invariant mass   */ Int_t imID   = task->CreateValue(AliRsnMiniValue::kInvMass, kFALSE);
  /* IM resolution    */ Int_t resID  = task->CreateValue(AliRsnMiniValue::kInvMassRes, kTRUE);
  /* transv. momentum */ Int_t ptID   = task->CreateValue(AliRsnMiniValue::kPt, kFALSE);
  /* centrality       */ Int_t centID = task->CreateValue(AliRsnMiniValue::kMult, kFALSE);
  /* pseudorapidity   */ Int_t etaID  = task->CreateValue(AliRsnMiniValue::kEta, kFALSE);
  /* rapidity         */ Int_t yID    = task->CreateValue(AliRsnMiniValue::kY, kFALSE);
  /* 1st daughter pt  */ Int_t fdpt   = task->CreateValue(AliRsnMiniValue::kFirstDaughterPt, kFALSE);
  /* 2nd daughter pt  */ Int_t sdpt   = task->CreateValue(AliRsnMiniValue::kSecondDaughterPt, kFALSE);
  /* 1st daughter p   */ Int_t fdp    = task->CreateValue(AliRsnMiniValue::kFirstDaughterP, kFALSE);
  /* 2nd daughter p   */ Int_t sdp    = task->CreateValue(AliRsnMiniValue::kSecondDaughterP, kFALSE);
  TString output = "SPARSE";

  // set daughter cuts
  //use default quality cuts std 2010 with crossed rows TPC
  Bool_t useCrossedRows = 1; Float_t nsigma = 3.0;
  AliRsnCutSetDaughterParticle * cutSetPi = new AliRsnCutSetDaughterParticle("cutPi", AliRsnCutSetDaughterParticle::kQualityStd2010, AliPID::kPion, nsigma, aodFilterBit, useCrossedRows);
  AliRsnCutSetDaughterParticle * cutSetKa = new AliRsnCutSetDaughterParticle("cutKa", AliRsnCutSetDaughterParticle::kQualityStd2010, AliPID::kKaon, nsigma, aodFilterBit, useCrossedRows);
  AliRsnCutSetDaughterParticle * cutSetPro = new AliRsnCutSetDaughterParticle("cutPro", AliRsnCutSetDaughterParticle::kQualityStd2010, AliPID::kProton, nsigma, aodFilterBit, useCrossedRows);
  Int_t icutPi = task->AddTrackCuts(cutSetPi);
  Int_t icutKa = task->AddTrackCuts(cutSetKa);
  Int_t icutPro = task->AddTrackCuts(cutSetPro);
  
  if (enableMonitor){
    Printf("======== Cut monitoring enabled");
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/AddMonitorOutput.C");
    AddMonitorOutput(isMC, cutSetPi->GetMonitorOutput(), monitorOpt.Data());
    AddMonitorOutput(isMC, cutSetKa->GetMonitorOutput()), monitorOpt.Data();
    AddMonitorOutput(isMC, cutSetPro->GetMonitorOutput()), monitorOpt.Data();
  }  

  Int_t icut1, icut2;
  if (d1==AliRsnDaughter::kProton) icut2 = icutPro;
  else if (d1==AliRsnDaughter::kKaon) icut2 = icutKa;
  else icut1 = icutPi;

  if (d2==AliRsnDaughter::kProton) icut2 = icutPro;
  else if (d2==AliRsnDaughter::kKaon) icut2 = icutKa;
  else icut2 = icutPi;

  //TRUE RECO PAIRS - MASS
  AliRsnMiniOutput * outtrue = task->CreateOutput(Form("true_%s", partname.Data()), output.Data(),"TRUE");
  outtrue->SetCutID(0, icut1);
  outtrue->SetCutID(1, icut2);
  outtrue->SetCharge(0, charge1);
  outtrue->SetCharge(1, charge2);
  outtrue->SetDaughter(0, d1);
  outtrue->SetDaughter(1, d2);
  outtrue->SetMotherPDG(pdgCode);
  outtrue->SetMotherMass(mass);
  outtrue->SetPairCuts(cutsPair);
  // axis X: invmass 
  outtrue->AddAxis(imID, nbins, masslow, massup);
  //axis Y: mother pt
  outtrue->AddAxis(ptID, 200, 0.0, 20.0); //default use mother pt
  // axis Z: centrality-multiplicity
  if (!isPP)
    outtrue->AddAxis(centID, 100, 0.0, 100.0);
  else 
    outtrue->AddAxis(centID, 400, 0.0, 400.0);

  //TRUE RECO PAIRS - MASS RESOLUTION
  AliRsnMiniOutput * outres = task->CreateOutput(Form("res_%s", partname.Data()), output.Data(),"TRUE");
  outres->SetCutID(0, icut1);
  outres->SetCutID(1, icut2);
  outres->SetCharge(0, charge1);
  outres->SetCharge(1, charge2);
  outres->SetDaughter(0, d1);
  outres->SetDaughter(1, d2);
  outres->SetMotherPDG(pdgCode);
  outres->SetMotherMass(mass);
  outres->SetPairCuts(cutsPair);
  // axis X: invmass resolution
  outres->AddAxis(resID, 200, -0.01, 0.01);
  //axis Y: mother pt
  outres->AddAxis(ptID, 200, 0.0, 20.0);
  // axis Z: centrality-multiplicity
  if (!isPP)
    outres->AddAxis(centID, 100, 0.0, 100.0);
  else 
    outres->AddAxis(centID, 400, 0.0, 400.0);
    
  //GENERATED PAIRS
  AliRsnMiniOutput * outm = task->CreateOutput(Form("mother_%s", partname.Data()), output.Data(),"MOTHER");
  outm->SetDaughter(0, d1);
  outm->SetDaughter(1, d2);
  outm->SetMotherPDG(pdgCode);
  outm->SetMotherMass(mass);
  outm->SetPairCuts(cutsPair);
  outm->AddAxis(imID, nbins, masslow, massup);
  outm->AddAxis(ptID, 200, 0.0, 20.0);
  if (!isPP){
    outm->AddAxis(centID, 100, 0.0, 100.0);
  }   else    { 
    outm->AddAxis(centID, 400, 0.0, 400.0);
  }

  return kTRUE;
}
