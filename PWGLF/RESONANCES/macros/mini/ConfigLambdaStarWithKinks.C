/***************************************************************************
Paraskevi.Ganoti@cern.ch - created on 2019 with the help of F. Bellini and 
other scripts in the repo.  
Configuration script for Lambda(1520) with Kaon identification via kinks analysis
****************************************************************************/

Bool_t ConfigLambdaStarWithKinks(AliRsnMiniAnalysisTask *task, 
		Bool_t                 isMC, 
		AliPIDResponse::EBeamType collSys = AliPIDResponse::kPBPB, //=0, kPPB=1, kPBPB=2
                AliRsnCutSet           *cutsPair,             //cuts on the pair
                AliRsnCutSet           *cutsPairY,             //cuts on the pair
		Bool_t                 enaMultSel = kTRUE,    //enable multiplicity axis
      		Float_t                masslow = 1.4,         //inv mass axis low edge 
		Float_t                massup = 2.2,           //inv mass axis upper edge 
		Int_t                  nbins = 800,           //inv mass axis n bins
		Float_t                nsigma = 3.0,          //nsigma of TPC PID cut
		Bool_t                 enableMonitor = kTRUE, //enable single track QA plots
		Float_t                kaonPIDCut=4.)        //nsigma V0 daughters
 
{
  //-----------------------
  //General 
  //-----------------------
  TString partname="LambdaStar";
//  RSNPID  d1 = AliRsnDaughter::kPion;
//  RSNPID  d2 = AliRsnDaughter::kLambda;
  Int_t   aodFilterBit = 0;

  //Additional options for monitoring plots
  TString monitorOpt = "NoSIGN";
  
  //-----------------------
  // CUTS
  //-----------------------
  //use default quality cuts std 2011 with crossed rows TPC
  Bool_t useCrossedRows = 1;

  AliRsnCutTrackQuality * fCutQuality = new AliRsnCutTrackQuality("CutQuality");
  fCutQuality->SetDefaults2011(useCrossedRows, kFALSE);
  
  AliRsnCutTOFMatch  *iCutTOFMatch     = new AliRsnCutTOFMatch("CutTOFMatch");
  AliRsnCutPIDNSigma *iCutTPCNSigma    = new AliRsnCutPIDNSigma("CutTPCNSigma", AliPID::kProton, AliRsnCutPIDNSigma::kTPC);//, AliRsnCutPIDNSigma::kTPCinnerP );
  AliRsnCutPIDNSigma *iCutTOFNSigma    = new AliRsnCutPIDNSigma("CutTOFNSigma", AliPID::kProton, AliRsnCutPIDNSigma::kTOF);//, AliRsnCutPIDNSigma::kP );
  AliRsnCutPIDNSigma *iCutTPCTOFNSigma = new AliRsnCutPIDNSigma("CutTPCTOFNSigma", AliPID::kProton, AliRsnCutPIDNSigma::kTPC);
  //for setting PID cuts without any selection in pT
  //iCutTPCNSigma->SinglePIDRange(fNsigmaTPC);
 
//for setting PID cuts in given pT ranges
  iCutTPCNSigma->AddPIDRange(4.00, 0.00, 0.25);
  iCutTPCNSigma->AddPIDRange(3.00, 0.25, 0.70);
  iCutTPCNSigma->AddPIDRange(2.00, 0.70, 1.10);
  //
  iCutTOFNSigma->AddPIDRange(0.00, 0.00, 0.80);  
  iCutTOFNSigma->AddPIDRange(3.00, 0.80, 1.e6);
  //
  iCutTPCTOFNSigma->SinglePIDRange(5.0);

 AliRsnCutSet * protonCutSet = new AliRsnCutSet("ProtonCutSet", AliRsnTarget::kDaughter);
 protonCutSet->AddCut(fCutQuality);
 protonCutSet->AddCut(iCutTOFMatch);
 protonCutSet->AddCut(iCutTPCNSigma);
 protonCutSet->AddCut(iCutTPCTOFNSigma);
 protonCutSet->AddCut(iCutTOFNSigma);

  // scheme:
      // quality & [ ( TOFmatch & TOF & TPCTOF ) || ( TPConly ) ]

 protonCutSet->SetCutScheme( Form(" %s & ( ( %s & %s & %s ) | ( %s ) )",
			 fCutQuality->GetName(),
			 iCutTOFMatch->GetName(), iCutTOFNSigma->GetName(), iCutTPCTOFNSigma->GetName(),
			 iCutTPCNSigma->GetName()) ) ;

  Int_t icutPr = task->AddTrackCuts(protonCutSet);
 
  //set daughter cuts
  Int_t icut1 = icutPr;

  // Kink

  Float_t crossedRows = 70;
  Float_t rowsbycluster = 0.8;
     
  AliESDtrackCuts *esdTrackCuts = new AliESDtrackCuts("kinkQuality");
  esdTrackCuts->SetEtaRange(-0.8,0.8);
  esdTrackCuts->SetPtRange(0.15,1.E10);
  esdTrackCuts->SetRequireTPCRefit();
  esdTrackCuts->SetRequireITSRefit();
  esdTrackCuts->SetAcceptKinkDaughters(1); //
  //esdTrackCuts->SetMinNCrossedRowsTPC(useCrossedRows);
  //esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(rowsbycluster);
  esdTrackCuts->SetMaxChi2PerClusterTPC(4);
  esdTrackCuts->SetMinNClustersTPC(30);
  esdTrackCuts->SetMaxDCAToVertexXY(0.3);
  esdTrackCuts->SetMaxDCAToVertexZ(2.5);
 
  AliRsnCutTrackQuality * fCutQualityKink = new AliRsnCutTrackQuality("KinkCutQuality");
  fCutQualityKink->SetESDtrackCuts(esdTrackCuts);
  AliRsnCutPIDkink *iCutkaonkink    = new AliRsnCutPIDkink("CutKaonKink", AliPID::kKaon, AliRsnCutPIDkink::kTPC);
  iCutkaonkink->AddPIDRange(2.5, 0.2, 15.);

  AliRsnCutSet * kaonCutSet = new AliRsnCutSet("kaonCutSet", AliRsnTarget::kDaughter);
  
  kaonCutSet->AddCut(fCutQualityKink); 
  kaonCutSet->AddCut(iCutkaonkink);

  //scheme (quality & kink Kaon)
  kaonCutSet->SetCutScheme( Form(" %s & %s", fCutQualityKink->GetName(), iCutkaonkink->GetName() ));

  // kaonCutSet->SetCutScheme( Form(" %s", fCutQualityKink->GetName() ));
   Int_t icutKaon = task->AddTrackCuts(kaonCutSet);
   Int_t icut2 = icutKaon;
 
  //QA plots 
  if (enableMonitor){
    Printf("======== Cut monitoring enabled");
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/AddMonitorOutput.C");
//    AddMonitorOutput(isMC, cutSetQuality->GetMonitorOutput(), monitorOpt.Data());
      AddMonitorOutput(isMC, protonCutSet->GetMonitorOutput(), monitorOpt.Data(), 0);
      AddMonitorOutput(isMC, kaonCutSet->GetMonitorOutput(), monitorOpt.Data(), 0);
  }  

  //-----------------------
  // OUTPUT
  //-----------------------
  /* invariant mass       */ Int_t imID   = task->CreateValue(AliRsnMiniValue::kInvMass, kFALSE);
  /* IM resolution        */ Int_t resID  = task->CreateValue(AliRsnMiniValue::kInvMassRes, kTRUE);
  /* transv. momentum     */ Int_t ptID   = task->CreateValue(AliRsnMiniValue::kPt, kFALSE);
  /* centrality           */ Int_t multID = task->CreateValue(AliRsnMiniValue::kMult, kFALSE);
  /* pseudorapidity       */ Int_t etaID  = task->CreateValue(AliRsnMiniValue::kEta, kFALSE);
  /* rapidity             */ Int_t yID    = task->CreateValue(AliRsnMiniValue::kY, kFALSE);
  /* 1st daughter pt      */ Int_t fdpt   = task->CreateValue(AliRsnMiniValue::kFirstDaughterPt, kFALSE);
  /* 2nd daughter pt      */ Int_t sdpt   = task->CreateValue(AliRsnMiniValue::kSecondDaughterPt, kFALSE);
  /* 1st daughter p       */ Int_t fdp    = task->CreateValue(AliRsnMiniValue::kFirstDaughterP, kFALSE);
  /* 2nd daughter p       */ Int_t sdp    = task->CreateValue(AliRsnMiniValue::kSecondDaughterP, kFALSE);
   /* gener. transv. mom.  */
  /* Asymmetry data       */ Int_t asymd   = task->CreateValue(AliRsnMiniValue::kAsym, kFALSE);
  /* Asymmetry            */ Int_t asym   = task->CreateValue(AliRsnMiniValue::kAsym, kTRUE);

  if (isMC) Int_t ptIDgen= task->CreateValue(AliRsnMiniValue::kPt, kTRUE);

  TString output[6] = {"HIST",      "HIST",          "HIST",      "HIST",          "HIST",        "HIST"          }; // or "SPARSE"
  TString name[6] =   {"Lambda",   "AntiLambda",   "Lambdamix", "AntiLambdamix",   "LambdaLSP",   "LambdamLSM"    };
  TString comp[6] =   {"PAIR",      "PAIR",         "MIX",         "MIX",           "PAIR",       "PAIR"          };
  Char_t charge1[6] = {  '+',         '-',          '+',            '-',             '+',          '-'            };
  Char_t charge2[6] = {  '-',         '+',          '-',            '+',             '+',          '-'            };
  Int_t ipdg [6] =    {  3124,        -3124,          3124,         -3124,           3124 ,        3124           };
  Double_t mass [6]=  { 1.5195,       1.5195,        1.5195,       1.5195,          1.5195,       1.5195          }; 
  Int_t cutID1 [6]= {icutPr, icutPr, icutPr, icutPr, icutPr, icutPr}; 
  Int_t cutID2 [6]= {icutKaon, icutKaon, icutKaon, icutKaon, icutKaon, icutKaon};

  //DATA 
  for (Int_t i = 0; i < 6; i++) {
    output[i] = "SPARSE";
    AliRsnMiniOutput *out = task->CreateOutput(Form("Lambda1520Kink_%s", name[i].Data()), output[i].Data(), comp[i].Data());
    out->SetCutID(0, cutID1[i]);
    out->SetCutID(1, cutID2[i]);
    out->SetDaughter(0, AliRsnDaughter::kProton);
    out->SetDaughter(1, AliRsnDaughter::kKaon);
    out->SetCharge(0, charge1[i]);
    out->SetCharge(1, charge2[i]);
    out->SetMotherPDG(ipdg[i]);
    out->SetMotherMass(mass[i]);
    out->SetPairCuts(cutsPair);
    // axis X: invmass 
    out->AddAxis(imID, nbins, masslow, massup);
    //axis Y: mother pt
    out->AddAxis(ptID, 150, 0.0, 15.0); //default use mother pt
    //axis Z: multiplicity
      if (enaMultSel) out->AddAxis(multID, 100, 0.0, 100.0);
  }

    AliRsnMiniOutput *outasmd = task->CreateOutput("hAsymmetryDAtaLambda", "HIST", "PAIR");
    outasmd->SetDaughter(0, AliRsnDaughter::kProton);
    outasmd->SetDaughter(1, AliRsnDaughter::kKaon);
    outasmd->SetCutID(0, icutPr);
    outasmd->SetCutID(1, icutKaon);
    outasmd->SetMotherPDG(3124);
    outasmd->SetMotherMass(1.5195);
    outasmd->SetPairCuts(cutsPair);
    outasmd->SetCharge(0, charge1[0]);
    outasmd->SetCharge(1, charge2[0]);
    outasmd->AddAxis(asymd, 200, -1.0, 1.0);
    if (enaMultSel) outasmd->AddAxis(multID, 100, 0.0, 100.0);   
 
  // AddMonitorOutput_LambdaPt(cutSetLs->GetMonitorOutput());

  if (isMC) {

    //TRUE RECO PAIRS - MASS
    output[0] = "SPARSE";
    // create output
    AliRsnMiniOutput *out = task->CreateOutput(Form("Lambda1520Kink_TrueMC_%s", name[0].Data()), output[0].Data(), "TRUE");
    // selection settings
    out->SetCutID(0, cutID1[0]);
    out->SetCutID(1, cutID2[0]);
    out->SetDaughter(0, AliRsnDaughter::kProton);
    out->SetDaughter(1, AliRsnDaughter::kKaon);
    out->SetCharge(0, charge1[0]);
    out->SetCharge(1, charge2[0]);
    out->SetMotherPDG(ipdg[0]);
    out->SetMotherMass(mass[0]);
    out->SetPairCuts(cutsPair);
    // binnings
    out->AddAxis(imID, 800, 1.4, 2.2);
    out->AddAxis(ptID, 150, 0.0, 15.0);
    
    if (enaMultSel) out->AddAxis(multID, 100, 0.0, 100.0);

    output[1] = "SPARSE";
    // create output
    AliRsnMiniOutput *out = task->CreateOutput(Form("ALambda1520Kink_TrueMC_%s", name[1].Data()), output[1].Data(), "TRUE");
    // selection settings
    out->SetCutID(0, cutID1[1]);
    out->SetCutID(1, cutID2[1]);
    out->SetDaughter(0, AliRsnDaughter::kProton);
    out->SetDaughter(1, AliRsnDaughter::kKaon);
    out->SetCharge(0, charge1[1]);
    out->SetCharge(1, charge2[1]);
    out->SetMotherPDG(ipdg[1]);
    out->SetMotherMass(mass[1]);
    out->SetPairCuts(cutsPair);
    // binnings
    out->AddAxis(imID, 800, 1.4, 2.2);
    out->AddAxis(ptID, 150, 0.0, 15.0);
    //out->AddAxis(lambdaDCA, 10, 0.0, 1.0);
    
    if (enaMultSel) out->AddAxis(multID, 100, 0.0, 100.0);
    
    //GENERATED PAIRS
    AliRsnMiniOutput *outasm = task->CreateOutput("hAsymmetryMCLambda", "HIST", "MOTHER");
    outasm->SetDaughter(0, AliRsnDaughter::kProton);
    outasm->SetDaughter(1, AliRsnDaughter::kKaon);
    outasm->SetMotherPDG(3124);
    outasm->SetMotherMass(1.5195);
    outasm->SetPairCuts(cutsPairY);
    outasm->AddAxis(asym, 200, -1.0, 1.0);
    if (enaMultSel) outasm->AddAxis(multID, 100, 0.0, 100.0);

    output[0] = "HIST";
    AliRsnMiniOutput * outm = task->CreateOutput(Form("motherLambda1520Kinkpt_%s", name[0].Data()), output[0].Data(),"MOTHER");
    outm->SetDaughter(0, AliRsnDaughter::kProton);
    outm->SetDaughter(1, AliRsnDaughter::kKaon);
    outm->SetMotherPDG(3124);
    outm->SetMotherMass(1.5195);
    outm->SetPairCuts(cutsPairY);
    outm->AddAxis(ptIDgen, 150, 0.0, 15.0);
    if (enaMultSel) outm->AddAxis(multID, 100, 0.0, 100.0);

    output[1] = "HIST";
    AliRsnMiniOutput * outm = task->CreateOutput(Form("Lambda1520Kinkpt_%s", name[1].Data()), output[1].Data(),"MOTHER");
    outm->SetDaughter(0, AliRsnDaughter::kProton);
    outm->SetDaughter(1, AliRsnDaughter::kKaon);
    outm->SetMotherPDG(-3124);
    outm->SetMotherMass(1.5195);
    outm->SetPairCuts(cutsPairY);
    outm->AddAxis(ptIDgen, 150, 0.0, 15.0);
    if (enaMultSel) outm->AddAxis(multID, 100, 0.0, 100.0);
 
  }

  return kTRUE;
}

