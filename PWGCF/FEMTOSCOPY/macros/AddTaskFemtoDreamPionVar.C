AliAnalysisTaskSE* AddTaskFemtoDreamPionVar(
    bool isMC=false, bool MCtemplatefit=false, float fSpherDown=0.7, float ptLow=0.5, float fdPhidEta=0.01,
    TString CentEst="kInt7", const char *cutVar = "0") {

  TString suffix = TString::Format("%s", cutVar);

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  if (!mgr)
  {
    printf("No analysis manager to connect to!\n");
    return nullptr;
  }
  if (!mgr->GetInputEventHandler()) {
    printf("This task requires an input event handler!\n");
    return nullptr;
  }

  AliFemtoDreamEventCuts *evtCuts=
      AliFemtoDreamEventCuts::StandardCutsRun2();
  //This sets the method we want to use to clean up events with negative or too
  //low multiplicity. Usually you use the matching multiplicity estiamtor in your
  //event collection
  // Not mention in AN oder Indico
  evtCuts->CleanUpMult(false,false,false,true);
  evtCuts->SetZVtxPosition(-10., 10.);
  // Only use those events where more than two primary tracks with |eta|<0.8 and pT>0.5 GeV/c see AN

  AliFemtoDreamTrackCuts *fTrackCutsPosPion=new AliFemtoDreamTrackCuts();
  AliFemtoDreamTrackCuts *fTrackCutsNegPion=new AliFemtoDreamTrackCuts();

  if (suffix == "1") {//pT low
  fSpherDown = 0.7;
  //Track Cuts are defined here
  //positive pions
  fdPhidEta=0.012;
  fTrackCutsPosPion->SetIsMonteCarlo(isMC);
  fTrackCutsPosPion->SetCutCharge(1);
  fTrackCutsPosPion->SetFilterBit(96);
  fTrackCutsPosPion->SetPtRange(0.11, 4.3);
  fTrackCutsPosPion->SetEtaRange(-0.82, 0.82);
  fTrackCutsPosPion->SetNClsTPC(78);
  fTrackCutsPosPion->SetDCAReCalculation(true);
  fTrackCutsPosPion->SetDCAVtxZ(0.33);
  fTrackCutsPosPion->SetDCAVtxXY(0.33);
  fTrackCutsPosPion->SetPID(AliPID::kPion, 0.48);
  fTrackCutsPosPion->SetRejLowPtPionsTOF(false);
  fTrackCutsPosPion->SetMinimalBooking(false);

  //The same things for negative pions
  fTrackCutsNegPion->SetIsMonteCarlo(isMC);
  fTrackCutsNegPion->SetCutCharge(-1);
  fTrackCutsNegPion->SetFilterBit(96);
  fTrackCutsNegPion->SetPtRange(0.11, 4.3);
  fTrackCutsNegPion->SetEtaRange(-0.82, 0.82);
  fTrackCutsNegPion->SetNClsTPC(78);
  fTrackCutsNegPion->SetDCAReCalculation(true);
  fTrackCutsNegPion->SetDCAVtxZ(0.33);
  fTrackCutsNegPion->SetDCAVtxXY(0.33);
  fTrackCutsNegPion->SetPID(AliPID::kPion, 0.48);
  fTrackCutsNegPion->SetRejLowPtPionsTOF(false);
  fTrackCutsNegPion->SetMinimalBooking(false);
  } else if (suffix == "2") {//pT high
  fSpherDown = 0.71;
  //Track Cuts are defined here
  //positive pions
  fdPhidEta=0.011;
  fTrackCutsPosPion->SetIsMonteCarlo(isMC);
  fTrackCutsPosPion->SetCutCharge(1);
  fTrackCutsPosPion->SetFilterBit(96);
  fTrackCutsPosPion->SetPtRange(0.12, 5.0);
  fTrackCutsPosPion->SetEtaRange(-0.8, 0.8);
  fTrackCutsPosPion->SetNClsTPC(77);
  fTrackCutsPosPion->SetDCAReCalculation(true);
  fTrackCutsPosPion->SetDCAVtxZ(0.26);
  fTrackCutsPosPion->SetDCAVtxXY(0.26);
  fTrackCutsPosPion->SetPID(AliPID::kPion, 0.46);
  fTrackCutsPosPion->SetRejLowPtPionsTOF(false);
  fTrackCutsPosPion->SetMinimalBooking(false);

  //The same things for negative pions
  fTrackCutsNegPion->SetIsMonteCarlo(isMC);
  fTrackCutsNegPion->SetCutCharge(-1);
  fTrackCutsNegPion->SetFilterBit(96);
  fTrackCutsNegPion->SetPtRange(0.14, 5.0);
  fTrackCutsNegPion->SetEtaRange(-0.8, 0.8);
  fTrackCutsNegPion->SetNClsTPC(77);
  fTrackCutsNegPion->SetDCAReCalculation(true);
  fTrackCutsNegPion->SetDCAVtxZ(0.26);
  fTrackCutsNegPion->SetDCAVtxXY(0.26);
  fTrackCutsNegPion->SetPID(AliPID::kPion, 0.46);
  fTrackCutsNegPion->SetRejLowPtPionsTOF(false);
  fTrackCutsNegPion->SetMinimalBooking(false);
  } else if ( suffix == "3" ) {//eta low
  fSpherDown = 0.73;
  //Track Cuts are defined here
  //positive pions
  fdPhidEta=0.010;
  fTrackCutsPosPion->SetIsMonteCarlo(isMC);
  fTrackCutsPosPion->SetCutCharge(1);
  fTrackCutsPosPion->SetFilterBit(96);
  fTrackCutsPosPion->SetPtRange(0.18, 4.3);
  fTrackCutsPosPion->SetEtaRange(-0.79, 0.79);
  fTrackCutsPosPion->SetNClsTPC(80);
  fTrackCutsPosPion->SetDCAReCalculation(true);
  fTrackCutsPosPion->SetDCAVtxZ(0.25);
  fTrackCutsPosPion->SetDCAVtxXY(0.25);
  fTrackCutsPosPion->SetPID(AliPID::kPion, 0.5);
  fTrackCutsPosPion->SetRejLowPtPionsTOF(false);
  fTrackCutsPosPion->SetMinimalBooking(false);

  //The same things for negative pions
  fTrackCutsNegPion->SetIsMonteCarlo(isMC);
  fTrackCutsNegPion->SetCutCharge(-1);
  fTrackCutsNegPion->SetFilterBit(96);
  fTrackCutsNegPion->SetPtRange(0.18, 4.3);
  fTrackCutsNegPion->SetEtaRange(-0.79, 0.79);
  fTrackCutsNegPion->SetNClsTPC(80);
  fTrackCutsNegPion->SetDCAReCalculation(true);
  fTrackCutsNegPion->SetDCAVtxZ(0.25);
  fTrackCutsNegPion->SetDCAVtxXY(0.25);
  fTrackCutsNegPion->SetPID(AliPID::kPion, 0.5);
  fTrackCutsNegPion->SetRejLowPtPionsTOF(false);
  fTrackCutsNegPion->SetMinimalBooking(false);
  } else if ( suffix == "4" ) {//eta high
  fSpherDown = 0.69;
  //Track Cuts are defined here
  //positive pions
  fdPhidEta=0.009;
  fTrackCutsPosPion->SetIsMonteCarlo(isMC);
  fTrackCutsPosPion->SetCutCharge(1);
  fTrackCutsPosPion->SetFilterBit(96);
  fTrackCutsPosPion->SetPtRange(0.11, 4.0);
  fTrackCutsPosPion->SetEtaRange(-0.9, 0.9);
  fTrackCutsPosPion->SetNClsTPC(90);
  fTrackCutsPosPion->SetDCAReCalculation(true);
  fTrackCutsPosPion->SetDCAVtxZ(0.33);
  fTrackCutsPosPion->SetDCAVtxXY(0.33);
  fTrackCutsPosPion->SetPID(AliPID::kPion, 0.54);
  fTrackCutsPosPion->SetRejLowPtPionsTOF(false);
  fTrackCutsPosPion->SetMinimalBooking(false);

  //The same things for negative pions
  fTrackCutsNegPion->SetIsMonteCarlo(isMC);
  fTrackCutsNegPion->SetCutCharge(-1);
  fTrackCutsNegPion->SetFilterBit(96);
  fTrackCutsNegPion->SetPtRange(0.11, 4.0);
  fTrackCutsNegPion->SetEtaRange(-0.9, 0.9);
  fTrackCutsNegPion->SetNClsTPC(90);
  fTrackCutsNegPion->SetDCAReCalculation(true);
  fTrackCutsNegPion->SetDCAVtxZ(0.33);
  fTrackCutsNegPion->SetDCAVtxXY(0.33);
  fTrackCutsNegPion->SetPID(AliPID::kPion, 0.54);
  fTrackCutsNegPion->SetRejLowPtPionsTOF(false);
  fTrackCutsNegPion->SetMinimalBooking(false);
  } else if ( suffix == "5" ) {//eta high conservative
  fSpherDown = 0.7;
  //Track Cuts are defined here
  //positive pions
  fdPhidEta=0.008;
  fTrackCutsPosPion->SetIsMonteCarlo(isMC);
  fTrackCutsPosPion->SetCutCharge(1);
  fTrackCutsPosPion->SetFilterBit(96);
  fTrackCutsPosPion->SetPtRange(0.12, 4.8);
  fTrackCutsPosPion->SetEtaRange(-0.85, 0.85);
  fTrackCutsPosPion->SetNClsTPC(88);
  fTrackCutsPosPion->SetDCAReCalculation(true);
  fTrackCutsPosPion->SetDCAVtxZ(0.26);
  fTrackCutsPosPion->SetDCAVtxXY(0.26);
  fTrackCutsPosPion->SetPID(AliPID::kPion, 0.45);
  fTrackCutsPosPion->SetRejLowPtPionsTOF(false);
  fTrackCutsPosPion->SetMinimalBooking(false);

  //The same things for negative pions
  fTrackCutsNegPion->SetIsMonteCarlo(isMC);
  fTrackCutsNegPion->SetCutCharge(-1);
  fTrackCutsNegPion->SetFilterBit(96);
  fTrackCutsNegPion->SetPtRange(0.12, 4.8);
  fTrackCutsNegPion->SetEtaRange(-0.85, 0.85);
  fTrackCutsNegPion->SetNClsTPC(88);
  fTrackCutsNegPion->SetDCAReCalculation(true);
  fTrackCutsNegPion->SetDCAVtxZ(0.26);
  fTrackCutsNegPion->SetDCAVtxXY(0.26);
  fTrackCutsNegPion->SetPID(AliPID::kPion, 0.45);
  fTrackCutsNegPion->SetRejLowPtPionsTOF(false);
  fTrackCutsNegPion->SetMinimalBooking(false);
  } else if ( suffix == "6" ) {//eta low conservative
  fSpherDown = 0.68;
  //Track Cuts are defined here
  //positive pions
  fdPhidEta=0.009;
  fTrackCutsPosPion->SetIsMonteCarlo(isMC);
  fTrackCutsPosPion->SetCutCharge(1);
  fTrackCutsPosPion->SetFilterBit(96);
  fTrackCutsPosPion->SetPtRange(0.14, 4.1);
  fTrackCutsPosPion->SetEtaRange(-0.77, 0.77);
  fTrackCutsPosPion->SetNClsTPC(78);
  fTrackCutsPosPion->SetDCAReCalculation(true);
  fTrackCutsPosPion->SetDCAVtxZ(0.31);
  fTrackCutsPosPion->SetDCAVtxXY(0.31);
  fTrackCutsPosPion->SetPID(AliPID::kPion, 0.54);
  fTrackCutsPosPion->SetRejLowPtPionsTOF(false);
  fTrackCutsPosPion->SetMinimalBooking(false);

  //The same things for negative pions
  fTrackCutsNegPion->SetIsMonteCarlo(isMC);
  fTrackCutsNegPion->SetCutCharge(-1);
  fTrackCutsNegPion->SetFilterBit(96);
  fTrackCutsNegPion->SetPtRange(0.14, 4.1);
  fTrackCutsNegPion->SetEtaRange(-0.77, 0.77);
  fTrackCutsNegPion->SetNClsTPC(78);
  fTrackCutsNegPion->SetDCAReCalculation(true);
  fTrackCutsNegPion->SetDCAVtxZ(0.31);
  fTrackCutsNegPion->SetDCAVtxXY(0.31);
  fTrackCutsNegPion->SetPID(AliPID::kPion, 0.54);
  fTrackCutsNegPion->SetRejLowPtPionsTOF(false);
  fTrackCutsNegPion->SetMinimalBooking(false);
  } else if ( suffix == "7" ) {//#TPC Cls low
  fSpherDown = 0.67;
  //Track Cuts are defined here
  //positive pions
  fdPhidEta=0.010;
  fTrackCutsPosPion->SetIsMonteCarlo(isMC);
  fTrackCutsPosPion->SetCutCharge(1);
  fTrackCutsPosPion->SetFilterBit(96);
  fTrackCutsPosPion->SetPtRange(0.11, 4.7);
  fTrackCutsPosPion->SetEtaRange(-0.87, 0.87);
  fTrackCutsPosPion->SetNClsTPC(76);
  fTrackCutsPosPion->SetDCAReCalculation(true);
  fTrackCutsPosPion->SetDCAVtxZ(0.35);
  fTrackCutsPosPion->SetDCAVtxXY(0.35);
  fTrackCutsPosPion->SetPID(AliPID::kPion, 0.46);
  fTrackCutsPosPion->SetRejLowPtPionsTOF(false);
  fTrackCutsPosPion->SetMinimalBooking(false);

  //The same things for negative pions
  fTrackCutsNegPion->SetIsMonteCarlo(isMC);
  fTrackCutsNegPion->SetCutCharge(-1);
  fTrackCutsNegPion->SetFilterBit(96);
  fTrackCutsNegPion->SetPtRange(0.11, 4.7);
  fTrackCutsNegPion->SetEtaRange(-0.87, 0.87);
  fTrackCutsNegPion->SetNClsTPC(76);
  fTrackCutsNegPion->SetDCAReCalculation(true);
  fTrackCutsNegPion->SetDCAVtxZ(0.35);
  fTrackCutsNegPion->SetDCAVtxXY(0.35);
  fTrackCutsNegPion->SetPID(AliPID::kPion, 0.46);
  fTrackCutsNegPion->SetRejLowPtPionsTOF(false);
  fTrackCutsNegPion->SetMinimalBooking(false);
  } else if ( suffix == "8" ) {//#TPC Cls high conservative
  fSpherDown = 0.67;
  //Track Cuts are defined here
  //positive pions
  fdPhidEta=0.011;
  fTrackCutsPosPion->SetIsMonteCarlo(isMC);
  fTrackCutsPosPion->SetCutCharge(1);
  fTrackCutsPosPion->SetFilterBit(96);
  fTrackCutsPosPion->SetPtRange(0.14, 4.0);
  fTrackCutsPosPion->SetEtaRange(-0.76, 0.76);
  fTrackCutsPosPion->SetNClsTPC(90);
  fTrackCutsPosPion->SetDCAReCalculation(true);
  fTrackCutsPosPion->SetDCAVtxZ(0.3);
  fTrackCutsPosPion->SetDCAVtxXY(0.3);
  fTrackCutsPosPion->SetPID(AliPID::kPion, 0.47);
  fTrackCutsPosPion->SetRejLowPtPionsTOF(false);
  fTrackCutsPosPion->SetMinimalBooking(false);

  //The same things for negative pions
  fTrackCutsNegPion->SetIsMonteCarlo(isMC);
  fTrackCutsNegPion->SetCutCharge(-1);
  fTrackCutsNegPion->SetFilterBit(96);
  fTrackCutsNegPion->SetPtRange(0.14, 4.0);
  fTrackCutsNegPion->SetEtaRange(-0.76, 0.76);
  fTrackCutsNegPion->SetNClsTPC(90);
  fTrackCutsNegPion->SetDCAReCalculation(true);
  fTrackCutsNegPion->SetDCAVtxZ(0.3);
  fTrackCutsNegPion->SetDCAVtxXY(0.3);
  fTrackCutsNegPion->SetPID(AliPID::kPion, 0.47);
  fTrackCutsNegPion->SetRejLowPtPionsTOF(false);
  fTrackCutsNegPion->SetMinimalBooking(false);
  } else if ( suffix == "9" ) {//#TPC Cls high
  fSpherDown = 0.68;
  //Track Cuts are defined here
  //positive pions
  fdPhidEta=0.012;
  fTrackCutsPosPion->SetIsMonteCarlo(isMC);
  fTrackCutsPosPion->SetCutCharge(1);
  fTrackCutsPosPion->SetFilterBit(96);
  fTrackCutsPosPion->SetPtRange(0.12, 4.2);
  fTrackCutsPosPion->SetEtaRange(-0.83, 0.83);
  fTrackCutsPosPion->SetNClsTPC(80);
  fTrackCutsPosPion->SetDCAReCalculation(true);
  fTrackCutsPosPion->SetDCAVtxZ(0.27);
  fTrackCutsPosPion->SetDCAVtxXY(0.27);
  fTrackCutsPosPion->SetPID(AliPID::kPion, 0.54);
  fTrackCutsPosPion->SetRejLowPtPionsTOF(false);
  fTrackCutsPosPion->SetMinimalBooking(false);

  //The same things for negative pions
  fTrackCutsNegPion->SetIsMonteCarlo(isMC);
  fTrackCutsNegPion->SetCutCharge(-1);
  fTrackCutsNegPion->SetFilterBit(96);
  fTrackCutsNegPion->SetPtRange(0.12, 4.2);
  fTrackCutsNegPion->SetEtaRange(-0.83, 0.83);
  fTrackCutsNegPion->SetNClsTPC(80);
  fTrackCutsNegPion->SetDCAReCalculation(true);
  fTrackCutsNegPion->SetDCAVtxZ(0.27);
  fTrackCutsNegPion->SetDCAVtxXY(0.27);
  fTrackCutsNegPion->SetPID(AliPID::kPion, 0.54);
  fTrackCutsNegPion->SetRejLowPtPionsTOF(false);
  fTrackCutsNegPion->SetMinimalBooking(false);
  } else if ( suffix == "10" ) {//#TPC Cls low conservative
  fSpherDown = 0.69;
  //Track Cuts are defined here
  //positive pions
  fdPhidEta=0.008;
  fTrackCutsPosPion->SetIsMonteCarlo(isMC);
  fTrackCutsPosPion->SetCutCharge(1);
  fTrackCutsPosPion->SetFilterBit(96);
  fTrackCutsPosPion->SetPtRange(0.13, 4.9);
  fTrackCutsPosPion->SetEtaRange(-0.88, 0.88);
  fTrackCutsPosPion->SetNClsTPC(75);
  fTrackCutsPosPion->SetDCAReCalculation(true);
  fTrackCutsPosPion->SetDCAVtxZ(0.31);
  fTrackCutsPosPion->SetDCAVtxXY(0.31);
  fTrackCutsPosPion->SetPID(AliPID::kPion, 0.46);
  fTrackCutsPosPion->SetRejLowPtPionsTOF(false);
  fTrackCutsPosPion->SetMinimalBooking(false);

  //The same things for negative pions
  fTrackCutsNegPion->SetIsMonteCarlo(isMC);
  fTrackCutsNegPion->SetCutCharge(-1);
  fTrackCutsNegPion->SetFilterBit(96);
  fTrackCutsNegPion->SetPtRange(0.13, 4.9);
  fTrackCutsNegPion->SetEtaRange(-0.88, 0.88);
  fTrackCutsNegPion->SetNClsTPC(75);
  fTrackCutsNegPion->SetDCAReCalculation(true);
  fTrackCutsNegPion->SetDCAVtxZ(0.31);
  fTrackCutsNegPion->SetDCAVtxXY(0.31);
  fTrackCutsNegPion->SetPID(AliPID::kPion, 0.46);
  fTrackCutsNegPion->SetRejLowPtPionsTOF(false);
  fTrackCutsNegPion->SetMinimalBooking(false);
  } else if ( suffix == "11" ) {//DCA XYZ low conservative
  fSpherDown = 0.7;
  //Track Cuts are defined here
  //positive pions
  fdPhidEta=0.009;
  fTrackCutsPosPion->SetIsMonteCarlo(isMC);
  fTrackCutsPosPion->SetCutCharge(1);
  fTrackCutsPosPion->SetFilterBit(96);
  fTrackCutsPosPion->SetPtRange(0.14, 4.0);
  fTrackCutsPosPion->SetEtaRange(-0.79, 0.79);
  fTrackCutsPosPion->SetNClsTPC(77);
  fTrackCutsPosPion->SetDCAReCalculation(true);
  fTrackCutsPosPion->SetDCAVtxZ(0.28);
  fTrackCutsPosPion->SetDCAVtxXY(0.28);
  fTrackCutsPosPion->SetPID(AliPID::kPion, 0.5);
  fTrackCutsPosPion->SetRejLowPtPionsTOF(false);
  fTrackCutsPosPion->SetMinimalBooking(false);

  //The same things for negative pions
  fTrackCutsNegPion->SetIsMonteCarlo(isMC);
  fTrackCutsNegPion->SetCutCharge(-1);
  fTrackCutsNegPion->SetFilterBit(96);
  fTrackCutsNegPion->SetPtRange(0.14, 4.0);
  fTrackCutsNegPion->SetEtaRange(-0.79, 0.79);
  fTrackCutsNegPion->SetNClsTPC(77);
  fTrackCutsNegPion->SetDCAReCalculation(true);
  fTrackCutsNegPion->SetDCAVtxZ(0.28);
  fTrackCutsNegPion->SetDCAVtxXY(0.28);
  fTrackCutsNegPion->SetPID(AliPID::kPion, 0.5);
  fTrackCutsNegPion->SetRejLowPtPionsTOF(false);
  fTrackCutsNegPion->SetMinimalBooking(false);
  } else if ( suffix == "12" ) {//DCA XYZ high conservative
  fSpherDown = 0.71;
  //Track Cuts are defined here
  //positive pions
  fdPhidEta=0.010;
  fTrackCutsPosPion->SetIsMonteCarlo(isMC);
  fTrackCutsPosPion->SetCutCharge(1);
  fTrackCutsPosPion->SetFilterBit(96);
  fTrackCutsPosPion->SetPtRange(0.13, 4.0);
  fTrackCutsPosPion->SetEtaRange(-0.8, 0.8);
  fTrackCutsPosPion->SetNClsTPC(76);
  fTrackCutsPosPion->SetDCAReCalculation(true);
  fTrackCutsPosPion->SetDCAVtxZ(0.32);
  fTrackCutsPosPion->SetDCAVtxXY(0.32);
  fTrackCutsPosPion->SetPID(AliPID::kPion, 0.5);
  fTrackCutsPosPion->SetRejLowPtPionsTOF(false);
  fTrackCutsPosPion->SetMinimalBooking(false);

  //The same things for negative pions
  fTrackCutsNegPion->SetIsMonteCarlo(isMC);
  fTrackCutsNegPion->SetCutCharge(-1);
  fTrackCutsNegPion->SetFilterBit(96);
  fTrackCutsNegPion->SetPtRange(0.13, 4.0);
  fTrackCutsNegPion->SetEtaRange(-0.8, 0.8);
  fTrackCutsNegPion->SetNClsTPC(76);
  fTrackCutsNegPion->SetDCAReCalculation(true);
  fTrackCutsNegPion->SetDCAVtxZ(0.32);
  fTrackCutsNegPion->SetDCAVtxXY(0.32);
  fTrackCutsNegPion->SetPID(AliPID::kPion, 0.5);
  fTrackCutsNegPion->SetRejLowPtPionsTOF(false);
  fTrackCutsNegPion->SetMinimalBooking(false);
  } else if ( suffix == "13" ) {//DCA XYZ high
  fSpherDown = 0.72;
  //Track Cuts are defined here
  //positive pions
  fdPhidEta=0.011;
  fTrackCutsPosPion->SetIsMonteCarlo(isMC);
  fTrackCutsPosPion->SetCutCharge(1);
  fTrackCutsPosPion->SetFilterBit(96);
  fTrackCutsPosPion->SetPtRange(0.14, 4.5);
  fTrackCutsPosPion->SetEtaRange(-0.84, 0.84);
  fTrackCutsPosPion->SetNClsTPC(75);
  fTrackCutsPosPion->SetDCAReCalculation(true);
  fTrackCutsPosPion->SetDCAVtxZ(0.35);
  fTrackCutsPosPion->SetDCAVtxXY(0.35);
  fTrackCutsPosPion->SetPID(AliPID::kPion, 0.53);
  fTrackCutsPosPion->SetRejLowPtPionsTOF(false);
  fTrackCutsPosPion->SetMinimalBooking(false);

  //The same things for negative pions
  fTrackCutsNegPion->SetIsMonteCarlo(isMC);
  fTrackCutsNegPion->SetCutCharge(-1);
  fTrackCutsNegPion->SetFilterBit(96);
  fTrackCutsNegPion->SetPtRange(0.14, 4.5);
  fTrackCutsNegPion->SetEtaRange(-0.84, 0.84);
  fTrackCutsNegPion->SetNClsTPC(75);
  fTrackCutsNegPion->SetDCAReCalculation(true);
  fTrackCutsNegPion->SetDCAVtxZ(0.35);
  fTrackCutsNegPion->SetDCAVtxXY(0.35);
  fTrackCutsNegPion->SetPID(AliPID::kPion, 0.53);
  fTrackCutsNegPion->SetRejLowPtPionsTOF(false);
  fTrackCutsNegPion->SetMinimalBooking(false);
  } else if ( suffix == "14" ) {//DCA XYZ low
  fSpherDown = 0.73;
  //Track Cuts are defined here
  //positive pions
  fdPhidEta=0.012;
  fTrackCutsPosPion->SetIsMonteCarlo(isMC);
  fTrackCutsPosPion->SetCutCharge(1);
  fTrackCutsPosPion->SetFilterBit(96);
  fTrackCutsPosPion->SetPtRange(0.11, 4.1);
  fTrackCutsPosPion->SetEtaRange(-0.78, 0.78);
  fTrackCutsPosPion->SetNClsTPC(76);
  fTrackCutsPosPion->SetDCAReCalculation(true);
  fTrackCutsPosPion->SetDCAVtxZ(0.25);
  fTrackCutsPosPion->SetDCAVtxXY(0.25);
  fTrackCutsPosPion->SetPID(AliPID::kPion, 0.52);
  fTrackCutsPosPion->SetRejLowPtPionsTOF(false);
  fTrackCutsPosPion->SetMinimalBooking(false);

  //The same things for negative pions
  fTrackCutsNegPion->SetIsMonteCarlo(isMC);
  fTrackCutsNegPion->SetCutCharge(-1);
  fTrackCutsNegPion->SetFilterBit(96);
  fTrackCutsNegPion->SetPtRange(0.11, 4.1);
  fTrackCutsNegPion->SetEtaRange(-0.78, 0.78);
  fTrackCutsNegPion->SetNClsTPC(76);
  fTrackCutsNegPion->SetDCAReCalculation(true);
  fTrackCutsNegPion->SetDCAVtxZ(0.25);
  fTrackCutsNegPion->SetDCAVtxXY(0.25);
  fTrackCutsNegPion->SetPID(AliPID::kPion, 0.52);
  fTrackCutsNegPion->SetRejLowPtPionsTOF(false);
  fTrackCutsNegPion->SetMinimalBooking(false);
  } else if ( suffix == "15" ) {//Sphericity low
  fSpherDown = 0.72;
  //Track Cuts are defined here
  //positive pions
  fdPhidEta=0.012;
  fTrackCutsPosPion->SetIsMonteCarlo(isMC);
  fTrackCutsPosPion->SetCutCharge(1);
  fTrackCutsPosPion->SetFilterBit(96);
  fTrackCutsPosPion->SetPtRange(0.13, 4.6);
  fTrackCutsPosPion->SetEtaRange(-0.74, 0.74);
  fTrackCutsPosPion->SetNClsTPC(86);
  fTrackCutsPosPion->SetDCAReCalculation(true);
  fTrackCutsPosPion->SetDCAVtxZ(0.25);
  fTrackCutsPosPion->SetDCAVtxXY(0.25);
  fTrackCutsPosPion->SetPID(AliPID::kPion, 0.47);
  fTrackCutsPosPion->SetRejLowPtPionsTOF(false);
  fTrackCutsPosPion->SetMinimalBooking(false);

  //The same things for negative pions
  fTrackCutsNegPion->SetIsMonteCarlo(isMC);
  fTrackCutsNegPion->SetCutCharge(-1);
  fTrackCutsNegPion->SetFilterBit(96);
  fTrackCutsNegPion->SetPtRange(0.13, 4.6);
  fTrackCutsNegPion->SetEtaRange(-0.74, 0.74);
  fTrackCutsNegPion->SetNClsTPC(86);
  fTrackCutsNegPion->SetDCAReCalculation(true);
  fTrackCutsNegPion->SetDCAVtxZ(0.25);
  fTrackCutsNegPion->SetDCAVtxXY(0.25);
  fTrackCutsNegPion->SetPID(AliPID::kPion, 0.47);
  fTrackCutsNegPion->SetRejLowPtPionsTOF(false);
  fTrackCutsNegPion->SetMinimalBooking(false);
  } else if ( suffix == "16" ) {//Sphericity high
  fSpherDown = 0.71;
  //Track Cuts are defined here
  //positive pions
  fdPhidEta=0.011;
  fTrackCutsPosPion->SetIsMonteCarlo(isMC);
  fTrackCutsPosPion->SetCutCharge(1);
  fTrackCutsPosPion->SetFilterBit(96);
  fTrackCutsPosPion->SetPtRange(0.13, 4.2);
  fTrackCutsPosPion->SetEtaRange(-0.78, 0.78);
  fTrackCutsPosPion->SetNClsTPC(90);
  fTrackCutsPosPion->SetDCAReCalculation(true);
  fTrackCutsPosPion->SetDCAVtxZ(0.33);
  fTrackCutsPosPion->SetDCAVtxXY(0.33);
  fTrackCutsPosPion->SetPID(AliPID::kPion, 0.45);
  fTrackCutsPosPion->SetRejLowPtPionsTOF(false);
  fTrackCutsPosPion->SetMinimalBooking(false);

  //The same things for negative pions
  fTrackCutsNegPion->SetIsMonteCarlo(isMC);
  fTrackCutsNegPion->SetCutCharge(-1);
  fTrackCutsNegPion->SetFilterBit(96);
  fTrackCutsNegPion->SetPtRange(0.13, 4.2);
  fTrackCutsNegPion->SetEtaRange(-0.78, 0.78);
  fTrackCutsNegPion->SetNClsTPC(90);
  fTrackCutsNegPion->SetDCAReCalculation(true);
  fTrackCutsNegPion->SetDCAVtxZ(0.33);
  fTrackCutsNegPion->SetDCAVtxXY(0.33);
  fTrackCutsNegPion->SetPID(AliPID::kPion, 0.45);
  fTrackCutsNegPion->SetRejLowPtPionsTOF(false);
  fTrackCutsNegPion->SetMinimalBooking(false);
  } else if ( suffix == "17" ) {//PID low
  fSpherDown = 0.7;
  //Track Cuts are defined here
  //positive pions
  fdPhidEta=0.011;
  fTrackCutsPosPion->SetIsMonteCarlo(isMC);
  fTrackCutsPosPion->SetCutCharge(1);
  fTrackCutsPosPion->SetFilterBit(96);
  fTrackCutsPosPion->SetPtRange(0.12, 4.5);
  fTrackCutsPosPion->SetEtaRange(-0.76, 0.76);
  fTrackCutsPosPion->SetNClsTPC(80);
  fTrackCutsPosPion->SetDCAReCalculation(true);
  fTrackCutsPosPion->SetDCAVtxZ(0.3);
  fTrackCutsPosPion->SetDCAVtxXY(0.3);
  fTrackCutsPosPion->SetPID(AliPID::kPion, 0.46);
  fTrackCutsPosPion->SetRejLowPtPionsTOF(false);
  fTrackCutsPosPion->SetMinimalBooking(false);

  //The same things for negative pions
  fTrackCutsNegPion->SetIsMonteCarlo(isMC);
  fTrackCutsNegPion->SetCutCharge(-1);
  fTrackCutsNegPion->SetFilterBit(96);
  fTrackCutsNegPion->SetPtRange(0.12, 4.5);
  fTrackCutsNegPion->SetEtaRange(-0.76, 0.76);
  fTrackCutsNegPion->SetNClsTPC(80);
  fTrackCutsNegPion->SetDCAReCalculation(true);
  fTrackCutsNegPion->SetDCAVtxZ(0.26);
  fTrackCutsNegPion->SetDCAVtxXY(0.26);
  fTrackCutsNegPion->SetPID(AliPID::kPion, 0.46);
  fTrackCutsNegPion->SetRejLowPtPionsTOF(false);
  fTrackCutsNegPion->SetMinimalBooking(false);
  } else if ( suffix == "18" ) {//PID high
  fSpherDown = 0.69;
  //Track Cuts are defined here
  //positive pions
  fdPhidEta=0.01;
  fTrackCutsPosPion->SetIsMonteCarlo(isMC);
  fTrackCutsPosPion->SetCutCharge(1);
  fTrackCutsPosPion->SetFilterBit(96);
  fTrackCutsPosPion->SetPtRange(0.14, 4.7);
  fTrackCutsPosPion->SetEtaRange(-0.86, 0.86);
  fTrackCutsPosPion->SetNClsTPC(77);
  fTrackCutsPosPion->SetDCAReCalculation(true);
  fTrackCutsPosPion->SetDCAVtxZ(0.27);
  fTrackCutsPosPion->SetDCAVtxXY(0.27);
  fTrackCutsPosPion->SetPID(AliPID::kPion, 0.49);
  fTrackCutsPosPion->SetRejLowPtPionsTOF(false);
  fTrackCutsPosPion->SetMinimalBooking(false);

  //The same things for negative pions
  fTrackCutsNegPion->SetIsMonteCarlo(isMC);
  fTrackCutsNegPion->SetCutCharge(-1);
  fTrackCutsNegPion->SetFilterBit(96);
  fTrackCutsNegPion->SetPtRange(0.14, 4.7);
  fTrackCutsNegPion->SetEtaRange(-0.86, 0.86);
  fTrackCutsNegPion->SetNClsTPC(77);
  fTrackCutsNegPion->SetDCAReCalculation(true);
  fTrackCutsNegPion->SetDCAVtxZ(0.27);
  fTrackCutsNegPion->SetDCAVtxXY(0.27);
  fTrackCutsNegPion->SetPID(AliPID::kPion, 0.49);
  fTrackCutsNegPion->SetRejLowPtPionsTOF(false);
  fTrackCutsNegPion->SetMinimalBooking(false);
  } else if ( suffix == "19" ) {//CPR low
  fSpherDown = 0.68;
  fdPhidEta=0.01;
  //Track Cuts are defined here
  //positive pions
  fTrackCutsPosPion->SetIsMonteCarlo(isMC);
  fTrackCutsPosPion->SetCutCharge(1);
  fTrackCutsPosPion->SetFilterBit(96);
  fTrackCutsPosPion->SetPtRange(0.11, 4.2);
  fTrackCutsPosPion->SetEtaRange(-0.75, 0.75);
  fTrackCutsPosPion->SetNClsTPC(77);
  fTrackCutsPosPion->SetDCAReCalculation(true);
  fTrackCutsPosPion->SetDCAVtxZ(0.26);
  fTrackCutsPosPion->SetDCAVtxXY(0.26);
  fTrackCutsPosPion->SetPID(AliPID::kPion, 0.47);
  fTrackCutsPosPion->SetRejLowPtPionsTOF(false);
  fTrackCutsPosPion->SetMinimalBooking(false);

  //The same things for negative pions
  fTrackCutsNegPion->SetIsMonteCarlo(isMC);
  fTrackCutsNegPion->SetCutCharge(-1);
  fTrackCutsNegPion->SetFilterBit(96);
  fTrackCutsNegPion->SetPtRange(0.11, 4.2);
  fTrackCutsNegPion->SetEtaRange(-0.75, 0.75);
  fTrackCutsNegPion->SetNClsTPC(77);
  fTrackCutsNegPion->SetDCAReCalculation(true);
  fTrackCutsNegPion->SetDCAVtxZ(0.26);
  fTrackCutsNegPion->SetDCAVtxXY(0.26);
  fTrackCutsNegPion->SetPID(AliPID::kPion, 0.47);
  fTrackCutsNegPion->SetRejLowPtPionsTOF(false);
  fTrackCutsNegPion->SetMinimalBooking(false);
  } else if ( suffix == "20" ) {//CPR high
  fSpherDown = 0.67;
  fdPhidEta=0.009;
  //Track Cuts are defined here
  //positive pions
  fTrackCutsPosPion->SetIsMonteCarlo(isMC);
  fTrackCutsPosPion->SetCutCharge(1);
  fTrackCutsPosPion->SetFilterBit(96);
  fTrackCutsPosPion->SetPtRange(0.12, 4.4);
  fTrackCutsPosPion->SetEtaRange(-0.83, 0.83);
  fTrackCutsPosPion->SetNClsTPC(84);
  fTrackCutsPosPion->SetDCAReCalculation(true);
  fTrackCutsPosPion->SetDCAVtxZ(0.35);
  fTrackCutsPosPion->SetDCAVtxXY(0.35);
  fTrackCutsPosPion->SetPID(AliPID::kPion, 0.54);
  fTrackCutsPosPion->SetRejLowPtPionsTOF(false);
  fTrackCutsPosPion->SetMinimalBooking(false);

  //The same things for negative pions
  fTrackCutsNegPion->SetIsMonteCarlo(isMC);
  fTrackCutsNegPion->SetCutCharge(-1);
  fTrackCutsNegPion->SetFilterBit(96);
  fTrackCutsNegPion->SetPtRange(0.12, 4.4);
  fTrackCutsNegPion->SetEtaRange(-0.83, 0.83);
  fTrackCutsNegPion->SetNClsTPC(84);
  fTrackCutsNegPion->SetDCAReCalculation(true);
  fTrackCutsNegPion->SetDCAVtxZ(0.35);
  fTrackCutsNegPion->SetDCAVtxXY(0.35);
  fTrackCutsNegPion->SetPID(AliPID::kPion, 0.54);
  fTrackCutsNegPion->SetRejLowPtPionsTOF(false);
  fTrackCutsNegPion->SetMinimalBooking(false);
  } else if (suffix == "21") {//pT high
  fSpherDown = 0.67;
  //Track Cuts are defined here
  //positive pions
  fdPhidEta=0.009;
  fTrackCutsPosPion->SetIsMonteCarlo(isMC);
  fTrackCutsPosPion->SetCutCharge(1);
  fTrackCutsPosPion->SetFilterBit(96);
  fTrackCutsPosPion->SetPtRange(0.14, 5.0);
  fTrackCutsPosPion->SetEtaRange(-0.84, 0.84);
  fTrackCutsPosPion->SetNClsTPC(87);
  fTrackCutsPosPion->SetDCAReCalculation(true);
  fTrackCutsPosPion->SetDCAVtxZ(0.26);
  fTrackCutsPosPion->SetDCAVtxXY(0.26);
  fTrackCutsPosPion->SetPID(AliPID::kPion, 0.51);
  fTrackCutsPosPion->SetRejLowPtPionsTOF(false);
  fTrackCutsPosPion->SetMinimalBooking(false);

  //The same things for negative pions
  fTrackCutsNegPion->SetIsMonteCarlo(isMC);
  fTrackCutsNegPion->SetCutCharge(-1);
  fTrackCutsNegPion->SetFilterBit(96);
  fTrackCutsNegPion->SetPtRange(0.14, 5.0);
  fTrackCutsNegPion->SetEtaRange(-0.84, 0.84);
  fTrackCutsNegPion->SetNClsTPC(87);
  fTrackCutsNegPion->SetDCAReCalculation(true);
  fTrackCutsNegPion->SetDCAVtxZ(0.26);
  fTrackCutsNegPion->SetDCAVtxXY(0.26);
  fTrackCutsNegPion->SetPID(AliPID::kPion, 0.51);
  fTrackCutsNegPion->SetRejLowPtPionsTOF(false);
  fTrackCutsNegPion->SetMinimalBooking(false);
  } else if ( suffix == "22" ) {//eta low
  fSpherDown = 0.68;
  //Track Cuts are defined here
  //positive pions
  fdPhidEta=0.008;
  fTrackCutsPosPion->SetIsMonteCarlo(isMC);
  fTrackCutsPosPion->SetCutCharge(1);
  fTrackCutsPosPion->SetFilterBit(96);
  fTrackCutsPosPion->SetPtRange(0.11, 4.3);
  fTrackCutsPosPion->SetEtaRange(-0.75, 0.75);
  fTrackCutsPosPion->SetNClsTPC(82);
  fTrackCutsPosPion->SetDCAReCalculation(true);
  fTrackCutsPosPion->SetDCAVtxZ(0.34);
  fTrackCutsPosPion->SetDCAVtxXY(0.34);
  fTrackCutsPosPion->SetPID(AliPID::kPion, 0.49);
  fTrackCutsPosPion->SetRejLowPtPionsTOF(false);
  fTrackCutsPosPion->SetMinimalBooking(false);

  //The same things for negative pions
  fTrackCutsNegPion->SetIsMonteCarlo(isMC);
  fTrackCutsNegPion->SetCutCharge(-1);
  fTrackCutsNegPion->SetFilterBit(96);
  fTrackCutsNegPion->SetPtRange(0.11, 4.3);
  fTrackCutsNegPion->SetEtaRange(-0.75, 0.75);
  fTrackCutsNegPion->SetNClsTPC(82);
  fTrackCutsNegPion->SetDCAReCalculation(true);
  fTrackCutsNegPion->SetDCAVtxZ(0.34);
  fTrackCutsNegPion->SetDCAVtxXY(0.34);
  fTrackCutsNegPion->SetPID(AliPID::kPion, 0.49);
  fTrackCutsNegPion->SetRejLowPtPionsTOF(false);
  fTrackCutsNegPion->SetMinimalBooking(false);
  } else if ( suffix == "23" ) {//eta high
  fSpherDown = 0.69;
  //Track Cuts are defined here
  //positive pions
  fdPhidEta=0.008;
  fTrackCutsPosPion->SetIsMonteCarlo(isMC);
  fTrackCutsPosPion->SetCutCharge(1);
  fTrackCutsPosPion->SetFilterBit(96);
  fTrackCutsPosPion->SetPtRange(0.12, 4.7);
  fTrackCutsPosPion->SetEtaRange(-0.9, 0.9);
  fTrackCutsPosPion->SetNClsTPC(79);
  fTrackCutsPosPion->SetDCAReCalculation(true);
  fTrackCutsPosPion->SetDCAVtxZ(0.26);
  fTrackCutsPosPion->SetDCAVtxXY(0.26);
  fTrackCutsPosPion->SetPID(AliPID::kPion, 0.47);
  fTrackCutsPosPion->SetRejLowPtPionsTOF(false);
  fTrackCutsPosPion->SetMinimalBooking(false);

  //The same things for negative pions
  fTrackCutsNegPion->SetIsMonteCarlo(isMC);
  fTrackCutsNegPion->SetCutCharge(-1);
  fTrackCutsNegPion->SetFilterBit(96);
  fTrackCutsNegPion->SetPtRange(0.12, 4.7);
  fTrackCutsNegPion->SetEtaRange(-0.9, 0.9);
  fTrackCutsNegPion->SetNClsTPC(79);
  fTrackCutsNegPion->SetDCAReCalculation(true);
  fTrackCutsNegPion->SetDCAVtxZ(0.26);
  fTrackCutsNegPion->SetDCAVtxXY(0.26);
  fTrackCutsNegPion->SetPID(AliPID::kPion, 0.47);
  fTrackCutsNegPion->SetRejLowPtPionsTOF(false);
  fTrackCutsNegPion->SetMinimalBooking(false);
  } else if ( suffix == "24" ) {//eta high conservative
  fSpherDown = 0.7;
  //Track Cuts are defined here
  //positive pions
  fdPhidEta=0.008;
  fTrackCutsPosPion->SetIsMonteCarlo(isMC);
  fTrackCutsPosPion->SetCutCharge(1);
  fTrackCutsPosPion->SetFilterBit(96);
  fTrackCutsPosPion->SetPtRange(0.13, 4.3);
  fTrackCutsPosPion->SetEtaRange(-0.85, 0.85);
  fTrackCutsPosPion->SetNClsTPC(76);
  fTrackCutsPosPion->SetDCAReCalculation(true);
  fTrackCutsPosPion->SetDCAVtxZ(0.3);
  fTrackCutsPosPion->SetDCAVtxXY(0.3);
  fTrackCutsPosPion->SetPID(AliPID::kPion, 0.5);
  fTrackCutsPosPion->SetRejLowPtPionsTOF(false);
  fTrackCutsPosPion->SetMinimalBooking(false);

  //The same things for negative pions
  fTrackCutsNegPion->SetIsMonteCarlo(isMC);
  fTrackCutsNegPion->SetCutCharge(-1);
  fTrackCutsNegPion->SetFilterBit(96);
  fTrackCutsNegPion->SetPtRange(0.13, 4.3);
  fTrackCutsNegPion->SetEtaRange(-0.85, 0.85);
  fTrackCutsNegPion->SetNClsTPC(76);
  fTrackCutsNegPion->SetDCAReCalculation(true);
  fTrackCutsNegPion->SetDCAVtxZ(0.3);
  fTrackCutsNegPion->SetDCAVtxXY(0.3);
  fTrackCutsNegPion->SetPID(AliPID::kPion, 0.5);
  fTrackCutsNegPion->SetRejLowPtPionsTOF(false);
  fTrackCutsNegPion->SetMinimalBooking(false);
  } else if ( suffix == "25" ) {//eta low conservative
  fSpherDown = 0.7;
  //Track Cuts are defined here
  //positive pions
  fdPhidEta=0.008;
  fTrackCutsPosPion->SetIsMonteCarlo(isMC);
  fTrackCutsPosPion->SetCutCharge(1);
  fTrackCutsPosPion->SetFilterBit(96);
  fTrackCutsPosPion->SetPtRange(0.14, 4.5);
  fTrackCutsPosPion->SetEtaRange(-0.75, 0.75);
  fTrackCutsPosPion->SetNClsTPC(87);
  fTrackCutsPosPion->SetDCAReCalculation(true);
  fTrackCutsPosPion->SetDCAVtxZ(0.32);
  fTrackCutsPosPion->SetDCAVtxXY(0.32);
  fTrackCutsPosPion->SetPID(AliPID::kPion, 0.48);
  fTrackCutsPosPion->SetRejLowPtPionsTOF(false);
  fTrackCutsPosPion->SetMinimalBooking(false);

  //The same things for negative pions
  fTrackCutsNegPion->SetIsMonteCarlo(isMC);
  fTrackCutsNegPion->SetCutCharge(-1);
  fTrackCutsNegPion->SetFilterBit(96);
  fTrackCutsNegPion->SetPtRange(0.14, 4.5);
  fTrackCutsNegPion->SetEtaRange(-0.75, 0.75);
  fTrackCutsNegPion->SetNClsTPC(87);
  fTrackCutsNegPion->SetDCAReCalculation(true);
  fTrackCutsNegPion->SetDCAVtxZ(0.32);
  fTrackCutsNegPion->SetDCAVtxXY(0.32);
  fTrackCutsNegPion->SetPID(AliPID::kPion, 0.48);
  fTrackCutsNegPion->SetRejLowPtPionsTOF(false);
  fTrackCutsNegPion->SetMinimalBooking(false);
  } else if ( suffix == "26" ) {//#TPC Cls low
  fSpherDown = 0.71;
  //Track Cuts are defined here
  //positive pions
  fdPhidEta=0.009;
  fTrackCutsPosPion->SetIsMonteCarlo(isMC);
  fTrackCutsPosPion->SetCutCharge(1);
  fTrackCutsPosPion->SetFilterBit(96);
  fTrackCutsPosPion->SetPtRange(0.14, 4.0);
  fTrackCutsPosPion->SetEtaRange(-0.77, 0.77);
  fTrackCutsPosPion->SetNClsTPC(90);
  fTrackCutsPosPion->SetDCAReCalculation(true);
  fTrackCutsPosPion->SetDCAVtxZ(0.3);
  fTrackCutsPosPion->SetDCAVtxXY(0.3);
  fTrackCutsPosPion->SetPID(AliPID::kPion, 0.46);
  fTrackCutsPosPion->SetRejLowPtPionsTOF(false);
  fTrackCutsPosPion->SetMinimalBooking(false);

  //The same things for negative pions
  fTrackCutsNegPion->SetIsMonteCarlo(isMC);
  fTrackCutsNegPion->SetCutCharge(-1);
  fTrackCutsNegPion->SetFilterBit(96);
  fTrackCutsNegPion->SetPtRange(0.14, 4.0);
  fTrackCutsNegPion->SetEtaRange(-0.77, 0.77);
  fTrackCutsNegPion->SetNClsTPC(90);
  fTrackCutsNegPion->SetDCAReCalculation(true);
  fTrackCutsNegPion->SetDCAVtxZ(0.3);
  fTrackCutsNegPion->SetDCAVtxXY(0.3);
  fTrackCutsNegPion->SetPID(AliPID::kPion, 0.46);
  fTrackCutsNegPion->SetRejLowPtPionsTOF(false);
  fTrackCutsNegPion->SetMinimalBooking(false);
  } else if ( suffix == "27" ) {//#TPC Cls high conservative
  fSpherDown = 0.71;
  //Track Cuts are defined here
  //positive pions
  fdPhidEta=0.009;
  fTrackCutsPosPion->SetIsMonteCarlo(isMC);
  fTrackCutsPosPion->SetCutCharge(1);
  fTrackCutsPosPion->SetFilterBit(96);
  fTrackCutsPosPion->SetPtRange(0.13, 4.3);
  fTrackCutsPosPion->SetEtaRange(-0.82, 0.82);
  fTrackCutsPosPion->SetNClsTPC(85);
  fTrackCutsPosPion->SetDCAReCalculation(true);
  fTrackCutsPosPion->SetDCAVtxZ(0.29);
  fTrackCutsPosPion->SetDCAVtxXY(0.29);
  fTrackCutsPosPion->SetPID(AliPID::kPion, 0.54);
  fTrackCutsPosPion->SetRejLowPtPionsTOF(false);
  fTrackCutsPosPion->SetMinimalBooking(false);

  //The same things for negative pions
  fTrackCutsNegPion->SetIsMonteCarlo(isMC);
  fTrackCutsNegPion->SetCutCharge(-1);
  fTrackCutsNegPion->SetFilterBit(96);
  fTrackCutsNegPion->SetPtRange(0.13, 4.3);
  fTrackCutsNegPion->SetEtaRange(-0.82, 0.82);
  fTrackCutsNegPion->SetNClsTPC(85);
  fTrackCutsNegPion->SetDCAReCalculation(true);
  fTrackCutsNegPion->SetDCAVtxZ(0.29);
  fTrackCutsNegPion->SetDCAVtxXY(0.29);
  fTrackCutsNegPion->SetPID(AliPID::kPion, 0.54);
  fTrackCutsNegPion->SetRejLowPtPionsTOF(false);
  fTrackCutsNegPion->SetMinimalBooking(false);
  } else if ( suffix == "28" ) {//#TPC Cls high
  fSpherDown = 0.72;
  //Track Cuts are defined here
  //positive pions
  fdPhidEta=0.01;
  fTrackCutsPosPion->SetIsMonteCarlo(isMC);
  fTrackCutsPosPion->SetCutCharge(1);
  fTrackCutsPosPion->SetFilterBit(96);
  fTrackCutsPosPion->SetPtRange(0.14, 4.8);
  fTrackCutsPosPion->SetEtaRange(-0.89, 0.89);
  fTrackCutsPosPion->SetNClsTPC(90);
  fTrackCutsPosPion->SetDCAReCalculation(true);
  fTrackCutsPosPion->SetDCAVtxZ(0.26);
  fTrackCutsPosPion->SetDCAVtxXY(0.26);
  fTrackCutsPosPion->SetPID(AliPID::kPion, 0.47);
  fTrackCutsPosPion->SetRejLowPtPionsTOF(false);
  fTrackCutsPosPion->SetMinimalBooking(false);

  //The same things for negative pions
  fTrackCutsNegPion->SetIsMonteCarlo(isMC);
  fTrackCutsNegPion->SetCutCharge(-1);
  fTrackCutsNegPion->SetFilterBit(96);
  fTrackCutsNegPion->SetPtRange(0.14, 4.8);
  fTrackCutsNegPion->SetEtaRange(-0.89, 0.89);
  fTrackCutsNegPion->SetNClsTPC(90);
  fTrackCutsNegPion->SetDCAReCalculation(true);
  fTrackCutsNegPion->SetDCAVtxZ(0.26);
  fTrackCutsNegPion->SetDCAVtxXY(0.26);
  fTrackCutsNegPion->SetPID(AliPID::kPion, 0.47);
  fTrackCutsNegPion->SetRejLowPtPionsTOF(false);
  fTrackCutsNegPion->SetMinimalBooking(false);
  } else if ( suffix == "29" ) {//#TPC Cls low conservative
  fSpherDown = 0.72;
  //Track Cuts are defined here
  //positive pions
  fdPhidEta=0.01;
  fTrackCutsPosPion->SetIsMonteCarlo(isMC);
  fTrackCutsPosPion->SetCutCharge(1);
  fTrackCutsPosPion->SetFilterBit(96);
  fTrackCutsPosPion->SetPtRange(0.13, 4.3);
  fTrackCutsPosPion->SetEtaRange(-0.76, 0.76);
  fTrackCutsPosPion->SetNClsTPC(75);
  fTrackCutsPosPion->SetDCAReCalculation(true);
  fTrackCutsPosPion->SetDCAVtxZ(0.33);
  fTrackCutsPosPion->SetDCAVtxXY(0.33);
  fTrackCutsPosPion->SetPID(AliPID::kPion, 0.5);
  fTrackCutsPosPion->SetRejLowPtPionsTOF(false);
  fTrackCutsPosPion->SetMinimalBooking(false);

  //The same things for negative pions
  fTrackCutsNegPion->SetIsMonteCarlo(isMC);
  fTrackCutsNegPion->SetCutCharge(-1);
  fTrackCutsNegPion->SetFilterBit(96);
  fTrackCutsNegPion->SetPtRange(0.13, 4.3);
  fTrackCutsNegPion->SetEtaRange(-0.76, 0.76);
  fTrackCutsNegPion->SetNClsTPC(75);
  fTrackCutsNegPion->SetDCAReCalculation(true);
  fTrackCutsNegPion->SetDCAVtxZ(0.33);
  fTrackCutsNegPion->SetDCAVtxXY(0.33);
  fTrackCutsNegPion->SetPID(AliPID::kPion, 0.5);
  fTrackCutsNegPion->SetRejLowPtPionsTOF(false);
  fTrackCutsNegPion->SetMinimalBooking(false);
  } else if ( suffix == "30" ) {//DCA XYZ low conservative
  fSpherDown = 0.73;
  //Track Cuts are defined here
  //positive pions
  fdPhidEta=0.011;
  fTrackCutsPosPion->SetIsMonteCarlo(isMC);
  fTrackCutsPosPion->SetCutCharge(1);
  fTrackCutsPosPion->SetFilterBit(96);
  fTrackCutsPosPion->SetPtRange(0.14, 4.8);
  fTrackCutsPosPion->SetEtaRange(-0.83, 0.83);
  fTrackCutsPosPion->SetNClsTPC(87);
  fTrackCutsPosPion->SetDCAReCalculation(true);
  fTrackCutsPosPion->SetDCAVtxZ(0.28);
  fTrackCutsPosPion->SetDCAVtxXY(0.28);
  fTrackCutsPosPion->SetPID(AliPID::kPion, 0.46);
  fTrackCutsPosPion->SetRejLowPtPionsTOF(false);
  fTrackCutsPosPion->SetMinimalBooking(false);

  //The same things for negative pions
  fTrackCutsNegPion->SetIsMonteCarlo(isMC);
  fTrackCutsNegPion->SetCutCharge(-1);
  fTrackCutsNegPion->SetFilterBit(96);
  fTrackCutsNegPion->SetPtRange(0.14, 4.8);
  fTrackCutsNegPion->SetEtaRange(-0.83, 0.83);
  fTrackCutsNegPion->SetNClsTPC(87);
  fTrackCutsNegPion->SetDCAReCalculation(true);
  fTrackCutsNegPion->SetDCAVtxZ(0.28);
  fTrackCutsNegPion->SetDCAVtxXY(0.28);
  fTrackCutsNegPion->SetPID(AliPID::kPion, 0.46);
  fTrackCutsNegPion->SetRejLowPtPionsTOF(false);
  fTrackCutsNegPion->SetMinimalBooking(false);
  } else if ( suffix == "31" ) {//DCA XYZ high conservative
  fSpherDown = 0.73;
  //Track Cuts are defined here
  //positive pions
  fdPhidEta=0.011;
  fTrackCutsPosPion->SetIsMonteCarlo(isMC);
  fTrackCutsPosPion->SetCutCharge(1);
  fTrackCutsPosPion->SetFilterBit(96);
  fTrackCutsPosPion->SetPtRange(0.11, 4.1);
  fTrackCutsPosPion->SetEtaRange(-0.77, 0.77);
  fTrackCutsPosPion->SetNClsTPC(76);
  fTrackCutsPosPion->SetDCAReCalculation(true);
  fTrackCutsPosPion->SetDCAVtxZ(0.32);
  fTrackCutsPosPion->SetDCAVtxXY(0.32);
  fTrackCutsPosPion->SetPID(AliPID::kPion, 0.54);
  fTrackCutsPosPion->SetRejLowPtPionsTOF(false);
  fTrackCutsPosPion->SetMinimalBooking(false);

  //The same things for negative pions
  fTrackCutsNegPion->SetIsMonteCarlo(isMC);
  fTrackCutsNegPion->SetCutCharge(-1);
  fTrackCutsNegPion->SetFilterBit(96);
  fTrackCutsNegPion->SetPtRange(0.11, 4.1);
  fTrackCutsNegPion->SetEtaRange(-0.77, 0.77);
  fTrackCutsNegPion->SetNClsTPC(76);
  fTrackCutsNegPion->SetDCAReCalculation(true);
  fTrackCutsNegPion->SetDCAVtxZ(0.32);
  fTrackCutsNegPion->SetDCAVtxXY(0.32);
  fTrackCutsNegPion->SetPID(AliPID::kPion, 0.54);
  fTrackCutsNegPion->SetRejLowPtPionsTOF(false);
  fTrackCutsNegPion->SetMinimalBooking(false);
  } else if ( suffix == "32" ) {//DCA XYZ high
  fSpherDown = 0.69;
  //Track Cuts are defined here
  //positive pions
  fdPhidEta=0.012;
  fTrackCutsPosPion->SetIsMonteCarlo(isMC);
  fTrackCutsPosPion->SetCutCharge(1);
  fTrackCutsPosPion->SetFilterBit(96);
  fTrackCutsPosPion->SetPtRange(0.11, 4.7);
  fTrackCutsPosPion->SetEtaRange(-0.84, 0.84);
  fTrackCutsPosPion->SetNClsTPC(77);
  fTrackCutsPosPion->SetDCAReCalculation(true);
  fTrackCutsPosPion->SetDCAVtxZ(0.35);
  fTrackCutsPosPion->SetDCAVtxXY(0.35);
  fTrackCutsPosPion->SetPID(AliPID::kPion, 0.46);
  fTrackCutsPosPion->SetRejLowPtPionsTOF(false);
  fTrackCutsPosPion->SetMinimalBooking(false);

  //The same things for negative pions
  fTrackCutsNegPion->SetIsMonteCarlo(isMC);
  fTrackCutsNegPion->SetCutCharge(-1);
  fTrackCutsNegPion->SetFilterBit(96);
  fTrackCutsNegPion->SetPtRange(0.11, 4.7);
  fTrackCutsNegPion->SetEtaRange(-0.84, 0.84);
  fTrackCutsNegPion->SetNClsTPC(77);
  fTrackCutsNegPion->SetDCAReCalculation(true);
  fTrackCutsNegPion->SetDCAVtxZ(0.35);
  fTrackCutsNegPion->SetDCAVtxXY(0.35);
  fTrackCutsNegPion->SetPID(AliPID::kPion, 0.46);
  fTrackCutsNegPion->SetRejLowPtPionsTOF(false);
  fTrackCutsNegPion->SetMinimalBooking(false);
  } else if ( suffix == "33" ) {//DCA XYZ low
  fSpherDown = 0.69;
  //Track Cuts are defined here
  //positive pions
  fdPhidEta=0.012;
  fTrackCutsPosPion->SetIsMonteCarlo(isMC);
  fTrackCutsPosPion->SetCutCharge(1);
  fTrackCutsPosPion->SetFilterBit(96);
  fTrackCutsPosPion->SetPtRange(0.12, 4.9);
  fTrackCutsPosPion->SetEtaRange(-0.89, 0.89);
  fTrackCutsPosPion->SetNClsTPC(75);
  fTrackCutsPosPion->SetDCAReCalculation(true);
  fTrackCutsPosPion->SetDCAVtxZ(0.25);
  fTrackCutsPosPion->SetDCAVtxXY(0.25);
  fTrackCutsPosPion->SetPID(AliPID::kPion, 0.47);
  fTrackCutsPosPion->SetRejLowPtPionsTOF(false);
  fTrackCutsPosPion->SetMinimalBooking(false);

  //The same things for negative pions
  fTrackCutsNegPion->SetIsMonteCarlo(isMC);
  fTrackCutsNegPion->SetCutCharge(-1);
  fTrackCutsNegPion->SetFilterBit(96);
  fTrackCutsNegPion->SetPtRange(0.12, 4.9);
  fTrackCutsNegPion->SetEtaRange(-0.89, 0.89);
  fTrackCutsNegPion->SetNClsTPC(75);
  fTrackCutsNegPion->SetDCAReCalculation(true);
  fTrackCutsNegPion->SetDCAVtxZ(0.25);
  fTrackCutsNegPion->SetDCAVtxXY(0.25);
  fTrackCutsNegPion->SetPID(AliPID::kPion, 0.47);
  fTrackCutsNegPion->SetRejLowPtPionsTOF(false);
  fTrackCutsNegPion->SetMinimalBooking(false);
  } else if ( suffix == "34" ) {//Sphericity low
  fSpherDown = 0.68;
  //Track Cuts are defined here
  //positive pions
  fdPhidEta=0.008;
  fTrackCutsPosPion->SetIsMonteCarlo(isMC);
  fTrackCutsPosPion->SetCutCharge(1);
  fTrackCutsPosPion->SetFilterBit(96);
  fTrackCutsPosPion->SetPtRange(0.14, 4.1);
  fTrackCutsPosPion->SetEtaRange(-0.77, 0.77);
  fTrackCutsPosPion->SetNClsTPC(78);
  fTrackCutsPosPion->SetDCAReCalculation(true);
  fTrackCutsPosPion->SetDCAVtxZ(0.29);
  fTrackCutsPosPion->SetDCAVtxXY(0.29);
  fTrackCutsPosPion->SetPID(AliPID::kPion, 0.53);
  fTrackCutsPosPion->SetRejLowPtPionsTOF(false);
  fTrackCutsPosPion->SetMinimalBooking(false);

  //The same things for negative pions
  fTrackCutsNegPion->SetIsMonteCarlo(isMC);
  fTrackCutsNegPion->SetCutCharge(-1);
  fTrackCutsNegPion->SetFilterBit(96);
  fTrackCutsNegPion->SetPtRange(0.14, 4.1);
  fTrackCutsNegPion->SetEtaRange(-0.77, 0.77);
  fTrackCutsNegPion->SetNClsTPC(78);
  fTrackCutsNegPion->SetDCAReCalculation(true);
  fTrackCutsNegPion->SetDCAVtxZ(0.29);
  fTrackCutsNegPion->SetDCAVtxXY(0.29);
  fTrackCutsNegPion->SetPID(AliPID::kPion, 0.53);
  fTrackCutsNegPion->SetRejLowPtPionsTOF(false);
  fTrackCutsNegPion->SetMinimalBooking(false);
  } else if ( suffix == "35" ) {//Sphericity low
  fSpherDown = 0.68;
  //Track Cuts are defined here
  //positive pions
  fdPhidEta=0.009;
  fTrackCutsPosPion->SetIsMonteCarlo(isMC);
  fTrackCutsPosPion->SetCutCharge(1);
  fTrackCutsPosPion->SetFilterBit(96);
  fTrackCutsPosPion->SetPtRange(0.11, 4.8);
  fTrackCutsPosPion->SetEtaRange(-0.76, 0.76);
  fTrackCutsPosPion->SetNClsTPC(79);
  fTrackCutsPosPion->SetDCAReCalculation(true);
  fTrackCutsPosPion->SetDCAVtxZ(0.34);
  fTrackCutsPosPion->SetDCAVtxXY(0.34);
  fTrackCutsPosPion->SetPID(AliPID::kPion, 0.55);
  fTrackCutsPosPion->SetRejLowPtPionsTOF(false);
  fTrackCutsPosPion->SetMinimalBooking(false);

  //The same things for negative pions
  fTrackCutsNegPion->SetIsMonteCarlo(isMC);
  fTrackCutsNegPion->SetCutCharge(-1);
  fTrackCutsNegPion->SetFilterBit(96);
  fTrackCutsNegPion->SetPtRange(0.11, 4.8);
  fTrackCutsNegPion->SetEtaRange(-0.76, 0.76);
  fTrackCutsNegPion->SetNClsTPC(79);
  fTrackCutsNegPion->SetDCAReCalculation(true);
  fTrackCutsNegPion->SetDCAVtxZ(0.34);
  fTrackCutsNegPion->SetDCAVtxXY(0.34);
  fTrackCutsNegPion->SetPID(AliPID::kPion, 0.55);
  fTrackCutsNegPion->SetRejLowPtPionsTOF(false);
  fTrackCutsNegPion->SetMinimalBooking(false);
  } else if ( suffix == "36" ) {//Sphericity low
  fSpherDown = 0.72;
  //Track Cuts are defined here
  //positive pions
  fdPhidEta=0.01;
  fTrackCutsPosPion->SetIsMonteCarlo(isMC);
  fTrackCutsPosPion->SetCutCharge(1);
  fTrackCutsPosPion->SetFilterBit(96);
  fTrackCutsPosPion->SetPtRange(0.14, 4.0);
  fTrackCutsPosPion->SetEtaRange(-0.83, 0.83);
  fTrackCutsPosPion->SetNClsTPC(77);
  fTrackCutsPosPion->SetDCAReCalculation(true);
  fTrackCutsPosPion->SetDCAVtxZ(0.3);
  fTrackCutsPosPion->SetDCAVtxXY(0.3);
  fTrackCutsPosPion->SetPID(AliPID::kPion, 0.45);
  fTrackCutsPosPion->SetRejLowPtPionsTOF(false);
  fTrackCutsPosPion->SetMinimalBooking(false);

  //The same things for negative pions
  fTrackCutsNegPion->SetIsMonteCarlo(isMC);
  fTrackCutsNegPion->SetCutCharge(-1);
  fTrackCutsNegPion->SetFilterBit(96);
  fTrackCutsNegPion->SetPtRange(0.14, 4.0);
  fTrackCutsNegPion->SetEtaRange(-0.83, 0.83);
  fTrackCutsNegPion->SetNClsTPC(77);
  fTrackCutsNegPion->SetDCAReCalculation(true);
  fTrackCutsNegPion->SetDCAVtxZ(0.3);
  fTrackCutsNegPion->SetDCAVtxXY(0.3);
  fTrackCutsNegPion->SetPID(AliPID::kPion, 0.45);
  fTrackCutsNegPion->SetRejLowPtPionsTOF(false);
  fTrackCutsNegPion->SetMinimalBooking(false);
  } else if ( suffix == "37" ) {//Sphericity high
  fSpherDown = 0.71;
  //Track Cuts are defined here
  //positive pions
  fdPhidEta=0.011;
  fTrackCutsPosPion->SetIsMonteCarlo(isMC);
  fTrackCutsPosPion->SetCutCharge(1);
  fTrackCutsPosPion->SetFilterBit(96);
  fTrackCutsPosPion->SetPtRange(0.11, 4.0);
  fTrackCutsPosPion->SetEtaRange(-0.86, 0.86);
  fTrackCutsPosPion->SetNClsTPC(89);
  fTrackCutsPosPion->SetDCAReCalculation(true);
  fTrackCutsPosPion->SetDCAVtxZ(0.27);
  fTrackCutsPosPion->SetDCAVtxXY(0.27);
  fTrackCutsPosPion->SetPID(AliPID::kPion, 0.53);
  fTrackCutsPosPion->SetRejLowPtPionsTOF(false);
  fTrackCutsPosPion->SetMinimalBooking(false);

  //The same things for negative pions
  fTrackCutsNegPion->SetIsMonteCarlo(isMC);
  fTrackCutsNegPion->SetCutCharge(-1);
  fTrackCutsNegPion->SetFilterBit(96);
  fTrackCutsNegPion->SetPtRange(0.11, 4.0);
  fTrackCutsNegPion->SetEtaRange(-0.86, 0.86);
  fTrackCutsNegPion->SetNClsTPC(89);
  fTrackCutsNegPion->SetDCAReCalculation(true);
  fTrackCutsNegPion->SetDCAVtxZ(0.27);
  fTrackCutsNegPion->SetDCAVtxXY(0.27);
  fTrackCutsNegPion->SetPID(AliPID::kPion, 0.53);
  fTrackCutsNegPion->SetRejLowPtPionsTOF(false);
  fTrackCutsNegPion->SetMinimalBooking(false);
  } else if ( suffix == "38" ) {//PID low
  fSpherDown = 0.67;
  //Track Cuts are defined here
  //positive pions
  fdPhidEta=0.012;
  fTrackCutsPosPion->SetIsMonteCarlo(isMC);
  fTrackCutsPosPion->SetCutCharge(1);
  fTrackCutsPosPion->SetFilterBit(96);
  fTrackCutsPosPion->SetPtRange(0.14, 5.0);
  fTrackCutsPosPion->SetEtaRange(-0.81, 0.81);
  fTrackCutsPosPion->SetNClsTPC(84);
  fTrackCutsPosPion->SetDCAReCalculation(true);
  fTrackCutsPosPion->SetDCAVtxZ(0.28);
  fTrackCutsPosPion->SetDCAVtxXY(0.28);
  fTrackCutsPosPion->SetPID(AliPID::kPion, 0.55);
  fTrackCutsPosPion->SetRejLowPtPionsTOF(false);
  fTrackCutsPosPion->SetMinimalBooking(false);

  //The same things for negative pions
  fTrackCutsNegPion->SetIsMonteCarlo(isMC);
  fTrackCutsNegPion->SetCutCharge(-1);
  fTrackCutsNegPion->SetFilterBit(96);
  fTrackCutsNegPion->SetPtRange(0.14, 5.0);
  fTrackCutsNegPion->SetEtaRange(-0.81, 0.81);
  fTrackCutsNegPion->SetNClsTPC(84);
  fTrackCutsNegPion->SetDCAReCalculation(true);
  fTrackCutsNegPion->SetDCAVtxZ(0.28);
  fTrackCutsNegPion->SetDCAVtxXY(0.28);
  fTrackCutsNegPion->SetPID(AliPID::kPion, 0.55);
  fTrackCutsNegPion->SetRejLowPtPionsTOF(false);
  fTrackCutsNegPion->SetMinimalBooking(false);
  } else if ( suffix == "39" ) {//PID high
  fSpherDown = 0.67;
  //Track Cuts are defined here
  //positive pions
  fdPhidEta=0.012;
  fTrackCutsPosPion->SetIsMonteCarlo(isMC);
  fTrackCutsPosPion->SetCutCharge(1);
  fTrackCutsPosPion->SetFilterBit(96);
  fTrackCutsPosPion->SetPtRange(0.13, 4.0);
  fTrackCutsPosPion->SetEtaRange(-0.8, 0.8);
  fTrackCutsPosPion->SetNClsTPC(75);
  fTrackCutsPosPion->SetDCAReCalculation(true);
  fTrackCutsPosPion->SetDCAVtxZ(0.32);
  fTrackCutsPosPion->SetDCAVtxXY(0.32);
  fTrackCutsPosPion->SetPID(AliPID::kPion, 0.45);
  fTrackCutsPosPion->SetRejLowPtPionsTOF(false);
  fTrackCutsPosPion->SetMinimalBooking(false);

  //The same things for negative pions
  fTrackCutsNegPion->SetIsMonteCarlo(isMC);
  fTrackCutsNegPion->SetCutCharge(-1);
  fTrackCutsNegPion->SetFilterBit(96);
  fTrackCutsNegPion->SetPtRange(0.13, 4.0);
  fTrackCutsNegPion->SetEtaRange(-0.8, 0.8);
  fTrackCutsNegPion->SetNClsTPC(75);
  fTrackCutsNegPion->SetDCAReCalculation(true);
  fTrackCutsNegPion->SetDCAVtxZ(0.32);
  fTrackCutsNegPion->SetDCAVtxXY(0.32);
  fTrackCutsNegPion->SetPID(AliPID::kPion, 0.45);
  fTrackCutsNegPion->SetRejLowPtPionsTOF(false);
  fTrackCutsNegPion->SetMinimalBooking(false);
  } else if ( suffix == "40" ) {//CPR low
  fSpherDown = 0.69;
  fdPhidEta=0.011;
  //Track Cuts are defined here
  //positive pions
  fTrackCutsPosPion->SetIsMonteCarlo(isMC);
  fTrackCutsPosPion->SetCutCharge(1);
  fTrackCutsPosPion->SetFilterBit(96);
  fTrackCutsPosPion->SetPtRange(0.14, 4.8);
  fTrackCutsPosPion->SetEtaRange(-0.78, 0.78);
  fTrackCutsPosPion->SetNClsTPC(79);
  fTrackCutsPosPion->SetDCAReCalculation(true);
  fTrackCutsPosPion->SetDCAVtxZ(0.27);
  fTrackCutsPosPion->SetDCAVtxXY(0.27);
  fTrackCutsPosPion->SetPID(AliPID::kPion, 0.48);
  fTrackCutsPosPion->SetRejLowPtPionsTOF(false);
  fTrackCutsPosPion->SetMinimalBooking(false);

  //The same things for negative pions
  fTrackCutsNegPion->SetIsMonteCarlo(isMC);
  fTrackCutsNegPion->SetCutCharge(-1);
  fTrackCutsNegPion->SetFilterBit(96);
  fTrackCutsNegPion->SetPtRange(0.14, 4.8);
  fTrackCutsNegPion->SetEtaRange(-0.78, 0.78);
  fTrackCutsNegPion->SetNClsTPC(79);
  fTrackCutsNegPion->SetDCAReCalculation(true);
  fTrackCutsNegPion->SetDCAVtxZ(0.27);
  fTrackCutsNegPion->SetDCAVtxXY(0.27);
  fTrackCutsNegPion->SetPID(AliPID::kPion, 0.48);
  fTrackCutsNegPion->SetRejLowPtPionsTOF(false);
  fTrackCutsNegPion->SetMinimalBooking(false);
  } else if ( suffix == "41" ) {//CPR high
  fSpherDown = 0.71;
  fdPhidEta=0.01;
  //Track Cuts are defined here
  //positive pions
  fTrackCutsPosPion->SetIsMonteCarlo(isMC);
  fTrackCutsPosPion->SetCutCharge(1);
  fTrackCutsPosPion->SetFilterBit(96);
  fTrackCutsPosPion->SetPtRange(0.14, 4.0);
  fTrackCutsPosPion->SetEtaRange(-0.79, 0.79);
  fTrackCutsPosPion->SetNClsTPC(84);
  fTrackCutsPosPion->SetDCAReCalculation(true);
  fTrackCutsPosPion->SetDCAVtxZ(0.26);
  fTrackCutsPosPion->SetDCAVtxXY(0.26);
  fTrackCutsPosPion->SetPID(AliPID::kPion, 0.5);
  fTrackCutsPosPion->SetRejLowPtPionsTOF(false);
  fTrackCutsPosPion->SetMinimalBooking(false);

  //The same things for negative pions
  fTrackCutsNegPion->SetIsMonteCarlo(isMC);
  fTrackCutsNegPion->SetCutCharge(-1);
  fTrackCutsNegPion->SetFilterBit(96);
  fTrackCutsNegPion->SetPtRange(0.14, 4.0);
  fTrackCutsNegPion->SetEtaRange(-0.79, 0.79);
  fTrackCutsNegPion->SetNClsTPC(84);
  fTrackCutsNegPion->SetDCAReCalculation(true);
  fTrackCutsNegPion->SetDCAVtxZ(0.26);
  fTrackCutsNegPion->SetDCAVtxXY(0.26);
  fTrackCutsNegPion->SetPID(AliPID::kPion, 0.5);
  fTrackCutsNegPion->SetRejLowPtPionsTOF(false);
  fTrackCutsNegPion->SetMinimalBooking(false);
  } else if ( suffix == "42" ) {//CPR low
  fSpherDown = 0.7;
  fdPhidEta=0.008;
  //Track Cuts are defined here
  //positive pions
  fTrackCutsPosPion->SetIsMonteCarlo(isMC);
  fTrackCutsPosPion->SetCutCharge(1);
  fTrackCutsPosPion->SetFilterBit(96);
  fTrackCutsPosPion->SetPtRange(0.13, 4.6);
  fTrackCutsPosPion->SetEtaRange(-0.77, 0.77);
  fTrackCutsPosPion->SetNClsTPC(80);
  fTrackCutsPosPion->SetDCAReCalculation(true);
  fTrackCutsPosPion->SetDCAVtxZ(0.31);
  fTrackCutsPosPion->SetDCAVtxXY(0.31);
  fTrackCutsPosPion->SetPID(AliPID::kPion, 0.5);
  fTrackCutsPosPion->SetRejLowPtPionsTOF(false);
  fTrackCutsPosPion->SetMinimalBooking(false);

  //The same things for negative pions
  fTrackCutsNegPion->SetIsMonteCarlo(isMC);
  fTrackCutsNegPion->SetCutCharge(-1);
  fTrackCutsNegPion->SetFilterBit(96);
  fTrackCutsNegPion->SetPtRange(0.13, 4.6);
  fTrackCutsNegPion->SetEtaRange(-0.77, 0.77);
  fTrackCutsNegPion->SetNClsTPC(80);
  fTrackCutsNegPion->SetDCAReCalculation(true);
  fTrackCutsNegPion->SetDCAVtxZ(0.31);
  fTrackCutsNegPion->SetDCAVtxXY(0.31);
  fTrackCutsNegPion->SetPID(AliPID::kPion, 0.5);
  fTrackCutsNegPion->SetRejLowPtPionsTOF(false);
  fTrackCutsNegPion->SetMinimalBooking(false);
  } else if ( suffix == "43" ) {//CPR high
  fSpherDown = 0.7;
  fdPhidEta=0.009;
  //Track Cuts are defined here
  //positive pions
  fTrackCutsPosPion->SetIsMonteCarlo(isMC);
  fTrackCutsPosPion->SetCutCharge(1);
  fTrackCutsPosPion->SetFilterBit(96);
  fTrackCutsPosPion->SetPtRange(0.12, 4.9);
  fTrackCutsPosPion->SetEtaRange(-0.88, 0.88);
  fTrackCutsPosPion->SetNClsTPC(77);
  fTrackCutsPosPion->SetDCAReCalculation(true);
  fTrackCutsPosPion->SetDCAVtxZ(0.34);
  fTrackCutsPosPion->SetDCAVtxXY(0.34);
  fTrackCutsPosPion->SetPID(AliPID::kPion, 0.46);
  fTrackCutsPosPion->SetRejLowPtPionsTOF(false);
  fTrackCutsPosPion->SetMinimalBooking(false);

  //The same things for negative pions
  fTrackCutsNegPion->SetIsMonteCarlo(isMC);
  fTrackCutsNegPion->SetCutCharge(-1);
  fTrackCutsNegPion->SetFilterBit(96);
  fTrackCutsNegPion->SetPtRange(0.12, 4.9);
  fTrackCutsNegPion->SetEtaRange(-0.88, 0.88);
  fTrackCutsNegPion->SetNClsTPC(77);
  fTrackCutsNegPion->SetDCAReCalculation(true);
  fTrackCutsNegPion->SetDCAVtxZ(0.34);
  fTrackCutsNegPion->SetDCAVtxXY(0.34);
  fTrackCutsNegPion->SetPID(AliPID::kPion, 0.46);
  fTrackCutsNegPion->SetRejLowPtPionsTOF(false);
  fTrackCutsNegPion->SetMinimalBooking(false);
  } else if ( suffix == "44" ) {//CPR high
  fSpherDown = 0.7;
  fdPhidEta=0.01;
  //Track Cuts are defined here
  //positive pions
  fTrackCutsPosPion->SetIsMonteCarlo(isMC);
  fTrackCutsPosPion->SetCutCharge(1);
  fTrackCutsPosPion->SetFilterBit(96);
  fTrackCutsPosPion->SetPtRange(0.14, 4.2);
  fTrackCutsPosPion->SetEtaRange(-0.84, 0.84);
  fTrackCutsPosPion->SetNClsTPC(78);
  fTrackCutsPosPion->SetDCAReCalculation(true);
  fTrackCutsPosPion->SetDCAVtxZ(0.32);
  fTrackCutsPosPion->SetDCAVtxXY(0.32);
  fTrackCutsPosPion->SetPID(AliPID::kPion, 0.51);
  fTrackCutsPosPion->SetRejLowPtPionsTOF(false);
  fTrackCutsPosPion->SetMinimalBooking(false);

  //The same things for negative pions
  fTrackCutsNegPion->SetIsMonteCarlo(isMC);
  fTrackCutsNegPion->SetCutCharge(-1);
  fTrackCutsNegPion->SetFilterBit(96);
  fTrackCutsNegPion->SetPtRange(0.14, 4.2);
  fTrackCutsNegPion->SetEtaRange(-0.84, 0.84);
  fTrackCutsNegPion->SetNClsTPC(78);
  fTrackCutsNegPion->SetDCAReCalculation(true);
  fTrackCutsNegPion->SetDCAVtxZ(0.32);
  fTrackCutsNegPion->SetDCAVtxXY(0.32);
  fTrackCutsNegPion->SetPID(AliPID::kPion, 0.51);
  fTrackCutsNegPion->SetRejLowPtPionsTOF(false);
  fTrackCutsNegPion->SetMinimalBooking(false);
  } else {//Default
  //Track Cuts are defined here
  //positive pions
  fTrackCutsPosPion->SetIsMonteCarlo(isMC);
  fTrackCutsPosPion->SetCutCharge(1);
  fTrackCutsPosPion->SetFilterBit(96);
  fTrackCutsPosPion->SetPtRange(0.14, 4.0);
  fTrackCutsPosPion->SetEtaRange(-0.8, 0.8);
  fTrackCutsPosPion->SetNClsTPC(80);
  fTrackCutsPosPion->SetDCAReCalculation(true);
  fTrackCutsPosPion->SetDCAVtxZ(0.3);
  fTrackCutsPosPion->SetDCAVtxXY(0.3);
  fTrackCutsPosPion->SetPID(AliPID::kPion, 0.5);
  fTrackCutsPosPion->SetRejLowPtPionsTOF(false);
  fTrackCutsPosPion->SetMinimalBooking(false);

  //The same things for negative pions
  fTrackCutsNegPion->SetIsMonteCarlo(isMC);
  fTrackCutsNegPion->SetCutCharge(-1);
  fTrackCutsNegPion->SetFilterBit(96);
  fTrackCutsNegPion->SetPtRange(0.14, 4.0);
  fTrackCutsNegPion->SetEtaRange(-0.8, 0.8);
  fTrackCutsNegPion->SetNClsTPC(80);
  fTrackCutsNegPion->SetDCAReCalculation(true);
  fTrackCutsNegPion->SetDCAVtxZ(0.3);
  fTrackCutsNegPion->SetDCAVtxXY(0.3);
  fTrackCutsNegPion->SetPID(AliPID::kPion, 0.5);
  fTrackCutsNegPion->SetRejLowPtPionsTOF(false);
  fTrackCutsNegPion->SetMinimalBooking(false);
  }

  evtCuts->SetSphericityCuts(fSpherDown, 1.0, ptLow);  


 
  //Now we define stuff we want for our Particle collection
  //Thanks, CINT - will not compile due to an illegal constructor
  //std::vector<int> PDGParticles ={2212,2212,3122,3122,3312,3312};
  //First we need to tell him about the particles we mix, from the
  //PDG code the mass is obtained.
  //Order must mymik the order as the particles are added to the PairCleaner
  std::vector<int> PDGParticles;
  PDGParticles.push_back(211); // pi+
  PDGParticles.push_back(-211); // pi-

  //We need to set the ZVtx bins
  std::vector<float> ZVtxBins;
  ZVtxBins.push_back(-10);
  ZVtxBins.push_back(-8);
  ZVtxBins.push_back(-6);
  ZVtxBins.push_back(-4);
  ZVtxBins.push_back(-2);
  ZVtxBins.push_back(0);
  ZVtxBins.push_back(2);
  ZVtxBins.push_back(4);
  ZVtxBins.push_back(6);
  ZVtxBins.push_back(8);
  ZVtxBins.push_back(10);
  //The Multiplicity bins are set here
  std::vector<int> MultBins;
  MultBins.push_back(0);
  MultBins.push_back(18);
  MultBins.push_back(30);

  //The next part is for the result histograms. The order of hist. is the following:
  //                Particle1     Particle2
  //Particle 1       Hist 1         Hist2
  //
  //Particle 2                      Hist3
  //The same way the values for binning, minimum and maximum k* range have to be set!
  //Number of bins
  std::vector<int> NBins;
  NBins.push_back(750);
  NBins.push_back(750);
  NBins.push_back(750);
  std::vector<float> kMin;
  //minimum k* value
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  //maximum k* value
  std::vector<float> kMax;
  kMax.push_back(3.);
  kMax.push_back(3.);
  kMax.push_back(3.);
  //pair rejection
  std::vector<bool> closeRejection;
  closeRejection.push_back(true); // pi+ pi+
  closeRejection.push_back(false); // pi+ pi- 
  closeRejection.push_back(true); // pi- pi-

  //QA plots for tracks
  std::vector<int> pairQA;
  pairQA.push_back(11); // pi+ pi+
  pairQA.push_back(11); // pi+ pi-
  pairQA.push_back(11); // pi- pi-

  //To put all this into the task we add it to our collection config object in
  //the following way:
  AliFemtoDreamCollConfig *config=new AliFemtoDreamCollConfig("Femto","Femto");
  config->SetZBins(ZVtxBins);
  config->SetMultBins(MultBins);
  //Do you want to have an explicit binning of the correlation function for each multiplicity
  //bin set above?
  config->SetMultBinning(true);
  config->SetPDGCodes(PDGParticles);
  config->SetNBinsHist(NBins);
  config->SetMinKRel(kMin);
  config->SetMaxKRel(kMax);
  config->SetClosePairRejection(closeRejection);
  config->SetDeltaEtaMax(fdPhidEta); // https://alice-notes.web.cern.ch/system/files/notes/analysis/616/2018-08-10-NotepPb.pdf
  config->SetDeltaPhiMax(fdPhidEta);
  config->SetExtendedQAPairs(pairQA);
  //Here we set the mixing depth.
  config->SetMixingDepth(10); // AN
  config->SetkTBinning(true);
  config->SetmTBinning(true);
  config->SetMinimalBookingME(false);
  config->SetdPhidEtaPlots(false);
  config->SetkTandMultBinning(true);
  config->SetdPhidEtaPlotsSmallK(false);
  config->SetPhiEtaBinnign(false);
  
  if (isMC) {
      config->SetMomentumResolution(true);
    } else {
      std::cout << "You are trying to request the Momentum Resolution without MC Info; fix it wont work! \n";
    }
    if (isMC) {
      config->SetPhiEtaBinnign(true);
    } else {
      std::cout << "You are trying to request the Eta Phi Plots without MC Info; fix it wont work! \n";
    }
  
  //now we create the task
  AliAnalysisTaskFemtoDreamPion *task=
      new AliAnalysisTaskFemtoDreamPion("FemtoDreamDefaultPion",isMC);
  //THIS IS VERY IMPORTANT ELSE YOU DONT PROCESS ANY EVENTS
  //kINT7 == Minimum bias
  //kHighMultV0 high multiplicity triggered by the V0 detector
  if(CentEst == "kInt7"){
	task->SetTrigger(AliVEvent::kINT7);
    task->SelectCollisionCandidates(AliVEvent::kINT7);
    std::cout << "Added kINT7 Trigger \n";
  } else if (CentEst == "kHM") {
    task->SelectCollisionCandidates(AliVEvent::kHighMultV0);
    std::cout << "Added kHighMultV0 Trigger \n";
  }else{
    std::cout << "=====================================================================" << std::endl;
    std::cout << "=====================================================================" << std::endl;
    std::cout << "Centrality Estimator not set, fix it else your Results will be empty!" << std::endl;
    std::cout << "=====================================================================" << std::endl;
    std::cout << "=====================================================================" << std::endl;
  }

  //Throw all our settings to the task
  task->SetEventCuts(evtCuts);
  task->SetTrackCutsPosPion(fTrackCutsPosPion);
  task->SetTrackCutsNegPion(fTrackCutsNegPion);
  task->SetCollectionConfig(config);
  task->SetIsMC(isMC);

  mgr->AddTask(task);

  TString file = AliAnalysisManager::GetCommonFileName();

  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();

  mgr->ConnectInput(task, 0, cinput);

  AliAnalysisDataContainer *coutputQA;
  TString addon = "";
  if (CentEst == "kInt7") {
  addon += "MB";
  } else if (CentEst == "kHM") {
  addon += "HM";
  }
  TString QAName = Form("%sResults%s", addon.Data(), suffix.Data());
  coutputQA = mgr->CreateContainer(
    QAName.Data(),
    TList::Class(),
    AliAnalysisManager::kOutputContainer,
    Form("%s:%s", file.Data(), QAName.Data()));
  mgr->ConnectOutput(task, 1, coutputQA);

  return task;
}
