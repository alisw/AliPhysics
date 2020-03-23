AliAnalysisTaskSE* AddTaskNanoAODFemtoDreamPion(
    bool isMC=false, bool MCtemplatefit=false, float fSpherDown=0.7, float fdPhidEta=0.01,
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

  if (suffix == "10") {//pT low
  //Track Cuts are defined here
  //positive pions
  fTrackCutsPosPion->SetIsMonteCarlo(isMC);
  fTrackCutsPosPion->SetCutCharge(1);
  fTrackCutsPosPion->SetFilterBit(96);
  fTrackCutsPosPion->SetPtRange(0.11, 4.0);
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
  fTrackCutsNegPion->SetPtRange(0.11, 4.0);
  fTrackCutsNegPion->SetEtaRange(-0.8, 0.8);
  fTrackCutsNegPion->SetNClsTPC(80);
  fTrackCutsNegPion->SetDCAReCalculation(true);
  fTrackCutsNegPion->SetDCAVtxZ(0.3);
  fTrackCutsNegPion->SetDCAVtxXY(0.3);
  fTrackCutsNegPion->SetPID(AliPID::kPion, 0.5);
  fTrackCutsNegPion->SetRejLowPtPionsTOF(false);
  fTrackCutsNegPion->SetMinimalBooking(false);
  } else if (suffix == "11") {//pT high
  //Track Cuts are defined here
  //positive pions
  fTrackCutsPosPion->SetIsMonteCarlo(isMC);
  fTrackCutsPosPion->SetCutCharge(1);
  fTrackCutsPosPion->SetFilterBit(96);
  fTrackCutsPosPion->SetPtRange(0.14, 5.0);
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
  fTrackCutsNegPion->SetPtRange(0.14, 5.0);
  fTrackCutsNegPion->SetEtaRange(-0.8, 0.8);
  fTrackCutsNegPion->SetNClsTPC(80);
  fTrackCutsNegPion->SetDCAReCalculation(true);
  fTrackCutsNegPion->SetDCAVtxZ(0.3);
  fTrackCutsNegPion->SetDCAVtxXY(0.3);
  fTrackCutsNegPion->SetPID(AliPID::kPion, 0.5);
  fTrackCutsNegPion->SetRejLowPtPionsTOF(false);
  fTrackCutsNegPion->SetMinimalBooking(false);
  } else if ( suffix == 20 ) {//eta low
  //Track Cuts are defined here
  //positive pions
  fTrackCutsPosPion->SetIsMonteCarlo(isMC);
  fTrackCutsPosPion->SetCutCharge(1);
  fTrackCutsPosPion->SetFilterBit(96);
  fTrackCutsPosPion->SetPtRange(0.14, 4.0);
  fTrackCutsPosPion->SetEtaRange(-0.6, 0.6);
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
  fTrackCutsNegPion->SetEtaRange(-0.6, 0.6);
  fTrackCutsNegPion->SetNClsTPC(80);
  fTrackCutsNegPion->SetDCAReCalculation(true);
  fTrackCutsNegPion->SetDCAVtxZ(0.3);
  fTrackCutsNegPion->SetDCAVtxXY(0.3);
  fTrackCutsNegPion->SetPID(AliPID::kPion, 0.5);
  fTrackCutsNegPion->SetRejLowPtPionsTOF(false);
  fTrackCutsNegPion->SetMinimalBooking(false);
  } else if ( suffix == 21 ) {//eta high
  //Track Cuts are defined here
  //positive pions
  fTrackCutsPosPion->SetIsMonteCarlo(isMC);
  fTrackCutsPosPion->SetCutCharge(1);
  fTrackCutsPosPion->SetFilterBit(96);
  fTrackCutsPosPion->SetPtRange(0.14, 4.0);
  fTrackCutsPosPion->SetEtaRange(-0.9, 0.9);
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
  fTrackCutsNegPion->SetEtaRange(-0.9, 0.9);
  fTrackCutsNegPion->SetNClsTPC(80);
  fTrackCutsNegPion->SetDCAReCalculation(true);
  fTrackCutsNegPion->SetDCAVtxZ(0.3);
  fTrackCutsNegPion->SetDCAVtxXY(0.3);
  fTrackCutsNegPion->SetPID(AliPID::kPion, 0.5);
  fTrackCutsNegPion->SetRejLowPtPionsTOF(false);
  fTrackCutsNegPion->SetMinimalBooking(false);
  } else if ( suffix == 22 ) {//eta high conservative
  //Track Cuts are defined here
  //positive pions
  fTrackCutsPosPion->SetIsMonteCarlo(isMC);
  fTrackCutsPosPion->SetCutCharge(1);
  fTrackCutsPosPion->SetFilterBit(96);
  fTrackCutsPosPion->SetPtRange(0.14, 4.0);
  fTrackCutsPosPion->SetEtaRange(-0.85, 0.85);
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
  fTrackCutsNegPion->SetEtaRange(-0.85, 0.85);
  fTrackCutsNegPion->SetNClsTPC(80);
  fTrackCutsNegPion->SetDCAReCalculation(true);
  fTrackCutsNegPion->SetDCAVtxZ(0.3);
  fTrackCutsNegPion->SetDCAVtxXY(0.3);
  fTrackCutsNegPion->SetPID(AliPID::kPion, 0.5);
  fTrackCutsNegPion->SetRejLowPtPionsTOF(false);
  fTrackCutsNegPion->SetMinimalBooking(false);
  } else if ( suffix == 23 ) {//eta low conservative
  //Track Cuts are defined here
  //positive pions
  fTrackCutsPosPion->SetIsMonteCarlo(isMC);
  fTrackCutsPosPion->SetCutCharge(1);
  fTrackCutsPosPion->SetFilterBit(96);
  fTrackCutsPosPion->SetPtRange(0.14, 4.0);
  fTrackCutsPosPion->SetEtaRange(-0.75, 0.75);
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
  fTrackCutsNegPion->SetEtaRange(-0.75, 0.75);
  fTrackCutsNegPion->SetNClsTPC(80);
  fTrackCutsNegPion->SetDCAReCalculation(true);
  fTrackCutsNegPion->SetDCAVtxZ(0.3);
  fTrackCutsNegPion->SetDCAVtxXY(0.3);
  fTrackCutsNegPion->SetPID(AliPID::kPion, 0.5);
  fTrackCutsNegPion->SetRejLowPtPionsTOF(false);
  fTrackCutsNegPion->SetMinimalBooking(false);
  } else if ( suffix == 30 ) {//#TPC Cls low
  //Track Cuts are defined here
  //positive pions
  fTrackCutsPosPion->SetIsMonteCarlo(isMC);
  fTrackCutsPosPion->SetCutCharge(1);
  fTrackCutsPosPion->SetFilterBit(96);
  fTrackCutsPosPion->SetPtRange(0.14, 4.0);
  fTrackCutsPosPion->SetEtaRange(-0.8, 0.8);
  fTrackCutsPosPion->SetNClsTPC(70);
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
  fTrackCutsNegPion->SetNClsTPC(70);
  fTrackCutsNegPion->SetDCAReCalculation(true);
  fTrackCutsNegPion->SetDCAVtxZ(0.3);
  fTrackCutsNegPion->SetDCAVtxXY(0.3);
  fTrackCutsNegPion->SetPID(AliPID::kPion, 0.5);
  fTrackCutsNegPion->SetRejLowPtPionsTOF(false);
  fTrackCutsNegPion->SetMinimalBooking(false);
  } else if ( suffix == 31 ) {//#TPC Cls high conservative
  //Track Cuts are defined here
  //positive pions
  fTrackCutsPosPion->SetIsMonteCarlo(isMC);
  fTrackCutsPosPion->SetCutCharge(1);
  fTrackCutsPosPion->SetFilterBit(96);
  fTrackCutsPosPion->SetPtRange(0.14, 4.0);
  fTrackCutsPosPion->SetEtaRange(-0.8, 0.8);
  fTrackCutsPosPion->SetNClsTPC(85);
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
  fTrackCutsNegPion->SetNClsTPC(85);
  fTrackCutsNegPion->SetDCAReCalculation(true);
  fTrackCutsNegPion->SetDCAVtxZ(0.3);
  fTrackCutsNegPion->SetDCAVtxXY(0.3);
  fTrackCutsNegPion->SetPID(AliPID::kPion, 0.5);
  fTrackCutsNegPion->SetRejLowPtPionsTOF(false);
  fTrackCutsNegPion->SetMinimalBooking(false);
  } else if ( suffix == 32 ) {//#TPC Cls high
  //Track Cuts are defined here
  //positive pions
  fTrackCutsPosPion->SetIsMonteCarlo(isMC);
  fTrackCutsPosPion->SetCutCharge(1);
  fTrackCutsPosPion->SetFilterBit(96);
  fTrackCutsPosPion->SetPtRange(0.14, 4.0);
  fTrackCutsPosPion->SetEtaRange(-0.8, 0.8);
  fTrackCutsPosPion->SetNClsTPC(90);
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
  fTrackCutsNegPion->SetNClsTPC(90);
  fTrackCutsNegPion->SetDCAReCalculation(true);
  fTrackCutsNegPion->SetDCAVtxZ(0.3);
  fTrackCutsNegPion->SetDCAVtxXY(0.3);
  fTrackCutsNegPion->SetPID(AliPID::kPion, 0.5);
  fTrackCutsNegPion->SetRejLowPtPionsTOF(false);
  fTrackCutsNegPion->SetMinimalBooking(false);
  } else if ( suffix == 33 ) {//#TPC Cls low conservative
  //Track Cuts are defined here
  //positive pions
  fTrackCutsPosPion->SetIsMonteCarlo(isMC);
  fTrackCutsPosPion->SetCutCharge(1);
  fTrackCutsPosPion->SetFilterBit(96);
  fTrackCutsPosPion->SetPtRange(0.14, 4.0);
  fTrackCutsPosPion->SetEtaRange(-0.8, 0.8);
  fTrackCutsPosPion->SetNClsTPC(75);
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
  fTrackCutsNegPion->SetNClsTPC(75);
  fTrackCutsNegPion->SetDCAReCalculation(true);
  fTrackCutsNegPion->SetDCAVtxZ(0.3);
  fTrackCutsNegPion->SetDCAVtxXY(0.3);
  fTrackCutsNegPion->SetPID(AliPID::kPion, 0.5);
  fTrackCutsNegPion->SetRejLowPtPionsTOF(false);
  fTrackCutsNegPion->SetMinimalBooking(false);
  } else if ( suffix == 40 ) {//DCA XYZ low conservative
  //Track Cuts are defined here
  //positive pions
  fTrackCutsPosPion->SetIsMonteCarlo(isMC);
  fTrackCutsPosPion->SetCutCharge(1);
  fTrackCutsPosPion->SetFilterBit(96);
  fTrackCutsPosPion->SetPtRange(0.14, 4.0);
  fTrackCutsPosPion->SetEtaRange(-0.8, 0.8);
  fTrackCutsPosPion->SetNClsTPC(80);
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
  fTrackCutsNegPion->SetEtaRange(-0.8, 0.8);
  fTrackCutsNegPion->SetNClsTPC(80);
  fTrackCutsNegPion->SetDCAReCalculation(true);
  fTrackCutsNegPion->SetDCAVtxZ(0.28);
  fTrackCutsNegPion->SetDCAVtxXY(0.28);
  fTrackCutsNegPion->SetPID(AliPID::kPion, 0.5);
  fTrackCutsNegPion->SetRejLowPtPionsTOF(false);
  fTrackCutsNegPion->SetMinimalBooking(false);
  } else if ( suffix == 41 ) {//DCA XYZ high conservative
  //Track Cuts are defined here
  //positive pions
  fTrackCutsPosPion->SetIsMonteCarlo(isMC);
  fTrackCutsPosPion->SetCutCharge(1);
  fTrackCutsPosPion->SetFilterBit(96);
  fTrackCutsPosPion->SetPtRange(0.14, 4.0);
  fTrackCutsPosPion->SetEtaRange(-0.8, 0.8);
  fTrackCutsPosPion->SetNClsTPC(80);
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
  fTrackCutsNegPion->SetPtRange(0.14, 4.0);
  fTrackCutsNegPion->SetEtaRange(-0.8, 0.8);
  fTrackCutsNegPion->SetNClsTPC(80);
  fTrackCutsNegPion->SetDCAReCalculation(true);
  fTrackCutsNegPion->SetDCAVtxZ(0.32);
  fTrackCutsNegPion->SetDCAVtxXY(0.32);
  fTrackCutsNegPion->SetPID(AliPID::kPion, 0.5);
  fTrackCutsNegPion->SetRejLowPtPionsTOF(false);
  fTrackCutsNegPion->SetMinimalBooking(false);
  } else if ( suffix == 42 ) {//DCA XYZ high
  //Track Cuts are defined here
  //positive pions
  fTrackCutsPosPion->SetIsMonteCarlo(isMC);
  fTrackCutsPosPion->SetCutCharge(1);
  fTrackCutsPosPion->SetFilterBit(96);
  fTrackCutsPosPion->SetPtRange(0.14, 4.0);
  fTrackCutsPosPion->SetEtaRange(-0.8, 0.8);
  fTrackCutsPosPion->SetNClsTPC(80);
  fTrackCutsPosPion->SetDCAReCalculation(true);
  fTrackCutsPosPion->SetDCAVtxZ(0.35);
  fTrackCutsPosPion->SetDCAVtxXY(0.35);
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
  fTrackCutsNegPion->SetDCAVtxZ(0.35);
  fTrackCutsNegPion->SetDCAVtxXY(0.35);
  fTrackCutsNegPion->SetPID(AliPID::kPion, 0.5);
  fTrackCutsNegPion->SetRejLowPtPionsTOF(false);
  fTrackCutsNegPion->SetMinimalBooking(false);
  } else if ( suffix == 43 ) {//DCA XYZ low
  //Track Cuts are defined here
  //positive pions
  fTrackCutsPosPion->SetIsMonteCarlo(isMC);
  fTrackCutsPosPion->SetCutCharge(1);
  fTrackCutsPosPion->SetFilterBit(96);
  fTrackCutsPosPion->SetPtRange(0.14, 4.0);
  fTrackCutsPosPion->SetEtaRange(-0.8, 0.8);
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
  fTrackCutsNegPion->SetPtRange(0.14, 4.0);
  fTrackCutsNegPion->SetEtaRange(-0.8, 0.8);
  fTrackCutsNegPion->SetNClsTPC(80);
  fTrackCutsNegPion->SetDCAReCalculation(true);
  fTrackCutsNegPion->SetDCAVtxZ(0.25);
  fTrackCutsNegPion->SetDCAVtxXY(0.25);
  fTrackCutsNegPion->SetPID(AliPID::kPion, 0.5);
  fTrackCutsNegPion->SetRejLowPtPionsTOF(false);
  fTrackCutsNegPion->SetMinimalBooking(false);
  } else if ( suffix == 50 ) {//Sphericity low
  fSpherDown = 0.6;
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
  } else if ( suffix == 51 ) {//Sphericity high
  fSpherDown = 0.8;
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
  } else if ( suffix == 60 ) {//PID low
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
  fTrackCutsPosPion->SetPID(AliPID::kPion, 0.45);
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
  fTrackCutsNegPion->SetPID(AliPID::kPion, 0.45);
  fTrackCutsNegPion->SetRejLowPtPionsTOF(false);
  fTrackCutsNegPion->SetMinimalBooking(false);
  } else if ( suffix == 61 ) {//PID high
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
  fTrackCutsPosPion->SetPID(AliPID::kPion, 0.55);
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
  fTrackCutsNegPion->SetPID(AliPID::kPion, 0.55);
  fTrackCutsNegPion->SetRejLowPtPionsTOF(false);
  fTrackCutsNegPion->SetMinimalBooking(false);
  } else if ( suffix == 70 ) {//CPR low
  fdPhidEta=0.008;
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
  } else if ( suffix == 71 ) {//CPR high
  fdPhidEta=0.012;
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

  evtCuts->SetSphericityCuts(fSpherDown, 1.0);  


  
  //Check proton-pi(-)
  //AliFemtoDreamTrackCuts* fTrackCutsProton = 
  //	AliFemtoDreamTrackCuts::PrimProtonCuts(isMC, false, true, false);

  //Now we define stuff we want for our Particle collection
  //Thanks, CINT - will not compile due to an illegal constructor
  //std::vector<int> PDGParticles ={2212,2212,3122,3122,3312,3312};
  //First we need to tell him about the particles we mix, from the
  //PDG code the mass is obtained.
  //Order must mymik the order as the particles are added to the PairCleaner
  std::vector<int> PDGParticles;
  PDGParticles.push_back(211); // pi+
  PDGParticles.push_back(-211); // pi-
 // PDGParticles.push_back(2212); // p

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
 // NBins.push_back(750);
 // NBins.push_back(750);
 // NBins.push_back(750);
  std::vector<float> kMin;
  //minimum k* value
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
 // kMin.push_back(0.);
 // kMin.push_back(0.);
 // kMin.push_back(0.);
  //maximum k* value
  std::vector<float> kMax;
  kMax.push_back(3.);
  kMax.push_back(3.);
  kMax.push_back(3.);
 // kMax.push_back(3.);
 // kMax.push_back(3.);
 // kMax.push_back(3.);
  //pair rejection
  std::vector<bool> closeRejection;
  closeRejection.push_back(true); // pi+ pi+
  closeRejection.push_back(false); // pi+ pi- 
  closeRejection.push_back(true); // pi- pi-
 // closeRejection.push_back(true); // pi+ p
 // closeRejection.push_back(false); // pi- p 
 // closeRejection.push_back(true); // p p

  if (suffix == "5") {
    //Deactivate the ClosePairRejection
    fdPhidEta=0.;
    closeRejection.clear();
    closeRejection.push_back(false); // pi+ pi+
    closeRejection.push_back(false); // pi+ pi-
    closeRejection.push_back(false); // pi- pi-
  //  closeRejection.push_back(false); // pi+ p
  //  closeRejection.push_back(false); // pi- p
  //  closeRejection.push_back(false); // p p
  }

  //QA plots for tracks
  std::vector<int> pairQA;
  pairQA.push_back(11); // pi+ pi+
  pairQA.push_back(11); // pi+ pi-
  pairQA.push_back(11); // pi- pi-
 // pairQA.push_back(11); // pi+ p
 // pairQA.push_back(11); // pi- p
 // pairQA.push_back(11); // p p

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
  AliAnalysisTaskNanoAODFemtoDreamPion *task=
      new AliAnalysisTaskNanoAODFemtoDreamPion("NanoAODFemtoDreamPion",isMC);
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
  //task->SetTrackCutsProton(fTrackCutsProton);
  task->SetCollectionConfig(config);

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
