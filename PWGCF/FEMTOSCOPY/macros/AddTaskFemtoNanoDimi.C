#ifndef __CINT__
#include <vector>
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskNanoBBar.h"
#include "AliFemtoDreamEventCuts.h"
#include "AliFemtoDreamTrackCuts.h"
#include "AliFemtoDreamCascadeCuts.h"
#include "AliFemtoDreamCollConfig.h"
#endif


AliAnalysisTaskSE *AddTaskFemtoNanoDimi(bool fullBlastQA = false,//1
									 bool isMC = false,				                        //2
									 int fFilterBit = 128,			                      //3
									 TString triggerData = "kInt7",	                  //4
                   TString selectSB = "SL1",                    //4
                   TString mixmethod = "0",
									 TString CutVar = "0"
								 ) {                   //5

  TString suffix = TString::Format("%s_%s_%s",selectSB.Data(),mixmethod.Data(),CutVar.Data());

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskFemtoNanoDimi()", "No analysis manager found.");
    return 0x0;
  }

  // ================== GetInputEventHandler =============================
  AliVEventHandler *inputHandler = mgr->GetInputEventHandler();
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();

  //========= Init subtasks and start analysis ============================
  // Event Cuts
  AliFemtoDreamEventCuts *evtCuts = AliFemtoDreamEventCuts::StandardCutsRun2();
//  void CleanUpMult(bool SPD, bool v0A, bool v0C, bool RefMult) {
  evtCuts->CleanUpMult(false, false, false, true);
  evtCuts->SetDoSpherocityCuts(false);
  evtCuts->SetDoSphericityCuts(false);

  // Proton Cuts
  AliFemtoDreamTrackCuts *TrackCuts = AliFemtoDreamTrackCuts::PrimProtonCuts(
		  isMC, true, false, true);//DCAplots,CombSigma,ContribSplitting
  TrackCuts->SetFilterBit(fFilterBit);
  TrackCuts->SetCutCharge(1);

  AliFemtoDreamTrackCuts *AntiTrackCuts =
      AliFemtoDreamTrackCuts::PrimProtonCuts(isMC, true, false, true);
  AntiTrackCuts->SetFilterBit(fFilterBit);
  AntiTrackCuts->SetCutCharge(-1);

  //Lambda Cuts
  float SidebandLow=0;
  float SidebandUp=0;
  float LambdaWindow=0;
  //SL = sideband left
  //SE = sideband right
  //S# = signal with # MeV around the peak
  if(selectSB=="SL1"){
    SidebandLow = 1.08;
    SidebandUp = 1.103;
  } else if (selectSB=="SL2"){
    SidebandLow = 1.085;
    SidebandUp = 1.103;
  } else if (selectSB=="SR1"){
    SidebandLow = 1.129;
    SidebandUp = 1.155;
  } else if (selectSB=="SR2"){
    SidebandLow = 1.129;
    SidebandUp = 1.2;
  } else if (selectSB=="SL3"){
    SidebandLow = 1.090;
    SidebandUp = 1.103;
  } else if (selectSB=="SL4"){
    SidebandLow = 1.095;
    SidebandUp = 1.108;
  } else if (selectSB=="SR3"){
    SidebandLow = 1.129;
    SidebandUp = 1.145;
  } else if (selectSB=="SR4"){
    SidebandLow = 1.129;
    SidebandUp = 1.140;
  } else if (selectSB=="SR5"){
    SidebandLow = 1.124;
    SidebandUp = 1.140;
  } else if (selectSB=="SR6"){
    SidebandLow = 1.124;
    SidebandUp = 1.135;
  } else if (selectSB=="S40"){//for IMS
	  LambdaWindow=0.04;
  } else if (selectSB=="S4"){//default (peak width c.a. 1.3, so this is c.a. 3 sigma)
	  LambdaWindow=0.004;
  } else if (selectSB=="S2p5"){//2.5 MeV, which is just below 2 sigma of the peak
	  LambdaWindow=0.0025;
  } else{
	LambdaWindow=0.004;
  }

  AliFemtoDreamv0Cuts *v0Cuts = AliFemtoDreamv0Cuts::LambdaCuts(isMC, true, true);
  if(LambdaWindow) v0Cuts->SetCutInvMass(LambdaWindow);
  else v0Cuts->SetCutWindow(SidebandLow,SidebandUp);
  AliFemtoDreamTrackCuts *Posv0Daug = AliFemtoDreamTrackCuts::DecayProtonCuts(isMC, true, false);//PileUpRej, false
  AliFemtoDreamTrackCuts *Negv0Daug = AliFemtoDreamTrackCuts::DecayPionCuts(isMC, true, false);
  v0Cuts->SetPosDaugterTrackCuts(Posv0Daug);
  v0Cuts->SetNegDaugterTrackCuts(Negv0Daug);
  v0Cuts->SetPDGCodePosDaug(2212);  //Proton
  v0Cuts->SetPDGCodeNegDaug(211);  //Pion
  v0Cuts->SetPDGCodev0(3122);  //Lambda

  AliFemtoDreamv0Cuts *Antiv0Cuts = AliFemtoDreamv0Cuts::LambdaCuts(isMC, true, true);
  if(LambdaWindow) Antiv0Cuts->SetCutInvMass(LambdaWindow);
  else Antiv0Cuts->SetCutWindow(SidebandLow,SidebandUp);
  AliFemtoDreamTrackCuts *PosAntiv0Daug = AliFemtoDreamTrackCuts::DecayPionCuts(isMC, true, false);
  PosAntiv0Daug->SetCutCharge(1);
  AliFemtoDreamTrackCuts *NegAntiv0Daug =
      AliFemtoDreamTrackCuts::DecayProtonCuts(isMC, true, false);
  NegAntiv0Daug->SetCutCharge(-1);

  Antiv0Cuts->SetPosDaugterTrackCuts(PosAntiv0Daug);
  Antiv0Cuts->SetNegDaugterTrackCuts(NegAntiv0Daug);
  Antiv0Cuts->SetPDGCodePosDaug(211);  //Pion
  Antiv0Cuts->SetPDGCodeNegDaug(2212);  //Proton
  Antiv0Cuts->SetPDGCodev0(-3122);  //Lambda

    //Cascade Cuts

  AliFemtoDreamCascadeCuts* CascadeCuts = AliFemtoDreamCascadeCuts::XiCuts(
      isMC, false);
  CascadeCuts->SetXiCharge(-1);
    AliFemtoDreamTrackCuts *XiNegCuts = AliFemtoDreamTrackCuts::Xiv0PionCuts(
      isMC, true, false);
  XiNegCuts->SetCheckTPCRefit(false);  //for nanos this is already done while prefiltering
  AliFemtoDreamTrackCuts *XiPosCuts = AliFemtoDreamTrackCuts::Xiv0ProtonCuts(
      isMC, true, false);
  XiPosCuts->SetCheckTPCRefit(false);  //for nanos this is already done while prefiltering
  AliFemtoDreamTrackCuts *XiBachCuts = AliFemtoDreamTrackCuts::XiBachPionCuts(
      isMC, true, false);
  XiBachCuts->SetCheckTPCRefit(false);  //for nanos this is already done while prefiltering

  CascadeCuts->Setv0Negcuts(XiNegCuts);
  CascadeCuts->Setv0PosCuts(XiPosCuts);
  CascadeCuts->SetBachCuts(XiBachCuts);
  CascadeCuts->SetPDGCodeCasc(3312);
  CascadeCuts->SetPDGCodev0(3122);
  CascadeCuts->SetPDGCodePosDaug(2212);
  CascadeCuts->SetPDGCodeNegDaug(-211);
  CascadeCuts->SetPDGCodeBach(-211);

    AliFemtoDreamCascadeCuts* AntiCascadeCuts = AliFemtoDreamCascadeCuts::XiCuts(
      isMC, false);
  AntiCascadeCuts->SetXiCharge(1);
  AliFemtoDreamTrackCuts *AntiXiNegCuts =
      AliFemtoDreamTrackCuts::Xiv0ProtonCuts(isMC, true, false);
  AntiXiNegCuts->SetCutCharge(-1);
  AntiXiNegCuts->SetCheckTPCRefit(false);  //for nanos this is already done while prefiltering
  AliFemtoDreamTrackCuts *AntiXiPosCuts = AliFemtoDreamTrackCuts::Xiv0PionCuts(
      isMC, true, false);
  AntiXiPosCuts->SetCutCharge(1);
  AntiXiPosCuts->SetCheckTPCRefit(false);  //for nanos this is already done while prefiltering
  AliFemtoDreamTrackCuts *AntiXiBachCuts =
      AliFemtoDreamTrackCuts::XiBachPionCuts(isMC, true, false);
  AntiXiBachCuts->SetCutCharge(1);
  AntiXiBachCuts->SetCheckTPCRefit(false);  //for nanos this is already done while prefiltering

  AntiCascadeCuts->Setv0Negcuts(AntiXiNegCuts);
  AntiCascadeCuts->Setv0PosCuts(AntiXiPosCuts);
  AntiCascadeCuts->SetBachCuts(AntiXiBachCuts);
  AntiCascadeCuts->SetPDGCodeCasc(-3312);
  AntiCascadeCuts->SetPDGCodev0(-3122);
  AntiCascadeCuts->SetPDGCodePosDaug(211);
  AntiCascadeCuts->SetPDGCodeNegDaug(-2212);
  AntiCascadeCuts->SetPDGCodeBach(211);


	CascadeCuts->SetMinimalBooking(true);
	AntiCascadeCuts->SetMinimalBooking(true);
  if (!fullBlastQA) {
    evtCuts->SetMinimalBooking(true);
    TrackCuts->SetMinimalBooking(true);
    AntiTrackCuts->SetMinimalBooking(true);
    v0Cuts->SetMinimalBooking(true);
    Antiv0Cuts->SetMinimalBooking(true);
  }


  AliFemtoDreamCollConfig *config = new AliFemtoDreamCollConfig("Femto",
                                                                "Femto", false);
  // Femto Collection
  std::vector<int> PDGParticles;
  PDGParticles.push_back(2212);//p
  PDGParticles.push_back(2212);
  PDGParticles.push_back(3122);//Lambda
  PDGParticles.push_back(3122);
  PDGParticles.push_back(3312);//Cascade
  PDGParticles.push_back(3312);

  std::vector<int> NBins;
  std::vector<float> kMin;
  std::vector<float> kMax;
  std::vector<int> pairQA;
  std::vector<bool> closeRejection;

	//pairs (corrected):
	//pp								0
	//p bar p						1
	//p La							2
	//p bar La					3
	//p Xi							4
	//p bar Xi					5
	//bar p bar p				6
	//bar p La					7
	//bar p bar La			8
	//bar p Xi					9
	//bar p bar Xi			10
	//La La							11
	//La bar La					12
	//La Xi							13
	//La bar Xi					14
	//bar La bar La			15
	//bar La Xi					16
	//bar La bar Xi			17
	//Xi Xi							18
	//Xi bar Xi					19
	//Xi bar Xi bar			20
  const int nPairs = 21;
  for (int i = 0; i < nPairs; ++i) {
    pairQA.push_back(0);
    closeRejection.push_back(false);
    NBins.push_back(384);
    kMin.push_back(0.);
    kMax.push_back(4.608);
  }
  //pairQA[0] = 11;
	//pairQA[6] = 11;
  pairQA[2] = 12;
  pairQA[8] = 12;

	closeRejection[0] = true;  // pp
  closeRejection[6] = true;  // barp barp

  closeRejection[18] = true;  // Xi Xi
  closeRejection[20] = true;  // barXi barXi

  config->SetPDGCodes(PDGParticles);
  config->SetNBinsHist(NBins);
  config->SetMinKRel(kMin);
  config->SetMaxKRel(kMax);
  config->SetClosePairRejection(closeRejection);
  config->SetDeltaEtaMax(0.012);
  config->SetDeltaPhiMax(0.012);
  config->SetExtendedQAPairs(pairQA);

//Setting the configurations of the mixing methods to run on trains:

  if(mixmethod == "0"){
    config->SetMixingDepth(10);
    config->SetUseEventMixing(true);
  } else if (mixmethod == "1a"){
    config->SetUseEventMixing(false);
    config->SetUsePhiSpinning(true);
    config->SetControlMethod(AliFemtoDreamCollConfig::kCorrelatedPhi);
    config->SetCorrelationRange(0.1);
    config->SetSpinningDepth(1);
  } else if (mixmethod == "1b"){
    config->SetUseEventMixing(false);
    config->SetUsePhiSpinning(true);
    config->SetControlMethod(AliFemtoDreamCollConfig::kCorrelatedPhi);
    config->SetCorrelationRange(0.2);
    config->SetSpinningDepth(1);
  } else if (mixmethod == "1c"){
    config->SetUseEventMixing(false);
    config->SetUsePhiSpinning(true);
    config->SetControlMethod(AliFemtoDreamCollConfig::kCorrelatedPhi);
    config->SetCorrelationRange(0.3);
    config->SetSpinningDepth(1);
  } else if (mixmethod == "2"){
    config->SetUseEventMixing(false);
    config->SetUsePhiSpinning(true);
    config->SetControlMethod(AliFemtoDreamCollConfig::kStravinsky);
    config->SetSpinningDepth(1);
  } else if (mixmethod == "3"){
    config->SetUseEventMixing(false);
    config->SetUsePhiSpinning(true);
    config->SetControlMethod(AliFemtoDreamCollConfig::kPhiSpin);
    config->SetSpinningDepth(1);
  }

	bool FetmoBoyzSyst = false;
	bool DimiSyst = false;
	bool DimiSpecial = false;
	int iCutVar = CutVar.Atoi();
	int iCPA = 0;
	int iV0d = 0;
	//normal syst
	if(iCutVar<100){
		FetmoBoyzSyst = true;
	}
	//dimi syst
	else if(iCutVar<200){
		DimiSyst = true;
	}
	//special cuts. Last digit is CPA cut, second to last is the v0 daugter cuts
	//the ones to check: 1000,1001,1002,1003,1004   maybe 1010,1011,1012,1013,1014
	//1000 is the default, so we can skip it
	//all of this should be done both for S4 and S40
	//for S40 -> de we need all data ?
	else if(iCutVar>=1000){
		DimiSpecial = true;
		iCPA = (iCutVar)%10;
		iV0d = (iCutVar/10)%10;
	}

	if(FetmoBoyzSyst){
		if (CutVar == "1") {
	    TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
	    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);

	    TrackCuts->SetNClsTPC(90);
	    AntiTrackCuts->SetNClsTPC(90);

	    v0Cuts->SetCutCPA(0.995);
	    Antiv0Cuts->SetCutCPA(0.995);

	    Posv0Daug->SetNClsTPC(80);
	    Negv0Daug->SetNClsTPC(80);
	    PosAntiv0Daug->SetNClsTPC(80);
	    NegAntiv0Daug->SetNClsTPC(80);

	    v0Cuts->SetCutDCADaugToPrimVtx(0.06);
	    Antiv0Cuts->SetCutDCADaugToPrimVtx(0.06);

	  } else if (CutVar == "2") {
	    TrackCuts->SetPtRange(0.6, 4.05);
	    AntiTrackCuts->SetPtRange(0.6, 4.05);

	    TrackCuts->SetEtaRange(-0.77, 0.77);
	    AntiTrackCuts->SetEtaRange(-0.77, 0.77);

	    config->SetDeltaEtaMax(0.019);
	    config->SetDeltaPhiMax(0.019);

	    Posv0Daug->SetPID(AliPID::kProton, 999.9, 4);
	    Negv0Daug->SetPID(AliPID::kPion, 999.9, 4);
	    PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, 4);
	    NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, 4);

	  } else if (CutVar == "3") {
	    v0Cuts->SetCutCPA(0.995);
	    Antiv0Cuts->SetCutCPA(0.995);

	    Posv0Daug->SetEtaRange(-0.77, 0.77);
	    Negv0Daug->SetEtaRange(-0.77, 0.77);
	    PosAntiv0Daug->SetEtaRange(-0.77, 0.77);
	    NegAntiv0Daug->SetEtaRange(-0.77, 0.77);

	    v0Cuts->SetCutDCADaugTov0Vtx(1.2);
	    Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);

	  } else if (CutVar == "4") {
	    TrackCuts->SetEtaRange(-0.77, 0.77);
	    AntiTrackCuts->SetEtaRange(-0.77, 0.77);

	    TrackCuts->SetNClsTPC(70);
	    AntiTrackCuts->SetNClsTPC(70);

	    Posv0Daug->SetNClsTPC(80);
	    Negv0Daug->SetNClsTPC(80);
	    PosAntiv0Daug->SetNClsTPC(80);
	    NegAntiv0Daug->SetNClsTPC(80);

	    v0Cuts->SetCutDCADaugToPrimVtx(0.06);
	    Antiv0Cuts->SetCutDCADaugToPrimVtx(0.06);

	  } else if (CutVar == "5") {
	    TrackCuts->SetPtRange(0.4, 4.05);
	    AntiTrackCuts->SetPtRange(0.4, 4.05);

	    TrackCuts->SetEtaRange(-0.77, 0.77);
	    AntiTrackCuts->SetEtaRange(-0.77, 0.77);

	    TrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
	    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);

	    config->SetDeltaEtaMax(0.019);
	    config->SetDeltaPhiMax(0.019);

	    v0Cuts->SetCutCPA(0.995);
	    Antiv0Cuts->SetCutCPA(0.995);

	    Posv0Daug->SetEtaRange(-0.83, 0.83);
	    Negv0Daug->SetEtaRange(-0.83, 0.83);
	    PosAntiv0Daug->SetEtaRange(-0.83, 0.83);
	    NegAntiv0Daug->SetEtaRange(-0.83, 0.83);

	    v0Cuts->SetCutDCADaugTov0Vtx(1.2);
	    Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);

	  } else if (CutVar == "6") {
	    TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
	    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);

	    TrackCuts->SetNClsTPC(90);
	    AntiTrackCuts->SetNClsTPC(90);

	    v0Cuts->SetCutCPA(0.995);
	    Antiv0Cuts->SetCutCPA(0.995);

	    Posv0Daug->SetPID(AliPID::kProton, 999.9, 4);
	    Negv0Daug->SetPID(AliPID::kPion, 999.9, 4);
	    PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, 4);
	    NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, 4);

	    Posv0Daug->SetEtaRange(-0.83, 0.83);
	    Negv0Daug->SetEtaRange(-0.83, 0.83);
	    PosAntiv0Daug->SetEtaRange(-0.83, 0.83);
	    NegAntiv0Daug->SetEtaRange(-0.83, 0.83);

	  } else if (CutVar == "7") {
	    TrackCuts->SetEtaRange(-0.85, 0.85);
	    AntiTrackCuts->SetEtaRange(-0.85, 0.85);

	    TrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
	    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);

	    config->SetDeltaEtaMax(0.019);
	    config->SetDeltaPhiMax(0.019);

	    v0Cuts->SetCutCPA(0.995);
	    Antiv0Cuts->SetCutCPA(0.995);

	    Posv0Daug->SetEtaRange(-0.83, 0.83);
	    Negv0Daug->SetEtaRange(-0.83, 0.83);
	    PosAntiv0Daug->SetEtaRange(-0.83, 0.83);
	    NegAntiv0Daug->SetEtaRange(-0.83, 0.83);

	    v0Cuts->SetCutDCADaugTov0Vtx(1.2);
	    Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);

	    v0Cuts->SetCutDCADaugToPrimVtx(0.06);
	    Antiv0Cuts->SetCutDCADaugToPrimVtx(0.06);

	  } else if (CutVar == "8") {
	    TrackCuts->SetPtRange(0.4, 4.05);
	    AntiTrackCuts->SetPtRange(0.4, 4.05);

	    TrackCuts->SetEtaRange(-0.77, 0.77);
	    AntiTrackCuts->SetEtaRange(-0.77, 0.77);

	    TrackCuts->SetNClsTPC(90);
	    AntiTrackCuts->SetNClsTPC(90);

	    config->SetDeltaEtaMax(0.019);
	    config->SetDeltaPhiMax(0.019);

	    Posv0Daug->SetEtaRange(-0.83, 0.83);
	    Negv0Daug->SetEtaRange(-0.83, 0.83);
	    PosAntiv0Daug->SetEtaRange(-0.83, 0.83);
	    NegAntiv0Daug->SetEtaRange(-0.83, 0.83);

	  } else if (CutVar == "9") {
	    Posv0Daug->SetNClsTPC(80);
	    Negv0Daug->SetNClsTPC(80);
	    PosAntiv0Daug->SetNClsTPC(80);
	    NegAntiv0Daug->SetNClsTPC(80);

	    Posv0Daug->SetEtaRange(-0.77, 0.77);
	    Negv0Daug->SetEtaRange(-0.77, 0.77);
	    PosAntiv0Daug->SetEtaRange(-0.77, 0.77);
	    NegAntiv0Daug->SetEtaRange(-0.77, 0.77);

	    v0Cuts->SetCutDCADaugTov0Vtx(1.2);
	    Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);

	  } else if (CutVar == "10") {
	    TrackCuts->SetPtRange(0.4, 4.05);
	    AntiTrackCuts->SetPtRange(0.4, 4.05);

	    TrackCuts->SetEtaRange(-0.85, 0.85);
	    AntiTrackCuts->SetEtaRange(-0.85, 0.85);

	    Posv0Daug->SetPID(AliPID::kProton, 999.9, 4);
	    Negv0Daug->SetPID(AliPID::kPion, 999.9, 4);
	    PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, 4);
	    NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, 4);

	    Posv0Daug->SetEtaRange(-0.77, 0.77);
	    Negv0Daug->SetEtaRange(-0.77, 0.77);
	    PosAntiv0Daug->SetEtaRange(-0.77, 0.77);
	    NegAntiv0Daug->SetEtaRange(-0.77, 0.77);

	    v0Cuts->SetCutDCADaugToPrimVtx(0.06);
	    Antiv0Cuts->SetCutDCADaugToPrimVtx(0.06);

	  } else if (CutVar == "11") {
	    TrackCuts->SetEtaRange(-0.77, 0.77);
	    AntiTrackCuts->SetEtaRange(-0.77, 0.77);

	    TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
	    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);

	    Posv0Daug->SetNClsTPC(80);
	    Negv0Daug->SetNClsTPC(80);
	    PosAntiv0Daug->SetNClsTPC(80);
	    NegAntiv0Daug->SetNClsTPC(80);

	    Posv0Daug->SetEtaRange(-0.77, 0.77);
	    Negv0Daug->SetEtaRange(-0.77, 0.77);
	    PosAntiv0Daug->SetEtaRange(-0.77, 0.77);
	    NegAntiv0Daug->SetEtaRange(-0.77, 0.77);

	    v0Cuts->SetCutDCADaugTov0Vtx(1.2);
	    Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);

	  } else if (CutVar == "12") {
	    TrackCuts->SetPtRange(0.4, 4.05);
	    AntiTrackCuts->SetPtRange(0.4, 4.05);

	    TrackCuts->SetEtaRange(-0.85, 0.85);
	    AntiTrackCuts->SetEtaRange(-0.85, 0.85);

	    TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
	    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);

	    config->SetDeltaEtaMax(0.019);
	    config->SetDeltaPhiMax(0.019);

	    Posv0Daug->SetPID(AliPID::kProton, 999.9, 4);
	    Negv0Daug->SetPID(AliPID::kPion, 999.9, 4);
	    PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, 4);
	    NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, 4);

	    Posv0Daug->SetEtaRange(-0.77, 0.77);
	    Negv0Daug->SetEtaRange(-0.77, 0.77);
	    PosAntiv0Daug->SetEtaRange(-0.77, 0.77);
	    NegAntiv0Daug->SetEtaRange(-0.77, 0.77);

	  } else if (CutVar == "13") {
	    TrackCuts->SetEtaRange(-0.85, 0.85);
	    AntiTrackCuts->SetEtaRange(-0.85, 0.85);

	    TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
	    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);

	    TrackCuts->SetNClsTPC(90);
	    AntiTrackCuts->SetNClsTPC(90);

	    Posv0Daug->SetEtaRange(-0.83, 0.83);
	    Negv0Daug->SetEtaRange(-0.83, 0.83);
	    PosAntiv0Daug->SetEtaRange(-0.83, 0.83);
	    NegAntiv0Daug->SetEtaRange(-0.83, 0.83);

	    v0Cuts->SetCutDCADaugTov0Vtx(1.2);
	    Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);

	    v0Cuts->SetCutDCADaugToPrimVtx(0.06);
	    Antiv0Cuts->SetCutDCADaugToPrimVtx(0.06);

	  } else if (CutVar == "14") {
	    TrackCuts->SetPtRange(0.4, 4.05);
	    AntiTrackCuts->SetPtRange(0.4, 4.05);

	    TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
	    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);

	    config->SetDeltaEtaMax(0.019);
	    config->SetDeltaPhiMax(0.019);

	    Posv0Daug->SetNClsTPC(80);
	    Negv0Daug->SetNClsTPC(80);
	    PosAntiv0Daug->SetNClsTPC(80);
	    NegAntiv0Daug->SetNClsTPC(80);

	    Posv0Daug->SetEtaRange(-0.83, 0.83);
	    Negv0Daug->SetEtaRange(-0.83, 0.83);
	    PosAntiv0Daug->SetEtaRange(-0.83, 0.83);
	    NegAntiv0Daug->SetEtaRange(-0.83, 0.83);

	  } else if (CutVar == "15") {
	    TrackCuts->SetPtRange(0.6, 4.05);
	    AntiTrackCuts->SetPtRange(0.6, 4.05);

	    TrackCuts->SetEtaRange(-0.85, 0.85);
	    AntiTrackCuts->SetEtaRange(-0.85, 0.85);

	    TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
	    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);

	    TrackCuts->SetNClsTPC(90);
	    AntiTrackCuts->SetNClsTPC(90);

	    v0Cuts->SetCutCPA(0.995);
	    Antiv0Cuts->SetCutCPA(0.995);

	    Posv0Daug->SetNClsTPC(80);
	    Negv0Daug->SetNClsTPC(80);
	    PosAntiv0Daug->SetNClsTPC(80);
	    NegAntiv0Daug->SetNClsTPC(80);

	    v0Cuts->SetCutDCADaugTov0Vtx(1.2);
	    Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);

	    v0Cuts->SetCutDCADaugToPrimVtx(0.06);
	    Antiv0Cuts->SetCutDCADaugToPrimVtx(0.06);

	  } else if (CutVar == "16") {
	    TrackCuts->SetPtRange(0.4, 4.05);
	    AntiTrackCuts->SetPtRange(0.4, 4.05);

	    TrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
	    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);

	    TrackCuts->SetNClsTPC(70);
	    AntiTrackCuts->SetNClsTPC(70);

	    config->SetDeltaEtaMax(0.019);
	    config->SetDeltaPhiMax(0.019);

	    Posv0Daug->SetNClsTPC(80);
	    Negv0Daug->SetNClsTPC(80);
	    PosAntiv0Daug->SetNClsTPC(80);
	    NegAntiv0Daug->SetNClsTPC(80);

	    Posv0Daug->SetEtaRange(-0.77, 0.77);
	    Negv0Daug->SetEtaRange(-0.77, 0.77);
	    PosAntiv0Daug->SetEtaRange(-0.77, 0.77);
	    NegAntiv0Daug->SetEtaRange(-0.77, 0.77);

	    v0Cuts->SetCutDCADaugToPrimVtx(0.06);
	    Antiv0Cuts->SetCutDCADaugToPrimVtx(0.06);

	  } else if (CutVar == "17") {
	    TrackCuts->SetPtRange(0.4, 4.05);
	    AntiTrackCuts->SetPtRange(0.4, 4.05);

	    TrackCuts->SetEtaRange(-0.85, 0.85);
	    AntiTrackCuts->SetEtaRange(-0.85, 0.85);

	    TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
	    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);

	    config->SetDeltaEtaMax(0.019);
	    config->SetDeltaPhiMax(0.019);

	    v0Cuts->SetCutCPA(0.995);
	    Antiv0Cuts->SetCutCPA(0.995);

	    Posv0Daug->SetEtaRange(-0.77, 0.77);
	    Negv0Daug->SetEtaRange(-0.77, 0.77);
	    PosAntiv0Daug->SetEtaRange(-0.77, 0.77);
	    NegAntiv0Daug->SetEtaRange(-0.77, 0.77);

	    v0Cuts->SetCutDCADaugToPrimVtx(0.06);
	    Antiv0Cuts->SetCutDCADaugToPrimVtx(0.06);

	  } else if (CutVar == "18") {
	    TrackCuts->SetPtRange(0.4, 4.05);
	    AntiTrackCuts->SetPtRange(0.4, 4.05);

	    TrackCuts->SetEtaRange(-0.85, 0.85);
	    AntiTrackCuts->SetEtaRange(-0.85, 0.85);

	    TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
	    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);

	    TrackCuts->SetNClsTPC(90);
	    AntiTrackCuts->SetNClsTPC(90);

	    Posv0Daug->SetNClsTPC(80);
	    Negv0Daug->SetNClsTPC(80);
	    PosAntiv0Daug->SetNClsTPC(80);
	    NegAntiv0Daug->SetNClsTPC(80);

	    Posv0Daug->SetEtaRange(-0.77, 0.77);
	    Negv0Daug->SetEtaRange(-0.77, 0.77);
	    PosAntiv0Daug->SetEtaRange(-0.77, 0.77);
	    NegAntiv0Daug->SetEtaRange(-0.77, 0.77);

	  } else if (CutVar == "19") {
	    TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
	    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);

	    TrackCuts->SetNClsTPC(70);
	    AntiTrackCuts->SetNClsTPC(70);

	    config->SetDeltaEtaMax(0.019);
	    config->SetDeltaPhiMax(0.019);

	    v0Cuts->SetCutCPA(0.995);
	    Antiv0Cuts->SetCutCPA(0.995);

	    Posv0Daug->SetNClsTPC(80);
	    Negv0Daug->SetNClsTPC(80);
	    PosAntiv0Daug->SetNClsTPC(80);
	    NegAntiv0Daug->SetNClsTPC(80);

	    v0Cuts->SetCutDCADaugTov0Vtx(1.2);
	    Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);

	  } else if (CutVar == "20") {
	    TrackCuts->SetPtRange(0.4, 4.05);
	    AntiTrackCuts->SetPtRange(0.4, 4.05);

	    TrackCuts->SetNClsTPC(90);
	    AntiTrackCuts->SetNClsTPC(90);

	    config->SetDeltaEtaMax(0.019);
	    config->SetDeltaPhiMax(0.019);

	    Posv0Daug->SetNClsTPC(80);
	    Negv0Daug->SetNClsTPC(80);
	    PosAntiv0Daug->SetNClsTPC(80);
	    NegAntiv0Daug->SetNClsTPC(80);

	    Posv0Daug->SetEtaRange(-0.83, 0.83);
	    Negv0Daug->SetEtaRange(-0.83, 0.83);
	    PosAntiv0Daug->SetEtaRange(-0.83, 0.83);
	    NegAntiv0Daug->SetEtaRange(-0.83, 0.83);

	  } else if (CutVar == "21") {

	    TrackCuts->SetEtaRange(-0.85, 0.85);
	    AntiTrackCuts->SetEtaRange(-0.85, 0.85);

	    TrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
	    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);

	    TrackCuts->SetNClsTPC(70);
	    AntiTrackCuts->SetNClsTPC(70);

	    config->SetDeltaEtaMax(0.019);
	    config->SetDeltaPhiMax(0.019);

	    Posv0Daug->SetNClsTPC(80);
	    Negv0Daug->SetNClsTPC(80);
	    PosAntiv0Daug->SetNClsTPC(80);
	    NegAntiv0Daug->SetNClsTPC(80);

	    Posv0Daug->SetEtaRange(-0.77, 0.77);
	    Negv0Daug->SetEtaRange(-0.77, 0.77);
	    PosAntiv0Daug->SetEtaRange(-0.77, 0.77);
	    NegAntiv0Daug->SetEtaRange(-0.77, 0.77);

	  } else if (CutVar == "22") {
	    TrackCuts->SetPtRange(0.4, 4.05);
	    AntiTrackCuts->SetPtRange(0.4, 4.05);

	    TrackCuts->SetEtaRange(-0.77, 0.77);
	    AntiTrackCuts->SetEtaRange(-0.77, 0.77);

	    TrackCuts->SetNClsTPC(70);
	    AntiTrackCuts->SetNClsTPC(70);

	    config->SetDeltaEtaMax(0.019);
	    config->SetDeltaPhiMax(0.019);

	    v0Cuts->SetCutCPA(0.995);
	    Antiv0Cuts->SetCutCPA(0.995);

	    Posv0Daug->SetPID(AliPID::kProton, 999.9, 4);
	    Negv0Daug->SetPID(AliPID::kPion, 999.9, 4);
	    PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, 4);
	    NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, 4);

	  } else if (CutVar == "23") {
	    TrackCuts->SetEtaRange(-0.77, 0.77);
	    AntiTrackCuts->SetEtaRange(-0.77, 0.77);

	    TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
	    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);

	    Posv0Daug->SetPID(AliPID::kProton, 999.9, 4);
	    Negv0Daug->SetPID(AliPID::kPion, 999.9, 4);
	    PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, 4);
	    NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, 4);

	    Posv0Daug->SetEtaRange(-0.83, 0.83);
	    Negv0Daug->SetEtaRange(-0.83, 0.83);
	    PosAntiv0Daug->SetEtaRange(-0.83, 0.83);
	    NegAntiv0Daug->SetEtaRange(-0.83, 0.83);

	    v0Cuts->SetCutDCADaugTov0Vtx(1.2);
	    Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);

	  } else if (CutVar == "24") {
	    TrackCuts->SetEtaRange(-0.77, 0.77);
	    AntiTrackCuts->SetEtaRange(-0.77, 0.77);

	    TrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
	    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);

	    TrackCuts->SetNClsTPC(90);
	    AntiTrackCuts->SetNClsTPC(90);

	    config->SetDeltaEtaMax(0.019);
	    config->SetDeltaPhiMax(0.019);

	    v0Cuts->SetCutCPA(0.995);
	    Antiv0Cuts->SetCutCPA(0.995);

	    Posv0Daug->SetPID(AliPID::kProton, 999.9, 4);
	    Negv0Daug->SetPID(AliPID::kPion, 999.9, 4);
	    PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, 4);
	    NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, 4);

	    Posv0Daug->SetEtaRange(-0.83, 0.83);
	    Negv0Daug->SetEtaRange(-0.83, 0.83);
	    PosAntiv0Daug->SetEtaRange(-0.83, 0.83);
	    NegAntiv0Daug->SetEtaRange(-0.83, 0.83);

	    //XI

	  } else if (CutVar == "25") {
	    TrackCuts->SetPtRange(0.4, 4.05);
	    AntiTrackCuts->SetPtRange(0.4, 4.05);

	    TrackCuts->SetEtaRange(-0.83, 0.83);
	    AntiTrackCuts->SetEtaRange(-0.83, 0.83);

	    config->SetDeltaEtaMax(0.019);
	    config->SetDeltaPhiMax(0.019);

	    Posv0Daug->SetPID(AliPID::kProton, 999.9, 4);
	    Negv0Daug->SetPID(AliPID::kPion, 999.9, 4);
	    PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, 4);
	    NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, 4);

	    Posv0Daug->SetEtaRange(-0.77, 0.77);
	    Negv0Daug->SetEtaRange(-0.77, 0.77);
	    PosAntiv0Daug->SetEtaRange(-0.77, 0.77);
	    NegAntiv0Daug->SetEtaRange(-0.77, 0.77);

	    v0Cuts->SetCutDCADaugToPrimVtx(0.06);
	    Antiv0Cuts->SetCutDCADaugToPrimVtx(0.06);

	  } else if (CutVar == "26") {
	    TrackCuts->SetPtRange(0.6, 4.05);
	    AntiTrackCuts->SetPtRange(0.6, 4.05);

	    TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
	    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);

	    config->SetDeltaEtaMax(0.019);
	    config->SetDeltaPhiMax(0.019);

	    Posv0Daug->SetEtaRange(-0.77, 0.77);
	    Negv0Daug->SetEtaRange(-0.77, 0.77);
	    PosAntiv0Daug->SetEtaRange(-0.77, 0.77);
	    NegAntiv0Daug->SetEtaRange(-0.77, 0.77);

	  } else if (CutVar == "27") {
	    TrackCuts->SetEtaRange(-0.77, 0.77);
	    AntiTrackCuts->SetEtaRange(-0.77, 0.77);

	    TrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
	    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);

	    TrackCuts->SetNClsTPC(90);
	    AntiTrackCuts->SetNClsTPC(90);

	    Posv0Daug->SetPID(AliPID::kProton, 999.9, 4);
	    Negv0Daug->SetPID(AliPID::kPion, 999.9, 4);
	    PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, 4);
	    NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, 4);

	    Posv0Daug->SetEtaRange(-0.77, 0.77);
	    Negv0Daug->SetEtaRange(-0.77, 0.77);
	    PosAntiv0Daug->SetEtaRange(-0.77, 0.77);
	    NegAntiv0Daug->SetEtaRange(-0.77, 0.77);

	    v0Cuts->SetCutDCADaugTov0Vtx(1.2);
	    Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);

	  } else if (CutVar == "28") {
	    TrackCuts->SetNClsTPC(90);
	    AntiTrackCuts->SetNClsTPC(90);

	    v0Cuts->SetCutCPA(0.995);
	    Antiv0Cuts->SetCutCPA(0.995);

	    v0Cuts->SetCutDCADaugToPrimVtx(0.06);
	    Antiv0Cuts->SetCutDCADaugToPrimVtx(0.06);

	  } else if (CutVar == "29") {
	    TrackCuts->SetPtRange(0.4, 4.05);
	    AntiTrackCuts->SetPtRange(0.4, 4.05);

	    TrackCuts->SetEtaRange(-0.77, 0.77);
	    AntiTrackCuts->SetEtaRange(-0.77, 0.77);

	    TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
	    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);

	    config->SetDeltaEtaMax(0.019);
	    config->SetDeltaPhiMax(0.019);

	    Posv0Daug->SetNClsTPC(80);
	    Negv0Daug->SetNClsTPC(80);
	    PosAntiv0Daug->SetNClsTPC(80);
	    NegAntiv0Daug->SetNClsTPC(80);

	    v0Cuts->SetCutDCADaugToPrimVtx(0.06);
	    Antiv0Cuts->SetCutDCADaugToPrimVtx(0.06);
	  } else if (CutVar == "30") {
	    TrackCuts->SetEtaRange(-0.85, 0.85);
	    AntiTrackCuts->SetEtaRange(-0.85, 0.85);

	    config->SetDeltaEtaMax(0.019);
	    config->SetDeltaPhiMax(0.019);

	    v0Cuts->SetCutCPA(0.995);
	    Antiv0Cuts->SetCutCPA(0.995);

	    Posv0Daug->SetPID(AliPID::kProton, 999.9, 4);
	    Negv0Daug->SetPID(AliPID::kPion, 999.9, 4);
	    PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, 4);
	    NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, 4);

	    Posv0Daug->SetNClsTPC(80);
	    Negv0Daug->SetNClsTPC(80);
	    PosAntiv0Daug->SetNClsTPC(80);
	    NegAntiv0Daug->SetNClsTPC(80);

	    Posv0Daug->SetEtaRange(-0.83, 0.83);
	    Negv0Daug->SetEtaRange(-0.83, 0.83);
	    PosAntiv0Daug->SetEtaRange(-0.83, 0.83);
	    NegAntiv0Daug->SetEtaRange(-0.83, 0.83);

	    v0Cuts->SetCutDCADaugTov0Vtx(1.2);
	    Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);
	  } else if (CutVar == "31") {
	    TrackCuts->SetPtRange(0.4, 4.05);
	    AntiTrackCuts->SetPtRange(0.4, 4.05);

	    TrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
	    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);

	    TrackCuts->SetNClsTPC(70);
	    AntiTrackCuts->SetNClsTPC(70);

	    config->SetDeltaEtaMax(0.019);
	    config->SetDeltaPhiMax(0.019);

	    v0Cuts->SetCutCPA(0.995);
	    Antiv0Cuts->SetCutCPA(0.995);

	    Posv0Daug->SetEtaRange(-0.83, 0.83);
	    Negv0Daug->SetEtaRange(-0.83, 0.83);
	    PosAntiv0Daug->SetEtaRange(-0.83, 0.83);
	    NegAntiv0Daug->SetEtaRange(-0.83, 0.83);

	    v0Cuts->SetCutDCADaugTov0Vtx(1.2);
	    Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);
	  } else if (CutVar == "32") {

	    TrackCuts->SetEtaRange(-0.77, 0.77);
	    AntiTrackCuts->SetEtaRange(-0.77, 0.77);

	    TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
	    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);

	    TrackCuts->SetNClsTPC(90);
	    AntiTrackCuts->SetNClsTPC(90);

	    v0Cuts->SetCutCPA(0.995);
	    Antiv0Cuts->SetCutCPA(0.995);

	    Posv0Daug->SetNClsTPC(80);
	    Negv0Daug->SetNClsTPC(80);
	    PosAntiv0Daug->SetNClsTPC(80);
	    NegAntiv0Daug->SetNClsTPC(80);

	    v0Cuts->SetCutDCADaugToPrimVtx(0.06);
	    Antiv0Cuts->SetCutDCADaugToPrimVtx(0.06);
	  } else if (CutVar == "33") {
	    TrackCuts->SetEtaRange(-0.77, 0.77);
	    AntiTrackCuts->SetEtaRange(-0.77, 0.77);

	    config->SetDeltaEtaMax(0.019);
	    config->SetDeltaPhiMax(0.019);

	    v0Cuts->SetCutCPA(0.995);
	    Antiv0Cuts->SetCutCPA(0.995);

	    Posv0Daug->SetNClsTPC(80);
	    Negv0Daug->SetNClsTPC(80);
	    PosAntiv0Daug->SetNClsTPC(80);
	    NegAntiv0Daug->SetNClsTPC(80);

	    Posv0Daug->SetEtaRange(-0.77, 0.77);
	    Negv0Daug->SetEtaRange(-0.77, 0.77);
	    PosAntiv0Daug->SetEtaRange(-0.77, 0.77);
	    NegAntiv0Daug->SetEtaRange(-0.77, 0.77);
	  } else if (CutVar == "34") {
	    TrackCuts->SetEtaRange(-0.77, 0.77);
	    AntiTrackCuts->SetEtaRange(-0.77, 0.77);

	    TrackCuts->SetNClsTPC(70);
	    AntiTrackCuts->SetNClsTPC(70);

	    config->SetDeltaEtaMax(0.019);
	    config->SetDeltaPhiMax(0.019);

	    Posv0Daug->SetNClsTPC(80);
	    Negv0Daug->SetNClsTPC(80);
	    PosAntiv0Daug->SetNClsTPC(80);
	    NegAntiv0Daug->SetNClsTPC(80);

	    Posv0Daug->SetEtaRange(-0.77, 0.77);
	    Negv0Daug->SetEtaRange(-0.77, 0.77);
	    PosAntiv0Daug->SetEtaRange(-0.77, 0.77);
	    NegAntiv0Daug->SetEtaRange(-0.77, 0.77);

	  } else if (CutVar == "35") {
	    TrackCuts->SetPtRange(0.4, 4.05);
	    AntiTrackCuts->SetPtRange(0.4, 4.05);

	    TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
	    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);

	    TrackCuts->SetNClsTPC(70);
	    AntiTrackCuts->SetNClsTPC(70);

	    config->SetDeltaEtaMax(0.019);
	    config->SetDeltaPhiMax(0.019);

	    v0Cuts->SetCutCPA(0.995);
	    Antiv0Cuts->SetCutCPA(0.995);

	    Posv0Daug->SetNClsTPC(80);
	    Negv0Daug->SetNClsTPC(80);
	    PosAntiv0Daug->SetNClsTPC(80);
	    NegAntiv0Daug->SetNClsTPC(80);

	    v0Cuts->SetCutDCADaugTov0Vtx(1.2);
	    Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);
	  } else if (CutVar == "36") {
	    TrackCuts->SetPtRange(0.6, 4.05);
	    AntiTrackCuts->SetPtRange(0.6, 4.05);

	    TrackCuts->SetEtaRange(-0.85, 0.85);
	    AntiTrackCuts->SetEtaRange(-0.85, 0.85);

	    TrackCuts->SetNClsTPC(70);
	    AntiTrackCuts->SetNClsTPC(70);

	    config->SetDeltaEtaMax(0.019);
	    config->SetDeltaPhiMax(0.019);

	    Posv0Daug->SetPID(AliPID::kProton, 999.9, 4);
	    Negv0Daug->SetPID(AliPID::kPion, 999.9, 4);
	    PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, 4);
	    NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, 4);

	    Posv0Daug->SetNClsTPC(80);
	    Negv0Daug->SetNClsTPC(80);
	    PosAntiv0Daug->SetNClsTPC(80);
	    NegAntiv0Daug->SetNClsTPC(80);

	    Posv0Daug->SetEtaRange(-0.83, 0.83);
	    Negv0Daug->SetEtaRange(-0.83, 0.83);
	    PosAntiv0Daug->SetEtaRange(-0.83, 0.83);
	    NegAntiv0Daug->SetEtaRange(-0.83, 0.83);

	    v0Cuts->SetCutDCADaugTov0Vtx(1.2);
	    Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);
	  } else if (CutVar == "37") {
	    TrackCuts->SetEtaRange(-0.85, 0.85);
	    AntiTrackCuts->SetEtaRange(-0.85, 0.85);

	    TrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
	    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);

	    TrackCuts->SetNClsTPC(90);
	    AntiTrackCuts->SetNClsTPC(90);

	    config->SetDeltaEtaMax(0.019);
	    config->SetDeltaPhiMax(0.019);

	    v0Cuts->SetCutCPA(0.995);
	    Antiv0Cuts->SetCutCPA(0.995);

	    Posv0Daug->SetEtaRange(-0.83, 0.83);
	    Negv0Daug->SetEtaRange(-0.83, 0.83);
	    PosAntiv0Daug->SetEtaRange(-0.83, 0.83);
	    NegAntiv0Daug->SetEtaRange(-0.83, 0.83);

	    v0Cuts->SetCutDCADaugTov0Vtx(1.2);
	    Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);

	  } else if (CutVar == "38") {
	    TrackCuts->SetEtaRange(-0.77, 0.77);
	    AntiTrackCuts->SetEtaRange(-0.77, 0.77);

	    TrackCuts->SetNClsTPC(90);
	    AntiTrackCuts->SetNClsTPC(90);

	    config->SetDeltaEtaMax(0.019);
	    config->SetDeltaPhiMax(0.019);

	    v0Cuts->SetCutCPA(0.995);
	    Antiv0Cuts->SetCutCPA(0.995);

	    Posv0Daug->SetNClsTPC(80);
	    Negv0Daug->SetNClsTPC(80);
	    PosAntiv0Daug->SetNClsTPC(80);
	    NegAntiv0Daug->SetNClsTPC(80);

	    Posv0Daug->SetEtaRange(-0.77, 0.77);
	    Negv0Daug->SetEtaRange(-0.77, 0.77);
	    PosAntiv0Daug->SetEtaRange(-0.77, 0.77);
	    NegAntiv0Daug->SetEtaRange(-0.77, 0.77);
	  } else if (CutVar == "39") {
	    TrackCuts->SetPtRange(0.4, 4.05);
	    AntiTrackCuts->SetPtRange(0.4, 4.05);

	    TrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
	    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);

	    TrackCuts->SetNClsTPC(90);
	    AntiTrackCuts->SetNClsTPC(90);

	    config->SetDeltaEtaMax(0.019);
	    config->SetDeltaPhiMax(0.019);

	    Posv0Daug->SetPID(AliPID::kProton, 999.9, 4);
	    Negv0Daug->SetPID(AliPID::kPion, 999.9, 4);
	    PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, 4);
	    NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, 4);

	    Posv0Daug->SetNClsTPC(80);
	    Negv0Daug->SetNClsTPC(80);
	    PosAntiv0Daug->SetNClsTPC(80);
	    NegAntiv0Daug->SetNClsTPC(80);

	    Posv0Daug->SetEtaRange(-0.77, 0.77);
	    Negv0Daug->SetEtaRange(-0.77, 0.77);
	    PosAntiv0Daug->SetEtaRange(-0.77, 0.77);
	    NegAntiv0Daug->SetEtaRange(-0.77, 0.77);

	    v0Cuts->SetCutDCADaugTov0Vtx(1.2);
	    Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);
	  } else if (CutVar == "40") {
	    TrackCuts->SetPtRange(0.6, 4.05);
	    AntiTrackCuts->SetPtRange(0.6, 4.05);

	    TrackCuts->SetNClsTPC(70);
	    AntiTrackCuts->SetNClsTPC(70);

	    Posv0Daug->SetPID(AliPID::kProton, 999.9, 4);
	    Negv0Daug->SetPID(AliPID::kPion, 999.9, 4);
	    PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, 4);
	    NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, 4);

	    v0Cuts->SetCutDCADaugTov0Vtx(1.2);
	    Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);
	  } else if (CutVar == "41") {
	    TrackCuts->SetPtRange(0.4, 4.05);
	    AntiTrackCuts->SetPtRange(0.4, 4.05);

	    TrackCuts->SetEtaRange(-0.77, 0.77);
	    AntiTrackCuts->SetEtaRange(-0.77, 0.77);

	    TrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
	    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);

	    config->SetDeltaEtaMax(0.019);
	    config->SetDeltaPhiMax(0.019);

	    v0Cuts->SetCutCPA(0.995);
	    Antiv0Cuts->SetCutCPA(0.995);

	    Posv0Daug->SetEtaRange(-0.77, 0.77);
	    Negv0Daug->SetEtaRange(-0.77, 0.77);
	    PosAntiv0Daug->SetEtaRange(-0.77, 0.77);
	    NegAntiv0Daug->SetEtaRange(-0.77, 0.77);

	    v0Cuts->SetCutDCADaugTov0Vtx(1.2);
	    Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);
	  } else if (CutVar == "42") {
	    TrackCuts->SetPtRange(0.6, 4.05);
	    AntiTrackCuts->SetPtRange(0.6, 4.05);

	    TrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);
	    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 2.5);

	    v0Cuts->SetCutCPA(0.995);
	    Antiv0Cuts->SetCutCPA(0.995);

	    Posv0Daug->SetPID(AliPID::kProton, 999.9, 4);
	    Negv0Daug->SetPID(AliPID::kPion, 999.9, 4);
	    PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, 4);
	    NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, 4);

	    Posv0Daug->SetNClsTPC(80);
	    Negv0Daug->SetNClsTPC(80);
	    PosAntiv0Daug->SetNClsTPC(80);
	    NegAntiv0Daug->SetNClsTPC(80);

	    v0Cuts->SetCutDCADaugTov0Vtx(1.2);
	    Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);
	  } else if (CutVar == "43") {
	    TrackCuts->SetPtRange(0.6, 4.05);
	    AntiTrackCuts->SetPtRange(0.6, 4.05);

	    TrackCuts->SetEtaRange(-0.85, 0.85);
	    AntiTrackCuts->SetEtaRange(-0.85, 0.85);

	    TrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
	    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);

	    config->SetDeltaEtaMax(0.019);
	    config->SetDeltaPhiMax(0.019);

	    Posv0Daug->SetEtaRange(-0.83, 0.83);
	    Negv0Daug->SetEtaRange(-0.83, 0.83);
	    PosAntiv0Daug->SetEtaRange(-0.83, 0.83);
	    NegAntiv0Daug->SetEtaRange(-0.83, 0.83);
	  } else if (CutVar == "44") {
	    TrackCuts->SetEtaRange(-0.85, 0.85);
	    AntiTrackCuts->SetEtaRange(-0.85, 0.85);

	    TrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);
	    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 3.5);

	    config->SetDeltaEtaMax(0.019);
	    config->SetDeltaPhiMax(0.019);

	    v0Cuts->SetCutCPA(0.995);
	    Antiv0Cuts->SetCutCPA(0.995);

	    Posv0Daug->SetPID(AliPID::kProton, 999.9, 4);
	    Negv0Daug->SetPID(AliPID::kPion, 999.9, 4);
	    PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, 4);
	    NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, 4);

	    Posv0Daug->SetNClsTPC(80);
	    Negv0Daug->SetNClsTPC(80);
	    PosAntiv0Daug->SetNClsTPC(80);
	    NegAntiv0Daug->SetNClsTPC(80);

	    Posv0Daug->SetEtaRange(-0.83, 0.83);
	    Negv0Daug->SetEtaRange(-0.83, 0.83);
	    PosAntiv0Daug->SetEtaRange(-0.83, 0.83);
	    NegAntiv0Daug->SetEtaRange(-0.83, 0.83);

	    v0Cuts->SetCutDCADaugToPrimVtx(0.06);
	    Antiv0Cuts->SetCutDCADaugToPrimVtx(0.06);
	  }
	}
	else if(DimiSyst){
		const int NumObservables = 11;
		//last digit is the nsigma, first digit is bool smallest nsigma
		const double ptmin_p[] = {0.5,0.4,0.6};//0
		const double eta_p[] = {0.8,0.77,0.85};//1
		const double nsig_p[] = {3,2.5,3.5};//2
		const double ncls_p[] = {80,70,90};//3
		//const double pt_V0[] = {0.3,0.24,0.36};//4
		///apparantly, not included before.
		const double pt_V0[] = {0.3,0.3,0.3};//4
		const double cpa_v0[] = {0.999,0.9988,0.9992};//5
		const int V0d_var[] = {05,04};//6
		const double ncls_V0d[] = {70,80};//7
		const double eta_V0[] = {0.8,0.77,0.83};//8
		const double dca_ppi[] = {1.5,1.2};//9
		const double dca_V0d_PV[] = {0.05,0.06};//10

		const int SEED = iCutVar;
		TRandom3 rangen(SEED);
		//const int NumVar = 48;
		//vary only part of observables. For prob 4/11=0.36, we will end up with 4 simultanious variations on average
		const double ProbToVary = 4./11.;

		//for(int iVar=0; iVar<NumVar; iVar++){
			for(int iObs=0; iObs<NumObservables; iObs++){
				int iv=0;
				//non-default variation
				if(iCutVar!=100&&rangen.Uniform()<ProbToVary){
					if(iObs==6||iObs==7||iObs==9||iObs==10) iv = rangen.Integer(1)+1;
					else iv = rangen.Integer(2)+1;
				}
				switch (iObs) {
					case 0:
						TrackCuts->SetPtRange(ptmin_p[iv], 4.05);
						AntiTrackCuts->SetPtRange(ptmin_p[iv], 4.05);
						//printf("ptmin_p[%u]=%f\n",iv,ptmin_p[iv]);
						break;
					case 1:
						TrackCuts->SetEtaRange(-eta_p[iv],eta_p[iv]);
						AntiTrackCuts->SetEtaRange(-eta_p[iv], eta_p[iv]);
						//printf("eta_p[%u]=%f\n",iv,eta_p[iv]);
						break;
					case 2:
						TrackCuts->SetPID(AliPID::kProton, 0.75, nsig_p[iv]);
						AntiTrackCuts->SetPID(AliPID::kProton, 0.75, nsig_p[iv]);
						//printf("nsig_p[%u]=%f\n",iv,nsig_p[iv]);
						break;
					case 3:
						TrackCuts->SetNClsTPC(ncls_p[iv]);
						AntiTrackCuts->SetNClsTPC(ncls_p[iv]);
						//printf("ncls_p[%u]=%f\n",iv,ncls_p[iv]);
						break;
					case 4:
						v0Cuts->SetPtRange(pt_V0[iv], 999.9);
						Antiv0Cuts->SetPtRange(pt_V0[iv], 999.9);
						//printf("pt_V0[%u]=%f\n",iv,pt_V0[iv]);
						break;
					case 5:
						v0Cuts->SetCutCPA(cpa_v0[iv]);
						Antiv0Cuts->SetCutCPA(cpa_v0[iv]);
						//printf("cpa_v0[%u]=%f\n",iv,cpa_v0[iv]);
						break;
					case 6:
						Posv0Daug->SetCutSmallestSig(V0d_var[iv]/10);
						Negv0Daug->SetCutSmallestSig(V0d_var[iv]/10);
						PosAntiv0Daug->SetCutSmallestSig(V0d_var[iv]/10);
						NegAntiv0Daug->SetCutSmallestSig(V0d_var[iv]/10);

						Posv0Daug->SetPID(AliPID::kProton, 999.9, V0d_var[iv]%10);
						Negv0Daug->SetPID(AliPID::kPion, 999.9, V0d_var[iv]%10);
						PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, V0d_var[iv]%10);
						NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, V0d_var[iv]%10);

						//printf("V0d_var[%u]=%i\n",iv,V0d_var[iv]);
						break;
					case 7:
						Posv0Daug->SetNClsTPC(ncls_V0d[iv]);
						Negv0Daug->SetNClsTPC(ncls_V0d[iv]);
						PosAntiv0Daug->SetNClsTPC(ncls_V0d[iv]);
						NegAntiv0Daug->SetNClsTPC(ncls_V0d[iv]);
						//printf("ncls_V0d[%u]=%f\n",iv,ncls_V0d[iv]);
						break;
					case 8:
						Posv0Daug->SetEtaRange(-eta_V0[iv], eta_V0[iv]);
						Negv0Daug->SetEtaRange(-eta_V0[iv], eta_V0[iv]);
						PosAntiv0Daug->SetEtaRange(-eta_V0[iv], eta_V0[iv]);
						NegAntiv0Daug->SetEtaRange(-eta_V0[iv], eta_V0[iv]);
						//printf("eta_V0[%u]=%f\n",iv,eta_V0[iv]);
						break;
					case 9:
						v0Cuts->SetCutDCADaugTov0Vtx(dca_ppi[iv]);
						Antiv0Cuts->SetCutDCADaugTov0Vtx(dca_ppi[iv]);
						//printf("dca_ppi[%u]=%f\n",iv,dca_ppi[iv]);
						break;
					case 10:
						v0Cuts->SetCutDCADaugToPrimVtx(dca_V0d_PV[iv]);
						Antiv0Cuts->SetCutDCADaugToPrimVtx(dca_V0d_PV[iv]);
						//printf("dca_V0d_PV[%u]=%f\n",iv,dca_V0d_PV[iv]);
						break;
					default:
						break;
				}//switch
			}//iObs
		//}//iVar
	}//DimiSyst
	else if(DimiSpecial){

		//printf("iCPA=%i; iV0d=%i\n",iCPA,iV0d);

		if(iCPA==0)			{v0Cuts->SetCutCPA(0.99);		Antiv0Cuts->SetCutCPA(0.99);}
		else if(iCPA==1){v0Cuts->SetCutCPA(0.995);	Antiv0Cuts->SetCutCPA(0.995);}
		else if(iCPA==2){v0Cuts->SetCutCPA(0.999);	Antiv0Cuts->SetCutCPA(0.999);}
		else if(iCPA==3){v0Cuts->SetCutCPA(0.9995);	Antiv0Cuts->SetCutCPA(0.9995);}
		else if(iCPA==4){v0Cuts->SetCutCPA(0.9999);	Antiv0Cuts->SetCutCPA(0.9999);}

		//the 4 and 5 were used to do super high purity, without biasing the fraction of primaries
		if(iV0d==0){
			Posv0Daug->SetCutSmallestSig(false);
			Negv0Daug->SetCutSmallestSig(false);
			PosAntiv0Daug->SetCutSmallestSig(false);
			NegAntiv0Daug->SetCutSmallestSig(false);

			Posv0Daug->SetPID(AliPID::kProton, 999.9, 5);
			Negv0Daug->SetPID(AliPID::kPion, 999.9, 5);
			PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, 5);
			NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, 5);
		}
		else if(iV0d==1){
			Posv0Daug->SetCutSmallestSig(false);
			Negv0Daug->SetCutSmallestSig(false);
			PosAntiv0Daug->SetCutSmallestSig(false);
			NegAntiv0Daug->SetCutSmallestSig(false);

			Posv0Daug->SetPID(AliPID::kProton, 999.9, 3);
			Negv0Daug->SetPID(AliPID::kPion, 999.9, 3);
			PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, 3);
			NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, 3);
		}
		else if(iV0d==2){
			Posv0Daug->SetCutSmallestSig(true);
			Negv0Daug->SetCutSmallestSig(true);
			PosAntiv0Daug->SetCutSmallestSig(true);
			NegAntiv0Daug->SetCutSmallestSig(true);

			Posv0Daug->SetPID(AliPID::kProton, 999.9, 5);
			Negv0Daug->SetPID(AliPID::kPion, 999.9, 5);
			PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, 5);
			NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, 5);
		}
		else if(iV0d==3){
			Posv0Daug->SetCutSmallestSig(true);
			Negv0Daug->SetCutSmallestSig(true);
			PosAntiv0Daug->SetCutSmallestSig(true);
			NegAntiv0Daug->SetCutSmallestSig(true);

			Posv0Daug->SetPID(AliPID::kProton, 999.9, 3);
			Negv0Daug->SetPID(AliPID::kPion, 999.9, 3);
			PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, 3);
			NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, 3);
		}
		else if(iV0d==4){
			Posv0Daug->SetCutSmallestSig(true);
			Negv0Daug->SetCutSmallestSig(true);
			PosAntiv0Daug->SetCutSmallestSig(true);
			NegAntiv0Daug->SetCutSmallestSig(true);

			Posv0Daug->SetPID(AliPID::kProton, 999.9, 3);
			Negv0Daug->SetPID(AliPID::kPion, 999.9, 3);
			PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, 3);
			NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, 3);

			v0Cuts->SetCutDCADaugTov0Vtx(0.1);
			Antiv0Cuts->SetCutDCADaugTov0Vtx(0.1);

		}
		else if(iV0d==5){
			Posv0Daug->SetCutSmallestSig(true);
			Negv0Daug->SetCutSmallestSig(true);
			PosAntiv0Daug->SetCutSmallestSig(true);
			NegAntiv0Daug->SetCutSmallestSig(true);

			Posv0Daug->SetPID(AliPID::kProton, 999.9, 3);
			Negv0Daug->SetPID(AliPID::kPion, 999.9, 3);
			PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, 3);
			NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, 3);

			v0Cuts->SetCutDCADaugTov0Vtx(0.05);
			Antiv0Cuts->SetCutDCADaugTov0Vtx(0.05);
		}
	}

  config->SetMultiplicityEstimator(AliFemtoDreamEvent::kRef08);

  std::vector<int> MultBins;
  MultBins.push_back(0);
  MultBins.push_back(4);
  MultBins.push_back(8);
  MultBins.push_back(12);
  MultBins.push_back(16);
  MultBins.push_back(20);
  MultBins.push_back(24);
  MultBins.push_back(28);
  MultBins.push_back(32);
  MultBins.push_back(36);
  MultBins.push_back(40);
  MultBins.push_back(44);
  MultBins.push_back(48);
  MultBins.push_back(52);
  MultBins.push_back(56);
  MultBins.push_back(60);
  MultBins.push_back(64);
  MultBins.push_back(68);
  MultBins.push_back(72);
  MultBins.push_back(76);
  MultBins.push_back(80);
  MultBins.push_back(84);
  MultBins.push_back(88);
  MultBins.push_back(92);
  MultBins.push_back(96);
  MultBins.push_back(100);

  config->SetMultBins(MultBins);

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

  config->SetZBins(ZVtxBins);

  config->SetMultBinning(true);
  config->SetmTBinning(true);

  config->SetdPhidEtaPlotsSmallK(false);
  config->SetdPhidEtaPlots(false);
  config->SetPhiEtaBinnign(false);

  if (fullBlastQA) {
    config->SetkTBinning(true);
    config->SetPtQA(true);
    config->SetMassQA(true);
  }

  if (!fullBlastQA) {
    config->SetMinimalBookingME(true);
    config->SetMinimalBookingSample(true);
  }

  if (isMC) {
    config->SetMomentumResolution(true);//kstar true vs. kstar reco
  } else {
    std::cout
        << "You are trying to request the Momentum Resolution without MC Info; fix it wont work! \n";
  }

  AliAnalysisTaskNanoBBar* task = new AliAnalysisTaskNanoBBar("femtoDimi",isMC);

  if (!fullBlastQA) {
    task->SetRunTaskLightWeight(true);
  }
  if(triggerData=="kINT7"){
	  task->SelectCollisionCandidates(AliVEvent::kINT7);
  }else if(triggerData=="kHM"){
	  task->SelectCollisionCandidates(AliVEvent::kHighMultV0);
  }
  task->SetEventCuts(evtCuts);
  task->SetProtonCuts(TrackCuts);
  task->SetAntiProtonCuts(AntiTrackCuts);
  task->Setv0Cuts(v0Cuts);
  task->SetAntiv0Cuts(Antiv0Cuts);
  task->SetXiCuts(CascadeCuts);
  task->SetAntiXiCuts(AntiCascadeCuts);
  task->SetCorrelationConfig(config);
  task->SetUseDumpster(false);
  mgr->AddTask(task);

  TString addon = "";
  if (triggerData == "kINT7") {
    addon += "MBDimi";
  } else if (triggerData == "kHM") {
    addon += "HMDimi";
  }

  TString file = AliAnalysisManager::GetCommonFileName();

  mgr->ConnectInput(task, 0, cinput);

  TString QAName = Form("%sQA%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputQA = mgr->CreateContainer(
		  QAName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), QAName.Data()));
  mgr->ConnectOutput(task, 1, coutputQA);

  TString EvtCutsName = Form("%sEvtCuts%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *coutputEvtCuts = mgr->CreateContainer(
      EvtCutsName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), EvtCutsName.Data()));
  mgr->ConnectOutput(task, 2, coutputEvtCuts);

  TString TrackCutsName = Form("%sTrackCuts%s", addon.Data(), suffix.Data());
  AliAnalysisDataContainer *couputTrkCuts = mgr->CreateContainer(
      TrackCutsName.Data(), TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), TrackCutsName.Data()));
  mgr->ConnectOutput(task, 3, couputTrkCuts);

  TString AntiTrackCutsName = Form("%sAntiTrackCuts%s", addon.Data(),
                                   suffix.Data());
  AliAnalysisDataContainer *coutputAntiTrkCuts = mgr->CreateContainer(
      AntiTrackCutsName.Data(), TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), AntiTrackCutsName.Data()));
  mgr->ConnectOutput(task, 4, coutputAntiTrkCuts);

  AliAnalysisDataContainer *coutputv0Cuts;
  TString v0CutsName = Form("%sv0Cuts%s", addon.Data(), suffix.Data());
  coutputv0Cuts = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      v0CutsName.Data(),
      TList::Class(), AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), v0CutsName.Data()));
  mgr->ConnectOutput(task, 5, coutputv0Cuts);

  AliAnalysisDataContainer *coutputAntiv0Cuts;
  TString Antiv0CutsName = Form("%sAntiv0Cuts%s", addon.Data(), suffix.Data());
  coutputAntiv0Cuts = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      Antiv0CutsName.Data(),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), Antiv0CutsName.Data()));
  mgr->ConnectOutput(task, 6, coutputAntiv0Cuts);

  AliAnalysisDataContainer *coutputXiCuts;
  TString XiCutsName = Form("%sXiCuts%s", addon.Data(), suffix.Data());
  coutputXiCuts = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      XiCutsName.Data(),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), XiCutsName.Data()));
  mgr->ConnectOutput(task, 7, coutputXiCuts);

  AliAnalysisDataContainer *coutputAntiXiCuts;
  TString AntiXiCutsName = Form("%sAntiXiCuts%s", addon.Data(), suffix.Data());
  coutputAntiXiCuts = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      AntiXiCutsName.Data(),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), AntiXiCutsName.Data()));
  mgr->ConnectOutput(task, 8, coutputAntiXiCuts);

  AliAnalysisDataContainer *coutputResults;
  TString ResultsName = Form("%sResults%s", addon.Data(), suffix.Data());
  coutputResults = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      ResultsName.Data(),
      TList::Class(), AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), ResultsName.Data()));
  mgr->ConnectOutput(task, 9, coutputResults);

  AliAnalysisDataContainer *coutputResultsQA;
  TString ResultsQAName = Form("%sResultsQA%s", addon.Data(), suffix.Data());
  coutputResultsQA = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      ResultsQAName.Data(),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), ResultsQAName.Data()));
  mgr->ConnectOutput(task, 10, coutputResultsQA);

  AliAnalysisDataContainer *coutputResultsSample;
  TString ResultsSampleName = Form("%sResultsSample%s", addon.Data(),
                                   suffix.Data());
  coutputResultsSample = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      ResultsSampleName.Data(),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), ResultsSampleName.Data()));
  mgr->ConnectOutput(task, 11, coutputResultsSample);

  AliAnalysisDataContainer *coutputResultsSampleQA;
  TString ResultsSampleQAName = Form("%sResultsSampleQA%s", addon.Data(),
                                     suffix.Data());
  coutputResultsSampleQA = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      ResultsSampleQAName.Data(),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), ResultsSampleQAName.Data()));
  mgr->ConnectOutput(task, 12, coutputResultsSampleQA);

  AliAnalysisDataContainer *coutputDumpster;
  TString DumpsterName = Form("%sDumpster%s", addon.Data(), suffix.Data());
  coutputDumpster = mgr->CreateContainer(
      //@suppress("Invalid arguments") it works ffs
      DumpsterName.Data(),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), DumpsterName.Data()));
  mgr->ConnectOutput(task, 13, coutputDumpster);

   if (isMC) {
    AliAnalysisDataContainer *coutputTrkCutsMC;
    TString TrkCutsMCName = Form("%sTrkCutsMC%s",addon.Data(),suffix.Data());
    coutputTrkCutsMC = mgr->CreateContainer(
        //@suppress("Invalid arguments") it works ffs
        TrkCutsMCName.Data(),
        TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), TrkCutsMCName.Data()));
    mgr->ConnectOutput(task, 14, coutputTrkCutsMC);

    AliAnalysisDataContainer *coutputAntiTrkCutsMC;
    TString AntiTrkCutsMCName = Form("%sAntiTrkCutsMC%s",addon.Data(),suffix.Data());
    coutputAntiTrkCutsMC = mgr->CreateContainer(
        //@suppress("Invalid arguments") it works ffs
        AntiTrkCutsMCName.Data(),
        TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), AntiTrkCutsMCName.Data()));
    mgr->ConnectOutput(task, 15, coutputAntiTrkCutsMC);

    AliAnalysisDataContainer *coutputv0CutsMC;
    TString v0CutsMCName = Form("%sv0CutsMC%s",addon.Data(),suffix.Data());
    coutputv0CutsMC = mgr->CreateContainer(
        //@suppress("Invalid arguments") it works ffs
        v0CutsMCName.Data(),
        TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), v0CutsMCName.Data()));
    mgr->ConnectOutput(task, 16, coutputv0CutsMC);

    AliAnalysisDataContainer *coutputAntiv0CutsMC;
    TString Antiv0CutsMCName = Form("%sAntiv0CutsMC%s",addon.Data(),suffix.Data());
    coutputAntiv0CutsMC = mgr->CreateContainer(
        //@suppress("Invalid arguments") it works ffs
        Antiv0CutsMCName.Data(),
        TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), Antiv0CutsMCName.Data()));
    mgr->ConnectOutput(task, 17, coutputAntiv0CutsMC);

    AliAnalysisDataContainer *coutputXiCutsMC;
    TString XiCutsMCName = Form("%sXiCutsMC%s",addon.Data(),suffix.Data());
    coutputXiCutsMC = mgr->CreateContainer(
        //@suppress("Invalid arguments") it works ffs
        XiCutsMCName.Data(),
        TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), XiCutsMCName.Data()));
    mgr->ConnectOutput(task, 18, coutputXiCutsMC);

    AliAnalysisDataContainer *coutputAntiXiCutsMC;
    TString AntiXiCutsMCName = Form("%sAntiXiCutsMC%s",addon.Data(),suffix.Data());
    coutputAntiXiCutsMC = mgr->CreateContainer(
        //@suppress("Invalid arguments") it works ffs
        AntiXiCutsMCName.Data(),
        TList::Class(),
        AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), AntiXiCutsMCName.Data()));
    mgr->ConnectOutput(task, 19, coutputAntiXiCutsMC);

   }

  return task;
}
