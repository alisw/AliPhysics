
/*********************************************************************
 *                                                                   *
 * ConfigFemtoAnalysis.C - configuration macro for the femtoscopic   *
 * analysis, meant as a QA process for two-particle effects          *
 *                                                                   *
 * Author: Adam Kisiel (Adam.Kisiel@cern.ch)                         *
 *********************************************************************
 *    K+K- in PbPb@2.76TeV                                           *
 * Update: Konstantin.Mikhaylov@cern.ch                              *
 *         KpmHBT16  0010  03-Mar-2016                               *
 *         TOF PID from 0.45 GeV/c                                   *
 *         SetNsigmaTPC400_450(1.0) and other set to standard        *
 *         no cut on mom and no on gamma                             *
 *         no cut on phi* eta (set to 0)                             *
 *         10 MeV/q-bin in CF                                        *
 *********************************************************************/

#if !defined(__CINT__) || defined(__MAKECINT_)
#include "AliFemtoManager.h"
#include "AliFemtoEventReaderESDChain.h"
#include "AliFemtoEventReaderESDChainKine.h"
#include "AliFemtoEventReaderAODChain.h"
#include "AliFemtoSimpleAnalysis.h"
#include "AliFemtoBasicEventCut.h"
#include "AliFemtoESDTrackCut.h"
//#include "AliFemtoKKTrackCut.h"
#include "AliFemtoKpm45TrackCut.h"
#include "AliFemtoCorrFctn.h"
#include "AliFemtoCutMonitorParticleYPt.h"
#include "AliFemtoCutMonitorParticleVertPos.h"
#include "AliFemtoCutMonitorParticleMomRes.h"
#include "AliFemtoCutMonitorParticlePID.h"
#include "AliFemtoCutMonitorEventMult.h"
#include "AliFemtoCutMonitorEventVertex.h"
#include "AliFemtoShareQualityTPCEntranceSepPairCut.h"
#include "AliFemtoPairCutAntiGammaAlpha.h"
#include "AliFemtoPairCutRadialDistance.h"
#include "AliFemtoPairCutRadialDistanceKK.h"
#include "AliFemtoQinvCorrFctn.h"
#include "AliFemtoShareQualityCorrFctn.h"
#include "AliFemtoTPCInnerCorrFctn.h"
#include "AliFemtoVertexMultAnalysis.h"
#include "AliFemtoCorrFctn3DSpherical.h"
#include "AliFemtoChi2CorrFctn.h"
#include "AliFemtoCorrFctnTPCNcls.h"
#include "AliFemtoBPLCMS3DCorrFctn.h"
#include "AliFemtoCorrFctn3DLCMSSym.h"
#include "AliFemtoModelBPLCMSCorrFctn.h"
#include "AliFemtoModelCorrFctn3DSpherical.h"
#include "AliFemtoModelGausLCMSFreezeOutGenerator.h"
#include "AliFemtoModelGausRinvFreezeOutGenerator.h"
#include "AliFemtoModelManager.h"
#include "AliFemtoModelWeightGeneratorBasic.h"
#include "AliFemtoModelWeightGeneratorLednicky.h"
#include "AliFemtoCorrFctnDirectYlm.h"
#include "AliFemtoModelCorrFctnDirectYlm.h"
#include "AliFemtoModelCorrFctnSource.h"
#include "AliFemtoModelCorrFctn.h"
#include "AliFemtoKTPairCut.h"
//... double ratio ...
#include "AliFemtoCorrFctnNonIdDR.h"
#endif

//________________________________________________________________________
AliFemtoManager* ConfigFemtoAnalysis() {

  double PionMass = 0.13956995;
  double KaonMass = 0.493677;
  const int cMu=1;
  const int cKt=2;

  //-------Single track cuts------------------------------------------------->
  double DCAxy=0.3;//cm // our standard is 0.20 cm; super narrow was 0.015cm
  double DCAz =0.3;//cm // our standard is 0.15 cm;
  //-------Single track cuts-------------------------------------------------<
  //=======Double track cuts=================================================>
  //Dhevan's : PhiStarDifferenceMinimum=0.06; EtaDifferenceMinimum=0.02;
  double PhiStarDifferenceMinimum=0.02; //[radian]
  double EtaDifferenceMinimum=0.04; //[radian]
  //=======Double track cuts=================================================<

  // Switches for QA analyses
  int runmults[cMu] = {1};//, 0, 0, 0};
  int multbins[cMu+1] = {0, 900};//, 300, 500, 900};

  int runch[2] = {1, 0};//K+-
  const char *chrgs[2] = { "Kp", "Km"};
  
  int runktdep = 1;
  //double ktrng[cKt+1] = {0.2, 0.5, 1.0};//orig  
  double ktrng[cKt+1] = {0.2, 0.5, 1.5};
  //double ktrng[cKt+1] = {0.2, 1.0};
  
  

  int run3d = 0;
  int runshlcms = 0;

  int runtype = 2; // Types 0 - global, 1 - ITS only, 2 - TPC Inner
  int isrealdata = 1;

  //   AliFemtoEventReaderESDChainKine* Reader=new AliFemtoEventReaderESDChainKine();
  //   Reader->SetConstrained(true);
  //   Reader->SetUseTPCOnly(false);

  double shqmax;
  double shqmaxSH;
  int nbinssh = 200;// 10 MeV/q-bin //orig 100, 20 MeV per bin

  if (runshlcms) shqmaxSH = 0.25;
  shqmax = 2.0;//K+-
  
  

  AliFemtoEventReaderAODChain *Reader = new AliFemtoEventReaderAODChain();
    Reader->SetFilterBit(7);
/////    Reader->SetCentralityPreSelection(0.0001, 900);
    Reader->SetUseMultiplicity(AliFemtoEventReaderAOD::kCentrality);

    Reader->SetReadMC(kTRUE);
    Reader->SetDCAglobalTrack(kTRUE);//proverit'!!!
    //
    Reader->SetKaonAnalysis(kTRUE);// new set to choose Kaon by PDG code (7Nov2016)


/*
AliFemtoEventReaderAOD *Reader = new AliFemtoEventReaderAODMultSelection();
//AliFemtoEventReaderAOD *Reader = new AliFemtoEventReaderAOD();
    Reader->SetFilterBit(7);
    Reader->SetEPVZERO(kTRUE);
    Reader->SetUseMultiplicity(AliFemtoEventReaderAOD::kCentrality);
/////    Reader->SetCentralityFlattening(kFALSE);
    Reader->SetReadV0(0);
    // rdr->SetPrimaryVertexCorrectionTPCPoints(kTRUE);
    Reader->SetDCAglobalTrack(kTRUE);
*/


    /*AliFemtoModelGausLCMSFreezeOutGenerator *tFreeze = new AliFemtoModelGausLCMSFreezeOutGenerator();
    tFreeze->SetSizeOut(3.0*TMath::Sqrt(2.0));         
    tFreeze->SetSizeSide(3.0*TMath::Sqrt(2.0));
    tFreeze->SetSizeLong(3.0*TMath::Sqrt(2.0));*/

    //Generate freeze-out coordinates as a 3D gaussian sphere in PRF
    AliFemtoModelGausRinvFreezeOutGenerator *tFreeze = new AliFemtoModelGausRinvFreezeOutGenerator();
    tFreeze->SetSizeInv(5.0*TMath::Sqrt(2.0));//r_0=3fm it should be time to sqrt(2) !!!!KM
    
    
    AliFemtoModelWeightGeneratorLednicky *tWeight = new AliFemtoModelWeightGeneratorLednicky();
    tWeight->SetPairType(AliFemtoModelWeightGenerator::KaonPlusKaonPlus());
    //tWeight->SetCoulOff();
    tWeight->SetCoulOn();
    tWeight->SetQuantumOn();
    tWeight->SetStrongOff();
//    tWeight->SetStrongOn();
    tWeight->Set3BodyOff();

/////////////    tWeight->SetKpKmModelType(14,1);//(Martin-f0, Achasov2-a0,0->phi off;1->phi on)
 

  AliFemtoModelManager *tModelManager = new AliFemtoModelManager();
  tModelManager->AcceptFreezeOutGenerator(tFreeze);
  tModelManager->AcceptWeightGenerator(tWeight);
  tModelManager->CreateCopyHiddenInfo(kTRUE);

  AliFemtoManager* Manager=new AliFemtoManager();
  Manager->SetEventReader(Reader);

  AliFemtoVertexMultAnalysis    *anetaphitpc[20];
  AliFemtoBasicEventCut         *mecetaphitpc[20];
  AliFemtoCutMonitorEventMult   *cutPassEvMetaphitpc[20];
  AliFemtoCutMonitorEventMult   *cutFailEvMetaphitpc[20];
  AliFemtoCutMonitorEventVertex *cutPassEvVetaphitpc[20];
  AliFemtoCutMonitorEventVertex *cutFailEvVetaphitpc[20];
  AliFemtoKKTrackCut           *dtc1etaphitpc[20];
  AliFemtoKKTrackCut           *dtc2etaphitpc[20];
  //AliFemtoKpm45TrackCut           *dtc1etaphitpc[20];
  //AliFemtoKpm45TrackCut           *dtc2etaphitpc[20];
  //  AliFemtoESDTrackCut           *dtc1etaphitpc[20];
  //  AliFemtoESDTrackCut           *dtc2etaphitpc[20];
  AliFemtoCutMonitorParticleYPt *cutPass1YPtetaphitpc[20];
  AliFemtoCutMonitorParticleYPt *cutFail1YPtetaphitpc[20];
  AliFemtoCutMonitorParticlePID *cutPass1PIDetaphitpc[20];
  AliFemtoCutMonitorParticlePID *cutFail1PIDetaphitpc[20];
  AliFemtoCutMonitorParticleYPt *cutPass2YPtetaphitpc[20];
  AliFemtoCutMonitorParticleYPt *cutFail2YPtetaphitpc[20];
  AliFemtoCutMonitorParticlePID *cutPass2PIDetaphitpc[20];
  AliFemtoCutMonitorParticlePID *cutFail2PIDetaphitpc[20];
  ////AliFemtoPairCutAntiGammaAlpha      *sqpcetaphitpc[20];//K+-
  //    AliFemtoShareQualityTPCEntranceSepPairCut      *sqpcetaphitpc[20];
  //AliFemtoPairCutRadialDistance      *sqpcetaphitpc[20];//AliFemto dphi* cut
  //~~~
  AliFemtoPairCutRadialDistanceKK      *sqpcetaphitpc[20];//Dhevan's dphi* cut
  //
  //AliFemtoCorrFctnDirectYlm     *cylmetaphitpc[20];
  AliFemtoCorrFctnDEtaDPhi      *cdedpetaphi[20*10];//20->20*10 due to kT
  //AliFemtoChi2CorrFctn          *cchiqinvetaphitpc[20];
  AliFemtoKTPairCut             *ktpcuts[20*8];
  //AliFemtoCorrFctnDirectYlm     *cylmkttpc[20*8];
  AliFemtoQinvCorrFctn          *cqinvkttpc[20*8];
  //AliFemtoCorrFctn3DLCMSSym     *cq3dlcmskttpc[20*8];
  AliFemtoCorrFctnTPCNcls       *cqinvnclstpc[20];
  AliFemtoShareQualityCorrFctn  *cqinvsqtpc[20*10];
  AliFemtoChi2CorrFctn          *cqinvchi2tpc[20];
  AliFemtoTPCInnerCorrFctn      *cqinvinnertpc[20*10];
  AliFemtoCorrFctnGammaMonitorAlpha  *cgamma[20*10];
  //... double ratio ...
  // AliFemtoCorrFctnNonIdDR       *cfdourat[20*10];
  AliFemtoModelCorrFctn   *cqinvkttpcmodel[20*8];
 // AliFemtoCorrFctnGammaMonitorAlpha  *cgamma[20*10];


  // *** Third QA task - HBT analysis with all pair cuts off, TPC only ***
  // *** Begin K+K- analysis ***
  int aniter = 0;

  for (int imult=0; imult<cMu/*4*/; imult++) {
    if (runmults[imult]) {
      for (int ichg=0; ichg<1/*K+-*/; ichg++) {//one loop
	if (runch[ichg]) {
	  aniter = ichg*5+imult;

	  anetaphitpc[aniter] = new AliFemtoVertexMultAnalysis(4, -8.0, 8.0, 5, multbins[imult], multbins[imult+1]);
	  anetaphitpc[aniter]->SetNumEventsToMix(3);
	  anetaphitpc[aniter]->SetMinSizePartCollection(1);

	  mecetaphitpc[aniter] = new AliFemtoBasicEventCut();
	  mecetaphitpc[aniter]->SetEventMult(0,100000);
	  mecetaphitpc[aniter]->SetVertZPos(-8.0,8.0);
	  //    mecetaphitpc->SetAcceptBadVertex(kTRUE);
	  
	  //Multiplicity---->
	  cutPassEvMetaphitpc[aniter] = new AliFemtoCutMonitorEventMult(Form("cutPass%stpcM%i", chrgs[ichg], imult));
	  cutFailEvMetaphitpc[aniter] = new AliFemtoCutMonitorEventMult(Form("cutFail%stpcM%i", chrgs[ichg], imult));
	  mecetaphitpc[aniter]->AddCutMonitor(cutPassEvMetaphitpc[aniter], cutFailEvMetaphitpc[aniter]);
	  //Multiplicity----<
	  
	  //Vertex---------->
	  cutPassEvVetaphitpc[aniter] = new AliFemtoCutMonitorEventVertex(Form("cutPass%stpcM%i", chrgs[ichg], imult));
	  cutFailEvVetaphitpc[aniter] = new AliFemtoCutMonitorEventVertex(Form("cutFail%stpcM%i", chrgs[ichg], imult));
	  mecetaphitpc[aniter]->AddCutMonitor(cutPassEvVetaphitpc[aniter], cutFailEvVetaphitpc[aniter]);
	  //Vertex----------<
	  
	  dtc1etaphitpc[aniter] = new AliFemtoKKTrackCut();
	  dtc2etaphitpc[aniter] = new AliFemtoKKTrackCut();//K+-
	  //dtc1etaphitpc[aniter] = new AliFemtoKpm45TrackCut();
	  //dtc2etaphitpc[aniter] = new AliFemtoKpm45TrackCut();//K+-
	  //--- K+- --->
	  dtc1etaphitpc[aniter]->SetCharge(1.0);
	  dtc2etaphitpc[aniter]->SetCharge(1.0);
	  //--- K+- ---<
	    
	  dtc1etaphitpc[aniter]->SetPt(0.14,1.5); //0.14,1.5); //--- K+-
	  dtc1etaphitpc[aniter]->SetEta(-0.8,0.8);//--- K+-
	  dtc2etaphitpc[aniter]->SetPt(0.14,1.5); //0.14,1.5);  //--- K+-
	  dtc2etaphitpc[aniter]->SetEta(-0.8,0.8); //--- K+-
	  
	  //--- Choose Kaon as Most Probable (switch on all cuts: TPC, TOF)---
	  dtc1etaphitpc[aniter]->SetMass(KaonMass);
	  dtc1etaphitpc[aniter]->SetMostProbableKaon();
	  dtc2etaphitpc[aniter]->SetMass(KaonMass);
	  dtc2etaphitpc[aniter]->SetMostProbableKaon();
//------------------- November 2013 -----------------------------------< 
         // new cuts to remove electron (do not take into analysis if 400<p<500) 
	 dtc1etaphitpc[aniter]->SetNsigmaTPCle250(2.0);
	 dtc1etaphitpc[aniter]->SetNsigmaTPC250_400(2.0);
	 dtc1etaphitpc[aniter]->SetNsigmaTPC400_450(1.0);//cut on e+e- orig(2.0);
	 dtc1etaphitpc[aniter]->SetNsigmaTPC450_500(2.0);//cut on e+e- orig(2.0);
	 dtc1etaphitpc[aniter]->SetNsigmaTPCge500(3.0);  
	 // new cuts are stronger, better separation of pion in TOF 
	 // when momentum is greater then 800 MeV/c
	 dtc1etaphitpc[aniter]->SetNsigmaTOF500_800(2.0);
	 dtc1etaphitpc[aniter]->SetNsigmaTOF800_1000(1.5);
	 dtc1etaphitpc[aniter]->SetNsigmaTOFge1000(1.0);
	 //K+ contour (default p0=9999, p1=0, p2=0):  
	 //dtc1etaphitpc[aniter]->SetKaonTCPdEdxContourP0( 196.25);//december 2014
	 //dtc1etaphitpc[aniter]->SetKaonTCPdEdxContourP1(-420.00);//december 2014
	 //dtc1etaphitpc[aniter]->SetKaonTCPdEdxContourP2( 300.00);//december 2014

	 dtc2etaphitpc[aniter]->SetNsigmaTPCle250(2.0);
	 dtc2etaphitpc[aniter]->SetNsigmaTPC250_400(2.0);
	 dtc2etaphitpc[aniter]->SetNsigmaTPC400_450(1.0);//cut on e+e- orig(2.0);
	 dtc2etaphitpc[aniter]->SetNsigmaTPC450_500(2.0);//cut on e+e- orig(2.0);
	 dtc2etaphitpc[aniter]->SetNsigmaTPCge500(3.0);    
	 // new cuts are stronger, better separation of pion in TOF 
	 // when momentum is greater then 800 MeV/c
	 dtc2etaphitpc[aniter]->SetNsigmaTOF500_800(2.0);
	 dtc2etaphitpc[aniter]->SetNsigmaTOF800_1000(1.5);
	 dtc2etaphitpc[aniter]->SetNsigmaTOFge1000(1.0);
	 //K+ contour (default p0=9999, p1=0, p2=0):  
	 //dtc2etaphitpc[aniter]->SetKaonTCPdEdxContourP0( 196.25);//december 2014
	 //dtc2etaphitpc[aniter]->SetKaonTCPdEdxContourP1(-420.00);//december 2014
	 //dtc2etaphitpc[aniter]->SetKaonTCPdEdxContourP2( 300.00);//december 2014
	 //------------------- November 2013 ----------------------------------->
	  
	  // *** Track quality cuts ***
	

  cutPassEvMetaphitpc[aniter] = new AliFemtoCutMonitorEventMult(Form("cutPass%stpcM%i", chrgs[ichg], imult));
	  cutFailEvMetaphitpc[aniter] = new AliFemtoCutMonitorEventMult(Form("cutFail%stpcM%i", chrgs[ichg], imult));
	  mecetaphitpc[aniter]->AddCutMonitor(cutPassEvMetaphitpc[aniter], cutFailEvMetaphitpc[aniter]);
	  
	  cutPassEvVetaphitpc[aniter] = new AliFemtoCutMonitorEventVertex(Form("cutPass%stpcM%i", chrgs[ichg], imult));
	  cutFailEvVetaphitpc[aniter] = new AliFemtoCutMonitorEventVertex(Form("cutFail%stpcM%i", chrgs[ichg], imult));
	  mecetaphitpc[aniter]->AddCutMonitor(cutPassEvVetaphitpc[aniter], cutFailEvVetaphitpc[aniter]);
	  
	  dtc1etaphitpc[aniter] = new AliFemtoKKTrackCut();
//	  dtc1etaphitpc[aniter] = new AliFemtoESDTrackCut();
	  //     dtc1etaphitpc[aniter]->SetPidProbPion(0.0,1.001);
	  //     dtc1etaphitpc[aniter]->SetPidProbMuon(0.0,1.0);
	  //     dtc1etaphitpc[aniter]->SetPidProbKaon(0.0,1.0);
	  //     dtc1etaphitpc[aniter]->SetPidProbProton(0.0,1.0);
	  if (ichg == 0)
	    dtc1etaphitpc[aniter]->SetCharge(1.0);
	  else if (ichg == 1)
	    dtc1etaphitpc[aniter]->SetCharge(-1.0);
	    
	  dtc1etaphitpc[aniter]->SetPt(0.14,1.5);
	  //	  dtc1etaphitpc[aniter]->SetEta(-1.2,1.2);
	  dtc1etaphitpc[aniter]->SetEta(-0.8,0.8); //0.5
	  // 	//    dtc1etaphitpc[aniter]->SetEta(-0.5,0.5);
///	  dtc1etaphitpc[aniter]->SetMass(PionMass);
	  dtc1etaphitpc[aniter]->SetMass(KaonMass);
	  
	  
	  ////	  dtc1etaphitpc[aniter]->SetminTPCncls(80);
		  
///////   ----!!!!!!	   
	  dtc1etaphitpc[aniter]->SetMostProbableKaon();  //!!!!!!
	  //------------------- November 2013 -----------------------------------< 
	  //New class in AliFemo: PWGCF/FEMTOSCOPY/AliFemtoUser/AliFemtoKKTrackCut.cxx
	  dtc1etaphitpc[aniter]->SetNsigmaTPCle250(2.0);
	  dtc1etaphitpc[aniter]->SetNsigmaTPC250_400(2.0);
	  dtc1etaphitpc[aniter]->SetNsigmaTPC400_450(2.0);
	  dtc1etaphitpc[aniter]->SetNsigmaTPC450_500(2.0);
	  dtc1etaphitpc[aniter]->SetNsigmaTPCge500(3.0);
	  // new cuts are stronger, better separation of pion in TOF 
	  // when momentum is greater then 800 MeV/c
	  dtc1etaphitpc[aniter]->SetNsigmaTOF500_800(2.0);
	  dtc1etaphitpc[aniter]->SetNsigmaTOF800_1000(1.5);
	  dtc1etaphitpc[aniter]->SetNsigmaTOFge1000(1.0);
	  //------------------- November 2013 ----------------------------------->

////	  	  dtc1etaphitpc[aniter]->SetMostProbablePion();
	  // 	// Track quality cuts
	  if (runtype == 0) {
	    dtc1etaphitpc[aniter]->SetStatus(AliESDtrack::kTPCrefit|AliESDtrack::kITSrefit);
	    //	    dtc1etaphitpc[aniter]->SetStatus(AliESDtrack::kTPCrefit);
	    //    dtc1etaphitpc[aniter]->SetStatus(AliESDtrack::kITSrefit);
	    dtc1etaphitpc[aniter]->SetminTPCncls(80);
	    dtc1etaphitpc[aniter]->SetRemoveKinks(kTRUE);
	    dtc1etaphitpc[aniter]->SetLabel(kFALSE);
	    //    dtc1etaphitpc[aniter]->SetMaxITSChiNdof(6.0);
	    dtc1etaphitpc[aniter]->SetMaxTPCChiNdof(4.0);
	    dtc1etaphitpc[aniter]->SetMaxImpactXY(DCAxy);
	    //Poland: dtc1etaphitpc[aniter]->SetMaxImpactXYPtDep(0.0182, 0.0350, -1.01);
	    dtc1etaphitpc[aniter]->SetMaxImpactZ(DCAz);
	    //      dtc1etaphitpc[aniter]->SetMaxSigmaToVertex(6.0);
	  }
	  else if (runtype == 1) {
	    //      dtc1etaphitpc[aniter]->SetStatus(AliESDtrack::kTPCrefit|AliESDtrack::kITSrefit);
	    //    dtc1etaphitpc[aniter]->SetStatus(AliESDtrack::kTPCrefit);
	    //	    dtc1etaphitpc[aniter]->SetStatus(AliESDtrack::kITSrefit|AliESDtrack::kITSpureSA);
	    //      dtc1etaphitpc[aniter]->SetminTPCncls(70);
	    dtc1etaphitpc[aniter]->SetStatus(AliESDtrack::kITSrefit);
	    dtc1etaphitpc[aniter]->SetRemoveKinks(kTRUE);
	    dtc1etaphitpc[aniter]->SetLabel(kFALSE);
	    //    dtc1etaphitpc[aniter]->SetMaxITSChiNdof(6.0);
	    //      dtc1etaphitpc[aniter]->SetMaxTPCChiNdof(6.0);
	    dtc1etaphitpc[aniter]->SetMaxImpactXY(DCAxy);
	    dtc1etaphitpc[aniter]->SetMaxImpactZ(DCAz);
	    //      dtc1etaphitpc[aniter]->SetMaxSigmaToVertex(6.0);
	  }
	  else if (runtype == 2) {
	    //	    dtc1etaphitpc[aniter]->SetStatus(AliESDtrack::kTPCrefit|AliESDtrack::kITSrefit);
	    dtc1etaphitpc[aniter]->SetStatus(AliESDtrack::kTPCin);
	    //	    dtc1etaphitpc[aniter]->SetStatus(AliESDtrack::kTPCrefit);
	    //    dtc1etaphitpc[aniter]->SetStatus(AliESDtrack::kITSrefit);
	    dtc1etaphitpc[aniter]->SetminTPCncls(80); //was "0"
	    dtc1etaphitpc[aniter]->SetRemoveKinks(kTRUE);
	    dtc1etaphitpc[aniter]->SetLabel(kFALSE);
	    //    dtc1etaphitpc[aniter]->SetMaxITSChiNdof(6.0);
	    dtc1etaphitpc[aniter]->SetMaxTPCChiNdof(4.0);
	    dtc1etaphitpc[aniter]->SetMaxImpactXY(DCAxy);
	    //dtc1etaphitpc[aniter]->SetMaxImpactXYPtDep(0.0182, 0.0350, -1.01);
	    dtc1etaphitpc[aniter]->SetMaxImpactZ(DCAz);  //3.0
	    //      dtc1etaphitpc[aniter]->SetMaxSigmaToVertex(6.0);
	  }

	  

	  cutPass1YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutPass1%stpcM%i", chrgs[ichg], imult), 0.493677);
	  cutFail1YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutFail1%stpcM%i", chrgs[ichg], imult), 0.493677);
	  dtc1etaphitpc[aniter]->AddCutMonitor(cutPass1YPtetaphitpc[aniter], cutFail1YPtetaphitpc[aniter]);
	  
	  cutPass1PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cutPass1%stpcM%i", chrgs[ichg], imult),1);
	  cutFail1PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cutFail1%stpcM%i", chrgs[ichg], imult),1);
	  dtc1etaphitpc[aniter]->AddCutMonitor(cutPass1PIDetaphitpc[aniter], cutFail1PIDetaphitpc[aniter]);

	  
	 // sqpcetaphitpc[aniter] = new AliFemtoPairCutAntiGamma();
	//  sqpcetaphitpc[aniter] = new AliFemtoShareQualityTPCEntranceSepPairCut();
	  
          if (ichg < 2) {
//	    sqpcetaphitpc[aniter] = new AliFemtoPairCutRadialDistance();//AliFemto dphi* cut
    	  sqpcetaphitpc[aniter] = new AliFemtoPairCutRadialDistanceKK();  //Dhevan's dphi* cut
	  if (runtype == 0) {
	    sqpcetaphitpc[aniter]->SetShareQualityMax(1.0);
	    sqpcetaphitpc[aniter]->SetShareFractionMax(0.05);
	    sqpcetaphitpc[aniter]->SetRemoveSameLabel(kFALSE);
	    // sqpcetaphitpc[aniter]->SetMaxEEMinv(0.0);
	    // sqpcetaphitpc[aniter]->SetMaxThetaDiff(0.0);
	    	    //ml sqpcetaphitpc[aniter]->SetTPCEntranceSepMinimum(1.5);
//ml	    sqpcetaphitpc[aniter]->SetRadialDistanceMinimum(0.12, 0.03);
//ml	    sqpcetaphitpc[aniter]->SetEtaDifferenceMinimum(0.02);

        
	   //--------- AK cuts ----------->>>>
        /*
            sqpcetaphitpc[aniter]->SetMinimumRadius(1.2); //0.8
            sqpcetaphitpc[aniter]->SetPhiStarMin(kFALSE);
     	    sqpcetaphitpc[aniter]->SetPhiStarDifferenceMinimum(PhiStarDifferenceMinimum);
            sqpcetaphitpc[aniter]->SetEtaDifferenceMinimum(EtaDifferenceMinimum);
           //sqpcetaphitpc[aniter]->SetMinimumRadius(0.8);//not need for AliFemtoPairCutRadialDistanceKK()
           */
           
	   //---------  eta-phi*  -----------<<<

	  
     	    sqpcetaphitpc[aniter]->SetPhiStarDifferenceMinimum(PhiStarDifferenceMinimum);
            sqpcetaphitpc[aniter]->SetEtaDifferenceMinimum(EtaDifferenceMinimum);
	 
	    
	  }
	  else if (runtype == 1) {
	    sqpcetaphitpc[aniter]->SetShareQualityMax(1.0);
	    sqpcetaphitpc[aniter]->SetShareFractionMax(1.05);
	    sqpcetaphitpc[aniter]->SetRemoveSameLabel(kFALSE);
	    // sqpcetaphitpc[aniter]->SetMaxEEMinv(0.002);
	    // sqpcetaphitpc[aniter]->SetMaxThetaDiff(0.008);
	    	    //ml sqpcetaphitpc[aniter]->SetTPCEntranceSepMinimum(5.0);
//	    sqpcetaphitpc[aniter]->SetRadialDistanceMinimum(1.2, 0.03);
//	    sqpcetaphitpc[aniter]->SetEtaDifferenceMinimum(0.02);

//  void SetPhiStarDifferenceMinimum(double dtpc);
//  void SetEtaDifferenceMinimum(double etpc);
//  void SetMinimumRadius(double minrad);
//  void SetMagneticFieldSign(int magsign);

	   //--------- AK cuts ----------->>>>
/*
            sqpcetaphitpc[aniter]->SetMinimumRadius(1.2); //0.8
            sqpcetaphitpc[aniter]->SetPhiStarMin(kFALSE);
     	    sqpcetaphitpc[aniter]->SetPhiStarDifferenceMinimum(PhiStarDifferenceMinimum);
            sqpcetaphitpc[aniter]->SetEtaDifferenceMinimum(EtaDifferenceMinimum);
            */
           //sqpcetaphitpc[aniter]->SetMinimumRadius(0.8);//not need for AliFemtoPairCutRadialDistanceKK()
	   //--------- km:  eta-phi* Dhevan's custs -----------<<<

         /////////sqpcetaphitpc[aniter]->SetMagneticFieldSign(1);
     //sqpcetaphitpc[aniter]->SetMagneticFieldSign(1.0);

	    sqpcetaphitpc[aniter]->SetPhiStarDifferenceMinimum(PhiStarDifferenceMinimum);
            sqpcetaphitpc[aniter]->SetEtaDifferenceMinimum(EtaDifferenceMinimum);


	  }
	  else if (runtype == 2) {
            sqpcetaphitpc[aniter]->SetDataType(AliFemtoPairCut::kAOD);
	    sqpcetaphitpc[aniter]->SetShareQualityMax(1.0);
	    sqpcetaphitpc[aniter]->SetShareFractionMax(0.05);
	    sqpcetaphitpc[aniter]->SetRemoveSameLabel(kFALSE);
	    // sqpcetaphitpc[aniter]->SetMaxEEMinv(0.0);
	    // sqpcetaphitpc[aniter]->SetMaxThetaDiff(0.0);
	    	//ml    sqpcetaphitpc[aniter]->SetTPCEntranceSepMinimum(0.0);
//ml	    sqpcetaphitpc[aniter]->SetRadialDistanceMinimum(1.2, 0.045);
//ml	    sqpcetaphitpc[aniter]->SetEtaDifferenceMinimum(0.016);

	   //--------- km:  eta-phi* Dhevan's custs ----------->>>>
        
        /*     
            sqpcetaphitpc[aniter]->SetMinimumRadius(1.2); //0.8
            sqpcetaphitpc[aniter]->SetPhiStarMin(kFALSE);
            sqpcetaphitpc[aniter]->SetPhiStarDifferenceMinimum(PhiStarDifferenceMinimum);//0.045); // 0.012 - pions, 0.017 - kaons, 0.018
            sqpcetaphitpc[aniter]->SetEtaDifferenceMinimum(EtaDifferenceMinimum);//0.01);
        */
            
	    sqpcetaphitpc[aniter]->SetPhiStarDifferenceMinimum(PhiStarDifferenceMinimum);
            sqpcetaphitpc[aniter]->SetEtaDifferenceMinimum(EtaDifferenceMinimum);

	  }
         }        

	  
	  
	  anetaphitpc[aniter]->SetEventCut(mecetaphitpc[aniter]);
	  anetaphitpc[aniter]->SetFirstParticleCut(dtc1etaphitpc[aniter]);//K+
	  anetaphitpc[aniter]->SetSecondParticleCut(dtc2etaphitpc[aniter]);//K-
	  anetaphitpc[aniter]->SetPairCut(sqpcetaphitpc[aniter]);


	  if (runktdep) {
	    int ktm;
	    for (int ikt=0; ikt<cKt/*8*/; ikt++) {
	      ktm = aniter*cKt/*8*/ + ikt;
	      ktpcuts[ktm] = new AliFemtoKTPairCut(ktrng[ikt], ktrng[ikt+1]);
	      
              cqinvkttpc[ktm] = new AliFemtoQinvCorrFctn(Form("cqinv%stpcM%ikT%i", chrgs[ichg], imult, ikt),nbinssh,0.0,(imult>6)?shqmax*2.5:shqmax);

	      cqinvkttpc[ktm]->SetPairSelectionCut(ktpcuts[ktm]);//add kT bins
	      anetaphitpc[aniter]->AddCorrFctn(cqinvkttpc[ktm]);// add CF histos

	      //~~~
	      cqinvsqtpc[ktm] = new AliFemtoShareQualityCorrFctn(Form("cqinvsq%stpcM%ikT%i", chrgs[ichg], imult, ikt),nbinssh,0.0,shqmax);
	      cqinvsqtpc[ktm]->SetPairSelectionCut(ktpcuts[ktm]);
	      anetaphitpc[aniter]->AddCorrFctn(cqinvsqtpc[ktm]);

           // model--------------
           
             cqinvkttpcmodel[ktm] = new AliFemtoModelCorrFctn(Form("cqinvModel%stpcM%ikT%i", chrgs[ichg], imult, ikt),100, 0.0,0.2);

	      cqinvkttpcmodel[ktm]->SetPairSelectionCut(ktpcuts[ktm]);//add kT bins
	      cqinvkttpcmodel[ktm]->SetKaonPDG(kTRUE);//Special MC analysis for K selected by PDG code -->
              cqinvkttpcmodel[ktm]->ConnectToManager(tModelManager);
	      anetaphitpc[aniter]->AddCorrFctn(cqinvkttpcmodel[ktm]);// add CF histos
           

	    }
	  }
	  
	  // cdedpetaphi[aniter] = new AliFemtoCorrFctnDEtaDPhi(Form("cdedp%stpcM%i", chrgs[ichg], imult),24, 24);
	  // anetaphitpc[aniter]->AddCorrFctn(cdedpetaphi[aniter]);
	  
	  Manager->AddAnalysis(anetaphitpc[aniter]);	
	}
      }
    }
  }
  // *** End Kaon-Kaon (positive) analysis

  return Manager;
}                         
                      
