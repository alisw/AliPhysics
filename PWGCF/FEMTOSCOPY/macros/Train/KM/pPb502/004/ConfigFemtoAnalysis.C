// Konstantin Mikhaylov: 30-March-2023
// 004: PWGCF/FEMTOSCOPY/macros/Train/KM/pPb502/004/ConfigFemtoAnalysis.C
// second wagon with AliFemtoKpm45TrackCut
// instead of Lena's wagon AliFemtoKKTrackCut.h
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#if !defined(__CINT__) || defined(__MAKECINT_)
#include "AliFemtoManager.h"
//#include "AliFemtoEventReaderESDChain.h"
//#include "AliFemtoEventReaderESDChainKine.h"
//#include "AliFemtoEventReaderAODChain.h"
#include "AliFemtoEventReaderAOD.h"
#include "AliFemtoEventReaderAODMultSelection.h"
#include "AliFemtoSimpleAnalysis.h"
#include "AliFemtoBasicEventCut.h"
#include "AliFemtoSphericityEventCut.h"
//#include "AliFemtoSpherocityEventCut.h"
#include "AliFemtoESDTrackCut.h"
#include "AliFemtoKpm45TrackCut.h"
//#include "AliFemtoKKTrackCut.h"
//#include "AliFemtoKKTrackCutTest.h"
#include "AliFemtoCorrFctn.h"
#include "AliFemtoCutMonitorParticleYPt.h"
#include "AliFemtoCutMonitorParticleVertPos.h"
#include "AliFemtoCutMonitorParticleMomRes.h"
#include "AliFemtoCutMonitorParticlePID.h"
#include "AliFemtoCutMonitorEventMult.h"
#include "AliFemtoCutMonitorEventVertex.h"
#include "AliFemtoCutMonitorEventSphericity.h"
#include "AliFemtoShareQualityTPCEntranceSepPairCut.h"
#include "AliFemtoPairCutAntiGamma.h"
#include "AliFemtoPairCutRadialDistance.h"
#include "AliFemtoQinvCorrFctn.h"
#include "AliFemtoShareQualityCorrFctn.h"
#include "AliFemtoPairCutRadialDistanceKKdist.h"
#include "AliFemtoTPCInnerCorrFctn.h"
#include "AliFemtoVertexMultAnalysis.h"
#include "AliFemtoCorrFctn3DSpherical.h"
#include "AliFemtoChi2CorrFctn.h"
#include "AliFemtoCorrFctnTPCNcls.h"
#include "AliFemtoBPLCMS3DCorrFctn.h"
#include "AliFemtoCorrFctn3DLCMSSym.h"
#include "AliFemtoBPLCMS3DCorrFctnKK.h"
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
#include "AliFemtoCutMonitorParticlePtPDG.h"
#include "AliFemtoKTPairCut.h"
#include "AliFemtoCutMonitorCollections.h"
#endif

//________________________________________________________________________
AliFemtoManager* ConfigFemtoAnalysis() {

  double PionMass = 0.13956995;
  double KaonMass = 0.493677;
  
  Double_t DCAxy=0.135;//cm [0.135 - main]
  Double_t DCAz =0.130;//cm [0.130 - main]
	
  //multiplicity bins
  int runmults[3] = {1, 1, 1};
  double multbins[4] = {0.01, 200, 400, 900};
  //int runmults[6] = {1, 1, 1, 1, 1, 1};
  //int multbins[7] = {0.01, 50, 100, 200, 400, 600, 800};

  int runch[2] = {1, 1};
  const char *chrgs[2] = { "Kp", "Km"};

  int runktdep = 1;
  double ktrng[3] = {0.2, 0.5, 1.0};
  //double ktrng[9] = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.3};

  int run3d = 0;
  //int runshlcms = 1;
  int runshlcms = 0;
  double shqmax;
  int nbinssh = 200;
  //int nbinssh = 100;

  //if (runshlcms) shqmax = 2.0;
  if (runshlcms) shqmax = 0.25;
  else shqmax = 2.0;

  ///AliFemtoEventReaderAODChain *Reader = new AliFemtoEventReaderAODChain();
  ///Reader->SetUseMultiplicity(AliFemtoEventReaderAODChain::kCentrality);
  ///Reader->SetFilterBit(7);
  //Reader->SetCentralityPreSelection(0.001, 900);
  //Reader->SetDCAglobalTrack(kTRUE);
  //Reader->SetpA2013(kTRUE);

  //AliFemtoEventReaderAOD *Reader = new AliFemtoEventReaderAODMultSelection();
  AliFemtoEventReaderAODMultSelection *Reader = new AliFemtoEventReaderAODMultSelection();
  Reader->SetFilterBit(8);
  Reader->SetReadV0(1);
  Reader->SetEPVZERO(kTRUE);
  Reader->SetCentralityFlattening(kTRUE);
  Reader->SetUseMultiplicity(AliFemtoEventReaderAODMultSelection::kCentralityV0A);
  Reader->SetDCAglobalTrack(kTRUE);
  //Reader->SetUseMultiplicity(AliFemtoEventReaderAOD::kCentrality);
  /*Reader->SetCentralityFlattening(kFALSE);
  Reader->SetReadV0(0);
  Reader->SetDCAglobalTrack(kTRUE);*/


  /*AliFemtoEventReaderAOD *Reader = new AliFemtoEventReaderAODMultSelection();
  Reader->SetFilterBit(7);
  Reader->SetEPVZERO(kTRUE);
  Reader->SetUseMultiplicity(AliFemtoEventReaderAOD::kCentrality);
  Reader->SetCentralityFlattening(kFALSE);
  Reader->SetReadV0(0);//from ML*/
  
  AliFemtoManager* Manager=new AliFemtoManager();
  Manager->SetEventReader(Reader);

  AliFemtoVertexMultAnalysis    *anetaphitpc[20];
  AliFemtoBasicEventCut         *mecetaphitpc[20];
  //AliFemtoSphericityEventCut         *mecetaphitpc[20];
  //AliFemtoSpherocityEventCut         *mecetaphitpc[20];
  AliFemtoCutMonitorEventMult   *cutPassEvMetaphitpc[20];
  AliFemtoCutMonitorEventMult   *cutFailEvMetaphitpc[20];
  AliFemtoCutMonitorEventVertex *cutPassEvVetaphitpc[20];
  AliFemtoCutMonitorEventVertex *cutFailEvVetaphitpc[20];
  AliFemtoCutMonitorCollections   *cutPassColletaphitpc[20];
  AliFemtoCutMonitorCollections   *cutFailColletaphitpc[20];
  //AliFemtoKKTrackCut           *dtc1etaphitpc[20];
  //AliFemtoKKTrackCut           *dtc2etaphitpc[20];
  AliFemtoKpm45TrackCut           *dtc1etaphitpc[20];
  AliFemtoKpm45TrackCut           *dtc2etaphitpc[20];
  //AliFemtoESDTrackCut           *dtc1etaphitpc[20];
  //AliFemtoESDTrackCut           *dtc2etaphitpc[20];
  //AliFemtoKKTrackCutTest        *dtc1etaphitpc[20];
  //AliFemtoKKTrackCutTest        *dtc2etaphitpc[20];
  AliFemtoCutMonitorParticleYPt *cutPass1YPtetaphitpc[20];
  AliFemtoCutMonitorParticleYPt *cutFail1YPtetaphitpc[20];
  AliFemtoCutMonitorParticlePID *cutPass1PIDetaphitpc[20];
  AliFemtoCutMonitorParticlePID *cutFail1PIDetaphitpc[20];
  AliFemtoCutMonitorParticleYPt *cutPass2YPtetaphitpc[20];
  AliFemtoCutMonitorParticleYPt *cutFail2YPtetaphitpc[20];
  AliFemtoCutMonitorParticlePID *cutPass2PIDetaphitpc[20];
  AliFemtoCutMonitorParticlePID *cutFail2PIDetaphitpc[20];
  AliFemtoCutMonitorEventSphericity *cutPassEvSpher[20];
  AliFemtoCutMonitorEventSphericity *cutFailEvSpher[20];
  //AliFemtoPairCutAntiGamma      *sqpcetaphitpc[20];
  //AliFemtoShareQualityPairCut      *sqpcetaphitpc[20];
  //AliFemtoPairCutRadialDistance      *sqpcetaphitpc[20];
  AliFemtoPairCutRadialDistanceKKdist  *sqpcetaphitpc[20];
  AliFemtoCorrFctnDirectYlm     *cylmetaphitpc[20];
  AliFemtoCorrFctnDEtaDPhi      *cdedpetaphi[20];
  AliFemtoChi2CorrFctn          *cchiqinvetaphitpc[20];
  //AliFemtoKTPairCut             *ktpcuts[20*7];
  //AliFemtoCorrFctnDirectYlm     *cylmkttpc[20*7];
  //AliFemtoQinvCorrFctn          *cqinvkttpc[20*7];
  //AliFemtoCorrFctn3DLCMSSym     *cq3dlcmskttpc[20*7];
  AliFemtoKTPairCut             *ktpcuts[20*8];
  AliFemtoCorrFctnDirectYlm     *cylmkttpc[20*8];
  AliFemtoQinvCorrFctn          *cqinvkttpc[20*8];
  //AliFemtoCorrFctn3DLCMSSym     *cq3dlcmskttpc[20*8];
  AliFemtoBPLCMS3DCorrFctnKK     *cq3dlcmskttpc[20*8];
  AliFemtoCorrFctnTPCNcls       *cqinvnclstpc[20];
  AliFemtoShareQualityCorrFctn  *cqinvsqtpc[20*10];
  AliFemtoChi2CorrFctn          *cqinvchi2tpc[20];
  AliFemtoTPCInnerCorrFctn      *cqinvinnertpc[20*10];

  // *** Begin pion-pion analysis ***
  int aniter = 0;

  for (int imult=0; imult<3; imult++) {
    if (runmults[imult]) {
      for (int ichg=0; ichg<1/*2*/; ichg++) {//'1' for K+K-
	if (runch[ichg]) {
	  aniter = ichg*3+imult;

	  anetaphitpc[aniter] = new AliFemtoVertexMultAnalysis(10, -10.0, 10.0, 4, multbins[imult], multbins[imult+1]);
	  anetaphitpc[aniter]->SetNumEventsToMix(3);
	  anetaphitpc[aniter]->SetMinSizePartCollection(1);
	  anetaphitpc[aniter]->SetVerboseMode(kFALSE);
      
	  mecetaphitpc[aniter] = new AliFemtoBasicEventCut();
	  mecetaphitpc[aniter]->SetEventMult(0,10000);
	  mecetaphitpc[aniter]->SetVertZPos(-10,10);
	  
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
	  
	  cutPassColletaphitpc[aniter] = new AliFemtoCutMonitorCollections(Form("cutPass%stpcM%i", chrgs[ichg], imult));
          cutFailColletaphitpc[aniter] = new AliFemtoCutMonitorCollections(Form("cutFail%stpcM%i", chrgs[ichg], imult));
          mecetaphitpc[aniter]->AddCutMonitor(cutPassColletaphitpc[aniter], cutFailColletaphitpc[aniter]);
          
          cutPassEvSpher[aniter] = new AliFemtoCutMonitorEventSphericity(Form("cutPass%stpcM%i", chrgs[ichg], imult));
          cutFailEvSpher[aniter] = new AliFemtoCutMonitorEventSphericity(Form("cutFail%stpcM%i", chrgs[ichg], imult));
          mecetaphitpc[aniter]->AddCutMonitor(cutPassEvSpher[aniter], cutFailEvSpher[aniter]);


	  //dtc1etaphitpc[aniter] = new AliFemtoKKTrackCut();//K+
	  //dtc2etaphitpc[aniter] = new AliFemtoKKTrackCut();//K-
	  dtc1etaphitpc[aniter] = new AliFemtoKpm45TrackCut();//K+
	  dtc2etaphitpc[aniter] = new AliFemtoKpm45TrackCut();//K-

	  //--- K+- --->
	  dtc1etaphitpc[aniter]->SetCharge(1.0);
	  dtc2etaphitpc[aniter]->SetCharge(-1.0);
	  //--- K+- ---<

	   dtc1etaphitpc[aniter]->SetPt(0.14,1.5);//K+
	   dtc1etaphitpc[aniter]->SetEta(-0.8,0.8);//K+
	   dtc2etaphitpc[aniter]->SetPt(0.14,1.5);//K-
	   dtc2etaphitpc[aniter]->SetEta(-0.8,0.8);//K-

	  //PID method
	   //dtc1etaphitpc[aniter]->SetMass(PionMass);
	   dtc1etaphitpc[aniter]->SetMass(KaonMass);//K+
	   //dtc1etaphitpc[aniter]->SetMostProbablePion();
	   dtc1etaphitpc[aniter]->SetMostProbableKaon();//K+
	   //dtc1etaphitpc[aniter]->SetPIDMethod(AliFemtoESDTrackCut::kContour);
	   dtc2etaphitpc[aniter]->SetMass(KaonMass);//K-
	   dtc2etaphitpc[aniter]->SetMostProbableKaon();//K-

	   
//------------------- November 2013 -----------------------------------< 
	  // new cuts to remove electron (do not take into analysis if 400<p<500)
	   //K+++++++++:
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
	   //K----------:
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
	 //K- contour (default p0=9999, p1=0, p2=0):  
	 //dtc2etaphitpc[aniter]->SetKaonTCPdEdxContourP0( 196.25);//december 2014
	 //dtc2etaphitpc[aniter]->SetKaonTCPdEdxContourP1(-420.00);//december 2014
	 //dtc2etaphitpc[aniter]->SetKaonTCPdEdxContourP2( 300.00);//december 2014
	  //------------------- November 2013 ----------------------------------->
	   //Track quality cuts
	   //K+:
	   dtc1etaphitpc[aniter]->SetminTPCncls(80);
	   dtc1etaphitpc[aniter]->SetRemoveKinks(kTRUE);
	   dtc1etaphitpc[aniter]->SetMaxTPCChiNdof(4.0);
	   //dtc1etaphitpc[aniter]->SetMaxITSChiNdof(36);	  
	   dtc1etaphitpc[aniter]->SetLabel(kFALSE);
	   dtc1etaphitpc[aniter]->SetMaxImpactZ(DCAz);//main
	   dtc1etaphitpc[aniter]->SetMaxImpactXY(DCAxy);//main
	  //dtc1etaphitpc[aniter]->SetMaxImpactZ(0.14);
	  //dtc1etaphitpc[aniter]->SetMaxImpactXY(0.145);
	  //dtc1etaphitpc[aniter]->SetMaxImpactXYPtDep(0.0105, 0.0350, -1.1);
	  //dtc1etaphitpc[aniter]->SetMaxImpactXYPtDep(0.0182, 0.0350, -1.01);
	  //dtc1etaphitpc[aniter]->SetMaxSigmaToVertex(6.0);
	   //K-:
	   dtc2etaphitpc[aniter]->SetminTPCncls(80);
	   dtc2etaphitpc[aniter]->SetRemoveKinks(kTRUE);
	   dtc2etaphitpc[aniter]->SetMaxTPCChiNdof(4.0);
	   //dtc2etaphitpc[aniter]->SetMaxITSChiNdof(36);	  
	   dtc2etaphitpc[aniter]->SetLabel(kFALSE);
	   dtc2etaphitpc[aniter]->SetMaxImpactZ(DCAz);//main
	   dtc2etaphitpc[aniter]->SetMaxImpactXY(DCAxy);//main
	  //dtc2etaphitpc[aniter]->SetMaxImpactZ(0.14);
	  //dtc2etaphitpc[aniter]->SetMaxImpactXY(0.145);
	  //dtc2etaphitpc[aniter]->SetMaxImpactXYPtDep(0.0105, 0.0350, -1.1);
	  //dtc2etaphitpc[aniter]->SetMaxImpactXYPtDep(0.0182, 0.0350, -1.01);
	  //dtc2etaphitpc[aniter]->SetMaxSigmaToVertex(6.0);

	  //cutPass1YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutPass1%stpcM%i", chrgs[ichg], imult), 0.13957);
	  //cutFail1YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutFail1%stpcM%i", chrgs[ichg], imult), 0.13957);
	  cutPass1YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutPass1%stpcM%i", chrgs[ichg], imult), 0.493677);
	  cutFail1YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutFail1%stpcM%i", chrgs[ichg], imult), 0.493677);
	  dtc1etaphitpc[aniter]->AddCutMonitor(cutPass1YPtetaphitpc[aniter], cutFail1YPtetaphitpc[aniter]);
	  
	  //cutPass1PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cutPass1%stpcM%i", chrgs[ichg], imult),0);
	  //cutFail1PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cutFail1%stpcM%i", chrgs[ichg], imult),0);
	  cutPass1PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cutPass1%stpcM%i", chrgs[ichg], imult),1);
	  cutFail1PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cutFail1%stpcM%i", chrgs[ichg], imult),1);
	  dtc1etaphitpc[aniter]->AddCutMonitor(cutPass1PIDetaphitpc[aniter], cutFail1PIDetaphitpc[aniter]);
	  
	  //dtc1etaphitpc[aniter]->SetStatus(AliESDtrack::kTPCin); ///my

 
	    
	  //sqpcetaphitpc[aniter] = new AliFemtoPairCutAntiGamma();
	  //sqpcetaphitpc[aniter] = new AliFemtoPairCutRadialDistance();
	  //sqpcetaphitpc[aniter] = new AliFemtoShareQualityTPCEntranceSepPairCut();
	  //sqpcetaphitpc[aniter] = new AliFemtoPairCutRadialDistance();

	  //sqpcetaphitpc[aniter] = new AliFemtoShareQualityPairCut();
	  //sqpcetaphitpc[aniter]->SetShareQualityMax(1.0);
	  //sqpcetaphitpc[aniter]->SetShareFractionMax(0.05);

	  sqpcetaphitpc[aniter] = new AliFemtoPairCutRadialDistanceKKdist();
	  sqpcetaphitpc[aniter]->SetShareQualityMax(1.0);
	  sqpcetaphitpc[aniter]->SetShareFractionMax(0.05);
	  sqpcetaphitpc[aniter]->SetAverageSeparation(3.0);

	  sqpcetaphitpc[aniter]->SetRemoveSameLabel(kFALSE);
	  
	 // sqpcetaphitpc[aniter]->SetPhiStarDifferenceMinimum(0.04);
         // sqpcetaphitpc[aniter]->SetEtaDifferenceMinimum(0.02);
	    
	  //sqpcetaphitpc[aniter]->SetEtaDifferenceMinimum(0.016);
	  //sqpcetaphitpc[aniter]->SetPhiStarDifferenceMinimum(0.02);
	  //runtype==0
	  // sqpcetaphitpc[aniter]->SetMaxEEMinv(0.0);
	  // sqpcetaphitpc[aniter]->SetMaxThetaDiff(0.0);
	  // sqpcetaphitpc[aniter]->SetTPCEntranceSepMinimum(1.5);
	  // sqpcetaphitpc[aniter]->SetPhiStarDistanceMinimum(0.03);
	  // sqpcetaphitpc[aniter]->SetRadialDistanceMinimum(0.12, 0.03);
	  //sqpcetaphitpc[aniter]->SetEtaDifferenceMinimum(0.02);
	  //runtype==1
	  //	    sqpcetaphitpc[aniter]->SetTPCEntranceSepMinimum(5.0);
	  //	    sqpcetaphitpc[aniter]->SetPhiStarDistanceMinimum(0.03);
	  //sqpcetaphitpc[aniter]->SetRadialDistanceMinimum(1.2, 0.03);
	  //sqpcetaphitpc[aniter]->SetEtaDifferenceMinimum(0.02);
	  //runtype==2
	  //	    sqpcetaphitpc[aniter]->SetTPCEntranceSepMinimum(1.0);
	  //	    sqpcetaphitpc[aniter]->SetPhiStarDistanceMinimum(0.03);
	  //sqpcetaphitpc[aniter]->SetRadialDistanceMinimum(1.2, 0.045);
	  //sqpcetaphitpc[aniter]->SetEtaDifferenceMinimum(0.016);
	  //sqpcetaphitpc[aniter]->SetPhiStarDifferenceMinimum(0.02);
	  
	  anetaphitpc[aniter]->SetEventCut(mecetaphitpc[aniter]);
	  anetaphitpc[aniter]->SetFirstParticleCut(dtc1etaphitpc[aniter]);//K+
	  anetaphitpc[aniter]->SetSecondParticleCut(dtc2etaphitpc[aniter]);//K-//anetaphitpc[aniter]->SetSecondParticleCut(dtc1etaphitpc[aniter]);
	  anetaphitpc[aniter]->SetPairCut(sqpcetaphitpc[aniter]);
	  
	  //Correlation functions

	  //Spherical harmonics (without kT bins)
	  //cylmetaphitpc[aniter] = new AliFemtoCorrFctnDirectYlm(Form("cylm%stpcM%i", chrgs[ichg], imult),3,nbinssh,0.0,shqmax,runshlcms);
	  //anetaphitpc[aniter]->AddCorrFctn(cylmetaphitpc[aniter]);

	  //Qinv (without kT bins)-->
	  cqinvkttpc[aniter] = new AliFemtoQinvCorrFctn(Form("cqinv%stpcM%i", chrgs[ichg], imult),nbinssh,0.0,shqmax);
	  anetaphitpc[aniter]->AddCorrFctn(cqinvkttpc[aniter]);
	  //Qinv (without kT bins)--<
	  
	  // cqinvnclstpc[aniter] = new AliFemtoCorrFctnTPCNcls(Form("cqinvncls%stpcM%i", chrgs[ichg], imult),nbinssh,0.0,shqmax);
	  // anetaphitpc[aniter]->AddCorrFctn(cqinvnclstpc[aniter]);

	  // cqinvchi2tpc[aniter] = new AliFemtoChi2CorrFctn(Form("cqinvchi2%stpcM%i", chrgs[ichg], imult),nbinssh,0.0,shqmax);
	  // anetaphitpc[aniter]->AddCorrFctn(cqinvchi2tpc[aniter]);

	  if (runktdep) {
	    int ktm;
	    for (int ikt=0; ikt<2; ikt++) {
	      ktm = aniter*2 + ikt;
	      ktpcuts[ktm] = new AliFemtoKTPairCut(ktrng[ikt], ktrng[ikt+1]);
	      	      
	      cqinvkttpc[ktm] = new AliFemtoQinvCorrFctn(Form("cqinv%stpcM%ikT%i", chrgs[ichg], imult, ikt),nbinssh,0.0, shqmax);
	      cqinvkttpc[ktm]->SetPairSelectionCut(ktpcuts[ktm]);
	      anetaphitpc[aniter]->AddCorrFctn(cqinvkttpc[ktm]);

	      cqinvsqtpc[ktm] = new AliFemtoShareQualityCorrFctn(Form("cqinvsq%stpcM%ikT%i", chrgs[ichg], imult, ikt),nbinssh,0.0,shqmax);
	      cqinvsqtpc[ktm]->SetPairSelectionCut(ktpcuts[ktm]);
	      anetaphitpc[aniter]->AddCorrFctn(cqinvsqtpc[ktm]);

	      cqinvinnertpc[ktm] = new AliFemtoTPCInnerCorrFctn(Form("cqinvinner%stpcM%ikT%i", chrgs[ichg], imult, ikt),nbinssh,0.0,shqmax);
	      cqinvinnertpc[ktm]->SetPairSelectionCut(ktpcuts[ktm]);
	      cqinvinnertpc[ktm]->SetRadius(1.2);
	      anetaphitpc[aniter]->AddCorrFctn(cqinvinnertpc[ktm]);

	    }
	  }
	  
	  /*cdedpetaphi[aniter] = new AliFemtoCorrFctnDEtaDPhi(Form("cdedp%stpcM%i", chrgs[ichg], imult),39, 39);
	    anetaphitpc[aniter]->AddCorrFctn(cdedpetaphi[aniter]);*/
	  
	  Manager->AddAnalysis(anetaphitpc[aniter]);	
	}
      }
    }
  }
  // *** End K+K- analysis

  return Manager;
}                         
                      
