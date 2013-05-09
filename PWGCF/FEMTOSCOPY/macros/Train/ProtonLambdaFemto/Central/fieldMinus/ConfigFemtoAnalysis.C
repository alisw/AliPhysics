/*********************************************************************
 *                                                                   *
 * ConfigFemtoAnalysis.C - configuration macro for the femtoscopic   *
 * analysis, meant as a QA process for two-particle effects          *
 *                                                                   *
 * Author: Adam Kisiel (Adam.Kisiel@cern.ch)                         *
 *                                                                   *
 *********************************************************************/

#if !defined(__CINT__) || defined(__MAKECINT_)
#include "AliFemtoManager.h"
#include "AliFemtoEventReaderESDChain.h"
#include "AliFemtoEventReaderESDChainKine.h"
#include "AliFemtoSimpleAnalysis.h"
#include "AliFemtoBasicEventCut.h"
#include "AliFemtoESDTrackCut.h"
#include "AliFemtoCorrFctn.h"
#include "AliFemtoCutMonitorParticleYPt.h"
#include "AliFemtoCutMonitorParticleVertPos.h"
#include "AliFemtoCutMonitorParticleMomRes.h"
#include "AliFemtoCutMonitorParticlePID.h"
#include "AliFemtoCutMonitorEventMult.h"
#include "AliFemtoCutMonitorEventVertex.h"
#include "AliFemtoShareQualityTPCEntranceSepPairCut.h"
#include "AliFemtoPairCutAntiGamma.h"
#include "AliFemtoQinvCorrFctn.h"
#include "AliFemtoShareQualityCorrFctn.h"
#include "AliFemtoTPCInnerCorrFctn.h"
#include "AliFemtoVertexMultAnalysis.h"
#include "AliFemtoCorrFctn3DSpherical.h"
#include "AliFemtoChi2CorrFctn.h"
#include "AliFemtoCorrFctnTPCNcls.h"
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
#include "AliFemtoEnumeration.h"
#include "AliFemtoV0PairCut.h"
#include "AliFemtoV0TrackPairCut.h"
#include "AliFemtoV0TrackCut.h"
#include "AliFemtoEventReaderAODChain.h"
#include "AliFemtoCorrFctnNonIdDR.h"
#include "AliFemtoCutMonitorCollections.h"
#include "AliFemtoPairCutRadialDistance.h"
#include "AliFemtoAvgSepCorrFctn.h"
#endif

//________________________________________________________________________
AliFemtoManager* ConfigFemtoAnalysis() {

  double PionMass = 0.13956995;
  double KaonMass = 0.493677;
  double ProtonMass = 0.938272;
  double LambdaMass = 1.115683;

  double psi = TMath::Pi()/2.;
  double psid = TMath::Pi()/6.;

  // Switches for QA analyses
  int runmults[10] = {1, 1, 0, 0, 0, 0, 0, 0, 0, 0};
  int multbins[11] = {0.00001, 100, 200, 300, 400, 500, 600, 700, 800, 900, 900};
  int runch[10] = {1, 1, 1, 1, 1, 1, 1, 0, 0, 0};
  const char *chrgs[10] = { "V0LL", "V0ALAL", "V0LAL", "V0PL", "V0APL", "V0PAL", "V0APAL","PP","PAP","APAP" };

 
  int runktdep = 0;
  double ktrng[3] = {0.01, 0.7, 100};

  int numOfMultBins = 10;
  int numOfChTypes = 10;
  int numOfkTBins = 2;


  int runshlcms = 1;

  int runtype = 2; // Types 0 - global, 1 - ITS only, 2 - TPC Inner
 
  double shqmax;
  int nbinssh = 200;

  shqmax = 1.0;

  AliFemtoEventReaderAODChain* Reader=new AliFemtoEventReaderAODChain();
  Reader->SetFilterBit(7);
  Reader->SetCentralityPreSelection(0.00001, 910);
  Reader->SetReadV0(1); //Read V0
  Reader->SetMagneticFieldSign(-1.0); //-1->field1, 1->field3

  AliFemtoManager* Manager=new AliFemtoManager();
  Manager->SetEventReader(Reader);


  AliFemtoVertexMultAnalysis    *anetaphitpc[320];
  AliFemtoBasicEventCut         *mecetaphitpc[320];
  AliFemtoCutMonitorEventMult   *cutPassEvMetaphitpc[320];
  AliFemtoCutMonitorEventMult   *cutFailEvMetaphitpc[320];
  AliFemtoCutMonitorEventVertex *cutPassEvVetaphitpc[320];
  AliFemtoCutMonitorEventVertex *cutFailEvVetaphitpc[320];
  AliFemtoCutMonitorCollections   *cutPassColletaphitpc[320];
  AliFemtoCutMonitorCollections   *cutFailColletaphitpc[320];
  AliFemtoV0TrackCut           *dtc1etaphitpc[320];
  AliFemtoV0TrackCut           *dtc2etaphitpc[320];
  AliFemtoESDTrackCut           *dtc3etaphitpc[320];
  AliFemtoESDTrackCut           *dtc4etaphitpc[320];
  AliFemtoCutMonitorParticleYPt *cutPass1YPtetaphitpc[320];
  AliFemtoCutMonitorParticleYPt *cutFail1YPtetaphitpc[320];
  AliFemtoCutMonitorParticlePID *cutPass1PIDetaphitpc[320];
  AliFemtoCutMonitorParticlePID *cutFail1PIDetaphitpc[320];
  AliFemtoCutMonitorParticleYPt *cutPass2YPtetaphitpc[320];
  AliFemtoCutMonitorParticleYPt *cutFail2YPtetaphitpc[320];
  AliFemtoCutMonitorParticlePID *cutPass2PIDetaphitpc[320];
  AliFemtoCutMonitorParticlePID *cutFail2PIDetaphitpc[320];
  AliFemtoCutMonitorV0          *cutPass1V0[320];
  AliFemtoCutMonitorV0          *cutFail1V0[320];
  AliFemtoCutMonitorV0          *cutPass2V0[320];
  AliFemtoCutMonitorV0          *cutFail2V0[320];
  AliFemtoV0PairCut             *sqp1cetaphitpc[320];
  AliFemtoV0TrackPairCut        *sqp2cetaphitpc[320];
  AliFemtoV0TrackPairCut        *sqp3cetaphitpc[320];
  AliFemtoV0TrackPairCut        *sqp4cetaphitpc[320];
  //AliFemtoPairCutAntiGamma      *sqp3cetaphitpc[320];
  AliFemtoPairCutRadialDistance *sqp5cetaphitpc[320];
  AliFemtoCorrFctnDirectYlm     *cylmetaphitpc[320];
  AliFemtoCorrFctnDEtaDPhi      *cdedpetaphi[320*6];
  AliFemtoChi2CorrFctn          *cchiqinvetaphitpc[320];
  AliFemtoKTPairCut             *ktpcuts[320*6];
  AliFemtoCorrFctnDirectYlm     *cylmkttpc[320*6];
  AliFemtoQinvCorrFctn          *cqinvkttpc[320*6];
  AliFemtoAvgSepCorrFctn        *avgsepcorr[320*6];
  //AliFemtoBPLCMS3DCorrFctn      *cq3dlcmskttpc[320*6];
  AliFemtoCorrFctnNonIdDR       *cnonidtpc[320*6];
  AliFemtoCorrFctnTPCNcls       *cqinvnclstpc[320];
  AliFemtoShareQualityCorrFctn  *cqinvsqtpc[320];
  AliFemtoTPCInnerCorrFctn      *cqinvtitpc[320];
  //AliFemtoChi2CorrFctn          *cqinvchi2tpc[320];
  AliFemtoQinvCorrFctn          *cqinvtpc[100*3];


  // *** Third QA task - HBT analysis with all pair cuts off, TPC only ***
  // *** Begin pion-pion (positive) analysis ***
  int aniter = 0;

  for (int imult=0; imult<numOfMultBins; imult++) {
    if (runmults[imult]) {
      for (int ichg=0; ichg<numOfChTypes; ichg++) {
	if (runch[ichg]) {
	      aniter = imult * numOfChTypes + ichg;



	      if(ichg==3 || ichg==4 || ichg==5 || ichg==6)
		runshlcms = 0;
	      else
		runshlcms = 1;

	      anetaphitpc[aniter] = new AliFemtoVertexMultAnalysis(10, -10.0, 10.0, 4, multbins[imult], multbins[imult+1]);
	      anetaphitpc[aniter]->SetNumEventsToMix(10);
	      anetaphitpc[aniter]->SetMinSizePartCollection(1);
	      anetaphitpc[aniter]->SetVerboseMode(kFALSE);

	      mecetaphitpc[aniter] = new AliFemtoBasicEventCut();
	      mecetaphitpc[aniter]->SetEventMult(0,100000);
	      mecetaphitpc[aniter]->SetVertZPos(-10,10);
	  
	      // mecetaphitpc->SetAcceptBadVertex(kTRUE);

	      cutPassEvMetaphitpc[aniter] = new AliFemtoCutMonitorEventMult(Form("cutPass%stpcM%i", chrgs[ichg], imult));
	      cutFailEvMetaphitpc[aniter] = new AliFemtoCutMonitorEventMult(Form("cutFail%stpcM%i", chrgs[ichg], imult));
	      mecetaphitpc[aniter]->AddCutMonitor(cutPassEvMetaphitpc[aniter], cutFailEvMetaphitpc[aniter]);
	  
	      cutPassEvVetaphitpc[aniter] = new AliFemtoCutMonitorEventVertex(Form("cutPass%stpcM%i", chrgs[ichg], imult));
	      cutFailEvVetaphitpc[aniter] = new AliFemtoCutMonitorEventVertex(Form("cutFail%stpcM%i", chrgs[ichg], imult));
	      mecetaphitpc[aniter]->AddCutMonitor(cutPassEvVetaphitpc[aniter], cutFailEvVetaphitpc[aniter]);

	      cutPassColletaphitpc[aniter] = new AliFemtoCutMonitorCollections(Form("cutPass%stpcM%i", chrgs[ichg], imult));
	      cutFailColletaphitpc[aniter] = new AliFemtoCutMonitorCollections(Form("cutFail%stpcM%i", chrgs[ichg], imult));
	      mecetaphitpc[aniter]->AddCutMonitor(cutPassColletaphitpc[aniter], cutFailColletaphitpc[aniter]);

	      //V0 first particle cut -> Lambda ichg 0, 2, 6
	      dtc1etaphitpc[aniter] = new AliFemtoV0TrackCut();
	      dtc1etaphitpc[aniter]->SetMass(LambdaMass);
	      dtc1etaphitpc[aniter]->SetEta(0.8); //0.8
	      dtc1etaphitpc[aniter]->SetPt(0.4,5.0); //0.4,100
	      dtc1etaphitpc[aniter]->SetEtaDaughters(0.8); //0.8
	      dtc1etaphitpc[aniter]->SetPtPosDaughter(0.7,5.0); //0.5
	      dtc1etaphitpc[aniter]->SetPtNegDaughter(0.16,5.0); //0.16
	      dtc1etaphitpc[aniter]->SetTPCnclsDaughters(80); //80
	      dtc1etaphitpc[aniter]->SetNdofDaughters(4.0); //4.0
	      dtc1etaphitpc[aniter]->SetStatusDaughters(AliESDtrack::kTPCrefit/* | AliESDtrack::kITSrefit*/);
	      dtc1etaphitpc[aniter]->SetOnFlyStatus(kFALSE);
	      dtc1etaphitpc[aniter]->SetParticleType(0);
	      dtc1etaphitpc[aniter]->SetMaxDcaV0Daughters(0.4); //1.5 Jai, 0.6 //0.4
	      dtc1etaphitpc[aniter]->SetMaxDcaV0(1.0); //5.0
	      dtc1etaphitpc[aniter]->SetMinDaughtersToPrimVertex(0.1); //0.01
	      dtc1etaphitpc[aniter]->SetMaxCosPointingAngle(0.998); //0.99 - Jai //0.998
	      dtc1etaphitpc[aniter]->SetInvariantMassLambda(1.112683,1.118683);


	      //V0 second particle cut -> AntiLambda ichg 1, 3, 4, 5
	      dtc2etaphitpc[aniter] = new AliFemtoV0TrackCut();
	      dtc2etaphitpc[aniter]->SetMass(LambdaMass);
	      dtc2etaphitpc[aniter]->SetEta(0.8);
	      dtc2etaphitpc[aniter]->SetPt(0.4,5.0);
	      dtc2etaphitpc[aniter]->SetEtaDaughters(0.8);
	      dtc2etaphitpc[aniter]->SetPtPosDaughter(0.16,5.0);
	      dtc2etaphitpc[aniter]->SetPtNegDaughter(0.7,5.0);
	      dtc2etaphitpc[aniter]->SetTPCnclsDaughters(80);
	      dtc2etaphitpc[aniter]->SetNdofDaughters(4.0); //4.0
	      dtc2etaphitpc[aniter]->SetStatusDaughters(AliESDtrack::kTPCrefit/* | AliESDtrack::kITSrefit*/);
	      dtc2etaphitpc[aniter]->SetOnFlyStatus(kFALSE); //kTRUE
	      dtc2etaphitpc[aniter]->SetParticleType(1);
	      dtc2etaphitpc[aniter]->SetMaxDcaV0Daughters(0.4); //1.5 Jai, 0.6
	      dtc2etaphitpc[aniter]->SetMaxDcaV0(1.0);
	      dtc2etaphitpc[aniter]->SetMinDaughtersToPrimVertex(0.1);
	      dtc2etaphitpc[aniter]->SetMaxCosPointingAngle(0.998); //0.99 - Jai
	      dtc2etaphitpc[aniter]->SetInvariantMassLambda(1.112683,1.118683);

	      //ESD first particle cut -> Proton 3, 5; AntiProton 4, 6, 7, 8
	      dtc3etaphitpc[aniter] = new AliFemtoESDTrackCut();
	      dtc3etaphitpc[aniter]->SetMostProbableProton();
	      dtc3etaphitpc[aniter]->SetMass(ProtonMass);
	      dtc3etaphitpc[aniter]->SetCharge(1.0);
	      dtc3etaphitpc[aniter]->SetPt(0.7,5.0);
	      dtc3etaphitpc[aniter]->SetEta(-0.8,0.8);

	      // Track quality cuts	 
	      dtc3etaphitpc[aniter]->SetStatus(AliESDtrack::kTPCin);
	      dtc3etaphitpc[aniter]->SetminTPCncls(80);
	      dtc3etaphitpc[aniter]->SetRemoveKinks(kTRUE);
	      dtc3etaphitpc[aniter]->SetLabel(kFALSE);
	      dtc3etaphitpc[aniter]->SetMaxTPCChiNdof(4.0);
	      //dtc3etaphitpc[aniter]->SetMaxImpactXY(0.2);
	      dtc3etaphitpc[aniter]->SetMaxImpactZ(2.0);
	      dtc3etaphitpc[aniter]->SetMaxImpactXYPtDep(0.018, 0.035, -1.01);

	      //ESD first particle cut -> Proton 3, 5; AntiProton 4, 6, 8, 9
	      dtc4etaphitpc[aniter] = new AliFemtoESDTrackCut();
	      dtc4etaphitpc[aniter]->SetMostProbableProton();
	      dtc4etaphitpc[aniter]->SetMass(ProtonMass);
	      dtc4etaphitpc[aniter]->SetCharge(-1.0);
	      dtc4etaphitpc[aniter]->SetPt(0.7,5.0);
	      dtc4etaphitpc[aniter]->SetEta(-0.8,0.8);

	      // Track quality cuts	 
	      dtc4etaphitpc[aniter]->SetStatus(AliESDtrack::kTPCin);
	      dtc4etaphitpc[aniter]->SetminTPCncls(80);
	      dtc4etaphitpc[aniter]->SetRemoveKinks(kTRUE);
	      dtc4etaphitpc[aniter]->SetLabel(kFALSE);
	      dtc4etaphitpc[aniter]->SetMaxTPCChiNdof(4.0);
	      //dtc4etaphitpc[aniter]->SetMaxImpactXY(0.2);
	      dtc4etaphitpc[aniter]->SetMaxImpactZ(2.0);
	      dtc4etaphitpc[aniter]->SetMaxImpactXYPtDep(0.018, 0.035, -1.01);

	 
	      //V0 monitor
	      cutPass1V0[aniter] = new AliFemtoCutMonitorV0(Form("cutPass1%stpcM%i", chrgs[ichg], imult));
	      cutFail1V0[aniter] = new AliFemtoCutMonitorV0(Form("cutFail1%stpcM%i", chrgs[ichg], imult));
	      //dtc1etaphitpc[aniter]->AddCutMonitor(cutPass1V0[aniter], cutFail1V0[aniter]);

	      cutPass2V0[aniter] = new AliFemtoCutMonitorV0(Form("cutPass2%stpcM%i", chrgs[ichg], imult));
	      cutFail2V0[aniter] = new AliFemtoCutMonitorV0(Form("cutFail2%stpcM%i", chrgs[ichg], imult));
	      //dtc2etaphitpc[aniter]->AddCutMonitor(cutPass2V0[aniter], cutFail2V0[aniter]);
 
	      cutPass1YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutPass1%stpcM%i", chrgs[ichg], imult), ProtonMass);
	      cutFail1YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutFail1%stpcM%i", chrgs[ichg], imult), ProtonMass);
	      //if(ichg==7 || ichg==8 || ichg==9)
	      //dtc3etaphitpc[aniter]->AddCutMonitor(cutPass1YPtetaphitpc[aniter], cutFail1YPtetaphitpc[aniter]);

	      cutPass2YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutPass2%stpcM%i", chrgs[ichg], imult), ProtonMass);
	      cutFail2YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutFail2%stpcM%i", chrgs[ichg], imult), ProtonMass);
	      //if(ichg==7 || ichg==8 || ichg==9)
	      //dtc4etaphitpc[aniter]->AddCutMonitor(cutPass2YPtetaphitpc[aniter], cutFail2YPtetaphitpc[aniter]);
	  
	      cutPass1PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cutPass1%stpcM%i", chrgs[ichg], imult),2);
	      cutFail1PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cutFail1%stpcM%i", chrgs[ichg], imult),2);
	      //dtc1etaphitpc[aniter]->AddCutMonitor(cutPass1PIDetaphitpc[aniter], cutFail1PIDetaphitpc[aniter]);

	      cutPass1PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cutPass1%stpcM%i", chrgs[ichg], imult),2);
	      cutFail1PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cutFail1%stpcM%i", chrgs[ichg], imult),2);
	      //if(ichg == 4 || ichg == 5 || ichg == 6 || ichg == 7 || ichg == 8 || ichg == 9)
	      //dtc3etaphitpc[aniter]->AddCutMonitor(cutPass1PIDetaphitpc[aniter], cutFail1PIDetaphitpc[aniter]);


	      sqp1cetaphitpc[aniter] = new AliFemtoV0PairCut();
	      sqp1cetaphitpc[aniter]->SetDataType(AliFemtoPairCut::kAOD);
	      sqp1cetaphitpc[aniter]->SetTPCEntranceSepMinimum(0.00001);
	      sqp1cetaphitpc[aniter]->SetTPCExitSepMinimum(-1.);
	      sqp1cetaphitpc[aniter]->SetMinAvgSeparation(0,3);
	      sqp1cetaphitpc[aniter]->SetMinAvgSeparation(1,0);
	      sqp1cetaphitpc[aniter]->SetMinAvgSeparation(2,0);
	      sqp1cetaphitpc[aniter]->SetMinAvgSeparation(3,3);

	      sqp2cetaphitpc[aniter] = new AliFemtoV0TrackPairCut(); //lambda-proton
	      sqp2cetaphitpc[aniter]->SetShareQualityMax(1.0); //between V0 daughter and track
	      sqp2cetaphitpc[aniter]->SetShareFractionMax(0.05);
	      sqp2cetaphitpc[aniter]->SetTPCOnly(kTRUE);
	      sqp2cetaphitpc[aniter]->SetDataType(AliFemtoPairCut::kAOD);
	      sqp2cetaphitpc[aniter]->SetTPCEntranceSepMinimum(0.00001);
	      sqp2cetaphitpc[aniter]->SetTPCExitSepMinimum(-1.);
	      sqp2cetaphitpc[aniter]->SetKstarCut(0.04,AliFemtoV0TrackPairCut::kLambda,AliFemtoV0TrackPairCut::kProton); //0 - lambda, 2 - proton
	      sqp2cetaphitpc[aniter]->SetMinAvgSeparation(0,11); //0 - track-pos, 1 - track-neg
	      sqp2cetaphitpc[aniter]->SetMinAvgSeparation(1,0);

	      sqp3cetaphitpc[aniter] = new AliFemtoV0TrackPairCut(); //antilambda-antiproton
	      sqp3cetaphitpc[aniter]->SetShareQualityMax(1.0); //between V0 daughter and track
	      sqp3cetaphitpc[aniter]->SetShareFractionMax(0.05);
	      sqp3cetaphitpc[aniter]->SetTPCOnly(kTRUE);
	      sqp3cetaphitpc[aniter]->SetDataType(AliFemtoPairCut::kAOD);
	      sqp3cetaphitpc[aniter]->SetTPCEntranceSepMinimum(0.00001);
	      sqp3cetaphitpc[aniter]->SetTPCExitSepMinimum(-1.);
	      sqp3cetaphitpc[aniter]->SetKstarCut(0.04,AliFemtoV0TrackPairCut::kAntiLambda,AliFemtoV0TrackPairCut::kAntiProton); //1 - antilambda, 3 - antiproton
	      sqp3cetaphitpc[aniter]->SetMinAvgSeparation(0,0); //0 - track-pos, 1 - track-neg
	      sqp3cetaphitpc[aniter]->SetMinAvgSeparation(1,11);

	      sqp4cetaphitpc[aniter] = new AliFemtoV0TrackPairCut(); //lambda-antiproton, antilambda-proton
	      sqp4cetaphitpc[aniter]->SetShareQualityMax(1.0); //between V0 daughter and track
	      sqp4cetaphitpc[aniter]->SetShareFractionMax(0.05);
	      sqp4cetaphitpc[aniter]->SetTPCOnly(kTRUE);
	      sqp4cetaphitpc[aniter]->SetDataType(AliFemtoPairCut::kAOD);
	      sqp4cetaphitpc[aniter]->SetTPCEntranceSepMinimum(0.00001);
	      sqp4cetaphitpc[aniter]->SetTPCExitSepMinimum(-1.);
	      //SetMinAvgSeparation w if'ach ponizej

	      //sqp3cetaphitpc[aniter] = new AliFemtoPairCutAntiGamma();
	      sqp5cetaphitpc[aniter] = new AliFemtoPairCutRadialDistance();
	      sqp5cetaphitpc[aniter]->SetMagneticFieldSign(-1.0); //-1->field1, 1->field3
	      sqp5cetaphitpc[aniter]->SetPhiStarDifferenceMinimum(0.012);
	      sqp5cetaphitpc[aniter]->SetEtaDifferenceMinimum(0.017);
	      sqp5cetaphitpc[aniter]->SetShareQualityMax(1.0);
	      sqp5cetaphitpc[aniter]->SetShareFractionMax(0.05);
	      sqp5cetaphitpc[aniter]->SetRemoveSameLabel(kFALSE);
	      sqp5cetaphitpc[aniter]->SetMaxEEMinv(0.002);
	      sqp5cetaphitpc[aniter]->SetMaxThetaDiff(0.008);
	      sqp5cetaphitpc[aniter]->SetDataType(AliFemtoPairCut::kAOD);
	      sqp5cetaphitpc[aniter]->SetTPCEntranceSepMinimum(0.0001);
          
	      avgsepcorr[aniter] = new AliFemtoAvgSepCorrFctn(Form("Avgsep%stpcM%i", chrgs[ichg], imult),5000,0,500);
	  

	      if(ichg == 0) //V0LL
		{
		  anetaphitpc[aniter]->SetEventCut(mecetaphitpc[aniter]);
		  anetaphitpc[aniter]->SetFirstParticleCut(dtc1etaphitpc[aniter]);
		  anetaphitpc[aniter]->SetSecondParticleCut(dtc1etaphitpc[aniter]);
		  anetaphitpc[aniter]->SetPairCut(sqp1cetaphitpc[aniter]);
		  avgsepcorr[aniter]->SetPairType(AliFemtoAvgSepCorrFctn::kV0s);
		}
	      else if(ichg == 1) //V0ALAL
		{
		  anetaphitpc[aniter]->SetEventCut(mecetaphitpc[aniter]);
		  anetaphitpc[aniter]->SetFirstParticleCut(dtc2etaphitpc[aniter]);
		  anetaphitpc[aniter]->SetSecondParticleCut(dtc2etaphitpc[aniter]);
		  anetaphitpc[aniter]->SetPairCut(sqp1cetaphitpc[aniter]);
		  avgsepcorr[aniter]->SetPairType(AliFemtoAvgSepCorrFctn::kV0s);
		}
	      else if(ichg == 2) //VOLAL
		{
		  anetaphitpc[aniter]->SetEventCut(mecetaphitpc[aniter]);
		  anetaphitpc[aniter]->SetFirstParticleCut(dtc1etaphitpc[aniter]);
		  anetaphitpc[aniter]->SetSecondParticleCut(dtc2etaphitpc[aniter]);
		  anetaphitpc[aniter]->SetPairCut(sqp1cetaphitpc[aniter]);
		  avgsepcorr[aniter]->SetPairType(AliFemtoAvgSepCorrFctn::kV0s);
		}
	      else if(ichg == 3) //V0PL
		{
		  anetaphitpc[aniter]->SetEventCut(mecetaphitpc[aniter]);
		  anetaphitpc[aniter]->SetFirstParticleCut(dtc1etaphitpc[aniter]);
		  anetaphitpc[aniter]->SetSecondParticleCut(dtc3etaphitpc[aniter]);
		  anetaphitpc[aniter]->SetPairCut(sqp2cetaphitpc[aniter]);
		  avgsepcorr[aniter]->SetPairType(AliFemtoAvgSepCorrFctn::kTrackV0);
		}
	      else if(ichg == 4) //V0APL
		{
		  sqp4cetaphitpc[aniter]->SetMinAvgSeparation(0,0); //0 - track-pos, 1 - track-neg
		  sqp4cetaphitpc[aniter]->SetMinAvgSeparation(1,11);
		  anetaphitpc[aniter]->SetEventCut(mecetaphitpc[aniter]);
		  anetaphitpc[aniter]->SetFirstParticleCut(dtc1etaphitpc[aniter]);
		  anetaphitpc[aniter]->SetSecondParticleCut(dtc4etaphitpc[aniter]);
		  anetaphitpc[aniter]->SetPairCut(sqp4cetaphitpc[aniter]);
		  avgsepcorr[aniter]->SetPairType(AliFemtoAvgSepCorrFctn::kTrackV0);
		}
	      else if(ichg == 5) //V0PAL
		{
		  sqp4cetaphitpc[aniter]->SetMinAvgSeparation(0,11); //0 - track-pos, 1 - track-neg
		  sqp4cetaphitpc[aniter]->SetMinAvgSeparation(1,0);
		  anetaphitpc[aniter]->SetEventCut(mecetaphitpc[aniter]);
		  anetaphitpc[aniter]->SetFirstParticleCut(dtc2etaphitpc[aniter]);
		  anetaphitpc[aniter]->SetSecondParticleCut(dtc3etaphitpc[aniter]);
		  anetaphitpc[aniter]->SetPairCut(sqp4cetaphitpc[aniter]);
		  avgsepcorr[aniter]->SetPairType(AliFemtoAvgSepCorrFctn::kTrackV0);
		}
	      else if(ichg == 6) //V0APAL
		{
		  anetaphitpc[aniter]->SetEventCut(mecetaphitpc[aniter]);
		  anetaphitpc[aniter]->SetFirstParticleCut(dtc2etaphitpc[aniter]);
		  anetaphitpc[aniter]->SetSecondParticleCut(dtc4etaphitpc[aniter]);
		  anetaphitpc[aniter]->SetPairCut(sqp3cetaphitpc[aniter]);
		  avgsepcorr[aniter]->SetPairType(AliFemtoAvgSepCorrFctn::kTrackV0);
		}
	      else if(ichg ==7) //PP
		{
		  anetaphitpc[aniter]->SetEventCut(mecetaphitpc[aniter]);
		  anetaphitpc[aniter]->SetFirstParticleCut(dtc3etaphitpc[aniter]);
		  anetaphitpc[aniter]->SetSecondParticleCut(dtc3etaphitpc[aniter]);
		  anetaphitpc[aniter]->SetPairCut(sqp5cetaphitpc[aniter]);
		  avgsepcorr[aniter]->SetPairType(AliFemtoAvgSepCorrFctn::kTracks);
		}
	      else if(ichg ==8) //PAP
		{
		  anetaphitpc[aniter]->SetEventCut(mecetaphitpc[aniter]);
		  anetaphitpc[aniter]->SetFirstParticleCut(dtc3etaphitpc[aniter]);
		  anetaphitpc[aniter]->SetSecondParticleCut(dtc4etaphitpc[aniter]);
		  anetaphitpc[aniter]->SetPairCut(sqp5cetaphitpc[aniter]);
		  avgsepcorr[aniter]->SetPairType(AliFemtoAvgSepCorrFctn::kTracks);
		}
	      else if(ichg ==9) //APAP
		{
		  anetaphitpc[aniter]->SetEventCut(mecetaphitpc[aniter]);
		  anetaphitpc[aniter]->SetFirstParticleCut(dtc4etaphitpc[aniter]);
		  anetaphitpc[aniter]->SetSecondParticleCut(dtc4etaphitpc[aniter]);
		  anetaphitpc[aniter]->SetPairCut(sqp5cetaphitpc[aniter]);
		  avgsepcorr[aniter]->SetPairType(AliFemtoAvgSepCorrFctn::kTracks);
		}
	  
	      anetaphitpc[aniter]->AddCorrFctn(avgsepcorr[aniter]);

	      //cylmetaphitpc[aniter] = new AliFemtoCorrFctnDirectYlm(Form("cylm%stpcM%i", chrgs[ichg], imult),3,nbinssh,0.0,shqmax,runshlcms);
	      //anetaphitpc[aniter]->AddCorrFctn(cylmetaphitpc[aniter]);
	  
	      //cqinvnclstpc[aniter] = new AliFemtoCorrFctnTPCNcls(Form("cqinvncls%stpcM%i", chrgs[ichg], imult),nbinssh,0.0,shqmax);
	      //anetaphitpc[aniter]->AddCorrFctn(cqinvnclstpc[aniter]);

	      //cdedpetaphi[aniter] = new AliFemtoCorrFctnDEtaDPhi(Form("cdedp%stpcM%i", chrgs[ichg], imult),45, 45);
	      //anetaphitpc[aniter]->AddCorrFctn(cdedpetaphi[aniter]);

          
	      //cqinvtitpc[aniter] = new AliFemtoTPCInnerCorrFctn(Form("cqinvti%stpcM%i", chrgs[ichg], imult),nbinssh,0.0,shqmax);
	      //anetaphitpc[aniter]->AddCorrFctn(cqinvtitpc[aniter]);
	  
	      //if(ichg==3 || ichg==4 || ichg==5 || ichg==6 || ichg==8) { //PL, APL, PAL, APAL, PAP
	      //cnonidtpc[aniter] = new AliFemtoCorrFctnNonIdDR(Form("cnonid%stpcM%i", chrgs[ichg], imult), nbinssh, 0.0,shqmax);
	      //anetaphitpc[aniter]->AddCorrFctn(cnonidtpc[aniter]);
	      //}
	      //else {
	      //cqinvtpc[aniter] = new AliFemtoQinvCorrFctn(Form("cqinv%stpcM%i", chrgs[ichg], imult),nbinssh,0.0,shqmax);
	      //anetaphitpc[aniter]->AddCorrFctn(cqinvtpc[aniter]);
	      //}
	  
	      cylmkttpc[aniter] = new AliFemtoCorrFctnDirectYlm(Form("cylm%stpcM%i", chrgs[ichg], imult),3, nbinssh, 0.0,shqmax, runshlcms);
	      anetaphitpc[aniter]->AddCorrFctn(cylmkttpc[aniter]);
	  
	      if (runktdep)
		{
		  int ktm;
		  for (int ikt=0; ikt<numOfkTBins; ikt++) {
		    ktm = aniter*numOfkTBins + ikt;
		    ktpcuts[ktm] = new AliFemtoKTPairCut(ktrng[ikt], ktrng[ikt+1]);
	      
		    cylmkttpc[ktm] = new AliFemtoCorrFctnDirectYlm(Form("cylm%stpcM%ikT%i", chrgs[ichg], imult, ikt),3,nbinssh, 0.0,shqmax,runshlcms);
		    cylmkttpc[ktm]->SetPairSelectionCut(ktpcuts[ktm]);
		    anetaphitpc[aniter]->AddCorrFctn(cylmkttpc[ktm]);

		    //cdedpetaphi[ktm] = new AliFemtoCorrFctnDEtaDPhi(Form("cdedp%stpcM%ikT%i", chrgs[ichg], imult,ikt),45, 45);
		    // cdedpetaphi[ktm]->SetPairSelectionCut(ktpcuts[ktm]);
		    //anetaphitpc[aniter]->AddCorrFctn(cdedpetaphi[ktm]);

		    //if(ichg==3 || ichg==4 || ichg==5 || ichg==6 || ichg==8) { //PL, APL, PAL, APAL, PAP
		    //cnonidtpc[ktm] = new AliFemtoCorrFctnNonIdDR(Form("cnonid%stpcM%ikT%i", chrgs[ichg], imult,ikt),nbinssh, 0.0,shqmax);
		    //cnonidtpc[ktm]->SetPairSelectionCut(ktpcuts[ktm]);     
		    //anetaphitpc[aniter]->AddCorrFctn(cnonidtpc[ktm]);
		    //}

		    //else{ 
		    //cqinvkttpc[ktm] = new AliFemtoQinvCorrFctn(Form("cqinv%stpcM%ikT%i", chrgs[ichg], imult, ikt),nbinssh,0.0,shqmax);
		    //cqinvkttpc[ktm]->SetPairSelectionCut(ktpcuts[ktm]);
		    //anetaphitpc[aniter]->AddCorrFctn(cqinvkttpc[ktm]);
             
		    //}
     
		  }
		}
	  
	      Manager->AddAnalysis(anetaphitpc[aniter]);	
 	}
      }
    }
  }
  // *** End of analysis

  return Manager;
}                         
                      
