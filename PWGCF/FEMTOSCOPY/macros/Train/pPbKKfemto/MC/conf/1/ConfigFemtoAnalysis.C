
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
#include "AliFemtoPairCutAntiGamma.h"
#include "AliFemtoPairCutRadialDistance.h"
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
#include "AliFemtoCutMonitorParticlePtPDG.h"
#include "AliFemtoKTPairCut.h"
#include "AliFemtoCutMonitorCollections.h"
#endif

//________________________________________________________________________
AliFemtoManager* ConfigFemtoAnalysis() {

  double PionMass = 0.13956995;
  double KaonMass = 0.493677;
	
  //multiplicity bins
  //int runmults[10] = {1, 1, 1, 1, 0, 0, 0, 0, 0, 0};
  //int multbins[11] = {0.01, 200, 400, 600, 1000, 950, 500, 600, 700, 800, 900};
  int runmults[3] = {1, 1, 1};
  int multbins[4] = {0.01, 200, 400, 900};

  int runch[2] = {1, 0};
  //const char *chrgs[2] = { "pip", "pim" };
  const char *chrgs[2] = { "Kp", "Km"};

  int runktdep = 1;
  //double ktrng[8] = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0};
  double ktrng[3] = {0.2, 0.5, 1.0};

  int run3d = 0; // Do 3D cartesian analysis?
  //int runshlcms = 1;
  int runshlcms = 0;

  //PhysicsSelection set only in runBatch with trigger + Physics Selection Task


  double shqmax;
  //int nbinssh = 200;
  int nbinssh = 200;

  //if (runshlcms) shqmax = 2.0;
  if (runshlcms) shqmax = 0.25;
  else shqmax = 2.0;

  AliFemtoEventReaderAODChain *Reader = new AliFemtoEventReaderAODChain();
  //Reader->SetUseMultiplicity(AliFemtoEventReaderAODChain::kCentrality);
  Reader->SetFilterBit(7);
  Reader->SetCentralityPreSelection(0.0001, 900);
  Reader->SetReadMC(kTRUE);
  Reader->SetDCAglobalTrack(kTRUE);

  AliFemtoModelGausLCMSFreezeOutGenerator *tFreeze = new AliFemtoModelGausLCMSFreezeOutGenerator();
  tFreeze->SetSizeOut(1.7*TMath::Sqrt(1.5));         
  tFreeze->SetSizeSide(1.7*TMath::Sqrt(1.5));
  tFreeze->SetSizeLong(1.7*TMath::Sqrt(1.5));

  AliFemtoModelWeightGeneratorLednicky *tWeight = new AliFemtoModelWeightGeneratorLednicky();
  tWeight->SetPairType(AliFemtoModelWeightGenerator::KaonPlusKaonMinus());
  tWeight->SetCoulOff();
  tWeight->SetQuantumOff();
  tWeight->SetStrongOff();
  tWeight->Set3BodyOff();

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
  AliFemtoCutMonitorCollections   *cutPassColletaphitpc[20];
  AliFemtoCutMonitorCollections   *cutFailColletaphitpc[20];
  //AliFemtoKKTrackCut           *dtc1etaphitpc[20];
  //AliFemtoKKTrackCut           *dtc2etaphitpc[20];
  AliFemtoKpm45TrackCut           *dtc1etaphitpc[20];
  AliFemtoKpm45TrackCut           *dtc2etaphitpc[20];
  //AliFemtoESDTrackCut           *dtc1etaphitpc[20];
  //AliFemtoESDTrackCut           *dtc2etaphitpc[20];
  AliFemtoCutMonitorParticleYPt *cutPass1YPtetaphitpc[20];
  AliFemtoCutMonitorParticleYPt *cutFail1YPtetaphitpc[20];
  AliFemtoCutMonitorParticleYPt *cutPass2YPtetaphitpc[20];
  AliFemtoCutMonitorParticleYPt *cutFail2YPtetaphitpc[20];
  AliFemtoCutMonitorParticlePID *cutPass1PIDetaphitpc[20];
  AliFemtoCutMonitorParticlePID *cutFail1PIDetaphitpc[20];
  AliFemtoCutMonitorParticlePID *cutPass2PIDetaphitpc[20];
  AliFemtoCutMonitorParticlePID *cutFail2PIDetaphitpc[20];
  AliFemtoPairCutAntiGamma      *sqpcetaphitpc[20];
  //AliFemtoShareQualityPairCut      *sqpcetaphitpc[20];
  //AliFemtoPairCutRadialDistance      *sqpcetaphitpc[20];
  //AliFemtoCorrFctnDirectYlm     *cylmetaphitpc[20];
  //AliFemtoCorrFctnDEtaDPhi      *cdedpetaphi[20];
  AliFemtoChi2CorrFctn          *cchiqinvetaphitpc[20];
  //AliFemtoKTPairCut             *ktpcuts[20*7];
  //AliFemtoCorrFctnDirectYlm     *cylmkttpc[20*7];
  //AliFemtoQinvCorrFctn          *cqinvkttpc[20*7];
  //AliFemtoCorrFctn3DLCMSSym     *cq3dlcmskttpc[20*7];
  AliFemtoKTPairCut             *ktpcuts[20*8];
  AliFemtoCorrFctnDirectYlm     *cylmkttpc[20*8];
  AliFemtoQinvCorrFctn          *cqinvkttpc[20*8];
  AliFemtoCorrFctn3DLCMSSym     *cq3dlcmskttpc[20*8];
  AliFemtoCorrFctnTPCNcls       *cqinvnclstpc[20];
  AliFemtoShareQualityCorrFctn  *cqinvsqtpc[20*10];
  AliFemtoChi2CorrFctn          *cqinvchi2tpc[20];
  AliFemtoTPCInnerCorrFctn      *cqinvinnertpc[20*10];
  AliFemtoCorrFctnGammaMonitor  *cgamma[20*10];
  AliFemtoModelCorrFctn   *cqinvkttpcmodel[20*8];

  // *** Begin pion-pion analysis ***
  int aniter = 0;

  for (int imult=0; imult<3; imult++) {
    if (runmults[imult]) {
      for (int ichg=0; ichg<1; ichg++) {
	if (runch[ichg]) {
	  aniter = ichg*3+imult;

	  anetaphitpc[aniter] = new AliFemtoVertexMultAnalysis(10, -8.0, 8.0, 4, multbins[imult], multbins[imult+1]);
	  anetaphitpc[aniter]->SetNumEventsToMix(5);
	  anetaphitpc[aniter]->SetMinSizePartCollection(1);

	  mecetaphitpc[aniter] = new AliFemtoBasicEventCut();
	  mecetaphitpc[aniter]->SetEventMult(0,10000);
	  mecetaphitpc[aniter]->SetVertZPos(-8.0,8.0);

	  cutPassEvMetaphitpc[aniter] = new AliFemtoCutMonitorEventMult(Form("cutPass%stpcM%i", chrgs[ichg], imult));
	  cutFailEvMetaphitpc[aniter] = new AliFemtoCutMonitorEventMult(Form("cutFail%stpcM%i", chrgs[ichg], imult));
	  mecetaphitpc[aniter]->AddCutMonitor(cutPassEvMetaphitpc[aniter], cutFailEvMetaphitpc[aniter]);
	  
	  cutPassEvVetaphitpc[aniter] = new AliFemtoCutMonitorEventVertex(Form("cutPass%stpcM%i", chrgs[ichg], imult));
	  cutFailEvVetaphitpc[aniter] = new AliFemtoCutMonitorEventVertex(Form("cutFail%stpcM%i", chrgs[ichg], imult));
	  mecetaphitpc[aniter]->AddCutMonitor(cutPassEvVetaphitpc[aniter], cutFailEvVetaphitpc[aniter]);
	  
	  cutPassColletaphitpc[aniter] = new AliFemtoCutMonitorCollections(Form("cutPass%stpcM%i", chrgs[ichg], imult));
          cutFailColletaphitpc[aniter] = new AliFemtoCutMonitorCollections(Form("cutFail%stpcM%i", chrgs[ichg], imult));
          mecetaphitpc[aniter]->AddCutMonitor(cutPassColletaphitpc[aniter], cutFailColletaphitpc[aniter]);

	  //dtc1etaphitpc[aniter] = new AliFemtoESDTrackCut();
	  //dtc1etaphitpc[aniter] = new AliFemtoKKTrackCut();
	  //dtc2etaphitpc[aniter] = new AliFemtoKKTrackCut();
	  dtc1etaphitpc[aniter] = new AliFemtoKpm45TrackCut();
	  dtc2etaphitpc[aniter] = new AliFemtoKpm45TrackCut();//K+-

	  if (ichg == 0){
	    dtc1etaphitpc[aniter]->SetCharge(1.0);
	    dtc2etaphitpc[aniter]->SetCharge(-1.0);
	  }
	  else if (ichg == 1){
	    dtc1etaphitpc[aniter]->SetCharge(-1.0);
	    dtc2etaphitpc[aniter]->SetCharge(1.0);
	  }

	  //dtc1etaphitpc[aniter]->SetPt(0.12,4.0);
	  //dtc1etaphitpc[aniter]->SetEta(-1.2,1.2);
	  dtc1etaphitpc[aniter]->SetPt(0.14,1.5);
	  dtc1etaphitpc[aniter]->SetEta(-0.8,0.8);
	  dtc2etaphitpc[aniter]->SetPt(0.14,1.5);
	  dtc2etaphitpc[aniter]->SetEta(-0.8,0.8);

	  //PID method
	  //dtc1etaphitpc[aniter]->SetMass(PionMass);
	  //dtc1etaphitpc[aniter]->SetMostProbablePion();
	  dtc1etaphitpc[aniter]->SetMass(KaonMass);
	  dtc1etaphitpc[aniter]->SetMostProbableKaon();
	  dtc2etaphitpc[aniter]->SetMass(KaonMass);
	  dtc2etaphitpc[aniter]->SetMostProbableKaon();
	  //dtc1etaphitpc[aniter]->SetPIDMethod(AliFemtoESDTrackCut::kContour);

//------------------- November 2013 -----------------------------------< 
	  // new cuts to remove electron (do not take into analysis if 400<p<500) 
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

	 dtc2etaphitpc[aniter]->SetNsigmaTPCle250(2.0);
	 dtc2etaphitpc[aniter]->SetNsigmaTPC250_400(2.0);
	 dtc2etaphitpc[aniter]->SetNsigmaTPC400_450(1.0);
	 dtc2etaphitpc[aniter]->SetNsigmaTPC450_500(2.0);
	 dtc2etaphitpc[aniter]->SetNsigmaTPCge500(3.0);    
	 // new cuts are stronger, better separation of pion in TOF 
	 // when momentum is greater then 800 MeV/c
	 dtc2etaphitpc[aniter]->SetNsigmaTOF500_800(2.0);
	 dtc2etaphitpc[aniter]->SetNsigmaTOF800_1000(1.5);
	 dtc2etaphitpc[aniter]->SetNsigmaTOFge1000(1.0);
	  //------------------- November 2013 -------------------------->

	  //Track quality cuts
	  /////dtc1etaphitpc[aniter]->SetStatus(AliESDtrack::kTPCrefit|AliESDtrack::kITSrefit);
	  //dtc1etaphitpc[aniter]->SetminTPCncls(80);
	  //dtc2etaphitpc[aniter]->SetminTPCncls(80);
	  dtc1etaphitpc[aniter]->SetRemoveKinks(kTRUE);
	  dtc2etaphitpc[aniter]->SetRemoveKinks(kTRUE);

	  //dtc1etaphitpc[aniter]->SetMaxTPCChiNdof(4.0);
	  /////dtc1etaphitpc[aniter]->SetMaxITSChiNdof(36);	  
	  dtc1etaphitpc[aniter]->SetLabel(kFALSE);
	  //dtc2etaphitpc[aniter]->SetMaxTPCChiNdof(4.0);
	  dtc2etaphitpc[aniter]->SetLabel(kFALSE);
	  
	  //primary particles: hits in ITS + DCA cut
	  //dtc1etaphitpc[aniter]->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
	  //				 AliESDtrackCuts::kAny);
	  //dtc1etaphitpc[aniter]->SetMaxImpactZ(3.0);
	  //dtc1etaphitpc[aniter]->SetMaxImpactXY(2.4);
	  //dtc2etaphitpc[aniter]->SetMaxImpactZ(3.0);
	  //dtc2etaphitpc[aniter]->SetMaxImpactXY(2.4);
	  /////dtc1etaphitpc[aniter]->SetMaxImpactXYPtDep(0.0105, 0.0350, -1.1);
	  //dtc1etaphitpc[aniter]->SetMaxImpactXYPtDep(0.0182, 0.0350, -1.01);
	  //dtc1etaphitpc[aniter]->SetMaxSigmaToVertex(6.0);

	  //cutPass1YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutPass1%stpcM%i", chrgs[ichg], imult), 0.13957);
	  //cutFail1YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutFail1%stpcM%i", chrgs[ichg], imult), 0.13957);
	  cutPass1YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutPass1%stpcM%i", chrgs[0], imult), 0.493677);
	  cutFail1YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutFail1%stpcM%i", chrgs[0], imult), 0.493677);
	  dtc1etaphitpc[aniter]->AddCutMonitor(cutPass1YPtetaphitpc[aniter], cutFail1YPtetaphitpc[aniter]);

	  cutPass2YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutPass2%stpcM%i", chrgs[1], imult), 0.493677);
	  cutFail2YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutFail2%stpcM%i", chrgs[1], imult), 0.493677);
	  dtc2etaphitpc[aniter]->AddCutMonitor(cutPass2YPtetaphitpc[aniter], cutFail2YPtetaphitpc[aniter]);

	  
	  //cutPass1PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cutPass1%stpcM%i", chrgs[ichg], imult),0);
	  //cutFail1PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cutFail1%stpcM%i", chrgs[ichg], imult),0);
	  cutPass1PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cutPass1%stpcM%i", chrgs[ichg], imult),1);
	  cutFail1PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cutFail1%stpcM%i", chrgs[ichg], imult),1);
	  dtc1etaphitpc[aniter]->AddCutMonitor(cutPass1PIDetaphitpc[aniter], cutFail1PIDetaphitpc[aniter]);


	  sqpcetaphitpc[aniter] = new AliFemtoPairCutAntiGamma();
	  //sqpcetaphitpc[aniter] = new AliFemtoShareQualityPairCut();
	  //sqpcetaphitpc[aniter] = new AliFemtoShareQualityTPCEntranceSepPairCut();
	  //sqpcetaphitpc[aniter] = new AliFemtoPairCutRadialDistance();
	  //sqpcetaphitpc[aniter]->SetShareQualityMax(1.0);
	  //sqpcetaphitpc[aniter]->SetShareFractionMax(0.05);
	  sqpcetaphitpc[aniter]->SetRemoveSameLabel(kFALSE);
	  //sqpcetaphitpc[aniter]->SetEtaDifferenceMinimum(0.016);
	  //sqpcetaphitpc[aniter]->SetPhiStarDifferenceMinimum(0.02);
	  //runtype==0

	  //sqpcetaphitpc[aniter]->SetMaxEEMinv(0.01);
	  //sqpcetaphitpc[aniter]->SetMaxThetaDiff(0.1);
	  //sqpcetaphitpc[aniter]->SetTPCEntranceSepMinimum(0.00001);

	  sqpcetaphitpc[aniter]->SetMaxEEMinv(0.0);
	  sqpcetaphitpc[aniter]->SetMaxThetaDiff(0.0);
	  sqpcetaphitpc[aniter]->SetTPCEntranceSepMinimum(0.0);

	  //sqpcetaphitpc[aniter]->SetTPCEntranceSepMinimum(0.00001);

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
	  anetaphitpc[aniter]->SetFirstParticleCut(dtc1etaphitpc[aniter]);
	  anetaphitpc[aniter]->SetSecondParticleCut(dtc2etaphitpc[aniter]);
	  anetaphitpc[aniter]->SetPairCut(sqpcetaphitpc[aniter]);
	  
	  //Correlation functions

	  //Qinv (without kT bins)
	  //cqinvkttpc[aniter] = new AliFemtoQinvCorrFctn(Form("cqinv%stpcM%i", chrgs[ichg], imult),nbinssh,0.0,shqmax);
	  //anetaphitpc[aniter]->AddCorrFctn(cqinvkttpc[aniter]);

	  if (runktdep) {
	    int ktm;
	    for (int ikt=0; ikt<2; ikt++) {
	      ktm = aniter*2 + ikt;
	      ktpcuts[ktm] = new AliFemtoKTPairCut(ktrng[ikt], ktrng[ikt+1]);
	      
	      cqinvkttpc[ktm] = new AliFemtoQinvCorrFctn(Form("cqinv%stpcM%ikT%i", chrgs[ichg], imult, ikt),nbinssh,0.0, shqmax);
	      cqinvkttpc[ktm]->SetPairSelectionCut(ktpcuts[ktm]);
	      anetaphitpc[aniter]->AddCorrFctn(cqinvkttpc[ktm]);

	      // model--------------
	      cqinvkttpcmodel[ktm] = new AliFemtoModelCorrFctn(Form("cqinvModel%stpcM%ikT%i", chrgs[ichg], imult, ikt),nbinssh,0.0,(imult>6)?shqmax*2.5:shqmax);
	      cqinvkttpcmodel[ktm]->SetPairSelectionCut(ktpcuts[ktm]);//add kT bins
              cqinvkttpcmodel[ktm]->ConnectToManager(tModelManager);
	      anetaphitpc[aniter]->AddCorrFctn(cqinvkttpcmodel[ktm]);// add CF histos

	      cqinvsqtpc[ktm] = new AliFemtoShareQualityCorrFctn(Form("cqinvsq%stpcM%ikT%i", chrgs[ichg], imult, ikt),nbinssh,0.0,shqmax);
	      cqinvsqtpc[ktm]->SetPairSelectionCut(ktpcuts[ktm]);
	      anetaphitpc[aniter]->AddCorrFctn(cqinvsqtpc[ktm]);

	      cqinvinnertpc[ktm] = new AliFemtoTPCInnerCorrFctn(Form("cqinvinner%stpcM%ikT%i", chrgs[ichg], imult, ikt),nbinssh,0.0,shqmax);
	      cqinvinnertpc[ktm]->SetPairSelectionCut(ktpcuts[ktm]);
	      cqinvinnertpc[ktm]->SetRadius(1.6);
	      anetaphitpc[aniter]->AddCorrFctn(cqinvinnertpc[ktm]);
         
              cgamma[aniter] = new AliFemtoCorrFctnGammaMonitor(Form("cgammaM%ikT%i", imult, ikt),200,200);
              anetaphitpc[aniter]->AddCorrFctn(cgamma[aniter]);

	    }
	  }
	  
	  //cdedpetaphi[aniter] = new AliFemtoCorrFctnDEtaDPhi(Form("cdedp%stpcM%i", chrgs[ichg], imult),39, 39);
	  //anetaphitpc[aniter]->AddCorrFctn(cdedpetaphi[aniter]);
	  
	  Manager->AddAnalysis(anetaphitpc[aniter]);	
	}
      }
    }
  }
  // *** End K+K- analysis

  return Manager;
}
