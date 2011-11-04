
/*********************************************************************
 *                                                                   *
 * Configfemtoanalysis.C - configuration macro for the femtoscopic   *
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
#include "AliFemtoCorrFctnNonIdDR.h"
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
#endif

//________________________________________________________________________
AliFemtoManager* ConfigFemtoAnalysis() {

  double PionMass = 0.13956995;
  double KaonMass = 0.493677;
  double ProtonMass = 0.938272013;
	

  int runmults[10] = {1, 1, 1, 1, 1, 1, 1, 1, 0, 0};
  int multbins[11] = {0.001, 50, 100, 200, 300, 400, 500, 600, 700, 800, 900};
  
  int runch[3] = {1, 0, 0};
  const char *chrgs[3] = { "PP", "APAP", "PAP" };
  

  int runktdep = 1;
  double ktrng[3] = {0.01, 0.7, 100.0};
  //  double ktrng[8] = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0};

  int numOfMultBins = 10;  
  int numOfChTypes = 3;
  int numOfkTbins = 2;

  int runqinv = 1;
  int runshlcms = 1;// 0:PRF(PAP), 1:LCMS(PP,APAP)

  int runtype = 2; // Types 0 - global, 1 - ITS only, 2 - TPC Inner
  int isrealdata = 1;

  int gammacut = 1;
  
  double shqmax = 0.5;
  int nbinssh = 100;

  AliFemtoEventReaderAODChain *Reader = new AliFemtoEventReaderAODChain();
  Reader->SetFilterBit(7);
  Reader->SetCentralityPreSelection(0.001, 710);

  AliFemtoManager* Manager = new AliFemtoManager();
  Manager->SetEventReader(Reader);

  AliFemtoVertexMultAnalysis    *anetaphitpc[10*3];
  AliFemtoBasicEventCut         *mecetaphitpc[10*3];
  AliFemtoCutMonitorEventMult   *cutPassEvMetaphitpc[10*3];
  AliFemtoCutMonitorEventMult   *cutFailEvMetaphitpc[10*3];
  AliFemtoCutMonitorEventVertex *cutPassEvVetaphitpc[10*3];
  AliFemtoCutMonitorEventVertex *cutFailEvVetaphitpc[10*3];
  AliFemtoESDTrackCut           *dtc1etaphitpc[10*3];
  AliFemtoESDTrackCut           *dtc2etaphitpc[10*3];
  AliFemtoCutMonitorParticleYPt *cutPass1YPtetaphitpc[10*3];
  AliFemtoCutMonitorParticleYPt *cutFail1YPtetaphitpc[10*3];
  AliFemtoCutMonitorParticlePID *cutPass1PIDetaphitpc[10*3];
  AliFemtoCutMonitorParticlePID *cutFail1PIDetaphitpc[10*3];
  AliFemtoCutMonitorParticleYPt *cutPass2YPtetaphitpc[10*3];
  AliFemtoCutMonitorParticleYPt *cutFail2YPtetaphitpc[10*3];
  AliFemtoCutMonitorParticlePID *cutPass2PIDetaphitpc[10*3];
  AliFemtoCutMonitorParticlePID *cutFail2PIDetaphitpc[10*3];
  //   AliFemtoPairCutAntiGamma      *sqpcetaphitpcdiff[10*3];
  //   AliFemtoShareQualityTPCEntranceSepPairCut      *sqpcetaphitpcsame[10*3];
  AliFemtoPairCutAntiGamma      *sqpcetaphitpc[10*3];
  //  AliFemtoPairCutRadialDistance      *sqpcetaphitpc[10*3];
  //  AliFemtoChi2CorrFctn          *cchiqinvetaphitpc[20*2];
  AliFemtoKTPairCut             *ktpcuts[10*3*2];
  AliFemtoCorrFctnDirectYlm     *cylmkttpc[10*3];
  //AliFemtoCorrFctnDirectYlm     *cylmetaphitpc[10*3];
  AliFemtoQinvCorrFctn          *cqinvkttpc[10*3*2];
  AliFemtoQinvCorrFctn          *cqinvtpc[10*3];
  AliFemtoCorrFctnNonIdDR       *ckstartpc[10*3];
  AliFemtoCorrFctnDEtaDPhi      *cdedpetaphi[10*3];

  //   AliFemtoCorrFctn3DLCMSSym     *cq3dlcmskttpc[20*2];
  //   AliFemtoCorrFctnTPCNcls       *cqinvnclstpc[20];
  //   AliFemtoShareQualityCorrFctn  *cqinvsqtpc[20*10];
  //   AliFemtoChi2CorrFctn          *cqinvchi2tpc[20];
  //   AliFemtoTPCInnerCorrFctn      *cqinvinnertpc[20*10];
  
  // *** Third QA task - HBT analysis with all pair cuts off, TPC only ***
  // *** Begin pion-pion (positive) analysis ***
  int aniter = 0;  

  for (int imult = 0; imult < numOfMultBins; imult++) {
    if (runmults[imult]) {

      for (int ichg = 0; ichg < numOfChTypes; ichg++) {
	if (runch[ichg]) {

	  aniter = ichg * numOfMultBins + imult;
	  anetaphitpc[aniter] = new AliFemtoVertexMultAnalysis(8, -8.0, 8.0, 4, multbins[imult], multbins[imult+1]);
	  anetaphitpc[aniter]->SetNumEventsToMix(10);
	  anetaphitpc[aniter]->SetMinSizePartCollection(1);

	  mecetaphitpc[aniter] = new AliFemtoBasicEventCut();
	  mecetaphitpc[aniter]->SetEventMult(0.001,100000);
	  mecetaphitpc[aniter]->SetVertZPos(-8,8);

	  if (isrealdata)
	    mecetaphitpc[aniter]->SetAcceptOnlyPhysics(kTRUE);
	  	  
	  cutPassEvMetaphitpc[aniter] = new AliFemtoCutMonitorEventMult(Form("cutPass%stpcM%i", chrgs[ichg], imult));
	  cutFailEvMetaphitpc[aniter] = new AliFemtoCutMonitorEventMult(Form("cutFail%stpcM%i", chrgs[ichg], imult));
	  mecetaphitpc[aniter]->AddCutMonitor(cutPassEvMetaphitpc[aniter], cutFailEvMetaphitpc[aniter]);
	  
	  cutPassEvVetaphitpc[aniter] = new AliFemtoCutMonitorEventVertex(Form("cutPass%stpcM%i", chrgs[ichg], imult));
	  cutFailEvVetaphitpc[aniter] = new AliFemtoCutMonitorEventVertex(Form("cutFail%stpcM%i", chrgs[ichg], imult));
	  mecetaphitpc[aniter]->AddCutMonitor(cutPassEvVetaphitpc[aniter], cutFailEvVetaphitpc[aniter]);
	  
	  dtc1etaphitpc[aniter] = new AliFemtoESDTrackCut();
	  dtc2etaphitpc[aniter] = new AliFemtoESDTrackCut();

	  if (ichg == 0)
	    dtc1etaphitpc[aniter]->SetCharge(1.0);
	  else if (ichg == 1)
	    dtc1etaphitpc[aniter]->SetCharge(-1.0);
	  else if (ichg == 2) {
    	    dtc1etaphitpc[aniter]->SetCharge(-1.0);
	    dtc2etaphitpc[aniter]->SetCharge(1.0);
	  }
	  
	  dtc1etaphitpc[aniter]->SetEta(-0.8,0.8);
	  dtc1etaphitpc[aniter]->SetMass(ProtonMass);	  
	  dtc1etaphitpc[aniter]->SetMostProbableProton();

	  if (ichg == 2) {	    
	    dtc2etaphitpc[aniter]->SetEta(-0.8,0.8);
	    dtc2etaphitpc[aniter]->SetMass(ProtonMass);	  
	    dtc2etaphitpc[aniter]->SetMostProbableProton();
	  }
	  
	  // Track quality cuts

	  if (runtype == 0) {
	    dtc1etaphitpc[aniter]->SetStatus(AliESDtrack::kTPCrefit|AliESDtrack::kITSrefit);
	    //	    dtc1etaphitpc[aniter]->SetStatus(AliESDtrack::kTPCrefit);
	    //    dtc1etaphitpc[aniter]->SetStatus(AliESDtrack::kITSrefit);
	    dtc1etaphitpc[aniter]->SetminTPCncls(80);
	    dtc1etaphitpc[aniter]->SetRemoveKinks(kTRUE);
	    dtc1etaphitpc[aniter]->SetLabel(kFALSE);
	    //    dtc1etaphitpc[aniter]->SetMaxITSChiNdof(6.0);
	    dtc1etaphitpc[aniter]->SetMaxTPCChiNdof(4.0);
	    dtc1etaphitpc[aniter]->SetMaxImpactXY(0.2);
	    //            dtc1etaphitpc[aniter]->SetMaxImpactXYPtDep(0.0182, 0.0350, -1.01);
	    dtc1etaphitpc[aniter]->SetMaxImpactZ(0.15);
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
	    dtc1etaphitpc[aniter]->SetMaxImpactXY(0.2);
	    dtc1etaphitpc[aniter]->SetMaxImpactZ(0.25);
	    //      dtc1etaphitpc[aniter]->SetMaxSigmaToVertex(6.0);
	  }
	  else if (runtype == 2) {
	    dtc1etaphitpc[aniter]->SetStatus(AliESDtrack::kTPCin);
	    dtc1etaphitpc[aniter]->SetminTPCncls(80);
	    dtc1etaphitpc[aniter]->SetRemoveKinks(kTRUE);
	    dtc1etaphitpc[aniter]->SetLabel(kFALSE);
	    dtc1etaphitpc[aniter]->SetMaxTPCChiNdof(4.0);
	    dtc1etaphitpc[aniter]->SetMaxImpactXY(2.4);
	    dtc1etaphitpc[aniter]->SetMaxImpactZ(3.0);

	    if (ichg == 2) {
	      dtc2etaphitpc[aniter]->SetStatus(AliESDtrack::kTPCin);
	      dtc2etaphitpc[aniter]->SetminTPCncls(80);
	      dtc2etaphitpc[aniter]->SetRemoveKinks(kTRUE);
	      dtc2etaphitpc[aniter]->SetLabel(kFALSE);
	      dtc2etaphitpc[aniter]->SetMaxTPCChiNdof(4.0);
	      dtc2etaphitpc[aniter]->SetMaxImpactXY(2.4);
	      dtc2etaphitpc[aniter]->SetMaxImpactZ(3.0);
	    }
	    
	  }
	  
	  cutPass1YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutPass1%stpcM%i", chrgs[ichg], imult),ProtonMass);
	  cutFail1YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutFail1%stpcM%i", chrgs[ichg], imult),ProtonMass);
	  dtc1etaphitpc[aniter]->AddCutMonitor(cutPass1YPtetaphitpc[aniter], cutFail1YPtetaphitpc[aniter]);
	  
	  cutPass1PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cutPass1%stpcM%i", chrgs[ichg], imult),2);//0-pion,1-kaon,2-proton
	  cutFail1PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cutFail1%stpcM%i", chrgs[ichg], imult),2);
	  dtc1etaphitpc[aniter]->AddCutMonitor(cutPass1PIDetaphitpc[aniter], cutFail1PIDetaphitpc[aniter]);

	  if (ichg == 2){
	    cutPass2PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cutPass2%stpcM%i", chrgs[ichg], imult),2);//0-pion,1-kaon,2-proton
	    cutFail2PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cutFail2%stpcM%i", chrgs[ichg], imult),2);
	    dtc2etaphitpc[aniter]->AddCutMonitor(cutPass2PIDetaphitpc[aniter], cutFail2PIDetaphitpc[aniter]);
	  }

	  sqpcetaphitpc[aniter] = new AliFemtoPairCutAntiGamma();

	  if (runtype == 0) {
	    sqpcetaphitpc[aniter]->SetShareQualityMax(1.0);
	    sqpcetaphitpc[aniter]->SetShareFractionMax(0.05);
	    sqpcetaphitpc[aniter]->SetRemoveSameLabel(kFALSE);
	    // sqpcetaphitpc[aniter]->SetMaxEEMinv(0.0);
	    // sqpcetaphitpc[aniter]->SetMaxThetaDiff(0.0);
	    //	    sqpcetaphitpc[aniter]->SetTPCEntranceSepMinimum(1.5);
	    //sqpcetaphitpc[aniter]->SetRadialDistanceMinimum(0.12, 0.03);
	    //	    sqpcetaphitpc[aniter]->SetEtaDifferenceMinimum(0.02);
	  }
	  else if (runtype == 1) {
	    sqpcetaphitpc[aniter]->SetShareQualityMax(1.0);
	    sqpcetaphitpc[aniter]->SetShareFractionMax(1.05);
	    sqpcetaphitpc[aniter]->SetRemoveSameLabel(kFALSE);
	    // sqpcetaphitpc[aniter]->SetMaxEEMinv(0.002);
	    // sqpcetaphitpc[aniter]->SetMaxThetaDiff(0.008);
	    //	    sqpcetaphitpc[aniter]->SetTPCEntranceSepMinimum(5.0);
	    //sqpcetaphitpc[aniter]->SetRadialDistanceMinimum(1.2, 0.03);
	    //	    sqpcetaphitpc[aniter]->SetEtaDifferenceMinimum(0.02);
	  }
	  else if (runtype == 2) {
	    sqpcetaphitpc[aniter]->SetUseAOD(kTRUE);
	    sqpcetaphitpc[aniter]->SetShareQualityMax(1.0);
	    sqpcetaphitpc[aniter]->SetShareFractionMax(0.05);
	    sqpcetaphitpc[aniter]->SetRemoveSameLabel(kFALSE);

	    if (gammacut == 0) {	      
	      sqpcetaphitpc[aniter]->SetMaxEEMinv(0.0);
	      sqpcetaphitpc[aniter]->SetMaxThetaDiff(0.0);
	    }
	    else if (gammacut == 1) { 
	      sqpcetaphitpc[aniter]->SetMaxEEMinv(0.002);
	      sqpcetaphitpc[aniter]->SetMaxThetaDiff(0.008);
	    }

	    //phi-star cut - values from Johana
	    // 	    sqpcetaphitpc[aniter]->SetPhiStarDifferenceMinimum(0.012);
	    // 	    sqpcetaphitpc[aniter]->SetEtaDifferenceMinimum(0.017);

	  }
	  
	  anetaphitpc[aniter]->SetEventCut(mecetaphitpc[aniter]);

	  if (ichg == 2) {
	    anetaphitpc[aniter]->SetFirstParticleCut(dtc1etaphitpc[aniter]);
	    anetaphitpc[aniter]->SetSecondParticleCut(dtc2etaphitpc[aniter]);
	  }
	  else {
	    anetaphitpc[aniter]->SetFirstParticleCut(dtc1etaphitpc[aniter]);
	    anetaphitpc[aniter]->SetSecondParticleCut(dtc1etaphitpc[aniter]);
	  }
	    
	  anetaphitpc[aniter]->SetPairCut(sqpcetaphitpc[aniter]);
	  

	  if (ichg == 2) {
	    ckstartpc[aniter] = new AliFemtoCorrFctnNonIdDR(Form("ckstar%stpcM%i", chrgs[ichg], imult),nbinssh,0.0,shqmax);
	    anetaphitpc[aniter]->AddCorrFctn(ckstartpc[aniter]);
	  }
	  else {  
	    cqinvtpc[aniter] = new AliFemtoQinvCorrFctn(Form("cqinv%stpcM%i", chrgs[ichg], imult),nbinssh,0.0,shqmax);
	    anetaphitpc[aniter]->AddCorrFctn(cqinvtpc[aniter]);
	  }

          cylmkttpc[aniter] = new AliFemtoCorrFctnDirectYlm(Form("cylm%stpcM%i", chrgs[ichg], imult),2,nbinssh, 0.0,shqmax,runshlcms);
          anetaphitpc[aniter]->AddCorrFctn(cylmkttpc[aniter]);

	  
	  if (runktdep) {
	    int ktm;
	    for (int ikt=0; ikt<numOfkTbins; ikt++) {

	      ktm = aniter * numOfkTbins + ikt;
	      ktpcuts[ktm] = new AliFemtoKTPairCut(ktrng[ikt], ktrng[ikt+1]);
	      
	      //	      cylmkttpc[ktm] = new AliFemtoCorrFctnDirectYlm(Form("cylm%stpcM%ikT%i", chrgs[ichg], imult, ikt),3,
	      //							     nbinssh, 0.0,
	      //							     (imult>6)?shqmax*2.5:shqmax,
	      //							     runshlcms);
	      //	      cylmkttpc[ktm]->SetPairSelectionCut(ktpcuts[ktm]);
	      //	      anetaphitpc[aniter]->AddCorrFctn(cylmkttpc[ktm]);
	      
	      cqinvkttpc[ktm] = new AliFemtoQinvCorrFctn(Form("cqinv%stpcM%ikT%i", chrgs[ichg], imult, ikt),nbinssh,0.0,(imult>6)?shqmax*2.5:shqmax);
	      cqinvkttpc[ktm]->SetPairSelectionCut(ktpcuts[ktm]);
	      anetaphitpc[aniter]->AddCorrFctn(cqinvkttpc[ktm]);

	      // 	      cqinvsqtpc[ktm] = new AliFemtoShareQualityCorrFctn(Form("cqinvsq%stpcM%ikT%i", chrgs[ichg], imult, ikt),nbinssh,0.0,shqmax);
	      // 	      cqinvsqtpc[ktm]->SetPairSelectionCut(ktpcuts[ktm]);
	      // 	      anetaphitpc[aniter]->AddCorrFctn(cqinvsqtpc[ktm]);

	      // 	      cqinvinnertpc[ktm] = new AliFemtoTPCInnerCorrFctn(Form("cqinvinner%stpcM%ikT%i", chrgs[ichg], imult, ikt),nbinssh,0.0,shqmax);
	      // 	      cqinvinnertpc[ktm]->SetPairSelectionCut(ktpcuts[ktm]);
	      // 	      cqinvinnertpc[ktm]->SetRadius(1.2);
	      // 	      anetaphitpc[aniter]->AddCorrFctn(cqinvinnertpc[ktm]);

	      // 	      if (run3d) {
	      // 		cq3dlcmskttpc[ktm] = new AliFemtoCorrFctn3DLCMSSym(Form("cq3d%stpcM%ikT%i", chrgs[ichg], imult, ikt),60,(imult>3)?((imult>6)?((imult>7)?0.6:0.4):0.25):0.15);
	      // 		cq3dlcmskttpc[ktm]->SetPairSelectionCut(ktpcuts[ktm]);
	      // 		anetaphitpc[aniter]->AddCorrFctn(cq3dlcmskttpc[ktm]);
	      // 	      }
	    }
	  }
	  
	 // cdedpetaphi[aniter] = new AliFemtoCorrFctnDEtaDPhi(Form("cdedp%stpcM%i", chrgs[ichg], imult),24, 24);
	 // anetaphitpc[aniter]->AddCorrFctn(cdedpetaphi[aniter]);
	  
	  Manager->AddAnalysis(anetaphitpc[aniter]);	
	}
      }
    }
  }
  // *** End pion-pion (positive) analysis

  return Manager;
}                         
