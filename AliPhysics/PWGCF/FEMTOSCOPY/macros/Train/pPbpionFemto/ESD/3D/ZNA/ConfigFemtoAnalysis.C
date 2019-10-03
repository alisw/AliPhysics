
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
  int runmults[10] = {1, 1, 1, 1, 0, 0, 0, 0, 0, 0};
  int multbins[11] = {0.01, 200, 400, 600, 900, 950, 500, 600, 700, 800, 900};

  int runch[2] = {1, 1};
  const char *chrgs[2] = { "pip", "pim" };
  
  int runktdep = 1;
  double ktrng[8] = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0};

  int run3d = 1; // Do 3D cartesian analysis?
  int runshlcms = 1;

  //PhysicsSelection set only in runBatch with trigger + Physics Selection Task


  double shqmax;
  int nbinssh = 200;

  if (runshlcms) shqmax = 2.0;
  else shqmax = 0.9;

  //AliFemtoEventReaderAODChain *Reader = new AliFemtoEventReaderAODChain();
  //Reader->SetFilterBit(4);
  //Reader->SetCentralityPreSelection(500, 950);

  //AliFemtoEventReaderESDChainKine* Reader=new AliFemtoEventReaderESDChainKine();
  //Reader->SetConstrained(true);
  // Reader->SetUseTPCOnly(false);

  AliFemtoEventReaderESDChain *Reader = new AliFemtoEventReaderESDChain();
  Reader->SetUseMultiplicity(AliFemtoEventReaderESDChain::kCentralityZNA);
  Reader->SetReadTrackType(AliFemtoEventReaderESDChain::kGlobal);


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
  AliFemtoESDTrackCut           *dtc1etaphitpc[20];
  AliFemtoESDTrackCut           *dtc2etaphitpc[20];
  AliFemtoCutMonitorParticleYPt *cutPass1YPtetaphitpc[20];
  AliFemtoCutMonitorParticleYPt *cutFail1YPtetaphitpc[20];
  AliFemtoCutMonitorParticlePID *cutPass1PIDetaphitpc[20];
  AliFemtoCutMonitorParticlePID *cutFail1PIDetaphitpc[20];
  AliFemtoCutMonitorParticleYPt *cutPass2YPtetaphitpc[20];
  AliFemtoCutMonitorParticleYPt *cutFail2YPtetaphitpc[20];
  AliFemtoCutMonitorParticlePID *cutPass2PIDetaphitpc[20];
  AliFemtoCutMonitorParticlePID *cutFail2PIDetaphitpc[20];
  //  AliFemtoPairCutAntiGamma      *sqpcetaphitpc[20];
  AliFemtoShareQualityTPCEntranceSepPairCut      *sqpcetaphitpc[20];
  //AliFemtoPairCutRadialDistance      *sqpcetaphitpc[20];
  AliFemtoCorrFctnDirectYlm     *cylmetaphitpc[20];
  AliFemtoCorrFctnDEtaDPhi      *cdedpetaphi[20];
  AliFemtoChi2CorrFctn          *cchiqinvetaphitpc[20];
  AliFemtoKTPairCut             *ktpcuts[20*7];
  AliFemtoCorrFctnDirectYlm     *cylmkttpc[20*7];
  AliFemtoQinvCorrFctn          *cqinvkttpc[20*7];
  AliFemtoCorrFctn3DLCMSSym     *cq3dlcmskttpc[20*7];
  AliFemtoCorrFctnTPCNcls       *cqinvnclstpc[20];
  AliFemtoShareQualityCorrFctn  *cqinvsqtpc[20*10];
  AliFemtoChi2CorrFctn          *cqinvchi2tpc[20];
  AliFemtoTPCInnerCorrFctn      *cqinvinnertpc[20*10];

  // *** Begin pion-pion analysis ***
  int aniter = 0;

  for (int imult=0; imult<10; imult++) {
    if (runmults[imult]) {
      for (int ichg=0; ichg<2; ichg++) {
	if (runch[ichg]) {
	  aniter = ichg*10+imult;

	  anetaphitpc[aniter] = new AliFemtoVertexMultAnalysis(10, -10.0, 10.0, 4, multbins[imult], multbins[imult+1]);
	  anetaphitpc[aniter]->SetNumEventsToMix(5);
	  anetaphitpc[aniter]->SetMinSizePartCollection(1);

	  mecetaphitpc[aniter] = new AliFemtoBasicEventCut();
	  mecetaphitpc[aniter]->SetEventMult(0,10000);
	  mecetaphitpc[aniter]->SetVertZPos(-10,10);

	  cutPassEvMetaphitpc[aniter] = new AliFemtoCutMonitorEventMult(Form("cutPass%stpcM%i", chrgs[ichg], imult),500);
	  cutFailEvMetaphitpc[aniter] = new AliFemtoCutMonitorEventMult(Form("cutFail%stpcM%i", chrgs[ichg], imult),500);
	  mecetaphitpc[aniter]->AddCutMonitor(cutPassEvMetaphitpc[aniter], cutFailEvMetaphitpc[aniter]);
	  
	  cutPassEvVetaphitpc[aniter] = new AliFemtoCutMonitorEventVertex(Form("cutPass%stpcM%i", chrgs[ichg], imult));
	  cutFailEvVetaphitpc[aniter] = new AliFemtoCutMonitorEventVertex(Form("cutFail%stpcM%i", chrgs[ichg], imult));
	  mecetaphitpc[aniter]->AddCutMonitor(cutPassEvVetaphitpc[aniter], cutFailEvVetaphitpc[aniter]);
	  
	  cutPassColletaphitpc[aniter] = new AliFemtoCutMonitorCollections(Form("cutPass%stpcM%i", chrgs[ichg], imult));
          cutFailColletaphitpc[aniter] = new AliFemtoCutMonitorCollections(Form("cutFail%stpcM%i", chrgs[ichg], imult));
          mecetaphitpc[aniter]->AddCutMonitor(cutPassColletaphitpc[aniter], cutFailColletaphitpc[aniter]);

	  dtc1etaphitpc[aniter] = new AliFemtoESDTrackCut();

	  if (ichg == 0)
	    dtc1etaphitpc[aniter]->SetCharge(1.0);
	  else if (ichg == 1)
	    dtc1etaphitpc[aniter]->SetCharge(-1.0);

	  dtc1etaphitpc[aniter]->SetPt(0.12,4.0);
	  dtc1etaphitpc[aniter]->SetEta(-1.2,1.2);
	  
	  //PID method
	  dtc1etaphitpc[aniter]->SetMass(PionMass);
	  dtc1etaphitpc[aniter]->SetMostProbablePion();
	  //dtc1etaphitpc[aniter]->SetPIDMethod(AliFemtoESDTrackCut::kContour);

	  //Track quality cuts
	  dtc1etaphitpc[aniter]->SetStatus(AliESDtrack::kTPCrefit|AliESDtrack::kITSrefit);
	  dtc1etaphitpc[aniter]->SetminTPCncls(50);
	  dtc1etaphitpc[aniter]->SetRemoveKinks(kTRUE);


	  dtc1etaphitpc[aniter]->SetMaxTPCChiNdof(4.0);
	  dtc1etaphitpc[aniter]->SetMaxITSChiNdof(36);	  
	  dtc1etaphitpc[aniter]->SetLabel(kFALSE);
	  
	  //primary particles: hits in ITS + DCA cut
	  dtc1etaphitpc[aniter]->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
	  				 AliESDtrackCuts::kAny);
	  dtc1etaphitpc[aniter]->SetMaxImpactZ(2.0);
	  //dtc1etaphitpc[aniter]->SetMaxImpactXY(2.4);
	  dtc1etaphitpc[aniter]->SetMaxImpactXYPtDep(0.0105, 0.0350, -1.1);
	  //dtc1etaphitpc[aniter]->SetMaxImpactXYPtDep(0.0182, 0.0350, -1.01);
	  //dtc1etaphitpc[aniter]->SetMaxSigmaToVertex(6.0);

	  cutPass1YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutPass1%stpcM%i", chrgs[ichg], imult), 0.13957);
	  cutFail1YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutFail1%stpcM%i", chrgs[ichg], imult), 0.13957);
	  dtc1etaphitpc[aniter]->AddCutMonitor(cutPass1YPtetaphitpc[aniter], cutFail1YPtetaphitpc[aniter]);
	  
	  cutPass1PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cutPass1%stpcM%i", chrgs[ichg], imult),0);
	  cutFail1PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cutFail1%stpcM%i", chrgs[ichg], imult),0);
	  dtc1etaphitpc[aniter]->AddCutMonitor(cutPass1PIDetaphitpc[aniter], cutFail1PIDetaphitpc[aniter]);
	  
	  //sqpcetaphitpc[aniter] = new AliFemtoPairCutAntiGamma();
	  sqpcetaphitpc[aniter] = new AliFemtoShareQualityTPCEntranceSepPairCut();
	  //sqpcetaphitpc[aniter] = new AliFemtoPairCutRadialDistance();
	  sqpcetaphitpc[aniter]->SetShareQualityMax(1.0);
	  sqpcetaphitpc[aniter]->SetShareFractionMax(0.05);
	  sqpcetaphitpc[aniter]->SetRemoveSameLabel(kFALSE);
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
	  anetaphitpc[aniter]->SetFirstParticleCut(dtc1etaphitpc[aniter]);
	  anetaphitpc[aniter]->SetSecondParticleCut(dtc1etaphitpc[aniter]);
	  anetaphitpc[aniter]->SetPairCut(sqpcetaphitpc[aniter]);
	  
	  //Correlation functions

	  //Spherical harmonics (without kT bins)
	  //cylmetaphitpc[aniter] = new AliFemtoCorrFctnDirectYlm(Form("cylm%stpcM%i", chrgs[ichg], imult),3,nbinssh,0.0,shqmax,runshlcms);
	  //anetaphitpc[aniter]->AddCorrFctn(cylmetaphitpc[aniter]);

	  //Qinv (without kT bins)
	  //cqinvkttpc[aniter] = new AliFemtoQinvCorrFctn(Form("cqinv%stpcM%i", chrgs[ichg], imult),nbinssh,0.0,shqmax);
	  //anetaphitpc[aniter]->AddCorrFctn(cqinvkttpc[aniter]);

	  //3D cartesian (without kT bins)
	  if(run3d){
	    cq3dlcmskttpc[aniter] = new AliFemtoCorrFctn3DLCMSSym(Form("cq3d%stpcM%i", chrgs[ichg], imult),100,0.6);
	    anetaphitpc[aniter]->AddCorrFctn(cq3dlcmskttpc[aniter]);
	  }
	  
	  // cqinvnclstpc[aniter] = new AliFemtoCorrFctnTPCNcls(Form("cqinvncls%stpcM%i", chrgs[ichg], imult),nbinssh,0.0,shqmax);
	  // anetaphitpc[aniter]->AddCorrFctn(cqinvnclstpc[aniter]);

	  // cqinvchi2tpc[aniter] = new AliFemtoChi2CorrFctn(Form("cqinvchi2%stpcM%i", chrgs[ichg], imult),nbinssh,0.0,shqmax);
	  // anetaphitpc[aniter]->AddCorrFctn(cqinvchi2tpc[aniter]);

	  if (runktdep) {
	    int ktm;
	    for (int ikt=0; ikt<7; ikt++) {
	      ktm = aniter*7 + ikt;
	      ktpcuts[ktm] = new AliFemtoKTPairCut(ktrng[ikt], ktrng[ikt+1]);
	      
  	      //cylmkttpc[ktm] = new AliFemtoCorrFctnDirectYlm(Form("cylm%stpcM%ikT%i", chrgs[ichg], imult, ikt),3,nbinssh, 0.0, shqmax, runshlcms);
  	      //cylmkttpc[ktm]->SetPairSelectionCut(ktpcuts[ktm]);
  	      //anetaphitpc[aniter]->AddCorrFctn(cylmkttpc[ktm]);
	      
	      //cqinvkttpc[ktm] = new AliFemtoQinvCorrFctn(Form("cqinv%stpcM%ikT%i", chrgs[ichg], imult, ikt),nbinssh,0.0, shqmax);
	      //cqinvkttpc[ktm]->SetPairSelectionCut(ktpcuts[ktm]);
	      //anetaphitpc[aniter]->AddCorrFctn(cqinvkttpc[ktm]);

	      //cqinvsqtpc[ktm] = new AliFemtoShareQualityCorrFctn(Form("cqinvsq%stpcM%ikT%i", chrgs[ichg], imult, ikt),nbinssh,0.0,shqmax);
	      //cqinvsqtpc[ktm]->SetPairSelectionCut(ktpcuts[ktm]);
	      //anetaphitpc[aniter]->AddCorrFctn(cqinvsqtpc[ktm]);

	      //cqinvinnertpc[ktm] = new AliFemtoTPCInnerCorrFctn(Form("cqinvinner%stpcM%ikT%i", chrgs[ichg], imult, ikt),nbinssh,0.0,shqmax);
	      //cqinvinnertpc[ktm]->SetPairSelectionCut(ktpcuts[ktm]);
	      //cqinvinnertpc[ktm]->SetRadius(1.2);
	      //anetaphitpc[aniter]->AddCorrFctn(cqinvinnertpc[ktm]);

	      if (run3d) {
		//		cq3dlcmskttpc[ktm] = new AliFemtoCorrFctn3DLCMSSym(Form("cq3d%stpcM%ikT%i", chrgs[ichg], imult, ikt),60,(imult>3)?((imult>6)?((imult>7)?0.6:0.4):0.25):0.15);
		cq3dlcmskttpc[ktm] = new AliFemtoCorrFctn3DLCMSSym(Form("cq3d%stpcM%ikT%i", chrgs[ichg], imult, ikt),100,0.6);
		cq3dlcmskttpc[ktm]->SetPairSelectionCut(ktpcuts[ktm]);
		anetaphitpc[aniter]->AddCorrFctn(cq3dlcmskttpc[ktm]);
	      }
	    }
	  }
	  
	  //cdedpetaphi[aniter] = new AliFemtoCorrFctnDEtaDPhi(Form("cdedp%stpcM%i", chrgs[ichg], imult),39, 39);
	  //anetaphitpc[aniter]->AddCorrFctn(cdedpetaphi[aniter]);
	  
	  Manager->AddAnalysis(anetaphitpc[aniter]);	
	}
      }
    }
  }
  // *** End pion-pion analysis

  return Manager;
}
