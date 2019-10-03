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

//_
AliFemtoManager* ConfigFemtoAnalysis() {


  double PionMass = 0.13956995;
  double KaonMass = 0.493677;
  double ProtonMass = 0.938272013;

  int runmults[10] = {1, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  int multbins[11] = {0.001, 50, 100, 200, 300, 400, 500, 600, 700, 800, 900};
  
  int runch[4] = {1, 1, 1, 1};
  const char *chrgs[4] = { "PPKP", "PMKM", "PPKM","PMKP"};

 int runktdep = 1;
  double ktrng[4] = {0.2, 0.3, 0.4, 0.5};
  //  double ktrng[8] = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0};

  int numOfMultBins = 10;  
  int numOfChTypes = 4;
  int numOfkTbins = 3;

  int runqinv = 1;
  int run3d = 1; // Do 3D cartesian analysis?
  int runshlcms = 1;// 0:PRF(PAP), 1:LCMS(PP,APAP)

  int runtype = 2; // Types 0 - global, 1 - ITS only, 2 - TPC Inner
  int isrealdata = 1;

  int gammacut = 1;
  double shqmax;
  if (runshlcms) shqmax = 2.0;
  else shqmax = 0.9;

  int nbinssh = 200;

  AliFemtoEventReaderAODChain *Reader = new AliFemtoEventReaderAODChain();
  Reader->SetFilterBit(7);
  Reader->SetCentralityPreSelection(0.001, 950);

  //AliFemtoEventReaderESDChainKine* Reader=new AliFemtoEventReaderESDChainKine();
  //Reader->SetConstrained(true);
  // Reader->SetUseTPCOnly(false);

  //AliFemtoEventReaderESDChain *Reader = new AliFemtoEventReaderESDChain();
  //Reader->SetUseMultiplicity(AliFemtoEventReaderESDChain::kV0Centrality);
  //Reader->SetReadTrackType(AliFemtoEventReaderESDChain::kGlobal);


  AliFemtoManager* Manager=new AliFemtoManager();
  Manager->SetEventReader(Reader);

  AliFemtoVertexMultAnalysis    *anetaphitpc[40];
  AliFemtoBasicEventCut         *mecetaphitpc[40];
  AliFemtoCutMonitorEventMult   *cutPassEvMetaphitpc[40];
  AliFemtoCutMonitorEventMult   *cutFailEvMetaphitpc[40];
  AliFemtoCutMonitorEventVertex *cutPassEvVetaphitpc[40];
  AliFemtoCutMonitorEventVertex *cutFailEvVetaphitpc[40];
  AliFemtoCutMonitorCollections   *cutPassColletaphitpc[40];
  AliFemtoCutMonitorCollections   *cutFailColletaphitpc[40];
  AliFemtoESDTrackCut           *dtc1etaphitpc[40];
  AliFemtoESDTrackCut           *dtc2etaphitpc[40];
  AliFemtoCutMonitorParticleYPt *cutPass1YPtetaphitpc[40];
  AliFemtoCutMonitorParticleYPt *cutFail1YPtetaphitpc[40];
  AliFemtoCutMonitorParticlePID *cutPass1PIDetaphitpc[40];
  AliFemtoCutMonitorParticlePID *cutFail1PIDetaphitpc[40];
  AliFemtoCutMonitorParticleYPt *cutPass2YPtetaphitpc[40];
  AliFemtoCutMonitorParticleYPt *cutFail2YPtetaphitpc[40];
  AliFemtoCutMonitorParticlePID *cutPass2PIDetaphitpc[40];
  AliFemtoCutMonitorParticlePID *cutFail2PIDetaphitpc[40];
  //  AliFemtoPairCutAntiGamma      *sqpcetaphitpc[20];
  //AliFemtoShareQualityTPCEntranceSepPairCut      *sqpcetaphitpc[20];
  AliFemtoPairCutRadialDistance      *sqpcetaphitpc[40];
  AliFemtoPairCutRadialDistanceLM      *sqpcetaphitpcRD[40];
  AliFemtoCorrFctnDirectYlm     *cylmetaphitpc[40];
  AliFemtoCorrFctnDEtaDPhi      *cdedpetaphi[40];
  AliFemtoChi2CorrFctn          *cchiqinvetaphitpc[40];
  AliFemtoKTPairCut             *ktpcuts[40*7];
  AliFemtoCorrFctnDirectYlm     *cylmkttpc[40*7];
  AliFemtoQinvCorrFctn          *cqinvkttpc[40*7];
  AliFemtoCorrFctn3DLCMSSym     *cq3dlcmskttpc[40*7];
  AliFemtoCorrFctn3DSpherical   *cq3dspherical[40*7];
  AliFemtoCorrFctnTPCNcls       *cqinvnclstpc[40];
  AliFemtoShareQualityCorrFctn  *cqinvsqtpc[40*10];
  AliFemtoChi2CorrFctn          *cqinvchi2tpc[40];
  AliFemtoTPCInnerCorrFctn      *cqinvinnertpc[40*10];

  // *** Begin pion-pion (positive) analysis ***
  int aniter = 0;  

  for (int imult = 0; imult < numOfMultBins; imult++) {
    if (runmults[imult]) {

      for (int ichg = 0; ichg < numOfChTypes; ichg++) {
	if (runch[ichg]) {

	  aniter = ichg * numOfMultBins + imult;
	  anetaphitpc[aniter] = new AliFemtoVertexMultAnalysis(10, -10.0, 10.0, 4, multbins[imult], multbins[imult+1]);
	  anetaphitpc[aniter]->SetNumEventsToMix(5);
	  anetaphitpc[aniter]->SetMinSizePartCollection(1);

	  mecetaphitpc[aniter] = new AliFemtoBasicEventCut();
	  mecetaphitpc[aniter]->SetEventMult(0.001,100000);
	  mecetaphitpc[aniter]->SetVertZPos(-10,10);

	  //if (isrealdata)
	  //  mecetaphitpc[aniter]->SetAcceptOnlyPhysics(kTRUE);
	  	  
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
	  dtc2etaphitpc[aniter] = new AliFemtoESDTrackCut();

	  if (ichg == 0)
	    {
	      dtc1etaphitpc[aniter]->SetCharge(1.0);
	      dtc2etaphitpc[aniter]->SetCharge(1.0);
	    }
	  else if (ichg == 1)
	    {
	      dtc1etaphitpc[aniter]->SetCharge(-1.0);
	      dtc2etaphitpc[aniter]->SetCharge(-1.0);
	    }
	  else if (ichg == 2) 
	    {
	      dtc1etaphitpc[aniter]->SetCharge(1.0);
	      dtc2etaphitpc[aniter]->SetCharge(-1.0);
	    }
	  else if (ichg == 3) 
	    {
	      dtc1etaphitpc[aniter]->SetCharge(-1.0);
	      dtc2etaphitpc[aniter]->SetCharge(1.0);
	    }

	  dtc1etaphitpc[aniter]->SetPt(0.14,1.5);
	  dtc1etaphitpc[aniter]->SetEta(-0.8,0.8);
	  dtc1etaphitpc[aniter]->SetMass(PionMass);	  
	  dtc1etaphitpc[aniter]->SetMostProbablePion();	 
 
	  dtc2etaphitpc[aniter]->SetPt(0.14,1.5);
          dtc2etaphitpc[aniter]->SetEta(-0.8,0.8);
	  dtc2etaphitpc[aniter]->SetMass(KaonMass);	  
	  dtc2etaphitpc[aniter]->SetMostProbableKaon();
	  
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
            
            dtc2etaphitpc[aniter]->SetStatus(AliESDtrack::kTPCin);
	    dtc2etaphitpc[aniter]->SetminTPCncls(80);
	    dtc2etaphitpc[aniter]->SetRemoveKinks(kTRUE);
	    dtc2etaphitpc[aniter]->SetLabel(kFALSE);
	    dtc2etaphitpc[aniter]->SetMaxTPCChiNdof(4.0);
	    dtc2etaphitpc[aniter]->SetMaxImpactXY(2.4);
	    dtc2etaphitpc[aniter]->SetMaxImpactZ(3.0);
             }

	  cutPass1YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutPass1%stpcM%i", chrgs[ichg], imult),PionMass);
	  cutFail1YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutFail1%stpcM%i", chrgs[ichg], imult),PionMass);
	  dtc1etaphitpc[aniter]->AddCutMonitor(cutPass1YPtetaphitpc[aniter], cutFail1YPtetaphitpc[aniter]);
	  
	  cutPass1PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cutPass1%stpcM%i", chrgs[ichg], imult),0);//0-pion,1-kaon,2-proton
	  cutFail1PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cutFail1%stpcM%i", chrgs[ichg], imult),0);
	  dtc1etaphitpc[aniter]->AddCutMonitor(cutPass1PIDetaphitpc[aniter], cutFail1PIDetaphitpc[aniter]);

//================================================================
	      cutPass2YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutPass2%stpcM%i", chrgs[ichg], imult),KaonMass);
	      cutFail2YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutFail2%stpcM%i", chrgs[ichg], imult),KaonMass);
	      dtc2etaphitpc[aniter]->AddCutMonitor(cutPass2YPtetaphitpc[aniter], cutFail2YPtetaphitpc[aniter]);

	      cutPass2PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cutPass2%stpcM%i", chrgs[ichg], imult),1);//0-pion,1-kaon,2-proton
	      cutFail2PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cutFail2%stpcM%i", chrgs[ichg], imult),1);
	      dtc2etaphitpc[aniter]->AddCutMonitor(cutPass2PIDetaphitpc[aniter], cutFail2PIDetaphitpc[aniter]);


	 sqpcetaphitpc[aniter] = new AliFemtoPairCutRadialDistance();
          sqpcetaphitpcRD[aniter] = new AliFemtoPairCutRadialDistanceLM();

	/*  if( runtype==0){
	   sqpcetaphitpc[aniter]->SetMaxEEMinv(0.0);
	   sqpcetaphitpc[aniter]->SetMaxThetaDiff(0.0);
	   sqpcetaphitpc[aniter]->SetTPCEntranceSepMinimum(1.5);
	  // sqpcetaphitpc[aniter]->SetPhiStarDistanceMinimum(0.03);
	  // sqpcetaphitpc[aniter]->SetRadialDistanceMinimum(0.12, 0.03);
	  //sqpcetaphitpc[aniter]->SetEtaDifferenceMinimum(0.02);
        }

	  else if (runtype==1){
	  	    sqpcetaphitpc[aniter]->SetTPCEntranceSepMinimum(5.0);
	  	    sqpcetaphitpc[aniter]->SetPhiStarDistanceMinimum(0.03);
	  //sqpcetaphitpc[aniter]->SetRadialDistanceMinimum(1.2, 0.03);
	  //sqpcetaphitpc[aniter]->SetEtaDifferenceMinimum(0.02);
}
	else if ( runtype==2 ){
	  sqpcetaphitpc[aniter]->SetShareQualityMax(1.0);
	  sqpcetaphitpc[aniter]->SetShareFractionMax(0.05);
	  sqpcetaphitpc[aniter]->SetRemoveSameLabel(kFALSE);
          // sqpcetaphitpc[aniter]->SetMinimumRadius(1.6);
	  //sqpcetaphitpc[aniter]->SetEtaDifferenceMinimum(0.02);
	  //sqpcetaphitpc[aniter]->SetPhiStarDifferenceMinimum(0.045);  
          //sqpcetaphitpcRD[aniter]->SetMinimumRadius(1.6);
	  //sqpcetaphitpcRD[aniter]->SetEtaDifferenceMinimum(0.02);
	  //sqpcetaphitpcRD[aniter]->SetPhiStarDifferenceMinimum(0.045);
	  }
*/
	  anetaphitpc[aniter]->SetEventCut(mecetaphitpc[aniter]);
	  anetaphitpc[aniter]->SetFirstParticleCut(dtc1etaphitpc[aniter]);
	  anetaphitpc[aniter]->SetSecondParticleCut(dtc2etaphitpc[aniter]);
	 anetaphitpc[aniter]->SetPairCut(sqpcetaphitpc[aniter]);
	  
	  

	  if (runktdep) {
	    int ktm;
	    for (int ikt=0; ikt<numOfkTbins; ikt++) {
	      ktm = aniter*numOfkTbins + ikt;
	      ktpcuts[ktm] = new AliFemtoKTPairCut(ktrng[ikt], ktrng[ikt+1]);
	      
              cqinvkttpc[ktm] = new AliFemtoQinvCorrFctn(Form("cqinv%stpcM%ikT%i", chrgs[ichg], imult, ikt),nbinssh,0.0, shqmax);
	      cqinvkttpc[ktm]->SetPairSelectionCut(ktpcuts[ktm]);
	      anetaphitpc[aniter]->AddCorrFctn(cqinvkttpc[ktm]);

	      if (run3d) {

		cq3dlcmskttpc[ktm] = new AliFemtoCorrFctn3DLCMSSym(Form("cq3d%stpcM%ikT%i", chrgs[ichg], imult, ikt),100,0.5);
		cq3dlcmskttpc[ktm]->SetPairSelectionCut(ktpcuts[ktm]);
		anetaphitpc[aniter]->AddCorrFctn(cq3dlcmskttpc[ktm]);

              /* cq3dspherical[ktm] = new AliFemtoCorrFctn3DSpherical(Form("cq3d%stpcM%ikT%i",chrgs[ichg], imult, ikt), nbinssh,0.0,shqmax,nbinssh,nbinssh);
               cq3dspherical[ktm]->SetPairSelectionCut(ktpcuts[ktm]);
               anetaphitpc[aniter]->AddCorrFctn(cq3dspherical[ktm]);*/

	      }
	    }
	  }
	  
	
	  
	  Manager->AddAnalysis(anetaphitpc[aniter]);	
	}
      }
    }
  }
  // *** End pion-pion analysis

  return Manager;
}




















