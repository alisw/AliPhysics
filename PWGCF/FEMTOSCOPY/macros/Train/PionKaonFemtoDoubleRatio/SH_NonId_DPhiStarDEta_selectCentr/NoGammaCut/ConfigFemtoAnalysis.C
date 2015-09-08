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
#include "AliFemtoCorrFctnNonIdDR.h"
#include "AliFemtoCorrFctnDPhiStarDEta.h"
#endif

//_
AliFemtoManager* ConfigFemtoAnalysis(int runcentrality0, int runcentrality1, int runcentrality2, int runcentrality3, int runcentrality4,int runcentrality5, int runcentrality6, int runcentrality7, int runcentrality8, int runcentrality9, int runSHCorrFctn, int runNonIdCorrFctn) {


  double PionMass = 0.13956995;
  double KaonMass = 0.493677;
  double ProtonMass = 0.938272013;

  const int numOfMultBins = 10;  
  const int numOfChTypes = 4;
  const int numOfkTbins = 7;

  int runmults[numOfMultBins] = {runcentrality0, runcentrality1, runcentrality2, runcentrality3, runcentrality4, runcentrality5, runcentrality6, runcentrality7, runcentrality8, runcentrality9};
  int multbins[numOfMultBins + 1] = {0.001, 50, 100, 200, 300, 400, 500, 600, 700, 800, 900};
  
  int runch[numOfChTypes] = {1, 1, 1, 1};
  const char *chrgs[numOfChTypes] = { "PIpKp", "PImKm", "PIpKm","PImKp"};

  int runktdep = 0;
  double ktrng[numOfkTbins + 1] = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0};



  int runqinv = 1;
  int run3d = 0; // Do 3D cartesian analysis?
  int runshlcms = 0;// 0:PRF(PAP), 1:LCMS(PP,APAP)

  int runtype = 2; // Types 0 - global, 1 - ITS only, 2 - TPC Inner
  int isrealdata = 1;

  int gammacut = 1;
  double shqmax = 0.5;
  //if (runshlcms) shqmax = 2.0;
  //else shqmax = 0.9;

  int nbinssh = 100;

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

  const int size = numOfMultBins * numOfChTypes;
  AliFemtoVertexMultAnalysis    *anetaphitpc[size];
  AliFemtoBasicEventCut         *mecetaphitpc[size];
  AliFemtoCutMonitorEventMult   *cutPassEvMetaphitpc[size];
  AliFemtoCutMonitorEventMult   *cutFailEvMetaphitpc[size];
  AliFemtoCutMonitorEventVertex *cutPassEvVetaphitpc[size];
  AliFemtoCutMonitorEventVertex *cutFailEvVetaphitpc[size];
  AliFemtoCutMonitorCollections   *cutPassColletaphitpc[size];
  AliFemtoCutMonitorCollections   *cutFailColletaphitpc[size];
  AliFemtoESDTrackCut           *dtc1etaphitpc[size];
  AliFemtoESDTrackCut           *dtc2etaphitpc[size];
  AliFemtoCutMonitorParticleYPt *cutPass1YPtetaphitpc[size];
  AliFemtoCutMonitorParticleYPt *cutFail1YPtetaphitpc[size];
  AliFemtoCutMonitorParticlePID *cutPass1PIDetaphitpc[size];
  AliFemtoCutMonitorParticlePID *cutFail1PIDetaphitpc[size];
  AliFemtoCutMonitorParticleYPt *cutPass2YPtetaphitpc[size];
  AliFemtoCutMonitorParticleYPt *cutFail2YPtetaphitpc[size];
  AliFemtoCutMonitorParticlePID *cutPass2PIDetaphitpc[size];
  AliFemtoCutMonitorParticlePID *cutFail2PIDetaphitpc[size];
  AliFemtoPairCutAntiGamma      *sqpcetaphitpc[size];
  //AliFemtoShareQualityTPCEntranceSepPairCut      *sqpcetaphitpc[20];
  //AliFemtoPairCutRadialDistance      *sqpcetaphitpc[size];
  AliFemtoPairCutRadialDistanceLM      *sqpcetaphitpcRD[size];
  AliFemtoCorrFctnDirectYlm     *cylmetaphitpc[size];
  AliFemtoCorrFctnDEtaDPhi      *cdedpetaphi[size];
  AliFemtoChi2CorrFctn          *cchiqinvetaphitpc[size];
  AliFemtoKTPairCut             *ktpcuts[size*numOfkTbins];
  AliFemtoCorrFctnDirectYlm     *cylmkttpc[size*numOfkTbins];
  AliFemtoQinvCorrFctn          *cqinvkttpc[size*numOfkTbins];
  AliFemtoCorrFctn3DLCMSSym     *cq3dlcmskttpc[size*numOfkTbins];
  AliFemtoCorrFctn3DSpherical   *cq3dspherical[size*numOfkTbins];
  AliFemtoCorrFctnTPCNcls       *cqinvnclstpc[size];
  AliFemtoChi2CorrFctn          *cqinvchi2tpc[size];

  AliFemtoModelGausLCMSFreezeOutGenerator *gausLCMSFreezeOutGenerator[size];
  AliFemtoModelWeightGeneratorBasic      *weightGeneratorLednicky[size];
  AliFemtoModelManager          *tModelManager[size];
  AliFemtoCorrFctnNonIdDR       *cnonidtpc[size*numOfkTbins];
  AliFemtoCorrFctnDPhiStarDEta  *cdphistardeta08[size*numOfkTbins];
  AliFemtoCorrFctnDPhiStarDEta  *cdphistardeta12[size*numOfkTbins];
  AliFemtoCorrFctnDPhiStarDEta  *cdphistardeta16[size*numOfkTbins];
  AliFemtoCorrFctnDPhiStarDEta  *cdphistardeta20[size*numOfkTbins];


  // *** Begin pion-kaon analysis ***
  int aniter = 0;  

  for (int imult = 0; imult < numOfMultBins; imult++) {
    if (runmults[imult]) {

      for (int ichg = 0; ichg < numOfChTypes; ichg++) {
	if (runch[ichg]) {
	  
	  //Iterator:
	  aniter = ichg * numOfMultBins + imult;


	  //Mix events with respect to the z position of the primary vertex and event total multipliticy:
	  anetaphitpc[aniter] = new AliFemtoVertexMultAnalysis(10, -10.0, 10.0, 4, multbins[imult], multbins[imult+1]);
	  anetaphitpc[aniter]->SetNumEventsToMix(5);
	  anetaphitpc[aniter]->SetMinSizePartCollection(1);

	  //Select basic cuts:
	  mecetaphitpc[aniter] = new AliFemtoBasicEventCut();
	  mecetaphitpc[aniter]->SetEventMult(0.001,100000);
	  mecetaphitpc[aniter]->SetVertZPos(-10,10);

	  //Study the multiplicity distribution:
	  cutPassEvMetaphitpc[aniter] = new AliFemtoCutMonitorEventMult(Form("cutPass%stpcM%i", chrgs[ichg], imult),500);
	  cutFailEvMetaphitpc[aniter] = new AliFemtoCutMonitorEventMult(Form("cutFail%stpcM%i", chrgs[ichg], imult),500);
	  mecetaphitpc[aniter]->AddCutMonitor(cutPassEvMetaphitpc[aniter], cutFailEvMetaphitpc[aniter]);
	  //Study the distribution and error of the primary vertex:
	  cutPassEvVetaphitpc[aniter] = new AliFemtoCutMonitorEventVertex(Form("cutPass%stpcM%i", chrgs[ichg], imult));
	  cutFailEvVetaphitpc[aniter] = new AliFemtoCutMonitorEventVertex(Form("cutFail%stpcM%i", chrgs[ichg], imult));
	  mecetaphitpc[aniter]->AddCutMonitor(cutPassEvVetaphitpc[aniter], cutFailEvVetaphitpc[aniter]);
	  //Study the multiplicity distribution:
          cutPassColletaphitpc[aniter] = new AliFemtoCutMonitorCollections(Form("cutPass%stpcM%i", chrgs[ichg], imult));
          cutFailColletaphitpc[aniter] = new AliFemtoCutMonitorCollections(Form("cutFail%stpcM%i", chrgs[ichg], imult));
          mecetaphitpc[aniter]->AddCutMonitor(cutPassColletaphitpc[aniter], cutFailColletaphitpc[aniter]);


	  //Basic track cut:
	  dtc1etaphitpc[aniter] = new AliFemtoESDTrackCut();
	  dtc2etaphitpc[aniter] = new AliFemtoESDTrackCut();

	  //Set charge of particles:
	  if (ichg == 0) {
	    dtc1etaphitpc[aniter]->SetCharge(1.0);
	    dtc2etaphitpc[aniter]->SetCharge(1.0);
	  }
	  else if (ichg == 1) {
	    dtc1etaphitpc[aniter]->SetCharge(-1.0);
	    dtc2etaphitpc[aniter]->SetCharge(-1.0);
	  }
	  else if (ichg == 2) {
	    dtc1etaphitpc[aniter]->SetCharge(1.0);
	    dtc2etaphitpc[aniter]->SetCharge(-1.0);
	  }
	  else if (ichg == 3) {
	    dtc1etaphitpc[aniter]->SetCharge(-1.0);
	    dtc2etaphitpc[aniter]->SetCharge(1.0);
	  }

	  //Set particle 1:
	  dtc1etaphitpc[aniter]->SetPt(0.14,1.5);
	  dtc1etaphitpc[aniter]->SetEta(-0.8,0.8);
	  dtc1etaphitpc[aniter]->SetMass(PionMass);	  
	  dtc1etaphitpc[aniter]->SetMostProbablePion();	 
 
	  //Set particle 2:
	  dtc2etaphitpc[aniter]->SetPt(0.14,1.5);
          dtc2etaphitpc[aniter]->SetEta(-0.8,0.8);
	  dtc2etaphitpc[aniter]->SetMass(KaonMass);	  
	  dtc2etaphitpc[aniter]->SetMostProbableKaon();
	  

	  // Track quality cuts

	  /*if (runtype == 0) {
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
	    }*/

	  //============PION============

	  //The cut monitor for particles to study the difference between reconstructed and true momentum: 
	  cutPass1YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutPass1%stpcM%i", chrgs[ichg], imult),PionMass);
	  cutFail1YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutFail1%stpcM%i", chrgs[ichg], imult),PionMass);
	  dtc1etaphitpc[aniter]->AddCutMonitor(cutPass1YPtetaphitpc[aniter], cutFail1YPtetaphitpc[aniter]);
	  //The cut monitor for particles to study various aspects of the PID determination:
	  cutPass1PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cutPass1%stpcM%i", chrgs[ichg], imult),0);//0-pion,1-kaon,2-proton
	  cutFail1PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cutFail1%stpcM%i", chrgs[ichg], imult),0);
	  dtc1etaphitpc[aniter]->AddCutMonitor(cutPass1PIDetaphitpc[aniter], cutFail1PIDetaphitpc[aniter]);

	  //============KAON============

	  cutPass2YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutPass2%stpcM%i", chrgs[ichg], imult),KaonMass);
	  cutFail2YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutFail2%stpcM%i", chrgs[ichg], imult),KaonMass);
	  dtc2etaphitpc[aniter]->AddCutMonitor(cutPass2YPtetaphitpc[aniter], cutFail2YPtetaphitpc[aniter]);

	  cutPass2PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cutPass2%stpcM%i", chrgs[ichg], imult),1);//0-pion,1-kaon,2-proton
	  cutFail2PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cutFail2%stpcM%i", chrgs[ichg], imult),1);
	  dtc2etaphitpc[aniter]->AddCutMonitor(cutPass2PIDetaphitpc[aniter], cutFail2PIDetaphitpc[aniter]);

	  //A pair cut which checks for some pair qualities that attempt to identify slit/doubly reconstructed tracks and also selects pairs based on their separation at the entrance to the TPC
	  //sqpcetaphitpc[aniter] = new AliFemtoPairCutRadialDistance();
          //sqpcetaphitpcRD[aniter] = new AliFemtoPairCutRadialDistanceLM();

	  sqpcetaphitpc[aniter] = new AliFemtoPairCutAntiGamma();
	  sqpcetaphitpc[aniter]->SetDataType(AliFemtoPairCut::kKine);

	  if(gammacut == 1) {
	    sqpcetaphitpc[aniter]->SetMaxEEMinv(0.002);
	    sqpcetaphitpc[aniter]->SetMaxThetaDiff(0.008);
	  }

	  if( runtype==0){
	    sqpcetaphitpc[aniter]->SetMaxEEMinv(0.0);
	    sqpcetaphitpc[aniter]->SetMaxThetaDiff(0.0);
	    sqpcetaphitpc[aniter]->SetTPCEntranceSepMinimum(1.5);
	    // sqpcetaphitpc[aniter]->SetPhiStarDistanceMinimum(0.03);
	    // sqpcetaphitpc[aniter]->SetRadialDistanceMinimum(0.12, 0.03);
	    //sqpcetaphitpc[aniter]->SetEtaDifferenceMinimum(0.02);
	  }

	  else if (runtype==1){
	    //sqpcetaphitpc[aniter]->SetTPCEntranceSepMinimum(5.0);
	    //sqpcetaphitpc[aniter]->SetPhiStarDistanceMinimum(0.03);

	    //sqpcetaphitpc[aniter]->SetRadialDistanceMinimum(1.2, 0.03);
	    //sqpcetaphitpc[aniter]->SetEtaDifferenceMinimum(0.02);
	  }
	  else if ( runtype==2 ){
	    sqpcetaphitpc[aniter]->SetShareQualityMax(1.0);
	    sqpcetaphitpc[aniter]->SetShareFractionMax(0.05);
	    sqpcetaphitpc[aniter]->SetRemoveSameLabel(kFALSE);
	    sqpcetaphitpc[aniter]->SetTPCEntranceSepMinimum(0.0);
	    // sqpcetaphitpc[aniter]->SetMinimumRadius(1.6);
	    //sqpcetaphitpc[aniter]->SetEtaDifferenceMinimum(0.02);
	    //sqpcetaphitpc[aniter]->SetPhiStarDifferenceMinimum(0.045);  
	    //sqpcetaphitpcRD[aniter]->SetMinimumRadius(1.6);
	    //sqpcetaphitpcRD[aniter]->SetEtaDifferenceMinimum(0.02);
	    //sqpcetaphitpcRD[aniter]->SetPhiStarDifferenceMinimum(0.045);
	  }

	  anetaphitpc[aniter]->SetEventCut(mecetaphitpc[aniter]);
	  anetaphitpc[aniter]->SetFirstParticleCut(dtc1etaphitpc[aniter]);
	  anetaphitpc[aniter]->SetSecondParticleCut(dtc2etaphitpc[aniter]);
	  anetaphitpc[aniter]->SetPairCut(sqpcetaphitpc[aniter]);
	  
	  
	  //Correlation functions

	  //Spherical harmonics (without kT bins)
	  if(runSHCorrFctn == 1) {
	    cylmetaphitpc[aniter] = new AliFemtoCorrFctnDirectYlm(Form("cylm%stpcM%i", chrgs[ichg], imult),3,nbinssh,0.0,shqmax,runshlcms);
	    anetaphitpc[aniter]->AddCorrFctn(cylmetaphitpc[aniter]);
	  }
	  
	  //Qinv (without kT bins)
	  //cqinvkttpc[aniter] = new AliFemtoQinvCorrFctn(Form("cqinv%stpcM%i", chrgs[ichg], imult),nbinssh,0.0, shqmax);
	  //anetaphitpc[aniter]->AddCorrFctn(cqinvkttpc[aniter]);

	  // NonId (without kT bins)
	  // Correlation function for non-identical particles uses k* as a function variable. Stores the correlation 
	  // function separately for positive and negative signs of k* projections into out, side and long directions, 
	  // enabling the calculations of double ratios
	  if(runNonIdCorrFctn == 1) {
	    cnonidtpc[aniter] = new AliFemtoCorrFctnNonIdDR(Form("cnonid%stpcM%i", chrgs[ichg], imult), nbinssh, 0.0, shqmax);
	    anetaphitpc[aniter]->AddCorrFctn(cnonidtpc[aniter]);
	  }
	  
	  // DPhiStarDEta (without bins)
	  // Correlation function for two particle correlations which uses dPhi* and dEta as a function variables.
	  cdphistardeta08[aniter] = new AliFemtoCorrFctnDPhiStarDEta(Form("cdphistardeta08%stpcM%i", chrgs[ichg], imult), 0.8, 51, -0.05, 0.05, 51, -0.4, 0.4);
	  anetaphitpc[aniter]->AddCorrFctn(cdphistardeta08[aniter]);
	  cdphistardeta12[aniter] = new AliFemtoCorrFctnDPhiStarDEta(Form("cdphistardeta12%stpcM%i", chrgs[ichg], imult), 1.2, 51, -0.05, 0.05, 51, -0.4, 0.4);
	  anetaphitpc[aniter]->AddCorrFctn(cdphistardeta12[aniter]);
	  cdphistardeta16[aniter] = new AliFemtoCorrFctnDPhiStarDEta(Form("cdphistardeta16%stpcM%i", chrgs[ichg], imult), 1.6, 51, -0.05, 0.05, 51, -0.4, 0.4);
	  anetaphitpc[aniter]->AddCorrFctn(cdphistardeta16[aniter]);
	  cdphistardeta20[aniter] = new AliFemtoCorrFctnDPhiStarDEta(Form("cdphistardeta20%stpcM%i", chrgs[ichg], imult), 2.0, 51, -0.05, 0.05, 51, -0.4, 0.4);
	  anetaphitpc[aniter]->AddCorrFctn(cdphistardeta20[aniter]);

	  if (runktdep) {
	    int ktm;
	    for (int ikt=0; ikt<numOfkTbins; ikt++) {
	      ktm = aniter*numOfkTbins + ikt;
	      ktpcuts[ktm] = new AliFemtoKTPairCut(ktrng[ikt], ktrng[ikt+1]);
	      

  	      cylmkttpc[ktm] = new AliFemtoCorrFctnDirectYlm(Form("cylm%stpcM%ikT%i", chrgs[ichg], imult, ikt),3,nbinssh, 0.0,shqmax, runshlcms);
  	      cylmkttpc[ktm]->SetPairSelectionCut(ktpcuts[ktm]);
  	      anetaphitpc[aniter]->AddCorrFctn(cylmkttpc[ktm]);

              //cqinvkttpc[ktm] = new AliFemtoQinvCorrFctn(Form("cqinv%stpcM%ikT%i", chrgs[ichg], imult, ikt),nbinssh,0.0, shqmax);
	      //cqinvkttpc[ktm]->SetPairSelectionCut(ktpcuts[ktm]);
	      //anetaphitpc[aniter]->AddCorrFctn(cqinvkttpc[ktm]);

	      //cnonidtpc[ktm] = new AliFemtoCorrFctnNonIdDR(Form("cnonid%stpcM%ikT%i", chrgs[ichg], imult, ikt), nbinssh, 0.0, shqmax);
	      //cnonidtpc[ktm]->SetPairSelectionCut(ktpcuts[ktm]);
	      //anetaphitpc[aniter]->AddCorrFctn(cnonidtpc[ktm]);

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
  // *** End pion-kaon analysis

  return Manager;
}




















