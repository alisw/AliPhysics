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
#include "AliFemtoCutMonitorParticlePIDBeta.h"
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
#include "AliFemtoPairCutRadialDistanceAsymmetric.h"
#include "AliFemtoKKTrackCut.h"
#include "AliFemtoKKTrackCutFull.h"
#include "AliFemtoPairCutMergedFraction.h"
#include "AliFemtoCorrFctnDPhiStarKStarMergedFraction.h"
#include "AliFemtoCorrFctnDPhiStarKStarAverageMergedPointsFraction.h"
#include "AliFemtoBetaTPairCut.h"
#include "AliFemtoCutMonitorPairBetaT.h"
#endif
AliFemtoEventReaderAODMultSelection* GetReader2015(bool mcAnalysis);
AliFemtoEventReaderAODChain* GetReader2011(bool mcAnalysis);

AliFemtoEventReaderAODMultSelection* GetReader2015(bool mcAnalysis)
{
  AliFemtoEventReaderAODMultSelection* Reader = new AliFemtoEventReaderAODMultSelection();
  Reader->SetFilterMask(96);
  //Reader->SetReadV0(1);
  Reader->SetUseMultiplicity(AliFemtoEventReaderAOD::kCentrality);
  Reader->SetEPVZERO(kTRUE);
  Reader->SetCentralityFlattening(kTRUE);
  //Reader->SetReadCascade(kTRUE);
  Reader->SetPrimaryVertexCorrectionTPCPoints(kTRUE);

  Reader->SetUseAliEventCuts(kTRUE);
  Reader->SetTrackPileUpRemoval(kTRUE);
  
  if(mcAnalysis) Reader->SetReadMC(kTRUE);
  
  return Reader;
}

AliFemtoEventReaderAODChain* GetReader2011(bool mcAnalysis)
{
  AliFemtoEventReaderAODChain* Reader = new AliFemtoEventReaderAODChain();
  Reader->SetFilterBit(7);
  //Reader->SetReadV0(1);
  Reader->SetUseMultiplicity(AliFemtoEventReaderAOD::kCentrality);
  Reader->SetEPVZERO(kTRUE);
  Reader->SetCentralityFlattening(kTRUE);
  //Reader->SetReadCascade(kTRUE);
  Reader->SetPrimaryVertexCorrectionTPCPoints(kTRUE);
  if(mcAnalysis) Reader->SetReadMC(kTRUE);
  
  return Reader;
}

//_
AliFemtoManager* ConfigFemtoAnalysis(int runcentrality0, int runcentrality1, int runcentrality2, int runcentrality3, int runcentrality4,int runcentrality5, int runcentrality6, int runSHCorrFctn, int runNonIdCorrFctn, int paircutantigammaon, int paircutmergedfractionon, double distance, double fraction1, int runDPhiStarKStarMergedFraction, int runDPhiStarKStarAverageMergedPointsFraction, int runDPhiStarDEta, int turnOnMonitors, int turnOnBetaTMonitor, int runbetatdep, int runbetatylm, int runbetatnonid,int year) {


  double PionMass = 0.13957018;//0.13956995;
  double KaonMass = 0.493677;
  double ProtonMass = 0.938272013;
  double DeuteronMass =1.8756;

  const int numOfMultBins = 7;  
  const int numOfChTypes = 4;
  const int numOfkTbins = 2;

  int runmults[numOfMultBins] = {runcentrality0, runcentrality1, runcentrality2, runcentrality3, runcentrality4, runcentrality5, runcentrality6};
  int multbins[numOfMultBins + 1] = {0, 50, 100, 200, 300, 400, 500, 900};
  
  int runch[numOfChTypes] = {1, 1, 1, 1};
  const char *chrgs[numOfChTypes] = { "KpDp", "KmDm", "KpDm","KmDp"};

  //int runktdep = 1;
  double ktrng[numOfkTbins + 1] = {0.5, 0.85, 1.0};
  
  int runshlcms = 0;// 0:PRF(PAP), 1:LCMS(PP,APAP)

  int runtype = 2; // Types 0 - global, 1 - ITS only, 2 - TPC Inner

  int gammacut = 1;
  double shqmax = 1.0;
  //if (runshlcms) shqmax = 2.0;
  //else shqmax = 0.9;

  int nbinssh = 200;

  // create analysis managers
  AliFemtoManager* Manager=new AliFemtoManager();
  AliFemtoModelManager *modelMgr = new AliFemtoModelManager();
 // int year=2015;
  // add event reader
  if(year==2015)
  {
    AliFemtoEventReaderAODMultSelection* Reader2015 = GetReader2015(kFALSE);
    Manager->SetEventReader(Reader2015);
  }
  else if(year==2011)
  {
    AliFemtoEventReaderAODChain* Reader2011 = GetReader2011(kFALSE);
    Manager->SetEventReader(Reader2011);
  }
  

  /*
  AliFemtoEventReaderAODChain *Reader = new AliFemtoEventReaderAODChain();
  if(mask==96){
     Reader->SetFilterMask(96);}//globaltrack
  else if(mask==7){
     Reader->SetFilterBit(7);}//tpconlytrack
  else if(mask==768){
     Reader->SetFilterMask(768);}//hybridtrack

  Reader->SetCentralityPreSelection(0, 950);

  AliFemtoManager* Manager=new AliFemtoManager();
  Manager->SetEventReader(Reader);
  */




  const int size = numOfMultBins * numOfChTypes;
  AliFemtoVertexMultAnalysis    *anetaphitpc[size];
  AliFemtoBasicEventCut         *mecetaphitpc[size];
  AliFemtoCutMonitorEventMult   *cutPassEvMetaphitpc[size];
  AliFemtoCutMonitorEventMult   *cutFailEvMetaphitpc[size];
  AliFemtoCutMonitorEventVertex *cutPassEvVetaphitpc[size];
  AliFemtoCutMonitorEventVertex *cutFailEvVetaphitpc[size];
  AliFemtoCutMonitorCollections   *cutPassColletaphitpc[size];
  AliFemtoCutMonitorCollections   *cutFailColletaphitpc[size];
  //AliFemtoESDTrackCut           *dtc2etaphitpc[size];
  AliFemtoKKTrackCutFull            *dtc1etaphitpc[size];
  AliFemtoESDTrackCut           *dtc2etaphitpc[size];
  AliFemtoCutMonitorParticleYPt *cutPass1YPtetaphitpc[size];
  AliFemtoCutMonitorParticleYPt *cutFail1YPtetaphitpc[size];
  AliFemtoCutMonitorParticlePIDBeta *cutPass1PIDetaphitpc[size];
  AliFemtoCutMonitorParticlePIDBeta *cutFail1PIDetaphitpc[size];
  AliFemtoCutMonitorParticleYPt *cutPass2YPtetaphitpc[size];
  AliFemtoCutMonitorParticleYPt *cutFail2YPtetaphitpc[size];
  AliFemtoCutMonitorParticlePIDBeta *cutPass2PIDetaphitpc[size];
  AliFemtoCutMonitorParticlePIDBeta *cutFail2PIDetaphitpc[size];
  AliFemtoPairCutAntiGamma      *sqpcetaphitpc[size];
  //AliFemtoShareQualityTPCEntranceSepPairCut      *sqpcetaphitpc[20];
  //AliFemtoPairCutRadialDistance      *sqpcetaphitpc[size];
  //AliFemtoPairCutRadialDistanceAsymmetric *sqpcetaphitpc[size];
  //AliFemtoPairCutMergedFraction *sqpcetaphitpcmf[size];
  //AliFemtoPairCutRadialDistanceLM      *sqpcetaphitpcRD[size];
  AliFemtoCorrFctnDirectYlm     *cylmetaphitpc[size];
  //AliFemtoCorrFctnDEtaDPhi      *cdedpetaphi[size];
  //AliFemtoChi2CorrFctn          *cchiqinvetaphitpc[size];
  AliFemtoCorrFctnDPhiStarKStarMergedFraction *dphistarkstarmftpc[size];
  AliFemtoCorrFctnDPhiStarKStarAverageMergedPointsFraction *dphistarkstarampftpc[size];
  //AliFemtoKTPairCut             *ktpcuts[size*numOfkTbins];
  AliFemtoBetaTPairCut          *ktpcuts[size*numOfkTbins];
  AliFemtoCutMonitorPairBetaT   *cutpassbetatcutmonitor[size*numOfkTbins];
  AliFemtoCutMonitorPairBetaT   *cutfailbetatcutmonitor[size*numOfkTbins];
  //AliFemtoPairCutMergedFraction *ktpcuts[size*numOfkTbins];
  AliFemtoCorrFctnDirectYlm     *cylmkttpc[size*numOfkTbins];
  //AliFemtoQinvCorrFctn          *cqinvkttpc[size*numOfkTbins];
  AliFemtoCorrFctn3DLCMSSym     *cq3dlcmskttpc[size*numOfkTbins];
  //AliFemtoCorrFctn3DSpherical   *cq3dspherical[size*numOfkTbins];
  //AliFemtoCorrFctnTPCNcls       *cqinvnclstpc[size];
  //AliFemtoChi2CorrFctn          *cqinvchi2tpc[size];

  //AliFemtoModelGausLCMSFreezeOutGenerator *gausLCMSFreezeOutGenerator[size];
  //AliFemtoModelWeightGeneratorBasic      *weightGeneratorLednicky[size];
  //AliFemtoModelManager          *tModelManager[size];
  AliFemtoCorrFctnNonIdDR       *cnonidtpc[size];
  AliFemtoCorrFctnDPhiStarDEta  *cdphistardeta08[size];
  AliFemtoCorrFctnDPhiStarDEta  *cdphistardeta12[size];
  AliFemtoCorrFctnDPhiStarDEta  *cdphistardeta16[size];
  AliFemtoCorrFctnDPhiStarDEta  *cdphistardeta20[size];
  AliFemtoCorrFctnNonIdDR       *cnonidkttpc[size*numOfkTbins];


  // *** Begin kaon-deuteron analysis ***
  int aniter = 0;  

  for (int imult = 0; imult < numOfMultBins; imult++) {
    if (runmults[imult]) {

      for (int ichg = 0; ichg < numOfChTypes; ichg++) {
	if (runch[ichg]) {
	  
	  //Iterator:
	  aniter = ichg * numOfMultBins + imult;


	  //Mix events with respect to the z position of the primary vertex and event total multipliticy:
	  anetaphitpc[aniter] = new AliFemtoVertexMultAnalysis(10, -10.0, 10.0, 4, multbins[imult], multbins[imult+1]);
	  anetaphitpc[aniter]->SetNumEventsToMix(10);
	  anetaphitpc[aniter]->SetMinSizePartCollection(1);
	  anetaphitpc[aniter]->SetVerboseMode(kFALSE);
	  
	  //Select basic cuts:
	  mecetaphitpc[aniter] = new AliFemtoBasicEventCut();
	  mecetaphitpc[aniter]->SetEventMult(0,100000);
	  mecetaphitpc[aniter]->SetVertZPos(-10,10);

	  //Study the multiplicity distribution:
	  if(turnOnMonitors == 1) {
	    cutPassEvMetaphitpc[aniter] = new AliFemtoCutMonitorEventMult(Form("cutPass%stpcM%i", chrgs[ichg], imult),500);
	    cutFailEvMetaphitpc[aniter] = new AliFemtoCutMonitorEventMult(Form("cutFail%stpcM%i", chrgs[ichg], imult),500);
	    mecetaphitpc[aniter]->AddCutMonitor(cutPassEvMetaphitpc[aniter], cutFailEvMetaphitpc[aniter]);
	  }
	  //Study the distribution and error of the primary vertex:
	  if(turnOnMonitors == 1) {
	    cutPassEvVetaphitpc[aniter] = new AliFemtoCutMonitorEventVertex(Form("cutPass%stpcM%i", chrgs[ichg], imult));
	    cutFailEvVetaphitpc[aniter] = new AliFemtoCutMonitorEventVertex(Form("cutFail%stpcM%i", chrgs[ichg], imult));
	    mecetaphitpc[aniter]->AddCutMonitor(cutPassEvVetaphitpc[aniter], cutFailEvVetaphitpc[aniter]);
	  }
	  //Study the multiplicity distribution:
	  if(turnOnMonitors == 1) {
	    cutPassColletaphitpc[aniter] = new AliFemtoCutMonitorCollections(Form("cutPass%stpcM%i", chrgs[ichg], imult));
	    cutFailColletaphitpc[aniter] = new AliFemtoCutMonitorCollections(Form("cutFail%stpcM%i", chrgs[ichg], imult));
	    mecetaphitpc[aniter]->AddCutMonitor(cutPassColletaphitpc[aniter], cutFailColletaphitpc[aniter]);
	  }



	  //Basic track cut for kaons:
	  dtc1etaphitpc[aniter] = new AliFemtoKKTrackCutFull();
	  dtc1etaphitpc[aniter]->SetNsigmaTPCle250(2.0);
	  dtc1etaphitpc[aniter]->SetNsigmaTPC250_400(2.0);
	  dtc1etaphitpc[aniter]->SetNsigmaTPC400_450(1.0);
	  dtc1etaphitpc[aniter]->SetNsigmaTPC450_500(3.0);
	  dtc1etaphitpc[aniter]->SetNsigmaTOF450_500(2.0);
	  dtc1etaphitpc[aniter]->UseNsigmaTOF450_500(true);
	  dtc1etaphitpc[aniter]->SetNsigmaTPCge500(3.0);
	  dtc1etaphitpc[aniter]->SetNsigmaTOF500_800(2.0);
	  dtc1etaphitpc[aniter]->SetNsigmaTOF800_1000(1.5);
	  dtc1etaphitpc[aniter]->SetNsigmaTOFge1000(1.0);

	  //Basic track cut for deuterons:
	  dtc2etaphitpc[aniter] = new AliFemtoESDTrackCut();
	  dtc2etaphitpc[aniter]->SetNsigmaTPCTOF(false); //no TOF Nsigma cut!
	  dtc2etaphitpc[aniter]->SetNsigma(2.0);
	  dtc2etaphitpc[aniter]->SetNsigmaMass(-1);// -1 (if no tof mass cut) or 0.7 (the +/- val of mass squared around PDG mass

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
	  dtc1etaphitpc[aniter]->SetPt(0.19,1.5);
          dtc1etaphitpc[aniter]->SetEta(-0.8,0.8);
	  dtc1etaphitpc[aniter]->SetMass(KaonMass);	  
	  dtc1etaphitpc[aniter]->SetMostProbableKaon();
	  
	  //Set particle 2:
	  dtc2etaphitpc[aniter]->SetPt(1.0,4.0);
          dtc2etaphitpc[aniter]->SetEta(-0.8,0.8);
	  dtc2etaphitpc[aniter]->SetMass(DeuteronMass);	  
	  dtc2etaphitpc[aniter]->SetMostProbableDeuteron();


	  //============KAON============

	  //The cut monitor for particles to study the difference between reconstructed and true momentum: 
	  if(turnOnMonitors == 1) { 
	    cutPass1YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutPass1%stpcM%i", chrgs[ichg], imult),KaonMass);
	    cutFail1YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutFail1%stpcM%i", chrgs[ichg], imult),KaonMass);
	    dtc1etaphitpc[aniter]->AddCutMonitor(cutPass1YPtetaphitpc[aniter], cutFail1YPtetaphitpc[aniter]);
	  }
	  //The cut monitor for particles to study various aspects of the PID determination:
	  if(turnOnMonitors == 1) {
	    cutPass1PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePIDBeta(Form("cutPass1%stpcM%i", chrgs[ichg], imult),1);//0-pion,1-kaon,2-proton, 3-deuteron
	    cutFail1PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePIDBeta(Form("cutFail1%stpcM%i", chrgs[ichg], imult),1);
	    dtc1etaphitpc[aniter]->AddCutMonitor(cutPass1PIDetaphitpc[aniter], cutFail1PIDetaphitpc[aniter]);
	  }
	  
	  //============DEUTERON============

	  if(turnOnMonitors == 1) { 
	    cutPass2YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutPass2%stpcM%i", chrgs[ichg], imult),DeuteronMass);
	    cutFail2YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutFail2%stpcM%i", chrgs[ichg], imult),DeuteronMass);
	    dtc2etaphitpc[aniter]->AddCutMonitor(cutPass2YPtetaphitpc[aniter], cutFail2YPtetaphitpc[aniter]);
	  }
	  if(turnOnMonitors == 1) {
	    cutPass2PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePIDBeta(Form("cutPass2%stpcM%i", chrgs[ichg], imult),3,-20000.0,20000.0);//0-pion,1-kaon,2-proton,3-deuteron
	    cutFail2PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePIDBeta(Form("cutFail2%stpcM%i", chrgs[ichg], imult),3,-20000.0,20000.0);
	    dtc2etaphitpc[aniter]->AddCutMonitor(cutPass2PIDetaphitpc[aniter], cutFail2PIDetaphitpc[aniter]);
	  }
	  
	  //A pair cut which checks for some pair qualities that attempt to identify slit/doubly reconstructed tracks and also selects pairs based on their separation at the entrance to the TPC
	  //sqpcetaphitpc[aniter] = new AliFemtoPairCutRadialDistance();
          //sqpcetaphitpcRD[aniter] = new AliFemtoPairCutRadialDistanceLM();

	  if (paircutmergedfractionon == 1) 
	    sqpcetaphitpc[aniter] = new AliFemtoPairCutMergedFraction(distance, fraction1, 0.01, 0.8, 2.5);
	  else
	    sqpcetaphitpc[aniter] = new AliFemtoPairCutAntiGamma();
	  
	  if(paircutantigammaon == 1)
	    sqpcetaphitpc[aniter]->SetDataType(AliFemtoPairCut::kAOD);
	  else
	    sqpcetaphitpc[aniter]->SetDataType(AliFemtoPairCut::kKine);

	  if(gammacut == 1) {
	    sqpcetaphitpc[aniter]->SetMaxEEMinv(0.002);
	    sqpcetaphitpc[aniter]->SetMaxThetaDiff(0.008);
	  }

	  if( runtype==0){
	    sqpcetaphitpc[aniter]->SetMaxEEMinv(0.0);
	    sqpcetaphitpc[aniter]->SetMaxThetaDiff(0.0);
	    sqpcetaphitpc[aniter]->SetTPCEntranceSepMinimum(1.5);
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
	  }
	  
	  // Study the betaT distribution:
	  if(turnOnBetaTMonitor == 1) {
	  cutpassbetatcutmonitor[aniter] = new AliFemtoCutMonitorPairBetaT(Form("cutPass%stpcM%i", chrgs[ichg], imult), 100, 0.0, 1.0, KaonMass, DeuteronMass);
	  cutfailbetatcutmonitor[aniter] = new AliFemtoCutMonitorPairBetaT(Form("cutFail%stpcM%i", chrgs[ichg], imult), 100, 0.0, 1.0, KaonMass, DeuteronMass);
	  sqpcetaphitpc[aniter]->AddCutMonitor(cutpassbetatcutmonitor[aniter], cutfailbetatcutmonitor[aniter]);
	  }
	  
	  anetaphitpc[aniter]->SetEnablePairMonitors(true);
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

	  // NonId (without kT bins)
	  // Correlation function for non-identical particles uses k* as a function variable. Stores the correlation 
	  // function separately for positive and negative signs of k* projections into out, side and long directions, 
	  // enabling the calculations of double ratios
	  if(runNonIdCorrFctn == 1) {
	    cnonidtpc[aniter] = new AliFemtoCorrFctnNonIdDR(Form("cnonid%stpcM%i", chrgs[ichg], imult), nbinssh, 0.0, shqmax);
	    anetaphitpc[aniter]->AddCorrFctn(cnonidtpc[aniter]);
	  }

	  if(runDPhiStarKStarMergedFraction == 1) {
	    dphistarkstarmftpc[aniter] = new AliFemtoCorrFctnDPhiStarKStarMergedFraction(Form("cdphistarkstarmergedfraction%stpcM%iD%lfF%lf", chrgs[ichg], imult, distance, fraction1), 0.8, 2.5, distance, fraction1, 0.01, 51, 0.0, 0.5, 127, -1.0, 1.0);
	    anetaphitpc[aniter]->AddCorrFctn(dphistarkstarmftpc[aniter]);
	  }
	  
	  if(runDPhiStarKStarAverageMergedPointsFraction == 1) {
	    dphistarkstarampftpc[aniter] = new AliFemtoCorrFctnDPhiStarKStarAverageMergedPointsFraction(Form("cdphistarkstaraveragemergedpointsfraction%stpcM%iD%lf", chrgs[ichg], imult, distance), 0.8, 2.5, distance, 0.01, 51, 0.0, 0.5, 127, -0.5, 0.5);
	    anetaphitpc[aniter]->AddCorrFctn(dphistarkstarampftpc[aniter]);
	  }
	  // DPhiStarDEta (without bins)
	  // Correlation function for two particle correlations which uses dPhi* and dEta as a function variables.
	  if(runDPhiStarDEta == 1) {
	  cdphistardeta08[aniter] = new AliFemtoCorrFctnDPhiStarDEta(Form("cdphistardeta08%stpcM%i", chrgs[ichg], imult), 0.8, 51, -0.05, 0.05, 127, -1.0, 1.0);
	  anetaphitpc[aniter]->AddCorrFctn(cdphistardeta08[aniter]);
	  cdphistardeta12[aniter] = new AliFemtoCorrFctnDPhiStarDEta(Form("cdphistardeta12%stpcM%i", chrgs[ichg], imult), 1.2, 51, -0.05, 0.05, 127, -1.0, 1.0);
	  anetaphitpc[aniter]->AddCorrFctn(cdphistardeta12[aniter]);
	  cdphistardeta16[aniter] = new AliFemtoCorrFctnDPhiStarDEta(Form("cdphistardeta16%stpcM%i", chrgs[ichg], imult), 1.6, 51, -0.05, 0.05, 127, -1.0, 1.0);
	  anetaphitpc[aniter]->AddCorrFctn(cdphistardeta16[aniter]);
	  cdphistardeta20[aniter] = new AliFemtoCorrFctnDPhiStarDEta(Form("cdphistardeta20%stpcM%i", chrgs[ichg], imult), 2.0, 51, -0.05, 0.05, 127, -1.0, 1.0);
	  anetaphitpc[aniter]->AddCorrFctn(cdphistardeta20[aniter]);
	  }

	  if (runbetatdep) {
	    int ktm;
	    for (int ikt=0; ikt<numOfkTbins; ikt++) {
	      ktm = aniter*numOfkTbins + ikt;
	      //ktpcuts[ktm] = new AliFemtoKTPairCut(ktrng[ikt], ktrng[ikt+1]);
	      ktpcuts[ktm] = new AliFemtoBetaTPairCut(ktrng[ikt], ktrng[ikt + 1], KaonMass, DeuteronMass);
	      
	      if (runbetatylm) {
		cylmkttpc[ktm] = new AliFemtoCorrFctnDirectYlm(Form("cylm%stpcM%iD%lfF%lfbetat%d", chrgs[ichg], imult, distance, fraction1, ikt),3,nbinssh, 0.0,shqmax, runshlcms);
		cylmkttpc[ktm]->SetPairSelectionCut(ktpcuts[ktm]);
		anetaphitpc[aniter]->AddCorrFctn(cylmkttpc[ktm]);
	      }
	      
	      if (runbetatnonid) {
		cnonidkttpc[ktm] = new AliFemtoCorrFctnNonIdDR(Form("cnonid%stpcM%ibetaT%i", chrgs[ichg], imult,ikt), nbinssh, 0.0, shqmax);
		cnonidkttpc[ktm]->SetPairSelectionCut(ktpcuts[ktm]);
		anetaphitpc[aniter]->AddCorrFctn(cnonidkttpc[ktm]);
	      }
	      
	    }
	  }
	  
	
	  
	  Manager->AddAnalysis(anetaphitpc[aniter]);	
	}
      }
    }
  }
  // *** End kaon-deuteron analysis

  return Manager;
}
