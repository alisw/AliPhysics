 
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
#include "AliFemtoSphericityEventCut.h"
#include "AliFemtoSpherocityEventCut.h"
#include "AliFemtoBasicEventCut.h" 
#include "AliFemtoESDTrackCut.h" 
#include "AliFemtoESDTrackCutPlusJets25.h" 
//#include "AliFemtoKKTrackCut.h"
#include "AliFemtoKpm45TrackCut.h"
#include "AliFemtoCorrFctn.h"
#include "AliFemtoCutMonitorParticleYPt.h"
#include "AliFemtoCutMonitorParticleVertPos.h"
#include "AliFemtoCutMonitorParticleMomRes.h"
#include "AliFemtoCutMonitorParticlePID.h"
#include "AliFemtoCutMonitorEventMult.h"
#include "AliFemtoCutMonitorEventVertex.h"
#include "AliFemtoCutMonitorEventSphericity.h"
#include "AliFemtoCutMonitorDphiDeta1.h"
#include "AliFemtoShareQualityTPCEntranceSepPairCut.h"
#include "AliFemtoPairCutAntiGamma.h"
#include "AliFemtoPairCutRadialDistance.h"
#include "AliFemtoPairCutRadialDistanceKK.h"
#include "AliFemtoPairCutRadialDistanceKKdist.h"
#include "AliFemtoQinvCorrFctn.h"
#include "AliFemtoShareQualityCorrFctn.h"
#include "AliFemtoTPCInnerCorrFctn.h"
#include "AliFemtoVertexMultAnalysis.h"
#include "AliFemtoCorrFctn3DSpherical.h"
#include "AliFemtoChi2CorrFctn.h"
#include "AliFemtoCorrFctnTPCNcls.h"
#include "AliFemtoBPLCMS3DCorrFctn.h"
#include "AliFemtoBPLCMS3DCorrFctnKK.h"
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
#include "AliFemtoKTPairCut.h"
#include "AliFemtoCutMonitorParticlePtPDG.h";
#include "AliFemtoAvgSepCorrFctn.h";
//#include "AliFemtoCorrFctnDEtaDPhiKK2.h";
#include "AliFemtoCorrFctnDEtaDPhi.h";


#endif

//________________________________________________________________________
AliFemtoManager* ConfigFemtoAnalysis() {

  double PionMass = 0.13956995;
  double KaonMass = 0.493677;
 ///// const int cMu=4;
 //////// const int cKt=3;
  
  
  const int cMu=3;
  const int cKt=7; //8
  
  
  

  //-------Single track cuts------------------------------------------------->
  double DCAxy=0.3;//2.4;// cm // our standard is 0.20 cm; super narrow was 0.015cm
  double DCAz =0.3;//3.0;// cm // our standard is 0.15 cm;
  //-------Single track cuts-------------------------------------------------<
  //=======Double track cuts=================================================>
  //Dhevan's : PhiStarDifferenceMinimum=0.06; EtaDifferenceMinimum=0.02;
  //standart
  //double PhiStarDifferenceMinimum=0.017; //[radian]
 // double EtaDifferenceMinimum=0.015; //[radian]
 //for test
   //double PhiStarDifferenceMinimum=0.03; //[radian]
 // double EtaDifferenceMinimum=0.02; //[radian]
  // double PhiStarDifferenceMinimum=0.04; //[radian]
  // double EtaDifferenceMinimum=0.02; //[radian]
  
  
  
  //=======Double track cuts=================================================<

  // Switches for QA analyses
 
  int runmults[4] = {1, 1, 1, 0};
 // int multbins[5] = {0, 900, 300, 500, 900};
//old pp
   int multbins[5] = {1, 18, 30, 100, 1000};
  //int multbins[5] = {1, 100000, 30, 100, 1000};
  //.................................................

  int runch[2] = {1, 1};
  const char *chrgs[2] = { "Kp", "Km"};
  
  
  int runktdep = 1;

 //double ktrng[3] = {0.15, 0.5, 1.2};
  //double ktrng[3] = {0.15, 6.0, 1.3};
//double ktrng[3] = {2.0, 4.0, 2.0};
 

double ktrng[8] = {0.15, 0.3, 0.5, 0.7, 0.9, 1.2, 3.0, 5.0};
 
 // double ktrng[3] = {0.15, 1.2};

// double ktrng[8] = {0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 2.0};

  int run3d = 0;
  int runshlcms = 0;

  int runtype = 0; // Types 0 - global, 1 - ITS only, 2 - TPC Inner
  int isrealdata = 1;

  //   AliFemtoEventReaderESDChainKine* Reader=new AliFemtoEventReaderESDChainKine();
  //   Reader->SetConstrained(true);
  //   Reader->SetUseTPCOnly(false);

  double shqmax;
  double shqmaxSH;
  int nbinssh = 100;
  int nbinsh3D = 160;

//ml  if (runshlcms) shqmax = 0.25;
//  else shqmax = 0.9;
  

  if (runshlcms) shqmaxSH = 0.25;
  shqmax = 0.9;
  

 

AliFemtoEventReaderAODChain *Reader = new AliFemtoEventReaderAODChain();
//AliFemtoEventReaderAOD *Reader = new AliFemtoEventReaderAODMultSelection();


   Reader->SetFilterMask(96);
  // Reader->SetFilterMask(128);
 //  Reader->SetFilterMask(768);
   
   
//    Reader->SetFilterBit(128);
    
     Reader->SetUseMultiplicity(AliFemtoEventReaderAOD::kRefComb08);
    
    Reader->SetCentralityFlattening(kFALSE);
    Reader->SetReadV0(0);
    // rdr->SetPrimaryVertexCorrectionTPCPoints(kTRUE);
    Reader->SetDCAglobalTrack(kTRUE);

    
//    Reader->SetReadMC(kTRUE);
    //  Reader->SetPionAnalysis(kTRUE);

    
 
//    Reader->SetUseMVPlpSelection(kTRUE);
//    Reader->SetTrackPileUpRemoval(kTRUE);
//    Reader->SetUseAliEventCuts(kTRUE);
    
    
    
  AliFemtoManager* Manager=new AliFemtoManager();
  Manager->SetEventReader(Reader);

  AliFemtoVertexMultAnalysis    *anetaphitpc[20];
 ///////// AliFemtoBasicEventCut         *mecetaphitpc[20];
  AliFemtoSphericityEventCut    *mecetaphitpc[20];
 // AliFemtoSpherocityEventCut    *mecetaphitpc[20];

  AliFemtoCutMonitorEventMult   *cutPassEvMetaphitpc[20];
  AliFemtoCutMonitorEventMult   *cutFailEvMetaphitpc[20];
  // AliFemtoCutMonitorEventVertex *cutPassEvVetaphitpc[20];
  // AliFemtoCutMonitorEventVertex *cutFailEvVetaphitpc[20];
  AliFemtoCutMonitorEventSphericity  *cutPassEvSpher[20];
  AliFemtoCutMonitorEventSphericity  *cutFailEvSpher[20];

  AliFemtoCutMonitorDphiDeta1  *cutPassDphiDeta[20];
  AliFemtoCutMonitorDphiDeta1  *cutFailDphiDeta[20];


  AliFemtoCutMonitorCollections   *cutPassColletaphitpc[20];
  AliFemtoCutMonitorCollections   *cutFailColletaphitpc[20];
 

 // AliFemtoKpm45TrackCut           *dtc1etaphitpc[20];
 // AliFemtoKpm45TrackCut           *dtc2etaphitpc[20];

 // AliFemtoESDTrackCut           *dtc1etaphitpc[20];
 // AliFemtoESDTrackCut           *dtc2etaphitpc[20];
  

  AliFemtoESDTrackCutPlusJets25           *dtc1etaphitpc[20];
  AliFemtoESDTrackCutPlusJets25           *dtc2etaphitpc[20];
  
  
  
//    AliFemtoKKTrackCutTest           *dtc1etaphitpc[20];
//    AliFemtoKKTrackCutTest           *dtc2etaphitpc[20];


  AliFemtoCutMonitorParticleYPt *cutPass1YPtetaphitpc[20];
  AliFemtoCutMonitorParticleYPt *cutFail1YPtetaphitpc[20];
 
  AliFemtoCutMonitorParticlePID *cutPass1PIDetaphitpc[20];
  AliFemtoCutMonitorParticlePID *cutFail1PIDetaphitpc[20];
   
 // AliFemtoCutMonitorParticlePtPDG *cutPass1PIDetaphitpc[20];
 // AliFemtoCutMonitorParticlePtPDG *cutFail1PIDetaphitpc[20];
  
  
  AliFemtoCutMonitorParticleYPt *cutPass2YPtetaphitpc[20];
  AliFemtoCutMonitorParticleYPt *cutFail2YPtetaphitpc[20];
  AliFemtoCutMonitorParticlePID *cutPass2PIDetaphitpc[20];
  AliFemtoCutMonitorParticlePID *cutFail2PIDetaphitpc[20];


//AliFemtoPairCutAntiGamma      *sqpcetaphitpc[20];
//    AliFemtoShareQualityTPCEntranceSepPairCut      *sqpcetaphitpc[20];
 // AliFemtoPairCutRadialDistance      *sqpcetaphitpc[20];//AliFemto dphi* cut
 //// AliFemtoPairCutRadialDistanceKK      *sqpcetaphitpc[20];//Dhevan's dphi* cut
  AliFemtoPairCutRadialDistanceKKdist      *sqpcetaphitpc[20];//Dhevan's dphi* cut
 
  AliFemtoCorrFctnDirectYlm     *cylmetaphitpc[20];
  //AliFemtoCorrFctnDEtaDPhiKK2      *cdedpetaphi[20*10];//20->20*10 due to kT
  AliFemtoCorrFctnDEtaDPhi     *cdedpetaphi[20*10];//20->20*10 due to kT
  
  AliFemtoChi2CorrFctn          *cchiqinvetaphitpc[20];
  AliFemtoKTPairCut             *ktpcuts[20*8];
  AliFemtoCorrFctnDirectYlm     *cylmkttpc[20*8];
  AliFemtoQinvCorrFctn          *cqinvkttpc[20*8];
 // AliFemtoCorrFctn3DLCMSSym     *cq3dlcmskttpc[20*8];
  AliFemtoBPLCMS3DCorrFctnKK  *cq3dlcmskttpc[20*8];
 // AliFemtoBPLCMS3DCorrFctn  *cq3dlcmskttpc[20*8];
 
  AliFemtoCorrFctnTPCNcls       *cqinvnclstpc[20];
  AliFemtoShareQualityCorrFctn  *cqinvsqtpc[20*10];
  AliFemtoChi2CorrFctn          *cqinvchi2tpc[20];
  //AliFemtoTPCInnerCorrFctn      *cqinvinnertpc[20*10];
  AliFemtoAvgSepCorrFctn  *cqinvinnertpc[20*10];
  

  // *** Third QA task - HBT analysis with all pair cuts off, TPC only ***
  // *** Begin Kaon-Kaon (positive) analysis ***
  int aniter = 0;

  bool verbose=false;

  for (int imult=0; imult<cMu/*4*/; imult++) {
    if (runmults[imult]) {
      for (int ichg=0; ichg<2; ichg++) {
	if (runch[ichg]) {
	  aniter = ichg*cMu+imult; //0, 1(ich=0) ,2,3


	  anetaphitpc[aniter] = new AliFemtoVertexMultAnalysis(10, -10.0, 10.0, 8, multbins[imult], multbins[imult+1]);
	  anetaphitpc[aniter]->SetNumEventsToMix(10);
	  anetaphitpc[aniter]->SetMinSizePartCollection(1);
	  anetaphitpc[aniter]->SetVerboseMode(verbose);

/*
	  mecetaphitpc[aniter] = new AliFemtoBasicEventCut();
	  mecetaphitpc[aniter]->SetEventMult(0.01,100000);
//	  mecetaphitpc[aniter]->SetVertZPos(-10.0,10.0);
	  mecetaphitpc[aniter]->SetVertZPos(-10.0,10.0);
	  */

	  
	  mecetaphitpc[aniter] = new AliFemtoSphericityEventCut();
	    //	  mecetaphitpc[aniter] = new AliFemtoSpherocityEventCut();
          mecetaphitpc[aniter]->SetEventMult(1,100000);
          mecetaphitpc[aniter]->SetVertZPos(-10.0,10.0);
          mecetaphitpc[aniter]->SetStMin(0.0);
          mecetaphitpc[aniter]->SetStMax(0.05);

	  /* //was in aliroot 5.03.76
	  if (isrealdata)
	     mecetaphitpc[aniter]->SetAcceptOnlyPhysics(kTRUE);
	  */
	  //    mecetaphitpc->SetAcceptBadVertex(kTRUE);
	  
	  cutPassEvMetaphitpc[aniter] = new AliFemtoCutMonitorEventMult(Form("cutPass%stpcM%i", chrgs[ichg], imult));
	  cutFailEvMetaphitpc[aniter] = new AliFemtoCutMonitorEventMult(Form("cutFail%stpcM%i", chrgs[ichg], imult));
	  mecetaphitpc[aniter]->AddCutMonitor(cutPassEvMetaphitpc[aniter], cutFailEvMetaphitpc[aniter]);
	  
	  //	  cutPassEvVetaphitpc[aniter] = new AliFemtoCutMonitorEventVertex(Form("cutPass%stpcM%i", chrgs[ichg], imult));
	  //	  cutFailEvVetaphitpc[aniter] = new AliFemtoCutMonitorEventVertex(Form("cutFail%stpcM%i", chrgs[ichg], imult));
	  //	  mecetaphitpc[aniter]->AddCutMonitor(cutPassEvVetaphitpc[aniter], cutFailEvVetaphitpc[aniter]);
	  
	  

	   cutPassEvSpher[aniter] = new AliFemtoCutMonitorEventSphericity(Form("cutPass%stpcM%i", chrgs[ichg], imult));
	   cutFailEvSpher[aniter] = new AliFemtoCutMonitorEventSphericity(Form("cutFail%stpcM%i", chrgs[ichg], imult));
	   mecetaphitpc[aniter]->AddCutMonitor(cutPassEvSpher[aniter], cutFailEvSpher[aniter]);
	  
	  
	   cutPassDphiDeta[aniter] = new AliFemtoCutMonitorDphiDeta1(Form("cutPass%stpcM%i", chrgs[ichg], imult));
	   cutFailDphiDeta[aniter] = new AliFemtoCutMonitorDphiDeta1(Form("cutFail%stpcM%i", chrgs[ichg], imult));
	   mecetaphitpc[aniter]->AddCutMonitor(cutPassDphiDeta[aniter], cutFailDphiDeta[aniter]);	  
	  
	  
	  //Study the collection multiplicity distribution
//	  cutPassColletaphitpc[aniter] = new AliFemtoCutMonitorCollections(Form("cutPass%stpcM%i", chrgs[ichg], imult));
//	  cutFailColletaphitpc[aniter] = new AliFemtoCutMonitorCollections(Form("cutFail%stpcM%i", chrgs[ichg], imult));
//	  mecetaphitpc[aniter]->AddCutMonitor(cutPassColletaphitpc[aniter], cutFailColletaphitpc[aniter]);
	  
	  
	  
//	  dtc1etaphitpc[aniter] = new AliFemtoKpm45TrackCut()
//           dtc1etaphitpc[aniter] = new AliFemtoKKTrackCutTest();

	   // dtc1etaphitpc[aniter] = new AliFemtoESDTrackCut();
	   dtc1etaphitpc[aniter] = new AliFemtoESDTrackCutPlusJets25();
           
           dtc1etaphitpc[aniter]->SetNsigmaTPCTOF(true);
	   dtc1etaphitpc[aniter]->AddCutMonitor(cutPassDphiDeta[aniter]);
           dtc1etaphitpc[aniter]->SetNsigma(3.0);
           
           
	  //     dtc1etaphitpc[aniter]->SetPidProbPion(0.0,1.001);
	  //     dtc1etaphitpc[aniter]->SetPidProbMuon(0.0,1.0);
	  //     dtc1etaphitpc[aniter]->SetPidProbKaon(0.0,1.0);
	  //     dtc1etaphitpc[aniter]->SetPidProbProton(0.0,1.0);
	  if (ichg == 0)
	    dtc1etaphitpc[aniter]->SetCharge(1.0);
	  else if (ichg == 1)
	    dtc1etaphitpc[aniter]->SetCharge(-1.0);
	    
	  dtc1etaphitpc[aniter]->SetPt(0.15,5.0);
	  //	  dtc1etaphitpc[aniter]->SetEta(-1.2,1.2);
	  dtc1etaphitpc[aniter]->SetEta(-0.8,0.8); //0.5
	  // 	//    dtc1etaphitpc[aniter]->SetEta(-0.5,0.5);
	  dtc1etaphitpc[aniter]->SetMass(PionMass);
//	  dtc1etaphitpc[aniter]->SetMass(KaonMass);
	   

	  dtc1etaphitpc[aniter]->SetMostProbablePion();  //!!!!!!

	  	 
	  // 	// Track quality cuts
	  if (runtype == 0) {
	    /*
	    dtc1etaphitpc[aniter]->SetStatus(AliESDtrack::kTPCrefit|AliESDtrack::kITSrefit);
	    //	    dtc1etaphitpc[aniter]->SetStatus(AliESDtrack::kTPCrefit);
	    //    dtc1etaphitpc[aniter]->SetStatus(AliESDtrack::kITSrefit);
	    dtc1etaphitpc[aniter]->SetminTPCncls(80);
	    dtc1etaphitpc[aniter]->SetRemoveKinks(kTRUE);
	    dtc1etaphitpc[aniter]->SetLabel(kFALSE);
	    //    dtc1etaphitpc[aniter]->SetMaxITSChiNdof(6.0);
	    dtc1etaphitpc[aniter]->SetMaxTPCChiNdof(4.0);
	    //Poland: 
            dtc1etaphitpc[aniter]->SetMaxImpactXYPtDep(0.0182, 0.0350, -1.01);
	    //      dtc1etaphitpc[aniter]->SetMaxSigmaToVertex(6.0);
	    */
           dtc1etaphitpc[aniter]->SetMaxImpactXY(DCAxy);
           dtc1etaphitpc[aniter]->SetMaxImpactZ(DCAz);

	  }
	  else if (runtype == 1) {

	  }
	  else if (runtype == 2) {


	  }

	  

	  cutPass1YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutPass1%stpcM%i", chrgs[ichg], imult), 0.139569);
	  cutFail1YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutFail1%stpcM%i", chrgs[ichg], imult), 0.139569);
	  dtc1etaphitpc[aniter]->AddCutMonitor(cutPass1YPtetaphitpc[aniter], cutFail1YPtetaphitpc[aniter]);
	  
	  cutPass1PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cutPass1%stpcM%i", chrgs[ichg], imult),0);
	  cutFail1PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cutFail1%stpcM%i", chrgs[ichg], imult),0);
	  dtc1etaphitpc[aniter]->AddCutMonitor(cutPass1PIDetaphitpc[aniter], cutFail1PIDetaphitpc[aniter]);

//	  cutPass1PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePtPDG(Form("cutPass1%stpcM%i", chrgs[ichg], imult),0);
//	  cutFail1PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePtPDG(Form("cutFail1%stpcM%i", chrgs[ichg], imult),0);
//	  dtc1etaphitpc[aniter]->AddCutMonitor(cutPass1PIDetaphitpc[aniter], cutFail1PIDetaphitpc[aniter]);

	  
	 // sqpcetaphitpc[aniter] = new AliFemtoPairCutAntiGamma();
	//  sqpcetaphitpc[aniter] = new AliFemtoShareQualityTPCEntranceSepPairCut();
	  
          if (ichg < 2) {
	    sqpcetaphitpc[aniter] = new AliFemtoPairCutRadialDistanceKKdist();//AliFemto dphi* cut
//    	  sqpcetaphitpc[aniter] = new AliFemtoPairCutRadialDistanceKK();  //Dhevan's dphi* cut
	  if (runtype == 0) {
	    sqpcetaphitpc[aniter]->SetShareQualityMax(1.0);
	    sqpcetaphitpc[aniter]->SetShareFractionMax(0.05);
	    sqpcetaphitpc[aniter]->SetRemoveSameLabel(kFALSE);
            sqpcetaphitpc[aniter]->SetAverageSeparation(0.0); //0.8
	  }
	  else if (runtype == 1) {
	    sqpcetaphitpc[aniter]->SetShareQualityMax(1.0);
	    sqpcetaphitpc[aniter]->SetShareFractionMax(1.05);
	    sqpcetaphitpc[aniter]->SetRemoveSameLabel(kFALSE);
            sqpcetaphitpc[aniter]->SetAverageSeparation(0.0); //0.8
	  }
	  else if (runtype == 2) {
            sqpcetaphitpc[aniter]->SetDataType(AliFemtoPairCut::kAOD);
	    sqpcetaphitpc[aniter]->SetShareQualityMax(1.0);
	    sqpcetaphitpc[aniter]->SetShareFractionMax(0.05);
	    sqpcetaphitpc[aniter]->SetRemoveSameLabel(kFALSE);
            sqpcetaphitpc[aniter]->SetAverageSeparation(0.0); //0.8
	  }
         }        

	  
	  anetaphitpc[aniter]->SetEventCut(mecetaphitpc[aniter]);
	  anetaphitpc[aniter]->SetFirstParticleCut(dtc1etaphitpc[aniter]);
	  anetaphitpc[aniter]->SetSecondParticleCut(dtc1etaphitpc[aniter]);
	  anetaphitpc[aniter]->SetPairCut(sqpcetaphitpc[aniter]);
	  
	


	  if (runktdep) {
	    int ktm;
	    for (int ikt=0; ikt<cKt/*8*/; ikt++) {
	      ktm = aniter*cKt/*8*/ + ikt;
	      ktpcuts[ktm] = new AliFemtoKTPairCut(ktrng[ikt], ktrng[ikt+1]);
	      
//              cqinvkttpc[ktm] = new AliFemtoQinvCorrFctn(Form("cqinv%stpcM%ikT%i", chrgs[ichg], imult, ikt),nbinssh,0.0,(imult>6)?shqmax*2.5:shqmax);

	      cqinvkttpc[ktm] = new AliFemtoQinvCorrFctn(Form("cqinv%stpcM%ikT%i", chrgs[ichg], imult, ikt),nbinssh,0.0,2.0);
	      cqinvkttpc[ktm]->SetPairSelectionCut(ktpcuts[ktm]);
	      anetaphitpc[aniter]->AddCorrFctn(cqinvkttpc[ktm]);
	      
	      
	      
        //      cdedpetaphi[ktm] = new AliFemtoCorrFctnDEtaDPhi(Form("cdetadphi%stpcM%ikT%i", chrgs[ichg], imult, ikt), 45,45);
	//      cdedpetaphi[ktm]->SetPairSelectionCut(ktpcuts[ktm]);
        //      anetaphitpc[aniter]->AddCorrFctn(cdedpetaphi[ktm]);
	      
	      

	      cqinvsqtpc[ktm] = new AliFemtoShareQualityCorrFctn(Form("cqinvsq%stpcM%ikT%i", chrgs[ichg], imult, ikt),nbinssh,0.0,shqmax);
	      cqinvsqtpc[ktm]->SetPairSelectionCut(ktpcuts[ktm]);
	      anetaphitpc[aniter]->AddCorrFctn(cqinvsqtpc[ktm]);

/*
	      cqinvinnertpc[ktm] = new AliFemtoTPCInnerCorrFctn(Form("cqinvinner%stpcM%ikT%i", chrgs[ichg], imult, ikt),nbinssh,0.0,shqmax);
	      cqinvinnertpc[ktm]->SetPairSelectionCut(ktpcuts[ktm]);
	      cqinvinnertpc[ktm]->SetRadius(1.2);
	      anetaphitpc[aniter]->AddCorrFctn(cqinvinnertpc[ktm]);
*/	      
	      //cqinvinnertpc[ktm] = new AliFemtoAvgSepCorrFctn(0,50,0.0,100.);
              //anetaphitpc[aniter]->AddCorrFctn(cqinvinnertpc[ktm]);
              
              
              


	      if (run3d) {
	//	cq3dlcmskttpc[ktm] = new AliFemtoCorrFctn3DLCMSSym(Form("cq3d%stpcM%ikT%i", chrgs[ichg], imult, ikt),60,(imult>3)?((imult>6)?((imult>7)?0.6:0.4):0.25):0.15);

        //    AliFemtoBPLCMS3DCorrFctn *cq3dallpiptpc = new AliFemtoBPLCMS3DCorrFctn("cq3dallpiptpc",100,-1.5,1.5);
                                 
//our
	cq3dlcmskttpc[ktm] = new AliFemtoBPLCMS3DCorrFctnKK(Form("cq3d%stpcM%ikT%i", chrgs[ichg], imult, ikt),nbinsh3D,-0.3,0.3);
//	cq3dlcmskttpc[ktm] = new AliFemtoBPLCMS3DCorrFctn(Form("cq3d%stpcM%ikT%i", chrgs[ichg], imult, ikt),nbinsh3D,-0.3,0.3);
//	cq3dlcmskttpc[ktm] = new AliFemtoCorrFctn3DLCMSSym(Form("cq3d%stpcM%ikT%i", chrgs[ichg], imult, ikt),nbinsh3D,0.5);
		cq3dlcmskttpc[ktm]->SetPairSelectionCut(ktpcuts[ktm]);
		anetaphitpc[aniter]->AddCorrFctn(cq3dlcmskttpc[ktm]);
	      }
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
                      
