// Konstantin Mikhaylov: 20-DEC-2024
// 012: PWGCF/FEMTOSCOPY/macros/Train/KM/pPb502/011/ConfigFemtoAnalysis.C
//          ____  ::: KM 012->KK 13_DEC_2024
// Test for K+K+ to select spherical evens 0.7<S<1
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
//#include "AliFemtoESDTrackCut.h" //pion-pion
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
#include "AliFemtoPairCutRadialDistanceKKdist.h" //KK
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
  //const char *chrgs[2] = { "pip", "pim"};
  const char *chrgs[2] = { "Kp", "Km"};

  int runktdep = 1;
  const int cktbins=2; //corresponds to  number of kT ranges
  double ktrng[cktbins+1] = {0.2, 0.5, 1.0};
  //double ktrng[9] = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.3};

  int run3d = 1;//0;
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
  ///////////////////////////////////////////////////////////////////////////////////////////////
  //AliFemtoBasicEventCut         *mecetaphitpc[20];//COMMENT this line in case of sphericity!!!!
  AliFemtoSphericityEventCut    *mecetaphitpc[20];//sphericity->()
  //AliFemtoSpherocityEventCut         *mecetaphitpc[20];
  AliFemtoCutMonitorEventMult   *cutPassEvMetaphitpc[20];
  AliFemtoCutMonitorEventMult   *cutFailEvMetaphitpc[20];
  AliFemtoCutMonitorEventVertex *cutPassEvVetaphitpc[20];
  AliFemtoCutMonitorEventVertex *cutFailEvVetaphitpc[20];
  AliFemtoCutMonitorCollections   *cutPassColletaphitpc[20];
  AliFemtoCutMonitorCollections   *cutFailColletaphitpc[20];
  //AliFemtoKKTrackCut           *dtc1etaphitpc[20];
  //AliFemtoKKTrackCut           *dtc2etaphitpc[20];
     AliFemtoKpm45TrackCut           *dtc1etaphitpc[20];//K+K- or KK(++,--)
     AliFemtoKpm45TrackCut           *dtc2etaphitpc[20];//K+K-
  //AliFemtoESDTrackCut           *dtc1etaphitpc[20];//pipi
  //AliFemtoESDTrackCut           *dtc2etaphitpc[20];//pipi
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
  AliFemtoCutMonitorEventSphericity *cutPassEvSpher[20];////sphericity->()??KM
  AliFemtoCutMonitorEventSphericity *cutFailEvSpher[20];////sphericity->()??KM
  //AliFemtoPairCutAntiGamma      *sqpcetaphitpc[20];
  //AliFemtoShareQualityPairCut      *sqpcetaphitpc[20];//pipi not for K+K-
  //AliFemtoPairCutRadialDistance      *sqpcetaphitpc[20];
  AliFemtoPairCutRadialDistanceKKdist  *sqpcetaphitpc[20];//K+K- of KK(++,--)
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
      for (int ichg=0; ichg<1/*2*/; ichg++) {//'1' for K+K- or K+K+ ; 2 KK(++,--)
	if (runch[ichg]) {
	  aniter = ichg*3+imult;

	  anetaphitpc[aniter] = new AliFemtoVertexMultAnalysis(10, -10.0, 10.0, 4, multbins[imult], multbins[imult+1]);
	  anetaphitpc[aniter]->SetNumEventsToMix(3);
	  anetaphitpc[aniter]->SetMinSizePartCollection(1);
	  anetaphitpc[aniter]->SetVerboseMode(kFALSE);
      
	  mecetaphitpc[aniter] = new AliFemtoBasicEventCut();
	  mecetaphitpc[aniter]->SetEventMult(0,10000);
	  mecetaphitpc[aniter]->SetVertZPos(-10,10);
	  //sphericity->() ===>
	  mecetaphitpc[aniter] = new AliFemtoSphericityEventCut();
          mecetaphitpc[aniter]->SetEventMult(0.01,100000);
          mecetaphitpc[aniter]->SetVertZPos(-10.0,10.0);
          mecetaphitpc[aniter]->SetStMin(0.7);//spherical 
          mecetaphitpc[aniter]->SetStMax(1.0);//events 0.7<S<1
	  //sphericity->() ===<
	  
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
	  
          //sphericity->() ===>
	  cutPassEvSpher[aniter] = new AliFemtoCutMonitorEventSphericity(Form("cutPass%stpcM%i", chrgs[ichg], imult));
          cutFailEvSpher[aniter] = new AliFemtoCutMonitorEventSphericity(Form("cutFail%stpcM%i", chrgs[ichg], imult));
          mecetaphitpc[aniter]->AddCutMonitor(cutPassEvSpher[aniter], cutFailEvSpher[aniter]);
	  //sphericity->() ===<

	  //PID for kaons or pions --->
	  //dtc1etaphitpc[aniter] = new AliFemtoESDTrackCut();//pipi(++,--)
	  dtc1etaphitpc[aniter] = new AliFemtoKpm45TrackCut();//KK(++,--)
	  if (ichg == 0)
	    dtc1etaphitpc[aniter]->SetCharge(1.0);
	  else if (ichg == 1)
	  dtc1etaphitpc[aniter]->SetCharge(-1.0);
	  //p_T, |eta| ---------------------------->
	  dtc1etaphitpc[aniter]->SetPt(0.14/*0.13*/,1.5);//0.14 from Elena ??
	  dtc1etaphitpc[aniter]->SetEta(-0.8,0.8);
	  //p_T, |eta| ----------------------------<
	  //Kaon or pion 
	  dtc1etaphitpc[aniter]->SetMass(KaonMass);//->SetMass(PionMass);
	  dtc1etaphitpc[aniter]->SetMostProbableKaon();//->SetMostProbablePion();
	  //----Kch special PID since 2013 run1 --------->
	  dtc1etaphitpc[aniter]->SetNsigmaTPCle250(2.0);
	  dtc1etaphitpc[aniter]->SetNsigmaTPC250_400(2.0);
	  dtc1etaphitpc[aniter]->SetNsigmaTPC400_450(1.0);//cut on e+e- orig(2.0);
	  dtc1etaphitpc[aniter]->SetNsigmaTPC450_500(2.0);//cut on e+e- orig(2.0);
	  dtc1etaphitpc[aniter]->SetNsigmaTPCge500(3.0);  
	  dtc1etaphitpc[aniter]->SetNsigmaTOF500_800(2.0);
	  dtc1etaphitpc[aniter]->SetNsigmaTOF800_1000(1.5);
	  dtc1etaphitpc[aniter]->SetNsigmaTOFge1000(1.0);
	  //----Kch special PID since 2013 run1 ---------<
	  //Track quality cuts: single Kaon 
	  dtc1etaphitpc[aniter]->SetminTPCncls(80);
	  dtc1etaphitpc[aniter]->SetRemoveKinks(kTRUE);
	  dtc1etaphitpc[aniter]->SetMaxTPCChiNdof(4.0);
	  //dtc1etaphitpc[aniter]->SetMaxITSChiNdof(36);	  
	  dtc1etaphitpc[aniter]->SetLabel(kFALSE);
	  dtc1etaphitpc[aniter]->SetMaxImpactZ(DCAxy/*0.3*/);//0.13 - kaons
	  dtc1etaphitpc[aniter]->SetMaxImpactXY(DCAz/*0.3*/);//0.135 - kaons

	  //PID for kaons or pions ---<


	  cutPass1YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutPass1%stpcM%i", chrgs[ichg], imult), 0.493677);
	  cutFail1YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutFail1%stpcM%i", chrgs[ichg], imult), 0.493677);
	  dtc1etaphitpc[aniter]->AddCutMonitor(cutPass1YPtetaphitpc[aniter], cutFail1YPtetaphitpc[aniter]);
	  
	  cutPass1PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cutPass1%stpcM%i", chrgs[ichg], imult),1);
	  cutFail1PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cutFail1%stpcM%i", chrgs[ichg], imult),1);
	  dtc1etaphitpc[aniter]->AddCutMonitor(cutPass1PIDetaphitpc[aniter], cutFail1PIDetaphitpc[aniter]);
	  /* //pi-pi 
	  sqpcetaphitpc[aniter] = new AliFemtoShareQualityPairCut();
	  sqpcetaphitpc[aniter]->SetShareQualityMax(1.0);
	  sqpcetaphitpc[aniter]->SetShareFractionMax(0.05);
	  sqpcetaphitpc[aniter]->SetRemoveSameLabel(kFALSE);
	  */
	  //K-K two track cuts: ------------------------------------------->
	  sqpcetaphitpc[aniter] = new AliFemtoPairCutRadialDistanceKKdist();
	  sqpcetaphitpc[aniter]->SetShareQualityMax(1.0);
	  sqpcetaphitpc[aniter]->SetShareFractionMax(0.05);
	  sqpcetaphitpc[aniter]->SetAverageSeparation(3.0);
	  sqpcetaphitpc[aniter]->SetRemoveSameLabel(kFALSE);
	  //K-K two track cuts: -------------------------------------------<
	  

	  anetaphitpc[aniter]->SetEventCut(mecetaphitpc[aniter]);
	  anetaphitpc[aniter]->SetFirstParticleCut(dtc1etaphitpc[aniter]);
	  anetaphitpc[aniter]->SetSecondParticleCut(dtc1etaphitpc[aniter]);// to be checked ????
	  anetaphitpc[aniter]->SetPairCut(sqpcetaphitpc[aniter]);
	  
	  //Correlation functions

	  //Spherical harmonics (without kT bins)
	  //cylmetaphitpc[aniter] = new AliFemtoCorrFctnDirectYlm(Form("cylm%stpcM%i", chrgs[ichg], imult),3,nbinssh,0.0,shqmax,runshlcms);
	  //anetaphitpc[aniter]->AddCorrFctn(cylmetaphitpc[aniter]);

	  //Qinv (without kT bins)-->
	  cqinvkttpc[aniter] = new AliFemtoQinvCorrFctn(Form("cqinv%stpcM%i", chrgs[ichg], imult),nbinssh,0.0,shqmax);
	  anetaphitpc[aniter]->AddCorrFctn(cqinvkttpc[aniter]);
	  //Qinv (without kT bins)--<
	  
	  //3D cartesian (without kT bins) --->
	  if(run3d){
	    //cq3dlcmskttpc[aniter] = new AliFemtoCorrFctn3DLCMSSym(Form("cq3d%stpcM%i", chrgs[ichg], imult),60,0.5);
	    cq3dlcmskttpc[aniter] = new AliFemtoBPLCMS3DCorrFctnKK(Form("cq3d%stpcM%ikTall", chrgs[ichg], imult),150,-3.0,3.0);
	    anetaphitpc[aniter]->AddCorrFctn(cq3dlcmskttpc[aniter]);
	    
	  }
	  //3D cartesian (without kT bins) <---

	  // cqinvnclstpc[aniter] = new AliFemtoCorrFctnTPCNcls(Form("cqinvncls%stpcM%i", chrgs[ichg], imult),nbinssh,0.0,shqmax);
	  // anetaphitpc[aniter]->AddCorrFctn(cqinvnclstpc[aniter]);

	  // cqinvchi2tpc[aniter] = new AliFemtoChi2CorrFctn(Form("cqinvchi2%stpcM%i", chrgs[ichg], imult),nbinssh,0.0,shqmax);
	  // anetaphitpc[aniter]->AddCorrFctn(cqinvchi2tpc[aniter]);

	  if (runktdep) {
	    int ktm;
	    for (int ikt=0; ikt<cktbins/*2*/; ikt++) {
	      ktm = aniter*cktbins + ikt;
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

	      //3d for kT bins ---->
	      if (run3d) {
		//cq3dlcmskttpc[ktm] = new AliFemtoCorrFctn3DLCMSSym(Form("cq3d%stpcM%ikT%i", chrgs[ichg], imult, ikt),60,(imult>3)?((imult>6)?((imult>7)?0.6:0.4):0.25):0.15);
		//cq3dlcmskttpc[ktm] = new AliFemtoCorrFctn3DLCMSSym(Form("cq3d%stpcM%ikT%i", chrgs[ichg], imult, ikt),60,0.5);
		cq3dlcmskttpc[ktm] = new AliFemtoBPLCMS3DCorrFctnKK(Form("cq3d%stpcM%ikT%i", chrgs[ichg], imult, ikt),150,-3.0,3.0);
		cq3dlcmskttpc[ktm]->SetPairSelectionCut(ktpcuts[ktm]);
		anetaphitpc[aniter]->AddCorrFctn(cq3dlcmskttpc[ktm]);
	      }
	      //3d for kT bins ----<
	    }
	  }
	  
	  /*cdedpetaphi[aniter] = new AliFemtoCorrFctnDEtaDPhi(Form("cdedp%stpcM%i", chrgs[ichg], imult),39, 39);
	    anetaphitpc[aniter]->AddCorrFctn(cdedpetaphi[aniter]);*/
	  
	  Manager->AddAnalysis(anetaphitpc[aniter]);	
	}
      }
    }
  }
  // *** End pipi analysis

  return Manager;
}                         
                      
