//Analog PWGCF/FEMTOSCOPY/macros/Train/KpKmpp13TeV/ConfigFemtoAnalysis.C
// but with FB 128 and w/o Sphericity selection
// AOD data files, Run2
//
//Konstantin 9 June 2020 modification to K+K- 
//June 2020: No sphericity cut, filterbit 128.
// KchKch the same statistics is in PWGCF/FEMTOSCOPY/macros/Train/malinina/May2020/KK_13TeV_dist/ConfigFemtoAnalysis.C (Ludmila)
 
#if !defined(__CINT__) || defined(__MAKECINT_) 
#include "AliFemtoManager.h" 
#include "AliFemtoEventReaderESDChain.h" 
#include "AliFemtoEventReaderESDChainKine.h" 
#include "AliFemtoEventReaderAODChain.h" 
#include "AliFemtoSimpleAnalysis.h" 
#include "AliFemtoSphericityEventCut.h"
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
#include "AliFemtoCutMonitorEventSphericity.h"
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
#endif

//________________________________________________________________________
AliFemtoManager* ConfigFemtoAnalysis() {

  double PionMass = 0.13956995;
  double KaonMass = 0.493677;
 ///// const int cMu=4;
 //////// const int cKt=3;
  
  
  const int cMu=3;
  const int cKt=2;
  
  
  

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

//test Tom
//  int runmults[4] = {1, 0, 0, 0};
  
 // int multbins[5] = {0, 900, 300, 500, 900};
//old pp
  int multbins[5] = {1, 18, 30, 100, 1000};
 
 //test Tom
 // int multbins[5] = {1, 10000, 30, 100, 1000};
  
  
  
  //.................................................

  int runch[2] = {1, 1};
  const char *chrgs[2] = { "Kp", "Km"};
  
  
  int runktdep = 1;

  double ktrng[3] = {0.15, 0.5, 1.2};

// double ktrng[8] = {0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 2.0};

  int run3d = 0;
  int runshlcms = 0;

  int runtype = 2; // Types 0 - global, 1 - ITS only, 2 - TPC Inner
  int isrealdata = 1;

  //   AliFemtoEventReaderESDChainKine* Reader=new AliFemtoEventReaderESDChainKine();
  //   Reader->SetConstrained(true);
  //   Reader->SetUseTPCOnly(false);

  double shqmax;
  double shqmaxSH;
  int nbinssh = 100;
  int nbinsh3D = 160;

  if (runshlcms) shqmaxSH = 0.25;
  shqmax = 0.9;
  



AliFemtoEventReaderAODChain *Reader = new AliFemtoEventReaderAODChain();
//AliFemtoEventReaderAOD *Reader = new AliFemtoEventReaderAODMultSelection();
//  Reader->SetFilterMask(96);
  Reader->SetFilterBit(128);
//    Reader->SetFilterBit(7);
    
////    Reader->SetEPVZERO(kTRUE);
//    Reader->SetUseMultiplicity(AliFemtoEventReaderAOD::kCentrality);

//    Reader->SetNoCentrality(kTRUE);
//    Reader->SetUseMultiplicity(AliFemtoEventReaderAOD::kTPCOnlyRef);
//    Reader->SetUseMultiplicity(AliFemtoEventReaderAOD::kReference);
//    Reader->SetUseMultiplicity(AliFemtoEventReaderAOD::kGlobalCount);
     Reader->SetUseMultiplicity(AliFemtoEventReaderAOD::kRefComb08);
    
    Reader->SetCentralityFlattening(kFALSE);
    Reader->SetReadV0(0);
    // rdr->SetPrimaryVertexCorrectionTPCPoints(kTRUE);
    Reader->SetDCAglobalTrack(kTRUE);

//pile-up suppression
    Reader->SetUseMVPlpSelection(kTRUE);
    Reader->SetTrackPileUpRemoval(kTRUE);
    Reader->SetUseAliEventCuts(kTRUE);
    Reader->SetDCAglobalTrack(kTRUE);
    Reader->SetMinPlpContribSPD(3.0);
    Reader->SetIsPileUpEvent(kTRUE);


  AliFemtoManager* Manager=new AliFemtoManager();
  Manager->SetEventReader(Reader);

  AliFemtoVertexMultAnalysis    *anetaphitpc[20];
  AliFemtoBasicEventCut         *mecetaphitpc[20];
  //AliFemtoSphericityEventCut    *mecetaphitpc[20];
  AliFemtoCutMonitorEventMult   *cutPassEvMetaphitpc[20];
  AliFemtoCutMonitorEventMult   *cutFailEvMetaphitpc[20];
  AliFemtoCutMonitorEventVertex *cutPassEvVetaphitpc[20];
  AliFemtoCutMonitorEventVertex *cutFailEvVetaphitpc[20];
  AliFemtoCutMonitorCollections   *cutPassColletaphitpc[20];
  AliFemtoCutMonitorCollections   *cutFailColletaphitpc[20];
  AliFemtoKpm45TrackCut           *dtc1etaphitpc[20];
  AliFemtoKpm45TrackCut           *dtc2etaphitpc[20];
  AliFemtoCutMonitorParticleYPt *cutPass1YPtetaphitpc[20];
  AliFemtoCutMonitorParticleYPt *cutFail1YPtetaphitpc[20];
  AliFemtoCutMonitorParticlePID *cutPass1PIDetaphitpc[20];
  AliFemtoCutMonitorParticlePID *cutFail1PIDetaphitpc[20];
  AliFemtoCutMonitorParticleYPt *cutPass2YPtetaphitpc[20];
  AliFemtoCutMonitorParticleYPt *cutFail2YPtetaphitpc[20];
  AliFemtoCutMonitorParticlePID *cutPass2PIDetaphitpc[20];
  AliFemtoCutMonitorParticlePID *cutFail2PIDetaphitpc[20];
  AliFemtoCutMonitorEventSphericity  *cutPassEvSpher[20];
  AliFemtoCutMonitorEventSphericity  *cutFailEvSpher[20];

  AliFemtoPairCutRadialDistanceKKdist      *sqpcetaphitpc[20];//Dhevan's dphi* cut
 
  AliFemtoCorrFctnDirectYlm     *cylmetaphitpc[20];
  AliFemtoCorrFctnDEtaDPhi      *cdedpetaphi[20*10];//20->20*10 due to kT
  AliFemtoChi2CorrFctn          *cchiqinvetaphitpc[20];
  AliFemtoKTPairCut             *ktpcuts[20*8];
  AliFemtoCorrFctnDirectYlm     *cylmkttpc[20*8];
  AliFemtoQinvCorrFctn          *cqinvkttpc[20*8];
  AliFemtoBPLCMS3DCorrFctnKK  *cq3dlcmskttpc[20*8];
 
  AliFemtoCorrFctnTPCNcls       *cqinvnclstpc[20];
  AliFemtoShareQualityCorrFctn  *cqinvsqtpc[20*10];
  AliFemtoChi2CorrFctn          *cqinvchi2tpc[20];
  AliFemtoTPCInnerCorrFctn      *cqinvinnertpc[20*10];

  // *** Third QA task - HBT analysis with all pair cuts off, TPC only ***
  // *** Begin Kaon-Kaon (positive) analysis ***
  int aniter = 0;

  bool verbose=false;

  for (int imult=0; imult<cMu/*4*/; imult++) {
    if (runmults[imult]) {
      for (int ichg=0; ichg<1/*K+-*/; ichg++) {
	if (runch[ichg]) {
	  aniter = ichg*cMu+imult; //0, 1(ich=0) ,2,3


	  anetaphitpc[aniter] = new AliFemtoVertexMultAnalysis(10, -10.0, 10.0, 8, multbins[imult], multbins[imult+1]);
	  anetaphitpc[aniter]->SetNumEventsToMix(10);
	  anetaphitpc[aniter]->SetMinSizePartCollection(1);
	  anetaphitpc[aniter]->SetVerboseMode(verbose);
	  //FB 128:
	  mecetaphitpc[aniter] = new AliFemtoBasicEventCut();
	  mecetaphitpc[aniter]->SetEventMult(0.01,100000);
	  mecetaphitpc[aniter]->SetVertZPos(-10.0,10.0);

	  /*
	  mecetaphitpc[aniter] = new AliFemtoSphericityEventCut();
          mecetaphitpc[aniter]->SetEventMult(0.01,100000);
          mecetaphitpc[aniter]->SetVertZPos(-10.0,10.0);
	  //Sphericity selection:
          mecetaphitpc[aniter]->SetStMin(0.7);
          mecetaphitpc[aniter]->SetStMax(1.0);
	  */

	  //Multiplicity---->	  
	  cutPassEvMetaphitpc[aniter] = new AliFemtoCutMonitorEventMult(Form("cutPass%stpcM%i", chrgs[ichg], imult));
	  cutFailEvMetaphitpc[aniter] = new AliFemtoCutMonitorEventMult(Form("cutFail%stpcM%i", chrgs[ichg], imult));
	  mecetaphitpc[aniter]->AddCutMonitor(cutPassEvMetaphitpc[aniter], cutFailEvMetaphitpc[aniter]);
	  //Sphericity---->
 	  cutPassEvSpher[aniter] = new AliFemtoCutMonitorEventSphericity(Form("cutPass%stpcM%i", chrgs[ichg], imult));
	  cutFailEvSpher[aniter] = new AliFemtoCutMonitorEventSphericity(Form("cutFail%stpcM%i", chrgs[ichg], imult));
	  mecetaphitpc[aniter]->AddCutMonitor(cutPassEvSpher[aniter], cutFailEvSpher[aniter]);

	  	  
	  
	  //Main modification to K+K- starts here
	  dtc1etaphitpc[aniter] = new AliFemtoKpm45TrackCut();//K+
	  dtc2etaphitpc[aniter] = new AliFemtoKpm45TrackCut();//K-
	  //--- K+- --->
	  dtc1etaphitpc[aniter]->SetCharge(1.0);
	  dtc2etaphitpc[aniter]->SetCharge(-1.0);
	  //--- K+- ---<
	  //--- Choose Kaon as Most Probable (switch on all cuts: TPC, TOF)---
	  dtc1etaphitpc[aniter]->SetMass(KaonMass);
	  dtc1etaphitpc[aniter]->SetMostProbableKaon();
	  dtc2etaphitpc[aniter]->SetMass(KaonMass);
	  dtc2etaphitpc[aniter]->SetMostProbableKaon();
	  // new cuts to remove electron (do not take into analysis if 400<p<500) 
	  //K+
	  dtc1etaphitpc[aniter]->SetNsigmaTPCle250(2.0);
	  dtc1etaphitpc[aniter]->SetNsigmaTPC250_400(2.0);
	  dtc1etaphitpc[aniter]->SetNsigmaTPC400_450(1.0);//cut on e+e- orig(2.0);
	  dtc1etaphitpc[aniter]->SetNsigmaTPC450_500(2.0);//cut on e+e- orig(2.0);
	  dtc1etaphitpc[aniter]->SetNsigmaTPCge500(3.0);  
	  dtc1etaphitpc[aniter]->SetNsigmaTOF500_800(2.0);
	  dtc1etaphitpc[aniter]->SetNsigmaTOF800_1000(1.5);
	  dtc1etaphitpc[aniter]->SetNsigmaTOFge1000(1.0);
	  //K-
	  dtc2etaphitpc[aniter]->SetNsigmaTPCle250(2.0);
	  dtc2etaphitpc[aniter]->SetNsigmaTPC250_400(2.0);
	  dtc2etaphitpc[aniter]->SetNsigmaTPC400_450(1.0);//cut on e+e- orig(2.0);
	  dtc2etaphitpc[aniter]->SetNsigmaTPC450_500(2.0);//cut on e+e- orig(2.0);
	  dtc2etaphitpc[aniter]->SetNsigmaTPCge500(3.0);    
	  dtc2etaphitpc[aniter]->SetNsigmaTOF500_800(2.0);
	  dtc2etaphitpc[aniter]->SetNsigmaTOF800_1000(1.5);
	  dtc2etaphitpc[aniter]->SetNsigmaTOFge1000(1.0);

	  // *** Track quality cuts ***
	  //K+
	  dtc1etaphitpc[aniter]->SetMaxImpactXY(DCAxy);
	  dtc1etaphitpc[aniter]->SetMaxImpactZ(DCAz); 
	  //K- 
	  dtc2etaphitpc[aniter]->SetMaxImpactXY(DCAxy);
	  dtc2etaphitpc[aniter]->SetMaxImpactZ(DCAz);
	    
	  //New class to calculate track distance(possible not needed to K+K-)
	  sqpcetaphitpc[aniter] = new AliFemtoPairCutRadialDistanceKKdist();
	  sqpcetaphitpc[aniter]->SetDataType(AliFemtoPairCut::kAOD);
	  sqpcetaphitpc[aniter]->SetShareQualityMax(1.0);
	  sqpcetaphitpc[aniter]->SetShareFractionMax(0.05);
	  sqpcetaphitpc[aniter]->SetRemoveSameLabel(kFALSE);
	  sqpcetaphitpc[aniter]->SetAverageSeparation(12.0); //0.8
	  //AliFemtoVertexMultAnalysis additional settings
	  anetaphitpc[aniter]->SetEventCut(mecetaphitpc[aniter]);
	  anetaphitpc[aniter]->SetFirstParticleCut(dtc1etaphitpc[aniter]);//K+
	  anetaphitpc[aniter]->SetSecondParticleCut(dtc2etaphitpc[aniter]);//K-
	  anetaphitpc[aniter]->SetPairCut(sqpcetaphitpc[aniter]);
	  
	  //Qinv (without kT bins)-->
	  //cqinvkttpc[aniter] = new AliFemtoQinvCorrFctn(Form("cqinv%stpcM%i", chrgs[ichg], imult),nbinssh,0.0,shqmax);
	  //anetaphitpc[aniter]->AddCorrFctn(cqinvkttpc[aniter]);
	  //Qinv (without kT bins)--<

	  if (runktdep) {
	    int ktm;
	    for (int ikt=0; ikt<cKt/*8*/; ikt++) {
	      ktm = aniter*cKt/*8*/ + ikt;
	      ktpcuts[ktm] = new AliFemtoKTPairCut(ktrng[ikt], ktrng[ikt+1]);	      
	      cqinvkttpc[ktm] = new AliFemtoQinvCorrFctn(Form("cqinv%stpcM%ikT%i", chrgs[ichg], imult, ikt),nbinssh,0.0,2.0);
	      cqinvkttpc[ktm]->SetPairSelectionCut(ktpcuts[ktm]);
	      anetaphitpc[aniter]->AddCorrFctn(cqinvkttpc[ktm]);
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
                      
