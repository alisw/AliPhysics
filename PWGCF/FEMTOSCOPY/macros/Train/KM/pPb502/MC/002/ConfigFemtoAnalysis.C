//
//  December 2023: Konstantin.Mikhaylov@cern.ch
// MC  PWGCF/FEMTOSCOPY/macros/Train/KM/pPb502/MC/KM_pPb502_MC_002/ConfigFemtoAnalysis.C
// --> MonteCarlo wagon for DPMJET model:
// --> Add martix for Momentum Resolution as well as in PbPb /cern/kmikhail/Grid/GIT/AliPhysics/PWGCF/FEMTOSCOPY/macros/Train/KpmHBT_MC/Feb2018/ConfigFemtoAnalysis.C
//
#if !defined(__CINT__) || defined(__MAKECINT_)
#include "AliFemtoManager.h"
#include "AliFemtoEventReaderAOD.h"
#include "AliFemtoEventReaderAODMultSelection.h"
#include "AliFemtoSimpleAnalysis.h"
#include "AliFemtoBasicEventCut.h"
#include "AliFemtoSphericityEventCut.h"
#include "AliFemtoESDTrackCut.h"
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
#include "AliFemtoPairCutRadialDistanceKKdist.h"
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


AliFemtoManager* ConfigFemtoAnalysis() {
  double PionMass = 0.13956995;
  double KaonMass = 0.493677;
  Double_t DCAxy=0.135;//cm [0.135 - main]
  Double_t DCAz =0.130;//cm [0.130 - main]

  //multiplicity bins
  int runmults[3] = {1, 1, 1};
  int multbins[4] = {0, 200, 400, 900};

  int runch[2] = {1, 1};
  const char *chrgs[2] = { "Kp", "Km"};

  int runktdep = 1;
  double ktrng[3] = {0.2, 0.5, 1.0};
  //double ktrng[9] = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.3};
  int run3d = 0;//no 3d
  int runshlcms = 0;//1;
  double shqmax;
  int nbinssh = 200;
  if (runshlcms) shqmax = 0.25;
  else shqmax = 2.0;//CF at to 2 GeV/c

  AliFemtoEventReaderAODMultSelection *Reader = new AliFemtoEventReaderAODMultSelection();
  Reader->SetFilterBit(8);
  Reader->SetReadV0(1);
  Reader->SetEPVZERO(kTRUE);
  Reader->SetCentralityFlattening(kTRUE);
  Reader->SetUseMultiplicity(AliFemtoEventReaderAODMultSelection::kCentrality);
  Reader->SetDCAglobalTrack(kTRUE);

  //Model function K+K- as well as in PbPb^ PWGCF/FEMTOSCOPY/macros/Train/KpmHBT_MC/Feb2018/ConfigFemtoAnalysis.C ------>
   /*AliFemtoModelGausLCMSFreezeOutGenerator *tFreeze = new AliFemtoModelGausLCMSFreezeOutGenerator();
    tFreeze->SetSizeOut(3.0*TMath::Sqrt(2.0));         
    tFreeze->SetSizeSide(3.0*TMath::Sqrt(2.0));
    tFreeze->SetSizeLong(3.0*TMath::Sqrt(2.0));*/

    //Generate freeze-out coordinates as a 3D gaussian sphere in PRF
    AliFemtoModelGausRinvFreezeOutGenerator *tFreeze = new AliFemtoModelGausRinvFreezeOutGenerator();
    //tFreeze->SetSizeInv(3.0*TMath::Sqrt(2.0));//r_0=3fm it should be time to sqrt(2) !!!!KM
    //tFreeze->SetSizeInv(4.0*TMath::Sqrt(2.0));//r_0=4fm it should be time to sqrt(2) !!!!KM
    tFreeze->SetSizeInv(5.0*TMath::Sqrt(2.0));//r_0=5fm it should be time to sqrt(2) !!!!KM
    //tFreeze->SetSizeInv(6.0*TMath::Sqrt(2.0));//r_0=6fm it should be time to sqrt(2) !!!!KM

    //Feb 2018: Switch off all FSI and phi-meson
    AliFemtoModelWeightGeneratorLednicky *tWeight = new AliFemtoModelWeightGeneratorLednicky();
    tWeight->SetPairType(AliFemtoModelWeightGenerator::KaonPlusKaonMinus());
    tWeight->SetCoulOff();
    //tWeight->SetCoulOn();
    tWeight->SetQuantumOff();
    tWeight->SetStrongOff();
    //tWeight->SetStrongOn();
    tWeight->Set3BodyOff();

    //tWeight->SetKpKmModelType(14,1);//(Martin-f0, Achasov2-a0,0->phi off;1->phi on)
    //No phi:
    tWeight->SetKpKmModelType(14,0);//(Martin-f0, Achasov2-a0,0->phi off;1->phi on)
 

  AliFemtoModelManager *tModelManager = new AliFemtoModelManager();
  tModelManager->AcceptFreezeOutGenerator(tFreeze);
  tModelManager->AcceptWeightGenerator(tWeight);
  tModelManager->CreateCopyHiddenInfo(kTRUE);
  
  //Model function K+K- as well as in PbPb^ PWGCF/FEMTOSCOPY/macros/Train/KpmHBT_MC/Feb2018/ConfigFemtoAnalysis.C ------<
  
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
  AliFemtoKpm45TrackCut *dtc1etaphitpc[20];
  AliFemtoKpm45TrackCut *dtc2etaphitpc[20];
  AliFemtoCutMonitorParticleYPt *cutPass1YPtetaphitpc[20];
  AliFemtoCutMonitorParticleYPt *cutFail1YPtetaphitpc[20];
  AliFemtoCutMonitorParticlePID *cutPass1PIDetaphitpc[20];
  AliFemtoCutMonitorParticlePID *cutFail1PIDetaphitpc[20];
  AliFemtoCutMonitorParticleYPt *cutPass2YPtetaphitpc[20];
  AliFemtoCutMonitorParticleYPt *cutFail2YPtetaphitpc[20];
  AliFemtoCutMonitorParticlePID *cutPass2PIDetaphitpc[20];
  AliFemtoCutMonitorParticlePID *cutFail2PIDetaphitpc[20];
  AliFemtoCutMonitorEventSphericity *cutPassEvSpher[20];
  AliFemtoCutMonitorEventSphericity *cutFailEvSpher[20];

  AliFemtoPairCutRadialDistanceKKdist  *sqpcetaphitpc[20];
  AliFemtoCorrFctnDirectYlm     *cylmetaphitpc[20];
  AliFemtoCorrFctnDEtaDPhi      *cdedpetaphi[20];
  AliFemtoChi2CorrFctn          *cchiqinvetaphitpc[20];
  AliFemtoKTPairCut             *ktpcuts[20*8];
  AliFemtoCorrFctnDirectYlm     *cylmkttpc[20*8];
  AliFemtoQinvCorrFctn          *cqinvkttpc[20*8];
  AliFemtoBPLCMS3DCorrFctnKK     *cq3dlcmskttpc[20*8];
  AliFemtoCorrFctnTPCNcls       *cqinvnclstpc[20];
  AliFemtoShareQualityCorrFctn  *cqinvsqtpc[20*10];
  AliFemtoChi2CorrFctn          *cqinvchi2tpc[20];
  AliFemtoTPCInnerCorrFctn      *cqinvinnertpc[20*10];

  //...Begin K+K- analysis
  int aniter = 0;

  for (int imult=0; imult<3; imult++) {
    if (runmults[imult]) {
      for (int ichg=0; ichg<1/*2*/; ichg++) { //1 -> K+K- only 
	if (runch[ichg]) {
	  aniter = ichg*3+imult;
	  anetaphitpc[aniter] = new AliFemtoVertexMultAnalysis(10, -10.0, 10.0, 4, multbins[imult], multbins[imult+1]);
	  anetaphitpc[aniter]->SetNumEventsToMix(3);
	  anetaphitpc[aniter]->SetMinSizePartCollection(1);
	  anetaphitpc[aniter]->SetVerboseMode(kFALSE);

	  mecetaphitpc[aniter] = new AliFemtoBasicEventCut();
	  mecetaphitpc[aniter]->SetEventMult(0,10000);
	  mecetaphitpc[aniter]->SetVertZPos(-10,10);
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
	  
          cutPassEvSpher[aniter] = new AliFemtoCutMonitorEventSphericity(Form("cutPass%stpcM%i", chrgs[ichg], imult));
          cutFailEvSpher[aniter] = new AliFemtoCutMonitorEventSphericity(Form("cutFail%stpcM%i", chrgs[ichg], imult));
          mecetaphitpc[aniter]->AddCutMonitor(cutPassEvSpher[aniter], cutFailEvSpher[aniter]);

	  dtc1etaphitpc[aniter] = new AliFemtoKpm45TrackCut();//K+
	  dtc2etaphitpc[aniter] = new AliFemtoKpm45TrackCut();//K-
	  //--- K+- --->
	  dtc1etaphitpc[aniter]->SetCharge(1.0);
	  dtc2etaphitpc[aniter]->SetCharge(-1.0);
	  //--- K+- ---<
	   dtc1etaphitpc[aniter]->SetPt(0.14,1.5);//K+
	   dtc1etaphitpc[aniter]->SetEta(-0.8,0.8);//K+
	   dtc2etaphitpc[aniter]->SetPt(0.14,1.5);//K-
	   dtc2etaphitpc[aniter]->SetEta(-0.8,0.8);//K-
	   //PID method
	   //dtc1etaphitpc[aniter]->SetMass(PionMass);dtc1etaphitpc[aniter]->SetMostProbablePion();
	   dtc1etaphitpc[aniter]->SetMass(KaonMass);//K+
	   dtc1etaphitpc[aniter]->SetMostProbableKaon();//K+
	   dtc2etaphitpc[aniter]->SetMass(KaonMass);//K-
	   dtc2etaphitpc[aniter]->SetMostProbableKaon();//K-
	   
	  // new cuts to remove electron (do not take into analysis if 400<p<500)
	   //K+++++++++:
	   dtc1etaphitpc[aniter]->SetNsigmaTPCle250(2.0);
	   dtc1etaphitpc[aniter]->SetNsigmaTPC250_400(2.0);
	   dtc1etaphitpc[aniter]->SetNsigmaTPC400_450(1.0);//cut on e+e- orig(2.0);
	   dtc1etaphitpc[aniter]->SetNsigmaTPC450_500(2.0);//cut on e+e- orig(2.0);
	   dtc1etaphitpc[aniter]->SetNsigmaTPCge500(3.0);  
	   // new cuts are stronger, better separation of pion in TOF 
	   // when momentum is greater then 800 MeV/c
	   dtc1etaphitpc[aniter]->SetNsigmaTOF500_800(2.0);
	   dtc1etaphitpc[aniter]->SetNsigmaTOF800_1000(1.5);
	   dtc1etaphitpc[aniter]->SetNsigmaTOFge1000(1.0);
	   //K----------:
	 dtc2etaphitpc[aniter]->SetNsigmaTPCle250(2.0);
	 dtc2etaphitpc[aniter]->SetNsigmaTPC250_400(2.0);
	 dtc2etaphitpc[aniter]->SetNsigmaTPC400_450(1.0);//cut on e+e- orig(2.0);
	 dtc2etaphitpc[aniter]->SetNsigmaTPC450_500(2.0);//cut on e+e- orig(2.0);
	 dtc2etaphitpc[aniter]->SetNsigmaTPCge500(3.0);    
	 // new cuts are stronger, better separation of pion in TOF 
	 // when momentum is greater then 800 MeV/c
	 dtc2etaphitpc[aniter]->SetNsigmaTOF500_800(2.0);
	 dtc2etaphitpc[aniter]->SetNsigmaTOF800_1000(1.5);
	 dtc2etaphitpc[aniter]->SetNsigmaTOFge1000(1.0);
	 //Track quality cuts
	 //K+:
	 dtc1etaphitpc[aniter]->SetminTPCncls(80);
	 dtc1etaphitpc[aniter]->SetRemoveKinks(kTRUE);
	 dtc1etaphitpc[aniter]->SetMaxTPCChiNdof(4.0);
	 dtc1etaphitpc[aniter]->SetLabel(kFALSE);
	 dtc1etaphitpc[aniter]->SetMaxImpactZ(DCAz);//main
	 dtc1etaphitpc[aniter]->SetMaxImpactXY(DCAxy);//main
	 //K-:
	 dtc2etaphitpc[aniter]->SetminTPCncls(80);
	 dtc2etaphitpc[aniter]->SetRemoveKinks(kTRUE);
	 dtc2etaphitpc[aniter]->SetMaxTPCChiNdof(4.0);
	 dtc2etaphitpc[aniter]->SetLabel(kFALSE);
	 dtc2etaphitpc[aniter]->SetMaxImpactZ(DCAz);//main
	 dtc2etaphitpc[aniter]->SetMaxImpactXY(DCAxy);//main

	 cutPass1YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutPass1%stpcM%i", chrgs[ichg], imult), 0.493677);
	 cutFail1YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutFail1%stpcM%i", chrgs[ichg], imult), 0.493677);
	 dtc1etaphitpc[aniter]->AddCutMonitor(cutPass1YPtetaphitpc[aniter], cutFail1YPtetaphitpc[aniter]);
	 
	 cutPass1PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cutPass1%stpcM%i", chrgs[ichg], imult),1);
	 cutFail1PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cutFail1%stpcM%i", chrgs[ichg], imult),1);
	 dtc1etaphitpc[aniter]->AddCutMonitor(cutPass1PIDetaphitpc[aniter], cutFail1PIDetaphitpc[aniter]);

	 sqpcetaphitpc[aniter] = new AliFemtoPairCutRadialDistanceKKdist();
	 sqpcetaphitpc[aniter]->SetShareQualityMax(1.0);
	 sqpcetaphitpc[aniter]->SetShareFractionMax(0.05);
	 sqpcetaphitpc[aniter]->SetAverageSeparation(3.0);
	 sqpcetaphitpc[aniter]->SetRemoveSameLabel(kFALSE);
	 
	 anetaphitpc[aniter]->SetEventCut(mecetaphitpc[aniter]);
	 anetaphitpc[aniter]->SetFirstParticleCut(dtc1etaphitpc[aniter]);//K+
	 anetaphitpc[aniter]->SetSecondParticleCut(dtc2etaphitpc[aniter]);//K-;
	 anetaphitpc[aniter]->SetPairCut(sqpcetaphitpc[aniter]);

	 //...Correlation functions
	 //...Qinv (without kT bins)--->
	 cqinvkttpc[aniter] = new AliFemtoQinvCorrFctn(Form("cqinv%stpcM%i", chrgs[ichg], imult),nbinssh,0.0,shqmax);
	 anetaphitpc[aniter]->AddCorrFctn(cqinvkttpc[aniter]);
	 //...Qinv (without kT bins)---<
	 if (runktdep) {
	   int ktm;
	   for (int ikt=0; ikt<2; ikt++) {
	     ktm = aniter*2 + ikt;
	     ktpcuts[ktm] = new AliFemtoKTPairCut(ktrng[ikt], ktrng[ikt+1]);
	     //...C(qinv,kT)--->
	     cqinvkttpc[ktm] = new AliFemtoQinvCorrFctn(Form("cqinv%stpcM%ikT%i", chrgs[ichg], imult, ikt),nbinssh,0.0, shqmax);
	     cqinvkttpc[ktm]->SetPairSelectionCut(ktpcuts[ktm]);
	     anetaphitpc[aniter]->AddCorrFctn(cqinvkttpc[ktm]);
	     //...C(qinv,kT)---<
	     
	     // model CF for MR MTX QgenQrec-------------->
	     cqinvkttpcmodel[ktm] = new AliFemtoModelCorrFctnKK(Form("cqinvModel%stpcM%ikT%i", chrgs[ichg], imult, ikt),nbinssh,0.0,shqmax);//0.0,(imult>6)?shqmax*2.5:shqmax);

	      cqinvkttpcmodel[ktm]->SetPairSelectionCut(ktpcuts[ktm]);//add kT bins
	      cqinvkttpcmodel[ktm]->SetKaonPDG(kTRUE);//Special MC analysis for K selected by PDG code -->
              cqinvkttpcmodel[ktm]->ConnectToManager(tModelManager);
	      anetaphitpc[aniter]->AddCorrFctn(cqinvkttpcmodel[ktm]);// add CF histos
	     // model CF for MR MTX QgenQrec--------------<
 
	     cqinvsqtpc[ktm] = new AliFemtoShareQualityCorrFctn(Form("cqinvsq%stpcM%ikT%i", chrgs[ichg], imult, ikt),nbinssh,0.0,shqmax);
	     cqinvsqtpc[ktm]->SetPairSelectionCut(ktpcuts[ktm]);
	     anetaphitpc[aniter]->AddCorrFctn(cqinvsqtpc[ktm]);

	     cqinvinnertpc[ktm] = new AliFemtoTPCInnerCorrFctn(Form("cqinvinner%stpcM%ikT%i", chrgs[ichg], imult, ikt),nbinssh,0.0,shqmax);
	     cqinvinnertpc[ktm]->SetPairSelectionCut(ktpcuts[ktm]);
	     cqinvinnertpc[ktm]->SetRadius(1.2);
	     anetaphitpc[aniter]->AddCorrFctn(cqinvinnertpc[ktm]);
	   }
	 }
	 Manager->AddAnalysis(anetaphitpc[aniter]);
	}
      }
    }
  }

  return Manager;
}
