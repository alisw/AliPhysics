
/*********************************************************************
 *                                                                   *
 * ConfigFemtoAnalysis.C - configuration macro for the femtoscopic   *
 * analysis, meant as a QA process for two-particle effects          *
 *                                                                   *
 * Author: Adam Kisiel (Adam.Kisiel@cern.ch)                         *
 *********************************************************************
 *    K+K- in PbPb@2.76TeV                                           *
 * Update: Konstantin.Mikhaylov@cern.ch                              *
 *         KpmHBT16  0010  03-Mar-2016                               *
 *         TOF PID from 0.45 GeV/c                                   *
 *         SetNsigmaTPC400_450(1.0) and other set to standard        *
 *         no cut on mom and no on gamma                             *
 *         no cut on phi* eta (set to 0)                             *
 *         10 MeV/q-bin in CF                                        *
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
#include "AliFemtoPairCutAntiGammaAlpha.h"
#include "AliFemtoPairCutRadialDistance.h"
#include "AliFemtoPairCutRadialDistanceKK.h"
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
#include "AliFemtoModelCorrFctn.h"
#include "AliFemtoKTPairCut.h"
//... double ratio ...
#include "AliFemtoCorrFctnNonIdDR.h"
#endif

//________________________________________________________________________
AliFemtoManager* ConfigFemtoAnalysis3fm() {

  double PionMass = 0.13956995;
  double KaonMass = 0.493677;
  const int cMu=1;
  //const int cKt=2;//test
  const int cKt=8;//test

  //-------Single track cuts------------------------------------------------->
  //March 2017, CERN, do strict DCA cuts (one oder slaller)
  double DCAxy=0.24;//2.4;//cm // our standard is 0.20 cm; super narrow was 0.015cm
  double DCAz =0.3;//3.0;//cm // our standard is 0.15 cm;
  //wide DCA
  //double DCAxy=2.4;//cm // our standard is 0.20 cm; super narrow was 0.015cm
  //double DCAz =3.0;//cm // our standard is 0.15 cm;
  //-------Single track cuts-------------------------------------------------<
  //=======Double track cuts=================================================>
  //Dhevan's : PhiStarDifferenceMinimum=0.06; EtaDifferenceMinimum=0.02;
  double PhiStarDifferenceMinimum=0.;//0.02; //[radian]
  double EtaDifferenceMinimum=0.;//0.04; //[radian]
  //=======Double track cuts=================================================<

  // Switches for QA analyses
  int runmults[cMu] = {1};//, 0, 0, 0};
  int multbins[cMu+1] = {0, 900};//, 300, 500, 900};

  int runch[2] = {1, 0};//K+-
  const char *chrgs[2] = { "Kp", "Km"};
  
  int runktdep = 1;
  //double ktrng[cKt+1] = {0.2, 0.5, 1.0};//orig  
  //double ktrng[cKt+1] = {0.2, 0.5, 1.5};
  double ktrng[cKt+1] = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.3};
  
  

  int run3d = 0;
  int runshlcms = 0;

  int runtype = 2; // Types 0 - global, 1 - ITS only, 2 - TPC Inner
  int isrealdata = 1;

  //   AliFemtoEventReaderESDChainKine* Reader=new AliFemtoEventReaderESDChainKine();
  //   Reader->SetConstrained(true);
  //   Reader->SetUseTPCOnly(false);

  double shqmax;
  double shqmaxSH;
  int nbinssh = 200;// 10 MeV/q-bin //orig 100, 20 MeV per bin

  if (runshlcms) shqmaxSH = 0.25;
  shqmax = 2.0;//K+-
  
;

  AliFemtoEventReaderAODChain *Reader = new AliFemtoEventReaderAODChain();
    Reader->SetFilterBit(7);
    Reader->SetCentralityPreSelection(0.0001, 900);
    Reader->SetReadMC(kTRUE);
    Reader->SetDCAglobalTrack(kTRUE);//proverit'!!!
    //
    Reader->SetKaonAnalysis(kTRUE);// new set to choose Kaon by PDG code (7Nov2016)


    /*AliFemtoModelGausLCMSFreezeOutGenerator *tFreeze = new AliFemtoModelGausLCMSFreezeOutGenerator();
    tFreeze->SetSizeOut(3.0*TMath::Sqrt(2.0));         
    tFreeze->SetSizeSide(3.0*TMath::Sqrt(2.0));
    tFreeze->SetSizeLong(3.0*TMath::Sqrt(2.0));*/

    //Generate freeze-out coordinates as a 3D gaussian sphere in PRF
    AliFemtoModelGausRinvFreezeOutGenerator *tFreeze = new AliFemtoModelGausRinvFreezeOutGenerator();
    tFreeze->SetSizeInv(5.0*TMath::Sqrt(2.0));//r_0=3fm it should be time to sqrt(2) !!!!KM
    
    
    AliFemtoModelWeightGeneratorLednicky *tWeight = new AliFemtoModelWeightGeneratorLednicky();
    tWeight->SetPairType(AliFemtoModelWeightGenerator::KaonPlusKaonMinus());
    //tWeight->SetCoulOff();
    tWeight->SetCoulOn();
    tWeight->SetQuantumOff();
    //tWeight->SetStrongOff();
    tWeight->SetStrongOn();
    tWeight->Set3BodyOff();

    tWeight->SetKpKmModelType(14,1);//(Martin-f0, Achasov2-a0,0->phi off;1->phi on)
    //tWeight->SetKpKmModelType(22,1);//(Antonelli-f0, Antonelli-a0,0->phi off;1->phi on)
    tWeight->SetSphere();
    tWeight->SetSphere();

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
  AliFemtoKKTrackCut           *dtc1etaphitpc[20];
  AliFemtoKKTrackCut           *dtc2etaphitpc[20];
  //AliFemtoKpm45TrackCut           *dtc1etaphitpc[20];
  //AliFemtoKpm45TrackCut           *dtc2etaphitpc[20];
  //  AliFemtoESDTrackCut           *dtc1etaphitpc[20];
  //  AliFemtoESDTrackCut           *dtc2etaphitpc[20];
  AliFemtoCutMonitorParticleYPt *cutPass1YPtetaphitpc[20];
  AliFemtoCutMonitorParticleYPt *cutFail1YPtetaphitpc[20];
  AliFemtoCutMonitorParticlePID *cutPass1PIDetaphitpc[20];
  AliFemtoCutMonitorParticlePID *cutFail1PIDetaphitpc[20];
  AliFemtoCutMonitorParticleYPt *cutPass2YPtetaphitpc[20];
  AliFemtoCutMonitorParticleYPt *cutFail2YPtetaphitpc[20];
  AliFemtoCutMonitorParticlePID *cutPass2PIDetaphitpc[20];
  AliFemtoCutMonitorParticlePID *cutFail2PIDetaphitpc[20];
  AliFemtoPairCutAntiGammaAlpha      *sqpcetaphitpc[20];//K+-
  //    AliFemtoShareQualityTPCEntranceSepPairCut      *sqpcetaphitpc[20];
  //AliFemtoPairCutRadialDistance      *sqpcetaphitpc[20];//AliFemto dphi* cut
  //~~~
  ////////////////////AliFemtoPairCutRadialDistanceKK      *sqpcetaphitpc[20];//Dhevan's dphi* cut
  //
  //AliFemtoCorrFctnDirectYlm     *cylmetaphitpc[20];
  AliFemtoCorrFctnDEtaDPhi      *cdedpetaphi[20*10];//20->20*10 due to kT
  //AliFemtoChi2CorrFctn          *cchiqinvetaphitpc[20];
  AliFemtoKTPairCut             *ktpcuts[20*8];
  //AliFemtoCorrFctnDirectYlm     *cylmkttpc[20*8];
  AliFemtoQinvCorrFctn          *cqinvkttpc[20*8];
  //AliFemtoCorrFctn3DLCMSSym     *cq3dlcmskttpc[20*8];
  AliFemtoCorrFctnTPCNcls       *cqinvnclstpc[20];
  AliFemtoShareQualityCorrFctn  *cqinvsqtpc[20*10];
  AliFemtoChi2CorrFctn          *cqinvchi2tpc[20];
  AliFemtoTPCInnerCorrFctn      *cqinvinnertpc[20*10];
  AliFemtoCorrFctnGammaMonitorAlpha  *cgamma[20*10];
  //... double ratio ...
  // AliFemtoCorrFctnNonIdDR       *cfdourat[20*10];
  AliFemtoModelCorrFctn   *cqinvkttpcmodel[20*8];
 // AliFemtoCorrFctnGammaMonitorAlpha  *cgamma[20*10];


  // *** Third QA task - HBT analysis with all pair cuts off, TPC only ***
  // *** Begin K+K- analysis ***
  int aniter = 0;

  for (int imult=0; imult<cMu/*4*/; imult++) {
    if (runmults[imult]) {
      for (int ichg=0; ichg<1/*K+-*/; ichg++) {//one loop
	if (runch[ichg]) {
	  aniter = ichg*5+imult;

	  anetaphitpc[aniter] = new AliFemtoVertexMultAnalysis(4, -8.0, 8.0, 5, multbins[imult], multbins[imult+1]);
	  anetaphitpc[aniter]->SetNumEventsToMix(3);
	  anetaphitpc[aniter]->SetMinSizePartCollection(1);

	  mecetaphitpc[aniter] = new AliFemtoBasicEventCut();
	  mecetaphitpc[aniter]->SetEventMult(0,100000);
	  mecetaphitpc[aniter]->SetVertZPos(-8.0,8.0);
	  //    mecetaphitpc->SetAcceptBadVertex(kTRUE);
	  
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
	  
	  dtc1etaphitpc[aniter] = new AliFemtoKKTrackCut();
	  dtc2etaphitpc[aniter] = new AliFemtoKKTrackCut();//K+-
	  //dtc1etaphitpc[aniter] = new AliFemtoKpm45TrackCut();
	  //dtc2etaphitpc[aniter] = new AliFemtoKpm45TrackCut();//K+-
	  //--- K+- --->
	  dtc1etaphitpc[aniter]->SetCharge(1.0);
	  dtc2etaphitpc[aniter]->SetCharge(-1.0);
	  //--- K+- ---<
	    
	  dtc1etaphitpc[aniter]->SetPt(0.14,1.5); //0.14,1.5); //--- K+-
	  dtc1etaphitpc[aniter]->SetEta(-0.8,0.8);//--- K+-
	  dtc2etaphitpc[aniter]->SetPt(0.14,1.5); //0.14,1.5);  //--- K+-
	  dtc2etaphitpc[aniter]->SetEta(-0.8,0.8); //--- K+-
	  
	  //--- Choose Kaon as Most Probable (switch on all cuts: TPC, TOF)---
	  dtc1etaphitpc[aniter]->SetMass(KaonMass);
	  dtc1etaphitpc[aniter]->SetMostProbableKaon();
	  dtc2etaphitpc[aniter]->SetMass(KaonMass);
	  dtc2etaphitpc[aniter]->SetMostProbableKaon();
//------------------- November 2013 -----------------------------------< 
         // new cuts to remove electron (do not take into analysis if 400<p<500) 
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
	 //K+ contour (default p0=9999, p1=0, p2=0):  
	 //dtc1etaphitpc[aniter]->SetKaonTCPdEdxContourP0( 196.25);//december 2014
	 //dtc1etaphitpc[aniter]->SetKaonTCPdEdxContourP1(-420.00);//december 2014
	 //dtc1etaphitpc[aniter]->SetKaonTCPdEdxContourP2( 300.00);//december 2014

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
	 //K+ contour (default p0=9999, p1=0, p2=0):  
	 //dtc2etaphitpc[aniter]->SetKaonTCPdEdxContourP0( 196.25);//december 2014
	 //dtc2etaphitpc[aniter]->SetKaonTCPdEdxContourP1(-420.00);//december 2014
	 //dtc2etaphitpc[aniter]->SetKaonTCPdEdxContourP2( 300.00);//december 2014
	 //------------------- November 2013 ----------------------------------->
	  
	  // *** Track quality cuts ***
	 //K+
	 dtc1etaphitpc[aniter]->SetStatus(AliESDtrack::kTPCin);
	 dtc1etaphitpc[aniter]->SetminTPCncls(80); 
	 dtc1etaphitpc[aniter]->SetRemoveKinks(kTRUE);
	 dtc1etaphitpc[aniter]->SetLabel(kFALSE);
	 dtc1etaphitpc[aniter]->SetMaxTPCChiNdof(4.0);
	 dtc1etaphitpc[aniter]->SetMaxImpactXY(DCAxy);
	 dtc1etaphitpc[aniter]->SetMaxImpactZ(DCAz); 
	 //K- 
	 dtc2etaphitpc[aniter]->SetStatus(AliESDtrack::kTPCin);
	 dtc2etaphitpc[aniter]->SetminTPCncls(80); 
	 dtc2etaphitpc[aniter]->SetRemoveKinks(kTRUE);
	 dtc2etaphitpc[aniter]->SetLabel(kFALSE);
	 dtc2etaphitpc[aniter]->SetMaxTPCChiNdof(4.0);
	 dtc2etaphitpc[aniter]->SetMaxImpactXY(DCAxy);
	 dtc2etaphitpc[aniter]->SetMaxImpactZ(DCAz);  
	  
	  //K+
	  cutPass1YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutPass1%stpcM%i", chrgs[0], imult), 0.493677);
	  cutFail1YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutFail1%stpcM%i", chrgs[0], imult), 0.493677);
	  dtc1etaphitpc[aniter]->AddCutMonitor(cutPass1YPtetaphitpc[aniter], cutFail1YPtetaphitpc[aniter]);
	  //K- 
	  cutPass2YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutPass2%stpcM%i", chrgs[1], imult), 0.493677);
	  cutFail2YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutFail2%stpcM%i", chrgs[1], imult), 0.493677);
	  dtc2etaphitpc[aniter]->AddCutMonitor(cutPass2YPtetaphitpc[aniter], cutFail2YPtetaphitpc[aniter]);
	  
	  //K+
	  cutPass1PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cutPass1%stpcM%i", chrgs[0], imult),1);
	  cutFail1PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cutFail1%stpcM%i", chrgs[0], imult),1);
	  dtc1etaphitpc[aniter]->AddCutMonitor(cutPass1PIDetaphitpc[aniter], cutFail1PIDetaphitpc[aniter]);

	  //K-
	  cutPass2PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cutPass2%stpcM%i", chrgs[1], imult),1);
	  cutFail2PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cutFail2%stpcM%i", chrgs[1], imult),1);
	  dtc1etaphitpc[aniter]->AddCutMonitor(cutPass2PIDetaphitpc[aniter], cutFail2PIDetaphitpc[aniter]);

	  
	 // AliFemtoCorrFctnGammaMonitorAlpha
	   sqpcetaphitpc[aniter] = new AliFemtoPairCutAntiGammaAlpha();
	 
	  //  sqpcetaphitpc[aniter] = new AliFemtoShareQualityTPCEntranceSepPairCut();
	  
          //if (ichg < 2) {
	    //sqpcetaphitpc[aniter] = new AliFemtoPairCutRadialDistance();//AliFemto dphi* cut

//	    sqpcetaphitpc[aniter] = new AliFemtoPairCutRadialDistanceKK();  //Dhevan's dphi* cut

//~~~	  sqpcetaphitpc[aniter] = new AliFemtoPairCutAntiGamma();

	  sqpcetaphitpc[aniter]->SetRemoveSameLabel(kFALSE);
          //e+e-e+e-e+e-e+e-e+e-e+e-e+e-e+e- Set cut on gumma conversion e+e-
	  sqpcetaphitpc[aniter]->SetMaxEEMinv(0.05);
	  sqpcetaphitpc[aniter]->SetMaxAlphaDiff(0.9998);
	  sqpcetaphitpc[aniter]->SetTPCEntranceSepMinimum(0.00001);
//~~~	  sqpcetaphitpc[aniter]->SetMaxEEMinv(0.0);
//~~~	  sqpcetaphitpc[aniter]->SetMaxThetaDiff(0.0);
//~~~	  sqpcetaphitpc[aniter]->SetTPCEntranceSepMinimum(0.0);

     	   //--------- km:  eta-phi* Dhevan's custs ----------->>>>
     //	    sqpcetaphitpc[aniter]->SetPhiStarDifferenceMinimum(PhiStarDifferenceMinimum);
    //        sqpcetaphitpc[aniter]->SetEtaDifferenceMinimum(EtaDifferenceMinimum);
           //sqpcetaphitpc[aniter]->SetMinimumRadius(0.8);//not need for AliFemtoPairCutRadialDistanceKK()
           

	  
	  anetaphitpc[aniter]->SetEventCut(mecetaphitpc[aniter]);
	  anetaphitpc[aniter]->SetFirstParticleCut(dtc1etaphitpc[aniter]);//K+
	  anetaphitpc[aniter]->SetSecondParticleCut(dtc2etaphitpc[aniter]);//K-
	  anetaphitpc[aniter]->SetPairCut(sqpcetaphitpc[aniter]);

	  //Qinv (without kT bins)-->
	  // cqinvkttpc[aniter] = new AliFemtoQinvCorrFctn(Form("cqinv%stpcM%i", chrgs[ichg], imult),nbinssh,0.0,shqmax);
	  //  anetaphitpc[aniter]->AddCorrFctn(cqinvkttpc[aniter]);

	  //theoretical Qinv
	//    cqinvkttpcmodel[aniter] = new AliFemtoModelCorrFctn(Form("cqinvModel%stpcM%i", chrgs[ichg], imult),nbinssh,0.0,shqmax);
	    
	//    cqinvkttpcmodel[aniter]->ConnectToManager(tModelManager);
	//    anetaphitpc[aniter]->AddCorrFctn(cqinvkttpcmodel[aniter]);

      

	  if (runktdep) {
	    int ktm;
	    for (int ikt=0; ikt<cKt/*8*/; ikt++) {
	      ktm = aniter*cKt/*8*/ + ikt;
	      ktpcuts[ktm] = new AliFemtoKTPairCut(ktrng[ikt], ktrng[ikt+1]);
	      
              cqinvkttpc[ktm] = new AliFemtoQinvCorrFctn(Form("cqinv%stpcM%ikT%i", chrgs[ichg], imult, ikt),nbinssh,0.0,(imult>6)?shqmax*2.5:shqmax);

	      cqinvkttpc[ktm]->SetPairSelectionCut(ktpcuts[ktm]);//add kT bins
	      anetaphitpc[aniter]->AddCorrFctn(cqinvkttpc[ktm]);// add CF histos

	      //~~~
	      //cqinvsqtpc[ktm] = new AliFemtoShareQualityCorrFctn(Form("cqinvsq%stpcM%ikT%i", chrgs[ichg], imult, ikt),nbinssh,0.0,shqmax);
	      //cqinvsqtpc[ktm]->SetPairSelectionCut(ktpcuts[ktm]);
	      //  anetaphitpc[aniter]->AddCorrFctn(cqinvsqtpc[ktm]);

           // model--------------
           
          cqinvkttpcmodel[ktm] = new AliFemtoModelCorrFctn(Form("cqinvModel%stpcM%ikT%i", chrgs[ichg], imult, ikt),nbinssh,0.0,(imult>6)?shqmax*2.5:shqmax);

	      cqinvkttpcmodel[ktm]->SetPairSelectionCut(ktpcuts[ktm]);//add kT bins
	      cqinvkttpcmodel[ktm]->SetKaonPDG(kTRUE);//Special MC analysis for K selected by PDG code -->
              cqinvkttpcmodel[ktm]->ConnectToManager(tModelManager);
	      anetaphitpc[aniter]->AddCorrFctn(cqinvkttpcmodel[ktm]);// add CF histos
           
	    


              cgamma[aniter] = new AliFemtoCorrFctnGammaMonitorAlpha(Form("cgammaM%ikT%i", imult, ikt),200,200);
              anetaphitpc[aniter]->AddCorrFctn(cgamma[aniter]);



	      //... double ratio ...>
	      /*
	      cfdourat[ktm] = new 
		AliFemtoCorrFctnNonIdDR(Form("cfKstr%stpcM%ikT%i", chrgs[ichg], imult, ikt), 60, 0, 0.3);//AliFemtoCorrFctnNonIdDR(char* title, const int& nbins, const float& QinvLo, const float& QinvHi)
	      cfdourat[ktm]->SetPairSelectionCut(ktpcuts[ktm]);
	      anetaphitpc[aniter]->AddCorrFctn(cfdourat[ktm]);// add CF histos
	      //... double ratio ...<

	      */

	      //if (run3d) {
	      //cq3dlcmskttpc[ktm] = new AliFemtoCorrFctn3DLCMSSym(Form("cq3d%stpcM%ikT%i", chrgs[ichg], imult, ikt),60,(imult>3)?((imult>6)?((imult>7)?0.6:0.4):0.25):0.15);
	      //		//cq3dlcmskttpc[ktm] = new AliFemtoCorrFctn3DLCMSSym(Form("cq3d%stpcM%ikT%i", chrgs[ichg], imult, ikt),50,0.5);
	      //cq3dlcmskttpc[ktm]->SetPairSelectionCut(ktpcuts[ktm]);
	      //anetaphitpc[aniter]->AddCorrFctn(cq3dlcmskttpc[ktm]);
	      //}
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
                      
