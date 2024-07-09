//---------------------------------------------------------------->
// MC PURITY Pb-Pb 5.02 TeV                                          
// ConfigFemtoAnalysis.C 
// July 2024, Konstantin.Mikhaylov@cern.ch                           
//
//  few multiplicities and one kT to find out MC purity
//----------------------------------------------------------------<

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
#include "AliFemtoCutMonitorParticlePtPDG.h"
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
#include "AliFemtoModelCorrFctnTrueQ3D.h"
#include "AliFemtoModelCorrFctn3DKKGR.h"
#endif

//________________________________________________________________________
AliFemtoManager* ConfigFemtoAnalysis() {

  double PionMass = 0.13956995;
  double KaonMass = 0.493677;
  const int cMu=3;//orig //Paper PbPb 8 mults.
  const int cKt=1;  //4

  //-------Single track cuts------------------------------------------------->
  double DCAxy=0.3;//cm // our standard is 0.20 cm; super narrow was 0.015cm
  double DCAz =0.3;//cm // our standard is 0.15 cm;
  //-------Single track cuts-------------------------------------------------<
  //=======Double track cuts=================================================>
  //Dhevan's : PhiStarDifferenceMinimum=0.06; EtaDifferenceMinimum=0.02;
  // double PhiStarDifferenceMinimum=0.02; //[radian]
  //double EtaDifferenceMinimum=0.04; //[radian]

  //resolution study
  double PhiStarDifferenceMinimum=0.0; //[radian]
  double EtaDifferenceMinimum=0.0; //[radian]


  //=======Double track cuts=================================================<

  // Switches for QA analyses
  int runmults[cMu] = {1,1};//, 0, 0, 0};
  //pPb Paper: 0−20, 20-40, 40−90 [%]
  //int multbins[cMu+1] = {0, 50, 900};//test of 2 multiplicities
  int multbins[cMu+1] = {0, 200, 400, 900};

  int runch[2] = {1, 1};//K+-
  const char *chrgs[2] = { "Kp", "Km"};
  
  int runktdep = 0;//1;
  //double ktrng[cKt+1] = {0.2, 0.5, 1.0};//orig  
//  double ktrng[cKt+1] = {0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.5 (!)}; not 1.3
  //double ktrng[cKt+1] = {0.2, 0.4, 0.6, 0.8, 1.2};
  double ktrng[cKt+1] = {0.2, 10.0};
  
  

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
  shqmax = 0.4;//K+-
  
  

  AliFemtoEventReaderAODChain *Reader = new AliFemtoEventReaderAODChain();
    Reader->SetFilterBit(128);
    Reader->SetUseMultiplicity(AliFemtoEventReaderAOD::kCentrality);
    Reader->SetReadMC(kTRUE);
    Reader->SetDCAglobalTrack(kTRUE);//proverit'!!!s

  AliFemtoManager* Manager=new AliFemtoManager();
    Manager->SetEventReader(Reader);

  AliFemtoVertexMultAnalysis    *anetaphitpc[20];
  AliFemtoBasicEventCut         *mecetaphitpc[20];
  AliFemtoCutMonitorEventMult   *cutPassEvMetaphitpc[20];
  AliFemtoCutMonitorEventMult   *cutFailEvMetaphitpc[20];
  AliFemtoCutMonitorEventVertex *cutPassEvVetaphitpc[20];
  AliFemtoCutMonitorEventVertex *cutFailEvVetaphitpc[20];
  //AliFemtoKKTrackCut           *dtc1etaphitpc[20];
  //AliFemtoKKTrackCut           *dtc2etaphitpc[20];
  AliFemtoKpm45TrackCut           *dtc1etaphitpc[20];
  AliFemtoKpm45TrackCut           *dtc2etaphitpc[20];
  //  AliFemtoESDTrackCut           *dtc1etaphitpc[20];
  //  AliFemtoESDTrackCut           *dtc2etaphitpc[20];
  AliFemtoCutMonitorParticleYPt *cutPass1YPtetaphitpc[20];
  AliFemtoCutMonitorParticleYPt *cutFail1YPtetaphitpc[20];
  AliFemtoCutMonitorParticlePtPDG *cutPass1PIDetaphitpc[20];
  AliFemtoCutMonitorParticlePtPDG *cutFail1PIDetaphitpc[20];
  AliFemtoPairCutRadialDistanceKK      *sqpcetaphitpc[20];//Dhevan's dphi* cut


  // *** Third QA task - HBT analysis with all pair cuts off, TPC only ***
  // *** Begin KK analysis ***
  int aniter = 0;
  bool verbose=false;
  for (int imult=0; imult<cMu/*4*/; imult++) {
    if (runmults[imult]) {
      for (int ichg=0; ichg<2/*K+-*/; ichg++) {//one loop
	  	if (runch[ichg]){
			aniter = ichg*cMu+imult;
			anetaphitpc[aniter] = new AliFemtoVertexMultAnalysis(4, -8.0, 8.0, 5, multbins[imult], multbins[imult+1]);
			anetaphitpc[aniter]->SetNumEventsToMix(20);
			anetaphitpc[aniter]->SetMinSizePartCollection(1);
			anetaphitpc[aniter]->SetVerboseMode(verbose);
			
			mecetaphitpc[aniter] = new AliFemtoBasicEventCut();
			mecetaphitpc[aniter]->SetEventMult(0,100000);
			mecetaphitpc[aniter]->SetVertZPos(-8.0,8.0);
			
			cutPassEvMetaphitpc[aniter] = new AliFemtoCutMonitorEventMult(Form("cutPass%stpcM%i", chrgs[ichg], imult));
			cutFailEvMetaphitpc[aniter] = new AliFemtoCutMonitorEventMult(Form("cutFail%stpcM%i", chrgs[ichg], imult));
			mecetaphitpc[aniter]->AddCutMonitor(cutPassEvMetaphitpc[aniter], cutFailEvMetaphitpc[aniter]);
			
			cutPassEvVetaphitpc[aniter] = new AliFemtoCutMonitorEventVertex(Form("cutPass%stpcM%i", chrgs[ichg], imult));
			cutFailEvVetaphitpc[aniter] = new AliFemtoCutMonitorEventVertex(Form("cutFail%stpcM%i", chrgs[ichg], imult));
			mecetaphitpc[aniter]->AddCutMonitor(cutPassEvVetaphitpc[aniter], cutFailEvVetaphitpc[aniter]);
			
			//dtc1etaphitpc[aniter] = new AliFemtoKKTrackCut();
			dtc1etaphitpc[aniter] = new AliFemtoKpm45TrackCut();
			dtc2etaphitpc[aniter] = new AliFemtoKpm45TrackCut();
			//	  dtc1etaphitpc[aniter] = new AliFemtoESDTrackCut();
			//     dtc1etaphitpc[aniter]->SetPidProbPion(0.0,1.001);
			//     dtc1etaphitpc[aniter]->SetPidProbMuon(0.0,1.0);
			//     dtc1etaphitpc[aniter]->SetPidProbKaon(0.0,1.0);
			//     dtc1etaphitpc[aniter]->SetPidProbProton(0.0,1.0);
			if (ichg == 0)
				dtc1etaphitpc[aniter]->SetCharge(1.0);
			else if (ichg == 1)
				dtc1etaphitpc[aniter]->SetCharge(-1.0);
				
			dtc1etaphitpc[aniter]->SetPt(0.14,1.5);
			//	  dtc1etaphitpc[aniter]->SetEta(-1.2,1.2);
			dtc1etaphitpc[aniter]->SetEta(-0.8,0.8); //0.5
			// 	//    dtc1etaphitpc[aniter]->SetEta(-0.5,0.5);
			//dtc1etaphitpc[aniter]->SetMass(PionMass);
			dtc1etaphitpc[aniter]->SetMass(KaonMass);
			
			////	  dtc1etaphitpc[aniter]->SetminTPCncls(80);
			///////   ----!!!!!!	   
			dtc1etaphitpc[aniter]->SetMostProbableKaon();  //!!!!!!
			//dtc1etaphitpc[aniter]->SetMostProbablePion();
			//------------------- November 2013 -----------------------------------< 
			//New class in AliFemo: PWGCF/FEMTOSCOPY/AliFemtoUser/AliFemtoKKTrackCut.cxx
			
			dtc1etaphitpc[aniter]->SetNsigmaTPCle250(2.0);
			dtc1etaphitpc[aniter]->SetNsigmaTPC250_400(2.0);
			dtc1etaphitpc[aniter]->SetNsigmaTPC400_450(1.0);
			dtc1etaphitpc[aniter]->SetNsigmaTPC450_500(2.0);
			dtc1etaphitpc[aniter]->SetNsigmaTPCge500(3.0);
			// new cuts are stronger, better separation of pion in TOF 
			// when momentum is greater then 800 MeV/c
			dtc1etaphitpc[aniter]->SetNsigmaTOF500_800(2.0);
			dtc1etaphitpc[aniter]->SetNsigmaTOF800_1000(1.5);
			dtc1etaphitpc[aniter]->SetNsigmaTOFge1000(1.0);
			
			//------------------- November 2013 ----------------------------------->
			
			////	  	  dtc1etaphitpc[aniter]->SetMostProbablePion();
			// 	// Track quality cuts
			if (runtype == 0) {
				dtc1etaphitpc[aniter]->SetStatus(AliESDtrack::kTPCrefit|AliESDtrack::kITSrefit);
				//	    dtc1etaphitpc[aniter]->SetStatus(AliESDtrack::kTPCrefit);
				//    dtc1etaphitpc[aniter]->SetStatus(AliESDtrack::kITSrefit);
				dtc1etaphitpc[aniter]->SetminTPCncls(80);
				dtc1etaphitpc[aniter]->SetRemoveKinks(kTRUE);
				dtc1etaphitpc[aniter]->SetLabel(kFALSE);
				//    dtc1etaphitpc[aniter]->SetMaxITSChiNdof(6.0);
				dtc1etaphitpc[aniter]->SetMaxTPCChiNdof(4.0);
				dtc1etaphitpc[aniter]->SetMaxImpactXY(DCAxy);
				//Poland: dtc1etaphitpc[aniter]->SetMaxImpactXYPtDep(0.0182, 0.0350, -1.01);
				dtc1etaphitpc[aniter]->SetMaxImpactZ(DCAz);
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
				dtc1etaphitpc[aniter]->SetMaxImpactXY(DCAxy);
				dtc1etaphitpc[aniter]->SetMaxImpactZ(DCAz);
				//      dtc1etaphitpc[aniter]->SetMaxSigmaToVertex(6.0);
			}
			else if (runtype == 2) {
				//	    dtc1etaphitpc[aniter]->SetStatus(AliESDtrack::kTPCrefit|AliESDtrack::kITSrefit);
				dtc1etaphitpc[aniter]->SetStatus(AliESDtrack::kTPCin);
				//	    dtc1etaphitpc[aniter]->SetStatus(AliESDtrack::kTPCrefit);
				//    dtc1etaphitpc[aniter]->SetStatus(AliESDtrack::kITSrefit);
				dtc1etaphitpc[aniter]->SetminTPCncls(80); //was "0"
				dtc1etaphitpc[aniter]->SetRemoveKinks(kTRUE);
				dtc1etaphitpc[aniter]->SetLabel(kFALSE);
				//    dtc1etaphitpc[aniter]->SetMaxITSChiNdof(6.0);
				dtc1etaphitpc[aniter]->SetMaxTPCChiNdof(4.0);
				dtc1etaphitpc[aniter]->SetMaxImpactXY(DCAxy);
				//dtc1etaphitpc[aniter]->SetMaxImpactXYPtDep(0.0182, 0.0350, -1.01);
				dtc1etaphitpc[aniter]->SetMaxImpactZ(DCAz);  //3.0
				//      dtc1etaphitpc[aniter]->SetMaxSigmaToVertex(6.0);
			}
			
			cutPass1YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutPass1%stpcM%i", chrgs[ichg], imult), 0.493677);
			cutFail1YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutFail1%stpcM%i", chrgs[ichg], imult), 0.493677);
			dtc1etaphitpc[aniter]->AddCutMonitor(cutPass1YPtetaphitpc[aniter], cutFail1YPtetaphitpc[aniter]);
			//purity: AliFemtoCutMonitorParticlePtPD -> histograms hard coded pT range: 100,0.1,2.0 
			cutPass1PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePtPDG(Form("cutPass1%stpcM%i", chrgs[ichg], imult),1);
			cutFail1PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePtPDG(Form("cutFail1%stpcM%i", chrgs[ichg], imult),1);
			dtc1etaphitpc[aniter]->AddCutMonitor(cutPass1PIDetaphitpc[aniter], cutFail1PIDetaphitpc[aniter]);

			if (ichg < 2)
				sqpcetaphitpc[aniter] = new AliFemtoPairCutRadialDistanceKK();  //Dhevan's dphi* cut
			
			
			anetaphitpc[aniter]->SetEventCut(mecetaphitpc[aniter]);
			anetaphitpc[aniter]->SetFirstParticleCut(dtc1etaphitpc[aniter]);//K+
			anetaphitpc[aniter]->SetSecondParticleCut(dtc1etaphitpc[aniter]);//K-
			anetaphitpc[aniter]->SetPairCut(sqpcetaphitpc[aniter]);
			
			Manager->AddAnalysis(anetaphitpc[aniter]);
		}
	  }
    }
  }
  // *** End Kaon-Kaon (positive) analysis

  return Manager;
}                         
                      
