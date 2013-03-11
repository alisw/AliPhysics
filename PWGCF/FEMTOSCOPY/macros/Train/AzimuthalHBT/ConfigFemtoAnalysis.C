
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
#include "AliFemtoPairCutRadialDistanceLM.h"
#include "AliFemtoQinvCorrFctn.h"
#include "AliFemtoShareQualityCorrFctn.h"
#include "AliFemtoTPCInnerCorrFctn.h"
#include "AliFemtoAnalysisAzimuthalPbPb.h"
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
#endif

//________________________________________________________________________
AliFemtoManager* ConfigFemtoAnalysis() {
	
	double PionMass = 0.13956995;
	double KaonMass = 0.493677;
	
	int runmults[10] = {0, 0, 0, 1, 1, 1, 1, 0, 0, 0};
	int multbins[11] = {0, 50, 100, 200, 300, 400, 500, 700, 800, 900, 990};
	
	int runch[2] = {1, 1};
	const char *chrgs[2] = { "pip", "pim" };
	
	int runktdep = 1;
	double ktrng[5] = {0.2, 0.3, 0.4, 0.5, 0.7};
	int phirange[7] = {-15, 15, 45, 75, 105, 135, 165};
	
	int runtype = 2; // Types 0 - global, 1 - ITS only, 2 - TPC Inner
	int isrealdata = 0;
	
	AliFemtoEventReaderAODChain *Reader = new AliFemtoEventReaderAODChain();
	Reader->SetFilterBit(7);
	Reader->SetCentralityPreSelection(180, 701);
	
	AliFemtoManager* Manager=new AliFemtoManager();
	Manager->SetEventReader(Reader);
	
	AliFemtoAnalysisAzimuthalPbPb    *ana[2][10];
	AliFemtoBasicEventCut         *mecetaphitpc[2][10];
	AliFemtoCutMonitorEventMult   *cutPassEvMetaphitpc[2][10];
	AliFemtoCutMonitorEventMult   *cutFailEvMetaphitpc[2][10];
	AliFemtoCutMonitorEventVertex *cutPassEvVetaphitpc[2][10];
	AliFemtoCutMonitorEventVertex *cutFailEvVetaphitpc[2][10];
	AliFemtoESDTrackCut           *dtc1etaphitpc[2][10];
	AliFemtoESDTrackCut           *dtc2etaphitpc[2][10];
	//   AliFemtoCutMonitorParticleYPt *cutPass1YPtetaphitpc[20];
	//   AliFemtoCutMonitorParticleYPt *cutFail1YPtetaphitpc[20];
	//   AliFemtoCutMonitorParticlePID *cutPass1PIDetaphitpc[20];
	//   AliFemtoCutMonitorParticlePID *cutFail1PIDetaphitpc[20];
	//   AliFemtoCutMonitorParticleYPt *cutPass2YPtetaphitpc[20];
	//   AliFemtoCutMonitorParticleYPt *cutFail2YPtetaphitpc[20];
	//   AliFemtoCutMonitorParticlePID *cutPass2PIDetaphitpc[20];
	//   AliFemtoCutMonitorParticlePID *cutFail2PIDetaphitpc[20];
	//  AliFemtoPairCutAntiGamma      *sqpcetaphitpc[20];
	//  AliFemtoShareQualityTPCEntranceSepPairCut      *sqpcetaphitpc[20];
	AliFemtoPairCutRadialDistanceLM      *sqpcetaphitpc[2][10];
	AliFemtoPairCutRadialDistanceLM      *sqpcetaphitpcRD[2][10];
	//   AliFemtoCorrFctnDirectYlm     *cylmetaphitpc[20];
	//   AliFemtoCorrFctnDEtaDPhi      *cdedpetaphi[20];
	//   AliFemtoChi2CorrFctn          *cchiqinvetaphitpc[20];
	AliFemtoKTPairCut             *ktpaircut[2][10][6][8];
	//   AliFemtoCorrFctnDirectYlm     *cylmkttpc[20*7];
	//   AliFemtoQinvCorrFctn          *cqinvkttpc[20*7];
	AliFemtoBPLCMS3DCorrFctn     *cq3dlcmskttpc[2][10][6][8];
	//   AliFemtoCorrFctnTPCNcls       *cqinvnclstpc[20];
	//   AliFemtoShareQualityCorrFctn  *cqinvsqtpc[20*10];
	//   AliFemtoChi2CorrFctn          *cqinvchi2tpc[20];
	//   AliFemtoTPCInnerCorrFctn      *cqinvinnertpc[20*10];
	
	
	for (int imult=0; imult<10; imult++) {
		if (runmults[imult]) {
			for (int ichg=0; ichg<2; ichg++) {
				if (runch[ichg]) {
					
					ana[ichg][imult] = new AliFemtoAnalysisAzimuthalPbPb(4, -8.0, 8.0, 5, multbins[imult], multbins[imult+1],6);
					ana[ichg][imult]->SetNumEventsToMix(3);
					ana[ichg][imult]->SetMinSizePartCollection(4);
					ana[ichg][imult]->SetEPhistname(Form("hist%i%i",ichg,imult));
					
					mecetaphitpc[ichg][imult] = new AliFemtoBasicEventCut();
					mecetaphitpc[ichg][imult]->SetEventMult(0.001,100000);
					mecetaphitpc[ichg][imult]->SetVertZPos(-8,8);
					
					cutPassEvMetaphitpc[ichg][imult] = new AliFemtoCutMonitorEventMult(Form("cutPass%stpcM%i", chrgs[ichg], imult));
					cutFailEvMetaphitpc[ichg][imult] = new AliFemtoCutMonitorEventMult(Form("cutFail%stpcM%i", chrgs[ichg], imult));
					mecetaphitpc[ichg][imult]->AddCutMonitor(cutPassEvMetaphitpc[ichg][imult], cutFailEvMetaphitpc[ichg][imult]);
					
					cutPassEvVetaphitpc[ichg][imult] = new AliFemtoCutMonitorEventVertex(Form("cutPass%stpcM%i", chrgs[ichg], imult));
					cutFailEvVetaphitpc[ichg][imult] = new AliFemtoCutMonitorEventVertex(Form("cutFail%stpcM%i", chrgs[ichg], imult));
					mecetaphitpc[ichg][imult]->AddCutMonitor(cutPassEvVetaphitpc[ichg][imult], cutFailEvVetaphitpc[ichg][imult]);
					
					dtc1etaphitpc[ichg][imult] = new AliFemtoESDTrackCut();
					
					if (ichg == 0)
						dtc1etaphitpc[ichg][imult]->SetCharge(1.0);
					else if (ichg == 1)
						dtc1etaphitpc[ichg][imult]->SetCharge(-1.0);
					
					dtc1etaphitpc[ichg][imult]->SetPt(0.15,2.0);
					dtc1etaphitpc[ichg][imult]->SetEta(-0.8,0.8);
					// 	  dtc1etaphitpc->SetPIDMethod(1);
					dtc1etaphitpc[ichg][imult]->SetMass(PionMass);
					//	  dtc1etaphitpc->SetminTPCncls(80);
					
					dtc1etaphitpc[ichg][imult]->SetMostProbablePion();
					
					//	    dtc1etaphitpc->SetStatus(AliESDtrack::kTPCrefit|AliESDtrack::kITSrefit);
					dtc1etaphitpc[ichg][imult]->SetStatus(AliESDtrack::kTPCin);
					//	    dtc1etaphitpc->SetStatus(AliESDtrack::kTPCrefit);
					//    dtc1etaphitpc->SetStatus(AliESDtrack::kITSrefit);
					dtc1etaphitpc[ichg][imult]->SetminTPCncls(0);
					dtc1etaphitpc[ichg][imult]->SetRemoveKinks(kTRUE);
					dtc1etaphitpc[ichg][imult]->SetLabel(kFALSE);
					//    dtc1etaphitpc->SetMaxITSChiNdof(6.0);
					dtc1etaphitpc[ichg][imult]->SetMaxTPCChiNdof(4.0);
					dtc1etaphitpc[ichg][imult]->SetMaxImpactXY(0.20);
					//            dtc1etaphitpc->SetMaxImpactXYPtDep(0.0182, 0.0350, -1.01);
					dtc1etaphitpc[ichg][imult]->SetMaxImpactZ(0.15);
					//      dtc1etaphitpc->SetMaxSigmaToVertex(6.0);
					
					// sqpcetaphitpc = new AliFemtoPairCutAntiGamma();
					//	  sqpcetaphitpc = new AliFemtoShareQualityTPCEntranceSepPairCut();
					sqpcetaphitpc[ichg][imult] = new AliFemtoPairCutRadialDistanceLM();
					sqpcetaphitpcRD[ichg][imult] = new AliFemtoPairCutRadialDistanceLM();
					
					sqpcetaphitpc[ichg][imult]->SetShareQualityMax(1.0);
					sqpcetaphitpc[ichg][imult]->SetShareFractionMax(1);
					sqpcetaphitpc[ichg][imult]->SetRemoveSameLabel(kFALSE);
					sqpcetaphitpc[ichg][imult]->SetMinimumRadius(0.8);
					sqpcetaphitpc[ichg][imult]->SetPhiStarDifferenceMinimum(0.0);
					sqpcetaphitpc[ichg][imult]->SetEtaDifferenceMinimum(0.0);
					
					sqpcetaphitpcRD[ichg][imult]->SetShareQualityMax(1.0);
					sqpcetaphitpcRD[ichg][imult]->SetShareFractionMax(0.05);
					sqpcetaphitpcRD[ichg][imult]->SetRemoveSameLabel(kFALSE);
					sqpcetaphitpcRD[ichg][imult]->SetMinimumRadius(0.8);
					sqpcetaphitpcRD[ichg][imult]->SetPhiStarDifferenceMinimum(0.012);
					sqpcetaphitpcRD[ichg][imult]->SetEtaDifferenceMinimum(0.017);
					
					ana[ichg][imult]->SetEventCut(mecetaphitpc[ichg][imult]);
					ana[ichg][imult]->SetFirstParticleCut(dtc1etaphitpc[ichg][imult]);
					ana[ichg][imult]->SetSecondParticleCut(dtc1etaphitpc[ichg][imult]);
					ana[ichg][imult]->SetPairCut(sqpcetaphitpc[ichg][imult]);
					ana[ichg][imult]->SetPairCutRD(sqpcetaphitpcRD[ichg][imult]);
					
					
					for (int ikt=0; ikt<4; ikt++){
						for (int iphi=0; iphi<6; iphi++){
							ktpaircut[ichg][imult][ikt][iphi] = new AliFemtoKTPairCut(ktrng[ikt],ktrng[ikt+1]);
							ktpaircut[ichg][imult][ikt][iphi]->SetPhiRange(phirange[iphi],phirange[iphi+1]);
							
							// 		cylmkt[ich][imult][ikt][iphi] = new AliFemtoCorrFctnDirectYlm(Form("cylm%imult%ikT%iRP%i", ich, imult, ikt, iphi),2,90,0.0,0.3,1);
							// 		cylmkt[ich][imult][ikt][iphi]->SetPairSelectionCut(ktpaircut[ich][imult][ikt][iphi]);
							// 		ana[ich][imult]->AddCorrFctn(cylmkt[ich][imult][ikt][iphi]);
							
							// 		cqinvkt[ich][imult][ikt][iphi] = new AliFemtoQinvCorrFctn(Form("cqinv%imult%ikT%iRP%i", ich, imult, ikt, iphi),20,0.0,0.2);
							// 		cqinvkt[ich][imult][ikt][iphi]->SetPairSelectionCut(ktpaircut[ich][imult][ikt][iphi]);
							// 		ana[ich][imult]->AddCorrFctn(cqinvkt[ich][imult][ikt][iphi]);
							
							cq3dlcmskttpc[ichg][imult][ikt][iphi] = new AliFemtoBPLCMS3DCorrFctn(Form("cq3d%imult%ikT%iRP%i", ichg, imult, ikt, iphi),30,-0.15,0.15);
							cq3dlcmskttpc[ichg][imult][ikt][iphi]->SetPairSelectionCut(ktpaircut[ichg][imult][ikt][iphi]);
							ana[ichg][imult]->AddCorrFctn(cq3dlcmskttpc[ichg][imult][ikt][iphi]);
						}
					}
					Manager->AddAnalysis(ana[ichg][imult]);	
				}
			}
		}
	}
	return Manager;
}                         

