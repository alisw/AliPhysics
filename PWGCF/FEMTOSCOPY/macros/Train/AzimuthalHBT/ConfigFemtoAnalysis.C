
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
#include "AliFemtoCutMonitorCollections.h"

#endif

//________________________________________________________________________
AliFemtoManager* ConfigFemtoAnalysis(int getbinwidth=20) {
	
	double PionMass = 0.13956995;
	double KaonMass = 0.493677;
	
	// int runmults[10] = {1, 1, 0, 0, 0, 0, 0, 0, 0, 0};
	int runmults[10] = {1, 1, 1, 1, 1, 1, 1, 0, 0, 0};
	
	//int runmults[10] = {1, 1, 1, 1, 1, 1, 1, 0, 0, 0};
	
	
	int multbins[11] = {0, 50, 100, 200, 300, 400, 500, 700, 800, 850, 900};
	
	int runch[2] = {1, 1};
	const char *chrgs[2] = { "pip", "pim" };
	
	int runktdep = 1;
	double ktrng[5] = {0.2, 0.3, 0.4, 0.5, 0.7};
	//	double ktrng[7] = {0.2,0.3, 0.4, 0.5,0.6,0.7,0.8};
	
	//  int phirange[7] = {-15, 15, 45, 75, 105, 135, 165};
	int phirange[10] = {-15, 5, 25,  45,65,85,105,125,145,165};
	
	
	int runtype = 2; // Types 0 - global, 1 - ITS only, 2 - TPC Inner
	int isrealdata = 0;
	
	AliFemtoEventReaderAODChain *Reader = new AliFemtoEventReaderAODChain();
	Reader->SetFilterBit(7);
	Reader->SetCentralityPreSelection(0.000001, 701);
	
	
	AliFemtoManager* Manager=new AliFemtoManager();
	Manager->SetEventReader(Reader);
	
	AliFemtoAnalysisAzimuthalPbPb   *ana[2][10];
	AliFemtoVertexMultAnalysis    *anetaphitpc[320];
	
	AliFemtoBasicEventCut         *mecetaphitpc[320];
	AliFemtoCutMonitorEventMult   *cutPassEvMetaphitpc[320];
	AliFemtoCutMonitorEventMult   *cutFailEvMetaphitpc[320];
	AliFemtoCutMonitorCollections *cutPassColletaphitpc[320];
	AliFemtoCutMonitorCollections *cutFailColletaphitpc[320];
	
	AliFemtoCutMonitorEventVertex *cutPassEvVetaphitpc[320];
	AliFemtoCutMonitorEventVertex *cutFailEvVetaphitpc[320];
	
	AliFemtoESDTrackCut           *dtc1etaphitpc[2][10];
	AliFemtoESDTrackCut           *dtc2etaphitpc[2][10];
	
	AliFemtoPairCutRadialDistanceLM      *sqpcetaphitpc[2][10];
	AliFemtoPairCutRadialDistanceLM      *sqpcetaphitpcRD[2][10];
	AliFemtoPairCutRadialDistanceLM      *sqpc3etaphitpc[320];
	AliFemtoPairCutRadialDistanceLM      *sqpc3etaphitpcRD[320];
	
	AliFemtoKTPairCut             *ktpaircut[2][10][6][8];
	
	AliFemtoBPLCMS3DCorrFctn     *cq3dlcmskttpc[2][10][6][8];
	
	AliFemtoESDTrackCut           *dtc3etaphitpc[320];
	AliFemtoESDTrackCut           *dtc4etaphitpc[320];
	
	AliFemtoPairCutAntiGamma      *sqp3cetaphitpc[320];
	AliFemtoCorrFctnDirectYlm     *cylmetaphitpc[320];
	AliFemtoCorrFctnDEtaDPhi      *cdedpetaphi[320*6];
	AliFemtoChi2CorrFctn          *cchiqinvetaphitpc[320];
	AliFemtoKTPairCut             *ktpcuts[320*6];
	AliFemtoCorrFctnDirectYlm     *cylmkttpc[320*6];
	AliFemtoQinvCorrFctn          *cqinvkttpc[320*6];
	AliFemtoCorrFctnTPCNcls       *cqinvnclstpc[320];
	AliFemtoShareQualityCorrFctn  *cqinvsqtpc[320];
	AliFemtoTPCInnerCorrFctn      *cqinvtitpc[320];
	AliFemtoQinvCorrFctn          *cqinvtpc[100*3];
	
	int numOfMultBins=10;
	int aniter = 0;
	int isrealdata = 0; //AOD = 0 (doesn't change anything)
	
	
	
	for (int imult=0; imult<10; imult++) {
		if (runmults[imult]) {
			for (int ichg=0; ichg<2; ichg++) {
				if (runch[ichg]) {
					aniter = ichg*numOfMultBins+imult;
					
					anetaphitpc[aniter] = new AliFemtoVertexMultAnalysis(8, -8.0, 8.0, 4, multbins[imult], multbins[imult+1]);
					anetaphitpc[aniter]->SetNumEventsToMix(3);
					anetaphitpc[aniter]->SetMinSizePartCollection(4);
					
					ana[ichg][imult] = new AliFemtoAnalysisAzimuthalPbPb(4, -8.0, 8.0, 4, multbins[imult], multbins[imult+1],9);
					
					ana[ichg][imult]->SetNumEventsToMix(3);
					ana[ichg][imult]->SetMinSizePartCollection(4);
					ana[ichg][imult]->SetEPhistname(Form("hist%i%i",ichg,imult));
					
					mecetaphitpc[aniter] = new AliFemtoBasicEventCut();
					mecetaphitpc[aniter]->SetEventMult(10,100000);  //remove 0 events
					mecetaphitpc[aniter]->SetVertZPos(-8,8);
					
					cutPassEvMetaphitpc[aniter] = new AliFemtoCutMonitorEventMult(Form("cutPass%stpcM%i", chrgs[ichg], imult));
					cutFailEvMetaphitpc[aniter] = new AliFemtoCutMonitorEventMult(Form("cutFail%stpcM%i", chrgs[ichg], imult));
					mecetaphitpc[aniter]->AddCutMonitor(cutPassEvMetaphitpc[aniter], cutFailEvMetaphitpc[aniter]);
					
					cutPassEvVetaphitpc[aniter] = new AliFemtoCutMonitorEventVertex(Form("cutPass%stpcM%i", chrgs[ichg], imult));
					cutFailEvVetaphitpc[aniter] = new AliFemtoCutMonitorEventVertex(Form("cutFail%stpcM%i", chrgs[ichg], imult));
					mecetaphitpc[aniter]->AddCutMonitor(cutPassEvVetaphitpc[aniter], cutFailEvVetaphitpc[aniter]);
					
					cutPassColletaphitpc[aniter] = new AliFemtoCutMonitorCollections(Form("cutPass%stpcM%i", chrgs[ichg], imult));
					cutFailColletaphitpc[aniter] = new AliFemtoCutMonitorCollections(Form("cutFail%stpcM%i", chrgs[ichg], imult));
					mecetaphitpc[aniter]->AddCutMonitor(cutPassColletaphitpc[aniter], cutFailColletaphitpc[aniter]);
					
					
					//ESD first particle cut -> Pion 3, 5; AntiPion 4, 6
					dtc3etaphitpc[aniter] = new AliFemtoESDTrackCut();
					dtc3etaphitpc[aniter]->SetMostProbablePion();
					dtc3etaphitpc[aniter]->SetMass(PionMass);
					if(ichg == 0 || ichg == 1)
						dtc3etaphitpc[aniter]->SetCharge(1.0);
					else
						dtc3etaphitpc[aniter]->SetCharge(-1.0);
					
					dtc3etaphitpc[aniter]->SetPt(0.15,2.0);
					dtc3etaphitpc[aniter]->SetEta(-0.8,0.8);
					
					
					// Track quality cuts	 
					dtc3etaphitpc[aniter]->SetStatus(AliESDtrack::kTPCin);
					dtc3etaphitpc[aniter]->SetminTPCncls(80);
					dtc3etaphitpc[aniter]->SetRemoveKinks(kTRUE);
					dtc3etaphitpc[aniter]->SetLabel(kFALSE);
					dtc3etaphitpc[aniter]->SetMaxTPCChiNdof(4.0);
					dtc3etaphitpc[aniter]->SetMaxImpactXY(2.4);
					dtc3etaphitpc[aniter]->SetMaxImpactZ(3.0);
					
					//ESD first particle cut -> Pion 3, 5; AntiPion 4, 6, 8, 9
					dtc4etaphitpc[aniter] = new AliFemtoESDTrackCut();
					dtc4etaphitpc[aniter]->SetMostProbablePion();
					dtc4etaphitpc[aniter]->SetMass(PionMass);
					dtc4etaphitpc[aniter]->SetCharge(-1.0);
					dtc4etaphitpc[aniter]->SetPt(0.15,2.0);
					dtc4etaphitpc[aniter]->SetEta(-0.8,0.8);
					
					// Track quality cuts
					dtc4etaphitpc[aniter]->SetStatus(AliESDtrack::kTPCin);
					dtc4etaphitpc[aniter]->SetminTPCncls(80);
					dtc4etaphitpc[aniter]->SetRemoveKinks(kTRUE);
					dtc4etaphitpc[aniter]->SetLabel(kFALSE);
					dtc4etaphitpc[aniter]->SetMaxTPCChiNdof(4.0);
					dtc4etaphitpc[aniter]->SetMaxImpactXY(2.4);
					dtc4etaphitpc[aniter]->SetMaxImpactZ(3.0);
					
					
					
					dtc1etaphitpc[ichg][imult] = new AliFemtoESDTrackCut();
					
					if (ichg == 0)
						dtc1etaphitpc[ichg][imult]->SetCharge(1.0);
					else if (ichg == 1)
						dtc1etaphitpc[ichg][imult]->SetCharge(-1.0);
					
					dtc1etaphitpc[ichg][imult]->SetPt(0.15,2.0);
					dtc1etaphitpc[ichg][imult]->SetEta(-0.8,0.8);
					dtc1etaphitpc[ichg][imult]->SetMass(PionMass);
					dtc1etaphitpc[ichg][imult]->SetMostProbablePion();
					dtc1etaphitpc[ichg][imult]->SetStatus(AliESDtrack::kTPCin);
					dtc1etaphitpc[ichg][imult]->SetminTPCncls(80);
					dtc1etaphitpc[ichg][imult]->SetRemoveKinks(kTRUE);
					dtc1etaphitpc[ichg][imult]->SetLabel(kFALSE);
					dtc1etaphitpc[ichg][imult]->SetMaxTPCChiNdof(4.0);
					dtc1etaphitpc[ichg][imult]->SetMaxImpactXY(2.4);
					dtc1etaphitpc[ichg][imult]->SetMaxImpactZ(3.0);
					
					sqpcetaphitpc[ichg][imult] = new AliFemtoPairCutRadialDistanceLM();
					sqpcetaphitpcRD[ichg][imult] = new AliFemtoPairCutRadialDistanceLM();
					
					sqpc3etaphitpc[aniter]=new AliFemtoPairCutRadialDistanceLM();
					sqpc3etaphitpcRD[aniter]=new AliFemtoPairCutRadialDistanceLM();
					
					
					sqpcetaphitpc[ichg][imult]->SetShareQualityMax(1.0);
					sqpcetaphitpc[ichg][imult]->SetShareFractionMax(1);
					sqpcetaphitpc[ichg][imult]->SetRemoveSameLabel(kFALSE);
					sqpcetaphitpc[ichg][imult]->SetMinimumRadius(0.8);
					sqpcetaphitpc[ichg][imult]->SetPhiStarDifferenceMinimum(0.0);
					sqpcetaphitpc[ichg][imult]->SetEtaDifferenceMinimum(0.0);
					
					sqpc3etaphitpc[aniter]->SetShareQualityMax(1.0);
					sqpc3etaphitpc[aniter]->SetShareFractionMax(1);
					sqpc3etaphitpc[aniter]->SetRemoveSameLabel(kFALSE);
					sqpc3etaphitpc[aniter]->SetMinimumRadius(0.8);
					sqpc3etaphitpc[aniter]->SetPhiStarDifferenceMinimum(0.0);
					sqpc3etaphitpc[aniter]->SetEtaDifferenceMinimum(0.0);
					
					sqpcetaphitpcRD[ichg][imult]->SetShareQualityMax(1.0);
					sqpcetaphitpcRD[ichg][imult]->SetShareFractionMax(0.05);
					sqpcetaphitpcRD[ichg][imult]->SetRemoveSameLabel(kFALSE);
					sqpcetaphitpcRD[ichg][imult]->SetMinimumRadius(0.8);
					sqpcetaphitpcRD[ichg][imult]->SetPhiStarDifferenceMinimum(0.012);
					sqpcetaphitpcRD[ichg][imult]->SetEtaDifferenceMinimum(0.017);
					
					sqpc3etaphitpcRD[aniter]->SetShareQualityMax(1.0);
					sqpc3etaphitpcRD[aniter]->SetShareFractionMax(0.05);
					sqpc3etaphitpcRD[aniter]->SetRemoveSameLabel(kFALSE);
					sqpc3etaphitpcRD[aniter]->SetMinimumRadius(0.8);
					sqpc3etaphitpcRD[aniter]->SetPhiStarDifferenceMinimum(0.012);
					sqpc3etaphitpcRD[aniter]->SetEtaDifferenceMinimum(0.017);
					
					if(ichg == 0) //p-p
					{
						anetaphitpc[aniter]->SetEventCut(mecetaphitpc[aniter]);
						anetaphitpc[aniter]->SetFirstParticleCut(dtc3etaphitpc[aniter]);
						anetaphitpc[aniter]->SetSecondParticleCut(dtc3etaphitpc[aniter]);
						anetaphitpc[aniter]->SetPairCut(sqpc3etaphitpc[aniter]);
						anetaphitpc[aniter]->SetPairCut(sqpc3etaphitpcRD[aniter]);
						
					}
					else if(ichg == 1) //ap-ap
					{
						anetaphitpc[aniter]->SetEventCut(mecetaphitpc[aniter]);
						anetaphitpc[aniter]->SetFirstParticleCut(dtc4etaphitpc[aniter]);
						anetaphitpc[aniter]->SetSecondParticleCut(dtc4etaphitpc[aniter]);
						anetaphitpc[aniter]->SetPairCut(sqpc3etaphitpc[aniter]);
						anetaphitpc[aniter]->SetPairCut(sqpc3etaphitpcRD[aniter]);
						
					}
					
					ana[ichg][imult]->SetEventCut(mecetaphitpc[aniter]);
					ana[ichg][imult]->SetFirstParticleCut(dtc1etaphitpc[ichg][imult]);
					ana[ichg][imult]->SetSecondParticleCut(dtc1etaphitpc[ichg][imult]);
					ana[ichg][imult]->SetPairCut(sqpcetaphitpc[ichg][imult]);
					ana[ichg][imult]->SetPairCutRD(sqpcetaphitpcRD[ichg][imult]);
					
					cdedpetaphi[aniter] = new AliFemtoCorrFctnDEtaDPhi(Form("cdedp%stpcM%i", chrgs[ichg], imult),45, 45);
					anetaphitpc[aniter]->AddCorrFctn(cdedpetaphi[aniter]);
					
					cqinvtpc[aniter] = new AliFemtoQinvCorrFctn(Form("cqinv%stpcM%i", chrgs[ichg], imult),100,0.0,0.5);
					anetaphitpc[aniter]->AddCorrFctn(cqinvtpc[aniter]);
					
					for (int ikt=0; ikt<4; ikt++){
						for (int iphi=0; iphi<9; iphi++){
							
							ktpaircut[ichg][imult][ikt][iphi] = new AliFemtoKTPairCut(ktrng[ikt],ktrng[ikt+1]);
							ktpaircut[ichg][imult][ikt][iphi]->SetPhiRange(phirange[iphi],phirange[iphi+1]);
							
							
							//cq3dlcmskttpc[ichg][imult][ikt][iphi] = new AliFemtoBPLCMS3DCorrFctn(Form("cq3d%imult%ikT%iRP%i", ichg, imult, ikt, iphi),30,-0.15,0.15);
							cq3dlcmskttpc[ichg][imult][ikt][iphi] = new AliFemtoBPLCMS3DCorrFctn(Form("cq3d%imult%ikT%iRP%i", ichg, imult, ikt, iphi),getbinwidth,-0.2,0.2);
							
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

