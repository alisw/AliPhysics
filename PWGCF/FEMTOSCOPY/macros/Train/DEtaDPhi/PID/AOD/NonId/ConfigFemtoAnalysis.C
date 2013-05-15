/*********************************************************************
 *																							 *
 * Configfemtoanalysis.C - configuration macro for the femtoscopic	 *
 * analysis, meant as a QA process for two-particle effects				 *
 *																							 *
 * Author: Adam Kisiel (Adam.Kisiel@cern.ch)									 *
 *																							 *
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
//#include "AliFemtoCutMonitorParticleVertPos.h"
//#include "AliFemtoCutMonitorParticleMomRes.h"
#include "AliFemtoCutMonitorParticlePID.h"
//#include "AliFemtoCutMonitorEventMult.h"
//#include "AliFemtoCutMonitorEventVertex.h"
#include "AliFemtoShareQualityTPCEntranceSepPairCut.h"
#include "AliFemtoPairCutAntiGamma.h"
#include "AliFemtoPairCutRadialDistance.h"
#include "AliFemtoQinvCorrFctn.h"
//#include "AliFemtoTPCInnerCorrFctn.h"
#include "AliFemtoVertexMultAnalysis.h"
//#include "AliFemtoCutMonitorParticlePtPDG.h"
#include "AliFemtoPairCutPt.h"
//#include "AliFemtoCorrFctnMassInvMonitor.h"
#endif

//________________________________________________________________________
AliFemtoManager* ConfigFemtoAnalysis() {

	double PionMass = 0.13956995;
	double KaonMass = 0.493677;
	double ProtonMass = 0.938272013;

	const int numOfMultBins = 5;	
	const int numOfChTypes = 12;
	const int numOfkTbins = 2;

	int runmults[numOfMultBins] = {1, 1, 1, 0, 1};
	int multbins[numOfMultBins+1] = {0, 20, 50, 150, 2, 150};
	
	int runch[numOfChTypes] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
	const char *chrgs[numOfChTypes] = { "PIpKp", "PImKm", "PIpKm", "PImKp", "PIpPp","PImPm","PIpPm","PImPp", "KpPp","KmPm","KpPm","KmPp"};
	
	int runktdep = 0;
	double ktrng[numOfkTbins+1] = {0.0, 0.7, 100.0};

	int runtype = 0; // Types 0 - global, 1 - ITS only, 2 - TPC Inner	//global tracks ->mfit ITS+TPC

	int gammacut = 0;	// cut na ee z gamma 
	

	//AliFemtoEventReaderESDChainKine* Reader=new AliFemtoEventReaderESDChainKine();
	//Reader->SetUseMultiplicity(AliFemtoEventReaderESDChainKine::kReferenceITSTPC);
	//Reader->SetMagneticFieldSign(1.0);
	//AliFemtoEventReaderESDChain *Reader = new AliFemtoEventReaderESDChain();
	//Reader->SetUseMultiplicity(AliFemtoEventReaderESDChain::kReferenceITSTPC);

	AliFemtoEventReaderAODChain *Reader = new AliFemtoEventReaderAODChain();
	Reader->SetFilterBit(0);
	Reader->SetUseMultiplicity(AliFemtoEventReaderAOD::kReference);

	// //Reader->SetCentralityPreSelection(0.001, 910);

	AliFemtoManager* Manager = new AliFemtoManager();
	Manager->SetEventReader(Reader);

	AliFemtoVertexMultAnalysis		*anetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoBasicEventCut				 *mecetaphitpc[numOfMultBins*numOfChTypes];
/*	AliFemtoCutMonitorEventMult	 *cutPassEvMetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoCutMonitorEventMult	 *cutFailEvMetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoCutMonitorEventVertex *cutPassEvVetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoCutMonitorEventVertex *cutFailEvVetaphitpc[numOfMultBins*numOfChTypes];
*/	AliFemtoESDTrackCut					 *dtc1etaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoESDTrackCut					 *dtc2etaphitpc[numOfMultBins*numOfChTypes];

	//	 AliFemtoPairCutAntiGamma			*sqpcetaphitpcdiff[numOfMultBins*numOfChTypes];
	//	 AliFemtoShareQualityTPCEntranceSepPairCut			*sqpcetaphitpcsame[numOfMultBins*numOfChTypes];
	AliFemtoPairCutAntiGamma			*sqpcetaphitpc[numOfMultBins*numOfChTypes];
	//	AliFemtoPairCutRadialDistance			*sqpcetaphitpc[numOfMultBins*numOfChTypes];
	//	AliFemtoChi2CorrFctn					*cchiqinvetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoPairCutPt						 *ktpcuts[numOfMultBins*numOfChTypes];
	AliFemtoQinvCorrFctn					*cqinvkttpc[numOfMultBins*numOfChTypes];
	AliFemtoQinvCorrFctn					*cqinvtpc[numOfMultBins*numOfChTypes];
	AliFemtoCorrFctnDEtaDPhi			*cdedpetaphi[numOfMultBins*numOfChTypes];
	//AliFemtoCorrFctnMassInvMonitor			*cMinvMonitor[numOfMultBins*numOfChTypes];

	
	// *** Third QA task - HBT analysis with all pair cuts off, TPC only ***
	// *** Begin pion-pion (positive) analysis ***
	int aniter = 0;	

	for (int imult = 0; imult < numOfMultBins; imult++)
	{
		if (runmults[imult])
		{
			for (int ichg = 0; ichg < numOfChTypes; ichg++)
			{
				if (runch[ichg])
				{

					aniter = ichg * numOfMultBins + imult;
					anetaphitpc[aniter] = new AliFemtoVertexMultAnalysis(8, -8.0, 8.0, 6, multbins[imult], multbins[imult+1]);
					anetaphitpc[aniter]->SetNumEventsToMix(10);
					anetaphitpc[aniter]->SetMinSizePartCollection(1);

					mecetaphitpc[aniter] = new AliFemtoBasicEventCut();
					mecetaphitpc[aniter]->SetEventMult(0.001,100000);
					mecetaphitpc[aniter]->SetVertZPos(-8,8);//cm
				
			/*		cutPassEvMetaphitpc[aniter] = new AliFemtoCutMonitorEventMult(Form("cutPass%stpcM%i", chrgs[ichg], imult));
					cutFailEvMetaphitpc[aniter] = new AliFemtoCutMonitorEventMult(Form("cutFail%stpcM%i", chrgs[ichg], imult));
					mecetaphitpc[aniter]->AddCutMonitor(cutPassEvMetaphitpc[aniter], cutFailEvMetaphitpc[aniter]);
		
					cutPassEvVetaphitpc[aniter] = new AliFemtoCutMonitorEventVertex(Form("cutPass%stpcM%i", chrgs[ichg], imult));
					cutFailEvVetaphitpc[aniter] = new AliFemtoCutMonitorEventVertex(Form("cutFail%stpcM%i", chrgs[ichg], imult));
					mecetaphitpc[aniter]->AddCutMonitor(cutPassEvVetaphitpc[aniter], cutFailEvVetaphitpc[aniter]);
		*/
					dtc1etaphitpc[aniter] = new AliFemtoESDTrackCut();
					dtc2etaphitpc[aniter] = new AliFemtoESDTrackCut();

					dtc1etaphitpc[aniter]->SetEta(-1.0,1.0);
					dtc2etaphitpc[aniter]->SetEta(-1.0,1.0);
							
					if (ichg == 0) // PIpKp
					{
						dtc1etaphitpc[aniter]->SetCharge(1.0);
						dtc1etaphitpc[aniter]->SetPt(0.2,20);
						dtc1etaphitpc[aniter]->SetMass(PionMass);		
						dtc1etaphitpc[aniter]->SetMostProbablePion();
						dtc2etaphitpc[aniter]->SetCharge(1.0);
						dtc2etaphitpc[aniter]->SetPt(0.3,20);
						dtc2etaphitpc[aniter]->SetMass(KaonMass);		
						dtc2etaphitpc[aniter]->SetMostProbableKaon();
					}		
					if (ichg == 1) // PImKm
					{
						dtc1etaphitpc[aniter]->SetCharge(-1.0);
						dtc1etaphitpc[aniter]->SetPt(0.2,20);
						dtc1etaphitpc[aniter]->SetMass(PionMass);		
						dtc1etaphitpc[aniter]->SetMostProbablePion();
						dtc2etaphitpc[aniter]->SetCharge(-1.0);
						dtc2etaphitpc[aniter]->SetPt(0.3,20);
						dtc2etaphitpc[aniter]->SetMass(KaonMass);		
						dtc2etaphitpc[aniter]->SetMostProbableKaon();
					}		
					if (ichg == 2) //PIpKm
					{
						dtc1etaphitpc[aniter]->SetCharge(1.0);
						dtc1etaphitpc[aniter]->SetPt(0.2,20);
						dtc1etaphitpc[aniter]->SetMass(PionMass);		
						dtc1etaphitpc[aniter]->SetMostProbablePion();
						dtc2etaphitpc[aniter]->SetCharge(-1.0);
						dtc2etaphitpc[aniter]->SetPt(0.3,20);
						dtc2etaphitpc[aniter]->SetMass(KaonMass);		
						dtc2etaphitpc[aniter]->SetMostProbableKaon();
					}
					if (ichg == 3) // PImKp
					{
						dtc1etaphitpc[aniter]->SetCharge(-1.0);
						dtc1etaphitpc[aniter]->SetPt(0.2,20);
						dtc1etaphitpc[aniter]->SetMass(PionMass);		
						dtc1etaphitpc[aniter]->SetMostProbablePion();
						dtc2etaphitpc[aniter]->SetCharge(1.0);
						dtc2etaphitpc[aniter]->SetPt(0.3,20);
						dtc2etaphitpc[aniter]->SetMass(KaonMass);		
						dtc2etaphitpc[aniter]->SetMostProbableKaon();
					}
					if (ichg == 4) // PIpPp
					{
						dtc1etaphitpc[aniter]->SetCharge(1.0);
						dtc1etaphitpc[aniter]->SetPt(0.2,20);
						dtc1etaphitpc[aniter]->SetMass(PionMass);		
						dtc1etaphitpc[aniter]->SetMostProbablePion();
						dtc2etaphitpc[aniter]->SetCharge(1.0);
						dtc2etaphitpc[aniter]->SetPt(0.5,20);
						dtc2etaphitpc[aniter]->SetMass(ProtonMass);		
						dtc2etaphitpc[aniter]->SetMostProbableProton();
					}
					if (ichg == 5) // PImPm
					{
						dtc1etaphitpc[aniter]->SetCharge(-1.0);
						dtc1etaphitpc[aniter]->SetPt(0.2,20);
						dtc1etaphitpc[aniter]->SetMass(PionMass);		
						dtc1etaphitpc[aniter]->SetMostProbablePion();
						dtc2etaphitpc[aniter]->SetCharge(-1.0);
						dtc2etaphitpc[aniter]->SetPt(0.5,20);
						dtc2etaphitpc[aniter]->SetMass(ProtonMass);		
						dtc2etaphitpc[aniter]->SetMostProbableProton();
					}
					if (ichg == 6) // PIpPm
					{
						dtc1etaphitpc[aniter]->SetCharge(1.0);
						dtc1etaphitpc[aniter]->SetPt(0.2,20);
						dtc1etaphitpc[aniter]->SetMass(PionMass);		
						dtc1etaphitpc[aniter]->SetMostProbablePion();
						dtc2etaphitpc[aniter]->SetCharge(-1.0);
						dtc2etaphitpc[aniter]->SetPt(0.5,20);
						dtc2etaphitpc[aniter]->SetMass(ProtonMass);		
						dtc2etaphitpc[aniter]->SetMostProbableProton();
					}
					if (ichg == 7) // PImPp
					{
						dtc1etaphitpc[aniter]->SetCharge(-1.0);
						dtc1etaphitpc[aniter]->SetPt(0.2,20);
						dtc1etaphitpc[aniter]->SetMass(PionMass);		
						dtc1etaphitpc[aniter]->SetMostProbablePion();
						dtc2etaphitpc[aniter]->SetCharge(1.0);
						dtc2etaphitpc[aniter]->SetPt(0.5,20);
						dtc2etaphitpc[aniter]->SetMass(ProtonMass);		
						dtc2etaphitpc[aniter]->SetMostProbableProton();
					}
					if (ichg == 8) // KpPp
					{
						dtc1etaphitpc[aniter]->SetCharge(1.0);
						dtc1etaphitpc[aniter]->SetPt(0.3,20);
						dtc1etaphitpc[aniter]->SetMass(KaonMass);		
						dtc1etaphitpc[aniter]->SetMostProbableKaon();
						dtc2etaphitpc[aniter]->SetCharge(1.0);
						dtc2etaphitpc[aniter]->SetPt(0.5,20);
						dtc2etaphitpc[aniter]->SetMass(ProtonMass);		
						dtc2etaphitpc[aniter]->SetMostProbableProton();
					}
					if (ichg == 9) // KpPp
					{
						dtc1etaphitpc[aniter]->SetCharge(1.0);
						dtc1etaphitpc[aniter]->SetPt(0.3,20);
						dtc1etaphitpc[aniter]->SetMass(KaonMass);		
						dtc1etaphitpc[aniter]->SetMostProbableKaon();
						dtc2etaphitpc[aniter]->SetCharge(1.0);
						dtc2etaphitpc[aniter]->SetPt(0.5,20);
						dtc2etaphitpc[aniter]->SetMass(ProtonMass);		
						dtc2etaphitpc[aniter]->SetMostProbableProton();
					}
					if (ichg == 10) // KpPp
					{
						dtc1etaphitpc[aniter]->SetCharge(1.0);
						dtc1etaphitpc[aniter]->SetPt(0.3,20);
						dtc1etaphitpc[aniter]->SetMass(KaonMass);		
						dtc1etaphitpc[aniter]->SetMostProbableKaon();
						dtc2etaphitpc[aniter]->SetCharge(1.0);
						dtc2etaphitpc[aniter]->SetPt(0.5,20);
						dtc2etaphitpc[aniter]->SetMass(ProtonMass);		
						dtc2etaphitpc[aniter]->SetMostProbableProton();
					}
					if (ichg ==11) // KpPp
					{
						dtc1etaphitpc[aniter]->SetCharge(1.0);
						dtc1etaphitpc[aniter]->SetPt(0.3,20);
						dtc1etaphitpc[aniter]->SetMass(KaonMass);		
						dtc1etaphitpc[aniter]->SetMostProbableKaon();
						dtc2etaphitpc[aniter]->SetCharge(1.0);
						dtc2etaphitpc[aniter]->SetPt(0.5,20);
						dtc2etaphitpc[aniter]->SetMass(ProtonMass);		
						dtc2etaphitpc[aniter]->SetMostProbableProton();
					}

				// Track quality cuts

					if (runtype == 0)
					{
						dtc1etaphitpc[aniter]->SetStatus(AliESDtrack::kTPCrefit|AliESDtrack::kITSrefit);
						dtc1etaphitpc[aniter]->SetminTPCncls(70);
						dtc1etaphitpc[aniter]->SetRemoveKinks(kTRUE);
						dtc1etaphitpc[aniter]->SetLabel(kFALSE);
						//	dtc1etaphitpc[aniter]->SetMaxITSChiNdof(6.0);
						dtc1etaphitpc[aniter]->SetMaxTPCChiNdof(4.0);
						dtc1etaphitpc[aniter]->SetMaxImpactXYPtDep(0.018, 0.035, -1.01); 	//	DCA xy
						dtc1etaphitpc[aniter]->SetMaxImpactZ(2);	//DCA Z
						//	dtc1etaphitpc[aniter]->SetMaxSigmaToVertex(6.0);
						dtc2etaphitpc[aniter]->SetStatus(AliESDtrack::kTPCrefit|AliESDtrack::kITSrefit);
						dtc2etaphitpc[aniter]->SetminTPCncls(70);
						dtc2etaphitpc[aniter]->SetRemoveKinks(kTRUE);
						dtc2etaphitpc[aniter]->SetLabel(kFALSE);
						//	dtc2etaphitpc[aniter]->SetMaxITSChiNdof(6.0);
						dtc2etaphitpc[aniter]->SetMaxTPCChiNdof(4.0);
						dtc2etaphitpc[aniter]->SetMaxImpactXYPtDep(0.018, 0.035, -1.01) ;
						dtc2etaphitpc[aniter]->SetMaxImpactZ(2);
						//	dtc2etaphitpc[aniter]->SetMaxSigmaToVertex(6.0);
					}
					else if (runtype == 1)
					{
						;
					}
					else if (runtype == 2)
					{
						;
					}

					sqpcetaphitpc[aniter] = new AliFemtoPairCutAntiGamma();

					if (runtype == 0)
					{
						sqpcetaphitpc[aniter]->SetShareQualityMax(1.0);		// two track cuts on splitting and merging
						sqpcetaphitpc[aniter]->SetShareFractionMax(0.05);	//  ile moga miec wspolnych klastrow
						sqpcetaphitpc[aniter]->SetRemoveSameLabel(kFALSE);
						// sqpcetaphitpc[aniter]->SetMaxEEMinv(0.0);
						// sqpcetaphitpc[aniter]->SetMaxThetaDiff(0.0);
						//	sqpcetaphitpc[aniter]->SetTPCEntranceSepMinimum(1.5);
						// sqpcetaphitpc[aniter]->SetRadialDistanceMinimum(0.12, 0.03);
						//	sqpcetaphitpc[aniter]->SetEtaDifferenceMinimum(0.02);
					}
					else if (runtype == 1)
					{
						
					}
					else if (runtype == 2)
					{
						
					}
		
					anetaphitpc[aniter]->SetEventCut(mecetaphitpc[aniter]);

					anetaphitpc[aniter]->SetFirstParticleCut(dtc1etaphitpc[aniter]);
					anetaphitpc[aniter]->SetSecondParticleCut(dtc2etaphitpc[aniter]);
						//cMinvMonitor[aniter] = new AliFemtoCorrFctnMassInvMonitor(Form("%sM%i",chrgs[ichg], imult),500,ProtonMass);
						//anetaphitpc[aniter]->AddCorrFctn(cMinvMonitor[aniter]);
					
					anetaphitpc[aniter]->SetPairCut(sqpcetaphitpc[aniter]);
					cdedpetaphi[aniter] = new AliFemtoCorrFctnDEtaDPhi(Form("cdedp%stpcM%i", chrgs[ichg], imult),35, 35);
					anetaphitpc[aniter]->AddCorrFctn(cdedpetaphi[aniter]);

					if (runktdep)
					{
						int ktm;
						for (int ikt=0; ikt<numOfkTbins; ikt++)
						{
							ktm = aniter * numOfkTbins + ikt;
							ktpcuts[ktm] = new AliFemtoPairCutPt(ktrng[ikt], ktrng[ikt+1]);
				
							//cqinvkttpc[ktm] = new AliFemtoQinvCorrFctn(Form("cqinv%stpcM%ikT%i", chrgs[ichg], imult, ikt),nbinssh,0.0,(imult>6)?shqmax*2.5:shqmax);
							//cqinvkttpc[ktm]->SetPairSelectionCut(ktpcuts[ktm]);
							//anetaphitpc[aniter]->AddCorrFctn(cqinvkttpc[ktm]);

							cdedpetaphi[ktm] = new AliFemtoCorrFctnDEtaDPhi(Form("cdedp%stpcM%ipT%i", chrgs[ichg], imult,ikt),35, 35);
							cdedpetaphi[ktm]->SetPairSelectionCut(ktpcuts[ktm]);
							anetaphitpc[aniter]->AddCorrFctn(cdedpetaphi[ktm]);
						}
					}		
					Manager->AddAnalysis(anetaphitpc[aniter]);	
				}
			}
		}
	}

	return Manager;
}												 
