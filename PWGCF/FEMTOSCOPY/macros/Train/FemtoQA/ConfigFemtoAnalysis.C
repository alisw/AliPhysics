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
#include "AliFemtoCutMonitorParticleVertPos.h"
#include "AliFemtoCutMonitorParticleMomRes.h"
#include "AliFemtoCutMonitorParticlePID.h"
#include "AliFemtoCutMonitorEventMult.h"
#include "AliFemtoCutMonitorEventVertex.h"
#include "AliFemtoShareQualityTPCEntranceSepPairCut.h"
#include "AliFemtoPairCutAntiGamma.h"
#include "AliFemtoPairCutRadialDistance.h"
#include "AliFemtoQinvCorrFctn.h"
#include "AliFemtoCorrFctnNonIdDR.h"
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
#include "AliFemtoPairCutPt.h"
#endif

//________________________________________________________________________
AliFemtoManager* ConfigFemtoAnalysis() {

	double PionMass = 0.13956995;
	double KaonMass = 0.493677;
	double ProtonMass = 0.938272013;
	
	const int numOfMultBins = 4;	
	const int numOfChTypes = 3;
	const int numOfkTbins = 4;

	int runmults[numOfMultBins] = {1,1,1,1};
	int multbins[numOfMultBins+1] = {0.001, 100, 200, 400, 600};
	
	int runch[numOfChTypes] = {1, 1, 1};
	const char *chrgs[numOfChTypes] = { "pip", "pim", "pippim" };
	
	int runktdep = 1;
	double ktrng[numOfkTbins+1] = {0.2, 0.4, 0.6, 0.8, 1.2};

	int gammacut = 0;	// cut na ee for gamma 
	
	double shqmax = 0.5; 
	int nbinssh = 200;

	//AliFemtoEventReaderESDChain *Reader = new AliFemtoEventReaderESDChain();
	//Reader->SetUseMultiplicity(AliFemtoEventReaderESDChain::kGlobalCount);
	
	AliFemtoEventReaderAODChain *Reader = new AliFemtoEventReaderAODChain();
	Reader->SetFilterBit(5);
	Reader->SetUseMultiplicity(AliFemtoEventReaderAOD::kCentrality);

	AliFemtoManager* Manager = new AliFemtoManager();
	Manager->SetEventReader(Reader);

	AliFemtoVertexMultAnalysis		*anetaphitpc[320];
	AliFemtoBasicEventCut			 *mecetaphitpc[320];
	AliFemtoCutMonitorEventMult	 *cutPassEvMetaphitpc[320];
	AliFemtoCutMonitorEventMult	 *cutFailEvMetaphitpc[320];
	AliFemtoCutMonitorEventVertex  *cutPassEvVetaphitpc[320];
	AliFemtoCutMonitorEventVertex  *cutFailEvVetaphitpc[320];

	AliFemtoESDTrackCut		*dtc5etaphitpc[320];
	AliFemtoESDTrackCut		*dtc6etaphitpc[320];

	AliFemtoTPCInnerCorrFctn      *PhiStarEtaetaphitpc[320];

	AliFemtoCutMonitorParticleYPt *cutPass6YPtetaphitpc[320];
	AliFemtoCutMonitorParticleYPt *cutFail6YPtetaphitpc[320];
	AliFemtoCutMonitorParticlePID *cutPass6PIDetaphitpc[320];
	AliFemtoCutMonitorParticlePID *cutFail6PIDetaphitpc[320];
	//	 AliFemtoShareQualityTPCEntranceSepPairCut	*sqpcetaphitpcsame[320];
	AliFemtoPairCutAntiGamma			*sqpcetaphitpc[320];
	//	AliFemtoPairCutRadialDistance			*sqpcetaphitpc[320];
	//	AliFemtoChi2CorrFctn				*cchiqinvetaphitpc[320];
	AliFemtoKTPairCut				 *ktpcuts[320];
	AliFemtoQinvCorrFctn				*cqinvkttpc[320];
	AliFemtoQinvCorrFctn				*cqinvtpc[320];
	AliFemtoCorrFctnDEtaDPhi			*cdedpetaphi[320];


	
	// *** Third QA task - HBT analysis ***
	// *** Begin pion-pion analysis ***
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

					//*** Event monitors ***
			
					cutPassEvMetaphitpc[aniter] = new AliFemtoCutMonitorEventMult(Form("cutPass%stpcM%i", chrgs[ichg], imult));
					cutFailEvMetaphitpc[aniter] = new AliFemtoCutMonitorEventMult(Form("cutFail%stpcM%i", chrgs[ichg], imult));
					mecetaphitpc[aniter]->AddCutMonitor(cutPassEvMetaphitpc[aniter], cutFailEvMetaphitpc[aniter]);
		
					cutPassEvVetaphitpc[aniter] = new AliFemtoCutMonitorEventVertex(Form("cutPass%stpcM%i", chrgs[ichg], imult));
					cutFailEvVetaphitpc[aniter] = new AliFemtoCutMonitorEventVertex(Form("cutFail%stpcM%i", chrgs[ichg], imult));
					mecetaphitpc[aniter]->AddCutMonitor(cutPassEvVetaphitpc[aniter], cutFailEvVetaphitpc[aniter]);

					//*** Cuts ***

					dtc5etaphitpc[aniter] = new AliFemtoESDTrackCut();
					dtc6etaphitpc[aniter] = new AliFemtoESDTrackCut();
  
					dtc5etaphitpc[aniter]->SetCharge(1.0);
					dtc6etaphitpc[aniter]->SetCharge(-1.0);

				
					dtc5etaphitpc[aniter]->SetEta(-1.0,1.0);
					dtc5etaphitpc[aniter]->SetPt(0.2,20);
					dtc5etaphitpc[aniter]->SetMass(PionMass);		
					dtc5etaphitpc[aniter]->SetMostProbablePion();

					dtc6etaphitpc[aniter]->SetEta(-1.0,1.0);
					dtc6etaphitpc[aniter]->SetPt(0.2,20);
					dtc6etaphitpc[aniter]->SetMass(PionMass);		
					dtc6etaphitpc[aniter]->SetMostProbablePion();

					// Track quality cuts					
					
					dtc5etaphitpc[aniter]->SetStatus(AliESDtrack::kTPCrefit|AliESDtrack::kITSrefit);
					dtc5etaphitpc[aniter]->SetminTPCncls(70);
					dtc5etaphitpc[aniter]->SetRemoveKinks(kTRUE);
					dtc5etaphitpc[aniter]->SetLabel(kFALSE);
					//	dtc5etaphitpc[aniter]->SetMaxITSChiNdof(6.0);
					dtc5etaphitpc[aniter]->SetMaxTPCChiNdof(4.0);
					dtc5etaphitpc[aniter]->SetMaxImpactXYPtDep(0.018, 0.035, -1.01) ;
					//	dtc5etaphitpc[aniter]->SetMaxImpactXYPtDep(0.0182, 0.0350, -1.01);
					dtc5etaphitpc[aniter]->SetMaxImpactZ(2);
					//	dtc5etaphitpc[aniter]->SetMaxSigmaToVertex(6.0);
					
                                             
					dtc6etaphitpc[aniter]->SetStatus(AliESDtrack::kTPCrefit|AliESDtrack::kITSrefit);
					dtc6etaphitpc[aniter]->SetminTPCncls(70);
					dtc6etaphitpc[aniter]->SetRemoveKinks(kTRUE);
					dtc6etaphitpc[aniter]->SetLabel(kFALSE);
					//      dtc6etaphitpc[aniter]->SetMaxITSChiNdof(6.0);
					dtc6etaphitpc[aniter]->SetMaxTPCChiNdof(4.0);
					dtc6etaphitpc[aniter]->SetMaxImpactXYPtDep(0.018, 0.035, -1.01) ;
					//      dtc6etaphitpc[aniter]->SetMaxImpactXYPtDep(0.0182, 0.0350, -1.01);
					dtc6etaphitpc[aniter]->SetMaxImpactZ(2);
					//      dtc6etaphitpc[aniter]->SetMaxSigmaToVertex(6.0);
                                             
					//*** Single track monitors ***
				
					cutPass6YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutPass%stpcM%i", chrgs[ichg], imult),PionMass);
					cutFail6YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutFail%stpcM%i", chrgs[ichg], imult),PionMass);
					dtc5etaphitpc[aniter]->AddCutMonitor(cutPass6YPtetaphitpc[aniter], cutFail6YPtetaphitpc[aniter]);				   

					cutPass6PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cutPass%stpcM%i", chrgs[ichg], imult),0);//0-pion,1-kaon,2-proton
					cutFail6PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cutFail%stpcM%i", chrgs[ichg], imult),0);
					dtc5etaphitpc[aniter]->AddCutMonitor(cutPass6PIDetaphitpc[aniter], cutFail6PIDetaphitpc[aniter]);


					//*** Two-track cuts monitors ***
					sqpcetaphitpc[aniter] = new AliFemtoPairCutAntiGamma();
		
					sqpcetaphitpc[aniter]->SetShareQualityMax(1.0);	// two track cuts on splitting and merging
					sqpcetaphitpc[aniter]->SetShareFractionMax(0.05);
					sqpcetaphitpc[aniter]->SetRemoveSameLabel(kFALSE);

					if (gammacut == 0)
					  {
					    sqpcetaphitpc[aniter]->SetMaxEEMinv(0.0);
					    sqpcetaphitpc[aniter]->SetMaxThetaDiff(0.0);
					  }
					else if (gammacut == 1)
					  { 
					    sqpcetaphitpc[aniter]->SetMaxEEMinv(0.002);
					    sqpcetaphitpc[aniter]->SetMaxThetaDiff(0.008);
					  }
				
		
					anetaphitpc[aniter]->SetEventCut(mecetaphitpc[aniter]);

					if(ichg==0)//pip
					{
						anetaphitpc[aniter]->SetFirstParticleCut(dtc5etaphitpc[aniter]);
						anetaphitpc[aniter]->SetSecondParticleCut(dtc5etaphitpc[aniter]);
					}
					if(ichg==1)//pim
					{
						anetaphitpc[aniter]->SetFirstParticleCut(dtc6etaphitpc[aniter]);
						anetaphitpc[aniter]->SetSecondParticleCut(dtc6etaphitpc[aniter]);
					
					}
					if(ichg==2)//pip pim
					{
						anetaphitpc[aniter]->SetFirstParticleCut(dtc5etaphitpc[aniter]);
						anetaphitpc[aniter]->SetSecondParticleCut(dtc6etaphitpc[aniter]);
					}

					anetaphitpc[aniter]->SetPairCut(sqpcetaphitpc[aniter]);

					// *** Correlation functions ***

					//Deta-Dphi correlation function
					cdedpetaphi[aniter] = new AliFemtoCorrFctnDEtaDPhi(Form("cdedp%stpcM%i", chrgs[ichg], imult),35, 35);
				
					//qinv correlation function
					cqinvtpc[aniter] = new AliFemtoQinvCorrFctn(Form("cqinv%stpcM%i", chrgs[ichg], imult),nbinssh,0.0,2.0);
					anetaphitpc[aniter]->AddCorrFctn(cdedpetaphi[aniter]);

					//Phi*-Eta monitor
					PhiStarEtaetaphitpc[aniter] = new AliFemtoTPCInnerCorrFctn(Form("PhistarEta%stpcM%i", chrgs[ichg], imult),nbinssh,0.0,2.0);
					anetaphitpc[aniter]->AddCorrFctn(PhiStarEtaetaphitpc[aniter]);

					if (runktdep)
					{
						int ktm;
						for (int ikt=0; ikt<numOfkTbins; ikt++)
						{
							ktm = aniter * numOfkTbins + ikt;
							ktpcuts[ktm] = new AliFemtoKTPairCut(ktrng[ikt], ktrng[ikt+1]);

							cqinvkttpc[ktm] = new AliFemtoQinvCorrFctn(Form("cqinv%stpcM%ikT%i", chrgs[ichg], imult, ikt),nbinssh,0.0,2.0);
							cqinvkttpc[ktm]->SetPairSelectionCut(ktpcuts[ktm]);
							anetaphitpc[aniter]->AddCorrFctn(cqinvkttpc[ktm]);

							PhiStarEtaetaphitpc[aniter] = new AliFemtoTPCInnerCorrFctn(Form("PhistarEta%stpcM%ikT%i", chrgs[ichg], imult, ikt),nbinssh,0.0,2.0);
							anetaphitpc[aniter]->AddCorrFctn(PhiStarEtaetaphitpc[aniter]);

						}
					}		
					Manager->AddAnalysis(anetaphitpc[aniter]);	
				}
			}
		}
	}
	// *** End pion-pion (positive) analysis

	return Manager;
}												 
