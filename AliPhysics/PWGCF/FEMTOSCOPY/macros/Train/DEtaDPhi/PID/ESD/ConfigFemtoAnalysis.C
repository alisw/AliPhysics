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
	
	const int numOfMultBins = 5;	
	const int numOfChTypes = 13;
	const int numOfkTbins = 2;

	int runmults[numOfMultBins] = {1, 1, 1, 0, 1};
	int multbins[numOfMultBins+1] = {2, 20, 50,150,2,150};
	
	int runch[numOfChTypes] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
	const char *chrgs[numOfChTypes] = { "PP", "aPaP", "PaP", "KpKp", "KmKm", "KpKm", "PIpPIp", "PImPIm", "PIpPIm", "all", "plus", "minus", "mixed" };
	
	int runktdep = 0;
	double ktrng[numOfkTbins+1] = {0.0, 0.7, 100.0};

	int runqinv = 1;
	int runshlcms = 1;// 0:PRF(PAP), 1:LCMS(PP,APAP)

	int runtype = 0; // Types 0 - global, 1 - ITS only, 2 - TPC Inner	//global tracks ->mfit ITS+TPC
	int isrealdata = 1;

	int gammacut = 1;	// cut na ee z gamma 
	
	double shqmax = 0.5; 
	int nbinssh = 100;

	AliFemtoEventReaderESDChain *Reader = new AliFemtoEventReaderESDChain();
	Reader->SetUseMultiplicity(AliFemtoEventReaderESDChain::kReferenceITSTPC);
	
	//AliFemtoEventReaderAODChain *Reader = new AliFemtoEventReaderAODChain();
	//Reader->SetFilterBit(7);
	// //Reader->SetCentralityPreSelection(0.001, 910);

	AliFemtoManager* Manager = new AliFemtoManager();
	Manager->SetEventReader(Reader);

	AliFemtoVertexMultAnalysis		*anetaphitpc[320];
	AliFemtoBasicEventCut				 *mecetaphitpc[320];
	AliFemtoCutMonitorEventMult	 *cutPassEvMetaphitpc[320];
	AliFemtoCutMonitorEventMult	 *cutFailEvMetaphitpc[320];
	AliFemtoCutMonitorEventVertex *cutPassEvVetaphitpc[320];
	AliFemtoCutMonitorEventVertex *cutFailEvVetaphitpc[320];
	AliFemtoESDTrackCut					 *dtc1etaphitpc[320];
	AliFemtoESDTrackCut					 *dtc2etaphitpc[320];
	AliFemtoESDTrackCut					 *dtc3etaphitpc[320];
	AliFemtoESDTrackCut					 *dtc4etaphitpc[320];
	AliFemtoESDTrackCut					 *dtc5etaphitpc[320];
	AliFemtoESDTrackCut					 *dtc6etaphitpc[320];
	AliFemtoESDTrackCut                                      *dtc7etaphitpc[320];
        AliFemtoESDTrackCut                                      *dtc8etaphitpc[320];
        AliFemtoESDTrackCut                                      *dtc9etaphitpc[320];
	AliFemtoCutMonitorParticleYPt *cutPass1YPtetaphitpc[320];
	AliFemtoCutMonitorParticleYPt *cutFail1YPtetaphitpc[320];
	AliFemtoCutMonitorParticlePID *cutPass1PIDetaphitpc[320];
	AliFemtoCutMonitorParticlePID *cutFail1PIDetaphitpc[320];
	AliFemtoCutMonitorParticleYPt *cutPass2YPtetaphitpc[320];
	AliFemtoCutMonitorParticleYPt *cutFail2YPtetaphitpc[320];
	AliFemtoCutMonitorParticlePID *cutPass2PIDetaphitpc[320];
	AliFemtoCutMonitorParticlePID *cutFail2PIDetaphitpc[320];
	AliFemtoCutMonitorParticleYPt *cutPass3YPtetaphitpc[320];
	AliFemtoCutMonitorParticleYPt *cutFail3YPtetaphitpc[320];
	AliFemtoCutMonitorParticlePID *cutPass3PIDetaphitpc[320];
	AliFemtoCutMonitorParticlePID *cutFail3PIDetaphitpc[320];
	AliFemtoCutMonitorParticleYPt *cutPass4YPtetaphitpc[320];
	AliFemtoCutMonitorParticleYPt *cutFail4YPtetaphitpc[320];
	AliFemtoCutMonitorParticlePID *cutPass4PIDetaphitpc[320];
	AliFemtoCutMonitorParticlePID *cutFail4PIDetaphitpc[320];
	AliFemtoCutMonitorParticleYPt *cutPass5YPtetaphitpc[320];
	AliFemtoCutMonitorParticleYPt *cutFail5YPtetaphitpc[320];
	AliFemtoCutMonitorParticlePID *cutPass5PIDetaphitpc[320];
	AliFemtoCutMonitorParticlePID *cutFail5PIDetaphitpc[320];
	AliFemtoCutMonitorParticleYPt *cutPass6YPtetaphitpc[320];
	AliFemtoCutMonitorParticleYPt *cutFail6YPtetaphitpc[320];
	AliFemtoCutMonitorParticlePID *cutPass6PIDetaphitpc[320];
	AliFemtoCutMonitorParticlePID *cutFail6PIDetaphitpc[320];
	//	 AliFemtoPairCutAntiGamma			*sqpcetaphitpcdiff[320];
	//	 AliFemtoShareQualityTPCEntranceSepPairCut			*sqpcetaphitpcsame[320];
	AliFemtoPairCutAntiGamma			*sqpcetaphitpc[320];
	//	AliFemtoPairCutRadialDistance			*sqpcetaphitpc[320];
	//	AliFemtoChi2CorrFctn					*cchiqinvetaphitpc[320];
	AliFemtoPairCutPt						 *ktpcuts[320];
	AliFemtoQinvCorrFctn					*cqinvkttpc[320];
	AliFemtoQinvCorrFctn					*cqinvtpc[320];
	AliFemtoCorrFctnDEtaDPhi			*cdedpetaphi[320];


	
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

					//if (isrealdata)mecetaphitpc[aniter]->SetAcceptOnlyPhysics(kTRUE);
				
					cutPassEvMetaphitpc[aniter] = new AliFemtoCutMonitorEventMult(Form("cutPass%stpcM%i", chrgs[ichg], imult));
					cutFailEvMetaphitpc[aniter] = new AliFemtoCutMonitorEventMult(Form("cutFail%stpcM%i", chrgs[ichg], imult));
					mecetaphitpc[aniter]->AddCutMonitor(cutPassEvMetaphitpc[aniter], cutFailEvMetaphitpc[aniter]);
		
					cutPassEvVetaphitpc[aniter] = new AliFemtoCutMonitorEventVertex(Form("cutPass%stpcM%i", chrgs[ichg], imult));
					cutFailEvVetaphitpc[aniter] = new AliFemtoCutMonitorEventVertex(Form("cutFail%stpcM%i", chrgs[ichg], imult));
					mecetaphitpc[aniter]->AddCutMonitor(cutPassEvVetaphitpc[aniter], cutFailEvVetaphitpc[aniter]);
		
					dtc1etaphitpc[aniter] = new AliFemtoESDTrackCut();
					dtc2etaphitpc[aniter] = new AliFemtoESDTrackCut();
					dtc3etaphitpc[aniter] = new AliFemtoESDTrackCut();
					dtc4etaphitpc[aniter] = new AliFemtoESDTrackCut();
					dtc5etaphitpc[aniter] = new AliFemtoESDTrackCut();
					dtc6etaphitpc[aniter] = new AliFemtoESDTrackCut();
                                        dtc7etaphitpc[aniter] = new AliFemtoESDTrackCut();
                                        dtc8etaphitpc[aniter] = new AliFemtoESDTrackCut();
                                        dtc9etaphitpc[aniter] = new AliFemtoESDTrackCut();

					dtc1etaphitpc[aniter]->SetCharge(1.0);
					dtc2etaphitpc[aniter]->SetCharge(-1.0);
					dtc3etaphitpc[aniter]->SetCharge(1.0);
					dtc4etaphitpc[aniter]->SetCharge(-1.0);
					dtc5etaphitpc[aniter]->SetCharge(1.0);
					dtc6etaphitpc[aniter]->SetCharge(-1.0);
					dtc8etaphitpc[aniter]->SetCharge(1.0);
					dtc9etaphitpc[aniter]->SetCharge(-1.0);
							

					dtc1etaphitpc[aniter]->SetEta(-1.0,1.0);
					dtc1etaphitpc[aniter]->SetPt(0.5,20);
					dtc1etaphitpc[aniter]->SetMass(ProtonMass);		
					dtc1etaphitpc[aniter]->SetMostProbableProton();
					dtc2etaphitpc[aniter]->SetEta(-1.0,1.0);
					dtc2etaphitpc[aniter]->SetPt(0.5,20);
					dtc2etaphitpc[aniter]->SetMass(ProtonMass);		
					dtc2etaphitpc[aniter]->SetMostProbableProton();
					dtc3etaphitpc[aniter]->SetEta(-1.0,1.0);
					dtc3etaphitpc[aniter]->SetPt(0.3,20);
					dtc3etaphitpc[aniter]->SetMass(KaonMass);
					dtc3etaphitpc[aniter]->SetMostProbableKaon();
					dtc4etaphitpc[aniter]->SetEta(-1.0,1.0);
					dtc4etaphitpc[aniter]->SetPt(0.3,20);
					dtc4etaphitpc[aniter]->SetMass(KaonMass);		
					dtc4etaphitpc[aniter]->SetMostProbableKaon();
					dtc5etaphitpc[aniter]->SetEta(-1.0,1.0);
					dtc5etaphitpc[aniter]->SetPt(0.2,20);
					dtc5etaphitpc[aniter]->SetMass(PionMass);		
					dtc5etaphitpc[aniter]->SetMostProbablePion();
					dtc6etaphitpc[aniter]->SetEta(-1.0,1.0);
					dtc6etaphitpc[aniter]->SetPt(0.2,20);
					dtc6etaphitpc[aniter]->SetMass(PionMass);		
					dtc6etaphitpc[aniter]->SetMostProbablePion();
                                        dtc7etaphitpc[aniter]->SetEta(-1.0,1.0);
                                        dtc7etaphitpc[aniter]->SetPt(0.2,20);
                                        dtc8etaphitpc[aniter]->SetEta(-1.0,1.0);
                                        dtc8etaphitpc[aniter]->SetPt(0.2,20);
                                        dtc9etaphitpc[aniter]->SetEta(-1.0,1.0);
                                        dtc9etaphitpc[aniter]->SetPt(0.2,20);

				// Track quality cuts

					if (runtype == 0)
					{
						if(ichg==0 || ichg==2)
						{
							dtc1etaphitpc[aniter]->SetStatus(AliESDtrack::kTPCrefit|AliESDtrack::kITSrefit);
							dtc1etaphitpc[aniter]->SetminTPCncls(70);
							dtc1etaphitpc[aniter]->SetRemoveKinks(kTRUE);
							dtc1etaphitpc[aniter]->SetLabel(kFALSE);
							//	dtc1etaphitpc[aniter]->SetMaxITSChiNdof(6.0);
							dtc1etaphitpc[aniter]->SetMaxTPCChiNdof(4.0);	// pisac
							dtc1etaphitpc[aniter]->SetMaxImpactXYPtDep(0.018, 0.035, -1.01); 	//	DCA xy
							//	dtc1etaphitpc[aniter]->SetMaxImpactXYPtDep(0.0182, 0.0350, -1.01);
							dtc1etaphitpc[aniter]->SetMaxImpactZ(2);	//DCA Z
							//	dtc1etaphitpc[aniter]->SetMaxSigmaToVertex(6.0);
						}
						if(ichg==1 || ichg==2)
						{
							dtc2etaphitpc[aniter]->SetStatus(AliESDtrack::kTPCrefit|AliESDtrack::kITSrefit);
							dtc2etaphitpc[aniter]->SetminTPCncls(70);
							dtc2etaphitpc[aniter]->SetRemoveKinks(kTRUE);
							dtc2etaphitpc[aniter]->SetLabel(kFALSE);
							//	dtc2etaphitpc[aniter]->SetMaxITSChiNdof(6.0);
							dtc2etaphitpc[aniter]->SetMaxTPCChiNdof(4.0);
							dtc2etaphitpc[aniter]->SetMaxImpactXYPtDep(0.018, 0.035, -1.01) ;
							//	dtc2etaphitpc[aniter]->SetMaxImpactXYPtDep(0.0182, 0.0350, -1.01);
							dtc2etaphitpc[aniter]->SetMaxImpactZ(2);
							//	dtc2etaphitpc[aniter]->SetMaxSigmaToVertex(6.0);
						}
						if(ichg==3 || ichg==5)
						{
							dtc3etaphitpc[aniter]->SetStatus(AliESDtrack::kTPCrefit|AliESDtrack::kITSrefit);
							dtc3etaphitpc[aniter]->SetminTPCncls(70);
							dtc3etaphitpc[aniter]->SetRemoveKinks(kTRUE);
							dtc3etaphitpc[aniter]->SetLabel(kFALSE);
							//	dtc3etaphitpc[aniter]->SetMaxITSChiNdof(6.0);
							dtc3etaphitpc[aniter]->SetMaxTPCChiNdof(4.0);
							dtc3etaphitpc[aniter]->SetMaxImpactXYPtDep(0.018, 0.035, -1.01) ;
							//	dtc3etaphitpc[aniter]->SetMaxImpactXYPtDep(0.0182, 0.0350, -1.01);
							dtc3etaphitpc[aniter]->SetMaxImpactZ(2);
							//	dtc3etaphitpc[aniter]->SetMaxSigmaToVertex(6.0);
						}
						if(ichg==4 || ichg==5)
						{
							dtc4etaphitpc[aniter]->SetStatus(AliESDtrack::kTPCrefit|AliESDtrack::kITSrefit);
							dtc4etaphitpc[aniter]->SetminTPCncls(70);
							dtc4etaphitpc[aniter]->SetRemoveKinks(kTRUE);
							dtc4etaphitpc[aniter]->SetLabel(kFALSE);
							//	dtc4etaphitpc[aniter]->SetMaxITSChiNdof(6.0);
							dtc4etaphitpc[aniter]->SetMaxTPCChiNdof(4.0);
							dtc4etaphitpc[aniter]->SetMaxImpactXYPtDep(0.018, 0.035, -1.01) ;
							//	dtc4etaphitpc[aniter]->SetMaxImpactXYPtDep(0.0182, 0.0350, -1.01);
							dtc4etaphitpc[aniter]->SetMaxImpactZ(2);
							//	dtc4etaphitpc[aniter]->SetMaxSigmaToVertex(6.0);
						}
						if(ichg==6 || ichg==8)
						{
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
						}
                                                if(ichg==7 || ichg ==8)
                                                {
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
                                                }

						if(ichg==9)
						{
							dtc7etaphitpc[aniter]->SetStatus(AliESDtrack::kTPCrefit|AliESDtrack::kITSrefit);
							dtc7etaphitpc[aniter]->SetminTPCncls(70);
							dtc7etaphitpc[aniter]->SetRemoveKinks(kTRUE);
							dtc7etaphitpc[aniter]->SetLabel(kFALSE);
							//	dtc6etaphitpc[aniter]->SetMaxITSChiNdof(6.0);
							dtc7etaphitpc[aniter]->SetMaxTPCChiNdof(4.0);
							dtc7etaphitpc[aniter]->SetMaxImpactXYPtDep(0.018, 0.035, -1.01) ;
							//	dtc6etaphitpc[aniter]->SetMaxImpactXYPtDep(0.0182, 0.0350, -1.01);
							dtc7etaphitpc[aniter]->SetMaxImpactZ(2);
							//	dtc6etaphitpc[aniter]->SetMaxSigmaToVertex(6.0);
						}
                                                if(ichg==10 || ichg==12)
                                                {
                                                        dtc8etaphitpc[aniter]->SetStatus(AliESDtrack::kTPCrefit|AliESDtrack::kITSrefit);
                                                        dtc8etaphitpc[aniter]->SetminTPCncls(70);
                                                        dtc8etaphitpc[aniter]->SetRemoveKinks(kTRUE);
                                                        dtc8etaphitpc[aniter]->SetLabel(kFALSE);
                                                        //      dtc6etaphitpc[aniter]->SetMaxITSChiNdof(6.0);
                                                        dtc8etaphitpc[aniter]->SetMaxTPCChiNdof(4.0);
                                                        dtc8etaphitpc[aniter]->SetMaxImpactXYPtDep(0.018, 0.035, -1.01) ;
                                                        //      dtc6etaphitpc[aniter]->SetMaxImpactXYPtDep(0.0182, 0.0350, -1.01);
                                                        dtc8etaphitpc[aniter]->SetMaxImpactZ(2);
                                                        //      dtc6etaphitpc[aniter]->SetMaxSigmaToVertex(6.0);
                                                }

                                                if(ichg==11 || ichg==12)
                                                {
                                                        dtc9etaphitpc[aniter]->SetStatus(AliESDtrack::kTPCrefit|AliESDtrack::kITSrefit);
                                                        dtc9etaphitpc[aniter]->SetminTPCncls(70);
                                                        dtc9etaphitpc[aniter]->SetRemoveKinks(kTRUE);
                                                        dtc9etaphitpc[aniter]->SetLabel(kFALSE);
                                                        //      dtc6etaphitpc[aniter]->SetMaxITSChiNdof(6.0);
                                                        dtc9etaphitpc[aniter]->SetMaxTPCChiNdof(4.0);
                                                        dtc9etaphitpc[aniter]->SetMaxImpactXYPtDep(0.018, 0.035, -1.01) ;
                                                        //      dtc6etaphitpc[aniter]->SetMaxImpactXYPtDep(0.0182, 0.0350, -1.01);
                                                        dtc9etaphitpc[aniter]->SetMaxImpactZ(2);
                                                        //      dtc6etaphitpc[aniter]->SetMaxSigmaToVertex(6.0);
                                                }

					}
					else if (runtype == 1)
					{
						;
					}
					else if (runtype == 2)
					{
						;
					}
					if(ichg==0 || ichg==2)
					{
						cutPass1YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutPass1%stpcM%i", chrgs[ichg], imult),ProtonMass);
						cutFail1YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutFail1%stpcM%i", chrgs[ichg], imult),ProtonMass);
						dtc1etaphitpc[aniter]->AddCutMonitor(cutPass1YPtetaphitpc[aniter], cutFail1YPtetaphitpc[aniter]);
						cutPass1PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cutPass1%stpcM%i", chrgs[ichg], imult),2);//0-pion,1-kaon,2-proton
						cutFail1PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cutFail1%stpcM%i", chrgs[ichg], imult),2);
						dtc1etaphitpc[aniter]->AddCutMonitor(cutPass1PIDetaphitpc[aniter], cutFail1PIDetaphitpc[aniter]);
					}
					if(ichg==1 || ichg==2)
					{
						cutPass2YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutPass2%stpcM%i", chrgs[ichg], imult),ProtonMass);
						cutFail2YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutFail2%stpcM%i", chrgs[ichg], imult),ProtonMass);
						dtc2etaphitpc[aniter]->AddCutMonitor(cutPass2YPtetaphitpc[aniter], cutFail2YPtetaphitpc[aniter]);

						cutPass2PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cutPass2%stpcM%i", chrgs[ichg], imult),2);//0-pion,1-kaon,2-proton
						cutFail2PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cutFail2%stpcM%i", chrgs[ichg], imult),2);
						dtc2etaphitpc[aniter]->AddCutMonitor(cutPass2PIDetaphitpc[aniter], cutFail2PIDetaphitpc[aniter]);
					}

					if(ichg==3 || ichg==5)
					{
						cutPass3YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutPass3%stpcM%i", chrgs[ichg], imult),KaonMass);
						cutFail3YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutFail3%stpcM%i", chrgs[ichg], imult),KaonMass);
						dtc3etaphitpc[aniter]->AddCutMonitor(cutPass3YPtetaphitpc[aniter], cutFail3YPtetaphitpc[aniter]);
						cutPass3PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cutPass3%stpcM%i", chrgs[ichg], imult),1);//0-pion,1-kaon,2-proton
						cutFail3PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cutFail3%stpcM%i", chrgs[ichg], imult),1);
						dtc3etaphitpc[aniter]->AddCutMonitor(cutPass3PIDetaphitpc[aniter], cutFail3PIDetaphitpc[aniter]);
					}
					if(ichg==4 || ichg==5)
					{
						cutPass4YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutPass4%stpcM%i", chrgs[ichg], imult),KaonMass);
						cutFail4YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutFail4%stpcM%i", chrgs[ichg], imult),KaonMass);
						dtc4etaphitpc[aniter]->AddCutMonitor(cutPass4YPtetaphitpc[aniter], cutFail4YPtetaphitpc[aniter]);

						cutPass4PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cutPass4%stpcM%i", chrgs[ichg], imult),1);//0-pion,1-kaon,2-proton
						cutFail4PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cutFail4%stpcM%i", chrgs[ichg], imult),1);
						dtc4etaphitpc[aniter]->AddCutMonitor(cutPass4PIDetaphitpc[aniter], cutFail4PIDetaphitpc[aniter]);
					}

					if(ichg==6 || ichg==8)
					{
						cutPass5YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutPass5%stpcM%i", chrgs[ichg], imult),PionMass);
						cutFail5YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutFail5%stpcM%i", chrgs[ichg], imult),PionMass);
						dtc5etaphitpc[aniter]->AddCutMonitor(cutPass5YPtetaphitpc[aniter], cutFail5YPtetaphitpc[aniter]);
						cutPass5PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cutPass5%stpcM%i", chrgs[ichg], imult),0);//0-pion,1-kaon,2-proton
						cutFail5PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cutFail5%stpcM%i", chrgs[ichg], imult),0);
						dtc5etaphitpc[aniter]->AddCutMonitor(cutPass5PIDetaphitpc[aniter], cutFail5PIDetaphitpc[aniter]);
					}
					if(ichg==7 || ichg==8)
					{
						cutPass6YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutPass6%stpcM%i", chrgs[ichg], imult),PionMass);
						cutFail6YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutFail6%stpcM%i", chrgs[ichg], imult),PionMass);
						dtc6etaphitpc[aniter]->AddCutMonitor(cutPass6YPtetaphitpc[aniter], cutFail6YPtetaphitpc[aniter]);

						cutPass6PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cutPass6%stpcM%i", chrgs[ichg], imult),0);//0-pion,1-kaon,2-proton
						cutFail6PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cutFail6%stpcM%i", chrgs[ichg], imult),0);
						dtc6etaphitpc[aniter]->AddCutMonitor(cutPass6PIDetaphitpc[aniter], cutFail6PIDetaphitpc[aniter]);
					}
					if(ichg>8)
					  {
					    cutPass6YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutPass%stpcM%i", chrgs[ichg], imult),PionMass);
					    cutFail6YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutFail%stpcM%i", chrgs[ichg], imult),PionMass);
					    if(ichg==9) dtc7etaphitpc[aniter]->AddCutMonitor(cutPass6YPtetaphitpc[aniter], cutFail6YPtetaphitpc[aniter]);
					    if(ichg==10 || ichg==12) dtc8etaphitpc[aniter]->AddCutMonitor(cutPass6YPtetaphitpc[aniter], cutFail6YPtetaphitpc[aniter]);
					    if(ichg==11 || ichg==12) dtc9etaphitpc[aniter]->AddCutMonitor(cutPass6YPtetaphitpc[aniter], cutFail6YPtetaphitpc[aniter]);

					    cutPass6PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cutPass%stpcM%i", chrgs[ichg], imult),0);//0-pion,1-kaon,2-proton
					    cutFail6PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cutFail%stpcM%i", chrgs[ichg], imult),0);
		    			    if(ichg==9) dtc7etaphitpc[aniter]->AddCutMonitor(cutPass6PIDetaphitpc[aniter], cutFail6PIDetaphitpc[aniter]);
					    if(ichg==10 || ichg==12) dtc8etaphitpc[aniter]->AddCutMonitor(cutPass6PIDetaphitpc[aniter], cutFail6PIDetaphitpc[aniter]);
					    if(ichg==11 || ichg==12) dtc9etaphitpc[aniter]->AddCutMonitor(cutPass6PIDetaphitpc[aniter], cutFail6PIDetaphitpc[aniter]);
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
						sqpcetaphitpc[aniter]->SetShareQualityMax(1.0);
						sqpcetaphitpc[aniter]->SetShareFractionMax(1.05);
						sqpcetaphitpc[aniter]->SetRemoveSameLabel(kFALSE);
						// sqpcetaphitpc[aniter]->SetMaxEEMinv(0.002);
						// sqpcetaphitpc[aniter]->SetMaxThetaDiff(0.008);
						//	sqpcetaphitpc[aniter]->SetTPCEntranceSepMinimum(5.0);
						// sqpcetaphitpc[aniter]->SetRadialDistanceMinimum(1.2, 0.03);
						//	sqpcetaphitpc[aniter]->SetEtaDifferenceMinimum(0.02);
					}
					else if (runtype == 2)
					{
						sqpcetaphitpc[aniter]->SetDataType(kESD);
						sqpcetaphitpc[aniter]->SetShareQualityMax(1.0);
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
					}
		
					anetaphitpc[aniter]->SetEventCut(mecetaphitpc[aniter]);

					if(ichg==0)
					{
						anetaphitpc[aniter]->SetFirstParticleCut(dtc1etaphitpc[aniter]);
						anetaphitpc[aniter]->SetSecondParticleCut(dtc1etaphitpc[aniter]);
					}
					if(ichg==1)
					{
						anetaphitpc[aniter]->SetFirstParticleCut(dtc2etaphitpc[aniter]);
						anetaphitpc[aniter]->SetSecondParticleCut(dtc2etaphitpc[aniter]);
					}
					if(ichg==2)
					{
						anetaphitpc[aniter]->SetFirstParticleCut(dtc1etaphitpc[aniter]);
						anetaphitpc[aniter]->SetSecondParticleCut(dtc2etaphitpc[aniter]);
					}
					if(ichg==3)
					{
						anetaphitpc[aniter]->SetFirstParticleCut(dtc3etaphitpc[aniter]);
						anetaphitpc[aniter]->SetSecondParticleCut(dtc3etaphitpc[aniter]);
					}
					if(ichg==4)
					{
						anetaphitpc[aniter]->SetFirstParticleCut(dtc4etaphitpc[aniter]);
						anetaphitpc[aniter]->SetSecondParticleCut(dtc4etaphitpc[aniter]);
					
					}
					if(ichg==5)
					{
						anetaphitpc[aniter]->SetFirstParticleCut(dtc3etaphitpc[aniter]);
						anetaphitpc[aniter]->SetSecondParticleCut(dtc4etaphitpc[aniter]);
					}
					if(ichg==6)
					{
						anetaphitpc[aniter]->SetFirstParticleCut(dtc5etaphitpc[aniter]);
						anetaphitpc[aniter]->SetSecondParticleCut(dtc5etaphitpc[aniter]);
					}
					if(ichg==7)
					{
						anetaphitpc[aniter]->SetFirstParticleCut(dtc6etaphitpc[aniter]);
						anetaphitpc[aniter]->SetSecondParticleCut(dtc6etaphitpc[aniter]);
					
					}
					if(ichg==8)
					{
						anetaphitpc[aniter]->SetFirstParticleCut(dtc5etaphitpc[aniter]);
						anetaphitpc[aniter]->SetSecondParticleCut(dtc6etaphitpc[aniter]);
					}

                                        if(ichg==9)
                                        {
					  anetaphitpc[aniter]->SetFirstParticleCut(dtc7etaphitpc[aniter]);
					  anetaphitpc[aniter]->SetSecondParticleCut(dtc7etaphitpc[aniter]);
                                        }

                                        if(ichg==10)
                                        {
                                                anetaphitpc[aniter]->SetFirstParticleCut(dtc8etaphitpc[aniter]);
						anetaphitpc[aniter]->SetSecondParticleCut(dtc8etaphitpc[aniter]);
                                        }

                                        if(ichg==11)
                                        {
                                                anetaphitpc[aniter]->SetFirstParticleCut(dtc9etaphitpc[aniter]);
						anetaphitpc[aniter]->SetSecondParticleCut(dtc9etaphitpc[aniter]);
                                         
                                        }
                                        if(ichg==12)
                                        {
                                                anetaphitpc[aniter]->SetFirstParticleCut(dtc8etaphitpc[aniter]);
                                                anetaphitpc[aniter]->SetSecondParticleCut(dtc9etaphitpc[aniter]);
                                        }



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
	// *** End pion-pion (positive) analysis

	return Manager;
}												 
