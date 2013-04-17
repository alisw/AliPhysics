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
#include "AliFemtoSimpleAnalysis.h"
#include "AliFemtoBasicEventCut.h"
#include "AliFemtoESDTrackCut.h"
#include "AliFemtoCorrFctn.h"
#include "AliFemtoShareQualityTPCEntranceSepPairCut.h"
#include "AliFemtoPairCutAntiGamma.h"
#include "AliFemtoPairCutRadialDistance.h"
#include "AliFemtoShareQualityCorrFctn.h"
#include "AliFemtoTPCInnerCorrFctn.h"
#include "AliFemtoVertexMultAnalysis.h"
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
	int multbins[numOfMultBins+1] = {0, 20, 50, 150, 2, 150};
	
	int runch[numOfChTypes] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
	const char *chrgs[numOfChTypes] = { "PP", "aPaP", "PaP", "KpKp", "KmKm", "KpKm", "PIpPIp", "PImPIm", "PIpPIm", "all", "plus", "minus", "mixed" };
	
	int runktdep = 0;
	double ktrng[numOfkTbins+1] = {0.0, 0.7, 100.0};

	int runqinv = 1;
	int runshlcms = 1;// 0:PRF(PAP), 1:LCMS(PP,APAP)

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
	AliFemtoCutMonitorEventVertex *cutFailEvVetaphitpc[numOfMultBins*numOfChTypes];*/
	AliFemtoESDTrackCut					 *dtc1etaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoESDTrackCut					 *dtc2etaphitpc[numOfMultBins*numOfChTypes];
/*	AliFemtoCutMonitorParticleYPt *cutPass1YPtetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoCutMonitorParticleYPt *cutFail1YPtetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoCutMonitorParticlePID *cutPass1PIDetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoCutMonitorParticlePID *cutFail1PIDetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoCutMonitorParticleYPt *cutPass2YPtetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoCutMonitorParticleYPt *cutFail2YPtetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoCutMonitorParticlePID *cutPass2PIDetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoCutMonitorParticlePID *cutFail2PIDetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoCutMonitorParticleYPt *cutPass3YPtetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoCutMonitorParticleYPt *cutFail3YPtetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoCutMonitorParticlePID *cutPass3PIDetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoCutMonitorParticlePID *cutFail3PIDetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoCutMonitorParticleYPt *cutPass4YPtetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoCutMonitorParticleYPt *cutFail4YPtetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoCutMonitorParticlePID *cutPass4PIDetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoCutMonitorParticlePID *cutFail4PIDetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoCutMonitorParticleYPt *cutPass5YPtetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoCutMonitorParticleYPt *cutFail5YPtetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoCutMonitorParticlePID *cutPass5PIDetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoCutMonitorParticlePID *cutFail5PIDetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoCutMonitorParticleYPt *cutPass6YPtetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoCutMonitorParticleYPt *cutFail6YPtetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoCutMonitorParticlePID *cutPass6PIDetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoCutMonitorParticlePID *cutFail6PIDetaphitpc[numOfMultBins*numOfChTypes];*/
	//AliFemtoCutMonitorParticlePtPDG *cutPassPtPDGetaphitpc[numOfMultBins*numOfChTypes];
	//AliFemtoCutMonitorParticlePtPDG *cutFailPtPDGetaphitpc[numOfMultBins*numOfChTypes];

	//	 AliFemtoShareQualityTPCEntranceSepPairCut			*sqpcetaphitpcsame[numOfMultBins*numOfChTypes];
	AliFemtoPairCutAntiGamma			*sqpcetaphitpc[numOfMultBins*numOfChTypes];
	//	AliFemtoPairCutRadialDistance			*sqpcetaphitpc[numOfMultBins*numOfChTypes];
	//	AliFemtoChi2CorrFctn					*cchiqinvetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoPairCutPt						 *ktpcuts[numOfMultBins*numOfChTypes];
//	AliFemtoQinvCorrFctn					*cqinvkttpc[numOfMultBins*numOfChTypes];
//	AliFemtoQinvCorrFctn					*cqinvtpc[numOfMultBins*numOfChTypes];
	AliFemtoCorrFctnDEtaDPhi			*cdedpetaphi[numOfMultBins*numOfChTypes];
	//AliFemtoCorrFctnMassInvMonitor			*cMinvMonitor[numOfMultBins*numOfChTypes];

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
				
					/*cutPassEvMetaphitpc[aniter] = new AliFemtoCutMonitorEventMult(Form("cutPass%stpcM%i", chrgs[ichg], imult));
					cutFailEvMetaphitpc[aniter] = new AliFemtoCutMonitorEventMult(Form("cutFail%stpcM%i", chrgs[ichg], imult));
					mecetaphitpc[aniter]->AddCutMonitor(cutPassEvMetaphitpc[aniter], cutFailEvMetaphitpc[aniter]);
		
					cutPassEvVetaphitpc[aniter] = new AliFemtoCutMonitorEventVertex(Form("cutPass%stpcM%i", chrgs[ichg], imult));
					cutFailEvVetaphitpc[aniter] = new AliFemtoCutMonitorEventVertex(Form("cutFail%stpcM%i", chrgs[ichg], imult));
					mecetaphitpc[aniter]->AddCutMonitor(cutPassEvVetaphitpc[aniter], cutFailEvVetaphitpc[aniter]);
		*/
					dtc1etaphitpc[aniter] = new AliFemtoESDTrackCut();
					dtc2etaphitpc[aniter] = new AliFemtoESDTrackCut();

					dtc1etaphitpc[aniter]->SetCharge(1.0);
					dtc2etaphitpc[aniter]->SetCharge(-1.0);

					dtc1etaphitpc[aniter]->SetEta(-1.0,1.0);

					if( ichg == 0 || ichg == 1 || ichg == 2)
					{
						dtc1etaphitpc[aniter]->SetPt(0.5,20);
						dtc1etaphitpc[aniter]->SetMass(ProtonMass);		
						dtc1etaphitpc[aniter]->SetMostProbableProton();
						dtc2etaphitpc[aniter]->SetPt(0.5,20);
						dtc2etaphitpc[aniter]->SetMass(ProtonMass);		
						dtc2etaphitpc[aniter]->SetMostProbableProton();
					}
					if( ichg == 3 || ichg == 4 || ichg == 5)
					{
						dtc1etaphitpc[aniter]->SetPt(0.3,20);
						dtc1etaphitpc[aniter]->SetMass(KaonMass);
						dtc1etaphitpc[aniter]->SetMostProbableKaon();
						dtc2etaphitpc[aniter]->SetPt(0.3,20);
						dtc2etaphitpc[aniter]->SetMass(KaonMass);
						dtc2etaphitpc[aniter]->SetMostProbableKaon();
					}
					if( ichg == 6 || ichg == 7 || ichg == 8)
					{
						dtc5etaphitpc[aniter]->SetPt(0.2,20);
						dtc5etaphitpc[aniter]->SetMass(PionMass);		
						dtc5etaphitpc[aniter]->SetMostProbablePion();
						dtc6etaphitpc[aniter]->SetPt(0.2,20);
						dtc6etaphitpc[aniter]->SetMass(PionMass);		
						dtc6etaphitpc[aniter]->SetMostProbablePion();
		         }
					if( ichg == 9 )
					{
						dtc3etaphitpc[aniter] = new AliFemtoESDTrackCut();
						dtc3etaphitpc[aniter]->SetEta(-1.0,1.0);
	            	dtc3etaphitpc[aniter]->SetPt(0.2,20);
		         }
					if( ichg == 10 || ichg == 11 || ichg == 12)
					{
		         	dtc1etaphitpc[aniter]->SetPt(0.2,20);
		         	dtc2etaphitpc[aniter]->SetPt(0.2,20);
					}
				// Track quality cuts

					if (runtype == 0)
					{
				//		if(ichg==0 || ichg==2)
						{
							dtc1etaphitpc[aniter]->SetStatus(AliESDtrack::kTPCrefit|AliESDtrack::kITSrefit);
							dtc1etaphitpc[aniter]->SetminTPCncls(70);
							dtc1etaphitpc[aniter]->SetRemoveKinks(kTRUE);
							dtc1etaphitpc[aniter]->SetLabel(kFALSE);
							//	dtc1etaphitpc[aniter]->SetMaxITSChiNdof(6.0);
							dtc1etaphitpc[aniter]->SetMaxTPCChiNdof(4.0);
							dtc1etaphitpc[aniter]->SetMaxImpactXYPtDep(0.018, 0.035, -1.01); 	//	DCA xy
							//	dtc1etaphitpc[aniter]->SetMaxImpactXYPtDep(0.0182, 0.0350, -1.01);
							dtc1etaphitpc[aniter]->SetMaxImpactZ(2);	//DCA Z
							//	dtc1etaphitpc[aniter]->SetMaxSigmaToVertex(6.0);
						}
			//			if(ichg==1 || ichg==2)
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
						if(ichg==9)
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
					}
					else if (runtype == 1)
					{
						;
					}
					else if (runtype == 2)
					{
						;
					}
/*
					cutPassPtPDGetaphitpc[aniter] = new AliFemtoCutMonitorParticlePtPDG(Form("cutPass%stpcM%i", chrgs[ichg], imult),PionMass);
					cutFailPtPDGetaphitpc[aniter] = new AliFemtoCutMonitorParticlePtPDG(Form("cutFail%stpcM%i", chrgs[ichg], imult),PionMass);
					if(ichg==0 || ichg==2)
					  {
					    dtc1etaphitpc[aniter]->AddCutMonitor(cutPassPtPDGetaphitpc[aniter],cutFailPtPDGetaphitpc[aniter]);
					    cutPass1YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutPass1%stpcM%i", chrgs[ichg], imult),ProtonMass);
					    cutFail1YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutFail1%stpcM%i", chrgs[ichg], imult),ProtonMass);
					    dtc1etaphitpc[aniter]->AddCutMonitor(cutPass1YPtetaphitpc[aniter], cutFail1YPtetaphitpc[aniter]);
					    cutPass1PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cutPass1%stpcM%i", chrgs[ichg], imult),2);//0-pion,1-kaon,2-proton
					    cutFail1PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cutFail1%stpcM%i", chrgs[ichg], imult),2);
					    dtc1etaphitpc[aniter]->AddCutMonitor(cutPass1PIDetaphitpc[aniter], cutFail1PIDetaphitpc[aniter]);
					  }
					else if(ichg==1)
					  {
					    dtc2etaphitpc[aniter]->AddCutMonitor(cutPassPtPDGetaphitpc[aniter],cutFailPtPDGetaphitpc[aniter]);
					    cutPass2YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutPass2%stpcM%i", chrgs[ichg], imult),ProtonMass);
					    cutFail2YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutFail2%stpcM%i", chrgs[ichg], imult),ProtonMass);
					    dtc2etaphitpc[aniter]->AddCutMonitor(cutPass2YPtetaphitpc[aniter], cutFail2YPtetaphitpc[aniter]);

					    cutPass2PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cutPass2%stpcM%i", chrgs[ichg], imult),2);//0-pion,1-kaon,2-proton
					    cutFail2PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cutFail2%stpcM%i", chrgs[ichg], imult),2);
					    dtc2etaphitpc[aniter]->AddCutMonitor(cutPass2PIDetaphitpc[aniter], cutFail2PIDetaphitpc[aniter]);
					  }

					else if(ichg==3 || ichg==5)
					  {
					    dtc3etaphitpc[aniter]->AddCutMonitor(cutPassPtPDGetaphitpc[aniter],cutFailPtPDGetaphitpc[aniter]);
					    cutPass3YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutPass3%stpcM%i", chrgs[ichg], imult),KaonMass);
					    cutFail3YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutFail3%stpcM%i", chrgs[ichg], imult),KaonMass);
					    dtc3etaphitpc[aniter]->AddCutMonitor(cutPass3YPtetaphitpc[aniter], cutFail3YPtetaphitpc[aniter]);
					    cutPass3PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cutPass3%stpcM%i", chrgs[ichg], imult),1);//0-pion,1-kaon,2-proton
					    cutFail3PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cutFail3%stpcM%i", chrgs[ichg], imult),1);
					    dtc3etaphitpc[aniter]->AddCutMonitor(cutPass3PIDetaphitpc[aniter], cutFail3PIDetaphitpc[aniter]);
					  }
					else if(ichg==4)
					  {
					    dtc4etaphitpc[aniter]->AddCutMonitor(cutPassPtPDGetaphitpc[aniter],cutFailPtPDGetaphitpc[aniter]);
					    cutPass4YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutPass4%stpcM%i", chrgs[ichg], imult),KaonMass);
					    cutFail4YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutFail4%stpcM%i", chrgs[ichg], imult),KaonMass);
					    dtc4etaphitpc[aniter]->AddCutMonitor(cutPass4YPtetaphitpc[aniter], cutFail4YPtetaphitpc[aniter]);

					    cutPass4PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cutPass4%stpcM%i", chrgs[ichg], imult),1);//0-pion,1-kaon,2-proton
					    cutFail4PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cutFail4%stpcM%i", chrgs[ichg], imult),1);
					    dtc4etaphitpc[aniter]->AddCutMonitor(cutPass4PIDetaphitpc[aniter], cutFail4PIDetaphitpc[aniter]);
					  }

					else if(ichg==6 || ichg==8)
					  {
					    dtc5etaphitpc[aniter]->AddCutMonitor(cutPassPtPDGetaphitpc[aniter],cutFailPtPDGetaphitpc[aniter]);
					    cutPass5YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutPass5%stpcM%i", chrgs[ichg], imult),PionMass);
					    cutFail5YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutFail5%stpcM%i", chrgs[ichg], imult),PionMass);
					    dtc5etaphitpc[aniter]->AddCutMonitor(cutPass5YPtetaphitpc[aniter], cutFail5YPtetaphitpc[aniter]);
					    cutPass5PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cutPass5%stpcM%i", chrgs[ichg], imult),0);//0-pion,1-kaon,2-proton
					    cutFail5PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cutFail5%stpcM%i", chrgs[ichg], imult),0);
					    dtc5etaphitpc[aniter]->AddCutMonitor(cutPass5PIDetaphitpc[aniter], cutFail5PIDetaphitpc[aniter]);
					  }
					else if(ichg==7)
					  {
					    dtc6etaphitpc[aniter]->AddCutMonitor(cutPassPtPDGetaphitpc[aniter],cutFailPtPDGetaphitpc[aniter]);
					    cutPass6YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutPass6%stpcM%i", chrgs[ichg], imult),PionMass);
					    cutFail6YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutFail6%stpcM%i", chrgs[ichg], imult),PionMass);
					    dtc6etaphitpc[aniter]->AddCutMonitor(cutPass6YPtetaphitpc[aniter], cutFail6YPtetaphitpc[aniter]);

					    cutPass6PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cutPass6%stpcM%i", chrgs[ichg], imult),0);//0-pion,1-kaon,2-proton
					    cutFail6PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cutFail6%stpcM%i", chrgs[ichg], imult),0);
					    dtc6etaphitpc[aniter]->AddCutMonitor(cutPass6PIDetaphitpc[aniter], cutFail6PIDetaphitpc[aniter]);
					  }
						else
						  {
							 if(ichg==9) dtc7etaphitpc[aniter]->AddCutMonitor(cutPassPtPDGetaphitpc[aniter],cutFailPtPDGetaphitpc[aniter]);
							 if(ichg==10) dtc8etaphitpc[aniter]->AddCutMonitor(cutPassPtPDGetaphitpc[aniter],cutFailPtPDGetaphitpc[aniter]);
							 if(ichg==11) dtc9etaphitpc[aniter]->AddCutMonitor(cutPassPtPDGetaphitpc[aniter],cutFailPtPDGetaphitpc[aniter]);
							 cutPass6YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutPass%stpcM%i", chrgs[ichg], imult),PionMass);
							 cutFail6YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutFail%stpcM%i", chrgs[ichg], imult),PionMass);
							 if(ichg==9) dtc7etaphitpc[aniter]->AddCutMonitor(cutPass6YPtetaphitpc[aniter], cutFail6YPtetaphitpc[aniter]);
							 if(ichg==10) dtc8etaphitpc[aniter]->AddCutMonitor(cutPass6YPtetaphitpc[aniter], cutFail6YPtetaphitpc[aniter]);
							 if(ichg==11) dtc9etaphitpc[aniter]->AddCutMonitor(cutPass6YPtetaphitpc[aniter], cutFail6YPtetaphitpc[aniter]);

							 cutPass6PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cutPass%stpcM%i", chrgs[ichg], imult),0);//0-pion,1-kaon,2-proton
							 cutFail6PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cutFail%stpcM%i", chrgs[ichg], imult),0);
				 			    if(ichg==9) dtc7etaphitpc[aniter]->AddCutMonitor(cutPass6PIDetaphitpc[aniter], cutFail6PIDetaphitpc[aniter]);
							 if(ichg==10) dtc8etaphitpc[aniter]->AddCutMonitor(cutPass6PIDetaphitpc[aniter], cutFail6PIDetaphitpc[aniter]);
							 if(ichg==11) dtc9etaphitpc[aniter]->AddCutMonitor(cutPass6PIDetaphitpc[aniter], cutFail6PIDetaphitpc[aniter]);
						  }
*/
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
				
					if(ichg==0) //PP
					{
						anetaphitpc[aniter]->SetFirstParticleCut(dtc1etaphitpc[aniter]);
						anetaphitpc[aniter]->SetSecondParticleCut(dtc1etaphitpc[aniter]);
						//cMinvMonitor[aniter] = new AliFemtoCorrFctnMassInvMonitor(Form("%sM%i",chrgs[ichg], imult),500,ProtonMass);
						//anetaphitpc[aniter]->AddCorrFctn(cMinvMonitor[aniter]);
					}
					if(ichg==1) //aPaP
					{
						anetaphitpc[aniter]->SetFirstParticleCut(dtc2etaphitpc[aniter]);
						anetaphitpc[aniter]->SetSecondParticleCut(dtc2etaphitpc[aniter]);
				//		cMinvMonitor[aniter] = new AliFemtoCorrFctnMassInvMonitor(Form("%sM%i",chrgs[ichg], imult),500,ProtonMass);
						//anetaphitpc[aniter]->AddCorrFctn(cMinvMonitor[aniter]);
					}
					if(ichg==2) //PaP
					{
						anetaphitpc[aniter]->SetFirstParticleCut(dtc1etaphitpc[aniter]);
						anetaphitpc[aniter]->SetSecondParticleCut(dtc2etaphitpc[aniter]);
						//cMinvMonitor[aniter] = new AliFemtoCorrFctnMassInvMonitor(Form("%sM%i",chrgs[ichg], imult),500,ProtonMass);
						//anetaphitpc[aniter]->AddCorrFctn(cMinvMonitor[aniter]);
					}
					if(ichg==3) //KpKp
					{
						anetaphitpc[aniter]->SetFirstParticleCut(dtc1etaphitpc[aniter]);
						anetaphitpc[aniter]->SetSecondParticleCut(dtc1etaphitpc[aniter]);
						//cMinvMonitor[aniter] = new AliFemtoCorrFctnMassInvMonitor(Form("%sM%i",chrgs[ichg], imult),500,KaonMass);
						//anetaphitpc[aniter]->AddCorrFctn(cMinvMonitor[aniter]);
					}
					if(ichg==4) //KmKm
					{
						anetaphitpc[aniter]->SetFirstParticleCut(dtc2etaphitpc[aniter]);
						anetaphitpc[aniter]->SetSecondParticleCut(dtc2etaphitpc[aniter]);
						//cMinvMonitor[aniter] = new AliFemtoCorrFctnMassInvMonitor(Form("%sM%i",chrgs[ichg], imult),500,KaonMass);
						//anetaphitpc[aniter]->AddCorrFctn(cMinvMonitor[aniter]);
					
					}
					if(ichg==5) //KpKm
					{
						anetaphitpc[aniter]->SetFirstParticleCut(dtc1etaphitpc[aniter]);
						anetaphitpc[aniter]->SetSecondParticleCut(dtc2etaphitpc[aniter]);
						//cMinvMonitor[aniter] = new AliFemtoCorrFctnMassInvMonitor(Form("%sM%i",chrgs[ichg], imult),500,KaonMass);
						//anetaphitpc[aniter]->AddCorrFctn(cMinvMonitor[aniter]);
					}
					if(ichg==6)//PIpPIp
					{
						anetaphitpc[aniter]->SetFirstParticleCut(dtc1etaphitpc[aniter]);
						anetaphitpc[aniter]->SetSecondParticleCut(dtc1etaphitpc[aniter]);
						//cMinvMonitor[aniter] = new AliFemtoCorrFctnMassInvMonitor(Form("%sM%i",chrgs[ichg], imult),500,PionMass);
						//anetaphitpc[aniter]->AddCorrFctn(cMinvMonitor[aniter]);
					}
					if(ichg==7)//PImPIm
					{
						anetaphitpc[aniter]->SetFirstParticleCut(dtc2etaphitpc[aniter]);
						anetaphitpc[aniter]->SetSecondParticleCut(dtc2etaphitpc[aniter]);
						//cMinvMonitor[aniter] = new AliFemtoCorrFctnMassInvMonitor(Form("%sM%i",chrgs[ichg], imult),500,PionMass);
						//anetaphitpc[aniter]->AddCorrFctn(cMinvMonitor[aniter]);
					
					}
					if(ichg==8)//PIpPIm
					{
						anetaphitpc[aniter]->SetFirstParticleCut(dtc1etaphitpc[aniter]);
						anetaphitpc[aniter]->SetSecondParticleCut(dtc2etaphitpc[aniter]);
						//cMinvMonitor[aniter] = new AliFemtoCorrFctnMassInvMonitor(Form("%sM%i",chrgs[ichg], imult),500,PionMass);
						//anetaphitpc[aniter]->AddCorrFctn(cMinvMonitor[aniter]);
					}

		             if(ichg==9)//all
		             {
						 	anetaphitpc[aniter]->SetFirstParticleCut(dtc3etaphitpc[aniter]);
						 	anetaphitpc[aniter]->SetSecondParticleCut(dtc3etaphitpc[aniter]);
                   }

                   if(ichg==10)//plus
                   {
                   	anetaphitpc[aniter]->SetFirstParticleCut(dtc1etaphitpc[aniter]);
							anetaphitpc[aniter]->SetSecondParticleCut(dtc1etaphitpc[aniter]);
                   }

                   if(ichg==11)//minus
                   {
                    	anetaphitpc[aniter]->SetFirstParticleCut(dtc2etaphitpc[aniter]);
							anetaphitpc[aniter]->SetSecondParticleCut(dtc2etaphitpc[aniter]);
                   }
                   if(ichg==12)//mixed
                   {
			            anetaphitpc[aniter]->SetFirstParticleCut(dtc1etaphitpc[aniter]);
			            anetaphitpc[aniter]->SetSecondParticleCut(dtc2etaphitpc[aniter]);
			          //  cMinvMonitor[aniter] = new AliFemtoCorrFctnMassInvMonitor(Form("%sM%i",chrgs[ichg], imult),500,PionMass);
			          //  anetaphitpc[aniter]->AddCorrFctn(cMinvMonitor[aniter]);
						 }

					anetaphitpc[aniter]->SetPairCut(sqpcetaphitpc[aniter]);
					cdedpetaphi[aniter] = new AliFemtoCorrFctnDEtaDPhi(Form("cdedp%stpcM%i", chrgs[ichg], imult),35, 35);
					anetaphitpc[aniter]->AddCorrFctn(cdedpetaphi[aniter]);

					if (runktdep)
					{
						
					}		
					Manager->AddAnalysis(anetaphitpc[aniter]);	
				}
			}
		}
	}

	return Manager;
}												 
