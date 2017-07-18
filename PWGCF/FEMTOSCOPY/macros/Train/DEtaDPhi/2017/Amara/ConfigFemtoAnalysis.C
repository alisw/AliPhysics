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
	const int numOfChTypes = 20;
	const int numOfpTbins = 5;


	int runmults[numOfMultBins] = {0, 0, 0, 0, 1};
	int multbins[numOfMultBins+1] = {2, 20, 50,150,2,1500};
	
	int runch[numOfChTypes] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
	const char *chrgs[numOfChTypes] = { "PP", "aPaP", "PaP", "KpKp", "KmKm", "KpKm", "PIpPIp", "PImPIm", "PIpPIm", "all", "plus", "minus", "mixed",  "V0PL","V0PAL","V0APL","V0APAL","V0LL","V0ALAL","V0LAL" };
	
	
	
	int runptdep = 1;
	double ptrng[numOfpTbins+1] = {0.0, 0, 0, 0, 0, 0};
	double ptrngAll[numOfpTbins+1] = {0.0, 1.0, 2.0, 3.0, 4.0, 100.0};
	double ptrngPion[numOfpTbins+1] = {0.0, 0.8, 1.2, 1.4, 2.5, 100.0};
	double ptrngKaon[numOfpTbins+1] = {0.0, 1.5, 2.5, 3.5, 100.0, 0};
	double ptrngProton[numOfpTbins+1] = {0.0, 2.75, 100, 0, 0, 0};

	AliFemtoEventReaderKinematicsChain* Reader=new AliFemtoEventReaderKinematicsChain();

	AliFemtoManager* Manager = new AliFemtoManager();
	Manager->SetEventReader(Reader);



	AliFemtoVertexMultAnalysis	*anetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoBasicEventCut		*mecetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoCutMonitorEventMult	*cutPassEvMetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoCutMonitorEventMult	*cutFailEvMetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoCutMonitorEventVertex   *cutPassEvVetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoCutMonitorEventVertex   *cutFailEvVetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoMCTrackCut		*dtc1etaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoMCTrackCut		*dtc2etaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoMCTrackCut		*dtc3etaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoCutMonitorParticleYPt   *cutPass1YPtetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoCutMonitorParticleYPt   *cutFail1YPtetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoCutMonitorParticlePID   *cutPass1PIDetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoCutMonitorParticlePID   *cutFail1PIDetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoCutMonitorParticleYPt   *cutPass2YPtetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoCutMonitorParticleYPt   *cutFail2YPtetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoCutMonitorParticlePID   *cutPass2PIDetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoCutMonitorParticlePID   *cutFail2PIDetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoCutMonitorParticleYPt   *cutPass3YPtetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoCutMonitorParticleYPt   *cutFail3YPtetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoCutMonitorParticlePID   *cutPass3PIDetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoCutMonitorParticlePID   *cutFail3PIDetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoPairCutAntiGamma	*sqpcetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoPairCutPt		*ptpcuts[numOfMultBins*numOfChTypes*numOfpTbins];
	AliFemtoQinvCorrFctn		*cqinvkttpc[numOfMultBins*numOfChTypes*numOfpTbins];
	AliFemtoQinvCorrFctn		*cqinvtpc[numOfMultBins*numOfChTypes];
	AliFemtoCorrFctnDEtaDPhiSimple	*cdedpetaphi[numOfMultBins*numOfChTypes];
	AliFemtoCorrFctnDEtaDPhiSimple	*cdedpetaphiPt[numOfMultBins*numOfChTypes*numOfpTbins];
	AliFemtoCorrFctnNonIdDR         *cnonidtpc[numOfMultBins*numOfChTypes];

	
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
					int multmix = 5;
					if(imult == 4) multmix = 30; 
					anetaphitpc[aniter] = new AliFemtoVertexMultAnalysis(10, -10.0, 10.0, multmix, multbins[imult], multbins[imult+1]);
					anetaphitpc[aniter]->SetNumEventsToMix(10);
					anetaphitpc[aniter]->SetMinSizePartCollection(1);

					//*** Event cut ***
					mecetaphitpc[aniter] = new AliFemtoBasicEventCut();
					mecetaphitpc[aniter]->SetEventMult(0.001,100000);
					mecetaphitpc[aniter]->SetVertZPos(-10,10);//cm

					//****** event monitors **********	
					cutPassEvMetaphitpc[aniter] = new AliFemtoCutMonitorEventMult(Form("cutPass%stpcM%i", chrgs[ichg], imult));
					cutFailEvMetaphitpc[aniter] = new AliFemtoCutMonitorEventMult(Form("cutFail%stpcM%i", chrgs[ichg], imult));
					mecetaphitpc[aniter]->AddCutMonitor(cutPassEvMetaphitpc[aniter], cutFailEvMetaphitpc[aniter]);
		
					cutPassEvVetaphitpc[aniter] = new AliFemtoCutMonitorEventVertex(Form("cutPass%stpcM%i", chrgs[ichg], imult));
					cutFailEvVetaphitpc[aniter] = new AliFemtoCutMonitorEventVertex(Form("cutFail%stpcM%i", chrgs[ichg], imult));
					mecetaphitpc[aniter]->AddCutMonitor(cutPassEvVetaphitpc[aniter], cutFailEvVetaphitpc[aniter]);


					// ***** single particle track cuts *********
					dtc1etaphitpc[aniter] = new AliFemtoMCTrackCut();
					dtc2etaphitpc[aniter] = new AliFemtoMCTrackCut();
					dtc3etaphitpc[aniter] = new AliFemtoMCTrackCut();

					dtc1etaphitpc[aniter]->SetCharge(1.0);
					dtc2etaphitpc[aniter]->SetCharge(-1.0);

					dtc1etaphitpc[aniter]->SetEta(-0.8,0.8);
					dtc2etaphitpc[aniter]->SetEta(-0.8,0.8);
					dtc3etaphitpc[aniter]->SetEta(-0.8,0.8);


                                        if (ichg == 0 ||ichg == 1 ||ichg == 2)//protons 0-2
                                          {
                                            dtc1etaphitpc[aniter]->SetPt(0.5,2.5);
                                            dtc2etaphitpc[aniter]->SetPt(0.5,2.5);
                                           
					    dtc1etaphitpc[aniter]->SetPDG(2212);
					    dtc2etaphitpc[aniter]->SetPDG(2212);
                                          }

					if (ichg == 3 ||ichg == 4 ||ichg == 5)//kaons 3-5
                                          {
                                            dtc1etaphitpc[aniter]->SetPt(0.3,2.5);
                                            dtc2etaphitpc[aniter]->SetPt(0.3,2.5);

					    dtc1etaphitpc[aniter]->SetPDG(321);
					    dtc2etaphitpc[aniter]->SetPDG(321);

                                          }
                                        if (ichg == 6 ||ichg == 7 ||ichg == 8)//pions 6-8
                                          {
                                            dtc1etaphitpc[aniter]->SetPt(0.2,2.5);
					    dtc2etaphitpc[aniter]->SetPt(0.2,2.5);

					    dtc1etaphitpc[aniter]->SetPDG(211);
					    dtc2etaphitpc[aniter]->SetPDG(211);	
                                          }
                                        if (ichg == 9)//all
                                          {
                                            dtc3etaphitpc[aniter]->SetPt(0.2,2.5);
                                          }
                                        if (ichg == 10 ||ichg == 11 ||ichg == 12)//plus,minus,mixed
                                          {
                                            dtc1etaphitpc[aniter]->SetPt(0.2,2.5);
                                            dtc2etaphitpc[aniter]->SetPt(0.2,2.5);
                                          }
					if(ichg >= 13 && ichg <=16){
					  dtc1etaphitpc[aniter]->SetPt(0.6,2.5); //lambda
					  dtc1etaphitpc[aniter]->SetCharge(0.0);
					  if(ichg == 14 || ichg ==16) dtc1etaphitpc[aniter]->SetPDG(-3122);
					  else dtc1etaphitpc[aniter]->SetPDG(3122);
					  
					  dtc2etaphitpc[aniter]->SetPt(0.5,2.5);//proton
					  dtc2etaphitpc[aniter]->SetPDG(2212);
					  if(ichg == 15 || ichg ==16) dtc2etaphitpc[aniter]->SetCharge(-1.0);
					  else dtc2etaphitpc[aniter]->SetCharge(1.0);
					  
					}
					if (ichg == 17 ||ichg == 18 ||ichg == 19)//lambdas
                                          {
                                            dtc1etaphitpc[aniter]->SetPt(0.6,2.5);
					    dtc2etaphitpc[aniter]->SetPt(0.6,2.5);
					    dtc1etaphitpc[aniter]->SetCharge(0.0);
					    dtc2etaphitpc[aniter]->SetCharge(0.0);
					    dtc1etaphitpc[aniter]->SetPDG(3122);
					    dtc2etaphitpc[aniter]->SetPDG(-3122);	
					  }

				
					//**************** track Monitors ***************

					
					if(1)
					  {
					    cutPass3YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutPass%stpcM%i", chrgs[ichg], imult),PionMass);
					    cutFail3YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutFail%stpcM%i", chrgs[ichg], imult),PionMass);
					    if(ichg==9) dtc3etaphitpc[aniter]->AddCutMonitor(cutPass3YPtetaphitpc[aniter], cutFail3YPtetaphitpc[aniter]);
					    if(ichg==0||ichg==3||ichg==6||ichg==10) dtc1etaphitpc[aniter]->AddCutMonitor(cutPass3YPtetaphitpc[aniter], cutFail3YPtetaphitpc[aniter]);
					   
					  }
					 
					//******** Two - track cuts ************
					sqpcetaphitpc[aniter] = new AliFemtoPairCutAntiGamma();
					sqpcetaphitpc[aniter]->SetDataType(AliFemtoPairCut::kKine);				      
					
		
					//***** Setting cuts ***********

		
					// setting event cut
					anetaphitpc[aniter]->SetEventCut(mecetaphitpc[aniter]);
					//setting single track cuts
					if(ichg==0 || ichg==3  || ichg==6 || ichg==10 || ichg==17) //positive like-sign
					  {
					    anetaphitpc[aniter]->SetFirstParticleCut(dtc1etaphitpc[aniter]);
					    anetaphitpc[aniter]->SetSecondParticleCut(dtc1etaphitpc[aniter]);
					  }
					if(ichg==1 || ichg==4 || ichg==7 || ichg==11 || ichg == 18)//negative like-sign
					  {
					    anetaphitpc[aniter]->SetFirstParticleCut(dtc2etaphitpc[aniter]);
					    anetaphitpc[aniter]->SetSecondParticleCut(dtc2etaphitpc[aniter]);
					  }
					if(ichg==2 || ichg==5 || ichg==8 || ichg==12 || (ichg >= 13 && ichg <=16) || ichg == 19)//unlike-sign + proton/lambda
					  {
					    anetaphitpc[aniter]->SetFirstParticleCut(dtc1etaphitpc[aniter]);
					    anetaphitpc[aniter]->SetSecondParticleCut(dtc2etaphitpc[aniter]);
					  }
					if(ichg==9) //all
					  {
					    anetaphitpc[aniter]->SetFirstParticleCut(dtc3etaphitpc[aniter]);
					    anetaphitpc[aniter]->SetSecondParticleCut(dtc3etaphitpc[aniter]);
					  }
				      

					//setting two-track cuts
					anetaphitpc[aniter]->SetPairCut(sqpcetaphitpc[aniter]);


					//**** Correlation functions *******


					if(ichg >= 13)
					  cdedpetaphi[aniter] = new AliFemtoCorrFctnDEtaDPhiSimple(Form("cdedp%stpcM%i", chrgs[ichg], imult),23, 23);
					else
					  cdedpetaphi[aniter] = new AliFemtoCorrFctnDEtaDPhiSimple(Form("cdedp%stpcM%i", chrgs[ichg], imult),29, 29);
					anetaphitpc[aniter]->AddCorrFctn(cdedpetaphi[aniter]);





					//**** Pt dependence *******
					if (runptdep)
					  {
					    if(ichg<=2){
					      for(int pit=0;pit<=numOfpTbins;pit++)
						ptrng[pit]=ptrngProton[pit];
					    }
					    else if(ichg>2 && ichg<6){
					      for(int pit=0;pit<=numOfpTbins;pit++)
						ptrng[pit]=ptrngKaon[pit];
					    }
					    else if(ichg>=6 && ichg<=8){
					      for(int pit=0;pit<=numOfpTbins;pit++)
						ptrng[pit]=ptrngPion[pit];
					    }
					    else if(ichg>=9){
					      for(int pit=0;pit<=numOfpTbins;pit++)
						ptrng[pit]=ptrngAll[pit];
					    }

					    int ptm;
					    for (int ipt=0; ipt<numOfpTbins; ipt++)
					      {
						if(ptrng[ipt+1]==0) continue;
						ptm = aniter * numOfpTbins + ipt;
						ptpcuts[ptm] = new AliFemtoPairCutPt(ptrng[ipt], ptrng[ipt+1]);
				

						cdedpetaphiPt[ptm] = new AliFemtoCorrFctnDEtaDPhiSimple(Form("cdedp%stpcM%ipT%i", chrgs[ichg], imult,ipt),23, 23);
						cdedpetaphiPt[ptm]->SetPairSelectionCut(ptpcuts[ptm]);
						anetaphitpc[aniter]->AddCorrFctn(cdedpetaphiPt[ptm]);

					      }
					  }
					

		
					Manager->AddAnalysis(anetaphitpc[aniter]);	
				}
			}
		}
	}
	return Manager;
}												 
