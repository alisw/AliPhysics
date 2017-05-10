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
#include "AliFemtoCorrFctnDirectYlm.h"
#include "AliFemtoCutMonitorParticlePtPDG.h"
#include "AliFemtoKTPairCut.h"
#include "AliFemtoPairCutPt.h"
#endif

//________________________________________________________________________
AliFemtoManager* ConfigFemtoAnalysis(const char* params) {

	double PionMass = 0.13956995;
	double KaonMass = 0.493677;
	double ProtonMass = 0.938272013;
	
	const int numOfMultBins = 5;	
	const int numOfChTypes = 20;
	const int numOfkTbins = 5;

	char *parameter[20];
	if(strlen(params)!=0)
	  {
	    parameter[0] = strtok(params, ","); // Splits spaces between words in params
	    cout<<"Parameter [0] (filterbit):"<<parameter[0]<<endl; // Writes first parameter
	    parameter[1] = strtok(NULL, ",");
	    cout<<"Parameter [1] (ktdep):"<<parameter[1]<<" "<<endl;
	    parameter[2] = strtok(NULL, ",");
	    cout<<"Parameter [2] (multdep):"<<parameter[2]<<" "<<endl;	  
	    parameter[3] = strtok(NULL, ",");
	    cout<<"Parameter [3]: (MinPlpContribSPD)"<<parameter[3]<<" "<<endl;
	    parameter[4] = strtok(NULL, ",");
	    cout<<"Parameter [4]: (multbino)"<<parameter[4]<<" "<<endl;
	    parameter[5] = strtok(NULL, ",");
	    cout<<"Parameter [5]: (zvertbino)"<<parameter[5]<<" "<<endl;
	    parameter[6] = strtok(NULL, ",");
	    cout<<"Parameter [6]: (ifGlobalTracks=true/false)"<<parameter[6]<<" "<<endl;
	    parameter[7] = strtok(NULL, ",");
	    cout<<"Parameter [7]: (shareQuality)"<<parameter[7]<<" "<<endl;
	    parameter[8] = strtok(NULL, ",");
	    cout<<"Parameter [8]: (shareFraction)"<<parameter[8]<<" "<<endl;
	    parameter[9] = strtok(NULL, ",");
	    cout<<"Parameter [9]: (ifElectronRejection)"<<parameter[9]<<" "<<endl;
	    parameter[10] = strtok(NULL, ",");
	    cout<<"Parameter [10]: (nSigma)"<<parameter[10]<<" "<<endl;
	    parameter[11] = strtok(NULL, ",");
	    cout<<"Parameter [12]: (etaMin)"<<parameter[11]<<" "<<endl;
	    parameter[12] = strtok(NULL, ",");
	    cout<<"Parameter [12]: (etaMax)"<<parameter[12]<<" "<<endl;
	    parameter[13] = strtok(NULL, ",");
	    cout<<"Parameter [13]: (ispileup)"<<parameter[13]<<" "<<endl;
	    parameter[14] = strtok(NULL, ",");
	    cout<<"Parameter [14]: (max pT)"<<parameter[14]<<" "<<endl;
	    parameter[15] = strtok(NULL, ",");
	    cout<<"Parameter [15]: (SetMostProbable 1)"<<parameter[15]<<" "<<endl;
	    parameter[16] = strtok(NULL, ",");
	    cout<<"Parameter [16]: (SetMostProbable 2)"<<parameter[16]<<" "<<endl;
	    parameter[17] = strtok(NULL, ",");
	    cout<<"Parameter [17]: (SetMostProbable 3)"<<parameter[17]<<" "<<endl;
	    parameter[18] = strtok(NULL, ",");
	    cout<<"Parameter [18]: (FILE no)"<<parameter[18]<<" "<<endl;
	  }
	int filterbit = atoi(parameter[0]); //96 / 768 / 128 
	int runktdep = atoi(parameter[1]); //0
	int runmultdep = atoi(parameter[2]); //0
	int minPlpContribSPD = atoi(parameter[3]); //3
	int multbino = atoi(parameter[4]); //30
	int zvertbino = atoi(parameter[5]); //10
	Bool_t ifGlobalTracks=kFALSE; if(atoi(parameter[6]))ifGlobalTracks=kTRUE;//kTRUE 
	double shareQuality = atof(parameter[7]); //0.00
	double shareFraction = atof(parameter[8]); //0.05
	bool ifElectronRejection = atoi(parameter[9]); //true
	double nSigmaVal = atof(parameter[10]); //3.0
	double nEtaMin = atof(parameter[11]); //-0.8
	double nEtaMax = atof(parameter[12]);  //0.8
	bool ifIsPileUp = atoi(parameter[13]); //true
	double maxPt = atof(parameter[14]);  //4.0

	int setMostProb1 = atoi(parameter[15]);
	int setMostProb2 = atoi(parameter[16]);
	int setMostProb3 = atoi(parameter[17]);

	int fileNo = atoi(parameter[18]);
	char* fileName[300];
	//if(fileNo==0)
	//  strcpy(fileName,"alien:///alice/cern.ch/user/m/majanik/2014/DEtaDPhi/Trains/Corrections/Train7Light/1Dmap_FB96_MCOnly_DoubleCounting.root");
	//cout<<"Filename: "<<Form("%s",fileName)<<endl;

	//printf("*** Connect to AliEn ***\n");
	//TGrid::Connect("alien://");

	int runmults[numOfMultBins] = {0, 0, 0, 0, 1};
	int multbins[numOfMultBins+1] = {2, 20, 50,150,2,1500};
	
	int runch[numOfChTypes] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
	const char *chrgs[numOfChTypes] = { "PP", "aPaP", "PaP", "KpKp", "KmKm", "KpKm", "PIpPIp", "PImPIm", "PIpPIm", "all", "plus", "minus", "mixed",  "V0PL","V0PAL","V0APL","V0APAL","V0LL","V0ALAL","V0LAL" };
	
	//int runktdep = 1;
	double ktrng[numOfkTbins+1] = {0.0, 0, 0, 0, 0, 0};
	double ktrngAll[numOfkTbins+1] = {0.0, 1.0, 2.0, 3.0, 4.0, 100.0};
	double ktrngPion[numOfkTbins+1] = {0.0, 0.8, 1.2, 1.4, 2.5, 100.0};
	double ktrngKaon[numOfkTbins+1] = {0.0, 1.5, 2.5, 3.5, 100.0, 0};
	double ktrngProton[numOfkTbins+1] = {0.0, 2.75, 100, 0, 0, 0};

	int runqinv = 1;
	int runshlcms = 1;// 0:PRF(PAP), 1:LCMS(PP,APAP)

	int runtype = 0; // Types 0 - global, 1 - ITS only, 2 - TPC Inner	//global tracks ->mfit ITS+TPC
	int owncuts = 1; 
	int owndca = 1;

	int gammacut = 0;	// cut na ee z gamma 
	
	double shqmax = 1.0; 
	int nbinssh = 100;

	// AliFemtoEventReaderESDChain *Reader = new AliFemtoEventReaderESDChain();
	// Reader->SetUseMultiplicity(AliFemtoEventReaderESDChain::kReferenceITSTPC);

	//AliFemtoEventReaderAODChain *Reader = new AliFemtoEventReaderAODChain();
	//Reader->SetFilterMask(96);
	//Reader->SetDCAglobalTrack(kTRUE); //false for FB7, true for the rest //we do not use DCA at all
	//Reader->SetUseMultiplicity(AliFemtoEventReaderAOD::kReference);

	// Reader->SetMinPlpContribSPD(3);
	// Reader->SetIsPileUpEvent(kTRUE);



	AliFemtoEventReaderKinematicsChain* Reader=new AliFemtoEventReaderKinematicsChain();
	Reader->RemoveWeakDecaysManually(kFALSE);
	Reader->ReadOnlyPrimaries(kFALSE);
	Reader->ReadOnlyPrimariesV0(kFALSE);
	
	
	AliFemtoManager* Manager = new AliFemtoManager();
	Manager->SetEventReader(Reader);

	AliFemtoVertexMultAnalysis		*anetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoBasicEventCut				 *mecetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoCutMonitorEventMult	 *cutPassEvMetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoCutMonitorEventMult	 *cutFailEvMetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoCutMonitorEventVertex *cutPassEvVetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoCutMonitorEventVertex *cutFailEvVetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoMCTrackCut			 *dtc1etaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoMCTrackCut			 *dtc2etaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoMCTrackCut			 *dtc3etaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoCutMonitorParticleYPt *cutPass1YPtetaphitpc[numOfMultBins*numOfChTypes];
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
	//	 AliFemtoShareQualityTPCEntranceSepPairCut			*sqpcetaphitpcsame[numOfMultBins*numOfChTypes];
	//AliFemtoPairCutAntiGamma			*sqpcetaphitpc[numOfMultBins*numOfChTypes];
	//AliFemtoPairCutRadialDistance			*sqpcetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoPairCutAntiGamma			*sqpcetaphitpc[numOfMultBins*numOfChTypes];
	//	AliFemtoChi2CorrFctn					*cchiqinvetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoPairCutPt					*ktpcuts[numOfMultBins*numOfChTypes*numOfkTbins];
	AliFemtoQinvCorrFctn					*cqinvkttpc[numOfMultBins*numOfChTypes*numOfkTbins];
	AliFemtoQinvCorrFctn					*cqinvtpc[numOfMultBins*numOfChTypes];
	AliFemtoCorrFctnDEtaDPhi			*cdedpetaphi[numOfMultBins*numOfChTypes];
	AliFemtoCorrFctnDEtaDPhi			*cdedpetaphiPt[numOfMultBins*numOfChTypes*numOfkTbins];
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
					anetaphitpc[aniter] = new AliFemtoVertexMultAnalysis(zvertbino, -10.0, 10.0, multbino, multbins[imult], multbins[imult+1]);
					anetaphitpc[aniter]->SetNumEventsToMix(10);
					anetaphitpc[aniter]->SetMinSizePartCollection(1);
					anetaphitpc[aniter]->SetVerboseMode(kFALSE);//~~~~~~~~~~~~~~~~

					//*** Event cut ***
					mecetaphitpc[aniter] = new AliFemtoBasicEventCut();
					mecetaphitpc[aniter]->SetEventMult(0.001,100000);
					mecetaphitpc[aniter]->SetVertZPos(-10,10);//cm

					//****** event monitors **********	
					cutPassEvMetaphitpc[aniter] = new AliFemtoCutMonitorEventMult(Form("cutPass%stpcM%i", chrgs[ichg], imult));
					cutFailEvMetaphitpc[aniter] = new AliFemtoCutMonitorEventMult(Form("cutFail%stpcM%i", chrgs[ichg], imult));
					mecetaphitpc[aniter]->AddCutMonitor(cutPassEvMetaphitpc[aniter], cutFailEvMetaphitpc[aniter]);
		
					//cutPassEvVetaphitpc[aniter] = new AliFemtoCutMonitorEventVertex(Form("cutPass%stpcM%i", chrgs[ichg], imult));
					//cutFailEvVetaphitpc[aniter] = new AliFemtoCutMonitorEventVertex(Form("cutFail%stpcM%i", chrgs[ichg], imult));
					//mecetaphitpc[aniter]->AddCutMonitor(cutPassEvVetaphitpc[aniter], cutFailEvVetaphitpc[aniter]);


					// ***** single particle track cuts *********
					dtc1etaphitpc[aniter] = new AliFemtoMCTrackCut();
					dtc2etaphitpc[aniter] = new AliFemtoMCTrackCut();
					dtc3etaphitpc[aniter] = new AliFemtoMCTrackCut();

					dtc1etaphitpc[aniter]->SetCharge(1.0);
					dtc2etaphitpc[aniter]->SetCharge(-1.0);

					dtc1etaphitpc[aniter]->SetEta(nEtaMin,nEtaMax);
					dtc2etaphitpc[aniter]->SetEta(nEtaMin,nEtaMax);
					dtc3etaphitpc[aniter]->SetEta(nEtaMin,nEtaMax);

					// dtc1etaphitpc[aniter]->SetElectronRejection(true);	
					// dtc2etaphitpc[aniter]->SetElectronRejection(true);
					// dtc3etaphitpc[aniter]->SetElectronRejection(true); 

                                        if (ichg == 0 ||ichg == 1 ||ichg == 2)//protons 0-2
                                          {
                                            dtc1etaphitpc[aniter]->SetPt(0.5,maxPt);
                                            dtc2etaphitpc[aniter]->SetPt(0.5,maxPt);
                                           
					    dtc1etaphitpc[aniter]->SetPDG(2212);
					    dtc2etaphitpc[aniter]->SetPDG(2212);
                                          }

					if (ichg == 3 ||ichg == 4 ||ichg == 5)//kaons 3-5
                                          {
                                            dtc1etaphitpc[aniter]->SetPt(0.3,maxPt);
                                            dtc2etaphitpc[aniter]->SetPt(0.3,maxPt);

					    dtc1etaphitpc[aniter]->SetPDG(321);
					    dtc2etaphitpc[aniter]->SetPDG(321);

                                          }
                                        if (ichg == 6 ||ichg == 7 ||ichg == 8)//pions 6-8
                                          {
                                            dtc1etaphitpc[aniter]->SetPt(0.2,maxPt);
					    dtc2etaphitpc[aniter]->SetPt(0.2,maxPt);

					    dtc1etaphitpc[aniter]->SetPDG(211);
					    dtc2etaphitpc[aniter]->SetPDG(211);	
                                          }
                                        if (ichg == 9)//all
                                          {
                                            dtc3etaphitpc[aniter]->SetPt(0.2,maxPt);
                                          }
                                        if (ichg == 10 ||ichg == 11 ||ichg == 12)//plus,minus,mixed
                                          {
                                            dtc1etaphitpc[aniter]->SetPt(0.2,maxPt);
                                            dtc2etaphitpc[aniter]->SetPt(0.2,maxPt);
                                          }
					if(ichg >= 13 && ichg <=16){
					  dtc1etaphitpc[aniter]->SetPt(0.6,maxPt); //lambda
					  dtc1etaphitpc[aniter]->SetCharge(0.0);
					  if(ichg == 14 || ichg ==16) dtc1etaphitpc[aniter]->SetPDG(-3122);
					  else dtc1etaphitpc[aniter]->SetPDG(3122);
					  
					  dtc2etaphitpc[aniter]->SetPt(0.5,maxPt);//proton
					  dtc2etaphitpc[aniter]->SetPDG(2212);
					  if(ichg == 15 || ichg ==16) dtc2etaphitpc[aniter]->SetCharge(-1.0);
					  else dtc2etaphitpc[aniter]->SetCharge(1.0);
					  
					}
					if (ichg == 17 ||ichg == 18 ||ichg == 19)//lambdas
                                          {
                                            dtc1etaphitpc[aniter]->SetPt(0.6,maxPt);
					    dtc2etaphitpc[aniter]->SetPt(0.6,maxPt);
					    dtc1etaphitpc[aniter]->SetCharge(0.0);
					    dtc2etaphitpc[aniter]->SetCharge(0.0);
					    dtc1etaphitpc[aniter]->SetPDG(3122);
					    dtc2etaphitpc[aniter]->SetPDG(-3122);	
					  }
				
					//**************** track Monitors ***************

					
					if(1)//ichg>8)
					  {
					    cutPass3YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutPass%stpcM%i", chrgs[ichg], imult),PionMass);
					    cutFail3YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutFail%stpcM%i", chrgs[ichg], imult),PionMass);
					    if(ichg==9) dtc3etaphitpc[aniter]->AddCutMonitor(cutPass3YPtetaphitpc[aniter], cutFail3YPtetaphitpc[aniter]);
					    if(ichg==0||ichg==3||ichg==6||ichg==10) dtc1etaphitpc[aniter]->AddCutMonitor(cutPass3YPtetaphitpc[aniter], cutFail3YPtetaphitpc[aniter]);
					    //if(ichg==1||ichg==4||ichg==7||ichg==11) dtc2etaphitpc[aniter]->AddCutMonitor(cutPass3YPtetaphitpc[aniter], cutFail3YPtetaphitpc[aniter]);

					                              
					    // cutPass1PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cutPass%stpcM%i", chrgs[ichg], imult),0);
					    // cutFail1PIDetaphitpc[aniter]  = new AliFemtoCutMonitorParticlePID(Form("cutFail%stpcM%i", chrgs[ichg], imult),0);
					    // if(ichg==0||ichg==3||ichg==6||ichg==10) dtc1etaphitpc[aniter]->AddCutMonitor(cutPass1PIDetaphitpc[aniter], cutFail1PIDetaphitpc[aniter]);


					  }
					 
					//******** Two - track cuts ************
					sqpcetaphitpc[aniter] = new AliFemtoPairCutAntiGamma();

					//sqpcetaphitpc[aniter] = new AliFemtoPairCutRadialDistance();
					sqpcetaphitpc[aniter]->SetDataType(AliFemtoPairCut::kKine);
					
					/*
					sqpcetaphitpc[aniter]->SetShareQualityMax(0.0);		// two track cuts on splitting and merging  //1- wylaczany 0 -wlaczany   
					sqpcetaphitpc[aniter]->SetShareFractionMax(0.05);	//  ile moga miec wspolnych klastrow //1 - wylaczany, 0.05 - wlaczany
					sqpcetaphitpc[aniter]->SetRemoveSameLabel(kFALSE);
					sqpcetaphitpc[aniter]->SetMaximumRadius(0.82);
					sqpcetaphitpc[aniter]->SetMinimumRadius(0.8);
					sqpcetaphitpc[aniter]->SetPhiStarDifferenceMinimum(0.02);
					sqpcetaphitpc[aniter]->SetEtaDifferenceMinimum(0.02);

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
					*/
					// sqpcetaphitpc[aniter]->SetTPCEntranceSepMinimum(1.5);
					// sqpcetaphitpc[aniter]->SetRadialDistanceMinimum(0.12, 0.03);
					// sqpcetaphitpc[aniter]->SetEtaDifferenceMinimum(0.02);
					  
					
		

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
					  cdedpetaphi[aniter] = new AliFemtoCorrFctnDEtaDPhi(Form("cdedp%stpcM%i", chrgs[ichg], imult),23, 23);
					else
					  cdedpetaphi[aniter] = new AliFemtoCorrFctnDEtaDPhi(Form("cdedp%stpcM%i", chrgs[ichg], imult),29, 29);
					anetaphitpc[aniter]->AddCorrFctn(cdedpetaphi[aniter]);



					/*
					if(ichg==0 || ichg==1 || ichg==3 || ichg==4 || ichg==6 || ichg==7 || ichg==9 || ichg==10 || ichg==11 || ichg==12 || ichg==17 || ichg==18) //PP, aPaP, LL, ALAL
                                        {
					  cqinvtpc[aniter] = new AliFemtoQinvCorrFctn(Form("cqinv%stpcM%i", chrgs[ichg], imult),nbinssh,0.0,shqmax); //femto qinv, for identical mass particles
					  anetaphitpc[aniter]->AddCorrFctn(cqinvtpc[aniter]);
					}
					else //PaP, PL, APL, PAL, APAL, LAL
                                        {
					  cnonidtpc[aniter] = new AliFemtoCorrFctnNonIdDR(Form("cnonid%stpcM%i", chrgs[ichg], imult), nbinssh, 0.0,shqmax); //for non-identical partcles
					  anetaphitpc[aniter]->AddCorrFctn(cnonidtpc[aniter]);
                                        }
					*/



					

					
					if (runktdep)
					{
 					if(ichg<=2){
						for(int kit=0;kit<=numOfkTbins;kit++)
						    ktrng[kit]=ktrngProton[kit];
					  }
					else if(ichg>2&&ichg<6){
					       for(int kit=0;kit<=numOfkTbins;kit++)
					            ktrng[kit]=ktrngKaon[kit];
                                          }
					else if(ichg>=6&&ichg<=8){
					    for(int kit=0;kit<=numOfkTbins;kit++)
					       ktrng[kit]=ktrngPion[kit];
                                          }
					else if(ichg>=9){
					    for(int kit=0;kit<=numOfkTbins;kit++)
					       ktrng[kit]=ktrngAll[kit];
                                          }

						int ktm;
						for (int ikt=0; ikt<numOfkTbins; ikt++)
						{
							if(ktrng[ikt+1]==0) continue;
							ktm = aniter * numOfkTbins + ikt;
							ktpcuts[ktm] = new AliFemtoPairCutPt(ktrng[ikt], ktrng[ikt+1]);
				
							//cqinvkttpc[ktm] = new AliFemtoQinvCorrFctn(Form("cqinv%stpcM%ikT%i", chrgs[ichg], imult, ikt),nbinssh,0.0,(imult>6)?shqmax*2.5:shqmax);
							//cqinvkttpc[ktm]->SetPairSelectionCut(ktpcuts[ktm]);
							//anetaphitpc[aniter]->AddCorrFctn(cqinvkttpc[ktm]);

							cdedpetaphiPt[ktm] = new AliFemtoCorrFctnDEtaDPhi(Form("cdedp%stpcM%ipT%i", chrgs[ichg], imult,ikt),23, 23);
							cdedpetaphiPt[ktm]->SetPairSelectionCut(ktpcuts[ktm]);
							anetaphitpc[aniter]->AddCorrFctn(cdedpetaphiPt[ktm]);

						}
					}		
					Manager->AddAnalysis(anetaphitpc[aniter]);	
				}
			}
		}
	}
	return Manager;
}												 
