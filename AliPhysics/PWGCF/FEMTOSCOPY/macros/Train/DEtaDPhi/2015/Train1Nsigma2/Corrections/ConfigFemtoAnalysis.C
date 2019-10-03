/********************************************************************* *
 * Configfemtoanalysis.C - configuration macro for the femtoscopic
 * analysis, meant as a QA process for two-particle effects	
 *								
 * Author: Adam Kisiel (Adam.Kisiel@cern.ch)		
 *						
 *********************************************************************/
#if !defined(__CINT__) || defined(__MAKECINT_)
#include "AliFemtoManager.h"
#include "AliFemtoEventReaderESDChain.h"
#include "AliFemtoEventReaderESDChainKine.h"
#include "AliFemtoEventReaderAODChain.h"
#include "AliFemtoSimpleAnalysis.h"
#include "AliFemtoBasicEventCut.h"
#include "AliFemtoMJTrackCut.h"
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
#include "AliFemtoCorrFctnDEtaDPhiCorrections.h"
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
AliFemtoManager* ConfigFemtoAnalysis(const char* params) {

	double PionMass = 0.13956995;
	double KaonMass = 0.493677;
	double ProtonMass = 0.938272013;
	
	const int numOfMultBins = 5;	
	const int numOfChTypes = 13;
	const int numOfkTbins = 5;

	char *parameter[21];
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
	    parameter[19] = strtok(NULL, ",");
	    cout<<"Parameter [19]: (monitors)"<<parameter[19]<<" "<<endl;
	    parameter[20] = strtok(NULL, ",");
	    cout<<"Parameter [20]: (nSigma2)"<<parameter[20]<<" "<<endl;
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
	if(fileNo==0)
	  strcpy(fileName,"alien:///alice/cern.ch/user/m/majanik/2014/DEtaDPhi/Trains/Corrections/Train7Light/1Dmap_FB96_MCOnly_DoubleCounting.root");
	else if(fileNo==1)
	  strcpy(fileName,"alien:///alice/cern.ch/user/m/majanik/2014/DEtaDPhi/Trains/Corrections/Train7Light/1Dmap_FB96_MCDCA_DoubleCounting.root");
	else if(fileNo==2)
	  strcpy(fileName,"alien:///alice/cern.ch/user/m/majanik/2014/DEtaDPhi/Trains/Corrections/Train7Light/1Dmap_FB96_MCOnly_NoDoubleCounting.root");
	else if(fileNo==3)
	  strcpy(fileName,"alien:///alice/cern.ch/user/m/majanik/2014/DEtaDPhi/Trains/Corrections/Train7Light/1Dmap_FB96_MCDCA_NoDoubleCounting.root");
	else if(fileNo==4)
	  strcpy(fileName,"alien:///alice/cern.ch/user/m/majanik/2014/DEtaDPhi/Trains/Corrections/Train7Light/1Dmap_FB96_MCOnly_Exclusive.root");
	else if(fileNo==5)
	  strcpy(fileName,"alien:///alice/cern.ch/user/m/majanik/2014/DEtaDPhi/Trains/Corrections/Train7Light/1Dmap_FB96_MCDCA_Exclusive.root");
	else if(fileNo==6)
	  strcpy(fileName,"alien:///alice/cern.ch/user/m/majanik/2014/DEtaDPhi/Trains/Corrections/Train7Light/1Dmap_FB96_MCOnly_DoubleCountingNsigma2.root");
	else if(fileNo==7)
	  strcpy(fileName,"alien:///alice/cern.ch/user/m/majanik/2014/DEtaDPhi/Trains/Corrections/Train7Light/1Dmap_FB96_MCDCA_DoubleCountingNsigma2.root");
	else if(fileNo==8)
	  strcpy(fileName,"alien:///alice/cern.ch/user/m/majanik/2014/DEtaDPhi/Trains/Corrections/Train7Light/1Dmap_FB768_MCOnly_NoDoubleCounting.root");
	else if(fileNo==9)
	  strcpy(fileName,"alien:///alice/cern.ch/user/m/majanik/2014/DEtaDPhi/Trains/Corrections/Train7Light/1Dmap_FB768_MCDCA_NoDoubleCounting.root");
	else if(fileNo==10)
	  strcpy(fileName,"alien:///alice/cern.ch/user/m/majanik/2014/DEtaDPhi/Trains/Corrections/Train7Light/1Dmap_FB768_MCOnly_DoubleCounting.root");
	else if(fileNo==11)
	  strcpy(fileName,"alien:///alice/cern.ch/user/m/majanik/2014/DEtaDPhi/Trains/Corrections/Train7Light/1Dmap_FB768_MCDCA_DoubleCounting.root");
	else if(fileNo==12)
	  strcpy(fileName,"alien:///alice/cern.ch/user/m/majanik/2014/DEtaDPhi/Trains/Corrections/Train7Light/4Dmap_FB96_MCOnly_NoDoubleCounting.root");
	else if(fileNo==13)
	  strcpy(fileName,"alien:///alice/cern.ch/user/m/majanik/2014/DEtaDPhi/Trains/Corrections/Train7Light/4Dmap_FB96_MCDCA_NoDoubleCounting.root");
	else if(fileNo==14)
	  strcpy(fileName,"alien:///alice/cern.ch/user/m/majanik/2014/DEtaDPhi/Trains/Corrections/Train7Light/4Dmap_FB768_MCOnly_DoubleCounting.root");
	else if(fileNo==15)
	  strcpy(fileName,"alien:///alice/cern.ch/user/m/majanik/2014/DEtaDPhi/Trains/Corrections/Train7Light/4Dmap_FB768_MCDCA_DoubleCounting.root");
	else if(fileNo==16)
	  strcpy(fileName,"alien:///alice/cern.ch/user/m/majanik/2015/DEtaDPhi/Trains/Corrections/CorrectionFiles/1Dmap_FB96_MCOnly_ExclusivePIDNsigmaHalf2.root");
	else if(fileNo==17)
	  strcpy(fileName,"alien:///alice/cern.ch/user/m/majanik/2015/DEtaDPhi/Trains/Corrections/CorrectionFiles/AnalysisNote/1Dmap_FB96_MCOnlyLHC10f6a_Exclusive2Rej3.root");
	else if(fileNo==18)
	  strcpy(fileName,"alien:///alice/cern.ch/user/m/majanik/2015/DEtaDPhi/Trains/Corrections/CorrectionFiles/AnalysisNote/1Dmap_FB768_MCOnlyLHC10f6a_Exclusive2Rej3.root");
	else if(fileNo==19)
	  strcpy(fileName,"alien:///alice/cern.ch/user/m/majanik/2015/DEtaDPhi/Trains/Corrections/CorrectionFiles/AnalysisNote/1Dmap_FB128_MCOnlyLHC10f6a_Exclusive2Rej3.root");
	else if(fileNo==20)
	  strcpy(fileName,"alien:///alice/cern.ch/user/m/majanik/2015/DEtaDPhi/Trains/Corrections/CorrectionFiles/AnalysisNote/1Dmap_FB96_MCOnlyLHC10f6a_Exclusive2.root");
	else if(fileNo==21)
	  strcpy(fileName,"alien:///alice/cern.ch/user/m/majanik/2015/DEtaDPhi/Trains/Corrections/CorrectionFiles/AnalysisNote/1Dmap_FB96_MCOnlyLHC10f6a_NoDouble2.root");

	cout<<"Filename: "<<Form("%s",fileName)<<endl;


	Bool_t ifMonitors=kFALSE; if(atoi(parameter[19]))ifMonitors=kTRUE;//kTRUE 
	double nSigmaVal2 = atof(parameter[20]); //2.0 or 3.0

	printf("*** Connect to AliEn ***\n");
	TGrid::Connect("alien://");

	int runmults[numOfMultBins] = {0, 0, 0, 0, 1};
	if(runmultdep)	  {runmults[0]=1; runmults[1]=1; runmults[2]=1;	  }
	int multbins[numOfMultBins+1] = {2, 20, 50,150,2,150};
	
	int runch[numOfChTypes] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1};
	const char *chrgs[numOfChTypes] = { "PP", "aPaP", "PaP", "KpKp", "KmKm", "KpKm", "PIpPIp", "PImPIm", "PIpPIm", "all", "plus", "minus", "mixed" };
	
	double ktrng[numOfkTbins+1] = {0.0, 0, 0, 0, 0, 0};
	double ktrngAll[numOfkTbins+1] = {0.0, 1.0, 2.0, 3.0, 4.0, 100.0};
	double ktrngPion[numOfkTbins+1] = {0.0, 0.8, 1.2, 1.4, 2.5, 100.0};
	double ktrngKaon[numOfkTbins+1] = {0.0, 2.5, 3.5, 100.0, 0, 0};
	double ktrngProton[numOfkTbins+1] = {0.0, 2.75, 100, 0, 0, 0};

	int runqinv = 1;
	int runshlcms = 1;// 0:PRF(PAP), 1:LCMS(PP,APAP)

	int runtype = 0; // Types 0 - global, 1 - ITS only, 2 - TPC Inner	//global tracks ->mfit ITS+TPC
	int owncuts = 0; 
	int owndca = 0;

	int gammacut = 0;	// cut na ee z gamma 
	
	double shqmax = 0.5; 
	int nbinssh = 100;

	//AliFemtoEventReaderESDChain *Reader = new AliFemtoEventReaderESDChain();
	//Reader->SetUseMultiplicity(AliFemtoEventReaderESDChain::kGlobalCount);


	AliFemtoEventReaderAODChain *Reader = new AliFemtoEventReaderAODChain();
	Reader->SetFilterMask(filterbit);
	Reader->SetDCAglobalTrack(ifGlobalTracks); //false for FB7, true for the rest //we do not use DCA at all
	Reader->SetUseMultiplicity(AliFemtoEventReaderAOD::kReference);
	Reader->SetMinPlpContribSPD(minPlpContribSPD);
	Reader->SetIsPileUpEvent(ifIsPileUp);

	AliFemtoManager* Manager = new AliFemtoManager();
	Manager->SetEventReader(Reader);

	AliFemtoVertexMultAnalysis		*anetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoBasicEventCut				 *mecetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoCutMonitorEventMult	 *cutPassEvMetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoCutMonitorEventMult	 *cutFailEvMetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoCutMonitorEventVertex *cutPassEvVetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoCutMonitorEventVertex *cutFailEvVetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoMJTrackCut			 *dtc1etaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoMJTrackCut			 *dtc2etaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoMJTrackCut			 *dtc3etaphitpc[numOfMultBins*numOfChTypes];
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
	AliFemtoShareQualityPairCut			*sqpcetaphitpc[numOfMultBins*numOfChTypes];
	//	AliFemtoChi2CorrFctn					*cchiqinvetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoPairCutPt					*ktpcuts[numOfMultBins*numOfChTypes*numOfkTbins];
	AliFemtoQinvCorrFctn					*cqinvkttpc[numOfMultBins*numOfChTypes*numOfkTbins];
	AliFemtoQinvCorrFctn					*cqinvtpc[numOfMultBins*numOfChTypes];
	AliFemtoCorrFctnDEtaDPhiCorrections			*cdedpetaphi[numOfMultBins*numOfChTypes];
	AliFemtoCorrFctnDEtaDPhiCorrections			*cdedpetaphiPt[numOfMultBins*numOfChTypes*numOfkTbins];


	
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
					dtc1etaphitpc[aniter] = new AliFemtoMJTrackCut();
					dtc2etaphitpc[aniter] = new AliFemtoMJTrackCut();
					dtc1etaphitpc[aniter]->SetNsigmaTPCTOF(kTRUE);
					dtc2etaphitpc[aniter]->SetNsigmaTPCTOF(kTRUE);
					dtc1etaphitpc[aniter]->SetNsigma(nSigmaVal);
					dtc2etaphitpc[aniter]->SetNsigma(nSigmaVal);
					dtc1etaphitpc[aniter]->SetNsigma2(nSigmaVal2);
					dtc2etaphitpc[aniter]->SetNsigma2(nSigmaVal2);
					//dtc3etaphitpc[aniter]->SetNsigma(3.0);

					dtc1etaphitpc[aniter]->SetCharge(1.0);
					dtc2etaphitpc[aniter]->SetCharge(-1.0);

					dtc1etaphitpc[aniter]->SetEta(nEtaMin,nEtaMax);
					dtc2etaphitpc[aniter]->SetEta(nEtaMin,nEtaMax);


					dtc1etaphitpc[aniter]->SetElectronRejection(ifElectronRejection);	
					dtc2etaphitpc[aniter]->SetElectronRejection(ifElectronRejection);



                                        if (ichg == 0 ||ichg == 1 ||ichg == 2)//protons 0-2
                                          {
                                            dtc1etaphitpc[aniter]->SetPt(0.5,maxPt);
                                            dtc2etaphitpc[aniter]->SetPt(0.5,maxPt);
                                           
					    dtc1etaphitpc[aniter]->SetMass(ProtonMass);		
					    dtc1etaphitpc[aniter]->SetMostProbable(setMostProb1);//cut on Nsigma in pT not p
					    dtc2etaphitpc[aniter]->SetMass(ProtonMass);		
					    dtc2etaphitpc[aniter]->SetMostProbable(setMostProb1);//cut on Nsigma in pT not p
                                          }

					if (ichg == 3 ||ichg == 4 ||ichg == 5)//kaons 3-5
                                          {
                                            dtc1etaphitpc[aniter]->SetPt(0.3,maxPt);
                                            dtc2etaphitpc[aniter]->SetPt(0.3,maxPt);
 					    dtc1etaphitpc[aniter]->SetMass(KaonMass);
					    dtc1etaphitpc[aniter]->SetMostProbable(setMostProb2);//cut on Nsigma in pT not p
					    dtc2etaphitpc[aniter]->SetMass(KaonMass);
					    dtc2etaphitpc[aniter]->SetMostProbable(setMostProb2);//cut on Nsigma in pT not p

                                          }
                                        if (ichg == 6 ||ichg == 7 ||ichg == 8)//pions 6-8
                                          {
                                            dtc1etaphitpc[aniter]->SetPt(0.2,maxPt);
					    dtc2etaphitpc[aniter]->SetPt(0.2,maxPt);

					    dtc1etaphitpc[aniter]->SetMass(PionMass);		
					    dtc1etaphitpc[aniter]->SetMostProbable(setMostProb3);//cut on Nsigma in pT not p
					    dtc2etaphitpc[aniter]->SetMass(PionMass);		
					    dtc2etaphitpc[aniter]->SetMostProbable(setMostProb3);//cut on Nsigma in pT not p
                                          }
                                        if (ichg == 9)//all
                                          {

					    dtc3etaphitpc[aniter] = new AliFemtoMJTrackCut();
					    dtc3etaphitpc[aniter]->SetNsigmaTPCTOF(kTRUE);
					    dtc3etaphitpc[aniter]->SetEta(nEtaMin,nEtaMax);
					    dtc3etaphitpc[aniter]->SetElectronRejection(ifElectronRejection); 
                                            dtc3etaphitpc[aniter]->SetPt(0.2,maxPt);
                                          }
                                        if (ichg == 10 ||ichg == 11 ||ichg == 12)//plus,minus,mixed
                                          {
                                            dtc1etaphitpc[aniter]->SetPt(0.2,maxPt);
                                            dtc2etaphitpc[aniter]->SetPt(0.2,maxPt);
                                          }

					//****** DCA ******

					if(owndca){
					  dtc1etaphitpc[aniter]->SetMaxImpactXYPtDep(0.018, 0.035, -1.01); 	//	DCA xy
					  //dtc1etaphitpc[aniter]->SetMaxImpactXYPtDep(0.0182, 0.0350, -1.01);
					  dtc1etaphitpc[aniter]->SetMaxImpactZ(2);	//DCA Z
					  dtc2etaphitpc[aniter]->SetMaxImpactXYPtDep(0.018, 0.035, -1.01); 	//	DCA xy
					  dtc2etaphitpc[aniter]->SetMaxImpactZ(2);	//DCA Z
					  if (ichg == 9){dtc3etaphitpc[aniter]->SetMaxImpactXYPtDep(0.018, 0.035, -1.01); 	//	DCA xy
					  dtc3etaphitpc[aniter]->SetMaxImpactZ(2);}	//DCA Z
					}
					//****** Track quality cuts ******

					if(owncuts){

					  dtc1etaphitpc[aniter]->SetStatus(AliESDtrack::kTPCrefit|AliESDtrack::kITSrefit);
					  dtc1etaphitpc[aniter]->SetminTPCncls(70);
					  dtc1etaphitpc[aniter]->SetRemoveKinks(kTRUE);
					  dtc1etaphitpc[aniter]->SetLabel(kFALSE);
					  //	dtc1etaphitpc[aniter]->SetMaxITSChiNdof(6.0);
					  dtc1etaphitpc[aniter]->SetMaxTPCChiNdof(4.0);	// pisac
					  //dtc1etaphitpc[aniter]->SetMaxSigmaToVertex(6.0);

					  dtc2etaphitpc[aniter]->SetStatus(AliESDtrack::kTPCrefit|AliESDtrack::kITSrefit);
					  dtc2etaphitpc[aniter]->SetminTPCncls(70);
					  dtc2etaphitpc[aniter]->SetRemoveKinks(kTRUE);
					  dtc2etaphitpc[aniter]->SetLabel(kFALSE);
					  //	dtc2etaphitpc[aniter]->SetMaxITSChiNdof(6.0);
					  dtc2etaphitpc[aniter]->SetMaxTPCChiNdof(4.0);
					  //	dtc2etaphitpc[aniter]->SetMaxSigmaToVertex(6.0);
					  if (ichg == 9){
					    dtc3etaphitpc[aniter]->SetStatus(AliESDtrack::kTPCrefit|AliESDtrack::kITSrefit);
					    dtc3etaphitpc[aniter]->SetminTPCncls(70);
					    dtc3etaphitpc[aniter]->SetRemoveKinks(kTRUE);
					    dtc3etaphitpc[aniter]->SetLabel(kFALSE);
					    //	dtc3etaphitpc[aniter]->SetMaxITSChiNdof(6.0);
					    dtc3etaphitpc[aniter]->SetMaxTPCChiNdof(4.0);
					    //	dtc3etaphitpc[aniter]->SetMaxSigmaToVertex(6.0);
					  }
	
					}
					//**************** track Monitors ***************

					if(ifMonitors)//ichg>8)
					  {
					    //FULL
					    if(ichg<2 || ichg==3||ichg==4 || ichg==6|| ichg==7||ichg==9||ichg==10||ichg==11){ 
					    //if(ichg==0 || ichg==3 || ichg==6 || ichg==10){ 
					      cutPass3YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutPass%stpcM%i", chrgs[ichg], imult),PionMass);
					      cutFail3YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutFail%stpcM%i", chrgs[ichg], imult),PionMass);
					    }
					    if(ichg==9) dtc3etaphitpc[aniter]->AddCutMonitor(cutPass3YPtetaphitpc[aniter], cutFail3YPtetaphitpc[aniter]);
					    if(ichg==0||ichg==3||ichg==6||ichg==10)  dtc1etaphitpc[aniter]->AddCutMonitor(cutPass3YPtetaphitpc[aniter], cutFail3YPtetaphitpc[aniter]);
					    if(ichg==1||ichg==4||ichg==7||ichg==11) dtc2etaphitpc[aniter]->AddCutMonitor(cutPass3YPtetaphitpc[aniter], cutFail3YPtetaphitpc[aniter]);


					    if(ichg<2){ //PP, PaP
					      cutPass3PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cutPass%stpcM%i", chrgs[ichg], imult),2);//0-pion,1-kaon,2-proton
					      cutFail3PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cutFail%stpcM%i", chrgs[ichg], imult),2);
					      if(ichg==0) dtc1etaphitpc[aniter]->AddCutMonitor(cutPass3PIDetaphitpc[aniter], cutFail3PIDetaphitpc[aniter]);
					      if(ichg==1) dtc2etaphitpc[aniter]->AddCutMonitor(cutPass3PIDetaphitpc[aniter], cutFail3PIDetaphitpc[aniter]);
		    			    }
					    else if(ichg>=3 && ichg<=4){//KpKp, KmKm
					      cutPass3PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cutPass%stpcM%i", chrgs[ichg], imult),1);//0-pion,1-kaon,2-proton
					      cutFail3PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cutFail%stpcM%i", chrgs[ichg], imult),1);
					      if(ichg==3) dtc1etaphitpc[aniter]->AddCutMonitor(cutPass3PIDetaphitpc[aniter], cutFail3PIDetaphitpc[aniter]);
					      if(ichg==4) dtc2etaphitpc[aniter]->AddCutMonitor(cutPass3PIDetaphitpc[aniter], cutFail3PIDetaphitpc[aniter]);
					    }
					    else if(ichg==6||ichg==10||ichg==7||ichg==11||ichg==9){ //pions, all
					      cutPass3PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cutPass%stpcM%i", chrgs[ichg], imult),0);//0-pion,1-kaon,2-proton
					      cutFail3PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cutFail%stpcM%i", chrgs[ichg], imult),0);
					      if(ichg==6||ichg==10) dtc1etaphitpc[aniter]->AddCutMonitor(cutPass3PIDetaphitpc[aniter], cutFail3PIDetaphitpc[aniter]);
					      if(ichg==7||ichg==11) dtc2etaphitpc[aniter]->AddCutMonitor(cutPass3PIDetaphitpc[aniter], cutFail3PIDetaphitpc[aniter]);
					      if(ichg==9) dtc3etaphitpc[aniter]->AddCutMonitor(cutPass3PIDetaphitpc[aniter], cutFail3PIDetaphitpc[aniter]);
					    }

					  }
					 
					//******** Two - track cuts ************
					//sqpcetaphitpc[aniter] = new AliFemtoPairCutAntiGamma();
					sqpcetaphitpc[aniter] = new AliFemtoShareQualityPairCut();
				
					sqpcetaphitpc[aniter]->SetShareQualityMax(shareQuality);		// two track cuts on splitting and merging  //1- wylaczany 0 -wlaczany   
					sqpcetaphitpc[aniter]->SetShareFractionMax(shareFraction);	//  ile moga miec wspolnych klastrow //1 - wylaczany, 0.05 - wlaczany
					sqpcetaphitpc[aniter]->SetRemoveSameLabel(kFALSE);
					// sqpcetaphitpc[aniter]->SetMaximumRadius(0.82);
					// sqpcetaphitpc[aniter]->SetMinimumRadius(0.8);
					// sqpcetaphitpc[aniter]->SetPhiStarDifferenceMinimum(0.02);
					// sqpcetaphitpc[aniter]->SetEtaDifferenceMinimum(0.02);
					/*
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
					if(ichg==0 || ichg==3  || ichg==6 || ichg==10) //positive like-sign
					{
						anetaphitpc[aniter]->SetFirstParticleCut(dtc1etaphitpc[aniter]);
						anetaphitpc[aniter]->SetSecondParticleCut(dtc1etaphitpc[aniter]);
					}
					if(ichg==1 || ichg==4 || ichg==7 || ichg==11)//negative like-sign
					{
						anetaphitpc[aniter]->SetFirstParticleCut(dtc2etaphitpc[aniter]);
						anetaphitpc[aniter]->SetSecondParticleCut(dtc2etaphitpc[aniter]);
					}
					if(ichg==2 || ichg==5 || ichg==8 || ichg==12)//unlike-sign
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

					cdedpetaphi[aniter] = new AliFemtoCorrFctnDEtaDPhiCorrections(Form("cdedp%stpcM%i", chrgs[ichg], imult),29, 29);
					cdedpetaphi[aniter]->SetDoFullAnalysis(kFALSE);
					if(ichg==0 || ichg==1 || ichg==2)
					  cdedpetaphi[aniter]->LoadCorrectionTabFromROOTFile1D(Form("%s",fileName), AliFemtoCorrFctnDEtaDPhiCorrections::kProton, AliFemtoCorrFctnDEtaDPhiCorrections::kProton/*,1,1,1,0*/);
					//cdedpetaphi[aniter]->LoadCorrectionTabFromROOTFile(Form("%s",fileName), AliFemtoCorrFctnDEtaDPhiCorrections::kProton, AliFemtoCorrFctnDEtaDPhiCorrections::kProton,1,1,1,1);
					else if(ichg==3 || ichg==4 || ichg==5)
					  cdedpetaphi[aniter]->LoadCorrectionTabFromROOTFile1D(Form("%s",fileName), AliFemtoCorrFctnDEtaDPhiCorrections::kKaon, AliFemtoCorrFctnDEtaDPhiCorrections::kKaon/*,1,1,1,0*/);
					  //cdedpetaphi[aniter]->LoadCorrectionTabFromROOTFile(Form("%s",fileName), AliFemtoCorrFctnDEtaDPhiCorrections::kKaon, AliFemtoCorrFctnDEtaDPhiCorrections::kKaon,1,1,1,1);	
					 
					else if(ichg==6 || ichg==7 || ichg==8)
					  cdedpetaphi[aniter]->LoadCorrectionTabFromROOTFile1D(Form("%s",fileName), AliFemtoCorrFctnDEtaDPhiCorrections::kPion, AliFemtoCorrFctnDEtaDPhiCorrections::kPion/*,1,1,1,0*/);
					  //cdedpetaphi[aniter]->LoadCorrectionTabFromROOTFile(Form("%s",fileName), AliFemtoCorrFctnDEtaDPhiCorrections::kPion, AliFemtoCorrFctnDEtaDPhiCorrections::kPion,1,1,1,1);	
					
					else
					  cdedpetaphi[aniter]->LoadCorrectionTabFromROOTFile1D(Form("%s",fileName), AliFemtoCorrFctnDEtaDPhiCorrections::kAll, AliFemtoCorrFctnDEtaDPhiCorrections::kAll/*,1,1,1,0*/);
					  //cdedpetaphi[aniter]->LoadCorrectionTabFromROOTFile(Form("%s",fileName), AliFemtoCorrFctnDEtaDPhiCorrections::kAll, AliFemtoCorrFctnDEtaDPhiCorrections::kAll,1,1,1,1);
					
			     

	     
					anetaphitpc[aniter]->AddCorrFctn(cdedpetaphi[aniter]);

		
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

							cdedpetaphiPt[ktm] = new AliFemtoCorrFctnDEtaDPhiCorrections(Form("cdedp%stpcM%ipT%i", chrgs[ichg], imult,ikt),29, 29);
							cdedpetaphiPt[ktm]->SetPairSelectionCut(ktpcuts[ktm]);

							if(ichg==0 || ichg==1 || ichg==2)
							  cdedpetaphiPt[ktm]->LoadCorrectionTabFromROOTFile1D(Form("%s",fileName), AliFemtoCorrFctnDEtaDPhiCorrections::kProton, AliFemtoCorrFctnDEtaDPhiCorrections::kProton/*,1,1,1,1*/);
							else if(ichg==3 || ichg==4 || ichg==5)
							  cdedpetaphiPt[ktm]->LoadCorrectionTabFromROOTFile1D(Form("%s",fileName), AliFemtoCorrFctnDEtaDPhiCorrections::kKaon, AliFemtoCorrFctnDEtaDPhiCorrections::kKaon);
							else if(ichg==6 || ichg==7 || ichg==8)
							  cdedpetaphiPt[ktm]->LoadCorrectionTabFromROOTFile1D(Form("%s",fileName), AliFemtoCorrFctnDEtaDPhiCorrections::kPion, AliFemtoCorrFctnDEtaDPhiCorrections::kPion);
							else
							  cdedpetaphiPt[ktm]->LoadCorrectionTabFromROOTFile1D(Form("%s",fileName), AliFemtoCorrFctnDEtaDPhiCorrections::kAll, AliFemtoCorrFctnDEtaDPhiCorrections::kAll);
				

							/*if(ichg==0 || ichg==1 || ichg==2)
							  cdedpetaphiPt[ktm]->LoadCorrectionTabFromFile("alien:///alice/cern.ch/user/l/lgraczyk/2014/DEtaDPhi/CorrectionTables/pTab.txt","alien:///alice/cern.ch/user/l/lgraczyk/2014/DEtaDPhi/CorrectionTables/protonCorrTab.txt");
							else if(ichg==3 || ichg==4 || ichg==5)
							  cdedpetaphiPt[ktm]->LoadCorrectionTabFromFile("alien:///alice/cern.ch/user/l/lgraczyk/2014/DEtaDPhi/CorrectionTables/pTab.txt","alien:///alice/cern.ch/user/l/lgraczyk/2014/DEtaDPhi/CorrectionTables/kaonCorrTab.txt");	
							else if(ichg==6 || ichg==7 || ichg==8)
							  cdedpetaphiPt[ktm]->LoadCorrectionTabFromFile("alien:///alice/cern.ch/user/l/lgraczyk/2014/DEtaDPhi/CorrectionTables/pTab.txt","alien:///alice/cern.ch/user/l/lgraczyk/2014/DEtaDPhi/CorrectionTables/pionCorrTab.txt");
							else
							cdedpetaphiPt[ktm]->LoadCorrectionTabFromFile("alien:///alice/cern.ch/user/l/lgraczyk/2014/DEtaDPhi/CorrectionTables/pTab.txt","alien:///alice/cern.ch/user/l/lgraczyk/2014/DEtaDPhi/CorrectionTables/allCorrTab.txt"); 	*/		
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
