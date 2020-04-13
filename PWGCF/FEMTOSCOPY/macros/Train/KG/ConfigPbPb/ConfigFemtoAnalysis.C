/**********************************************************************																							 *
 * Configfemtoanalysis.C - configuration macro for the femtoscopic	  *
 * analysis, meant as a QA process for two-particle effects			  *
 *																	  *
 * Author: Adam Kisiel (Adam.Kisiel@cern.ch)						  *
 *																	  *
 *********************************************************************/
#if !defined(__CINT__) || defined(__MAKECINT_)
#include "AliFemtoManager.h"
#include "AliFemtoEventReaderESDChain.h"
#include "AliFemtoEventReaderESDChainKine.h"
#include "AliFemtoEventReaderAODMultSelection.h"
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
#include "AliFemtoCorrFctnDEtaDPhi.h"
#include "AliFemtoCorrFctnDEtaDPhiSimpleWithCorrections.h"
#include "AliFemtoCorrFctnDYDPhiSimple.h"
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
#include "AliFemtoCorrFctnDPhiStarDEta.h"

#include "AliFemtoV0PairCut.h"
#include "AliFemtoV0TrackPairCut.h"
#include "AliFemtoV0TrackCut.h"
#endif

//________________________________________________________________________
AliFemtoManager* ConfigFemtoAnalysis(const char* params) {

	double PionMass = 0.13956995;
	double KaonMass = 0.493677;
	double ProtonMass = 0.938272013;
	double LambdaMass = 1.115683;
	
	const int numOfMultBins = 1;	
	const int numOfChTypes = 16; //13
	const int numOfkTbins = 1;

	bool performSharedDaughterCut = false;
	bool enablePairMonitors = true;

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
	int filterbit = atoi(parameter[0]); //128  (768 / 128) 
	int runktdep = atoi(parameter[1]); //0
	int runmultdep = atoi(parameter[2]); //0
	int minPlpContribSPD = atoi(parameter[3]); //3
	int multbino = atoi(parameter[4]); //200
	int zvertbino = atoi(parameter[5]); //10
	Bool_t ifGlobalTracks=kFALSE; if(atoi(parameter[6]))ifGlobalTracks=kTRUE;//kFALSE 
	double shareQuality = atof(parameter[7]); //0.00
	double shareFraction = atof(parameter[8]); //0.05
	bool ifElectronRejection = atoi(parameter[9]); //true
	double nSigmaVal = atof(parameter[10]); //2.0
	double nEtaMin = atof(parameter[11]); //-0.8
	double nEtaMax = atof(parameter[12]);  //0.8
	bool ifIsPileUp = atoi(parameter[13]); //true
	double maxPt = atof(parameter[14]);  //4.0

	int setMostProb1 = atoi(parameter[15]);
	int setMostProb2 = atoi(parameter[16]);
	int setMostProb3 = atoi(parameter[17]);

	Bool_t ifMonitors=kFALSE; if(atoi(parameter[19]))ifMonitors=kTRUE;//kTRUE 
	Bool_t ifV0Monitors = kFALSE;//TRUE;
	double nSigmaVal2 = atof(parameter[20]); //3.0 (or 2.0)

	printf("*** Connect to AliEn ***\n");
	TGrid::Connect("alien://");

	int runmults[numOfMultBins] = {1};
	int multbins[numOfMultBins+1] = {0, 200};
	
	int runch[numOfChTypes] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0}; // 1 - wlacza czastki do analizy
	const char *chrgs[numOfChTypes] = { "PP", "aPaP", "PaP", "KpKp", "KmKm", "KpKm", "PIpPIp", "PImPIm", "PIpPIm", "V0LL", "V0ALAL", "V0LAL", "all", "plus", "minus", "mixed"};
	

	int runqinv = 1;


	int runtype = 0; // Types 0 - global, 1 - ITS only, 2 - TPC Inner	//global tracks ->mfit ITS+TPC
	int owncuts = 0; 
	int owndca = 0;

	int gammacut = 0;	// cut na ee z gamma 
	
	double shqmax = 2.0; 
	int nbinssh = 20;

	//AliFemtoEventReaderESDChain *Reader = new AliFemtoEventReaderESDChain();
	//Reader->SetUseMultiplicity(AliFemtoEventReaderESDChain::kGlobalCount);


	AliFemtoEventReaderAODMultSelection *Reader = new AliFemtoEventReaderAODMultSelection();
	Reader->SetFilterMask(filterbit);
	Reader->SetDCAglobalTrack(ifGlobalTracks); //false for FB7, true for the rest //we do not use DCA at all
	Reader->SetUseMultiplicity(AliFemtoEventReaderAOD::kCentrality);
	Reader->SetMinPlpContribSPD(minPlpContribSPD);
	Reader->SetIsPileUpEvent(ifIsPileUp);
	Reader->SetReadV0(kFALSE);
	//Reader->SetReadMC(kTRUE);
	
	AliFemtoManager* Manager = new AliFemtoManager();
	Manager->SetEventReader(Reader);

	AliFemtoVertexMultAnalysis	*anetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoBasicEventCut		*mecetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoCutMonitorEventMult	*cutPassEvMetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoCutMonitorEventMult	*cutFailEvMetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoCutMonitorEventVertex   *cutPassEvVetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoCutMonitorEventVertex   *cutFailEvVetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoMJTrackCut		*dtc1etaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoMJTrackCut		*dtc2etaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoMJTrackCut		*dtc3etaphitpc[numOfMultBins*numOfChTypes];
	
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
	
	//	 AliFemtoShareQualityTPCEntranceSepPairCut			*sqpcetaphitpcsame[numOfMultBins*numOfChTypes];
	// AliFemtoPairCutAntiGamma	*sqpcetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoPairCutRadialDistance			*sqpcetaphitpc[numOfMultBins*numOfChTypes];
	//AliFemtoShareQualityPairCut			*sqpcetaphitpc[numOfMultBins*numOfChTypes];
	
	//	AliFemtoChi2CorrFctn					*cchiqinvetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoPairCutPt               *ktpcuts[numOfMultBins*numOfChTypes*numOfkTbins];
	AliFemtoQinvCorrFctn		*cqinvkttpc[numOfMultBins*numOfChTypes*numOfkTbins];
	AliFemtoQinvCorrFctn		*cqinvtpc[numOfMultBins*numOfChTypes];
	AliFemtoCorrFctnDEtaDPhiSimpleWithCorrections *cdedpetaphi[numOfMultBins*numOfChTypes];
	//AliFemtoCorrFctnDEtaDPhi         *cdedpetaphi[numOfMultBins*numOfChTypes];
	AliFemtoCorrFctnDEtaDPhiSimple	*cdedpetaphinocorr[numOfMultBins*numOfChTypes];

	//	AliFemtoCorrFctnDYDPhiSimpleWithCorrections	*cdydpyphinocorr[numOfMultBins*numOfChTypes];
	AliFemtoCorrFctnDEtaDPhiCorrections *cdedpetaphiPt[numOfMultBins*numOfChTypes*numOfkTbins];

	AliFemtoCorrFctnDPhiStarDEta  *cdphistardeta08[numOfMultBins**numOfkTbins];
	AliFemtoCorrFctnDPhiStarDEta  *cdphistardeta12[numOfMultBins**numOfkTbins];
	AliFemtoCorrFctnDPhiStarDEta  *cdphistardeta16[numOfMultBins**numOfkTbins];
	AliFemtoCorrFctnDPhiStarDEta  *cdphistardeta20[numOfMultBins**numOfkTbins];
	
        AliFemtoCorrFctnNonIdDR         *cnonidtpc[numOfMultBins*numOfChTypes];

	

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
					anetaphitpc[aniter] = new AliFemtoVertexMultAnalysis(7, -7.0, 7.0, multbino, multbins[imult], multbins[imult+1]);
					anetaphitpc[aniter]->SetNumEventsToMix(5);
					anetaphitpc[aniter]->SetMinSizePartCollection(1);
					anetaphitpc[aniter]->SetVerboseMode(kFALSE);//~~~~~~~~~~~~~~~~

					//*** Event cut ***
					mecetaphitpc[aniter] = new AliFemtoBasicEventCut();
					mecetaphitpc[aniter]->SetEventMult(0.001,100000);
					mecetaphitpc[aniter]->SetVertZPos(-10,10);//cm

					//****** event monitors **********	
					cutPassEvMetaphitpc[aniter] = new AliFemtoCutMonitorEventMult(Form("cutPass%stpcM%i", chrgs[ichg], imult), 2000, 20000.5);
					cutFailEvMetaphitpc[aniter] = new AliFemtoCutMonitorEventMult(Form("cutFail%stpcM%i", chrgs[ichg], imult), 2000, 20000.5);
					mecetaphitpc[aniter]->AddCutMonitor(cutPassEvMetaphitpc[aniter], cutFailEvMetaphitpc[aniter]);

					/*
					cutPassEvVetaphitpc[aniter] = new AliFemtoCutMonitorEventVertex(Form("cutPass%stpcM%i", chrgs[ichg], imult));
					cutFailEvVetaphitpc[aniter] = new AliFemtoCutMonitorEventVertex(Form("cutFail%stpcM%i", chrgs[ichg], imult));
					mecetaphitpc[aniter]->AddCutMonitor(cutPassEvVetaphitpc[aniter], cutFailEvVetaphitpc[aniter]);
					*/

					// ***** single particle track cuts *********
					dtc1etaphitpc[aniter] = new AliFemtoMJTrackCut();
					dtc1etaphitpc[aniter]->SetCharge(1.0);
					dtc1etaphitpc[aniter]->SetEta(nEtaMin,nEtaMax);
					dtc1etaphitpc[aniter]->SetNsigma(nSigmaVal); //w zaleznosci od pedow jest brane albo nsigma1 albo nsigma2
					dtc1etaphitpc[aniter]->SetNsigma2(nSigmaVal2);

					dtc1etaphitpc[aniter]->SetNsigmaTPCTOF(kTRUE);//bierzemy jednoczesnie sigma z TPC i TOF jednoczesnie
					dtc1etaphitpc[aniter]->SetElectronRejection(ifElectronRejection);
					if (ichg == 0 || ichg == 1 ||ichg == 2) //protons
					  {
					    dtc1etaphitpc[aniter]->SetPt(0.5,maxPt); 
					    dtc1etaphitpc[aniter]->SetMass(ProtonMass);		
					    dtc1etaphitpc[aniter]->SetMostProbable(setMostProb1);//cut on Nsigma in pT not p
					    //czy dana czastka jest P, K czy PI - prawdopodobienstwo, krzywe Bettego-Blocha. 
					    //zestaw ustawien dla PID, jezeli dana czastka wpada w nsigma dla wiecej niz jednego typu to odrzucamy ja z analizy
					  }
					if (ichg == 3 ||ichg == 4 ||ichg == 5)//kaons 3-5
						{
					    dtc1etaphitpc[aniter]->SetPt(0.3,maxPt);
					    dtc1etaphitpc[aniter]->SetMass(KaonMass);
					    dtc1etaphitpc[aniter]->SetMostProbable(setMostProb2);//cut on Nsigma in pT not p
						}
					if (ichg == 6 ||ichg == 7 ||ichg == 8)//pions 6-8
						{
					    dtc1etaphitpc[aniter]->SetPt(0.2,maxPt);
					    dtc1etaphitpc[aniter]->SetMass(PionMass);		
					    dtc1etaphitpc[aniter]->SetMostProbable(setMostProb3);//cut on Nsigma in pT not p - PID - rodzaj metody		   
						}
					if (ichg == 13 ||ichg == 14 ||ichg == 15)//plus,minus,mixed
						{
					    dtc1etaphitpc[aniter]->SetPt(0.2,maxPt);
						}
	      
					dtc2etaphitpc[aniter] = new AliFemtoMJTrackCut();
					dtc2etaphitpc[aniter]->SetCharge(-1.0);
					dtc2etaphitpc[aniter]->SetEta(nEtaMin,nEtaMax);
					dtc2etaphitpc[aniter]->SetNsigma(nSigmaVal); //ustawia PID
					dtc2etaphitpc[aniter]->SetNsigma2(nSigmaVal2);//
					dtc2etaphitpc[aniter]->SetNsigmaTPCTOF(kTRUE);//
					dtc2etaphitpc[aniter]->SetElectronRejection(ifElectronRejection);
					if (ichg == 0 || ichg == 1 ||ichg == 2) //protons
					  {
					    dtc2etaphitpc[aniter]->SetPt(0.5,maxPt); 
					    dtc2etaphitpc[aniter]->SetMass(ProtonMass);		
					    dtc2etaphitpc[aniter]->SetMostProbable(setMostProb1);//cut on Nsigma in pT not p
					  }
					  
					if (ichg == 3 ||ichg == 4 ||ichg == 5)//kaons 3-5
					  {
					    dtc2etaphitpc[aniter]->SetPt(0.3,maxPt);
					    dtc2etaphitpc[aniter]->SetMass(KaonMass);
					    dtc2etaphitpc[aniter]->SetMostProbable(setMostProb2);//cut on Nsigma in pT not p
					  }
					if (ichg == 6 ||ichg == 7 ||ichg == 8)//pions 6-8
					  {
					    dtc2etaphitpc[aniter]->SetPt(0.2,maxPt);
					    dtc2etaphitpc[aniter]->SetMass(PionMass);		
					    dtc2etaphitpc[aniter]->SetMostProbable(setMostProb3);//cut on Nsigma in pT not p			   
					  }
					if (ichg == 13 ||ichg == 14 ||ichg == 15)//plus,minus,mixed
					{
					    dtc2etaphitpc[aniter]->SetPt(0.2,maxPt);
					}

					if (ichg == 12)//all
					{

					    dtc3etaphitpc[aniter] = new AliFemtoMJTrackCut();
					    dtc3etaphitpc[aniter]->SetNsigmaTPCTOF(kTRUE);
					    dtc3etaphitpc[aniter]->SetEta(nEtaMin,nEtaMax);
					    dtc3etaphitpc[aniter]->SetElectronRejection(ifElectronRejection); 
						dtc3etaphitpc[aniter]->SetPt(0.2,maxPt);
					}

					
					

					
					//**************** track Monitors ***************

					if(ifMonitors)//ichg>8)
					  {

					   

					    //FULL
						if(ichg<2 || ichg==3 || ichg==4 || ichg==6 || ichg==7 /*|| ichg == 9 || ichg == 10*/ || ichg==12 || ichg==13 || ichg==14){ 
					      cutPass3YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutPass%stpcM%i", chrgs[ichg], imult),PionMass);
					      cutFail3YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutFail%stpcM%i", chrgs[ichg], imult),PionMass);
					    }
					    if(ichg==12) dtc3etaphitpc[aniter]->AddCutMonitor(cutPass3YPtetaphitpc[aniter], cutFail3YPtetaphitpc[aniter]);
					    if(ichg==0||ichg==3||ichg==6||ichg==13)  dtc1etaphitpc[aniter]->AddCutMonitor(cutPass3YPtetaphitpc[aniter], cutFail3YPtetaphitpc[aniter]);
					    if(ichg==1||ichg==4||ichg==7||ichg==14) dtc2etaphitpc[aniter]->AddCutMonitor(cutPass3YPtetaphitpc[aniter], cutFail3YPtetaphitpc[aniter]);
			
						
						
//monitor rysuje tylko histogramy!!!
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
					  
					    
					    else if(ichg==6||ichg==13||ichg==7||ichg==14||ichg==12){ //pions, all
					      cutPass3PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cutPass%stpcM%i", chrgs[ichg], imult),0);//0-pion,1-kaon,2-proton
					      cutFail3PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cutFail%stpcM%i", chrgs[ichg], imult),0);
					      if(ichg==6||ichg==13) dtc1etaphitpc[aniter]->AddCutMonitor(cutPass3PIDetaphitpc[aniter], cutFail3PIDetaphitpc[aniter]);
					      if(ichg==7||ichg==14) dtc2etaphitpc[aniter]->AddCutMonitor(cutPass3PIDetaphitpc[aniter], cutFail3PIDetaphitpc[aniter]);
					      if(ichg==12) dtc3etaphitpc[aniter]->AddCutMonitor(cutPass3PIDetaphitpc[aniter], cutFail3PIDetaphitpc[aniter]);
					    }



					  }
					 
					//******** Two - track cuts ************
					sqpcetaphitpc[aniter] = new AliFemtoPairCutRadialDistance();
					sqpcetaphitpc[aniter]->SetShareQualityMax(shareQuality);
					sqpcetaphitpc[aniter]->SetShareFractionMax(shareFraction);
					sqpcetaphitpc[aniter]->SetRemoveSameLabel(kFALSE);
					sqpcetaphitpc[aniter]->SetMinimumRadius(1.6);
					sqpcetaphitpc[aniter]->SetEtaDifferenceMinimum(0.02);
					sqpcetaphitpc[aniter]->SetPhiStarDifferenceMinimum(0.045);


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
					// sqpcetaphitpc[aniter]->SetTPCEntranceSepMinimum(1.5);
					// sqpcetaphitpc[aniter]->SetRadialDistanceMinimum(0.12, 0.03);
					// sqpcetaphitpc[aniter]->SetEtaDifferenceMinimum(0.02);

					

					
					//***** Setting cuts ***********
					// setting event cut
					anetaphitpc[aniter]->SetEventCut(mecetaphitpc[aniter]);
					//setting single track cuts
					if(ichg==0 || ichg==3  || ichg==6 ||ichg==13) //positive like-sign
					{
					  anetaphitpc[aniter]->SetPairCut(sqpcetaphitpc[aniter]);
					  anetaphitpc[aniter]->SetFirstParticleCut(dtc1etaphitpc[aniter]);
					  anetaphitpc[aniter]->SetSecondParticleCut(dtc1etaphitpc[aniter]);
					}
					if(ichg==1 || ichg==4 || ichg==7 || ichg==14)//negative like-sign
					{
					  anetaphitpc[aniter]->SetPairCut(sqpcetaphitpc[aniter]);
					  anetaphitpc[aniter]->SetFirstParticleCut(dtc2etaphitpc[aniter]);
					  anetaphitpc[aniter]->SetSecondParticleCut(dtc2etaphitpc[aniter]);
					}
					if(ichg==2 || ichg==5 || ichg==8 || ichg==15)//unlike-sign
					{
					  anetaphitpc[aniter]->SetPairCut(sqpcetaphitpc[aniter]);
					  anetaphitpc[aniter]->SetFirstParticleCut(dtc1etaphitpc[aniter]);
					  anetaphitpc[aniter]->SetSecondParticleCut(dtc2etaphitpc[aniter]);
					}
					if(ichg==12) //all
                    {
					  anetaphitpc[aniter]->SetPairCut(sqpcetaphitpc[aniter]);
					  anetaphitpc[aniter]->SetFirstParticleCut(dtc3etaphitpc[aniter]);
					  anetaphitpc[aniter]->SetSecondParticleCut(dtc3etaphitpc[aniter]);
					}
					
					


				//**** Correlation functions *******	
					//***without corrections*****
					if(ichg >= 13)
					  cdedpetaphinocorr[aniter] = new AliFemtoCorrFctnDEtaDPhiSimple(Form("cdedpnocorr%stpcM%i", chrgs[ichg], imult),23, 23);
					else
					  cdedpetaphinocorr[aniter] = new AliFemtoCorrFctnDEtaDPhiSimple(Form("cdedpnocorr%stpcM%i", chrgs[ichg], imult),35, 35);
					anetaphitpc[aniter]->AddCorrFctn(cdedpetaphinocorr[aniter]);
				


					//***with corrections****
					if(ichg >= 13)
					  cdedpetaphi[aniter] = new AliFemtoCorrFctnDEtaDPhiSimpleWithCorrections(Form("cdedp%stpcM%i", chrgs[ichg], imult),23, 23);
					else
					  cdedpetaphi[aniter] = new AliFemtoCorrFctnDEtaDPhiSimpleWithCorrections(Form("cdedp%stpcM%i", chrgs[ichg], imult),35, 35);

					
					if(ichg==0)
					  cdedpetaphi[aniter]->SetParticleTypes(AliFemtoCorrFctnDEtaDPhiSimpleWithCorrections::kProton, AliFemtoCorrFctnDEtaDPhiSimpleWithCorrections::kProton);
					else if(ichg==1)
					  cdedpetaphi[aniter]->SetParticleTypes(AliFemtoCorrFctnDEtaDPhiSimpleWithCorrections::kProtonMinus, AliFemtoCorrFctnDEtaDPhiSimpleWithCorrections::kProtonMinus);
					else if(ichg==2)
					  cdedpetaphi[aniter]->SetParticleTypes(AliFemtoCorrFctnDEtaDPhiSimpleWithCorrections::kProton, AliFemtoCorrFctnDEtaDPhiSimpleWithCorrections::kProtonMinus);
					else if(ichg==3)
					  cdedpetaphi[aniter]->SetParticleTypes(AliFemtoCorrFctnDEtaDPhiSimpleWithCorrections::kKaon, AliFemtoCorrFctnDEtaDPhiSimpleWithCorrections::kKaon);
					else if(ichg==4)
					  cdedpetaphi[aniter]->SetParticleTypes(AliFemtoCorrFctnDEtaDPhiSimpleWithCorrections::kKaonMinus, AliFemtoCorrFctnDEtaDPhiSimpleWithCorrections::kKaonMinus);
					else if(ichg==5)
					  cdedpetaphi[aniter]->SetParticleTypes(AliFemtoCorrFctnDEtaDPhiSimpleWithCorrections::kKaon, AliFemtoCorrFctnDEtaDPhiSimpleWithCorrections::kKaonMinus);
					else if(ichg==6)
					  cdedpetaphi[aniter]->SetParticleTypes(AliFemtoCorrFctnDEtaDPhiSimpleWithCorrections::kPion, AliFemtoCorrFctnDEtaDPhiSimpleWithCorrections::kPion);
					else if(ichg==7)
					  cdedpetaphi[aniter]->SetParticleTypes(AliFemtoCorrFctnDEtaDPhiSimpleWithCorrections::kPionMinus, AliFemtoCorrFctnDEtaDPhiSimpleWithCorrections::kPionMinus);
					else if(ichg==8)
					  cdedpetaphi[aniter]->SetParticleTypes(AliFemtoCorrFctnDEtaDPhiSimpleWithCorrections::kPion, AliFemtoCorrFctnDEtaDPhiSimpleWithCorrections::kPionMinus);
					else if(ichg==10 || ichg==11 || ichg==12)
					  cdedpetaphi[aniter]->SetParticleTypes(AliFemtoCorrFctnDEtaDPhiSimpleWithCorrections::kAll, AliFemtoCorrFctnDEtaDPhiSimpleWithCorrections::kAll);					
					else if (ichg==13 || ichg==20 || ichg==27)
					  cdedpetaphi[aniter]->SetParticleTypes(AliFemtoCorrFctnDEtaDPhiSimpleWithCorrections::kLambda, AliFemtoCorrFctnDEtaDPhiSimpleWithCorrections::kProton);
					else if (ichg==14 || ichg==21 || ichg==28)
					  cdedpetaphi[aniter]->SetParticleTypes(AliFemtoCorrFctnDEtaDPhiSimpleWithCorrections::kLambdaMinus, AliFemtoCorrFctnDEtaDPhiSimpleWithCorrections::kProton);
					else if (ichg==15 || ichg==22 || ichg==29)
					  cdedpetaphi[aniter]->SetParticleTypes(AliFemtoCorrFctnDEtaDPhiSimpleWithCorrections::kLambda, AliFemtoCorrFctnDEtaDPhiSimpleWithCorrections::kProtonMinus);
					else if (ichg==16 || ichg==23 || ichg==30)
					  cdedpetaphi[aniter]->SetParticleTypes(AliFemtoCorrFctnDEtaDPhiSimpleWithCorrections::kLambdaMinus, AliFemtoCorrFctnDEtaDPhiSimpleWithCorrections::kProtonMinus);
					


					
					anetaphitpc[aniter]->AddCorrFctn(cdedpetaphi[aniter]);
						       							
					Manager->AddAnalysis(anetaphitpc[aniter]);
				}
			}
		}
	}
					    
	return Manager;
}								 
