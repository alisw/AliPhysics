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
	int multbins[numOfMultBins+1] = {0, 100};
	
	int runch[numOfChTypes] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0}; // 1 - wlacza czastki do analizy
	const char *chrgs[numOfChTypes] = { "PP", "aPaP", "PaP", "KpKp", "KmKm", "KpKm", "PIpPIp", "PImPIm", "PIpPIm", "V0LL", "V0ALAL", "V0LAL", "all", "plus", "minus", "mixed"};
	

	int runqinv = 1;


	int runtype = 0; // Types 0 - global, 1 - ITS only, 2 - TPC Inner	//global tracks ->mfit ITS+TPC
	int owncuts = 0; 
	int owndca = 0;

	int gammacut = 1;	// cut na ee z gamma 
	
	double shqmax = 2.0; 
	int nbinssh = 20;

	//AliFemtoEventReaderESDChain *Reader = new AliFemtoEventReaderESDChain();
	//Reader->SetUseMultiplicity(AliFemtoEventReaderESDChain::kGlobalCount);


	AliFemtoEventReaderAODChain *Reader = new AliFemtoEventReaderAODChain();
	Reader->SetFilterMask(filterbit);
	Reader->SetDCAglobalTrack(ifGlobalTracks); //false for FB7, true for the rest //we do not use DCA at all
	Reader->SetUseMultiplicity(AliFemtoEventReaderAOD::kCentrality);
	Reader->SetMinPlpContribSPD(minPlpContribSPD);
	Reader->SetIsPileUpEvent(ifIsPileUp);
	Reader->SetReadV0(kTRUE);
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
	AliFemtoV0TrackCut              *dtc4etaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoV0TrackCut              *dtc5etaphitpc[numOfMultBins*numOfChTypes];
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
	AliFemtoCutMonitorV0            *cutPass1V0[numOfMultBins*numOfChTypes];
	AliFemtoCutMonitorV0            *cutFail1V0[numOfMultBins*numOfChTypes];
	AliFemtoCutMonitorV0            *cutPass2V0[numOfMultBins*numOfChTypes];
	AliFemtoCutMonitorV0            *cutFail2V0[numOfMultBins*numOfChTypes];
	//	 AliFemtoShareQualityTPCEntranceSepPairCut			*sqpcetaphitpcsame[numOfMultBins*numOfChTypes];
	AliFemtoPairCutAntiGamma	*sqpcetaphitpc[numOfMultBins*numOfChTypes];
	//AliFemtoPairCutRadialDistance			*sqpcetaphitpc[numOfMultBins*numOfChTypes];
	//AliFemtoShareQualityPairCut			*sqpcetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoV0PairCut               *sqp1cetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoV0TrackPairCut          *sqp2cetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoV0TrackPairCut          *sqp3cetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoV0TrackPairCut          *sqp4cetaphitpc[numOfMultBins*numOfChTypes];
	//	AliFemtoChi2CorrFctn					*cchiqinvetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoPairCutPt               *ktpcuts[numOfMultBins*numOfChTypes*numOfkTbins];
	AliFemtoQinvCorrFctn		*cqinvkttpc[numOfMultBins*numOfChTypes*numOfkTbins];
	AliFemtoQinvCorrFctn		*cqinvtpc[numOfMultBins*numOfChTypes];
	//AliFemtoCorrFctnDEtaDPhiSimpleWithCorrections *cdedpetaphi[numOfMultBins*numOfChTypes];
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
					anetaphitpc[aniter] = new AliFemtoVertexMultAnalysis(zvertbino, -10.0, 10.0, multbino, multbins[imult], multbins[imult+1]);
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

					//*********V0 cuts********************
//czastka rozpadajaca sie na plus i minus: K0 short; Lamdba (proton i pion) i anty 
					//V0 first particle cut
					dtc4etaphitpc[aniter] = new AliFemtoV0TrackCut();
					dtc4etaphitpc[aniter]->SetMass(LambdaMass);
					dtc4etaphitpc[aniter]->SetEta(0.8); //0.8					
					dtc4etaphitpc[aniter]->SetPt(0.7,maxPt);
					dtc4etaphitpc[aniter]->SetEtaDaughters(0.8);
					dtc4etaphitpc[aniter]->SetPtPosDaughter(0.5,4.0);//pos - dodatnia c?rka :) - od 0.5 bo ponizej duzo z materialu
					dtc4etaphitpc[aniter]->SetPtNegDaughter(0.16,4.0);// pionow jest duzo mniej niz z materialu
					dtc4etaphitpc[aniter]->SetTPCnclsDaughters(70);// ilosc clustrow - punkty ktore zostaly uzyte do rekonstrukcji trackow
					dtc4etaphitpc[aniter]->SetNdofDaughters(4.0); //4.0 - chi^2/ndf jakosc fitu tracku
					dtc4etaphitpc[aniter]->SetStatusDaughters(AliESDtrack::kTPCrefit/* | AliESDtrack::kITSrefit*/);
					dtc4etaphitpc[aniter]->SetOnFlyStatus(kFALSE); //kTRUE - algorytm szukania V0 - drugi offFly - wlaczony
					dtc4etaphitpc[aniter]->SetParticleType(0);//lambda
					dtc4etaphitpc[aniter]->SetMaxDcaV0Daughters(1.0); //0.5 cm
					//dtc4etaphitpc[aniter]->SetMaxDcaV0(1.0); //0.5 - odleglosc danego V0 od primary vertex
					dtc4etaphitpc[aniter]->SetMinDaughtersToPrimVertex(0.06, 0.06); //0.05 - odleglosc od Primary vertex
					dtc4etaphitpc[aniter]->SetMaxCosPointingAngle(0.99); //0.9993 
					//dtc4etaphitpc[aniter]->SetMaxV0DecayLength(100.0); //60 - odleglosc od PV do secondary vertex
					dtc4etaphitpc[aniter]->SetInvariantMassLambda(LambdaMass-0.0038,LambdaMass+0.0038);
					dtc4etaphitpc[aniter]->SetInvariantMassRejectK0s(0.48,0.515); //odrzucamy czastki majace masy blisko kaonow
					dtc4etaphitpc[aniter]->SetRadiusV0Min(0.5); //minimalna odleglosc od primary vertex do secondary
					dtc4etaphitpc[aniter]->SetNsigmaPosDaughter(5.0);
					dtc4etaphitpc[aniter]->SetNsigmaNegDaughter(5.0);// dla TPC
					dtc4etaphitpc[aniter]->SetRequireTOFPion(false);//czy uzywamy TOF dla identyfikacji c?rek
					dtc4etaphitpc[aniter]->SetRequireTOFProton(false);
				      
					//V0 second particle cut
					dtc5etaphitpc[aniter] = new AliFemtoV0TrackCut();
					dtc5etaphitpc[aniter]->SetMass(LambdaMass);
					dtc5etaphitpc[aniter]->SetEta(0.8);
					dtc5etaphitpc[aniter]->SetPt(0.7,maxPt);
					dtc5etaphitpc[aniter]->SetEtaDaughters(0.8);
					dtc5etaphitpc[aniter]->SetPtPosDaughter(0.16,4.0);
					dtc5etaphitpc[aniter]->SetPtNegDaughter(0.3,4.0);
					dtc5etaphitpc[aniter]->SetTPCnclsDaughters(70);
					dtc5etaphitpc[aniter]->SetNdofDaughters(4.0); //4.0
					dtc5etaphitpc[aniter]->SetStatusDaughters(AliESDtrack::kTPCrefit/* | AliESDtrack::kITSrefit*/);
					dtc5etaphitpc[aniter]->SetOnFlyStatus(kFALSE); //kTRUE
					dtc5etaphitpc[aniter]->SetParticleType(1);
					dtc5etaphitpc[aniter]->SetMaxDcaV0Daughters(1.0); //0.5
					//dtc5etaphitpc[aniter]->SetMaxDcaV0(1.0); //0.5
					dtc5etaphitpc[aniter]->SetMinDaughtersToPrimVertex(0.06, 0.06); //0.05
					dtc5etaphitpc[aniter]->SetMaxCosPointingAngle(0.99); //0.9993
					//dtc5etaphitpc[aniter]->SetMaxV0DecayLength(100.0); //60
					dtc5etaphitpc[aniter]->SetInvariantMassLambda(LambdaMass-0.0038,LambdaMass+0.0038);
					dtc5etaphitpc[aniter]->SetInvariantMassRejectK0s(0.48,0.515);
					dtc5etaphitpc[aniter]->SetRadiusV0Min(0.5);
					dtc5etaphitpc[aniter]->SetNsigmaPosDaughter(5.0);
					dtc5etaphitpc[aniter]->SetNsigmaNegDaughter(5.0);
					dtc5etaphitpc[aniter]->SetRequireTOFPion(false);
					dtc5etaphitpc[aniter]->SetRequireTOFProton(false);

					
					

					//****** DCA ******
					/*
					if(owndca){
					  dtc1etaphitpc[aniter]->SetMaxImpactXYPtDep(0.018, 0.035, -1.01); 	//	DCA xy - parametry zalezno?ci DCA od pt
					  //dtc1etaphitpc[aniter]->SetMaxImpactXYPtDep(0.0182, 0.0350, -1.01);
					  dtc1etaphitpc[aniter]->SetMaxImpactZ(2);	//DCA Z
					  dtc2etaphitpc[aniter]->SetMaxImpactXYPtDep(0.018, 0.035, -1.01); 	//	DCA xy
					  dtc2etaphitpc[aniter]->SetMaxImpactZ(2);	//DCA Z
					  if (ichg == 9){dtc3etaphitpc[aniter]->SetMaxImpactXYPtDep(0.018, 0.035, -1.01); 	//	DCA xy
					  dtc3etaphitpc[aniter]->SetMaxImpactZ(2);}	//DCA Z - cm
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
					  dtc2etaphitpc[aniter]->SetRemoveKinks(kTRUE); //wyrzucamy tracki "?amane"
					  dtc2etaphitpc[aniter]->SetLabel(kFALSE); //Label - numer czastki
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
					*/
					//**************** track Monitors ***************

					if(ifMonitors)//ichg>8)
					  {

					    if(ifV0Monitors){
					      // //V0 monitors (memory leak problems?)
					      cutPass1V0[aniter] = new AliFemtoCutMonitorV0(Form("cutPass1%stpcM%i", chrgs[ichg], imult));
					      cutFail1V0[aniter] = new AliFemtoCutMonitorV0(Form("cutFail1%stpcM%i", chrgs[ichg], imult));
					      dtc4etaphitpc[aniter]->AddCutMonitor(cutPass1V0[aniter], cutFail1V0[aniter]);
	  
					      cutPass2V0[aniter] = new AliFemtoCutMonitorV0(Form("cutPass2%stpcM%i", chrgs[ichg], imult));
					      cutFail2V0[aniter] = new AliFemtoCutMonitorV0(Form("cutFail2%stpcM%i", chrgs[ichg], imult));
					      dtc5etaphitpc[aniter]->AddCutMonitor(cutPass2V0[aniter], cutFail2V0[aniter]);

					    
				
					    }

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
					sqpcetaphitpc[aniter] = new AliFemtoPairCutAntiGamma();// dziedziczy po ShareQualityPairCut
					//sqpcetaphitpc[aniter] = new AliFemtoShareQualityPairCut();
					sqpcetaphitpc[aniter]->SetShareQualityMax(shareQuality);	// two track cuts on splitting and merging  //1- wylaczany 0 -wlaczany   
					sqpcetaphitpc[aniter]->SetShareFractionMax(shareFraction);	//  ile moga miec wspolnych klastrow //1 - wylaczany, 0.05 - wlaczany
					sqpcetaphitpc[aniter]->SetRemoveSameLabel(kFALSE);
				
					
					if (gammacut == 0)
					  {
					    sqpcetaphitpc[aniter]->SetMaxEEMinv(0.0);//Masa inv e+ e-
					    sqpcetaphitpc[aniter]->SetMaxThetaDiff(0.0);
					  }
					else if (gammacut == 1)
					  {
					    sqpcetaphitpc[aniter]->SetMaxEEMinv(0.002);
					    sqpcetaphitpc[aniter]->SetMaxThetaDiff(0.008); //czy leca w podobnym kierunku
					  }				     
	
					  
					//V0 two-track cuts
			
					
					sqp1cetaphitpc[aniter] = new AliFemtoV0PairCut();
					sqp1cetaphitpc[aniter]->SetDataType(AliFemtoPairCut::kAOD);//wazne z czego czytamy
					sqp1cetaphitpc[aniter]->SetTPCEntranceSepMinimum(0.00001);
					sqp1cetaphitpc[aniter]->SetTPCExitSepMinimum(-1.);
					sqp1cetaphitpc[aniter]->SetMinAvgSeparation(0,3);//odleglosci c?rek
					sqp1cetaphitpc[aniter]->SetMinAvgSeparation(1,0);
					sqp1cetaphitpc[aniter]->SetMinAvgSeparation(2,0);
					sqp1cetaphitpc[aniter]->SetMinAvgSeparation(3,3);


					
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
					
					if(ichg == 9) //V0LL
					  {
					    //anetaphitpc[aniter]->SetEventCut(mecetaphitpc[aniter]);
					    anetaphitpc[aniter]->SetFirstParticleCut(dtc4etaphitpc[aniter]);
					    anetaphitpc[aniter]->SetSecondParticleCut(dtc4etaphitpc[aniter]);
					    anetaphitpc[aniter]->SetPairCut(sqp1cetaphitpc[aniter]);
					    //avgsepcorr[aniter]->SetPairType(AliFemtoAvgSepCorrFctn::kV0s);
					  }
					if(ichg == 10) //V0ALAL
					  {
					    anetaphitpc[aniter]->SetV0SharedDaughterCut(performSharedDaughterCut); //???
					    //anetaphitpc[aniter]->SetEventCut(mecetaphitpc[aniter]);
					    anetaphitpc[aniter]->SetFirstParticleCut(dtc5etaphitpc[aniter]);
					    anetaphitpc[aniter]->SetSecondParticleCut(dtc5etaphitpc[aniter]);
					    anetaphitpc[aniter]->SetPairCut(sqp1cetaphitpc[aniter]);
					    //avgsepcorr[aniter]->SetPairType(AliFemtoAvgSepCorrFctn::kV0s);
					  }
					if(ichg == 11) //VOLAL
					  {
					    anetaphitpc[aniter]->SetV0SharedDaughterCut(performSharedDaughterCut);
					    //anetaphitpc[aniter]->SetEventCut(mecetaphitpc[aniter]);
					    anetaphitpc[aniter]->SetFirstParticleCut(dtc4etaphitpc[aniter]);
					    anetaphitpc[aniter]->SetSecondParticleCut(dtc5etaphitpc[aniter]);
					    anetaphitpc[aniter]->SetPairCut(sqp1cetaphitpc[aniter]);
					    //avgsepcorr[aniter]->SetPairType(AliFemtoAvgSepCorrFctn::kV0s);
					  }


					 //**** Correlation functions *******
					//cqinvtpc[aniter] = new AliFemtoQinvCorrFctn(Form("Qinv_%s_M%i", chrgs[ichg], imult), nbinssh, 0.0, shqmax);
					//anetaphitpc[aniter]->AddCorrFctn(cqinvtpc[aniter]);

					//***without corrections*****
					cdedpetaphinocorr[aniter] = new AliFemtoCorrFctnDEtaDPhiSimple(Form("cdedpnocorr%stpcM%i", chrgs[ichg], imult),35, 35);
					anetaphitpc[aniter]->AddCorrFctn(cdedpetaphinocorr[aniter]);
										

				

					
		
					Manager->AddAnalysis(anetaphitpc[aniter]);
				}
			}
		}
	}
					    
	return Manager;
}								 
