/********************************************************************* *																							 *
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
	
	const int numOfMultBins = 5;	
	const int numOfChTypes = 34; //13
	const int numOfkTbins = 5;

	bool performSharedDaughterCut = true;
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

	Bool_t ifMonitors=kFALSE; if(atoi(parameter[19]))ifMonitors=kTRUE;//kTRUE 
	double nSigmaVal2 = atof(parameter[20]); //2.0 or 3.0

	printf("*** Connect to AliEn ***\n");
	TGrid::Connect("alien://");

	int runmults[numOfMultBins] = {0, 0, 0, 0, 1};
	if(runmultdep)	  {runmults[0]=1; runmults[1]=1; runmults[2]=1;	  }
	int multbins[numOfMultBins+1] = {2, 20, 50,150,2,150};
	
	int runch[numOfChTypes] = {/*protons*/1, 1, 1, /* kaons */ 1, 1, 1, /* pions */ 1, 1, 1, /* no PID */ 0, 0, 0, 0,/*other*/ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, /*p-lam */ 1, 1, 1, 1, /* lambdas */ 1, 1, 1};
	const char *chrgs[numOfChTypes] = { "PP", "aPaP", "PaP", "KpKp", "KmKm", "KpKm", "PIpPIp", "PImPIm", "PIpPIm", "all", "plus", "minus", "mixed", "V0PLlowPt","V0PALlowPt","V0APLlowPt","V0APALlowPt","V0LLlowPt","V0LALlowPt","V0ALALlowPt", "V0PLhighPt","V0PALhighPt","V0APLhighPt","V0APALhighPt","V0LLhighPt","V0LALhighPt","V0ALALhighPt", "V0PL","V0PAL","V0APL","V0APAL","V0LL","V0LAL","V0ALAL" };
	
	
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

	int gammacut = 1;	// cut na ee z gamma 
	
	double shqmax = 1.0; 
	int nbinssh = 100;

	//AliFemtoEventReaderESDChain *Reader = new AliFemtoEventReaderESDChain();
	//Reader->SetUseMultiplicity(AliFemtoEventReaderESDChain::kGlobalCount);


	AliFemtoEventReaderAODChain *Reader = new AliFemtoEventReaderAODChain();
	Reader->SetFilterMask(filterbit);
	Reader->SetDCAglobalTrack(ifGlobalTracks); //false for FB7, true for the rest //we do not use DCA at all
	Reader->SetUseMultiplicity(AliFemtoEventReaderAOD::kReference);
	Reader->SetMinPlpContribSPD(minPlpContribSPD);
	Reader->SetIsPileUpEvent(ifIsPileUp);
	Reader->SetReadV0(kTRUE);

	AliFemtoManager* Manager = new AliFemtoManager();
	Manager->SetEventReader(Reader);

	AliFemtoVertexMultAnalysis	*anetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoBasicEventCut		*mecetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoCutMonitorEventMult	*cutPassEvMetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoCutMonitorEventMult	*cutFailEvMetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoCutMonitorEventVertex   *cutPassEvVetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoCutMonitorEventVertex   *cutFailEvVetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoCutMonitorCollections   *cutPassColletaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoCutMonitorCollections   *cutFailColletaphitpc[numOfMultBins*numOfChTypes];
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
	AliFemtoCorrFctnDEtaDPhiSimpleWithCorrections *cdedpetaphi[numOfMultBins*numOfChTypes];
	AliFemtoCorrFctnDEtaDPhiSimple	*cdedpetaphinocorr[numOfMultBins*numOfChTypes];
	//	AliFemtoCorrFctnDYDPhiSimpleWithCorrections	*cdydpyphinocorr[numOfMultBins*numOfChTypes];
	AliFemtoCorrFctnDEtaDPhiCorrections *cdedpetaphiPt[numOfMultBins*numOfChTypes*numOfkTbins];

	AliFemtoCorrFctnDPhiStarDEta  *cdphistardeta08[numOfMultBins**numOfkTbins];
	AliFemtoCorrFctnDPhiStarDEta  *cdphistardeta12[numOfMultBins**numOfkTbins];
	AliFemtoCorrFctnDPhiStarDEta  *cdphistardeta16[numOfMultBins**numOfkTbins];
	AliFemtoCorrFctnDPhiStarDEta  *cdphistardeta20[numOfMultBins**numOfkTbins];
	
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
					mecetaphitpc[aniter]->SetVertZPos(0,4);//cm

					//****** event monitors **********	
					cutPassEvMetaphitpc[aniter] = new AliFemtoCutMonitorEventMult(Form("cutPass%stpcM%i", chrgs[ichg], imult));
					cutFailEvMetaphitpc[aniter] = new AliFemtoCutMonitorEventMult(Form("cutFail%stpcM%i", chrgs[ichg], imult));
					mecetaphitpc[aniter]->AddCutMonitor(cutPassEvMetaphitpc[aniter], cutFailEvMetaphitpc[aniter]);
		
					//Study the collection multiplicity distribution
					cutPassColletaphitpc[aniter] = new AliFemtoCutMonitorCollections(Form("cutPass%stpcM%i", chrgs[ichg], imult));
					cutFailColletaphitpc[aniter] = new AliFemtoCutMonitorCollections(Form("cutFail%stpcM%i", chrgs[ichg], imult));
					mecetaphitpc[aniter]->AddCutMonitor(cutPassColletaphitpc[aniter], cutFailColletaphitpc[aniter]);
					
					// ***** single particle track cuts *********
					dtc1etaphitpc[aniter] = new AliFemtoMJTrackCut();
					dtc1etaphitpc[aniter]->SetCharge(1.0);
					dtc1etaphitpc[aniter]->SetEta(nEtaMin,nEtaMax);
					dtc1etaphitpc[aniter]->SetNsigma(nSigmaVal);
					dtc1etaphitpc[aniter]->SetNsigma2(nSigmaVal2);
					dtc1etaphitpc[aniter]->SetNsigmaTPCTOF(kTRUE);
					dtc1etaphitpc[aniter]->SetElectronRejection(ifElectronRejection);
					if (ichg == 0 || ichg == 1 ||ichg == 2 || ichg == 13 || ichg == 14 || ichg == 15 || ichg == 16 || ichg == 20 || ichg == 21 || ichg == 22 || ichg == 23 || ichg == 27 || ichg == 28 || ichg == 29 || ichg == 30) //protons
					  {
					    dtc1etaphitpc[aniter]->SetPt(0.5,maxPt); 
					    dtc1etaphitpc[aniter]->SetMass(ProtonMass);		
					    dtc1etaphitpc[aniter]->SetMostProbable(setMostProb1);//cut on Nsigma in pT not p
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
					    dtc1etaphitpc[aniter]->SetMostProbable(setMostProb3);//cut on Nsigma in pT not p			   
					  }
                                        if (ichg == 10 ||ichg == 11 ||ichg == 12)//plus,minus,mixed
                                          {
					    dtc1etaphitpc[aniter]->SetPt(0.2,maxPt);
					  }
	      
					dtc2etaphitpc[aniter] = new AliFemtoMJTrackCut();
					dtc2etaphitpc[aniter]->SetCharge(-1.0);
					dtc2etaphitpc[aniter]->SetEta(nEtaMin,nEtaMax);
					dtc2etaphitpc[aniter]->SetNsigma(nSigmaVal);
					dtc2etaphitpc[aniter]->SetNsigma2(nSigmaVal2);
					dtc2etaphitpc[aniter]->SetNsigmaTPCTOF(kTRUE);
					dtc2etaphitpc[aniter]->SetElectronRejection(ifElectronRejection);
					if (ichg == 0 || ichg == 1 ||ichg == 2 || ichg == 13 || ichg == 14 || ichg == 15 || ichg == 16 || ichg == 20 || ichg == 21 || ichg == 22 || ichg == 23 || ichg == 27 || ichg == 28 || ichg == 29 || ichg == 30) //protons
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
                                        if (ichg == 10 ||ichg == 11 ||ichg == 12)//plus,minus,mixed
                                          {
					    dtc2etaphitpc[aniter]->SetPt(0.2,maxPt);
					  }

                                        if (ichg == 9)//all
                                          {

					    dtc3etaphitpc[aniter] = new AliFemtoMJTrackCut();
					    dtc3etaphitpc[aniter]->SetNsigmaTPCTOF(kTRUE);
					    dtc3etaphitpc[aniter]->SetEta(nEtaMin,nEtaMax);
					    dtc3etaphitpc[aniter]->SetElectronRejection(ifElectronRejection); 
                                            dtc3etaphitpc[aniter]->SetPt(0.2,maxPt);
                                          }

					//*********V0 cuts********************
					//V0 first particle cut
					dtc4etaphitpc[aniter] = new AliFemtoV0TrackCut();
					dtc4etaphitpc[aniter]->SetMass(LambdaMass);
					dtc4etaphitpc[aniter]->SetEta(0.8); //0.8
					if(ichg>=20 && ichg<=26)
					  dtc4etaphitpc[aniter]->SetPt(1.4,maxPt); //0.5,5.0
					if(ichg>=13 && ichg<=19)
					  dtc4etaphitpc[aniter]->SetPt(0.6,1.4);
					if(ichg>=27 && ichg<=33)
					  dtc4etaphitpc[aniter]->SetPt(0.6,maxPt);
					dtc4etaphitpc[aniter]->SetEtaDaughters(0.8);
					dtc4etaphitpc[aniter]->SetPtPosDaughter(0.3,4.0);
					dtc4etaphitpc[aniter]->SetPtNegDaughter(0.16,4.0);
					dtc4etaphitpc[aniter]->SetTPCnclsDaughters(70);
					dtc4etaphitpc[aniter]->SetNdofDaughters(4.0); //4.0
					dtc4etaphitpc[aniter]->SetStatusDaughters(AliESDtrack::kTPCrefit/* | AliESDtrack::kITSrefit*/);
					dtc4etaphitpc[aniter]->SetOnFlyStatus(kFALSE); //kTRUE
					dtc4etaphitpc[aniter]->SetParticleType(0);
					dtc4etaphitpc[aniter]->SetMaxDcaV0Daughters(1.0); //0.5
					dtc4etaphitpc[aniter]->SetMaxDcaV0(0.6); //0.5
					dtc4etaphitpc[aniter]->SetMinDaughtersToPrimVertex(0.06, 0.06); //0.05
					dtc4etaphitpc[aniter]->SetMaxCosPointingAngle(0.99); //0.9993
					dtc4etaphitpc[aniter]->SetMaxV0DecayLength(60.0); //60
					dtc4etaphitpc[aniter]->SetInvariantMassLambda(LambdaMass-0.0038,LambdaMass+0.0038);
					dtc4etaphitpc[aniter]->SetInvariantMassRejectK0s(0.48,0.515);
					dtc4etaphitpc[aniter]->SetRadiusV0Min(0.5);
					dtc4etaphitpc[aniter]->SetNsigmaPosDaughter(5.0);
					dtc4etaphitpc[aniter]->SetNsigmaNegDaughter(5.0);
					dtc4etaphitpc[aniter]->SetRequireTOFPion(false);
					dtc4etaphitpc[aniter]->SetRequireTOFProton(false);
				      
					//V0 second particle cut
					dtc5etaphitpc[aniter] = new AliFemtoV0TrackCut();
					dtc5etaphitpc[aniter]->SetMass(LambdaMass);
					dtc5etaphitpc[aniter]->SetEta(0.8);
					if(ichg>=20 && ichg<=26)
					  dtc5etaphitpc[aniter]->SetPt(1.4,maxPt);
					if(ichg>=13 && ichg<=19)
					  dtc5etaphitpc[aniter]->SetPt(0.6,1.4);
					if(ichg>=27 && ichg<=33)
					  dtc5etaphitpc[aniter]->SetPt(0.6,maxPt);
					dtc5etaphitpc[aniter]->SetEtaDaughters(0.8);
					dtc5etaphitpc[aniter]->SetPtPosDaughter(0.16,4.0);
					dtc5etaphitpc[aniter]->SetPtNegDaughter(0.3,4.0);
					dtc5etaphitpc[aniter]->SetTPCnclsDaughters(70);
					dtc5etaphitpc[aniter]->SetNdofDaughters(4.0); //4.0
					dtc5etaphitpc[aniter]->SetStatusDaughters(AliESDtrack::kTPCrefit/* | AliESDtrack::kITSrefit*/);
					dtc5etaphitpc[aniter]->SetOnFlyStatus(kFALSE); //kTRUE
					dtc5etaphitpc[aniter]->SetParticleType(1);
					dtc5etaphitpc[aniter]->SetMaxDcaV0Daughters(1.0); //0.5
					dtc5etaphitpc[aniter]->SetMaxDcaV0(0.6); //0.5
					dtc5etaphitpc[aniter]->SetMinDaughtersToPrimVertex(0.06, 0.06); //0.05
					dtc5etaphitpc[aniter]->SetMaxCosPointingAngle(0.99); //0.9993
					dtc5etaphitpc[aniter]->SetMaxV0DecayLength(60.0); //60
					dtc5etaphitpc[aniter]->SetInvariantMassLambda(LambdaMass-0.0038,LambdaMass+0.0038);
					dtc5etaphitpc[aniter]->SetInvariantMassRejectK0s(0.48,0.515);
					dtc5etaphitpc[aniter]->SetRadiusV0Min(0.5);
					dtc5etaphitpc[aniter]->SetNsigmaPosDaughter(5.0);
					dtc5etaphitpc[aniter]->SetNsigmaNegDaughter(5.0);
					dtc5etaphitpc[aniter]->SetRequireTOFPion(false);
					dtc5etaphitpc[aniter]->SetRequireTOFProton(false);

					
					

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

					    if(0){
					      // //V0 monitors (memory leak problems?)
					      cutPass1V0[aniter] = new AliFemtoCutMonitorV0(Form("cutPass1%stpcM%i", chrgs[ichg], imult));
					      cutFail1V0[aniter] = new AliFemtoCutMonitorV0(Form("cutFail1%stpcM%i", chrgs[ichg], imult));
					      dtc4etaphitpc[aniter]->AddCutMonitor(cutPass1V0[aniter], cutFail1V0[aniter]);
	  
					      cutPass2V0[aniter] = new AliFemtoCutMonitorV0(Form("cutPass2%stpcM%i", chrgs[ichg], imult));
					      cutFail2V0[aniter] = new AliFemtoCutMonitorV0(Form("cutFail2%stpcM%i", chrgs[ichg], imult));
					      dtc5etaphitpc[aniter]->AddCutMonitor(cutPass2V0[aniter], cutFail2V0[aniter]);

					      anetaphitpc[aniter]->SetV0SharedDaughterCut(performSharedDaughterCut);
					      anetaphitpc[aniter]->SetEnablePairMonitors(enablePairMonitors);
				
					    }

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
					sqpcetaphitpc[aniter] = new AliFemtoPairCutAntiGamma();
					//sqpcetaphitpc[aniter] = new AliFemtoShareQualityPairCut();
					sqpcetaphitpc[aniter]->SetShareQualityMax(shareQuality);	// two track cuts on splitting and merging  //1- wylaczany 0 -wlaczany   
					sqpcetaphitpc[aniter]->SetShareFractionMax(shareFraction);	//  ile moga miec wspolnych klastrow //1 - wylaczany, 0.05 - wlaczany
					sqpcetaphitpc[aniter]->SetRemoveSameLabel(kFALSE);
					// sqpcetaphitpc[aniter]->SetMaximumRadius(0.82);
					// sqpcetaphitpc[aniter]->SetMinimumRadius(0.8);
					// sqpcetaphitpc[aniter]->SetPhiStarDifferenceMinimum(0.02);
					// sqpcetaphitpc[aniter]->SetEtaDifferenceMinimum(0.02);
					
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
					  
					//V0 two-track cuts
			
					
					sqp1cetaphitpc[aniter] = new AliFemtoV0PairCut();
					sqp1cetaphitpc[aniter]->SetDataType(AliFemtoPairCut::kAOD);
					sqp1cetaphitpc[aniter]->SetTPCEntranceSepMinimum(0.00001);
					sqp1cetaphitpc[aniter]->SetTPCExitSepMinimum(-1.);
					sqp1cetaphitpc[aniter]->SetMinAvgSeparation(0,3);
					sqp1cetaphitpc[aniter]->SetMinAvgSeparation(1,0);
					sqp1cetaphitpc[aniter]->SetMinAvgSeparation(2,0);
					sqp1cetaphitpc[aniter]->SetMinAvgSeparation(3,3);

					sqp2cetaphitpc[aniter] = new AliFemtoV0TrackPairCut(); //lambda-proton
					sqp2cetaphitpc[aniter]->SetShareQualityMax(1.0); //between V0 daughter and track
					sqp2cetaphitpc[aniter]->SetShareFractionMax(0.05);
					//sqp2cetaphitpc[aniter]->SetTPCOnly(kTRUE);
					sqp2cetaphitpc[aniter]->SetDataType(AliFemtoPairCut::kAOD);
					sqp2cetaphitpc[aniter]->SetTPCEntranceSepMinimum(0.00001);
					sqp2cetaphitpc[aniter]->SetTPCExitSepMinimum(-1.);
					sqp2cetaphitpc[aniter]->SetKstarCut(0.04,AliFemtoV0TrackPairCut::kLambda,AliFemtoV0TrackPairCut::kProton); //0 - lambda, 2 - proton
					sqp2cetaphitpc[aniter]->SetMinAvgSeparation(0,11); //0 - track-pos, 1 - track-neg
					sqp2cetaphitpc[aniter]->SetMinAvgSeparation(1,0);

					sqp3cetaphitpc[aniter] = new AliFemtoV0TrackPairCut(); //antilambda-antiproton
					sqp3cetaphitpc[aniter]->SetShareQualityMax(1.0); //between V0 daughter and track
					sqp3cetaphitpc[aniter]->SetShareFractionMax(0.05);
					//sqp3cetaphitpc[aniter]->SetTPCOnly(kTRUE);
					sqp3cetaphitpc[aniter]->SetDataType(AliFemtoPairCut::kAOD);
					sqp3cetaphitpc[aniter]->SetTPCEntranceSepMinimum(0.00001);
					sqp3cetaphitpc[aniter]->SetTPCExitSepMinimum(-1.);
					sqp3cetaphitpc[aniter]->SetKstarCut(0.04,AliFemtoV0TrackPairCut::kAntiLambda,AliFemtoV0TrackPairCut::kAntiProton); //1 - antilambda, 3 - antiproton
					sqp3cetaphitpc[aniter]->SetMinAvgSeparation(0,0); //0 - track-pos, 1 - track-neg
					sqp3cetaphitpc[aniter]->SetMinAvgSeparation(1,11);

					sqp4cetaphitpc[aniter] = new AliFemtoV0TrackPairCut(); //lambda-antiproton, antilambda-proton
					sqp4cetaphitpc[aniter]->SetShareQualityMax(1.0); //between V0 daughter and track
					sqp4cetaphitpc[aniter]->SetShareFractionMax(0.05);
					//sqp4cetaphitpc[aniter]->SetTPCOnly(kTRUE);
					sqp4cetaphitpc[aniter]->SetDataType(AliFemtoPairCut::kAOD);
					sqp4cetaphitpc[aniter]->SetTPCEntranceSepMinimum(0.00001);
					sqp4cetaphitpc[aniter]->SetTPCExitSepMinimum(-1.);
					//SetMinAvgSeparation w if'ach ponizej
					
					//***** Setting cuts ***********
					// setting event cut
					anetaphitpc[aniter]->SetEventCut(mecetaphitpc[aniter]);
					//setting single track cuts
					if(ichg==0 || ichg==3  || ichg==6 || ichg==10) //positive like-sign
					{
					  anetaphitpc[aniter]->SetPairCut(sqpcetaphitpc[aniter]);
					  anetaphitpc[aniter]->SetFirstParticleCut(dtc1etaphitpc[aniter]);
					  anetaphitpc[aniter]->SetSecondParticleCut(dtc1etaphitpc[aniter]);
					}
					if(ichg==1 || ichg==4 || ichg==7 || ichg==11)//negative like-sign
					{
					  anetaphitpc[aniter]->SetPairCut(sqpcetaphitpc[aniter]);
					  anetaphitpc[aniter]->SetFirstParticleCut(dtc2etaphitpc[aniter]);
					  anetaphitpc[aniter]->SetSecondParticleCut(dtc2etaphitpc[aniter]);
					}
					if(ichg==2 || ichg==5 || ichg==8 || ichg==12)//unlike-sign
					{
					  anetaphitpc[aniter]->SetPairCut(sqpcetaphitpc[aniter]);
					  anetaphitpc[aniter]->SetFirstParticleCut(dtc1etaphitpc[aniter]);
					  anetaphitpc[aniter]->SetSecondParticleCut(dtc2etaphitpc[aniter]);
					}
					if(ichg==9) //all
                                        {
					  anetaphitpc[aniter]->SetPairCut(sqpcetaphitpc[aniter]);
					  anetaphitpc[aniter]->SetFirstParticleCut(dtc3etaphitpc[aniter]);
					  anetaphitpc[aniter]->SetSecondParticleCut(dtc3etaphitpc[aniter]);
                                        }
					if(ichg == 17 || ichg == 24 || ichg == 31) //V0LL
					  {
					    anetaphitpc[aniter]->SetV0SharedDaughterCut(performSharedDaughterCut);
					    //anetaphitpc[aniter]->SetEventCut(mecetaphitpc[aniter]);
					    anetaphitpc[aniter]->SetFirstParticleCut(dtc4etaphitpc[aniter]);
					    anetaphitpc[aniter]->SetSecondParticleCut(dtc4etaphitpc[aniter]);
					    anetaphitpc[aniter]->SetPairCut(sqp1cetaphitpc[aniter]);
					    //avgsepcorr[aniter]->SetPairType(AliFemtoAvgSepCorrFctn::kV0s);
					  }
					if(ichg == 19 || ichg == 26 || ichg == 33) //V0ALAL
					  {
					    anetaphitpc[aniter]->SetV0SharedDaughterCut(performSharedDaughterCut);
					    //anetaphitpc[aniter]->SetEventCut(mecetaphitpc[aniter]);
					    anetaphitpc[aniter]->SetFirstParticleCut(dtc5etaphitpc[aniter]);
					    anetaphitpc[aniter]->SetSecondParticleCut(dtc5etaphitpc[aniter]);
					    anetaphitpc[aniter]->SetPairCut(sqp1cetaphitpc[aniter]);
					    //avgsepcorr[aniter]->SetPairType(AliFemtoAvgSepCorrFctn::kV0s);
					  }
					if(ichg == 18 || ichg == 25 || ichg == 32) //VOLAL
					  {
					    anetaphitpc[aniter]->SetV0SharedDaughterCut(performSharedDaughterCut);
					    //anetaphitpc[aniter]->SetEventCut(mecetaphitpc[aniter]);
					    anetaphitpc[aniter]->SetFirstParticleCut(dtc4etaphitpc[aniter]);
					    anetaphitpc[aniter]->SetSecondParticleCut(dtc5etaphitpc[aniter]);
					    anetaphitpc[aniter]->SetPairCut(sqp1cetaphitpc[aniter]);
					    //avgsepcorr[aniter]->SetPairType(AliFemtoAvgSepCorrFctn::kV0s);
					  }
					if(ichg == 13 || ichg == 20 || ichg == 27) //V0PL
					  {
					    anetaphitpc[aniter]->SetV0SharedDaughterCut(performSharedDaughterCut);
					    anetaphitpc[aniter]->SetEventCut(mecetaphitpc[aniter]);
					    anetaphitpc[aniter]->SetFirstParticleCut(dtc4etaphitpc[aniter]);
					    anetaphitpc[aniter]->SetSecondParticleCut(dtc1etaphitpc[aniter]);
					    anetaphitpc[aniter]->SetPairCut(sqp2cetaphitpc[aniter]);
					    //avgsepcorr[aniter]->SetPairType(AliFemtoAvgSepCorrFctn::kTrackV0);
					  }
					if(ichg == 15 || ichg == 22 || ichg == 29) //V0APL
					  {
					    sqp4cetaphitpc[aniter]->SetMinAvgSeparation(0,0); //0 - track-pos, 1 - track-neg
					    sqp4cetaphitpc[aniter]->SetMinAvgSeparation(1,11);
					    anetaphitpc[aniter]->SetV0SharedDaughterCut(performSharedDaughterCut);
					    //anetaphitpc[aniter]->SetEventCut(mecetaphitpc[aniter]);
					    anetaphitpc[aniter]->SetFirstParticleCut(dtc4etaphitpc[aniter]);
					    anetaphitpc[aniter]->SetSecondParticleCut(dtc2etaphitpc[aniter]);
					    anetaphitpc[aniter]->SetPairCut(sqp4cetaphitpc[aniter]);
					    //avgsepcorr[aniter]->SetPairType(AliFemtoAvgSepCorrFctn::kTrackV0);
					  }
					if(ichg == 14 || ichg == 21 || ichg == 28) //V0PAL
					  {
					    anetaphitpc[aniter]->SetV0SharedDaughterCut(performSharedDaughterCut);
					    sqp4cetaphitpc[aniter]->SetMinAvgSeparation(0,11); //0 - track-pos, 1 - track-neg
					    sqp4cetaphitpc[aniter]->SetMinAvgSeparation(1,0);
					    //anetaphitpc[aniter]->SetEventCut(mecetaphitpc[aniter]);
					    anetaphitpc[aniter]->SetFirstParticleCut(dtc5etaphitpc[aniter]);
					    anetaphitpc[aniter]->SetSecondParticleCut(dtc1etaphitpc[aniter]);
					    anetaphitpc[aniter]->SetPairCut(sqp4cetaphitpc[aniter]);
					    //avgsepcorr[aniter]->SetPairType(AliFemtoAvgSepCorrFctn::kTrackV0);
					  }
					if(ichg == 16 || ichg == 23 || ichg == 30) //V0APAL
					  {
					    anetaphitpc[aniter]->SetV0SharedDaughterCut(performSharedDaughterCut);
					    //anetaphitpc[aniter]->SetEventCut(mecetaphitpc[aniter]);
					    anetaphitpc[aniter]->SetFirstParticleCut(dtc5etaphitpc[aniter]);
					    anetaphitpc[aniter]->SetSecondParticleCut(dtc2etaphitpc[aniter]);
					    anetaphitpc[aniter]->SetPairCut(sqp3cetaphitpc[aniter]);
					    //avgsepcorr[aniter]->SetPairType(AliFemtoAvgSepCorrFctn::kTrackV0);
					  }



					//**** Correlation functions *******	
					//***without corrections*****
					if(ichg >= 13)
					  cdedpetaphinocorr[aniter] = new AliFemtoCorrFctnDEtaDPhiSimple(Form("cdedpnocorr%stpcM%i", chrgs[ichg], imult),23, 23);
					else
					  cdedpetaphinocorr[aniter] = new AliFemtoCorrFctnDEtaDPhiSimple(Form("cdedpnocorr%stpcM%i", chrgs[ichg], imult),29, 29);
					anetaphitpc[aniter]->AddCorrFctn(cdedpetaphinocorr[aniter]);
					if(ichg==0 || ichg==1 || ichg==31 || ichg==33) //PP, aPaP, LL, ALAL
                                        {
					  cqinvtpc[aniter] = new AliFemtoQinvCorrFctn(Form("cqinv%stpcM%i", chrgs[ichg], imult),nbinssh,0.0,shqmax); //femto qinv, for identical mass particles
					  anetaphitpc[aniter]->AddCorrFctn(cqinvtpc[aniter]);
					}
					if(ichg==2 || ichg==27 || ichg==28 || ichg==29 || ichg==30 || ichg==32) //PaP, PL, APL, PAL, APAL, LAL
                                        {
					  cnonidtpc[aniter] = new AliFemtoCorrFctnNonIdDR(Form("cnonid%stpcM%i", chrgs[ichg], imult), nbinssh, 0.0,shqmax); //for non-identical partcles
					  anetaphitpc[aniter]->AddCorrFctn(cnonidtpc[aniter]);
                                        }


					//***with corrections****
					if(ichg >= 13)
					  cdedpetaphi[aniter] = new AliFemtoCorrFctnDEtaDPhiSimpleWithCorrections(Form("cdedp%stpcM%i", chrgs[ichg], imult),23, 23);
					else
					  cdedpetaphi[aniter] = new AliFemtoCorrFctnDEtaDPhiSimpleWithCorrections(Form("cdedp%stpcM%i", chrgs[ichg], imult),29, 29);

					
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
					else if (ichg==17 || ichg==24 || ichg==31)
					  cdedpetaphi[aniter]->SetParticleTypes(AliFemtoCorrFctnDEtaDPhiSimpleWithCorrections::kLambda, AliFemtoCorrFctnDEtaDPhiSimpleWithCorrections::kLambda);
					else if (ichg==18 || ichg==25 || ichg==32)
					  cdedpetaphi[aniter]->SetParticleTypes(AliFemtoCorrFctnDEtaDPhiSimpleWithCorrections::kLambda, AliFemtoCorrFctnDEtaDPhiSimpleWithCorrections::kLambdaMinus);
					else if (ichg==19 || ichg==26 || ichg==33)
					  cdedpetaphi[aniter]->SetParticleTypes(AliFemtoCorrFctnDEtaDPhiSimpleWithCorrections::kLambdaMinus, AliFemtoCorrFctnDEtaDPhiSimpleWithCorrections::kLambdaMinus);


					
					anetaphitpc[aniter]->AddCorrFctn(cdedpetaphi[aniter]);
						       							
					Manager->AddAnalysis(anetaphitpc[aniter]);
				}
			}
		}
	}
					    
	return Manager;
}								 
