 /* Configfemtoanalysis.C - configuration macro for the femtoscopic	 *
 /********************************************************************* *																							 *
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
 #include "AliFemtoEventReaderNanoAODChain.h"
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
 #include "AliFemtoPairCutRadialDistance.h"
 #include "AliFemtoCorrFctnDEtaDPhiStar.h"

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
	double XiMass = 1.32171;
  double OmegaMass = 1.67245;


	const int numOfMultBins = 5;
	const int numOfChTypes = 42; //34 + 4 xi + 4 Omega
	const int numOfkTbins = 5;

	bool performSharedDaughterCut = true;
	bool enablePairMonitors = true;

	char *par = new char[strlen(params)+1];
        strcpy(par,params);

	char *parameter[21];
	if(strlen(params)!=0)
	  {
	    parameter[0] = strtok(par, ","); // Splits spaces between words in params
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
	int multbins[numOfMultBins+1] = {2, 20, 50,150,2,1500};

	int runch[numOfChTypes] = {/*protons*/1, 1, 1, /* kaons */ 1, 1, 1, /* pions */ 1, 1, 1, /* no PID */ 0, 0, 0, 0,/*V0 pT*/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, /*p-lam */ 1, 1, 1, 1, /* lambdas */ 1, 1, 1,/* p-Xi */ 1, 1, 1, 1,/* p-Omega*/ 1, 1, 1, 1};
	const char *chrgs[numOfChTypes] = { "PP", "aPaP", "PaP", "KpKp", "KmKm", "KpKm", "PIpPIp", "PImPIm", "PIpPIm", "all", "plus", "minus", "mixed", "V0PLlowPt","V0PALlowPt","V0APLlowPt","V0APALlowPt","V0LLlowPt","V0LALlowPt","V0ALALlowPt", "V0PLhighPt","V0PALhighPt","V0APLhighPt","V0APALhighPt","V0LLhighPt","V0LALhighPt","V0ALALhighPt", "V0PL","V0PAL","V0APL","V0APAL","V0LL","V0LAL","V0ALAL","PXim","aPXim" ,"PXip","aPXip", ,"POmegam","aPOmegam" ,"POmegap","aPOmegap"};


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


	AliFemtoEventReaderNanoAODChain *Reader = new AliFemtoEventReaderNanoAODChain();
	Reader->SetFilterMask(filterbit);
	Reader->SetCovMatPresent(false);
	Reader->SetDCAglobalTrack(1); //false for FB7, true for the rest //we do not use DCA at all
	Reader->SetUseMultiplicity("V0M");

  Reader->SetReadV0(kTRUE);
  Reader->SetReadCascade(kTRUE);

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
	AliFemtoXiTrackCut           *tXiCut[numOfMultBins*numOfChTypes];
	AliFemtoXiTrackCut           *tAXiCut[numOfMultBins*numOfChTypes];
  AliFemtoXiTrackCut           *tOmegaCut[numOfMultBins*numOfChTypes];
  AliFemtoXiTrackCut           *tAOmegaCut[numOfMultBins*numOfChTypes];
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
	AliFemtoCutMonitorXi             *cutPass1Xi[numOfMultBins*numOfChTypes];
	AliFemtoCutMonitorXi             *cutFail1Xi[numOfMultBins*numOfChTypes];
	AliFemtoCutMonitorXi             *cutPass2Xi[numOfMultBins*numOfChTypes];
	AliFemtoCutMonitorXi             *cutFail2Xi[numOfMultBins*numOfChTypes];
  AliFemtoCutMonitorXi             *cutPass1Omega[numOfMultBins*numOfChTypes];
  AliFemtoCutMonitorXi             *cutFail1Omega[numOfMultBins*numOfChTypes];
  AliFemtoCutMonitorXi             *cutPass2Omega[numOfMultBins*numOfChTypes];
  AliFemtoCutMonitorXi             *cutFail2Omega[numOfMultBins*numOfChTypes];
	AliFemtoPairCutRadialDistance  *sqpcetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoV0PairCut               *sqp1cetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoV0TrackPairCut          *sqp2cetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoV0TrackPairCut          *sqp3cetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoV0TrackPairCut          *sqp4cetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoXiTrackPairCut           *tXiTrackPairCut[numOfMultBins*numOfChTypes];
  AliFemtoXiTrackPairCut           *tOmegaTrackPairCut[numOfMultBins*numOfChTypes];

//	AliFemtoChi2CorrFctn					*cchiqinvetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoPairCutPt               *ktpcuts[numOfMultBins*numOfChTypes*numOfkTbins];
	AliFemtoQinvCorrFctn		*cqinvkttpc[numOfMultBins*numOfChTypes*numOfkTbins];
	AliFemtoQinvCorrFctn		*cqinvtpc[numOfMultBins*numOfChTypes];
	AliFemtoCorrFctnDEtaDPhiSimpleWithCorrections *cdedpetaphi[numOfMultBins*numOfChTypes];
	AliFemtoCorrFctnDEtaDPhiSimple	*cdedpetaphinocorr[numOfMultBins*numOfChTypes];
	//	AliFemtoCorrFctnDYDPhiSimpleWithCorrections	*cdydpyphinocorr[numOfMultBins*numOfChTypes];
	AliFemtoCorrFctnDEtaDPhiCorrections *cdedpetaphiPt[numOfMultBins*numOfChTypes*numOfkTbins];

	AliFemtoCorrFctnDPhiStarDEta  *cdphistardeta08[5000];
	AliFemtoCorrFctnDPhiStarDEta  *cdphistardeta12[500];
	AliFemtoCorrFctnDPhiStarDEta  *cdphistardeta16[500];
	AliFemtoCorrFctnDPhiStarDEta  *cdphistardeta20[500];

        AliFemtoCorrFctnNonIdDR         *cnonidtpc[numOfMultBins*numOfChTypes];
        AliFemtoCorrFctnDEtaDPhiStar  *cdetadphistar[numOfMultBins*numOfChTypes];


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
					anetaphitpc[aniter]->SetNumEventsToMix(5);
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
					if (ichg == 0 || ichg == 1 ||ichg == 2 || ichg == 13 || ichg == 14 || ichg == 15 || ichg == 16 || ichg == 20 || ichg == 21 || ichg == 22 || ichg == 23 || ichg == 27 || ichg == 28 || ichg == 29 || ichg == 30 || ichg==34 || ichg==36 || ichg == 38 || ichg == 40) //protons
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
					if (ichg == 0 || ichg == 1 ||ichg == 2 || ichg == 13 || ichg == 14 || ichg == 15 || ichg == 16 || ichg == 20 || ichg == 21 || ichg == 22 || ichg == 23 || ichg == 27 || ichg == 28 || ichg == 29 || ichg == 30 || ichg==35 || ichg==37|| ichg == 39 || ichg == 41) //antiprotons
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
					dtc4etaphitpc[aniter] = new AliFemtoV0TrackCut(); //Lambda
					dtc4etaphitpc[aniter]->SetMass(LambdaMass);
					dtc4etaphitpc[aniter]->SetEta(0.8); //0.8
          dtc4etaphitpc[aniter]->SetNanoAODAnalysis(kTRUE);
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
					dtc4etaphitpc[aniter]->SetOnFlyStatus(kFALSE); //kTRUE
					dtc4etaphitpc[aniter]->SetParticleType(0);
					dtc4etaphitpc[aniter]->SetMaxDcaV0Daughters(1.0); //0.5
					dtc4etaphitpc[aniter]->SetMaxDcaV0(0.6); //0.5
					dtc4etaphitpc[aniter]->SetMinDaughtersToPrimVertex(0.06, 0.06); //0.05
					dtc4etaphitpc[aniter]->SetMaxCosPointingAngle(0.99); //0.9993
					dtc4etaphitpc[aniter]->SetMaxV0DecayLength(60.0); //60
					//dtc4etaphitpc[aniter]->SetInvariantMassLambda(LambdaMass-0.0038,LambdaMass+0.0038);
          dtc4etaphitpc[aniter]->SetSidebandAnalysis(kTRUE);
          dtc4etaphitpc[aniter]->SetInvariantMassLambdaSideband(0.905, 1.105, 1.275, 1.475 );
					dtc4etaphitpc[aniter]->SetInvariantMassRejectK0s(0.48,0.515);
					dtc4etaphitpc[aniter]->SetRadiusV0Min(0.5);
					dtc4etaphitpc[aniter]->SetNsigmaPosDaughter(5.0);
					dtc4etaphitpc[aniter]->SetNsigmaNegDaughter(5.0);
					dtc4etaphitpc[aniter]->SetRequireTOFPion(false);
					dtc4etaphitpc[aniter]->SetRequireTOFProton(false);

					//V0 second particle cut
					dtc5etaphitpc[aniter] = new AliFemtoV0TrackCut(); //Antilambda
					dtc5etaphitpc[aniter]->SetMass(LambdaMass);
					dtc5etaphitpc[aniter]->SetEta(0.8);
          dtc5etaphitpc[aniter]->SetNanoAODAnalysis(kTRUE);
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
					dtc5etaphitpc[aniter]->SetOnFlyStatus(kFALSE); //kTRUE
					dtc5etaphitpc[aniter]->SetParticleType(1);
					dtc5etaphitpc[aniter]->SetMaxDcaV0Daughters(1.0); //0.5
					dtc5etaphitpc[aniter]->SetMaxDcaV0(0.6); //0.5
					dtc5etaphitpc[aniter]->SetMinDaughtersToPrimVertex(0.06, 0.06); //0.05
					dtc5etaphitpc[aniter]->SetMaxCosPointingAngle(0.99); //0.9993
					dtc5etaphitpc[aniter]->SetMaxV0DecayLength(60.0); //60
					//dtc5etaphitpc[aniter]->SetInvariantMassLambda(LambdaMass-0.0038,LambdaMass+0.0038);
          dtc5etaphitpc[aniter]->SetSidebandAnalysis(kTRUE);
          dtc5etaphitpc[aniter]->SetInvariantMassAntiLambdaSideband(0.905, 1.105, 1.275, 1.475 );
					dtc5etaphitpc[aniter]->SetInvariantMassRejectK0s(0.48,0.515);
					dtc5etaphitpc[aniter]->SetRadiusV0Min(0.5);
					dtc5etaphitpc[aniter]->SetNsigmaPosDaughter(5.0);
					dtc5etaphitpc[aniter]->SetNsigmaNegDaughter(5.0);
					dtc5etaphitpc[aniter]->SetRequireTOFPion(false);
					dtc5etaphitpc[aniter]->SetRequireTOFProton(false);

					//Cascade cuts
					//from J. Buxton
					//xi cut
					//NOTE: the SetMass call actually is important
					//      This should be set to the mass of the particle of interest, here the Xi
					//      Be sure to not accidentally set it again in the Lambda cuts (for instance, when copy/pasting the lambda cuts from above!)

					//Xi -> Lam Pi-
          if(ichg == 34 || ichg == 35 || ichg == 36 || ichg == 37){

  					tXiCut[aniter] = new AliFemtoXiTrackCut();
  					// %%%%%%%%%%%%%%%%%%%%%%%% Version 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  					tXiCut[aniter]->SetParticleTypeXi(AliFemtoXiTrackCut::kXiMinus);  //kXiMinus = 0 //sprawdza Nsigma<3 czy bachelor jest pionem
  					tXiCut[aniter]->SetChargeXi(-1);
  					tXiCut[aniter]->SetPtXi(0.5,100);
  					tXiCut[aniter]->SetEtaXi(0.8);
  					tXiCut[aniter]->SetMass(XiMass);
  	        tXiCut[aniter]->SetNanoAODAnalysis(kTRUE);
  				//	tXiCut[aniter]->SetInvariantMassXi(XiMass-0.005,XiMass+0.005); //++ bylo 006
          tXiCut[aniter]->SetSidebandAnalysis(kTRUE);
          tXiCut[aniter]->SetInvariantMassXiSideband(XiMass - 0.35, XiMass-0.15, XiMass+0.15, XiMass+0.35 );

  					tXiCut[aniter]->SetMinCosPointingAngleXi(0.97); //++ bylo 0.99
  					tXiCut[aniter]->SetMaxDecayLengthXi(100.);
  					tXiCut[aniter]->SetMaxDcaXi(100);
  					tXiCut[aniter]->SetInvariantMassRejectOmega(1.667,1.677);//++ NEW: omega rejection od 1.667 do 1.677 !


  					//XiDaughters
  					tXiCut[aniter]->SetMaxDcaXiDaughters(1.6);//++ bylo 0.3
  					tXiCut[aniter]->SetRadiusXiMin(0.8); //++ NEW!
  					tXiCut[aniter]->SetRadiusXiMax(200); //++ NEW!

  					tXiCut[aniter]->SetIgnoreOnFlyStatus(kTRUE);
  					//Bachelor cuts (here = PiM)
  					tXiCut[aniter]->SetMinDcaXiBac(0.05); //++ bylo 0.03
  					tXiCut[aniter]->SetEtaBac(0.8);
  					tXiCut[aniter]->SetTPCnclsBac(70); //++
  					tXiCut[aniter]->SetPtBac(0.3,100.);//++
  					//++ brakuje cut-u z bachelor do V0 (<1.3)!!!!++++++++++


  					//Lambda cuts (regular V0)
  					tXiCut[aniter]->SetParticleType(AliFemtoV0TrackCut::kLambda); //0=lambda
  					tXiCut[aniter]->SetMinDcaV0(0.07); //++ bylo 0.1
  					tXiCut[aniter]->SetInvariantMassLambda(LambdaMass-0.005,LambdaMass+0.005);
  					tXiCut[aniter]->SetMinCosPointingAngle(0.97); //++ bylo 0.998
  					tXiCut[aniter]->SetEta(0.8);
  					tXiCut[aniter]->SetPt(0.0,100);
  					//tXiCut[aniter]->SetOnFlyStatus(kFALSE);
  					tXiCut[aniter]->SetMaxV0DecayLength(100.);
  					tXiCut[aniter]->SetRadiusV0Min(1.4); //++ NEW!
  					tXiCut[aniter]->SetRadiusV0Max(200); //++ NEW!

  					//Lambda daughter cuts
  					tXiCut[aniter]->SetMinDaughtersToPrimVertex(0.04,0.04); //++   pierwsza pos, druga neg, bylio (0.1,0.1); albo  pion 0.04, proton 0.03
  					tXiCut[aniter]->SetMaxDcaV0Daughters(1.6); //++ bylo 0.8
  					tXiCut[aniter]->SetEtaDaughters(0.8);    //++
  					tXiCut[aniter]->SetPtPosDaughter(0.3,99); //++
  					tXiCut[aniter]->SetPtNegDaughter(0.3,99); //++
  					tXiCut[aniter]->SetTPCnclsDaughters(70); //++
  				//	tXiCut[aniter]->SetNanoAODAnalysis(kTRUE);  //yes or no?
  					tXiCut[aniter]->SetNsigmaPosDaughter(5.0); //++
  					tXiCut[aniter]->SetNsigmaNegDaughter(5.0); //++
  					tXiCut[aniter]->SetRequireTOFProton(kFALSE); //++
  					// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  					tXiCut[aniter]->SetMinvPurityAidHistoXi("XiPurityAid","XiMinvBeforeFinalCut",100,XiMass-0.035,XiMass+0.035);
  					tXiCut[aniter]->SetMinvPurityAidHistoV0("LambdaPurityAid","LambdaMinvBeforeFinalCut",100,LambdaMass-0.035,LambdaMass+0.035);


  					//antiXi cut
  					//NOTE: the SetMass call actually is important
  					//      This should be set to the mass of the particle of interest, here the Xi
  					//      Be sure to not accidentally set it again in the Lambda cuts (for instance, when copy/pasting the lambda cuts from above!)

  					//AXi -> ALam Pi+

  					tAXiCut[aniter] = new AliFemtoXiTrackCut();

  					// %%%%%%%%%%%%%%%%%%%%%%%% Version 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  					tAXiCut[aniter]->SetParticleTypeXi(AliFemtoXiTrackCut::kXiPlus); //kXiPlus = 1 //sprawdza Nsigma<3 czy bachelor jest pionem
  					tAXiCut[aniter]->SetChargeXi(1);
  					tAXiCut[aniter]->SetPtXi(0.5,100);
  					tAXiCut[aniter]->SetEtaXi(0.8);
  					tAXiCut[aniter]->SetMass(XiMass);
  					//tAXiCut[aniter]->SetInvariantMassXi(XiMass-0.005,XiMass+0.005);
            tXiCut[aniter]->SetSidebandAnalysis(kTRUE);
            tXiCut[aniter]->SetInvariantMassXiSideband(XiMass - 0.35, XiMass-0.15, XiMass+0.15, XiMass+0.35 );
            
  					tAXiCut[aniter]->SetMinCosPointingAngleXi(0.97);
  					tAXiCut[aniter]->SetMaxDecayLengthXi(100.0);
  					tAXiCut[aniter]->SetMaxDcaXi(100);
  					tAXiCut[aniter]->SetInvariantMassRejectOmega(1.667,1.677);//++ NEW: omega rejection od 1.667 do 1.677 !


  					tAXiCut[aniter]->SetIgnoreOnFlyStatus(kTRUE);
  					//XiDaughters
  					tAXiCut[aniter]->SetMaxDcaXiDaughters(1.6);
  					tAXiCut[aniter]->SetRadiusXiMin(0.8); //++ NEW!
  					tAXiCut[aniter]->SetRadiusXiMax(200); //++ NEW!

  					//Bachelor cuts (here = PiP)
  					tAXiCut[aniter]->SetMinDcaXiBac(0.05);
  					tAXiCut[aniter]->SetEtaBac(0.8);
  					tAXiCut[aniter]->SetTPCnclsBac(70);
  					tAXiCut[aniter]->SetPtBac(0.3,100);
  					tAXiCut[aniter]->SetNanoAODAnalysis(kTRUE); //yes or no?

  					//AntiLambda cuts (regular V0)
  					tAXiCut[aniter]->SetParticleType(AliFemtoV0TrackCut::kAntiLambda); //1=anti-lambda
  					tAXiCut[aniter]->SetMinDcaV0(0.07);
  					tAXiCut[aniter]->SetInvariantMassLambda(LambdaMass-0.005,LambdaMass+0.005);
  					tAXiCut[aniter]->SetMinCosPointingAngle(0.97);
  					tAXiCut[aniter]->SetEta(0.8);
  					tAXiCut[aniter]->SetPt(0.,100);
  					//tAXiCut[aniter]->SetOnFlyStatus(kFALSE);  //CHECK kTRUE STATUS AS WELL?
  					tAXiCut[aniter]->SetMaxV0DecayLength(100.);
  					tAXiCut[aniter]->SetRadiusV0Min(1.4); //++ NEW!
  					tAXiCut[aniter]->SetRadiusV0Max(200); //++ NEW!
  					//Lambda daughter cuts
  					tAXiCut[aniter]->SetMinDaughtersToPrimVertex(0.04,0.04);
  					tAXiCut[aniter]->SetMaxDcaV0Daughters(1.6);
  					tAXiCut[aniter]->SetEtaDaughters(0.8);
  					tAXiCut[aniter]->SetPtPosDaughter(0.3,99); //0.16 for pions
  					tAXiCut[aniter]->SetPtNegDaughter(0.3,99); //0.5 for anti-protons
  					tAXiCut[aniter]->SetTPCnclsDaughters(70);
  					//tAXiCut[aniter]->SetStatusDaughters(AliESDtrack::kTPCrefit);  //yes or no?
  					tAXiCut[aniter]->SetNsigmaPosDaughter(5.0); //++
  					tAXiCut[aniter]->SetNsigmaNegDaughter(5.0); //++
  					tAXiCut[aniter]->SetRequireTOFProton(kFALSE); //++

  					// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  					tAXiCut[aniter]->SetMinvPurityAidHistoXi("AXiPurityAid","AXiMinvBeforeFinalCut",100,XiMass-0.035,XiMass+0.035);
  					tAXiCut[aniter]->SetMinvPurityAidHistoV0("AntiLambdaPurityAid","AntiLambdaMinvBeforeFinalCut",100,LambdaMass-0.035,LambdaMass+0.035);
        }

        if(ichg == 38 || ichg == 39 || ichg == 40 || ichg == 41){
        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OMEGA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            tOmegaCut[aniter] = new AliFemtoXiTrackCut();
            tOmegaCut[aniter]->SetParticleTypeXi(AliFemtoXiTrackCut::kOmegaMinus);  //kXiMinus = 0 //sprawdza Nsigma<3 czy bachelor jest pionem
            tOmegaCut[aniter]->SetChargeXi(-1);
            tOmegaCut[aniter]->SetPtXi(0.5,100);
            tOmegaCut[aniter]->SetEtaXi(0.8);
            tOmegaCut[aniter]->SetMass(OmegaMass);
            tOmegaCut[aniter]->SetNanoAODAnalysis(kTRUE);
            //tOmegaCut[aniter]->SetInvariantMassOmega(OmegaMass-0.0575,OmegaMass+0.0575); //++ bylo 006
            tOmegaCut[aniter]->SetSidebandAnalysis(kTRUE);
            tOmegaCut[aniter]->SetInvariantMassOmegaSideband(1.642, 1.658, 1.684, 1.7 );

            tOmegaCut[aniter]->SetMinCosPointingAngleXi(0.97); //++ bylo 0.99
            tOmegaCut[aniter]->SetMaxDecayLengthXi(100.);
            tOmegaCut[aniter]->SetMaxDcaXi(100);
            tOmegaCut[aniter]->SetInvariantMassRejectXi(XiMass - 0.005,XiMass + 0.005);//++ NEW: omega rejection od 1.667 do 1.677 !

            //XiDaughters
            tOmegaCut[aniter]->SetMaxDcaXiDaughters(1.6);//++ bylo 0.3
            tOmegaCut[aniter]->SetRadiusXiMin(0.8); //++ NEW!
            tOmegaCut[aniter]->SetRadiusXiMax(200); //++ NEW!

            tOmegaCut[aniter]->SetIgnoreOnFlyStatus(kTRUE);
            //Bachelor cuts (here = PiM)
            tOmegaCut[aniter]->SetMinDcaXiBac(0.03); //++ bylo 0.03
            tOmegaCut[aniter]->SetEtaBac(0.8);
            tOmegaCut[aniter]->SetTPCnclsBac(70); //++
            tOmegaCut[aniter]->SetPtBac(0.3,100.);//++
            //++ brakuje cut-u z bachelor do V0 (<1.3)!!!!++++++++++

            //Lambda cuts (regular V0)
            tOmegaCut[aniter]->SetParticleType(AliFemtoV0TrackCut::kLambda); //0=lambda
            tOmegaCut[aniter]->SetMinDcaV0(0.07); //++ bylo 0.1
            tOmegaCut[aniter]->SetInvariantMassLambda(LambdaMass-0.005,LambdaMass+0.005);
            tOmegaCut[aniter]->SetMinCosPointingAngle(0.97); //++ bylo 0.998
            tOmegaCut[aniter]->SetEta(0.8);
            tOmegaCut[aniter]->SetPt(0.0,100);
            //tOmegaCut[aniter]->SetOnFlyStatus(kFALSE);
            tOmegaCut[aniter]->SetMaxV0DecayLength(100.);
            tOmegaCut[aniter]->SetRadiusV0Min(1.4); //++ NEW!
            tOmegaCut[aniter]->SetRadiusV0Max(200); //++ NEW!

            //Lambda daughter cuts
            tOmegaCut[aniter]->SetMinDaughtersToPrimVertex(0.04,0.04); //++   pierwsza pos, druga neg, bylio (0.1,0.1); albo  pion 0.04, proton 0.03
            tOmegaCut[aniter]->SetMaxDcaV0Daughters(1.6); //++ bylo 0.8
            tOmegaCut[aniter]->SetEtaDaughters(0.8);    //++
            tOmegaCut[aniter]->SetPtPosDaughter(0.3,99); //++
            tOmegaCut[aniter]->SetPtNegDaughter(0.3,99); //++
            tOmegaCut[aniter]->SetTPCnclsDaughters(70); //++
            //	tOmegaCut[aniter]->SetNanoAODAnalysis(kTRUE);  //yes or no?
            tOmegaCut[aniter]->SetNsigmaPosDaughter(5.0); //++
            tOmegaCut[aniter]->SetNsigmaNegDaughter(5.0); //++
            tOmegaCut[aniter]->SetRequireTOFProton(kFALSE); //++
            // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            tOmegaCut[aniter]->SetMinvPurityAidHistoXi("OmegaPurityAid","OmegaMinvBeforeFinalCut",100,OmegaMass-0.035,OmegaMass+0.035);
            tOmegaCut[aniter]->SetMinvPurityAidHistoV0("LambdaPurityAid","LambdaMinvBeforeFinalCut",100,LambdaMass-0.035,LambdaMass+0.035);


            //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Anty Omega %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            tAOmegaCut[aniter] = new AliFemtoXiTrackCut();

            tAOmegaCut[aniter]->SetParticleTypeXi(AliFemtoXiTrackCut::kOmegaPlus); //kXiPlus = 1 //sprawdza Nsigma<3 czy bachelor jest pionem
            tAOmegaCut[aniter]->SetChargeXi(1);
            tAOmegaCut[aniter]->SetPtXi(0.5,100);
            tAOmegaCut[aniter]->SetEtaXi(0.8);
            tAOmegaCut[aniter]->SetMass(OmegaMass);
          //  tAOmegaCut[aniter]->SetInvariantMassOmega(OmegaMass-0.0575,OmegaMass+0.0575);
            tAOmegaCut[aniter]->SetSidebandAnalysis(kTRUE);
            tAOmegaCut[aniter]->SetInvariantMassOmegaSideband(1.642, 1.658, 1.684, 1.7 );

            tAOmegaCut[aniter]->SetMinCosPointingAngleXi(0.97);
            tAOmegaCut[aniter]->SetMaxDecayLengthXi(100.0);
            tAOmegaCut[aniter]->SetMaxDcaXi(100);
            tAOmegaCut[aniter]->SetInvariantMassRejectXi(XiMass - 0.005, XiMass + 0.005);//++ NEW: omega rejection od 1.667 do 1.677 !

            tAOmegaCut[aniter]->SetIgnoreOnFlyStatus(kTRUE);
            //XiDaughters
            tAOmegaCut[aniter]->SetMaxDcaXiDaughters(1.6);
            tAOmegaCut[aniter]->SetRadiusXiMin(0.8); //++ NEW!
            tAOmegaCut[aniter]->SetRadiusXiMax(200); //++ NEW!

            //Bachelor cuts (here = KP)
            tAOmegaCut[aniter]->SetMinDcaXiBac(0.05);
            tAOmegaCut[aniter]->SetEtaBac(0.8);
            tAOmegaCut[aniter]->SetTPCnclsBac(70);
            tAOmegaCut[aniter]->SetPtBac(0.3,100);
            tAOmegaCut[aniter]->SetNanoAODAnalysis(kTRUE); //yes or no?

            //AntiLambda cuts (regular V0)
            tAOmegaCut[aniter]->SetParticleType(AliFemtoV0TrackCut::kAntiLambda); //1=anti-lambda
            tAOmegaCut[aniter]->SetMinDcaV0(0.07);
            tAOmegaCut[aniter]->SetInvariantMassLambda(LambdaMass-0.005,LambdaMass+0.005);
            tAOmegaCut[aniter]->SetMinCosPointingAngle(0.97);
            tAOmegaCut[aniter]->SetEta(0.8);
            tAOmegaCut[aniter]->SetPt(0.,100);
            //tAOmegaCut[aniter]->SetOnFlyStatus(kFALSE);  //CHECK kTRUE STATUS AS WELL?
            tAOmegaCut[aniter]->SetMaxV0DecayLength(100.);
            tAOmegaCut[aniter]->SetRadiusV0Min(1.4); //++ NEW!
            tAOmegaCut[aniter]->SetRadiusV0Max(200); //++ NEW!
            //Lambda daughter cuts
            tAOmegaCut[aniter]->SetMinDaughtersToPrimVertex(0.04,0.04);
            tAOmegaCut[aniter]->SetMaxDcaV0Daughters(1.6);
            tAOmegaCut[aniter]->SetEtaDaughters(0.8);
            tAOmegaCut[aniter]->SetPtPosDaughter(0.3,99); //0.16 for pions
            tAOmegaCut[aniter]->SetPtNegDaughter(0.3,99); //0.5 for anti-protons
            tAOmegaCut[aniter]->SetTPCnclsDaughters(70);
            //tAOmegaCut[aniter]->SetStatusDaughters(AliESDtrack::kTPCrefit);  //yes or no?
            tAOmegaCut[aniter]->SetNsigmaPosDaughter(5.0); //++
            tAOmegaCut[aniter]->SetNsigmaNegDaughter(5.0); //++
            tAOmegaCut[aniter]->SetRequireTOFProton(kFALSE); //++

            // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            tAOmegaCut[aniter]->SetMinvPurityAidHistoXi("AOmegaPurityAid","AOmegaMinvBeforeFinalCut",100,OmegaMass-0.035,OmegaMass+0.035);
            tAOmegaCut[aniter]->SetMinvPurityAidHistoV0("AntiLambdaPurityAid","AntiLambdaMinvBeforeFinalCut",100,LambdaMass-0.035,LambdaMass+0.035);
          }


					//**************** track Monitors ***************

					if(ifMonitors)//ichg>8)
					  {

					    if(ichg>=17 && ichg<=33){
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
					    else if(ichg==34||ichg==35||ichg==36||ichg==37) //PXim, aPXim, PXip, aPXip
					      {
						     cutPass1Xi[aniter] = new AliFemtoCutMonitorXi(Form("cutPass%stpcM%i", chrgs[ichg], imult));
						     cutFail1Xi[aniter] = new AliFemtoCutMonitorXi(Form("cutFail%stpcM%i", chrgs[ichg], imult));
						     if(ichg==34||ichg==35) tXiCut[aniter]->AddCutMonitor(cutPass1Xi[aniter],cutFail1Xi[aniter]);
						     if(ichg==36||ichg==37) tAXiCut[aniter]->AddCutMonitor(cutPass1Xi[aniter],cutFail1Xi[aniter]);
					      }
                else if(ichg==38||ichg==39||ichg==40||ichg==41) //POm, aPOm, POp, aPOp
                  {
              cutPass1Omega[aniter] = new AliFemtoCutMonitorXi(Form("cutPass%stpcM%i", chrgs[ichg], imult));
              cutFail1Omega[aniter] = new AliFemtoCutMonitorXi(Form("cutFail%stpcM%i", chrgs[ichg], imult));
              if(ichg==38||ichg==40) tOmegaCut[aniter]->AddCutMonitor(cutPass1Omega[aniter],cutFail1Omega[aniter]);
              if(ichg==39||ichg==41) tAOmegaCut[aniter]->AddCutMonitor(cutPass1Omega[aniter],cutFail1Omega[aniter]);
                  }



					  }

					//******** Two - track cuts ************
					sqpcetaphitpc[aniter] = new AliFemtoPairCutRadialDistance();
					sqpcetaphitpc[aniter]->SetRemoveSameLabel(kFALSE);
					sqpcetaphitpc[aniter]->SetMinimumRadius(1.6);
					sqpcetaphitpc[aniter]->SetEtaDifferenceMinimum(0.02);
					sqpcetaphitpc[aniter]->SetPhiStarDifferenceMinimum(0.045);
          sqpcetaphitpc[aniter]->SetNanoAODAnalysis(kTRUE);

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

  				//V0 two-track cuts

					sqp1cetaphitpc[aniter] = new AliFemtoV0PairCut();
          sqp1cetaphitpc[aniter]->SetNanoAODAnalysis(kTRUE);
					sqp1cetaphitpc[aniter]->SetDataType(AliFemtoPairCut::kAOD);
					sqp1cetaphitpc[aniter]->SetTPCEntranceSepMinimum(0.00001);
					sqp1cetaphitpc[aniter]->SetTPCExitSepMinimum(-1.);

					sqp2cetaphitpc[aniter] = new AliFemtoV0TrackPairCut(); //lambda-proton
          sqp2cetaphitpc[aniter]->SetNanoAODAnalysis(kTRUE);
					sqp2cetaphitpc[aniter]->SetShareQualityMax(1.0); //between V0 daughter and track
					sqp2cetaphitpc[aniter]->SetShareFractionMax(0.05);
					sqp2cetaphitpc[aniter]->SetTPCOnly(kTRUE);
					sqp2cetaphitpc[aniter]->SetDataType(AliFemtoPairCut::kAOD);
					sqp2cetaphitpc[aniter]->SetTPCEntranceSepMinimum(0.00001);
					sqp2cetaphitpc[aniter]->SetTPCExitSepMinimum(-1.);
					sqp2cetaphitpc[aniter]->SetKstarCut(0.04,AliFemtoV0TrackPairCut::kLambda,AliFemtoV0TrackPairCut::kProton); //0 - lambda, 2 - proton


					sqp3cetaphitpc[aniter] = new AliFemtoV0TrackPairCut(); //antilambda-antiproton
          sqp3cetaphitpc[aniter]->SetNanoAODAnalysis(kTRUE);
					sqp3cetaphitpc[aniter]->SetShareQualityMax(1.0); //between V0 daughter and track
					sqp3cetaphitpc[aniter]->SetShareFractionMax(0.05);
					sqp3cetaphitpc[aniter]->SetTPCOnly(kTRUE);
					sqp3cetaphitpc[aniter]->SetDataType(AliFemtoPairCut::kAOD);
					sqp3cetaphitpc[aniter]->SetTPCEntranceSepMinimum(0.00001);
					sqp3cetaphitpc[aniter]->SetTPCExitSepMinimum(-1.);
					sqp3cetaphitpc[aniter]->SetKstarCut(0.04,AliFemtoV0TrackPairCut::kAntiLambda,AliFemtoV0TrackPairCut::kAntiProton); //1 - antilambda, 3 - antiproton

					sqp4cetaphitpc[aniter] = new AliFemtoV0TrackPairCut(); //lambda-antiproton, antilambda-proton
          sqp4cetaphitpc[aniter]->SetNanoAODAnalysis(kTRUE);
					sqp4cetaphitpc[aniter]->SetShareQualityMax(1.0); //between V0 daughter and track
					sqp4cetaphitpc[aniter]->SetShareFractionMax(0.05);
					sqp4cetaphitpc[aniter]->SetTPCOnly(kTRUE);
					sqp4cetaphitpc[aniter]->SetDataType(AliFemtoPairCut::kAOD);
					sqp4cetaphitpc[aniter]->SetTPCEntranceSepMinimum(0.00001);
					sqp4cetaphitpc[aniter]->SetTPCExitSepMinimum(-1.);
					//SetMinAvgSeparation w if'ach ponizej


					tXiTrackPairCut[aniter] = new AliFemtoXiTrackPairCut(); //xi-proton, all combinations
	        tOmegaTrackPairCut[aniter] = new AliFemtoXiTrackPairCut();


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
					    anetaphitpc[aniter]->SetEventCut(mecetaphitpc[aniter]);
					    anetaphitpc[aniter]->SetFirstParticleCut(dtc4etaphitpc[aniter]);
					    anetaphitpc[aniter]->SetSecondParticleCut(dtc4etaphitpc[aniter]);
					    anetaphitpc[aniter]->SetPairCut(sqp1cetaphitpc[aniter]);

					  }
					if(ichg == 19 || ichg == 26 || ichg == 33) //V0ALAL
					  {
					    anetaphitpc[aniter]->SetV0SharedDaughterCut(performSharedDaughterCut);
					    anetaphitpc[aniter]->SetEventCut(mecetaphitpc[aniter]);
					    anetaphitpc[aniter]->SetFirstParticleCut(dtc5etaphitpc[aniter]);
					    anetaphitpc[aniter]->SetSecondParticleCut(dtc5etaphitpc[aniter]);
					    anetaphitpc[aniter]->SetPairCut(sqp1cetaphitpc[aniter]);
					  }
					if(ichg == 18 || ichg == 25 || ichg == 32) //VOLAL
					  {
					    anetaphitpc[aniter]->SetV0SharedDaughterCut(performSharedDaughterCut);
					    anetaphitpc[aniter]->SetEventCut(mecetaphitpc[aniter]);
					    anetaphitpc[aniter]->SetFirstParticleCut(dtc4etaphitpc[aniter]);
					    anetaphitpc[aniter]->SetSecondParticleCut(dtc5etaphitpc[aniter]);
					    anetaphitpc[aniter]->SetPairCut(sqp1cetaphitpc[aniter]);
					  }
					if(ichg == 13 || ichg == 20 || ichg == 27) //V0PL
					  {
					    anetaphitpc[aniter]->SetV0SharedDaughterCut(performSharedDaughterCut);
					    anetaphitpc[aniter]->SetEventCut(mecetaphitpc[aniter]);
					    anetaphitpc[aniter]->SetFirstParticleCut(dtc4etaphitpc[aniter]);
					    anetaphitpc[aniter]->SetSecondParticleCut(dtc1etaphitpc[aniter]);
					    anetaphitpc[aniter]->SetPairCut(sqp2cetaphitpc[aniter]);
					  }
					if(ichg == 15 || ichg == 22 || ichg == 29) //V0APL
					  {
					     anetaphitpc[aniter]->SetV0SharedDaughterCut(performSharedDaughterCut);
					    anetaphitpc[aniter]->SetEventCut(mecetaphitpc[aniter]);
					    anetaphitpc[aniter]->SetFirstParticleCut(dtc4etaphitpc[aniter]);
					    anetaphitpc[aniter]->SetSecondParticleCut(dtc2etaphitpc[aniter]);
					    anetaphitpc[aniter]->SetPairCut(sqp4cetaphitpc[aniter]);
					  }
					if(ichg == 14 || ichg == 21 || ichg == 28) //V0PAL
					  {
					    anetaphitpc[aniter]->SetV0SharedDaughterCut(performSharedDaughterCut);
					     anetaphitpc[aniter]->SetEventCut(mecetaphitpc[aniter]);
					    anetaphitpc[aniter]->SetFirstParticleCut(dtc5etaphitpc[aniter]);
					    anetaphitpc[aniter]->SetSecondParticleCut(dtc1etaphitpc[aniter]);
					    anetaphitpc[aniter]->SetPairCut(sqp4cetaphitpc[aniter]);
            }
					if(ichg == 16 || ichg == 23 || ichg == 30) //V0APAL
					  {
					    anetaphitpc[aniter]->SetV0SharedDaughterCut(performSharedDaughterCut);
					    anetaphitpc[aniter]->SetEventCut(mecetaphitpc[aniter]);
					    anetaphitpc[aniter]->SetFirstParticleCut(dtc5etaphitpc[aniter]);
					    anetaphitpc[aniter]->SetSecondParticleCut(dtc2etaphitpc[aniter]);
					    anetaphitpc[aniter]->SetPairCut(sqp3cetaphitpc[aniter]);
					   }
					if(ichg == 34) //PXim
					  {
					    anetaphitpc[aniter]->SetFirstParticleCut(tXiCut[aniter]);
					    anetaphitpc[aniter]->SetSecondParticleCut(dtc1etaphitpc[aniter]);
					    anetaphitpc[aniter]->SetPairCut(tXiTrackPairCut[aniter]);
					  }
					if(ichg == 35) //aPXim
					  {
					    anetaphitpc[aniter]->SetFirstParticleCut(tXiCut[aniter]);
					    anetaphitpc[aniter]->SetSecondParticleCut(dtc2etaphitpc[aniter]);
					    anetaphitpc[aniter]->SetPairCut(tXiTrackPairCut[aniter]);
					  }
					if(ichg == 36) //PXip
					  {
					    anetaphitpc[aniter]->SetFirstParticleCut(tAXiCut[aniter]);
					    anetaphitpc[aniter]->SetSecondParticleCut(dtc1etaphitpc[aniter]);
					    anetaphitpc[aniter]->SetPairCut(tXiTrackPairCut[aniter]);
					  }
					if(ichg == 37) //aPXip
					  {
					    anetaphitpc[aniter]->SetFirstParticleCut(tAXiCut[aniter]);
					    anetaphitpc[aniter]->SetSecondParticleCut(dtc2etaphitpc[aniter]);
					    anetaphitpc[aniter]->SetPairCut(tXiTrackPairCut[aniter]);
					  }
            if(ichg == 38) //POm
              {
                anetaphitpc[aniter]->SetFirstParticleCut(tOmegaCut[aniter]);
                anetaphitpc[aniter]->SetSecondParticleCut(dtc1etaphitpc[aniter]);
                anetaphitpc[aniter]->SetPairCut(tOmegaTrackPairCut[aniter]);
              }
            if(ichg == 39) //aPOm
              {
                anetaphitpc[aniter]->SetFirstParticleCut(tOmegaCut[aniter]);
                anetaphitpc[aniter]->SetSecondParticleCut(dtc2etaphitpc[aniter]);
                anetaphitpc[aniter]->SetPairCut(tOmegaTrackPairCut[aniter]);
              }
            if(ichg == 40) //POp
              {
                anetaphitpc[aniter]->SetFirstParticleCut(tAOmegaCut[aniter]);
                anetaphitpc[aniter]->SetSecondParticleCut(dtc1etaphitpc[aniter]);
                anetaphitpc[aniter]->SetPairCut(tOmegaTrackPairCut[aniter]);
              }
            if(ichg == 41) //aPOp
              {
                anetaphitpc[aniter]->SetFirstParticleCut(tAOmegaCut[aniter]);
                anetaphitpc[aniter]->SetSecondParticleCut(dtc2etaphitpc[aniter]);
                anetaphitpc[aniter]->SetPairCut(tOmegaTrackPairCut[aniter]);
              }


					//**** Correlation functions *******
					//***without corrections*****
					//if(ichg >=38)
					  //cdedpetaphinocorr[aniter] = new AliFemtoCorrFctnDEtaDPhiSimple(Form("cdedpnocorr%stpcM%i", chrgs[ichg], imult),11, 11);
				 if(ichg >= 13 || ichg < 42)
					  cdedpetaphinocorr[aniter] = new AliFemtoCorrFctnDEtaDPhiSimple(Form("cdedpnocorr%stpcM%i", chrgs[ichg], imult),23, 23);
					else
					  cdedpetaphinocorr[aniter] = new AliFemtoCorrFctnDEtaDPhiSimple(Form("cdedpnocorr%stpcM%i", chrgs[ichg], imult),29, 29);

					anetaphitpc[aniter]->AddCorrFctn(cdedpetaphinocorr[aniter]);

					if(ichg==0 || ichg==1 || ichg==31 || ichg==33) //PP, aPaP, LL, ALAL
                                        {
					  cqinvtpc[aniter] = new AliFemtoQinvCorrFctn(Form("cqinv%stpcM%i", chrgs[ichg], imult),nbinssh,0.0,shqmax); //femto qinv, for identical mass particles
					  anetaphitpc[aniter]->AddCorrFctn(cqinvtpc[aniter]);
					}
					if(ichg==2 || ichg==27 || ichg==28 || ichg==29 || ichg==30 || ichg==32 || ichg==34 || ichg == 35 || ichg == 36 || ichg == 37) //PaP, PL, APL, PAL, APAL, LAL, PXim, aPXim, PXip, aPXip
                                        {
					  cnonidtpc[aniter] = new AliFemtoCorrFctnNonIdDR(Form("cnonid%stpcM%i", chrgs[ichg], imult), nbinssh, 0.0,shqmax); //for non-identical partcles
					  anetaphitpc[aniter]->AddCorrFctn(cnonidtpc[aniter]);
                                        }


					//***with corrections****
					/*
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
					*/

					Manager->AddAnalysis(anetaphitpc[aniter]);
				}
			}
		}
	}

	return Manager;
}
