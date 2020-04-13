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
#include "AliFemtoCorrFctnDEtaDPhiSimple"

#include "AliFemtoV0PairCut.h"
#include "AliFemtoV0TrackPairCut.h"
#include "AliFemtoV0TrackCut.h"

#include "AliFemtoCorrFctnInvMass.h"
#endif

//________________________________________________________________________
AliFemtoManager* ConfigFemtoAnalysis(const char* params) {

	double PionMass = 0.13956995;
	double KaonMass = 0.493677;
	double ProtonMass = 0.938272013;
	double LambdaMass = 1.115683;
	double XiMass = 1.32171;

	
	const int numOfMultBins = 10;	
	const int numOfChTypes = 3;
	const int numOfkTbins = 6;

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
	    cout<<"Parameter [15]: (SetMostProbable)"<<parameter[15]<<" "<<endl;
	    parameter[16] = strtok(NULL, ",");
	    cout<<"Parameter [16]: (monitors)"<<parameter[16]<<" "<<endl;
	    parameter[17] = strtok(NULL, ",");
	    cout<<"Parameter [17]: (nSigma2)"<<parameter[17]<<" "<<endl;
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

	int setMostProb = atoi(parameter[15]);

	Bool_t ifMonitors=kFALSE; if(atoi(parameter[16]))ifMonitors=kTRUE;//kTRUE 
	double nSigmaVal2 = atof(parameter[17]); //2.0 or 3.0

	printf("*** Connect to AliEn ***\n");
	TGrid::Connect("alien://");

	int runmults[numOfMultBins] = { 1, 1, 0, 0, 0, 0, 0, 0, 0, 0};
	if(runmultdep)	  {runmults[0]=1; runmults[1]=1; runmults[2]=1;	  }
	int multbins[numOfMultBins+1] = {2, 11, 16, 22, 29, 36, 44, 57, 150,2,150};

	int runch[numOfChTypes] = { 1, 1, 1};
	const char *chrgs[numOfChTypes] = { "PIpPIp", "PImPIm", "PIpPIm"};
	

	
	
	double ktrng[numOfkTbins+1] = {0.13, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7};


	int runqinv = 1;
	int runshlcms = 1;// 0:PRF(PAP), 1:LCMS(PP,APAP)
	int run3d = 0;

	int runtype = 0; // Types 0 - global, 1 - ITS only, 2 - TPC Inner	//global tracks ->mfit ITS+TPC
	int owncuts = 0; 
	int owndca = 0;

	int gammacut = 1;	// cut na ee z gamma 
	
	double shqmax = 5.0; 
	int nbinssh = 500;

	//AliFemtoEventReaderESDChain *Reader = new AliFemtoEventReaderESDChain();
	//Reader->SetUseMultiplicity(AliFemtoEventReaderESDChain::kGlobalCount);


	AliFemtoEventReaderAODChain *Reader = new AliFemtoEventReaderAODChain();
	Reader->SetFilterMask(filterbit);
	Reader->SetDCAglobalTrack(ifGlobalTracks); //false for FB7, true for the rest //we do not use DCA at all
	Reader->SetUseMultiplicity(AliFemtoEventReaderAOD::kReference);
	Reader->SetMinPlpContribSPD(minPlpContribSPD);
	Reader->SetIsPileUpEvent(ifIsPileUp);
	Reader->SetReadV0(kTRUE);
	Reader->SetReadCascade(kTRUE);
	
	Reader->SetUseMVPlpSelection(kTRUE);
	Reader->SetTrackPileUpRemoval(kTRUE);
	
	Reader->SetUseAliEventCuts(kTRUE);

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
	AliFemtoPairCutAntiGamma	*sqpcetaphitpc[numOfMultBins*numOfChTypes];
	//AliFemtoPairCutRadialDistance			*sqpcetaphitpc[numOfMultBins*numOfChTypes];
	//AliFemtoShareQualityPairCut			*sqpcetaphitpc[numOfMultBins*numOfChTypes];


	
	//	AliFemtoChi2CorrFctn					*cchiqinvetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoKTPairCut               *ktpcuts[numOfMultBins*numOfChTypes*numOfkTbins];
	AliFemtoQinvCorrFctn		*cqinvkttpc[numOfMultBins*numOfChTypes*numOfkTbins];
	AliFemtoQinvCorrFctn		*cqinvtpc[numOfMultBins*numOfChTypes];
	AliFemtoCorrFctnDirectYlm *cylmetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoCorrFctnDirectYlm * cylmkttpc[numOfMultBins*numOfChTypes*numOfkTbins];
	AliFemtoCorrFctn3DLCMSSym *cq3dlcmskttpc[numOfMultBins*numOfChTypes];
	AliFemtoCorrFctnDEtaDPhiSimple	*cdedpetaphinocorr[numOfMultBins*numOfChTypes];
	  

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
					dtc1etaphitpc[aniter]->SetPt(0.2,maxPt);
					dtc1etaphitpc[aniter]->SetMass(PionMass);		
					dtc1etaphitpc[aniter]->SetMostProbable(setMostProb);//cut on Nsigma in pT not p			   
			
	      
					dtc2etaphitpc[aniter] = new AliFemtoMJTrackCut();
					dtc2etaphitpc[aniter]->SetCharge(-1.0);
					dtc2etaphitpc[aniter]->SetEta(nEtaMin,nEtaMax);
					dtc2etaphitpc[aniter]->SetNsigma(nSigmaVal);
					dtc2etaphitpc[aniter]->SetNsigma2(nSigmaVal2);
					dtc2etaphitpc[aniter]->SetNsigmaTPCTOF(kTRUE);
					dtc2etaphitpc[aniter]->SetElectronRejection(ifElectronRejection);
					dtc2etaphitpc[aniter]->SetPt(0.2,maxPt);
					dtc2etaphitpc[aniter]->SetMass(PionMass);		
					dtc2etaphitpc[aniter]->SetMostProbable(setMostProb);//cut on Nsigma in pT not p			   
		

	

					

					//****** DCA ******

					if(owndca){
					  dtc1etaphitpc[aniter]->SetMaxImpactXYPtDep(0.018, 0.035, -1.01); 	//	DCA xy
					  //dtc1etaphitpc[aniter]->SetMaxImpactXYPtDep(0.0182, 0.0350, -1.01);
					  dtc1etaphitpc[aniter]->SetMaxImpactZ(2);	//DCA Z
					  dtc2etaphitpc[aniter]->SetMaxImpactXYPtDep(0.018, 0.035, -1.01); 	//	DCA xy
					  dtc2etaphitpc[aniter]->SetMaxImpactZ(2);	//DCA Z
					
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
				
	
					}
					//**************** track Monitors ***************

					if(ifMonitors)
					  {

					    cutPass1YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutPass%stpcM%i", chrgs[ichg], imult),PionMass);
					    cutFail1YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutFail%stpcM%i", chrgs[ichg], imult),PionMass);

					    cutPass2YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutPass%stpcM%i", chrgs[ichg], imult),PionMass);
					    cutFail2YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutFail%stpcM%i", chrgs[ichg], imult),PionMass);

					    cutPass1PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cutPass%stpcM%i", chrgs[ichg], imult),0);//0-pion,1-kaon,2-proton
					    cutFail1PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cutFail%stpcM%i", chrgs[ichg], imult),0);

					    cutPass2PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cutPass%stpcM%i", chrgs[ichg], imult),0);//0-pion,1-kaon,2-proton
					    cutFail2PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cutFail%stpcM%i", chrgs[ichg], imult),0);

					    dtc1etaphitpc[aniter]->AddCutMonitor(cutPass1YPtetaphitpc[aniter], cutFail1YPtetaphitpc[aniter]);
					    dtc1etaphitpc[aniter]->AddCutMonitor(cutPass1PIDetaphitpc[aniter], cutFail1PIDetaphitpc[aniter]);

					    dtc2etaphitpc[aniter]->AddCutMonitor(cutPass2YPtetaphitpc[aniter], cutFail2YPtetaphitpc[aniter]);
					    dtc2etaphitpc[aniter]->AddCutMonitor(cutPass2PIDetaphitpc[aniter], cutFail2PIDetaphitpc[aniter]);


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
	
					
					//***** Setting cuts ***********
					// setting event cut
					anetaphitpc[aniter]->SetEventCut(mecetaphitpc[aniter]);
					//setting single track cuts
					if(ichg==0) //positive like-sign
					{
					  anetaphitpc[aniter]->SetPairCut(sqpcetaphitpc[aniter]);
					  anetaphitpc[aniter]->SetFirstParticleCut(dtc1etaphitpc[aniter]);
					  anetaphitpc[aniter]->SetSecondParticleCut(dtc1etaphitpc[aniter]);
					}
					if( ichg==1 )//negative like-sign
					{
					  anetaphitpc[aniter]->SetPairCut(sqpcetaphitpc[aniter]);
					  anetaphitpc[aniter]->SetFirstParticleCut(dtc2etaphitpc[aniter]);
					  anetaphitpc[aniter]->SetSecondParticleCut(dtc2etaphitpc[aniter]);
					}
					if(ichg==2)//unlike-sign
					{
					  anetaphitpc[aniter]->SetPairCut(sqpcetaphitpc[aniter]);
					  anetaphitpc[aniter]->SetFirstParticleCut(dtc1etaphitpc[aniter]);
					  anetaphitpc[aniter]->SetSecondParticleCut(dtc2etaphitpc[aniter]);
					}
				

					//**** Correlation functions *******	
					cqinvtpc[aniter] = new AliFemtoQinvCorrFctn(Form("cqinv%stpcM%i", chrgs[ichg], imult),nbinssh,0.0,shqmax); //femto qinv, for identical mass particles
					anetaphitpc[aniter]->AddCorrFctn(cqinvtpc[aniter]);

					 cdedpetaphinocorr[aniter] = new AliFemtoCorrFctnDEtaDPhiSimple(Form("cdedpnocorr%stpcM%i", chrgs[ichg], imult),29, 29);
					 anetaphitpc[aniter]->AddCorrFctn(cdedpetaphinocorr[aniter]);
				
					//Spherical harmonics (without kT bins)
					//cylmetaphitpc[aniter] = new AliFemtoCorrFctnDirectYlm(Form("cylm%stpcM%i", chrgs[ichg], imult),3,nbinssh,0.0,shqmax,runshlcms);
					//anetaphitpc[aniter]->AddCorrFctn(cylmetaphitpc[aniter]);

					//if(run3d){
					//cq3dlcmskttpc[aniter] = new AliFemtoCorrFctn3DLCMSSym(Form("cq3d%stpcM%i", chrgs[ichg], imult),60,0.5);
					//anetaphitpc[aniter]->AddCorrFctn(cq3dlcmskttpc[aniter]);
					//}

					if (runktdep) {
					  int ktm;
					  for (int ikt=0; ikt<6; ikt++) {
					    ktm = aniter*6 + ikt;
					    ktpcuts[ktm] = new AliFemtoKTPairCut(ktrng[ikt], ktrng[ikt+1]);
					    
					    //cylmkttpc[ktm] = new AliFemtoCorrFctnDirectYlm(Form("cylm%stpcM%ikT%i", chrgs[ichg], imult, ikt),3, nbinssh, 0.0, shqmax, runshlcms);
					    //cylmkttpc[ktm]->SetPairSelectionCut(ktpcuts[ktm]);
					    //anetaphitpc[aniter]->AddCorrFctn(cylmkttpc[ktm]);
					    
					    cqinvkttpc[ktm] = new AliFemtoQinvCorrFctn(Form("cqinv%stpcM%ikT%i", chrgs[ichg], imult, ikt),nbinssh,0.0, shqmax);
					    cqinvkttpc[ktm]->SetPairSelectionCut(ktpcuts[ktm]);
					    anetaphitpc[aniter]->AddCorrFctn(cqinvkttpc[ktm]);
					    
					  }
					}

						       							
					Manager->AddAnalysis(anetaphitpc[aniter]);
				}
			}
		}
	}
					    
	return Manager;
}								 
