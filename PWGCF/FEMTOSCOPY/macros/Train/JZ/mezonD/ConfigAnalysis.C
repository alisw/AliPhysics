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
#include "AliFemtoEventReaderAODMultSelection.h" //zmiana
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
//#include "AliFemtoCorrFctnDPhiStarDEta.h"
#include "AliFemtoCorrFctnDEtaDPhiStar.h"
#include "AliFemtoCorrFctnInvMass.h"

#include "AliFemtoV0PairCut.h"
#include "AliFemtoV0TrackPairCut.h"
#include "AliFemtoV0TrackCut.h"
#endif

#include <stdio.h>
#include <string.h>

//________________________________________________________________________
AliFemtoManager* ConfigFemtoAnalysis(const char* params) {

	double PionMass = 0.13956995;
	double KaonMass = 0.493677;
	double D0Mass = 1.86483;

	const int numOfMultBins = 1;
	const int numOfChTypes = 4;
	const int numOfkTbins = 1;

	bool performSharedDaughterCut = false;
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

	int setMostProb1 = atoi(parameter[15]); //21 (protons)
	int setMostProb2 = atoi(parameter[16]); //20 (kaons)
	int setMostProb3 = atoi(parameter[17]); //19 (pions)

	Bool_t ifMonitors=kFALSE; if(atoi(parameter[19]))ifMonitors=kTRUE;//kTRUE
	double nSigmaVal2 = atof(parameter[20]); //3.0 (or 2.0)

	printf("*** Connect to AliEn ***\n");
	TGrid::Connect("alien://");

	int runmults[numOfMultBins] = {1};
	int multbins[numOfMultBins+1] = {0, 1000};

	int runch[numOfChTypes] = {1, 1, 1, 1};
	const char *chrgs[numOfChTypes] = { "KpPi", "KmPi", "KPip", "KPim"};

	int runqinv = 1;

	int runtype = 0; // Types 0 - global, 1 - ITS only, 2 - TPC Inner	//global tracks ->mfit ITS+TPC
	int owncuts = 1;
	int owndca = 1;

	int gammacut = 1;	// cut na ee z gamma

	double shqmax = 2.0;
	int nbinssh = 20;

	//AliFemtoEventReaderESDChain *Reader = new AliFemtoEventReaderESDChain();
	//Reader->SetUseMultiplicity(AliFemtoEventReaderESDChain::kGlobalCount);


	AliFemtoEventReaderAODMultSelection  *Reader = new AliFemtoEventReaderAODMultSelection(); //tutaj byla zmiana
	Reader->SetFilterMask(filterbit);
	Reader->SetDCAglobalTrack(ifGlobalTracks); //false for FB7, true for the rest //we do not use DCA at all
	Reader->SetUseMultiplicity(AliFemtoEventReaderAOD::kCentrality); //zmiana z  kReference
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
	AliFemtoCutMonitorParticleYPt   *cutPass1YPtetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoCutMonitorParticleYPt   *cutFail1YPtetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoCutMonitorParticlePID   *cutPass1PIDetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoCutMonitorParticlePID   *cutFail1PIDetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoCutMonitorParticleYPt   *cutPass2YPtetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoCutMonitorParticleYPt   *cutFail2YPtetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoCutMonitorParticlePID   *cutPass2PIDetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoCutMonitorParticlePID   *cutFail2PIDetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoPairCutRadialDistance			*sqpcetaphitpc[numOfMultBins*numOfChTypes];
	AliFemtoPairCutPt               *ktpcuts[numOfMultBins*numOfkTbins];
	AliFemtoCorrFctnInvMass					*fInvMass[numOfMultBins*numOfChTypes];



	std::cout<<numOfMultBins<<" "<<numOfChTypes<<std::endl;
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
					mecetaphitpc[aniter]->SetVertZPos(-7,7);//cm

					//****** event monitors **********
					cutPassEvMetaphitpc[aniter] = new AliFemtoCutMonitorEventMult(Form("cutPass%stpcM%i", "KPi", imult), 2000, 20000.5);
					cutFailEvMetaphitpc[aniter] = new AliFemtoCutMonitorEventMult(Form("cutFail%stpcM%i", "KPi", imult), 2000, 20000.5);
					mecetaphitpc[aniter]->AddCutMonitor(cutPassEvMetaphitpc[aniter], cutFailEvMetaphitpc[aniter]);


					// ***** single particle track cuts *********

					// ***** First particle (kaons)
					dtc1etaphitpc[aniter] = new AliFemtoMJTrackCut();
					if (ichg == 0)
					{
						dtc1etaphitpc[aniter]->SetCharge(1.0);
					}
					else if (ichg == 1)
					{
						dtc1etaphitpc[aniter]->SetCharge(-1.0);
					}
					if (ichg == 2 || ichg == 3)
					{
						dtc1etaphitpc[aniter]->SetCharge(0);
					}
					dtc1etaphitpc[aniter]->SetEta(nEtaMin,nEtaMax);
					dtc1etaphitpc[aniter]->SetNsigma(nSigmaVal); //w zaleznosci od pedow jest brane albo nsigma1 albo nsigma2
					dtc1etaphitpc[aniter]->SetNsigma2(nSigmaVal2);

					dtc1etaphitpc[aniter]->SetNsigmaTPCTOF(kTRUE);//bierzemy jednoczesnie sigma z TPC i TOF jednoczesnie
					dtc1etaphitpc[aniter]->SetElectronRejection(ifElectronRejection);

					dtc1etaphitpc[aniter]->SetPt(0.3,maxPt);
					dtc1etaphitpc[aniter]->SetMass(KaonMass);
					dtc1etaphitpc[aniter]->SetMostProbable(setMostProb2);//cut on Nsigma in pT not p


					// ***** Second particle (pions)
					dtc2etaphitpc[aniter] = new AliFemtoMJTrackCut();
					if (ichg == 0 | ichg == 1)
					{
						dtc2etaphitpc[aniter]->SetCharge(0);
					}
					else if (ichg == 2)
					{
						dtc2etaphitpc[aniter]->SetCharge(1.0);
					}
					if (ichg == 2)
					{
						dtc2etaphitpc[aniter]->SetCharge(-1.0);
					}

					dtc2etaphitpc[aniter]->SetEta(nEtaMin,nEtaMax);
					dtc2etaphitpc[aniter]->SetNsigma(nSigmaVal); //ustawia PID
					dtc2etaphitpc[aniter]->SetNsigma2(nSigmaVal2);//
					dtc2etaphitpc[aniter]->SetNsigmaTPCTOF(kTRUE);//
					dtc2etaphitpc[aniter]->SetElectronRejection(ifElectronRejection);
					dtc2etaphitpc[aniter]->SetPt(0.2,maxPt);
					dtc2etaphitpc[aniter]->SetMass(PionMass);
					dtc2etaphitpc[aniter]->SetMostProbable(setMostProb3);//cut on Nsigma in pT not p


					//****** DCA ******

					if(owndca){
					  dtc1etaphitpc[aniter]->SetMaxImpactXYPtDep(1, 0, 0); 	//	DCA xy - parametry zalezno?ci DCA od pt
					  //dtc1etaphitpc[aniter]->SetMaxImpactXYPtDep(0.0182, 0.0350, -1.01);
					  dtc1etaphitpc[aniter]->SetMaxImpactZPtDep(1, 0, -1.01);	//DCA Z
					  dtc2etaphitpc[aniter]->SetMaxImpactXYPtDep(1, 0, 0); 	//	DCA xy
					  dtc2etaphitpc[aniter]->SetMaxImpactZPtDep(1, 0, -1.01);	//DCA Z

						//****** Track quality cuts ******
				  	if(owncuts){
				  		//dtc1etaphitpc[aniter]->SetStatus(AliESDtrack::kTPCrefit|AliESDtrack::kITSrefit);
							dtc1etaphitpc[aniter]->SetminTPCncls(70);
							dtc1etaphitpc[aniter]->SetRemoveKinks(kTRUE);
							dtc1etaphitpc[aniter]->SetLabel(kFALSE);
							//	dtc1etaphitpc[aniter]->SetMaxITSChiNdof(6.0);
							dtc1etaphitpc[aniter]->SetMaxTPCChiNdof(4.0);	// pisac
							dtc1etaphitpc[aniter]->SetMaxSigmaToVertex(6.0);

							//dtc2etaphitpc[aniter]->SetStatus(AliESDtrack::kTPCrefit|AliESDtrack::kITSrefit);
							dtc2etaphitpc[aniter]->SetminTPCncls(70);
							dtc2etaphitpc[aniter]->SetRemoveKinks(kTRUE); //wyrzucamy tracki "?amane"
							dtc2etaphitpc[aniter]->SetLabel(kFALSE); //Label - numer czastki
							//	dtc2etaphitpc[aniter]->SetMaxITSChiNdof(6.0);
							dtc2etaphitpc[aniter]->SetMaxTPCChiNdof(4.0);
							dtc2etaphitpc[aniter]->SetMaxSigmaToVertex(6.0);
						}
					}
					//**************** track Monitors ***************

					if(ifMonitors)//ichg>8)
					{
				//monitor rysuje tylko histogramy!!
						cutPass1PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cutPass%stpcM%i", "K", imult),1);//0-pion,1-kaon,2-proton
						cutFail1PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cutFail%stpcM%i", "K", imult),1);
						dtc1etaphitpc[aniter]->AddCutMonitor(cutPass1PIDetaphitpc[aniter], cutFail1PIDetaphitpc[aniter]);

						cutPass2PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cutPass%stpcM%i", "Pi", imult),0);//0-pion,1-kaon,2-proton
						cutFail2PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cutFail%stpcM%i", "Pi", imult),0);
						dtc2etaphitpc[aniter]->AddCutMonitor(cutPass2PIDetaphitpc[aniter], cutFail2PIDetaphitpc[aniter]);
					}

					sqpcetaphitpc[aniter] = new AliFemtoPairCutRadialDistance();
				  sqpcetaphitpc[aniter]->SetRemoveSameLabel(kFALSE);
				  sqpcetaphitpc[aniter]->SetMinimumRadius(1.6);
				  sqpcetaphitpc[aniter]->SetEtaDifferenceMinimum(0.02);
				  sqpcetaphitpc[aniter]->SetPhiStarDifferenceMinimum(0.045);


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


					//***** Setting cuts ***********
					// setting event cut
					anetaphitpc[aniter]->SetEventCut(mecetaphitpc[aniter]);
					//setting single track cuts
					anetaphitpc[aniter]->SetPairCut(sqpcetaphitpc[aniter]);
					anetaphitpc[aniter]->SetFirstParticleCut(dtc1etaphitpc[aniter]);
					anetaphitpc[aniter]->SetSecondParticleCut(dtc2etaphitpc[aniter]);



					//**** Correlation functions *******
					fInvMass[aniter] = new AliFemtoCorrFctnInvMass(Form("%s",chrgs[ichg]), 100, 1.750,1.950,KaonMass,PionMass);
					anetaphitpc[aniter]->AddCorrFctn(fInvMass[aniter]);

					//*** calculating invariant mass *****
					Manager->AddAnalysis(anetaphitpc[aniter]);
				}
			}
		}
	}

	return Manager;
}
