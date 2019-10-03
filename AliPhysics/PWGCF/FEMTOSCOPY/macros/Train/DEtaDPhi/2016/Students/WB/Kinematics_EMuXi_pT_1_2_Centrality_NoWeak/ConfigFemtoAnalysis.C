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
AliFemtoManager* ConfigFemtoAnalysis(const char* params) {

	double PionMass = 0.13956995; //masy
	double KaonMass = 0.493677;
	double ProtonMass = 0.9382720813;
	double LambdaMass = 1.115683;
   
	const int numOfMultBins = 1; //biny w krotniosci
	const int numOfChTypes = 13; // liczba mozliwych par
	const int numOfpTbins = 1; // biny w pedzie poprzecznym

   int runch[numOfChTypes] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0}; // 1 - wlacza czastki do analizy: wlaczone PIpPIm , all
   const char *chrgs[numOfChTypes] = { "EE", "EpEp", "EEp", "XIpXIp", "XImXIm", "XIpXIm", "MUpMUp", "MUmMUm", "MUpMUm", "all", "plus", "minus", "mixed"};
   // pary czastek ktore analizujemy PP - proton proton; KpKp - K plus
   int runktdep = 0; // czy chcemy wlaczac analize w zaleznosci od pedow poprzecznych; 0 - nie ma, 1 - jest
   double ktrng[numOfpTbins+1] = {0.0, 100.0}; // przedzialy sum pedow poprzecznych

   int gammacut = 1;   // cut for e+e- coming from gamma //wlaczenie cuty na konkretne pary czastek


//Ustalamy z ktora analiza mamy do czynienia
   //  ***Reader ESD***
   //AliFemtoEventReaderESDChain *Reader = new AliFemtoEventReaderESDChain();
   //Reader->SetUseMultiplicity(AliFemtoEventReaderESDChain::kReferenceITSTPC); //sposob liczenia multiplicity - ustawione zeby liczyc track po tracku, po ustawieniu konkretnych cutow
//dla olow-olow po centralnosci


   //  ***Reader AOD***		
   //AliFemtoEventReaderAODChain *Reader = new AliFemtoEventReaderAODChain();
   //Reader->SetFilterMask(96);		//wpisane DCA w filter bity				
   //Reader->SetDCAglobalTrack(kTRUE); //false for FB7, true for the rest //we do not use DCA (distance of closest approach - na przyklad zeby sprawdzic czy 2 czastki pochodza od jednej czastki - matki, primary vertex - punkt zderzenia - z niego pochodza wszystkie pierwotne czastki) at all
   //Reader->SetUseMultiplicity(AliFemtoEventReaderAOD::kReference);
   //Reader->SetMinPlpContribSPD(3); //PileUp - zderzenia nachodza na siebie, ten cut dodatkowo wyrzuca zderzenia pylab, ktorych nie wyrzucaja trigery
   //Reader->SetIsPileUpEvent(kTRUE);

   //  ***Reader Kinematics (Monte Carlo data)***
   AliFemtoEventReaderKinematicsChain* Reader=new AliFemtoEventReaderKinematicsChain();
   Reader->SetUseMultiplicity(AliFemtoEventReaderKinematicsChain::kImpactParameter);
   Reader->RemoveWeakDecaysManually(kTRUE);


   AliFemtoManager* Manager = new AliFemtoManager(); //typowy dla alifemto
   Manager->SetEventReader(Reader);
	// * Monitory * - histogramy ktore beda wypelniane podczas analizy - wlasciwie klasy generujace histogramy
   AliFemtoVertexMultAnalysis      *anetaphitpc[numOfMultBins*numOfChTypes]; //korelujemy eventy tylko podobne w wertexie i multiplicity
   AliFemtoBasicEventCut             *mecetaphitpc[numOfMultBins*numOfChTypes]; //podstawowy cut na zderzenia
   AliFemtoCutMonitorEventMult    *cutPassEvMetaphitpc[numOfMultBins*numOfChTypes];
   AliFemtoCutMonitorEventMult    *cutFailEvMetaphitpc[numOfMultBins*numOfChTypes];
   AliFemtoCutMonitorEventVertex *cutPassEvVetaphitpc[numOfMultBins*numOfChTypes];
   AliFemtoCutMonitorEventVertex *cutFailEvVetaphitpc[numOfMultBins*numOfChTypes];
   AliFemtoMCTrackCut          *dtc1etaphitpc[numOfMultBins*numOfChTypes]; // ograniczenia dla konkretnych czastek
   AliFemtoMCTrackCut          *dtc2etaphitpc[numOfMultBins*numOfChTypes];
   AliFemtoMCTrackCut          *dtc3etaphitpc[numOfMultBins*numOfChTypes];
   //AliFemtoCutMonitorParticleYPt *cutPass1YPtetaphitpc[320];
   //AliFemtoCutMonitorParticleYPt *cutFail1YPtetaphitpc[320];
   //AliFemtoCutMonitorParticlePID *cutPass1PIDetaphitpc[numOfMultBins*numOfChTypes];
   //AliFemtoCutMonitorParticlePID *cutFail1PIDetaphitpc[numOfMultBins*numOfChTypes];
   //AliFemtoCutMonitorParticleYPt *cutPass2YPtetaphitpc[320];
   //AliFemtoCutMonitorParticleYPt *cutFail2YPtetaphitpc[320];
   //AliFemtoCutMonitorParticlePID *cutPass2PIDetaphitpc[320];
   //AliFemtoCutMonitorParticlePID *cutFail2PIDetaphitpc[320];
   AliFemtoCutMonitorParticleYPt *cutPass3YPtetaphitpc[numOfMultBins*numOfChTypes];
   AliFemtoCutMonitorParticleYPt *cutFail3YPtetaphitpc[numOfMultBins*numOfChTypes];
   //AliFemtoCutMonitorParticlePID *cutPass3PIDetaphitpc[320];
   //AliFemtoCutMonitorParticlePID *cutFail3PIDetaphitpc[320];
   //    AliFemtoShareQualityTPCEntranceSepPairCut         *sqpcetaphitpcsame[320]; //sprawdza jak duzo trackow jest dzielone przez czastki
   AliFemtoPairCutAntiGamma         *sqpcetaphitpc[numOfMultBins*numOfChTypes]; //ograniczenia na pary czastek
   //AliFemtoPairCutRadialDistance         *sqpcetaphitpc[320];
   //   AliFemtoChi2CorrFctn               *cchiqinvetaphitpc[320];
   AliFemtoPairCutPt                   *ktpcuts[numOfMultBins*numOfChTypes*numOfpTbins];
   AliFemtoQinvCorrFctn               *cqinvkttpc[numOfMultBins*numOfChTypes*numOfpTbins];
   AliFemtoQinvCorrFctn               *cqinvtpc[numOfMultBins*numOfChTypes*numOfpTbins];
   AliFemtoCorrFctnDEtaDPhi         *cdedpetaphi[numOfMultBins*numOfChTypes*numOfpTbins];


   
   //***Delta Eta Delta Phi analysis for identified systems (chg: 0-8)***
   
   int aniter = 0;
   int runmults[numOfMultBins] = {1}; // 1 - wlaczony przedzial
   int multbins[numOfMultBins+1] = {0, 523}; // definiujemy przedzialy krotnosci // "Impact Parameter * 100" (523)
    // petla po wszystkich czastkach 
   for (int imult = 0; imult < numOfMultBins; imult++) // petla po krotnosciach
   {
      if (runmults[imult]) // dla kazdego wlaczonego multiplicity
      {
         for (int ichg = 0; ichg < numOfChTypes; ichg++)
         {
            if (runch[ichg]) //dla kazdych wlaczonych typow czastek
            {

               aniter = ichg * numOfMultBins + imult;
               //miksujemy tylko po podobnych krotnosciach
               int multmix = 5;
              
               anetaphitpc[aniter] = new AliFemtoVertexMultAnalysis(10, -10.0, 10.0, multmix, multbins[imult], multbins[imult+1]);
               anetaphitpc[aniter]->SetNumEventsToMix(5); //zwiekszamy statystyke w mianowkiku sygnalu korelacji
               anetaphitpc[aniter]->SetMinSizePartCollection(1); //co najmniej jedna czastka na event
               anetaphitpc[aniter]->SetVerboseMode(kFALSE);

               //*** Event cut ***
               mecetaphitpc[aniter] = new AliFemtoBasicEventCut();
               mecetaphitpc[aniter]->SetEventMult(0.001,100000); //(min, max)
               mecetaphitpc[aniter]->SetVertZPos(-10,10);//cm *akceptujemy zderenia ktore mialy miejsce w przedziale od -10cm do 10 cm od srodka detektora

               //****** event monitors **********   
               cutPassEvMetaphitpc[aniter] = new AliFemtoCutMonitorEventMult(Form("cutPass%stpcM%i", chrgs[ichg], imult), 2000, 20000.5);
               cutFailEvMetaphitpc[aniter] = new AliFemtoCutMonitorEventMult(Form("cutFail%stpcM%i", chrgs[ichg], imult), 2000, 20000.5);
               mecetaphitpc[aniter]->AddCutMonitor(cutPassEvMetaphitpc[aniter], cutFailEvMetaphitpc[aniter]);
      
               cutPassEvVetaphitpc[aniter] = new AliFemtoCutMonitorEventVertex(Form("cutPass%stpcM%i", chrgs[ichg], imult));
               cutFailEvVetaphitpc[aniter] = new AliFemtoCutMonitorEventVertex(Form("cutFail%stpcM%i", chrgs[ichg], imult));
               mecetaphitpc[aniter]->AddCutMonitor(cutPassEvVetaphitpc[aniter], cutFailEvVetaphitpc[aniter]);


               // ***** single particle track cuts *********
               dtc1etaphitpc[aniter] = new AliFemtoMCTrackCut();
               dtc2etaphitpc[aniter] = new AliFemtoMCTrackCut();
               dtc3etaphitpc[aniter] = new AliFemtoMCTrackCut();

               

               dtc1etaphitpc[aniter]->SetEta(-0.8,0.8);// psedopospiesznosc 0 - czastki wyprodukowane prostopadle do osi detektora
               dtc2etaphitpc[aniter]->SetEta(-0.8,0.8);// limit 0.9; powyzej czatka nie zostaje zarejestrowana
               dtc3etaphitpc[aniter]->SetEta(-0.8,0.8);
               

				if (ichg == 0 ||ichg == 1 ||ichg == 2)//electrons 0-2
				{
					dtc1etaphitpc[aniter]->SetCharge(-1.0);//ustawiamy ograniczenia na ladunek, zostawi dodatnie
					dtc2etaphitpc[aniter]->SetCharge(1.0);
					
					dtc1etaphitpc[aniter]->SetPt(1,2); // w zaleznosci od czastki rozne pedy
					dtc2etaphitpc[aniter]->SetPt(1,2);
											   
					dtc1etaphitpc[aniter]->SetPDG(11); //numer rodzaju czastek
					dtc2etaphitpc[aniter]->SetPDG(11);
					/*
					cutPass3YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutPass%stpcM%i", chrgs[ichg], imult),ProtonMass);
					cutFail3YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutFail%stpcM%i", chrgs[ichg], imult),ProtonMass);
					
					if(ichg==0) 
						dtc1etaphitpc[aniter]->AddCutMonitor(cutPass3YPtetaphitpc[aniter], cutFail3YPtetaphitpc[aniter]);
					if(ichg==1) 
						dtc2etaphitpc[aniter]->AddCutMonitor(cutPass3YPtetaphitpc[aniter], cutFail3YPtetaphitpc[aniter]);
					 */
				}
				if (ichg == 3 ||ichg == 4 ||ichg == 5)//xis 3-5
				{
					dtc1etaphitpc[aniter]->SetCharge(-1.0);//ustawiamy ograniczenia na ladunek, zostawi dodatnie
					dtc2etaphitpc[aniter]->SetCharge(1.0);
					
					dtc1etaphitpc[aniter]->SetPt(1,2);
					dtc2etaphitpc[aniter]->SetPt(1,2);

					dtc1etaphitpc[aniter]->SetPDG(3312);
					dtc2etaphitpc[aniter]->SetPDG(3312);
					/*
					cutPass3YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutPass%stpcM%i", chrgs[ichg], imult),KaonMass);
					cutFail3YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutFail%stpcM%i", chrgs[ichg], imult),KaonMass);
					
					if(ichg==3) 
						dtc1etaphitpc[aniter]->AddCutMonitor(cutPass3YPtetaphitpc[aniter], cutFail3YPtetaphitpc[aniter]);
					if(ichg==4) 
						dtc2etaphitpc[aniter]->AddCutMonitor(cutPass3YPtetaphitpc[aniter], cutFail3YPtetaphitpc[aniter]);
					*/
				}
				if (ichg == 6 ||ichg == 7 ||ichg == 8)//muons 6-8
				{
					dtc1etaphitpc[aniter]->SetCharge(-1.0);//ustawiamy ograniczenia na ladunek, zostawi dodatnie
					dtc2etaphitpc[aniter]->SetCharge(1.0);
					
					dtc1etaphitpc[aniter]->SetPt(1,2);
					dtc2etaphitpc[aniter]->SetPt(1,2);

					dtc1etaphitpc[aniter]->SetPDG(13);
					dtc2etaphitpc[aniter]->SetPDG(13); 
					/*
					cutPass3YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutPass%stpcM%i", chrgs[ichg], imult),PionMass);
					cutFail3YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutFail%stpcM%i", chrgs[ichg], imult),PionMass);  
				
					if(ichg==6) 
						dtc1etaphitpc[aniter]->AddCutMonitor(cutPass3YPtetaphitpc[aniter], cutFail3YPtetaphitpc[aniter]);
					if(ichg==7) 
						dtc2etaphitpc[aniter]->AddCutMonitor(cutPass3YPtetaphitpc[aniter], cutFail3YPtetaphitpc[aniter]);
					*/
				}
				
                if (ichg == 9)//all
                {
					
					dtc3etaphitpc[aniter]->SetPt(1,2);//min 0.12 TPC jeszcze wyciaga
                }
                if (ichg == 10 ||ichg == 11 ||ichg == 12)//plus,minus,mixed
                {
					dtc1etaphitpc[aniter]->SetCharge(1.0);//ustawiamy ograniczenia na ladunek, zostawi dodatnie
					dtc2etaphitpc[aniter]->SetCharge(-1.0);
					dtc1etaphitpc[aniter]->SetPt(1,2);
					dtc2etaphitpc[aniter]->SetPt(1,2);
                }
				
				
				
               //**************** track Monitors ***************
               /*if(1)//ichg>8)
                 {
                   cutPass3YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutPass%stpcM%i", chrgs[ichg], imult),PionMass);
                   cutFail3YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cutFail%stpcM%i", chrgs[ichg], imult),PionMass);
                   if(ichg==0||ichg==3||ichg==6) dtc1etaphitpc[aniter]->AddCutMonitor(cutPass3YPtetaphitpc[aniter], cutFail3YPtetaphitpc[aniter]);
                   if(ichg==1||ichg==4||ichg==7) dtc2etaphitpc[aniter]->AddCutMonitor(cutPass3YPtetaphitpc[aniter], cutFail3YPtetaphitpc[aniter]);
                 }*/
                
               //******** Two - track cuts ************
				sqpcetaphitpc[aniter] = new AliFemtoPairCutAntiGamma();
               //sqpcetaphitpc[aniter] = new AliFemtoPairCutRadialDistance();
               sqpcetaphitpc[aniter]->SetDataType(AliFemtoPairCut::kKine);
               
               
               //***** Setting cuts ***********
               // setting event cut
               anetaphitpc[aniter]->SetEventCut(mecetaphitpc[aniter]);
               //setting single track cuts
               if(ichg == 0 || ichg == 3 || ichg == 6) //positive like-sign
               {
                  anetaphitpc[aniter]->SetFirstParticleCut(dtc1etaphitpc[aniter]);
                  anetaphitpc[aniter]->SetSecondParticleCut(dtc1etaphitpc[aniter]);
               }
               if(ichg == 1 || ichg == 4 || ichg == 7)//negative like-sign
               {
                  anetaphitpc[aniter]->SetFirstParticleCut(dtc2etaphitpc[aniter]);
                  anetaphitpc[aniter]->SetSecondParticleCut(dtc2etaphitpc[aniter]);
               }
               if(ichg == 2 || ichg == 5 || ichg == 8)//unlike-sign
               {
                  anetaphitpc[aniter]->SetFirstParticleCut(dtc1etaphitpc[aniter]);
                  anetaphitpc[aniter]->SetSecondParticleCut(dtc2etaphitpc[aniter]);
               }
				if(ichg==9) //all
               {
					  anetaphitpc[aniter]->SetFirstParticleCut(dtc3etaphitpc[aniter]);
					  anetaphitpc[aniter]->SetSecondParticleCut(dtc3etaphitpc[aniter]);
               }
               if(ichg==10) //positive like-sign
					{
						anetaphitpc[aniter]->SetFirstParticleCut(dtc1etaphitpc[aniter]);
						anetaphitpc[aniter]->SetSecondParticleCut(dtc1etaphitpc[aniter]);
					}
				if(ichg==11)//negative like-sign
				{
					anetaphitpc[aniter]->SetFirstParticleCut(dtc2etaphitpc[aniter]);
					anetaphitpc[aniter]->SetSecondParticleCut(dtc2etaphitpc[aniter]);
				}
				if(ichg==12)//unlike-sign
				{
					anetaphitpc[aniter]->SetFirstParticleCut(dtc1etaphitpc[aniter]);
					anetaphitpc[aniter]->SetSecondParticleCut(dtc2etaphitpc[aniter]);
				}
               
               //setting two-track cuts
               anetaphitpc[aniter]->SetPairCut(sqpcetaphitpc[aniter]);

               //**** Correlation functions *******

               cdedpetaphi[aniter] = new AliFemtoCorrFctnDEtaDPhi(Form("DEtaDPhi_%s_M%i", chrgs[ichg], imult),35, 35);// torzenie funkcji korelacji
               anetaphitpc[aniter]->AddCorrFctn(cdedpetaphi[aniter]);
               cqinvtpc[aniter] = new AliFemtoQinvCorrFctn(Form("Qinv_%s_M%i", chrgs[ichg], imult), 20, 0, 2);
               anetaphitpc[aniter]->AddCorrFctn(cqinvtpc[aniter]);

               if (runktdep)
               {
                  int ktm;
                  for (int ikt=0; ikt<numOfpTbins; ikt++)// dla kazdego przedzialu funkcje liczone oddzielnie
                  {
                     ktm = aniter * numOfpTbins + ikt;
                     ktpcuts[ktm] = new AliFemtoPairCutPt(ktrng[ikt], ktrng[ikt+1]);

                     cdedpetaphi[ktm] = new AliFemtoCorrFctnDEtaDPhi(Form("cdedp%stpcM%ipT%i", chrgs[ichg], imult,ikt),35, 35);
                     cdedpetaphi[ktm]->SetPairSelectionCut(ktpcuts[ktm]);
                     anetaphitpc[aniter]->AddCorrFctn(cdedpetaphi[ktm]);

                  }
               }      
               Manager->AddAnalysis(anetaphitpc[aniter]);   
            }
         }
      }
   }
   return Manager;// manager glowny program wykonujacy analize
}

