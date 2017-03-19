#if !defined(__CINT__) || defined(__MAKECINT_)
#include "AliFemtoManager.h"
#include "AliFemtoEventReaderESDChain.h"
#include "AliFemtoEventReaderESDChainKine.h"
#include "AliFemtoEventReaderAODChain.h"
#include "AliFemtoSimpleAnalysis.h"
#include "AliFemtoBasicEventCut.h"
#include "AliFemtoESDTrackCut.h"
#include "AliFemtoCorrFctn.h"
#include "AliFemtoCutMonitorParticleYPtWithWeights.h"
#include "AliFemtoCutMonitorParticleVertPos.h"
#include "AliFemtoCutMonitorParticleMomRes.h"
#include "AliFemtoCutMonitorParticlePID.h"
#include "AliFemtoCutMonitorEventMult.h"
#include "AliFemtoCutMonitorEventVertex.h"
#include "AliFemtoShareQualityTPCEntranceSepPairCut.h"
#include "AliFemtoPairCutAntiGamma.h"
#include "AliFemtoPairCutRadialDistance.h"
#include "AliFemtoModelCorrFctnQinv.h"
#include "AliFemtoQinvCorrFctn.h"
#include "AliFemtoCorrFctnDEtaDPhi.h"
#include "AliFemtoModelCorrFctnWithWeights.h"
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
#include "AliFemtoModelWeightGenerator.h"
#include "AliFemtoModelWeightGeneratorBasic.h"
#include "AliFemtoModelWeightGeneratorLednicky.h"
#include "AliFemtoCorrFctnDirectYlm.h"
#include "AliFemtoModelCorrFctnDirectYlm.h"
#include "AliFemtoModelCorrFctnSource.h"
#include "AliFemtoCutMonitorParticlePtPDG.h"
#include "AliFemtoKTPairCut.h"
#include "AliFemtoPairCutPt.h"
#include "AliFemtoModelCorrFctnDEtaDPhiWithWeights.h"
#include "AliFemtoAvgSepCorrFctn.h"
#include "AliFemtoModelGausRinvFreezeOutGenerator.h"
#include "AliFemtoQinvCorrFctnWithWeights.h"
#include "AliFemtoModelCorrFctnDEtaDPhiWithWeights.h"
#include "AliFemtoModelCorrFctn.h"
#include "AliFemtoModelCorrFctnDEtaDPhiRM.h"
#include "AliFemtoCutMonitorPairKT.h"
#include "AliFemtoPairCutMInv.h"

#endif

//________________________________________________________________________
AliFemtoManager* ConfigFemtoAnalysis(const char* params) 
{
//masy czastek w GeV
   double PionMass = 0.13956995;
   double KaonMass = 0.493677;
   double ProtonMass = 0.938272013;
   double LambdaMass = 1.115683;
   double RhoMass = 0.77549;
   double PhiMass = 1.019445;

   //double *part_mases[4] = {ProtonMass, KaonMass, PionMass, LambdaMass};

/************** MULTIPLICITY *************/   
   const int numOfMultBins = 5;
   int runmults[numOfMultBins] = {1, 0, 0, 0, 0};    // 1 wlaczony przedzial krotnosci, 0 wylaczony
   int multbins[numOfMultBins+1] = {2, 20000, 0, 5000, 10001, 15000}; //definiujemy przedzialy krotnosci 
/************** OTHER SETTINGS *************/   
   int gammacut = 1;   // cut for e+e- coming from gamma
   int runktdep = 1;
   int resonancecut = 0;
/************** PARTICLE TYPES FOR ANALYSIS *************/   
	//ilosc par dla ktorych liczymy korelacje
   const int numOfPairTypes = 16;
	//typy analizowanych par czastek, all - bierzemy (dowolna, dowolna), plus - (dodatnia, dodatnia), mixed (plus, minus) Ostateczna ilosc wykresow = liczba par * liczba binow multiplicity * [liczba binow pedow + 1] (+1 bo chcemy zazwyczaj wykres bez binowania)
  const char* filterTypes[8] = {"p+", "p-","K+", "K-", "pi+", "pi-", "lambda", "lambdab"};
   const char *pair_types[numOfPairTypes] = { "PP" /*0*/, "aPaP", "PaP", "KpKp", "KmKm", "KpKm", "PIpPIp", "PImPIm", "PIpPIm", \
   "LamLam" /*9*/, "aLam_aLam", "Lam_aLam", "all" /*12*/, "plus", "minus", "mixed"};
	// 1 wlaczaja nam pary do analizy, 0 wylaczaja
  int runch[numOfPairTypes]= {0/*PP*/, 0/*aPaP*/, 0/*PaP*/, 0/*KpKp*/, 0/*KmKm*/, 0/*KpKm*/, 1,/*PIpPIp*/, 1/*PImPIm*/, 1/*PIpPIm*/,\
    0 /*"LamLam"*/, 0 /*"aLam_aLam"*/, 0 /*"Lam_aLam"*/, 0/*all*/, 0/*plus*/, 0/*minus*/, 0/*mixed*/};
  
  int weightSetting = 3;
  if(strlen(params)!=0)
    {
      weightSetting = atoi(strtok(params, ","));
      runktdep = atoi(strtok(NULL, ",")); 
      resonancecut =  atoi(strtok(NULL, ","));
      for(int ll=0; ll<numOfPairTypes; ll++)
        {
          if(ll==0) cout<<"\nSELECTED PAIRS: ";
          runch[ll]= atoi(strtok(NULL, ","));
          cout<<runch[ll];
          if(ll==numOfPairTypes-1) cout<<endl;
        }
    }
    else
    {
      runch[12]=1;
    }
  
  
  ///////////////////////////////////////////
  
	// TGrid::Connect("alien://");
	//  TFile *filterFile = TFile::Open(filterPath);


  ///////////////////////////////////////////
  // TH2D *filter1 = new TH2D();
  // TH2D *filter2 = new TH2D();

  
   

/************** K_T SETTINGS *************/   
   const int numOfkTbins = 11;
   double kt_bins[numOfkTbins+1] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 100.0};   


   //  ***Reader Kinematics (Monte Carlo data)***
   AliFemtoEventReaderKinematicsChain* Reader=new AliFemtoEventReaderKinematicsChain();

/************** MANAGER *************/   

   //ten maanger zbiera analizy typowe dla alifemto, tylko
   AliFemtoManager* Manager = new AliFemtoManager();
   Manager->SetEventReader(Reader);

   
	   AliFemtoModelWeightGeneratorLednicky *tWeight = new AliFemtoModelWeightGeneratorLednicky();
     if(weightSetting==0)
        {
          tWeight->SetDefaultCalcPar(); // Default is CoulOn, QuantumOn, StrongOn, 3BodyOff, Square, T0ApproxOff
        }
      else if(weightSetting==1)
      {
          tWeight->SetQuantumOn();
          tWeight->SetCoulOff();
          tWeight-> SetStrongOff();
          tWeight->Set3BodyOff();
          tWeight->SetSquare();
          tWeight->SetT0ApproxOff();
      }
      else if(weightSetting==2)
      {
          tWeight->SetQuantumOff();
          tWeight->SetCoulOff();
          tWeight-> SetStrongOn();
          tWeight->Set3BodyOff();
          tWeight->SetSquare();
          tWeight->SetT0ApproxOff();
      }
      else if(weightSetting==3)
      {
          tWeight->SetQuantumOff();
          tWeight->SetCoulOn();
          tWeight-> SetStrongOff();
          tWeight->Set3BodyOff();
          tWeight->SetSquare();
          tWeight->SetT0ApproxOff();
      }
      else if(weightSetting==4)
      {
          tWeight->SetQuantumOff();
          tWeight->SetCoulOff();
          tWeight-> SetStrongOff();
          tWeight->Set3BodyOn();
          tWeight->SetSquare();
          tWeight->SetT0ApproxOff();
      }
       else if(weightSetting==5)
      {
          tWeight->SetQuantumOn();
          tWeight->SetCoulOn();
          tWeight-> SetStrongOff();
          tWeight->Set3BodyOff();
          tWeight->SetSquare();
          tWeight->SetT0ApproxOff();
      }
       else if(weightSetting==6)
      {
          tWeight->SetQuantumOn();
          tWeight->SetCoulOn();
          tWeight-> SetStrongOn();
          tWeight->Set3BodyOff();
          tWeight->SetSquare();
          tWeight->SetT0ApproxOff();
      }
	    else 
      {
          tWeight->SetQuantumOff();
          tWeight->SetCoulOff();
          tWeight-> SetStrongOff();
          tWeight->Set3BodyOff();
          // tWeight->SetSquare();
          tWeight->SetT0ApproxOff();
      }

	   AliFemtoModelManager *tModelManager = new AliFemtoModelManager(); 
	   AliFemtoModelGausRinvFreezeOutGenerator *tFreezepim = new AliFemtoModelGausRinvFreezeOutGenerator();
	   tFreezepim->SetSizeInv(2.0); 
		 tModelManager->AcceptFreezeOutGenerator(tFreezepim); 
		tModelManager->CreateCopyHiddenInfo(kFALSE); 





           		  

/************** CREATING MONITORS AND CUTS *************/   

// tzw monitory, histogramy lub wykresy ktore sa uzupelniane podczas analizy
//liczby w nawiasach to wielkosc tablicy wskaznikow
//to sa tak na prawde klasy, monitory generuja rozklady jedno-dwu-czastkowe

   AliFemtoVertexMultAnalysis      *analysis[640]; //glowna analiza, korelujemy tylko czastki ktore sa podobne w z-vertex(polozenie primary\
   vertex na osi z) i multiplicity - porownujemy tylko zderzenia ktory zaszly w podobnym miejscu i mialy podobne multiplicity
   AliFemtoBasicEventCut             *basic_event_cut[640]; //podstawowy cut na zderzenia (moga tez byc cuty na tracki i na pary czastek)
   AliFemtoCutMonitorEventMult    *cutPassEvM[640]; //standardowe do multiplicity, zwykle histogramy z liczba eventow
   AliFemtoCutMonitorEventMult    *cutFailEvM[640];
   AliFemtoCutMonitorEventVertex *cutPassEvV[640];
   AliFemtoCutMonitorEventVertex *cutFailEvV[640];
	//cut to ograniczenie na analizowane czastki, pojedyncze tracki
   AliFemtoMCTrackCut          *part1_mc_track_cut[640];
   AliFemtoMCTrackCut          *part2_mc_track_cut[640]; //nazwy niekoniecznie odpowiadaja faktycznie uzywanym detektorom
   AliFemtoMCTrackCut          *part3_mc_track_cut[640];
   //AliFemtoCutMonitorParticleYPtWithWeights *cutPass1YPtetaphitpc[640]; //monitory na pojedyncze czastki
   //AliFemtoCutMonitorParticleYPtWithWeights *cutFail1YPtetaphitpc[640];
   AliFemtoCutMonitorParticlePID *cutPass1PID[640]; //generuja wykresy PID
   AliFemtoCutMonitorParticlePID *cutFail1PID[640];
   //AliFemtoCutMonitorParticleYPtWithWeights *cutPass2YPtetaphitpc[640];
   //AliFemtoCutMonitorParticleYPtWithWeights *cutFail2YPtetaphitpc[640];
   AliFemtoCutMonitorParticlePID *cutPass2PID[640]; //generuja wykresy PID na pary
   AliFemtoCutMonitorParticlePID *cutFail2PID[640];

   AliFemtoCutMonitorParticleYPt *cutPass3YPt_normal[640]; //generuje wykresy rapidity od pedu
   AliFemtoCutMonitorParticleYPt *cutFail3YPt_normal[640];

   // AliFemtoCutMonitorParticleYPtWithWeights *cutPass3YPt[640]; //generuje wykresy rapidity od pedu
   // AliFemtoCutMonitorParticleYPtWithWeights *cutFail3YPt[640];
   // AliFemtoCutMonitorParticleYPtWithWeights *cutPass3YPtFilter[640]; //generuje wykresy rapidity od pedu
   // AliFemtoCutMonitorParticleYPtWithWeights *cutFail3YPtFilter[640];

   AliFemtoCutMonitorParticlePID *cutPass3PID[640]; //Pass - przeszly cuty, Fail - nie przeszly TRACKcutow (przeszly event cuty)
   AliFemtoCutMonitorParticlePID *cutFail3PID[640];

   //    AliFemtoShareQualityTPCEntranceSepPairCut         *antigamma_cutsame[640]; //ShareQualityPairCut - stnadardowy cut na pary czastek, \
   sprawdza jak duzo klastrow jest sharowanych przez dana pare trackow, jezeli jest ich duzo to dwie czastki sa najpewniej jedna \
   TPCEntranceSep - do sprawdzamy na wejsciu czy odleglosc miedzy trackami nie jest za mala
   //cut na pary czastek, uzywa share-quality, sprawdza czy nie mamy do czynienia z e+ i e- z kreacji par
   AliFemtoPairCutAntiGamma         *antigamma_cut[640];
   //AliFemtoPairCutRadialDistance         *antigamma_cut[640];  //tu liczy sie odlegosc katowa a nie zwykla
   //   AliFemtoChi2CorrFctn               *cchiqinvetaphitpc[640]; //?
   // AliFemtoPairCutPt                   *kt_sum_cuts[640]; //cut na sume pedow poprzecznych 
   // AliFemtoQinvCorrFctn               *cQinv_kt[640]; //Qinv najbardziej podst femtoskopowa funkcja
   
   AliFemtoQinvCorrFctn               *cQinv[640];
   AliFemtoCorrFctnDEtaDPhi         *cDetaDphi[640];
   // AliFemtoQinvCorrFctnWithWeights               *cQinvWB[640];
   // AliFemtoCorrFctnDEtaDPhiWithWeights         *cDetaDphiWB[640];
   AliFemtoModelCorrFctn            *cQinvModel[640];
   AliFemtoModelCorrFctn            *cQinv_kt[1280];
    AliFemtoModelCorrFctnDEtaDPhiRM    *cDetaDphiModel[640];
   // AliFemtoModelCorrFctnWithWeights              *cQinvModelWithWeights[640];
   // AliFemtoModelCorrFctnDEtaDPhiWithWeights         *cDetaDphiModelWithWeights[640];
    AliFemtoPairCutPt                   *kt_sum_cuts[1280]; //cut na sume pedow poprzecznych 
   AliFemtoKTPairCut              *kt_cut_bins[1280];
   AliFemtoKTPairCut              *kt_cut[1280];
      AliFemtoCutMonitorPairKT     *ktDist[640];
    AliFemtoPairCutMInv   *invmasscuts[640];





/************** MAIN PART OF ANALYSIS *************/   

   
   	//***Delta Eta Delta Phi analysis for identified systems (chg: 0-8)***
	// wielka petla po wszystkich czastkach
	//imainloop - zmienna iteracyjna
  	int imainloop = 0;

	//petla po krotnosciach
   for (int imult = 0; imult < numOfMultBins; imult++)
   {
      if (runmults[imult])
      {
		// po typach czastek
         for (int itype = 0; itype < numOfPairTypes; itype++)
         {
            if (runch[itype]) //sprawdzamy czy dany typ czastek jest brany pod uwage do analizy
            {
              
  	          
                 	if(itype==0 || itype==1)
                  tWeight->SetPairType(AliFemtoModelWeightGenerator::ProtonProton());
                else if(itype==2)
                  tWeight->SetPairType(AliFemtoModelWeightGenerator::ProtonAntiproton());
                else if(itype==3 || itype ==4)
                  tWeight->SetPairType(AliFemtoModelWeightGenerator::KaonPlusKaonPlus());
                else if(itype==5)
                  tWeight->SetPairType(AliFemtoModelWeightGenerator::KaonPlusKaonMinus());
                else if(itype==6 || itype ==7)
                  tWeight->SetPairType(AliFemtoModelWeightGenerator::PionPlusPionPlus());
                else if(itype==8)
                  tWeight->SetPairType(AliFemtoModelWeightGenerator::PionPlusPionMinus());
                else if(itype==9)
                  tWeight->SetPairType(AliFemtoModelWeightGenerator::LambdaLambda());
                else if(itype==10)
                  tWeight->SetPairType(AliFemtoModelWeightGenerator::AntilambdaAntilambda());
                else if(itype==11)
                  tWeight->SetPairType(AliFemtoModelWeightGenerator::LambdaAntilambda());
                else
                   tWeight->SetPairType(AliFemtoModelWeightGenerator::PairTypeNone()); 

                   tModelManager->AcceptWeightGenerator(tWeight);
               
           	   
               imainloop = itype * numOfMultBins + imult; //struktura danych(histogramow w tablicy): PPM0, PPM1, PPM2 ... aPaPM0, aPaPM1...
               int multmix = 5; //miksujemy tylko eventy o podobnych krotnosciach
               // if(imult == 4) multmix = 30; //whole multiplicity
               analysis[imainloop] = new AliFemtoVertexMultAnalysis(10, -10.0, 10.0, multmix, multbins[imult], multbins[imult+1]); //tworzenie analizy
               analysis[imainloop]->SetNumEventsToMix(10); //zwiekszamy statystyke mianownika f korelacyjnej, okreslamy tu z ilu eventow (miksujemy)
               analysis[imainloop]->SetMinSizePartCollection(1); //przynajmniej jedna czastka musi przejsc nasze cuty
               analysis[imainloop]->SetVerboseMode(kFALSE); 
				//ograniczenie na event = jedna kolizje
               //*** Event cut ***
               basic_event_cut[imainloop] = new AliFemtoBasicEventCut(); //tworzymy cut na multiplicity
               basic_event_cut[imainloop]->SetEventMult(1,100000); //cut na multiplicity dolny i gorny
            	//odcinamy zderzenie ktore dzieja sie w odleglosci wiekszej niz 10 cm od srodka detektora na osi z
               basic_event_cut[imainloop]->SetVertZPos(-10,10);//cm

               //****** event monitors **********   //najpierw tworzymy cuta, zaraz potem monitory
               cutPassEvM[imainloop] = new AliFemtoCutMonitorEventMult(Form("cutPass_%s_M%i", pair_types[itype], imult), 2000, 20000.5);
               cutFailEvM[imainloop] = new AliFemtoCutMonitorEventMult(Form("cutFail_%s_M%i", pair_types[itype], imult), 2000, 20000.5);
               basic_event_cut[imainloop]->AddCutMonitor(cutPassEvM[imainloop], cutFailEvM[imainloop]);
               //do wykresow Z-vertex
               cutPassEvV[imainloop] = new AliFemtoCutMonitorEventVertex(Form("cutPass_%s_M%i", pair_types[itype], imult));
               cutFailEvV[imainloop] = new AliFemtoCutMonitorEventVertex(Form("cutFail_%s_M%i", pair_types[itype], imult));
               basic_event_cut[imainloop]->AddCutMonitor(cutPassEvV[imainloop], cutFailEvV[imainloop]);


				//_________Cuty na poj czastki_______________________
               // ***** single particle track cuts *********
               part1_mc_track_cut[imainloop] = new AliFemtoMCTrackCut();
               part2_mc_track_cut[imainloop] = new AliFemtoMCTrackCut();
               part3_mc_track_cut[imainloop] = new AliFemtoMCTrackCut();
           		//ograniczenie na ladunek
               part1_mc_track_cut[imainloop]->SetCharge(1.0); //zostawia tylko dodatnie
               part2_mc_track_cut[imainloop]->SetCharge(-1.0);
            	//ograniczamy pseudopospiesznosc, tym wieksza im bardziej czastka leci w kierunku wiazki, 0 jesli leci prostopadle
            	//wartosci ponizej sa podyktowane geometria detektora
               part1_mc_track_cut[imainloop]->SetEta(-0.8, 0.8); //zostawia tylko w takim zakresie ety
               part2_mc_track_cut[imainloop]->SetEta(-0.8, 0.8);
               part3_mc_track_cut[imainloop]->SetEta(-0.8, 0.8);

               	//_______Cuty dla konkretnych typow czastek_________
            	//uzupelniane inne ograniczenia w zaleznosci od typu czastki
               if (itype == 0 ||itype == 1 ||itype == 2)//protons 0-2
               {
                  part1_mc_track_cut[imainloop]->SetPt(0.5,4); //Pt - ped poprzeczny
                  part2_mc_track_cut[imainloop]->SetPt(0.5,4);
                             
                  part1_mc_track_cut[imainloop]->SetPDG(2212); //PDG - numer czastki z bookletu, tutaj okreslamy ze to sa wlasnie protony
                  part2_mc_track_cut[imainloop]->SetPDG(2212);
               }
               if (itype == 3 ||itype == 4 ||itype == 5)//kaons 3-5
               {
                  part1_mc_track_cut[imainloop]->SetPt(0.3, 4);
                  part2_mc_track_cut[imainloop]->SetPt(0.3, 4);

                  part1_mc_track_cut[imainloop]->SetPDG(321);
                  part2_mc_track_cut[imainloop]->SetPDG(321);

               }
               if (itype == 6 ||itype == 7 ||itype == 8)//pions 6-8
               {
                  part1_mc_track_cut[imainloop]->SetPt(0.2,4);
                  part2_mc_track_cut[imainloop]->SetPt(0.2,4);

                  part1_mc_track_cut[imainloop]->SetPDG(211);
                  part2_mc_track_cut[imainloop]->SetPDG(211);   
                }
                if(itype == 9)  //lambdas 9-11
                {
                  part1_mc_track_cut[imainloop]->SetPt(0.7,4);
                  part1_mc_track_cut[imainloop]->SetPDG(3122);
                  part1_mc_track_cut[imainloop]->SetCharge(0.0);
                 
                }
            	if(itype == 10)
            	{
                  part2_mc_track_cut[imainloop]->SetCharge(0.0);
                  part2_mc_track_cut[imainloop]->SetPt(0.7,4);
                  part2_mc_track_cut[imainloop]->SetPDG(-3122);

            	}
               if(itype == 11)
                {
                  part1_mc_track_cut[imainloop]->SetPt(0.7,4.0);
                  part2_mc_track_cut[imainloop]->SetPt(0.7,4.0);

                  part1_mc_track_cut[imainloop]->SetCharge(0.0);
                  part2_mc_track_cut[imainloop]->SetCharge(0.0);
                  part1_mc_track_cut[imainloop]->SetPDG(3122);
                  part2_mc_track_cut[imainloop]->SetPDG(-3122);
                }
                if (itype == 12)//all
                {
                  part3_mc_track_cut[imainloop]->SetPt(0.12,4.0);
                }
                if (itype == 13 ||itype == 14 ||itype == 15)//plus,minus,mixed
                {
                  part1_mc_track_cut[imainloop]->SetPt(0.12,4.0);
                  part2_mc_track_cut[imainloop]->SetPt(0.12,4.0);

                }
         

               //**************** track Monitors ***************
               //tutaj jest jak wlaczalibysmy monitory na pojedyncze czastki
               //mozna odkomentowac, powinno dzialac
               if(itype<12)
                 {
                  double mass = 0.0;
                  if(itype<3)
                     mass = ProtonMass;
                  else if(itype>=3 & itype<6)
                     mass = KaonMass;
                  else if(itype>=6 & itype<9)
                     mass = PionMass;
                  else if(itype>=9 & itype<12)
                     mass = LambdaMass;

                	              
		                cutPass3YPt_normal[imainloop] = new AliFemtoCutMonitorParticleYPt(Form("cutPass_%s_M%i", pair_types[itype], imult), mass);
	                  cutFail3YPt_normal[imainloop] = new AliFemtoCutMonitorParticleYPt(Form("cutFail_%s_M%i", pair_types[itype], imult), mass);
	                  if(itype%3 ==0 ) part1_mc_track_cut[imainloop]->AddCutMonitor(cutPass3YPt_normal[imainloop], cutFail3YPt_normal[imainloop]);
	                  else if(itype%3 == 1) part2_mc_track_cut[imainloop]->AddCutMonitor(cutPass3YPt_normal[imainloop], cutFail3YPt_normal[imainloop]);
	              
                 }

				//_____________Cuty na pary________________________
               //******** Two - track cuts ************
        	
               
               //cut below is onlyc reated to store monitor
                 kt_cut[imainloop] = new AliFemtoKTPairCut(-1.0, 100.0); //loose cut
                 ktDist[imainloop] = new AliFemtoCutMonitorPairKT(Form("ktDist_%s", pair_types[itype]), 0.0, 10.0);           

                 if(resonancecut && itype == 8)
                 {
                 	// cout<<"pion resonancecut~"<<endl;
                 	invmasscuts[imainloop] = new AliFemtoPairCutMInv(PionMass, PionMass, 0.95*RhoMass, 1.05*RhoMass);
	                invmasscuts[imainloop] -> AddCutMonitor(ktDist[imainloop]);
	                analysis[imainloop] -> SetPairCut(invmasscuts[imainloop]);
             	 }
             	 else if(resonancecut && itype == 5)
                 {
                 	// cout<<"kaon resonancecut~"<<endl;
					invmasscuts[imainloop] = new AliFemtoPairCutMInv(KaonMass, KaonMass, 0.95*PhiMass, 1.05*PhiMass);
	                invmasscuts[imainloop] -> AddCutMonitor(ktDist[imainloop]);
	                analysis[imainloop] -> SetPairCut(invmasscuts[imainloop]);
             	 }
             	 else  analysis[imainloop] -> SetPairCut(kt_cut[imainloop]);
              
               if(itype >= 12)
               	{
               		//wyciecie par elektron pozyton
	              	antigamma_cut[imainloop] = new AliFemtoPairCutAntiGamma();
	                antigamma_cut[imainloop]->SetDataType(AliFemtoPairCut::kKine); //kinematics 
               		analysis[imainloop]->SetPairCut(antigamma_cut[imainloop]);
               	}

               analysis[imainloop]->SetEnablePairMonitors(kTRUE);
               
               //_____________Cuty na eventy______________________
               //***** Setting cuts ***********
               // setting event cut
               analysis[imainloop]->SetEventCut(basic_event_cut[imainloop]);


               //setting single track cuts
               if(itype==0 || itype==3 || itype==6) //positive like-sign
               {
                  analysis[imainloop]->SetFirstParticleCut(part1_mc_track_cut[imainloop]);
                  analysis[imainloop]->SetSecondParticleCut(part1_mc_track_cut[imainloop]);
               }
               if(itype==1 || itype==4 || itype==7)//negative like-sign
               {
                  analysis[imainloop]->SetFirstParticleCut(part2_mc_track_cut[imainloop]);
                  analysis[imainloop]->SetSecondParticleCut(part2_mc_track_cut[imainloop]);
               }
               if(itype==2 || itype==5 || itype==8)//unlike-sign
               {
                  analysis[imainloop]->SetFirstParticleCut(part1_mc_track_cut[imainloop]);
                  analysis[imainloop]->SetSecondParticleCut(part2_mc_track_cut[imainloop]);
               }
  	           if(itype == 9) 
  	            {
  	            	analysis[imainloop]->SetFirstParticleCut(part1_mc_track_cut[imainloop]);
                    analysis[imainloop]->SetSecondParticleCut(part1_mc_track_cut[imainloop]);
  	            }
  	            if(itype == 10) 
                 {
                    analysis[imainloop]->SetFirstParticleCut(part2_mc_track_cut[imainloop]);
                    analysis[imainloop]->SetSecondParticleCut(part2_mc_track_cut[imainloop]);
                 }
                 if(itype == 11) 
  	            {
  	            	analysis[imainloop]->SetFirstParticleCut(part1_mc_track_cut[imainloop]);
                    analysis[imainloop]->SetSecondParticleCut(part2_mc_track_cut[imainloop]);
  	            }    
      				  if(itype==12) //all
      	           	{	
      					  analysis[imainloop]->SetFirstParticleCut(part3_mc_track_cut[imainloop]);
      					  analysis[imainloop]->SetSecondParticleCut(part3_mc_track_cut[imainloop]);
      	           	}
      	           if(itype==13) //positive like-sign
         				{
         					analysis[imainloop]->SetFirstParticleCut(part1_mc_track_cut[imainloop]);
         					analysis[imainloop]->SetSecondParticleCut(part1_mc_track_cut[imainloop]);
         				}
      				  if(itype==14)//negative like-sign
         				{
         					analysis[imainloop]->SetFirstParticleCut(part2_mc_track_cut[imainloop]);
         					analysis[imainloop]->SetSecondParticleCut(part2_mc_track_cut[imainloop]);
         				}
         				if(itype==15)//unlike-sign
         				{
         					analysis[imainloop]->SetFirstParticleCut(part1_mc_track_cut[imainloop]);
         					analysis[imainloop]->SetSecondParticleCut(part2_mc_track_cut[imainloop]);
         				}
                 
               cQinvModel[imainloop] = new AliFemtoModelCorrFctn(Form("cQinv_Model_%s_M%i", pair_types[itype], imult), 20, 0, 2); //sprawdzic, czy wart arg sa ok
               if(resonancecut && (itype==8 || itype==5)) cQinvModel[imainloop]->SetPairSelectionCut(invmasscuts[imainloop]);
               cQinvModel[imainloop] -> ConnectToManager(tModelManager); 
               analysis[imainloop]->AddCorrFctn(cQinvModel[imainloop]);
               // cQinvModel[imainloop]->Report();

               // cDetaDphiModel[imainloop] = new AliFemtoModelCorrFctnDEtaDPhiRM(Form("cdedp_Model_%s_M%i", pair_types[itype], imult), 35, 35, PionMass, PionMass);
               // if(resonancecut && (itype==8 || itype==5)) cDetaDphiModel[imainloop]->SetPairSelectionCut(invmasscuts[imainloop]);
               // cDetaDphiModel[imainloop] -> ConnectToManager(tModelManager); 
               // analysis[imainloop]->AddCorrFctn(cDetaDphiModel[imainloop]);
 
                if (runktdep)
               {
                  int ktm;
                  for (int ikt=0; ikt<numOfkTbins; ikt++)
                  {
                     ktm = imainloop * numOfkTbins + ikt;
                      kt_cut_bins[ktm] = new AliFemtoKTPairCut(kt_bins[ikt], kt_bins[ikt+1]);
                      ktDist[ktm] = new AliFemtoCutMonitorPairKT(Form("ktDist_%s", pair_types[itype]), kt_bins[ikt], kt_bins[ikt+1]);
                      kt_cut_bins[ktm] -> AddCutMonitor(ktDist[ktm]);
                     
                     cQinv_kt[ktm] = new AliFemtoModelCorrFctn(Form("cQinv_Model_%s_M%i_kT%i", pair_types[itype], imult, ikt), 20, 0, 2);
                     cQinv_kt[ktm]->SetPairSelectionCut(kt_cut_bins[ktm]);
                     cQinv_kt[ktm]->ConnectToManager(tModelManager); 
                     analysis[imainloop]->AddCorrFctn(cQinv_kt[ktm]);
                  }
               }

               Manager->AddAnalysis(analysis[imainloop]);   
            }
         }
      }
   }
   return Manager;
}

