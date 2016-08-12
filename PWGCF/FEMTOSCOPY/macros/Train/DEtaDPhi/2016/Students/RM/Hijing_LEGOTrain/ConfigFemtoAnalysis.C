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
//masy czastek w GeV
   double PionMass = 0.13956995;
   double KaonMass = 0.493677;
   double ProtonMass = 0.938272013;
   double LambdaMass = 1.115683;

   //double *part_mases[4] = {ProtonMass, KaonMass, PionMass, LambdaMass};

/************** MULTIPLICITY *************/   
	//ilosc binow krotnosci czastek, krotnosc (MULTIPLICITY) - ilosc czastek wyprodukowanych w zderzeniu; bierzemy biny np. 0-10, 11-20 itd...; multiplicity uzywa sie wymiennie z centralnoscia (albo jedno albo drugie),
	//zazwyczaj tworzy osobne histogramy dla czastek z danego binu multiplicity, piszemy wtedy M0, M1 ...
	//suma wszystkich wejsc na histogramie multiplicity daje liczbe wszystkich zderzen wzietych pod uwage
	//calka daje sume czastek we wszystkich zderzeniach
   const int numOfMultBins = 5;
   int runmults[numOfMultBins] = {1, 0, 0, 0, 0};    // 1 wlaczony przedzial krotnosci, 0 wylaczony
   int multbins[numOfMultBins+1] = {2, 20000, 0, 5000, 10001, 15000}; //definiujemy przedzialy krotnosci 

/************** PARTICLE TYPES FOR ANALYSIS *************/   
	//ilosc par dla ktorych liczymy korelacje
   const int numOfPairTypes = 16;
	//typy analizowanych par czastek, all - bierzemy (dowolna, dowolna), plus - (dodatnia, dodatnia), mixed (plus, minus) Ostateczna ilosc wykresow = liczba par * liczba binow multiplicity * [liczba binow pedow + 1] (+1 bo chcemy zazwyczaj wykres bez binowania)
   const char *pair_types[numOfPairTypes] = { "PP" /*0*/, "aPaP", "PaP", "KpKp", "KmKm", "KpKm", "PIpPIp", "PImPIm", "PIpPIm", \
   "LamLam" /*9*/, "aLam_aLam", "Lam_aLam", "all" /*12*/, "plus", "minus", "mixed"};
	// 1 wlaczaja nam pary do analizy, 0 wylaczaja
   int runch[numOfPairTypes] = {1/*PP*/, 1/*aPaP*/, 1/*PaP*/, 1/*KpKp*/, 1/*KmKm*/, 1/*KpKm*/, 1,/*PIpPIp*/, 1/*PImPIm*/, 1/*PIpPIm*/,\
    1 /*"LamLam"*/, 1 /*"aLam_aLam"*/, 1 /*"Lam_aLam"*/, 1/*all*/, 1/*plus*/, 1/*minus*/, 1/*mixed*/};

/************** K_T SETTINGS *************/   
	//ilosc binow przy pedzie; kt = (p1 - p2)/2 W TEJ ANALIZIE TO MOZE BYC pt = |pt1| + |pt2|
	//jak wybierzemy pare, to dzielimy pedy na cztery biny i pozniej mozemy robic histogram dla danej pary o zadanym zakresie pedow
   const int numOfkTbins = 4;
	//przedzialy kt, przedzialy kt - granice dla binow, liczba binow jest zdefiniiowana wyzej
   double kt_bins[numOfkTbins+1] = {0.0, 0.75, 1.5, 2.25, 100.0};   

/************** OTHER SETTINGS *************/   
	//jesli 1 to analizowane sa funkcje korelacyjne w roznych przedzialach sredniego pedu poprzecznego kt
   int runktdep = 0; //run kT dependance - jak 0 to nie ma analizy

   //wlaczenie konkretnego cutu na dwie czastki, ma wyrzucac elektrony i pozytony pochodzace od fotonow
   int gammacut = 1;   // cut for e+e- coming from gamma

/************** READER CHOICE *************/   

//wybieramy readera pliku poprzez comment/uncomment
//dosyc istotne, ustalamy tu z jaka analiza mamy do czynienia, z jakimi plikami pracujemy

   //  ***Reader ESD***
   //AliFemtoEventReaderESDChain *Reader = new AliFemtoEventReaderESDChain();
   //sposob liczenia multiplicity, ustawione aby liczyc track po tracku dla zdefiniowanych dobrze trackow, zeby dowiedziec sie jakie dokladnie czastki \
   brane pod uwage nalezy wejsc w szczegoly READERa \
   dla zdezen Pb-Pb wybiera sie po centrality - centralnosc zostala przyporzadkowana eventowi w trakcie rekonstrukcji
   //Reader->SetUseMultiplicity(AliFemtoEventReaderESDChain::kReferenceITSTPC);

   //  ***Reader AOD***
   //AliFemtoEventReaderAODChain *Reader = new AliFemtoEventReaderAODChain();
   //w AOD sa rozne rodzaje/zestawy trackow - filtry, ktore usuwaja czesci trackow, zostawiajac tylko interesujace (np. pochodzace tylko z niektorych detektorow)
   //Reader->SetFilterMask(96); //maska - pozwala wybrac zestaw filtrow, np wykluczajacych sie; 96 = 64+32 dwa filterbity
   //DCA - distance of closest approach - najblizsza mozliwa odleglosc, zazwyczaj miedzy trackiem a primary vertex(punkt zdezenia), dla analizy VIZIRo liczymy DCA \
   zeby zobaczyc czy 2 czastki wychodza z tego samego punktu
   //FB7 - filter bit 7
   //ponizszy seter mowi, zeby ped brac obliczenia dot DCA z global trackow (po wszystkich detektorach, ale faktycznie ITS i TPC)
   //Reader->SetDCAglobalTrack(kTRUE); //false for FB7, true for the rest //we do not use DCA at all - niektore filterbity maja juz wpisane DCA, 96 juz ma
  
   //sposob liczenia multiplicity   
   //Reader->SetUseMultiplicity(AliFemtoEventReaderAOD::kReference); //jak w ESD, tylko dla AOD, czastka po czastce
   //Pilap - nakladanie sie sygnalow z kilku eventow (trigger sobie nie zawsze radzi)
   //to sa cuty zeby poradzic sobie z pilap'em
   //Reader->SetMinPlpContribSPD(3); 
   //Reader->SetIsPileUpEvent(kTRUE);  //jezeli event jest oznaczony flaga pilap-event to odrzuc taki event

   //  ***Reader Kinematics (Monte Carlo data)***
   AliFemtoEventReaderKinematicsChain* Reader=new AliFemtoEventReaderKinematicsChain();

/************** MANAGER *************/   

   //ten maanger zbiera analizy typowe dla alifemto, tylko
   AliFemtoManager* Manager = new AliFemtoManager();
   Manager->SetEventReader(Reader);

/************** CREATING MONITORS AND CUTS *************/   

// tzw monitory, histogramy lub wykresy ktore sa uzupelniane podczas analizy
//liczby w nawiasach to wielkosc tablicy wskaznikow
//to sa tak na prawde klasy, monitory generuja rozklady jedno-dwu-czastkowe

   AliFemtoVertexMultAnalysis      *analysis[640]; //glowna analiza, korelujemy tylko czastki ktore sa podobne w z-vertex(polozenie primary\
   vertex na osi z) i multiplicity - porownujemy tylko zderzenia ktory zaszly w podobnym miejscu i mialy podobne multiplicity
   AliFemtoBasicEventCut             *basic_event_cut[320]; //podstawowy cut na zderzenia (moga tez byc cuty na tracki i na pary czastek)
   AliFemtoCutMonitorEventMult    *cutPassEvM[320]; //standardowe do multiplicity, zwykle histogramy z liczba eventow
   AliFemtoCutMonitorEventMult    *cutFailEvM[320];
   AliFemtoCutMonitorEventVertex *cutPassEvV[320];
   AliFemtoCutMonitorEventVertex *cutFailEvV[320];
	//cut to ograniczenie na analizowane czastki, pojedyncze tracki
   AliFemtoMCTrackCut          *part1_mc_track_cut[320];
   AliFemtoMCTrackCut          *part2_mc_track_cut[320]; //nazwy niekoniecznie odpowiadaja faktycznie uzywanym detektorom
   AliFemtoMCTrackCut          *part3_mc_track_cut[320];
   //AliFemtoCutMonitorParticleYPt *cutPass1YPtetaphitpc[320]; //monitory na pojedyncze czastki
   //AliFemtoCutMonitorParticleYPt *cutFail1YPtetaphitpc[320];
   AliFemtoCutMonitorParticlePID *cutPass1PID[320]; //generuja wykresy PID
   AliFemtoCutMonitorParticlePID *cutFail1PID[320];
   //AliFemtoCutMonitorParticleYPt *cutPass2YPtetaphitpc[320];
   //AliFemtoCutMonitorParticleYPt *cutFail2YPtetaphitpc[320];
   AliFemtoCutMonitorParticlePID *cutPass2PID[320]; //generuja wykresy PID na pary
   AliFemtoCutMonitorParticlePID *cutFail2PID[320];

   AliFemtoCutMonitorParticleYPt *cutPass3YPt[320]; //generuje wykresy rapidity od pedu
   AliFemtoCutMonitorParticleYPt *cutFail3YPt[320];

   AliFemtoCutMonitorParticlePID *cutPass3PID[320]; //Pass - przeszly cuty, Fail - nie przeszly TRACKcutow (przeszly event cuty)
   AliFemtoCutMonitorParticlePID *cutFail3PID[320];
   //    AliFemtoShareQualityTPCEntranceSepPairCut         *antigamma_cutsame[320]; //ShareQualityPairCut - stnadardowy cut na pary czastek, \
   sprawdza jak duzo klastrow jest sharowanych przez dana pare trackow, jezeli jest ich duzo to dwie czastki sa najpewniej jedna \
   TPCEntranceSep - do sprawdzamy na wejsciu czy odleglosc miedzy trackami nie jest za mala
   //cut na pary czastek, uzywa share-quality, sprawdza czy nie mamy do czynienia z e+ i e- z kreacji par
   AliFemtoPairCutAntiGamma         *antigamma_cut[320];
   //AliFemtoPairCutRadialDistance         *antigamma_cut[320];  //tu liczy sie odlegosc katowa a nie zwykla
   //   AliFemtoChi2CorrFctn               *cchiqinvetaphitpc[320]; //?
   AliFemtoPairCutPt                   *kt_sum_cuts[640]; //cut na sume pedow poprzecznych 
   AliFemtoQinvCorrFctn               *cQinv_kt[320]; //Qinv najbardziej podst femtoskopowa funkcja
   AliFemtoQinvCorrFctn               *cQinv[320];
	//komponenty funkcji korelacyjnej
   AliFemtoCorrFctnDEtaDPhi         *cDetaDphi[640];


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

               imainloop = itype * numOfMultBins + imult; //struktura danych(histogramow w tablicy): PPM0, PPM1, PPM2 ... aPaPM0, aPaPM1...
               int multmix = 5; //miksujemy tylko eventy o podobnych krotnosciach
               if(imult == 4) multmix = 30; //whole multiplicity
               analysis[imainloop] = new AliFemtoVertexMultAnalysis(10, -10.0, 10.0, multmix, multbins[imult], multbins[imult+1]); //tworzenie analizy
               analysis[imainloop]->SetNumEventsToMix(10); //zwiekszamy statystyke mianownika f korelacyjnej, okreslamy tu z ilu eventow (miksujemy)
               analysis[imainloop]->SetMinSizePartCollection(1); //przynajmniej jedna czastka musi przejsc nasze cuty
				//ograniczenie na event = jedna kolizje
               //*** Event cut ***
               basic_event_cut[imainloop] = new AliFemtoBasicEventCut(); //tworzymy cut na multiplicity
               basic_event_cut[imainloop]->SetEventMult(0.001,100000); //cut na multiplicity dolny i gorny
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
               part1_mc_track_cut[imainloop]->SetEta(-0.8,0.8); //zostawia tylko w takim zakresie ety
               part2_mc_track_cut[imainloop]->SetEta(-0.8,0.8);
               part3_mc_track_cut[imainloop]->SetEta(-0.8,0.8);

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
                  part1_mc_track_cut[imainloop]->SetPt(0.3,4);
                  part2_mc_track_cut[imainloop]->SetPt(0.3,4);

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
                  part1_mc_track_cut[imainloop]->SetPt(0.7,4);
                  part2_mc_track_cut[imainloop]->SetPt(0.7,4);
                  part1_mc_track_cut[imainloop]->SetCharge(0.0);
                  part2_mc_track_cut[imainloop]->SetCharge(0.0);

                  part1_mc_track_cut[imainloop]->SetPDG(3122);
                  part2_mc_track_cut[imainloop]->SetPDG(-3122);
                }
                if (itype == 12)//all
                {
                  part3_mc_track_cut[imainloop]->SetPt(0.12,4);
                }
                if (itype == 13 ||itype == 14 ||itype == 15)//plus,minus,mixed
                {
                  part1_mc_track_cut[imainloop]->SetPt(0.12,4);
                  part2_mc_track_cut[imainloop]->(0.12,4);
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
                  cutPass3YPt[imainloop] = new AliFemtoCutMonitorParticleYPt(Form("cutPass_%s_M%i", pair_types[itype], imult), mass);
                  cutFail3YPt[imainloop] = new AliFemtoCutMonitorParticleYPt(Form("cutFail_%s_M%i", pair_types[itype], imult), mass);
                  if(itype%3 ==0 ) part1_mc_track_cut[imainloop]->AddCutMonitor(cutPass3YPt[imainloop], cutFail3YPt[imainloop]);
                  else if(itype%3 == 1) part2_mc_track_cut[imainloop]->AddCutMonitor(cutPass3YPt[imainloop], cutFail3YPt[imainloop]);
                 }

				//_____________Cuty na pary________________________
               //******** Two - track cuts ************
        		//wyciecie par elektron pozyton
               antigamma_cut[imainloop] = new AliFemtoPairCutAntiGamma();
               //antigamma_cut[imainloop] = new AliFemtoPairCutRadialDistance();
               antigamma_cut[imainloop]->SetDataType(AliFemtoPairCut::kKine); //kinematics 
               
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
           
               //setting two-track cuts
               analysis[imainloop]->SetPairCut(antigamma_cut[imainloop]);

               //______________Funkcje korelacyjne___________________
               //**** Correlation functions *******
               //dodajemy funkcje korelacyjna
               cQinv[imainloop] = new AliFemtoQinvCorrFctn(Form("cQinv_%s_M%i", pair_types[itype], imult), 20, 0, 2); //sprawdzic, czy wart arg sa ok
               analysis[imainloop]->AddCorrFctn(cQinv[imainloop]);

               cDetaDphi[imainloop] = new AliFemtoCorrFctnDEtaDPhi(Form("cdedp_%s_M%i", pair_types[itype], imult),35, 35);
               analysis[imainloop]->AddCorrFctn(cDetaDphi[imainloop]);

               //jesli 1 to analizowane sa funkcje korelacyjne w roznych przedzialach sredniego pedu poprzecznego kt
               if (runktdep)
               {
                  //liczenie dla pedu poprzecznego, jesli to wlaczylismy
                  //generuje dodatkowe histogramy
                  int ktm;
                  for (int ikt=0; ikt<numOfkTbins; ikt++)
                  {
                     ktm = imainloop * numOfkTbins + ikt;
                     kt_sum_cuts[ktm] = new AliFemtoPairCutPt(kt_bins[ikt], kt_bins[ikt+1]);

                     cQinv[ktm] = new AliFemtoQinvCorrFctn(Form("cQinv_%s_M%ipT%i", pair_types[itype], imult, ikt), 20, 0, 2);
                     cQinv[ktm]->SetPairSelectionCut(kt_sum_cuts[ktm]);
                     analysis[imainloop]->AddCorrFctn(cQinv[ktm]);

                     cDetaDphi[ktm] = new AliFemtoCorrFctnDEtaDPhi(Form("cdedp_%s_M%ipT%i", pair_types[itype], imult,ikt),35, 35);
                     cDetaDphi[ktm]->SetPairSelectionCut(kt_sum_cuts[ktm]);
                     analysis[imainloop]->AddCorrFctn(cDetaDphi[ktm]);

                  }
               }

            	//dodajemy do Managera(glownego programu) funkcje z analiza
               Manager->AddAnalysis(analysis[imainloop]);   
            }
         }
      }
   }
   return Manager;
}

