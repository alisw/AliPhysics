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
		#include "AliFemtoModelCorrFctnQinv.h"
		#include "AliFemtoCorrFctnDEtaDPhi.h"
		#include "AliFemtoModelCorrFctn.h"
		#include "AliFemtoModelWeightGenerator.h"
		#include "AliFemtoModelCorrFctnDEtaDPhi.h"
		#include "AliFemtoAvgSepCorrFctn.h"
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
		AliFemtoManager* ConfigFemtoAnalysis() {

		   double PionMass = 0.13956995;
		   double KaonMass = 0.493677;
		   double ProtonMass = 0.938272013;
		   double LambdaMass = 1.115683;

		   //lambda0 (3122), antylambda0 (-3122) -> 3 nowe korelacje: lambda0-lambda0, antylambda-antylambda, lambda0-antylambda


		   // 0-10% centralno??, to znaczy 10% wszystkich zderze?, kt?re s? zdarzeniami o najwi?kszej multiplicity
		   // multiplicity- jeden ze sposobow jakim mozemy mierzyc centralnosc

		   // bierzemy bin i wypelniamy inne histogramy
		   // inaczej wygladaja femtoskopowe funkcje korelacyjne dla zderzen o niskiej i wysokiej krotnosci
		   // w jednym binie liczba czastek ze zderzenia - moze byc kilka zderzen o tej samej liczbie czastek
		   // suma po wejsciach histogramu daje liczbe czastek we wszystkich zderzeniach
		   
		   const int numOfMultBins = 5; // ilosc binow krotnosci czastek - multiplicity  = ilosc czastek wyprodukowanych w zderzeniu
		   // bardzo ciezko zdefiniowac multiplicity, bo nie wiemy, ile bylo w zderzeniu, tylko ile nam zlapal detektor
		   const int numOfPairBins = 16; // biny dla par dla ktorych liczymy korelacje (liczba par)
		   const int numOfpTbins = 4; // roznica pedow czastek
		   // w tej analizie to ptSum = suma pedow pT dwoch czastek, kt - wektorowo, a tutaj -moduly

		   // pedy dzielimy na 4 biny po wybraniu pary i mozemy robic histogram dla danej pary o zadanym zakresie pedow

		   // wylaczone sa teraz mixed, plus, minus
		   int pairs[numOfPairBins] = {1,1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}; // wlaczamy konkretne pary
		   // ostateczna ilosc wykresow = liczba par*liczba binowmultiplicity*(liczba binow pedow+1(bo chcemy wykres bez binowania w pedzie))
		   const char *pairTypes[numOfPairBins] = { "PP", "aPaP", "PaP", "KpKp", "KmKm", "KpKm", "PIpPIp", "PImPIm", "PIpPIm", "all", "plus", "minus", "mixed","LL","LaL","aLaL"};
		   // mixed - jedna dodatnia, druga - ujemna, plus - dodatnie hadrony, all - dowolna z dowolna,minus - obie ujemne
		   
		   int runpTdep = 0; // chcemy wlaczac analize w pT, czy nie - true/false
		   double ptrng[numOfpTbins+1] = {0.0, 0.75, 1.5, 2.25, 100.0}; // tablica = granice dla binow kT

		   // wlaczenie cutu na konkretne czastki
		   int gammacut = 1;   // cut for e+e- coming from gamma

		   //dosc istotna czesc
		   // ustalamy, z ktora analiza mamy do czynienia: esd, aod, kinematics

		   //  ***Reader ESD***
		   //AliFemtoEventReaderESDChain *Reader = new AliFemtoEventReaderESDChain();

		   // kreference - wazne!, ustawiamy, czy chcemy pod wzgledem krotnosci, czy centralnosci
		   // kluster - punkt, ktory ma zrekonstruowany ped
		   // liczy truck po trucku, ale tylko dla zdefiniowanych dobrze truckow
		   // ustawiamy konkretne cuty na czastke i liczymy po tej(konkretnej) liczbie czastek
		   // dla Pb - Pb- po centralnosci - w momencie rekonstrukcji zostaly uznane centralnosci dla poszczegolnych eventow
		   // centralnosc wylicza sie najpierw dla wszystkich eventow, zeby miec referencje i powiedziec,jaka centralnosc ma dany event

		   // istnieja trucki ktore sa tylko TPC
		   // sa takie co maja wejscia w ITD, jesli chcemy badac czastki pierwotne

		   //Reader->SetUseMultiplicity(AliFemtoEventReaderESDChain::kReferenceITSTPC);

		   //  ***Reader AOD***
		   //AliFemtoEventReaderAODChain *Reader = new AliFemtoEventReaderAODChain();
		   //Reader->SetFilterMask(96); // maska jest juz zbiorem truckow z odpowiednimi cut'ami
		   // maska ustawia kilka filter bitow na raz
		   // 96 rozpisane binarnie 64+32 - sa to 2 filterbity, ktore wypelniaja dziure w eta
		   // DCA = distance of closest approach (najblizsza mozliwa odleglosc pomiedzy czyms a czyms - zazwyczaj pomiedzy punktem a primary vertexem-punktem zderzenia w pierwotnym zderzeniu)
		   // wszystkie wytworzone czastki pochodza z primary vertexu
		   // np.liczy,y DCA dla 2 czastek po rekonstrukcji i sprawdzamy jak blisko siebie byly czastki

		   // FB7 - filterbit 7
		   //Reader->SetDCAglobalTrack(kTRUE); //false for FB7, true for the rest //we do not use DCA at all
		   // mowi, zeby obliczenia dotyczace DCA brac z global truckow - ITS i TPC
		   // global trucki - w teorii najbardziej podstawowe ktore zawieraja wszystkie info
		   // ale w filtrze 96 juz jest ustawiony filterbit, wiec jego juz nie uzywamy 

		   // tak samo znowu ustawiamy sposob liczenia czastka po czastce
		   //Reader->SetUseMultiplicity(AliFemtoEventReaderAOD::kReference);
		   //Reader->SetMinPlpContribSPD(3); // pochodzi PileUp- poprzednie eventy ktore wczesniej sie odbyly nachodza na nasza probke i wychodzi mieszanka
		   // wiekszosc trigger wycina, ale nie wszystko
		   // sa cuty, ktore dodatkowo wyrzucaja taki pileup
		   // eventy pileUp'owe sa oznaczone czasem i flaga
		   //Reader->SetIsPileUpEvent(kTRUE);

		   // wszystko sa to ustawienia reader'a

		   //  ***Reader Kinematics (Monte Carlo data)***
		   AliFemtoEventReaderKinematicsChain* Reader=new AliFemtoEventReaderKinematicsChain();
		   Reader->IsMisalignment(kTRUE);

		   // typowy manager dla AliFemto, ktory zbiera analizy typowe dla Alifemto, jest on tylko dla femtoskopii, a tamte byly ogolne
		   AliFemtoManager* Manager = new AliFemtoManager();
		   Manager->SetEventReader(Reader);

		      /************** WEIGHT GENERATOR ***********/
	   AliFemtoModelWeightGeneratorLednicky *tWeight = new AliFemtoModelWeightGeneratorLednicky();
	   AliFemtoModelWeightGeneratorBasic *tWeightBasic = new AliFemtoModelWeightGeneratorBasic();
	   // tWeight->SetPairType(AliFemtoModelWeightGenerator::ProtonProton());
	   tWeight ->SetCoulOn();
	   tWeight -> SetQuantumOn();
	   tWeight -> SetStrongOn();
	   tWeight -> Set3BodyOn();

	   AliFemtoModelManager *tModelManager = new AliFemtoModelManager(); 
	 //   AliFemtoModelGausLCMSFreezeOutGenerator *tFreezepim = new AliFemtoModelGausLCMSFreezeOutGenerator();
	 //  	tFreezepim->SetSizeOut(1.8*TMath::Sqrt(2.0));                                                
	 //  	tFreezepim->SetSizeSide(1.3*TMath::Sqrt(2.0));                                               
	 //  	tFreezepim->SetSizeLong(1.6*TMath::Sqrt(2.0));      
		// tModelManager->AcceptFreezeOutGenerator(tFreezepim); 
		tModelManager->CreateCopyHiddenInfo(kFALSE); 


		   // korelujemy czastki, ktore sa tylko podobne w multiplicity i z - vertex
		   // zwykle mamy cut, ktore bierze w przedziale odleg?o?ci od -10cm do 10cm
		   // z - vertex = polozenie primary vertexu na osi z
		   // VertexMult - wymuszamy na naszych zdarzeniach, zeby byly podobne w z-vertexie
		   // czyli analizujemy tylko zdarzenia ktore maja podobna multiplicity i sa w takim peczku wokol 0 w danym przedziale osi z
		   AliFemtoVertexMultAnalysis      *analysis[640];
		   AliFemtoBasicEventCut             *eventCuts[320];
		   // cut na eventy - na konkretne eventy
		   AliFemtoCutMonitorEventMult    *cutPassEvMult[320];
		   AliFemtoCutMonitorEventMult    *cutFailEvMult[320];
		   AliFemtoCutMonitorEventVertex *cutPassEvVertex[320];
		   AliFemtoCutMonitorEventVertex *cutFailEvVertex[320];
		   // cut na pojedyncza czastke
		   AliFemtoMCTrackCut          *dtruckCut1[320]; 
		   AliFemtoMCTrackCut          *dtruckCut2[320];
		   AliFemtoMCTrackCut          *dtruckCut3[320];
		   // AliFemtoCutMonitorParticleYPt *cutPass1YPtetaphitpc[320];
		   // AliFemtoCutMonitorParticleYPt *cutFail1YPtetaphitpc[320];

		   //Monitory - generuja rozklady glownie 1-czastkowe, ale tez 2-czastkowe - sa to histogramy, ktore beda podczas analizy
		   //  te klasy zawieraja klasy predefiniowanych histogramow
		   // trzeba obiekt zapisac do TH1D, a potem do drzewa rootowego
		   AliFemtoCutMonitorParticlePID *cutPass1PID[320];
		   AliFemtoCutMonitorParticlePID *cutFail1PID[320];
		   //AliFemtoCutMonitorParticleYPt *cutPass2YPtetaphitpc[320];
		   //AliFemtoCutMonitorParticleYPt *cutFail2YPtetaphitpc[320];

		   // generuja zestaw histogramow dla pid konkretnych, argument = liczba binow
		   // nazwa nie ma znaczenia, pass/fail tylko sa wazne - daja info o czastkach, ktore przeszly cuty albo nie przeszly cutow
		   // wszystkie cuty na trucki sa traktowane razem
		   // dla 1 czastki mamy tylko 1 truckcut'a
		   // najpierw event cut, potem truck cut, potem pair cut
		   AliFemtoCutMonitorParticlePID *cutPass2PID[320];
		   AliFemtoCutMonitorParticlePID *cutFail2PID[320];

		   // monitory na histogram rapidity Y od pT dla poszczegolnych czastek
		   AliFemtoCutMonitorParticleYPt *cutPassYPt[320];
		   AliFemtoCutMonitorParticleYPt *cutFailYPt[320];
		   AliFemtoCutMonitorParticlePID *cutPass3PID[320];
		   AliFemtoCutMonitorParticlePID *cutFail3PID[320];


		   // pair cuty  probuja sie sprzeciwic- splitting(2 wiazki dzieli - jedna jako 2 widzi), merging - 2 jako jedna widzi
		   // share quality - sprawdzamy jak duzo 2 takie same trucki maja klustrow jednakowych
		   // TPCEntranceSeparation - patrzymy, czy odleglosc miedzy 2 truckami jest za mala - bo wtedy moze byc wlasnie merging albo splitting
		   //    AliFemtoShareQualityTPCEntranceSepPairCut         *sqpcetaphitpcsame[320];

		   // cut na pary czastek
		   // zamiast sprawdzac odleglosc w TPC uzywa rowniez ShareQuality
		   AliFemtoPairCutAntiGamma         *pairCutAntiGamma[320];
		   //AliFemtoPairCutRadialDistance         *pairCutAntiGamma[320]; // "fi star" - liczymy odleglosc katowa, a nie zwykla
		   //   AliFemtoChi2CorrFctn               *cchiqinvetaphitpc[320]; // ?
		   AliFemtoPairCutPt                   *pTCuts[640]; // cut na sume pedow porpzecznych (dla 2 czastek) - bedziemy uzywac, jak bysmy chcieli wlaczyc binowanie w ptsum
		   AliFemtoQinvCorrFctn               *qinv_pTCorrFcn[320];



		   AliFemtoQinvCorrFctn               *qinvCorrFcn[320]; // moja tablica na Qinv


		   AliFemtoCorrFctnDEtaDPhi         *cdEtadPhi[640];

		   	//komponenty funkcji korelacyjnej
	   AliFemtoCorrFctnDEtaDPhi         *cDetaDphi[640];
	   AliFemtoModelCorrFctn              *cQinvModel[320];
	   AliFemtoModelCorrFctnDEtaDPhi         *cDetaDphiModel[640];


		   
		   //***Delta Eta Delta Phi analysis for identified systems (chg: 0-8)***
		   
		   int corrFcnIter = 0; // sluzy do iteracji liczby 0-320 po kolei wszystkie funkcje korelacyjne
		   int runmults[numOfMultBins] = {1, 0, 0, 0, 0}; // wlaczamy przedzialy krotnosci
		   int multbins[numOfMultBins+1] = {0, 20000, 2000, 1500, 1800, 2500}; // definiujemy przedzialy krotnosci ; dla Pb-Pb to sa tysiace
		    
		   for (int imult = 0; imult < numOfMultBins; imult++)
		   {
		      if (runmults[imult]) // jesli wlaczony typ czastek w multiplicity
		      {
		         for (int pairBinIter = 0; pairBinIter < numOfPairBins; pairBinIter++)
		         {
		            if (pairs[pairBinIter])  // jesli wlaczony po typach czastek
		            {
		            	if(pairBinIter==0)
	                  tWeight->SetPairType(AliFemtoModelWeightGenerator::ProtonProton());
	              	else if(pairBinIter==2)
	              		tWeight->SetPairType(AliFemtoModelWeightGenerator::ProtonAntiproton());
	              	else if(pairBinIter==3)
	              		tWeight->SetPairType(AliFemtoModelWeightGenerator::KaonPlusKaonPlus());
	              	else if(pairBinIter==5)
	              		tWeight->SetPairType(AliFemtoModelWeightGenerator::KaonPlusKaonMinus());
	              	else if(pairBinIter==6)
	              		tWeight->SetPairType(AliFemtoModelWeightGenerator::PionPlusPionPlus());
	              	else if(pairBinIter==8)
	              		tWeight->SetPairType(AliFemtoModelWeightGenerator::PionPlusPionMinus());
	              	else
	                   tWeight->SetPairType(AliFemtoModelWeightGenerator::PairTypeNone()); 

	               tModelManager->AcceptWeightGenerator(tWeight);

		               corrFcnIter = pairBinIter * numOfMultBins + imult; // wskazuje na to, ktory bin jest wlaczony -> wskazuje na obiekt przechowujacy funkcje korelacyjna

		               // patrzymy na wykres i ustawiamy odpowiednia liczbe multmix - liczba, na jaka dzielimy kazdy z binow w multiplicity, zeby brac
		               // czastki o podobnym multiplicity
		               int multmix = 100; // mixujemy tylko eventy o podobnych krotnosciach - zalezy od tego, jak szeroki jest bin multiplicity
		               // tej linijki na razie nie wykorzystuje; chodzi o to, zeby ustawic maksymalna wartosc z przedzialow multiplicity
		               if(imult == 4) multmix = 30; //whole multiplicity
		               analysis[corrFcnIter] = new AliFemtoVertexMultAnalysis(10, -10.0, 10.0, multmix, multbins[imult], multbins[imult+1]); // wlacza obliczenia w danym binie
		               // (UInt_t binsVe00rtex,Double_t minVertex,Double_t maxVertex,UInt_t binsMult,Double_t minMult, Double_t maxMult)

		               // oglada sie wykresy, ktore wychodza; jesli trojkat nie jest plaski, to zwiekszamy, w innym wypadku - nie warto
		               // dla statystyki wystarczy 10, bo juz mamy trojkat
		               analysis[corrFcnIter]->SetNumEventsToMix(5);// ustawiamy z iloma eventami, z ktorymi mieszamy ten konkretny event
		               // jak wezmiemy wiecej, to zwieksza nam statystyka, ale i tak mamy trojkat - za duzo nie trzeba, bo nie pomaga, tylko zajmuje czas
		               analysis[corrFcnIter]->SetMinSizePartCollection(1); // w danej analizie musi byc co najmniej 1 czastka, ktora przejdzie cuta, inaczej nic nie bedzie w detektorze zarejestrowane
		               analysis[corrFcnIter]->SetVerboseMode(kFALSE);
		               //*** Event cut ***
		               eventCuts[corrFcnIter] = new AliFemtoBasicEventCut(); // tworzymy po raz 1-szy event cut
		               eventCuts[corrFcnIter]->SetEventMult(0.001,100000); // tnie na multiplicity - a naprawde wyrzuca eventy, ktore maja 0 w multiplicity, args:min, max multiplicity
		               eventCuts[corrFcnIter]->SetVertZPos(-10,10);//cm // odleg?o?? w detektorze od 0 wzd?u? osi, gdzie zachodzi zderzenie

		               //****** event monitors **********   
		               // monitory do event cuta
		               // -szy zdefiniowany dla multiplicity - w katalogu pwgcf ->alifemto
		               // MC - mamy czastki nienaladowane rowniez(uwaga na roznice z eksperymentem)
		               cutPassEvMult[corrFcnIter] = new AliFemtoCutMonitorEventMult(Form("cutPass%sM%i", pairTypes[pairBinIter], imult), 2000, 20000.5);
		               cutFailEvMult[corrFcnIter] = new AliFemtoCutMonitorEventMult(Form("cutFail%sM%i", pairTypes[pairBinIter], imult), 2000, 20000.5);
		               eventCuts[corrFcnIter]->AddCutMonitor(cutPassEvMult[corrFcnIter], cutFailEvMult[corrFcnIter]);
		      
		               cutPassEvVertex[corrFcnIter] = new AliFemtoCutMonitorEventVertex(Form("cutPass%sM%i", pairTypes[pairBinIter], imult));
		               cutFailEvVertex[corrFcnIter] = new AliFemtoCutMonitorEventVertex(Form("cutFail%sM%i", pairTypes[pairBinIter], imult));
		               eventCuts[corrFcnIter]->AddCutMonitor(cutPassEvVertex[corrFcnIter], cutFailEvVertex[corrFcnIter]);

		              

		               // ***** single particle track cuts *********
		               dtruckCut1[corrFcnIter] = new AliFemtoMCTrackCut();
		               dtruckCut2[corrFcnIter] = new AliFemtoMCTrackCut();
		               dtruckCut3[corrFcnIter] = new AliFemtoMCTrackCut();
		          

		               dtruckCut1[corrFcnIter]->SetCharge(1.0); /// ustawiamy konkretny ladunek dla cutow na czastki
		               dtruckCut2[corrFcnIter]->SetCharge(-1.0);

		               // cuty na czastki, a nic nie zwiazane z funkcja korelacyjna
		               dtruckCut1[corrFcnIter]->SetEta(-0.8,0.8); // plateau obciete, dlatego 0.8 - przechodza takie ktore sie zalapuja o max przekatna detektora
		               dtruckCut2[corrFcnIter]->SetEta(-0.8,0.8);
		               dtruckCut3[corrFcnIter]->SetEta(-0.8,0.8);


		               // ustawiamy to dla cuty - ze bedziemy analizowali tylko protony o takim pedzie, takiej eta i takim ladunku
		               if (pairBinIter == 0 ||pairBinIter == 1 ||pairBinIter == 2)//protons 0-2
		               {
		               	  
		                  dtruckCut1[corrFcnIter]->SetPt(0.5,4); // ped poprzeczny 02. piony, 0.3 kaony, 0.5protony
		                  dtruckCut2[corrFcnIter]->SetPt(0.5,4); // prootony z materialow sa wybijane, maja nizsze pedy, datego je wycinamy
		                  // zawsze mozemy porownywac tylko czastki o tych samych pedach!!!!
		                                           
		                  dtruckCut1[corrFcnIter]->SetPDG(2212);
		                  dtruckCut2[corrFcnIter]->SetPDG(2212);
		               }
		               if (pairBinIter == 3 ||pairBinIter == 4 ||pairBinIter == 5)//kaons 3-5
		               {


		                  dtruckCut1[corrFcnIter]->SetPt(0.3,4); // 1 rodzaj czastek
		                  dtruckCut2[corrFcnIter]->SetPt(0.3,4); // 2 rodzaj czastek

		                  dtruckCut1[corrFcnIter]->SetPDG(321);
		                  dtruckCut2[corrFcnIter]->SetPDG(321);

		               }
		               if (pairBinIter == 6 ||pairBinIter == 7 ||pairBinIter == 8)//pions 6-8
		               {


		                  dtruckCut1[corrFcnIter]->SetPt(0.2,4);
		                  dtruckCut2[corrFcnIter]->SetPt(0.2,4);

		                  dtruckCut1[corrFcnIter]->SetPDG(211);
		                  dtruckCut2[corrFcnIter]->SetPDG(211);   
		                }
		                if (pairBinIter == 9)//all
		                {
		                  dtruckCut3[corrFcnIter]->SetPt(0.12,4); // TPC wyciaga min 0.12
		                }
		                if (pairBinIter == 10 ||pairBinIter == 11 ||pairBinIter == 12)//plus,minus,mixed
		                {
		                  dtruckCut1[corrFcnIter]->SetPt(0.12,4);
		                  dtruckCut2[corrFcnIter]->SetPt(0.12,4);
		                }

		                // zrobic 3 ify
		                if (pairBinIter == 13 ||pairBinIter == 14 ||pairBinIter == 15)//lambda0-lambda0, antylambda-antylambda, lambda0- antylambda0
		                {

		                	// dla lambd trzeba ustawic ladunek na zero i PDG w odroznieniu od innych czastek ma byc ujemny dla antylambd
		                  dtruckCut1[corrFcnIter]->SetCharge(0.0);
		                  dtruckCut2[corrFcnIter] -> SetCharge(0.0);

		                  dtruckCut1[corrFcnIter]->SetPDG(3122);
		                  dtruckCut2[corrFcnIter]->SetPDG(-3122);   

		                  dtruckCut1[corrFcnIter]->SetPt(0.7,4);
		                  dtruckCut2[corrFcnIter]->SetPt(0.7,4);
		                }



		                // monitory na pojedyncze trucki, to by wypelnialo dodatkowe histogramy
		                // prawdopodobnie mozna odkomentowac
		               //**************** track Monitors ***************
		       


		                if (pairBinIter==0 || pairBinIter==1){

		                  cutPassYPt[corrFcnIter] = new AliFemtoCutMonitorParticleYPt(Form("cutPass%sM%i", pairTypes[pairBinIter], imult),ProtonMass);
		                  cutFailYPt[corrFcnIter] = new AliFemtoCutMonitorParticleYPt(Form("cutFail%sM%i", pairTypes[pairBinIter], imult),ProtonMass);

		                }
		                else if (pairBinIter==3 || pairBinIter==4){

		               	  cutPassYPt[corrFcnIter] = new AliFemtoCutMonitorParticleYPt(Form("cutPass%sM%i", pairTypes[pairBinIter], imult),KaonMass);
		                  cutFailYPt[corrFcnIter] = new AliFemtoCutMonitorParticleYPt(Form("cutFail%sM%i", pairTypes[pairBinIter], imult),KaonMass);
		                }
		                else  if (pairBinIter==6 || pairBinIter==7){

		               	  cutPassYPt[corrFcnIter] = new AliFemtoCutMonitorParticleYPt(Form("cutPass%sM%i", pairTypes[pairBinIter], imult),PionMass);
		                  cutFailYPt[corrFcnIter] = new AliFemtoCutMonitorParticleYPt(Form("cutFail%sM%i", pairTypes[pairBinIter], imult),PionMass);
		                }
		                else if (pairBinIter==13 || pairBinIter==15){

		                  cutPassYPt[corrFcnIter] = new AliFemtoCutMonitorParticleYPt(Form("cutPass%sM%i", pairTypes[pairBinIter], imult),LambdaMass);
		                  cutFailYPt[corrFcnIter] = new AliFemtoCutMonitorParticleYPt(Form("cutFail%sM%i", pairTypes[pairBinIter], imult),LambdaMass);


		               }
		               // else{

		               // 	   cutPassYPt[corrFcnIter] = new AliFemtoCutMonitorParticleYPt(Form("cutPass%sM%i", pairTypes[pairBinIter], imult),PionMass);
		               //     cutFailYPt[corrFcnIter] = new AliFemtoCutMonitorParticleYPt(Form("cutFail%sM%i", pairTypes[pairBinIter], imult),PionMass);
		               // }


		                   if(pairBinIter==0||pairBinIter==3||pairBinIter==6 || pairBinIter==13) dtruckCut1[corrFcnIter]->AddCutMonitor(cutPassYPt[corrFcnIter], cutFailYPt[corrFcnIter]);
		                   if(pairBinIter==1||pairBinIter==4||pairBinIter==7 || pairBinIter==15) dtruckCut2[corrFcnIter]->AddCutMonitor(cutPassYPt[corrFcnIter], cutFailYPt[corrFcnIter]);
		               //  }


		               //******** Two - track cuts ************
		               pairCutAntiGamma[corrFcnIter] = new AliFemtoPairCutAntiGamma();
		               //pairCutAntiGamma[corrFcnIter] = new AliFemtoPairCutRadialDistance();
		               pairCutAntiGamma[corrFcnIter]->SetDataType(AliFemtoPairCut::kKine); //kinematics
		               

		               //dla MC nie ma na szczesci splitting, merging,itd
		               
		               //***** Setting cuts ***********
		               // setting event cut
		               analysis[corrFcnIter]->SetEventCut(eventCuts[corrFcnIter]);
		               //setting single track cuts
		                // dla kazdego typu czastek ustawiamy cuty
		               if(pairBinIter==0 || pairBinIter==3 || pairBinIter==6) //positive like-sign
		               {
		                 // ustawiamy czastki o tym samy znaku (jednakowo zdefiniowane)
		                  analysis[corrFcnIter]->SetFirstParticleCut(dtruckCut1[corrFcnIter]);
		                  analysis[corrFcnIter]->SetSecondParticleCut(dtruckCut1[corrFcnIter]);
		               }
		               if(pairBinIter==1 || pairBinIter==4 || pairBinIter==7)//negative like-sign
		               {
		                  analysis[corrFcnIter]->SetFirstParticleCut(dtruckCut2[corrFcnIter]);
		                  analysis[corrFcnIter]->SetSecondParticleCut(dtruckCut2[corrFcnIter]);
		               }
		               if(pairBinIter==2 || pairBinIter==5 || pairBinIter==8)//unlike-sign // rozne-znaki
		               {
		                  analysis[corrFcnIter]->SetFirstParticleCut(dtruckCut1[corrFcnIter]);
		                  analysis[corrFcnIter]->SetSecondParticleCut(dtruckCut2[corrFcnIter]);
		               }
							if(pairBinIter==9) //all
		               {
							  analysis[corrFcnIter]->SetFirstParticleCut(dtruckCut3[corrFcnIter]);
							  analysis[corrFcnIter]->SetSecondParticleCut(dtruckCut3[corrFcnIter]);
		               }
		               if(pairBinIter==10) //positive like-sign
							{
								analysis[corrFcnIter]->SetFirstParticleCut(dtruckCut1[corrFcnIter]);
								analysis[corrFcnIter]->SetSecondParticleCut(dtruckCut1[corrFcnIter]);
							}
							if(pairBinIter==11)//negative like-sign
							{
								analysis[corrFcnIter]->SetFirstParticleCut(dtruckCut2[corrFcnIter]);
								analysis[corrFcnIter]->SetSecondParticleCut(dtruckCut2[corrFcnIter]);
							}
							if(pairBinIter==12)//unlike-sign
							{
								analysis[corrFcnIter]->SetFirstParticleCut(dtruckCut1[corrFcnIter]);
								analysis[corrFcnIter]->SetSecondParticleCut(dtruckCut2[corrFcnIter]);
							}
		                  if(pairBinIter==13)//unlike-sign
		               {
		                  analysis[corrFcnIter]->SetFirstParticleCut(dtruckCut1[corrFcnIter]);
		                  analysis[corrFcnIter]->SetSecondParticleCut(dtruckCut1[corrFcnIter]);
		               }
		                  if(pairBinIter==14)//unlike-sign
		               {
		                  analysis[corrFcnIter]->SetFirstParticleCut(dtruckCut1[corrFcnIter]);
		                  analysis[corrFcnIter]->SetSecondParticleCut(dtruckCut2[corrFcnIter]);
		               }
		                  if(pairBinIter==15)//unlike-sign
		               {
		                  analysis[corrFcnIter]->SetFirstParticleCut(dtruckCut2[corrFcnIter]);
		                  analysis[corrFcnIter]->SetSecondParticleCut(dtruckCut2[corrFcnIter]);
		               }
		               
		               //setting two-track cuts
		               analysis[corrFcnIter]->SetPairCut(pairCutAntiGamma[corrFcnIter]);

		               //**** Correlation functions *******

		               // dodawane sa funkcje korelacyjne
		               // trzeba dodac zeby wyswietlal histogram q_invariant
		               // dodajemy funkcje korelacyjne, dopisac 2 linijki dla Qinv
		               cdEtadPhi[corrFcnIter] = new AliFemtoCorrFctnDEtaDPhi(Form("cdedp%sM%i", pairTypes[pairBinIter], imult),35, 35);
		               analysis[corrFcnIter]->AddCorrFctn(cdEtadPhi[corrFcnIter]);

		               qinvCorrFcn[corrFcnIter] = new AliFemtoQinvCorrFctn(Form("qinv%sM%i", pairTypes[pairBinIter], imult),40,0,2.0); // Qinv_Low, Qinv_High
		               analysis[corrFcnIter] -> AddCorrFctn(qinvCorrFcn[corrFcnIter]);

	      //______________Funkcje korelacyjne___________________
	               //**** Correlation functions *******
	               //dodajemy funkcje korelacyjna
	               // cQinv[corrFcnIter ] = new AliFemtoQinvCorrFctn(Form("NumcQinv_%s_M%i", pair_types[itype], imult), 20, 0, 2); //sprawdzic, czy wart arg sa ok
	               // analysis[corrFcnIter ]->AddCorrFctn(cQinv[corrFcnIter ]);

	              //  cDetaDphi[corrFcnIter ] = new AliFemtoCorrFctnDEtaDPhi(Form("cdedp_%s_M%i", pair_types[itype], imult),35, 35);
	              //  analysis[corrFcnIter ]->AddCorrFctn(cDetaDphi[corrFcnIter ]);
	              // //Model correlation functions 
	               cQinvModel[corrFcnIter ] = new AliFemtoModelCorrFctn(Form("cQinv_Model_%s_M%i", pairTypes[pairBinIter], imult), 20, 0, 2); //sprawdzic, czy wart arg sa ok
	               // cQinvModel[corrFcnIter ] -> SetPairType(AliFemtoAvgSepCorrFctn::kTracks);
	               // cQinvModel[corrFcnIter ] -> SetExpectedPDGCodes(2122,2122);
	               cQinvModel[corrFcnIter ] -> ConnectToManager(tModelManager); 
	               analysis[corrFcnIter ]->AddCorrFctn(cQinvModel[corrFcnIter ]);
	               cQinvModel[corrFcnIter ]->Report();

	             
	               cDetaDphiModel[corrFcnIter ] = new AliFemtoModelCorrFctnDEtaDPhi(Form("cdedp_Model_%s_M%i", pairTypes[pairBinIter], imult),35, 35);
	               cDetaDphiModel[corrFcnIter ] -> ConnectToManager(tModelManager); 
	               analysis[corrFcnIter ]->AddCorrFctn(cDetaDphiModel[corrFcnIter ]);


		               if (runpTdep)
		               {




	         

		                  int ktm;
		                  for (int ipt=0; ipt<numOfpTbins; ipt++)
		                  {
		                     ktm = corrFcnIter * numOfpTbins + ipt;
		                     pTCuts[ktm] = new AliFemtoPairCutPt(ptrng[ipt], ptrng[ipt+1]);

		                     cdEtadPhi[ktm] = new AliFemtoCorrFctnDEtaDPhi(Form("qinv%sM%ipT%i", pairTypes[pairBinIter], imult,ipt),35, 35);
		                     cdEtadPhi[ktm]->SetPairSelectionCut(pTCuts[ktm]);
		                     analysis[corrFcnIter]->AddCorrFctn(cdEtadPhi[ktm]);


		                     qinvCorrFcn[ktm] = new AliFemtoQinvCorrFctn(Form("qinv%sM%ipT%i", pairTypes[pairBinIter], imult,ipt),40,0,2.0);
		                     qinvCorrFcn[ktm]->SetPairSelectionCut(pTCuts[ktm]);
		                     analysis[corrFcnIter]->AddCorrFctn(qinvCorrFcn[ktm]);


		                  }
		               }      
		               Manager->AddAnalysis(analysis[corrFcnIter]);   
		            }
		         }
		      }
		   }
		   return Manager;
		}

