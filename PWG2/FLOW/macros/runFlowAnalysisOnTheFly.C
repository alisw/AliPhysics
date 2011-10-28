// Settings for the simulation of events 'on the fly': 
//  a) Determine how many events you want to create;
//  b) Set random or same seed for random generator;
//  c) Determine multiplicites of events;
//  d) Parametrize the phi distribution;
//   d1) Enable/disable uniform event-wise fluctuations of v2; 
//   d2) Enable/diable pt dependence of v2;  
//  e) Parametrize the pt distribution;
//  f) Determine how many times each sampled particle will be taken (simulating nonflow);
//  g) Configure detector's acceptance;
//  h) Decide which flow analysis methods you will use;
//  i) Define simple cuts for Reference Particle (RP) selection;
//  j) Define simple cuts for Particle of Interest (POI) selection;
//  k) Define the ranges for two subevents separated with eta gap (needed only for SP method);
//  l) Enable/disable usage of particle weights.

// a) Determine how many events you want to create:
Int_t iNevts = 1000; // total statistics

// b) Set random or same seed for random generator:
Bool_t bSameSeed = kFALSE; // if kTRUE, the created events are the same when re-doing flow analysis 'on the fly'   

// c) Determine multiplicites of events:
//    Remark 1: Multiplicity M for each event is sampled uniformly from interval iMinMult <= M < iMaxMult;
//    Remark 2: For constant M of e.g. 500 for each event, set iMinMult = 500 and iMaxMult = 501.
Int_t iMinMult = 500; // uniformly sampled multiplicity is >= iMinMult
Int_t iMaxMult = 501; // uniformly sampled multiplicity is < iMaxMult 

// d) Parametrize the phi distribution:
//    Remark 1: Hardwired is Fourier-like distribution f(phi) = (1/2pi)(1+sum_{n=1}^{4} 2v_n cos[n(phi-rp)]),             
//              where reaction plane (rp) is sampled uniformly for each event from interval [0,2pi]
Double_t dV1 = 0.0; // constant harmonic v1
Double_t dV2 = 0.05; // constant harmonic v2 
Double_t dV3 = 0.0; // constant harmonic v3
Double_t dV4 = 0.0; // constant harmonic v4
//    Remark 2: By default all harmonics are constant for each event and for each particle. However, for v2 
//              the uniform event-wise fluctuations or pt dependence can be enabled:
//  d1) Enable/disable uniform event-wise fluctuations of v2: 
Bool_t bUniformFluctuationsV2 = kFALSE; // enable uniform event-wise flow fluctuations (set than also dMinV2 and dMaxV2 bellow)
Double_t dMinV2 = 0.04; // lower boundary on v2, when bUniformFluctuationsV2 = kTRUE
Double_t dMaxV2 = 0.06; // upper boundary on v2, when bUniformFluctuationsV2 = kTRUE
//  d2) Enable/disable pt dependence of v2: 
Bool_t bPtDependentV2 = kFALSE; // enable pt dependence of v2 (set then also dV2vsPtMax and dV2vsPtCutOff bellow) 
Double_t dV2vsPtCutOff = 2.0; // up to pt = dV2vsPtCutOff v2 is growing linearly as a function of pt
Double_t dV2vsPtMax = 0.20; // for pt >= dV2vsPtCutOff, v2(pt) = dV2vsPtMax 

// e) Parametrize the pt distribution:
//    Remark: Hardwired is Boltzmann distribution f(pt) = pt*exp[-sqrt(dMass^2+pt^2)/dT] 
Double_t dMass = 0.13957; // mass in GeV/c^2 (e.g. m_{pions} = 0.13957)
Double_t dTemperature = 0.44; // "temperature" in GeV/c (increase this parameter to get more high pt particles) 

// f) Determine how many times each sampled particle will be taken in the analysis (simulating nonflow):
Int_t nTimes = 1; // e.g. for nTimes = 2, strong 2-particle nonflow correlations are introduced 

// g) Configure detector's acceptance:
Bool_t uniformAcceptance = kTRUE; // if kTRUE: detectors has uniform azimuthal acceptance.
                                  // if kFALSE: you will simulate detector with non-uniform acceptance in one or 
                                  // two sectors. For each sector you specify phiMin, phiMax and probability p. 
                                  // Then all particles emitted in direction phiMin < phi < phiMax will be taken 
                                  // with probability p. If p = 0, that sector is completely blocked. Set bellow 
                                  // phiMin1, phiMax1, p1 for the first sector and phiMin2, phiMax2, p2 for the second 
                                  // sector. If you set phiMin2 = phiMax2 = p2 = 0, only first non-uniform sector is 
                                  // simulated.
// 1st non-uniform sector:
Double_t phiMin1 = 60; // first non-uniform sector starts at this azimuth (in degrees)
Double_t phiMax1 = 120; // first non-uniform sector ends at this azimuth (in degrees)
Double_t p1 = 0.5; // probablitity that particles emitted in [phiMin1,phiMax1] are taken
// 2nd non-uniform sector:
Double_t phiMin2 = 0.; // first non-uniform sector starts at this azimuth (in degrees)
Double_t phiMax2 = 0.; // first non-uniform sector ends at this azimuth (in degrees)
Double_t p2 = 0.; // probablitity that particles emitted in [phiMin2,phiMax2] are taken

// h) Decide which flow analysis methods you will use:
Bool_t MCEP     = kTRUE; // Monte Carlo Event Plane
Bool_t SP       = kTRUE; // Scalar Product (a.k.a 'flow analysis with eta gaps')
Bool_t GFC      = kTRUE; // Generating Function Cumulants
Bool_t QC       = kTRUE; // Q-cumulants
Bool_t FQD      = kTRUE; // Fitted q-distribution
Bool_t LYZ1SUM  = kTRUE; // Lee-Yang Zero (sum generating function), first pass over the data
Bool_t LYZ1PROD = kTRUE; // Lee-Yang Zero (product generating function), first pass over the data
Bool_t LYZ2SUM  = kFALSE; // Lee-Yang Zero (sum generating function), second pass over the data
Bool_t LYZ2PROD = kFALSE; // Lee-Yang Zero (product generating function), second pass over the data
Bool_t LYZEP    = kFALSE; // Lee-Yang Zero Event Plane
Bool_t MH       = kFALSE; // Mixed Harmonics (used for strong parity violation studies) 
Bool_t NL       = kFALSE; // Nested Loops (neeed for debugging, only for developers)

// i) Define simple cuts for Reference Particle (RP) selection:
Double_t ptMinRP = 0.0; // in GeV
Double_t ptMaxRP = 10.0; // in GeV
Double_t etaMinRP = -1.;
Double_t etaMaxRP = 1.;
Double_t phiMinRP = 0.0; // in degrees
Double_t phiMaxRP = 360.0; // in degrees
Bool_t bUseChargeRP = kFALSE; // if kFALSE, RPs with both sign of charges are taken
Int_t chargeRP = 1; // +1 or -1

// j) Define simple cuts for Particle of Interest (POI) selection:
Double_t ptMinPOI = 0.0; // in GeV
Double_t ptMaxPOI = 10.0; // in GeV
Double_t etaMinPOI = -1.; // 
Double_t etaMaxPOI = 1.;
Double_t phiMinPOI = 0.0; // in degrees
Double_t phiMaxPOI = 360.0; // in degrees
Bool_t bUseChargePOI = kFALSE; // if kFALSE, POIs with both sign of charges are taken
Int_t chargePOI = -1; // +1 or -1

// k) Define the ranges for two subevents separated with eta gap (needed only for SP method):
Double_t etaMinA = -0.8; // minimum eta of subevent A
Double_t etaMaxA = -0.5; // maximum eta of subevent A
Double_t etaMinB = 0.5; // minimum eta of subevent B
Double_t etaMaxB = 0.8; // maximum eta of subevent B 

// l) Enable/disable usage of particle weights:
Bool_t usePhiWeights = kFALSE; // phi weights
Bool_t usePtWeights  = kFALSE; // pt weights 
Bool_t useEtaWeights = kFALSE; // eta weights

enum anaModes {mLocal,mLocalSource,mLocalPAR};
// mLocal: Analyze data on your computer using aliroot
// mLocalPAR: Analyze data on your computer using root + PAR files
// mLocalSource: Analyze data on your computer using root + source files
                                          
#include "TStopwatch.h"
#include "TObjArray"
#include "Riostream.h"
#include "TFile.h"

int runFlowAnalysisOnTheFly(Int_t mode=mLocal)
{
 // Beging analysis 'on the fly'.
 
 // a) Formal necessities....;
 // b) Initialize the flow event maker 'on the fly';
 // c) If enabled, access particle weights from external file; 
 // d) Configure the flow analysis methods;
 // e) Simple cuts for RPs;
 // f) Simple cuts for POIs;
 // g) Create and analyse events 'on the fly'; 
 // h) Create the output file and directory structure for the final results of all methods; 
 // i) Calculate and store the final results of all methods.
 
 // a) Formal necessities....:
 CheckUserSettings();
 WelcomeMessage();
 TStopwatch timer;
 timer.Start(); 
 LoadLibraries(mode);
 
 // b) Initialize the flow event maker 'on the fly':
 UInt_t uiSeed = 0; // if uiSeed is 0, the seed is determined uniquely in space and time via TUUID
 if(bSameSeed){uiSeed = 44;} 
 AliFlowEventSimpleMakerOnTheFly* eventMakerOnTheFly = new AliFlowEventSimpleMakerOnTheFly(uiSeed);
 eventMakerOnTheFly->SetMinMult(iMinMult);
 eventMakerOnTheFly->SetMaxMult(iMaxMult); 
 eventMakerOnTheFly->SetMass(dMass);
 eventMakerOnTheFly->SetTemperature(dTemperature);
 eventMakerOnTheFly->SetV1(dV1);
 eventMakerOnTheFly->SetV2(dV2);
 eventMakerOnTheFly->SetV3(dV3);
 eventMakerOnTheFly->SetV4(dV4); 
 if(bUniformFluctuationsV2)
 {
  eventMakerOnTheFly->SetUniformFluctuationsV2(bUniformFluctuationsV2); 
  eventMakerOnTheFly->SetMinV2(dMinV2);
  eventMakerOnTheFly->SetMaxV2(dMaxV2);
 }
 if(bPtDependentV2) 
 {
  eventMakerOnTheFly->SetPtDependentV2(bPtDependentV2);
  eventMakerOnTheFly->SetV2vsPtCutOff(dV2vsPtCutOff);
  eventMakerOnTheFly->SetV2vsPtMax(dV2vsPtMax);
 } 
 eventMakerOnTheFly->SetSubeventEtaRange(etaMinA,etaMaxA,etaMinB,etaMaxB); 
 eventMakerOnTheFly->SetNTimes(nTimes); 
 if(!uniformAcceptance)
 {
  eventMakerOnTheFly->SetUniformAcceptance(kFALSE);
  eventMakerOnTheFly->SetFirstSectorPhiMin(phiMin1);
  eventMakerOnTheFly->SetFirstSectorPhiMax(phiMax1);
  eventMakerOnTheFly->SetFirstSectorProbability(p1);
  eventMakerOnTheFly->SetSecondSectorPhiMin(phiMin2);
  eventMakerOnTheFly->SetSecondSectorPhiMax(phiMax2);
  eventMakerOnTheFly->SetSecondSectorProbability(p2);
 } 
 eventMakerOnTheFly->Init();

 // c) If enabled, access particle weights from external file: 
 TFile *fileWithWeights = NULL;
 TList *listWithWeights = NULL; 
 if(usePhiWeights||usePtWeights||useEtaWeights) 
 {
  fileWithWeights = TFile::Open("weights.root","READ");
  if(fileWithWeights) 
  {
   listWithWeights = (TList*)fileWithWeights->Get("weights");
  }
  else
  {
   cout << " WARNING: the file <weights.root> with weights from the previous run was not found."<<endl;
   break;
  }    
 } // end of if(usePhiWeights||usePtWeights||useEtaWeights) 

 // d) Configure the flow analysis methods:
 AliFlowAnalysisWithQCumulants *qc = NULL;
 AliFlowAnalysisWithCumulants *gfc = NULL;
 AliFlowAnalysisWithFittingQDistribution *fqd = NULL;
 AliFlowAnalysisWithLeeYangZeros *lyz1sum  = NULL;
 AliFlowAnalysisWithLeeYangZeros *lyz1prod = NULL;
 AliFlowAnalysisWithLeeYangZeros *lyz2sum  = NULL;
 AliFlowAnalysisWithLeeYangZeros *lyz2prod = NULL;
 AliFlowAnalysisWithLYZEventPlane *lyzep = NULL;
 AliFlowAnalysisWithScalarProduct *sp = NULL;
 AliFlowAnalysisWithMixedHarmonics *mh = NULL;
 AliFlowAnalysisWithNestedLoops *nl = NULL;
 AliFlowAnalysisWithMCEventPlane *mcep = NULL;   
 // MCEP = monte carlo event plane
 if(MCEP) 
 {
  //AliFlowAnalysisWithMCEventPlane *mcep = new AliFlowAnalysisWithMCEventPlane();
  mcep = new AliFlowAnalysisWithMCEventPlane();
  mcep->SetHarmonic(2); // default is v2
  mcep->Init();
 } // end of if(MCEP)
 // SP = Scalar Product 
 if(SP) 
 {
  sp = new AliFlowAnalysisWithScalarProduct();
  if(listWithWeights){sp->SetWeightsList(listWithWeights);}
  sp->SetUsePhiWeights(usePhiWeights);  
  sp->SetHarmonic(2);
  sp->SetApplyCorrectionForNUA(kFALSE);
  sp->Init();
 } // end of if(SP)
 // QC = Q-cumulants  
 if(QC) 
 { 
  qc = new AliFlowAnalysisWithQCumulants();
  if(listWithWeights){qc->SetWeightsList(listWithWeights);}
  if(usePhiWeights){qc->SetUsePhiWeights(usePhiWeights);}
  if(usePtWeights){qc->SetUsePtWeights(usePtWeights);}
  if(useEtaWeights){qc->SetUseEtaWeights(useEtaWeights);}
  qc->SetHarmonic(2);
  qc->SetCalculateDiffFlow(kTRUE);
  qc->SetCalculate2DDiffFlow(kFALSE); // vs (pt,eta)
  qc->SetApplyCorrectionForNUA(kFALSE);
  qc->SetFillMultipleControlHistograms(kFALSE);     
  qc->SetMultiplicityWeight("combinations"); // default (other supported options are "unit" and "multiplicity")
  qc->SetCalculateCumulantsVsM(kFALSE);
  qc->SetCalculateAllCorrelationsVsM(kFALSE); // calculate all correlations in mixed harmonics "vs M"
  qc->SetnBinsMult(10000);
  qc->SetMinMult(0);
  qc->SetMaxMult(10000);      
  qc->SetBookOnlyBasicCCH(kFALSE); // book only basic common control histograms
  qc->SetCalculateDiffFlowVsEta(kTRUE); // if you set kFALSE only differential flow vs pt is calculated
  qc->Init();  
 } // end of if(QC)
 // GFC = Generating Function Cumulants 
 if(GFC) 
 {
  gfc = new AliFlowAnalysisWithCumulants();
  if(listWithWeights){gfc->SetWeightsList(listWithWeights);}
  if(usePhiWeights){gfc->SetUsePhiWeights(usePhiWeights);}
  if(usePtWeights){gfc->SetUsePtWeights(usePtWeights);}
  if(useEtaWeights){gfc->SetUseEtaWeights(useEtaWeights);}
  // calculation vs multiplicity:
  gfc->SetCalculateVsMultiplicity(kFALSE);   
  gfc->SetnBinsMult(10000);
  gfc->SetMinMult(0);
  gfc->SetMaxMult(10000);   
  gfc->SetHarmonic(2);
  // tuning of interpolating parameters:
  gfc->SetTuneParameters(kFALSE);
  Double_t r0[10] = {1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7}; // up to 10 values allowed
  for(Int_t r=0;r<10;r++){gfc->SetTuningR0(r0[r],r);}
  gfc->Init();
 } // end of if(GFC) 
 // FQD = Fitting q-distribution 
 if(FQD) 
 {
  fqd = new AliFlowAnalysisWithFittingQDistribution();
  if(listWithWeights){fqd->SetWeightsList(listWithWeights);}
  if(usePhiWeights){fqd->SetUsePhiWeights(usePhiWeights);} 
  fqd->SetHarmonic(2); 
  fqd->Init();
 } // end of if(FQD)
 // LYZ1 = Lee-Yang Zeroes first run
 if(LYZ1SUM) 
 {
  lyz1sum = new AliFlowAnalysisWithLeeYangZeros();
  lyz1sum->SetFirstRun(kTRUE);
  lyz1sum->SetUseSum(kTRUE);
  lyz1sum->Init();
 } // end of if(LYZ1SUM)
 if(LYZ1PROD) 
 {
  lyz1prod = new AliFlowAnalysisWithLeeYangZeros();
  lyz1prod->SetFirstRun(kTRUE);
  lyz1prod->SetUseSum(kFALSE);
  lyz1prod->Init();
 } // end of if(LYZ1PROD)
 // LYZ2 = Lee-Yang Zeroes second run
 if(LYZ2SUM) 
 {
  lyz2sum = new AliFlowAnalysisWithLeeYangZeros();
  // read the input file from the first run 
  TString inputFileNameLYZ2SUM = "outputLYZ1SUManalysis.root" ;
  TFile* inputFileLYZ2SUM = new TFile(inputFileNameLYZ2SUM.Data(),"READ");
  if(!inputFileLYZ2SUM || inputFileLYZ2SUM->IsZombie()) 
  { 
   cerr<<" ERROR: To run LYZ2SUM you need the output file from LYZ1SUM. This file is not there! Please run LYZ1SUM first."<<endl;
   break; 
  } else 
    { 
     TList* inputListLYZ2SUM = (TList*)inputFileLYZ2SUM->Get("cobjLYZ1SUM");  
     if(!inputListLYZ2SUM){cout<<"Input list LYZ2SUM is NULL pointer!"<<endl;break;}
     else 
     {
      cout<<"LYZ2SUM input file/list read..."<<endl;
      lyz2sum->SetFirstRunList(inputListLYZ2SUM);
      lyz2sum->SetFirstRun(kFALSE);
      lyz2sum->SetUseSum(kTRUE);
      lyz2sum->Init();
     }
   }
 } // end of if(LYZ2SUM)
 if(LYZ2PROD) 
 {
  lyz2prod = new AliFlowAnalysisWithLeeYangZeros();
  // read the input file from the first run 
  TString inputFileNameLYZ2PROD = "outputLYZ1PRODanalysis.root" ;
  TFile* inputFileLYZ2PROD = new TFile(inputFileNameLYZ2PROD.Data(),"READ");
  if(!inputFileLYZ2PROD || inputFileLYZ2PROD->IsZombie()) 
  { 
   cerr<<" ERROR: To run LYZ2PROD you need the output file from LYZ1PROD. This file is not there! Please run LYZ1PROD first."<<endl;
   break; 
  } else 
    { 
     TList* inputListLYZ2PROD = (TList*)inputFileLYZ2PROD->Get("cobjLYZ1PROD");  
     if(!inputListLYZ2PROD){cout<<"Input list LYZ2PROD is NULL pointer!"<<endl;break;}
     else 
     {
      cout<<"LYZ2PROD input file/list read..."<<endl;
      lyz2prod->SetFirstRunList(inputListLYZ2PROD);
      lyz2prod->SetFirstRun(kFALSE);
      lyz2prod->SetUseSum(kFALSE);
      lyz2prod->Init();
     }
   }
 } // end of if(LYZ2PROD) 
 // LYZEP = Lee-Yang Zeroes event plane
 if(LYZEP) 
 {
  AliFlowLYZEventPlane *ep = new AliFlowLYZEventPlane() ;
  AliFlowAnalysisWithLYZEventPlane* lyzep = new AliFlowAnalysisWithLYZEventPlane();
   // read the input file from the second lyz run 
   TString inputFileNameLYZEP = "outputLYZ2SUManalysis.root" ;
   TFile* inputFileLYZEP = new TFile(inputFileNameLYZEP.Data(),"READ");
   if(!inputFileLYZEP || inputFileLYZEP->IsZombie()) { 
     cerr << " ERROR: To run LYZEP you need the output file from LYZ2SUM. This file is not there! Please run LYZ2SUM first." << endl ; 
     break;
   }
   else { 
     TList* inputListLYZEP = (TList*)inputFileLYZEP->Get("cobjLYZ2SUM");  
     if (!inputListLYZEP) {cout<<"Input list LYZEP is NULL pointer!"<<endl; break;}
     else {
       cout<<"LYZEP input file/list read..."<<endl;
       ep   ->SetSecondRunList(inputListLYZEP);
       lyzep->SetSecondRunList(inputListLYZEP);
       ep   ->Init();
       lyzep->Init();
     }
   }
 }
 // Mixed Harmonics:
 if(MH) 
 {
  mh = new AliFlowAnalysisWithMixedHarmonics();
  mh->SetHarmonic(1); // integer n in expression cos[n(2phi1-phi2-phi3)] = v2n*vn^2
  mh->SetMinMultiplicity(100); 
  mh->SetNoOfMultipicityBins(5);  
  mh->SetMultipicityBinWidth(200);   
  mh->Init(); 
 } // end of if(MH)
 // NL = Nested Loops:
 if(NL) 
 {
  nl = new AliFlowAnalysisWithNestedLoops();
  nl->Init();
 } // end of if(NL) 
 
 // e) Simple cuts for RPs: 
 AliFlowTrackSimpleCuts *cutsRP = new AliFlowTrackSimpleCuts();
 cutsRP->SetPtMax(ptMaxRP);
 cutsRP->SetPtMin(ptMinRP);
 cutsRP->SetEtaMax(etaMaxRP);
 cutsRP->SetEtaMin(etaMinRP);
 cutsRP->SetPhiMax(phiMaxRP*TMath::Pi()/180.);
 cutsRP->SetPhiMin(phiMinRP*TMath::Pi()/180.);
 if(bUseChargeRP){cutsRP->SetCharge(chargeRP);}

 // f) Simple cuts for POIs: 
 AliFlowTrackSimpleCuts *cutsPOI = new AliFlowTrackSimpleCuts();
 cutsPOI->SetPtMax(ptMaxPOI);
 cutsPOI->SetPtMin(ptMinPOI);
 cutsPOI->SetEtaMax(etaMaxPOI);
 cutsPOI->SetEtaMin(etaMinPOI);
 cutsPOI->SetPhiMax(phiMaxPOI*TMath::Pi()/180.);
 cutsPOI->SetPhiMin(phiMinPOI*TMath::Pi()/180.);
 if(bUseChargePOI){cutsPOI->SetCharge(chargePOI);}
                                       
 // g) Create and analyse events 'on the fly':
 for(Int_t i=0;i<iNevts;i++) 
 {   
  // Creating the event 'on the fly':
  AliFlowEventSimple *event = eventMakerOnTheFly->CreateEventOnTheFly(cutsRP,cutsPOI); 
  // Passing the created event to flow analysis methods:
  if(MCEP){mcep->Make(event);}
  if(QC){qc->Make(event);}
  if(GFC){gfc->Make(event);}
  if(FQD){fqd->Make(event);}
  if(LYZ1SUM){lyz1sum->Make(event);}
  if(LYZ1PROD){lyz1prod->Make(event);}
  if(LYZ2SUM){lyz2sum->Make(event);}
  if(LYZ2PROD){lyz2prod->Make(event);}
  if(LYZEP){lyzep->Make(event,ep);}
  if(SP){sp->Make(event);}
  if(MH){mh->Make(event);}
  if(NL){nl->Make(event);}
  delete event;
 } // end of for(Int_t i=0;i<iNevts;i++)

 // h) Create the output file and directory structure for the final results of all methods: 
 TString outputFileName = "AnalysisResults.root";  
 TFile *outputFile = new TFile(outputFileName.Data(),"RECREATE");
 const Int_t nMethods = 12;
 TString method[nMethods] = {"MCEP","SP","GFC","QC","FQD","LYZ1SUM","LYZ1PROD","LYZ2SUM","LYZ2PROD","LYZEP","MH","NL"};
 TDirectoryFile *dirFileFinal[nMethods] = {NULL};
 TString fileName[nMethods]; 
 for(Int_t i=0;i<nMethods;i++)
 {
  fileName[i]+="output";
  fileName[i]+=method[i].Data();
  fileName[i]+="analysis";
  dirFileFinal[i] = new TDirectoryFile(fileName[i].Data(),fileName[i].Data());
 } 
 
 // i) Calculate and store the final results of all methods:
 if(MCEP){mcep->Finish();mcep->WriteHistograms(dirFileFinal[0]);}
 if(SP){sp->Finish();sp->WriteHistograms(dirFileFinal[1]);}
 if(GFC){gfc->Finish();gfc->WriteHistograms(dirFileFinal[2]);}
 if(QC){qc->Finish();qc->WriteHistograms(dirFileFinal[3]);}
 if(FQD){fqd->Finish();fqd->WriteHistograms(dirFileFinal[4]);}
 if(LYZ1SUM){lyz1sum->Finish();lyz1sum->WriteHistograms(dirFileFinal[5]);}
 if(LYZ1PROD){lyz1prod->Finish();lyz1prod->WriteHistograms(dirFileFinal[6]);}
 if(LYZ2SUM){lyz2sum->Finish();lyz2sum->WriteHistograms(dirFileFinal[7]);}
 if(LYZ2PROD){lyz2prod->Finish();lyz2prod->WriteHistograms(dirFileFinal[8]);}
 if(LYZEP){lyzep->Finish();lyzep->WriteHistograms(dirFileFinal[9]);}
 if(MH){mh->Finish();mh->WriteHistograms(dirFileFinal[10]);}
 if(NL){nl->Finish();nl->WriteHistograms(dirFileFinal[11]);}
 
 outputFile->Close();
 delete outputFile;
 
 cout<<endl;
 cout<<endl;
 cout<<" ---- LANDED SUCCESSFULLY ---- "<<endl;
 cout<<endl; 
 
 timer.Stop();
 cout << endl;
 timer.Print();

} // end of int runFlowAnalysisOnTheFly(Int_t mode=mLocal)

void SetupPar(char* pararchivename)
{
  //Load par files, create analysis libraries
  //For testing, if par file already decompressed and modified
  //classes then do not decompress.
 
  TString cdir(Form("%s", gSystem->WorkingDirectory() )) ; 
  TString parpar(Form("%s.par", pararchivename)) ; 
  if ( gSystem->AccessPathName(parpar.Data()) ) {
    gSystem->ChangeDirectory(gSystem->Getenv("ALICE_ROOT")) ;
    TString processline(Form(".! make %s", parpar.Data())) ; 
    gROOT->ProcessLine(processline.Data()) ;
    gSystem->ChangeDirectory(cdir) ; 
    processline = Form(".! mv /tmp/%s .", parpar.Data()) ;
    gROOT->ProcessLine(processline.Data()) ;
  } 
  if ( gSystem->AccessPathName(pararchivename) ) {  
    TString processline = Form(".! tar xvzf %s",parpar.Data()) ;
    gROOT->ProcessLine(processline.Data());
  }
  
  TString ocwd = gSystem->WorkingDirectory();
  gSystem->ChangeDirectory(pararchivename);
  
  // check for BUILD.sh and execute
  if (!gSystem->AccessPathName("PROOF-INF/BUILD.sh")) {
    printf("*******************************\n");
    printf("*** Building PAR archive    ***\n");
    cout<<pararchivename<<endl;
    printf("*******************************\n");
    
    if (gSystem->Exec("PROOF-INF/BUILD.sh")) {
      Error("runProcess","Cannot Build the PAR Archive! - Abort!");
      return -1;
    }
  }
  // check for SETUP.C and execute
  if (!gSystem->AccessPathName("PROOF-INF/SETUP.C")) {
    printf("*******************************\n");
    printf("*** Setup PAR archive       ***\n");
    cout<<pararchivename<<endl;
    printf("*******************************\n");
    gROOT->Macro("PROOF-INF/SETUP.C");
  }
  
  gSystem->ChangeDirectory(ocwd.Data());
  printf("Current dir: %s\n", ocwd.Data());
}

void CheckUserSettings()
{
 // Check if user settings make sense before taking off.

 if(iNevts <= 0)
 {
  printf("\n WARNING: nEvts <= 0 !!!! Please check your settings before taking off.\n\n");
  exit(0);
 } 
 if(iMinMult < 0.)
 {
  printf("\n WARNING: iMinMult < 0 !!!! Please check your settings before taking off.\n\n");
  exit(0);
 }
 if(iMaxMult <= 0.)
 {
  printf("\n WARNING: iMaxMult <= 0 !!!! Please check your settings before taking off.\n\n");
  exit(0);
 }
 if(iMinMult >= iMaxMult)
 {
  printf("\n WARNING: iMinMult >= iMaxMult !!!! Please check your settings before taking off.\n\n");
  exit(0);
 }
 if(dMass < 0.)
 {
  printf("\n WARNING: dMass < 0 !!!! Please check your settings before taking off.\n\n");
  exit(0);
 }
 if(dTemperature <= 1e-44)
 {
  printf("\n WARNING: dTemperature <= 0 !!!! Please check your settings before taking off.\n\n");
  exit(0);
 }
 if(TMath::Abs(dV1) > 0.5)
 {
  printf("\n WARNING: |dV1| > 0.5 !!!! Please check your settings before taking off.\n\n");
  exit(0);
 }
 if(TMath::Abs(dV2) > 0.5)
 {
  printf("\n WARNING: |dV2| > 0.5 !!!! Please check your settings before taking off.\n\n");
  exit(0);
 }
 if(TMath::Abs(dV3) > 0.5)
 {
  printf("\n WARNING: |dV3| > 0.5 !!!! Please check your settings before taking off.\n\n");
  exit(0);
 }
 if(TMath::Abs(dV4) > 0.5)
 {
  printf("\n WARNING: |dV4| > 0.5 !!!! Please check your settings before taking off.\n\n");
  exit(0);
 }
 if(LYZ1SUM && LYZ2SUM)
 {
  cout<<" WARNING: You cannot run LYZ1 and LYZ2 at the same time! LYZ2 needs the output from LYZ1."<<endl; 
  exit(0); 
 }
 if(LYZ1PROD && LYZ2PROD)  
 {
  cout<<" WARNING: You cannot run LYZ1 and LYZ2 at the same time! LYZ2 needs the output from LYZ1."<<endl; 
  exit(0); 
 }
 if(LYZ2SUM && LYZEP) 
 {
  cout<<" WARNING: You cannot run LYZ2 and LYZEP at the same time! LYZEP needs the output from LYZ2."<<endl; 
  exit(0); 
 }
 if(LYZ1SUM && LYZEP) 
 {
  cout<<" WARNING: You cannot run LYZ1 and LYZEP at the same time! LYZEP needs the output from LYZ2."<<endl; 
  exit(0); 
 }

 if(!uniformAcceptance && phiMin1 > phiMax1)
 {
  cout<<" WARNING: You must have phiMin1 < phiMax1 !!!!"<<endl;
  exit(0);
 }
 if(!uniformAcceptance && !((TMath::Abs(phiMin2) < 1.e-44) && (TMath::Abs(phiMax2) < 1.e-44) && (TMath::Abs(p2) < 1.e-44)) 
    && (phiMin2 < phiMax1 || phiMin2 > phiMax2))
 {
  cout<<" WARNING: You must have phiMin2 > phiMax1 and phiMin2 < phiMax2 !!!!"<<endl;
  exit(0);
 }
 if((phiMin1 < 0 || phiMin1 > 360) || (phiMax1 < 0 || phiMax1 > 360) || 
    (phiMin2 < 0 || phiMin2 > 360) || (phiMax2 < 0 || phiMax2 > 360) )
 {
  cout<<" WARNING: You must take azimuthal angles from interval [0,360] !!!!"<<endl;
  exit(0);
 }
 if((p1 < 0 || p1 > 1) || (p2 < 0 || p2 > 1))
 {
  cout<<" WARNING: you must take p1 and p2 from interval [0,1] !!!!"<<endl;
  exit(0);
 }
 if(bPtDependentV2 && bUniformFluctuationsV2)
 {
  cout<<" WARNING: Uniform fluctuations not supported for pt denependent v2 !!!!"<<endl;
  exit(0);
 }

} // end of void CheckUserSettings()

void WelcomeMessage()
{
 // Welcome.

 cout<<endl;
 cout<<endl;
 cout<<"      ---- ARE YOU READY TO FLY ? ----      "<<endl;
 cout<<endl;
 
 gSystem->Sleep(1544);
 
 cout<<endl;
 cout<<" ---- BEGIN FLOW ANALYSIS 'ON THE FLY' ---- "<<endl;
 cout<<endl;
 cout<<endl;

 gSystem->Sleep(1544);

} // end of void WelcomeMessage()

void LoadLibraries(const anaModes mode) {
  
  //--------------------------------------
  // Load the needed libraries most of them already loaded by aliroot
  //--------------------------------------
  //gSystem->Load("libTree");
  gSystem->Load("libGeom");
  gSystem->Load("libVMC");
  gSystem->Load("libXMLIO");
  gSystem->Load("libPhysics");
  
  //----------------------------------------------------------
  // >>>>>>>>>>> Local mode <<<<<<<<<<<<<< 
  //----------------------------------------------------------
  if (mode==mLocal) {
    //--------------------------------------------------------
    // If you want to use already compiled libraries 
    // in the aliroot distribution
    //--------------------------------------------------------
    gSystem->Load("libSTEERBase");
    gSystem->Load("libESD");
    gSystem->Load("libAOD");
    gSystem->Load("libANALYSIS");
    gSystem->Load("libANALYSISalice");
    gSystem->Load("libCORRFW");
    cerr<<"libCORRFW loaded..."<<endl;
    gSystem->Load("libPWG2flowCommon");
    cerr<<"libPWG2flowCommon loaded..."<<endl;
    gSystem->Load("libPWG2flowTasks");
    cerr<<"libPWG2flowTasks loaded..."<<endl;
  }
  
  else if (mode == mLocalPAR) {
    //--------------------------------------------------------
    //If you want to use root and par files from aliroot
    //--------------------------------------------------------  
     //If you want to use root and par files from aliroot
    //--------------------------------------------------------  
    SetupPar("STEERBase");
    SetupPar("ESD");
    SetupPar("AOD");
    SetupPar("ANALYSIS");
    SetupPar("ANALYSISalice");
    SetupPar("PWG2AOD");
    SetupPar("CORRFW");
    SetupPar("PWG2flowCommon");
    cerr<<"PWG2flowCommon.par loaded..."<<endl;
    SetupPar("PWG2flowTasks");
    cerr<<"PWG2flowTasks.par loaded..."<<endl;
  }
  
  //---------------------------------------------------------
  // <<<<<<<<<< Source mode >>>>>>>>>>>>
  //---------------------------------------------------------
  else if (mode==mLocalSource) {
 
    // In root inline compile

    // Constants  
    gROOT->LoadMacro("AliFlowCommon/AliFlowCommonConstants.cxx+");
    gROOT->LoadMacro("AliFlowCommon/AliFlowLYZConstants.cxx+");
    
    // Flow event
    gROOT->LoadMacro("AliFlowCommon/AliFlowVector.cxx+"); 
    gROOT->LoadMacro("AliFlowCommon/AliFlowTrackSimple.cxx+");    
    gROOT->LoadMacro("AliFlowCommon/AliFlowTrackSimpleCuts.cxx+");    
    gROOT->LoadMacro("AliFlowCommon/AliFlowEventSimple.cxx+");
    
    // Output histosgrams
    gROOT->LoadMacro("AliFlowCommon/AliFlowCommonHist.cxx+");
    gROOT->LoadMacro("AliFlowCommon/AliFlowCommonHistResults.cxx+");
    gROOT->LoadMacro("AliFlowCommon/AliFlowLYZHist1.cxx+");
    gROOT->LoadMacro("AliFlowCommon/AliFlowLYZHist2.cxx+");
    
    // Functions needed for various methods
    gROOT->LoadMacro("AliFlowCommon/AliCumulantsFunctions.cxx+");
    gROOT->LoadMacro("AliFlowCommon/AliFlowLYZEventPlane.cxx+");
    
    // Flow Analysis code for various methods
    gROOT->LoadMacro("AliFlowCommon/AliFlowAnalysisWithMCEventPlane.cxx+"); 
    gROOT->LoadMacro("AliFlowCommon/AliFlowAnalysisWithScalarProduct.cxx+");
    gROOT->LoadMacro("AliFlowCommon/AliFlowAnalysisWithLYZEventPlane.cxx+");
    gROOT->LoadMacro("AliFlowCommon/AliFlowAnalysisWithLeeYangZeros.cxx+");
    gROOT->LoadMacro("AliFlowCommon/AliFlowAnalysisWithCumulants.cxx+");
    gROOT->LoadMacro("AliFlowCommon/AliFlowAnalysisWithQCumulants.cxx+"); 
    gROOT->LoadMacro("AliFlowCommon/AliFlowAnalysisWithFittingQDistribution.cxx+");
    gROOT->LoadMacro("AliFlowCommon/AliFlowAnalysisWithMixedHarmonics.cxx+"); 
    gROOT->LoadMacro("AliFlowCommon/AliFlowAnalysisWithNestedLoops.cxx+");
    
    // Class to fill the FlowEvent on the fly (generate Monte Carlo events)
    gROOT->LoadMacro("AliFlowCommon/AliFlowEventSimpleMakerOnTheFly.cxx+");   
    
    cout << "finished loading macros!" << endl;  
    
  }  
  
}


