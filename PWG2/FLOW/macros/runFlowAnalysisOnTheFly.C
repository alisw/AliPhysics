#include "TStopwatch.h"
#include "TObjArray"
#include "Riostream.h"
#include "TFile.h"

//--------------------------------------------------------------------------------------
// RUN SETTINGS
// flow analysis method can be: (set to kTRUE or kFALSE)
Bool_t MCEP     = kTRUE;
Bool_t SP       = kTRUE;
Bool_t GFC      = kTRUE;
Bool_t QC       = kTRUE;
Bool_t FQD      = kTRUE;
Bool_t LYZ1SUM  = kTRUE;
Bool_t LYZ1PROD = kTRUE;
Bool_t LYZ2SUM  = kFALSE;
Bool_t LYZ2PROD = kFALSE;
Bool_t LYZEP    = kFALSE;
Bool_t MH       = kTRUE; // mixed harmonics 
Bool_t NL       = kFALSE; // nested loops
Bool_t MCEP_AH  = kFALSE; // MCEP in another harmonic 
//--------------------------------------------------------------------------------------

// Weights 
// use weights for Q-vector:
Bool_t usePhiWeights = kFALSE; // phi weights (correction for non-uniform azimuthal acceptance)
Bool_t usePtWeights  = kFALSE; // pt weights 
Bool_t useEtaWeights = kFALSE; // eta weights

// Define the range for eta subevents
Double_t minA = -0.9;
Double_t maxA = -0.01;
Double_t minB = 0.01;
Double_t maxB = 0.9;

// Define simple cuts for RP selection:
Double_t ptMinRP = 0.0; // in GeV
Double_t ptMaxRP = 10.0; // in GeV
Double_t etaMinRP  = -1.;
Double_t etaMaxRP  = 1.;
Double_t phiMinRP  = 0.0; // in degrees
Double_t phiMaxRP  = 360.0; // in degrees
Int_t pidRP = -1; // to be improved (supported eventually)

// Define simple cuts for POI selection:
Double_t ptMinPOI = 0.0; // in GeV
Double_t ptMaxPOI = 10.0; // in GeV
Double_t etaMinPOI  = -1.;
Double_t etaMaxPOI  = 1.;
Double_t phiMinPOI  = 0.0; // in degrees
Double_t phiMaxPOI  = 360.0; // in degrees
Int_t pidPOI = -1; // to be improved (supported eventually)

// Parameters for the simulation of events 'on the fly': 
//===SEED========================================================
// use always the same seed for random generators. 
// usage of same seed (kTRUE) is relevant in two cases:
// 1.) If you want to use LYZ method to calcualte differential flow;
// 2.) If you want to use phi weights for GFC, QC and FQD
Bool_t bSameSeed = kFALSE;  

//===NONFLOW=============================================
// number of times to use each track (to simulate nonflow)
Int_t iLoops = 1; 
Double_t phiRange = 0.0; // If the original track with azimuthal angle phi is splitted, 
                         // the azimuthal angle of splitted track is sampled uniformly 
                         // from (phi-phiRange, phi+phiRange). If phiRange = 0.0, the
                         // azimuthal angle of splitted track is the same as azimuthal  
                         // angle of original track. (Remark: phiRange's unit is degree.)
Double_t ptRange = 0.0; // If the original track with momentum pt is splitted, 
                        // the momentum of splitted track is sampled uniformly 
                        // from (pt-ptRange, pt+ptRange). If ptRange = 0.0, the
                        // momentum of splitted track is the same as momentum 
                        // of original track. (Remark: ptRange's unit is GeV.)
Double_t etaRange = 0.0; // If the original track with pseudorapidity eta is splitted, 
                         // the pseudorapidity of splitted track is sampled uniformly 
                         // from (eta-etaRange, eta+etaRange). If etaRange = 0.0, the
                         // pseudorapidity of splitted track will be the same as the
                         // pseudorapidity of original track.
// in addition one can simulate nonflow only in a certain detector's sector 
// ranging from nonflowSectorMin to nonflowSectorMax. Outside this sector the
// tracks will not be splitted. Use the following two settings only if iLoops>1:                         
Double_t nonflowSectorMin = 0.0; // detector's sector in which tracks will be splitted starts at this angle (in degrees)                         
Double_t nonflowSectorMax = 360.0; // detector's sector in which tracks will be splitted ends at this angle (in degrees)                        

//===GLAUBER MODEL================================================================
Bool_t bUseGlauberModel = kFALSE; // 1.) if kTRUE = multiplicity and constant flow harmonics are 
                                  //                determined event-by-event from Glauber model
                                  //                (pt and/or eta dependence of flow harmonics not supported;
                                  //                settings for multiplicity and flow harmonics bellow are irrelevant) 
                                  // 2.) if kFALSE = multiplicity and flow harmonics are determined
                                  //                 independently event-by-event with bellow settings
                                  //                 (pt and/or eta dependence of flow harmonics supported)

//===FLOW HARMONICS===============================================================
// harmonics V1, V2, V4... are constants or functions of pt and eta:         
Bool_t bPtDependentHarmonicV1 = kFALSE; 
Bool_t bEtaDependentHarmonicV1 = kFALSE;            
Bool_t bPtDependentHarmonicV2 = kFALSE; 
Bool_t bEtaDependentHarmonicV2 = kFALSE; 
Bool_t bPtDependentHarmonicV4 = kFALSE; 
Bool_t bEtaDependentHarmonicV4 = kFALSE; 
// 1.) if you use constant harmonics (bPtDependentHarmonic* = kFALSE, bEtaDependentHarmonic* = kFALSE)
//     you can additionally select if V2 will be sampled e-b-e from Gaussian or from uniform distribution
//     (constant harmonics V1 and V4 will be always sampled from Gaussian):  
Bool_t bConstantV2IsSampledFromGauss = kTRUE;
 //  1a.) if kTRUE = elliptic flow of RPs is sampled e-b-e from Gaussian distribution with
 //                  mean = dV2RP and spread = dV2SpreadRP (analogously for V1 and V4).
 //  1b.) if kFALSE = elliptic flow of RPs is sampled e-b-e uniformly from 
 //                   interval [dMinV2RP,dMaxV2RP] (not supported for V1 and V4).
 //  1c.) for a fixed elliptic flow e-b-e use Gaussian with zero spread or use uniform with dMinV2RP=dMaxV2RP 
 // V2:
 Double_t dV2RP = 0.05;       // elliptic flow of RPs (if sampled from Gaussian) 
 Double_t dV2SpreadRP = 0.0;  // elliptic flow spread of RPs (if sampled from Gaussian)
 Double_t dMinV2RP = 0.04;    // minimal elliptic flow of RPs (if sampled uniformly)
 Double_t dMaxV2RP = 0.06;    // maximal elliptic flow of RPs (if sampled uniformly)
 // V1:
 Double_t dV1RP = 0.0; // directed flow of RPs 
 Double_t dV1SpreadRP = 0.0; // directed flow spread of RPs
 // V4:
 Double_t dV4RP = 0.0; // harmonic V4 of RPs (to be improved: name needed) 
 Double_t dV4SpreadRP = 0.0; // harmonic V4's spread of RPs (to be improved: name needed)
// 2.) if you use (pt,eta) dependent harmonic V1 (bPtDependentHarmonicV1 = kTRUE or bEtaDependentHarmonicV1 = kTRUE): 
//  2a.) V1(pt) is linear up to pt = dV1PtCutOff and for pt > dV1PtCutOff it is constant, V1(pt) = dV1vsPtEtaMax:
//  2b.) V1(eta) is determined from formula V1(eta) = -eta
 Double_t dV1vsPtEtaMax = 0.10; // maximum value of V1(pt,eta) (for pt >= dV1PtCutOff and at eta = 0)
 Double_t dV1PtCutOff = 2.0; // up to this pt value V1(pt) is growing linearly versus pt
// 3.) if you use (pt,eta) dependent harmonic V2 (bPtDependentHarmonicV2 = kTRUE or bEtaDependentHarmonicV2 = kTRUE): 
//  2a.) V2(pt) is linear up to pt = dV2PtCutOff and for pt > dV2PtCutOff it is constant, V2(pt) = dV2vsPtEtaMax:
//  2b.) V2(eta) is Gaussian centered at midrapidity (eta=0) with V2(0) = dV2vsPtEtaMax and spread = dV2vsEtaSpread:
 Double_t dV2vsPtEtaMax = 0.20; // maximum value of V2(pt,eta) (for pt >= dV2PtCutOff and at eta = 0)
 Double_t dV2PtCutOff = 2.0; // up to this pt value V2(pt) is growing linearly versus pt
 Double_t dV2vsEtaSpread = 0.75; // V2(eta) is Gaussian centered at midrapidity (eta=0) with spread = dV2vsEtaSpread
// 4.) if you use (pt,eta) dependent harmonic V4 (bPtDependentHarmonicV4 = kTRUE or bEtaDependentHarmonicV4 = kTRUE):
//     V4(pt,eta) is determined as V4(pt,eta) = pow(V2(pt,eta),2.) 

//===MULTIPLICITY===============================================================
// 1.) if kTRUE  = multiplicitiy of RPs is sampled e-b-e from Gaussian distribution with
//                 mean = iMultiplicityOfRP and spread = dMultiplicitySpreadOfRP
// 2.) if kFALSE = multiplicitiy of RPs is sampled e-b-e uniformly from 
//                 interval [iMinMultOfRP,iMaxMultOfRP]
// 3.) for a fixed multiplicity use Gaussian with zero spread or use uniform with iMinMult=iMaxMult
Bool_t bMultDistrOfRPsIsGauss = kTRUE; 
                    
Int_t iMultiplicityOfRP = 500;        // mean multiplicity of RPs (if sampled from Gaussian)
Double_t dMultiplicitySpreadOfRP = 0; // multiplicity spread of RPs (if sampled from Gaussian)
Int_t iMinMultOfRP = 400;              // minimal multiplicity of RPs(if sampled uniformly)
Int_t iMaxMultOfRP = 600;             // maximal multiplicity of RPs (if sampled uniformly)
                    
//===DETECTOR ACCEPTANCE===============================================================

// 1.) if kTRUE = detectors has uniform azimuthal acceptance
// 2.) if kFALSE = you will simulate detector with non-uniform acceptance in one or two sectors. 
//                 For each of two sectors you specify phi_min, phi_max and probability p. Then all particles 
//                 going in direction phi_min < phi < phi_max will be taken with probability p. If p = 0, that
//                 sector is blocked. Set bellow phimin1, phimax1, p1 for the first sector and phimin2, phimax2, p2 
//                 for the second sector. If you set phimin2 = phimax2 = p2 = 0, only first non-uniform sector is 
//                 simulated.
Bool_t uniformAcceptance = kTRUE;
                                                                         
// settings for non-uniform acceptance:
// Remark: set the angles in degrees from interval [0,360] and probability from interval [0,1]

// 1st non-uniform sector:
Double_t phimin1 = 60;  // first non-uniform sector starts at this azimuth
Double_t phimax1 = 120; // first non-uniform sector ends at this azimuth
Double_t p1 = 0.5;     // e.g. if p1 = 0 all particles emitted in phimin1 < phi < phimax1 are blocked
                        // e.g. if p1 = 0.5 half of the particles emitted in phimin1 < phi < phimax1 are blocked

// 2nd non-uniform sector (Remark: if you do NOT want to simulate this sector, set phimin2 = phimax2 = p2 = 0):                 
Double_t phimin2 = 0.0; // second non-uniform sector starts at this azimuth (make sure phimin2 > phimax1 !!!!)
Double_t phimax2 = 0.0; // second non-uniform sector ends at this azimuth
Double_t p2 = 0.0;


//===MODIFYING Pt SPECTRA===============================================================
Double_t dTemperatureOfRP = 0.44; // 'temperature' of RPs in GeV/c (increase this parameter to get more high pt RPs) 

enum anaModes {mLocal,mLocalSource,mLocalPAR};
// mLocal: Analyze data on your computer using aliroot
// mLocalPAR: Analyze data on your computer using root + PAR files
// mLocalSource: Analyze data on your computer using root + source files
                                          
int runFlowAnalysisOnTheFly(Int_t nEvts=1000, Int_t mode=mLocal)
{
 CheckUserSettings();
 TStopwatch timer;
 timer.Start();
 
 cout<<endl;
 cout<<endl;
 cout<<"      ---- ARE YOU READY TO FLY ? ----      "<<endl;
 cout<<endl;
 
 cout<<endl;
 cout<<" ---- BEGIN FLOW ANALYSIS 'ON THE FLY' ---- "<<endl;
 cout<<endl;
 cout<<endl;
 
 LoadLibraries(mode);

 // Initialize the seed for random generator
 UInt_t sseed = 0;
 
 if(bSameSeed) 
 {
  sseed = 44; // the default constant value for seed for random generators
 } 

 if(!bSameSeed)
 {
  TTimeStamp dt;
  sseed = dt.GetNanoSec()/1000;
 }
 
 cout<<endl;
 cout<<"Seed for the random generators is "<<sseed<<endl;
 cout<<endl;
 
 //---------------------------------------------------------------------------------------
 // If the weights are used: 
 TFile *fileWithWeights = NULL;
 TList *listWithWeights = NULL;
  
 if(usePhiWeights||usePtWeights||useEtaWeights) {
   fileWithWeights = TFile::Open("weights.root","READ");
   if(fileWithWeights) {
     listWithWeights = (TList*)fileWithWeights->Get("weights");
   }
   else
     {cout << " WARNING: the file <weights.root> with weights from the previous run was not found."<<endl;
      break;
     }    
 }
 
 //---------------------------------------------------------------------------------------
 // Initialize the flowevent maker
 AliFlowEventSimpleMakerOnTheFly* eventMakerOnTheFly = new AliFlowEventSimpleMakerOnTheFly(sseed);
 eventMakerOnTheFly->Init();
 //set the range for eta subevents
 eventMakerOnTheFly->SetSubeventEtaRange(minA,maxA,minB,maxB); 
  
 //---------------------------------------------------------------------------------------
 // Initialize all the flow methods for default analysis:  
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
 AliFlowAnalysisWithMCEventPlane *mcep_AH = NULL;   

 // MCEP = monte carlo event plane
 if (MCEP) {
   AliFlowAnalysisWithMCEventPlane *mcep = new AliFlowAnalysisWithMCEventPlane();
   //mcep->SetHarmonic(2); // default is v2
   //mcep->SetFlowOfResonances(kTRUE);
   mcep->Init();
 }

 // MCEP = monte carlo event plane in another harmonic: 
 if(MCEP_AH)
 {
  AliFlowAnalysisWithMCEventPlane *mcep_ah = new AliFlowAnalysisWithMCEventPlane();
  mcep_ah->SetHarmonic(1);
  mcep_ah->Init();
 }
 
 // Mixed harmonics:
 if(MH) 
 {
  AliFlowAnalysisWithMixedHarmonics *mh = new AliFlowAnalysisWithMixedHarmonics();
  mh->SetCorrelatorInteger(1); // integer n in expression cos[n(2phi1-phi2-phi3)] = v2n*vn^2
  mh->SetMinMultiplicity(100); 
  mh->SetNoOfMultipicityBins(5);  
  mh->SetMultipicityBinWidth(200);   
  mh->Init(); 
 }
 
 // NL = nested loops:
 if(NL) {
   AliFlowAnalysisWithNestedLoops *nl = new AliFlowAnalysisWithNestedLoops();
   nl->Init();
 }

 // QC = Q-cumulants  
 if(QC) { 
   AliFlowAnalysisWithQCumulants* qc = new AliFlowAnalysisWithQCumulants();
   if(listWithWeights) qc->SetWeightsList(listWithWeights);
   if(usePhiWeights) qc->SetUsePhiWeights(usePhiWeights);
   if(usePtWeights) qc->SetUsePtWeights(usePtWeights);
   if(useEtaWeights) qc->SetUseEtaWeights(useEtaWeights);
   // qc->SetHarmonic(2); // default is v2
   // qc->SetApplyCorrectionForNUA(kTRUE); // default
   // qc->SetCalculate2DFlow(kFALSE); // default
   // qc->SetMultiplicityWeight("combinations"); // default
   // qc->SetMultiplicityWeight("unit");
   // qc->SetMultiplicityWeight("multiplicity");  
   qc->Init();  
 }
  
 // GFC = Generating Function Cumulants 
 if(GFC) {
   AliFlowAnalysisWithCumulants* gfc = new AliFlowAnalysisWithCumulants();
   if(listWithWeights) gfc->SetWeightsList(listWithWeights);
   if(usePhiWeights) gfc->SetUsePhiWeights(usePhiWeights);
   if(usePtWeights) gfc->SetUsePtWeights(usePtWeights);
   if(useEtaWeights) gfc->SetUseEtaWeights(useEtaWeights);
   gfc->Init();
 }
 
 // FQD = Fitting q-distribution 
 if(FQD) {
   AliFlowAnalysisWithFittingQDistribution* fqd = new AliFlowAnalysisWithFittingQDistribution();
   if(listWithWeights) fqd->SetWeightsList(listWithWeights);
   if(usePhiWeights) fqd->SetUsePhiWeights(usePhiWeights);  
   fqd->Init();
 }
 
 // SP = Scalar Product 
 if(SP) {
   AliFlowAnalysisWithScalarProduct* sp = new AliFlowAnalysisWithScalarProduct();
   if(listWithWeights) sp->SetWeightsList(listWithWeights);
   sp->SetUsePhiWeights(usePhiWeights);  
   sp->Init();
 }

 // LYZ1 = Lee-Yang Zeroes first run
 if(LYZ1SUM) {
   AliFlowAnalysisWithLeeYangZeros* lyz1sum = new AliFlowAnalysisWithLeeYangZeros();
   lyz1sum->SetFirstRun(kTRUE);
   lyz1sum->SetUseSum(kTRUE);
   lyz1sum->Init();
 }
 if(LYZ1PROD) {
   AliFlowAnalysisWithLeeYangZeros* lyz1prod = new AliFlowAnalysisWithLeeYangZeros();
   lyz1prod->SetFirstRun(kTRUE);
   lyz1prod->SetUseSum(kFALSE);
   lyz1prod->Init();
 }
 // LYZ2 = Lee-Yang Zeroes second run
 if(LYZ2SUM) {
   AliFlowAnalysisWithLeeYangZeros* lyz2sum = new AliFlowAnalysisWithLeeYangZeros();
   // read the input file from the first run 
   TString inputFileNameLYZ2SUM = "outputLYZ1SUManalysis.root" ;
   TFile* inputFileLYZ2SUM = new TFile(inputFileNameLYZ2SUM.Data(),"READ");
   if(!inputFileLYZ2SUM || inputFileLYZ2SUM->IsZombie()) { 
     cerr << " ERROR: To run LYZ2SUM you need the output file from LYZ1SUM. This file is not there! Please run LYZ1SUM first." << endl ;
     break; 
   }
   else { 
     TList* inputListLYZ2SUM = (TList*)inputFileLYZ2SUM->Get("cobjLYZ1SUM");  
     if (!inputListLYZ2SUM) {cout<<"Input list LYZ2SUM is NULL pointer!"<<endl; break;}
     else {
       cout<<"LYZ2SUM input file/list read..."<<endl;
       lyz2sum->SetFirstRunList(inputListLYZ2SUM);
       lyz2sum->SetFirstRun(kFALSE);
       lyz2sum->SetUseSum(kTRUE);
       lyz2sum->Init();
     }
   }
 }
 if(LYZ2PROD) {
   AliFlowAnalysisWithLeeYangZeros* lyz2prod = new AliFlowAnalysisWithLeeYangZeros();
   // read the input file from the first run 
   TString inputFileNameLYZ2PROD = "outputLYZ1PRODanalysis.root" ;
   TFile* inputFileLYZ2PROD = new TFile(inputFileNameLYZ2PROD.Data(),"READ");
   if(!inputFileLYZ2PROD || inputFileLYZ2PROD->IsZombie()) { 
     cerr << " ERROR: To run LYZ2PROD you need the output file from LYZ1PROD. This file is not there! Please run LYZ1PROD first." << endl ;
     break; 
   }
   else { 
     TList* inputListLYZ2PROD = (TList*)inputFileLYZ2PROD->Get("cobjLYZ1PROD");  
     if (!inputListLYZ2PROD) {cout<<"Input list LYZ2PROD is NULL pointer!"<<endl; break;}
     else {
       cout<<"LYZ2PROD input file/list read..."<<endl;
       lyz2prod->SetFirstRunList(inputListLYZ2PROD);
       lyz2prod->SetFirstRun(kFALSE);
       lyz2prod->SetUseSum(kFALSE);
       lyz2prod->Init();
     }
   }
 }

 // LYZEP = Lee-Yang Zeroes event plane
 if(LYZEP) {
   AliFlowLYZEventPlane* ep = new AliFlowLYZEventPlane() ;
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
 //---------------------------------------------------------------------------------------
 
 // set the global event parameters:
 eventMakerOnTheFly->SetNoOfLoops(iLoops);
 eventMakerOnTheFly->SetPhiRange(phiRange*TMath::Pi()/180.);
 eventMakerOnTheFly->SetPtRange(ptRange);
 eventMakerOnTheFly->SetEtaRange(etaRange);
 eventMakerOnTheFly->SetNonflowSectorMin(nonflowSectorMin*TMath::Pi()/180.);
 eventMakerOnTheFly->SetNonflowSectorMax(nonflowSectorMax*TMath::Pi()/180.);
 eventMakerOnTheFly->SetUseGlauberModel(bUseGlauberModel);
 
 if(!bUseGlauberModel)
 {
  if(bMultDistrOfRPsIsGauss)
  {
   eventMakerOnTheFly->SetMultDistrOfRPsIsGauss(bMultDistrOfRPsIsGauss);
   eventMakerOnTheFly->SetMultiplicityOfRP(iMultiplicityOfRP);
   eventMakerOnTheFly->SetMultiplicitySpreadOfRP(dMultiplicitySpreadOfRP);
  } else
    {
     eventMakerOnTheFly->SetMultDistrOfRPsIsGauss(bMultDistrOfRPsIsGauss);
     eventMakerOnTheFly->SetMinMultOfRP(iMinMultOfRP);
     eventMakerOnTheFly->SetMaxMultOfRP(iMaxMultOfRP); 
    }
 } // end of if(!bUseGlauberModel)
  
 eventMakerOnTheFly->SetTemperatureOfRP(dTemperatureOfRP);
 eventMakerOnTheFly->SetPtDependentHarmonicV1(bPtDependentHarmonicV1);
 eventMakerOnTheFly->SetEtaDependentHarmonicV1(bEtaDependentHarmonicV1);
 eventMakerOnTheFly->SetPtDependentHarmonicV2(bPtDependentHarmonicV2);
 eventMakerOnTheFly->SetEtaDependentHarmonicV2(bEtaDependentHarmonicV2);
 eventMakerOnTheFly->SetPtDependentHarmonicV4(bPtDependentHarmonicV4);
 eventMakerOnTheFly->SetEtaDependentHarmonicV4(bEtaDependentHarmonicV4);
 // V1:
 if(!(bPtDependentHarmonicV1 || bEtaDependentHarmonicV1))
 {
  if(!bUseGlauberModel)
  {
   eventMakerOnTheFly->SetV1RP(dV1RP);
   eventMakerOnTheFly->SetV1SpreadRP(dV1SpreadRP);  
  } 
 } else // (pt,eta) dependent V1
   {
    eventMakerOnTheFly->SetV1vsPtEtaMax(dV1vsPtEtaMax);
    if(bPtDependentHarmonicV1)
    {
     eventMakerOnTheFly->SetV1PtCutOff(dV1PtCutOff);  
    }
   }
 // V2:
 if(!(bPtDependentHarmonicV2 || bEtaDependentHarmonicV2)) // constant V2
 { 
  if(!bUseGlauberModel)
  {
   eventMakerOnTheFly->SetConstantV2IsSampledFromGauss(bConstantV2IsSampledFromGauss);
   if(bConstantV2IsSampledFromGauss) // Gauss
   {
    eventMakerOnTheFly->SetV2RP(dV2RP);
    eventMakerOnTheFly->SetV2SpreadRP(dV2SpreadRP);  
   } else // uniform
     {
      eventMakerOnTheFly->SetMinV2RP(dMinV2RP);
      eventMakerOnTheFly->SetMaxV2RP(dMaxV2RP);  
     }
  } // end of if(!bUseGlauberModel)  
 } else // (pt,eta) dependent V2
   {
    eventMakerOnTheFly->SetV2vsPtEtaMax(dV2vsPtEtaMax);
    if(bPtDependentHarmonicV2)
    {
     eventMakerOnTheFly->SetV2PtCutOff(dV2PtCutOff);  
    }
    if(bEtaDependentHarmonicV2)
    {
     eventMakerOnTheFly->SetV2vsEtaSpread(dV2vsEtaSpread);      
    }
   }  
 // V4:
 if(!(bPtDependentHarmonicV4 || bEtaDependentHarmonicV4))
 {
  if(!bUseGlauberModel)
  {
   eventMakerOnTheFly->SetV4RP(dV4RP);
   eventMakerOnTheFly->SetV4SpreadRP(dV4SpreadRP);  
  }
 } 
 
 // non-uniform acceptance:
 if(!uniformAcceptance)
 {
  eventMakerOnTheFly->SetFirstSectorPhiMin(phimin1);
  eventMakerOnTheFly->SetFirstSectorPhiMax(phimax1);
  eventMakerOnTheFly->SetFirstSectorProbability(p1);
  eventMakerOnTheFly->SetSecondSectorPhiMin(phimin2);
  eventMakerOnTheFly->SetSecondSectorPhiMax(phimax2);
  eventMakerOnTheFly->SetSecondSectorProbability(p2);
 }
 
 // simple cuts for RPs and POIs:
 AliFlowTrackSimpleCuts *cutsRP = new AliFlowTrackSimpleCuts();
 cutsRP->SetPtMax(ptMaxRP);
 cutsRP->SetPtMin(ptMinRP);
 cutsRP->SetEtaMax(etaMaxRP);
 cutsRP->SetEtaMin(etaMinRP);
 cutsRP->SetPhiMax(phiMaxRP*TMath::Pi()/180.);
 cutsRP->SetPhiMin(phiMinRP*TMath::Pi()/180.);
 cutsRP->SetPID(pidRP);
  
 AliFlowTrackSimpleCuts *cutsPOI = new AliFlowTrackSimpleCuts();
 cutsPOI->SetPtMax(ptMaxPOI);
 cutsPOI->SetPtMin(ptMinPOI);
 cutsPOI->SetEtaMax(etaMaxPOI);
 cutsPOI->SetEtaMin(etaMinPOI);
 cutsPOI->SetPhiMax(phiMaxPOI*TMath::Pi()/180.);
 cutsPOI->SetPhiMin(phiMinPOI*TMath::Pi()/180.);
 cutsPOI->SetPID(pidPOI);
                                       
 //---------------------------------------------------------------------------------------  
 // create and analyze events 'on the fly':

 for(Int_t i=0;i<nEvts;i++) {   
   // creating the event with above settings:
   AliFlowEventSimple *event = eventMakerOnTheFly->CreateEventOnTheFly(cutsRP,cutsPOI); 
   
   // analyzing the created event 'on the fly':
   // do flow analysis for various methods:
   if(MCEP)    mcep->Make(event);
   if(QC)      qc->Make(event);
   if(GFC)     gfc->Make(event);
   if(FQD)     fqd->Make(event);
   if(LYZ1SUM) lyz1sum->Make(event);
   if(LYZ1PROD)lyz1prod->Make(event);
   if(LYZ2SUM) lyz2sum->Make(event);
   if(LYZ2PROD)lyz2prod->Make(event);
   if(LYZEP)   lyzep->Make(event,ep);
   if(SP)      sp->Make(event);
   if(MH)      mh->Make(event);
   if(NL)      nl->Make(event);
   if(MCEP_AH) mcep_ah->Make(event);
   
   delete event;
 } // end of for(Int_t i=0;i<nEvts;i++)
 //---------------------------------------------------------------------------------------  

 //---------------------------------------------------------------------------------------  
 // create a new file which will hold the final results of all methods:
 TString outputFileName = "AnalysisResults.root";  
 TFile *outputFile = new TFile(outputFileName.Data(),"RECREATE");
 // create a new file for each method wich will hold list with final results:
 const Int_t nMethods = 13;
 TString method[nMethods] = {"MCEP","SP","GFC","QC","FQD","LYZ1SUM","LYZ1PROD","LYZ2SUM","LYZ2PROD","LYZEP","MH","NL","MCEP_AH"};
 TDirectoryFile *dirFileFinal[nMethods] = {NULL};
 TString fileName[nMethods]; 
 for(Int_t i=0;i<nMethods;i++)
 {
  // form a file name for each method:
  fileName[i]+="output";
  fileName[i]+=method[i].Data();
  fileName[i]+="analysis";
  dirFileFinal[i] = new TDirectoryFile(fileName[i].Data(),fileName[i].Data());
 } 
 
 // calculating and storing the final results of default flow analysis:
 if(MCEP)    {mcep->Finish();    mcep->WriteHistograms(dirFileFinal[0]);}
 if(SP)      {sp->Finish();      sp->WriteHistograms(dirFileFinal[1]);}
 if(GFC)     {gfc->Finish();     gfc->WriteHistograms(dirFileFinal[2]);}
 if(QC)      {qc->Finish();      qc->WriteHistograms(dirFileFinal[3]);}
 if(FQD)     {fqd->Finish();     fqd->WriteHistograms(dirFileFinal[4]);}
 if(LYZ1SUM) {lyz1sum->Finish(); lyz1sum->WriteHistograms(dirFileFinal[5]);}
 if(LYZ1PROD){lyz1prod->Finish();lyz1prod->WriteHistograms(dirFileFinal[6]);}
 if(LYZ2SUM) {lyz2sum->Finish(); lyz2sum->WriteHistograms(dirFileFinal[7]);}
 if(LYZ2PROD){lyz2prod->Finish();lyz2prod->WriteHistograms(dirFileFinal[8]);}
 if(LYZEP)   {lyzep->Finish();   lyzep->WriteHistograms(dirFileFinal[9]);}
 if(MH)      {mh->Finish();      mh->WriteHistograms(dirFileFinal[10]);}
 if(NL)      {nl->Finish();      nl->WriteHistograms(dirFileFinal[11]);}
 if(MCEP_AH) {mcep_ah->Finish(); mcep_ah->WriteHistograms(dirFileFinal[12]);}
 //---------------------------------------------------------------------------------------  
 
 outputFile->Close();
 delete outputFile;
 
 cout<<endl;
 cout<<endl;
 cout<<" ---- LANDED SUCCESSFULLY ---- "<<endl;
 cout<<endl; 
 
 timer.Stop();
 cout << endl;
 timer.Print();
}

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
 // Check if user settings make sense before taking off:
  
 if (LYZ1SUM && LYZ2SUM)  {cout<<"WARNING: you cannot run LYZ1 and LYZ2 at the same time! LYZ2 needs the output from LYZ1.  "<<endl; exit(); }
 if (LYZ1PROD && LYZ2PROD)  {cout<<"WARNING: you cannot run LYZ1 and LYZ2 at the same time! LYZ2 needs the output from LYZ1.  "<<endl; exit(); }
 if (LYZ2SUM && LYZEP) {cout<<"WARNING: you cannot run LYZ2 and LYZEP at the same time! LYZEP needs the output from LYZ2."<<endl; exit(); }
 if (LYZ1SUM && LYZEP) {cout<<"WARNING: you cannot run LYZ1 and LYZEP at the same time! LYZEP needs the output from LYZ2."<<endl; exit(); }

 if(!uniformAcceptance && phimin1 > phimax1)
 {
  cout<<"WARNING: you must have phimin1 < phimax1 !!!!"<<endl;
  break;
 }

 if (!uniformAcceptance && !((phimin2 == 0.) && (phimax2 == 0.) && (p2 == 0.)) && (phimin2 < phimax1 || phimin2 > phimax2))
 {
  cout<<"WARNING: you must have phimin2 > phimax1 and phimin2 < phimax2 !!!!"<<endl;
  break;
 }
 
 if((phimin1 < 0 || phimin1 > 360) || (phimax1 < 0 || phimax1 > 360) || 
    (phimin2 < 0 || phimin2 > 360) || (phimax2 < 0 || phimax2 > 360) )
 {
  cout<<"WARNING: you must take azimuthal angles from interval [0,360] !!!!"<<endl;
  break;
 }
 
 if((p1 < 0 || p1 > 1) || (p2 < 0 || p2 > 1))
 {
  cout<<"WARNING: you must take p1 and p2 from interval [0,1] !!!!"<<endl;
  break;
 }
 
 if(bUseGlauberModel)
 {
  if(bPtDependentHarmonicV1||bEtaDependentHarmonicV1||
     bPtDependentHarmonicV2||bEtaDependentHarmonicV2||
     bPtDependentHarmonicV4||bEtaDependentHarmonicV4)
  {   
   cout<<"WARNING: When using Glauber model pt and/or eta dependence of flow harmonics is not supported !!!!"<<endl;
   cout<<"         Set all booleans bPtDependentHarmonic* to kFALSE in the macro."<<endl;
   exit(0); 
  }
 }

 if(bPtDependentHarmonicV4 && !bPtDependentHarmonicV2))
 {
  cout<<"WARNING: V4(pt,eta) is determined as pow(V2(pt,eta),2.) !!!!"<<endl;
  cout<<"         You must also set bPtDependentHarmonicV2 = kTRUE in the macro."<<endl;
  exit(0);
 }
   
 if(bEtaDependentHarmonicV4 && !bEtaDependentHarmonicV2))
 {
  cout<<"WARNING: V4(pt,eta) is determined as pow(V2(pt,eta),2.) !!!!"<<endl;
  cout<<"         You must also set bEtaDependentHarmonicV2 = kTRUE in the macro."<<endl;
  exit(0); 
 }
  
 if(bPtDependentHarmonicV4 && (bEtaDependentHarmonicV4 != bEtaDependentHarmonicV2))
 {
  cout<<"WARNING: V4(pt,eta) is determined as pow(V2(pt,eta),2.) !!!!"<<endl;
  cout<<"         You must also have bEtaDependentHarmonicV2 = bEtaDependentHarmonicV4 in the macro."<<endl;
  exit(0); 
 }
 
 if(bEtaDependentHarmonicV4 && (bPtDependentHarmonicV4 != bPtDependentHarmonicV2))
 {
  cout<<"WARNING: V4(pt,eta) is determined as pow(V2(pt,eta),2.) !!!!"<<endl;
  cout<<"         You must also have bPtDependentHarmonicV2 = bPtDependentHarmonicV4 in the macro."<<endl;
  exit(0); 
 }
  
} // end of void CheckUserSettings()

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
    gROOT->LoadMacro("AliFlowCommon/AliFlowCumuConstants.cxx+");
    
    // Flow event
    gROOT->LoadMacro("AliFlowCommon/AliFlowVector.cxx+"); 
    gROOT->LoadMacro("AliFlowCommon/AliFlowTrackSimple.cxx+");
    gROOT->LoadMacro("AliFlowCommon/AliFlowEvent.cxx+");
    gROOT->LoadMacro("AliFlowCommon/AliFlowEventSimple.cxx+");
    
    // Cuts
    gROOT->LoadMacro("AliFlowCommon/AliFlowTrackSimpleCuts.cxx+");    
    
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


