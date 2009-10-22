#include "TStopwatch.h"
#include "TObjArray"
#include "Riostream.h"
#include "TFile.h"

//--------------------------------------------------------------------------------------
// RUN SETTINGS
// flow analysis method can be: (set to kTRUE or kFALSE)
Bool_t SP       = kTRUE;
Bool_t LYZ1SUM  = kTRUE;
Bool_t LYZ1PROD = kTRUE;
Bool_t LYZ2SUM  = kFALSE;
Bool_t LYZ2PROD = kFALSE;
Bool_t LYZEP    = kFALSE;
Bool_t GFC      = kTRUE;
Bool_t QC       = kTRUE;
Bool_t FQD      = kTRUE;
Bool_t MCEP     = kTRUE;
//--------------------------------------------------------------------------------------

// Weights 
// use weights for Q-vector:
Bool_t usePhiWeights = kFALSE; // phi weights (correction for non-uniform azimuthal acceptance)
Bool_t usePtWeights  = kFALSE; // pt weights 
Bool_t useEtaWeights = kFALSE; // eta weights

// Run same flow analysis method but with different settings/aims
// You will have to label each setting/aim with your own label (see examples bellow): 
Bool_t GFC_Additional_Analysis = kFALSE;
Bool_t QC_Additional_Analysis  = kFALSE;
Bool_t FQD_Additional_Analysis = kFALSE;

// Parameters for the simulation of events 'on the fly': 
Bool_t bSameSeed = kFALSE; // use always the same seed for random generators. 
                           // usage of same seed (kTRUE) is relevant in two cases:
                           // 1.) If you want to use LYZ method to calcualte differential flow;
                           // 2.) If you want to use phi weights for GFC, QC and FQD
                           
Bool_t bConstantHarmonics = kTRUE; // harmonics V1, V2, V4... are constant (kTRUE) or functions of pt and eta (kFALSE)

Int_t iLoops = 1; // number of times to use each track (to simulate nonflow)

Bool_t bMultDistrOfRPsIsGauss = kTRUE; // 1.) if kTRUE  = multiplicitiy of RPs is sampled e-b-e from Gaussian distribution with
                                        //                 mean = iMultiplicityOfRP and spread = dMultiplicitySpreadOfRP
                                        // 2.) if kFALSE = multiplicitiy of RPs is sampled e-b-e uniformly from 
                                        //                 interval [iMinMultOfRP,iMaxMultOfRP]
                                        // 3.) for a fixed multiplicity use Gaussian with zero spread or use uniform with iMinMult=iMaxMult
                                        
Bool_t bV2DistrOfRPsIsGauss = kTRUE; // 1.) if kTRUE  = elliptic flow of RPs is sampled e-b-e from Gaussian distribution with
                                      //                 mean = dV2RP and spread = dV2SpreadRP
                                      // 2.) if kFALSE = elliptic flow of RPs is sampled e-b-e uniformly from 
                                      //                 interval [dMinV2RP,dMaxV2RP]
                                      // 3.) for a fixed elliptic flow use Gaussian with zero spread or use uniform with dMinV2RP=dMaxV2RP

Bool_t uniformAcceptance = kTRUE; // 1.) if kTRUE = detectors has uniform azimuthal acceptance
                                  // 2.) if kFALSE = you will simulate detector with non-uniform acceptance in one or two sectors. 
                                  //                 For each of two sectors you specify phi_min, phi_max and probability p. Then all particles 
                                  //                 going in direction phi_min < phi < phi_max will be taken with probability p. If p = 0, that
                                  //                 sector is blocked. Set bellow phimin1, phimax1, p1 for the first sector and phimin2, phimax2, p2 
                                  //                 for the second sector. If you set phimin2 = phimax2 = p2 = 0, only first non-uniform sector is 
                                  //                 simulated.
                                                                                                                                                                                                                                                          
Int_t iMultiplicityOfRP = 500;        // mean multiplicity of RPs (if sampled from Gaussian)
Double_t dMultiplicitySpreadOfRP = 0; // multiplicity spread of RPs (if sampled from Gaussian)
Int_t iMinMultOfRP = 400;             // minimal multiplicity of RPs (if sampled uniformly)
Int_t iMaxMultOfRP = 600;             // maximal multiplicity of RPs (if sampled uniformly)

Double_t dTemperatureOfRP = 0.44; // 'temperature' of RPs in GeV/c (increase this parameter to get more high pt RPs) 

//......................................................................................  
// if you use (pt,eta) dependent harmonics (bConstantHarmonics = kFALSE):
Double_t dPtCutOff = 2.0; // V2(pt) is linear up to pt = 2 GeV and for pt > 2 GeV it is constant: V2(pt) = dVRPMax
Double_t dV2RPMax = 0.20; // maximum value of V2(pt) for pt >= 2GeV
//...................................................................................... 

//......................................................................................  
// if you use constant harmonics (bConstantHarmonics = kTRUE) (i.e. no pt dependence):
Double_t dV2RP = 0.05;       // elliptic flow of RPs (if sampled from Gaussian)
Double_t dV2SpreadRP = 0.0;  // elliptic flow spread of RPs (if sampled from Gaussian)
Double_t dMinV2RP = 0.04;    // minimal elliptic flow of RPs (if sampled uniformly)
Double_t dMaxV2RP = 0.06;    // maximal elliptic flow of RPs (if sampled uniformly)

Double_t dV1RP = 0.0; // directed flow of RPs
Double_t dV1SpreadRP = 0.0; // directed flow spread of RPs

Double_t dV4RP = 0.0; // harmonic V4 of RPs (to be improved: name needed)
Double_t dV4SpreadRP = 0.0; // harmonic V4's spread of RPs (to be improved: name needed)
//......................................................................................  

//......................................................................................  
// settings for non-uniform acceptance:
// Remark: set the angles in degrees from interval [0,360] and probability from interval [0,1]

// 1st non-uniform sector:
Double_t phimin1 = 60;  // first non-uniform sector starts at this azimuth
Double_t phimax1 = 120; // first non-uniform sector ends at this azimuth
Double_t p1 = 0.33;        // e.g. if p1 = 0 all particles emitted in phimin1 < phi < phimax1 are blocked
                        // e.g. if p1 = 0.5 half of the particles emitted in phimin1 < phi < phimax1 are blocked

// 2nd non-uniform sector (Remark: if you do NOT want to simulate this sector, set phimin2 = phimax2 = p2 = 0):                 
Double_t phimin2 = 0.0; // second non-uniform sector starts at this azimuth (make sure phimin2 > phimax1 !!!!)
Double_t phimax2 = 0.0; // second non-uniform sector ends at this azimuth
Double_t p2 = 0.0;
//......................................................................................  

enum anaModes {mLocal,mLocalSource,mLocalPAR};
// mLocal: Analyze data on your computer using aliroot
// mLocalPAR: Analyze data on your computer using root + PAR files
// mLocalSource: Analyze data on your computer using root + source files
                                          
int runFlowAnalysisOnTheFly(Int_t mode=mLocal, Int_t nEvts=440)
{
 TStopwatch timer;
 timer.Start();
 
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
 AliFlowAnalysisWithMCEventPlane *mcep = NULL;   

 // MCEP = monte carlo event plane
 if (MCEP) {
   AliFlowAnalysisWithMCEventPlane *mcep = new AliFlowAnalysisWithMCEventPlane();
   // mcep->SetHarmonic(2); // default is v2
   mcep->Init();
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
   qc->SetEvaluateNestedLoopsForIntFlow(kFALSE);
   qc->SetEvaluateNestedLoopsForDiffFlow(kFALSE);
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
 
 //---------------------------------------------------------------------------------------
 // Initialize all the flow methods for additional analysis with different settings/aims
 // Label each setting/aim with different label !!!!
 
 // GFC:    
 TString gfcDefaultName = "outputGFCanalysis"; 
 // 1.) GFC analysis for elliptic flow with r0 = 1.5: 
 AliFlowAnalysisWithCumulants *gfc_1;
 TString gfcAnalysisLabels_1 = "_r0_1.5"; // all histograms and output file name will have this label
 TString gfcOutputFileName_1;
 gfcOutputFileName_1 = gfcDefaultName.Data();
 gfcOutputFileName_1 += gfcAnalysisLabels_1.Data();
 gfcOutputFileName_1 += ".root";
 if(GFC_Additional_Analysis)
 {
  gfc_1 = new AliFlowAnalysisWithCumulants();
  //gfc_1->SetAnalysisLabel(gfcAnalysisLabels_1.Data());
  if(listWithWeights) gfc_1->SetWeightsList(listWithWeights);
  if(usePhiWeights) gfc_1->SetUsePhiWeights(usePhiWeights);
  if(usePtWeights) gfc_1->SetUsePtWeights(usePtWeights);
  if(useEtaWeights) gfc_1->SetUseEtaWeights(useEtaWeights);
  gfc_1->Init();
 }
 
 // QC:
 TString qcDefaultName = "outputQCanalysis"; 
 // 1.) QC analysis for directed flow: 
 AliFlowAnalysisWithQCumulants *qc_1;
 TString qcAnalysisLabels_1 = "_v1"; // all histograms and output file name will have this label
 TString qcOutputFileName_1;
 qcOutputFileName_1 = qcDefaultName.Data();
 qcOutputFileName_1 += qcAnalysisLabels_1.Data();
 qcOutputFileName_1 += ".root";
 if(QC_Additional_Analysis)
 {
  qc_1 = new AliFlowAnalysisWithQCumulants();
  //qc_1->SetAnalysisLabel(qcAnalysisLabels_1->Data());
  if(listWithWeights) qc_1->SetWeightsList(listWithWeights);
  if(usePhiWeights) qc_1->SetUsePhiWeights(usePhiWeights);
  qc_1->Init();
 }
 
 // FQD:
 TString fqdDefaultName = "outputFQDanalysis";
 // 1.) FQD fitting with fixed sigma:
 AliFlowAnalysisWithFittingQDistribution *fqd_1;
 TString fqdAnalysisLabels_1 = "_fixedSigma"; // all histograms and output file name will have this label
 TString fqdOutputFileName_1;
 fqdOutputFileName_1 = fqdDefaultName.Data();
 fqdOutputFileName_1 += fqdAnalysisLabels_1.Data();
 fqdOutputFileName_1 += ".root";
 if(FQD_Additional_Analysis)
 { 
  fqd_1 = new AliFlowAnalysisWithFittingQDistribution();
  //fqd_1->SetAnalysisLabel(fqdAnalysisLabels_1->Data());
  if(listWithWeights) fqd_1->SetWeightsList(listWithWeights);
  if(usePhiWeights) fqd_1->SetUsePhiWeights(usePhiWeights);
  fqd_1->Init();
 }
 //---------------------------------------------------------------------------------------
  
 // set the global event parameters: 
 eventMakerOnTheFly->SetNoOfLoops(iLoops);
 
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
  
 eventMakerOnTheFly->SetTemperatureOfRP(dTemperatureOfRP);

 eventMakerOnTheFly->SetV1RP(dV1RP);
 eventMakerOnTheFly->SetV1SpreadRP(dV1SpreadRP);  
 eventMakerOnTheFly->SetV4RP(dV4RP);
 eventMakerOnTheFly->SetV4SpreadRP(dV4SpreadRP);  
 
 // constant harmonic V2:
 if(bConstantHarmonics)
 { 
  eventMakerOnTheFly->SetUseConstantHarmonics(bConstantHarmonics);
  if(bV2DistrOfRPsIsGauss)
  {
   eventMakerOnTheFly->SetV2DistrOfRPsIsGauss(bV2DistrOfRPsIsGauss);
   eventMakerOnTheFly->SetV2RP(dV2RP);
   eventMakerOnTheFly->SetV2SpreadRP(dV2SpreadRP);  
  } else
    {
     eventMakerOnTheFly->SetV2DistrOfRPsIsGauss(bV2DistrOfRPsIsGauss);
     eventMakerOnTheFly->SetMinV2RP(dMinV2RP);
     eventMakerOnTheFly->SetMaxV2RP(dMaxV2RP);  
    }
 }
 
 // (pt,eta) dependent harmonic V2:
 if(!bConstantHarmonics)
 {
  eventMakerOnTheFly->SetUseConstantHarmonics(bConstantHarmonics);
  eventMakerOnTheFly->SetV2RPMax(dV2RPMax);
  eventMakerOnTheFly->SetPtCutOff(dPtCutOff);  
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
       
 //---------------------------------------------------------------------------------------  
 // create and analyze events 'on the fly':

 for(Int_t i=0;i<nEvts;i++) {   
   // creating the event with above settings:
   AliFlowEventSimple *event = eventMakerOnTheFly->CreateEventOnTheFly(); 
   
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
   
   if(GFC_Additional_Analysis)
   {
    // r0 = 1.5:
    gfc_1->Make(event);
   }
   if(QC_Additional_Analysis)
   {
    // v1:
    qc_1->Make(event);
   }
   if(FQD_Additional_Analysis)
   {
    // fixed sigma:
    fqd_1->Make(event);
   }
   
   delete event;
 } // end of for(Int_t i=0;i<nEvts;i++)
 //---------------------------------------------------------------------------------------  



 //---------------------------------------------------------------------------------------  
 // calculating and storing the final results of default flow analysis:
 if(MCEP)    {mcep->Finish();    mcep->WriteHistograms("outputMCEPanalysis.root");}
 if(SP)      {sp->Finish();      sp->WriteHistograms("outputSPanalysis.root");}
 if(QC)      {qc->Finish();      qc->WriteHistograms("outputQCanalysis.root");}
 if(GFC)     {gfc->Finish();     gfc->WriteHistograms("outputGFCanalysis.root");}
 if(FQD)     {fqd->Finish();     fqd->WriteHistograms("outputFQDanalysis.root");}
 if(LYZ1SUM) {lyz1sum->Finish(); lyz1sum->WriteHistograms("outputLYZ1SUManalysis.root");}
 if(LYZ1PROD){lyz1prod->Finish();lyz1prod->WriteHistograms("outputLYZ1PRODanalysis.root");}
 if(LYZ2SUM) {lyz2sum->Finish(); lyz2sum->WriteHistograms("outputLYZ2SUManalysis.root");}
 if(LYZ2PROD){lyz2prod->Finish();lyz2prod->WriteHistograms("outputLYZ2PRODanalysis.root");}
 if(LYZEP)   {lyzep->Finish();   lyzep->WriteHistograms("outputLYZEPanalysis.root");}
 //---------------------------------------------------------------------------------------  
 
 //---------------------------------------------------------------------------------------  
 // calculating and storing the final results of flow analysis with different settings/aims:
 if(GFC_Additional_Analysis)
 {
  // r0 = 1.5:
  gfc_1->Finish();
  gfc_1->WriteHistograms(gfcOutputFileName_1.Data());
 }
 if(QC_Additional_Analysis)
 {
  // v1:
  qc_1->Finish();
  qc_1->WriteHistograms(qcOutputFileName_1.Data());
 }
 if(FQD_Additional_Analysis)
 {
  // fixed sigma:
  fqd_1->Finish();
  fqd_1->WriteHistograms(fqdOutputFileName_1.Data());
 }
 //---------------------------------------------------------------------------------------
 
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

void LoadLibraries(const anaModes mode) {
  
  //--------------------------------------
  // Load the needed libraries most of them already loaded by aliroot
  //--------------------------------------
  gSystem->Load("libTree.so");
  gSystem->Load("libGeom.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libXMLIO.so");
  gSystem->Load("libPhysics.so");
  
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
    gSystem->Load("libCORRFW.so");
    cerr<<"libCORRFW.so loaded..."<<endl;
    gSystem->Load("libPWG2flowCommon.so");
    cerr<<"libPWG2flowCommon.so loaded..."<<endl;
    gSystem->Load("libPWG2flowTasks.so");
    cerr<<"libPWG2flowTasks.so loaded..."<<endl;
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
    
    // Class to fill the FlowEvent on the fly (generate Monte Carlo events)
    gROOT->LoadMacro("AliFlowCommon/AliFlowEventSimpleMakerOnTheFly.cxx+");   
    
    cout << "finished loading macros!" << endl;  
    
  }  
  
}


