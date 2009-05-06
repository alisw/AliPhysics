#include "TStopwatch.h"
#include "TObjArray"
#include "Riostream.h"
#include "TFile.h"

//--------------------------------------------------------------------------------------
// RUN SETTINGS
// flow analysis method can be: (set to kTRUE or kFALSE)
Bool_t SP    = kTRUE;
Bool_t LYZ1  = kTRUE;
Bool_t LYZ2  = kFALSE;
Bool_t LYZEP = kFALSE;
Bool_t GFC   = kTRUE;
Bool_t QC    = kTRUE;
Bool_t FQD   = kTRUE;
Bool_t MCEP  = kTRUE; 
//--------------------------------------------------------------------------------------

Bool_t bSameSeed = kFALSE; // use always the same seed for random generators
Bool_t bConstantHarmonics = kFALSE; // harmonics V1, V2, V4... are constant (kTRUE) or functions of pt and eta (kFALSE)

// Set the event parameters:
Int_t iLoops = 1; // number of times to use each track (to simulate nonflow)

Int_t iMultiplicityOfRP = 500; // multiplicity of RPs
Double_t dMultiplicitySpreadOfRP = 0; // multiplicity spread of RPs
Double_t dTemperatureOfRP = 0.44; // 'temperature' of RPs in GeV/c (increase this parameter to get more high pt RPs) 

//......................................................................................  
// if you use (pt,eta) dependent harmonics (bConstantHarmonics = kFALSE):
Double_t dPtCutOff = 2.0; // V2(pt) is linear up to pt = 2 GeV and for pt > 2 GeV it is constant: V2(pt) = dVRPMax
Double_t dV2RPMax = 0.20; // maximum value of V2(pt) for pt >= 2GeV
//...................................................................................... 

//......................................................................................  
// if you use constant harmonics (bConstantHarmonics = kTRUE):
Double_t dV2RP = 0.05; // elliptic flow of RPs
Double_t dV2SpreadRP = 0.; // elliptic flow spread of RPs

Double_t dV1RP = 0.0; // directed flow of RPs
Double_t dV1SpreadRP = 0.0; // directed flow spread of RPs

Double_t dV4RP = 0.0; // harmonic V4 of RPs (to be improved: name needed)
Double_t dV4SpreadRP = 0.0; // harmonic V4's spread of RPs (to be improved: name needed)
//......................................................................................  

enum anaModes {mLocal,mLocalSource,mLocalPAR};
// mLocal: Analyze data on your computer using aliroot
// mLocalPAR: Analyze data on your computer using root + PAR files
// mLocalSource: Analyze data on your computer using root + source files
                                          
int runFlowAnalysisOnTheFly(Int_t mode=mLocal, Int_t nEvts=100)
{
 TStopwatch timer;
 timer.Start();
 
 if (LYZ1 && LYZ2)  {cout<<"WARNING: you cannot run LYZ1 and LYZ2 at the same time! LYZ2 needs the output from LYZ1.  "<<endl; exit(); }
 if (LYZ2 && LYZEP) {cout<<"WARNING: you cannot run LYZ2 and LYZEP at the same time! LYZEP needs the output from LYZ2."<<endl; exit(); }
 if (LYZ1 && LYZEP) {cout<<"WARNING: you cannot run LYZ1 and LYZEP at the same time! LYZEP needs the output from LYZ2."<<endl; exit(); }
 
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
 
 //---------------------------------------------------------------------------------------
 // Initialize the flowevent maker
 AliFlowEventSimpleMakerOnTheFly* eventMakerOnTheFly = new AliFlowEventSimpleMakerOnTheFly(sseed);
 eventMakerOnTheFly->Init();
  
 //---------------------------------------------------------------------------------------
 // Initialize all the flow methods:  
 AliFlowAnalysisWithQCumulants    *qc    = NULL;
 AliFlowAnalysisWithCumulants     *gfc   = NULL;
 AliFittingQDistribution          *fqd   = NULL;
 AliFlowAnalysisWithLeeYangZeros  *lyz1  = NULL;
 AliFlowAnalysisWithLeeYangZeros  *lyz2  = NULL;
 AliFlowAnalysisWithLYZEventPlane *lyzep = NULL;
 AliFlowAnalysisWithScalarProduct *sp    = NULL;
 AliFlowAnalysisWithMCEventPlane  *mcep  = NULL;   

 // MCEP = monte carlo event plane
 if (MCEP) {
   AliFlowAnalysisWithMCEventPlane *mcep = new AliFlowAnalysisWithMCEventPlane();
   mcep->Init();
 }

  // QC = Q-cumulants  
 if(QC) { 
   AliFlowAnalysisWithQCumulants* qc = new AliFlowAnalysisWithQCumulants();
   qc->Init();
 }
  
 // GFC = Generating Function Cumulants 
 if(GFC) {
   AliFlowAnalysisWithCumulants* gfc = new AliFlowAnalysisWithCumulants();
   gfc->Init();
 }
 
 // FQD = Fitting q-distribution 
 if(FQD) {
   AliFittingQDistribution* fqd = new AliFittingQDistribution();
   fqd->Init();
 }
 
 // SP = Scalar Product 
 if(SP) {
   AliFlowAnalysisWithScalarProduct* sp = new AliFlowAnalysisWithScalarProduct();
   sp->Init();
 }

 // LYZ1 = Lee-Yang Zeroes first run
 if(LYZ1) {
   AliFlowAnalysisWithLeeYangZeros* lyz1 = new AliFlowAnalysisWithLeeYangZeros();
   lyz1->SetFirstRun(kTRUE);
   lyz1->SetUseSum(kTRUE);
   lyz1->Init();
 }

 // LYZ2 = Lee-Yang Zeroes second run
 if(LYZ2) {
   AliFlowAnalysisWithLeeYangZeros* lyz2 = new AliFlowAnalysisWithLeeYangZeros();
   // read the input file from the first run 
   TString inputFileNameLYZ2 = "outputLYZ1analysis.root" ;
   TFile* inputFileLYZ2 = new TFile(inputFileNameLYZ2.Data(),"READ");
   if(!inputFileLYZ2 || inputFileLYZ2->IsZombie()) { 
     cerr << " ERROR: NO First Run file... " << endl ;
     break; 
   }
   else { 
     TList* inputListLYZ2 = (TList*)inputFileLYZ2->Get("cobjLYZ1");  
     if (!inputListLYZ2) {cout<<"Input list is NULL pointer!"<<endl; break;}
     else {
       cout<<"LYZ2 input file/list read..."<<endl;
       lyz2->SetFirstRunList(inputListLYZ2);
       lyz2->SetFirstRun(kFALSE);
       lyz2->SetUseSum(kTRUE);
       lyz2->Init();
     }
   }
 }
 
 // LYZEP = Lee-Yang Zeroes event plane
 if(LYZEP) {
   AliFlowLYZEventPlane* ep = new AliFlowLYZEventPlane() ;
   AliFlowAnalysisWithLYZEventPlane* lyzep = new AliFlowAnalysisWithLYZEventPlane();
   // read the input file from the second lyz run 
   TString inputFileNameLYZEP = "outputLYZ2analysis.root" ;
   TFile* inputFileLYZEP = new TFile(inputFileNameLYZEP.Data(),"READ");
   if(!inputFileLYZEP || inputFileLYZEP->IsZombie()) { 
     cerr << " ERROR: NO Second Run file... " << endl ; 
     break;
   }
   else { 
     TList* inputListLYZEP = (TList*)inputFileLYZEP->Get("cobjLYZ2");  
     if (!inputListLYZEP) {cout<<"Input list is NULL pointer!"<<endl; break;}
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
 eventMakerOnTheFly->SetMultiplicityOfRP(iMultiplicityOfRP);
 eventMakerOnTheFly->SetMultiplicitySpreadOfRP(dMultiplicitySpreadOfRP);
 eventMakerOnTheFly->SetTemperatureOfRP(dTemperatureOfRP);

 eventMakerOnTheFly->SetV1RP(dV1RP);
 eventMakerOnTheFly->SetV1SpreadRP(dV1SpreadRP);  
 eventMakerOnTheFly->SetV4RP(dV4RP);
 eventMakerOnTheFly->SetV4SpreadRP(dV4SpreadRP);  
 
 // constant harmonic V2:
 if(bConstantHarmonics)
 { 
  eventMakerOnTheFly->SetUseConstantHarmonics(bConstantHarmonics);
  eventMakerOnTheFly->SetV2RP(dV2RP);
  eventMakerOnTheFly->SetV2SpreadRP(dV2SpreadRP);  
 }
 // (pt,eta) dependent harmonic V2:
 if(!bConstantHarmonics)
 {
  eventMakerOnTheFly->SetUseConstantHarmonics(bConstantHarmonics);
  eventMakerOnTheFly->SetV2RPMax(dV2RPMax);
  eventMakerOnTheFly->SetPtCutOff(dPtCutOff);  
 }
       
 //---------------------------------------------------------------------------------------  
 // create and analyze events 'on the fly':

 for(Int_t i=0;i<nEvts;i++) {   
   // creating the event with above settings:
   AliFlowEventSimple *event = eventMakerOnTheFly->CreateEventOnTheFly(); 
   
   // analyzing the created event 'on the fly':
   // do flow analysis for various methods:
   if(MCEP) mcep->Make(event);
   if(QC) qc->Make(event);
   if(GFC) gfc->Make(event);
   if(FQD) fqd->Make(event);
   if(LYZ1) lyz1->Make(event);
   if(LYZ2) lyz2->Make(event);
   if(LYZEP) lyzep->Make(event,ep);
   if(SP) sp->Make(event);
   
   delete event;
 } // end of for(Int_t i=0;i<nEvts;i++)
 //---------------------------------------------------------------------------------------  



 //---------------------------------------------------------------------------------------  
 // calculating and storing the final results of flow analysis
 if(MCEP) {mcep->Finish(); mcep->WriteHistograms("outputMCEPanalysis.root");}
 if(SP) {sp->Finish(); sp->WriteHistograms("outputSPanalysis.root");}
 if(QC) {qc->Finish(); qc->WriteHistograms("outputQCanalysis.root");}
 if(GFC) {gfc->Finish(); gfc->WriteHistograms("outputGFCanalysis.root");}
 if(FQD) {fqd->Finish(); fqd->WriteHistograms("outputFQDanalysis.root");}
 if(LYZ1) {lyz1->Finish(); lyz1->WriteHistograms("outputLYZ1analysis.root");}
 if(LYZ2) {lyz2->Finish(); lyz2->WriteHistograms("outputLYZ2analysis.root");}
 if(LYZEP) {lyzep->Finish(); lyzep->WriteHistograms("outputLYZEPanalysis.root");}
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
    gROOT->LoadMacro("AliFlowCommon/AliFittingFunctionsForQDistribution.cxx+");
    gROOT->LoadMacro("AliFlowCommon/AliFlowLYZEventPlane.cxx+");
    
    // Flow Analysis code for various methods
    gROOT->LoadMacro("AliFlowCommon/AliFlowAnalysisWithMCEventPlane.cxx+"); 
    gROOT->LoadMacro("AliFlowCommon/AliFlowAnalysisWithScalarProduct.cxx+");
    gROOT->LoadMacro("AliFlowCommon/AliFlowAnalysisWithLYZEventPlane.cxx+");
    gROOT->LoadMacro("AliFlowCommon/AliFlowAnalysisWithLeeYangZeros.cxx+");
    gROOT->LoadMacro("AliFlowCommon/AliFlowAnalysisWithCumulants.cxx+");
    gROOT->LoadMacro("AliFlowCommon/AliFlowAnalysisWithQCumulants.cxx+"); 
    gROOT->LoadMacro("AliFlowCommon/AliFittingQDistribution.cxx+");
    
    // Class to fill the FlowEvent on the fly (generate Monte Carlo events)
    gROOT->LoadMacro("AliFlowCommon/AliFlowEventSimpleMakerOnTheFly.cxx+");   
    
    cout << "finished loading macros!" << endl;  
    
  }  
  
}


