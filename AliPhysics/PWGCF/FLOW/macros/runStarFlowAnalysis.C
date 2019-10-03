//example script on what to do with the star events
//run e.g. like this:
//                    root runStarFlowAnalysis.C
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
Bool_t MH       = kFALSE; // mixed harmonics 
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


void  runStarFlowAnalysis(const char* inputDataFiles="/Users/snelling/alice_data/jthomas/testData/")
{
  gSystem->Load("libTree");
  gSystem->Load("libVMC");
  gSystem->Load("libPhysics");
  gSystem->Load("libPWGflowBase");

  Int_t maxNumberOfEvents = 1000000;

  TStopwatch timer;
  timer.Start();
  
  if (LYZ1SUM && LYZ2SUM) {cout<<"WARNING: you cannot run LYZ1 and LYZ2 at the same time! LYZ2 needs the output from LYZ1."<<endl; exit(1); }
  if (LYZ1PROD && LYZ2PROD) {cout<<"WARNING: you cannot run LYZ1 and LYZ2 at the same time! LYZ2 needs the output from LYZ1."<<endl; exit(1); }
  if (LYZ2SUM && LYZEP) {cout<<"WARNING: you cannot run LYZ2 and LYZEP at the same time! LYZEP needs the output from LYZ2."<<endl; exit(1); }
  if (LYZ1SUM && LYZEP) {cout<<"WARNING: you cannot run LYZ1 and LYZEP at the same time! LYZEP needs the output from LYZ2."<<endl; exit(1); }



  //define reference particles
  AliStarTrackCuts* rpCuts = AliStarTrackCuts::StandardCuts();
  rpCuts->SetPtMin(0.05);
  rpCuts->SetPtMax(10.);

  //define particles of interest
  AliStarTrackCuts* poiCuts = AliStarTrackCuts::StandardCuts();
  poiCuts->SetPtMin(0.05);
  poiCuts->SetPtMax(10.);

  //define event cuts
  AliStarEventCuts* starEventCuts = AliStarEventCuts::StandardCuts();
  starEventCuts-> SetCentralityIDMax(3);
  starEventCuts-> SetCentralityIDMin(3);

  //if the weights are used: 
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
  AliFlowAnalysisWithMCEventPlane *mcep_ah = NULL;   

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
  mh->SetHarmonic(1); // integer n in expression cos[n(2phi1-phi2-phi3)] = v2n*vn^2
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
   qc->SetnBinsMult(10000);
   qc->SetMinMult(0);
   qc->SetMaxMult(10000);      
   qc->Init();  
 }
  
 // GFC = Generating Function Cumulants 
 if(GFC) {
   AliFlowAnalysisWithCumulants* gfc = new AliFlowAnalysisWithCumulants();
   if(listWithWeights) gfc->SetWeightsList(listWithWeights);
   if(usePhiWeights) gfc->SetUsePhiWeights(usePhiWeights);
   if(usePtWeights) gfc->SetUsePtWeights(usePtWeights);
   if(useEtaWeights) gfc->SetUseEtaWeights(useEtaWeights);
   // calculation vs multiplicity:
   gfc->SetCalculateVsMultiplicity(kFALSE);   
   gfc->SetnBinsMult(10000);
   gfc->SetMinMult(0);
   gfc->SetMaxMult(10000);   
   // tuning of interpolating parameters:
   gfc->SetTuneParameters(kFALSE);
   Double_t r0[10] = {1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7}; // up to 10 values allowed
   for(Int_t r=0;r<10;r++) {gfc->SetTuningR0(r0[r],r);}
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


  Int_t i=0;
  AliStarEventReader starReader(inputDataFiles) ;
  while ( starReader.GetNextEvent() )                                // Get next event
  {
    AliStarEvent* starEvent = starReader.GetEvent();
    if ( !starEventCuts->PassesCuts(starEvent) ) continue;              // Test if the event is good

    AliFlowEventSimple* flowEvent = new AliFlowEventStar(starEvent,rpCuts,poiCuts);  // make a flow event from a star event (aka "the magic")
    flowEvent->TagSubeventsInEta(minA, maxA, minB, maxB );

    /////analysis here////////////////

    // do flow analysis for various methods
    if(MCEP)    mcep->Make(flowEvent);
    if(QC)      qc->Make(flowEvent);
    if(GFC)     gfc->Make(flowEvent);
    if(FQD)     fqd->Make(flowEvent);
    if(LYZ1SUM) lyz1sum->Make(flowEvent);
    if(LYZ1PROD)lyz1prod->Make(flowEvent);
    if(LYZ2SUM) lyz2sum->Make(flowEvent);
    if(LYZ2PROD)lyz2prod->Make(flowEvent);
    if(LYZEP)   lyzep->Make(flowEvent,ep);
    if(SP)      sp->Make(flowEvent);	      
    if(MH)      mh->Make(flowEvent);
    if(NL)      nl->Make(flowEvent);
    if(MCEP_AH) mcep_ah->Make(flowEvent);

    //////////////////////////////////

    flowEvent->Print();

    delete flowEvent;

    i++;
    if (i>maxNumberOfEvents) break;
  }

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
  
  delete rpCuts;
  delete poiCuts;
  delete starEventCuts;
}
