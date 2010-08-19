//example script on what to do with the star events
//run e.g. like this:
//                    root runStarFlowAnalysis.C

//--------------------------------------------------------------------------------------
// Run flow analysis on local data with custom FlowEvent maker
// RUN SETTINGS
//flow analysis method can be: (set to kTRUE or kFALSE)
Bool_t SP       = kTRUE;
Bool_t LYZ1SUM  = kTRUE;
Bool_t LYZ1PROD = kTRUE;
Bool_t LYZ2SUM  = kFALSE; 
Bool_t LYZ2PROD = kFALSE;
Bool_t LYZEP    = kFALSE; 
Bool_t GFC      = kTRUE;
Bool_t QC       = kTRUE;
Bool_t FQD      = kTRUE;
Bool_t MH       = kTRUE; 
Bool_t NL       = kFALSE; 
Bool_t MCEP     = kFALSE; //does not work yet 24/12/08
//--------------------------------------------------------------------------------------

// Weights 
// Use weights for Q vector
Bool_t usePhiWeights = kFALSE; //Phi (correction for non-uniform azimuthal acceptance)
Bool_t usePtWeights  = kFALSE; //v'(pt) (differential flow in pt)
Bool_t useEtaWeights = kFALSE; //v'(eta) (differential flow in eta)


void  runStarFlowAnalysis()
{
  gSystem->Load("libTree.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libPWG2flowCommon");

  Int_t maxNumberOfEvents = 1000000;

  TStopwatch timer;
  timer.Start();
  
  if (LYZ1SUM && LYZ2SUM) {cout<<"WARNING: you cannot run LYZ1 and LYZ2 at the same time! LYZ2 needs the output from LYZ1."<<endl; exit(1); }
  if (LYZ1PROD && LYZ2PROD) {cout<<"WARNING: you cannot run LYZ1 and LYZ2 at the same time! LYZ2 needs the output from LYZ1."<<endl; exit(1); }
  if (LYZ2SUM && LYZEP) {cout<<"WARNING: you cannot run LYZ2 and LYZEP at the same time! LYZEP needs the output from LYZ2."<<endl; exit(1); }
  if (LYZ1SUM && LYZEP) {cout<<"WARNING: you cannot run LYZ1 and LYZEP at the same time! LYZEP needs the output from LYZ2."<<endl; exit(1); }



  //define reference particles
  AliStarTrackCuts* rpCuts = AliStarTrackCuts::StandardCuts();

  //define particles of interest
  AliStarTrackCuts* poiCuts = AliStarTrackCuts::StandardCuts();
  poiCuts->SetPtMin(0.05);

  //define event cuts
  AliStarEventCuts* starEventCuts = AliStarEventCuts::StandardCuts();
  starEventCuts-> SetCentralityIDMax(8);
  starEventCuts-> SetCentralityIDMin(8);

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

  //flow methods:  
  AliFlowAnalysisWithQCumulants *qc = NULL;
  AliFlowAnalysisWithCumulants *gfc = NULL;
  AliFlowAnalysisWithFittingQDistribution *fqd = NULL;
  AliFlowAnalysisWithLeeYangZeros *lyz1sum = NULL;
  AliFlowAnalysisWithLeeYangZeros *lyz1prod = NULL;
  AliFlowAnalysisWithLeeYangZeros *lyz2sum = NULL;
  AliFlowAnalysisWithLeeYangZeros *lyz2prod = NULL;
  AliFlowAnalysisWithLYZEventPlane *lyzep = NULL;
  AliFlowAnalysisWithScalarProduct *sp = NULL;
  AliFlowAnalysisWithMCEventPlane *mcep = NULL;     
  AliFlowAnalysisWithMixedHarmonics *mh = NULL;
  AliFlowAnalysisWithNestedLoops *nl = NULL;

  //MCEP = monte carlo event plane
  if (MCEP) {
    AliFlowAnalysisWithMCEventPlane *mcep = new AliFlowAnalysisWithMCEventPlane();
    mcep->Init();
  }

  //QC = Q-cumulants  
  if(QC) { 
    AliFlowAnalysisWithQCumulants* qc = new AliFlowAnalysisWithQCumulants();
    if(listWithWeights) qc->SetWeightsList(listWithWeights);
    if(usePhiWeights) qc->SetUsePhiWeights(usePhiWeights);
    if(usePtWeights) qc->SetUsePtWeights(usePtWeights);
    if(useEtaWeights) qc->SetUseEtaWeights(useEtaWeights);
    qc->Init();
  }
  
  //GFC = Generating Function Cumulants 
  if(GFC) {
    AliFlowAnalysisWithCumulants* gfc = new AliFlowAnalysisWithCumulants();
    if(listWithWeights) gfc->SetWeightsList(listWithWeights);
    if(usePhiWeights) gfc->SetUsePhiWeights(usePhiWeights);
    if(usePtWeights) gfc->SetUsePtWeights(usePtWeights);
    if(useEtaWeights) gfc->SetUseEtaWeights(useEtaWeights);
    gfc->Init();
  }
  
  //FQD = Fitting q-distribution 
  if(FQD) {
    AliFlowAnalysisWithFittingQDistribution* fqd = new AliFlowAnalysisWithFittingQDistribution();
    if(listWithWeights) fqd->SetWeightsList(listWithWeights);
    if(usePhiWeights) fqd->SetUsePhiWeights(usePhiWeights);
    fqd->Init();
  }

  //SP = Scalar Product 
  if(SP) {
    AliFlowAnalysisWithScalarProduct* sp = new AliFlowAnalysisWithScalarProduct();
    if(usePhiWeights) sp->SetUsePhiWeights(usePhiWeights);
    sp->Init();
  }

  //LYZ1 = Lee-Yang Zeroes first run
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
  //LYZ2 = Lee-Yang Zeroes second run
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
      if (!inputListLYZ2SUM) {cout<<"SUM Input list is NULL pointer!"<<endl; break;}
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
      if (!inputListLYZ2PROD) {cout<<"PROD Input list is NULL pointer!"<<endl; break;}
      else {
	cout<<"LYZ2PROD input file/list read..."<<endl;
	lyz2prod->SetFirstRunList(inputListLYZ2PROD);
	lyz2prod->SetFirstRun(kFALSE);
	lyz2prod->SetUseSum(kTRUE);
	lyz2prod->Init();
      }
    }
  }
 //LYZEP = Lee-Yang Zeroes event plane
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
  // MH = Mixed Harmonics:  
  if(MH) { 
    AliFlowAnalysisWithMixedHarmonics* mh = new AliFlowAnalysisWithMixedHarmonics();
    if(listWithWeights) mh->SetWeightsList(listWithWeights);
    //if(usePhiWeights) mh->SetUsePhiWeights(usePhiWeights); // to be improved (enabled)
    //if(usePtWeights) mh->SetUsePtWeights(usePtWeights); // to be improved (enabled)
    //if(useEtaWeights) mh->SetUseEtaWeights(useEtaWeights); // to be improved (enabled)
    mh->Init();
  }
  // NL = Nested Loops:  
  if(NL) { 
    AliFlowAnalysisWithNestedLoops* nl = new AliFlowAnalysisWithNestedLoops();
    if(listWithWeights) nl->SetWeightsList(listWithWeights);
    //if(usePhiWeights) nl->SetUsePhiWeights(usePhiWeights); // to be improved (enabled)
    //if(usePtWeights) nl->SetUsePtWeights(usePtWeights); // to be improved (enabled)
    //if(useEtaWeights) nl->SetUseEtaWeights(useEtaWeights); // to be improved (enabled)
    nl->Init();
  }

  //------------------------------------------------------------------------


  Int_t i=0;
  AliStarEventReader starReader("/Users/snelling/alice_data/jthomas/testData/") ;
  while ( starReader.GetNextEvent() )                                // Get next event
  {
    AliStarEvent* starEvent = starReader.GetEvent();
    if ( !starEventCuts->PassesCuts(starEvent) ) continue;              // Test if the event is good

    AliFlowEventSimple* flowEvent = new AliFlowEventStar(starEvent,rpCuts,poiCuts);  // make a flow event from a star event (aka "the magic")

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

    //////////////////////////////////

    //starEvent->Print("all");
    //flowEvent->Print();

    delete flowEvent;

    i++;
    if (i>maxNumberOfEvents) break;
  }

  // create a new file which will hold the final results of all methods:
  TString outputFileName = "AnalysisResults.root";  
  TFile *outputFile = new TFile(outputFileName.Data(),"RECREATE");
  // create a new file for each method wich will hold list with final results:
  const Int_t nMethods = 12;
  TString method[nMethods] = {"MCEP","SP","GFC","QC","FQD","LYZ1SUM","LYZ1PROD","LYZ2SUM","LYZ2PROD","LYZEP","MH","NL"};
  TDirectoryFile *dirFileFinal[nMethods] = {NULL};
  TString fileNameMethod[nMethods]; 
  for(Int_t i=0;i<nMethods;i++)
    {
      // form a file name for each method:
      fileNameMethod[i]+="output";
      fileNameMethod[i]+=method[i].Data();
      fileNameMethod[i]+="analysis";
      dirFileFinal[i] = new TDirectoryFile(fileNameMethod[i].Data(),fileNameMethod[i].Data());
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
  //---------------------------------------------------------------------------------------  
  
  outputFile->Close();
  delete outputFile;
  
  cout << endl;
  cout << " Fini ... " << endl;
  cout << endl;
  
  timer.Stop();
  cout << endl;
  timer.Print();
  
  delete rpCuts;
  delete poiCuts;
  delete starEventCuts;
}
