//--------------------------------------------------------------------------
// Base macro for submitting muon Resolution analysis.
//
// In case it is not run with full aliroot, it needs the following libraries:
//  - libSTEERBase.so
//  - libESD.so
//  - libAOD.so
//  - libANALYSIS.so
//  - libANALYSISalice.so
//  - libGui.so
//  - libMinuit.so
//  - libProofPlayer.so
//  - libXMLParser.so
//  - libRAWDatabase.so
//  - libCDB.so
//  - libSTEER.so
//  - libMUONcore.so
//  - libMUONmapping.so
//  - libMUONcalib.so
//  - libMUONgeometry.so
//  - libMUONtrigger.so
//  - libMUONraw.so
//  - libMUONbase.so
//  - libMUONrec.so
//  - libCORRFW.so
//  - libPWG3muondep.so
//
// The macro reads ESDs and store outputs in standard output file (AnalysisResults.root)
// It also needs to load magnetic field, mapping, geometry (+alignment), and recoonstruction parameters from the OCDB
//
// Author: Philippe Pillot - SUBATECH Nantes
//--------------------------------------------------------------------------

TString aliroot="VO_ALICE@AliRoot::v4-19-15-AN";

enum {kLocal, kInteractif_xml, kInteractif_ESDList, kProof};

// global declaration needed in Proof mode, certainly because of the interpretor
TStopwatch* localTimer;
TMultiGraph* mgClusterResXVsStep;
TMultiGraph* mgClusterResYVsStep;
Double_t clusterResNB[10];
Double_t clusterResB[10];
Double_t clusterResNBErr[10];
Double_t clusterResBErr[10];

void RunMuonResolution(TString smode = "local", TString inputFileName = "AliESDs.root", Int_t nSteps = 10,
		       Bool_t selectPhysics = kFALSE, Bool_t matchTrig = kFALSE, Double_t minMomentum = 0.,
		       Bool_t correctForSystematics = kTRUE, Int_t extrapMode = 1)
{
  /// Compute the cluster resolution by studying cluster-track residual, deconvoluting from track resolution
  /// if extrapMode == 0: extrapolate from the closest cluster
  /// if extrapMode == 1: extrapolate from the previous cluster except between stations 2-3-4
  /// if correctForSystematics == kTRUE: the systematic shifts of the residuals is included in the resolution
  
  // timer start...
  localTimer = new TStopwatch;
  
  // check parameters
  nSteps = TMath::Max(nSteps,1);
  if (extrapMode != 0 && extrapMode != 1) {
    Error("RunMuonResolution","incorrect extrapolation mode!");
    return;
  }
  
  // Check runing mode
  Int_t mode = GetMode(smode, inputFileName);
  if(mode < 0){
    Error("RunMuonResolution","Please provide either an ESD root file a collection of ESDs or a dataset.");
    return;
  }
  
  gSystem->Load("libVMC");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  
  // Load additional libraries
  gSystem->Load("libProofPlayer");
  TString extraLibs="Physics:Minuit:XMLParser:Gui:RAWDatabase:CDB:STEER:MUONcore:MUONmapping:MUONcalib:MUONgeometry:MUONtrigger:MUONraw:MUONbase:MUONrec:CORRFW:PWG3muondep";
  TObjArray* libs = extraLibs.Tokenize(":");
  for (Int_t i = 0; i < libs->GetEntriesFast(); i++)
    gSystem->Load(Form("lib%s",static_cast<TObjString*>(libs->UncheckedAt(i))->GetName()));
  
  // check for old output file to removed
  char remove = '';
  if (!gSystem->Exec("ls chamberResolution_step*[0-9].root")) {
    cout<<"above files must be removed from the current directory. Delete? [y=yes, n=no] "<<flush;
    while (remove != 'y' && remove != 'n') cin>>remove;
    if (remove == 'y') gSystem->Exec("rm -f chamberResolution_step*[0-9].root");
    else {
      Error("RunMuonResolution","cannot proceed with these files there otherwise results will be mixed up!");
      return;
    }
  }
  
  // OCDB access
  if(mode == kLocal) AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  else if (mode != kProof) {
    if (!TGrid::Connect("alien://")) return;
    if(mode == kInteractif_ESDList || !gSystem->AccessPathName("ConfigureCuts.C"))
      AliCDBManager::Instance()->SetDefaultStorage("raw://");
  }
  
  // Create input object
  TObject* inputObj = 0x0;
  if (mode == kProof) inputObj = new TObjString(inputFileName);
  else inputObj = CreateChain(mode, inputFileName);
  if (!inputObj) return;
  
  // set starting chamber resolution (if -1 they will be loaded from recoParam in the task)
  for (Int_t i = 0; i < AliMUONConstants::NTrackingCh(); i++) {
    clusterResNB[i] = -1.;
    clusterResB[i] = -1.;
    clusterResNBErr[i] = 0.;
    clusterResBErr[i] = 0.;
  }
  
  // output graphs
  mgClusterResXVsStep = new TMultiGraph("mgClusterResXVsStep","cluster X-resolution versus step;step;#sigma_{X} (cm)");
  mgClusterResYVsStep = new TMultiGraph("mgClusterResYVsStep","cluster Y-resolution versus step;step;#sigma_{Y} (cm)");
  TGraphErrors* gClusterResXVsStep[10];
  TGraphErrors* gClusterResYVsStep[10];
  for (Int_t i = 0; i < AliMUONConstants::NTrackingCh(); i++) {
    gClusterResXVsStep[i] = new TGraphErrors(nSteps+1);
    gClusterResXVsStep[i]->SetName(Form("gResX_ch%d",i+1));
    gClusterResXVsStep[i]->SetMarkerStyle(kFullDotMedium);
    gClusterResXVsStep[i]->SetMarkerColor(i+1+i/9);
    mgClusterResXVsStep->Add(gClusterResXVsStep[i],"lp");
    
    gClusterResYVsStep[i] = new TGraphErrors(nSteps+1);
    gClusterResYVsStep[i]->SetName(Form("gResY_ch%d",i+1));
    gClusterResYVsStep[i]->SetMarkerStyle(kFullDotMedium);
    gClusterResYVsStep[i]->SetMarkerColor(i+1+i/9);
    mgClusterResYVsStep->Add(gClusterResYVsStep[i],"lp");
  }
  
  // loop over step
  for (Int_t iStep = 0; iStep < nSteps; iStep++) {
    cout<<"step "<<iStep+1<<"/"<<nSteps<<endl;
    
    // Connect to proof if needed and prepare environment
    if (mode == kProof) {
      
      // set general environment and close previous session
      if (iStep == 0) {
	gEnv->SetValue("XSec.GSI.DelegProxy","2");
	TProof::AddEnvVar("ALIROOT_EXTRA_LIBS", extraLibs.Data());
      } else gProof->Close("s");
      
      // connect
      if (gSystem->Getenv("alien_API_USER") == NULL) TProof::Open("alice-caf.cern.ch","workers=20");
      else TProof::Open(Form("%s@alice-caf.cern.ch",gSystem->Getenv("alien_API_USER")),"workers=20");
      if (!gProof) return;
      
      // set environment and compile task on workers
      gProof->EnablePackage(aliroot.Data());
      
      // prepare OCDB access on workers
      gProof->Exec("AliCDBManager::Instance()->SetDefaultStorage(\"raw://\")");
      
    }
    
    // run one step
    run(mode, iStep, inputObj, extrapMode, correctForSystematics, minMomentum, matchTrig, selectPhysics, clusterResNB, clusterResB);
    
    // fill graph with starting resolutions from the task at first step
    if (iStep == 0) for (Int_t i = 0; i < AliMUONConstants::NTrackingCh(); i++) {
      gClusterResXVsStep[i]->SetPoint(0, 0, clusterResNB[i]);
      gClusterResXVsStep[i]->SetPointError(0, 0., clusterResNBErr[i]);
      gClusterResYVsStep[i]->SetPoint(0, 0, clusterResB[i]);
      gClusterResYVsStep[i]->SetPointError(0, 0., clusterResBErr[i]);
    }
    
    // read the chamber resolution from the output file
    if (!GetChamberResolution(iStep, clusterResNB, clusterResB, clusterResNBErr, clusterResBErr)) return;
    
    // fill graphs with computed resolutions
    for (Int_t i = 0; i < AliMUONConstants::NTrackingCh(); i++) {
      gClusterResXVsStep[i]->SetPoint(iStep+1, iStep+1, clusterResNB[i]);
      gClusterResXVsStep[i]->SetPointError(iStep+1, 0., clusterResNBErr[i]);
      gClusterResYVsStep[i]->SetPoint(iStep+1, iStep+1, clusterResB[i]);
      gClusterResYVsStep[i]->SetPointError(iStep+1, 0., clusterResBErr[i]);
    }
    
  }
  
  // copy final results in results.root file
  gSystem->Exec(Form("cp chamberResolution_step%d.root results.root", nSteps-1));
  
  // display convergence
  TCanvas* convergence = new TCanvas("convergence","convergence");
  convergence->Divide(1,2);
  convergence->cd(1);
  mgClusterResXVsStep->Draw("ap");
  convergence->cd(2);
  mgClusterResYVsStep->Draw("ap");
  
  // save convergence plots
  TFile* outFile = TFile::Open("results.root","UPDATE");
  if (!outFile || !outFile->IsOpen()) return;
  outFile->cd();
  mgClusterResXVsStep->Write();
  mgClusterResYVsStep->Write();
  convergence->Write();
  outFile->Close();
  
  // print results
  printf("\nchamber resolution:\n");
  printf(" - non-bending:");
  for (Int_t i = 0; i < AliMUONConstants::NTrackingCh(); i++) printf((i==0)?" %5.3f":", %5.3f",clusterResNB[i]);
  printf("\n -     bending:");
  for (Int_t i = 0; i < AliMUONConstants::NTrackingCh(); i++) printf((i==0)?" %6.4f":", %6.4f",clusterResB[i]);
  printf("\n\n");
  
  // ...timer stop
  localTimer->Stop();
  localTimer->Print();
  
}

//______________________________________________________________________________
void run(Int_t mode, Int_t iStep, TObject* input, Int_t extrapMode, Bool_t correctForSystematics,
	 Double_t minMomentum, Bool_t matchTrig, Bool_t selectPhysics, Double_t clusterResNB[10], Double_t clusterResB[10])
{
  /// launch the analysis with these parameters
  
  // Create the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("MuonResolutionAnalysis");
  
  // ESD input handler
  AliESDInputHandler* esdH = new AliESDInputHandler();
  esdH->SetReadFriends(kFALSE);
  mgr->SetInputEventHandler(esdH);
  
  // event selection
  if (selectPhysics) {
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
    AliPhysicsSelectionTask* physicsSelection = AddTaskPhysicsSelection();
    if(!physicsSelection) {
      Error("RunMuonResolution","AliPhysicsSelectionTask not created!");
      return;
    }
    physicsSelection->GetPhysicsSelection()->SetUseMuonTriggers();
  }
  
  // Muon Resolution analysis
  gROOT->LoadMacro("$ALICE_ROOT/PWG3/muondep/AddTaskMuonResolution.C");
  AliAnalysisManager::SetCommonFileName(Form("chamberResolution_step%d.root", iStep));
  AliAnalysisTaskMuonResolution* muonResolution = AddTaskMuonResolution(selectPhysics, matchTrig, minMomentum, correctForSystematics, extrapMode);
  if(!muonResolution) {
    Error("RunMuonResolution","AliAnalysisTaskMuonResolution not created!");
    return;
  }
  muonResolution->SetStartingResolution(clusterResNB, clusterResB);
  
  // Enable debug printouts
  //mgr->SetDebugLevel(2);
  
  // start local analysis
  if (mgr->InitAnalysis()) {
    mgr->PrintStatus();
    if (mode == kProof) mgr->StartAnalysis("proof", Form("%s",static_cast<TObjString*>(input)->GetName()));
    else mgr->StartAnalysis("local", static_cast<TChain*>(input));
  }
  
  // save the summary canvases
  if (muonResolution->GetCanvases()) {
    TFile* outFile = TFile::Open(Form("chamberResolution_step%d.root", iStep),"UPDATE");
    if (outFile && outFile->IsOpen()) {
      muonResolution->GetCanvases()->Write();
      outFile->Close();
    }
  }
  
  // return starting chamber resolution from the task
  muonResolution->GetStartingResolution(clusterResNB, clusterResB);
  
  // clean memory
  delete mgr;
  TObject::SetObjectStat(kFALSE);
}

//______________________________________________________________________________
Bool_t GetChamberResolution(Int_t iStep, Double_t clusterResNB[10], Double_t clusterResB[10], Double_t clusterResNBErr[10], Double_t clusterResBErr[10])
{
  /// read the chamber resolution from the output file
  
  TFile* outFile = TFile::Open(Form("chamberResolution_step%d.root", iStep),"READ");
  
  if (!outFile || !outFile->IsOpen()) {
    Error("GetChamberResolution","output file does not exist!");
    return kFALSE;
  }
  
  TObjArray* summary = static_cast<TObjArray*>(outFile->FindObjectAny("ChamberRes"));
  TGraphErrors* gCombinedResidualXPerChSigma = (summary) ? static_cast<TGraphErrors*>(summary->FindObject("gCombinedResidualXPerChSigma")) : 0x0;
  TGraphErrors* gCombinedResidualYPerChSigma = (summary) ? static_cast<TGraphErrors*>(summary->FindObject("gCombinedResidualYPerChSigma")) : 0x0;
  
  if (!gCombinedResidualXPerChSigma || !gCombinedResidualYPerChSigma) {
    Error("GetChamberResolution","resolution graphs do not exist!");
    return kFALSE;
  }
  
  Double_t dummy;
  for (Int_t i = 0; i < AliMUONConstants::NTrackingCh(); i++) {
    gCombinedResidualXPerChSigma->GetPoint(i, dummy, clusterResNB[i]);
    gCombinedResidualYPerChSigma->GetPoint(i, dummy, clusterResB[i]);
    clusterResNBErr[i] = gCombinedResidualXPerChSigma->GetErrorY(i);
    clusterResBErr[i] = gCombinedResidualYPerChSigma->GetErrorY(i);
  }
  
  outFile->Close();
  
  return kTRUE;
}

//______________________________________________________________________________
Int_t GetMode(TString smode, TString input)
{
  if (smode == "local") {
    if ( input.EndsWith(".xml") ) return kInteractif_xml;
    else if ( input.EndsWith(".txt") ) return kInteractif_ESDList;
    else if ( input.EndsWith(".root") ) return kLocal;    
  } else if (smode == "proof") return kProof;
  return -1;
}

//______________________________________________________________________________
TChain* CreateChainFromCollection(const char *xmlfile)
{
  // Create a chain from the collection of tags.
  TAlienCollection* coll = TAlienCollection::Open(xmlfile);
  if (!coll) {
    Error("CreateChainFromCollection", Form("Cannot create an AliEn collection from %s!", xmlfile));
    return NULL;
  }
  
  TGridResult* tagResult = coll->GetGridResult("",kFALSE,kFALSE);
  AliTagAnalysis *tagAna = new AliTagAnalysis("ESD");
  tagAna->ChainGridTags(tagResult);
  
  AliRunTagCuts      *runCuts = new AliRunTagCuts();
  AliLHCTagCuts      *lhcCuts = new AliLHCTagCuts();
  AliDetectorTagCuts *detCuts = new AliDetectorTagCuts();
  AliEventTagCuts    *evCuts  = new AliEventTagCuts();
  
  // Check if the cuts configuration file was provided
  if (!gSystem->AccessPathName("ConfigureCuts.C")) {
    gROOT->LoadMacro("ConfigureCuts.C");
    ConfigureCuts(runCuts, lhcCuts, detCuts, evCuts);
  }
  
  TChain *chain = tagAna->QueryTags(runCuts, lhcCuts, detCuts, evCuts);
  if (!chain || !chain->GetNtrees()) return NULL;
  chain->ls();
  return chain;
}

//______________________________________________________________________________
TChain* CreateChainFromFile(const char *rootfile)
{
  // Create a chain using the root file.
  TChain* chain = new TChain("esdTree");
  chain->Add(rootfile);
  if (!chain->GetNtrees()) return NULL;
  chain->ls();
  return chain;
}

//______________________________________________________________________________
TChain* CreateChainFromESDList(const char *esdList)
{
  // Create a chain using tags from the run list.
  TChain* chain = new TChain("esdTree");
  ifstream inFile(esdList);
  TString inFileName;
  if (inFile.is_open()) {
    while (! inFile.eof() ) {
      inFileName.ReadLine(inFile,kFALSE);
      if(!inFileName.EndsWith(".root")) continue;
      chain->Add(inFileName.Data());
    }
  }
  inFile.close();
  if (!chain->GetNtrees()) return NULL;
  chain->ls();
  return chain;
}

//______________________________________________________________________________
TChain* CreateChain(Int_t mode, TString input)
{
  printf("*******************************\n");
  printf("*** Getting the Chain       ***\n");
  printf("*******************************\n");
  if(mode == kInteractif_xml) return CreateChainFromCollection(input.Data());
  else if (mode == kInteractif_ESDList) return CreateChainFromESDList(input.Data());
  else if (mode == kLocal) return CreateChainFromFile(input.Data());
  else return NULL;
}

