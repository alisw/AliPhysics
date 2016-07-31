//--------------------------------------------------------------------------
// Base macro for submitting muon Resolution analysis.
//
// It needs to load magnetic field, mapping, geometry (+alignment), and reconstruction parameters from the OCDB
// It is done by the task and OCDB storage locations can be parameterized in MuonResolution.C (default is raw://)
//
// The task reads ESDs
// Intermediate results are stored in a file chamberResolution_step<#step>.root
// Final results are stored in the file results.root
//
// Author: Philippe Pillot - SUBATECH Nantes
//--------------------------------------------------------------------------

//______________________________________________________________________________
void RunMuonResolution(TString smode = "local", TString inputFileName = "AliESDs.root",
		       TString rootVersion = "", TString aliphysicsVersion = "vAN-20160524-1", Int_t nSteps = 5,
		       Bool_t selectPhysics = kTRUE, Bool_t selectTrigger = kTRUE, Bool_t matchTrig = kTRUE,
		       Bool_t applyAccCut = kTRUE, Bool_t applyPDCACut = kTRUE, Double_t minMomentum = 0., Double_t minPt = 0.,
                       Bool_t isMC = kFALSE, Bool_t correctForSystematics = kTRUE, Int_t extrapMode = 1,
                       Bool_t shiftHalfCh = kFALSE, Bool_t shiftDE = kFALSE, Int_t nevents = 1234567890)
{
  /// Compute the cluster resolution by studying cluster-track residual, deconvoluting from track resolution
  /// - smode = "local" or "caf"
  /// - inputFileName = an ESD root file or a list of ESDs or a collection of ESDs or a dataset in proof mode
  /// - rootVersion = version of root package to enable on AAF (only used in proof mode)
  /// - alirootVersion = version of aliroot package to enable on AAF (only used in proof mode)
  /// - aliphysicsVersion = version of aliphysics package to enable on AAF (only used in proof mode)
  /// - nSteps = number of times to task is run (at each step it starts with the chamber resolution obtained in the previous one)
  /// - selectPhysics : apply or not the physics selection
  /// - selectTrigger : select only muon trigger events or not (the type of trigger can eventually be changed)
  /// - matchTrigger : select only the tracks matching the trigger or not
  /// - applyAccCut : select only the tracks passing the Rabs and the eta cut or not
  /// - minMomentum : select only the tracks with a total momentum above this value
  /// - if correctForSystematics == kTRUE: the systematic shifts of the residuals is included in the resolution
  /// - if extrapMode == 0: extrapolate from the closest cluster
  /// - if extrapMode == 1: extrapolate from the previous cluster except between stations 2-3-4
  /// - nevents = maximum number of processed events
  
  if (smode == "saf3" && gSystem->GetFromPipe("hostname") != "nansafmaster3.in2p3.fr") {
    
    // run on SAF3
    if (!RunAnalysisOnSAF3(aliphysicsVersion, inputFileName)) return;
    
    // draw the results locally
    TFile* outFile = TFile::Open("results.root","READ");
    if (!outFile || !outFile->IsOpen()) return;
    outFile->FindObjectAny("convergenceRes")->Draw();
    outFile->FindObjectAny("convergenceShift")->Draw();
    outFile->Close();
    
  } else {
    
    gROOT->ProcessLine(".include $ALICE_ROOT/include");
    gROOT->ProcessLine(".include $ALICE_PHYSICS/include");
    
    // compile analysis macro locally
    if (smode == "saf3") gROOT->LoadMacro("MuonResolution.C++g");
    else gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/MUON/dep/MuonResolution.C++g");
    MuonResolution(smode, inputFileName, rootVersion, aliphysicsVersion, nSteps, selectPhysics, selectTrigger, matchTrig, applyAccCut, applyPDCACut, minMomentum, minPt, isMC, correctForSystematics, extrapMode, shiftHalfCh, shiftDE, nevents);
    
  }
  
}

//______________________________________________________________________________
Bool_t RunAnalysisOnSAF3(TString aliphysicsVersion, TString dataset)
{
  /// Run analysis on SAF3
  
  // --- mount nansafmaster3 ---
  TString saf3dir = gSystem->ExpandPathName("$HOME/saf3");
  if (gSystem->AccessPathName(saf3dir.Data())) gSystem->Exec(Form("mkdir %s", saf3dir.Data()));
  if (gSystem->AccessPathName(Form("%s/.vaf", saf3dir.Data()))) {
    Int_t ret = gSystem->Exec(Form("sshfs -o ssh_command=\"gsissh -p1975\" nansafmaster3.in2p3.fr: %s", saf3dir.Data()));
    if (ret != 0) {
      cout<<"mounting of saf3 folder failed"<<endl;
      return kFALSE;
    }
  }
  
  // --- create the executable to run on SAF3 ---
  CreateSAF3Executable(dataset);
  
  // --- copy files needed for this analysis ---
  if (!CopyFileOnSAF3(dataset)) {
    cout << "cp problem" << endl;
    return kFALSE;
  }
  
  // --- change the AliPhysics version on SAF3 ---
  gSystem->Exec(Form("sed -i '' 's/VafAliPhysicsVersion.*/VafAliPhysicsVersion=\"%s\"/g' $HOME/saf3/.vaf/vaf.conf", aliphysicsVersion.Data()));
  
  // --- enter SAF3 and run analysis ---
  TString analysisLocation = gSystem->pwd();
  analysisLocation.ReplaceAll(Form("%s/", gSystem->Getenv("HOME")), "");
  gSystem->Exec(Form("gsissh -p 1975 -t -Y nansafmaster3.in2p3.fr 'cd %s; ~/saf3-enter \"\" \"./runAnalysis.sh 2>&1 | tee runAnalysis.log; exit\"'", analysisLocation.Data()));
  
  // --- copy analysis results (assuming analysis run smootly) ---
  gSystem->Exec(Form("cp -p %s/%s/*.root .", saf3dir.Data(), analysisLocation.Data()));
  
  return kTRUE;
  
}

//______________________________________________________________________________
Bool_t CopyFileOnSAF3(TString dataset)
{
  /// Copy files needed for this analysis from current to saf3 directory
  
  TString saf3dir = gSystem->ExpandPathName("$HOME/saf3");
  if (gSystem->AccessPathName(Form("%s/.vaf", saf3dir.Data()))) {
    cout<<"saf3 folder is not mounted"<<endl;
    cout<<"please retry as it can take some time to mount it"<<endl;
    return kFALSE;
  }
  
  TString remoteLocation = gSystem->pwd();
  remoteLocation.ReplaceAll(gSystem->Getenv("HOME"),saf3dir.Data());
  if (gSystem->AccessPathName(remoteLocation.Data())) gSystem->Exec(Form("mkdir -p %s", remoteLocation.Data()));
  
  gSystem->Exec(Form("cp -p $ALICE_PHYSICS/PWGPP/MUON/dep/RunMuonResolution.C %s/RunMuonResolution.C", remoteLocation.Data()));
  gSystem->Exec(Form("cp -p $ALICE_PHYSICS/PWGPP/MUON/dep/MuonResolution.C %s/MuonResolution.C", remoteLocation.Data()));
  gSystem->Exec(Form("cp -p $ALICE_PHYSICS/PWGPP/MUON/dep/AddTaskMuonResolution.C %s/AddTaskMuonResolution.C", remoteLocation.Data()));
  //gSystem->Exec(Form("cp -p $ALICE_PHYSICS/PARfiles/PWGPPMUONdep.par %s/PWGPPMUONdep.par", remoteLocation.Data()));
  gSystem->Exec(Form("cp runAnalysis.sh %s/runAnalysis.sh", remoteLocation.Data()));
  
  if (dataset.EndsWith(".txt")) {
    gSystem->Exec(Form("cat %s | awk {'print $1\";Mode=cache\"}' > datasetSaf3.txt", dataset.Data()));
    gSystem->Exec(Form("cp datasetSaf3.txt %s/datasetSaf3.txt", remoteLocation.Data()));
  } else if (dataset.EndsWith(".root"))
    gSystem->Exec(Form("cp %s %s/%s", dataset.Data(), remoteLocation.Data(), gSystem->BaseName(dataset.Data())));
  
  return kTRUE;
  
}

//______________________________________________________________________________
void CreateSAF3Executable(TString dataset)
{
  /// Create the executable to run on SAF3
  ofstream outFile("runAnalysis.sh");
  outFile << "#!/bin/bash" << endl;
  outFile << "vafctl start" << endl;
  Int_t nWorkers = 88;
  outFile << "nWorkers=" << nWorkers << endl;
  outFile << "let \"nWorkers -= `pod-info -n`\"" << endl;
  outFile << "echo \"requesting $nWorkers additional workers\"" << endl;
  outFile << "vafreq $nWorkers" << endl;
  outFile << "vafwait " << nWorkers << endl;
  TString macro = gSystem->GetFromPipe("tail -n 1 $HOME/.root_hist | sed 's/(.*)//g;s/^\.x\ //g;s:^.*/::g'");
  TString arg = gSystem->GetFromPipe("tail -n 1 $HOME/.root_hist | sed 's/.*(/(/g'");
  if (dataset.EndsWith(".txt")) arg.ReplaceAll(dataset.Data(), "datasetSaf3.txt");
  else if (dataset.EndsWith(".root")) arg.ReplaceAll(dataset.Data(), gSystem->BaseName(dataset.Data()));
  outFile << "root -b -q '" << macro.Data() << arg.Data() << "'" << endl;
  outFile << "vafctl stop" << endl;
  outFile.close();
  gSystem->Exec("chmod u+x runAnalysis.sh");
}

