//--------------------------------------------------------------------------
// Macro compiled and launch by RunMuonResolution.C for submitting muon Resolution analysis locally or on CAF.
// See RunMuonResolution.C for more details
//
// Author: Philippe Pillot - SUBATECH Nantes
//--------------------------------------------------------------------------

#if !defined(__CINT__) || defined(__MAKECINT__)
// ROOT includes
#include <fstream>
#include <TString.h>
#include <TStopwatch.h>
#include <TMultiGraph.h>
#include <TSystem.h>
#include <TChain.h>
#include <TGraphErrors.h>
#include <TProof.h>
#include <TList.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TGrid.h>
#include <TEnv.h>
#include <TROOT.h>
#include <TAxis.h>
#include <THashList.h>
#include <TFileCollection.h>
#include <TGridCollection.h>
#include <TGridResult.h>

// STEER includes
#include "AliLog.h"
#include "AliCDBManager.h"
#include "AliAnalysisManager.h"
#include "AliESDInputHandler.h"
#include "AliTagAnalysis.h"
#include "AliRunTagCuts.h"
#include "AliLHCTagCuts.h"
#include "AliDetectorTagCuts.h"
#include "AliEventTagCuts.h"
#include "AliAnalysisDataContainer.h"
#include "AliAnalysisAlien.h"

// PHYSICS includes
#include "AliPhysicsSelectionTask.h"
#include "AliPhysicsSelection.h"
#include "AliMultSelectionTask.h"
#include "AliAnalysisTaskMuonResolution.h"

// MUON includes
#include "AliMpCDB.h"
#include "AliMpDetElement.h"
#include "AliMpDDLStore.h"
#include "AliMUONCalibParamND.h"
#include "AliMUON2DMap.h"
#include "AliMUONTrackerData.h"
#include "AliMUONPainterDataRegistry.h"
#include "AliMUONTrackerDataWrapper.h"

#include "AliMuonEventCuts.h"
#include "AliMuonTrackCuts.h"

#include "AddTaskMuonResolution.C"

#endif

enum {kLocal, kInteractif_xml, kInteractif_ESDList, kProof, kGrid, kTerminate};
Int_t nDE = 200;

Bool_t  Resume(Int_t mode, Int_t &firstStep, Double_t clusterResNB[10], Double_t clusterResB[10],
	       Double_t clusterResNBErr[10], Double_t clusterResBErr[10],
	       Bool_t shiftHalfCh, Double_t halfChShiftNB[20], Double_t halfChShiftB[20],
	       Double_t halfChShiftNBErr[20], Double_t halfChShiftBErr[20],
	       Bool_t shiftDE, Double_t deShiftNB[200], Double_t deShiftB[200],
	       TGraphErrors* clusterResXVsStep[10], TGraphErrors* clusterResYVsStep[10],
	       TGraphErrors* halfChShiftXVsStep[20], TGraphErrors* halfChShiftYVsStep[20]);
void    LoadAlirootOnProof(TString& aaf, TString rootVersion, TString aliphysicsVersion, Int_t iStep);
void CreateAlienHandler(TString runMode, TString& aliphysicsVersion, TString& runListName,
                        TString &dataDir, TString &dataPattern, TString &outDir, Int_t iStep,
                        TString runFormat, Int_t maxFilesPerJob, Int_t maxMergeFiles, Int_t maxMergeStages);
AliAnalysisTaskMuonResolution* CreateAnalysisTrain(Int_t mode, Int_t iStep, Bool_t selectPhysics, Bool_t selectTrigger,
                                                   Bool_t matchTrig, Bool_t applyAccCut, Bool_t applyPDCACut,
                                                   Double_t minMomentum, Double_t minPt, Bool_t isMC, Bool_t correctForSystematics,
						   Int_t extrapMode, Double_t clusterResNB[10], Double_t clusterResB[10],
						   Bool_t shiftHalfCh, Double_t halfChShiftNB[20], Double_t halfChShiftB[20],
						   Bool_t shiftDE, Double_t deShiftNB[200], Double_t deShiftB[200]);
Bool_t  GetChamberResolution(Int_t iStep, Double_t clusterResNB[10], Double_t clusterResB[10],
			     Double_t clusterResNBErr[10], Double_t clusterResBErr[10]);
Bool_t  AddHalfChShift(Int_t iStep, Double_t halfChShiftNB[20], Double_t halfChShiftB[20], Double_t halfChShiftNBErr[20], Double_t halfChShiftBErr[20]);
Bool_t  AddDEShift(Int_t iStep, Double_t deShiftNB[200], Double_t deShiftB[200]);
void    AddMCHViews(TString smode, TFile* file);
AliMUONTrackerData* ConvertGraph(TGraphErrors& g, const char* name);
Int_t   GetMode(TString smode, TString input);
TChain* CreateChainFromCollection(const char *xmlfile);
TChain* CreateChainFromFile(const char *rootfile);
TChain* CreateChainFromESDList(const char *esdList);
TChain* CreateChain(Int_t mode, TString input);

//______________________________________________________________________________
void MuonResolution(TString smode, TString inputFileName, Int_t nSteps,
                    TString rootVersion, TString aliphysicsVersion,
                    TString dataDir, TString dataPattern, TString runFormat, TString outDir,
                    Int_t maxFilesPerJob, Int_t maxMergeFiles, Int_t maxMergeStages,
		    Bool_t selectPhysics, Bool_t selectTrigger, Bool_t matchTrig, Bool_t applyAccCut, Bool_t applyPDCACut,
		    Double_t minMomentum, Double_t minPt, Bool_t isMC, Bool_t correctForSystematics, Int_t extrapMode,
		    Bool_t shiftHalfCh, Bool_t shiftDE, Int_t nevents)
{
  /// Compute the cluster resolution by studying cluster-track residual, deconvoluting from track resolution
  
  // timer start...
  TStopwatch* localTimer = new TStopwatch;
  
  // check parameters
  nSteps = TMath::Max(nSteps,1);
  if (extrapMode != 0 && extrapMode != 1) {
    Error("MuonResolution","incorrect extrapolation mode!");
    return;
  }
  
  // Check runing mode
  Int_t mode = GetMode(smode, inputFileName);
  if(mode < 0){
    Error("MuonResolution","Please provide either an ESD root file, a list of ESDs, a collection of ESDs or a dataset.");
    return;
  }
  
  // set starting chamber resolution (if -1 they will be loaded from recoParam in the task)
  Double_t clusterResNB[10] = {-1., -1., -1., -1., -1., -1., -1., -1., -1., -1.};
  Double_t clusterResB[10] = {-1., -1., -1., -1., -1., -1., -1., -1., -1., -1.};
  Double_t clusterResNBErr[10] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
  Double_t clusterResBErr[10] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
  Double_t halfChShiftNB[20] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
  Double_t halfChShiftB[20] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
  Double_t halfChShiftNBErr[20] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
  Double_t halfChShiftBErr[20] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
  Double_t deShiftNB[200];
  Double_t deShiftB[200];
  for (Int_t i=0; i<200; i++) {
    deShiftNB[i] = 0.;
    deShiftB[i] = 0.;
  }
  
  // output graphs
  TMultiGraph* mgClusterResXVsStep = new TMultiGraph("mgClusterResXVsStep","cluster X-resolution versus step;step;#sigma_{X} (cm)");
  TMultiGraph* mgClusterResYVsStep = new TMultiGraph("mgClusterResYVsStep","cluster Y-resolution versus step;step;#sigma_{Y} (cm)");
  TGraphErrors* clusterResXVsStep[10];
  TGraphErrors* clusterResYVsStep[10];
  for (Int_t i = 0; i < 10; i++) {
    clusterResXVsStep[i] = new TGraphErrors(nSteps+1);
    clusterResXVsStep[i]->SetName(Form("gResX_ch%d",i+1));
    clusterResXVsStep[i]->SetMarkerStyle(kFullDotMedium);
    clusterResXVsStep[i]->SetMarkerColor(i+1+i/9);
    mgClusterResXVsStep->Add(clusterResXVsStep[i],"lp");
    
    clusterResYVsStep[i] = new TGraphErrors(nSteps+1);
    clusterResYVsStep[i]->SetName(Form("gResY_ch%d",i+1));
    clusterResYVsStep[i]->SetMarkerStyle(kFullDotMedium);
    clusterResYVsStep[i]->SetMarkerColor(i+1+i/9);
    mgClusterResYVsStep->Add(clusterResYVsStep[i],"lp");
  }
  TMultiGraph* mgHalfChShiftXVsStep = new TMultiGraph("mgHalfChShiftXVsStep","half-chamber displacement in X direction versus step;step;#Delta_{X} (cm)");
  TMultiGraph* mgHalfChShiftYVsStep = new TMultiGraph("mgHalfChShiftYVsStep","half-chamber displacement in Y direction versus step;step;#Delta_{Y} (cm)");
  TGraphErrors* halfChShiftXVsStep[20];
  TGraphErrors* halfChShiftYVsStep[20];
  for (Int_t i = 0; i < 20; i++) {
    halfChShiftXVsStep[i] = new TGraphErrors(nSteps+1);
    halfChShiftXVsStep[i]->SetName(Form("gShiftX_hch%d",i+1));
    halfChShiftXVsStep[i]->SetMarkerStyle(kFullDotMedium);
    halfChShiftXVsStep[i]->SetMarkerColor(i+1+i/9+i/18);
    mgHalfChShiftXVsStep->Add(halfChShiftXVsStep[i],"lp");
    halfChShiftXVsStep[i]->SetPoint(0, 0, halfChShiftNB[i]);
    halfChShiftXVsStep[i]->SetPointError(0, 0., halfChShiftNBErr[i]);
    
    halfChShiftYVsStep[i] = new TGraphErrors(nSteps+1);
    halfChShiftYVsStep[i]->SetName(Form("gShiftY_hch%d",i+1));
    halfChShiftYVsStep[i]->SetMarkerStyle(kFullDotMedium);
    halfChShiftYVsStep[i]->SetMarkerColor(i+1+i/9+i/18);
    mgHalfChShiftYVsStep->Add(halfChShiftYVsStep[i],"lp");
    halfChShiftYVsStep[i]->SetPoint(0, 0, halfChShiftB[i]);
    halfChShiftYVsStep[i]->SetPointError(0, 0., halfChShiftBErr[i]);
  }
  
  // check for old output files
  Int_t firstStep = 0;
  char remove = '\0';
  if (!gSystem->Exec("ls chamberResolution_step*[0-9].root")) {
    cout<<"Above files already exist in the current directory. [d=delete, r=resume, e=exit] "<<flush;
    while (remove != 'd' && remove != 'r' && remove != 'e') cin>>remove;
    if (remove == 'y') gSystem->Exec("rm -f chamberResolution_step*[0-9].root");
    else if (remove == 'r' && !Resume(mode, firstStep, clusterResNB, clusterResB, clusterResNBErr, clusterResBErr,
				      shiftHalfCh, halfChShiftNB, halfChShiftB, halfChShiftNBErr, halfChShiftBErr,
				      shiftDE, deShiftNB, deShiftB, clusterResXVsStep, clusterResYVsStep,
				      halfChShiftXVsStep, halfChShiftYVsStep)) return;
    else if (remove == 'e') return;
  }
  
  // Create input object
  TObject* inputObj = 0x0;
  if (mode != kGrid && mode != kTerminate) {
    if (mode == kProof) {
      if (inputFileName.EndsWith(".root")) {
        TFile *inFile = TFile::Open(inputFileName.Data(),"READ");
        if (inFile && inFile->IsOpen()) {
          inputObj = dynamic_cast<TFileCollection*>(inFile->FindObjectAny("dataset"));
          inFile->Close();
        }
      } else inputObj = new TObjString(inputFileName);
    } else inputObj = CreateChain(mode, inputFileName);
    if (!inputObj) return;
  }

  // loop over step
  for (Int_t iStep = firstStep; iStep < nSteps; iStep++) {
    cout<<"step "<<iStep+1<<"/"<<nSteps<<endl;
    
    // create the analysis train
    AliAnalysisTaskMuonResolution *muonResolution = CreateAnalysisTrain(mode, iStep, selectPhysics, selectTrigger, matchTrig,
				                        applyAccCut, applyPDCACut, minMomentum, minPt, isMC,
                                                        correctForSystematics, extrapMode, clusterResNB, clusterResB,
                                                        shiftHalfCh, halfChShiftNB, halfChShiftB, shiftDE, deShiftNB, deShiftB);
    if (!muonResolution) return;
    
    // prepare proof or grid environment
    if (mode == kProof) LoadAlirootOnProof(smode, rootVersion, aliphysicsVersion, iStep-firstStep);
    else if (mode == kGrid || mode == kTerminate) {
      if (!gGrid) TGrid::Connect("alien://");
      TString resultDir = Form("%s/%s/step%d", gGrid->GetHomeDirectory(), outDir.Data(), iStep);
      if ((smode == "submit" || smode == "full") && AliAnalysisAlien::DirectoryExists(resultDir.Data())) {
        cout << endl << "Output directory alien://" << resultDir <<" already exist." << endl;
        cout << "Do you want to continue? [Y/n] " << endl;
        TString reply = "";
        reply.Gets(stdin,kTRUE);
        reply.ToLower();
        if (reply.BeginsWith("n")) return;
      }
      CreateAlienHandler(smode, aliphysicsVersion, inputFileName, dataDir, dataPattern, outDir, iStep,
                         runFormat, maxFilesPerJob, maxMergeFiles, maxMergeStages);
    }
    
    // start analysis
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (mgr->InitAnalysis()) {
      mgr->PrintStatus();
      if (mode == kGrid) {
        mgr->StartAnalysis("grid");
        if (smode != "terminate") return; // stop here if we don't have yet the result of this step
      } else if (mode == kTerminate) mgr->StartAnalysis("grid terminate");
      else if (mode == kProof) {
        if (inputObj->IsA() == TFileCollection::Class())
          mgr->StartAnalysis("proof", static_cast<TFileCollection*>(inputObj), nevents);
        else mgr->StartAnalysis("proof", static_cast<TObjString*>(inputObj)->GetName(), nevents);
      } else mgr->StartAnalysis("local", static_cast<TChain*>(inputObj), nevents);
    }
    
    // save the summary canvases and mchview display
    if (muonResolution->GetCanvases()) {
      TFile* outFile = TFile::Open(Form("chamberResolution_step%d.root", iStep),"UPDATE");
      if (outFile && outFile->IsOpen()) {
	outFile->cd();
	muonResolution->GetCanvases()->Write();
	AddMCHViews(smode, outFile);
	outFile->Close();
	delete outFile;
      }
    }
    
    // fill graphs with starting resolutions from the task at very first step
    if (iStep == 0) {
      muonResolution->GetStartingResolution(clusterResNB, clusterResB);
      for (Int_t i = 0; i < 10; i++) {
	clusterResXVsStep[i]->SetPoint(0, 0, clusterResNB[i]);
	clusterResXVsStep[i]->SetPointError(0, 0., clusterResNBErr[i]);
	clusterResYVsStep[i]->SetPoint(0, 0, clusterResB[i]);
	clusterResYVsStep[i]->SetPointError(0, 0., clusterResBErr[i]);
      }
    }
    
    // read the chamber resolution from the output file
    if (!GetChamberResolution(iStep, clusterResNB, clusterResB, clusterResNBErr, clusterResBErr)) return;
    
    // fill graphs with computed resolutions
    for (Int_t i = 0; i < 10; i++) {
      clusterResXVsStep[i]->SetPoint(iStep+1, iStep+1, clusterResNB[i]);
      clusterResXVsStep[i]->SetPointError(iStep+1, 0., clusterResNBErr[i]);
      clusterResYVsStep[i]->SetPoint(iStep+1, iStep+1, clusterResB[i]);
      clusterResYVsStep[i]->SetPointError(iStep+1, 0., clusterResBErr[i]);
    }
    
    // get the half-chamber displacements currently used and add the new measurements from the output file
    muonResolution->GetHalfChShift(halfChShiftNB, halfChShiftB);
    if (!AddHalfChShift(iStep, halfChShiftNB, halfChShiftB, halfChShiftNBErr, halfChShiftBErr)) return;
    
    // fill graphs with computed displacements
    for (Int_t i = 0; i < 20; i++) {
      halfChShiftXVsStep[i]->SetPoint(iStep+1, iStep+1, halfChShiftNB[i]);
      halfChShiftXVsStep[i]->SetPointError(iStep+1, 0., halfChShiftNBErr[i]);
      halfChShiftYVsStep[i]->SetPoint(iStep+1, iStep+1, halfChShiftB[i]);
      halfChShiftYVsStep[i]->SetPointError(iStep+1, 0., halfChShiftBErr[i]);
    }
    
    // get the DE displacements currently used and add the new measurements from the output file
    muonResolution->GetDEShift(deShiftNB, deShiftB);
    if (!AddDEShift(iStep, deShiftNB, deShiftB)) return;
    
    // clean memory
    mgr->SetAnalysisType(AliAnalysisManager::kLocalAnalysis); // to make sure that all output containers are deleted
    delete mgr;
    TObject::SetObjectStat(kFALSE);
    
    // in grid mode, stop here unless it is the final step
    if ((mode == kGrid || mode == kTerminate) && iStep != nSteps-1) return;
    
  }
  
  // copy final results in results.root file
  gSystem->Exec(Form("cp chamberResolution_step%d.root results.root", nSteps-1));
  
  // display convergence of cluster resolution
  TCanvas* convergence1 = new TCanvas("convergenceRes","convergence of cluster resolution");
  convergence1->Divide(1,2);
  convergence1->cd(1);
  mgClusterResXVsStep->Draw("ap");
  convergence1->cd(2);
  mgClusterResYVsStep->Draw("ap");
  
  // display convergence of half-chamber displacements
  TCanvas* convergence2 = new TCanvas("convergenceShift","convergence of half-chamber displacements");
  convergence2->Divide(1,2);
  convergence2->cd(1);
  mgHalfChShiftXVsStep->Draw("ap");
  convergence2->cd(2);
  mgHalfChShiftYVsStep->Draw("ap");
  
  // save convergence plots
  TFile* outFile = TFile::Open("results.root","UPDATE");
  if (!outFile || !outFile->IsOpen()) return;
  outFile->cd();
  mgClusterResXVsStep->Write();
  mgClusterResYVsStep->Write();
  convergence1->Write();
  mgHalfChShiftXVsStep->Write();
  mgHalfChShiftYVsStep->Write();
  convergence2->Write();
  outFile->Close();
  delete outFile;
  
  // print final half-chamber displacements
  printf("\nhalf-chamber total displacements:\n");
  printf(" - non-bending:");
  for (Int_t i = 0; i < 20; i++) printf((i==0)?" %6.4f":", %6.4f", halfChShiftNB[i]);
  printf("\n -     bending:");
  for (Int_t i = 0; i < 20; i++) printf((i==0)?" %6.4f":", %6.4f", halfChShiftB[i]);
  printf("\n\n");
  
  // print final DE displacements
  printf("\nDE total displacements:\n");
  printf(" - non-bending:");
  for (Int_t i = 0; i < nDE; i++) printf((i==0)?" %6.4f":", %6.4f", deShiftNB[i]);
  printf("\n -     bending:");
  for (Int_t i = 0; i < nDE; i++) printf((i==0)?" %6.4f":", %6.4f", deShiftB[i]);
  printf("\n\n");
  
  // ...timer stop
  localTimer->Stop();
  localTimer->Print();
  delete localTimer;
  
}

//______________________________________________________________________________
Bool_t Resume(Int_t mode, Int_t &firstStep, Double_t clusterResNB[10], Double_t clusterResB[10],
	      Double_t clusterResNBErr[10], Double_t clusterResBErr[10],
	      Bool_t shiftHalfCh, Double_t halfChShiftNB[20], Double_t halfChShiftB[20],
	      Double_t halfChShiftNBErr[20], Double_t halfChShiftBErr[20],
	      Bool_t shiftDE, Double_t deShiftNB[200], Double_t deShiftB[200],
	      TGraphErrors* clusterResXVsStep[10], TGraphErrors* clusterResYVsStep[10],
	      TGraphErrors* halfChShiftXVsStep[20], TGraphErrors* halfChShiftYVsStep[20])
{
  /// resume analysis from desired step
  /// do not remove any step in terminate-only mode
  
  while (kTRUE) {
    
    // Get the step to restart from
    cout<<"From which step (included) you want to resume? [#, e=exit] "<<flush;
    TString step = "";
    do {step.Gets(stdin,kTRUE);} while (!step.IsDigit() && step != "e");
    if (step == "e") return kFALSE;
    firstStep = step.Atoi();
    
    // restart from scratch if requested
    if (firstStep == 0 && mode != kTerminate) {
      gSystem->Exec("rm -f chamberResolution_step*[0-9].root");
      return kTRUE;
    }
    
    // look for results from the previous step
    if (firstStep != 0 && gSystem->AccessPathName(Form("chamberResolution_step%d.root", firstStep-1))) {
      cout<<"No result found from the previous step ("<<firstStep-1<<"). Unable to resume from step "<<firstStep<<endl;
      continue;
    }
    
    // fill graph with starting resolutions
    for (Int_t i = 0; i < 10; i++) {
      clusterResXVsStep[i]->SetPoint(0, 0, clusterResNB[i]);
      clusterResXVsStep[i]->SetPointError(0, 0., clusterResNBErr[i]);
      clusterResYVsStep[i]->SetPoint(0, 0, clusterResB[i]);
      clusterResYVsStep[i]->SetPointError(0, 0., clusterResBErr[i]);
    }
    
    // loop over previous steps
    Bool_t missingInfo = kFALSE;
    for (Int_t iStep = 0; iStep < firstStep; iStep++) {
      
      // read the chamber resolution from the output file
      if (!GetChamberResolution(iStep, clusterResNB, clusterResB, clusterResNBErr, clusterResBErr) && iStep == firstStep-1) {
	missingInfo = kTRUE;
	break;
      }
      
      // fill graphs with computed resolutions
      for (Int_t i = 0; i < 10; i++) {
	clusterResXVsStep[i]->SetPoint(iStep+1, iStep+1, clusterResNB[i]);
	clusterResXVsStep[i]->SetPointError(iStep+1, 0., clusterResNBErr[i]);
	clusterResYVsStep[i]->SetPoint(iStep+1, iStep+1, clusterResB[i]);
	clusterResYVsStep[i]->SetPointError(iStep+1, 0., clusterResBErr[i]);
      }
      
      // reset the half-chamber displacements if not used and add the new measurements from the output file
      if (!shiftHalfCh) for (Int_t i=0; i<20; i++) {
	halfChShiftNB[i] = 0.; halfChShiftB[i] = 0.;
	halfChShiftNBErr[i] = 0.; halfChShiftBErr[i] = 0.;
      }
      if (!AddHalfChShift(iStep, halfChShiftNB, halfChShiftB, halfChShiftNBErr, halfChShiftBErr) && shiftHalfCh) {
	missingInfo = kTRUE;
	break;
      }
      
      // fill graphs with computed displacements
      for (Int_t i = 0; i < 20; i++) {
	halfChShiftXVsStep[i]->SetPoint(iStep+1, iStep+1, halfChShiftNB[i]);
	halfChShiftXVsStep[i]->SetPointError(iStep+1, 0., halfChShiftNBErr[i]);
	halfChShiftYVsStep[i]->SetPoint(iStep+1, iStep+1, halfChShiftB[i]);
	halfChShiftYVsStep[i]->SetPointError(iStep+1, 0., halfChShiftBErr[i]);
      }
      
      // add the new measurements of DE displacements from the output file if in use
      if (shiftDE && !AddDEShift(iStep, deShiftNB, deShiftB)) {
	missingInfo = kTRUE;
	break;
      }
      
    }
    
    // check if missing important results from previous steps
    if (missingInfo) continue;
    
    // keep previous steps and remove the others
    if (mode != kTerminate) {
      gSystem->Exec("mkdir __TMP__");
      for (Int_t iStep = 0; iStep < firstStep; iStep++)
        if (!gSystem->AccessPathName(Form("chamberResolution_step%d.root", iStep)))
          gSystem->Exec(Form("mv chamberResolution_step%d.root __TMP__", iStep));
      gSystem->Exec("rm -f chamberResolution_step*[0-9].root");
      gSystem->Exec("mv __TMP__/chamberResolution_step*[0-9].root .");
      gSystem->Exec("rm -rf __TMP__");
    }
    
    return kTRUE;
  }
  
}

//______________________________________________________________________________
void LoadAlirootOnProof(TString& aaf, TString rootVersion, TString aliphysicsVersion, Int_t iStep)
{
  /// Load aliroot packages and set environment on Proof
  
  // set general environment and close previous session
  if (iStep == 0) gEnv->SetValue("XSec.GSI.DelegProxy","2");
  else gProof->Close("s");
  
  // connect
  if (aaf == "saf3") TProof::Open("pod://");
  else {
    TString location = (aaf == "caf") ? "alice-caf.cern.ch" : "nansafmaster2.in2p3.fr"; //"localhost:1093"
    TString nWorkers = (aaf == "caf") ? "workers=80" : ""; //"workers=3x"
    TString user = (gSystem->Getenv("alien_API_USER") == NULL) ? "" : Form("%s@",gSystem->Getenv("alien_API_USER"));
    TProof::Mgr(Form("%s%s",user.Data(), location.Data()))->SetROOTVersion(Form("VO_ALICE@ROOT::%s",rootVersion.Data()));
    TProof::Open(Form("%s%s/?N",user.Data(), location.Data()), nWorkers.Data());
  }
  if (!gProof) return;
  
  // set environment and load libraries on workers
  TList* list = new TList();
  list->Add(new TNamed("ALIROOT_MODE", "base"));
  list->Add(new TNamed("ALIROOT_ENABLE_ALIEN", "1"));
  if (aaf == "saf3") {
    TString home = gSystem->Getenv("HOME");
    gProof->UploadPackage(Form("%s/AliceVaf.par", home.Data()));
    gProof->EnablePackage(Form("%s/AliceVaf.par", home.Data()), list, (iStep!=0));
  } else gProof->EnablePackage(Form("VO_ALICE@AliPhysics::%s",aliphysicsVersion.Data()), list, kTRUE);
//  gProof->UploadPackage("$ALICE_PHYSICS/PARfiles/PWGPPMUONdep.par");
//  gProof->EnablePackage("$ALICE_PHYSICS/PARfiles/PWGPPMUONdep.par", kTRUE);
//  gProof->UploadPackage("PWGPPMUONdep.par");
//  gProof->EnablePackage("PWGPPMUONdep.par", kTRUE);
  
}

//______________________________________________________________________________
void CreateAlienHandler(TString runMode, TString& aliphysicsVersion, TString& runListName,
                        TString &dataDir, TString &dataPattern, TString &outDir, Int_t iStep,
                        TString runFormat, Int_t maxFilesPerJob, Int_t maxMergeFiles, Int_t maxMergeStages)
{
  // Configure the alien plugin
  AliAnalysisAlien *plugin = new AliAnalysisAlien();
  
  // If the run mode is merge, run in mode terminate to merge via jdl
  // If the run mode is terminate, disable the mergin via jdl
  Bool_t merge = kTRUE;
  if (runMode.Contains("terminate")) merge = kFALSE;
  else if (runMode == "merge") runMode = "terminate";
  
  // Set the run mode (can be "full", "test", "offline", "submit" or "terminate")
  plugin->SetRunMode(runMode.Data());
  
  // Set the number of input files in test mode
  plugin->SetNtestFiles(2);
  
  // Set versions of used packages
  plugin->SetAPIVersion("V1.1x");
  if (!aliphysicsVersion.IsNull()) plugin->SetAliPhysicsVersion(aliphysicsVersion.Data());
  
  // Declare input data to be processed
  plugin->SetGridDataDir(dataDir.Data());
  plugin->SetDataPattern(dataPattern.Data());
  ifstream inFile(runListName.Data());
  TString currRun;
  if (inFile.is_open())
  {
    while (! inFile.eof() )
    {
      currRun.ReadLine(inFile,kTRUE); // Read line
      if(currRun.IsNull()) continue;
      plugin->AddRunNumber(Form(runFormat.Data(), currRun.Atoi()));
    }
  }
  inFile.close();
  
  // Define alien work directory where all files will be copied. Relative to alien $HOME
  plugin->SetGridWorkingDir(outDir.Data());
  
  // Declare alien output directory. Relative to working directory
  plugin->SetGridOutputDir(Form("step%d",iStep));
  
  // Set the ouput directory of each masterjob to the run number
  plugin->SetOutputToRunNo();
  
  // Declare all libraries (other than the default ones for the framework)
  plugin->SetAdditionalRootLibs("libGui.so libProofPlayer.so libXMLParser.so");
  
  // Optionally add include paths
  plugin->AddIncludePath("-I$ALICE_ROOT/include");
  plugin->AddIncludePath("-I$ALICE_PHYSICS/include");
  plugin->AddIncludePath("-I.");
  
  // Optionally add packages
//  plugin->EnablePackage("PWGPPMUONdep.par");
  
  // Optionally set a name for the generated analysis macro (default MyAnalysis.C)
  plugin->SetAnalysisMacro("Resolution.C");
  plugin->SetExecutable("Resolution.sh");
  
  // Optionally set time to live (default 30000 sec)
  plugin->SetTTL(30000);
  
  // Optionally set input format (default xml-single)
  plugin->SetInputFormat("xml-single");
  
  // Optionally modify the name of the generated JDL (default analysis.jdl)
  plugin->SetJDLName("Resolution.jdl");
  
  // Optionally modify job price (default 1)
  plugin->SetPrice(1);
  
  // Optionally modify split mode (default 'se')
  plugin->SetSplitMode("se");
  
  // Optionally modify the maximum number of files per job
  plugin->SetSplitMaxInputFileNumber(maxFilesPerJob);
  
  // Merge via JDL
  if (merge) plugin->SetMergeViaJDL(kTRUE);
  
  // Optionally set the maximum number of files merged together in one stage
  plugin->SetMaxMergeFiles(maxMergeFiles);
  
  // Optionally set the maximum number of merging stages
  plugin->SetMaxMergeStages(maxMergeStages);
  
  // Exclude given output file(s) from registration/merging
  plugin->SetRegisterExcludes("AnalysisResults.root EventStat_temp.root");
  
  // Save the log files
  plugin->SetKeepLogs();
  
  AliAnalysisManager::GetAnalysisManager()->SetGridHandler(plugin);
}

//______________________________________________________________________________
AliAnalysisTaskMuonResolution* CreateAnalysisTrain(Int_t mode, Int_t iStep, Bool_t selectPhysics, Bool_t selectTrigger,
                                                   Bool_t matchTrig, Bool_t applyAccCut, Bool_t applyPDCACut,
                                                   Double_t minMomentum, Double_t minPt, Bool_t isMC, Bool_t correctForSystematics,
                                                   Int_t extrapMode, Double_t clusterResNB[10], Double_t clusterResB[10],
                                                   Bool_t shiftHalfCh, Double_t halfChShiftNB[20], Double_t halfChShiftB[20],
                                                   Bool_t shiftDE, Double_t deShiftNB[200], Double_t deShiftB[200])
{
  /// create the analysis train and configure it
  
  // Create the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("MuonResolutionAnalysis");
  //mgr->SetNSysInfo(100);
  //mgr->SetDebugLevel(3);
  
  // ESD input handler
  AliESDInputHandler* esdH = new AliESDInputHandler();
  esdH->SetReadFriends(kFALSE);
  esdH->SetInactiveBranches("*");
  esdH->SetActiveBranches("MuonTracks MuonClusters MuonPads AliESDRun. AliESDHeader. AliMultiplicity. AliESDFMD. AliESDVZERO. SPDVertex. PrimaryVertex. AliESDZDC.");
  mgr->SetInputEventHandler(esdH);
  
  // event selection
  UInt_t eventSelectionMask = 0;
  if (selectPhysics) {
    gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
    AliPhysicsSelectionTask* physicsSelection = reinterpret_cast<AliPhysicsSelectionTask*>(gROOT->ProcessLineSync(TString::Format("AddTaskPhysicsSelection(%d)", isMC)));
    if (!physicsSelection) {
      Error("CreateAnalysisTrain","AliPhysicsSelectionTask not created!");
      return 0x0;
    }
    //if (!isMC) physicsSelection->GetPhysicsSelection()->SetUseBXNumbers(kFALSE); // Needed to merge runs with different running scheme
    eventSelectionMask |= AliMuonEventCuts::kPhysicsSelected;
  }
  if (selectTrigger) eventSelectionMask |= AliMuonEventCuts::kSelectedTrig;
  
  // multiplicity/centrality selection
  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");
  AliMultSelectionTask *mult = reinterpret_cast<AliMultSelectionTask*>(gROOT->ProcessLineSync("AddTaskMultSelection()"));
  if (!mult) {
    Error("CreateAnalysisTrain","AliMultSelectionTask not created!");
    return 0x0;
  }
//  mult->SetAlternateOADBforEstimators("LHC15o");
  
  // Muon Resolution analysis
  TString outputFileName = Form("chamberResolution_step%d.root", iStep);
  AliAnalysisManager::SetCommonFileName(outputFileName.Data());
  AliAnalysisTaskMuonResolution *muonResolution = AddTaskMuonResolution(minMomentum, minPt, correctForSystematics, extrapMode);
  if (!muonResolution) {
    Error("CreateAnalysisTrain","AliAnalysisTaskMuonResolution not created!");
    return 0x0;
  }
  /*if (mode == kLocal) muonResolution->SetDefaultStorage("local://$ALIROOT_OCDB_ROOT/OCDB");
  else muonResolution->SetDefaultStorage("raw://");*/
  muonResolution->SetDefaultStorage("raw://");
//  muonResolution->SetDefaultStorage("local:///cvmfs/alice-ocdb.cern.ch/calibration/data/2017/OCDB");
  if (mode != kProof) muonResolution->ShowProgressBar();
  muonResolution->PrintClusterRes(kTRUE, kTRUE);
  muonResolution->SetStartingResolution(clusterResNB, clusterResB);
  muonResolution->RemoveMonoCathodClusters(kTRUE, kFALSE);
//  muonResolution->FitResiduals(kFALSE);
//  muonResolution->ImproveTracks(kTRUE);
//  muonResolution->ReAlign("", 1, -1, "alien://folder=/alice/cern.ch/user/h/hupereir/CDB/LHC17g_realign");
//  muonResolution->ReAlign("", 1, 0, "");
//  muonResolution->SetAlignStorage("", 1);
//  muonResolution->SetMuonSign(-1);
//  muonResolution->UseMCLabel();
  
  if (shiftHalfCh) {
    muonResolution->SetHalfChShift(halfChShiftNB, halfChShiftB);
    muonResolution->ShiftHalfCh();
    muonResolution->PrintHalfChShift();
  }
  if (shiftDE) {
    muonResolution->SetDEShift(deShiftNB, deShiftB);
    muonResolution->ShiftDE();
    muonResolution->PrintDEShift();
  }
  
  if (eventSelectionMask != 0) {
    AliMuonEventCuts eventCuts("muEventCuts", "muEventCuts");
    eventCuts.SetFilterMask(eventSelectionMask);
    if (selectPhysics) eventCuts.SetPhysicsSelectionMask(AliVEvent::kAny);
//    if (selectPhysics) eventCuts.SetPhysicsSelectionMask(AliVEvent::kINT7inMUON);
//    if (selectPhysics) eventCuts.SetPhysicsSelectionMask(AliVEvent::kMuonUnlikeLowPt7);
    if (selectTrigger) eventCuts.SetTrigClassPatterns(eventCuts.GetDefaultTrigClassPatterns());
    muonResolution->SetMuonEventCuts(eventCuts);
  }
  
  UInt_t trackSelectionMask = 0;
  if (matchTrig) trackSelectionMask |= AliMuonTrackCuts::kMuMatchLpt;
  if (applyAccCut) trackSelectionMask |= AliMuonTrackCuts::kMuEta | AliMuonTrackCuts::kMuThetaAbs;
  if (applyPDCACut) trackSelectionMask |= AliMuonTrackCuts::kMuPdca;
  if (trackSelectionMask != 0) {
    AliMuonTrackCuts trackCuts("muTrackCuts", "muTrackCuts");
    trackCuts.SetAllowDefaultParams();
    if (isMC) trackCuts.SetIsMC();
    trackCuts.SetFilterMask(trackSelectionMask);
    muonResolution->SetMuonTrackCuts(trackCuts);
  }
  
  return muonResolution;
  
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
  for (Int_t i = 0; i < 10; i++) {
    gCombinedResidualXPerChSigma->GetPoint(i, dummy, clusterResNB[i]);
    gCombinedResidualYPerChSigma->GetPoint(i, dummy, clusterResB[i]);
    clusterResNBErr[i] = gCombinedResidualXPerChSigma->GetErrorY(i);
    clusterResBErr[i] = gCombinedResidualYPerChSigma->GetErrorY(i);
  }
  
  outFile->Close();
  //delete outFile;
  
  return kTRUE;
}

//______________________________________________________________________________
Bool_t AddHalfChShift(Int_t iStep, Double_t halfChShiftNB[20], Double_t halfChShiftB[20], Double_t halfChShiftNBErr[20], Double_t halfChShiftBErr[20])
{
  /// read the chamber resolution from the output file
  
  TFile* outFile = TFile::Open(Form("chamberResolution_step%d.root", iStep),"READ");
  
  if (!outFile || !outFile->IsOpen()) {
    Error("AddHalfChShift","output file does not exist!");
    return kFALSE;
  }
  
  TObjArray* summary = static_cast<TObjArray*>(outFile->FindObjectAny("ChamberRes"));
  TGraphErrors* gResidualXPerHalfChMean = (summary) ? static_cast<TGraphErrors*>(summary->FindObject("gResidualXPerHalfChMean_ClusterIn")) : 0x0;
  TGraphErrors* gResidualYPerHalfChMean = (summary) ? static_cast<TGraphErrors*>(summary->FindObject("gResidualYPerHalfChMean_ClusterIn")) : 0x0;
  
  if (!gResidualXPerHalfChMean || !gResidualYPerHalfChMean) {
    Error("AddHalfChShift","half-chamber shift graphs do not exist!");
    return kFALSE;
  }
  
  Double_t dummy, dx, dy;
  for (Int_t i = 0; i < 20; i++) {
    gResidualXPerHalfChMean->GetPoint(i, dummy, dx);
    halfChShiftNB[i] += dx;
    halfChShiftNBErr[i] = gResidualXPerHalfChMean->GetErrorY(i);
    gResidualYPerHalfChMean->GetPoint(i, dummy, dy);
    halfChShiftB[i] += dy;
    halfChShiftBErr[i] = gResidualYPerHalfChMean->GetErrorY(i);
  }
  
  outFile->Close();
  //delete outFile;
  
  return kTRUE;
}

//______________________________________________________________________________
Bool_t AddDEShift(Int_t iStep, Double_t deShiftNB[200], Double_t deShiftB[200])
{
  /// read the chamber resolution from the output file
  
  TFile* outFile = TFile::Open(Form("chamberResolution_step%d.root", iStep),"READ");
  
  if (!outFile || !outFile->IsOpen()) {
    Error("AddDEShift","output file does not exist!");
    return kFALSE;
  }
  
  TObjArray* summary = static_cast<TObjArray*>(outFile->FindObjectAny("ChamberRes"));
  TGraphErrors* gResidualXPerDEMean = (summary) ? static_cast<TGraphErrors*>(summary->FindObject("gResidualXPerDEMean_ClusterIn")) : 0x0;
  TGraphErrors* gResidualYPerDEMean = (summary) ? static_cast<TGraphErrors*>(summary->FindObject("gResidualYPerDEMean_ClusterIn")) : 0x0;
  
  if (!gResidualXPerDEMean || !gResidualYPerDEMean) {
    Error("AddDEShift","DE shift graphs do not exist!");
    return kFALSE;
  }
  
  Double_t dummy, dx, dy;
  nDE = gResidualXPerDEMean->GetN();
  for (Int_t i = 0; i < nDE; i++) {
    gResidualXPerDEMean->GetPoint(i, dummy, dx);
    deShiftNB[i] += dx;
    gResidualYPerDEMean->GetPoint(i, dummy, dy);
    deShiftB[i] += dy;
  }
  
  outFile->Close();
  //delete outFile;
  
  return kTRUE;
}

//______________________________________________________________________________
void AddMCHViews(TString smode, TFile* file)
{
  /// Get from the file the graphs containing data per DE, convert them into mchview objects and save them
  
  if (  ! AliMpDDLStore::Instance(false) )
  {
    Warning("AddMCHViews","mapping was not loaded. Loading it from OCDB");
    if (smode == "saf3") AliCDBManager::Instance()->SetDefaultStorage("raw://");
    else AliCDBManager::Instance()->SetDefaultStorage("local://$ALIROOT_OCDB_ROOT/OCDB");
    AliCDBManager::Instance()->SetRun(999999999);
  }
  
  AliMpCDB::LoadAll();
  
  TObjArray* summary = static_cast<TObjArray*>(file->FindObjectAny("ChamberRes"));
  if (!summary) {
    Error("AddMCHViews","resolution graphs do not exist!");
    return;
  }
  
  TGraphErrors* g = 0x0;
  g = static_cast<TGraphErrors*>(summary->FindObject("gCombinedResidualXPerDESigma"));
  if (g) {
    file->cd();
    AliMUONTrackerData* data = ConvertGraph(*g, "resoX");
    data->Write();
    delete data;
  }
  
  g = static_cast<TGraphErrors*>(summary->FindObject("gCombinedResidualYPerDESigma"));
  if (g) {
    file->cd();
    AliMUONTrackerData* data = ConvertGraph(*g, "resoY");
    data->Write();
    delete data;
  }
  
  g = static_cast<TGraphErrors*>(summary->FindObject("gResidualXPerDEMean_ClusterOut"));
  if (g) {
    file->cd();
    AliMUONTrackerData* data = ConvertGraph(*g, "shiftX");
    data->Write();
    delete data;
  }
  
  g = static_cast<TGraphErrors*>(summary->FindObject("gResidualYPerDEMean_ClusterOut"));
  if (g) {
    file->cd();
    AliMUONTrackerData* data = ConvertGraph(*g, "shiftY");
    data->Write();
    delete data;
  }
}

//______________________________________________________________________________
AliMUONTrackerData* ConvertGraph(TGraphErrors& g, const char* name)
{
  /// Convert graph containing data per DE into mchview object
  
  AliMUON2DMap deValues(kFALSE);
  
  for ( Int_t i = 0 ; i < g.GetN(); ++i ) 
  {
    double y = g.GetY()[i];
    double ey = g.GetEY()[i];
    int detElemId;
    sscanf(g.GetXaxis()->GetBinLabel(i+1),"%d",&detElemId);
    
    AliMpDetElement* de = AliMpDDLStore::Instance()->GetDetElement(detElemId);
    
    AliMUONVCalibParam* param = new AliMUONCalibParamND(5, 1, detElemId, 0);
    
    Double_t sumn = 1000.0;
    Double_t sumw = sumn*y;
    Double_t sumw2 = (sumn-1)*ey*ey+sumw*sumw/sumn;
    
    param->SetValueAsDouble(0,0,sumw);
    param->SetValueAsDouble(0,1,sumw2);
    param->SetValueAsDouble(0,2,sumn);
    param->SetValueAsDouble(0,3,de->NofChannels());
    param->SetValueAsDouble(0,4,1);
    
    deValues.Add(param);
  }
  
  AliMUONTrackerData* data = new AliMUONTrackerData(name,name,deValues,1);
  data->SetDimensionName(0,name);
  
  return data;
}

//______________________________________________________________________________
Int_t GetMode(TString smode, TString input)
{
  if (smode == "local") {
    if ( input.EndsWith(".xml") ) return kInteractif_xml;
    else if ( input.EndsWith(".txt") ) return kInteractif_ESDList;
    else if ( input.EndsWith(".root") ) return kLocal;    
  } else if (smode == "caf" || smode == "saf" || smode == "saf3") return kProof;
  else if ((smode == "test" || smode == "offline" || smode == "submit" || smode == "full" ||
            smode == "merge" || smode == "terminate") && input.EndsWith(".txt")) return kGrid;
  else if (smode == "terminateonly") return kTerminate;
  return -1;
}

//______________________________________________________________________________
TChain* CreateChainFromCollection(const char *xmlfile)
{
  // Create a chain from the collection of tags.
  if (!TGrid::Connect("alien://")) return NULL;
  
  TGridCollection* coll = gGrid->OpenCollection(xmlfile);
  if (!coll) {
    Error("CreateChainFromCollection", "Cannot create the AliEn collection");
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
  if (!gSystem->AccessPathName("ConfigureCuts.C"))
    gROOT->ProcessLine(Form(".x ConfigureCuts.C((AliRunTagCuts*)%p, (AliLHCTagCuts*)%p, (AliDetectorTagCuts*)%p,"
			    " (AliEventTagCuts*)%p)", runCuts, lhcCuts, detCuts, evCuts));
  
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

