#include <TFile.h>
#include <TList.h>
#include <TTree.h>
#include <TChain.h>
#include <TClonesArray.h>
#include <TObjString.h>
#include <TString.h>
#include <TROOT.h>

#include "AliHighPtDeDxData.h"
#include "AliHighPtDeDxMc.h"
#include "AliHighPtDeDxCalib.h"

#include <AliXRDPROOFtoolkit.h>

#include "DebugClasses.C"
#include "my_tools.C"
#include "my_functions.C"

#include <iostream>
#include <fstream>
#include <string>

using namespace std;

/*

  In this version the start step for calibrations is 2. We do not correct for
  the Ncl dependence as this is still uncertain of this should be done or not.

  It is good to make a small test first because the MIP peak might have to be
  adjusted. With pass0 the default 40-60 window should be good, but

  To run calibrations:
  ====================

  Use AliRoot because of AliXRDPROOFtoolkit:
   gSystem->AddIncludePath("-I$ALICE_ROOT/TPC")
   gSystem->AddIncludePath("-I../lib")
   gSystem->AddIncludePath("-I../grid")
   gSystem->AddIncludePath("-I../macros")
   gROOT->SetMacroPath(".:../macros:../grid:../lib/")
  .L libMyDeDxAnalysis.so 
  .L my_functions.C+
  .L my_tools.C+
  .L DebugClasses.C+
  .L run_code.C+

  // Step 1: create calibrations

  CreateCalib("7tev_b.dat", kFALSE, "7tev_b.root", 0, 2)  

  // This function is not used anymore
  // CreateCalibV0("7tev_b_v0.dat", kFALSE, "7tev_b_v0_loose.root", 0)


  // Step 2 and 3 are found in calibrate_de_dx.C


  // Step 5: create data
  CreateOutput("7tev_b.dat", kFALSE, "7tev_b.root", 0, "fitparameters/7tev_b.root")
  CreateOutputV0("7tev_b_v0.dat", kFALSE, "7tev_b_v0_loose.root", 0, "fitparameters/7tev_b.root")




  // Test examples
  // Step 1
  CreateCalib("7tev_b_test.dat", kFALSE, "7tev_b_test.root", 0, 2)  
  // Step 5
  CreateOutput("7tev_b_test.dat", kFALSE, "7tev_b_test.root", 0, "fitparameters/7tev_b_test.root")
  CreateOutputV0("7tev_b_test_v0.dat", kFALSE, "7tev_b_test_v0_loose.root", 0, "fitparameters/7tev_b_test.root")

 */


const Int_t nPtBins = 68;
//Double_t xBins[nPtBins+1] = {0., 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1 , 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 9.0, 10.0};
Double_t xBins[nPtBins+1] = {0. ,  0.05, 0.1,  0.15, 0.2,  0.25, 0.3,  0.35, 0.4,  0.45,
			     0.5,  0.55, 0.6,  0.65, 0.7,  0.75, 0.8,  0.85, 0.9,  0.95,
			     1.0,  1.1 , 1.2,  1.3 , 1.4,  1.5 , 1.6,  1.7 , 1.8,  1.9 ,
			     2.0,  2.2 , 2.4,  2.6 , 2.8,  3.0 , 3.2,  3.4 , 3.6,  3.8 ,
			     4.0,  4.5 , 5.0,  5.5 , 6.0,  6.5 , 7.0,  8.0 , 9.0,  10.0,
			     11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 18.0, 20.0, 22.0, 24.0,
			     26.0, 28.0, 30.0, 32.0, 34.0, 36.0, 40.0, 45.0, 50.0 };

TFile* dataOutFile = 0;
TFile* mcOutFile   = 0;

TF1* piFunc = 0;
TF1* kFunc  = 0;
TF1* pFunc = 0;
TF1* eFunc = 0;
TF1* sigmaFunc = 0;

void CreateOutput(const Char_t* dataFileName, Bool_t isMc, 
		  const Char_t* outfilename, Int_t maxEvents,
		  const Char_t* fitFileName=0);
void AddObject(TList* list, Int_t filter, Bool_t phiCut, Int_t run, 
	       Bool_t analyzeMc, Bool_t etaCut = kFALSE, Bool_t etaAbs = kFALSE, 
	       Int_t etaLow=-8, Int_t etaHigh=8);
void CreateOutputV0(const Char_t* dataFileName, Bool_t isMc, 
		    const Char_t* outFileName, Int_t maxEvents,
		    const Char_t* fitFileName=0);
void AddObjectV0(TList* list, const Char_t* baseName, Bool_t phiCut, Int_t run, 
		 Bool_t analyzeMc, Bool_t etaCut = kFALSE, Bool_t etaAbs = kFALSE, 
	       Int_t etaLow=-8, Int_t etaHigh=8);

void CreateCalib(const Char_t* dataFileName, Bool_t isMc, 
		 const Char_t* outfilename, Int_t maxEvents, Int_t startStep);
void AddCalibObject(TList* list, Int_t filter, Bool_t phiCut, Int_t run, 
		    Bool_t analyzeMc, Bool_t etaCut = kFALSE, Bool_t etaAbs = kFALSE, 
		    Int_t etaLow=-8, Int_t etaHigh=8);
void CreateCalibV0(const Char_t* dataFileName, Bool_t isMc, 
		   const Char_t* outFileName, Int_t maxEvents);
void AddCalibObjectV0(TList* list, const Char_t* baseName, Bool_t phiCut, 
		      Int_t run, Bool_t analyzeMc);
Int_t CalculateFilter(DeDxTrack* track);

//____________________________________________________________________________
void CreateOutput(const Char_t* dataFileName, Bool_t isMc, 
		  const Char_t* outFileName, Int_t maxEvents,
		  const Char_t* fitFileName)
{

  TF1* fDeDxVsEtaNeg = 0;
  TF1* fDeDxVsEtaPos = 0;
  TF1* fDeDxVsNcl    = 0;

  if(fitFileName) {

    TFile* fitFile = FindFileFresh(fitFileName);
    if(!fitFile)
      return;
    DeDxFitInfo* fitPar = (DeDxFitInfo*)fitFile->Get("fitInfo");
    fitPar->Print();
  

    if(!fitPar->calibFileName.IsNull()) {

      cout << "Setting calibFile: " << fitPar->calibFileName << endl;
      TFile* calibFile = FindFileFresh(fitPar->calibFileName);
      if(!calibFile)
	return;
      AliHighPtDeDxCalib* calib = (AliHighPtDeDxCalib*)GetObject(calibFile, 1, kTRUE, 0);
      calib->Print();
      fDeDxVsEtaNeg = calib->GetDeDxVsEtaNeg();
      fDeDxVsEtaPos = calib->GetDeDxVsEtaPos();
      fDeDxVsNcl    = calib->GetDeDxVsNcl();
    }
  
    fixMIP      = fitPar->MIP;
    fixPlateau  = fitPar->plateau;

    Double_t dedxPar[6]  = {0, 0, 0, 0, 0, 0};
    Double_t sigmaPar[6] = {0, 0, 0, 0, 0, 0};
    
    dedxPar[0] = fitPar->optionDeDx;
    for(Int_t i = 0; i < fitPar->nDeDxPar; i++) {
      dedxPar[i+1] = fitPar->parDeDx[i];
    }

    sigmaPar[0] = fitPar->optionSigma;
    for(Int_t i = 0; i < fitPar->nSigmaPar; i++) {
      sigmaPar[i+1] = fitPar->parSigma[i];
    }

    piFunc = new TF1("piFunc", FitFunc, 0, 100, fitPar->nDeDxPar+1);
    piFunc->SetParameters(dedxPar);

    kFunc = new TF1("kFunc", FitFunc, 0, 100, fitPar->nDeDxPar+1);
    kFunc->SetParameters(dedxPar);
    kFunc->SetParameter(0, kFunc->GetParameter(0)+10);

    pFunc = new TF1("pFunc", FitFunc, 0, 100, fitPar->nDeDxPar+1);
    pFunc->SetParameters(dedxPar);
    pFunc->SetParameter(0, pFunc->GetParameter(0)+20);

    eFunc = new TF1("eFunc", FitFunc, 0, 100, fitPar->nDeDxPar+1);
    eFunc->SetParameters(dedxPar);
    eFunc->SetParameter(0, eFunc->GetParameter(0)+30);

    sigmaFunc = new TF1("sigmaFunc", SigmaFunc, 0, 100, fitPar->nSigmaPar+1); 
    sigmaFunc->SetParameters(sigmaPar);
  }

  CreateDir("data");
  dataOutFile = new TFile(Form("data/%s", outFileName), "RECREATE");
  if(isMc) {
    
    CreateDir("mc");
    mcOutFile = new TFile(Form("mc/%s", outFileName), "RECREATE");
  }

  // Create output objects
  dataOutFile->cd();
  TList* runList = new TList();
  runList->SetOwner(kTRUE);
  runList->SetBit(TObject::kSingleKey);
  
  TList* analysisList = new TList();
  analysisList->SetOwner(kFALSE);
  
  //        TList    filter, phi cut, run, isMc
  AddObject(analysisList, 1, kTRUE,  0, isMc); 
  if(kTRUE) {
    AddObject(analysisList, 1, kTRUE,  0, isMc, kTRUE, kFALSE, -8, -6); 
    AddObject(analysisList, 1, kTRUE,  0, isMc, kTRUE, kFALSE, -6, -4); 
    AddObject(analysisList, 1, kTRUE,  0, isMc, kTRUE, kFALSE, -4, -2); 
    AddObject(analysisList, 1, kTRUE,  0, isMc, kTRUE, kFALSE, -2,  0); 
    AddObject(analysisList, 1, kTRUE,  0, isMc, kTRUE, kFALSE, 0, 2); 
    AddObject(analysisList, 1, kTRUE,  0, isMc, kTRUE, kFALSE, 2, 4); 
    AddObject(analysisList, 1, kTRUE,  0, isMc, kTRUE, kFALSE, 4, 6); 
    AddObject(analysisList, 1, kTRUE,  0, isMc, kTRUE, kFALSE, 6, 8); 
  }
  AddObject(analysisList, 1, kTRUE,  0, isMc, kTRUE, kFALSE, -8, 0); 
  AddObject(analysisList, 1, kTRUE,  0, isMc, kTRUE, kFALSE, 0, 8); 
  AddObject(analysisList, 1, kTRUE,  0, isMc, kFALSE, kTRUE, 0, 4); 
  AddObject(analysisList, 1, kTRUE,  0, isMc, kFALSE, kTRUE, 4, 8); 
  AddObject(analysisList, 1, kFALSE, 0, isMc); 
  AddObject(analysisList, 2, kTRUE,  0, isMc); 
  AddObject(analysisList, 2, kFALSE, 0, isMc); 
  
  TTree* Tree = 0;
  
  if(strstr(dataFileName, ".dat")) {
    
    AliXRDPROOFtoolkit tool;
    TChain* chain = tool.MakeChain(dataFileName,"tree", 0, 1000);
    //    chain->Lookup();
    Tree = chain;
  } else {
    TFile* dataFile = FindFileFresh(dataFileName);
    if(!dataFile)
      return;
    
    Tree = (TTree*)dataFile->Get("tree");
  }

  TClonesArray* trackArray = 0;
  TClonesArray* mcTrackArray = 0;
  DeDxEvent* event = 0;
  Tree->SetBranchAddress("event", &event);
  Tree->SetBranchAddress("trackGlobalPar"  , &trackArray);
  if(isMc)
    Tree->SetBranchAddress("trackMC"  , &mcTrackArray);

  Int_t nEvents = Tree->GetEntries();
  cout << "Number of events: " << nEvents << endl;
  
  if(maxEvents>0 && maxEvents < nEvents) {
    
    nEvents = maxEvents;
    cout << "N events was reduced to: " << maxEvents << endl;
  }

  Int_t currentRun = 0;
  Int_t nBad = 0;
  TIter* iter = new TIter(analysisList);
  
  for(Int_t n = 0; n < nEvents; n++) {
    
    Tree->GetEntry(n);
    
    if((n+1)%1000000==0)
      cout << "Event: " << n+1 << "/" << nEvents << endl;

    // A few bad entries
    if(event->run == -1) {
      nBad++;
      continue;
    }

    if(event->run != currentRun) {
      
      cout << "New run: " << event->run << endl;
      currentRun = event->run;
      
      // Check if run objects exist
      TObjString* runString = new TObjString(Form("%d", currentRun));
      if(!runList->FindObject(runString->GetString().Data())) {
	
	runList->Add(runString);
	
	//        TList    filter, phi cut, run, isMc
	AddObject(analysisList, 1, kTRUE,  currentRun, isMc); 
	AddObject(analysisList, 1, kFALSE, currentRun, isMc); 
	// AddObject(analysisList, 2, kTRUE,  currentRun, isMc); 
	// AddObject(analysisList, 2, kFALSE, currentRun, isMc); 
	
	// Is this really necessary?
	delete iter;
	iter = new TIter(analysisList);

      } else {

	delete runString;
      }
    }

    // iterate over analysis list
    iter->Reset();
    AliHighPtDeDxBase* analysis = 0;
    while ((analysis = dynamic_cast<AliHighPtDeDxBase*> (iter->Next()))) {
      
      // First we set the global properties
      // If I want to use a narrower vtx I have here to redefine my 
      // vtxstatus according to the new vtx range
      analysis->SetEventRun(event->run);
      analysis->SetEventMag(event->mag);
      analysis->SetEventVtxStatus(event->vtxstatus);
      analysis->SetEventTrigger(event->trig);

      Int_t vtxstatusmc = 0;
      if(TMath::Abs(event->zvtxMC) < 10.0)
	vtxstatusmc = 1;
      analysis->SetEventVtxStatusMc(vtxstatusmc);

      if(!analysis->EventAccepted()) // only checks for runs currently
	continue;

      // There is a small prolem in making this really nice, because we need
      // also event info from events that are rejected by the vtx cut to
      // correct our data. That is why we for bnow only have the run cut there
      
      analysis->FillEventInfo();
      
      // First we do the special MC part
      if(isMc) {
	
	AliHighPtDeDxMc* mcAnalysis = dynamic_cast<AliHighPtDeDxMc*> (analysis);
	if(mcAnalysis) {
	  if(vtxstatusmc) {
	    
	    // loop over mc tracks
	    const Int_t nMcTracks = mcTrackArray->GetEntries();
	    
	    for(Int_t i = 0; i < nMcTracks; i++) {
	      
	      DeDxTrackMC* trackMC = (DeDxTrackMC*)mcTrackArray->At(i);
	      
	      mcAnalysis->SetTrackPtMc(trackMC->ptMC);
	      mcAnalysis->SetTrackEtaMc(trackMC->etaMC);
	      mcAnalysis->SetTrackChargeMc(trackMC->qMC);
	      mcAnalysis->SetTrackPidMc(trackMC->pidMC);
	      if(mcAnalysis->TrackAcceptedMc()) {
		// if(trackMC->ptMC<2.0)
		//   mcAnalysis->FillTrackInfoMc(100.0);
		// else
		  mcAnalysis->FillTrackInfoMc(1.0);
	      }
	    }
	  }
	}
      }
      
      // The trig==1 is always true for real data, but not for MC data
      if(!event->trig)
	continue;
	
      if(event->vtxstatus<1) // only fill tracks for events with vtx inside cuts
	continue;
      
      const Int_t nTracks = trackArray->GetEntries();
      
      for(Int_t i = 0; i < nTracks; i++) {
	
	DeDxTrack* track = (DeDxTrack*)trackArray->At(i);
	
	Double_t eta  = track->eta;
	Double_t dedx = track->dedx;
	Int_t ncl     = track->ncl;
	if(fDeDxVsEtaNeg) {

	  dedx *= 50.0 / fDeDxVsNcl->Eval(ncl);

	  if(eta < 0) 
	    dedx *= 50.0 / fDeDxVsEtaNeg->Eval(eta);
	  else
	    dedx *= 50.0 / fDeDxVsEtaPos->Eval(eta);
	}
	
	analysis->SetTrackCharge(track->q);
	analysis->SetTrackP(track->p);
	analysis->SetTrackPt(track->pt);
	//
	Int_t filter = CalculateFilter(track);
	analysis->SetTrackFilter(filter);
	analysis->SetTrackPhi(track->phi);
	analysis->SetTrackDeDx(dedx);
	analysis->SetTrackNcl(ncl);
	analysis->SetTrackEta(eta);
	analysis->SetTrackPidMc(track->pid);

	if(analysis->TrackAccepted()) {
	  // if(track->pt<2.0)
	  //   analysis->FillTrackInfo(100.0);
	  // else
	    analysis->FillTrackInfo(1.0);
	}
      }
    }
  }
  
  if(isMc) {

    TList* runListMc = (TList*)runList->Clone();
    mcOutFile->cd();
    runListMc->Write();

    // I need to somehow specify autowrite, but I am not sure how, so for now
    // this a bit stupid workaround...
    iter->Reset();
    AliHighPtDeDxBase* analysis = 0;
    while ((analysis = dynamic_cast<AliHighPtDeDxBase*> (iter->Next()))) {
      
      AliHighPtDeDxMc* mcAnalysis = dynamic_cast<AliHighPtDeDxMc*> (analysis);
      if(mcAnalysis) {
	mcAnalysis->Write();
      }
    }
    mcOutFile->Close();
    delete mcOutFile;
    mcOutFile = 0;
  }

  dataOutFile->cd();
  runList->Write("runList");
  iter->Reset();
  AliHighPtDeDxBase* analysis = 0;
  while ((analysis = dynamic_cast<AliHighPtDeDxBase*> (iter->Next()))) {
    
    AliHighPtDeDxData* dataAnalysis = dynamic_cast<AliHighPtDeDxData*> (analysis);
    if(dataAnalysis) {
      dataAnalysis->Write();
    }
  }
  dataOutFile->Close();
  delete dataOutFile;
  dataOutFile = 0;

  cout << "Nbad (runno == -1) : " << nBad << endl;
}

//___________________________________________________________________________
void AddObject(TList* list, Int_t filter, Bool_t phiCut, Int_t run, 
	       Bool_t analyzeMc, Bool_t etaCut, Bool_t etaAbs, 
	       Int_t etaLow, Int_t etaHigh)
{
  if(etaCut && etaAbs) {

    Fatal("AddObject", "You cannot have a cut on both abs and signed eta");
    return;
  }
  TString objectName("filter");
  objectName += filter;
  if(phiCut)
    objectName += "phicut";
  if(run) {
    objectName += "_";
    objectName += run;
  }
  if(etaCut) {
    objectName += "eta";
    objectName += etaLow;
    objectName += etaHigh;
  }
  if(etaAbs) {
    objectName += "etaabs";
    objectName += etaLow;
    objectName += etaHigh;
  }

  dataOutFile->cd();
  AliHighPtDeDxData* data = new AliHighPtDeDxData(objectName.Data(), "Analysis data");
  if(run) {
    data->SetUseRunCut(kTRUE);
    data->SetRun(run);
  }

  list->Add(data);
  data->SetIsMc(analyzeMc);
  data->SetUseFilterCut(kTRUE);
  data->SetFilter(filter);
  data->SetUsePhiCut(phiCut);
  if(phiCut) {
    data->SetPhiCutLow(data->GetStandardPhiCutLow());
    data->SetPhiCutHigh(data->GetStandardPhiCutHigh());
  }
  if(etaCut)
    data->SetUseEtaCut(kTRUE);
  if(etaAbs)
    data->SetUseEtaCutAbs(kTRUE);
  if(etaCut || etaAbs) {
    data->SetEtaLow(Double_t(etaLow)/10.0);
    data->SetEtaHigh(Double_t(etaHigh)/10.0);
  }

  data->Init(nPtBins, xBins);
  data->Print();

  if(piFunc && kFunc && pFunc && eFunc && sigmaFunc) {

    data->SetPionDeDxFunction(piFunc);
    data->SetKaonDeDxFunction(kFunc);
    data->SetProtonDeDxFunction(pFunc);
    data->SetElectronDeDxFunction(eFunc);
    data->SetSigmaDeDxFunction(sigmaFunc);
  }

  if(analyzeMc) {
    // create the correction class also
    
    mcOutFile->cd();
    AliHighPtDeDxMc* mc = new AliHighPtDeDxMc(objectName.Data(), "Analysis mc");
    if(run) {
      mc->SetUseRunCut(kTRUE);
      mc->SetRun(run);
    }
    
    list->Add(mc);
    mc->SetUseFilterCut(kTRUE);
    mc->SetFilter(filter);
    mc->SetUsePhiCut(phiCut);
    if(phiCut) {
      mc->SetPhiCutLow(mc->GetStandardPhiCutLow());
      mc->SetPhiCutHigh(mc->GetStandardPhiCutHigh());
    }
    mc->Init(nPtBins, xBins);
    mc->Print();
  }
}

//____________________________________________________________________________
void CreateOutputV0(const Char_t* dataFileName, Bool_t isMc, 
		    const Char_t* outFileName, Int_t maxEvents,
		    const Char_t* fitFileName)
{

  TF1* fDeDxVsEtaNeg = 0;
  TF1* fDeDxVsEtaPos = 0;
  TF1* fDeDxVsNcl    = 0;

  if(fitFileName) {

    TFile* fitFile = FindFileFresh(fitFileName);
    if(!fitFile)
      return;
    DeDxFitInfo* fitPar = (DeDxFitInfo*)fitFile->Get("fitInfo");
    fitPar->Print();
  

    if(!fitPar->calibFileName.IsNull()) {

      cout << "Setting calibFile: " << fitPar->calibFileName << endl;
      TFile* calibFile = FindFileFresh(fitPar->calibFileName);
      if(!calibFile)
	return;
      AliHighPtDeDxCalib* calib = (AliHighPtDeDxCalib*)GetObject(calibFile, 1, kTRUE, 0);
      calib->Print();
      fDeDxVsEtaNeg = calib->GetDeDxVsEtaNeg();
      fDeDxVsEtaPos = calib->GetDeDxVsEtaPos();
      fDeDxVsNcl    = calib->GetDeDxVsNcl();
    }
  
    fixMIP      = fitPar->MIP;
    fixPlateau  = fitPar->plateau;

    Double_t dedxPar[6]  = {0, 0, 0, 0, 0, 0};
    Double_t sigmaPar[6] = {0, 0, 0, 0, 0, 0};
    
    dedxPar[0] = fitPar->optionDeDx;
    for(Int_t i = 0; i < fitPar->nDeDxPar; i++) {
      dedxPar[i+1] = fitPar->parDeDx[i];
    }

    sigmaPar[0] = fitPar->optionSigma;
    for(Int_t i = 0; i < fitPar->nSigmaPar; i++) {
      sigmaPar[i+1] = fitPar->parSigma[i];
    }

    piFunc = new TF1("piFunc", FitFunc, 0, 100, fitPar->nDeDxPar+1);
    piFunc->SetParameters(dedxPar);

    kFunc = new TF1("kFunc", FitFunc, 0, 100, fitPar->nDeDxPar+1);
    kFunc->SetParameters(dedxPar);
    kFunc->SetParameter(0, kFunc->GetParameter(0)+10);

    pFunc = new TF1("pFunc", FitFunc, 0, 100, fitPar->nDeDxPar+1);
    pFunc->SetParameters(dedxPar);
    pFunc->SetParameter(0, pFunc->GetParameter(0)+20);

    eFunc = new TF1("eFunc", FitFunc, 0, 100, fitPar->nDeDxPar+1);
    eFunc->SetParameters(dedxPar);
    eFunc->SetParameter(0, eFunc->GetParameter(0)+30);

    sigmaFunc = new TF1("sigmaFunc", SigmaFunc, 0, 100, fitPar->nSigmaPar+1); 
    sigmaFunc->SetParameters(sigmaPar);
  }

  CreateDir("data");
  dataOutFile = new TFile(Form("data/%s", outFileName), "RECREATE");

  // Create output objects
  dataOutFile->cd();
  TList* runList = new TList();
  runList->SetOwner(kTRUE);
  runList->SetBit(TObject::kSingleKey);
  
  TList* analysisList = new TList();
  analysisList->SetOwner(kFALSE);
  
  //               TList      v0type, phi cut, run, isMc
  AddObjectV0(analysisList, "lambda", kTRUE,  0, isMc); 
  AddObjectV0(analysisList, "lambda", kTRUE,  0, isMc, kTRUE, kFALSE, -8, 0); 
  AddObjectV0(analysisList, "lambda", kTRUE,  0, isMc, kTRUE, kFALSE, 0, 8); 
  AddObjectV0(analysisList, "lambda", kTRUE,  0, isMc, kFALSE, kTRUE, 0, 4); 
  AddObjectV0(analysisList, "lambda", kTRUE,  0, isMc, kFALSE, kTRUE, 4, 8); 
  AddObjectV0(analysisList, "kaon", kTRUE,  0, isMc); 
  AddObjectV0(analysisList, "kaon", kTRUE,  0, isMc, kTRUE, kFALSE, -8, 0); 
  AddObjectV0(analysisList, "kaon", kTRUE,  0, isMc, kTRUE, kFALSE, 0, 8); 
  AddObjectV0(analysisList, "kaon", kTRUE,  0, isMc, kFALSE, kTRUE, 0, 4); 
  AddObjectV0(analysisList, "kaon", kTRUE,  0, isMc, kFALSE, kTRUE, 4, 8); 
  
  TTree* Tree = 0;
  
  if(strstr(dataFileName, ".dat")) {
    
    AliXRDPROOFtoolkit tool;
    TChain* chain = tool.MakeChain(dataFileName,"tree", 0, 1000);
    //    chain->Lookup();
    Tree = chain;
  } else {
    TFile* dataFile = FindFileFresh(dataFileName);
    if(!dataFile)
      return;
    
    Tree = (TTree*)dataFile->Get("tree");
  }

  TClonesArray* v0Array = 0;
  DeDxEvent* event = 0;
  Tree->SetBranchAddress("event", &event);
  Tree->SetBranchAddress("v0GlobalPar"  , &v0Array);

  Int_t nEvents = Tree->GetEntries();
  cout << "Number of events: " << nEvents << endl;
  
  if(maxEvents>0 && maxEvents < nEvents) {
    
    nEvents = maxEvents;
    cout << "N events was reduced to: " << maxEvents << endl;
  }

  Int_t currentRun = 0;
  Int_t nBad = 0;
  TIter* iter = new TIter(analysisList);
  
  for(Int_t n = 0; n < nEvents; n++) {
    
    Tree->GetEntry(n);
    
    if((n+1)%1000000==0)
      cout << "Event: " << n+1 << "/" << nEvents << endl;

    if(event->run == -1) {
      nBad++;
      continue;
    }
    // if(event->run == 126437)
    //   continue;
    
    if(event->run != currentRun) {
      
      cout << "New run: " << event->run << endl;
      currentRun = event->run;
      
      // Check if run objects exist
      TObjString* runString = new TObjString(Form("%d", currentRun));
      if(!runList->FindObject(runString->GetString().Data())) {
	
	runList->Add(runString);
	
	//               TList      v0type, phi cut, run, isMc
	AddObjectV0(analysisList, "lambda", kTRUE,  currentRun, isMc); 
	AddObjectV0(analysisList, "kaon", kTRUE,  currentRun, isMc); 
	
	// Is this really necessary?
	delete iter;
	iter = new TIter(analysisList);

      } else {

	delete runString;
      }
    }

    // iterate over analysis list
    iter->Reset();
    AliHighPtDeDxBase* analysis = 0;
    while ((analysis = dynamic_cast<AliHighPtDeDxBase*> (iter->Next()))) {
      
      // First we set the global properties
      // If I want to use a narrower vtx I have here to redefine my 
      // vtxstatus according to the new vtx range
      analysis->SetEventRun(event->run);
      analysis->SetEventMag(event->mag);
      analysis->SetEventVtxStatus(event->vtxstatus);
      analysis->SetEventTrigger(event->trig);

      Int_t vtxstatusmc = 0;
      if(TMath::Abs(event->zvtxMC) < 10.0)
	vtxstatusmc = 1;
      analysis->SetEventVtxStatusMc(vtxstatusmc);

      if(!analysis->EventAccepted()) // only checks for runs currently
	continue;
      // There is a small prolem in making this really nice, because we need
      // also event info from events that are rejected by the vtx cut to
      // correct our data. That is why we for bnow only have the run cut there
      
      analysis->FillEventInfo();
      
      
      // The trig==1 is always true for real data, but not for MC data
      if(!event->trig)
	continue;
	
      if(event->vtxstatus<1) // only fill tracks for events with vtx inside cuts
	continue;
      
      const Int_t nV0s = v0Array->GetEntries();
      
      for(Int_t i = 0; i < nV0s; i++) {
	
	DeDxV0* v0 = (DeDxV0*)v0Array->At(i);
	
	// if(v0->ptrack.filter != 2)
	//   continue;

	if(v0->status != 0)
	  continue;

	if(v0->dmassG < 0.1)
	  continue;

	// if(v0->decayr < 5.0)
	//   continue;

	// if(v0->decayr > 40.0)
	//   continue;
	
	const Double_t dmassK  = TMath::Abs(v0->dmassK0);
	const Double_t dmassL  = TMath::Abs(v0->dmassL);
	const Double_t dmassAL = TMath::Abs(v0->dmassAL);
	
	Bool_t fillPos = kFALSE;
	Bool_t fillNeg = kFALSE;
					   
	
	if(strstr(analysis->GetName(), "lambda")) {
	  
	  if(dmassK<0.01)
	    continue;

	  if(dmassL<0.01&&dmassL>0.01)
	    fillPos = kTRUE;

	  if(dmassAL<0.01&&dmassL)
	    fillNeg = kTRUE;
	} else { // kaons

	  if(dmassL<0.01 || dmassAL<0.01)
	    continue;

	  if(dmassK<0.01) {
	    fillPos = kTRUE;
	    fillNeg = kTRUE;
	  }
	}

	for(Int_t j = 0; j < 2; j++) {

	  DeDxTrack* track = 0;
	  
	  if(j==0) {

	    if(fillNeg)
	      track = &(v0->ntrack);
	    else
	      continue;
	  } else {

	    if(fillPos)
	      track = &(v0->ptrack);
	    else
	      continue;
	  }

	  Double_t eta  = track->eta;
	  Double_t dedx = track->dedx;
	  Int_t ncl     = track->ncl;
	  if(fDeDxVsEtaNeg) {
	    
	    dedx *= 50.0 / fDeDxVsNcl->Eval(ncl);
	    
	    if(eta < 0) 
	      dedx *= 50.0 / fDeDxVsEtaNeg->Eval(eta);
	    else
	      dedx *= 50.0 / fDeDxVsEtaPos->Eval(eta);
	  }
	  
	  analysis->SetTrackCharge(track->q);
	  analysis->SetTrackP(track->p);
	  analysis->SetTrackPt(track->pt);
	  // NB! Not used
	  analysis->SetTrackFilter(track->filter);
	  analysis->SetTrackPhi(track->phi);
	  analysis->SetTrackDeDx(dedx);
	  analysis->SetTrackNcl(ncl);
	  analysis->SetTrackBeta(track->beta);
	  analysis->SetTrackPidMc(track->pid);
	  analysis->SetTrackEta(eta);
	  
	  if(analysis->TrackAccepted()) {
	    analysis->FillTrackInfo(1.0);
	  }
	}
      }
    }
  }
  
  dataOutFile->cd();
  runList->Write("runList");
  iter->Reset();
  AliHighPtDeDxBase* analysis = 0;
  while ((analysis = dynamic_cast<AliHighPtDeDxBase*> (iter->Next()))) {
    
    AliHighPtDeDxData* dataAnalysis = dynamic_cast<AliHighPtDeDxData*> (analysis);
    if(dataAnalysis) {
      dataAnalysis->Write();
    }
  }
  dataOutFile->Close();
  delete dataOutFile;
  dataOutFile = 0;

  cout << "Nbad (runno == -1) : " << nBad << endl;
}


//___________________________________________________________________________
void AddObjectV0(TList* list, const Char_t* baseName, Bool_t phiCut, Int_t run, 
		 Bool_t analyzeMc, Bool_t etaCut, Bool_t etaAbs, 
		 Int_t etaLow, Int_t etaHigh)
{
  if(etaCut && etaAbs) {

    Fatal("AddObject", "You cannot have a cut on both abs and signed eta");
    return;
  }

  TString objectName(baseName);
  if(phiCut)
    objectName += "phicut";
  if(run) {
    objectName += "_";
    objectName += run;
  }
  if(etaCut) {
    objectName += "eta";
    objectName += etaLow;
    objectName += etaHigh;
  }
  if(etaAbs) {
    objectName += "etaabs";
    objectName += etaLow;
    objectName += etaHigh;
  }
  dataOutFile->cd();
  AliHighPtDeDxData* data = new AliHighPtDeDxData(objectName.Data(), "Analysis data");
  if(run) {
    data->SetUseRunCut(kTRUE);
    data->SetRun(run);
  }
  
  list->Add(data);
  data->SetIsMc(analyzeMc);
  data->SetUseFilterCut(kFALSE);

  data->SetUsePhiCut(phiCut);
  if(phiCut) {
    data->SetPhiCutLow(data->GetStandardPhiCutLow());
    data->SetPhiCutHigh(data->GetStandardPhiCutHigh());
  }
  if(etaCut)
    data->SetUseEtaCut(kTRUE);
  if(etaAbs)
    data->SetUseEtaCutAbs(kTRUE);
  if(etaCut || etaAbs) {
    data->SetEtaLow(Double_t(etaLow)/10.0);
    data->SetEtaHigh(Double_t(etaHigh)/10.0);
  }

  data->Init(nPtBins, xBins);
  data->Print();

  if(piFunc && kFunc && pFunc && eFunc && sigmaFunc) {

    data->SetPionDeDxFunction(piFunc);
    data->SetKaonDeDxFunction(kFunc);
    data->SetProtonDeDxFunction(pFunc);
    data->SetElectronDeDxFunction(eFunc);
    data->SetSigmaDeDxFunction(sigmaFunc);
  }
}

//____________________________________________________________________________
void CreateCalib(const Char_t* dataFileName, Bool_t isMc, 
		 const Char_t* outFileName, Int_t maxEvents, Int_t startStep)
{
  if(startStep < 1 || startStep > 3) {

    cout << "Start step should be 1,2, or 3 - 2 is recommended" << endl;
  }

  if(startStep == 1) {

    for(Int_t i = 0; i < 10; i++)
      cout << "Step 1 calibration is depreceated! - 2 is recommended" << endl;
  }
  
  if(startStep<3) {

    CreateDir("calib_eta");
    dataOutFile = new TFile(Form("calib_eta/%s", outFileName), "RECREATE");
  } else {

    CreateDir("no_calib_eta");
    dataOutFile = new TFile(Form("no_calib_eta/%s", outFileName), "RECREATE");
  }

  // Create output objects
  dataOutFile->cd();
  TList* runList = new TList();
  runList->SetOwner(kTRUE);
  runList->SetBit(TObject::kSingleKey);
  
  TList* calibList = new TList();
  calibList->SetOwner(kFALSE);
  
  //        TList    filter, phi cut, run, isMc
  AddCalibObject(calibList, 1, kTRUE,  0, isMc); 
  AddCalibObject(calibList, 2, kTRUE,  0, isMc); 
  AddCalibObject(calibList, 1, kTRUE,  0, isMc, kTRUE, kFALSE, -8, 0); 
  AddCalibObject(calibList, 1, kTRUE,  0, isMc, kTRUE, kFALSE, 0, 8); 
  AddCalibObject(calibList, 1, kTRUE,  0, isMc, kFALSE, kTRUE, 0, 4); 
  AddCalibObject(calibList, 1, kTRUE,  0, isMc, kFALSE, kTRUE, 4, 8); 
  //  AddCalibObject(calibList, 1, kFALSE, 0, isMc); 
  //  AddCalibObject(calibList, 2, kTRUE,  0, isMc); 
  //  AddCalibObject(calibList, 2, kFALSE, 0, isMc); 
  
  TTree* Tree = 0;
  
  if(strstr(dataFileName, ".dat")) {
    
    AliXRDPROOFtoolkit tool;
    TChain* chain = tool.MakeChain(dataFileName,"tree", 0, 1000);
    //    chain->Lookup();
    Tree = chain;
  } else {
    TFile* dataFile = FindFileFresh(dataFileName);
    if(!dataFile)
      return;
    
    Tree = (TTree*)dataFile->Get("tree");
  }

  TClonesArray* trackArray = 0;
  TClonesArray* mcTrackArray = 0;
  DeDxEvent* event = 0;
  Tree->SetBranchAddress("event", &event);
  Tree->SetBranchAddress("trackGlobalPar"  , &trackArray);
  if(isMc)
    Tree->SetBranchAddress("trackMC"  , &mcTrackArray);

  Int_t nEvents = Tree->GetEntries();
  cout << "Number of events: " << nEvents << endl;
  
  if(maxEvents>0 && maxEvents < nEvents) {
    
    nEvents = maxEvents;
    cout << "N events was reduced to: " << maxEvents << endl;
  }

  Int_t currentRun = 0;
  Int_t nBad = 0;
  TIter* iter = new TIter(calibList);
  AliHighPtDeDxCalib* calib = 0;

  for(Int_t step = startStep; step <= 3; step++) {
    
    iter->Reset();
    while ((calib = dynamic_cast<AliHighPtDeDxCalib*> (iter->Next()))) {
      
      calib->SetStep(step);
      // Here you can set the MIP interval
      // default is 40-60
      if(step==startStep) {
	calib->SetDeDxMIPMin(25);
	calib->SetDeDxMIPMax(55);
      } else {
	calib->SetDeDxMIPMin(40);
	calib->SetDeDxMIPMax(60);
      }
      calib->Init(step, nPtBins, xBins);
    }
    
    for(Int_t n = 0; n < nEvents; n++) {
      
      Tree->GetEntry(n);
      
      if((n+1)%1000000==0)
	cout << "Event: " << n+1 << "/" << nEvents << endl;
      
      if(event->run == -1) {
	nBad++;
	continue;
      }
      // if(event->run == 126437)
      //   continue;
      
      if(event->run != currentRun) {
	
	cout << "New run: " << event->run << endl;
	currentRun = event->run;
	
	// Check if run objects exist
	TObjString* runString = new TObjString(Form("%d", currentRun));
	if(!runList->FindObject(runString->GetString().Data())) {
	  
	  runList->Add(runString);
	  
	  //        TList    filter, phi cut, run, isMc
	  AddCalibObject(calibList, 1, kTRUE,  currentRun, isMc); 
	  // AddCalibObject(calibList, 1, kFALSE, currentRun, isMc); 
	  // AddCalibObject(calibList, 2, kTRUE,  currentRun, isMc); 
	  // AddCalibObject(calibList, 2, kFALSE, currentRun, isMc); 
	  
	  // Is this really necessary?
	  delete iter;
	  iter = new TIter(calibList);

	  iter->Reset();
	  while ((calib = dynamic_cast<AliHighPtDeDxCalib*> (iter->Next()))) {
	    
	    calib->SetStep(step);
	    if(step==startStep) {
	      calib->SetDeDxMIPMin(25);
	      calib->SetDeDxMIPMax(55);
	    } else {
	      calib->SetDeDxMIPMin(40);
	      calib->SetDeDxMIPMax(60);
	    }
	    calib->Init(step, nPtBins, xBins);
	  }
	  
	} else {

	  delete runString;
	}
      }
      
      // iterate over calib list
      iter->Reset();
      while ((calib = dynamic_cast<AliHighPtDeDxCalib*> (iter->Next()))) {
	
	// First we set the global properties
	// If I want to use a narrower vtx I have here to redefine my 
	// vtxstatus according to the new vtx range
	calib->SetEventRun(event->run);
	calib->SetEventMag(event->mag);
	calib->SetEventVtxStatus(event->vtxstatus);
	calib->SetEventTrigger(event->trig);
	
	if(!calib->EventAccepted()) // only checks for runs currently
	  continue;

	calib->FillEventInfo();
	
	// The trig==1 is always true for real data, but not for MC data
	if(!event->trig)
	  continue;
	
	if(event->vtxstatus<1) // only fill tracks for events with vtx inside cuts
	  continue;
	
	const Int_t nTracks = trackArray->GetEntries();
	
	for(Int_t i = 0; i < nTracks; i++) {
	  
	  DeDxTrack* track = (DeDxTrack*)trackArray->At(i);
	  
	  calib->SetTrackCharge(track->q);
	  calib->SetTrackP(track->p);
	  calib->SetTrackPt(track->pt);
	  Int_t filter = CalculateFilter(track);
	  calib->SetTrackFilter(filter);
	  calib->SetTrackPhi(track->phi);
	  calib->SetTrackDeDx(track->dedx);
	  calib->SetTrackNcl(track->ncl);
	  calib->SetTrackBeta(track->beta);
	  calib->SetTrackPidMc(track->pid);
	  calib->SetTrackEta(track->eta);
	  
	  if(calib->TrackAccepted()) {
	    // if(track->pt<2.0)
	    //   calib->FillTrackInfo(100.0);
	    // else
	      calib->FillTrackInfo(1.0);
	  }
	}
      }
    }

    if(step == 2) { // perform eta calibration

      cout << "Starting Eta cal" << endl;
      iter->Reset();
      while ((calib = dynamic_cast<AliHighPtDeDxCalib*> (iter->Next()))) {
	
	calib->PerformEtaCal();
      }
      cout << "Finishing Eta cal" << endl;
    } else if(step == 1) { // perform ncl calibration
      

      cout << "Starting Ncl cal" << endl;
      iter->Reset();
      while ((calib = dynamic_cast<AliHighPtDeDxCalib*> (iter->Next()))) {
	
	calib->PerformNclCal();
      }
      cout << "Finished Ncl cal" << endl;
    }
  }

  cout << "Nbad (runno == -1) : " << nBad << endl;

  dataOutFile->cd();
  runList->Write("runList");
  iter->Reset();
  cout << "Writing calibration data" << endl;
  while ((calib = dynamic_cast<AliHighPtDeDxCalib*> (iter->Next()))) {
    
    cout << "Writing calibration object " << calib->GetName() << endl;
    calib->Write();
  }
  cout << "Finsihed writing calibration data" << endl;
  dataOutFile->Close();
  delete dataOutFile;
  dataOutFile = 0;
}

//___________________________________________________________________________
void AddCalibObject(TList* list, Int_t filter, Bool_t phiCut, Int_t run, 
		    Bool_t analyzeMc, Bool_t etaCut, Bool_t etaAbs, 
		    Int_t etaLow, Int_t etaHigh)
{
  if(etaCut && etaAbs) {

    Fatal("AddCalibObject", "You cannot have a cut on both abs and signed eta");
    return;
  }
  TString objectName("filter");
  objectName += filter;
  if(phiCut)
    objectName += "phicut";
  if(run) {
    objectName += "_";
    objectName += run;
  }
  if(etaCut) {
    objectName += "eta";
    objectName += etaLow;
    objectName += etaHigh;
  }
  if(etaAbs) {
    objectName += "etaabs";
    objectName += etaLow;
    objectName += etaHigh;
  }

  //  dataOutFile->cd();
  AliHighPtDeDxCalib* data = new AliHighPtDeDxCalib(objectName.Data(), "Calib data");
  if(run) {
    data->SetUseRunCut(kTRUE);
    data->SetRun(run);
  }

  list->Add(data);
  data->SetIsMc(analyzeMc);
  data->SetUseFilterCut(kTRUE);
  data->SetFilter(filter);
  data->SetUsePhiCut(phiCut);
  if(phiCut) {
    data->SetPhiCutLow(data->GetStandardPhiCutLow());
    data->SetPhiCutHigh(data->GetStandardPhiCutHigh());
  }
  if(etaCut)
    data->SetUseEtaCut(kTRUE);
  if(etaAbs)
    data->SetUseEtaCutAbs(kTRUE);
  if(etaCut || etaAbs) {
    data->SetEtaLow(Double_t(etaLow)/10.0);
    data->SetEtaHigh(Double_t(etaHigh)/10.0);
  }
  if(run) {
    data->SetUseRunCut(kTRUE);
    data->SetRun(run);
  }

  data->Init(nPtBins, xBins);
  data->Print();
}

//____________________________________________________________________________
void CreateCalibV0(const Char_t* dataFileName, Bool_t isMc, 
		   const Char_t* outFileName, Int_t maxEvents)
{
  CreateDir("no_calib_eta");
  dataOutFile = new TFile(Form("no_calib_eta/%s", outFileName), "RECREATE");
  
  // Create output objects
  dataOutFile->cd();
  TList* runList = new TList();
  runList->SetOwner(kTRUE);
  runList->SetBit(TObject::kSingleKey);
  
  TList* calibList = new TList();
  calibList->SetOwner(kFALSE);
  
  //        TList    filter, phi cut, run, isMc
  AddCalibObjectV0(calibList, "lambda", kTRUE,  0, isMc); 
  //  AddCalibObject(calibList, 1, kFALSE, 0, isMc); 
  AddCalibObjectV0(calibList, "kaon", kTRUE,  0, isMc); 
  //  AddCalibObject(calibList, 2, kFALSE, 0, isMc); 
  
  TTree* Tree = 0;
  
  if(strstr(dataFileName, ".dat")) {
    
    AliXRDPROOFtoolkit tool;
    TChain* chain = tool.MakeChain(dataFileName,"tree", 0, 1000);
    //    chain->Lookup();
    Tree = chain;
  } else {
    TFile* dataFile = FindFileFresh(dataFileName);
    if(!dataFile)
      return;
    
    Tree = (TTree*)dataFile->Get("tree");
  }
  
  TClonesArray* v0Array = 0;
  DeDxEvent* event = 0;
  Tree->SetBranchAddress("event", &event);
  Tree->SetBranchAddress("v0"  , &v0Array);
  
  Int_t nEvents = Tree->GetEntries();
  cout << "Number of events: " << nEvents << endl;
  
  if(maxEvents>0 && maxEvents < nEvents) {
    
    nEvents = maxEvents;
    cout << "N events was reduced to: " << maxEvents << endl;
  }
  
  Int_t currentRun = 0;
  Int_t nBad = 0;
  TIter* iter = new TIter(calibList);
  AliHighPtDeDxCalib* calib = 0;
  
  iter->Reset();
  while ((calib = dynamic_cast<AliHighPtDeDxCalib*> (iter->Next()))) {
    
    calib->SetStep(3);
  }
  
  for(Int_t n = 0; n < nEvents; n++) {
    
    Tree->GetEntry(n);
    
    if((n+1)%1000000==0)
      cout << "Event: " << n+1 << "/" << nEvents << endl;
    
    if(event->run == -1) {
      nBad++;
      continue;
    }
    // if(event->run == 126437)
    //   continue;
    
    if(event->run != currentRun) {
      
      cout << "New run: " << event->run << endl;
      currentRun = event->run;
      
      // Check if run objects exist
      TObjString* runString = new TObjString(Form("%d", currentRun));
      if(!runList->FindObject(runString->GetString().Data())) {
	
	runList->Add(runString);
	
	//        TList    filter, phi cut, run, isMc
	AddCalibObjectV0(calibList, "lambda", kTRUE,  currentRun, isMc); 
	AddCalibObjectV0(calibList, "kaon", kTRUE,  currentRun, isMc); 
	
	// Is this really necessary?
	delete iter;
	iter = new TIter(calibList);
	
	iter->Reset();
	while ((calib = dynamic_cast<AliHighPtDeDxCalib*> (iter->Next()))) {
	  
	  calib->SetStep(3);
	}
	
	} else {
	
	delete runString;
      }
    }
    
    // iterate over calib list
    iter->Reset();
    while ((calib = dynamic_cast<AliHighPtDeDxCalib*> (iter->Next()))) {
      
      // First we set the global properties
      // If I want to use a narrower vtx I have here to redefine my 
      // vtxstatus according to the new vtx range
      calib->SetEventRun(event->run);
      calib->SetEventMag(event->mag);
      calib->SetEventVtxStatus(event->vtxstatus);
      calib->SetEventTrigger(event->trig);
      
      if(!calib->EventAccepted()) // only checks for runs currently
	continue;

      calib->FillEventInfo();
	
      // The trig==1 is always true for real data, but not for MC data
      if(!event->trig)
	continue;
      
      if(event->vtxstatus<1) // only fill tracks for events with vtx inside cuts
	continue;
      
      const Int_t nV0s = v0Array->GetEntries();
      
      for(Int_t i = 0; i < nV0s; i++) {
	
	DeDxV0* v0 = (DeDxV0*)v0Array->At(i);
	
	// if(v0->ptrack.filter != 2)
	//   continue;

	if(v0->status != 0)
	  continue;

	if(v0->dmassG < 0.1)
	  continue;

	// if(v0->decayr < 5.0)
	//   continue;

	// if(v0->decayr > 40.0)
	//   continue;
	
	const Double_t dmassK  = TMath::Abs(v0->dmassK0);
	const Double_t dmassL  = TMath::Abs(v0->dmassL);
	const Double_t dmassAL = TMath::Abs(v0->dmassAL);
	
	Bool_t fillPos = kFALSE;
	Bool_t fillNeg = kFALSE;
					   

	if(strstr(calib->GetName(), "lambda")) {
	  
	  if(dmassK<0.01)
	    continue;

	  if(dmassL<0.01)
	    fillPos = kTRUE;

	  if(dmassAL<0.01)
	    fillNeg = kTRUE;
	} else { // kaons

	  if(dmassL<0.01 || dmassAL<0.01)
	    continue;

	  if(dmassK<0.01) {
	    fillPos = kTRUE;
	    fillNeg = kTRUE;
	  }
	}

	for(Int_t j = 0; j < 2; j++) {

	  DeDxTrack* track = 0;
	  
	  if(j==0) {

	    if(fillNeg)
	      track = &(v0->ntrack);
	    else
	      continue;
	  } else {

	    if(fillPos)
	      track = &(v0->ptrack);
	    else
	      continue;
	  }
	  
	  calib->SetTrackCharge(track->q);
	  calib->SetTrackP(track->p);
	  calib->SetTrackPt(track->pt);
	  // NB! Filter is not used for these tracks!
	  calib->SetTrackFilter(track->filter);
	  calib->SetTrackPhi(track->phi);
	  calib->SetTrackDeDx(track->dedx);
	  calib->SetTrackNcl(track->ncl);
	  calib->SetTrackBeta(track->beta);
	  calib->SetTrackPidMc(track->pid);
	  calib->SetTrackEta(track->eta);
	  
	  if(calib->TrackAccepted()) {
	    calib->FillTrackInfo(1.0);
	  }
	}
      }
    }
  }
  dataOutFile->cd();
  runList->Write("runList");
  iter->Reset();
  while ((calib = dynamic_cast<AliHighPtDeDxCalib*> (iter->Next()))) {
    
    calib->Write();
  }
  dataOutFile->Close();
  delete dataOutFile;
  dataOutFile = 0;

  cout << "Nbad (runno == -1) : " << nBad << endl;
}

//___________________________________________________________________________
void AddCalibObjectV0(TList* list, const Char_t* baseName, Bool_t phiCut, 
		      Int_t run, Bool_t analyzeMc)
{
  TString objectName(baseName);
  if(phiCut)
    objectName += "phicut";
  if(run) {
    objectName += "_";
    objectName += run;
  }
  dataOutFile->cd();
  AliHighPtDeDxCalib* data = new AliHighPtDeDxCalib(objectName.Data(), "Calib data");
  if(run) {
    data->SetUseRunCut(kTRUE);
    data->SetRun(run);
  }
  
  list->Add(data);
  data->SetIsMc(analyzeMc);
  data->SetUseFilterCut(kFALSE);
  data->SetUsePhiCut(phiCut);
  if(phiCut) {
    data->SetPhiCutLow(data->GetStandardPhiCutLow());
    data->SetPhiCutHigh(data->GetStandardPhiCutHigh());
  }
  data->Init(nPtBins, xBins);
  data->Print();
}

//___________________________________________________________________________
Int_t CalculateFilter(DeDxTrack* track)
{
  Int_t filter = 0;
  if(track->filterset3)
    filter += 1;
  if(track->filterset2 && !track->filterset3)
    filter += 2;
  if(track->filterset1)
    filter += 4;
  return filter;
}
