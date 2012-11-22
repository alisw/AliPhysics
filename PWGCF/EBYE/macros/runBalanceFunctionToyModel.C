//=========Total multiplicity=========//
Double_t nTotalMultiplicityMean = 1000.;
Double_t nTotalMultiplicitySigma = 10.;

//=========Net charge=========//
Double_t nNetChargeMean = 50.0;
Double_t nNetChargeSigma = 3.0;

//==============Particles and spectra==============//
Double_t gAllChargesTemperature = 0.11; //in GeV
Double_t gPionPercentage = 0.8;
Double_t gPionTemperature = 0.1; //in GeV
Double_t gKaonPercentage = 0.12;
Double_t gKaonTemperature = 0.12; //in GeV
Double_t gProtonPercentage = 0.08;
Double_t gProtonTemperature = 0.2; //in GeV
//==============Particles and spectra==============//

//==============Flow values==============//
Double_t gDirectedFlow = 0.0;
Double_t gEllipticFlow = 0.07;
Double_t gTriangularFlow = 0.0;
Double_t gQuandrangularFlow = 0.0;
Double_t gPentangularFlow = 0.0;
//==============Flow values==============//

//=========Acceptance definition=========//
Double_t gEtaMin = -1.0;
Double_t gEtaMax = 1.0;
Double_t gPtMin = 0.1;
Double_t gPtMax = 100.0;
//=========Acceptance definition=========//

//=========Acceptance filter=========//
Bool_t kUseAcceptanceFilter = kFALSE;
const char *gAcceptanceFilterFile = "efficiencyALICE.root";
//=========Acceptance filter=========//

//=========Dynamical Correlations=========//
Bool_t kUseDynamicalCorrelations = kFALSE;
Double_t gDynamicalCorrelationsPercentage = 0.1;
//=========Dynamical Correlations=========//

Bool_t kUseDebug = kFALSE;

// Run macro used for the toy model analysis
// Author: Panos.Christakoglou@nikhef.nl

//______________________________________________________________________________
void runBalanceFunctionToyModel(Int_t nEvents = 10,
				Bool_t kUseAllCharges = kTRUE) {
  TStopwatch timer;
  timer.Start();

  // load libraries
  gSystem->Load("libCore.so");        
  gSystem->Load("libGeom.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libTree.so");
  gSystem->Load("libSTEERBase.so");
  gSystem->Load("libESD.so");
  gSystem->Load("libAOD.so");
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  gSystem->Load("libEventMixing.so");
  gSystem->Load("libPWGCFebye.so");
  
  //configure the bf objects
  gROOT->LoadMacro("configBalanceFunctionAnalysis.C");
  AliBalance *bf  = GetBalanceFunctionObject("MC");
  AliBalance *bfs = GetBalanceFunctionObject("MC",kTRUE);
  
  //Configure the toy model object
  AliAnalysisTaskToyModel *toyModelAnalysis = new AliAnalysisTaskToyModel();
  if(kUseDebug) toyModelAnalysis->SetDebugFlag();
  toyModelAnalysis->SetAnalysisObject(bf);
  toyModelAnalysis->SetShufflingObject(bfs);
  toyModelAnalysis->SetTotalMultiplicity(nTotalMultiplicityMean,nTotalMultiplicitySigma);
  toyModelAnalysis->SetNetCharge(nNetChargeMean,nNetChargeSigma);
  toyModelAnalysis->SetKinematicsCutsMC(gPtMin,gPtMax,gEtaMin,gEtaMax);

  if(kUseAllCharges) {
    toyModelAnalysis->SetSpectraTemperatureForAllCharges(gAllChargesTemperature);
    toyModelAnalysis->SetDirectedFlowForAllCharges(gDirectedFlow);
    toyModelAnalysis->SetEllipticFlowForAllCharges(gEllipticFlow);
    toyModelAnalysis->SetTriangularFlowForAllCharges(gTriangularFlow);
    toyModelAnalysis->SetQuandrangularFlowForAllCharges(gQuandrangularFlow);
    toyModelAnalysis->SetPentangularFlowForAllCharges(gPentangularFlow);
  }
  else {
    //Pions
    toyModelAnalysis->SetPionPercentage(gPionPercentage);
    toyModelAnalysis->SetSpectraTemperatureForPions(gPionTemperature);
    toyModelAnalysis->SetDirectedFlowForPions(gDirectedFlow);
    toyModelAnalysis->SetEllipticFlowForPions(gEllipticFlow);
    toyModelAnalysis->SetTriangularFlowForPions(gTriangularFlow);
    toyModelAnalysis->SetQuandrangularFlowForPions(gQuandrangularFlow);
    toyModelAnalysis->SetPentangularFlowForPions(gPentangularFlow);

    //Kaons
    toyModelAnalysis->SetKaonPercentage(gKaonPercentage);
    toyModelAnalysis->SetSpectraTemperatureForKaons(gKaonTemperature);
    toyModelAnalysis->SetDirectedFlowForKaons(gDirectedFlow);
    toyModelAnalysis->SetEllipticFlowForKaons(gEllipticFlow);
    toyModelAnalysis->SetTriangularFlowForKaons(gTriangularFlow);
    toyModelAnalysis->SetQuandrangularFlowForKaons(gQuandrangularFlow);
    toyModelAnalysis->SetPentangularFlowForKaons(gPentangularFlow);

    //Protons
    toyModelAnalysis->SetProtonPercentage(gProtonPercentage);
    toyModelAnalysis->SetSpectraTemperatureForProtons(gProtonTemperature);
    toyModelAnalysis->SetDirectedFlowForProtons(gDirectedFlow);
    toyModelAnalysis->SetEllipticFlowForProtons(gEllipticFlow);
    toyModelAnalysis->SetTriangularFlowForProtons(gTriangularFlow);
    toyModelAnalysis->SetQuandrangularFlowForProtons(gQuandrangularFlow);
    toyModelAnalysis->SetPentangularFlowForProtons(gPentangularFlow);
  }

  //Dynamical correlations
  if(kUseDynamicalCorrelations) 
    toyModelAnalysis->SetCorrelationPercentage(gDynamicalCorrelationsPercentage);

  //Acceptance filter
  if(kUseAcceptanceFilter) {
    TFile *gParamFile = TFile::Open(gAcceptanceFilterFile);
    if((!gParamFile) || (!gParamFile->IsOpen())) {
      Printf("File %s not found!!!",acceptanceFilename);
      return;
    }

    TString gParamName;
    for(Int_t iCentrality = 0; iCentrality < numberOfCentralityBins; iCentrality++) {
      gParamName = "gParamCentrality0";//centrality 0-5%
      TF1 *gParameterization = dynamic_cast<TF1 *>(gParamFile->Get(gParamName.Data()));
    }
    toyModelAnalysis->SetAcceptanceParameterization(gParameterization);
  }

  toyModelAnalysis->Init();
  toyModelAnalysis->CreateOutputObjects();
  toyModelAnalysis->Run(nEvents);
  toyModelAnalysis->FinishOutput();

  // Print real and CPU time used for analysis:  
  timer.Stop();
  timer.Print();
}

