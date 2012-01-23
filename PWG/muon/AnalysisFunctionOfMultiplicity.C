/*************************************************************************************************************

This macro analyses the production of SingleMuon and J/Psi as a function of the collision multiplicity in the SPD.
It reads and analyses the output of the Analysis Task PWG3/muon/AliAnalysisTaskMuonCollisionMultiplicity.
Refer to the corresponding files to know what the output of this task is.

Thismacro use a number of input files whose names are hard-coded
These are :
pTRanges.txt : different ranges in pT in which the study is done.
etaRnages.txt : different ranges in eta in which the study is done.


*************************************************************************************************************/




// ROOT includes
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <THnSparse.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TFile.h>
#include <TList.h>
#include <Riostream.h>
#include <TMath.h>

// RooFit includes
#include <RooRealVar.h>
#include <RooAbsReal.h>
#include <RooArgSet.h>
#include <RooCBShape.h>
#include <RooGaussian.h>
#include <RooExponential.h>
#include <RooAddPdf.h>
#include <RooDataHist.h>
#include <RooFitResult.h>
#include <RooPlot.h>

// std includes
#include <vector>




// These functions are in charge of doing all the analysis
void ProduceTriggerGraph(TFile *inputFile, TFile *outputFile, std::vector<Double_t> multiplicityRanges, Bool_t applyZCut, Bool_t applyPileUpCut);

void ProduceSingleMuonRawGraph(TFile *inputFile, TFile *outputFile, std::vector<Double_t> multiplicityRanges, 
			       Bool_t applyZCut, Bool_t applyPileUpCut, Bool_t applyMatchTrigCut, Bool_t applyThetaAbsCut, Bool_t applyPDCACut);

void AnalysisSingleMuon(TFile *outputFile);

void SingleMuonYieldGraph(TFile *outputFile, TGraphErrors *CINT1B);

void SingleMuonYieldOverMeanMultGraph(TFile *outputFile, TGraphErrors *CINT1B);

void SingleMuonYieldNormalisedGraph(TFile *outputFile, TGraphErrors *CMUS1B);

void ProduceFitResults(TFile *inputFile, TFile *outputFile, std::vector<Double_t> multiplicityRanges, 
		       Bool_t applyZCut, Bool_t applyPileUpCut, Double_t numberMatchTrig, Bool_t applyThetaAbsCut, Bool_t applyPDCACut);

void FitInvariantMass(TList *fitList, TH1D *invariantMass, Bool_t areParametersFixed, TH1D *referenceParameters);

void ProduceDimuonGraph(TFile *outputFile, std::vector<Double_t> multiplicityRanges);


///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
void AnalysisFunctionOfMultiplicity(Bool_t doSingleMuonAnalysis = kTRUE, Bool_t doJPsiAnalysis = kTRUE, TString inputName = "./Input/MUON.Multiplicity.THnSparse.ESDs.root", TString multiplicityFileName = "multiplicityJPsi.txt", TString outputName = "./Output/MultiplicityResults.LHC10e.JPsi.ESDs.root")
{
  // The main function of the macro

  // Open the input file and create the output file
  TFile *inputFile = TFile::Open(inputName, "read");
  TFile *outputFile = TFile::Open(outputName, "recreate");

  // First of all, we need to figure which conditions we are using
  // For the event : the cut on z_vertex and if we have pile up from the SPD
  // For the muons : the usual cuts for SMR  
  Bool_t applyZCut = kTRUE;
  Bool_t applyPileUpCut = kTRUE;
  Bool_t applyMatchTrigCut = kTRUE;
  Bool_t applyThetaAbsCut = kTRUE;
  Bool_t applyPDCACut = kTRUE;
  Double_t nMatchTrig = 1.0;


  // Now, we need the multiplicity ranges
  // We get it from an input file
  ifstream multiplicityFile(multiplicityFileName);
  Double_t multiplicityRange = 0.0;
  std::vector <Double_t> multiplicityRanges;
  while(multiplicityFile >> multiplicityRange)
    multiplicityRanges.push_back(multiplicityRange);

  // With that, we can produce the number of CINT1B in each bin
  // Specifficaly, we'll save two TGrpahErrors in the output :
  // - the graph containing the number of minimum bias event in the bin
  // - the graph containing the mean multiplicity in the bin
  ProduceTriggerGraph(inputFile, outputFile, multiplicityRanges, applyZCut, applyPileUpCut);

  // Now, the muon analysis
  if (doSingleMuonAnalysis)
    {
      // First, we produce the single muon raw graph
      // it is the graph giving the number of muons detected in different the multiplicity bins, and in different eta and pT domains
      ProduceSingleMuonRawGraph(inputFile, outputFile, multiplicityRanges, applyZCut, applyPileUpCut, applyMatchTrigCut, applyThetaAbsCut, applyPDCACut);

      // Raw graph are created, it is time to analyse them
      AnalysisSingleMuon(outputFile);
    }


  // Then, the J/Psi analysis
  if (doJPsiAnalysis)
    {
      ProduceFitResults(inputFile, outputFile, multiplicityRanges, applyZCut, applyPileUpCut, nMatchTrig, applyThetaAbsCut, applyPDCACut);
     ProduceDimuonGraph(outputFile, multiplicityRanges);
    }

  inputFile->Close();
  outputFile->Close();
  cout << "End" << endl;
}



///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
void ProduceTriggerGraph(TFile *inputFile, TFile *outputFile, std::vector<Double_t> multiplicityRanges, Bool_t applyZCut, Bool_t applyPileUpCut)
{
  // Compute the number of CINT1B and CMUS1B in each bin
  // Along the x axis, the points are placed at the mean multiplicity of the bin, with an error equal to the error on the mean multiplicity
  
  
  // First of all, we need to get the number of CINT1B as a function of the multiplicity in a TH1D
  // Retrieving the THnSparse in the input
  THnSparseD *inputCINT1B = (THnSparseD *) ((TList *) (inputFile->Get("Trigger;1"))->FindObject("CINT1B") );
  THnSparseD *inputCMUS1B = (THnSparseD *) ((TList *) (inputFile->Get("Trigger;1"))->FindObject("CMUS1B") );


  // Reminder of the architecture of this THnSparse
  // dimension 0 : multiplicity of the event
  // dimension 1 : is the vertex in the z range (0 for no, 1 for yes)?
  // dimension 2 : is it an event without pile up (0 for no, 1 for yes)?

  // Now, applying the cuts on z vertex and pile-up
  if (applyZCut)
    {
      inputCINT1B->GetAxis(1)->SetRangeUser(1.0, 2.0);
      inputCMUS1B->GetAxis(1)->SetRangeUser(1.0, 2.0);
    }
  if (applyPileUpCut)
    {
      inputCINT1B->GetAxis(2)->SetRangeUser(1.0, 2.0);
      inputCMUS1B->GetAxis(2)->SetRangeUser(1.0, 2.0);
    }
  
  // cuts applied, we can project the THnSparse on TH1D, because it is easier to handle
  TH1D *histoCINT1B = inputCINT1B->Projection(0, "E");
  TH1D *histoCMUS1B = inputCMUS1B->Projection(0, "E");
  
  
  // We can now create and fill the two TGraphErrors
  TGraphErrors *graphCINT1B = new TGraphErrors(multiplicityRanges.size() - 1);
  graphCINT1B->SetName("graphCINT1B");
  TGraphErrors *graphCMUS1B = new TGraphErrors(multiplicityRanges.size() - 1);
  graphCMUS1B->SetName("graphCMUS1B");
  
  // Graphs for plotting purposes
  TGraphAsymmErrors *plotCINT1B = new TGraphAsymmErrors(multiplicityRanges.size() - 1);
  plotCINT1B->SetName("plotCINT1B");
  TGraphAsymmErrors *plotCMUS1B = new TGraphAsymmErrors(multiplicityRanges.size() - 1);
  plotCMUS1B->SetName("plotCMUS1B");
  
  for (UInt_t iRange = 0; iRange < multiplicityRanges.size() - 1; iRange++)
    {
      // CINT1B first
      Int_t firstBin = histoCINT1B->GetXaxis()->FindBin(multiplicityRanges[iRange]);
      Int_t lastBin = histoCINT1B->GetXaxis()->FindBin(multiplicityRanges[iRange+1]);
      histoCINT1B->GetXaxis()->SetRange(firstBin, lastBin);

      Double_t total = 0.0;
      Double_t totalErr = 0.0;
      total = histoCINT1B->IntegralAndError(firstBin, lastBin, totalErr);

      Double_t mean = histoCINT1B->GetMean();
      Double_t meanErr = histoCINT1B->GetMeanError();
      
      graphCINT1B->SetPoint(iRange, mean, total);
      graphCINT1B->SetPointError(iRange, meanErr, totalErr);

      plotCINT1B->SetPoint(iRange, mean, total);
      plotCINT1B->SetPointError(iRange, 
				mean - multiplicityRanges[iRange], multiplicityRanges[iRange+1] - mean,
				totalErr, totalErr);

      // CMUS1B second
      firstBin = histoCMUS1B->GetXaxis()->FindBin(multiplicityRanges[iRange]);
      lastBin = histoCMUS1B->GetXaxis()->FindBin(multiplicityRanges[iRange+1]);
      histoCMUS1B->GetXaxis()->SetRange(firstBin, lastBin);

      total = histoCMUS1B->IntegralAndError(firstBin, lastBin, totalErr);

      mean = histoCMUS1B->GetMean();
      meanErr = histoCMUS1B->GetMeanError();
      
      graphCMUS1B->SetPoint(iRange, mean, total);
      graphCMUS1B->SetPointError(iRange, meanErr, totalErr);

      plotCMUS1B->SetPoint(iRange, mean, total);
      plotCMUS1B->SetPointError(iRange, 
				mean - multiplicityRanges[iRange], multiplicityRanges[iRange+1] - mean,
				totalErr, totalErr);
    }

  // Saving the graphs
  outputFile->WriteTObject(plotCINT1B);
  outputFile->WriteTObject(plotCMUS1B);
  outputFile->WriteTObject(graphCINT1B);
  outputFile->WriteTObject(graphCMUS1B);

  delete graphCINT1B;
  delete histoCINT1B;
  delete inputCINT1B;
  delete graphCMUS1B;
  delete histoCMUS1B;
  delete inputCMUS1B;
}


void ProduceSingleMuonRawGraph(TFile *inputFile, TFile *outputFile, std::vector<Double_t> multiplicityRanges, 
			       Bool_t applyZCut, Bool_t applyPileUpCut, Bool_t applyMatchTrigCut, Bool_t applyThetaAbsCut, Bool_t applyPDCACut)
{
  // Produce the raw graph for Single Muon that will be used all along the analysis
  // There are several graph :
  //  - over all the range in both eta and pT
  //  - Single Muons Reference : pT > 1 GeV
  //  - Heavy Flavor Muons : pT > 4 GeV
  //  - bin by bin in pT and over all the range in eta
  //  - bin by bin in eta and with single muon reference
  //  and right now, there is not enough stat to do it bin by bin in both pT and eta
  

  // Retrieve the THnSparse and the Minimum bias histo 
  THnSparseD *inputSingleMuon = (THnSparseD *) ((TList *) (inputFile->Get("SingleMuon;1"))->FindObject("muonCMUS1B") );
  TGraphErrors *CINT1B = (TGraphErrors *) (outputFile->Get("graphCINT1B;1") );

  // Reminder of the architecture of this THnSparse
  // dimension 0  : multiplicity of the event
  // dimension 1  : is the vertex in the z range (0 for no, 1 for yes)?
  // dimension 2  : is it an event without pile up (0 for no, 1 for yes)?
  // dimension 3  : does the muon match the trigger (0 for no, 1 for yes)?
  // dimension 4  : theta_abs of the muon
  // dimension 5  : eta of the muon
  // dimension 6  : p DCA_x of the muon
  // dimension 7  : p DCA_y of the muon
  // dimension 8  : p DCA of the muon
  // dimension 9  : p of the muon
  // dimension 10 : pT of the muon



  // get the pT and eta ranges from files.
  // Beware, the names of the files are hard-coded, so make sure to have them in your folder.
  // I might change this in the future, but I felt the main function had enough parameters already.
  ifstream pTFile ("pTRanges.txt");
  Double_t pTRange = 0.0;
  std::vector <Double_t> pTRanges;
  while(pTFile >> pTRange)
    pTRanges.push_back(pTRange);

  ifstream etaFile ("etaRanges.txt");
  Double_t etaRange = 0.0;
  std::vector <Double_t> etaRanges;
  while(etaFile >> etaRange)
    etaRanges.push_back(etaRange);

 
  // Apply all the cuts
  // Beware, theta_abs and pDCA cuts are hard-coded
  if (applyZCut)
    inputSingleMuon->GetAxis(1)->SetRangeUser(1.0, 2.0);
  if (applyPileUpCut)
    inputSingleMuon->GetAxis(2)->SetRangeUser(1.0, 2.0);
  if (applyMatchTrigCut)
    inputSingleMuon->GetAxis(3)->SetRangeUser(1.0, 2.0);
  if (applyThetaAbsCut)
    inputSingleMuon->GetAxis(4)->SetRangeUser(2.0, 9.0);
  if (applyPDCACut)
    inputSingleMuon->GetAxis(8)->SetRangeUser(0.0, 450.0); // no cut yet, I need to check what are the usual values

  
  // First, we get all the raw histos
  TList *rawSingleMuon = new TList();
  rawSingleMuon->SetName("rawSingleMuon");

  // All integrated
  // There are some hard cuts 
  // - pT : 0.5 -> 8 GeV, because we are not confident on what we are seeing outside of these limits
  // - eta : -4.0, -> -2.5, acceptance of the spectrometer

  Double_t etaMinAbsolute = -4.0;
  Double_t etaMaxAbsolute = -2.5;
  inputSingleMuon->GetAxis(5)->SetRangeUser(etaMinAbsolute, etaMaxAbsolute);

  Double_t pTMinAbsolute = 0.5;
  Double_t pTMaxAbsolute = 8.0;
  inputSingleMuon->GetAxis(10)->SetRangeUser(pTMinAbsolute, pTMaxAbsolute);
  
  // Putting this in a TH3D, since it is easier to use and allow for computation of errors on the integral
  // TH3D is multiplicity, eta and pT
  TH3D *histoSingleMuon = inputSingleMuon->Projection(0, 5, 10, "E");

  // Declare all the histos
  TGraphAsymmErrors *allSingleMuon = new TGraphAsymmErrors(multiplicityRanges.size() - 1);
  allSingleMuon->SetName("allSingleMuon");
  TGraphAsymmErrors *rawSMR = new TGraphAsymmErrors(multiplicityRanges.size() - 1); // Single Muon Reference (pT > 1.0 GeV)
  rawSMR->SetName("rawSMR");
  TGraphAsymmErrors *rawHFM = new TGraphAsymmErrors(multiplicityRanges.size() - 1); // Heavy Flavor Muon (pT > 4.0 GeV)
  rawHFM->SetName("rawHFM");


  for (UInt_t iRange = 0; iRange < multiplicityRanges.size() - 1; iRange++)
    {
      Int_t firstBinX = histoSingleMuon->GetXaxis()->FindBin(multiplicityRanges[iRange]);
      Int_t lastBinX = histoSingleMuon->GetXaxis()->FindBin(multiplicityRanges[iRange+1]);
      Int_t firstBinY = histoSingleMuon->GetYaxis()->GetFirst();
      Int_t lastBinY = histoSingleMuon->GetYaxis()->GetLast();
      Int_t firstBinZ = histoSingleMuon->GetZaxis()->GetFirst();
      Int_t lastBinZ = histoSingleMuon->GetZaxis()->GetLast();

      Double_t nSingleMuon = 0.0;
      Double_t errSingleMuon = 0.0;
      nSingleMuon = histoSingleMuon->IntegralAndError(firstBinX, lastBinX, firstBinY, lastBinY, firstBinZ, lastBinZ, errSingleMuon);
      
      allSingleMuon->SetPoint(iRange, CINT1B->GetX()[iRange], nSingleMuon);
      allSingleMuon->SetPointError(iRange, 
				   CINT1B->GetX()[iRange] - multiplicityRanges[iRange],
				   multiplicityRanges[iRange+1] - CINT1B->GetX()[iRange],
				   errSingleMuon,
				   errSingleMuon);

      // SMR
      firstBinZ = histoSingleMuon->GetZaxis()->FindBin(1.0);
      lastBinZ = histoSingleMuon->GetZaxis()->GetLast();
      nSingleMuon = histoSingleMuon->IntegralAndError(firstBinX, lastBinX, firstBinY, lastBinY, firstBinZ, lastBinZ, errSingleMuon);
      
      rawSMR->SetPoint(iRange, CINT1B->GetX()[iRange], nSingleMuon);
      rawSMR->SetPointError(iRange, 
			    CINT1B->GetX()[iRange] - multiplicityRanges[iRange],
			    multiplicityRanges[iRange+1] - CINT1B->GetX()[iRange],
			    errSingleMuon,
			    errSingleMuon);

      // HFM
      firstBinZ = histoSingleMuon->GetZaxis()->FindBin(4.0);
      lastBinZ = histoSingleMuon->GetZaxis()->GetLast();
      nSingleMuon = histoSingleMuon->IntegralAndError(firstBinX, lastBinX, firstBinY, lastBinY, firstBinZ, lastBinZ, errSingleMuon);
      
      rawHFM->SetPoint(iRange, CINT1B->GetX()[iRange], nSingleMuon);
      rawHFM->SetPointError(iRange, 
			    CINT1B->GetX()[iRange] - multiplicityRanges[iRange],
			    multiplicityRanges[iRange+1] - CINT1B->GetX()[iRange],
			    errSingleMuon,
			    errSingleMuon);
    }

  TList *integratedSingleMuon = new TList();
  integratedSingleMuon->SetName("rawSingleMuonIntegrated");

  integratedSingleMuon->AddAt(allSingleMuon, 0);
  integratedSingleMuon->AddAt(rawSMR, 1);
  integratedSingleMuon->AddAt(rawHFM, 2);

  rawSingleMuon->AddAt(integratedSingleMuon, 0);
  //delete allSingleMuon;
  //delete rawSMR;
  //delete rawHFM;

  // Now, for the study in pT
  TList *rawSingleMuonPt = new TList();
  rawSingleMuonPt->SetName("rawSingleMuonPt");

  for (UInt_t ipT = 0; ipT < pTRanges.size()-1; ipT++)
    {
      TString nameGraph;
      nameGraph.Form("pTRange_%1f-%1f", pTRanges[ipT], pTRanges[ipT+1]);
      TGraphAsymmErrors *singleMuonPtRange = new TGraphAsymmErrors();
      singleMuonPtRange->SetName(nameGraph);

      for (UInt_t iRange = 0; iRange < multiplicityRanges.size() - 1; iRange++)
	{
	  Int_t firstBinX = histoSingleMuon->GetXaxis()->FindBin(multiplicityRanges[iRange]);
	  Int_t lastBinX = histoSingleMuon->GetXaxis()->FindBin(multiplicityRanges[iRange+1]);
	  Int_t firstBinY = histoSingleMuon->GetYaxis()->GetFirst();
	  Int_t lastBinY = histoSingleMuon->GetYaxis()->GetLast();
	  Int_t firstBinZ = histoSingleMuon->GetZaxis()->FindBin(pTRanges[ipT]);
	  Int_t lastBinZ = histoSingleMuon->GetZaxis()->FindBin(pTRanges[ipT]+1);

	  Double_t nSingleMuon = 0.0;
	  Double_t errSingleMuon = 0.0;
	  nSingleMuon = histoSingleMuon->IntegralAndError(firstBinX, lastBinX, firstBinY, lastBinY, firstBinZ, lastBinZ, errSingleMuon);
	  
	  singleMuonPtRange->SetPoint(iRange, CINT1B->GetX()[iRange], nSingleMuon);
	  singleMuonPtRange->SetPointError(iRange, 
					   CINT1B->GetX()[iRange] - multiplicityRanges[iRange],
					   multiplicityRanges[iRange+1] - CINT1B->GetX()[iRange],
					   errSingleMuon,
					   errSingleMuon);
	  
	}

      rawSingleMuonPt->AddAt(singleMuonPtRange, ipT);
      //delete singleMuonPtRange;
    }

  rawSingleMuon->AddAt(rawSingleMuonPt, 1);
 
  // At last, the study in eta
  // Both on SMR and HFM
  TList *rawSingleMuonEta = new TList();
  rawSingleMuonEta->SetName("rawSingleMuonEta");

  for (UInt_t iEta = 0; iEta < etaRanges.size()-1; iEta++)
    {
      TString nameGraph;
      nameGraph.Form("SMRetaRange_%1f-%1f", etaRanges[iEta], etaRanges[iEta+1]);
      TGraphAsymmErrors *SMREtaRange = new TGraphAsymmErrors();
      SMREtaRange->SetName(nameGraph);

      nameGraph.Form("HFMetaRange_%1f-%1f", etaRanges[iEta], etaRanges[iEta+1]);
      TGraphAsymmErrors *HFMEtaRange = new TGraphAsymmErrors();
      HFMEtaRange->SetName(nameGraph);

      for (UInt_t iRange = 0; iRange < multiplicityRanges.size() - 1; iRange++)
	{
	  // First, SMR
	  Int_t firstBinX = histoSingleMuon->GetXaxis()->FindBin(multiplicityRanges[iRange]);
	  Int_t lastBinX = histoSingleMuon->GetXaxis()->FindBin(multiplicityRanges[iRange+1]);
	  Int_t firstBinY = histoSingleMuon->GetYaxis()->FindBin(etaRanges[iEta]);
	  Int_t lastBinY = histoSingleMuon->GetYaxis()->FindBin(etaRanges[iEta]+1);
	  Int_t firstBinZ = histoSingleMuon->GetZaxis()->FindBin(1.0);
	  Int_t lastBinZ = histoSingleMuon->GetZaxis()->GetLast();

	  Double_t nSingleMuon = 0.0;
	  Double_t errSingleMuon = 0.0;
	  nSingleMuon = histoSingleMuon->IntegralAndError(firstBinX, lastBinX, firstBinY, lastBinY, firstBinZ, lastBinZ, errSingleMuon);
	  
	  SMREtaRange->SetPoint(iRange, CINT1B->GetX()[iRange], nSingleMuon);
	  SMREtaRange->SetPointError(iRange, 
				     CINT1B->GetX()[iRange] - multiplicityRanges[iRange],
				     multiplicityRanges[iRange+1] - CINT1B->GetX()[iRange],
				     errSingleMuon,
				     errSingleMuon);

	  // Then, HFM
	  firstBinZ = histoSingleMuon->GetZaxis()->FindBin(4.0);
	  lastBinZ = histoSingleMuon->GetZaxis()->GetLast();

	  nSingleMuon = histoSingleMuon->IntegralAndError(firstBinX, lastBinX, firstBinY, lastBinY, firstBinZ, lastBinZ, errSingleMuon);
	  
	  HFMEtaRange->SetPoint(iRange, CINT1B->GetX()[iRange], nSingleMuon);
	  HFMEtaRange->SetPointError(iRange, 
				     CINT1B->GetX()[iRange] - multiplicityRanges[iRange],
				     multiplicityRanges[iRange+1] - CINT1B->GetX()[iRange],
				     errSingleMuon,
				     errSingleMuon);

	  
	}

      rawSingleMuonEta->AddAt(SMREtaRange, iEta);
      rawSingleMuonEta->AddAt(HFMEtaRange, iEta + etaRanges.size()-1);
      //delete SMREtaRange;
      //delete HFMEtaRange;
    }

  rawSingleMuon->AddAt(rawSingleMuonEta, 2);

  outputFile->WriteTObject(rawSingleMuon);
  delete integratedSingleMuon;
  delete rawSingleMuonPt;
  delete rawSingleMuonEta;
}


void AnalysisSingleMuon(TFile *outputFile)
{
  // Analyse the raw data for single muons

  // First, we retrieve the CINT1B graph, as usual
  TGraphErrors *CINT1B = (TGraphErrors *) outputFile->Get("graphCINT1B;1");
  TGraphErrors *CMUS1B = (TGraphErrors *) outputFile->Get("graphCMUS1B;1");

  // Creation of all the yield graphs
  SingleMuonYieldGraph(outputFile, CINT1B);

  // Creation of the yield over mean mult graph
  SingleMuonYieldOverMeanMultGraph(outputFile, CINT1B);

  // Creation of the yield with the reference bin normalised to 1, for comparison
  SingleMuonYieldNormalisedGraph(outputFile, CMUS1B);
}




//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
void SingleMuonYieldGraph(TFile *outputFile, TGraphErrors *CINT1B)
{
  // This function use the raw graphs of outputfile and the CINT1B graph to create the yield graph
  // yield are raw divided by CINT1B
  TList *yieldSingleMuon = new TList();
  yieldSingleMuon->SetName("yieldSingleMuon");

  for (Int_t iList = 0; iList < ((TList *) outputFile->Get("rawSingleMuon;1"))->GetEntries(); iList++)
    {
      TList *currentList = (TList *) (( TList *) outputFile->Get("rawSingleMuon;1"))->At(iList);
      TList *newList = new TList();
      TString listName = currentList->GetName();
      listName.ReplaceAll("raw", 3, "yield", 5);
      newList->SetName(listName);

      for (Int_t iGraph = 0; iGraph < currentList->GetEntries(); iGraph++)
	{
	  TGraphAsymmErrors *rawGraph = (TGraphAsymmErrors *) currentList->At(iGraph);

	  TString yieldName = rawGraph->GetName();
	  yieldName.ReplaceAll("raw", 3, "yield", 5);
	  TGraphAsymmErrors *yieldGraph = new TGraphAsymmErrors(rawGraph->GetN());
	  yieldGraph->SetName(yieldName);

	  // fill the yield graph
	  for (Int_t iBin = 0; iBin < yieldGraph->GetN(); iBin++ )
	    {
	      if (CINT1B->GetY()[iBin] != 0.0)
		{
		  yieldGraph->SetPoint(iBin, 
				       rawGraph->GetX()[iBin], 
				       rawGraph->GetY()[iBin]/CINT1B->GetY()[iBin]);
		  if (rawGraph->GetY()[iBin] != 0.0)
		    { 
		      Double_t error = yieldGraph->GetY()[iBin] * 
			TMath::Sqrt((rawGraph->GetEYlow()[iBin]/rawGraph->GetY()[iBin])*(rawGraph->GetEYlow()[iBin]/rawGraph->GetY()[iBin]) + 
				    (CINT1B->GetEY()[iBin]/CINT1B->GetY()[iBin])*(CINT1B->GetEY()[iBin]/CINT1B->GetY()[iBin]));			
		      yieldGraph->SetPointError(iBin, 
					      rawGraph->GetEXlow()[iBin],
					      rawGraph->GetEXhigh()[iBin],
					      error, error); 
		    }
		}
	      
	      else
		yieldGraph->SetPoint(iBin, 0, 0);
	    }

	  newList->AddAt(yieldGraph, iGraph);
	}
      
      yieldSingleMuon->AddAt(newList, iList);
    }
  
  outputFile->WriteTObject(yieldSingleMuon);
}



//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
void SingleMuonYieldOverMeanMultGraph(TFile *outputFile, TGraphErrors *CINT1B)
{
  // This function use the yield graphs of outputfile to create the yield graph divided by the mean multiplicity in each bin
  // The CINT1B graph is necessary to get the error on the mean multiplicity
  TList *yieldOverMeanMultSingleMuon = new TList();
  yieldOverMeanMultSingleMuon->SetName("yieldOverMeanMultSingleMuon");

  for (Int_t iList = 0; iList < ((TList *) outputFile->Get("yieldSingleMuon;1"))->GetEntries(); iList++)
    {
      TList *currentList = (TList *) (( TList *) outputFile->Get("yieldSingleMuon;1"))->At(iList);
      TList *newList = new TList();
      TString listName = currentList->GetName();
      listName.ReplaceAll("yield", 5, "yieldOverMeanMult", 17);
      newList->SetName(listName);

      for (Int_t iGraph = 0; iGraph < currentList->GetEntries(); iGraph++)
	{
	  TGraphAsymmErrors *yieldGraph = (TGraphAsymmErrors *) currentList->At(iGraph);

	  TString yieldOverMeanMultName = yieldGraph->GetName();
	  yieldOverMeanMultName.ReplaceAll("yield", 5, "yieldOverMeanMult", 17);
	  TGraphAsymmErrors *yieldOverMeanMultGraph = new TGraphAsymmErrors(yieldGraph->GetN());
	  yieldOverMeanMultGraph->SetName(yieldOverMeanMultName);

	  // fill the yield graph
	  for (Int_t iBin = 0; iBin < yieldOverMeanMultGraph->GetN(); iBin++ )
	    {
	      if (CINT1B->GetX()[iBin] != 0.0)
		{
		  yieldOverMeanMultGraph->SetPoint(iBin, 
						   yieldGraph->GetX()[iBin], 
						   yieldGraph->GetY()[iBin]/CINT1B->GetX()[iBin]);
		  if (yieldGraph->GetY()[iBin] != 0.0)
		    {
		      Double_t error = yieldOverMeanMultGraph->GetY()[iBin]*
			TMath::Sqrt((yieldGraph->GetEYlow()[iBin]/yieldGraph->GetY()[iBin])*(yieldGraph->GetEYlow()[iBin]/yieldGraph->GetY()[iBin]) +
				    (CINT1B->GetEX()[iBin]/CINT1B->GetX()[iBin])*(CINT1B->GetEX()[iBin]/CINT1B->GetX()[iBin]));
		      yieldOverMeanMultGraph->SetPointError(iBin, 
							    yieldGraph->GetEXlow()[iBin],
							    yieldGraph->GetEXhigh()[iBin],
							    error, error); 
		    }
		}
	      
	      else
		yieldOverMeanMultGraph->SetPoint(iBin, 0, 0);
	    }
	  
	  newList->AddAt(yieldOverMeanMultGraph, iGraph);
	}
      
      yieldOverMeanMultSingleMuon->AddAt(newList, iList);
    }
  
  outputFile->WriteTObject(yieldOverMeanMultSingleMuon);
}


//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
void SingleMuonYieldNormalisedGraph(TFile *outputFile, TGraphErrors *CMUS1B)
{
  // This function normalise the yield over mean multiplicity so that the reference bin in multiplicity is equal to 1.
  // The reference bin is defined it as the bin with the maximum number of CMUS1B
  TList *yieldNormalisedSingleMuon = new TList();
  yieldNormalisedSingleMuon->SetName("yieldNormalisedSingleMuon");


  // The first step is to find the reference bin in multiplicity
  // we define it as the bin with the maximum number of CMUS1B
  
  Double_t maxCMUS1B = 0.0;
  Int_t referenceBin = 0;
  for (Int_t iBin = 0; iBin < CMUS1B->GetN(); iBin++)
    if (CMUS1B->GetX()[iBin] > 10.0)
      if (CMUS1B->GetY()[iBin] > maxCMUS1B)
	{
	  maxCMUS1B = CMUS1B->GetY()[iBin];
	  referenceBin = iBin;
	}

  // Now the loop, as usual
  // We don't propagate the error on the normalisation factor, since this is artificial
  for (Int_t iList = 0; iList < ((TList *) outputFile->Get("yieldOverMeanMultSingleMuon;1"))->GetEntries(); iList++)
    {
      TList *currentList = (TList *) (( TList *) outputFile->Get("yieldOverMeanMultSingleMuon;1"))->At(iList);
      TList *newList = new TList();
      TString listName = currentList->GetName();
      listName.ReplaceAll("yieldOverMeanMult", 17, "yieldNormalised", 15);
      newList->SetName(listName);

      for (Int_t iGraph = 0; iGraph < currentList->GetEntries(); iGraph++)
	{
	  TGraphAsymmErrors *yieldOverMeanMultGraph = (TGraphAsymmErrors *) currentList->At(iGraph);

	  TString yieldNormalisedName = yieldOverMeanMultGraph->GetName();
	  yieldNormalisedName.ReplaceAll("yieldOverMeanMult", 17, "yieldNormalised", 15);
	  TGraphAsymmErrors *yieldNormalisedGraph = new TGraphAsymmErrors(yieldOverMeanMultGraph->GetN());
	  yieldNormalisedGraph->SetName(yieldNormalisedName);

	  // fill the yield graph
	  for (Int_t iBin = 0; iBin < yieldNormalisedGraph->GetN(); iBin++ )
	    {
	      if (yieldOverMeanMultGraph->GetY()[iBin] != 0.0)
		{
		  yieldNormalisedGraph->SetPoint(iBin, 
						 yieldOverMeanMultGraph->GetX()[iBin], 
						 yieldOverMeanMultGraph->GetY()[iBin]/yieldOverMeanMultGraph->GetY()[referenceBin]);
		  if (yieldOverMeanMultGraph->GetY()[iBin] != 0.0)
		    yieldNormalisedGraph->SetPointError(iBin,
							yieldOverMeanMultGraph->GetEXlow()[iBin],
							yieldOverMeanMultGraph->GetEXhigh()[iBin],
							yieldOverMeanMultGraph->GetEYlow()[iBin]/yieldOverMeanMultGraph->GetY()[referenceBin],
							yieldOverMeanMultGraph->GetEYhigh()[iBin]/yieldOverMeanMultGraph->GetY()[referenceBin]); 
		}
	      
	      else
		yieldNormalisedGraph->SetPoint(iBin, 0, 0);
	    }

	  newList->AddAt(yieldNormalisedGraph, iGraph);
	}
      
      yieldNormalisedSingleMuon->AddAt(newList, iList);
    }
  
  outputFile->WriteTObject(yieldNormalisedSingleMuon);
}


/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
void ProduceFitResults(TFile *inputFile, TFile *outputFile, std::vector<Double_t> multiplicityRanges, 
		       Bool_t applyZCut, Bool_t applyPileUpCut, Double_t numberMatchTrig, Bool_t applyThetaAbsCut, Bool_t applyPDCACut)
{
  // Produce the raw graph for Dimuons that will be used all along the analysis
  // There are two graphs :
  //  - J/Psi
  //  - Background in the J/Psi range

  // Retrieve the THnSparse and the Minimum bias histo 
  THnSparseD *inputDimuon = (THnSparseD *) ((TList *) (inputFile->Get("Dimuon;1"))->FindObject("dimuonCMUS1B") );

  // Reminder of the structure of the dimuon histo
  // dimension 0  : multiplicity of the event
  // dimension 1  : is the vertex in the z range (0 for no, 1 for yes)?
  // dimension 2  : is it an event without pile up (0 for no, 1 for yes)?
  // dimension 3  : does the first muon match the trigger (0 for no, 1 for yes)?
  // dimension 4  : does the second muon match the trigger (0 for no, 1 for yes)?
  // dimension 5  : number of muons matching the trigger in the dimuon
  // dimension 6  : theta_abs of the first muon
  // dimension 7  : theta_abs of the second muon
  // dimension 8  : eta of the first muon
  // dimension 9  : eta of the second muon
  // dimension 10 : p DCA_x of the first muon
  // dimension 11 : p DCA_y of the first muon
  // dimension 12 : p DCA_x of the second muon
  // dimension 13 : p DCA_y of the second muon
  // dimension 14 : p of the first muon
  // dimension 15 : p of the second muon
  // dimension 16 : p of the dimuon
  // dimension 17 : pT of the dimuon
  // dimension 18 : invariant mass of the dimuon



  // get the pT and eta ranges from files.
  // Beware, the names of the files are hard-coded, so make sure to have them in your folder.
  // I might change this in the future, but I felt the main function had enough parameters already.
  ifstream pTFile ("pTRanges.txt");
  Double_t pTRange = 0.0;
  std::vector <Double_t> pTRanges;
  while(pTFile >> pTRange)
    pTRanges.push_back(pTRange);

  ifstream etaFile ("etaRanges.txt");
  Double_t etaRange = 0.0;
  std::vector <Double_t> etaRanges;
  while(etaFile >> etaRange)
    etaRanges.push_back(etaRange);

 
  // Apply all the cuts
  // Beware, theta_abs and pDCA cuts are hard-coded
  if (applyZCut)
    inputDimuon->GetAxis(1)->SetRangeUser(1.0, 2.0);
  if (applyPileUpCut)
    inputDimuon->GetAxis(2)->SetRangeUser(1.0, 2.0);
  inputDimuon->GetAxis(5)->SetRangeUser(numberMatchTrig, 3.0);
  if (applyThetaAbsCut)
    {
      inputDimuon->GetAxis(6)->SetRangeUser(2.0, 9.0);
      inputDimuon->GetAxis(7)->SetRangeUser(2.0, 9.0);
    }
  if (applyPDCACut)
    {// no cut yet, I need to check what are the usual values
      inputDimuon->GetAxis(12)->SetRangeUser(0.0, 450.0); 
      inputDimuon->GetAxis(15)->SetRangeUser(0.0, 450.0); 
    }


  // First, we get the invariant mass histos
  TList *dimuonFitResults = new TList();
  dimuonFitResults->SetName("dimuonFitResults");

  // There are some hard cuts 
  // - eta : -4.0, -> -2.5, acceptance of the spectrometer

  Double_t etaMinAbsolute = -4.0;
  Double_t etaMaxAbsolute = -2.5;
  inputDimuon->GetAxis(8)->SetRangeUser(etaMinAbsolute, etaMaxAbsolute);
  inputDimuon->GetAxis(9)->SetRangeUser(etaMinAbsolute, etaMaxAbsolute);

  // First, we get the invariant mass histo for the whole range in multiplicity
  inputDimuon->GetAxis(0)->SetRangeUser(multiplicityRanges[0], multiplicityRanges[multiplicityRanges.size()-1]);
  TH1D *invariantMassIntegrated = inputDimuon->Projection(18, "E");
  invariantMassIntegrated->SetName("invariantMassIntegrated");


  // Now for the fit on all the multiplicity
  // This is used to get a value for some parameters
  TList *fitAll = new TList();
  fitAll->SetName("AllMultiplicity");


  // Create the container for the reference parameters
  TH1D *referenceParameters = new TH1D("referenceParameters", "parameters", 8, 1.0, 9.0);
  referenceParameters->Sumw2();
  referenceParameters->GetXaxis()->SetBinLabel(1, "meanJPsi");
  referenceParameters->GetXaxis()->SetBinLabel(2, "sigmaJPsi");
  referenceParameters->GetXaxis()->SetBinLabel(3, "alphaJPsi");
  referenceParameters->GetXaxis()->SetBinLabel(4, "nJPsi");
  referenceParameters->GetXaxis()->SetBinLabel(5, "meanPsiPrime");
  referenceParameters->GetXaxis()->SetBinLabel(6, "sigmaPsiPrime");
  referenceParameters->GetXaxis()->SetBinLabel(7, "expBkg1");
  referenceParameters->GetXaxis()->SetBinLabel(8, "expBkg2");
  

  FitInvariantMass(fitAll, invariantMassIntegrated, kFALSE, referenceParameters);
  dimuonFitResults->AddAt(fitAll, 0);

  // Now, we make a fit for each range in multiplicity
  for (UInt_t iMult = 0; iMult < multiplicityRanges.size()-1; iMult++)
    {
      TString name;
      name.Form("Multiplicity_%d-%d", static_cast<Int_t>(multiplicityRanges[iMult]), static_cast<Int_t>(multiplicityRanges[iMult+1]));
      TList *fitRange = new TList();
      fitRange->SetName(name);

      inputDimuon->GetAxis(0)->SetRangeUser(multiplicityRanges[iMult], multiplicityRanges[iMult+1]);
      TH1D *invariantMassRange = inputDimuon->Projection(18, "E");
      invariantMassRange->SetName("invariantMass");

      FitInvariantMass(fitRange, invariantMassRange, kTRUE, referenceParameters);
      dimuonFitResults->AddAt(fitRange, iMult + 1);
    }

  outputFile->WriteTObject(dimuonFitResults);
}




/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void FitInvariantMass(TList *fitList, TH1D *invariantMass, Bool_t areParametersFixed, TH1D *referenceParameters)
{
  // Fit the given histo and put the results in outputFile
  // The fit is : Crystal Ball (J/Psi) + gaussian (Psi') + two exponentials (background)

  // The parameters histo architecture is the following :
  // bin 1 : mean J/Psi
  // bin 2 : sigma J/Psi
  // bin 3 : alpha J/Psi
  // bin 4 : n J/Psi
  // bin 5 : mean Psi'
  // bin 6 : sigma Psi'
  // bin 7 : exponential factor background 1
  // bin 8 : exponential factor background 2

  // The results histo architecture is the following :
  // bin 1 : J/Psi
  // bin 2 : Background J/Psi at 2 sigma
  // bin 3 : signal over background
  // bin 4 : chi2

  // There are some parameters that are not always free :
  // mean of J/Psi peak
  // sigma of J/Psi peak
  // alpha of J/Psi function (Crystall Ball function)
  // n of J/Psi function
  // mean of Psi' peak
  // sigma of Psi' peak
  // If there are fixed, their values are taken from the parameters histo


  // First, we have to create two histos, one that will contain the parameters of the fit, and the other the results
  TH1D *parameters = new TH1D("parameters", "parameters", 8, 1.0, 9.0);
  parameters->Sumw2();
  parameters->GetXaxis()->SetBinLabel(1, "meanJPsi");
  parameters->GetXaxis()->SetBinLabel(2, "sigmaJPsi");
  parameters->GetXaxis()->SetBinLabel(3, "alphaJPsi");
  parameters->GetXaxis()->SetBinLabel(4, "nJPsi");
  parameters->GetXaxis()->SetBinLabel(5, "meanPsiPrime");
  parameters->GetXaxis()->SetBinLabel(6, "sigmaPsiPrime");
  parameters->GetXaxis()->SetBinLabel(7, "expBkg1");
  parameters->GetXaxis()->SetBinLabel(8, "expBkg2");
  TH1D *results = new TH1D("results", "results", 4, 1.0, 5.0);
  results->Sumw2();
  results->GetXaxis()->SetBinLabel(1, "JPsi");
  results->GetXaxis()->SetBinLabel(2, "bkg");
  results->GetXaxis()->SetBinLabel(3, "SB");
  results->GetXaxis()->SetBinLabel(4, "chi2");


  // Declare observable M
  RooRealVar M("M","Dimuon invariant mass [GeV/c^2]", 2.0, 5.0);


  // Declare Cristal Ball for J/Psi
  RooRealVar mean_jpsi("mean_jpsi","mean_jpsi", 3.096, 2.8, 3.2);
  RooRealVar sigma_jpsi("sigma_jpsi","sigma_jpsi", 0.070, 0, 0.2);
  RooRealVar alpha_jpsi("alpha","alpha", 5, 0, 10);
  RooRealVar n_jpsi("n","n",5,0,10);
  RooCBShape jPsi("jpsi", "crystal ball PDF", M, mean_jpsi, sigma_jpsi, alpha_jpsi, n_jpsi);
  

  //Declare gaussian for Psi prime
  RooRealVar mean_psip("mean_psip", "mean_psip", 3.686, 3.6, 3.75);
  RooRealVar sigma_psip("sigma_psip","sigma_psip", 0.086);
  RooGaussian psiPrime("psip","Gaussian", M, mean_psip, sigma_psip);
  
  //Declare exponentials for background
  RooRealVar c1("c1", "c1", -10, -0.1);
  RooExponential bkg_exp1("bkg_exp1", "background", M, c1);
  
  RooRealVar c2("c2", "c2", -5.0, -1);
  RooExponential bkg_exp2("bkg_exp2", "background", M, c2);

  // Give a value and fix sigma, alpha and n in the Cristal Ball
  // The values correspond to a previous fit (ideally, with all the statistic)
  // I should automatize this one of these day
  if (areParametersFixed)
    {
      Double_t fixMeanJPsi = referenceParameters->GetBinContent(parameters->GetXaxis()->FindBin("meanJPsi"));
      mean_jpsi.setVal(fixMeanJPsi);
      mean_jpsi.setError(0.0);
      mean_jpsi.setRange(fixMeanJPsi, fixMeanJPsi);

      Double_t fixSigmaJPsi = referenceParameters->GetBinContent(parameters->GetXaxis()->FindBin("sigmaJPsi"));
      sigma_jpsi.setVal(fixSigmaJPsi);
      sigma_jpsi.setError(0.0);
      sigma_jpsi.setRange(fixSigmaJPsi, fixSigmaJPsi);

      Double_t fixAlphaJPsi = referenceParameters->GetBinContent(parameters->GetXaxis()->FindBin("alphaJPsi"));
      alpha_jpsi.setVal(fixAlphaJPsi);
      alpha_jpsi.setError(0.0);
      alpha_jpsi.setRange(fixAlphaJPsi, fixAlphaJPsi);

      Double_t fixNJPsi = referenceParameters->GetBinContent(parameters->GetXaxis()->FindBin("nJPsi"));
      n_jpsi.setVal(fixNJPsi);
      n_jpsi.setError(0.0);
      n_jpsi.setRange(fixNJPsi, fixNJPsi);

      Double_t fixMeanPsiPrime = referenceParameters->GetBinContent(parameters->GetXaxis()->FindBin("meanPsiPrime"));
      mean_psip.setVal(fixMeanPsiPrime);
      mean_psip.setError(0.0);
      mean_psip.setRange(fixMeanPsiPrime, fixMeanPsiPrime);

      Double_t fixSigmaPsiPrime = referenceParameters->GetBinContent(parameters->GetXaxis()->FindBin("sigmaPsiPrime"));
      sigma_psip.setVal(fixSigmaPsiPrime);
      sigma_psip.setError(0.0);
      sigma_psip.setRange(fixSigmaPsiPrime, fixSigmaPsiPrime);
    }


  // Sum the composite signal and background into an extended pdf nsig*sig+nbkg*bkg
  RooRealVar fitJPsi("fitJPsi", "number of J/Psi events", 1600 ,0.0, 15000);
  RooRealVar fitPsiPrime("fitPsiPrime", "number of Psi Prime events", 50, 0.0, 200);
  RooRealVar fitBkg1("fitBkg1", "number of background events", 5000, 0.0, 15000);
  RooRealVar fitBkg2("fitBkg2", "number of background events", 5000, 0.0, 15000);
  

  RooAddPdf *fitFunction = new RooAddPdf("model", "(CB+Gauss+2exp)", RooArgList(bkg_exp1, bkg_exp2, jPsi, psiPrime), RooArgList(fitBkg1, fitBkg2, fitJPsi, fitPsiPrime));

  // Define the histo to be fitted
  RooDataHist *invariantMassFit = new RooDataHist("invariantMassFit", "invariantMassFit", M, invariantMass);

  // Fit the invariant mass
  // This is the important part
  RooFitResult *fitResult = fitFunction->fitTo(*invariantMassFit, RooFit::Save());
  //RooFitResult *fitResult = fitFunction->chi2FitTo(*invariantMassFit, RooFit::Save());

  RooPlot *plot = M.frame();
  TString title = "fittedInvariantMass";
  plot->SetTitle(title);
  plot->SetName(title);

  invariantMassFit->plotOn(plot, RooFit::Name("invariantMassFit"));
  fitFunction->plotOn(plot, RooFit::Name("fitFunction"), RooFit::Range(2.0, 5.0));


  // Fill the histo with the final values of the parameters
  parameters->SetBinContent(parameters->GetXaxis()->FindBin("meanJPsi"), mean_jpsi.getVal());
  parameters->SetBinError(parameters->GetXaxis()->FindBin("meanJPsi"), mean_jpsi.getError());
  
  parameters->SetBinContent(parameters->GetXaxis()->FindBin("sigmaJPsi"), sigma_jpsi.getVal());
  parameters->SetBinError(parameters->GetXaxis()->FindBin("sigmaJPsi"), sigma_jpsi.getError());
  
  parameters->SetBinContent(parameters->GetXaxis()->FindBin("alphaJPsi"), alpha_jpsi.getVal());
  parameters->SetBinError(parameters->GetXaxis()->FindBin("alphaJPsi"), alpha_jpsi.getError());
  
  parameters->SetBinContent(parameters->GetXaxis()->FindBin("nJPsi"), n_jpsi.getVal());
  parameters->SetBinError(parameters->GetXaxis()->FindBin("nJPsi"), n_jpsi.getError());
  
  parameters->SetBinContent(parameters->GetXaxis()->FindBin("meanPsiPrime"), mean_psip.getVal());
  parameters->SetBinError(parameters->GetXaxis()->FindBin("meanPsiPrime"), mean_psip.getError());

  parameters->SetBinContent(parameters->GetXaxis()->FindBin("sigmaPsiPrime"), sigma_psip.getVal());
  parameters->SetBinError(parameters->GetXaxis()->FindBin("sigmaPsiPrime"), sigma_psip.getError());
  
  // Fill the histo with the reference parameters
  if (!areParametersFixed)
    {
      referenceParameters->SetBinContent(parameters->GetXaxis()->FindBin("meanJPsi"), mean_jpsi.getVal());
      referenceParameters->SetBinError(parameters->GetXaxis()->FindBin("meanJPsi"), mean_jpsi.getError());

      referenceParameters->SetBinContent(parameters->GetXaxis()->FindBin("sigmaJPsi"), sigma_jpsi.getVal());
      referenceParameters->SetBinError(parameters->GetXaxis()->FindBin("sigmaJPsi"), sigma_jpsi.getError());

      referenceParameters->SetBinContent(parameters->GetXaxis()->FindBin("alphaJPsi"), alpha_jpsi.getVal());
      referenceParameters->SetBinError(parameters->GetXaxis()->FindBin("alphaJPsi"), alpha_jpsi.getError());

      referenceParameters->SetBinContent(parameters->GetXaxis()->FindBin("nJPsi"), n_jpsi.getVal());
      referenceParameters->SetBinError(parameters->GetXaxis()->FindBin("nJPsi"), n_jpsi.getError());

      referenceParameters->SetBinContent(parameters->GetXaxis()->FindBin("meanPsiPrime"), mean_psip.getVal());
      referenceParameters->SetBinError(parameters->GetXaxis()->FindBin("meanPsiPrime"), mean_psip.getError());

      referenceParameters->SetBinContent(parameters->GetXaxis()->FindBin("sigmaPsiPrime"), sigma_psip.getVal());
      referenceParameters->SetBinError(parameters->GetXaxis()->FindBin("sigmaPsiPrime"), sigma_psip.getError());
    }

  parameters->SetBinContent(parameters->GetXaxis()->FindBin("expBkg1"), c1.getVal());
  parameters->SetBinError(parameters->GetXaxis()->FindBin("expBkg1"), c1.getError());
  
  parameters->SetBinContent(parameters->GetXaxis()->FindBin("expBkg2"), c2.getVal());
  parameters->SetBinError(parameters->GetXaxis()->FindBin("expBkg2"), c2.getError());
  

  // Define the range to compute the results
  // We choose two sigma
  // We also need two new variables to compute the background
  Double_t lowerBoundJPsi = mean_jpsi.getVal() - 2.0*sigma_jpsi.getVal();
  Double_t upperBoundJPsi = mean_jpsi.getVal() + 2.0*sigma_jpsi.getVal();
  M.setRange("JPsiRange" ,lowerBoundJPsi, upperBoundJPsi);

  RooAbsReal *fracBkgRange1 = bkg_exp1.createIntegral(M, M, "JPsiRange");
  RooAbsReal *fracBkgRange2 = bkg_exp2.createIntegral(M, M, "JPsiRange");

  // Get the results, and fill the results histo
  Double_t nJPsi = fitJPsi.getVal();
  Double_t errJPsi = fitJPsi.getError();
  results->SetBinContent(results->GetXaxis()->FindBin("JPsi"), nJPsi);
  results->SetBinError(results->GetXaxis()->FindBin("JPsi"), errJPsi);

  Double_t nBkg = (fitBkg1.getVal() * fracBkgRange1->getVal() + fitBkg2.getVal() * fracBkgRange2->getVal());
  Double_t errBkg = (fitBkg1.getError() * fracBkgRange1->getVal() + fitBkg2.getError() * fracBkgRange2->getVal());
  results->SetBinContent(results->GetXaxis()->FindBin("bkg"), nBkg);
  results->SetBinError(results->GetXaxis()->FindBin("bkg"), errBkg);

  Double_t SB = 0.0;
  Double_t errSB = 0.0;
  if (nJPsi != 0.0 && nBkg != 0.0)
    {
      SB = nJPsi/nBkg;
      errSB = SB * (errJPsi/nJPsi + errBkg/nBkg);
    }
  results->SetBinContent(results->GetXaxis()->FindBin("SB"), SB);
  results->SetBinError(results->GetXaxis()->FindBin("SB"), errSB);

  Int_t nDF = fitResult->floatParsFinal().getSize();
  Double_t chi2 = plot->chiSquare("fitFunction", "invariantMassFit", nDF);
  results->SetBinContent(results->GetXaxis()->FindBin("chi2"), chi2);
  results->SetBinError(results->GetXaxis()->FindBin("chi2"), 0.0);

  fitList->AddAt(invariantMass, 0);
  fitList->AddAt(plot, 1);
  fitList->AddAt(parameters, 2);
  fitList->AddAt(results, 3);
}


void ProduceDimuonGraph(TFile *outputFile, std::vector<Double_t> multiplicityRanges)
{
  // This function will create all the graphs for dimuons (J/Psi and Background)
  // - raw
  // - yield
  // - yield over mean mult
  // - yield over mean mult normalised

  // First, get the CINT1B and CMUS1B graphs
  TGraphErrors *CINT1B = (TGraphErrors *) outputFile->Get("graphCINT1B");
  TGraphErrors *CMUS1B = (TGraphErrors *) outputFile->Get("graphCMUS1B");

  TList *JPsiGraph = new TList();
  JPsiGraph->SetName("JPsiGraph");
  TList *bkgGraph = new TList();
  bkgGraph->SetName("bkgGraph");

  // Raw Graphs
  TGraphAsymmErrors *rawJPsiGraph = new TGraphAsymmErrors(multiplicityRanges.size()-1);
  rawJPsiGraph->SetName("rawJPsiGraph");
  TGraphAsymmErrors *rawBkgGraph = new TGraphAsymmErrors(multiplicityRanges.size()-1);
  rawBkgGraph->SetName("rawBkgGraph");

  for (UInt_t iMult = 0; iMult < multiplicityRanges.size()-1; iMult++)
    {
      TString name;
      name.Form("Multiplicity_%d-%d", static_cast<Int_t>(multiplicityRanges[iMult]), static_cast<Int_t>(multiplicityRanges[iMult+1]));
      TH1D *results = (TH1D *) ((TList *) ((TList *) outputFile->Get("dimuonFitResults;1"))->FindObject(name))->FindObject("results");

      // J/Psi
      Double_t nJPsi = results->GetBinContent(results->GetXaxis()->FindBin("JPsi"));
      Double_t errJPsi = results->GetBinError(results->GetXaxis()->FindBin("JPsi"));

      rawJPsiGraph->SetPoint(iMult, CINT1B->GetX()[iMult], nJPsi);
      rawJPsiGraph->SetPointError(iMult,
				  CINT1B->GetX()[iMult] - multiplicityRanges[iMult], 
				  multiplicityRanges[iMult+1] - CINT1B->GetX()[iMult],
				  errJPsi, errJPsi);

      // Background
      Double_t nBkg = results->GetBinContent(results->GetXaxis()->FindBin("bkg"));
      Double_t errBkg = results->GetBinError(results->GetXaxis()->FindBin("bkg"));

      rawJPsiGraph->SetPoint(iMult, CINT1B->GetX()[iMult], nBkg);
      rawJPsiGraph->SetPointError(iMult,
				  CINT1B->GetX()[iMult] - multiplicityRanges[iMult], 
				  multiplicityRanges[iMult+1] - CINT1B->GetX()[iMult],
				  errBkg, errBkg);
    }

  JPsiGraph->AddAt(rawJPsiGraph, 0);
  bkgGraph->AddAt(rawBkgGraph, 0);

  // Yield Graphs
  TGraphAsymmErrors *yieldJPsiGraph = new TGraphAsymmErrors(multiplicityRanges.size()-1);
  yieldJPsiGraph->SetName("yieldJPsiGraph");
  TGraphAsymmErrors *yieldBkgGraph = new TGraphAsymmErrors(multiplicityRanges.size()-1);
  yieldBkgGraph->SetName("yieldBkgGraph");
  
  for (UInt_t iMult = 0; iMult < multiplicityRanges.size()-1; iMult++)
    {
      if (CINT1B->GetY()[iMult] != 0.0)
	{
	  // J/Psi
	  yieldJPsiGraph->SetPoint(iMult, rawJPsiGraph->GetX()[iMult], rawJPsiGraph->GetY()[iMult]/CINT1B->GetY()[iMult]);
	  if (rawJPsiGraph->GetY()[iMult] != 0.0)
	    {
	      Double_t error = yieldJPsiGraph->GetY()[iMult]*
		TMath::Sqrt((rawJPsiGraph->GetEYlow()[iMult]/rawJPsiGraph->GetY()[iMult])*(rawJPsiGraph->GetEYlow()[iMult]/rawJPsiGraph->GetY()[iMult]) + 
			    (CINT1B->GetEY()[iMult]/CINT1B->GetY()[iMult])*(CINT1B->GetEY()[iMult]/CINT1B->GetY()[iMult]));
	      yieldJPsiGraph->SetPointError(iMult,
					    rawJPsiGraph->GetEXlow()[iMult],
					    rawJPsiGraph->GetEXhigh()[iMult],
					    error, error);
	    }
	  // Background
	  yieldBkgGraph->SetPoint(iMult, rawBkgGraph->GetX()[iMult], rawBkgGraph->GetY()[iMult]/CINT1B->GetY()[iMult]);
	  if (rawBkgGraph->GetY()[iMult] != 0.0)
	    {
	      Double_t error = yieldBkgGraph->GetY()[iMult]*
		TMath::Sqrt((rawBkgGraph->GetEYlow()[iMult]/rawBkgGraph->GetY()[iMult])*(rawBkgGraph->GetEYlow()[iMult]/rawBkgGraph->GetY()[iMult]) + 
			    (CINT1B->GetEY()[iMult]/CINT1B->GetY()[iMult])*(CINT1B->GetEY()[iMult]/CINT1B->GetY()[iMult]));
		
	      yieldBkgGraph->SetPointError(iMult,
					   rawBkgGraph->GetEXlow()[iMult],
					   rawBkgGraph->GetEXhigh()[iMult],
					   error, error);
	    }
	}
      
      else 
	{
	  yieldJPsiGraph->SetPoint(iMult, rawJPsiGraph->GetX()[iMult], 0);
	  yieldBkgGraph->SetPoint(iMult, rawBkgGraph->GetX()[iMult], 0);
	}
    }

  JPsiGraph->AddAt(yieldJPsiGraph, 1);
  bkgGraph->AddAt(yieldBkgGraph, 1);

  // Yield over Mean Mult graph 
  TGraphAsymmErrors *yieldOverMeanMultJPsiGraph = new TGraphAsymmErrors(multiplicityRanges.size()-1);
  yieldOverMeanMultJPsiGraph->SetName("yieldOverMeanMultJPsiGraph");
  TGraphAsymmErrors *yieldOverMeanMultBkgGraph = new TGraphAsymmErrors(multiplicityRanges.size()-1);
  yieldOverMeanMultBkgGraph->SetName("yieldOverMeanMultBkgGraph");
  
  for (UInt_t iMult = 0; iMult < multiplicityRanges.size()-1; iMult++)
    {
      if (CINT1B->GetX()[iMult] != 0.0)
	{
	  // J/Psi
	  yieldOverMeanMultJPsiGraph->SetPoint(iMult, yieldJPsiGraph->GetX()[iMult], yieldJPsiGraph->GetY()[iMult]/CINT1B->GetX()[iMult]);
	  if (yieldJPsiGraph->GetY()[iMult] != 0.0)
	    {
	      Double_t error = yieldOverMeanMultJPsiGraph->GetY()[iMult]*
		TMath::Sqrt((yieldJPsiGraph->GetEYlow()[iMult]/yieldJPsiGraph->GetY()[iMult])*(yieldJPsiGraph->GetEYlow()[iMult]/yieldJPsiGraph->GetY()[iMult]) +
			    (CINT1B->GetEX()[iMult]/CINT1B->GetX()[iMult])*(CINT1B->GetEX()[iMult]/CINT1B->GetX()[iMult]));
	      yieldOverMeanMultJPsiGraph->SetPointError(iMult,
							yieldJPsiGraph->GetEXlow()[iMult],
							yieldJPsiGraph->GetEXhigh()[iMult],
							error, error);
	    }
	  // Background
	  yieldOverMeanMultBkgGraph->SetPoint(iMult, yieldBkgGraph->GetX()[iMult], yieldBkgGraph->GetY()[iMult]/CINT1B->GetX()[iMult]);
	  if (yieldBkgGraph->GetY()[iMult] != 0.0)
	    {
	      Double_t error = yieldOverMeanMultBkgGraph->GetY()[iMult]*
		TMath::Sqrt((yieldBkgGraph->GetEYlow()[iMult]/yieldBkgGraph->GetY()[iMult])*(yieldBkgGraph->GetEYlow()[iMult]/yieldBkgGraph->GetY()[iMult]) +
			    (CINT1B->GetEX()[iMult]/CINT1B->GetX()[iMult])*(CINT1B->GetEX()[iMult]/CINT1B->GetX()[iMult]));
	      yieldOverMeanMultBkgGraph->SetPointError(iMult,
						       yieldBkgGraph->GetEXlow()[iMult],
						       yieldBkgGraph->GetEXhigh()[iMult],
						       error, error);
	    }
	}
      
      else 
	{
	  yieldOverMeanMultJPsiGraph->SetPoint(iMult, yieldJPsiGraph->GetX()[iMult], 0);
	  yieldOverMeanMultBkgGraph->SetPoint(iMult, yieldBkgGraph->GetX()[iMult], 0);
	}
    }

  JPsiGraph->AddAt(yieldOverMeanMultJPsiGraph, 2);
  bkgGraph->AddAt(yieldOverMeanMultBkgGraph, 2);

  // Yield over Mean Mult normalised to get the bin with the highest number of CMUS1B equal to 1 
  Double_t maxCMUS1B = 0.0;
  Int_t referenceBin = 0;
  for (Int_t iBin = 0; iBin < CMUS1B->GetN(); iBin++)
    if (CMUS1B->GetX()[iBin] > 10.0)
      if (CMUS1B->GetY()[iBin] > maxCMUS1B)
	{
	  maxCMUS1B = CMUS1B->GetY()[iBin];
	  referenceBin = iBin;
	}
  
  TGraphAsymmErrors *yieldNormalisedJPsiGraph = new TGraphAsymmErrors(multiplicityRanges.size()-1);
  yieldNormalisedJPsiGraph->SetName("yieldNormalisedJPsiGraph");
  TGraphAsymmErrors *yieldNormalisedBkgGraph = new TGraphAsymmErrors(multiplicityRanges.size()-1);
  yieldNormalisedBkgGraph->SetName("yieldNormalisedBkgGraph");

  for (UInt_t iMult = 0; iMult < multiplicityRanges.size()-1; iMult++)
    {
      // JPsi
      if (yieldOverMeanMultJPsiGraph->GetY()[referenceBin] != 0.0)
	{
	  yieldNormalisedJPsiGraph->SetPoint(iMult, 
					     yieldOverMeanMultJPsiGraph->GetX()[iMult], 
					     yieldOverMeanMultJPsiGraph->GetY()[iMult]/yieldOverMeanMultJPsiGraph->GetY()[referenceBin]);
	  yieldNormalisedJPsiGraph->SetPointError(iMult,
						  yieldOverMeanMultJPsiGraph->GetEXlow()[iMult],
						  yieldOverMeanMultJPsiGraph->GetEXhigh()[iMult],
						  yieldOverMeanMultJPsiGraph->GetEYlow()[iMult]/yieldOverMeanMultJPsiGraph->GetY()[referenceBin],
						  yieldOverMeanMultJPsiGraph->GetEYhigh()[iMult]/yieldOverMeanMultJPsiGraph->GetY()[referenceBin]);
	}
      
      // Bkg
      if (yieldOverMeanMultBkgGraph->GetY()[referenceBin] != 0.0)
	{
	  yieldNormalisedBkgGraph->SetPoint(iMult, 
					    yieldOverMeanMultBkgGraph->GetX()[iMult], 
					    yieldOverMeanMultBkgGraph->GetY()[iMult]/yieldOverMeanMultBkgGraph->GetY()[referenceBin]);
	  yieldNormalisedBkgGraph->SetPointError(iMult,
						 yieldOverMeanMultBkgGraph->GetEXlow()[iMult],
						 yieldOverMeanMultBkgGraph->GetEXhigh()[iMult],
						 yieldOverMeanMultBkgGraph->GetEYlow()[iMult]/yieldOverMeanMultBkgGraph->GetY()[referenceBin],
						 yieldOverMeanMultBkgGraph->GetEYhigh()[iMult]/yieldOverMeanMultBkgGraph->GetY()[referenceBin]);
	}

    }

  JPsiGraph->AddAt(yieldNormalisedJPsiGraph, 3);
  bkgGraph->AddAt(yieldNormalisedBkgGraph, 3);

  outputFile->WriteTObject(JPsiGraph);
  outputFile->WriteTObject(bkgGraph);
}
