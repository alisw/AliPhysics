/////////////////////////////////////////////////////////////////////////////////////////
// AODQAChecks.C
//
// Performs some run-by-run high-level QA checks on AOD data 
//    for runs marked with global quality "good" in the RCT
//    for the LHC10h period
//
// Note that the following macros must be in the current
//    directory: CheckEfficiencies.C, CheckIntegratedRawYield.C,
//    CheckMatchingEfficiency.C, CheckNSigmaStability.C,
//    FindOutliers.C, and PlotAndSave.C
//
// Written by John Groh, summer student (Summer 2012)
// Advisors: Michele Floris and Alexander Kalweit
/////////////////////////////////////////////////////////////////////////////////////////

// for cout
#include <iostream>

// for processing .txt file with run numbers
#include <fstream>
using namespace std;

// global constants
const Int_t nPart = 3;
TString Particle[nPart] = {"Pion", "Kaon", "Proton"};
const Int_t nCharge = 2;
TString Charge[nCharge] = {"Pos", "Neg"};
TString Sign[nCharge] = {"Plus", "Minus"};
Int_t Color[nPart] = {1,2,4};
Int_t Marker[nPart*nCharge] = {20,21,22,24,25,26};
TString Names[nPart*nCharge] = {"#pi^{+}",
				"K^{+}",
				"p",
				"#pi^{-}",
				"K^{-}",
				"#bar{p}"};

// possible cuts to use -   0    1    2    3
Double_t CentCutMin[4] =   {0,  30,  30,  30};
Double_t CentCutMax[4] =   {5,  40,  40,  40};
Double_t QvecCutMin[4] =   {0,   0,   0, 1.5};
Double_t QvecCutMax[4] = {100, 100, 0.4, 100};
Double_t EtaMin[4] =    {-0.8,-0.8,-0.8,-0.8};
Double_t EtaMax[4] =     {0.8, 0.8, 0.8, 0.8};
Double_t Nsigmapid = 3.;
UInt_t trkbit = 1024;

// fixed Pt values chosen for efficiency and matching efficiency plots (GeV)
const Float_t FixedPtEff = 0.4;
const Float_t FixedPtMatchEff = 0.9;

// inclusions for functions called
#include "CheckIntegratedRawYield.C"
#include "CheckNSigmaStability.C"
#include "CheckEfficiencies.C"
#include "CheckMatchingEfficiency.C"
#include "PlotAndSave.C"
#include "FindOutliers.C"

//--------------------
// Start of function -
//--------------------
// useMC selects whether to use Monte Carlo or Data - must be run seperately for each
// icut selects which cut to use (2 and 3 have low statistics and cause divide-by-zero errors, so don't use)
// nSigmaCut is the size of the deviation with which to determine outliers
void AODQAChecks(Bool_t useMC = 1, Int_t icut = 1, const Float_t nSigmaCut = 3)
{
  Printf("\n\n--- Running AODQAChecks() with useMC = %i ---\n\n",useMC);

  gSystem->Load("libCore");
  gSystem->Load("libPhysics");
  gSystem->Load("libTree");
  gSystem->Load("libMatrix");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libOADB");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libTender");
  gSystem->Load("libCORRFW");
  gSystem->Load("libPWGTools");
  gSystem->Load("libPWGLFspectra");
  gSystem->Load("libProof");
  gSystem->AddIncludePath("-I$ALICE_ROOT/include");

  TString fold = "AODQAChecks";

  // get number of runs used
  Int_t nRunsTemp = 0;
  Int_t dummy;
  ifstream runList("output/AODQAChecks/RunList.txt");
  while (!runList.eof())
    {
      runList >> dummy;
      nRunsTemp++;
    }
  runList.close();

  // I know this is silly, but you need a const for array sizes.
  // Also, it seems to want to read one more line than there is in the file
  const Int_t nRuns = nRunsTemp - 1;

  // fill an array with the run numbers and another with the file names
  runList.open("output/AODQAChecks/RunList.txt");
  Int_t runs[nRuns];
  TString dataFiles[nRuns];
  
  for (Int_t irun=0; irun<nRuns; irun++)
    {
      runList >> runs[irun];
      if (useMC) dataFiles[irun] = Form("output/%s/MC/AnalysisResults%i.root",fold.Data(),runs[irun]);
      else dataFiles[irun] = Form("output/%s/DATA/AnalysisResults%i.root",fold.Data(),runs[irun]);
    }
  runList.close();

  // choose which TDirectory in the .root file to use
  TString sname;
  if (useMC)
    sname = Form("OutputAODSpectraTask_MC_Cent%.0fto%.0f_QVec%.1fto%.1f_Eta%.1fto%.1f_%.1fSigmaPID_TrBit%d",CentCutMin[icut],CentCutMax[icut],QvecCutMin[icut],QvecCutMax[icut],EtaMin[icut],EtaMax[icut],Nsigmapid,trkbit);
  else
    sname = Form("OutputAODSpectraTask_Data_Cent%.0fto%.0f_QVec%.1fto%.1f_Eta%.1fto%.1f_%.1fSigmaPID_TrBit%d",CentCutMin[icut],CentCutMax[icut],QvecCutMin[icut],QvecCutMax[icut],EtaMin[icut],EtaMax[icut],Nsigmapid,trkbit);

  // create a .root file for output
  TFile * fout = new TFile(Form("results/%s/Res_%s_AODQAChecks.root",fold.Data(),sname.Data()),"RECREATE");
  TFile * file;
  TDirectoryFile * Dir;
  AliSpectraAODTrackCuts * tcuts;
  AliSpectraAODEventCuts * ecuts;
  AliSpectraAODHistoManager * hman;

  // histograms for stability of nsigma distribution plots
  TH1F * TPCnsigMeanTrendPion = new TH1F("TPCnsigMeanTrendPion","",nRuns,0,nRuns);
  TH1F * TPCnsigMeanTrendKaon = new TH1F("TPCnsigMeanTrendKaon","",nRuns,0,nRuns);
  TH1F * TPCnsigMeanTrendProton = new TH1F("TPCnsigMeanTrendProton","",nRuns,0,nRuns);
  TH1F * TPCnsigSigmaTrendPion = new TH1F("TPCnsigSigmaTrendPion","",nRuns,0,nRuns);
  TH1F * TPCnsigSigmaTrendKaon = new TH1F("TPCnsigSigmaTrendKaon","",nRuns,0,nRuns);
  TH1F * TPCnsigSigmaTrendProton = new TH1F("TPCnsigSigmaTrendProton","",nRuns,0,nRuns);
  TH1F * TOFnsigMeanTrendPion = new TH1F("TOFnsigMeanTrendPion","",nRuns,0,nRuns);
  TH1F * TOFnsigMeanTrendKaon = new TH1F("TOFnsigMeanTrendKaon","",nRuns,0,nRuns);
  TH1F * TOFnsigMeanTrendProton = new TH1F("TOFnsigMeanTrendProton","",nRuns,0,nRuns);
  TH1F * TOFnsigSigmaTrendPion = new TH1F("TOFnsigSigmaTrendPion","",nRuns,0,nRuns);
  TH1F * TOFnsigSigmaTrendKaon = new TH1F("TOFnsigSigmaTrendKaon","",nRuns,0,nRuns);
  TH1F * TOFnsigSigmaTrendProton = new TH1F("TOFnsigSigmaTrendProton","",nRuns,0,nRuns);

  // histograms for integrated raw yield stability plots
  TH1F * IntegRawYieldAll = new TH1F("IntegRawYieldAll","",nRuns,0,nRuns);
  TH1F * IntegRawYieldPiPlus = new TH1F("IntegRawYieldPiPlus","",nRuns,0,nRuns);
  TH1F * IntegRawYieldPiMinus = new TH1F("IntegRawYieldPiMinus","",nRuns,0,nRuns);
  TH1F * IntegRawYieldKPlus = new TH1F("IntegRawYieldKPlus","",nRuns,0,nRuns);  
  TH1F * IntegRawYieldKMinus = new TH1F("IntegRawYieldKMinus","",nRuns,0,nRuns);
  TH1F * IntegRawYieldProton = new TH1F("IntegRawYieldProton","",nRuns,0,nRuns);
  TH1F * IntegRawYieldAntiproton = new TH1F("IntegRawYieldAntiproton","",nRuns,0,nRuns);

  // objects for efficiency plots
  if (useMC)
    {
      TCanvas * cEfficienciesAllRuns = new TCanvas("cEfficienciesAllRuns","cEfficienciesAllRuns",50,25,700,500);
      cEfficienciesAllRuns->Divide(3,2);
      TH1F * EfficiencyPiPlus = new TH1F("EfficiencyPiPlus","",nRuns,0,nRuns);
      TH1F * EfficiencyKPlus = new TH1F("EfficiencyKPlus","",nRuns,0,nRuns);
      TH1F * EfficiencyProton = new TH1F("EfficiencyProton","",nRuns,0,nRuns);
      TH1F * EfficiencyPiMinus = new TH1F("EfficiencyPiMinus","",nRuns,0,nRuns);
      TH1F * EfficiencyKMinus = new TH1F("EfficiencyKMinus","",nRuns,0,nRuns);
      TH1F * EfficiencyAntiproton = new TH1F("EfficiencyAntiproton","",nRuns,0,nRuns);
    }
  else
    {
      TCanvas * cEfficienciesAllRuns = 0x0;
      TH1F * EfficiencyPiPlus = 0x0;
      TH1F * EfficiencyKPlus = 0x0;
      TH1F * EfficiencyProton = 0x0;
      TH1F * EfficiencyPiMinus = 0x0;
      TH1F * EfficiencyKMinus = 0x0;
      TH1F * EfficiencyAntiproton = 0x0;
    }


  // objects for matching efficiency plots
  TH1F * MatchEffPos = new TH1F("MatchEffPos","",nRuns,0,nRuns);
  MatchEffPos->Sumw2();
  TH1F * MatchEffNeg = new TH1F("MatchEffNeg","",nRuns,0,nRuns);
  MatchEffNeg->Sumw2();

  // overall loop over data files (1 per run)
  for (Int_t irun=0; irun<nRuns; irun++)
    {
      Printf("\n--- Processing data from run %i ---\n", runs[irun]);

      // open file and print all objects in file
      file = TFile::Open(dataFiles[irun].Data());     
      if (!file)
	{
	  Printf("\n\n!!! ERROR: Could not open file %s !!!\n\n",dataFiles[irun].Data());
	  continue;
	}
      file->Print();
      Printf("sname: %s", sname.Data());

      // choose the right directory, event, and histo managers
      Dir = (TDirectoryFile*)file->Get(Form("%s",sname.Data()));
      if (!Dir)
	{
	  Printf("\n\n!!! ERROR: Dir is a null pointer. Skipping to next file !!!\n\n");
	  continue;
	}
      ecuts = (AliSpectraAODEventCuts*)Dir->Get("Event Cuts");
      if (!ecuts)
	{
	  Printf("\n\n!!! ERROR: ecuts is a null pointer. Skipping to next file !!!\n\n");
	  continue;
	}      
      tcuts = (AliSpectraAODTrackCuts*)Dir->Get("Track Cuts");
      if (!tcuts)
	{
	  Printf("\n\n!!! ERROR: tcuts is a null pointer. Skipping to next file !!!\n\n");
	  continue;
	}      
      hman = (AliSpectraAODHistoManager*)Dir->Get("SpectraHistos");
      if (!hman)
	{
	  Printf("\n\n!!! ERROR: hman is a null pointer. Skipping to next file !!!\n\n");
	  continue;
	}

      //---------------------------------------------------------------------------------------
      // Find the trends in the means and sigmas of fits to the peaks of fixed-Pt projections -
      // of the nsigma distributions of the TPC and TOF as a function of the run number       -
      //---------------------------------------------------------------------------------------
      CheckNSigmaStability(hman,
			   TPCnsigMeanTrendPion,
			   TPCnsigMeanTrendKaon,
			   TPCnsigMeanTrendProton, 
			   TPCnsigSigmaTrendPion,
			   TPCnsigSigmaTrendKaon,
			   TPCnsigSigmaTrendProton,
			   TOFnsigMeanTrendPion,
			   TOFnsigMeanTrendKaon,
			   TOFnsigMeanTrendProton,
			   TOFnsigSigmaTrendPion,
			   TOFnsigSigmaTrendKaon,
			   TOFnsigSigmaTrendProton,
			   runs,
			   nRuns,
			   irun,
			   useMC);

      //---------------------------------------------------------------------------------
      // plot the stability of the integrated raw yield as a function of the run number -
      //---------------------------------------------------------------------------------
      CheckIntegratedRawYield(ecuts,
			      hman,
			      IntegRawYieldAll,
			      IntegRawYieldPiPlus,
			      IntegRawYieldPiMinus,
			      IntegRawYieldKPlus,
			      IntegRawYieldKMinus,
			      IntegRawYieldProton,
			      IntegRawYieldAntiproton,
			      runs,
			      irun,
			      nRuns,
			      Names,
			      useMC);
      
      //-----------------------------------------------------------------------------
      // find the trend in the MC correction factor as a function of the run number -
      //-----------------------------------------------------------------------------
      if (useMC) CheckEfficiencies(hman,
				   cEfficienciesAllRuns,
				   EfficiencyPiPlus,
				   EfficiencyKPlus,
				   EfficiencyProton,
				   EfficiencyPiMinus,
				   EfficiencyKMinus,
				   EfficiencyAntiproton,
				   FixedPtEff,
				   runs,
				   nRuns,
				   irun);

      //-------------------------------------------------------------------------------
      // for a fixed Pt, plot the matching efficiency as a function of the run number -
      //-------------------------------------------------------------------------------
      CheckMatchingEfficiency(tcuts,
			      MatchEffPos,
			      MatchEffNeg,
			      FixedPtMatchEff,
			      runs,
			      irun,
			      nRuns,
			      useMC);

      //---------------------------------------------------------------------------
      // save eta-phi distributions for each run on a seperate page of a pdf file -
      //---------------------------------------------------------------------------
      // canvas for temporarily drawing
      TCanvas * cEtaPhi = new TCanvas("cEtaPhi","cEtaPhi");
      TH2F * hEtaPhi = (TH2F*)((TH2F*)tcuts->GetHistoEtaPhiHighPt()->Clone("hEtaPhi"));
      hEtaPhi->Rebin2D(4);
      hEtaPhi->Scale(1./hEtaPhi->GetEntries());
      gPad->SetGridy();
      gPad->SetGridx();
      hEtaPhi->DrawClone("COLZ");

      // legend
      TLegend * lEtaPhi = new TLegend(0.4,0.85,.6,.95);
      lEtaPhi->SetFillColor(0);
      lEtaPhi->AddEntry(hEtaPhi,Form("Run %i",runs[irun]),"");
      if (useMC) lEtaPhi->AddEntry(hEtaPhi,"MC","");
      else lEtaPhi->AddEntry(hEtaPhi,"DATA","");
      lEtaPhi->DrawClone();

      // save to a pdf and close the temporary canvas
      if (useMC)
	{
	  if (irun == 0) cEtaPhi->SaveAs("Plots/MC/EtaPhi.pdf(","pdf");
	  else if (irun < nRuns-1) cEtaPhi->SaveAs("Plots/MC/EtaPhi.pdf","pdf");
	  else if (irun == nRuns-1) cEtaPhi->SaveAs("Plots/MC/EtaPhi.pdf)","pdf");
	  cEtaPhi->Close();
	}
      else
	{
	  if (irun == 0) cEtaPhi->SaveAs("Plots/DATA/EtaPhi.pdf(","pdf");
	  else if (irun < nRuns-1) cEtaPhi->SaveAs("Plots/DATA/EtaPhi.pdf","pdf");
	  else if (irun == nRuns-1) cEtaPhi->SaveAs("Plots/DATA/EtaPhi.pdf)","pdf");
	  cEtaPhi->Close();
	}
      
      // close the file and clean up before the next iteration
      file->Close();
      //delete file;
      //delete Dir;
      //delete hman;

      Printf("\n--- Done Processing data from run %i ---\n", runs[irun]);

    }

  //------------------------------------------------------------
  // Plot all the results and save them to the output file     -
  //------------------------------------------------------------
  PlotAndSave(runs,
	      nRuns,
	      useMC,
	      fout,
	      
	      TPCnsigMeanTrendPion,
	      TPCnsigMeanTrendKaon,
	      TPCnsigMeanTrendProton,
	      TPCnsigSigmaTrendPion,
	      TPCnsigSigmaTrendKaon,
	      TPCnsigSigmaTrendProton,
	      TOFnsigMeanTrendPion,
	      TOFnsigMeanTrendKaon,
	      TOFnsigMeanTrendProton,
	      TOFnsigSigmaTrendPion,
	      TOFnsigSigmaTrendKaon,
	      TOFnsigSigmaTrendProton,
	      
	      IntegRawYieldAll,
	      IntegRawYieldPiPlus,
	      IntegRawYieldKPlus,
	      IntegRawYieldProton,
	      IntegRawYieldPiMinus,
	      IntegRawYieldKMinus,
	      IntegRawYieldAntiproton,
	      
	      EfficiencyPiPlus,
	      EfficiencyKPlus,
	      EfficiencyProton,
	      EfficiencyPiMinus,
	      EfficiencyKMinus,
	      EfficiencyAntiproton,
	      
	      MatchEffPos,
	      MatchEffNeg);

  //------------------------------------------------------------------------------
  // Find the outliers from the efficiency, raw yield, matching efficiency, and  -
  // nsigma fit plots and print the corresponding run numbers to the console.    -
  // Also creates a .txt file in the results/ directory with info about outliers -
  //------------------------------------------------------------------------------
  FindOutliers(runs,
	       nRuns,
	       useMC,
	       icut,
	       nSigmaCut,
	       fout,
	       
	       TPCnsigMeanTrendPion,
	       TPCnsigMeanTrendKaon,
	       TPCnsigMeanTrendProton,
	       TPCnsigSigmaTrendPion,
	       TPCnsigSigmaTrendKaon,
	       TPCnsigSigmaTrendProton,
	       TOFnsigMeanTrendPion,
	       TOFnsigMeanTrendKaon,
	       TOFnsigMeanTrendProton,
	       TOFnsigSigmaTrendPion,
	       TOFnsigSigmaTrendKaon,
	       TOFnsigSigmaTrendProton,
	       
	       IntegRawYieldAll,
	       
	       EfficiencyPiPlus,
	       EfficiencyKPlus,
	       EfficiencyProton,
	       EfficiencyPiMinus,
	       EfficiencyKMinus,
	       EfficiencyAntiproton,
	       
	       MatchEffPos,
	       MatchEffNeg);

  
  // close the output file
  fout->Close();
  // also close the unused efficiency canvas
  if (useMC) cEfficienciesAllRuns->Close();
  
  Printf("\n\n--- End of AODQAChecks() ---\n");
}







