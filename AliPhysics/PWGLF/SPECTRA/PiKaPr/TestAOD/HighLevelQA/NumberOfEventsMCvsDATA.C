/////////////////////////////////////////////////////////
// NumberOfEventsMCvsDATA.C                            //
//                                                     //
// plots the ratio (# of MC events)/(# of DATA events) //
//    against the run number                           //
//                                                     //
// Written by John Groh                                //
/////////////////////////////////////////////////////////

// for cout
#include <iostream>
using namespace std;

// possible cuts to use -   0    1    2    3
Double_t CentCutMin[4] =   {0,  30,  30,  30};
Double_t CentCutMax[4] =   {5,  40,  40,  40};
Double_t QvecCutMin[4] =   {0,   0,   0, 1.5};
Double_t QvecCutMax[4] = {100, 100, 0.4, 100};
Double_t EtaMin[4] =    {-0.8,-0.8,-0.8,-0.8};
Double_t EtaMax[4] =     {0.8, 0.8, 0.8, 0.8};
Double_t Nsigmapid = 3.;
UInt_t trkbit = 1024;

//--------------------
// Start of function -
//--------------------
// icut selects which cut to use (2 and 3 have too low statistics and cause divide-by-zero errors, so don't use)
// nSigmaCut is the size of the deviation with which to determine outliers
void NumberOfEventsMCvsDATA(Int_t icut = 1, const Float_t nSigmaCut = 3)
{  
  // load libraries (I'm not sure I need all of these...)
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

  // fill an array with the run numbers and another two with the file names
  runList.open("output/AODQAChecks/RunList.txt");
  Int_t runs[nRuns];
  TString dataFilesDATA[nRuns];
  TString dataFilesMC[nRuns];

  for (Int_t irun=0; irun<nRuns; irun++)
    {
      runList >> runs[irun];
      dataFilesDATA[irun] = Form("output/AODQAChecks/DATA/AnalysisResults%i.root",runs[irun]);
      dataFilesMC[irun] = Form("output/AODQAChecks/MC/AnalysisResults%i.root",runs[irun]);
    }
  runList.close();

  // choose which TDirectories in the .root files to use
  TString snameDATA;
  TString snameMC;
  snameDATA = Form("OutputAODSpectraTask_Data_Cent%.0fto%.0f_QVec%.1fto%.1f_Eta%.1fto%.1f_%.1fSigmaPID_TrBit%d",CentCutMin[icut],CentCutMax[icut],QvecCutMin[icut],QvecCutMax[icut],EtaMin[icut],EtaMax[icut],Nsigmapid,trkbit);
  snameMC = Form("OutputAODSpectraTask_MC_Cent%.0fto%.0f_QVec%.1fto%.1f_Eta%.1fto%.1f_%.1fSigmaPID_TrBit%d",CentCutMin[icut],CentCutMax[icut],QvecCutMin[icut],QvecCutMax[icut],EtaMin[icut],EtaMax[icut],Nsigmapid,trkbit);

  TFile * fileDATA;
  TDirectoryFile * DirDATA;
  AliSpectraAODEventCuts * ecutsDATA;
  TFile * fileMC;
  TDirectoryFile * DirMC;
  AliSpectraAODEventCuts * ecutsMC;

  // histogram to be filled w/ the ratio (# of MC events)/(# of DATA events)
  TH1F * hNEventsRatio = new TH1F("hNEventsRatio","",nRuns,0,nRuns);

  // loop over all runs
  for (Int_t irun=0; irun<nRuns; irun++)
    {
      Printf("\n--- Processing run %i ---",runs[irun]);
      
      // open the file
      fileDATA = TFile::Open(dataFilesDATA[irun].Data());
      if (!fileDATA)
	{
	  Printf("\n!!! ERROR: Could not open DATA file %s !!!\n",dataFiles[irun].Data());
	  continue;
	}
      fileMC = TFile::Open(dataFilesMC[irun].Data());
      if (!fileMC)
	{
	  Printf("\n!!! ERROR: Could not open MC file %s !!!\n",dataFiles[irun].Data());
	  continue;
	}
      
      // choose the right directory and event cut objects
      DirDATA = (TDirectoryFile*)fileDATA->Get(Form("%s",snameDATA.Data()));
      if (!DirDATA)
	{
	  Printf("\n!!! ERROR: DirDATA is a null pointer. Skipping to next file !!!\n");
	  continue;
	}
      DirMC = (TDirectoryFile*)fileMC->Get(Form("%s",snameMC.Data()));
      if (!DirMC)
	{
	  Printf("\n!!! ERROR: DirMC is a null pointer. Skipping to next file !!!\n");
	  continue;
	}
      ecutsDATA = (AliSpectraAODEventCuts*)DirDATA->Get("Event Cuts");
      if (!ecutsDATA)
	{
	  Printf("\n!!! ERROR: ecutsDATA is a null pointer. Skipping to next file !!!\n");
	  continue;
	}      
      ecutsMC = (AliSpectraAODEventCuts*)DirMC->Get("Event Cuts");
      if (!ecutsMC)
	{
	  Printf("\n!!! ERROR: ecutsMC is a null pointer. Skipping to next file !!!\n");
	  continue;
	}      

      // calculate the ratio and fill the histogram
      Int_t nEventsDATA = ecutsDATA->NumberOfEvents();
      Int_t nEventsMC = ecutsMC->NumberOfEvents();
      Float_t ratio = (Float_t)nEventsMC / (Float_t)nEventsDATA;
      hNEventsRatio->SetBinContent(irun+1,ratio);
      Printf("# of MC events: %i\n# of DATA events: %i\nRatio: %.2f",nEventsMC,nEventsDATA,ratio);

    } // end loop over runs

  // draw it
  TCanvas * cNEventsRatio = new TCanvas("cNEventsRatio","cNEventsRatio");
  for (Int_t irun=0; irun<nRuns; irun++)
    hNEventsRatio->GetXaxis()->SetBinLabel(irun+1,Form("%i",runs[irun]));
  hNEventsRatio->GetYaxis()->SetTitle("(# of MC events)/(# of DATA events)");
  hNEventsRatio->GetYaxis()->SetTitleOffset(1.2);
  hNEventsRatio->SetMarkerStyle(21);
  hNEventsRatio->SetMarkerColor(kBlue);
  hNEventsRatio->SetStats(kFALSE);
  hNEventsRatio->DrawCopy("P");
}





