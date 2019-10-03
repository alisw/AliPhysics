#include <TFile.h>
#include <TH1.h>
#include <TMath.h>
#include <TGraphErrors.h>

#include "my_tools.C"

#include <iostream>

using namespace std;

/*
  Info:

  Method to convert charged pT histograms to have the same histogram binning
  as our histograms have.  They in general go to 100 GeV/c in pT while we only
  go to 50. In that range the binning is EXACTLY the same (and the code checks
  it).

  In charged spectra there is NO NORMALIZATION uncertainty. We add it by hand
  to the syst error.

  * Need at some point to get the final spectra
  * Need at some point to get the final normalization uncertainty

  
  To run:
  
  gSystem->AddIncludePath("-I../macros")
  gSystem->AddIncludePath("-I../lib")
  gROOT->SetMacroPath(".:../macros")
  .L my_tools.C+
  .L convertToMyBins.C+


  ConvertPPfile("charged_spectra/dNdPt_pp_7TeV.root",    "NOPID_pp_7TeV.root", 0.034);
  ConvertPPfile("charged_spectra/dNdPt_pp_2.76TeV.root", "NOPID_pp_2.76TeV.root", 0.034);
  ConvertPPfile("charged_spectra/dNdPt_pp_900GeV.root",  "NOPID_pp_900GeV.root", 0.034);




  Older methods that can be useful for AA:
  The method convertToMyBinsAA is for pp and AA spectra.
  The method convertToMyBinsRaa is for RAA.
  The method convertToMyBinsOld is for AA data when the centralities had to be merged.

 */



TH1D* Convert(TH1D* histIn, Double_t addSyst=0);
TH1D* ConvertGraph(TGraphErrors* graphIn);
void convertToMyBinsAA();
void convertToMyBinsRaa(const Char_t* fileName = "PbPb_RAA_2760GeV_200511.root");
void convertToMyBinsOld(const Char_t* fileName = "dndpt_all_2011-05-15.root");


void  ConvertPPfile(const Char_t* inFileName, const Char_t* outFileName, Double_t normSyst=0.0)
{
  TFile* fileDataCharged = FindFileFresh(inFileName);
  
  if(!fileDataCharged)
    return;
  
  TH1D* hPtHelp = (TH1D*)fileDataCharged->Get("dNdPt_stat");
  TH1D* hPtCharged = Convert(hPtHelp);
  hPtCharged->SetName("hDnDptCharged_PP");
  hPtCharged->SetMarkerStyle(29);

  hPtHelp = (TH1D*)fileDataCharged->Get("dNdPt_syst");
  TH1D* hPtChargedSyst = Convert(hPtHelp, normSyst);
  hPtChargedSyst->SetName("hDnDptChargedSyst_PP");
  hPtChargedSyst->SetMarkerStyle(29);

  TFile* fileOut = new TFile(outFileName, "RECREATE");
  hPtCharged->Write();
  hPtChargedSyst->Write();
  fileOut->Close();  
}

//________________________________________________________________________________________
TH1D* Convert(TH1D* histIn, Double_t addSyst)
{
  const Int_t nPtBins = 68;
  Double_t xBins[nPtBins+1] = {0. ,  0.05, 0.1,  0.15, 0.2,  0.25, 0.3,  0.35, 0.4,  0.45,
			       0.5,  0.55, 0.6,  0.65, 0.7,  0.75, 0.8,  0.85, 0.9,  0.95,
			       1.0,  1.1 , 1.2,  1.3 , 1.4,  1.5 , 1.6,  1.7 , 1.8,  1.9 ,
			       2.0,  2.2 , 2.4,  2.6 , 2.8,  3.0 , 3.2,  3.4 , 3.6,  3.8 ,
			       4.0,  4.5 , 5.0,  5.5 , 6.0,  6.5 , 7.0,  8.0 , 9.0,  10.0,
			       11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 18.0, 20.0, 22.0, 24.0,
			       26.0, 28.0, 30.0, 32.0, 34.0, 36.0, 40.0, 45.0, 50.0 };
  

  TH1D* hist = new TH1D("hist", "", nPtBins, xBins);
  hist->SetTitle(histIn->GetTitle());
  hist->GetXaxis()->SetTitle(histIn->GetXaxis()->GetTitle());
  hist->GetYaxis()->SetTitle(histIn->GetYaxis()->GetTitle());

  const Double_t deltapt = 0.0001;

  for(Int_t bin = 1; bin <= nPtBins; bin++) {

    // check bin size
    if(TMath::Abs(hist->GetXaxis()->GetBinLowEdge(bin) -
		  histIn->GetXaxis()->GetBinLowEdge(bin)) > deltapt) {
      cout << "pt edge low does not agree!!!!!!!!!!!!!!!" << endl;
      return 0;
    }
    if(TMath::Abs(hist->GetXaxis()->GetBinUpEdge(bin) -
		  histIn->GetXaxis()->GetBinUpEdge(bin)) > deltapt) {
      cout << "pt edge high does not agree!!!!!!!!!!!!!!!" << endl;
      return 0;
    }
    if(TMath::Abs(hist->GetXaxis()->GetBinCenter(bin) -
		  histIn->GetXaxis()->GetBinCenter(bin)) > deltapt) {
      cout << "pt center does not agree!!!!!!!!!!!!!!!" << endl;
      return 0;
    }


    hist->SetBinContent(bin, histIn->GetBinContent(bin)); 

    Double_t error = histIn->GetBinError(bin); 
    if(addSyst) {

      Double_t extra_error = addSyst*histIn->GetBinContent(bin);
      error = TMath::Sqrt(error*error + extra_error*extra_error);
    }
    hist->SetBinError(bin, error);
  }

  return hist;
}

//________________________________________________________________________________________
TH1D* ConvertGraph(TGraphErrors* graphIn)
{
  const Int_t nPtBins = 68;
  Double_t xBins[nPtBins+1] = {0. ,  0.05, 0.1,  0.15, 0.2,  0.25, 0.3,  0.35, 0.4,  0.45,
			       0.5,  0.55, 0.6,  0.65, 0.7,  0.75, 0.8,  0.85, 0.9,  0.95,
			       1.0,  1.1 , 1.2,  1.3 , 1.4,  1.5 , 1.6,  1.7 , 1.8,  1.9 ,
			       2.0,  2.2 , 2.4,  2.6 , 2.8,  3.0 , 3.2,  3.4 , 3.6,  3.8 ,
			       4.0,  4.5 , 5.0,  5.5 , 6.0,  6.5 , 7.0,  8.0 , 9.0,  10.0,
			       11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 18.0, 20.0, 22.0, 24.0,
			       26.0, 28.0, 30.0, 32.0, 34.0, 36.0, 40.0, 45.0, 50.0 };
  

  TH1D* hist = new TH1D("hist", "", nPtBins, xBins);
  // hist->SetTitle(graphIn->GetTitle());
  // hist->GetXaxis()->SetTitle(graphIn->GetXaxis()->GetTitle());
  // hist->GetYaxis()->SetTitle(graphIn->GetYaxis()->GetTitle());

  const Double_t deltapt = 0.0001;

  // Jacek does not fill the first 3 bins of his graph?
  for(Int_t bin = 4; bin <= nPtBins; bin++) {

    // check bin size
    if(TMath::Abs(hist->GetXaxis()->GetBinCenter(bin) -
		  graphIn->GetX()[bin-1]) > deltapt) {
      cout << "pt center does not agree!!!!!!!!!!!!!!!" << endl;
      return 0;
    }


    hist->SetBinContent(bin, graphIn->GetY()[bin-1]); 
    hist->SetBinError(bin, graphIn->GetEY()[bin-1]); 
  }

  return hist;
}

//________________________________________________________________________________________
void convertToMyBinsAA()
{
  const Int_t nCent = 7;
  const Int_t cent[nCent+1] = {-5, 0, 5, 10, 20, 40, 60, 80};
  const Int_t colors[nCent] = {kBlack, kRed+1, kOrange+7, kOrange, kGreen+2, kCyan+2, kBlue+1};  

  TH1D* hPtCharged[nCent] = {0, 0, 0, 0, 0, 0, 0};
  TH1D* hPtChargedSyst[nCent] = {0, 0, 0, 0, 0, 0, 0};

  for (Int_t i = 0; i < nCent; i++) {
    
    if(cent[i]==-5) {

      TFile* fileDataCharged = FindFileFresh("dNdPt_pp_2.76TeV.root");
      
      if(!fileDataCharged)
	return;
      
      TH1D* hPtHelp = (TH1D*)fileDataCharged->Get("dNdPt_stat");
      hPtCharged[i] = Convert(hPtHelp);
      hPtCharged[i]->SetName("hDnDptCharged_PP");
      hPtCharged[i]->SetMarkerColor(colors[i]);
      hPtCharged[i]->SetMarkerStyle(29);
      hPtHelp = (TH1D*)fileDataCharged->Get("dNdPt_syst");
      hPtChargedSyst[i] = Convert(hPtHelp, 0.034); // 3.4 % norm error
      /*
	Dear Peter,
	
	for 2.76TeV the systematic error for MB->INEL is 3.40%,, we decided to 
	drop the uncertainty from VdM scan since we don't really use it...
	
	Michael	
      */
      hPtChargedSyst[i]->SetName("hDnDptChargedSyst_PP");
      hPtChargedSyst[i]->SetMarkerColor(colors[i]);
      hPtChargedSyst[i]->SetMarkerStyle(29);
    } else {

      TFile* fileDataCharged = FindFileFresh(Form("dNdPt_PbPb_2.76ATeV_centrality_%d-%d.root", cent[i], cent[i+1]));
      
      if(!fileDataCharged)
	return;

      TH1D* hPtHelp = (TH1D*)fileDataCharged->Get("dNdPt_stat");
      hPtCharged[i] = Convert(hPtHelp);
      hPtCharged[i]->SetName(Form("hDnDptCharged_%d_%d", cent[i], cent[i+1]));
      hPtCharged[i]->SetMarkerColor(colors[i]);
      hPtCharged[i]->SetMarkerStyle(29);
      hPtHelp = (TH1D*)fileDataCharged->Get("dNdPt_syst");
      hPtChargedSyst[i] = Convert(hPtHelp);
      hPtChargedSyst[i]->SetName(Form("hDnDptChargedSyst_%d_%d", cent[i], cent[i+1]));
      hPtChargedSyst[i]->SetMarkerColor(colors[i]);
      hPtChargedSyst[i]->SetMarkerStyle(29);
    }
  }

  TFile* fileOut = new TFile("dNdPt_Charged.root", "RECREATE");
  for (Int_t i = 0; i < nCent; i++) {
    hPtCharged[i]->Write();
    hPtChargedSyst[i]->Write();
  }
  fileOut->Close();  
}

//________________________________________________________________________________________
void convertToMyBinsRaa(const Char_t* fileName)
{
  TFile* fileData = FindFileFresh(fileName);

  if(!fileData)
    return;
  
  TGraphErrors* hRaa_0_5 = (TGraphErrors*)fileData->Get("raa_c0_5_stat");
  TH1D* hRaaCharged_0_5 = ConvertGraph(hRaa_0_5);
  hRaaCharged_0_5->SetName("hRaaCharged_0_5");

  TGraphErrors* hRaa_5_10 = (TGraphErrors*)fileData->Get("raa_c5_10_stat");
  TH1D* hRaaCharged_5_10 = ConvertGraph(hRaa_5_10);
  hRaaCharged_5_10->SetName("hRaaCharged_5_10");

  TGraphErrors* hRaa_10_20 = (TGraphErrors*)fileData->Get("raa_c10_20_stat");
  TH1D* hRaaCharged_10_20 = ConvertGraph(hRaa_10_20);
  hRaaCharged_10_20->SetName("hRaaCharged_10_20");

  TGraphErrors* hRaa_20_40 = (TGraphErrors*)fileData->Get("raa_c20_40_stat");
  TH1D* hRaaCharged_20_40 = ConvertGraph(hRaa_20_40);
  hRaaCharged_20_40->SetName("hRaaCharged_20_40");

  TGraphErrors* hRaa_40_60 = (TGraphErrors*)fileData->Get("raa_c40_60_stat");
  TH1D* hRaaCharged_40_60 = ConvertGraph(hRaa_40_60);
  hRaaCharged_40_60->SetName("hRaaCharged_40_60");

  TGraphErrors* hRaa_60_80 = (TGraphErrors*)fileData->Get("raa_c60_80_stat");
  TH1D* hRaaCharged_60_80 = ConvertGraph(hRaa_60_80);
  hRaaCharged_60_80->SetName("hRaaCharged_60_80");

  TGraphErrors* hRaa_0_20 = (TGraphErrors*)fileData->Get("raa_c0_20_stat");
  TH1D* hRaaCharged_0_20 = ConvertGraph(hRaa_0_20);
  hRaaCharged_0_20->SetName("hRaaCharged_0_20");

  TGraphErrors* hRaa_40_80 = (TGraphErrors*)fileData->Get("raa_c40_80_stat");
  TH1D* hRaaCharged_40_80 = ConvertGraph(hRaa_40_80);
  hRaaCharged_40_80->SetName("hRaaCharged_40_80");


  TFile* fileOut = new TFile("Raa_Charged.root", "RECREATE");
  hRaaCharged_0_5->Write();
  hRaaCharged_5_10->Write();
  hRaaCharged_10_20->Write();
  hRaaCharged_20_40->Write();
  hRaaCharged_40_60->Write();
  hRaaCharged_60_80->Write();
  hRaaCharged_0_20->Write();
  hRaaCharged_40_80->Write();
  fileOut->Close();  
}

//________________________________________________________________________________________

void convertToMyBinsOld(const Char_t* fileName)
{
  TFile* fileData = FindFileFresh(fileName);

  if(!fileData)
    return;

  TH1D* hDnDPt_0_5 = (TH1D*)fileData->Get("dNdPt_PbPb_c0_5");
  TH1D* hDnDptCharged_0_5 = Convert(hDnDPt_0_5);
  hDnDptCharged_0_5->SetName("hDnDptCharged_0_5");

  TH1D* hDnDPt_5_10 = (TH1D*)fileData->Get("dNdPt_PbPb_c5_10");
  TH1D* hDnDptCharged_5_10 = Convert(hDnDPt_5_10);
  hDnDptCharged_5_10->SetName("hDnDptCharged_5_10");

  TH1D* hDnDPt_10_20 = (TH1D*)fileData->Get("dNdPt_PbPb_c10_20");
  TH1D* hDnDptCharged_10_20 = Convert(hDnDPt_10_20);
  hDnDptCharged_10_20->SetName("hDnDptCharged_10_20");

  TH1D* hDnDPt_20_40 = (TH1D*)fileData->Get("dNdPt_PbPb_c20_30");
  TH1D* hHelp =  (TH1D*)fileData->Get("dNdPt_PbPb_c30_40");
  hDnDPt_20_40->Add(hHelp);
  hDnDPt_20_40->Scale(1.0/2.0);
  TH1D* hDnDptCharged_20_40 = Convert(hDnDPt_20_40);
  hDnDptCharged_20_40->SetName("hDnDptCharged_20_40");

  TH1D* hDnDPt_40_60 = (TH1D*)fileData->Get("dNdPt_PbPb_c40_50");
  hHelp =  (TH1D*)fileData->Get("dNdPt_PbPb_c50_60");
  hDnDPt_40_60->Add(hHelp);
  hDnDPt_40_60->Scale(1.0/2.0);
  TH1D* hDnDptCharged_40_60 = Convert(hDnDPt_40_60);
  hDnDptCharged_40_60->SetName("hDnDptCharged_40_60");

  TH1D* hDnDPt_60_80 = (TH1D*)fileData->Get("dNdPt_PbPb_c60_70");
  hHelp =  (TH1D*)fileData->Get("dNdPt_PbPb_c70_80");
  hDnDPt_60_80->Add(hHelp);
  hDnDPt_60_80->Scale(1.0/2.0);
  TH1D* hDnDptCharged_60_80 = Convert(hDnDPt_60_80);
  hDnDptCharged_60_80->SetName("hDnDptCharged_60_80");


  TFile* fileOut = new TFile("dNdPt_Charged.root", "RECREATE");
  hDnDptCharged_0_5->Write();
  hDnDptCharged_5_10->Write();
  hDnDptCharged_10_20->Write();
  hDnDptCharged_20_40->Write();
  hDnDptCharged_40_60->Write();
  hDnDptCharged_60_80->Write();
  fileOut->Close();  
}
