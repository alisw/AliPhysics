#include <TCanvas.h>

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TGrid.h>
#include <TFile.h>
#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TList.h>
#include <TMath.h>
#include <TGraphErrors.h>
#include <TString.h>
#include "Riostream.h"
#include "stdio.h"
using namespace std;
#endif

//-----------------------------------------------------------------------------
// Function declaration
void QAFillEventSelection();
void QAFillOccupancy();
void QAFillClusters();
void QAFillTracks();
void QAFillNPi0();

void QAWriteEventSelection();
void QAWriteOccupancy();
void QAWriteClusters();
void QAWriteTracks();
void QAWriteNPi0();

TH1D *Pi0InvMass(TH2 *hMassPt, Float_t ptmin, Float_t ptmax);
TH1F *AddWriteTH1F(const char *name, const char *title, Double_t* values, Double_t* errors);
void Fit1Pi0(TH1D *hMass,
	     const Int_t polN,
	     const Double_t mFitMin,
	     const Double_t mFitMax,
	     Double_t &nPi0, Double_t &ePi0);
Double_t pi0massP1(Double_t *x, Double_t *par);
Double_t pi0massP2(Double_t *x, Double_t *par);
Double_t pi0mass  (Double_t *x, Double_t *par);
Double_t bgP2     (Double_t *x, Double_t *par);

//-----------------------------------------------------------------------------
// Global variabes
const Int_t kNCents = 1;
const Int_t kNPID = 8+4;
const char* kPIDNames[kNPID] = {"All", "Allwou", "Disp", "Disp2", "Dispwou", "CPV", "CPV2", "Both",
				"Allcore", "Dispcore", "CPVcore", "Bothcore"};

Int_t runIndex;
Int_t runNum[500];
TList *listHist;
Double_t run[500], erun[500];
// Event selection:
Double_t nTotal4[500];
Double_t rVtx[500], rVtxZ10[500], rVtxZ10Cent[500];
Double_t eVtx[500], eVtxZ10[500], eVtxZ10Cent[500];
// Occupancy:
Double_t nCellsM1[500], nCellsM2[500], nCellsM3[500];
// Clusters:
Double_t nCluEvent[500], eCluEvent[500];
Double_t cluEnergy[500], ecluEnergy[500];
Double_t nPhotPID[kNCents][kNPID][500], enPhotPID[kNCents][kNPID][500];
Double_t mEnPID[kNCents][kNPID][500],   emEnPID[kNCents][kNPID][500];
Double_t nPhotPIDHigh[kNCents][kNPID][500], enPhotPIDHigh[kNCents][kNPID][500];
Double_t mEnPIDHigh[kNCents][kNPID][500],   emEnPIDHigh[kNCents][kNPID][500];


// Tracks:
Double_t nTracks0[500], eTracks0[500];
Double_t nTracks1[500], eTracks1[500];
// Pi0:
Double_t mPi0[500], emPi0[500];
Double_t wPi0[500], ewPi0[500];
Double_t nPi0Event[500], ePi0Event[500];


//-----------------------------------------------------------------------------
void ExtractQA(const TString filelist="filelist.txt",
	       const TString runlist="runlist.txt")
{
  // Loop over per-run root files and run various QA functions

  //TGrid::Connect("alien://");

  ifstream in;
  in.open(filelist.Data());

  ifstream inRuns;
  inRuns.open(runlist.Data());

  char rootFileName[256];
  TFile *rootFile;

  Int_t runNumber=0;
  runIndex = 0;
  while (1) {
    in >> rootFileName;
    if (!in.good()) break;
    inRuns >> runNumber;
    if (!inRuns.good()) break;
    
    if(169094 == runNumber )
      continue;
    
    runNum[runIndex] = runNumber;

    printf("root file is %s, run # = %d\n",rootFileName,runNumber);
    // char *runNum = strtok(rootFileName+35,".");
    rootFile = TFile::Open(rootFileName,"read");
    listHist = (TList*)rootFile->Get("PHOSPi0Flow/PHOSPi0FlowCoutput1");

    run[runIndex]            = runIndex+1;

    QAFillEventSelection();
    QAFillOccupancy();
    QAFillClusters();
    QAFillTracks();
//     QAFillNPi0();

    listHist->Clear();
    rootFile->Close();
    runIndex++;
  }

  
  TFile *fileQA = TFile::Open("runQA.root","recreate");
  QAWriteEventSelection();
  QAWriteOccupancy();
  QAWriteClusters();
  QAWriteTracks();
//   QAWriteNPi0();
  fileQA         ->Close();

}

//-----------------------------------------------------------------------------
void QAFillEventSelection()
{
  TH1F *hSelEvents = (TH1F*)listHist->FindObject("hTotSelEvents");
  
  Double_t nTotal      = hSelEvents->GetBinContent(1);
  Int_t nVtx           = hSelEvents->GetBinContent(2);
  Int_t nVtxZ10        = hSelEvents->GetBinContent(3);
  Int_t nVtxZ10Cent    = hSelEvents->GetBinContent(4);
  
  if( ! nTotal )
    return;
  
  nTotal4[runIndex] = hSelEvents->GetBinContent(4);
  
  rVtx[runIndex]           = nVtx          /nTotal;
  rVtxZ10[runIndex]        = nVtxZ10       /nTotal;
  rVtxZ10Cent[runIndex]    = nVtxZ10Cent   /nTotal;
  
  
  //Change in logic (plus clarification through rewriting): t(1+p)/n/n -> p(1-p)/n:
  eVtx[runIndex]           = TMath::Sqrt(       rVtx[runIndex] * (1 - rVtx[runIndex])        /nTotal); 
  eVtxZ10[runIndex]        = TMath::Sqrt(    rVtxZ10[runIndex] * (1 - rVtxZ10[runIndex])     /nTotal);
  eVtxZ10Cent[runIndex]    = TMath::Sqrt(rVtxZ10Cent[runIndex] * (1 - rVtxZ10Cent[runIndex]) /nTotal);
  
}

//-----------------------------------------------------------------------------
void QAFillOccupancy()
{
  TH2D *hCellNXZM1 = (TH2D*)listHist->FindObject("hCellNXZM1");
  TH2D *hCellNXZM2 = (TH2D*)listHist->FindObject("hCellNXZM2");
  TH2D *hCellNXZM3 = (TH2D*)listHist->FindObject("hCellNXZM3");
  
  // Count cells in each module
  
  Int_t nCells1=0, nCells2=0, nCells3=0;
  for (Int_t cellX=1; cellX<=64; cellX++) {
    for (Int_t cellZ=1; cellZ<=56; cellZ++) {
      if (hCellNXZM1->GetBinContent(cellX,cellZ) > 0) nCells1++;
      if (hCellNXZM2->GetBinContent(cellX,cellZ) > 0) nCells2++;
      if (hCellNXZM3->GetBinContent(cellX,cellZ) > 0) nCells3++;
    }
  }
  nCellsM1[runIndex] = nCells1;
  nCellsM2[runIndex] = nCells2;
  nCellsM3[runIndex] = nCells3;
}

//-----------------------------------------------------------------------------
void QAFillClusters()
{
  TH1F *hTotSelEvents = (TH1F*)listHist->FindObject("hTotSelEvents");
  Double_t nEvents4 = hTotSelEvents->GetBinContent(4);
  if( ! nEvents4 ) return;

  TH2F *hCluEvsClu = (TH2F*)listHist->FindObject("hCluEvsClu");
  TH1D *hClusterE = hCluEvsClu->ProjectionX("hClusterE",4,41);

  cluEnergy[runIndex]  = hClusterE->GetMean();
  ecluEnergy[runIndex] = hClusterE->GetMeanError();
  
  nCluEvent[runIndex]  = hClusterE->Integral() / nEvents4;
  eCluEvent[runIndex]  = TMath::Sqrt( hClusterE->Integral() )/nEvents4;


  for(int cent = 0; cent < kNCents; ++cent) {
    for(int ipid = 0; ipid < kNPID; ++ipid) {
      TH1* hPhot = listHist->FindObject( Form("hPhot%s_cen%d", kPIDNames[ipid], cent) );

      hPhot->SetAxisRange(0., 100.);
      double nPhot = hPhot->Integral() /nEvents4;
      nPhotPID [cent][ipid][runIndex] = nPhot;
      enPhotPID[cent][ipid][runIndex] = TMath::Sqrt( nPhot /nEvents4 );
      mEnPID   [cent][ipid][runIndex] = hPhot->GetMean();
      emEnPID  [cent][ipid][runIndex] = hPhot->GetMeanError();

      hPhot->SetAxisRange(1., 100.);
      double nPhotHigh = hPhot->Integral() /nEvents4;
      nPhotPIDHigh [cent][ipid][runIndex] = nPhotHigh;
      enPhotPIDHigh[cent][ipid][runIndex] = TMath::Sqrt( nPhotHigh /nEvents4 );
      mEnPIDHigh   [cent][ipid][runIndex] = hPhot->GetMean();
      emEnPIDHigh  [cent][ipid][runIndex] = hPhot->GetMeanError();

    }
  }
}

//-----------------------------------------------------------------------------
void QAFillTracks()
{
  TH2F *hCenTrack = (TH2F*)listHist->FindObject("hCenTrack");
  TH1D *hTrackMult = hCenTrack->ProjectionY("hTrackMult", 1, hCenTrack->GetNbinsY());

  nTracks0[runIndex] = hTrackMult->GetMean();
  eTracks0[runIndex] = hTrackMult->GetMeanError();
  
//   hTrackMult->DrawCopy();
//  hTrackMult->SetAxisRange(5., hTrackMult->GetXaxis()->GetXmax());
//   new TCanvas;
//   hTrackMult->DrawCopy();

//   nTracks1[runIndex] = hTrackMult->GetMean();
//   eTracks1[runIndex] = hTrackMult->GetMeanError();
}

//-----------------------------------------------------------------------------
void QAFillNPi0()
{
  TCanvas* canv = new TCanvas;
  canv->Divide(3,3);
  canv->cd(1);
  
  TH1 *hTotSelEvents = (TH1*)listHist->FindObject("hTotSelEvents");
  TH2F *hReMassPt = (TH2F*)listHist->FindObject("hPi0All_cen0");
  TH2 *hMiMassPt = (TH2*)listHist->FindObject("hMiPi0All_cen0");
  TH1* hReMass = Pi0InvMass(hReMassPt,2.,10.);
  TH1* hMiMass = Pi0InvMass(hMiMassPt,2.,10.);
  hReMass->Rebin(2);
  hMiMass->Rebin(2);
  hReMass->SetAxisRange(0.,0.3);
  hMiMass->SetAxisRange(0.,0.3);
  
  // Draw
  hReMassPt->DrawCopy("colz");
  canv->cd(2);
  hReMass->DrawCopy();
  canv->cd(3);
  
  
  TH1* hReMiRatio = (TH1*)hReMass->Clone("hReMiRatio");
  TH1* hPi0SubBG  = (TH1*)hReMass->Clone("hPi0SubBG");
  hReMiRatio->Divide(hMiMass);
  hReMiRatio->DrawCopy();
  canv->cd(4);
  

  TF1 * fitM = new TF1("fitM",pi0massP2,0.,1.,6) ;
  fitM->SetLineColor(kRed);
  fitM->SetParName(0,"A") ;
  fitM->SetParName(1,"m_{0}") ;
  fitM->SetParName(2,"#sigma") ;
  fitM->SetParName(3,"a_{0}") ;
  fitM->SetParName(4,"a_{1}") ;
  fitM->SetParName(5,"a_{2}") ;

  TF1 * fitG = new TF1("fitG",pi0mass,0.,1.,3) ;
  fitG->SetLineColor(kRed);
  fitG->SetParName(0,"A") ;
  fitG->SetParName(1,"m_{0}") ;
  fitG->SetParName(2,"#sigma") ;

  TF1 * fitBG = new TF1("fitBG",bgP2,0.,1.,3) ;
  fitBG->SetLineColor(kBlue);
  fitBG->SetLineStyle(2);
  fitBG->SetParName(0,"a_{0}") ;
  fitBG->SetParName(1,"a_{1}") ;
  fitBG->SetParName(2,"a_{2}") ;

  fitM->SetParameters(0.1,0.136,0.007,0.0013,-0.0007, 0.0) ;
  fitM->SetParLimits(0,0.000,1.000) ;
  fitM->SetParLimits(1,0.130,0.14) ;
  fitM->SetParLimits(2,0.005,0.012) ;

  Double_t rangeMin=0.06 ;
  Double_t rangeMax=0.25 ;
  hReMiRatio->Fit(fitM,"NQ","",rangeMin,rangeMax) ;
  hReMiRatio->Fit(fitM,"MQ","",rangeMin,rangeMax) ;
  fitBG->SetParameters(fitM->GetParameter(3),
		       fitM->GetParameter(4),
		       fitM->GetParameter(5)); 
  hMiMass->Multiply(fitBG) ;
  hPi0SubBG ->Add(hMiMass,-1.) ;


  fitG->SetParameters(100.,fitM->GetParameter(1),fitM->GetParameter(2)) ;
  fitG->SetParLimits(0,0.000,1.e+5) ;
  fitG->SetParLimits(1,0.120,0.145) ;
  fitG->SetParLimits(2,0.005,0.012) ;
  hPi0SubBG->Fit(fitG,"Q","",rangeMin,rangeMax);
  Double_t pi0Peak  = fitG->GetParameter(1);
  Double_t pi0Sigma = fitG->GetParameter(2);
  Double_t epi0Peak  = fitG->GetParError(1);
  Double_t epi0Sigma = fitG->GetParError(2);
  Int_t iMin = hPi0SubBG->FindBin(pi0Peak - 3*pi0Sigma);
  Int_t iMax = hPi0SubBG->FindBin(pi0Peak + 3*pi0Sigma);

  Double_t nPi0,ePi0;
  
  nPi0 = hPi0SubBG->Integral(iMin,iMax);
  ePi0 = TMath::Sqrt(nPi0);
  // Fit1Pi0(hMass,2,0.05,0.20,nPi0,ePi0); 
  
  Double_t nEvents4 = hTotSelEvents->GetBinContent(4);
  mPi0     [runIndex] = pi0Peak;
  emPi0    [runIndex] = epi0Peak;
  wPi0     [runIndex] = pi0Sigma;
  ewPi0    [runIndex] = epi0Sigma;
  nPi0Event[runIndex] = nPi0/nEvents4;
  ePi0Event[runIndex] = ePi0/nEvents4;

  // printf("Npi0(%d,%d) = %.1f\n",iMin,iMax,nPi0);
  // TCanvas *c1 = new TCanvas("c1","c1",0,0,800,600);
  // c1->Divide(2,2);
  // c1->cd(1);
  // hReMass->DrawCopy();
  // c1->cd(2);
  // hMiMass->DrawCopy();
  // c1->cd(3);
  // hReMiRatio->DrawCopy();
  // fitBG->Draw("same");
  // c1->cd(4);
  // hPi0SubBG->DrawCopy();

}

//-----------------------------------------------------------------------------
void QAWriteEventSelection()
{
  TH1F *hNSelected       = new TH1F("hNSelected","N_{selected}",runIndex,0,runIndex);
  TH1F *grVtx           = new TH1F("grVtx","N_{vertex exists} / N_{total}",runIndex,0,runIndex);
  TH1F *grVtxZ10        = new TH1F("grVtxZ10","N_{vertex exists, |Z|<10 cm} / N_{total}",runIndex,0,runIndex);
  TH1F *grVtxZ10Cent    = new TH1F("grVtxZ10Cent","N_{vertex exists, |Z|<10 cm, centrality} / N_{total}",runIndex,0,runIndex);

  for (Int_t i=0; i<runIndex; i++) {
    hNSelected->SetBinContent(i+1, nTotal4[i]);

    grVtx          ->SetBinContent(i+1,rVtx[i]);
    grVtxZ10       ->SetBinContent(i+1,rVtxZ10[i]);
    grVtxZ10Cent   ->SetBinContent(i+1,rVtxZ10Cent[i]);

    grVtx          ->SetBinError(i+1,eVtx[i]);
    grVtxZ10       ->SetBinError(i+1,eVtxZ10[i]);
    grVtxZ10Cent   ->SetBinError(i+1,eVtxZ10Cent[i]);

    grVtx          ->GetXaxis()->SetBinLabel(i+1,Form("%d",runNum[i]));
    grVtxZ10       ->GetXaxis()->SetBinLabel(i+1,Form("%d",runNum[i]));
    grVtxZ10Cent   ->GetXaxis()->SetBinLabel(i+1,Form("%d",runNum[i]));
  }

  grVtx          ->LabelsOption("v");
  grVtxZ10       ->LabelsOption("v");
  grVtxZ10Cent   ->LabelsOption("v");

  grVtx          ->SetMarkerStyle(33);
  grVtxZ10       ->SetMarkerStyle(33);
  grVtxZ10Cent   ->SetMarkerStyle(33);

  hNSelected     ->Write();
  grVtx          ->Write();
  grVtxZ10       ->Write();
  grVtxZ10Cent   ->Write();
}

//-----------------------------------------------------------------------------
void QAWriteOccupancy()
{
  // Number of fired cells in each module

  TH1F *grNCellsM1      = new TH1F("grNCellsM1","N_{cells} in module 1",runIndex,0,runIndex);
  TH1F *grNCellsM2      = new TH1F("grNCellsM2","N_{cells} in module 2",runIndex,0,runIndex);
  TH1F *grNCellsM3      = new TH1F("grNCellsM3","N_{cells} in module 3",runIndex,0,runIndex);

  for (Int_t i=0; i<runIndex; i++) {
    grNCellsM1          ->SetBinContent(i+1,nCellsM1[i]);
    grNCellsM2          ->SetBinContent(i+1,nCellsM2[i]);
    grNCellsM3          ->SetBinContent(i+1,nCellsM3[i]);

    grNCellsM1          ->SetBinError(i+1,0);
    grNCellsM2          ->SetBinError(i+1,0);
    grNCellsM3          ->SetBinError(i+1,0);

    grNCellsM1          ->GetXaxis()->SetBinLabel(i+1,Form("%d",runNum[i]));
    grNCellsM2          ->GetXaxis()->SetBinLabel(i+1,Form("%d",runNum[i]));
    grNCellsM3          ->GetXaxis()->SetBinLabel(i+1,Form("%d",runNum[i]));
    
    grNCellsM1->GetYaxis()->SetTitle("N_{cells}");
    grNCellsM2->GetYaxis()->SetTitle("N_{cells}");
    grNCellsM3->GetYaxis()->SetTitle("N_{cells}");
  }

  grNCellsM1     ->LabelsOption("v");
  grNCellsM2     ->LabelsOption("v");
  grNCellsM3     ->LabelsOption("v");
  grNCellsM1     ->SetMarkerStyle(33);
  grNCellsM2     ->SetMarkerStyle(33);
  grNCellsM3     ->SetMarkerStyle(33);
  grNCellsM1     ->Write();
  grNCellsM2     ->Write();
  grNCellsM3     ->Write();
}

//-----------------------------------------------------------------------------
void QAWriteClusters()
{
  TH1F *grECluster = new TH1F("grECluster","#LTE_{cluster}#GT",runIndex,0,runIndex);
  for (Int_t i=0; i<runIndex; i++) {
    grECluster          ->SetBinContent(i+1,cluEnergy[i]);
    grECluster          ->SetBinError(i+1,ecluEnergy[i]);
    grECluster          ->GetXaxis()->SetBinLabel(i+1,Form("%d",runNum[i]));
  }

  grECluster ->LabelsOption("v");
  grECluster ->SetMarkerStyle(33);
  grECluster ->Write();

  TH1F *grNCluster = new TH1F("grNCluster","#LTN_{clusters}#GT",runIndex,0,runIndex);
  for (Int_t i=0; i<runIndex; i++) {
    grNCluster          ->SetBinContent(i+1,nCluEvent[i]);
    grNCluster          ->SetBinError(i+1,eCluEvent[i]);
    grNCluster          ->GetXaxis()->SetBinLabel(i+1,Form("%d",runNum[i]));
  }
  grNCluster ->LabelsOption("v");
  grNCluster ->SetMarkerStyle(33);
  grNCluster ->Write();

  for(int cent=0; cent<kNCents; ++cent) {
    for(int ipid = 0; ipid < kNPID; ++ipid) {
      TH1* hPhot = listHist->FindObject( Form("hPhot%s_cen%d", kPIDNames[ipid], cent) );
      AddWriteTH1F(Form("grNPhot%s_cen%d", kPIDNames[ipid], cent), Form("#LTN_{clusters}^{%s}#GT, c.bin=%d", kPIDNames[ipid], cent),
		   nPhotPID[cent][ipid], enPhotPID[cent][ipid]);
      AddWriteTH1F(Form("grEn%s_cen%d", kPIDNames[ipid], cent), Form("#LTE_{clusters}^{%s}#GT, c.bin=%d", kPIDNames[ipid], cent),
		   mEnPID[cent][ipid], emEnPID[cent][ipid]);
      AddWriteTH1F(Form("grNPhot%sHigh_cen%d", kPIDNames[ipid], cent), Form("#LTN_{clusters}^{%s}#GT, c.bin=%d", kPIDNames[ipid], cent),
		   nPhotPIDHigh[cent][ipid], enPhotPIDHigh[cent][ipid]);
      AddWriteTH1F(Form("grEn%sHigh_cen%d", kPIDNames[ipid], cent), Form("#LTE_{clusters}^{%s}#GT, c.bin=%d", kPIDNames[ipid], cent),
		   mEnPIDHigh[cent][ipid], emEnPIDHigh[cent][ipid]);
    }
  }
}

//-----------------------------------------------------------------------------
void QAWriteTracks()
{
  TH1F *grNTracks0 = new TH1F("grNTracks0","#LTN_{tracks}#GT",runIndex,0,runIndex);
  for (Int_t i=0; i<runIndex; i++) {
    grNTracks0          ->SetBinContent(i+1,nTracks0[i]);
    grNTracks0          ->SetBinError(i+1,eTracks0[i]);
    grNTracks0          ->GetXaxis()->SetBinLabel(i+1,Form("%d",runNum[i]));
  }
  grNTracks0 ->LabelsOption("v");
  grNTracks0 ->SetMarkerStyle(33);
  grNTracks0 ->Write();

//   TH1F *grNTracks1 = new TH1F("grNTracks1","#LTN_{tracks}#GT, N #geq 5",runIndex,0,runIndex);
//   for (Int_t i=0; i<runIndex; i++) {
//     grNTracks1          ->SetBinContent(i+1,nTracks0[i]);
//     grNTracks1          ->SetBinError(i+1,eTracks0[i]);
//     grNTracks1          ->GetXaxis()->SetBinLabel(i+1,Form("%d",runNum[i]));
//   }
//   grNTracks1 ->LabelsOption("v");
//   grNTracks1 ->SetMarkerStyle(33);
//   grNTracks1 ->Write();
}

//-----------------------------------------------------------------------------
void QAWriteNPi0()
{
  TH1F *grNPi0 = new TH1F("grNPi0","N(#pi^{0}) per event",runIndex,0,runIndex);
  for (Int_t i=0; i<runIndex; i++) {
    grNPi0          ->SetBinContent(i+1,nPi0Event[i]);
    grNPi0          ->SetBinError(i+1,ePi0Event[i]);
    grNPi0          ->GetXaxis()->SetBinLabel(i+1,Form("%d",runNum[i]));
  }
  grNPi0     ->LabelsOption("v");
  grNPi0     ->SetMarkerStyle(33);
  grNPi0     ->Write();

  TH1F *grMPi0 = new TH1F("grMPi0","M(#pi^{0})",runIndex,0,runIndex);
  for (Int_t i=0; i<runIndex; i++) {
    grMPi0          ->SetBinContent(i+1,mPi0[i]);
    grMPi0          ->SetBinError(i+1,emPi0[i]);
    grMPi0          ->GetXaxis()->SetBinLabel(i+1,Form("%d",runNum[i]));
  }
  grMPi0     ->LabelsOption("v");
  grMPi0     ->SetMarkerStyle(33);
  grMPi0     ->Write();

  TH1F *grWPi0 = new TH1F("grWPi0","#sigma(#pi^{0})",runIndex,0,runIndex);
  for (Int_t i=0; i<runIndex; i++) {
    grWPi0          ->SetBinContent(i+1,wPi0[i]);
    grWPi0          ->SetBinError(i+1,ewPi0[i]);
    grWPi0          ->GetXaxis()->SetBinLabel(i+1,Form("%d",runNum[i]));
  }
  grWPi0     ->LabelsOption("v");
  grWPi0     ->SetMarkerStyle(33);
  grWPi0     ->Write();
}

//-----------------------------------------------------------------------------
TH1D * Pi0InvMass(TH2 *hMassPt, Float_t ptmin, Float_t ptmax)
{
  Int_t bmin = hMassPt->GetYaxis()->FindBin(ptmin);
  Int_t bmax = hMassPt->GetYaxis()->FindBin(ptmax);
  
  TString name = hMassPt->GetName();
  name += "_M";
  TH1D* hMass = hMassPt->ProjectionX(name,bmin,bmax);
  char htitle[64];
  sprintf(htitle,"%.1f<p_{T}<%.1f GeV/c",ptmin,ptmax);
  hMass->SetTitle(htitle);
  
  return hMass;

}

//-----------------------------------------------------------------------------
TH1F *AddWriteTH1F(const char *name, const char *title, Double_t* values, Double_t* errors)
{
  TH1F* hist = new TH1F(name, title, runIndex, 0, runIndex);
  hist->SetContent(values);
  hist->SetError(errors);

  for(Int_t i=0; i<runIndex; i++) {
    hist->GetXaxis()->SetBinLabel( i+1, Form("%d",runNum[i]) );
  }

  hist->LabelsOption("v");
  hist->SetMarkerStyle(33);
  hist->Write();

  return hist;
}

//-----------------------------------------------------------------------------
void Fit1Pi0(TH1D *hMass,
	     const Int_t polN,
	     const Double_t mFitMin,
	     const Double_t mFitMax,
	     Double_t &nPi0, Double_t &ePi0)
{
  // This script takes a 2D histogram hMassPt with invariant mass vs
  // pt of cluster pairs:
  // - slices them along pt bins,
  // - fits 1D invariant mass specrta by Gaussian+polynomial
  // - calculates the number of pi0's as an integral of Gaussian
  // - calculates background under pi0 as an integral of polynomial
  // Author: Yuri Kharlov.

  TF1 *fitfun = 0;
  if      (polN == 1)
    fitfun = new TF1("fitfun",pi0massP1,0.1,0.7,5);
  else if (polN == 2)
    fitfun = new TF1("fitfun",pi0massP2,0.1,0.7,6);
  fitfun->SetLineColor(kRed);
  fitfun->SetLineWidth(2);

  fitfun->SetParName(0,"A");
  fitfun->SetParName(1,"m_{0}");
  fitfun->SetParName(2,"#sigma");
  fitfun->SetParName(3,"a_{0}");
  fitfun->SetParName(4,"a_{1}");
  if (polN == 2)
    fitfun->SetParName(5,"a_{2}");
  
  fitfun->SetParLimits(0,  1.000,1.e+5);
  fitfun->SetParLimits(1,  0.09,0.18);
  fitfun->SetParLimits(2,  0.003,0.020);

  Int_t   nM     = hMass->GetNbinsX();
  Float_t mMin   = hMass->GetXaxis()->GetXmin();
  Float_t mMax   = hMass->GetXaxis()->GetXmax();
  Float_t mStep  = (mMax-mMin)/nM;

  hMass->SetXTitle("M_{#gamma#gamma}, GeV/c");
  hMass->SetAxisRange(0.01,1.0);
  Float_t nPi0Max = hMass->Integral(hMass->FindBin(0.30),
                                    hMass->FindBin(0.40));
  for (Int_t iM=1; iM<nM; iM++)
    if (hMass->GetBinContent(iM) == 0) hMass->SetBinError(iM,1.0);
  hMass->SetMinimum(0.01);
  if      (polN == 1)
    fitfun->SetParameters(nPi0Max/4,0.135,0.010,1.,0.);
  else if (polN == 2)
    fitfun->SetParameters(nPi0Max/4,0.135,0.010,1.,0.,0.);
  hMass->Fit("fitfun","Q0","",mFitMin,mFitMax);

  nPi0 = fitfun->GetParameter(0)*fitfun->GetParameter(2)/
    mStep*TMath::Sqrt(TMath::TwoPi());
  ePi0 = TMath::Sqrt(nPi0);
}


//-----------------------------------------------------------------------------
Double_t pi0massP2(Double_t *x, Double_t *par)
{
  Double_t gaus = par[0] * TMath::Exp( -(x[0]-par[1])*(x[0]-par[1]) /
                                       (2*par[2]*par[2]) );
  Double_t back = par[3] + par[4]*x[0] + par[5]*x[0]*x[0];
  return gaus+back;
}

//-----------------------------------------------------------------------------
Double_t pi0massP1(Double_t *x, Double_t *par)
{
  Double_t gaus = par[0] * TMath::Exp( -(x[0]-par[1])*(x[0]-par[1]) /
                                       (2*par[2]*par[2]) );
  Double_t back = par[3] + par[4]*x[0];
  return gaus+back;
}

//-----------------------------------------------------------------------------
Double_t pi0mass(Double_t *x, Double_t *par)
{
  Double_t gaus = par[0] * TMath::Exp( -(x[0]-par[1])*(x[0]-par[1]) /
                                       (2*par[2]*par[2]) );
  return gaus;
}

//-----------------------------------------------------------------------------
Double_t bgP2(Double_t *x, Double_t *par)
{
  Double_t back = par[0] + par[1]*x[0] + par[2]*x[0]*x[0];
  return back;
}
