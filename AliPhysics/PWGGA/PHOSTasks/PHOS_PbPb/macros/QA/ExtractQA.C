#include <TCanvas.h>

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TGrid.h>
#include <TStyle.h>
#include <TRandom.h>
#include <TFile.h>
#include <TF1.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TList.h>
#include <TMath.h>
#include <TGraphErrors.h>
#include <TString.h>
#include "TFitResult.h"
#include "THashList.h"
#include "Riostream.h"
#include "stdio.h"
using namespace std;
#endif

//-----------------------------------------------------------------------------
// Function declaration
void QAFillEventSelection();
void QAFillOccupancy();
void QAFillClusters();
void QAFillRP();
void QAFillTracks();
void QAFillNPi0();

void QAWriteEventSelection();
void QAWriteOccupancy();
void QAWriteClusters();
void QAWriteRP();
void QAWriteTracks();
void QAWriteNPi0();

TH1D *Pi0InvMass(TH2 *hMassPt, Float_t ptmin, Float_t ptmax);
TH1F *AddWriteTH1F(const char *name, const char *title, Double_t* values, Double_t* errors);
const TF1* GetPeriodRPFit(const char* histName = "phiRP");
void AddNoise(TH1D* hist, double noise);
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
const Int_t kNEventsBin = 4;
const Int_t kNCents = 1;
const Int_t kNPID = 8+6;
const char* kPIDNames[kNPID] = {"All", "Allwou", "Disp", "Disp2", "Dispwou", "CPV", "CPV2", "Both",
				"Allcore", "Dispcore", "CPVcore", "Bothcore", "Both2core", "Disp2core"};
const char* fullMergeFileName = "AnalysisResults.root";

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
const float highPt = 1.;
Double_t nPhotPIDHigh[kNCents][kNPID][500], enPhotPIDHigh[kNCents][kNPID][500];
Double_t mEnPIDHigh[kNCents][kNPID][500],   emEnPIDHigh[kNCents][kNPID][500];
// Reaction Plane
const int nRPD = 3;// Reaction Plane Detector
enum RPD { V0A=0, V0C=1, TPC=2};
const char* rpdNames[nRPD] = {"V0A", "V0C", "TPC" };
const char* rpdId[nRPD] = {"V0A", "V0C", "" };
const char* flatId[2] = {"", "flat"};
//Double_t 
Double_t rpSE[nRPD][2][2][500];// [rpd][flat][val/err][run]
Double_t rpSin[nRPD][2][4][2][500];// [rpd][flat][parameter][val/err][run]
Double_t rpChi2[nRPD][2][2][500];// [rpd][flat][chi val/err][run]
// Tracks:
Double_t nTracks0[500], eTracks0[500];
Double_t nTracks1[500], eTracks1[500];
// Pi0:
Double_t mPi0[500], emPi0[500];
Double_t wPi0[500], ewPi0[500];
Double_t nPi0Event[500], ePi0Event[500];


//-----------------------------------------------------------------------------
void ExtractQA(const TString runFile="runFile.txt",
	       const TString outputFileName = "outputQA.root"
	       )
{
  // Loop over per-run root files and run various QA functions

  //TGrid::Connect("alien://");

  ifstream in;
  in.open(runFile.Data());

  TFile *rootFile;

  runIndex = 0;
  while (1) {
    char rootFileName[256];
    Int_t runNumber=0;
    in >> runNumber >> rootFileName;
    if (!in.good()) break;

    runNum[runIndex] = runNumber;

    printf("root file is %s, run # = %d\n",rootFileName,runNumber);
    // char *runNum = strtok(rootFileName+35,".");
    rootFile = TFile::Open(rootFileName,"read");
    listHist = (TList*)rootFile->Get("PHOSPi0Flow_kMB/PHOSPi0Flow_kMBCoutput1");

    run[runIndex]            = runIndex+1;
    QAFillEventSelection();
    QAFillOccupancy();
    QAFillClusters();
    QAFillRP();
    QAFillTracks();
    QAFillNPi0();

    listHist->Clear();
    rootFile->Close();
    runIndex++;
//     if(runIndex>3)
//       break;
  }


  TFile *fileQA = TFile::Open(outputFileName.Data(), "recreate");
  QAWriteEventSelection();
  QAWriteOccupancy();
  QAWriteClusters();
  QAWriteRP();
  QAWriteTracks();
  QAWriteNPi0();
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
      TObject* obj = listHist->FindObject( Form("hPhot%s_cen%d", kPIDNames[ipid], cent) );
      TH1* hPhot = dynamic_cast<TH1*> ( obj );

      hPhot->SetAxisRange(0., 100.);
      double nPhot = hPhot->Integral() /nEvents4;
      nPhotPID [cent][ipid][runIndex] = nPhot;
      enPhotPID[cent][ipid][runIndex] = TMath::Sqrt( nPhot /nEvents4 );
      mEnPID   [cent][ipid][runIndex] = hPhot->GetMean();
      emEnPID  [cent][ipid][runIndex] = hPhot->GetMeanError();

      hPhot->SetAxisRange(highPt, 100.);
      double nPhotHigh = hPhot->Integral() /nEvents4;
      nPhotPIDHigh [cent][ipid][runIndex] = nPhotHigh;
      enPhotPIDHigh[cent][ipid][runIndex] = TMath::Sqrt( nPhotHigh /nEvents4 );
      mEnPIDHigh   [cent][ipid][runIndex] = hPhot->GetMean();
      emEnPIDHigh  [cent][ipid][runIndex] = hPhot->GetMeanError();

    }
  }
}
//-----------------------------------------------------------------------------
void QAFillRP()
{
  //TH1F *hev         = (TH1F*)listHist->FindObject("hTotSelEvents") ;
  TH2F* phiRP       = (TH2F*)listHist->FindObject("phiRP");
  TH2F* phiRPV0A    = (TH2F*)listHist->FindObject("phiRPV0A");
  TH2F* phiRPV0C    = (TH2F*)listHist->FindObject("phiRPV0C");
  TH2F* phiRPV0Aflat= (TH2F*)listHist->FindObject("phiRPV0Aflat");
  TH2F* phiRPV0Cflat= (TH2F*)listHist->FindObject("phiRPV0Cflat");
  TH2F* phiRPflat   = (TH2F*)listHist->FindObject("phiRPflat");

  //int nEvents = hev->GetBinContent(kNEventsBin);

  TH1D* phiRP1[nRPD][2] = {{0}};
  phiRP1[V0A][0] = phiRPV0A->ProjectionX();
  phiRP1[V0C][0] = phiRPV0C->ProjectionX();
  phiRP1[TPC][0] = phiRP->ProjectionX();
  phiRP1[V0A][1] = phiRPV0Aflat->ProjectionX();
  phiRP1[V0C][1] = phiRPV0Cflat->ProjectionX();
  phiRP1[TPC][1] = phiRPflat->ProjectionX();

  for(int rpd=0; rpd<nRPD; ++rpd) {
    for(int flat=0; flat<2; ++flat) {
      TH1D* hist = phiRP1[rpd][flat];
      const int nEntries = hist->GetEntries();
      const int nBins = hist->GetNbinsX();
      if ( nEntries < 1000 ) {
	//Printf("skipping RP: %s, very low number of entries", hist->GetName());
	continue;
      }
      
      // rpSE
      const Double_t* array = & hist->GetArray()[1]; 
      Double_t mean = TMath::Mean(nBins, array);
      Double_t rms = TMath::RMS(nBins, array) / mean * TMath::Sqrt(double(nBins)/(nBins-1));
      rpSE[rpd][flat][0][runIndex] = rms; 
      rpSE[rpd][flat][1][runIndex] = 0;

      // rpSin
      const TF1* periodFunc = GetPeriodRPFit(Form("phiRP%s%s", rpdId[rpd], flatId[flat]));
      for(int param=1; param < 4; ++param) {
	TF1* func = (TF1*)periodFunc->Clone("newRPFunc");
	for(int param2 = 1; param2 < 4; ++param2)
	  if(param != param2 ) func->FixParameter(param2, periodFunc->GetParameter(param2));
	if(1 == param ) func->SetParLimits(1, 0, 1);
	if(2 == param ) func->SetParLimits(2, 0, 100);
	if(3 == param ) func->SetParLimits(3, periodFunc->GetParameter(3)-TMath::Pi(), periodFunc->GetParameter(3)+TMath::Pi());
	int error = hist->Fit(func, "0Q");
	if( error ) {
	  rpSin[rpd][flat][param][0][runIndex] = 0;
	  if(3 == param) rpSin[rpd][flat][3][0][runIndex] = -2*TMath::Pi();
	  rpSin[rpd][flat][param][1][runIndex] = 0;	  
	}
	else {
	  rpSin[rpd][flat][param][0][runIndex] = func->GetParameter(param);
	  if(3 == param) rpSin[rpd][flat][3][0][runIndex] = func->GetParameter(3) - periodFunc->GetParameter(3);
	  rpSin[rpd][flat][param][1][runIndex] = func->GetParError(param);
	}
	delete func;
      }

      // rpChi2
      if(flat)
	AddNoise(hist, 0.01); // stricter if flattened
      else
	AddNoise(hist, 0.1);
      TFitResultPtr result = hist->Fit("pol0", "NSQ");
      if( result.Get() && !int(result) ) {
	int ndf = result->Ndf();
	rpChi2[rpd][flat][0][runIndex] = result->Chi2()/ndf;
	rpChi2[rpd][flat][1][runIndex] = TMath::Sqrt(2*ndf)/ndf;
      } else {
	rpChi2[rpd][flat][0][runIndex] = 0;
	rpChi2[rpd][flat][1][runIndex] = 0;
      }
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
  TH1 *hTotSelEvents = (TH1*)listHist->FindObject("hTotSelEvents");
  Int_t nEvents4 = hTotSelEvents->GetBinContent(4);

  if( nEvents4 < 10000 ) {
    Printf(" -> skipping due to small number of selected events: %d", nEvents4);
    return;
  }


  TCanvas* canv = new TCanvas;
  canv->Divide(3,3);
  canv->cd(1);

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
  //hReMass->SetAxisRange(0.1, 0.2);
  hReMass->DrawCopy("E");
  canv->cd(3);


  TH1* hReMiRatio = (TH1*)hReMass->Clone("hReMiRatio");
  TH1* hPi0SubBG  = (TH1*)hReMass->Clone("hPi0SubBG");
  hReMiRatio->Divide(hMiMass);
  hReMiRatio->DrawCopy("E");
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

  fitM->SetParameters(0.003, 0.134, 0.007, 0.14, -0.02, -0.01) ;
  fitM->SetParLimits(0, 0., 1. ) ;
  fitM->SetParLimits(1, 0.1, 0.2 ) ;
  fitM->SetParLimits(2, 0., 0.05) ;

  Double_t rangeMin=0.06 ;
  Double_t rangeMax=0.25 ;
  TFitResultPtr mrp = hReMiRatio->Fit(fitM,"Q","",rangeMin,rangeMax) ;
  //mrp = hReMiRatio->Fit(fitM,"MQ","",rangeMin,rangeMax) ; // "M" always fail...
  int error = mrp;
  if( error % 1000) {
    Printf(" -> fit of fitM to hReMiRatio failed with error code %d", error % 1000);
    return;
  }
  else if( error )
    Printf("Warning: failure of 'improve result' of fit of fitM to hReMiRatio");
  fitBG->SetParameters(fitM->GetParameter(3),
		       fitM->GetParameter(4),
		       fitM->GetParameter(5));
  hReMiRatio->SetAxisRange(rangeMin, rangeMax);
  hReMiRatio->DrawCopy();
  canv->cd(5);


  hMiMass->Multiply(fitBG) ;
  hPi0SubBG ->Add(hMiMass,-1.) ;


  fitG->SetParameters(100.,fitM->GetParameter(1),fitM->GetParameter(2)) ;
  // fitG->SetParLimits(0,0.000,1.e+5) ;
  fitG->SetParLimits(1, rangeMin, rangeMax) ;
  fitG->SetParLimits(2, 0., rangeMax) ;
  mrp = hPi0SubBG->Fit(fitG,"Q","",rangeMin,rangeMax);
  if( (error=mrp) ) {
    Printf(" -> fit of fitG to hPi0SubBG failed with error code %d, skipping", error );
    return;
  }
  hPi0SubBG->SetAxisRange(rangeMin, rangeMax);
  hPi0SubBG->DrawCopy();



  Double_t pi0Peak  = fitG->GetParameter(1);
  Double_t pi0Sigma = fitG->GetParameter(2);
  Double_t epi0Peak  = fitG->GetParError(1);
  Double_t epi0Sigma = fitG->GetParError(2);
  Int_t iMin = hPi0SubBG->FindBin(pi0Peak - 3*pi0Sigma);
  Int_t iMax = hPi0SubBG->FindBin(pi0Peak + 3*pi0Sigma);

  Double_t nPi0,ePi0;


  nPi0 = hPi0SubBG->IntegralAndError(iMin, iMax, ePi0);
  //ePi0 = TMath::Sqrt(nPi0);
  // Fit1Pi0(hMass,2,0.05,0.20,nPi0,ePi0);

  mPi0     [runIndex] = pi0Peak;
  emPi0    [runIndex] = epi0Peak;
  wPi0     [runIndex] = pi0Sigma;
  ewPi0    [runIndex] = epi0Sigma;
  nPi0Event[runIndex] = nPi0/nEvents4;
  ePi0Event[runIndex] = ePi0/nEvents4;
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

  TString name, title;
  for(int cent=0; cent<kNCents; ++cent) {
    for(int ipid = 0; ipid < kNPID; ++ipid) {
      TObject* obj = listHist->FindObject( Form("hPhot%s_cen%d", kPIDNames[ipid], cent) );
      TH1* hPhot = dynamic_cast<TH1*> (obj);
      name = Form("grNPhot%s_cen%d", kPIDNames[ipid], cent);
      title = Form("#LTN_{clusters}^{%s}#GT, c.bin=%d", kPIDNames[ipid], cent);
      AddWriteTH1F(name, title, nPhotPID[cent][ipid], enPhotPID[cent][ipid]);
      name = Form("grEn%s_cen%d", kPIDNames[ipid], cent);
      title = Form("#LTE_{clusters}^{%s}#GT, c.bin=%d", kPIDNames[ipid], cent);
      AddWriteTH1F(name , title, mEnPID[cent][ipid], emEnPID[cent][ipid]);
      name = Form("grNPhot%sHigh_cen%d", kPIDNames[ipid], cent);
      title = Form("#LTN_{clusters}^{%s}#GT, c.bin=%d, p_T>%f", kPIDNames[ipid], cent, highPt);
      AddWriteTH1F(name, title, nPhotPIDHigh[cent][ipid], enPhotPIDHigh[cent][ipid]);
      name = Form("grEn%sHigh_cen%d", kPIDNames[ipid], cent);
      title = Form("#LTE_{clusters}^{%s}#GT, c.bin=%d, p_T>%f", kPIDNames[ipid], cent, highPt);
      AddWriteTH1F(name, title, mEnPIDHigh[cent][ipid], emEnPIDHigh[cent][ipid]);
    }
  }
}

void QAWriteRP()
{
  for(int rpd=0; rpd<nRPD; ++rpd) {
    for(int flat=0; flat<2; ++flat) {
      TString name;
      TString title;
      
      //rpSE
      name = Form("grSERP%s%s", rpdNames[rpd], flatId[flat]);
      title = Form("RMS (SE) of RP of %s, %s", rpdNames[rpd], flatId[flat]);
      TH1F* grSE = new TH1F(name.Data(), title.Data(), runIndex, 0, runIndex);
      for (Int_t i=0; i<runIndex; i++) {
	grSE          ->SetBinContent(i+1,  rpSE[rpd][flat][0][i]);
	grSE          ->SetBinError(i+1,    rpSE[rpd][flat][1][i]);
	grSE          ->GetXaxis()->SetBinLabel(i+1,Form("%d",runNum[i]));
      }
      grSE->GetYaxis()->SetTitle("Std. Err.");
      grSE->LabelsOption("v");
      grSE->SetMarkerStyle(33);
      grSE->Write();
      
      // rpSin
      for(int param=1; param < 4; ++param) {
	if( flat ) {
	  name = Form("grSinRP%sflat%d", rpdNames[rpd], param);
	  title = Form("Parameter [%d] of fit constraining others, flat, %s ", param, rpdNames[rpd]);
	} else {
	  name = Form("grSinRP%s%d", rpdNames[rpd], param);
	  title = Form("Parameter [%d] of fit constraining others, %s ", param, rpdNames[rpd]);
	}
	TH1F* grSin = new TH1F(name, title, runIndex, 0, runIndex);
	Printf(name.Data());
	for (Int_t i=0; i<runIndex; i++) {
	  //Printf("  %d chi2/ndf: %f, \terror:%f", i, rpSin[rpd][flat][i][0],rpSin[rpd][flat][i][1]);
	  grSin          ->SetBinContent(i+1,  rpSin[rpd][flat][param][0][i]);
	  grSin          ->SetBinError(i+1,    rpSin[rpd][flat][param][1][i]);
	  grSin          ->GetXaxis()->SetBinLabel(i+1,Form("%d",runNum[i]));
	}
	grSin->GetYaxis()->SetTitle(Form("p[%d]",param));
	if( 1 == param ) grSin->GetYaxis()->SetTitle("s_{1}");
	if( 2 == param ) grSin->GetYaxis()->SetTitle("#omega");
	if( 3 == param ) grSin->GetYaxis()->SetTitle("#psi - #psi_{P}");
	grSin->LabelsOption("v");
	grSin->SetMarkerStyle(33);
	grSin->Write();
      }
      // rpChi2
      if( flat ) {
	name = Form("grChi2RP%sflat", rpdNames[rpd]);
	title = Form("#Chi^2 / ndf to fit of pol0 of flattened RP of %s", rpdNames[rpd]);
      } else {
	name = Form("grChi2RP%s", rpdNames[rpd]);
	title = Form("#Chi^2 / ndf to fit of pol0 of raw RP of %s", rpdNames[rpd]);
      }
      TH1F* grChi2 = new TH1F(name, title, runIndex, 0, runIndex);
      Printf(name.Data());
      for (Int_t i=0; i<runIndex; i++) {
	//Printf("  %d chi2/ndf: %f, \terror:%f", i, rpChi2[rpd][flat][i][0],rpChi2[rpd][flat][i][1]);
	grChi2          ->SetBinContent(i+1,  rpChi2[rpd][flat][0][i]);
	grChi2          ->SetBinError(i+1,    rpChi2[rpd][flat][1][i]);
	grChi2          ->GetXaxis()->SetBinLabel(i+1,Form("%d",runNum[i]));
      }
      grChi2->GetYaxis()->SetTitle("#Chi^2 / ndf");
      grChi2->LabelsOption("v");
      grChi2->SetMarkerStyle(33);
      grChi2->Write();
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
  grNPi0->GetYaxis()->SetTitle("#LTN#GT");
  grNPi0     ->Write();

  TH1F *grMPi0 = new TH1F("grMPi0","M(#pi^{0})",runIndex,0,runIndex);
  for (Int_t i=0; i<runIndex; i++) {
    grMPi0          ->SetBinContent(i+1,mPi0[i]);
    grMPi0          ->SetBinError(i+1,emPi0[i]);
    grMPi0          ->GetXaxis()->SetBinLabel(i+1,Form("%d",runNum[i]));
  }
  grMPi0     ->LabelsOption("v");
  grMPi0     ->SetMarkerStyle(33);
  grMPi0->GetYaxis()->SetTitle("#LTM#GT");
  grMPi0     ->Write();

  TH1F *grWPi0 = new TH1F("grWPi0","#sigma(#pi^{0})",runIndex,0,runIndex);
  for (Int_t i=0; i<runIndex; i++) {
    grWPi0          ->SetBinContent(i+1,wPi0[i]);
    grWPi0          ->SetBinError(i+1,ewPi0[i]);
    grWPi0          ->GetXaxis()->SetBinLabel(i+1,Form("%d",runNum[i]));
  }
  grWPi0     ->LabelsOption("v");
  grWPi0     ->SetMarkerStyle(33);
  grWPi0->GetYaxis()->SetTitle("#LT#sigma#GT");
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

  hMass->Sumw2();

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

  if( TString(name).Contains("En") )
    hist->GetYaxis()->SetTitle("#LTE#GT");
  if( TString(name).Contains("NPhot") )
    hist->GetYaxis()->SetTitle("#LTN#GT");

  return hist;
}

const TF1* GetPeriodRPFit(const char* histName)
{
  TString name = Form("f1RPFit%s", histName);
  
  // if created earlier, find in list
  static THashList* list = new THashList;
  TF1* func = (TF1*) list->FindObject(name.Data());
  if( func ) return func;

  TFile* file = TFile::Open(fullMergeFileName, "read");
  TList* hlist = (TList*) file->Get("PHOSPi0Flow_kCentral/PHOSPi0Flow_kCentralCoutput1");
  TH2F* hist2d = (TH2F*) hlist->FindObject(histName);
  TH1D* hist = hist2d->ProjectionX();
  int nBins = hist->GetNbinsX();
  const Double_t* array = & hist->GetArray()[1]; 
  
  //TCanvas* canv = new TCanvas;
  func = new TF1(name.Data(), "[0]*(1. + [1]*sin([2]*x + [3]))", 0, TMath::Pi());
  double mean = TMath::Mean(nBins, array);
  double rms = TMath::RMS(nBins, array) / mean;
  func->SetParameters(mean, rms, 0.6*TMath::Pi(), -0.2*TMath::Pi());
  func->SetParLimits(0, 0, 2*mean);
  func->SetParLimits(1, 0, 1);
  func->SetParLimits(2, 0, 100);
  func->SetParLimits(3, -TMath::Pi(), TMath::Pi());
  func->SetParNames("s_{0}", "s_{1}", "#omega", "#psi");
  hist->GetXaxis()->SetTitle("#phi");
  TCanvas* canv = new TCanvas(name.Data());
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  int error = hist->Fit(func, "Q");
  if( ! name.Contains("flat") ) {
    gStyle->SetOptFit(1);
  } else {
    gStyle->SetOptFit(1);
  }
  hist->Draw("E");
  canv->SaveAs(Form("imgs/%s.pdf", name.Data()));
  if( error )
    Printf("GetPeriodRPFit::ERROR, could not fit %s, error:%d", histName, error);

  //hist->Draw();
  //func->Draw("same");
  //canv->SaveAs("rpfit.pdf");
  
  // Fit
  list->Add(func);
  return func;
}


//-----------------------------------------------------------------------------
void AddNoise(TH1D* hist, double noise)
{
  int nBins = hist->GetNbinsX();
  double mean = TMath::Mean(nBins, &((hist->GetArray())[1]) );
  //Printf("Adding noise to: %s", hist->GetName());
  for(int bin=1; bin<=nBins; ++bin) {
    hist->SetBinContent( bin, hist->GetBinContent(bin) + gRandom->Gaus(0, noise) );
    double err = hist->GetBinError(bin);
    hist->SetBinError( bin, TMath::Sqrt( err*err + noise*noise ) );
  }
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
