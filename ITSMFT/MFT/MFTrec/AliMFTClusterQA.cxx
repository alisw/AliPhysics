#include "TObject.h"
#include "AliRunLoader.h"
#include "AliRun.h"
#include "AliLoader.h"
#include "AliMFT.h"
#include "TClonesArray.h"
#include "AliMFTCluster.h"
#include "AliMFTSegmentation.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "AliLog.h"
#include "TString.h"
#include "TTree.h"

#include "AliMFTClusterQA.h"

//====================================================================================================================================================
//
// Class for the analysis of the MFT clusters (a.k.a. rec points). Few QA histograms are created
//
// Contact author: antonio.uras@cern.ch
//
//====================================================================================================================================================

ClassImp(AliMFTClusterQA)

//====================================================================================================================================================

AliMFTClusterQA::AliMFTClusterQA():
  TObject(),
  fMFTLoader(0),
  fRunLoader(0),
  fMFT(0),
  fNEvents(0),
  fEv(0),
  fFileOut(0),
  fReadDir(0),
  fOutDir(0)
{
  
  // default constructor

  for (Int_t iPlane=0; iPlane<AliMFTConstants::kNDisks; iPlane++) {
    fHistNClustersPerEvent[iPlane]     = 0;
    fHistNPixelsPerCluster[iPlane]     = 0;
    fHistClusterSizeX[iPlane]          = 0; 
    fHistClusterSizeY[iPlane]          = 0;
    fHistClusterRadialPosition[iPlane] = 0;
    fClusterScatterPlotXY[iPlane]      = 0;
  }

}

//====================================================================================================================================================

void AliMFTClusterQA::Init(Char_t *readDir, Char_t *outDir, Int_t nEventsToAnalyze) {

  fReadDir = readDir;
  fOutDir  = outDir;

  fRunLoader = AliRunLoader::Open(Form("%s/galice.root", fReadDir.Data()));
  gAlice = fRunLoader->GetAliRun();
  if (!gAlice) fRunLoader->LoadgAlice();
  fNEvents = fRunLoader->GetNumberOfEvents();
  if (nEventsToAnalyze>0) fNEvents = TMath::Min(fNEvents, nEventsToAnalyze);

  fMFT = (AliMFT*) gAlice->GetDetector("MFT"); 

  BookHistos();

  fMFTLoader = fRunLoader->GetDetectorLoader("MFT");
  fMFTLoader -> LoadRecPoints("READ");

}

//====================================================================================================================================================

Bool_t AliMFTClusterQA::LoadNextEvent() {

  if (fEv>=fNEvents) return kFALSE;
  AliDebug(1, Form("event %5d",fEv));
  
  fRunLoader->GetEvent(fEv);
  fEv++;

  if (!fMFTLoader->TreeR()->GetEvent()) return kTRUE;

  for (Int_t iPlane=0; iPlane<AliMFTConstants::kNDisks; iPlane++) {
    Int_t nClusters = fMFT->GetRecPointsList(iPlane)->GetEntries();
    fHistNClustersPerEvent[iPlane] -> Fill(nClusters);
    fClusterScatterPlotXY[iPlane]  -> Fill(0., 0.);    // "scaler" bin
    AliDebug(1,Form("nClusters = %5d", nClusters));
    for (Int_t iCluster=0; iCluster<nClusters; iCluster++) {
      AliMFTCluster *cluster = (AliMFTCluster*) (fMFT->GetRecPointsList(iPlane))->At(iCluster);
      fHistNPixelsPerCluster[iPlane] -> Fill(cluster->GetSize());
      fHistClusterSizeX[iPlane]      -> Fill(cluster->GetErrX()*1.e4);   // converted in microns
      fHistClusterSizeY[iPlane]      -> Fill(cluster->GetErrY()*1.e4);   // converted in microns
      fHistClusterRadialPosition[iPlane] -> Fill(TMath::Sqrt(cluster->GetX()*cluster->GetX()+cluster->GetY()*cluster->GetY()));
      fClusterScatterPlotXY[iPlane]      -> Fill(cluster->GetX(), cluster->GetY());
    }
  }

  return kTRUE;

}  

//====================================================================================================================================================

void AliMFTClusterQA::BookHistos() {

  fFileOut = new TFile(Form("%s/MFT.RecPoints.QA.root",fOutDir.Data()), "recreate");

  for (Int_t iPlane=0; iPlane<AliMFTConstants::kNDisks; iPlane++) {

    fHistNClustersPerEvent[iPlane] = new TH1D(Form("fHistNClustersPerEvent_Pl%02d",iPlane), 
					      Form("Number of clusters per event in Plane%02d",iPlane),
					      25000, -0.5, 24999.5);

    fHistNPixelsPerCluster[iPlane] = new TH1D(Form("fHistNPixelsPerCluster_Pl%02d",iPlane), 
					      Form("Number of pixels per cluster in Plane%02d",iPlane),
					      15, -0.5, 14.5);

    fHistClusterSizeX[iPlane]      = new TH1D(Form("fHistClusterSizeX_Pl%02d",iPlane), 
					      Form("#Deltax for clusters in Plane%02d",iPlane), 
					      100, 0., 100.);
    
    fHistClusterSizeY[iPlane]      = new TH1D(Form("fHistClusterSizeY_Pl%02d",iPlane), 
					      Form("#Deltay for clusters in Plane%02d",iPlane), 
					      100, 0., 100.);
    
    fHistNClustersPerEvent[iPlane] -> SetXTitle("N_{clusters} per Event");
    fHistNPixelsPerCluster[iPlane] -> SetXTitle("N_{pixels} per Cluster");
    fHistClusterSizeX[iPlane]      -> SetXTitle("#Deltax  [#mum]");
    fHistClusterSizeY[iPlane]      -> SetXTitle("#Deltay  [#mum]");

    fHistNClustersPerEvent[iPlane] -> Sumw2();
    fHistNPixelsPerCluster[iPlane] -> Sumw2();
    fHistClusterSizeX[iPlane]      -> Sumw2();
    fHistClusterSizeY[iPlane]      -> Sumw2();

    //------------------------------------------------------------

    Int_t rMax = 20.;//

    fHistClusterRadialPosition[iPlane] = new TH1D(Form("fHistClusterRadialPosition_Pl%02d",iPlane),
						  Form("Cluster radial position (Plane%02d)",iPlane),
						  rMax, 0, Double_t(rMax)/10.);

    fClusterScatterPlotXY[iPlane] = new TH2D(Form("fClusterScatterPlotXY_Pl%02d",iPlane),
					     Form("Cluster scatter plot (Plane%02d)",iPlane),
					     2*rMax+1, (-rMax-0.5)/10., (rMax+0.5)/10., 2*rMax+1, (-rMax-0.5)/10., (rMax+0.5)/10.);
    
    fHistClusterRadialPosition[iPlane] -> SetXTitle("R  [cm]");
    fClusterScatterPlotXY[iPlane]      -> SetXTitle("X  [cm]");
    fClusterScatterPlotXY[iPlane]      -> SetYTitle("Y  [cm]");

    fHistClusterRadialPosition[iPlane] -> Sumw2();
    fClusterScatterPlotXY[iPlane]      -> Sumw2();
    
  }
  
}

//====================================================================================================================================================

void AliMFTClusterQA::Terminate() {

  AliInfo("Writing QA histos...");

  // ---- equalize radial clusters distribution ----------------------

  for (Int_t iPlane=0; iPlane<AliMFTConstants::kNDisks; iPlane++) {
    for (Int_t iBin=1; iBin<=fHistClusterRadialPosition[iPlane]->GetNbinsX(); iBin++) {
      Double_t rMin = fHistClusterRadialPosition[iPlane]->GetBinLowEdge(iBin);        // in cm
      Double_t rMax = fHistClusterRadialPosition[iPlane]->GetBinWidth(iBin) + rMin;   // in cm
      Double_t area = 100.*TMath::Pi()*(rMax*rMax - rMin*rMin);                       // in mm^2
      Double_t density = fHistClusterRadialPosition[iPlane]->GetBinContent(iBin)/area;
      fHistClusterRadialPosition[iPlane]->SetBinContent(iBin, density);
      fHistClusterRadialPosition[iPlane]->SetBinError(iBin, fHistClusterRadialPosition[iPlane]->GetBinError(iBin)/area);
    }
    fHistClusterRadialPosition[iPlane] -> SetBinContent(1, fEv);      // "scaler" bin
  }

  // -----------------------------------------------------------------

  fFileOut->cd();

  for (Int_t iPlane=0; iPlane<AliMFTConstants::kNDisks; iPlane++) {
    fHistNClustersPerEvent[iPlane]     -> Write();
    fHistNPixelsPerCluster[iPlane]     -> Write();
    fHistClusterSizeX[iPlane]          -> Write();
    fHistClusterSizeY[iPlane]          -> Write();
    fHistClusterRadialPosition[iPlane] -> Write();
    fClusterScatterPlotXY[iPlane]      -> Write();
  }

  fFileOut -> Close();

}

//====================================================================================================================================================
