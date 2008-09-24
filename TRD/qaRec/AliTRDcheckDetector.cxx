#include <TFile.h>
#include <TH1F.h>
#include <TGraph.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TProfile.h>
#include <TROOT.h>

#include "AliTRDcluster.h"
#include "AliTRDseedV1.h"
#include "AliTRDtrackV1.h"
#include "TTreeStream.h"

#include "AliTRDtrackInfo/AliTRDtrackInfo.h"
#include "AliTRDcheckDetector.h"

#include <cstdio>

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  Reconstruction QA                                                     //
//                                                                        //
//  Task doing basic checks for tracking and detector performance         //
//                                                                        //
//  Authors:                                                              //
//    Anton Andronic <A.Andronic@gsi.de>                                  //
//    Markus Fasel <M.Fasel@gsi.de>                                       //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

//_______________________________________________________
AliTRDcheckDetector::AliTRDcheckDetector():
  AliTRDrecoTask("detChecker", "Basic Detector Checker")
  ,fPHSdetector(0x0)
  ,fPHSsector(0x0)
  ,fQCLdetector(0x0)
  ,fQCLsector(0x0)
  ,fQTdetector(0x0)
  ,fQTsector(0x0)
{
  //
  // Default constructor
  //
}

//_______________________________________________________
AliTRDcheckDetector::~AliTRDcheckDetector(){
  //
  // Destructor
  // 
  if(fPHSdetector) delete fPHSdetector;
  if(fPHSsector) delete fPHSsector;
  if(fQCLdetector) delete fQCLdetector;
  if(fQCLsector) delete fQCLsector;
  if(fQTdetector) delete fQTdetector;
  if(fQTsector) delete fQTsector;
}

//_______________________________________________________
void AliTRDcheckDetector::CreateOutputObjects(){
  //
  // Create Output Objects
  //
  fContainer = new TObjArray();
  TH1F *hNtrks = new TH1F("hNtrks", "Number of Tracks per event", 100, 0, 100);
  TH1F *hNcls = new TH1F("hNcls", "Nr. of clusters per track", 30, 0, 180);
  TH1F *hNtls = new TH1F("hNtls", "Nr. tracklets per track", 10, 0, 10);
  TH1F *hNclTls = new TH1F("hNclTls","Mean Number of clusters per tracklet", 30, 0, 30);
  TH1F *hChi2 = new TH1F("hChi2", "Chi2", 200, 0, 20);
  TH1F *hChi2n = new TH1F("hChi2n", "Norm. Chi2 (tracklets)", 50, 0, 5);
  TH1F *hSM = new TH1F("hSM", "Track Counts in Supermodule", 18, 0., 18.);
  TH1F *hQcl = new TH1F("hQcl1", "Cluster charge", 200, -1200, 1200);
  TProfile *hPH = new TProfile("hPH", "Average PH", 31, -0.5, 30.5);
  TH1F * hQT = new TH1F("hQT", "Total Charge Deposit", 6000, 0, 6000);
  TH1F * hEventsTrigger = new TH1F("hEventsTrigger", "Trigger Class", 20, 0, 20);
  // Register Histograms
  fContainer->Add(hNtrks);
  fContainer->Add(hEventsTrigger);
  fContainer->Add(hNcls);
  fContainer->Add(hNtls);
  fContainer->Add(hNclTls);
  fContainer->Add(hChi2);
  fContainer->Add(hChi2n);
  fContainer->Add(hSM);
  fContainer->Add(hPH);
  fContainer->Add(hQcl);
  fContainer->Add(hQT);

  // Initialise the PHS, cluster charge and total charge deposit histograms for each detector
  fPHSdetector = new TObjArray();
  fPHSdetector->SetName("fPHSdetector");
  fQCLdetector = new TObjArray();
  fQCLdetector->SetName("fQCLdetector");
  fQTdetector = new TObjArray();
  fQTdetector->SetName("fQTdetector");
  for(Int_t idet = 0; idet < kNDetectors; idet++){
    fPHSdetector->Add(new TProfile(Form("hPHSdet%d", idet), Form("Average Pulse Height in Detector %d", idet), 31, -0.5, 30.5));
    fQCLdetector->Add(new TH1F(Form("hQd%d",idet),Form("Qd%d",idet), 100,0,5000));
    fQTdetector->Add(new TH1F(Form("hQTd%d", idet), Form("Qtotal%d", idet), 6000, 0, 6000));
  }
  fContainer->Add(fPHSdetector);
  fContainer->Add(fQCLdetector);
  fContainer->Add(fQTdetector);
  
  // Initialise the PHS, cluster charge and total charge deposit histograms for each sector
  fPHSsector = new TObjArray();
  fPHSsector->SetName("fPHSsector");
  fQCLsector = new TObjArray();
  fQCLsector->SetName("fQCLsector");
  fQTsector = new TObjArray();
  fQTsector->SetName("fQTsector");
  for(Int_t isec = 0; isec < kNSectors; isec++){
    fPHSsector->Add(new TProfile(Form("hPHSsec%d", isec), Form("Average Pulse Height in Sector %d", isec), 31, -0.5, 30.5));
    fQCLsector->Add(new TH1F(Form("hQd%d",isec),Form("Qd%d",isec), 100,0,5000));
    fQTsector->Add(new TH1F(Form("hQTd%d", isec), Form("Qtotal%d", isec), 6000, 0, 6000));
  }
  fContainer->Add(fPHSsector);
  fContainer->Add(fQCLsector);
  fContainer->Add(fQTsector);
}

//_______________________________________________________
void AliTRDcheckDetector::Exec(Option_t *){
  //
  // Execution function
  // Filling TRD quality histos
  //
  Int_t nTracks = 0;		// Count the number of tracks per event
  Int_t triggermask = 0;
  AliTRDtrackInfo *fTrackInfo = 0x0;
  AliTRDtrackV1 *fTRDtrack = 0x0;
  AliTRDseedV1 *fTracklet = 0x0;
  AliTRDcluster *fTRDcluster = 0x0;
  for(Int_t iti = 0; iti < fTracks->GetEntriesFast(); iti++){
    fTrackInfo = dynamic_cast<AliTRDtrackInfo *>(fTracks->UncheckedAt(iti));
    if(!fTrackInfo || !(fTRDtrack = fTrackInfo->GetTRDtrack())) continue;
    if(!triggermask) triggermask = fTrackInfo->GetTriggerCluster();
    Int_t nclusters = fTRDtrack->GetNumberOfClusters();
    Int_t ntracklets = fTRDtrack->GetNumberOfTracklets();
    Float_t chi2 = fTRDtrack->GetChi2();
    // Fill Histograms
    dynamic_cast<TH1F *>(fContainer->UncheckedAt(kNclustersHist))->Fill(nclusters);
    dynamic_cast<TH1F *>(fContainer->UncheckedAt(kNtrackletsHist))->Fill(ntracklets);
    dynamic_cast<TH1F *>(fContainer->UncheckedAt(kChi2))->Fill(chi2);
    dynamic_cast<TH1F *>(fContainer->UncheckedAt(kChi2Normalized))->Fill(chi2/static_cast<Float_t>(ntracklets));
    // now loop over single tracklets
    if(ntracklets > 2){
      Int_t sector = -1;
      for(Int_t ilayer = 0; ilayer < kNLayers; ilayer++){
        if(!(fTracklet = fTRDtrack->GetTracklet(ilayer)) || !fTracklet->IsOK()) continue;
        Float_t nClustersTracklet = fTracklet->GetN();
        dynamic_cast<TH1F *>(fContainer->UncheckedAt(kNclusterTrackletHist))->Fill(nClustersTracklet);
        if(nClustersTracklet){
          Float_t Qtot = 0;
          Int_t detector = -1;
          for(Int_t itb = 0; itb < kNTimebins; itb++){
            if(!(fTRDcluster = fTracklet->GetClusters(itb))) continue;
            Int_t localtime        = fTRDcluster->GetLocalTimeBin();
            Double_t clusterCharge = fTRDcluster->GetQ();
            detector               = fTRDcluster->GetDetector();
            sector                 = static_cast<Int_t>(detector/kNDetectorsSector);
            Double_t absolute_charge = TMath::Abs(clusterCharge);
            Qtot += absolute_charge;
            // Fill Histograms
            // 1st. PHS
            dynamic_cast<TProfile *>(fContainer->UncheckedAt(kPulseHeight))->Fill(localtime, absolute_charge);
            dynamic_cast<TProfile *>(fPHSdetector->UncheckedAt(detector))->Fill(localtime, absolute_charge);
            dynamic_cast<TProfile *>(fPHSsector->UncheckedAt(sector))->Fill(localtime, absolute_charge);
            // 2nd. cluster charge
            dynamic_cast<TH1F *>(fContainer->UncheckedAt(kClusterCharge))->Fill(absolute_charge);
            dynamic_cast<TH1F *>(fQCLdetector->UncheckedAt(detector))->Fill(absolute_charge);
            dynamic_cast<TH1F *>(fQCLsector->UncheckedAt(sector))->Fill(absolute_charge);
          }
          // Fill total cluster charge histograms
          dynamic_cast<TH1F *>(fContainer->UncheckedAt(kChargeDeposit))->Fill(Qtot);
          dynamic_cast<TH1F *>(fQTdetector->UncheckedAt(detector))->Fill(Qtot);
          dynamic_cast<TH1F *>(fQTsector->UncheckedAt(sector))->Fill(Qtot);
        }
      }
      dynamic_cast<TH1F *>(fContainer->UncheckedAt(kNTracksSectorHist))->Fill(sector);
    }
    nTracks++;
  }
  dynamic_cast<TH1F *>(fContainer->UncheckedAt(kNTracksEventHist))->Fill(nTracks);
 	dynamic_cast<TH1F *>(fContainer->UncheckedAt(kNEventsTrigger))->Fill(triggermask);
  PostData(0, fContainer);
}

//_______________________________________________________
void AliTRDcheckDetector::Terminate(Option_t *){
  //
  // Terminate function
  //
}

//_______________________________________________________
Bool_t AliTRDcheckDetector::PostProcess(){
	//
	// Do Postprocessing (for the moment set the number of Reference histograms)
	//
	
	TH1 * histo = dynamic_cast<TH1F *>(fContainer->UncheckedAt(kNTracksEventHist));
	histo->GetXaxis()->SetTitle("Number of Tracks");
	histo->GetYaxis()->SetTitle("Events");
	// Calculate percentage of events containing tracks
	Int_t entries = (Int_t)histo->GetEntries();
  Float_t percTracks = entries > 0 ? (entries - histo->GetBinContent(1))/Float_t(entries) : 0.;
	printf("Percentage of Events containing TRD Tracks: %6.3f %%\n", 1.E2*percTracks);
	histo = dynamic_cast<TH1F *>(fContainer->UncheckedAt(kNclustersHist));
	histo->GetXaxis()->SetTitle("Number of Clusters");
	histo->GetYaxis()->SetTitle("Entries");
	histo = dynamic_cast<TH1F *>(fContainer->UncheckedAt(kNtrackletsHist));
	histo->GetXaxis()->SetTitle("Number of Tracklets");
	histo->GetYaxis()->SetTitle("Entries");
	histo = dynamic_cast<TH1F *>(fContainer->UncheckedAt(kNclusterTrackletHist));
	histo->GetXaxis()->SetTitle("Number of Clusters");
	histo->GetYaxis()->SetTitle("Entries");
	histo = dynamic_cast<TH1F *>(fContainer->UncheckedAt(kChi2));
	histo->GetXaxis()->SetTitle("#chi^2");
	histo->GetYaxis()->SetTitle("Entries");
	histo = dynamic_cast<TH1F *>(fContainer->UncheckedAt(kNTracksSectorHist));
	histo->GetXaxis()->SetTitle("Sector");
	histo->GetYaxis()->SetTitle("Number of Tracks");
	histo = dynamic_cast<TProfile *>(fContainer->UncheckedAt(kPulseHeight));
	histo->GetXaxis()->SetTitle("Time / 100ns");
	histo->GetYaxis()->SetTitle("Average Pulse Height (a. u.)");
	histo = dynamic_cast<TH1F *>(fContainer->UncheckedAt(kClusterCharge));
	histo->GetXaxis()->SetTitle("Cluster Charge (a.u.)");
	histo->GetYaxis()->SetTitle("Entries");
	histo = dynamic_cast<TH1F *>(fContainer->UncheckedAt(kChargeDeposit));
	histo->GetXaxis()->SetTitle("Charge Deposit (a.u.)");
	histo->GetYaxis()->SetTitle("Entries");
	// Calculate the Percentage of events containing tracks as function of the trigger cluster
	histo = dynamic_cast<TH1F *>(fContainer->UncheckedAt(kNEventsTrigger));
	Double_t nEvents = histo->Integral();
	TGraph *percentages = new TGraph(histo->GetNbinsX());
	for(Int_t ibin = 0; ibin < histo->GetNbinsX(); ibin++)
		percentages->SetPoint(ibin, ibin, histo->GetBinContent(histo->FindBin(ibin))/nEvents);
	percentages->GetXaxis()->SetTitle("Trigger Cluster");
	percentages->GetYaxis()->SetTitle("#%Events");
	percentages->SetMarkerColor(kBlue);
	percentages->SetMarkerStyle(22);
	fContainer->Add(percentages);
	fNRefFigures = 10;
	return kTRUE;
}

//_______________________________________________________
void AliTRDcheckDetector::GetRefFigure(Int_t ifig, Int_t &first, Int_t &last, Option_t *opt){
	//
	// Setting Reference Figures
	//
	switch(ifig){
		case 0:	first = last = kNTracksEventHist;
						break;
		case 1:	first = last = kNclustersHist;
						break;
		case 2:	first = last = kNtrackletsHist;
						break;
		case 3:	first = last = kNclusterTrackletHist;
						break;
		case 4:	first = last = kChi2;
						break;
		case 5:	first = last = kNTracksSectorHist;
						break;
		case 6:	first = last = kPulseHeight;
						break;
		case 7:	first = last = kClusterCharge;
						break;
		case 8:	first = last = kChargeDeposit;
						break;
		case 9: first = last = kPostProcessing;
						break;
		default: first = last = kNTracksEventHist;
						break;
	};
}

