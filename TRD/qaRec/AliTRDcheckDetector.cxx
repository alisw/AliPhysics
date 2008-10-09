#include <TFile.h>
#include <TH1F.h>
#include <TGraph.h>
#include <TMath.h>
#include <TMap.h>
#include <TObjArray.h>
#include <TObject.h>
#include <TObjString.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TROOT.h>

#include "AliTRDcluster.h"
#include "AliESDHeader.h"
#include "AliESDRun.h"
#include "AliTRDseedV1.h"
#include "AliTRDtrackV1.h"
#include "TTreeStream.h"

#include "AliTRDtrackInfo/AliTRDtrackInfo.h"
#include "AliTRDtrackInfo/AliTRDeventInfo.h"
#include "AliTRDcheckDetector.h"

#include <cstdio>
#include <iostream>

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
  AliTRDrecoTask("DetChecker", "Basic Detector Checker")
  ,fEventInfo(0x0)
  ,fTriggerNames(0x0)
{
  //
  // Default constructor
  //
  DefineInput(1,AliTRDeventInfo::Class());
}

//_______________________________________________________
AliTRDcheckDetector::~AliTRDcheckDetector(){
  //
  // Destructor
  // 
  if(fEventInfo) delete fEventInfo;
  if(fTriggerNames) delete fTriggerNames;
}

//_______________________________________________________
void AliTRDcheckDetector::ConnectInputData(Option_t *opt){
	//
	// Connect the Input data with the task
	//
	AliTRDrecoTask::ConnectInputData(opt);
	fEventInfo = dynamic_cast<AliTRDeventInfo *>(GetInputData(1));
}

//_______________________________________________________
void AliTRDcheckDetector::CreateOutputObjects(){
  //
  // Create Output Objects
  //
  fContainer = new TObjArray(25);
  // Register Histograms
  fContainer->Add(new TH1F("hNtrks", "Number of Tracks per event", 100, 0, 100));
  fContainer->Add(new TH1F("hEventsTriggerTracks", "Trigger Class (Tracks)", 100, 0, 100));
  fContainer->Add(new TH1F("hNcls", "Nr. of clusters per track", 181, -0.5, 180.5));
  fContainer->Add(new TH1F("hNtls", "Nr. tracklets per track", 7, -0.5, 6.5));
  fContainer->Add(new TH1F("hNclTls","Mean Number of clusters per tracklet", 31, -0.5, 30.5));
  fContainer->Add(new TH1F("hChi2", "Chi2", 200, 0, 20));
  fContainer->Add(new TH1F("hChi2n", "Norm. Chi2 (tracklets)", 50, 0, 5));
  fContainer->Add(new TH1F("hSM", "Track Counts in Supermodule", 18, -0.5, 17.5));
	// Detector signal on Detector-by-Detector basis
  fContainer->Add(new TProfile("hPHdetector", "Average PH", 31, -0.5, 30.5));
  fContainer->Add(new TH1F("hQclDetector", "Cluster charge", 200, 0, 1200));
  fContainer->Add(new TH1F("hQTdetector", "Total Charge Deposit", 6000, 0, 6000));
  fContainer->Add(new TH1F("hEventsTrigger", "Trigger Class", 100, 0, 100));
  
  fTriggerNames = new TMap();
}

//_______________________________________________________
void AliTRDcheckDetector::Exec(Option_t *){
  //
  // Execution function
  // Filling TRD quality histos
  //
  if(!HasMCdata() && fEventInfo->GetEventHeader()->GetEventType() != 7) return;	// For real data we select only physical events
  Int_t nTracks = 0;		// Count the number of tracks per event
  Int_t triggermask = fEventInfo->GetEventHeader()->GetTriggerMask();
  TString triggername =  fEventInfo->GetRunInfo()->GetFiredTriggerClasses(triggermask);
  if(fDebugLevel > 6)printf("Trigger cluster: %d, Trigger class: %s\n", triggermask, triggername.Data());
  dynamic_cast<TH1F *>(fContainer->UncheckedAt(kNEventsTrigger))->Fill(triggermask);
  AliTRDtrackInfo *fTrackInfo = 0x0;
  AliTRDtrackV1 *fTRDtrack = 0x0;
  AliTRDseedV1 *fTracklet = 0x0;
  AliTRDcluster *fTRDcluster = 0x0;
  for(Int_t iti = 0; iti < fTracks->GetEntriesFast(); iti++){
    fTrackInfo = dynamic_cast<AliTRDtrackInfo *>(fTracks->UncheckedAt(iti));
    if(!fTrackInfo || !(fTRDtrack = fTrackInfo->GetTRDtrack())) continue;
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
            dynamic_cast<TProfile *>(fContainer->UncheckedAt(kPulseHeight))->Fill(localtime, absolute_charge);
            dynamic_cast<TH1F *>(fContainer->UncheckedAt(kClusterCharge))->Fill(absolute_charge);
            if(fDebugLevel > 2){
            	(*fDebugStream) << "PulseHeight"
            		<< "Detector="	<< detector
            		<< "Sector="		<< sector
            		<< "Timebin="		<< localtime
            		<< "Charge="		<< absolute_charge
            		<< "\n";
            }
          }
          dynamic_cast<TH1F *>(fContainer->UncheckedAt(kChargeDeposit))->Fill(Qtot);
          if(fDebugLevel > 3){
          	(*fDebugStream) << "ChargeDeposit"
          		<< "Detector="	<< detector
          		<< "Sector="		<< sector
          		<< "QT="				<< Qtot
          		<< "\n";
          }
        }
      }
      dynamic_cast<TH1F *>(fContainer->UncheckedAt(kNTracksSectorHist))->Fill(sector);
    }
    nTracks++;
  }
  if(nTracks){
  	dynamic_cast<TH1F *>(fContainer->UncheckedAt(kNEventsTriggerTracks))->Fill(triggermask);
  	dynamic_cast<TH1F *>(fContainer->UncheckedAt(kNTracksEventHist))->Fill(nTracks);
  }
  if(triggermask <= 20 && !fTriggerNames->FindObject(Form("%d", triggermask))){
    fTriggerNames->Add(new TObjString(Form("%d", triggermask)), new TObjString(triggername));
    // also set the label for both histograms
    TH1 *histo = dynamic_cast<TH1F *>(fContainer->UncheckedAt(kNEventsTriggerTracks));
    histo->GetXaxis()->SetBinLabel(histo->FindBin(triggermask), triggername);
    histo = dynamic_cast<TH1F *>(fContainer->UncheckedAt(kNEventsTrigger));
    histo->GetXaxis()->SetBinLabel(histo->FindBin(triggermask), triggername);
  }
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
	
	TH1 * histo = 0x0;
	histo = dynamic_cast<TH1F *>(fContainer->UncheckedAt(kNTracksEventHist));
	histo->GetXaxis()->SetTitle("Number of Tracks");
	histo->GetYaxis()->SetTitle("Events");
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
	
	// Calculate the purity of the trigger clusters
	histo = dynamic_cast<TH1F *>(fContainer->UncheckedAt(kNEventsTrigger));
	TH1F *histoTracks = dynamic_cast<TH1F *>(fContainer->UncheckedAt(kNEventsTriggerTracks));
	histoTracks->Divide(histo);
	Float_t purities[20], val = 0;
	TString triggernames[20];
	Int_t nTriggerClasses = 0;
	for(Int_t ibin = 1; ibin <= histo->GetNbinsX(); ibin++){
		if((val = histoTracks->GetBinContent(ibin))){
			purities[nTriggerClasses] = val;
			triggernames[nTriggerClasses] = histoTracks->GetXaxis()->GetBinLabel(ibin);
			nTriggerClasses++;
		}
	}
	TH1F *hTriggerInf = new TH1F("fTriggerInf", "Trigger Information", TMath::Max(nTriggerClasses, 1), 0, TMath::Max(nTriggerClasses, 1));
	for(Int_t ibin = 1; ibin <= nTriggerClasses; ibin++){
		hTriggerInf->SetBinContent(ibin, purities[ibin-1]);
		hTriggerInf->GetXaxis()->SetBinLabel(ibin, triggernames[ibin-1].Data());
	}
	hTriggerInf->GetXaxis()->SetTitle("Trigger Cluster");
	hTriggerInf->GetYaxis()->SetTitle("Ratio");
	hTriggerInf->GetYaxis()->SetRangeUser(0,1);
//	hTriggerInf->SetMarkerColor(kBlue);
//	hTriggerInf->SetMarkerStyle(22);
	fContainer->Add(hTriggerInf);
	fNRefFigures = 10;
	return kTRUE;
}

//_______________________________________________________
void AliTRDcheckDetector::GetRefFigure(Int_t ifig, Int_t &first, Int_t &last, Option_t *opt){
	//
	// Setting Reference Figures
	//
	opt = "pl";
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
		case 9: first = last = kPurity;
						opt="bar";
						break;
		default: first = last = kNTracksEventHist;
						break;
	};
}

