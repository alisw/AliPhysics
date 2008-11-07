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

#include "AliLog.h"
#include "AliTRDcluster.h"
#include "AliESDHeader.h"
#include "AliESDRun.h"
#include "AliTRDgeometry.h"
#include "AliTRDseedV1.h"
#include "AliTRDtrackV1.h"
#include "AliTrackReference.h"
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
  InitFunctorList();
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
  OpenFile(0,"RECREATE");
  fContainer = Histos();
  if(!fTriggerNames) fTriggerNames = new TMap();
}

//_______________________________________________________
void AliTRDcheckDetector::Exec(Option_t *opt){
  //
  // Execution function
  // Filling TRD quality histos
  //
  if(!HasMCdata() && fEventInfo->GetEventHeader()->GetEventType() != 7) return;	// For real data we select only physical events
  AliTRDrecoTask::Exec(opt);  
  Int_t nTracks = 0;		// Count the number of tracks per event
  Int_t triggermask = fEventInfo->GetEventHeader()->GetTriggerMask();
  TString triggername =  fEventInfo->GetRunInfo()->GetFiredTriggerClasses(triggermask);
  if(fDebugLevel > 6)printf("Trigger cluster: %d, Trigger class: %s\n", triggermask, triggername.Data());
  dynamic_cast<TH1F *>(fContainer->UncheckedAt(kNEventsTrigger))->Fill(triggermask);
  for(Int_t iti = 0; iti < fTracks->GetEntriesFast(); iti++){
    if(!fTracks->UncheckedAt(iti)) continue;
    AliTRDtrackInfo *fTrackInfo = dynamic_cast<AliTRDtrackInfo *>(fTracks->UncheckedAt(iti));
    if(!fTrackInfo->GetTrack()) continue;
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
void AliTRDcheckDetector::GetRefFigure(Int_t ifig){
  //
  // Setting Reference Figures
  //
  TH1 *h = 0x0;
  switch(ifig){
  case 0:	
    ((TH1F*)fContainer->At(kNTracksEventHist))->Draw("pl");
    break;
  case 1:
    ((TH1F*)fContainer->At(kNclustersHist))->Draw("pl");
    break;
  case 2:
    h = (TH1F*)fContainer->At(kNtrackletsHist);
    if(!h->GetEntries()) break;
    h->Scale(100./h->Integral());
    h->GetXaxis()->SetRangeUser(.5, 6.5);
    h->SetFillColor(kGreen);
    h->SetBarOffset(.2);
    h->SetBarWidth(.6);
    h->Draw("bar1");
    break;
  case 3:
    ((TH1F*)fContainer->At(kNclusterTrackletHist))->Draw("pc");
    break;
  case 4:
    ((TH1F*)fContainer->At(kChi2))->Draw("");
    break;
  case 5:
    h = (TH1F*)fContainer->At(kNTracksSectorHist);
    if(!h->GetEntries()) break;
    h->Scale(100./h->Integral());
    h->SetFillColor(kGreen);
    h->SetBarOffset(.2);
    h->SetBarWidth(.6);
    h->Draw("bar1");
    break;
  case 6:
    h = (TH1F*)fContainer->At(kPulseHeight);
    h->SetMarkerStyle(24);
    h->Draw("e1");
    break;
  case 7:
    ((TH1F*)fContainer->At(kClusterCharge))->Draw("c");
    break;
  case 8:
    ((TH1F*)fContainer->At(kChargeDeposit))->Draw("c");
    break;
  case 9: 
    h=(TH1F*)fContainer->At(kPurity);
    h->SetBarOffset(.2);
    h->SetBarWidth(.6);
    h->SetFillColor(kGreen);
    h->Draw("bar1");
    break;
  default:
    ((TH1F*)fContainer->At(kNTracksEventHist))->Draw("pl");
    break;
  }
}

//_______________________________________________________
TObjArray *AliTRDcheckDetector::Histos(){
  //
  // Create QA histograms
  //
  if(fContainer) return fContainer;
  
  fContainer = new TObjArray(25);
  // Register Histograms
  fContainer->AddAt(new TH1F("hNtrks", "Number of Tracks per event", 100, 0, 100), kNTracksEventHist);
  fContainer->AddAt(new TH1F("hEventsTriggerTracks", "Trigger Class (Tracks)", 100, 0, 100), kNEventsTriggerTracks);
  fContainer->AddAt(new TH1F("hNcls", "Nr. of clusters per track", 181, -0.5, 180.5), kNclustersHist);
  fContainer->AddAt(new TH1F("hNtls", "Nr. tracklets per track", 7, -0.5, 6.5), kNtrackletsHist);
  fContainer->AddAt(new TH1F("hNclTls","Mean Number of clusters per tracklet", 31, -0.5, 30.5), kNclusterTrackletHist);
  fContainer->AddAt(new TH1F("hChi2", "Chi2", 200, 0, 20), kChi2);
  fContainer->AddAt(new TH1F("hChi2n", "Norm. Chi2 (tracklets)", 50, 0, 5), kChi2Normalized);
  fContainer->AddAt(new TH1F("hSM", "Track Counts in Supermodule", 18, -0.5, 17.5), kNTracksSectorHist);
  // Detector signal on Detector-by-Detector basis
  fContainer->AddAt(new TProfile("hPHdetector", "Average PH", 31, -0.5, 30.5), kPulseHeight);
  fContainer->AddAt(new TH1F("hQclDetector", "Cluster charge", 200, 0, 1200), kClusterCharge);
  fContainer->AddAt(new TH1F("hQTdetector", "Total Charge Deposit", 6000, 0, 6000), kChargeDeposit);
  fContainer->AddAt(new TH1F("hEventsTrigger", "Trigger Class", 100, 0, 100), kNEventsTrigger);

  return fContainer;
}

/*
* Plotting Functions
*/

//_______________________________________________________
TH1 *AliTRDcheckDetector::PlotMeanNClusters(const AliTRDtrackV1 *track){
  //
  // Plot the mean number of clusters per tracklet
  //
  if(track) fTrack = track;
  if(!fTrack){
    AliWarning("No Track defined.");
    return 0x0;
  }
  TH1 *h = 0x0;
  if(!(h = dynamic_cast<TH1F *>(fContainer->At(kNclusterTrackletHist)))){
    AliWarning("No Histogram defined.");
    return 0x0;
  }
  AliTRDseedV1 *tracklet = 0x0;
  for(Int_t itl = 0; itl < AliTRDgeometry::kNlayer; itl++){
    if(!(tracklet = fTrack->GetTracklet(itl)) || !tracklet->IsOK()) continue;
    h->Fill(tracklet->GetN());
  }
  return h;
}

//_______________________________________________________
TH1 *AliTRDcheckDetector::PlotNClusters(const AliTRDtrackV1 *track){
  //
  // Plot the number of clusters in one track
  //
  if(track) fTrack = track;
  if(!fTrack){
    AliWarning("No Track defined.");
    return 0x0;
  }
  TH1 *h = 0x0;
  if(!(h = dynamic_cast<TH1F *>(fContainer->At(kNclustersHist)))){
    AliWarning("No Histogram defined.");
    return 0x0;
  }
  
  Int_t nclusters = 0;
  AliTRDseedV1 *tracklet = 0x0;
  for(Int_t itl = 0; itl < AliTRDgeometry::kNlayer; itl++){
    if(!(tracklet = fTrack->GetTracklet(itl)) || !tracklet->IsOK()) continue;
    nclusters += tracklet->GetN();
    if(fDebugLevel > 2){
      Int_t crossing = tracklet->GetNChange();
      AliTRDcluster *c = 0x0;
      for(Int_t itime = 0; itime < kNTimeBins; itime++){
        if(!(c = tracklet->GetClusters(itime))) continue;
        break;
      }
      Int_t detector = c->GetDetector();
      Float_t sector = static_cast<Int_t>(detector/AliTRDgeometry::kNdets);
      Float_t theta = TMath::ATan(tracklet->GetZfit(1));
      Float_t phi = TMath::ATan(tracklet->GetYfit(1));
      Float_t momentum = 0.;
      Int_t pdg = 0;
      Int_t kinkIndex = fESD ? fESD->GetKinkIndex() : 0;
      UShort_t TPCncls = fESD ? fESD->GetTPCncls() : 0;
      if(fMC){
        if(fMC->GetTrackRef()) momentum = fMC->GetTrackRef()->P();
        pdg = fMC->GetPDG();
      }
      (*fDebugStream) << "NClusters"
        << "Detector="  << detector
        << "Sector="    << sector
        << "crossing="  << crossing
        << "momentum="	<< momentum
        << "pdg="				<< pdg
        << "theta="			<< theta
        << "phi="				<< phi
        << "kinkIndex="	<< kinkIndex
        << "TPCncls="		<< TPCncls
        << "nclusters=" << nclusters
        << "\n";
    }
  }
  h->Fill(nclusters);
  return h;
}

//_______________________________________________________
TH1 *AliTRDcheckDetector::PlotChi2(const AliTRDtrackV1 *track){
  //
  // Plot the chi2 of the track
  //
  if(track) fTrack = track;
  if(!fTrack){
    AliWarning("No Track defined.");
    return 0x0;
  }
  TH1 *h = 0x0;
  if(!(h = dynamic_cast<TH1F *>(fContainer->At(kChi2)))){
    AliWarning("No Histogram defined.");
    return 0x0;
  }
  h->Fill(fTrack->GetChi2());
  return h;
}

//_______________________________________________________
TH1 *AliTRDcheckDetector::PlotNormalizedChi2(const AliTRDtrackV1 *track){
  //
  // Plot the chi2 of the track
  //
  if(track) fTrack = track;
  if(!fTrack){
    AliWarning("No Track defined.");
    return 0x0;
  }
  TH1 *h = 0x0;
  if(!(h = dynamic_cast<TH1F *>(fContainer->At(kChi2Normalized)))){
    AliWarning("No Histogram defined.");
    return 0x0;
  }
  Int_t nTracklets = 0;
  AliTRDseedV1 *tracklet = 0x0;
  for(Int_t itl = 0; itl < AliTRDgeometry::kNlayer; itl++){
    if(!(tracklet = fTrack->GetTracklet(itl)) || !tracklet->IsOK()) continue;
    nTracklets++;
  }
  h->Fill(fTrack->GetChi2()/nTracklets);
  return h;
}


//_______________________________________________________
TH1 *AliTRDcheckDetector::PlotNTracklets(const AliTRDtrackV1 *track){
  //
  // Plot the number of tracklets
  //
  if(track) fTrack = track;
  if(!fTrack){
    AliWarning("No Track defined.");
    return 0x0;
  }
  TH1 *h = 0x0;
  if(!(h = dynamic_cast<TH1F *>(fContainer->At(kNtrackletsHist)))){
    AliWarning("No Histogram defined.");
    return 0x0;
  }
  Int_t nTracklets = 0;
  AliTRDseedV1 *tracklet = 0x0;
  for(Int_t itl = 0; itl < AliTRDgeometry::kNlayer; itl++){
    if(!(tracklet = fTrack->GetTracklet(itl)) || !tracklet->IsOK()) continue;
    nTracklets++;
  }
  h->Fill(nTracklets);
  return h;
}

//_______________________________________________________
TH1 *AliTRDcheckDetector::PlotPulseHeight(const AliTRDtrackV1 *track){
  //
  // Plot the average pulse height
  //
  if(track) fTrack = track;
  if(!fTrack){
    AliWarning("No Track defined.");
    return 0x0;
  }
  TProfile *h = 0x0;
  if(!(h = dynamic_cast<TProfile *>(fContainer->At(kPulseHeight)))){
    AliWarning("No Histogram defined.");
    return 0x0;
  }
  AliTRDseedV1 *tracklet = 0x0;
  AliTRDcluster *c = 0x0;
  for(Int_t itl = 0; itl < AliTRDgeometry::kNlayer; itl++){
    if(!(tracklet = fTrack->GetTracklet(itl)) || !tracklet->IsOK())continue;
    for(Int_t itime = 0; itime < kNTimeBins; itime++){
      if(!(c = tracklet->GetClusters(itime))) continue;
      Int_t localtime        = c->GetLocalTimeBin();
      Double_t absolute_charge = TMath::Abs(c->GetQ());
      h->Fill(localtime, absolute_charge);
      if(fDebugLevel > 3){
        Int_t crossing = tracklet->GetNChange();
        Int_t detector = c->GetDetector();
        Float_t sector = static_cast<Int_t>(detector/AliTRDgeometry::kNdets);
        Float_t theta = TMath::ATan(tracklet->GetZfit(1));
        Float_t phi = TMath::ATan(tracklet->GetYfit(1));
        Float_t momentum = 0.;
        Int_t pdg = 0;
        Int_t kinkIndex = fESD ? fESD->GetKinkIndex() : 0;
        UShort_t TPCncls = fESD ? fESD->GetTPCncls() : 0;
        if(fMC){
          if(fMC->GetTrackRef()) momentum = fMC->GetTrackRef()->P();
          pdg = fMC->GetPDG();
        }
        (*fDebugStream) << "PulseHeight"
          << "Detector="	<< detector
          << "Sector="		<< sector
          << "crossing="	<< crossing
          << "Timebin="		<< localtime
          << "Charge="		<< absolute_charge
          << "momentum="	<< momentum
          << "pdg="				<< pdg
          << "theta="			<< theta
          << "phi="				<< phi
          << "kinkIndex="	<< kinkIndex
          << "TPCncls="		<< TPCncls
          << "\n";
      }
    }
  }
  return h;
}

//_______________________________________________________
TH1 *AliTRDcheckDetector::PlotClusterCharge(const AliTRDtrackV1 *track){
  //
  // Plot the cluster charge
  //
  if(track) fTrack = track;
  if(!fTrack){
    AliWarning("No Track defined.");
    return 0x0;
  }
  TH1 *h = 0x0;
  if(!(h = dynamic_cast<TH1F *>(fContainer->At(kClusterCharge)))){
    AliWarning("No Histogram defined.");
    return 0x0;
  }
  AliTRDseedV1 *tracklet = 0x0;
  AliTRDcluster *c = 0x0;
  for(Int_t itl = 0; itl < AliTRDgeometry::kNlayer; itl++){
    if(!(tracklet = fTrack->GetTracklet(itl)) || !tracklet->IsOK())continue;
    for(Int_t itime = 0; itime < kNTimeBins; itime++){
      if(!(c = tracklet->GetClusters(itime))) continue;
      h->Fill(c->GetQ());
    }
  }
  return h;
}

//_______________________________________________________
TH1 *AliTRDcheckDetector::PlotChargeDeposit(const AliTRDtrackV1 *track){
  //
  // Plot the charge deposit per chamber
  //
  if(track) fTrack = track;
  if(!fTrack){
    AliWarning("No Track defined.");
    return 0x0;
  }
  TH1 *h = 0x0;
  if(!(h = dynamic_cast<TH1F *>(fContainer->At(kChargeDeposit)))){
    AliWarning("No Histogram defined.");
    return 0x0;
  }
  AliTRDseedV1 *tracklet = 0x0;
  AliTRDcluster *c = 0x0, *c1 = 0x0;	// c1 for the Debug Stream
  Double_t Qtot = 0;
  for(Int_t itl = 0x0; itl < AliTRDgeometry::kNlayer; itl++){
    if(!(tracklet = fTrack->GetTracklet(itl)) || !tracklet->IsOK()) continue;
    Qtot = 0;
    c1 = 0x0;
    for(Int_t itime = 0; itime < kNTimeBins; itime++){
      if(!(c = tracklet->GetClusters(itime))) continue;
      if(!c1) c1 = c;
      Qtot += TMath::Abs(c->GetQ());
    }
    h->Fill(Qtot);
    if(fDebugLevel > 3){
      Int_t crossing = tracklet->GetNChange();
      Int_t detector = c1->GetDetector();
      Float_t sector = static_cast<Int_t>(detector/AliTRDgeometry::kNdets);
      Float_t theta = TMath::ATan(tracklet->GetZfit(1));
      Float_t phi = TMath::ATan(tracklet->GetYfit(1));
      Float_t momentum = 0.;
      Int_t pdg = 0;
      Int_t kinkIndex = fESD ? fESD->GetKinkIndex() : 0;
      UShort_t TPCncls = fESD ? fESD->GetTPCncls() : 0;
      if(fMC){
	      if(fMC->GetTrackRef()) momentum = fMC->GetTrackRef()->P();
        pdg = fMC->GetPDG();
      }
      (*fDebugStream) << "ChargeDeposit"
        << "Detector="  << detector
        << "Sector="    << sector
        << "crossing="  << crossing
        << "momentum="	<< momentum
        << "pdg="				<< pdg
        << "theta="			<< theta
        << "phi="				<< phi
        << "kinkIndex="	<< kinkIndex
        << "TPCncls="		<< TPCncls
        << "QT="        << Qtot
        << "\n";
    }
  }
  return h;
}

//_______________________________________________________
TH1 *AliTRDcheckDetector::PlotTracksSector(const AliTRDtrackV1 *track){
  //
  // Plot the number of tracks per Sector
  //
  if(track) fTrack = track;
  if(!fTrack){
    AliWarning("No Track defined.");
    return 0x0;
  }
  TH1 *h = 0x0;
  if(!(h = dynamic_cast<TH1F *>(fContainer->At(kNTracksSectorHist)))){
    AliWarning("No Histogram defined.");
    return 0x0;
  }
  AliTRDseedV1 *tracklet = 0x0;
  AliTRDcluster *c = 0x0;
  Int_t sector = -1;
  for(Int_t itl = 0; itl < AliTRDgeometry::kNlayer; itl++){
    if(!(tracklet = fTrack->GetTracklet(itl)) || !tracklet->IsOK()) continue;
    for(Int_t itime = 0; itime < kNTimeBins; itime++){
      if(!(c = tracklet->GetClusters(itime))) continue;
      sector = static_cast<Int_t>(c->GetDetector()/AliTRDgeometry::kNdets);
    }
    break;
  }
  h->Fill(sector);
  return h;
}

