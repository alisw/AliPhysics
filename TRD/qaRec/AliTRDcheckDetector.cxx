#include <TAxis.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TGaxis.h>
#include <TGraph.h>
#include <TMath.h>
#include <TMap.h>
#include <TObjArray.h>
#include <TObject.h>
#include <TObjString.h>
#include <TPad.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TROOT.h>

#include "AliLog.h"
#include "AliTRDcluster.h"
#include "AliESDHeader.h"
#include "AliESDRun.h"
#include "AliESDtrack.h"
#include "AliTRDgeometry.h"
#include "AliTRDpadPlane.h"
#include "AliTRDSimParam.h"
#include "AliTRDseedV1.h"
#include "AliTRDtrackV1.h"
#include "AliTRDtrackerV1.h"
#include "AliTRDReconstructor.h"
#include "AliTrackReference.h"
#include "AliTrackPointArray.h"
#include "AliTracker.h"
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
  ,fReconstructor(0x0)
  ,fGeo(0x0)
{
  //
  // Default constructor
  //
  DefineInput(1,AliTRDeventInfo::Class());
  fReconstructor = new AliTRDReconstructor;
  fReconstructor->SetRecoParam(AliTRDrecoParam::GetLowFluxParam());
  fGeo = new AliTRDgeometry;
  InitFunctorList();
}

//_______________________________________________________
AliTRDcheckDetector::~AliTRDcheckDetector(){
  //
  // Destructor
  // 
  if(fEventInfo) delete fEventInfo;
  if(fTriggerNames) delete fTriggerNames;
  delete fReconstructor;
  delete fGeo;
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
  histo = dynamic_cast<TH1F *>(fContainer->UncheckedAt(kNTrackletsVsFindable));
  histo->GetXaxis()->SetTitle("Ratio Found/Findable Tracklets");
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
  fNRefFigures = 11;
  return kTRUE;
}

//_______________________________________________________
void AliTRDcheckDetector::GetRefFigure(Int_t ifig){
  //
  // Setting Reference Figures
  //
  TH1 *h = 0x0, *h1 = 0x0, *h2 = 0x0;
  TGaxis *axis = 0x0;
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
    h = (TH1F*)fContainer->At(kNTrackletsVsFindable);
    if(!h->GetEntries()) break;
    h->Scale(100./h->Integral());
    h->GetXaxis()->SetRangeUser(0.005, 1.005);
    h->SetFillColor(kGreen);
    h->SetBarOffset(.2);
    h->SetBarWidth(.6);
    h->Draw("bar1");
    break;
  case 4:
    ((TH1F*)fContainer->At(kNclusterTrackletHist))->Draw("pc");
    break;
  case 5:
    ((TH1F*)fContainer->At(kChi2))->Draw("");
    break;
  case 6:
    h = (TH1F*)fContainer->At(kNTracksSectorHist);
    if(!h->GetEntries()) break;
    h->Scale(100./h->Integral());
    h->SetFillColor(kGreen);
    h->SetBarOffset(.2);
    h->SetBarWidth(.6);
    h->Draw("bar1");
    break;
  case 7:
    h = (TH1F*)fContainer->At(kPulseHeight);
    h->SetMarkerStyle(24);
    h->SetMarkerColor(kBlack);
    h->SetLineColor(kBlack);
    h->Draw("e1");
    // copy the second histogram in a new one with the same x-dimension as the phs with respect to time
    h1 = (TH1F *)fContainer->At(kPulseHeightDistance);
    h2 = new TH1F("hphs1","Average PH", 31, -0.5, 30.5);
    for(Int_t ibin = h1->GetXaxis()->GetFirst(); ibin < h1->GetNbinsX(); ibin++) 
      h2->SetBinContent(ibin, h1->GetBinContent(ibin));
    h2->SetMarkerStyle(22);
    h2->SetMarkerColor(kBlue);
    h2->SetLineColor(kBlue);
    h2->Draw("e1same");
    gPad->Update();
    // create axis according to the histogram dimensions of the original second histogram
    axis = new TGaxis(gPad->GetUxmin(),
                      gPad->GetUymax(),
                      gPad->GetUxmax(),
                      gPad->GetUymax(),
                      -0.08, 4.88, 510,"-L");
    axis->SetLineColor(kBlue);
    axis->SetLabelColor(kBlue);
    axis->SetTextColor(kBlue);
    axis->SetTitle("x_{c}-x_{0} / cm");
    axis->Draw();
    break;
  case 8:
    ((TH1F*)fContainer->At(kClusterCharge))->Draw("c");
    break;
  case 9:
    ((TH1F*)fContainer->At(kChargeDeposit))->Draw("c");
    break;
  case 10: 
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
  fContainer->AddAt(new TH1F("hNtlsFindable", "Ratio of found/findable Tracklets" , 101, -0.005, 1.005), kNTrackletsVsFindable);
  fContainer->AddAt(new TH1F("hNclTls","Mean Number of clusters per tracklet", 31, -0.5, 30.5), kNclusterTrackletHist);
  fContainer->AddAt(new TH1F("hChi2", "Chi2", 200, 0, 20), kChi2);
  fContainer->AddAt(new TH1F("hChi2n", "Norm. Chi2 (tracklets)", 50, 0, 5), kChi2Normalized);
  fContainer->AddAt(new TH1F("hSM", "Track Counts in Supermodule", 18, -0.5, 17.5), kNTracksSectorHist);
  // Detector signal on Detector-by-Detector basis
  fContainer->AddAt(new TProfile("hPHdetector", "Average PH", 31, -0.5, 30.5), kPulseHeight);
  fContainer->AddAt(new TProfile("hPHdistance", "Average PH", 31, -0.08, 4.88), kPulseHeightDistance);
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
      for(Int_t itime = 0; itime < AliTRDtrackerV1::GetNTimeBins(); itime++){
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
  Int_t nTracklets = GetNTracklets(fTrack);
  h->Fill(nTracklets);
  if(fDebugLevel > 3){
    if(nTracklets == 1){
      // If we have one Tracklet, check in which layer this happens
      Int_t layer = -1;
      AliTRDseedV1 *tracklet = 0x0;
      for(Int_t il = 0; il < AliTRDgeometry::kNlayer; il++){
        if((tracklet = fTrack->GetTracklet(il)) && tracklet->IsOK()){layer =  il; break;}
      }
      (*fDebugStream) << "PlotNTracklets"
        << "Layer=" << layer
        << "\n";
    }
  }
  return h;
}

//_______________________________________________________
TH1 *AliTRDcheckDetector::PlotTrackletsVsFindable(const AliTRDtrackV1 *track){
  //
  // Plots the ratio of number of tracklets vs.
  // number of findable tracklets
  //
  // Findable tracklets are defined as track prolongation
  // to layer i does not hit the dead area +- epsilon
  //
  // In order to check whether tracklet hist active area in Layer i, 
  // the track is refitted and the fitted position + an uncertainty 
  // range is compared to the chamber border (also with a different
  // uncertainty)
  //
  // For the track fit two cases are distinguished:
  // If the track is a stand alone track (defined by the status bit 
  // encoding, then the track is fitted with the tilted Rieman model
  // Otherwise the track is fitted with the Kalman fitter in two steps:
  // Since the track parameters are give at the outer point, we first 
  // fit in direction inwards. Afterwards we fit again in direction outwards
  // to extrapolate the track to layers which are not reached by the track
  // For the Kalman model, the radial track points have to be shifted by
  // a distance epsilon in the direction that we want to fit
  //
  const Float_t epsilon = 0.01;   // dead area tolerance
  const Float_t epsilon_R = 1;    // shift in radial direction of the anode wire position (Kalman filter only)
  const Float_t delta_y = 0.7;    // Tolerance in the track position in y-direction
  const Float_t delta_z = 7.0;    // Tolerance in the track position in z-direction (Padlength)
  Double_t x_anode[AliTRDgeometry::kNlayer] = {300.2, 312.8, 325.4, 338.0, 350.6, 363.2}; // Take the default X0
 
  if(track) fTrack = track;
  if(!fTrack){
    AliWarning("No Track defined.");
    return 0x0;
  }
  TH1 *h = 0x0;
  if(!(h = dynamic_cast<TH1F *>(fContainer->At(kNTrackletsVsFindable)))){
    AliWarning("No Histogram defined.");
    return 0x0;
  }
  Int_t nFound = 0, nFindable = 0;
  Int_t stack = -1;
  Double_t ymin = 0., ymax = 0., zmin = 0., zmax = 0.;
  Double_t y = 0., z = 0.;
  AliTRDseedV1 *tracklet = 0x0;
  AliTRDpadPlane *pp;  
  for(Int_t il = 0; il < AliTRDgeometry::kNlayer; il++){
    if((tracklet = fTrack->GetTracklet(il)) && tracklet->IsOK()){
      tracklet->SetReconstructor(fReconstructor);
      nFound++;
    }
  }
  // 2 Different cases:
  // 1st stand alone: here we cannot propagate, but be can do a Tilted Rieman Fit
  // 2nd barrel track: here we propagate the track to the layers
  AliTrackPoint points[6];
  Float_t xyz[3];
  memset(xyz, 0, sizeof(Float_t) * 3);
  if(((fESD->GetStatus() & AliESDtrack::kTRDout) > 0) && !((fESD->GetStatus() & AliESDtrack::kTRDin) > 0)){
    // stand alone track
    for(Int_t il = 0; il < AliTRDgeometry::kNlayer; il++){
      xyz[0] = x_anode[il];
      points[il].SetXYZ(xyz);
    }
    AliTRDtrackerV1::FitRiemanTilt(const_cast<AliTRDtrackV1 *>(fTrack), 0x0, kTRUE, 6, points);
  } else {
    // barrel track
    //
    // 2 Steps:
    // -> Kalman inwards
    // -> Kalman outwards
    AliTRDtrackV1 copy_track(*fTrack);  // Do Kalman on a (non-constant) copy of the track
    AliTrackPoint points_inward[6], points_outward[6];
    for(Int_t il = AliTRDgeometry::kNlayer; il--;){
      // In order to avoid complications in the Kalman filter if the track points have the same radial
      // position like the tracklets, we have to shift the radial postion of the anode wire by epsilon
      // in the direction we want to go
      // The track points have to be in reverse order for the Kalman Filter inwards
      xyz[0] = x_anode[AliTRDgeometry::kNlayer - il - 1] - epsilon_R;
      points_inward[il].SetXYZ(xyz);
      xyz[0] = x_anode[il] + epsilon_R;
      points_outward[il].SetXYZ(xyz);
    }
    /*for(Int_t ipt = 0; ipt < AliTRDgeometry::kNlayer; ipt++)
      printf("%d. X = %f\n", ipt, points[ipt].GetX());*/
    // Kalman inwards
    AliTRDtrackerV1::FitKalman(&copy_track, 0x0, kFALSE, 6, points_inward);
    memcpy(points, points_inward, sizeof(AliTrackPoint) * 6); // Preliminary store the inward results in the Array points
    // Kalman outwards
    AliTRDtrackerV1::FitKalman(&copy_track, 0x0, kTRUE, 6, points_inward);
    memcpy(points, points_outward, sizeof(AliTrackPoint) * AliTRDgeometry::kNlayer);
  }
  for(Int_t il = 0; il < AliTRDgeometry::kNlayer; il++){
    y = points[il].GetY();
    z = points[il].GetZ();
    if((stack = fGeo->GetStack(z, il)) < 0) continue; // Not findable
    pp = fGeo->GetPadPlane(il, stack);
    ymin = pp->GetCol0() + epsilon;
    ymax = pp->GetColEnd() - epsilon; 
    zmin = pp->GetRowEnd() + epsilon; 
    zmax = pp->GetRow0() - epsilon;
    // ignore y-crossing (material)
    if((z + delta_z > zmin && z - delta_z < zmax) && (y + delta_y > ymin && y - delta_y < ymax)) nFindable++;
      if(fDebugLevel > 3){
        Double_t pos_tracklet[2] = {tracklet ? tracklet->GetYfit(0) : 0, tracklet ? tracklet->GetMeanz() : 0};
        Int_t hasTracklet = tracklet ? 1 : 0;
        (*fDebugStream)   << "GetFindableTracklets"
          << "layer="     << il
          << "ytracklet=" << pos_tracklet[0]
          << "ytrack="    << y
          << "ztracklet=" << pos_tracklet[1]
          << "ztrack="    << z
          << "tracklet="  << hasTracklet
          << "\n";
      }
  }
  
  h->Fill(nFindable > 0 ? TMath::Min(nFound/static_cast<Double_t>(nFindable), 1.) : 1);
  if(fDebugLevel > 2) AliInfo(Form("Findable[Found]: %d[%d|%f]", nFindable, nFound, nFound/static_cast<Float_t>(nFindable > 0 ? nFindable : 1)));
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
    for(Int_t itime = 0; itime < AliTRDtrackerV1::GetNTimeBins(); itime++){
      if(!(c = tracklet->GetClusters(itime))) continue;
      Int_t localtime        = c->GetLocalTimeBin();
      Double_t absolute_charge = TMath::Abs(c->GetQ());
      h->Fill(localtime, absolute_charge);
      if(fDebugLevel > 3){
        Int_t crossing = tracklet->GetNChange();
        Int_t detector = c->GetDetector();
        Double_t distance[2];
        GetDistanceToTracklet(distance, tracklet, c);
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
          << "dy="        << distance[0]
          << "dz="        << distance[1]
          << "c.="        << c
          << "\n";
      }
    }
  }
  return h;
}

//_______________________________________________________
TH1 *AliTRDcheckDetector::PlotPHSdistance(const AliTRDtrackV1 *track){
  //
  // Plots the average pulse height vs the distance from the anode wire
  // (plus const anode wire offset)
  //
  if(track) fTrack = track;
  if(!fTrack){
    AliWarning("No Track defined.");
    return 0x0;
  }
  TProfile *h = 0x0;
  if(!(h = dynamic_cast<TProfile *>(fContainer->At(kPulseHeightDistance)))){
    AliWarning("No Histogram defined.");
    return 0x0;
  }
  Int_t offset = AliTRDSimParam::Instance()->GetAnodeWireOffset();
  AliTRDseedV1 *tracklet = 0x0;
  AliTRDcluster *c = 0x0;
  Double_t distance = 0;
  for(Int_t itl = 0; itl < AliTRDgeometry::kNlayer; itl++){
    if(!(tracklet = fTrack->GetTracklet(itl)) || !(tracklet->IsOK())) continue;
    tracklet->ResetClusterIter();
    while((c = tracklet->NextCluster())){
      distance = tracklet->GetX0() - c->GetX() + offset;
      h->Fill(distance, TMath::Abs(c->GetQ()));
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
    for(Int_t itime = 0; itime < AliTRDtrackerV1::GetNTimeBins(); itime++){
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
  Int_t nTracklets = 0;
  if(fDebugLevel > 3)
    nTracklets = GetNTracklets(fTrack); // fill NTracklet to the Debug Stream
  for(Int_t itl = 0x0; itl < AliTRDgeometry::kNlayer; itl++){
    if(!(tracklet = fTrack->GetTracklet(itl)) || !tracklet->IsOK()) continue;
    Qtot = 0;
    c1 = 0x0;
    for(Int_t itime = 0; itime < AliTRDtrackerV1::GetNTimeBins(); itime++){
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
        << "nTracklets="<< nTracklets
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
    for(Int_t itime = 0; itime < AliTRDtrackerV1::GetNTimeBins(); itime++){
      if(!(c = tracklet->GetClusters(itime))) continue;
      sector = static_cast<Int_t>(c->GetDetector()/AliTRDgeometry::kNdets);
    }
    break;
  }
  h->Fill(sector);
  return h;
}

//_______________________________________________________
Int_t AliTRDcheckDetector::GetNTracklets(const AliTRDtrackV1 *track){
  //
  // Count the number of tracklets per track
  //
  if(!track){
    AliError("No track");
    return 0;
  }
  Int_t nTracklets = 0;
  AliTRDseedV1 *tracklet = 0x0;
  for(Int_t il = AliTRDgeometry::kNlayer; il--;){
    if((tracklet = track->GetTracklet(il)) && tracklet->IsOK()) nTracklets++;
  }
  return nTracklets;
}

//________________________________________________________
void AliTRDcheckDetector::SetRecoParam(AliTRDrecoParam *r)
{

  fReconstructor->SetRecoParam(r);
}

//________________________________________________________
void AliTRDcheckDetector::GetDistanceToTracklet(Double_t *dist, AliTRDseedV1 *tracklet, AliTRDcluster *c){
  dist[0] = c->GetY() - (tracklet->GetYfit(0) + tracklet->GetYfit(1)*(c->GetX() - tracklet->GetX0()));
  dist[1] = c->GetZ() - tracklet->GetMeanz();
}
