////////////////////////////////////////////////////////////////////////////
//                                                                        //
//                                                                        //
//  Basic checks for tracking and detector performance                    //
//  
//     For the moment (15.10.2009) the following checks are implemented    //
//       - Number of clusters/track
//       - Number of clusters/tracklet
//       - Number of tracklets/track from different sources (Barrel, StandAlone)
//       - Number of findable tracklets
//       - Number of tracks per event and TRD sector
//       - <PH>
//       - Chi2 distribution for tracks
//       - Charge distribution per cluster and tracklet
//       - Number of events and tracks per trigger source 
//       - Trigger purity
//       - Track and Tracklet propagation status
//
//  Authors:                                                              //
//    Anton Andronic <A.Andronic@gsi.de>                                  //
//    Alexandru Bercuci <A.Bercuci@gsi.de>                                //
//    Markus Fasel <M.Fasel@gsi.de>                                       //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include <TAxis.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH1I.h>
#include <TH2F.h>
#include <TF1.h>
#include <TGaxis.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TLinearFitter.h>
#include <TMath.h>
#include <TMap.h>
#include <TProfile2D.h>
#include <TObjArray.h>
#include <TObject.h>
#include <TObjString.h>

#include <TPad.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TROOT.h>
#include <TChain.h>

#include "AliLog.h"
#include "AliTRDcluster.h"
#include "AliESDHeader.h"
#include "AliESDRun.h"
#include "AliESDtrack.h"
#include "AliExternalTrackParam.h"
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

#include "info/AliTRDtrackInfo.h"
#include "info/AliTRDeventInfo.h"
#include "AliTRDcheckDET.h"

#include <cstdio>
#include <iostream>

ClassImp(AliTRDcheckDET)

//_______________________________________________________
AliTRDcheckDET::AliTRDcheckDET():
  AliTRDrecoTask()
  ,fEventInfo(NULL)
  ,fTriggerNames(NULL)
  ,fReconstructor(NULL)
  ,fGeo(NULL)
  ,fFlags(0)
{
  //
  // Default constructor
  //
  SetNameTitle("checkDET", "Basic TRD data checker");
}

//_______________________________________________________
AliTRDcheckDET::AliTRDcheckDET(char* name):
  AliTRDrecoTask(name, "Basic TRD data checker")
  ,fEventInfo(NULL)
  ,fTriggerNames(NULL)
  ,fReconstructor(NULL)
  ,fGeo(NULL)
  ,fFlags(0)
{
  //
  // Default constructor
  //
  DefineInput(2, AliTRDeventInfo::Class());

  fReconstructor = new AliTRDReconstructor;
  fReconstructor->SetRecoParam(AliTRDrecoParam::GetLowFluxParam());
  fGeo = new AliTRDgeometry;
  InitFunctorList();
}


//_______________________________________________________
AliTRDcheckDET::~AliTRDcheckDET(){
  //
  // Destructor
  // 
  if(fTriggerNames) delete fTriggerNames;
  delete fReconstructor;
  delete fGeo;
}

//_______________________________________________________
void AliTRDcheckDET::UserCreateOutputObjects(){
  //
  // Create Output Objects
  //
  if(!HasFunctorList()) InitFunctorList();
  fContainer = Histos();
  if(!fTriggerNames) fTriggerNames = new TMap();
}

//_______________________________________________________
void AliTRDcheckDET::UserExec(Option_t *opt){
  //
  // Execution function
  // Filling TRD quality histos
  //

  fEventInfo = dynamic_cast<AliTRDeventInfo *>(GetInputData(2));
  AliDebug(2, Form("EventInfo[%p] Header[%p]", (void*)fEventInfo, (void*)(fEventInfo?fEventInfo->GetEventHeader():NULL)));

  AliTRDrecoTask::UserExec(opt);  
  Int_t nTracks = 0;		// Count the number of tracks per event
  Int_t triggermask = fEventInfo->GetEventHeader()->GetTriggerMask();
  TString triggername =  fEventInfo->GetRunInfo()->GetFiredTriggerClasses(triggermask);
  AliDebug(6, Form("Trigger cluster: %d, Trigger class: %s\n", triggermask, triggername.Data()));
  dynamic_cast<TH1F *>(fContainer->UncheckedAt(kNeventsTrigger))->Fill(triggermask);
  for(Int_t iti = 0; iti < fTracks->GetEntriesFast(); iti++){
    if(!fTracks->UncheckedAt(iti)) continue;
    AliTRDtrackInfo *fTrackInfo = dynamic_cast<AliTRDtrackInfo *>(fTracks->UncheckedAt(iti));
    if(!fTrackInfo->GetTrack()) continue;
    nTracks++;
  }

  if(nTracks){
    dynamic_cast<TH1F *>(fContainer->UncheckedAt(kNeventsTriggerTracks))->Fill(triggermask);
    dynamic_cast<TH1F *>(fContainer->UncheckedAt(kNtracksEvent))->Fill(nTracks);
  }
  if(triggermask <= 20 && !fTriggerNames->FindObject(Form("%d", triggermask))){
    fTriggerNames->Add(new TObjString(Form("%d", triggermask)), new TObjString(triggername));
    // also set the label for both histograms
    TH1 *histo = dynamic_cast<TH1F *>(fContainer->UncheckedAt(kNeventsTriggerTracks));
    histo->GetXaxis()->SetBinLabel(histo->FindBin(triggermask), triggername);
    histo = dynamic_cast<TH1F *>(fContainer->UncheckedAt(kNeventsTrigger));
    histo->GetXaxis()->SetBinLabel(histo->FindBin(triggermask), triggername);
  }
  PostData(1, fContainer);
}


//_______________________________________________________
Bool_t AliTRDcheckDET::PostProcess(){
  //
  // Do Postprocessing (for the moment set the number of Reference histograms)
  //
  
  TH1 * h = NULL;
  
  // Calculate of the trigger clusters purity
  h = dynamic_cast<TH1F *>(fContainer->FindObject("hEventsTrigger"));
  TH1F *h1 = dynamic_cast<TH1F *>(fContainer->FindObject("hEventsTriggerTracks"));
  h1->Divide(h);
  Float_t purities[20], val = 0;
  TString triggernames[20];
  Int_t nTriggerClasses = 0;
  for(Int_t ibin = 1; ibin <= h->GetNbinsX(); ibin++){
    if((val = h1->GetBinContent(ibin))){
      purities[nTriggerClasses] = val;
      triggernames[nTriggerClasses] = h1->GetXaxis()->GetBinLabel(ibin);
      nTriggerClasses++;
    }
  }
  h = dynamic_cast<TH1F *>(fContainer->UncheckedAt(kTriggerPurity));
  TAxis *ax = h->GetXaxis();
  for(Int_t itrg = 0; itrg < nTriggerClasses; itrg++){
    h->Fill(itrg, purities[itrg]);
    ax->SetBinLabel(itrg+1, triggernames[itrg].Data());
  }
  ax->SetRangeUser(-0.5, nTriggerClasses+.5);
  h->GetYaxis()->SetRangeUser(0,1);

  // track status
  h=dynamic_cast<TH1F*>(fContainer->At(kTrackStatus));
  Float_t ok = h->GetBinContent(1);
  Int_t nerr = h->GetNbinsX();
  for(Int_t ierr=nerr; ierr--;){
    h->SetBinContent(ierr+1, ok>0.?1.e2*h->GetBinContent(ierr+1)/ok:0.);
  }
  h->SetBinContent(1, 0.);

  // tracklet status
  
  TObjArray *arr = dynamic_cast<TObjArray*>(fContainer->UncheckedAt(kTrackletStatus));
  for(Int_t ily = AliTRDgeometry::kNlayer; ily--;){
    h=dynamic_cast<TH1F*>(arr->At(ily));
    Float_t okB = h->GetBinContent(1);
    Int_t nerrB = h->GetNbinsX();
    for(Int_t ierr=nerrB; ierr--;){
      h->SetBinContent(ierr+1, okB>0.?1.e2*h->GetBinContent(ierr+1)/okB:0.);
    }
    h->SetBinContent(1, 0.);
  }

  fNRefFigures = 17;

  return kTRUE;
}

//_______________________________________________________
void AliTRDcheckDET::MakeSummary(){
  //
  // Create summary plots for TRD check DET
  // This function creates 2 summary plots:
  // - General Quantities
  // - PHS
  // The function will reuse GetRefFigure
  //
  
  TCanvas *cOut = new TCanvas(Form("summary%s1", GetName()), Form("Summary 1 for task %s", GetName()), 1024, 768);
  cOut->Divide(3,3);

  // Create figures using GetRefFigure
  cOut->cd(1); GetRefFigure(kFigNtracksEvent);  
  cOut->cd(2); GetRefFigure(kFigNtracksSector);
  cOut->cd(3); GetRefFigure(kFigNclustersTrack);
  cOut->cd(4); GetRefFigure(kFigNclustersTracklet);
  cOut->cd(5); GetRefFigure(kFigNtrackletsTrack);
  cOut->cd(6); GetRefFigure(kFigNTrackletsP);
  cOut->cd(7); GetRefFigure(kFigChargeCluster);
  cOut->cd(8); GetRefFigure(kFigChargeTracklet);
  cOut->SaveAs(Form("TRDsummary%s1.gif", GetName()));
  delete cOut;

  // Second Plot: PHS
  cOut = new TCanvas(Form("summary%s2", GetName()), Form("Summary 2 for task %s", GetName()), 1024, 512);
  cOut->cd(); GetRefFigure(kFigPH);
  cOut->SaveAs(Form("TRDsummary%s2.gif", GetName())); 
  delete cOut;

  // Third Plot: Mean Number of clusters as function of eta, phi and layer
   cOut = new TCanvas(Form("summary%s3", GetName()), Form("Summary 3 for task %s", GetName()), 1024, 768);
  cOut->cd(); MakePlotMeanClustersLayer();
  cOut->SaveAs(Form("TRDsummary%s3.gif", GetName())); 
  delete cOut;

}

//_______________________________________________________
Bool_t AliTRDcheckDET::GetRefFigure(Int_t ifig){
  //
  // Setting Reference Figures
  //
  gPad->SetLogy(0);
  gPad->SetLogx(0);
  TH1 *h = NULL; TObjArray *arr=NULL;
  TLegend *leg = NULL;
  Bool_t kFIRST(1);
  switch(ifig){
  case kFigNclustersTrack:
    (h=(TH1F*)fContainer->FindObject("hNcls"))->Draw("pl");
    PutTrendValue("NClustersTrack", h->GetMean());
    PutTrendValue("NClustersTrackRMS", h->GetRMS());
    return kTRUE;
  case kFigNclustersTracklet:
    (h =(TH1F*)fContainer->FindObject("hNclTls"))->Draw("pc");
    PutTrendValue("NClustersTracklet", h->GetMean());
    PutTrendValue("NClustersTrackletRMS", h->GetRMS());
    return kTRUE;
  case kFigNtrackletsTrack:
    h=MakePlotNTracklets();
    PutTrendValue("NTrackletsTrack", h->GetMean());
    PutTrendValue("NTrackletsTrackRMS", h->GetRMS());
    return kTRUE;
  case kFigNTrackletsP:
    MakePlotnTrackletsVsP();
    return kTRUE;
  case kFigNtrackletsCross:
    h = (TH1F*)fContainer->FindObject("hNtlsCross");
    if(!MakeBarPlot(h, kRed)) break;
    PutTrendValue("NTrackletsCross", h->GetMean());
    PutTrendValue("NTrackletsCrossRMS", h->GetRMS());
    return kTRUE;
  case kFigNtrackletsFindable:
    h = (TH1F*)fContainer->FindObject("hNtlsFindable");
    if(!MakeBarPlot(h, kGreen)) break;
    PutTrendValue("NTrackletsFindable", h->GetMean());
    PutTrendValue("NTrackletsFindableRMS", h->GetRMS());
    return kTRUE;
  case kFigNtracksEvent:
    (h = (TH1F*)fContainer->FindObject("hNtrks"))->Draw("pl");
    PutTrendValue("NTracksEvent", h->GetMean());
    PutTrendValue("NTracksEventRMS", h->GetRMS());
    return kTRUE;
  case kFigNtracksSector:
    h = (TH1F*)fContainer->FindObject("hNtrksSector");
    if(!MakeBarPlot(h, kGreen)) break;
    PutTrendValue("NTracksSector", h->Integral()/h->GetNbinsX());
    return kTRUE;
  case kFigTrackStatus:
    if(!(h=(TH1F *)fContainer->FindObject("hTrackStatus"))) break;
    h->GetXaxis()->SetRangeUser(0.5, -1);
    h->GetYaxis()->CenterTitle();
    h->Draw("c");
    PutTrendValue("TrackStatus", h->Integral());
    gPad->SetLogy(0);
    return kTRUE;
  case kFigTrackletStatus:
    if(!(arr = dynamic_cast<TObjArray*>(fContainer->At(kTrackletStatus)))) break;
    leg = new TLegend(.68, .7, .97, .97);
    leg->SetBorderSize(0);leg->SetFillStyle(0);
    leg->SetHeader("TRD layer");
    for(Int_t ily=AliTRDgeometry::kNlayer; ily--;){
      if(!(h=dynamic_cast<TH1F*>(arr->At(ily)))) continue;
      if(kFIRST){
        h->Draw("pl");
        h->GetXaxis()->SetRangeUser(0.5, -1);
        h->GetYaxis()->CenterTitle();
        kFIRST = kFALSE;
      } else h->Draw("samepl");
      leg->AddEntry(h, Form("ly = %d", ily), "l");
      PutTrendValue(Form("TrackletStatus%d", ily), h->Integral());
    }
    leg->Draw();
    gPad->SetLogy(0);
    return kTRUE;
  case kFigChi2:
    return kTRUE;
    MakePlotChi2();
    return kTRUE;
  case kFigPH:
    gPad->SetMargin(0.125, 0.015, 0.1, 0.1);
    MakePlotPulseHeight();
    gPad->SetLogy(0);
    return kTRUE;
  case kFigChargeCluster:
    (h = (TH1F*)fContainer->FindObject("hQcl"))->Draw("c");
    gPad->SetLogy(1);
    PutTrendValue("ChargeCluster", h->GetMaximumBin());
    PutTrendValue("ChargeClusterRMS", h->GetRMS());
    return kTRUE;
  case kFigChargeTracklet:
    (h=(TH1F*)fContainer->FindObject("hQtrklt"))->Draw("c");
    PutTrendValue("ChargeTracklet", h->GetMaximumBin());
    PutTrendValue("ChargeTrackletRMS", h->GetRMS());
    return kTRUE;
  case kFigNeventsTrigger:
    ((TH1F*)fContainer->FindObject("hEventsTrigger"))->Draw("");
    return kTRUE;
  case kFigNeventsTriggerTracks:
    ((TH1F*)fContainer->FindObject("hEventsTriggerTracks"))->Draw("");
    return kTRUE;
  case kFigTriggerPurity: 
    if(!MakeBarPlot((TH1F*)fContainer->FindObject("hTriggerPurity"), kGreen)) break;
    break;
  default:
    break;
  }
  AliInfo(Form("Reference plot [%d] missing result", ifig));
  return kFALSE;
}

//_______________________________________________________
TObjArray *AliTRDcheckDET::Histos(){
  //
  // Create QA histograms
  //
    
  if(fContainer) return fContainer;
  
  fContainer = new TObjArray(20);
  //fContainer->SetOwner(kTRUE);

  // Register Histograms
  TH1 * h = NULL;
  TAxis *ax = NULL;
  if(!(h = (TH1F *)gROOT->FindObject("hNcls"))){
    h = new TH1F("hNcls", "N_{clusters} / track", 181, -0.5, 180.5);
    h->GetXaxis()->SetTitle("N_{clusters}");
    h->GetYaxis()->SetTitle("Entries");
  } else h->Reset();
  fContainer->AddAt(h, kNclustersTrack);

  TObjArray *arr = new TObjArray(AliTRDgeometry::kNlayer);
  arr->SetOwner(kTRUE);  arr->SetName("clusters");
  fContainer->AddAt(arr, kNclustersLayer);
  for(Int_t ily=AliTRDgeometry::kNlayer; ily--;){
    if(!(h = (TProfile2D *)gROOT->FindObject(Form("hNcl%d", ily)))){
      h = new TProfile2D(Form("hNcl%d", ily), Form("Mean Number of clusters in Layer %d", ily), 100, -1.0, 1.0, 50, -1.1*TMath::Pi(), 1.1*TMath::Pi());
      h->GetXaxis()->SetTitle("#eta");
      h->GetYaxis()->SetTitle("#phi");
    } else h->Reset();
    arr->AddAt(h, ily);
  }

  if(!(h = (TH1F *)gROOT->FindObject("hNclTls"))){
    h = new TH1F("hNclTls","N_{clusters} / tracklet", 51, -0.5, 50.5);
    h->GetXaxis()->SetTitle("N_{clusters}");
    h->GetYaxis()->SetTitle("Entries");
  } else h->Reset();
  fContainer->AddAt(h, kNclustersTracklet);

  if(!(h = (TH1F *)gROOT->FindObject("hNtls"))){
    h = new TH1F("hNtls", "N_{tracklets} / track", AliTRDgeometry::kNlayer, 0.5, 6.5);
    h->GetXaxis()->SetTitle("N^{tracklet}");
    h->GetYaxis()->SetTitle("freq. [%]");
  } else h->Reset();
  fContainer->AddAt(h, kNtrackletsTrack);

  if(!(h = (TH1F *)gROOT->FindObject("htlsSTA"))){
    h = new TH1F("hNtlsSTA", "N_{tracklets} / track (Stand Alone)", AliTRDgeometry::kNlayer, 0.5, 6.5);
    h->GetXaxis()->SetTitle("N^{tracklet}");
    h->GetYaxis()->SetTitle("freq. [%]");
  }
  fContainer->AddAt(h, kNtrackletsSTA);

  // Binning for momentum dependent tracklet Plots
  const Int_t kNp(30);
  Float_t P=0.2;
  Float_t binsP[kNp+1], binsTrklt[AliTRDgeometry::kNlayer+1];
  for(Int_t i=0;i<kNp+1; i++,P+=(TMath::Exp(i*i*.001)-1.)) binsP[i]=P;
  for(Int_t il = 0; il <= AliTRDgeometry::kNlayer; il++) binsTrklt[il] = 0.5 + il;
  if(!(h = (TH1F *)gROOT->FindObject("htlsBAR"))){
    // Make tracklets for barrel tracks momentum dependent (if we do not exceed min and max values)
    h = new TH2F("hNtlsBAR", 
    "N_{tracklets} / track;p [GeV/c];N^{tracklet};freq. [%]", 
    kNp, binsP, AliTRDgeometry::kNlayer, binsTrklt);
  }
  fContainer->AddAt(h, kNtrackletsBAR);

  // 
  if(!(h = (TH1F *)gROOT->FindObject("hNtlsCross"))){
    h = new TH1F("hNtlsCross", "N_{tracklets}^{cross} / track", 7, -0.5, 6.5);
    h->GetXaxis()->SetTitle("n_{row cross}");
    h->GetYaxis()->SetTitle("freq. [%]");
  } else h->Reset();
  fContainer->AddAt(h, kNtrackletsCross);

  if(!(h = (TH1F *)gROOT->FindObject("hNtlsFindable"))){
    h = new TH1F("hNtlsFindable", "Found/Findable Tracklets" , 101, -0.005, 1.005);
    h->GetXaxis()->SetTitle("r [a.u]");
    h->GetYaxis()->SetTitle("Entries");
  } else h->Reset();
  fContainer->AddAt(h, kNtrackletsFindable);

  if(!(h = (TH1F *)gROOT->FindObject("hNtrks"))){
    h = new TH1F("hNtrks", "N_{tracks} / event", 100, 0, 100);
    h->GetXaxis()->SetTitle("N_{tracks}");
    h->GetYaxis()->SetTitle("Entries");
  } else h->Reset();
  fContainer->AddAt(h, kNtracksEvent);

  if(!(h = (TH1F *)gROOT->FindObject("hNtrksSector"))){
    h = new TH1F("hNtrksSector", "N_{tracks} / sector", AliTRDgeometry::kNsector, -0.5, 17.5);
    h->GetXaxis()->SetTitle("sector");
    h->GetYaxis()->SetTitle("freq. [%]");
  } else h->Reset();
  fContainer->AddAt(h, kNtracksSector);

  if(!(h = (TH1F*)gROOT->FindObject("hTrackStatus"))){
    const Int_t nerr = 7;
    h = new TH1F("hTrackStatus", "Track Status", nerr, -0.5, nerr-0.5);
    const Char_t *label[nerr] = {"OK", "PROL", "PROP", "AJST", "SNP", "TINI", "UPDT"};
    ax = h->GetXaxis();
    for(Int_t ierr=nerr; ierr--;) ax->SetBinLabel(ierr+1, label[ierr]);
    h->SetYTitle("Relative Error to Good [%]");
  }
  fContainer->AddAt(h, kTrackStatus);

  arr = new TObjArray(AliTRDgeometry::kNlayer);
  arr->SetOwner(kTRUE);  arr->SetName("TrackletStatus");
  fContainer->AddAt(arr, kTrackletStatus);
  for(Int_t ily=AliTRDgeometry::kNlayer; ily--;){
    if(!(h = (TH1F *)gROOT->FindObject(Form("hTrackletStatus%d", ily)))){
      const Int_t nerr = 8;
      h = new TH1F(Form("hTrackletStatus%d", ily), "Tracklet status", nerr, -0.5, nerr-0.5);
      h->SetLineColor(ily+1);
      const Char_t *label[nerr] = {"OK", "Geom", "Bound", "NoCl", "NoAttach", "NoClTr", "NoFit", "Chi2"};
      ax = h->GetXaxis();
      for(Int_t ierr=nerr; ierr--;) ax->SetBinLabel(ierr+1, label[ierr]);
      h->SetYTitle("Relative Error to Good [%]");
    } else h->Reset();
    arr->AddAt(h, ily);
  }

  // <PH> histos
  arr = new TObjArray(3);
  arr->SetOwner(kTRUE);  arr->SetName("<PH>");
  fContainer->AddAt(arr, kPH);
  if(!(h = (TH1F *)gROOT->FindObject("hPHt"))){
    h = new TProfile("hPHt", "<PH>", 31, -0.5, 30.5);
    h->GetXaxis()->SetTitle("Time / 100ns");
    h->GetYaxis()->SetTitle("<PH> [a.u]");
  } else h->Reset();
  arr->AddAt(h, 0);
  if(!(h = (TH1F *)gROOT->FindObject("hPHx")))
    h = new TProfile("hPHx", "<PH>", 31, -0.08, 4.88);
  else h->Reset();
  arr->AddAt(h, 1);
  if(!(h = (TH2F *)gROOT->FindObject("hPH2D"))){
    h = new TH2F("hPH2D", "Charge Distribution / time", 31, -0.5, 30.5, 100, 0, 1024);
    h->GetXaxis()->SetTitle("Time / 100ns");
    h->GetYaxis()->SetTitle("Charge / a.u.");
  } else h->Reset();
  arr->AddAt(h, 2);

  // Chi2 histos
  if(!(h = (TH2S*)gROOT->FindObject("hChi2"))){
    h = new TH2S("hChi2", "#chi^{2} per track", AliTRDgeometry::kNlayer, .5, AliTRDgeometry::kNlayer+.5, 100, 0, 50);
    h->SetXTitle("ndf");
    h->SetYTitle("#chi^{2}/ndf");
    h->SetZTitle("entries");
  } else h->Reset();
  fContainer->AddAt(h, kChi2);

  if(!(h = (TH1F *)gROOT->FindObject("hQcl"))){
    h = new TH1F("hQcl", "Q_{cluster}", 200, 0, 1200);
    h->GetXaxis()->SetTitle("Q_{cluster} [a.u.]");
    h->GetYaxis()->SetTitle("Entries");
  }else h->Reset();
  fContainer->AddAt(h, kChargeCluster);

  if(!(h = (TH1F *)gROOT->FindObject("hQtrklt"))){
    h = new TH1F("hQtrklt", "Q_{tracklet}", 6000, 0, 6000);
    h->GetXaxis()->SetTitle("Q_{tracklet} [a.u.]");
    h->GetYaxis()->SetTitle("Entries");
  }else h->Reset();
  fContainer->AddAt(h, kChargeTracklet);


  if(!(h = (TH1F *)gROOT->FindObject("hEventsTrigger")))
    h = new TH1F("hEventsTrigger", "Trigger Class", 100, 0, 100);
  else h->Reset();
  printf("Histos Adding \n");
  
  fContainer->AddAt(h, kNeventsTrigger);

  if(!(h = (TH1F *)gROOT->FindObject("hEventsTriggerTracks")))
    h = new TH1F("hEventsTriggerTracks", "Trigger Class (Tracks)", 100, 0, 100);
  else h->Reset();
  fContainer->AddAt(h, kNeventsTriggerTracks);

  if(!(h = (TH1F *)gROOT->FindObject("hTriggerPurity"))){
    h = new TH1F("hTriggerPurity", "Trigger Purity", 10, -0.5, 9.5);
    h->GetXaxis()->SetTitle("Trigger Cluster");
    h->GetYaxis()->SetTitle("freq.");
  } else h->Reset();
  fContainer->AddAt(h, kTriggerPurity);

  return fContainer;
}

/*
* Plotting Functions
*/

//_______________________________________________________
TH1 *AliTRDcheckDET::PlotTrackStatus(const AliTRDtrackV1 *track)
{
//
// Plot the track propagation status. The following errors are defined (see AliTRDtrackV1::ETRDtrackError)
//   PROL - track prolongation failure
//   PROP - track propagation failure
//   AJST - crossing sectors failure
//   SNP  - too large bending
//   TINI - tracklet initialization failure
//   UPDT - track position/covariance update failure 
//
// Performance plot looks as below:
//Begin_Html
//<img src="TRD/trackStatus.gif">
//End_Html
//
  if(track) fkTrack = track;
  if(!fkTrack){
    AliDebug(4, "No Track defined.");
    return NULL;
  }
  TH1 *h = NULL;
  if(!(h = dynamic_cast<TH1F *>(fContainer->At(kTrackStatus)))){
    AliWarning("No Histogram defined.");
    return NULL;
  }
  h->Fill(fkTrack->GetStatusTRD());
  return h;
}

//_______________________________________________________
TH1 *AliTRDcheckDET::PlotTrackletStatus(const AliTRDtrackV1 *track)
{
//
// Plot the tracklet propagation status. The following errors are defined for tracklet (see AliTRDtrackV1::ETRDlayerError)
//   Geom   - 
//   Bound  - tracklet too close to chamber walls
//   NoCl   - no clusters in the track roads
//   NoAttach - fail to attach clusters
//   NoClTr - fail to use clusters for fit
//   NoFit  - tracklet fit failled
//   Chi2   - chi2 tracklet-track over threshold
//
// Performance plot looks as below:
//Begin_Html
//<img src="TRD/trackletStatus.gif">
//End_Html
//
  if(track) fkTrack = track;
  if(!fkTrack){
    AliDebug(4, "No Track defined.");
    return NULL;
  }
  TObjArray *arr =NULL;
  if(!(arr = dynamic_cast<TObjArray*>(fContainer->At(kTrackletStatus)))){
    AliWarning("Histograms not defined.");
    return NULL;
  }

  TH1 *h = NULL;
  for(Int_t ily=AliTRDgeometry::kNlayer; ily--;){
    if(!(h = dynamic_cast<TH1F*>(arr->At(ily)))){
      AliWarning(Form("Missing histo for layer %d.", ily));
      continue;
    }
    h->Fill(fkTrack->GetStatusTRD(ily));
  }
  return h;
}

//_______________________________________________________
TH1 *AliTRDcheckDET::PlotNClustersTracklet(const AliTRDtrackV1 *track){
  //
  // Plot the mean number of clusters per tracklet
  //
  if(track) fkTrack = track;
  if(!fkTrack){
    AliDebug(4, "No Track defined.");
    return NULL;
  }
  AliExternalTrackParam *par = fkTrack->GetTrackIn() ? fkTrack->GetTrackIn() : fkTrack->GetTrackOut();
  TH1 *h = NULL;
  TProfile2D *hlayer = NULL;
  Double_t eta = 0., phi = 0.;
  if(!(h = dynamic_cast<TH1F *>(fContainer->At(kNclustersTracklet)))){
    AliWarning("No Histogram defined.");
    return NULL;
  }
  AliTRDseedV1 *tracklet = NULL;
  TObjArray *histosLayer = dynamic_cast<TObjArray *>(fContainer->At(kNclustersLayer));
  if(!histosLayer){
    AliWarning("No Histograms for single layer defined");
  }
  for(Int_t itl = 0; itl < AliTRDgeometry::kNlayer; itl++){
    if(!(tracklet = fkTrack->GetTracklet(itl)) || !tracklet->IsOK()) continue;
    h->Fill(tracklet->GetN2());
    if(histosLayer && par){
      if((hlayer = dynamic_cast<TProfile2D *>(histosLayer->At(itl)))){
        GetEtaPhiAt(par, tracklet->GetX0(), eta, phi);
        hlayer->Fill(eta, phi, tracklet->GetN2());
      }
    }
  }
  return h;
}

//_______________________________________________________
TH1 *AliTRDcheckDET::PlotNClustersTrack(const AliTRDtrackV1 *track){
  //
  // Plot the number of clusters in one track
  //
  if(track) fkTrack = track;
  if(!fkTrack){
    AliDebug(4, "No Track defined.");
    return NULL;
  }
  TH1 *h = NULL;
  if(!(h = dynamic_cast<TH1F *>(fContainer->At(kNclustersTrack)))){
    AliWarning("No Histogram defined.");
    return NULL;
  }
  
  Int_t nclusters = 0;
  AliTRDseedV1 *tracklet = NULL;
  AliExternalTrackParam *par = fkTrack->GetTrackOut() ? fkTrack->GetTrackOut() : fkTrack->GetTrackIn();
  if(!par) return NULL;
  Double_t momentumRec = par->P();
  for(Int_t itl = 0; itl < AliTRDgeometry::kNlayer; itl++){
    if(!(tracklet = fkTrack->GetTracklet(itl)) || !tracklet->IsOK()) continue;
    Int_t n(tracklet->GetN());
    nclusters += n;
    if(DebugLevel() > 2){
      Int_t crossing = Int_t(tracklet->IsRowCross());
      Int_t detector = tracklet->GetDetector();
      Float_t theta = TMath::ATan(tracklet->GetZref(1));
      Float_t phi = TMath::ATan(tracklet->GetYref(1));
      Float_t momentumMC = 0.;
      Int_t pdg = 0;
      Int_t kinkIndex = fkESD ? fkESD->GetKinkIndex() : 0;
      UShort_t nclsTPC = fkESD ? fkESD->GetTPCncls() : 0;
      if(fkMC){
        if(fkMC->GetTrackRef()) momentumMC = fkMC->GetTrackRef()->P();
        pdg = fkMC->GetPDG();
      }
      (*DebugStream()) << "NClustersTrack"
        << "Detector="  << detector
        << "crossing="  << crossing
        << "momentumMC="<< momentumMC
        << "momentumRec=" << momentumRec
        << "pdg="				<< pdg
        << "theta="			<< theta
        << "phi="				<< phi
        << "kinkIndex="	<< kinkIndex
        << "TPCncls="		<< nclsTPC
        << "TRDncls="   << n
        << "\n";
    }
  }
  h->Fill(nclusters);
  return h;
}


//_______________________________________________________
TH1 *AliTRDcheckDET::PlotNTrackletsTrack(const AliTRDtrackV1 *track){
  //
  // Plot the number of tracklets
  //
  if(track) fkTrack = track;
  if(!fkTrack){
    AliDebug(4, "No Track defined.");
    return NULL;
  }
  TH1 *h = NULL, *hSta = NULL; TH2 *hBarrel = NULL;
  if(!(h = dynamic_cast<TH1F *>(fContainer->At(kNtrackletsTrack)))){
    AliWarning("No Histogram defined.");
    return NULL;
  }
  Int_t nTracklets = fkTrack->GetNumberOfTracklets();
  h->Fill(nTracklets);
  if(!fkESD) return h;
  Int_t status = fkESD->GetStatus();

/*  printf("in/out/refit/pid: TRD[%d|%d|%d|%d]\n", status &AliESDtrack::kTRDin ? 1 : 0, status &AliESDtrack::kTRDout ? 1 : 0, status &AliESDtrack::kTRDrefit ? 1 : 0, status &AliESDtrack::kTRDpid ? 1 : 0);*/
  Double_t p = 0.;
  Int_t method = -1;    // to distinguish between stand alone and full barrel tracks in the debugging
  if((status & AliESDtrack::kTRDin) != 0){
    method = 1;
    // Full Barrel Track: Save momentum dependence
    if(!(hBarrel = dynamic_cast<TH2F *>(fContainer->At(kNtrackletsBAR)))){
      AliWarning("Method: Barrel.  Histogram not processed!");
      return NULL;
    }
    AliExternalTrackParam *par(fkTrack->GetTrackIn());
    if(!par){
      AliError("Input track params missing");
      return NULL;
    }
    p = par->P(); // p needed later in the debug streaming
    hBarrel->Fill(p, nTracklets);
  } else {
    // Stand alone Track: momentum dependence not usefull
    method = 0;
    if(!(hSta = dynamic_cast<TH1F *>(fContainer->At(kNtrackletsSTA)))) {
      AliWarning("Method: StandAlone.  Histogram not processed!");
      return NULL;
    }
    hSta->Fill(nTracklets);
  }

  if(DebugLevel() > 2){
    AliTRDseedV1 *tracklet = NULL;
    Int_t sector = -1, stack = -1, detector;
    for(Int_t itl = 0; itl < AliTRDgeometry::kNlayer; itl++){
      if(!(tracklet = fkTrack->GetTracklet(itl)) || !(tracklet->IsOK())) continue;
      detector = tracklet->GetDetector();
      sector = fGeo->GetSector(detector);
      stack = fGeo->GetStack(detector);
      break;
    }
    (*DebugStream()) << "NTrackletsTrack"
      << "Sector="      << sector
      << "Stack="        << stack 
      << "NTracklets="  << nTracklets
      << "Method="      << method
      << "p="           << p
      << "\n";
  }
  if(DebugLevel() > 3){
    AliTRDseedV1 *tracklet = NULL;
    for(Int_t il = 0; il < AliTRDgeometry::kNlayer; il++){
      if((tracklet = fkTrack->GetTracklet(il)) && tracklet->IsOK()){
        (*DebugStream()) << "NTrackletsLayer"
        << "Layer=" << il
        << "p=" << p
        << "\n";
      }
    }
  }
  return h;
}


//_______________________________________________________
TH1 *AliTRDcheckDET::PlotNTrackletsRowCross(const AliTRDtrackV1 *track){
  //
  // Plot the number of tracklets
  //
  if(track) fkTrack = track;
  if(!fkTrack){
    AliDebug(4, "No Track defined.");
    return NULL;
  }
  TH1 *h = NULL;
  if(!(h = dynamic_cast<TH1F *>(fContainer->At(kNtrackletsCross)))){
    AliWarning("No Histogram defined.");
    return NULL;
  }

  Int_t ncross = 0;
  AliTRDseedV1 *tracklet = NULL;
  for(Int_t il = 0; il < AliTRDgeometry::kNlayer; il++){
    if(!(tracklet = fkTrack->GetTracklet(il)) || !tracklet->IsOK()) continue;

    if(tracklet->IsRowCross()) ncross++;
  }
  h->Fill(ncross);
  return h;
}

//_______________________________________________________
TH1 *AliTRDcheckDET::PlotFindableTracklets(const AliTRDtrackV1 *track){
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
  const Float_t epsilonR = 1;    // shift in radial direction of the anode wire position (Kalman filter only)
  const Float_t deltaY = 0.7;    // Tolerance in the track position in y-direction
  const Float_t deltaZ = 7.0;    // Tolerance in the track position in z-direction (Padlength)
  Double_t xAnode[AliTRDgeometry::kNlayer] = {300.2, 312.8, 325.4, 338.0, 350.6, 363.2}; // Take the default X0
 
  if(track) fkTrack = track;
  if(!fkTrack){
    AliDebug(4, "No Track defined.");
    return NULL;
  }
  TH1 *h = NULL;
  if(!(h = dynamic_cast<TH1F *>(fContainer->At(kNtrackletsFindable)))){
    AliWarning("No Histogram defined.");
    return NULL;
  }
  Int_t nFound = 0, nFindable = 0;
  Int_t stack = -1;
  Double_t ymin = 0., ymax = 0., zmin = 0., zmax = 0.;
  Double_t y = 0., z = 0.;
  AliTRDseedV1 *tracklet = NULL;
  AliTRDpadPlane *pp;  
  for(Int_t il = 0; il < AliTRDgeometry::kNlayer; il++){
    if((tracklet = fkTrack->GetTracklet(il)) && tracklet->IsOK()){
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
  if(((fkESD->GetStatus() & AliESDtrack::kTRDout) > 0) && !((fkESD->GetStatus() & AliESDtrack::kTRDin) > 0)){
    // stand alone track
    for(Int_t il = 0; il < AliTRDgeometry::kNlayer; il++){
      xyz[0] = xAnode[il];
      points[il].SetXYZ(xyz);
    }
    AliTRDtrackerV1::FitRiemanTilt(const_cast<AliTRDtrackV1 *>(fkTrack), NULL, kTRUE, 6, points);
  } else {
    // barrel track
    //
    // 2 Steps:
    // -> Kalman inwards
    // -> Kalman outwards
    AliTRDtrackV1 copyTrack(*fkTrack);  // Do Kalman on a (non-constant) copy of the track
    AliTrackPoint pointsInward[6], pointsOutward[6];
    for(Int_t il = AliTRDgeometry::kNlayer; il--;){
      // In order to avoid complications in the Kalman filter if the track points have the same radial
      // position like the tracklets, we have to shift the radial postion of the anode wire by epsilon
      // in the direction we want to go
      // The track points have to be in reverse order for the Kalman Filter inwards
      xyz[0] = xAnode[AliTRDgeometry::kNlayer - il - 1] - epsilonR;
      pointsInward[il].SetXYZ(xyz);
      xyz[0] = xAnode[il] + epsilonR;
      pointsOutward[il].SetXYZ(xyz);
    }
    /*for(Int_t ipt = 0; ipt < AliTRDgeometry::kNlayer; ipt++)
      printf("%d. X = %f\n", ipt, points[ipt].GetX());*/
    // Kalman inwards
    AliTRDtrackerV1::FitKalman(&copyTrack, NULL, kFALSE, 6, pointsInward);
    memcpy(points, pointsInward, sizeof(AliTrackPoint) * 6); // Preliminary store the inward results in the Array points
    // Kalman outwards
    AliTRDtrackerV1::FitKalman(&copyTrack, NULL, kTRUE, 6, pointsInward);
    memcpy(points, pointsOutward, sizeof(AliTrackPoint) * AliTRDgeometry::kNlayer);
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
    if((z + deltaZ > zmin && z - deltaZ < zmax) && (y + deltaY > ymin && y - deltaY < ymax)) nFindable++;
      if(DebugLevel() > 3){
        Double_t posTracklet[2] = {tracklet ? tracklet->GetYfit(0) : 0, tracklet ? tracklet->GetZfit(0) : 0};
        Int_t hasTracklet = tracklet ? 1 : 0;
        (*DebugStream())   << "FindableTracklets"
          << "layer="     << il
          << "ytracklet=" << posTracklet[0]
          << "ytrack="    << y
          << "ztracklet=" << posTracklet[1]
          << "ztrack="    << z
          << "tracklet="  << hasTracklet
          << "\n";
      }
  }
  
  h->Fill(nFindable > 0 ? TMath::Min(nFound/static_cast<Double_t>(nFindable), 1.) : 1);
  AliDebug(2, Form("Findable[Found]: %d[%d|%f]", nFindable, nFound, nFound/static_cast<Float_t>(nFindable > 0 ? nFindable : 1)));
  return h;
}


//_______________________________________________________
TH1 *AliTRDcheckDET::PlotChi2(const AliTRDtrackV1 *track){
  //
  // Plot the chi2 of the track
  //
  if(track) fkTrack = track;
  if(!fkTrack){
    AliDebug(4, "No Track defined.");
    return NULL;
  }
  TH1 *h = NULL;
  if(!(h = dynamic_cast<TH2S*>(fContainer->At(kChi2)))) {
    AliWarning("No Histogram defined.");
    return NULL;
  }
  Int_t n = fkTrack->GetNumberOfTracklets();
  if(!n) return NULL;

  h->Fill(n, fkTrack->GetChi2()/n);
  return h;
}


//_______________________________________________________
TH1 *AliTRDcheckDET::PlotPHt(const AliTRDtrackV1 *track){
  //
  // Plot the average pulse height
  //
  if(track) fkTrack = track;
  if(!fkTrack){
    AliDebug(4, "No Track defined.");
    return NULL;
  }
  TProfile *h = NULL; TH2F *phs2D = NULL;
  if(!(h = dynamic_cast<TProfile *>(((TObjArray*)(fContainer->At(kPH)))->At(0)))){
    AliWarning("No Histogram defined.");
    return NULL;
  }
  if(!(phs2D = dynamic_cast<TH2F *>(((TObjArray*)(fContainer->At(kPH)))->At(2)))){
    AliWarning("2D Pulse Height histogram not defined. Histogramm cannot be filled");
  }
  AliTRDseedV1 *tracklet = NULL;
  AliTRDcluster *c = NULL;
  for(Int_t itl = 0; itl < AliTRDgeometry::kNlayer; itl++){
    if(!(tracklet = fkTrack->GetTracklet(itl)) || !tracklet->IsOK())continue;
    Int_t crossing = Int_t(tracklet->IsRowCross());
    Int_t detector = tracklet->GetDetector();
    tracklet->ResetClusterIter();
    while((c = tracklet->NextCluster())){
      if(!IsUsingClustersOutsideChamber() && !c->IsInChamber()) continue;
      Int_t localtime        = c->GetLocalTimeBin();
      Double_t absoluteCharge = TMath::Abs(c->GetQ());
      h->Fill(localtime, absoluteCharge);
      phs2D->Fill(localtime, absoluteCharge); 
      if(DebugLevel() > 3){
        Int_t inChamber = c->IsInChamber() ? 1 : 0;
        Double_t distance[2];
        GetDistanceToTracklet(distance, tracklet, c);
        Float_t theta = TMath::ATan(tracklet->GetZref(1));
        Float_t phi = TMath::ATan(tracklet->GetYref(1));
        AliExternalTrackParam *trdPar = fkTrack->GetTrackIn();
        Float_t momentumMC = 0, momentumRec = trdPar ? trdPar->P() : track->P(); // prefer Track Low
        Int_t pdg = 0;
        Int_t kinkIndex = fkESD ? fkESD->GetKinkIndex() : 0;
        UShort_t TPCncls = fkESD ? fkESD->GetTPCncls() : 0;
        if(fkMC){
          if(fkMC->GetTrackRef()) momentumMC = fkMC->GetTrackRef()->P();
          pdg = fkMC->GetPDG();
        }
        (*DebugStream()) << "PHt"
          << "Detector="	<< detector
          << "crossing="	<< crossing
          << "inChamber=" << inChamber
          << "Timebin="		<< localtime
          << "Charge="		<< absoluteCharge
          << "momentumMC="	<< momentumMC
          << "momentumRec="	<< momentumRec
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
TH1 *AliTRDcheckDET::PlotPHx(const AliTRDtrackV1 *track){
  //
  // Plots the average pulse height vs the distance from the anode wire
  // (plus const anode wire offset)
  //
  if(track) fkTrack = track;
  if(!fkTrack){
     AliDebug(4, "No Track defined.");
     return NULL;
  }
  TProfile *h = NULL;
  if(!(h = dynamic_cast<TProfile *>(((TObjArray*)(fContainer->At(kPH)))->At(1)))){
    AliWarning("No Histogram defined.");
    return NULL;
  }
  AliTRDseedV1 *tracklet(NULL);
  AliTRDcluster *c(NULL);
  Double_t xd(0.), dqdl(0.);
  TVectorD vq(AliTRDseedV1::kNtb), vxd(AliTRDseedV1::kNtb), vdqdl(AliTRDseedV1::kNtb);
  for(Int_t itl = 0; itl < AliTRDgeometry::kNlayer; itl++){
    if(!(tracklet = fkTrack->GetTracklet(itl)) || !(tracklet->IsOK())) continue;
    Int_t det(tracklet->GetDetector());
    Bool_t rc(tracklet->IsRowCross());
    for(Int_t ic(0); ic<AliTRDseedV1::kNtb; ic++){
      Bool_t kFIRST(kFALSE);
      if(!(c = tracklet->GetClusters(ic))){
         if(!(c = tracklet->GetClusters(AliTRDseedV1::kNtb+ic))) continue;
      } else kFIRST=kTRUE;
      if(!IsUsingClustersOutsideChamber() && !c->IsInChamber()) continue;
      xd = tracklet->GetX0() - c->GetX(); vxd[ic] = xd;
      dqdl=tracklet->GetdQdl(ic); vdqdl[ic] = dqdl;
      vq[ic]=c->GetQ();
      if(kFIRST && (c = tracklet->GetClusters(AliTRDseedV1::kNtb+ic))) vq[ic]+=c->GetQ();
      h->Fill(xd, dqdl);
    }
    if(DebugLevel() > 3){
      (*DebugStream()) << "PHx"
        << "det="  << det
        << "rc="   << rc
        << "xd="   << &vxd
        << "q="    << &vq
        << "dqdl=" << &vdqdl
        << "\n";
    }
  }  
  return h;
}

//_______________________________________________________
TH1 *AliTRDcheckDET::PlotChargeCluster(const AliTRDtrackV1 *track){
  //
  // Plot the cluster charge
  //
  if(track) fkTrack = track;
  if(!fkTrack){
    AliDebug(4, "No Track defined.");
    return NULL;
  }
  TH1 *h = NULL;
  if(!(h = dynamic_cast<TH1F *>(fContainer->At(kChargeCluster)))){
    AliWarning("No Histogram defined.");
    return NULL;
  }
  AliTRDseedV1 *tracklet = NULL;
  AliTRDcluster *c = NULL;
  for(Int_t itl = 0; itl < AliTRDgeometry::kNlayer; itl++){
    if(!(tracklet = fkTrack->GetTracklet(itl)) || !tracklet->IsOK())continue;
    for(Int_t ic(0); ic < AliTRDseedV1::kNtb; ic++){
      Bool_t kFIRST(kFALSE);
      if(!(c = tracklet->GetClusters(ic))) {
        if(!(c = tracklet->GetClusters(AliTRDseedV1::kNtb+ic))) continue;
      } else kFIRST = kTRUE;
      Float_t q(c->GetQ());
      if(kFIRST && (c = tracklet->GetClusters(AliTRDseedV1::kNtb+ic))) q+=c->GetQ();
      h->Fill(q);
    }
  }
  return h;
}

//_______________________________________________________
TH1 *AliTRDcheckDET::PlotChargeTracklet(const AliTRDtrackV1 *track){
  //
  // Plot the charge deposit per chamber
  //
  if(track) fkTrack = track;
  if(!fkTrack){
    AliDebug(4, "No Track defined.");
    return NULL;
  }
  TH1 *h = NULL;
  if(!(h = dynamic_cast<TH1F *>(fContainer->At(kChargeTracklet)))){
    AliWarning("No Histogram defined.");
    return NULL;
  }
  AliTRDseedV1 *tracklet = NULL;
  AliTRDcluster *c = NULL;
  Double_t qTot = 0;
  Int_t nTracklets =fkTrack->GetNumberOfTracklets();
  for(Int_t itl(0); itl < AliTRDgeometry::kNlayer; itl++){
    if(!(tracklet = fkTrack->GetTracklet(itl)) || !tracklet->IsOK()) continue;
    qTot = 0.;
    for(Int_t ic = AliTRDseedV1::kNclusters; ic--;){
      if(!(c = tracklet->GetClusters(ic))) continue;
      qTot += TMath::Abs(c->GetQ());
    }
    h->Fill(qTot);
    if(DebugLevel() > 3){
      Int_t crossing = (Int_t)tracklet->IsRowCross();
      Int_t detector = tracklet->GetDetector();
      Float_t theta = TMath::ATan(tracklet->GetZfit(1));
      Float_t phi = TMath::ATan(tracklet->GetYfit(1));
      Float_t momentum = 0.;
      Int_t pdg = 0;
      Int_t kinkIndex = fkESD ? fkESD->GetKinkIndex() : 0;
      UShort_t nclsTPC = fkESD ? fkESD->GetTPCncls() : 0;
      if(fkMC){
	      if(fkMC->GetTrackRef()) momentum = fkMC->GetTrackRef()->P();
        pdg = fkMC->GetPDG();
      }
      (*DebugStream()) << "ChargeTracklet"
        << "Detector="  << detector
        << "crossing="  << crossing
        << "momentum="	<< momentum
        << "nTracklets="<< nTracklets
        << "pdg="				<< pdg
        << "theta="			<< theta
        << "phi="				<< phi
        << "kinkIndex="	<< kinkIndex
        << "TPCncls="		<< nclsTPC
        << "QT="        << qTot
        << "\n";
    }
  }
  return h;
}

//_______________________________________________________
TH1 *AliTRDcheckDET::PlotNTracksSector(const AliTRDtrackV1 *track){
  //
  // Plot the number of tracks per Sector
  //
  if(track) fkTrack = track;
  if(!fkTrack){
    AliDebug(4, "No Track defined.");
    return NULL;
  }
  TH1 *h = NULL;
  if(!(h = dynamic_cast<TH1F *>(fContainer->At(kNtracksSector)))){
    AliWarning("No Histogram defined.");
    return NULL;
  }

  // TODO we should compare with
  // sector = Int_t(track->GetAlpha() / AliTRDgeometry::GetAlpha());

  AliTRDseedV1 *tracklet = NULL;
  Int_t sector = -1;
  for(Int_t itl = 0; itl < AliTRDgeometry::kNlayer; itl++){
    if(!(tracklet = fkTrack->GetTracklet(itl)) || !tracklet->IsOK()) continue;
    sector = static_cast<Int_t>(tracklet->GetDetector()/AliTRDgeometry::kNdets);
    break;
  }
  h->Fill(sector);
  return h;
}


//________________________________________________________
void AliTRDcheckDET::SetRecoParam(AliTRDrecoParam *r)
{

  fReconstructor->SetRecoParam(r);
}

//________________________________________________________
void AliTRDcheckDET::GetDistanceToTracklet(Double_t *dist, AliTRDseedV1 * const tracklet, AliTRDcluster * const c)
{
  Float_t x = c->GetX();
  dist[0] = c->GetY() - tracklet->GetYat(x);
  dist[1] = c->GetZ() - tracklet->GetZat(x);
}

//________________________________________________________
void AliTRDcheckDET::GetEtaPhiAt(AliExternalTrackParam *track, Double_t x, Double_t &eta, Double_t &phi){
  //
  // Get phi and eta at a given radial position
  // 
  AliExternalTrackParam workpar(*track);

  Double_t posLocal[3];
  Bool_t sucPos = workpar.GetXYZAt(x, fEventInfo->GetRunInfo()->GetMagneticField(), posLocal);
  Double_t sagPhi = sucPos ? TMath::ATan2(posLocal[1], posLocal[0]) : 0.;
  phi = sagPhi;
  eta = workpar.Eta();
}


//_______________________________________________________
TH1* AliTRDcheckDET::MakePlotChi2()
{
// Plot chi2/track normalized to number of degree of freedom 
// (tracklets) and compare with the theoretical distribution.
// 
// Alex Bercuci <A.Bercuci@gsi.de>

  return NULL;

  TH2S *h2 = (TH2S*)fContainer->At(kChi2);
  TF1 f("fChi2", "[0]*pow(x, [1]-1)*exp(-0.5*x)", 0., 50.);
  f.SetParLimits(1,1, 1e100);
  TLegend *leg = new TLegend(.7,.7,.95,.95);
  leg->SetBorderSize(1); leg->SetHeader("Tracklets per Track");
  TH1D *h1 = NULL;
  Bool_t kFIRST = kTRUE;
  for(Int_t il=1; il<=h2->GetNbinsX(); il++){
    h1 = h2->ProjectionY(Form("pyChi2%d", il), il, il);
    if(h1->Integral()<50) continue;
    h1->Scale(1./h1->Integral());
    h1->SetMarkerStyle(7);h1->SetMarkerColor(il);
    h1->SetLineColor(il);h1->SetLineStyle(2);
    f.SetParameter(1, .5*il);f.SetLineColor(il);
    h1->Fit(&f, "QW+", kFIRST ? "pc": "pcsame");
    leg->AddEntry(h1, Form("%d", il), "l");
    if(kFIRST){
      h1->GetXaxis()->SetRangeUser(0., 25.);
    }
    kFIRST = kFALSE;
  }
  leg->Draw();
  gPad->SetLogy();
  return h1;
}


//________________________________________________________
TH1* AliTRDcheckDET::MakePlotNTracklets(){
  //
  // Make nice bar plot of the number of tracklets in each method
  //
  TH2F *tmp = (TH2F *)fContainer->FindObject("hNtlsBAR");
  TH1D *hBAR = tmp->ProjectionY();
  TH1F *hSTA = (TH1F *)fContainer->FindObject("hNtlsSTA");
  TH1F *hCON = (TH1F *)fContainer->FindObject("hNtls");
  TLegend *leg = new TLegend(0.13, 0.75, 0.39, 0.89);
  leg->SetBorderSize(1);
  leg->SetFillColor(0);

  Float_t scale = hCON->Integral();
  hCON->Scale(100./scale);
  hCON->SetFillColor(kRed);hCON->SetLineColor(kRed);
  hCON->SetBarWidth(0.2);
  hCON->SetBarOffset(0.6);
  hCON->SetStats(kFALSE);
  hCON->GetYaxis()->SetRangeUser(0.,40.);
  hCON->GetYaxis()->SetTitleOffset(1.2);
  hCON->Draw("bar1"); leg->AddEntry(hCON, "Total", "f");
  hCON->SetMaximum(55.);

  hBAR->Scale(100./scale);
  hBAR->SetFillColor(kGreen);hBAR->SetLineColor(kGreen);
  hBAR->SetBarWidth(0.2);
  hBAR->SetBarOffset(0.2);
  hBAR->SetTitle("");
  hBAR->SetStats(kFALSE);
  hBAR->GetYaxis()->SetRangeUser(0.,40.);
  hBAR->GetYaxis()->SetTitleOffset(1.2);
  hBAR->Draw("bar1same"); leg->AddEntry(hBAR, "Barrel", "f");

  hSTA->Scale(100./scale);
  hSTA->SetFillColor(kBlue);hSTA->SetLineColor(kBlue);
  hSTA->SetBarWidth(0.2);
  hSTA->SetBarOffset(0.4);
  hSTA->SetTitle("");
  hSTA->SetStats(kFALSE);
  hSTA->GetYaxis()->SetRangeUser(0.,40.);
  hSTA->GetYaxis()->SetTitleOffset(1.2);
  hSTA->Draw("bar1same"); leg->AddEntry(hSTA, "Stand Alone", "f");
  leg->Draw();
  gPad->Update();
  return hCON;
}

//________________________________________________________
void AliTRDcheckDET::MakePlotnTrackletsVsP(){
  //
  // Plot abundance of tracks with number of tracklets as function of momentum
  //




  Color_t color[AliTRDgeometry::kNlayer] = {kBlue, kOrange, kBlack, kGreen, kCyan, kRed};
  TH1 *h(NULL); TGraphErrors *g[AliTRDgeometry::kNlayer];
  for(Int_t itl(0); itl<AliTRDgeometry::kNlayer; itl++){
    g[itl] = new TGraphErrors();
    g[itl]->SetLineColor(color[itl]);
    g[itl]->SetMarkerColor(color[itl]);
    g[itl]->SetMarkerStyle(20 + itl);
  }

  TH2 *hBar = (TH2F *)fContainer->FindObject("hNtlsBAR");
  TAxis *ax(hBar->GetXaxis());
  Int_t np(ax->GetNbins());
  for(Int_t ipBin(1); ipBin<np; ipBin++){
    h = hBar->ProjectionY("npBin", ipBin, ipBin);
    if(!Int_t(h->Integral())) continue;
    h->Scale(100./h->Integral());
    Float_t p(ax->GetBinCenter(ipBin)); 
    Float_t dp(ax->GetBinWidth(ipBin)); 
    Int_t ip(g[0]->GetN());
    for(Int_t itl(AliTRDgeometry::kNlayer); itl--;){
      g[itl]->SetPoint(ip, p, h->GetBinContent(itl+1));
      g[itl]->SetPointError(ip, dp/3.46, h->GetBinError(itl+1));
    }
  }

  TLegend *leg = new TLegend(0.76, 0.6, 1., 0.9);
  leg->SetBorderSize(0);
  leg->SetHeader("Tracklet/Track");
  leg->SetFillStyle(0);
  h = hBar->ProjectionX("npxBin"); h->Reset();
  h->SetTitle("");
  h->GetYaxis()->SetRangeUser(1., 99.);
  h->GetYaxis()->SetMoreLogLabels();
  h->GetYaxis()->CenterTitle();
  h->GetYaxis()->SetTitleOffset(1.2);
  h->SetYTitle("Prob. [%]");
  h->GetXaxis()->SetRangeUser(0.4, 12.);
  h->GetXaxis()->SetMoreLogLabels();
  h->GetXaxis()->CenterTitle();
  h->Draw("p");
  for(Int_t itl(AliTRDgeometry::kNlayer); itl--;){
    g[itl]->Draw("pc");
    leg->AddEntry(g[itl], Form("n = %d", itl+1),"pl");
  }

  leg->Draw();
  gPad->SetLogx();gPad->SetLogy();
}

//________________________________________________________
Bool_t AliTRDcheckDET::MakePlotPulseHeight(){
  //
  // Create Plot of the Pluse Height Spectrum
  //
  TCanvas *output = gPad->GetCanvas();
  output->Divide(2);
  output->cd(1);
  TH1 *h, *h1, *h2;
  TObjArray *arr = (TObjArray*)fContainer->FindObject("<PH>");
  h = (TH1F*)arr->At(0);
  h->SetMarkerStyle(24);
  h->SetMarkerColor(kBlack);
  h->SetLineColor(kBlack);
  h->GetYaxis()->SetTitleOffset(1.5);
  h->Draw("e1");
  // Trending for the pulse height: plateau value, slope and timebin of the maximum
  TLinearFitter fit(1,"pol1");
  Double_t time = 0.;
  for(Int_t itime = 10; itime <= 20; itime++){
    time = static_cast<Double_t>(itime);
    fit.AddPoint(&time, h->GetBinContent(itime + 1), h->GetBinError(itime + 1));
  }
  fit.Eval();
  Double_t plateau = fit.GetParameter(0) + 12 * fit.GetParameter(1);
  Double_t slope = fit.GetParameter(1);
  PutTrendValue("PHplateau", plateau);
  PutTrendValue("PHslope", slope);
  PutTrendValue("PHamplificationPeak", static_cast<Double_t>(h->GetMaximumBin()-1));
  AliDebug(1, Form("plateau %f, slope %f, MaxTime %f", plateau, slope, static_cast<Double_t>(h->GetMaximumBin()-1)));
//   copy the second histogram in a new one with the same x-dimension as the phs with respect to time
  h1 = (TH1F *)arr->At(1);
  h2 = new TH1F("hphs1","Average PH", 31, -0.5, 30.5);
  for(Int_t ibin = h1->GetXaxis()->GetFirst(); ibin < h1->GetNbinsX(); ibin++) 
    h2->SetBinContent(ibin, h1->GetBinContent(ibin));
  h2->SetMarkerStyle(22);
  h2->SetMarkerColor(kBlue);
  h2->SetLineColor(kBlue);
  h2->Draw("e1same");
  gPad->Update();
//   create axis according to the histogram dimensions of the original second histogram
  TGaxis *axis = new TGaxis(gPad->GetUxmin(),
                    gPad->GetUymax(),
                    gPad->GetUxmax(),
                    gPad->GetUymax(),
                    -0.08, 4.88, 510,"-L");
  axis->SetLineColor(kBlue);
  axis->SetLabelColor(kBlue);
  axis->SetTextColor(kBlue);
  axis->SetTitle("x_{0}-x_{c} [cm]");
  axis->Draw();

  output->cd(2);
  TH2 *ph2d = (TH2F *)arr->At(2);
  ph2d->GetYaxis()->SetTitleOffset(1.8);
  ph2d->SetStats(kFALSE);
  ph2d->Draw("colz");
  return kTRUE;
}

//________________________________________________________
void AliTRDcheckDET::MakePlotMeanClustersLayer(){
  //
  // Create Summary plot for the mean number of clusters per layer
  //
  TCanvas *output = gPad->GetCanvas();
  output->Divide(3,2);
  TObjArray *histos = (TObjArray *)fContainer->At(kNclustersLayer);
  if(!histos){
    AliWarning("Histos for each layer not found");
    return;
  }
  TProfile2D *hlayer = NULL;
  for(Int_t ily = 0; ily < AliTRDgeometry::kNlayer; ily++){
    hlayer = dynamic_cast<TProfile2D *>(histos->At(ily));
    output->cd(ily + 1);
    gPad->SetGrid(0,0);
    hlayer->Draw("colz");
  }
}

//________________________________________________________
Bool_t AliTRDcheckDET::MakeBarPlot(TH1 *histo, Int_t color){
  //
  // Draw nice bar plots
  //
  if(!histo->GetEntries()) return kFALSE;
  histo->Scale(100./histo->Integral());
  histo->SetFillColor(color);
  histo->SetBarOffset(.2);
  histo->SetBarWidth(.6);
  histo->Draw("bar1");
  return kTRUE;
}
