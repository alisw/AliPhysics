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

#include <TArrayD.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH1I.h>
#include <TH2F.h>
#include <TH3S.h>
#include <TH3F.h>
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
#include "AliTRDinfoGen.h"
#include "AliTRDcheckDET.h"
#include "AliTRDpwg1Helper.h"

#include <cstdio>
#include <iostream>

ClassImp(AliTRDcheckDET)

const Color_t AliTRDcheckDET::fkColorsCentrality[AliTRDeventInfo::kCentralityClasses] = {kTeal, kOrange, kGreen, kBlue ,kRed};

//_______________________________________________________
AliTRDcheckDET::AliTRDcheckDET():
  AliTRDrecoTask()
  ,fCentralityClass(-1)
  ,fTriggerNames(NULL)
  ,fFlags(0)
{
  //
  // Default constructor
  //
  SetNameTitle("TRDcheckDET", "Basic TRD data checker");
}

//_______________________________________________________
AliTRDcheckDET::AliTRDcheckDET(char* name):
  AliTRDrecoTask(name, "Basic TRD data checker")
  ,fCentralityClass(-1)
  ,fTriggerNames(NULL)
  ,fFlags(0)
{
  //
  // Default constructor
  //
  InitFunctorList();
}


//_______________________________________________________
AliTRDcheckDET::~AliTRDcheckDET(){
  //
  // Destructor
  // 
  if(fTriggerNames) delete fTriggerNames;
}

//_______________________________________________________
void AliTRDcheckDET::UserCreateOutputObjects(){
  //
  // Create Output Objects
  //
  AliTRDrecoTask::UserCreateOutputObjects();
  if(!fTriggerNames) fTriggerNames = new TMap();
}

//_______________________________________________________
void AliTRDcheckDET::UserExec(Option_t *opt){
  //
  // Execution function
  // Filling TRD quality histos
  //
  AliDebug(2, Form("EventInfo[%p] Header[%p]", (void*)fEvent, (void*)(fEvent?fEvent->GetEventHeader():NULL)));
  if(fEvent) fCentralityClass = fEvent->GetCentrality();
  else fCentralityClass = -1;  // Assume pp
 
  AliTRDrecoTask::UserExec(opt);  

  TH1F *histo(NULL); AliTRDtrackInfo *fTrackInfo(NULL); Int_t nTracks(0);		// Count the number of tracks per event
  for(Int_t iti = 0; iti < fTracks->GetEntriesFast(); iti++){
    if(!fTracks->UncheckedAt(iti)) continue;
    if(!(fTrackInfo = dynamic_cast<AliTRDtrackInfo *>(fTracks->UncheckedAt(iti)))) continue;
    if(!fTrackInfo->GetTrack()) continue;
    nTracks++;
  }
  if(nTracks)
    if((histo = dynamic_cast<TH1F *>(fContainer->UncheckedAt(kNtracksEvent)))) histo->Fill(nTracks);

  if(!fEvent->GetEventHeader()) return; // For trigger statistics event header is essential
  Int_t triggermask = fEvent->GetEventHeader()->GetTriggerMask();
  TString triggername =  fEvent->GetRunInfo()->GetFiredTriggerClasses(triggermask);
  AliDebug(6, Form("Trigger cluster: %d, Trigger class: %s\n", triggermask, triggername.Data()));
  if((histo = dynamic_cast<TH1F *>(fContainer->UncheckedAt(kNeventsTrigger)))) histo->Fill(triggermask);

  if(nTracks){
    if((histo = dynamic_cast<TH1F *>(fContainer->UncheckedAt(kNeventsTriggerTracks)))) histo->Fill(triggermask);
  }
  if(triggermask <= 20 && !fTriggerNames->FindObject(Form("%d", triggermask))){
    fTriggerNames->Add(new TObjString(Form("%d", triggermask)), new TObjString(triggername));
    // also set the label for both histograms
    if((histo = dynamic_cast<TH1F *>(fContainer->UncheckedAt(kNeventsTriggerTracks))))
      histo->GetXaxis()->SetBinLabel(histo->FindBin(triggermask), triggername);
    if((histo = dynamic_cast<TH1F *>(fContainer->UncheckedAt(kNeventsTrigger))))
      histo->GetXaxis()->SetBinLabel(histo->FindBin(triggermask), triggername);
  }
}


//_______________________________________________________
//_______________________________________________________
Bool_t AliTRDcheckDET::PostProcess(){
  //
  // Do Postprocessing (for the moment set the number of Reference histograms)
  //
  
  TH1F *h(NULL), *h1(NULL);

  // Calculate of the trigger clusters purity
  if((h  = dynamic_cast<TH1F *>(fContainer->FindObject("hEventsTrigger"))) &&
     (h1 = dynamic_cast<TH1F *>(fContainer->FindObject("hEventsTriggerTracks")))) {
    h1->Divide(h);
    Float_t purities[20], val = 0; memset(purities, 0, 20*sizeof(Float_t));
    TString triggernames[20];
    Int_t nTriggerClasses = 0;
    for(Int_t ibin = 1; ibin <= h->GetNbinsX(); ibin++){
      if((val = h1->GetBinContent(ibin))){
        purities[nTriggerClasses] = val;
        triggernames[nTriggerClasses] = h1->GetXaxis()->GetBinLabel(ibin);
        nTriggerClasses++;
      }
    }

    if((h = dynamic_cast<TH1F *>(fContainer->UncheckedAt(kTriggerPurity)))){
      TAxis *ax = h->GetXaxis();
      for(Int_t itrg = 0; itrg < nTriggerClasses; itrg++){
        h->Fill(itrg, purities[itrg]);
        ax->SetBinLabel(itrg+1, triggernames[itrg].Data());
      }
      ax->SetRangeUser(-0.5, nTriggerClasses+.5);
      h->GetYaxis()->SetRangeUser(0,1);
    }
  }

  // track status
  if((h=dynamic_cast<TH1F*>(fContainer->At(kTrackStatus)))){
    Float_t ok = h->GetBinContent(1);
    Int_t nerr = h->GetNbinsX();
    for(Int_t ierr=nerr; ierr--;){
      h->SetBinContent(ierr+1, ok>0.?1.e2*h->GetBinContent(ierr+1)/ok:0.);
    }
    h->SetBinContent(1, 0.);
  }
  // tracklet status
  
  TObjArray *arr(NULL);
  if(( arr = dynamic_cast<TObjArray*>(fContainer->UncheckedAt(kTrackletStatus)))){
    for(Int_t ily = AliTRDgeometry::kNlayer; ily--;){
      if(!(h=dynamic_cast<TH1F*>(arr->At(ily)))) continue;
      Float_t okB = h->GetBinContent(1);
      Int_t nerrB = h->GetNbinsX();
      for(Int_t ierr=nerrB; ierr--;){
        h->SetBinContent(ierr+1, okB>0.?1.e2*h->GetBinContent(ierr+1)/okB:0.);
      }
      h->SetBinContent(1, 0.);
    }
  }
  fNRefFigures = 17;

  return kTRUE;
}

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
  cOut->SaveAs("TRD_TrackQuality.gif");
  delete cOut;

  // Second Plot: PHS
  cOut = new TCanvas(Form("summary%s2", GetName()), Form("Summary 2 for task %s", GetName()), 1024, 512);
  cOut->cd(); GetRefFigure(kFigPH);
  cOut->SaveAs("TRD_PH.gif"); 
  delete cOut;

  // Third Plot: Mean Number of clusters as function of eta, phi and layer
   cOut = new TCanvas(Form("summary%s3", GetName()), Form("Summary 3 for task %s", GetName()), 1024, 768);
  cOut->cd(); MakePlotMeanClustersLayer();
  cOut->SaveAs("TRD_MeanNclusters.gif"); 
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
		MakePlotNclustersTrack();
		return kTRUE;
	case kFigNclustersTracklet:
		MakePlotNclustersTracklet();
		return kTRUE;
	case kFigNtrackletsTrack:
		h=MakePlotNTracklets();
		if(h){
			PutTrendValue("NTrackletsTrack", h->GetMean());
			PutTrendValue("NTrackletsTrackRMS", h->GetRMS());
		}
		return kTRUE;
	case kFigNTrackletsP:
		MakePlotnTrackletsVsP();
		return kTRUE;
	case kFigNtrackletsCross:
		h = ProjectCentrality((TH2*)fContainer->FindObject("hNtlsCross"), -1);
		if(!MakeBarPlot(h, kRed)) break;
		PutTrendValue("NTrackletsCross", h->GetMean());
		PutTrendValue("NTrackletsCrossRMS", h->GetRMS());
		return kTRUE;
	case kFigNtrackletsFindable:
		h = ProjectCentrality((TH2*)fContainer->FindObject("hNtlsFindable"), -1);
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
		MakePlotChi2();
		return kTRUE;
	case kFigPH:
		gPad->SetMargin(0.125, 0.015, 0.1, 0.1);
		MakePlotPulseHeight();
		gPad->SetLogy(0);
		return kTRUE;
	case kFigChargeCluster:
		h = ProjectCentrality((TH2*)fContainer->FindObject("hQcl"), -1);
		h->Draw("c");
		gPad->SetLogy(1);
		PutTrendValue("ChargeCluster", h->GetMaximumBin());
		PutTrendValue("ChargeClusterRMS", h->GetRMS());
		return kTRUE;
	case kFigChargeTracklet:
		MakePlotTrackletCharge();
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
  fContainer->SetOwner(kTRUE);

  // Register Histograms
  TH1 * h = NULL;
  TH2 * h2 = NULL;      // Pointer for two dimensional histograms
  TH3 * h3 = NULL;      // Pointer for tree dimensional histograms
  TAxis *ax = NULL;
  if(!(h2 = (TH2F *)gROOT->FindObject("hNcls"))){
    h2 = new TH2F("hNcls", "N_{clusters}/track;N_{clusters};Centrality;Entries", 181, -0.5, 180.5, AliTRDeventInfo::kCentralityClasses + 1, -1.5, AliTRDeventInfo::kCentralityClasses - 0.5);
  } else h2->Reset();
  fContainer->AddAt(h2, kNclustersTrack);
 
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

  if(!(h2 = (TH2F *)gROOT->FindObject("hNclTls"))){
    h2 = new TH2F("hNclTls","N_{clusters}/tracklet;N_{clusters};Entries", 51, -0.5, 50.5, AliTRDeventInfo::kCentralityClasses + 1, -1.5, AliTRDeventInfo::kCentralityClasses - 0.5);
  } else h2->Reset();
  fContainer->AddAt(h2, kNclustersTracklet);

  if(!(h2 = (TH2F *)gROOT->FindObject("hNtls"))){
    h2 = new TH2F("hNtls", "N_{tracklets}/track;N^{tracklet};Centrality;freq.[%]", AliTRDgeometry::kNlayer, 0.5, 6.5, AliTRDeventInfo::kCentralityClasses + 1, -1.5, AliTRDeventInfo::kCentralityClasses - 0.5);
  } else h2->Reset();
  fContainer->AddAt(h2, kNtrackletsTrack);

  if(!(h = (TH2F *)gROOT->FindObject("htlsSTA"))){
    h = new TH2F("hNtlsSTA", "#splitline{N_{tracklets}/track}{Stand Alone};N^{tracklet};Centrality;freq.[%]", AliTRDgeometry::kNlayer, 0.5, 6.5, AliTRDeventInfo::kCentralityClasses + 1, -1.5, AliTRDeventInfo::kCentralityClasses - 0.5);
  }
  fContainer->AddAt(h, kNtrackletsSTA);

  // Binning for momentum dependent tracklet Plots
  const Int_t kNp(30);
  Float_t p=0.2;
  Float_t binsP[kNp+1], binsTrklt[AliTRDgeometry::kNlayer+1], binsCent[AliTRDeventInfo::kCentralityClasses+2];
  for(Int_t i=0;i<kNp+1; i++,p+=(TMath::Exp(i*i*.001)-1.)) binsP[i]=p;
  for(Int_t il = 0; il <= AliTRDgeometry::kNlayer; il++) binsTrklt[il] = 0.5 + il;
  for(Int_t icent = -1; icent < AliTRDeventInfo::kCentralityClasses + 1; icent++) binsCent[icent+1] = icent - 0.5;
  if(!(h3 = (TH3F *)gROOT->FindObject("htlsBAR"))){
    // Make tracklets for barrel tracks momentum dependent (if we do not exceed min and max values)
    h3 = new TH3F("hNtlsBAR", 
    "N_{tracklets}/track;p [GeV/c];N^{tracklet};freq. [%]",
    kNp, binsP, AliTRDgeometry::kNlayer, binsTrklt, AliTRDeventInfo::kCentralityClasses + 1, binsCent);
  }
  fContainer->AddAt(h3, kNtrackletsBAR);

  // 
  if(!(h2 = (TH2F *)gROOT->FindObject("hNtlsCross"))){
    h2 = new TH2F("hNtlsCross", "N_{tracklets}^{cross}/track;n_{row cross};Centrality;freq.[%]", 7, -0.5, 6.5, AliTRDeventInfo::kCentralityClasses + 1, -1.5, AliTRDeventInfo::kCentralityClasses - 0.5);
  } else h2->Reset();
  fContainer->AddAt(h2, kNtrackletsCross);

  if(!(h2 = (TH2F *)gROOT->FindObject("hNtlsFindable"))){
    h2 = new TH2F("hNtlsFindable", "Found/Findable Tracklets;r[a.u];Centrality;Entries" , 101, -0.005, 1.005,  AliTRDeventInfo::kCentralityClasses + 1, -1.5, AliTRDeventInfo::kCentralityClasses - 0.5);
  } else h2->Reset();
  fContainer->AddAt(h2, kNtrackletsFindable);

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
  if(!(h3 = (TH3F *)gROOT->FindObject("hPHx"))){
    h3 = new TH3F("hPHx", "<PH>(x);x[mm];Centrality;Charge[a.u.]", 31, -0.08, 4.88, 100, 0, 1024, AliTRDeventInfo::kCentralityClasses + 1, -1.5, AliTRDeventInfo::kCentralityClasses - 0.5);
  } else h3->Reset();
  arr->AddAt(h3, 0);
  if(!(h3 = (TH3F *)gROOT->FindObject("hPHt"))){
    h3 = new TH3F("hPHt", "<PH>(t);time[100ns];Centrality;Charge[a.u.]", 31, -0.5, 30.5, 100, 0, 1024, AliTRDeventInfo::kCentralityClasses + 1, -1.5, AliTRDeventInfo::kCentralityClasses - 0.5);
  } else h3->Reset();
  arr->AddAt(h3, 1);

  // Chi2 histos
  if(!(h3 = (TH3S*)gROOT->FindObject("hChi2"))){
    h3 = new TH3S("hChi2", "#chi^{2}/track;ndf;#chi^{2}/ndf;Centrality", AliTRDgeometry::kNlayer, .5, AliTRDgeometry::kNlayer+.5, 100, 0, 50, AliTRDeventInfo::kCentralityClasses + 1, -1.5, AliTRDeventInfo::kCentralityClasses - 0.5);
  } else h3->Reset();
  fContainer->AddAt(h3, kChi2);

  if(!(h2 = (TH2F *)gROOT->FindObject("hQcl"))){
    h2 = new TH2F("hQcl", "Q_{cluster};Charge[a.u.];Centrality;Entries", 200, 0, 1200, AliTRDeventInfo::kCentralityClasses + 1, -1.5, AliTRDeventInfo::kCentralityClasses - 0.5);
  }else h2->Reset();
  fContainer->AddAt(h2, kChargeCluster);

  if(!(h2 = (TH2F *)gROOT->FindObject("hQtrklt"))){
    h2 = new TH2F("hQtrklt", "Q_{tracklet};Charge[a.u.];Centrality;Entries", 6000, 0, 6000, AliTRDeventInfo::kCentralityClasses + 1, -1.5, AliTRDeventInfo::kCentralityClasses - 0.5);
  }else h2->Reset();
  fContainer->AddAt(h2, kChargeTracklet);


  if(!(h = (TH1F *)gROOT->FindObject("hEventsTrigger")))
    h = new TH1F("hEventsTrigger", "Trigger Class", 100, 0, 100);
  else h->Reset();
  
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

//_______________________________________________________
TH1 *AliTRDcheckDET::ProjectCentrality(TH2 *hIn, Int_t centralityBin){
  //
  // Project histogram to a given centrality Bin
  //
  if(!hIn) return NULL;
  if(centralityBin >= AliTRDeventInfo::kCentralityClasses) centralityBin = -1;
  Int_t binMin = centralityBin > -1 ? centralityBin + 2 : 0, 
        binMax = centralityBin > -1 ? centralityBin + 2 : -1;
  return hIn->ProjectionX(Form("%s_%d", hIn->GetName(), centralityBin), binMin, binMax);
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
  TH2 *h = NULL;
  TProfile2D *hlayer = NULL;
  Double_t eta = 0., phi = 0.;
  if(!(h = dynamic_cast<TH2F *>(fContainer->At(kNclustersTracklet)))){
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
    h->Fill(tracklet->GetN2(), fCentralityClass);
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
  TH2 *h = NULL;
  if(!(h = dynamic_cast<TH2F *>(fContainer->At(kNclustersTrack)))){
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
        << "Centrality="<< fCentralityClass
        << "crossing="  << crossing
        << "momentumMC="<< momentumMC
        << "momentumRec="<< momentumRec
        << "pdg="				<< pdg
        << "theta="			<< theta
        << "phi="				<< phi
        << "kinkIndex="	<< kinkIndex
        << "TPCncls="		<< nclsTPC
        << "TRDncls="   << n
        << "\n";
    }
  }
  h->Fill(nclusters, fCentralityClass);
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
  TH2 *h = NULL, *hSta = NULL; TH3 *hBarrel = NULL;
  if(!(h = dynamic_cast<TH2F *>(fContainer->At(kNtrackletsTrack)))){
    AliWarning("No Histogram defined.");
    return NULL;
  }
  Int_t nTracklets = fkTrack->GetNumberOfTracklets();
  h->Fill(nTracklets, fCentralityClass);
  if(!fkESD) return h;
  Int_t status = fkESD->GetStatus();

/*  printf("in/out/refit/pid: TRD[%d|%d|%d|%d]\n", status &AliESDtrack::kTRDin ? 1 : 0, status &AliESDtrack::kTRDout ? 1 : 0, status &AliESDtrack::kTRDrefit ? 1 : 0, status &AliESDtrack::kTRDpid ? 1 : 0);*/
  Double_t p = 0.;
  Int_t method = -1;    // to distinguish between stand alone and full barrel tracks in the debugging
  if((status & AliESDtrack::kTRDin) != 0){
    method = 1;
    // Full Barrel Track: Save momentum dependence
    if(!(hBarrel = dynamic_cast<TH3F *>(fContainer->At(kNtrackletsBAR)))){
      AliWarning("Method: Barrel.  Histogram not processed!");
      return NULL;
    }
    AliExternalTrackParam *par(fkTrack->GetTrackIn());
    if(!par){
      AliError("Input track params missing");
      return NULL;
    }
    p = par->P(); // p needed later in the debug streaming
    hBarrel->Fill(p, nTracklets, fCentralityClass);
  } else {
    // Stand alone Track: momentum dependence not usefull
    method = 0;
    if(!(hSta = dynamic_cast<TH2F *>(fContainer->At(kNtrackletsSTA)))) {
      AliWarning("Method: StandAlone.  Histogram not processed!");
      return NULL;
    }
    hSta->Fill(nTracklets, fCentralityClass);
  }

  if(DebugLevel() > 2){
    AliTRDseedV1 *tracklet = NULL;
    AliTRDgeometry *geo(AliTRDinfoGen::Geometry());
    Int_t sector = -1, stack = -1, detector;
    for(Int_t itl = 0; itl < AliTRDgeometry::kNlayer; itl++){
      if(!(tracklet = fkTrack->GetTracklet(itl)) || !(tracklet->IsOK())) continue;
      detector = tracklet->GetDetector();
      sector = geo->GetSector(detector);
      stack = geo->GetStack(detector);
      break;
    }
    (*DebugStream()) << "NTrackletsTrack"
      << "Sector="      << sector
      << "Stack="       << stack
      << "Centrality="  << fCentralityClass
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
  TH2 *h = NULL;
  if(!(h = dynamic_cast<TH2F *>(fContainer->At(kNtrackletsCross)))){
    AliWarning("No Histogram defined.");
    return NULL;
  }

  Int_t ncross = 0;
  AliTRDseedV1 *tracklet = NULL;
  for(Int_t il = 0; il < AliTRDgeometry::kNlayer; il++){
    if(!(tracklet = fkTrack->GetTracklet(il)) || !tracklet->IsOK()) continue;

    if(tracklet->IsRowCross()) ncross++;
  }
  h->Fill(ncross, fCentralityClass);
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
  TH2 *h = NULL;
  if(!(h = dynamic_cast<TH2F *>(fContainer->At(kNtrackletsFindable)))){
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
      tracklet->SetReconstructor(AliTRDinfoGen::Reconstructor());
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
  AliTRDgeometry *geo(AliTRDinfoGen::Geometry());
  for(Int_t il = 0; il < AliTRDgeometry::kNlayer; il++){
    y = points[il].GetY();
    z = points[il].GetZ();
    if((stack = geo->GetStack(z, il)) < 0) continue; // Not findable
    pp = geo->GetPadPlane(il, stack);
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
  
  h->Fill((nFindable > 0 ? TMath::Min(nFound/static_cast<Double_t>(nFindable), 1.) : 1), fCentralityClass);
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
  TH3 *h = NULL;
  if(!(h = dynamic_cast<TH3S*>(fContainer->At(kChi2)))) {
    AliWarning("No Histogram defined.");
    return NULL;
  }
  Int_t n = fkTrack->GetNumberOfTracklets();
  if(!n) return NULL;

  h->Fill(n, fkTrack->GetChi2()/n, fCentralityClass);
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
  TH3F *h = NULL;
  if(!(h = dynamic_cast<TH3F *>(((TObjArray*)(fContainer->At(kPH)))->At(1)))){
    AliWarning("No Histogram defined.");
    return NULL;
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
      h->Fill(localtime, absoluteCharge, fCentralityClass);
      if(DebugLevel() > 3){
        Int_t inChamber = c->IsInChamber() ? 1 : 0;
        Double_t distance[2];
        GetDistanceToTracklet(distance, tracklet, c);
        Float_t theta = TMath::ATan(tracklet->GetZref(1));
        Float_t phi = TMath::ATan(tracklet->GetYref(1));
        AliExternalTrackParam *trdPar = fkTrack->GetTrackIn();
        Float_t momentumMC = 0, momentumRec = trdPar ? trdPar->P() : fkTrack->P(); // prefer Track Low
        Int_t pdg = 0;
        Int_t kinkIndex = fkESD ? fkESD->GetKinkIndex() : 0;
        UShort_t tpcNCLS = fkESD ? fkESD->GetTPCncls() : 0;
        if(fkMC){
          if(fkMC->GetTrackRef()) momentumMC = fkMC->GetTrackRef()->P();
          pdg = fkMC->GetPDG();
        }
        (*DebugStream()) << "PHt"
          << "Detector="	<< detector
          << "Centrality="<< fCentralityClass
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
          << "TPCncls="		<< tpcNCLS
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
  TH3 *h = NULL;
  if(!(h = dynamic_cast<TH3F *>(((TObjArray*)(fContainer->At(kPH)))->At(0)))){
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
      h->Fill(xd, dqdl, fCentralityClass);
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
  TH2 *h = NULL;
  if(!(h = dynamic_cast<TH2F *>(fContainer->At(kChargeCluster)))){
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
      h->Fill(q, fCentralityClass);
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
  TH2 *h = NULL;
  if(!(h = dynamic_cast<TH2F *>(fContainer->At(kChargeTracklet)))){
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
    h->Fill(qTot, fCentralityClass);
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
        << "Centrality="<< fCentralityClass
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
void AliTRDcheckDET::GetDistanceToTracklet(Double_t *dist, AliTRDseedV1 * const tracklet, AliTRDcluster * const c)
{
  Float_t x = c->GetX();
  dist[0] = c->GetY() - tracklet->GetYat(x);
  dist[1] = c->GetZ() - tracklet->GetZat(x);
}

//________________________________________________________
void AliTRDcheckDET::GetEtaPhiAt(const AliExternalTrackParam *track, Double_t x, Double_t &eta, Double_t &phi){
  //
  // Get phi and eta at a given radial position
  // 
  AliExternalTrackParam workpar(*track);

  Double_t posLocal[3];
  Bool_t sucPos = workpar.GetXYZAt(x, fEvent->GetRunInfo()->GetMagneticField(), posLocal);
  Double_t sagPhi = sucPos ? TMath::ATan2(posLocal[1], posLocal[0]) : 0.;
  phi = sagPhi;
  eta = workpar.Eta();
}


//_______________________________________________________
TH1* AliTRDcheckDET::MakePlotChi2() const
{
// Plot chi2/track normalized to number of degree of freedom 
// (tracklets) and compare with the theoretical distribution.
// 
// Alex Bercuci <A.Bercuci@gsi.de>

  return NULL;

/*  TH2S *h2 = (TH2S*)fContainer->At(kChi2);
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
  return h1;*/
}


//________________________________________________________
TH1* AliTRDcheckDET::MakePlotNTracklets(){
	//
	// Make nice bar plot of the number of tracklets in each method
	//
	TH3 *hTracklets3D = dynamic_cast<TH3F *>(fContainer->FindObject("hNtlsBAR"));
	if(!hTracklets3D){
		AliError("Tracklet Histogram not found");
		return NULL;
	}
	TH1 *hBAR = hTracklets3D->Project3D("y");
	hBAR->SetName("hBAR");
	hBAR->SetTitle("Number of Tracklets");
	hBAR->Scale(1./hBAR->Integral());
	hBAR->SetLineColor(kBlack);
	hBAR->SetLineWidth(2);
	hBAR->Draw();

	// Draw also centrality-dependent plots
	TH1 *hBarCent[AliTRDeventInfo::kCentralityClasses];
	Int_t nHistsCentrality = 0;
	for(Int_t icent = 0; icent < AliTRDeventInfo::kCentralityClasses; icent++){
		hTracklets3D->GetZaxis()->SetRange(icent+2, icent+2);
		hBarCent[icent] = hTracklets3D->Project3D("y");
		if(!(hBarCent[icent] && hBarCent[icent]->GetEntries())){
			delete hBarCent[icent];
			hBarCent[icent] = NULL;
			continue;
		}
		hBarCent[icent]->SetName(Form("hBarCent_%d", icent));
		hBarCent[icent]->SetTitle("Number of Tracklets");
		hBarCent[icent]->Scale(1./hBarCent[icent]->Integral());
		hBarCent[icent]->SetLineColor(fkColorsCentrality[icent]);
		hBarCent[icent]->Draw("same");
		nHistsCentrality++;
	}
	hTracklets3D->GetZaxis()->SetRange(0, hTracklets3D->GetNbinsZ());
	AliInfo(Form("Number of Centrality Classes: %d", nHistsCentrality));
	if(nHistsCentrality){
		TLegend *leg = new TLegend(0.5, 0.6, 0.89, 0.89);
		leg->SetBorderSize(0);
		leg->SetFillStyle(0);
		for(Int_t icent = 0; icent < AliTRDeventInfo::kCentralityClasses; icent++)
			if(hBarCent[icent]) leg->AddEntry(hBarCent[icent], Form("Centrality Class %d", icent), "l");
		leg->Draw();
		gPad->Update();
	}
	return hBAR;

	/*
	  TH2F *tmp = (TH2F *)fContainer->FindObject("hNtlsBAR");
	  TH1D *hBAR = tmp->ProjectionY();
	  TH1F *hSTA = (TH1F *)fContainer->FindObject("hNtlsSTA");
	  TH1F *hCON = (TH1F *)fContainer->FindObject("hNtls");
	  TLegend *leg = new TLegend(0.13, 0.75, 0.39, 0.89);
	  leg->SetBorderSize(1);
	  leg->SetFillColor(0);

	  Float_t scale = hCON->Integral();
	  if(scale) hCON->Scale(100./scale);
	  hCON->SetFillColor(kRed);hCON->SetLineColor(kRed);
	  hCON->SetBarWidth(0.2);
	  hCON->SetBarOffset(0.6);
	  hCON->SetStats(kFALSE);
	  hCON->GetYaxis()->SetRangeUser(0.,40.);
	  hCON->GetYaxis()->SetTitleOffset(1.2);
	  hCON->Draw("bar1"); leg->AddEntry(hCON, "Total", "f");
	  hCON->SetMaximum(55.);

	  if(scale) hBAR->Scale(100./scale);
	  hBAR->SetFillColor(kGreen);hBAR->SetLineColor(kGreen);
	  hBAR->SetBarWidth(0.2);
	  hBAR->SetBarOffset(0.2);
	  hBAR->SetTitle("");
	  hBAR->SetStats(kFALSE);
	  hBAR->GetYaxis()->SetRangeUser(0.,40.);
	  hBAR->GetYaxis()->SetTitleOffset(1.2);
	  hBAR->Draw("bar1same"); leg->AddEntry(hBAR, "Barrel", "f");

	  if(scale) hSTA->Scale(100./scale);
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
	 */
}

//________________________________________________________
void AliTRDcheckDET::MakePlotNclustersTrack(){
	//
	// Plot number of clusters
	// Put histos from all centrality classes into one pad
	//
	TH2 *hClusters = dynamic_cast<TH2 *>(fContainer->FindObject("hNcls"));
	if(!hClusters){
		AliError("Cluster histogram not found in the output");
	}
	TH1 *hAllCentrality = hClusters->ProjectionX("hNcls_all");
	hAllCentrality->SetTitle("Number of clusters/track");
	hAllCentrality->Scale(1./hAllCentrality->Integral());
	hAllCentrality->SetLineColor(kBlack);
	hAllCentrality->SetLineWidth(2);
	hAllCentrality->GetYaxis()->SetRangeUser(0, 0.02);
	hAllCentrality->Draw();
	PutTrendValue("NClustersTrack", hAllCentrality->GetMean());
	PutTrendValue("NClustersTrackRMS", hAllCentrality->GetRMS());

	// Now look at single centrality classes
	TH1 *hProjCentral[AliTRDeventInfo::kCentralityClasses];
	Int_t nHistsCentrality = 0;
	for(Int_t icent = 0; icent < AliTRDeventInfo::kCentralityClasses; icent++){
		hProjCentral[icent] = hClusters->ProjectionX(Form("hNcls_%d", icent), icent+1, icent+1);
		if(!hProjCentral[icent]->GetEntries()){
			delete hProjCentral[icent];
			hProjCentral[icent] = NULL;
			continue;
		}
		hProjCentral[icent]->Scale(1./hProjCentral[icent]->Integral());
		hProjCentral[icent]->SetTitle("Number of clusters/track");
		hProjCentral[icent]->SetLineColor(fkColorsCentrality[icent]);
		hProjCentral[icent]->GetYaxis()->SetRangeUser(0, 0.03);
		hProjCentral[icent]->Draw("same");
		nHistsCentrality++;
	}
	AliInfo(Form("%d centrality classes found", nHistsCentrality));
	if(nHistsCentrality){
		// Draw nice legend
		TLegend *leg = new TLegend(0.5, 0.6, 0.89, 0.89);
		leg->SetBorderSize(0);
		leg->SetFillStyle(0);
		for(Int_t icent = 0; icent < AliTRDeventInfo::kCentralityClasses; icent++){
			if(hProjCentral[icent]) leg->AddEntry(hProjCentral[icent], Form("Centrality Class %d", icent), "l");
		}
		leg->Draw();
		gPad->Update();
	}
}

//________________________________________________________
void AliTRDcheckDET::MakePlotNclustersTracklet(){
	//
	// Plot number of clusters for different centrality classes
	//
	TH2 *hClusters = dynamic_cast<TH2*>(fContainer->FindObject("hNclTls"));
	if(!hClusters){
		AliError("Histogram for the number of clusters per tracklet not available");
		return;
	}
	TH1 *hAllCentrality = hClusters->ProjectionX("hNclsTls_all");
	hAllCentrality->SetTitle("Number of clusters/track");
	hAllCentrality->Scale(1./hAllCentrality->Integral());
	hAllCentrality->SetLineColor(kBlack);
	hAllCentrality->SetLineWidth(2);
	hAllCentrality->GetYaxis()->SetRangeUser(0, 0.3);
	hAllCentrality->Draw("pc");
	PutTrendValue("NClustersTracklet", hAllCentrality->GetMean());
	PutTrendValue("NClustersTrackletRMS", hAllCentrality->GetRMS());

	// Now look at single centrality classes
	TH1 *hProjCentral[AliTRDeventInfo::kCentralityClasses];
	Int_t nHistsCentrality = 0;
	for(Int_t icent = 0; icent < AliTRDeventInfo::kCentralityClasses; icent++){
		hProjCentral[icent] = hClusters->ProjectionX(Form("hNclsTls_%d", icent), icent+1, icent+1);
		if(!hProjCentral[icent]->GetEntries()){
			delete hProjCentral[icent];
			hProjCentral[icent] = NULL;
			continue;
		}
		hProjCentral[icent]->Scale(1./hProjCentral[icent]->Integral());
		hProjCentral[icent]->SetTitle("Number of clusters/track");
		hProjCentral[icent]->SetLineColor(fkColorsCentrality[icent]);
		hProjCentral[icent]->GetYaxis()->SetRangeUser(0, 0.3);
		hProjCentral[icent]->Draw("pcsame");
		nHistsCentrality++;
	}
	AliInfo(Form("%d centrality classes found", nHistsCentrality));
	if(nHistsCentrality){
		// Draw nice legend
		TLegend *leg = new TLegend(0.5, 0.6, 0.89, 0.89);
		leg->SetBorderSize(0);
		leg->SetFillStyle(0);
		for(Int_t icent = 0; icent < AliTRDeventInfo::kCentralityClasses; icent++){
			if(hProjCentral[icent]) leg->AddEntry(hProjCentral[icent], Form("Centrality Class %d", icent), "l");
		}
		leg->Draw();
		gPad->Update();
	}
}

//________________________________________________________
void AliTRDcheckDET::MakePlotTrackletCharge(){
	//
	// Draw tracklet charge for different centrality classes
	//
	TH2 *hQt = dynamic_cast<TH2*>(fContainer->FindObject("hQtrklt"));
	if(!hQt){
		AliError("Histogram for tracklet charges not found");
		return;
	}
	// First Project all charhes
	TH1 *hQtAll = hQt->ProjectionX("hQtAll");
	hQtAll->SetTitle("Tracklet Charge");
	Double_t scalefactor = 0.7 * hQtAll->Integral() / hQtAll->GetMaximum();
	hQtAll->Scale(scalefactor/hQtAll->Integral());
	hQtAll->GetYaxis()->SetRangeUser(0., 1.);
	hQtAll->SetLineColor(kBlack);
	hQtAll->SetLineWidth(2);
	hQtAll->Draw("c");
	PutTrendValue("ChargeTracklet", hQtAll->GetMaximumBin());
	PutTrendValue("ChargeTrackletRMS", hQtAll->GetRMS());

	TH1 *hQtCent[AliTRDeventInfo::kCentralityClasses];
	Int_t nHistsCentrality = 0;
	for(Int_t icent = 0; icent < AliTRDeventInfo::kCentralityClasses; icent++){
		hQtCent[icent] = hQt->ProjectionX(Form("hQt_%d", icent), icent+1, icent+1);
		if(!hQtCent[icent]->GetEntries()){
			delete hQtCent[icent];
			continue;
		}
		hQtCent[icent]->SetTitle("Tracklet Charge");
		hQtCent[icent]->Scale(scalefactor/hQtCent[icent]->Integral());
		hQtCent[icent]->GetYaxis()->SetRangeUser(0., 1.);
		hQtCent[icent]->SetLineColor(fkColorsCentrality[icent]);
		hQtCent[icent]->Draw("csame");
		nHistsCentrality++;
	}
	AliInfo(Form("%d centrality classes found", nHistsCentrality));
	if(nHistsCentrality){
		TLegend *leg = new TLegend(0.5, 0.6, 0.89, 0.89);
		leg->SetBorderSize(0);
		leg->SetFillStyle(0);
		for(Int_t icent = 0; icent < AliTRDeventInfo::kCentralityClasses; icent++){
			if(hQtCent[icent]) leg->AddEntry(hQtCent[icent], Form("Centrality Class %d", icent), "l");
		}
		leg->Draw();
		gPad->Update();
	}
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

  TH3 *hBar3D = dynamic_cast<TH3F *>(fContainer->FindObject("hNtlsBAR"));
  if(!hBar3D){
	  AliError("Histogram for the number of tracklets vs p not available");
	  return;
  }
  TH2 *hBar = (TH2 *)hBar3D->Project3D("yx");
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

	TObjArray *arr = (TObjArray*)fContainer->FindObject("<PH>");
	//TH3F *hPhx = dynamic_cast<TH3F *>(arr->At(0)),
	TH3F *hPht = dynamic_cast<TH3F *>(arr->At(1));
	if(!hPht) return kFALSE;
	// Project centrality of the 2 histograms
	//TH2 *hProjCentX = dynamic_cast<TH2 *>(hPhx->Project3D("yx")),
	TH2 *hProjCentT = dynamic_cast<TH2 *>(hPht->Project3D("yx"));
  if(!hProjCentT) return kFALSE;
	//hProjCentX->SetName("hProjCentX");
	hProjCentT->SetName("hProjCentT");
	// Draw 2D histogram versus time on pad 2
	output->cd(2);
	hProjCentT->SetTitle("2D Pulse Height vs. Time");
	hProjCentT->GetYaxis()->SetTitleOffset(1.8);
  hProjCentT->GetYaxis()->SetRangeUser(0., 600.);
	hProjCentT->SetStats(kFALSE);
	hProjCentT->Draw("colz");

  // Fit each slice of the pulse height spectrum with a Landau function
  TGraphErrors *landaufit = new TGraphErrors;
  landaufit->SetTitle();
  landaufit->SetMarkerColor(kWhite);
  landaufit->SetLineColor(kWhite);
  landaufit->SetMarkerStyle(25);
  TH1 *projection;
  Int_t ntime = 0;
  for(Int_t it = 1; it <= hProjCentT->GetXaxis()->GetLast(); it++){
    projection = hProjCentT->ProjectionY("proj", it, it);
    if(!projection->GetEntries()){
      delete projection;
      continue;  
    }
    TF1 fitfun("fitfun", "landau", 0, 1000);
    fitfun.SetParLimits(0, 1e-1, 1e9); fitfun.SetParameter(0, 1000);
    fitfun.SetParLimits(1, 10, 200); fitfun.SetParameter(1, 50);
    fitfun.SetParLimits(2, 1e-1, 200); fitfun.SetParameter(2, 10);
    projection->Fit(&fitfun, "QN", "", 0, 1000); 
    landaufit->SetPoint(ntime, hProjCentT->GetXaxis()->GetBinCenter(it), fitfun.GetParameter(1));
    landaufit->SetPointError(ntime, hProjCentT->GetXaxis()->GetBinWidth(it)/2., fitfun.GetParError(1));
    ntime++;
    delete projection;
  }
  landaufit->Draw("lpesame");

	// Fill 1D PHS as function of time resp. radius (same binning of the 2 histograms)
	//TH1 *hPhsX = new TH1F("hPhsX", "Average PH vs X", hProjCentT->GetXaxis()->GetNbins(), hProjCentT->GetXaxis()->GetXbins()->GetArray()),
  Int_t nbinsT = hProjCentT->GetXaxis()->GetNbins();
	TH1	*hPhsT = new TH1F("hPhsT", "Average Ph vs Time; Time (100 ns); Average Pulse Height (a.u.)", nbinsT, -0.5, nbinsT - 0.5),
		*htmp = NULL;
	/*for(Int_t irad = 1; irad <= hProjCentX->GetXaxis()->GetNbins(); irad++){
		htmp = hProjCentX->ProjectionY("tmp", irad, irad);
		hPhsX->SetBinContent(irad, htmp->GetMean());
		hPhsX->SetBinError(irad, htmp->GetMeanError());
		delete htmp;
	}*/
	for(Int_t it = 1; it <= hProjCentT->GetXaxis()->GetNbins(); it++){
		htmp = hProjCentT->ProjectionY("tmp", it, it);
		hPhsT->SetBinContent(it, htmp->GetMean());
		hPhsT->SetBinError(it, htmp->GetMeanError());
		delete htmp;
	}
	output->cd(1);
	// Draw 1D histograms
	if(hPhsT->GetEntries()){
		hPhsT->SetMarkerStyle(24);
		hPhsT->SetMarkerColor(kBlack);
		hPhsT->SetLineColor(kBlack);
		hPhsT->GetYaxis()->SetTitleOffset(1.5);
		hPhsT->Draw("e1");
		// Now fit the PHS with respect to time to get plateau and slope
		// Trending for the pulse height: plateau value, slope and timebin of the maximum
		TLinearFitter fit(1,"pol1");
		Double_t time = 0.;
		for(Int_t itime = 10; itime <= 20; itime++){
			time = Double_t(itime);
			Double_t err(hPhsT->GetBinError(itime + 1));
			if(err>1.e-10) fit.AddPoint(&time, hPhsT->GetBinContent(itime + 1), err);
		}
		if(!fit.Eval()){
			Double_t plateau = fit.GetParameter(0) + 12 * fit.GetParameter(1);
			Double_t slope = fit.GetParameter(1);
			PutTrendValue("PHplateau", plateau);
			PutTrendValue("PHslope", slope);
			PutTrendValue("PHamplificationPeak", static_cast<Double_t>(hPhsT->GetMaximumBin()-1));
			AliDebug(1, Form("plateau %f, slope %f, MaxTime %f", plateau, slope, static_cast<Double_t>(hPhsT->GetMaximumBin()-1)));
		}
	}
	/*if(hPhsX->GetEntries()){
		hPhsX->SetMarkerStyle(22);
		hPhsX->SetMarkerColor(kBlue);
		hPhsX->SetLineColor(kBlue);
		hPhsX->GetYaxis()->SetTitleOffset(1.5);
		hPhsX->Draw("e1same");
		// create axis according to the histogram dimensions of the original second histogram
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
		gPad->Update();
	}*/

	// Centrality-dependent Pulse-Height Spectrum
	TH1 *hPhtCent[AliTRDeventInfo::kCentralityClasses];
	TH2 *hPtmp;
	Int_t nHistsCentrality = 0;
	for(Int_t icent = 0; icent < AliTRDeventInfo::kCentralityClasses; icent++){
		hPht->GetZaxis()->SetRange(icent+2,icent+2);
		hPtmp = dynamic_cast<TH2*>(hPht->Project3D("yx"));
		if(!(hPtmp && hPtmp->GetEntries())){
			hPhtCent[icent] = NULL;
			continue;
		}
		hPhtCent[icent] = new TH1F(Form("hPhtCent_%d", icent), "Average Ph vs Time", hPtmp->GetNbinsX(), hPtmp->GetXaxis()->GetXbins()->GetArray());
		for(Int_t it = 1; it <= hPtmp->GetNbinsX(); it++){
			htmp = hPtmp->ProjectionY("htmp", it, it);
			hPhtCent[icent]->SetBinContent(it, htmp->GetMean());
			hPhtCent[icent]->SetBinError(it, htmp->GetMeanError());
			delete htmp;
		}
		delete hPtmp;
		hPhtCent[icent]->SetMarkerStyle(24);
		hPhtCent[icent]->SetMarkerColor(fkColorsCentrality[icent]);
		hPhtCent[icent]->SetLineColor(fkColorsCentrality[icent]);
		hPhtCent[icent]->Draw("e1same");
		nHistsCentrality++;
	}
	hPht->GetZaxis()->SetRange(0, hPht->GetNbinsZ());
	if(nHistsCentrality){
		TLegend *leg = new TLegend(0.5, 0.6, 0.89, 0.89);
		leg->SetBorderSize(0);
		leg->SetFillStyle(0);
		for(Int_t icent = 0; icent < AliTRDeventInfo::kCentralityClasses; icent++){
			if(hPhtCent[icent]) leg->AddEntry(hPhtCent[icent], Form("Centrality Class %d", icent), "p");
		}
		leg->Draw();
		gPad->Update();
	}

	// delete 2D Projection of the PHS vs x since it is only used to calculate the 1D projection
	//delete hProjCentX;

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
    if(!(hlayer = dynamic_cast<TProfile2D *>(histos->At(ily)))) continue;
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
