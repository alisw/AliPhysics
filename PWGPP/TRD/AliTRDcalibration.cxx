
/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/

/////////////////////////////////////////////////////////////////////////////////
//                                                                             
// AliTRDcalibration                                                            
//                                                                             
// Task to run the calibration offline.
// Author:
//   R. Bailhache (rbailhache@ikf.uni-frankfurt.de, R.Bailhache@gsi.de)
//           
//////////////////////////////////////////////////////////////////////////////////


#include "Riostream.h"
#include "TChain.h"
#include "TTree.h"
#include "TProfile2D.h"
#include "TH2I.h"
#include "TH1F.h"
#include "TList.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TObject.h"
#include "TFile.h"
#include "TObjArray.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TGraphErrors.h"

#include "AliTRDrecoTask.h"
#include "AliAnalysisManager.h"

#include "AliESDInputHandler.h"
#include "AliTRDtrackV1.h"
#include "AliTRDseedV1.h"
#include "AliTRDcluster.h"
#include "info/AliTRDtrackInfo.h"
#include "AliTRDcalibDB.h"

#include "AliTRDCalibraFillHisto.h"
#include "AliTRDCalibraFit.h"
#include "AliTRDCalibraVdriftLinearFit.h"
#include "AliTRDCalibraMode.h"
#include "AliTRDCalibraVector.h"
#include "./Cal/AliTRDCalPad.h"
#include "./Cal/AliTRDCalDet.h"

#include "AliLog.h"

#include "AliTRDcalibration.h"


ClassImp(AliTRDcalibration)

//________________________________________________________________________
AliTRDcalibration::AliTRDcalibration() 
  :AliTRDrecoTask()
  ,fTrackInfo(0)
  ,ftrdTrack(0)
  ,fcl(0)
  ,fTRDCalibraFillHisto(0)
  ,fNbTRDTrack(0)
  ,fNbTRDTrackOffline(0)
  ,fNbTRDTrackStandalone(0)
  ,fNbTRDTracklet(0)
  ,fNbTRDTrackletOffline(0)
  ,fNbTRDTrackletStandalone(0)
  ,fNbTimeBin(0x0)
  ,fNbTimeBinOffline(0x0)
  ,fNbTimeBinStandalone(0x0)
  ,fNbClusters(0)
  ,fNbClustersOffline(0)
  ,fNbClustersStandalone(0)
  ,fPHSM(0)
  ,fCHSM(0)
  ,fPHSum(0)
  ,fCHSum(0)
  ,fDetSum(0)
  ,fDetSumVector(0)
  ,fHisto2d(kTRUE)
  ,fVector2d(kFALSE)
  ,fVdriftLinear(kTRUE)
  ,flow(0)
  ,fhigh(30)
  ,fNbTimeBins(0)
  ,ffillZero(kFALSE)
  ,fnormalizeNbOfCluster(kFALSE)
  ,fmaxCluster(0)
  ,fOfflineTracks(kFALSE)
  ,fStandaloneTracks(kFALSE)
  ,fCompressPerDetector(kFALSE)
  ,fGraph(0x0)
  ,fArrayCalib(0x0)
  ,fPostProcess(kFALSE)
{
  // Constructor
  
  fNRefFigures = 17;

  for(Int_t k = 0; k < 3; k++)
    {
      fNz[k]=0;
      fNrphi[k]=0;
    }

}  

AliTRDcalibration::AliTRDcalibration(char* name) 
  :AliTRDrecoTask(name, "Calibration on tracks")
  ,fTrackInfo(0)
  ,ftrdTrack(0)
  ,fcl(0)
  ,fTRDCalibraFillHisto(0)
  ,fNbTRDTrack(0)
  ,fNbTRDTrackOffline(0)
  ,fNbTRDTrackStandalone(0)
  ,fNbTRDTracklet(0)
  ,fNbTRDTrackletOffline(0)
  ,fNbTRDTrackletStandalone(0)
  ,fNbTimeBin(0x0)
  ,fNbTimeBinOffline(0x0)
  ,fNbTimeBinStandalone(0x0)
  ,fNbClusters(0)
  ,fNbClustersOffline(0)
  ,fNbClustersStandalone(0)
  ,fPHSM(0)
  ,fCHSM(0)
  ,fPHSum(0)
  ,fCHSum(0)
  ,fDetSum(0)
  ,fDetSumVector(0)
  ,fHisto2d(kTRUE)
  ,fVector2d(kFALSE)
  ,fVdriftLinear(kTRUE)
  ,flow(0)
  ,fhigh(30)
  ,fNbTimeBins(0)
  ,ffillZero(kFALSE)
  ,fnormalizeNbOfCluster(kFALSE)
  ,fmaxCluster(0)
  ,fOfflineTracks(kFALSE)
  ,fStandaloneTracks(kFALSE)
  ,fCompressPerDetector(kFALSE)
  ,fGraph(0x0)
  ,fArrayCalib(0x0)
  ,fPostProcess(kFALSE)
{
  // Constructor
  
  fNRefFigures = 17;

  for(Int_t k = 0; k < 3; k++)
    {
      fNz[k]=0;
      fNrphi[k]=0;
    }

}  

//________________________________________________________________________
AliTRDcalibration::~AliTRDcalibration() 
{
  // Default destructor
  if (AliAnalysisManager::GetAnalysisManager()->IsProofMode()) return;
  //
  if(fNbTRDTrack) delete fNbTRDTrack;
  if(fNbTRDTrackOffline) delete fNbTRDTrackOffline;
  if(fNbTRDTrackStandalone) delete fNbTRDTrackStandalone;
  if(fNbTRDTracklet) delete fNbTRDTracklet;
  if(fNbTRDTrackletOffline) delete fNbTRDTrackletOffline;
  if(fNbTRDTrackletStandalone) delete fNbTRDTrackletStandalone;
  if(fNbTimeBin) delete fNbTimeBin;
  if(fNbTimeBinOffline) delete fNbTimeBinOffline;
  if(fNbTimeBinStandalone) delete fNbTimeBinStandalone;
  if(fNbClusters) delete fNbClusters;
  if(fNbClustersOffline) delete fNbClustersOffline;
  if(fNbClustersStandalone) delete fNbClustersStandalone;
  if(fPHSM) delete fPHSM;
  if(fCHSM) delete fCHSM;
  if(fPHSum) delete fPHSum;
  if(fCHSum) delete fCHSum;
  if(fDetSum) delete fDetSum;
  if(fDetSumVector) delete fDetSumVector;
  if(fGraph){fGraph->Delete(); delete fGraph;}
  if(fArrayCalib){fArrayCalib->Delete(); delete fArrayCalib;}
   
}
//________________________________________________________________________
void AliTRDcalibration::UserCreateOutputObjects() 
{
  // Create output objects

  // Number of time bins
  if(fNbTimeBins==0) {
    AliTRDcalibDB *cal = AliTRDcalibDB::Instance();
    fNbTimeBins = cal->GetNumberOfTimeBinsDCS();
    if(fNbTimeBins <= 0){ 
      AliWarning(Form("No of TimeBins from DB [%d] use default [30]", fNbTimeBins));
      fNbTimeBins = 30;
    }
  }
  
  // instance calibration: what to calibrate
  fTRDCalibraFillHisto = AliTRDCalibraFillHisto::Instance();
  fTRDCalibraFillHisto->SetHisto2d(fHisto2d); // choose to use histograms
  fTRDCalibraFillHisto->SetVector2d(fVector2d); // choose to use vectors
  fTRDCalibraFillHisto->SetCH2dOn();  // choose to calibrate the gain
  fTRDCalibraFillHisto->SetPH2dOn();  // choose to calibrate the drift velocity
  fTRDCalibraFillHisto->SetPRF2dOn(); // choose to look at the PRF
  fTRDCalibraFillHisto->SetLinearFitterOn(fVdriftLinear); // Other possibility vdrift VDRIFT
  fTRDCalibraFillHisto->SetLinearFitterDebugOn(fVdriftLinear); // Other possibility vdrift

  // segmentation (should be per default the max and add at the end)
  for(Int_t k = 0; k < 3; k++){
    if(((fNz[k] != 10) && (fNrphi[k] != 10)) && ((fNz[k] != 100) && (fNrphi[k] != 100))) {
      fTRDCalibraFillHisto->SetNz(k,fNz[k]);                                    // Mode calibration
      fTRDCalibraFillHisto->SetNrphi(k,fNrphi[k]);                             // Mode calibration
    }
    else {
      if((fNz[k] == 100) && (fNrphi[k] == 100))  {
	if(fVector2d) AliInfo("The mode all together is not supported by the vector method");
	fTRDCalibraFillHisto->SetAllTogether(k);
      }
      if((fNz[k] == 10) && (fNrphi[k] == 10))  {
	if(fVector2d) AliInfo("The mode per supermodule is not supported by the vector method");
	fTRDCalibraFillHisto->SetPerSuperModule(k);
      }
    }
  }

  // Debug level
  fTRDCalibraFillHisto->SetDebugLevel(DebugLevel()); //debug stuff

  // Init the stuff
  fTRDCalibraFillHisto->Init2Dhistos(fNbTimeBins); // initialise the histos

  // cuts
  fTRDCalibraFillHisto->SetNumberClusters(flow); // At least flow clusters
  fTRDCalibraFillHisto->SetNumberClustersf(fhigh); // The more fhigh clusters
  fTRDCalibraFillHisto->SetFillWithZero(ffillZero); // Fill zeros
  fTRDCalibraFillHisto->SetNormalizeNbOfCluster(fnormalizeNbOfCluster); // For iterations

  // Add them to the container
  fContainer = new TObjArray();
  if(fHisto2d) {
    fContainer->Add(fTRDCalibraFillHisto->GetCH2d()); //TH2I
    fContainer->Add(fTRDCalibraFillHisto->GetPH2d()); //TProfile2D
    fContainer->Add(fTRDCalibraFillHisto->GetPRF2d()); //TProfile2D
  }
  if(fVdriftLinear) fContainer->Add(fTRDCalibraFillHisto->GetVdriftLinearFit()); // Other drift velocity 
  if(fVector2d) fContainer->Add(fTRDCalibraFillHisto->GetCalibraVector()); //calibra vector
      
  if(DebugLevel()) {
    
    // Init the debug histos
    fNbTRDTrack = new TH1F("TRDTrack","TRDTrack",500,0,500);
    fNbTRDTrack->Sumw2();
    fNbTRDTrackOffline = new TH1F("TRDTrackOffline","TRDTrackOffline",500,0,500);
    fNbTRDTrackOffline->Sumw2();
    fNbTRDTrackStandalone = new TH1F("TRDTrackStandalone","TRDTrackStandalone",500,0,500);
    fNbTRDTrackStandalone->Sumw2();
    //
    fNbTRDTracklet = new TH1F("TRDTracklet","TRDTracklet",540,0.,540.);
    fNbTRDTracklet->Sumw2();
    fNbTRDTrackletOffline = new TH1F("TRDTrackletOffline","TRDTrackletOffline",540,0.,540.);
    fNbTRDTrackletOffline->Sumw2();
    fNbTRDTrackletStandalone = new TH1F("TRDTrackletStandalone","TRDTrackletStandalone",540,0.,540.);
    fNbTRDTrackletStandalone->Sumw2();
    //
    fNbTimeBin = new TH1F("TimeBin","TimeBin",35,0,35);
    fNbTimeBin->Sumw2();
    fNbTimeBinOffline = new TH1F("TimeBinOffline","TimeBinOffline",35,0,35);
    fNbTimeBinOffline->Sumw2();
    fNbTimeBinStandalone = new TH1F("TimeBinStandalone","TimeBinStandalone",35,0,35);
    fNbTimeBinStandalone->Sumw2();
    //
    fNbClusters = new TH1F("NbClusters","",35,0,35);
    fNbClusters->Sumw2();
    fNbClustersOffline = new TH1F("NbClustersOffline","",35,0,35);
    fNbClustersOffline->Sumw2();
    fNbClustersStandalone = new TH1F("NbClustersStandalone","",35,0,35);
    fNbClustersStandalone->Sumw2();
    //
    fPHSM = new TProfile2D("PH2dSM","Nz10Nrphi10"
			    ,fNbTimeBins,-0.05,(Double_t)((fNbTimeBins-0.5)/10.0)
			    ,18,0,18);
    fPHSM->SetYTitle("Det/pad groups");
    fPHSM->SetXTitle("time [#mus]");
    fPHSM->SetZTitle("<PH> [a.u.]");
    fPHSM->SetStats(0);
    //
    fCHSM = new TH2I("CH2dSM","Nz10Nrphi10",50,0,300,18,0,18);
    fCHSM->SetYTitle("Det/pad groups");
    fCHSM->SetXTitle("charge deposit [a.u]");
    fCHSM->SetZTitle("counts");
    fCHSM->SetStats(0);
    fCHSM->Sumw2();
    //
    fPHSum = new TProfile2D("PH2dSum","Nz100Nrphi100"
			    ,fNbTimeBins,-0.05,(Double_t)((fNbTimeBins-0.5)/10.0)
			    ,1,0,1);
    fPHSum->SetYTitle("Det/pad groups");
    fPHSum->SetXTitle("time [#mus]");
    fPHSum->SetZTitle("<PH> [a.u.]");
    fPHSum->SetStats(0);
    //
    fCHSum = new TH2I("CH2dSum","Nz100Nrphi100",50,0,300,1,0,1);
    fCHSum->SetYTitle("Det/pad groups");
    fCHSum->SetXTitle("charge deposit [a.u]");
    fCHSum->SetZTitle("counts");
    fCHSum->SetStats(0);
    fCHSum->Sumw2();
    
    // Add them
    fContainer->Add(fNbTRDTrack);
    fContainer->Add(fNbTRDTrackOffline);
    fContainer->Add(fNbTRDTrackStandalone);
    fContainer->Add(fNbTRDTracklet);
    fContainer->Add(fNbTRDTrackletOffline);
    fContainer->Add(fNbTRDTrackletStandalone);
    fContainer->Add(fNbTimeBin);
    fContainer->Add(fNbTimeBinOffline);
    fContainer->Add(fNbTimeBinStandalone);
    fContainer->Add(fNbClusters);
    fContainer->Add(fNbClustersOffline);
    fContainer->Add(fNbClustersStandalone);
    fContainer->Add(fPHSM);
    fContainer->Add(fCHSM);
    fContainer->Add(fPHSum);
    fContainer->Add(fCHSum);

  }
  // Post output data
  PostData(1, fContainer);
}

//________________________________________________________________________
void AliTRDcalibration::UserExec(Option_t *) 
{
  //
  // Execute function where the reference data are filled
  //

  if(!fTracks) return;
  
  // In total
  Int_t nbTrdTracks = 0;
  // standalone
  Int_t nbTrdTracksStandalone = 0;
  // offline
  Int_t nbTrdTracksOffline = 0;
  

  //
  // Loop on track in the event
  //
  //printf("Total of %d\n",fTracks->GetEntriesFast());
  for(Int_t itrk=0; itrk < fTracks->GetEntriesFast(); itrk++){
    
    //printf("itrk %d\n",itrk);

    fTrackInfo = (AliTRDtrackInfo*)fTracks->UncheckedAt(itrk);
    ftrdTrack = fTrackInfo->GetTrack();
    if(!ftrdTrack) continue;

    nbTrdTracks++;
  
    fTRDCalibraFillHisto->UpdateHistogramsV1(ftrdTrack);

    if(DebugLevel()) {
      
      Bool_t standalonetracklet = kFALSE;  
      const AliTRDseedV1 *tracklet = 0x0;
      //
      // Loop on tracklet in the event
      //
      for(Int_t itr = 0; itr < 6; itr++){
	//printf("itr %d\n",itr);
 	if(!(tracklet = ftrdTrack->GetTracklet(itr))) continue;
	if(!tracklet->IsOK()) continue;
	// standalone
	if(tracklet->IsStandAlone()) standalonetracklet = kTRUE;
	Int_t nbclusters = 0;
	// For PH
	Double_t phtb[AliTRDseedV1::kNtb];
	memset(phtb, 0, AliTRDseedV1::kNtb*sizeof(Double_t));
	// For CH
	Double_t sum = 0.0;
	// normalisation
	Float_t normalisation = 6.67;
	Int_t detector = 0;
	Int_t sector = 0;
	Int_t crossrow = 0;
	// Check no shared clusters
	for(int icc=AliTRDseedV1::kNtb; icc<AliTRDseedV1::kNclusters; icc++){
	  if((fcl = tracklet->GetClusters(icc)))  crossrow = 1;
	}
	// Loop on clusters
	for(int ic=0; ic<AliTRDseedV1::kNtb; ic++){
	  if(!(fcl = tracklet->GetClusters(ic))) continue;
	  nbclusters++;
	  Int_t time = fcl->GetPadTime();
	  Float_t ch =  tracklet->GetdQdl(ic);
	  Float_t qcl = TMath::Abs(fcl->GetQ());
	  detector = fcl->GetDetector();	  
	  if((time>-1) && (time<fNbTimeBins)) phtb[time]=qcl;
	  sum += ch/normalisation;
	  fNbTimeBin->Fill(time);
	  if(tracklet->IsStandAlone()) fNbTimeBinStandalone->Fill(time);
	  else fNbTimeBinOffline->Fill(time);
	}
	sector = AliTRDgeometry::GetSector(detector);

	fNbTRDTracklet->Fill(detector);
	if(tracklet->IsStandAlone()) fNbTRDTrackletStandalone->Fill(detector);
	else fNbTRDTrackletOffline->Fill(detector);
	
	fNbClusters->Fill(nbclusters);
	if(tracklet->IsStandAlone())  fNbClustersStandalone->Fill(nbclusters);
	else  fNbClustersOffline->Fill(nbclusters);
	
	if((nbclusters > flow) && (nbclusters < fhigh)){
	  fCHSM->Fill(sum/20.0,sector+0.5);
	  fCHSum->Fill(sum/20.0,0.5);
	  for(int ic=0; ic<fNbTimeBins; ic++){
	    if(ffillZero) {
	      fPHSum->Fill((Double_t)(ic/10.0),0.5,(Double_t)phtb[ic]);
	      fPHSM->Fill((Double_t)(ic/10.0),sector+0.5,(Double_t)phtb[ic]);
	    }
	    else {
	      if(phtb[ic] > 0.0) {
		fPHSum->Fill((Double_t)(ic/10.0),0.0,(Double_t)phtb[ic]);
		fPHSM->Fill((Double_t)(ic/10.0),sector+0.5,(Double_t)phtb[ic]);
	      }
	    }
	  }
	}
      }
    
    if(standalonetracklet) nbTrdTracksStandalone++;
    else nbTrdTracksOffline++;
    
    }
    
  }
  
  if(DebugLevel()) {
    
    //Fill Histos
    fNbTRDTrack->Fill(nbTrdTracks);
    fNbTRDTrackStandalone->Fill(nbTrdTracksStandalone);
    fNbTRDTrackOffline->Fill(nbTrdTracksOffline);
    
  }
}  
    
//________________________________________________________________________
void AliTRDcalibration::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query

  //printf("terminate\n");

  if(fTRDCalibraFillHisto) fTRDCalibraFillHisto->DestroyDebugStreamer();
 
}
//________________________________________________________
Bool_t AliTRDcalibration::GetRefFigure(Int_t ifig)
{
  //
  // Draw filled histos
  //
  
  gStyle->SetPalette(1);
  gStyle->SetOptStat(1111);
  gStyle->SetPadBorderMode(0);
  gStyle->SetCanvasColor(10);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadRightMargin(0.13);

  if(!fContainer) return kFALSE;
  
  switch(ifig){
  case kNbTrack:{
    TCanvas *c0 = new TCanvas("c0","c0",10,10,510,510);
    TLegend *legNbTrack = new TLegend(.7, .7, .98, .98);
    legNbTrack->SetBorderSize(1);
    TH1F *h  = 0x0;
    TH1F *ho = 0x0;
    TH1F *hs = 0x0;
    if(!(h = (TH1F *)fContainer->FindObject("TRDTrack"))) break;
    if(!(ho = (TH1F *)fContainer->FindObject("TRDTrackOffline"))) break;
    if(!(hs = (TH1F *)fContainer->FindObject("TRDTrackStandalone"))) break;
    c0->cd();
    //gPad->SetLogy();
    gPad->SetGridy();
    gPad->SetGridx();
    h->Draw();
    ho->Draw("same");
    hs->Draw("same");
    legNbTrack->AddEntry(h, "all", "p");
    legNbTrack->AddEntry(ho, "offline", "p");
    legNbTrack->AddEntry(hs, "standalone", "p");
    legNbTrack->Draw("same");
    return kTRUE;
  }
  case kNbTracklet:{
    TLegend *legNbTracklet = new TLegend(.7, .7, .98, .98);
    legNbTracklet->SetBorderSize(1);
    TH1F *h = 0x0;
    TH1F *ho = 0x0;
    TH1F *hs = 0x0;
    if(!(h = (TH1F *)fContainer->FindObject("TRDTracklet"))) break;
    if(!(ho = (TH1F *)fContainer->FindObject("TRDTrackletOffline"))) break;
    if(!(hs = (TH1F *)fContainer->FindObject("TRDTrackletStandalone"))) break;
    h->Draw();
    ho->Draw("same");
    hs->Draw("same");
    legNbTracklet->AddEntry(h, "all", "p");
    legNbTracklet->AddEntry(ho, "offline", "p");
    legNbTracklet->AddEntry(hs, "standalone", "p");
    legNbTracklet->Draw("same");
    gPad->SetLogy();
    gPad->SetGridy();
    gPad->SetGridx();
    return kTRUE;
  }
  case kNbTimeBin:{
    TLegend *legNbTimeBin = new TLegend(.7, .7, .98, .98);
    legNbTimeBin->SetBorderSize(1);
    TH1F *h = 0x0;
    TH1F *ho = 0x0;
    TH1F *hs = 0x0;
    if(!(h = (TH1F *)fContainer->FindObject("TimeBin"))) break;
    if(!(ho = (TH1F *)fContainer->FindObject("TimeBinOffline"))) break;
    if(!(hs = (TH1F *)fContainer->FindObject("TimeBinStandalone"))) break;
    h->Draw();
    ho->Draw("same");
    hs->Draw("same");
    legNbTimeBin->AddEntry(h, "all", "p");
    legNbTimeBin->AddEntry(ho, "offline", "p");
    legNbTimeBin->AddEntry(hs, "standalone", "p");
    legNbTimeBin->Draw("same");
    //gPad->SetLogy();
    gPad->SetGridy();
    gPad->SetGridx();
    return kTRUE;
  }
  case kNbClusters:{
    TLegend *legNbClusters = new TLegend(.7, .7, .98, .98);
    legNbClusters->SetBorderSize(1);
    TH1F *h = 0x0;
    TH1F *ho = 0x0;
    TH1F *hs = 0x0;
    if(!(h = (TH1F *)fContainer->FindObject("NbClusters"))) break;
    if(!(ho = (TH1F *)fContainer->FindObject("NbClustersOffline"))) break;
    if(!(hs = (TH1F *)fContainer->FindObject("NbClustersStandalone"))) break;
    h->Draw();
    ho->Draw("same");
    hs->Draw("same");
    legNbClusters->AddEntry(h, "all", "p");
    legNbClusters->AddEntry(ho, "offline", "p");
    legNbClusters->AddEntry(hs, "standalone", "p");
    legNbClusters->Draw("same");
    gPad->SetLogy();
    gPad->SetGridy();
    gPad->SetGridx();
    return kTRUE;
  }
  case kPHSum:{
    TProfile2D *h = 0x0;
    if(!(h = (TProfile2D *)fContainer->FindObject("PH2dSum"))) break;
    TH1D *projh = h->ProjectionX("projh",1,1,"e");
    projh->Draw();
    gPad->SetGridy();
    gPad->SetGridx();
    return kTRUE;
  }
  case kCHSum:{
    TH2I *h = 0x0;
    if(!(h = (TH2I *)fContainer->FindObject("CH2dSum"))) break;
    TH1D *projh = h->ProjectionX("projhh",1,1,"e");
    projh->Draw();
    gPad->SetGridy();
    gPad->SetGridx();
    return kTRUE;
  }
  case kPH2D:{
    if(!fHisto2d) {
      AliInfo("Histo was not filled!");
      break;
    }
    TProfile2D *h = 0x0;
    if(!(h = (TProfile2D *)fContainer->FindObject("PH2d"))) break;
    h->Draw("lego");
    return kTRUE;
  }
  case kCH2D:{
    if(!fHisto2d) {
      AliInfo("Histo was not filled!");
      break;
    }
    TH2I *h = 0x0;
    if(!(h = (TH2I *)fContainer->FindObject("CH2d"))) break;
    h->Draw("lego");
    return kTRUE;
  }
  case kPRF2D:{
    if(!fHisto2d) {
      AliInfo("Histo was not filled!");
      break;
    }
    TProfile2D *h = 0x0;
    if(!(h = (TProfile2D *)fContainer->FindObject("PRF2d"))) break;
    h->Draw("lego");
    return kTRUE;
  }
  case kPH2DVector:{
    if(!fVector2d) {
      AliInfo("vector was not filled!");
      break;
    }
    AliTRDCalibraVector *v = 0x0;
    TGraphErrors *vdet = 0x0; 
    if(!(v = (AliTRDCalibraVector *)fContainer->FindObject("AliTRDCalibraVector"))) break;
    Int_t detectormax = -1;
    Int_t groupmax    = -1;
    if(!v->FindTheMaxEntries(1,detectormax,groupmax)) break;
    if(!(vdet = v->ConvertVectorPHTGraphErrors((Int_t)detectormax,groupmax,"plotPH2dVector"))) break;
    Int_t nbeentries = 0;
    TH1F *ko = v->CorrectTheError(vdet,nbeentries);
    ko->Draw();
    AliInfo(Form("There are %d entries in the detector %d and group %d",nbeentries,detectormax,groupmax));
    return kTRUE;
  }
case kCH2DVector:{
    if(!fVector2d) {
      AliInfo("vector was not filled!");
      break;
    }
    AliTRDCalibraVector *v = 0x0;
    TH1F *vdet = 0x0; 
    if(!(v = (AliTRDCalibraVector *)fContainer->FindObject("AliTRDCalibraVector"))) break;
    Int_t detectormax = -1;
    Int_t groupmax    = -1;
    if(!v->FindTheMaxEntries(0,detectormax,groupmax)) break;
    if(!(vdet = v->ConvertVectorCHHisto((Int_t)detectormax,groupmax,"plotCH2dVector"))) break;
    vdet->Draw();
    AliInfo(Form("The detectormax and groupmax are %d and %d",detectormax,groupmax));
    return kTRUE;
  }
  case kPRF2DVector:{
    if(!fVector2d) {
      AliInfo("vector was not filled!");
      break;
    }
    AliTRDCalibraVector *v = 0x0;
    TGraphErrors *vdet = 0x0; 
    if(!(v = (AliTRDCalibraVector *)fContainer->FindObject("AliTRDCalibraVector"))) break;
    Int_t detectormax  = -1;
    Int_t groupmax     = -1;
    Int_t nbeentries   = 0;
    if(!v->FindTheMaxEntries(2,detectormax,groupmax)) break;
    if(!(vdet = v->ConvertVectorPRFTGraphErrors((Int_t)detectormax,groupmax,"plotPRF2dVector"))) break;
    TH1F *ko = v->CorrectTheError(vdet,nbeentries);
    ko->Draw();
    AliInfo(Form("The detectormax and groupmax are %d and %d",detectormax,groupmax));
    return kTRUE;
  }
  case kLinearFitter:{
    if(!fVdriftLinear) {
      AliInfo("vdrift linear was not filled!");
      break;
    }
    AliTRDCalibraVdriftLinearFit *h = 0x0;
    TH2S *hdetector = 0x0; 
    if(!(h = (AliTRDCalibraVdriftLinearFit *)fContainer->FindObject("AliTRDCalibraVdriftLinearFit"))) break;
    Double_t entries[540];
    for(Int_t k = 0; k < 540; k++){
      entries[k] = 0.0;
      hdetector = 0x0;
      if(!(hdetector = (TH2S *)h->GetLinearFitterHisto(k,kFALSE))) continue;
      entries[k] = hdetector->GetEntries();
    }
    Double_t max = -10.0;
    Double_t detectormax = -1;
    for(Int_t k = 0; k < 540; k++){
      if(entries[k] > max) {
	max = entries[k];
	detectormax = k;
      }
    }
    hdetector = 0x0;
    if((TMath::Abs(max) <= 0.001) || (detectormax <0.0) || (detectormax >=540.0)) break;
    if(!(hdetector = (TH2S *)h->GetLinearFitterHisto((Int_t)detectormax,kFALSE))) break;
    AliInfo(Form("The detector with the maximum of entries is %f",detectormax));
    hdetector->Draw();
    return kTRUE;
  }
  case kGainFactor:{
    if(!fPostProcess){
      if(!PostProcess()) break;
    }
    TGraph *fgain = (TGraph *) fGraph->At(0);
    if(!fgain) break;
    fgain->Draw("ALP");
    return kTRUE;
  }
  case kVdriftT0Factor:{
    if(!fPostProcess){
      if(!PostProcess()) break;
    }
    TCanvas *c = new TCanvas("c","c",10,10,510,510);
    c->Divide(2,1);
    TGraph *fvd = (TGraph *) fGraph->At(1);
    if(fvd){
      c->cd(1);
      fvd->Draw("ALP");
    } 
    TGraph *ft0 = (TGraph *) fGraph->At(2);
    if(ft0){
      c->cd(2);
      ft0->Draw("ALP");
    } 
    return kTRUE;
  }
  case kVdriftLorentzAngleFactor:{
    if(!fVdriftLinear) {
      AliInfo("vdrift linear was not filled!");
      break;
    }
    if(!fPostProcess){
      if(!PostProcess()) break;
    }
    TCanvas *c = new TCanvas("c","c",10,10,510,510);
    c->Divide(2,1);
    TGraph *fvdl = (TGraph *) fGraph->At(3);
    if(fvdl){
      c->cd(1);
      fvdl->Draw("ALP");
    } 
    TGraph *flal = (TGraph *) fGraph->At(4);
    if(flal){
      c->cd(2);
      flal->Draw("ALP");
    } 
    return kTRUE;
  }
  case kPRFFactor:{
    if(!fPostProcess){
      if(!PostProcess()) break;
    }
    TGraph *fprf = (TGraph *) fGraph->At(5);
    if(!fprf) break;
    fprf->Draw("ALP");
    return kTRUE;
  }
  }
  
  return kFALSE;
  
}
//________________________________________________________________________
Bool_t AliTRDcalibration::PostProcess()
{
  // 
  // Fit the filled histos
  // Put the calibration object in fArrayCalib
  // 0 and 1 AliTRDCalDet and AliTRDCalPad gain
  // 2 and 3 AliTRDCalDet and AliTRDCalPad vdrift PH
  // 4 and 5 AliTRDCalDet and AliTRDCalPad timeoffset PH
  // 6 AliTRDCalPad PRF
  // 7 and 8 AliTRDCalDet and AliTRDCalPad vdrift second
  // 9 and 10 AliTRDCalDet and AliTRDCalPad lorentz angle second
  //

  if(!fArrayCalib){
    fArrayCalib = new TObjArray(11);
    fArrayCalib->SetOwner();
  }
  else {
    delete fArrayCalib;
    PostProcess();
  }
  
  if(!fGraph){
    fGraph = new TObjArray(6);
    fGraph->SetOwner();
  }
  else {
    delete fGraph;
    PostProcess();
  }

  Bool_t storage[3] = {kFALSE,kFALSE,kFALSE};

  // Objects for fitting
  AliTRDCalibraFit *calibra = AliTRDCalibraFit::Instance();
  calibra->SetDebugLevel(2); // 0 rien, 1 fitvoir, 2 debug files, 3 one detector  
  
  // Take the stuff
  if (!fContainer) {
    Printf("ERROR: list not available");
    return kFALSE;
  }

  if(fHisto2d && fVector2d) AliInfo("We will only look at histos. Set fHisto2d off if you don't want");
  AliTRDCalibraVector *calibraVector = 0x0;
  if(fVector2d) calibraVector = (AliTRDCalibraVector *) fContainer->FindObject("CalibraVector");
  //
  // GAIN TH2I
  //
  Bool_t pass = kFALSE; 
  AliTRDCalibraVector *vvect = 0x0;
  if(fHisto2d) {
    TH2I *histogain = (TH2I *) fContainer->FindObject("CH2d");  
    if(histogain) {
      histogain->SetDirectory(0);
      calibra->SetMinEntries(20); 
      if(fCompressPerDetector){
	if(AddStatsPerDetector(histogain)) pass = calibra->AnalyseCH(fCHSum);
      }
      else pass = calibra->AnalyseCH(histogain);
    }
  }
  else {
    if(fVector2d && calibraVector) {
      calibra->SetMinEntries(20); 
      if(fCompressPerDetector){
	if(!(vvect = calibraVector->AddStatsPerDetectorCH())) return kFALSE;
       	pass = calibra->AnalyseCH(vvect);
      }
      else pass = calibra->AnalyseCH(calibraVector);
    }
  }
  
  if(pass)
    {
      Int_t nbtg = 6*4*18*((Int_t) ((AliTRDCalibraMode *)calibra->GetCalibraMode())->GetDetChamb0(0))
	+ 6*  18*((Int_t) ((AliTRDCalibraMode *)calibra->GetCalibraMode())->GetDetChamb2(0));
      Int_t nbfit = calibra->GetNumberFit();  //Number of fits
      Int_t nbE   = calibra->GetNumberEnt();  //Number of detector mit entries
      // enough statistics
      if ((nbtg >                  0) && 
	  (nbfit        >= 0.001*nbE))
	{
	  // create the cal objects
	  calibra->PutMeanValueOtherVectorFit(1,kTRUE); 
	  TObjArray object           = calibra->GetVectorFit();
	  AliTRDCalDet *objgaindet   = calibra->CreateDetObjectGain(&object);
	  TObject *objgainpad        = calibra->CreatePadObjectGain();
	  // store
	  fArrayCalib->AddAt(objgaindet,0);
	  fArrayCalib->AddAt(objgainpad,1);
	  storage[0] = kTRUE;
	  // Make graph
	  TGraph *graph = 0x0;
	  if(FillGraphIndex(&object,graph)){ 
	    fGraph->AddAt(graph,0);
	  }
	}//if(enough statistics?)
      calibra->ResetVectorFit();
    }
  else return kFALSE;
  
  //
  // VDRIFT average pulse height
  //
  pass = kFALSE; 
  if(fHisto2d) {
    TProfile2D *histodriftvelocity = (TProfile2D *) fContainer->FindObject("PH2d");  
    if(histodriftvelocity) {
      histodriftvelocity->SetDirectory(0);  
      calibra->SetMinEntries(20*20);  
      if(fCompressPerDetector){
	if(AddStatsPerDetector(histodriftvelocity,1)) {
	  pass = calibra->AnalysePH(fDetSumVector);
	}
      }
      else pass = calibra->AnalysePH(histodriftvelocity);
    }
  }
  else {
    if(fVector2d && calibraVector) {
      calibra->SetMinEntries(20*20);  
      if(fCompressPerDetector){
	if(!(vvect = calibraVector->AddStatsPerDetectorPH())) return kFALSE;
       	pass = calibra->AnalysePH(vvect);
      }
      else pass = calibra->AnalysePH(calibraVector);  
    }
  }

  if(pass) {
    Int_t nbtg = 6*4*18*((Int_t) ((AliTRDCalibraMode *)calibra->GetCalibraMode())->GetDetChamb0(1))
      + 6*  18*((Int_t) ((AliTRDCalibraMode *)calibra->GetCalibraMode())->GetDetChamb2(1));
    Int_t nbfit  = calibra->GetNumberFit();
    Int_t nbE    = calibra->GetNumberEnt();
    // enough statistics
    if ((nbtg >                  0) && 
	(nbfit        >= 0.001*nbE))
      {
	// create the cal objects
	calibra->PutMeanValueOtherVectorFit(1,kTRUE);
	calibra->PutMeanValueOtherVectorFit2(1,kTRUE); 
	TObjArray object  = calibra->GetVectorFit();
	AliTRDCalDet *objdriftvelocitydet = calibra->CreateDetObjectVdrift(&object,kTRUE);
	TObject *objdriftvelocitypad      = calibra->CreatePadObjectVdrift();
	TObjArray objectt          = calibra->GetVectorFit2();
	AliTRDCalDet *objtime0det  = calibra->CreateDetObjectT0(&object,kTRUE);
	TObject *objtime0pad       = calibra->CreatePadObjectT0();
	// store
	fArrayCalib->AddAt(objdriftvelocitydet,2);
	fArrayCalib->AddAt(objdriftvelocitypad,3);
	//
	fArrayCalib->AddAt(objtime0det,4);
	fArrayCalib->AddAt(objtime0pad,5);
	// Make graph
	TGraph *graph = 0x0;
	if(FillGraphIndex(&object,graph)){ 
	  fGraph->AddAt(graph,1);
	}
	TGraph *graphh = 0x0;
	if(FillGraphIndex(&objectt,graphh)){ 
	  fGraph->AddAt(graphh,2);
	}
      }//if(enough statistics)
    calibra->ResetVectorFit();
  }
  else return kFALSE;
  
  //
  // PRF
  //
  pass = kFALSE; 
  if(fHisto2d) {
    TProfile2D *histoprf = (TProfile2D *) fContainer->FindObject("PRF2d");
    if (histoprf) {
      histoprf->SetDirectory(0);  
      calibra->SetMinEntries(600); 
      if(fCompressPerDetector){
    	if(AddStatsPerDetector(histoprf,2)) pass = calibra->AnalysePRFMarianFit(fDetSumVector);
      }
      else pass = calibra->AnalysePRFMarianFit(histoprf);
    }
  }
  else {
    if(fVector2d && calibraVector) {
      calibra->SetMinEntries(600);  
      if(fCompressPerDetector){
	if(!(vvect =calibraVector->AddStatsPerDetectorPRF())) return kFALSE;
	pass = calibra->AnalysePRFMarianFit(vvect);
      }
      else pass = calibra->AnalysePRFMarianFit(calibraVector);  
    }
  }
  
  if(pass){
    Int_t nbtg = 6*4*18*((Int_t) ((AliTRDCalibraMode *)calibra->GetCalibraMode())->GetDetChamb0(2))
      + 6*  18*((Int_t) ((AliTRDCalibraMode *)calibra->GetCalibraMode())->GetDetChamb2(2));
    Int_t nbfit        = calibra->GetNumberFit();
    Int_t nbE          = calibra->GetNumberEnt();
    // enough statistics
    if ((nbtg >                  0) && 
	(nbfit        >= 0.001*nbE)) {
      TObjArray object            = calibra->GetVectorFit();
      TObject *objPRFpad          = calibra->CreatePadObjectPRF(&object);
      // store
      fArrayCalib->AddAt(objPRFpad,6);
      // Make graph
      TGraph *graph = 0x0;
      if(FillGraphIndex(&object,graph)){ 
	fGraph->AddAt(graph,5);
      }
    }
    calibra->ResetVectorFit();
  }
  else return kFALSE;
  
  //
  // VDRIFT linear fit 
  //
  AliTRDCalibraVdriftLinearFit *vlinearfit = (AliTRDCalibraVdriftLinearFit *) fContainer->FindObject("LinearVdriftFit"); 
  if (vlinearfit) {
    calibra->SetMinEntries(20*20);     
    if(calibra->AnalyseLinearFitters(vlinearfit)) {
      
      Int_t nbtg = 6*4*18*((Int_t) ((AliTRDCalibraMode *)calibra->GetCalibraMode())->GetDetChamb0(2))
	+ 6*  18*((Int_t) ((AliTRDCalibraMode *)calibra->GetCalibraMode())->GetDetChamb2(2));
      Int_t nbfit  = calibra->GetNumberFit();
      Int_t nbE    = calibra->GetNumberEnt();
      // enough statistics
      if ((nbtg >                  0) && 
	  (nbfit        >= 0.001*nbE))
	{
	  // create the cal objects
	  calibra->PutMeanValueOtherVectorFit(1,kTRUE);
	  calibra->PutMeanValueOtherVectorFit2(1,kTRUE); 
	  TObjArray object  = calibra->GetVectorFit();
	  AliTRDCalDet *objdriftvelocitydet = calibra->CreateDetObjectVdrift(&object,kTRUE);
	  TObject *objdriftvelocitypad      = calibra->CreatePadObjectVdrift();
	  TObjArray objectt          = calibra->GetVectorFit2();
	  AliTRDCalDet *objtime0det  = calibra->CreateDetObjectT0(&object,kTRUE);
	  TObject *objtime0pad       = calibra->CreatePadObjectT0();
	  // store dummy
	  fArrayCalib->AddAt(objdriftvelocitydet,7);
	  fArrayCalib->AddAt(objdriftvelocitypad,8);
	  //
	  fArrayCalib->AddAt(objtime0det,9);
	  fArrayCalib->AddAt(objtime0pad,10);
	  // Make graph
	  TGraph *graph = 0x0;
	  if(FillGraphIndex(&object,graph)){ 
	    fGraph->AddAt(graph,3);
	  }
	  TGraph *graphh = 0x0;
	  if(FillGraphIndex(&objectt,graphh)){ 
	    fGraph->AddAt(graphh,4);
	  }
	}//if(enough statistics)
    }// if fit
    calibra->ResetVectorFit();
  }
  else return kFALSE;
  
  fPostProcess = kTRUE;
  
  return kTRUE;
  
}

//________________________________________________________________________
Bool_t AliTRDcalibration::FillGraphIndex(const TObjArray *vectora,TGraph *graph) const
{
  //
  // Fill one value (the first one) per detector
  //

  Int_t loop = (Int_t) vectora->GetEntriesFast();
  if(loop != 540) {
    AliInfo("The Vector Fit is not complete!");
    return kFALSE;
  }
  
  Double_t x[540];
  Double_t y[540];
  for (Int_t k = 0; k < loop; k++) {
    if(!vectora->At(k)){
      x[k] = -1.0;
      y[k] = -1.0;
      continue;
    }
    x[k]  = ((AliTRDCalibraFit::AliTRDFitInfo *) vectora->At(k))->GetDetector();
    y[k]  = ((Float_t *)((AliTRDCalibraFit::AliTRDFitInfo *) vectora->At(k))->GetCoef())[0];
  }

  if(!graph){
    graph = new TGraph(540,&x[0],&y[0]);
    graph->SetMarkerStyle(20);
  } else{ 
    graph->~TGraph();
    new(graph) TGraph(540,&x[0],&y[0]);
  }

  return kTRUE;

}
//________________________________________________________________________
Bool_t AliTRDcalibration::AddStatsPerDetector(const TH2I *ch) 
{
  //
  // Add statistic per detector
  //
  
  AliTRDCalibraMode calibMode = AliTRDCalibraMode();
  const char *name = ch->GetTitle();
  if(!SetNzFromTObject(name,0,&calibMode)) return 0x0;
  if(!SetNrphiFromTObject(name,0,&calibMode)) return 0x0;
  if(((calibMode.GetNz(0) == 100) && (calibMode.GetNrphi(0) == 100)) || ((calibMode.GetNz(0) == 10) && (calibMode.GetNrphi(0) == 10))) return kFALSE;

  Int_t    nybins  = ch->GetNbinsY();// groups number
  Int_t    nxbins  = ch->GetNbinsX();// number of bins X
  TAxis   *xaxis   = ch->GetXaxis();
  Double_t lowedge  = xaxis->GetBinLowEdge(1);
  Double_t upedge   = xaxis->GetBinUpEdge(nxbins);

  // number per chamber 2
  calibMode.ModePadCalibration(2,0);
  calibMode.ModePadFragmentation(0,2,0,0);
  calibMode.SetDetChamb2(0);
  Int_t perChamber2 = (Int_t) calibMode.GetDetChamb2(0);

  // number per other chamber
  calibMode.ModePadCalibration(0,0);
  calibMode.ModePadFragmentation(0,0,0,0);
  calibMode.SetDetChamb0(0);
  Int_t perChamber0 = (Int_t) calibMode.GetDetChamb0(0);

  if(nybins != (6*18*perChamber2+6*4*18*perChamber0)) return kFALSE;

  // Create Histo
  TString nname((const char *)ch->GetName());
  nname  += "PerDetector";
  TString title("Nz");
  title += 0;
  title += "Nrphi";
  title += 0;
  if(!fCHSum) fCHSum = new TH2I((const char *)nname,(const char *)title
				,nxbins,lowedge,upedge,540,0,540);
  else{ 
    fCHSum->~TH2I();
    new(fCHSum) TH2I((const Char_t *) nname,(const char *)title
		     ,nxbins,lowedge,upedge,540,0,540);
  }
  fCHSum->SetYTitle("Detector number");
  fCHSum->SetXTitle("charge deposit [a.u]");
  fCHSum->SetZTitle("counts");
  fCHSum->SetStats(0);
  fCHSum->Sumw2();

  Int_t counter = 0;
  
  for(Int_t det = 0; det < 540; det++){

    Int_t numberofgroup = 0;
    if(AliTRDgeometry::GetStack(det) == 2) numberofgroup = perChamber2;
    else numberofgroup = perChamber0;
    TH1I *projch = (TH1I *) ch->ProjectionX("projch",counter+1,counter+numberofgroup,(Option_t *)"e");
    projch->SetDirectory(0);
       
    for(Int_t nx = 0; nx <= nxbins; nx++) {
      fCHSum->SetBinContent(nx,det+1,projch->GetBinContent(nx));
      fCHSum->SetBinError(nx,det+1,projch->GetBinError(nx));
    }

    counter += numberofgroup;
    
    delete projch;

  }

  return kTRUE;

}
//_____________________________________________________________________________________________________________________
Bool_t AliTRDcalibration::AddStatsPerDetector(const TProfile2D *ph,Int_t i)
{
  //
  // Add statistic per detector
  //

  AliTRDCalibraMode calibMode = AliTRDCalibraMode();
  const char *name = ph->GetTitle();
  //printf("name %s\n",name);
  if(!SetNzFromTObject(name,0,&calibMode)) return kFALSE;
  if(!SetNrphiFromTObject(name,0,&calibMode)) return kFALSE;
  if(((calibMode.GetNz(0) == 100) && (calibMode.GetNrphi(0) == 100)) || ((calibMode.GetNz(0) == 10) && (calibMode.GetNrphi(0) == 10))) return kFALSE;
  //printf("Found mode Mz %d, Nrphi %d\n",calibMode.GetNz(0),calibMode.GetNrphi(0));  


  Int_t    nybins  = ph->GetNbinsY();// groups number
  Int_t    nxbins  = ph->GetNbinsX();// number of bins X
  TAxis   *xaxis = ph->GetXaxis();
  Double_t lowedge  = xaxis->GetBinLowEdge(1);
  Double_t upedge   = xaxis->GetBinUpEdge(nxbins);

  // number per chamber 2
  calibMode.ModePadCalibration(2,0);
  calibMode.ModePadFragmentation(0,2,0,0);
  calibMode.SetDetChamb2(0);
  Int_t perChamber2 = (Int_t) calibMode.GetDetChamb2(0);

  // number per other chamber
  calibMode.ModePadCalibration(0,0);
  calibMode.ModePadFragmentation(0,0,0,0);
  calibMode.SetDetChamb0(0);
  Int_t perChamber0 = (Int_t) calibMode.GetDetChamb0(0);

  if(nybins != (6*18*perChamber2+6*4*18*perChamber0)) return kFALSE;
  
  // Create calvector 
  TString nbname((const char *)ph->GetName());
  nbname  += "PerDetectorVector";
  if(!fDetSumVector) fDetSumVector = new AliTRDCalibraVector();
  else{ 
    fDetSumVector->~AliTRDCalibraVector();
    new(fDetSumVector) AliTRDCalibraVector();
  }
  if(i==1){
    fDetSumVector->SetTimeMax(nxbins);
  }
  if(i==2){
    fDetSumVector->SetNumberBinPRF(nxbins);
    fDetSumVector->SetPRFRange(TMath::Abs(lowedge));
  }
  fDetSumVector->SetDetCha0(i,1);
  fDetSumVector->SetDetCha2(i,1);
  fDetSumVector->SetNzNrphi(i,0,0);
  if(i==2) {
    Int_t nbg = GetNumberOfGroupsPRF((const char *)name);
    fDetSumVector->SetNbGroupPRF(nbg);
  }

  // Create Histo
  TString nname((const char *)ph->GetName());
  nname  += "PerDetector";
  TString title("Nz");
  title += 0;
  title += "Nrphi";
  title += 0;
  if(!fDetSum) fDetSum = new TH2D((const char *)nname,(const Char_t *) title
				,nxbins,lowedge,upedge,540,0,540);
  else{ 
    fDetSum->~TH2D();
    new(fDetSum) TH2D((const Char_t *) nname,(const Char_t *) title
		     ,nxbins,lowedge,upedge,540,0,540);
  }
  fDetSum->SetYTitle("Detector number");
  fDetSum->SetXTitle(xaxis->GetTitle());
  fDetSum->SetStats(0);
  
  Int_t counter = 0;

  for(Int_t det = 0; det < 540; det++){

    Int_t numberofgroup = 0;
    if(AliTRDgeometry::GetStack(det) == 2) numberofgroup = perChamber2;
    else numberofgroup = perChamber0;
    
    for(Int_t nx = 1; nx <= nxbins; nx++) {
      
      Double_t entries = 0.0;
      Double_t sumw2 = 0.0;
      Double_t sumw = 0.0;

      for(Int_t k = counter+1; k <= (counter+numberofgroup); k++){
	Int_t  binnumber = ph->GetBin(nx,k);
	entries += ph->GetBinEntries(binnumber);
	sumw2 += (ph->GetBinError(binnumber)*ph->GetBinError(binnumber)+ph->GetBinContent(binnumber)*ph->GetBinContent(binnumber))*ph->GetBinEntries(binnumber);
	sumw += ph->GetBinContent(binnumber)*ph->GetBinEntries(binnumber);
      }

      Double_t mean = 0.0;
      if(entries > 0.0) mean = sumw/entries;
      Double_t squaremean = 0.0;
      if(entries > 0.0) squaremean = sumw2/entries;
      Double_t errorf = squaremean - mean*mean;
      Double_t error = 0.0;
      if(entries > 0.0) error = TMath::Sqrt(TMath::Abs(errorf)/entries);
      
      fDetSum->SetBinContent(nx,det+1,mean);
      fDetSum->SetBinError(nx,det+1,error);

      if(i==1) fDetSumVector->FillVectorPH(det,0,nx-1,(Int_t)entries,(Float_t)mean,(Float_t)squaremean);
      if(i==2) fDetSumVector->FillVectorPRF(det,0,nx-1,(Int_t)entries,(Float_t)mean,(Float_t)squaremean);
      
    }
    
    counter += numberofgroup;

  }

  return kTRUE;

  
}
//_____________________________________________________________________________
Bool_t AliTRDcalibration::SetNrphiFromTObject(const char *name, Int_t i, AliTRDCalibraMode *calibMode) const
{
  //
  // Set the granularity from object
  //  
  
  const Char_t *patternrphi0 = "Nrphi0";
  const Char_t *patternrphi1 = "Nrphi1";
  const Char_t *patternrphi2 = "Nrphi2";
  const Char_t *patternrphi3 = "Nrphi3";
  const Char_t *patternrphi4 = "Nrphi4";
  const Char_t *patternrphi5 = "Nrphi5";
  const Char_t *patternrphi6 = "Nrphi6";

  
  const Char_t *patternrphi10 = "Nrphi10";
  const Char_t *patternrphi100 = "Nrphi100";
  const Char_t *patternz10 = "Nz10";
  const Char_t *patternz100 = "Nz100";

  // Nrphi mode
  if ((strstr(name,patternrphi100)) && (strstr(name,patternz100))) {
    calibMode->SetAllTogether(i);
    return kTRUE;
  }
  if ((strstr(name,patternrphi10)) && (strstr(name,patternz10))) {
    calibMode->SetPerSuperModule(i);
    return kTRUE;
  }
  
  if (strstr(name,patternrphi0)) {
    calibMode->SetNrphi(i ,0);
    return kTRUE;
  }
  if (strstr(name,patternrphi1)) {
    calibMode->SetNrphi(i, 1);
    return kTRUE;
  }
  if (strstr(name,patternrphi2)) {
    calibMode->SetNrphi(i, 2);
    return kTRUE;
  }
  if (strstr(name,patternrphi3)) {
    calibMode->SetNrphi(i, 3);
    return kTRUE;
  }
  if (strstr(name,patternrphi4)) {
    calibMode->SetNrphi(i, 4);
    return kTRUE;
  }
  if (strstr(name,patternrphi5)) {
    calibMode->SetNrphi(i, 5);
    return kTRUE;
  }
  if (strstr(name,patternrphi6)) {
    calibMode->SetNrphi(i, 6);
    return kTRUE;
  }
  
  calibMode->SetNrphi(i ,0);
  return kFALSE;
  
}
//_____________________________________________________________________________
Bool_t AliTRDcalibration::SetNzFromTObject(const char *name, Int_t i, AliTRDCalibraMode *calibMode) const
{
  //
  // Set fNz[i] of the AliTRDCalibraFit::Instance()
  // corresponding to the given TObject
  //

  // Some patterns
  const Char_t *patternz0    = "Nz0";
  const Char_t *patternz1    = "Nz1";
  const Char_t *patternz2    = "Nz2";
  const Char_t *patternz3    = "Nz3";
  const Char_t *patternz4    = "Nz4";

  const Char_t *patternrphi10 = "Nrphi10";
  const Char_t *patternrphi100 = "Nrphi100";
  const Char_t *patternz10 = "Nz10";
  const Char_t *patternz100 = "Nz100";

  if ((strstr(name,patternrphi100)) && (strstr(name,patternz100))) {
    calibMode->SetAllTogether(i);
    return kTRUE;
  }
  if ((strstr(name,patternrphi10)) && (strstr(name,patternz10))) {
    calibMode->SetPerSuperModule(i);
    return kTRUE;
  }
  if (strstr(name,patternz0)) {
    calibMode->SetNz(i, 0);
    return kTRUE;
  }
  if (strstr(name,patternz1)) {
    calibMode->SetNz(i ,1);
    return kTRUE;
  }
  if (strstr(name,patternz2)) {
    calibMode->SetNz(i ,2);
    return kTRUE;
  }
  if (strstr(name,patternz3)) {
    calibMode->SetNz(i ,3);
    return kTRUE;  
  }
  if (strstr(name,patternz4)) {
    calibMode->SetNz(i ,4);
    return kTRUE;
  }
 
  calibMode->SetNz(i ,0);
  return kFALSE;
}
//____________________________________________________________________________________________________
Int_t AliTRDcalibration::GetNumberOfGroupsPRF(const char* nametitle) const
{
  //
  // Get numberofgroupsprf
  //
  
  // Some patterns
  const Char_t *pattern0 = "Ngp0";
  const Char_t *pattern1 = "Ngp1";
  const Char_t *pattern2 = "Ngp2";
  const Char_t *pattern3 = "Ngp3";
  const Char_t *pattern4 = "Ngp4";
  const Char_t *pattern5 = "Ngp5";
  const Char_t *pattern6 = "Ngp6";

  // Nrphi mode
  if (strstr(nametitle,pattern0)) {
    return 0;
  }
  if (strstr(nametitle,pattern1)) {
    return 1;
  }
  if (strstr(nametitle,pattern2)) {
    return 2;
  }
  if (strstr(nametitle,pattern3)) {
    return 3;
  }
  if (strstr(nametitle,pattern4)) {
    return 4;
  }
  if (strstr(nametitle,pattern5)) {
    return 5;
  }
  if (strstr(nametitle,pattern6)){
    return 6;
  }
  else return -1;
}
