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

/* $Id$ */

//_________________________________________________________________________
// An analysis task to check the TRD data in simulated data
//
//*-- Sylwester Radomski
//////////////////////////////////////////////////////////////////////////////
// track type codding
//
// tpci = kTPCin
// tpco = kTPCout
// tpcz = kTPCout && !kTRDout
// trdo = kTRDout
// trdr = kTRDref
// trdz = kTRDout && !kTRDref
// 

#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h>
#include <TGaxis.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TROOT.h>
#include <TStyle.h>

#include "AliTRDQATask.h"
#include "AliESD.h"
#include "AliLog.h"

//______________________________________________________________________________
AliTRDQATask::AliTRDQATask(const char *name) : 
  AliAnalysisTask(name,""),  
  fChain(0),
  fESD(0)
{
  // Constructor.
  // Input slot #0 works with an Ntuple
  DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TH1 container
  DefineOutput(0,  TObjArray::Class()) ; 
}

//______________________________________________________________________________
void AliTRDQATask::Init(const Option_t *)
{
  // Initialisation of branch container and histograms 

  AliInfo(Form("*** Initialization of %s", GetName())) ; 
  
  // Get input data
  fChain = dynamic_cast<TChain *>(GetInputData(0)) ;
  if (!fChain) {
    AliError(Form("Input 0 for %s not found\n", GetName()));
    return ;
  }

  if (!fESD) {
    // One should first check if the branch address was taken by some other task
    char ** address = (char **)GetBranchAddress(0, "ESD");
    if (address) {
      fESD = (AliESD *)(*address); 
      AliInfo("Old ESD found.");
    }
    if (!fESD) {
      fESD = new AliESD();
      SetBranchAddress(0, "ESD", &fESD);  
      if (fESD) AliInfo("ESD connected.");
    }
  }
  // The output objects will be written to 
  TDirectory * cdir = gDirectory ; 
  OpenFile(0, Form("%s.root", GetName()), "RECREATE"); 
  if (cdir) 
    cdir->cd() ; 

  // create histograms 

  fNTracks     = new TH1D("ntracks", ";number of all tracks", 500, -0.5, 499.5); 
  fEventSize   = new TH1D("evSize", ";event size (MB)", 100, 0, 5);

  fTrackStatus = new TH1D("trackStatus", ";status bit", 32, -0.5, 31.5);
  fKinkIndex   = new TH1D("kinkIndex", ";kink index", 41, -20.5, 20.5);
  
  fParIn  = new TH1D("parIn", "Inner Plane", 2, -0.5, 1.5);
  fParOut = new TH1D("parOut", "Outer Plane", 2, -0.5, 1.5);

  fXIn  = new TH1D("xIn", ";X at the inner plane (cm)", 200, 50, 250);
  fXOut = new TH1D("xOut", ";X at the outer plane (cm)", 300, 50, 400);
  
  const int nNameAlpha = 4;
  const char *namesAlpha[nNameAlpha] = {"alphaTPCi", "alphaTPCo", "alphaTRDo", "alphaTRDr"};
  //TH1D *fAlpha[4];
  for(int i=0; i<nNameAlpha; i++) {
    fAlpha[i] = new TH1D(namesAlpha[i], "alpha", 100, -4, 4);
  }
  fSectorTRD = new TH1D ("sectorTRD", ";sector TRD", 20, -0.5, 19.5);


  // track parameters
  const int nbits = 6;
  const char *suf[nbits] = {"TPCi", "TPCo", "TPCz", "TRDo", "TRDr", "TRDz"};
  for(int i=0; i<nbits; i++) {
    fPt[i]      = new TH1D(Form("pt%s",suf[i]), ";p_{T} (GeV/c);entries TPC", 50, 0, 10);
    fTheta[i]   = new TH1D(Form("theta%s", suf[i]), "theta (rad)", 100, -4, 4); 
    fSigmaY[i]  = new TH1D(Form("sigmaY%s",suf[i]),  ";sigma Y (cm)", 200, 0, 1);
    fChi2[i]    = new TH1D(Form("Chi2%s", suf[i]), ";#chi2 / ndf", 100, 0, 10);
    fPlaneYZ[i] = new TH2D(Form("planeYZ%s", suf[i]), Form("%sy (cm);z (cm)", suf[i]), 
			   100, -60, 60, 500, -500, 500);
  }
  
  // efficiency
  fEffPt[0] = (TH1D*) fPt[0]->Clone(Form("eff_%s_%s", suf[0], suf[1]));
  fEffPt[1] = (TH1D*) fPt[0]->Clone(Form("eff_%s_%s", suf[1], suf[3]));
  fEffPt[2] = (TH1D*) fPt[0]->Clone(Form("eff_%s_%s", suf[3], suf[4]));
  fEffPt[3] = (TH1D*) fPt[0]->Clone(Form("eff_%s_%s", suf[1], suf[4]));

  for(int i=0; i<4; i++) {
    fEffPt[i]->Sumw2();
    fEffPt[i]->SetMarkerStyle(20);
    fEffPt[i]->SetMinimum(0);
    fEffPt[i]->SetMaximum(1.1);
  }

  // track features
  fClustersTRD[0] = new TH1D("clsTRDo", "TRDo;number of clusters", 130, -0.5, 129.5);;
  fClustersTRD[1] = new TH1D("clsTRDr", "TRDr;number of clusters", 130, -0.5, 129.5);;
  fClustersTRD[2] = new TH1D("clsTRDz", "TRDz;number of clusters", 130, -0.5, 129.5);;

  // for good refitted tracks only
  fTime    = new TH1D("time", ";time bin", 25, -0.5, 24.5);
  fBudget  = new TH1D("budget", ";material budget", 100, 0, 100);
  fQuality = new TH1D("quality", ";track quality", 100, 0, 1.1);
  fSignal  = new TH1D("signal", ";signal", 100, 0, 1e3);  
  
  // dEdX and PID
  fTrdSigMom = new TH2D("trdSigMom", ";momentum (GeV/c);signal", 100, 0, 3, 100, 0, 1e3);
  fTpcSigMom = new TH2D("tpcSigMom", ";momentum (GeV/c);signal", 100, 0, 3, 100, 0, 200);
  
  const char *pidName[6] = {"El", "Mu", "Pi", "K", "P", "Ch"};
  for(int i=0; i<6; i++) {
    
    // TPC
    fTpcPID[i] = new TH1D(Form("tpcPid%s",pidName[i]), pidName[i], 100, 0, 1.5);
    fTpcPID[i]->GetXaxis()->SetTitle("probability");
    
    fTpcSigMomPID[i] = new TH2D(Form("tpcSigMom%s",pidName[i]), "", 100, 0, 3, 100, 0, 200);
    fTpcSigMomPID[i]->SetTitle(Form("%s;momentum (GeV/c);signal",pidName[i])); 
    
    // TRD
    fTrdPID[i] = new TH1D(Form("trdPid%s",pidName[i]), pidName[i], 100, 0, 1.5);
    fTrdPID[i]->GetXaxis()->SetTitle("probability");
     
    fTrdSigMomPID[i] = new TH2D(Form("trdSigMom%s",pidName[i]), "", 100, 0, 3, 100, 0, 1e3);
    fTrdSigMomPID[i]->SetTitle(Form("%s;momentum (GeV/c);signal",pidName[i]));  
  }

  
  // create output container
  fOutputContainer = new TObjArray(150); 
  
  // register histograms to the container  
  TIter next(gDirectory->GetList());
  TObject *obj;
  int counter = 0;
  
  while (obj = next.Next()) {
    if (obj->InheritsFrom("TH1"))  fOutputContainer->AddAt(obj, counter++);
  }

  AliInfo(Form("Number of histograms = %d", counter));

 }

//______________________________________________________________________________
void AliTRDQATask::Exec(Option_t *) 
{
  // Process one event
  
   Long64_t entry = fChain->GetReadEntry() ;
  
  // Processing of one event 
   
  if (!fESD) {
    AliError("fESD is not connected to the input!") ; 
    return ; 
  }
  
  if ( !((entry-1)%100) ) 
    AliInfo(Form("%s ----> Processing event # %lld",  (dynamic_cast<TChain *>(fChain))->GetFile()->GetName(), entry)) ; 

  int nTracks = fESD->GetNumberOfTracks();
  fNTracks->Fill(nTracks); 

  // track loop
  for(int i=0; i<nTracks; i++) {
    
    AliESDtrack *track = fESD->GetTrack(i);
    const AliExternalTrackParam *paramOut = track->GetOuterParam();
    const AliExternalTrackParam *paramIn = track->GetInnerParam();

    fParIn->Fill(!!paramIn);
    if (!paramIn) continue;
    fXIn->Fill(paramIn->GetX());

    fParOut->Fill(!!paramOut);
    if (!paramOut) continue;
    fXOut->Fill(paramOut->GetX());
 
    int sector = GetSector(paramOut->GetAlpha());
    if (!CheckSector(sector)) continue;
    fSectorTRD->Fill(sector);

    fKinkIndex->Fill(track->GetKinkIndex(0));
    if (track->GetKinkIndex(0)) continue;    

    UInt_t u = 1;
    UInt_t status = track->GetStatus();
    for(int bit=0; bit<32; bit++) 
      if (u<<bit & status) fTrackStatus->Fill(bit);

    const int nbits = 6; 
    int bit[6] = {0,0,0,0,0,0};    
    bit[0] = status & AliESDtrack::kTPCin;
    bit[1] = status & AliESDtrack::kTPCout;
    bit[2] = (status & AliESDtrack::kTPCout) && !(status & AliESDtrack::kTRDout);
    bit[3] = status & AliESDtrack::kTRDout;
    bit[4] = status & AliESDtrack::kTRDrefit;
    bit[5] = (status & AliESDtrack::kTRDout) && !(status & AliESDtrack::kTRDrefit);

    
    // transverse momentum
    const double *val = track->GetParameter(); // parameters at the vertex
    double pt = 1./TMath::Abs(val[4]);

    for(int b=0; b<nbits; b++) {
      if (bit[b]) {
	fPt[b]->Fill(pt); 
	fTheta[b]->Fill(val[3]);
	fSigmaY[b]->Fill(TMath::Sqrt(paramOut->GetSigmaY2()));
	fChi2[b]->Fill(track->GetTRDchi2()/track->GetTRDncls());    
	fPlaneYZ[b]->Fill(paramOut->GetY(), paramOut->GetZ()); 
      }
    }

    // sectors
    if (bit[1]) {
      fAlpha[0]->Fill(paramIn->GetAlpha());
      fAlpha[1]->Fill(paramOut->GetAlpha());
    }
    
    if (bit[3]) fAlpha[2]->Fill(paramOut->GetAlpha());
    if (bit[4]) fAlpha[3]->Fill(paramOut->GetAlpha());

    // clusters
    for(int b=0; b<3; b++) 
      if (bit[3+b]) fClustersTRD[b]->Fill(track->GetTRDncls());

    // refitted only
    if (!bit[4]) continue;

    fQuality->Fill(track->GetTRDQuality());
    fBudget->Fill(track->GetTRDBudget());
    fSignal->Fill(track->GetTRDsignal());
	
    fTrdSigMom->Fill(track->GetP(), track->GetTRDsignal());
    fTpcSigMom->Fill(track->GetP(), track->GetTPCsignal());

    // PID only
    if (status & AliESDtrack::kTRDpid) {
      
      for(int l=0; l<6; l++) fTime->Fill(track->GetTRDTimBin(l));

      // fill pid histograms
      double trdr0 = 0, tpcr0 = 0;
      int trdBestPid = 5, tpcBestPid = 5;  // charged
      const double minPidValue =  0.9;

      double pp[5];
      track->GetTPCpid(pp); // ESD inconsequence

      for(int pid=0; pid<5; pid++) {
	
	trdr0 += track->GetTRDpid(pid);
	tpcr0 += pp[pid];
	
	fTrdPID[pid]->Fill(track->GetTRDpid(pid));
	fTpcPID[pid]->Fill(pp[pid]);
	
	if (track->GetTRDpid(pid) > minPidValue) trdBestPid = pid;
	if (pp[pid] > minPidValue) tpcBestPid = pid;
      }
      
      fTrdPID[5]->Fill(trdr0); // check unitarity
      fTrdSigMomPID[trdBestPid]->Fill(track->GetP(), track->GetTRDsignal());
      
      fTpcPID[5]->Fill(tpcr0); // check unitarity
      fTpcSigMomPID[tpcBestPid]->Fill(track->GetP(), track->GetTPCsignal());
    }
    
  }

  CalculateEff();
  PostData(0, fOutputContainer);
}

//______________________________________________________________________________
void AliTRDQATask::Terminate(Option_t *)
{
  // Processing when the event loop is ended
  AliInfo("TRD QA module");

  // create efficiency histograms
  
  CalculateEff();
  PostData(0, fOutputContainer);

  DrawESD() ; 
  DrawGeoESD() ; 
  //DrawConvESD() ; 
  DrawPidESD() ; 
}

//______________________________________________________________________________
int AliTRDQATask::GetSector(double alpha) 
{
  // Gets the sector number 

  double size = TMath::DegToRad() * 20.;
  int sector = (int)((alpha + TMath::Pi())/size);
  return sector;
}

//______________________________________________________________________________
int AliTRDQATask::CheckSector(int sector) 
{  
  // Checks the sector number
  const int nSec = 8;
  int sec[] = {2,3,5,6,11,12,13,15};
  
  for(int i=0; i<nSec; i++) 
    if (sector == sec[i]) return 1;
  
  return 0;
}

//______________________________________________________________________________
void AliTRDQATask::CalculateEff() 
{
  // calculates the efficiency
  
  for(int i=0; i<4; i++) fEffPt[i]->Reset();
  
  fEffPt[0]->Add(fPt[1]);
  fEffPt[0]->Divide(fPt[0]);
  
  fEffPt[1]->Add(fPt[3]);
  fEffPt[1]->Divide(fPt[1]);
  
  fEffPt[2]->Add(fPt[4]);
  fEffPt[2]->Divide(fPt[3]);
  
  fEffPt[3]->Add(fPt[4]);
  fEffPt[3]->Divide(fPt[1]);
}

//______________________________________________________________________________
void AliTRDQATask::DrawESD() 
{
  // Makes a few plots

  AliInfo("Plotting....") ; 
  
  TCanvas * cTRD = new TCanvas("cTRD", "TRD ESD Test", 400, 10, 600, 700) ;
  cTRD->Divide(6,3) ;

  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  
  TGaxis::SetMaxDigits(3);
  
  gStyle->SetLabelFont(52, "XYZ");
  gStyle->SetTitleFont(62, "XYZ");
  gStyle->SetPadRightMargin(0.02);

  // draw all 
  
  const int nplots = 18;
  const int nover[nplots] = {1,1,1,4,1,1,1,1,1,1,2,1,1,3,1,1,1,1};
  const int nnames = 24;
  const char *names[nnames] = {
    "ntracks", "kinkIndex", "trackStatus", 
    "ptTPCi", "ptTPCo", "ptTRDo", "ptTRDr",  "ptTPCz", "ptTRDz",
    "eff_TPCi_TPCo",  "eff_TPCo_TRDo", "eff_TRDo_TRDr",  "eff_TPCo_TRDr",
    "clsTRDo", "clsTRDr", "clsTRDz", 
    "alphaTPCi", "alphaTPCo", "alphaTRDo", "alphaTRDr", "sectorTRD",
    "time", "budget", "signal"
  };
  
  const int logy[nnames] = {
    1,1,1,
    1,1,1,
    0,0,0,0,
    1,1,
    0,0,0,0,0,
    0,1,1
  };

  int nhist=0;
  for(int i=0; i<nplots; i++) {
  cTRD->cd(i+1) ;
  
  //  new TCanvas(names[i], names[nhist], 500, 300);
    gPad->SetLogy(logy[i]);
      
    for(int j=0; j<nover[i]; j++) {
      TH1D *hist = dynamic_cast<TH1D*>(gDirectory->FindObject(names[nhist++]));
      if (!hist) continue;
      
      if (strstr(hist->GetName(), "eff")) {
	hist->SetMarkerStyle(20);
	hist->SetMinimum(0);
	hist->SetMaximum(1.2);
      }

      if (!j) hist->Draw();
      else hist->Draw("SAME");
    }
  }
  cTRD->Print("TRD_ESD.gif");
}

//______________________________________________________________________________
void AliTRDQATask::DrawGeoESD() 
{
  // Makes a few plots

  AliInfo("Plotting....") ; 

  TCanvas * cTRDGeo = new TCanvas("cTRDGeo", "TRD ESDGeo Test", 400, 10, 600, 700) ;
  cTRDGeo->Divide(4,2) ;

  gStyle->SetOptStat(0);
  TGaxis::SetMaxDigits(3);
  
  gStyle->SetLabelFont(52, "XYZ");
  gStyle->SetTitleFont(62, "XYZ");
  
  gStyle->SetPadTopMargin(0.06);
  gStyle->SetTitleBorderSize(0);
  
  // draw all   
  const int nnames = 7;
  const char *names[nnames] = {
    "xIn", "xOut",
    "planeYZTPCo", "planeYZTPCz", "planeYZTRDo", "planeYZTRDr", "planeYZTRDz",
  };
  
  const char *opt[nnames] = {
    "", "",
    "colz","colz", "colz", "colz", "colz"
  };
  
  const int logy[nnames] = {
    1,1,
    0,0,0,0,0
  };
  
  for(int i=0; i<nnames; i++) {
  cTRDGeo->cd(i+1) ;
    TH1D *hist = dynamic_cast<TH1D*>(gDirectory->FindObject(names[i]));
    if (!hist) continue;
    
    //if (i<2) new TCanvas(names[i], names[i], 500, 300);
    //else new TCanvas(names[i], names[i], 300, 900);
   
    gPad->SetLogy(logy[i]);
    if (strstr(opt[i],"colz")) gPad->SetRightMargin(0.1);
    
    hist->Draw(opt[i]);    
    AliInfo(Form("%s\t%d", names[i], hist->GetEntries()));
  }
  
  cTRDGeo->Print("TRD_Geo.gif");
}

//______________________________________________________________________________
void AliTRDQATask::DrawConvESD() 
{
  // Makes a few plots

  AliInfo("Plotting....") ; 
  TCanvas * cTRDConv = new TCanvas("cTRDConv", "TRD ESDConv Test", 400, 10, 600, 700) ;
  cTRDConv->Divide(3,2) ;

  gROOT->SetStyle("Plain");
  gROOT->ForceStyle();
  gStyle->SetPalette(1);
  
  TGaxis::SetMaxDigits(3);
  
  gStyle->SetLabelFont(52, "XYZ");
  gStyle->SetTitleFont(62, "XYZ");
  gStyle->SetPadRightMargin(0.02);

  const int nnames = 9;
  const int nplots = 5;
  const int nover[nplots] = {3,1,1,3,1}; 
  
  const char *names[nnames] = {
    "sigmaYTPCo","sigmaYTRDo", "sigmaYTRDr", "sigmaYTPCz", "sigmaYTRDz",
    "Chi2TPCo", "Chi2TRDo", "Chi2TRDr", "Chi2TRDz"
  };
  
  const char *opt[nplots] = {
    "", "", "","","",
  };
  
  const int logy[nplots] = {
    0,0,0,1,1
  };

  int nhist = 0;
  for(int i=0; i<nplots; i++) {
    cTRDConv->cd(i+1) ;
    //new TCanvas(names[i], names[i], 500, 300);
    gPad->SetLogy(logy[i]);
    if (strstr(opt[i],"colz")) gPad->SetRightMargin(0.1);
   
    for(int j=0; j<nover[i]; j++) {
      TH1D *hist = dynamic_cast<TH1D*>(gDirectory->FindObject(names[nhist++]));
      if (!j) hist->Draw(opt[i]);
      else hist->Draw("same");
    }

  }
    cTRDConv->Print("TRD_Conv.eps");
}

//______________________________________________________________________________
void AliTRDQATask::DrawPidESD() 
{
  // Makes a few plots

  AliInfo("Plotting....") ; 
  TCanvas * cTRDPid = new TCanvas("cTRDPid", "TRD ESDPid Test", 400, 10, 600, 700) ;
  cTRDPid->Divide(9,3) ;

  gROOT->SetStyle("Plain");
  gROOT->ForceStyle();
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  
  TGaxis::SetMaxDigits(3);
  
  gStyle->SetLabelFont(52, "XYZ");
  gStyle->SetTitleFont(62, "XYZ");

  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadRightMargin(0.02);

  // draw all 
 
  const int nnames = 27;
  const char *names[nnames] = {

    "signal", "trdSigMom", "tpcSigMom",

    "trdPidEl", "trdPidMu", "trdPidPi", "trdPidK", "trdPidP", "trdPidCh",
    "trdSigMomEl", "trdSigMomMu", "trdSigMomPi", "trdSigMomK", "trdSigMomP", "trdSigMomCh",
    
    "tpcPidEl", "tpcPidMu", "tpcPidPi", "tpcPidK", "tpcPidP", "tpcPidCh",
    "tpcSigMomEl", "tpcSigMomMu", "tpcSigMomPi", "tpcSigMomK", "tpcSigMomP", "tpcSigMomCh"

  };
  
  const char *opt[nnames] = {

    "", "colz", "colz",

    "", "", "", "", "", "" ,
    "colz", "colz", "colz", "colz", "colz", "colz",

    "", "", "", "", "", "" ,
    "colz", "colz", "colz", "colz", "colz", "colz" 
  };
  
  const int logy[nnames] = {

    0,0,0,

    1,1,1,1,1,1,
    0,0,0,0,0,0,

    1,1,1,1,1,1,
    0,0,0,0,0,0    
  };

  for(int i=0; i<nnames; i++) {
    cTRDPid->cd(i+1) ;

    TH1D *hist = dynamic_cast<TH1D*>(gDirectory->FindObject(names[i]));
    if (!hist) continue;
    
    //new TCanvas(names[i], names[i], 500, 300);
    gPad->SetLogy(logy[i]);
    if (strstr(opt[i],"colz")) gPad->SetRightMargin(0.1);
    
    if (strstr(names[i],"sigMom")) gPad->SetLogz(1);
    if (strstr(names[i],"SigMom")) gPad->SetLogz(1);

    hist->Draw(opt[i]);    
    AliInfo(Form("%s\t%d", names[i], hist->GetEntries()));
  }
   cTRDPid->Print("TRD_Pid.gif");
}
