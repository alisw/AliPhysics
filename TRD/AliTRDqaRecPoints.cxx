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

/* $Id: AliTRDqaRecPoints.cxx 23387 2008-01-17 17:25:16Z cblume $ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  Produces the data needed to calculate the quality assurance.          //
//  All data must be mergeable objects.                                   //
//                                                                        //
//  Author:                                                               //
//    Sylwester Radomski (radomski@physi.uni-heidelberg.de)               //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

// --- ROOT system ---
#include <TClonesArray.h>
#include <TFile.h> 
#include <TH1D.h> 
#include <TH2D.h>
#include <TH3D.h>
#include <TProfile.h>
#include <TF1.h>
#include <TCanvas.h>

// --- AliRoot header files ---
#include "AliESDEvent.h"
#include "AliLog.h"
#include "AliTRDcluster.h"
#include "AliTRDqaRecPoints.h"
#include "AliTRDgeometry.h"

#include "AliQAChecker.h"

ClassImp(AliTRDqaRecPoints)

//____________________________________________________________________________ 

AliTRDqaRecPoints::AliTRDqaRecPoints() : 
  TObject(),
  fnEvents(0),
  fHist(0),
  fnPad(0),
  fRef(0)
{
  //
  // Default constructor
  //

}

//____________________________________________________________________________ 

AliTRDqaRecPoints::AliTRDqaRecPoints(const AliTRDqaRecPoints &/*qa*/) : 
  TObject(),
  fnEvents(0),
  fHist(0),
  fnPad(0),
  fRef(0)
{
  //
  // Copy constructor
  //

}

//____________________________________________________________________________ 
void AliTRDqaRecPoints::Process(const char* filename)
{
  //
  // Detector specific actions at end of cycle
  //
  //TStopwatch watch;
  //watch.Start();
  
  AliInfo("End of TRD cycle");
  
  //if (task == AliQAv1::kRECPOINTS) {
  
  TH1D *hist = new TH1D("fitHist", "", 200, -0.5, 199.5);
  //fHist->Print();
  
  // fill detector map;
  for(int i=0; i<540; i++) {
    Double_t v = ((TH1D*)fHist->At(0))->GetBinContent(i+1);
    Int_t sm = i/30;
    Int_t det = i%30;
    
    TH2D *detMap = (TH2D*)fHist->At(87);
    Int_t bin = detMap->FindBin(sm, det);
    detMap->SetBinContent(bin, v);
  }
  
  
  // Rec points full chambers
  for (Int_t i=0; i<540; i++) {
    
    //AliInfo(Form("I = %d", i));
    
    //TH1D *h = ((TH2D*)fHist->At(1))->ProjectionY(Form("qaTRD_recPoints_amp_%d",i), i+1, i+1);
    hist->Reset();
    for(Int_t b=1; b<hist->GetXaxis()->GetNbins()-1; b++) {
      Double_t xvalue = hist->GetBinCenter(b);
      Int_t bin = ((TH2D*)fHist->At(1))->FindBin(i,xvalue);
      Double_t value =  ((TH2D*)fHist->At(1))->GetBinContent(bin);
      hist->SetBinContent(b, value);
    }
    
    //printf("Sum = %d %f\n", i, hist->GetSum());
    if (hist->GetSum() < 100) continue; // chamber not present
    
    hist->Fit("landau", "q0", "goff", 10, 180);
    TF1 *fit = hist->GetFunction("landau");
    ((TH1D*)fHist->At(12))->Fill(fit->GetParameter(1));
    ((TH1D*)fHist->At(13))->Fill(fit->GetParameter(2));
  }
  
  // time-bin by time-bin sm by sm
  for(Int_t i=0; i<18; i++) { // loop over super-modules
    
    for(Int_t j=0; j<35; j++) { // loop over time bins
      
      //TH1D *h =  ((TH3D*)fHist->At(10))->ProjectionZ(Form("ampTime_%d",i), i+1, i+1, j+1, j+1);     
      hist->Reset();
      for(Int_t b=1; b<hist->GetXaxis()->GetNbins()-1; b++) {
	Double_t xvalue = hist->GetBinCenter(b);
	Double_t svalue = 0;
	
	for(Int_t det=i*30; det<(i+1)*30; det++) { // loop over detectors
	  Int_t bin = ((TH3D*)fHist->At(10))->FindBin(det,j,xvalue);
	  Double_t value =  ((TH3D*)fHist->At(10))->GetBinContent(bin);
	  svalue += value;
	}
	//printf("v = %f\n", value);
	hist->SetBinContent(b, svalue);
      }
      
      if (hist->GetSum() < 100) continue;
      //printf("fitting %d %d %f\n", i, j, hist->GetSum());
      
      hist->Fit("landau", "q0", "goff", 10, 180);
      TF1 *fit = hist->GetFunction("landau");
      
      TH1D *h1 = (TH1D*)fHist->At(14+18+i);
      Int_t bin = h1->FindBin(j);
      // printf("%d %d %d\n", det, j, bin);
      h1->SetBinContent(bin, TMath::Abs(fit->GetParameter(1)));
    }
  }
  

  // time-bin by time-bin chamber by chamber

  for (Int_t i=0; i<540; i++) {
    
    //TH1D *test = ((TH3D*)fHist->At(10))->ProjectionZ(Form("ampTime_%d",i), i+1, i+1, 0, 35);     
    //if (test->GetSum() < 100) continue;
    
    //AliInfo(Form("fitting det = %d", i));
    
    for(Int_t j=0; j<35; j++) {
      
      //TH1D *h =  ((TH3D*)fHist->At(10))->ProjectionZ(Form("ampTime_%d",i), i+1, i+1, j+1, j+1);     
      hist->Reset();
      for(Int_t b=1; b<hist->GetXaxis()->GetNbins()-1; b++) {
	Double_t xvalue = hist->GetBinCenter(b);
	Int_t bin = ((TH3D*)fHist->At(10))->FindBin(i,j,xvalue);
	Double_t value =  ((TH3D*)fHist->At(10))->GetBinContent(bin);
	//printf("v = %f\n", value);
	hist->SetBinContent(b, value);
      }
      
      if (hist->GetSum() < 100) continue;
      //printf("fitting %d %d %f\n", i, j, hist->GetSum());
      
      hist->Fit("landau", "q0", "goff", 10, 180);
      TF1 *fit = hist->GetFunction("landau");
      
      Int_t sm = i/30;
      Int_t det = i%30;
      TH2D *h2 = (TH2D*)fHist->At(14+sm);
      Int_t bin = h2->FindBin(det,j);
      // printf("%d %d %d\n", det, j, bin);
      h2->SetBinContent(bin, TMath::Abs(fit->GetParameter(1)));
      h2->SetBinError(bin,fit->GetParError(1));
    }
  }
    
  if (hist) delete hist;
    
  
  TFile *outFile = new TFile(filename, "RECREATE");
  outFile->mkdir("TRD");
  gDirectory->cd("TRD");
  gDirectory->mkdir("RecPoints");
  gDirectory->cd("RecPoints");
  fHist->Write();

  if (fRef) {
    for(Int_t i=0; i<540; i++) {
      //fRefHist[i]->Scale(1./fnEvents);
      fRefHist[i]->Write();
    }
  }

  outFile->Close();

}

//____________________________________________________________________________ 

void AliTRDqaRecPoints::Init()
{
  //
  // Create Reconstructed Points histograms in RecPoints subdir
  //

  //const Int_t kNhist = 14 + 4 * 18 + 2;
  const Int_t kNhist = 88+2*18+1;
  TH1 *hist[kNhist];

  hist[0] = new TH1D("qaTRD_recPoints_det", ";Detector ID of the cluster", 540, -0.5, 539.5);
  hist[1] = new TH2D("qaTRD_recPoints_amp", ";Amplitude", 540, -0.5, 539, 200, -0.5, 199.5);
  hist[2] = new TH1D("qaTRD_recPoints_npad", ";Number of Pads", 12, -0.5, 11.5);

  hist[3] = new TH1D("qaTRD_recPoints_dist2", ";residuals [2pad]", 100, -1, 1);
  hist[4] = new TH1D("qaTRD_recPoints_dist3", ";residuals [3pad]", 100, -1, 1);
  hist[5] = new TH1D("qaTRD_recPoints_dist4", ";residuals [4pad]", 100, -1, 1);
  hist[6] = new TH1D("qaTRD_recPoints_dist5", ";residuals [5pad]", 100, -1, 1);

  hist[7] = new TH2D("qaTRD_recPoints_rowCol", ";row;col", 16, -0.5, 15.5, 145, -0.5, 144.5);
  hist[8] = new TH1D("qaTRD_recPoints_time", ";time bin", 35, -0.5, 34.5);
  hist[9] = new TH1D("qaTRD_recPoints_nCls", ";number of clusters", 500, -0.5, 499.5);

  hist[10] = new TH3D("qaTRD_recPoints_sigTime", ";chamber;time bin;signal", 
		      540, -0.5, 539.5, 35, -0.5, 34.5, 200, -0.5, 199.5);
  hist[11] = new TProfile("qaTRD_recPoints_prf", ";distance;center of gravity"
                         , 120, -0.6, 0.6, -1.2, 1.2, "");

  hist[12] = new TH1D("qaTRD_recPoints_ampMPV", ";amplitude MPV", 150, 0, 150);
  hist[13] = new TH1D("qaTRD_recPoints_ampSigma", ";amplitude Sigma", 200, 0, 200); 
  
  // chamber by chamber
  for(Int_t i=0; i<18; i++) {
    hist[14+i] = new TH2D(Form("qaTRD_recPoints_sigTime_sm%d",i), Form("sm%d;det;time bin"), 
			30, -0.5, 29.5, 35, -0.5, 34.5);
    hist[14+i]->SetMinimum(20);
    hist[14+i]->SetMaximum(40);
  }
 
  // time bin by time bin sm-by-sm
  for(Int_t i=0; i<18; i++) {
    hist[14+18+i] = new TH1D(Form("qaTRD_recPoints_sigTimeShape_sm%d", i), 
			     Form("sm%d;time bin;signal"),
			     35, -0.5, 34.5);

    hist[14+18+i]->SetMaximum(120);    
  }

  // str = 50
  for(Int_t i=0; i<18; i++) {
    hist[50+i] = new TH1D(Form("qaTRD_recPoints_nCls_sm%d",i),
			  Form("sm%d;time bin;number of clusters",i),
			  35, -0.5, 34.5);
  }

  // str = 68
  for(Int_t i=0; i<18; i++) {
    hist[68+i] = new TH1D(Form("qaTRD_recPoints_totalCharge_sm%d", i),
			  Form("sm%d;time bin;total charge", i),
			  35, -0.5, 34.5);
  }

  hist[86] = new TH1D("qaTRD_recPoints_signal", ";amplitude", 200, -0.5, 199.5);
  hist[87] = new TH2D("qaTRD_recPoints_detMap", ";sm;chamber", 18, -0.5, 17.5, 30, -0.5, 29.5);

  for(Int_t i=0; i<18; i++)
    hist[88+i] = new TH2D(Form("qaTRD_recPoints_XY_sm%d", i), 
			  Form("SM%d;Y;X",i), 240, -60, 60, 200, 290, 370);

  for(Int_t i=0; i<18; i++)
    hist[106+i] = new TH2D(Form("qaTRD_recPoints_XPad_sm%d", i), 
			  Form("SM%d;Y;X",i), 144, -0.5, 143.5, 200, 290, 370);

 
  hist[124] = new TH1D("qRef", "ref", 100, 0, 1e4);

  fHist = new TObjArray(200);
  for(Int_t i=0; i<kNhist; i++) {
    //hist[i]->Sumw2();
    fHist->AddAt(hist[i], i);
  }

  // reference histograms
  if (fRef) {
    for(Int_t i=0; i<540; i++) {
      fRefHist[i] = new TH2D(Form("refRecPoints_sm%d", i), "",
			     16, -0.5, 15.5, 144, -0.5, 143.5); //, 30, -0.5, 29.5);
    }
  } else {
    
    TFile *refFile = new TFile("outRef.root");
    refFile->cd("TRD/RecPoints");

    for(Int_t i=0; i<540; i++) {
      fRefHist[i] = (TH2D*)gDirectory->Get(Form("refRecPoints_sm%d", i));
    }
  }
}

//____________________________________________________________________________ 

void AliTRDqaRecPoints::AddEvent(TTree * clustersTree)
{
  //  
  // Makes data from RecPoints
  // 
  
  //  Info("MakeRecPoints", "making");

  Int_t nsize = Int_t(clustersTree->GetTotBytes() / (sizeof(AliTRDcluster))); 
  TObjArray *clusterArray = new TObjArray(nsize+1000); 

  TBranch *branch = clustersTree->GetBranch("TRDcluster");
  if (!branch) {
    AliError("Can't get the branch !");
    return;
  }
  branch->SetAddress(&clusterArray); 

  // Loop through all entries in the tree
  Int_t nEntries   = (Int_t) clustersTree->GetEntries();
  Int_t nbytes     = 0;
  AliTRDcluster *c = 0;
  Int_t nDet[540];
  for (Int_t i=0; i<540; i++) nDet[i] = 0;

  for (Int_t iEntry = 0; iEntry < nEntries; iEntry++) {    

    //printf("Entry = %d\n", iEntry);

    // Import the tree
    nbytes += clustersTree->GetEvent(iEntry);  

    // Get the number of points in the detector
    Int_t nCluster = clusterArray->GetEntriesFast();  

    // Loop through all TRD digits
    for (Int_t iCluster = 0; iCluster < nCluster; iCluster++) { 
      c = (AliTRDcluster *) clusterArray->UncheckedAt(iCluster);

      Int_t iDet = c->GetDetector();
      if (iDet < 0 || iDet > 539) continue;

      
      Int_t iSM = iDet / 30;
      //Int_t iStack = iDet % 30;
      Int_t nPad = c->GetNPads();

      if (fnPad && nPad != fnPad) continue;
      
      //if (iSM == 0 && iStack == 29) continue;
      //if (iSM == 8 && iStack == 11) continue;
      //if (iSM == 8 && iStack == 7)  continue;      
      
      Int_t padRow = c->GetPadRow();
      Int_t padCol = c->GetPadCol();
      //Int_t timeBin = c->GetPadTime();

      Double_t refQ = 0;

      if (fRef) {
	fRefHist[iDet]->Fill(padRow, padCol, c->GetQ());
      } else {
	Int_t bin = fRefHist[iDet]->FindBin(padRow, padCol);
	refQ = fRefHist[iDet]->GetBinContent(bin);
	//printf("bin = %d\n", bin);
      }

      ((TH1D*)fHist->At(124))->Fill(refQ);
      //printf("ref Q = %lf\n", refQ);
      
      Double_t charge = c->GetQ() - (refQ / (490. * 30));
      if (charge < 0) continue;

      if (charge > 20) {
	((TH2D*)fHist->At(88+iSM))->Fill(c->GetY(), c->GetX());
	((TH2D*)fHist->At(106+iSM))->Fill(c->GetPadCol(), c->GetX());
      }

      nDet[iDet]++;
      ((TH1D*)fHist->At(0))->Fill(iDet);
      ((TH1D*)fHist->At(86))->Fill(charge);
      ((TH1D*)fHist->At(1))->Fill(iDet, charge);
      ((TH1D*)fHist->At(2))->Fill(c->GetNPads());
      if (c->GetNPads() < 6)
	((TH1D*)fHist->At(1+c->GetNPads()))->Fill(c->GetCenter());
      
      //if (c->GetPadTime() < 5)
      ((TH2D*)fHist->At(7))->Fill(padRow, c->GetPadCol());
      ((TH1D*)fHist->At(8))->Fill(c->GetPadTime());

      ((TH3D*)fHist->At(10))->Fill(iDet, c->GetPadTime(), charge);
   
      ((TH1D*)fHist->At(50+iSM))->Fill(c->GetPadTime());
      ((TH1D*)fHist->At(68+iSM))->Fill(c->GetPadTime(), charge);

      // PRF for 2pad
      //if (c->GetNPads() == 2) {
      Short_t *sig = c->GetSignals();
      Double_t frac = -10;

      if (sig[0] == 0 && sig[1] == 0 && sig[2] == 0 && sig[5] == 0 && sig[6] == 0) 
	frac = 1. * sig[4] / (sig[3] + sig[4]);

      if (sig[0] == 0 && sig[1] == 0 && sig[4] == 0 && sig[5] == 0 && sig[6] == 0)
	frac = -1. * sig[2] / (sig[2] + sig[3]);

      if (frac > -10)  ((TProfile*)fHist->At(11))->Fill(c->GetCenter(), frac);
	
      //}
    }
  }

  for(Int_t i=0; i<540; i++) 
    if (nDet[i] > 0) ((TH1D*)fHist->At(9))->Fill(nDet[i]);

  delete clusterArray;

  /*
  TFile *outFile = new TFile("outQA.root", "RECREATE");
  outFile->mkdir("TRD");
  gDirectory->cd("TRD");
  gDirectory->mkdir("RecPoints");
  gDirectory->cd("RecPoints");
  fHist->Write();
  outFile->Close();
  */
}

//____________________________________________________________________________ 
