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

/*
$Log$
Revision 1.5.4.1  2002/06/10 15:26:12  hristov
Merged with v3-08-02

Revision 1.9  2002/05/13 09:53:08  hristov
Some frequent problems corrected: arrays with variable size have to be defined via the operator new, default values for the arguments have to be  used only in the header files, etc.

Revision 1.8  2002/05/08 18:21:40  kowal2
Now the code is blind to the binning used for pulls, efficiencies etc.

Revision 1.7  2002/04/10 16:30:07  kowal2
logs added

*/


/**************************************************************************
 *                                                                        *
 * This class builds AliTPCtrack objects from generated tracks to feed    *
 * ITS tracking (V2). The AliTPCtrack is built from its first hit in      *
 * the TPC. The track is assigned a Kalman-like covariance matrix         *
 * depending on its pT and pseudorapidity and track parameters are        *
 * smeared according to this covariance matrix.                           *
 * Output file contains sorted tracks, ready for matching with ITS.       *
 *                                                                        *
 * For details:                                                           *
 * http://www.pd.infn.it/alipd/talks/soft/adIII02/TPCtrackingParam.htm    *
 *                                                                        *
 * Test macro is: AliBarrelRec_TPCparam.C                                 *   
 *                                                                        *
 * 2002/04/28: Major upgrade of the class                                 *
 * - Covariance matrices and pulls are now separeted for pions, kaons     *
 *   and electrons.                                                       *
 * - A parameterization of dE/dx in the TPC has been included and it is   *
 *   used to assign a mass to each track according to a rough dE/dx PID.  *
 * - All the "numbers" have been moved to the file with the DataBase,     *
 *   they are read as objects of the class AliTPCkineGrid, and assigned   *
 *   to data memebers of the class AliTPCtrackerParam.                    *
 * - All the code necessary to create a BataBase has been included in     *
 *   class (see the macro AliTPCtrackingParamDB.C for the details).       *
 *                                                                        *
 *  Origin: Andrea Dainese, Padova - e-mail: andrea.dainese@pd.infn.it    *
 *                                                                        *
 **************************************************************************/
#include <TChain.h>
#include <TF1.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMatrixD.h>
#include <TNtuple.h>
#include <TSystem.h>
#include "alles.h"
#include "AliGausCorr.h"
#include "AliKalmanTrack.h"
#include "AliMagF.h"
#include "AliMagFCM.h"
#include "AliTPCkineGrid.h"
#include "AliTPCtrack.h"
#include "AliTPCtrackerParam.h"

Double_t RegFunc(Double_t *x,Double_t *par) {
// This is the function used to regularize the covariance matrix
  Double_t value = par[0]+par[1]/TMath::Power(x[0],par[2])/
                   TMath::Power(x[1],par[3]);
  return value;
}
Double_t FitRegFunc(Double_t *x,Double_t *par) {
// This is the function used for the fit of the covariance 
// matrix as a function of the momentum. 
// For cosl the average value 0.87 is used.
  Double_t value = par[0]+par[1]/TMath::Power(x[0],par[2])/
                   TMath::Power(0.87,par[3]);
  return value;
}

// structure for DB building
typedef struct {
  Int_t    pdg;
  Int_t    bin;
  Double_t r;
  Double_t p;
  Double_t pt;
  Double_t cosl;
  Double_t eta;
  Double_t dpt;
  Double_t dP0,dP1,dP2,dP3,dP4;
  Double_t c00,c10,c11,c20,c21,c22,c30,c31,c32,c33,c40,c41,c42,c43,c44;
  Double_t dEdx;
} COMPTRACK;
// cov matrix structure 
typedef struct {
  Double_t c00,c10,c11,c20,c21,c22,c30,c31,c32,c33,c40,c41,c42,c43,c44;
} COVMATRIX;

ClassImp(AliTPCtrackerParam)

//-----------------------------------------------------------------------------
AliTPCtrackerParam::AliTPCtrackerParam(const Int_t coll,const Double_t Bz)
{
//-----------------------------------------------------------------------------
// This is the class conctructor 
//-----------------------------------------------------------------------------

  fBz = Bz;             // value of the z component of L3 field (Tesla)
  fColl = coll;         // collision code (0: PbPb6000)
  fSelAndSmear = kTRUE; // by default selection and smearing are done

  if(fBz!=0.4) {
    cerr<<"AliTPCtrackerParam::AliTPCtrackerParam:  Invalid field!\n";
    cerr<<"      Available:  0.4"<<endl;
  }
  if(fColl!=0) { 
    cerr<<"AliTPCtrackerParam::AliTPCtrackerParam:  Invalid collision!\n";
    cerr<<"      Available:  0   ->   PbPb6000"<<endl; 
  }

  fDBfileName = gSystem->Getenv("ALICE_ROOT");  
  fDBfileName.Append("/TPC/CovMatrixDB_");
  if(fColl==0) fDBfileName.Append("PbPb6000");
  if(fBz==0.4) fDBfileName.Append("_B0.4T.root");
}
//-----------------------------------------------------------------------------
AliTPCtrackerParam::~AliTPCtrackerParam() 
{}
//-----------------------------------------------------------------------------

Int_t AliTPCtrackerParam::BuildTPCtracks(const TFile *inp, TFile *out, Int_t n)
{
//-----------------------------------------------------------------------------
// This function creates the TPC parameterized tracks
//-----------------------------------------------------------------------------

  TFile *fileDB=0;
  TTree *covTreePi[50];
  TTree *covTreeKa[50];
  TTree *covTreeEl[50];

  if(fSelAndSmear) {
    cerr<<"+++\n+++ Reading DataBase from:\n+++ "<<
      fDBfileName.Data()<<"\n+++\n"; 
    // Read paramters from DB file
    if(!ReadAllData(fDBfileName.Data())) {
      cerr<<"AliTPCtrackerParam::BuildTPCtracks: \
             Could not read data from DB\n\n"; return 1; 
    }
    // Read the trees with regularized cov. matrices from DB
    TString str;
    fileDB = TFile::Open(fDBfileName.Data());
    Int_t nBinsPi = fDBgridPi.GetTotBins();
    for(Int_t l=0; l<nBinsPi; l++) {
      str = "/CovMatrices/Pions/CovTreePi_bin";
      str+=l;
      covTreePi[l] = (TTree*)fileDB->Get(str.Data());
    }
    Int_t nBinsKa = fDBgridKa.GetTotBins();
    for(Int_t l=0; l<nBinsKa; l++) {
      str = "/CovMatrices/Kaons/CovTreeKa_bin";
      str+=l;
      covTreeKa[l] = (TTree*)fileDB->Get(str.Data());
    }
    Int_t nBinsEl = fDBgridEl.GetTotBins();
    for(Int_t l=0; l<nBinsEl; l++) {
      str = "/CovMatrices/Electrons/CovTreeEl_bin";
      str+=l;
      covTreeEl[l] = (TTree*)fileDB->Get(str.Data());
    }

  } else cerr<<"\n ! Creating ALL TRUE tracks at TPC 1st hit !\n\n";

  TFile *infile=(TFile*)inp;
  infile->cd();

  // Get gAlice object from file
  if(!(gAlice=(AliRun*)infile->Get("gAlice"))) {
    cerr<<"gAlice has not been found on galice.root !\n";
    return 1;
  }

  // Check if value in the galice file is equal to selected one (fBz)
  AliMagFCM *fiel = (AliMagFCM*)gAlice->Field();
  Double_t fieval=(Double_t)fiel->SolenoidField()/10.;
  printf("Magnetic field is %6.2f Tesla\n",fieval);
  if(fBz!=fieval) {
    cerr<<"AliTPCtrackerParam::BuildTPCtracks:  Invalid field!"<<endl;
    cerr<<"Field selected is: "<<fBz<<" T\n";
    cerr<<"Field found on file is: "<<fieval<<" T\n";
    return 0;
  }

  // Get TPC detector 
  AliTPC *TPC=(AliTPC*)gAlice->GetDetector("TPC");
  Int_t ver = TPC->IsVersion(); 
  cerr<<"+++ TPC version "<<ver<<" has been found !\n";
  AliTPCParam *digp=(AliTPCParam*)infile->Get("75x40_100x60");
  if(digp){
    delete digp;
    digp = new AliTPCParamSR();
  }
  else digp=(AliTPCParam*)infile->Get("75x40_100x60_150x60");
  
  if(!digp) { cerr<<"TPC parameters have not been found !\n"; return 1; }
  TPC->SetParam(digp);

  // Set the conversion constant between curvature and Pt
  AliKalmanTrack::SetConvConst(100/0.299792458/fBz);

  TParticle   *Part=0;
  AliTPChit   *tpcHit=0;
  AliTPCtrack *tpctrack=0;
  Double_t     hPx,hPy,hPz,hPt,hEta,xg,yg,zg,xl,yl,zl;
  Double_t     alpha;
  Float_t      cosAlpha,sinAlpha;
  Int_t        bin,label,pdg,charge;
  Int_t        tracks=0;
  Int_t        nParticles,nTracks,arrentr;
  Char_t       tname[100];
  //Int_t nSel=0,nAcc=0;

  // loop over first n events in file
  for(Int_t evt=0; evt<n; evt++){
    cerr<<"+++\n+++ Processing event "<<evt<<"\n+++\n";

    // tree for TPC tracks
    sprintf(tname,"TreeT_TPC_%d",evt);
    TTree *tracktree = new TTree(tname,"Tree with TPC tracks");
    tracktree->Branch("tracks","AliTPCtrack",&tpctrack,20000,0);

    // array for TPC tracks
    TObjArray tarray(20000);

    // get the particles stack
    nParticles = gAlice->GetEvent(evt);
    Bool_t * done = new Bool_t[nParticles];
    Int_t  * pdgCodes = new Int_t[nParticles];

    // loop on particles and store pdg codes
    for(Int_t l=0; l<nParticles; l++) {
      Part        = gAlice->Particle(l);
      pdgCodes[l] = Part->GetPdgCode();
      done[l]     = kFALSE;
    }

    // Get TreeH with hits
    TTree *TH=gAlice->TreeH(); 
    nTracks=(Int_t)TH->GetEntries();
    cerr<<"+++\n+++ Number of particles in event "<<evt<<":  "<<nParticles<<
      "\n+++\n+++ Number of \"primary tracks\"(entries in TreeH): "<<nTracks<<
      "\n+++\n\n";

    // loop over entries in TreeH
    for(Int_t l=0; l<nTracks; l++) {
      if(l%1000==0) cerr<<"  --- Processing primary track "
			<<l<<" of "<<nTracks<<" ---\r";
      TPC->ResetHits();
      TH->GetEvent(l);
      // Get FirstHit
      tpcHit=(AliTPChit*)TPC->FirstHit(-1);
      for( ; tpcHit; tpcHit=(AliTPChit*)TPC->NextHit() ) {
        if(tpcHit->fQ !=0.) continue;
        // Get particle momentum at hit
        hPx=tpcHit->X(); hPy=tpcHit->Y(); hPz=tpcHit->Z();
        hPt=TMath::Sqrt(hPx*hPx+hPy*hPy);
        // reject hits with Pt<mag*0.45 GeV/c
        if(hPt<(fBz*0.45)) continue;

        // Get track label
        label=tpcHit->Track();
        // check if this track has already been processed
        if(done[label]) continue;
        // PDG code & electric charge
        pdg = pdgCodes[label];
        if(pdg>200 || pdg==-11 || pdg==-13) { charge=1; }
        else if(pdg<-200 || pdg==11 || pdg==13) { charge=-1; }
        else continue;
        pdg = TMath::Abs(pdg);
        if(pdg>3000) pdg=211;
        if(fSelAndSmear) SetParticle(pdg);

        if((tpcHit=(AliTPChit*)TPC->NextHit())==0) break;
        if(tpcHit->fQ != 0.) continue;
        // Get global coordinates of hit
        xg=tpcHit->X(); yg=tpcHit->Y(); zg=tpcHit->Z();
        if(TMath::Sqrt(xg*xg+yg*yg)>90.) continue;

        // Get TPC sector, Alpha angle and local coordinates
        // cerr<<"Sector "<<tpcHit->fSector<<endl;
        digp->AdjustCosSin(tpcHit->fSector,cosAlpha,sinAlpha);
        alpha = TMath::ATan2(sinAlpha,cosAlpha);
        xl = xg*cosAlpha + yg*sinAlpha;
        yl =-xg*sinAlpha + yg*cosAlpha;
        zl = zg;
        //printf("Alpha %f   xl %f  yl %f  zl %f\n",Alpha,xl,yl,zl);

        // reject tracks which are not in the TPC acceptance
        if(TMath::Abs(zl+(244.-xl)*hPz/hPt)>252.) continue;

        hEta = -TMath::Log(TMath::Tan(0.25*TMath::Pi()-0.5*TMath::ATan(hPz/hPt)));  

        // Apply selection according to TPC efficiency
        //if(TMath::Abs(pdg)==211) nAcc++;
        if(fSelAndSmear && !SelectedTrack(hPt,hEta)) continue; 
        //if(TMath::Abs(pdg)==211) nSel++;

        // create AliTPCtrack object
        BuildTrack(alpha,xl,yl,zl,hPx,hPy,hPz,hPt,charge);
	
        if(fSelAndSmear) {
          bin = fDBgrid->GetBin(hPt,hEta);
          switch (pdg) {
          case 211:
            fCovTree = covTreePi[bin];
            break;
          case 321:
            fCovTree = covTreeKa[bin];
            break;
          case 2212:
            fCovTree = covTreePi[bin];
            break;
          case 11:
            fCovTree = covTreeEl[bin];
            break;
          case 13:
            fCovTree = covTreePi[bin];
            break;
          }
          // deal with covariance matrix and smearing of parameters
          CookTrack(hPt,hEta);

          // assign the track a dE/dx and make a rough PID
          CookdEdx(hPt,hEta);
        }

        // put track in array
        AliTPCtrack *iotrack = new AliTPCtrack(fTrack);
        iotrack->SetLabel(label);
        tarray.AddLast(iotrack);
        // Mark track as "done" and register the pdg code
        done[label]=kTRUE; 
        tracks++;
      }
 
    } // loop over entries in TreeH

    // sort array with TPC tracks (decreasing pT)
    tarray.Sort();

    arrentr = tarray.GetEntriesFast();
    for(Int_t l=0; l<arrentr; l++) {
      tpctrack=(AliTPCtrack*)tarray.UncheckedAt(l);
      tracktree->Fill();
    }

    // write the tree with tracks in the output file
    out->cd();
    tracktree->Write();

    delete tracktree;
    delete [] done;
    delete [] pdgCodes;

    printf("\n\n+++\n+++ Number of TPC tracks: %d\n+++\n",tracks);
    //cerr<<"Average Eff: "<<(Float_t)nSel/nAcc<<endl;

  } // loop on events

  if(fileDB) fileDB->Close();

  return 0;
}
//-----------------------------------------------------------------------------
void AliTPCtrackerParam::AnalyzedEdx(const Char_t *outName,Int_t pdg) {
//-----------------------------------------------------------------------------
// This function computes the dE/dx for pions, kaons, protons and electrons,
// in the [pT,eta] bins.
// Input file is CovMatrix_AllEvts.root.  
//-----------------------------------------------------------------------------

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(10001);

  Char_t *part="PIONS";
  Double_t ymax=500.;

  // create a chain with compared tracks
  TChain cmptrkchain("cmptrktree");
  cmptrkchain.Add("CovMatrix_AllEvts_1.root");
  cmptrkchain.Add("CovMatrix_AllEvts_2.root");
  cmptrkchain.Add("CovMatrix_AllEvts_3.root");

  COMPTRACK cmptrk; 
  cmptrkchain.SetBranchAddress("comptracks",&cmptrk);
  Int_t entries = (Int_t)cmptrkchain.GetEntries(); 
  cerr<<" Number of entries: "<<entries<<endl;

  InitializeKineGrid("DB","points");
  InitializeKineGrid("dEdx","bins");

  switch(pdg) {
  case 211:
    part = "PIONS";
    ymax = 100.;
    break;
  case 321:
    part = "KAONS";
    ymax = 300.;
    break;
  case 2212:
    part = "PROTONS";
    ymax = 500.;
    break;
  case 11:
    part = "ELECTRONS";
    ymax = 100.;
    break;
  }

  SetParticle(pdg);

  const Int_t nTotBins = fDBgrid->GetTotBins(); 

  cerr<<" Fit bins: "<<nTotBins<<endl;

  Int_t bin=0;
  Int_t * n = new Int_t[nTotBins];
  Double_t * p = new Double_t[nTotBins];
  Double_t * ep = new Double_t[nTotBins];
  Double_t * mean = new Double_t[nTotBins];
  Double_t * sigma = new Double_t[nTotBins];

  for(Int_t l=0; l<nTotBins; l++) {
    n[l] = 1; // set to 1 to avoid divisions by 0
    p[l] = mean[l] = sigma[l] = ep[l] = 0.; 
  }

  // loop on chain entries for the mean
  for(Int_t l=0; l<entries; l++) {
    cmptrkchain.GetEvent(l);
    if(TMath::Abs(cmptrk.pdg)!=pdg) continue;
    bin = fDBgrid->GetBin(cmptrk.pt,cmptrk.eta);
    p[bin] += cmptrk.p;
    mean[bin] += cmptrk.dEdx;
    n[bin]++;
  } // loop on chain entries

  for(Int_t l=0; l<nTotBins; l++) {
    p[l] /= n[l];
    mean[l] /= n[l];
    n[l] = 1; // set to 1 to avoid divisions by 0
  }

  // loop on chain entries for the sigma
  for(Int_t l=0; l<entries; l++) {
    cmptrkchain.GetEvent(l);
    if(TMath::Abs(cmptrk.pdg)!=pdg) continue;
    bin = fDBgrid->GetBin(cmptrk.pt,cmptrk.eta);
    if(cmptrk.p<1. && TMath::Abs(cmptrk.p-p[bin])>0.025) continue;
    n[bin]++;
    sigma[bin] += (cmptrk.dEdx-mean[bin])*(cmptrk.dEdx-mean[bin]);
  } // loop on chain entries
  
  for(Int_t l=0; l<nTotBins; l++) {
    sigma[l] = TMath::Sqrt(sigma[l]/n[l]);
  }

  // create the canvas
  TCanvas *canv = new TCanvas("canv","dEdx",0,0,900,700); 

  // create the graph for dEdx vs p
  TGraphErrors *gr = new TGraphErrors(nTotBins,p,mean,ep,sigma);
  TString title("  : dE/dx vs momentum"); title.Prepend(part);
  TH2F *frame = new TH2F("frame1",title.Data(),2,0.1,50,2,0,ymax);
  frame->SetXTitle("p [GeV/c]");
  frame->SetYTitle("dE/dx [a.u.]");
  canv->SetLogx();
  frame->Draw();
  gr->Draw("P");

  switch(pdg) {
  case 211:
    for(Int_t i=0; i<nTotBins; i++) {
      fdEdxMeanPi.SetParam(i,mean[i]);
      fdEdxRMSPi.SetParam(i,sigma[i]);
    }    
    break;
  case 321:
    for(Int_t i=0; i<nTotBins; i++) {
      fdEdxMeanKa.SetParam(i,mean[i]);
      fdEdxRMSKa.SetParam(i,sigma[i]);
    }    
    break;
  case 2212:
    for(Int_t i=0; i<nTotBins; i++) {
      fdEdxMeanPr.SetParam(i,mean[i]);
      fdEdxRMSPr.SetParam(i,sigma[i]);
    }    
    break;
  case 11:
    for(Int_t i=0; i<nTotBins; i++) {
      fdEdxMeanEl.SetParam(i,mean[i]);
      fdEdxRMSEl.SetParam(i,sigma[i]);
    }    
    break;
  }

  // write results to file
  WritedEdx(outName,pdg);

  delete [] n;
  delete [] p;
  delete [] ep;
  delete [] mean;
  delete [] sigma;
  
  return;
}
//-----------------------------------------------------------------------------
void AliTPCtrackerParam::AnalyzePulls(const Char_t *outName) {
//-----------------------------------------------------------------------------
// This function computes the pulls for pions, kaons and electrons,
// in the [pT,eta] bins.
// Input file is CovMatrix_AllEvts.root.
// Output file is pulls.root.  
//-----------------------------------------------------------------------------

  // create a chain with compared tracks
  TChain cmptrkchain("cmptrktree");
  cmptrkchain.Add("CovMatrix_AllEvts_1.root");
  cmptrkchain.Add("CovMatrix_AllEvts_2.root");
  cmptrkchain.Add("CovMatrix_AllEvts_3.root");

  COMPTRACK cmptrk; 
  cmptrkchain.SetBranchAddress("comptracks",&cmptrk);
  Int_t entries = (Int_t)cmptrkchain.GetEntries(); 
  cerr<<" Number of entries: "<<entries<<endl;

  Int_t thisPdg=0;
  Char_t hname[100], htitle[100];
  Int_t nTotBins,bin;

  AliTPCkineGrid  pulls[5];
  TH1F *hDum = new TH1F("name","title",100,-7.,7.);
  TF1 *g = new TF1("g","gaus");

  InitializeKineGrid("pulls","bins");
  InitializeKineGrid("DB","points");



  // loop on the particles Pi,Ka,El
  for(Int_t part=0; part<3; part++) {

    switch (part) {
    case 0: // pions
      thisPdg=211; 
      cerr<<" Processing pions ...\n";
      break;   
    case 1: // kaons
      thisPdg=321; 
      cerr<<" Processing kaons ...\n";
      break;
      
    case 2: // electrons
      thisPdg=11; 
      cerr<<" Processing electrons ...\n";
      break;
    }

    SetParticle(thisPdg);

    for(Int_t i=0;i<5;i++) {
      pulls[i].~AliTPCkineGrid(); 
      new(&pulls[i]) AliTPCkineGrid(*(fPulls+i));
    }
    nTotBins = fDBgrid->GetTotBins();
 
    // create histograms for the all the bins
    TH1F *hPulls0_=NULL;
    TH1F *hPulls1_=NULL;
    TH1F *hPulls2_=NULL;
    TH1F *hPulls3_=NULL;
    TH1F *hPulls4_=NULL;

    hPulls0_ = new TH1F[nTotBins]; 
    hPulls1_ = new TH1F[nTotBins]; 
    hPulls2_ = new TH1F[nTotBins]; 
    hPulls3_ = new TH1F[nTotBins]; 
    hPulls4_ = new TH1F[nTotBins]; 


    for(Int_t i=0; i<nTotBins; i++) {
      sprintf(hname,"hPulls0_%d",i);
      sprintf(htitle,"P0 pulls for bin %d",i);
      hDum->SetName(hname); hDum->SetTitle(htitle);
      hPulls0_[i] = *hDum;
      sprintf(hname,"hPulls1_%d",i);
      sprintf(htitle,"P1 pulls for bin %d",i);
      hDum->SetName(hname); hDum->SetTitle(htitle);
      hPulls1_[i] = *hDum;
      sprintf(hname,"hPulls2_%d",i);
      sprintf(htitle,"P2 pulls for bin %d",i);
      hDum->SetName(hname); hDum->SetTitle(htitle);
      hPulls2_[i] = *hDum;
      sprintf(hname,"hPulls3_%d",i);
      sprintf(htitle,"P3 pulls for bin %d",i);
      hDum->SetName(hname); hDum->SetTitle(htitle);
      hPulls3_[i] = *hDum;
      sprintf(hname,"hPulls4_%d",i);
      sprintf(htitle,"P4 pulls for bin %d",i);
      hDum->SetName(hname); hDum->SetTitle(htitle);
      hPulls4_[i] = *hDum;
    }

    // loop on chain entries 
    for(Int_t i=0; i<entries; i++) {
      cmptrkchain.GetEvent(i);
      if(TMath::Abs(cmptrk.pdg)!=thisPdg) continue;
      // fill histograms with the pulls
      bin = fDBgrid->GetBin(cmptrk.pt,cmptrk.eta); 
      hPulls0_[bin].Fill(cmptrk.dP0/TMath::Sqrt(cmptrk.c00));
      hPulls1_[bin].Fill(cmptrk.dP1/TMath::Sqrt(cmptrk.c11));
      hPulls2_[bin].Fill(cmptrk.dP2/TMath::Sqrt(cmptrk.c22));
      hPulls3_[bin].Fill(cmptrk.dP3/TMath::Sqrt(cmptrk.c33));
      hPulls4_[bin].Fill(cmptrk.dP4/TMath::Sqrt(cmptrk.c44));
    } // loop on chain entries

    // compute the sigma of the distributions
    for(Int_t i=0; i<nTotBins; i++) {
      if(hPulls0_[i].GetEntries()>1000) {
	g->SetRange(-3.*hPulls0_[i].GetRMS(),3.*hPulls0_[i].GetRMS());
	hPulls0_[i].Fit("g","R,Q,N");
	pulls[0].SetParam(i,g->GetParameter(2));
      } else pulls[0].SetParam(i,-1.);
      if(hPulls1_[i].GetEntries()>1000) {
	g->SetRange(-3.*hPulls1_[i].GetRMS(),3.*hPulls1_[i].GetRMS());
	hPulls1_[i].Fit("g","R,Q,N");
	pulls[1].SetParam(i,g->GetParameter(2));
      } else pulls[1].SetParam(i,-1.);
      if(hPulls2_[i].GetEntries()>1000) {
	g->SetRange(-3.*hPulls2_[i].GetRMS(),3.*hPulls2_[i].GetRMS());
	hPulls2_[i].Fit("g","R,Q,N");
	pulls[2].SetParam(i,g->GetParameter(2));
      } else pulls[2].SetParam(i,-1.);
      if(hPulls3_[i].GetEntries()>1000) {
	g->SetRange(-3.*hPulls3_[i].GetRMS(),3.*hPulls3_[i].GetRMS());
	hPulls3_[i].Fit("g","R,Q,N");
	pulls[3].SetParam(i,g->GetParameter(2));
      } else pulls[3].SetParam(i,-1.);
      if(hPulls4_[i].GetEntries()>1000) {
	g->SetRange(-3.*hPulls4_[i].GetRMS(),3.*hPulls4_[i].GetRMS());
	hPulls4_[i].Fit("g","R,Q,N");
	pulls[4].SetParam(i,g->GetParameter(2));
      } else pulls[4].SetParam(i,-1.);
    } // loop on bins


    switch (part) {
    case 0: // pions
      for(Int_t i=0;i<5;i++) {
	fPullsPi[i].~AliTPCkineGrid();
	new(&fPullsPi[i]) AliTPCkineGrid(pulls[i]);
      }
      break;
    case 1: // kaons
      for(Int_t i=0;i<5;i++) {
	fPullsKa[i].~AliTPCkineGrid();
	new(&fPullsKa[i]) AliTPCkineGrid(pulls[i]);
      }
      break;
    case 2: // electrons
      for(Int_t i=0;i<5;i++) {
	fPullsEl[i].~AliTPCkineGrid();
	new(&fPullsEl[i]) AliTPCkineGrid(pulls[i]);
      }
      break;
    }

    delete [] hPulls0_;
    delete [] hPulls1_;
    delete [] hPulls2_;
    delete [] hPulls3_;
    delete [] hPulls4_;
    
  } // loop on particle species

  // write pulls to file
  WritePulls(outName);


  return;
}
//-----------------------------------------------------------------------------
void AliTPCtrackerParam::BuildTrack(Double_t alpha,
		                    Double_t x,Double_t y,Double_t z,
		                    Double_t px,Double_t py,
		                    Double_t pz,Double_t pt,
		                    Int_t ch) {  
//-----------------------------------------------------------------------------
// This function uses GEANT info to set true track parameters
//-----------------------------------------------------------------------------
  Double_t xref = x;
  Double_t xx[5],cc[15];
  cc[0]=cc[2]=cc[5]=cc[9]=cc[14]=10.;
  cc[1]=cc[3]=cc[4]=cc[6]=cc[7]=cc[8]=cc[10]=cc[11]=cc[12]=cc[13]=0.;
  
  // Magnetic field
  TVector3 bfield(0.,0.,fBz);
  
  
  // radius [cm] of track projection in (x,y) 
  Double_t rho = pt*100./0.299792458/bfield.Z();
  // center of track projection in local reference frame
  TVector3 hmom,hpos;


  // position (local) and momentum (local) at the hit
  // in the bending plane (z=0)
  hpos.SetXYZ(x,y,0.);
  hmom.SetXYZ(px*TMath::Cos(alpha)+py*TMath::Sin(alpha),-px*TMath::Sin(alpha)+py*TMath::Cos(alpha),0.);
  TVector3 vrho = hmom.Cross(bfield);
  vrho *= ch;
  vrho.SetMag(rho);

  TVector3 vcenter = hpos+vrho;

  Double_t x0 = vcenter.X();

  // fX     = xref         X-coordinate of this track (reference plane)
  // fAlpha = Alpha        Rotation angle the local (TPC sector) 
  // fP0    = YL           Y-coordinate of a track
  // fP1    = ZG           Z-coordinate of a track
  // fP2    = C*x0         x0 is center x in rotated frame
  // fP3    = Tgl          tangent of the track momentum dip angle
  // fP4    = C            track curvature
  xx[0] = y;
  xx[1] = z;
  xx[3] = pz/pt;
  xx[4] = -ch/rho;
  xx[2] = xx[4]*x0;

  // create the object AliTPCtrack    
  AliTPCtrack track(0,xx,cc,xref,alpha);
  new(&fTrack) AliTPCtrack(track);

  return;
}
//-----------------------------------------------------------------------------
void AliTPCtrackerParam::CompareTPCtracks(
			   const Char_t* galiceName,
			   const Char_t* trkGeaName,
			   const Char_t* trkKalName,
			   const Char_t* covmatName,
			   const Char_t* tpceffName) const {
//-----------------------------------------------------------------------------
// This function compares tracks from TPC Kalman Filter V2 with 
// true tracks at TPC 1st hit. It gives:
//   - a tree with Kalman cov. matrix elements, resolutions, dEdx
//   - the efficiencies as a function of particle type, pT, eta
//-----------------------------------------------------------------------------

  TFile *kalFile    = TFile::Open(trkKalName);
  TFile *geaFile    = TFile::Open(trkGeaName);
  TFile *galiceFile = TFile::Open(galiceName);

  // particles from TreeK
  AliRun *gAlice = (AliRun*)galiceFile->Get("gAlice");
  const Int_t nparticles = gAlice->GetEvent(0);

  Int_t * kalLab = new Int_t[nparticles];
  for(Int_t i=0; i<nparticles; i++) kalLab[i] = -1; 
 

  // tracks from Kalman
  TTree *kaltree=(TTree*)kalFile->Get("TreeT_TPC_0");
  AliTPCtrack *kaltrack=new AliTPCtrack; 
  kaltree->SetBranchAddress("tracks",&kaltrack);
  Int_t kalEntries = (Int_t)kaltree->GetEntries();

  // tracks from 1st hit
  TTree *geatree=(TTree*)geaFile->Get("TreeT_TPC_0");
  AliTPCtrack *geatrack=new AliTPCtrack; 
  geatree->SetBranchAddress("tracks",&geatrack);
  Int_t geaEntries = (Int_t)geatree->GetEntries();

  cerr<<"+++\n+++ Number of tracks:  TPC Kalman  = "<<kalEntries<<endl<<"+++                    TPC 1st hit = "<<geaEntries<<endl<<"+++\n";

  // set pointers for TPC tracks: 
  // the entry number of the track labelled l is stored in kalLab[l]
  Int_t fake=0,mult=0;
  for (Int_t j=0; j<kalEntries; j++) {
    kaltree->GetEvent(j);
    if(kaltrack->GetLabel()>=0) { 
      if(kalLab[kaltrack->GetLabel()]!=-1) mult++; 
      kalLab[kaltrack->GetLabel()] = j;
    }
    else fake++;
  }
  
  cerr<<"+++ Number of fake tracks in TPC Kalman: "<<fake<<endl;
  cerr<<"+++ Number of multiply found tracks in TPC Kalman: "<<mult<<endl;

  TParticle *Part;
  Int_t    label;
  Double_t cc[15],dAlpha;
  //                 ---  [Pt,eta] binning ---
  //
  //               |eta|<0.3    0.3<=|eta|<0.6    0.6<=|eta|
  //      Pt<0.3       0              1                2 
  // 0.3<=Pt<0.5       3              4                5  
  // 0.5<=Pt<1.0       6              7                8 
  // 1.0<=Pt<1.5       9             10               11 
  // 1.5<=Pt<3.0      12             13               14 
  // 3.0<=Pt<5.0      15             16               17 
  // 5.0<=Pt<7.0      18             19               20 
  // 7.0<=Pt<15.      21             22               23 
  // 15.<=Pt          24             25               26
  // 
  Int_t    pi=0,ka=0,mu=0,el=0,pr=0;
  Int_t    geaPi[27],geaKa[27],geaPr[27],geaEl[27],geaMu[27];
  Int_t    kalPi[27],kalKa[27],kalPr[27],kalEl[27],kalMu[27];
  Float_t  effPi[27],effKa[27],effPr[27],effEl[27],effMu[27];
  Int_t bin;
  for(Int_t j=0; j<27; j++) {
    geaPi[j]=geaKa[j]=geaPr[j]=geaEl[j]=geaMu[j]=0;
    kalPi[j]=kalKa[j]=kalPr[j]=kalEl[j]=kalMu[j]=0;
    effPi[j]=effKa[j]=effPr[j]=effEl[j]=effMu[j]=-1.;
  }

  COMPTRACK cmptrk;
  // create the tree for comparison results
  TTree *cmptrktree = new TTree("cmptrktree","results of track comparison");
  cmptrktree->Branch("comptracks",&cmptrk,"pdg/I:bin:r/D:p:pt:cosl:eta:dpt:dP0:dP1:dP2:dP3:dP4:c00:c10:c11:c20:c21:c22:c30:c31:c32:c33:c40:c41:c42:c43:c44:dEdx");
  
  cerr<<"Doing track comparison...\n";
  // loop on tracks at 1st hit
  for(Int_t j=0; j<geaEntries; j++) {
    geatree->GetEvent(j);
    
    label = geatrack->GetLabel();
    Part = (TParticle*)gAlice->Particle(label);
    cmptrk.pdg = Part->GetPdgCode();
    cmptrk.eta = Part->Eta();
    cmptrk.r = TMath::Sqrt(Part->Vx()*Part->Vx()+Part->Vy()*Part->Vy());

    cmptrk.pt   = 1/TMath::Abs(geatrack->Get1Pt());
    cmptrk.cosl = TMath::Cos(TMath::ATan(geatrack->GetTgl()));
    cmptrk.p    = cmptrk.pt/cmptrk.cosl;


    bin = GetBin(cmptrk.pt,cmptrk.eta);
    cmptrk.bin = bin;
    if(abs(cmptrk.pdg)==211)  geaPi[bin]++;  
    if(abs(cmptrk.pdg)==321)  geaKa[bin]++;   
    if(abs(cmptrk.pdg)==2212) geaPr[bin]++;   
    if(abs(cmptrk.pdg)==11)   geaEl[bin]++;  
    if(abs(cmptrk.pdg)==13)   geaMu[bin]++; 


    // check if there is the corresponding track in the TPC kalman and get it
    if(kalLab[label]==-1) continue;
    kaltree->GetEvent(kalLab[label]);

    // and go on only if it has xref = 84.57 cm (inner pad row)
    if(kaltrack->GetX()>90.) continue;

    if(abs(cmptrk.pdg)==211)  { kalPi[bin]++; pi++; }
    if(abs(cmptrk.pdg)==321)  { kalKa[bin]++; ka++; }
    if(abs(cmptrk.pdg)==2212) { kalPr[bin]++; pr++; }
    if(abs(cmptrk.pdg)==11)   { kalEl[bin]++; el++; }
    if(abs(cmptrk.pdg)==13)   { kalMu[bin]++; mu++; }

    kaltrack->PropagateTo(geatrack->GetX());

    cmptrk.dEdx = kaltrack->GetdEdx();

    // compute errors on parameters
    dAlpha = kaltrack->GetAlpha()-geatrack->GetAlpha();
    if(TMath::Abs(dAlpha)>0.1) { cerr<<" ! WRONG SECTOR !\n"; continue; }

    cmptrk.dP0 = kaltrack->GetY()-geatrack->GetY();
    cmptrk.dP1 = kaltrack->GetZ()-geatrack->GetZ();
    cmptrk.dP2 = kaltrack->GetEta()-geatrack->GetEta();
    cmptrk.dP3 = kaltrack->GetTgl()-geatrack->GetTgl();
    cmptrk.dP4 = kaltrack->GetC()-geatrack->GetC();
    cmptrk.dpt = 1/kaltrack->Get1Pt()-1/geatrack->Get1Pt();

    // get covariance matrix
    // beware: lines 3 and 4 are inverted!
    kaltrack->GetCovariance(cc);

    cmptrk.c00 = cc[0];
    cmptrk.c10 = cc[1];
    cmptrk.c11 = cc[2];
    cmptrk.c20 = cc[3];
    cmptrk.c21 = cc[4];
    cmptrk.c22 = cc[5];
    cmptrk.c30 = cc[10];
    cmptrk.c31 = cc[11];
    cmptrk.c32 = cc[12];
    cmptrk.c33 = cc[14];
    cmptrk.c40 = cc[6];
    cmptrk.c41 = cc[7];
    cmptrk.c42 = cc[8];
    cmptrk.c43 = cc[13];
    cmptrk.c44 = cc[9];

    // fill tree
    cmptrktree->Fill();

  } // loop on tracks at TPC 1st hit

  cerr<<"+++\n+++ Number of compared tracks: "<<pi+ka+el+mu+pr<<endl;
  cerr<<"+++ Pions: "<<pi<<", Kaons: "<<ka<<", Protons : "<<pr<<", Electrons: "<<el<<", Muons: "<<mu<<endl;
  cerr<<"+++"<<endl;

  // Write tree to file
  TFile *outfile = new TFile(covmatName,"recreate");
  cmptrktree->Write();
  outfile->Close();


  // Write efficiencies to file
  FILE *effFile = fopen(tpceffName,"w");
  //fprintf(effFile,"%d\n",kalEntries);
  for(Int_t j=0; j<27; j++) {
    if(geaPi[j]>=5) effPi[j]=(Float_t)kalPi[j]/geaPi[j];
    if(geaKa[j]>=5) effKa[j]=(Float_t)kalKa[j]/geaKa[j];
    if(geaPr[j]>=5) effPr[j]=(Float_t)kalPr[j]/geaPr[j];
    if(geaEl[j]>=5) effEl[j]=(Float_t)kalEl[j]/geaEl[j];
    if(geaMu[j]>=5) effMu[j]=(Float_t)kalMu[j]/geaMu[j];
    fprintf(effFile,"%f  %f  %f  %f  %f\n",effPi[j],effKa[j],effPr[j],effEl[j],effMu[j]);
  }

  for(Int_t j=0; j<27; j++) {
    fprintf(effFile,"%d  %d  %d  %d  %d\n",geaPi[j],geaKa[j],geaPr[j],geaEl[j],geaMu[j]);
  }
  for(Int_t j=0; j<27; j++) {
    fprintf(effFile,"%d  %d  %d  %d  %d\n",kalPi[j],kalKa[j],kalPr[j],kalEl[j],kalMu[j]);
  }
  fclose(effFile);

  // delete AliRun object
  delete gAlice; gAlice=0;
  
  // close all input files
  kalFile->Close();
  geaFile->Close();
  galiceFile->Close();

  delete [] kalLab;

  return;
}
//-----------------------------------------------------------------------------
void AliTPCtrackerParam::CookdEdx(Double_t pt,Double_t eta) {
//-----------------------------------------------------------------------------
// This function assigns the track a dE/dx and makes a rough PID
//-----------------------------------------------------------------------------

  Double_t mean = fdEdxMean->GetValueAt(pt,eta);
  Double_t rms  = fdEdxRMS->GetValueAt(pt,eta);

  Double_t dEdx = gRandom->Gaus(mean,rms);

  fTrack.SetdEdx(dEdx);

  AliTPCtrackParam t(fTrack);

  //Very rough PID
  Double_t p = TMath::Sqrt(1.+t.GetTgl()*t.GetTgl())*pt;

  if (p<0.6) {
    if (dEdx < 39.+ 12./(p+0.25)/(p+0.25)) { 
      t.AssignMass(0.13957); new(&fTrack) AliTPCtrack(t); return;
    }
    if (dEdx < 39.+ 12./p/p) { 
      t.AssignMass(0.49368); new(&fTrack) AliTPCtrack(t); return;
    }
    t.AssignMass(0.93827); new(&fTrack) AliTPCtrack(t); return;
  }

  if (p<1.2) {
    if (dEdx < 39.+ 12./(p+0.25)/(p+0.25)) { 
      t.AssignMass(0.13957); new(&fTrack) AliTPCtrack(t); return;
    }
    t.AssignMass(0.93827); new(&fTrack) AliTPCtrack(t); return;
  }

  t.AssignMass(0.13957); new(&fTrack) AliTPCtrack(t); return;
}
//-----------------------------------------------------------------------------
void AliTPCtrackerParam::CookTrack(Double_t pt,Double_t eta) {
//-----------------------------------------------------------------------------
// This function deals with covariance matrix and smearing
//-----------------------------------------------------------------------------

  COVMATRIX covmat;
  Double_t  p,cosl;
  Double_t  trkKine[2],trkRegPar[4];    
  Double_t  xref,alpha,xx[5],xxsm[5],cc[15];
  Int_t     treeEntries;

  fCovTree->SetBranchAddress("matrix",&covmat);

  // get random entry from the tree
  treeEntries = (Int_t)fCovTree->GetEntries();
  fCovTree->GetEvent(gRandom->Integer(treeEntries));

  // get P and Cosl from track
  cosl = TMath::Cos(TMath::ATan(fTrack.GetTgl()));
  p    = 1./TMath::Abs(fTrack.Get1Pt())/cosl;

  trkKine[0] = p;
  trkKine[1] = cosl; 

  // get covariance matrix from regularized matrix    
  for(Int_t l=0;l<4;l++) trkRegPar[l] = (*fRegPar)(0,l);
  cc[0] = covmat.c00*RegFunc(trkKine,trkRegPar);
  cc[1] = covmat.c10;
  for(Int_t l=0;l<4;l++) trkRegPar[l] = (*fRegPar)(1,l);
  cc[2] = covmat.c11*RegFunc(trkKine,trkRegPar);
  for(Int_t l=0;l<4;l++) trkRegPar[l] = (*fRegPar)(2,l);
  cc[3] = covmat.c20*RegFunc(trkKine,trkRegPar);
  cc[4] = covmat.c21;
  for(Int_t l=0;l<4;l++) trkRegPar[l] = (*fRegPar)(3,l);
  cc[5] = covmat.c22*RegFunc(trkKine,trkRegPar);
  cc[6] = covmat.c30;
  for(Int_t l=0;l<4;l++) trkRegPar[l] = (*fRegPar)(4,l);
  cc[7] = covmat.c31*RegFunc(trkKine,trkRegPar);
  cc[8] = covmat.c32;
  for(Int_t l=0;l<4;l++) trkRegPar[l] = (*fRegPar)(5,l);
  cc[9] = covmat.c33*RegFunc(trkKine,trkRegPar);
  for(Int_t l=0;l<4;l++) trkRegPar[l] = (*fRegPar)(6,l);
  cc[10]= covmat.c40*RegFunc(trkKine,trkRegPar);
  cc[11]= covmat.c41;
  for(Int_t l=0;l<4;l++) trkRegPar[l] = (*fRegPar)(7,l);
  cc[12]= covmat.c42*RegFunc(trkKine,trkRegPar);
  cc[13]= covmat.c43;
  for(Int_t l=0;l<4;l++) trkRegPar[l] = (*fRegPar)(8,l);
  cc[14]= covmat.c44*RegFunc(trkKine,trkRegPar);
   
  TMatrixD covMatSmear(5,5);
    
  covMatSmear = GetSmearingMatrix(cc,pt,eta);

  // get track original parameters
  xref =fTrack.GetX();
  alpha=fTrack.GetAlpha();
  xx[0]=fTrack.GetY();
  xx[1]=fTrack.GetZ();
  xx[2]=fTrack.GetX()*fTrack.GetC()-fTrack.GetSnp();
  xx[3]=fTrack.GetTgl();
  xx[4]=fTrack.GetC();
    
  // use smearing matrix to smear the original parameters
  SmearTrack(xx,xxsm,covMatSmear);
    
  AliTPCtrack track(0,xxsm,cc,xref,alpha);
  new(&fTrack) AliTPCtrack(track);
  
  return; 
}
//-----------------------------------------------------------------------------
void AliTPCtrackerParam::DrawEffs(const Char_t* inName,Int_t pdg) {
//-----------------------------------------------------------------------------
// This function draws the TPC efficiencies in the [pT,eta] bins
//-----------------------------------------------------------------------------

  ReadEffs(inName);
  SetParticle(pdg);

  const Int_t n = fEff->GetPointsPt();
  Double_t * effsA = new Double_t[n];
  Double_t * effsB = new Double_t[n];
  Double_t * effsC = new Double_t[n];
  Double_t * pt = new Double_t[n];

  fEff->GetArrayPt(pt);
  for(Int_t i=0;i<n;i++) {
    effsA[i] = fEff->GetParam(i,0);
    effsB[i] = fEff->GetParam(i,1);
    effsC[i] = fEff->GetParam(i,2);
  }
  
  TGraph *grA = new TGraph(n,pt,effsA);
  TGraph *grB = new TGraph(n,pt,effsB);
  TGraph *grC = new TGraph(n,pt,effsC);

  grA->SetMarkerStyle(20); grA->SetMarkerColor(2); grA->SetMarkerSize(1.5);
  grB->SetMarkerStyle(21); grB->SetMarkerColor(3); grB->SetMarkerSize(1.5);
  grC->SetMarkerStyle(22); grC->SetMarkerColor(4); grC->SetMarkerSize(1.5);

  TString title("Distribution of the TPC efficiencies");
  switch (pdg) {
  case 211:
    title.Prepend("PIONS - ");
    break;
  case 321:
    title.Prepend("KAONS - ");
    break;
  case 2212:
    title.Prepend("PROTONS - ");
    break;
  case 11:
    title.Prepend("ELECTRONS - ");
    break;
  case 13:
    title.Prepend("MUONS - ");
    break;
  }


  gStyle->SetOptStat(0);
  TCanvas *c = new TCanvas("c","effs",0,0,900,700);
  c->SetLogx();
  c->SetGridx();
  c->SetGridy();

  TH2F *frame = new TH2F("frame",title.Data(),2,0.1,30,2,0,1);
  frame->SetXTitle("p_{T} [GeV/c]");
  frame->Draw();
  grA->Draw("P");
  grB->Draw("P");
  grC->Draw("P");

  TLegend *leg = new TLegend(0.2,0.2,0.4,0.4);
  leg->AddEntry(grA,"|#eta|<0.3","p");
  leg->AddEntry(grB,"0.3<|#eta|<0.6","p");
  leg->AddEntry(grC,"0.6<|#eta|<0.9","p");
  leg->SetFillColor(0);
  leg->Draw();

  delete [] effsA;
  delete [] effsB;
  delete [] effsC;
  delete [] pt;

  return;
}
//-----------------------------------------------------------------------------
void AliTPCtrackerParam::DrawPulls(const Char_t* inName,Int_t pdg,
				   Int_t par) {
//-----------------------------------------------------------------------------
// This function draws the pulls in the [pT,eta] bins
//-----------------------------------------------------------------------------

  ReadPulls(inName);
  SetParticle(pdg);

  const Int_t n = (fPulls+par)->GetPointsPt();
  Double_t * pullsA = new Double_t[n];
  Double_t * pullsB = new Double_t[n];
  Double_t * pullsC = new Double_t[n];
  Double_t * pt = new Double_t[n];
  (fPulls+par)->GetArrayPt(pt);  
  for(Int_t i=0;i<n;i++) {
    pullsA[i] = (fPulls+par)->GetParam(i,0);
    pullsB[i] = (fPulls+par)->GetParam(i,1);
    pullsC[i] = (fPulls+par)->GetParam(i,2);
  }

  TGraph *grA = new TGraph(n,pt,pullsA);
  TGraph *grB = new TGraph(n,pt,pullsB);
  TGraph *grC = new TGraph(n,pt,pullsC);

  grA->SetMarkerStyle(20); grA->SetMarkerColor(2); grA->SetMarkerSize(1.5);
  grB->SetMarkerStyle(21); grB->SetMarkerColor(3); grB->SetMarkerSize(1.5);
  grC->SetMarkerStyle(22); grC->SetMarkerColor(4); grC->SetMarkerSize(1.5);

  TString title("Distribution of the pulls: ");
  switch (pdg) {
  case 211:
    title.Prepend("PIONS - ");
    break;
  case 321:
    title.Prepend("KAONS - ");
    break;
  case 11:
    title.Prepend("ELECTRONS - ");
    break;
  }
  switch (par) {
  case 0:
    title.Append("y");
    break;
  case 1:
    title.Append("z");
    break;
  case 2:
    title.Append("    #eta");
    break;
  case 3:
    title.Append("tg    #lambda");
    break;
  case 4:
    title.Append("C");
    break;
  }

  gStyle->SetOptStat(0);
  TCanvas *c = new TCanvas("c","pulls",0,0,900,700);
  c->SetLogx();
  c->SetGridx();
  c->SetGridy();

  TH2F *frame = new TH2F("frame",title.Data(),2,0.1,30,2,0,2);
  frame->SetXTitle("p_{T} [GeV/c]");
  frame->Draw();
  grA->Draw("P");
  grB->Draw("P");
  grC->Draw("P");

  TLegend *leg = new TLegend(0.2,0.2,0.4,0.4);
  leg->AddEntry(grA,"|#eta|<0.3","p");
  leg->AddEntry(grB,"0.3<|#eta|<0.6","p");
  leg->AddEntry(grC,"0.6<|#eta|<0.9","p");
  leg->SetFillColor(0);
  leg->Draw();

  delete [] pullsA;
  delete [] pullsB;
  delete [] pullsC;
  delete [] pt;

  return;
}
//-----------------------------------------------------------------------------
Int_t AliTPCtrackerParam::GetBin(Double_t pt,Double_t eta) const {
//-----------------------------------------------------------------------------
// This function tells bin number in a grid (pT,eta) 
//-----------------------------------------------------------------------------
  if(TMath::Abs(eta)<0.3) {
    if(pt<0.3)            return 0;
    if(pt>=0.3 && pt<0.5) return 3;
    if(pt>=0.5 && pt<1.)  return 6;
    if(pt>=1. && pt<1.5)  return 9;
    if(pt>=1.5 && pt<3.)  return 12;
    if(pt>=3. && pt<5.)   return 15;
    if(pt>=5. && pt<7.)   return 18;
    if(pt>=7. && pt<15.)  return 21;
    if(pt>=15.)           return 24;
  }
  if(TMath::Abs(eta)>=0.3 && TMath::Abs(eta)<0.6) {
    if(pt<0.3)            return 1;
    if(pt>=0.3 && pt<0.5) return 4;
    if(pt>=0.5 && pt<1.)  return 7;
    if(pt>=1. && pt<1.5)  return 10;
    if(pt>=1.5 && pt<3.)  return 13;
    if(pt>=3. && pt<5.)   return 16;
    if(pt>=5. && pt<7.)   return 19;
    if(pt>=7. && pt<15.)  return 22;
    if(pt>=15.)           return 25;
  }
  if(TMath::Abs(eta)>=0.6) {
    if(pt<0.3)            return 2;
    if(pt>=0.3 && pt<0.5) return 5;
    if(pt>=0.5 && pt<1.)  return 8;
    if(pt>=1. && pt<1.5)  return 11;
    if(pt>=1.5 && pt<3.)  return 14;
    if(pt>=3. && pt<5.)   return 17;
    if(pt>=5. && pt<7.)   return 20;
    if(pt>=7. && pt<15.)  return 23;
    if(pt>=15.)           return 26;
  }

  return -1;

}
//-----------------------------------------------------------------------------
TMatrixD AliTPCtrackerParam::GetSmearingMatrix(Double_t* cc,Double_t pt,
					       Double_t eta) const {
//-----------------------------------------------------------------------------
// This function stretches the covariance matrix according to the pulls
//-----------------------------------------------------------------------------

  TMatrixD covMat(5,5);

  covMat(0,0)=cc[0];
  covMat(1,0)=cc[1];  covMat(0,1)=covMat(1,0);
  covMat(1,1)=cc[2];
  covMat(2,0)=cc[3];  covMat(0,2)=covMat(2,0);
  covMat(2,1)=cc[4];  covMat(1,2)=covMat(2,1);
  covMat(2,2)=cc[5];
  covMat(3,0)=cc[6];  covMat(0,3)=covMat(3,0);
  covMat(3,1)=cc[7];  covMat(1,3)=covMat(3,1);
  covMat(3,2)=cc[8];  covMat(2,3)=covMat(3,2);
  covMat(3,3)=cc[9];
  covMat(4,0)=cc[10]; covMat(0,4)=covMat(4,0);
  covMat(4,1)=cc[11]; covMat(1,4)=covMat(4,1);
  covMat(4,2)=cc[12]; covMat(2,4)=covMat(4,2);
  covMat(4,3)=cc[13]; covMat(3,4)=covMat(4,3);
  covMat(4,4)=cc[14];

  TMatrixD stretchMat(5,5);
  for(Int_t k=0;k<5;k++) {
    for(Int_t l=0;l<5;l++) {
      stretchMat(k,l)=0.;
    }
  }

  for(Int_t i=0;i<5;i++) stretchMat(i,i) = 
			   TMath::Sqrt((fPulls+i)->GetValueAt(pt,eta)); 

  TMatrixD mat(stretchMat,TMatrixD::kMult,covMat);
  TMatrixD covMatSmear(mat,TMatrixD::kMult,stretchMat);

  return covMatSmear;
}
//-----------------------------------------------------------------------------
void AliTPCtrackerParam::InitializeKineGrid(Option_t* which,Option_t* how) {
//-----------------------------------------------------------------------------
// This function initializes ([pt,eta] points) the data members AliTPCkineGrid
// which = "DB"     -> initialize fDBgrid... members
//         "eff"    -> initialize fEff... members
//         "pulls"  -> initialize fPulls... members
//         "dEdx"   -> initialize fdEdx... members
// how =   "points" -> initialize with the points of the grid
//         "bins"   -> initialize with the bins of the grid
//-----------------------------------------------------------------------------

  const char *points = strstr(how,"points");
  const char *bins   = strstr(how,"bins");
  const char *DB     = strstr(which,"DB");
  const char *eff    = strstr(which,"eff");
  const char *pulls  = strstr(which,"pulls");
  const char *dEdx   = strstr(which,"dEdx");


  Int_t nEta=0,nPtPi=0,nPtKa=0,nPtPr=0,nPtEl=0,nPtMu=0;

  Double_t etaPoints[2] = {0.3,0.6};
  Double_t etaBins[3]   = {0.15,0.45,0.75};
  Double_t ptPoints8[8] = {0.3,0.5,1.,1.5,3.,5.,7.,15.};
  Double_t ptPoints5[5] = {0.3,0.5,1.,1.5,5.};
  Double_t ptBins9[9]   = {0.244,0.390,0.676,1.190,2.36,4.,6.,10.,20.};
  Double_t ptBins6[6]   = {0.244,0.390,0.676,1.190,2.36,10.};

  Double_t *eta=0,*ptPi=0,*ptKa=0,*ptPr=0,*ptEl=0,*ptMu=0;

  if(points) {
    nEta  = 2;
    nPtPi = 8;
    nPtKa = 5;
    nPtPr = 4;
    nPtEl = 8;
    nPtMu = 2;
    eta  = etaPoints;
    ptPi = ptPoints8;
    ptKa = ptPoints5;
    ptPr = ptPoints8;
    ptEl = ptPoints8;
    ptMu = ptPoints8;
  }
  if(bins) {
    nEta  = 3;
    nPtPi = 9;
    nPtKa = 6;
    nPtPr = 5;
    nPtEl = 9;
    nPtMu = 3;
    eta  = etaBins;
    ptPi = ptBins9;
    ptKa = ptBins6;
    ptPr = ptBins9;
    ptEl = ptBins9;
    ptMu = ptBins9;
  }

  AliTPCkineGrid *dummy=0;

  if(DB) {    
    dummy = new AliTPCkineGrid(nPtPi,nEta,ptPi,eta);
    new(&fDBgridPi) AliTPCkineGrid(*dummy);
    delete dummy;
    dummy = new AliTPCkineGrid(nPtKa,nEta,ptKa,eta);
    new(&fDBgridKa) AliTPCkineGrid(*dummy);
    delete dummy;
    dummy = new AliTPCkineGrid(nPtEl,nEta,ptEl,eta);
    new(&fDBgridEl) AliTPCkineGrid(*dummy);
    delete dummy;
  }
  if(eff) {    
    dummy = new AliTPCkineGrid(nPtPi,nEta,ptPi,eta);
    new(&fEffPi) AliTPCkineGrid(*dummy);
    delete dummy;
    dummy = new AliTPCkineGrid(nPtKa,nEta,ptKa,eta);
    new(&fEffKa) AliTPCkineGrid(*dummy);
    delete dummy;
    dummy = new AliTPCkineGrid(nPtPr,nEta,ptPr,eta);
    new(&fEffPr) AliTPCkineGrid(*dummy);
    delete dummy;
    dummy = new AliTPCkineGrid(nPtEl,nEta,ptEl,eta);
    new(&fEffEl) AliTPCkineGrid(*dummy);
    delete dummy;
    dummy = new AliTPCkineGrid(nPtMu,nEta,ptMu,eta);
    new(&fEffMu) AliTPCkineGrid(*dummy);
    delete dummy;
  }
  if(pulls) {    
    dummy = new AliTPCkineGrid(nPtPi,nEta,ptPi,eta);
    for(Int_t i=0;i<5;i++) new(&fPullsPi[i]) AliTPCkineGrid(*dummy);
    delete dummy;
    dummy = new AliTPCkineGrid(nPtKa,nEta,ptKa,eta);
    for(Int_t i=0;i<5;i++) new(&fPullsKa[i]) AliTPCkineGrid(*dummy);
    delete dummy;
    dummy = new AliTPCkineGrid(nPtEl,nEta,ptEl,eta);
    for(Int_t i=0;i<5;i++) new(&fPullsEl[i]) AliTPCkineGrid(*dummy);
    delete dummy;
  }
  if(dEdx) {    
    dummy = new AliTPCkineGrid(nPtPi,nEta,ptPi,eta);
    new(&fdEdxMeanPi) AliTPCkineGrid(*dummy);
    new(&fdEdxRMSPi) AliTPCkineGrid(*dummy);
    delete dummy;
    dummy = new AliTPCkineGrid(nPtKa,nEta,ptKa,eta);
    new(&fdEdxMeanKa) AliTPCkineGrid(*dummy);
    new(&fdEdxRMSKa) AliTPCkineGrid(*dummy);
    delete dummy;
    dummy = new AliTPCkineGrid(nPtPi,nEta,ptPi,eta);
    new(&fdEdxMeanPr) AliTPCkineGrid(*dummy);
    new(&fdEdxRMSPr) AliTPCkineGrid(*dummy);
    delete dummy;
    dummy = new AliTPCkineGrid(nPtEl,nEta,ptEl,eta);
    new(&fdEdxMeanEl) AliTPCkineGrid(*dummy);
    new(&fdEdxRMSEl) AliTPCkineGrid(*dummy);
    delete dummy;
  }

  return;
}
//-----------------------------------------------------------------------------
void AliTPCtrackerParam::MakeDataBase() {
//-----------------------------------------------------------------------------
// This function creates the DB file and store in it:
//  - TPC Efficiencies for pi,ka,pr,el,mu     (AliTPCkineGrid class)
//  - Pulls for pi,ka,el                      (AliTPCkineGrid class)
//  - Regularization parameters for pi,ka,el  (TMatrixD class)
//  - dE/dx parameterization for pi,ka,pr,el  (AliTPCkineGrid class)  
//  - Regularized cov. matrices for pi,ka,el  (COVMATRIX structure)
//-----------------------------------------------------------------------------

  // define some file names
  Char_t *effFile  ="CovMatrixDB_partial.root";
  Char_t *pullsFile="CovMatrixDB_partial.root";
  Char_t *regPiFile="CovMatrixDB_partial.root";
  Char_t *regKaFile="CovMatrixDB_partial.root";
  Char_t *regElFile="CovMatrixDB_partial.root";
  Char_t *dEdxPiFile="dEdxPi.root";
  Char_t *dEdxKaFile="dEdxKa.root";
  Char_t *dEdxPrFile="dEdxPr.root";
  Char_t *dEdxElFile="dEdxEl.root";
  Char_t *cmFile1  ="CovMatrix_AllEvts_1.root";
  Char_t *cmFile2  ="CovMatrix_AllEvts_2.root";
  Char_t *cmFile3  ="CovMatrix_AllEvts_3.root";

  // store the effieciencies
  ReadEffs(effFile);
  WriteEffs(fDBfileName.Data());

  // store the pulls
  ReadPulls(pullsFile);
  WritePulls(fDBfileName.Data());

  //store the regularization parameters
  ReadRegParams(regPiFile,211);
  WriteRegParams(fDBfileName.Data(),211);
  ReadRegParams(regKaFile,321);
  WriteRegParams(fDBfileName.Data(),321);
  ReadRegParams(regElFile,11);
  WriteRegParams(fDBfileName.Data(),11);

  // store the dEdx parameters
  ReaddEdx(dEdxPiFile,211);
  WritedEdx(fDBfileName.Data(),211);
  ReaddEdx(dEdxKaFile,321);
  WritedEdx(fDBfileName.Data(),321);
  ReaddEdx(dEdxPrFile,2212);
  WritedEdx(fDBfileName.Data(),2212);
  ReaddEdx(dEdxElFile,11);
  WritedEdx(fDBfileName.Data(),11);


  //
  // store the regularized covariance matrices
  //
  InitializeKineGrid("DB","points");

  const Int_t nBinsPi = fDBgridPi.GetTotBins();
  const Int_t nBinsKa = fDBgridKa.GetTotBins();
  const Int_t nBinsEl = fDBgridEl.GetTotBins();


  // create the trees for cov. matrices
  // trees for pions
  TTree *CovTreePi_ = NULL;
  CovTreePi_ = new TTree[nBinsPi]; 
  // trees for kaons
  TTree *CovTreeKa_ = NULL;
  CovTreeKa_ = new TTree[nBinsKa]; 
  // trees for electrons
  TTree *CovTreeEl_ = NULL;
  CovTreeEl_ = new TTree[nBinsEl]; 

  Char_t hname[100], htitle[100];
  COVMATRIX covmat;

  
  for(Int_t i=0; i<nBinsPi; i++) {
    sprintf(hname,"CovTreePi_bin%d",i);
    sprintf(htitle,"Tree with cov matrix elements for bin %d",i);
    CovTreePi_[i].SetName(hname); CovTreePi_[i].SetTitle(htitle);
    CovTreePi_[i].Branch("matrix",&covmat,"c00/D:c10:c11:c20:c21:c22:c30:c31:c32:c33:c40:c41:c42:c43:c44",5000000);
  }
  for(Int_t i=0; i<nBinsKa; i++) {
    sprintf(hname,"CovTreeKa_bin%d",i);
    sprintf(htitle,"Tree with cov matrix elements for bin %d",i);
    CovTreeKa_[i].SetName(hname); CovTreeKa_[i].SetTitle(htitle);
    CovTreeKa_[i].Branch("matrix",&covmat,"c00/D:c10:c11:c20:c21:c22:c30:c31:c32:c33:c40:c41:c42:c43:c44",1000000);
  }
  for(Int_t i=0; i<nBinsEl; i++) {
    sprintf(hname,"CovTreeEl_bin%d",i);
    sprintf(htitle,"Tree with cov matrix elements for bin %d",i);
    CovTreeEl_[i].SetName(hname); CovTreeEl_[i].SetTitle(htitle);
    CovTreeEl_[i].Branch("matrix",&covmat,"c00/D:c10:c11:c20:c21:c22:c30:c31:c32:c33:c40:c41:c42:c43:c44",1000000);
  }
  
  // create the chain with the compared tracks
  TChain cmptrkchain("cmptrktree");
  cmptrkchain.Add(cmFile1);
  cmptrkchain.Add(cmFile2);
  cmptrkchain.Add(cmFile3);

  COMPTRACK cmptrk; 
  cmptrkchain.SetBranchAddress("comptracks",&cmptrk);
  Int_t entries = (Int_t)cmptrkchain.GetEntries(); 
  cerr<<" Number of entries: "<<entries<<endl;

  Int_t trkPdg,trkBin;
  Double_t trkKine[2],trkRegPar[4]; 
  Int_t * nPerBinPi = new Int_t[nBinsPi];
  for(Int_t k=0;k<nBinsPi;k++) nPerBinPi[k]=0;
  Int_t * nPerBinKa = new Int_t[nBinsKa];
  for(Int_t k=0;k<nBinsKa;k++) nPerBinKa[k]=0;
  Int_t * nPerBinEl = new Int_t[nBinsEl];
  for(Int_t k=0;k<nBinsEl;k++) nPerBinEl[k]=0;

  // loop on chain entries 
  for(Int_t l=0; l<entries; l++) {
    if(l % 10000 == 0) cerr<<"--- Processing track "<<l<<" of "<<entries<<" ---"<<endl;
    // get the event
    cmptrkchain.GetEvent(l);
    // get the pdg code
    trkPdg = TMath::Abs(cmptrk.pdg);
    // use only pions, kaons, electrons
    if(trkPdg!=211 && trkPdg!=321 && trkPdg!=11) continue;
    SetParticle(trkPdg);
    trkBin = fDBgrid->GetBin(cmptrk.pt,cmptrk.eta);
    //cerr<<cmptrk.pt<<"  "<<cmptrk.eta<<"  "<<trkBin<<endl;

    if(trkPdg==211 && nPerBinPi[trkBin]>=5000) continue;
    if(trkPdg==321 && nPerBinKa[trkBin]>=5000) continue;
    if(trkPdg==11  && nPerBinEl[trkBin]>=5000) continue;

    trkKine[0] = cmptrk.p;
    trkKine[1] = cmptrk.cosl;
    
    // get regularized covariance matrix
    for(Int_t k=0;k<4;k++) trkRegPar[k] = (*fRegPar)(0,k);
    covmat.c00 = cmptrk.c00/RegFunc(trkKine,trkRegPar);
    covmat.c10 = cmptrk.c10;
    for(Int_t k=0;k<4;k++) trkRegPar[k] = (*fRegPar)(1,k);
    covmat.c11 = cmptrk.c11/RegFunc(trkKine,trkRegPar);
    for(Int_t k=0;k<4;k++) trkRegPar[k] = (*fRegPar)(2,k);
    covmat.c20 = cmptrk.c20/RegFunc(trkKine,trkRegPar);
    covmat.c21 = cmptrk.c21;
    for(Int_t k=0;k<4;k++) trkRegPar[k] = (*fRegPar)(3,k);
    covmat.c22 = cmptrk.c22/RegFunc(trkKine,trkRegPar);
    covmat.c30 = cmptrk.c30;
    for(Int_t k=0;k<4;k++) trkRegPar[k] = (*fRegPar)(4,k);
    covmat.c31 = cmptrk.c31/RegFunc(trkKine,trkRegPar);
    covmat.c32 = cmptrk.c32;
    for(Int_t k=0;k<4;k++) trkRegPar[k] = (*fRegPar)(5,k);
    covmat.c33 = cmptrk.c33/RegFunc(trkKine,trkRegPar);
    for(Int_t k=0;k<4;k++) trkRegPar[k] = (*fRegPar)(6,k);
    covmat.c40 = cmptrk.c40/RegFunc(trkKine,trkRegPar);
    covmat.c41 = cmptrk.c41;
    for(Int_t k=0;k<4;k++) trkRegPar[k] = (*fRegPar)(7,k);
    covmat.c42 = cmptrk.c42/RegFunc(trkKine,trkRegPar);
    covmat.c43 = cmptrk.c43;
    for(Int_t k=0;k<4;k++) trkRegPar[k] = (*fRegPar)(8,k);
    covmat.c44 = cmptrk.c44/RegFunc(trkKine,trkRegPar);

    // fill the tree
    switch (trkPdg) {
    case 211: // pions
      CovTreePi_[trkBin].Fill();
      nPerBinPi[trkBin]++;
      break;
    case 321: // kaons
      CovTreeKa_[trkBin].Fill();
      nPerBinKa[trkBin]++;
      break;
    case 11: // electrons
      CovTreeEl_[trkBin].Fill();
      nPerBinEl[trkBin]++;
      break;
    }
  } // loop on chain entries

  // store all trees the DB file
  TFile *DBfile = new TFile(fDBfileName.Data(),"update");
  DBfile->mkdir("CovMatrices");
  gDirectory->cd("/CovMatrices");
  gDirectory->mkdir("Pions");
  gDirectory->mkdir("Kaons");
  gDirectory->mkdir("Electrons");
  // store pions
  gDirectory->cd("/CovMatrices/Pions");
  fDBgridPi.SetName("DBgridPi"); fDBgridPi.Write();
  for(Int_t i=0;i<nBinsPi;i++) CovTreePi_[i].Write();
  // store kaons
  gDirectory->cd("/CovMatrices/Kaons");
  fDBgridKa.SetName("DBgridKa"); fDBgridKa.Write();
  for(Int_t i=0;i<nBinsKa;i++) CovTreeKa_[i].Write();
  // store electrons
  gDirectory->cd("/CovMatrices/Electrons");
  fDBgridEl.SetName("DBgridEl"); fDBgridEl.Write();
  for(Int_t i=0;i<nBinsEl;i++) CovTreeEl_[i].Write();

  DBfile->Close();
  delete [] nPerBinPi;
  delete [] nPerBinKa;
  delete [] nPerBinEl;
  
  return;
}
//-----------------------------------------------------------------------------
void AliTPCtrackerParam::MergeEvents(Int_t evFirst,Int_t evLast) {
//-----------------------------------------------------------------------------
// This function: 1) merges the files from track comparison
//                   (beware: better no more than 100 events per file)
//                2) computes the average TPC efficiencies
//-----------------------------------------------------------------------------

  Char_t *outName="TPCeff.root";

  Int_t    nCol;
  Float_t effPi,effKa,effPr,effEl,effMu;
  Float_t avEffPi[27],avEffKa[27],avEffPr[27],avEffEl[27],avEffMu[27];
  Int_t    evtsPi[27],evtsKa[27],evtsPr[27],evtsEl[27],evtsMu[27];
  for(Int_t j=0; j<27; j++) {
    avEffPi[j]=avEffKa[j]=avEffPr[j]=avEffEl[j]=avEffMu[j]=0.;
    evtsPi[j]=evtsKa[j]=evtsPr[j]=evtsEl[j]=evtsMu[j]=0;
  }

  // create the chain for the tree of compared tracks
  TChain ch1("cmptrktree");
  TChain ch2("cmptrktree");
  TChain ch3("cmptrktree");


  for(Int_t j=evFirst; j<=evLast; j++) {
    cerr<<"Processing event "<<j<<endl;

    TString covName("CovMatrix.");
    TString effName("TPCeff.");
    covName+=j;
    effName+=j;
    covName.Append(".root");
    effName.Append(".dat");

    if(gSystem->AccessPathName(covName.Data(),kFileExists) ||
       gSystem->AccessPathName(effName.Data(),kFileExists)) continue;

    if(j<=100) ch1.Add(covName.Data());
    if(j>100 && j<=200) ch2.Add(covName.Data());
    if(j>200) ch3.Add(covName.Data());


    FILE *effIn = fopen(effName.Data(),"r");
    for(Int_t k=0; k<27; k++) {
      nCol = fscanf(effIn,"%f  %f  %f  %f  %f",&effPi,&effKa,&effPr,&effEl,&effMu);
      if(effPi>=0.) {avEffPi[k]+=effPi; evtsPi[k]++;}
      if(effKa>=0.) {avEffKa[k]+=effKa; evtsKa[k]++;}
      if(effPr>=0.) {avEffPr[k]+=effPr; evtsPr[k]++;}
      if(effEl>=0.) {avEffEl[k]+=effEl; evtsEl[k]++;}
      if(effMu>=0.) {avEffMu[k]+=effMu; evtsMu[k]++;}
    }
    fclose(effIn);

  } // loop on events

  // merge chain in one file
  TFile *covOut=0;
  covOut = new TFile("CovMatrix_AllEvts_1.root","recreate");
  ch1.Merge(covOut,1000000000);
  covOut->Close();
  delete covOut;
  covOut = new TFile("CovMatrix_AllEvts_2.root","recreate");
  ch2.Merge(covOut,1000000000);
  covOut->Close();
  delete covOut;
  covOut = new TFile("CovMatrix_AllEvts_3.root","recreate");
  ch3.Merge(covOut,1000000000);
  covOut->Close();
  delete covOut;


  // efficiencies 

  InitializeKineGrid("eff","bins");

  // pions
  for(Int_t j=0; j<27; j++) {
    if(evtsPi[j]==0) evtsPi[j]++;
    fEffPi.SetParam(j,(Double_t)avEffPi[j]/evtsPi[j]);
  }
 
  // kaons
  Int_t j=0;
  for(Int_t k=0; k<27; k++) {
    if(k>14 && k!=21 && k!=22 && k!=23) continue;
    if(evtsKa[k]==0) evtsKa[k]++;
    fEffKa.SetParam(j,(Double_t)avEffKa[k]/evtsKa[k]);
    j++;
  }

  // protons
  for(Int_t j=0; j<15; j++) {
    if(evtsPr[j]==0) evtsPr[j]++;
    fEffPr.SetParam(j,(Double_t)avEffPr[j]/evtsPr[j]);
  }

  // electrons
  for(Int_t j=0; j<27; j++) {
    if(evtsEl[j]==0) evtsEl[j]++;
    fEffEl.SetParam(j,(Double_t)avEffEl[j]/evtsEl[j]);
  }

  // muons
  for(Int_t j=0; j<9; j++) {
    if(evtsMu[j]==0) evtsMu[j]++;
    fEffMu.SetParam(j,(Double_t)avEffMu[j]/evtsMu[j]);
  }

  // write efficiencies to a file
  WriteEffs(outName);

  return;
}
//-----------------------------------------------------------------------------
Int_t AliTPCtrackerParam::ReadAllData(const Char_t* inName) {
//-----------------------------------------------------------------------------
// This function reads all parameters from the DB
//-----------------------------------------------------------------------------

  if(!ReadEffs(inName)) return 0;
  if(!ReadPulls(inName)) return 0;        
  if(!ReadRegParams(inName,211)) return 0; 
  if(!ReadRegParams(inName,321)) return 0; 
  if(!ReadRegParams(inName,11)) return 0; 
  if(!ReaddEdx(inName,211)) return 0;
  if(!ReaddEdx(inName,321)) return 0;
  if(!ReaddEdx(inName,2212)) return 0;
  if(!ReaddEdx(inName,11)) return 0;
  if(!ReadDBgrid(inName)) return 0;

  return 1;
}
//-----------------------------------------------------------------------------
Int_t AliTPCtrackerParam::ReaddEdx(const Char_t* inName,Int_t pdg) {
//-----------------------------------------------------------------------------
// This function reads the dEdx parameters from the DB
//-----------------------------------------------------------------------------

  if(gSystem->AccessPathName(inName,kFileExists)) { 
    cerr<<"AliTPCtrackerParam::ReaddEdx: "<<inName<<" not found\n"; 
    return 0; }

  TFile *inFile = TFile::Open(inName);
  switch (pdg) {
  case 211:
    fdEdxMean = (AliTPCkineGrid*)inFile->Get("/dEdx/Pions/dEdxMeanPi");
    fdEdxMeanPi.~AliTPCkineGrid();
    new(&fdEdxMeanPi) AliTPCkineGrid(*fdEdxMean);
    fdEdxRMS = (AliTPCkineGrid*)inFile->Get("/dEdx/Pions/dEdxRMSPi");
    fdEdxRMSPi.~AliTPCkineGrid();
    new(&fdEdxRMSPi) AliTPCkineGrid(*fdEdxRMS);
    break;
  case 321:
    fdEdxMean = (AliTPCkineGrid*)inFile->Get("/dEdx/Kaons/dEdxMeanKa");
    fdEdxMeanKa.~AliTPCkineGrid();
    new(&fdEdxMeanKa) AliTPCkineGrid(*fdEdxMean);
    fdEdxRMS = (AliTPCkineGrid*)inFile->Get("/dEdx/Kaons/dEdxRMSKa");
    fdEdxRMSKa.~AliTPCkineGrid();
    new(&fdEdxRMSKa) AliTPCkineGrid(*fdEdxRMS);
    break;
  case 2212:
    fdEdxMean = (AliTPCkineGrid*)inFile->Get("/dEdx/Protons/dEdxMeanPr");
    fdEdxMeanPr.~AliTPCkineGrid();
    new(&fdEdxMeanPr) AliTPCkineGrid(*fdEdxMean);
    fdEdxRMS = (AliTPCkineGrid*)inFile->Get("/dEdx/Protons/dEdxRMSPr");
    fdEdxRMSPr.~AliTPCkineGrid();
    new(&fdEdxRMSPr) AliTPCkineGrid(*fdEdxRMS);
    break;
  case 11:
    fdEdxMean = (AliTPCkineGrid*)inFile->Get("/dEdx/Electrons/dEdxMeanEl");
    fdEdxMeanEl.~AliTPCkineGrid();
    new(&fdEdxMeanEl) AliTPCkineGrid(*fdEdxMean);
    fdEdxRMS = (AliTPCkineGrid*)inFile->Get("/dEdx/Electrons/dEdxRMSEl");
    fdEdxRMSEl.~AliTPCkineGrid();
    new(&fdEdxRMSEl) AliTPCkineGrid(*fdEdxRMS);
    break;
  }
  inFile->Close();

  return 1;
}
//-----------------------------------------------------------------------------
Int_t AliTPCtrackerParam::ReadDBgrid(const Char_t* inName) {
//-----------------------------------------------------------------------------
// This function reads the kine grid from the DB
//-----------------------------------------------------------------------------

  if(gSystem->AccessPathName(inName,kFileExists)) {
    cerr<<"AliTPCtrackerParam::ReadCovMatrices: "<<inName<<" not found\n"; 
    return 0; }

  TFile *inFile = TFile::Open(inName);

  // first read the DB grid for the different particles
  fDBgrid = (AliTPCkineGrid*)inFile->Get("/CovMatrices/Pions/DBgridPi");
  fDBgridPi.~AliTPCkineGrid();
  new(&fDBgridPi) AliTPCkineGrid(*fDBgrid);
  fDBgrid = (AliTPCkineGrid*)inFile->Get("/CovMatrices/Kaons/DBgridKa");
  fDBgridKa.~AliTPCkineGrid();
  new(&fDBgridKa) AliTPCkineGrid(*fDBgrid);
  fDBgrid = (AliTPCkineGrid*)inFile->Get("/CovMatrices/Electrons/DBgridEl");
  fDBgridEl.~AliTPCkineGrid();
  new(&fDBgridEl) AliTPCkineGrid(*fDBgrid);

  inFile->Close();

  return 1;
}
//-----------------------------------------------------------------------------
Int_t AliTPCtrackerParam::ReadEffs(const Char_t* inName) {
//-----------------------------------------------------------------------------
// This function reads the TPC efficiencies from the DB
//-----------------------------------------------------------------------------

  if(gSystem->AccessPathName(inName,kFileExists)) {
    cerr<<"AliTPCtrackerParam::ReadEffs: "<<inName<<" not found\n"; 
    return 0; }

  TFile *inFile = TFile::Open(inName);

  fEff = (AliTPCkineGrid*)inFile->Get("/Efficiencies/Pions/EffPi");
  fEffPi.~AliTPCkineGrid();
  new(&fEffPi) AliTPCkineGrid(*fEff);
  fEff = (AliTPCkineGrid*)inFile->Get("/Efficiencies/Kaons/EffKa");
  fEffKa.~AliTPCkineGrid();
  new(&fEffKa) AliTPCkineGrid(*fEff);
  fEff = (AliTPCkineGrid*)inFile->Get("/Efficiencies/Protons/EffPr");
  fEffPr.~AliTPCkineGrid();
  new(&fEffPr) AliTPCkineGrid(*fEff);
  fEff = (AliTPCkineGrid*)inFile->Get("/Efficiencies/Electrons/EffEl");
  fEffEl.~AliTPCkineGrid();
  new(&fEffEl) AliTPCkineGrid(*fEff);
  fEff = (AliTPCkineGrid*)inFile->Get("/Efficiencies/Muons/EffMu");
  fEffMu.~AliTPCkineGrid();
  new(&fEffMu) AliTPCkineGrid(*fEff);

  inFile->Close();

  return 1;
}
//-----------------------------------------------------------------------------
Int_t AliTPCtrackerParam::ReadPulls(const Char_t* inName) {
//-----------------------------------------------------------------------------
// This function reads the pulls from the DB
//-----------------------------------------------------------------------------

  if(gSystem->AccessPathName(inName,kFileExists)) { 
    cerr<<"AliTPCtrackerParam::ReadPulls: "<<inName<<" not found\n"; 
    return 0; }

  TFile *inFile = TFile::Open(inName);

  for(Int_t i=0; i<5; i++) {
    TString pi("/Pulls/Pions/PullsPi_"); pi+=i;
    TString ka("/Pulls/Kaons/PullsKa_"); ka+=i;
    TString el("/Pulls/Electrons/PullsEl_"); el+=i;
    fPulls = (AliTPCkineGrid*)inFile->Get(pi.Data());
    fPullsPi[i].~AliTPCkineGrid();
    new(&fPullsPi[i]) AliTPCkineGrid(*fPulls);
    fPulls = (AliTPCkineGrid*)inFile->Get(ka.Data());
    fPullsKa[i].~AliTPCkineGrid();
    new(&fPullsKa[i]) AliTPCkineGrid(*fPulls);
    fPulls = (AliTPCkineGrid*)inFile->Get(el.Data());
    fPullsEl[i].~AliTPCkineGrid();
    new(&fPullsEl[i]) AliTPCkineGrid(*fPulls);
  }

  inFile->Close();

  return 1;
}
//-----------------------------------------------------------------------------
Int_t AliTPCtrackerParam::ReadRegParams(const Char_t* inName,Int_t pdg) {
//-----------------------------------------------------------------------------
// This function reads the regularization parameters from the DB
//-----------------------------------------------------------------------------

  if(gSystem->AccessPathName(inName,kFileExists)) { 
    cerr<<"AliTPCtrackerParam::ReadRegParams: "<<inName<<" not found\n"; 
    return 0; }

  TFile *inFile = TFile::Open(inName);
  switch (pdg) {
  case 211:
    fRegPar = (TMatrixD*)inFile->Get("/RegParams/Pions/RegPions");
    new(&fRegParPi) TMatrixD(*fRegPar);
    break;
  case 321:
    fRegPar = (TMatrixD*)inFile->Get("/RegParams/Kaons/RegKaons");
    new(&fRegParKa) TMatrixD(*fRegPar);
    break;
  case 11:
    fRegPar = (TMatrixD*)inFile->Get("/RegParams/Electrons/RegElectrons");
    new(&fRegParEl) TMatrixD(*fRegPar);
    break;
  }
  inFile->Close();

  return 1;
}
//-----------------------------------------------------------------------------
void AliTPCtrackerParam::RegularizeCovMatrix(const Char_t *outName,Int_t pdg) {
//-----------------------------------------------------------------------------
// This function regularizes the elements of the covariance matrix
// that show a momentum depence:
// c00, c11, c22, c33, c44, c20, c24, c40, c31 
// The regularization is done separately for pions, kaons, electrons:
// give "Pion","Kaon" or "Electron" as first argument.
//-----------------------------------------------------------------------------

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(10001);

  Int_t thisPdg=211;
  Char_t *part="Pions - ";

  InitializeKineGrid("DB","points");
  SetParticle(pdg);
  const Int_t fitbins = fDBgrid->GetBinsPt();
  cerr<<" Fit bins:  "<<fitbins<<endl;

  switch (pdg) {
  case 211: // pions
    thisPdg=211;  
    part="Pions - ";
    cerr<<" Processing pions ...\n";
    break;
  case 321: //kaons
    thisPdg=321; 
    part="Kaons - ";
    cerr<<" Processing kaons ...\n";
    break;
  case 11: // electrons
    thisPdg= 11; 
    part="Electrons - ";
    cerr<<" Processing electrons ...\n";
    break;
  }

  // create the chain with all compared tracks
  TChain cmptrkchain("cmptrktree");
  cmptrkchain.Add("CovMatrix_AllEvts_1.root");
  cmptrkchain.Add("CovMatrix_AllEvts_2.root");
  cmptrkchain.Add("CovMatrix_AllEvts_3.root");

  COMPTRACK cmptrk; 
  cmptrkchain.SetBranchAddress("comptracks",&cmptrk);
  Int_t entries = (Int_t)cmptrkchain.GetEntries(); 

    
  Int_t pbin;
  Int_t * n = new Int_t[fitbins];
  Int_t * n00 = new Int_t[fitbins];
  Int_t * n11 = new Int_t[fitbins];
  Int_t * n20 = new Int_t[fitbins];
  Int_t * n22 = new Int_t[fitbins];
  Int_t * n31 = new Int_t[fitbins];
  Int_t * n33 = new Int_t[fitbins];
  Int_t * n40 = new Int_t[fitbins];
  Int_t * n42 = new Int_t[fitbins];
  Int_t * n44 = new Int_t[fitbins];
  Double_t * p = new Double_t[fitbins];
  Double_t * ep = new Double_t[fitbins];
  Double_t * mean00 = new Double_t[fitbins];
  Double_t * mean11 = new Double_t[fitbins];
  Double_t * mean20 = new Double_t[fitbins];
  Double_t * mean22 = new Double_t[fitbins];
  Double_t * mean31 = new Double_t[fitbins];
  Double_t * mean33 = new Double_t[fitbins];
  Double_t * mean40 = new Double_t[fitbins];
  Double_t * mean42 = new Double_t[fitbins];
  Double_t * mean44 = new Double_t[fitbins];
  Double_t * sigma00 = new Double_t[fitbins];
  Double_t * sigma11 = new Double_t[fitbins];
  Double_t * sigma20 = new Double_t[fitbins];
  Double_t * sigma22 = new Double_t[fitbins];
  Double_t * sigma31 = new Double_t[fitbins];
  Double_t * sigma33 = new Double_t[fitbins];
  Double_t * sigma40 = new Double_t[fitbins];
  Double_t * sigma42 = new Double_t[fitbins];
  Double_t * sigma44 = new Double_t[fitbins];
  Double_t * rmean = new Double_t[fitbins];
  Double_t * rsigma = new Double_t[fitbins];
  Double_t fitpar[4];

  for(Int_t l=0; l<fitbins; l++) {
    n[l]=1;
    n00[l]=n11[l]=n20[l]=n22[l]=n31[l]=n33[l]=n40[l]=n42[l]=n44[l]=1;
    p[l ]=ep[l]=0.;
    mean00[l]=mean11[l]=mean20[l]=mean22[l]=mean31[l]=mean33[l]=mean40[l]=mean42[l]=mean44[l]=0.;
    sigma00[l]=sigma11[l]=sigma20[l]=sigma22[l]=sigma31[l]=sigma33[l]=sigma40[l]=sigma42[l]=sigma44[l]=0.;
  }

  // loop on chain entries for mean 
  for(Int_t l=0; l<entries; l++) {
    cmptrkchain.GetEvent(l);
    if(TMath::Abs(cmptrk.pdg)!=thisPdg) continue;
    pbin = (Int_t)fDBgrid->GetBin(cmptrk.pt,cmptrk.eta)/fDBgrid->GetBinsEta();
    n[pbin]++;
    p[pbin]+=cmptrk.p;
    mean00[pbin]+=cmptrk.c00;
    mean11[pbin]+=cmptrk.c11;
    mean20[pbin]+=cmptrk.c20;
    mean22[pbin]+=cmptrk.c22;
    mean31[pbin]+=cmptrk.c31;
    mean33[pbin]+=cmptrk.c33;
    mean40[pbin]+=cmptrk.c40;
    mean42[pbin]+=cmptrk.c42;
    mean44[pbin]+=cmptrk.c44;
  } // loop on chain entries

  for(Int_t l=0; l<fitbins; l++) {
    p[l]/=n[l];
    mean00[l]/=n[l];
    mean11[l]/=n[l];
    mean20[l]/=n[l];
    mean22[l]/=n[l];
    mean31[l]/=n[l];
    mean33[l]/=n[l];
    mean40[l]/=n[l];
    mean42[l]/=n[l];
    mean44[l]/=n[l];
  }

  // loop on chain entries for sigma
  for(Int_t l=0; l<entries; l++) {
    cmptrkchain.GetEvent(l);
    if(TMath::Abs(cmptrk.pdg)!=thisPdg) continue;
    pbin = (Int_t)fDBgrid->GetBin(cmptrk.pt,cmptrk.eta)/fDBgrid->GetBinsEta();
    if(TMath::Abs(cmptrk.c00-mean00[pbin])<0.4*mean00[pbin]) { n00[pbin]++;
      sigma00[pbin]+=(cmptrk.c00-mean00[pbin])*(cmptrk.c00-mean00[pbin]); }
    if(TMath::Abs(cmptrk.c11-mean11[pbin])<0.4*mean11[pbin]) { n11[pbin]++;
      sigma11[pbin]+=(cmptrk.c11-mean11[pbin])*(cmptrk.c11-mean11[pbin]); }
    if(TMath::Abs(cmptrk.c20-mean20[pbin])<0.4*mean20[pbin]) { n20[pbin]++;
      sigma20[pbin]+=(cmptrk.c20-mean20[pbin])*(cmptrk.c20-mean20[pbin]); }
    if(TMath::Abs(cmptrk.c22-mean22[pbin])<0.4*mean22[pbin]) { n22[pbin]++;
      sigma22[pbin]+=(cmptrk.c22-mean22[pbin])*(cmptrk.c22-mean22[pbin]); }
    if(TMath::Abs(cmptrk.c31-mean31[pbin])<-0.4*mean31[pbin]) { n31[pbin]++;
      sigma31[pbin]+=(cmptrk.c31-mean31[pbin])*(cmptrk.c31-mean31[pbin]); }
    if(TMath::Abs(cmptrk.c33-mean33[pbin])<0.4*mean33[pbin]) { n33[pbin]++;
      sigma33[pbin]+=(cmptrk.c33-mean33[pbin])*(cmptrk.c33-mean33[pbin]); }
    if(TMath::Abs(cmptrk.c40-mean40[pbin])<0.4*mean40[pbin]) { n40[pbin]++;
      sigma40[pbin]+=(cmptrk.c40-mean40[pbin])*(cmptrk.c40-mean40[pbin]); }
    if(TMath::Abs(cmptrk.c42-mean42[pbin])<0.4*mean42[pbin]) { n42[pbin]++;
      sigma42[pbin]+=(cmptrk.c42-mean42[pbin])*(cmptrk.c42-mean42[pbin]); }
    if(TMath::Abs(cmptrk.c44-mean44[pbin])<0.4*mean44[pbin]) { n44[pbin]++;
      sigma44[pbin]+=(cmptrk.c44-mean44[pbin])*(cmptrk.c44-mean44[pbin]); }
  } // loop on chain entries
 
  for(Int_t l=0; l<fitbins; l++) {
    sigma00[l] = TMath::Sqrt(sigma00[l]/n00[l]);
    sigma11[l] = TMath::Sqrt(sigma11[l]/n11[l]);
    sigma20[l] = TMath::Sqrt(sigma20[l]/n20[l]);
    sigma22[l] = TMath::Sqrt(sigma22[l]/n22[l]);
    sigma31[l] = TMath::Sqrt(sigma31[l]/n31[l]);
    sigma33[l] = TMath::Sqrt(sigma33[l]/n33[l]);
    sigma40[l] = TMath::Sqrt(sigma40[l]/n40[l]);
    sigma42[l] = TMath::Sqrt(sigma42[l]/n42[l]);
    sigma44[l] = TMath::Sqrt(sigma44[l]/n44[l]);
  }
  

  // Fit function
  TF1 *func = new TF1("FitRegFunc",FitRegFunc,0.23,50.,4);
  func->SetParNames("A_meas","A_scatt","B1","B2");

  // line to draw on the plots
  TLine *lin = new TLine(-1,1,1.69,1);
  lin->SetLineStyle(2);
  lin->SetLineWidth(2);

  // matrix used to store fit results
  TMatrixD fitRes(9,4);

  //    --- c00 ---

  // create the canvas
  TCanvas *canv00 = new TCanvas("canv00","c00",0,0,700,900); 
  canv00->Divide(1,2);
  // create the graph for cov matrix
  TGraphErrors *gr00 = new TGraphErrors(fitbins,p,mean00,ep,sigma00);
  TString title00("C(y,y)"); title00.Prepend(part);
  TH2F *frame00 = new TH2F("frame00",title00.Data(),2,0.1,50,2,0,5e-3);
  frame00->SetXTitle("p [GeV/c]");
  canv00->cd(1);  gPad->SetLogx();
  frame00->Draw();
  gr00->Draw("P");
  // Sets initial values for parameters
  func->SetParameters(1.6e-3,1.9e-4,1.5,0.);
  // Fit points in range defined by function
  gr00->Fit("FitRegFunc","R,Q");
  func->GetParameters(fitpar);
  for(Int_t i=0; i<4; i++) fitRes(0,i)=fitpar[i];
  for(Int_t l=0; l<fitbins; l++) {
    rmean[l]  = mean00[l]/FitRegFunc(&p[l],fitpar);
    rsigma[l] = sigma00[l]/FitRegFunc(&p[l],fitpar);
  }
  // create the graph the regularized cov. matrix
  TGraphErrors *gr00reg = new TGraphErrors(fitbins,p,rmean,ep,rsigma);
  TString regtitle00("C(y,y)/(A_meas+A_scatt/p^{B1}/cos^{B2} #lambda"); 
  regtitle00.Prepend(part);
  TH2F *frame00reg = new TH2F("frame00reg",regtitle00.Data(),2,0.1,50,2,0,2);
  frame00reg->SetXTitle("p [GeV/c]");
  canv00->cd(2); gPad->SetLogx();
  frame00reg->Draw();
  gr00reg->Draw("P");
  lin->Draw("same");


  //    --- c11 ---

  // create the canvas
  TCanvas *canv11 = new TCanvas("canv11","c11",0,0,700,900); 
  canv11->Divide(1,2);
  // create the graph for cov matrix
  TGraphErrors *gr11 = new TGraphErrors(fitbins,p,mean11,ep,sigma11);
  TString title11("C(z,z)"); title11.Prepend(part);
  TH2F *frame11 = new TH2F("frame11",title11.Data(),2,0.1,50,2,0,6e-3);
  frame11->SetXTitle("p [GeV/c]");
  canv11->cd(1);  gPad->SetLogx();
  frame11->Draw();
  gr11->Draw("P");
  // Sets initial values for parameters
  func->SetParameters(1.2e-3,8.1e-4,1.,0.);
  // Fit points in range defined by function
  gr11->Fit("FitRegFunc","R,Q");
  func->GetParameters(fitpar);
  for(Int_t i=0; i<4; i++) fitRes(1,i)=fitpar[i];
  for(Int_t l=0; l<fitbins; l++) {
    rmean[l]  = mean11[l]/FitRegFunc(&p[l],fitpar);
    rsigma[l] = sigma11[l]/FitRegFunc(&p[l],fitpar);
  }
  // create the graph the regularized cov. matrix
  TGraphErrors *gr11reg = new TGraphErrors(fitbins,p,rmean,ep,rsigma);
  TString regtitle11("C(z,z)/(A_meas+A_scatt/p^{B1}/cos^{B2} #lambda"); 
  regtitle11.Prepend(part);
  TH2F *frame11reg = new TH2F("frame11reg",regtitle11.Data(),2,0.1,50,2,0,2);
  frame11reg->SetXTitle("p [GeV/c]");
  canv11->cd(2); gPad->SetLogx();
  frame11reg->Draw();
  gr11reg->Draw("P");
  lin->Draw("same");


  //    --- c20 ---

  // create the canvas
  TCanvas *canv20 = new TCanvas("canv20","c20",0,0,700,900); 
  canv20->Divide(1,2);
  // create the graph for cov matrix
  TGraphErrors *gr20 = new TGraphErrors(fitbins,p,mean20,ep,sigma20);
  TString title20("C(#eta, y)"); title20.Prepend(part);
  TH2F *frame20 = new TH2F("frame20",title20.Data(),2,0.1,50,2,0,2.5e-4);
  frame20->SetXTitle("p [GeV/c]");
  canv20->cd(1);  gPad->SetLogx();
  frame20->Draw();
  gr20->Draw("P");
  // Sets initial values for parameters
  func->SetParameters(7.3e-5,1.2e-5,1.5,0.);
  // Fit points in range defined by function
  gr20->Fit("FitRegFunc","R,Q");
  func->GetParameters(fitpar);
  for(Int_t i=0; i<4; i++) fitRes(2,i)=fitpar[i];
  for(Int_t l=0; l<fitbins; l++) {
    rmean[l]  = mean20[l]/FitRegFunc(&p[l],fitpar);
    rsigma[l] = sigma20[l]/FitRegFunc(&p[l],fitpar);
  }
  // create the graph the regularized cov. matrix
  TGraphErrors *gr20reg = new TGraphErrors(fitbins,p,rmean,ep,rsigma);
  TString regtitle20("C(#eta, y)/(A_meas+A_scatt/p^{B1}/cos^{B2} #lambda"); 
  regtitle20.Prepend(part);
  TH2F *frame20reg = new TH2F("frame20reg",regtitle20.Data(),2,0.1,50,2,0,2);
  frame20reg->SetXTitle("p [GeV/c]");
  canv20->cd(2); gPad->SetLogx();
  frame20reg->Draw();
  gr20reg->Draw("P");
  lin->Draw("same");


  //    --- c22 ---

  // create the canvas
  TCanvas *canv22 = new TCanvas("canv22","c22",0,0,700,900); 
  canv22->Divide(1,2);
  // create the graph for cov matrix
  TGraphErrors *gr22 = new TGraphErrors(fitbins,p,mean22,ep,sigma22);
  TString title22("C(#eta, #eta)"); title22.Prepend(part);
  TH2F *frame22 = new TH2F("frame22",title22.Data(),2,0.1,50,2,0,3e-5);
  frame22->SetXTitle("p [GeV/c]");
  canv22->cd(1);  gPad->SetLogx();
  frame22->Draw();
  gr22->Draw("P");
  // Sets initial values for parameters
  func->SetParameters(5.2e-6,1.1e-6,2.,1.);
  // Fit points in range defined by function
  gr22->Fit("FitRegFunc","R,Q");
  func->GetParameters(fitpar);
  for(Int_t i=0; i<4; i++) fitRes(3,i)=fitpar[i];
  for(Int_t l=0; l<fitbins; l++) {
    rmean[l]  = mean22[l]/FitRegFunc(&p[l],fitpar);
    rsigma[l] = sigma22[l]/FitRegFunc(&p[l],fitpar);
  }
  // create the graph the regularized cov. matrix
  TGraphErrors *gr22reg = new TGraphErrors(fitbins,p,rmean,ep,rsigma);
  TString regtitle22("C(#eta, #eta)/(A_meas+A_scatt/p^{B1}/cos^{B2} #lambda"); 
  regtitle22.Prepend(part);
  TH2F *frame22reg = new TH2F("frame22reg",regtitle22.Data(),2,0.1,50,2,0,2);
  frame22reg->SetXTitle("p [GeV/c]");
  canv22->cd(2); gPad->SetLogx();
  frame22reg->Draw();
  gr22reg->Draw("P");
  lin->Draw("same");


  //    --- c31 ---

  // create the canvas
  TCanvas *canv31 = new TCanvas("canv31","c31",0,0,700,900); 
  canv31->Divide(1,2);
  // create the graph for cov matrix
  TGraphErrors *gr31 = new TGraphErrors(fitbins,p,mean31,ep,sigma31);
  TString title31("C(tg #lambda,z)"); title31.Prepend(part);
  TH2F *frame31 = new TH2F("frame31",title31.Data(),2,0.1,50,2,-2e-4,0);
  frame31->SetXTitle("p [GeV/c]");
  canv31->cd(1);  gPad->SetLogx();
  frame31->Draw();
  gr31->Draw("P");
  // Sets initial values for parameters
  func->SetParameters(-1.2e-5,-1.2e-5,1.5,3.);
  // Fit points in range defined by function
  gr31->Fit("FitRegFunc","R,Q");
  func->GetParameters(fitpar);
  for(Int_t i=0; i<4; i++) fitRes(4,i)=fitpar[i];
  for(Int_t l=0; l<fitbins; l++) {
    rmean[l]  = mean31[l]/FitRegFunc(&p[l],fitpar);
    rsigma[l] = -sigma31[l]/FitRegFunc(&p[l],fitpar);
  }
  // create the graph the regularized cov. matrix
  TGraphErrors *gr31reg = new TGraphErrors(fitbins,p,rmean,ep,rsigma);
  TString regtitle31("C(tg #lambda,z)/(A_meas+A_scatt/p^{B1}/cos^{B2} #lambda"); 
  regtitle31.Prepend(part);
  TH2F *frame31reg = new TH2F("frame31reg",regtitle31.Data(),2,0.1,50,2,0,2);
  frame31reg->SetXTitle("p [GeV/c]");
  canv31->cd(2); gPad->SetLogx();
  frame31reg->Draw();
  gr31reg->Draw("P");
  lin->Draw("same");


  //    --- c33 ---

  // create the canvas
  TCanvas *canv33 = new TCanvas("canv33","c33",0,0,700,900); 
  canv33->Divide(1,2);
  // create the graph for cov matrix
  TGraphErrors *gr33 = new TGraphErrors(fitbins,p,mean33,ep,sigma33);
  TString title33("C(tg #lambda,tg #lambda)"); title33.Prepend(part);
  TH2F *frame33 = new TH2F("frame33",title33.Data(),2,0.1,50,2,0,1e-5);
  frame33->SetXTitle("p [GeV/c]");
  canv33->cd(1);  gPad->SetLogx();
  frame33->Draw();
  gr33->Draw("P");
  // Sets initial values for parameters
  func->SetParameters(1.3e-7,4.6e-7,1.7,4.);
  // Fit points in range defined by function
  gr33->Fit("FitRegFunc","R,Q");
  func->GetParameters(fitpar);
  for(Int_t i=0; i<4; i++) fitRes(5,i)=fitpar[i];
  for(Int_t l=0; l<fitbins; l++) {
    rmean[l]  = mean33[l]/FitRegFunc(&p[l],fitpar);
    rsigma[l] = sigma33[l]/FitRegFunc(&p[l],fitpar);
  }
  // create the graph the regularized cov. matrix
  TGraphErrors *gr33reg = new TGraphErrors(fitbins,p,rmean,ep,rsigma);
  TString regtitle33("C(tg #lambda,tg #lambda)/(A_meas+A_scatt/p^{B1}/cos^{B2} #lambda"); 
  regtitle33.Prepend(part);
  TH2F *frame33reg = new TH2F("frame33reg",regtitle33.Data(),2,0.1,50,2,0,2);
  frame33reg->SetXTitle("p [GeV/c]");
  canv33->cd(2); gPad->SetLogx();
  frame33reg->Draw();
  gr33reg->Draw("P");
  lin->Draw("same");


  //    --- c40 ---

  // create the canvas
  TCanvas *canv40 = new TCanvas("canv40","c40",0,0,700,900); 
  canv40->Divide(1,2);
  // create the graph for cov matrix
  TGraphErrors *gr40 = new TGraphErrors(fitbins,p,mean40,ep,sigma40);
  TString title40("C(C,y)"); title40.Prepend(part);
  TH2F *frame40 = new TH2F("frame40",title40.Data(),2,0.1,50,2,0,1e-6);
  frame40->SetXTitle("p [GeV/c]");
  canv40->cd(1);  gPad->SetLogx();
  frame40->Draw();
  gr40->Draw("P");
  // Sets initial values for parameters
  func->SetParameters(4.e-7,4.4e-8,1.5,0.);
  // Fit points in range defined by function
  gr40->Fit("FitRegFunc","R,Q");
  func->GetParameters(fitpar);
  for(Int_t i=0; i<4; i++) fitRes(6,i)=fitpar[i];
  for(Int_t l=0; l<fitbins; l++) {
    rmean[l]  = mean40[l]/FitRegFunc(&p[l],fitpar);
    rsigma[l] = sigma40[l]/FitRegFunc(&p[l],fitpar);
  }
  // create the graph the regularized cov. matrix
  TGraphErrors *gr40reg = new TGraphErrors(fitbins,p,rmean,ep,rsigma);
  TString regtitle40("C(C,y)/(A_meas+A_scatt/p^{B1}/cos^{B2} #lambda"); 
  regtitle40.Prepend(part);
  TH2F *frame40reg = new TH2F("frame40reg",regtitle40.Data(),2,0.1,50,2,0,2);
  frame40reg->SetXTitle("p [GeV/c]");
  canv40->cd(2); gPad->SetLogx();
  frame40reg->Draw();
  gr40reg->Draw("P");
  lin->Draw("same");


  //    --- c42 ---

  // create the canvas
  TCanvas *canv42 = new TCanvas("canv42","c42",0,0,700,900); 
  canv42->Divide(1,2);
  // create the graph for cov matrix
  TGraphErrors *gr42 = new TGraphErrors(fitbins,p,mean42,ep,sigma42);
  TString title42("C(C, #eta)"); title42.Prepend(part);
  TH2F *frame42 = new TH2F("frame42",title42.Data(),2,0.1,50,2,0,2.2e-7);
  frame42->SetXTitle("p [GeV/c]");
  canv42->cd(1);  gPad->SetLogx();
  frame42->Draw();
  gr42->Draw("P");
  // Sets initial values for parameters
  func->SetParameters(3.e-8,8.2e-9,2.,1.);
  // Fit points in range defined by function
  gr42->Fit("FitRegFunc","R,Q");
  func->GetParameters(fitpar);
  for(Int_t i=0; i<4; i++) fitRes(7,i)=fitpar[i];
  for(Int_t l=0; l<fitbins; l++) {
    rmean[l]  = mean42[l]/FitRegFunc(&p[l],fitpar);
    rsigma[l] = sigma42[l]/FitRegFunc(&p[l],fitpar);
  }
  // create the graph the regularized cov. matrix
  TGraphErrors *gr42reg = new TGraphErrors(fitbins,p,rmean,ep,rsigma);
  TString regtitle42("C(C, #eta)/(A_meas+A_scatt/p^{B1}/cos^{B2} #lambda"); 
  regtitle42.Prepend(part);
  TH2F *frame42reg = new TH2F("frame42reg",regtitle42.Data(),2,0.1,50,2,0,2);
  frame42reg->SetXTitle("p [GeV/c]");
  canv42->cd(2); gPad->SetLogx();
  frame42reg->Draw();
  gr42reg->Draw("P");
  lin->Draw("same");


  //    --- c44 ---

  // create the canvas
  TCanvas *canv44 = new TCanvas("canv44","c44",0,0,700,900); 
  canv44->Divide(1,2);
  // create the graph for cov matrix
  TGraphErrors *gr44 = new TGraphErrors(fitbins,p,mean44,ep,sigma44);
  TString title44("C(C,C)"); title44.Prepend(part);
  TH2F *frame44 = new TH2F("frame44",title44.Data(),2,0.1,50,2,0,2e-9);
  frame44->SetXTitle("p [GeV/c]");
  canv44->cd(1);  gPad->SetLogx();
  frame44->Draw();
  gr44->Draw("P");
  // Sets initial values for parameters
  func->SetParameters(1.8e-10,5.8e-11,2.,3.);
  // Fit points in range defined by function
  gr44->Fit("FitRegFunc","R,Q");
  func->GetParameters(fitpar);
  for(Int_t i=0; i<4; i++) fitRes(8,i)=fitpar[i];
  for(Int_t l=0; l<fitbins; l++) {
    rmean[l]  = mean44[l]/FitRegFunc(&p[l],fitpar);
    rsigma[l] = sigma44[l]/FitRegFunc(&p[l],fitpar);
  }
  // create the graph the regularized cov. matrix
  TGraphErrors *gr44reg = new TGraphErrors(fitbins,p,rmean,ep,rsigma);
  TString regtitle44("C(C,C)/(A_meas+A_scatt/p^{B1}/cos^{B2} #lambda"); 
  regtitle44.Prepend(part);
  TH2F *frame44reg = new TH2F("frame44reg",regtitle44.Data(),2,0.1,50,2,0,2);
  frame44reg->SetXTitle("p [GeV/c]");
  canv44->cd(2); gPad->SetLogx();
  frame44reg->Draw();
  gr44reg->Draw("P");
  lin->Draw("same");




  switch (pdg) {
  case 211:
    new(&fRegParPi) TMatrixD(fitRes);
    break;
  case 321:
    new(&fRegParKa) TMatrixD(fitRes);
    break;
  case 11:
    new(&fRegParEl) TMatrixD(fitRes);
    break;
  }

  // write fit parameters to file
  WriteRegParams(outName,pdg);

  delete [] n;
  delete [] n00;
  delete [] n11;
  delete [] n20;
  delete [] n22;
  delete [] n31;
  delete [] n33;
  delete [] n40;
  delete [] n42;
  delete [] n44;
  delete [] p;
  delete [] ep;
  delete [] mean00;
  delete [] mean11;
  delete [] mean20;
  delete [] mean22;
  delete [] mean31;
  delete [] mean33;
  delete [] mean40;
  delete [] mean42;
  delete [] mean44;
  delete [] sigma00;
  delete [] sigma11;
  delete [] sigma20;
  delete [] sigma22;
  delete [] sigma31;
  delete [] sigma33;
  delete [] sigma40;
  delete [] sigma42;
  delete [] sigma44;
  delete [] rmean;
  delete [] rsigma;

  return;
}
//-----------------------------------------------------------------------------
Bool_t AliTPCtrackerParam::SelectedTrack(Double_t pt,Double_t eta) const {
//-----------------------------------------------------------------------------
// This function makes a selection according to TPC tracking efficiency
//-----------------------------------------------------------------------------

  Double_t eff=0.;
 
  eff = fEff->GetValueAt(pt,eta);

  if(gRandom->Rndm() < eff) return kTRUE;

  return kFALSE;
}
//-----------------------------------------------------------------------------
void AliTPCtrackerParam::SetParticle(Int_t pdg) {
//-----------------------------------------------------------------------------
// This function set efficiencies, pulls, etc... for the particle
// specie of the current particle
//-----------------------------------------------------------------------------

  switch (pdg) {
  case 211:
    fDBgrid = &fDBgridPi;
    fEff    = &fEffPi;
    fPulls  =  fPullsPi;
    fRegPar = &fRegParPi;
    fdEdxMean = &fdEdxMeanPi;
    fdEdxRMS  = &fdEdxRMSPi;
    break;
  case 321:
    fDBgrid = &fDBgridKa;
    fEff    = &fEffKa;
    fPulls  =  fPullsKa;
    fRegPar = &fRegParKa;
    fdEdxMean = &fdEdxMeanKa;
    fdEdxRMS  = &fdEdxRMSKa;
    break;
  case 2212:
    fDBgrid = &fDBgridPi;
    fEff    = &fEffPr;
    fPulls  =  fPullsPi;
    fRegPar = &fRegParPi;
    fdEdxMean = &fdEdxMeanPr;
    fdEdxRMS  = &fdEdxRMSPr;
    break;
  case 11:
    fDBgrid = &fDBgridEl;
    fEff    = &fEffEl;
    fPulls  =  fPullsEl;
    fRegPar = &fRegParEl;
    fdEdxMean = &fdEdxMeanEl;
    fdEdxRMS  = &fdEdxRMSEl;
    break;
  case 13:
    fDBgrid = &fDBgridPi;
    fEff    = &fEffMu;
    fPulls  =  fPullsPi;
    fRegPar = &fRegParPi;
    fdEdxMean = &fdEdxMeanPi;
    fdEdxRMS  = &fdEdxRMSPi;
    break;
  }

  return;
}
//-----------------------------------------------------------------------------
void AliTPCtrackerParam::SmearTrack(Double_t *xx,Double_t *xxsm,TMatrixD cov)
  const {
//-----------------------------------------------------------------------------
// This function smears track parameters according to streched cov. matrix
//-----------------------------------------------------------------------------
  AliGausCorr *corgen = new AliGausCorr(cov,5);
  TArrayD corr(5);
  corgen->GetGaussN(corr);
  delete corgen;
  corgen = 0;

  for(Int_t l=0;l<5;l++) {
    xxsm[l] = xx[l]+corr[l];
  }

  return;
}
//-----------------------------------------------------------------------------
Int_t AliTPCtrackerParam::WritedEdx(const Char_t *outName,Int_t pdg) {
//-----------------------------------------------------------------------------
// This function writes the dEdx parameters to the DB
//-----------------------------------------------------------------------------

  Option_t *opt;
  Char_t *dirName="Pions";
  Char_t *meanName="dEdxMeanPi";
  Char_t *RMSName="dEdxRMSPi";

  SetParticle(pdg);

  if(gSystem->AccessPathName(outName,kFileExists))  { opt="recreate"; 
  } else { opt="update"; }

  switch (pdg) {
  case 211:
    dirName="Pions";
    meanName="dEdxMeanPi";
    RMSName="dEdxRMSPi";
    break;
  case 321:
    dirName="Kaons";
    meanName="dEdxMeanKa";
    RMSName="dEdxRMSKa";
    break;
  case 2212:
    dirName="Protons";
    meanName="dEdxMeanPr";
    RMSName="dEdxRMSPr";
    break;
  case 11:
    dirName="Electrons";
    meanName="dEdxMeanEl";
    RMSName="dEdxRMSEl";
    break;
  }

  TFile *outFile = new TFile(outName,opt);
  if(!gDirectory->cd("/dEdx")) {
    outFile->mkdir("dEdx");
    gDirectory->cd("/dEdx"); 
  }
  TDirectory *dir2 = gDirectory->mkdir(dirName);
  dir2->cd();
  fdEdxMean->SetName(meanName); fdEdxMean->Write();
  fdEdxRMS->SetName(RMSName);  fdEdxRMS->Write();

  outFile->Close();
  delete outFile;


  return 1;
}
//-----------------------------------------------------------------------------
Int_t AliTPCtrackerParam::WriteEffs(const Char_t *outName) {
//-----------------------------------------------------------------------------
// This function writes the TPC efficiencies to the DB
//-----------------------------------------------------------------------------


  Option_t *opt;

  if(gSystem->AccessPathName(outName,kFileExists))  { opt="recreate"; 
  } else { opt="update"; }

  TFile *outFile = new TFile(outName,opt);

  outFile->mkdir("Efficiencies");
  gDirectory->cd("/Efficiencies");
  gDirectory->mkdir("Pions");
  gDirectory->mkdir("Kaons");
  gDirectory->mkdir("Protons");
  gDirectory->mkdir("Electrons");
  gDirectory->mkdir("Muons");

  gDirectory->cd("/Efficiencies/Pions");
  fEffPi.SetName("EffPi");
  fEffPi.Write();
  gDirectory->cd("/Efficiencies/Kaons");
  fEffKa.SetName("EffKa");
  fEffKa.Write();
  gDirectory->cd("/Efficiencies/Protons");
  fEffPr.SetName("EffPr");
  fEffPr.Write();
  gDirectory->cd("/Efficiencies/Electrons");
  fEffEl.SetName("EffEl");
  fEffEl.Write();
  gDirectory->cd("/Efficiencies/Muons");
  fEffMu.SetName("EffMu");
  fEffMu.Write();

  outFile->Close();

  delete outFile;

  return 1;
}
//-----------------------------------------------------------------------------
Int_t AliTPCtrackerParam::WritePulls(const Char_t *outName) {
//-----------------------------------------------------------------------------
// This function writes the pulls to the DB
//-----------------------------------------------------------------------------

  Option_t *opt;

  if(gSystem->AccessPathName(outName,kFileExists))  { opt="recreate"; 
  } else { opt="update"; }

  TFile *outFile = new TFile(outName,opt);

  outFile->mkdir("Pulls");
  gDirectory->cd("/Pulls");
  gDirectory->mkdir("Pions");
  gDirectory->mkdir("Kaons");
  gDirectory->mkdir("Electrons");

  for(Int_t i=0;i<5;i++) {
    TString pi("PullsPi_"); pi+=i;
    TString ka("PullsKa_"); ka+=i;
    TString el("PullsEl_"); el+=i;
    fPullsPi[i].SetName(pi.Data());
    fPullsKa[i].SetName(ka.Data());
    fPullsEl[i].SetName(el.Data());
    gDirectory->cd("/Pulls/Pions");
    fPullsPi[i].Write();
    gDirectory->cd("/Pulls/Kaons");
    fPullsKa[i].Write();
    gDirectory->cd("/Pulls/Electrons");
    fPullsEl[i].Write();
  }
  outFile->Close();
  delete outFile;

  return 1;
}
//-----------------------------------------------------------------------------
Int_t AliTPCtrackerParam::WriteRegParams(const Char_t *outName,Int_t pdg) {
//-----------------------------------------------------------------------------
// This function writes the regularization parameters to the DB
//-----------------------------------------------------------------------------

  Option_t *opt;
  Char_t *dirName="Pions";
  Char_t *keyName="RegPions";

  SetParticle(pdg);

  if(gSystem->AccessPathName(outName,kFileExists))  { opt="recreate"; 
  } else { opt="update"; }

  switch (pdg) {
  case 211:
    dirName="Pions";
    keyName="RegPions";
    break;
  case 321:
    dirName="Kaons";
    keyName="RegKaons";
    break;
  case 11:
    dirName="Electrons";
    keyName="RegElectrons";
    break;
  }

  TFile *outFile = new TFile(outName,opt);
  if(!gDirectory->cd("/RegParams")) {
    outFile->mkdir("RegParams");
    gDirectory->cd("/RegParams"); 
  }
  TDirectory *dir2 = gDirectory->mkdir(dirName);
  dir2->cd();
  fRegPar->Write(keyName);

  outFile->Close();
  delete outFile;


  return 1;
}









