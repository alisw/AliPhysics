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
 * Alice Internal Note 2003-011                                           *
 *                                                                        *
 * Test macro is: AliBarrelRec_TPCparam.C                                 *   
 *                                                                        *
 * 2002/10/01: Introduction of the database for pp collisions (B=0.4 T)   *
 * - Everything done separately for pions, kaons, protons, electrons and  *
 *   muons.                                                               *
 * - Now (only for pp) the tracks are built from the AliTrackReferences   *
 *   which contain the position and momentum of all tracks at R = 83 cm;  *
 *   This eliminates the loss of tracks due the dead zone of the TPC      *
 *   where the 1st hit is not created.                                    *
 * - In AliBarrelRec_TPCparam.C there many possible ways of determining   *
 *   the z coordinate of the primary vertex in pp events (from pixels,    *
 *   from ITS tracks, smearing according to resolution given by tracks.   *
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
 * 2006/03/16: Adapted to ESD input/output                                *
 *                                                                        *
 *  Origin: Andrea Dainese, Padova - e-mail: andrea.dainese@pd.infn.it    *
 *       adapted to ESD output by Marcello Lunardon, Padova               *
 **************************************************************************/
//  *
// This is a dummy comment
//
//
// *
//------- Root headers --------
#include <Riostream.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TF1.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH2F.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMatrixD.h>
#include <TParticle.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TFile.h>
//------ AliRoot headers ------
#include "AliGausCorr.h"
#include "AliTracker.h"
#include "AliMC.h"
#include "AliMagF.h"
#include "AliRun.h"
#include "AliRunLoader.h"
#include "AliTPC.h"
#include "AliTPCParamSR.h"
#include "AliTPCkineGrid.h"
#include "AliTPCtrack.h"
#include "AliTPCtrackerParam.h"
#include "AliTrackReference.h"
//-----------------------------

Double_t RegFunc(Double_t *x,Double_t *par) {
// This is the function used to regularize the covariance matrix
  Double_t value = par[0]+par[1]/TMath::Power(x[0],par[2]);
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
AliTPCtrackerParam::AliTPCtrackerParam(Int_t kcoll, Double_t kBz,
				       const char* evfoldname):TObject(),
    fEvFolderName(evfoldname),
    fBz(kBz),
    fColl(kcoll),
    fSelAndSmear(kTRUE),
    fDBfileName(""),
    fTrack(),
    fCovTree(0),
    fDBgrid(0),
    fDBgridPi(),
    fDBgridKa(),
    fDBgridPr(),
    fDBgridEl(),
    fDBgridMu(),
    fEff(0),
    fEffPi(),
    fEffKa(),
    fEffPr(),
    fEffEl(),
    fEffMu(),
    fPulls(0),
    fRegPar(0),
    fRegParPi(),
    fRegParKa(),
    fRegParPr(),
    fRegParEl(),
    fRegParMu(),
    fdEdxMean(0),
    fdEdxMeanPi(),
    fdEdxMeanKa(),
    fdEdxMeanPr(),
    fdEdxMeanEl(),
    fdEdxMeanMu(),
    fdEdxRMS(0),
    fdEdxRMSPi(),
    fdEdxRMSKa(),
    fdEdxRMSPr(),
    fdEdxRMSEl(),
    fdEdxRMSMu() 
{
//-----------------------------------------------------------------------------
// This is the class conctructor 
//-----------------------------------------------------------------------------

  // fBz = kBz;             // value of the z component of L3 field (Tesla)
  //  fColl = kcoll;         // collision code (0: PbPb6000; 1: pp)
  //  fSelAndSmear = kTRUE; // by default selection and smearing are done

  if(fBz!=0.4 && fBz!=0.5) {
    Fatal("AliTPCtrackerParam","AliTPCtrackerParam::AliTPCtrackerParam:  Invalid field!\n      Available:  0.4 or 0.5");
  }
  if(fColl!=0 && fColl!=1) {
    Fatal("AliTPCtrackerParam","AliTPCtrackerParam::AliTPCtrackerParam:  Invalid collision!\n      Available:  0   ->   PbPb6000\n                  1   ->   pp"); 
  }

  fDBfileName = gSystem->Getenv("ALICE_ROOT");  
  fDBfileName.Append("/TPC/CovMatrixDB_");
  //fDBfileName = "CovMatrixDB_";
  if(fColl==0) fDBfileName.Append("PbPb6000");
  if(fColl==1) fDBfileName.Append("pp");
  if(fBz==0.4) fDBfileName.Append("_B0.4T.root");
  // use same DB for 0.4 and 0.5 T; for 0.5 T, correction done in CookTrack()
  if(fBz==0.5) fDBfileName.Append("_B0.4T.root");
}
//-----------------------------------------------------------------------------
AliTPCtrackerParam::~AliTPCtrackerParam() {}
//____________________________________________________________________________
AliTPCtrackerParam::AliTPCtrackerParam( const AliTPCtrackerParam& p)
    :TObject(p),
    fEvFolderName(""),
    fBz(0.),
    fColl(0),
    fSelAndSmear(0),
    fDBfileName(""),
    fTrack(),
    fCovTree(0),
    fDBgrid(0),
    fDBgridPi(),
    fDBgridKa(),
    fDBgridPr(),
    fDBgridEl(),
    fDBgridMu(),
    fEff(0),
    fEffPi(),
    fEffKa(),
    fEffPr(),
    fEffEl(),
    fEffMu(),
    fPulls(0),
    fRegPar(0),
    fRegParPi(),
    fRegParKa(),
    fRegParPr(),
    fRegParEl(),
    fRegParMu(),
    fdEdxMean(0),
    fdEdxMeanPi(),
    fdEdxMeanKa(),
    fdEdxMeanPr(),
    fdEdxMeanEl(),
    fdEdxMeanMu(),
    fdEdxRMS(0),
    fdEdxRMSPi(),
    fdEdxRMSKa(),
    fdEdxRMSPr(),
    fdEdxRMSEl(),
    fdEdxRMSMu() 
{
  // dummy copy constructor
}
//----------------------------------------------------------------------------
AliTPCtrackerParam::AliTPCseedGeant::AliTPCseedGeant(
		    Double_t x,Double_t y,Double_t z,
		    Double_t px,Double_t py,Double_t pz,
		    Int_t lab)
                    :TObject(),
      fXg(x),
      fYg(y),
      fZg(z),
      fPx(px),
      fPy(py),
      fPz(pz),
      fAlpha(0.),
      fLabel(lab),
      fSector(0)
 
{
//----------------------------------------------------------------------------
// Constructor of the geant seeds
//----------------------------------------------------------------------------

      Double_t a = TMath::ATan2(y,x)*180./TMath::Pi();
      if(a<0) a += 360.;
      fSector = (Int_t)(a/20.);
      fAlpha = 10.+20.*fSector;
      fAlpha /= 180.;
      fAlpha *= TMath::Pi();
}
//-----------------------------------------------------------------------------
Int_t AliTPCtrackerParam::Init() {
//-----------------------------------------------------------------------------
// This function reads the DB from file
//-----------------------------------------------------------------------------

  if(fSelAndSmear) {
    printf("+++\n+++ Reading DataBase from:%s\n+++\n+++\n",fDBfileName.Data()); 
    // Read paramters from DB file
    if(!ReadAllData(fDBfileName.Data())) {
      printf("AliTPCtrackerParam::BuildTPCtracks: \
             Could not read data from DB\n\n"); return 1; 
    }
    
  } else printf("\n ! Creating ALL TRUE tracks at TPC inner radius !\n\n");


  // Check if value in the galice file is equal to selected one (fBz)
  AliMagF *fiel = (AliMagF*)gAlice->Field();
  Double_t fieval=TMath::Abs((Double_t)fiel->SolenoidField()/10.);
  printf("Magnetic field is %6.2f Tesla\n",fieval);
  if(fBz!=fieval) {
    printf("AliTPCtrackerParam::BuildTPCtracks:  Invalid field!");
    printf("Field selected is: %f T\n",fBz);
    printf("Field found on file is: %f T\n",fieval);
    return 1;
  }

  // Set the conversion constant between curvature and Pt
  AliTracker::SetFieldMap(fiel,kTRUE);

  return 0;
}
//-----------------------------------------------------------------------------
Int_t AliTPCtrackerParam::BuildTPCtracks(AliESD *event) {
//-----------------------------------------------------------------------------
// This function creates the TPC parameterized tracks and writes them
// as AliESDtrack objects in the ESD event
//-----------------------------------------------------------------------------


  if(!event) return -1;

  AliRunLoader* rl = AliRunLoader::GetRunLoader(fEvFolderName);
  if (rl == 0x0) {
    Error("BuildTPCtracks","Can not get Run Loader from event folder named %s",
	  fEvFolderName.Data());
    return 2;
  }
  AliLoader* tpcloader = rl->GetLoader("TPCLoader");
  if (tpcloader == 0x0) {
    Error("BuildTPCtracks","Can not get TPC Loader from Run Loader.");
    return 3;
  }
  
  // Get gAlice object from file  
  if(!gAlice) rl->LoadgAlice();
  //rl->LoadHeader();
  rl->LoadKinematics();
  tpcloader->LoadHits("read");
  
  if(!(gAlice=rl->GetAliRun())) {
    printf("Cannot get gAlice from Run Loader !\n");
    return 1;
  }

  // Get TPC detector 
  AliTPC *tpc=(AliTPC*)gAlice->GetDetector("TPC");

  rl->CdGAFile();
  AliTPCParam *digp=(AliTPCParam*)gDirectory->Get("75x40_100x60");
  if(digp){
    delete digp;
    digp = new AliTPCParamSR();
  }
  else digp=(AliTPCParam*)gDirectory->Get("75x40_100x60_150x60");
  
  if(!digp) { cerr<<"TPC parameters have not been found !\n"; return 1; }
  tpc->SetParam(digp);

  TParticle       *part=0;
  AliTPCseedGeant *seed=0;
  AliTPCtrack     *tpctrack=0;
  Double_t     sPt,sEta;
  Int_t        bin,label,pdg,charge;
  Int_t        tracks;
  Int_t        nParticles,nSeeds,arrentr;
  //Int_t nSel=0,nAcc=0;

  Int_t evt=event->GetEventNumber();
  
  tracks=0;

  // array for TPC tracks
  TObjArray tArray(20000);
  
  // array for TPC seeds with geant info
  TObjArray sArray(20000);
  
  // get the particles stack
  nParticles = (Int_t)gAlice->GetEvent(evt);
    
  Bool_t   *done     = new Bool_t[nParticles];
  Int_t    *pdgCodes = new Int_t[nParticles];
  
  // loop on particles and store pdg codes
  for(Int_t l=0; l<nParticles; l++) {
    part        = (TParticle*)gAlice->GetMCApp()->Particle(l);
    pdgCodes[l] = part->GetPdgCode();
    done[l]     = kFALSE;
  }

  printf("+++ Number of particles:  %d\n",nParticles);

  // Create the seeds for the TPC tracks at the inner radius of TPC
  if(fColl==0) {
    // Get TreeH with hits
    TTree *th = tpcloader->TreeH(); 
    MakeSeedsFromHits(tpc,th,sArray);
  } else {
    // Get TreeTR with track references
    rl->LoadTrackRefs();
    TTree *ttr = rl->TreeTR();
    if (!ttr) return -3;
    MakeSeedsFromRefs(ttr,sArray);
  }

  nSeeds = sArray.GetEntries();
  printf("+++ Number of seeds: %d\n",nSeeds);
    
  // loop over entries in sArray
  for(Int_t l=0; l<nSeeds; l++) {
    //if(l%1==0) printf("  --- Processing seed %d of %d ---\n",l,nSeeds);
    
    seed = (AliTPCseedGeant*)sArray.At(l);
    
    // Get track label
    label = seed->GetLabel();
    
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
    
    sPt  = seed->GetPt();
    sEta = seed->GetEta();
    
    // Apply selection according to TPC efficiency
    //if(TMath::Abs(pdg)==211) nAcc++;
    if(fSelAndSmear && !SelectedTrack(sPt,sEta)) continue; 
    //if(TMath::Abs(pdg)==211) nSel++;

    // create AliTPCtrack object
    BuildTrack(seed,charge);

    if(fSelAndSmear) {
      bin = fDBgrid->GetBin(sPt,sEta);
      switch (pdg) {
      case 211:
	//fCovTree = &(fCovTreePi[bin]);
	fCovTree = fCovTreePi[bin];
	break;
      case 321:
	//fCovTree = &(fCovTreeKa[bin]);
	fCovTree = fCovTreeKa[bin];
	break;
      case 2212:
	//fCovTree = &(fCovTreePr[bin]);
	fCovTree = fCovTreePr[bin];
	break;
      case 11:
	//fCovTree = &(fCovTreeEl[bin]);
	fCovTree = fCovTreeEl[bin];
	break;
      case 13:
	//fCovTree = &(fCovTreeMu[bin]);
	fCovTree = fCovTreeMu[bin];
	break;
      }
      // deal with covariance matrix and smearing of parameters
      CookTrack(sPt,sEta);

      // assign the track a dE/dx and make a rough PID
      CookdEdx(sPt,sEta);
    }
    
    // put track in array
    AliTPCtrack *iotrack = new AliTPCtrack(fTrack);
    iotrack->SetLabel(label);
    tArray.AddLast(iotrack);
    // Mark track as "done" and register the pdg code
    done[label] = kTRUE; 
    tracks++;
    
  } // loop over entries in sArray
  
  // sort array with TPC tracks (decreasing pT)
  tArray.Sort();
  
  // convert to AliESDtrack and write to AliESD
  arrentr = tArray.GetEntriesFast(); 
  Int_t idx;
  Double_t wgts[5];
  for(Int_t l=0; l<arrentr; l++) {
    tpctrack=(AliTPCtrack*)tArray.UncheckedAt(l);
    AliESDtrack ioESDtrack;
    ioESDtrack.UpdateTrackParams(tpctrack,AliESDtrack::kTPCin);
    // rough PID
    wgts[0]=0.; wgts[1]=0.; wgts[2]=0.; wgts[3]=0.; wgts[4]=0.;
    if(TMath::Abs(tpctrack->GetMass()-0.9)<0.1) {
      idx = 4; // proton
    } else if(TMath::Abs(tpctrack->GetMass()-0.5)<0.1) {
      idx = 3; // kaon
    } else {
      idx = 2; // pion
    }
    wgts[idx] = 1.;
    ioESDtrack.SetESDpid(wgts);
    event->AddTrack(&ioESDtrack);
  }
  
  
  delete [] done;
  delete [] pdgCodes;
  
  printf("+++ Number of TPC tracks: %d\n",tracks);
  //cerr<<"Average Eff: "<<(Float_t)nSel/nAcc<<endl;
  
  sArray.Delete();
  tArray.Delete();
  
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

  const char *part="PIONS";
  Double_t ymax=500.;

  /*
  // create a chain with compared tracks
  TChain *cmptrkchain = new ("cmptrktree");
  cmptrkchain.Add("CovMatrix_AllEvts.root");
  //cmptrkchain.Add("CovMatrix_AllEvts_1.root");
  //cmptrkchain.Add("CovMatrix_AllEvts_2.root");
  //cmptrkchain.Add("CovMatrix_AllEvts_3.root");
  */


  TFile *infile = TFile::Open("CovMatrix_AllEvts.root");
  TTree *cmptrktree = (TTree*)infile->Get("cmptrktree");

  COMPTRACK cmptrk; 
  cmptrktree->SetBranchAddress("comptracks",&cmptrk);
  Int_t entries = (Int_t)cmptrktree->GetEntries(); 
  cerr<<" Number of entries: "<<entries<<endl;

  InitializeKineGrid("DB");
  InitializeKineGrid("dEdx");

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
  case 13:
    part = "MUONS";
    ymax = 100.;
    break;
  }

  SetParticle(pdg);

  const Int_t knTotBins = fDBgrid->GetTotBins(); 

  cerr<<" Fit bins: "<<knTotBins<<endl;

  Int_t bin=0;
  Int_t        *n = new Int_t[knTotBins];
  Double_t     *p = new Double_t[knTotBins];
  Double_t    *ep = new Double_t[knTotBins];
  Double_t  *mean = new Double_t[knTotBins];
  Double_t *sigma = new Double_t[knTotBins];

  for(Int_t l=0; l<knTotBins; l++) {
    n[l] = 1; // set to 1 to avoid divisions by 0
    p[l] = mean[l] = sigma[l] = ep[l] = 0.; 
  }

  // loop on chain entries for the mean
  for(Int_t l=0; l<entries; l++) {
    cmptrktree->GetEvent(l);
    if(TMath::Abs(cmptrk.pdg)!=pdg) continue;
    bin = fDBgrid->GetBin(cmptrk.pt,cmptrk.eta);
    p[bin] += cmptrk.p;
    mean[bin] += cmptrk.dEdx;
    n[bin]++;
  } // loop on chain entries

  for(Int_t l=0; l<knTotBins; l++) {
    p[l] /= n[l];
    mean[l] /= n[l];
    n[l] = 1; // set to 1 to avoid divisions by 0
  }

  // loop on chain entries for the sigma
  for(Int_t l=0; l<entries; l++) {
    cmptrktree->GetEvent(l);
    if(TMath::Abs(cmptrk.pdg)!=pdg) continue;
    bin = fDBgrid->GetBin(cmptrk.pt,cmptrk.eta);
    if(cmptrk.p<1. && TMath::Abs(cmptrk.p-p[bin])>0.025) continue;
    n[bin]++;
    sigma[bin] += (cmptrk.dEdx-mean[bin])*(cmptrk.dEdx-mean[bin]);
  } // loop on chain entries
  
  for(Int_t l=0; l<knTotBins; l++) {
    sigma[l] = TMath::Sqrt(sigma[l]/n[l]);
  }

  // create the canvas
  TCanvas *canv = new TCanvas("canv","dEdx",0,0,900,700); 

  // create the graph for dEdx vs p
  TGraphErrors *gr = new TGraphErrors(knTotBins,p,mean,ep,sigma);
  TString title("  : dE/dx vs momentum"); title.Prepend(part);
  TH2F *frame = new TH2F("frame1",title.Data(),2,0.1,50,2,0,ymax);
  frame->SetXTitle("p [GeV/c]");
  frame->SetYTitle("dE/dx [a.u.]");
  canv->SetLogx();
  frame->Draw();
  gr->Draw("P");

  switch(pdg) {
  case 211:
    for(Int_t i=0; i<knTotBins; i++) {
      fdEdxMeanPi.SetParam(i,mean[i]);
      fdEdxRMSPi.SetParam(i,sigma[i]);
    }    
    break;
  case 321:
    for(Int_t i=0; i<knTotBins; i++) {
      fdEdxMeanKa.SetParam(i,mean[i]);
      fdEdxRMSKa.SetParam(i,sigma[i]);
    }    
    break;
  case 2212:
    for(Int_t i=0; i<knTotBins; i++) {
      fdEdxMeanPr.SetParam(i,mean[i]);
      fdEdxRMSPr.SetParam(i,sigma[i]);
    }    
    break;
  case 11:
    for(Int_t i=0; i<knTotBins; i++) {
      fdEdxMeanEl.SetParam(i,mean[i]);
      fdEdxRMSEl.SetParam(i,sigma[i]);
    }    
    break;
  case 13:
    for(Int_t i=0; i<knTotBins; i++) {
      fdEdxMeanMu.SetParam(i,mean[i]);
      fdEdxRMSMu.SetParam(i,sigma[i]);
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

  /*
  // create a chain with compared tracks
  TChain *cmptrkchain = new ("cmptrktree");
  cmptrkchain.Add("CovMatrix_AllEvts.root");
  //cmptrkchain.Add("CovMatrix_AllEvts_1.root");
  //cmptrkchain.Add("CovMatrix_AllEvts_2.root");
  //cmptrkchain.Add("CovMatrix_AllEvts_3.root");
  */


  TFile *infile = TFile::Open("CovMatrix_AllEvts.root");
  TTree *cmptrktree = (TTree*)infile->Get("cmptrktree");

  COMPTRACK cmptrk; 
  cmptrktree->SetBranchAddress("comptracks",&cmptrk);
  Int_t entries = (Int_t)cmptrktree->GetEntries(); 
  cerr<<" Number of entries: "<<entries<<endl;

  Int_t thisPdg=0;
  Char_t hname[100], htitle[100];
  Int_t nTotBins,bin;

  AliTPCkineGrid  pulls[5];
  TH1F *hDum = new TH1F("name","title",100,-7.,7.);
  TF1 *g = new TF1("g","gaus");

  InitializeKineGrid("pulls");
  InitializeKineGrid("DB");



  // loop on the particles Pi,Ka,Pr,El,Mu
  for(Int_t part=0; part<5; part++) {

    switch (part) {
    case 0: // pions
      thisPdg=211; 
      cerr<<" Processing pions ...\n";
      break;   
    case 1: // kaons
      thisPdg=321; 
      cerr<<" Processing kaons ...\n";
      break;
    case 2: // protons
      thisPdg=2212; 
      cerr<<" Processing protons ...\n";
      break;
    case 3: // electrons
      thisPdg=11; 
      cerr<<" Processing electrons ...\n";
      break;
    case 4: // muons
      thisPdg=13; 
      cerr<<" Processing muons ...\n";
      break;
    }

    SetParticle(thisPdg);

    for(Int_t i=0;i<5;i++) {
      pulls[i].~AliTPCkineGrid(); 
      new(&pulls[i]) AliTPCkineGrid(*(fPulls+i));
    }
    nTotBins = fDBgrid->GetTotBins();
    cerr<<"nTotBins = "<<nTotBins<<endl; 

    // create histograms for the all the bins
    TH1F *hPulls0=NULL;
    TH1F *hPulls1=NULL;
    TH1F *hPulls2=NULL;
    TH1F *hPulls3=NULL;
    TH1F *hPulls4=NULL;

    hPulls0 = new TH1F[nTotBins]; 
    hPulls1 = new TH1F[nTotBins]; 
    hPulls2 = new TH1F[nTotBins]; 
    hPulls3 = new TH1F[nTotBins]; 
    hPulls4 = new TH1F[nTotBins]; 


    for(Int_t i=0; i<nTotBins; i++) {
      sprintf(hname,"hPulls0%d",i);
      sprintf(htitle,"P0 pulls for bin %d",i);
      hDum->SetName(hname); hDum->SetTitle(htitle);
      hPulls0[i] = *hDum;
      sprintf(hname,"hPulls1%d",i);
      sprintf(htitle,"P1 pulls for bin %d",i);
      hDum->SetName(hname); hDum->SetTitle(htitle);
      hPulls1[i] = *hDum;
      sprintf(hname,"hPulls2%d",i);
      sprintf(htitle,"P2 pulls for bin %d",i);
      hDum->SetName(hname); hDum->SetTitle(htitle);
      hPulls2[i] = *hDum;
      sprintf(hname,"hPulls3%d",i);
      sprintf(htitle,"P3 pulls for bin %d",i);
      hDum->SetName(hname); hDum->SetTitle(htitle);
      hPulls3[i] = *hDum;
      sprintf(hname,"hPulls4%d",i);
      sprintf(htitle,"P4 pulls for bin %d",i);
      hDum->SetName(hname); hDum->SetTitle(htitle);
      hPulls4[i] = *hDum;
    }

    // loop on chain entries 
    for(Int_t i=0; i<entries; i++) {
      cmptrktree->GetEvent(i);
      if(TMath::Abs(cmptrk.pdg)!=thisPdg) continue;
      // fill histograms with the pulls
      bin = fDBgrid->GetBin(cmptrk.pt,cmptrk.eta);
      //cerr<<" pt "<<cmptrk.pt<<"   eta "<<cmptrk.eta<<"   bin "<<bin<<endl; 
      hPulls0[bin].Fill(cmptrk.dP0/TMath::Sqrt(cmptrk.c00));
      hPulls1[bin].Fill(cmptrk.dP1/TMath::Sqrt(cmptrk.c11));
      hPulls2[bin].Fill(cmptrk.dP2/TMath::Sqrt(cmptrk.c22));
      hPulls3[bin].Fill(cmptrk.dP3/TMath::Sqrt(cmptrk.c33));
      hPulls4[bin].Fill(cmptrk.dP4/TMath::Sqrt(cmptrk.c44));
    } // loop on chain entries

    // compute the sigma of the distributions
    for(Int_t i=0; i<nTotBins; i++) {
      if(hPulls0[i].GetEntries()>10) {
	g->SetRange(-3.*hPulls0[i].GetRMS(),3.*hPulls0[i].GetRMS());
	hPulls0[i].Fit("g","R,Q,N");
	pulls[0].SetParam(i,g->GetParameter(2));
      } else pulls[0].SetParam(i,-1.);
      if(hPulls1[i].GetEntries()>10) {
	g->SetRange(-3.*hPulls1[i].GetRMS(),3.*hPulls1[i].GetRMS());
	hPulls1[i].Fit("g","R,Q,N");
	pulls[1].SetParam(i,g->GetParameter(2));
      } else pulls[1].SetParam(i,-1.);
      if(hPulls2[i].GetEntries()>10) {
	g->SetRange(-3.*hPulls2[i].GetRMS(),3.*hPulls2[i].GetRMS());
	hPulls2[i].Fit("g","R,Q,N");
	pulls[2].SetParam(i,g->GetParameter(2));
      } else pulls[2].SetParam(i,-1.);
      if(hPulls3[i].GetEntries()>10) {
	g->SetRange(-3.*hPulls3[i].GetRMS(),3.*hPulls3[i].GetRMS());
	hPulls3[i].Fit("g","R,Q,N");
	pulls[3].SetParam(i,g->GetParameter(2));
      } else pulls[3].SetParam(i,-1.);
      if(hPulls4[i].GetEntries()>10) {
	g->SetRange(-3.*hPulls4[i].GetRMS(),3.*hPulls4[i].GetRMS());
	hPulls4[i].Fit("g","R,Q,N");
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
    case 2: // protons
      for(Int_t i=0;i<5;i++) {
	fPullsPr[i].~AliTPCkineGrid();
	new(&fPullsPr[i]) AliTPCkineGrid(pulls[i]);
      }
      break;
    case 3: // electrons
      for(Int_t i=0;i<5;i++) {
	fPullsEl[i].~AliTPCkineGrid();
	new(&fPullsEl[i]) AliTPCkineGrid(pulls[i]);
      }
      break;
    case 4: // muons
      for(Int_t i=0;i<5;i++) {
	fPullsMu[i].~AliTPCkineGrid();
	new(&fPullsMu[i]) AliTPCkineGrid(pulls[i]);
	//cerr<<" mu  pulls "<<i<<"  "<<fPullsMu[i].GetParam(0)<<endl;
      }
      break;
    }

    delete [] hPulls0;
    delete [] hPulls1;
    delete [] hPulls2;
    delete [] hPulls3;
    delete [] hPulls4;
    
  } // loop on particle species

  // write pulls to file
  WritePulls(outName);


  return;
}
//-----------------------------------------------------------------------------
void AliTPCtrackerParam::AnalyzeResolutions(Int_t pdg) {
//-----------------------------------------------------------------------------
// This function computes the resolutions:
// - dy 
// - dC 
// - dPt/Pt
// as a function of Pt
// Input file is CovMatrix_AllEvts.root.  
//-----------------------------------------------------------------------------

  /*
  // create a chain with compared tracks
  TChain *cmptrkchain = new ("cmptrktree");
  cmptrkchain.Add("CovMatrix_AllEvts.root");
  //cmptrkchain.Add("CovMatrix_AllEvts_1.root");
  //cmptrkchain.Add("CovMatrix_AllEvts_2.root");
  //cmptrkchain.Add("CovMatrix_AllEvts_3.root");
  */


  TFile *infile = TFile::Open("CovMatrix_AllEvts.root");
  TTree *cmptrktree = (TTree*)infile->Get("cmptrktree");

  COMPTRACK cmptrk; 
  cmptrktree->SetBranchAddress("comptracks",&cmptrk);
  Int_t entries = (Int_t)cmptrktree->GetEntries(); 
  cerr<<" Number of entries: "<<entries<<endl;


  Int_t bin = 0;

  InitializeKineGrid("DB");
  InitializeKineGrid("eff");

  SetParticle(pdg);

  const Int_t knPtBins = fEff->GetPointsPt();
  cerr<<"knPtBins = "<<knPtBins<<endl; 
  Double_t *dP0     = new Double_t[knPtBins];
  Double_t *dP4     = new Double_t[knPtBins];
  Double_t *dPtToPt = new Double_t[knPtBins];
  Double_t *pt      = new Double_t[knPtBins];
  fEff->GetArrayPt(pt);


  TH1F *hDumP0 = new TH1F("nameP0","dy",100,-0.3,0.3);
  TH1F *hDumP4 = new TH1F("nameP4","dC",100,-0.0005,0.0005);
  TH1F *hDumPt = new TH1F("namePt","dp_{T}/p_{T}",100,-0.5,0.5);

  TF1 *g = new TF1("g","gaus");

  // create histograms for the all the bins
  TH1F *hP0=NULL;
  TH1F *hP4=NULL;
  TH1F *hPt=NULL;

  hP0 = new TH1F[knPtBins]; 
  hP4 = new TH1F[knPtBins]; 
  hPt = new TH1F[knPtBins]; 

  for(Int_t i=0; i<knPtBins; i++) {
    hP0[i] = *hDumP0;
    hP4[i] = *hDumP4;
    hPt[i] = *hDumPt;
  }

  // loop on chain entries 
  for(Int_t i=0; i<entries; i++) {
    cmptrktree->GetEvent(i);
    if(TMath::Abs(cmptrk.pdg)!=pdg) continue;
    // fill histograms with the residuals
    bin = (Int_t)fDBgrid->GetBin(cmptrk.pt,cmptrk.eta)/fDBgrid->GetBinsEta();
    //cerr<<" pt "<<cmptrk.pt<<"   eta "<<cmptrk.eta<<"   bin "<<bin<<endl; 
    hP0[bin].Fill(cmptrk.dP0);
    hP4[bin].Fill(cmptrk.dP4);
    hPt[bin].Fill(cmptrk.dpt/cmptrk.pt);
  } // loop on chain entries


  TCanvas *cP0res = new TCanvas("cP0res","cP0res",0,0,1200,700);
  cP0res->Divide(5,2);
  TCanvas *cP4res = new TCanvas("cP4res","cP4res",0,0,1200,700);
  cP4res->Divide(5,2);
  TCanvas *cPtres = new TCanvas("cPtres","cPtres",0,0,1200,700);
  cPtres->Divide(5,2);

  // Draw histograms
  for(Int_t i=0; i<knPtBins; i++) {
    cP0res->cd(i+1); hP0[i].Draw();
    cP4res->cd(i+1); hP4[i].Draw();
    cPtres->cd(i+1); hPt[i].Draw();
  }


  // compute the sigma of the distributions
  for(Int_t i=0; i<knPtBins; i++) {
    if(hP0[i].GetEntries()>10) {
      g->SetRange(-3.*hP0[i].GetRMS(),3.*hP0[i].GetRMS());
      hP0[i].Fit("g","R,Q,N");
      dP0[i] = g->GetParameter(2);
    } else dP0[i] = 0.;
    if(hP4[i].GetEntries()>10) {
      g->SetRange(-3.*hP4[i].GetRMS(),3.*hP4[i].GetRMS());
      hP4[i].Fit("g","R,Q,N");
      dP4[i] = g->GetParameter(2);
    } else dP4[i] = 0.;
    if(hPt[i].GetEntries()>10) {
      g->SetRange(-3.*hPt[i].GetRMS(),3.*hPt[i].GetRMS());
      hPt[i].Fit("g","R,Q,N");
      dPtToPt[i] = 100.*g->GetParameter(2);
    } else dPtToPt[i] = 0.;
  } // loop on bins

  
  TGraph *grdP0 = new TGraph(knPtBins,pt,dP0);
  TGraph *grdP4 = new TGraph(knPtBins,pt,dP4);
  TGraph *grdPtToPt = new TGraph(knPtBins,pt,dPtToPt);

  grdP0->SetMarkerStyle(20); grdP0->SetMarkerColor(2); grdP0->SetMarkerSize(1.5);
  grdP4->SetMarkerStyle(21); grdP4->SetMarkerColor(3); grdP4->SetMarkerSize(1.5);
  grdPtToPt->SetMarkerStyle(22); grdPtToPt->SetMarkerColor(4); grdPtToPt->SetMarkerSize(1.5);

  // Plot Results
  gStyle->SetOptStat(0);
  TCanvas *c1 = new TCanvas("c1","dP0",0,0,900,700);
  c1->SetLogx();
  c1->SetGridx();
  c1->SetGridy();

  TH2F *frame1 = new TH2F("frame1","y resolution VS p_{T} in TPC",2,0.1,30,2,0,0.1);
  frame1->SetXTitle("p_{T} [GeV/c]");
  frame1->SetYTitle("#sigma(y) [cm]");
  frame1->Draw();
  grdP0->Draw("P");


  TCanvas *c2 = new TCanvas("c2","dP4",0,0,900,700);
  c2->SetLogx();
  c2->SetGridx();
  c2->SetGridy();

  TH2F *frame2 = new TH2F("frame2","C resolution VS p_{T} in TPC",2,0.1,30,2,0,0.0001);
  frame2->SetXTitle("p_{T} [GeV/c]");
  frame2->SetYTitle("#sigma(C) [1/cm]");
  frame2->Draw();
  grdP4->Draw("P");

  TCanvas *c3 = new TCanvas("c3","dPtToPt",0,0,900,700);
  c3->SetLogx();
  c3->SetLogy();
  c3->SetGridx();
  c3->SetGridy();

  TH2F *frame3 = new TH2F("frame3","Relative p_{T} resolution VS p_{T} in TPC",2,0.1,30,2,0.1,30.);
  frame3->SetXTitle("p_{T} [GeV/c]");
  frame3->SetYTitle("dp_{T}/p_{T} (%)");
  frame3->Draw();
  grdPtToPt->Draw("P");


  delete [] dP0;
  delete [] dP4;
  delete [] dPtToPt;
  delete [] pt;

  
  delete [] hP0;
  delete [] hP4;
  delete [] hPt;
  
  return;
}
//-----------------------------------------------------------------------------
void AliTPCtrackerParam::BuildTrack(AliTPCseedGeant *s,Int_t ch) {  
//-----------------------------------------------------------------------------
// This function uses GEANT info to set true track parameters
//-----------------------------------------------------------------------------
  Double_t xref = s->GetXL();
  Double_t xx[5],cc[15];
  cc[0]=cc[2]=cc[5]=cc[9]=cc[14]=10.;
  cc[1]=cc[3]=cc[4]=cc[6]=cc[7]=cc[8]=cc[10]=cc[11]=cc[12]=cc[13]=0.;
  
  // Magnetic field
  TVector3 bfield(0.,0.,-fBz);
  
  
  // radius [cm] of track projection in (x,y) 
  Double_t rho = s->GetPt()*100./0.299792458/TMath::Abs(bfield.Z());
  // center of track projection in local reference frame
  TVector3 sMom,sPos;


  // position (local) and momentum (local) at the seed
  // in the bending plane (z=0)
  sPos.SetXYZ(s->GetXL(),s->GetYL(),0.);
  sMom.SetXYZ(s->GetPx()*TMath::Cos(s->GetAlpha())+s->GetPy()*TMath::Sin(s->GetAlpha()),-s->GetPx()*TMath::Sin(s->GetAlpha())+s->GetPy()*TMath::Cos(s->GetAlpha()),0.);
  TVector3 vrho = sMom.Cross(bfield);
  vrho *= ch;
  vrho.SetMag(rho);

  TVector3 vcenter = sPos+vrho;

  Double_t x0 = vcenter.X();

  // fX     = xref         X-coordinate of this track (reference plane)
  // fAlpha = Alpha        Rotation angle the local (TPC sector) 
  // fP0    = YL           Y-coordinate of a track
  // fP1    = ZG           Z-coordinate of a track
  // fP2    = sin(phi)     sine of the (local) azimuthal angle
  // fP3    = Tgl          tangent of the track momentum dip angle
  // fP4    = C            track curvature
  xx[0] = s->GetYL();
  xx[1] = s->GetZL();
  xx[2] = ch/rho*(xref-x0);
  xx[3] = s->GetPz()/s->GetPt();
  xx[4] = ch/rho;

  // create the object AliTPCtrack    
  AliTPCtrack track(xref,s->GetAlpha(),xx,cc,0);
  new(&fTrack) AliTPCtrack(track);

  return;
}
//-----------------------------------------------------------------------------
void AliTPCtrackerParam::CompareTPCtracks(
			   Int_t nEvents,
			   const Char_t* galiceName,
			   const Char_t* trkGeaName,
			   const Char_t* trkKalName,
			   const Char_t* covmatName,
			   const Char_t* tpceffasciiName,
			   const Char_t* tpceffrootName) {
//-----------------------------------------------------------------------------
// This function compares tracks from TPC Kalman Filter V2 with 
// geant tracks at TPC 1st hit. It gives:
//   - a tree with Kalman cov. matrix elements, resolutions, dEdx
//   - the efficiencies as a function of particle type, pT, eta
//-----------------------------------------------------------------------------

  TFile *kalFile    = TFile::Open(trkKalName);
  TFile *geaFile    = TFile::Open(trkGeaName);
  TFile *galiceFile = TFile::Open(galiceName);

  // get the AliRun object
  AliRun *gAlice = (AliRun*)galiceFile->Get("gAlice");


  // create the tree for comparison results
  COMPTRACK cmptrk;
  TTree *cmptrktree = new TTree("cmptrktree","results of track comparison");
  cmptrktree->Branch("comptracks",&cmptrk,"pdg/I:bin:r/D:p:pt:cosl:eta:dpt:dP0:dP1:dP2:dP3:dP4:c00:c10:c11:c20:c21:c22:c30:c31:c32:c33:c40:c41:c42:c43:c44:dEdx");

  InitializeKineGrid("eff");
  InitializeKineGrid("DB");
  Int_t   effBins = fEffPi.GetTotPoints();
  Int_t effBinsPt = fEffPi.GetPointsPt();
  Double_t    *pt = new Double_t[effBinsPt];
  fEffPi.GetArrayPt(pt);

  TParticle *part;
  Double_t ptgener;
  Bool_t   usethis;
  Int_t    label;
  Double_t dAlpha;
  Int_t    pi=0,ka=0,mu=0,el=0,pr=0;
  Int_t   *geaPi = new Int_t[effBins];
  Int_t   *geaKa = new Int_t[effBins];
  Int_t   *geaPr = new Int_t[effBins];
  Int_t   *geaEl = new Int_t[effBins];
  Int_t   *geaMu = new Int_t[effBins];
  Int_t   *kalPi = new Int_t[effBins];
  Int_t   *kalKa = new Int_t[effBins];
  Int_t   *kalPr = new Int_t[effBins];
  Int_t   *kalEl = new Int_t[effBins];
  Int_t   *kalMu = new Int_t[effBins];
  Float_t *effPi = new Float_t[effBins];
  Float_t *effKa = new Float_t[effBins];
  Float_t *effPr = new Float_t[effBins];
  Float_t *effEl = new Float_t[effBins];
  Float_t *effMu = new Float_t[effBins];
  Int_t bin;
  for(Int_t j=0; j<effBins; j++) {
    geaPi[j]=geaKa[j]=geaPr[j]=geaEl[j]=geaMu[j]=0;
    kalPi[j]=kalKa[j]=kalPr[j]=kalEl[j]=kalMu[j]=0;
    effPi[j]=effKa[j]=effPr[j]=effEl[j]=effMu[j]=-1.;
  }

  Char_t tname[100];

  // loop on events in file
  for(Int_t evt=0; evt<nEvents; evt++) {  
    cerr<<"\n  --- Reading tracks for event "<<evt<<" ---\n\n";
    sprintf(tname,"TreeT_TPC_%d",evt);
    
    // particles from TreeK
    const Int_t knparticles = gAlice->GetEvent(evt);

    Int_t *kalLab = new Int_t[knparticles];
    for(Int_t i=0; i<knparticles; i++) kalLab[i] = -1; 
 

    // tracks from Kalman
    TTree *kaltree=(TTree*)kalFile->Get(tname);
    if(!kaltree) continue;
    AliTPCtrack *kaltrack=new AliTPCtrack; 
    kaltree->SetBranchAddress("tracks",&kaltrack);
    Int_t kalEntries = (Int_t)kaltree->GetEntries();

    // tracks from 1st hit
    TTree *geatree=(TTree*)geaFile->Get(tname);
    if(!geatree) continue;
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
    
    /*
    // Read the labels of the seeds
    char sname[100];
    Int_t sLabel,ncol;
    Bool_t *hasSeed = new Bool_t[knparticles];
    for(Int_t i=0; i<knparticles; i++) hasSeed[i] = kFALSE; 
    sprintf(sname,"seedLabels.%d.dat",evt);
    FILE *seedFile = fopen(sname,"r");
    while(1) {
      ncol = fscanf(seedFile,"%d",&sLabel);
      if(ncol<1) break;
      if(sLabel>0) hasSeed[sLabel]=kTRUE;
    }
    fclose(seedFile);
    */

    cerr<<"Doing track comparison...\n";
    // loop on tracks at 1st hit
    for(Int_t j=0; j<geaEntries; j++) {
      geatree->GetEvent(j);
      
      label = geatrack->GetLabel();
      part = (TParticle*)gAlice->GetMCApp()->Particle(label);
      
      // use only injected tracks with fixed values of pT
      ptgener = part->Pt();
      usethis = kFALSE;
      for(Int_t l=0; l<fEffPi.GetPointsPt(); l++) {
	if(TMath::Abs(ptgener-pt[l])<0.01) usethis = kTRUE;
      }
      if(!usethis) continue;

      // check if it has the seed
      //if(!hasSeed[label]) continue;

      /*
      // check if track is entirely contained in a TPC sector
      Bool_t out = kFALSE;
      for(Int_t l=0; l<10; l++) {
	Double_t x = 85. + (250.-85.)*(Double_t)l/9.;
	//cerr<<"x "<<x<<"    X "<<geatrack->GetX()<<endl;
	Double_t y = geatrack->GetY() + (
        TMath::Sqrt(1-(geatrack->GetC()*geatrack->GetX()-geatrack->GetEta())*
                      (geatrack->GetC()*geatrack->GetX()-geatrack->GetEta()))
       -TMath::Sqrt(1-(geatrack->GetC()*x-geatrack->GetEta())*
		      (geatrack->GetC()*x-geatrack->GetEta()))
        )/geatrack->GetC();
	y = TMath::Abs(y);
	//cerr<<"y "<<y<<"    Y "<<geatrack->GetY()<<endl;
	if(y > 0.8*x*TMath::Tan(TMath::Pi()/18)) { out = kTRUE; break; }
      }
      if(out) continue;      
      */

      cmptrk.pdg = part->GetPdgCode();
      cmptrk.eta = part->Eta();
      cmptrk.r = TMath::Sqrt(part->Vx()*part->Vx()+part->Vy()*part->Vy());
      
      cmptrk.pt   = 1/TMath::Abs(geatrack->Get1Pt());
      cmptrk.cosl = TMath::Cos(TMath::ATan(geatrack->GetTgl()));
      cmptrk.p    = cmptrk.pt/cmptrk.cosl;
    

      bin = fDBgridPi.GetBin(cmptrk.pt,cmptrk.eta);
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
      
      kaltrack->PropagateTo(geatrack->GetX(),1.e-9,0.);
      
      cmptrk.dEdx = kaltrack->GetdEdx();
      
      // compute errors on parameters
      dAlpha = kaltrack->GetAlpha()-geatrack->GetAlpha();
      if(TMath::Abs(dAlpha)>0.1) { cerr<<" ! WRONG SECTOR !\n"; continue; }
      
      cmptrk.dP0 = kaltrack->GetY()-geatrack->GetY();
      cmptrk.dP1 = kaltrack->GetZ()-geatrack->GetZ();
      cmptrk.dP2 = kaltrack->GetSnp()-geatrack->GetSnp();
      cmptrk.dP3 = kaltrack->GetTgl()-geatrack->GetTgl();
      cmptrk.dP4 = kaltrack->GetC()-geatrack->GetC();
      cmptrk.dpt = 1/kaltrack->Get1Pt()-1/geatrack->Get1Pt();
    
      // get covariance matrix
      // beware: lines 3 and 4 in the matrix are inverted!
      //kaltrack->GetCovariance(cc);

      cmptrk.c00 = kaltrack->GetSigmaY2();
      cmptrk.c10 = kaltrack->GetSigmaZY();
      cmptrk.c11 = kaltrack->GetSigmaZ2();
      cmptrk.c20 = kaltrack->GetSigmaSnpY();
      cmptrk.c21 = kaltrack->GetSigmaSnpY();
      cmptrk.c22 = kaltrack->GetSigmaSnp2();
      cmptrk.c30 = kaltrack->GetSigmaTglY();
      cmptrk.c31 = kaltrack->GetSigmaTglZ();
      cmptrk.c32 = kaltrack->GetSigmaTglSnp();
      cmptrk.c33 = kaltrack->GetSigmaTgl2();
      cmptrk.c40 = kaltrack->GetSigma1PtY();
      cmptrk.c41 = kaltrack->GetSigma1PtZ();
      cmptrk.c42 = kaltrack->GetSigma1PtSnp();
      cmptrk.c43 = kaltrack->GetSigma1PtTgl();
      cmptrk.c44 = kaltrack->GetSigma1Pt2();
    
      // fill tree
      cmptrktree->Fill();

    } // loop on tracks at TPC 1st hit

    delete [] kalLab;
    //delete [] hasSeed;

  } // end loop on events in file


  cerr<<"+++\n+++ Number of compared tracks: "<<pi+ka+el+mu+pr<<endl;
  cerr<<"+++ Pions: "<<pi<<", Kaons: "<<ka<<", Protons : "<<pr<<", Electrons: "<<el<<", Muons: "<<mu<<endl;
  cerr<<"+++"<<endl;

  // Write tree to file
  TFile *outfile = new TFile(covmatName,"recreate");
  cmptrktree->Write();
  outfile->Close();


  // Write efficiencies to ascii file
  FILE *effFile = fopen(tpceffasciiName,"w");
  //fprintf(effFile,"%d\n",kalEntries);
  for(Int_t j=0; j<effBins; j++) {
    if(geaPi[j]>=100) effPi[j]=(Float_t)kalPi[j]/geaPi[j];
    if(geaKa[j]>=100) effKa[j]=(Float_t)kalKa[j]/geaKa[j];
    if(geaPr[j]>=100) effPr[j]=(Float_t)kalPr[j]/geaPr[j];
    if(geaEl[j]>=500) effEl[j]=(Float_t)kalEl[j]/geaEl[j];
    if(geaMu[j]>=100) effMu[j]=(Float_t)kalMu[j]/geaMu[j];
    fprintf(effFile,"%f  %f  %f  %f  %f\n",effPi[j],effKa[j],effPr[j],effEl[j],effMu[j]);
  }

  for(Int_t j=0; j<effBins; j++) {
    fprintf(effFile,"%d  %d  %d  %d  %d\n",geaPi[j],geaKa[j],geaPr[j],geaEl[j],geaMu[j]);
  }
  for(Int_t j=0; j<effBins; j++) {
    fprintf(effFile,"%d  %d  %d  %d  %d\n",kalPi[j],kalKa[j],kalPr[j],kalEl[j],kalMu[j]);
  }
  fclose(effFile);

  // Write efficiencies to root file
  for(Int_t j=0; j<effBins; j++) {
    fEffPi.SetParam(j,(Double_t)effPi[j]);
    fEffKa.SetParam(j,(Double_t)effKa[j]);
    fEffPr.SetParam(j,(Double_t)effPr[j]);
    fEffEl.SetParam(j,(Double_t)effEl[j]);
    fEffMu.SetParam(j,(Double_t)effMu[j]);
  }
  WriteEffs(tpceffrootName);

  // delete AliRun object
  delete gAlice; gAlice=0;
  
  // close all input files
  kalFile->Close();
  geaFile->Close();
  galiceFile->Close();

  delete [] pt;
  delete [] effPi;
  delete [] effKa;
  delete [] effPr;
  delete [] effEl;
  delete [] effMu;
  delete [] geaPi;
  delete [] geaKa;
  delete [] geaPr;
  delete [] geaEl;
  delete [] geaMu;
  delete [] kalPi;
  delete [] kalKa;
  delete [] kalPr;
  delete [] kalEl;
  delete [] kalMu;

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

  Double_t massPi = (Double_t)TDatabasePDG::Instance()->GetParticle(211)->Mass();
  Double_t massKa = (Double_t)TDatabasePDG::Instance()->GetParticle(321)->Mass();
  Double_t massPr = (Double_t)TDatabasePDG::Instance()->GetParticle(2212)->Mass();

  if (p<0.6) {
    if (dEdx < 39.+ 12./(p+0.25)/(p+0.25)) { 
      t.AssignMass(massPi); new(&fTrack) AliTPCtrack(t); return;
    }
    if (dEdx < 39.+ 12./p/p) { 
      t.AssignMass(massKa); new(&fTrack) AliTPCtrack(t); return;
    }
    t.AssignMass(massPr); new(&fTrack) AliTPCtrack(t); return;
  }

  if (p<1.2) {
    if (dEdx < 39.+ 12./(p+0.25)/(p+0.25)) { 
      t.AssignMass(massPi); new(&fTrack) AliTPCtrack(t); return;
    }
    t.AssignMass(massPr); new(&fTrack) AliTPCtrack(t); return;
  }

  t.AssignMass(massPi); new(&fTrack) AliTPCtrack(t); return;
}
//-----------------------------------------------------------------------------
void AliTPCtrackerParam::CookTrack(Double_t pt,Double_t eta) {
//-----------------------------------------------------------------------------
// This function deals with covariance matrix and smearing
//-----------------------------------------------------------------------------

  COVMATRIX covmat;
  Double_t  p,cosl;
  Double_t  trkKine[1],trkRegPar[3];    
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

  // get covariance matrix from regularized matrix    
  for(Int_t l=0;l<3;l++) trkRegPar[l] = (*fRegPar)(0,l);
  cc[0] = covmat.c00*RegFunc(trkKine,trkRegPar);
  cc[1] = covmat.c10;
  for(Int_t l=0;l<3;l++) trkRegPar[l] = (*fRegPar)(1,l);
  cc[2] = covmat.c11*RegFunc(trkKine,trkRegPar);
  for(Int_t l=0;l<3;l++) trkRegPar[l] = (*fRegPar)(2,l);
  cc[3] = covmat.c20*RegFunc(trkKine,trkRegPar);
  cc[4] = covmat.c21;
  for(Int_t l=0;l<3;l++) trkRegPar[l] = (*fRegPar)(3,l);
  cc[5] = covmat.c22*RegFunc(trkKine,trkRegPar);
  cc[6] = covmat.c30;
  for(Int_t l=0;l<3;l++) trkRegPar[l] = (*fRegPar)(4,l);
  cc[7] = covmat.c31*RegFunc(trkKine,trkRegPar);
  cc[8] = covmat.c32;
  for(Int_t l=0;l<3;l++) trkRegPar[l] = (*fRegPar)(5,l);
  cc[9] = covmat.c33*RegFunc(trkKine,trkRegPar);
  for(Int_t l=0;l<3;l++) trkRegPar[l] = (*fRegPar)(6,l);
  cc[10]= covmat.c40*RegFunc(trkKine,trkRegPar);
  cc[11]= covmat.c41;
  for(Int_t l=0;l<3;l++) trkRegPar[l] = (*fRegPar)(7,l);
  cc[12]= covmat.c42*RegFunc(trkKine,trkRegPar);
  cc[13]= covmat.c43;
  for(Int_t l=0;l<3;l++) trkRegPar[l] = (*fRegPar)(8,l);
  cc[14]= covmat.c44*RegFunc(trkKine,trkRegPar);  


  TMatrixD covMatSmear(5,5);    
  covMatSmear = GetSmearingMatrix(cc,pt,eta);

  // get track original parameters
  xref =fTrack.GetX();
  alpha=fTrack.GetAlpha();
  xx[0]=fTrack.GetY();
  xx[1]=fTrack.GetZ();
  xx[2]=fTrack.GetSnp();
  xx[3]=fTrack.GetTgl();
  xx[4]=fTrack.GetC();
    
  // use smearing matrix to smear the original parameters
  xxsm[0]=xref;
  SmearTrack(xx,xxsm,covMatSmear);
    
  AliTPCtrack track(xref,alpha,xxsm,cc,0);
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

  const Int_t kn = fEff->GetPointsPt();
  Double_t *effsA = new Double_t[kn];
  Double_t *effsB = new Double_t[kn];
  Double_t *effsC = new Double_t[kn];
  Double_t *pt    = new Double_t[kn];

  fEff->GetArrayPt(pt);
  for(Int_t i=0;i<kn;i++) {
    effsA[i] = fEff->GetParam(i,0);
    effsB[i] = fEff->GetParam(i,1);
    effsC[i] = fEff->GetParam(i,2);
  }
  
  TGraph *grA = new TGraph(kn,pt,effsA);
  TGraph *grB = new TGraph(kn,pt,effsB);
  TGraph *grC = new TGraph(kn,pt,effsC);

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

  const Int_t kn = (fPulls+par)->GetPointsPt();
  Double_t *pullsA = new Double_t[kn];
  Double_t *pullsB = new Double_t[kn];
  Double_t *pullsC = new Double_t[kn];
  Double_t *pt     = new Double_t[kn];
  (fPulls+par)->GetArrayPt(pt);  
  for(Int_t i=0;i<kn;i++) {
    pullsA[i] = (fPulls+par)->GetParam(i,0);
    pullsB[i] = (fPulls+par)->GetParam(i,1);
    pullsC[i] = (fPulls+par)->GetParam(i,2);
  }

  TGraph *grA = new TGraph(kn,pt,pullsA);
  TGraph *grB = new TGraph(kn,pt,pullsB);
  TGraph *grC = new TGraph(kn,pt,pullsC);

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

  for(Int_t i=0;i<5;i++) {
    stretchMat(i,i) = (fPulls+i)->GetValueAt(pt,eta); 
    if(stretchMat(i,i)==0.) stretchMat(i,i) = 1.;
  }

  TMatrixD mat(stretchMat,TMatrixD::kMult,covMat);
  TMatrixD covMatSmear(mat,TMatrixD::kMult,stretchMat);

  return covMatSmear;
}
//-----------------------------------------------------------------------------
void AliTPCtrackerParam::InitializeKineGrid(Option_t* which) {
//-----------------------------------------------------------------------------
// This function initializes ([pt,eta] points) the data members AliTPCkineGrid
// which = "DB"     -> initialize fDBgrid... members
//         "eff"    -> initialize fEff... members
//         "pulls"  -> initialize fPulls... members
//         "dEdx"   -> initialize fdEdx... members
//-----------------------------------------------------------------------------

  const char *db     = strstr(which,"DB");
  const char *eff    = strstr(which,"eff");
  const char *pulls  = strstr(which,"pulls");
  const char *dEdx   = strstr(which,"dEdx");


  Int_t nEta=0, nPt=0;

  Double_t etaPoints[2] = {0.3,0.6};
  Double_t etaBins[3]   = {0.15,0.45,0.75};

  Double_t ptPoints[9]  = {0.4,0.6,0.8,1.2,1.7,3.,5.,8.,15.};
  Double_t ptBins[10]  = {0.3,0.5,0.7,1.,1.5,2.,4.,6.,10.,20.};


  Double_t *eta=0,*pt=0;

  if(db) {
    nEta = 2;
    nPt  = 9;
    eta  = etaPoints;
    pt   = ptPoints;
  } else {
    nEta = 3;
    nPt  = 10;
    eta  = etaBins;
    pt   = ptBins;
  }

  AliTPCkineGrid *dummy=0;

  if(db) {    
    dummy = new AliTPCkineGrid(nPt,nEta,pt,eta);
    new(&fDBgridPi) AliTPCkineGrid(*dummy);
    new(&fDBgridKa) AliTPCkineGrid(*dummy);
    new(&fDBgridPr) AliTPCkineGrid(*dummy);
    new(&fDBgridEl) AliTPCkineGrid(*dummy);
    new(&fDBgridMu) AliTPCkineGrid(*dummy);
    delete dummy;
  }
  if(eff) {    
    dummy = new AliTPCkineGrid(nPt,nEta,pt,eta);
    new(&fEffPi) AliTPCkineGrid(*dummy);
    new(&fEffKa) AliTPCkineGrid(*dummy);
    new(&fEffPr) AliTPCkineGrid(*dummy);
    new(&fEffEl) AliTPCkineGrid(*dummy);
    new(&fEffMu) AliTPCkineGrid(*dummy);
    delete dummy;
  }
  if(pulls) {    
    dummy = new AliTPCkineGrid(nPt,nEta,pt,eta);
    for(Int_t i=0;i<5;i++) new(&fPullsPi[i]) AliTPCkineGrid(*dummy);
    for(Int_t i=0;i<5;i++) new(&fPullsKa[i]) AliTPCkineGrid(*dummy);
    for(Int_t i=0;i<5;i++) new(&fPullsPr[i]) AliTPCkineGrid(*dummy);
    for(Int_t i=0;i<5;i++) new(&fPullsEl[i]) AliTPCkineGrid(*dummy);
    for(Int_t i=0;i<5;i++) new(&fPullsMu[i]) AliTPCkineGrid(*dummy);
    delete dummy;
  }
  if(dEdx) {    
    dummy = new AliTPCkineGrid(nPt,nEta,pt,eta);
    new(&fdEdxMeanPi) AliTPCkineGrid(*dummy);
    new(&fdEdxRMSPi) AliTPCkineGrid(*dummy);
    new(&fdEdxMeanKa) AliTPCkineGrid(*dummy);
    new(&fdEdxRMSKa) AliTPCkineGrid(*dummy);
    new(&fdEdxMeanPr) AliTPCkineGrid(*dummy);
    new(&fdEdxRMSPr) AliTPCkineGrid(*dummy);
    new(&fdEdxMeanEl) AliTPCkineGrid(*dummy);
    new(&fdEdxRMSEl) AliTPCkineGrid(*dummy);
    new(&fdEdxMeanMu) AliTPCkineGrid(*dummy);
    new(&fdEdxRMSMu) AliTPCkineGrid(*dummy);
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
  const char *effFile   ="TPCeff.root";
  const char *pullsFile ="pulls.root";
  const char *regPiFile ="regPi.root";
  const char *regKaFile ="regKa.root";
  const char *regPrFile ="regPr.root";
  const char *regElFile ="regEl.root";
  const char *regMuFile ="regMu.root";
  const char *dEdxPiFile="dEdxPi.root";
  const char *dEdxKaFile="dEdxKa.root";
  const char *dEdxPrFile="dEdxPr.root";
  const char *dEdxElFile="dEdxEl.root";
  const char *dEdxMuFile="dEdxMu.root";
  const char *cmFile    ="CovMatrix_AllEvts.root";
  /*
  Char_t *cmFile1  ="CovMatrix_AllEvts_1.root";
  Char_t *cmFile2  ="CovMatrix_AllEvts_2.root";
  Char_t *cmFile3  ="CovMatrix_AllEvts_3.root";
  */

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
  ReadRegParams(regPrFile,2212);
  WriteRegParams(fDBfileName.Data(),2212);
  ReadRegParams(regElFile,11);
  WriteRegParams(fDBfileName.Data(),11);
  ReadRegParams(regMuFile,13);
  WriteRegParams(fDBfileName.Data(),13);

  // store the dEdx parameters
  ReaddEdx(dEdxPiFile,211);
  WritedEdx(fDBfileName.Data(),211);
  ReaddEdx(dEdxKaFile,321);
  WritedEdx(fDBfileName.Data(),321);
  ReaddEdx(dEdxPrFile,2212);
  WritedEdx(fDBfileName.Data(),2212);
  ReaddEdx(dEdxElFile,11);
  WritedEdx(fDBfileName.Data(),11);
  ReaddEdx(dEdxMuFile,13);
  WritedEdx(fDBfileName.Data(),13);


  //
  // store the regularized covariance matrices
  //
  InitializeKineGrid("DB");

  const Int_t knBinsPi = fDBgridPi.GetTotBins();
  const Int_t knBinsKa = fDBgridKa.GetTotBins();
  const Int_t knBinsPr = fDBgridPr.GetTotBins();
  const Int_t knBinsEl = fDBgridEl.GetTotBins();
  const Int_t knBinsMu = fDBgridMu.GetTotBins();


  // create the trees for cov. matrices
  // trees for pions
  TTree *covTreePi1 = NULL;
  covTreePi1 = new TTree[knBinsPi]; 
  // trees for kaons
  TTree *covTreeKa1 = NULL;
  covTreeKa1 = new TTree[knBinsKa]; 
  // trees for protons
  TTree *covTreePr1 = NULL;
  covTreePr1 = new TTree[knBinsPr]; 
  // trees for electrons
  TTree *covTreeEl1 = NULL;
  covTreeEl1 = new TTree[knBinsEl]; 
  // trees for muons
  TTree *covTreeMu1 = NULL;
  covTreeMu1 = new TTree[knBinsMu]; 

  Char_t hname[100], htitle[100];
  COVMATRIX covmat;

  
  for(Int_t i=0; i<knBinsPi; i++) {
    sprintf(hname,"CovTreePi_bin%d",i);
    sprintf(htitle,"Tree with cov matrix elements for bin %d",i);
    covTreePi1[i].SetName(hname); covTreePi1[i].SetTitle(htitle);
    covTreePi1[i].Branch("matrix",&covmat,"c00/D:c10:c11:c20:c21:c22:c30:c31:c32:c33:c40:c41:c42:c43:c44",5000000);
  }
  for(Int_t i=0; i<knBinsKa; i++) {
    sprintf(hname,"CovTreeKa_bin%d",i);
    sprintf(htitle,"Tree with cov matrix elements for bin %d",i);
    covTreeKa1[i].SetName(hname); covTreeKa1[i].SetTitle(htitle);
    covTreeKa1[i].Branch("matrix",&covmat,"c00/D:c10:c11:c20:c21:c22:c30:c31:c32:c33:c40:c41:c42:c43:c44",1000000);
  }
  for(Int_t i=0; i<knBinsPr; i++) {
    sprintf(hname,"CovTreePr_bin%d",i);
    sprintf(htitle,"Tree with cov matrix elements for bin %d",i);
    covTreePr1[i].SetName(hname); covTreePr1[i].SetTitle(htitle);
    covTreePr1[i].Branch("matrix",&covmat,"c00/D:c10:c11:c20:c21:c22:c30:c31:c32:c33:c40:c41:c42:c43:c44",1000000);
  }
  for(Int_t i=0; i<knBinsEl; i++) {
    sprintf(hname,"CovTreeEl_bin%d",i);
    sprintf(htitle,"Tree with cov matrix elements for bin %d",i);
    covTreeEl1[i].SetName(hname); covTreeEl1[i].SetTitle(htitle);
    covTreeEl1[i].Branch("matrix",&covmat,"c00/D:c10:c11:c20:c21:c22:c30:c31:c32:c33:c40:c41:c42:c43:c44",1000000);
  }
  for(Int_t i=0; i<knBinsMu; i++) {
    sprintf(hname,"CovTreeMu_bin%d",i);
    sprintf(htitle,"Tree with cov matrix elements for bin %d",i);
    covTreeMu1[i].SetName(hname); covTreeMu1[i].SetTitle(htitle);
    covTreeMu1[i].Branch("matrix",&covmat,"c00/D:c10:c11:c20:c21:c22:c30:c31:c32:c33:c40:c41:c42:c43:c44",1000000);
  }

  /*  
  // create the chain with the compared tracks
  TChain *cmptrktree = new TChain("cmptrktree");
  cmptrkchain.Add(cmFile1);
  cmptrkchain.Add(cmFile2);
  cmptrkchain.Add(cmFile3);
  */

  TFile *infile = TFile::Open(cmFile);
  TTree *cmptrktree = (TTree*)infile->Get("cmptrktree");
 
  COMPTRACK cmptrk; 
  cmptrktree->SetBranchAddress("comptracks",&cmptrk);
  Int_t entries = (Int_t)cmptrktree->GetEntries(); 
  cerr<<" Number of entries: "<<entries<<endl;

  Int_t trkPdg,trkBin;
  Double_t trkKine[1],trkRegPar[3]; 
  Int_t *nPerBinPi = new Int_t[knBinsPi];
  for(Int_t k=0;k<knBinsPi;k++) nPerBinPi[k]=0;
  Int_t *nPerBinKa = new Int_t[knBinsKa];
  for(Int_t k=0;k<knBinsKa;k++) nPerBinKa[k]=0;
  Int_t *nPerBinMu = new Int_t[knBinsMu];
  for(Int_t k=0;k<knBinsMu;k++) nPerBinMu[k]=0;
  Int_t *nPerBinEl = new Int_t[knBinsEl];
  for(Int_t k=0;k<knBinsEl;k++) nPerBinEl[k]=0;
  Int_t *nPerBinPr = new Int_t[knBinsPr];
  for(Int_t k=0;k<knBinsPr;k++) nPerBinPr[k]=0;

  // loop on chain entries 
  for(Int_t l=0; l<entries; l++) {
    if(l % 10000 == 0) cerr<<"--- Processing track "<<l<<" of "<<entries<<" ---"<<endl;
    // get the event
    cmptrktree->GetEvent(l);
    // get the pdg code
    trkPdg = TMath::Abs(cmptrk.pdg);
    // use only pions, kaons, protons, electrons, muons
    if(trkPdg!=211 && trkPdg!=321 && trkPdg!=2212 && trkPdg!=11 && trkPdg!=13) continue;
    SetParticle(trkPdg);
    trkBin = fDBgrid->GetBin(cmptrk.pt,cmptrk.eta);
    //cerr<<cmptrk.pt<<"  "<<cmptrk.eta<<"  "<<trkBin<<endl;

    if(trkPdg==211 && nPerBinPi[trkBin]>=5000) continue;
    if(trkPdg==321 && nPerBinKa[trkBin]>=5000) continue;
    if(trkPdg==2212 && nPerBinPr[trkBin]>=5000) continue;
    if(trkPdg==11  && nPerBinEl[trkBin]>=5000) continue;
    if(trkPdg==13 && nPerBinMu[trkBin]>=5000) continue;

    trkKine[0] = cmptrk.p;
    
    // get regularized covariance matrix
    for(Int_t k=0;k<3;k++) trkRegPar[k] = (*fRegPar)(0,k);
    covmat.c00 = cmptrk.c00/RegFunc(trkKine,trkRegPar);
    covmat.c10 = cmptrk.c10;
    for(Int_t k=0;k<3;k++) trkRegPar[k] = (*fRegPar)(1,k);
    covmat.c11 = cmptrk.c11/RegFunc(trkKine,trkRegPar);
    for(Int_t k=0;k<3;k++) trkRegPar[k] = (*fRegPar)(2,k);
    covmat.c20 = cmptrk.c20/RegFunc(trkKine,trkRegPar);
    covmat.c21 = cmptrk.c21;
    for(Int_t k=0;k<3;k++) trkRegPar[k] = (*fRegPar)(3,k);
    covmat.c22 = cmptrk.c22/RegFunc(trkKine,trkRegPar);
    covmat.c30 = cmptrk.c30;
    for(Int_t k=0;k<3;k++) trkRegPar[k] = (*fRegPar)(4,k);
    covmat.c31 = cmptrk.c31/RegFunc(trkKine,trkRegPar);
    covmat.c32 = cmptrk.c32;
    for(Int_t k=0;k<3;k++) trkRegPar[k] = (*fRegPar)(5,k);
    covmat.c33 = cmptrk.c33/RegFunc(trkKine,trkRegPar);
    for(Int_t k=0;k<3;k++) trkRegPar[k] = (*fRegPar)(6,k);
    covmat.c40 = cmptrk.c40/RegFunc(trkKine,trkRegPar);
    covmat.c41 = cmptrk.c41;
    for(Int_t k=0;k<3;k++) trkRegPar[k] = (*fRegPar)(7,k);
    covmat.c42 = cmptrk.c42/RegFunc(trkKine,trkRegPar);
    covmat.c43 = cmptrk.c43;
    for(Int_t k=0;k<3;k++) trkRegPar[k] = (*fRegPar)(8,k);
    covmat.c44 = cmptrk.c44/RegFunc(trkKine,trkRegPar);

    // fill the tree
    switch (trkPdg) {
    case 211: // pions
      covTreePi1[trkBin].Fill();
      nPerBinPi[trkBin]++;
      break;
    case 321: // kaons
      covTreeKa1[trkBin].Fill();
      nPerBinKa[trkBin]++;
      break;
    case 2212: // protons
      covTreePr1[trkBin].Fill();
      nPerBinPr[trkBin]++;
      break;
    case 11: // electrons
      covTreeEl1[trkBin].Fill();
      nPerBinEl[trkBin]++;
      break;
    case 13: // muons
      covTreeMu1[trkBin].Fill();
      nPerBinMu[trkBin]++;
      break;
    }
  } // loop on chain entries

  // store all trees the DB file
  TFile *dbfile = new TFile(fDBfileName.Data(),"update");
  dbfile->mkdir("CovMatrices");
  gDirectory->cd("/CovMatrices");
  gDirectory->mkdir("Pions");
  gDirectory->mkdir("Kaons");
  gDirectory->mkdir("Protons");
  gDirectory->mkdir("Electrons");
  gDirectory->mkdir("Muons");
  // store pions
  gDirectory->cd("/CovMatrices/Pions");
  fDBgridPi.SetName("DBgridPi"); fDBgridPi.Write();
  for(Int_t i=0;i<knBinsPi;i++) covTreePi1[i].Write();
  // store kaons
  gDirectory->cd("/CovMatrices/Kaons");
  fDBgridKa.SetName("DBgridKa"); fDBgridKa.Write();
  for(Int_t i=0;i<knBinsKa;i++) covTreeKa1[i].Write();
  // store kaons
  gDirectory->cd("/CovMatrices/Protons");
  fDBgridPr.SetName("DBgridPr"); fDBgridPr.Write();
  for(Int_t i=0;i<knBinsPr;i++) covTreePr1[i].Write();
  // store electrons
  gDirectory->cd("/CovMatrices/Electrons");
  fDBgridEl.SetName("DBgridEl"); fDBgridEl.Write();
  for(Int_t i=0;i<knBinsEl;i++) covTreeEl1[i].Write();
  // store kaons
  gDirectory->cd("/CovMatrices/Muons");
  fDBgridMu.SetName("DBgridMu"); fDBgridMu.Write();
  for(Int_t i=0;i<knBinsMu;i++) covTreeMu1[i].Write();

  dbfile->Close();
  delete [] nPerBinPi;
  delete [] nPerBinKa;
  delete [] nPerBinPr;
  delete [] nPerBinEl;
  delete [] nPerBinMu; 
 
  return;
}
//-----------------------------------------------------------------------------
void AliTPCtrackerParam::MakeSeedsFromHits(AliTPC *tpc,TTree *th,
					   TObjArray &seedArray) const {
//-----------------------------------------------------------------------------
// This function makes the seeds for tracks from the 1st hits in the TPC
//-----------------------------------------------------------------------------

  Double_t xg,yg,zg,px,py,pz,pt;
  Int_t label;
  Int_t nTracks=(Int_t)th->GetEntries();

  cout<<"+++ Number of \"primary tracks\"(entries in TreeH): "<<nTracks<<
         "\n";
   
  AliTPChit *tpcHit=0;

  // loop over entries in TreeH
  for(Int_t l=0; l<nTracks; l++) {
    if(l%1000==0) cerr<<"  --- Processing primary track "
		      <<l<<" of "<<nTracks<<" ---\r";
    tpc->ResetHits();
    th->GetEvent(l);
    // Get FirstHit
    tpcHit=(AliTPChit*)tpc->FirstHit(-1);
    for( ; tpcHit; tpcHit=(AliTPChit*)tpc->NextHit() ) {
      if(tpcHit->fQ !=0.) continue;
      // Get particle momentum at hit
      px=tpcHit->X(); py=tpcHit->Y(); pz=tpcHit->Z();

      pt=TMath::Sqrt(px*px+py*py);
      // reject hits with Pt<mag*0.45 GeV/c
      if(pt<(fBz*0.45)) continue;

      // Get track label
      label=tpcHit->Track();
      
      if((tpcHit=(AliTPChit*)tpc->NextHit())==0) break;
      if(tpcHit->fQ != 0.) continue;
      // Get global coordinates of hit
      xg=tpcHit->X(); yg=tpcHit->Y(); zg=tpcHit->Z();
      if(TMath::Sqrt(xg*xg+yg*yg)>90.) continue;

      // create the seed
      AliTPCseedGeant *ioseed = new AliTPCseedGeant(xg,yg,zg,px,py,pz,label);

      // reject tracks which are not in the TPC acceptance
      if(!ioseed->InTPCAcceptance()) { delete ioseed; continue; }

      // put seed in array
      seedArray.AddLast(ioseed);
    }
    
  } // loop over entries in TreeH

  return;
}
//-----------------------------------------------------------------------------
void AliTPCtrackerParam::MakeSeedsFromRefs(TTree *ttr,
					   TObjArray &seedArray) const {
//-----------------------------------------------------------------------------
// This function makes the seeds for tracks from the track references
//-----------------------------------------------------------------------------

  Double_t xg,yg,zg,px,py,pz,pt;
  Int_t label;
  Int_t nnn,k;

  TClonesArray *tkRefArray = new TClonesArray("AliTrackReference");

  TBranch *b =(TBranch*)ttr->GetBranch("TPC");
  if(!b) {cerr<<"TPC branch of TreeTR not found"<<endl; return; }
  b->SetAddress(&tkRefArray);
  Int_t nTkRef = (Int_t)b->GetEntries();
  cerr<<"+++ Number of entries in TreeTR(TPC): "<<nTkRef<<"\n";

  // loop on track references
  for(Int_t l=0; l<nTkRef; l++){
    //if(l%1000==0) cerr<<"  --- Processing primary track "
    //	      <<l<<" of "<<nTkRef<<" ---\r";
    if(!b->GetEvent(l)) continue;
    nnn = tkRefArray->GetEntriesFast();

    if(nnn <= 0) continue;   
    for(k=0; k<nnn; k++) {
      AliTrackReference *tkref = (AliTrackReference*)tkRefArray->UncheckedAt(k);

      xg = tkref->X();
      yg = tkref->Y();
      zg = tkref->Z();
      px = tkref->Px();
      py = tkref->Py();
      pz = tkref->Pz();
      label = tkref->GetTrack();

      pt=TMath::Sqrt(px*px+py*py);
      // reject hits with Pt<mag*0.45 GeV/c
      if(pt<(fBz*0.45)) continue;

      // create the seed
      AliTPCseedGeant *ioseed = new AliTPCseedGeant(xg,yg,zg,px,py,pz,label);

      // reject if not at the inner part of TPC
      if(TMath::Abs(ioseed->GetXL()-83.8) > 0.2) { 
	//cerr<<"not at TPC inner part: XL = "<<ioseed->GetXL()<<endl;
	delete ioseed; continue; 
      }

      // reject tracks which are not in the TPC acceptance
      if(!ioseed->InTPCAcceptance()) { 
	delete ioseed; continue; 
      }

      // put seed in array
      seedArray.AddLast(ioseed);

    }
    

  } // loop on track references

  delete tkRefArray;


  return;
}
//-----------------------------------------------------------------------------
void AliTPCtrackerParam::MergeEvents(Int_t evFirst,Int_t evLast) {
//-----------------------------------------------------------------------------
// This function: 1) merges the files from track comparison
//                   (beware: better no more than 100 events per file)
//                2) computes the average TPC efficiencies
//-----------------------------------------------------------------------------

  const char *outName="TPCeff.root";

  // merge files with tracks
  cerr<<" ******** MERGING FILES **********\n\n";

  // create the chain for the tree of compared tracks
  TChain ch1("cmptrktree");
  TChain ch2("cmptrktree");
  TChain ch3("cmptrktree");

  for(Int_t j=evFirst; j<=evLast; j++) {
    cerr<<"Processing event "<<j<<endl;

    TString covName("CovMatrix.");
    covName+=j;
    covName.Append(".root");

    if(gSystem->AccessPathName(covName.Data(),kFileExists)) continue;
       

    if(j<=100) ch1.Add(covName.Data());
    if(j>100 && j<=200) ch2.Add(covName.Data());
    if(j>200) ch3.Add(covName.Data());

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
  cerr<<" ***** EFFICIENCIES ******\n\n";

  ReadEffs("TPCeff.1.root");

  Int_t n = fEffPi.GetTotPoints();
  Double_t *avEffPi = new Double_t[n];
  Double_t *avEffKa = new Double_t[n];
  Double_t *avEffPr = new Double_t[n];
  Double_t *avEffEl = new Double_t[n];
  Double_t *avEffMu = new Double_t[n];
  Int_t    *evtsPi  = new Int_t[n];
  Int_t    *evtsKa  = new Int_t[n];
  Int_t    *evtsPr  = new Int_t[n];
  Int_t    *evtsEl  = new Int_t[n];
  Int_t    *evtsMu  = new Int_t[n];

  for(Int_t j=0; j<n; j++) {
    avEffPi[j]=avEffKa[j]=avEffPr[j]=avEffEl[j]=avEffMu[j]=0.;
    evtsPi[j]=evtsKa[j]=evtsPr[j]=evtsEl[j]=evtsMu[j]=0;
  }

  for(Int_t j=evFirst; j<=evLast; j++) {
    cerr<<"Processing event "<<j<<endl;

    TString effName("TPCeff.");
    effName+=j;
    effName.Append(".root");
       
    if(gSystem->AccessPathName(effName.Data(),kFileExists)) continue;

    ReadEffs(effName.Data());

    for(Int_t k=0; k<n; k++) {
      if(fEffPi.GetParam(k)>=0.) {avEffPi[k]+=fEffPi.GetParam(k); evtsPi[k]++;}
      if(fEffKa.GetParam(k)>=0.) {avEffKa[k]+=fEffKa.GetParam(k); evtsKa[k]++;}
      if(fEffPr.GetParam(k)>=0.) {avEffPr[k]+=fEffPr.GetParam(k); evtsPr[k]++;}
      if(fEffEl.GetParam(k)>=0.) {avEffEl[k]+=fEffEl.GetParam(k); evtsEl[k]++;}
      if(fEffMu.GetParam(k)>=0.) {avEffMu[k]+=fEffMu.GetParam(k); evtsMu[k]++;}
    }

  } // loop on events

  // compute average efficiencies
  for(Int_t j=0; j<n; j++) {
    if(evtsPi[j]==0) evtsPi[j]++;
    fEffPi.SetParam(j,(Double_t)avEffPi[j]/evtsPi[j]);
    if(evtsKa[j]==0) evtsKa[j]++;
    fEffKa.SetParam(j,(Double_t)avEffKa[j]/evtsKa[j]);
    if(evtsPr[j]==0) evtsPr[j]++;
    fEffPr.SetParam(j,(Double_t)avEffPr[j]/evtsPr[j]);
    if(evtsEl[j]==0) evtsEl[j]++;
    fEffEl.SetParam(j,(Double_t)avEffEl[j]/evtsEl[j]);
    if(evtsMu[j]==0) evtsMu[j]++;
    fEffMu.SetParam(j,(Double_t)avEffMu[j]/evtsMu[j]);
  }
 
  // write efficiencies to a file
  WriteEffs(outName);

  delete [] avEffPi;
  delete [] avEffKa;
  delete [] avEffPr;
  delete [] avEffEl;
  delete [] avEffMu;
  delete [] evtsPi;
  delete [] evtsKa;
  delete [] evtsPr;
  delete [] evtsEl;
  delete [] evtsMu;

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
  if(!ReadRegParams(inName,2212)) return 0; 
  if(!ReadRegParams(inName,11)) return 0; 
  if(!ReadRegParams(inName,13)) return 0; 
  if(!ReaddEdx(inName,211)) return 0;
  if(!ReaddEdx(inName,321)) return 0;
  if(!ReaddEdx(inName,2212)) return 0;
  if(!ReaddEdx(inName,11)) return 0;
  if(!ReaddEdx(inName,13)) return 0;
  if(!ReadDBgrid(inName)) return 0;
  if(!ReadCovTrees(inName)) return 0;

  return 1;
}
//-----------------------------------------------------------------------------
Int_t AliTPCtrackerParam::ReadCovTrees(const Char_t* inName) {
//-----------------------------------------------------------------------------
// This function reads the covariance matrix trees from the DB
//-----------------------------------------------------------------------------

  if(gSystem->AccessPathName(inName,kFileExists)) { 
    cerr<<"AliTPCtrackerParam::ReaddEdx: "<<inName<<" not found\n"; 
    return 0; }

  TString str;

  TFile *inFile = TFile::Open(inName);


  Int_t nBinsPi = fDBgridPi.GetTotBins();
  for(Int_t l=0; l<nBinsPi; l++) {
    str = "/CovMatrices/Pions/CovTreePi_bin";
    str+=l;
    fCovTreePi[l] = (TTree*)inFile->Get(str.Data());
    //fCovTree = (TTree*)inFile->Get(str.Data());
    //fCovTreePi[l] = new TTree(*fCovTree);
  }
  Int_t nBinsKa = fDBgridKa.GetTotBins();
  for(Int_t l=0; l<nBinsKa; l++) {
    str = "/CovMatrices/Kaons/CovTreeKa_bin";
    str+=l;
    fCovTreeKa[l] = (TTree*)inFile->Get(str.Data());
    //fCovTree = (TTree*)inFile->Get(str.Data());
    //fCovTreeKa[l] = new TTree(*fCovTree);
  }
  Int_t nBinsPr = fDBgridPr.GetTotBins();
  for(Int_t l=0; l<nBinsPr; l++) {
    if(fColl==0) {      
      str = "/CovMatrices/Pions/CovTreePi_bin";
    } else {
      str = "/CovMatrices/Protons/CovTreePr_bin";
    }
    str+=l;
    fCovTreePr[l] = (TTree*)inFile->Get(str.Data());
    //fCovTree = (TTree*)inFile->Get(str.Data());
    //fCovTreePr[l] = new TTree(*fCovTree);
  }
  Int_t nBinsEl = fDBgridEl.GetTotBins();
  for(Int_t l=0; l<nBinsEl; l++) {
    str = "/CovMatrices/Electrons/CovTreeEl_bin";
    str+=l;
    fCovTreeEl[l] = (TTree*)inFile->Get(str.Data());
    //fCovTree = (TTree*)inFile->Get(str.Data());
    //fCovTreeEl[l] = new TTree(*fCovTree);
  }
  Int_t nBinsMu = fDBgridMu.GetTotBins();
  for(Int_t l=0; l<nBinsMu; l++) {
    if(fColl==0) {      
      str = "/CovMatrices/Pions/CovTreePi_bin";
    } else {
      str = "/CovMatrices/Muons/CovTreeMu_bin";
    }
    str+=l;
    fCovTreeMu[l] = (TTree*)inFile->Get(str.Data());
    //fCovTree = (TTree*)inFile->Get(str.Data());
    //fCovTreeMu[l] = new TTree(*fCovTree);
  }

  //inFile->Close();

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
  case 13:
    if(fColl==0) {
    fdEdxMean = (AliTPCkineGrid*)inFile->Get("/dEdx/Pions/dEdxMeanPi");
    fdEdxRMS = (AliTPCkineGrid*)inFile->Get("/dEdx/Pions/dEdxRMSPi");
    } else {
      fdEdxMean = (AliTPCkineGrid*)inFile->Get("/dEdx/Muons/dEdxMeanMu");
      fdEdxRMS = (AliTPCkineGrid*)inFile->Get("/dEdx/Muons/dEdxRMSMu");
    }
    fdEdxMeanMu.~AliTPCkineGrid();
    new(&fdEdxMeanMu) AliTPCkineGrid(*fdEdxMean);
    fdEdxRMSMu.~AliTPCkineGrid();
    new(&fdEdxRMSMu) AliTPCkineGrid(*fdEdxRMS);
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
  if(fColl==0) {
    fDBgrid = (AliTPCkineGrid*)inFile->Get("/CovMatrices/Pions/DBgridPi");
  } else {
    fDBgrid = (AliTPCkineGrid*)inFile->Get("/CovMatrices/Protons/DBgridPr");
  }
  fDBgridPr.~AliTPCkineGrid();
  new(&fDBgridPr) AliTPCkineGrid(*fDBgrid);
  fDBgrid = (AliTPCkineGrid*)inFile->Get("/CovMatrices/Electrons/DBgridEl");
  fDBgridEl.~AliTPCkineGrid();      
  new(&fDBgridEl) AliTPCkineGrid(*fDBgrid);
  if(fColl==0) {
    fDBgrid = (AliTPCkineGrid*)inFile->Get("/CovMatrices/Pions/DBgridPi");
  } else {
    fDBgrid = (AliTPCkineGrid*)inFile->Get("/CovMatrices/Muons/DBgridMu");
  }
  fDBgridMu.~AliTPCkineGrid();
  new(&fDBgridMu) AliTPCkineGrid(*fDBgrid);

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
    TString pr("/Pulls/Protons/PullsPr_"); pr+=i;
    TString el("/Pulls/Electrons/PullsEl_"); el+=i;
    TString mu("/Pulls/Muons/PullsMu_"); mu+=i;

    fPulls = (AliTPCkineGrid*)inFile->Get(pi.Data());
    fPullsPi[i].~AliTPCkineGrid();
    new(&fPullsPi[i]) AliTPCkineGrid(*fPulls);

    fPulls = (AliTPCkineGrid*)inFile->Get(ka.Data());
    fPullsKa[i].~AliTPCkineGrid();
    new(&fPullsKa[i]) AliTPCkineGrid(*fPulls);

    if(fColl==0) {
      fPulls = (AliTPCkineGrid*)inFile->Get(pi.Data());
    } else {
      fPulls = (AliTPCkineGrid*)inFile->Get(pr.Data());
    }
    fPullsPr[i].~AliTPCkineGrid();
    new(&fPullsPr[i]) AliTPCkineGrid(*fPulls);

    fPulls = (AliTPCkineGrid*)inFile->Get(el.Data());
    fPullsEl[i].~AliTPCkineGrid();
    new(&fPullsEl[i]) AliTPCkineGrid(*fPulls);

    if(fColl==0) {
      fPulls = (AliTPCkineGrid*)inFile->Get(pi.Data());
    } else {
      fPulls = (AliTPCkineGrid*)inFile->Get(mu.Data());
    }
    fPullsMu[i].~AliTPCkineGrid();
    new(&fPullsMu[i]) AliTPCkineGrid(*fPulls);
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
  // Introduced change to "reverse" the TMatrixD read from file.
  // Needed because storage mode of TMatrixD changed from column-wise
  // to rwo-wise in ROOT.
  //
  // A.Dainese 03/06/05 

  TMatrixD dummyMatrix(9,3);

  TFile *inFile = TFile::Open(inName);
  switch (pdg) {
  case 211:
    fRegPar = (TMatrixD*)inFile->Get("/RegParams/Pions/RegPions");
    new(&fRegParPi) TMatrixD(*fRegPar);
    //printf("before\n");
    //for(Int_t jj=0;jj<9;jj++) printf("%13.10f   %13.10f  %13.10f\n",fRegParPi(jj,0),fRegParPi(jj,1),fRegParPi(jj,2));
    dummyMatrix = fRegParPi;  
    fRegParPi(0,0) = dummyMatrix(0,0);
    fRegParPi(0,1) = dummyMatrix(0,1);
    fRegParPi(0,2) = dummyMatrix(0,2);
    fRegParPi(1,0) = dummyMatrix(3,0);
    fRegParPi(1,1) = dummyMatrix(1,1);
    fRegParPi(1,2) = dummyMatrix(1,2);
    fRegParPi(2,0) = dummyMatrix(6,0);
    fRegParPi(2,1) = dummyMatrix(3,2);
    fRegParPi(2,2) = dummyMatrix(2,2);
    fRegParPi(3,0) = dummyMatrix(1,0);
    fRegParPi(3,1) = dummyMatrix(4,0);
    fRegParPi(3,2) = dummyMatrix(7,0);
    fRegParPi(4,0) = dummyMatrix(3,1);
    fRegParPi(4,1) = dummyMatrix(4,1);
    fRegParPi(4,2) = dummyMatrix(7,1);
    fRegParPi(5,0) = dummyMatrix(6,1);
    fRegParPi(5,1) = dummyMatrix(4,2);
    fRegParPi(5,2) = dummyMatrix(7,2);
    fRegParPi(6,0) = dummyMatrix(2,0);
    fRegParPi(6,1) = dummyMatrix(5,0);
    fRegParPi(6,2) = dummyMatrix(8,0);
    fRegParPi(7,0) = dummyMatrix(2,1);
    fRegParPi(7,1) = dummyMatrix(5,1);
    fRegParPi(7,2) = dummyMatrix(8,1);
    fRegParPi(8,0) = dummyMatrix(6,2);
    fRegParPi(8,1) = dummyMatrix(5,2);
    fRegParPi(8,2) = dummyMatrix(8,2);
    //printf("after\n");
    //for(Int_t jj=0;jj<9;jj++) printf("%13.10f   %13.10f  %13.10f\n",fRegParPi(jj,0),fRegParPi(jj,1),fRegParPi(jj,2));
    break;
  case 321:
    fRegPar = (TMatrixD*)inFile->Get("/RegParams/Kaons/RegKaons");
    new(&fRegParKa) TMatrixD(*fRegPar);
    dummyMatrix = fRegParKa;  
    fRegParKa(0,0) = dummyMatrix(0,0);
    fRegParKa(0,1) = dummyMatrix(0,1);
    fRegParKa(0,2) = dummyMatrix(0,2);
    fRegParKa(1,0) = dummyMatrix(3,0);
    fRegParKa(1,1) = dummyMatrix(1,1);
    fRegParKa(1,2) = dummyMatrix(1,2);
    fRegParKa(2,0) = dummyMatrix(6,0);
    fRegParKa(2,1) = dummyMatrix(3,2);
    fRegParKa(2,2) = dummyMatrix(2,2);
    fRegParKa(3,0) = dummyMatrix(1,0);
    fRegParKa(3,1) = dummyMatrix(4,0);
    fRegParKa(3,2) = dummyMatrix(7,0);
    fRegParKa(4,0) = dummyMatrix(3,1);
    fRegParKa(4,1) = dummyMatrix(4,1);
    fRegParKa(4,2) = dummyMatrix(7,1);
    fRegParKa(5,0) = dummyMatrix(6,1);
    fRegParKa(5,1) = dummyMatrix(4,2);
    fRegParKa(5,2) = dummyMatrix(7,2);
    fRegParKa(6,0) = dummyMatrix(2,0);
    fRegParKa(6,1) = dummyMatrix(5,0);
    fRegParKa(6,2) = dummyMatrix(8,0);
    fRegParKa(7,0) = dummyMatrix(2,1);
    fRegParKa(7,1) = dummyMatrix(5,1);
    fRegParKa(7,2) = dummyMatrix(8,1);
    fRegParKa(8,0) = dummyMatrix(6,2);
    fRegParKa(8,1) = dummyMatrix(5,2);
    fRegParKa(8,2) = dummyMatrix(8,2);
    break;
  case 2212:
    if(fColl==0) {
      fRegPar = (TMatrixD*)inFile->Get("/RegParams/Pions/RegPions");
    } else {
      fRegPar = (TMatrixD*)inFile->Get("/RegParams/Protons/RegProtons");
    }
    new(&fRegParPr) TMatrixD(*fRegPar);
    dummyMatrix = fRegParPr;  
    fRegParPr(0,0) = dummyMatrix(0,0);
    fRegParPr(0,1) = dummyMatrix(0,1);
    fRegParPr(0,2) = dummyMatrix(0,2);
    fRegParPr(1,0) = dummyMatrix(3,0);
    fRegParPr(1,1) = dummyMatrix(1,1);
    fRegParPr(1,2) = dummyMatrix(1,2);
    fRegParPr(2,0) = dummyMatrix(6,0);
    fRegParPr(2,1) = dummyMatrix(3,2);
    fRegParPr(2,2) = dummyMatrix(2,2);
    fRegParPr(3,0) = dummyMatrix(1,0);
    fRegParPr(3,1) = dummyMatrix(4,0);
    fRegParPr(3,2) = dummyMatrix(7,0);
    fRegParPr(4,0) = dummyMatrix(3,1);
    fRegParPr(4,1) = dummyMatrix(4,1);
    fRegParPr(4,2) = dummyMatrix(7,1);
    fRegParPr(5,0) = dummyMatrix(6,1);
    fRegParPr(5,1) = dummyMatrix(4,2);
    fRegParPr(5,2) = dummyMatrix(7,2);
    fRegParPr(6,0) = dummyMatrix(2,0);
    fRegParPr(6,1) = dummyMatrix(5,0);
    fRegParPr(6,2) = dummyMatrix(8,0);
    fRegParPr(7,0) = dummyMatrix(2,1);
    fRegParPr(7,1) = dummyMatrix(5,1);
    fRegParPr(7,2) = dummyMatrix(8,1);
    fRegParPr(8,0) = dummyMatrix(6,2);
    fRegParPr(8,1) = dummyMatrix(5,2);
    fRegParPr(8,2) = dummyMatrix(8,2);
    break;
  case 11:
    fRegPar = (TMatrixD*)inFile->Get("/RegParams/Electrons/RegElectrons");
    new(&fRegParEl) TMatrixD(*fRegPar);
    dummyMatrix = fRegParEl;  
    fRegParEl(0,0) = dummyMatrix(0,0);
    fRegParEl(0,1) = dummyMatrix(0,1);
    fRegParEl(0,2) = dummyMatrix(0,2);
    fRegParEl(1,0) = dummyMatrix(3,0);
    fRegParEl(1,1) = dummyMatrix(1,1);
    fRegParEl(1,2) = dummyMatrix(1,2);
    fRegParEl(2,0) = dummyMatrix(6,0);
    fRegParEl(2,1) = dummyMatrix(3,2);
    fRegParEl(2,2) = dummyMatrix(2,2);
    fRegParEl(3,0) = dummyMatrix(1,0);
    fRegParEl(3,1) = dummyMatrix(4,0);
    fRegParEl(3,2) = dummyMatrix(7,0);
    fRegParEl(4,0) = dummyMatrix(3,1);
    fRegParEl(4,1) = dummyMatrix(4,1);
    fRegParEl(4,2) = dummyMatrix(7,1);
    fRegParEl(5,0) = dummyMatrix(6,1);
    fRegParEl(5,1) = dummyMatrix(4,2);
    fRegParEl(5,2) = dummyMatrix(7,2);
    fRegParEl(6,0) = dummyMatrix(2,0);
    fRegParEl(6,1) = dummyMatrix(5,0);
    fRegParEl(6,2) = dummyMatrix(8,0);
    fRegParEl(7,0) = dummyMatrix(2,1);
    fRegParEl(7,1) = dummyMatrix(5,1);
    fRegParEl(7,2) = dummyMatrix(8,1);
    fRegParEl(8,0) = dummyMatrix(6,2);
    fRegParEl(8,1) = dummyMatrix(5,2);
    fRegParEl(8,2) = dummyMatrix(8,2);
    break;
  case 13:
    if(fColl==0) {
      fRegPar = (TMatrixD*)inFile->Get("/RegParams/Pions/RegPions");
    } else {
      fRegPar = (TMatrixD*)inFile->Get("/RegParams/Muons/RegMuons");
    }
    new(&fRegParMu) TMatrixD(*fRegPar);
    dummyMatrix = fRegParMu;  
    fRegParMu(0,0) = dummyMatrix(0,0);
    fRegParMu(0,1) = dummyMatrix(0,1);
    fRegParMu(0,2) = dummyMatrix(0,2);
    fRegParMu(1,0) = dummyMatrix(3,0);
    fRegParMu(1,1) = dummyMatrix(1,1);
    fRegParMu(1,2) = dummyMatrix(1,2);
    fRegParMu(2,0) = dummyMatrix(6,0);
    fRegParMu(2,1) = dummyMatrix(3,2);
    fRegParMu(2,2) = dummyMatrix(2,2);
    fRegParMu(3,0) = dummyMatrix(1,0);
    fRegParMu(3,1) = dummyMatrix(4,0);
    fRegParMu(3,2) = dummyMatrix(7,0);
    fRegParMu(4,0) = dummyMatrix(3,1);
    fRegParMu(4,1) = dummyMatrix(4,1);
    fRegParMu(4,2) = dummyMatrix(7,1);
    fRegParMu(5,0) = dummyMatrix(6,1);
    fRegParMu(5,1) = dummyMatrix(4,2);
    fRegParMu(5,2) = dummyMatrix(7,2);
    fRegParMu(6,0) = dummyMatrix(2,0);
    fRegParMu(6,1) = dummyMatrix(5,0);
    fRegParMu(6,2) = dummyMatrix(8,0);
    fRegParMu(7,0) = dummyMatrix(2,1);
    fRegParMu(7,1) = dummyMatrix(5,1);
    fRegParMu(7,2) = dummyMatrix(8,1);
    fRegParMu(8,0) = dummyMatrix(6,2);
    fRegParMu(8,1) = dummyMatrix(5,2);
    fRegParMu(8,2) = dummyMatrix(8,2);
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
  const char *part="Pions - ";

  InitializeKineGrid("DB");
  SetParticle(pdg);
  const Int_t kfitbins = fDBgrid->GetBinsPt();
  cerr<<" Fit bins:  "<<kfitbins<<endl;

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
  case 2212: //protons
    thisPdg=2212; 
    part="Protons - ";
    cerr<<" Processing protons ...\n";
    break;
  case 11: // electrons
    thisPdg= 11; 
    part="Electrons - ";
    cerr<<" Processing electrons ...\n";
    break;
  case 13: // muons
    thisPdg= 13; 
    part="Muons - ";
    cerr<<" Processing muons ...\n";
    break;
  }


  /*
  // create a chain with compared tracks
  TChain *cmptrkchain = new ("cmptrktree");
  cmptrkchain.Add("CovMatrix_AllEvts.root");
  //cmptrkchain.Add("CovMatrix_AllEvts_1.root");
  //cmptrkchain.Add("CovMatrix_AllEvts_2.root");
  //cmptrkchain.Add("CovMatrix_AllEvts_3.root");
  */


  TFile *infile = TFile::Open("CovMatrix_AllEvts.root");
  TTree *cmptrktree = (TTree*)infile->Get("cmptrktree");

  COMPTRACK cmptrk; 
  cmptrktree->SetBranchAddress("comptracks",&cmptrk);
  Int_t entries = (Int_t)cmptrktree->GetEntries(); 

    
  Int_t pbin;
  Int_t    *n       = new Int_t[kfitbins];
  Int_t    *n00     = new Int_t[kfitbins];
  Int_t    *n11     = new Int_t[kfitbins];
  Int_t    *n20     = new Int_t[kfitbins];
  Int_t    *n22     = new Int_t[kfitbins];
  Int_t    *n31     = new Int_t[kfitbins];
  Int_t    *n33     = new Int_t[kfitbins];
  Int_t    *n40     = new Int_t[kfitbins];
  Int_t    *n42     = new Int_t[kfitbins];
  Int_t    *n44     = new Int_t[kfitbins];
  Double_t *p       = new Double_t[kfitbins];
  Double_t *ep      = new Double_t[kfitbins];
  Double_t *mean00  = new Double_t[kfitbins];
  Double_t *mean11  = new Double_t[kfitbins];
  Double_t *mean20  = new Double_t[kfitbins];
  Double_t *mean22  = new Double_t[kfitbins];
  Double_t *mean31  = new Double_t[kfitbins];
  Double_t *mean33  = new Double_t[kfitbins];
  Double_t *mean40  = new Double_t[kfitbins];
  Double_t *mean42  = new Double_t[kfitbins];
  Double_t *mean44  = new Double_t[kfitbins];
  Double_t *sigma00 = new Double_t[kfitbins];
  Double_t *sigma11 = new Double_t[kfitbins];
  Double_t *sigma20 = new Double_t[kfitbins];
  Double_t *sigma22 = new Double_t[kfitbins];
  Double_t *sigma31 = new Double_t[kfitbins];
  Double_t *sigma33 = new Double_t[kfitbins];
  Double_t *sigma40 = new Double_t[kfitbins];
  Double_t *sigma42 = new Double_t[kfitbins];
  Double_t *sigma44 = new Double_t[kfitbins];
  Double_t *rmean   = new Double_t[kfitbins];
  Double_t *rsigma  = new Double_t[kfitbins];
  Double_t fitpar[3];

  for(Int_t l=0; l<kfitbins; l++) {
    n[l]=1;
    n00[l]=n11[l]=n20[l]=n22[l]=n31[l]=n33[l]=n40[l]=n42[l]=n44[l]=1;
    p[l ]=ep[l]=0.;
    mean00[l]=mean11[l]=mean20[l]=mean22[l]=mean31[l]=mean33[l]=mean40[l]=mean42[l]=mean44[l]=0.;
    sigma00[l]=sigma11[l]=sigma20[l]=sigma22[l]=sigma31[l]=sigma33[l]=sigma40[l]=sigma42[l]=sigma44[l]=0.;
  }

  // loop on chain entries for mean 
  for(Int_t l=0; l<entries; l++) {
    cmptrktree->GetEvent(l);
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

  for(Int_t l=0; l<kfitbins; l++) {
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
    cmptrktree->GetEvent(l);
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
 
  for(Int_t l=0; l<kfitbins; l++) {
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
  TF1 *func = new TF1("RegFunc",RegFunc,0.23,50.,3);
  func->SetParNames("A_meas","A_scatt","B");

  // line to draw on the plots
  TLine *lin = new TLine(-1,1,1.69,1);
  lin->SetLineStyle(2);
  lin->SetLineWidth(2);

  // matrix used to store fit results
  TMatrixD fitRes(9,3);

  //    --- c00 ---

  // create the canvas
  TCanvas *canv00 = new TCanvas("canv00","c00",0,0,700,900); 
  canv00->Divide(1,2);
  // create the graph for cov matrix
  TGraphErrors *gr00 = new TGraphErrors(kfitbins,p,mean00,ep,sigma00);
  TString title00("C(y,y)"); title00.Prepend(part);
  TH2F *frame00 = new TH2F("frame00",title00.Data(),2,0.1,50,2,0,5e-3);
  frame00->SetXTitle("p [GeV/c]");
  canv00->cd(1);  gPad->SetLogx();
  frame00->Draw();
  gr00->Draw("P");
  // Sets initial values for parameters
  func->SetParameters(1.6e-3,1.9e-4,1.5);
  // Fit points in range defined by function
  gr00->Fit("RegFunc","R,Q");
  func->GetParameters(fitpar);
  for(Int_t i=0; i<3; i++) fitRes(0,i)=fitpar[i];
  for(Int_t l=0; l<kfitbins; l++) {
    rmean[l]  = mean00[l]/RegFunc(&p[l],fitpar);
    rsigma[l] = sigma00[l]/RegFunc(&p[l],fitpar);
  }
  // create the graph the regularized cov. matrix
  TGraphErrors *gr00reg = new TGraphErrors(kfitbins,p,rmean,ep,rsigma);
  TString regtitle00("C(y,y)/(A_meas+A_scatt/p^{B})"); 
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
  TGraphErrors *gr11 = new TGraphErrors(kfitbins,p,mean11,ep,sigma11);
  TString title11("C(z,z)"); title11.Prepend(part);
  TH2F *frame11 = new TH2F("frame11",title11.Data(),2,0.1,50,2,0,6e-3);
  frame11->SetXTitle("p [GeV/c]");
  canv11->cd(1);  gPad->SetLogx();
  frame11->Draw();
  gr11->Draw("P");
  // Sets initial values for parameters
  func->SetParameters(1.2e-3,8.1e-4,1.);
  // Fit points in range defined by function
  gr11->Fit("RegFunc","R,Q");
  func->GetParameters(fitpar);
  for(Int_t i=0; i<3; i++) fitRes(1,i)=fitpar[i];
  for(Int_t l=0; l<kfitbins; l++) {
    rmean[l]  = mean11[l]/RegFunc(&p[l],fitpar);
    rsigma[l] = sigma11[l]/RegFunc(&p[l],fitpar);
  }
  // create the graph the regularized cov. matrix
  TGraphErrors *gr11reg = new TGraphErrors(kfitbins,p,rmean,ep,rsigma);
  TString regtitle11("C(z,z)/(A_meas+A_scatt/p^{B})"); 
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
  TGraphErrors *gr20 = new TGraphErrors(kfitbins,p,mean20,ep,sigma20);
  TString title20("C(#eta, y)"); title20.Prepend(part);
  TH2F *frame20 = new TH2F("frame20",title20.Data(),2,0.1,50,2,0,2.5e-4);
  frame20->SetXTitle("p [GeV/c]");
  canv20->cd(1);  gPad->SetLogx();
  frame20->Draw();
  gr20->Draw("P");
  // Sets initial values for parameters
  func->SetParameters(7.3e-5,1.2e-5,1.5);
  // Fit points in range defined by function
  gr20->Fit("RegFunc","R,Q");
  func->GetParameters(fitpar);
  for(Int_t i=0; i<3; i++) fitRes(2,i)=fitpar[i];
  for(Int_t l=0; l<kfitbins; l++) {
    rmean[l]  = mean20[l]/RegFunc(&p[l],fitpar);
    rsigma[l] = sigma20[l]/RegFunc(&p[l],fitpar);
  }
  // create the graph the regularized cov. matrix
  TGraphErrors *gr20reg = new TGraphErrors(kfitbins,p,rmean,ep,rsigma);
  TString regtitle20("C(#eta, y)/(A_meas+A_scatt/p^{B})"); 
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
  TGraphErrors *gr22 = new TGraphErrors(kfitbins,p,mean22,ep,sigma22);
  TString title22("C(#eta, #eta)"); title22.Prepend(part);
  TH2F *frame22 = new TH2F("frame22",title22.Data(),2,0.1,50,2,0,3e-5);
  frame22->SetXTitle("p [GeV/c]");
  canv22->cd(1);  gPad->SetLogx();
  frame22->Draw();
  gr22->Draw("P");
  // Sets initial values for parameters
  func->SetParameters(5.2e-6,1.1e-6,2.);
  // Fit points in range defined by function
  gr22->Fit("RegFunc","R,Q");
  func->GetParameters(fitpar);
  for(Int_t i=0; i<3; i++) fitRes(3,i)=fitpar[i];
  for(Int_t l=0; l<kfitbins; l++) {
    rmean[l]  = mean22[l]/RegFunc(&p[l],fitpar);
    rsigma[l] = sigma22[l]/RegFunc(&p[l],fitpar);
  }
  // create the graph the regularized cov. matrix
  TGraphErrors *gr22reg = new TGraphErrors(kfitbins,p,rmean,ep,rsigma);
  TString regtitle22("C(#eta, #eta)/(A_meas+A_scatt/p^{B})"); 
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
  TGraphErrors *gr31 = new TGraphErrors(kfitbins,p,mean31,ep,sigma31);
  TString title31("C(tg #lambda,z)"); title31.Prepend(part);
  TH2F *frame31 = new TH2F("frame31",title31.Data(),2,0.1,50,2,-2e-4,0);
  frame31->SetXTitle("p [GeV/c]");
  canv31->cd(1);  gPad->SetLogx();
  frame31->Draw();
  gr31->Draw("P");
  // Sets initial values for parameters
  func->SetParameters(-1.2e-5,-1.2e-5,1.5);
  // Fit points in range defined by function
  gr31->Fit("RegFunc","R,Q");
  func->GetParameters(fitpar);
  for(Int_t i=0; i<3; i++) fitRes(4,i)=fitpar[i];
  for(Int_t l=0; l<kfitbins; l++) {
    rmean[l]  = mean31[l]/RegFunc(&p[l],fitpar);
    rsigma[l] = -sigma31[l]/RegFunc(&p[l],fitpar);
  }
  // create the graph the regularized cov. matrix
  TGraphErrors *gr31reg = new TGraphErrors(kfitbins,p,rmean,ep,rsigma);
  TString regtitle31("C(tg #lambda,z)/(A_meas+A_scatt/p^{B})"); 
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
  TGraphErrors *gr33 = new TGraphErrors(kfitbins,p,mean33,ep,sigma33);
  TString title33("C(tg #lambda,tg #lambda)"); title33.Prepend(part);
  TH2F *frame33 = new TH2F("frame33",title33.Data(),2,0.1,50,2,0,1e-5);
  frame33->SetXTitle("p [GeV/c]");
  canv33->cd(1);  gPad->SetLogx();
  frame33->Draw();
  gr33->Draw("P");
  // Sets initial values for parameters
  func->SetParameters(1.3e-7,4.6e-7,1.7);
  // Fit points in range defined by function
  gr33->Fit("RegFunc","R,Q");
  func->GetParameters(fitpar);
  for(Int_t i=0; i<3; i++) fitRes(5,i)=fitpar[i];
  for(Int_t l=0; l<kfitbins; l++) {
    rmean[l]  = mean33[l]/RegFunc(&p[l],fitpar);
    rsigma[l] = sigma33[l]/RegFunc(&p[l],fitpar);
  }
  // create the graph the regularized cov. matrix
  TGraphErrors *gr33reg = new TGraphErrors(kfitbins,p,rmean,ep,rsigma);
  TString regtitle33("C(tg #lambda,tg #lambda)/(A_meas+A_scatt/p^{B})"); 
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
  TGraphErrors *gr40 = new TGraphErrors(kfitbins,p,mean40,ep,sigma40);
  TString title40("C(C,y)"); title40.Prepend(part);
  TH2F *frame40 = new TH2F("frame40",title40.Data(),2,0.1,50,2,0,1e-6);
  frame40->SetXTitle("p [GeV/c]");
  canv40->cd(1);  gPad->SetLogx();
  frame40->Draw();
  gr40->Draw("P");
  // Sets initial values for parameters
  func->SetParameters(4.e-7,4.4e-8,1.5);
  // Fit points in range defined by function
  gr40->Fit("RegFunc","R,Q");
  func->GetParameters(fitpar);
  for(Int_t i=0; i<3; i++) fitRes(6,i)=fitpar[i];
  for(Int_t l=0; l<kfitbins; l++) {
    rmean[l]  = mean40[l]/RegFunc(&p[l],fitpar);
    rsigma[l] = sigma40[l]/RegFunc(&p[l],fitpar);
  }
  // create the graph the regularized cov. matrix
  TGraphErrors *gr40reg = new TGraphErrors(kfitbins,p,rmean,ep,rsigma);
  TString regtitle40("C(C,y)/(A_meas+A_scatt/p^{B})"); 
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
  TGraphErrors *gr42 = new TGraphErrors(kfitbins,p,mean42,ep,sigma42);
  TString title42("C(C, #eta)"); title42.Prepend(part);
  TH2F *frame42 = new TH2F("frame42",title42.Data(),2,0.1,50,2,0,2.2e-7);
  frame42->SetXTitle("p [GeV/c]");
  canv42->cd(1);  gPad->SetLogx();
  frame42->Draw();
  gr42->Draw("P");
  // Sets initial values for parameters
  func->SetParameters(3.e-8,8.2e-9,2.);
  // Fit points in range defined by function
  gr42->Fit("RegFunc","R,Q");
  func->GetParameters(fitpar);
  for(Int_t i=0; i<3; i++) fitRes(7,i)=fitpar[i];
  for(Int_t l=0; l<kfitbins; l++) {
    rmean[l]  = mean42[l]/RegFunc(&p[l],fitpar);
    rsigma[l] = sigma42[l]/RegFunc(&p[l],fitpar);
  }
  // create the graph the regularized cov. matrix
  TGraphErrors *gr42reg = new TGraphErrors(kfitbins,p,rmean,ep,rsigma);
  TString regtitle42("C(C, #eta)/(A_meas+A_scatt/p^{B})"); 
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
  TGraphErrors *gr44 = new TGraphErrors(kfitbins,p,mean44,ep,sigma44);
  TString title44("C(C,C)"); title44.Prepend(part);
  TH2F *frame44 = new TH2F("frame44",title44.Data(),2,0.1,50,2,0,2e-9);
  frame44->SetXTitle("p [GeV/c]");
  canv44->cd(1);  gPad->SetLogx();
  frame44->Draw();
  gr44->Draw("P");
  // Sets initial values for parameters
  func->SetParameters(1.8e-10,5.8e-11,2.);
  // Fit points in range defined by function
  gr44->Fit("RegFunc","R,Q");
  func->GetParameters(fitpar);
  for(Int_t i=0; i<3; i++) fitRes(8,i)=fitpar[i];
  for(Int_t l=0; l<kfitbins; l++) {
    rmean[l]  = mean44[l]/RegFunc(&p[l],fitpar);
    rsigma[l] = sigma44[l]/RegFunc(&p[l],fitpar);
  }
  // create the graph the regularized cov. matrix
  TGraphErrors *gr44reg = new TGraphErrors(kfitbins,p,rmean,ep,rsigma);
  TString regtitle44("C(C,C)/(A_meas+A_scatt/p^{B})"); 
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
  case 2212:
    new(&fRegParPr) TMatrixD(fitRes);
    break;
  case 11:
    new(&fRegParEl) TMatrixD(fitRes);
    break;
  case 13:
    new(&fRegParMu) TMatrixD(fitRes);
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
    fDBgrid = &fDBgridPr;
    fEff    = &fEffPr;
    fPulls  =  fPullsPr;
    fRegPar = &fRegParPr;
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
    fDBgrid = &fDBgridMu;
    fEff    = &fEffMu;
    fPulls  =  fPullsMu;
    fRegPar = &fRegParMu;
    fdEdxMean = &fdEdxMeanMu;
    fdEdxRMS  = &fdEdxRMSMu;
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
  Double_t xref=xxsm[0]; xxsm[0]=0;

  AliGausCorr *corgen = new AliGausCorr(cov,5);
  TArrayD corr(5);
  corgen->GetGaussN(corr);
  // check on fP2(ESD) = fX*fP4-fP2(TPC)
  // with fX=xref (not smeared), fP4=xx[4]+corr[4] e fP2=xx[2]+corr[2];
  // if fP2(ESD)^2 > 1 -> resmear...
  Double_t fp2esd=xref*(xx[4]+corr[4])-(xx[2]+corr[2]);
  while ( (fp2esd*fp2esd) > 1.0 ) {
    Warning("SmearTrack()","Badly smeared track, retrying...");
    corgen->GetGaussN(corr);
    fp2esd=xref*(xx[4]+corr[4])-(xx[2]+corr[2]);
  }

  delete corgen;
  corgen = 0;

  for(Int_t l=0;l<5;l++) xxsm[l] = xx[l]+corr[l];

  return;
}
//-----------------------------------------------------------------------------
Int_t AliTPCtrackerParam::WritedEdx(const Char_t *outName,Int_t pdg) {
//-----------------------------------------------------------------------------
// This function writes the dEdx parameters to the DB
//-----------------------------------------------------------------------------

  Option_t *opt;
  const char *dirName="Pions";
  const char *meanName="dEdxMeanPi";
  const char *rmsName="dEdxRMSPi";

  SetParticle(pdg);

  if(gSystem->AccessPathName(outName,kFileExists))  { opt="recreate"; 
  } else { opt="update"; }

  switch (pdg) {
  case 211:
    dirName="Pions";
    meanName="dEdxMeanPi";
    rmsName="dEdxRMSPi";
    break;
  case 321:
    dirName="Kaons";
    meanName="dEdxMeanKa";
    rmsName="dEdxRMSKa";
    break;
  case 2212:
    dirName="Protons";
    meanName="dEdxMeanPr";
    rmsName="dEdxRMSPr";
    break;
  case 11:
    dirName="Electrons";
    meanName="dEdxMeanEl";
    rmsName="dEdxRMSEl";
    break;
  case 13:
    dirName="Muons";
    meanName="dEdxMeanMu";
    rmsName="dEdxRMSMu";
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
  fdEdxRMS->SetName(rmsName);  fdEdxRMS->Write();

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
  gDirectory->mkdir("Protons");
  gDirectory->mkdir("Electrons");
  gDirectory->mkdir("Muons");

  for(Int_t i=0;i<5;i++) {
    TString pi("PullsPi_"); pi+=i;
    TString ka("PullsKa_"); ka+=i;
    TString pr("PullsPr_"); pr+=i;
    TString el("PullsEl_"); el+=i;
    TString mu("PullsMu_"); mu+=i;
    fPullsPi[i].SetName(pi.Data());
    fPullsKa[i].SetName(ka.Data());
    fPullsPr[i].SetName(pr.Data());
    fPullsEl[i].SetName(el.Data());
    fPullsMu[i].SetName(mu.Data());
    gDirectory->cd("/Pulls/Pions");
    fPullsPi[i].Write();
    gDirectory->cd("/Pulls/Kaons");
    fPullsKa[i].Write();
    gDirectory->cd("/Pulls/Protons");
    fPullsPr[i].Write();
    gDirectory->cd("/Pulls/Electrons");
    fPullsEl[i].Write();
    gDirectory->cd("/Pulls/Muons");
    fPullsMu[i].Write();
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
  const char *dirName="Pions";
  const char *keyName="RegPions";

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
  case 2212:
    dirName="Protons";
    keyName="RegProtons";
    break;
  case 11:
    dirName="Electrons";
    keyName="RegElectrons";
    break;
  case 13:
    dirName="Muons";
    keyName="RegMuons";
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
//-----------------------------------------------------------------------------
//*****************************************************************************
//-----------------------------------------------------------------------------























