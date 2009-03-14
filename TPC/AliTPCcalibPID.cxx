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


Basic calibration and QA class for the PID of the TPC based on dEdx.

The corrections which are done on the cluster level (angle,drift time etc.) are here checked by plotting the dE/dx for selected tracks as a function of these variables. In addition the Bethe-Bloch-Parameterisation is fitted to the distribution and the corresponding parameters are being extracted.

1.a) Fast equalization for different gain in pad regions:

TFile f("CalibObjects.root")
AliTPCcalibPID *cal = f->Get("TPCCalib")->FindObject("calibPID")

cal->GetHistLongMediumRatio()->Projection(0)->Fit("gaus")
cal->GetHistShortMediumRatio()->Projection(0)->Fit("gaus")

1.b) Update OCDB:
.x $ALICE_ROOT/TPC/macros/ConfigOCDB.C
AliTPCClusterParam * parcl = AliTPCcalibDB::Instance()->GetClusterParam();
(*parcl->fQpadMnorm)[ipad] = oldvalue*corrFactor

Int_t runNumber = 0;
AliCDBMetaData *metaData= new AliCDBMetaData();
metaData->SetObjectClassName("AliTPCClusterParam");
metaData->SetResponsible("Alexander Kalweit");
metaData->SetBeamPeriod(1);
metaData->SetAliRootVersion("05-23-02"); 
metaData->SetComment("October runs calibration");
AliCDBId id1("TPC/Calib/ClusterParam", runNumber, AliCDBRunRange::Infinity());
gStorage = AliCDBManager::Instance()->GetStorage("local://$ALICE_ROOT/OCDB");
gStorage->Put(parcl, id1, metaData);
  

2.) Checks should be done with particles on the Fermi-Plateau and whoch leave a long track:

cal->GetHistQmax()->GetAxis(6)->SetRangeUser(120,160)
cal->GetHistQmax()->GetAxis(4)->SetRangeUser(20,50)

variables to check:

cal->GetHistQmax()->Projection(0,1)->Draw() z-dep.
cal->GetHistQmax()->Projection(0,2)->Draw() snp-dep.
cal->GetHistQmax()->Projection(0,3)->Draw() tgl-dep.



    Comments to be written here:
    1. What do we calibrate.
    2. How to interpret results
    3. Simple example
    4. Analysis using debug streamers.    

Double_t initialParam[] = {3.81470,15.2608,3.25306e-07,2.51791,2.71012}


Send comments etc. to: A.Kalweit@gsi.de, marian.ivanov@cern.ch
*/


#include "Riostream.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TList.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TF1.h"
#include "TVectorD.h"
#include "TProfile.h"

#include "AliTPCcalibDB.h"
#include "AliTPCclusterMI.h"
#include "AliTPCClusterParam.h"
#include "AliTPCseed.h"
#include "AliESDVertex.h"
#include "AliESDEvent.h"
#include "AliESDfriend.h"
#include "AliESDInputHandler.h"
#include "AliAnalysisManager.h"
#include "AliTPCParam.h"

#include "AliLog.h"

#include "AliTPCcalibPID.h"

#include "TTreeStream.h"


ClassImp(AliTPCcalibPID)


AliTPCcalibPID::AliTPCcalibPID() 
  :AliTPCcalibBase(),
   fMIP(0),
   fLowerTrunc(0),
   fUpperTrunc(0),
   fUseShapeNorm(0),
   fUsePosNorm(0),
   fUsePadNorm(0),
   fIsCosmic(0),  
   fHistNTracks(0),
   fClusters(0),   
   fPileUp(0),
   fLandau(0),
   fDeDxQmax(0),
   fDeDxQtot(0),
   fDeDxRatio(0),
   fDeDxShortMediumRatio(0),
   fDeDxLongMediumRatio(0)
{  
  //
  // Empty default cosntructor
  //
  AliInfo("Default Constructor");  
}


AliTPCcalibPID::AliTPCcalibPID(const Text_t *name, const Text_t *title) 
  :AliTPCcalibBase(),
   fMIP(0),
   fLowerTrunc(0),
   fUpperTrunc(0),
   fUseShapeNorm(0),
   fUsePosNorm(0),
   fUsePadNorm(0),
   fIsCosmic(0),  
   fHistNTracks(0),
   fClusters(0),   
   fPileUp(0),
   fLandau(0),
   fDeDxQmax(0),
   fDeDxQtot(0),
   fDeDxRatio(0),
   fDeDxShortMediumRatio(0),
   fDeDxLongMediumRatio(0)
{
  //
  //
  //  
  SetName(name);
  SetTitle(title);
  //
  fMIP = 40.;
  fLowerTrunc = 0;
  fUpperTrunc = 0.6;
  //
  fUseShapeNorm = kTRUE;
  fUsePosNorm = kFALSE;
  fUsePadNorm = 1;
  //
  fIsCosmic  = kTRUE;
  //
  //                  dE/dx,  z, phi, theta,    p,  bg, ncls
  Int_t binsQA[7]    = {150, 10,  10,    10,   50,   1,  8};
  Double_t xminQA[7] = {0.,  0,    0,     0, 0.01, 0.1, 60};
  Double_t xmaxQA[7] = {10., 1,  0.6,   1.5,   50, 500, 160};
  //
  Double_t xminRa[7] = {0.5, 0,    0,     0, 0.01, 0.1, 60};
  Double_t xmaxRa[7] = { 2., 1,  0.6,   1.5,   50, 500, 160};
  
  // z,sin(phi),tan(theta),p,betaGamma,ncls
  fDeDxQmax  = new THnSparseS("fDeDxQmax","Qmax;(z,sin(phi),tan(theta),p,betaGamma,ncls); TPC signal Qmax (a.u.)",7,binsQA,xminQA,xmaxQA);
  fDeDxQtot  = new THnSparseS("fDeDxQtot","Qmax;(z,sin(phi),tan(theta),p,betaGamma,ncls); TPC signal Qmax (a.u.)",7,binsQA,xminQA,xmaxQA);
  fDeDxRatio = new THnSparseS("fDeDxRatio","Qmax;(z,sin(phi),tan(theta),p,betaGamma,ncls); TPC signal Qmax (a.u.)",7,binsQA,xminRa,xmaxRa);
  fDeDxShortMediumRatio = new THnSparseS("fDeDxQmax","Qmax;(z,sin(phi),tan(theta),p,betaGamma,ncls); TPC signal Qmax (a.u.)",7,binsQA,xminRa,xmaxRa);
  fDeDxLongMediumRatio  = new THnSparseS("fDeDxQmax","Qmax;(z,sin(phi),tan(theta),p,betaGamma,ncls); TPC signal Qmax (a.u.)",7,binsQA,xminRa,xmaxRa);
  BinLogX(fDeDxQmax,4); BinLogX(fDeDxQmax,5);
  BinLogX(fDeDxQtot,4); BinLogX(fDeDxQtot,5);
  BinLogX(fDeDxRatio,4); BinLogX(fDeDxRatio,5);
  BinLogX(fDeDxShortMediumRatio,4); BinLogX(fDeDxShortMediumRatio,5);
  BinLogX(fDeDxLongMediumRatio,4); BinLogX(fDeDxLongMediumRatio,5);
  //
  fHistNTracks = new TH1F("ntracks","Number of Tracks per Event; number of tracks per event; number of tracks",1001,-0.5,1000.5);
  fClusters = new TH1F("signal","Number of Clusters per track; number of clusters per track n_{cl}; counts",40,0,160);
  fPileUp = new TH2F("PileUp","timing plots.; #Delta L ; dEdx signal ",400,-1,1,400,0,200);
  fLandau = new TH2F("Landau","Landau.; pad type ; cluster charge ",3,-0.5,2.5,400,0,1000);
  //
  AliInfo("Non Default Constructor");  
  //
}


AliTPCcalibPID::~AliTPCcalibPID(){
  //
  //
  //
}



void AliTPCcalibPID::Process(AliESDEvent *event) {
  //
  //
  //
  if (!event) {
    Printf("ERROR: ESD not available");
    return;
  }  

  Int_t ntracks=event->GetNumberOfTracks(); 
  fHistNTracks->Fill(ntracks);
  
  AliESDfriend *ESDfriend=static_cast<AliESDfriend*>(event->FindListObject("AliESDfriend"));
  if (!ESDfriend) {
   Printf("ERROR: ESDfriend not available");
   return;
  }  
  //
  // track loop
  //
  for (Int_t i=0;i<ntracks;++i) {

    AliESDtrack *track = event->GetTrack(i);
    if (!track) continue;
    
    
    const AliExternalTrackParam * trackIn = track->GetInnerParam();
    const AliExternalTrackParam * trackOut = track->GetOuterParam();
    if (!trackIn) continue;
    if (!trackOut) continue;

    // calculate necessary track parameters
    //Double_t meanP = 0.5*(trackIn->GetP() + trackOut->GetP());
    Double_t meanP = trackIn->GetP();
    Double_t d = trackIn->GetLinearD(0,0);
    Int_t NclsDeDx = track->GetTPCNcls();

    //if (meanP > 0.7 || meanP < 0.2) continue;
     if (NclsDeDx < 60) continue;     

    // exclude tracks which do not look like primaries or are simply too short or on wrong sectors

    //if (TMath::Abs(trackIn->GetSnp()) > 3*0.4) continue;
    //if (TMath::Abs(trackIn->GetZ()) > 150) continue;   
    //if (seed->CookShape(1) > 1) continue;
    //if (TMath::Abs(trackIn->GetY()) > 20) continue;
    //if (TMath::Abs(d)>20) continue;   // distance to the 0,0; select only tracks which cross chambers under proper angle
    //if (TMath::Abs(trackIn->GetSnp()) > 0.6) continue;
    //if (TMath::Abs(trackOut->GetSnp()) > 0.2) continue;
    if (TMath::Abs(trackIn->GetAlpha()+0.872665)<0.01 || TMath::Abs(trackOut->GetAlpha()+0.872665)<0.01) continue;  // Funny sector !
    
    // Get seeds
    AliESDfriendTrack *friendTrack = ESDfriend->GetTrack(i);
    if (!friendTrack) continue;
    TObject *calibObject;
    AliTPCseed *seed = 0;
    for (Int_t l=0;(calibObject=friendTrack->GetCalibObject(l));++l) {
      if ((seed=dynamic_cast<AliTPCseed*>(calibObject))) break;
    }    

    if (seed) {
      if (meanP > 30 && TMath::Abs(trackIn->GetSnp()) < 0.2) fClusters->Fill(track->GetTPCNcls());
      // calculate dEdx
      // (Double_t low, Double_t up, Int_t type, Int_t i1, Int_t i2, Bool_t shapeNorm,Bool_t posNorm, Int_t padNorm, Int_t returnVal)
      Double_t TPCsignalTot    = (1/fMIP)*seed->CookdEdxNorm(fLowerTrunc,fUpperTrunc,0,0,159,fUseShapeNorm,fUsePosNorm,fUsePadNorm,0);
      Double_t TPCsignalMax    = (1/fMIP)*seed->CookdEdxNorm(fLowerTrunc,fUpperTrunc,1,0,159,fUseShapeNorm,fUsePosNorm,fUsePadNorm,0);
      //
      Double_t TPCsignalShort  = seed->CookdEdxNorm(fLowerTrunc,fUpperTrunc,1,0,63,fUseShapeNorm,fUsePosNorm,fUsePadNorm,0);
      Double_t TPCsignalMedium = seed->CookdEdxNorm(fLowerTrunc,fUpperTrunc,1,63,63+64,fUseShapeNorm,fUsePosNorm,fUsePadNorm,0);
      Double_t TPCsignalLong   = seed->CookdEdxNorm(fLowerTrunc,fUpperTrunc,1,63+64,159,fUseShapeNorm,fUsePosNorm,fUsePadNorm,0);
      //
      //Double_t driftMismatch = seed->GetDriftTimeMismatch(0,0.95,0.5);
      Double_t driftMismatch = 0;
      Double_t drift = 1 - (TMath::Abs(trackIn->GetZ()) + TMath::Abs(trackOut->GetZ()))/500.;        
      
      // particle identification
      Double_t mass = 0.105658369;// muon
      
      if (meanP > 30) {
	fPileUp->Fill(driftMismatch,TPCsignalMax);
	//
	fLandau->Fill(0.1,0.6);
      }
      //var. of interest, z,sin(phi),tan(theta),p,betaGamma,ncls
      Double_t snp = TMath::Abs(trackIn->GetSnp());
      Double_t tgl = TMath::Abs(trackIn->GetTgl());
      //
      Double_t vecQmax[7] = {TPCsignalMax,drift,snp,tgl,meanP,meanP/mass,NclsDeDx};
      Double_t vecQtot[7] = {TPCsignalTot,drift,snp,tgl,meanP,meanP/mass,NclsDeDx};
      Double_t vecRatio[7] = {TPCsignalMax/TPCsignalTot,drift,snp,tgl,meanP,meanP/mass,NclsDeDx};
      Double_t vecShortMediumRatio[7] = {TPCsignalShort/TPCsignalMedium,drift,snp,tgl,meanP,meanP/mass,NclsDeDx};
      Double_t vecLongMediumRatio[7] = {TPCsignalLong/TPCsignalMedium,drift,snp,tgl,meanP,meanP/mass,NclsDeDx};
          
      fDeDxQmax->Fill(vecQmax); 
      fDeDxQtot->Fill(vecQtot); 
      fDeDxRatio->Fill(vecRatio); 
      fDeDxShortMediumRatio->Fill(vecShortMediumRatio); 
      fDeDxLongMediumRatio->Fill(vecLongMediumRatio);       
      
    }
    
  }
  
  
}



void AliTPCcalibPID::Analyze() {


  return;

}


Long64_t AliTPCcalibPID::Merge(TCollection *li) {

  TIterator* iter = li->MakeIterator();
  AliTPCcalibPID* cal = 0;

  while ((cal = (AliTPCcalibPID*)iter->Next())) {
    if (!cal->InheritsFrom(AliTPCcalibPID::Class())) {
      Error("Merge","Attempt to add object of class %s to a %s", cal->ClassName(), this->ClassName());
      return -1;
    }
    
    if (cal->GetHistNTracks()) fHistNTracks->Add(cal->GetHistNTracks());
    if (cal->GetHistClusters()) fClusters->Add(cal->GetHistClusters());
    if (cal->GetHistPileUp()) fPileUp->Add(cal->GetHistPileUp());
    if (cal->GetHistLandau()) fLandau->Add(cal->GetHistLandau());
    //
    if (cal->GetHistQmax()) fDeDxQmax->Add(cal->GetHistQmax());
    if (cal->GetHistQtot()) fDeDxQtot->Add(cal->GetHistQtot());
    if (cal->GetHistRatio()) fDeDxRatio->Add(cal->GetHistRatio());
    if (cal->GetHistShortMediumRatio()) fDeDxShortMediumRatio->Add(cal->GetHistShortMediumRatio());
    if (cal->GetHistLongMediumRatio()) fDeDxLongMediumRatio->Add(cal->GetHistLongMediumRatio());
  }
  
  return 0;
  
}




void AliTPCcalibPID::MakeReport() {
  //
  // 1. standard dEdx picture
  TCanvas *cDeDx = new TCanvas("cDeDx","cDeDx");
  GetHistQmax()->Projection(0,4)->Draw("colz");
  return;
}



void AliTPCcalibPID::BinLogX(THnSparse *h, Int_t axisDim) {

  // Method for the correct logarithmic binning of histograms

  TAxis *axis = h->GetAxis(axisDim);
  int bins = axis->GetNbins();

  Double_t from = axis->GetXmin();
  Double_t to = axis->GetXmax();
  Double_t *new_bins = new Double_t[bins + 1];
   
  new_bins[0] = from;
  Double_t factor = pow(to/from, 1./bins);
  
  for (int i = 1; i <= bins; i++) {
   new_bins[i] = factor * new_bins[i-1];
  }
  axis->Set(bins, new_bins);
  delete new_bins;
  
}

