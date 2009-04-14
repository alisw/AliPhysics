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

cal->GetHistRatioMaxTot()->Projection(0)->Fit("gaus")
cal->GetHistRatioTracklet()->Projection(0)->Fit("gaus")

1.b) Update OCDB:
.x $ALICE_ROOT/TPC/macros/ConfigOCDB.C
AliTPCClusterParam * parcl = AliTPCcalibDB::Instance()->GetClusterParam();
*parcl->fQpadMnorm)[ipad] = oldvalue*corrFactor

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
   fDeDxRatioMaxTot(0),
   fDeDxRatioQmax(0),
   fDeDxRatioQtot(0),
   fDeDxRatioTruncQtot(0),
   fDeDxRatioTruncQmax(0)
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
   fDeDxRatioMaxTot(0),
   fDeDxRatioQmax(0), 
   fDeDxRatioQtot(0) ,
   fDeDxRatioTruncQtot(0),
   fDeDxRatioTruncQmax(0)
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
  fUsePosNorm = 0;
  fUsePadNorm = 1;
  //
  fIsCosmic  = kTRUE;
  //
  //                  dE/dx,  z, phi, theta,    p,  bg, ncls, tracklet type
  Int_t binsQA[8]    = {300, 20,  10,    20,   50, 50,  8,    5};
  Double_t xminQA[8] = {0.2,  -1,    0,  -1.5, 0.01, 0.1, 60,   0};
  Double_t xmaxQA[8] = {10.,  1,  0.9,   1.5,   50, 500, 160,  5};
  //
  //
  //
  //                  dE/dx,  z, phi, theta, dEdx, dEdx*dl, ncls, tracklet type
  Int_t    binsRA[9] = {100, 20,  10,    20,   25,  25,      4,    5};
  Double_t xminRa[9] = {0.5, -1,    0,  -1.5,  0.2, 0.2,     60,    0};
  Double_t xmaxRa[9] = {2,  1,  0.9,   1.5,    5,   5,    160,    5};
  
  // z,sin(phi),tan(theta),p,betaGamma,ncls
  fDeDxQmax  = new THnSparseS("fDeDxQmax","Qmax;(z,sin(phi),tan(theta),p,betaGamma,ncls,type); TPC signal Qmax (a.u.)",8,binsQA,xminQA,xmaxQA);
  fDeDxQtot  = new THnSparseS("fDeDxQtot","Q_{tot};(z,sin(phi),tan(theta),p,betaGamma,ncls,type); TPC signal Qmax (a.u.)",8,binsQA,xminQA,xmaxQA);
  //
  // ratio histograms
  //
  fDeDxRatioMaxTot = new THnSparseS("fDeDxRatioMaxTot","Q_{max}/Q_{tot};(z,sin(phi),tan(theta),dEdx,dEdx*dl,ncls,type); TPC signal Qmax/Qtot (a.u.)",8,binsRA,xminRa,xmaxRa);
  fDeDxRatioQmax = new THnSparseS("fDeDxRatioQmax","Q_{mtracklet}/Q_{mtrack};(z,sin(phi),tan(theta),dEdx,dEdx*dl,ncls,type,qtupe); TPC signal Tracklet/Track (a.u.)",8,binsRA,xminRa,xmaxRa);
  fDeDxRatioQtot = new THnSparseS("fDeDxRatioQtot","Q_{ttracklet}/Q_{ttrack};(z,sin(phi),tan(theta),dEdx,dEdx*dl,ncls,type,qtupe); TPC signal Tracklet/Track (a.u.)",8,binsRA,xminRa,xmaxRa);
  fDeDxRatioTruncQmax = new THnSparseS("fDeDxRatioTrunQmax","Q_{max}/Q_{maxtrunc};(z,sin(phi),tan(theta),dEdx,dEdx*dl,ncls,type,qtupe); TPC signal Full/Trunc. (a.u.)",8,binsRA,xminRa,xmaxRa);


  fDeDxRatioTruncQtot = new THnSparseS("fDeDxRatioTruncQtot","Q_{tot}/Q_{tottrunc};(z,sin(phi),tan(theta),dEdx,dEdx*dl,ncls,type,qtupe); TPC signal Full/Trunc (a.u.)",8,binsRA,xminRa,xmaxRa);


  BinLogX(fDeDxQmax,4); BinLogX(fDeDxQmax,5);BinLogX(fDeDxQmax,0);
  BinLogX(fDeDxQtot,4); BinLogX(fDeDxQtot,5);BinLogX(fDeDxQmax,0);
  //
  BinLogX(fDeDxRatioMaxTot,4); BinLogX(fDeDxRatioMaxTot,5);
  BinLogX(fDeDxRatioQmax,4); BinLogX(fDeDxRatioQmax,5);
  BinLogX(fDeDxRatioQtot,4); BinLogX(fDeDxRatioQtot,5);
  BinLogX(fDeDxRatioTruncQmax,4); BinLogX(fDeDxRatioTruncQmax,5);
  BinLogX(fDeDxRatioTruncQtot,4); BinLogX(fDeDxRatioTruncQtot,5);
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
  delete fHistNTracks;            //  histogram showing number of ESD tracks per event
  delete fClusters;               //  histogram showing the number of clusters per track
  delete fPileUp;                 //  histogram which shows correlation between time mismatch and dEdx signal
  delete fLandau;                 //  histogran which shows Landau distribution for the three pad geometries
  //
  delete fDeDxQmax;               //  histogram which shows dEdx (Qmax) as a function of z,sin(phi),tan(theta),p,betaGamma
  delete fDeDxQtot;               //  histogram which shows dEdx (Qtot) as a function of z,sin(phi),tan(theta),p,betaGamma
  //
  // ratio histograms
  //
  delete fDeDxRatioMaxTot;              //  histogram which shows dEdx ratio (Qmax/Qtot) as a function of z,sin(phi),tan(theta),dEdx,dEdx*dl
  delete fDeDxRatioQmax;   // dEdx ratio (tracklet/track) as a function of z,sin(phi),tan(theta),dEdx,dEdx*dl
  delete fDeDxRatioQtot;   // dEdx ratio (tracklet/track) as a function of z,sin(phi),tan(theta),dEdx,dEdx*dl
  delete fDeDxRatioTruncQtot;   // dEdx ratio (tracklet/track) as a function of z,sin(phi),tan(theta),dEdx,dEdx*dl
  delete fDeDxRatioTruncQmax;   // dEdx ratio (tracklet/track) as a function of z,sin(phi),tan(theta),dEdx,dEdx*d
}



void AliTPCcalibPID::Process(AliESDEvent *event) {
  //
  //
  //
  if (!event) {
    Printf("ERROR: ESD not available");
    return;
  }  
  const Int_t kMinClustersTracklet = 25;
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
    
    
    AliExternalTrackParam * trackIn  = (AliExternalTrackParam *)track->GetInnerParam();
    AliExternalTrackParam * trackOut = (AliExternalTrackParam *)track->GetOuterParam();
    if (!trackIn) continue;
    if (!trackOut) continue;

    // calculate necessary track parameters
    //Double_t meanP = 0.5*(trackIn->GetP() + trackOut->GetP());
    Double_t meanP = trackIn->GetP();
    //Double_t d = trackIn->GetLinearD(0,0);
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
      Double_t TPCsignalMax    = (1/fMIP)*seed->CookdEdxNorm(fLowerTrunc,fUpperTrunc,1,0,159,fUseShapeNorm,fUsePosNorm,fUsePadNorm,0);
      //
      //
      //Double_t driftMismatch = seed->GetDriftTimeMismatch(0,0.95,0.5);
      Double_t driftMismatch = 0;
      //      Double_t drift = 1 - (TMath::Abs(trackIn->GetZ()) + TMath::Abs(trackOut->GetZ()))/500.;        
      
      // particle identification
      Double_t mass = 0.105658369;// muon
      
      if (meanP > 30) {
	fPileUp->Fill(driftMismatch,TPCsignalMax);
	//
	fLandau->Fill(0.1,0.6);
      }
      //var. of interest, z,sin(phi),tan(theta),p,betaGamma,ncls
      Double_t snpIn   = TMath::Abs(trackIn->GetSnp());
      Double_t snpOut  = TMath::Abs(trackOut->GetSnp());
      Double_t tglIn   = trackIn->GetTgl();
      Double_t tglOut  = trackOut->GetTgl();
      Double_t driftIn = trackIn->GetZ()/250.;
      Double_t driftOut= trackIn->GetZ()/250.;
      //
      // dEdx in region
      //
      Float_t dEdxTot[5],dEdxTotFull[5];
      Float_t dEdxMax[5],dEdxMaxFull[5];
      Float_t   ncl[5];
      for (Int_t itype=0; itype<5;itype++){
	Int_t row0=0, row1 =159;
	if (itype==1) {row0=0;      row1 = 63;};
	if (itype==2) {row0= 64;    row1=63+64;}
	if (itype==3) {row0= 63+64; row1=159;}
	if (fUsePosNorm==0){
	  dEdxTot[itype]= (1/fMIP)*seed->CookdEdxNorm(fLowerTrunc,fUpperTrunc,0,row0,row1,fUseShapeNorm,fUsePosNorm,fUsePadNorm,0);
	  dEdxMax[itype]= (1/fMIP)*seed->CookdEdxNorm(fLowerTrunc,fUpperTrunc,1,row0,row1,fUseShapeNorm,fUsePosNorm,fUsePadNorm,0);
	  // non trucated dEdx
	  dEdxTotFull[itype]= (1/fMIP)*seed->CookdEdxNorm(0.0,0.99,0,row0,row1,fUseShapeNorm,fUsePosNorm,fUsePadNorm,0);
	  dEdxMaxFull[itype]= (1/fMIP)*seed->CookdEdxNorm(0.0,0.99,1,row0,row1,fUseShapeNorm,fUsePosNorm,fUsePadNorm,0);
	}else{
	  dEdxTot[itype]= (1/fMIP)*seed->CookdEdxAnalytical(fLowerTrunc,fUpperTrunc,0,row0,row1);
	  dEdxMax[itype]= (1/fMIP)*seed->CookdEdxAnalytical(fLowerTrunc,fUpperTrunc,1,row0,row1);
	  // non trucated dEdx
	  dEdxTotFull[itype]= (1/fMIP)*seed->CookdEdxAnalytical(0.0,0.99,0,row0,row1);
	  dEdxMaxFull[itype]= (1/fMIP)*seed->CookdEdxAnalytical(0.0,0.99,1,row0,row1);
	}

	ncl[itype]=(seed->CookdEdxAnalytical(fLowerTrunc,fUpperTrunc,1,row0,row1,2));
      }
      //
      //
      //
      Float_t wmeanTot=0, wmeanMax=0, sumW=0;
      Double_t length[3] = {0.75,1,1.5};
      //       //
      
      for (Int_t ipad=0; ipad<3; ipad++){
	if (ncl[1+ipad]<3) continue;
	Double_t weight = Double_t(ncl[1+ipad])*TMath::Sqrt(length[ipad]);
	wmeanTot+=dEdxTot[1+ipad]*weight;
	wmeanMax+=dEdxMax[1+ipad]*weight;
	sumW+=weight;
      }
      if (sumW>0){
	dEdxTot[4]= wmeanTot/sumW;
	dEdxMax[4]= wmeanMax/sumW;	
      }
      for (Int_t itype=0;itype<5;itype++){
	//
	Float_t snp=(TMath::Abs(snpIn)+TMath::Abs(snpOut))*0.5;
	Float_t tgl=(tglIn+tglOut)*0.5;
	Float_t drift = (driftIn+driftOut)*0.5;
	if (itype==1) {snp = TMath::Abs(snpIn); tgl = tglIn; drift= driftIn;};
	if (itype==3) {snp = TMath::Abs(snpOut); tgl = tglOut;drift=driftOut;};
	if (ncl[itype]<kMinClustersTracklet) continue;
	Float_t deltaL = TMath::Sqrt(1+snp*snp+tgl*tgl);
	//
	Double_t vecQmax[8] = {dEdxMax[itype],drift,snp,tgl,meanP,meanP/mass,NclsDeDx, itype};
	Double_t vecQtot[8] = {dEdxTot[itype],drift,snp,tgl,meanP,meanP/mass,NclsDeDx, itype};
	//
	//
	//
	Double_t ratioMaxTot           = (dEdxTot[itype]>0)  ? dEdxMax[itype]/dEdxTot[itype]:0;
	Double_t ratioTrackletTot      = (dEdxTot[0]>0)      ? dEdxTot[itype]/dEdxTot[0]:0;
	Double_t ratioTrackletMax      = (dEdxMax[0]>0)      ? dEdxMax[itype]/dEdxMax[0]:0;
	Double_t ratioTrackletTruncTot = (dEdxTot[0]>0)      ? dEdxTotFull[itype]/dEdxTot[itype]:0;
	Double_t ratioTrackletTruncMax = (dEdxMax[0]>0)      ? dEdxMaxFull[itype]/dEdxMax[itype]:0;

	Double_t vecRatioMaxTot[8]      = {ratioMaxTot,      drift,snp,tgl,dEdxTot[0],  dEdxTot[0]*deltaL,NclsDeDx,itype};
	Double_t vecRatioTrackletTot[8] = {ratioTrackletTot, drift,snp,tgl,dEdxTot[0],  dEdxTot[0]*deltaL,NclsDeDx,itype};	
	Double_t vecRatioTrackletMax[8] = {ratioTrackletMax, drift,snp,tgl,dEdxMax[0],  dEdxMax[0]*deltaL,NclsDeDx,itype};	
	Double_t vecRatioTrackletTruncTot[8] = {ratioTrackletTruncTot, drift,snp,tgl,dEdxTot[0],  dEdxTot[0]*deltaL,NclsDeDx,itype};	
	Double_t vecRatioTrackletTruncMax[8] = {ratioTrackletTruncMax, drift,snp,tgl,dEdxMax[0],  dEdxMax[0]*deltaL,NclsDeDx,itype};	
	fDeDxQmax->Fill(vecQmax); 
	fDeDxQtot->Fill(vecQtot); 
	fDeDxRatioMaxTot->Fill(vecRatioMaxTot); 
	fDeDxRatioQmax->Fill(vecRatioTrackletTot); 
	fDeDxRatioQtot->Fill(vecRatioTrackletMax); 
	fDeDxRatioTruncQmax->Fill(vecRatioTrackletTruncTot); 
	fDeDxRatioTruncQtot->Fill(vecRatioTrackletTruncMax); 
	//
	TTreeSRedirector * cstream =  GetDebugStreamer();
	if (cstream){
	  TVectorD vQT(9,vecQtot);
	  TVectorD vQM(9,vecQmax);
	  TVectorD vQMT(9, vecRatioMaxTot);
	  TVectorD vQRT(9,vecRatioTrackletTot);
	  TVectorD vQRM(9,vecRatioTrackletMax);
	  TVectorD vQRTT(9,vecRatioTrackletTruncTot);
	  TVectorD vQRTM(9,vecRatioTrackletTruncMax);
	  TVectorF  vNcl(5,ncl);
	  (*cstream) << "dEdx" <<
	    "run="<<fRun<<              //  run number
	    "event="<<fEvent<<          //  event number
	    "time="<<fTime<<            //  time stamp of event
	    "trigger="<<fTrigger<<      //  trigger
	    "mag="<<fMagF<<             //  magnetic field	      
	    "track.="<<seed<<           //  original seed
	    "trin.="<<trackIn<<           //  inner param
	    "trout.="<<trackOut<<         //  outer param
	    "ncl.="<<&vNcl<<            //  number of clusters
	    "itype="<<itype<<
	    //
	    "vQT.="<<&vQT<<             // trunc mean - total charge
	    "vQM.="<<&vQM<<	        // trunc mean - max charge 
	    //
	    "vQMT.="<<&vQMT<<           // ratio max/tot      
	    "vQRT.="<<&vQRT<<           // ratio tracklet/track - total charge
	    "vQRM.="<<&vQRM<<           // ratio tracklet/track - max charge
	    "vQRTT.="<<&vQRTT<<         // ratio trunc/full     - total charge
	    "vQRTM.="<<&vQRTM<<         // ratio trunc/full     - max charge
	    "\n";
	}
      }   
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
    if (cal->GetHistRatioMaxTot()) fDeDxRatioMaxTot->Add(cal->GetHistRatioMaxTot());
    if (cal->GetHistRatioQmax()) fDeDxRatioQmax->Add(cal->GetHistRatioQmax());
    if (cal->GetHistRatioQtot()) fDeDxRatioQtot->Add(cal->GetHistRatioQtot());
    if (cal->GetHistRatioTruncQmax()) fDeDxRatioTruncQmax->Add(cal->GetHistRatioTruncQmax());
    if (cal->GetHistRatioTruncQtot()) fDeDxRatioTruncQtot->Add(cal->GetHistRatioTruncQtot());
  }
  
  return 0;
  
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





void AliTPCcalibPID::MakeReport(const char *outputpath) {
  //
  // Make a standard QA plots
  //
  for (Int_t ipad=0;ipad<4;ipad++){
    DrawRatioTot(ipad,outputpath);
    DrawRatioMax(ipad,outputpath);
  }
  DrawRatiodEdx(0.5,3,outputpath);
  DrawResolBGQtot(140,160,1,40,outputpath);
  DrawResolBGQmax(140,160,1,40,outputpath);
  return;
}

void  AliTPCcalibPID::DrawRatioTot(Int_t ipad, const char* outputpath){
  //
  // Draw default ratio histogram for given pad type
  // ipad - 0 - short pads
  //        1 - medium pads
  //        2 - long pads
  //
  //  Int_t kmicolors[10]={1,2,3,4,6,7,8,9,10,11};
  Int_t kmimarkers[10]={21,22,23,24,25,26,27,28,29,30};
  AliTPCcalibPID * pid = this;
  TCanvas *canvas= new TCanvas(Form("QtotRatio%d)",ipad),Form("QtotRatio%d)",ipad),600,800);
  canvas->Divide(3,2);
  pid->GetHistRatioQtot()->GetAxis(7)->SetRange(ipad+2,ipad+2);
  canvas->cd(1);
  TH1 * his0 = pid->GetHistRatioQtot()->Projection(0);
  his0->SetXTitle("dEdx_{tracklet}/dEdx_{track}");
  his0->SetYTitle("");  
  his0->SetTitle(Form("dEdx_{tracklet}/dEdx_{track} type %d",ipad));
  his0->Draw();
  //
  canvas->cd(2);
  TH1 * hprofz = (TH1*) (pid->GetHistRatioQtot()->Projection(0,1)->ProfileX());
  hprofz->SetXTitle("drift length");
  hprofz->SetYTitle("dEdx_{tracklet}/dEdx_{track}");
  hprofz->SetTitle(Form("dEdx_{tracklet}/dEdx_{track} type %d",ipad));
  hprofz->SetMarkerStyle(kmimarkers[0]);
  hprofz->SetMaximum(1.1);
  hprofz->SetMinimum(0.9);
  hprofz->Draw();
  //
  canvas->cd(3);
  TH1 * hprofphi = (TH1*) (pid->GetHistRatioQtot()->Projection(0,2)->ProfileX());
  hprofphi->SetXTitle("sin(#phi)");
  hprofphi->SetYTitle("dEdx_{tracklet}/dEdx_{track}");
  hprofphi->SetTitle(Form("dEdx_{tracklet}/dEdx_{track} type %d",ipad));    
  hprofphi->SetMarkerStyle(kmimarkers[0]);  
  hprofz->SetTitle(Form("dEdx_{tracklet}/dEdx_{track} type %d",ipad));
  hprofphi->SetMaximum(1.1);
  hprofphi->SetMinimum(0.9);
  hprofphi->Draw();
  //
  canvas->cd(4);
  TH1 * hproftheta = (TH1*) (pid->GetHistRatioQtot()->Projection(0,3)->ProfileX());
  hproftheta->SetXTitle("tan(#theta)");
  hproftheta->SetYTitle("dEdx_{tracklet}/dEdx_{track}");
  hproftheta->SetTitle(Form("dEdx_{tracklet}/dEdx_{track} type %d",ipad));
  hproftheta->SetMarkerStyle(kmimarkers[0]);
  hproftheta->SetMaximum(1.1);
  hproftheta->SetMinimum(0.9);
  hproftheta->Draw();
  //
  canvas->cd(5);
  TH1 * hprofdedx = (TH1*) (pid->GetHistRatioQtot()->Projection(0,4)->ProfileX());
  hprofdedx->SetXTitle("dEdx_{track}");
  hprofdedx->SetYTitle("dEdx_{tracklet}/dEdx_{track}");
  hprofdedx->SetTitle(Form("dEdx_{tracklet}/dEdx_{track} type %d",ipad));
  hprofdedx->SetMarkerStyle(kmimarkers[0]);
  hprofdedx->SetMaximum(1.1);
  hprofdedx->SetMinimum(0.9);
  hprofdedx->Draw();

  canvas->cd(6);
  TH1 * hprofdedxL = (TH1*) (pid->GetHistRatioQtot()->Projection(0,5)->ProfileX());
  hprofdedxL->SetXTitle("dEdx_{track}#Delta_{x}");
  hprofdedxL->SetYTitle("dEdx_{tracklet}/dEdx_{track}");
  hprofdedxL->SetTitle(Form("dEdx_{tracklet}/dEdx_{track} type %d",ipad));
  hprofdedxL->SetMarkerStyle(kmimarkers[0]);
  hprofdedxL->SetMaximum(1.1);
  hprofdedxL->SetMinimum(0.9);
  hprofdedxL->Draw();



  canvas->SaveAs(Form("%s/QtotRatioType%d.eps",outputpath,ipad));
  canvas->SaveAs(Form("%s/QtotRatioType%d.gif",outputpath,ipad));
}

void  AliTPCcalibPID::DrawRatioMax(Int_t ipad, const char* outputpath){
  //
  // Draw default ration histogram for given pad type
  // ipad - 0 - short pads
  //        1 - medium pads
  //        2 - long pads
  //
  //  Int_t kmicolors[10]={1,2,3,4,6,7,8,9,10,11};
  Int_t kmimarkers[10]={21,22,23,24,25,26,27,28,29,30};
  AliTPCcalibPID * pid = this;
  TCanvas *canvas= new TCanvas(Form("QmaxRatio%d)",ipad),Form("QmaxRatio%d)",ipad),600,800);
  canvas->Divide(3,2);
  pid->GetHistRatioQmax()->GetAxis(7)->SetRange(ipad+2,ipad+2);
  canvas->cd(1);
  TH1 * his0 = pid->GetHistRatioQmax()->Projection(0);
  his0->SetXTitle("dEdx_{tracklet}/dEdx_{track}");
  his0->SetYTitle("");  
  his0->SetTitle(Form("dEdx_{tracklet}/dEdx_{track} type %d",ipad));
  his0->Draw();
  //
  canvas->cd(2);
  TH1 * hprofz = (TH1*) (pid->GetHistRatioQmax()->Projection(0,1)->ProfileX());
  hprofz->SetXTitle("drift length");
  hprofz->SetYTitle("dEdx_{tracklet}/dEdx_{track}");
  hprofz->SetTitle(Form("dEdx_{tracklet}/dEdx_{track} type %d",ipad));
  hprofz->SetMarkerStyle(kmimarkers[0]);
  hprofz->SetMaximum(1.1);
  hprofz->SetMinimum(0.9);
  hprofz->Draw();
  //
  canvas->cd(3);
  TH1 * hprofphi = (TH1*) (pid->GetHistRatioQmax()->Projection(0,2)->ProfileX());
  hprofphi->SetXTitle("sin(#phi)");
  hprofphi->SetYTitle("dEdx_{tracklet}/dEdx_{track}");
  hprofphi->SetTitle(Form("dEdx_{tracklet}/dEdx_{track} type %d",ipad));    
  hprofphi->SetMarkerStyle(kmimarkers[0]);  
  hprofphi->SetTitle(Form("dEdx_{tracklet}/dEdx_{track} type %d",ipad));
  hprofphi->SetMaximum(1.1);
  hprofphi->SetMinimum(0.9);
  hprofphi->Draw();
  //
  canvas->cd(4);
  TH1 * hproftheta = (TH1*) (pid->GetHistRatioQmax()->Projection(0,3)->ProfileX());
  hproftheta->SetXTitle("tan(#theta)");
  hproftheta->SetYTitle("dEdx_{tracklet}/dEdx_{track}");
  hproftheta->SetTitle(Form("dEdx_{tracklet}/dEdx_{track} type %d",ipad));
  hproftheta->SetMarkerStyle(kmimarkers[0]);
  hproftheta->SetMaximum(1.1);
  hproftheta->SetMinimum(0.9);
  hproftheta->Draw();

  canvas->cd(5);
  TH1 * hprofdedx = (TH1*) (pid->GetHistRatioQmax()->Projection(0,4)->ProfileX());
  hprofdedx->SetXTitle("dEdx_{track}");
  hprofdedx->SetYTitle("dEdx_{tracklet}/dEdx_{track}");
  hprofdedx->SetTitle(Form("dEdx_{tracklet}/dEdx_{track} type %d",ipad));
  hprofdedx->SetMarkerStyle(kmimarkers[0]);
  hprofdedx->SetMaximum(1.1);
  hprofdedx->SetMinimum(0.9);
  hprofdedx->Draw();

  canvas->cd(6);
  TH1 * hprofdedxL = (TH1*) (pid->GetHistRatioQmax()->Projection(0,5)->ProfileX());
  hprofdedxL->SetXTitle("dEdx_{track}#Delta_{x}");
  hprofdedxL->SetYTitle("dEdx_{tracklet}/dEdx_{track}");
  hprofdedxL->SetTitle(Form("dEdx_{tracklet}/dEdx_{track} type %d",ipad));
  hprofdedxL->SetMarkerStyle(kmimarkers[0]);
  hprofdedxL->SetMaximum(1.1);
  hprofdedxL->SetMinimum(0.9);
  hprofdedxL->Draw();


  canvas->SaveAs(Form("%s/QmaxRatioType%d.eps",outputpath,ipad));
  canvas->SaveAs(Form("%s/QmaxRatioType%d.gif",outputpath,ipad));
}



void AliTPCcalibPID::DrawRatiodEdx(Float_t demin, Float_t demax, const char* outputpath){
  //
  //
  //
  //
  //  Int_t kmicolors[10]={1,2,3,4,6,7,8,9,10,11};
  Int_t kmimarkers[10]={21,22,23,24,25,26,27,28,29,30};
  AliTPCcalibPID * pid = this;
  TCanvas *canvas= new TCanvas("QRatiodEdx","QRatiodEdx",600,800);
  canvas->Divide(2,4);
  canvas->SetLogx(kTRUE);
  TH1 * hprofP=0;
  for (Int_t ipad=0;ipad<4;ipad++){
    pid->GetHistRatioQmax()->GetAxis(7)->SetRange(ipad+2,ipad+2);
    pid->GetHistRatioQtot()->GetAxis(7)->SetRange(ipad+2,ipad+2);
    pid->GetHistRatioQmax()->GetAxis(5)->SetRangeUser(demin,demax);
    pid->GetHistRatioQtot()->GetAxis(5)->SetRangeUser(demin,demax);

    canvas->cd(ipad*2+1)->SetLogx(kFALSE);
    hprofP  = (TH1*) (pid->GetHistRatioQmax()->Projection(0,5)->ProfileX());
    hprofP->SetXTitle("dEdx_{track}");
    hprofP->SetYTitle("dEdx_{tracklet}/dEdx_{track}");
    hprofP->SetTitle(Form("Q_{max} dEdx_{tracklet}/dEdx_{track} type %d",ipad));
    hprofP->SetMarkerStyle(kmimarkers[0]);     
    hprofP->SetMaximum(1.1);
    hprofP->SetMinimum(0.9);
    hprofP->Draw();
    // pad  Tot
    canvas->cd(ipad*2+2)->SetLogx(kFALSE);
    hprofP  = (TH1*) (pid->GetHistRatioQtot()->Projection(0,5)->ProfileX());
    hprofP->SetXTitle("dEdx_{track}");
    hprofP->SetYTitle("dEdx_{tracklet}/dEdx_{track}");
    hprofP->SetTitle(Form("Q_{tot} dEdx_{tracklet}/dEdx_{track} type %d",ipad));
    hprofP->SetMarkerStyle(kmimarkers[0]);     
    hprofP->SetMaximum(1.1);
    hprofP->SetMinimum(0.9);
    hprofP->Draw();
  }
  //
  //
  canvas->SaveAs(Form("%s/QratiodEdx.eps",outputpath));
  canvas->SaveAs(Form("%s/QratiodEdx.gif",outputpath));
}



void AliTPCcalibPID::DrawResolBGQtot(Int_t minClusters, Int_t maxClusters, Float_t minp, Float_t maxp, const char *outputpath, Bool_t resol){
  //
  // 
  //
  AliTPCcalibPID * pid = this;
  Int_t kmimarkers[10]={21,22,23,24,25,26,27,28,29,30};

  TObjArray arr;
  TH2 * his   =0; 
  TH1 * hmean =0;
  TH1 * hsigma=0;
  //
  // set cut
  pid->GetHistQtot()->GetAxis(6)->SetRangeUser(minClusters,maxClusters);
  pid->GetHistQtot()->GetAxis(5)->SetRangeUser(1,10000);
  pid->GetHistQtot()->GetAxis(4)->SetRangeUser(minp,maxp);
  TCanvas *canvas= new TCanvas("dEdxResolQ_{Tot}","dEdxResolQ_{Tot}",800,600);
  canvas->Divide(2,3);
  //
  //
  //
  for (Int_t ipad=0;ipad<5;ipad++){
    canvas->cd(1+ipad)->SetLogx(kTRUE);
    if (ipad<4) pid->GetHistQtot()->GetAxis(7)->SetRange(ipad+2,ipad+2);
    if (ipad==4) pid->GetHistQtot()->GetAxis(7)->SetRange(1,1);
    his = (TH2*)(pid->GetHistQtot()->Projection(0,5));
    his->FitSlicesY(0,0,-1,0,"QNR",&arr);
    if (resol){
      hmean = (TH1*)arr.At(1);
      hsigma = (TH1*)arr.At(2)->Clone();    
      hsigma->SetMaximum(0.11);
      hsigma->SetMinimum(0.4);    
      hsigma->Divide(hmean);
    }else{
      hsigma = (TH1*)arr.At(1)->Clone();
      hsigma->SetMaximum(2);
      hsigma->SetMinimum(0.5);
    }
    arr.SetOwner(kTRUE);
    arr.Delete();
    delete his;

    hsigma->SetXTitle("#beta#gamma");
    hsigma->SetYTitle("#sigma_{dEdx}/dEdx");
    hsigma->SetTitle(Form("#sigma_{dEdx}/dEdx_{tot} Pad %d",ipad));
    hsigma->SetName(Form("#sigma_{dEdx}/dEdx_{tot} Pad %d",ipad));
    if (ipad==4) {
      hsigma->SetTitle(Form("#sigma_{dEdx}/dEdx_{tot} Full"));
      hsigma->SetName(Form("#sigma_{dEdx}/dEdx_{tot} Full"));
    }
    hsigma->SetMarkerStyle(kmimarkers[0]);
    hsigma->Draw();
  }
  if (resol){
    canvas->SaveAs(Form("%s/dEdxResolTot.eps",outputpath));
    canvas->SaveAs(Form("%s/dEdxResolTot.gif",outputpath));
  }else {
    canvas->SaveAs(Form("%s/dEdxBGTot.eps",outputpath));
    canvas->SaveAs(Form("%s/dEdxBGTot.gif",outputpath));
  }
}

void AliTPCcalibPID::DrawResolBGQmax(Int_t minClusters, Int_t maxClusters, Float_t minp, Float_t maxp, const char *outputpath, Bool_t resol){
  //
  // Int_t minClusters=140;  Int_t maxClusters=200; const char *outputpath="picPID06"
  //
  Int_t kmimarkers[10]={21,22,23,24,25,26,27,28,29,30};
  AliTPCcalibPID * pid = this;
  TObjArray arr;
  TH2 * his   =0; 
  TH1 * hmean =0;
  TH1 * hsigma=0;
  //
  // set cut
  pid->GetHistQmax()->GetAxis(6)->SetRangeUser(minClusters,maxClusters);
  pid->GetHistQmax()->GetAxis(4)->SetRangeUser(minp,maxp);
  pid->GetHistQmax()->GetAxis(5)->SetRangeUser(1,10000);
  TCanvas *canvas= new TCanvas("dEdxResolQ_{Max}","dEdxResolQ_{Max}",800,600);
  canvas->Divide(2,3);
  //
  //
  //
  for (Int_t ipad=0;ipad<5;ipad++){
    canvas->cd(1+ipad)->SetLogx(kTRUE);
    if (ipad<4) pid->GetHistQmax()->GetAxis(7)->SetRange(ipad+2,ipad+2);
    if (ipad==4) pid->GetHistQmax()->GetAxis(7)->SetRange(1,1);
    his = (TH2*)(pid->GetHistQmax()->Projection(0,5));
    his->FitSlicesY(0,0,-1,0,"QNR",&arr);
    if (resol){
      hmean = (TH1*)arr.At(1);
      hsigma = (TH1*)arr.At(2)->Clone();
      hsigma->SetMaximum(0.11);
      hsigma->SetMinimum(0.4);
      hsigma->Divide(hmean);
    }else{
      hsigma = (TH1*)arr.At(1)->Clone();      
      hsigma->SetMaximum(2.);
      hsigma->SetMinimum(0.5);
    }
    arr.SetOwner(kTRUE);
    arr.Delete();
    delete his;
    hsigma->SetXTitle("#beta#gamma");
    hsigma->SetYTitle("#sigma_{dEdx}/dEdx");
    hsigma->SetTitle(Form("#sigma_{dEdx}/dEdx_{max} Pad %d",ipad));
    hsigma->SetName(Form("#sigma_{dEdx}/dEdx_{max} Pad %d",ipad));
    if (ipad==4){
      hsigma->SetTitle(Form("#sigma_{dEdx}/dEdx_{max} Full"));
      hsigma->SetName(Form("#sigma_{dEdx}/dEdx_{max} Full"));
    }
    hsigma->SetMarkerStyle(kmimarkers[0]);
    hsigma->Draw();
  }

  if (resol){
    canvas->SaveAs(Form("%s/dEdxResolMax.eps",outputpath));
    canvas->SaveAs(Form("%s/dEdxResolMax.gif",outputpath));
  }else {
    canvas->SaveAs(Form("%s/dEdxBGMax.eps",outputpath));
    canvas->SaveAs(Form("%s/dEdxBGMax.gif",outputpath));
  }


}



void AliTPCcalibPID::DumpTree(THnSparse * hndim, const char * outname){
  //
  // make a tree
  // the tree will be written to the file - outname
  //
  //
  

  TTreeSRedirector *pcstream = new TTreeSRedirector(outname);
  //
  //
  Double_t position[10];
  Double_t value; 
  Int_t *bins = new Int_t[10];
  //
  //
  const Float_t rmsy0=0.5, rmsz0=1;
  //
  Int_t ndim = hndim->GetNdimensions();
  //
  AliTPCParam * param   = AliTPCcalibDB::Instance()->GetParameters(); 
  for (Long64_t i = 0; i < hndim->GetNbins(); ++i) {
    value = hndim->GetBinContent(i, bins);
    for (Int_t idim = 0; idim < ndim; idim++) {
      position[idim] = hndim->GetAxis(idim)->GetBinCenter(bins[idim]);
    }      
    Int_t sector=36;
    Int_t row=0;
    if (TMath::Abs(position[7]-1.5)<0.1) sector=0;  // inner sector
    if (TMath::Abs(position[7]-3.5)<0.1) row=96;    // long pads

    Double_t ty = TMath::Tan(TMath::ASin(position[2]));
    Double_t tz = position[3]*TMath::Sqrt(1+ty*ty);
    Double_t padLength= param->GetPadPitchLength(sector,row);
    Double_t padWidth = param->GetPadPitchWidth(sector);
    Double_t zwidth   = param->GetZWidth();
    Double_t zlength=250-TMath::Abs(position[1])*250.;
    // diffusion in pad, time bin  units
    Double_t diffT=TMath::Sqrt(zlength)*param->GetDiffT()/padWidth;
    Double_t diffL=TMath::Sqrt(zlength)*param->GetDiffL()/zwidth;
    //
    // transform angular effect to pad units
    Double_t pky   = ty*padLength/padWidth;
    Double_t pkz   = tz*padLength/zwidth;
    //
    //
    Double_t sy = TMath::Sqrt(rmsy0*rmsy0+diffT*diffT+pky*pky/12.);
    Double_t sz = TMath::Sqrt(rmsz0*rmsz0+diffL*diffL+pkz*pkz/12.);
    Int_t ipad = TMath::Nint(position[7]-1.5);
     //
    (*pcstream)<<"Dump"<<
      "bincont="<<value<<      // bin content
      "val="<<position[0]<<    // parameter difference
      "ipad="<<ipad<<          // pad type
      "dr="<<position[1]<<     //drift
      "ty="<<ty<<              //phi
      "tz="<<tz<<              //theta      
      "p2="<<position[2]<<      //p2
      "p3="<<position[3]<<     //p3
      "p="<<position[4]<<      //p
      "bg="<<position[5]<<     //bg
      "ncl="<<position[6]<<    //ncl
      "type="<<position[7]<<   //tracklet
      "tot="<<position[8]<<    //is tot 
      "sy="<<sy<<              // sigma y
      "sz="<<sz<<              // sigma z
      "pky="<<pky<<            // ky - angle in pad units
      "pkz="<<pkz<<            // kz - angle in pad units
      "diy="<<diffT<<
      "diz="<<diffL<<
      "\n";
  }
  delete pcstream;
}


void AliTPCcalibPID::DumpTrees(){
  //
  // dump the content of histogram to the tree
  //
  AliTPCcalibPID *pid = this;
  DumpTree(pid->GetHistQtot(),"dumpQtot.root");
  DumpTree(pid->GetHistQmax(),"dumpQmax.root");
  DumpTree(pid->GetHistRatioQtot(),"dumpRatioQtot.root");
  DumpTree(pid->GetHistRatioQmax(),"dumpRatioQmax.root");
}



