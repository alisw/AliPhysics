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
  Class for histogramming of cluster residuals
  and fitting of the unlinearities. To be used only for data without
  magnetic field. The laser tracks should be also rejected.
  //

  Track fitting:
  The track is fitted using linear and parabolic model 
  The edge clusters are removed from the fit.
  Edge pad-row - first and last 15 padrows removed from the fit

  Unlinearities fitting:
  Unlinearities at the edge aproximated using two exponential decays.

  Model:
  dz = dz0(r,z) +dr(r,z)*tan(theta) 
  dy =          +dr(r,z)*tan(phi)

   
  .x ~/NimStyle.C
  gSystem->Load("libANALYSIS");
  gSystem->Load("libTPCcalib");
  TFile fcalib("CalibObjects.root");
  TObjArray * array = (TObjArray*)fcalib.Get("TPCCalib");
  AliTPCcalibUnlinearity * calibUnlin = ( AliTPCcalibUnlinearity *)array->FindObject("calibUnlinearity");
  //
  
*/


#include "TLinearFitter.h"

#include "Riostream.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "THnSparse.h"
#include "TList.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TF1.h"
#include "TVectorD.h"
#include "TProfile.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TCut.h"


#include "AliTPCclusterMI.h"
#include "AliTPCseed.h"
#include "AliESDVertex.h"
#include "AliESDEvent.h"
#include "AliESDfriend.h"
#include "AliESDInputHandler.h"
#include "AliAnalysisManager.h"
#include "AliTracker.h"
#include "AliMagF.h"
#include "AliTPCCalROC.h"

#include "AliLog.h"


#include "TTreeStream.h"
#include "AliTPCTracklet.h"
#include "TTimeStamp.h"
#include "AliTPCcalibDB.h"
#include "AliDCSSensorArray.h"
#include "AliDCSSensor.h"
#include "TLinearFitter.h"
#include "AliTPCClusterParam.h"
#include "TStatToolkit.h"


#include "AliTPCcalibUnlinearity.h"
#include "AliTPCPointCorrection.h"


ClassImp(AliTPCcalibUnlinearity)

AliTPCcalibUnlinearity::AliTPCcalibUnlinearity():
  AliTPCcalibBase(),
  fInit(kFALSE),
  fDiffHistoLineY(0),
  fDiffHistoParY(0),
  fDiffHistoLineZ(0),
  fDiffHistoParZ(0),
  fFittersOutR(38),
  fFittersOutZ(38),
  fParamsOutR(38),
  fParamsOutZ(38),
  fErrorsOutR(38),
  fErrorsOutZ(38),
  fDistRPHIPlus(74),
  fDistRPHIMinus(74),
  fFitterQuadrantY(16*38),
  fFitterQuadrantPhi(16*38)
{
  //
  // Defualt constructor
  //
}

AliTPCcalibUnlinearity::AliTPCcalibUnlinearity(const Text_t *name, const Text_t *title):
  AliTPCcalibBase(name,title),
  fInit(kFALSE),
  fDiffHistoLineY(0),
  fDiffHistoParY(0),  
  fDiffHistoLineZ(0),
  fDiffHistoParZ(0),  
  fFittersOutR(38),
  fFittersOutZ(38),
  fParamsOutR(38),
  fParamsOutZ(38),
  fErrorsOutR(38),
  fErrorsOutZ(38),
  fDistRPHIPlus(74),
  fDistRPHIMinus(74),
  fFitterQuadrantY(16*38),
  fFitterQuadrantPhi(16*38)
{
  //
  // Non default constructor
  //
  MakeHisto();
}

AliTPCcalibUnlinearity::~AliTPCcalibUnlinearity(){
  //
  //
  //
  if (fDiffHistoLineZ) delete fDiffHistoLineY;
  if (fDiffHistoParZ)  delete fDiffHistoParY;
  if (fDiffHistoLineZ) delete fDiffHistoLineZ;
  if (fDiffHistoParZ)  delete fDiffHistoParZ;
}





void AliTPCcalibUnlinearity::Process(AliTPCseed *seed) {
  //
  // 
  //  
  if (!fInit) {
    Init();   // work around for PROOF - initialize firs time used
  }
  const Int_t  kdrow=10;
  const Int_t kMinCluster=35;
  if (!seed) return;
  if (TMath::Abs(fMagF)>0.01) return;   // use only data without field
  Int_t counter[72];
  for (Int_t i=0;i<72;i++) counter[i]=0;
  for (Int_t irow=kdrow;irow<159-kdrow;irow++) {
    AliTPCclusterMI *cluster=seed->GetClusterPointer(irow);
    if (!cluster) continue;
    if (cluster->GetType()<0) continue;
    counter[cluster->GetDetector()]++;
  }
  //
  for (Int_t is=0; is<72;is++){
    if (counter[is]>kMinCluster) ProcessDiff(seed, is);
    if (counter[is]>kMinCluster) {
      AlignOROC(seed, is);
    }
  }
}


void AliTPCcalibUnlinearity::ProcessDiff(AliTPCseed *track, Int_t isec){
  //
  // process differnce of the cluster in respect with linear and parabolic fit
  //
  const  Double_t kXmean=134;
  const Int_t     kdrow=10;
  const Float_t kSagitaCut = 1;
  static TLinearFitter fy1(2,"pol1");
  static TLinearFitter fy2(3,"pol2");
  static TLinearFitter fz1(2,"pol1");
  static TLinearFitter fz2(3,"pol2");
  //
  static TVectorD py1(2);
  static TVectorD py2(3);
  static TVectorD pz1(2);
  static TVectorD pz2(3);  
  //  
  fy1.ClearPoints();
  fy2.ClearPoints();
  fz1.ClearPoints();
  fz2.ClearPoints();
  //
  for (Int_t irow=kdrow;irow<159-kdrow;irow++) {
    AliTPCclusterMI *c=track->GetClusterPointer(irow);
    if (!c) continue;
    if (c->GetDetector()!=isec) continue;
    if (c->GetType()<0) continue;
    Double_t dx = c->GetX()-kXmean;
    Double_t x[2]={dx, dx*dx};
    fy1.AddPoint(x,c->GetY());
    fy2.AddPoint(x,c->GetY());
    fz1.AddPoint(x,c->GetZ());
    fz2.AddPoint(x,c->GetZ());
  }
  fy1.Eval(); fy1.GetParameters(py1);
  fy2.Eval(); fy2.GetParameters(py2);
  fz1.Eval(); fz1.GetParameters(pz1);
  fz2.Eval(); fz2.GetParameters(pz2);
  //
  
  for (Int_t irow=0;irow<159;irow++) {
    AliTPCclusterMI *c=track->GetClusterPointer(irow);
    if (!c) continue;
    if (c->GetDetector()!=isec) continue;
    Double_t dx = c->GetX()-kXmean;
    Double_t y1 = py1[0]+py1[1]*dx;
    Double_t y2 = py2[0]+py2[1]*dx+py2[2]*dx*dx;
    //
    Double_t z1 = pz1[0]+pz1[1]*dx;
    Double_t z2 = pz2[0]+pz2[1]*dx+pz2[2]*dx*dx;
    //
    Double_t edgeMinus= -y1+c->GetX()*TMath::Tan(TMath::Pi()/18.);
    Double_t edgePlus =  y1+c->GetX()*TMath::Tan(TMath::Pi()/18.);
    Int_t    npads = AliTPCROC::Instance()->GetNPads(isec,c->GetRow());
    Float_t    dpad  = TMath::Min(c->GetPad(), npads-1- c->GetPad());
    Float_t dy =c->GetY()-y1;
    //
    // Corrections 
    //    
    AliTPCPointCorrection * corr =  AliTPCPointCorrection::Instance();
    Double_t corrclY=0, corrtrY=0, corrR=0;
    corrclY =  corr->RPhiCOGCorrection(isec,c->GetRow(), c->GetPad(),
				       c->GetY(),c->GetY(), c->GetZ(), py1[1],c->GetMax(),2.5); 
    corrtrY =  corr->RPhiCOGCorrection(isec,c->GetRow(), c->GetPad(),
				       c->GetY(),y1, c->GetZ(), py1[1], c->GetMax(),2.5);     
    corrR   = corr->CorrectionOutR0(kFALSE,kFALSE,c->GetX(),c->GetY(),c->GetZ(),isec);
    //
    //
    //
    if (fStreamLevel>1){
      TTreeSRedirector *cstream = GetDebugStreamer();
      if (cstream){
	(*cstream)<<"Diff"<<
	  "run="<<fRun<<              //  run number
	  "event="<<fEvent<<          //  event number
	  "time="<<fTime<<            //  time stamp of event
	  "trigger="<<fTrigger<<      //  trigger
	  "mag="<<fMagF<<             //  magnetic field	  
	  "isec="<<isec<<             // current sector
	  "npads="<<npads<<           // number of pads at given sector
	  "dpad="<<dpad<<             // distance to edge pad
	  //
	  "crR="<<corrR<<              // radial correction
	  "cclY="<<corrclY<<            // edge effect correction cl
	  "ctrY="<<corrtrY<<            // edge effect correction using track
	  //
	  "Cl.="<<c<<
	  "y1="<<y1<<
	  "y2="<<y2<<
	  "z1="<<z1<<
	  "z2="<<z2<<
	  "py1.="<<&py1<<
	  "py2.="<<&py2<<
	  "pz1.="<<&pz1<<
	  "pz2.="<<&pz2<<
	  "eP="<<edgePlus<<
	  "eM="<<edgeMinus<<
	  "\n";
      }
    }
    if (TMath::Min(edgeMinus,edgePlus)<6){
      AddPointRPHI(isec,c->GetX(),c->GetY(),c->GetZ(),y1,z1,py1[1],pz1[1],1);
      TTreeSRedirector *cstream = GetDebugStreamer();
      if (TMath::Abs(dy)<2 && cstream){
	(*cstream)<<"rphi"<<
	  "isec="<<isec<<             // current sector
	  "npads="<<npads<<           // number of pads at given sector
	  "dpad="<<dpad<<             // distance to edge pad
	  "eP="<<edgePlus<<           // distance to edge 
	  "eM="<<edgeMinus<<          // distance to edge minus
	  //
	  "crR="<<corrR<<              // radial correction
	  "cclY="<<corrclY<<            // edge effect correction cl
	  "ctrY="<<corrtrY<<            // edge effect correction using track
	  //
	  "dy="<<dy<<                 // dy
	  "Cl.="<<c<<
	  "y1="<<y1<<
	  "y2="<<y2<<
	  "z1="<<z1<<
	  "z2="<<z2<<
	  "py1.="<<&py1<<
	  "pz1.="<<&pz1<<
	  "\n";
      }
    }
    if (c->GetType()<0) continue;  // don't use edge rphi cluster
    //
    //
    Double_t x[10];
    x[0]=c->GetDetector();
    x[1]=c->GetRow();
    x[2]=c->GetY()/c->GetX();
    x[3]=c->GetZ();
    //
    // apply sagita cut
    //
    if (TMath::Abs(py2[2]*150*150/4)<kSagitaCut && 
	TMath::Abs(pz2[2]*150*150/4)<kSagitaCut){
      //
      x[4]=py1[1];
      x[5]=c->GetY()-y1;
      fDiffHistoLineY->Fill(x);
      x[5]=c->GetY()-y2;
      //fDiffHistoParY->Fill(x);
      //
      x[4]=pz1[1];
      x[5]=c->GetY()-z1;
      fDiffHistoLineZ->Fill(x);
      x[5]=c->GetY()-z2;
      //fDiffHistoParZ->Fill(x);      
      // Apply sagita cut
      AddPoint(isec,c->GetX(),c->GetY(),c->GetZ(),y1,z1,py1[1],pz1[1],1);
    }
	     
  }
}


void AliTPCcalibUnlinearity::AlignOROC(AliTPCseed *track, Int_t isec){
  //
  // fit the tracklet in OROC sepatately in Quadrant
  //
  //
  if (isec<36) return; 
  const Int_t kMinClusterF=35;
  const Int_t kMinClusterQ=10;
  //
  const  Int_t     kdrow1 =8;       // rows to skip at the end      
  const  Int_t     kdrow0 =2;        // rows to skip at beginning  
  const  Float_t   kedgey=3;
  const  Float_t   kMaxDist=0.5;
  const  Float_t   kMaxCorrY=0.1;
  const  Float_t   kPRFWidth = 0.6;   //cut 2 sigma of PRF
  //
  //
  AliTPCROC * roc = AliTPCROC::Instance();
  AliTPCPointCorrection * corr =  AliTPCPointCorrection::Instance();
  const  Double_t kXmean=roc->GetPadRowRadii(isec,53);
  //
  // full track fit parameters
  // 
  TLinearFitter fyf(2,"pol1");
  TLinearFitter fzf(2,"pol1");
  TVectorD pyf(2), peyf(2),pzf(2), pezf(2);
  Int_t nf=0;
  //
  // make full fit as reference
  //
  for (Int_t iter=0; iter<2; iter++){
    fyf.ClearPoints();
    for (Int_t irow=kdrow0;irow<159-kdrow1;irow++) {
      AliTPCclusterMI *c=track->GetClusterPointer(irow);
      if (!c) continue;
      if (c->GetDetector()!=isec) continue;
      if (c->GetRow()<kdrow0) continue;
      Double_t dx = c->GetX()-kXmean;
      Double_t x[2]={dx, dx*dx};
      if (iter==0 &&c->GetType()<0) continue;
      if (iter==1){	
	Double_t yfit  =  pyf[0]+pyf[1]*dx;
	Double_t dedge =  c->GetX()*TMath::Tan(TMath::Pi()/18.)-TMath::Abs(yfit);
	if (TMath::Abs(c->GetY()-yfit)>kMaxDist) continue;
	if (dedge<kedgey) continue;
	Double_t corrtrY =  
	  corr->RPhiCOGCorrection(isec,c->GetRow(), c->GetPad(),
				  c->GetY(),yfit, c->GetZ(), pyf[1], c->GetMax(),2.5);
	if (TMath::Abs(corrtrY)>kMaxCorrY) continue;
      }
      fyf.AddPoint(x,c->GetY(),0.1);
      fzf.AddPoint(x,c->GetZ(),0.1);
    }
    nf = fyf.GetNpoints();
    if (nf<kMinClusterF) return;   // not enough points - skip 
    fyf.Eval(); 
    fyf.GetParameters(pyf); 
    fyf.GetErrors(peyf);
    fzf.Eval(); 
    fzf.GetParameters(pzf); 
    fzf.GetErrors(pezf);
  }
  //
  // Make Fitters and params for 3 fitters
  //
  TLinearFitter *fitters[4];
  Int_t npoints[4];
  TVectorD params[4];
  TVectorD errors[4];
  Double_t chi2C[4];
  for (Int_t i=0;i<4;i++) {
    fitters[i] = new TLinearFitter(2,"pol1");
    npoints[i]=0;
    params[i].ResizeTo(2);
    errors[i].ResizeTo(2);
  }
  //
  // Update fitters
  //
  for (Int_t irow=kdrow0;irow<159-kdrow1;irow++) {
    AliTPCclusterMI *c=track->GetClusterPointer(irow);
    if (!c) continue;
    if (c->GetDetector()!=isec) continue;
    if (c->GetRow()<kdrow0) continue;      
    Double_t dx = c->GetX()-kXmean;
    Double_t x[2]={dx, dx*dx};
    Double_t yfit  =  pyf[0]+pyf[1]*dx;
    Double_t dedge =  c->GetX()*TMath::Tan(TMath::Pi()/18.)-TMath::Abs(yfit);
    if (TMath::Abs(c->GetY()-yfit)>kMaxDist) continue;
    if (dedge<kedgey) continue;
    Double_t corrtrY =  
      corr->RPhiCOGCorrection(isec,c->GetRow(), c->GetPad(),
			      c->GetY(),yfit, c->GetZ(), pyf[1], c->GetMax(),2.5);
    if (TMath::Abs(corrtrY)>kMaxCorrY) continue;  
    if (dx<0){
      if (yfit<-kPRFWidth) fitters[0]->AddPoint(x,c->GetY(),0.1);
      if (yfit>kPRFWidth) fitters[1]->AddPoint(x,c->GetY(),0.1);
    }
    if (dx>0){
      if (yfit<-kPRFWidth) fitters[2]->AddPoint(x,c->GetY(),0.1);
      if (yfit>kPRFWidth) fitters[3]->AddPoint(x,c->GetY(),0.1);
    }
  }

  for (Int_t i=0;i<4;i++) {
    npoints[i] = fitters[i]->GetNpoints();
    if (npoints[i]>=kMinClusterQ){
      fitters[i]->Eval();
      Double_t chi2Fac = TMath::Sqrt(fitters[i]->GetChisquare()/(fitters[i]->GetNpoints()-2));
      chi2C[i]=chi2Fac;
      fitters[i]->GetParameters(params[i]);
      fitters[i]->GetErrors(errors[i]);
      // renormalize errors
      errors[i][0]*=chi2Fac;
      errors[i][1]*=chi2Fac;
    }
  }
  //
  // Fill fitters
  //
  for (Int_t i0=0;i0<4;i0++){
    for (Int_t i1=0;i1<4;i1++){
      if (i0==i1) continue;
      if(npoints[i0]<kMinClusterQ) continue;
      if(npoints[i1]<kMinClusterQ) continue;
      Int_t index0=i0*4+i1;
      Double_t xy[1];
      Double_t xfi[1];
      xy[0] = pyf[1];   
      xfi[0] = pyf[1];
      //
      Double_t dy = params[i1][0]-params[i0][0];
      Double_t sy = TMath::Sqrt(errors[i1][0]*errors[i1][0]+errors[i0][0]*errors[i0][0]);
      Double_t dphi = params[i1][1]-params[i0][1];
      Double_t sphi = TMath::Sqrt(errors[i1][1]*errors[i1][1]+errors[i0][1]*errors[i0][1]);
      //
      Int_t sector = isec-36;
      //
      if (TMath::Abs(dy/(sy+0.1))>5.) continue;          // 5 sigma cut
      if (TMath::Abs(dphi/(sphi+0.004))>5.) continue;    // 5 sigma cut 
      TLinearFitter * fitterY,*fitterPhi;
      fitterY = (TLinearFitter*) fFitterQuadrantY.At(16*sector+index0);
      if (TMath::Abs(dy/(sy+0.1))<5.)  fitterY->AddPoint(xy,dy,sy);
      fitterY = (TLinearFitter*) fFitterQuadrantY.At(16*36+index0);
      if (TMath::Abs(dy/(sy+0.1))<5.)  fitterY->AddPoint(xy,dy,sy);
      //
      fitterPhi = (TLinearFitter*) fFitterQuadrantPhi.At(16*sector+index0);
      if (TMath::Abs(dphi/(sphi+0.001))<5.)  fitterPhi->AddPoint(xfi,dphi,sphi);
      fitterPhi = (TLinearFitter*) fFitterQuadrantPhi.At(16*36+index0);
      if (TMath::Abs(dphi/(sphi+0.001))<5.)  fitterPhi->AddPoint(xfi,dphi,sphi);
    }
  }
  //
  // Dump debug information
  //
  if (fStreamLevel>0){    
    TTreeSRedirector *cstream = GetDebugStreamer();      
    Int_t sector = isec-36;
    if (cstream){      
      for (Int_t i0=0;i0<4;i0++){
	for (Int_t i1=i0;i1<4;i1++){
	  if (i0==i1) continue;
	  if(npoints[i0]<kMinClusterQ) continue;
	  if(npoints[i1]<kMinClusterQ) continue;
	  (*cstream)<<"fitQ"<<
	    "run="<<fRun<<              //  run number
	    "event="<<fEvent<<          //  event number
	    "time="<<fTime<<            //  time stamp of event
	    "trigger="<<fTrigger<<      //  trigger
	    "mag="<<fMagF<<             //  magnetic field	  
	    "sec="<<sector<<            // current sector from 0 to 36
	    "isec="<<isec<<             //  current sector
	    // full fit
	    "nf="<<nf<<                 //  total number of points
	    "pyf.="<<&pyf<<             //  full OROC fit y
	    "pzf.="<<&pzf<<             //  full OROC fit z
	    // quadrant fit
	    "i0="<<i0<<                 // quadrant number
	    "i1="<<i1<<                 
	    "n0="<<npoints[i0]<<        // number of points
	    "n1="<<npoints[i1]<<
	    "p0.="<<&params[i0]<<       // parameters
	    "p1.="<<&params[i1]<< 
	    "e0.="<<&errors[i0]<<       // errors
	    "e1.="<<&errors[i1]<<
	    "chi0="<<chi2C[i0]<<       // chi2s
	    "chi1="<<chi2C[i1]<<

	    "\n";
	}    
      }
    }
  }
}





void AliTPCcalibUnlinearity::MakeHisto(){
  //
  //
  //
  TString  axisName[10];
  Double_t xmin[10],  xmax[10];
  Int_t    nbins[10];
  //
  //
  axisName[0]="sector";
  xmin[0]=0; xmax[0]=72; nbins[0]=72;
  //          
  axisName[1]="row";
  xmin[1]=0; xmax[1]=96; nbins[1]=96;
  //
  axisName[2]="tphi";            // tan phi - ly/lx 
  xmin[2]=-0.17; xmax[2]=0.17; nbins[2]=5;
  //  
  axisName[3]="lz";              // global z        
  xmin[3]=-250; xmax[3]=250; nbins[3]=10;
  //
  axisName[4]="k";              // tangent - in angle
  xmin[4]=-0.5; xmax[4]=0.5; nbins[4]=10;
  //
  //
  axisName[5]="delta";          // delta
  xmin[5]=-0.5; xmax[5]=0.5; nbins[5]=100;
  //
  //
  fDiffHistoLineY = new THnSparseF("hnDistHistoLineY","hnDistHistoLineY",6,nbins,xmin,xmax);
  fDiffHistoParY = new THnSparseF("hnDistHistoParY","hnDistHistoParY",6,nbins,xmin,xmax);
  fDiffHistoLineZ = new THnSparseF("hnDistHistoLineZ","hnDistHistoLineZ",6,nbins,xmin,xmax);
  fDiffHistoParZ = new THnSparseF("hnDistHistoParZ","hnDistHistoParZ",6,nbins,xmin,xmax);
 //  for (Int_t i=0; i<8;i++){
//     fDiffHistoLineY->GetAxis(i)->SetName(axisName[i].Data());
//     fDiffHistoParY->GetAxis(i)->SetName(axisName[i].Data());
//     fDiffHistoLineY->GetAxis(i)->SetTitle(axisName[i].Data());
//     fDiffHistoParY->GetAxis(i)->SetTitle(axisName[i].Data());
//     fDiffHistoLineZ->GetAxis(i)->SetName(axisName[i].Data());
//     fDiffHistoParZ->GetAxis(i)->SetName(axisName[i].Data());
//     fDiffHistoLineZ->GetAxis(i)->SetTitle(axisName[i].Data());
//     fDiffHistoParZ->GetAxis(i)->SetTitle(axisName[i].Data());
//   }
  
  //
  // 
  //
  char hname[1000];
  TH2F * his=0;
  for (Int_t isec=0;isec<74;isec++){
    sprintf(hname,"DeltaRPhiPlus_Sector%d",isec);
    his = new TH2F(hname,hname,100,1,4,100,-0.5,0.5);
    his->SetDirectory(0);
    fDistRPHIPlus.AddAt(his,isec);
    sprintf(hname,"DeltaRPhiMinus_Sector%d",isec);
    his = new TH2F(hname,hname,100,1,4,100,-0.5,0.5);
    his->SetDirectory(0);
    fDistRPHIMinus.AddAt(his,isec);
  }
}


void AliTPCcalibUnlinearity::Terminate(){
  //
  // Terminate function
  // call base terminate + Eval of fitters
  //
  Info("AliTPCcalibUnlinearities","Terminate");
  EvalFitters();
  DumpTree();
  AliTPCcalibBase::Terminate();
}

Long64_t AliTPCcalibUnlinearity::Merge(TCollection *list) {
  //
  // merge results
  //
  TIterator* iter = list->MakeIterator();
  AliTPCcalibUnlinearity* cal = 0;
  //
  while ((cal = (AliTPCcalibUnlinearity*)iter->Next())) {
    if (!cal->InheritsFrom(AliTPCcalibUnlinearity::Class())) {
      return -1;
    }    
    Add(cal);
  }
  return 0;  
}

void    AliTPCcalibUnlinearity::Add(AliTPCcalibUnlinearity * calib){
  //
  //
  //
  if (!fInit) Init();
  if (fDiffHistoLineY && calib->fDiffHistoLineY){
    fDiffHistoLineY->Add(calib->fDiffHistoLineY);
  }
  if (fDiffHistoParY && calib->fDiffHistoParY){
    fDiffHistoParY->Add(calib->fDiffHistoParY);
  }
  if (fDiffHistoLineZ && calib->fDiffHistoLineZ){
    fDiffHistoLineZ->Add(calib->fDiffHistoLineZ);
  }
  if (fDiffHistoParZ && calib->fDiffHistoParZ){
    fDiffHistoParZ->Add(calib->fDiffHistoParZ);
  }
  for (Int_t isec=0;isec<38;isec++){
    TLinearFitter *f0r = (TLinearFitter*)fFittersOutR.At(isec);
    TLinearFitter *f1r = (TLinearFitter*)(calib->fFittersOutR.At(isec));
    if (f0r&&f1r) f0r->Add(f1r);
    TLinearFitter *f0z = (TLinearFitter*)fFittersOutZ.At(isec);
    TLinearFitter *f1z = (TLinearFitter*)(calib->fFittersOutZ.At(isec));
    if (f0z&&f1z) f0z->Add(f1z);
  }

  for (Int_t isec=0;isec<16*38;isec++){
    TLinearFitter *f0y = (TLinearFitter*)fFitterQuadrantY.At(isec);
    TLinearFitter *f1y = (TLinearFitter*)(calib->fFitterQuadrantY.At(isec));
    if (f0y&&f1y) f0y->Add(f1y);
    TLinearFitter *f0phi = (TLinearFitter*)fFitterQuadrantPhi.At(isec);
    TLinearFitter *f1phi = (TLinearFitter*)(calib->fFitterQuadrantPhi.At(isec));
    if (f0phi&&f1phi) f0phi->Add(f1phi);
  }
  
  
  for (Int_t isec=0;isec<74;isec++){
    TH2F * h0p= (TH2F*)(fDistRPHIPlus.At(isec));
    TH2F * h1p= (TH2F*)(calib->fDistRPHIPlus.At(isec));
    if (h0p&&h1p) h0p->Add(h1p);
    TH2F * h0m= (TH2F*)(fDistRPHIMinus.At(isec));
    TH2F * h1m= (TH2F*)(calib->fDistRPHIMinus.At(isec));
    if (h0m&&h1m) h0m->Add(h1m);
  }

 
}



void AliTPCcalibUnlinearity::DumpTree(const char *fname){
  //
  //
  // 
  TTreeSRedirector *cstream =new TTreeSRedirector(fname);
  if (!cstream) return;
  //  
  THnSparse *his=0;  
  Double_t position[10];
  Double_t value; 
  Int_t *bins = new Int_t[10];
  //
  for (Int_t ihis=0; ihis<2; ihis++){
    his =  (ihis==0)? fDiffHistoLineY:fDiffHistoParY;
    if (!his) continue;
    //
    Int_t ndim = his->GetNdimensions();
    //
    for (Long64_t i = 0; i < his->GetNbins(); ++i) {
      value = his->GetBinContent(i, bins);
      for (Int_t idim = 0; idim < ndim; idim++) {
	position[idim] = his->GetAxis(idim)->GetBinCenter(bins[idim]);
      }      
      (*cstream)<<"Resol"<<
	"ihis="<<ihis<<
	"bincont="<<value;
      for (Int_t idim=0;idim<ndim;idim++){
	(*cstream)<<"Resol"<<Form("%s=", his->GetAxis(idim)->GetName())<<position[idim];
      }
      (*cstream)<<"Resol"<<
	"\n";      
    }
  }
  delete cstream;
}


void AliTPCcalibUnlinearity::AddPoint(Int_t sector, Double_t cx, Double_t cy, Double_t cz, Double_t ty, Double_t tz,  Double_t ky, Double_t kz, Int_t npoints){
  //
  //
  // Simple distortion fit in outer sectors
  //
  // Model - 2 exponential decay toward the center of chamber
  //       - Decay length 10 and 5 cm
  //       - Decay amplitude - Z dependent 
  //
 
  Double_t side   = (-1+((sector%36)<18)*2); // A side +1   C side -1
  Double_t dzt    = (cz-tz)*side;
  Double_t radius = TMath::Sqrt(cx*cx+cy*cy);  
  AliTPCROC * roc = AliTPCROC::Instance();
  Double_t xout = roc->GetPadRowRadiiUp(roc->GetNRows(37)-1);
  Double_t dout = xout-radius;
  if (dout>30) return;
  //drift
  Double_t dr   = 0.5 - TMath::Abs(cz/250.);
  Double_t dy = cy-ty;
  Double_t eout10 = TMath::Exp(-dout/10.);
  Double_t eout5 = TMath::Exp(-dout/5.);
  //
  Double_t eoutp  = (eout10+eout5)*0.5;
  Double_t eoutm  = (eout10-eout5)*0.5;

  //
  Double_t xxxR[6], xxxZ[6], xxxRZ[6];
  //TString fstring="";
  xxxZ[0]=eoutp;                //fstring+="eoutp++";  
  xxxZ[1]=eoutm;               //fstring+="eoutm++";  
  xxxZ[2]=dr*eoutp;            //fstring+="dr*eoutp++";  
  xxxZ[3]=dr*eoutm;            //fstring+="dr*eoutm++";  
  xxxZ[4]=dr*dr*eoutp;         //fstring+="dr*dr*eoutp++";  
  xxxZ[5]=dr*dr*eoutm;         //fstring+="dr*dr*eoutm++";    
  //
  xxxR[0]=ky*eoutp;                //fstring+="eoutp++";  
  xxxR[1]=ky*eoutm;               //fstring+="eoutm++";  
  xxxR[2]=ky*dr*eoutp;            //fstring+="dr*eoutp++";  
  xxxR[3]=ky*dr*eoutm;            //fstring+="dr*eoutm++";  
  xxxR[4]=ky*dr*dr*eoutp;         //fstring+="dr*dr*eoutp++";  
  xxxR[5]=ky*dr*dr*eoutm;         //fstring+="dr*dr*eoutm++";    
  //
  xxxRZ[0]=kz*eoutp;                //fstring+="eoutp++";  
  xxxRZ[1]=kz*eoutm;               //fstring+="eoutm++";  
  xxxRZ[2]=kz*dr*eoutp;            //fstring+="dr*eoutp++";  
  xxxRZ[3]=kz*dr*eoutm;            //fstring+="dr*eoutm++";  
  xxxRZ[4]=kz*dr*dr*eoutp;         //fstring+="dr*dr*eoutp++";  
  xxxRZ[5]=kz*dr*dr*eoutm;         //fstring+="dr*dr*eoutm++";    
  //
  TLinearFitter * fitter=0;
  Double_t err=0.1;
  for (Int_t ip=0; ip<npoints; ip++){
    //
    fitter = (TLinearFitter*)fFittersOutR.At(sector%36);
    fitter->AddPoint(xxxR,dy,err);
    //fitter->AddPoint(xxxRZ,dz);
    fitter = (TLinearFitter*)fFittersOutZ.At(sector%36);
    fitter->AddPoint(xxxZ,dzt,err);
    //
    if (side>0){
      fitter = (TLinearFitter*)fFittersOutR.At(36);
      fitter->AddPoint(xxxR,dy,err);
      //fitter->AddPoint(xxxRZ,dz);
      fitter = (TLinearFitter*)fFittersOutZ.At(36);
      fitter->AddPoint(xxxZ,dzt,err);
    }
    if (side<0){
      fitter = (TLinearFitter*)fFittersOutR.At(37);
      fitter->AddPoint(xxxR,dy,err);
      //fitter->AddPoint(xxxRZ,dz);
      fitter = (TLinearFitter*)fFittersOutZ.At(37);
      fitter->AddPoint(xxxZ,dzt,err);
    }
  }
}
void AliTPCcalibUnlinearity::AddPointRPHI(Int_t sector, Double_t cx, Double_t cy, Double_t cz, Double_t ty, Double_t tz,  Double_t /*ky*/, Double_t /*kz*/, Int_t npoints){
  //
  //
  //
  Float_t kMaxDz = 0.5;  // cut on z in cm
  Double_t dy = cy-ty;
  Double_t dz = cz-tz;
  if (TMath::Abs(dz)>kMaxDz) return;
  Double_t edgeMinus= -ty+cx*TMath::Tan(TMath::Pi()/18.);
  Double_t edgePlus =  ty+cx*TMath::Tan(TMath::Pi()/18.);
  //
  TH2F * his =0;
  his = (TH2F*)fDistRPHIPlus.At(sector);
  his->Fill(edgePlus,dy,npoints);
  his = (TH2F*)((sector<36)? fDistRPHIPlus.At(72):fDistRPHIPlus.At(73));
  his->Fill(edgePlus,dy,npoints);
  his = (TH2F*)fDistRPHIMinus.At(sector);
  his->Fill(edgeMinus,dy,npoints);
  his = (TH2F*)((sector<36)? fDistRPHIMinus.At(72):fDistRPHIMinus.At(73));
  his->Fill(edgeMinus,dy,npoints);
}



void AliTPCcalibUnlinearity::Init(){
  //
  //
  //
  //
  MakeHisto();
  // Make outer fitters
  TLinearFitter * fitter = 0;
  for (Int_t ifit=0; ifit<38; ifit++){
    fitter = new TLinearFitter(7,"hyp6");
    fitter->StoreData(kFALSE);
    fFittersOutR.AddAt(fitter,ifit);
    fitter = new TLinearFitter(7,"hyp6");
    fitter->StoreData(kFALSE);
    fFittersOutZ.AddAt(fitter,ifit);
  }
  for (Int_t ifit=0; ifit<16*38;ifit++){
    fitter = new TLinearFitter(2,"hyp1");
    fitter->StoreData(kFALSE);
    fFitterQuadrantY.AddAt(fitter,ifit);
    fitter = new TLinearFitter(2,"hyp1");
    fitter->StoreData(kFALSE);
    fFitterQuadrantPhi.AddAt(fitter,ifit);    
  }  
  fInit= kTRUE;
}

void AliTPCcalibUnlinearity::EvalFitters(){
  //
  // 
  // Evaluate fitters
  // 
  TLinearFitter * fitter = 0;
  TVectorD vec;
  for (Int_t ifit=0; ifit<38; ifit++){
    fitter=(TLinearFitter*)fFittersOutR.At(ifit);
    if (fitter) {
      fitter->Eval();
      fitter->GetParameters(vec);
      fParamsOutR.AddAt(vec.Clone(),ifit);
      fitter->GetErrors(vec);
      fErrorsOutR.AddAt(vec.Clone(),ifit);
    }
    fitter=(TLinearFitter*)fFittersOutZ.At(ifit);
    if (fitter) {
      fitter->Eval();
      fitter->GetParameters(vec);
      fParamsOutZ.AddAt(vec.Clone(),ifit);
      fitter->GetErrors(vec);
      fErrorsOutZ.AddAt(vec.Clone(),ifit);
    }
  }
}





void AliTPCcalibUnlinearity::ProcessTree(TTree * tree, Long64_t nmaxPoints){
  //
  //
  //
  //  TTree * tree = chainUnlinP;
  AliTPCcalibUnlinearity *calib = this;
  //
  AliTPCclusterMI * cl=0;
  Double_t ty,tz;
  TVectorD *vy=0, *vz=0;
  TVectorD *vy2=0, *vz2=0;
  Long64_t nentries = tree->GetEntries();
  if (nmaxPoints>0) nentries = TMath::Min(nentries,nmaxPoints);
  //
  {
  for (Long64_t i=0; i<nentries; i++){
    tree->SetBranchAddress("Cl.",&cl);
    tree->SetBranchAddress("y1",&ty);
    tree->SetBranchAddress("z1",&tz);
    tree->SetBranchAddress("py1.",&vy);
    tree->SetBranchAddress("pz1.",&vz);
    tree->SetBranchAddress("py2.",&vy2);
    tree->SetBranchAddress("pz2.",&vz2);
    //
    tree->GetEntry(i);
    if (!cl) continue;
    if (i%10000==0) printf("%d\n",(Int_t)i);
    Int_t row= cl->GetRow();
    if (cl->GetDetector()>36) row+=64;
    calib->AddPointRPHI(cl->GetDetector(),cl->GetX(),cl->GetY(),cl->GetZ(),
		      ty,tz,(*vy)[1],(*vz)[1],1);

    if (cl->GetType()<0) continue; 
    Double_t dy = cl->GetY()-ty;
    Double_t dz = cl->GetZ()-tz;
    //
    //
    if (TMath::Abs(dy)>0.25) continue;
    if (TMath::Abs(dz)>0.25) continue;
    
    if (TMath::Abs((*vy2)[2]*150*150/4)>1) continue;
    if (TMath::Abs((*vz2)[2]*150*150/4)>1) continue;
      // Apply sagita cut
    if (cl->GetType()>=0) {
      calib->AddPoint(cl->GetDetector(),cl->GetX(),cl->GetY(),cl->GetZ(),
		      ty,tz,(*vy)[1],(*vz)[1],1);
    }
  }
  }
}


void  AliTPCcalibUnlinearity::MakeQPosNormAll(TTree * chainUnlinD, AliTPCClusterParam * param, Int_t maxPoints){
  //
  // Make position corrections
  // for the moment Only using debug streamer 
  // chainUnlinD  - debug tree
  // param     - parameters to be updated
  // maxPoints - maximal number of points using for fit
  // verbose   - print info flag
  //
  // Current implementation - only using debug streamers
  // 
  
  /*    
  //Defaults
  Int_t maxPoints=100000;
  */
  /*
    Usage: 
    //0. Load libraries
    gSystem->Load("libANALYSIS");
    gSystem->Load("libSTAT");
    gSystem->Load("libTPCcalib");
    

    //1. Load Parameters to be modified:
    //e.g:
    
    AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
    AliCDBManager::Instance()->SetRun(0) 
    AliTPCClusterParam * param = AliTPCcalibDB::Instance()->GetClusterParam();
    //
    //2. Load the debug streamer
    gSystem->AddIncludePath("-I$ALICE_ROOT/TPC/macros");
    gROOT->LoadMacro("$ALICE_ROOT/TPC/macros/AliXRDPROOFtoolkit.cxx+")
    AliXRDPROOFtoolkit tool;  
    TChain * chainUnlin = tool.MakeChain("unlin.txt","Resol",0,10200);
    chainUnlin->Lookup();
    TChain * chainUnlinD = tool.MakeChain("unlin.txt","Diff",0,10200);
    chainUnlinD->Lookup();

    //3. Do fits and store results
    // 
    AliTPCcalibUnlinearity::MakeQPosNormAll(chainUnlinD,param,200000,0) ;
    TFile f("paramout.root","recreate");
    param->Write("clusterParam");
    f.Close();
    //4. Verification
    TFile f2("paramout.root");
    AliTPCClusterParam *param2 = (AliTPCClusterParam*)f2.Get("clusterParam");
    param2->SetInstance(param2);
    chainUnlinD->Draw("fitZ0:AliTPCClusterParam::SPosCorrection(1,0,Cl.fPad,Cl.fTimeBin,Cl.fZ,Cl.fSigmaY2,Cl.fSigmaZ2,Cl.fMax)","Cl.fDetector<36","",10000); // should be line 
    
   */


  TStatToolkit toolkit;
  Double_t chi2;
  TVectorD fitParamY0;
  TVectorD fitParamY1;
  TVectorD fitParamZ0;
  TVectorD fitParamZ1;
  TMatrixD covMatrix;
  Int_t npoints;
  
  chainUnlinD->SetAlias("dp","(-1+(Cl.fZ>0)*2)*((Cl.fPad-int(Cl.fPad))-0.5)");
  chainUnlinD->SetAlias("dt","(-1+(Cl.fZ>0)*2)*((Cl.fTimeBin-0.66-int(Cl.fTimeBin-0.66))-0.5)");
  chainUnlinD->SetAlias("sp","(sin(dp*pi)-dp*pi)");
  chainUnlinD->SetAlias("st","(sin(dt)-dt)");
  chainUnlinD->SetAlias("dq","(-1+(Cl.fZ>0)*2)/Cl.fMax");
  //
  chainUnlinD->SetAlias("di","sqrt(1.-abs(Cl.fZ/250.))");
  //
  //
  //
  TCut cutE("Cl.fType>=0");
  TCut cutDy("abs(Cl.fY-y1)<0.4");
  TCut cutDz("abs(Cl.fZ-z1)<0.4");
  TCut cutY2("abs(py2.fElements[2])*150^2/4<1");
  TCut cutZ2("abs(pz2.fElements[2])*150^2/4<1");
  TCut cutAll = cutE+cutDy+cutDz+cutY2+cutZ2;
  //
  TCut cutA("Cl.fZ>0");
  TCut cutC("Cl.fZ<0");

  TString fstringY="";  
  //
  fstringY+="(dp)++";            //1
  fstringY+="(dp)*di++";         //2
  fstringY+="(sp)++";            //3
  fstringY+="(sp)*di++";         //4
  fstringY+="(dq)++";            //5
  TString fstringZ="";  
  fstringZ+="(dt)++";            //1
  fstringZ+="(dt)*di++";         //2
  fstringZ+="(st)++";            //3
  fstringZ+="(st)*di++";         //4
  fstringZ+="(dq)++";            //5
  //
  // Z corrections
  //
  TString *strZ0 = toolkit.FitPlane(chainUnlinD,"(Cl.fZ-z2):1",fstringZ.Data(), "Cl.fDetector<36"+cutAll, chi2,npoints,fitParamZ0,covMatrix,-1,0,maxPoints);
  printf("Z0 - chi2/npoints = %f\n",TMath::Sqrt(chi2/npoints));
  strZ0->Tokenize("++")->Print();
  param->fPosZcor[0] = (TVectorD*) fitParamZ0.Clone();
  //
  TString *strZ1 = toolkit.FitPlane(chainUnlinD,"(Cl.fZ-z2):1",fstringZ.Data(), "Cl.fDetector>36"+cutAll, chi2,npoints,fitParamZ1,covMatrix,-1,0,maxPoints);
  //
  printf("Z1 - chi2/npoints = %f\n",TMath::Sqrt(chi2/npoints));
  strZ1->Tokenize("++")->Print();
  param->fPosZcor[1] = (TVectorD*) fitParamZ1.Clone();
  param->fPosZcor[2] = (TVectorD*) fitParamZ1.Clone();
  //
  // Y corrections
  //   
  TString *strY0 = toolkit.FitPlane(chainUnlinD,"(Cl.fY-y2):1",fstringY.Data(), "Cl.fDetector<36"+cutAll, chi2,npoints,fitParamY0,covMatrix,-1,0,maxPoints);
  printf("Y0 - chi2/npoints = %f\n",TMath::Sqrt(chi2/npoints)); 
  strY0->Tokenize("++")->Print();
  param->fPosYcor[0] = (TVectorD*) fitParamY0.Clone();
  

  TString *strY1 = toolkit.FitPlane(chainUnlinD,"(Cl.fY-y2):1",fstringY.Data(), "Cl.fDetector>36"+cutAll, chi2,npoints,fitParamY1,covMatrix,-1,0,maxPoints);
  //
  printf("Y1 - chi2/npoints = %f\n",TMath::Sqrt(chi2/npoints));  
  strY1->Tokenize("++")->Print();
  param->fPosYcor[1] = (TVectorD*) fitParamY1.Clone();
  param->fPosYcor[2] = (TVectorD*) fitParamY1.Clone();
  //
  //
  //
  chainUnlinD->SetAlias("fitZ0",strZ0->Data());
  chainUnlinD->SetAlias("fitZ1",strZ1->Data());
  chainUnlinD->SetAlias("fitY0",strY0->Data());
  chainUnlinD->SetAlias("fitY1",strY1->Data());
}







