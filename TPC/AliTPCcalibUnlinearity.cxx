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
#include "AliTPCcalibLaser.h"
#include "AliDCSSensorArray.h"
#include "AliDCSSensor.h"
#include "TLinearFitter.h"
#include "AliTPCClusterParam.h"
#include "TStatToolkit.h"


#include "AliTPCcalibUnlinearity.h"

AliTPCcalibUnlinearity* AliTPCcalibUnlinearity::fgInstance = 0;

ClassImp(AliTPCcalibUnlinearity)

AliTPCcalibUnlinearity::AliTPCcalibUnlinearity():
  AliTPCcalibBase(),
  fDiffHistoLine(0),
  fDiffHistoPar(0),
  fFittersOutR(38),
  fFittersOutZ(38),
  fParamsOutR(38),
  fParamsOutZ(38),
  fErrorsOutR(38),
  fErrorsOutZ(38)
{
  //
  // Defualt constructor
  //
}

AliTPCcalibUnlinearity::AliTPCcalibUnlinearity(const Text_t *name, const Text_t *title):
  AliTPCcalibBase(name,title),
  fDiffHistoLine(0),
  fDiffHistoPar(0),  
  fFittersOutR(38),
  fFittersOutZ(38),
  fParamsOutR(38),
  fParamsOutZ(38),
  fErrorsOutR(38),
  fErrorsOutZ(38)
{
  //
  // Non default constructor
  //
  MakeFitters();
}

AliTPCcalibUnlinearity::~AliTPCcalibUnlinearity(){
  //
  //
  //
  if (fDiffHistoLine) delete fDiffHistoLine;
  if (fDiffHistoPar)  delete fDiffHistoPar;
}


AliTPCcalibUnlinearity* AliTPCcalibUnlinearity::Instance()
{
  //
  // Singleton implementation
  // Returns an instance of this class, it is created if neccessary
  //
  if (fgInstance == 0){
    fgInstance = new AliTPCcalibUnlinearity();
  }
  return fgInstance;
}




void AliTPCcalibUnlinearity::Process(AliTPCseed *seed) {
  //
  // 
  //  
  const Int_t  kdrow=10;
  const Int_t kMinCluster=35;
  if (!fDiffHistoLine) MakeHisto();    
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
  }
}


void AliTPCcalibUnlinearity::ProcessDiff(AliTPCseed *track, Int_t isec){
  //
  // process differnce of the cluster in respect with linear and parabolic fit
  //
  const  Double_t kXmean=134;
  const Int_t     kdrow=10;
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
    if (c->GetType()<0) continue;
    Double_t dx = c->GetX()-kXmean;
    Double_t y1 = py1[0]+py1[1]*dx;
    Double_t y2 = py2[0]+py2[1]*dx+py2[2]*dx*dx;
    //
    Double_t z1 = pz1[0]+pz1[1]*dx;
    Double_t z2 = pz2[0]+pz2[1]*dx+pz2[2]*dx*dx;
    //
    //
    Double_t x[10];
    x[0]=isec;
    x[1]=irow;
    x[2]=c->GetY();
    x[3]=c->GetZ();
    x[4]=py1[1];
    x[5]=pz1[1];
    x[6]=py2[2]*150*150/4; // sagita 150 cm
    //
    x[7]=c->GetY()-y1;
    x[8]=c->GetZ()-z1;
    fDiffHistoLine->Fill(x);
    x[7]=c->GetY()-y2;
    x[8]=c->GetZ()-z2;
    fDiffHistoPar->Fill(x);   

    if (TMath::Abs(py2[2]*150*150/4)<1 && TMath::Abs(pz2[2]*150*150/4)<1){
      // Apply sagita cut
      AddPoint(isec,irow,c->GetZ()-z1, c->GetY()-y1,
	       py1[1],pz1[1],1-TMath::Abs(c->GetZ()/250.),1);
    }
	     
    if (fStreamLevel>1){
      TTreeSRedirector *cstream = GetDebugStreamer();
      if (cstream){
	(*cstream)<<"Diff"<<
	  "run="<<fRun<<              //  run number
	  "event="<<fEvent<<          //  event number
	  "time="<<fTime<<            //  time stamp of event
	  "trigger="<<fTrigger<<      //  trigger
	  "mag="<<fMagF<<             //  magnetic field	  
	  "isec="<<isec<<
	  "Cl.="<<c<<
	  "y1="<<y1<<
	  "y2="<<y2<<
	  "z1="<<z1<<
	  "z2="<<z2<<
	  "py1.="<<&py1<<
	  "py2.="<<&py2<<
	  "pz1.="<<&pz1<<
	  "pz2.="<<&pz2<<
	  "\n";
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
  xmin[1]=0; xmax[1]=159; nbins[1]=159;
  //
  axisName[2]="ly";
  xmin[2]=-50; xmax[2]=50; nbins[2]=10;
  //
  axisName[3]="lz";
  xmin[3]=-250; xmax[3]=250; nbins[3]=50;
  //
  axisName[4]="p2";
  xmin[4]=-0.6; xmax[4]=0.6; nbins[4]=12;
  //
  axisName[5]="p3";
  xmin[5]=-0.6; xmax[5]=0.6; nbins[5]=12;
  //
  axisName[6]="p4";
  xmin[6]=-2; xmax[6]=2; nbins[6]=20;
  //
  axisName[7]="dy";
  xmin[7]=-0.5; xmax[7]=0.5; nbins[7]=100;
  //
  axisName[8]="dz";
  xmin[8]=-0.5; xmax[8]=0.5; nbins[8]=100;
  //
  fDiffHistoLine = new THnSparseF("DistHistoLine","DistHistoLine",9,nbins,xmin,xmax);
  fDiffHistoPar = new THnSparseF("DistHistoPar","DistHistoPar",9,nbins,xmin,xmax);
  for (Int_t i=0; i<9;i++){
    fDiffHistoLine->GetAxis(i)->SetName(axisName[i].Data());
    fDiffHistoPar->GetAxis(i)->SetName(axisName[i].Data());
    fDiffHistoLine->GetAxis(i)->SetTitle(axisName[i].Data());
    fDiffHistoPar->GetAxis(i)->SetTitle(axisName[i].Data());
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


void AliTPCcalibUnlinearity::DumpTree(){
  //
  //
  // 
  if (fStreamLevel==0) return;
  TTreeSRedirector *cstream = GetDebugStreamer();
  if (!cstream) return;
  //  
  THnSparse *his=0;  
  Double_t position[10];
  Double_t value; 
  Int_t *bins = new Int_t[10];
  //
  for (Int_t ihis=0; ihis<2; ihis++){
    his =  (ihis==0)? fDiffHistoLine:fDiffHistoPar;
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
}


void AliTPCcalibUnlinearity::AddPoint(Int_t sector, Int_t row, Double_t dz, Double_t dy, Double_t p2, Double_t p3, Double_t dr, Int_t npoints){
  //
  //
  // Simple distortion fit in outer sectors
  //
  // Model - 2 exponential decay toward the center of chamber
  //       - Decay length 10 and 5 cm
  //       - Decay amplitude - Z dependent 
  //
  /*
  chainUnlin->SetAlias("side","(-1+((sector%36)<18)*2)");
  chainUnlin->SetAlias("sa","sin(pi*sector*20/180)");
  chainUnlin->SetAlias("ca","cos(pi*sector*20/180)");
  chainUnlin->SetAlias("dzt","dz*side");
  chainUnlin->SetAlias("dr","(1-abs(lz/250))");
  chainUnlin->SetAlias("dout","(159-row)*1.5");
  chainUnlin->SetAlias("din","row*0.75");
  chainUnlin->SetAlias("eout10","exp(-(dout)/10.)");
  chainUnlin->SetAlias("eout5","exp(-(dout)/5.)");  
  */

  Double_t side   = (-1+((sector%36)<18)*2); // A side +1   C side -1
  Double_t dzt    = dz*side;
  Double_t dout   = (159-row)*1.5;  // distance to the last pad row - (valid only for long pads)
  if (dout>30) return; // use only edge pad rows
  dr-=0.5;             // dr shifted to the middle - reduce numerical instabilities

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
  xxxR[0]=p2*eoutp;                //fstring+="eoutp++";  
  xxxR[1]=p2*eoutm;               //fstring+="eoutm++";  
  xxxR[2]=p2*dr*eoutp;            //fstring+="dr*eoutp++";  
  xxxR[3]=p2*dr*eoutm;            //fstring+="dr*eoutm++";  
  xxxR[4]=p2*dr*dr*eoutp;         //fstring+="dr*dr*eoutp++";  
  xxxR[5]=p2*dr*dr*eoutm;         //fstring+="dr*dr*eoutm++";    
  //
  xxxRZ[0]=p3*eoutp;                //fstring+="eoutp++";  
  xxxRZ[1]=p3*eoutm;               //fstring+="eoutm++";  
  xxxRZ[2]=p3*dr*eoutp;            //fstring+="dr*eoutp++";  
  xxxRZ[3]=p3*dr*eoutm;            //fstring+="dr*eoutm++";  
  xxxRZ[4]=p3*dr*dr*eoutp;         //fstring+="dr*dr*eoutp++";  
  xxxRZ[5]=p3*dr*dr*eoutm;         //fstring+="dr*dr*eoutm++";    
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


void AliTPCcalibUnlinearity::MakeFitters(){
  //
  //
  //
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

Double_t AliTPCcalibUnlinearity::GetDr(Int_t sector, Double_t dout, Double_t dr){
  //
  //
  //
  TVectorD * vec = GetParamOutR(sector);
  if (!vec) return 0;
  dr-=0.5;        // dr=0 at the middle drift length
  Double_t eout10 = TMath::Exp(-dout/10.);
  Double_t eout5 = TMath::Exp(-dout/5.);		    
  Double_t eoutp  = (eout10+eout5)*0.5;
  Double_t eoutm  = (eout10-eout5)*0.5;

  Double_t result=0;
  result+=(*vec)[1]*eoutp;
  result+=(*vec)[2]*eoutm;
  result+=(*vec)[3]*eoutp*dr;
  result+=(*vec)[4]*eoutm*dr;
  result+=(*vec)[5]*eoutp*dr*dr;
  result+=(*vec)[6]*eoutm*dr*dr;
  return result;
}


Double_t AliTPCcalibUnlinearity::GetGDr(Int_t stype,Float_t gx, Float_t gy,Float_t gz){
  //
  //
  //
  Double_t alpha    = TMath::ATan2(gy,gx);
  if (alpha<0)  alpha+=TMath::Pi()*2;
  Int_t lsector = Int_t(18*alpha/(TMath::Pi()*2));
  alpha = (lsector+0.5)*TMath::Pi()/9.;
  //
  Double_t r[3];
  r[0]=gx; r[1]=gy;
  Double_t cs=TMath::Cos(-alpha), sn=TMath::Sin(-alpha), x=r[0];
  r[0]=x*cs - r[1]*sn; r[1]=x*sn + r[1]*cs;
  //
  AliTPCROC * roc = AliTPCROC::Instance();
  Double_t xout = roc->GetPadRowRadiiUp(roc->GetNRows(37)-1);
  Double_t dout = xout-r[0];
  if (dout<0) return 0;
  Double_t dr   = 1-TMath::Abs(gz/250.);
  if (gz<0) lsector+=18;
  if (stype==0) lsector = (gz>0) ? 36:37;
  if (stype<0) return lsector;  // 
  return GetDr(lsector,dout,dr);
}


Double_t AliTPCcalibUnlinearity::GetDz(Int_t sector, Double_t dout, Double_t dr){
  //
  //
  //
  TVectorD * vec = GetParamOutZ(sector);
  if (!vec) return 0;
  dr-=0.5; // 0 at the middle
  Double_t eout10 = TMath::Exp(-dout/10.);
  Double_t eout5 = TMath::Exp(-dout/5.);
  Double_t eoutp  = (eout10+eout5)*0.5;
  Double_t eoutm  = (eout10-eout5)*0.5;

  
  Double_t result=0;
  result+=(*vec)[1]*eoutp;
  result+=(*vec)[2]*eoutm;
  result+=(*vec)[3]*eoutp*dr;
  result+=(*vec)[4]*eoutm*dr;
  result+=(*vec)[5]*eoutp*dr*dr;
  result+=(*vec)[6]*eoutm*dr*dr;
  return result;
}


Double_t      AliTPCcalibUnlinearity::SGetDr(Int_t sector, Double_t dout, Double_t dr){
  //
  //
  //
  return fgInstance->GetDr(sector,dout,dr);
}
Double_t      AliTPCcalibUnlinearity::SGetGDr(Int_t stype,Float_t gx, Float_t gy,Float_t gz){
  //
  //
  //
  return fgInstance->GetGDr(stype,gx,gy,gz);
}
Double_t      AliTPCcalibUnlinearity::SGetDz(Int_t sector, Double_t dout, Double_t dr){
  //
  //
  //
  return fgInstance->GetDz(sector,dout,dr);
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
    if (cl->GetType()<0) continue; 
    Double_t dy = cl->GetY()-ty;
    Double_t dz = cl->GetZ()-tz;
    Double_t dr = 1-TMath::Abs(cl->GetZ()/250);
    //
    //
    if (TMath::Abs(dy)>0.25) continue;
    if (TMath::Abs(dz)>0.25) continue;
    
    if (TMath::Abs((*vy2)[2]*150*150/4)>1) continue;
    if (TMath::Abs((*vz2)[2]*150*150/4)>1) continue;
      // Apply sagita cut
    calib->AddPoint(cl->GetDetector(), row, dz, dy,
		    (*vy)[1],(*vz)[1], dr, 1);
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







