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
  Only Linear and parabolic fit used
  To be used for laser or for data without field
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

#include "AliTPCclusterMI.h"
#include "AliTPCseed.h"
#include "AliESDVertex.h"
#include "AliESDEvent.h"
#include "AliESDfriend.h"
#include "AliESDInputHandler.h"
#include "AliAnalysisManager.h"
#include "AliTracker.h"
#include "AliMagFMaps.h"
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


#include "AliTPCcalibUnlinearity.h"


ClassImp(AliTPCcalibUnlinearity)

AliTPCcalibUnlinearity::AliTPCcalibUnlinearity():
  AliTPCcalibBase(),
  fDiffHistoLine(0),
  fDiffHistoPar(0),
  fFittersOutR(38),
  fFittersOutZ(38)
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
  fFittersOutZ(38)
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


void AliTPCcalibUnlinearity::Process(AliTPCseed *seed) {
  //
  // 
  //  
  const Int_t  kdrow=10;
  const Int_t kMinCluster=35;
  if (!fDiffHistoLine) MakeHisto();    
  if (!seed) return;
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
    Float_t dx = c->GetX()-kXmean;
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
    Float_t dx = c->GetX()-kXmean;
    Float_t y1 = py1[0]+py1[1]*dx;
    Float_t y2 = py2[0]+py2[1]*dx+py2[2]*dx*dx;
    //
    Float_t z1 = pz1[0]+pz1[1]*dx;
    Float_t z2 = pz2[0]+pz2[1]*dx+pz2[2]*dx*dx;
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

    if (TMath::Abs(py2[2]*150*150/4)<1 && TMath::Abs(py2[2]*150*150/4)<1){
      // Apply sagita cut
      AddPoint(isec,irow,c->GetZ()-z1, c->GetY()-y1,
	       py1[1],pz1[1],1-TMath::Abs(c->GetZ()/250.),1);
    }
	     
    if (fStreamLevel>1){
      TTreeSRedirector *cstream = GetDebugStreamer();
      if (cstream){
	(*cstream)<<"Diff"<<
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


void AliTPCcalibUnlinearity::AddPoint(Int_t sector, Int_t row, Float_t dz, Float_t dy, Float_t p2, Float_t p3, Float_t dr, Int_t npoints){
  //
  //
  // Simple distortion fit in outer sectors
  //
  // Model - 2 exponential decay toward the center of chamber
  //       - Decay length 10 and 20 cm
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
  chainUnlin->SetAlias("eout20","exp(-(dout)/20.)");  
  */
  Float_t side   = (-1+((sector%36)<18)*2); // A side +1   C side -1
  Float_t dzt    = dz*side;
  Float_t dout   = (159-row)*1.5;  // distance to the last pad row - (valid only for long pads)
  Float_t eout10 = TMath::Exp(-dout/10.);
  Float_t eout20 = TMath::Exp(-dout/20.);
  //
  Double_t xxxR[6], xxxZ[6], xxxRZ[6];
  //TString fstring="";
  xxxZ[0]=eout20;                //fstring+="eout20++";  
  xxxZ[1]=eout10;               //fstring+="eout10++";  
  xxxZ[1]=dr*eout20;            //fstring+="dr*eout20++";  
  xxxZ[1]=dr*eout10;            //fstring+="dr*eout10++";  
  xxxZ[1]=dr*dr*eout20;         //fstring+="dr*dr*eout20++";  
  xxxZ[1]=dr*dr*eout10;         //fstring+="dr*dr*eout10++";    
  //
  xxxR[0]=p2*eout20;                //fstring+="eout20++";  
  xxxR[1]=p2*eout10;               //fstring+="eout10++";  
  xxxR[1]=p2*dr*eout20;            //fstring+="dr*eout20++";  
  xxxR[1]=p2*dr*eout10;            //fstring+="dr*eout10++";  
  xxxR[1]=p2*dr*dr*eout20;         //fstring+="dr*dr*eout20++";  
  xxxR[1]=p2*dr*dr*eout10;         //fstring+="dr*dr*eout10++";    
  //
  xxxRZ[0]=p3*eout20;                //fstring+="eout20++";  
  xxxRZ[1]=p3*eout10;               //fstring+="eout10++";  
  xxxRZ[1]=p3*dr*eout20;            //fstring+="dr*eout20++";  
  xxxRZ[1]=p3*dr*eout10;            //fstring+="dr*eout10++";  
  xxxRZ[1]=p3*dr*dr*eout20;         //fstring+="dr*dr*eout20++";  
  xxxRZ[1]=p3*dr*dr*eout10;         //fstring+="dr*dr*eout10++";    
  //
  TLinearFitter * fitter=0;
  for (Int_t ip=0; ip<npoints; ip++){
    //
    fitter = (TLinearFitter*)fFittersOutR.At(sector%36);
    fitter->AddPoint(xxxR,dy);
    fitter->AddPoint(xxxRZ,dz);
    fitter = (TLinearFitter*)fFittersOutZ.At(sector%36);
    fitter->AddPoint(xxxZ,dzt);
    //
    if (side>0){
      fitter = (TLinearFitter*)fFittersOutR.At(36);
      fitter->AddPoint(xxxR,dy);
      fitter->AddPoint(xxxRZ,dz);
      fitter = (TLinearFitter*)fFittersOutZ.At(36);
      fitter->AddPoint(xxxZ,dzt);
    }
    if (side<0){
      fitter = (TLinearFitter*)fFittersOutR.At(36);
      fitter->AddPoint(xxxR,dy);
      fitter->AddPoint(xxxRZ,dz);
      fitter = (TLinearFitter*)fFittersOutZ.At(36);
      fitter->AddPoint(xxxZ,dzt);
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
  for (Int_t ifit=0; ifit<38; ifit++){
    fitter=(TLinearFitter*)fFittersOutR.At(ifit);
    if (fitter) fitter->Eval();
    fitter=(TLinearFitter*)fFittersOutZ.At(ifit);
    if (fitter) fitter->Eval();
  }
}


void AliTPCcalibUnlinearity::ProcessTree(TTree * tree, Int_t nmaxPoints){
  //
  //
  //
  //  TTree * tree = chainUnlinP;
  AliTPCcalibUnlinearity *calib = this;
  //
  AliTPCclusterMI * cl=0;
  Float_t ty,tz;
  TVectorD *vy=0, *vz=0;
  Int_t nentries = TMath::Min(Int_t(tree->GetEntries()), nmaxPoints);
  //
  {
  for (Int_t i=0; i<nentries; i++){
    tree->SetBranchAddress("Cl.",&cl);
    tree->SetBranchAddress("y1",&ty);
    tree->SetBranchAddress("z1",&tz);
    tree->SetBranchAddress("py1.",&vy);
    tree->SetBranchAddress("pz1.",&vz);
    //
    tree->GetEntry(i);
    if (!cl) continue;
    calib->AddPoint(cl->GetDetector(), cl->GetRow(), cl->GetZ()-tz, cl->GetY()-ty,
		    (*vy)(1),(*vz)(1), 1-TMath::Abs(cl->GetZ()/250),1);
  }
  }
}
