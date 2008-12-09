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




*/



#include "AliExternalComparison.h"
#include "AliExternalTrackParam.h"
#include "AliTrackReference.h"
#include "THnSparse.h"
#include "TParticle.h"
#include "TMatrix.h"
#include "TParticlePDG.h"
#include "TParticle.h"
#include "TTreeStream.h"
#include "TAxis.h"

ClassImp(AliExternalComparison)


AliExternalComparison::AliExternalComparison() 
  :TNamed(),
   fResolHistos(0),
   fPullHistos(0),
   fRangeMatrix(0),
   fCutMatrix(0)
{  
  //
  // Default constructor
  //
}


AliExternalComparison::AliExternalComparison(const Text_t *name, const Text_t *title)
  :TNamed(name,title),
   fResolHistos(0),
   fPullHistos(0),
   fRangeMatrix(0),
   fCutMatrix(0)
{
  //
  // Non default cosntructor
  //
}

AliExternalComparison::AliExternalComparison(const AliExternalComparison& comp)
  :TNamed(comp.fName,comp.fTitle),
   fResolHistos(new TObjArray(*(comp.fResolHistos))),
   fPullHistos(new TObjArray(*(comp.fPullHistos))),
   fRangeMatrix(new TMatrixD(*(comp.fRangeMatrix))),
   fCutMatrix(new TMatrixD(*(comp.fCutMatrix)))
{
  //
  // copy constructor
  //
}

AliExternalComparison& AliExternalComparison::operator=(const AliExternalComparison&comp)
{
  //
  //
  //
  SetName(comp.GetName());
  SetTitle(comp.GetTitle());
  fResolHistos=new TObjArray(*(comp.fResolHistos));
  fPullHistos=new TObjArray(*(comp.fPullHistos));
  fRangeMatrix=new TMatrixD(*(comp.fRangeMatrix));
  fCutMatrix=new TMatrixD(*(comp.fCutMatrix));
  return *this;
}


AliExternalComparison::~AliExternalComparison(){
  //
  //
  //
  if (fResolHistos) fResolHistos->Delete();
  if (fPullHistos)  fPullHistos->Delete();
  delete fResolHistos;             // resolution histogram
  delete fPullHistos;              // pull       histogram
  delete fRangeMatrix;            // range matrix
  delete fCutMatrix;            // range matrix
}
 

void AliExternalComparison::Analyze() {
  //
  //
  //
}


Long64_t AliExternalComparison::Merge(TCollection *li) {
  //
  //
  //   
  TIterator* iter = li->MakeIterator();
  AliExternalComparison* comp = 0;
  TString strName(GetName());
  while ((comp = (AliExternalComparison*)iter->Next())) {
    if (!comp->InheritsFrom(AliExternalComparison::Class())) {
      return -1;
    }
    if (strName.CompareTo(comp->GetName())!=0) return -1;
    // add histograms here...
    Add(comp);
  }
  return 0;  
}

void  AliExternalComparison:: Add(AliExternalComparison*comp){
  //
  // Add comparison object
  //
  if (!fResolHistos) return;
  for (Int_t i=0;i<5;i++){
    THnSparse * h0 = (THnSparse*)fResolHistos->At(i);
    THnSparse * h1 = (THnSparse*)comp->fResolHistos->At(i);
    if (h0&&h1) h0->Add(h1);
    h0 = (THnSparse*)fPullHistos->At(i);
    h1 = (THnSparse*)comp->fPullHistos->At(i);
    if (h0&&h1) h0->Add(h1);
  }
}



void   AliExternalComparison::Process(const AliExternalTrackParam *param0, const AliExternalTrackParam *param1){
  //
  // Process- fill histogram with residuals
  //    Tracks has to be in the same local X 
  //    
  if (!AcceptPair(param0,param1)) return;
  //
  if (!fResolHistos) MakeHistos(); 
  //
  //
  const Double_t *p0 = param0->GetParameter();
  const Double_t *p1 = param1->GetParameter();
  const Double_t *c0 = param0->GetCovariance();
  const Double_t *c1 = param1->GetCovariance();
  Double_t xyz[3];
  param1->GetXYZ(xyz);
  //
  Double_t vec[10];
  vec[0] = p1[0];
  vec[1] = p1[1];
  vec[2] = p1[2];
  vec[3] = p1[3];
  vec[4] = TMath::Sqrt(TMath::Abs(p1[4]))* ((p1[4]>0)? 1:-1);
  vec[5] = TMath::Sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]);
  vec[6] = TMath::ATan2(xyz[1],xyz[0]);
  //
  //
  THnSparse * his = 0;
  //
  vec[7] = p0[0]-p1[0];
  his = (THnSparse*)fResolHistos->At(0);
  his->Fill(vec);
  vec[7]=(p0[0]-p1[0])/TMath::Sqrt(c0[0]+c1[0]);
  his = (THnSparse*)fPullHistos->At(0);
  his->Fill(vec);
  //
  //
  vec[7] = p0[1]-p1[1];
  his = (THnSparse*)fResolHistos->At(1);
  his->Fill(vec);
  vec[7]=(p0[1]-p1[1])/TMath::Sqrt(c0[2]+c1[2]);
  his = (THnSparse*)fPullHistos->At(1);
  his->Fill(vec);
  //
  vec[7] = p0[2]-p1[2];
  his = (THnSparse*)fResolHistos->At(2);
  his->Fill(vec);
  vec[7]=(p0[2]-p1[2])/TMath::Sqrt(c0[5]+c1[5]);
  his = (THnSparse*)fPullHistos->At(2);
  his->Fill(vec);
  //
  vec[7] = p0[3]-p1[3];
  his = (THnSparse*)fResolHistos->At(3);
  his->Fill(vec);
  vec[7]=(p0[3]-p1[3])/TMath::Sqrt(c0[9]+c1[9]);
  his = (THnSparse*)fPullHistos->At(3);
  his->Fill(vec);
  //
  vec[7] = p0[4]-p1[4];
  his = (THnSparse*)fResolHistos->At(4);
  his->Fill(vec);
  vec[7]=(p0[4]-p1[4])/TMath::Sqrt(c0[14]+c1[14]);
  his = (THnSparse*)fPullHistos->At(4);
  his->Fill(vec);
}

void  AliExternalComparison::Process(const AliExternalTrackParam *param0, TParticle *part){
  //
  //
  //
  AliExternalTrackParam *param1 = MakeExternalParam(part);
  param1->Rotate(param0->GetAlpha());
  if (param0&&param1) Process(param0, param1);
}
void AliExternalComparison::Process(const AliExternalTrackParam *param0, TParticle *part, const AliTrackReference *ref){
  //
  //
  //
  AliExternalTrackParam *param1 = MakeExternalParam(part,ref);
  param1->Rotate(param0->GetAlpha());
  if (param0&&param1) Process(param0, param1);
}

THnSparse * AliExternalComparison::GetHisto(Int_t ivar, Int_t type) {
  //
  //
  //
  if (!fResolHistos) return 0;
  if (ivar<0 || ivar>4)  return 0;
  if (type==0) return (THnSparse*)fResolHistos->At(ivar);
  if (type==1) return (THnSparse*)fPullHistos->At(ivar);
  return 0;
}

void AliExternalComparison::MakeHistos(){
  //
  //
  // 
  if (!fRangeMatrix) SetDefaultRange();
  TMatrixD & mat = *fRangeMatrix;  
  //
  Double_t xmin[8], xmax[8];
  Double_t xpmin[8], xpmax[8];
  Int_t nbins[10];
  for (Int_t i=0;i<7;i++){
    xmin[i]=mat(i,0);
    xmax[i]=mat(i,1);
    xpmin[i]=mat(i,0);
    xpmax[i]=mat(i,1);
    nbins[i]=TMath::Nint(mat(i,2));
  }
  //
  fResolHistos = new TObjArray(5);
  fPullHistos = new TObjArray(5);
  //
  for (Int_t idelta=0;idelta<5;idelta++){
    xmin[7]=mat(idelta+7,0); 
    xmax[7]=mat(idelta+7,1); 
    nbins[7]=TMath::Nint(mat(idelta+7,2));
    fResolHistos->AddAt(new THnSparseF(Form("Resol%s",GetName()),Form("Resol%s",GetTitle()),
				     8,nbins,xmin,xmax),idelta);
    xmin[7]=-10;
    xmax[7]=10;
    fPullHistos->AddAt(new THnSparseF(Form("Pull%s",GetName()),Form("Pull%s",GetTitle()),
				      8,nbins,xmin,xmax),idelta);
  }
  //
}


void  AliExternalComparison::SetDefaultRange(Float_t scale,Float_t arm, Int_t nbins){
  //
  // set default binning
  // scale - define position range
  // arm  - Level arm
  // Binning:
  // 7 variables
  // 5 deltas

  if (!fRangeMatrix) fRangeMatrix = new TMatrixD(12,3);
  TMatrixD & mat = *fRangeMatrix;  
  //
  // variables
  //
  mat(0,0)=-50;    mat(0,1)=50;   mat(0,2)=50;  // P0 -y            range
  mat(1,0)=-250;   mat(1,1)=250;  mat(1,2)=50;  // P1 -Z            range
  mat(2,0)=-0.99;  mat(2,1)=0.99; mat(2,2)=20;  // P2 -snp          range
  mat(3,0)=-1.5;   mat(3,1)=1.5;  mat(3,2)=20;  // P3 -tantheta     range
  mat(4,0)=-3;     mat(4,1)=3;    mat(4,2)=50;   // sqrt(P4) -sqrt(1/pt) range  
  //
  mat(5,0)= 0;     mat(5,1)=250;  mat(5,2)=50;   // R  range  
  mat(6,0)= -TMath::Pi();   mat(6,1)=TMath::Pi();    mat(6,2)=18*4;   // fi range  
  //
  // Resolution
  //
  mat(7,0)=-scale;      mat(7,1)=scale;   mat(7,2)=nbins;       // P0 -y            range
  mat(8,0)=-scale;      mat(8,1)=scale;   mat(8,2)=nbins;       // P1 -Z            range
  mat(9,0)=-scale/arm;  mat(9,1)=scale/arm; mat(9,2)=nbins;     // P2 -snp          range
  mat(10,0)=-scale/arm; mat(10,1)=scale/arm;  mat(10,2)=nbins;  // P3 -tantheta     range
  mat(11,0)=-1;         mat(11,1)=1;   mat(11,2)=nbins;         // sqrt(P4) -sqrt(1/pt) range
  //
}

void  AliExternalComparison::SetDefaultCuts(){
  //

  if (!fCutMatrix) fCutMatrix = new TMatrixD(10,1);
  TMatrixD & mat = *fCutMatrix;  
  //
  mat(0,0)=10;     //dP0
  mat(1,0)=20;     //dP1
  mat(2,0)=0.05;   //dP2
  mat(3,0)=0.05;   //dP3
  mat(4,0)=1;      //dP4
  //
  mat(5,0)=10;      //dnP0
  mat(6,0)=50;      //dnP1
  mat(7,0)=10;      //dnP2
  mat(8,0)=10;      //dnP3
  mat(9,0)=10;     //dnP4
  //
}

void AliExternalComparison::SetResolRange(Int_t param, Float_t min, Float_t max, Int_t nbins){
  //
  //
  //
  if (!fRangeMatrix) SetDefaultRange();
  TMatrixD & mat = *fRangeMatrix; 
  if (param<0) return;
  if (param>4) return;
  mat(7+param,0)=min;      mat(7+param,1)=max;   mat(7+param,2)=nbins;
}

void AliExternalComparison::SetDistCut(Float_t dP0, Float_t dP1,Float_t dP2,Float_t dP3, Float_t dP4){
  //
  // Set diff cuts
  //
  if (!fCutMatrix) SetDefaultCuts();
  TMatrixD & mat = *fCutMatrix;  
  mat(0,0)=dP0;     //dP0
  mat(1,0)=dP1;     //dP1
  mat(2,0)=dP2;   //dP2
  mat(3,0)=dP3;   //dP3
  mat(4,0)=dP4;      //dP4

  
}
void AliExternalComparison::SetPullDistCut(Float_t dnP0, Float_t dnP1,Float_t dnP2,Float_t dnP3, Float_t dnP4){
  //
  //
  //
  if (!fCutMatrix) SetDefaultCuts();
  TMatrixD & mat = *fCutMatrix;  
  mat(5,0)=dnP0;     //dnP0
  mat(6,0)=dnP1;     //dnP1
  mat(7,0)=dnP2;     //dnP2
  mat(8,0)=dnP3;     //dnP3
  mat(9,0)=dnP4;     //dnP4
}

Bool_t   AliExternalComparison::AcceptPair(const AliExternalTrackParam *param0, const AliExternalTrackParam *param1){  
  //
  //
  //
  Bool_t isOK =kTRUE;
  if (!fCutMatrix) SetDefaultCuts();
  TMatrixD & mat = *fCutMatrix;  
  //
  const Double_t *p0 = param0->GetParameter();
  const Double_t *p1 = param1->GetParameter();
  const Double_t *c0 = param0->GetCovariance();
  const Double_t *c1 = param1->GetCovariance();
  //
  if (TMath::Abs(p0[0]-p1[0])>mat(0,0)) return kFALSE;
  if (TMath::Abs(p0[1]-p1[1])>mat(1,0)) return kFALSE;
  if (TMath::Abs(p0[2]-p1[2])>mat(2,0)) return kFALSE;
  if (TMath::Abs(p0[3]-p1[3])>mat(3,0)) return kFALSE;
  if (TMath::Abs(p0[4]-p1[4])>mat(4,0)) return kFALSE;
  if ((c0[0]+c1[0])<0) return kFALSE;
  if ((c0[2]+c1[2])<0) return kFALSE;
  if ((c0[5]+c1[5])<0) return kFALSE;
  if ((c0[9]+c1[9])<0) return kFALSE;
  if ((c0[14]+c1[14])<0) return kFALSE;
  
  if (TMath::Abs((p0[0]-p1[0])/TMath::Sqrt(c0[0]+c1[0]))>mat(5,0)) return kFALSE;
  if (TMath::Abs((p0[1]-p1[1])/TMath::Sqrt(c0[2]+c1[2]))>mat(6,0)) return kFALSE;
  if (TMath::Abs((p0[2]-p1[2])/TMath::Sqrt(c0[5]+c1[5]))>mat(7,0)) return kFALSE;
  if (TMath::Abs((p0[3]-p1[3])/TMath::Sqrt(c0[9]+c1[9]))>mat(8,0)) return kFALSE;
  if (TMath::Abs((p0[4]-p1[4])/TMath::Sqrt(c0[14]+c1[14]))>mat(9,0)) return kFALSE;
  return isOK;
}

  
AliExternalTrackParam *AliExternalComparison::MakeExternalParam(TParticle *part, const AliTrackReference *ref){
  //
  // make a Kalman track from track reference
  //
  Double_t xyz[3]={ref->X(),ref->Y(),ref->Z()};
  Double_t pxyz[3]={ref->Px(),ref->Py(),ref->Pz()};
  Double_t cv[21];
  for (Int_t i=0; i<21;i++) cv[i]=0;
  if (!part->GetPDG()) return 0;
  AliExternalTrackParam * param = new AliExternalTrackParam(xyz,pxyz,cv,TMath::Nint(part->GetPDG()->Charge()/3.));
  return param;
}


AliExternalTrackParam * AliExternalComparison::MakeExternalParam(TParticle *part){
  //
  //
  //
  Double_t xyz[3]={part->Vx(),part->Vy(),part->Vz()};
  Double_t pxyz[3]={part->Px(),part->Py(),part->Pz()};
  Double_t cv[21];
  for (Int_t i=0; i<21;i++) cv[i]=0;
  if (!part->GetPDG()) return 0;
  AliExternalTrackParam * param = new AliExternalTrackParam(xyz,pxyz,cv,TMath::Nint(part->GetPDG()->Charge()/3.));
  return param;
}




void AliExternalComparison::MakeComparisonTree(const char * outname){
  //
  // make a comparison tree
  // the tree will be written to the file - outname
  //
  /*
    TFile f("comp.root");
    AliExternalComparison *comp =  (AliExternalComparison*)f.Get("compAlign");
    outname="compTree.root";
  */
  AliExternalComparison *comp = this;
  TTreeSRedirector *pcstream = new TTreeSRedirector(outname);
  //
  THnSparse *his=0;  
  Double_t position[10];
  Double_t value; 
  Int_t *bins = new Int_t[10];
  //
  //
  for (Int_t isPull=0; isPull<2; isPull++){
    for (Int_t ires=0;ires<5; ires++){
      his = comp->GetHisto(ires,isPull);
      if (!his) continue;
      //
      Int_t ndim = his->GetNdimensions();
      //
      for (Long64_t i = 0; i < his->GetNbins(); ++i) {
	value = his->GetBinContent(i, bins);
	for (Int_t idim = 0; idim < ndim; idim++) {
	  position[idim] = his->GetAxis(idim)->GetBinCenter(bins[idim]);
	}      
	(*pcstream)<<"Resol"<<
	  "isPull="<<isPull<<     // normalized error indication
	  "index="<<ires<<        // parameter difference index
	  "bincont="<<value<<     // bin content
	  "val="<<position[7]<<   // parameter difference
	  "p0="<<position[0]<<    //P0
	  "p1="<<position[1]<<    //P1
	  "p2="<<position[2]<<    //P2
	  "p3="<<position[3]<<    //P3
	  "p4="<<position[4]<<    //P4
	  "R="<<position[5]<<     //Radius
	  "phi="<<position[6]<<     //Radius
	  "\n";
      }
    }
  }
  delete pcstream;
}
