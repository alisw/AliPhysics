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
  Calculation of the Electric field:
  Sollution of laplace equation in cartezian system, with boundary condition.
  Se details:
  http://web.mit.edu/6.013_book/www/chapter5/5.10.html


*/

/* $Id: AliTPCEfield.cxx 41275 2010-05-16 22:23:06Z marian $ */

#include "TTreeStream.h"
#include "TMath.h"
#include "TLinearFitter.h"
#include "TRandom.h"
#include "AliTPCEfield.h"

ClassImp(AliTPCEfield)


AliTPCEfield* AliTPCEfield::fgInstance=0;

AliTPCEfield::AliTPCEfield():
  TNamed(),  
  fScale(0),
  fMaxFreq(0),
  fIs2D(kTRUE),
  fWorkspace(0)   // file with trees, pictures ...
{
  //
  for (Int_t i=0; i<3; i++){
    fMin[i]=0; fMax[i]=0;
  }
  fgInstance =this;
}

AliTPCEfield::AliTPCEfield(const char* name, Int_t maxFreq, Bool_t is2D, Bool_t useLinear):
  TNamed(name,name),  
  fScale(0),
  fMaxFreq(maxFreq),
  fIs2D(is2D),
  fUseLinear(useLinear),
  fWorkspace(0)   // file with trees, pictures ...
{
  //
  for (Int_t i=0; i<3; i++){
    fMin[i]=0; fMax[i]=0;
  }
  fWorkspace=new TTreeSRedirector(Form("%s.root",name));
  MakeFitFunctions(maxFreq);
  fgInstance =this;
}

void AliTPCEfield::MakeFitFunctions(Int_t maxFreq){
  //
  // fit functions = f(x,y,z) = fx(x)*fy(y)*fz(z)
  // function can be of following types: 
  // 0 - constant
  // 1 - linear
  // 2 - hypx
  // 3 - hypy
  // 4 - hypz
  //
  Int_t nfunctions=0;
  if (fIs2D)     nfunctions = 1+(maxFreq)*8;
  if (!fIs2D)    nfunctions = 1+(maxFreq)*8;
  if (fUseLinear) nfunctions+=2;
  fFitFunctions  = new TMatrixD(nfunctions,4);
  fFitParam      = new TVectorD(nfunctions);
  fFitCovar      = new TMatrixD(nfunctions,nfunctions);
  TMatrixD &fitF = *fFitFunctions;
  // constant function
  Int_t counter=1;
  fitF(0,0)=0;
  //
  // linear functions in one cordinate - constan in others
  if (fUseLinear) 
    for (Int_t ifx=0; ifx<=1; ifx++)
    for (Int_t ify=0; ify<=1; ify++)
      for (Int_t ifz=0; ifz<=1; ifz++){
	if (fIs2D && ifz>0) continue;
	if ((ifx+ify+ifz)==0) continue;
	if ((ifx+ify+ifz)>1) continue;
	fitF(counter,0)= 1;         // function type 
	fitF(counter,1)= ifx;       // type of x function
	fitF(counter,2)= ify;       // type of y function
	fitF(counter,3)= ifz;       // type of z function
	counter++;
      }
  //
  if (fIs2D){
    for (Int_t ihyp=0; ihyp<2; ihyp++)      
	for (Int_t ifx=-maxFreq; ifx<=maxFreq; ifx++)
	  for (Int_t ify=-maxFreq; ify<=maxFreq; ify++){
	    if (TMath::Abs(ify)!=TMath::Abs(ifx)) continue;	    
	    if (ifx==0) continue;
	    if (ify==0) continue;
	    fitF(counter,0)= 2+ihyp;    // function type 
	    fitF(counter,1)= ifx;       // type of x function  - + sinus - cosin
	    fitF(counter,2)= ify;       // type of y function
	    fitF(counter,3)= 0;         // type of y function
	    counter++;
	  }
  }

}



AliTPCEfield::~AliTPCEfield() {
  //
  // Destructor
  //
  if (fWorkspace) delete fWorkspace;
}

void AliTPCEfield::SetRange(Double_t x0, Double_t x1, Double_t y0, Double_t y1, Double_t z0,Double_t z1){
  //
  // Set the ranges - coordinates are rescaled in order to use proper
  // cos,sin expansion in scaled space
  //
  fMin[0]=x0; fMax[0]=x1;
  fMin[1]=y0; fMax[1]=y1;
  fMin[2]=z0; fMax[2]=z1;
  if (fIs2D) fScale=0.5*(TMath::Abs(x1-x0)+TMath::Abs(y1-y0));
  if (!fIs2D) fScale=0.5*(TMath::Abs(x1-x0)+TMath::Abs(y1-y0)+TMath::Abs(z1-z0));
}


void AliTPCEfield::AddBoundaryLine(Double_t x0,Double_t y0,Double_t z0,  Double_t v0, Double_t x1, Double_t y1, Double_t z1,Double_t v1, Int_t id, Int_t npoints){
  //
  // Add a e field boundary line
  // From point (x0,y0) to point (x1,y1)
  // Linear decrease of potential is assumed
  // Boundary can be identified using boundary ID
  // The line is written into tree Boundary
  // 
  Double_t deltaX = (x1-x0);
  Double_t deltaY = (y1-y0);
  Double_t deltaZ = (z1-z0);
  Double_t deltaV = (v1-v0);  
  for (Int_t ipoint=0; ipoint<npoints; ipoint++){
    Double_t bpoint=gRandom->Rndm();
    Double_t x = x0+deltaX*bpoint;
    Double_t y = y0+deltaY*bpoint;
    Double_t z = z0+deltaZ*bpoint;
    Double_t v = v0+deltaV*bpoint;    
    (*fWorkspace)<<"Boundary"<<
      "x="<<x<<       // x coordinate
      "y="<<y<<       // y coordinate
      "z="<<z<<       // z coordinate
      "v="<<v<<       // potential
      "id="<<id<<     // boundary ID
      "\n";    
  }
}

TTree * AliTPCEfield::GetTree(const char * tname){
  //  
  //
  //
  return ((*fWorkspace)<<tname).GetTree();
}

Double_t AliTPCEfield::Field(Int_t ftype,  Double_t ifx, Double_t ify, Double_t ifz, Double_t x, Double_t y, Double_t z){
  //
  // Field component in 
  // f frequency
  Double_t fx=1,fy=1,fz=1;
  const Double_t kEps=0.01;
  //
  if (ftype==0) return 1;
  if (ftype==1) {
    if (TMath::Nint(ifx)==1) return x;
    if (TMath::Nint(ify)==1) return y;
    if (TMath::Nint(ifz)==1) return z;
  }
  Double_t pi = TMath::Pi();
  if (ifx>kEps)  fx = (ftype==2) ? SinHNorm(ifx*pi*x,TMath::Abs(ifx*pi)): TMath::Sin(ifx*pi*x);
  if (ifx<-kEps) fx = (ftype==2) ? CosHNorm(ifx*pi*x,TMath::Abs(ifx*pi)): TMath::Cos(ifx*pi*x);
  //
  if (ify>kEps)  fy = (ftype==3) ? SinHNorm(ify*pi*y,TMath::Abs(ify*pi)): TMath::Sin(ify*pi*y);
  if (ify<-kEps) fy = (ftype==3) ? CosHNorm(ify*pi*y,TMath::Abs(ify*pi)): TMath::Cos(ify*pi*y);
  //
  if (ifz>kEps)  fz = (ftype==4) ? SinHNorm(ifz*pi*z,TMath::Abs(ifz*pi)): TMath::Sin(ifz*pi*z);
  if (ifz<-kEps) fz = (ftype==4) ? CosHNorm(ifz*pi*z,TMath::Abs(ifz*pi)): TMath::Cos(ifz*pi*z);
  (*fWorkspace)<<"eval"<<
    "x="<<x<<
    "y="<<y<<
    "z="<<z<<
    "ifx="<<ifx<<
    "ify="<<ify<<
    "ifz="<<ifz<<
    "fx="<<fx<<
    "fy="<<fy<<
    "fz="<<fz<<
    "ftype="<<ftype<<
    "\n";
  return fx*fy*fz;
}


Double_t AliTPCEfield::FieldDn(Int_t ftype, Double_t ifx, Double_t ify, Double_t ifz, Int_t dn, Double_t x, Double_t y, Double_t z){
  //
  //
  //
 //
  // Field component in 
  // f frequency
  Double_t fx=1,fy=1,fz=1;
  const Double_t kEps=0.01;
  //
  if (ftype==0) return 0.;
  if (ftype==1) {
    Double_t value=0;
    if (TMath::Nint(ifx)==1 &&dn==0) value=1.;
    if (TMath::Nint(ify)==1 &&dn==1) value=1.;
    if (TMath::Nint(ifz)==1 &&dn==2) value=1.;
    return value;    
  }
  Double_t pi = TMath::Pi();
  if (ifx>kEps)  fx = (ftype==2) ? SinHNorm(ifx*pi*x,TMath::Abs(ifx*pi)): TMath::Sin(ifx*pi*x);
  if (ifx<-kEps) fx = (ftype==2) ? CosHNorm(ifx*pi*x,TMath::Abs(ifx*pi)): TMath::Cos(ifx*pi*x);
  //
  if (ify>kEps)  fy = (ftype==3) ? SinHNorm(ify*pi*y,TMath::Abs(ify*pi)): TMath::Sin(ify*pi*y);
  if (ify<-kEps) fy = (ftype==3) ? CosHNorm(ify*pi*y,TMath::Abs(ify*pi)): TMath::Cos(ify*pi*y);
  //
  if (ifz>kEps)  fz = (ftype==4) ? SinHNorm(ifz*pi*z,TMath::Abs(ifz*pi)): TMath::Sin(ifz*pi*z);
  if (ifz<-kEps) fz = (ftype==4) ? CosHNorm(ifz*pi*z,-TMath::Abs(ifz*pi)): TMath::Cos(ifz*pi*z);
  
  if (dn==0){
    if (ifx>kEps)  fx = (ftype==2) ? CosHNorm(ifx*pi*x,TMath::Abs(ifx*pi)): TMath::Cos(ifx*pi*x);
    if (ifx<-kEps) fx = (ftype==2) ? SinHNorm(ifx*pi*x,TMath::Abs(ifx*pi)): -TMath::Sin(ifx*pi*x);
    fx*=ifx*pi;
  }
  if (dn==1){
    if (ify>kEps)  fy = (ftype==3) ? CosHNorm(ify*pi*y,TMath::Abs(ify*pi)): TMath::Cos(ify*pi*y);
    if (ify<-kEps) fy = (ftype==3) ? SinHNorm(ify*pi*y,TMath::Abs(ify*pi)): -TMath::Sin(ify*pi*y);
    fy*=ify*pi;
  }
  if (dn==2){
    if (ifz>kEps)  fz = (ftype==4) ? CosHNorm(ifz*pi*z,TMath::Abs(ifz*pi)): TMath::Cos(ifz*pi*z);
    if (ifz<-kEps) fz = (ftype==4) ? SinHNorm(ifz*pi*z,TMath::Abs(ifz*pi)): -TMath::Sin(ifz*pi*z);
    fz*=ifz*pi;
  }

  return fx*fy*fz;
}




Double_t AliTPCEfield::EvalField(Int_t ifun, Double_t x, Double_t y, Double_t z, Int_t type){
  //
  // Evaluate function ifun at position gx amd gy
  // type == 0 - field
  //      == 1 - Ex
  //      == 2 - Ey
  //      == 3 - Ez
  TMatrixD &mat    = *fFitFunctions;
  Int_t     fid   = TMath::Nint(mat(ifun,0));
  Double_t   ifx   = (mat(ifun,1));
  Double_t   ify   = (mat(ifun,2));
  Double_t   ifz   = (mat(ifun,3));
  //
  if (type==0) return Field(fid,ifx,ify,ifz, x, y,z);
  if (type>0)  return FieldDn(fid,ifx,ify,ifz,type-1, x, y,z);
  return 0;
}

Double_t AliTPCEfield::Eval(Double_t x, Double_t y, Double_t z, Int_t type){
  //
  // Evaluate function ifun at position gx amd gy
  // type == 0 - field
  //      == 1 - Ex
  //      == 2 - Ey
  //      == 3 - Ez
  Double_t value=0;   
  Double_t lx= 2.*(x-(fMin[0]+fMax[0])*0.5)/fScale;
  Double_t ly= 2.*(y-(fMin[1]+fMax[1])*0.5)/fScale;
  Double_t lz= 2.*(z-(fMin[2]+fMax[2])*0.5)/fScale;
  //
  Int_t nfun=fFitFunctions->GetNrows();
  for (Int_t ifun=0; ifun<nfun; ifun++){
    if (type==0) value+=(*fFitParam)[ifun]*EvalField(ifun,lx,ly,lz,type);
    if (type>0)  value+=2*(*fFitParam)[ifun]*EvalField(ifun,lx,ly,lz,type)/fScale;
  }
  return value;
}

Double_t AliTPCEfield::EvalS(Double_t x, Double_t y, Double_t z,  Int_t type){
  //
  // static evaluation - possible to use it in the TF1 
  //
  return fgInstance->Eval(x,y,z,type);
}

void AliTPCEfield::FitField(){
  //
  // Fit the e field
  // Minimize chi2 residuals at the boundary points 
  // ?Tempoary sollution - integrals can be calculated analytically -
  //
  Int_t nfun=fFitFunctions->GetNrows();
  Double_t *fun =new Double_t[nfun];
  fFitter= new TLinearFitter(nfun, Form("hyp%d", nfun-1));
  //
  TTree * tree = GetTree("Boundary");
  Int_t npoints = tree->GetEntries();
  Int_t   *indexes = new Int_t[npoints]; 
  Double_t *rindex  = new Double_t[npoints]; 
  //

  Double_t x=0, y=0, z=0, v=0;
  tree->SetBranchAddress("x",&x);
  tree->SetBranchAddress("y",&y);
  tree->SetBranchAddress("z",&z);
  tree->SetBranchAddress("v",&v);
  TMatrixD valMatrix(npoints,4);
  for (Int_t ipoint=0; ipoint<npoints; ipoint++){
    tree->GetEntry(ipoint);
    valMatrix(ipoint,0)=2.*(x-(fMin[0]+fMax[0])*0.5)/fScale;
    valMatrix(ipoint,1)=2.*(y-(fMin[1]+fMax[1])*0.5)/fScale;
    valMatrix(ipoint,2)=2.*(z-(fMin[2]+fMax[2])*0.5)/fScale;
    valMatrix(ipoint,3)=v;    
    rindex[ipoint]=gRandom->Rndm();
  }
  TMath::Sort(npoints,rindex,indexes);
  for (Int_t jpoint=0; jpoint<npoints; jpoint++){
    Int_t ipoint=indexes[jpoint];
    //
    Double_t lx= valMatrix(ipoint,0);
    Double_t ly= valMatrix(ipoint,1);
    Double_t lz= valMatrix(ipoint,2);
    Double_t value = valMatrix(ipoint,3);
    for (Int_t ifun=1; ifun<nfun; ifun++){      
      Double_t ffun=EvalField(ifun,lx,ly,lz,0);
      fun[ifun-1]=ffun;	
    }
    fFitter->AddPoint(fun,value,1);    
  }
  fFitter->Eval();
  fFitter->GetCovarianceMatrix(*fFitCovar);
  fFitter->GetParameters(*fFitParam); 

  (*fWorkspace)<<"fitGlobal"<<
    "covar.="<<fFitCovar<<        // covariance matrix
    "param.="<<fFitParam<<        // fit parameters
    "fun.="<<fFitFunctions<<      // type of fit functions
    "\n";

  for (Int_t ifun=1; ifun<nfun; ifun++){
    
  }


}


TMatrixD* AliTPCEfield::MakeCorrelation(TMatrixD &matrix){
  //
  //
  //
  Int_t nrows = matrix.GetNrows();
  TMatrixD * mat = new TMatrixD(nrows,nrows);
  for (Int_t irow=0; irow<nrows; irow++)
    for (Int_t icol=0; icol<nrows; icol++){
      (*mat)(irow,icol)= matrix(irow,icol)/TMath::Sqrt(matrix(irow,irow)*matrix(icol,icol));
    }
  return mat;
}




void AliTPCEfield::DumpField(Double_t gridSize, Double_t step){
  //
  //
  //
  Double_t stepSize=0.001*fScale/fMaxFreq;
  //
  for (Double_t x = fMin[0]+stepSize; x<=fMax[0]-stepSize; x+=gridSize){
    for (Double_t y = fMin[1]+stepSize; y<=fMax[1]-stepSize; y+=gridSize)
      for (Double_t z = fMin[2]; z<=fMax[2]; z+=gridSize){
	//
	//
	Double_t v  =  Eval(x,y,z,0);
	Double_t ex =  Eval(x,y,z,1);
	Double_t ey =  Eval(x,y,z,2);
	Double_t ez =  Eval(x,y,z,3);
	Double_t dexdx =  (Eval(x,y,z,1)-Eval(x-stepSize,y,z,1))/stepSize;  // numerical derivative
	Double_t deydy =  (Eval(x,y,z,2)-Eval(x,y-stepSize,z,2))/stepSize;
	Double_t dezdz =  (Eval(x,y,z,3)-Eval(x,y,z-stepSize,3))/stepSize;
	(*fWorkspace)<<"dumpField"<<
	  "x="<<x<<  // position
	  "y="<<y<<
	  "z="<<z<<
	  "v="<<v<<  // potential
	  "ex="<<ex<<  // Efield
	  "ey="<<ey<<
	  "ez="<<ez<<
	  "dexdx="<<dexdx<<  //gradient of e field
	  "deydy="<<deydy<<
	  "dezdz="<<dezdz<<
	  "\n";
      }
  }
  //
  //
  for (Double_t x = fMin[0]+stepSize; x<=fMax[0]-stepSize; x+=gridSize){
    Double_t sumEx=0;
    Double_t sumEy=0;
    Double_t sumExEy=0;
    for (Double_t y = fMin[1]+0.2; y<=fMax[1]-stepSize; y+=step)
      for (Double_t z = fMin[2]; z<=fMax[2]; z+=gridSize){
	//
	//
	Double_t v  =  Eval(x,y+step*0.5,z,0);
	Double_t ex =  Eval(x,y+step*0.5,z,1);
	Double_t ey =  Eval(x,y+step*0.5,z,2);
	Double_t ez =  Eval(x,y+step*0.5,z,3);
	sumEx+=ex*step;
	sumEy+=ey*step;
	sumExEy+=(ex/ey)*step;
	(*fWorkspace)<<"dumpDistortion"<<
	  "x="<<x<<  // position
	  "y="<<y<<
	  "z="<<z<<
	  "v="<<v<<  // potential
	  "ex="<<ex<<  // Efield
	  "ey="<<ey<<
	  "ez="<<ez<<
	  //                 
	  "sEx="<<sumEx<<       // field x integral
	  "sEy="<<sumEy<<       // field y integral
	  "sExEy="<<sumExEy<<   // tan integral
	  "\n";
      }
  }
}


void MakeTPC2DExample(AliTPCEfield *field){
  //
  /*
    .L  $ALICE_ROOT/TPC/AliTPCEfield.cxx++
    AliTPCEfield *field =  new AliTPCEfield("field",20, kTRUE,kTRUE);
    MakeTPC2DExample(field)
    field->FitField()
    sqrt(field->fFitter.GetChisquare()/field->fFitter.GetNpoints())

    TF2 f2("f2","AliTPCEfield::EvalS(x,y,0,0)",90,245,0,250);
    f2->SetNpx(100);     f2->SetNpy(100); f2->Draw("colz");

    TF2 f2x("f2x","AliTPCEfield::EvalS(x,y,0,1)",90,240,0,240);
    f2x->SetNpx(100);     f2x->SetNpy(100);  f2x->Draw("surf2");

    TF2 f2y("f2y","AliTPCEfield::EvalS(x,y,0,2)",90,240,0,240);
    f2y->SetNpx(100);     f2y->SetNpy(100);  f2y->Draw("surf2");

    field->MakeCorrelation(*(field->fFitCovar)).Print()
    Double_t index[100000];
    for (Int_t i=0; i<field->fFitCovar->GetNrows(); i++) index[i]=i;
    TGraph gr(field->fFitCovar->GetNrows(), index, field->fFitParam->GetMatrixArray());
    gr->Draw("alp");

    field->GetTree()->Draw("AliTPCEfield::EvalS(x,y,0,0):v");


    TF2 f2xdy("f2xdy","AliTPCEfield::EvalS(x,y,0,1)/AliTPCEfield::EvalS(x,y,0,2)",90,240,0,240);
    f2xdy->SetNpx(100);     f2xdy->SetNpy(100);  f2xdy->Draw("colz");

  */

  Double_t p0[4];
  Double_t p1[4];
  Double_t xmin=85, xmax=245;
  Double_t ymin=0, ymax=250, deltaY=0.1*ymax;
  Double_t vup=1;
  Int_t npoints=1000;
  field->SetRange(xmin, xmax,ymin,ymax,0,0);
  // upper part
  p0[0]=xmin; p0[1]=ymax+deltaY; p0[2]=0; p0[3]=vup;
  p1[0]=xmax; p1[1]=ymax-deltaY; p1[2]=0; p1[3]=vup;
  field->AddBoundaryLine(p0[0],p0[1], p0[2],p0[3],p1[0],p1[1], p1[2],p1[3],1,npoints);
  //left
  p0[0]=xmin; p0[1]=ymin+deltaY; p0[2]=0; p0[3]=0;
  p1[0]=xmin; p1[1]=ymax+deltaY; p1[2]=0; p1[3]=vup;
  field->AddBoundaryLine(p0[0],p0[1], p0[2],p0[3],p1[0],p1[1], p1[2],p1[3],2,npoints);
  //right
  p0[0]=xmax; p0[1]=ymin-deltaY; p0[2]=0; p0[3]=0;
  p1[0]=xmax; p1[1]=ymax-deltaY; p1[2]=0; p1[3]=vup;
  field->AddBoundaryLine(p0[0],p0[1], p0[2],p0[3],p1[0],p1[1], p1[2],p1[3],3,npoints);
  //
  //ROC
  p0[0]=xmin; p0[1]=deltaY; p0[2]=0; p0[3]=-0;
  p1[0]=xmax; p1[1]=-deltaY; p1[2]=0; p1[3]=0;
  field->AddBoundaryLine(p0[0],p0[1], p0[2],p0[3],p1[0],p1[1], p1[2],p1[3],4,npoints);		 
}



