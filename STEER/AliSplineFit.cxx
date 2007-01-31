/**************************************************************************
 * Copyright(c) 2006-07, ALICE Experiment at CERN, All rights reserved. *
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

//-------------------------------------------------------------------------
//                Implementation of the AliSplineFit class
//   The class performs a spline fit on an incoming TGraph. The graph is 
//   divided into several parts (identified by knots between each part). 
//   Spline fits are performed on each part. According to user parameters,
//   the function, first and second derivative are requested to be continuous 
//   at each knot.
//        Origin: Marian Ivanov, CERN, Marian.Ivanov@cern.ch
//   Adjustments by Haavard Helstrup,  Haavard.Helstrup@cern.ch
//-------------------------------------------------------------------------


#include "AliSplineFit.h"
 
ClassImp(AliSplineFit)

TLinearFitter*
AliSplineFit::fitterStatic()
{
  static TLinearFitter* fit = new TLinearFitter(4,"pol3","");
  return fit;
}

AliSplineFit::AliSplineFit() :
  fBDump(kFALSE),
  fGraph (0),
  fNmin (0),
  fSigma (0),
  fMaxDelta (0),
  fN0 (0),
  fParams (0),
  fCovars (0),
  fIndex (0),
  fN    (0),
  fChi2  (0.0),
  fX   (0),
  fY0  (0),
  fY1  (0),
  fChi2I (0)
  //
  // Default constructor
  //
{ }



AliSplineFit::AliSplineFit(const AliSplineFit& source) :
  TObject(source),
  fBDump (source.fBDump),
  fGraph (source.fGraph),
  fNmin  (source.fNmin),
  fSigma (source.fSigma),
  fMaxDelta (source.fMaxDelta),
  fN0    (source.fN0),
  fN     (source.fN),
  fChi2  (source.fChi2)
{
//
//  Copy constructor
//
  fIndex = new Int_t[fN0];
  fParams = new TClonesArray("TVectorD",fN0);
  fCovars = new TClonesArray("TMatrixD",fN0);
  fParams = (TClonesArray*)source.fParams->Clone();
  fCovars = (TClonesArray*)source.fCovars->Clone();
  for (Int_t i=0; i<fN0; i++) fIndex[i] = source.fIndex[i];
  
  fX     = new Double_t[fN];
  fY0    = new Double_t[fN];
  fY1    = new Double_t[fN];
  fChi2I = new Double_t[fN];
  for (Int_t i=0; i<fN; i++){
    fX[i]  = source.fX[i];
    fY0[i] = source.fY0[i];
    fY1[i] = source.fY1[i];
  }
}
AliSplineFit& AliSplineFit::operator=(const AliSplineFit& source){
//
//  assignment operator
//    
  if (&source == this) return *this;

//
// reassign memory as previous fit could have a different size
//

  if ( fN0 != source.fN0) { 

    delete fParams;
    delete fCovars;
    delete []fIndex;

    fN0 = source.fN0;
    fIndex = new Int_t[fN0];
    fParams = new TClonesArray("TVectorD",fN0);
    fCovars = new TClonesArray("TMatrixD",fN0);
  }
  if ( fN != source.fN) { 

    delete []fX;
    delete []fY0;
    delete []fY1;
    delete []fChi2I;
    fN = source.fN;  
    fX     = new Double_t[fN];
    fY0    = new Double_t[fN];
    fY1    = new Double_t[fN];
    fChi2I = new Double_t[fN];
  }

// use copy constructor (without reassigning memory) to copy values
  
  new (this) AliSplineFit(source);
  
  return *this;
}

  
AliSplineFit::~AliSplineFit(){
  //
  // destructor. Don't delete fGraph, as this normally comes as input parameter
  //
  delete []fX;
  delete []fY0;
  delete []fY1;
  delete []fChi2I;
  delete fParams;
  delete fCovars;
  delete []fIndex;
}

Double_t   AliSplineFit::Eval(Double_t x, Int_t deriv) const{
  //
  // evaluate value at x
  //   deriv = 0: function value
  //         = 1: first derivative
  //         = 2: 2nd derivative
  //         = 3: 3rd derivative
  //
  //  a2 = -(3*a0 -3*b0 + 2*a1*dx +b1*dx)/(dx*dx)
  //  a3 = -(-2*a0+2*b0 - a1*dx - b1*dx)/(dx*dx*dx) 

  Int_t index = TMath::BinarySearch(fN,fX,x);
  if (index<0) index =0;
  if (index>fN-2) index =fN-2;
  //
  Double_t dx   = x-fX[index];
  Double_t dxc  = fX[index+1]-fX[index];
  Double_t y0   = fY0[index];
  Double_t y1   = fY1[index];
  Double_t y01  = fY0[index+1];
  Double_t y11  = fY1[index+1];
  Double_t y2   = -(3.*y0-3.*y01+2*y1*dxc+y11*dxc)/(dxc*dxc);
  Double_t y3   = -(-2.* y0 + 2*y01 -  y1*dxc - y11*dxc) /(dxc*dxc*dxc);
  Double_t val = y0+y1*dx+y2*dx*dx+y3*dx*dx*dx;
  if (deriv==1) val = y1+2.*y2*dx+3.*y3*dx*dx;
  if (deriv==2) val = 2.*y2+6.*y3*dx;
  if (deriv==3) val = 6*y3;
  return val;
}


TGraph * AliSplineFit::GenerGraph(Int_t npoints, Double_t fraction, Double_t s1, Double_t s2, Double_t s3, Int_t der){
  //
  // generate random graph 
  // xrange 0,1
  // yrange 0,1
  // s1, s2, s3 -  sigma of derivative
  // fraction   -  

  Double_t *value = new Double_t[npoints];
  Double_t *time  = new Double_t[npoints];
  Double_t d0=0, d1=0,d2=0,d3=0;
  value[0] = d0;
  time[0]  = 0;
  for(Int_t i=1; i<npoints; i++){
    Double_t dtime = 1./npoints;
    Double_t dd1 = dtime;
    Double_t dd2 = dd1*dd1;
    Double_t dd3 = dd2*dd1;
    d0 += d1*dd1 + d2*dd2/2. + d3*dd3/6.;
    d1 += d2*dd1 +d3*dd2/2;
    d2 += d3*dd1;
    value[i] = d0;
    time[i]  = time[i-1]+dtime;
    d1 =(1.-fraction)*d1+fraction*(gRandom->Exp(s1))*(gRandom->Rndm()-0.5);
    d2 =(1.-fraction)*d2+fraction*(gRandom->Exp(s2))*(gRandom->Rndm()-0.5);
    d3 =(1.-fraction)*d3+fraction*(gRandom->Exp(s3))*(gRandom->Rndm()-0.5);
    if (gRandom->Rndm()<fraction) d3 =(1.-fraction)*d3+fraction*(gRandom->BreitWigner(0,s3));
  }
  Double_t dmean = (value[npoints-1]-value[0])/(time[npoints-1]-time[0]);
  Double_t min = value[0];
  Double_t max = value[0];
  for (Int_t i=0; i<npoints; i++){
    value[i]  = value[i]-dmean*(time[i]-time[0]); 
    if (value[i]<min) min=value[i];
    if (value[i]>max) max=value[i];
  }

  for (Int_t i=0; i<npoints; i++){
    value[i]  = (value[i]-min)/(max-min); 
  }
  if (der==1) for (Int_t i=1; i<npoints; i++){
    value[i-1]  =  (value[i]-value[i-1])/(time[i]-time[i-1]);
  }

  TGraph * graph = new TGraph(npoints,time,value);
 
  delete [] value;
  delete [] time; 
  return graph;  
}


TGraph * AliSplineFit::GenerNoise(TGraph * graph0, Double_t sigma0){
  //
  // add noise to graph
  //

  Int_t npoints=graph0->GetN();
  Double_t *value = new Double_t[npoints];
  Double_t *time  = new Double_t[npoints];
  for(Int_t i=0; i<npoints; i++){
    time[i]  = graph0->GetX()[i];
    value[i] = graph0->GetY()[i]+gRandom->Gaus(0,sigma0);
  }
  TGraph * graph = new TGraph(npoints,time,value);

  delete [] value;
  delete [] time;
  return graph;  
}
 

TGraph * AliSplineFit::MakeGraph(Double_t xmin, Double_t xmax, Int_t npoints, Int_t deriv) const {
  //
  // if npoints<=0 draw derivative
  //

  TGraph *graph =0;
  if (npoints<=0) {
    if (deriv<=0) return new TGraph(fN,fX,fY0);
    if (deriv==1) return new TGraph(fN,fX,fY1);
    if (deriv>2) return new TGraph(fN-1,fX,fChi2I);
  }
  Double_t * x = new Double_t[npoints+1];
  Double_t * y = new Double_t[npoints+1];
  for (Int_t ip=0; ip<=npoints; ip++){
    x[ip] = xmin+ (xmax-xmin)*(Double_t(ip)/Double_t(npoints));
    y[ip] = Eval(x[ip],deriv);
  }

  graph = new TGraph(npoints,x,y);
  delete [] x;
  delete [] y;
  return graph;
}

TGraph * AliSplineFit::MakeDiff(TGraph * graph0) const {
  //
  // Make graph of difference to reference graph 
  //
  
  Int_t npoints=graph0->GetN();
  TGraph *graph =0;
  Double_t * x = new Double_t[npoints];
  Double_t * y = new Double_t[npoints];
  for (Int_t ip=0; ip<npoints; ip++){
    x[ip] = graph0->GetX()[ip];
    y[ip] = Eval(x[ip],0)-graph0->GetY()[ip];
  }
  graph = new TGraph(npoints,x,y);
  delete [] x;
  delete [] y;
  return graph;
}


TH1F * AliSplineFit::MakeDiffHisto(TGraph * graph0) const {
  //
  // Make histogram of difference to reference graph 
  //

  Int_t npoints=graph0->GetN();
  Float_t min=1e+39,max=-1e+39;
  for (Int_t ip=0; ip<npoints; ip++){
    Double_t x = graph0->GetX()[ip];
    Double_t y = Eval(x,0)-graph0->GetY()[ip];
    if (ip==0) {
      min = y;
      max = y;
    }else{
      if (y<min) min=y;
      if (y>max) max=y;
    }    
  }

  TH1F *his = new TH1F("hdiff","hdiff", 100, min, max);
  for (Int_t ip=0; ip<npoints; ip++){
    Double_t x = graph0->GetX()[ip];
    Double_t y = Eval(x,0)-graph0->GetY()[ip];
    his->Fill(y);
  }

  return his;
}



void AliSplineFit::InitKnots(TGraph * graph, Int_t min, Int_t iter, Double_t maxDelta){
  //
  // initialize knots + estimate sigma of noise + make initial parameters 
  //
  //

  const Double_t kEpsilon = 1.e-7;
  fGraph  = graph;
  fNmin   = min;
  fMaxDelta = maxDelta;
  Int_t npoints = fGraph->GetN();
  fN0           = (npoints/fNmin)+1;
  Float_t delta = Double_t(npoints)/Double_t(fN0-1);

  fParams = new TClonesArray("TVectorD",fN0);
  fCovars = new TClonesArray("TMatrixD",fN0);
  fIndex  = new Int_t[fN0];
  TLinearFitter fitterLocal(4,"pol3");  // local fitter
  Double_t sigma2 =0;


  Double_t yMin=graph->GetY()[0];
  Double_t yMax=graph->GetY()[0];

  for (Int_t iKnot=0; iKnot<fN0; iKnot++){
    Int_t index0 = TMath::Nint(Double_t(iKnot)*Double_t(delta));
    Int_t index1 = TMath::Min(TMath::Nint(Double_t(iKnot+1)*Double_t(delta)),npoints-1);
    Int_t indexM = (iKnot>0) ? fIndex[iKnot-1]:index0;
    fIndex[iKnot]=TMath::Min(index0, npoints-1);
    Float_t startX =graph->GetX()[fIndex[iKnot]];

    for (Int_t ipoint=indexM; ipoint<index1; ipoint++){
      Double_t dxl   =graph->GetX()[ipoint]-startX;
      Double_t  y    = graph->GetY()[ipoint];
      if (y<yMin) yMin=y;
      if (y>yMax) yMax=y;
      fitterLocal.AddPoint(&dxl,y,1);
    }

    fitterLocal.Eval();
    sigma2 += fitterLocal.GetChisquare()/Double_t((index1-indexM)-4.);
    TMatrixD   * covar = new ((*fCovars)[iKnot]) TMatrixD(4,4);
    TVectorD   * param = new ((*fParams)[iKnot]) TVectorD(4);
    fitterLocal.GetParameters(*param);
    fitterLocal.GetCovarianceMatrix(*covar);
    fitterLocal.ClearPoints();
  }
  fSigma  =TMath::Sqrt(sigma2/Double_t(fN0));   // mean sigma
  Double_t tDiff = ((yMax-yMin)+TMath::Abs(yMax)+TMath::Abs(yMin))*kEpsilon;
  fSigma += tDiff+fMaxDelta/TMath::Sqrt(npoints);
  fMaxDelta +=tDiff;
  for (Int_t iKnot=0; iKnot<fN0; iKnot++){
    TMatrixD & cov = *((TMatrixD*)fCovars->At(iKnot));
    cov*=fSigma*fSigma;
  }  
  OptimizeKnots(iter);

  fN = 0;
  for (Int_t iKnot=0; iKnot<fN0; iKnot++) if (fIndex[iKnot]>=0) fN++;
  fX  = new Double_t[fN];
  fY0 = new Double_t[fN];
  fY1 = new Double_t[fN];
  fChi2I = new Double_t[fN];
  Int_t iKnot=0;
  for (Int_t i=0; i<fN0; i++){
    if (fIndex[i]<0) continue; 
    if (iKnot>=fN) {
      printf("AliSplineFit::InitKnots: Knot number > Max knot number\n");
      break;
    }
    TVectorD   * param = (TVectorD*) fParams->At(i);
    fX[iKnot]  = fGraph->GetX()[fIndex[i]];
    fY0[iKnot] = (*param)(0);
    fY1[iKnot] = (*param)(1);    
    fChi2I[iKnot] = 0;
    iKnot++;
  }
}


Int_t AliSplineFit::OptimizeKnots(Int_t nIter){
  //
  //
  //
  const Double_t kMaxChi2= 5;
  Int_t nKnots=0;
  TTreeSRedirector cstream("SplineIter.root");
  for (Int_t iIter=0; iIter<nIter; iIter++){
    if (fBDump) cstream<<"Fit"<<
      "iIter="<<iIter<<
      "fit.="<<this<<
      "\n";
    nKnots=2;
    for (Int_t iKnot=1; iKnot<fN0-1; iKnot++){
      if (fIndex[iKnot]<0) continue;   //disabled knot
      Double_t chi2 = CheckKnot(iKnot); 
      Double_t startX = fGraph->GetX()[fIndex[iKnot]]; 
      if (fBDump) {
        TMatrixD   * covar = (TMatrixD*)fCovars->At(iKnot);
        TVectorD   * param = (TVectorD*)fParams->At(iKnot);
        cstream<<"Chi2"<<
	 "iIter="<<iIter<<
	 "iKnot="<<iKnot<<
	 "chi2="<<chi2<<
	 "x="<<startX<<
	 "param="<<param<<
	 "covar="<<covar<<
	 "\n";
      }
      if (chi2>kMaxChi2) { nKnots++;continue;}
      fIndex[iKnot]*=-1;
      Int_t iPrevious=iKnot-1;
      Int_t iNext    =iKnot+1;
      while (fIndex[iPrevious]<0) iPrevious--;
      while (fIndex[iNext]<0) iNext++;
      RefitKnot(iPrevious);
      RefitKnot(iNext);
      iKnot++;
      while (iKnot<fN0-1&& fIndex[iKnot]<0) iKnot++;
    }
  }
  return nKnots;
}


Bool_t   AliSplineFit::RefitKnot(Int_t iKnot){
  //
  //
  //

  Int_t iPrevious=(iKnot>0)  ?iKnot-1: 0;
  Int_t iNext    =(iKnot<fN0)?iKnot+1: fN0-1;
  while (iPrevious>0&&fIndex[iPrevious]<0) iPrevious--;
  while (iNext<fN0&&fIndex[iNext]<0) iNext++;
  if (iPrevious<0) iPrevious=0;
  if (iNext>=fN0) iNext=fN0-1;
  
  Double_t startX = fGraph->GetX()[fIndex[iKnot]]; 
  AliSplineFit::fitterStatic()->ClearPoints();
  Int_t indPrev = fIndex[iPrevious];
  Int_t indNext = fIndex[iNext];
  Double_t *graphX = fGraph->GetX();
  Double_t *graphY = fGraph->GetY();

  // make arrays for points to fit (to save time)

  Int_t nPoints = indNext-indPrev;
  Double_t *xPoint = new Double_t[3*nPoints];
  Double_t *yPoint = &xPoint[nPoints];
  Double_t *ePoint = &xPoint[2*nPoints];
  Int_t indVec=0;
  for (Int_t iPoint=indPrev; iPoint<indNext; iPoint++, indVec++){
    Double_t dxl   = graphX[iPoint]-startX;
    Double_t  y    = graphY[iPoint];
    xPoint[indVec] = dxl;
    yPoint[indVec] = y;
    ePoint[indVec] =  fSigma;
//    ePoint[indVec] =  fSigma+TMath::Abs(y)*kEpsilon;
//    AliSplineFit::fitterStatic.AddPoint(&dxl,y,fSigma+TMath::Abs(y)*kEpsilon);
  }
  AliSplineFit::fitterStatic()->AssignData(nPoints,1,xPoint,yPoint,ePoint);
  AliSplineFit::fitterStatic()->Eval();

//  delete temporary arrays

  delete [] xPoint; 
  
  TMatrixD   * covar = (TMatrixD*)fCovars->At(iKnot);
  TVectorD   * param = (TVectorD*)fParams->At(iKnot);
  AliSplineFit::fitterStatic()->GetParameters(*param);
  AliSplineFit::fitterStatic()->GetCovarianceMatrix(*covar);
  return 0;
}


Float_t AliSplineFit::CheckKnot(Int_t iKnot){
  //
  //
  //

  Int_t iPrevious=iKnot-1;
  Int_t iNext    =iKnot+1;
  while (fIndex[iPrevious]<0) iPrevious--;
  while (fIndex[iNext]<0) iNext++;
  TVectorD &pPrevious = *((TVectorD*)fParams->At(iPrevious));
  TVectorD &pNext     = *((TVectorD*)fParams->At(iNext));
  TVectorD &pKnot     = *((TVectorD*)fParams->At(iKnot));
  TMatrixD &cPrevious = *((TMatrixD*)fCovars->At(iPrevious));
  TMatrixD &cNext     = *((TMatrixD*)fCovars->At(iNext));
  TMatrixD &cKnot     = *((TMatrixD*)fCovars->At(iKnot));
  Double_t xPrevious = fGraph->GetX()[fIndex[iPrevious]];
  Double_t xNext     = fGraph->GetX()[fIndex[iNext]];
  Double_t xKnot     = fGraph->GetX()[fIndex[iKnot]];

  // extra variables introduced to save processing time

  Double_t dxc  = xNext-xPrevious;
  Double_t invDxc = 1./dxc;
  Double_t invDxc2 = invDxc*invDxc;
  TMatrixD  tPrevious(4,4);
  TMatrixD  tNext(4,4);

  tPrevious(0,0) = 1;    tPrevious(1,1) = 1;
  tPrevious(2,0) = -3.*invDxc2;
  tPrevious(2,1) = -2.*invDxc;
  tPrevious(3,0) =  2.*invDxc2*invDxc;
  tPrevious(3,1) =  1.*invDxc2;
  tNext(2,0)     =  3.*invDxc2;      tNext(2,1)     = -1*invDxc;
  tNext(3,0)     = -2.*invDxc2*invDxc;  tNext(3,1)     =  1.*invDxc2;
  TMatrixD  tpKnot(4,4);
  TMatrixD  tpNext(4,4);
  Double_t dx = xKnot-xPrevious;
  tpKnot(0,0) = 1;      tpKnot(1,1) = 1;        tpKnot(2,2) = 1;        tpKnot(3,3) = 1;
  tpKnot(0,1) = dx;     tpKnot(0,2) = dx*dx;    tpKnot(0,3) = dx*dx*dx;
  tpKnot(1,2) = 2.*dx;  tpKnot(1,3) = 3.*dx*dx;
  tpKnot(2,3) = 3.*dx;
  Double_t dxn = xNext-xPrevious;
  tpNext(0,0) = 1;       tpNext(1,1) = 1;        tpNext(2,2) = 1;        tpNext(3,3) = 1;
  tpNext(0,1) = dxn;     tpNext(0,2) = dxn*dxn;    tpNext(0,3) = dxn*dxn*dxn;
  tpNext(1,2) = 2.*dxn;  tpNext(1,3) = 3.*dxn*dxn;
  tpNext(2,3) = 3.*dxn;

  //
  // matrix and vector at previous
  //

  TVectorD  sPrevious = tPrevious*pPrevious+tNext*pNext;
  TVectorD  sKnot     = tpKnot*sPrevious;
  TVectorD  sNext     = tpNext*sPrevious;  
  
  TMatrixD csPrevious00(tPrevious, TMatrixD::kMult,cPrevious);
  csPrevious00 *= tPrevious.T();
  TMatrixD csPrevious01(tNext,TMatrixD::kMult,cNext);
  csPrevious01*=tNext.T();
  TMatrixD  csPrevious(csPrevious00,TMatrixD::kPlus,csPrevious01);
  TMatrixD  csKnot(tpKnot,TMatrixD::kMult,csPrevious);
  csKnot*=tpKnot.T();
  TMatrixD  csNext(tpNext,TMatrixD::kMult,csPrevious);
  csNext*=tpNext.T();

  TVectorD dPrevious = pPrevious-sPrevious;
  TVectorD dKnot     = pKnot-sKnot;
  TVectorD dNext     = pNext-sNext;
  //
  //
  TMatrixD prec(4,4);
  prec(0,0) = (fMaxDelta*fMaxDelta);
  prec(1,1) = prec(0,0)*invDxc2;
  prec(2,2) = prec(1,1)*invDxc2;
  prec(3,3) = prec(2,2)*invDxc2;

//   prec(0,0) = (fMaxDelta*fMaxDelta);
//   prec(1,1) = (fMaxDelta*fMaxDelta)/(dxc*dxc);
//   prec(2,2) = (fMaxDelta*fMaxDelta)/(dxc*dxc*dxc*dxc);
//   prec(3,3) = (fMaxDelta*fMaxDelta)/(dxc*dxc*dxc*dxc*dxc*dxc);

  csPrevious+=cPrevious;
  csPrevious+=prec;
  csPrevious.Invert(); 
  Double_t  chi2P     = dPrevious*(csPrevious*dPrevious);

  csKnot+=cKnot;
  csKnot+=prec;
  csKnot.Invert();
  Double_t  chi2K     = dKnot*(csKnot*dKnot);
  
  csNext+=cNext;
  csNext+=prec;
  csNext.Invert();
  Double_t  chi2N     = dNext*(csNext*dNext);    

  return (chi2P+chi2K+chi2N)/8.;


}
                  
void AliSplineFit::SplineFit(Int_t nder){
  //
  // Cubic spline fit of graph
  // 
  // nder
  // nder<0  - no continuity requirement
  //     =0  - continous  0 derivative
  //     =1  - continous  1 derivative
  //     >1  - continous  2 derivative  
  //
  if (!fGraph) return;
  TGraph * graph = fGraph;
  if (nder>1) nder=2;
  Int_t nknots  = fN;
  Int_t npoints = graph->GetN(); 
  //
  //
  // spline fit
  // each knot 4 parameters  
  // 
  TMatrixD       *pmatrix = 0;
  TVectorD       *pvalues = 0;  
  if (nder>1){
    pmatrix = new TMatrixD(4*(nknots-1)+3*(nknots-2), 4*(nknots-1)+3*(nknots-2));
    pvalues = new TVectorD(4*(nknots-1)+3*(nknots-2));
  }
  if (nder==1){
    pmatrix = new TMatrixD(4*(nknots-1)+2*(nknots-2), 4*(nknots-1)+2*(nknots-2));
    pvalues = new TVectorD(4*(nknots-1)+2*(nknots-2));
  }
  if (nder==0){
    pmatrix = new TMatrixD(4*(nknots-1)+1*(nknots-2), 4*(nknots-1)+1*(nknots-2));
    pvalues = new TVectorD(4*(nknots-1)+1*(nknots-2));
  }
  if (nder<0){
    pmatrix = new TMatrixD(4*(nknots-1)+0*(nknots-2), 4*(nknots-1)+0*(nknots-2));
    pvalues = new TVectorD(4*(nknots-1)+0*(nknots-2));
  }
  
  
  TMatrixD &matrix = *pmatrix;
  TVectorD &values = *pvalues;
  Int_t    current = 0;
//
//  defined extra variables (current4 etc.) to save processing time.
//  fill normal matrices, then copy to sparse matrix.
//  
  Double_t *graphX = graph->GetX();
  Double_t *graphY = graph->GetY();
  for (Int_t ip=0;ip<npoints;ip++){
    if (current<nknots-2&&graphX[ip]>fX[current+1]) current++;
    Double_t xmiddle = (fX[current+1]+fX[current])*0.5;
    Double_t x1 = graphX[ip]- xmiddle;
    Double_t x2 = x1*x1;
    Double_t x3 = x2*x1;
    Double_t x4 = x2*x2;
    Double_t x5 = x3*x2;
    Double_t x6 = x3*x3;
    Double_t y  = graphY[ip];
    Int_t current4 = 4*current;

    matrix(current4  , current4  )+=1;
    matrix(current4  , current4+1)+=x1;
    matrix(current4  , current4+2)+=x2;
    matrix(current4  , current4+3)+=x3;
    //
    matrix(current4+1, current4  )+=x1;
    matrix(current4+1, current4+1)+=x2;
    matrix(current4+1, current4+2)+=x3;
    matrix(current4+1, current4+3)+=x4;
    //
    matrix(current4+2, current4  )+=x2;
    matrix(current4+2, current4+1)+=x3;
    matrix(current4+2, current4+2)+=x4;
    matrix(current4+2, current4+3)+=x5;
    //
    matrix(current4+3, current4  )+=x3;
    matrix(current4+3, current4+1)+=x4;
    matrix(current4+3, current4+2)+=x5;
    matrix(current4+3, current4+3)+=x6;
    //
    values(current4  ) += y;
    values(current4+1) += y*x1;
    values(current4+2) += y*x2;
    values(current4+3) += y*x3;
  }
  //
  // constraint 0
  //
  Int_t offset =4*(nknots-1)-1;
  if (nder>=0) for (Int_t iknot = 1; iknot<nknots-1; iknot++){

    Double_t dxm  =  (fX[iknot]-fX[iknot-1])*0.5;
    Double_t dxp  = -(fX[iknot+1]-fX[iknot])*0.5;
    Double_t dxm2 = dxm*dxm;
    Double_t dxp2 = dxp*dxp;
    Double_t dxm3 = dxm2*dxm;
    Double_t dxp3 = dxp2*dxp;
    Int_t iknot4  = 4*iknot;
    Int_t iknot41 = 4*(iknot-1);
    Int_t offsKnot = offset+iknot;
    //
    // condition on knot
    //
    // a0[i] = a0m[i-1]  + a1m[i-1]*dxm + a2m[i-1]*dxm^2 + a3m[i-1]*dxm^3
    // a0[i] = a0m[i-0]  + a1m[i-0]*dxp + a2m[i-0]*dxp^2 + a3m[i-0]*dxp^3
    // (a0m[i-1]  + a1m[i-1]*dxm + a2m[i-1]*dxm^2 + a3m[i-1]*dxm^3) -
    // (a0m[i-0]  + a1m[i-0]*dxp + a2m[i-0]*dxp^2 + a3m[i-0]*dxp^3)  = 0
    
    matrix(offsKnot, iknot41  )=1;
    matrix(offsKnot, iknot4   )=-1;

    matrix(offsKnot, iknot41+1)=dxm;
    matrix(offsKnot, iknot4 +1)=-dxp;

    matrix(offsKnot, iknot41+2)=dxm2;
    matrix(offsKnot, iknot4 +2)=-dxp2;

    matrix(offsKnot, iknot41+3)=dxm3;
    matrix(offsKnot, iknot4 +3)=-dxp3;

    matrix(iknot41  , offsKnot)=1;
    matrix(iknot41+1, offsKnot)=dxm;
    matrix(iknot41+2, offsKnot)=dxm2;
    matrix(iknot41+3, offsKnot)=dxm3;
    matrix(iknot4  , offsKnot)=-1;
    matrix(iknot4+1, offsKnot)=-dxp;
    matrix(iknot4+2, offsKnot)=-dxp2;
    matrix(iknot4+3, offsKnot)=-dxp3;
  }
  //
  // constraint 1
  //
  offset =4*(nknots-1)-1+(nknots-2);
  if (nder>=1)for (Int_t iknot = 1; iknot<nknots-1; iknot++){

    Double_t dxm  =  (fX[iknot]-fX[iknot-1])*0.5;
    Double_t dxp  = -(fX[iknot+1]-fX[iknot])*0.5;
    Double_t dxm2 = dxm*dxm;
    Double_t dxp2 = dxp*dxp;
    Int_t iknot4  = 4*iknot;
    Int_t iknot41 = 4*(iknot-1);
    Int_t offsKnot = offset+iknot;
    //
    // condition on knot derivation
    //
    // a0d[i] =  a1m[i-1] + 2*a2m[i-1]*dxm + 3*a3m[i-1]*dxm^2
    // a0d[i] =  a1m[i-0] + 2*a2m[i-0]*dxp + 3*a3m[i-0]*dxp^2
    
    //
    matrix(offsKnot, iknot41+1)= 1;
    matrix(offsKnot, iknot4 +1)=-1;

    matrix(offsKnot, iknot41+2)= 2.*dxm;
    matrix(offsKnot, iknot4 +2)=-2.*dxp;

    matrix(offsKnot, iknot41+3)= 3.*dxm2;
    matrix(offsKnot, iknot4 +3)=-3.*dxp2;

    matrix(iknot41+1, offsKnot)=1;
    matrix(iknot41+2, offsKnot)=2.*dxm;
    matrix(iknot41+3, offsKnot)=3.*dxm2;

    matrix(iknot4+1, offsKnot)=-1.;
    matrix(iknot4+2, offsKnot)=-2.*dxp;
    matrix(iknot4+3, offsKnot)=-3.*dxp2;
  }
  //
  // constraint 2
  //
  offset =4*(nknots-1)-1+2*(nknots-2);
  if (nder>=2) for (Int_t iknot = 1; iknot<nknots-1; iknot++){

    Double_t dxm  =  (fX[iknot]-fX[iknot-1])*0.5;
    Double_t dxp  = -(fX[iknot+1]-fX[iknot])*0.5;
    Int_t iknot4  = 4*iknot;
    Int_t iknot41 = 4*(iknot-1);
    Int_t offsKnot = offset+iknot;
    //
    // condition on knot second derivative
    //
    // a0dd[i] =  2*a2m[i-1] + 6*a3m[i-1]*dxm
    // a0dd[i] =  2*a2m[i-0] + 6*a3m[i-0]*dxp    
    //
    //
    matrix(offsKnot, iknot41+2)= 2.;
    matrix(offsKnot, iknot4 +2)=-2.;

    matrix(offsKnot, iknot41+3)= 6.*dxm;
    matrix(offsKnot, iknot4 +3)=-6.*dxp;

    matrix(iknot41+2, offsKnot)=2.;
    matrix(iknot41+3, offsKnot)=6.*dxm;

    matrix(iknot4+2, offsKnot)=-2.;
    matrix(iknot4+3, offsKnot)=-6.*dxp;
  }
 
// sparse matrix to do fit
  
  TMatrixDSparse smatrix(matrix);
  TDecompSparse svd(smatrix,0);
  Bool_t ok;
  const TVectorD results = svd.Solve(values,ok);

  for (Int_t iknot = 0; iknot<nknots-1; iknot++){

    Double_t dxm  =  -(fX[iknot+1]-fX[iknot])*0.5;

    fY0[iknot] = results(4*iknot)+ results(4*iknot+1)*dxm+results(4*iknot+2)*dxm*dxm+
      results(4*iknot+3)*dxm*dxm*dxm;

    fY1[iknot] = results(4*iknot+1)+2.*results(4*iknot+2)*dxm+
      3*results(4*iknot+3)*dxm*dxm;
  }
  Int_t   iknot2= nknots-1;
  Int_t   iknot = nknots-2;
  Double_t dxm   =  (fX[iknot2]-fX[iknot2-1])*0.5;

  fY0[iknot2] = results(4*iknot)+ results(4*iknot+1)*dxm+results(4*iknot+2)*dxm*dxm+
    results(4*iknot+3)*dxm*dxm*dxm;

  fY1[iknot2] = results(4*iknot+1)+2.*results(4*iknot+2)*dxm+
      3*results(4*iknot+3)*dxm*dxm;

  delete  pmatrix;
  delete  pvalues;

}





void AliSplineFit::MakeKnots0(TGraph * graph, Double_t maxdelta, Int_t minpoints){
  //
  // make knots  - restriction max distance and minimum points
  //

  Int_t npoints  = graph->GetN();
  Double_t *xknots = new Double_t[npoints];
  Int_t nknots =0;
  Int_t ipoints =0;
  //
  // generate knots
  //
  for (Int_t ip=0;ip<npoints;ip++){
    if (graph->GetX()[ip]-xknots[nknots-1]>maxdelta && ipoints>minpoints){
      xknots[nknots] = graph->GetX()[ip];
      ipoints=1;
      nknots++;
    }
    ipoints++;
  }
  if (npoints-ipoints>minpoints){
    xknots[nknots] = graph->GetX()[npoints-1];
    nknots++;
  }else{
    xknots[nknots-1] = graph->GetX()[npoints-1];
  }

  fN = nknots;
  fX = new Double_t[nknots];
  fY0 = new Double_t[nknots];
  fY1 = new Double_t[nknots];
  fChi2I= new Double_t[nknots];
  for (Int_t i=0; i<nknots; i++) fX[i]= xknots[i];  
  delete [] xknots;
}




void AliSplineFit::MakeSmooth(TGraph * graph, Float_t ratio, char * type){
  //
  // Interface to GraphSmooth
  //

  TGraphSmooth smooth;
  Int_t    npoints2 = TMath::Nint(graph->GetN()*ratio);  
  TGraph * graphT0 = smooth.SmoothKern(graph,type,ratio);
  if (!graphT0) return;
  TGraph  graphT1(npoints2);
  for (Int_t ipoint=0; ipoint<npoints2; ipoint++){
    Int_t pointS = TMath::Nint(ipoint/ratio);
    if (ipoint==npoints2-1) pointS=graph->GetN()-1;
    graphT1.SetPoint(ipoint, graphT0->GetX()[pointS] , graphT0->GetY()[pointS]);
  }  
  TSpline3 spline2("spline", &graphT1);
  Update(&spline2, npoints2);
}


void AliSplineFit::Update(TSpline3 *spline, Int_t nknots){
  //
  //
  //

  fN = nknots;
  fX = new Double_t[nknots];
  fY0 = new Double_t[nknots];
  fY1 = new Double_t[nknots];
  Double_t d0, d1;
  fChi2I= 0;
  for (Int_t i=0; i<nknots; i++) {
    spline->GetCoeff(i,fX[i],fY0[i], fY1[i],d0,d1);
  }
}




void AliSplineFit::Test(Int_t npoints, Int_t ntracks, Float_t snoise){  
  //
  // test function
  //

  AliSplineFit fit;
  AliSplineFit fitS;
  TGraph * graph0=0;
  TGraph * graph1=0;
  
  TTreeSRedirector *pcstream = new TTreeSRedirector("TestSmooth.root");
  for (Int_t i=0; i<ntracks; i++){
    graph0 = AliSplineFit::GenerGraph(npoints,0.05,0,0,1,0);  
    graph1 = AliSplineFit::GenerNoise(graph0,snoise);  
    fit.InitKnots(graph1, 10,10, 0.00);
    TGraph *d0 = fit.MakeDiff(graph0);
    TGraph *g0 = fit.MakeGraph(0,1,1000,0);
    fit.SplineFit(2);
    TH1F * h2 = fit.MakeDiffHisto(graph0);
    TGraph *d2 = fit.MakeDiff(graph0);
    TGraph *g2 = fit.MakeGraph(0,1,1000,0);
    fit.SplineFit(1);
    TH1F * h1 = fit.MakeDiffHisto(graph0);
    TGraph *d1 = fit.MakeDiff(graph0);
    TGraph *g1 = fit.MakeGraph(0,1,1000,0);

    Float_t ratio = Float_t(fit.fN)/Float_t(npoints);
    fitS.MakeSmooth(graph1,ratio,"box");
    TGraph *dS = fitS.MakeDiff(graph0);
    TGraph *gS = fit.MakeGraph(0,1,1000,0);

    TH1F * hS = fitS.MakeDiffHisto(graph0);
    Double_t mean2  = h2->GetMean();
    Double_t sigma2 = h2->GetRMS();
    Double_t mean1  = h1->GetMean();
    Double_t sigma1 = h1->GetRMS();
    Double_t meanS  = hS->GetMean();
    Double_t sigmaS = hS->GetRMS();
    char fname[100];
    if (fit.fN<20){
      sprintf(fname,"pol%d",fit.fN);
    }else{
      sprintf(fname,"pol%d",19);
    }
    TF1 fpol("fpol",fname);
    graph1->Fit(&fpol);
    TGraph dpol(*graph1);
    TGraph gpol(*graph1);
    for (Int_t ipoint=0; ipoint<graph1->GetN(); ipoint++){
      dpol.GetY()[ipoint]= graph0->GetY()[ipoint]-
	fpol.Eval(graph0->GetX()[ipoint]);
      gpol.GetY()[ipoint]= fpol.Eval(graph0->GetX()[ipoint]);
    }
    (*pcstream)<<"Test"<<
      "Event="<<i<<
      "Graph0.="<<graph0<<
      "Graph1.="<<graph1<<
      "G0.="<<g0<<
      "G1.="<<g1<<
      "G2.="<<g2<<
      "GS.="<<gS<<
      "GP.="<<&gpol<<
      "D0.="<<d0<<
      "D1.="<<d1<<
      "D2.="<<d2<<
      "DS.="<<dS<<
      "DP.="<<&dpol<<
      "Npoints="<<fit.fN<<
      "Mean1="<<mean1<<
      "Mean2="<<mean2<<
      "MeanS="<<meanS<<
      "Sigma1="<<sigma1<<
      "Sigma2="<<sigma2<<
      "SigmaS="<<sigmaS<<
      "\n";
    
    delete graph0;
    delete graph1;
    delete g1;
    delete g2;
    delete gS;
    delete h1;
    delete h2;
    delete hS;
  }
  delete pcstream;    
}
