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

//-----------------------------------------------------------------
//   Implementation of the derived class for track residuals
//   based on linear chi2 minimization 
//  The minimization relies on the fact that the alignment parameters     
//   (angles and translations) are small.                                  
//   TLinearFitter used for minimization
//   Possibility to fix Paramaters
//   FixParameter()ReleaseParameter();
//   Possibility to define fraction of outliers to be skipped
//
//   marian.ivanov@cern.ch
//
//-----------------------------------------------------------------

#include <TMath.h>
#include <TMinuit.h>
#include <TGeoMatrix.h>

#include "AliLog.h"
#include "AliAlignObj.h"
#include "AliTrackPointArray.h"
#include "AliTrackResidualsLinear.h"
#include "AliAlignObj.h"
#include "TLinearFitter.h"
#include  "TDecompSVD.h"

ClassImp(AliTrackResidualsLinear)

//______________________________________________________________________________
AliTrackResidualsLinear::AliTrackResidualsLinear():
  AliTrackResiduals(),
  fFitter(0),
  fFraction(-1),
  fChi2Orig(0)
{
  // Default constructor
  for (Int_t ipar=0; ipar<6; ipar++){
    fBFixed[ipar] = kFALSE;
    fFixed[ipar]  = 0;;
    fParams[ipar]  = 0;
  }  
  for (Int_t icov=0; icov<36; icov++){ fCovar[icov]=0;}
}

//______________________________________________________________________________
AliTrackResidualsLinear::AliTrackResidualsLinear(Int_t ntracks):
  AliTrackResiduals(ntracks),
  fFitter(new TLinearFitter(6,"hyp6")),
  fFraction(-1),
  fChi2Orig(0)
{
  // Constructor
  for (Int_t ipar=0; ipar<6; ipar++){
    fBFixed[ipar] = kFALSE;
    fFixed[ipar]  = 0;
    fParams[ipar]  = 0;
  }
  for (Int_t icov=0; icov<36; icov++){ fCovar[icov]=0;}
}
 
//______________________________________________________________________________
AliTrackResidualsLinear::AliTrackResidualsLinear(const AliTrackResidualsLinear &res):
  AliTrackResiduals(res),
  fFitter(new TLinearFitter(*(res.fFitter))),
  fFraction(res.fFraction),
  fChi2Orig(res.fChi2Orig)
{
  // Copy constructor
  //..
  for (Int_t ipar=0; ipar<6; ipar++){
    fBFixed[ipar]  = res.fBFixed[ipar];
    fFixed[ipar]   = res.fFixed[ipar];
    fParams[ipar]  = res.fParams[ipar];
  }
  for (Int_t icov=0; icov<36; icov++){ fCovar[icov]= res.fCovar[icov];}
  fChi2Orig = res.fChi2Orig;
}

//______________________________________________________________________________
AliTrackResidualsLinear &AliTrackResidualsLinear::operator= (const AliTrackResidualsLinear& res)
{
  // Assignment operator
 ((AliTrackResiduals *)this)->operator=(res);
 return *this;
}
//______________________________________________________________________________
AliTrackResidualsLinear::~AliTrackResidualsLinear()
{
  //
  //
  //
  delete fFitter;
}


//______________________________________________________________________________
Bool_t AliTrackResidualsLinear::Minimize()
{
  // Implementation of fast linear Chi2 minimizer
  // based on TLinear fitter
  //
  if (!fFitter) fFitter = new TLinearFitter(6,"hyp6");
  fFitter->StoreData(kTRUE);
  fFitter->ClearPoints();
  fChi2Orig = 0;
  AliTrackPoint p1,p2;
  for (Int_t itrack = 0; itrack < fLast; itrack++) {
    if (!fVolArray[itrack] || !fTrackArray[itrack]) continue;
    for (Int_t ipoint = 0; ipoint < fVolArray[itrack]->GetNPoints(); ipoint++) {
      fVolArray[itrack]->GetPoint(p1,ipoint);
      fTrackArray[itrack]->GetPoint(p2,ipoint);
      AddPoints(p1,p2);
    }
  }
  Bool_t isOK = Update();
  if (!isOK) return isOK;
  //
  TGeoHMatrix  matrix;
  fAlignObj->GetMatrix(matrix);
  return isOK;
}



//______________________________________________________________________________
void AliTrackResidualsLinear::AddPoints(AliTrackPoint &p, AliTrackPoint &pprime)
{
  //
  // add points to linear fitter - option with correlation betwee measurement in different dimensions
  // p1      - local point
  // pprime  - track extrapolation point
  //
  Float_t xyz[3],xyzp[3];
  Float_t cov[6],covp[6];
  p.GetXYZ(xyz,cov); pprime.GetXYZ(xyzp,covp);
  //
  //
  TMatrixD mcov(3,3);   // local point covariance
  mcov(0,0) = cov[0]; mcov(0,1) = cov[1]; mcov(0,2) = cov[2];
  mcov(1,0) = cov[1]; mcov(1,1) = cov[3]; mcov(1,2) = cov[4];
  mcov(2,0) = cov[2]; mcov(2,1) = cov[4]; mcov(2,2) = cov[5];
  TMatrixD mcovp(3,3); //  extrapolation point covariance  
  mcovp(0,0) = covp[0]; mcovp(0,1) = covp[1]; mcovp(0,2) = covp[2];
  mcovp(1,0) = covp[1]; mcovp(1,1) = covp[3]; mcovp(1,2) = covp[4];
  mcovp(2,0) = covp[2]; mcovp(2,1) = covp[4]; mcovp(2,2) = covp[5];
  mcov+=mcovp;
  //mcov.Invert();
  if (!mcov.IsValid()) return; 
  TMatrixD mcovBack = mcov;  // for debug purposes
  //
  // decompose matrix
  //
  TDecompSVD svd(mcov);              // mcov  = svd.fU * covDiagonal * svd.fV.Invert   
  if (!svd.Decompose()) return;      // decomposition failed
  TMatrixD   matrixV = svd.GetV();   // transformation matrix to diagonalize covariance matrix
  Double_t   covDiagonal[3] = {svd.GetSig()[0],svd.GetSig()[1],svd.GetSig()[2]};    // diagonalized covariance matrix
  //
  // residual vector 
  TMatrixD  deltaR(3,1);
  deltaR(0,0) = (xyzp[0]-xyz[0]); 
  deltaR(1,0) = (xyzp[1]-xyz[1]);
  deltaR(2,0) = (xyzp[2]-xyz[2]);   
  //
  // parametrization matrix
  //
  TMatrixD        mparam(3,6);
  mparam(0,0) = 1;      mparam(1,0) = 0;       mparam(2,0) = 0;            // xshift
  mparam(0,1) = 0;      mparam(1,1) = 1;       mparam(2,1) = 0;            // yshift
  mparam(0,2) = 0;      mparam(1,2) = 0;       mparam(2,2) = 1;            // zshift
  mparam(0,3) = 0;      mparam(1,3) =-xyz[2];  mparam(2,3) = xyz[1];       // x rotation
  mparam(0,4) = xyz[2]; mparam(1,4) = 0;       mparam(2,4) =-xyz[0];       // y rotation
  mparam(0,5) =-xyz[1]; mparam(1,5) = xyz[0];  mparam(2,5) = 0;            // z rotation
  //
  
  TMatrixD  deltaT(matrixV, TMatrixD::kTransposeMult, deltaR);   // tranformed delta
  TMatrixD  mparamT(matrixV,TMatrixD::kTransposeMult, mparam);   // tranformed linear transformation
  if (AliLog::GetDebugLevel("","AliTrackResidualsLinear")>2){    
    //
    // debug part
    //
    //   covDiag = U^-1 * mcov * V      -- diagonalization of covariance matrix

    TMatrixD   matrixU = svd.GetU();                      // transformation matrix to diagonalize covariance matrix
    TMatrixD   matrixUI= svd.GetU(); 
    matrixUI.Invert();
    //
    TMatrixD   test0   = matrixUI*matrixV;                // test matrix - should be unit matrix
    TMatrixD   test1   = matrixUI*mcovBack*matrixV;       // test matrix - diagonal - should be diagonal with covDiagonal on diag
    TMatrixD   test2   = matrixU.T()*matrixV;             // test ortogonality - shoul be unit
    printf("Test matrix 2 - should be unit\n");
    test2.Print();
    printf("Test matrix 0 - should be unit\n"); 
    test0.Print();
    printf("Test matrix 1 - should be diagonal\n"); 
    test1.Print();
    printf("Diagonal matrix\n"); 
    svd.GetSig().Print();
    printf("Original param matrix\n"); 
    mparam.Print();
    printf("Rotated  param matrix\n"); 
    mparamT.Print();
    //
    //
    printf("Trans Matrix:\n");
    matrixV.Print();
    printf("Delta Orig\n");
    deltaR.Print();
    printf("Delta Rotated");
    deltaT.Print();
    //
    //    
  }
  //
  Double_t sumChi2=0;
  for (Int_t idim = 1; idim<3; idim++){
    Double_t yf;     // input points to fit in TLinear fitter
    Double_t xf[6];  // input points to fit
    yf = deltaT(idim,0);
    for (Int_t ipar =0; ipar<6; ipar++) xf[ipar] = mparamT(idim,ipar); 
    if (covDiagonal[idim]>0.){
      fFitter->AddPoint(xf,yf, TMath::Sqrt(covDiagonal[idim]));
      // accumulate chi2
      Double_t chi2       = (yf*yf)/covDiagonal[idim];
      fChi2Orig += (yf*yf)/covDiagonal[idim];  
      if (chi2>100 && AliLog::GetDebugLevel("","AliTrackResidualsLinear")>1){
	printf("Too big chi2- %f\n",chi2);
	printf("Delta Orig\n");
	deltaR.Print();
	printf("Delta Rotated");
	deltaT.Print();
	matrixV.Print();
	printf("Too big chi2 - End\n");	
      }
      sumChi2+=chi2;
    }
    else{
      printf("Bug\n");
    }
  }
  if (AliLog::GetDebugLevel("","AliTrackResidualsLinear")>1){
    TMatrixD matChi0=(mcov.Invert()*deltaR);
    TMatrixD matChi2=deltaR.T()*matChi0;
    printf("Chi2:\t%f\t%f", matChi2(0,0), sumChi2);
  }

  fNdf +=2;
}

//______________________________________________________________________________
Bool_t AliTrackResidualsLinear::Update()
{
  // Find the alignment parameters
  // using TLinear fitter  + fill data containers
  // 
  //
  fFitter->Eval();
  //
  // TLinear fitter put as first parameter offset - fixing parameter shifted by one
  //
  fFitter->FixParameter(0);
  for (Int_t ipar =0; ipar<6; ipar++){
    if (fBFixed[ipar])  fFitter->FixParameter(ipar+1,fFixed[ipar]);
  }
  if (fFraction>0.5) {
    fFitter->EvalRobust(fFraction);
  }else{
    fFitter->Eval();
  }
  //
  fFitter->ReleaseParameter(0);
  for (Int_t ipar=0; ipar<6; ipar++) {
    if (fBFixed[ipar])  fFitter->ReleaseParameter(ipar+1);
  }
    

  //
  fChi2 = fFitter->GetChisquare();
  fNdf -= 6;
  TVectorD vector(7);
  fFitter->GetParameters(vector);
  fParams[0] = vector[1];
  fParams[1] = vector[2];
  fParams[2] = vector[3];  
  fParams[3] = vector[4];
  fParams[4] = vector[5];
  fParams[5] = vector[6];
  TMatrixD covar(7,7);
  fFitter->GetCovarianceMatrix(covar);
  for (Int_t i0=0; i0 <6; i0++)
    for (Int_t j0=0; j0 <6; j0++){
      fCovar[i0*6+j0] = covar(i0+1,j0+1);
    }
  //
  fAlignObj->SetPars(fParams[0], fParams[1], fParams[2],
		     TMath::RadToDeg()*fParams[3],
		     TMath::RadToDeg()*fParams[4],
		     TMath::RadToDeg()*fParams[5]);
  return kTRUE;
}
