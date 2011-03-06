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

// Implemenatation of the K-Means Clustering Algorithm
// See http://en.wikipedia.org/wiki/K-means_clustering and references therein.
//
// This particular implementation is the so called Soft K-means algorithm.
// It has been modified to work on the cylindrical topology in eta-phi space.
//
// Author: Andreas Morsch (CERN)
// andreas.morsch@cern.ch
 
#include "AliKMeansClustering.h"
#include <TMath.h>
#include <TRandom.h>
#include <TH1F.h>

ClassImp(AliKMeansClustering)

Double_t AliKMeansClustering::fBeta = 10.;

 
Int_t AliKMeansClustering::SoftKMeans(Int_t k, Int_t n, const Double_t* x, const Double_t* y, Double_t* mx, Double_t* my , Double_t* rk )
{
    //
    // The soft K-means algorithm
    //
    Int_t i,j;
    //
    // (1) Initialisation of the k means

    for (i = 0; i < k; i++) {
	mx[i] = 2. * TMath::Pi() * gRandom->Rndm();
	my[i] = -1. + 2. * gRandom->Rndm();
    }

    //
    // (2a) The responsibilities
    Double_t** r = new Double_t*[n]; // responsibilities
    for (j = 0; j < n; j++) {r[j] = new Double_t[k];}
    //
    // (2b) Normalisation
    Double_t* nr   = new Double_t[n];
    // (3) Iterations
    Int_t nit = 0;
    
    while(1) {
	nit++;
      //
      // Assignment step
      //
      for (j = 0; j < n; j++) {
	nr[j] = 0.;
	for (i = 0; i < k; i++) {
	  r[j][i] = TMath::Exp(- fBeta * d(mx[i], my[i], x[j], y[j]));
	  nr[j] += r[j][i];
	} // mean i
      } // data point j
	
      for (j = 0; j < n; j++) {
	for (i = 0; i < k; i++) {
	  r[j][i] /=  nr[j];
	} // mean i
      } // data point j
      
	//
	// Update step
      Double_t di = 0;
      
      for (i = 0; i < k; i++) {
	  Double_t oldx = mx[i];
	  Double_t oldy = my[i];
	  
	  mx[i] = x[0];
	  my[i] = y[0];
	  rk[i] = r[0][i];
	
	for (j = 1; j < n; j++) {
	    Double_t xx =  x[j];
//
// Here we have to take into acount the cylinder topology where phi is defined mod 2xpi
// If two coordinates are separated by more than pi in phi one has to be shifted by +/- 2 pi

	    Double_t dx = mx[i] - x[j];
	    if (dx >  TMath::Pi()) xx += 2. * TMath::Pi();
	    if (dx < -TMath::Pi()) xx -= 2. * TMath::Pi();
	    mx[i] = mx[i] * rk[i] + r[j][i] * xx;
	    my[i] = my[i] * rk[i] + r[j][i] * y[j];
	    rk[i] += r[j][i];
	    mx[i] /= rk[i];
	    my[i] /= rk[i];	    
	    if (mx[i] > 2. * TMath::Pi()) mx[i] -= 2. * TMath::Pi();
	    if (mx[i] < 0.              ) mx[i] += 2. * TMath::Pi();
	} // Data
	di += d(mx[i], my[i], oldx, oldy);
      } // means 
	//
	// ending condition
      if (di < 1.e-8 || nit > 1000) break;
    } // while

// Clean-up    
    delete[] nr;
    for (j = 0; j < n; j++) delete[] r[j];
    delete[] r;
// 
    return (nit < 1000);
    
}

Int_t AliKMeansClustering::SoftKMeans2(Int_t k, Int_t n, Double_t* x, Double_t* y, Double_t* mx, Double_t* my , Double_t* sigma2, Double_t* rk )
{
    //
    // The soft K-means algorithm
    //
    Int_t i,j;
    //
    // (1) Initialisation of the k means using k-means++ recipe
    // 
     OptimalInit(k, n, x, y, mx, my);
    //
    // (2a) The responsibilities
    Double_t** r = new Double_t*[n]; // responsibilities
    for (j = 0; j < n; j++) {r[j] = new Double_t[k];}
    //
    // (2b) Normalisation
    Double_t* nr = new Double_t[n];
    //
    // (2c) Weights 
    Double_t* pi = new Double_t[k];
    //
    //
    // (2d) Initialise the responsibilties and weights
    for (j = 0; j < n; j++) {
      nr[j] = 0.;
      for (i = 0; i < k; i++) {
	r[j][i] = TMath::Exp(- fBeta * d(mx[i], my[i], x[j], y[j]));
	nr[j] += r[j][i];
      } // mean i
    } // data point j
    
    for (i = 0; i < k; i++) {
      rk[i]    = 0.;
      sigma2[i] = 1./fBeta;
 
      for (j = 0; j < n; j++) {
	r[j][i] /=  nr[j];
	rk[i] += r[j][i];
      } // mean i
      pi[i] = rk[i] / Double_t(n);
    } // data point j
    // (3) Iterations
    Int_t nit = 0;

    while(1) {
	nit++;
      //
      // Assignment step
      //
      for (j = 0; j < n; j++) {
	nr[j] = 0.;
	for (i = 0; i < k; i++) {
	  r[j][i] = pi[i] * TMath::Exp(- d(mx[i], my[i], x[j], y[j]) / sigma2[i] ) 
	    / (2. * sigma2[i] * TMath::Pi() * TMath::Pi());
	  nr[j] += r[j][i];
	} // mean i
      } // data point j
	
      for (i = 0; i < k; i++) {
	for (j = 0; j < n; j++) {
	  r[j][i] /=  nr[j];
	} // mean i
      } // data point j
      
	//
	// Update step
      Double_t di = 0;
      
      for (i = 0; i < k; i++) {
	  Double_t oldx = mx[i];
	  Double_t oldy = my[i];
	  
	  mx[i] = x[0];
	  my[i] = y[0];
	  rk[i] = r[0][i];
	for (j = 1; j < n; j++) {
	    Double_t xx =  x[j];
//
// Here we have to take into acount the cylinder topology where phi is defined mod 2xpi
// If two coordinates are separated by more than pi in phi one has to be shifted by +/- 2 pi

	    Double_t dx = mx[i] - x[j];
	    if (dx >  TMath::Pi()) xx += 2. * TMath::Pi();
	    if (dx < -TMath::Pi()) xx -= 2. * TMath::Pi();
	    if (r[j][i] > 1.e-15) {
	      mx[i] = mx[i] * rk[i] + r[j][i] * xx;
	      my[i] = my[i] * rk[i] + r[j][i] * y[j];
	      rk[i] += r[j][i];
	      mx[i] /= rk[i];
	      my[i] /= rk[i];	
	    }    
	    if (mx[i] > 2. * TMath::Pi()) mx[i] -= 2. * TMath::Pi();
	    if (mx[i] < 0.              ) mx[i] += 2. * TMath::Pi();
	} // Data
	di += d(mx[i], my[i], oldx, oldy);

      } // means 
      //
      // Sigma
      for (i = 0; i < k; i++) {
	sigma2[i] = 0.;
	for (j = 0; j < n; j++) {
	  sigma2[i] += r[j][i] * d(mx[i], my[i], x[j], y[j]);
	} // Data
	sigma2[i] /= rk[i];
	if (sigma2[i] < 0.0025) sigma2[i] = 0.0025;
      } // Clusters    
      //
      // Fractions
      for (i = 0; i < k; i++) pi[i] = rk[i] / Double_t(n);
      //
// ending condition
      if (di < 1.e-8 || nit > 1000) break;
    } // while

// Clean-up    
    delete[] nr;
    delete[] pi;
    for (j = 0; j < n; j++) delete[] r[j];
    delete[] r;
// 
    return (nit < 1000);
}

Int_t AliKMeansClustering::SoftKMeans3(Int_t k, Int_t n, Double_t* x, Double_t* y, Double_t* mx, Double_t* my , 
				       Double_t* sigmax2, Double_t* sigmay2, Double_t* rk )
{
    //
    // The soft K-means algorithm
    //
    Int_t i,j;
    //
    // (1) Initialisation of the k means using k-means++ recipe
    // 
     OptimalInit(k, n, x, y, mx, my);
    //
    // (2a) The responsibilities
    Double_t** r = new Double_t*[n]; // responsibilities
    for (j = 0; j < n; j++) {r[j] = new Double_t[k];}
    //
    // (2b) Normalisation
    Double_t* nr = new Double_t[n];
    //
    // (2c) Weights 
    Double_t* pi = new Double_t[k];
    //
    //
    // (2d) Initialise the responsibilties and weights
    for (j = 0; j < n; j++) {
      nr[j] = 0.;
      for (i = 0; i < k; i++) {

	r[j][i] = TMath::Exp(- fBeta * d(mx[i], my[i], x[j], y[j]));
	nr[j] += r[j][i];
      } // mean i
    } // data point j
    
    for (i = 0; i < k; i++) {
      rk[i]    = 0.;
      sigmax2[i] = 1./fBeta;
      sigmay2[i] = 1./fBeta;
 
      for (j = 0; j < n; j++) {
	r[j][i] /=  nr[j];
	rk[i] += r[j][i];
      } // mean i
      pi[i] = rk[i] / Double_t(n);
    } // data point j
    // (3) Iterations
    Int_t nit = 0;

    while(1) {
	nit++;
      //
      // Assignment step
      //
      for (j = 0; j < n; j++) {
	nr[j] = 0.;
	for (i = 0; i < k; i++) {

	  Double_t dx = TMath::Abs(mx[i]-x[j]);
	  if (dx > TMath::Pi()) dx = 2. * TMath::Pi() - dx;
	  Double_t dy = TMath::Abs(my[i]-y[j]);
	  r[j][i] = pi[i] * TMath::Exp(-0.5 *  (dx * dx / sigmax2[i] + dy * dy / sigmay2[i])) 
	    / (2. * TMath::Sqrt(sigmax2[i] * sigmay2[i]) * TMath::Pi() * TMath::Pi());
	  nr[j] += r[j][i];
	} // mean i
      } // data point j
	
      for (i = 0; i < k; i++) {
	for (j = 0; j < n; j++) {
	  r[j][i] /=  nr[j];
	} // mean i
      } // data point j
      
	//
	// Update step
      Double_t di = 0;
      
      for (i = 0; i < k; i++) {
	  Double_t oldx = mx[i];
	  Double_t oldy = my[i];
	  
	  mx[i] = x[0];
	  my[i] = y[0];
	  rk[i] = r[0][i];
	for (j = 1; j < n; j++) {
	    Double_t xx =  x[j];
//
// Here we have to take into acount the cylinder topology where phi is defined mod 2xpi
// If two coordinates are separated by more than pi in phi one has to be shifted by +/- 2 pi

	    Double_t dx = mx[i] - x[j];
	    if (dx >  TMath::Pi()) xx += 2. * TMath::Pi();
	    if (dx < -TMath::Pi()) xx -= 2. * TMath::Pi();
	    if (r[j][i] > 1.e-15) {
	      mx[i] = mx[i] * rk[i] + r[j][i] * xx;
	      my[i] = my[i] * rk[i] + r[j][i] * y[j];
	      rk[i] += r[j][i];
	      mx[i] /= rk[i];
	      my[i] /= rk[i];	
	    }    
	    if (mx[i] > 2. * TMath::Pi()) mx[i] -= 2. * TMath::Pi();
	    if (mx[i] < 0.              ) mx[i] += 2. * TMath::Pi();
	} // Data
	di += d(mx[i], my[i], oldx, oldy);

      } // means 
      //
      // Sigma
      for (i = 0; i < k; i++) {
	sigmax2[i] = 0.;
	sigmay2[i] = 0.;

	for (j = 0; j < n; j++) {
	  Double_t dx = TMath::Abs(mx[i]-x[j]);
	  if (dx > TMath::Pi()) dx = 2. * TMath::Pi() - dx;
	  Double_t dy = TMath::Abs(my[i]-y[j]);
	  sigmax2[i] += r[j][i] * dx * dx;
	  sigmay2[i] += r[j][i] * dy * dy;
	} // Data
	sigmax2[i] /= rk[i];
	sigmay2[i] /= rk[i];
	if (sigmax2[i] < 0.0025) sigmax2[i] = 0.0025;
	if (sigmay2[i] < 0.0025) sigmay2[i] = 0.0025;
      } // Clusters    
      //
      // Fractions
      for (i = 0; i < k; i++) pi[i] = rk[i] / Double_t(n);
      //
// ending condition
      if (di < 1.e-8 || nit > 1000) break;
    } // while

// Clean-up    
    delete[] nr;
    delete[] pi;
    for (j = 0; j < n; j++) delete[] r[j];
    delete[] r;
// 
    return (nit < 1000);
}

Double_t AliKMeansClustering::d(Double_t mx, Double_t my, Double_t x, Double_t y)
{
    //
    // Distance definition 
    // Quasi - Euclidian on the eta-phi cylinder
    
    Double_t dx = TMath::Abs(mx-x);
    if (dx > TMath::Pi()) dx = 2. * TMath::Pi() - dx;
    
    return (0.5*(dx * dx + (my - y) * (my - y)));
}



void AliKMeansClustering::OptimalInit(Int_t k, Int_t n, const Double_t* x, const Double_t* y, Double_t* mx, Double_t* my)
{
  //  
  // Optimal initialisation using the k-means++ algorithm
  // http://en.wikipedia.org/wiki/K-means%2B%2B
  //
  // k-means++ is an algorithm for choosing the initial values for k-means clustering in statistics and machine learning. 
  // It was proposed in 2007 by David Arthur and Sergei Vassilvitskii as an approximation algorithm for the NP-hard k-means problem---
  // a way of avoiding the sometimes poor clusterings found by the standard k-means algorithm.
  //
  //
  TH1F d2("d2", "", n, -0.5, Float_t(n)-0.5);
  d2.Reset();

  // (1) Chose first center as a random point among the input data.
  Int_t ir = Int_t(Float_t(n) * gRandom->Rndm());
  mx[0] = x[ir];
  my[0] = y[ir];

  // (2) Iterate
  Int_t icl = 1;
  while(icl < k)
    {
      // find min distance to existing clusters
      for (Int_t j = 0; j < n; j++) {
	Double_t dmin = 1.e10;
	for (Int_t i = 0; i < icl; i++) {
	  Double_t dij = d(mx[i], my[i], x[j], y[j]);
	  if (dij < dmin) dmin = dij;
	} // clusters
	d2.Fill(Float_t(j), dmin);
      } // data points
      // select a new cluster from data points with probability ~d2
      ir = Int_t(d2.GetRandom() + 0.5);
      mx[icl] = x[ir];
      my[icl] = y[ir];
      icl++;
    } // icl
}


ClassImp(AliKMeansResult)


    
AliKMeansResult::AliKMeansResult(Int_t k):
    TObject(),
    fK(k),
    fMx    (new Double_t[k]),
    fMy    (new Double_t[k]),
    fSigma2(new Double_t[k]),
    fRk    (new Double_t[k]),
    fTarget(new Double_t[k]),
    fInd   (new Int_t[k])
{
// Constructor
}

AliKMeansResult::AliKMeansResult(const AliKMeansResult &res):
  TObject(res),
  fK(res.GetK()),
  fMx(new Double_t[res.GetK()]),
  fMy(new Double_t[res.GetK()]),
  fSigma2(new Double_t[res.GetK()]),
  fRk(new Double_t[res.GetK()]),
  fTarget(new Double_t[res.GetK()]),
  fInd(new Int_t[res.GetK()])
{
  // Copy constructor
  for (Int_t i = 0; i <fK; i++) {
    fMx[i]     = (res.GetMx())    [i];
    fMy[i]     = (res.GetMy())    [i];
    fSigma2[i] = (res.GetSigma2())[i];
    fRk[i]     = (res.GetRk())    [i];
    fTarget[i] = (res.GetTarget())[i];
    fInd[i]    = (res.GetInd())   [i];
  }
}

AliKMeansResult& AliKMeansResult::operator=(const AliKMeansResult& res)
{
  //
  // Assignment operator
  if (this != &res) {
    fK = res.GetK();
    for (Int_t i = 0; i <fK; i++) {
      fMx[i]     = (res.GetMx())    [i];
      fMy[i]     = (res.GetMy())    [i];
      fSigma2[i] = (res.GetSigma2())[i];
      fRk[i]     = (res.GetRk())    [i];
      fTarget[i] = (res.GetTarget())[i];
      fInd[i]    = (res.GetInd())   [i];
    }
  }
  return *this;
}


AliKMeansResult::~AliKMeansResult()
{
// Destructor
    delete[] fMx;
    delete[] fMy;
    delete[] fSigma2;
    delete[] fRk;
    delete[] fInd;
    delete[] fTarget;
}

void AliKMeansResult::Sort()
{
  // Build target array and sort
  // Sort clusters
  for (Int_t i = 0; i < fK; i++) {
    if (fRk[i] > 2.9) {
      fTarget[i] = fRk[i] / fSigma2[i];
    }
    else fTarget[i] = 0.;
  }
    
  TMath::Sort(fK, fTarget, fInd);
}

void AliKMeansResult::Sort(Int_t n, const Double_t* x, const Double_t* y)
{
  // Build target array and sort
  for (Int_t i = 0; i < fK; i++)
    {
      Int_t nc = 0;
      for (Int_t j = 0; j < n; j++)
	{
	  if (2. * AliKMeansClustering::d(fMx[i], fMy[i], x[j], y[j])  <  2.28 * fSigma2[i]) nc++;
	}

      if (nc > 2) {
	fTarget[i] = Double_t(nc) / (2.28 * fSigma2[i]);
      } else {
	fTarget[i] = 0.;
      }
    }

  TMath::Sort(fK, fTarget, fInd);
}

void AliKMeansResult::CopyResults(const AliKMeansResult* res)
{
  fK = res->GetK();
  for (Int_t i = 0; i <fK; i++) {
    fMx[i]     = (res->GetMx())    [i];
    fMy[i]     = (res->GetMy())    [i];
    fSigma2[i] = (res->GetSigma2())[i];
    fRk[i]     = (res->GetRk())    [i];
    fTarget[i] = (res->GetTarget())[i];
    fInd[i]    = (res->GetInd())   [i];
  }
}
