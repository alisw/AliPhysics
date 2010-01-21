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

ClassImp(AliKMeansClustering)

Double_t AliKMeansClustering::fBeta = 10.;
 
void AliKMeansClustering::SoftKMeans(Int_t k, Int_t n, Double_t* x, Double_t* y, Double_t* mx, Double_t* my , Double_t* rk )
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
    Double_t* nr = new Double_t[n];

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
}

Double_t AliKMeansClustering::d(Double_t mx, Double_t my, Double_t x, Double_t y)
{
    //
    // Distance definition 
    // Quasi - Euclidian on the eta-phi cylinder
    
    Double_t dx = TMath::Abs(mx-x);
    if (dx > TMath::Pi()) dx = 2. * TMath::Pi() - dx;
    
    return (dx * dx + (my - y) * (my - y));
}
