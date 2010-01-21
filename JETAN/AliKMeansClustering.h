#ifndef ALIKMEANSCLUSTERING_H
#define ALIKMEANSCLUSTERING_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Implemenatation of the K-Means Clustering Algorithm
// http://en.wikipedia.org/wiki/K-means_clustering
// This particular implementation is the so called Soft K-means algorithm.
// It has been modified to work on the cylindrical topology in eta-phi space.
//
// Author: Andreas Morsch (CERN)
// andreas.morsch@cern.ch

#include <TObject.h>
 
class AliKMeansClustering : public TObject
{
 public:
  AliKMeansClustering()          {}
  virtual ~AliKMeansClustering() {}
  
  static void SoftKMeans(Int_t k, Int_t n, Double_t* x, Double_t* y, Double_t* mx, Double_t* my , Double_t* rk );
  static void SetBeta(Double_t beta) {fBeta = beta;}
  static Double_t d(Double_t mx, Double_t my, Double_t x, Double_t y);
protected:
  static Double_t fBeta;
  
  ClassDef(AliKMeansClustering, 1)
};
 
#endif
