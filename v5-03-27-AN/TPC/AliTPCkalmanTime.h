#ifndef ALITPCKALMANTIME_H
#define ALITPCKALMANTIME_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "TNamed.h"
#include "TMatrixD.h"
class TTreeSRedirector;

class AliTPCkalmanTime: public TNamed{
public:
  AliTPCkalmanTime();
  AliTPCkalmanTime(Double_t time, Double_t xoff, Double_t k, Double_t sigmaxoff, Double_t sigmak);
  void Propagate(Double_t time, Double_t sigma,  TTreeSRedirector *debug=0);
  void Update(Double_t x, Double_t xerr, Double_t ptratio, TTreeSRedirector *debug=0);
  static void TestMC(const char * fname);
public:
  void Init(Double_t time, Double_t xoff, Double_t k, Double_t sigmaxoff, Double_t sigmak);
  TMatrixD * fState;           // state vector
  TMatrixD * fCovariance;      // covariance
  Double_t   fTime;            // current time
private:
  AliTPCkalmanTime&  operator=(const AliTPCkalmanTime&);// not implemented
  AliTPCkalmanTime(const AliTPCkalmanTime&); //not implemented
  ClassDef(AliTPCkalmanTime,1);
};

#endif

