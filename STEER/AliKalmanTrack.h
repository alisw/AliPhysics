#ifndef ALIKALMANTRACK_H
#define ALIKALMANTRACK_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//-------------------------------------------------------------------------
//                          Class AliKalmanTrack
//
//         Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch 
//-------------------------------------------------------------------------

#include <TObject.h>

class AliCluster;

class AliKalmanTrack : public TObject {
public:
  AliKalmanTrack() {fN=0; fChi2=0; fLab=-3141593;}
  AliKalmanTrack(const AliKalmanTrack& t);
  virtual ~AliKalmanTrack() {}
  Int_t Compare(const TObject *o) const;
  void SetLabel(Int_t lab) {fLab=lab;} 

  Double_t GetPredictedChi2(const AliCluster *cluster) const;
  Bool_t IsSortable() const {return kTRUE;}
  Int_t GetLabel() const {return fLab;}
  void GetCovariance(Double_t cov[15]) const;
  Double_t GetChi2() const {return fChi2;}
  Int_t GetNumberOfClusters() const {return fN;}

 virtual Double_t GetPt() const=0;
 virtual Double_t GetP()  const=0;
 virtual void GetPxPyPz(Double_t &px, Double_t &py, Double_t &pz) const=0;
 virtual Int_t PropagateTo(Double_t xr,Double_t x0,Double_t rho,Double_t pm)=0;
 virtual void Update(const AliCluster* c, Double_t chi2, UInt_t i)=0;

protected: 
  Int_t fLab;             // track label

  Double_t fP0;           // track parameter
  Double_t fP1;           // track parameter
  Double_t fP2;           // track parameter
  Double_t fP3;           // track parameter
  Double_t fP4;           // track parameter

  Double_t fC00;                         // covariance
  Double_t fC10, fC11;                   // matrix
  Double_t fC20, fC21, fC22;             // of the
  Double_t fC30, fC31, fC32, fC33;       // track
  Double_t fC40, fC41, fC42, fC43, fC44; // parameters

  Double_t fChi2;         // total chi2 value for this track
  Short_t fN;             // number of clusters 

  ClassDef(AliKalmanTrack,1)    // Reconstructed track
};

#endif


