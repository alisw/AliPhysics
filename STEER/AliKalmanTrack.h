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
  AliKalmanTrack() {fLab=-3141593; fChi2=0; fN=0;}
  AliKalmanTrack(const AliKalmanTrack &t) {fLab=t.fLab;fChi2=t.fChi2;fN=t.fN;}
  virtual ~AliKalmanTrack(){};
  void SetLabel(Int_t lab) {fLab=lab;}

  Bool_t   IsSortable() const {return kTRUE;}
  Int_t    GetLabel()   const {return fLab;}
  Double_t GetChi2()    const {return fChi2;}
  Int_t    GetNumberOfClusters() const {return fN;}
  virtual Int_t GetClusterIndex(Int_t i) const { //reserved for AliTracker
    printf("AliKalmanTrack::GetClusterIndex(Int_t i) must be overloaded !\n");
    return 0;
  } 

  virtual Int_t Compare(const TObject *o) const=0; 

  virtual void GetExternalParameters(Double_t &xr, Double_t x[5]) const=0;
  virtual void GetExternalCovariance(Double_t cov[15]) const=0;

  virtual Double_t GetPredictedChi2(const AliCluster *cluster) const=0;
  virtual 
  Int_t PropagateTo(Double_t xr,Double_t x0,Double_t rho,Double_t pm)=0;
  virtual Int_t Update(const AliCluster* c, Double_t chi2, UInt_t i)=0;

protected:
  void SetChi2(Double_t chi2) {fChi2=chi2;} 
  void SetNumberOfClusters(Int_t n) {fN=n;} 

private: 
  Int_t fLab;             // track label
  Double_t fChi2;         // total chi2 value for this track
  Int_t fN;               // number of associated clusters

  ClassDef(AliKalmanTrack,1)    // Reconstructed track
};

#endif


