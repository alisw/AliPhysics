#ifndef ALIKALMANTRACK_H
#define ALIKALMANTRACK_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//-------------------------------------------------------------------------
//                          Class AliKalmanTrack
//      fixed the interface for the derived reconstructed track classes 
//            Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch 
//-------------------------------------------------------------------------

#include <TObject.h>
#include "AliLog.h"
#include "AliPID.h"

class AliCluster;

class AliKalmanTrack : public TObject {
public:
  AliKalmanTrack();
  AliKalmanTrack(const AliKalmanTrack &t);

  virtual ~AliKalmanTrack(){};
  void SetLabel(Int_t lab) {fLab=lab;}
  void SetFakeRatio(Float_t ratio) {fFakeRatio=ratio;}

  Bool_t   IsSortable() const {return kTRUE;}
  Int_t    GetLabel()   const {return fLab;}
  Float_t    GetFakeRatio()   const {return fFakeRatio;}
  Double_t GetChi2()    const {return fChi2;}
  Double_t GetMass()    const {return fMass;}
  Int_t    GetNumberOfClusters() const {return fN;}
  virtual Int_t GetClusterIndex(Int_t) const { //reserved for AliTracker
    AliWarning("Method must be overloaded !\n");
    return 0;
  } 
  virtual Double_t GetPIDsignal() const {
    AliWarning("Method must be overloaded !\n");
    return 0.;
  }

  virtual Double_t GetAlpha() const {
    AliWarning("Method must be overloaded !\n");
    return 0.;
  }
  virtual Double_t GetSigmaY2() const {
    AliWarning("Method must be overloaded !\n");
    return 0.;
  }
  virtual Double_t GetSigmaZ2() const {
    AliWarning("Method must be overloaded !\n");
    return 0.;
  }

  virtual Int_t Compare(const TObject *) const {return 0;} 

  virtual void GetExternalParameters(Double_t&/*xr*/,Double_t/*x*/[5]) const=0;
  virtual void GetExternalCovariance(Double_t /*cov*/[15]) const = 0;

  virtual Double_t GetPredictedChi2(const AliCluster *) const = 0;
  virtual Int_t PropagateTo(Double_t/*xr*/,Double_t/*x0*/,Double_t/*rho*/) = 0;

  virtual Int_t Update(const AliCluster*, Double_t /*chi2*/, UInt_t) = 0;

  static Double_t MeanMaterialBudget(Double_t *start, Double_t *end, Double_t *mparam);
 
  // Time integration (S.Radomski@gsi.de)
  void   StartTimeIntegral();
  void SetIntegratedLength(Double_t l) {fIntegratedLength=l;}
  void SetIntegratedTimes(const Double_t *times);

  Bool_t IsStartedTimeIntegral() const {return fStartTimeIntegral;}
  void     AddTimeStep(Double_t length);
  void GetIntegratedTimes(Double_t *times) const;
  Double_t GetIntegratedTime(Int_t pdg) const;
  Double_t GetIntegratedLength() const {return fIntegratedLength;}
  virtual void GetXYZ(Float_t r[3]) const = 0;

  static Double_t GetConvConst();
  void SetNumberOfClusters(Int_t n) {fN=n;} 
protected:
  void     SaveLocalConvConst();
  Double_t GetLocalConvConst() const {return fLocalConvConst;}

  void External2Helix(Double_t helix[6]) const;

  void SetChi2(Double_t chi2) {fChi2=chi2;} 
  void SetMass(Double_t mass) {fMass=mass;}

  Int_t fLab;             // track label
  Float_t fFakeRatio;     // fake ratio
  Double_t fChi2;         // total chi2 value for this track
  Double_t fMass;         // mass hypothesis
  Int_t fN;               // number of associated clusters

private:
  Double_t fLocalConvConst;   //local conversion "curvature(1/cm) -> pt(GeV/c)"

  // variables for time integration (S.Radomski@gsi.de)
  Bool_t  fStartTimeIntegral;       // indicator wether integrate time
  Double_t fIntegratedTime[AliPID::kSPECIES];       // integrated time
  Double_t fIntegratedLength;        // integrated length
  
  ClassDef(AliKalmanTrack,5)    // Reconstructed track
};

#endif


