#ifndef ALITPCTRACK_H
#define ALITPCTRACK_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
/* $Id$ */

//-------------------------------------------------------
//                    TPC Track Class
//
//   Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch 
//
// The track parameterization is fixed in the following way:                   
//      param0:   local Y-coordinate of a track (cm)                
//      param1:   local Z-coordinate of a track (cm)                
//      param2:   local sine of the track momentum azimuth angle    
//      param3:   tangent of the track momentum dip angle           
//      param4:   1/pt (1/(GeV/c))                                  
//
//-------------------------------------------------------

#include <AliKalmanTrack.h>
#include <TMath.h>

#include "AliTPCreco.h"
#include "AliExternalTrackParam.h"
class AliESDtrack;
class AliESDVertex;

//_____________________________________________________________________________
class AliTPCtrack : public AliKalmanTrack {
public:
  AliTPCtrack();
  AliTPCtrack(Double_t x, Double_t alpha, const Double_t p[5], 
              const Double_t cov[15], Int_t index); 
  AliTPCtrack(const AliESDtrack& t);
  AliTPCtrack(const AliTPCtrack& t);
  virtual ~AliTPCtrack() {}

  Int_t Compare(const TObject *o) const;

  void SetdEdx(Double_t dedx) {fdEdx=dedx;}
  Double_t GetdEdx()  const {return fdEdx;}
  Double_t GetPIDsignal()  const {return GetdEdx();}

  Int_t GetClusterIndex(Int_t i) const {return fIndex[i];}
  void  SetClusterIndex(Int_t i, Int_t idx) {fIndex[i]=idx;}

  Double_t GetC() const {return AliExternalTrackParam::GetC(GetBz());}

  Double_t GetPredictedChi2(const AliCluster *cluster) const;
  Bool_t PropagateTo(Double_t xr, Double_t rho=0.9e-3, Double_t x0=28.94);
  Bool_t Update(const AliCluster *c, Double_t chi2, Int_t i);
  Bool_t Rotate(Double_t alpha) {
    return AliExternalTrackParam::Rotate(GetAlpha()+alpha);
  }
  
  Bool_t PropagateToVertex(const AliESDVertex *v, 
                           Double_t rho=1.2e-3, Double_t x0=36.66);
  void ResetClusters() {SetNumberOfClusters(0); SetChi2(0.);}
  void UpdatePoints();//update points 
  Float_t* GetPoints() {return fPoints;}

  Float_t Density(Int_t row0, Int_t row1); //calculate cluster density
  Float_t Density2(Int_t row0, Int_t row1); //calculate cluster density
  Double_t GetD(Double_t x=0, Double_t y=0) const {
     return AliExternalTrackParam::GetD(x,y,GetBz());
  }
  AliExternalTrackParam & GetReference(){ return fReference;}
  void UpdateReference(){ new (&fReference) AliExternalTrackParam(*this);}
  Int_t   GetKinkIndex(Int_t i) const{ return fKinkIndexes[i];}
  Int_t*  GetKinkIndexes() { return fKinkIndexes;}
  Int_t   GetV0Index(Int_t i) const{ return fV0Indexes[i];}
  Int_t*  GetV0Indexes() { return fV0Indexes;}

  void SetFirstPoint(Int_t f) {fFirstPoint=f;}
  void SetLastPoint(Int_t f) {fLastPoint=f;}
  void SetRemoval(Int_t f) {fRemoval=f;}
  void SetLab2(Int_t f) {fLab2=f;}
  void SetKinkIndex(Int_t i, Int_t idx) {fKinkIndexes[i]=idx;}
  void SetBConstrain(Bool_t b) {fBConstrain=b;}
  void SetNFoundable(Int_t n) {fNFoundable=n;}
  void SetNShared(Int_t s) {fNShared=s;}

  Int_t GetFirstPoint() const {return fFirstPoint;}
  Int_t GetLastPoint() const {return fLastPoint;}
  Int_t GetRemoval() const {return fRemoval;}
  Int_t GetLab2() const {return fLab2;}
  Bool_t GetBConstrain() const {return fBConstrain;}
  Int_t GetNShared() const {return fNShared;}
  Int_t GetNFoundable() const {return fNFoundable;}

protected: 
  Double_t GetBz() const;
  Double_t fdEdx;           // dE/dx

  Int_t fIndex[kMaxRow];       // indices of associated clusters 
  Float_t fPoints[4];            //first, max dens row  end points of the track and max density
  // MI addition
  Float_t fSdEdx;           // sigma of dedx 
  //
  Int_t   fNFoundable;      //number of foundable clusters - dead zone taken to the account
  Bool_t  fBConstrain;   // indicate seeding with vertex constrain
  Int_t   fLastPoint;     // last  cluster position     
  Int_t   fFirstPoint;    // first cluster position
  Int_t fRemoval;         // removal factor
  Int_t fTrackType;       // track type - 0 - normal - 1 - kink -  2 -V0  3- double found
  Int_t fLab2;            // index of corresponding track (kink, V0, double)
  Int_t   fNShared;       // number of shared points 
  AliExternalTrackParam   fReference; // track parameters at the middle of the chamber
  Float_t  fKinkPoint[12];      //radius, of kink,  dfi and dtheta
  Int_t    fKinkIndexes[3];     // kink indexes - minus = mother + daughter
  Int_t    fV0Indexes[3];     // kink indexes - minus = mother + daughter

  ClassDef(AliTPCtrack,3)   // Time Projection Chamber reconstructed tracks
};

#endif

