#ifndef ALIITSTRACKV2_H
#define ALIITSTRACKV2_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//-------------------------------------------------------------------------
//                       ITS Track Class
//
//        Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch 
//     dEdx analysis by: Boris Batyunya, JINR, Boris.Batiounia@cern.ch
//-------------------------------------------------------------------------

#include <AliKalmanTrack.h>
#include "AliITSRecoParam.h"
#include "AliITSgeomTGeo.h"

class AliESDtrack;

//_____________________________________________________________________________
class AliITStrackV2 : public AliKalmanTrack {
public:
  AliITStrackV2();
  AliITStrackV2(AliESDtrack& t,Bool_t c=kFALSE) throw (const Char_t *);
  AliITStrackV2(const AliITStrackV2& t);

  Bool_t CorrectForMeanMaterial(Double_t xOverX0, Double_t xTimesRho,
				Bool_t anglecorr=kFALSE) {
    return AliExternalTrackParam::CorrectForMeanMaterial(xOverX0,xTimesRho,GetMass(),anglecorr);
  }
  Bool_t CorrectForMaterial(Double_t d, Double_t x0=AliITSRecoParam::GetX0Air()) {
    // deprecated: use CorrectForMeanMaterial instead
    return AliExternalTrackParam::CorrectForMaterial(d,x0,GetMass());
  }
  Bool_t PropagateTo(Double_t xr, Double_t d, Double_t x0=AliITSRecoParam::GetX0Air());
  Bool_t PropagateToTGeo(Double_t xToGo, Int_t nstep, Double_t &xOverX0, Double_t &xTimesRho);
  Bool_t PropagateToTGeo(Double_t xToGo, Int_t nstep=1) {
    Double_t dummy1,dummy2; return PropagateToTGeo(xToGo,nstep,dummy1,dummy2);
  }
  Double_t GetPredictedChi2(const AliCluster *cluster) const;
  Bool_t Update(const AliCluster *cl, Double_t chi2, Int_t i);

  Bool_t PropagateToVertex(const AliESDVertex *v,Double_t d=0.,Double_t x0=0.);
  Bool_t Propagate(Double_t alpha, Double_t xr);
  Bool_t MeanBudgetToPrimVertex(Double_t xyz[3], Double_t step, Double_t &d) const;
  Bool_t Improve(Double_t x0,Double_t xyz[3],Double_t ers[3]);

  void SetdEdx(Double_t dedx) {fdEdx=dedx;}
  void SetSampledEdx(Float_t q, Int_t i);
  Float_t GetSampledEdx(Int_t i) const {return fdEdxSample[i];}
  void CookdEdx(Double_t low=0., Double_t up=0.51);
  void SetDetectorIndex(Int_t i) {SetLabel(i);}
  void ResetClusters() { SetChi2(0.); SetNumberOfClusters(0); }
  void UpdateESDtrack(ULong_t flags) const;
  
  AliESDtrack *GetESDtrack() const {return fESDtrack;}

  Int_t GetDetectorIndex() const {return GetLabel();}
  Double_t GetdEdx() const {return fdEdx;}
  Double_t GetPIDsignal() const {return GetdEdx();}
  Double_t GetC() const {return AliExternalTrackParam::GetC(GetBz());}
  Double_t GetD(Double_t x, Double_t y) const {
    return AliExternalTrackParam::GetD(x,y,GetBz());
  }
  void GetDZ(Double_t xv, Double_t yv, Double_t zv, Float_t dz[2]) const {
    return AliExternalTrackParam::GetDZ(xv,yv,zv,GetBz(),dz);
  }

  Bool_t GetGlobalXYZat(Double_t xloc,Double_t &x,Double_t &y,Double_t &z) const;
  Bool_t GetPhiZat(Double_t r,Double_t &phi,Double_t &z) const;
  Bool_t GetLocalXat(Double_t r,Double_t &xloc) const;

  Int_t Compare(const TObject *o) const;
  Int_t GetClusterIndex(Int_t i) const {return fIndex[i];}
  Bool_t Invariant() const;

  void  SetExtraCluster(Int_t i, Int_t idx) {fIndex[AliITSgeomTGeo::kNLayers+i]=idx;}
  Int_t GetExtraCluster(Int_t i) const {return fIndex[AliITSgeomTGeo::kNLayers+i];}

protected:
  Double_t GetBz() const ;
  Double_t fdEdx;            // dE/dx

  static const Int_t fgkWARN; //! used for debugging purposes
  Float_t fdEdxSample[4];    // array of dE/dx samples b.b.

  Int_t fIndex[2*AliITSgeomTGeo::kNLayers]; // indices of associated clusters 

  AliESDtrack *fESDtrack;    //! pointer to the connected ESD track

private:
  AliITStrackV2 &operator=(const AliITStrackV2 &tr);
  ClassDef(AliITStrackV2,7)  //ITS reconstructed track
};

inline
void AliITStrackV2::SetSampledEdx(Float_t q, Int_t i) {
  //----------------------------------------------------------------------
  // This function stores dEdx sample corrected for the track segment length 
  // Origin: Boris Batyunya, JINR, Boris.Batiounia@cern.ch
  //----------------------------------------------------------------------
  if (i<0) return;
  if (i>3) return;
  Double_t s=GetSnp(), t=GetTgl();
  q *= TMath::Sqrt((1-s*s)/(1+t*t));
  fdEdxSample[i]=q;
}

#endif


