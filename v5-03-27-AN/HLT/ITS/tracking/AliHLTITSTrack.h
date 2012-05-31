#ifndef ALIHLTITSTRACK1_H
#define ALIHLTITSTRACK1_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


#include "AliKalmanTrack.h"
#include "AliITSRecoParam.h"
#include "AliITSgeomTGeo.h"

class AliESDtrack;
class AliESDVertex;
class AliTracker;
class AliESDtrack;

//_____________________________________________________________________________
class AliHLTITSTrack : public AliKalmanTrack 
{
 public:

  AliHLTITSTrack();
  AliHLTITSTrack(AliESDtrack& t,Bool_t c=kFALSE) throw (const Char_t *);
  AliHLTITSTrack(AliExternalTrackParam& t ) throw (const Char_t *);
  AliHLTITSTrack(const AliHLTITSTrack& t);
  AliHLTITSTrack &operator=(const AliHLTITSTrack& t);

  Int_t TPCtrackId(){ return fTPCtrackId;}
  void SetTPCtrackId( Int_t v ){ fTPCtrackId = v;}

  Int_t GetProlongationFast(Double_t alpha, Double_t xr,Double_t &y, Double_t &z);


  Float_t GetExpQ() const {return fExpQ;}
  void SetExpQ(Float_t f) {fExpQ=f;}

  Double_t GetPredictedChi2(const AliCluster* c) const;
  Double_t GetPredictedChi2(Double_t cy, Double_t cz, Double_t cerr2Y, Double_t cerr2Z) const;


  Bool_t CorrectForMeanMaterial(Double_t xOverX0, Double_t xTimesRho,
				Bool_t anglecorr=kFALSE) {
    return AliExternalTrackParam::CorrectForMeanMaterial(xOverX0,xTimesRho,GetMass(),anglecorr);
  }

  Bool_t PropagateTo(Double_t xr, Double_t d, Double_t x0=AliITSRecoParam::GetX0Air());

  Bool_t PropagateToTGeo(Double_t xToGo, Int_t nstep, Double_t &xOverX0, Double_t &xTimesRho, Bool_t addTime=kTRUE);
  Bool_t PropagateToTGeo(Double_t xToGo, Int_t nstep=1, Bool_t addTime=kTRUE) {
    Double_t dummy1,dummy2; return PropagateToTGeo(xToGo,nstep,dummy1,dummy2,addTime);
  }

  Bool_t Update(const AliCluster *cl, Double_t chi2, Int_t i);
  
  Bool_t Propagate(Double_t alpha, Double_t xr);
  Bool_t Propagate(Double_t xr) { return Propagate(GetAlpha(),xr); }
  
  void ResetClusters();
  
  using AliExternalTrackParam::GetC;
  Double_t GetC() const {return AliExternalTrackParam::GetC(GetBz());}
  Double_t GetD(Double_t x, Double_t y) const {
    return AliExternalTrackParam::GetD(x,y,GetBz());
  }
 
  Bool_t GetGlobalXYZat(Double_t xloc,Double_t &x,Double_t &y,Double_t &z) const;
  Bool_t GetPhiZat(Double_t r,Double_t &phi,Double_t &z) const;
  //Bool_t GetLocalXat(Double_t r,Double_t &xloc) const;
  Bool_t GetLocalXPhiZat(Double_t r,Double_t &xloc, double &phi, double &z ) const ;

  Int_t GetClusterIndex(Int_t i) const {return fIndex[i];}
  void SetClusterIndex(Int_t i, Int_t index ) { fIndex[i] = index;}

  Bool_t GetYZAtPhiX( double phi, double x,
		      double &y, double&z, double &snp, double cov[3] ) const;

  Bool_t GetLocalYZat(Double_t xloc, Double_t &y, Double_t &z) const ;

 protected:


  Int_t fIndex[2*AliITSgeomTGeo::kNLayers]; // indices of associated clusters 

  Float_t fExpQ;            // expected Q
  Int_t  fTPCtrackId;

  ClassDef(AliHLTITSTrack,0)   //HLT ITS tracker
};





#endif


