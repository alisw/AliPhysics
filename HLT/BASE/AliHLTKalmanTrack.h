//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

#ifndef _ALIHLTKALMANTRACK_H_
#define _ALIHLTKALMANTRACK_H_

#include "AliKalmanTrack.h"

class AliHLTExternalTrackParam;

/**
 * @class AliHLTKalmanTrack
 * The class is only used for copy AliHLTExternalTrackParam tracks to esd tracks
 */
class AliHLTKalmanTrack :public AliKalmanTrack
{
 public:
  AliHLTKalmanTrack(): AliKalmanTrack(){}
  AliHLTKalmanTrack( const AliHLTExternalTrackParam &t);
  
  Double_t GetPredictedChi2(const AliCluster * /*c*/) const { return 0; }
  Bool_t PropagateTo(Double_t /*xr*/, Double_t /*x0*/, Double_t /*rho*/){ return 0; }
  Bool_t Update(const AliCluster* /*c*/, Double_t /*chi2*/, Int_t /*index*/){ return 0; }
};


#endif
