//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTGLOBALBARRELTRACK_H
#define ALIHLTGLOBALBARRELTRACK_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTGlobalBarrelTrack.h
    @author Matthias Richter
    @date   2009-06-24
    @brief  An AliKalmanTrack implementation for global HLT barrel tracks.
*/

#include "AliKalmanTrack.h"
#include "AliHLTDataTypes.h"
#include "AliHLTExternalTrackParam.h"
#include <vector>
using namespace std;

/**
 * @class AliHLTGlobalBarrelTrack
 * Representation of global HLT barrel tracks.
 *
 * @ingroup alihlt_global_components
 */
class AliHLTGlobalBarrelTrack : public AliKalmanTrack
{
 public:
  /** standard constructor */
  AliHLTGlobalBarrelTrack();
  /** copy constructor */
  AliHLTGlobalBarrelTrack(const AliHLTGlobalBarrelTrack& t);
  /** copy constructor */
  AliHLTGlobalBarrelTrack(const AliHLTExternalTrackParam& p);

  /** assignment operator */
  template <class c>
  AliHLTGlobalBarrelTrack& operator=(const c& t);
  /** destructor */
  ~AliHLTGlobalBarrelTrack();

  /// inherited from AliKalmanTrack
  Int_t GetClusterIndex(Int_t i) const { 
    return (i<(int)fPoints.size()) ?fPoints[i] :0;
  } 

  /// inherited from AliKalmanTrack, dummy implementation
  virtual Int_t GetNumberOfTracklets() const {return 0;}
  /// inherited from AliKalmanTrack, dummy implementation
  virtual Int_t GetTrackletIndex(Int_t) const {return -1;} 
  /// inherited from AliKalmanTrack, dummy implementation
  virtual Double_t GetPIDsignal() const {return 0.;}

  /// Get the x position of the last assigned point
  Double_t GetLastPointX() const {return fLastX;}
  /// Get the y position of the last assigned point
  Double_t GetLastPointY() const {return fLastY;}
  /// return Track ID
  Int_t TrackID() const {return fTrackID;}
  /// return Track ID, inherited for AliExternalTrackParam
  Int_t GetID() const {return fTrackID;}

  /// Get the number of associated points
  UInt_t GetNumberOfPoints() const;

  /// Get the list of associated points
  const UInt_t* GetPoints() const;

  /// Set the list of associated points
  int SetPoints(const UInt_t* pArray, UInt_t arraySize);

  static int ConvertTrackDataArray(const AliHLTTracksData* pTracks, unsigned sizeInByte, vector<AliHLTGlobalBarrelTrack> &tgtArray);

  /// inherited from AliKalmanTrack, dummy implementation
  Double_t GetPredictedChi2(const AliCluster*) const {return 0.0;}

  /// inherited from AliKalmanTrack, dummy implementation
  Bool_t PropagateTo(Double_t, Double_t, Double_t) {return kFALSE;}

  /// inherited from AliKalmanTrack, dummy implementation
  Bool_t Update(const AliCluster*, Double_t, Int_t) {return kFALSE;}

  /// Inherited from TObject, prints the track parameters
  virtual void Print(Option_t* option = "") const;

  Double_t GetPathLengthTo( Double_t x, Double_t b ) const;

 protected:

 private:

  /// array of points
  vector<UInt_t> fPoints; //

  /// x position of the last assigned point
  Double_t fLastX; //
  /// y position of the last assigned point
  Double_t fLastY; //

  /// track Id for identification during the reconstruction
  Int_t   fTrackID; //

  ClassDef(AliHLTGlobalBarrelTrack, 0)
};
#endif
