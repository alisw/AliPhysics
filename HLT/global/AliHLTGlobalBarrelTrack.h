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
  /** assignment operator */
  AliHLTGlobalBarrelTrack& operator=(const AliHLTGlobalBarrelTrack& t);
  /** assignment operator */
  AliHLTGlobalBarrelTrack& operator=(const AliHLTExternalTrackParam& extp);
  /** destructor */
  ~AliHLTGlobalBarrelTrack();

  /// Get the x position of the last assigned point
  Double_t GetLastPointX() const {return fLastX;}
  /// Get the y position of the last assigned point
  Double_t GetLastPointY() const {return fLastY;}

  /// Get the number of associated points
  UInt_t GetNumberOfPoints() const;

  /// Get the list of associated points
  const UInt_t* GetPoints() const;

  /// Set the list of associated points
  int SetPoints(const UInt_t* pArray, UInt_t arraySize);

  static int ConvertTrackDataArray(const AliHLTTracksData* pTracks, int sizeInByte, vector<AliHLTGlobalBarrelTrack> &tgtArray);

  /// dummy function required by AliKalmanTrack
  Double_t GetPredictedChi2(const AliCluster*) const {return 0.0;}

  /// dummy function required by AliKalmanTrack
  Bool_t PropagateTo(Double_t, Double_t, Double_t) {return kFALSE;}

  /// dummy function required by AliKalmanTrack
  Bool_t Update(const AliCluster*, Double_t, Int_t) {return kFALSE;}

 protected:

 private:
  /// array of points
  vector<UInt_t> fPoints; //

  /// x position of the last assigned point
  Double_t fLastX; //
  /// y position of the last assigned point
  Double_t fLastY; //

  ClassDef(AliHLTGlobalBarrelTrack, 0)
};
#endif
