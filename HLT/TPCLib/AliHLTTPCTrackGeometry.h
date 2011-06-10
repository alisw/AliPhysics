//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTTPCTRACKGEOMETRY_H
#define ALIHLTTPCTRACKGEOMETRY_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/// @file   AliHLTTPCTrackGeometry.h
/// @author Matthias Richter
/// @date   2011-05-20
/// @brief  Desciption of a track by a sequence of track points
///

#include "AliHLTTrackGeometry.h"

class AliHLTGlobalBarrelTrack;

/**
 * @class AliHLTTPCTrackGeometry
 * Description of tracks in the TPC geometry.
 * This implementation describes track points of the track crossing with
 * the pad-row plane (yz). The radial distance of the plane is the x coordinate.
 * 
 * The 32bit ids of the planes follow the same coding as the TPC clusters.
 * - bit 25-31:  7 bit slice number
 * - bit 22-24:  3 bit partition number
 * - bit  0-21: 22 bit row number
 * @ingroup alihlt-tpc
 */
class AliHLTTPCTrackGeometry : public AliHLTTrackGeometry
{
 public:
  /// standard constructor
  AliHLTTPCTrackGeometry();
  /// copy constructor
  AliHLTTPCTrackGeometry(const AliHLTTPCTrackGeometry&);
  /// assignment operator
  AliHLTTPCTrackGeometry& operator=(const AliHLTTPCTrackGeometry&);

  /// destructor
  ~AliHLTTPCTrackGeometry();

  /// alpha of the plane
  virtual float GetPlaneAlpha(AliHLTUInt32_t planeId) const;
  /// radial distance from global {0,0,0}
  virtual float GetPlaneR(AliHLTUInt32_t planeId) const;
  /// theta of the plane
  virtual float GetPlaneTheta(AliHLTUInt32_t planeId) const ;
  /// check bounds in u and v coordinate
  virtual bool CheckBounds(AliHLTUInt32_t planeId, float u, float v) const;

  // track interface

  /// calculate the track points, expects the global magnetic field to be initialized
  virtual int CalculateTrackPoints(const AliHLTExternalTrackParam& track);

  /// calculate the track points, expects the global magnetic field to be initialized
  virtual int CalculateTrackPoints(AliHLTGlobalBarrelTrack& track);

  /// find the track point which can be associated to a spacepoint with coordinates and id
  virtual int FindMatchingTrackPoint(AliHLTUInt32_t spacepointId, float spacepoint[3], AliHLTUInt32_t& planeId);

 private:
  /// calculate the track points, expects the global magnetic field to be initialized
  int CalculateTrackPoints(AliHLTGlobalBarrelTrack& track, int firstpadrow, int step);

  ClassDef(AliHLTTPCTrackGeometry, 0)
};
#endif
