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
#include <vector>

class AliHLTGlobalBarrelTrack;
class AliHLTDataDeflater;

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

  /// register track points in the index grid
  virtual int RegisterTrackPoints(AliHLTTrackGrid* pGrid) const;

  /// fill track points to index grid
  virtual int FillTrackPoints(AliHLTTrackGrid* pGrid) const;

  virtual int Write(const AliHLTGlobalBarrelTrack& track,
	    AliHLTSpacePointContainer* pSpacePoints,
	    AliHLTDataDeflater* pDeflater,
	    AliHLTUInt8_t* outputPtr,
	    AliHLTUInt32_t size,
	    vector<AliHLTUInt32_t>* writtenClusterIds=NULL,
	    const char* option="") const;

  virtual int WriteAssociatedClusters(AliHLTSpacePointContainer* pSpacePoints,
				      AliHLTDataDeflater* pDeflater,
				      vector<AliHLTUInt32_t>* writtenClusterIds=NULL,
				      const char* option="") const;

  int InitDriftTimeTransformation(float mA, float nA, float mC, float nC) {
    fDriftTimeFactorA=mA; fDriftTimeOffsetA=nA; fDriftTimeFactorC=mC; fDriftTimeOffsetC=nC; return 0;
  }

  struct AliHLTTPCTrackBlock {
    AliHLTUInt16_t   fSize; //! size in byte of the complete track block
    AliHLTUInt8_t    fSlice; //! slice number -> rotation angle of local coordinates
    AliHLTUInt8_t    fReserved; //! reserved field to fill 32bit
    AliHLTFloat32_t  fX; //! first X
    AliHLTFloat32_t  fY; //! first Y
    AliHLTFloat32_t  fZ; //! first Z
    AliHLTFloat32_t  fSinPsi; // local sine of the track momentum azimuthal angle
    AliHLTFloat32_t  fTgl;    // tangent of the track momentum dip angle
    AliHLTFloat32_t  fq1Pt;   // 1/pt (1/(GeV/c))
  };

  /// create a collection of all points
  virtual AliHLTSpacePointContainer* ConvertToSpacePoints() const {return ConvertToSpacePoints(false);}
  virtual AliHLTSpacePointContainer* ConvertToSpacePoints(bool bAssociated) const;

  /// get raw track point of id
  const AliHLTTrackPoint* GetRawTrackPoint(AliHLTUInt32_t id) const;
  /// get raw track point of id
  AliHLTTrackPoint* GetRawTrackPoint(AliHLTUInt32_t id);

  int FillRawResidual(int coordinate, TH2* histo, AliHLTSpacePointContainer* points) const;

  const vector<AliHLTTrackGeometry::AliHLTTrackPoint>& GetRawPoints() const {return fRawTrackPoints;}

 private:
  /// calculate the track points, expects the global magnetic field to be initialized
  int CalculateTrackPoints(AliHLTGlobalBarrelTrack& track, int firstpadrow, int step);

  vector<AliHLTTrackPoint> fRawTrackPoints; // list of points in raw coordinates

  float fDriftTimeFactorA; //! drift time A side
  float fDriftTimeOffsetA; //! drift time A side
  float fDriftTimeFactorC; //! drift time C side
  float fDriftTimeOffsetC; //! drift time C side

  ClassDef(AliHLTTPCTrackGeometry, 0)
};
#endif
