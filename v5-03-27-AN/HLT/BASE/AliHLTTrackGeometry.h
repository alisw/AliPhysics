//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTTRACKGEOMETRY_H
#define ALIHLTTRACKGEOMETRY_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/// @file   AliHLTTrackGeometry.h
/// @author Matthias Richter
/// @date   2011-05-20
/// @brief  Desciption of a track by a sequence of track points
///

#include <vector>
#include <TObject.h>
#include "AliHLTLogging.h"
#include "AliHLTDataTypes.h"
#include "AliHLTStdIncludes.h"
#include "AliHLTExternalTrackParam.h"
#include "AliHLTIndexGrid.h"

class AliHLTGlobalBarrelTrack;
class AliHLTSpacePointContainer;
class TH2;

/**
 * @class AliHLTTrackGeometry
 * Description of a track by a sequence of track points given by the
 * intersection of the track with a set of planes.
 *
 * Each plane describes a local 2D coordinate system with origin at the
 * norm vector going through global {0,0,0}. The local coordinates u and v
 * correspond to global y and z at zero alpha and theta. In that case, r
 * corresponds to global x.
 *
 * Each plane is identified by
 * - rotation angle alpha in the global xy plane
 * - rotation angle theta
 * - unique 32bit AliHLTUInt32_t plane id
 *
 * The geometry of planes can be defined in child classes defining the
 * plane id.
 * @ingroup alihlt_base
 */
class AliHLTTrackGeometry : public TObject, public AliHLTLogging
{
 public:
  /// standard constructor
  AliHLTTrackGeometry();
  /// copy constructor
  AliHLTTrackGeometry(const AliHLTTrackGeometry&);
  /// assignment operator
  AliHLTTrackGeometry& operator=(const AliHLTTrackGeometry&);

  /// destructor
  ~AliHLTTrackGeometry();

  typedef AliHLTIndexGrid<int, AliHLTUInt32_t> AliHLTTrackGrid;

  /// set the track id
  void SetTrackId(int trackId) {fTrackId=trackId;}
  /// get the track id
  int GetTrackId() const {return fTrackId;}

  enum {
    kLower = 0,
    kUpper = 1,
    kBoundsDimension
  };

  /**
   * class AliHLTTrackPlane
   * Helper class to describe a plane
   */
  class AliHLTTrackPlane {
  public:
    AliHLTTrackPlane(AliHLTUInt32_t id, float alpha, float r, float theta) 
      : fId(id), fAlpha(alpha), fR(r), fTheta(theta), fBoundsU(), fBoundsV() {
      fBoundsU[kLower]=-5.0; fBoundsU[kUpper]=5.0; fBoundsV[kLower]=-5.0; fBoundsV[kUpper]=5.0;
    }
    virtual ~AliHLTTrackPlane() {}

    /// id of the plane
    AliHLTUInt32_t GetId() const {return fId;}
    /// alpha of the plane
    float GetPlaneAlpha() const {return fAlpha;}
    /// radial distance from global {0,0,0}
    float GetPlaneR() const {return fR;}
    /// theta of the plane
    float GetPlaneTheta() const {return fTheta;}

    // set bounds of u coordinate
    void SetBoundsU(const float bounds[kBoundsDimension]) {
      fBoundsU[kLower]=bounds[kLower]; fBoundsU[kUpper]=bounds[kUpper];
    }
    // set bounds of v coordinate
    void SetBoundsV(const float bounds[kBoundsDimension]) {
      fBoundsV[kLower]=bounds[kLower]; fBoundsV[kUpper]=bounds[kUpper];
    }

    bool CheckBounds(float u, float v) const {
      return (u>=fBoundsU[kLower]) && (u<=fBoundsU[kUpper]) 
	&& (v>=fBoundsV[kLower]) && (v<=fBoundsV[kUpper]);
    }

  private:
    // standard constructor prohibited
    AliHLTTrackPlane();

    AliHLTUInt32_t fId; // unique id
    float fAlpha;       // alpha of the plane
    float fR;           // radial distance from global {0,0,0}
    float fTheta;       // theta of the plane
    float fBoundsU[kBoundsDimension];  // bounds u coordinate
    float fBoundsV[kBoundsDimension];  // bounds v coordinate
  };

  class AliHLTTrackPlaneYZ : public AliHLTTrackPlane {
  public:
    AliHLTTrackPlaneYZ(AliHLTUInt32_t id, float alpha, float r) 
      : AliHLTTrackPlane(id, alpha, r, 0.0) {}
    ~AliHLTTrackPlaneYZ() {}

  private:
    // standard constructor prohibited
    AliHLTTrackPlaneYZ();
  };

  struct AliHLTTrackSpacepoint {
    AliHLTTrackSpacepoint(AliHLTUInt32_t id, float dU=-1000., float dV=-1000.)
      : fId(id), fdU(dU), fdV(dV) {}

    int SetResidual(int coordinate, float value) {
      if (coordinate==0) fdU=value;
      else if (coordinate==1) fdV=value;
      return 0;
    }

    float GetResidual(int coordinate) const {
      if (coordinate==0) return fdU;
      else if (coordinate==1) return fdV;
      return -1000.;
    }

    AliHLTUInt32_t fId; // space point id
    float fdU;          // residual of the spacepoint
    float fdV;          // residual of the spacepoint
  };

  class AliHLTTrackPoint {
  public:
    // constructor
    AliHLTTrackPoint(AliHLTUInt32_t id, float u, float v)
      : fId(id), fU(u), fV(v), fSpacePoints() {}
    // copy constructor
    AliHLTTrackPoint(const AliHLTTrackPoint& src)
      : fId(src.fId), fU(src.fU), fV(src.fV), fSpacePoints(src.fSpacePoints) {}
    // assignment operator
    AliHLTTrackPoint& operator=(const AliHLTTrackPoint& src) {
      if (this!=&src) {fId=src.fId; fU=src.fU; fV=src.fV; fSpacePoints=src.fSpacePoints;}
      return *this;
    }
    ~AliHLTTrackPoint() {}

    bool operator==(const AliHLTTrackPoint& point) const {
      return point.fId==fId;
    }
    bool operator==(AliHLTUInt32_t id) const {
      return id==fId;
    }

    /// id of the plane
    AliHLTUInt32_t GetId() const {return fId;}
    /// u coordinate
    float GetU() const {return fU;}
    /// u coordinate
    float GetV() const {return fV;}

    /// check associate space point
    bool HaveAssociatedSpacePoint() const {
      return fSpacePoints.size()>0;
    }

    /// associate a space point
    int AddAssociatedSpacePoint(AliHLTUInt32_t spacepointId, float dU=-1000., float dV=-1000.) {
      fSpacePoints.push_back(AliHLTTrackSpacepoint(spacepointId, dU, dV));
      return 0;
    }

    const vector<AliHLTTrackSpacepoint>& GetSpacepoints() const {return fSpacePoints;}
    vector<AliHLTTrackSpacepoint>& GetSpacepoints() {return fSpacePoints;}

    int SetResidual(AliHLTUInt32_t id, int coordinate, float value) {
      for (unsigned i=0; i<fSpacePoints.size(); i++) {
	if (fSpacePoints[i].fId!=id) continue;
	return fSpacePoints[i].SetResidual(coordinate, value);
      }
      return -ENOENT;
    }

  private:
    // standard constructor prohibited
    AliHLTTrackPoint();

    AliHLTUInt32_t fId; // unique id
    float fU;           // u coordinate
    float fV;           // v coordinate
    vector<AliHLTTrackSpacepoint> fSpacePoints;
  };

  // interface to plane description

  /// alpha of the plane
  virtual float GetPlaneAlpha(AliHLTUInt32_t planeId) const = 0;
  /// radial distance from global {0,0,0}
  virtual float GetPlaneR(AliHLTUInt32_t planeId) const = 0;
  /// theta of the plane
  virtual float GetPlaneTheta(AliHLTUInt32_t planeId) const = 0;
  /// check bounds in u and v coordinate
  virtual bool CheckBounds(AliHLTUInt32_t planeId, float u, float v) const =0;

  // track interface

  /// calculate the track points, expects the global magnetic field to be initialized
  virtual int CalculateTrackPoints(const AliHLTExternalTrackParam& track) = 0;

  /// calculate the track points, expects the global magnetic field to be initialized
  virtual int CalculateTrackPoints(AliHLTGlobalBarrelTrack& track) = 0;

  /// associate all space points of the container to the calculated track points
  int AssociateSpacePoints(AliHLTSpacePointContainer& points);
  /// associate specified space points of the container to the calculated track points
  int AssociateSpacePoints(const AliHLTUInt32_t* trackpoints, AliHLTUInt32_t nofPoints, AliHLTSpacePointContainer& points);
  int AssociateUnusedSpacePoints(AliHLTSpacePointContainer& points);

  /// find the track point which can be associated to a spacepoint with coordinates and id
  virtual int FindMatchingTrackPoint(AliHLTUInt32_t spacepointId, float spacepoint[3], AliHLTUInt32_t& planeId) = 0;

  /// register track points in the index grid
  virtual int RegisterTrackPoints(AliHLTTrackGrid* pGrid) const;

  /// fill track points to index grid
  virtual int FillTrackPoints(AliHLTTrackGrid* pGrid) const;

  /// get track point of id
  const AliHLTTrackPoint* GetTrackPoint(AliHLTUInt32_t id) const;
  /// get track point of id
  AliHLTTrackPoint* GetTrackPoint(AliHLTUInt32_t id);

  /// create a collection of all points
  virtual AliHLTSpacePointContainer* ConvertToSpacePoints() const {return ConvertToSpacePoints(false);}
  virtual AliHLTSpacePointContainer* ConvertToSpacePoints(bool bAssociated) const;

  /// set the spacepoint associated with a track point
  /// @param  planeId       track point
  /// @param  spacepointId  space point id to be associated with track point
  /// @param  status        status flag
  /// @return 0 on success, -ENOENT planeId not found
  int SetAssociatedSpacePoint(AliHLTUInt32_t planeId, AliHLTUInt32_t spacepointId, int status, float fdU=0.0, float fdV=0.0);

  /// get the spacepoint associated with a track point
  /// @param  planeId       id of the track point
  /// @param  spacepointId  target to get the spacepoint data
  /// @return status flag if found, -ENOENT planeId not found, -ENODATA no associated spacepoint found
  int GetAssociatedSpacePoint(AliHLTUInt32_t planeId, AliHLTUInt32_t& spacepointId) const;

  // services

  int FillResidual(int coordinate, TH2* histo) const;

  void SetVerbosity(int verbosity) {fVerbosity=verbosity;}
  int GetVerbosity() const {return fVerbosity;}

  /// inherited from TObject: clear the object and reset pointer references
  virtual void Clear(Option_t * /*option*/ ="");

  /// inherited from TObject
  virtual void Print(Option_t *option="") const;

  virtual void Print(ostream& out, Option_t *option="") const;

  /// Inherited from TObject, draw the track points
  virtual void Draw(Option_t *option="");

 protected:
  int AddTrackPoint(const AliHLTTrackPoint& point, AliHLTUInt32_t selectionMask=kAliHLTVoidDataSpec);

  const vector<AliHLTTrackPoint>& TrackPoints() const {return fTrackPoints;}

 private:
  vector<AliHLTTrackPoint> fTrackPoints; // list of points
  vector<AliHLTUInt32_t> fSelectionMasks; // selection masks

  int fTrackId; // track id
  int fVerbosity; // verbosity

  ClassDef(AliHLTTrackGeometry, 0)
};

ostream& operator<<(ostream &out, const AliHLTTrackGeometry& c);

#endif
