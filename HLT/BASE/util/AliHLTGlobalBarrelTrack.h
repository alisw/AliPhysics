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

class TClonesArray;
class AliHLTSpacePointContainer;
class AliHLTTrackGeometry;

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
  /** copy constructor */
  AliHLTGlobalBarrelTrack(const AliExternalTrackParam& p);

  /// assignment operator
  /// the standard assignment operator for AliHLTGlobalBarrelTrack is in principle
  /// covered by the template definition, however, compiler does not seem to recognize
  /// correctly -> effC++ warning, 
  AliHLTGlobalBarrelTrack& operator=(const AliHLTGlobalBarrelTrack& t) {
    if (this==&t) return *this;
    this->~AliHLTGlobalBarrelTrack(); new (this) AliHLTGlobalBarrelTrack(t);
    return *this;
  }
  template <class c>
  AliHLTGlobalBarrelTrack& operator=(const c& t) {
    this->~AliHLTGlobalBarrelTrack(); new (this) AliHLTGlobalBarrelTrack(t);
    return *this;
  }

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

  /// Set the space point data
  void SetSpacePointContainer(AliHLTSpacePointContainer* points) {fSpacePoints=points;}
  /// Set track point container
  void SetTrackGeometry(AliHLTTrackGeometry* points);
  /// Get track point container
  AliHLTTrackGeometry* GetTrackGeometry() const {return fTrackPoints;}

  /// associate the track space points to the calculated track points
  int AssociateSpacePoints(AliHLTTrackGeometry* trackpoints, AliHLTSpacePointContainer& spacepoints) const;

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

  /// Inherited from TObject, draw the track
  virtual void Draw(Option_t *option="");

  int DrawProjXYSpacePoints(Option_t *option, const AliHLTSpacePointContainer* fSpacePoints, const float scale, float center[2]);
  int DrawProjXYTrack(Option_t *option, const float scale, float center[2]);
  int DrawProjXYTrackPoints(Option_t *option, const float scale, const float center[2], int firstpadrow, int step, float lastpoint[2]);

  Double_t GetPathLengthTo( Double_t x, Double_t b ) const;

  static int ReadTracks(const char* filename, TClonesArray& tgt, AliHLTComponentDataType dt=kAliHLTVoidDataType, unsigned specification=kAliHLTVoidDataSpec);
  static int ReadTrackList(const char* listfile, TClonesArray& tgt, AliHLTComponentDataType dt=kAliHLTVoidDataType, unsigned specification=kAliHLTVoidDataSpec);

  /// calculate crossing point with a plane parallel to z axis, at distance x
  /// and phi
  int CalculateCrossingPoint(float xPlane, float phiPlane, float& u, float& v);

  /// calculate and set internal helix parameters using the global magnetic field
  int CalculateHelixParams();
  /// calculate and set internal helix parameters
  int CalculateHelixParams(float bfield);

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

  /// helix parameters
  Float_t fHelixRadius; //
  Float_t fHelixCenterX; //
  Float_t fHelixCenterY; //

  /// the space points assigned to the track
  AliHLTSpacePointContainer* fSpacePoints; //!
  /// the track points according to a geometry
  AliHLTTrackGeometry* fTrackPoints; //!

  ClassDef(AliHLTGlobalBarrelTrack, 0)
};
#endif
