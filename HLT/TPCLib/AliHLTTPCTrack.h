// XEmacs -*-C++-*-
// $Id$
// Original: AliHLTTrack.h,v 1.18 2005/03/31 04:48:58 cvetan 

#ifndef ALIHLTTPCTRACK_H
#define ALIHLTTPCTRACK_H

//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/// @file   AliHLTTPCTrack.h
/// @author Anders Vestbo, Uli Frankenfeld, maintained by Matthias Richter
/// @date   
/// @brief  HLT TPC track base class (conformal mapping)
///

#include "AliTPCtrack.h"

class AliHLTTPCVertex;
struct AliHLTTPCSpacePointData;

/**
 * @class AliHLTTPCTrack
 * This class implements the representation of a TPC track, used by the
 * HLT conformal mapping track finder. <br>
 * It was originally separated from the offline TPC track class, but in
 * order to adjust the output format to the offline ESD, AliHLTTPCTrack
 * now inherits from AliHLTtrack.
 *
 * @ingroup alihlt_tpc
 */
class AliHLTTPCTrack : public AliTPCtrack {

 public:
  
  AliHLTTPCTrack();
  virtual ~AliHLTTPCTrack();
  
  /**
   * Copy track parameters.
   * @param track   pointer to source track
   */
  virtual void Copy(AliHLTTPCTrack* track);
  using AliTPCtrack::Copy;  //TODO: Check if "virtual void Copy(TObject*)" does what it is supposed to do.

  /**
   * Compare two tracks by the number of hits
   * @return 0 if equal number of hits, 
   *         1 if this > track
   *        -1 if this < track
   */
  virtual Int_t Compare(const AliHLTTPCTrack *track) const;
  using AliTPCtrack::Compare;  //TODO: Check if "virtual Int_t Compare(TObject*)" does what it is supposed to do.

  /**
   * Fit the assigned spacepoints to a helix.
   * The function sets teh track parameters.
   */
  virtual void CalculateHelix();
  
  Bool_t CalculateReferencePoint(Double_t angle,Double_t radius=132);//Calculate Reference Point
  Bool_t CalculateEdgePoint(Double_t angle);//Calculate crossing point with line
  Bool_t CalculatePoint(Double_t xplane);   //Calculate crossing point with X-plane
  Bool_t IsPoint() {return fIsPoint;}
  Double_t GetCrossingAngle(Int_t padrow,Int_t slice=-1);
  Bool_t GetCrossingPoint(Int_t padrow,Float_t *xyz);
  Double_t GetDistance(Double_t /*x0*/,Double_t /*x1*/){return 0;}
  void UpdateToFirstPoint();
  void GetClosestPoint(AliHLTTPCVertex *vertex,Double_t &closestX,Double_t &closestY,Double_t &closestZ);
  void Rotate(Int_t slice,Bool_t tolocal=kFALSE);
  Bool_t IsLocal() const {return fIsLocal;}
  virtual void Print(Option_t* option = "") const;
  using AliTPCtrack::Print;

  // getter
  Double_t GetFirstPointX() const {return fFirstPoint[0];}
  Double_t GetFirstPointY() const {return fFirstPoint[1];}
  Double_t GetFirstPointZ() const {return fFirstPoint[2];}
  Double_t GetLastPointX() const {return fLastPoint[0];}
  Double_t GetLastPointY() const {return fLastPoint[1];}
  Double_t GetLastPointZ() const {return fLastPoint[2];}

  Double_t GetPointPsi() const {return fPointPsi;}
  Double_t GetPointX() const {return fPoint[0];}
  Double_t GetPointY() const {return fPoint[1];}
  Double_t GetPointZ() const {return fPoint[2];}

  Double_t GetPt() const {return fPt;}
  Double_t GetTgl() const {return fTanl;}
  Double_t GetPsi() const {return fPsi;}
  Double_t GetPhi0() const {return fPhi0;}
  Double_t GetR0() const {return fR0;}
  Double_t GetZ0() const {return fFirstPoint[2];}
  Float_t GetPID() const {return fPID;}

  Double_t GetPterr() const {return fPterr;}
  Double_t GetPsierr() const {return fPsierr;}
  Double_t GetTglerr() const {return fTanlerr;}
  Double_t GetZ0err() const {return fZ0err;}
  Double_t GetY0err() const {return fY0err;}

  Double_t GetKappa() const {return fKappa;}
  Double_t GetRadius() const {return fRadius;}
  Double_t GetCenterX() const {return fCenterX;}
  Double_t GetCenterY() const {return fCenterY;}

  Int_t GetNHits() const {return fNHits;}
  Int_t   GetNumberOfPoints()   const {return fNHits;}
  Bool_t  ComesFromMainVertex() const {return fFromMainVertex;}
    
  Double_t GetPx() const {return fPt*cos(fPsi);}
  Double_t GetPy() const {return fPt*sin(fPsi);}
  Double_t GetPz() const {return fPt*fTanl;}
  
  Double_t GetP() const;
  Double_t GetPseudoRapidity() const;
  Double_t GetRapidity() const;
  
  Int_t GetCharge() const {return fQ;}
  Int_t GetMCid() const {return fMCid;}
  Double_t GetLength() const {return fLength;}
  Double_t GetLengthXY() const ;
  Double_t GetLengthTot() const;
  
  Int_t GetFirstRow() const {return fRowRange[0];}
  Int_t GetLastRow()  const {return fRowRange[1];}
  Int_t GetSector()   const {return fSector;}

  UInt_t *GetHitNumbers() {return fHitNumbers;}
  Int_t GetId(){ return fId; }

  // setter   
  void SetPID(Float_t pid) {fPID=pid;}  
  void SetMCid(Int_t f) {fMCid = f;}
  void SetFirstPoint(Double_t f,Double_t g,Double_t h) {fFirstPoint[0]=f; fFirstPoint[1]=g; fFirstPoint[2]=h;}
  void SetLastPoint(Double_t f,Double_t g,Double_t h) {fLastPoint[0]=f; fLastPoint[1]=g; fLastPoint[2]=h;}
  void SetHits(Int_t nhits,UInt_t *hits);
  void SetPhi0(Double_t f) {fPhi0 = f;}
  void SetPsi(Double_t f) {fPsi = f;}
  void SetR0(Double_t f) {fR0 = f;}
  void SetTgl(Double_t f) {fTanl =f;}
  void SetZ0(Double_t f) {fFirstPoint[2] = f;}
  void SetPt(Double_t f) {fPt = f;}
  void SetLength(Double_t f) {fLength = f;}
  void SetPterr(Double_t f) {fPterr = f;}
  void SetPsierr(Double_t f) {fPsierr = f;}
  void SetZ0err(Double_t f) {fZ0err = f;}
  void SetY0err(Double_t f) {fY0err = f;}  
  void SetTglerr(Double_t f) {fTanlerr = f;}
  void SetKappa(Double_t f) {fKappa = f;}
  void SetNHits(Int_t f) {fNHits = f;}
  void SetRowRange(Int_t f,Int_t g) {fRowRange[0]=f; fRowRange[1]=g;}
  void SetSector(Int_t f) {fSector = f;}
  void SetRadius(Double_t f) {fRadius = f;}
  void SetCenterX(Double_t f) {fCenterX = f;}
  void SetCenterY(Double_t f) {fCenterY = f;}
  void SetCharge(Int_t f) {fQ = f;}
  void SetId( Int_t f ) { fId = f; }

  void ComesFromMainVertex(Bool_t f) {fFromMainVertex = f;}

  /**
   * Convert all track parameters to the format of AliKalmanTrack
   * The AliKalmanTrack class implements the track parametrization for offline ITS, TPC
   * and TRD tracking. The function calculates and sets the parameters of the
   * parent class (Note: AliHLTTPCTrack inherits from AliTPCtrack and thus
   * AliKalmanTrack).
   */
  int Convert2AliKalmanTrack();

  /**
   * Check the structure members to be within reasonable limits.
   */
  int CheckConsistency();

  /**
   * Check consistency of a double member
   */
  int CheckDoubleMember(double* pMember, double def, const char* name) const;

 private:

  Int_t fNHits; //Number of hits
  Int_t fMCid;  //Assigned id from MC data.

  Double_t fKappa;   // Signed curvature (projected to a circle)
  Double_t fRadius;  // Radius of the helix (projected to a circle)
  Double_t fCenterX; // x coordinate of the center of the helix (projected to a circle)
  Double_t fCenterY; // y coordinate of the center of the helix (projected to a circle)
  Bool_t   fFromMainVertex; // true if tracks origin is the main vertex, otherwise false
  
  Int_t fRowRange[2]; //Subsector where this track was build
  Int_t fSector;      //Sector # where  this track was build

  //data from momentum fit
  Int_t    fQ;    //charge measured fit
    
  //track parameters:
  Double_t fTanl; //tan of dipangle
  Double_t fPsi;  //azimuthal angle of the momentum 
  Double_t fPt;   //transverse momentum
  Double_t fLength; //length of track (s)
  
  Double_t fPterr;   //error in pt
  Double_t fPsierr;  //error in psi
  Double_t fZ0err;   //error in first point
  Double_t fY0err;   //error in first point
  Double_t fTanlerr; //error in tanl

  Double_t fPhi0; //azimuthal angle of the first point
  Double_t fR0;   //radius of the first point
  Double_t fZ0;   //z coordinate of the first point (fFirstPoint[2])

  Double_t fFirstPoint[3]; //first point
  Double_t fLastPoint[3];  //last point
  Double_t fPoint[3]; //point
  Double_t fPointPsi; //azimuthal angle of the momentum at Point

  Bool_t fIsPoint;  //Helix crosses the X-plane
  Bool_t fIsLocal; //Track given in local coordinates.

  Float_t fPID; //pid 
  static const int fgkHitArraySize=159; // size of hit array
  UInt_t fHitNumbers[fgkHitArraySize]; //Array of hit numbers for this track
  Int_t fId; // unique ID of the track

  Bool_t IsPoint(Bool_t ispoint) {fIsPoint = ispoint;return fIsPoint;}

  ClassDef(AliHLTTPCTrack,3) //Base track class
};
#endif
