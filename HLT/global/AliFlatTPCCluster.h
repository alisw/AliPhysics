#ifndef ALIFLATTPCCLUSTER_H
#define ALIFLATTPCCLUSTER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               *
 * Primary Authors : Sergey Gorbunov, Jochen Thaeder, Chiara Zampolli     */

/**
 * >> Flat structure representing a TPC cluster <<
 */

#include "Rtypes.h"
#include "AliVMisc.h"
#include "AliTPCclusterMI.h"
#include "AliComplexCluster.h"

class AliFlatTPCCluster
{
  public:

 AliFlatTPCCluster() : fPad(0.), fTimeBin(0.), fX(0.), fY(0.), fZ(0.), fSector(0), fPadRow(0), fSigmaY2(0.), fSigmaZ2(0.), fCharge(0), fQMax(0), fTrackAngleY(0), fTrackAngleZ(0) {}

  AliFlatTPCCluster(AliVConstructorReinitialisationFlag ); // do nothing
 
  void Reinitialize(){} // do nothing

  ~AliFlatTPCCluster() {}

  void SetPad(Float_t pad)         {fPad = pad;}
  void SetTimeBin(Float_t time)    {fTimeBin = time;}
  void SetX(Float_t x)             {fX = x;}
  void SetY(Float_t y)             {fY = y;}
  void SetZ(Float_t z)             {fZ = z;}
  void SetSector(UShort_t sector)  {fSector = sector;}
  void SetPadRow(UShort_t padrow)  {fPadRow = padrow;}
  void SetSigmaY2(Float_t sigmaY2) {fSigmaY2 = sigmaY2;}
  void SetSigmaZ2(Float_t sigmaZ2) {fSigmaZ2 = sigmaZ2;}
  void SetCharge(UShort_t charge)  {fCharge = charge;}
  void SetQMax(UShort_t qmax)      {fQMax = qmax;}
  void SetTrackAngleY( Float_t angY ) {fTrackAngleY = angY;}  
  void SetTrackAngleZ( Float_t angZ ) {fTrackAngleZ = angZ;}
 
  Float_t  GetPad()     const      {return fPad;}
  Float_t  GetTimeBin() const      {return fTimeBin;}
  Float_t  GetX()       const      {return fX;}
  Float_t  GetY()       const      {return fY;}
  Float_t  GetZ()       const      {return fZ;}
  UShort_t GetSector()  const      {return fSector;}
  UShort_t GetPadRow()  const      {return fPadRow;}
  Float_t  GetSigmaY2() const      {return fSigmaY2;}
  Float_t  GetSigmaZ2() const      {return fSigmaZ2;}
  UShort_t GetCharge()  const      {return fCharge;}
  UShort_t GetQMax()    const      {return fQMax;}
  
  Float_t GetTrackAngleY() const {return fTrackAngleY;}
  Float_t GetTrackAngleZ() const {return fTrackAngleZ;}

  void SetTPCCluster( const AliTPCclusterMI *c, const AliTPCTrackerPoints::Point *p );
  void GetTPCCluster( AliTPCclusterMI *c, AliTPCTrackerPoints::Point *p  ) const;

  private:
  
  Float_t fPad;     // Pad coordinate in local coordinates  
  Float_t fTimeBin; // Time coordinate in local coordinates  
  Float_t fX;       // X coordinate in local coordinates
  Float_t fY;       // Y coordinate in local coordinates
  Float_t fZ;       // Z coordinate in local coordinates
  UChar_t fSector;  // TPC sector
  UChar_t fPadRow;  // Pad row number withing the sector
  Float_t fSigmaY2; // error (former width) of the clusters
  Float_t fSigmaZ2; // error (former width) of the clusters
  UInt_t  fCharge;  // total charge of cluster
  UInt_t  fQMax;    // QMax of cluster
  Float_t fTrackAngleY; // tracker point angle Y
  Float_t fTrackAngleZ; // tracker point angle Z
};

#pragma GCC diagnostic ignored "-Weffc++" 
inline AliFlatTPCCluster::AliFlatTPCCluster(AliVConstructorReinitialisationFlag ){}
#pragma GCC diagnostic warning "-Weffc++" 

inline void AliFlatTPCCluster::SetTPCCluster( const AliTPCclusterMI *c, const AliTPCTrackerPoints::Point *p  )
{
  if( !c ) return;
  SetPad( c->GetPad() );
  SetTimeBin( c->GetTimeBin() );
  SetX( c->GetX() );
  SetY( c->GetY() );
  SetZ( c->GetZ() );
  SetSector( c->GetDetector() );
  SetPadRow( c->GetRow() );
  SetSigmaY2( c->GetSigmaY2() );
  SetSigmaZ2( c->GetSigmaZ2() );
  SetCharge( c->GetQ() );
  SetQMax( c->GetMax() );
  if( p ){
    SetTrackAngleY( p->GetAngleY() );  
    SetTrackAngleZ( p->GetAngleZ() );  
  } else {
    SetTrackAngleY( 0. );  
    SetTrackAngleZ( 0. );  
  }
}

inline void AliFlatTPCCluster::GetTPCCluster( AliTPCclusterMI *c, AliTPCTrackerPoints::Point *p ) const
{
  if( c ){
    c->SetPad( GetPad() );
    c->SetTimeBin( GetTimeBin() );
    c->SetX( GetX() );
    c->SetY( GetY() );
    c->SetZ( GetZ() );
    c->SetDetector( GetSector() );
    c->SetRow( GetPadRow() );
    c->SetSigmaY2( GetSigmaY2() );
    c->SetSigmaZ2( GetSigmaZ2() );
    c->SetQ( GetCharge() );
    c->SetMax( GetQMax() );
  }
  if( p ){
    p->SetAngleY( GetTrackAngleY() );
    p->SetAngleZ( GetTrackAngleZ() );
  }
}

#endif
