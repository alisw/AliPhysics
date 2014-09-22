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

class AliFlatTPCCluster
{
  public:

  AliFlatTPCCluster() : fX(0.), fY(0.), fZ(0.), fSector(0), fPadRow(0), fSigmaY2(0.), fSigmaZ2(0.), fCharge(0), fQMax(0) {}

  AliFlatTPCCluster(AliVConstructorReinitialisationFlag ); // do nothing
 
  void Reinitialize(){} // do nothing

  ~AliFlatTPCCluster() {}

  void SetX(Float_t x)             {fX = x;}
  void SetY(Float_t y)             {fY = y;}
  void SetZ(Float_t z)             {fZ = z;}
  void SetSector(UShort_t sector)  {fSector = sector;}
  void SetPadRow(UShort_t padrow)  {fPadRow = padrow;}
  void SetSigmaY2(Float_t sigmaY2) {fSigmaY2 = sigmaY2;}
  void SetSigmaZ2(Float_t sigmaZ2) {fSigmaZ2 = sigmaZ2;}
  void SetCharge(UShort_t charge)  {fCharge = charge;}
  void SetQMax(UShort_t qmax)      {fQMax = qmax;}
  
 
  Float_t  GetX()       const      {return fX;}
  Float_t  GetY()       const      {return fY;}
  Float_t  GetZ()       const      {return fZ;}
  UShort_t GetSector()  const      {return fSector;}
  UShort_t GetPadRow()  const      {return fPadRow;}
  Float_t  GetSigmaY2() const      {return fSigmaY2;}
  Float_t  GetSigmaZ2() const      {return fSigmaZ2;}
  UShort_t GetCharge()  const      {return fCharge;}
  UShort_t GetQMax()    const      {return fQMax;}

  void SetTPCCluster( const AliTPCclusterMI *c );
  void GetTPCCluster( AliTPCclusterMI *c ) const;

  private:

  Float_t fX;       // X coordinate in local coordinates
  Float_t fY;       // Y coordinate in local coordinates
  Float_t fZ;       // Z coordinate in local coordinates
  UChar_t fSector;  // TPC sector
  UChar_t fPadRow;  // Pad row number withing the sector
  Float_t fSigmaY2; // error (former width) of the clusters
  Float_t fSigmaZ2; // error (former width) of the clusters
  UInt_t  fCharge;  // total charge of cluster
  UInt_t  fQMax;    // QMax of cluster
  
};

#pragma GCC diagnostic ignored "-Weffc++" 
inline AliFlatTPCCluster::AliFlatTPCCluster(AliVConstructorReinitialisationFlag ){}
#pragma GCC diagnostic warning "-Weffc++" 

inline void AliFlatTPCCluster::SetTPCCluster( const AliTPCclusterMI *c )
{
  SetX( c->GetX() );
  SetY( c->GetY() );
  SetZ( c->GetZ() );
  SetSector( c->GetDetector() );
  SetPadRow( c->GetRow() );
  SetSigmaY2( c->GetSigmaY2() );
  SetSigmaZ2( c->GetSigmaZ2() );
  SetCharge( c->GetQ() );
  SetQMax( c->GetMax() );
}
 
inline void AliFlatTPCCluster::GetTPCCluster( AliTPCclusterMI *c ) const
{
  if( !c ) return;
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

#endif
