#ifndef ALIFLATTPCCLUSTER_H
#define ALIFLATTPCCLUSTER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               *
 * Primary Authors : Sergey Gorbunov, Jochen Thaeder, Chiara Zampolli     */

/**
 * >> Flat structure representing a TPC cluster <<
 */

#include "Rtypes.h"
#include "AliVVcluster.h"
#include "AliFlatESDMisc.h"

class AliFlatTPCCluster: public AliVVcluster
{
 friend class AliFlatESDTrack;
  public:
  void SetX(Float_t x)             {fX = x;}
  void SetY(Float_t y)             {fY = y;}
  void SetZ(Float_t z)             {fZ = z;}
  void SetPadRow(Short_t padrow)   {fPadRow = padrow;}
  void SetSigmaY2(Float_t sigmaY2) {fSigmaY2 = sigmaY2;}
  void SetSigmaZ2(Float_t sigmaZ2) {fSigmaZ2 = sigmaZ2;}
  void SetCharge(UShort_t charge)  {fCharge = charge;}
  void SetQMax(UShort_t qmax)      {fQMax = qmax;}

  Float_t  GetX()       const      {return fX;}
  Float_t  GetY()       const      {return fY;}
  Float_t  GetZ()       const      {return fZ;}
  UShort_t GetPadRow()  const      {return fPadRow;}
  Float_t  GetSigmaY2() const      {return fSigmaY2;}
  Float_t  GetSigmaZ2() const      {return fSigmaZ2;}
  UShort_t GetCharge()  const      {return fCharge;}
  UShort_t GetQMax()    const      {return fQMax;}

  AliFlatTPCCluster() 
  : fX(0.), fY(0.), fZ(0.), fPadRow(0), fSigmaY2(0.), fSigmaZ2(0.), fCharge(0), fQMax(0) {}

  static Bool_t SortClusters(const AliFlatTPCCluster &first, const AliFlatTPCCluster &second){
    // Method to sort two clusters according to pad row
    Int_t padrowfirst  = first.GetPadRow();  
    Int_t padrowsecond = second.GetPadRow();
    return (padrowfirst < padrowsecond);
  }
  
  private:
  AliFlatTPCCluster(AliFlatESDSpecialConstructorFlag) {}
  virtual ~AliFlatTPCCluster() {}
  Float_t fX;       // X coordinate in local coordinates
  Float_t fY;       // Y coordinate in local coordinates
  Float_t fZ;       // Z coordinate in local coordinates
  UChar_t fPadRow;  // Pad row number
  Float_t fSigmaY2; // error (former width) of the clusters
  Float_t fSigmaZ2; // error (former width) of the clusters
  UInt_t  fCharge;  // total charge of cluster
  UInt_t  fQMax;    // QMax of cluster
  
  
  
};

#endif
