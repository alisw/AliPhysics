#ifndef AliPHOSDATreeCluster_H
#define AliPHOSDATreeCluster_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

// --
// --
// Implementation for TTree output in PHOS DA
// for calibrating energy by pi0 and MIP.
// --
// -- Author: Hisayuki Torii (Hiroshima Univ.)
// --

#include "AliPHOSDATreeDigit.h"
#include <iosfwd>

class AliPHOSDATreeCluster{

  friend std::ostream& operator<<(std::ostream& out,const AliPHOSDATreeCluster& cluster);

 public:

  //AliPHOSDATreeCluster():fEnergy(0),fX(0),fY(0),fZ(0),fNDigits(0),fDigits(0){/**/};
  //AliPHOSDATreeCluster(float energy,float x,float y,float z):fEnergy(energy),fX(x),fY(y),fZ(z),fNDigits(0),fDigits(0){/**/};
  //void Set(float energy,float x,float y,float z){fEnergy=energy; fX=x; fY=y; fZ=z;};
  //float GetX(){ return fX; };
  //float GetY(){ return fY; };
  //float GetZ(){ return fZ; };
  AliPHOSDATreeCluster():fEnergy(0),fRow(-100),fCol(-100),fNDigits(0),fDigits(0){/**/};
  AliPHOSDATreeCluster(float energy,float row,float col):fEnergy(energy),fRow(row),fCol(col),fNDigits(0),fDigits(0){/**/};
  virtual ~AliPHOSDATreeCluster(){ delete[] fDigits; };
  AliPHOSDATreeCluster(const AliPHOSDATreeCluster& cluster);
  AliPHOSDATreeCluster& operator=(const AliPHOSDATreeCluster& cluster);
  void Set(float energy,float row,float col){fEnergy=energy; fRow=row; fCol=col; };
  void SetEnergy(float energy){fEnergy=energy;};
  float GetEnergy() const{ return fEnergy; };
  float GetRow() const{ return fRow; };
  float GetCol() const{ return fCol; };
  bool CalculateProperty();
  int GetNDigits() const{ return fNDigits; };
  AliPHOSDATreeDigit& GetDigit(int ndigit){
    return fDigits[ndigit];
  };
  AliPHOSDATreeDigit& GetMaxDigit(){
    return fDigits[0];
  };
  bool Append(AliPHOSDATreeDigit& digit);
  bool Append(AliPHOSDATreeCluster& cluster);
  bool IsNeighbor(const AliPHOSDATreeDigit& digit) const;
  bool IsNeighbor(const AliPHOSDATreeCluster& cluster) const;
  void Reset();
  void Print(Option_t *option="") const;

 private:

  float fEnergy;               // Energy in GeV
  float fRow;                  // PHOS Internal Coordinates, 0 - 63
  float fCol;                  // PHOS Internal Coordinates, 0 - 55
  //float fX, fY, fZ;
  int fNDigits;                // Number of digits
  AliPHOSDATreeDigit* fDigits; //[fNDigits]

  ClassDef(AliPHOSDATreeCluster,1) // Simple Cluster Structure for PHOS DA
};

#endif
