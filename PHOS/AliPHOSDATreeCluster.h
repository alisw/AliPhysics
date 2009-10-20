// --
// --
// Implementation for TTree output in PHOS DA
// for calibrating energy by pi0 and MIP.
// --
// -- Author: Hisayuki Torii (Hiroshima Univ.)
// --

#ifndef AliPHOSDATreeCluster_H
#define AliPHOSDATreeCluster_H

#include <iostream>
#include "AliPHOSDATreeDigit.h"

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
  virtual ~AliPHOSDATreeCluster(){ if(fNDigits>0) delete[] fDigits; };
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
    if( ndigit >= 0 && ndigit < fNDigits ) return fDigits[ndigit];
    else std::cout<<" AliPHOSDATreeCluster::GetDigit("<<ndigit<<")::Error. Out of range > "<<fNDigits<<std::endl;
  };
  AliPHOSDATreeDigit& GetMaxDigit(){
    if( fNDigits >= 0 ) return fDigits[0];
    else std::cout<<" AliPHOSDATreeCluster::GetMaxDigit()::Warning No digit information."<<std::endl;
  };
  bool Append(AliPHOSDATreeDigit& digit);
  bool Append(AliPHOSDATreeCluster& cluster);
  bool IsNeighbor(const AliPHOSDATreeDigit& digit) const;
  bool IsNeighbor(const AliPHOSDATreeCluster& cluster) const;
  void Reset();
  void Print(char* opt="");

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
