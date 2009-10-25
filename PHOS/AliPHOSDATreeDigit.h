#ifndef AliPHOSDATreeDigit_H
#define AliPHOSDATreeDigit_H
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

class AliPHOSDATreeDigit {

  friend std::ostream& operator<<(std::ostream& out,const AliPHOSDATreeDigit& digit);

 public:

  AliPHOSDATreeDigit():fEnergy(0),fAbsId(-1){/**/};
  AliPHOSDATreeDigit(float energy,int row,int col):fEnergy(energy),fAbsId(row*56+col){/**/};
  AliPHOSDATreeDigit(float energy,short absid):fEnergy(energy),fAbsId(absid){/**/};
  virtual ~AliPHOSDATreeDigit(){/**/};
  void Set(float energy,int row,int col){fEnergy=energy; fAbsId = row*56 + col;};
  void Set(float energy,short absid){fEnergy=energy; fAbsId = absid;};
  void SetEnergy(float energy){ fEnergy=energy; };
  float GetEnergy() const{return fEnergy;};
  int GetCol() const{return fAbsId % 56; };     // 0 - 55 : col( z ) number
  int GetRow() const{return (int) fAbsId/56; }; // 0 - 63 : row( x ) number
  short GetAbsId() const{return fAbsId; };
  short GetDigitId() const{return fAbsId; };
  bool IsValid() const { if( fAbsId>=0 && fAbsId < 3584 ) return true; return false; };
  bool IsNeighbor(const AliPHOSDATreeDigit &digit) const;
  void Print(Option_t *opt="") const;

  //AliPHOSDATreeDigit(float energy,int mod,int row,int col):fEnergy(energy),fAbsId(mod*3584+row*56+col){/**/};
  //void Set(float energy,int mod,int row,int col){
  //fEnergy=energy; fAbsId = mod * 3584 + row*56 + col;};
  //int GetMod(){return (int)(fAbsId/3584); };      // 0 -  4 : module number
  //int GetCol(){return (fAbsId%3584) % 56; };      // 0 - 55 : col( z ) number
  //int GetRow(){return (int)(fAbsId%3584) / 56; }; // 0 - 63 : row( x ) number
  //bool IsValid(){ if( fAbsId>=0 && fAbsId < 17920 ) return true; return false; };

 private:

  float fEnergy; // Energy in GeV
  short fAbsId;  // col + row * 56 [+mod*3584 in some cases]

  ClassDef(AliPHOSDATreeDigit,1) // Simple Digit Structure for PHOS DA
};

#endif
