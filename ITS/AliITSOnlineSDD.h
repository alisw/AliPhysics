#ifndef ALIITSONLINESDD_H
#define ALIITSONLINESDD_H


///////////////////////////////////////////////////////////////////
//                                                               //
// Base class for SDD detector algorithms                        //
// Origin: F.Prino, Torino, prino@to.infn.it                     //
//                                                               //
///////////////////////////////////////////////////////////////////

#include<TObject.h>

class AliITSOnlineSDD : public TObject {

 public:
  AliITSOnlineSDD();
  AliITSOnlineSDD(Int_t nddl, Int_t ncarlos, Int_t sid);
  virtual ~AliITSOnlineSDD(){};

  void SetDDL(Int_t nd){fDDL=nd;}
  void SetCarlos(Int_t nc){fCarlos=nc;}
  void SetDetectorSide(Int_t sid){fSide=sid;}
  void SetFirstGoodTB(Int_t itb=1){fFirstGoodTB=itb;}
  void SetLastGoodTB(Int_t itb=127){fLastGoodTB=itb;}

  Int_t GetDDL() const {return fDDL;}
  Int_t GetCarlos() const {return fCarlos;}
  Int_t GetDetectorSide() const {return fSide;}
  Int_t GetFirstGoodTB() const {return fFirstGoodTB;}
  Int_t GetLastGoodTB() const {return fLastGoodTB;}

 protected:
  static const Int_t fgkNAnodes = 256; // number of anodes in each half-module
  Int_t fDDL;         // SDD DDL number (from 0 to 24)
  Int_t fCarlos;      // carlos number inside DDL (from 0 to 11)
  Int_t fSide;        // detector side (0-1)
  Int_t fFirstGoodTB; // first good time bin (to exclude time bin 0)
  Int_t fLastGoodTB;  // last good time bin (to exclude time bin 255)

  ClassDef(AliITSOnlineSDD,3);
};
#endif
