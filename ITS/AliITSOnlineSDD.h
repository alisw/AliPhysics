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
  AliITSOnlineSDD(Int_t mod, Int_t sid);
  virtual ~AliITSOnlineSDD(){};

  void SetModule(Int_t mod){fModuleId=mod;}
  void SetDetectorSide(Int_t sid){fSide=sid;}
  void SetFirstGoodTB(Int_t itb=1){fFirstGoodTB=itb;}
  void SetLastGoodTB(Int_t itb=254){fLastGoodTB=itb;}

  Int_t GetModuleId() const {return fModuleId;}
  Int_t GetDetectorSide() const {return fSide;}
  Int_t GetFirstGoodTB() const {return fFirstGoodTB;}
  Int_t GetLastGoodTB() const {return fLastGoodTB;}

 protected:
  static const Int_t fgkNAnodes = 256; // number of anodes in each half-module
  Int_t fModuleId;    // module number from 0 to 255
  Int_t fSide;        // detector side (0-1)
  Int_t fFirstGoodTB; // first good time bin (to exclude time bin 0)
  Int_t fLastGoodTB;  // last good time bin (to exclude time bin 255)

  ClassDef(AliITSOnlineSDD,2);
};
#endif
