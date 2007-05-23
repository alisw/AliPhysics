#ifndef ALIITSONLINESDD_H
#define ALIITSONLINESDD_H


///////////////////////////////////////////////////////////////////
//                                                               //
// Base class for SDD detector algorithms                        //
// Origin: F.Prino, Torino, prino@to.infn.it                     //
//                                                               //
///////////////////////////////////////////////////////////////////

#include<TObject.h>
#include<TMath.h>

class AliITSOnlineSDD : public TObject {

 public:
  AliITSOnlineSDD();
  AliITSOnlineSDD(Int_t mod, Int_t sid);
  virtual ~AliITSOnlineSDD(){};

  void SetModule(Int_t mod){fModuleId=mod;}
  void SetDetectorSide(Int_t sid){fSide=sid;}

  Int_t GetModuleId() const {return fModuleId;}
  Int_t GetDetectorSide() const {return fSide;}

 protected:
  static const Int_t fgkNAnodes = 256;
  Int_t fModuleId; // module number from 0 to 255
  Int_t fSide;     // detector side (0-1)

  ClassDef(AliITSOnlineSDD,1);
};
#endif
