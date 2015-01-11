#ifndef ALIITS_FONOISESPD_H
#define ALIITS_FONOISESPD_H

/////////////////////////////////////////////////////////////////////
// Author: Henrik Tydesjo                                          //
//                                                                 //
// This class is used to store Fast-OR noise values in OCDB.       //
// One value per pixel chip.                                       //
// The values are the probability that a pixel chip will generate  //
// a fast-OR signal independently (originating from noise).        //
//                                                                 //
/////////////////////////////////////////////////////////////////////

#include <TObject.h>
#include <TError.h>

class AliITSFONoiseSPD : public TObject {

 public:
  AliITSFONoiseSPD();
  AliITSFONoiseSPD(const AliITSFONoiseSPD& foNoi);
  virtual ~AliITSFONoiseSPD();
  AliITSFONoiseSPD& operator=(const AliITSFONoiseSPD& foNoi);

  virtual void    ResetValues();
  virtual void    SetChipNoise(UInt_t eq, UInt_t hs, UInt_t chip, Float_t value);
  virtual Float_t GetChipNoise(UInt_t eq, UInt_t hs, UInt_t chip) const;
  virtual Float_t GetNoise(UInt_t eq, UInt_t hs, UInt_t chip, UInt_t /*col*/, UInt_t /*row*/) const
    {return GetChipNoise(eq,hs,chip);}

 protected:
  Float_t fChipNoise[20][6][10]; // noise values per chip

  ClassDef(AliITSFONoiseSPD,1)      // FO Noise Base
};

#endif
