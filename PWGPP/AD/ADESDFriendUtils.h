// -*- C++ -*-
// $Id$

/*  See cxx source for full Copyright notice */

//-----------------------------------------------------------------
//            AD ESD friend utilities
//-----------------------------------------------------------------

#include <TObject.h>

class AliADCalibData;
class AliESDADfriend;

class ADESDFriendUtils : public TObject {
public:
  ADESDFriendUtils();
  virtual ~ADESDFriendUtils();

  void    Init(Int_t runNumber);
  void    Update(const AliESDADfriend*);

  Float_t GetADCPedSub(Int_t ch, Int_t bc) const { return fADCPedSub[ch][bc]; }
  Bool_t  IsPileUp(Int_t ch, Float_t thr=20.0f) const;

protected:

private:
  // not implemented:
  ADESDFriendUtils(const ADESDFriendUtils&);
  ADESDFriendUtils& operator=(const ADESDFriendUtils&);

  AliADCalibData *fCalibData;         //! pedestals from OCDB
  Float_t         fADCPedSub[16][21]; //! pedestal subtracted ADC value

  ClassDef(ADESDFriendUtils, 1);
} ;
