#ifndef PMDcell_H
#define PMDcell_H
//-----------------------------------------------------//
//                                                     //
//  Date   : August 05 2003                            //
//                                                     //
//  Store cell/track info which is used to assign      //
//  the correct track number to a multiple hit cell    //
//                                                     //
//-----------------------------------------------------//

#include "Riostream.h"
#include "Rtypes.h"
#include "TObject.h"
#include "TClonesArray.h"

class AliPMDcell : public TObject
{
  
 protected:
  Int_t   fTrNumber, fSMNumber, fXpos, fYpos;
  Float_t fEdep;

 public:
  AliPMDcell();
  AliPMDcell(Int_t /* trnumber */, Int_t /* smnumber */,
	      Int_t /* xpos */, Int_t /* ypos */, Float_t /* edep */);
  AliPMDcell(AliPMDcell *pmdcell) {*this = *pmdcell;}
  
  virtual ~AliPMDcell();

  Int_t   GetTrackNumber() const;
  Int_t   GetSMNumber() const;
  Int_t   GetX() const;
  Int_t   GetY() const;
  Float_t GetEdep() const;
  
  ClassDef(AliPMDcell,1)
};

#endif
