#ifndef PMDdigit_H
#define PMDdigit_H
//-----------------------------------------------------//
//                                                     //
//  Date   : August 05 2003                            //
//                                                     //
//  Store digits  for PMD                              //
//                                                     //
//-----------------------------------------------------//

#include "Riostream.h"
#include "Rtypes.h"
#include "TObject.h"
#include "TClonesArray.h"

class AliPMDdigit : public TObject
{
  
 protected:
  Int_t   fTrNumber, fDet, fSMNumber, fCellNumber;
  Float_t fADC;

 public:
  AliPMDdigit();
  AliPMDdigit(Int_t /* trnumber */, Int_t /* det */, Int_t /* smnumber */,
	      Int_t /* cellnumber */, Float_t /* adc */);
  AliPMDdigit(AliPMDdigit *pmddigit) {*this = *pmddigit;}
  
  virtual ~AliPMDdigit();

  Int_t   GetTrackNumber() const;
  Int_t   GetDetector() const;
  Int_t   GetSMNumber() const;
  Int_t   GetCellNumber() const;
  Float_t GetADC() const;
  
  ClassDef(AliPMDdigit,1)
};

#endif
