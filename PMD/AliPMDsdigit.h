#ifndef PMDsdigit_H
#define PMDsdigit_H
//-----------------------------------------------------//
//                                                     //
//                                                     //
//  Date   : August 05 2003                            //
//  used to store the info into TreeS                  //
//                                                     //
//-----------------------------------------------------//

#include "Riostream.h"
#include "Rtypes.h"
#include "TObject.h"
#include "TClonesArray.h"

class AliPMDsdigit : public TObject
{
  
 protected:

  Int_t   fTrNumber, fDet, fSMN, fCellNumber;
  Float_t fEdep;

 public:
  AliPMDsdigit();
  AliPMDsdigit(Int_t /* trnumber */, Int_t /* det */, Int_t /* smn */,
	       Int_t /* cellnumber */, Float_t /* edep */);
  AliPMDsdigit(AliPMDsdigit *pmdsdigit) {*this = *pmdsdigit;}
  
  virtual ~AliPMDsdigit();

  Int_t   GetTrackNumber() const;
  Int_t   GetDetector() const;
  Int_t   GetSMNumber() const;
  Int_t   GetCellNumber() const;
  Float_t GetCellEdep() const;
  
  ClassDef(AliPMDsdigit,1)
};

#endif
