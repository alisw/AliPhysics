//-----------------------------------------------------//
//                                                     //
//  Date   : August 05 2003                            //
//                                                     //
//  Store digits for ALICE-PMD                         //
//                                                     //
//-----------------------------------------------------//

#include "AliPMDdigit.h"
#include <stdio.h>

ClassImp(AliPMDdigit)

AliPMDdigit::AliPMDdigit()
{
  fTrNumber   = 0;
  fDet        = 0;
  fSMNumber   = 0;
  fCellNumber = 0;
  fADC        = 0.;
}

AliPMDdigit::AliPMDdigit(Int_t trnumber, Int_t det, Int_t smnumber, 
			 Int_t cellnumber, Float_t adc)
{
  fTrNumber   = trnumber;
  fDet        = det;
  fSMNumber   = smnumber;
  fCellNumber = cellnumber;
  fADC        = adc;
}
AliPMDdigit::~AliPMDdigit()
{

}
Int_t AliPMDdigit::GetTrackNumber() const
{
  return fTrNumber;
}
Int_t AliPMDdigit::GetDetector() const
{
  return fDet;
}
Int_t AliPMDdigit::GetSMNumber() const
{
  return fSMNumber;
}
Int_t AliPMDdigit::GetCellNumber() const
{
  return fCellNumber;
}
Float_t AliPMDdigit::GetADC() const
{
  return fADC;
}

