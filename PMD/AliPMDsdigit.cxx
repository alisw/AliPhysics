//-----------------------------------------------------//
//                                                     //
//  Date   : August 05 2003                            //
//  used to store the info into TreeS                  //
//                                                     //
//-----------------------------------------------------//

#include "AliPMDsdigit.h"
#include <stdio.h>

ClassImp(AliPMDsdigit)

AliPMDsdigit::AliPMDsdigit()
{
  fTrNumber   = 0;
  fDet        = 0;
  fSMN        = 0;
  fCellNumber = 0;
  fEdep       = 0.;
}

AliPMDsdigit::AliPMDsdigit(Int_t trnumber, Int_t det, Int_t smn,
			   Int_t cellnumber, Float_t edep)
{
  fTrNumber   = trnumber;
  fDet        = det;
  fSMN        = smn;
  fCellNumber = cellnumber;
  fEdep       = edep;
}
AliPMDsdigit::~AliPMDsdigit()
{

}
Int_t AliPMDsdigit::GetTrackNumber() const
{
  return fTrNumber;
}
Int_t AliPMDsdigit::GetDetector() const
{
  return fDet;
}
Int_t AliPMDsdigit::GetSMNumber() const
{
  return fSMN;
}
Int_t AliPMDsdigit::GetCellNumber() const
{
  return fCellNumber;
}

Float_t AliPMDsdigit::GetCellEdep() const
{
  return fEdep;
}
