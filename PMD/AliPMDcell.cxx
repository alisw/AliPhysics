//-----------------------------------------------------//
//                                                     //
//  Date   : August 05 2003                            //
//                                                     //
//  Store cell/track info which is used to assign      //
//  the correct track number to a multiple hit cell    //
//                                                     //
//-----------------------------------------------------//

#include "AliPMDcell.h"

ClassImp(AliPMDcell)

AliPMDcell::AliPMDcell()
{
  fTrNumber = 0;
  fSMNumber = 0;
  fXpos     = 0;
  fYpos     = 0;
  fEdep     = 0.;
}

AliPMDcell::AliPMDcell(Int_t trnumber, Int_t smnumber, 
			 Int_t xpos, Int_t ypos, Float_t edep)
{
  fTrNumber = trnumber;
  fSMNumber = smnumber;
  fXpos     = xpos;
  fYpos     = ypos;
  fEdep     = edep;

}
AliPMDcell::~AliPMDcell()
{

}

Int_t AliPMDcell::GetTrackNumber() const
{
  return fTrNumber;
}
Int_t AliPMDcell::GetSMNumber() const
{
  return fSMNumber;
}
Int_t AliPMDcell::GetX() const
{
  return fXpos;
}
Int_t AliPMDcell::GetY() const
{
  return fYpos;
}

Float_t AliPMDcell::GetEdep() const
{
  return fEdep;
}

