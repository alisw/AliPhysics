/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD calibration class for TRD DCS GTU parameters                         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliTRDCalDCSGTU.h"

ClassImp(AliTRDCalDCSGTU)

//_____________________________________________________________________________
AliTRDCalDCSGTU::AliTRDCalDCSGTU()
  :TNamed()
  ,fDCSID(0)
  ,fSMMask()
  ,fStackMask()
  ,fLinkMask()
{
  //
  // AliTRDCalDCSGTU default constructor
  //

  // link not active: 0, link active: 1, bit not set: 2
  for (Int_t i=0;i<18;i++) {
    fSMMask[i] = 2;
    for (Int_t j=0;j<5;j++) {
      fStackMask[i][j] = 2;
      for (Int_t k=0;k<12;k++) {
	fLinkMask[i][j][k] = 2;
      }
    }
  }

}

//_____________________________________________________________________________
AliTRDCalDCSGTU::AliTRDCalDCSGTU(const char *name, const char *title)
  :TNamed(name,title)
  ,fDCSID(0)
  ,fSMMask()
  ,fStackMask()
  ,fLinkMask()
{
  //
  // AliTRDCalDCSGTU constructor
  //

  // link not active: 0, link active: 1, bit not set: 2
  for(int i=0;i<18;i++) {
    fSMMask[i] = 2;
    for(int j=0;j<5;j++) {
      fStackMask[i][j] = 2;
      for(int k=0;k<12;k++) {
	fLinkMask[i][j][k] = 2;
      }
    }
  }

}

//_____________________________________________________________________________
Int_t AliTRDCalDCSGTU::SetSMMask(const char* smmask)
{
  // if something goes wrong here, the errorcode is 10x
  TString smMaskStr = smmask;
  // return false in case of wrong length
  if (smMaskStr.Length() != 18) return 101;

  for (Int_t i=0;i<18;i++) {
    TString bit = smMaskStr[i];
    if (!bit.IsDigit()) return 102; // must be 0 or 1 -> a digit!
    fSMMask[i] = bit.Atoi();
  }
  return 0;
}

//_____________________________________________________________________________
Int_t AliTRDCalDCSGTU::SetLinkMask(Int_t smid, Int_t stid, const char* lkmask)
{
  // if something goes wrong here, the errorcode is 11x
  if (smid == 99 || stid == 99) return 111; // 99: missing assignment

  TString lkMaskStr = lkmask;
  // return false in case of wrong length
  if (lkMaskStr.Length() != 12) return 112;

  for (Int_t i=0;i<12;i++) {
    TString bit = lkMaskStr[i];
    if (!bit.IsDigit()) return 113; // must be 0 or 1
    fLinkMask[smid][stid][i] = bit.Atoi();
  }
  return 0;
}

