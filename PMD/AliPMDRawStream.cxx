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

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
///
/// This class provides access to PMD digits in raw data.
///
/// It loops over all PMD digits in the raw data given by the AliRawReader.
/// The Next method goes to the next digit. If there are no digits left
/// it returns kFALSE.
/// Several getters provide information about the current digit.
///
///////////////////////////////////////////////////////////////////////////////

#include "AliPMDRawStream.h"
#include "AliRawReader.h"

ClassImp(AliPMDRawStream)


//_____________________________________________________________________________
AliPMDRawStream::AliPMDRawStream(AliRawReader* rawReader) :
  fRawReader(rawReader),
  fModule(-1),
  fPrevModule(-1),
  fMCM(-1),
  fChannel(-1),
  fRow(-1),
  fColumn(-1),
  fSignal(-1),
  fDetector(-1),
  fSMN(-1)
{
// create an object to read PMD raw digits

  fRawReader->Select(12);
}

//_____________________________________________________________________________
AliPMDRawStream::AliPMDRawStream(const AliPMDRawStream& stream) :
  TObject(stream),
  fRawReader(NULL),
  fModule(-1),
  fPrevModule(-1),
  fMCM(-1),
  fChannel(-1),
  fRow(-1),
  fColumn(-1),
  fSignal(-1),
  fDetector(-1),
  fSMN(-1)
{
// copy constructor

  Fatal("AliPMDRawStream", "copy constructor not implemented");
}

//_____________________________________________________________________________
AliPMDRawStream& AliPMDRawStream::operator = (const AliPMDRawStream& 
					      /* stream */)
{
// assignment operator

  Fatal("operator =", "assignment operator not implemented");
  return *this;
}

//_____________________________________________________________________________
AliPMDRawStream::~AliPMDRawStream()
{
// destructor

}


//_____________________________________________________________________________
Bool_t AliPMDRawStream::Next()
{
// read the next raw digit
// returns kFALSE if there is no digit left

  fPrevModule = fModule;

  UInt_t data;
  if (!fRawReader->ReadNextInt(data)) return kFALSE;

  fSignal  = data & 0x0FFF;
  fChannel = (data >> 12) & 0x003F;
  fMCM     = (data >> 18) & 0x07FF;

  Int_t  iddl  = fRawReader->GetDDLID();
  Int_t  ium;
  GetRowCol(iddl, fMCM, fChannel, ium, fRow, fColumn);

  if (iddl < 4)
    {
      fModule = iddl*6 + ium;
      fDetector = 0;
      fSMN = iddl*6 + ium;
    }
  else if (iddl == 4)
    {
      if (ium < 6)
	{
	  fModule = 24 + ium;
	  fSMN    = ium;
	}
      else if (ium >= 6)
	{
	  fModule = 30 + ium;
	  fSMN    = 6 + ium;
	}
      fDetector = 1;
    }
  else if (iddl == 5)
    {

      if (ium < 6)
	{
	  fModule = 30 + ium;
	  fSMN    = 6 + ium;
	}
      else if (ium >= 6)
	{
	  fModule = 36 + ium;
	  fSMN    = 12 + ium;
	}
      fDetector = 1;
    }


  return kTRUE;
}


//_____________________________________________________________________________
void AliPMDRawStream::GetRowCol(Int_t ddlno, UInt_t mcmno, UInt_t chno,
				Int_t &um, Int_t &row, Int_t &col) const
{
// decode: ddlno, mcmno, chno -> um, row, col


  Int_t remmcm = 0;
  Int_t divmcm = 0;

  static const UInt_t kCh[64] = { 21, 25, 29, 28, 17, 24, 20, 16,
				  12, 13, 8, 4, 0, 1, 9, 5,
				  10, 6, 2, 3, 14, 7, 11, 15,
				  19, 18, 23, 27, 31, 30, 22, 26,
				  53, 57, 61, 60, 49, 56, 52, 48,
				  44, 45, 40, 36, 32, 33, 41, 37,
				  42, 38, 34, 35, 46, 39, 43, 47,
				  51, 50, 55, 59, 63, 62, 54, 58 };

  if (ddlno == 0 || ddlno == 1)
    {
      um  = mcmno/72;
      Int_t mcmnonew = mcmno - 72*um;
      Int_t rowcol  = kCh[chno];
      Int_t irownew = rowcol/4;
      Int_t icolnew = rowcol%4;
      
      remmcm  = mcmnonew%12;
      divmcm  = mcmnonew/12;
      
      row = 16*divmcm + irownew;
      col =  4*remmcm + icolnew;
      // This obtatined row and col (0,0) are the top left corner.
      // Needs transformation to get the Geant (0,0)

    }
  else   if (ddlno == 2 || ddlno == 3)
    {
      um  = mcmno/72;
      Int_t mcmnonew = mcmno - 72*um;
      Int_t rowcol  = kCh[chno];
      Int_t irownew = rowcol/4;
      Int_t icolnew = rowcol%4;
      
      remmcm  = mcmnonew%24;
      divmcm  = mcmnonew/24;
      
      row = 16*divmcm + irownew;
      col =  4*remmcm + icolnew;
      // This obtatined row and col (0,0) are the top left corner.
      // Needs transformation to get the Geant (0,0)

    }
  else if (ddlno == 4 || ddlno == 5)
    {
      um  = mcmno/72;
      Int_t mcmnonew = mcmno - 72*um;
      Int_t rowcol  = kCh[chno];
      Int_t irownew = rowcol/4;
      Int_t icolnew = rowcol%4;
      
      if (um < 6)
	{
	  remmcm  = mcmnonew%12;
	  divmcm  = mcmnonew/12;
	}
      else if (um >= 6)
	{
	  remmcm  = mcmnonew%24;
	  divmcm  = mcmnonew/24;
	}

      row = 16*divmcm + irownew;
      col =  4*remmcm + icolnew;
      // This obtatined row and col (0,0) are the top left corner.
      // Needs transformation to get the Geant (0,0)

    }

}
