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
  ConvertDDL2SMN(iddl, ium, fSMN, fModule, fDetector);
  TransformH2S(fSMN, fRow, fColumn);

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
//_____________________________________________________________________________
void AliPMDRawStream::ConvertDDL2SMN(Int_t iddl, Int_t ium, Int_t &smn,
				    Int_t &module, Int_t &detector) const
{
  // This converts the DDL number to Module Number which runs from 0-47
  // Serial module number in one detector which runs from 0-23
  // Also gives the detector number (0:PRE plane, 1:CPV plane)
  if (iddl < 4)
    {
      module = iddl*6 + ium;
      detector = 0;
      smn = iddl*6 + ium;
    }
  else if (iddl == 4)
    {
      if (ium < 6)
	{
	  module = 24 + ium;
	  smn    = ium;
	}
      else if (ium >= 6)
	{
	  module = 30 + ium;
	  smn    = 6 + ium;
	}
      detector = 1;
    }
  else if (iddl == 5)
    {
      
      if (ium < 6)
	{
	  module = 30 + ium;
	  smn    = 6 + ium;
	}
      else if (ium >= 6)
	{
	  module = 36 + ium;
	  smn    = 12 + ium;
	}
      detector = 1;
    }
}
//_____________________________________________________________________________
void AliPMDRawStream::TransformH2S(Int_t smn, Int_t &row, Int_t &col) const
{
  // Transform the Hardware (0,0) coordinate to Software (0,0) coordinate
  // and also writes in the digit form
  // i.e., For SuperModule 1 &2, instead of 96x48 it is 48x96
  // For Supermodule 3 & 4, 48x96

  Int_t irownew1 = 0;
  Int_t icolnew1 = 0;
  Int_t irownew  = 0;
  Int_t icolnew  = 0;
  // Transform all the (0,0) coordinates to the geant frame
  if(smn < 6)
    {
      irownew1 = 95 - row;
      icolnew1 = col;
    }
  else if(smn >= 6 && smn < 12)
    {
      irownew1 = row;
      icolnew1 = 47 - col;
    }
  else if(smn >= 12 && smn < 18)
    {
      irownew1 = 47 - row;
      icolnew1 = col;
    }
  else if(smn >= 18 && smn < 24)
    {
      irownew1 = row;
      icolnew1 = 95 - col;
    }
  
  // for smn < 12          : row = 96, column = 48
  // for smn>= 12 and < 24 : row = 48, column = 96
  // In order to make it uniform dimension, smn < 12 are inverted
  // i.e., row becomes column and column becomes row
  // for others it remains same
  // This is further inverted back while calculating eta and phi
  if(smn < 12)
    {
      // SupeModule 1 and 2 : Rows are inverted to columns and vice versa
      // and at the time of calculating the eta,phi it is again reverted
      // back
      irownew = icolnew1;
      icolnew = irownew1;
    }
  else if( smn >= 12 && smn < 24)
    {
      irownew = irownew1;
      icolnew = icolnew1;
    }

  row = irownew;
  col = icolnew;
}
//_____________________________________________________________________________
