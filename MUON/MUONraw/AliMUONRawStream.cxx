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


//-----------------------------------------------------------------------------
/// \class AliMUONRawStream
/// This base class to MUON raw stream
///
///
/// \author Christian Finck
//-----------------------------------------------------------------------------

#include "AliMUONRawStream.h"

#include "AliRawReader.h"


/// \cond CLASSIMP
ClassImp(AliMUONRawStream)
/// \endcond

//___________________________________________
AliMUONRawStream::AliMUONRawStream()
 : TObject(),
   fRawReader(),
   fEnableErrorLogger(kFALSE)
{
  ///
  /// Default ctor for monitoring purposes
  ///
  
  
}

//_________________________________________________________________
AliMUONRawStream::AliMUONRawStream(AliRawReader* rawReader)
:  TObject(),
   fRawReader(rawReader),
   fEnableErrorLogger(kFALSE)
{
  ///
  /// ctor with AliRawReader as argument
  /// for reconstruction purpose
  ///
  
  
}

//___________________________________
AliMUONRawStream::~AliMUONRawStream()
{
  ///
  /// clean up
  ///
}

//_________________________________________________________________
void AliMUONRawStream::Swap(UInt_t* buffer, Int_t size) const
{
  /// swap from little to big endian
 
    RawWord *word, temp;
    word = (RawWord *) buffer;

    for (int i = 0 ; i < size; i++) {
      temp = *(((RawWord*)buffer)+i);
      word->fB1 = temp.fB4;
      word->fB2 = temp.fB3;
      word->fB3 = temp.fB2;
      word->fB4 = temp.fB1;
      word++;
    }

}                
