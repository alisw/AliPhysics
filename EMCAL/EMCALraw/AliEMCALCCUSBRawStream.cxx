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

#include "AliEMCALCCUSBRawStream.h"

// AliRoot system
#include "AliRawReader.h"
#include "AliLog.h"

/// \cond CLASSIMP
ClassImp(AliEMCALCCUSBRawStream) ;
/// \endcond

///
/// Constructor
//_______________________________________________________________________
AliEMCALCCUSBRawStream::AliEMCALCCUSBRawStream(AliRawReader* rawReader) :
  fRawReader(rawReader),
  fData(0),
  fHeader(0),
  fOptHeader(0),
  fEventLength(0),
  fEOBuffer(0)
{
  fRawReader = rawReader;
  
  fRawReader->Reset();
  fRawReader->SelectEquipment(1, 1, 1);
  
  for(Int_t i = 0 ; i < fgkNTDC; i++ )          fTDC[i] = 0 ;             
  for(Int_t i = 0 ; i < fgkNQDC; i++ )          fQDC[i] = 0 ;                         
  for(Int_t i = 0 ; i < fgkNScalerCCUSB; i++ )  fScalerCCUSB[i]  = 0;          
  for(Int_t i = 0 ; i < fgkNScalerLecroy; i++ ) fScalerLecroy[i] = 0;       
}

///
/// Read the next raw digit.
/// \return kFALSE if there is no digit left
//_______________________________________________________________________
Bool_t AliEMCALCCUSBRawStream::Next()
{
  
  if ( fEOBuffer == 0xFFFF ) { fEOBuffer = 0; return kFALSE; }
  
  if (!fRawReader->ReadNextInt((UInt_t&) fHeader)) 
  {
    AliError("No header");
    return kFALSE;
  }
  
  if (!fRawReader->ReadNextInt((UInt_t&) fOptHeader)) 
  {
    AliError("No optional header");
    return kFALSE;
  }
  
  if (!fRawReader->ReadNextInt((UInt_t&) fEventLength)) 
  {
    AliError("No event length");
    return kFALSE;
  }
  
  for (Int_t i = 0; i < fgkNScalerCCUSB; i++) 
  {  
    if (!fRawReader->ReadNext((UChar_t*)&fData,8)) 
    {
      AliError("Internal CC-USB scaler issing");
      return kFALSE;
    }
    
    fScalerCCUSB[i] = fData;
  }
  
  for (Int_t i = 0; i < fgkNScalerLecroy; i++) 
  {  
    if (!fRawReader->ReadNext((UChar_t*)&fData,8)) 
    {
      AliError("Lecroy scaler missing");
      return kFALSE;
    }
    
    fScalerLecroy[i] = fData;
  }
  
  for (Int_t i = 0; i < fgkNTDC; i++) 
  {  
    if (!fRawReader->ReadNextInt(fData)) 
    {
      AliError("Incomplete TDC equipment");
      return kFALSE;
    }
    
    fTDC[i] = fData;
  }
  
  for (Int_t i = 0; i < fgkNQDC; i++) 
  {  
    if (!fRawReader->ReadNextInt(fData)) 
    {
      AliError("Incomplete QDC equipment");
      return kFALSE;
    }
    
    fQDC[i] = fData;
  }
  
  if ( !fRawReader->ReadNextInt((UInt_t&) fEOBuffer) ) 
  {
    AliError("No end of buffer");
    return kFALSE;
  }
  
  return kTRUE;
}

