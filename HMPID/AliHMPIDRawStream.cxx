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
///
/// This is a class for reading the HMPID raw data
/// The format of the raw data corresponds to the one
/// which was documented by Paolo Martinengo.
///
///////////////////////////////////////////////////////////////////////////////

#include "AliHMPIDRawStream.h"
#include "AliRawReader.h"
#include "AliLog.h"
#include "AliDAQ.h"

ClassImp(AliHMPIDRawStream)

//_____________________________________________________________________________
AliHMPIDRawStream::AliHMPIDRawStream(AliRawReader* rawReader) :
  fDDLNumber(-1),
  fRawReader(rawReader),
  fData(NULL),
  fPosition(-1)
{
  // Constructor
  Init();

  fRawReader->Reset();
  fRawReader->Select("HMPID");
}

//_____________________________________________________________________________
AliHMPIDRawStream::~AliHMPIDRawStream()
{
  // destructor
}

//_____________________________________________________________________________
void AliHMPIDRawStream::Init()
{
  // Initalize the container
  // with the pad charges
  for(Int_t i = 0; i < kNRows; i++)
    for(Int_t j = 0; j < kNDILOGICAdd; j++)
      for(Int_t k = 0; k < kNPadAdd; k++)
	fCharge[i][j][k] = -1;
}

//_____________________________________________________________________________
void AliHMPIDRawStream::Reset()
{
  // reset raw stream params

  // Reinitalize the containers
  Init();

  fDDLNumber = -1;
  fPosition = -1;
  fData = NULL;

  if (fRawReader) fRawReader->Reset();
}

//_____________________________________________________________________________
Bool_t AliHMPIDRawStream::Next()
{
  // read next DDL raw data from the VZERO raw data stream
  // return kFALSE in case of error or no data left

  do {
    if (!fRawReader->ReadNextData(fData)) return kFALSE;
  } while (fRawReader->GetDataSize() == 0);

  if (fRawReader->GetDataSize() != 47148) {
    fRawReader->AddFatalErrorLog(kRawDataSizeErr,Form("size %d != 47148",fRawReader->GetDataSize()));
    AliWarning(Form("Wrong HMPID raw data size: %d, expected 47148 bytes!",fRawReader->GetDataSize()));
    return kTRUE;
  }
  fDDLNumber = fRawReader->GetDDLID();
  fPosition = 0;

  Init();

  // Look over rows
  for(Int_t iRow = 1; iRow < kNRows; iRow++) {
    // Read row marker
    UInt_t rowMarker = GetNextWord() & 0x1ffffff;
    if (rowMarker != 0x1ea32a8) {
      fRawReader->AddMajorErrorLog(kRowMarkerErr);
      AliWarning(Form("Wrong row marker %x for row %d, expected 0xx1ea32a8!",rowMarker,iRow));
      return kTRUE;
    }
    UInt_t dilogic = 0, row = 0;
    for(Int_t iDILOGIC = 1; iDILOGIC < kNDILOGICAdd; iDILOGIC++) {
      // Read pad charges
      for(Int_t iPad = 0; iPad < kNPadAdd; iPad++) {
	UInt_t data = GetNextWord();
	row = (data >> 22) & 0x1f;
	if (row < 1 || row >= kNRows) {
	  fRawReader->AddMajorErrorLog(kWrongRowErr,Form("row %d",row));
	  AliWarning(Form("Wrong row index: %d, expected (1 -> %d)!",row,kNRows));
	  row = iRow;
	}
	dilogic = (data >> 18) & 0xf;
	if (dilogic < 1 || dilogic >= kNDILOGICAdd) {
	  fRawReader->AddMajorErrorLog(kWrongDilogicErr,Form("dil %d",dilogic));
	  AliWarning(Form("Wrong DILOGIC index: %d, expected (1 -> %d)!",dilogic,kNDILOGICAdd));
	  dilogic = iDILOGIC;
	}
	UInt_t pad = (data >> 12) & 0x3f;
	if (pad >= kNPadAdd) {
	  fRawReader->AddMajorErrorLog(kWrongPadErr,Form("pad %d",pad));
	  AliWarning(Form("Wrong pad index: %d, expected (0 -> %d)!",pad,kNPadAdd));
	  pad = iPad;
	}
	fCharge[row][dilogic][pad] = data & 0xfff;
      }
      // Now read the end-of-event word
      UInt_t eOfEvent = GetNextWord() & 0xfffffff;
      if (!((eOfEvent >> 27) & 0x1)) {
      	fRawReader->AddMajorErrorLog(kEoEFlagErr);
      	AliWarning(Form("Missing end-of-event flag! (%x)",eOfEvent));
      	return kTRUE;
      }
      UInt_t wc = eOfEvent & 0x7f;
      if (wc != 48) {
	fRawReader->AddMajorErrorLog(kEoESizeErr,Form("eoe size=%d",wc));
	AliWarning(Form("Wrong end-of-event word-count:%d, expected 48!",wc));
	return kTRUE;
      }
      UInt_t da = (eOfEvent >> 18) & 0xf;
      if (da != dilogic) {
	fRawReader->AddMajorErrorLog(kEoEDILOGICErr,Form("eoe dil %d != %d",da,dilogic));
	AliWarning(Form("Wrong DILOGIC address found in end-of-event: %d, expected %d!",da,dilogic));
	return kTRUE;
      }
      UInt_t ca = (eOfEvent >> 22) & 0x1f;
      if (ca != row) {
	fRawReader->AddMajorErrorLog(kEoERowErr,Form("eoe row %d != %d",ca,row));
	AliWarning(Form("Wrong row index found in end-of-event: %d, expected %d!",ca,row));
	return kTRUE;
      }
    }

    // Read the segment marker
    // One maker per 8 rows
    if (iRow%8 == 0) {
      UInt_t segWord = GetNextWord();
      if ((segWord >> 8) != 0xab0f58) {
	fRawReader->AddMajorErrorLog(kBadSegWordErr);
	AliWarning(Form("Wrong segment word signature: %x, expected 0xab0f58!",(segWord >> 8)));
	return kTRUE;
      }
      if ((segWord & 0xff) != (((UInt_t)iRow + 7) / 8)) {
	fRawReader->AddMajorErrorLog(kWrongSegErr,Form("seg %d != %d",segWord & 0xff,(iRow + 7) / 8));
	AliWarning(Form("Segment index (%d) does not correspond to the one expected from row index (%d)!",segWord & 0xff,(iRow + 7) / 8));
	return kTRUE;
      }
    }
  }    
  return kTRUE;
}

//_____________________________________________________________________________
UInt_t AliHMPIDRawStream::GetNextWord()
{
  // This method returns the next 32 bit word
  // inside the raw data payload.
  // The method is supposed to be endian (platform)
  // independent.
  if (!fData || fPosition < 0) AliFatal("Raw data payload buffer is not yet initialized !");

  UInt_t word = 0;
  word |= fData[fPosition++];
  word |= fData[fPosition++] << 8;
  word |= fData[fPosition++] << 16;
  word |= fData[fPosition++] << 24;

  return word;
}

//_____________________________________________________________________________
Short_t AliHMPIDRawStream::GetCharge(Int_t row, Int_t dilogic, Int_t pad) const
{
  // The method returns the charge collected
  // in a particular channel
  // Return -1 in case the charge from the channels
  // has not been read or invalid arguments

  if (row < 1 || row >= kNRows) {
    AliError(Form("Wrong row index %d!",row));
    return 0;
  }

  if (dilogic < 1 || dilogic >= kNDILOGICAdd) {
    AliError(Form("Wrong DILOGIC address %d!",dilogic));
    return 0;
  }
  
  if (pad >= kNPadAdd) {
    AliError(Form("Wrong pad index %d!",pad));
    return 0;
  }

  return fCharge[row][dilogic][pad];
}
