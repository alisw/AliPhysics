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
/// This is a base class for reading raw data digits in Altro format.
/// The class is able to read the RCU v3 and above formats.
/// The main difference between the format V3 and older ones is in
/// the coding of the 10-bit Altro payload words. In V3 3 10-bit words
/// are coded in one 32-bit word. The bits 30 and 31 are used to identify
/// the payload, altro header and RCU trailer contents.
///
///
/// cvetan.cheshkov@cern.ch 1/04/2009
///////////////////////////////////////////////////////////////////////////////

#include "AliAltroRawStreamV3.h"
#include "AliRawReader.h"
#include "AliLog.h"

ClassImp(AliAltroRawStreamV3)


//_____________________________________________________________________________
AliAltroRawStreamV3::AliAltroRawStreamV3(AliRawReader* rawReader) :
  fIsShortDataHeader(kFALSE),
  fDDLNumber(-1),
  fRCUId(-1),
  fHWAddress(-1),
  fRawReader(rawReader),
  fData(NULL),
  fPosition(-1),
  fCount(-1),
  fStartTimeBin(-1),
  fBunchLength(-1),
  fBadChannel(kFALSE),
  fPayloadSize(-1),
  fBunchDataPointer(NULL),
  fBunchDataIndex(-1),
  fRCUTrailerData(NULL),
  fRCUTrailerSize(0),
  fFECERRA(0),
  fFECERRB(0),
  fERRREG2(0),
  fERRREG3(0),
  fERRREG4(0),
  fActiveFECsA(0),
  fActiveFECsB(0),
  fAltroCFG1(0),
  fAltroCFG2(0)
{
  // Constructor
  // Create an object to read Altro raw digits in
  // RCU version 3 and beyond format
  for(Int_t i = 0; i < kMaxNTimeBins; i++) fBunchData[i] = 0;
}

//_____________________________________________________________________________
AliAltroRawStreamV3::~AliAltroRawStreamV3()
{
// destructor
// nothing to do
}

//_____________________________________________________________________________
AliAltroRawStreamV3::AliAltroRawStreamV3(const AliAltroRawStreamV3& stream) :
  TObject(stream),
  fIsShortDataHeader(stream.fIsShortDataHeader),
  fDDLNumber(stream.fDDLNumber),
  fRCUId(stream.fRCUId),
  fHWAddress(stream.fHWAddress),
  fRawReader(stream.fRawReader),
  fData(stream.fData),
  fPosition(stream.fPosition),
  fCount(stream.fCount),
  fStartTimeBin(stream.fStartTimeBin),
  fBunchLength(stream.fBunchLength),
  fBadChannel(stream.fBadChannel),
  fPayloadSize(stream.fPayloadSize),
  fBunchDataPointer(stream.fBunchDataPointer),
  fBunchDataIndex(stream.fBunchDataIndex),
  fRCUTrailerData(stream.fRCUTrailerData),
  fRCUTrailerSize(stream.fRCUTrailerSize),
  fFECERRA(stream.fFECERRA),
  fFECERRB(stream.fFECERRB),
  fERRREG2(stream.fERRREG2),
  fERRREG3(stream.fERRREG3),
  fERRREG4(stream.fERRREG4),
  fActiveFECsA(stream.fActiveFECsA),
  fActiveFECsB(stream.fActiveFECsB),
  fAltroCFG1(stream.fAltroCFG1),
  fAltroCFG2(stream.fAltroCFG2)
{
  // Copy constructor
  // Copy the bunch data array
  for(Int_t i = 0; i < kMaxNTimeBins; i++) fBunchData[i] = stream.fBunchData[i];
}

//_____________________________________________________________________________
AliAltroRawStreamV3& AliAltroRawStreamV3::operator = (const AliAltroRawStreamV3& stream)
{
  // assignment operator
  // ... 
  if(&stream == this) return *this;

  fIsShortDataHeader = stream.fIsShortDataHeader;
  fDDLNumber         = stream.fDDLNumber;
  fRCUId             = stream.fRCUId;
  fHWAddress         = stream.fHWAddress;
  fRawReader         = stream.fRawReader;
  fData              = stream.fData;
  fPosition          = stream.fPosition;
  fCount             = stream.fCount;
  fStartTimeBin      = stream.fStartTimeBin;
  fBunchLength       = stream.fBunchLength;
  fBadChannel        = stream.fBadChannel;
  fPayloadSize       = stream.fPayloadSize;
  fBunchDataPointer  = stream.fBunchDataPointer;
  fBunchDataIndex    = stream.fBunchDataIndex;
  fRCUTrailerData    = stream.fRCUTrailerData;
  fRCUTrailerSize    = stream.fRCUTrailerSize;
  fFECERRA           = stream.fFECERRA;
  fFECERRB           = stream.fFECERRB;
  fERRREG2           = stream.fERRREG2;
  fERRREG3           = stream.fERRREG3;
  fERRREG4           = stream.fERRREG4;
  fActiveFECsA       = stream.fActiveFECsA;
  fActiveFECsB       = stream.fActiveFECsB;
  fAltroCFG1         = stream.fAltroCFG1;
  fAltroCFG2         = stream.fAltroCFG2;

  for(Int_t i = 0; i < kMaxNTimeBins; i++) fBunchData[i] = stream.fBunchData[i];

  return *this;
}

//_____________________________________________________________________________
void AliAltroRawStreamV3::Reset()
{
// Complete reset of raw stream params
// Reset of the raw-reader as well

  fDDLNumber = fRCUId = fHWAddress = -1;
  fPosition = fCount = -1;
  fBunchLength = fStartTimeBin = -1;
  fBadChannel = kFALSE;
  fPayloadSize = -1;
  fBunchDataPointer = NULL;
  fBunchDataIndex = -1;

  fRCUTrailerData = NULL;
  fRCUTrailerSize = 0;

  fFECERRA = fFECERRB = fERRREG2 = fERRREG3 = fERRREG4 = fActiveFECsA = fActiveFECsB = fAltroCFG1 = fAltroCFG2 = 0;

  if (fRawReader) fRawReader->Reset();

}

//_____________________________________________________________________________
Bool_t AliAltroRawStreamV3::NextDDL()
{
// Read the next DDL payload (CDH + RCU trailer)
// Updates the information which is coming from these
// two sources

  fPosition = 0;
  // Get next DDL payload
  // return wtih false in case no more data payloads
  // are found
  do {
    if (!fRawReader->ReadNextData(fData)) return kFALSE;
  } while (fRawReader->GetDataSize() == 0);

  fDDLNumber = fRawReader->GetDDLID();

  UChar_t rcuVer = fRawReader->GetBlockAttributes();

  if (!ReadRCUTrailer(rcuVer)) return kFALSE;

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliAltroRawStreamV3::NextChannel()
{
  // Read the next Altro channel from the
  // raw-data stream
  // Updates the channel hardware address member and
  // channel data size. Sets the error flag in case
  // RCU signals readout error in this channel
  fCount = -1;
  fBadChannel = kFALSE;
  fBunchDataIndex = 0;

  UInt_t word = 0;
  do {
    word = Get32bitWord(fPosition++);
    if (fPosition > fPayloadSize) return kFALSE;
  }
  while ((word >> 30) != 1);

  // check for readout errors
  fBadChannel = (word >> 29) & 0x1;

  // extract channel payload and hw address
  fCount = (word >> 16) & 0x3FF; 
  fHWAddress = word & 0xFFF;

  // Now unpack the altro data
  // Revert the order of the samples
  // inside the bunch so that the
  // first time is first in the samples
  // array
  Int_t isample = 0;
  Int_t nwords = (fCount+2)/3;
  for (Int_t iword = 0; iword < nwords; iword++) {
    word = Get32bitWord(fPosition++);
    if ((word >> 30) != 0) {
      // Unexpected end of altro channel payload
      AliWarning(Form("Unexpected end of payload in altro channel payload! Address=0x%x, word=0x%x",
		      fHWAddress,word));
      fRawReader->AddMinorErrorLog(kAltroPayloadErr,Form("hw=0x%x",fHWAddress));
      fCount = -1;
      fPosition--;
      return kFALSE;
    }
    fBunchData[isample++] = (word >> 20) & 0x3FF;
    fBunchData[isample++] = (word >> 10) & 0x3FF;
    fBunchData[isample++] = word & 0x3FF;
  }  

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliAltroRawStreamV3::NextBunch()
{
  // Read next altro bunch from the
  // raw-data stream.
  // Updates the start/end time-bins
  // and the array with altro samples
  fBunchLength = fStartTimeBin = -1;
  fBunchDataPointer = NULL;

  if ((fBunchDataIndex >= fCount) || fBadChannel) return kFALSE;

  fBunchLength = fBunchData[fBunchDataIndex];
  if ((fBunchDataIndex + fBunchLength) > fCount) {
    // Too long bunch detected
    AliWarning(Form("Too long bunch detected in Address=0x%x ! Expected <= %d 10-bit words, found %d !",
		    fHWAddress,fCount-fBunchDataIndex,fBunchLength));
    fRawReader->AddMinorErrorLog(kAltroBunchHeadErr,Form("hw=0x%x",fHWAddress));
    fCount = fBunchLength = -1;
    return kFALSE;
  }
  fBunchDataIndex++;
  fBunchLength -= 2;

  fStartTimeBin = fBunchData[fBunchDataIndex++];

  fBunchDataPointer = &fBunchData[fBunchDataIndex];

  fBunchDataIndex += fBunchLength;

  return kTRUE;
}

//_____________________________________________________________________________
UInt_t AliAltroRawStreamV3::Get32bitWord(Int_t index) const
{
  // This method returns the 32 bit word at a given
  // position inside the raw data payload.
  // The 'index' points to the beginning of the word.
  // The method is supposed to be endian (platform)
  // independent.

  index = (index << 2); 
  UInt_t word  = 0;
  word |= fData[index++];
  word |= fData[index++] << 8;
  word |= fData[index++] << 16;
  word |= fData[index++] << 24;

  return word;
}

///_____________________________________________________________________________
Bool_t AliAltroRawStreamV3::ReadRCUTrailer(UChar_t rcuVer)
{
  // Read the RCU trailer according
  // to the RCU formware version
  // specified in CDH
  // Cross-check with version found in the
  // trailer

  fRCUTrailerData = NULL;
  fRCUTrailerSize = 0;
  fPayloadSize = -1;

  Int_t index = fRawReader->GetDataSize()/4;

  // First read 32-bit word with the
  // trailer size (7 bits), RCU ID (9 bits) and
  // RCU firmware version (8 bits?)
  // The two major bit should be 11 (identifies
  // the end of the trailer)
  UInt_t word = Get32bitWord(--index);

  if ((word >> 30) != 3) {
    fRawReader->AddFatalErrorLog(kRCUTrailerErr,"");
    AliError("Last RCU trailer word not found!");
    return kFALSE;
  }

  UChar_t ver = (word >> 16) & 0xFF;
  if (ver != rcuVer) {
    // Wrong RCU firmware version detected
    fRawReader->AddMajorErrorLog(kRCUVerErr,Form("%d!=%d",
						 ver,rcuVer));
    AliDebug(1,Form("Wrong RCU firmware version detected: %d != %d",
		    ver,rcuVer));
  }

  fRCUId = (Int_t)((word >> 7) & 0x1FF);
  Int_t trailerSize = (word & 0x7F);

  // Now read the beginning of the trailer
  // where the payload size is written
  if (trailerSize < 2) {
    fRawReader->AddMajorErrorLog(kRCUTrailerErr,Form("tr=%d bytes",
						     trailerSize*4));
    AliWarning(Form("Invalid trailer size found (%d bytes) !",
		    trailerSize*4));
    return kFALSE;
  }

  trailerSize -= 2;
  fRCUTrailerSize = trailerSize*4;

  for (; trailerSize > 0; trailerSize--) {
    word = Get32bitWord(--index);
    if ((word >> 30) != 2) {
      fRawReader->AddMinorErrorLog(kRCUTrailerErr,"missing 10");
      AliWarning("Missing RCU trailer identifier pattern!");
      continue;
    }
    Int_t parCode = (word >> 26) & 0xF;
    Int_t parData = word & 0x3FFFFFF;
    switch (parCode) {
    case 1:
      // ERR_REG1
      fFECERRA = ((parData >> 13) & 0x1FFF) << 7;
      fFECERRB = ((parData & 0x1FFF)) << 7;
      break;
    case 2:
      // ERR_REG2
      fERRREG2 = parData & 0x1FF;
      break;
    case 3:
      // ERR_REG3
      fERRREG3 = parData & 0xFFF;
      break;
    case 4:
      // ERR_REG4
      fERRREG4 = parData & 0xFFF;
      break;
    case 5:
      // FEC_RO_A
      fActiveFECsA = parData & 0xFFFF;
      break;
    case 6:
      // FEC_RO_B
      fActiveFECsB = parData & 0xFFFF;
      break;
    case 7:
      // RDO_CFG1
      fAltroCFG1 = parData & 0xFFFFF;
      break;
    case 8:
      // RDO_CFG2
      fAltroCFG2 = parData & 0x1FFFFFF;
     break;
    default:
      fRawReader->AddMinorErrorLog(kRCUTrailerErr,"undef word");
      AliWarning(Form("Undefined parameter code %d, ignore it !",
		      parCode));
      break;
    }
  }

  if (index < 1) {
    fRawReader->AddMajorErrorLog(kRCUTrailerErr,Form("tr=%d raw=%d bytes",
						     fRCUTrailerSize+8,
						     fRawReader->GetDataSize()));
    AliWarning(Form("Invalid trailer size found (%d bytes) ! The size is bigger than the raw data size (%d bytes)!",
		    fRCUTrailerSize,
		    fRawReader->GetDataSize()));
  }

  fRCUTrailerData = fData + index*4;

  // Now read the payload size
  // (First word in the RCU trailer)
  fPayloadSize = Get32bitWord(--index) & 0x3FFFFFF;

  if ((fRawReader->GetDataSize() - fRCUTrailerSize - 8) != (fPayloadSize*4)) {
    fRawReader->AddMajorErrorLog(kRCUTrailerSizeErr,Form("h=%d tr=%d rcu=%d bytes",
							 fRawReader->GetDataSize(),
							 fRCUTrailerSize+8,
							 fPayloadSize*4));
    AliWarning(Form("Inconsistent raw data size ! Raw data size - %d bytes (from CDH), RCU trailer - %d bytes, raw data size (from RCU trailer) - %d bytes !",
		    fRawReader->GetDataSize(),
		    fRCUTrailerSize+8,
		    fPayloadSize*4));
  }

  return kTRUE;
}

//_____________________________________________________________________________
Int_t AliAltroRawStreamV3::GetBranch() const
{
  // The method provides the RCU branch index (0 or 1)
  // for the current hardware address.
  // In case the hardware address has not been yet
  // initialized, the method returns -1
  if (fHWAddress == -1) return -1;

  return ((fHWAddress >> 11) & 0x1);
}

//_____________________________________________________________________________
Int_t AliAltroRawStreamV3::GetFEC() const
{
  // The method provides the front-end card index
  // for the current hardware address.
  // In case the hardware address has not been yet
  // initialized, the method returns -1
  if (fHWAddress == -1) return -1;

  return ((fHWAddress >> 7) & 0xF);
}

//_____________________________________________________________________________
Int_t AliAltroRawStreamV3::GetAltro() const
{
  // The method provides the altro chip index
  // for the current hardware address.
  // In case the hardware address has not been yet
  // initialized, the method returns -1
  if (fHWAddress == -1) return -1;

  return ((fHWAddress >> 4) & 0x7);
}

//_____________________________________________________________________________
Int_t AliAltroRawStreamV3::GetChannel() const
{
  // The method provides the channel index
  // for the current hardware address.
  // In case the hardware address has not been yet
  // initialized, the method returns -1
  if (fHWAddress == -1) return -1;

  return (fHWAddress & 0xF);
}

//_____________________________________________________________________________
Bool_t AliAltroRawStreamV3::GetRCUTrailerData(UChar_t*& data) const
{
  // Return a pointer to the RCU trailer
  // data. Should be called always after
  // the RCU trailer was already processed
  // in the GetPosition() method
  if (!fRCUTrailerSize || !fRCUTrailerData) {
    AliError("No valid RCU trailer data is found !");
    data = NULL;
    return kFALSE;
  }

  data = fRCUTrailerData;

  return kTRUE;
}

//_____________________________________________________________________________
void AliAltroRawStreamV3::PrintRCUTrailer() const
{
  // Prints the contents of
  // the RCU trailer data
  printf("RCU trailer:\n===========\n");
  printf("FECERRA: 0x%x\nFECERRB: 0x%x\n",fFECERRA,fFECERRB);
  printf("ERRREG2: 0x%x\n",fERRREG2);
  printf("#channels skipped due to address mismatch: %d\n",GetNChAddrMismatch());
  printf("#channels skipped due to bad block length: %d\n",GetNChLengthMismatch());
  printf("Active FECs (branch A): 0x%x\nActive FECs (branch B): 0x%x\n",fActiveFECsA,fActiveFECsB);
  printf("Baseline corr: 0x%x\n",GetBaselineCorr());
  printf("Number of presamples: %d\nNumber of postsamples: %d\n",GetNPresamples(),GetNPostsamples());
  printf("Second baseline corr: %d\n",GetSecondBaselineCorr());
  printf("GlitchFilter: %d\n",GetGlitchFilter());
  printf("Number of non-ZS postsamples: %d\nNumber of non-ZS presamples: %d\n",GetNNonZSPostsamples(),GetNNonZSPresamples());
  printf("Number of ALTRO buffers: %d\n",GetNAltroBuffers());
  printf("Number of pretrigger samples: %d\n",GetNPretriggerSamples());
  printf("Number of samples per channel: %d\n",GetNSamplesPerCh());
  printf("Sparse readout: %d\n",GetSparseRO());
  printf("Sampling time: %e s\n",GetTSample());
  printf("L1 Phase: %e s\n",GetL1Phase());
  printf("AltroCFG1: 0x%x\nAltroCFG2: 0x%x\n",GetAltroCFG1(),GetAltroCFG2());
  printf("===========\n");
}

//_____________________________________________________________________________
void AliAltroRawStreamV3::SelectRawData(Int_t detId)
{
  // Select the raw data for specific
  // detector id
  AliDebug(1,Form("Selecting raw data for detector %d",detId));
  fRawReader->Select(detId);
}

//_____________________________________________________________________________
void AliAltroRawStreamV3::SelectRawData(const char *detName)
{
  // Select the raw data for specific
  // detector name
  AliDebug(1,Form("Selecting raw data for detector %s",detName));
  fRawReader->Select(detName);
}

//_____________________________________________________________________________
Double_t AliAltroRawStreamV3::GetTSample() const
{
  // Returns the sampling time
  // in seconds. In case the rcu trailer
  // was note read, return an invalid number (0)

  if (!fRCUTrailerData) return 0.;

  const Double_t kLHCTimeSample = 25.0e-9; // LHC clocks runs at 40 MHz
  UChar_t fq = (fAltroCFG2 >> 5) & 0xF;
  Double_t tSample;
  switch (fq) {
  case 0:
    // 20 MHz
    tSample = 2.0*kLHCTimeSample;
    break;
  case 1:
    // 10 Mhz
    tSample = 4.0*kLHCTimeSample;
    break;
  case 2:
    // 5 MHz
    tSample = 8.0*kLHCTimeSample;
    break;
  default:
    AliWarning(Form("Invalid sampling frequency value %d !",
		      fq));
    tSample = 0.;
    break;
  }

  return tSample;
}

//_____________________________________________________________________________
Double_t AliAltroRawStreamV3::GetL1Phase() const
{
  // Returns the L1 phase w.r.t to the
  // LHC clock
  if (!fRCUTrailerData) return 0.;

  const Double_t kLHCTimeSample = 25.0e-9; // LHC clocks runs at 40 MHz
  Double_t phase = ((Double_t)(fAltroCFG2 & 0x1F))*kLHCTimeSample;

  Double_t tSample = GetTSample();
  if (phase >= tSample) {
    AliWarning(Form("Invalid L1 trigger phase (%f >= %d) !",
		    phase,tSample));
    phase = 0.;
  }

  return phase;
}

//_____________________________________________________________________________
void AliAltroRawStreamV3::AddMappingErrorLog(const char *message)
{
  // Signal a mapping error
  // The method can be used by the TPC,PHOS,EMCAL,FMD raw stream
  // classes in order to log an error related to bad altro mapping

  if (fRawReader) fRawReader->AddMinorErrorLog(kBadAltroMapping,message);
}
