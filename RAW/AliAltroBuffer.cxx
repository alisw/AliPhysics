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

// Interface to the Altro format
// to read and write digits
// To be used in Alice Data Challenges 
// and in the compression of the RAW data

#include "AliAltroBuffer.h"
#include "AliAltroMapping.h"
#include "AliRawDataHeaderSim.h"
#include "AliLog.h"
#include "AliFstream.h"
//#include <stdlib.h>


ClassImp(AliAltroBuffer)

//_____________________________________________________________________________
AliAltroBuffer::AliAltroBuffer(const char* fileName, const AliAltroMapping *mapping):
  fShift(0),
  fCurrentCell(0),
  fFreeCellBuffer(16),
  fVerbose(0),
  fFile(NULL),
  fDataHeaderPos(0),
  fMapping(mapping)
{
  //the buffer is cleaned 
  for (Int_t i = 0; i < 5; i++) fBuffer[i] = 0;

  //open the output file
  fFile = new AliFstream(fileName);

}

//_____________________________________________________________________________
AliAltroBuffer::~AliAltroBuffer()
{
// destructor

  //Flush out the Buffer content at the end only if Buffer wasn't completely filled
  Flush();
  if (fVerbose) Info("~AliAltroBuffer", "File Created");

  delete fFile;

}

//_____________________________________________________________________________
AliAltroBuffer::AliAltroBuffer(const AliAltroBuffer& source):
  TObject(source),
  fShift(source.fShift),
  fCurrentCell(source.fCurrentCell),
  fFreeCellBuffer(source.fFreeCellBuffer),
  fVerbose(source.fVerbose),
  fFile(NULL),
  fDataHeaderPos(source.fDataHeaderPos),
  fMapping(source.fMapping)
{
// Copy Constructor

  Fatal("AliAltroBuffer", "copy constructor not implemented");
}

//_____________________________________________________________________________
AliAltroBuffer& AliAltroBuffer::operator = (const AliAltroBuffer& /*source*/)
{
//Assigment operator

  Fatal("operator =", "assignment operator not implemented");
  return *this;
}

//_____________________________________________________________________________
void AliAltroBuffer::Flush()
{
// Flushes the Buffer content 
  if (fFreeCellBuffer != 16) {
    Int_t temp = fFreeCellBuffer;
    for (Int_t i = 0; i < temp; i++){
      FillBuffer(0x2AA);
    }//end for
  }//end if
}

//_____________________________________________________________________________
void AliAltroBuffer::FillBuffer(Int_t val)
{
//Fills the Buffer with 16 ten bits words and write into a file 

  if ((val > 0x3FF) || (val < 0)) {
    Error("FillBuffer", "Value out of range (10 bits): %d", val);
    val = 0x3FF;
  }
  fFreeCellBuffer--;

  fBuffer[fCurrentCell] |= (val << fShift);
  fShift += 10;

  if (fShift > 32) {
    fCurrentCell++;
    fShift -= 32;
    fBuffer[fCurrentCell] |= (val >> (10 - fShift));
  }

  if (fShift == 32) {
    //Buffer is written into a file
    fFile->WriteBuffer((char*)fBuffer, sizeof(UInt_t)*5);
    //Buffer is empty
    for (Int_t j = 0; j < 5; j++) fBuffer[j] = 0;
    fShift = 0;
    fCurrentCell = 0;
    fFreeCellBuffer = 16;
  }
}


//_____________________________________________________________________________
void AliAltroBuffer::WriteDummyTrailer(Int_t wordsNumber, Int_t padNumber,
				       Int_t rowNumber, Int_t secNumber)
{
//Writes a trailer of 40 bits

   Int_t num = fFreeCellBuffer % 4;
   for(Int_t i = 0; i < num; i++) {
     FillBuffer(0x2AA);
   }//end for
   FillBuffer(wordsNumber);
   FillBuffer(padNumber);
   FillBuffer(rowNumber);
   FillBuffer(secNumber);
}

//_____________________________________________________________________________
void AliAltroBuffer::WriteTrailer(Int_t wordsNumber, Int_t padNumber,
				  Int_t rowNumber, Int_t secNumber)
{
//Writes a trailer of 40 bits

  if (!fMapping) {
    AliError("No ALTRO mapping information is loaded! Filling a dummy trailer!");
    return WriteDummyTrailer(wordsNumber,padNumber,
			     rowNumber,secNumber);
  }

  Short_t hwAddress = fMapping->GetHWAddress(rowNumber,padNumber,secNumber);
  if (hwAddress == -1)
    AliFatal(Form("No hardware (ALTRO) adress found for these pad-row (%d) and pad (%d) indeces !",rowNumber,padNumber));
  WriteTrailer(wordsNumber,hwAddress);
}

//_____________________________________________________________________________
void AliAltroBuffer::WriteTrailer(Int_t wordsNumber, Short_t hwAddress)
{
//Writes a trailer of 40 bits using
//a given hardware adress
  Int_t num = fFreeCellBuffer % 4;
  for(Int_t i = 0; i < num; i++) {
    FillBuffer(0x2AA);
  }//end for
  Int_t temp;
  temp = hwAddress & 0x3FF;
  FillBuffer(temp);

  temp = (wordsNumber << 6) & 0x3FF;
  temp |= (0xA << 2);
  temp |= ((hwAddress >> 10) & 0x3);
  FillBuffer(temp);

  temp = 0xA << 6;
  temp |= ((wordsNumber & 0x3FF) >> 4);
  FillBuffer(temp);

  temp = 0x2AA;
  FillBuffer(temp);
}

//_____________________________________________________________________________
void AliAltroBuffer::WriteChannel(Int_t padNumber, Int_t rowNumber, 
				  Int_t secNumber,
				  Int_t nTimeBins, const Int_t* adcValues,
				  Int_t threshold)
{
  //Write all ADC values and the trailer of a channel
  Int_t nWords = WriteBunch(nTimeBins,adcValues,threshold);
  // write the trailer
  if (nWords) WriteTrailer(nWords, padNumber, rowNumber, secNumber);
}

//_____________________________________________________________________________
void AliAltroBuffer::WriteChannel(Short_t hwAddress,
				  Int_t nTimeBins, const Int_t* adcValues,
				  Int_t threshold)
{
  //Write all ADC values and the trailer of a channel
  Int_t nWords = WriteBunch(nTimeBins,adcValues,threshold);
  // write the trailer
  if (nWords) WriteTrailer(nWords, hwAddress);
}

//_____________________________________________________________________________
Int_t AliAltroBuffer::WriteBunch(Int_t nTimeBins, const Int_t* adcValues,
				 Int_t threshold)
{
  //Write all ADC values
  //Return number of words written

  Int_t nWords = 0;
  Int_t timeBin = -1;
  Int_t bunchLength = 0;

  // loop over time bins
  for (Int_t iTime = 0; iTime < nTimeBins; iTime++) {
    if (adcValues[iTime] >= threshold) { // ADC value above threshold
      FillBuffer(adcValues[iTime]);
      nWords++;
      timeBin = iTime;
      bunchLength++;

    } else if (timeBin >= 0) {  // end of bunch
      FillBuffer(timeBin);
      FillBuffer(bunchLength + 2);
      nWords += 2;
      timeBin = -1;
      bunchLength = 0;
    }
  }

  if (timeBin >= 0) {  // end of bunch
    FillBuffer(timeBin);
    FillBuffer(bunchLength + 2);
    nWords += 2;
  }

  return nWords;
}

//_____________________________________________________________________________
void AliAltroBuffer::WriteDataHeader(Bool_t dummy, Bool_t compressed)
{
//Write a (dummy or real) DDL data header, 
//set the compression bit if compressed

  AliRawDataHeaderSim header;
  if (dummy) {
    //if size=0 it means that this data header is a dummy data header
    fDataHeaderPos = fFile->Tellp();
    fFile->WriteBuffer((char*)(&header), sizeof(header));
  } else {
    WriteRCUTrailer(0);
    UInt_t currentFilePos = fFile->Tellp();
    fFile->Seekp(fDataHeaderPos);
    header.fSize = currentFilePos-fDataHeaderPos;
    header.SetAttribute(0);  // valid data
    if (compressed) header.SetAttribute(1); 
    fFile->WriteBuffer((char*)(&header), sizeof(header));
    fFile->Seekp(currentFilePos);
  }
}

//_____________________________________________________________________________
void AliAltroBuffer::WriteRCUTrailer(Int_t rcuId)
{
  // Writes the RCU trailer
  // rcuId the is serial number of the corresponding
  // RCU. The basic format of the trailer can be
  // found in the RCU manual.
  // This method should be called at the end of
  // raw data writing.

  UInt_t currentFilePos = fFile->Tellp();
  UInt_t size = currentFilePos-fDataHeaderPos;
  size -= sizeof(AliRawDataHeader);
  
  if ((size % 5) != 0) {
    AliFatal(Form("The current raw data payload is not a mutiple of 5 (%d) ! Can not write the RCU trailer !",size));
    return;
  }

  // Now put the size in unit of number of 40bit words
  size /= 5;
  fFile->WriteBuffer((char *)(&size),sizeof(UInt_t));

  // Now several not yet full defined fields
  // In principle they are supposed to contain
  // information about the sampling frequency,
  // L1 phase, list of 'dead' FECs, etc.
  //  UInt_t buffer[n];
  //  fFile->WriteBuffer((char *)(buffer),sizeof(UInt_t)*n);
  
  //  Now the RCU identifier and size of the trailer
  //  FOr the moment the triler size is 2 32-bit words
  UInt_t buffer = 2;
  buffer |= ((rcuId & 0x3FF) << 22);
  fFile->WriteBuffer((char *)(&buffer),sizeof(UInt_t));

}
