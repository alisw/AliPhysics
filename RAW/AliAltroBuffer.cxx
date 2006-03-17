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
// Author: D.Favretto

#include "AliAltroBuffer.h"
#include "AliAltroMapping.h"
#include "AliRawDataHeader.h"
#include "AliLog.h"
#include <Riostream.h>
#include <stdlib.h>


ClassImp(AliAltroBuffer)

//_____________________________________________________________________________
AliAltroBuffer::AliAltroBuffer(const char* fileName, Int_t flag, const AliAltroMapping *mapping):
  fShift(0),
  fCurrentCell(0),
  fFreeCellBuffer(0),
  fFlag(flag),
  fVerbose(0),
  fFile(NULL),
  fMaskBackward(0xFF),
  fFilePosition(0),
  fFileEnd(0),
  fDataHeaderPos(0),
  fEndingFillWords(0),
  fMapping(mapping)
{
//if flag = 1 the actual object is used in the write mode
//if flag = 0 the actual object is used in the read mode

  //the buffer is cleaned 
  for (Int_t i = 0; i < 5; i++) fBuffer[i] = 0;

  if (flag) {
    fFreeCellBuffer = 16;
    fShift = 32; 
    //open the output file
#ifndef __DECCXX
    fFile = new fstream(fileName, ios::binary|ios::out);
#else
    fFile = new fstream(fileName, ios::out);
#endif
  } else {
    //open the input file
#ifndef __DECCXX
    fFile = new fstream(fileName, ios::binary|ios::in);
#else
    fFile = new fstream(fileName, ios::in);
#endif
    if (!fFile) {
      Error("AliAltroBuffer", "File doesn't exist: %s", fileName);
      return;
    }
    fShift = 0;
    //To get the file dimension (position of the last element in term of bytes)
    fFile->seekg(0, ios::end);
    fFilePosition = fFile->tellg();
    fFileEnd = fFilePosition;
    fFile->seekg(0);
  }

}

//_____________________________________________________________________________
AliAltroBuffer::~AliAltroBuffer()
{
// destructor

  if (fFlag) {
    //Flush out the Buffer content at the end only if Buffer wasn't completely filled
    Flush();
    if (fVerbose) Info("~AliAltroBuffer", "File Created");
  }//end if
  fFile->close();
  delete fFile;

}

//_____________________________________________________________________________
AliAltroBuffer::AliAltroBuffer(const AliAltroBuffer& source):
  TObject(source),
  fShift(source.fShift),
  fCurrentCell(source.fCurrentCell),
  fFreeCellBuffer(source.fFreeCellBuffer),
  fFlag(source.fFlag),
  fVerbose(source.fVerbose),
  fFile(NULL),
  fMaskBackward(source.fMaskBackward),
  fFilePosition(source.fFilePosition),
  fFileEnd(source.fFileEnd),
  fDataHeaderPos(source.fDataHeaderPos),
  fEndingFillWords(source.fEndingFillWords),
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
Int_t AliAltroBuffer::GetNext()
{
//It reads a 10 bits word in forward dicection from the Buffer.
//A new Buffer is read from the file only when Buffer is empty.
//If there aren't elements anymore -1 is returned otherwise 
//the next element is returned

  UInt_t mask = 0xFFC00000;
  UInt_t temp;
  UInt_t value;
  if (!fShift) {
    if (fFile->tellg() >= (Int_t)fFileEnd) return -1;
    if (fFile->read((char*)fBuffer, sizeof(UInt_t)*5)) {
      fCurrentCell = 0;
      fShift = 22;
      value = fBuffer[fCurrentCell] & mask;
      value = value >> 22;
      fBuffer[fCurrentCell] = fBuffer[fCurrentCell] << 10;
      return value;      
    } else {
      return -1;
    }
  } else {
    if (fShift >= 10) {
      value = fBuffer[fCurrentCell] & mask;
      value = value >> 22;
      fShift -= 10;
      fBuffer[fCurrentCell] = fBuffer[fCurrentCell] << 10;
    } else {
      value = fBuffer[fCurrentCell] & mask;
      fCurrentCell++;
      temp = fBuffer[fCurrentCell];
      temp = temp >> fShift;
      temp = temp & mask;
      value = value | temp;
      value = value >> 22;
      fBuffer[fCurrentCell] = fBuffer[fCurrentCell] << (10-fShift);
      fShift += 22;
    }
    return value;
  }//end else
}

//_____________________________________________________________________________
Int_t AliAltroBuffer::GetNextBackWord()
{
//It reads a 10 bits word in backward dicection from the Buffer.
//A new Buffer is read from the file only when Buffer is empty.
//If there aren't elements anymore -1 is returned otherwise 
//the next element is returned

  UInt_t mask = 0x3FF;
  UInt_t temp;
  UInt_t value;
  if (!fShift) {
    if (fFilePosition > fDataHeaderPos){
      fFilePosition -= sizeof(UInt_t)*5;
      fFile->seekg(fFilePosition);
      fFile->read((char*)fBuffer, sizeof(UInt_t)*5);
      
      fCurrentCell = 4;
      fShift = 22;
      fMaskBackward = 0xFF;
      value = fBuffer[fCurrentCell] & mask;
      fBuffer[fCurrentCell] = fBuffer[fCurrentCell] >> 10;
      return value;
    } else {
      fFile->seekg(fDataHeaderPos);
      return -1;
    }
  } else {
    if (fShift >= 10) {
      value = fBuffer[fCurrentCell] & mask;
      fShift -= 10;
      fBuffer[fCurrentCell] = fBuffer[fCurrentCell] >> 10;
    } else {
      value = fBuffer[fCurrentCell];
      fCurrentCell--;
      temp = fBuffer[fCurrentCell] & mask;
      temp = temp & fMaskBackward;
      fMaskBackward = fMaskBackward >> 2;
      temp = temp << fShift;
      value = value | temp;
      fBuffer[fCurrentCell] = fBuffer[fCurrentCell] >> (10-fShift);
      fShift = 22 + fShift;
    }
    return value;
  }//end else
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
  if (fShift < 10) {
    Int_t temp = val;
    val = val >> (10-fShift);
    fBuffer[fCurrentCell] |= val;
    fCurrentCell++;
    fShift += 32;
    val = temp;
  }
  fShift -= 10;
  val = val << fShift;
  fBuffer[fCurrentCell] |= val;
  if (!fShift) {
    //Buffer is written into a file
    fFile->write((char*)fBuffer, sizeof(UInt_t)*5);
    //Buffer is empty
    for (Int_t j = 0; j < 5; j++) fBuffer[j] = 0;
    fShift = 32;
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

  Short_t hwAdress = fMapping->GetHWAdress(rowNumber,padNumber,secNumber);
  if (hwAdress == -1)
    AliFatal(Form("No hardware (ALTRO) adress found for these pad-row (%d) and pad (%d) indeces !",rowNumber,padNumber));
  WriteTrailer(wordsNumber,hwAdress);
}

//_____________________________________________________________________________
void AliAltroBuffer::WriteTrailer(Int_t wordsNumber, Short_t hwAdress)
{
//Writes a trailer of 40 bits using
//a given hardware adress
  Int_t num = fFreeCellBuffer % 4;
  for(Int_t i = 0; i < num; i++) {
    FillBuffer(0x2AA);
  }//end for
  Int_t temp;
  temp = 0x2AA;
  FillBuffer(temp);
  temp = 0xA << 6;
  temp |= ((wordsNumber & 0x3FF) >> 4);
  FillBuffer(temp);
  temp = (wordsNumber << 6) & 0x3FF;
  temp |= (0xA << 2);

  temp |= (hwAdress >> 10) & 0x3;
  FillBuffer(temp);
  temp = hwAdress & 0x3FF;
  FillBuffer(temp);
}

//_____________________________________________________________________________
Bool_t AliAltroBuffer::ReadDummyTrailer(Int_t& wordsNumber, Int_t& padNumber,
					Int_t& rowNumber, Int_t& secNumber)
{
//Read a dummy trailer of 40 bits in the forward reading mode

  wordsNumber = GetNext();
  if (wordsNumber == -1) return kFALSE;
  padNumber = GetNext();
  if (padNumber == -1) return kFALSE;
  rowNumber = GetNext();
  if (rowNumber == -1) return kFALSE;
  secNumber = GetNext();
  if (secNumber == -1) return kFALSE;
  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliAltroBuffer::ReadTrailer(Int_t& wordsNumber, Int_t& padNumber,
				   Int_t& rowNumber, Int_t& secNumber)
{
//Read a trailer of 40 bits in the forward reading mode
  if (!fMapping) {
    AliError("No ALTRO mapping information is loaded! Reading a dummy trailer!");
    return ReadDummyTrailer(wordsNumber,padNumber,
			    rowNumber,secNumber);
  }

  Short_t hwAdress;
  if (!ReadTrailer(wordsNumber,hwAdress)) return kFALSE;
  rowNumber = fMapping->GetPadRow(hwAdress);
  padNumber = fMapping->GetPad(hwAdress);
  secNumber = fMapping->GetSector(hwAdress);

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliAltroBuffer::ReadTrailer(Int_t& wordsNumber, Short_t& hwAdress)
{
//Read a trailer of 40 bits in the forward reading mode

  Int_t temp = GetNext();
  if (temp != 0x2AA)
    AliFatal(Form("Incorrect trailer found ! Expecting 0x2AA but found %x !",temp));

  temp = GetNext();
  if ((temp >> 6) != 0xA)
    AliFatal(Form("Incorrect trailer found ! Expecting 0xA but found %x !",temp >> 6));
  wordsNumber = (temp << 4) & 0x3FF;

  temp = GetNext();
  wordsNumber |= (temp >> 6);
  if (((temp >> 2) & 0xF) != 0xA)
    AliFatal(Form("Incorrect trailer found ! Expecting second 0xA but found %x !",temp >> 6));
  hwAdress = (temp & 0x3) << 10;

  temp = GetNext();
  hwAdress |= temp;

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliAltroBuffer::ReadDummyTrailerBackward(Int_t& wordsNumber, Int_t& padNumber,
						Int_t& rowNumber, Int_t& secNumber)
{
//Read a trailer of 40 bits in the backward reading mode

  Int_t temp;
  fEndingFillWords = 0;
  do {
    temp = GetNextBackWord();
    fEndingFillWords++;
    if (temp == -1) return kFALSE;
  } while (temp == 0x2AA);  
  fEndingFillWords--;
  secNumber = temp;
  rowNumber = GetNextBackWord();
  if (rowNumber == -1) return kFALSE;
  padNumber = GetNextBackWord();
  if (padNumber == -1) return kFALSE;
  wordsNumber = GetNextBackWord();
  if (wordsNumber == -1) return kFALSE;
  return kTRUE;
} 

//_____________________________________________________________________________
Bool_t AliAltroBuffer::ReadTrailerBackward(Int_t& wordsNumber, Int_t& padNumber,
					   Int_t& rowNumber, Int_t& secNumber)
{
//Read a trailer of 40 bits in the backward reading mode
  if (!fMapping) {
    AliError("No ALTRO mapping information is loaded! Reading a dummy trailer!");
    return ReadDummyTrailerBackward(wordsNumber,padNumber,
				    rowNumber,secNumber);
  }

  Short_t hwAdress;
  if (!ReadTrailerBackward(wordsNumber,hwAdress)) return kFALSE;
  rowNumber = fMapping->GetPadRow(hwAdress);
  padNumber = fMapping->GetPad(hwAdress);
  secNumber = fMapping->GetSector(hwAdress);

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliAltroBuffer::ReadTrailerBackward(Int_t& wordsNumber, Short_t& hwAdress)
{
//Read a trailer of 40 bits in the backward reading mode

  Int_t temp;
  fEndingFillWords = 0;
  do {
    temp = GetNextBackWord();
    fEndingFillWords++;
    if (temp == -1) return kFALSE;
  } while (temp == 0x2AA);  
  fEndingFillWords--;

  hwAdress = temp;

  temp = GetNextBackWord();
  hwAdress |= (temp & 0x3) << 10;
  if (((temp >> 2) & 0xF) != 0xA)
    AliFatal(Form("Incorrect trailer found ! Expecting second 0xA but found %x !",temp >> 6));
  wordsNumber = (temp >> 6);

  temp = GetNextBackWord();
  wordsNumber |= (temp << 4) & 0x3FF;
  if ((temp >> 6) != 0xA)
    AliFatal(Form("Incorrect trailer found ! Expecting 0xA but found %x !",temp >> 6));
  
  temp = GetNextBackWord();
  if (temp != 0x2AA)
    AliFatal(Form("Incorrect trailer found ! Expecting 0x2AA but found %x !",temp));

  return kTRUE;
} 

//_____________________________________________________________________________
void AliAltroBuffer::WriteChannel(Int_t padNumber, Int_t rowNumber, 
				  Int_t secNumber,
				  Int_t nTimeBins, const Int_t* adcValues,
				  Int_t threshold)
{
//Write all ADC values and the trailer of a channel

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

  // write the trailer
  WriteTrailer(nWords, padNumber, rowNumber, secNumber);
}

//_____________________________________________________________________________
void AliAltroBuffer::ReadChannel(Int_t padNumber, Int_t rowNumber, 
				 Int_t secNumber,
				 Int_t& nTimeBins, Int_t* adcValues)
{
//Read all ADC values and the trailer of a channel

  Int_t wordsNumber;
  if (!ReadTrailer(wordsNumber,padNumber,
		   rowNumber,secNumber)) return;

  if (wordsNumber < 0) return;
  // Number of fill words 
  Int_t nFillWords;
  if ((wordsNumber % 4) == 0)
    nFillWords = 0;
  else
    nFillWords = 4 - wordsNumber % 4;
  // Read the fill words 
  for (Int_t i = 0; i < nFillWords; i++) {
    Int_t temp = GetNext();
    if (temp != 0x2AA) 
      AliFatal(Form("Invalid fill word, expected 0x2AA, but got %X", temp));
  }

  // Decoding
  Int_t lastWord =  wordsNumber;
  nTimeBins = -1;
  while (lastWord > 0) { 
    Int_t l =  GetNext(); 
    if (l < 0) AliFatal(Form("Bad bunch length (%d) !", l));
    Int_t t =  GetNext(); 
    if (t < 0) AliFatal(Form("Bad bunch time (%d) !", t));
    lastWord -= 2;
    if (nTimeBins == -1) nTimeBins = t + 1;
    for (Int_t i = 2; i < l; i++) {
      Int_t amp = GetNext();
      if (amp < 0) AliFatal(Form("Bad adc value (%X) !", amp));
      adcValues[t - (i-2)] = amp;
      lastWord--;
    }
  }

} 

//_____________________________________________________________________________
void AliAltroBuffer::WriteDataHeader(Bool_t dummy, Bool_t compressed)
{
//Write a (dummy or real) DDL data header, 
//set the compression bit if compressed

  AliRawDataHeader header;
  if (dummy) {
    //if size=0 it means that this data header is a dummy data header
    fDataHeaderPos = fFile->tellp();
    fFile->write((char*)(&header), sizeof(header));
  } else {
    UInt_t currentFilePos = fFile->tellp();
    fFile->seekp(fDataHeaderPos);
    header.fSize = currentFilePos-fDataHeaderPos;
    header.SetAttribute(0);  // valid data
    if (compressed) header.SetAttribute(1); 
    fFile->write((char*)(&header), sizeof(header));
    fFile->seekp(currentFilePos);
  }
}

//_____________________________________________________________________________
Bool_t AliAltroBuffer::ReadDataHeader()
{
//Read the DDL data header at the beginning of the file, 
//returns true in case of valid data

  AliRawDataHeader header;
  UInt_t currentPos = fFile->tellp();
  fFile->seekp(0);
  if (!fFile->read((char*)(&header), sizeof(header))) return kFALSE;
  fDataHeaderPos = fFile->tellp();
  fFile->seekp(currentPos);
  return header.TestAttribute(0);
}

