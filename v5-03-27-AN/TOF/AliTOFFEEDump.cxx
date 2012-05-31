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

/*
  author: Roberto Preghenella (preghenella@bo.infn.it)
*/
///////////////////////////////////////////////////////////////
//                                                           //
//   This classes provide the object to store the full dump  //
//   of TOF FEE configuration database.                      //
//                                                           //
///////////////////////////////////////////////////////////////

#include "AliTOFFEEDump.h"
#include <string.h>
#include <iostream>
#include <fstream>
#include "TSystem.h"
#include "AliLog.h"

ClassImp(AliTOFFEEDump)

//_______________________________________________________________

AliTOFFEEDump::AliTOFFEEDump() :
  TObject(),
  fSize(0),
  fData(NULL)
{
  /* default constructor */
}

#if 0
//_______________________________________________________________

AliTOFFEEDump::AliTOFFEEDump(const AliTOFFEEDump &source) :
  TObject(source),
  fSize(source.fSize),
  fData(NULL)
{
  /* copy constructor */
  
  /* check size */
  if (fSize == 0) return;

  /* allocate and copy data */
  fData = new UChar_t[fSize];
  memcpy(fData, source.fData, fSize);
}

//_______________________________________________________________

AliTOFFEEDump &
AliTOFFEEDump::operator=(const AliTOFFEEDump &source)
{
  /* operator= */
  
  /* check source and destination size */
  if (source.fSize == 0 || fSize != source.fSize) return *this;

  /* copy data */
  memcpy(fData, source.fData, fSize);
  return *this;
}
#endif

//_______________________________________________________________

AliTOFFEEDump::~AliTOFFEEDump()
{
  /* default destructor */

  if (fData) delete [] fData;
}

//_______________________________________________________________

Bool_t
AliTOFFEEDump::operator!=(const AliTOFFEEDump &source)
{
  /* operator!= */
  
  /* check size */
  if (fSize != source.fSize) return kTRUE;

  /* check data */
  if (memcmp(fData, source.fData, fSize) != 0) return kTRUE;

  return kFALSE;
}

//_______________________________________________________________

Bool_t
AliTOFFEEDump::ReadFromFile(const Char_t *filename)
{
  /* read from file */

  /* open file */
  Char_t *expandedFileName = gSystem->ExpandPathName(filename);
  std::ifstream is;
  is.open(expandedFileName, std::ios::binary);
  if (!is.is_open()) {
    AliError(Form("error while opening TOF FEE dump file: %s", filename));
    return kFALSE;
  }
  AliInfo(Form("TOF FEE dump file opened: %s", filename));

  /* get file size */
  Int_t begin = is.tellg();
  is.seekg(0, std::ios::end); /* end */
  Int_t end = is.tellg();
  Int_t size = end - begin;
  is.seekg(0, std::ios::beg); /* rewind file */
  if (size <= 0) {
    AliError(Form("error while getting TOF FEE dump file size: %d", size));
    return kFALSE;
  }
  AliInfo(Form("got TOF FEE dump file size: %d", size));

  /* check previous allocation */
  if (fData) {
    AliWarning("data already allocated, old data will be overwritten");
    delete [] fData;
  }

  /* allocate and read data */
  fSize = size;
  fData = new UChar_t[fSize];
  is.read((Char_t *)fData, fSize);
  AliInfo(Form("TOF FEE dump file stored"));

  /* close file */
  is.close();

  return kTRUE;
}

//_______________________________________________________________

void
AliTOFFEEDump::DumpData() {
  /* dump data */

  printf("*** TOF FEE dump data ***\n");
  printf("data size = %d bytes\n", fSize);
  printf("*************************\n");
  Int_t nwords = fSize / 4;
  UInt_t *data = (UInt_t *)fData;
  for (Int_t iword = 0; iword < nwords; iword++) {
    if (iword != 0 && iword % 4 == 0) printf("\n");
    printf("%08x ", data[iword]);
  }
  printf("\n*************************\n");
  
}
