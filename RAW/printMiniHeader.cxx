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
//                                                                           //
// this program prints the mini headers in a raw data file                   //
//                                                                           //
// type "printMiniHeader -?" to get a description of the arguments           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <stdlib.h>
#include <string.h>
#include <TError.h>
#include <TSystem.h>
#include "AliMiniHeader.h"


//_____________________________________________________________________________
int main(int argc, char** argv)
{
  if ((argc < 2) || (argc > 2) || (strcmp(argv[1], "-?") == 0) || 
      (strcmp(argv[1], "-h") == 0) || (strcmp(argv[1], "--help") == 0)) {
    printf("\nUsage: printMiniHeader fileName\n\n");
    printf("         fileName   : name of the raw data file\n\n");
    return 0;
  }

  // open the raw data file
  FILE* file = fopen(argv[1], "rb");
  if (!file) {
    ::Fatal("readMiniHeader", "could not open file %s", argv[1]);
  }
  fseek(file, 0, SEEK_END);
  Int_t size = ftell(file);
  fseek(file, 0, SEEK_SET);

  // read the mini headers
  AliMiniHeader header;
  while (!feof(file)) {
    if (!fread(&header, sizeof(header), 1, file)) {
      if (!feof(file)) {
	::Error("readMiniHeader", "could not read mini header at position %d", 
		ftell(file));
      }
      break;
    }

    // check the magic word
    if ((header.fMagicWord[2] != 0x12) ||
	(header.fMagicWord[1] != 0x34) ||
	(header.fMagicWord[0] != 0x56)) {
      ::Error("readMiniHeader", "the magic word is wrong %2x%2x%2x", 
	      header.fMagicWord[2], header.fMagicWord[1], 
	      header.fMagicWord[0]);
      break;
    }

    // check the size
    Int_t bytesLeft = size-ftell(file);
    if ((Int_t)header.fSize > bytesLeft) {
      ::Error("readMiniHeader", 
	      "the raw data size exceeds the file size by %d bytes",
	      header.fSize - bytesLeft);
    }

    // print the header information
    printf("size        : %d\n", header.fSize);
    printf("detector ID : %d\n", header.fDetectorID);
    printf("DDL ID      : %d\n", header.fDDLID);
    printf("compressed  : %d\n", header.fCompressionFlag);
    printf("version     : %d\n", header.fVersion);
    printf("\n");

    // go to next mini header
    if ((Int_t)header.fSize > bytesLeft) break;
    fseek(file, header.fSize, SEEK_CUR);
  }

  fclose(file);

  return 0;
}
