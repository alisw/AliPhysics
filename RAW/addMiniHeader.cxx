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
// this program adds a mini header to a binary file                          //
//                                                                           //
// type "addMiniHeader -?" to get a description of the arguments             //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <stdlib.h>
#include <string.h>
#include <sys/sendfile.h>
#include <TError.h>
#include <TSystem.h>
#include "AliMiniHeader.h"


//_____________________________________________________________________________
int main(int argc, char** argv)
{
  if ((argc < 2) || (argc > 4) || (strcmp(argv[1], "-?") == 0) || 
      (strcmp(argv[1], "-h") == 0) || (strcmp(argv[1], "--help") == 0)) {
    printf("\nUsage: addMiniHeader fileName [detectorID [DDLID]]\n\n");
    printf("         fileName   : name of the raw data file\n");
    printf("         detectorID : number of the detector (default 0xFF)\n");
    printf("         DDLID      : number of the DDL (default 0xFFFF)\n\n");
    return 0;
  }

  // prepare the mini header
  AliMiniHeader header;
  header.fMagicWord[2] = 0x12;
  header.fMagicWord[1] = 0x34;
  header.fMagicWord[0] = 0x56;
  header.fDetectorID = 0xFF;
  header.fDDLID = 0xFFFF;
  header.fCompressionFlag = 0;
  header.fVersion = 1;

  // open the raw data file
  FILE* file = fopen(argv[1], "rb");
  if (!file) {
    ::Fatal("addMiniHeader", "could not open file %s", argv[1]);
  }

  // set the correct raw data size in the mini header
  fseek(file, 0, SEEK_END);
  header.fSize = ftell(file);
  fseek(file, 0, SEEK_SET);

  // set the detector and DDL id if they are given
  if (argc > 2) {
    header.fDetectorID = atoi(argv[2]);
  }
  if (argc > 3) {
    header.fDDLID = atoi(argv[3]);
  }

  // open a temporary file and write the mini header to it
  char tmpFileName[] = "tmpXXXXXX";
  Int_t tmpID = mkstemp(tmpFileName);
  if (tmpID < 0) {
    ::Fatal("addMiniHeader", "could not open temporary file");
  }
  FILE* tmp = fdopen(tmpID, "w");
  fwrite(&header, sizeof(header), 1, tmp);

  // copy the raw data to the temporary file
  char buffer[256];
  while (Int_t size = fread(buffer, 1, 256, file)) {
    fwrite(buffer, size, 1, tmp);
  }
  fclose(file);
  fclose(tmp);

  // replace the original file by the temporary one
  gSystem->Rename(tmpFileName, argv[1]);

  return 0;
}
