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

//
// Macro to generate manu list per buspatch for station 3, 4 and 5
// Christian Finck, Subatech
//

#if !defined(__CINT__) || defined(__MAKECINT__)

#include "AliMpSegmentation.h"
#include "AliMpVSegmentation.h"
#include "AliMpDEManager.h"
#include "AliMpBusPatch.h"
#include "AliMpIntPair.h"

#include <TArrayI.h>
#include <TList.h>

#endif

void testManuListSt345(Char_t* fileNamePre = "BusPatchToManu")
{
   
  Int_t planeOffset = 1024;
  Int_t nBusPatch = 0;
  Int_t nManu = 0;

  Char_t fileName[255];
  AliMpBusPatch* busPatchManager = new AliMpBusPatch();
  busPatchManager->ReadBusPatchFile();

  // loop over DDL station 345
  for (Int_t iDDL = 8; iDDL < 20; ++iDDL) {

    sprintf(fileName,"%s%d%s", fileNamePre, iDDL, ".dat");
 
    FILE* fp = fopen(fileName,"w");
    printf("DDL # %d\n", iDDL);
    fprintf(fp, "DDL # %d\n", iDDL);
    TArrayI deArray = busPatchManager->GetDeInDDL(iDDL);

    // list of DE in the given DDL
    for (Int_t detElemId = 0; detElemId < deArray.GetSize(); ++detElemId) {

      // list of bus patch for a given DE
      TArrayI* busArray = busPatchManager->GetBusfromDE(deArray[detElemId]);

      // list of manu per DE
      TList manuList;
      for ( Int_t cath = 0; cath <=1 ; ++cath ) {
	const AliMpVSegmentation* seg 
	  = AliMpSegmentation::Instance()->GetMpSegmentation(deArray[detElemId],cath);
        
	TArrayI manus;

	seg->GetAllElectronicCardIDs(manus);
          
	// filling class
	for ( Int_t im = 0; im < manus.GetSize(); ++im ) {

	  AliMpIntPair* manu = new AliMpIntPair((manus[im] & 0x3FF), cath, kTRUE);// remove the offfset
	  manuList.Add(manu);
	}        
      }
      manuList.Sort();
      //      printf("Number of manu %d\n", manuList.GetEntries());

      // writing into files/stdout
      Int_t posPrev = -1;
      Int_t manuIdPrev = 0;
      for (Int_t iEntry = 0; iEntry < manuList.GetEntries(); iEntry++) {

	AliMpIntPair* manuPtr = (AliMpIntPair*)manuList.At(iEntry);
	Int_t pos = manuPtr->GetFirst()/100;
	Int_t manuId;
	nManu++;

	if (manuPtr->GetSecond())
	  manuId = manuPtr->GetFirst() + planeOffset;
	else
	  manuId = manuPtr->GetFirst();

	if (pos != posPrev) {
	  if (posPrev != -1) {
	    printf("%d;\n", manuIdPrev);
	    fprintf(fp,"%d;\n", manuIdPrev);

	  }
	  printf("%d %d %d-", deArray[detElemId], 
		 AliMpBusPatch::GetLocalBusID(busArray->At(pos),iDDL), manuId);
	  fprintf(fp,"%d %d-", AliMpBusPatch::GetLocalBusID(busArray->At(pos),iDDL), manuId);
	  nBusPatch++;

	} else if ( manuId != manuIdPrev+1) {
	  printf("%d; ", manuIdPrev);
	  printf("%d-", manuId);
	  fprintf(fp,"%d; ", manuIdPrev);
	  fprintf(fp,"%d-", manuId);

	}
     
	posPrev = pos;
	manuIdPrev = manuId;
	if (iEntry ==  manuList.GetEntries()-1) {
	  printf("%d;\n", manuId);
	  fprintf(fp,"%d;\n", manuId);
	}
      }
    }
    fclose(fp);
  }
  fflush(stdout);
  printf("Number of buspatches %d and manus %d\n", nBusPatch, nManu);

  delete busPatchManager;
}

