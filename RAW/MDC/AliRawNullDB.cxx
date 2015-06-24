// @(#) $Id$
// Author: Fons Rademakers  26/11/99

/**************************************************************************
 * Copyright(c) 1998-2003, ALICE Experiment at CERN, All rights reserved. *
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

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// AliRawNullDB                                                         //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "AliRawNullDB.h"


ClassImp(AliRawNullDB)


//______________________________________________________________________________
AliRawNullDB::AliRawNullDB(AliRawEventV2 *event,
			   AliESDEvent *esd,
			   Int_t compress,
			   const char* fileName,Int_t basketsize, Long64_t autoflush)
   : AliRawDB(event, esd, compress, fileName, basketsize, autoflush)
{
   // Create a new raw DB that will wrtie to /dev/null.

}

//______________________________________________________________________________
const char *AliRawNullDB::GetFileName() const
{
   // Return /dev/null as filename.

   return "/dev/null";
}

//______________________________________________________________________________
Long64_t AliRawNullDB::Close()
{
   // Close raw RFIO DB.

   if (!fRawDB) return 0;

   if (!fRawDB->IsOpen()) return 0;

   fRawDB->cd();

   // Write the tree.
   Bool_t error = kFALSE;
   if (fTree)
     if (fTree->Write() == 0)
       error = kTRUE;
   if (fESDTree)
     if (fESDTree->Write() == 0)
       error = kTRUE;

   // Close DB, this also deletes the fTree
   fRawDB->Close();

   fTree = NULL;

   Long64_t filesize = fRawDB->GetEND();

   delete fRawDB;
   fRawDB = 0;
   if(!error)
     return filesize;
   else
     return -1;
}
