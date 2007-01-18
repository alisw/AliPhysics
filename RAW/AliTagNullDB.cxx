// @(#)alimdc:$Name$:$Id$
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
// AliTagNullDB                                                         //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "AliTagNullDB.h"


ClassImp(AliTagNullDB)


//______________________________________________________________________________
AliTagNullDB::AliTagNullDB(AliRawEventTag *eventTag) :
   AliTagDB(eventTag, "/dev/null")
{
   // Create tag db writing to /dev/null.

}

//______________________________________________________________________________
const char *AliTagNullDB::GetFileName() const
{
   // Return /dev/null as filename.

   return "/dev/null";
}

//______________________________________________________________________________
void AliTagNullDB::Close()
{
   // Close null tag DB.

   if (!fTagDB) return;

   fTagDB->cd();

   // Write the tree.
   fTree->Write();

   // Close DB, this also deletes the fTree
   fTagDB->Close();

   delete fTagDB;
   fTagDB = 0;
}
