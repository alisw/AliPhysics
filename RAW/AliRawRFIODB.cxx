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
// AliRawRFIODB                                                         //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <TSystem.h>
#include <TUrl.h>

#include "AliRawRFIODB.h"


ClassImp(AliRawRFIODB)


//______________________________________________________________________________
AliRawRFIODB::AliRawRFIODB(AliRawEvent *event,
			   AliESD *esd,
			   Int_t compress,
			   const char* fileName)
   : AliRawDB(event, esd, compress, fileName)
{
   // Create a new raw DB that will be accessed via RFIO.

   static int init = 0;
   // Set STAGE_POOL environment variable to current host
   if (!init) {
      // THESE ENVIRONMENT SYMBOLS ARE NOW DEFINED BY THE ALICE DATE SETUP
      // THEREFORE WE SHALL NOT USE ANY HARDCODED VALUES BUT RATHER USE
      // WHATEVER HAS BEEN SET IN THE DATE SITE
      //gSystem->Setenv("STAGE_POOL", "lcg00");
      //gSystem->Setenv("STAGE_HOST", "stage013");

      // however for sanity we check if they are really set
      if (!gSystem->Getenv("STAGE_POOL"))
         Error("AliRawRFIODB", "STAGE_POOL not set");
      if (!gSystem->Getenv("STAGE_HOST"))
         Error("AliRawRFIODB", "STAGE_HOST not set");
      init = 1;
   }

   if (fRawDB) fRawDB->UseCache(50, 0x200000);  //0x100000 = 1MB)
}

//______________________________________________________________________________
const char *AliRawRFIODB::GetFileName() const
{
   // Return filename based on hostname and date and time. This will make
   // each file unique. Also the directory will be made unique for each
   // day by adding the date to the fs. Assumes there is always enough
   // space on the device.

   static TString fname;

   TString fs = fFS1;
   TDatime dt;

   // make a new subdirectory for each day
   fs += "/adc-";
   fs += dt.GetDate();

   Long_t id, size, flags, time;
   if (gSystem->GetPathInfo(fs, &id, &size, &flags, &time) == 1) {
      // directory does not exist, create it
      if (gSystem->mkdir(fs, kTRUE) == -1) {
         Error("GetFileName", "cannot create dir %s, using %s", fs.Data(),
               fFS1.Data());
         fs = fFS1;
      }
   }
   // FIXME: should check if fs is a directory

   TString hostname = gSystem->HostName();
   Int_t pos;
   if ((pos = hostname.Index(".")) != kNPOS)
      hostname.Remove(pos);

   fname = fs + "/" + hostname + "_";
   fname += dt.GetDate();
   fname += "_";
   fname += dt.GetTime();
   fname += ".root";

   return fname;
}

//______________________________________________________________________________
void AliRawRFIODB::Close()
{
   // Close raw RFIO DB.

   if (!fRawDB) return;

   fRawDB->cd();

   // Write the tree.
   fTree->Write();
   if (fESDTree) fESDTree->Write();

   // Close DB, this also deletes the fTree
   fRawDB->Close();

   if (fDeleteFiles) {
      TUrl u(fRawDB->GetName());
      gSystem->Exec(Form("rfrm %s", u.GetFile()));
   }

   delete fRawDB;
   fRawDB = 0;
}
