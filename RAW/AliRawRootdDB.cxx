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
// AliRawRootdDB                                                        //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <TSystem.h>

#include "AliMDC.h"

#include "AliRawRootdDB.h"


ClassImp(AliRawRootdDB)


//______________________________________________________________________________
AliRawRootdDB::AliRawRootdDB(AliRawEvent *event, Double_t maxsize, Int_t compress)
   : AliRawDB(event, maxsize, compress, kFALSE)
{
   // Create a new raw DB that will be accessed via rootd daemon.

   if (!Create())
      MakeZombie();
   else
      fRawDB->UseCache(50, 0x200000);  //0x100000 = 1MB)
}

//______________________________________________________________________________
const char *AliRawRootdDB::GetFileName() const
{
   // Return filename based on hostname and date and time. This will make
   // each file unique. Also the directory will be made unique for each
   // day by adding the date to the fs. Assumes there is always enough
   // space on the device.

   static TString fname;

   TString fs = AliMDC::RootdFS();
   TDatime dt;

#if 0
   // make a new subdirectory for each day
   fs += "/adc-";
   fs += dt.GetDate();

   Long_t id, size, flags, time;
   if (gSystem->GetPathInfo(fs, &id, &size, &flags, &time) == 1) {
      // directory does not exist, create it
      if (gSystem->mkdir(fs, kTRUE) == -1) {
         Error("GetFileName", "cannot create dir %s, using %s", fs.Data(),
               AliMDC::RootdFS());
         fs = AliMDC::RootdFS();
      }
   }
   // FIXME: should check if fs is a directory
#endif

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
void AliRawRootdDB::Close()
{
   // Close raw rootd DB.

   if (!fRawDB) return;

   fRawDB->cd();

   // Write the tree.
   fTree->Write();

   // Close DB, this also deletes the fTree
   fRawDB->Close();

#if 0
   // can use services of TFTP
   if (AliMDC::DeleteFiles())
      gSystem->Exec(Form("rfrm %s", fRawDB->GetName()));
#endif

   delete fRawDB;
   fRawDB = 0;
}
