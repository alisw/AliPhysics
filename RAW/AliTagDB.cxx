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
// AliTagDB                                                             //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <TSystem.h>

#include "AliMDC.h"

#include "AliTagDB.h"


ClassImp(AliTagDB)


//______________________________________________________________________________
AliTagDB::AliTagDB(AliRawEventHeader *header, Double_t maxsize, Bool_t create)
{
   // Create tag DB.

   fHeader   = header;
   fMaxSize  = maxsize;

   if (create) {
      if (!Create())
         MakeZombie();
   }
}

//______________________________________________________________________________
AliTagDB::AliTagDB(const AliTagDB& tagDB): TObject(tagDB)
{
// copy constructor

  Fatal("AliTagDB", "copy constructor not implemented");
}

//______________________________________________________________________________
AliTagDB& AliTagDB::operator = (const AliTagDB& /*tagDB*/)
{
// assignment operator

  Fatal("operator =", "assignment operator not implemented");
  return *this;
}

//______________________________________________________________________________
Bool_t AliTagDB::Create()
{
   // Create a new tag DB.

   fTagDB = new TFile(GetFileName(), "RECREATE",
                      Form("ALICE MDC%d tag DB", AliMDC::kMDC), 1);
   if (fTagDB->IsZombie()) {
      Error("Create", "error opening tag DB");
      fTagDB = 0;
      return kFALSE;
   }

   // Create ROOT Tree object container
   fTree = new TTree("TAG", Form("ALICE MDC%d header data tree", AliMDC::kMDC));
   fTree->SetAutoSave(100000000);  // autosave when 100 Mbyte written

   Int_t bufsize = 32000;
   Int_t split   = 1;
   fTree->Branch("header", "AliRawEventHeader", &fHeader, bufsize, split);

   return kTRUE;
}

//______________________________________________________________________________
void AliTagDB::Close()
{
   // Close tag DB.

   if (!fTagDB) return;

   fTagDB->cd();

   // Write the tree.
   fTree->Write();

   // Close DB, this also deletes the fTree
   fTagDB->Close();

   if (AliMDC::DeleteFiles())
      gSystem->Unlink(fTagDB->GetName());

   delete fTagDB;
   fTagDB = 0;
}

//______________________________________________________________________________
Bool_t AliTagDB::NextFile()
{
   // Close te current file and open a new one.
   // Returns kFALSE in case opening failed.

   Close();

   if (!Create()) return kFALSE;
   return kTRUE;
}

//______________________________________________________________________________
Float_t AliTagDB::GetCompressionFactor() const
{
   // Return compression factor.

   if (fTree->GetZipBytes() == 0.)
      return 1.0;
   else
      return fTree->GetTotBytes()/fTree->GetZipBytes();
}

//______________________________________________________________________________
const char *AliTagDB::GetFileName() const
{
   // Return filename based on hostname and date and time. This will make
   // each file unique. The tags will be stored in the /data1/tags directory.

   static char fname[64];
   const char *fs = AliMDC::TagDBFS();

   // check that fs exists (crude check fails if fs is a file)
   gSystem->MakeDirectory(fs);

   char hostname[64];

   strcpy(hostname, gSystem->HostName());

   char *s;
   if ((s = strchr(hostname, '.')))
      *s = 0;

   TDatime dt;

   sprintf(fname, "%s/%s_%d_%d.root", fs, hostname, dt.GetDate(), dt.GetTime());

   return fname;
}
