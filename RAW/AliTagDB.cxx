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

#include <errno.h>

#include <TSystem.h>

#include "AliESD.h"

#include "AliRawDB.h"
#include "AliRawEventHeaderBase.h"
#include "AliTagDB.h"


ClassImp(AliTagDB)


//______________________________________________________________________________
AliTagDB::AliTagDB(AliRawEventHeaderBase *header, const char* fileName) :
  fTagDB(NULL),
  fTree(NULL),
  fHeader(header),
  fMaxSize(-1),
  fFS(""),
  fDeleteFiles(kFALSE)
{
   // Create tag DB.

   if (fileName) {
      if (!Create(fileName))
         MakeZombie();
   }
}

//______________________________________________________________________________
Bool_t AliTagDB::Create(const char* fileName)
{
   // Create a new tag DB.

   const char *name = fileName;
   if (!name) name = GetFileName();
   fTagDB = new TFile(name, "RECREATE",
                      Form("ALICE tag DB (%s)", AliRawDB::GetAliRootTag()), 1);
   if (fTagDB->IsZombie()) {
      Error("Create", "error opening tag DB");
      fTagDB = 0;
      return kFALSE;
   }
   // Put wide read-write permissions
   if(gSystem->Chmod(name,438)) {
     Error("Create", "can't set permissions for tag DB file");
     fTagDB = 0;
     return kFALSE;
   }

   // Create ROOT Tree object container
   fTree = new TTree("TAG", Form("ALICE header data tree (%s)", AliRawDB::GetAliRootTag()));
   fTree->SetAutoSave(100000000);  // autosave when 100 Mbyte written

   Int_t bufsize = 32000;
   Int_t split   = 1;
   const char *headername = fHeader->GetName();
   fTree->Branch("header", headername, &fHeader, bufsize, split);

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

   if (fDeleteFiles)
      gSystem->Unlink(fTagDB->GetName());

   delete fTagDB;
   fTagDB = 0;
}

//______________________________________________________________________________
Bool_t AliTagDB::NextFile(const char* fileName)
{
   // Close te current file and open a new one.
   // Returns kFALSE in case opening failed.

   Close();

   if (!Create(fileName)) return kFALSE;
   return kTRUE;
}

//______________________________________________________________________________
void AliTagDB::SetFS(const char* fs)
{
// set the file system location

  fFS = fs;
  if (fs) {
    gSystem->ResetErrno();
    gSystem->MakeDirectory(fs);
    if (gSystem->GetErrno() && gSystem->GetErrno() != EEXIST) {
      SysError("SetFS", "mkdir %s", fs);
    }
  }
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
   const char *fs = fFS;

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
