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
// AliTagDB                                                             //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <errno.h>

#include <TSystem.h>
#include <TTimeStamp.h>
#include <TBranch.h>

#include "AliESD.h"

#include "AliRawDB.h"
#include "AliRawEventTag.h"
#include "AliTagDB.h"
#include "AliRawEventHeaderBase.h"


ClassImp(AliTagDB)


//______________________________________________________________________________
AliTagDB::AliTagDB(AliRawEventTag *eventTag, const char* fileName) :
  fTagDB(NULL),
  fTree(NULL),
  fEventTag(eventTag),
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

static void BranchResetBit(TBranch *b) 
{
  // Reset MapObject on this branch and all the sub-branches

  b->ResetBit( kBranchObject | kBranchAny ); // Or in newer ROOT: b->ResetBit( kMapObject )
  TIter next( b->GetListOfBranches() );
  TBranch *sub = 0;
  while ( (sub = (TBranch*)next() ) ) {
    BranchResetBit( sub );
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
   fTree = new TTree("T", Form("ALICE raw-data tag tree (%s)", AliRawDB::GetAliRootTag()));
   fTree->SetAutoSave(100000000);  // autosave when 100 Mbyte written

   Int_t bufsize = 32000;
   Int_t split   = 1;
   const char *tagname = fEventTag->GetName();
   TBranch * b = fTree->Branch("TAG", tagname, &fEventTag, bufsize, split);
   BranchResetBit(b);

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

   static TString fname;
   const char *fs = fFS;

   // check that fs exists (crude check fails if fs is a file)
   gSystem->MakeDirectory(fs);

   // Get the run number
   Int_t runNumber = -1;
   if (fEventTag) {
     AliRawEventHeaderBase *header = fEventTag->GetHeader();
     if (header) runNumber = header->Get("RunNb");
   }

   TString hostname;
   hostname.Form("%s",gSystem->HostName());
   if ( hostname.First('.') > 0 ) hostname.Resize(hostname.First('.'));

   TTimeStamp ts;

   fname.Form("%s/Run%d.%s_%d_%d_%d.RAW.tag.root", fs, runNumber, hostname.Data(),
	      ts.GetDate(), ts.GetTime(), ts.GetNanoSec());

   return fname.Data();
}
