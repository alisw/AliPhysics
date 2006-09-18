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
// AliRunDB                                                             //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <TSystem.h>
#include <TFile.h>
#include <TString.h>
#include <TDatime.h>
#include <TSQLServer.h>
#include <TSQLResult.h>
#include <TGrid.h>

#include "AliStats.h"
#include "AliRawDB.h"

#include "AliRunDB.h"


ClassImp(AliRunDB)


//______________________________________________________________________________
AliRunDB::AliRunDB(const char* localFS, Bool_t rdbms,
		   const char* alienHost, const char* alienDir) :
  fRunDB(NULL),
  fRDBMS(rdbms),
  fAlienHost(alienHost),
  fAlienDir(alienDir)
{
   // Open run database, and get or create tree.

   if (!localFS) return;

   // Get hostname
   char hostname[64], filename[64];

   // check that fs exists (crude check fails if fs is a file)
   gSystem->MakeDirectory(localFS);

   // Put wide read-write permissions
   if(gSystem->Chmod(localFS,1023)) {
     Error("AliRunDB","can't set permissions for run DB directory");
     return;
   }

   strcpy(hostname, gSystem->HostName());

   char *s;
   if ((s = strchr(hostname, '.')))
      *s = 0;

   sprintf(filename, "%s/%s_rundb.root", localFS, hostname);

   if (!gSystem->AccessPathName(filename, kFileExists))
      fRunDB = new TFile(filename, "UPDATE");
   else
     fRunDB = new TFile(filename, "CREATE", Form("ALICE Run DB (%s)", AliRawDB::GetAliRootTag()));

   // Put wide read-write permissions
   if(gSystem->Chmod(filename,438)) {
     Error("AliRunDB","can't set permissions for run DB file");
     return;
   }
}

//______________________________________________________________________________
void AliRunDB::Update(AliStats *stats)
{
  UpdateLocal(stats);
  UpdateRDBMS(stats);
  UpdateAliEn(stats);
}

//______________________________________________________________________________
void AliRunDB::UpdateLocal(AliStats *stats)
{
   // Add stats object to database.

   if (!stats || !fRunDB) return;

   TDirectory *ds = gDirectory;
   fRunDB->cd();

   char sname[64];
   char *s = (char*)strrchr(stats->GetFileName(), '/');
   if (s) {
      s++;
      strcpy(sname, s);
   } else
      strcpy(sname, stats->GetFileName());
   s = strchr(sname, '.');
   if (s) *s = 0;

   stats->Write(sname);

   ds->cd();
}

//______________________________________________________________________________
void AliRunDB::UpdateRDBMS(AliStats *stats)
{
   // Add stats object to central MySQL DB.

   if (!stats || !fRDBMS) return;

   char sql[4096];
   char bt[25], et[25];

   strcpy(bt, stats->GetBeginTime().AsSQLString());
   strcpy(et, stats->GetEndTime().AsSQLString());

   sprintf(sql, "INSERT INTO mdccatalog VALUES (0, '%s', %d, "
           "%d, %d, %d, %d, %d, %d, %.2f, '%s', '%s', '%s')",
           stats->GetFileName(), (int)stats->GetFileSize(), stats->GetEvents(),
           stats->GetFirstRun(), stats->GetFirstEvent(), stats->GetLastRun(),
           stats->GetLastEvent(), stats->GetCompressionMode(),
           stats->GetCompressionFactor(), stats->GetFilterState() ? "on" : "off",
           bt, et);

   // open connection to MySQL server on pcsalo
//    TSQLServer *db = TSQLServer::Connect("mysql://pcsalo.cern.ch/mdc", "alice", "amdc");

//    if (!db || db->IsZombie()) {
//       Error("UpdateRDBMS", "failed to connect to MySQL server on pcsalo");
//       printf("%s\n", sql);
//       delete db;
//       return;
//    }

//    TSQLResult *res = db->Query(sql);

//    if (!res) {
//       Error("UpdateRDBMS", "insert into mdccatalog failed");
//       printf("%s\n", sql);
//    }

//    delete res;
//    delete db;
}

//______________________________________________________________________________
void AliRunDB::UpdateAliEn(AliStats *stats)
{
   // Record file in AliEn catalog.

   if (!stats || fAlienHost.IsNull()) return;

   TGrid *g = TGrid::Connect(fAlienHost, "");

   //Protection in case root is compiled without AliEn support
   if(!g) {
      Error("UpdateAliEn", "ROOT compiled without AliEn support");
      return;
   }

   TString lfn = fAlienDir;
   TDatime dt;

   // make a subdirectory for each day
   lfn += "/adc-";
   lfn += dt.GetDate();

   // check if directory exists, if not create it
#if ROOT_VERSION_CODE < ROOT_VERSION(5,0,0)
   Grid_ResultHandle_t res = 0;
   if (!(res = g->OpenDir(lfn))) {
      // directory does not exist, create it
      if (g->Mkdir(lfn) == -1) {
         Error("UpdateAliEn", "cannot create directory %s", lfn.Data());
         lfn = fAlienDir;
      }
   }
   if (res) g->CloseResult(res);
#else
   Error("UpdateAliEn", "needs to be ported to new TGrid");
#endif

   lfn += "/";
   lfn += gSystem->BaseName(stats->GetFileName());

#if ROOT_VERSION_CODE < ROOT_VERSION(5,0,0)
   Int_t result = g->AddFile(lfn, stats->GetFileName(),
			     (int)stats->GetFileSize());
   if (result == -1) {
      Error("UpdateAliEn", "error adding file to AliEn catalog");
      printf("AliEn: AddFile(%s, %s, %d)\n", lfn.Data(), stats->GetFileName(),
             (int)stats->GetFileSize());
   }
#else
   Error("UpdateAliEn", "needs to be ported to new TGrid");
#endif

   delete g;
}

//______________________________________________________________________________
void AliRunDB::Close()
{
   // Close run database.

   if (fRunDB) fRunDB->Close();
   delete fRunDB;
}

