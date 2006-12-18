// @(#)alimdc:$Name$:$Id$
// Author: Fons Rademakers  26/11/99

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// alimdc                                                               //
//                                                                      //
// Main program used to create application that reads a data stream     //
// from the DATE DAQ system and that creates a ROOT database.           //
//                                                                      //
// Written by: Fons Rademakers, 1/4/99.                                 //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdlib.h>
#include <fcntl.h>
#include <errno.h>

#include <TROOT.h>
#include <TSystem.h>
#include <TError.h>

#include <AliLog.h>

#ifdef USE_SMI
extern "C" {
   #include <smirtl.h>
}
#endif

#include "AliMDC.h"

#ifdef __APPLE__
// avoid loading pythia and pdf
#include <Hepevt.h>
#endif

//______________________________________________________________________________
static void AliMDCErrorHandler(int level, Bool_t abort, const char *location,
                               const char *msg)
{
   // The default error handler function. It prints the message on stderr and
   // if abort is set it aborts the application. Comapared to the default
   // ROOT error handler this one also prints the date and time in front
   // of each message.

   if (level < gErrorIgnoreLevel)
      return;

   const char *type = 0;

   if (level >= kInfo)
      type = "Info";
   if (level >= kWarning)
      type = "Warning";
   if (level >= kError)
      type = "Error";
   if (level >= kBreak)
      type = "\n *** Break ***";
   if (level >= kSysError)
      type = "SysError";
   if (level >= kFatal)
      type = "Fatal";

   TDatime dt;

   if (level >= kBreak && level < kSysError)
      fprintf(stderr, "%s: %s %s\n", dt.AsSQLString(), type, msg);
   else if (!location || strlen(location) == 0)
      fprintf(stderr, "%s: %s: %s\n", dt.AsSQLString(), type, msg);
   else
      fprintf(stderr, "%s: %s in <%s>: %s\n", dt.AsSQLString(), type, location,
              msg);

   fflush(stderr);
   if (abort) {
      fprintf(stderr, "aborting\n");
      fflush(stderr);
      if (gSystem) {
         gSystem->StackTrace();
         gSystem->Abort();
      } else
         ::abort();
   }
}

#ifdef USE_SMI
static void SMI_handle_command()
{
   // Handle SMI commands

   char action[64], param[64];
   int n_params;

   smi_get_action(action, &n_params);
   if (n_params >= 1) {
      smi_get_par_value("PARAM", param);
   } else {
      strcpy(param, "");
   }
   if (strcmp(action, "STOP") == 0) {
      if (AliMDC::Instance()) AliMDC::Instance()->SetStopLoop();
   }
   smi_set_state("RUNNING");
}
#endif

//______________________________________________________________________________
static void Usage(const char *prognam)
{
#ifdef USE_SMI
      fprintf(stderr, "Usage: %s <sminame> <dbsize> <tagdbsize> <filter> <compmode> [date_file]\n",
              prognam);
      fprintf(stderr, " <sminame> = name used by SMI\n");
#else
      fprintf(stderr, "Usage: %s <dbsize> <tagdbsize> <filter> <compmode> [date_file]\n",
              prognam);
#endif
      fprintf(stderr, " <dbsize> = maximum raw DB size (in bytes)\n");
      fprintf(stderr, "    (precede by - to delete raw and tag databases on close)\n");
      fprintf(stderr, " <tagdbsize> = maximum tag DB size (in bytes, 0 for no tag DB)\n");
      fprintf(stderr, " <filter> = state of 3rd level filter (0: off, 1: transparent, 2: on)\n");
      fprintf(stderr, " <compmode> = compression level (see TFile)\n");
      fprintf(stderr, "    (precede by - to use RFIO, -0 is RFIO and 0 compression)\n");
      fprintf(stderr, "    (precede by + to use rootd, +0 is rootd and 0 compression)\n");
      fprintf(stderr, "    (precede by %% to use Castor/rootd, %%0 is Castor/rootd and 0 compression)\n");
      fprintf(stderr, "    (precede by @ to use /dev/null as sink)\n");
      fprintf(stderr, " [date_file] = optional input file (default reads from DATE EventBuffer)\n");
      fprintf(stderr, "    (precede with - for endless loop on same file (use SIGUSR1 to stop)\n");
}

//______________________________________________________________________________
int main(int argc, char **argv)
{
   // Convert a DATE data stream to a ROOT DB.

   // Set ROOT in batch mode
   gROOT->SetBatch();

   // Set custom error handler
   AliLog::SetHandleRootMessages(kFALSE);
   SetErrorHandler(AliMDCErrorHandler);

   // Default file system locations
#ifdef USE_EB
   const char* rawDBFS[2] = { "/data1/mdc", "/data2/mdc" };
   const char* tagDBFS    = "/data1/mdc/tags";
   const char* rfioFS     = "rfio:/castor/cern.ch/lcg/dc5";
   const char* castorFS   = "castor:/castor/cern.ch/lcg/dc5";
#else
   const char* rawDBFS[2] = { "/tmp/mdc1", "/tmp/mdc2" };
   const char* tagDBFS    = "/tmp/mdc1/tags";
   TString user(gSystem->Getenv("USER")[0] + TString("/") + 
		gSystem->Getenv("USER"));
   TString rfioStr("rfio:/castor/cern.ch/user/" + user);
   const char* rfioFS     = rfioStr.Data();
   TString castorStr("castor:/castor/cern.ch/user/" + user);
   const char* castorFS   = castorStr.Data();
#endif
   const char* rootdFS    = "root://localhost//tmp/mdc1";

   // User defined file system locations
   if (gSystem->Getenv("ALIMDC_RAWDB1")) 
     rawDBFS[0] = gSystem->Getenv("ALIMDC_RAWDB1");
   if (gSystem->Getenv("ALIMDC_RAWDB2")) 
     rawDBFS[1] = gSystem->Getenv("ALIMDC_RAWDB2");
   if (gSystem->Getenv("ALIMDC_TAGDB")) 
     tagDBFS = gSystem->Getenv("ALIMDC_TAGDB");
   if (gSystem->Getenv("ALIMDC_RFIO")) 
     rfioFS = gSystem->Getenv("ALIMDC_RFIO");
   if (gSystem->Getenv("ALIMDC_CASTOR")) 
     castorFS = gSystem->Getenv("ALIMDC_CASTOR");
   if (gSystem->Getenv("ALIMDC_ROOTD")) 
     rootdFS = gSystem->Getenv("ALIMDC_ROOTD");

   // Handle command line arguments
   if ((argc == 2 && (!strcmp(argv[1], "-?") || !strcmp(argv[1], "-help"))) ||
#ifdef USE_SMI
       argc > 7 || argc < 6) {
#else
       argc > 6 || argc < 5) {
#endif
      Usage(argv[0]);
      return 1;
   }

   Int_t iarg = 1;
#ifdef USE_SMI
   char smiobj[128];
   strcpy(smiobj, argv[iarg]);
   smi_attach(smiobj, SMI_handle_command);
   smi_volatile();
   smi_set_state("RUNNING");
   iarg++;
#endif

   AliMDC::EWriteMode wmode = AliMDC::kLOCAL;
   Int_t    filterMode = 0;
   Bool_t   useLoop = kFALSE;
   Bool_t   delFiles = kFALSE;
   Int_t    compress;
   Double_t maxFileSize;
   Double_t maxTagSize;
   const char* fs1 = NULL;
   const char* fs2 = NULL;

   // no special arg checking so don't make errors
   if (argv[iarg][0] == '-') {
      delFiles = kTRUE;
      maxFileSize = atoi(argv[iarg]+1);
   } else
      maxFileSize = atoi(argv[iarg]);
   if (maxFileSize < 1000 || maxFileSize > 2.e9) {
      Error(argv[0], "unreasonable file size %f\n", maxFileSize);
      return 1;
   }
   iarg++;

   maxTagSize = atoi(argv[iarg]);
   if (maxTagSize > 0 && (maxTagSize < 1000 || maxTagSize > 2.e9)) {
      Error(argv[0], "unreasonable tag file size %f\n", maxTagSize);
      return 1;
   }
   if (maxTagSize == 0) tagDBFS = NULL;
   iarg++;

   filterMode = atoi(argv[iarg]);
   if (filterMode < 0 || filterMode > 2) {
      Error(argv[0], "unreasonable filter mode %d\n", filterMode);
      return 1;
   }
   iarg++;

   if (argv[iarg][0] == '-') {
      wmode = AliMDC::kRFIO;
      compress = atoi(argv[iarg]+1);
      fs1 = rfioFS;
   } else if (argv[iarg][0] == '+') {
      wmode = AliMDC::kROOTD;
      compress = atoi(argv[iarg]+1);
      fs1 = rootdFS;
   } else if (argv[iarg][0] == '%') {
      wmode = AliMDC::kCASTOR;
      compress = atoi(argv[iarg]+1);
      fs1 = castorFS;
   } else if (argv[iarg][0] == '@') {
      wmode = AliMDC::kDEVNULL;
      compress = atoi(argv[iarg]+1);
   } else {
      compress = atoi(argv[iarg]);
      fs1 = rawDBFS[0];
      fs2 = rawDBFS[1];
   }
   if (compress > 9) {
      Error(argv[0], "unreasonable compression mode %d\n", compress);
      return 1;
   }
   iarg++;

   char* file = NULL;
   if (iarg < argc) {
      file = argv[iarg];
      if (argv[iarg][0] == '-') {
         useLoop = kTRUE;
         file = argv[iarg]+1;
      }
   }

   // Create MDC processor object and process input stream
   AliMDC mdcproc(compress, delFiles, AliMDC::EFilterMode(filterMode), 
		  maxTagSize, tagDBFS);

   Int_t result = mdcproc.Run(file, useLoop, wmode, maxFileSize, fs1, fs2);

   if (result == 0)
      Info(argv[0], "normal termination of run");
   else
      Error(argv[0], "error termination of run, status: %d", result);
   return result;
}
