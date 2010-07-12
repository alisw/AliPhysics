#include "AliDimIntNotifier.h"
#include "AliOnlineReco.h"

#include <TRint.h>

#include <TSQLServer.h>
#include <TSQLResult.h>
#include <TSQLRow.h>

#include <TTimeStamp.h>

int main(int argc, char **argv)
{
  AliDimIntNotifier::SetMainThreadId();

  Bool_t test = argc > 1 && strcmp("-test", argv[1]) == 0;

  TRint app("App", &argc, argv);

  AliOnlineReco *win = new AliOnlineReco;
  win->MapWindow();

  TString autoRun(gSystem->Getenv("ONLINERECO_AUTORUN"));
  if (autoRun == "1" || autoRun.CompareTo("true", TString::kIgnoreCase) == 0)
  {
    win->SetAutoRunMode(kTRUE);
  }

  if (test)
  {
    win->SetTestMode();

    win->GetSOR(0)->infoHandlerTest(2214);
    win->GetSOR(0)->infoHandlerTest(2215);
    win->GetSOR(0)->infoHandlerTest(2224);
    win->GetSOR(0)->infoHandlerTest(2244);

    printf("win = (AliOnlineReco*) 0x%lx\n", (unsigned long)win);
  }
  else
  {
    TString baseDir = gSystem->Getenv("ONLINERECO_BASE_DIR");
    if (baseDir.IsNull())
    {
      printf("ERROR: ONLINERECO_BASE_DIR is not set. Exiting...");
      return 0;
    }

    const char *dbHost = "aldaqdb";
    Int_t dbPort = 3306;
    const char *dbName = "LOGBOOK";
    const char *user = "logbook";
    const char *password = "alice";

    TSQLServer* server = TSQLServer::Connect(Form("mysql://%s:%d/%s", dbHost, dbPort, dbName), user, password);
    if (!server) {
      printf("ERROR: Could not connect to DAQ Logbook");
      return 0;
    }
    TString sqlQuery;
    TTimeStamp ts;
    sqlQuery.Form("SELECT run FROM logbook WHERE DAQ_time_start > %u AND DAQ_time_end IS NULL AND partition REGEXP 'PHYSICS.*'",
		  (UInt_t)ts.GetSec()-86400);
    TSQLResult* result = server->Query(sqlQuery);
    if (!result)
    {
      printf("ERROR: Can't execute query <%s>!", sqlQuery.Data());
      return 0;
    }
    if (result->GetRowCount() == 0)
    {
      printf("No active physics runs found");
    }
    else
    {
      for (Int_t iRow = 0; iRow < result->GetRowCount(); iRow++)
      {
	TSQLRow* row = result->Next();
	TString runStr = row->GetField(0);
	if (runStr.IsDigit())
	  win->StartOfRun(runStr.Atoi());
	delete row;
      }
    }
    delete result;
  }

  app.Run(kTRUE);
  return 0;
}
