#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AliCDBManager.h"
#include "AliGRPManager.h"
#include "AliGRPObject.h"
#include "AliRunInfo.h"
#include "AliTimeStamp.h"
#include <ctime>
#include <fstream>
#include <iostream>

#endif

using namespace std;

// Macro to produce the start and end timestamp of runs
// Runs can be given from an input file (comma separated runs) or directly as input to the function e.g.
// aliroot -l -b -q 'GetStartAndEndOfRunTimestamp.C("280235")'
// aliroot -l -b -q 'GetStartAndEndOfRunTimestamp.C("list_of_runs.txt")'
// Timestamps are written to the specified file
// Duplicate entries are not checked
// Author: Nicolo' Jacazio

TString ProcessRun(const AliGRPObject* grpObj, time_t& start, time_t& end);

void GetStartAndEndOfRunTimestamp(TString runlist = "280235", TString outputfile = "/tmp/RunStartTimes.txt")
{
  TObjArray* run_arr = 0x0;
  if (runlist.Contains(".")) { // Input file with runs mode
    run_arr = new TObjArray();
    Printf("Using input file runlist %s", runlist.Data());
    std::ifstream infile(runlist);
    std::string line;
    while (std::getline(infile, line)) {
      TString l = line;
      l = l.Strip(TString::kBoth, ' ');
      TString separator = ",";
      TObjArray* line_arr = l.Tokenize(separator);
      for (Int_t i = 0; i < line_arr->GetEntries(); i++) {
        TString n = line_arr->At(i)->GetName();
        n = n.Strip(TString::kBoth, ' ');
        Printf("Adding run %i %s to the list of retrieval", run_arr->GetEntries(), n.Data());
        if (!n.IsDigit()) {
          Printf("Error: %s is not digit!", n.Data());
          return;
        }
        run_arr->AddLast(new TNamed(n, n));
      }
    }
  } else { // Run list mode
    run_arr = runlist.Tokenize(" ");
    for (Int_t i = 0; i < run_arr->GetEntries(); i++) {
      Printf("Added run %i %s to the list of retrieval", i, run_arr->At(i)->GetName());
    }
  }
  //
  Printf("\nConnecting to OCDB");
  AliCDBManager* CDBman = AliCDBManager::Instance();
  CDBman->SetDefaultStorage("raw://");
  for (Int_t i = 0; i < run_arr->GetEntries(); i++) {
    TString run = run_arr->At(i)->GetName();
    CDBman->SetRun(run.Atoi());
    AliGRPManager grpMan;
    Bool_t status = grpMan.ReadGRPEntry();
    const AliGRPObject* grpObj = grpMan.GetGRPData();
    time_t start, end;
    TString s = ProcessRun(grpObj, start, end);
    Printf("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
    Printf("> Run %s %s <", run.Data(), s.Data());
    Printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
    // Writing to file
    ofstream out_file;
    out_file.open(outputfile, ios::out | ios::app);
    out_file << "Run " << run << " starts at " << start << " and ends at " << end << " (" << s << ")\n";
    out_file.close();
  }
  Printf("Timestamps written to '%s'", outputfile.Data());
}

TString ProcessRun(const AliGRPObject* grpObj, time_t& start, time_t& end)
{
  start = grpObj->GetTimeStart();
  std::tm* start_time = std::localtime(&start);
  char buffer_start[32];
  std::strftime(buffer_start, sizeof(buffer_start), "%a, %d.%m.%Y %H:%M:%S", start_time);
  end = grpObj->GetTimeEnd();
  std::tm* end_time = std::localtime(&end);
  char buffer_end[32];
  std::strftime(buffer_end, sizeof(buffer_end), "%a, %d.%m.%Y %H:%M:%S", end_time);
  TString s = Form("started at %s ended at %s", buffer_start, buffer_end);
  return s;
}
