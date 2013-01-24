#ifndef __CINT__
# include <TString.h>
# include <TGrid.h>
# include <TObjArray.h>
# include <TSystem.h>
# include <TError.h>
# include <TEnv.h>
# include <fstream>
# include <cstdio>
#else
class TString;
#endif

void ShowJobState(const TString& status, Int_t jid=-1)
{
  TString out;
  switch (status[0]) { 
  case 'D': out = "DONE"; break;
  case 'E': out = "ERROR"; break;
  case 'R': out = "RUNNING"; break; 
  case 'W': out = "WAITING"; break;
  case 'O': out = "WAITING_QUOTA"; break;
  case 'A': out = "ASSIGNED"; break;
  case 'S': out = "STARTED" ; break;
  case 'I': out = "INSERTING"; break;
  case 'K': out = "KILLED"; break;
  default:  out = "UNKNOWN"; break;
  }
  if (status[1] != '\0' && 
      (status[0] != 'O' || status[0] != 'S')) { 
    out.Append("_");
    switch (status[1]) { 
    case 'S': out.Append(status[0] == 'E' ? "SUBMIT" : "SPLIT"); break;
    case 'X': out.Append("EXPIRED"); break;
    case 'A': out.Append("ASSIGNING"); break;
    case 'E': out.Append("EXECUTING"); break;
    case 'V': out.Append("VALIDATING"); break;
    case 'd': 
      if (status[0] == 'I') {
	out = "INTERACTIVE_IDLE";
	break;
      } // Fall through on else 
    case 'a': 
      if (status[0] == 'I') {
	out = "INTERACTIVE_USED";
	break;
      } // Fall through on else
    default:  out.Append("UNKNOWN"); break;
    }
  }
  if (jid > 0) Printf("%-9d %s", jid, out.Data());
  else         Printf("%s", out.Data());
}

void GridMonitor(Int_t jobId=0)
{
  gSystem->RedirectOutput("/dev/null");
  gEnv->SetValue("XSec.GSI.DelegProxy", "2");
  TGrid::Connect("alien:///");
  gSystem->RedirectOutput(0);
  if (!gGrid) { 
    Error("GridMonitor", "Failed to connect");
    return;
  }

  TString fn("gridMonitor");
  FILE* fp = gSystem->TempFileName(fn);
  gSystem->RedirectOutput(fn);
  gGrid->Command("ps -Ax");
  gGrid->Stdout();
  gSystem->RedirectOutput(0);
  gGrid->Stderr();
  fclose(fp);

  std::ifstream in(fn.Data());
  bool found = false;

  while (!in.eof()) { 
    TString l;
    l.ReadLine(in);
    if (in.bad()) break;

    TObjArray* tokens = l.Tokenize(" \t");
    if (tokens->GetEntriesFast() < 2) break;

    TString    user   = tokens->At(0)->GetName(); // user name 
    TString    sjid   = tokens->At(1)->GetName(); // Job ID
    TString    stat   = tokens->At(2)->GetName(); // State 
    
    Int_t      jid    = sjid.Atoi();
    
    if (jobId > 0) {
      if (jid != jobId) continue;
      found = true;
    }
    ShowJobState(stat, jobId > 0 ? -1 : jid);
    if (found) break;
  }
  in.close();
}

      
      
  
