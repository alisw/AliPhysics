/**
 * @file   GridWatch.C
 * @author Christian Holm Christensen <cholm@master.hehi.nbi.dk>
 * @date   Thu Jan 24 23:06:08 2013
 * 
 * @brief Script to watch master jobs and automatically submit
 * terminate jobs
 * 
 * 
 * @ingroup pwglf_forward_trains_helper
 */
#ifndef __CINT__
# include <TString.h>
# include <TGrid.h>
# include <TSystem.h>
# include <TObjArray.h>
# include <iostream>
# include <fstream>
# include <TError.h>
# include <TDatime.h>
# include <TEnv.h>
#else 
class TString;
#endif
#include <TArrayI.h>

/** 
 * Create token name 
 * 
 * @param name   Name
 * @param ext    Extension
 * @param merge  Merge state or not 
 * 
 * @return Formatted string
 * @ingroup pwglf_forward_trains_helper
 */
TString TokenName(const TString& name, 
		   const TString& ext, 
		   Bool_t merge=false)
{
  return TString::Format("%s%s.%s", name.Data(), 
			 (merge ? "_merge" : ""), ext.Data());
}
/** 
 * Check if we have a particular file 
 * 
 * @param name   Base name 
 * @param ext    Extension 
 * @param merge  Merging stage or not 
 * 
 * @return true if file exits
 * @ingroup pwglf_forward_trains_helper
 */
Bool_t CheckTokens(const TString& name, 
		   const TString& ext, 
		   Bool_t merge=false)
{
  // TSystem::AccessPathName return false if file is there 
  return !gSystem->AccessPathName(TokenName(name, ext, merge));
}
/** 
 * Remove a token file 
 * 
 * @param name   Base name 
 * @param ext    Extension 
 * @param merge  Merging stage or not 
 * @ingroup pwglf_forward_trains_helper
 */
void RemoveTokens(const TString& name, 
		  const TString& ext, 
		  Bool_t merge=false)
{
  gSystem->Unlink(TokenName(name, ext, merge));
}

/** 
 * Read one line of text from file and return tokens.
 * 
 * @param name   Base name 
 * @param ext    Extension
 * @param merge  If true append "_merge" to name
 * 
 * @return Array of tokens or null
 *
 * @ingroup pwglf_forward_trains_helper
 */
TObjArray* ReadTokens(const TString& name, 
		      const TString& ext, 
		      bool merge=false) 
{
  TString fn = TString::Format("%s%s.%s", name.Data(), 
			       (merge ? "_merge" : ""), ext.Data());
  std::ifstream in(fn.Data());
  if (!in) { 
    Error("ReadTokens", "Failed to open %s", fn.Data());
    return 0;
  }
  TString ln;
  ln.ReadLine(in);
  in.close();
  
  if (ln.IsNull()) return 0;
  return ln.Tokenize(" \t");
}
  
/** 
 * Read list of job IDs from file 
 * 
 * @param name   Base name 
 * @param merge  If true append "_merge" to name
 * 
 * @return Array of job IDs or null
 *
 * @ingroup pwglf_forward_trains_helper
 */
TObjArray* ReadJobIDs(const TString& name, bool merge=false)
{
  return ReadTokens(name, "jobid", merge);
}

/** 
 * Read list of job stages from file 
 * 
 * @param name   Base name 
 * @param merge  If true append "_merge" to name
 * 
 * @return Array of job stages or null
 *
 * @ingroup pwglf_forward_trains_helper
 */
TObjArray* ReadStages(const TString& name, bool merge=false)
{
  return ReadTokens(name, "stage", merge);
}

/** 
 * Parse the job IDs into an array of integers 
 * 
 * @param jobIds List of jobs
 * @param ret    Return array
 * 
 * @return true on success
 *
 * @ingroup pwglf_forward_trains_helper
 */
Bool_t ParseJobIDs(const TObjArray* jobIds, TArrayI& ret)
{
  if (!jobIds) return false;

  Int_t n = jobIds->GetEntries();
  ret.Set(n);
  ret.Reset(-1);

  for (Int_t i = 0; i < n; i++) { 
    TObjString*    id = static_cast<TObjString*>(jobIds->At(i));
    const TString& s  = id->String();
    ret.SetAt(s.Atoi(), i);
  }
  return true;
}

/** 
 * Parse string representing status and return human-readable string 
 * 
 * @param status Return from ps command
 * @param out    Output
 * 
 * @return true on success
 *
 * @ingroup pwglf_forward_trains_helper
 */
Bool_t ParseState(const TString& status, TString& out)
{
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
  default:  out = "UNKNOWN"; return false;
  }
  if (status[1] != '\0' && 
      (status[0] != 'O' || status[0] != 'S')) { 
    out.Append("_");
    switch (status[1]) { 
    case 'S': out.Append(status[0] == 'E' ? "SUBMIT" : "SPLIT"); break;
    case 'X': out.Append("EXPIRED"); break;
    case 'A': out.Append("ASSIGNING"); break;
    case 'E': out.Append("EXECUTING"); break;
    case 'V': 
      if (status[0] == 'S') out = "SAVING"; 
      else out.Append("VALIDATING"); 
      break;
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
    default:  out.Append("UNKNOWN"); return false;
    }
    if (status[2] != '\0') { 
      switch (status[2]) {
      case 'V': if (status[0] == 'E') out.ReplaceAll("SUBMIT", "SAVING");
	break;
      default: out.Append("_UNKNOWN");
      }
    }
  }
  return true;
}  

/** 
 * Get the job state 
 * 
 * @param jobId Job status 
 * @param out   Output status string 
 * 
 * @return true on success
 *
 * @ingroup pwglf_forward_trains_helper
 */
Bool_t GetJobState(Int_t jobId, TString& out) 
{
  out = "MISSING";

  TString fn("gridMonitor");
  FILE* fp = gSystem->TempFileName(fn);

  gSystem->RedirectOutput(fn);
  gGrid->Command("ps -Ax");
  gGrid->Stdout();
  gSystem->RedirectOutput(0);
  gGrid->Stderr();

  fclose(fp);

  std::ifstream in(fn.Data());

  while (!in.eof()) { 
    TString l;
    l.ReadLine(in);
    if (in.bad()) break;

    TObjArray* tokens = l.Tokenize(" \t");
    if (tokens->GetEntries() < 2) break;

    //TString  user   = tokens->At(0)->GetName(); 
    TString    sjid   = tokens->At(1)->GetName(); // Job ID
    TString    stat   = tokens->At(2)->GetName(); // State 
    Int_t      jid    = sjid.Atoi();
    
    if (jid != jobId) continue;

    ParseState(stat, out);
    break;
  }

  in.close();
  gSystem->Unlink(fn);

  return true;
}

/** 
 * Get the job states
 * 
 * @param jobs   List of job IDs
 * @param states On return the states
 * 
 * @return true on success
 *
 * @ingroup pwglf_forward_trains_helper
 */
Bool_t GetJobStates(const TArrayI& jobs, TObjArray& states)
{
  Int_t n = jobs.GetSize();
  states.Expand(n);
  for (Int_t i = 0; i < n; i++) {
    TObjString* s = static_cast<TObjString*>(states.At(i));
    if (!s) states.AddAt(s = new TObjString(""), i);
    s->SetString("MISSING");
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

  while (!in.eof()) {
    TString l;
    l.ReadLine(in);
    if (in.bad()) break;
    if (l.IsNull()) continue;

    TObjArray* tokens = l.Tokenize(" \t");
    if (tokens->GetEntries() < 3) { 
      Warning("GetJobStates", "Got too few tokens (%d): %s", 
	      tokens->GetEntries(), l.Data());
      tokens->Print();
      break;
    }

    //TString  user   = tokens->At(0)->GetName(); 
    TString    sjid   = tokens->At(1)->GetName(); // Job ID
    TString    stat   = tokens->At(2)->GetName(); // State 
    Int_t      jid    = sjid.Atoi();
    
    for (Int_t i = 0; i < n; i++) { 
      if (jid != jobs.At(i)) continue;
      TObjString* s = static_cast<TObjString*>(states.At(i));
      TString out;
      if (!ParseState(stat, out)) continue;
      s->SetString(out);
    }
  }

  in.close();
  gSystem->Unlink(fn);

  return true;
}


/** 
 * Wait of jobs to finish 
 * 
 * @param jobs    List of jobs
 * @param stages  Stages
 * @param delay   Delay for check
 * 
 * @return true on success, false otherwise
 *
 * @ingroup pwglf_forward_trains_helper
 */
Bool_t WaitForJobs(TArrayI& jobs, TObjArray* stages, Int_t delay)
{
  Bool_t stopped = false;
  TFileHandler h(0, 0x1);
  do { 
    Bool_t allDone = true;
    TDatime t;
    Printf("--- %4d/%02d/%02d %02d:%02d:%02d [Press enter to pause] ---", 
	   t.GetYear(), t.GetMonth(), t.GetDay(), 
	   t.GetHour(), t.GetMinute(), t.GetSecond());

    TObjArray states;
    GetJobStates(jobs, states);

    Int_t missing = 0;
    Int_t total   = jobs.GetSize();
    // Bool_t allAccounted = false;
    for (Int_t i = 0; i < total; i++) { 
      Int_t job = jobs.At(i);

      if (job < 0) continue;

      TObjString* obj = static_cast<TObjString*>(states.At(i));
      const TString& state = obj->String();
      
      if (state.BeginsWith("ERROR_"))
	jobs.SetAt(-1, i);
      else if (state.EqualTo("MISSING")) 
	missing++;
      else if (!state.EqualTo("DONE")) 
	allDone = false;
      

      Printf(" %d(%s)=%s", job, stages->At(i)->GetName(), state.Data());
      
    }
    if (allDone) break;
    if (missing >= total) {
      Error("GetJobStates", "Info on all jobs missing");
      break;
    }
    if (gSystem->Select(&h, 1000*delay)) {
      // Got input on std::cin 
      std::string l;
      std::getline(std::cin, l);
      std::cout << "Do you want to terminate now [yN]? " << std::flush;
      std::getline(std::cin, l);
      if (l[0] == 'y' || l[0] == 'Y') { 
	stopped = true;
	break;
      }
    }

    // gSystem->Sleep(1000*delay);
  } while (true);

  return true;
}
/** 
 * Watch Grid for termination of main job, and submit merging jobs as needed. 
 * 
 * @param name   Name of the job
 * @param delay  Delay between updates in seconds
 *
 * @ingroup pwglf_forward_trains_helper
 */
void GridWatch(const TString& name, UShort_t delay=5*60)
{
  gEnv->SetValue("XSec.GSI.DelegProxy", "2");
  TGrid::Connect("alien:///");
  if (!gGrid) { 
    Error("GridWatch", "Failed to connect to the Grid");
    return;
  }

  TObjArray* jobIDs = ReadJobIDs(name, false);
  TObjArray* stages = ReadStages(name, false);

  if (!jobIDs || !stages) return;

  TArrayI jobs;
  if (!ParseJobIDs(jobIDs, jobs)) return;

  if (!(CheckTokens(name, "jobid", true) && 
	CheckTokens(name, "stage", true))) 
    WaitForJobs(jobs, stages, delay);

  delete jobIDs;
  delete stages;

  // return;
  do {
    if (!CheckTokens(name, "jobid", true) && 
	!CheckTokens(name, "stage", true)) {
      Printf("Now executing terminate");
      gSystem->Exec("aliroot -l -b -q Terminate.C");
      gSystem->Sleep(10*1000);
    }
    Printf("Reading job ids");

    jobIDs = ReadJobIDs(name, true);
    stages = ReadStages(name, true);
    
    if (!ParseJobIDs(jobIDs, jobs)) {
      Error("GridWatch", "Failed to parse job ids %s", 
	    TokenName(name,"jobid",true).Data());
      return;
    }

    WaitForJobs(jobs, stages, delay);
    
    Bool_t allFinal = true;
    for (Int_t i = 0; i < jobs.GetSize(); i++) {
      if (jobs.At(i) < 0) continue;

      const TString& s = static_cast<TObjString*>(stages->At(i))->String();
      if (!s.BeginsWith("final_")) allFinal = false;
    }
    
    delete jobIDs;
    delete stages;

    Printf("All jobs in final stage");
    if (allFinal) break;

    RemoveTokens(name, "jobid", true);
    RemoveTokens(name, "stage", true);
  } while (true);

  Printf("Finished");
}
//
// EOF
//

