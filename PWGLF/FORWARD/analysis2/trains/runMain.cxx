/**
 * @file   runMain.cxx
 * @author Christian Holm Christensen <cholm@master.hehi.nbi.dk>
 * @date   Fri Jun  1 13:56:29 2012
 * 
 * @brief  
 * 
 * @ingroup pwglf_forward_trains
 * 
 */
/** 
 * @defgroup pwglf_forward_trains_run Execute train
 *
 * Code for program to run trains 
 *
 * @ingroup pwglf_forward_trains
 */
#include "FORWARD/analysis2/trains/TrainSetup.C"
#include <TSystem.h>
#include <TROOT.h>
#include <TString.h>

/** 
 * Build the train setup script 
 * 
 * @param script Script name 
 * @param useTmp Whether to use temporary directory
 * @param useDbg Build class with debug symbols
 * 
 * @return true on success
 * 
 * @ingroup pwglf_forward_trains_run
 */
Bool_t BuildTrain(const char* script, 
		  Bool_t useTmp=false,
		  Bool_t useDbg=false)
{
  // --- Load needed libraries ---------------------------------------
  gROOT->LoadClass("TAlien",               "libRAliEn");
  gROOT->LoadClass("TVirtualMC",           "libVMC");
  gROOT->LoadClass("AliVEvent",            "libSTEERBase");
  gROOT->LoadClass("AliESDEvent",          "libESD");
  gROOT->LoadClass("AliAODEvent",          "libAOD");
  gROOT->LoadClass("AliAnalysisManager",   "libANALYSIS");
  gROOT->LoadClass("AliAnalysisTaskSE",    "libANALYSISalice");

  // --- Setup script path -------------------------------------------
  const char* aliPath   = gSystem->ExpandPathName("$ALICE_ROOT");
  const char* fwd2Path  = 
    gSystem->ExpandPathName("$ALICE_ROOT/PWGLF/FORWARD/analysis2");
  const char* macroPath = gROOT->GetMacroPath();
  gROOT->SetMacroPath(Form(".:%s:%s/trains:%s/scripts",
			   macroPath,fwd2Path,fwd2Path));
  
  // --- Setup include path ------------------------------------------
  gSystem->AddIncludePath(Form("-I%s", fwd2Path));
  gSystem->AddIncludePath(Form("-I%s/trains", fwd2Path));
  gSystem->AddIncludePath(Form("-I%s", aliPath));
  gSystem->AddIncludePath(Form("-I%s/include", aliPath));

  // --- Check that we can find the script ---------------------------
  TString path = gSystem->Which(gROOT->GetMacroPath(), script);
  if (path.IsNull()) { 
    path = gSystem->Which(gROOT->GetMacroPath(), Form("%s.C", script));
    if (path.IsNull()) { 
      Error("BuildTrain", "Cannot find script %s or %s.C in %s", 
	    script, script, gROOT->GetMacroPath());
      return false;
    }
  }
  // --- Possibly make a temporary copy ------------------------------
  if (useTmp) { 
    TString tmp("Train");
    FILE*   fp = gSystem->TempFileName(tmp, ".");
    fclose(fp);
    gSystem->Unlink(tmp);
    tmp.Append(".C");
    
    gSystem->CopyFile(path, tmp);
    path = tmp;
  }

  // --- Now compile the thing ---------------------------------------
  // Info("BuildTrain", "Building macro %s (ACLiC mode: +%s)", path.Data(),
  //      useDbg ? "+g" : "");
  Int_t ret = gROOT->LoadMacro(Form("%s+%s", path.Data(), useDbg ? "+g" : ""));

  // --- If we did a temporary file, remove stuff --------------------
  if (useTmp) {
    gSystem->Unlink(path);
    path.ReplaceAll(".C",   "_C.d");  gSystem->Unlink(path);
    path.ReplaceAll("_C.d", "_C.so"); gSystem->Unlink(path);
  }

  // --- Warning in case of problems ---------------------------------
  if (ret != 0) 
    Warning("BuildTrain", "Failed to build %s (%s)", path.Data(), script);

  return ret == 0;  
}
/** 
 * Make the train setup object 
 * 
 * @param trainClass Class-name of the train setup 
 * @param trainName  Name of the train. 
 * 
 * @return Newly allocated object or null
 * 
 * @ingroup pwglf_forward_trains_run
 */
TrainSetup* MakeTrain(const char* trainClass, const char* trainName)
{
  // Info("MakeTrain", "Making an object of class %s named '%s'", 
  //      trainClass, trainName);
  Long_t ret = gROOT->ProcessLine(Form("new %s(\"%s\")",trainClass,trainName));
  if (!ret) 
    Warning("MakeTrain", "Failed to make object of class %s", trainClass);

  return reinterpret_cast<TrainSetup*>(ret);
}

/** 
 * Append an element to a list 
 * 
 * @param str String to append to 
 * @param val Value to append
 * @param sep Separator
 * 
 * @ingroup pwglf_forward_trains_run
 */
void AppendTo(TString& str, const TString& val, const char* sep=",")
{
  if (!str.IsNull()) str.Append(sep);
  // Info("", "Appending %s to %s", val.Data(), str.Data());
  str.Append(val);
}
/** 
 * Append directory to header and script search path
 * 
 * @param dir Directory
 * 
 * @ingroup pwglf_forward_trains_run
 */
void AppendPath(const char* dir)
{
  gROOT->SetMacroPath(Form("%s:%s",gROOT->GetMacroPath(), dir));
  gSystem->AddIncludePath(Form("-I%s", dir));
}
/** 
 * Print a fake option description.  Used for options specific to this
 * program.
 * 
 * @param o    Output stream 
 * @param opt  Option (including meta argument)
 * @param desc Option description.
 * 
 * @ingroup pwglf_forward_trains_run
 */
void PrintFakeOption(std::ostream& o, const char* opt, const char* desc)
{
  o << "  --" << std::left << std::setw(30) << opt << " " << desc << std::endl;
}

/** 
 * Print usage information 
 * 
 * @param progname Program name 
 * @param o        Output stream
 * @param r        Optional runner. 
 * 
 * @ingroup pwglf_forward_trains_run
 */
void Usage(const char* progname, std::ostream& o, TrainSetup::Runner* r)
{
  o << "Usage: " << progname << " --class=CLASS --name=NAME [OPTIONS]\n\n"
    << "PROGRAM OPTIONS:\n";
  PrintFakeOption(o, "class=CLASS", "Train class");
  PrintFakeOption(o, "name=NAME",   "Name of train");
  PrintFakeOption(o, "events=NUMBER", "Max number of events to analyse");
  PrintFakeOption(o, "run=RUN",       "Comma separated list of runs");
  PrintFakeOption(o, "tmp",           "Copy code to temporary directory");
  PrintFakeOption(o, "debug",         "Create debugging symbols");
  PrintFakeOption(o, "include=DIRECTORY", "Append dir to macro/header path");
  if (!r) { 
    o << "\n"
      << "Additional options may be defined by the specific train setup\n"
      << std::endl;
    return;
  }
  o << "\nTRAIN OPTIONS:" << std::endl;
  r->PrintHelp(o, "--");
  o << "\n"
    << "The options --run and --include may be passed multiple times\n"
    << "Options --debug and --tmp toggle, so pass only once\n"
    << "For any other option, the last value passed is used\n"
    << std::endl;
}    

/** 
 * Run a given analysis train 
 * 
 * @verbatim
Usage: runTrain --class=CLASS --name=NAME [OPTIONS]

OPTIONS:
  --class=CLASS                    Train class
  --name=NAME                      Name of train
  --run=RUN                        Comma separated list of runs
  --tmp                            Copy class script to temporary directory
  --help                           Show this help
  --par                            Use PAR files (PROOF and Grid)
  --mc                             Assume simulation input
  --debug                          Execute in debugger
  --type=AOD|ESD                   Type of train
  --mode=LOCAL|PROOF|GRID          Execution mode
  --oper=TEST|TERMINATE|FULL|INIT  Operation mode
  --date=YYYY-MM-DD HH:MM:SS       Set date string
  --cluster=HOST                   PROOF cluster
  --dataSet=NAME                   Data set (PROOF only)
  --dataDir=DIRECTORY              Data directory
  --pass=NUMBER                    ESD Pass (grid only)
  --verb=NUMBER                    Verbosity
  --root=TAG                       ROOT version (Grid)
  --aliroot=TAG                    AliROOT version (Grid)
  --alien=TAG                      AliEn API version (Grid)
  --overwrite                      Allow overwrite
  --include=DIRECTORY              Append dir to macro/header path
  @endverbatim
 * 
 * Individual classes derived from TrainSetup may define additional
 * (or less) options.  Pass the option @b --help to see the list of
 * available options
 *
 * Usually, the train defines the execution mode (ESD or AOD) and the
 * option @b --type should not be given
 *
 * If a date is set (option @b --date) either to a string of the form 
 * @c YYYY-MM-DD HH:MM:SS or to the string @c now, then the date will be 
 * appended to the job name. 
 * 
 * A sub-directory named according to the name passed will be made,
 * and relevant files will copied there.
 * 
 * @param argc Count of command line arguments, including program name
 * @param argv Vector of command line arguments, including program name
 *  
 * @return 0 on success, 1 otherwise 
 * 
 * @ingroup pwglf_forward_trains_run
 */
int main(int argc, char** argv)
{
  TString trainOpts;
  TString trainClass;
  TString trainName; 
  TString trainRuns; 
  Int_t   trainEvents = -1;
  Bool_t  trainTmp = false;
  Bool_t  trainDbg = false;
  Bool_t  progHelp = false;

  for (int i = 1; i < argc; i++) { 
    // std::cout << "Arg # " << i << ": \"" << argv[i] << "\"" << std::endl;
    if (argv[i][0] == '-') { 
      if (argv[i][1] == '-') { // Long options 
	TString arg(argv[i]);
	arg.ReplaceAll("\"'", "");
	Int_t   eq = arg.Index("=");
	TString val;
	if (eq != kNPOS) val = arg(eq+1,arg.Length()-eq-1);
	if      (arg.BeginsWith("--class"))   trainClass  = val;
	else if (arg.BeginsWith("--name"))    trainName   = val;
	else if (arg.BeginsWith("--events"))  trainEvents = val.Atoi();
	else if (arg.BeginsWith("--run"))     AppendTo(trainRuns, val);
	else if (arg.BeginsWith("--tmp"))     trainTmp = !trainTmp;
	else if (arg.BeginsWith("--debug"))   trainDbg = !trainDbg;
	else if (arg.BeginsWith("--include")) AppendPath(val);
	else if (arg.BeginsWith("--help"))    progHelp=true; 
	else    AppendTo(trainOpts, arg(2,arg.Length()-2)); 
      }
      else {
	switch (argv[i][1]) { 
	case 'h': progHelp    = true; break;
	case 'c': trainClass  = argv[++i]; break;
	case 'n': trainName   = argv[++i]; break;
	case 'N': trainEvents = strtol(argv[++i],  NULL, 0); break;
	case 't': trainTmp    = !trainTmp; break;
	case 'd': trainDbg    = !trainDbg; break;
	case 'r': AppendTo(trainRuns, argv[++i]); break;
	case 'I': AppendPath(argv[++i]); break;
	default: 
	  Error("runMain", "Unknown option %s", argv[i]);
	  return 1;
	} // switch 
      } // Long options 
    } // Options 
  }

  if (trainClass.IsNull() || trainName.IsNull()) { 
    if (progHelp) {
      Usage(argv[0], std::cout, 0);
      return 0;
    }
    Error(argv[0], "No train-setup class name or train name given");
    return 1;
  }
  if (trainDbg) AppendTo(trainOpts, "debug");
  if (progHelp) AppendTo(trainOpts, "help");

  // Build the code 
  if (!BuildTrain(trainClass, trainTmp, trainDbg)) return 1;
  
  // Make our object 
  TrainSetup* trainSetup = MakeTrain(trainClass, trainName);
  if (!trainSetup) return 1;

  TrainSetup::Runner r(*trainSetup);
  if (!r.Init(trainOpts)) return 1;
  if (r.IsHelpAsked()) { 
    Usage(argv[0], std::cout, &r);
    return 0;
  }

  return r.Run(trainRuns, trainEvents, true) ? 0 : 1;
}
//
// EOF
//

  
  
  
