/*
  To compile, do 

  rootcint -f OptionDict.C -c Option.C 
  g++ `root-config --cflags --libs` \
    -lVMC -lGeom -lMinuit -lXMLIO -lTree -lTreePlayer \
    -I$ALICE_ROOT/include -L$ALICE_ROOT/lib/tgt_${ALICE_TARGET} \
    -lSTEERBase -lESD -lAOD -lANALYSIS -lOADB -lANALYSISalice \
    trainMain.cxx -o runTrain
*/
 
#include "AvailableSoftware.C" 
#include "ChainBuilder.C" 
#include "ParUtilities.C"
#include "OutputUtilities.C" 
#include "Option.C"
#include "Helper.C" 
#include "LocalHelper.C" 
#include "ProofHelper.C" 
#include "LiteHelper.C"
#include "AAFHelper.C" 
#include "PluginHelper.C"
#include "AAFPluginHelper.C" 
#include "GridHelper.C"
#include "TrainSetup.C"

#include <TGApplication.h>
#include <TRint.h>
#include <TROOT.h>
#include <TList.h>
#include <TObjString.h>
#include <TString.h>

#include <iostream>
#include <iomanip>

/** 
 * Custom timer to do a deferred start after the application 
 * has been started 
 */
struct Deferred : public TTimer
{
  Deferred(const TString& name, const TString& cls, 
	   const TCollection* opts, bool spawn)
    : TTimer(1000, false), 
      fName(name), 
      fClass(cls), 
      fOptions(opts),
      fSpawn(spawn)
  {
    Start(1000, true);
  }
  Deferred(const Deferred& o)
    : TTimer(), 
      fName(o.fName),
      fClass(o.fClass),
      fOptions(o.fOptions),
      fSpawn(o.fSpawn)
  {}
  Deferred& operator=(const Deferred& o) 
  { 
    if (&o == this) return *this;
    fName    = o.fName;
    fClass   = o.fClass;
    fOptions = o.fOptions;
    fSpawn   = o.fSpawn;
    return *this; 
  }
  Bool_t Notify()
  {
    // gSystem->RemoveTimer(this);
    //Info("Notify","Will run train setup: %s (%s)",fName.Data(),fClass.Data());
    return TrainSetup::Main(fName, fClass, fOptions, true, fSpawn);
  }
  TString fName;
  TString fClass;
  const TCollection* fOptions;
  Bool_t  fSpawn;
};

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
  o << "  --" << std::left << std::setw(30) << opt << " " << desc 
    << std::right << std::endl;
}

void ProgramUsage(const char* progname, std::ostream& o)
{
  o << "Usage: " << progname 
    << " --class=CLASS --name=NAME --url=URI [OPTIONS]\n\n"
    << progname << " specific options:\n";
  PrintFakeOption(o, "class=CLASS",       "Train class");
  PrintFakeOption(o, "name=NAME",         "Name of train");
  PrintFakeOption(o, "include=DIRECTORY", "Append dir to macro/header path");
  PrintFakeOption(o, "batch",             "Batch mode");
  PrintFakeOption(o, "spawn",             "Spawn interactive ROOT shell");
}
/** 
 * Print usage information 
 * 
 * @param progname Program name 
 * @param o        Output stream
 * 
 * @ingroup pwglf_forward_trains_run
 */
void Usage(const char* progname, std::ostream& o)
{
  ProgramUsage(progname, o);
  PrintFakeOption(o, "url=URI",           "Execution URI");
  o << "\nAlternatively to using --url=URI, one can use\n";
  PrintFakeOption(o, "where=BASE_URI",    "Set protocol, user, host, "
		  "and port URI");
  PrintFakeOption(o, "file=FILE_OR_PATH", "File/path part of URI");
  PrintFakeOption(o, "options=OPTIONS",   "Query options for URI");
  PrintFakeOption(o, "anchor=ANCHOR",     "Query anchor for URI");
}

int
main(int argc, char** argv)
{
  TList optList;
  TString name;
  TString cls;
  TString where;
  TString file;
  TString opts;
  TString anchor;
  Bool_t  batch   = false;
  Bool_t  help    = false;
  Bool_t  urlSeen = false;
  Bool_t  spawn   = false;
  // --- Parse options -----------------------------------------------
  for (int i = 1; i < argc; i++) { 
    if (argv[i][0] == '-' && argv[i][1] == '-') { 
      TString arg(argv[i]);
      TString val("");
      arg.ReplaceAll("\"'", "");
      Int_t   eq = arg.Index("=");
      if (eq != kNPOS) val = arg(eq+1, arg.Length()-eq-1);
      if      (arg.BeginsWith("--class"))   cls  = val;
      else if (arg.BeginsWith("--name"))    name = val;
      else if (arg.BeginsWith("--include")) AppendPath(val);
      else if (arg.BeginsWith("--batch"))   batch  = true;
      else if (arg.BeginsWith("--help"))    help   = true;
      else if (arg.BeginsWith("--where"))   where  = val;
      else if (arg.BeginsWith("--file"))    file   = val;
      else if (arg.BeginsWith("--opts"))    opts   = val;
      else if (arg.BeginsWith("--anchor"))  anchor = val;
      else if (arg.BeginsWith("--spawn"))   spawn  = true;
      else {
	if (arg.BeginsWith("--url")) urlSeen = true;
	optList.Add(new TObjString(&(argv[i][2])));
      }
    }
  }
  // --- Initial check or URI/WHERE ----------------------------------
  if (!where.IsNull()) {
    if (urlSeen) {
      Error("main", "option --url and --where mutually exclusive");
      return 1;
    }
    TUrl u(where);
    if (!file.IsNull())   u.SetFile(file);
    if (!opts.IsNull())   u.SetOptions(opts);
    if (!anchor.IsNull()) u.SetAnchor(anchor);
    optList.Add(new TObjString(Form("url=%s", u.GetUrl())));
  }

  // --- check for help ----------------------------------------------
  if (help) {
    if (cls.IsNull()) {
      Usage(argv[0], std::cout);
      return 0;
    }
    ProgramUsage(argv[0], std::cout);
    optList.Add(new TObjString("help"));
  }

  // --- Check name and class ----------------------------------------
  if (name.IsNull()) { 
    Error("main", "No name specified");
    return 1;
  }
  if (cls.IsNull()) { 
    Error("main", "No class specified");
    return 1;
  }

  // --- Setup script path -------------------------------------------
  const char* aliPath  = gSystem->ExpandPathName("$ALICE_ROOT");
  const char* fwdPath  = gSystem->ExpandPathName("$ALICE_ROOT/PWGLF/FORWARD/");
  AppendPath(aliPath);
  AppendPath(Form("%s/include",          aliPath));
  AppendPath(Form("%s/trains",           fwdPath));
  AppendPath(Form("%s/analysis2",        fwdPath));
  AppendPath(Form("%s/analysis2/trains", fwdPath));

  // --- Set-up Application ------------------------------------------
  TApplication* app = 0;
  gROOT->SetBatch(batch);
  if (spawn) {
    // Info("main", "Creating interpreter application");
    TRint* rint = new TRint("runTrain", 0, 0, 0, 0, true);
    rint->SetPrompt("runTrain[%3d]> ");
    app = rint;
  }
  else if (!batch) {
    // Info("main", "Creating GUI application");
    app = new TGApplication("runTrain", 0, 0);
  }
  if (app && !batch) app->InitializeGraphics();
  
  // --- run, possibly in a timer ------------------------------------
  Bool_t ret = true;
  if (!app) 
    ret = TrainSetup::Main(name, cls, &optList);
  else {
    // optList.ls();
    new Deferred(name, cls, &optList, spawn);
    gApplication->Run();
  }

  // --- Return ------------------------------------------------------
  return ret ? 0 : 1;
}
//
// EOF
//

	
    
