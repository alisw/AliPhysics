/**
 * @file   RunTrain.C
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Tue Oct 16 17:49:49 2012
 * 
 * @brief  Script to run a train 
 *
 * @ingroup pwglf_forward_trains
 */
/** 
 * Build a script 
 * 
 * @param name     Name of script 
 * @param verbose  Whether to be verbose 
 * @param force    Whether to force re-compilation 
 * @param debug    Whether to enable debug symbols  
 * 
 * @return true on success 
 *
 * @ingroup pwglf_forward_trains
 */
Bool_t
BuildScript(const char* name, Bool_t verbose, Bool_t force, Bool_t debug)
{
  const char opt[] = { (force ? '+' : debug ? 'g' : '\0'), 
		       (force && debug ? 'g' : '\0'), '\0' };
  if (verbose) Printf("Building %s ...",name);
  Int_t error;
  Int_t ret = gROOT->LoadMacro(Form("%s.C+%s", name, opt), &error);
  if (ret < 0 || error) { 
    Error("BuildScript", "Failed to build %s: %d", error);
    return false;
  }
  return true;
  
}

/** 
 * Build all helper classes 
 * 
 * @param verbose Whether to be verbose 
 * @param force   Whether to force re-builds
 * @param debug   Whether to enable debug symbols 
 * @param all     Whether to build all helpers 
 *
 * @return true on success 
 *
 * @ingroup pwglf_forward_trains
 */
Bool_t
BuildHelpers(Bool_t verbose, Bool_t force, Bool_t debug,
	     Bool_t all=false)
{
  gSystem->AddIncludePath("-I$ALICE_ROOT/include");
  gSystem->AddIncludePath("-I$ALICE_ROOT/PWGLF/FORWARD/trains/");
  gROOT->SetMacroPath(Form("%s:$ALICE_ROOT/PWGLF/FORWARD/trains",
			   gROOT->GetMacroPath()));
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  const char* scripts[] = { "AvailableSoftware", 
			    "ChainBuilder", 
			    "ParUtilities",
			    "OutputUtilities", 
			    "Option",
			    "Helper", 
			    "TrainSetup",
			    (all ? "LocalHelper" : 0), 
			    "ProofHelper", 
			    "LiteHelper", 
			    "AAFHelper", 
			    "PluginHelper", 
			    "AAFPluginHelper", 
			    "GridHelper", 
			    0 };
  const char** ptr = scripts;
  while ((*ptr)) {
    if (!BuildScript(*ptr, verbose, force, debug)) return false;
    ptr++;
  }
  return true;
}

/** 
 * Show the usage 
 * 
 *
 * @ingroup pwglf_forward_trains
 */
void PlainUsage()
{
  std::cout << "Usage: .x RunTrain.C(NAME,CLASS,OPTIONS)\n\n"
	    << "  NAME       Name of train (free form)\n"
	    << "  CLASS      Name of class implementing TrainSetup\n"
	    << "  OPTIONS    Comma separated list of options\n"
	    << std::endl;
}

/** 
 * Function to run a train.  
 * 
 * @param name  Name of the train. 
 * @param cls   class name of train setup
 * @param uri   Exection URI  
 * @param opts  Optons 
 * 
 * @return true on success
 *
 * @ingroup pwglf_forward_trains
 */
Bool_t RunTrain(const TString& name, const TString& cls, 
		const TUrl& uri,     const TString& opts)
{
  // Check for help 
  if (name.IsNull() || name.EqualTo("help", TString::kIgnoreCase) || 
      cls.IsNull()  || cls.EqualTo("help", TString::kIgnoreCase) || 
      !uri.IsValid()) {
    PlainUsage();
    return true;
  }
  
  Bool_t verb = opts.Contains("verbose");
  // Build our helpers 
  if (!BuildHelpers(verb, false, true)) return false;

  // Tokenize options 
  if (!opts.EndsWith(",")) opts.Append(",");
  opts.Append("url="); 
  opts.Append(uri.GetUrl());
  TObjArray* optList = opts.Tokenize(",");
  return TrainSetup::Main(name, cls, optList, false);
}
/*
 * EOF
 */	
    
