/**
 * @file   BaseConfig.C
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Wed Oct 15 12:52:58 2014
 * 
 * @brief  Base classes for configurations shared amoung steps. 
 * 
 * 
 */
//====================================================================
/** 
 * Base class for detector configuration. By default, everything is on
 * except ACORDE.
 */
struct VirtualDetCfg 
{
  virtual Bool_t UseABSO()   const { return true;  }
  virtual Bool_t UseACORDE() const { return false; }
  virtual Bool_t UseDIPO()   const { return true;  }
  virtual Bool_t UseEMCAL()  const { return true;  }
  virtual Bool_t UseFMD()    const { return true;  }
  virtual Bool_t UseFRAME()  const { return true;  }
  virtual Bool_t UseHALL()   const { return true;  }
  virtual Bool_t UseITS()    const { return true;  }
  virtual Bool_t UseMAG()    const { return true;  }
  virtual Bool_t UseMUON()   const { return true;  }
  virtual Bool_t UsePHOS()   const { return true;  }
  virtual Bool_t UsePIPE()   const { return true;  }
  virtual Bool_t UsePMD()    const { return true;  }
  virtual Bool_t UseHMPID()  const { return true;  }
  virtual Bool_t UseSHIL()   const { return true;  }
  virtual Bool_t UseT0()     const { return true;  }
  virtual Bool_t UseTOF()    const { return true;  }
  virtual Bool_t UseTPC()    const { return true;  }
  virtual Bool_t UseTRD()    const { return true;  }
  virtual Bool_t UseVZERO()  const { return true;  }
  virtual Bool_t UseZDC()    const { return true;  }
  virtual void Print() 
  {
    Printf("ABSO:   %3s", UseABSO()    ? "yes" : "no");
    Printf("ACORDE: %3s", UseACORDE()  ? "yes" : "no");
    Printf("DIPO:   %3s", UseDIPO()    ? "yes" : "no");
    Printf("EMCAL:  %3s", UseEMCAL()   ? "yes" : "no");
    Printf("FMD:    %3s", UseFMD()     ? "yes" : "no");
    Printf("FRAME:  %3s", UseFRAME()   ? "yes" : "no");
    Printf("HALL:   %3s", UseHALL()    ? "yes" : "no");
    Printf("ITS:    %3s", UseITS()     ? "yes" : "no");
    Printf("MAG:    %3s", UseMAG()     ? "yes" : "no");
    Printf("MUON:   %3s", UseMUON()    ? "yes" : "no");
    Printf("PHOS:   %3s", UsePHOS()    ? "yes" : "no");
    Printf("PIPE:   %3s", UsePIPE()    ? "yes" : "no");
    Printf("PMD:    %3s", UsePMD()     ? "yes" : "no");
    Printf("HMPID:  %3s", UseHMPID()   ? "yes" : "no");
    Printf("SHIL:   %3s", UseSHIL()    ? "yes" : "no");
    Printf("T0:     %3s", UseT0()      ? "yes" : "no");
    Printf("TOF:    %3s", UseTOF()     ? "yes" : "no");
    Printf("TPC:    %3s", UseTPC()     ? "yes" : "no");
    Printf("TRD:    %3s", UseTRD()     ? "yes" : "no");
    Printf("VZERO:  %3s", UseVZERO()   ? "yes" : "no");
    Printf("ZDC:    %3s", UseZDC()     ? "yes" : "no");
  }
  /** 
   * Get the string of enabled detectors for local reconstruction.
   * 
   * @param enable On return, contains string of enable detectors
   */
  void GetRecoString(TString& enable) const
  {
    if (UseITS())	Append2Str(enable, "ITS"); 
    if (UseTPC())	Append2Str(enable, "TPC"); 
    if (UseTRD())	Append2Str(enable, "TRD"); 
    if (UseTOF())	Append2Str(enable, "TOF"); 
    if (UsePHOS())	Append2Str(enable, "PHOS"); 
    if (UseHMPID())	Append2Str(enable, "HMPID"); 
    if (UseEMCAL())	Append2Str(enable, "EMCAL"); 
    if (UseMUON())	Append2Str(enable, "MUON"); 
    if (UseFMD())	Append2Str(enable, "FMD"); 
    if (UseZDC())	Append2Str(enable, "ZDC"); 
    if (UsePMD())	Append2Str(enable, "PMD"); 
    if (UseT0())	Append2Str(enable, "T0"); 
    if (UseVZERO())	Append2Str(enable, "VZERO");
  }
  /** 
   * Get the string of detectors for which we should make Summable
   * Digits
   * 
   * @param sDigits On returm contains the string of enable detectors
   */
  void GetSDigitString(TString& sDigits) const 
  {
    if (UseTRD())	Append2Str(sDigits, "TRD"); 
    if (UseTOF())	Append2Str(sDigits, "TOF"); 
    if (UsePHOS())	Append2Str(sDigits, "PHOS"); 
    if (UseHMPID())	Append2Str(sDigits, "HMPID"); 
    if (UseEMCAL())	Append2Str(sDigits, "EMCAL"); 
    if (UseMUON())	Append2Str(sDigits, "MUON"); 
    if (UseFMD())	Append2Str(sDigits, "FMD"); 
    if (UseZDC())	Append2Str(sDigits, "ZDC"); 
    if (UsePMD())	Append2Str(sDigits, "PMD"); 
    if (UseT0())	Append2Str(sDigits, "T0"); 
    if (UseVZERO())	Append2Str(sDigits, "VZERO");
  }
  /** 
   * Get the sting of detectors for which we should do the hit to
   * digit conversion directly.
   * 
   * @param fromHits On returm contains the string of enable detectors
   */
  void GetHits2DigitsString(TString& fromHits) const
  {
    if (UseITS())	Append2Str(fromHits, "ITS");
    if (UseTPC())	Append2Str(fromHits, "TPC");
  }
  /** 
   * Append a C style string to a string, possibly adding a space before
   * 
   * @param str     Where to append
   * @param append  What to append
   */
  static void Append2Str(TString& str, const char* append)
  {
    if (!str.IsNull()) str.Append(" ");
    str.Append(append);
  }
  /** 
   * Overload this member function to return the string "*OCDB*" to
   * load geometry from OCDB.  If this returns a null pointer, the
   * geometry will be constructed.  Alternatively, the one can specify
   * a file nanme (possibly on remote storage).
   * 
   * @return File name, "*OCDB*" string, or null 
   */
  virtual const char* GeometrySource() const
  {
    return 0;
  }
};
/** Global variable */
VirtualDetCfg* detCfg = 0;

//====================================================================
/**
 * Base class for the OCDG configration 
 */
struct VirtualOCDBCfg
{
  /** 
   * This member function must return the default prefix. 
   * 
   * @return Prefix of OCDB specific storages
   */
  virtual const char* Prefix() const { return ""; }
  /** 
   * This member function should define the real setup. 
   * 
   * @param forSim Whether we're setting up for simulations or not 
   */
  virtual void Init(Bool_t forSim) 
  {
    ::Fatal("VirtualOCDBConfig", "Dummy init called - redefine!");
  }
  /** 
   * Set the specific storage for a given key (possibly wild-carded). 
   * 
   * @param key    Key 
   * @param ideal  Whether it is residual or ideal
   */
  void AddStore(const char*    key, 
		Bool_t         ideal)
  {
    AliCDBManager* cdb = AliCDBManager::Instance();
    const char* prefix = Prefix();
    TString     path   = Form("alien://Folder=/alice/simulation/%s/%s",
			      prefix, !ideal ? "Residual" : "Ideal");
    ::Info("AddStore", "%s -> %s", key, path.Data());
    cdb->SetSpecificStorage(key, path);
  }
};

/** Global variable */
VirtualOCDBCfg* ocdbCfg = 0;

//====================================================================
/** 
 * Event generator configuration 
 * 
 */
struct VirtualEGCfg 
{
  TString runType;
  VirtualEGCfg() : runType("") {}
  virtual ~VirtualEGCfg() {}
  virtual Bool_t IsLego() const { return false; }
    /** 
   * Set the default generator based on the beam type 
   *
   * - p-p PYTHIA
   * - p-A or A-p DPMJet
   * - A-A Hijing 
   */
  static const char* DeduceRunType()
  {
    if      (grp->IsPP())                return "pythia";
    else if (grp->IsPA() || grp->IsAP()) return "dpmjet";
    else if (grp->IsAA())                return "hijing";
    return "hijing";
  }

  static void LoadLibrary(const TString& name,
			  const TString& cls="")
  {
    // If we're looking for a specific class, check that first, and
    // if available, do nothing;
    if (!cls.IsNull() && gROOT->GetClass(cls)) return;

    // Now check the list of loaded and linekd libraries 
    TString libs(gSystem->GetLibraries("", "SD"));

    // IF already in the list, do nothing 
    if (libs.Contains(name)) return;

    // Otherwise load the library 
    gSystem->Load(name);
  }
  /**
   * Load the general libraries needed
   *
   */
  static void LoadGen(const TString& runType) {
    LoadLibrary("liblhapdf","AliStructFuncType"); // Parton density functions
    LoadLibrary("libEGPythia6","TPythia6");       // TGenerator interface
    if (!runType.EqualTo("hydjet", TString::kIgnoreCase))
      LoadPythia(false);
  }

  /**
   * Load the pythia libraries
   *
   * @param vers Optional version post-fix
   */
  static void LoadPythia(Bool_t gen=true, const char* vers="6.4.21")
  {
    if (gen) LoadGen("");
    char m = vers[0];
    if (gROOT->GetClass(Form("AliPythia6%c", m))) return;

    TString v(vers); 
    v.ReplaceAll(".", "_");
    LoadLibrary("libmicrocern");
    LoadLibrary(Form("libpythia%s",v.Data()));
    LoadLibrary(Form("libAliPythia%c", m));
  }
  /**
   * Load HIJING libraries
   */
  static void LoadHijing()
  {
    LoadPythia();
    if (gROOT->GetClass("THijing")) return;
    LoadLibrary("libhijing");
    LoadLibrary("libTHijing");
    AliPDG::AddParticlesToPdgDataBase();
  }
  /**
   * Load HydJet libraries
   */
  static void LoadHydjet()
  {
    LoadLibrary("libTUHKMgen","TUHKMgen");
  }
  static void LoadEposLHC()
  {
    LoadLibrary("libEPOSLHC","TEposLHC");
  }
  static void LoadEpos()
  {
    LoadLibrary("libEPOS","TEpos");
  }
  static void LoadTherminator()
  {
    LoadLibrary("libTTherminator","TTherminator");
  }
  /**
   * Load DPMJet libraries
   */
  static void LoadDpmjet()
  {
    LoadPythia();
    if (gROOT->GetClass("TDPMjet")) return;
    LoadLibrary("libdpmjet");
    LoadLibrary("libTDPMjet");
  }
  /**
   * Load AMPT libraries
   */
  static void LoadAmpt()
  {
    LoadPythia();
    if (gROOT->GetClass("TAmpt")) return;
    LoadLibrary("libampt");
    LoadLibrary("libTAmpt");
  }
  /**
   * Make the generator
   *
   * @param rt    Event generator identifier 
   * @param b1    Least impact parameter 
   * @param b2    Largest impact parameter 
   * @param smear If true, smear interaction per event 
   *
   * @return Point to newly allocated generator or null
   * 
   */
  AliGenerator* MakeGenerator(const TString& rt, 
			      Float_t b1, 
			      Float_t b2, 
			      Bool_t smear=true)
  {
    if (rt.IsNull()) { 
      ::Fatal("MakeGenerator", "No EG spec given");
      return 0;
    }

    TString runType = rt;
    runType.ToLower();

    AliGenerator* g = CreateGenerator(runType,b1,b2);
    if (g && smear) {
      // Origin and sigmaZ will be supplied from OCDB.
      // We need to set the XY smearing though
      // Numbers from Ruben 
      g->SetSigma(0.0025, 0.0029, 0);
      g->SetVertexSmear(AliGenerator::kPerEvent);
    }

    return g;
  }
  /**
   * Make our decayer
   *
   * @param rt The EG to use 
   *
   * @return Newly allocated decayer or null
   */
  TVirtualMCDecayer* MakeDecayer(const TString& rt)
  {
    if (rt.IsNull()) { 
      ::Fatal("MakeGenerator", "No EG spec given");
      return 0;
    }

    TString runType = rt;
    rt.ToLower();

    TVirtualMCDecayer* decayer = CreateDecayer(runType);

    if (decayer) decayer->Init();
    return decayer;
  }
protected:
  /** 
   * Create the generator.  This function must be defined in a derived class. 
   * 
   * @param runType The generator ID (all lower case)
   * @param b1      Least impact parameter 
   * @param b2      Largest impact parameter 
   * 
   * @return Must return a pointer to a new AliGenerator or null
   */
  virtual AliGenerator* CreateGenerator(const TString& runType, 
					Float_t b1, 
					Float_t b2) = 0;
  /** 
   * Create the decayer.  This function must be defined in a derived class. 
   * 
   * @param runType The generator ID (all lower case)
   * 
   * @return Must return a pointer to a new TVirtualMCDecayer or null
   */
  virtual TVirtualMCDecayer* CreateDecayer(const TString& runType) = 0;

};
/** Global variable */
VirtualEGCfg* egCfg = 0;

//====================================================================
/**
 * Base class for trains 
 * 
 */
struct VirtualTrain 
{


  /** 
   * Run this train 
   * 
   * @param run 
   * @param xmlFile 
   * @param stage 
   * @param cdb 
   * 
   * @return 
   */
  Bool_t Run(UInt_t      run, 
	     const char* xmlFile = "wn.xml", 
	     Int_t       stage   = 0, 
	     const char* cdb     = "raw://")
  {
    // --- Load configuration script ---------------------------------
    LoadConfig();
    
    // --- Set-up for CDB access through Grid ------------------------
    TString cdbString(cdb);
    if (cdbString.Contains("raw://")) {
      TGrid::Connect("alien://");
      if (!gGrid || !gGrid->IsConnected()) {
	::Error("Run", "No grid connection");
	return false;
      }  
    }
    
    // --- Some environment variables --------------------------------
    // Temp dir is here, and compilation is here too 
    gSystem->Setenv("TMPDIR", gSystem->pwd());
    gSystem->SetBuildDir(gSystem->pwd(), kTRUE);

    // --- Now load common libraries ---------------------------------
    LoadBaseLibraries();
    
    // --- Now create and configure manager --------------------------
    AliAnalysisManager *mgr  = new AliAnalysisManager(GetName(), 
						      "Production train");
    mgr->SetRunFromPath(grp->run);
    
    // --- Create ESD input handler ------------------------------------
    AliESDInputHandlerRP *esdHandler = new AliESDInputHandlerRP();
    mgr->SetInputEventHandler(esdHandler);
    if (UseFriends()) {
      esdHandler->SetReadFriends(kTRUE);
      esdHandler->SetActiveBranches("ESDfriend");
    }

    // --- Monte Carlo handler -----------------------------------------
    if (UseMC()) {
      AliMCEventHandler* mcHandler = new AliMCEventHandler();
      mgr->SetMCtruthEventHandler(mcHandler);
      mcHandler->SetPreReadMode(1);
      mcHandler->SetReadTR(true);
    }
    // --- AOD output handler ----------------------------------------
    if (MakeAOD()) {
      AliAODHandler* aodHandler   = new AliAODHandler();
      aodHandler->SetOutputFileName("AliAOD.root");
      mgr->SetOutputEventHandler(aodHandler);
    }
    
    // --- Call user routine for adding tasks ------------------------
    if (!AddTasks()) return false;
    
    // --- Check if we are to merge ----------------------------------
    if (stage > 0) 
      return Merge(xmlfile, stage);

    // --- Otherwise run the train -----------------------------------
    TChain* chain = CreateChain();
    if (!chain) return false;

    TStopwatch timer;
    timer.Start();
    if (!mgr->InitAnalysis()) {
      ::Error("Run", "Failed to initialize the train");
      return false;
    }

    mgr->PrintStatus();
    mgr->SetSkipTerminate(kTRUE);
    mgr->StartAnalysis("local", chain);
    timer.Print();
    
  }
  /** 
   * Merge requested files 
   * 
   * @param dir    Output directory 
   * @param stage  Stage 
   * 
   * @return true on success 
   */
  Bool_t Merge(const char* dir, Int_t stage)
  {

    TStopwatch    timer;     
    timer.Start();
    TString       outputDir     = dir;
    Bool_t        final         = outputDir.Contains("Stage");
    TCollection*  outputFiles   = GetFilesToMerge(stage, final);
    if (!outputFiles) { 
      ::Warning("Merge", "Nothing to merge");
      return true;
    }
    TIter       iter(outputFiles);
    TObjString* str           = 0;
    Bool_t      merged        = kTRUE;
    while((str = static_cast<TObjString*>(iter()))) {
      TString& outputFile = str->GetString();
      // Skip already merged outputs
      if (!gSystem->AccessPathName(outputFile)) {
	::Warning("Merge","Output file <%s> found. Not merging again.",
		  outputFile.Data());
	continue;
      }
      merged = AliAnalysisAlien::MergeOutput(outputFile, 
					     outputDir, 
					     10, 
					     stage);
      if (merged) continue; 
      
      ::Error("Merge", "Cannot merge %s\n", outputFile.Data());
    }
    // --- possible merge file information files ---------------------
    if (MergeFileInfo()) { 
      TString infolog = "fileinfo.log";
      AliAnalysisAlien::MergeInfo(infolog, dir); 
    }

    // --- If not final stage, get out here --------------------------
    if (!final) { 
      ValidateOutput();
      timer.Print();
      return true;
    }
    
    // --- set up and run termiante ----------------------------------
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    mgr->SetSkipTerminate(kFALSE);
    if (!mgr->InitAnalysis()) {
      ::Error("Merge", "Failed to initialize the train");
      return false;
    }
    
    mgr->PrintStatus();
    mgr->StartAnalysis("gridterminate", (TTree*)0);
    ValidateOutput();
    timer.Print();

    return true;
  }
	       
  /** 
   * Load a library/module 
   * 
   * @param module Library/module name 
   * 
   * @return true on success
   */
  Bool_t LoadLibrary(const char *module)
  {
    // Load a module library in a given mode. Reports success.
    Int_t result = 0;
    TString mod(module);
    ::Info("LoadLibrary", "Loading %s", module);
    gROOT->IncreaseDirLevel();

    if (mod.IsNull()) {
      ::Error("LoadLibrary", "Empty module name");
      gROOT->DecreaseDirLevel();
      return kFALSE;
    }

    // If a library is specified, just load it
    if (mod.EndsWith(".so")) {
      mod.Remove(mod.Index(".so"));
      ::Info("LoadLibrary", "Loading .so: %s", mod.Data()); 
      result = gSystem->Load(mod);
      if (result < 0) {
	::Error("oadLibrary", "Could not load library %s", module);
      }
      gROOT->DecreaseDirLevel();      
      return (result >= 0);
    }
    // Check if the library is already loaded
    if (strlen(gSystem->GetLibraries(module, "", kFALSE)) > 0) {
      ::Info("LoadLibrary", "Module %s.so already loaded", module);
      gROOT->DecreaseDirLevel();      
      return kTRUE;
    }

    ::Info("LoadLibrary", "Trying to load lib%s.so", module);
    result = gSystem->Load(Form("lib%s", module));
    if (result < 0)
      ::Error("LoadLibrary", "Could not load module %s", module);
    ::Info("LoadLibrary", "Module %s, successfully loaded", module);
    gROOT->DecreaseDirLevel();      
    return (result >= 0);
  }
  /** 
   * Load common libraries 
   * 
   * @return true on sucess 
   */
  virtual Bool_t LoadBaseLibraries()
  {
    // Load common analysis libraries.
    if (!gSystem->Getenv("ALICE_PHYSICS")) {
      ::Error("LoadBaseLibraries", 
	      "Analysis trains requires that analysis libraries are "
	      "compiled with a local AliRoot");
      return false;
    }

    Bool_t success = true;
    // Load framework classes. Par option ignored here.
    success &= LoadLibrary("libSTEERBase.so");
    success &= LoadLibrary("libESD.so");
    success &= LoadLibrary("libAOD.so");
    success &= LoadLibrary("libANALYSIS.so");
    success &= LoadLibrary("libOADB.so");
    success &= LoadLibrary("libANALYSISalice.so");
    success &= LoadLibrary("libESDfilter.so");
    success &= LoadLibrary("libCORRFW.so");
    success &= LoadLibrary("libPWGPP.so");
    gROOT->ProcessLine(".include $ALICE_PHYSICS/include");
    if (success) {
      ::Info("LoadBaseLibraries", 
	     "Load common libraries:    SUCCESS");
      ::Info("LoadBaseLibraries", 
	     "Include path for Aclic compilation:\n%s",
	     gSystem->GetIncludePath());
    } else {
      ::Info("LoadBaseLibraries", 
	     "Load common libraries:    FAILED");
    }
    return success;
  }
  /** 
   * Create the input chain
   * 
   * @return Pointer to newly allocated train 
   */
  TChain* CreateChain()
  {
    if (gSystem->AccessPathName("AliESDs.root")) {
      ::Error("CreateChain", 
	      "File: AliESDs.root not in ./data dir");
      return 0;
    }
    
    // Create the input chain
    TChain* chain = new TChain("esdTree");
    chain->Add("AliESDs.root");
    if (!chain->GetNtrees()) {
      delete chain;
      chain = 0;
    }

    return chain;
  }
  /** 
   * Helper function to make @c outputs_valid file 
   * 
   */
  void ValidateOutput()
  {
    std::ofstream out;
    out.open("outputs_valid", ios::out);
    out.close();    
  }  

  /** 
   * @{ 
   * @name Functions to overload 
   */
  /** 
   * Load the configuration script. Override to load specific script.
   */
  virtual void LoadConfig() {};
  /** 
   * Override to set a name of the analysis manager 
   * 
   * @return Name of analysis manager 
   */
  virtual const char* GetName() const { return "dummy"; }
  /** 
   * Override to return true if friends are needed. 
   * 
   * @return false
   */
  virtual Bool_t UseFriends() const { return false; }
  /** 
   * Override to return true if MC info is needed
   * 
   * @return false
   */
  virtual Bool_t UseMC() const { return false; }
  /** 
   * Override to return true if AODs should be made 
   * 
   * @return false
   */
  virtual Bool_t MakeAOD() const { return false; }
  /**
   * User rountine for adding tasks. Override to add tasks to the
   * train.
   *
   * @return true
   */
  virtual Bool_t AddTasks() const { return true; }
  /** 
   * Override to return true to merge file information files. 
   * 
   * @return false
   */
  virtual Bool_t MergeFileInfo() const { return false; }
  /** 
   * Return the list of ouput files (TObjString objects)
   *
   * @param stage Merge stage 
   * @param final Final merging (also terminate)
   *
   * @return Pointer to TCollection. 
   */
  virtual TCollection* GetFilesToMerge(Int_t stage, Bool_t final) const 
  { 
    return 0; 
  }
};




//====================================================================
/**
 * A function so that we can do TROOT::Macro.  Does nothing but print a message.
 *
 */
void BaseConfig()
{
  Info("", "Defined base classes for configuration");
}
//
// EOF
//

