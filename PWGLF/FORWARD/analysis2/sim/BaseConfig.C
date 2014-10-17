/**
 * @file   BaseConfig.C
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Wed Oct 15 12:52:58 2014
 * 
 * @brief  Base classes for configurations shared amoung steps. 
 * 
 * 
 */
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
  
};
/** Global variable */
VirtualDetCfg* detCfg = 0;

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

struct VirtualEGCfg 
{
  TString runType;
  VirtualEGCfg() : runType("") {}
  virtual ~VirtualEGCfg() {}
  /**
   * Load the general libraries needed
   *
   */
  static void LoadGen(const TString& runType) {
    if (!gROOT->GetClass("AliStructFuncType"))
      gSystem->Load("liblhapdf");      // Parton density functions
    if (!gROOT->GetClass("TPythia6"))
      gSystem->Load("libEGPythia6");   // TGenerator interface
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
    gSystem->Load(Form("libpythia%s", vers));
    gSystem->Load(Form("libAliPythia%c", m));
  }
  /**
   * Load HIJING libraries
   */
  static void LoadHijing()
  {
    LoadPythia();
    if (gROOT->GetClass("THijing")) return;
    gSystem->Load("libhijing");
    gSystem->Load("libTHijing");
    AliPDG::AddParticlesToPdgDataBase();
  }
  /**
   * Load HydJet libraries
   */
  static void LoadHydjet()
  {
    if (gROOT->GetClass("TUHKMgen")) return;
    gSystem->Load("libTUHKMgen");
  }
  /**
   * Load DPMJet libraries
   */
  static void LoadDpmjet()
  {
    LoadPythia();
    if (gROOT->GetClass("TDPMjet")) return;
    gSystem->Load("libdpmjet");
    gSystem->Load("libTDPMjet");
  }
  /**
   * Load AMPT libraries
   */
  static void LoadAmpt()
  {
    LoadPythia();
    if (gROOT->GetClass("TAmpt")) return;
    gSystem->Load("libampt");
    gSystem->Load("libTAmpt");
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
    if (g && smear) g->SetVertexSmear(AliGenerator::kPerEvent);

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

