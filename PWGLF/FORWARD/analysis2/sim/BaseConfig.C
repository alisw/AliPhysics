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

// 
// EOF
//

