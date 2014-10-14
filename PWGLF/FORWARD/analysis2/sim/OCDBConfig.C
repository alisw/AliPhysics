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
VirtualOCDBCfg* ocdbCfg = 0;

/** 
 * Specific implementation.  Note, this requires that GRP.C has been
 * loaded and exectuted before calling Init.
 */
struct OCDBCfg : public VirtualOCDBCfg
{
  const char* Prefix() const { return "2008/v4-15-Release"; }
  void Init(Bool_t forSim)
  {
    Bool_t is10h = grp->period.EqualTo("LHC10h");
    
    // --- ITS  (1 Total) ----------------------------------------------
    AddStore("ITS/Align/Data",		forSim);
    if (!forSim) 
      AddStore("ITS/Align/SPDSparseDead",false);
  
    // --- MUON (1 object) ---------------------------------------------
    AddStore("MUON/Align/Data",		forSim); 

    // ---- TPC (6 total) ----------------------------------------------
    AddStore("TPC/Calib/TimeGain",	forSim);
    AddStore("TPC/Calib/ClusterParam",	forSim);
    AddStore("TPC/Calib/AltroConfig",	forSim);
    AddStore("TPC/Calib/Correction",	forSim);
    AddStore("TPC/Align/Data",		forSim);
    AddStore("TPC/Calib/TimeDrift",	forSim);
    AddStore("TPC/Calib/RecoParam",	(forSim && !is10h));
    
    // --- ZDC for 2010 the following is needed ------------------------
    // (https://savannah.cern.ch/task/?func=detailitem&item_id=33180#comment46)
    if (is10h) AddStore("ZDC/Align/Data",true); 
  }
};


void OCDBConfig()
{
  ::Info("OCDBConfig", "Creating OCDB configuration");
  ocdbCfg = new OCDBCfg;
}



