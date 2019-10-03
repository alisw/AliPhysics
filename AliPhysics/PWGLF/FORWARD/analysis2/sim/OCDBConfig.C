/**
 * @file   OCDBConfig.C
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Wed Oct 15 13:19:34 2014
 * 
 * @brief  Particular setup of specific storages
 * 
 * This is used by Simulate.C, Reconstruct.C
 */
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
    AddStore("TPC/Calib/RecoParam",	false);
    
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
// 
// EOF
// 
