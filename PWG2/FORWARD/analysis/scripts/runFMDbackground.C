void runFMDbackground() {

  // Mandatory fields for the background correction
  TString runMode   = "full";
  TString anaSource = "";
  TString addLibs   = "libPWG2forward.so";
  TString anaType   = "background";
  
  // Optional fields - Must however be asserted a value for the analysis to run 
  TString dataDir  = "/alice/sim/LHC10f6/126437";
  TString anaName  = "Correction";
  TString colSys   = "p-p";
  Float_t cmsNNGeV = 900;
  Float_t bkG      = 5;

  
  /** 
   * Run an FMD corrections job. 
   * 
   * @param runMode     Running mode (full, test, terminate, submit, offline)
   * @param anaType     What to do (background, collector, sharing)
   * @param dataDir     Input data directory 
   * @param anaSource   Analysis source (if any)
   * @param addLibs     Libraries to add 
   * @param anaName     Analysis name 
   * @param colSys      Collision system (p-p, Pb-Pb, A-A)
   * @param cmsNNGeV    Center of mass energy per nucleon in GeV
   * @param bkG         Magnetic field in kilo gaus 
   * @param aliceTag    AliROOT tag to use 
   * @param rootTag     ROOT tag to use 
   * @param apiTag      API tag to use 
   */
  gROOT->LoadMacro("runFMDjob.C");
  runFMDjob(runMode, 
	    anaType, 
	    dataDir, 
	    anaSource, 
	    addLibs, 
	    anaName, 
	    colSys, 
	    cmsNNGeV, 
	    bkG);
}
  
    
