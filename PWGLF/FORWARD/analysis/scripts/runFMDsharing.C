void runFMDsharing() 
{
  // Mandatory fields for the collector analysis
  TString runMode   = "full";
  TString anaSource = "";
  TString anaType   = "sharing";
  TString dataDir   = "/alice/sim/LHC10e13/";
  TString anaName   = "FMDCollector900GeVPythia";
  TString colSys    = "p-p";
  Float_t cmsNNGeV  = 900;
  Float_t bkG       = 5;

  
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
	    anaName, 
	    colSys, 
	    cmsNNGeV, 
	    bkG);
}
//
// EOF
//
