/**
 * @file   dndeta/tracklets3/Post.C
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Tue Sep 27 08:58:09 2016
 * 
 * @brief  Do the post processing - steering macro 
 * 
 * 
 * @ingroup pwglf_forward_tracklets
 */
/** 
 * @{
 * @name Post-processing steering code 
 *
 * @ingroup pwglf_forward_tracklets
 */ 
/** 
 * Print out usage information 
 * 
 * @param o Output stream
 *
 * @ingroup pwglf_forward_tracklets
 */
void Usage(std::ostream& o)
{
  o << "Usage: Post(SIM,REAL,OUT,PROC,VIZ,N)\n"
    << "\n"
    << "  SIM    Simulation data file \n"
    << "  REAL   Real data file \n"
    << "  OUTPUT (optional) Output name \n"
    << "  PROC   (optional) processing options \n"
    << "  VIZ    (optional) visualisation options \n"
    << "  N      (optional) max number of centrality bins \n"
    << std::endl;
  o << "Processing options:\n"
    << "\n"
    << "  0x0001   Do scaling by unity\n"
    << "  0x0002   Do scaling by full average\n"
    << "  0x0004   Do scaling by eta differential\n"
    << "  0x0008   Do scaling by fully differential\n"
    << "  0x0010   Correct for decay of strange to secondary \n"
    << "  0x0020   Correct for fake tracklets\n"
    << "  0x0040   Correct for non-flat centrality in each bin\n"
    << "  0x1000   MC closure test\n"
    << std::endl;
  o << "Visualization options:\n"
    << "\n"
    << "  0x0001   Draw general information\n"
    << "  0x0002   Draw parameters\n"
    << "  0x0004   Draw weights\n"
    << "  0x0008   Draw dNch/deta\n"
    << "  0x0010   Draw alphas\n"
    << "  0x0020   Draw delta information\n"
    << "  0x0040   Draw backgrounds\n"
    << "  0x0100   Whether to make a PDF\n"
    << "  0x0200   Whether to pause after each plot\n"
    << "  0x0400   Draw in landscape\n"
    << "  0x0800   Alternative markers\n"
    << std::endl;
}

/** 
 * Format input file name and return short name 
 * 
 * @param inp  Input 
 * @param shrt On return, the short name 
 * 
 * @return The input file name 
 *
 * @ingroup pwglf_forward_tracklets
 */
const TString& FormatInput(const char* inp, TString& shrt)
{
  static TString tmp;
  tmp = "";
  Long_t flags;
  if (gSystem->GetPathInfo(inp, 0, (Long_t*)0, &flags, 0) != 0) {
    Warning("FormatInput", "Cannot stat %s", inp);
    return tmp;
  }
  Info("FormatInput", "Input=%s stat=0x%x", inp, flags);
  if (flags & 0x1) {
    shrt = inp;
    tmp.Form("%s/tracklet_dndeta.root", inp);
  }
  else {
    shrt = gSystem->DirName(inp);
    tmp  = inp;
  }
  return tmp;        
}
/** 
 * 
 * 
 * @param sim    Simulation data file 
 * @param real   Real data file 
 * @param output (optional) Output name 
 * @param proc   (optional) processing options 
 * @param viz    (optional) visualisation options 
 * @param n      (optional) max number of centrality bins 
 *
 * Processing options: 
 *
 * - 0x0001   Do scaling by unity
 * - 0x0002   Do scaling by full average
 * - 0x0004   Do scaling by eta differential
 * - 0x0008   Do scaling by fully differential
 * - 0x0010   Correct for decay of strange to secondary 
 * - 0x0020   Correct for fake tracklets
 * - 0x0040   Correct for non-flat centrality in each bin
 * - 0x1000   MC closure test
 *
 * Visualization options: 
 *
 * - 0x0001   Draw general information
 * - 0x0002   Draw parameters
 * - 0x0004   Draw weights
 * - 0x0008   Draw dNch/deta
 * - 0x0010   Draw alphas
 * - 0x0020   Draw delta information
 * - 0x0040   Draw backgrounds
 * - 0x0100   Whether to make a PDF
 * - 0x0200   Whether to pause after each plot
 * - 0x0400   Draw in landscape
 * - 0x0800   Alternative markers
 * 
 * @ingroup pwglf_forward_tracklets
 * @relates AliTrackletdNdeta2
 */
void Post(const char* sim,
	  const char* real,
	  const char* output=0,
	  UInt_t      proc=0x2,
	  UInt_t      viz=0x32f,
	  UInt_t      n=10,
	  UInt_t      sNN=5023)
{
  if (TString(sim) .Contains("help",TString::kIgnoreCase) ||
      TString(real).Contains("help",TString::kIgnoreCase)) {
    Usage(std::cout);
    return;
  }
  // Set the path to the code 
  TString fwd = "$ALICE_PHYSICS/PWGLF/FORWARD/analysis2/dndeta/tracklets3";
  if (gSystem->Getenv("ANA_SRC")) fwd = "$ANA_SRC/dndeta/tracklets3";

  // Compile the needed code 
  gSystem->AddIncludePath(Form("-I%s",fwd.Data()));
  gROOT->LoadMacro(Form("%s/AliTrackletAODUtils.C++g",fwd.Data()));
  gROOT->LoadMacro(Form("%s/AliTrackletdNdeta2.C++g",fwd.Data()));

  // Set inputs and output
  TString realShrt, simShrt;  
  TString realFile = FormatInput(real, realShrt);
  TString simFile  = FormatInput(sim,  simShrt);
  TString outFile(output && output[0] != '\0' ? output :
		  Form("%s_%s",realShrt.Data(),simShrt.Data()));
  if (proc & 0x1) outFile.Append("_unit");
  if (proc & 0x2) outFile.Append("_const");
  if (proc & 0x4) outFile.Append("_eta");
  if (proc & 0x8) outFile.Append("_etaipz");
  gSystem->mkdir(outFile,true);
  Printf("===========================================\n"
	 " Real data file:        %s\n"
	 " Simulation data file:  %s\n"
	 " Output directory:      %s\n"
	 "===========================================",
	 realFile.Data(), simFile.Data(), outFile.Data());

  // Create the object to do the post processing, and run it
  AliTrackletdNdeta2* p = new AliTrackletdNdeta2;  
  p->Run(proc,viz,n,realFile,simFile,outFile);

  // Extract a GSE
  Printf("Extracting GraphSysErr object(s)");
  gROOT->LoadMacro(Form("%s/ExtractGSE2.C",fwd.Data()));
  ExtractGSE2(outFile,sNN);

  Printf("All output stored in %s", outFile.Data());
}
/* @} */
// EOF
