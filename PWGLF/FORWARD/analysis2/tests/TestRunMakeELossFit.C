/** 
 * Run the energy loss fit finder and generate corrections output file 
 * 
 * @param sys       Collision system 
 * @param cms       Center of mass energy per nucleon in GeV
 * @param field     Magnetic field 
 * @param mc        Whether this is for Monte-Carlo data
 * @param filename  Input file name 
 *
 * @ingroup pwg2_forward_scripts_tests
 *
 * @depcrecated 
 * The class AliFMDELossFitter automatically generates the
 * AliFMDCorrELossFit object.
 *
 */
void
TestRunMakeELossFit(UShort_t    sys, 
		    UShort_t    cms, 
		    Short_t     field, 
		    Bool_t      mc=false,
		    const char* filename="forward_eloss.root")
{
  std::cout << "Loading libraries ..." << std::endl;
  gROOT->Macro("$ALICE_ROOT/PWG2/FORWARD/analysis2/scripts/LoadLibs.C");

  std::cout << "Loading compile script ..." << std::endl;
  gROOT->LoadMacro("$ALICE_ROOT/PWG2/FORWARD/analysis2/scripts/Compile.C");
 
  std::cout << "Compiling MakeELossFit.C script ..." << std::endl;
  Compile("$ALICE_ROOT/PWG2/FORWARD/analysis2/scripts/TestMakeELossFit.C"); 

  std::cout << "Making MakeELossFit object (sys=" << sys 
	    << ", cms=" << cms << ", field=" << field << ", mc=" << mc 
	    << ")" << std::endl;
  MakeELossFit mef(sys, cms, field, mc, filename); 

  std::cout << "Running maker ..." << std::endl;
  mef.Run();
}
//
// EOF
//
