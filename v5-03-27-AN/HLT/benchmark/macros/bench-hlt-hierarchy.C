// $Id$
/**
 * Benchmark macro for transportation of fake data blocks through an HLT
 * hierarchy.
 *
 * Usage:
 * <pre>
 *   aliroot -b -q bench-hlt-hierarchy.C | tee bench-hlt-hierarchy.log
 * </pre>
 * Options:
 * <pre>
 *   aliroot -b -q bench-hlt-hierarchy.C'(2, "root raw")'
 * </pre>
 * First argument allows to specifiy the number of events, the second one
 * is a string defining the output, "root" for a statistics in a root file
 * and "raw" for just the raw data block. Can be combined. Default is "root".
 *
 * @ingroup alihlt_benchmark
 * @author Matthias.Richter@ift.uib.no
 */
void bench_hlt_hierarchy(int events=100, const char* option="root")
{
  TString options=option;
  AliHLTSystem gHLT(0x7c);
  gHLT.LoadComponentLibraries("libAliHLTUtil.so");

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // define the analysis chain to be run
  //
  const int iNofLevel0=6;
  const int iNofLevel1=18;
  const int iNofLevel2=2;
  TString level3Input;
  TString arg, component;
  for (int level2=0; level2<iNofLevel2; level2++) {
    TString level2Input;
    for (int level1=0; level1<iNofLevel1; level1++) {
      TString level1Input;
      for (int level0=0; level0<iNofLevel0; level0++) {
	// level0 publisher component
	arg.Form("-datatype DUMMYDAT FAKE -size 10000 -range 1000 -modulo 10 -decrement 1000");
	component.Form("L0_%d_%02d_%d", level2, level1, level0);
	AliHLTConfiguration l0conf(component.Data(), "DataGenerator", NULL , arg.Data());

	if (level1Input.Length()>0) level1Input+=" ";
	level1Input+=component;
      }

      // level 1 components
      arg.Form("-datatype DUMMYDAT FAKE -multiplier 0.2");
      component.Form("L1_%d_%02d", level2, level1);
      AliHLTConfiguration l1conf(component.Data(), "DataGenerator", level1Input.Data(), arg.Data());
      if (level2Input.Length()>0) level2Input+=" ";
      level2Input+=component;
    }

    // level 2 components
    arg.Form("-datatype DUMMYDAT FAKE -multiplier 1.0");
    component.Form("L2_%d", level2);
    AliHLTConfiguration l2conf(component.Data(),"DataGenerator",level2Input.Data(), arg.Data());
    if (level3Input.Length()>0) level3Input+=" ";
    level3Input+=component;
  }

  arg.Form("-datatype DUMMYDAT FAKE -multiplier 1.0");
  AliHLTConfiguration l3conf("L3","DataGenerator",level3Input.Data(),arg.Data());
  AliHLTConfiguration statroot("statroot", "StatisticsCollector"   , "L3", "-file HLT.statistics.root -publish 0");
  AliHLTConfiguration statraw("statraw", "FileWriter"   , "L3", "-datafile HLT.statistics.raw -concatenate-events -concatenate-blocks -datatype COMPSTAT PRIV");

  if (options.Contains("root"))
      gHLT.BuildTaskList("statroot");
  if (options.Contains("raw"))
      gHLT.BuildTaskList("statraw");
  gHLT.Run(events);
}
