// $Id$
/**
 * Helper macro to format a block of AliHLTComponentStatistics entries.
 * The block is usually created by attaching a file writer to a chain,
 * writing only COMPSTAT:PRIV data blocks.
 *
 * The macro translates the block into HLTruns a stand-alone chain
 * Usage:
 * <pre>
 *   aliroot -b -q format-statistics.C | tee format-statistics.log
 * </pre>
 *
 *
 * @ingroup alihlt_benchmark
 * @author Matthias.Richter@ift.uib.no
 */
void format_statistics(const char* infile, const char* outfile="HLT.statistics.root")
{
  AliHLTSystem gHLT;
  gHLT.LoadComponentLibraries("libAliHLTUtil.so");
  AliHLTConfiguration publisher("publisher", "FilePublisher", NULL, "-datatype 'COMPSTAT' 'PRIV' -datafile HLT.statistics.raw");

  AliHLTConfiguration sink1("sink1", "StatisticsCollector", "publisher", "-file out.root -publish 0 -arraysize 200000");

  gHLT.BuildTaskList("sink1");
  gHLT.Run();
}

void format_statistics()
{
  cerr << "===============================================================" << endl;
  cerr << "usage:" << endl;
  cerr << "  aliroot -b -q -l format-statistics.C'(\"infile\", \"outfile\")'" << endl << endl;
  cerr << "please provide input, e.g. \"HLT.statistics.raw\"" << endl;
  cerr << "default output file is \"HLT.statistics.root\"" << endl;
  cerr << "===============================================================" << endl;
}
