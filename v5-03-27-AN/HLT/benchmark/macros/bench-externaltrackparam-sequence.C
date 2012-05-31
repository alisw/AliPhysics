// $Id$
/*
 * Benchmark macro for transportation of AliExternalTrackParam arrays.
 *
 * Usage:
 * <pre>
 *   aliroot -b -q bench-externaltrackparam.C | tee bench-externaltrackparam.log
 * </pre>
 *
 *
 * @ingroup alihlt_benchmark
 * @author Matthias.Richter@ift.uib.no
 */
void bench_externaltrackparam_sequence(int events=100)
{
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // init the HLT system in order to define the analysis chain below
  //
  AliHLTSystem* gHLT=AliHLTPluginBase::GetInstance();
  gHLT->SetGlobalLoggingLevel(0x7c);

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // define the analysis chain to be run
  //
  int verbosity=0;
  int levels=11;
  TString lastPublisher;
  TString arg;
  for (int pub=0; pub<levels; pub++) {
    TString publisher;
    // publishers of AliExternalTrackParam arrays
    publisher.Form("PUB_%02d", pub);
    arg="-minsize 9000 -maxsize 10000";
    //arg="-minsize 1 -maxsize 2";
    if (pub<10) {
      // first 10 ones are TClonesArrays with the corresponding compression
      arg+=" -tclonesarray";
      //arg+=" -tobjarray";
      arg+=" -object-compression="; arg+=pub;
    } else {
      // last one is a c-array without compression
      arg+=" -carray";
    }

    arg+=" -rangeoffset -1000 -rangemodulo "; arg+=events/10;
    arg+=" -verbosity "; arg+=verbosity;
    AliHLTConfiguration publisherconf(publisher.Data(), "BenchmarkAliExternalTrackParam", lastPublisher.Data(), arg.Data());
    lastPublisher=publisher;
  }

  arg=" -verbosity "; arg+=verbosity;
  AliHLTConfiguration dumpconf("sink1", "BenchmarkAliExternalTrackParam"   , lastPublisher.Data(), "-verbosity 1");

  AliHLTConfiguration statconf("stat", "StatisticsCollector", /*dumpInput.Data()*/"sink1", "");

  AliHLTConfiguration writer("statwriter", "ROOTFileWriter"   , "stat", "-datafile HLT.statistics.root -concatenate-events -overwrite");

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // Init and run the chain
  //
  gHLT->LoadComponentLibraries("libAliHLTBenchmark.so libAliHLTUtil.so");
  gHLT->BuildTaskList("statwriter");
  gHLT->Run(events);
}
