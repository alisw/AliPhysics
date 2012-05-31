// $Id$
/**
 * @file conf-sample.C
 * @brief A sample configuration macro for HLT chains in AliRoot.
 *
 * The macro defines a simple HLT analysis chain consisting of multiple
 * data publishers with a processor each (just a dummy component copying
 * the data blocks), and a common data sink.
 *
 * \b Note: The file publisher needs a file to read, either you replace
 * \em some-data.dat with the path of an existing file or just create a
 * dummy file in the current working directory. Futhermore, there has to
 * be at least one simulated event since AliReconstruction relies on a
 * couple of filesin the folder.
 *
 * Usage: from the aliroot prompt
 * <pre>
 {
   AliReconstruction rec;                 // the reconstruction instance
   rec.SetInput("./");                    // to be independent of galice.root
   rec.SetLoadAlignFromCDB(kFALSE);
   rec.SetFillTriggerESD(kFALSE);
   rec.SetRunQA(":");
   rec.SetRunVertexFinder(kFALSE);
   rec.SetRunLocalReconstruction("HLT");  // run local rec only for HLT 
   rec.SetRunTracking("");                // switch off tracking
   rec.SetFillESD("HLT");                 // 
   rec.SetOption("HLT", "libAliHLTSample.so libAliHLTUtil.so "
                        "config=$ALICE_ROOT/HLT/exa/conf-sample.C "
                        "chains=sink");
   //rec.SetEventRange(0,0);
   rec.Run();
 }
 * </pre>
 *
 * @author Matthias.Richter@ift.uib.no
 * @ingroup alihlt_tutorial
 */
{
  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////
  //
  // the configuration
  const int nofPublishers=5;
  TString writerInput;

  for (int i=0; i<nofPublishers; i++) {
    TString publisherName;
    TString processorName;
    TString arg;
    publisherName.Form("fp_%d", i);
    arg.Form("-datatype 'DUMMYDAT' 'SMPL' -datafile some-data.dat "
	     "-dataspec %d", i);
    // The AliHLTFilePublisher (component Id \em 'FilePublisher' provides
    // the given file (see AliHLTFilePublisher for more options) to the
    // subsequent components in the chain.
    AliHLTConfiguration publisher(publisherName.Data(), "FilePublisher",
				  NULL, arg.Data());

    processorName.Form("cp_%d", i);
    // The AliHLTDummyComponent (Id \em 'Dummy') just forwards a certain
    // fraction of the input to the output or just repeats the input data
    // if percentage > 100
    AliHLTConfiguration copy(processorName.Data(), "Dummy", publisherName.Data(),
			     "output_percentage 80");

    // add all processors to the input of the data sink
    if (!writerInput.IsNull()) writerInput+=" ";
    writerInput+=processorName;
  }

  // The AliHLTFileWriter (Id 'FileWriter') is a data sink. It writes
  // all incoming data blocks to files. Several options available.
  AliHLTConfiguration sink("sink", "FileWriter", writerInput.Data(), NULL);
}
