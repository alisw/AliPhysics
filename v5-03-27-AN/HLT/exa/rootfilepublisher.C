/**
 * @file rootfilepublisher.C
 * @brief Macro for testing ROOT-file publishing and writing
 *
 * This macro is a testing/example macro of how to use the RootFilePublisher 
 * (AliHLTRootFilePublisher) and RootFileWriter (AliHLTRootFileWriter). 
 * It defines only two component in the chain, the publisher, which 
 * publishes the content of the root files according to the selection. 
 * Be aware there can be several root objects in one root file.
 *
 * The file myRootFile.root has to be changed with the actual file.
 * In this example, an ESD file is used, which contains a "esdTree" and
 * a "HLTesdTree" object. Only the "HLTesdTree" object is selected. The 
 * data blocks, which will are published here, have the data type of 
 * {ROOTTOBJ,"***"} and specification of 0x00000000. They also have to 
 * been adopted to the actual content of the root file. 
 *
 * @author thaeder@kip.uni-heidelberg.de
 * @ingroup alihlt_tutorial
 */

/** RootFilePublisher test macro
 *  @param nEvents Number of events which should be processed
 */
void rootfilepublisher(Int_t nEvents=1) {

  TString writerInput;
  TString arg;

  AliHLTSystem gHLT;
  gHLT.LoadComponentLibraries("libAliHLTUtil.so");

  // -- Root publisher configuration
  // !!! myRootFile.root has to be exchanged with an existing one.
  arg.Form("-objectname HLTesdTree -datatype 'ESD_TREE' 'TPC '-dataspec 0x00000000 -datafile myRootFile.root");

  // -- The AliHLTRootFilePublisher (Id 'ROOTFilePublisher') is a data source. 
  //    It provides the given files to the subsequent components in the chain.
  //    see AliHLTRootFilePublisher for more options
  AliHLTConfiguration RootPublisher("RootPublisher", "ROOTFilePublisher", NULL, arg.Data() );
  if (!writerInput.IsNull()) writerInput+=" ";
  writerInput+="RootPublisher";

  // -
  // -- Processing Components can be put in here
  // - 

  // -- The AliHLTRootFileWriter (Id 'ROOTFileWriter') is a data sink. It writes
  // all incoming data blocks to files. Several options available.
  AliHLTConfiguration RootWriter("RootWriter", "ROOTFileWriter", writerInput.Data(),"-datafile event");

  // -- Here you specify the top most configuration of the chain. The
  // configuration depends on all the parents. The task lisy is build
  // according to that.
  gHLT.BuildTaskList("RootWriter");
  gHLT.Run(nEvents);
}
