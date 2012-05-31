/**
 * @file ESDMCEventPublisher.C
 * @brief Macro for testing AliESDEvent and AliMCEvent publishing and writing
 *
 * This macro is a testing/example macro of how to use the ESDMCEventPublisher 
 * (AliHLTESDMCEventPublisherComponent) and RootFileWriter (AliHLTRootFileWriter). 
 * It defines only two component in the chain, the publisher, which 
 * publishes the content of the root files according to the selection. 
 * Be aware there can be several root objects in one root file.
 *
 * The file datapath has to be changed with an actual one. It must contain 
 * the files: <br>
 *   - AliESDs.root<br>
 *   - Kinematics.root<br>
 *   - galice.root<br>
 *   - TrackRefs.root<br>
 *
 * Entry type can be one, all or some of :<br>
 *   - ESD<br>
 *   - HLTESD<br>
 *   - MC<br>
 *
 * For more descriptions, especially the used datatypes and specification:
 * @see AliHLTESDMCEventPublisherComponent
 *
 * @author thaeder@kip.uni-heidelberg.de
 * @ingroup alihlt_tutorial
 */

/** ESDMCEventPublisher test macro
 *  @param nEvents Number of events which should be processed
 */
void ESDMCEventPublisher(Int_t nEvents=1) {

  TString writerInput;
  TString arg;

  AliHLTSystem gHLT;
  gHLT.LoadComponentLibraries("libAliHLTUtil.so");

  // -- Root publisher configuration
  // !!! myDataPath has to be exchanged with an existing one.
  arg.Form("-entrytype ESD -entrytype HLTESD -entrytype MC -dataspec 0x0000001F -datapath /opt/HLT/analysis/HLT-HEAD_2008-09-24/exa -datapath /home/jthaeder/jet/data/test");
  //  arg.Form("-entrytype ESD -entrytype HLTESD -entrytype MC -dataspec 0x0000001F -datapath mydatapath0 -datapath mydatapath1");

  // -- The AliHLTESDMCEventPublisher (Id 'ESDMCEventPublisher') is a data source. 
  //    It provides AliESDEvents and AliMCEvents out of the given datapaths to 
  //    the subsequent components in the chain.
  //    see AliHLTESDMCEventPublisherComponent for more options 
  AliHLTConfiguration ESDMCEventPublisher("ESDMCEventPublisher", "ESDMCEventPublisher", NULL, arg.Data() );
  if (!writerInput.IsNull()) writerInput+=" ";
  writerInput+="ESDMCEventPublisher";
  
  // -
  // -- Processing Components can be put in here
  // - 

  // -- The AliHLTRootFileWriter (Id 'ROOTFileWriter') is a data sink. It writes
  // all incoming data blocks to files. Several options available.
  AliHLTConfiguration rootWriter("RootWriter", "ROOTFileWriter", writerInput.Data(),"-datafile event");

  // -- Here you specify the top most configuration of the chain. The
  // configuration depends on all the parents. The task lisy is build
  // according to that.
  gHLT.BuildTaskList("RootWriter");
  gHLT.Run(nEvents);
}
