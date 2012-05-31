// $Id$
/**
 * @file sampleAliHLTOUTHandlerEsdBranch.C
 * @brief Example macro to illustrate function of AliHLTOUTHandlerEsdBranch.
 * @note  Macro requires simulated raw data, please run e.g. ppbench setup.
 *
 * HLTOUT handlers inherit from class AliHLTOUTHandler and are utilized to
 * process data blocks in the HLTOUT. As a specific case HLTOUT can contain
 * data blocks/streamed objects designated for storage in custom ESD branches.
 * AliHLTOUTHandlerEsdBranch provides an easy way to implement merging of
 * such blocks/objects into the hltEsd.
 *
 * Usage: Simulate a data sample (e.g. ppbench) and run this macro in the
 * folder for raw data reconstruction. 
 * <pre>
 * aliroot -b -q -l sampleAliHLTOUTHandlerEsdBranch.C \
 *            | tee sampleAliHLTOUTHandlerEsdBranch.log
 * </pre>
 *
 * The macro contains 3 parts:
 * - the hltEsd layout object is adjusted in order to allow for the
 *   new branch
 * - an object is generated, streamed and saved to file
 * - a small HLT chain is run
 * 
 * The HLT chain is as simple as just publishing the binary object file with
 * the appropriate data type. The published object becomes part of the HLTOUT
 * and the handler of the AliHLTAgentSample is called as it will announce
 * the ability to process data blocks of this specific type from HLTOUT.
 *
 * @author Matthias.Richter@ift.uib.no
 * @ingroup alihlt_tutorial
 */
void sampleAliHLTOUTHandlerEsdBranch()
{
  // require a simulated data sample, e.g. the ppbench setup
  // we need raw data in order to run a custom HLT chain
  TString datainput;
  if (!gSystem->AccessPathName("raw.root")) datainput="raw.root";
  else if (!gSystem->AccessPathName("raw0")) datainput="./";

  TString grpstorage;
  if (!gSystem->AccessPathName("GRP/GRP/Data")) grpstorage="local://$PWD";
  else if (!gSystem->AccessPathName("../GRP/GRP/Data")) grpstorage="local://$PWD/..";
  
  if (datainput.IsNull() || grpstorage.IsNull()) {
    cout << "sampleAliHLTOUTHandlerEsdBranch.C: illustrate AliHLTOUTHandlerEsdBranch" << endl;
    cout << " Usage: aliroot -b -q -l sampleAliHLTOUTHandlerEsdBranch.C" << endl;
    cout << "" << endl;
    cout << " Macro requieres a simulated data sample, please run e.g. ppbench setup" << endl;
    cout << " and run in the folder with the raw data" << endl;
    if (datainput.IsNull()) cout << "Error: no raw data input found" << endl;
    if (grpstorage.IsNull()) cout << "Error: no GRP entry found in neither ./ nor ../" << endl;
    return;
  }

  if (!gSystem->AccessPathName("galice.root")) {
    cout << "AliReconstruction from raw data requires to delete the galice.root file" << endl;
    cout << "of previous reconstruction cycles. IMPORTANT: do not delete the file in" << endl;
    cout << "the simulation folder but run in a subfolder like e.g. recraw" << endl;
    return;
  }

  // Basic settings
  TString command;

  TString objectFileName=gSystem->TempDirectory();
  objectFileName+="/";objectFileName+="inputobject.dat";

  // Create the input object
  TClonesArray* inputObject=new TClonesArray("TNamed");
  inputObject->SetName("MyPrivateBranch");

  // Init OCDB
  AliCDBManager * man = AliCDBManager::Instance();
  man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  // set specific storage for GRP entry
  man->SetSpecificStorage("GRP/GRP/Data", grpstorage);

  // Adjust hltEsd layout
  //
  // in order to allow custom branches for the HLT ESD a layout object is
  // stored in OCDB. This object will be generated automatically in the future
  // according to the online setup. It can be adapted manually to enable specific
  // branches. In this example we add the input object (still an empty TClonesArray)
  // to the layout and store the object.
  const char* esdLayoutEntry="HLT/ConfigHLT/esdLayout";
  if (gSystem->AccessPathName(esdLayoutEntry)) {
    command="mkdir -p "; command+=esdLayoutEntry;
    gSystem->Exec(command);
  }
  command="cp $ALICE_ROOT/OCDB/"; command+=esdLayoutEntry; command+="/Run* "; command+=esdLayoutEntry;
  gSystem->Exec(command);
  man->SetSpecificStorage(esdLayoutEntry, "local://$PWD");
  man->SetRun(0);

  AliCDBEntry* esdLayoutObject=man->Get(esdLayoutEntry);
  AliESDEvent* esdLayout=(AliESDEvent*)esdLayoutObject->GetObject();
  if (!esdLayout->FindListObject(inputObject->GetName())) {
    esdLayout->AddObject(inputObject);
  }
  man->Put(esdLayout, esdLayoutObject->GetId(), esdLayoutObject->GetMetaData());
  man->UnloadFromCache(esdLayoutEntry);

  // now add some data to the input object, stream it and save to file
  new ((*inputObject)[inputObject->GetEntriesFast()]) TNamed("data1","some data");
  new ((*inputObject)[inputObject->GetEntriesFast()]) TNamed("data2","some more data");
  AliHLTMessage* pMsg=AliHLTMessage::Stream(inputObject);
  if (!pMsg) {
    cerr << "failed to stream input object" << endl;
    return;
  }
  ofstream objectFile(objectFileName);
  objectFile.write(pMsg->Buffer(), pMsg->Length());
  objectFile.close();
  delete pMsg;
  if (gSystem->AccessPathName(objectFileName)) {
    cerr << "failed to write input file " << objectFileName << endl;
    return;
  }

  //////////////////////////////////////////////////////////////////////////////////////
  //
  // Reconstruction settings
  AliReconstruction rec;

  rec.SetRunReconstruction("HLT ITS TPC");
  rec.SetInput(datainput);
  rec.SetWriteESDfriend(kFALSE);

  // QA options
  rec.SetRunQA(":") ;
  //rec.SetQARefDefaultStorage("local://$ALICE_ROOT/QAref") ;

  //////////////////////////////////////////////////////////////////////////////////////
  //
  // setup the HLT system
  AliHLTSystem* pHLT=AliHLTPluginBase::GetInstance();

  // define the HLT chain to be run in AliReconstruction
  // arguments:
  //  1) id of the configuartion, later used to refer to this configuration
  //  2) id of the component to run
  //  3) parents, no parents in this case
  //  4) optional component arguments
  TString arguments;
  arguments.Form("-datatype ROOTTOBJ SMPL -datafile %s", objectFileName.Data());
  AliHLTConfiguration publisher("Object-Publisher", "FilePublisher", "", arguments.Data());

  // set option for the HLT module in AliReconstruction
  // arguments
  //  - ignore-hltout : ignore the HLTOUT payload from the HLT DDLs
  //  - libraries to be used as plugins
  //  - loglevel=0x7c : Info, Important, Warning, Error, Fatal
  //  - chains=Object-Publisher : chains to be run
  rec.SetOption("HLT",
		"ignore-hltout " 
		"libAliHLTSample.so "
		"loglevel=0x7c "
		"chains=Object-Publisher "
		);

  rec.SetRunPlaneEff(kFALSE);

  // switch off cleanESD
  rec.SetCleanESD(kFALSE);

  AliLog::Flush();
  rec.Run();

}
