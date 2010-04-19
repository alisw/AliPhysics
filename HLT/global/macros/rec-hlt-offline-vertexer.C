// $Id$
/**
 * @file rec-hlt-offline-vertexer.C
 * @brief Example macro to run the AliHLTGlobalOfflineVertexerComponent in
 * AliReconstruction.
 *
 * The macro supports two running modes for tests of the
 * 'GlobalOfflineVertexer' HLT component:
 * 1) subscribes to the output of the default HLT reconstruction chain,
 * extracts the ESD from the input, and runs the AliVertexerTracks.
 * 2) run on already reconstructed ESD, either "esdTree" or "HLTesdTree"
 *
 * <pre>
 * 1) Usage: aliroot -b -q -l \
 *     rec-hlt-offline-vertexer.C'("file", "cdb", minEvent, maxEvent)'
 * 2) Usage: aliroot -b -q -l \
 *     rec-hlt-offline-vertexer.C'("esdfile", "tree", "cdb", nofEvents)'
 *
 * Examples:
 *     rec-hlt-offline-vertexer.C'("alien:///alice/data/2009/.../....root")' 
 *     rec-hlt-offline-vertexer.C'("raw.root", "local://$ALICE_ROOT/OCDB", minEvent, MaxEvent)'
 *     rec-hlt-offline-vertexer.C'("raw.root", "esdTree", "local://$ALICE_ROOT/OCDB")'
 *
 * Defaults
 *     cdb="raw://"  -> take OCDB from GRID
 *     minEvent=-1   -> no lower event selection
 *     maxEvent=-1   -> no upper event selection
 *     nofEvents=-1  -> all events
 *
 * </pre>
 *
 * The input file can be a file on the grid, indicated by the tag
 * 'alien://'. By default also the OCDB is set to the GRID.
 * If either the file or the OCDB is taken from the GRID, the macro
 * connects to the Grid in the beginning.
 *
 * It is always a good idea to use the OCDB from the
 * Grid as this will contain all the necessary objects and the latest
 * calibration. The special URI 'raw://' is most advisable as it selects
 * the storage automatically from the run number. Other options are e.g.
 * - "alien://folder=/alice/data/2010/OCDB"
 * - "local://$ALICE_ROOT/OCDB"
 *
 * Note: You need a valid GRID token, if you want to access files directly
 * from the Grid, use 'alien-token-init' of your alien installation.
 *
 * @author Matthias.Richter@ift.uib.no
 * @ingroup alihlt_global
 */
void rec_hlt_offline_vertexer(const char *filename,
			      const char *cdbURI,
			      int minEvent=-1,
			      int maxEvent=-1)
{
  // connect to the GRID if we use a file or OCDB from the GRID
  TString struri=cdbURI;
  TString strfile=filename;
  if (struri.BeginsWith("raw://") ||
      strfile.Contains("://") && !strfile.Contains("local://")) {
    TGrid::Connect("alien");
  }

  // Set the CDB storage location
  AliCDBManager * man = AliCDBManager::Instance();
  man->SetDefaultStorage(cdbURI);
  if (struri.BeginsWith("local://")) {
    // set specific storage for GRP entry
    // search in the working directory and one level above, the latter
    // follows the standard simulation setup like e.g. in test/ppbench
    if (!gSystem->AccessPathName("GRP/GRP/Data")) {
      man->SetSpecificStorage("GRP/GRP/Data", "local://$PWD");
    } else if (!gSystem->AccessPathName("../GRP/GRP/Data")) {
      man->SetSpecificStorage("GRP/GRP/Data", "local://$PWD/..");      
    }
  }

  //////////////////////////////////////////////////////////////////////////////////////
  //
  // Reconstruction settings
  AliReconstruction rec;

  if (minEvent>=0 || maxEvent>minEvent) {
    if (minEvent<0) minEvent=0;
    if (maxEvent<minEvent) maxEvent=minEvent;
    rec.SetEventRange(minEvent,maxEvent);
  }
  rec.SetRunReconstruction("HLT ITS TPC");
  rec.SetWriteESDfriend(kFALSE);
  rec.SetInput(filename);

  // QA options
  rec.SetRunQA(":") ;
  //rec.SetQARefDefaultStorage("local://$ALICE_ROOT/QAref") ;

  //////////////////////////////////////////////////////////////////////////////////////
  //
  // setup the HLT system
  AliHLTSystem* pHLT=AliHLTPluginBase::GetInstance();

  // define a configuration for the GlobalOfflineVertexer component
  // a histgram component (GlobalVertexerHisto) creates histograms from
  // the output of the vertexer component. Histograms are saved in file
  // offlvertex.root
  //
  // The same chain is set up for the GlobalVertexer component, output
  // saved in glblvertex.root.
  //
  // arguments of a configuration:
  //  1) id of the configuartion, later used to refer to this configuration
  //  2) id of the component to run
  //  3) parents, here the configuration 'GLOBAL-esd-converter' of the libAliHLTGlobal
  //  4) optional component arguments
  AliHLTConfiguration offlvertexer("Offline-Vertexer", "GlobalOfflineVertexer", "GLOBAL-esd-converter", "");
  AliHLTConfiguration offlhistogram("Offline-Histogram", "GlobalVertexerHisto", "Offline-Vertexer", "");
  AliHLTConfiguration offlwriter("Offline-Writer", "ROOTFileWriter", "Offline-Histogram", "-datafile offlvertex.root -overwrite -concatenate-events");

  AliHLTConfiguration glblvertexer("Global-Vertexer", "GlobalVertexer", "GLOBAL-esd-converter", "");
  AliHLTConfiguration glblhistogram("Global-Histogram", "GlobalVertexerHisto", "Global-Vertexer", "");
  AliHLTConfiguration glblwriter("Global-Writer", "ROOTFileWriter", "Global-Histogram", "-datafile glblvertex.root -overwrite -concatenate-events");

  // set option for the HLT module in AliReconstruction
  // arguments
  //  - ignore-hltout : ignore the HLTOUT payload from the HLT DDLs
  //  - libraries to be used as plugins
  //  - loglevel=0x79 : Important, Warning, Error, Fatal
  //  - chains=Offline-Vertexer : chains to be run
  rec.SetOption("HLT",
		"ignore-hltout " 
		"libAliHLTUtil.so libAliHLTGlobal.so "
		"loglevel=0x79 "
		"chains=Offline-Writer,Global-Writer "
		);

  rec.SetRunPlaneEff(kFALSE);

  // switch off cleanESD
  rec.SetCleanESD(kFALSE);

  AliLog::Flush();
  rec.Run();

}

// run default chain and connect vertexer components to global ESD
// using the raw OCDB
void rec_hlt_offline_vertexer(const char *filename,
		       int minEvent=-1,
		       int maxEvent=-1)
{
  rec_hlt_offline_vertexer(filename, "raw://", minEvent, maxEvent);
}

// run on ESD
void rec_hlt_offline_vertexer(const char *esdfilename,
			      const char *treename,
			      const char *cdbURI,
			      int nofEvents=-1)
{
  // check the name of the tree
  TString strtree=treename;
  if (strtree.CompareTo("esdTree")==0) strtree="ESD";
  else if (strtree.CompareTo("HLTesdTree")==0) strtree="HLTESD";
  else {
    cerr << "invalid treename '" << treename << "', supported 'esdTree' and 'HLTesdTree'" << endl;
    return;
  }

  // connect to the GRID if we use a file or OCDB from the GRID
  TString struri=cdbURI;
  TString strfile=esdfilename;
  if (struri.BeginsWith("raw://") ||
      strfile.Contains("://") && !strfile.Contains("local://")) {
    TGrid::Connect("alien");
  }

  // open the ESD file and get the event count
  if (!strfile.EndsWith("/")) strfile+="/";
  strfile+="AliESDs.root";
  TFile* esdfile=TFile::Open(strfile);
  if (!esdfile || esdfile->IsZombie()) {
    cerr << "can not open file " << strfile << endl;
    return;
  }

  // get number of events
  TTree* pTree=NULL;
  esdfile->GetObject(treename, pTree);
  if (!pTree) {
    cerr << "can not find " << treename << " in file " << strfile << endl;
    return;
  }
  if (pTree->GetEntries()<=0) {
    cerr << "empty tree " << treename << " in file " << strfile << endl;
    return;
  }
  
  AliESDEvent* esd=new AliESDEvent;
  esd->ReadFromTree(pTree);
  pTree->GetEntry(0);

  if (nofEvents<0 || nofEvents>pTree->GetEntries())
    nofEvents=pTree->GetEntries();

  // Set the CDB storage location
  AliCDBManager * man = AliCDBManager::Instance();
  man->SetDefaultStorage(cdbURI);
  man->SetRun(esd->GetRunNumber());
  if (struri.BeginsWith("local://")) {
    // set specific storage for GRP entry
    // search in the working directory and one level above, the latter
    // follows the standard simulation setup like e.g. in test/ppbench
    if (!gSystem->AccessPathName("GRP/GRP/Data")) {
      man->SetSpecificStorage("GRP/GRP/Data", "local://$PWD");
    } else if (!gSystem->AccessPathName("../GRP/GRP/Data")) {
      man->SetSpecificStorage("GRP/GRP/Data", "local://$PWD/..");      
    }
  }

  //////////////////////////////////////////////////////////////////////////////////////
  //
  // setup the HLT system
  AliHLTSystem* pHLT=AliHLTPluginBase::GetInstance();

  TString arguments, vertexerinput;

  // ESD publisher component
  arguments=" -datapath "; arguments+=esdfilename;
  arguments+=" -entrytype "; arguments+=strtree;
  vertexerinput="ESD-publisher";
  AliHLTConfiguration esdpublisher(vertexerinput.Data(), "ESDMCEventPublisher", "", arguments.Data());

  AliHLTConfiguration offlvertexer("Offline-Vertexer", "GlobalOfflineVertexer", vertexerinput.Data(), "");
  AliHLTConfiguration offlhistogram("Offline-Histogram", "GlobalVertexerHisto", "Offline-Vertexer", "");
  AliHLTConfiguration offlwriter("Offline-Writer", "ROOTFileWriter", "Offline-Histogram", "-datafile offlvertex.root -overwrite -concatenate-events");

  AliHLTConfiguration glblvertexer("Global-Vertexer", "GlobalVertexer", vertexerinput.Data(), "");
  AliHLTConfiguration glblhistogram("Global-Histogram", "GlobalVertexerHisto", "Global-Vertexer", "");
  AliHLTConfiguration glblwriter("Global-Writer", "ROOTFileWriter", "Global-Histogram", "-datafile glblvertex.root -overwrite -concatenate-events");

  // set option for the HLT system
  // arguments
  //  - libraries to be used as plugins
  //  - loglevel=0x79 : Important, Warning, Error, Fatal
  pHLT->ScanOptions("libAliHLTUtil.so libAliHLTGlobal.so "
		    "loglevel=0x79 "
		    );

  pHLT->BuildTaskList("Offline-Writer");
  pHLT->BuildTaskList("Global-Writer");
  pHLT->Run(nofEvents);
}

// help
void rec_hlt_offline_vertexer()
{
  cout << "rec_hlt_offline_vertexer.C: The macro supports two running modes of the " << endl; 
  cout << "                            'GlobalOfflineVertexer' HLT component" << endl;
  cout << " Run in AliReconstruction" << endl;
  cout << "   Usage: aliroot -b -q -l \\" << endl;
  cout << "     rec-hlt-offline-vertexer.C'(\"file\", \"cdb\", minEvent, maxEvent)'" << endl;
  cout << "" << endl;
  cout << " Run on already reconstructed ESD, either \"esdTree\" or \"HLTesdTree\"" << endl;
  cout << "   Usage: aliroot -b -q -l \\" << endl;
  cout << "     rec-hlt-offline-vertexer.C'(\"esdfile\", \"tree\", \"cdb\", nofEvents)'" << endl;
  cout << "     three first arguments mandatory" << endl;
  cout << "" << endl;
  cout << " Examples:" << endl;
  cout << "     rec-hlt-offline-vertexer.C'(\"alien:///alice/data/2009/.../....root\")' " << endl;
  cout << "     rec-hlt-offline-vertexer.C'(\"raw.root\", \"local://$ALICE_ROOT/OCDB\", minEvent, MaxEvent)'" << endl;
  cout << "     rec-hlt-offline-vertexer.C'(\"raw.root\", \"esdTree\", \"local://$ALICE_ROOT/OCDB\")'" << endl;
  cout << "" << endl;
  cout << " Defaults" << endl;
  cout << "     cdb=\"raw://\"  -> take OCDB from GRID" << endl;
  cout << "     minEvent=-1   -> no lower event selection" << endl;
  cout << "     maxEvent=-1   -> no upper event selection" << endl;
  cout << "     nofEvents=-1  -> all events" << endl;
}
