/**
 * @file testCosmicsTrigger.C
 * @brief Test macro for the cosmics trigger
 *
 * Usage:
 * <pre>
 *   aliroot -b -q testCosmicsTrigger.C'("./","HLTesdTree","local://$ALICE_ROOT/OCDB", nOfEvents)' | tee cosmics-trigger.log
 *   aliroot -b -q testCosmicsTrigger.C'("./","esdTree","local://$ALICE_ROOT/OCDB", nOfEvents)' | tee cosmics-trigger.log
 *   aliroot -b -q testCosmicsTrigger.C'("alien:///alice/data/...","HLTesdTree","local://$ALICE_ROOT/OCDB", nOfEvents)' | tee cosmics-trigger.log
 * </pre>
 *
 * The macro takes as input an ESD file, as the component works on ESD tracks. 
 * The user can speficy which tree will be used. 
 *
 * @author Kalliopi Kanaki (Kalliopi.Kanaki@ift.uib.no)
 */

void testCosmicsTrigger(const char *esdfilename,
		        const char *treename,
			const char *cdbURI,
		        int nofEvents=-1
		        )
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
    cerr << "cannot open file " << strfile << endl;
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
  AliCDBManager *man = AliCDBManager::Instance();
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
  AliHLTSystem *pHLT = AliHLTPluginBase::GetInstance();
  
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libTPCcalib");  
  
  TString arguments, triggerinput;

  // ESD publisher component
  arguments=" -datapath "; arguments+=esdfilename;
  arguments+=" -entrytype "; arguments+=strtree;
  triggerinput="ESD-publisher";

  AliHLTConfiguration esdpublisher(triggerinput.Data(), "ESDMCEventPublisher", "", arguments.Data());

  AliHLTConfiguration cosmicstr("cosmics-trigger","CosmicsTrigger", triggerinput.Data(),"");
  //AliHLTConfiguration globaltriggerconf("global-trigger", "HLTGlobalTrigger", "cosmics-trigger" , "");
    
  // set option for the HLT system
  // arguments
  //  - libraries to be used as plugins
  //  - loglevel=0x79 : Important, Warning, Error, Fatal
  pHLT->ScanOptions("libAliHLTUtil.so libAliHLTMUON.so libAliHLTTRD.so libAliHLTGlobal.so libAliHLTTrigger.so loglevel=0x79");

  pHLT->BuildTaskList("cosmics-trigger");
  //pHLT->BuildTaskList("global-trigger");
  pHLT->Run(nofEvents);
}

void testCosmicsTrigger(){
   cout << "=====================" << endl;
   cout << "  " << endl;
   cout << "Usage: aliroot -b -q testCosmicsTrigger.C\'(\"./\",\"HLTesdTree\",\"local://$ALICE_ROOT/OCDB\", nOfEvents)\' | tee cosmics-trigger.log" << endl;
   cout << "  " << endl;
   cout << "=====================" << endl;


}
