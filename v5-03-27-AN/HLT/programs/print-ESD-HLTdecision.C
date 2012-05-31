// $Id$
/**
 * @file print-ESD-HLTdecision.C
 * @brief Print HLT decisions per event of ESD
 *
 * <pre>
 * Usage: aliroot -b -q print-ESD-HLTdecision.C
 * </pre>
 *
 * The input file can be a file on Grid like e.g.
 * "alien:///alice/data/2009/LHC09d/000104321/ESDs/pass5/09000104321018.30/AliESDs.root"
 * In that case you need a valid token in order to connect to the Grid.
 * Use 'alien-token-init' from your alien installation.
 *
 * @author Matthias.Richter@ift.uib.no
 * @ingroup alihlt_programs
 */
int print_ESD_HLTdecision(const char* esdFileName="AliESDs.root",
			   int minEvent=0, int maxEvent=-1)
{

  TString strfile=esdFileName;
  if (strfile.Contains("://") && !strfile.Contains("local://")) {
    TGrid::Connect("alien");
  }

  TFile* esdFile=NULL;
  if (esdFileName && esdFileName[0]!=0) esdFile=TFile::Open(esdFileName);
  if (!esdFileName || esdFileName[0]==0 || !esdFile || esdFile->IsZombie()) {
    if (esdFileName && esdFileName[0]!=0)
      cerr << "can not open esd file " << esdFileName << endl;
    else
      cerr << "print-ESD-HLTdecision.C Print HLT decisions per event of ESD" << endl;
    cerr << "===============================================================" << endl;
    cerr << "usage: aliroot -b -q -l print-ESD-HLTdecision.C'(\"AliESDs.root\" " << endl;
    cerr << "                                                 minEvent, maxEvent)'" << endl << endl;
    cerr << "  Parameter:" << endl;
    cerr << "       esdFileName      default \"AliESDs.root\"" << endl;
    cerr << "       minEvent         first event (optional)" << endl;
    cerr << "       maxEvent         last event (optional)" << endl << endl;
    cerr << "  The input file can be a file on Grid like e.g." << endl;
    cerr << "  \"alien:///alice/data/2009/LHC09d/000104321/ESDs/pass5/09000104321018.30/AliESDs.root\"" << endl;
    cerr << "===============================================================" << endl;
    return -1;
  }
  TTree* pTree=NULL;
  esdFile->GetObject("esdTree", pTree);
  if (!pTree) {
    cerr << "can not find ESD tree" << endl;
    return -1;
  }

  AliESDEvent* esd=new AliESDEvent;
  esd->CreateStdContent();
  esd->ReadFromTree(pTree);

  TTree* pHLTTree=NULL;
  esdFile->GetObject("HLTesdTree", pHLTTree);
  if (!pHLTTree) {
    cerr << "can not find HLT ESD tree" << endl;
    return -1;
  }

  if (pTree->GetEntries() != pHLTTree->GetEntries()) {
    cerr << "entries differ: ESD tree " << pTree->GetEntries() << "    HLT ESD tree " << pHLTTree->GetEntries() << endl;
  }
  
  AliESDEvent* HLTesd=new AliESDEvent;
  HLTesd->CreateStdContent();
  HLTesd->ReadFromTree(pHLTTree);

  for (int event=minEvent; 
       event<pTree->GetEntries() && (maxEvent<minEvent || event<=maxEvent);
       event++) {
    pTree->GetEvent(event);
    pHLTTree->GetEvent(event);
    cout << "=====================  event " << event << "  ==========================" << endl;
    cout << "\t ESD: # in file: " << esd->GetEventNumberInFile() << "   time stamp (UCT): " << esd->GetTimeStamp() << endl;
    cout << "\t orbit no (offline/hlt): " <<  esd->GetOrbitNumber() << "/" << HLTesd->GetOrbitNumber() << endl;
    cout << "\t bunch crossing (offline/hlt): " << esd->GetBunchCrossNumber() << "/" << HLTesd->GetBunchCrossNumber() << endl;
    cout << "\t tracks (offline/hlt):\t " << esd->GetNumberOfTracks() << "/"<< HLTesd->GetNumberOfTracks() << endl;
    cout << "\t hltESD content:" << endl;
    cout << "\t    fired triggers:\t" << HLTesd->GetFiredTriggerClasses() << endl;
    cout << "\t    trigger mask:\t 0x" << hex << HLTesd->GetTriggerMask() << dec << endl;
    cout << "\t    time stamp (UCT):\t " << HLTesd->GetTimeStamp() << endl;
    TObject* decision=HLTesd->GetHLTTriggerDecision();
    if (decision) decision->Print();
    else cout << "   no HLT decision found" << endl;
  }
  return 0;
}
