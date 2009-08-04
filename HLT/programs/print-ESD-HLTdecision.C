// $Id$
/**
 * @file print-ESD-HLTdecision.C
 * @brief Print HLT decisions per event of ESD
 *
 * <pre>
 * Usage: aliroot -b -q print-ESD-HLTdecision.C
 * </pre>
 *
 * @author Matthias.Richter@ift.uib.no
 * @ingroup alihlt_programs
 */
void print_ESD_HLTdecision(const char* esdFileName="AliESDs.root",
			   int minEvent=0, int maxEvent=-1)
{
  TFile* esdFile=NULL;
  if (esdFileName && esdFileName[0]!=0) esdFile=new TFile(esdFileName);
  if (!esdFileName || esdFileName[0]==0 || esdFile->IsZombie()) {
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
    cerr << "       maxEvent         last event (optional)" << endl;
    cerr << "===============================================================" << endl;
    return -1;
  }
  TTree* pTree=NULL;
  esdFile->GetObject("HLTesdTree", pTree);
  if (!pTree) {
    cout << "no HLT ESD tree found, trying ESD tree" << endl;
    esdFile->GetObject("esdTree", pTree);
  }
  if (!pTree) {
    cerr << "can not find (HLT) ESD tree" << endl;
    return -1;
  }
  
  AliESDEvent* esd=new AliESDEvent;
  esd->CreateStdContent();
  esd->ReadFromTree(pTree);
  for (int event=minEvent; 
       event<pTree->GetEntries() && (maxEvent<minEvent || event<=maxEvent);
       event++) {
    pTree->GetEvent(event);
    cout << "=====================  event " << event << "  ==========================" << endl;
    TObject* decision=esd->GetHLTTriggerDecision();
    if (decision) decision->Print();
    else cout << "   no HLT decision found" << endl;
  }
}
