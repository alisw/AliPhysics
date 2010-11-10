// $Id$
/**
 * @file print-RAW-HLTdecision.C
 * @brief Print HLT decisions per event of RAW data
 * The uses the RawReader interface for data access and can be used
 * with different types of data as long as a RawReader is implemented.
 *
 * <pre>
 * Usage: aliroot -b -q print-RAW-HLTdecision.C'(input, minEvent, maxEvent)'
 * 
 * For help: aliroot -b -q print-RAW-HLTdecision.C
 * </pre>
 *
 * The input file can be a file on Grid like e.g.
 * "alien:///alice/data/2009/LHC09d/000104321/raw/09000104321018.30.root"
 * In that case you need a valid token in order to connect to the Grid.
 * Use 'alien-token-init' from your alien installation.
 *
 * @author Matthias.Richter@ift.uib.no
 * @ingroup alihlt_programs
 */
void print_RAW_HLTdecision(const char* rawFileName,
			   int minEvent=0, int maxEvent=-1)
{
  AliHLTLogging log;
  log.SetGlobalLoggingLevel(0x7c);

  TString strfile=rawFileName;
  if (strfile.Contains("://") && !strfile.Contains("local://")) {
    TGrid::Connect("alien");
  }

  AliRawReader* rawReader=NULL;
  if (rawFileName && rawFileName[0]!=0) rawReader=AliRawReader::Create(rawFileName);
  if (!rawFileName || rawFileName[0]==0 || !rawReader) {
    if (rawFileName && rawFileName[0]!=0)
      cerr << "can not open raw reader " << rawFileName << endl;
    else
      cerr << "print-RAW-HLTdecision.C Print HLT decisions per event of RAW data" << endl;
    cerr << "===============================================================" << endl;
    cerr << "usage: aliroot -b -q -l print-RAW-HLTdecision.C'(rawFileName" << endl;
    cerr << "                                                 minEvent, maxEvent)'" << endl << endl;
    cerr << "  Parameter:" << endl;
    cerr << "       rawFileName      e.g \"raw.root\", \"./\"" << endl;
    cerr << "       minEvent         first event (optional)" << endl;
    cerr << "       maxEvent         last event (optional)" << endl;
    cerr << "===============================================================" << endl;
    return;
  }

  rawReader->RewindEvents();
  int event=0;
  if (!rawReader->NextEvent()) {
    cout << "no events found in " << rawFileName << endl;
    return;
  }

  AliHLTOUT* pHLTOUT=AliHLTOUT::New(rawReader);
  if (!pHLTOUT) {
    cerr << "can not create HLTOUT instance" << endl;
    return;
  }

  do {
    if (minEvent>=0 && event<minEvent) continue;
    cout << "=====================  event " << event << "  ==========================" << endl;
    if (pHLTOUT->Init()<0) {
      cerr << "failed to initialize HLTOUT for event " << event << endl;
      break;
    }
    bool found=false;
    
    // in the original implementation, the GlobalTriggerDecision has been
    // sent with data type kAliHLTDataTypeTObject, thats why we check for
    // both
    if (pHLTOUT->SelectFirstDataBlock(kAliHLTDataTypeGlobalTrigger)>=0 ||
        pHLTOUT->SelectFirstDataBlock(kAliHLTDataTypeTObject)>=0) {
      do {
        TObject* decision = pHLTOUT->GetDataObject();
        if (decision) {
	  if (decision->IsA() == AliHLTGlobalTriggerDecision::Class()) {
	    cout << "HLT Global Trigger: " << decision->GetOption() << "   " << decision->GetTitle() << endl;
	    decision->Print();
	    found=true;
	  }
	  pHLTOUT->ReleaseDataObject(decision);
        }
      } while (!found && pHLTOUT->SelectNextDataBlock());
    }

    if (!found) {
      cout << "   no HLT decision found" << endl;
    }
    pHLTOUT->Reset();
    event++;
  } while (rawReader->NextEvent() && (maxEvent<0 || event<maxEvent));
}

void print_RAW_HLTdecision()
{
  // print usage
  print_RAW_HLTdecision(NULL);
}
