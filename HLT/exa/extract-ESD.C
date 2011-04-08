// $Id$
/**
 * @file extract-ESD.C
 * @brief Tool to extract the HLT ESD from the raw data.
 *
 * Usage:
 * <pre>
 *   aliroot -b -q 'extract-ESD.C("raw.root", cdb, nOfEvents)' | tee extract-ESD.log
 *   aliroot -b -q 'extract-ESD.C("alien:///alice/data/2010/.../.root", "raw://", nOfEvents)' | tee extract-ESD.log
 * </pre>
 *
 * @author Kalliopi.Kanaki@ift.uib.no
 * @ingroup alihlt_tutorial
 */
void extract_ESD(const char* input, const char *cdbURI, int nofEvents=-1){

  // setup of the RawReader
  if(!input){
    cerr << "invalid path" << endl;
    cerr << "usage: aliroot -b -q extract-ESD.C'(\"raw.root\", \"cdb\", nOfEvents)'" << endl;
    return;
  }

  TString strfile = input;
  TString struri  = cdbURI;
  
  if(strfile.Contains("://") && !strfile.Contains("local://")) TGrid::Connect("alien");
  
  // Set the CDB storage location
  AliCDBManager *man = AliCDBManager::Instance();
  man->SetDefaultStorage(cdbURI);
  //man->SetRun(146018);

  AliRawReader* pRawReader=AliRawReader::Create(input);
  if (!pRawReader){
    cout << "cannot open RawReader for file " << input << endl;
    return;
  }
  
  if (!pRawReader->NextEvent()){
    cerr << "no events available" << endl;
    return;
  }
  
  if(struri.Contains("local") || strfile.Contains("alien")) man->SetRun(pRawReader->GetRunNumber());

  // setup of the HLT system
  AliHLTSystem *pHLT = AliHLTPluginBase::GetInstance();
  if(!pHLT) cerr << "fatal error: cannot get HLT instance" << endl;

  TString arg;
  arg.Form("-typeid ALIESDV0");
  AliHLTConfiguration publisher("hltout-publisher", "AliHLTOUTPublisher" , NULL, arg.Data());

  arg.Form("-treename HLTesdTree");
  AliHLTConfiguration collector("sink1","EsdCollector","hltout-publisher", arg.Data());

  // the reconstructor setup
  AliHLTReconstructor hltRec;
  hltRec.SetOption("libAliHLTUtil.so loglevel=0x7c chains=sink1 ignore-hltout ignore-ctp");
  if (hltRec.Init()<0) {
    cerr << "initialization of reconstructor failed" << endl;
    return;
  }

  // the reconstruction loop
  Int_t event=0;
  UChar_t* pData=NULL;
  pRawReader->RewindEvents();
  while(pRawReader->NextEvent() && (nofEvents<0 || event<nofEvents)){
    cout << "=======================================================" << endl;
    cout << "event " << event << endl;
    cout << "-------------------------------------------------------" << endl;
    pRawReader->Reset();
    hltRec.Reconstruct(pRawReader, NULL);
    event++;
  }
}

void extract_ESD(const char* input, int nofEvents=-1){
  extract_ESD(input, "raw://", nofEvents);
}

void extract_ESD(){
  cerr << "===============================================================" << endl;
  cerr << "usage: aliroot -b -q extract-ESD.C'(\"raw.root\", \"cdb\", nOfEvents)'"<< endl;
  cerr <<"       aliroot -b -q extract-ESD.C'(\"alien:///alice/data/2010/LHC10../../.root\", \"raw://\", nOfEvents)'"<< endl;
  cerr << "===============================================================" << endl;
}
