#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TFile.h"

#include "AliRawReaderRoot.h"
#include "AliTRDrawStream.h"
#include "AliTRDdigitsManager.h"
#include "AliESDTrdTrack.h"
#endif

void AliTRDreadRaw(TString filename, Int_t firstEvent = 0, Int_t nEvents = 100,
	      Bool_t readDigits = kTRUE, Bool_t readAtOnce = kTRUE, Int_t debugLevel = -1)
{
  // this macro demonstrates the usage of the TRD rawstream
  // including access to error codes, event statistics, tracklets and tracks

  Bool_t dump = kFALSE; // dump data from specified MCMs
  Bool_t oldReadoutOrder = kTRUE; // set correct readout order

  // create raw reader to retrieve data
  AliRawReader *reader = new AliRawReaderRoot(filename.Data());
  reader->Select("TRD");
  reader->GotoEvent(firstEvent);

  // create raw stream
  AliTRDrawStream *rawStream = new AliTRDrawStream(reader);

  // enable debug output for the rawStream
  if (debugLevel > 0) {
    AliLog::SetPrintLocation(kFALSE);
    AliLog::SetClassDebugLevel("AliTRDrawStream", debugLevel);
  }

  // create the digits manager if wanted
  AliTRDdigitsManager *digMgr = 0x0;
  if (readDigits) {
     digMgr = new AliTRDdigitsManager();
     digMgr->CreateArrays();
     rawStream->SetDigitsManager(digMgr);
  }

  // setup MCM readout order and enable error message in case of wrong order
  if (!oldReadoutOrder) {
    for (Int_t iMcm = 0; iMcm < 16; iMcm++)
      AliTRDrawStream::SetMCMReadoutPos(iMcm, iMcm);
  }
  AliTRDrawStream::SetErrorDebugLevel(AliTRDrawStream::kPosUnexp, 0);

  // if you want you can dump data from individual MCMs
  // add further MCMs as needed
  if (dump) {
    rawStream->SetDumpMCM(330, 0, 10);
  }

  Int_t iEvent = firstEvent;

  // create the array to hold the tracklets
  TClonesArray *trklArray = new TClonesArray("AliTRDtrackletWord", 500);

  // create the array to hold the tracks
  TClonesArray *trkArray = new TClonesArray("AliESDTrdTrack", 50);

  // open a file to hold output data (for testing only)
  TFile *f = TFile::Open("raw-out.root", "RECREATE");
  TTree *trackingTree = new TTree("trackingTree", "tree with tracklets and tracks");
  trackingTree->Branch("event", &iEvent);
  trackingTree->Branch("trkl", &trklArray);
  trackingTree->Branch("trkl", &trkArray);

  TTree *eventStats = new TTree("stats", "event statistics");
  AliTRDrawStream::AliTRDrawStats *stats = rawStream->GetStats();
  eventStats->Branch("event", &iEvent);
  eventStats->Branch("stats", &stats);

  // loop over events
  while (reader->NextEvent()) {
    iEvent++;

    // clear from previous event
    if (digMgr) {
      for (Int_t iDet = 0; iDet < 540; iDet++)
	digMgr->ClearArrays(iDet);
    }
    stats->ClearStats();

    if (readAtOnce) {
      rawStream->ReadEvent();
    }
    else {
      Int_t det;

      while (rawStream->NextDDL()) {
	while ((det = rawStream->NextChamber(digMgr)) > -1)
	;
      }
    }

    // store the output
    trackingTree->Fill();
    eventStats->Fill();

    if (iEvent >= (firstEvent + nEvents))
      break;
  }

  // retrieve error tree and write output to file
  TTree *t = rawStream->GetErrorTree();
  f->WriteTObject(t);
  f->WriteTObject(trackingTree);
  eventStats->Write();

  f->Close();
}
