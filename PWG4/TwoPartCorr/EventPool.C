// $Id:$

#include "EventPool.h"

bool debug = true;

void EventPool::PrintInfo()
{
  cout << " --- --- --- " << endl;
  cout << Form("%20s: %d", "Pool capacity", fMixDepth) << endl;
  cout << Form("%20s: %d", "Current depth", Depth()) << endl;
  cout << Form("%20s: %d-%d", "Mult. range", fMultMin, fMultMax) << endl;
  cout << Form("%20s: %.1f-%.1f", "Z-vtx range", fZvtxMin, fZvtxMax) << endl;

  return;
}

Bool_t EventPool::IsPoolReady()
{
  return (Int_t)fEvents.size()==fMixDepth;
}

Bool_t EventPool::EventMatchesBin(Int_t mult, Short_t zvtx)
{
  // N.B. Lower bin limit included; upper limit excluded.
  Bool_t multOK = (mult >= fMultMin && mult < fMultMax);
  Bool_t zvtxOK = (zvtx >= fZvtxMin && zvtx < fZvtxMax);
  return (multOK && zvtxOK);
}

Int_t EventPool::Depth()
{
  return fEvents.size();
}
Int_t EventPool::TracksInPool()
{
  int ntrk=0;
  for (int i=0; i<(int)fEvents.size(); i++) {
    ntrk += fNTracksInEvent.at(i);
  }
  return ntrk;
}

Int_t EventPool::SetEventMultRange(int multMin, int multMax)
{
  fMultMin = multMin;
  fMultMax = multMax;
  return 0;
}
Int_t EventPool::SetEventZvtxRange(int zvtxMin, int zvtxMax)
{
  fZvtxMin = zvtxMin;
  fZvtxMax = zvtxMax;
  return 0;
}

Int_t EventPool::UpdatePool(int iEvent, MyHeader* ev, TClonesArray* trk)
{
  // Initialize at any chosen starting event
  if (!fTracks) fTracks = new TClonesArray("MyPart", 1000);

  fMult = trk->GetEntries();
  fZvtx = ev->fVz; // Short_t. TODO: check--should fVc be used instead?

  if (!EventMatchesBin(fMult, fZvtx)) {
    fWasUpdated = false;
    return -1;
  }

  // Should see evsize = trsize (= fMixDepth once full).
  Int_t evsize = fEvents.size();
  Int_t ntsize = fNTracksInEvent.size();

  if (evsize != ntsize) 
    cout << "WARNING:  Event array and track counter array sizes do not match."
	 << " evsize = " << evsize
	 << " ntsize = " << ntsize
	 << endl;

  Bool_t firstReachedCapacity = false;
  if (evsize == fMixDepth - 1) firstReachedCapacity = true;

  // Remove 0th element before appending this event
  if (evsize >= fMixDepth) {
    // TODO: find out - does popping delete the elements from memory?
    fNTracksInEvent.pop_front(); // remove first int
    fEvents.pop_front();        // remove first track array 
    fEventIndex.pop_front();
  }

  // TODO: Clone() is expensive. Maybe a better way available?
  fTracks = (TClonesArray*) trk->Clone();

  fNTracksInEvent.push_back(ev->fNSelTracks);
  fEvents.push_back(fTracks);
  fEventIndex.push_back(iEvent);

  if (firstReachedCapacity) {
    cout << "Pool " << MultBinIndex() << ", " << ZvtxBinIndex() 
	 << " ready at event "<< iEvent;
    PrintInfo();
  }

  fWasUpdated = true;

  bool print = false;
  if (debug && print) {
    cout << " Event " << fEventIndex.back();
    cout << " PoolDepth = " << Depth();
    cout << " NTracks = " << NTracksInCurrentEvent();
    // for (int i=1; i<Depth(); i++) {
    //   cout << " " << NTracksInEvent(iEvent-i);
    // }
    cout << " TracksInPool = " << TracksInPool();
    //    cout << endl;
  }

  return 0;
}

MyPart* EventPool::GetRandomTrack()
{
  MyPart* trk = 0;
  TClonesArray* tca = 0;
  
 // index of random event in the pool
  UInt_t ranEvt = gRandom->Integer(fEvents.size()-1);
  tca = fEvents.at(ranEvt);

 // index of random track in event
  UInt_t ranTrk = gRandom->Integer(tca->GetEntries());
  
  trk = (MyPart*)tca->At(ranTrk);
  
  /*
  if (debug) {
    // loop over the event buffer
    for (int i=0; i<(int)fEvents.size(); i++) {
      TClonesArray* tca = fEvents.at(i);
      cout << tca->GetEntries() << " ";
    }
    cout << endl;
    // Get first track inside each
    for (int i=0; i<(int)fEvents.size(); i++) {
      TClonesArray* tca = fEvents.at(i);
      trk = (MyPart*)tca->At(0);
      cout << trk->Pt() << " ";
    }
    cout << endl;
  }
  */

  return trk;
}

Int_t EventPool::NTracksInCurrentEvent()
{
  return fNTracksInEvent.back();
}

// Important!! This fn. will break if selective filling is implemented
// (i.e. pool is not updated for every single event)
// TODO: Implement internal counters to fix this.
Int_t EventPool::NTracksInEvent(int iEvent)
{
  Int_t n = -1;
  Int_t curEvent = fEventIndex.back();
  Int_t offset = curEvent - iEvent;
  // pos = position of iEvent in rolling buffer
  int pos = fEventIndex.size() - offset - 1;

  if (offset==0) // (iEvent == curEvent)
    n = fNTracksInEvent.back();
  else if (offset < 0 || iEvent < 0) {
    n = 0;//    cout << " No info for event " << iEvent << " ";
  }
  else if (offset > 0 && offset <= (int)fEventIndex.size()) {
    n = fNTracksInEvent.at(pos);
  }
  else
    cout << "Event info no longer in memory" << endl;
  return n;
}


