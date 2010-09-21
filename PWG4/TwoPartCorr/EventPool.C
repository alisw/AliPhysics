// $Id$

#include "EventPool.h"

ClassImp(EventPool)

void EventPool::PrintInfo() const
{
  cout << " --- --- --- " << endl;
  cout << Form("%20s: %d", "Pool capacity", fMixDepth) << endl;
  cout << Form("%20s: %d", "Current depth", Depth()) << endl;
  cout << Form("%20s: %d to %d", "Mult. range", fMultMin, fMultMax) << endl;
  cout << Form("%20s: %.1f to %.1f", "Z-vtx range", fZvtxMin, fZvtxMax) << endl;

  return;
}

Bool_t EventPool::EventMatchesBin(Int_t mult, Short_t zvtx) const
{
  // N.B. Lower bin limit included; upper limit excluded.

  Bool_t multOK = (mult >= fMultMin && mult < fMultMax);
  Bool_t zvtxOK = (zvtx >= fZvtxMin && zvtx < fZvtxMax);
  return (multOK && zvtxOK);
}

Int_t EventPool::TracksInPool() const
{
  Int_t ntrk=0;
  for (Int_t i=0; i<(Int_t)fEvents.size(); ++i) {
    ntrk += fNTracksInEvent.at(i);
  }
  return ntrk;
}

Int_t EventPool::SetEventMultRange(Int_t multMin, Int_t multMax)
{
  fMultMin = multMin;
  fMultMax = multMax;
  return 0;
}

Int_t EventPool::SetEventZvtxRange(Int_t zvtxMin, Int_t zvtxMax)
{
  fZvtxMin = zvtxMin;
  fZvtxMax = zvtxMax;
  return 0;
}

Int_t EventPool::UpdatePool(Int_t iEvent, const MyHeader *ev, TClonesArray *trk)
{
  // Initialize at any chosen starting event
  if (!fTracks) 
    fTracks = new TClonesArray("MyPart", 1000);

  fMult = trk->GetEntries();
  fZvtx = ev->fVz;

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

  bool print = 1;
  if (fDebug && print) {
    cout << " Event " << fEventIndex.back();
    cout << " PoolDepth = " << Depth();
    cout << " NTracks = " << NTracksInCurrentEvent();
    cout << " TracksInPool = " << TracksInPool();
  }

  return 0;
}

MyPart* EventPool::GetRandomTrack() const
{
  MyPart* trk = 0;
  TClonesArray* tca = 0;
  UInt_t ranEvt = gRandom->Integer(fEvents.size());
  tca = fEvents.at(ranEvt);
  UInt_t ranTrk = gRandom->Integer(tca->GetEntries());
  trk = (MyPart*)tca->At(ranTrk);
  return trk;
}

// Important!! This fn. will break if selective filling is implemented
// (i.e. pool is not updated for every single event)
// TODO: Implement internal counters to fix this.
Int_t EventPool::NTracksInEvent(Int_t iEvent) const
{
  Int_t n = -1;
  Int_t curEvent = fEventIndex.back();
  Int_t offset = curEvent - iEvent;
  // pos = position of iEvent in rolling buffer
  Int_t pos = fEventIndex.size() - offset - 1;

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
