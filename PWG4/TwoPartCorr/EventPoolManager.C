// $Id$

#include "EventPoolManager.h"

ClassImp(GenericEventPool)

void GenericEventPool::PrintInfo() const
{
  cout << " --- --- --- " << endl;
  cout << Form("%20s: %d", "Pool capacity", fMixDepth) << endl;
  cout << Form("%20s: %d", "Current depth", Depth()) << endl;
  cout << Form("%20s: %d to %d", "Mult. range", fMultMin, fMultMax) << endl;
  cout << Form("%20s: %.1f to %.1f", "Z-vtx range", fZvtxMin, fZvtxMax) << endl;

  return;
}

Bool_t GenericEventPool::EventMatchesBin(Int_t mult, Short_t zvtx) const
{
  // N.B. Lower bin limit included; upper limit excluded.

  Bool_t multOK = (mult >= fMultMin && mult < fMultMax);
  Bool_t zvtxOK = (zvtx >= fZvtxMin && zvtx < fZvtxMax);
  return (multOK && zvtxOK);
}

Int_t GenericEventPool::TracksInPool() const
{
  Int_t ntrk=0;
  for (Int_t i=0; i<(Int_t)fEvents.size(); ++i) {
    ntrk += fNTracksInEvent.at(i);
  }
  return ntrk;
}

Int_t GenericEventPool::SetEventMultRange(Int_t multMin, Int_t multMax)
{
  fMultMin = multMin;
  fMultMax = multMax;
  return 0;
}

Int_t GenericEventPool::SetEventZvtxRange(Int_t zvtxMin, Int_t zvtxMax)
{
  fZvtxMin = zvtxMin;
  fZvtxMax = zvtxMax;
  return 0;
}

Int_t GenericEventPool::UpdatePool(Int_t iEvent, const MyHeader *ev, TObjArray *trk)
{
  // Initialize at any chosen starting event

  Int_t mult = trk->GetEntries();
  Double_t zvtx = ev->fVz;

  if (!EventMatchesBin(mult, zvtx)) {
    fWasUpdated = false;
    return -1;
  }

  fMult = mult;
  fZvtx = zvtx;

  // Should see evsize = trsize (= fMixDepth once full).
  Int_t evsize = fEvents.size();
  Int_t ntsize = fNTracksInEvent.size();

  if (evsize != ntsize) 
    cout << "WARNING:  Event array and track counter array sizes do not match."
	 << " evsize = " << evsize
	 << " ntsize = " << ntsize
	 << endl;

  Bool_t firstReachedCapacity = false;
  if (evsize == fMixDepth - 1) 
    firstReachedCapacity = true;

  // Remove 0th element before appending this event
  if (evsize >= fMixDepth) {
    TObjArray *fa = fEvents.front();
    delete fa;
    fEvents.pop_front();         // remove first track array 
    fNTracksInEvent.pop_front(); // remove first int
    fEventIndex.pop_front();
  }

  fNTracksInEvent.push_back(mult);
  fEvents.push_back(trk);
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

  return fEvents.size();
}

TObject* GenericEventPool::GetRandomTrack() const
{
  UInt_t ranEvt = gRandom->Integer(fEvents.size());
  TObjArray *tca = fEvents.at(ranEvt);
  UInt_t ranTrk = gRandom->Integer(tca->GetEntries());
  TObject *trk = (TObject*)tca->At(ranTrk);
  return trk;
}

Int_t GenericEventPool::NTracksInEvent(Int_t iEvent) const
{
  Int_t n = -1;
  Int_t curEvent = fEventIndex.back();
  Int_t offset = curEvent - iEvent;
  Int_t pos = fEventIndex.size() - offset - 1;

  if (offset==0)
    n = fNTracksInEvent.back();
  else if (offset < 0 || iEvent < 0) {
    n = 0;
  }
  else if (offset > 0 && offset <= (int)fEventIndex.size()) {
    n = fNTracksInEvent.at(pos);
  }
  else
    cout << "Event info no longer in memory" << endl;
  return n;
}

ClassImp(EventPoolManager)

Int_t EventPoolManager::InitEventPools(Int_t depth, 
                                       Int_t nMultBins, Double_t *multbin, 
                                       Int_t nZvtxBins, Double_t *zvtxbin)
{
  // First assign EventPoolManager members
  fNMultBins = nMultBins;
  fNZvtxBins = nZvtxBins;

  for (Int_t iM=0; iM<nMultBins; iM++) {
    std::vector<GenericEventPool*> evp;
    for (Int_t iZ=0; iZ<nZvtxBins; iZ++) {
      evp.push_back(new GenericEventPool(depth, 
                                         multbin[iM], multbin[iM+1], 
                                         zvtxbin[iZ], zvtxbin[iZ+1] ));
    }
    fEvPool.push_back(evp);
  }
  
  for (Int_t iM=0; iM<nMultBins; iM++) {
    for (Int_t iZ=0; iZ<nZvtxBins; iZ++) {
      fEvPool.at(iM).at(iZ)->SetMultBinIndex(iM);
      fEvPool.at(iM).at(iZ)->SetZvtxBinIndex(iZ);
    }
  }
  
  if (0) {
    cout << "fEvPool outer size: " << fEvPool.size() << endl;
    for (Int_t iM=0; iM<nMultBins; iM++) {
      for (Int_t iZ=0; iZ<nZvtxBins; iZ++) {
	if(fEvPool.at(iM).at(iZ)) {
	  cout << "multiplicity bin: " << iM;
	  cout << ", z-vertex bin: " << iZ;
	  fEvPool.at(iM).at(iZ)->PrintInfo();
	}
      }
    }
  }
  
  return fEvPool.size();
}

GenericEventPool *EventPoolManager::GetEventPool(Int_t iMult, Int_t iZvtx) const
{
  if (iMult < 0 || iMult >= fNMultBins) 
    return 0x0;
  if (iZvtx < 0 || iZvtx >= fNZvtxBins) 
    return 0x0;
  return fEvPool.at(iMult).at(iZvtx);
}

Int_t EventPoolManager::UpdatePools(Int_t iEvent, const MyHeader* ev, TObjArray *trk)
{
  for (Int_t iM=0; iM<fNMultBins; iM++) {
    for (Int_t iZ=0; iZ<fNZvtxBins; iZ++) {
      if (fEvPool.at(iM).at(iZ)->UpdatePool(iEvent, ev, trk)>-1)
        break;
    }
  }  
  return 0;
}
