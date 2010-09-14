// $Id:$

#ifndef EventPool_h
#define EventPool_h

#include <vector>
#include <deque>
#include "Rtypes.h"
#include <Riostream.h>
#include <TClonesArray.h>
#include "TFile.h"
#include "TMath.h"
#include "TRandom.h"
#include "TSystem.h"
#include "TreeClasses.h"

class EventPool : public TObject
{
public:
  EventPool(Int_t d) : fMixDepth(d), fMultMin(-999), fMultMax(-999), 
    fZvtxMin(-999), fZvtxMax(-999) {;}
  EventPool(Int_t d, Int_t multMin, Int_t multMax, 
	    Double_t zvtxMin, Double_t zvtxMax) : 
    fMixDepth(d), fMultMin(multMin), fMultMax(multMax), 
    fZvtxMin(zvtxMin), fZvtxMax(zvtxMax) {;}
  ~EventPool(){;}

  // TODO: implement cleanup
  // Int_t Drain();
  Int_t SetEventMultRange(int multMin, int multMax);
  Int_t SetEventZvtxRange(int zvtxMin, int zvtxMax);
  
  Int_t UpdatePool(int iEvent, MyHeader* ev, TClonesArray* trk);
  Bool_t WasUpdated() { return fWasUpdated; }
  Bool_t EventMatchesBin(Int_t mult, Short_t zvtx);
  Bool_t IsPoolReady();              // Contains fMixDepth events? 
  Int_t Depth();                     // Number of events in pool
  Int_t TracksInPool();              // Total tracks in pool
  Int_t NTracksInEvent(int iEvent);
  Int_t NTracksInCurrentEvent();
  MyPart* GetRandomTrack();
  void PrintInfo();
  Int_t MultBinIndex() { return fMultBinIndex; }
  Int_t ZvtxBinIndex() { return fZvtxBinIndex; }
  void SetMultBinIndex(Int_t iM) { fMultBinIndex = iM; }
  void SetZvtxBinIndex(Int_t iZ) { fZvtxBinIndex = iZ; }

protected:
  TClonesArray* fTracks; // Copy of trk array. Refreshes each event.

  //  Use as ring buffers
  deque<TClonesArray*> fEvents; // The guts of the class. Holds TObjArrays of MyParts.
  deque<int> fNTracksInEvent;
  deque<int> fEventIndex;

  Int_t fMixDepth;             // Number of evts. to mix with
  Int_t fMultMin, fMultMax;    // Track multiplicity bin range
  Double_t fZvtxMin, fZvtxMax; // Event z-vertex bin range
  Int_t fMult;                 // Tracks in current event
  Short_t fZvtx;               // Current z-vertex
  Bool_t fWasUpdated;          // Evt. succesfully passed selection?

  Int_t fMultBinIndex;
  Int_t fZvtxBinIndex;

  ClassDef(EventPool,1)
};
#endif



