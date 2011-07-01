#ifndef KiddiePoolClasses_h
#define KiddiePoolClasses_h

#include <vector>
#include <deque>
#include <Rtypes.h>
#include <Riostream.h>
#include <TObject.h>
#include <TFile.h>
#include <TMath.h>
#include <TSystem.h>

class MiniTrack;

// Low-memory event mixing classes:
// - KiddiePoolManager: Maintain 2D vector (z,cent) of KiddiePools.
// - KiddiePool: collections of sub-events in z and mult/cent bins.
// - MiniTrack: super-lightweight track "struct" class.
//
// Stores a buffer of tracks that updates continuously. The track type
// contained by the pools can be anything inheriting from
// TObject. Pools are updated based on maintaining a minimum fixed
// number of tracks. Multiplicity/centrality and z-vertex bins must be
// passed in at initialization. For example of implementation, see
// $ALICE_ROOT/PWG4/JetTasks/AliAnalysisTaskPhiCorrelations.cxx
//
// Authors: A. Adare and C. Loizides

typedef std::vector< MiniTrack > MiniEvent;

class MiniTrack
{
 public:
  MiniTrack() : fPt(0), fEta(0), fPhi(0), fSign(0) {};
  MiniTrack(Double_t pt, Double_t eta, Double_t phi, Int_t sign) :
   fPt(pt), fEta(eta), fPhi(phi), fSign(sign>=0?1:0) {};
 
  Float_t Eta()  const { return fEta; }
  Float_t Phi()  const { return fPhi; }
  Float_t Pt()   const { return fPt;  }
  Int_t   Sign() const { if (fSign) return 1; return -1; }

 protected: 
  Float_t fPt;
  Float_t fEta;
  Float_t fPhi;
  Bool_t  fSign;
};

class KiddiePool : public TObject
{
 public:
 KiddiePool(Int_t d) : 
   fEvents(0),
   fNTracksInEvent(0),
   fEventIndex(0),
   fMixDepth(d), 
   fMultMin(-999), 
   fMultMax(+999), 
   fZvtxMin(-999), 
   fZvtxMax(+999), 
   fWasUpdated(0), 
   fMultBinIndex(0), 
    fZvtxBinIndex(0), 
   fDebug(0), 
   fTargetTrackDepth(0),
   fFirstFilled(0),
   fNTimes(0) {;}
  
 KiddiePool(Int_t d, Double_t multMin, Double_t multMax, 
	    Double_t zvtxMin, Double_t zvtxMax) : 
   fEvents(0),
   fNTracksInEvent(0),
   fEventIndex(0),
   fMixDepth(d), 
   fMultMin(multMin), 
   fMultMax(multMax), 
   fZvtxMin(zvtxMin),
   fZvtxMax(zvtxMax),
   fWasUpdated(0),
    fMultBinIndex(0), 
   fZvtxBinIndex(0),
   fDebug(0),
   fTargetTrackDepth(0),
   fFirstFilled(0),
   fNTimes(0) {;}
  ~KiddiePool();
  
  Bool_t      EventMatchesBin(Int_t mult,    Double_t zvtx) const;
  Bool_t      EventMatchesBin(Double_t mult, Double_t zvtx) const;
  Bool_t      IsReady()                    const { return NTracksInPool() >= fTargetTrackDepth; }
  Bool_t      IsFirstReady()               const { return fFirstFilled;   }
  Int_t       GetNTimes()                  const { return fNTimes;        }
  Int_t       GetCurrentNEvents()          const { return fEvents.size(); }
  Int_t       GlobalEventIndex(Int_t j)    const;
  MiniEvent*  GetEvent(Int_t i)            const;
  Int_t       MultBinIndex()               const { return fMultBinIndex; }
  Int_t       NTracksInEvent(Int_t iEvent) const;
  Int_t       NTracksInCurrentEvent()      const { if (fNTracksInEvent.empty()) return 0;  
                                                   return fNTracksInEvent.back(); }
  void        PrintInfo()                  const;
  Int_t       NTracksInPool()              const;
  Bool_t      WasUpdated()                 const { return fWasUpdated; }
  Int_t       ZvtxBinIndex()               const { return fZvtxBinIndex; }
  void        SetDebug(Bool_t b)                 { fDebug = b; }
  void        SetTargetTrackDepth(Int_t d)       { fTargetTrackDepth = d; }
  Int_t       SetEventMultRange(Int_t    multMin, Int_t multMax);
  Int_t       SetEventMultRange(Double_t multMin, Double_t multMax);
  Int_t       SetEventZvtxRange(Double_t zvtxMin, Double_t zvtxMax);
  void        SetMultBinIndex(Int_t iM)          { fMultBinIndex = iM; }
  void        SetZvtxBinIndex(Int_t iZ)          { fZvtxBinIndex = iZ; }
  Int_t       UpdatePool(MiniEvent* miniEvt);
  
protected:
  deque< MiniEvent* >   fEvents;              //Pool of events storing Minitracks
  deque<int>            fNTracksInEvent;      //Tracks in event
  deque<int>            fEventIndex;          //Original event index
  Int_t                 fMixDepth;            //Number of evts. to mix with
  Double_t              fMultMin, fMultMax;   //Track multiplicity bin range
  Double_t              fZvtxMin, fZvtxMax;   //Event z-vertex bin range
  Bool_t                fWasUpdated;          //Evt. succesfully passed selection?
  Int_t                 fMultBinIndex;        //Multiplicity bin
  Int_t                 fZvtxBinIndex;        //Zvertex bin
  Int_t                 fDebug;               //If 1 then debug on
  Int_t                 fTargetTrackDepth;    //Number of tracks, once full
  Bool_t                fFirstFilled;         //Init to false
  Int_t                 fNTimes;              //Number of times init. condition reached

  ClassDef(KiddiePool,1) // Event pool class
};

class KiddiePoolManager : public TObject
{
public:
  KiddiePoolManager() : 
    fDebug(0),
    fNMultBins(0), 
    fNZvtxBins(0),
    fEvPool(0),
    fTargetTrackDepth(0) {;}
  KiddiePoolManager(Int_t maxEvts, Int_t minNTracks,
		    Int_t nMultBins, Double_t *multbins,
		    Int_t nZvtxBins, Double_t *zvtxbins);

  ~KiddiePoolManager();
  
  // First uses bin indices, second uses the variables themselves.
  KiddiePool *GetEventPool(Int_t iMult, Int_t iZvtx) const;
  KiddiePool *GetEventPool(Int_t    centVal, Double_t zvtxVal) const;
  KiddiePool *GetEventPool(Double_t centVal, Double_t zvtxVal) const;

  Int_t       InitEventPools(Int_t depth, 
			     Int_t nmultbins, Double_t *multbins, 
			     Int_t nzvtxbins, Double_t *zvtxbins);
  void        SetTargetTrackDepth(Int_t d) { fTargetTrackDepth = d;}
  Int_t       UpdatePools(MiniEvent*);
  void        SetDebug(Bool_t b) { fDebug = b; }

protected:
  Int_t      fDebug;                         // If 1 then debug on
  Int_t      fNMultBins;                     // number mult bins
  Int_t      fNZvtxBins;                     // number vertex bins
  vector< vector<KiddiePool*> > fEvPool;     // pool in bins of [fMultBin][fZvtxBin]
  Int_t      fTargetTrackDepth;              // Required track size, same for all pools.

  ClassDef(KiddiePoolManager,1)
};
#endif
