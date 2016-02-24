#ifndef AliEventPoolManager_h
#define AliEventPoolManager_h

#include <vector>
#include <deque>
#include <Rtypes.h>
#include <TObjArray.h>
#include "AliLog.h"
#include <AliVParticle.h>

// Generic event mixing classes
//
// Stores a buffer of tracks that updates continuously. The track type
// contained by the pools can be anything inheriting from
// TObject. Pools are updated based on maintaining a minimum fixed
// number of tracks. Multiplicity/centrality and z-vertex bins must be
// passed in at initialization. For example of implementation, see
// $ALICE_ROOT/PWGCF/Correlations/DPhi/AliAnalysisTaskPhiCorrelations.cxx
//
// Authors: A. Adare and C. Loizides

using std::deque;

class AliEventPool : public TObject
{
 public:
 AliEventPool() 
   : fEvents(0),
    fNTracksInEvent(0),
    fEventIndex(0),
    fMixDepth(1000), 
    fMultMin(-999), 
    fMultMax(+999), 
    fZvtxMin(-999), 
    fZvtxMax(+999), 
    fPsiMin(-999), 
    fPsiMax(+999), 
    fPtMin(-9999), 
    fPtMax(+9999), 
    fWasUpdated(0), 
    fMultBinIndex(0), 
    fZvtxBinIndex(0), 
    fPsiBinIndex(0), 
    fPtBinIndex(0),
    fDebug(0), 
    fTargetTrackDepth(0),
    fFirstFilled(0),
    fLockFlag(0),
    fSaveFlag(0),
    fNTimes(0),
    fTargetFraction(1),
    fTargetEvents(0)  {;} // default constructor needed for correct saving

 AliEventPool(Int_t d) 
   : fEvents(0),
    fNTracksInEvent(0),
    fEventIndex(0),
    fMixDepth(d), 
    fMultMin(-999), 
    fMultMax(+999), 
    fZvtxMin(-999), 
    fZvtxMax(+999), 
    fPsiMin(-999), 
    fPsiMax(+999), 
    fPtMin(-9999), 
    fPtMax(+9999), 
    fWasUpdated(0), 
    fMultBinIndex(0), 
    fZvtxBinIndex(0), 
    fPsiBinIndex(0), 
    fPtBinIndex(0),
    fDebug(0), 
    fTargetTrackDepth(0),
    fFirstFilled(0),
    fLockFlag(0),
    fSaveFlag(0),
    fNTimes(0),
    fTargetFraction(1),
    fTargetEvents(0)  {;}
  

 AliEventPool(Int_t d, Double_t multMin, Double_t multMax, 
        Double_t zvtxMin, Double_t zvtxMax,
        Double_t psiMin=-999., Double_t psiMax=999.,
        Double_t ptMin=-9999., Double_t ptMax=9999.) 
   : fEvents(0),
    fNTracksInEvent(0),
    fEventIndex(0),
    fMixDepth(d), 
    fMultMin(multMin), 
    fMultMax(multMax), 
    fZvtxMin(zvtxMin),
    fZvtxMax(zvtxMax),
    fPsiMin(psiMin),
    fPsiMax(psiMax),
    fPtMin(ptMin), 
    fPtMax(ptMax), 
    fWasUpdated(0),
    fMultBinIndex(0), 
    fZvtxBinIndex(0),
    fPsiBinIndex(0),
    fPtBinIndex(0),
    fDebug(0),
    fTargetTrackDepth(0),
    fFirstFilled(0),
    fLockFlag(0),
    fSaveFlag(0),
    fNTimes(0),
    fTargetFraction(1),
    fTargetEvents(0) {;}
  
  ~AliEventPool() {;}
  
  Bool_t      EventMatchesBin(Int_t mult,    Double_t zvtx, Double_t psi=0., Double_t pt=0.) const;
  Bool_t      EventMatchesBin(Double_t mult, Double_t zvtx, Double_t psi=0., Double_t pt=0.) const;
  Bool_t      IsReady()                    const { return IsReady(NTracksInPool(), GetCurrentNEvents()); }
  Bool_t      IsFirstReady()               const { return fFirstFilled;   }
  Int_t       GetNTimes()                  const { return fNTimes;        }
  Int_t       GetCurrentNEvents()          const { return fEvents.size(); }
  Int_t       GlobalEventIndex(Int_t j)    const;
  TObject    *GetRandomTrack()             const;
  TObjArray  *GetRandomEvent()             const;
  TObjArray  *GetEvent(Int_t i)            const;
  Int_t       MultBinIndex()               const { return fMultBinIndex; }
  Int_t       NTracksInEvent(Int_t iEvent) const;
  Int_t       NTracksInCurrentEvent()      const { return fNTracksInEvent.back(); }
  void        PrintInfo()                  const;
  Int_t       PsiBinIndex()                const { return fPsiBinIndex; }
  Int_t       PtBinIndex()                 const { return fPtBinIndex; }
  Int_t       NTracksInPool()              const;
  Bool_t      WasUpdated()                 const { return fWasUpdated; }
  Int_t       ZvtxBinIndex()               const { return fZvtxBinIndex; }
  void        SetDebug(Bool_t b)                 { fDebug = b; }
  void        SetTargetTrackDepth(Int_t d, Float_t fraction = 1.0) { fTargetTrackDepth = d; fTargetFraction = fraction; }
  void        SetTargetEvents(Int_t ev)          { fTargetEvents = ev; }
  Int_t       SetEventMultRange(Int_t    multMin, Int_t multMax);
  Int_t       SetEventMultRange(Double_t multMin, Double_t multMax);
  Int_t       SetEventZvtxRange(Double_t zvtxMin, Double_t zvtxMax);
  Int_t       SetEventPsiRange(Double_t psiMin, Double_t psiMax);
  Int_t       SetEventPtRange(Double_t ptMin, Double_t ptMax);
  void        SetMultBinIndex(Int_t iM) { fMultBinIndex = iM; }
  void        SetZvtxBinIndex(Int_t iZ) { fZvtxBinIndex = iZ; }
  void        SetPsiBinIndex(Int_t iP) { fPsiBinIndex  = iP; }
  void        SetPtBinIndex(Int_t iPt) { fPtBinIndex = iPt; }
  void        SetLockFlag(Bool_t val) { fLockFlag = val; }
  Bool_t      GetLockFlag() { return fLockFlag; }
  void        SetSaveFlag(Bool_t val) { fSaveFlag = val; }
  Bool_t      GetSaveFlag() { return fSaveFlag; }
  Double_t    GetPtMin() { return fPtMin; }
  Double_t    GetPtMax() { return fPtMax; }
  Double_t    GetPsiMin() { return fPsiMin; }
  Double_t    GetPsiMax() { return fPsiMax; }
  Double_t    GetMultMin() { return fMultMin; }
  Double_t    GetMultMax() { return fMultMax; }
  Double_t    GetZvtxMin() { return fZvtxMin; }
  Double_t    GetZvtxMax() { return fZvtxMax; }

  Int_t       UpdatePool(TObjArray *trk);
  Long64_t    Merge(TCollection* hlist);
//  deque<TObjArray*> GetEvents() { return fEvents; }
  void        Clear();

protected:
  Bool_t      IsReady(Int_t tracks, Int_t events) const { return (tracks >= fTargetFraction * fTargetTrackDepth) || ((fTargetEvents > 0) && (events >= fTargetEvents)); }
  
  deque<TObjArray*>     fEvents;              //Holds TObjArrays of MyTracklets
  deque<int>            fNTracksInEvent;      //Tracks in event
  deque<int>            fEventIndex;          //Original event index
  Int_t                 fMixDepth;            //Number of evts. to mix with
  Double_t              fMultMin, fMultMax;   //Track multiplicity bin range
  Double_t              fZvtxMin, fZvtxMax;   //Event z-vertex bin range
  Double_t              fPsiMin, fPsiMax;     //Event plane angle (Psi) bin range
  Double_t              fPtMin, fPtMax;       //Particle pt bin range
  Bool_t                fWasUpdated;          //Evt. succesfully passed selection?
  Int_t                 fMultBinIndex;        //Multiplicity bin
  Int_t                 fZvtxBinIndex;        //Zvertex bin
  Int_t                 fPsiBinIndex;         //Event plane angle (Psi) bin
  Int_t                 fPtBinIndex;          //Particle pt bin
  Int_t                 fDebug;               //If 1 then debug on
  Int_t                 fTargetTrackDepth;    //Number of tracks, once full
  Bool_t                fFirstFilled;         //Init to false
  Bool_t                fLockFlag;            //if locked, no update is allowed. Useful for external pools
  Bool_t                fSaveFlag;            //flag whether to save the pool to the output file or not
  Int_t                 fNTimes;              //Number of times init. condition reached
  Float_t               fTargetFraction;      //fraction of fTargetTrackDepth at which pool is ready (default: 1.0)
  Int_t                 fTargetEvents;        //if non-zero: number of filled events after which pool is ready regardless of fTargetTrackDepth (default: 0)

  ClassDef(AliEventPool,4) // Event pool class
};

class AliEventPoolManager : public TObject
{
public:
  AliEventPoolManager() 
    : fDebug(0),
    fNMultBins(0), 
    fNZvtxBins(0),
    fNPsiBins(0),
    fNPtBins(0),
    fMultBins(),
    fZvtxBins(),
    fPsiBins(),
    fPtBins(),
    fEvPool(0),
    fTargetTrackDepth(0) {}
  AliEventPoolManager(Int_t maxEvts, Int_t minNTracks,
          Int_t nMultBins, Double_t *multbins,
          Int_t nZvtxBins, Double_t *zvtxbins);
  
  AliEventPoolManager(Int_t maxEvts, Int_t minNTracks,
          Int_t nMultBins, Double_t *multbins,
          Int_t nZvtxBins, Double_t *zvtxbins,
          Int_t nPsiBins, Double_t *psibins);

  AliEventPoolManager(Int_t maxEvts, Int_t minNTracks,
          Int_t nMultBins, Double_t *multbins,
          Int_t nZvtxBins, Double_t *zvtxbins,
          Int_t nPsiBins, Double_t *psibins,
          Int_t nPtBins, Double_t *ptbins);

  AliEventPoolManager(Int_t maxEvts, Int_t minNTracks, const char* binning);


  ~AliEventPoolManager() {;}
  Long64_t    Merge(TCollection* hlist);

  // First uses bin indices, second uses the variables themselves.
  AliEventPool *GetEventPool(Int_t iMult, Int_t iZvtx, Int_t iPsi=0, Int_t iPt=0) const;
  AliEventPool *GetEventPool(Int_t centVal, Double_t zvtxVal, Double_t psiVal=0., Int_t iPt=0) const;
  AliEventPool *GetEventPool(Double_t centVal, Double_t zvtxVal, Double_t psiVal=0., Int_t iPt=0) const;

  Int_t       InitEventPools(Int_t depth, 
                Int_t nMultBins, Double_t *multbin, 
                Int_t nZvtxBins, Double_t *zvtxbin, 
                Int_t nPsiBins, Double_t *psibin,
                Int_t nPtBins, Double_t *ptbin);
  
  void        SetTargetTrackDepth(Int_t d) { fTargetTrackDepth = d;} // Same as for G.E.P. class
  Int_t       UpdatePools(TObjArray *trk);
  void        SetDebug(Bool_t b) { fDebug = b; }
  void        SetTargetValues(Int_t trackDepth, Float_t fraction, Int_t events);
  Int_t       GetNumberOfAllBins() {return fNPtBins*fNMultBins*fNZvtxBins*fNPsiBins;}
  Int_t       GetNumberOfPtBins() {return fNPtBins;}
  Int_t       GetNumberOfMultBins() {return fNMultBins;}
  Int_t       GetNumberOfZVtxBins() {return fNZvtxBins;}
  Int_t       GetNumberOfPsiBins() {return fNPsiBins;}

  void        Validate();
  void        ClearPools();
  void        ClearPools(Double_t minCent, Double_t maxCent,  Double_t minZvtx, Double_t maxZvtx, Double_t minPsi, Double_t maxPsi, Double_t minPt, Double_t maxPt);
  void        SetSaveFlag(Double_t minCent, Double_t maxCent,  Double_t minZvtx, Double_t maxZvtx, Double_t minPsi, Double_t maxPsi, Double_t minPt, Double_t maxPt);

 protected:
  Int_t      fDebug;                                    // If 1 then debug on
  Int_t      fNMultBins;                                // number mult bins
  Int_t      fNZvtxBins;                                // number vertex bins
  Int_t      fNPsiBins;                                 // number Event plane angle (Psi) bins
  Int_t      fNPtBins;                                  // number pt bins

  std::vector<Double_t> fMultBins;                      // mult bins
  std::vector<Double_t> fZvtxBins;                      // vertex bins
  std::vector<Double_t> fPsiBins;                       // Event plane angle (Psi) bins
  std::vector<Double_t> fPtBins;                        // pt bins

  std::vector<AliEventPool*> fEvPool;                   // pool in bins of [fNMultBin][fNZvtxBin][fNPsiBin][fNPtBins]
  Int_t      fTargetTrackDepth;                         // Required track size, same for all pools.

  Int_t       GetBinIndex(Int_t iMult, Int_t iZvtx, Int_t iPsi, Int_t iPt) const {return fNZvtxBins*fNPsiBins*fNPtBins*iMult + fNPsiBins*fNPtBins*iZvtx + fNPtBins*iPsi + iPt;}
  Double_t*   GetBinning(const char* configuration, const char* tag, Int_t& nBins) const;

  ClassDef(AliEventPoolManager,3)
};


#endif
