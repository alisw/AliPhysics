#ifndef AliToyMCReconstruction_H
#define AliToyMCReconstruction_H

#include <TString.h>
#include <TObject.h>
#include <TClonesArray.h>
#include <TExMap.h>

class TTree;

class TTreeSRedirector;
class AliExternalTrackParam;
class AliTPCtrackerSector;
class AliToyMCEvent;
class AliTPCCorrection;
class AliTPCseed;
class AliTPCtrackerRow;
class AliToyMCTrack;
class AliTPCclusterMI;
class AliRieman;
class AliTrackPoint;
class AliTPCParam;

class AliToyMCReconstruction : public TObject {
public:
  AliToyMCReconstruction();
  virtual ~AliToyMCReconstruction();
  
  enum ECorrType {
    kNoCorrection = 0,
    kTPCCenter,
    kAverageEta,
    kIdeal,
    kPreliminaryEta // NOT TO USE (needs fixing!!! Not yet in full code!!!)
  };

  enum EDet {
    kITS=0,
    kTPC,
    kTRD
  };
  
  void RunReco(const char* file, Int_t nmaxEv=-1);
  void RunRecoAllClusters(const char* file, Int_t nmaxEv=-1);
  void RunRecoAllClustersStandardTracking(const char* file, Int_t nmaxEv=-1);

  void RunFullTracking(const char* file, Int_t nmaxEv=-1);
  
  // reconstruction settings
  void      SetRecoSettings(Bool_t idealTracking, Int_t clusterType, ECorrType correctionType, Int_t seedingRow=130, Int_t seedingDist=10)
                           { fIdealTracking=idealTracking; fClusterType=clusterType; fSeedingRow=seedingRow, fSeedingDist=seedingDist, fCorrectionType=correctionType; }
  
  void      SetClusterType(Int_t type)  { fClusterType = type;    }
  Int_t     GetClusterType()  const     { return fClusterType;    }
  
  void      SetSeedingRow(Int_t row)    { fSeedingRow = row;      }
  Int_t     GetSeedingRow()  const      { return fSeedingRow;     }
  
  void      SetSeedingDist(Int_t dist)  { fSeedingDist = dist;    }
  Int_t     GetSeedingDist()  const     { return fSeedingDist;    }
  
  void      SetCorrType(ECorrType type) { fCorrectionType = type; }
  ECorrType GetCorrectionType() const   { return fCorrectionType; }

  void   SetUseMaterialBudget(Bool_t mat) { fUseMaterial = mat;   }
  Bool_t GetUseMaterialBudget() const   { return fUseMaterial;    }

  void   SetIdealTracking(Bool_t tr)    { fIdealTracking = tr;    }
  Bool_t GetIdealTracking()  const      { return fIdealTracking;  }
  
  void   SetTree(TTree *tree) { fTree=tree; }
  TTree* GetTree() const { return fTree; }

  AliExternalTrackParam* GetSeedFromTrack(const AliToyMCTrack * const tr, Bool_t forceSeed=kFALSE);
  AliExternalTrackParam* GetSeedFromTrackIdeal(const AliToyMCTrack * const tr, EDet det );
  AliExternalTrackParam* GetFittedTrackFromSeed(const AliToyMCTrack *tr, const AliExternalTrackParam *seed);
  AliExternalTrackParam* GetFittedTrackFromSeedAllClusters(const AliToyMCTrack *tr, const AliExternalTrackParam *seed, Int_t &nClus);
  AliExternalTrackParam* GetTrackRefit(const AliToyMCTrack * const tr, EDet det);

  AliToyMCTrack *ConvertTPCSeedToToyMCTrack(const AliTPCseed &seed);
  AliExternalTrackParam* GetRefittedTrack(const AliTPCseed &seed);
  
  AliTPCclusterMI* FindMiddleCluster(const AliTPCclusterMI *clInner, const AliTPCclusterMI *clOuter,
                                     Double_t roady, Double_t roadz,
                                     AliRieman &seedFit);

  void AddMiddleClusters(AliTPCseed *seed,
                         const AliTPCclusterMI *clInner, const AliTPCclusterMI *clOuter,
                         Double_t roady, Double_t roadz,
                         Int_t &nTotalClusters, AliRieman &seedFit);
  Int_t MakeSeeds2(TObjArray * arr, Int_t sec, Int_t i1, Int_t i2);
  void MakeSeeds(TObjArray * arr, Int_t sec, Int_t i1, Int_t i2);

  void SetRieman(const AliTPCseed &seed, AliRieman &rieman);
  void CopyRieman(const AliRieman &from, AliRieman &to);
  
  AliExternalTrackParam* ClusterToTrackAssociation(const AliTPCseed *seed, Int_t trackID, Int_t &nClus);
  void ClusterToTrackAssociation(AliTPCseed &seed);
  void AssociateClusters(AliTPCseed &seed, Int_t firstRow, Int_t lastRow, Bool_t direction);
  
  void InitSpaceCharge();

  void SetLongT0seed(Bool_t l) { fLongT0seed=l; }
  Bool_t GetLongT0seed() const { return fLongT0seed; }
  
  static TTree* ConnectTrees(const char* files);
  
  Double_t GetVDrift() const;
  Double_t GetZLength(Int_t roc) const;

  void InitStreamer(TString addName, Int_t level=1);

  void ConnectInputFile(const char* file, Int_t nmaxEv=-1);
  void Cleanup();

  void DumpTracksToTree(const char* file);
  
// private:
  AliToyMCReconstruction(const AliToyMCReconstruction &rec);
  AliToyMCReconstruction& operator= (AliToyMCReconstruction& rec);

  void SetTrackPointFromCluster(const AliTPCclusterMI *cl, AliTrackPoint &p);
  void ClusterToSpacePoint(const AliTPCclusterMI *cl, Float_t xyz[3]);

  Int_t LoadInnerSectors();
  Int_t LoadOuterSectors();
  
  Int_t GetSector(AliExternalTrackParam *track);
  void FillSectorStructure(Int_t maxev);
  void FillSectorStructureAC();

  void SetupTrackMaps();

  void CookLabel(AliTPCseed *seed, Double_t fraction, Int_t info[5]=0);

  void DumpSeedInfo(TObjArray *arr);
  void DumpTrackInfo(TObjArray *arr);
  void DumpSeedInfo(const AliToyMCTrack *toyTrack, AliTPCseed *seed);
  

  void MarkClustersUsed(AliTPCseed *seed);
  void ResetClustersZtoTime(AliTPCseed *seed);
  
  // reco settings
  Int_t  fSeedingRow;            // first row used for seeding
  Int_t  fSeedingDist;           // distance of seeds
  Int_t  fClusterType;           // cluster type to use
  ECorrType fCorrectionType;     // type of space point correction
  Bool_t fDoTrackFit;            // do track fitting
  Bool_t fUseMaterial;           // use material budget in tracking
  Bool_t fIdealTracking;         // use ideal coordinates for tracking

  Int_t  fNmaxEvents;            // maximum number of events

  // current reconstruction info
  Double_t fTime0;               // current time0 used for reconstruction
  Bool_t   fCreateT0seed;        // if current seed is the T0 seed
  Bool_t   fLongT0seed;          // if we should use a t0 seed including all clusters in the seed range
  
  TTreeSRedirector *fStreamer;   // debug streamer
  TFile *fInputFile;             // input file
  TTree *fTree;                  // input tree with ToyMC info
  AliToyMCEvent *fEvent;         // input event

  AliTPCParam *fTPCParam;            // tpc reco parameters
  AliTPCCorrection *fTPCCorrection; // space charge

  const Int_t fkNSectorInner;        //number of inner sectors
  AliTPCtrackerSector *fInnerSectorArray;  //array of inner sectors
  const Int_t fkNSectorOuter;        //number of outer sectors
  AliTPCtrackerSector *fOuterSectorArray;  //array of outer sectors

  TClonesArray fAllClusters;     //Array keeping all clusters for free seeding

  TExMap fMapTrackEvent;          // map global track number -> event number
  TExMap fMapTrackTrackInEvent;   // map global track number -> track in event

  Bool_t fIsAC;                     // if we are tracking with sector arrays running from 0-36 rather than 0-18
   
  ClassDef(AliToyMCReconstruction,0)
};


#endif
