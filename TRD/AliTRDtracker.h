#ifndef ALITRDTRACKER_H
#define ALITRDTRACKER_H   

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */ 

/* $Id$ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  The TRD tracker                                                       //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#ifndef ROOT_TObjArray
#include "TObjArray.h" 
#endif

#ifndef ALITRACKER_H
#include "AliTracker.h" 
#endif

#ifndef ALITRDPROPAGATIONLAYER_H
#include "AliTRDpropagationLayer.h"
#endif

class TFile;
class TTree;
class TH1D;
class TH2D;
class TParticle;
class TParticlePDG;
class TTreeSRedirector;

class AliTRDgeometry;
class AliTRDtrack;
class AliTRDtracklet;
class AliTRDcluster;
class AliTRDseed;
class AliESDEvent;
class AliTRDpropagationLayer;
class AliTRDReconstructor;

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  The standard TRD tracker                                                 //  
//  Based on Kalman filltering approach                                      //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

class AliTRDtracker : public AliTracker { 

 public:

  void InitLogHists();
  void SaveLogHists();

  enum { kMaxLayersPerSector   = 1000
       , kMaxTimeBinIndex      = 216
       , kTrackingSectors      = 18   };
  
  AliTRDtracker(AliTRDReconstructor *rec = 0x0);
  AliTRDtracker(const AliTRDtracker &t);
  AliTRDtracker(const TFile *in, AliTRDReconstructor *rec = 0x0);
  virtual         ~AliTRDtracker(); 
  AliTRDtracker   &operator=(const AliTRDtracker &/*t*/) { return *this;          } 
  
  void             SetReconstructor(AliTRDReconstructor *rec) {fReconstructor = rec;}
  void             SetAddTRDseeds()               { fAddTRDseeds = kTRUE;         }
  void             SetNoTilt()                    { fNoTilt      = kTRUE;         }
  
  Int_t            GetTimeBinsPerPlane() const    { return fTimeBinsPerPlane;     }   
  Double_t         GetMaxChi2() const             { return fgkMaxChi2;            }
  Float_t          GetLabelFraction() const       { return fgkLabelFraction;      }
  Float_t          GetMinClustersInTrack() const  { return fgkMinClustersInTrack; }
  Int_t            GetLastPlane(AliTRDtrack *track);
  Double_t         GetTiltFactor(const AliTRDcluster *c);
  virtual Bool_t   GetTrackPoint(Int_t index, AliTrackPoint& p) const;
  Double_t         GetX(Int_t sec, Int_t plane, Int_t localTB) const;
  Double_t         GetX(Int_t sec, Int_t pl) const
                                                  { return fTrSec[sec]->GetLayer(pl)->GetX();            }
  Int_t            GetGlobalTimeBin(Int_t sec, Int_t plane, Int_t localTB) const 
                                                  { return fTrSec[sec]->CookTimeBinIndex(plane,localTB); }
  Double_t         GetLayerNumber(Int_t sec, Double_t x) const 
                                                  { return fTrSec[sec]->GetLayerNumber(x);               }
  AliCluster      *GetCluster(Int_t index) const  { if (index >= fNclusters) return NULL; 
                                                    return (AliCluster *) fClusters->UncheckedAt(index); }
  
  static Int_t     Freq(Int_t n, const Int_t *inlist, Int_t *outlist, Bool_t down);    
  Int_t            Clusters2Tracks(AliESDEvent *event);
  Int_t            PropagateBack(AliESDEvent *event);
  Int_t            RefitInward(AliESDEvent *event);
  
  virtual void     CookLabel(AliKalmanTrack *t, Float_t wrong) const;
  
  Int_t            LocalToGlobalID(Int_t lid);
  Int_t            GlobalToLocalID(Int_t gid);
  
  Int_t            LoadClusters(TTree *cTree);
  void             UnloadClusters();
  virtual void     UseClusters(const AliKalmanTrack *t, Int_t from = 0) const;  
  Int_t            ReadClusters(TObjArray *array, TTree *in) const;
  AliTRDcluster   *GetCluster(AliTRDtrack *track, Int_t plane, Int_t timebin, UInt_t &index);
  Int_t            FindClusters(Int_t sector, Int_t t0, Int_t t1, AliTRDtrack *track
                              , Int_t *clusters, AliTRDtracklet &tracklet);  
  
 protected:

  Bool_t           AdjustSector(AliTRDtrack *track); 
  AliTRDtrack     *RegisterSeed(AliTRDseed *seeds, Double_t *params);
  Int_t            FollowBackProlongation(AliTRDtrack &t);
  //void             MakeSeedsMI(Int_t inner, Int_t outer, AliESDEvent *esd = 0);
  
 protected:

  //__________________________________________________________________________________________________
  class AliTRDtrackingSector {
    
   public:
    
    AliTRDtrackingSector(AliTRDgeometry* geo, Int_t gs);
    AliTRDtrackingSector(const AliTRDtrackingSector &/*t*/);
    ~AliTRDtrackingSector();
    
    AliTRDtrackingSector &operator=(const AliTRDtrackingSector &/*t*/) { return *this; }
    
    Int_t    GetNumberOfLayers() const             { return fN; }
    Int_t    GetNumberOfTimeBins() const;
    Int_t    GetLayerNumber(Double_t x) const;
    Int_t    GetInnerTimeBin() const;
    Int_t    GetOuterTimeBin() const;
    Int_t    GetLayerNumber(Int_t tb) const        { return fTimeBinIndex[tb];   }
    Double_t GetX(Int_t pl) const                  { return fLayers[pl]->GetX(); }
    AliTRDpropagationLayer* GetLayer(Int_t i)      { return fLayers[i];          }
    Int_t    GetSector() const {return fGeomSector;}	
    
    void     MapTimeBinLayers();
    Int_t    Find(Double_t x) const; 
    void     InsertLayer(AliTRDpropagationLayer *pl);
    Int_t    CookTimeBinIndex(Int_t plane, Int_t localTB) const;     
    
   private:
    
    Int_t                   fN;                              // Total number of layers
    AliTRDgeometry         *fGeom;                           // Geometry
    AliTRDpropagationLayer *fLayers[kMaxLayersPerSector];    // Layers   
    Int_t                   fTimeBinIndex[kMaxTimeBinIndex]; // Time bin index
    Int_t                   fGeomSector;                     // Sector# in AliTRDgeometry
    
  };
  
 protected:
  AliTRDReconstructor     *fReconstructor;
  AliTRDgeometry          *fGeom;                          // Pointer to TRD geometry
  AliTRDtrackingSector    *fTrSec[kTrackingSectors];       // Array of tracking sectors;    
  Int_t                    fNclusters;                     // Number of clusters in TRD 
  TObjArray               *fClusters;                      // List of clusters for all sectors
  Int_t                    fNseeds;                        // Number of track seeds  
  TObjArray               *fSeeds;                         // List of track seeds
  Int_t                    fNtracks;                       // Number of reconstructed tracks 
  TObjArray               *fTracks;                        // List of reconstructed tracks   
  Int_t                    fTimeBinsPerPlane;              // Timebins per plane in track prolongation 
  
  static const Double_t    fgkMaxChi2;                     // Max increment in track chi2 
  static const Float_t     fgkMinClustersInTrack;          // Min number of clusters in track
  static const Float_t     fgkLabelFraction;               // Min fraction of same label
  static const Double_t    fgkMaxSnp;                      // Maximal snp for tracking
  static const Double_t    fgkMaxStep;                     // Maximal step for tracking  
  
  Bool_t                   fAddTRDseeds;                   // Something else
  Bool_t                   fNoTilt;                        // No tilt, or what?

  // Histograms
  TH1D                    *fHBackfit;                      // QA histogram
  TH1D                    *fHClSearch;                     // QA histogram
  TH1D                    *fHRefit;                        // QA histogram
  TH1D                    *fHX;                            // QA histogram
  TH1D                    *fHNCl;                          // QA histogram
  TH1D                    *fHNClTrack;                     // QA histogram
  TH1D                    *fHFindCl[4];                    // QA histogram
  TH1D                    *fHMinYPos;                      // QA histogram
  TH1D                    *fHMinYNeg;                      // QA histogram
  TH1D                    *fHMinZ;                         // QA histogram
  TH2D                    *fHMinD;                         // QA histogram
  TH1D                    *fHDeltaX;                       // QA histogram
  TH1D                    *fHXCl;                          // QA histogram

 protected:

  Int_t            FollowProlongation(AliTRDtrack &t);
  Int_t            PropagateToX(AliTRDtrack &t, Double_t xToGo, Double_t maxStep);
  Double_t         ExpectedSigmaY2(Double_t r, Double_t tgl, Double_t pt) const;
  Double_t         ExpectedSigmaZ2(Double_t r, Double_t tgl) const;

 
  TTreeSRedirector        *fDebugStreamer;                 //!Debug streamer
  
  ClassDef(AliTRDtracker, 4)                               // TRD tracker
    
};


#endif 
