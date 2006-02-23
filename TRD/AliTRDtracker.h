#ifndef ALITRDTRACKER_H
#define ALITRDTRACKER_H   

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */ 

#include "AliTracker.h" 
#include "TObjArray.h" 

class TFile;
class TTree;
class TParticle;
class TParticlePDG;

class AliTRDgeometry;
class AliTRDparameter;
class AliTRDtrack;
class AliTRDtracklet;
class AliTRDcluster;
class AliESD;
class TTreeSRedirector;

const unsigned kMaxLayersPerSector = 1000;  
const unsigned kMaxTimeBinIndex = 216;  // (30 drift + 6 ampl) * 6 planes  
const unsigned kMaxClusterPerTimeBin = 7000; 
const unsigned kZones = 5; 
const Int_t    kTrackingSectors = 18; 

class AliTRDseed : public TObject{
  friend class AliTRDtracker;
 public:
  AliTRDseed();                 // default constructor
  ~AliTRDseed(){};              // default constructor
  static  void EvaluateUni(Int_t nvectors, Double_t *data, Double_t &mean, Double_t &sigma, Int_t hh);
  static  Float_t   FitRiemanTilt(AliTRDseed * seed, Bool_t error);
  void           UseClusters(); // use clusters
  void           Update();      // update information - without tilt correction
  void           CookLabels();  // cook label
  void           UpdateUsed();
  void           Reset();       // reset seed
  Bool_t         isOK(){return fN2>8;}
 private:
  Float_t        fTilt;         // tilting angle
  Float_t        fPadLength;    // pad length
  Float_t        fX0;           // x0 position
  Float_t        fX[25];        // !x position
  Float_t        fY[25];        // !y position
  Float_t        fZ[25];        // !z position
  Int_t          fIndexes[25];  // !indexes
  AliTRDcluster *fClusters[25]; // !clusters
  Bool_t         fUsable[25];   // !indication  - usable cluster
  Float_t        fYref[2];      // reference y
  Float_t        fZref[2];      // reference z
  Float_t        fYfit[2];      // y fit position +derivation
  Float_t        fYfitR[2];      // y fit position +derivation
  Float_t        fZfit[2];      // z fit position
  Float_t        fZfitR[2];      // z fit position
  Float_t        fSigmaY;        // "robust" sigma in Y - constant fit
  Float_t        fSigmaY2;       // "robust" sigma in Y - line fit
  Float_t        fMeanz;         // mean vaue of z
  Float_t        fZProb;         // max probbable z
  Int_t          fLabels[2];    // labels
  Int_t          fN;            // number of associated clusters
  Int_t          fN2;            // number of not crossed
  Int_t          fNUsed;        // number of used clusters
  Int_t          fFreq;         // freq
  Int_t          fNChange;      // change z counter
  Float_t        fMPads;        // mean number of pads per cluster
  // global
  //
  Float_t        fC;            // curvature
  Float_t        fCC;           // curvature with constrain
  Float_t        fChi2;         // global chi2
  Float_t        fChi2Z;        // global chi2
 private:
  ClassDef(AliTRDseed,1)  
};


class AliTRDtracker : public AliTracker { 

 public:


  AliTRDtracker();
  AliTRDtracker(const TFile *in);
  virtual ~AliTRDtracker(); 
  static        Int_t  Freq(Int_t n, const Int_t *inlist, Int_t *outlist, Bool_t down);    
  Int_t         Clusters2Tracks(AliESD* event);
  Int_t         PropagateBack(AliESD* event);
  Int_t         RefitInward(AliESD* event);
  
  Int_t         LoadClusters(TTree *cTree);
  void          UnloadClusters();
  AliCluster   *GetCluster(Int_t index) const { if (index >= fNclusters) return NULL; 
                                                return (AliCluster*) fClusters->UncheckedAt(index); };
  Bool_t        GetTrackPoint(Int_t index, AliTrackPoint& p) const;
  virtual void  CookLabel(AliKalmanTrack *t,Float_t wrong) const;
  virtual void  UseClusters(const AliKalmanTrack *t, Int_t from=0) const;  
  
  void          SetAddTRDseeds() { fAddTRDseeds = kTRUE; }
  void          SetNoTilt() { fNoTilt = kTRUE; }

  Double_t      GetTiltFactor(const AliTRDcluster* c);

  Int_t         ReadClusters(TObjArray *array, TTree *in) const;
  Int_t         CookSectorIndex(Int_t gs) const { return kTrackingSectors - 1 - gs; }
  AliTRDcluster * GetCluster(AliTRDtrack * track, Int_t plane, Int_t timebin, UInt_t &index);
  Int_t         GetLastPlane(AliTRDtrack * track); //return last updated plane
  Int_t FindClusters(Int_t sector, Int_t t0, Int_t t1, AliTRDtrack * track, Int_t *clusters, AliTRDtracklet& tracklet);

  Int_t    GetTimeBinsPerPlane()   const {return fTimeBinsPerPlane;}   
  Double_t GetMaxChi2()       const {return fgkMaxChi2;}
  Float_t  GetLabelFraction() const {return fgkLabelFraction;}

  Float_t  GetMinClustersInTrack() const {return fgkMinClustersInTrack;}

  // x <-> timebin conversions useful in analysis macros
  Double_t GetX(Int_t sec, Int_t plane, Int_t localTB) const;
  Double_t GetX(Int_t sec, Int_t pl) const { 
    return fTrSec[sec]->GetLayer(pl)->GetX(); }
  Int_t GetGlobalTimeBin(Int_t sec, Int_t plane, Int_t localTB) const {
    return fTrSec[sec]->CookTimeBinIndex(plane,localTB); }
  Double_t GetLayerNumber(Int_t sec, Double_t x) const {
    return fTrSec[sec]->GetLayerNumber(x); }

  class AliTRDpropagationLayer {
   // *****************  internal class *******************
   public: 
     AliTRDpropagationLayer(Double_t x, Double_t dx, Double_t rho, 
                            Double_t x0, Int_t tbIndex); 

     ~AliTRDpropagationLayer() { 
       if(fTimeBinIndex >= 0) { delete[] fClusters; delete[] fIndex; }
     }
     void InsertCluster(AliTRDcluster *c, UInt_t index);
     operator       Int_t() const {return fN;}
     AliTRDcluster* operator[](Int_t i) {return fClusters[i];}
     UInt_t         GetIndex(Int_t i) const {return fIndex[i];} 
     Double_t       GetX() const { return fX; }
     Double_t       GetdX() const { return fdX; }
     Double_t       GetRho() const { return fRho; }
     Double_t       GetX0() const { return fX0; }
     Int_t          GetTimeBinIndex() const { return fTimeBinIndex; }     
     Bool_t         GetPropagationParameters(Double_t y, Double_t z,
                                Double_t &dx, Double_t &rho, Double_t &x0, 
					     Bool_t &lookForCluster) const;
     Int_t          GetZone( Double_t z) const;
     Int_t          Find(Float_t y) const; 
     Int_t          FindNearestCluster(Float_t y, Float_t z, Float_t maxroad, Float_t maxroadz) const;
     void           SetZmax(Int_t cham, Double_t center, Double_t w)
                      { fZc[cham] = center;  fZmax[cham] = w; }
     void           SetZ(Double_t* center, Double_t *w, Double_t *wsensitive);
     void           SetHoles(Bool_t* holes);
     void           SetYmax(Double_t w, Double_t wsensitive) { fYmax = w; fYmaxSensitive = wsensitive; }
     Double_t       GetYmax() const { return fYmax; }
     Double_t       GetZmax(Int_t c) const { return fZmax[c]; }
     Double_t       GetZc(Int_t c) const { return fZc[c]; }
     
     void           SetHole(Double_t Zmax, Double_t Ymax,
                            Double_t rho = 1.29e-3, Double_t x0 = 36.66,
                            Double_t Yc = 0, Double_t Zc = 0);
                            
     Bool_t         IsSensitive() const {return (fTimeBinIndex>=0)? kTRUE: kFALSE;}
                       
     void    Clear() {for(Int_t i=0; i<fN; i++) fClusters[i] = NULL; fN = 0;}
     Bool_t  IsHole(Int_t zone) const  { return fIsHole[zone];}              
   private:     

     Int_t         fN;          // this is fN
     Int_t         fSec;        // sector mumber
     AliTRDcluster **fClusters; // array of pointers to clusters
     UInt_t        *fIndex;     // array of cluster indexes
     Double_t       fX;         // x coordinate of the middle plane
     Double_t       fdX;        // radial thickness of the time bin
     Double_t       fRho;       // default density of the material 
     Double_t       fX0;        // default radiation length 
     Int_t          fTimeBinIndex;  // plane * F(local_tb)  
     Double_t       fZc[kZones];  // Z position of the center for 5 active areas
     Double_t       fZmax[kZones]; // half of active area length in Z
     Double_t       fZmaxSensitive[kZones]; //sensitive area for detection Z     
     Bool_t         fIsHole[kZones]; //is hole in given sector       
     Double_t       fYmax;        // half of active area length in Y
     Double_t       fYmaxSensitive;        // half of active area length in Y

     Bool_t         fHole;        // kTRUE if there is a hole in the layer
     Double_t       fHoleZc;      // Z of the center of the hole 
     Double_t       fHoleZmax;    // half of the hole length in Z
     Double_t       fHoleYc;      // Y of the center of the hole 
     Double_t       fHoleYmax;    // half of the hole length in Y 
     Double_t       fHoleRho;     // density of the gas in the hole 
     Double_t       fHoleX0;      // radiation length of the gas in the hole 
   };

   class AliTRDtrackingSector {
   public:
     AliTRDtrackingSector(AliTRDgeometry* geo, Int_t gs, AliTRDparameter* par);
     ~AliTRDtrackingSector() { for(Int_t i=0; i<fN; i++) delete fLayers[i]; }
     Int_t    GetNumberOfLayers() const { return fN; }
     Int_t    GetNumberOfTimeBins() const;
     Double_t GetX(Int_t pl) const { return fLayers[pl]->GetX(); }
     void     MapTimeBinLayers();
     Int_t    GetLayerNumber(Double_t x) const;
     Int_t    GetInnerTimeBin() const;
     Int_t    GetOuterTimeBin() const;
     Int_t    GetLayerNumber(Int_t tb) const {return fTimeBinIndex[tb];}
     Int_t    Find(Double_t x) const; 
     void     InsertLayer(AliTRDpropagationLayer* pl);
     //     AliTRDpropagationLayer* operator[](Int_t i) { return fLayers[i]; }
     AliTRDpropagationLayer* GetLayer(Int_t i) { return fLayers[i]; }
     Int_t    CookTimeBinIndex(Int_t plane, Int_t localTB) const;     

   private:
     Int_t                     fN;      // total number of layers
     AliTRDgeometry            *fGeom;  // geometry
     AliTRDparameter           *fPar;   // parameter
     AliTRDpropagationLayer    *fLayers[kMaxLayersPerSector];   // layers   
     Int_t                     fTimeBinIndex[kMaxTimeBinIndex]; // time bin index
     Int_t                     fGeomSector;   // sector # in AliTRDgeometry
   };

 protected:

  friend class AliTRDtracker::AliTRDtrackingSector;

  AliTRDgeometry     *fGeom;            // Pointer to TRD geometry
  AliTRDparameter    *fPar;             // Pointer to TRD parameter

  AliTRDtrackingSector *fTrSec[kTrackingSectors];  // array of tracking sectors;    
  Int_t            fNclusters;        // Number of clusters in TRD 
  TObjArray        *fClusters;        // List of clusters for all sectors

  Int_t            fNseeds;           // Number of track seeds  
  TObjArray        *fSeeds;           // List of track seeds
   
  Int_t            fNtracks;          // Number of reconstructed tracks 
  TObjArray        *fTracks;          // List of reconstructed tracks   


  Int_t       fTimeBinsPerPlane;      // number of sensitive timebins per plane
                                      // in track prolongation             

  static const Double_t fgkMaxChi2;   // max increment in track chi2 
        
  static const Float_t  fgkMinClustersInTrack; // min number of clusters in track
                                               // out of total timebins

  static const Float_t  fgkLabelFraction;  // min fraction of same label
  //
  static const Double_t  fgkDriftCorrection; // correction coeficients for drift velocity
  static const Double_t  fgkOffset;         // correction coeficients 
  static const Double_t  fgkOffsetX;         // correction coeficients offset in X 
  static const Double_t  fgkCoef;           // correction coeficients
  static const Double_t  fgkMean;           // correction coeficients
  static const Double_t  fgkExB;            // correction coeficients

  Bool_t                fAddTRDseeds;      // Something else

  Bool_t                fNoTilt;           // No tilt, or what?
  Bool_t                fHoles[5][18];     // holes
  
  Bool_t AdjustSector(AliTRDtrack *track); 
 
 private:
  AliTRDtrack *   RegisterSeed(AliTRDseed * seeds, Double_t *params);
  void  MakeSeedsMI(Int_t inner, Int_t outer, AliESD *esd=0);

  Int_t         FollowProlongation(AliTRDtrack& t, Int_t rf);
  Int_t         FollowBackProlongation(AliTRDtrack& t);
  Int_t         Refit(AliTRDtrack& t, Int_t rf);
  void          CookdEdxTimBin(AliTRDtrack& t);  

  Int_t         PropagateToTPC(AliTRDtrack& t);
  Int_t         PropagateToOuterPlane(AliTRDtrack& t, Double_t x);

  Double_t      ExpectedSigmaY2(Double_t r, Double_t tgl, Double_t pt) const;
  Double_t      ExpectedSigmaZ2(Double_t r, Double_t tgl) const;
  TTreeSRedirector *fDebugStreamer;     //!debug streamer
  ClassDef(AliTRDtracker,2)           // manager base class  

};

#endif 
