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
class AliTRDcluster;
class AliTRDmcTrack;
class AliBarrelTrack;
class AliESD;

const unsigned kMaxLayersPerSector = 1000;  
const unsigned kMaxTimeBinIndex = 216;  // (30 drift + 6 ampl) * 6 planes  
const unsigned kMaxClusterPerTimeBin = 7000; 
const unsigned kZones = 5; 
const Int_t    kTrackingSectors = 18; 

class AliTRDtracker : public AliTracker { 

 public:

  AliTRDtracker();
  AliTRDtracker(const TFile *in);
  virtual ~AliTRDtracker(); 

  Int_t         Clusters2Tracks(AliESD* event);
  Int_t         PropagateBack(AliESD* event);
  Int_t         RefitInward(AliESD* event);

  Int_t         LoadClusters(TTree *cTree);
  void          UnloadClusters();
  AliCluster   *GetCluster(Int_t index) const { if (index >= fNclusters) return NULL; 
                                                return (AliCluster*) fClusters->UncheckedAt(index); };
  virtual void  CookLabel(AliKalmanTrack *t,Float_t wrong) const;
  virtual void  UseClusters(const AliKalmanTrack *t, Int_t from=0) const;  
  
  void          SetAddTRDseeds() { fAddTRDseeds = kTRUE; }
  void          SetNoTilt() { fNoTilt = kTRUE; }

  Double_t      GetTiltFactor(const AliTRDcluster* c);

  Int_t         ReadClusters(TObjArray *array, TTree *in) const;
  Int_t         CookSectorIndex(Int_t gs) const { return kTrackingSectors - 1 - gs; }
  AliTRDcluster * GetCluster(AliTRDtrack * track, Int_t plane, Int_t timebin);
  Int_t         GetLastPlane(AliTRDtrack * track); //return last updated plane

  Float_t  GetSeedGap()       const {return fgkSeedGap;}   
  Int_t    GetMaxGap()        const {return fMaxGap;}   
  Int_t    GetTimeBinsPerPlane()   const {return fTimeBinsPerPlane;}   
  Float_t  GetSeedStep()      const {return fgkSeedStep;}
  Float_t  GetSeedDepth()     const {return fgkSeedDepth;}
  Float_t  GetSkipDepth()     const {return fgkSkipDepth;}
  Double_t GetMaxChi2()       const {return fgkMaxChi2;}
  Float_t  GetMaxSeedC()      const {return fgkMaxSeedC;}
  Float_t  GetMaxSeedTan()    const {return fgkMaxSeedTan;}
  Double_t GetSeedErrorSY()   const {return fgkSeedErrorSY;}
  Double_t GetSeedErrorSY3()  const {return fgkSeedErrorSY3;}
  Double_t GetSeedErrorSZ()   const {return fgkSeedErrorSZ;}
  Float_t  GetLabelFraction() const {return fgkLabelFraction;}
  Float_t  GetWideRoad()      const {return fgkWideRoad;}

  Float_t  GetMinClustersInTrack() const {return fgkMinClustersInTrack;}
  Float_t  GetMinClustersInSeed()  const {return fgkMinClustersInSeed;} 
  Float_t  GetMaxSeedDeltaZ()      const {return fgkMaxSeedDeltaZ;}
  Float_t  GetMaxSeedVertexZ()     const {return fgkMaxSeedVertexZ;}

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
     void           GetPropagationParameters(Double_t y, Double_t z,
                                Double_t &dx, Double_t &rho, Double_t &x0, 
                                Bool_t &lookForCluster) const;
     Int_t          GetZone( Double_t z) const;
     Int_t          Find(Double_t y) const; 
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
     Float_t  GetTzeroShift() const { return fTzeroShift; }   
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
     Float_t                   fTzeroShift;   // T0 shift in cm
     Int_t                     fGeomSector;   // sector # in AliTRDgeometry
   };

 protected:

  AliTRDgeometry     *fGeom;            // Pointer to TRD geometry
  AliTRDparameter    *fPar;             // Pointer to TRD parameter

  AliTRDtrackingSector *fTrSec[kTrackingSectors];  // array of tracking sectors;    
  Int_t            fNclusters;        // Number of clusters in TRD 
  TObjArray        *fClusters;        // List of clusters for all sectors

  Int_t            fNseeds;           // Number of track seeds  
  TObjArray        *fSeeds;           // List of track seeds
   
  Int_t            fNtracks;          // Number of reconstructed tracks 
  TObjArray        *fTracks;          // List of reconstructed tracks   

  Float_t          fSY2corr;          // Correction coefficient for
                                      // cluster SigmaY2 

  Float_t          fSZ2corr;          // Correction coefficient for
                                      // cluster SigmaZ2 

  static const Float_t  fgkSeedGap;   // Distance between inner and outer
                                      // time bin in seeding 
                                      // (fraction of all time bins) 
  
  static const Float_t  fgkSeedStep;  // Step in iterations
  static const Float_t  fgkSeedDepth; // Fraction of TRD allocated for seeding
  static const Float_t  fgkSkipDepth; // Fraction of TRD which can be skipped
                                      // in track prolongation             
  Int_t       fTimeBinsPerPlane;      // number of sensitive timebins per plane
  Int_t       fMaxGap;                // max gap (in time bins) in the track  
                                      // in track prolongation             

  static const Double_t fgkMaxChi2;   // max increment in track chi2 
        
  static const Float_t  fgkMinClustersInTrack; // min number of clusters in track
                                               // out of total timebins

  static const Float_t  fgkMinFractionOfFoundClusters; // min found clusters 
                                                       // out of expected  

  static const Float_t  fgkMinClustersInSeed;  // min fraction of clusters in seed
  static const Float_t  fgkMaxSeedDeltaZ;   // max dZ in MakeSeeds
  static const Float_t  fgkMaxSeedDeltaZ12; // max abs(z1-z2) in MakeSeeds
  static const Float_t  fgkMaxSeedC;       // max initial curvature in MakeSeeds
  static const Float_t  fgkMaxSeedTan;     // max initial Tangens(lambda) in MakeSeeds
  static const Float_t  fgkMaxSeedVertexZ; // max vertex Z in MakeSeeds
  static const Double_t fgkSeedErrorSY;    // sy parameter in MakeSeeds
  static const Double_t fgkSeedErrorSY3;   // sy3 parameter in MakeSeeds
  static const Double_t fgkSeedErrorSZ;    // sz parameter in MakeSeeds
  static const Float_t  fgkLabelFraction;  // min fraction of same label
  static const Float_t  fgkWideRoad;       // max road width in FindProlongation

  Bool_t                fVocal;            // Whatever...
  Bool_t                fAddTRDseeds;      // Something else

  Bool_t                fNoTilt;           // No tilt, or what?
  Bool_t                fHoles[5][18];     // holes
  
  Bool_t AdjustSector(AliTRDtrack *track); 
 
  
  // Barrel tracks [SR, 03.04.2003]

  static const Int_t fgkFirstPlane;   // Id of the first (innermost) reference plane 
  static const Int_t fgkLastPlane;    // Id of the last (outermost) reference plane
  
  void SetBarrelTree(const char *mode);
  void StoreBarrelTrack(AliTRDtrack *ps, Int_t refPlane, Int_t isIn);
  
  TFile *fBarrelFile;                // Some kind of barrel file
  TTree *fBarrelTree;                // And a barrel tree
  TClonesArray *fBarrelArray;        // Wow, there even an array for that barrel
  AliBarrelTrack *fBarrelTrack;      // And, finally, the track

 private:

  virtual void  MakeSeeds(Int_t inner, Int_t outer, Int_t turn);

  Int_t         FollowProlongation(AliTRDtrack& t, Int_t rf);
  Int_t         FollowBackProlongation(AliTRDtrack& t);
  Int_t         Refit(AliTRDtrack& t, Int_t rf);

  Int_t         PropagateToTPC(AliTRDtrack& t);
  Int_t         PropagateToOuterPlane(AliTRDtrack& t, Double_t x);

  void          SetSY2corr(Float_t w)    {fSY2corr = w;}
  void          SetSZ2corr(Float_t w)    {fSZ2corr = w;}
  Double_t      ExpectedSigmaY2(Double_t r, Double_t tgl, Double_t pt) const;
  Double_t      ExpectedSigmaZ2(Double_t r, Double_t tgl) const;

  ClassDef(AliTRDtracker,1)           // manager base class  

};

#endif 
