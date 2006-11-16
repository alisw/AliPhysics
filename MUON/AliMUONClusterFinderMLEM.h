#ifndef ALIMUONCLUSTERFINDERMLEM_H
#define ALIMUONCLUSTERFINDERMLEM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/// \ingroup rec
/// \class AliMUONClusterFinderMLEM
/// \brief Cluster finder in MUON arm of ALICE
///
//  Author Alexander Zinchenko, JINR Dubna; Laurent Aphecetche, SUBATECH
//

class TH2D;
class TClonesArray;
class TMinuit;
class TStopwatch;

#ifndef ROOT_TObjArray
#  include "TObjArray.h"
#endif
#ifndef ROOT_TVector2
#  include "TVector2.h"
#endif

class AliMUONPad;

#include "AliMUONVClusterFinder.h"

class AliMUONClusterSplitterMLEM;

class AliMUONClusterFinderMLEM : public AliMUONVClusterFinder
{
public:
  AliMUONClusterFinderMLEM(Bool_t plot=kFALSE); // Constructor
  virtual ~AliMUONClusterFinderMLEM(); // Destructor

  virtual Bool_t Prepare(const AliMpVSegmentation* segmentations[2],
                         TClonesArray* digits[2]);
  
  virtual AliMUONCluster* NextCluster();
  
  virtual void Print(Option_t* opt="") const;

  virtual void Paint(Option_t* opt="");

private:
  AliMUONClusterFinderMLEM(const AliMUONClusterFinderMLEM& rhs);
  AliMUONClusterFinderMLEM& operator=(const AliMUONClusterFinderMLEM& rhs);

  Bool_t WorkOnPreCluster();

  /// Check precluster to simplify it (if possible), and return the simplified cluster
  AliMUONCluster* CheckPrecluster(const AliMUONCluster& cluster); 
  AliMUONCluster* CheckPreclusterTwoCathodes(AliMUONCluster* cluster); 
  AliMUONCluster* CheckPreclusterOneCathode(AliMUONCluster* cluster); 
  
  /// Checks whether a pad and a pixel have an overlapping area.
  Bool_t Overlap(const AliMUONPad& pad, const AliMUONPad& pixel); 
  
  /// build array of pixels
  void BuildPixArray(AliMUONCluster& cluster); 
  void BuildPixArrayOneCathode(AliMUONCluster& cluster); 
  void BuildPixArrayTwoCathodes(AliMUONCluster& cluster); 

  void RemovePixel(Int_t i);
  
  AliMUONPad* Pixel(Int_t i) const;
  
  void   AdjustPixel(AliMUONCluster& cluster, Float_t width, Int_t ixy); // adjust size of small pixels
  void   AdjustPixel(Double_t wxmin, Double_t wymin); // adjust size of large pixels

  Bool_t MainLoop(AliMUONCluster& cluster, Int_t iSimple); // repeat MLEM algorithm until pixels become sufficiently small
  
  void   Mlem(AliMUONCluster& cluster, Double_t *coef, Double_t *probi, Int_t nIter); // use MLEM for cluster finding
  
  void   FindCOG(TH2D *mlem, Double_t *xyc); // find COG position around maximum bin
  Int_t  FindNearest(AliMUONPad *pixPtr0); // find nearest neighbouring pixel to the given one

  Int_t FindLocalMaxima(TObjArray *pixArray, Int_t *localMax, Double_t *maxVal); // find local maxima 
  void  FlagLocalMax(TH2D *hist, Int_t i, Int_t j, Int_t *isLocalMax); // flag local max
  void  FindCluster(AliMUONCluster& cluster, Int_t *localMax, Int_t iMax); // find cluster around local max
  void  AddVirtualPad(AliMUONCluster& cluster); // add virtual pads for some clusters (if necessary)
  
  void  PadsInXandY(AliMUONCluster& cluster, Int_t &nInX, Int_t &nInY) const; // get number of pads in X and Y

  /// Process simple cluster
  void Simple(AliMUONCluster& cluster); 
  
  void Neighbours(Int_t cath, Int_t ix0, Int_t iy0, 
                  Int_t& nn,
                  Int_t* xList, Int_t* yList);

  void Plot(const char* outputfile);
    
  void ComputeCoefficients(AliMUONCluster& cluster, 
                           Double_t* coef, Double_t* probi);
  
  void CheckOverlaps();
  
  TStopwatch* Timer(Int_t i) const;
  
private:
    
  // Some constants
  static const Int_t fgkDim = 10000; ///< array size
  static const Double_t fgkZeroSuppression; ///< average zero suppression value
  static const Double_t fgkSaturation; ///< average saturation level
  static const Double_t fgkDistancePrecision; ///< used to check overlaps and so on
  static const TVector2 fgkIncreaseSize; ///< idem
  static const TVector2 fgkDecreaseSize; ///< idem
  
  AliMUONVClusterFinder* fPreClusterFinder; ///!< the pre-clustering worker
  AliMUONCluster* fPreCluster; ///<! current pre-cluster
  TObjArray fClusterList; ///!< clusters corresponding to the current pre-cluster
  
  Int_t fEventNumber; ///!< current event being processed
  Int_t fDetElemId; ///!< current DE being processed
  Int_t fClusterNumber; ///!< current cluster number
  
  const AliMpVSegmentation *fSegmentation[2]; //!< new segmentation
  
  Float_t    fZpad;             //!< z-coordinate of the hit
  Int_t      fReco;             //!< !=0 if run reco with writing of reconstructed clusters 
  Int_t fCathBeg;               //!< starting cathode (for combined cluster / track reco)
  Int_t fPadBeg[2];             //!< starting pads (for combined cluster / track reco)
  
  static     TMinuit* fgMinuit; //!< Fitter
  
  TObjArray* fPixArray; //!< collection of pixels
  Int_t fDebug; //!< debug level
  Bool_t fPlot; //!< whether we should plot thing (for debug only, quite slow!)
  
  TObjArray* fTimers; //!< internal timers
  enum ETimer { kMainLoop, kCheckPreCluster, kLast };
  
  AliMUONClusterSplitterMLEM* fSplitter; //!< helper class to go from pixel arrays to clusters
  Int_t fNClusters; //!< total number of clusters
  Int_t fNAddVirtualPads; //!< number of clusters for which we added virtual pads
  
  ClassDef(AliMUONClusterFinderMLEM,0) // cluster finder in MUON arm of ALICE
};

#endif
