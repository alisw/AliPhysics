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
class TMinuit;

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
  AliMUONClusterFinderMLEM(Bool_t plot, AliMUONVClusterFinder* clusterFinder); // Constructor
  virtual ~AliMUONClusterFinderMLEM(); // Destructor

  virtual Bool_t NeedSegmentation() const { return kTRUE; }
  
  using AliMUONVClusterFinder::Prepare;

  virtual Bool_t Prepare(Int_t detElemId,
                         TObjArray* pads[2],
                         const AliMpArea& area,
                         const AliMpVSegmentation* segmentations[2]);
  
  virtual AliMUONCluster* NextCluster();
  
  virtual void SetChargeHints(Double_t lowestPadCharge, Double_t lowestClusterCharge);
  
  virtual void Print(Option_t* opt="") const;

  virtual void Paint(Option_t* opt="");

  // Status flags for pads

               /// Return pad "basic" state flag
  static Int_t GetZeroFlag()       { return fgkZero; }
               /// Return do not kill flag
  static Int_t GetMustKeepFlag()   { return fgkMustKeep; }
               /// Return should be used for fit flag
  static Int_t GetUseForFitFlag()  { return fgkUseForFit; }
               /// Return processing is over flag
  static Int_t GetOverFlag()       { return fgkOver; }
               /// Return modified pad charge flag
  static Int_t GetModifiedFlag()   { return fgkModified; }
               /// Return coupled pad flag
  static Int_t GetCoupledFlag()    { return fgkCoupled; }
  
private:
  /// Not implemented
  AliMUONClusterFinderMLEM(const AliMUONClusterFinderMLEM& rhs);
  /// Not implemented
  AliMUONClusterFinderMLEM& operator=(const AliMUONClusterFinderMLEM& rhs);

  Bool_t WorkOnPreCluster();

  /// Check precluster to simplify it (if possible), and return the simplified cluster
  AliMUONCluster* CheckPrecluster(const AliMUONCluster& cluster); 
  AliMUONCluster* CheckPreclusterTwoCathodes(AliMUONCluster* cluster); 
  
  /// Checks whether a pad and a pixel have an overlapping area.
  Bool_t Overlap(const AliMUONPad& pad, const AliMUONPad& pixel); 
  
  /// build array of pixels
  void BuildPixArray(AliMUONCluster& cluster); 
  void BuildPixArrayOneCathode(AliMUONCluster& cluster); 
  void PadOverHist(Int_t idir, Int_t ix0, Int_t iy0, AliMUONPad *pad,
		   TH2D *hist1, TH2D *hist2);

  void RemovePixel(Int_t i);
  
  AliMUONPad* Pixel(Int_t i) const;
  
  Bool_t MainLoop(AliMUONCluster& cluster, Int_t iSimple); // repeat MLEM algorithm until pixels become sufficiently small
  
  void   Mlem(AliMUONCluster& cluster, const Double_t *coef, Double_t *probi, Int_t nIter); // use MLEM for cluster finding
  
  void   FindCOG(Double_t *xyc); // find COG position around maximum bin
  Int_t  FindNearest(const AliMUONPad *pixPtr0); // find nearest neighbouring pixel to the given one

  Int_t FindLocalMaxima(TObjArray *pixArray, Int_t *localMax, Double_t *maxVal); // find local maxima 
  void  FlagLocalMax(TH2D *hist, Int_t i, Int_t j, Int_t *isLocalMax); // flag local max
  void  FindCluster(AliMUONCluster& cluster, const Int_t *localMax, Int_t iMax); // find cluster around local max
  void  AddVirtualPad(AliMUONCluster& cluster); // add virtual pads for some clusters (if necessary)
  
  void  PadsInXandY(AliMUONCluster& cluster, Int_t &nInX, Int_t &nInY) const; // get number of pads in X and Y

  /// Process simple cluster
  void Simple(AliMUONCluster& cluster); 
  
  void Plot(const char* outputfile);
    
  void ComputeCoefficients(AliMUONCluster& cluster, 
                           Double_t* coef, Double_t* probi);
  
  void CheckOverlaps();
  void AddBinSimple(TH2D *mlem, Int_t ic, Int_t jc);
  void MaskPeaks(Int_t mask);

private:
  // Status flags for pads
  static const Int_t fgkZero; ///< pad "basic" state
  static const Int_t fgkMustKeep; ///< do not kill (for pixels)
  static const Int_t fgkUseForFit; ///< should be used for fit
  static const Int_t fgkOver; ///< processing is over
  static const Int_t fgkModified; ///< modified pad charge 
  static const Int_t fgkCoupled; ///< coupled pad  
      
  // Some constants
  static const Double_t fgkDistancePrecision; ///< used to check overlaps and so on
  static const TVector2 fgkIncreaseSize; ///< idem
  static const TVector2 fgkDecreaseSize; ///< idem
  
  AliMUONVClusterFinder* fPreClusterFinder; //!<! the pre-clustering worker
  AliMUONCluster* fPreCluster; //!<! current pre-cluster
  TObjArray fClusterList; //!<! clusters corresponding to the current pre-cluster
  
  Int_t fEventNumber; //!<! current event being processed
  Int_t fDetElemId; //!<! current DE being processed
  Int_t fClusterNumber; //!<! current cluster number
  
  const AliMpVSegmentation *fkSegmentation[2]; //!<! new segmentation
  
  //Int_t fCathBeg;               //!<! starting cathode (for combined cluster / track reco)
  //Int_t fPadBeg[2];             //!<! starting pads (for combined cluster / track reco)
  
  //static     TMinuit* fgMinuit; //!<! Fitter
  TH2D *fHistMlem; //!<! histogram for MLEM procedure
  TH2D *fHistAnode; //!<! histogram for local maxima search
  
  TObjArray* fPixArray; //!<! collection of pixels
  Int_t fDebug; //!<! debug level
  Bool_t fPlot; //!<! whether we should plot thing (for debug only, quite slow!)
  
  AliMUONClusterSplitterMLEM* fSplitter; //!<! helper class to go from pixel arrays to clusters
  Int_t fNClusters; //!<! total number of clusters
  Int_t fNAddVirtualPads; //!<! number of clusters for which we added virtual pads
  
  Double_t fLowestPixelCharge; //!<! see AliMUONRecoParam
  Double_t fLowestPadCharge; //!<! see AliMUONRecoParam
  Double_t fLowestClusterCharge; //!<! see AliMUONRecoParam
  
  ClassDef(AliMUONClusterFinderMLEM,0) // cluster finder in MUON arm of ALICE
};

#endif
