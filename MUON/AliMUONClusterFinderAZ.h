#ifndef ALIMUONCLUSTERFINDERAZ_H
#define ALIMUONCLUSTERFINDERAZ_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/// \ingroup rec
/// \class AliMUONClusterFinderAZ
/// \brief Cluster finder in MUON arm of ALICE
///
/// \author Alexander Zinchenko, JINR Dubna

#ifndef ALIMUONVCLUSTERFINDER_H
# include "AliMUONVClusterFinder.h"
#endif

class TH2D;
class TClonesArray;
class TMinuit;

class AliMpVSegmentation;
class AliMUONPixel;
//class AliMUONClusterDrawAZ;
#include "TMatrixDfwd.h"
class AliMUONVDigit;
class AliMUONRawCluster;
class AliMUONVDigitStore;
class AliMUONMathieson;
class AliMpVSegmentation;

class AliMUONClusterFinderAZ : public AliMUONVClusterFinder 
{
public:
  AliMUONClusterFinderAZ(Bool_t draw = 0); // Constructor
  virtual ~AliMUONClusterFinderAZ(); // Destructor

  virtual AliMUONCluster* NextCluster();

  virtual Bool_t Prepare(const AliMpVSegmentation* segmentations[2],
                         const AliMUONVDigitStore& digitStore);

  void     FindRawClusters(Int_t ch); // the same interface as for old cluster finder
  void     EventLoop(Int_t ch = 0);
  Bool_t   TestTrack(Int_t t) const; // test if track was selected
  
  /// Return the number of pads in the cluster on the given cathode
  Int_t    GetNPads(Int_t cath) const { return fnPads[cath]; }
  /// Return pad information \todo add more details
  Int_t    GetIJ(Int_t indx, Int_t iPad) const { return fPadIJ[indx][iPad]; }
  /// Return pad information \todo add more details
  Float_t  GetXyq(Int_t indx, Int_t iPad) const { return fXyq[indx][iPad]; }
  /// Return the flag for used pads
  Bool_t GetUsed(Int_t cath, Int_t dig) const { return fUsed[cath][dig]; }

  /// Mark used digits
  void SetUsed(Int_t cath, Int_t dig) { fUsed[cath][dig] = kTRUE; } 
  /// Unmark digits
  void SetUnused(Int_t cath, Int_t dig) { fUsed[cath][dig] = kFALSE; } 
  /// Set reco flag
  void SetReco(Int_t iReco) { fReco = iReco; } 
  /// Start  \todo add more details
  void SetStart(Int_t iCath, Int_t iPad) { fCathBeg = iCath; fPadBeg[0] = fPadBeg[1] = 0; fPadBeg[fCathBeg] = iPad; } 

  void ResetRawClusters();
  
  Float_t ChargeIntegration(Double_t x, Double_t y, 
                            Double_t padX, Double_t padY,
                            Double_t padDX, Double_t padDY);
    
private:
  // Some constants
  static const Int_t fgkDim = 10000; ///< array size
  static const Double_t fgkCouplMin; ///< threshold on coupling 
  static const Double_t fgkZeroSuppression; ///< average zero suppression value
  static const Double_t fgkSaturation; ///< average saturation level

  static  AliMUONClusterFinderAZ* fgClusterFinder; ///< the ClusterFinderAZ instance

  Int_t      fnPads[2];         //!< number of pads in the cluster on 2 cathodes
  Float_t    fXyq[7][fgkDim];   //!< pad information \todo add more details
  UInt_t     fDigitId[fgkDim];  //!< digit id of the pads (to find back the digit using the digitstore)
  Int_t      fPadIJ[4][fgkDim]; //!< pad information \todo add more details
  const AliMpVSegmentation *fSegmentation[2]; //!< new segmentation
  Int_t      fNpar;             //!< number of fit parameters
  Double_t   fQtot;             //!< total cluster charge
  Int_t      fReco;             //!< !=0 if run reco with writing of reconstructed clusters 
  Int_t fCathBeg;               //!< starting cathode (for combined cluster / track reco)
  Int_t fPadBeg[2];             //!< starting pads (for combined cluster / track reco)

  static     TMinuit* fgMinuit; //!< Fitter
  Bool_t     fUsed[2][fgkDim]; //!< flags for used pads
//  AliMUONClusterDrawAZ *fDraw; //!< drawing object 
  TObjArray* fPixArray; //!< collection of pixels
  Int_t fnCoupled; //!< number of coupled clusters in precluster
  Int_t fDebug; //!< debug level

  TClonesArray*           fRawClusters;        //!< array of cluster per ch.
  
  const AliMUONVDigitStore* fDigitStore; //!< digit store we're working on
  
  Int_t fDetElemId; //!< detection element id we're currently working on
  Int_t fChamberId; //!< chamber corresponding the fDetElemId
  
  AliMUONMathieson* fMathieson; //!< mathieson function to be used
  
  Int_t fCurrentCluster; //!< current cluster
  
  // Functions

  /// Not implemented
  AliMUONClusterFinderAZ(const AliMUONClusterFinderAZ& rhs);
  /// Not implemented
  AliMUONClusterFinderAZ& operator=(const AliMUONClusterFinderAZ& rhs);

  void   AddPad(Int_t cath, AliMUONVDigit& digit); // add a pad to the cluster
  Bool_t Overlap(Int_t cath, const AliMUONVDigit& dig); // check if the pad from one cathode overlaps with a pad in the cluster on the other cathode
  Bool_t Overlap(Float_t *xy1, Int_t iPad, Float_t *xy12, Int_t iSkip); // check if pads xy1 and iPad overlap and return overlap area
  Bool_t CheckPrecluster(Int_t *nShown); // check precluster to simplify it (if possible)
  void   BuildPixArray(); // build array of pixels
  void   AdjustPixel(Float_t width, Int_t ixy); // adjust size of small pixels
  void   AdjustPixel(Float_t wxmin, Float_t wymin); // adjust size of large pixels
  Bool_t MainLoop(Int_t iSimple); // repeat MLEM algorithm until pixels become sufficiently small
  void   Mlem(Double_t *coef, Double_t *probi, Int_t nIter); // use MLEM for cluster finding
  void   FindCOG(TH2D *mlem, Double_t *xyc); // find COG position around maximum bin
  Int_t  FindNearest(AliMUONPixel *pixPtr0); // find nearest neighbouring pixel to the given one
  void   Split(TH2D *mlem, Double_t *coef); // steering function for pixels
  void   AddBin(TH2D *mlem, Int_t ic, Int_t jc, Int_t mode, Bool_t* used, TObjArray *pix); // add a bin to the cluster
  TObject* BinToPix(TH2D *mlem, Int_t jc, Int_t ic); // hist. bin-to-pixel
  void   AddCluster(Int_t ic, Int_t nclust, TMatrixD *aijcluclu, Bool_t *used, Int_t *clustNumb, Int_t &nCoupled); // add a cluster to the group of coupled clusters
  Double_t MinGroupCoupl(Int_t nCoupled, Int_t *clustNumb, TMatrixD *aijcluclu, Int_t *minGroup); // find group of cluster with min. coupling to others
  Int_t  SelectPad(Int_t nCoupled, Int_t nForFit, Int_t *clustNumb, Int_t *clustFit, TMatrixD *aijcluclu); //select pads for fit
  void   Merge(Int_t nForFit, Int_t nCoupled, Int_t *clustNumb, Int_t *clustFit, TObjArray **clusters, TMatrixD *aijcluclu, TMatrixD *aijclupad); // merge clusters
  Int_t  Fit(Int_t iSimple, Int_t nfit, Int_t *clustFit, TObjArray **clusters, Double_t *parOk); // do the fitting 
  void  UpdatePads(Int_t nfit, Double_t *par); // subtract fitted charges from pads
  void  AddRawCluster(Double_t x, Double_t y, Double_t qTot, Double_t fmin, Int_t nfit, Int_t *tracks, Double_t sigx, Double_t sigy, Double_t dist); // add new reconstructed cluster
  Int_t FindLocalMaxima(TObjArray *pixArray, Int_t *localMax, Double_t *maxVal); // find local maxima 
  void  FlagLocalMax(TH2D *hist, Int_t i, Int_t j, Int_t *isLocalMax); // flag local max
  void  FindCluster(Int_t *localMax, Int_t iMax); // find cluster around local max
  void  AddVirtualPad(); // add virtual pads for some clusters (if necessary)
  void  PadsInXandY(Int_t &nInX, Int_t &nInY); // get number of pads in X and Y
  // This function is used for fitting
  void  Fcn1(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
  void Simple(); // process simple cluster

  void Errors(AliMUONRawCluster *clus); // correct coordinates and eval. errors
  void Errors(Int_t ny, Int_t nx, Int_t iby, Int_t ibx, Double_t fmin,
	      Double_t wy, Double_t wx, Int_t iover, 
	      Double_t dyc, Double_t dxc, Double_t qtot, 
	      Double_t &yrec, Double_t &xrec, Double_t &erry, Double_t &errx);

  virtual void Print(Option_t* opt="") const;
  void PrintPixel(Int_t i) const;
  void PrintPad(Int_t i) const;
    
  void Used(Int_t indx, Bool_t value);

  /// Dummy method for overloading warnings
  void FindCluster(int, int, int, AliMUONRawCluster&) {return;}
  /// Dummy method for overloading warnings
  void FindLocalMaxima(AliMUONRawCluster*) {return;}
  /// Dummy method for overloading warnings
  void Split(AliMUONRawCluster*) {return;}
  /// Dummy method for overloading warnings
  void AddRawCluster(AliMUONRawCluster&) {return;}

ClassDef(AliMUONClusterFinderAZ,0) // cluster finder in MUON arm of ALICE
};

#endif
