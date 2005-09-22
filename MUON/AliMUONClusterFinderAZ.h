#ifndef ALIMUONCLUSTERFINDERAZ_H
#define ALIMUONCLUSTERFINDERAZ_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// Revision of includes 07/05/2004

/// \ingroup rec
/// \class AliMUONClusterFinderAZ
/// \brief Cluster finder in MUON arm of ALICE

#include "AliMUONClusterFinderVS.h"

class TH2F;
class TH2D;
class TClonesArray;
class TMinuit;
class TMatrixD;

class AliSegmentation;
class AliMUONResponse;
class AliMUONPixel;

class AliMUONClusterFinderAZ : public AliMUONClusterFinderVS 
{
public:
  AliMUONClusterFinderAZ(Bool_t draw = 0, Int_t iReco = 1);// Constructor
  virtual ~AliMUONClusterFinderAZ(); // Destructor

  void     FindRawClusters(); // the same interface as for old cluster finder
  void     EventLoop(Int_t nev, Int_t ch); // first event 
  Bool_t   TestTrack(Int_t t) const; // test if track was selected
 
protected:
  AliMUONClusterFinderAZ(const AliMUONClusterFinderAZ& rhs);
  AliMUONClusterFinderAZ& operator=(const AliMUONClusterFinderAZ& rhs);

 private:
  // Some constants
  static const Int_t fgkDim = 10000; // array size
  static const Double_t fgkCouplMin; // threshold on coupling 

  static  AliMUONClusterFinderAZ* fgClusterFinder; // the ClusterFinderAZ instance

  Int_t      fnPads[2];        // ! number of pads in the cluster on 2 cathodes
  Float_t    fXyq[7][fgkDim];    // ! pad information
  Int_t      fPadIJ[2][fgkDim];  // ! pad information
  //AZ AliSegmentation *fSegmentation[2]; // ! old segmentation
  AliMUONGeometrySegmentation *fSegmentation[2]; // ! new segmentation
  AliMUONResponse *fResponse;// ! response
  Float_t    fZpad;            // ! z-coordinate of the hit
  Int_t      fNpar;            // ! number of fit parameters
  Double_t   fQtot;            // ! total cluster charge
  Int_t      fReco;            // ! =1 if run reco with writing of reconstructed clusters 

  static     TMinuit* fgMinuit; // ! Fitter
  Bool_t     fUsed[2][fgkDim]; // ! flags for used pads
  TH2F*      fHist[4]; // ! histograms
  TClonesArray *fMuonDigits; // ! pointer to digits
  Bool_t     fDraw; // ! draw flag
  Int_t      fnMu; // ! number of muons passing thru the selected area
  Double_t   fxyMu[2][7]; // ! muon information
  TObjArray* fPixArray; // ! collection of pixels
  Int_t fnCoupled; // ! number of coupled clusters in precluster
  Int_t fDebug; // ! debug level

  // Functions

  void   ModifyHistos(void); // modify histograms
  void   AddPad(Int_t cath, Int_t digit); // add a pad to the cluster
  //AZ Bool_t Overlap(Int_t cath, TObject *dig); // check if the pad from one cathode overlaps with a pad in the cluster on the other cathode
  Bool_t Overlap(Int_t cath, AliMUONDigit *dig); // check if the pad from one cathode overlaps with a pad in the cluster on the other cathode
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
  Int_t  Fit(Int_t nfit, Int_t *clustFit, TObjArray **clusters, Double_t *parOk); // do the fitting 
  void  UpdatePads(Int_t nfit, Double_t *par); // subtract fitted charges from pads
  void  AddRawCluster(Double_t x, Double_t y, Double_t qTot, Double_t fmin, Int_t nfit, Int_t *tracks, Double_t sigx, Double_t sigy, Double_t dist); // add new reconstructed cluster
  Int_t FindLocalMaxima(Int_t *localMax, Double_t *maxVal); // find local maxima 
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
  void DrawCluster(Int_t nev0, Int_t ch0); // draw precluster
  Int_t Next(Int_t &nev0, Int_t &ch0); // commands for drawing

  // Dummy methods for overloading warnings
  void FindCluster(int, int, int, AliMUONRawCluster&) {return;}
  void FindLocalMaxima(AliMUONRawCluster*) {return;}
  void Split(AliMUONRawCluster*) {return;}
  void AddRawCluster(AliMUONRawCluster&) {return;}

ClassDef(AliMUONClusterFinderAZ,0) // cluster finder in MUON arm of ALICE
};

#endif
