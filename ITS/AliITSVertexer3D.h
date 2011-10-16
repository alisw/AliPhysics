#ifndef ALIITSVERTEXER3D_H
#define ALIITSVERTEXER3D_H

#include<AliITSVertexer.h>

///////////////////////////////////////////////////////////////////
//                                                               //
// Class for primary vertex finding  (3D reconstruction)         //
//                                                               //
///////////////////////////////////////////////////////////////////

/* $Id$ */

#include <TClonesArray.h>
#include <TRandom3.h>
#include <AliESDVertex.h>
#include <TH3F.h>
#include <TBits.h>

class AliITSVertexer3D : public AliITSVertexer {

 public:

  AliITSVertexer3D();
  virtual ~AliITSVertexer3D();
  virtual AliESDVertex* FindVertexForCurrentEvent(TTree *itsClusterTree);
  void FindVertex3DIterative();
  void FindVertex3DIterativeMM();
  void FindVertex3D(TTree *itsClusterTree);
  AliESDVertex GetVertex3D() const {return fVert3D;}
  //  Double_t *Get3DPeak() {return f3DPeak;}
  virtual void PrintStatus() const;
  static Bool_t DistBetweenVertices(AliESDVertex &a, AliESDVertex &b, Double_t test, Double_t &dist);
  void SetWideFiducialRegion(Double_t dz = 40.0, Double_t dr=2.5){
    SetCoarseMaxRCut(dr);
    SetZCutDiamond(dz);
  }
  void SetNarrowFiducialRegion(Double_t dz = 0.5, Double_t dr=0.5){
    SetMaxRCut(dr);
    SetMaxZCut(dz);
  }
  void SetDeltaPhiCuts(Double_t dphiloose=0.5, Double_t dphitight=0.025){
    SetCoarseDiffPhiCut(dphiloose);
    SetDiffPhiMax(dphitight);
  }
  void SetCoarseDiffPhiCut(Double_t dphi = 0.5){fCoarseDiffPhiCut=dphi;}
  void SetFineDiffPhiCut(Double_t dphi = 0.05){fFineDiffPhiCut=dphi;}
  void SetCutOnPairs(Double_t cp = 0.15){fCutOnPairs = cp;}
  void SetCoarseMaxRCut(Double_t rad = 2.5){fCoarseMaxRCut=rad;}
  void SetMaxRCut(Double_t rad = 0.5){fMaxRCut=rad;}
  void SetMaxRCutAlgo2(Double_t rad = 0.2){fMaxRCut2=rad;}
  void SetZCutDiamond(Double_t zcut = 40.0){fZCutDiamond=zcut;}
  void SetMaxZCut(Double_t dz = 0.5){fMaxZCut=dz;}
  void SetDCACut(Double_t dca=0.1){fDCAcut=dca;} 
  void SetDiffPhiMax(Double_t pm = 0.025){fDiffPhiMax = pm;}
  void SetMeanPSelTracks(Double_t pGeV=0.875){fMeanPSelTrk = pGeV;}
  void SetMeanPtSelTracks(Double_t ptGeV=0.630){fMeanPtSelTrk = ptGeV;}
  void SetMeanPPtSelTracks(Double_t fieldTesla);
  void SetMinDCAforPileup(Double_t mindist=0.1) {fDCAforPileup=mindist;}
  void SetDeltaPhiforPileup(Double_t dphi=0.01) {fDiffPhiforPileup=dphi;}
  void SetPileupAlgo(UShort_t optalgo=1){fPileupAlgo=optalgo;}
  void SetBinSizeR(Double_t siz=0.1){fBinSizeR=siz;}
  void SetBinSizeZ(Double_t siz=0.8){fBinSizeZ=siz;}
  void SetHighMultAlgo(UChar_t n){
    if(n<2) fHighMultAlgo=n;
    else AliError("Only algos 0 and 1 implemented");
  }
  void SetHighMultDownscalingAlgo(){fHighMultAlgo=0;}
  void SetHighMultTracesAlgo(){fHighMultAlgo=1;}

  void SetMaxNumOfClustersForHighMult(Int_t ncl){fMaxNumOfCl=ncl;}
  void SetMaxNumOfClustersForDownScale(Int_t ncl){fMaxNumOfClForDownScale=ncl;}
  void SetMaxNumOfClustersForRebin(Int_t ncl){fMaxNumOfClForRebin=ncl;}
  Int_t GetMaxNumOfClustersForHighMult() const {return fMaxNumOfCl;}
  Int_t GetMaxNumOfClustersForDownScale() const {return fMaxNumOfClForDownScale;}
  Int_t GetMaxNumOfClustersForRebin() const {return fMaxNumOfClForRebin;}

protected:
  AliITSVertexer3D(const AliITSVertexer3D& vtxr);
  AliITSVertexer3D& operator=(const AliITSVertexer3D& /* vtxr */);
  Int_t FindTracklets(TTree *itsClusterTree, Int_t optCuts);
  Int_t Prepare3DVertex(Int_t optCuts);
  Int_t Prepare3DVertexPbPb();
  void ResetVert3D();
  void FindPeaks(TH3F* histo, Double_t *peak, Int_t &nOfTracklets, Int_t &nOfTimes);
  void PileupFromZ();
  void MarkUsedClusters();
  Int_t RemoveTracklets();
  void  FindOther3DVertices(TTree *itsClusterTree);
  Double_t GetFraction(Int_t itr) const;

  enum {kMaxCluPerMod=250};
  enum {kMaxPileupVertices=10};
  TClonesArray fLines;      //! array of tracklets
  AliESDVertex fVert3D;        // 3D Vertex
  Double_t fCoarseDiffPhiCut; // loose cut on DeltaPhi for RecPoint matching
  Double_t fFineDiffPhiCut; // tight value of DeltaPhi for RP matching (2nd method) 
  Double_t fCutOnPairs; //cut on distance between pairs of tracklets 
  Double_t fCoarseMaxRCut; // cut on tracklet DCA to Z axis
  Double_t fMaxRCut;     // cut on tracklet DCA to beam axis
  Double_t fMaxRCut2;    // cut on tracklet DCA to beam axis - algo2
  Double_t fZCutDiamond;   // cut on +-Z of the diamond
  Double_t fMaxZCut;   // cut on Z distance from estimated vertex
  Double_t fDCAcut; // cut on tracklet to tracklet and tracklet to vertex DCA
  Double_t fDiffPhiMax;     // Maximum delta phi allowed among corr. pixels
  Double_t fMeanPSelTrk; // GeV, mean P for tracks with dphi<0.01 rad
  Double_t fMeanPtSelTrk; // GeV, mean Pt for tracks with dphi<0.01 rad
  TBits   fUsedCluster;  // flag for used clusters in vertex calculation
  TH1F *fZHisto;           //! histogram with coarse z distribution
  Double_t  fDCAforPileup;  // Minimum DCA to 1st vertex for pileup tracklets 
  Double_t  fDiffPhiforPileup;  // Cut on delta phi for pileup 
  Double_t  fBinSizeR;      // Histo3D bin size along radius
  Double_t  fBinSizeZ;      // Histo3D bin size along z
  UShort_t fPileupAlgo;    // Algo for pileup identification
                           // 0->VertexerZ pileup algo
                           // 1->Unused RecPoints algo
  Int_t fMaxNumOfCl;          // max n. of clusters on L1 or L2 for high mult definition
  Int_t fMaxNumOfClForRebin;  // max n. of clusters on L1 or L2 for rebin
  Int_t fMaxNumOfClForDownScale;  // max n. of clusters on L1 or L2 for downscale
  Int_t  fNRecPLay1;       // number of rec ponts on SPD layer 1
  Int_t  fNRecPLay2;       // number of rec ponts on SPD layer 2
  Float_t f3DBinSize;           // Size of the 3D bins
  Bool_t fDoDownScale;     // Control downscaling of tracklets in high mult
  TRandom3 *fGenerForDownScale; // randomnumber generator fordownscaling
  Double_t f3DPeak[3];           // TH3F peak coords
  UChar_t fHighMultAlgo;    // algorithm used for high mult. events
  Bool_t fSwitchAlgorithm; // Switch between two algoritms in testing phase

  static const Int_t fgkMaxNumOfClDefault;      // Default max n. of clusters for downscale
  static const Int_t fgkMaxNumOfClRebinDefault; // Default max n. of clusters for rebin
  static const Int_t fgkMaxNumOfClDownscaleDefault; // Default max n. of clusters for rebin
  static const Float_t fgk3DBinSizeDefault;  // Default 3D bins size

  ClassDef(AliITSVertexer3D,15);

};

#endif
