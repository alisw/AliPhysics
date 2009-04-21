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
#include <AliESDVertex.h>
#include <TH3F.h>
#include <TBits.h>

class AliITSVertexer3D : public AliITSVertexer {

 public:

  AliITSVertexer3D();
  virtual ~AliITSVertexer3D();
  virtual AliESDVertex* FindVertexForCurrentEvent(TTree *itsClusterTree);
  void FindVertex3DIterative();
  void FindVertex3D(TTree *itsClusterTree);
  AliESDVertex GetVertex3D() const {return fVert3D;}
  virtual void PrintStatus() const;
  static Bool_t DistBetweenVertices(AliESDVertex &a, AliESDVertex &b, Double_t test, Double_t &dist);
  void SetWideFiducialRegion(Float_t dz = 20.0, Float_t dr=2.5){
    SetCoarseMaxRCut(dr);
    SetZCutDiamond(dz);
  }
  void SetNarrowFiducialRegion(Float_t dz = 0.5, Float_t dr=0.5){
    SetMaxRCut(dr);
    SetMaxZCut(dz);
  }
  void SetDeltaPhiCuts(Float_t dphiloose=0.5, Float_t dphitight=0.01){
    SetCoarseDiffPhiCut(dphiloose);
    SetDiffPhiMax(dphitight);
  }
  void SetCoarseDiffPhiCut(Float_t dphi = 0.5){fCoarseDiffPhiCut=dphi;}
  void SetFineDiffPhiCut(Float_t dphi = 0.05){fFineDiffPhiCut=dphi;}
  void SetCutOnPairs(Float_t cp = 0.1){fCutOnPairs = cp;}
  void SetCoarseMaxRCut(Float_t rad = 2.5){fCoarseMaxRCut=rad;}
  void SetMaxRCut(Float_t rad = 0.5){fMaxRCut=rad;}
  void SetZCutDiamond(Float_t zcut = 20.0){fZCutDiamond=zcut;}
  void SetMaxZCut(Float_t dz = 0.5){fMaxZCut=dz;}
  void SetDCACut(Float_t dca=0.1){fDCAcut=dca;} 
  void SetDiffPhiMax(Float_t pm = 0.01){fDiffPhiMax = pm;}
  void SetMeanPSelTracks(Float_t pGeV=0.875){fMeanPSelTrk = pGeV;}
  void SetMeanPtSelTracks(Float_t ptGeV=0.630){fMeanPtSelTrk = ptGeV;}
  void SetMeanPPtSelTracks(Float_t fieldTesla);
  void SetMinDCAforPileup(Float_t mindist=0.1) {fDCAforPileup=mindist;}
  void SetPileupAlgo(UShort_t optalgo=0){fPileupAlgo=optalgo;}

protected:
  AliITSVertexer3D(const AliITSVertexer3D& vtxr);
  AliITSVertexer3D& operator=(const AliITSVertexer3D& /* vtxr */);
  Int_t FindTracklets(TTree *itsClusterTree, Int_t optCuts);
  Int_t Prepare3DVertex(Int_t optCuts);
  void ResetVert3D();
  void FindPeaks(TH3F* histo, Double_t *peak, Int_t &nOfTracklets, Int_t &nOfTimes);
  void PileupFromZ();
  void MarkUsedClusters();
  Int_t RemoveTracklets();
  void  FindOther3DVertices(TTree *itsClusterTree);

  enum {kMaxCluPerMod=250};

  TClonesArray fLines;      //! array of tracklets
  AliESDVertex fVert3D;        // 3D Vertex
  Float_t fCoarseDiffPhiCut; // loose cut on DeltaPhi for RecPoint matching
  Float_t fFineDiffPhiCut; // tight value of DeltaPhi for RP matching (2nd method) 
  Float_t fCutOnPairs; //cut on distance between pairs of tracklets 
  Float_t fCoarseMaxRCut; // cut on tracklet DCA to Z axis
  Float_t fMaxRCut; // cut on tracklet DCA to beam axis
  Float_t fZCutDiamond;   // cut on +-Z of the diamond
  Float_t fMaxZCut;   // cut on Z distance from estimated vertex
  Float_t fDCAcut; // cut on tracklet to tracklet and tracklet to vertex DCA
  Float_t fDiffPhiMax;     // Maximum delta phi allowed among corr. pixels
  Float_t fMeanPSelTrk; // GeV, mean P for tracks with dphi<0.01 rad
  Float_t fMeanPtSelTrk; // GeV, mean Pt for tracks with dphi<0.01 rad
  TBits   fUsedCluster;  // flag for used clusters in vertex calculation
  TH1F *fZHisto;           //! histogram with coarse z distribution
  Float_t  fDCAforPileup;  // Minimum DCA to 1st vertex for pileup tracklets 
  UShort_t fPileupAlgo;    // Algo for pileup identification
                           // 0->VertexerZ pileup algo
                           // 1->Unused RecPoints algo

  ClassDef(AliITSVertexer3D,8);

};

#endif
