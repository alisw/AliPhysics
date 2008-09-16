#ifndef ALIITSMULTRECONSTRUCTOR_H
#define ALIITSMULTRECONSTRUCTOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/////////////////////////////////////////////////////////////////////////
// 
// AliITSMultReconstructor - find clusters in the pixels (theta and
// phi) and tracklets.
// 
// These can be used to extract charged particles multiplcicity from the ITS.
//
// A tracklet consist of two ITS clusters, one in the first pixel
// layer and one in the second. The clusters are associates if the 
// differencies in Phi (azimuth) and Zeta (longitudinal) are inside 
// a fiducial volume. In case of multiple candidates it is selected the
// candidate with minimum distance in Phi. 
// The boolean fOnlyOneTrackletPerC2 allows to control if two clusters 
// in layer 2 can be associated to the same cluster in layer 1 or not.
//
/////////////////////////////////////////////////////////////////////////

#include "TObject.h"

class TTree;
class TH1F;
class TH2F; 

class AliITSgeom;

class AliITSMultReconstructor : public TObject 
{
public:
  AliITSMultReconstructor();
  virtual ~AliITSMultReconstructor();

  void Reconstruct(TTree* tree, Float_t* vtx, Float_t* vtxRes);
  void LoadClusterFiredChips(TTree* tree);
  void FlagClustersInOverlapRegions(Int_t ic1,Int_t ic2);

  // Following members are set via AliITSRecoParam
  void SetOnlyOneTrackletPerC2(Bool_t b = kTRUE) {fOnlyOneTrackletPerC2 = b;}
  void SetPhiWindow(Float_t w=0.08) {fPhiWindow=w;}
  void SetZetaWindow(Float_t w=1.) {fZetaWindow=w;}
  void SetRemoveClustersFromOverlaps(Bool_t b = kFALSE) {fRemoveClustersFromOverlaps = b;}
  void SetPhiOverlapCut(Float_t w=0.005) {fPhiOverlapCut=w;}
  void SetZetaOverlapCut(Float_t w=0.05) {fZetaOverlapCut=w;}

  Int_t GetNClustersLayer1() const {return fNClustersLay1;}
  Int_t GetNClustersLayer2() const {return fNClustersLay2;}
  Int_t GetNTracklets() const {return fNTracklets;}
  Int_t GetNSingleClusters() const {return fNSingleCluster;}
  Short_t GetNFiredChips(Int_t layer) const {return fNFiredChips[layer];}

  Float_t* GetClusterLayer1(Int_t n) {return fClustersLay1[n];}
  Float_t* GetClusterLayer2(Int_t n) {return fClustersLay2[n];}
  Float_t* GetTracklet(Int_t n) {return fTracklets[n];}
  Float_t* GetCluster(Int_t n) {return fSClusters[n];}

  void SetHistOn(Bool_t b=kFALSE) {fHistOn=b;}
  void SaveHists();

protected:
  AliITSMultReconstructor(const AliITSMultReconstructor& mr);
  AliITSMultReconstructor& operator=(const AliITSMultReconstructor& mr);

  
  Float_t**     fClustersLay1;               // clusters in the 1st layer of ITS 
  Float_t**     fClustersLay2;               // clusters in the 2nd layer of ITS 
  Int_t*        fDetectorIndexClustersLay1;  // module index for clusters 1st ITS layer
  Int_t*        fDetectorIndexClustersLay2;  // module index for clusters 2nd ITS layer
  Bool_t*       fOverlapFlagClustersLay1;    // flag for clusters in the overlap regions 1st ITS layer
  Bool_t*       fOverlapFlagClustersLay2;    // flag for clusters in the overlap regions 2nd ITS layer 


  Float_t**     fTracklets;            // tracklets 
  Float_t**     fSClusters;            // single clusters (unassociated)
  Bool_t*       fAssociationFlag;      // flag for the associations 
  
  Int_t         fNClustersLay1;        // Number of clusters (Layer1)
  Int_t         fNClustersLay2;        // Number of clusters (Layer2)
  Int_t         fNTracklets;           // Number of tracklets
  Int_t         fNSingleCluster;       // Number of unassociated clusters
  Short_t       fNFiredChips[2];       // Number of fired chips in the two SPD layers
 
  // Following members are set via AliITSRecoParam
  Bool_t        fOnlyOneTrackletPerC2;         // Allow only one tracklet per cluster in the outer layer
  Float_t       fPhiWindow;                    // Search window in phi
  Float_t       fZetaWindow;                   // Search window in eta
  Bool_t        fRemoveClustersFromOverlaps;   // Option to skip clusters in the overlaps
  Float_t       fPhiOverlapCut;                // Fiducial window in phi for overlap cut
  Float_t       fZetaOverlapCut;               // Fiducial window in eta for overlap cut

  Bool_t        fHistOn;               // Option to define and fill the histograms 

  TH1F*         fhClustersDPhiAcc;     // Phi2 - Phi1 for tracklets 
  TH1F*         fhClustersDThetaAcc;   // Theta2 - Theta1 for tracklets 
  TH1F*         fhClustersDZetaAcc;    // z2 - z1projected for tracklets 
  TH1F*         fhClustersDPhiAll;     // Phi2 - Phi1 all the combinations 
  TH1F*         fhClustersDThetaAll;   // Theta2 - Theta1 all the combinations
  TH1F*         fhClustersDZetaAll;    // z2 - z1projected all the combinations
 
  TH2F*         fhDPhiVsDThetaAll;     // 2D plot for all the combinations  
  TH2F*         fhDPhiVsDThetaAcc;     // same plot for tracklets 
  TH2F*         fhDPhiVsDZetaAll;      // 2d plot for all the combination 
  TH2F*         fhDPhiVsDZetaAcc;      // same plot for tracklets 

  TH1F*         fhetaTracklets;        // Pseudorapidity distr. for tracklets 
  TH1F*         fhphiTracklets;        // Azimuthal (Phi) distr. for tracklets  
  TH1F*         fhetaClustersLay1;     // Pseudorapidity distr. for Clusters L. 1
  TH1F*         fhphiClustersLay1;     // Azimuthal (Phi) distr. for Clusters L. 1 


  void LoadClusterArrays(TTree* tree);

  ClassDef(AliITSMultReconstructor,3)
};

#endif
