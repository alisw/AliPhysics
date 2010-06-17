#ifndef ALIITSMULTRECONSTRUCTOR_H
#define ALIITSMULTRECONSTRUCTOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//_________________________________________________________________________
// 
//        Implementation of the ITS-SPD trackleter class
//
// It retrieves clusters in the pixels (theta and phi) and finds tracklets.
// These can be used to extract charged particle multiplicity from the ITS.
//
// A tracklet consists of two ITS clusters, one in the first pixel layer and 
// one in the second. The clusters are associated if the differences in 
// Phi (azimuth) and Theta (polar angle) are within fiducial windows.
// In case of multiple candidates the candidate with minimum
// distance is selected. 
//_________________________________________________________________________
#include "AliTrackleter.h"

class TBits;
class TTree;
class TH1F;
class TH2F; 
class AliITSDetTypeRec;
class AliITSgeom;
class AliESDEvent;
class AliESDtrack;
class AliVertex;
class AliMultiplicity;

class AliITSMultReconstructor : public AliTrackleter
{
public:
  //
  enum {kClTh,kClPh,kClZ,kClMC0,kClMC1,kClMC2,kClNPar};
  enum {kTrTheta,kTrPhi,kTrDPhi,kTrDTheta,kTrLab1,kTrLab2,kClID1,kClID2,kTrNPar};
  enum {kSCTh,kSCPh,kSCID,kSCNPar};
  enum {kITSTPC,kITSSAP,kITSTPCBit=BIT(kITSTPC),kITSSAPBit=BIT(kITSSAP)}; // RS
  AliITSMultReconstructor();
  virtual ~AliITSMultReconstructor();

  void Reconstruct(AliESDEvent* esd, TTree* treeRP);
  void Reconstruct(TTree* tree, Float_t* vtx, Float_t* vtxRes);   // old reconstructor invocation
  void FindTracklets(const Float_t* vtx); 
  void LoadClusterFiredChips(TTree* tree);
  void FlagClustersInOverlapRegions(Int_t ic1,Int_t ic2);
  void FlagTrackClusters(const AliESDtrack* track);
  void FlagIfPrimary(AliESDtrack* track, const AliVertex* vtx);
  void ProcessESDTracks();
  
  void CreateMultiplicityObject();
  //
  // Following members are set via AliITSRecoParam
  void SetPhiWindow(Float_t w=0.08) {fPhiWindow=w;}
  void SetThetaWindow(Float_t w=0.025) {fThetaWindow=w;}
  void SetPhiShift(Float_t w=0.0045) {fPhiShift=w;}
  void SetRemoveClustersFromOverlaps(Bool_t b = kFALSE) {fRemoveClustersFromOverlaps = b;}
  void SetPhiOverlapCut(Float_t w=0.005) {fPhiOverlapCut=w;}
  void SetZetaOverlapCut(Float_t w=0.05) {fZetaOverlapCut=w;}

  Int_t GetNClustersLayer1() const {return fNClustersLay1;}
  Int_t GetNClustersLayer2() const {return fNClustersLay2;}
  Int_t GetNTracklets() const {return fNTracklets;}
  Int_t GetNSingleClusters() const {return fNSingleCluster;}
  Short_t GetNFiredChips(Int_t layer) const {return fNFiredChips[layer];}

  Float_t* GetClusterLayer1(Int_t n) {return &fClustersLay1[n*kClNPar];}
  Float_t* GetClusterLayer2(Int_t n) {return &fClustersLay2[n*kClNPar];}

  Float_t* GetTracklet(Int_t n) {return fTracklets[n];}
  Float_t* GetCluster(Int_t n) {return fSClusters[n];}

  void SetHistOn(Bool_t b=kFALSE) {fHistOn=b;}
  void SaveHists();

  AliITSDetTypeRec *GetDetTypeRec() const {return fDetTypeRec;}
  void SetDetTypeRec(AliITSDetTypeRec *ptr){fDetTypeRec = ptr;}

protected:
  AliITSMultReconstructor(const AliITSMultReconstructor& mr);
  AliITSMultReconstructor& operator=(const AliITSMultReconstructor& mr);
  AliITSDetTypeRec* fDetTypeRec;            //! pointer to DetTypeRec
  AliESDEvent*      fESDEvent;              //! pointer to ESD event
  TTree*            fTreeRP;                //! ITS recpoints

  Char_t*       fUsedClusLay1;               // RS: flag of clusters usage in ESD tracks: 0=unused, bit0=TPC/ITS+ITSSA, bit1=ITSSA_Pure
  Char_t*       fUsedClusLay2;               // RS: flag of clusters usage in ESD tracks: 0=unused, bit0=TPC/ITS+ITSSA, bit1=ITSSA_Pure

  Float_t*      fClustersLay1;               // clusters in the 1st layer of ITS 
  Float_t*      fClustersLay2;               // clusters in the 2nd layer of ITS 
  Int_t*        fDetectorIndexClustersLay1;  // module index for clusters 1st ITS layer
  Int_t*        fDetectorIndexClustersLay2;  // module index for clusters 2nd ITS layer
  Bool_t*       fOverlapFlagClustersLay1;    // flag for clusters in the overlap regions 1st ITS layer
  Bool_t*       fOverlapFlagClustersLay2;    // flag for clusters in the overlap regions 2nd ITS layer 

  Float_t**     fTracklets;            // tracklets 
  Float_t**     fSClusters;            // single clusters (unassociated)
  
  Int_t         fNClustersLay1;        // Number of clusters (Layer1)
  Int_t         fNClustersLay2;        // Number of clusters (Layer2)
  Int_t         fNTracklets;           // Number of tracklets
  Int_t         fNSingleCluster;       // Number of unassociated clusters
  Short_t       fNFiredChips[2];       // Number of fired chips in the two SPD layers
 
  // Following members are set via AliITSRecoParam
  Float_t       fPhiWindow;                    // Search window in phi
  Float_t       fThetaWindow;                  // Search window in theta
  Float_t       fPhiShift;                     // Phi shift reference value (at 0.5 T) 
  Bool_t        fRemoveClustersFromOverlaps;   // Option to skip clusters in the overlaps
  Float_t       fPhiOverlapCut;                // Fiducial window in phi for overlap cut
  Float_t       fZetaOverlapCut;               // Fiducial window in eta for overlap cut

  Bool_t        fHistOn;               // Option to define and fill the histograms 

  TH1F*         fhClustersDPhiAcc;     // Phi2 - Phi1 for tracklets 
  TH1F*         fhClustersDThetaAcc;   // Theta2 - Theta1 for tracklets 
  TH1F*         fhClustersDPhiAll;     // Phi2 - Phi1 all the combinations 
  TH1F*         fhClustersDThetaAll;   // Theta2 - Theta1 all the combinations
 
  TH2F*         fhDPhiVsDThetaAll;     // 2D plot for all the combinations  
  TH2F*         fhDPhiVsDThetaAcc;     // same plot for tracklets 

  TH1F*         fhetaTracklets;        // Pseudorapidity distr. for tracklets 
  TH1F*         fhphiTracklets;        // Azimuthal (Phi) distr. for tracklets  
  TH1F*         fhetaClustersLay1;     // Pseudorapidity distr. for Clusters L. 1
  TH1F*         fhphiClustersLay1;     // Azimuthal (Phi) distr. for Clusters L. 1 


  void LoadClusterArrays(TTree* tree);

  ClassDef(AliITSMultReconstructor,7)
};

#endif
