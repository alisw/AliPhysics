#ifndef ALIITSUVERTEXER_H
#define ALIITSUVERTEXER_H
//#define MC_CHECK // comment out to enable MC checks for debugging

#include "AliVertexer.h"

/* Copyright(c) 2009-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
/////////////////////////////////////////////////////////////////////////////
//     Class for the reconstruction of the primary vertices using ITSU     //
/////////////////////////////////////////////////////////////////////////////

class TClonesArray;
class AliESDVertex;

class AliITSUVertexer : public AliVertexer {
 public:
  
  // Constructors and destructors
  AliITSUVertexer(Double_t phicut=0.005,Double_t zcut=0.002,Double_t paircut=0.04, Double_t clustercut=0.8, Int_t clcontrib=5);
  virtual ~AliITSUVertexer();

  // Public methods
  virtual AliESDVertex* GetAllVertices(Int_t& nVert) const { nVert=fNoVertices; return fVertices; };
  virtual AliESDVertex* FindVertexForCurrentEvent(TTree *); 
  virtual void PrintStatus() const;

  // Getters
  UInt_t   GetNoLines()       const { return fNoLines;    }
  UShort_t GetNumOfVertices() const { return fNoVertices; }

  // Setters
  void SetPhiCut(Double_t phicut) { fPhiCut=phicut; }
  void SetZCut(Double_t zcut)     { fZCut=zcut; }

  #ifdef MC_CHECK
  // Debug + MC truth
  UInt_t* GetParticleId(UInt_t &num) const { num=fGoodLines; return fParticleId; }
  UInt_t GetGoodLines() const { return fGoodLines; }
  UInt_t GetGoodLinesPhi() const { return fGoodLinesPhi; }
  UInt_t GetLinesPhi() const { return fLinesPhi; }
  #endif

 protected:
  // Methods
  AliITSUVertexer(AliITSUVertexer&);
  AliITSUVertexer& operator=(const AliITSUVertexer& other);
  void AddToCluster(UInt_t line,Bool_t weight=kFALSE,Int_t cl=-1);
  void CleanAndOrderClusters();
  void Clusterize(UInt_t l1, UInt_t l2, Bool_t weight=kFALSE);
  void ComputeClusterCentroid(UInt_t cl);
  void FindTracklets();
  void FindVerticesForCurrentEvent();
  Int_t MatchPoints(UShort_t layer, Double_t anchor, Double_t *p0=0x0, Double_t *p1=0x0);
  void MoveLabels(Short_t start, Short_t end);
  void Reset();
  void SortClusters();

  // Data members
  Int_t fClusterContribCut;
  Double_t fClusterCut;
  Int_t *fClusterIndex[3];           // AliITSUClusterPix index 
  Double_t *fClusterPhi[3];          // Phi of clusters
  TClonesArray *fClusters[3];        //! array of pointers to TClonesArray of AliITSUClusterPix not owned by this class
  TClonesArray fLines;               //! array of tracklets
  TClonesArray fLinesClusters;       // array of vertex candidates
  UInt_t fLinesPhi;                  // number of tracklets built by using the first two layers
  UInt_t fNoClusters;                // number of clusters
  UInt_t fNoLines;                   // number of tracklets
  UShort_t fNoVertices;              // number of vertices
  Double_t fPairCut;                 // cut on pair
  Double_t fPhiCut;                  // cut on deltaphi for cluster matching among first two layers
  Double_t fZCut;                    // cut on deltatheta for cluster matching among first two layers and the third one
  Bool_t *fUsedClusters[3];          // flag for used clusters in tracklet formation
  Short_t *fUsedLines;               // flag for used lines
  AliESDVertex *fVertices;           // array of vertices

  #ifdef MC_CHECK
  // MC truth methods
  Bool_t CheckMC(UInt_t,UInt_t,UInt_t);

  // MC truth data members
  UInt_t fGoodLines;
  UInt_t fGoodLinesPhi;
  UInt_t *fParticleId;  
  #endif

  ClassDef(AliITSUVertexer,1)
};

#endif
