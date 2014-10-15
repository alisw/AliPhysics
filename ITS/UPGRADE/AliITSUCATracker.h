//-*- Mode: C++ -*-
// ************************************************************************
// This file is property of and copyright by the ALICE ITSU Project       *
// ALICE Experiment at CERN, All rights reserved.                         *
// See cxx source for full Copyright notice                               *
//                                                                        *
//*************************************************************************

#ifndef ALIITSUCATRACKER_H
#define ALIITSUCATRACKER_H

#define _TUNING_

#include "AliTracker.h"
#include "AliITSUGeomTGeo.h"
#include <TClonesArray.h>
#include <vector>
#include "AliITSUMatLUT.h"
#include "AliITSUAux.h"
#include "AliExternalTrackParam.h"
#include "AliITSUCATrackingStation.h"
#include "AliITSUCACell.h"
#include "AliITSUTrackCooked.h"

class AliITSUReconstructor;
class AliITSURecoDet;
class AliITSURecoLayer;
class AliITSUTrackCooked;

class TTree;
class AliCluster;
class AliESDEvent;

typedef struct AliITSUCATrackingStation::ClsInfo ClsInfo_t;

using std::vector;

//__________________________________________________________________________________________________
class Doublets {
public:
  Doublets(int xx = 0, int yy = 0, float tL = 0.f, float ph = 0.f)
  : x((unsigned short)xx)
  , y((unsigned short)yy)
  , tanL(tL)
  , phi(ph) {}
  unsigned short x,y;
  float tanL, phi;
};

//__________________________________________________________________________________________________
class AliITSUCATracker : public AliTracker {
public:
  AliITSUCATracker(AliITSUReconstructor* rec=0);
  virtual ~AliITSUCATracker();

  // These functions must be implemented
  Int_t Clusters2Tracks(AliESDEvent *event);
  Int_t PropagateBack(AliESDEvent *event);
  Int_t RefitInward(AliESDEvent *event);
  Int_t LoadClusters(TTree *ct);
  void UnloadClusters();
  AliCluster *GetCluster(Int_t index) const;

  // Possibly, other public functions
  void     Init(AliITSUReconstructor* rec);
  Double_t RefitTrack(AliExternalTrackParam* trc, Int_t clInfo[2 * AliITSUAux::kMaxLayers], Double_t rDest, Int_t stopCond);
  Bool_t   PropagateSeed(AliExternalTrackParam *seed, Double_t xToGo, Double_t mass, Double_t maxStep=1.0, Bool_t matCorr=kTRUE);
  Double_t GetMaterialBudget(const double* pnt0, const double* pnt1, double& x2x0, double& rhol) const;
  Bool_t   GoToEntranceToLayer(AliExternalTrackParam* seed, AliITSURecoLayer* lr, Int_t dir, Bool_t check=kFALSE);
  Bool_t   GoToExitFromLayer(AliExternalTrackParam* seed, AliITSURecoLayer* lr, Int_t dir, Bool_t check=kTRUE);
  Bool_t   TransportToLayerX(AliExternalTrackParam* seed, Int_t lFrom, Int_t lTo, Double_t xStop);
  Bool_t   TransportToLayer(AliExternalTrackParam* seed, Int_t lFrom, Int_t lTo, Double_t rLim=-1);

  void SetChi2Cut(float cut) { fChi2Cut = cut; }
  void SetPhiCut(float cut) { fPhiCut = cut; }
  void SetZCut(float cut) { fZCut = cut; }

#ifdef _TUNING_
  bool                            fGood;
  TH1F *                          fGoodCombChi2[5];
  TH1F *                          fFakeCombChi2[5];
  TH1F *                          fGoodCombN[4];
  TH1F *                          fFakeCombN[4];
  TH1F *                          fGDZ[6];
  TH1F *                          fGDXY[6];
  TH1F *                          fFDZ[6];
  TH1F *                          fFDXY[6];
  TH1F *                          fGDCAZ[5];
  TH1F *                          fGDCAXY[5];
  TH1F *                          fFDCAZ[5];
  TH1F *                          fFDCAXY[5];
  TH1F *                          fTan;
  TH1F *                          fTanF;
  TH1F *                          fPhi;
  TH1F *                          fPhiF;
  TH1F *                          fNEntries;
#endif
  //
protected:
  AliITSUCATracker(const AliITSUCATracker&);

  //  void MakeTriplets();
  void CellsTreeTraversal(vector<AliITSUCARoad> &roads, const int &iD, const int &doubl);
  Bool_t InitTrackParams(AliITSUTrackCooked &track, int points[]);
  void FindTracksCA(int iteration);
  void GlobalFit();
  void MergeTracks( vector<AliITSUTrackCooked> &vec, bool flags[] );
  float FilterSeed(AliITSUCACell &c1, AliITSUCACell &c2, int lrStart);
  bool CellParams(int l, ClsInfo_t* c1, ClsInfo_t* c2, ClsInfo_t* c3, float &curv, float np[3]);

  Bool_t RefitAt(Double_t xx, AliITSUTrackCooked *t, const AliITSUTrackCooked *c);
private:
  AliITSUCATracker &operator=(const AliITSUCATracker &tr);
  void SetLabel(AliITSUTrackCooked &t, Float_t wrong);

  // Data members

  // classes for interfacing the geometry, materials etc.
  AliITSUReconstructor*           fReconstructor;  // ITS global reconstructor
  AliITSURecoDet*                 fITS;            // interface to ITS, borrowed from reconstructor
  AliITSUMatLUT*                  fMatLUT;         // material lookup table
  Bool_t                          fUseMatLUT;      //! use material lookup table rather than TGeo
  Double_t                        fCurrMass;       // assumption about particle mass


  // Internal tracker arrays, layers, modules, etc
  AliITSUCATrackingStation        fLayer[7];
  vector<bool>                    fUsedClusters[7];
  Float_t                         fChi2Cut;
  Float_t                         fPhiCut;
  Float_t                         fZCut;
  vector<Doublets>                fDoublets[6];
  vector<AliITSUCACell>           fCells[5];
  TClonesArray                   *fCandidates[4];

  static const Double_t           fgkToler;        // tracking tolerance
  static const Double_t           fgkChi2Cut;      // chi2 cut during track merging
  static const int                fgkNumberOfIterations;
  static const float              fgkR[7];
  //
  ClassDef(AliITSUCATracker,2)   //ITSU stand-alone tracker
};

#endif // ALIITSUCATRACKER_H
