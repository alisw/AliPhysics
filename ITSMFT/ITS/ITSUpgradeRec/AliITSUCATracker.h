//-*- Mode: C++ -*-
// ************************************************************************
// This file is property of and copyright by the ALICE ITSU Project       *
// ALICE Experiment at CERN, All rights reserved.                         *
// See cxx source for full Copyright notice                               *
//                                                                        *
//*************************************************************************

#ifndef ALIITSUCATRACKER_H
#define ALIITSUCATRACKER_H

//#define _TUNING_

#include <vector>

#include <TClonesArray.h>
#include "AliITSUCACell.h"
#include "AliITSUCATrackingStation.h"
typedef struct AliITSUCATrackingStation::ClsInfo ClsInfo_t;

#include "AliITSUTrackerGlo.h"

class AliITSUReconstructor;
class AliITSURecoDet;
class AliITSURecoLayer;
class AliITSUTrackCooked;

class TTree;
class AliCluster;
class AliESDEvent;

#ifdef _TUNING_
#include <TH1F.h>
#endif

using std::vector;

//__________________________________________________________________________________________________
class Doublets {
public:
  Doublets(int xx = 0, int yy = 0, float tL = 0.f, float ph = 0.f)
  : x((int)xx)
  , y((int)yy)
  , tanL(tL)
  , phi(ph) {}
  int x,y;
  float tanL, phi;
};

//__________________________________________________________________________________________________
class AliITSUCATracker : public AliITSUTrackerGlo {
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
  Double_t GetMaterialBudget(const double* p0, const double* p1, double& x2x0, double& rhol) const;
  Bool_t   GetSAonly() const { return fSAonly; }
  void     SetChi2Cut(float cut) { fChi2Cut = cut; }
  void     SetPhiCut(float cut) { fPhiCut = cut; }
  void     SetSAonly(Bool_t sa = kTRUE) { fSAonly = sa; }
  void     SetZCut(float cut) { fZCut = cut; }

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
  void ResetHistos();
#endif
  //
protected:
  bool   CellParams(int l, ClsInfo_t*__restrict__ c1, ClsInfo_t*__restrict__ c2,
                    ClsInfo_t*__restrict__ c3, float &curv, float np[3]);
  void   CellsTreeTraversal(vector<AliITSUCARoad> &roads, const int &iD, const int &doubl);
  void   FindTracksCA(int iteration);
  void   MakeCells(int iteration);
  Bool_t RefitAt(Double_t xx, AliITSUTrackCooked*__restrict__ t,
                 const AliITSUTrackCooked*__restrict__ c);
  void   SetCuts(int it);
  void   SetLabel(AliITSUTrackCooked &t, Float_t wrong);
  
private:
  AliITSUCATracker(const AliITSUCATracker&);
  AliITSUCATracker &operator=(const AliITSUCATracker &tr);

  // Data members

  // classes for interfacing the geometry, materials etc
  // Internal tracker arrays, layers, modules, etc
  AliITSUCATrackingStation        fLayer[7];
  vector<bool>                    fUsedClusters[7];
  Float_t                         fChi2Cut;
  Float_t                         fPhiCut;
  Float_t                         fZCut;
  vector<Doublets>                fDoublets[6];
  vector<AliITSUCACell>           fCells[5];
  TClonesArray                   *fCandidates[4];
  Bool_t                          fSAonly;             // kTRUE if the standalone tracking only

  
  // Cuts
  float fCPhi;
  float fCDTanL;
  float fCDPhi;
  float fCZ;
  float fCDCAz[5];
  float fCDCAxy[5];
  float fCDN[4];
  float fCDP[4];
  float fCDZ[6];

  static const Double_t           fgkChi2Cut;      // chi2 cut during track merging
  static const int                fgkNumberOfIterations;
  static const float              fgkR[7];
  //
  ClassDef(AliITSUCATracker,2)   //ITSU stand-alone tracker
};

#endif // ALIITSUCATRACKER_H
