#ifndef ALIITSUTRACKERSA_H
#define ALIITSUTRACKERSA_H

#define __DEBUG__
#ifdef __DEBUG__ 
#include <TString.h>
#endif

//-------------------------------------------------------------------------
//                   The stand-alone ITSU tracker
//     It reads AliITSUClusterPix clusters and writes the tracks to the ESD
//-------------------------------------------------------------------------

#include "AliTracker.h"
#include "AliITSUGeomTGeo.h"
#include <TClonesArray.h>
#include <vector>
#include "AliITSUTrackerSAaux.h"   // Structs and other stuff 
#include "AliITSUMatLUT.h"
#include "AliITSUAux.h"
#include "AliExternalTrackParam.h"
class AliITSUReconstructor;
class AliITSURecoDet;
class AliITSURecoLayer;

class TTree;
class AliCluster;
class AliESDEvent;

using std::vector;

//-------------------------------------------------------------------------
class AliITSUTrackerSA : public AliTracker {
public:
  AliITSUTrackerSA(AliITSUReconstructor* rec=0);
  virtual ~AliITSUTrackerSA();

  // These functions must be implemented 
  Int_t Clusters2Tracks(AliESDEvent *event);
  Int_t PropagateBack(AliESDEvent *event);
  Int_t RefitInward(AliESDEvent *event);
  Int_t LoadClusters(TTree *ct);
  void UnloadClusters();
  AliCluster *GetCluster(Int_t index) const;

  // Possibly, other public functions
  void     Init(AliITSUReconstructor* rec);
  Double_t RefitTrack(AliExternalTrackParam* trc, Int_t clInfo[2*AliITSUAux::kMaxLayers], Double_t rDest, Int_t stopCond);  
  Bool_t   PropagateSeed(AliExternalTrackParam *seed, Double_t xToGo, Double_t mass, Double_t maxStep=1.0, Bool_t matCorr=kTRUE);
  Double_t GetMaterialBudget(const double* pnt0, const double* pnt1, double& x2x0, double& rhol) const;
  Bool_t   GoToEntranceToLayer(AliExternalTrackParam* seed, AliITSURecoLayer* lr, Int_t dir, Bool_t check=kFALSE);
  Bool_t   GoToExitFromLayer(AliExternalTrackParam* seed, AliITSURecoLayer* lr, Int_t dir, Bool_t check=kTRUE);
  Bool_t   TransportToLayerX(AliExternalTrackParam* seed, Int_t lFrom, Int_t lTo, Double_t xStop);  
  Bool_t   TransportToLayer(AliExternalTrackParam* seed, Int_t lFrom, Int_t lTo, Double_t rLim=-1);
  //
protected:
  AliITSUTrackerSA(const AliITSUTrackerSA&);

  void MakeDoublets();
  //  void MakeTriplets();
  void CandidatesTreeTraversal( vector<trackC> &vec, int &iD, int &doubl, int &nCand);
  Bool_t InitTrackParams(trackC &track);
  void CASelection(AliESDEvent *ev);
  void GlobalFit();
  void ChiSquareSelection();
  // Other protected functions
  // (Sorting, labeling, calculations of "roads", etc)  
  static Double_t Curvature(Double_t x1,Double_t y1,Double_t x2,Double_t y2,Double_t x3,Double_t y3);

#ifdef __DEBUG__ 
  void PrintInfo(TString opt);
#endif

private:
  AliITSUTrackerSA &operator=(const AliITSUTrackerSA &tr);

  // Data members

  // classes for interfacing the geometry, materials etc.
  AliITSUReconstructor*           fReconstructor;  // ITS global reconstructor 
  AliITSURecoDet*                 fITS;            // interface to ITS, borrowed from reconstructor
  AliITSUMatLUT*                  fMatLUT;         // material lookup table
  Bool_t                          fUseMatLUT;      //! use material lookup table rather than TGeo
  Double_t                        fCurrMass;       // assumption about particle mass


  // Internal tracker arrays, layers, modules, etc
  vector<itsCluster> fClusters[7];
  TClonesArray *fClustersTC[7];
  vector<nPlets> fDoublets[6];
  Int_t *fIndex[7];
  Int_t fNClusters[7];
  Int_t fNDoublets[6];
  Float_t fPhiCut;
  Float_t fRPhiCut;
  Float_t fZCut;

  //
  static const Double_t           fgkToler;        // tracking tolerance
  //
  ClassDef(AliITSUTrackerSA,1)   //ITSU stand-alone tracker
};

#endif
