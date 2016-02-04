#ifndef ALIITSTRACKERUPGRADE_H
#define ALIITSTRACKERUPGRADE_H 

/* Copyright(c) 1998-2003, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
////////////////////////////////////////////////////
//  ITS Upgrade Stand alone tracker class         //
//  Authors: A.Mastroserio C.Terrevoli            //
//  e-mail:annalisa.mastroserio@cern.ch           //
//         cristina.terrevoli@ba.infn.it          //
//  tracks are saved as AliITStrackV2 objects     //
////////////////////////////////////////////////////



#include <TClonesArray.h>
#include "AliITSRecPointU.h"
#include "AliITStrackerMI.h"
#include "AliITSlayerUpgrade.h"
#include "AliITStrackU.h"

class AliITSclusterTable;
class AliITStrackU;
class AliITSsegmentationUpgrade;
class AliESDVertex;
class AliESDEvent;
class AliITSVertexer;
class TTree;
class TArrayD;

class AliITStrackerUpgrade : public AliITStrackerMI {


 public:

  AliITStrackerUpgrade();
  AliITStrackerUpgrade(Int_t nLay);
  virtual ~AliITStrackerUpgrade();  
  virtual Int_t Clusters2Tracks(AliESDEvent *event);  
  Int_t FindTracks(AliESDEvent* event,Bool_t useAllClusters=kFALSE);
  Int_t PropagateBack(AliESDEvent *event);
  Int_t CorrectForPipeMaterial(AliITStrackU *t, TString direction="inward");  
  Int_t RefitInward(AliESDEvent *event);
  AliCluster *GetCluster(Int_t index) const;

  AliITStrackV2* FitTrack(AliITStrackU* tr,Double_t* primaryVertex,Bool_t onePoint=kFALSE);
  void StoreTrack(AliITStrackV2 *t,AliESDEvent *event, Bool_t pureSA) /*const*/; 
  Int_t FindTrackLowChiSquare() const;
  Int_t LoadClusters(TTree *clusTree);
  void SetVertex(AliESDVertex *vtx){fVert = vtx;}
  void SetClusterTree(TTree * itscl){fITSclusters = itscl;}

  void SetOutwardFinding() {fInwardFlag=kFALSE;}
  void SetInwardFinding() {fInwardFlag=kTRUE;}
  void SetOuterStartLayer(Int_t osl = 0) {
    if(osl>(fNLayers-2)) AliWarning("Minimum Number of layers is 2. OuterStartLayer set to Nlayers-2");
    fOuterStartLayer = TMath::Min(fNLayers-2,osl);
  }
  Int_t GetOuterStartLayer() const {return fOuterStartLayer;}
  void SetInnerStartLayer(Int_t isl = 5) {
    if(isl<1) AliWarning("Minimum Number of layers is 2. InnerStartLayer set to 1");
    fInnerStartLayer = TMath::Max(1,isl);
  }
  Int_t GetInnerStartLayer() const {return fInnerStartLayer;}
  void SetSAFlag(Bool_t fl){fITSStandAlone=fl;}  // StandAlone flag setter
  Bool_t GetSAFlag() const {return fITSStandAlone;} // StandAlone flag gett
  void SetFixedWindowSizes(Int_t n=46, Double_t *phi=0, Double_t *lam=0);
  void SetCalculatedWindowSizes(Int_t n=10, Float_t phimin=0.002, Float_t phimax=0.0145, Float_t lambdamin=0.003, Float_t lambdamax=0.008);
  void UnloadClusters();
  void SetMinNPoints(Int_t np){fMinNPoints=np;}
  Int_t GetMinNPoints() const {return fMinNPoints;}
  void SetMinimumChargeSDDSSD(Float_t minq=0.){fMinQ=minq;}
  enum {kSAflag=0x8000}; //flag to mark clusters used in the SA tracker

  void SetNlayers(Int_t nlay) {fNLayers = nlay;}

 protected:

  //Initialization
  void Init();
  void  CreateLayers(Int_t iEvent);
  void ResetForFinding();
  void ResetTrackToFollow(const AliITStrackU &t) {
    fTrackToFollow.~AliITStrackU();
    Bool_t trackMI = kTRUE;
    new(&fTrackToFollow) AliITStrackU(t,trackMI);
  }


  void UpdatePoints();
  Bool_t SetFirstPoint(Int_t lay, Int_t clu, Double_t* primaryVertex);
  static Double_t Curvature(Double_t x1,Double_t y1,Double_t x2,Double_t y2,
			    Double_t x3,Double_t y3);

  Double_t ChoosePoint(Double_t p1, Double_t p2, Double_t pp); 

  static Int_t   FindIntersection(Float_t a1, Float_t b1, Float_t c1, Float_t c2, 
				  Float_t& x1,Float_t& y1, Float_t& x2, Float_t& y2);
  static Int_t   FindEquation(Float_t x1, Float_t y1, Float_t x2, Float_t y2, 
			      Float_t x3, Float_t y3,Float_t& a, Float_t& b, 
			      Float_t& c);
 
  Int_t FindLabel(AliITStrackV2* track) const;
 
  Int_t SearchClusters(Int_t layer,Double_t phiwindow,Double_t lambdawindow, 
                       AliITStrackU* trs,Double_t zvertex,Int_t flagp); 

  void GetCoorAngles(AliITSRecPointU* cl,Double_t &phi,Double_t &lambda,Double_t &x,Double_t &y,Double_t &z,Double_t* vertex);
  void GetCoorErrors(AliITSRecPointU* cl,Float_t &sx,Float_t &sy, Float_t &sz);
  AliITSclusterTable* GetClusterCoord(Int_t layer,Int_t n) const {return (AliITSclusterTable*)fCluCoord[layer]->UncheckedAt(n);}
  void RemoveClusterCoord(Int_t layer, Int_t n) {fCluCoord[layer]->RemoveAt(n);fCluCoord[layer]->Compress();}
  Bool_t RefitAtBase(Double_t x, AliITStrackU *track,
		     const Int_t *clusters);
  Int_t UpdateMI(AliITStrackU* track, const AliITSRecPointU* cl,Double_t chi2,Int_t layer) const;
  Int_t CorrectForLayerMaterial(AliITStrackU *t, Int_t layerindex, Double_t oldGlobXYZ[3], TString direction="inward");
  Double_t GetPredictedChi2MI(AliITStrackU* track, const AliITSRecPointU *cluster,Int_t layer);
  static Int_t GetError(Int_t layer,const AliITSRecPointU *cl,
                        Float_t tgl,Float_t tgphitr,Float_t expQ,
                        Float_t &erry,Float_t &errz,Float_t &covyz,
                        Bool_t addMisalErr=kTRUE);
  
  static Int_t GetErrorOrigRecPoint(const AliITSRecPointU*cl,
                                    Float_t &erry,Float_t &errz,Float_t &covyz);
  //  static void GetNTeor(Int_t layer,const AliITSRecPointU* cl,
  //                       Float_t tgl,Float_t tgphitr,
  //                       Float_t &ny,Float_t &nz);

  static const Int_t fgMaxNLayer = 8; //max number of layers in ITSUpgrade
  Int_t fNLayers;//number of layer in ITSUpgrade
  Double_t fPhiEstimate; //Estimation of phi angle on next layer
  Bool_t fITSStandAlone; //Tracking is performed in the ITS alone if kTRUE
  Float_t fPoint1[2];   //! coord. of 1-st point to evaluate the curvature
  Float_t fPoint2[2];   //! coord. of 2-nd point to evaluate the curvature
  Float_t fPoint3[2];   //! coord. of 3-rd point to evaluate the curvature
  Float_t fPointc[2];   //! current point coord (for curvature eval.)
  Double_t fLambdac;    //! current value of the Lambda angle in the window
  Double_t fPhic;       //! current value of the Phi angle in the window
  Float_t fCoef1;       //! param. of the equation of the circ. approx a layer
  Float_t fCoef2;       //! param. of the equation of the circ. approx a layer
  Float_t fCoef3;       //! param. of the equation of the circ. approx a layer
  Int_t fNloop;         //  Number of iterqations on phi and lambda windows
  Double_t *fPhiWin;    // phi window sizes
  Double_t *fLambdaWin; // lambda window sizes
  AliESDVertex *fVert;        //! primary vertex
  AliITSVertexer *fVertexer;  //! vertexer 
  TClonesArray *fListOfTracks;   //! container for found tracks 
  TClonesArray *fListOfUTracks; //! container for found SA tracks 
  TTree *fITSclusters;        //! pointer to ITS tree of clusters
  Bool_t fInwardFlag;       // set to kTRUE for inward track finding
  Int_t fOuterStartLayer;     // Outward search for tracks with <6 points: outer layer to start from
  Int_t fInnerStartLayer;     // Inward search for tracks with <6 points: inner layer to start from
  Int_t fMinNPoints;        // minimum number of clusters for a track
  Float_t fMinQ;              // lower cut on cluster charge (SDD and SSD)
  AliITStrackU fTrackToFollow;          

  AliITSlayerUpgrade** fLayers;
  AliITSsegmentationUpgrade *fSegmentation;

  TClonesArray** fCluLayer; //! array with clusters 
  TClonesArray** fCluCoord; //! array with cluster info

  private:
  AliITStrackerUpgrade(const AliITStrackerUpgrade& tracker);
  AliITStrackerUpgrade &operator=(const AliITStrackerUpgrade &tr);

  ClassDef(AliITStrackerUpgrade,2)
    };

#endif


