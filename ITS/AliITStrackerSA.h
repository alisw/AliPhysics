#ifndef ALIITSTRACKERSA_H
#define ALIITSTRACKERSA_H 



#include "AliITSgeomTGeo.h"
#include "AliITStrackerMI.h"

/* Copyright(c) 1998-2003, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


////////////////////////////////////////////////////
//  Stand alone tracker class                     //
//  Origin:  Elisabetta Crescio                   //
//  e-mail:  crescio@to.infn.it                   //
//  adapted for cosmics by A.Dainese              //
////////////////////////////////////////////////////

class AliITSclusterTable;
class AliITStrackSA;
class AliESDVertex;
class AliESDEvent;
class AliITSVertexer;
class TTree;
class TArrayD;

class AliITStrackerSA : public AliITStrackerMI {


 public:

  AliITStrackerSA();
  AliITStrackerSA(const Char_t *geom);
  AliITStrackerSA(const Char_t *geom,AliESDVertex *vert);
  AliITStrackerSA(const Char_t *geom,AliITSVertexer *vertexer);
  virtual ~AliITStrackerSA();  
  virtual Int_t Clusters2Tracks(AliESDEvent *event);
  Int_t FindTracks(AliESDEvent* event, Bool_t useAllClusters=kFALSE);

  AliITStrackV2* FitTrack(AliITStrackSA* tr,Double_t* primaryVertex,Bool_t onePoint=kFALSE);
  void StoreTrack(AliITStrackV2 *t,AliESDEvent *event, Bool_t pureSA) const; 
  Int_t FindTrackLowChiSquare() const;
  Int_t LoadClusters(TTree *cf) {Int_t rc=AliITStrackerMI::LoadClusters(cf); SetClusterTree(cf); return rc;}
  void SetVertex(AliESDVertex *vtx){fVert = vtx;}
  void SetClusterTree(TTree * itscl){fITSclusters = itscl;}

  void SetOutwardFinding() {fInwardFlag=kFALSE;}
  void SetInwardFinding() {fInwardFlag=kTRUE;}
  void SetOuterStartLayer(Int_t osl = 0) {
    if(osl>(AliITSgeomTGeo::GetNLayers()-2)) AliWarning("Minimum Number of layers is 2. OuterStartLayer set to Nlayers-2");
    fOuterStartLayer = TMath::Min(AliITSgeomTGeo::GetNLayers()-2,osl);
  }
  Int_t GetOuterStartLayer() const {return fOuterStartLayer;}
  void SetInnerStartLayer(Int_t isl = 5) {
    if(isl<1) AliWarning("Minimum Number of layers is 2. InnerStartLayer set to 1");
    fInnerStartLayer = TMath::Max(1,isl);
  }
  Int_t GetInnerStartLayer() const {return fInnerStartLayer;}

  void SetSAFlag(Bool_t fl){fITSStandAlone=fl;}  // StandAlone flag setter
  Bool_t GetSAFlag() const {return fITSStandAlone;} // StandAlone flag getter
  void SetFixedWindowSizes(Int_t n=46, Double_t *phi=0, Double_t *lam=0);
  void SetCalculatedWindowSizes(Int_t n=10, Double_t phimin=0.002, Double_t phimax=0.0145, Double_t lambdamin=0.003, Double_t lambdamax=0.008);

  void SetMinNPoints(Int_t np){fMinNPoints=np;}
  Int_t GetMinNPoints() const {return fMinNPoints;}
  void SetMinimumChargeSDDSSD(Double_t minq=0.){fMinQ=minq;}
  enum {kSAflag=0x8000}; //flag to mark clusters used in the SA tracker
 protected:

  //Initialization
  void Init();
  void ResetForFinding();
  void UpdatePoints();
  Bool_t SetFirstPoint(Int_t lay, Int_t clu, Double_t* primaryVertex);
  static Double_t Curvature(Double_t x1,Double_t y1,Double_t x2,Double_t y2,
                     Double_t x3,Double_t y3);

  Double_t ChoosePoint(Double_t p1, Double_t p2, Double_t pp); 

  static Int_t   FindIntersection(Double_t a1, Double_t b1, Double_t c1, Double_t c2, 
                           Double_t& x1,Double_t& y1, Double_t& x2, Double_t& y2);
  static Int_t   FindEquation(Double_t x1, Double_t y1, Double_t x2, Double_t y2, 
                       Double_t x3, Double_t y3,Double_t& a, Double_t& b, 
                       Double_t& c);
 
  Int_t FindLabel(AliITStrackV2* track);
  Int_t SearchClusters(Int_t layer,Double_t phiwindow,Double_t lambdawindow, 
                       AliITStrackSA* trs,Double_t zvertex,Int_t flagp); 

  void GetCoorAngles(AliITSRecPoint* cl,Double_t &phi,Double_t &lambda,Double_t &x,Double_t &y,Double_t &z,Double_t* vertex);
  void GetCoorErrors(AliITSRecPoint* cl,Double_t &sx,Double_t &sy, Double_t &sz);

  AliITSclusterTable* GetClusterCoord(Int_t layer,Int_t n) const {return (AliITSclusterTable*)fCluCoord[layer]->UncheckedAt(n);}
  void RemoveClusterCoord(Int_t layer, Int_t n) {fCluCoord[layer]->RemoveAt(n);fCluCoord[layer]->Compress();}


  Double_t fPhiEstimate; //Estimation of phi angle on next layer
  Bool_t fITSStandAlone; //Tracking is performed in the ITS alone if kTRUE
  Double_t fPoint1[2];   //! coord. of 1-st point to evaluate the curvature
  Double_t fPoint2[2];   //! coord. of 2-nd point to evaluate the curvature
  Double_t fPoint3[2];   //! coord. of 3-rd point to evaluate the curvature
  Double_t fPointc[2];   //! current point coord (for curvature eval.)
  Double_t fLambdac;    //! current value of the Lambda angle in the window
  Double_t fPhic;       //! current value of the Phi angle in the window
  Double_t fCoef1;       //! param. of the equation of the circ. approx a layer
  Double_t fCoef2;       //! param. of the equation of the circ. approx a layer
  Double_t fCoef3;       //! param. of the equation of the circ. approx a layer
  Int_t fNloop;         //  Number of iterqations on phi and lambda windows
  Double_t *fPhiWin;    // phi window sizes
  Double_t *fLambdaWin; // lambda window sizes
  AliESDVertex *fVert;        //! primary vertex
  AliITSVertexer *fVertexer;  //! vertexer 
  TClonesArray *fListOfTracks;   //! container for found tracks 
  TClonesArray *fListOfSATracks; //! container for found SA tracks 
  TTree *fITSclusters;        //! pointer to ITS tree of clusters
  Bool_t fInwardFlag;       // set to kTRUE for inward track finding
  Int_t fOuterStartLayer;     // Outward search for tracks with <6 points: outer layer to start from
  Int_t fInnerStartLayer;     // Inward search for tracks with <6 points: inner layer to start from
  Int_t fMinNPoints;        // minimum number of clusters for a track
  Double_t fMinQ;              // lower cut on cluster charge (SDD and SSD)

  TClonesArray** fCluLayer; //! array with clusters 
  TClonesArray** fCluCoord; //! array with cluster info

 private:
  AliITStrackerSA(const AliITStrackerSA& tracker);
  AliITStrackerSA& operator=(const AliITStrackerSA& source);

  ClassDef(AliITStrackerSA,11)
};

#endif


