#ifndef ALIITSTRACKERSA_H
#define ALIITSTRACKERSA_H 



#include "AliITStrackerV2.h"

/* Copyright(c) 1998-2003, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////////
//  Stand alone tracker class                     //
//  Origin:  Elisabetta Crescio                   //
//  e-mail:  crescio@to.infn.it                   //
////////////////////////////////////////////////////

class AliITSclusterTable;
class AliITStrackSA;
class AliESDVertex;
class AliITSVertexer;
class TTree;

class AliITStrackerSA : public AliITStrackerV2 {


 public:

  AliITStrackerSA();
  AliITStrackerSA(AliITSgeom *geom);
  AliITStrackerSA(AliITSgeom *geom,AliESDVertex *vert);
  AliITStrackerSA(AliITSgeom *geom,AliITSVertexer *vertexer);
  AliITStrackerSA(AliITStrackerSA& tracker);
  virtual ~AliITStrackerSA();  
  virtual Int_t Clusters2Tracks(AliESD *event){Int_t rc = AliITStrackerV2::Clusters2Tracks(event); if(!rc) rc=FindTracks(event); return rc;}
  Int_t FindTracks(AliESD* event);
  void  FindTracks(TTree *out,Int_t evnumber=0);
  AliITStrackV2* FitTrack(AliITStrackSA* tr,Double_t* primaryVertex,
                          Double_t *errorprimvert);

  AliITStrackV2* FindTrackLowChiSquare(TObjArray* tracklist, Int_t dim) const;
  Int_t LoadClusters(TTree *cf) {Int_t rc=AliITStrackerV2::LoadClusters(cf); SetClusterTree(cf);SetSixPoints(kTRUE); return rc;}
  void SetVertex(AliESDVertex *vtx){fVert = vtx;}
  void SetClusterTree(TTree * itscl){fITSclusters = itscl;}
  void SetSixPoints(Bool_t sp = kTRUE){fSixPoints = sp;}
  Bool_t GetSixPoints() const {return fSixPoints;}
  void SetWindowSizes(Int_t n=46, Double_t *phi=0, Double_t *lam=0);
  void UseFoundTracksV2(Int_t evnum,TTree* treev2);
  void UseFoundTracksV2(AliESD *event);

  enum {kSAflag=0x8000}; //flag to mark clusters used in the SA tracker

 protected:

  // copy constructor (NO copy allowed: the constructor is protected
  // to avoid misuse)
  AliITStrackerSA(const AliITStrackerSA& trkr);
  // assignment operator (NO assignment allowed)
  AliITStrackerSA& operator=(const AliITStrackerSA& /* trkr */);
  //Initialization
  void Init();
  void ResetForFinding();
  void UpdatePoints();

  static Double_t Curvature(Double_t x1,Double_t y1,Double_t x2,Double_t y2,
                     Double_t x3,Double_t y3);

  Double_t ChoosePoint(Double_t p1, Double_t p2, Double_t pp); 

  static Int_t   FindIntersection(Float_t a1, Float_t b1, Float_t c1, Float_t c2, 
                           Float_t& x1,Float_t& y1, Float_t& x2, Float_t& y2);
  static Int_t   FindEquation(Float_t x1, Float_t y1, Float_t x2, Float_t y2, 
                       Float_t x3, Float_t y3,Float_t& a, Float_t& b, 
                       Float_t& c);
 
  static Int_t FindLabel(Int_t l1, Int_t l2, Int_t l3, Int_t l4, Int_t l5, Int_t l6);
  static Int_t Label(Int_t gl1, Int_t gl2, Int_t gl3, Int_t gl4, Int_t gl5, 
              Int_t gl6,Int_t gl7, Int_t gl8, Int_t gl9, Int_t gl10,Int_t gl11,
              Int_t gl12, Int_t gl13, Int_t gl14,Int_t gl15, Int_t gl16, 
              Int_t gl17, Int_t gl18, Int_t numberofpoints=6);
 
  Int_t SearchClusters(Int_t layer,Double_t phiwindow,Double_t lambdawindow, 
                       AliITStrackSA* trs,Double_t zvertex,Int_t flagp); 
 
  Double_t fPhiEstimate; //Estimation of phi angle on next layer
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
  AliITSgeom *fGeom;          //! ITS geometry
  TObjArray *fListOfTracks;   //! container for found tracks 
  TTree *fITSclusters;        //! pointer to ITS tree of clusters
  Bool_t fSixPoints;          // If true 6/6 points are required (default). 5/6 otherwise
  AliITSclusterTable* fTable; //  table with clusters
  ClassDef(AliITStrackerSA,1)
};

#endif

