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
class AliITSVertex;
class AliITSVertexer;

class AliITStrackerSA : public AliITStrackerV2 {


 public:

  AliITStrackerSA();
  AliITStrackerSA(AliITSgeom *geom,AliITSVertex *vert);
  AliITStrackerSA(AliITSgeom *geom,AliITSVertexer *vertexer);
  AliITStrackerSA(AliITStrackerSA& tracker);
  virtual ~AliITStrackerSA();  
  void     FindTracks(TTree *clusterTree, TTree *out,Int_t evnumber=0,
                      char *opt="6/6");
  void     FindTracks(TTree *clusterTree, AliESD* event, Int_t evnumber=0,
		      char *opt="6/6");
  AliITStrackV2* FitTrack(AliITStrackSA* tr,Double_t* primaryVertex,
                          Double_t *errorprimvert,char *opt="6/6");

  AliITStrackV2* FindTrackLowChiSquare(TObjArray* tracklist, Int_t dim) const;
  void SetVertex(AliITSVertex *vtx){fVert = vtx;}
  void SetWindowSizes(Int_t n=46, Double_t *phi=0, Double_t *lam=0);
  void UseFoundTracksV2(Int_t evnum,TTree* treev2, TTree* clustertree);
  void UseFoundTracksV2(Int_t evnum,AliESD *event, TTree* clustertree);

 protected:

  // copy constructor (NO copy allowed: the constructor is protected
  // to avoid misuse)
  AliITStrackerSA(const AliITStrackerSA& trkr);
  // assignment operator (NO assignment allowed)
  AliITStrackerSA& operator=(const AliITStrackerSA& /* trkr */);
  //Initialization
  void Init();
  Int_t     GetFlagLoadedClusters() const {return fFlagLoad;}
   
  void     ResetForFinding();
  void     SetFlagLoadedClusters(Int_t d) {fFlagLoad=d;}

  void     UpdatePoints();

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
                       AliITStrackSA* trs,Double_t zvertex,Int_t flagp, AliITSclusterTable* table); 
 

  Double_t fPhiEstimate; //Estimation of phi angle on next layer
  Double_t fPoint1[2];   //! coord. of 1-st point to evaluate the curvature
  Double_t fPoint2[2];   //! coord. of 2-nd point to evaluate the curvature
  Double_t fPoint3[2];   //! coord. of 3-rd point to evaluate the curvature
  Double_t fPointc[2];   //! current point coord (for curvature eval.)
  Double_t fLambdac;    //! current value of the Lambda angle in the window
  Double_t fPhic;       //! current value of the Phi angle in the window
  Float_t fCoef1;       //! param. of the equation of the circ. approx a layer
  Float_t fCoef2;       //! param. of the equation of the circ. approx a layer
  Float_t fCoef3;       //! param. of the equation of the circ. approx a layer
  Int_t fNloop;         //  Number of iterqations on phi and lambda windows
  Double_t *fPhiWin;    // phi window sizes
  Double_t *fLambdaWin; // lambda window sizes
  AliITSVertex *fVert;        //! primary vertex
  AliITSVertexer *fVertexer;  //! vertexer 
  AliITSgeom *fGeom;          //! ITS geometry
  Int_t fFlagLoad;            //  flag for loaded clusters (1==already loaded)
  AliITSclusterTable* fTable; //  table with clusters
  ClassDef(AliITStrackerSA,1)
};

#endif

