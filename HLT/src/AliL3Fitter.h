// @(#) $Id$

#ifndef ALIL3_Fitter
#define ALIL3_Fitter

#include "AliL3RootTypes.h"

class AliL3Track;
class AliL3Vertex;
class AliL3SpacePointData;

class AliL3Fitter {

 private:
  AliL3Track *fTrack; //!
  AliL3Vertex *fVertex; //!
  Bool_t fVertexConstraint;
  AliL3SpacePointData *fClusters[36][6]; //!
  UInt_t fNcl[36][6];
  
 public:
  AliL3Fitter();
  AliL3Fitter(AliL3Vertex *vertex,Bool_t vertexconstraint=kTRUE);
  virtual ~AliL3Fitter();
  
  void LoadClusters(Char_t *path,Int_t event=0,Bool_t sp=kFALSE);
  void SortTrackClusters(AliL3Track *track);
  Int_t FitHelix(AliL3Track *track);
  Int_t FitCircle();
  Int_t FitLine();
  void NoVertex() {fVertexConstraint=kFALSE;}
  
  ClassDef(AliL3Fitter,1) //HLT fit class
};

#endif
