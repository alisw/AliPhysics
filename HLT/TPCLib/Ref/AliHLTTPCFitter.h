// @(#) $Id$

#ifndef ALIHLTTPC_Fitter
#define ALIHLTTPC_Fitter

#include "AliHLTTPCRootTypes.h"

class AliHLTTPCTrack;
class AliHLTTPCVertex;
class AliHLTTPCSpacePointData;

class AliHLTTPCFitter {

 private:
  AliHLTTPCTrack *fTrack; //!
  AliHLTTPCVertex *fVertex; //!
  Bool_t fVertexConstraint;
  AliHLTTPCSpacePointData *fClusters[36][6]; //!
  UInt_t fNcl[36][6];
  
 public:
  AliHLTTPCFitter();
  AliHLTTPCFitter(AliHLTTPCVertex *vertex,Bool_t vertexconstraint=kTRUE);
  virtual ~AliHLTTPCFitter();
  
  void LoadClusters(Char_t *path,Int_t event=0,Bool_t sp=kFALSE);
  void SortTrackClusters(AliHLTTPCTrack *track);
  void UpdateTrack(AliHLTTPCTrack *track);
  Int_t FitHelix(AliHLTTPCTrack *track);
  Int_t FitCircle();
  Int_t FitLine();
  void NoVertex() {fVertexConstraint=kFALSE;}
  
  ClassDef(AliHLTTPCFitter,1) //HLT fit class
};

#endif
