// @(#) $Id$

#ifndef ALIL3FitterH
#define ALIL3FitterH

//_____________________________________________________________
// AliL3Fitter
//
// Fit class HLT
//
// Author: Anders Vestbo <mailto:vestbo@fi.uib.no>

class AliL3Track;
class AliL3Vertex;
class AliL3SpacePointData;

class AliL3Fitter {

  public:
  AliL3Fitter();
  AliL3Fitter(AliL3Vertex *vertex,Bool_t vertexconstraint=kTRUE);
  virtual ~AliL3Fitter();
  
  void LoadClusters(Char_t *path,Int_t event=0,Bool_t sp=kFALSE);
  void SortTrackClusters(AliL3Track *track) const;
  Int_t FitHelix(AliL3Track *track);
  Int_t FitCircle();
  Int_t FitLine();
  void NoVertex() {fVertexConstraint=kFALSE;}
 
 private:
  AliL3Track *fTrack; //!                    actual track
  AliL3Vertex *fVertex; //!                  vertex info
  Bool_t fVertexConstraint; //               include vertex constraint
  AliL3SpacePointData *fClusters[36][6]; //! clusters
  UInt_t fNcl[36][6]; //                     cluster numbers
 
  ClassDef(AliL3Fitter,1) //HLT fit class
};

#endif
