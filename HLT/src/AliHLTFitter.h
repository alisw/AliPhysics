// @(#) $Id$

#ifndef ALIL3FITTER_H
#define ALIL3FITTER_H

//_____________________________________________________________
// AliHLTFitter
//
// Fit class HLT
//
// Author: Anders Vestbo <mailto:vestbo@fi.uib.no>

class AliHLTTrack;
class AliHLTVertex;
class AliHLTSpacePointData;

class AliHLTFitter {

  public:
  AliHLTFitter();
  AliHLTFitter(AliHLTVertex *vertex,Bool_t vertexconstraint=kTRUE);
  virtual ~AliHLTFitter();
  
  void LoadClusters(Char_t *path,Int_t event=0,Bool_t sp=kFALSE);
  void SortTrackClusters(AliHLTTrack *track) const;
  Int_t FitHelix(AliHLTTrack *track);
  Int_t FitCircle();
  Int_t FitLine();
  void NoVertex() {fVertexConstraint=kFALSE;}
 
 private:
  AliHLTTrack *fTrack; //!                    actual track
  AliHLTVertex *fVertex; //!                  vertex info
  Bool_t fVertexConstraint; //               include vertex constraint
  AliHLTSpacePointData *fClusters[36][6]; //! clusters
  UInt_t fNcl[36][6]; //                     cluster numbers
 
  ClassDef(AliHLTFitter,1) //HLT fit class
};

typedef AliHLTFitter AliL3Fitter; // for backward compatibility

#endif
