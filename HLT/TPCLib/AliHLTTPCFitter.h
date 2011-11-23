// $Id$
// Original: AliHLTFitter.h,v 1.7 2004/07/05 09:02:18 loizides 

#ifndef ALIHLTTPCFITTER_H
#define ALIHLTTPCFITTER_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/// @file   AliHLTTPCFitter.h
/// @author Anders Vestbo, maintained by Matthias Richter
/// @date   
/// @brief  Fit class HLT for helix
///

class AliHLTTPCTrack;
class AliHLTTPCVertex;
struct AliHLTTPCSpacePointData;

/** 
 * @class AliHLTTPCFitter
 * Fit class HLT for helix
 *
 * @ingroup alihlt_tpc
 */
class AliHLTTPCFitter {

  public:
  AliHLTTPCFitter();
  AliHLTTPCFitter(AliHLTTPCVertex *vertex,Bool_t vertexconstraint=kTRUE);
  virtual ~AliHLTTPCFitter();
  
  void SortTrackClusters(AliHLTTPCTrack *track) const;
  Int_t FitHelix(AliHLTTPCTrack *track);
  Int_t FitCircle();
  Int_t FitLine();
  void NoVertex() {fVertexConstraint=kFALSE;}
 
 private:
  /** copy constructor prohibited */
  AliHLTTPCFitter(const AliHLTTPCFitter& src);
  /** assignment operator prohibited */
  AliHLTTPCFitter& operator=(const AliHLTTPCFitter& src);

  AliHLTTPCTrack *fTrack; //!                    actual track
  AliHLTTPCVertex *fVertex; //!                  vertex info
  Bool_t fVertexConstraint; //               include vertex constraint
  AliHLTTPCSpacePointData *fClusters[36][6]; //! clusters
  UInt_t fNcl[36][6]; //                     cluster numbers
 
  ClassDef(AliHLTTPCFitter,0) //HLT fit class
};

#endif
