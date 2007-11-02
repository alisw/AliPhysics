// @(#) $Id$
// Original: AliHLTTrackMerger.h,v 1.6 2005/04/19 04:29:01 cvetan 
#ifndef ALIHLTTPCTRACKMERGER_H
#define ALIHLTTPCTRACKMERGER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTTPCTrackMerger.h
    @author Uli Frankenfeld, maintained by Matthias Richter
    @date   
    @brief  The HLT TPC track segment merger
*/

#ifndef __CINT__ 
#include "AliHLTTPCMerger.h"
#endif

#include "AliHLTTPCRootTypes.h"

class AliHLTTPCMerger;

/** 
 * @class AliHLTTPCTrackMerger
 * The HLTTPC track segment merger
 *
 *   This class is responsible for the merging of the HLT tracks
 *   between TPC sectors and readout patches
 */
class AliHLTTPCTrackMerger : public AliHLTTPCMerger {

 public:
  AliHLTTPCTrackMerger();
  AliHLTTPCTrackMerger(Int_t nsubsectors);
  /** destructor */
  virtual ~AliHLTTPCTrackMerger();

  void SetRows(Int_t *row);
  void InitSector(Int_t sector,Int_t subsector);
  void SlowMerge();
  void Merge();  //Loop over tracks from different subsectors
  void InterMerge();

 private:
  /** copy constructor prohibited */
  AliHLTTPCTrackMerger(const AliHLTTPCTrackMerger&);
  /** assignment operator prohibited */
  AliHLTTPCTrackMerger& operator=(const AliHLTTPCTrackMerger&);

  Int_t fSubSector;//Index of the readout patch inside given TPC sector
  Int_t fNSubSector;//Number of readout patches inside given TPC sector
  Int_t *fRowMin;//!
  Int_t *fRowMax;//!
  Bool_t fSlow;//Slow or fast merging
  void SlowMerge(AliHLTTPCTrackArray *mergedtrack,AliHLTTPCTrackArray *tracksin,AliHLTTPCTrackArray *tracksout,Double_t xval);
  Int_t Merge(AliHLTTPCTrackArray *mergedtrack,AliHLTTPCTrackArray *tracksin,AliHLTTPCTrackArray *tracksout);
  
  ClassDef(AliHLTTPCTrackMerger,1) //Track merging class 
};

#endif
