/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

///////////////////////////////////////////////
// Class to handle track quality in track pairs
// The track pair object to use is AliCFPair
// author : renaud.vernet@cern.ch
///////////////////////////////////////////////



#ifndef ALICFPAIRQUALITYCUTS_H
#define ALICFPAIRQUALITYCUTS_H

#include "AliCFCutBase.h"
#include "AliCFTrackQualityCuts.h"

class AliESDEvent;

class AliCFPairQualityCuts : public AliCFCutBase
{
 public :
  AliCFPairQualityCuts() ;
  AliCFPairQualityCuts(Char_t* name, Char_t* title) ;
  AliCFPairQualityCuts(const AliCFPairQualityCuts& c) ;
  AliCFPairQualityCuts& operator=(const AliCFPairQualityCuts& c) ;
  virtual ~AliCFPairQualityCuts() {delete fCutNeg; delete fCutPos; }

  Bool_t IsSelected(TObject* obj) ; 
  Bool_t IsSelected(TList* /*list*/) {return kTRUE;}
  virtual void SetMinNClusterTPC (UInt_t nClusNeg, UInt_t nClusPos) 
  {fCutNeg->SetMinNClusterTPC(nClusNeg); fCutPos->SetMinNClusterTPC(nClusPos);}
  virtual void SetMinNClusterITS (UInt_t nClusNeg, UInt_t nClusPos) 
  {fCutNeg->SetMinNClusterITS(nClusNeg); fCutPos->SetMinNClusterITS(nClusPos);}
  virtual void SetMaxChi2PerClusterTPC(Double32_t chi2Neg, Double32_t chi2Pos) 
  {fCutNeg->SetMaxChi2PerClusterTPC(chi2Neg); fCutPos->SetMaxChi2PerClusterTPC(chi2Pos);}
  virtual void SetMaxChi2PerClusterITS(Double32_t chi2Neg, Double32_t chi2Pos) 
  {fCutNeg->SetMaxChi2PerClusterITS(chi2Neg); fCutPos->SetMaxChi2PerClusterITS(chi2Pos);}
  virtual void SetMaxCovDiagonalElements(Double32_t* neg/*[5]*/, Double32_t* pos/*[5]*/) {
    fCutNeg->SetMaxCovDiagonalElements(neg[0],neg[1],neg[2],neg[3],neg[4]); 
    fCutPos->SetMaxCovDiagonalElements(pos[0],pos[1],pos[2],pos[3],pos[4]); }
  virtual void SetStatus(ULong_t statusNeg, ULong_t statusPos) 
  {fCutNeg->SetStatus(statusNeg); fCutPos->SetStatus(statusPos);}

 private :
  AliCFTrackQualityCuts *fCutNeg ; // quality cut on negative daughter
  AliCFTrackQualityCuts *fCutPos ; // quality cut on positive daughter

  ClassDef(AliCFPairQualityCuts,2);
};

#endif
