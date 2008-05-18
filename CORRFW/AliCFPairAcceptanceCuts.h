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


///////////////////////////////////////////////////////////////////////////
//          ----   CORRECTION FRAMEWORK   ----
// Class to cut on the number of AliTrackReference's 
// for each detector. Applies on pair of tracks (AliCFPair)
///////////////////////////////////////////////////////////////////////////
// author : R. Vernet (renaud.vernet@cern.ch)
///////////////////////////////////////////////////////////////////////////


#ifndef ALICFPAIRACCEPTANCECUTS_H
#define ALICFPAIRACCEPTANCECUTS_H

#include "AliCFAcceptanceCuts.h"
#include "AliCFCutBase.h"

class AliMCEvent;
class TBits;

class AliCFPairAcceptanceCuts : public AliCFCutBase
{
 public :
  AliCFPairAcceptanceCuts() ;
  AliCFPairAcceptanceCuts(const Char_t* name, const Char_t* title) ;
  AliCFPairAcceptanceCuts(const AliCFPairAcceptanceCuts& c) ;
  AliCFPairAcceptanceCuts& operator=(const AliCFPairAcceptanceCuts& c) ;
  virtual ~AliCFPairAcceptanceCuts() {delete fCutNeg; delete fCutPos; }
  Bool_t IsSelected(TObject* obj) ;
  Bool_t IsSelected(TList* /*list*/) {return kTRUE;}
  virtual void SetEvtInfo(TObject *mcInfo) ;
  virtual void SetMinNHitITS  (Int_t nHitNeg, Int_t nHitPos) {fCutNeg->SetMinNHitITS (nHitNeg); fCutPos->SetMinNHitITS (nHitPos);}
  virtual void SetMinNHitTPC  (Int_t nHitNeg, Int_t nHitPos) {fCutNeg->SetMinNHitTPC (nHitNeg); fCutPos->SetMinNHitTPC (nHitPos);}
  virtual void SetMinNHitTRD  (Int_t nHitNeg, Int_t nHitPos) {fCutNeg->SetMinNHitTRD (nHitNeg); fCutPos->SetMinNHitTRD (nHitPos);}
  virtual void SetMinNHitTOF  (Int_t nHitNeg, Int_t nHitPos) {fCutNeg->SetMinNHitTOF (nHitNeg); fCutPos->SetMinNHitTOF (nHitPos);}
  virtual void SetMinNHitMUON (Int_t nHitNeg, Int_t nHitPos) {fCutNeg->SetMinNHitMUON(nHitNeg); fCutPos->SetMinNHitMUON(nHitPos);}

  enum {
    kNCuts=2
  };

 protected:
  AliMCEvent          *fMCInfo ; // global event information
  AliCFAcceptanceCuts *fCutNeg ; // acceptance cut on negative daughter
  AliCFAcceptanceCuts *fCutPos ; // acceptance cut on positive daughter
  TBits               *fBitmap ; // cut bitmap    

 private:
  void SelectionBitMap(TObject* obj);

  ClassDef(AliCFPairAcceptanceCuts,1);
};

#endif
