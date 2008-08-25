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

///////////////////////////////////////////////////////
// Class to handle primariness criteria in track pairs
// The track pair object to use is AliCFPair
// author : renaud.vernet@cern.ch
///////////////////////////////////////////////////////

#ifndef ALICFPAIRISPRIMARYCUTS_H
#define ALICFPAIRISPRIMARYCUTS_H

#include "AliCFCutBase.h"
#include "AliCFTrackIsPrimaryCuts.h"

class AliESDEvent;

class AliCFPairIsPrimaryCuts : public AliCFCutBase
{
 public :
  AliCFPairIsPrimaryCuts() ;
  AliCFPairIsPrimaryCuts(Char_t* name, Char_t* title) ;
  AliCFPairIsPrimaryCuts(const AliCFPairIsPrimaryCuts& c) ;
  AliCFPairIsPrimaryCuts& operator=(const AliCFPairIsPrimaryCuts& c) ;
  virtual ~AliCFPairIsPrimaryCuts() {delete fCutNeg; delete fCutPos; }

  Bool_t IsSelected(TObject* obj) ;
  Bool_t IsSelected(TList* /*list*/) {return kTRUE;}
  virtual void SetMaxNSigmaToVertex(Double32_t neg, Double32_t pos) 
  {fCutNeg->SetMaxNSigmaToVertex(neg); fCutPos->SetMaxNSigmaToVertex(pos);}
  void SetRequireSigmaToVertex(Bool_t b1, Bool_t b2)
  {fCutNeg->SetRequireSigmaToVertex(b1); fCutPos->SetRequireSigmaToVertex(b2);}
  void SetAcceptKinkDaughters(Bool_t b1, Bool_t b2)
  {fCutNeg->SetAcceptKinkDaughters(b1); fCutPos->SetAcceptKinkDaughters(b2);}
  void SetAODType(Char_t typeNeg, Char_t typePos)
  {fCutNeg->SetAODType(typeNeg); fCutPos->SetAODType(typePos);}

  ClassDef(AliCFPairIsPrimaryCuts,2);

 private :
  AliCFTrackIsPrimaryCuts *fCutNeg ; // isprimary cut on negative daughter
  AliCFTrackIsPrimaryCuts *fCutPos ; // isprimary cut on positive daughter

};

#endif
