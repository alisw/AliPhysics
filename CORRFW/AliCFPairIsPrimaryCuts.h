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
  AliCFPairIsPrimaryCuts(const Char_t* name, const Char_t* title) ;
  AliCFPairIsPrimaryCuts(const AliCFPairIsPrimaryCuts& c) ;
  AliCFPairIsPrimaryCuts& operator=(const AliCFPairIsPrimaryCuts& c) ;
  virtual ~AliCFPairIsPrimaryCuts() {delete fCutNeg; delete fCutPos; }

  Bool_t IsSelected(TObject* obj) ;
  Bool_t IsSelected(TList* /*list*/) {return kTRUE;}
  virtual AliCFTrackIsPrimaryCuts* GetNegCut() const {return fCutNeg;}
  virtual AliCFTrackIsPrimaryCuts* GetPosCut() const {return fCutPos;}

  ClassDef(AliCFPairIsPrimaryCuts,2);

 private :
  AliCFTrackIsPrimaryCuts *fCutNeg ; // isprimary cut on negative daughter
  AliCFTrackIsPrimaryCuts *fCutPos ; // isprimary cut on positive daughter

};

#endif
