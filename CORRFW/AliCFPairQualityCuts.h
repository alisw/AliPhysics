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
  AliCFPairQualityCuts(const Char_t* name, const Char_t* title) ;
  AliCFPairQualityCuts(const AliCFPairQualityCuts& c) ;
  AliCFPairQualityCuts& operator=(const AliCFPairQualityCuts& c) ;
  virtual ~AliCFPairQualityCuts() {delete fCutNeg; delete fCutPos; }

  Bool_t IsSelected(TObject* obj) ; 
  Bool_t IsSelected(TList* /*list*/) {return kTRUE;}

  virtual AliCFTrackQualityCuts* GetNegCut() const {return fCutNeg;}
  virtual AliCFTrackQualityCuts* GetPosCut() const {return fCutPos;}

 private :
  AliCFTrackQualityCuts *fCutNeg ; // quality cut on negative daughter
  AliCFTrackQualityCuts *fCutPos ; // quality cut on positive daughter

  ClassDef(AliCFPairQualityCuts,2);
};

#endif
