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


#include "AliCFPairIsPrimaryCuts.h"
#include "AliCFPair.h"

ClassImp(AliCFPairIsPrimaryCuts)

//__________________________________________________________________________________
AliCFPairIsPrimaryCuts::AliCFPairIsPrimaryCuts() :
  AliCFCutBase(),
  fCutNeg(new AliCFTrackIsPrimaryCuts()),
  fCutPos(new AliCFTrackIsPrimaryCuts())
{
  //
  // Default constructor
  //
}
//__________________________________________________________________________________
AliCFPairIsPrimaryCuts::AliCFPairIsPrimaryCuts(Char_t* name, Char_t* title) :
  AliCFCutBase(name,title),
  fCutNeg(new AliCFTrackIsPrimaryCuts(name,title)),
  fCutPos(new AliCFTrackIsPrimaryCuts(name,title))
{
  //
  // Named Constructor
  //
}
//__________________________________________________________________________________
AliCFPairIsPrimaryCuts::AliCFPairIsPrimaryCuts(const AliCFPairIsPrimaryCuts& c) :
  AliCFCutBase(c),
  fCutNeg(c.fCutNeg),
  fCutPos(c.fCutPos)
{
  //
  // copy constructor
  //
}
//__________________________________________________________________________________
AliCFPairIsPrimaryCuts& AliCFPairIsPrimaryCuts::operator=(const AliCFPairIsPrimaryCuts& c)
{
  //
  // Assignment operator
  //
  if (this != &c) {
    AliCFCutBase::operator=(c) ;
    fCutNeg = c.fCutNeg ;
    fCutPos = c.fCutPos ;
  }
  return *this;
}

//__________________________________________________________________________________
Bool_t AliCFPairIsPrimaryCuts::IsSelected(TObject* obj) {
  //
  // loops over decisions of single cuts and returns if the track is accepted
  //

  if (!obj) return kFALSE ;
  TString className(obj->ClassName());
  if (className.CompareTo("AliCFPair") != 0) {
    Error("IsSelected","obj must point to a AliCFPair !");
    return kFALSE ;
  }

  AliCFPair* pair = dynamic_cast<AliCFPair*>(obj);

  AliVParticle* tneg = pair->GetNeg();
  AliVParticle* tpos = pair->GetPos();
  if (!tneg || !tpos) return kFALSE ;

  if ( ! fCutNeg->IsSelected((TObject*)tneg) || ! fCutPos->IsSelected((TObject*)tpos) ) return kFALSE ;

  return kTRUE ;
}
