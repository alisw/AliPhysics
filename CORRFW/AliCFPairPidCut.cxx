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

////////////////////////////////////////////////
// Class to define PID cuts on a pair of tracks
// The track pair object to use is AliCFPair
// author : renaud.vernet@cern.ch
////////////////////////////////////////////////

#include "AliCFPairPidCut.h"
#include "AliCFPair.h"

ClassImp(AliCFPairPidCut)

//__________________________________________________________________________________
AliCFPairPidCut::AliCFPairPidCut() :
  AliCFCutBase(),
  fCutNeg(new AliCFTrackCutPid()),
  fCutPos(new AliCFTrackCutPid())
{
  //
  // Default constructor
  //
}
//__________________________________________________________________________________
AliCFPairPidCut::AliCFPairPidCut(const Char_t* name, const Char_t* title) :
  AliCFCutBase(name,title),
  fCutNeg(new AliCFTrackCutPid(name,title)),
  fCutPos(new AliCFTrackCutPid(name,title))
{
  //
  // Constructor
  //
}
//__________________________________________________________________________________
AliCFPairPidCut::AliCFPairPidCut(const AliCFPairPidCut& c) :
  AliCFCutBase(c),
  fCutNeg(c.fCutNeg),
  fCutPos(c.fCutPos)
{
  //
  // copy constructor
  //
}
//__________________________________________________________________________________
AliCFPairPidCut& AliCFPairPidCut::operator=(const AliCFPairPidCut& c)
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
Bool_t AliCFPairPidCut::IsSelected(TObject* obj) {
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
  if ( ! fCutNeg->IsSelected(tneg) || ! fCutPos->IsSelected(tpos) ) return kFALSE ;
  return kTRUE ;
}
