 /**************************************************************************
 * Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
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
#include <iostream>
#include "AliEmcalTrackSelResultPtr.h"
#include "AliVTrack.h"

/// \cond CLASSIMP
ClassImp(AliEmcalTrackSelResultPtr)
/// \endcond

AliEmcalTrackSelResultPtr::AliEmcalTrackSelResultPtr() :
  TObject(),
  fTrack(nullptr),
  fSelectionResult(false),
  fFlag(0)
{

}

AliEmcalTrackSelResultPtr::AliEmcalTrackSelResultPtr(AliVTrack *trk, Bool_t selectionStatus, ULong_t flag) :
  TObject(),
  fTrack(trk),
  fSelectionResult(selectionStatus),
  fFlag(flag)
{

}

AliEmcalTrackSelResultPtr::AliEmcalTrackSelResultPtr(const AliEmcalTrackSelResultPtr &ref) :
  TObject(ref),
  fTrack(ref.fTrack),
  fSelectionResult(ref.fSelectionResult),
  fFlag(ref.fFlag)
{

}

AliEmcalTrackSelResultPtr::AliEmcalTrackSelResultPtr(AliEmcalTrackSelResultPtr &&ref) :
  TObject(ref),
  fTrack(ref.fTrack),
  fSelectionResult(ref.fSelectionResult),
  fFlag(ref.fFlag)
{
  ref.fTrack = nullptr;
}

AliEmcalTrackSelResultPtr &AliEmcalTrackSelResultPtr::operator =(const AliEmcalTrackSelResultPtr &ref){
  TObject::operator =(ref);
  if(this != &ref){
    fTrack = ref.fTrack;
    fSelectionResult = ref.fSelectionResult;
    fFlag = ref.fFlag;
  }
  return *this;
}

AliEmcalTrackSelResultPtr &AliEmcalTrackSelResultPtr::operator =(AliEmcalTrackSelResultPtr &&ref){
  TObject::operator =(ref);
  if(this != &ref){
    fTrack = ref.fTrack;
    fSelectionResult = ref.fSelectionResult;
    fFlag = ref.fFlag;

    delete ref.fTrack;
  }
  return *this;
}

Bool_t AliEmcalTrackSelResultPtr::operator ==(const AliEmcalTrackSelResultPtr &other) const {
  return fTrack == other.fTrack;
}

Bool_t AliEmcalTrackSelResultPtr::operator <(const AliEmcalTrackSelResultPtr &other) const {
  return fTrack < other.fTrack;
}

Bool_t AliEmcalTrackSelResultPtr::IsEqual(const TObject *o) const {
  const AliEmcalTrackSelResultPtr *otherobj = static_cast<const AliEmcalTrackSelResultPtr *>(o);
  if(!otherobj) return false;
  return *this == *otherobj;
}

Int_t AliEmcalTrackSelResultPtr::Compare(const TObject *o) const {
  const AliEmcalTrackSelResultPtr *otherobj = static_cast<const AliEmcalTrackSelResultPtr *>(o);
  if(!otherobj) return 1;
  if (*this == *otherobj) return 0;
  if (*this < *otherobj) return -1;
  return 1;
}

AliVTrack *AliEmcalTrackSelResultPtr::operator *() const {
  return fTrack;
}

AliVTrack *AliEmcalTrackSelResultPtr::operator ->() const {
  return fTrack;
}

void AliEmcalTrackSelResultPtr::PrintStream(std::ostream &stream) const {
  stream  << "Track selection result for track with address " << fTrack
          << ": Selection status: " << (fSelectionResult ? "true" : "false")
          << ", flag: " << fFlag;
}

std::ostream &operator<<(std::ostream &stream, const AliEmcalTrackSelResultPtr &o){
  o.PrintStream(stream);
  return stream;
}
