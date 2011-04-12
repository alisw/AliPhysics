/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
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

//*****************************************************
//   Class AliEventplane
//   author: Alberica Toia, Johanna Gramling
//*****************************************************
/// A container for the event plane stored in AOD in ESD
 
#include "AliEventplane.h"
#include "TVector2.h"
#include "AliVTrack.h"
#include "TObjArray.h"
#include "TArrayF.h"

ClassImp(AliEventplane)

AliEventplane::AliEventplane() : TNamed("Eventplane", "Eventplane"),
  fQVector(0),
  fQContributionX(0),
  fQContributionY(0),
  fEventplaneQ(0),
  fQsub1(0),
  fQsub2(0),
  fQsubRes(0)
{
  /// constructor
  fQContributionX = new TArrayF();
  fQContributionY = new TArrayF();
}

AliEventplane::AliEventplane(const AliEventplane& ep) : 
  TNamed(),
  fQVector(0),
  fQContributionX(0),
  fQContributionY(0),
  fEventplaneQ(0),
  fQsub1(0),
  fQsub2(0),
  fQsubRes(0)
{
  /// Copy constructor
  ((AliEventplane &) ep).CopyEP(*this);
}

AliEventplane& AliEventplane::operator=(const AliEventplane& ep)
{
  /// Assignment operator
  if (this!=&ep)
    ((AliEventplane &) ep).CopyEP(*this);

  return *this;
}

void AliEventplane::CopyEP(AliEventplane& ep) const
{ // copy function

  AliEventplane& target = (AliEventplane &) ep;
  if (fQContributionX)
      target.fQContributionX = fQContributionX;
  if (fQContributionY)
      target.fQContributionY = fQContributionY;
  if (fEventplaneQ)
      target.fEventplaneQ = fEventplaneQ;
  if (fQVector)
      target.fQVector = dynamic_cast<TVector2*> (fQVector->Clone());
  if (fQsub1)
      target.fQsub1 = dynamic_cast<TVector2*> (fQsub1->Clone());
  if (fQsub2)
      target.fQsub2 = dynamic_cast<TVector2*> (fQsub2->Clone());
  if (fQsubRes)
      target.fQsubRes = fQsubRes;
}

AliEventplane::~AliEventplane()
{
  /// destructor
  if (fQContributionX){
      delete fQContributionX;
      fQContributionX = 0;
  }
  if (fQContributionY){
      delete fQContributionY;
      fQContributionY = 0;
  }
  if (fQVector){
      delete fQVector;
      fQVector = 0;
  }
  if (fQsub1){
      delete fQsub1;
      fQsub1 = 0;
  }
    if (fQsub2){
      delete fQsub2;
      fQsub2 = 0;
  }
}

TVector2* AliEventplane::GetQVector()
{
  return fQVector;
}

Double_t AliEventplane::GetEventplane(const char *x)
{
  TString method = x;
  if(method.CompareTo("Q")==0)      return fEventplaneQ;
  return -1;
}

TVector2* AliEventplane::GetQsub1()
{
  return fQsub1;
}

TVector2* AliEventplane::GetQsub2()
{
  return fQsub2;
}

Double_t AliEventplane::GetQsubRes()
{
  return fQsubRes;
}

Bool_t AliEventplane::IsEventInEventplaneClass(Double_t a, Double_t b, const char *x)
{
  TString method = x;
  if ((method.CompareTo("Q")==0) && (fEventplaneQ >=a && fEventplaneQ < b)) return kTRUE;
  else return kFALSE;
}

Double_t AliEventplane::GetQContributionX(AliVTrack* track)
{ 
  return fQContributionX->GetAt(track->GetID());
}

Double_t AliEventplane::GetQContributionY(AliVTrack* track)
{ 
  return fQContributionY->GetAt(track->GetID());
}
