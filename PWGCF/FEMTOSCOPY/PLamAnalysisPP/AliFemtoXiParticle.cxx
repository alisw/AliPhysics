
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


#include "AliFemtoXiParticle.h"

ClassImp(AliFemtoXiParticle)

AliFemtoXiParticle::AliFemtoXiParticle() :
  fMomentum(),
  fPt(0),
  fMass(0),
  fDaughterID1(-9999),
  fDaughterID2(-9999),
  fBachID(-9999),
  fXitag(kTRUE),
  fPointing(0)
{
  //Default constructor
}
//_____________________________________________________________________________
/*
AliFemtoXiParticle::AliFemtoXiParticle(const AliFemtoXiParticle &obj) :
  fMomentum(),
  fPt(obj.fPt),
  fMass(obj.fMass),
  fDaughterID1(obj.fDaughterID1),
  fDaughterID2(obj.fDaughterID2),
  fBachID(obj.fBachID),
  fXitag(obj.fXitag),
  fPointing(obj.fPointing)
{
  // copy constructor
}
*/
//_____________________________________________________________________________
AliFemtoXiParticle &AliFemtoXiParticle::operator=(const AliFemtoXiParticle &obj)
{
  //Assignment operator
  if(this == &obj) return *this;
  
  fMomentum = obj.fMomentum;
  fPt = obj.fPt;
  fMass = obj.fMass;
  fDaughterID1 = obj.fDaughterID1;
  fDaughterID2 = obj.fDaughterID2;
  fBachID = obj.fBachID;
  fXitag = obj.fXitag;
  fPointing = obj.fPointing;

 return (*this);
}
//_____________________________________________________________________________
AliFemtoXiParticle::~AliFemtoXiParticle()
{
  // Destructor
}
//_____________________________________________________________________________
