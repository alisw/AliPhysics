
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
#include "AliFemtoLambdaParticle.h"

ClassImp(AliFemtoLambdaParticle)

AliFemtoLambdaParticle::AliFemtoLambdaParticle() :
  fMomentum(),
  fMomentumMC(),
  fMomentumMCMother(),
  fPDGCode(0),
  fPDGCodeMother(0),
  fPt(0),
  fMass(0),
  fDaughterID1(-9999),
  fDaughterID2(-9999),
  fV0tag(kTRUE),
  fPointing(0),
  fMomentumPosDaughter(),
  fMomentumNegDaughter(),
  fPhiPosdaughter(0),
  fPhiStarPosdaughter(),
  fEtaPosdaughter(0),
  fPhiNegdaughter(0),
  fPhiStarNegdaughter(),
  fEtaNegdaughter(0),
  fReal(kFALSE)
{
  //Default constructor
}
//_____________________________________________________________________________
//AliFemtoLambdaParticle::AliFemtoLambdaParticle(const AliFemtoLambdaParticle &obj) :
//  fMomentum(),
//  fMomentumMC(),
//  fMomentumMCMother(),
//  fPDGCode(obj.fPDGCode),
//  fPDGCodeMother(obj.fPDGCodeMother),
//  fPt(obj.fPt),
//  fMass(obj.fMass),
//  fDaughterID1(obj.fDaughterID1),
//  fDaughterID2(obj.fDaughterID2),
//  fV0tag(obj.fV0tag),
//  fPointing(obj.fPointing),
//  fPhiPosdaughter(obj.fPhiPosdaughter),
//  fEtaPosdaughter(obj.fEtaPosdaughter),
//  fPhiNegdaughter(obj.fPhiNegdaughter),
//  fEtaNegdaughter(obj.fEtaNegdaughter)
//{
//  // copy constructor
//}

//_____________________________________________________________________________
AliFemtoLambdaParticle &AliFemtoLambdaParticle::operator=(const AliFemtoLambdaParticle &obj)
{
  //Assignment operator
  if(this == &obj) return *this;
  
  fMomentum = obj.fMomentum;
  fMomentumMC = obj.fMomentumMC;
  fMomentumMCMother = obj.fMomentumMCMother;
  fPDGCode = obj.fPDGCode;
  fPDGCodeMother = obj.fPDGCodeMother;
  fPt = obj.fPt;
  fMass = obj.fMass;
  fDaughterID1 = obj.fDaughterID1;
  fDaughterID2 = obj.fDaughterID2;
  fV0tag = obj.fV0tag;
  fPointing = obj.fPointing;

  fMomentumPosDaughter = obj.fMomentumPosDaughter;
  fMomentumNegDaughter = obj.fMomentumNegDaughter;
  fPhiPosdaughter = obj.fPhiPosdaughter;
  fEtaPosdaughter = obj.fEtaPosdaughter;
  fPhiNegdaughter = obj.fPhiNegdaughter;
  fEtaNegdaughter = obj.fEtaNegdaughter;
  fReal = obj.fReal;


  for(int i=0;i<9;i++)
   {
     fPhiStarPosdaughter[i] = obj.fPhiStarPosdaughter[i];
     fPhiStarNegdaughter[i] = obj.fPhiStarNegdaughter[i];
     fPositionPosTPC[i] = obj.fPositionPosTPC[i];
     fPositionNegTPC[i] = obj.fPositionNegTPC[i];
   }
  
  
 return (*this);
}
//_____________________________________________________________________________
AliFemtoLambdaParticle::~AliFemtoLambdaParticle()
{
  // Destructor
}
//_____________________________________________________________________________
