
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
  fEtaPosdaughter(0),
  fPhiNegdaughter(0),
  fEtaNegdaughter(0),
  fReal(kFALSE)
{
  //Default constructor
}
//_____________________________________________________________________________
/*
AliFemtoLambdaParticle::AliFemtoLambdaParticle(const AliFemtoLambdaParticle &obj) :
  fMomentum(),
  fMomentumMC(),
  fMomentumMCMother(),
  fPDGCode(obj.fPDGCode),
  fPDGCodeMother(obj.fPDGCodeMother),
  fPt(obj.fPt),
  fMass(obj.fMass),
  fDaughterID1(obj.fDaughterID1),
  fDaughterID2(obj.fDaughterID2),
  fV0tag(obj.fV0tag),
  fPointing(obj.fPointing),
  fPhiPosdaughter(obj.fPhiPosdaughter),
  fEtaPosdaughter(obj.fEtaPosdaughter),
  fPhiNegdaughter(obj.fPhiNegdaughter),
  fEtaNegdaughter(obj.fEtaNegdaughter)
{
  // copy constructor
}
*/
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
     for(int j=0;j<3;j++)
       {
	 fPosDaughPosTPC[i][j] = obj.fPosDaughPosTPC[i][j];
	 fNegDaughPosTPC[i][j] = obj.fNegDaughPosTPC[i][j];
       }
   }
  
  
 return (*this);
}
//_____________________________________________________________________________
AliFemtoLambdaParticle::~AliFemtoLambdaParticle()
{
  // Destructor
}
//_____________________________________________________________________________
