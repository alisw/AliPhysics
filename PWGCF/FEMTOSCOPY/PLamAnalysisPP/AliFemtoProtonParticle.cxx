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


#include "AliFemtoProtonParticle.h"

ClassImp(AliFemtoProtonParticle)

AliFemtoProtonParticle::AliFemtoProtonParticle() :
  fMomentum(),
  fMomentumMC(),
  fMomentumMCMother(),
  fMomentumMCMotherParton(),
  fPDGCode(0),
  fPDGCodeMother(0),
  fPDGCodePartonMother(0),
  fPartonMotherLabel(0),
  fPt(0),
  fID(0),
  fPhi(0),
  fPhistar(),
  fEta(0),
  fReal(kFALSE),
  fProtonTag(kTRUE)
{

  fMomentum.SetXYZ(0.,0.,0.);
  fMomentumMC.SetXYZ(0.,0.,0.);
  fMomentumMCMother.SetXYZ(0.,0.,0.);
  fMomentumMCMotherParton.SetXYZ(0.,0.,0.);

  //Default constructor
}
//_____________________________________________________________________________
//AliFemtoProtonParticle::AliFemtoProtonParticle(const AliFemtoProtonParticle &obj) :
//  fMomentum(),
//  fMomentumMC(),
//  fMomentumMCMother(),
//  fPDGCode(obj.fPDGCode),
//  fPDGCodeMother(obj.fPDGCodeMother),
//  fPt(obj.fPt),
//  fID(obj.fID),
//  fPhi(obj.fPhi),
//  fEta(obj.fEta)
//{
//  // copy constructor
//}
//_____________________________________________________________________________
AliFemtoProtonParticle &AliFemtoProtonParticle::operator=(const AliFemtoProtonParticle &obj)
{
 //Assignment operator
 if(this == &obj) return *this;

 fMomentum = obj.fMomentum;
 fMomentumMC = obj.fMomentumMC;
 fMomentumMCMother = obj.fMomentumMCMother;
 fMomentumMCMotherParton = obj.fMomentumMCMotherParton;
 fPDGCode = obj.fPDGCode;
 fPDGCodeMother = obj.fPDGCodeMother;
 fPDGCodePartonMother = obj.fPDGCodePartonMother;
 fPartonMotherLabel = obj.fPartonMotherLabel;
 fReal = obj.fReal;
 fProtonTag = obj.fProtonTag;

 fPt = obj.fPt;
 fID = obj.fID;
 fPhi = obj.fPhi;
 fEta = obj.fEta;

 for(int i=0;i<9;i++)//nine different TPC radii, ok its not good to hard code the number
   {
     fPhistar[i] = obj.fPhistar[i];
     fPositionTPC[i] = obj.fPositionTPC[i];
     /*
     for(int j=0;j<3;j++)
       {
	 fPrimPosTPC[i][j] = obj.fPrimPosTPC[i][j];
       }
     */
   }

 return (*this);
}
//_____________________________________________________________________________
AliFemtoProtonParticle::~AliFemtoProtonParticle()
{
  // Destructor
}
//_____________________________________________________________________________
