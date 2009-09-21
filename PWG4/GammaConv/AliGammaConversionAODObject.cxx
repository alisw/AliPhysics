/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: Ana Marin, Kathrin Koch, Kenneth Aamodt                        *
 * Version 1.0                                                            *
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
//--------------------------------------------- 
// Class containing the aod information we need
//---------------------------------------------
////////////////////////////////////////////////

#include "AliGammaConversionAODObject.h"
//#include "AliAODv0.h"
#include "AliStack.h"
#include "AliESDEvent.h"
//#include "AliESDtrack.h"
#include "TParticle.h"

class AliAODv0;
class AliESDtrack;

using namespace std;

ClassImp(AliGammaConversionAODObject)



AliGammaConversionAODObject::AliGammaConversionAODObject() :
  TObject(),
  fPx(0),
  fPy(0),
  fPz(0),
  fLabel1(-1),
  fLabel2(-1),
  fMCStack(NULL),
  fESDEvent(NULL)
{
	
}


AliGammaConversionAODObject::AliGammaConversionAODObject(const AliGammaConversionAODObject & original) :
  TObject(original),
  fPx(original.fPx),
  fPy(original.fPy),
  fPz(original.fPz),
  fLabel1(original.fLabel1),
  fLabel2(original.fLabel2),
  fMCStack(original.fMCStack),
  fESDEvent(original.fESDEvent)
{
	
}


AliGammaConversionAODObject & AliGammaConversionAODObject::operator = (const AliGammaConversionAODObject & /*source*/)
{
  // assignment operator
  return *this;
}

Int_t AliGammaConversionAODObject::GetGammaMCLabel() const{
  // returns the MC label of the gamma (if both electrons have the same mother)
  Int_t iResult = -1;
  if(fMCStack != NULL){
    Int_t mcLabel1= GetElectronMCLabel1();
    Int_t mcLabel2= GetElectronMCLabel2();
    if(mcLabel1>=0 && mcLabel2>=0){
      TParticle *electron1 = fMCStack->Particle(mcLabel1);
      TParticle *electron2 = fMCStack->Particle(mcLabel2);
      if(electron1->GetMother(0) == electron2->GetMother(0)){
	iResult = electron1->GetMother(0);
      }
    }
  }
  return iResult;
}

Int_t AliGammaConversionAODObject::GetElectronUniqueID() const{
  // returns the unique id of the electrons if they have the same mother and unique id
  Int_t iResult = -1;
  if(fMCStack != NULL){
    Int_t mcLabel1= GetElectronMCLabel1();
    Int_t mcLabel2= GetElectronMCLabel2();
    if(mcLabel1>=0 && mcLabel2>=0){
      TParticle *electron1 = fMCStack->Particle(mcLabel1);
      TParticle *electron2 = fMCStack->Particle(mcLabel2);
      if(electron1->GetMother(0) == electron2->GetMother(0)){
	if(electron1->GetUniqueID() == electron2->GetUniqueID()){
	  iResult = (Int_t)electron1->GetUniqueID();
	}
      }
    }
  }
  return iResult;
}

Int_t AliGammaConversionAODObject::GetElectronMCLabel1() const{
  //returns the MC label of the first electron
  Int_t iResult=-1;
  if(fESDEvent != NULL){
    if(fLabel1>=0){
      iResult = (fESDEvent->GetTrack(fLabel1))->GetLabel();
    }
  }
  return iResult;
}

Int_t AliGammaConversionAODObject::GetElectronMCLabel2() const{
  //returns the MC label of the first electron
  Int_t iResult=-1;
  if(fESDEvent != NULL){
    if(fLabel2>=0){
      iResult = (fESDEvent->GetTrack(fLabel2))->GetLabel();
    }
  }
  return iResult;
}
