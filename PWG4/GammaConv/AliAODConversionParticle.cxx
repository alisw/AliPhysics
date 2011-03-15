/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: Ana Marin, Kathrin Koch, Kenneth Aamodt, Svein Lindal          *
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

#include "AliAODConversionParticle.h"
#include "AliStack.h"
#include "AliESDEvent.h"
#include "TParticle.h"
#include <iostream>

class AliAODv0;
class AliESDtrack;

using namespace std;

ClassImp(AliAODConversionParticle)



AliAODConversionParticle::AliAODConversionParticle() :
  AliAODPhoton(),
  fChi2(0),
  fIMass(0), 
  fTagged(kFALSE),
  fMCStack(NULL),
  fESDEvent(NULL)
{
  //Default constructor
  fLabel[0] = -1;
  fLabel[1] = -1;
}

AliAODConversionParticle::AliAODConversionParticle(TLorentzVector & momentum) :
  AliAODPhoton(momentum),
  fChi2(-1),
  fIMass(-1), 
  fTagged(kFALSE),
  fMCStack(NULL),
  fESDEvent(NULL)
{
  //Default constructor
  fLabel[0] = -1;
  fLabel[1] = -1;
}

AliAODConversionParticle::AliAODConversionParticle(AliGammaConversionAODObject *obj):
  AliAODPhoton(obj->Px(),obj->Py(),obj->Pz(),obj->E()),
  fChi2(0),
  fIMass(0), 
  fTagged(kFALSE),
  fMCStack(NULL),
  fESDEvent(NULL)
{

    fIMass=obj->IMass();
 
 //  fIMass=obj->M();
    //Default constructor
  fLabel[0] = obj->GetLabel1();
  fLabel[1] = obj->GetLabel2();

}

AliAODConversionParticle::AliAODConversionParticle(AliKFParticle * gammakf,Int_t label1,Int_t label2):
  AliAODPhoton(gammakf->Px(),gammakf->Py(),gammakf->Pz(), gammakf->E()),
  fChi2(0),
  fIMass(0), 
  fTagged(kFALSE),
  fMCStack(NULL),
  fESDEvent(NULL)
{
  
     fLabel[0] = label1;
     fLabel[1] = label2;

     fIMass=gammakf->GetMass();

    
}



AliAODConversionParticle::AliAODConversionParticle(AliAODConversionParticle *y1,AliAODConversionParticle *y2):
  AliAODPhoton(y1->Px()+y2->Px(),y1->Py()+y2->Py(),y1->Pz()+y2->Pz(),y1->E()+y2->E()),
  fChi2(0),
  fIMass(0), 
  fTagged(kFALSE),
  fMCStack(NULL),
  fESDEvent(NULL)
{
    fIMass=M();
    fLabel[0] = -1;
    fLabel[1] = -1;
}

AliAODConversionParticle::AliAODConversionParticle(const AliAODConversionParticle & original) :
  AliAODPhoton(original),
  fChi2(original.fChi2),
  fIMass(original.fIMass),
  fTagged(original.fTagged),
  fMCStack(original.fMCStack),
  fESDEvent(original.fESDEvent)
{
  //Copy constructor
  fLabel[0] = original.fLabel[0];
  fLabel[1] = original.fLabel[1];
}


AliAODConversionParticle & AliAODConversionParticle::operator = (const AliAODConversionParticle & /*source*/)
{
  // assignment operator
  return *this;
}


Int_t AliAODConversionParticle::GetGammaMCLabel() const{
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

Int_t AliAODConversionParticle::GetElectronUniqueID(Int_t mcLabel1, Int_t mcLabel2) const{
  Int_t iResult = -1;
  if(fMCStack != NULL){
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

Int_t AliAODConversionParticle::GetElectronUniqueID() const{
  // returns the unique id of the electrons if they have the same mother and unique id
  if(fMCStack != NULL){
    return GetElectronUniqueID(GetElectronMCLabel1(), GetElectronMCLabel2());
  }  
  return -1;
}


Int_t AliAODConversionParticle::GetMCLabel(Int_t label) const {
  //returns the MC label of the first electron
  Int_t iResult=-1;
  if(fESDEvent != NULL){
    if(label>=0){
      iResult = fESDEvent->GetTrack(label)->GetLabel();
    }
  }
  return iResult;
}

Int_t AliAODConversionParticle::GetElectronMCLabel1() const{
  //returns the MC label of the first electron
  if(fLabel[0] >=0) {
    return GetMCLabel(fLabel[0]);
  }
  
  return -1;
}

Int_t AliAODConversionParticle::GetElectronMCLabel2() const{
  //returns the MC label of the first electron
  if(fLabel[1] >=0) {
    return GetMCLabel(fLabel[1]);
  }
  
  return -1;
}

///_________________________________________________________
Bool_t AliAODConversionParticle::IsMySpawn(const Int_t trackId, const Int_t nSpawn, const Int_t * const spawnList) const {
  for(int iSpawn = 0; iSpawn < nSpawn; iSpawn++) {
    //cout << spawnList[iSpawn] << endl;;
    if(spawnList[iSpawn] == trackId) return kTRUE;
  }
  return kFALSE;
}

