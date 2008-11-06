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
///////////////////////////////////////////////////////////////////////////
// The class AliCFManager is designed to handle inside the 
// task the basic components for a particle-level correction of single 
// particle analysis
// The class provides methods to set lists of cuts and to loop over them 
// for several different selection steps to be then used for
// efficiency calculation.
// prototype version by S.Arcelli silvia.arcelli@cern.ch
///////////////////////////////////////////////////////////////////////////
#include "AliCFCutBase.h"
#include "AliCFManager.h"

ClassImp(AliCFManager)

//_____________________________________________________________________________
AliCFManager::AliCFManager() : 
  TNamed(),
  fNStepEvt(0),
  fNStepPart(0),
  fEvtContainer(0x0),
  fPartContainer(0x0),
  fEvtCutList(0x0),
  fPartCutList(0x0)
{ 
  //
  // ctor
  //
}
//_____________________________________________________________________________
AliCFManager::AliCFManager(Char_t* name, Char_t* title) : 
  TNamed(name,title),
  fNStepEvt(0),
  fNStepPart(0),
  fEvtContainer(0x0),
  fPartContainer(0x0),
  fEvtCutList(0x0),
  fPartCutList(0x0)
{ 
   //
   // ctor
   //
}
//_____________________________________________________________________________
AliCFManager::AliCFManager(const AliCFManager& c) : 
  TNamed(c),
  fNStepEvt(c.fNStepEvt),
  fNStepPart(c.fNStepPart),
  fEvtContainer(c.fEvtContainer),
  fPartContainer(c.fPartContainer),
  fEvtCutList(c.fEvtCutList),
  fPartCutList(c.fPartCutList)
{ 
   //
   //copy ctor
   //
}
//_____________________________________________________________________________
AliCFManager& AliCFManager::operator=(const AliCFManager& c)
{
  //
  // Assignment operator
  //
  if (this != &c) {
    TNamed::operator=(c) ;
  }
  
  this->fNStepEvt=c.fNStepEvt;
  this->fNStepPart=c.fNStepPart;
  this->fEvtContainer=c.fEvtContainer;
  this->fPartContainer=c.fPartContainer;
  this->fEvtCutList=c.fEvtCutList;
  this->fPartCutList=c.fPartCutList;
  return *this ;
}

//_____________________________________________________________________________
AliCFManager::~AliCFManager() {
   //
   //dtor
   //
}

//_____________________________________________________________________________
Bool_t AliCFManager::CheckParticleCuts(Int_t isel, TObject *obj, const TString  &selcuts) const {
  //
  // check whether object obj passes particle-level selection isel
  //

  if(isel>=fNStepPart){
    AliWarning(Form("Selection index out of Range! isel=%i, max. number of selections= %i", isel,fNStepPart));
    return kTRUE;
  }
  if(!fPartCutList[isel])return kTRUE;
  TObjArrayIter iter(fPartCutList[isel]);
  AliCFCutBase *cut = 0;
  while ( (cut = (AliCFCutBase*)iter.Next()) ) {
    TString cutName=cut->GetName();
    Bool_t checkCut=CompareStrings(cutName,selcuts);
    if(checkCut && !cut->IsSelected(obj)) return kFALSE;   
  }
  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliCFManager::CheckEventCuts(Int_t isel, TObject *obj, const TString  &selcuts) const{
  //
  // check whether object obj passes event-level selection isel
  //

  if(isel>=fNStepEvt){
    AliWarning(Form("Selection index out of Range! isel=%i, max. number of selections= %i", isel,fNStepEvt));
      return kTRUE;
  }
  if(!fEvtCutList[isel])return kTRUE;
  TObjArrayIter iter(fEvtCutList[isel]);
  AliCFCutBase *cut = 0;
  while ( (cut = (AliCFCutBase*)iter.Next()) ) {
    TString cutName=cut->GetName();
    Bool_t checkCut=CompareStrings(cutName,selcuts);
    if(checkCut && !cut->IsSelected(obj)) return kFALSE;   
  }
  return kTRUE;
}

//_____________________________________________________________________________
void  AliCFManager::SetEventInfo(TObject *obj) const {

  //Particle level cuts

  if (!fPartCutList) {
    AliWarning("No particle cut list found");
  }
  else {
    for(Int_t isel=0;isel<fNStepPart; isel++){
      if(!fPartCutList[isel])continue;  
      TObjArrayIter iter(fPartCutList[isel]);
      AliCFCutBase *cut = 0;
      while ( (cut = (AliCFCutBase*)iter.Next()) ) {
	cut->SetEvtInfo(obj);
      }    
    }
  }
  
  //Event level cuts 
  
  if (!fEvtCutList) {
    AliWarning("No event cut list found");
  }
  else {
    for(Int_t isel=0;isel<fNStepEvt; isel++){
      if(!fEvtCutList[isel])continue;  
      TObjArrayIter iter(fEvtCutList[isel]);
      AliCFCutBase *cut = 0;
      while ( (cut = (AliCFCutBase*)iter.Next()) ) {
	cut->SetEvtInfo(obj);
      }   
    }
  }
}

//_____________________________________________________________________________
Bool_t AliCFManager::CompareStrings(const TString  &cutname,const TString  &selcuts) const{
  //
  // compare two strings
  //

  if(selcuts.Contains("all"))return kTRUE;
  if ( selcuts.CompareTo(cutname) == 0 ||
       selcuts.BeginsWith(cutname+" ") ||
       selcuts.EndsWith(" "+cutname) ||
       selcuts.Contains(" "+cutname+" "))  return kTRUE; 
  return kFALSE;
}


//_____________________________________________________________________________
void AliCFManager::SetEventCutsList(Int_t isel, TObjArray* array) {
  //
  //Setter for event-level selection cut list at selection step isel
  //

  if (!fEvtContainer) {
    AliWarning("No event container defined, you may need to set it first!"); 
  }

  Int_t nstep = fNStepEvt;

  if (!fEvtCutList) fEvtCutList = new TObjArray*[nstep] ;
  if (isel >= nstep) {
    AliWarning(Form("Selection index out of Range! isel=%i, max. number of selections= %i", isel,nstep));
    return;
  }
  fEvtCutList[isel] = array;
}

//_____________________________________________________________________________
void AliCFManager::SetParticleCutsList(Int_t isel, TObjArray* array) {
  //
  //Setter for particle-level selection cut list at selection step isel
  //

  if (!fPartContainer) {
    AliWarning("No particle container defined, you may need to set it first!"); 
  }
  
  Int_t nstep = fNStepPart ;
  
  if (!fPartCutList) fPartCutList = new TObjArray*[nstep] ;
  if (isel >= nstep) {
    AliWarning(Form("Selection index out of Range! isel=%i, max. number of selections= %i", isel,nstep));
    return;
  }
  fPartCutList[isel] = array;
}
