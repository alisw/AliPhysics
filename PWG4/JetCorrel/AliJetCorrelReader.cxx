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
/* $Id: $ */

//______________________________________________________________________________________
// Class for input (ESD or AOD) reading and filling of Trigger&Associated particle lists
//-- Author: Paul Constantin

#include "AliJetCorrelReader.h"

using namespace std;

ClassImp(AliJetCorrelReader)

AliJetCorrelReader::AliJetCorrelReader() :
  jcESD(NULL), fSelector(NULL), fWriter(NULL){
  // constructor
}

AliJetCorrelReader::~AliJetCorrelReader(){
  // destructor
}

void AliJetCorrelReader::Init(AliJetCorrelSelector * const s, AliJetCorrelWriter * const w){
  // initialization method
  fSelector = s;
  fWriter = w;
}

Float_t AliJetCorrelReader::GetMultiplicity(){
  // event multiplicity
  if(!jcESD){
    std::cerr<<"AliJetCorrelReader::GetVertex() - ERROR : jcESD not set!"<<std::endl; 
    exit(-1);
  }
  // return jcESD->GetNumberOfTracks(); // ESD no of global tracks
  const AliMultiplicity* m = jcESD->GetMultiplicity(); // SPD no of tracklets
  return m->GetNumberOfTracklets();
}

Float_t AliJetCorrelReader::GetVertex(){
  // event vertex
  if(!jcESD){
    std::cerr<<"AliJetCorrelReader::GetVertex() - ERROR : jcESD not set!"<<std::endl; 
    exit(-1);
  }
  return jcESD->GetPrimaryVertex()->GetZ();
}

Bool_t AliJetCorrelReader::VtxOutPipe(){
  // returns true if vertex R >= beam pipe
  Float_t xVtx2 = jcESD->GetPrimaryVertex()->GetX()*jcESD->GetPrimaryVertex()->GetX();
  Float_t yVtx2 = jcESD->GetPrimaryVertex()->GetY()*jcESD->GetPrimaryVertex()->GetY();
  if(TMath::Sqrt(xVtx2+yVtx2)>3) return kTRUE;
  return kFALSE;
}

void AliJetCorrelReader::FillLists(CorrelList_t *list1, CorrelList_t *list2){
  // fills the trigger&associated particle lists
  PartType_t partType1 = list1->PartID();
  PartType_t partType2 = list2->PartID();
  Bool_t filled1 = list1->Filled();
  Bool_t filled2 = list2->Filled();
  // when needed, fill both lists simultaneously for speed:
  if(!filled1 && !filled2 && partType1==partType2){
    switch(partType1){
    case hadron:
      FillESDTrackLists(list1,list2);
      break;
    default:
      std::cerr<<"AliJetCorrelReader::FillLists() - ERROR: type not implemented!"<<std::endl;
      exit(-1);
    }
    list1->SetFilled(kTRUE);
    list2->SetFilled(kTRUE);
  } else {
    if(!filled1){
      FillList(list1);
      list1->SetFilled(kTRUE);
    }
    if(!filled2){
      FillList(list2);
      list2->SetFilled(kTRUE);
    }
  }
}

void AliJetCorrelReader::FillList(CorrelList_t *list){
  // calls the appropriate Fill method for the list's particle type
  PartType_t partType = list->PartID();
  switch(partType){
  case hadron:
    FillESDTrackList(list);
    break;
  case electron: 
    FillESDTrackList(list);
    break;
  case dielectron:
    FillESDDielectronList(list);
    break;
  case photon:
    FillESDPhotonList(list);
    break;
  case diphoton:
    FillESDDiphotonList(list);
    break;
  default:
    std::cerr<<"AliJetCorrelReader::FillList() - ERROR: type not implemented!"<<std::endl;
    exit(-1);
  }
}

void AliJetCorrelReader::FillESDTrackLists(CorrelList_t *list1, CorrelList_t *list2){
  // fills trigg&assoc lists simultaneously with ESD tracks
  PartType_t partType = list1->PartID(); // by definition the two lists store same particle

  UInt_t nTracks = jcESD->GetNumberOfTracks() ;
  if(nTracks<1) return;
  for(register UInt_t i=0; i<nTracks; i++){
    AliESDtrack *track = (AliESDtrack*)jcESD->GetTrack(i);

    Float_t pT = track->Pt();
    if(pT<fSelector->MinAssocPt()) continue;
    if(fSelector->GenQA()) fWriter->FillTrackQA(track,0);
    if(fSelector->LowQualityTrack(track)) continue;
    if(!fSelector->PassPID(track,partType)) continue;
    if(fSelector->GenQA()) fWriter->FillTrackQA(track,1);

    // fill CorrelList_t object with CorrelTrack_t if two-track cuts (TPC entrance) are used and
    // with CorrelParticle_t if not (base particle uses much less memory). Pair (ghost) cuts need
    // to be applied for both real/mixed pairs, but it's not clear yet whether this is necessary.
    CorrelTrack_t *hadr = new CorrelTrack_t;
    const AliExternalTrackParam* tpcEntry = track->GetInnerParam();
    hadr->SetTPCEntry(tpcEntry->GetX(),tpcEntry->GetY(),tpcEntry->GetZ());
//     CorrelParticle_t *hadr = new CorrelParticle_t; // comment out above 4 lines first
    hadr->SetPt(pT*track->Charge());
    hadr->SetPhi(track->Phi());
    hadr->SetEta(track->Eta());
    hadr->SetMass(track->GetMass());
    hadr->SetID(partType);

    if(list1->PoolID()==assocs && fSelector->AssocBin(pT)>=0) list1->Push(hadr->Copy());
    if(list1->PoolID()==triggs && fSelector->TriggBin(pT)>=0) list1->Push(hadr->Copy());
    if(list2->PoolID()==assocs && fSelector->AssocBin(pT)>=0) list2->Push(hadr->Copy());
    if(list2->PoolID()==triggs && fSelector->TriggBin(pT)>=0) list2->Push(hadr->Copy());
    delete hadr;
  } // ESD track loop
}

void AliJetCorrelReader::FillESDTrackList(CorrelList_t *list){
  // this method is called for: (1) associated hadrons, when trigger is not hadron;
  // (2) electrons to be used in dielectron reconstruction. Assoc pT cuts apply then...
  PartType_t partType = list->PartID();

  UInt_t nTracks = jcESD->GetNumberOfTracks();
  if(nTracks<1) return;
  for(register UInt_t i=0; i<nTracks; i++){
    AliESDtrack *track = (AliESDtrack*)jcESD->GetTrack(i);

    Float_t pT = track->Pt();
    if(pT<fSelector->MinAssocPt()) continue; 
    if(fSelector->GenQA()) fWriter->FillTrackQA(track,0);
    if(fSelector->LowQualityTrack(track)) continue;
    if(fSelector->GenQA()) fWriter->FillTrackQA(track,1);
    if(!fSelector->PassPID(track,partType)) continue;
    
    // fill CorrelList_t object with CorrelKFTrack_t if AliKFParticle is used for di-electrons and 
    // with CorrelParticle_t if this is done via TLorentzVector in CorrelRecoParent_t (less memory).
    // AliKFParticle allows for a vertex cut on the reconstructed di-electron
    if(partType==electron && fSelector->UseAliKF()){
      CorrelKFTrack_t *elec = new CorrelKFTrack_t;
      const AliExternalTrackParam* tPar = track->GetConstrainedParam();
      elec->SetParam(tPar->GetParameter());
      elec->SetCovar(tPar->GetCovariance());
      elec->SetPt(pT*track->Charge());
      elec->SetPhi(track->Phi());
      elec->SetEta(track->Eta());
      elec->SetMass(track->GetMass());
      elec->SetID(partType);
      list->Push(elec);
    } else {
      CorrelParticle_t *hadr = new CorrelParticle_t;
      hadr->SetPt(pT*track->Charge());
      hadr->SetPhi(track->Phi());
      hadr->SetEta(track->Eta());
      hadr->SetMass(track->GetMass());
      hadr->SetID(partType);
      list->Push(hadr);
    }   
  }
}

void AliJetCorrelReader::FillESDPhotonList(CorrelList_t *list){
  // TBI
  std::cerr<<"WARNING : FillESDPhotonList() not emplemented yet. Doing nothing..."<<std::endl;
  std::cerr<<"Photon list size:"<<list->Size()<<std::endl;
}

void AliJetCorrelReader::FillESDDiphotonList(CorrelList_t* list){
  // makes a diphoton list (see above warning!)
  CorrelList_t *fPhotonList = new CorrelList_t;
  fPhotonList->Label(photon,assocs,0); // event number unimportant here
  FillESDPhotonList(fPhotonList);
  FillParentList(list, fPhotonList);
  fWriter->FillParentNtuple(list);
  delete fPhotonList;
}

void AliJetCorrelReader::FillESDDielectronList(CorrelList_t* list){
  // makes a dielectron list
  CorrelList_t *fElectronList = new CorrelList_t;
  fElectronList->Label(electron,assocs,0); // event number unimportant here
  FillESDTrackList(fElectronList);
  FillParentList(list, fElectronList);
  fWriter->FillParentNtuple(list);
  delete fElectronList;
}

void AliJetCorrelReader::FillParentList(CorrelList_t *ParentList, CorrelList_t *ChildList){
  // makes a list of parent particles from a list of children of same type
  if(ChildList->Size()<2) return;

  CorrelListIter_t iterChild1, iterChild2;
  iterChild1 = ChildList->Head();
  while(!iterChild1.HasEnded()){
    CorrelParticle_t *child1 = iterChild1.Data(); iterChild1.Move();
    iterChild2 = iterChild1;
    while(!iterChild2.HasEnded()){
      CorrelParticle_t *child2 = iterChild2.Data(); iterChild2.Move();
      CorrelRecoParent_t *parent = new CorrelRecoParent_t;
      parent->SetEvent(jcESD);
      Bool_t goodParent = parent->Reconstruct(child1, child2, fSelector->UseAliKF());
      Bool_t inPtRange = (ParentList->PoolID()==assocs && fSelector->AssocBin(parent->Pt())>=0) ||
	(ParentList->PoolID()==triggs && fSelector->TriggBin(parent->Pt())>=0);
      if(goodParent && inPtRange) ParentList->Push(parent);
    } // 2nd particle loop
  } // 1st particle loop
}
