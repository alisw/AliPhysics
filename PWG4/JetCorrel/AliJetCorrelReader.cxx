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

//__________________________________________________________________________
// Class for input (ESD or AOD) reading.
// At the moment only ESD input is really implemented, AOD to be added later.
// Its products are the Trigger&Associated particle lists
//-- Author: Paul Constantin

#include "AliJetCorrelReader.h"

using namespace std;

ClassImp(AliJetCorrelReader)

AliJetCorrelReader::AliJetCorrelReader() :
  fjcESD(NULL), fSelector(NULL), fWriter(NULL){
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

Float_t AliJetCorrelReader::GetMultiplicity() const {
  // event multiplicity
  if(!fjcESD){
    std::cerr<<"AliJetCorrelReader::GetVertex() - ERROR : fjcESD not set!"<<std::endl; 
    exit(-1);
  }
  // return fjcESD->GetNumberOfTracks(); // ESD no of global tracks
  const AliMultiplicity* m = fjcESD->GetMultiplicity(); // SPD no of tracklets
  return m->GetNumberOfTracklets();
}

Float_t AliJetCorrelReader::GetVertex() const {
  // event vertex
  if(!fjcESD){
    std::cerr<<"AliJetCorrelReader::GetVertex() - ERROR : fjcESD not set!"<<std::endl; 
    exit(-1);
  }
  return fjcESD->GetPrimaryVertex()->GetZ();
}

Bool_t AliJetCorrelReader::VtxOutPipe() const {
  // returns true if vertex R >= beam pipe
  Float_t xVtx2 = fjcESD->GetPrimaryVertex()->GetX()*fjcESD->GetPrimaryVertex()->GetX();
  Float_t yVtx2 = fjcESD->GetPrimaryVertex()->GetY()*fjcESD->GetPrimaryVertex()->GetY();
  if(TMath::Sqrt(xVtx2+yVtx2)>3) return kTRUE;
  return kFALSE;
}

void AliJetCorrelReader::FillLists(CorrelList_t *trigList, CorrelList_t *assoList){
  // fills the trigger&associated particle lists, in this order!
  cPartType_t trigType = trigList->PartID();
  cPartType_t assoType = assoList->PartID();
  Bool_t filledTrig = trigList->Filled();
  Bool_t filledAsso = assoList->Filled();
  // when needed, fill both lists simultaneously for speed:
  if(!filledTrig && !filledAsso && trigType==assoType){
    switch(trigType){
    case t_hadron:
      FillESDTrackLists(trigList,assoList); // trigger list first!
      break;
    default:
      std::cerr<<"AliJetCorrelReader::FillLists() - ERROR: type not implemented!"<<std::endl;
      exit(-1);
    }
    trigList->SetFilled(kTRUE);
    assoList->SetFilled(kTRUE);
  } else {
    if(!filledTrig){
      FillList(trigList,kTRUE);
      trigList->SetFilled(kTRUE);
    }
    if(!filledAsso){
      FillList(assoList,kFALSE);
      assoList->SetFilled(kTRUE);
    }
  }
}

void AliJetCorrelReader::FillList(CorrelList_t *list, Bool_t isTrigg){
  // calls the appropriate Fill method for the list's particle type
  cPartType_t partType = list->PartID();
  switch(partType){
  case t_hadron:
    FillESDTrackList(list,isTrigg);
    break;
  case t_electron: 
    FillESDTrackList(list,isTrigg);
    break;
  case t_dielectron:
    FillESDDielectronList(list,isTrigg);
    break;
  case t_photon:
    FillESDPhotonList(list,isTrigg);
    break;
  case t_diphoton:
    FillESDDiphotonList(list,isTrigg);
    break;
  default:
    std::cerr<<"AliJetCorrelReader::FillList() - ERROR: type not implemented!"<<std::endl;
    exit(-1);
  }
}

void AliJetCorrelReader::FillESDTrackLists(CorrelList_t *trigList, CorrelList_t *assoList){
  // fills trigg&assoc lists simultaneously with ESD tracks
  cPartType_t partType = trigList->PartID(); // by definition the two lists store same particle

  UInt_t nTracks = fjcESD->GetNumberOfTracks() ;
  if(nTracks<1) return;
  for(register UInt_t i=0; i<nTracks; i++){
    AliESDtrack *track = (AliESDtrack*)fjcESD->GetTrack(i);

    Float_t pT = track->Pt();
    if(pT<fSelector->MinLowBin(t_asso)) continue;
    if(fSelector->GenQA()) fWriter->FillTrackQA(track,0);
    if(fSelector->LowQualityTrack(track)) continue;
    if(!fSelector->PassPID(track,partType)) continue;
    if(fSelector->GenQA()) fWriter->FillTrackQA(track,1);

    // fill CorrelList_t object with CorrelTrack_t if two-track cuts (TPC entrance) are used and
    // with CorrelParticle_t if not (base particle uses much less memory). Pair (ghost) cuts need
    // to be applied for both real/mixed pairs, but it's not clear yet whether this is necessary.
//     CorrelTrack_t *hadr = new CorrelTrack_t; // if uncommenting these 3 lines, comment out next line
//     const AliExternalTrackParam* tpcEntry = track->GetInnerParam();
//     hadr->SetTPCEntry(tpcEntry->GetX(),tpcEntry->GetY(),tpcEntry->GetZ());
    CorrelParticle_t *hadr = new CorrelParticle_t;
    hadr->SetPt(pT*track->Charge());
    hadr->SetPhi(track->Phi());
    hadr->SetEta(track->Eta());
    hadr->SetMass(track->GetMass());
    hadr->SetID(partType);

    if(fSelector->GetBin(t_trig,pT)>=0) trigList->Push(hadr->Copy());
    if(fSelector->GetBin(t_asso,pT)>=0) assoList->Push(hadr->Copy());
    delete hadr;
  } // ESD track loop
}

void AliJetCorrelReader::FillESDTrackList(CorrelList_t *list, Bool_t isTrigg){
  // this method is called for: (1) associated hadrons, when trigger is not hadron;
  // (2) electrons to be used in dielectron reconstruction. Assoc pT cuts apply then...
  if(isTrigg)
    {std::cerr<<"AliJetCorrelReader::FillESDTrackList - ERROR: not meant for triggers!"<<std::endl; exit(-1);}
  cPartType_t partType = list->PartID();

  UInt_t nTracks = fjcESD->GetNumberOfTracks();
  if(nTracks<1) return;
  for(register UInt_t i=0; i<nTracks; i++){
    AliESDtrack *track = (AliESDtrack*)fjcESD->GetTrack(i);

    Float_t pT = track->Pt();
    if(pT<fSelector->MinLowBin(t_asso)) continue; 
    if(fSelector->GenQA()) fWriter->FillTrackQA(track,0);
    if(fSelector->LowQualityTrack(track)) continue;
    if(fSelector->GenQA()) fWriter->FillTrackQA(track,1);
    if(!fSelector->PassPID(track,partType)) continue;
    
    // fill CorrelList_t object with CorrelKFTrack_t if AliKFParticle is used for di-electrons and 
    // with CorrelParticle_t if this is done via TLorentzVector in CorrelRecoParent_t (less memory).
    // AliKFParticle allows for a vertex cut on the reconstructed di-electron
    if(partType==t_electron && fSelector->UseAliKF()){
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

void AliJetCorrelReader::FillESDPhotonList(CorrelList_t *list, Bool_t isTrigg){
  // TBI
  std::cerr<<"WARNING : FillESDPhotonList() not emplemented yet. Doing nothing..."<<std::endl;
  std::cerr<<"Photon list size:"<<list->Size()<<std::endl;
}

void AliJetCorrelReader::FillESDDiphotonList(CorrelList_t* list, Bool_t isTrigg){
  // makes a diphoton list (see above warning!)
  CorrelList_t *fPhotonList = new CorrelList_t;
  fPhotonList->Label(t_photon,0); // event number unimportant here
  FillESDPhotonList(fPhotonList, isTrigg);
  FillParentList(list, fPhotonList, isTrigg);
  fWriter->FillParentNtuple(list);
  delete fPhotonList;
}

void AliJetCorrelReader::FillESDDielectronList(CorrelList_t* list, Bool_t isTrigg){
  // makes a dielectron list
  CorrelList_t *fElectronList = new CorrelList_t;
  fElectronList->Label(t_electron,0); // event number unimportant here
  FillESDTrackList(fElectronList, isTrigg);
  FillParentList(list, fElectronList, isTrigg);
  fWriter->FillParentNtuple(list);
  delete fElectronList;
}

void AliJetCorrelReader::FillParentList(CorrelList_t *ParentList, CorrelList_t *ChildList, Bool_t isTrigg){
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
      parent->SetEvent(fjcESD);
      Bool_t goodParent = parent->Reconstruct(child1, child2, fSelector->UseAliKF());
      Bool_t inPtRange = (!isTrigg && fSelector->GetBin(t_asso,parent->Pt())>=0) ||
	(isTrigg && fSelector->GetBin(t_trig,parent->Pt())>=0);
      if(goodParent && inPtRange) ParentList->Push(parent);
    } // 2nd particle loop
  } // 1st particle loop
}
