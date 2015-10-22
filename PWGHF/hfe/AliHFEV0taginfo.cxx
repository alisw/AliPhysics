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
//
// Class AliHFEV0taginfo
// Creates list of tracks with V0 information
//
// Author:
//   Jan Wagner <J.Wagner@gsi.de>
//




#include <iostream>
#include <TClass.h>
#include <TList.h>
#include <TMath.h>


#include "AliLog.h"
#include "AliAODEvent.h"
#include "AliAODv0.h"
#include "AliESDEvent.h"
#include "AliESDv0.h"

#include "AliHFEV0taginfo.h"


ClassImp(AliHFEV0taginfo)
ClassImp(AliHFEV0taginfo::AliHFEV0tag)

//___________________________________________________________________
AliHFEV0taginfo::AliHFEV0taginfo():
    TNamed(), 
    fIsAODana(false),
    fTaggedTracks(NULL),
    fV0finder(NULL),
    fAODV0finder(NULL)
{
    //
    // default constructor
    //
}
//___________________________________________________________________
AliHFEV0taginfo::AliHFEV0taginfo(const char* name):
    TNamed(name, ""), 
    fIsAODana(kFALSE),
    fTaggedTracks(NULL),
    fV0finder(NULL),
    fAODV0finder(NULL)
{
    //
    // constructor
    //

    fTaggedTracks = new TList();
    if(fTaggedTracks){
        fTaggedTracks->SetOwner();
    }
    fV0finder = new AliESDv0KineCuts();
    fAODV0finder = new AliAODv0KineCuts();
}

//__________________________________________________________________
AliHFEV0taginfo::AliHFEV0taginfo(const AliHFEV0taginfo &ref):
  TNamed(ref),
    fIsAODana(ref.fIsAODana),
    fTaggedTracks(NULL),
    fV0finder(ref.fV0finder),
    fAODV0finder(ref.fAODV0finder)
{
  //
  // Copy constructor
  // creates a new object with new (empty) containers
  //
    fTaggedTracks = new TList();
    if(fTaggedTracks){
        fTaggedTracks->SetOwner();
    }


    AliHFEV0tag *tmp = NULL;
    for(Int_t ien = 0; ien < ref.fTaggedTracks->GetEntries(); ien++){
        tmp = static_cast<AliHFEV0tag *>(ref.fTaggedTracks->At(ien));
        fTaggedTracks->Add(new AliHFEV0tag(tmp->GetTrackID(),tmp->GetPinfo(),tmp->GetProdR()));
    }
}

//__________________________________________________________________
AliHFEV0taginfo &AliHFEV0taginfo::operator=(const AliHFEV0taginfo &ref){
    //
    // Assignment operator
    // Cleanup old object, create a new one with new containers inside
    //
    if(this == &ref) return *this;
    this->~AliHFEV0taginfo(); 
    TNamed::operator=(ref);
    fIsAODana = ref.fIsAODana;
    fTaggedTracks = new TList();
    AliHFEV0tag *tmp = NULL;
    for(Int_t ien = 0; ien < ref.fTaggedTracks->GetEntries(); ien++){
        tmp = static_cast<AliHFEV0tag *>(ref.fTaggedTracks->At(ien));
        fTaggedTracks->Add(new AliHFEV0tag(tmp->GetTrackID(),tmp->GetPinfo(),tmp->GetProdR()));
    }
    fV0finder=ref.fV0finder;
    fAODV0finder=ref.fAODV0finder;
    return *this;
}
//___________________________________________________________________
AliHFEV0taginfo::~AliHFEV0taginfo(){

    //
    // Destructor
    //
    delete fTaggedTracks;
    delete fV0finder;
    delete fAODV0finder;
    AliDebug(6, "DESTRUCTOR");
}

//________________________________________________________________________________
// loops over V0s in event and fills fTaggedTracks with V0 tracks
void AliHFEV0taginfo::TagV0Tracks(AliVEvent *fEvent){

    if (!fEvent) return;

    const Int_t nTracks = fEvent->GetNumberOfTracks();
    if(nTracks < 2) return;

    const Int_t nV0s = fEvent->GetNumberOfV0s();
    if(nV0s < 1) return;
    AliDebug(3,Form("%d V0s found!",nV0s));
  
    if(fEvent->IsA() == AliESDEvent::Class()){
        AliDebug(4, "ESD part");
        AliESDEvent *esdevent = static_cast<AliESDEvent *>(fEvent);
        fV0finder->SetEvent(esdevent);

        for(Int_t i=0; i<nV0s; ++i){
            Int_t pdgP = 0;
            Int_t pdgN = 0;
            AliESDv0 *fV0 = esdevent->GetV0(i);
            if(!fV0) continue;
            if(fV0finder->ProcessV0(fV0,pdgP,pdgN)){
                AliDebug(5,Form("V0 has: pos pdg: %d, neg pdg: %d",pdgP,pdgN));
                AddTrack(fV0->GetPindex(),pdgP,TMath::Sqrt(fV0->Xv()*fV0->Xv()+fV0->Yv()*fV0->Yv()));
                AddTrack(fV0->GetNindex(),pdgN,TMath::Sqrt(fV0->Xv()*fV0->Xv()+fV0->Yv()*fV0->Yv()));
            }
        }
    } else if(fEvent->IsA() == AliAODEvent::Class()){
        AliDebug(4,"AOD part");
        AliAODEvent *aodevent = static_cast<AliAODEvent *>(fEvent);
        fAODV0finder->SetEvent(aodevent);

        for(Int_t i=0; i<nV0s; ++i){
            Int_t pdgP = 0;
            Int_t pdgN = 0;
            AliAODv0 *fV0 = aodevent->GetV0(i);
            if(!fV0) continue;
            if(fAODV0finder->ProcessV0(fV0,pdgP,pdgN)){
                AliDebug(5,Form("V0 has: pos pdg: %d, neg pdg: %d",pdgP,pdgN));
                AddTrack(fV0->GetPosID(),pdgP,fV0->RadiusV0());
                AddTrack(fV0->GetNegID(),pdgN,fV0->RadiusV0());
            }
        }
    }
}

//________________________________________________________________________________
//Translates the pdg code to AliPID enum and adds track to tagged list
void AliHFEV0taginfo::AddTrack(Int_t TrackID, Int_t pdgCode, Double_t prodR){

    if(TrackID<0) return;
    AliPID::EParticleType Pinfo;
    switch (TMath::Abs(pdgCode)){
        case  11:
            Pinfo = AliPID::kElectron;
            break;
        case  211:
            Pinfo = AliPID::kPion;
            break;
        case  2212:
            Pinfo = AliPID::kProton;
            break;
        default:
            return;
    }
    fTaggedTracks->Add(new AliHFEV0tag(TrackID, Pinfo, prodR));
    AliDebug(4,Form("Added new Track ID: %d with PID: %d, #entry: %d",TrackID, Pinfo, fTaggedTracks->GetEntries()));
}


//________________________________________________________________________________
//check for V0 information from track ID 
//returns AliPID::kUnknown if track ID not found
AliPID::EParticleType AliHFEV0taginfo::GetV0Info(Int_t trackID){

    AliHFEV0tag test(trackID, AliPID::kUnknown,0);
    AliHFEV0tag *result = dynamic_cast<AliHFEV0tag *>(fTaggedTracks->FindObject(&test));
    if(!result){ 
        AliDebug(6, Form("Could not find track ID %d", trackID));
        return AliPID::kUnknown;
    }
    return result->GetPinfo();
}

//________________________________________________________________________________
//check for V0 daughter production vertex from track ID
//returns -0.1 if track ID not found
Float_t AliHFEV0taginfo::GetV0ProdR(Int_t trackID){

    AliHFEV0tag test(trackID, AliPID::kUnknown, 0);
    AliHFEV0tag *result = dynamic_cast<AliHFEV0tag *>(fTaggedTracks->FindObject(&test));
    if(!result){
        AliDebug(6, Form("Could not find track ID %d", trackID));
        return -0.1;
    }
    return result->GetProdR();
}

//________________________________________________________________________________
//resets the fTaggedTracks TList
void AliHFEV0taginfo::Reset(){
    fTaggedTracks->Delete();
}


//___________________________________________________________________
AliHFEV0taginfo::AliHFEV0tag::AliHFEV0tag():
    TObject(), 
    fTrackID(0),
    fPinfo(AliPID::kUnknown),
    fProdR(0)
{
    // default constructor
}
//___________________________________________________________________
AliHFEV0taginfo::AliHFEV0tag::AliHFEV0tag(Int_t TrackID, AliPID::EParticleType Pinfo, Double_t ProdR):
    TObject(), 
    fTrackID(TrackID),
    fPinfo(Pinfo),
    fProdR(ProdR)
{
}

//____________________________________________________________
AliHFEV0taginfo::AliHFEV0tag::AliHFEV0tag(const AliHFEV0tag &ref):
    TObject(ref),
    fTrackID(ref.fTrackID),
    fPinfo(ref.fPinfo),
    fProdR(ref.fProdR)
{
    // Copy constructor
}

//____________________________________________________________
AliHFEV0taginfo::AliHFEV0tag &AliHFEV0taginfo::AliHFEV0tag::operator=(const AliHFEV0tag &ref){
    // Assignment operator
    if(this != &ref){
        TObject::operator=(ref);

        fTrackID = ref.fTrackID;
        fPinfo = ref.fPinfo;
        fProdR = ref.fProdR;
    }
    return *this;
}

//___________________________________________________________________
AliHFEV0taginfo::AliHFEV0tag::~AliHFEV0tag(){
    //
    // Destructor
    //
    AliDebug(6, "DESTRUCTOR");
}

//Set track ID and particle info 
//___________________________________________________________________
void AliHFEV0taginfo::AliHFEV0tag::SetTrack(Int_t trackID, AliPID::EParticleType Pinfo){
    fTrackID = trackID;
    fPinfo = Pinfo;
}

//Set track ID and production verxtex
//___________________________________________________________________
void AliHFEV0taginfo::AliHFEV0tag::SetProdR(Int_t trackID, Double_t prodR){
    fTrackID = trackID;
    fProdR = prodR;
}

//____________________________________________________________
Bool_t AliHFEV0taginfo::AliHFEV0tag::IsEqual(const TObject *ref) const {
    //
    // Check for equality  of track ID
    //
    const AliHFEV0tag *refObj = dynamic_cast<const AliHFEV0tag *>(ref);
    if(!refObj) return kFALSE;
    return (fTrackID == refObj->GetTrackID());
}
//____________________________________________________________
Int_t AliHFEV0taginfo::AliHFEV0tag::Compare(const TObject *ref) const{
    //
    // Compares two objects
    // Order:
    //   First compare track ID then particle info
    //
    const AliHFEV0tag *refObj = static_cast<const AliHFEV0tag *>(ref);
    if(fTrackID < refObj->GetTrackID()) return -1;
    else if(fTrackID > refObj->GetTrackID()) return 1;
    else{
        if(fPinfo < refObj->GetPinfo()) return -1;
        else if(fPinfo > refObj->GetPinfo()) return 1;
        else return 0;
    }
}


