#include "AliVTrack.h"
#include "AliVVertex.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliDalitzData.h"
#include "AliDalitzAODESD.h"
#include "AliAODHeader.h"
ClassImp( AliDalitzData )
//-----------------------------------------------------------------------------------------------
    AliDalitzData::AliDalitzData():
    fESDEvent(0),
    fAODEvent(0),
    fIsESD(kTRUE)
        {
    
        }
    
       AliDalitzData::AliDalitzData(AliESDEvent* lESDEvent):
    fESDEvent(0),
    fAODEvent(0),
    fIsESD(kTRUE)
       {
        fESDEvent=lESDEvent;
    };
        AliDalitzData::AliDalitzData(AliAODEvent* lAODEvent):
    fESDEvent(0),
    fAODEvent(0),
    fIsESD(kFALSE)
        {
        fAODEvent=lAODEvent;
    };
    AliDalitzData::~AliDalitzData(){}

      AliDalitzAODESD* AliDalitzData::GetTrack(Int_t i){
         if (fIsESD==kTRUE){
                AliDalitzAODESD* esdtrack = new AliDalitzAODESD((AliESDtrack*)fESDEvent->GetTrack(i));
                esdtrack->ComputeImpactParameter();
                return esdtrack;}
         else { AliDalitzAODESD* aodtrack= new AliDalitzAODESD((AliAODTrack*)fAODEvent->GetTrack(i));
                aodtrack->ComputeImpactParameter(fAODEvent->GetPrimaryVertex(),fAODEvent->GetMagneticField());
                return aodtrack;}
     }

    Int_t AliDalitzData::GetNumberOfTracks(){
         if (fIsESD==kTRUE) return fESDEvent->GetNumberOfTracks();
         else return fAODEvent->GetNumberOfTracks();
    }
    Int_t AliDalitzData::GetNumberOfTrackletsG(){
         if (fIsESD==kTRUE) return fESDEvent->GetMultiplicity()->GetNumberOfTracklets();
         else return fAODEvent->GetTracklets()->GetNumberOfTracklets();
    }
    Int_t AliDalitzData::GetNumberOfITSClustersG(Int_t i){
         if (fIsESD==kTRUE) return fESDEvent->GetMultiplicity()->GetNumberOfITSClusters(i);
         else return ((AliAODHeader*)fAODEvent->GetHeader())->GetNumberOfITSClusters(i);
    }
    const AliVVertex* AliDalitzData::GetPrimaryVertex(){
        if (fIsESD==kTRUE) return fESDEvent->GetPrimaryVertex();
        else return fAODEvent->GetPrimaryVertex();
    }
     Int_t AliDalitzData::GetNumberOfV0s(){
         if (fIsESD==kTRUE) return fESDEvent->GetNumberOfV0s();
         else return fAODEvent->GetNumberOfV0s();
     }
     
     
