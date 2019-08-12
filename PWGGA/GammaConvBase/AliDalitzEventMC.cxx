#include "AliVTrack.h"
#include "AliVVertex.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliDalitzEventMC.h"
#include "AliDalitzAODESDMC.h"
ClassImp( AliDalitzEventMC )
//-----------------------------------------------------------------------------------------------
  AliDalitzEventMC::AliDalitzEventMC():
    fESDEventMC(0),
    fAODEvent(0),
    fIsESDMC(kTRUE)
        {
    
        }
    
       AliDalitzEventMC::AliDalitzEventMC(AliMCEvent* lESDMCEvent):
    fESDEventMC(0),
    fAODEvent(0),
    fIsESDMC(kTRUE)
       {
         fESDEventMC=lESDMCEvent;
    };
        AliDalitzEventMC::AliDalitzEventMC(AliAODEvent* lAODMCEvent):
     fESDEventMC(0),
    fAODEvent(0),
    fIsESDMC(kFALSE)
        {
            fAODEvent=lAODMCEvent;
    };
    
        AliDalitzEventMC::~AliDalitzEventMC(){
        }
        
        

     AliDalitzAODESDMC* AliDalitzEventMC::Particle(Int_t i){
         if (fIsESDMC==kTRUE){
                AliDalitzAODESDMC* esdparticle= new AliDalitzAODESDMC((TParticle*)fESDEventMC->Particle(i));
                return esdparticle;}
                
        else {
            
              TClonesArray *AODMCTrackArray = dynamic_cast<TClonesArray*>(fAODEvent->FindListObject(AliAODMCParticle::StdBranchName()));
              
              AliAODMCParticle *aodparticle0 = (AliAODMCParticle*) AODMCTrackArray->At(i);
              AliDalitzAODESDMC* aodparticle1= new AliDalitzAODESDMC(aodparticle0);
              return aodparticle1;
        }
     }
     
