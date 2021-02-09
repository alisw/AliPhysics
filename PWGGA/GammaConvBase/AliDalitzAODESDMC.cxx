#include "AliVTrack.h"
#include "AliDalitzAODESDMC.h"
#include "TParticle.h"
#include "AliAODMCParticle.h"
#include "AliLog.h"
#include "TMath.h"

ClassImp( AliDalitzAODESDMC )
//-----------------------------------------------------------------------------------------------

AliDalitzAODESDMC::~AliDalitzAODESDMC() 
{
}

    Int_t AliDalitzAODESDMC::GetMotherG(){
        if (fIsESDMC==kTRUE) return tparticle->GetMother(0);
        else return aliaodmcparticle->GetMother();
    }
    Int_t AliDalitzAODESDMC::GetFirstDaughterG(){
        if (fIsESDMC==kTRUE) return tparticle->GetFirstDaughter();
        else return aliaodmcparticle->GetDaughterFirst();
    }
    Int_t AliDalitzAODESDMC::GetLastDaughterG(){
        if (fIsESDMC==kTRUE) return tparticle->GetLastDaughter();
        else return aliaodmcparticle->GetDaughterLast();
    }
    Int_t AliDalitzAODESDMC::GetNDaughtersG(){
        if (fIsESDMC==kTRUE) return tparticle->GetNDaughters();
        else return aliaodmcparticle->GetNDaughters();
    }
    //Int_t AliDalitzAODESDMC::GetFirstG(){
    //    if (fIsESD==kTRUE) return tparticle->GetMother(0);
    //    else return aliaodmcparaticle->GetMother();
    //}
    Int_t AliDalitzAODESDMC::GetPdgCodeG(){
        if (fIsESDMC==kTRUE) return tparticle->GetPdgCode();
        else return aliaodmcparticle->GetPdgCode();
    }
    Double_t AliDalitzAODESDMC::GetRatioVxyG(){
        if (fIsESDMC==kTRUE) return tparticle->R();
        else {
        Double_t vv[3];
        aliaodmcparticle->XvYvZv(vv);
        return TMath::Sqrt(vv[0]*vv[0]+vv[1]*vv[1]);
        }
    }
    Double_t AliDalitzAODESDMC::PtG(){
        if (fIsESDMC==kTRUE) return tparticle->Pt();
        else return aliaodmcparticle->Pt();
    }
    Double_t AliDalitzAODESDMC::PxG(){
        if (fIsESDMC==kTRUE) return tparticle->Px();
        else return aliaodmcparticle->Px();
    }
    Double_t AliDalitzAODESDMC::PyG(){
        if (fIsESDMC==kTRUE) return tparticle->Py();
        else return aliaodmcparticle->Py();
    }
    Double_t AliDalitzAODESDMC::PzG(){
        if (fIsESDMC==kTRUE) return tparticle->Pz();
        else return aliaodmcparticle->Pz();
    }
    Double_t AliDalitzAODESDMC::EtaG(){
        if (fIsESDMC==kTRUE) return tparticle->Eta();
        else return aliaodmcparticle->Eta();
    }
     Int_t AliDalitzAODESDMC::GetUniqueIDG(){
        if (fIsESDMC==kTRUE) return tparticle->GetUniqueID();
        else return aliaodmcparticle->GetMCProcessCode();
    }
    Double_t AliDalitzAODESDMC::VertexOnZ(){
        if (fIsESDMC==kTRUE) return tparticle->Vz();
        else return aliaodmcparticle->Xv();
    }
    Double_t AliDalitzAODESDMC::EnergyG(){
        if (fIsESDMC==kTRUE) return tparticle->Energy();
        else return aliaodmcparticle->E();
    }
    
