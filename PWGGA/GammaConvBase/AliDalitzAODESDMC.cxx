#include "AliVTrack.h"
#include "AliDalitzAODESDMC.h"
// #include "TParticle.h"
#include "AliAODMCParticle.h"
#include "AliLog.h"
#include "TMath.h"

ClassImp( AliDalitzAODESDMC )
//-----------------------------------------------------------------------------------------------

AliDalitzAODESDMC::~AliDalitzAODESDMC() 
{
}

    Int_t AliDalitzAODESDMC::GetMotherG(){
        if (fIsESDMC==kTRUE) return alimcparticle->GetMother();
        else return aliaodmcparticle->GetMother();
    }
    Int_t AliDalitzAODESDMC::GetFirstDaughterG(){
        if (fIsESDMC==kTRUE) return alimcparticle->GetDaughterFirst();
        else return aliaodmcparticle->GetDaughterFirst();
    }
    Int_t AliDalitzAODESDMC::GetLastDaughterG(){
        if (fIsESDMC==kTRUE) return alimcparticle->GetDaughterLast();
        else return aliaodmcparticle->GetDaughterLast();
    }
    Int_t AliDalitzAODESDMC::GetNDaughtersG(){
        if (fIsESDMC==kTRUE) return alimcparticle->GetNDaughters();
        else return aliaodmcparticle->GetNDaughters();
    }
    //Int_t AliDalitzAODESDMC::GetFirstG(){
    //    if (fIsESD==kTRUE) return alimcparticle->GetMother();
    //    else return aliaodmcparaticle->GetMother();
    //}
    Int_t AliDalitzAODESDMC::GetPdgCodeG(){
        if (fIsESDMC==kTRUE) return alimcparticle->PdgCode();
        else return aliaodmcparticle->GetPdgCode();
    }
    Double_t AliDalitzAODESDMC::GetRatioVxyG(){
        if (fIsESDMC==kTRUE) return TMath::Sqrt(TMath::Power(alimcparticle->Xv(),2)+TMath::Power(alimcparticle->Yv(),2));
        else {
            Double_t vv[3];
            aliaodmcparticle->XvYvZv(vv);
            return TMath::Sqrt(vv[0]*vv[0]+vv[1]*vv[1]);
        }
    }
    Double_t AliDalitzAODESDMC::PtG(){
        if (fIsESDMC==kTRUE) return alimcparticle->Pt();
        else return aliaodmcparticle->Pt();
    }
    Double_t AliDalitzAODESDMC::PxG(){
        if (fIsESDMC==kTRUE) return alimcparticle->Px();
        else return aliaodmcparticle->Px();
    }
    Double_t AliDalitzAODESDMC::PyG(){
        if (fIsESDMC==kTRUE) return alimcparticle->Py();
        else return aliaodmcparticle->Py();
    }
    Double_t AliDalitzAODESDMC::PzG(){
        if (fIsESDMC==kTRUE) return alimcparticle->Pz();
        else return aliaodmcparticle->Pz();
    }
    Double_t AliDalitzAODESDMC::EtaG(){
        if (fIsESDMC==kTRUE) return alimcparticle->Eta();
        else return aliaodmcparticle->Eta();
    }
     Int_t AliDalitzAODESDMC::GetUniqueIDG(){
        if (fIsESDMC==kTRUE) return alimcparticle->GetUniqueID();
        else return aliaodmcparticle->GetMCProcessCode();
    }
    Double_t AliDalitzAODESDMC::VertexOnZ(){
        if (fIsESDMC==kTRUE) return alimcparticle->Zv();
        else return aliaodmcparticle->Zv();
    }
    Double_t AliDalitzAODESDMC::EnergyG(){
        if (fIsESDMC==kTRUE) return alimcparticle->E();
        else return aliaodmcparticle->E();
    }
    
