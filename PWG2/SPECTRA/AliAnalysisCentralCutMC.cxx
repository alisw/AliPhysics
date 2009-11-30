/*************************************************************************
* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
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

//  ----------------------------------------------------
//  MC particle level cuts for azimuthal isotropic
//  expansion in highly central collisions analysis 
//  author: Cristian Andrei
//          acristian@niham.nipne.ro
//  ----------------------------------------------------

#include <TParticle.h>

#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliMCParticle.h"

#include "AliAnalysisCentralCutMC.h"

class TObject;

//____________________________________________________________________
ClassImp(AliAnalysisCentralCutMC)

//____________________________________________________________________
AliAnalysisCentralCutMC::AliAnalysisCentralCutMC(const Char_t* name, const Char_t* title) 
    :AliAnalysisCuts(name,title)
    ,fOnlyPrim(0)
    ,fMCEvt(0)
    ,fPDGCode(0)
    ,fPtMin(0)
    ,fPtMax(0)
    ,fEtaMin(0)
    ,fEtaMax(0)
{
    //constructor
    
    fOnlyPrim = kFALSE; //implicit consider si part secundare
    
    SetPDGCode();
    SetPtRange();
    SetEtaRange();

	printf("AliAnalysisCentralCutMC::Constructor\n");
}

//____________________________________________________________________
AliAnalysisCentralCutMC::~AliAnalysisCentralCutMC(){
// Destructor

    if(fMCEvt) delete fMCEvt;

	printf("AliAnalysisCentralCutMC::Destructor\n");
}

//___________________________________________________________________________
Bool_t AliAnalysisCentralCutMC::IsSelected(TObject* const obj){
// Check if the particle passes the cuts

    AliMCParticle *part = dynamic_cast<AliMCParticle *>(obj);

    if(!part){
		printf("AliAnalysisCentralCutMC:IsSelected ->Can't get particle!\n");
		return kFALSE;
    }
    
    if(!fMCEvt){
		printf("AliAnalysisCentralCutMC:IsSelected ->Can't get MCEvent!\n");
		return kFALSE;    
    }
    
    AliStack *stack = fMCEvt->Stack();
		if(!stack){
		printf("AliAnalysisCentralCutMC:IsSelected ->Can't get Stack!\n");
		return kFALSE;
    }

    Double_t pt = part->Pt();

    Double_t eta = part->Eta();

    
    if(fOnlyPrim){
		if(!IsPrimary(part, stack)) return kFALSE;
    }
    
    if(fPDGCode != 0){
		if(!CheckPDG(part, fPDGCode)) return kFALSE;
    }	

    if((pt < fPtMin) || (pt > fPtMax)) return kFALSE;
    
    if((eta < fEtaMin) || (eta > fEtaMax)) return kFALSE;


    return kTRUE;
}

//___________________________________________________________________________
Bool_t AliAnalysisCentralCutMC::IsPrimary(AliMCParticle* const part, AliStack* const stack){
// Check if the particle is primary

    Int_t index = part->GetLabel();
    
    TParticle* p = stack->Particle(index);
    if(!p){
	printf("AliAnalysisCentralCutMC:IsPrimary -> Can't get TParticle!\n");
	return kFALSE;    
    }

    Int_t ist  = p->GetStatusCode();

    if (ist > 1) return kFALSE;

//     Int_t pdg = TMath::Abs(p->GetPdgCode());

    if (index < stack->GetNprimary()) {
    	return kTRUE;
    }
    else return kFALSE;

}


void AliAnalysisCentralCutMC::ReceiveEvt(TObject* mcEvent) {
// Receive the event send from the task
// The event is needed in order to get the Stack

    if (!mcEvent){
        printf("Pointer to MC Event is null! \n");
        return;
    }

    fMCEvt = dynamic_cast<AliMCEvent*> (mcEvent);
    if(!fMCEvt){
		printf("AliAnalysisCentralCutMC:ReceiveEvt -> Can't get fMCEvt!\n");
		return;
    }

}


Bool_t AliAnalysisCentralCutMC::CheckPDG(AliMCParticle* const mcPart, Int_t const pdg) {
// Checks if the particle is of the wanted type

    TParticle* part = mcPart->Particle();

    if(!part){
		printf("AliAnalysisCentralCutMC:IsPrimary -> Can't get TParticle!\n");
		return kFALSE;
    }

    Int_t pdgCode = part->GetPdgCode();

    if (pdgCode != pdg ) return kFALSE;
    
    return kTRUE;
    
}
