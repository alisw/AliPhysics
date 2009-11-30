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

// ---------------------------------------------------
//  MC event level cuts for azimuthal isotropic
//  expansion in highly central collisions analysis
//  author: Cristian ANDREI
//          acristian@niham.nipne.ro
// ----------------------------------------------------


#include <TParticle.h>

#include "AliMCEvent.h"
#include "AliStack.h"

#include "AliAnalysisCentralCutEvtMC.h"


class TObject;

//____________________________________________________________________
ClassImp(AliAnalysisCentralCutEvtMC)

//____________________________________________________________________
AliAnalysisCentralCutEvtMC::AliAnalysisCentralCutEvtMC(const Char_t* name, const Char_t* title) 
    :AliAnalysisCuts(name,title)
    ,fReqMult(kFALSE) 
    ,fReqDir(kFALSE)
    ,fReqDirUnit(kFALSE)
    ,fMultMin(0)
    ,fMultMax(0)
    ,fDirMin(0)
    ,fDirMax(0)
    ,fDirUMin(0)
    ,fDirUMax(0)

{
    //constructor
	printf("AliAnalysisCentralCutEvtMC::Constructor\n");
}

AliAnalysisCentralCutEvtMC::~AliAnalysisCentralCutEvtMC(){
//destructor
printf("AliAnalysisCentralCutEvtMC::Destructor\n");
}

//___________________________________________________________________________
Bool_t AliAnalysisCentralCutEvtMC::IsSelected(TObject *obj){
// Checks whether the event passes the cuts

    if (!obj){
        printf("Pointer to MC Event is null! \n");
        return kFALSE;
    }

    AliMCEvent *mcEvent = dynamic_cast<AliMCEvent*> (obj);
		if(!mcEvent){
		printf("AliAnalysisCentralCutEvtMC -> Can't get MCEvt!\n");
		return kFALSE;
    }


    if(fReqMult){
		Int_t mult = CalcMult(mcEvent);
		if((mult<fMultMin)||(mult>fMultMax)){ 
			return kFALSE;
		}
    }

    if(fReqDir){
		Double_t dir = CalcDir(mcEvent);
		if((dir<fDirMin)||(dir>fDirMax)){
			return kFALSE;
		}
    }

    if(fReqDirUnit){
		Double_t dirU = CalcDirUnit(mcEvent);
		if((dirU<fDirUMin)||(dirU>fDirUMax)){
			return kFALSE;
		}
    }

    return kTRUE;
}

//___________________________________________________________________________
Int_t AliAnalysisCentralCutEvtMC::CalcMult(AliMCEvent* const mcEvent) {
// Computes the multiplicity of the event

    AliStack *stack=mcEvent->Stack();
    if(!stack) return -1;

    Int_t nPrim = stack->GetNprimary();

    Int_t charged = 0;
    Double_t eta;

    for(Int_t i = 0; i < nPrim; i++){//track loop -> compute multiplicity

		TParticle *particle=stack->Particle(i);
		if(!particle){
			continue;
		}
	
		eta = particle->Eta();
	
		if((eta > 0.5)||(eta < -0.5)){         //------------mid-rapidity
			continue;
		}
	
		TParticlePDG *particlePDG = particle->GetPDG(0); 
		if(!particlePDG){
			continue;
		}
		
		if(TMath::Abs(particlePDG->Charge()) < 3){
			continue;
		}
		
		charged++;  

    }//end track loop

    return charged;

}

//____________________________________________________________________________
Double_t AliAnalysisCentralCutEvtMC::CalcDir(AliMCEvent* const mcEvent) {
// computes the directivity of the event
	
    AliStack *stack=mcEvent->Stack();
    if(!stack) return -1;

    Int_t nPrim = stack->GetNprimary();

    Int_t goodtrack = 0;

    Double_t eta;
    Double_t pt;

    Double_t dir;
    Double_t px,py;
    Double_t sumaPt = 0;
    Double_t sumaPx = 0;
    Double_t sumaPy = 0;

    for(Int_t i = 0; i < nPrim; i++){//track loop -> compute directivity

		TParticle *particle=stack->Particle(i);
		if(!particle){
			continue;
		}
		
		eta = particle->Eta();
	
		if((eta > 1.9)||(eta < 0.0)){   // half of the SPD coverage -> directivity 
			continue;
		}
		
		TParticlePDG *particlePDG = particle->GetPDG(0); 
		if(!particlePDG){
			continue;
		}
		
		if(TMath::Abs(particlePDG->Charge()) < 3){
			continue;
		}
	
		px = particle->Px();
		py = particle->Py();
	
		pt = particle->Pt();
	
		sumaPx = sumaPx + px;
		sumaPy = sumaPy + py;
	
		sumaPt = sumaPt + pt;
	
		goodtrack++;
	 
    }//end track loop


    if(sumaPt < 0.0000001){
		return -1;
    }

    dir = (sqrt(pow(sumaPx,2)+pow(sumaPy,2)))/sumaPt;

    return dir;
}

//__________________________________________________________________
Double_t AliAnalysisCentralCutEvtMC::CalcDirUnit(AliMCEvent* const mcEvent) {
// computes the directivity using only the SPD unity vectors

    AliStack *stack=mcEvent->Stack();
    if(!stack){
		return -1;
	}

    Int_t nPrim = stack->GetNprimary();

    Int_t goodtrack = 0;

    Double_t eta;
    Double_t pt;

    Double_t dirU;
    Double_t px,py,pxU,pyU;
    Double_t sumaPxU = 0;
    Double_t sumaPyU = 0;


    for(Int_t i = 0; i < nPrim; i++){//track loop -> compute directivity

		TParticle *particle=stack->Particle(i);
		if(!particle){
			continue;
		}
		
		eta = particle->Eta();
	
		if((eta > 1.9)||(eta < 0.0)){ // half of the SPD coverage -> directivity 
			continue;
		}
		
		TParticlePDG *particlePDG = particle->GetPDG(0); 
		if(!particlePDG){
			continue;
		}
		
		if(TMath::Abs(particlePDG->Charge()) < 3){
			continue;
		}
	
		px = particle->Px();
		py = particle->Py();
		pt = particle->Pt();
	
		pxU = px/pt;
		pyU = py/pt;
	
		sumaPxU = sumaPxU + pxU;
		sumaPyU = sumaPyU + pyU;	
	
		goodtrack++;
	 
    }//end track loop

    if(goodtrack == 0){
		return -1;
	}

    dirU = (sqrt(pow(sumaPxU,2)+pow(sumaPyU,2)))/goodtrack;


    return dirU;
}
