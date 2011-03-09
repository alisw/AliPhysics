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

// -----------------------------------------------
// particle level cuts for azimuthal isotropic
// expansion in highly central collisions analysis
// author: Cristian Andrei
//         acristian@niham.nipne.ro
// ------------------------------------------------

#include <TF1.h>
#include <TFile.h>

#include "AliAnalysisCentralCutESD.h"
#include "AliESDtrack.h"


class TObject;
class TParticle;


//____________________________________________________________________
ClassImp(AliAnalysisCentralCutESD)

//____________________________________________________________________
AliAnalysisCentralCutESD::AliAnalysisCentralCutESD(const Char_t* name, const Char_t* title) 
    :AliAnalysisCuts(name,title)
    ,fReqPID(kFALSE)
    ,fReqCharge(kFALSE)
    ,fPartType(kPiPlus)
    ,fPIDtype("Custom")
    ,fPriorsFunc(kFALSE)
    ,fPartPriors()
    ,fElectronFunction(0)
    ,fMuonFunction(0)
    ,fPionFunction(0)
    ,fKaonFunction(0)
    ,fProtonFunction(0) 
{
// Constructor
// Initialize the priors
    fPartPriors[0] = 0.01;
    fPartPriors[1] = 0.01;
    fPartPriors[2] = 0.85;
    fPartPriors[3] = 0.1;
    fPartPriors[4] = 0.05;

    if(fPriorsFunc){
	TFile *f = TFile::Open("$ALICE_ROOT/PWG2/data/PriorProbabilities.root ");
	if(!f){
	    printf("Can't open PWG2 prior probabilities file!\n Exiting ...\n");
	    return;
	}
	fElectronFunction = (TF1 *)f->Get("fitElectrons");
	fMuonFunction = (TF1 *)f->Get("fitMuons");
	fPionFunction = (TF1 *)f->Get("fitPions");
	fKaonFunction = (TF1 *)f->Get("fitKaons");
	fProtonFunction = (TF1 *)f->Get("fitProtons");
    }

}

AliAnalysisCentralCutESD::~AliAnalysisCentralCutESD() {
// Destructor
// Delete the created priors


	if(fElectronFunction) delete fElectronFunction;
	if(fMuonFunction) delete fMuonFunction;
	if(fPionFunction) delete fPionFunction;
	if(fKaonFunction) delete fKaonFunction;
	if(fProtonFunction) delete fProtonFunction;

}


Bool_t AliAnalysisCentralCutESD::IsSelected(TObject *obj){
// Checks if a particle passes the cuts

    AliESDtrack *track = dynamic_cast<AliESDtrack *>(obj);

    if(!track){
		printf("AliAnalysisCentralCutESD:IsSelected ->Can't get track!\n");
		return kFALSE;
    }

    if(fReqCharge){
		if(!IsCharged(track)) return kFALSE;
    } 

    if(fReqPID){
		if(!IsA(track, fPartType)) return kFALSE;
    }

    return kTRUE;
}


Double_t AliAnalysisCentralCutESD::GetPriors(Int_t i, Double_t p) {
//Return the a priori probs

Double_t priors=0;
    if(fPriorsFunc) {
	if(i == 0) priors = fElectronFunction->Eval(p);
	if(i == 1) priors = fMuonFunction->Eval(p);
	if(i == 2) priors = fPionFunction->Eval(p);
	if(i == 3) priors = fKaonFunction->Eval(p);
	if(i == 4) priors = fProtonFunction->Eval(p);
    }
    else {
	priors = fPartPriors[i];
    }

  return priors;
}



Bool_t AliAnalysisCentralCutESD::IsA(AliESDtrack *track, PDG_t reqPartType){
// Determines the type of the particle
    Int_t charge;

    if(reqPartType < 0){
      charge = -1;
    }
    else{
      charge = 1;
    }

    Double_t probability[5] = {0.0,0.0,0.0,0.0,0.0};
    Double_t w[5] = {0.0,0.0,0.0,0.0,0.0};

    Long64_t partType = 0;

    Double_t p = track->P();

    track->GetESDpid(probability);

    Double_t s = 0.0;


    for(Int_t i = 0; i < AliPID::kSPECIES; i++){ 
		s += probability[i]*GetPriors(i,p);
    }
    
    if(!s < 0.000001) {
		for(Int_t i = 0; i < AliPID::kSPECIES; i++){ 
	    	w[i] = probability[i]*GetPriors(i,p)/s;
		}
    }


    if(fPIDtype.Contains("Bayesian")) {
		partType = TMath::LocMax(AliPID::kSPECIES,w);
		if(partType<0.) return kFALSE;
    }

    else if(fPIDtype.Contains("Custom")){
		for(Int_t i=0;i<AliPID::kSPECIES;i++)   {
			if(w[i]>0.9){
				partType = i;
			}
		}
    }

    else{
		printf("Unknown PID method!\n");
		return kFALSE;
    }
    if(partType<0.) return kFALSE;
    if((AliPID::ParticleCode(partType)) != reqPartType){
		return kFALSE;
    }

	if(track->Charge() != charge) return kFALSE;

    return kTRUE;

}

Bool_t AliAnalysisCentralCutESD::IsCharged(AliESDtrack* const track) const{

    if(track->Charge() == 0) return kFALSE;

    return kTRUE;

}
