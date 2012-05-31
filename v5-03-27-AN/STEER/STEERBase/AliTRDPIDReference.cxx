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
// Container class for the reference distributions for TRD PID
// The class contains the reference distributions and the momentum steps
// the references are taken at. Mapping is done inside. To derive references,
// the functions GetUpperReference and GetLowerReference return the next
// reference distribution object and the momentum step above respectively below
// the tracklet momentum.
//
// Authors:
//    Markus Fasel <M.Fasel@gsi.de>
//
#include "TObjArray.h"

#include "AliLog.h"

#include "AliTRDPIDReference.h"

ClassImp(AliTRDPIDReference)

//____________________________________________________________
AliTRDPIDReference::AliTRDPIDReference():
TNamed(),
fRefContainer(NULL),
fMomentumBins()
{
	//
	// Dummy constructor
	//
	SetBit(kIsOwner, kTRUE);
}

//____________________________________________________________
AliTRDPIDReference::AliTRDPIDReference(const Char_t *name):
		  TNamed(name, "TRD PID References"),
		  fRefContainer(NULL),
		  fMomentumBins()
{
	//
	// Default constructor
	//
	SetBit(kIsOwner, kTRUE);
}

//____________________________________________________________
AliTRDPIDReference::AliTRDPIDReference(const AliTRDPIDReference &ref):
		  TNamed(ref),
		  fRefContainer(ref.fRefContainer),
		  fMomentumBins(ref.fMomentumBins)
{
	//
	// Copy constructor
	// Only copies poiters, object is not the owner of the references
	//
	SetBit(kIsOwner, kFALSE);
}

//____________________________________________________________
AliTRDPIDReference &AliTRDPIDReference::operator=(const AliTRDPIDReference &ref){
	//
	// Assginment operator
	// Only copies poiters, object is not the owner of the references
	//
	if(this != &ref){
		TNamed::operator=(ref);
		if(TestBit(kIsOwner) && fRefContainer) delete fRefContainer;
		fRefContainer = ref.fRefContainer;
		fMomentumBins = ref.fMomentumBins;
		SetBit(kIsOwner, kFALSE);
	}
	return *this;
}

//____________________________________________________________
AliTRDPIDReference::~AliTRDPIDReference(){
	//
	// Destructor
	// references are deleted if the object is the owner
	//
	if(fRefContainer && TestBit(kIsOwner)) delete fRefContainer;
}

//____________________________________________________________
void AliTRDPIDReference::SetNumberOfMomentumBins(Int_t nBins, Float_t *momenta){
	//
	// Set the momentum binning
	//
	if(fRefContainer) fRefContainer->Clear();
	else{
		fRefContainer = new TObjArray;
		fRefContainer->SetOwner();
	}
	fRefContainer->Expand(nBins * AliPID::kSPECIES);
	fMomentumBins.Set(nBins,momenta);
}

//____________________________________________________________
void AliTRDPIDReference::AddReference(TObject *ref, AliPID::EParticleType spec, Int_t pbin){
	//
	// Add a new reference distribution for a given species and a givem
	// momentum value to the reference container
	// The reference distribution is associated with the momentum value defined for the
	// given momentum bin
	//
	if(!fRefContainer){
		AliError("Reference Container not initialized");
		return;
	}
	if(pbin > fMomentumBins.GetSize()){
		AliError("Pbin overflow");
		return;
	}
	AliDebug(1, Form("Adding object with address %p to position %d", ref, spec * fMomentumBins.GetSize() + pbin));
	fRefContainer->AddAt(ref, spec * fMomentumBins.GetSize() + pbin);
}

//____________________________________________________________
TObject *AliTRDPIDReference::GetUpperReference(AliPID::EParticleType spec, Float_t p, Float_t &pUpper) const{
	//
	// Get the next reference associated with a momentum larger than the requested momentum
	// In case no next upper reference is found, NULL is returned
	// The momentum value the reference is associated to is stored in the reference pUpper
	//
	Int_t bin = -1;
	pUpper = 20;
	for(Int_t ip = 0; ip < fMomentumBins.GetSize(); ip++){
		AliDebug(10, Form("Bin %d, p = %.1f", ip, fMomentumBins[ip]));
		if(p < fMomentumBins[ip]){
			AliDebug(10, "Bin found");
			bin = ip;
			break;
		}
	}
	AliDebug(2, Form("p = %.1f, bin = %d\n", p, bin));
	if(bin >= 0) {
		pUpper = fMomentumBins[bin];
		return fRefContainer->At(spec * fMomentumBins.GetSize() + bin);
	}
	else return NULL;
}

//____________________________________________________________
TObject *AliTRDPIDReference::GetLowerReference(AliPID::EParticleType spec, Float_t p, Float_t &pLower) const{
	//
	// Get the next reference associated with a momentum smaller than the requested momentum
	// In case no next lower reference is found, NULL is returned
	// The momentum value the reference is associated to is stored in the reference pLower
	//
	Int_t bin = -1;
	pLower = 0;
	for(Int_t ip = fMomentumBins.GetSize() - 1; ip >= 0; ip--){
		AliDebug(10, Form("Bin %d, p = %.1f", ip, fMomentumBins[ip]));
		if(p > fMomentumBins[ip]){
			AliDebug(10, "Bin found");
			bin = ip;
			break;
		}
	}
	AliDebug(2, Form("p = %.1f, bin = %d\n", p, bin));
	if(bin >= 0){
		pLower = fMomentumBins[bin];
		return fRefContainer->At(spec * fMomentumBins.GetSize() + bin);
	}
	else return NULL;
}

//____________________________________________________________
void AliTRDPIDReference::Print(const Option_t*) const{
	//
	// Print content of the PID reference container
	//
	printf("Number of Momentum Bins: %d\n", GetNumberOfMomentumBins());
	printf("=====================================\n");
	for(Int_t ip = 0; ip < GetNumberOfMomentumBins(); ip++){
		printf("Bin %d: p = %.1f\n", ip, fMomentumBins[ip]);
	}
	printf("=====================================\n");
	if(fRefContainer){
		printf("Content of the reference container:\n");
		for(Int_t ip = 0; ip < GetNumberOfMomentumBins(); ip++){
			printf("[");
			for(Int_t is = 0; is < 5; is++) printf("%p|", fRefContainer->At(is * fMomentumBins.GetSize() + ip));
			printf("]\n");
		}
	}
}
