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
//    Daniel Lohner <Daniel.Lohner@cern.ch>

#include "AliLog.h"

#include "AliTRDPIDResponseObject.h"

#ifndef AliTRDPIDREFERENCE_H
#include "AliTRDPIDReference.h"
#endif

#ifndef AliTRDPIDPARAMS_H
#include "AliTRDPIDParams.h"
#endif


ClassImp(AliTRDPIDResponseObject)

//____________________________________________________________
AliTRDPIDResponseObject::AliTRDPIDResponseObject():
    TNamed(),
    fNSlicesQ0(4)
{
    //
    // Dummy constructor
    //
    SetBit(kIsOwner, kTRUE);

    for(Int_t method=0;method<AliTRDPIDResponse::kNMethod;method++){
	fkPIDParams[method]=NULL;
	fkPIDReference[method]=NULL;
    }
}

//____________________________________________________________
AliTRDPIDResponseObject::AliTRDPIDResponseObject(const Char_t *name):
TNamed(name, "TRD PID Response Object"),
fNSlicesQ0(4)
{
	//
	// Default constructor
	//
	SetBit(kIsOwner, kTRUE);

	for(Int_t method=0;method<AliTRDPIDResponse::kNMethod;method++){
	    fkPIDParams[method]=NULL;
	    fkPIDReference[method]=NULL;
	}
}

//____________________________________________________________
AliTRDPIDResponseObject::AliTRDPIDResponseObject(const AliTRDPIDResponseObject &ref):
TNamed(ref),
fNSlicesQ0(ref.fNSlicesQ0)
{
    //
    // Copy constructor
    // Only copies pointers, object is not the owner of the references
    //
    SetBit(kIsOwner, kFALSE);

    for(Int_t method=0;method<AliTRDPIDResponse::kNMethod;method++){
	fkPIDParams[method]=ref.fkPIDParams[method];       // new Object is not owner, copy only pointer
	fkPIDReference[method]=ref.fkPIDReference[method];    // new Object is not owner, copy only pointer
    }
}
//____________________________________________________________
AliTRDPIDResponseObject &AliTRDPIDResponseObject::operator=(const AliTRDPIDResponseObject &ref){
	//
	// Assginment operator
	// Only copies poiters, object is not the owner of the references
	//
	if(this != &ref){
	    TNamed::operator=(ref);
	    fNSlicesQ0=ref.fNSlicesQ0;
	    for(Int_t method=0;method<AliTRDPIDResponse::kNMethod;method++){
		if(TestBit(kIsOwner) && fkPIDParams[method])delete fkPIDParams[method];
                if(TestBit(kIsOwner) && fkPIDReference[method])delete fkPIDReference[method];

		fkPIDParams[method]=ref.fkPIDParams[method];       // new Object is not owner, copy only pointer
		fkPIDReference[method]=ref.fkPIDReference[method];    // new Object is not owner, copy only pointer
	    }
	    SetBit(kIsOwner, kFALSE);
	}
	return *this;
}

//____________________________________________________________
AliTRDPIDResponseObject::~AliTRDPIDResponseObject(){
	//
	// Destructor
	// references are deleted if the object is the owner
	//
    if(fkPIDParams && TestBit(kIsOwner)) delete fkPIDParams;
    if(fkPIDReference && TestBit(kIsOwner)) delete fkPIDReference;

}

//____________________________________________________________
void AliTRDPIDResponseObject::SetPIDParams(AliTRDPIDParams *params,AliTRDPIDResponse::ETRDPIDMethod method){

    if(Int_t(method)>=Int_t(AliTRDPIDResponse::kNMethod)||Int_t(method)<0){
	AliError("Method does not exist");
	return;
    }
    if(fkPIDParams[method]){
	delete fkPIDParams[method];
        fkPIDParams[method]=NULL;
    }

    fkPIDParams[method]=new AliTRDPIDParams(*params);
}

//____________________________________________________________
void AliTRDPIDResponseObject::SetPIDReference(AliTRDPIDReference *reference,AliTRDPIDResponse::ETRDPIDMethod method){

    if(Int_t(method)>=Int_t(AliTRDPIDResponse::kNMethod)||Int_t(method)<0){
        AliError("Method does not exist");
	return;
    }
    if(fkPIDReference[method]){
	delete fkPIDReference[method];
	fkPIDReference[method]=NULL;
    }
    fkPIDReference[method]=new AliTRDPIDReference(*reference);
}

//____________________________________________________________
TObject *AliTRDPIDResponseObject::GetUpperReference(AliPID::EParticleType spec, Float_t p, Float_t &pUpper,AliTRDPIDResponse::ETRDPIDMethod method) const{

    if(Int_t(method)>=Int_t(AliTRDPIDResponse::kNMethod)||Int_t(method)<0){
	AliError("Method does not exist");
	return NULL;
    }
   
    if(fkPIDReference[method]){
	return fkPIDReference[method]->GetUpperReference(spec,p,pUpper);
    }
    return NULL;
}


//____________________________________________________________
TObject *AliTRDPIDResponseObject::GetLowerReference(AliPID::EParticleType spec, Float_t p, Float_t &pLower,AliTRDPIDResponse::ETRDPIDMethod method) const{

    if(Int_t(method)>=Int_t(AliTRDPIDResponse::kNMethod)||Int_t(method)<0){
	AliError("Method does not exist");
	return NULL;
    }

     if(fkPIDReference[method]){
	 return fkPIDReference[method]->GetLowerReference(spec,p,pLower);
     }
    return NULL;
}

//____________________________________________________________
Bool_t AliTRDPIDResponseObject::GetThresholdParameters(Int_t ntracklets, Double_t efficiency, Double_t *params,AliTRDPIDResponse::ETRDPIDMethod method) const{

    if(Int_t(method)>=Int_t(AliTRDPIDResponse::kNMethod)||Int_t(method)<0){
	AliError("Method does not exist");
	return kFALSE;
    }

    if(fkPIDParams[method]){
	return fkPIDParams[method]->GetThresholdParameters(ntracklets,efficiency,params);
    }
    return kFALSE;
}

//____________________________________________________________
Int_t AliTRDPIDResponseObject::GetNumberOfMomentumBins(AliTRDPIDResponse::ETRDPIDMethod method) const{

    if(Int_t(method)>=Int_t(AliTRDPIDResponse::kNMethod)||Int_t(method)<0){
	AliError("Method does not exist");
	return 0;
    }

    if(fkPIDReference[method]){
	return fkPIDReference[method]->GetNumberOfMomentumBins();
    }
    return 0;
}

//____________________________________________________________
void AliTRDPIDResponseObject::Print(const Option_t* opt) const{
	//
	// Print content of the PID object
	//
    printf("Content of AliTRDPIDResponseObject \n\n");
   
    for(Int_t method=0;method<AliTRDPIDResponse::kNMethod;method++){
	if(fkPIDReference[method])fkPIDReference[method]->Print(opt);
	if(fkPIDParams[method])printf("+ Threshold Parameters \n");
    }
}
