/**************************************************************************
 * Copyright(c) 2007-2010, ALICE Experiment at CERN, All rights reserved. *
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

//========================================================================
//
// This class generates misalignment for TPC. In particular it defines
// the misalignment in the three canonical scenarios: "ideal", "residual"
// and "full".
// It is meant to be run standalone or from the steering macro
// $ALICE_ROOT/macros/MakeAlignmentObjs.C
// looping on the detectors.
//
//========================================================================

#include "AliTPCMisAligner.h"
#include "AliGeomManager.h"
#include "TClonesArray.h"
#include "TGeoManager.h"
#include "TRandom.h"
#include "AliAlignObjParams.h"
#include "AliLog.h"

ClassImp(AliTPCMisAligner)

    //_______________________________________________________________________________________
AliTPCMisAligner::AliTPCMisAligner() : AliMisAligner()
{
    //
    // dummy constructor
    //
}

//_______________________________________________________________________________________
TClonesArray* AliTPCMisAligner::MakeAlObjsArray() {
    // builds and returns the array of alignment objects
    // according to the spcified misalignment scenario
    // ("ideal", "residual" or "full").
    //
    if(!AliGeomManager::GetGeometry())
    {
	AliError("No geometry loaded into AliGeomManager! Returning null pointer!");
	return 0;
    }

    TClonesArray *array = new TClonesArray("AliAlignObjParams",73);
    TClonesArray &alobj = *array;

    gRandom->SetSeed(4357);
    Int_t j = 0;
    // misalignment of the whole TPC according to survey
    Double_t dx=-0.159, dy=-0.05, dz=0.034, dpsi=-0.00183, dtheta=0.01835, dphi=0.02865;
    new(alobj[j++]) AliAlignObjParams("ALIC_1/TPC_M_1", 0, dx, dy, dz, dpsi, dtheta, dphi, kTRUE);
    AliAlignObjParams* alObjTPC = (AliAlignObjParams*) array->UncheckedAt(0);
    alObjTPC->ApplyToGeometry();

    if(TString(GetMisalType())=="ideal")
    {

	dx=0., dy=0., dz=0., dpsi=0., dtheta=0., dphi=0.;
	for (Int_t iLayer = AliGeomManager::kTPC1; iLayer <= AliGeomManager::kTPC2; iLayer++) {
	    for (Int_t iModule = 0; iModule < AliGeomManager::LayerSize(iLayer); iModule++) {

		UShort_t volid = AliGeomManager::LayerToVolUID(iLayer,iModule);
		const char *symname = AliGeomManager::SymName(volid);
		new(alobj[j++]) AliAlignObjParams(symname, volid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE);
	    }
	}

    }else if(TString(GetMisalType())=="residual" || TString(GetMisalType())=="full"){
	// sigma translation = 0.1 mm
	// sigma rotation = 0.1 mrad
	// RS = local
	
	Double_t sigmatr=0.01;
	Double_t sigmarot = 0.006;

	for (Int_t iLayer = AliGeomManager::kTPC1; iLayer <= AliGeomManager::kTPC2; iLayer++) {
	    for (Int_t iModule = 0; iModule < AliGeomManager::LayerSize(iLayer); iModule++) {

		dx = gRandom->Gaus(0,sigmatr);
		dy = gRandom->Gaus(0,sigmatr);
		dz = gRandom->Gaus(0,sigmatr);
		dpsi = gRandom->Gaus(0,sigmarot);
		dtheta = gRandom->Gaus(0,sigmarot);
		dphi = gRandom->Gaus(0,sigmarot);

		UShort_t volid = AliGeomManager::LayerToVolUID(iLayer,iModule);
		const char *symname = AliGeomManager::SymName(volid);
		new(alobj[j++]) AliAlignObjParams(symname, volid, dx, dy, dz, dpsi, dtheta, dphi, kFALSE);
	    }
	}
    }else{
	AliError(Form("\"%s\" is not a valid identifier for misalignment types. Exiting ...",GetMisalType()));
	return 0;
    }

    return array;
}

//_______________________________________________________________________________________
AliCDBMetaData* AliTPCMisAligner::GetCDBMetaData() const {
    // Returns the comment and responsible for the
    // AliCDBMetaData to be associated with the OCDB entry
    // containing the TPC array of misalignment objects
    //
    AliCDBMetaData* md = new AliCDBMetaData();
    md->SetResponsible("");

    if(TString(GetMisalType())=="ideal")
	md->SetComment("Alignment objects for TPC ideal misalignment");
    if(TString(GetMisalType())=="residual")
	md->SetComment("Alignment objects for TPC residual misalignment");
    if(TString(GetMisalType())=="full")
	md->SetComment("Alignment objects for TPC full misalignment");

    return md;
}
