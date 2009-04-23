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
// This class generates misalignment for VZERO. In particular it defines
// the misalignment in the three canonical scenarios: "ideal", "residual"
// and "full".
// It is meant to be run standalone or from the steering macro
// $ALICE_ROOT/macros/MakeAlignmentObjs.C
// looping on the detectors.
//
//========================================================================

#include "AliVZEROMisAligner.h"
#include "AliGeomManager.h"
#include "AliMathBase.h"
#include "TClonesArray.h"
#include "TRandom.h"
#include "AliAlignObjParams.h"
#include "AliLog.h"

ClassImp(AliVZEROMisAligner)

    //_______________________________________________________________________________________
AliVZEROMisAligner::AliVZEROMisAligner() : AliMisAligner()
{
    //
    // dummy constructor
    //
}

//_______________________________________________________________________________________
TClonesArray* AliVZEROMisAligner::MakeAlObjsArray() {
    // builds and returns the array of alignment objects
    // according to the spcified misalignment scenario
    // ("ideal", "residual" or "full").
    //
    TClonesArray *array = new TClonesArray("AliAlignObjParams",2);
    TClonesArray &alobj = *array;

    Double_t dx, dy, dz, dpsi, dtheta, dphi;
    gRandom->SetSeed(4321);
    Double_t sigmatr; // max shift in cm
    Double_t sigmarot; // max rot in degrees

    TString v0alignable[2]={"VZERO/V0C", "VZERO/V0A"};

    Int_t iIndex=0; // VZERO is not indexed
    AliGeomManager::ELayerID iLayer = AliGeomManager::kInvalidLayer;
    UShort_t volid = AliGeomManager::LayerToVolUID(iLayer,iIndex);

    if(TString(GetMisalType())=="ideal")
    {

	for(Int_t ii=0; ii<2; ii++)
	    new(alobj[ii]) AliAlignObjParams(v0alignable[ii].Data(), volid, 0., 0., 0., 0., 0., 0., kTRUE);

    }else if(TString(GetMisalType())=="residual" || TString(GetMisalType())=="full"){

	if(!AliGeomManager::GetGeometry())
	{
	    AliError("No geometry loaded into AliGeomManager! Returning null pointer!");
	    return 0;
	}

	sigmatr = 0.1;
	sigmarot = 0.5;

	for(Int_t ii=0; ii<2; ii++)
	{
	    dx = AliMathBase::TruncatedGaus(0.,sigmatr, 3*sigmatr);
	    dy = AliMathBase::TruncatedGaus(0.,sigmatr, 3*sigmatr);
	    dz = AliMathBase::TruncatedGaus(0.,sigmatr, 3*sigmatr);
	    dpsi   = AliMathBase::TruncatedGaus(0.,sigmarot, 3*sigmarot);
	    dtheta = AliMathBase::TruncatedGaus(0.,sigmarot, 3*sigmarot);
	    dphi   = AliMathBase::TruncatedGaus(0.,sigmarot, 3*sigmarot);
	    new(alobj[ii]) AliAlignObjParams(v0alignable[ii].Data(), volid, dx, dy, dz, dpsi, dtheta, dphi, kFALSE);
	}

    }else{
	AliError(Form("\"%s\" is not a valid identifier for misalignment types. Exiting ...",GetMisalType()));
	return 0;
    }

    return array;
}

//_______________________________________________________________________________________
AliCDBMetaData* AliVZEROMisAligner::GetCDBMetaData() const {
    // Returns the comment and responsible for the
    // AliCDBMetaData to be associated with the OCDB entry
    // containing the VZERO array of misalignment objects
    //
    AliCDBMetaData* md = new AliCDBMetaData();
    md->SetResponsible("Brigitte Cheynis");

    if(TString(GetMisalType())=="ideal")
	md->SetComment("Alignment objects for VZERO ideal misalignment");
    if(TString(GetMisalType())=="residual")
	md->SetComment("Alignment objects for VZERO residual misalignment");
    if(TString(GetMisalType())=="full")
	md->SetComment("Alignment objects for VZERO full misalignment");

    return md;
}
