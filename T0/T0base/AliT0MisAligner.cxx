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
// This class generates misalignment for T0. In particular it defines
// the misalignment in the three canonical scenarios: "ideal", "residual"
// and "full".
// It is meant to be run standalone or from the steering macro
// $ALICE_ROOT/macros/MakeAlignmentObjs.C
// looping on the detectors.
//
//========================================================================

#include "AliT0MisAligner.h"
#include "AliGeomManager.h"
#include "TClonesArray.h"
#include "TRandom.h"
#include "AliAlignObjParams.h"
#include "AliLog.h"

ClassImp(AliT0MisAligner)

    //_______________________________________________________________________________________
AliT0MisAligner::AliT0MisAligner() : AliMisAligner()
{
    //
    // dummy constructor
    //
}

//_______________________________________________________________________________________
TClonesArray* AliT0MisAligner::MakeAlObjsArray() {
    // builds and returns the array of alignment objects
    // according to the spcified misalignment scenario
    // ("ideal", "residual" or "full").
    //
    TClonesArray *array = new TClonesArray("AliAlignObjParams",4);
    TClonesArray &alobj = *array;

    Double_t dx,dy,dz,dpsi,dtheta,dphi;
    gRandom->SetSeed(4321);
    Double_t sigmatr = 0.006; // sigma for shifts in cm
    Double_t sigmarot = 0.001; // sigma for tilts in degrees

    TString symName[2] = {"/ALIC_1/0STR_1","/ALIC_1/0STL_1"};

    Int_t iIndex=0;
    AliGeomManager::ELayerID iLayer = AliGeomManager::kInvalidLayer;
    UShort_t volid = AliGeomManager::LayerToVolUID(iLayer,iIndex);

    if(TString(GetMisalType())=="ideal")
    {
	dx=0., dy=0., dz=0.;
	dpsi=0., dtheta=0., dphi=0.;
	for (Int_t imod=0; imod<2; imod++)
	{
	    new(alobj[imod]) AliAlignObjParams(symName[imod].Data(), volid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE);
	}
    }else if(TString(GetMisalType())=="residual" || TString(GetMisalType())=="full")
    {

	for (Int_t imod=0; imod<2; imod++)
	{
	    dx = gRandom->Gaus(0.,sigmatr);
	    dy = gRandom->Gaus(0.,sigmatr);
	    dz = gRandom->Gaus(0.,sigmatr);
	    dpsi = gRandom->Gaus(0.,sigmarot);
	    dtheta = gRandom->Gaus(0.,sigmarot);
	    dphi = gRandom->Gaus(0.,sigmarot);
	    new(alobj[imod]) AliAlignObjParams(symName[imod].Data(), volid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE);
	}

    }else{
	AliError(Form("\"%s\" is not a valid identifier for misalignment types. Exiting ...",GetMisalType()));
	return 0;
    }

    return array;
}

//_______________________________________________________________________________________
AliCDBMetaData* AliT0MisAligner::GetCDBMetaData() const {
    // Returns the comment and responsible for the
    // AliCDBMetaData to be associated with the OCDB entry
    // containing the T0 array of misalignment objects
    //
    AliCDBMetaData* md = new AliCDBMetaData();
    md->SetResponsible("Tomasz Malkiewicz");

    if(TString(GetMisalType())=="ideal")
	md->SetComment("Alignment objects for T0 ideal misalignment");
    if(TString(GetMisalType())=="residual")
	md->SetComment("Alignment objects for T0 residual misalignment");
    if(TString(GetMisalType())=="full")
	md->SetComment("Alignment objects for T0 full misalignment");

    return md;
}
