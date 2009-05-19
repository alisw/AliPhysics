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
// This class generates misalignment for HMPID. In particular it defines
// the misalignment in the three canonical scenarios: "ideal", "residual"
// and "full".
// It is meant to be run standalone or from the steering macro
// $ALICE_ROOT/macros/MakeAlignmentObjs.C
// looping on the detectors.
//
//========================================================================

#include "AliHMPIDMisAligner.h"
#include "AliGeomManager.h"
#include "TClonesArray.h"
#include "TRandom.h"
#include "TMath.h"
#include "AliAlignObjMatrix.h"
#include "AliLog.h"

ClassImp(AliHMPIDMisAligner)

//_______________________________________________________________________________________
AliHMPIDMisAligner::AliHMPIDMisAligner() : AliMisAligner()
{
    //
    // dummy constructor
    //
}

//_______________________________________________________________________________________
TClonesArray* AliHMPIDMisAligner::MakeAlObjsArray() {
    // builds and returns the array of alignment objects
    // according to the spcified misalignment scenario
    // ("ideal", "residual" or "full").
    //
    if(!AliGeomManager::GetGeometry())
    {
	AliError("No geometry loaded into AliGeomManager! Returning null pointer!");
	return 0;
    }

    TClonesArray *array = new TClonesArray("AliAlignObjMatrix",7);

    AliGeomManager::ELayerID idHMPID =  AliGeomManager::kHMPID;
    Double_t sigmaTrans, sigmaRot;
    Double_t dX=0., dY=0., dZ=0., dPsi=0., dTheta=0., dPhi=0.; //displacements

    // TRandom *pRnd   = new TRandom(4357);
    gRandom->SetSeed(4357);

    if(TString(GetMisalType())=="ideal")
    {

	for (Int_t iCh = 0; iCh < 7; iCh++) {
	    new((*array)[iCh]) AliAlignObjMatrix(AliGeomManager::SymName(idHMPID,iCh),
		    AliGeomManager::LayerToVolUID(idHMPID,iCh),dX,dY,dZ,dPsi,dTheta,dPhi,kTRUE);
	}

    }else if( TString(GetMisalType())=="residual" || TString(GetMisalType())=="full"){

	sigmaTrans=0.1; // 1mm
	sigmaRot=0.001*180/TMath::Pi(); // 1 mrad

	for (Int_t iCh = 0; iCh < 7; iCh++) {
	    dX     = (gRandom->Uniform()-0.5)*sigmaTrans;
	    dY     = (gRandom->Uniform()-0.5)*sigmaTrans;
	    dZ     = (gRandom->Uniform()-0.5)*sigmaTrans;
	    dPsi   = (gRandom->Uniform()-0.5)*sigmaRot;
	    dTheta = (gRandom->Uniform()-0.5)*sigmaRot;
	    dPhi   = (gRandom->Uniform()-0.5)*sigmaRot;
	    //    dX     = (pRnd->Uniform()-0.5)*sigmaTrans;    dY     = (pRnd->Uniform()-0.5)*sigmaTrans;    dZ     = (pRnd->Uniform()-0.5)*sigmaTrans;
	    //    dPsi   = (pRnd->Uniform()-0.5)*sigmaRot;    dTheta = (pRnd->Uniform()-0.5)*sigmaRot;    dPhi   = (pRnd->Uniform()-0.5)*sigmaRot;
	    new((*array)[iCh]) AliAlignObjMatrix(AliGeomManager::SymName(idHMPID,iCh),
		    AliGeomManager::LayerToVolUID(idHMPID,iCh),dX,dY,dZ,dPsi,dTheta,dPhi,kTRUE);
	}
    }else{
	AliError(Form("\"%s\" is not a valid identifier for misalignment types. Exiting ...",GetMisalType()));
	return 0;
    }

    return array;
}

//_______________________________________________________________________________________
AliCDBMetaData* AliHMPIDMisAligner::GetCDBMetaData() const {
    // Returns the comment and responsible for the
    // AliCDBMetaData to be associated with the OCDB entry
    // containing the HMPID array of misalignment objects
    //
    AliCDBMetaData* md = new AliCDBMetaData();
    md->SetResponsible("HMPID expert means nothing");

    if(TString(GetMisalType())=="ideal")
	md->SetComment("Alignment objects for HMPID ideal misalignment");
    if(TString(GetMisalType())=="residual")
	md->SetComment("Alignment objects for HMPID residual misalignment");
    if(TString(GetMisalType())=="full")
        md->SetComment("Full misalignment objects for HMPID produced with sigmaTrans=1mm and sigmaRot=1mrad");

    return md;
}
