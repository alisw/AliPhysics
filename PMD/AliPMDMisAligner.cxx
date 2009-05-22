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
// This class generates misalignment for PMD. In particular it defines
// the misalignment in the three canonical scenarios: "ideal", "residual"
// and "full".
// It is meant to be run standalone or from the steering macro
// $ALICE_ROOT/macros/MakeAlignmentObjs.C
// looping on the detectors.
//
// Macro to randomly displace the 4 sectors of the PMD
// in each plane. Each sector (to be misaligned) 
// of PMD houses the following :
// (a) 6 modules of preshower plane
// (b) 6 modules of veto plane
// (c) The FEE boards on back plane of each module
// (d) 6 modules of convertor plates
// The clustering is done module-wise
// The actual amount displacement will be provided
// by the survey data and has to be converted into
// displacement in x,y,z,theta, phi and psi 


// Now specify the path of the module to be misaligned
// as followed in the PMD geant

//   _____________
//  |    |        |
//  | 1  |   3    |
//  |    |________|
//  |____|___|    |
//  |        | 2  |
//  |   4    |    |
//  |________|____|
  
// Misalignment Matrix is expected to be
// same for sectors 1 and 4 
// and for the sectors 2 and 3
// As these will be mounted on the same
// Steel plate 
//========================================================================

#include "AliPMDMisAligner.h"
#include "AliGeomManager.h"
#include "TClonesArray.h"
#include "TRandom.h"
#include "AliAlignObjParams.h"
#include "AliLog.h"

ClassImp(AliPMDMisAligner)

    //_______________________________________________________________________________________
AliPMDMisAligner::AliPMDMisAligner() : AliMisAligner()
{
    //
    // dummy constructor
    //
}

//_______________________________________________________________________________________
TClonesArray* AliPMDMisAligner::MakeAlObjsArray() {
    // builds and returns the array of alignment objects
    // according to the spcified misalignment scenario
    // ("ideal", "residual" or "full").
    //
    TClonesArray *array = new TClonesArray("AliAlignObjParams",4);
    TClonesArray &alobj = *array;

    Double_t max_trans=0.1; // maximun shifts in X,Y,Z  in centimeters
    Double_t max_rot=0.1;   // maximum shifts in angles in degrees

    const char *Sector1="PMD/Sector1"; 
    const char *Sector2="PMD/Sector2"; 
    const char *Sector3="PMD/Sector3"; 
    const char *Sector4="PMD/Sector4"; 

    //Sectors 1 and 4
    Double_t dx14, dy14, dz14;          // Misalignment in X,Y and Z
    Double_t dpsi14, dtheta14, dphi14; //  Angular displacements
    //Sectors 2 and 3
    Double_t dx23, dy23, dz23;          // Misalignment in X,Y and Z
    Double_t dpsi23, dtheta23, dphi23; //  Angular displacements

    UShort_t iIndex=0; // PMD is not indexed
    AliGeomManager::ELayerID iLayer = AliGeomManager::kInvalidLayer;
    UShort_t volid = AliGeomManager::LayerToVolUID(iLayer,iIndex);
  
    gRandom->SetSeed(4357);
  
    if(TString(GetMisalType())=="ideal")
    {

	Double_t dx=0., dy=0., dz=0., dpsi=0., dtheta=0., dphi=0.;

	for(Int_t i=0; i<4; i++){
	    TString snSector(Form("PMD/Sector%d",i+1));
	    new(alobj[i]) AliAlignObjParams(snSector.Data(), volid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE);
	}

    }else if(TString(GetMisalType())=="residual" || TString(GetMisalType())=="full"){
	if(!AliGeomManager::GetGeometry())
	{
	    AliError("No geometry loaded into AliGeomManager! Returning null pointer!");
	    return 0;
	}

	// For sectors 1 and 4
	// Translation displacement
	dx14     = (gRandom->Uniform()-0.5)*max_trans;
	dy14     = (gRandom->Uniform()-0.5)*max_trans;
	dz14     = (gRandom->Uniform()-0.5)*max_trans;
	//Rotation angles
	dpsi14   = (gRandom->Uniform()-0.5)*max_rot;
	dtheta14 = (gRandom->Uniform()-0.5)*max_rot;
	dphi14   = (gRandom->Uniform()-0.5)*max_rot;

	// For sectors 2 and 3
	// Translation displacement
	dx23     = (gRandom->Uniform()-0.5)*max_trans;
	dy23     = (gRandom->Uniform()-0.5)*max_trans;
	dz23     = (gRandom->Uniform()-0.5)*max_trans;
	//Rotation angles
	dpsi23   = (gRandom->Uniform()-0.5)*max_rot;
	dtheta23 = (gRandom->Uniform()-0.5)*max_rot;
	dphi23   = (gRandom->Uniform()-0.5)*max_rot;

	new(alobj[0]) AliAlignObjParams(Sector1, volid, dx14, dy14, dz14, dpsi14, dtheta14, dphi14, kFALSE);
	new(alobj[1]) AliAlignObjParams(Sector2, volid, dx14, dy14, dz14, dpsi14, dtheta14, dphi14, kFALSE);
	new(alobj[2]) AliAlignObjParams(Sector3, volid, dx23, dy23, dz23, dpsi23, dtheta23, dphi23, kFALSE);
	new(alobj[3]) AliAlignObjParams(Sector4, volid, dx23, dy23, dz23, dpsi23, dtheta23, dphi23, kFALSE);

    }else{
	AliError(Form("\"%s\" is not a valid identifier for misalignment types. Exiting ...",GetMisalType()));
	return 0;
    }

    return array;
}

//_______________________________________________________________________________________
AliCDBMetaData* AliPMDMisAligner::GetCDBMetaData() const {
    // Returns the comment and responsible for the
    // AliCDBMetaData to be associated with the OCDB entry
    // containing the PMD array of misalignment objects
    //
    AliCDBMetaData* md = new AliCDBMetaData();
    md->SetResponsible("");

    if(TString(GetMisalType())=="ideal")
	md->SetComment("Alignment objects for PMD ideal misalignment");
    if(TString(GetMisalType())=="residual")
	md->SetComment("Alignment objects for PMD residual misalignment");
    if(TString(GetMisalType())=="full")
	md->SetComment("Alignment objects for PMD full misalignment");

    return md;
}
