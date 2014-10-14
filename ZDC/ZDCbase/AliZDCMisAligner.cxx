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
// This class generates misalignment for ZDC. In particular it defines
// the misalignment in the three canonical scenarios: "ideal", "residual"
// and "full".
// It is meant to be run standalone or from the steering macro
// $ALICE_ROOT/macros/MakeAlignmentObjs.C
// looping on the detectors.
//
//========================================================================

#include "AliZDCMisAligner.h"
#include "AliGeomManager.h"
#include "TClonesArray.h"
#include "AliAlignObjParams.h"
#include "AliLog.h"

ClassImp(AliZDCMisAligner)

    //_______________________________________________________________________________________
AliZDCMisAligner::AliZDCMisAligner() : AliMisAligner()
{
    //
    // dummy constructor
    //
}

//_______________________________________________________________________________________
TClonesArray* AliZDCMisAligner::MakeAlObjsArray() {
    // builds and returns the array of alignment objects
    // according to the spcified misalignment scenario
    // ("ideal", "residual" or "full").
    //
    TClonesArray *array = new TClonesArray("AliAlignObjParams",10);
    TClonesArray &alobj = *array;

    Double_t dx,dy,dz,dpsi,dtheta,dphi;
    if(TString(GetMisalType())=="ideal")
    {
	dx=0., dy=0., dz=0.;
	dpsi=0., dtheta=0., dphi=0.;
    }else if(TString(GetMisalType())=="residual" || TString(GetMisalType())=="full")
    {
	dx=0., dy=0.05, dz=0.;
	dpsi=0., dtheta=0., dphi=0.;
    }else{
	AliError(Form("\"%s\" is not a valid identifier for misalignment types. Exiting ...",GetMisalType()));
	return 0;
    }

    const char *zdcCn="ZDC/NeutronZDC_C";
    const char *zdcCp="ZDC/ProtonZDC_C";
    const char *zdcAn="ZDC/NeutronZDC_A";
    const char *zdcAp="ZDC/ProtonZDC_A";

    UShort_t iIndex=0;
    AliGeomManager::ELayerID iLayer = AliGeomManager::kInvalidLayer;
    UShort_t volid = AliGeomManager::LayerToVolUID(iLayer,iIndex);

    new(alobj[0]) AliAlignObjParams(zdcCn, volid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE);
    new(alobj[1]) AliAlignObjParams(zdcCp, volid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE);
    new(alobj[2]) AliAlignObjParams(zdcAn, volid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE);
    new(alobj[3]) AliAlignObjParams(zdcAp, volid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE);

    return array;
}

//_______________________________________________________________________________________
AliCDBMetaData* AliZDCMisAligner::GetCDBMetaData() const {
    // Returns the comment and responsible for the
    // AliCDBMetaData to be associated with the OCDB entry
    // containing the ZDC array of misalignment objects
    //
    AliCDBMetaData* md = new AliCDBMetaData();
    md->SetResponsible("Chiara Oppedisano");

    if(TString(GetMisalType())=="ideal")
	md->SetComment("Alignment objects for ZDC ideal misalignment");
    if(TString(GetMisalType())=="residual")
	md->SetComment("Alignment objects for ZDC residual misalignment");
    if(TString(GetMisalType())=="full")
	md->SetComment("Alignment objects for ZDC full misalignment");

    return md;
}
