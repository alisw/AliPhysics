/* *************************************************************************
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

/*
$Log$
Revision 1.8  2000/07/12 08:56:32  fca
Coding convention correction and warning removal

Revision 1.7  2000/06/28 21:27:45  morsch
Most coding rule violations corrected.
Still to do: Split the file (on file per class) ? Avoid the global variables.
Copy constructors and assignment operators (dummy ?)

Revision 1.6  2000/04/14 11:07:46  morsch
Correct volume to medium assignment in case several media are asigned to the
same material.

Revision 1.5  2000/03/20 15:11:03  fca
Mods to make the code compile on HP

Revision 1.4  2000/01/18 16:12:08  morsch
Bug in calculation of number of volume divisions and number of positionings corrected
Browser for Material and Media properties added

Revision 1.3  1999/11/14 14:31:14  fca
Correct small error and remove compilation warnings on HP

Revision 1.2  1999/11/10 16:53:35  fca
The new geometry viewer from A.Morsch

*/

/* 
 *  Version: 0
 *  Written by Andreas Morsch
 *  
 * 
 *
 * For questions critics and suggestions to this part of the code
 * contact andreas.morsch@cern.ch
 * 
 **************************************************************************/

#include "AliDrawVolume.h"

ClassImp(AliDrawVolume)

AliDrawVolume::AliDrawVolume(char* name)
{
// Constructor
    fName   = name;
    fTheta  = 30;
    fPhi    = 30;
    fPsi    = 0;
    fU      = 10;
    fV      = 10;
    fUscale = 0.01;
    fVscale = 0.01;
    fHide=0;
    fShadow=0;
    fFill=1;
    fSeen=1;
    fClip=0;
    fClipXmin=0.;
    fClipXmax=2000.;
    fClipYmin=0.;
    fClipYmax=2000.;
    fClipZmin=0.;
    fClipZmax=2000.;
}

char* AliDrawVolume::Name()
{
//
// Return volume name
    return fName;
}

    
void AliDrawVolume::Streamer(TBuffer &)
{
// Dummy Streamer
;
}



void AliDrawVolume::Draw(Option_t *)
{
// Wraps the geant Gdraw
    gMC->Gsatt(fName,"seen", fSeen);
    
    if (fHide) {
	gMC->Gdopt("hide", "on");
    } else {
	gMC->Gdopt("hide", "off");
    }

    if (fShadow) {
	gMC->Gdopt("shad", "on");
	gMC->Gsatt("*", "fill", fFill);
    } else {
	gMC->Gdopt("shad", "off");
    }

    	gMC->SetClipBox(".");
    if (fClip) {
	gMC->SetClipBox("*", fClipXmin, fClipXmax, 
			fClipYmin, fClipYmax, fClipZmin, fClipZmax);
    } else {
	gMC->SetClipBox(".");
    }
    

    gMC->Gdraw(fName, fTheta, fPhi, fPsi, fU, fV, fUscale, fVscale);
    THIGZ *higz = (THIGZ*)gROOT->GetListOfCanvases()->FindObject("higz");
    if (higz) higz->Update();
}

void AliDrawVolume::DrawSpec()
{
// Wraps the Geant DrawSpec
    gMC->Gsatt(fName,"seen", fSeen);
    
    if (fHide) {
	gMC->Gdopt("hide", "on");
    } else {
	gMC->Gdopt("hide", "off");
    }

    if (fShadow) {
	gMC->Gdopt("shad", "on");
	gMC->Gsatt("*", "fill", fFill);
    } else {
	gMC->Gdopt("shad", "off");
    }

    gMC->SetClipBox(".");
    if (fClip) {
	gMC->SetClipBox("*", fClipXmin, fClipXmax, fClipYmin, fClipYmax, fClipZmin, fClipZmax);
    } else {
	gMC->SetClipBox(".");
    }
    

    ((TGeant3*) gMC)->DrawOneSpec(fName);
    THIGZ *higz = (THIGZ*)gROOT->GetListOfCanvases()->FindObject("higz");
    if (higz) higz->Update();
}

void AliDrawVolume::SetParam(Int_t ip, Float_t param)
{
// Set drawing parameters
    switch (ip) {
    case kTheta:
	fTheta=param;
	break;
    case kPhi:
	fPhi=param;
	break;
    case kPsi:
	fPsi=param;
	break;
    case kU:
	fU=param;
	break;
    case kV:
	fV=param;
	break;
    case kUscale:
	fUscale=param;
	break;
    case kVscale:
	fVscale=param;
	break;
    case kHide:
	fHide=Int_t(param);
	break;
    case kShadow:
	fShadow=Int_t(param);
	break;
    case kFill:
	fFill=Int_t(param);
	break;
    case kSeen:
	fSeen=Int_t(param);
	break;
    case kClip:
	fClip=Int_t(param);
	break;
    case kClipXmin:
	fClipXmin=param;
	break;
    case kClipXmax:
	fClipXmax=param;
	break;
    case kClipYmin:
	fClipYmin=param;
	break;
    case kClipYmax:
	fClipYmax=param;
	break;
    case kClipZmin:
	fClipZmin=param;
	break;
    case kClipZmax:
	fClipZmax=param;
	break;
    }
}

Float_t  AliDrawVolume::GetParam(Int_t ip)
{
// Get drawing parameters
    switch (ip) {
    case kTheta:
	return fTheta;
    case kPhi:
	return fPhi;
    case kPsi:
	return fPsi;
    case kU:
	return fU;
    case kV:
	return fV;
    case kUscale:
	return fUscale;
    case kVscale:
	return fVscale;
    case kHide:
	return Float_t(fHide);
    case kShadow:
	return Float_t(fShadow);
    case kFill:
	return Float_t(fFill);
    case kSeen:
	return Float_t(fSeen);
    case kClip:
	return Float_t(fClip);
    case kClipXmin:
	return fClipXmin;
    case kClipXmax:
	return fClipXmax;
    case kClipYmin:
	return fClipYmin;
    case kClipYmax:
	return fClipYmax;
    case kClipZmin:
	return fClipZmin;
    case kClipZmax:
	return fClipZmax;
    default:
	return 0.;
    }
    return 0.;
}


