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
*/

/*
Old Logs: AliDrawVolume.cxx,v $
Revision 1.1  2000/07/13 16:19:10  fca
Mainly coding conventions + some small bug fixes

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

#include "AliG3Volume.h"
#include <TArrayF.h>
#include "TMaterial.h"
#include "TShape.h"
#include "TTUBE.h"
#include "TBRIK.h"
#include "TTRD1.h"
#include "TTRD2.h"
#include "TTRAP.h"
#include "TTUBS.h"
#include "TCONE.h"
#include "TCONS.h"
#include "TSPHE.h"
#include "TPARA.h"
#include "TPGON.h"
#include "TPCON.h"
#include "TTUBS.h"
#include "TELTU.h"
#include "THYPE.h"
#include "TGTRA.h"
#include "TCTUB.h"

ClassImp(AliG3Volume)

    AliG3Volume::AliG3Volume(const char* name) 
	: TNamed(name, " ")
{
// Constructor

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
    fNParam = 0;
    fPosition.Set(3);
    fCopies = new TList;
    fNCopies = 0;
    fPosp = kFALSE;
    fAxis = -1;
}

AliG3Volume::AliG3Volume(const AliG3Volume& volume)
    : TNamed(volume)
{
// Constructor

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
    fAxis = -1;
//
    fIdVolume   = volume.GetIdVolume();
    fIdCopy     = volume.GetIdCopy();
    fIdMedium   = volume.Medium();
    fIdMaterial = volume.Material();
    fPosition   = volume.Position(0);
    fShape      = volume.Shape();
    fRotMatrix  = volume.RotMatrix();
    TArrayF par;
    volume.Parameters(0, fParameters);
    fNCopies    = volume.NCopies();
    fPosp       = volume.Posp();
    
    fNParam     = volume.NParam();
    fCopies     = volume.Copies();
    fRotMatrix  = volume.RotMatrix();
    volume.Division(fNdiv, fAxis, fStartC, fStep);
    
}




void AliG3Volume::Draw(Option_t *)
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

void AliG3Volume::DrawSpec()
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

void AliG3Volume::SetParam(Int_t ip, Float_t param)
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

Float_t  AliG3Volume::GetParam(Int_t ip)
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

void  AliG3Volume::AddCopy(AliG3Volume* volume)
{
    volume->SetIdMaterial(Material());
    fCopies->Add(volume);
    fNCopies++;
}

AliG3Volume* AliG3Volume::Copy(Int_t i)
{
    return (AliG3Volume*) fCopies->At(i);
}


TArrayF AliG3Volume::Position(Int_t i) const
{
//
// Get position for volume copy i 
//
    if (i==0) {
	return fPosition;
    } else {
	return ((AliG3Volume*) fCopies->At(i-1))->Position(0);
    }
}

void AliG3Volume::SetPosition(Float_t x, Float_t y, Float_t z)
{
//
// Set position
//
  fPosition[0] = x;
  fPosition[1] = y;
  fPosition[2] = z;
}


void  AliG3Volume::SetParameters(Int_t np, Float_t* param) 
{
//
// Set parameters 
//
    fParameters.Set(np);
    for (Int_t j=0; j<np; j++) fParameters[j]=param[j];
    fNParam = np;
}
    

void  AliG3Volume::Parameters(Int_t i, TArrayF& param) const
{
//
// Get parameters for volume copy i 
//
    TArrayF p;
    if (i==0) {
	p = fParameters;
    } else {
	((AliG3Volume*) (fCopies->At(i-1)))->Parameters(0, p);
    }
    Int_t np = fNParam;
    param.Set(np);
    for (Int_t j=0; j<np; j++) {
	param[j] = p.At(j);
    }
}
    

void AliG3Volume::CreateTShape(char* nameV, TMaterial* mat)
{
//
// Create a root volume from G3 volume
//
    Int_t   ip;
    
    Float_t kRadDeg = 180./TMath::Pi();
    Float_t theta, phi, alpha1, alpha2;
    Float_t p1, p2;
    
    TShape* nShape=0;
    const char* tmp = mat->GetName();
    char nameM[21];
    strncpy(nameM, tmp, 20);
    nameM[20]='\0';
    switch(fShape)
    {
    case 1:
// BOX
	nShape = new TBRIK(nameV,"BRIK",nameM,fParameters[0], fParameters[1], fParameters[2]);
	break;
	
    case 2:
// TRD1	  
	nShape = new TTRD1(nameV, "TRD1", nameM, fParameters[0], fParameters[1], fParameters[2],
			   fParameters[3]);
	break;
	
    case 3:
// TRD2
	nShape = new TTRD2(nameV, "TRD2", nameM, fParameters[0], fParameters[1], fParameters[2], 
			   fParameters[3], fParameters[4]);
	break;
	
    case 4:
// TRAP
	p1     = fParameters[1];
	p2     = fParameters[2];
	
	theta  = TMath::ATan(TMath::Sqrt(p1*p1+p2*p2))*kRadDeg;
	phi    = TMath::ATan2(p2,p1)*kRadDeg;
	alpha1 = fParameters[6 ]*kRadDeg;
	alpha2 = fParameters[10]*kRadDeg;

	if (theta < 0.) theta+=180.;
	
	nShape =  new TTRAP(nameV, "TRAP", nameM, fParameters[0], 
			    theta, phi,
			    fParameters[3], fParameters[4], fParameters[5], 
			    alpha1, 
			    fParameters[7], fParameters[8], fParameters[9], 
			    alpha2);
	break;
	
    case 5:
// TUBE
	nShape =  new TTUBE(nameV,"TUBE",nameM,fParameters[0], fParameters[1], fParameters[2]);
	break;
	
    case 6:
// TUBS
	nShape =  new TTUBS(nameV,"TUBS",nameM,fParameters[0], fParameters[1], fParameters[2], 
			    fParameters[3], fParameters[4]);
	break;
	
    case 7:
// CONE
	nShape =  new TCONE(nameV, "CONE", nameM, fParameters[0], fParameters[1], fParameters[2], 
			    fParameters[3], fParameters[4]);
	break;
	
    case 8:
	
// CONS
	nShape =  new TCONS(nameV, "CONS", nameM, fParameters[0], fParameters[1], fParameters[2], 
			    fParameters[3], fParameters[4], fParameters[5], fParameters[6]);
	break;
	
    case 9:
// SPHE
	
	nShape =  new TSPHE(nameV, "SPHE", nameM, fParameters[0], fParameters[1], fParameters[2], 
			    fParameters[3], fParameters[4], fParameters[5]);
	break;
	
    case 10:
// PARA
	alpha1 = fParameters[3]*kRadDeg;
	p1     = fParameters[4];
	p2     = fParameters[5];
	theta  = TMath::ATan(TMath::Sqrt(p1*p1+p2*p2))*kRadDeg;
	phi    = TMath::ATan2(p2,p1)*kRadDeg;

	nShape =  new TPARA(nameV, "PARA", nameM, fParameters[0], fParameters[1], fParameters[2], 
				    alpha1, theta, phi);
	break;
	
    case 11:
// PGON
	nShape =  new TPGON(nameV, "PGON", nameM, fParameters[0], fParameters[1], Int_t(fParameters[2]), 
			    Int_t(fParameters[3]));
	for (ip=0; ip<Int_t(fParameters[3]); ip++) {
	    ((TPGON*) nShape)->DefineSection(ip, fParameters[4+3*ip], fParameters[4+3*ip+1], 
					     fParameters[4+3*ip+2]);
	}
	break;
	
    case 12:
// PCON
	nShape = new TPCON(nameV, "PCON", nameM, fParameters[0], fParameters[1], Int_t(fParameters[2]));
	for (ip=0; ip<Int_t(fParameters[2]); ip++) {
	    ((TPCON*) nShape)->DefineSection(ip, fParameters[3+3*ip], fParameters[3+3*ip+1], 
					     fParameters[3+3*ip+2]);
	}
	break;
	
    case 13:
// ELTU  
	nShape = new TELTU(nameV,"ELTU",nameM,fParameters[0], fParameters[1], fParameters[2]);
	break;
	
    case 14:
// HYPE
	nShape = new THYPE(nameV,"HYPE",nameM,fParameters[0], fParameters[1], fParameters[2], 
			   fParameters[3]);
	break;
	
    case 15:
// GTRA
	nShape = new TGTRA(nameV, "GTRA", nameM, fParameters[0], fParameters[1], fParameters[2], 
			   fParameters[3], fParameters[4], fParameters[5], fParameters[6], 
			   fParameters[7], fParameters[8], fParameters[9], fParameters[10], 
			   fParameters[11]);
	break;
	
    case 16:
// CTUB
	nShape = new TCTUB(nameV, "CTUB", nameM, fParameters[0], fParameters[1], fParameters[2], 
			   fParameters[3], fParameters[4], fParameters[5], fParameters[6], 
			   fParameters[7], fParameters[8], fParameters[9], fParameters[10]);
	break;
    default:
	break;
    }
    if (nShape) {
	Float_t density = mat->GetDensity();
	if (density < 0.01) {
	    nShape->SetVisibility(0);
	} else {
	    nShape->SetVisibility(1);
	}
	
	Int_t color = Int_t(density/20.*100.);
	nShape->SetLineColor(color);
    }
}

void   AliG3Volume::SetDivision(Int_t ndiv, Int_t axis, Float_t start, Float_t step)
{
    fNdiv   = ndiv;
    fAxis   = axis;
    fStartC = start;
    fStep   = step;
}

void   AliG3Volume::Division(Int_t& ndiv, Int_t& axis, Float_t& start, Float_t& step) const
{
    ndiv  = fNdiv;
    axis  = fAxis;
    start = fStartC;
    step  = fStep;
}
