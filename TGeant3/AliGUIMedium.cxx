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

#include "AliGUIMedium.h"

ClassImp(AliGUIMedium)

AliGUIMedium::AliGUIMedium()
{ 
// constructor
    fId=-1;
    fName = 0; 
}

AliGUIMedium::AliGUIMedium(Int_t imed, Int_t imat, char* name, Int_t isvol, 
			   Int_t ifield,
			   Float_t fieldm, Float_t tmaxfd, Float_t stemax, Float_t deemax,
			   Float_t epsil, Float_t stmin)
{
// constructor
    fId=imed;
    fIdMat=imat;
    fName=name;
    fIsvol=isvol;
    fIfield=ifield;
    fFieldm=fieldm;
    fTmaxfd=tmaxfd;
    fStemax=stemax;
    fDeemax=deemax;
    fEpsil=epsil;
    fStmin=stmin;
}

void AliGUIMedium::Dump()
{
// Dummy dump
    ;
}

Int_t AliGUIMedium::Id()
{
// return medium id
    return fId;
}

char*  AliGUIMedium::Name()
{
// return medium name
    return fName;
}

Float_t AliGUIMedium::GetPar(Int_t ipar)
{ 
// Get parameter number ipar
    Float_t p;
    if (ipar < 23) {
	p= fPars[ipar-1];   
    } else if(ipar >=23 && ipar <27) {
	p= fPars[ipar-1+3];   
    } else {
	p= fPars[ipar-1+4];   
    }
    
    return p;
}
 
void AliGUIMedium::Streamer(TBuffer &)
{
// dummy streamer
;
}















