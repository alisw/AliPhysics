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
Revision 1.1  2001/07/09 11:41:46  morsch
AliGUIMedium, AliGUIMaterial and AliDrawVolume obselete. Development continues
on AliG3Material, AliG3Medium and AliG3Volume.

*/

/*
Old Logs: AliGUIMaterial.cxx,v $
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

#include "AliG3Material.h"

ClassImp(AliG3Material)
AliG3Material::AliG3Material(char* name, char* title,
			       Float_t a, Float_t z, Float_t dens, Float_t radl, Float_t intl):
    TMaterial(name, title, a, z, dens, radl, intl)
{
    fId=-1;
}


void AliG3Material::Dump()
{
// Dump material information
    printf("\n *****************************************");
    printf("\n Material Number:   %10d", fId);
    printf("\n %s", GetName());
    printf("\n Mass   Number:     %10.2f", fA);    
    printf("\n Charge Number:     %10.2f", fZ);
    printf("\n Density:           %10.2f", fDensity);
    printf("\n Radiation  Length: %10.2f", fRadLength);
    printf("\n Absorption Length: %10.2f", fInterLength);        	
}















