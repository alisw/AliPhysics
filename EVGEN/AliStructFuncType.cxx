/**************************************************************************
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
Revision 1.2.4.1  2002/11/26 16:57:23  hristov
Merging NewIO with v3-09-04

Revision 1.2  2002/10/14 14:55:35  hristov
Merging the VirtualMC branch to the main development branch (HEAD)

Revision 1.1.2.1  2002/07/24 08:56:29  alibrary
Updating EVGEN on TVirtulaMC

Revision 1.1  2002/07/22 10:19:12  morsch
Interface to pdfset and structa added.

*/


#include "AliStructFuncType.h"

#ifndef WIN32
# define structa  structa_ 
# define pdfset   pdfset_
# define type_of_call
#else
# define structa  STRUCTA
# define pdfset   PDFSET
# define type_of_call _stdcall
#endif


extern "C" {
    void type_of_call pdfset(char parm[20][20], Double_t value[20]);
    
    void type_of_call structa(Double_t& xx, Double_t& qq, Double_t& a,
			       Double_t& upv, Double_t& dnv, Double_t& usea,
			       Double_t& dsea,
			       Double_t& str, Double_t& chm, Double_t& bot,
			       Double_t& top, Double_t& gl);
}

ClassImp(AliStructFuncType)

void AliStructFuncType::PdfSet(char parm[20][20], Double_t value[20])
{
    pdfset(parm, value);
}

void AliStructFuncType::StructA(Double_t xx, Double_t qq, Double_t a,
				Double_t& upv, Double_t& dnv, Double_t& usea,
				Double_t& dsea,
				Double_t& str, Double_t& chm, Double_t& bot,
				Double_t& top, Double_t& gl)
{
    structa(xx, qq, a, upv, dnv, usea, dsea, str, chm, bot, top, gl);
}


