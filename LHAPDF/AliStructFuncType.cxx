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

/* $Id$ */

// Helper class to interface pdflib and the TPythia 
// the c++ interface for Pythia
// The pdf codes used  in pdflib are mapped
// to a more user friendly enumeration type.
// Author: andreas.morsch@cern.ch

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


Int_t AliStructFuncType::PDFsetIndex(StrucFunc_t pdf)
{
// PDF set index
    Int_t pdfSetNumber[12] = {
	19170,
	19150,
	19070,
	19050,
	80060,
	10040,
	10100,
	10050,
	10041,
	10042,
        10800,
	11000
    };
    return pdfSetNumber[pdf];
}

TString AliStructFuncType::PDFsetName(StrucFunc_t pdf)
{
// PDF Set Name
    TString pdfsetName[12]   = {
	"cteq4l.LHgrid",
	"cteq4m.LHgrid",   
	"cteq5l.LHgrid",
	"cteq5m.LHgrid",
	"GRV98lo.LHgrid",   
	"cteq6.LHpdf",
	"cteq61.LHpdf",
	"cteq6m.LHpdf",
	"cteq6l.LHpdf", 
	"cteq6ll.LHpdf",
	"CT10.LHgrid", 
	"CT10nlo.LHgrid"
    };
    return pdfsetName[pdf];
}
