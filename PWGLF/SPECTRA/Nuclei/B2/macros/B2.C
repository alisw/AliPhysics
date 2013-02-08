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

// coalescence parameter
// author: Eulogio Serradilla <eulogio.serradilla@cern.ch>

#include <TSystem.h>
#include <TROOT.h>
#include <TFileMerger.h>
#include <TString.h>

#include "AliLnB2.h"

Int_t B2(  const TString& pSpectra   = "~/alice/output/Proton-lhc10d-Spectra.root"
         , const TString& dSpectra   = "~/alice/output/Deuteron-lhc10d-Spectra.root"
         , const TString& ptag       = "lhc10d"
         , const TString& dtag       = "lhc10d"
         , const TString& outputfile = "~/alice/output/lhc10d-B2.root"
         , const TString& otag       = "lhc10d"
         , const Bool_t   draw       = 1)
{
//
// coalescence parameter
//
	const Int_t kNpart      = 2;
	const TString kPrefix[] = { "", "Anti" };
	const Int_t kCharge[]   = { 1, -1 };
	
	TFileMerger m;
	
	for(Int_t i=0; i<kNpart; ++i)
	{
		TString b2file = kPrefix[i] + "B2.root";
		
		AliLnB2 b2(pSpectra, ptag, dSpectra, dtag, b2file, otag, 2, kCharge[i]);
		
		b2.Run();
		
		m.AddFile(b2file.Data(),0);
	}
	
	m.OutputFile(outputfile.Data());
	m.Merge();
	
	gSystem->Exec("rm -f B2.root AntiB2.root");
	
	// draw
	
	if(!draw) return 0;
	
	for(Int_t i=0; i<kNpart; ++i)
	{
		gROOT->ProcessLine(Form(".x DrawB2.C+g(\"%s\",\"%s\",\"%s\")", outputfile.Data(), otag.Data(), kPrefix[i].Data()));
	}
	
	return 0;
}
