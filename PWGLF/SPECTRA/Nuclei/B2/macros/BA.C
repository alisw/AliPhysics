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

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TFileMerger.h>
#include <TString.h>
#include "AliLnBA.h"
#endif

Int_t BA(  const TString& pSpectra   = "~/alice/output/Proton-lhc10d-Spectra.root"
         , const TString& dSpectra   = "~/alice/output/Deuteron-lhc10d-Spectra.root"
         , const TString& ptag       = "lhc10d"
         , const TString& dtag       = "lhc10d"
         , const TString& outputfile = "~/alice/output/lhc10d-B2.root"
         , const TString& otag       = "lhc10d"
         , const Bool_t   draw       = 1)
{
//
// coalescence parameter for nucleus and antinucleus
//
	using namespace std;
	
	Int_t A = 2;
	Int_t Z = 1;
	
	if(nucleus == "Deuteron")
	{
		A = 2;
		Z = 1;
	}
	else if(nucleus == "Triton")
	{
		A = 3;
		Z = 1;
	}
	else if(nucleus == "He3")
	{
		A = 3;
		Z = 2;
	}
	else
	{
		cerr << "only valid names: Deuteron, Triton and He3" << endl;
		return 1;
	}
	
	const Int_t kNpart = 2;
	const Int_t kZ[kNpart] = { Z, -Z };
	const TString kB2File[kNpart] = {"BA.root", "AntiBA.root" };
	
	TFileMerger m;
	
	for(Int_t i=0; i<kNpart; ++i)
	{
		AliLnBA b2(pSpectra, ptag, dSpectra, dtag, kB2File[i], otag, A, kZ[i]);
		
		b2.Run();
		
		m.AddFile(kB2File[i].Data(),0);
	}
	
	m.OutputFile(outputfile.Data());
	m.Merge();
	
	gSystem->Exec(Form("rm -f %s %s", kB2File[0].Data(), kB2File[1].Data()));
	
	if(!draw) return 0;
	
	const TString kNucleus[kNpart] = { nucleus, Form("Anti%s",nucleus.Data())};
	
	for(Int_t i=0; i<kNpart; ++i)
	{
		gROOT->ProcessLine(Form(".x DrawBA.C+g(\"%s\",\"%s\",\"%s\")", outputfile.Data(), otag.Data(), kNucleus[i].Data()));
	}
	
	return 0;
}
