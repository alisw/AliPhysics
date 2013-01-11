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

// driver for rebuilding the correction file
// author: Eulogio Serradilla <eulogio.serradilla@cern.ch>

#include <Riostream.h>
#include <TSystem.h>
#include <TString.h>
#include <TROOT.h>
#include <TFileMerger.h>

#include "AliLnUnfolding.h"
#include "AliLnFakeTracks.h"
#include "AliLnSecondaries.h"
#include "AliLnEfficiency.h"
#include "AliLnCorr.h"

ClassImp(AliLnCorr)

AliLnCorr::AliLnCorr(const TString& particle, const TString& dataFilename, const TString& simuFilename, const TString& simuFixFilename, const TString& outputFilename, const TString& otag)
: TObject()
, fOutputFilename(outputFilename)
{
//
// constructor
//

	fUnfolding   = new AliLnUnfolding(particle, simuFilename, Form("Unfolding-%s",outputFilename.Data()), otag);
	fFakeTracks  = new AliLnFakeTracks(particle, simuFilename, Form("FakeTracks-%s",outputFilename.Data()), otag);
	fSecondaries = new AliLnSecondaries(particle, dataFilename, simuFilename, Form("Secondaries-%s",outputFilename.Data()), otag);
	fEfficiency  = new AliLnEfficiency(particle, simuFixFilename, Form("Efficiency-%s",outputFilename.Data()), otag); // simufix for efficiencies
}

AliLnCorr::~AliLnCorr()
{
//
// destructor
//
	delete fUnfolding;
	delete fFakeTracks;
	delete fSecondaries;
	delete fEfficiency;
}

Int_t AliLnCorr::Exec()
{
//
// rebuild correction file
//
	fUnfolding->Exec();
	fFakeTracks->Exec();
	fSecondaries->Exec();
	fEfficiency->Exec();
	
	// merge the root files
	
	TString output1 = *fUnfolding->GetOutputFilename();
	TString output2 = *fFakeTracks->GetOutputFilename();
	TString output3 = *fSecondaries->GetOutputFilename();
	TString output4 = *fEfficiency->GetOutputFilename();
	
	TFileMerger m;
	
	m.AddFile(output1.Data(),0);
	m.AddFile(output2.Data(),0);
	m.AddFile(output3.Data(),0);
	m.AddFile(output4.Data(),0);
	
	m.OutputFile(fOutputFilename.Data());
	
	m.Merge();
	
	// remove tmp files
	gSystem->Exec(Form("rm -f %s %s %s %s", output1.Data(), output2.Data(),  output3.Data(),  output4.Data()));
	
	gSystem->Exec(Form("mv debug-%s debug-%s",output3.Data(),fOutputFilename.Data()));
	
	return 0;
}
