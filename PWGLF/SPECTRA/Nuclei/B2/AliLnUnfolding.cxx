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

// unfolding correction for the pt
// author: Eulogio Serradilla <eulogio.serradilla@cern.ch>

//#include <Riostream.h>
#include <TFile.h>
#include <TSystem.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TString.h>
#include <TF1.h>

#include "AliLnUnfolding.h"
#include "B2.h"

ClassImp(AliLnUnfolding)

AliLnUnfolding::AliLnUnfolding(const TString& particle, const TString& simuFilename, const TString& outputFilename, const TString& otag)
: TObject()
, fParticle(particle)
, fSimuFilename(simuFilename)
, fOutputFilename(outputFilename)
, fOutputTag(otag)
{
//
// constructor
//
}

AliLnUnfolding::~AliLnUnfolding()
{
//
// destructor
//
}

Int_t AliLnUnfolding::Exec()
{
//
// extract the detector response from the simulation for particle/antiparticle
//
	using namespace std;
	
	TFile* fsimu = new TFile(fSimuFilename.Data(), "read");
	if(fsimu->IsZombie()) exit(1);
	
	TFile* foutput = new TFile(fOutputFilename.Data(), "recreate");
	if(fOutputTag != "")
	{
		foutput->mkdir(fOutputTag.Data());
		foutput->cd(fOutputTag.Data());
	}
	
	TH2D* hResponseMtx = (TH2D*)FindObj(fsimu, fParticle + "_Prim_Response_Matrix");
	hResponseMtx->SetName(Form("%s_Response_Matrix",fParticle.Data()));
	hResponseMtx->Write();
	
	TH1D* hMeasuredPt = (TH1D*)FindObj(fsimu, fParticle + "_PID_Pt");
	hMeasuredPt->SetName(Form("%s_Measured_Pt", fParticle.Data()));
	hMeasuredPt->Reset();
	hMeasuredPt->Write();
	
	TH1D* hTruePt = (TH1D*)FindObj(fsimu, fParticle + "_Sim_PID_Prim_Pt");
	hTruePt->SetName(Form("%s_True_Pt", fParticle.Data()));
	hTruePt->Reset();
	hTruePt->Write();
	
	delete foutput;
	delete fsimu;
	
	return 0;
}
