/**************************************************************************
 * Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
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

// Container class for histograms needed in the analysis.

#include "AliJFlowHistos.h"
#include "AliJFlowBaseTask.h"
#include <TGrid.h>
#include <TPRegexp.h>

//Double_t AliJFlowHistos::CentBin[NCENT+1] = {0, 5, 10, 20, 30, 40, 50, 60};
Double_t AliJFlowHistos::CentBin[NCENT+1] = {0,0.001,0.01,0.1,0.5,1,2,3,4,5, 10, 20, 30, 40, 50, 60};
UInt_t AliJFlowHistos::NCentBin = sizeof(AliJFlowHistos::CentBin)/sizeof(AliJFlowHistos::CentBin[0])-1;
Double_t AliJFlowHistos::pttJacek[74] = {0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95,1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.2, 3.4, 3.6, 3.8, 4, 4.5, 5, 5.5, 6, 6.5, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 40, 45, 50, 60, 70, 80, 90, 100};
UInt_t AliJFlowHistos::NpttJacek = sizeof(AliJFlowHistos::pttJacek)/sizeof(AliJFlowHistos::pttJacek[0])-1;

//______________________________________________________________________________
AliJFlowHistos::AliJFlowHistos() :
	fHMG(NULL),
	fHistCentBin(),
	fBin_DetSet(),
	fBin_hh(),
	fh_pt(),
	fh_eta(),
	fh_phi(),
	fh_EP(),
	fhEPCorrInHar()
{
	
}

//______________________________________________________________________________
AliJFlowHistos::AliJFlowHistos(const AliJFlowHistos& obj) :
	fHMG(obj.fHMG),
	fHistCentBin(obj.fHistCentBin),
	fBin_DetSet(obj.fBin_DetSet),
	fBin_hh(obj.fBin_hh),
	fh_pt(obj.fh_pt),
	fh_eta(obj.fh_eta),
	fh_phi(obj.fh_phi),
	fh_EP(obj.fh_EP),
	fhEPCorrInHar(obj.fhEPCorrInHar)
{
	// copy constructor
}

//______________________________________________________________________________
AliJFlowHistos& AliJFlowHistos::operator=(const AliJFlowHistos& obj){
	// copy constructor
	return *this;
}

//______________________________________________________________________________
AliJFlowHistos::~AliJFlowHistos() {
	// destructor
	delete fHMG;
}


//______________________________________________________________________________
void AliJFlowHistos::CreateEventTrackHistos(){
	// Create basic event histograms
	fHMG = new AliJHistManager("AliJFlowHistManager","jflow");
	// set AliJBin here //
	fHistCentBin .Set("CentBin","CentBin","Cent:%d",AliJBin::kSingle).SetBin(NCentBin);
	fBin_DetSet .Set("DetSet","DetSet","DetSet:%d", AliJBin::kSingle).SetBin(AliJFlowBaseTask::D_COUNT);
	fBin_hh .Set("NHH","NHH","NHH:%d", AliJBin::kSingle).SetBin(2); // n=2 and 3

	fh_pt
		<< TH1D("hChargedPtJacek", "", AliJFlowHistos::NpttJacek, AliJFlowHistos::pttJacek)
		<< fHistCentBin << fBin_DetSet
		<< "END" ;

	fh_eta
		<< TH1D("h_eta", "h_eta", 300, -10, 10 )
		<< fHistCentBin << fBin_DetSet
		<< "END" ;
	fh_phi
		<< TH1D("h_phi", "h_phi", 100, -10, 10)
		<< fHistCentBin << fBin_DetSet
		<< "END" ;
	fh_EP
		<< TH1D("h_EP", "h_EP", 100, -10, 10)
		<< fHistCentBin << fBin_DetSet
		<< fBin_hh
		<< "END" ;
	fhEPCorrInHar
		<< TH1D("hEPCorrInHar", "hEPCorrInHar", 500, -TMath::Pi()*2, TMath::Pi()*2)
		<< fHistCentBin << fBin_DetSet
		<< fBin_hh
		<< "END" ;
	//fhEPCorr2D
	//	<< TH2D( "hEPCorr2D", "", 100,-TMath::Pi()/2, TMath::Pi()/2, 1000, -TMath::Pi()/2, TMath::Pi()/2)
	//	<< fHistCentBin << fBin_DetSet
	//	<< fBin_hh
	//	<< "END" ;

}

int AliJFlowHistos::GetCentralityClass(Double_t fCent){
	for(UInt_t iCbin = 0; iCbin < NCentBin; iCbin++){
		if(fCent > CentBin[iCbin] && fCent < CentBin[iCbin+1])
			return iCbin;
	}
	return -1;
}

