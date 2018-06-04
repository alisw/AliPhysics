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

#include "AliJCDijetHistos.h"
#include <TGrid.h>
#include <TPRegexp.h>

Double_t AliJCDijetHistos::CentBin[NCENT+1] = {0, 5, 10, 20, 30, 40, 50, 60};
UInt_t AliJCDijetHistos::NCentBin = sizeof(AliJCDijetHistos::CentBin)/sizeof(AliJCDijetHistos::CentBin[0])-1;
Double_t AliJCDijetHistos::pttJacek[74+16] = {0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95,1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.2, 3.4, 3.6, 3.8, 4, 4.5, 5, 5.5, 6, 6.5, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 40, 45, 50, 60, 70, 80, 90, 100, 110, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 340, 380, 420, 460, 500};
UInt_t AliJCDijetHistos::NpttJacek = sizeof(AliJCDijetHistos::pttJacek)/sizeof(AliJCDijetHistos::pttJacek[0])-1;

//______________________________________________________________________________
AliJCDijetHistos::AliJCDijetHistos() :
	fHMG(NULL),
	fHistCentBin(),
	fh_pt(),
	fh_eta(),
	fh_phi()
{
	
}

//______________________________________________________________________________
AliJCDijetHistos::AliJCDijetHistos(const AliJCDijetHistos& obj) :
	fHMG(obj.fHMG),
	fHistCentBin(obj.fHistCentBin),
	fh_pt(obj.fh_pt),
	fh_eta(obj.fh_eta),
	fh_phi(obj.fh_phi)
{
	// copy constructor
}

//______________________________________________________________________________
AliJCDijetHistos& AliJCDijetHistos::operator=(const AliJCDijetHistos& obj){
	// copy constructor
	return *this;
}

//______________________________________________________________________________
AliJCDijetHistos::~AliJCDijetHistos() {
	// destructor
	delete fHMG;
}


//______________________________________________________________________________
void AliJCDijetHistos::CreateEventTrackHistos(){
	// Create basic event histograms
	fHMG = new AliJHistManager("AliJCDijetHistManager","jcdijet");
	// set AliJBin here //
	fHistCentBin .Set("CentBin","CentBin","Cent:%d",AliJBin::kSingle).SetBin(NCentBin);

    int NBINS=170;
    double LogBinsX[NBINS+1], LimL=0.1, LimH=1000;
    double logBW = (log(LimH)-log(LimL))/NBINS;
    for(int ij=0;ij<=NBINS;ij++) LogBinsX[ij]=LimL*exp(ij*logBW);

    // ============= CHARGED PARTICLE HISTOS ============= 
	fh_pt
		<< TH1D("hPtJacek", "hPtJacek", AliJCDijetHistos::NpttJacek, AliJCDijetHistos::pttJacek)
        //<< TH1D("hChargedPt","Charged particle pt",NBINS, LogBinsX )
		<< fHistCentBin
		<< "END" ;

	fh_eta
		<< TH1D("h_eta", "h_eta", 300, -10, 10 )
		<< fHistCentBin
		<< "END" ;

	fh_phi
		<< TH1D("h_phi", "h_phi", 100, -10, 10)
		<< fHistCentBin
		<< "END" ;

    // ============= JET HISTOS ============= 
	fh_jetPt
		<< TH1D("hJetPtJacek", "hJetPtJacek", AliJCDijetHistos::NpttJacek, AliJCDijetHistos::pttJacek)
        //<< TH1D("hJetPt","Charged particle pt",NBINS, LogBinsX )
		<< fHistCentBin
		<< "END" ;

	fh_jetEta
		<< TH1D("h_jetEta", "h_jetEta", 300, -10, 10 )
		<< fHistCentBin
		<< "END" ;

	fh_jetPhi
		<< TH1D("h_jetPhi", "h_jetPhi", 100, -10, 10)
		<< fHistCentBin
		<< "END" ;

    // ============= DIJET HISTOS ============= 
	fh_DijetInvM
		<< TH1D("h_DijetInvM", "h_DijetInvM", NBINS, LogBinsX)
		<< fHistCentBin
		<< "END" ;

	fh_DijetPtPair
		<< TH1D("h_DijetPtPair", "h_DijetPtPair", AliJCDijetHistos::NpttJacek, AliJCDijetHistos::pttJacek )
		<< fHistCentBin
		<< "END" ;

	fh_DijetDeltaPhi
		<< TH1D("h_DijetDeltaPhi", "h_DijetDeltaPhi", 100, 0, 10)
		<< fHistCentBin
		<< "END" ;

	fh_DijetInvMDeltaPhiCut
		<< TH1D("h_DijetInvMDeltaPhiCut", "h_DijetInvMDeltaPhiCut", NBINS, LogBinsX)
		<< fHistCentBin
		<< "END" ;
}

int AliJCDijetHistos::GetCentralityClass(Double_t fCent){
	for(UInt_t iCbin = 0; iCbin < NCentBin; iCbin++){
		if(fCent > CentBin[iCbin] && fCent < CentBin[iCbin+1])
			return iCbin;
	}
	return -1;
}

