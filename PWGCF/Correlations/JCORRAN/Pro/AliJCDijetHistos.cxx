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
		<< TH1D("h_ptJacek", "h_ptJacek", AliJCDijetHistos::NpttJacek, AliJCDijetHistos::pttJacek)
        //<< TH1D("hChargedPt","Charged particle pt",NBINS, LogBinsX )
		<< fHistCentBin
		<< "END" ;

	fh_eta
		<< TH1D("h_eta", "h_eta", 100, -1.0, 1.0 )
		<< fHistCentBin
		<< "END" ;

	fh_phi
		<< TH1D("h_phi", "h_phi", 100, -TMath::Pi(), TMath::Pi())
		<< fHistCentBin
		<< "END" ;

	fh_etaPhi
		<< TH2D("h_etaPhi", "h_etaPhi", 100, -1.0, 1.0, 100, -TMath::Pi(), TMath::Pi())
		<< fHistCentBin
		<< "END" ;

    // ============= JET HISTOS ============= 
	fh_jetPt
		<< TH1D("h_jetPtJacek", "h_jetPtJacek", AliJCDijetHistos::NpttJacek, AliJCDijetHistos::pttJacek)
        //<< TH1D("hJetPt","Charged particle pt",NBINS, LogBinsX )
		<< fHistCentBin
		<< "END" ;

	fh_jetEta
		<< TH1D("h_jetEta", "h_jetEta", 100, -1.0, 1.0)
		<< fHistCentBin
		<< "END" ;

	fh_jetPhi
		<< TH1D("h_jetPhi", "h_jetPhi", 100, -TMath::Pi(), TMath::Pi())
		<< fHistCentBin
		<< "END" ;

	fh_jetEtaPhi
		<< TH2D("h_jetEtaPhi", "h_jetEtaPhi", 100, -1.0, 1.0, 100, -TMath::Pi(), TMath::Pi())
		<< fHistCentBin
		<< "END" ;

	fh_rho
		<< TH1D("h_rho", "h_rho", NBINS, LogBinsX)
		<< fHistCentBin
		<< "END" ;

	fh_rhom
		<< TH1D("h_rhom", "h_rhom", NBINS, LogBinsX)
		<< fHistCentBin
		<< "END" ;

	fh_jetArea
		<< TH1D("h_jetArea", "h_jetArea", AliJCDijetHistos::NpttJacek, AliJCDijetHistos::pttJacek)
		<< fHistCentBin
		<< "END" ;

	fh_jetAreaRho
		<< TH1D("h_jetAreaRho", "h_jetAreaRho", AliJCDijetHistos::NpttJacek, AliJCDijetHistos::pttJacek)
		<< fHistCentBin
		<< "END" ;

	fh_corrJetPt
		<< TH1D("h_corrJetPtJacek", "h_corrJetPtJacek", AliJCDijetHistos::NpttJacek, AliJCDijetHistos::pttJacek)
		<< fHistCentBin
		<< "END" ;

	fh_corrJetEta
		<< TH1D("h_corrJetEta", "h_corrJetEta", 100, -1.0, 1.0 )
		<< fHistCentBin
		<< "END" ;

	fh_corrJetPhi
		<< TH1D("h_corrJetPhi", "h_corrJetPhi", 100, -TMath::Pi(), TMath::Pi())
		<< fHistCentBin
		<< "END" ;

    // ============= DIJET HISTOS ============= 
	fh_dijetInvM
		<< TH1D("h_dijetInvM", "h_dijetInvM", NBINS, LogBinsX)
		<< fHistCentBin
		<< "END" ;

	fh_dijetPtPair
		<< TH1D("h_dijetPtPair", "h_dijetPtPair", AliJCDijetHistos::NpttJacek, AliJCDijetHistos::pttJacek )
		<< fHistCentBin
		<< "END" ;

	fh_dijetDeltaPhi
		<< TH1D("h_dijetDeltaPhi", "h_dijetDeltaPhi", 100, 0, 10)
		<< fHistCentBin
		<< "END" ;

	fh_dijetInvMDeltaPhiCut
		<< TH1D("h_dijetInvMDeltaPhiCut", "h_dijetInvMDeltaPhiCut", NBINS, LogBinsX)
		<< fHistCentBin
		<< "END" ;

	fh_corrDijetInvM
		<< TH1D("h_corrDijetInvM", "h_corrDijetInvM", NBINS, LogBinsX)
		<< fHistCentBin
		<< "END" ;

	fh_corrDijetPtPair
		<< TH1D("h_corrDijetPtPair", "h_corrDijetPtPair", AliJCDijetHistos::NpttJacek, AliJCDijetHistos::pttJacek )
		<< fHistCentBin
		<< "END" ;

	fh_corrDijetDeltaPhi
		<< TH1D("h_corrDijetDeltaPhi", "h_corrDijetDeltaPhi", 100, 0, 10)
		<< fHistCentBin
		<< "END" ;

	fh_corrDijetInvMDeltaPhiCut
		<< TH1D("h_corrDijetInvMDeltaPhiCut", "h_corrDijetInvMDeltaPhiCut", NBINS, LogBinsX)
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

