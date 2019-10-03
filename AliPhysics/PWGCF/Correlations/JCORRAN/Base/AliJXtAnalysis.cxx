
// ROOT Classes
#include <TChain.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TGraphErrors.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TVector.h>
#include <TVector3.h>
#include <TVectorT.h>
#include <TComplex.h>

// AliRoot classes
#include <AliAnalysisManager.h>
#include <AliInputEventHandler.h>

// Jyvaskyla Classes
#include "AliJXtAnalysis.h"
#include "AliJBaseTrack.h"
#include "AliJHistManager.h"
#include "AliJEfficiency.h"

ClassImp(AliJXtAnalysis)

//________________________________________________________________________
AliJXtAnalysis::AliJXtAnalysis() 
	: AliAnalysisTaskSE(), 
	AnaEntry(0),
	fInputList(0x0),
	fVertex(),
	fCent(0),
	fDebugLevel(0),
	fNCent(0),
	fCBin(0),
	fCentBin(),
	fNJacek(0),
	fPttJacek(),
	fnBinsPt(0),
	fLogBinsPt(),
	fnBinsXt(0),
	fLogBinsXt(),
	fsqrts(7000.),
	fetaRange(0.8),
	fisolR(0.4),
	fzvertexRange(10.0),
	fMinIsolPt(2.0),
	fEfficiency(0x0),
	fEfficiencyIsolated(0x0),
	fEffMode(1),
	fTrackFilterBit(0),
	fEffFilterBit(0),
	fHMG(NULL),
	fHistCentBin(),
	fVertexBin(),
 	fh_cent(),
	fh_ntracks(),
	fh_effCorr(),
 	fh_effCorrIsolated(),
	fh_vertex(),
 	fh_ptJacek(),
 	fh_charPt(),
	fh_invCharPt(),
	fh_charPtNoEff(),
	fh_charEta(),
	fh_charPhi(),
	fh_charXt(),
	fh_invCharXt(),
	fh_isolCharPt(),
	fh_isolInvCharPt(),
	fh_isolCharPtNoEff(),
	fh_isolCharEta(),
	fh_isolCharPhi(),
	fh_isolCharXt(),
	fh_isolInvCharXt(),
	fh_coneActivity(),
	fh_coneActivityIsolated()
{

	const int NCent = 9;
	double CentBin[NCent+1] = {0, 5, 10, 20, 30, 40, 50, 60, 80, 100};
	fNCent = NCent;
	fCentBin = new double[fNCent+1];
	for(int ic=0; ic<=NCent; ic++){
		fCentBin[ic] = CentBin[ic];
	}	

	// pt bins to check pt dist copied from AliJHistos
	const int nJacek = 73;
	double pttJacek[nJacek+1] = {0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.2, 3.4, 3.6, 3.8, 4, 4.5, 5, 5.5, 6, 6.5, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 40, 45, 50, 60, 70, 80, 90, 100};
	
	fNJacek = nJacek;
	fPttJacek = new double[fNJacek+1] ;
	for(int i=0; i<= fNJacek; i++){
		fPttJacek[i] = pttJacek[i];
	}
	
	// === pT binning
	fnBinsPt=150;
	fLogBinsPt = new double[fnBinsPt+1];
	double limL=0.1, limH=100, logBW = (log(limH)-log(limL))/fnBinsPt;
	for(int ij=0;ij<=fnBinsPt;ij++) fLogBinsPt[ij]=limL*exp(ij*logBW);
	
	// === xT binning
	fnBinsXt=200;
	fLogBinsXt = new double[fnBinsXt+1];
	double xTlimL = 1e-5, xTlimH = 0.02, xTlogBW = (log(xTlimH)-log(xTlimL))/fnBinsXt;
	for(int ij=0;ij<=fnBinsXt;ij++) fLogBinsXt[ij]=xTlimL*exp(ij*xTlogBW);

	// Constructor
}
//________________________________________________________________________
AliJXtAnalysis::AliJXtAnalysis(const char *name) 
	: AliAnalysisTaskSE(name), 
	AnaEntry(0),
	fInputList(0x0),
	fVertex(),
	fCent(0),
	fDebugLevel(0),
	fNCent(0),
	fCBin(0),
	fCentBin(),
	fNJacek(0),
	fLogBinsPt(),
	fnBinsXt(0),
	fLogBinsXt(),
	fsqrts(7000.),
	fetaRange(0.8),
	fisolR(0.4),
	fzvertexRange(10.0),
	fMinIsolPt(2.0),
        fEfficiency(0x0),
	fEfficiencyIsolated(0x0),
	fEffMode(1),
	fTrackFilterBit(0),
	fEffFilterBit(0),
	fHMG(NULL),
	fHistCentBin(),
	fVertexBin(),
	fh_cent(),
	fh_ntracks(),
	fh_effCorr(),
	fh_effCorrIsolated(),
	fh_vertex(),
	fh_ptJacek(),
	fh_charPt(),
	fh_invCharPt(),
	fh_charPtNoEff(),
	fh_charEta(),
	fh_charPhi(),
	fh_charXt(),
	fh_invCharXt(),
	fh_isolCharPt(),
	fh_isolInvCharPt(),
	fh_isolCharPtNoEff(),
	fh_isolCharEta(),
	fh_isolCharPhi(),
	fh_isolCharXt(),
	fh_isolInvCharXt(),
	fh_coneActivity(),
	fh_coneActivityIsolated()
{
	const int NCent = 9;
	double CentBin[NCent+1] = {0, 5, 10, 20, 30, 40, 50, 60, 80, 100};
	fNCent = NCent;
	fCentBin = new double[fNCent+1];
	for(int ic=0; ic<=NCent; ic++){
		fCentBin[ic] = CentBin[ic];
	}	


	// pt bins to check pt dist copied from AliJHistos
	const int nJacek = 73;
	double pttJacek[nJacek+1] = {0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.2, 3.4, 3.6, 3.8, 4, 4.5, 5, 5.5, 6, 6.5, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 40, 45, 50, 60, 70, 80, 90, 100};
	
	fNJacek = nJacek;
	fPttJacek = new double[fNJacek+1] ;
	for(int i=0; i<= fNJacek; i++){
		fPttJacek[i] = pttJacek[i];
	}
	
	// === pT binning
	fnBinsPt=150;
	fLogBinsPt = new double[fnBinsPt+1];
	double limL=0.1, limH=100, logBW = (log(limH)-log(limL))/fnBinsPt;
	for(int ij=0;ij<=fnBinsPt;ij++) fLogBinsPt[ij]=limL*exp(ij*logBW);
	
	// === xT binning
	fnBinsXt=200;
	fLogBinsXt = new double[fnBinsXt+1];
	double xTlimL = 1e-5, xTlimH = 0.02, xTlogBW = (log(xTlimH)-log(xTlimL))/fnBinsXt;
	for(int ij=0;ij<=fnBinsXt;ij++) fLogBinsXt[ij]=xTlimL*exp(ij*xTlogBW);
	
	// Constructor
	// Define input and output slots here
}

//________________________________________________________________________
AliJXtAnalysis::AliJXtAnalysis(const AliJXtAnalysis& a):
	AliAnalysisTaskSE(a.GetName()),
	AnaEntry(a.AnaEntry),
	fInputList(a.fInputList),
	fVertex(a.fVertex),
	fCent(a.fCent),
	fDebugLevel(a.fDebugLevel),
	fNCent(a.fNCent),
	fCBin(a.fCBin),
	fCentBin(a.fCentBin),
	fNJacek(a.fNJacek),
	fPttJacek(a.fPttJacek),
	fnBinsPt(a.fnBinsPt),
	fLogBinsPt(a.fLogBinsPt),
	fnBinsXt(a.fnBinsXt),
	fLogBinsXt(a.fLogBinsXt),
	fsqrts(a.fsqrts),
	fetaRange(a.fetaRange),
	fisolR(a.fisolR),
	fzvertexRange(a.fzvertexRange),
	fMinIsolPt(a.fMinIsolPt),
	fEfficiency(a.fEfficiency),
	fEfficiencyIsolated(a.fEfficiencyIsolated),
	fEffMode(a.fEffMode),
	fTrackFilterBit(a.fTrackFilterBit),
	fEffFilterBit(a.fEffFilterBit),
	fHMG(a.fHMG),
	fHistCentBin(a.fHistCentBin),
	fVertexBin(a.fVertexBin),
	fh_cent(a.fh_cent),
	fh_ntracks(a.fh_ntracks),
	fh_effCorr(a.fh_effCorr),
	fh_effCorrIsolated(a.fh_effCorrIsolated),
	fh_vertex(a.fh_vertex),
	fh_ptJacek(a.fh_ptJacek),
	fh_charPt(a.fh_charPt),
	fh_invCharPt(a.fh_invCharPt),
	fh_charPtNoEff(a.fh_charPtNoEff),
	fh_charEta(a.fh_charEta),
	fh_charPhi(a.fh_charPhi),
	fh_charXt(a.fh_charXt),
	fh_invCharXt(a.fh_invCharXt),
	fh_isolCharPt(a.fh_isolCharPt),
	fh_isolInvCharPt(a.fh_isolInvCharPt),
	fh_isolCharPtNoEff(a.fh_isolCharPtNoEff),
	fh_isolCharEta(a.fh_isolCharEta),
	fh_isolCharPhi(a.fh_isolCharPhi),
	fh_isolCharXt(a.fh_isolCharXt),
	fh_isolInvCharXt(a.fh_isolInvCharXt),
	fh_coneActivity(a.fh_coneActivity),
	fh_coneActivityIsolated(a.fh_coneActivityIsolated)
{
	//copy constructor
}
//________________________________________________________________________
AliJXtAnalysis& AliJXtAnalysis::operator = (const AliJXtAnalysis& ap){
	// assignment operator
	this->~AliJXtAnalysis();
	new(this) AliJXtAnalysis(ap);
	return *this;
}
//________________________________________________________________________
void AliJXtAnalysis::Init(){
}
//________________________________________________________________________
void AliJXtAnalysis::UserCreateOutputObjects(){
	
	// Create histograms
	// Called once
	AnaEntry = 0;
	// need to fill to book a histo
	
	fHMG = new AliJHistManager("AliJXtHistManager");
	
	cout << "DEBUG (sami):" << fNCent << endl;
	for(int i = 0; i < 200; i++) cout << "   " << fLogBinsXt[i];
	cout << "   " << endl;
		
	// set AliJBin here // 
	fHistCentBin .Set("CentBin","CentBin","Cent:%d",AliJBin::kSingle).SetBin(fNCent);
	fVertexBin .Set("Vtx","Vtx","Vtx:%d", AliJBin::kSingle).SetBin(3);

	// set AliJTH1D here //
	fh_cent
		<< TH1D("h_cent","h_cent", 100, 0, 100) 
		<< "END" ;

	fh_vertex
		<< TH1D("h_vertex","h_vertex", 400, -20, 20)
		<< fVertexBin
		<< "END" ;
	
	fh_ntracks 
		<< TH1D("h_tracks", "h_tracks", 1000, 0, 30000)
		<< fHistCentBin
		<< "END" ;

	fh_ptJacek
		<< TH1D("hChargedPtJacek", "", fNJacek, fPttJacek)
		<< fHistCentBin 
		<< "END" ; 

	fh_charPt
		<< TH1D("hChargedPt", "hChargedPt", fnBinsPt, fLogBinsPt)
		<< fHistCentBin
		<< "END" ;
	
	fh_invCharPt
		<< TH1D("hInvariantChargedPt", "hInvariantChargedPt", fnBinsPt, fLogBinsPt)
		<< fHistCentBin
		<< "END" ;
	
	fh_charPtNoEff
		<< TH1D("hChargedPtNoEff", "hChargedPtNoEff", fnBinsPt, fLogBinsPt)
		<< fHistCentBin
		<< "END" ;
	
	fh_charEta
		<< TH1D("hChargedEta", "hChargedEta", 100, -1.0, 1.0 )
		<< fHistCentBin
		<< "END" ;
	
	fh_charPhi
		<< TH1D("hChargedPhi", "hChargedPhi", 128, -3.2, 3.2 )
		<< fHistCentBin
		<< "END" ;
	
	fh_charXt
		<< TH1D("hChargedXt", "hChargetXt", fnBinsXt, fLogBinsXt )
		<< fHistCentBin
		<< "END" ;
	
	fh_invCharXt
		<< TH1D("hInvariantChargedXt", "hInvariantChargedXt", fnBinsXt, fLogBinsXt )
		<< fHistCentBin
		<< "END" ;
	
	fh_isolCharPt
		<< TH1D("hIsolatedChargedPt", "hIsolatedChargedPt", fnBinsPt, fLogBinsPt)
		<< fHistCentBin
		<< "END" ;
	
	fh_isolInvCharPt
		<< TH1D("hIsolatedInvariantChargedPt", "hIsolatedInvariantChargedPt", fnBinsPt, fLogBinsPt)
		<< fHistCentBin
		<< "END" ;
	
	fh_isolCharPtNoEff
		<< TH1D("hIsolatedChargedPtNoEff", "hIsolatedChargedPtNoEff", fnBinsPt, fLogBinsPt)
		<< fHistCentBin
		<< "END" ;
	
	fh_isolCharEta
		<< TH1D("hIsolatedChargedEta", "hIsolatedChargedEta", 100, -1.0, 1.0 )
		<< fHistCentBin
		<< "END" ;
	
	fh_isolCharPhi
		<< TH1D("hIsolatedChargedPhi", "hIsolatedChargedPhi", 128, -3.2, 3.2 )
		<< fHistCentBin
		<< "END" ;
	
	fh_isolCharXt
		<< TH1D("hIsolatedChargedXt", "hIsolatedChargetXt", fnBinsXt, fLogBinsXt )
		<< fHistCentBin
		<< "END" ;
	
	fh_isolInvCharXt
		<< TH1D("hIsolatedInvariantChargedXt", "hIsolatedInvariantChargedXt", fnBinsXt, fLogBinsXt )
		<< fHistCentBin
		<< "END" ;
	
	//AliJTH1D set done.
	fh_effCorr
		<< TProfile("pEffCorr", "pEffCorr", fnBinsPt, fLogBinsPt)
		<< "END" ;
	
	fh_effCorrIsolated
		<< TProfile("pEffCorrIsolated", "pEffCorrIsolated", fnBinsPt, fLogBinsPt)
		<< "END" ;
	
	fh_coneActivity
		<< TProfile("pConeActivity", "pConeActivity", fnBinsPt, fLogBinsPt)
		<< "END" ;
	
	fh_coneActivityIsolated
		<< TProfile("pConeActivityIsolated", "pConeActivityIsolated", fnBinsPt, fLogBinsPt)
		<< "END" ;
	
	fHMG->Print();
	fHMG->WriteConfig();

}
//________________________________________________________________________
AliJXtAnalysis::~AliJXtAnalysis() {
	delete fInputList;
	delete fHMG;
	delete []fPttJacek;
	delete []fCentBin;
}

//________________________________________________________________________
void AliJXtAnalysis::UserExec(Option_t *) {
	// Main loop
	// Called for each event
	// DO ANALYSIS WORK HERE // 
	
	// tighter z-vertex cut
	double zVert = fVertex[2];	
        if( zVert < -fzvertexRange || zVert > fzvertexRange ) return;
	
	// QA-histos for vertex, accepted events
	for(int iaxis=0; iaxis<3; iaxis++){
		fh_vertex[iaxis]->Fill( fVertex[iaxis] );
	}

	// find Centrality
	double inputCent = fCent;
	
	for(int iCbin=0; iCbin<fNCent; iCbin++){
		if( inputCent >= fCentBin[iCbin] && inputCent < fCentBin[iCbin+1] ){fCBin = iCbin;};
	}
	if(fCBin == -1){return;};
	
	int trk_number = fInputList->GetEntriesFast();
	fh_ntracks[fCBin]->Fill( trk_number ) ;
	fh_cent->Fill( inputCent ) ; 
	
         Long64_t ntracks = fInputList->GetEntriesFast();
         for( Long64_t it=0; it<ntracks; it++){
		// load track
		AliJBaseTrack *track = (AliJBaseTrack*)fInputList->At(it);
		
		// eta cut
		double eta = track->Eta();
		if( eta < -fetaRange || eta > fetaRange ) continue;
			
		// Basic variables
		double pT = track->Pt();
		double xT = 2.*pT/fsqrts;
		double phi = track->Phi();
		
		// Efficiency correction
		double effCorr = 1./fEfficiency->GetCorrection(pT, fEffFilterBit, inputCent);
		FillEfficiencyCheckHistogram(pT,effCorr);
		
		// Fill controll histograms. Note: still in full ALICE acceptance
		FillControlHistograms(pT, effCorr, fCBin);
		
		// Restric to fiducial acceptance and fill inclusive histograms
		if( eta < -fetaRange + fisolR || eta > fetaRange - fisolR ) continue; // Fiducial eta cut
		FillInclusiveHistograms(pT, xT, eta, phi, effCorr, fCBin);
		
		// BEGIN ISOLATION:
		if( pT < fMinIsolPt ) continue; // minimum pT cut
		
		double sum = 0.0;
		TLorentzVector lvTrack = TLorentzVector( track->Px(), track->Py(), track->Pz(), track->P() );
		for( Long64_t ia=0; ia<ntracks; ia++){
			if (ia == it ) continue; // do not count the track itself
			AliJBaseTrack *assocTrack = (AliJBaseTrack*)fInputList->At(ia);
			TLorentzVector lvAssoc = TLorentzVector( assocTrack->Px(), assocTrack->Py(), assocTrack->Pz(), assocTrack->P() );
			if( lvTrack.DeltaR( lvAssoc ) < fisolR ) sum += lvAssoc.Pt();
		}
		
		// Cone activity
		FillInclusiveConeActivities(pT,sum);
		
		// TODO: final criteria will come from isolation efficiency class
		double isolThreshold = 0.1 * pT;
		
		if( sum < isolThreshold ){
			if( fDebugLevel > 0 ){
				cout << "DEBUG: isolated particle found " << endl;
				cout << "fCentBin = " << fCentBin << ", pT = " << pT << ", xT = " << xT << " and eta = " << eta << endl;
			}
			
			FillIsolatedConeActivities(pT, sum);

			// isolated efficiency
			double effCorr = 1./fEfficiencyIsolated->GetCorrection(pT, fEffFilterBit, inputCent);
			FillIsolatedEfficiencyCheckHistogram(pT, effCorr);
			
			// Fill histos
			FillIsolatedHistograms(pT, xT, eta, phi, effCorr, fCBin);
		}
         }

	AnaEntry++;
	return;
}
//________________________________________________________________________
void AliJXtAnalysis::FillEfficiencyCheckHistogram(double pT, double effCorr){
	// to check the efficiency
	fh_effCorr->Fill(pT,1./effCorr);
}
//________________________________________________________________________
void AliJXtAnalysis::FillIsolatedEfficiencyCheckHistogram(double pT, double effCorr){
	// to check the efficiency
	fh_effCorrIsolated->Fill(pT,1./effCorr);
}
//________________________________________________________________________
void AliJXtAnalysis::FillControlHistograms(double pT, double effCorr, int centBin){
	// To compare with the published data, here the same binning as compared published
	fh_ptJacek[centBin]->Fill(pT, effCorr);
}
//________________________________________________________________________
void AliJXtAnalysis::FillInclusiveHistograms(double pT, double xT, double eta, double phi, double effCorr, int centBin){
	// Fill all inclusive histograms
	fh_charEta[centBin]->Fill(eta, effCorr);
	fh_charPhi[centBin]->Fill(phi, effCorr);
	fh_charPt[centBin]->Fill(pT, effCorr);
	fh_charPtNoEff[centBin]->Fill(pT);
	fh_charXt[centBin]->Fill(xT, effCorr);
	fh_invCharPt[centBin]->Fill(pT, pT > 0. ? effCorr/pT : 0.); // Can you do 1/pT here???
	fh_invCharXt[centBin]->Fill(xT, xT > 0. ? effCorr/xT : 0.); // Can you do 1/xT here???
}
//______________________________________________________________________________
void AliJXtAnalysis::FillIsolatedHistograms(double pT, double xT, double eta, double phi, double effCorr, int centBin){
	// Fill all inclusive histograms
	fh_isolCharEta[centBin]->Fill(eta, effCorr);
	fh_isolCharPhi[centBin]->Fill(phi, effCorr);
	fh_isolCharPt[centBin]->Fill(pT, effCorr);
	fh_isolCharPtNoEff[centBin]->Fill(pT);
	fh_isolCharXt[centBin]->Fill(xT, effCorr);
	fh_isolInvCharPt[centBin]->Fill(pT, pT > 0. ? effCorr/pT : 0.); // Can you do 1/pT here???
	fh_isolInvCharXt[centBin]->Fill(xT, xT > 0. ? effCorr/xT : 0.); // Can you do 1/xT here???
}
//______________________________________________________________________________
void AliJXtAnalysis::FillInclusiveConeActivities(double pT, double sumPt){
	fh_coneActivity->Fill(pT,sumPt);
}
//______________________________________________________________________________
void AliJXtAnalysis::FillIsolatedConeActivities(double pT, double sumPt){
	fh_coneActivityIsolated->Fill(pT,sumPt);
}
//______________________________________________________________________________
void AliJXtAnalysis::Terminate(Option_t *) 
{
	cout<<"Successfully Finished"<<endl;
}
//________________________________________________________________________
;
