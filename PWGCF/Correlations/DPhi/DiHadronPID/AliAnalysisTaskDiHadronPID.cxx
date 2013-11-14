/************************************************************************* 
* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. * 
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

// -----------------------------------------------------------------------
//  This analysis task fills histograms with PID information of tracks 
//  associated to a high p_T trigger.
// -----------------------------------------------------------------------
//  Author: Misha Veldhoen (misha.veldhoen@cern.ch)

#include <iostream>

// Basic Includes
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "THn.h"
#include "TFile.h"
#include "TChain.h"
#include "TObject.h"
#include "TObjArray.h"

// Manager/ Handler
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"

// Event pool includes.
#include "AliEventPoolManager.h"

// PID includes.
#include "AliPIDResponse.h"

// AOD includes.
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAODHandler.h"
#include "AliAODVertex.h"
#include "AliAODInputHandler.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"

// Additional includes.
#include "AliTrackDiHadronPID.h"
#include "AliAODTrackCutsDiHadronPID.h"
#include "AliAODEventCutsDiHadronPID.h"
#include "AliHistToolsDiHadronPID.h"
#include "AliFunctionsDiHadronPID.h"

// AnalysisTask headers.
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskDiHadronPID.h"

using namespace std;

ClassImp(AliAnalysisTaskDiHadronPID);

// -----------------------------------------------------------------------
AliAnalysisTaskDiHadronPID::AliAnalysisTaskDiHadronPID():
	AliAnalysisTaskSE(),
	fPIDResponse(0x0),
	fEventCuts(0x0),
	fTrackCutsTrigger(0x0),
	fTrackCutsAssociated(0x0),
	fPoolMgr(0x0),
	fTriggerTracks(0x0),
	fAssociatedTracks(0x0),
	fCurrentAODEvent(0x0),
	fOutputList(0x0),
	fPtSpectrumTOFbins(0x0),
	fCorrelationsTOFbins(0x0),
	fMixedEventsTOFbins(0x0),
	fPtSpectrumTOFTPCbins(0x0),
	fCorrelationsTOFTPCbins(0x0),
	fMixedEventsTOFTPCbins(0x0),
	fMixedEventsTOFTPCbinsPID(0x0),	
	fTOFhistos(0x0),
	fTOFmismatch(0x0),
	fTOFPtAxis(0x0),
	fTOFTPChistos(0x0),
	fTOFTPCmismatch(0x0),
	fTOFTPCPtAxis(0x0),
	fNDEtaBins(32),
	fNDPhiBins(32),	
	fMinNEventsForMixing(5),
	fPoolTrackDepth(2000),
	fPoolSize(1000),
	fMixEvents(kTRUE),
	fMixTriggers(kFALSE),
	fCalculateMismatch(kTRUE),
	fT0Fill(0x0),
	fLvsEta(0x0),
	fLvsEtaProjections(0x0),	
	fMakeTOFcorrelations(kTRUE),
	fMakeTOFTPCcorrelationsPi(kFALSE),
	fMakeTOFTPCcorrelationsKa(kFALSE),
	fMakeTOFTPCcorrelationsPr(kFALSE),	
	fTOFIntervalFactorTOFTPC(1.),
	fExtendPtAxis(kFALSE)

{

	//
	// Default Constructor.
	//

	if (fDebug > 0) {AliInfo("AliAnalysisTaskDiHadronPID Default Constructor.");}		

}

// -----------------------------------------------------------------------
AliAnalysisTaskDiHadronPID::AliAnalysisTaskDiHadronPID(const char* name):
	AliAnalysisTaskSE(name),
	fPIDResponse(0x0),
	fEventCuts(0x0),
	fTrackCutsTrigger(0x0),
	fTrackCutsAssociated(0x0),
	fPoolMgr(0x0),
	fTriggerTracks(0x0),
	fAssociatedTracks(0x0),	
	fCurrentAODEvent(0x0),
	fOutputList(0x0),
	fPtSpectrumTOFbins(0x0),
	fCorrelationsTOFbins(0x0),
	fMixedEventsTOFbins(0x0),
	fPtSpectrumTOFTPCbins(0x0),
	fCorrelationsTOFTPCbins(0x0),
	fMixedEventsTOFTPCbins(0x0),
	fMixedEventsTOFTPCbinsPID(0x0),			
	fTOFhistos(0x0),
	fTOFmismatch(0x0),
	fTOFPtAxis(0x0),
	fTOFTPChistos(0x0),
	fTOFTPCmismatch(0x0),
	fTOFTPCPtAxis(0x0),	
	fNDEtaBins(32),
	fNDPhiBins(32),
	fMinNEventsForMixing(5),
	fPoolTrackDepth(2000),
	fPoolSize(1000),
	fMixEvents(kTRUE),	
	fMixTriggers(kFALSE),
	fCalculateMismatch(kTRUE),
	fT0Fill(0x0),
	fLvsEta(0x0),
	fLvsEtaProjections(0x0),
	fMakeTOFcorrelations(kTRUE),
	fMakeTOFTPCcorrelationsPi(kFALSE),
	fMakeTOFTPCcorrelationsKa(kFALSE),
	fMakeTOFTPCcorrelationsPr(kFALSE),	
	fTOFIntervalFactorTOFTPC(1.),
	fExtendPtAxis(kFALSE)

{

	//
	// Named Constructor.
	//

	if (fDebug > 0) {cout << Form("File: %s, Line: %i, Function: %s",__FILE__,__LINE__,__func__) << endl;}

	DefineInput(0,TChain::Class());
	DefineOutput(1,TList::Class());

}

// -----------------------------------------------------------------------
AliAnalysisTaskDiHadronPID::~AliAnalysisTaskDiHadronPID() {

	//
	// Destructor.
	//

	if (fDebug > 0) {cout << Form("File: %s, Line: %i, Function: %s",__FILE__,__LINE__,__func__) << endl;}

	if (fPoolMgr) {delete fPoolMgr; fPoolMgr = 0x0;}
	if (fOutputList) {delete fOutputList; fOutputList = 0x0;}

}

// -----------------------------------------------------------------------
void AliAnalysisTaskDiHadronPID::UserCreateOutputObjects() {

	//
	// Create Output objects.
	//

	if (fDebug > 0) {cout << Form("File: %s, Line: %i, Function: %s",__FILE__,__LINE__,__func__) << endl;}

	AliAnalysisManager* manager = AliAnalysisManager::GetAnalysisManager();
	if (!manager) {AliFatal("Could not obtain analysis manager.");}	
	AliInputEventHandler* inputHandler = dynamic_cast<AliInputEventHandler*> (manager->GetInputEventHandler());
	if (!inputHandler) {AliFatal("Could not obtain input handler.");}	

	// Getting the pointer to the PID response object.
	fPIDResponse = inputHandler->GetPIDResponse();	
	if (!fPIDResponse) {AliFatal("Could not obtain PID response.");}

	// For now we don't bin in multiplicity for pp.
	TArrayD* centralityBins = 0x0;
	if (fEventCuts->GetIsPbPb()) {
		Double_t tmp[] = {0., 1., 2., 3., 4., 5., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100.1 };
		centralityBins = new TArrayD(15, tmp);
	} else {
		Double_t tmp[] = {0.,1.};
		centralityBins = new TArrayD(2, tmp);
	}

	Int_t nZvtxBins  = 7;
	Double_t vertexBins[] = {-7., -5., -3., -1., 1., 3., 5., 7.};

	fPoolMgr = new AliEventPoolManager(fPoolSize, fPoolTrackDepth, centralityBins->GetSize(), centralityBins->GetArray(), nZvtxBins, (Double_t*) vertexBins);
    
	delete centralityBins;

	// Create the output list.
	fOutputList = new TList();
	fOutputList->SetOwner(kTRUE);

	// Creating all requested histograms locally.
	fEventCuts->CreateHistos();
	fOutputList->Add(fEventCuts);

	fTrackCutsTrigger->CreateHistos();
	fOutputList->Add(fTrackCutsTrigger);

	fTrackCutsAssociated->CreateHistos();
	fOutputList->Add(fTrackCutsAssociated);

	TString speciesname[] = {"Pion","Kaon","Proton"};

	// Create TOF correlations histograms (DPhi,DEta,TOF).
	if (fMakeTOFcorrelations) {

		// Get the pT axis for the TOF PID correlations.
		Double_t* ptaxis = fTrackCutsAssociated->GetPtAxisPID();
		Int_t nptbins = fTrackCutsAssociated->GetNPtBinsPID();

		// Create Pt spectrum histogram.
		fPtSpectrumTOFbins = new TH1F("fPtSpectrumTOFbins","p_{T} Spectrum;p_{T} (GeV/c);Count",nptbins,ptaxis);
		fOutputList->Add(fPtSpectrumTOFbins);

		// Create unidentified correlations histogram.
		fCorrelationsTOFbins = AliHistToolsDiHadronPID::MakeHist3D("fCorrelationsTOFbins","Correlations;#Delta#phi;#Delta#eta;p_{T} (GeV/c)",
			fNDPhiBins,-TMath::Pi()/2.,3.*TMath::Pi()/2.,
			fNDEtaBins,-1.6,1.6,
			nptbins, ptaxis);
		fOutputList->Add(fCorrelationsTOFbins);

		// Create unidentified mixed events histogram.
		fMixedEventsTOFbins = AliHistToolsDiHadronPID::MakeHist3D("fMixedEventsTOFbins","Mixed Events;#Delta#phi;#Delta#eta;p_{T} (GeV/c)",
			fNDPhiBins,-TMath::Pi()/2.,3.*TMath::Pi()/2.,
			fNDEtaBins,-1.6,1.6,
			nptbins, ptaxis);
		fOutputList->Add(fMixedEventsTOFbins);

		// Create TOFPtaxis.
		fTOFPtAxis = new TAxis(nptbins, ptaxis);
		fTOFPtAxis->SetName("fTOFPtAxis");
		fTOFPtAxis->SetTitle("p_{T} GeV/c");

		// Create PID histograms.
		fTOFhistos = new TObjArray(3);
		fTOFhistos->SetOwner(kTRUE);	
		fTOFhistos->SetName("CorrelationsTOF");

		if (fCalculateMismatch) {
			fTOFmismatch = new TObjArray(3);
			fTOFmismatch->SetOwner(kTRUE);
			fTOFmismatch->SetName("MismatchTOF");
		}

		for (Int_t iSpecies = 0; iSpecies < 3; iSpecies++) {

			TObjArray* TOFhistosTmp = new TObjArray(fTOFPtAxis->GetNbins());
			TOFhistosTmp->SetOwner(kTRUE);
			TOFhistosTmp->SetName(speciesname[iSpecies].Data());

			TObjArray* TOFmismatchTmp = 0x0;
			if (fCalculateMismatch) {
				TOFmismatchTmp = new TObjArray(fTOFPtAxis->GetNbins());
				TOFmismatchTmp->SetOwner(kTRUE);
				TOFmismatchTmp->SetName(speciesname[iSpecies].Data());
			}

			for (Int_t iBinPt = 1; iBinPt < (fTOFPtAxis->GetNbins() + 1); iBinPt++) {

				Int_t iPtClass = fTrackCutsAssociated->GetPtClass(iBinPt);
				if (iPtClass == -1) {AliFatal("Not valid pT class."); continue;}

				Int_t NBinsTOF = fTrackCutsAssociated->GetNTOFbins(iPtClass,iSpecies);
				Double_t TOFmin = fTrackCutsAssociated->GetTOFmin(iPtClass,iSpecies);
				Double_t TOFmax = fTrackCutsAssociated->GetTOFmax(iPtClass,iSpecies);

				//cout << "ptbin: "<< iBinPt << " class: " << iPtClass << " TOFBins: " << NBinsTOF << " min: " << TOFmin << " max: " << TOFmax << endl; 

				// Correlation histogram.
				TH3F* htmp = new TH3F(Form("fCorrelationsTOF_%i",iBinPt),
					Form("%5.3f < p_{T} < %5.3f; #Delta#phi; #Delta#eta; t_{TOF} (ps)", fTOFPtAxis->GetBinLowEdge(iBinPt), fTOFPtAxis->GetBinUpEdge(iBinPt)), 
					fNDPhiBins, -TMath::Pi()/2., 3.*TMath::Pi()/2.,
					fNDEtaBins, -1.6, 1.6, NBinsTOF, TOFmin, TOFmax);
				htmp->SetDirectory(0);

				TOFhistosTmp->Add(htmp);

				if (fCalculateMismatch) {
					// Mismatch histogram.
					TH1F* htmp2 = new TH1F(Form("fMismatchTOF_%i",iBinPt),
						Form("%5.3f < p_{T} < %5.3f; t_{TOF} (ps)", fTOFPtAxis->GetBinLowEdge(iBinPt), fTOFPtAxis->GetBinUpEdge(iBinPt)),
						NBinsTOF, TOFmin, TOFmax);
					htmp2->SetDirectory(0);

					TOFmismatchTmp->Add(htmp2);	
				}	
			}

			fTOFhistos->Add(TOFhistosTmp);
			if (fCalculateMismatch) {fTOFmismatch->Add(TOFmismatchTmp);}

		}

		fOutputList->Add(fTOFhistos);
		if (fCalculateMismatch) {fOutputList->Add(fTOFmismatch);}

	}

	// Create TOF/TPC correlation histograms. (DPhi,DEta,TOF,TPC).
	if (fMakeTOFTPCcorrelationsPi || fMakeTOFTPCcorrelationsKa || fMakeTOFTPCcorrelationsPr) {

		Double_t ptarrayTOFTPC[16] = {2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 
									  2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 
									  4.2, 4.6, 5.0};
	  	const Int_t nptbins = 15;
		Double_t ptarrayTOFTPCext[26] = {1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 
									  2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 
									  2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 
									  4.2, 4.6, 5.0};
		const Int_t nptbinsext = 25;

		fTOFTPCPtAxis = new TAxis(fExtendPtAxis ? nptbinsext : nptbins, fExtendPtAxis ? ptarrayTOFTPCext : ptarrayTOFTPC);
		fTOFTPCPtAxis->SetName("fTOFTPCPtAxis");
		fTOFTPCPtAxis->SetTitle("p_{T} GeV/c");

		// Create Pt spectrum histogram.
		fPtSpectrumTOFTPCbins = new TH1F("fPtSpectrumTOFTPCbins","p_{T} Spectrum;p_{T} (GeV/c);Count",nptbins,ptarrayTOFTPC);
		fOutputList->Add(fPtSpectrumTOFTPCbins);

		// Create unidentified correlations histogram.
		fCorrelationsTOFTPCbins = AliHistToolsDiHadronPID::MakeHist3D("fCorrelationsTOFTPCbins","Correlations;#Delta#phi;#Delta#eta;p_{T} (GeV/c)",
			fNDPhiBins,-TMath::Pi()/2.,3.*TMath::Pi()/2.,
			fNDEtaBins,-1.6,1.6,fExtendPtAxis ? nptbinsext : nptbins, fExtendPtAxis ? ptarrayTOFTPCext : ptarrayTOFTPC);
		fOutputList->Add(fCorrelationsTOFTPCbins);

		// Create unidentified mixed events histogram.
		fMixedEventsTOFTPCbins = AliHistToolsDiHadronPID::MakeHist3D("fMixedEventsTOFTPCbins","Mixed Events;#Delta#phi;#Delta#eta;p_{T} (GeV/c)",
			fNDPhiBins,-TMath::Pi()/2.,3.*TMath::Pi()/2.,
			fNDEtaBins,-1.6,1.6,fExtendPtAxis ? nptbinsext : nptbins, fExtendPtAxis ? ptarrayTOFTPCext : ptarrayTOFTPC);
		fOutputList->Add(fMixedEventsTOFTPCbins);

		fTOFTPChistos = new TObjArray(3);
		fTOFTPChistos->SetOwner(kTRUE);
		fTOFTPChistos->SetName("CorrelationsTOFTPC");

		if (fCalculateMismatch) {
			fTOFTPCmismatch = new TObjArray(3);
			fTOFTPCmismatch->SetOwner(kTRUE);
			fTOFTPCmismatch->SetName("MismatchTOFTPC");
		}

		fMixedEventsTOFTPCbinsPID = new TObjArray(3);
		fMixedEventsTOFTPCbinsPID->SetOwner(kTRUE);
		fMixedEventsTOFTPCbinsPID->SetName("MixedEventsTOFTPC");

		for (Int_t iSpecies = 0; iSpecies < 3; iSpecies++) {

			// Create Mixed events with PID.
			TH3F* mixedeventsPID = AliHistToolsDiHadronPID::MakeHist3D(Form("fMixedEventsTOFTPC%s", speciesname[iSpecies].Data()),
			Form("Mixed Events %s;#Delta#phi;#Delta#eta;p_{T} (GeV/c)", speciesname[iSpecies].Data()),
			fNDPhiBins,-TMath::Pi()/2.,3.*TMath::Pi()/2.,
			fNDEtaBins,-1.6,1.6,fExtendPtAxis ? nptbinsext : nptbins, fExtendPtAxis ? ptarrayTOFTPCext : ptarrayTOFTPC);
			fMixedEventsTOFTPCbinsPID->Add(mixedeventsPID);
		
			// Create the directory structure Pion, Kaon, Proton, regardless
			// of wether the histograms are created (to keep the order.)
			TObjArray* TOFTPChistosTmp = new TObjArray(fTOFTPCPtAxis->GetNbins());
			TOFTPChistosTmp->SetOwner(kTRUE);
			TOFTPChistosTmp->SetName(speciesname[iSpecies].Data());

			TObjArray* TOFTPCmismatchTmp = 0x0;
			if (fCalculateMismatch) { 
				TOFTPCmismatchTmp = new TObjArray(fTOFTPCPtAxis->GetNbins());
				TOFTPCmismatchTmp->SetOwner(kTRUE);
				TOFTPCmismatchTmp->SetName(speciesname[iSpecies].Data());	
			}

			// Only Create the TOF/TPC histograms when requested.
			Bool_t MakeTOFTPCcorrelations[3] = {fMakeTOFTPCcorrelationsPi, fMakeTOFTPCcorrelationsKa, fMakeTOFTPCcorrelationsPr};
			if (MakeTOFTPCcorrelations[iSpecies]) {
				for (Int_t iBinPt = 1; iBinPt < (fTOFTPCPtAxis->GetNbins() + 1); iBinPt++) {
			
					// Approximate resolutions of TOF and TPC detector.
					const Double_t sTOFest = 110.;
					const Double_t sTPCest = 4.5;
					
					// Set range +/- 5 sigma of main peak. (+ 10 sigma for TOF max, for mismatches.)
					Double_t TOFmin = -5. * sTOFest;
					Double_t TOFmax = 10. * sTOFest;
					Double_t TPCmin = -4. * sTPCest;
					Double_t TPCmax = 4. * sTPCest;

					Double_t TOFexp = AliFunctionsDiHadronPID::TOFExpTime(fTOFTPCPtAxis->GetBinLowEdge(iBinPt), 0.4,  AliFunctionsDiHadronPID::M(iSpecies));
					Double_t TPCexp = AliFunctionsDiHadronPID::TPCExpdEdX(fTOFTPCPtAxis->GetBinLowEdge(iBinPt), 0.4,  AliFunctionsDiHadronPID::M(iSpecies));

					for (Int_t jSpecies = 0; jSpecies < 3; jSpecies++) {

						if (iSpecies == jSpecies) {continue;}

						Double_t TOFexpOther = AliFunctionsDiHadronPID::TOFExpTime(fTOFTPCPtAxis->GetBinLowEdge(iBinPt), 0.4,  AliFunctionsDiHadronPID::M(jSpecies));
						Double_t TPCexpOther = AliFunctionsDiHadronPID::TPCExpdEdX(fTOFTPCPtAxis->GetBinLowEdge(iBinPt), 0.4,  AliFunctionsDiHadronPID::M(jSpecies));
					
						// If any peak is within +/- 7 sigma, then also add this peak.
						if ( (TMath::Abs(TOFexp - TOFexpOther) < 7. * sTOFest) ||
							 (TMath::Abs(TPCexp - TPCexpOther) < 7. * sTPCest) ) {

							TOFmin = TMath::Min(TOFmin, (TOFexpOther - TOFexp - 2. * sTOFest) );
							TOFmax = TMath::Max(TOFmax, (TOFexpOther - TOFexp + 10. * sTOFest) );
							TPCmin = TMath::Min(TPCmin, (TPCexpOther - TPCexp - 2. * sTPCest) );
							TPCmax = TMath::Max(TPCmax, (TPCexpOther - TPCexp + 2. * sTPCest) );						

						}

					}

					// With the standard TOF range, fitting the deuterons and the TOF mismatches is
					// hard. This flag doubles the range of the TOF axis in the TOF/TPC histograms,
					// while leaving the resolution the same. Turning on this flag will greatly increase
					// the memory consumption of the task, to the point that it's probably too much 
					// to save a Buffer with all three species included.
					Double_t TOFreach = TOFmax - TOFmin;
					TOFmax += (TOFreach * (fTOFIntervalFactorTOFTPC - 1.));
					Int_t TOFbins = (Int_t)(60. * fTOFIntervalFactorTOFTPC);

					Int_t NBinsTOFTPC[4] = {32, 32, TOFbins, 40};
					Double_t minTOFTPC[4] = {-TMath::Pi()/2., -1.6, TOFmin, TPCmin};
					Double_t maxTOFTPC[4] = {3.*TMath::Pi()/2., 1.6, TOFmax, TPCmax};

					THnF* htmp = new THnF(Form("fCorrelationsTOFTPC_%i",iBinPt),
						Form("%5.3f < p_{T} < %5.3f", fTOFTPCPtAxis->GetBinLowEdge(iBinPt), fTOFTPCPtAxis->GetBinUpEdge(iBinPt)), 
						4, NBinsTOFTPC, minTOFTPC, maxTOFTPC);

					(htmp->GetAxis(0))->SetTitle("#Delta#phi");
					(htmp->GetAxis(1))->SetTitle("#Delta#eta");
					(htmp->GetAxis(2))->SetTitle("t_{TOF} (ps)");
					(htmp->GetAxis(3))->SetTitle("dE/dx (a.u.)");

					TOFTPChistosTmp->Add(htmp);

					if (fCalculateMismatch) { 
						// Mismatch histogram.
						TH2F* htmp2 = new TH2F(Form("fMismatchTOFTPC_%i",iBinPt),
							Form("%5.3f < p_{T} < %5.3f; t_{TOF} (ps); dE/dx (a.u.)", fTOFTPCPtAxis->GetBinLowEdge(iBinPt), fTOFTPCPtAxis->GetBinUpEdge(iBinPt)), 
							NBinsTOFTPC[2], TOFmin, TOFmax, NBinsTOFTPC[3], TPCmin, TPCmax);
						htmp2->SetDirectory(0);

						TOFTPCmismatchTmp->Add(htmp2);
		
					}
				} // End loop over pT bins.
			} // End species if.

			fTOFTPChistos->Add(TOFTPChistosTmp);
			if (fCalculateMismatch) {fTOFTPCmismatch->Add(TOFTPCmismatchTmp);}

		}

		fOutputList->Add(fTOFTPChistos);
		if (fCalculateMismatch) {fOutputList->Add(fTOFTPCmismatch);}
		fOutputList->Add(fMixedEventsTOFTPCbinsPID);

	}

	// Load external TOF histograms if flag is set.
	if (fCalculateMismatch) {LoadExtMismatchHistos();}

	PostData(1,fOutputList);

}

// -----------------------------------------------------------------------
void AliAnalysisTaskDiHadronPID::LocalInit() {

	//
	// Initialize on the this computer. 
	//

	if (fDebug > 0) {cout << Form("File: %s, Line: %i, Function: %s",__FILE__,__LINE__,__func__) << endl;}

}

// -----------------------------------------------------------------------
void AliAnalysisTaskDiHadronPID::UserExec(Option_t*) {

	//
	// Main Loop.
	//

	if (fDebug > 0) {cout << Form("File: %s, Line: %i, Function: %s",__FILE__,__LINE__,__func__) << endl;}

	// Input Current Event.
	fCurrentAODEvent = dynamic_cast<AliAODEvent*>(InputEvent());
	if (!fCurrentAODEvent) AliFatal("No Event Found!");

	if (!fEventCuts->IsSelected(fCurrentAODEvent)) {return;}

	// Fill the global tracks array. - NOT NEEDED I THINK, since we're not using
	// bit 1<<7 for the associated tracks!

	// Let the track cut objects know that a new event will start.
	fTrackCutsTrigger->StartNewEvent();
	fTrackCutsAssociated->StartNewEvent();

	// Create arrays for trigger/associated particles.
	fTriggerTracks = new TObjArray(350);
	fTriggerTracks->SetOwner(kTRUE);

	fAssociatedTracks = new TObjArray(3500);
	fAssociatedTracks->SetOwner(kTRUE);

	for (Int_t iTrack = 0; iTrack < fCurrentAODEvent->GetNumberOfTracks(); iTrack++) {

		AliAODTrack* track = (AliAODTrack*)fCurrentAODEvent->GetTrack(iTrack);
		AliTrackDiHadronPID* pidtrack = new AliTrackDiHadronPID(track,0x0,0x0,fPIDResponse);
		pidtrack->ForgetAboutPointers();
		pidtrack->SetDebugLevel(fDebug);

		Double_t rndhittime = -1.e21;
		if (fCalculateMismatch) rndhittime = GenerateRandomHit(pidtrack->Eta());

		// Fill the trigger/associated tracks array.
		if (fTrackCutsTrigger->IsSelectedData(pidtrack,rndhittime)) {fTriggerTracks->AddLast(pidtrack);}
		else if (fTrackCutsAssociated->IsSelectedData(pidtrack,rndhittime)) {
			
			fAssociatedTracks->AddLast(pidtrack);

			// Fill p_T spectrum.
			if (fPtSpectrumTOFbins) fPtSpectrumTOFbins->Fill(pidtrack->Pt());
			if (fPtSpectrumTOFTPCbins) fPtSpectrumTOFTPCbins->Fill(pidtrack->Pt());

			// Fill mismatch histograms with associateds.
			if (fCalculateMismatch && (rndhittime > -1.e20)) {

				Double_t apt = pidtrack->Pt();

				if (fMakeTOFcorrelations) {

					for (Int_t iSpecies = 0; iSpecies < 3; iSpecies++) {

						TObjArray* atmp = (TObjArray*)fTOFmismatch->At(iSpecies);
						Int_t ptbin = fTOFPtAxis->FindBin(apt);

						// Only fill if histogram exists in fTOFmismatch.
						if ( !(ptbin < 1) && !(ptbin > fTOFPtAxis->GetNbins()) ) {

							TH1F* htmp = (TH1F*)atmp->At(ptbin - 1);
							htmp->Fill(rndhittime - pidtrack->GetTOFsignalExpected(iSpecies));

						}
					}
				}

				Bool_t MakeTOFTPCcorrelations[3] = {fMakeTOFTPCcorrelationsPi, fMakeTOFTPCcorrelationsKa, fMakeTOFTPCcorrelationsPr};
				if (fMakeTOFTPCcorrelationsPi || fMakeTOFTPCcorrelationsKa || fMakeTOFTPCcorrelationsPr) { 

					for (Int_t iSpecies = 0; iSpecies < 3; iSpecies++) {

						if (!MakeTOFTPCcorrelations[iSpecies]) {continue;}

						TObjArray* atmp = (TObjArray*)fTOFTPCmismatch->At(iSpecies);
						Int_t ptbin = fTOFTPCPtAxis->FindBin(apt);

						// Only fill if histogram exists in fTOFTPCmismatch.
						if ( !(ptbin < 1) && !(ptbin > fTOFTPCPtAxis->GetNbins()) ) {

							TH2F* htmp = (TH2F*)atmp->At(ptbin - 1);
							htmp->Fill(rndhittime - pidtrack->GetTOFsignalExpected(iSpecies), pidtrack->GetTPCsignalMinusExpected(iSpecies));

						}
					}
				}

			}

		} 
		else {delete pidtrack;}

	}

	// Fill Correlation histograms.
	for (Int_t iTrigger = 0; iTrigger < fTriggerTracks->GetEntriesFast(); iTrigger++) {
		AliTrackDiHadronPID* triggertrack = (AliTrackDiHadronPID*)fTriggerTracks->At(iTrigger);

		for (Int_t iAssociated = 0; iAssociated < fAssociatedTracks->GetEntriesFast(); iAssociated++) {
			AliTrackDiHadronPID* associatedtrack = (AliTrackDiHadronPID*)fAssociatedTracks->At(iAssociated);

			Double_t DPhi = triggertrack->Phi() - associatedtrack->Phi();
			if (DPhi < -TMath::Pi()/2.) {DPhi += 2.*TMath::Pi();}
			if (DPhi > 3.*TMath::Pi()/2.) {DPhi -= 2.*TMath::Pi();}

			Double_t DEta = triggertrack->Eta() - associatedtrack->Eta();
			if (fCorrelationsTOFbins) fCorrelationsTOFbins->Fill(DPhi,DEta,associatedtrack->Pt());
			if (fCorrelationsTOFTPCbins) fCorrelationsTOFTPCbins->Fill(DPhi,DEta,associatedtrack->Pt());

			Double_t apt = associatedtrack->Pt();

			// Fill TOF correlations.
			if (fMakeTOFcorrelations) {
				
				for (Int_t iSpecies = 0; iSpecies < 3; iSpecies++) {

					TObjArray* atmp = (TObjArray*)fTOFhistos->At(iSpecies);
					Int_t ptbin = fTOFPtAxis->FindBin(apt);

					// Only fill if histogram exists in fTOFhistos.
					if ( !(ptbin < 1) && !(ptbin > fTOFPtAxis->GetNbins()) ) {

						TH3F* htmp = (TH3F*)atmp->At(ptbin - 1);
						htmp->Fill(DPhi, DEta, associatedtrack->GetTOFsignalMinusExpected(iSpecies));

					}
				}
			}

			// Fill TOF/ TPC Correlations.
			Bool_t MakeTOFTPCcorrelations[3] = {fMakeTOFTPCcorrelationsPi, fMakeTOFTPCcorrelationsKa, fMakeTOFTPCcorrelationsPr};
			if (fMakeTOFTPCcorrelationsPi || fMakeTOFTPCcorrelationsKa || fMakeTOFTPCcorrelationsPr) { 

				for (Int_t iSpecies = 0; iSpecies < 3; iSpecies++) {

					if (!MakeTOFTPCcorrelations[iSpecies]) {continue;}

					TObjArray* atmp = (TObjArray*)fTOFTPChistos->At(iSpecies);
					Int_t ptbin = fTOFTPCPtAxis->FindBin(apt);

					// Only fill if histogram exists in fTOFhistos.
					if ( !(ptbin < 1) && !(ptbin > fTOFTPCPtAxis->GetNbins()) ) {

						THnF* htmp = (THnF*)atmp->At(ptbin - 1);
						Double_t TOFTPCfill[4] = {DPhi, DEta, 
							associatedtrack->GetTOFsignalMinusExpected(iSpecies), associatedtrack->GetTPCsignalMinusExpected(iSpecies)}; 

						htmp->Fill(TOFTPCfill);

					}				
				}	
			}		
		}
	}

	//cout<<"Triggers: "<<fTriggerTracks->GetEntriesFast()<<" Associateds: "<<fAssociatedTracks->GetEntriesFast()<<endl;	

	// Determine vtxz of current event.
	AliAODVertex* currentprimaryvertex = fCurrentAODEvent->GetPrimaryVertex();
	Double_t vtxz = currentprimaryvertex->GetZ();

	// Determine centrality of current event (for PbPb).
	AliEventPool* poolin = 0x0;
	Float_t percentile = -1.;
	if (fEventCuts->GetIsPbPb()) {
		TString centralityestimator = fEventCuts->GetCentralityEstimator();
		AliCentrality* currentcentrality = fCurrentAODEvent->GetCentrality();
		percentile = currentcentrality->GetCentralityPercentile(centralityestimator.Data());

		poolin = fPoolMgr->GetEventPool(percentile, vtxz); 
		if (!poolin) {AliFatal(Form("No pool found for centrality = %f, vtxz = %f", percentile, vtxz));}
	} else {
		poolin = fPoolMgr->GetEventPool(0.5, vtxz);	// There are no multiplicity bins for pp yet.  
		if (!poolin) {AliFatal(Form("No pool found for vtxz = %f", vtxz));}
	}

	// TObjArray* fGlobalTracksArray; 

	// Give a print out of the pool manager's contents.
	if (fDebug > 0) PrintPoolManagerContents();

	// Mix events if there are enough events in the pool.
	if (poolin->GetCurrentNEvents() >= fMinNEventsForMixing) {
		//{cout << "Mixing Events." << endl;}

		// Loop over all events in the event pool.
		for (Int_t iMixEvent = 0; iMixEvent < poolin->GetCurrentNEvents(); iMixEvent++) {
	    	TObjArray* mixtracks = poolin->GetEvent(iMixEvent);

	    	// Mix either the triggers or the associateds.
	    	if (fMixTriggers) {

				// Loop over all associateds in this event.
				for (Int_t iAssociated = 0; iAssociated < fAssociatedTracks->GetEntriesFast(); iAssociated++) {
					AliTrackDiHadronPID* associatedtrack = (AliTrackDiHadronPID*)fAssociatedTracks->At(iAssociated);

		    		// Loop over all mixed tracks.
		    		for (Int_t iMixTrack = 0; iMixTrack < mixtracks->GetEntriesFast(); iMixTrack++) {
		    			AliTrackDiHadronPID* mixtrack = (AliTrackDiHadronPID*)mixtracks->At(iMixTrack);
							
						Double_t DPhi = mixtrack->Phi() - associatedtrack->Phi();
						if (DPhi < -TMath::Pi()/2.) {DPhi += 2.*TMath::Pi();}
						if (DPhi > 3.*TMath::Pi()/2.) {DPhi -= 2.*TMath::Pi();}

						Double_t DEta = mixtrack->Eta() - associatedtrack->Eta();
						if (fMixedEventsTOFbins) fMixedEventsTOFbins->Fill(DPhi,DEta,associatedtrack->Pt());
						if (fMixedEventsTOFTPCbins) fMixedEventsTOFTPCbins->Fill(DPhi,DEta,associatedtrack->Pt());

						// Fill the mixed event histograms with a 1 sigma PID cut.
						if (fMixedEventsTOFTPCbinsPID) {

							for (Int_t iSpecies = 0; iSpecies < 3; iSpecies++) {

								TH3F* mixedeventhist = (TH3F*)fMixedEventsTOFTPCbinsPID->At(iSpecies);

								// Check the nSigma of the associated tracks.
								Double_t nSigmaTOFTPC = TMath::Sqrt( 
								associatedtrack->GetNumberOfSigmasTOF(iSpecies) * associatedtrack->GetNumberOfSigmasTOF(iSpecies) +
								associatedtrack->GetNumberOfSigmasTPC(iSpecies) * associatedtrack->GetNumberOfSigmasTPC(iSpecies));

								if (nSigmaTOFTPC < 1.) {mixedeventhist->Fill(DPhi,DEta,associatedtrack->Pt());}

							}
						}

		    		}
		   		}

		   	} else {

				// Loop over all triggers in this event.
				for (Int_t iTrigger = 0; iTrigger < fTriggerTracks->GetEntriesFast(); iTrigger++) {
					AliTrackDiHadronPID* triggertrack = (AliTrackDiHadronPID*)fTriggerTracks->At(iTrigger);

		    		// Loop over all mixed tracks.
		    		for (Int_t iMixTrack = 0; iMixTrack < mixtracks->GetEntriesFast(); iMixTrack++) {
		    			AliTrackDiHadronPID* mixtrack = (AliTrackDiHadronPID*)mixtracks->At(iMixTrack);
							
						Double_t DPhi = triggertrack->Phi() - mixtrack->Phi();
						if (DPhi < -TMath::Pi()/2.) {DPhi += 2.*TMath::Pi();}
						if (DPhi > 3.*TMath::Pi()/2.) {DPhi -= 2.*TMath::Pi();}

						Double_t DEta = triggertrack->Eta() - mixtrack->Eta();
						if (fMixedEventsTOFbins) fMixedEventsTOFbins->Fill(DPhi,DEta,mixtrack->Pt());
						if (fMixedEventsTOFTPCbins) fMixedEventsTOFTPCbins->Fill(DPhi,DEta,mixtrack->Pt());
		    		
						// Fill the mixed event histograms with a 1 sigma PID cut.
						if (fMixedEventsTOFTPCbinsPID) {

							for (Int_t iSpecies = 0; iSpecies < 3; iSpecies++) {

								TH3F* mixedeventhist = (TH3F*)fMixedEventsTOFTPCbinsPID->At(iSpecies);

								// Check the nSigma of the associated tracks.
								Double_t nSigmaTOFTPC = TMath::Sqrt( 
								mixtrack->GetNumberOfSigmasTOF(iSpecies) * mixtrack->GetNumberOfSigmasTOF(iSpecies) +
								mixtrack->GetNumberOfSigmasTPC(iSpecies) * mixtrack->GetNumberOfSigmasTPC(iSpecies));

								if (nSigmaTOFTPC < 1.) {mixedeventhist->Fill(DPhi,DEta,mixtrack->Pt());}

							}
						}

		    		}
		   		}

		   	} // End if  	
	   	}
	}	

	// Update the event pool.
	AliEventPool* poolout = 0x0;
	if (fEventCuts->GetIsPbPb()) {
		poolout = fPoolMgr->GetEventPool(percentile, vtxz); // Get the buffer associated with the current centrality and z-vtx
		if (!poolout) AliFatal(Form("No pool found for centrality = %f, vtx_z = %f", percentile, vtxz));
	} else {
		poolout = fPoolMgr->GetEventPool(0.5, vtxz); // Get the buffer associated with the current centrality and z-vtx
		if (!poolout) AliFatal(Form("No pool found for vtx_z = %f", vtxz));
	}


	// Q: is it a problem that the fAssociatedTracks array can be bigger than the number of tracks inside?
	if (fMixTriggers) {
		poolout->UpdatePool(fTriggerTracks);
		fAssociatedTracks->Delete();
		delete fAssociatedTracks;
	}
	else {
		poolout->UpdatePool(fAssociatedTracks);
		fTriggerTracks->Delete();
		delete fTriggerTracks;
	}

	fTriggerTracks = 0x0;
	fAssociatedTracks = 0x0;

	// Tell the track cut object that the event is done.
	fTrackCutsTrigger->EventIsDone(0);
	fTrackCutsAssociated->EventIsDone(0);

	PostData(1,fOutputList);

}

// -----------------------------------------------------------------------
void AliAnalysisTaskDiHadronPID::SelectCollisionCandidates(UInt_t offlineTriggerMask) {

	// Overrides the method defined in AliAnalysisTaskSE. This is needed because
	// the event selection is not done in the task, but in the AliAODEventCutsDiHadronPID class.

	if (fDebug > 0) {cout << Form("File: %s, Line: %i, Function: %s",__FILE__,__LINE__,__func__) << endl;}
	if (!fEventCuts) {cout << Form("%s -> ERROR: No AliAODEventCutsDiHadronPID class created for the analysis...",__func__) << endl; return;}

	//fOfflineTriggerMask = offlineTriggerMask;
	fEventCuts->SetTrigger(offlineTriggerMask);

}

// -----------------------------------------------------------------------
void AliAnalysisTaskDiHadronPID::SetDebugLevel(Int_t level) {

	// Also propagates this setting to the track and event cuts.
	if (fDebug > 0) {cout << Form("File: %s, Line: %i, Function: %s",__FILE__,__LINE__,__func__) << endl;}

	fDebug = level;

	if (fEventCuts) {fEventCuts->SetDebugLevel(level);}
	if (fTrackCutsTrigger) {fTrackCutsTrigger->SetDebugLevel(level);}
	if (fTrackCutsAssociated) {fTrackCutsAssociated->SetDebugLevel(level);}

} 

// -----------------------------------------------------------------------
Bool_t AliAnalysisTaskDiHadronPID::LoadExtMismatchHistos() {

	//
	// Attempting to load a root file containing information needed
	// to generate random TOF hits.
 	//

	if (fDebug > 0) {cout << Form("File: %s, Line: %i, Function: %s",__FILE__,__LINE__,__func__) << endl;}

	// Opening external TOF file.
	if (fDebug > 0) cout<<"Trying to open TOFmismatchHistos.root ..."<<endl;
	TFile* fin = 0x0;
	fin = TFile::Open("alien:///alice/cern.ch/user/m/mveldhoe/rootfiles/TOFmismatchHistos.root");
	if (!fin) {
		AliWarning("Couln't open TOFmismatchHistos, will not calculate mismatches...");
		fCalculateMismatch = kFALSE;
		return kFALSE;
	}

	// Check if the required histograms are present.
	TH1F* tmp1 = (TH1F*)fin->Get("hNewT0Fill");
	if (!tmp1) {
		AliWarning("Couln't find hNewT0Fill, will not calculate mismatches...");
		fCalculateMismatch = kFALSE;
		return kFALSE;	
	}
	TH2F* tmp2 = (TH2F*)fin->Get("hLvsEta");
	if (!tmp2) {
		AliWarning("Couln't find hLvsEta, will not calculate mismatches...");
		fCalculateMismatch = kFALSE;
		return kFALSE;	
	}	

	// Make a deep copy of the files in the histogram.
	fT0Fill = (TH1F*)tmp1->Clone("fT0Fill");
	fLvsEta = (TH2F*)tmp2->Clone("fLvsEta");

	// Close the external file.
	AliInfo("Closing external file.");
	fin->Close();

	// Creating a TObjArray for LvsEta projections.
	const Int_t nbinseta = fLvsEta->GetNbinsX();
	fLvsEtaProjections = new TObjArray(nbinseta);
	fLvsEtaProjections->SetOwner(kTRUE);

	// Making the projections needed (excluding underflow/ overflow).
	for (Int_t iEtaBin = 1; iEtaBin < (nbinseta + 1); iEtaBin++) {
		TH1F* tmp = (TH1F*)fLvsEta->ProjectionY(Form("LvsEtaProjection_%i",iEtaBin),iEtaBin,iEtaBin);
		tmp->SetDirectory(0);
		fLvsEtaProjections->AddAt(tmp,iEtaBin - 1);
	}

	return kTRUE;

}

// -----------------------------------------------------------------------
Double_t AliAnalysisTaskDiHadronPID::GenerateRandomHit(Double_t eta) {

	//
	// Returns a random TOF time.
	//

	if (fDebug > 0) {cout << Form("File: %s, Line: %i, Function: %s",__FILE__,__LINE__,__func__) << endl;}

	// Default (error) value:
	Double_t rndhittime = -1.e21;

	// TOF mismatch flag is not turned on.
	if (!fCalculateMismatch) {
		AliFatal("Called GenerateRandomHit() method, but flag fCalculateMismatch not set.");
		return rndhittime;
	}

	// TOF doesn't extend much further than 0.8.
	if (TMath::Abs(eta) > 0.8) {
		if (fDebug) {AliInfo("Tried to get a random hit for a track with eta > 0.8.");}
		return rndhittime;
	}

	// Finding the bin of the eta.
	TAxis* etaAxis = fLvsEta->GetXaxis();
	Int_t etaBin = etaAxis->FindBin(eta);
	if (etaBin == 0 || (etaBin == etaAxis->GetNbins() + 1)) {return rndhittime;}

	const TH1F* lengthDistribution = (const TH1F*)fLvsEtaProjections->At(etaBin - 1);

	if (!lengthDistribution) {
		AliFatal("length Distribution not found.");
		return rndhittime;
	}

	Double_t currentRndLength = lengthDistribution->GetRandom(); // in cm.

	// Similar to Roberto's code.
	Double_t currentRndTime = currentRndLength / (TMath::C() * 1.e2 / 1.e12);
	Double_t t0fill = -1.26416e+04;
	rndhittime = fT0Fill->GetRandom() - t0fill + currentRndTime;

	return rndhittime;

}

// -----------------------------------------------------------------------
void AliAnalysisTaskDiHadronPID::PrintPoolManagerContents() {

	//
	// Prints out the current contents of the event pool manager.
	//

	// Determine the number of pools in the pool manager.
	AliEventPool* poolin = fPoolMgr->GetEventPool(0,0);
	Int_t NPoolsCentrality = 0;
	while (poolin) {
		NPoolsCentrality++;
		poolin = fPoolMgr->GetEventPool(NPoolsCentrality,0);
	} 

	poolin = fPoolMgr->GetEventPool(0,0);
	Int_t NPoolsVtxZ = 0;	
	while (poolin) {
		NPoolsVtxZ++;
		poolin = fPoolMgr->GetEventPool(0,NPoolsVtxZ);
	} 

	// Loop over all Pools in the matrix of the pool manager.
	cout<<" Pool manager contents: (Nevt,NTrack)"<<endl;
	for (Int_t iCentrality = 0; iCentrality < NPoolsCentrality; iCentrality++) {
		cout<<Form("Centrality Bin: %2i --> ", iCentrality);

		for (Int_t iVtxZ = 0; iVtxZ < NPoolsVtxZ; iVtxZ++) {

			poolin = fPoolMgr->GetEventPool(iCentrality, iVtxZ);

			cout<<Form("(%2i,%4i) ",poolin->GetCurrentNEvents(), poolin->NTracksInPool());

		}

		cout<<endl;
	}

}

// -----------------------------------------------------------------------
void AliAnalysisTaskDiHadronPID::Terminate(Option_t*) {;

	//
	// Called when task is done.
	//

	if (fDebug > 0) {cout << Form("File: %s, Line: %i, Function: %s",__FILE__,__LINE__,__func__) << endl;}

	delete fT0Fill;
	fT0Fill = 0x0;
	delete fLvsEta;
	fLvsEta = 0x0;
	delete fLvsEtaProjections;
	fLvsEtaProjections = 0x0;

}
