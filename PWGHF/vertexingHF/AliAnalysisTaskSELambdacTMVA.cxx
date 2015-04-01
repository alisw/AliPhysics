/**************************************************************************
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

/* $Id: AliAnalysisTaskSELambdacTMVA.cxx$ */

//*************************************************************************
// AliAnalysisTaskSE for the Lambdac candidates, for TMVA analysis,
// and checks on MC generated and reconstructed Lambdac
// 
// Modified from AliAnalysisTaskSELambdac
// Authors: Jaime Norman (jaime.norman@cern.ch)
//          Marcel Figueredo (marcel.figueredo@cern.ch)
//*************************************************************************


#include <TClonesArray.h>
#include <TNtuple.h>
#include <TCanvas.h>
#include <TList.h>
#include <TString.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TDatabasePDG.h>

#include <AliAnalysisDataSlot.h>
#include <AliAnalysisDataContainer.h>
#include "AliAnalysisManager.h"
#include "AliAODHandler.h"
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliAODRecoDecayHF3Prong.h"
#include "AliAnalysisVertexingHF.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskSELambdacTMVA.h"
#include "AliKFParticle.h"
#include "AliAODPidHF.h"
#include "AliRDHFCutsLctopKpi.h"
#include "AliRDHFCuts.h"
#include "AliKFVertex.h"
#include "AliESDVertex.h"
//#include "AliAODpidUtil.h"
#include "AliAODPid.h"
#include "AliInputEventHandler.h"
#include "AliPID.h"
#include "AliNormalizationCounter.h"
#include "AliVertexingHFUtils.h"
#include "AliCFVertexingHF3Prong.h"
#include "AliCFTaskVertexingHF.h"
#include "AliLog.h" 

ClassImp(AliAnalysisTaskSELambdacTMVA)


	//________________________________________________________________________
	AliAnalysisTaskSELambdacTMVA::AliAnalysisTaskSELambdacTMVA():
		AliAnalysisTaskSE(),
		fOutput(0), 
		fHistNEvents(0),
		fHistNEventsRejTM(0),
		fhSelectBit(0),
		fNtupleLambdac(0),
		fIsLc(0),
		fIsLcResonant(0),
		fRDCutsAnalysis(0),
		fListCuts(0),
		fFillNtuple(0),
		fReadMC(kFALSE),
		fMCPid(kFALSE),
		fRealPid(kFALSE),
		fResPid(kTRUE),
		fUseKF(kFALSE),
		fAnalysis(kFALSE),
		fVHF(0),
		fLcCut(kFALSE),
		fLcPIDCut(kFALSE),    
		fNentries(0),
		fPIDResponse(0),
		fCounter(0)

{
	//
	// Default constructor
	//

	for(Int_t i=0;i<11;i++) { 
		fhNBkgNI[i]=0x0;
		fhNLc[i]=0x0;
		fhNLcc[i]=0x0;
		fhNLcNonRc[i]=0x0;
		fhNLcL1520c[i]=0x0;
		fhNLcKstarc[i]=0x0;
		fhNLcDeltac[i]=0x0;
		fhNLcb[i]=0x0;
		fhNLcNonRb[i]=0x0;
		fhNLcL1520b[i]=0x0;
		fhNLcKstarb[i]=0x0;
		fhNLcDeltab[i]=0x0;

		fhPtEtaBkgNI[i]=0x0;
		fhPtEtaLc[i]=0x0;
		fhPtEtaLcc[i]=0x0;
		fhPtEtaLcNonRc[i]=0x0;
		fhPtEtaLcL1520c[i]=0x0;
		fhPtEtaLcKstarc[i]=0x0;
		fhPtEtaLcDeltac[i]=0x0;
		fhPtEtaLcb[i]=0x0;
		fhPtEtaLcNonRb[i]=0x0;
		fhPtEtaLcL1520b[i]=0x0;
		fhPtEtaLcKstarb[i]=0x0;
		fhPtEtaLcDeltab[i]=0x0;

		fhPtYBkgNI[i]=0x0;
		fhPtYLc[i]=0x0;
		fhPtYLcc[i]=0x0;
		fhPtYLcNonRc[i]=0x0;
		fhPtYLcL1520c[i]=0x0;
		fhPtYLcKstarc[i]=0x0;
		fhPtYLcDeltac[i]=0x0;
		fhPtYLcb[i]=0x0;
		fhPtYLcNonRb[i]=0x0;
		fhPtYLcL1520b[i]=0x0;
		fhPtYLcKstarb[i]=0x0;
		fhPtYLcDeltab[i]=0x0;

		fhPtPhiBkgNI[i]=0x0;
		fhPtPhiLc[i]=0x0;
		fhPtPhiLcc[i]=0x0;
		fhPtPhiLcNonRc[i]=0x0;
		fhPtPhiLcL1520c[i]=0x0;
		fhPtPhiLcKstarc[i]=0x0;
		fhPtPhiLcDeltac[i]=0x0;
		fhPtPhiLcb[i]=0x0;
		fhPtPhiLcNonRb[i]=0x0;
		fhPtPhiLcL1520b[i]=0x0;
		fhPtPhiLcKstarb[i]=0x0;
		fhPtPhiLcDeltab[i]=0x0;

	}
}

//________________________________________________________________________
AliAnalysisTaskSELambdacTMVA::AliAnalysisTaskSELambdacTMVA(const char *name,Int_t fillNtuple,AliRDHFCutsLctopKpi *lccutsana):
	AliAnalysisTaskSE(name),
	fOutput(0), 
	fHistNEvents(0),
	fHistNEventsRejTM(0),
	fhSelectBit(0),
	fNtupleLambdac(0),
	fIsLc(0),
	fIsLcResonant(0),
	fRDCutsAnalysis(lccutsana),
	fListCuts(0),
	fFillNtuple(fillNtuple),
	fReadMC(kFALSE),
	fMCPid(kFALSE),
	fRealPid(kFALSE),
	fResPid(kTRUE),
	fUseKF(kFALSE),
	fAnalysis(kFALSE),
	fVHF(0),
	fLcCut(kFALSE),
	fLcPIDCut(kFALSE),    
	fNentries(0),
	fPIDResponse(0),
	fCounter(0)
{
	//
	// Default constructor
	// Output slot #1 writes into a TList container
	//
	for(Int_t i=0;i<11;i++) { 
		fhNBkgNI[i]=0x0;
		fhNLc[i]=0x0;
		fhNLcc[i]=0x0;
		fhNLcNonRc[i]=0x0;
		fhNLcL1520c[i]=0x0;
		fhNLcKstarc[i]=0x0;
		fhNLcDeltac[i]=0x0;
		fhNLcb[i]=0x0;
		fhNLcNonRb[i]=0x0;
		fhNLcL1520b[i]=0x0;
		fhNLcKstarb[i]=0x0;
		fhNLcDeltab[i]=0x0;

		fhPtEtaBkgNI[i]=0x0;
		fhPtEtaLc[i]=0x0;
		fhPtEtaLcc[i]=0x0;
		fhPtEtaLcNonRc[i]=0x0;
		fhPtEtaLcL1520c[i]=0x0;
		fhPtEtaLcKstarc[i]=0x0;
		fhPtEtaLcDeltac[i]=0x0;
		fhPtEtaLcb[i]=0x0;
		fhPtEtaLcNonRb[i]=0x0;
		fhPtEtaLcL1520b[i]=0x0;
		fhPtEtaLcKstarb[i]=0x0;
		fhPtEtaLcDeltab[i]=0x0;

		fhPtYBkgNI[i]=0x0;
		fhPtYLc[i]=0x0;
		fhPtYLcc[i]=0x0;
		fhPtYLcNonRc[i]=0x0;
		fhPtYLcL1520c[i]=0x0;
		fhPtYLcKstarc[i]=0x0;
		fhPtYLcDeltac[i]=0x0;
		fhPtYLcb[i]=0x0;
		fhPtYLcNonRb[i]=0x0;
		fhPtYLcL1520b[i]=0x0;
		fhPtYLcKstarb[i]=0x0;
		fhPtYLcDeltab[i]=0x0;

		fhPtPhiBkgNI[i]=0x0;
		fhPtPhiLc[i]=0x0;
		fhPtPhiLcc[i]=0x0;
		fhPtPhiLcNonRc[i]=0x0;
		fhPtPhiLcL1520c[i]=0x0;
		fhPtPhiLcKstarc[i]=0x0;
		fhPtPhiLcDeltac[i]=0x0;
		fhPtPhiLcb[i]=0x0;
		fhPtPhiLcNonRb[i]=0x0;
		fhPtPhiLcL1520b[i]=0x0;
		fhPtPhiLcKstarb[i]=0x0;
		fhPtPhiLcDeltab[i]=0x0;

	}

	DefineOutput(1,TList::Class());  //My private output
	DefineOutput(2,TList::Class());
	DefineOutput(3,TH1F::Class());
	DefineOutput(4,AliNormalizationCounter::Class());
	if (fFillNtuple) {
		// Output slot #2 writes into a TNtuple container
		DefineOutput(5,TNtuple::Class());  //My private output
	}
}


//________________________________________________________________________
AliAnalysisTaskSELambdacTMVA::~AliAnalysisTaskSELambdacTMVA()
{
	//
	// Destructor
	//
	
	if (fOutput) {
		delete fOutput;
		fOutput = 0;
	}


	if (fVHF) {
		delete fVHF;
		fVHF = 0;
	}

	if(fRDCutsAnalysis){
		delete fRDCutsAnalysis;
		fRDCutsAnalysis = 0;
	}

	if (fListCuts) {
		delete fListCuts;
		fListCuts = 0;
	}
	if (fNentries){
		delete fNentries;
		fNentries = 0;
	}
	/*
		 if (fUtilPid){
		 delete fUtilPid;
		 fUtilPid = 0;
		 }
		 */
	if (fPIDResponse) {
		delete  fPIDResponse;
	}
	if(fCounter){
		delete fCounter;
		fCounter = 0;
	}

}  



//_________________________________________________________________
void AliAnalysisTaskSELambdacTMVA::Init()
{
	//
	// Initialization
	//

	if (fDebug > 1) printf("AnalysisTaskSELambdac::Init() \n");

	fListCuts=new TList();
	fListCuts->SetOwner();

	fListCuts->Add(new AliRDHFCutsLctopKpi(*fRDCutsAnalysis));
	PostData(2,fListCuts);
	return;
}

//________________________________________________________________________
void AliAnalysisTaskSELambdacTMVA::UserCreateOutputObjects()
{
	//
	// Create the output container
	//
	
	if (fDebug > 1) printf("AnalysisTaskSELambdac::UserCreateOutputObjects() \n");

	// Several histograms are more conveniently managed in a TList
	fOutput = new TList();
	fOutput->SetOwner();
	fOutput->SetName("OutputHistos");

	TString hisname,histitle;

	//Lc bit QA
	fhSelectBit			= new TH1F("hSelectBit","hSelectBit",5,-0.5,5.5);
	fhSelectBit->GetXaxis()->SetBinLabel(2,"All");
	fhSelectBit->GetXaxis()->SetBinLabel(3,"SelectionMap");
	fhSelectBit->GetXaxis()->SetBinLabel(4,"!LcCut");
	fhSelectBit->GetXaxis()->SetBinLabel(5,"!LcPID");
	fOutput->Add(fhSelectBit);

	fHistNEvents = new TH1F("fHistNEvents", "Number of processed events; ; Events",3,-1.5,1.5);
	fHistNEvents->Sumw2();
	fHistNEvents->SetMinimum(0);
	fOutput->Add(fHistNEvents);

	TString stepnames[11] = {"GeneratedLimAcc","Generated","GeneratedAcc","Reco3Prong","LcBit","IsSelectedTracks","IsInFidAcc","PtRange","IsSelectedCandidate","IsSelectedPID","IsSelectedNtuple"};
	for(Int_t i=0;i<11;i++) { // histograms for efficiency cross check
		//Lc vs Pt histograms
		hisname.Form("hNBkgNI%i",i);
		histitle.Form("N Bkg not injected %s",stepnames[i].Data());
		fhNBkgNI[i] = new TH1F(hisname.Data(),histitle.Data(),100,0,20);
		hisname.Form("hNLc%i",i);
		histitle.Form("N Lc %s",stepnames[i].Data());
		fhNLc[i] = new TH1F(hisname.Data(),histitle.Data(),100,0,20);
		hisname.Form("hNLcc%i",i);
		histitle.Form("N Lc from c %s",stepnames[i].Data());
		fhNLcc[i] = new TH1F(hisname.Data(),histitle.Data(),100,0,20);
		hisname.Form("hNLcNonRc%i",i);
		histitle.Form("N Lc non resonant from c %s",stepnames[i].Data());
		fhNLcNonRc[i] = new TH1F(hisname.Data(),histitle.Data(),100,0,20);
		hisname.Form("hNLcL1520c%i",i);
		histitle.Form("N Lc -> L(1520) + p from c %s",stepnames[i].Data());
		fhNLcL1520c[i] = new TH1F(hisname.Data(),histitle.Data(),100,0,20);
		hisname.Form("hNLcKstarc%i",i);
		histitle.Form("N Lc -> K* + pi from c %s",stepnames[i].Data());
		fhNLcKstarc[i] = new TH1F(hisname.Data(),histitle.Data(),100,0,20);
		hisname.Form("hNLcDeltac%i",i);
		histitle.Form("N Lc -> Delta + K from c %s",stepnames[i].Data());
		fhNLcDeltac[i] = new TH1F(hisname.Data(),histitle.Data(),100,0,20);
		hisname.Form("hNLcb%i",i);
		histitle.Form("N Lc from b %s",stepnames[i].Data());
		fhNLcb[i] = new TH1F(hisname.Data(),histitle.Data(),100,0,20);
		hisname.Form("hNLcNonRb%i",i);
		histitle.Form("N Lc non resonant from b %s",stepnames[i].Data());
		fhNLcNonRb[i] = new TH1F(hisname.Data(),histitle.Data(),100,0,20);
		hisname.Form("hNLcL1520b%i",i);
		histitle.Form("N Lc -> L(1520) + p from b %s",stepnames[i].Data());
		fhNLcL1520b[i] = new TH1F(hisname.Data(),histitle.Data(),100,0,20);
		hisname.Form("hNLcKstarb%i",i);
		histitle.Form("N Lc -> K* + pi from b %s",stepnames[i].Data());
		fhNLcKstarb[i] = new TH1F(hisname.Data(),histitle.Data(),100,0,20);
		hisname.Form("hNLcDeltab%i",i);
		histitle.Form("N Lc -> Delta + K from b %s",stepnames[i].Data());
		fhNLcDeltab[i] = new TH1F(hisname.Data(),histitle.Data(),100,0,20);
		fOutput->Add(fhNBkgNI[i]);
		fOutput->Add(fhNLc[i]);
		fOutput->Add(fhNLcc[i]);
		fOutput->Add(fhNLcNonRc[i]);
		fOutput->Add(fhNLcL1520c[i]);
		fOutput->Add(fhNLcKstarc[i]);
		fOutput->Add(fhNLcDeltac[i]);
		fOutput->Add(fhNLcb[i]);
		fOutput->Add(fhNLcNonRb[i]);
		fOutput->Add(fhNLcL1520b[i]);
		fOutput->Add(fhNLcKstarb[i]);
		fOutput->Add(fhNLcDeltab[i]);

		//Pt vs eta histograms
		hisname.Form("hPtEtaBkgNI%i",i);
		histitle.Form("Pt vs #eta Bkg not injected %s",stepnames[i].Data());
		fhPtEtaBkgNI[i] = new TH2F(hisname.Data(),histitle.Data(),20,0,20,200,-10,10);
		hisname.Form("hPtEtaLc%i",i);
		histitle.Form("Pt vs #eta Lc %s",stepnames[i].Data());
		fhPtEtaLc[i] = new TH2F(hisname.Data(),histitle.Data(),20,0,20,200,-10,10);
		hisname.Form("hPtEtaLcc%i",i);
		histitle.Form("Pt vs #eta Lc from c %s",stepnames[i].Data());
		fhPtEtaLcc[i] = new TH2F(hisname.Data(),histitle.Data(),20,0,20,200,-10,10);
		hisname.Form("hPtEtaLcNonRc%i",i);
		histitle.Form("Pt vs #eta Lc non resonant from c %s",stepnames[i].Data());
		fhPtEtaLcNonRc[i] = new TH2F(hisname.Data(),histitle.Data(),20,0,20,200,-10,10);
		hisname.Form("hPtEtaLcL1520c%i",i);
		histitle.Form("Pt vs #eta Lc -> L(1520) + p from c %s",stepnames[i].Data());
		fhPtEtaLcL1520c[i] = new TH2F(hisname.Data(),histitle.Data(),20,0,20,200,-10,10);
		hisname.Form("hPtEtaLcKstarc%i",i);
		histitle.Form("Pt vs #eta Lc -> K* + pi from c %s",stepnames[i].Data());
		fhPtEtaLcKstarc[i] = new TH2F(hisname.Data(),histitle.Data(),20,0,20,200,-10,10);
		hisname.Form("hPtEtaLcDeltac%i",i);
		histitle.Form("Pt vs #eta Lc -> Delta + K from c %s",stepnames[i].Data());
		fhPtEtaLcDeltac[i] = new TH2F(hisname.Data(),histitle.Data(),20,0,20,200,-10,10);
		hisname.Form("hPtEtaLcb%i",i);
		histitle.Form("Pt vs #eta Lc from b %s",stepnames[i].Data());
		fhPtEtaLcb[i] = new TH2F(hisname.Data(),histitle.Data(),20,0,20,200,-10,10);
		hisname.Form("hPtEtaLcNonRb%i",i);
		histitle.Form("Pt vs #eta Lc non resonant from b %s",stepnames[i].Data());
		fhPtEtaLcNonRb[i] = new TH2F(hisname.Data(),histitle.Data(),20,0,20,200,-10,10);
		hisname.Form("hPtEtaLcL1520b%i",i);
		histitle.Form("Pt vs #eta Lc -> L(1520) + p from b %s",stepnames[i].Data());
		fhPtEtaLcL1520b[i] = new TH2F(hisname.Data(),histitle.Data(),20,0,20,200,-10,10);
		hisname.Form("hPtEtaLcKstarb%i",i);
		histitle.Form("Pt vs #eta Lc -> K* + pi from b %s",stepnames[i].Data());
		fhPtEtaLcKstarb[i] = new TH2F(hisname.Data(),histitle.Data(),20,0,20,200,-10,10);
		hisname.Form("hPtEtaLcDeltab%i",i);
		histitle.Form("Pt vs #eta Lc -> Delta + K from b %s",stepnames[i].Data());
		fhPtEtaLcDeltab[i] = new TH2F(hisname.Data(),histitle.Data(),20,0,20,200,-10,10);
		fOutput->Add(fhPtEtaBkgNI[i]);
		fOutput->Add(fhPtEtaLc[i]);
		fOutput->Add(fhPtEtaLcc[i]);
		fOutput->Add(fhPtEtaLcNonRc[i]);
		fOutput->Add(fhPtEtaLcL1520c[i]);
		fOutput->Add(fhPtEtaLcKstarc[i]);
		fOutput->Add(fhPtEtaLcDeltac[i]);
		fOutput->Add(fhPtEtaLcb[i]);
		fOutput->Add(fhPtEtaLcNonRb[i]);
		fOutput->Add(fhPtEtaLcL1520b[i]);
		fOutput->Add(fhPtEtaLcKstarb[i]);
		fOutput->Add(fhPtEtaLcDeltab[i]);

		//Pt vs Y histograms
		hisname.Form("hPtYBkgNI%i",i);
		histitle.Form("Pt vs Y Bkg not injected %s",stepnames[i].Data());
		fhPtYBkgNI[i] = new TH2F(hisname.Data(),histitle.Data(),20,0,20,200,-10,10);
		hisname.Form("hPtYLc%i",i);
		histitle.Form("Pt vs Y Lc %s",stepnames[i].Data());
		fhPtYLc[i] = new TH2F(hisname.Data(),histitle.Data(),20,0,20,200,-10,10);
		hisname.Form("hPtYLcc%i",i);
		histitle.Form("Pt vs Y Lc from c %s",stepnames[i].Data());
		fhPtYLcc[i] = new TH2F(hisname.Data(),histitle.Data(),20,0,20,200,-10,10);
		hisname.Form("hPtYLcNonRc%i",i);
		histitle.Form("Pt vs Y Lc non resonant from c %s",stepnames[i].Data());
		fhPtYLcNonRc[i] = new TH2F(hisname.Data(),histitle.Data(),20,0,20,200,-10,10);
		hisname.Form("hPtYLcL1520c%i",i);
		histitle.Form("Pt vs Y Lc -> L(1520) + p from c %s",stepnames[i].Data());
		fhPtYLcL1520c[i] = new TH2F(hisname.Data(),histitle.Data(),20,0,20,200,-10,10);
		hisname.Form("hPtYLcKstarc%i",i);
		histitle.Form("Pt vs Y Lc -> K* + pi from c %s",stepnames[i].Data());
		fhPtYLcKstarc[i] = new TH2F(hisname.Data(),histitle.Data(),20,0,20,200,-10,10);
		hisname.Form("hPtYLcDeltac%i",i);
		histitle.Form("Pt vs Y Lc -> Delta + K from c %s",stepnames[i].Data());
		fhPtYLcDeltac[i] = new TH2F(hisname.Data(),histitle.Data(),20,0,20,200,-10,10);
		hisname.Form("hPtYLcb%i",i);
		histitle.Form("Pt vs Y Lc from b %s",stepnames[i].Data());
		fhPtYLcb[i] = new TH2F(hisname.Data(),histitle.Data(),20,0,20,200,-10,10);
		hisname.Form("hPtYLcNonRb%i",i);
		histitle.Form("Pt vs Y Lc non resonant from b %s",stepnames[i].Data());
		fhPtYLcNonRb[i] = new TH2F(hisname.Data(),histitle.Data(),20,0,20,200,-10,10);
		hisname.Form("hPtYLcL1520b%i",i);
		histitle.Form("Pt vs Y Lc -> L(1520) + p from b %s",stepnames[i].Data());
		fhPtYLcL1520b[i] = new TH2F(hisname.Data(),histitle.Data(),20,0,20,200,-10,10);
		hisname.Form("hPtYLcKstarb%i",i);
		histitle.Form("Pt vs Y Lc -> K* + pi from b %s",stepnames[i].Data());
		fhPtYLcKstarb[i] = new TH2F(hisname.Data(),histitle.Data(),20,0,20,200,-10,10);
		hisname.Form("hPtYLcDeltab%i",i);
		histitle.Form("Pt vs Y Lc -> Delta + K from b %s",stepnames[i].Data());
		fhPtYLcDeltab[i] = new TH2F(hisname.Data(),histitle.Data(),20,0,20,200,-10,10);
		fOutput->Add(fhPtYBkgNI[i]);
		fOutput->Add(fhPtYLc[i]);
		fOutput->Add(fhPtYLcc[i]);
		fOutput->Add(fhPtYLcNonRc[i]);
		fOutput->Add(fhPtYLcL1520c[i]);
		fOutput->Add(fhPtYLcKstarc[i]);
		fOutput->Add(fhPtYLcDeltac[i]);
		fOutput->Add(fhPtYLcb[i]);
		fOutput->Add(fhPtYLcNonRb[i]);
		fOutput->Add(fhPtYLcL1520b[i]);
		fOutput->Add(fhPtYLcKstarb[i]);
		fOutput->Add(fhPtYLcDeltab[i]);

		//Pt vs phi histograms
		hisname.Form("hPtPhiBkgNI%i",i);
		histitle.Form("Pt vs #phi Bkg not injected %s",stepnames[i].Data());
		fhPtPhiBkgNI[i] = new TH2F(hisname.Data(),histitle.Data(),20,0,20,70,0,7);
		hisname.Form("hPtPhiLc%i",i);
		histitle.Form("Pt vs #phi Lc %s",stepnames[i].Data());
		fhPtPhiLc[i] = new TH2F(hisname.Data(),histitle.Data(),20,0,20,70,0,7);
		hisname.Form("hPtPhiLcc%i",i);
		histitle.Form("Pt vs #phi Lc from c %s",stepnames[i].Data());
		fhPtPhiLcc[i] = new TH2F(hisname.Data(),histitle.Data(),20,0,20,70,0,7);
		hisname.Form("hPtPhiLcNonRc%i",i);
		histitle.Form("Pt vs #phi Lc non resonant from c %s",stepnames[i].Data());
		fhPtPhiLcNonRc[i] = new TH2F(hisname.Data(),histitle.Data(),20,0,20,70,0,7);
		hisname.Form("hPtPhiLcL1520c%i",i);
		histitle.Form("Pt vs #phi Lc -> L(1520) + p from c %s",stepnames[i].Data());
		fhPtPhiLcL1520c[i] = new TH2F(hisname.Data(),histitle.Data(),20,0,20,70,0,7);
		hisname.Form("hPtPhiLcKstarc%i",i);
		histitle.Form("Pt vs #phi Lc -> K* + pi from c %s",stepnames[i].Data());
		fhPtPhiLcKstarc[i] = new TH2F(hisname.Data(),histitle.Data(),20,0,20,70,0,7);
		hisname.Form("hPtPhiLcDeltac%i",i);
		histitle.Form("Pt vs #phi Lc -> Delta + K from c %s",stepnames[i].Data());
		fhPtPhiLcDeltac[i] = new TH2F(hisname.Data(),histitle.Data(),20,0,20,70,0,7);
		hisname.Form("hPtPhiLcb%i",i);
		histitle.Form("Pt vs #phi Lc from b %s",stepnames[i].Data());
		fhPtPhiLcb[i] = new TH2F(hisname.Data(),histitle.Data(),20,0,20,70,0,7);
		hisname.Form("hPtPhiLcNonRb%i",i);
		histitle.Form("Pt vs #phi Lc non resonant from b %s",stepnames[i].Data());
		fhPtPhiLcNonRb[i] = new TH2F(hisname.Data(),histitle.Data(),20,0,20,70,0,7);
		hisname.Form("hPtPhiLcL1520b%i",i);
		histitle.Form("Pt vs #phi Lc -> L(1520) + p from b %s",stepnames[i].Data());
		fhPtPhiLcL1520b[i] = new TH2F(hisname.Data(),histitle.Data(),20,0,20,70,0,7);
		hisname.Form("hPtPhiLcKstarb%i",i);
		histitle.Form("Pt vs #phi Lc -> K* + pi from b %s",stepnames[i].Data());
		fhPtPhiLcKstarb[i] = new TH2F(hisname.Data(),histitle.Data(),20,0,20,70,0,7);
		hisname.Form("hPtPhiLcDeltab%i",i);
		histitle.Form("Pt vs #phi Lc -> Delta + K from b %s",stepnames[i].Data());
		fhPtPhiLcDeltab[i] = new TH2F(hisname.Data(),histitle.Data(),20,0,20,70,0,7);
		fOutput->Add(fhPtPhiBkgNI[i]);
		fOutput->Add(fhPtPhiLc[i]);
		fOutput->Add(fhPtPhiLcc[i]);
		fOutput->Add(fhPtPhiLcNonRc[i]);
		fOutput->Add(fhPtPhiLcL1520c[i]);
		fOutput->Add(fhPtPhiLcKstarc[i]);
		fOutput->Add(fhPtPhiLcDeltac[i]);
		fOutput->Add(fhPtPhiLcb[i]);
		fOutput->Add(fhPtPhiLcNonRb[i]);
		fOutput->Add(fhPtPhiLcL1520b[i]);
		fOutput->Add(fhPtPhiLcKstarb[i]);
		fOutput->Add(fhPtPhiLcDeltab[i]);
	}

	// fhChi2 = new TH1F("fhChi2", "Chi2",100,0.,10.);
	// fhChi2->Sumw2();
	// fOutput->Add(fhChi2);

	fNentries=new TH1F("fNentries", "n Events/Candidates QA", 16,0.5,16.5);

	//Event and candidate QA - entries at each step
	fNentries->GetXaxis()->SetBinLabel(1,"nEventsRejTM");
	fNentries->GetXaxis()->SetBinLabel(2,"nEventsNoVtx");
	fNentries->GetXaxis()->SetBinLabel(3,"nEventsRejCutPileup");
	fNentries->GetXaxis()->SetBinLabel(4,"nLcGen");
	fNentries->GetXaxis()->SetBinLabel(5,"nLcGenFidAcc");
	fNentries->GetXaxis()->SetBinLabel(6,"nCandReco3Prong");
	fNentries->GetXaxis()->SetBinLabel(7,"nCandLcBit");
	fNentries->GetXaxis()->SetBinLabel(8,"nCandIsSelTracks");
	fNentries->GetXaxis()->SetBinLabel(9,"nCandIsInFidAcc");
	fNentries->GetXaxis()->SetBinLabel(10,"ptbin=-1");
	fNentries->GetXaxis()->SetBinLabel(11,"nCandIsSelCand");
	fNentries->GetXaxis()->SetBinLabel(12,"PID=0");
	fNentries->GetXaxis()->SetBinLabel(13,"PID=1");
	fNentries->GetXaxis()->SetBinLabel(14,"PID=2");
	fNentries->GetXaxis()->SetBinLabel(15,"PID=3");
	fNentries->GetXaxis()->SetBinLabel(16,"nLcSelected");
	fNentries->GetXaxis()->SetNdivisions(1,kFALSE);

	AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
	AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
	fPIDResponse = inputHandler->GetPIDResponse();

	if(fRDCutsAnalysis->GetIsUsePID()){
		fRDCutsAnalysis->GetPidHF()->SetPidResponse(fPIDResponse);
		fRDCutsAnalysis->GetPidpion()->SetPidResponse(fPIDResponse);
		fRDCutsAnalysis->GetPidprot()->SetPidResponse(fPIDResponse);
		fRDCutsAnalysis->GetPidHF()->SetOldPid(kFALSE);
		fRDCutsAnalysis->GetPidpion()->SetOldPid(kFALSE);
		fRDCutsAnalysis->GetPidprot()->SetOldPid(kFALSE);
	}


	PostData(1,fOutput);
	PostData(3,fNentries);

	TString normName="NormalizationCounter";
	AliAnalysisDataContainer *cont = GetOutputSlot(4)->GetContainer();
	if(cont)normName=(TString)cont->GetName(); 
	fCounter = new AliNormalizationCounter(normName.Data());
	fCounter->Init();
	PostData(4,fCounter);
	if (fFillNtuple) {
		//basic ntuple
		if(fFillNtuple==1)       fNtupleLambdac = new TNtuple("fNtupleLambdac", "Lc", "isLcBkg:InvMasspKpi:InvMasspiKp:PtTr0:PtTr1:PtTr2:PtLc:CosP:DecayL:DecayLSig:Dist12:SigVert:DCA:DecayLXY:DecayLXYSig:isLcResonant:selectionPID:Tr0Ppi:Tr0PK:Tr0Pp:Tr1Ppi:Tr1PK:Tr1Pp:Tr2Ppi:Tr2PK:Tr2Pp");
		//NSigma PID
		else if(fFillNtuple==2)  fNtupleLambdac = new TNtuple("fNtupleLambdac", "Lc", "isLcBkg:InvMasspKpi:InvMasspiKp:PtTr0:PtTr1:PtTr2:PtLc:CosP:DecayL:DecayLSig:Dist12:SigVert:DCA:DecayLXY:DecayLXYSig:isLcResonant:selectionPID:Tr0Ppi:Tr0PK:Tr0Pp:Tr1Ppi:Tr1PK:Tr1Pp:Tr2Ppi:Tr2PK:Tr2Pp:Tr0NSigmapi:Tr0NSigmaK:Tr0NSigmap:Tr1NSigmapi:Tr1NSigmaK:Tr1NSigmap:Tr2NSigmapi:Tr2NSigmaK:Tr2NSigmap");
		//2 prong decay products
		else if(fFillNtuple==3)  fNtupleLambdac = new TNtuple("fNtupleLambdac", "Lc", "isLcBkg:InvMasspKpi:InvMasspiKp:PtTr0:PtTr1:PtTr2:PtLc:CosP:DecayL:DecayLSig:Dist12:SigVert:DCA:DecayLXY:DecayLXYSig:isLcResonant:selectionPID:Tr0Ppi:Tr0PK:Tr0Pp:Tr1Ppi:Tr1PK:Tr1Pp:Tr2Ppi:Tr2PK:Tr2Pp:Tr0NSigmapi:Tr0NSigmaK:Tr0NSigmap:Tr1NSigmapi:Tr1NSigmaK:Tr1NSigmap:Tr2NSigmapi:Tr2NSigmaK:Tr2NSigmap:InvMasspK:InvMassKpi:InvMassppi:InvMassKp:InvMasspiK:InvMasspip");
		else AliFatal("Invalid fill ntuple argument");
		PostData(5,fNtupleLambdac);
	}

	return;
}

//________________________________________________________________________
void AliAnalysisTaskSELambdacTMVA::UserExec(Option_t */*option*/)
{
	//
	// Execute analysis for current event:
	// heavy flavor candidates association to MC truth
	//

	AliAODEvent *aod = dynamic_cast<AliAODEvent*> (InputEvent());
	//tmp
	fHistNEvents->Fill(0); // count event
	// Post the data already here


	TClonesArray *array3Prong = 0;
	TClonesArray *arrayLikeSign =0;
	if(!aod && AODEvent() && IsStandardAOD()) {
		// In case there is an AOD handler writing a standard AOD, use the AOD 
		// event in memory rather than the input (ESD) event.    
		aod = dynamic_cast<AliAODEvent*> (AODEvent());
		// in this case the braches in the deltaAOD (AliAOD.VertexingHF.root)
		// have to taken from the AOD event hold by the AliAODExtension
		AliAODHandler* aodHandler = (AliAODHandler*) 
			((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());
		if(aodHandler->GetExtensions()) {
			AliAODExtension *ext = (AliAODExtension*)aodHandler->GetExtensions()->FindObject("AliAOD.VertexingHF.root");
			AliAODEvent *aodFromExt = ext->GetAOD();
			array3Prong=(TClonesArray*)aodFromExt->GetList()->FindObject("Charm3Prong");
			arrayLikeSign=(TClonesArray*)aodFromExt->GetList()->FindObject("LikeSign3Prong");
		}
	} else if(aod) {
		array3Prong=(TClonesArray*)aod->GetList()->FindObject("Charm3Prong");
		arrayLikeSign=(TClonesArray*)aod->GetList()->FindObject("LikeSign3Prong");
	}

	if(!aod) return;

	TClonesArray *arrayMC=0;
	AliAODMCHeader *mcHeader=0;

	// load MC particles
	if(fReadMC){

		arrayMC =  (TClonesArray*)aod->GetList()->FindObject(AliAODMCParticle::StdBranchName());
		if(!arrayMC) {
			AliError("AliAnalysisTaskSELambdacTMVA::UserExec: MC particles branch not found!\n");
			return;
		}


		// load MC header
		mcHeader =  (AliAODMCHeader*)aod->GetList()->FindObject(AliAODMCHeader::StdBranchName());
		if(!mcHeader) {
			AliError("AliAnalysisTaskSELambdacTMVA::UserExec: MC header branch not found!\n");
			return;
		}
	}


	if(!array3Prong || !aod) {
		AliError("AliAnalysisTaskSELambdacTMVA::UserExec: Charm3Prong branch not found!\n");
		return;
	}
	if(!arrayLikeSign) {
		AliDebug(2,"AliAnalysisTaskSELambdacTMVA::UserExec: LikeSign3Prong branch not found!\n");
		//  return;
	}

	//Trigger mask = 0 rejection for LHC13d3
	Int_t runnumber = aod->GetRunNumber();
	if (aod->GetTriggerMask() == 0 && (runnumber >= 195344 && runnumber <= 195677)){
		Int_t nentriesTM = arrayMC->GetEntriesFast(); 
		AliDebug(2,Form("Event rejected because of null trigger mask, n entries = %i",nentriesTM));
		fNentries->Fill(1);
		return;
	}

	// fix for temporary bug in ESDfilter
	// the AODs with null vertex pointer didn't pass the PhysSel
	AliAODVertex *vtx1 = (AliAODVertex*)aod->GetPrimaryVertex();
	if(!vtx1 || TMath::Abs(aod->GetMagneticField())<0.001) {
		fNentries->Fill(2);
		return;
	}
	fCounter->StoreEvent(aod,fRDCutsAnalysis,fReadMC);
	TString trigclass=aod->GetFiredTriggerClasses();
	//if(trigclass.Contains("C0SMH-B-NOPF-ALLNOTRD") || trigclass.Contains("C0SMH-B-NOPF-ALL")) 
	Bool_t isEvSelAnCuts=fRDCutsAnalysis->IsEventSelected(aod);
	if(!isEvSelAnCuts){ //reject if not selected with analysis cut
		if(fRDCutsAnalysis->GetWhyRejection()==1) // rejected for pileup
			fNentries->Fill(3);
		return;
	}


	//Bool_t isThereA3prongWithGoodTracks = kFALSE;
	//Bool_t isThereA3ProngLcKine = kFALSE;
	//Bool_t isThereA3ProngLcKineANDpid = kFALSE;
	//Bool_t isThereA3ProngLcMC = kFALSE;
	//Bool_t isThereA3ProngCyes = kFALSE;
	//Bool_t isThereA3ProngByes = kFALSE;
	//Bool_t isThereA3ProngJPsi = kFALSE;

	Int_t n3Prong = array3Prong->GetEntriesFast();
	Int_t nSelectedloose[1]={0};
	Int_t nSelectedtight[1]={0};

	//MC Generated level
	//loop over MC particles to find Lc generated
	Int_t pdgcode = 0;
	Int_t pdgcodemother = 0;
	Bool_t isInFidAcc = kFALSE; 
	Bool_t isInAcc = kFALSE;
	if(fReadMC) {
		for (Int_t iPart=0; iPart<arrayMC->GetEntriesFast(); iPart++) { 
			fIsLc=0;
			fIsLcResonant=0;
			AliAODMCParticle* mcPart = dynamic_cast<AliAODMCParticle*>(arrayMC->At(iPart));
			if (!mcPart){
				AliError("Failed casting particle from MC array!, Skipping particle");
				continue;
			}
			isInFidAcc = fRDCutsAnalysis->IsInFiducialAcceptance(mcPart->Pt(),mcPart->Y());
			fCandidateVars[0] = mcPart->Pt();
			fCandidateVars[1] = mcPart->Eta();
			fCandidateVars[2] = mcPart->Y();
			fCandidateVars[3] = mcPart->Phi();
			//Check whether particle is Lc
			pdgcode = TMath::Abs(mcPart->GetPdgCode());
			if(pdgcode==4122) {
				AliDebug(2,"Found Lc! now check mother");
				fIsLc=1;
				AliVertexingHFUtils *util = new AliVertexingHFUtils();
				Int_t imother=mcPart->GetMother();
				if(imother>0) { //found mother 
					AliAODMCParticle* mcPartMother = dynamic_cast<AliAODMCParticle*>(arrayMC->At(imother));
					if(!mcPart){
						AliError("Failed casting mother particle, Skipping particle");
						continue;
					}
					pdgcodemother = TMath::Abs(mcPartMother->GetPdgCode());
					if(pdgcodemother == 5122){
						AliDebug(2,"Lc comes from Lb");
						fIsLc=2;
					}
				}
				//Check daughters
				SetLambdacDaugh(mcPart,arrayMC,isInAcc);
				if(fIsLcResonant>=1) {
					AliDebug(2,"Lc has p K pi in final state");
					fNentries->Fill(4);
					if(isInFidAcc){
						fNentries->Fill(5);
						if((TMath::Abs(mcPart->Y()) < 0.5)) {//limited acceptance
							AliDebug(2,"Lc in limited acceptance");
							FillEffHists(kGeneratedLimAcc);
						}
						FillEffHists(kGenerated); //MC generated ---> Should be within fiducial acceptance
						if(isInAcc) FillEffHists(kGeneratedAcc); //all daughters in acceptance
					}
				}
				else if(fIsLcResonant==0) AliDebug(2,"no p K pi final state");
				else AliError(Form("Not pKpi or background - should not happen! fIsLcResonant = %i",fIsLcResonant));
				delete util;
			}
		}
	}
	

	//Reconstruction level
	//loop over 3 prong candidates
	for (Int_t i3Prong = 0; i3Prong < n3Prong; i3Prong++) {
		AliAODRecoDecayHF3Prong *d = (AliAODRecoDecayHF3Prong*)array3Prong->UncheckedAt(i3Prong);
		//(AliAODEvent *aod,AliAODRecoDecayHF3Prong *part,TClonesArray *arrayMC)
		SetIsLc(d,arrayMC);
		fCandidateVars[0] = d->Pt();
		fCandidateVars[1] = d->Eta();
		fCandidateVars[2] = d->Y(4122);
		fCandidateVars[3] = d->Phi();

		FillEffHists(kReco3Prong);
		fNentries->Fill(6);

		Bool_t unsetvtx=kFALSE;
		if(!d->GetOwnPrimaryVtx()){
			d->SetOwnPrimaryVtx(vtx1);
			unsetvtx=kTRUE;
		}

		//Filter bit selection and QA:
		fhSelectBit->Fill(1);
		if(d->GetSelectionMap()){
			fhSelectBit->Fill(2);
			if(!d->HasSelectionBit(AliRDHFCuts::kLcCuts))  fhSelectBit->Fill(3);
			if(!d->HasSelectionBit(AliRDHFCuts::kLcPID))   fhSelectBit->Fill(4);
			if(fLcCut&&!d->HasSelectionBit(AliRDHFCuts::kLcCuts))		continue;
			if(fLcPIDCut&&!d->HasSelectionBit(AliRDHFCuts::kLcPID))		continue;
		}
		FillEffHists(kLcBit);
		fNentries->Fill(7);

		//track selection
		Int_t isSelectedTracks = fRDCutsAnalysis->IsSelected(d,AliRDHFCuts::kTracks,aod);
		if(!isSelectedTracks) continue;
		FillEffHists(kIsSelectedTracks);
		fNentries->Fill(8);

		//isThereA3prongWithGoodTracks=kTRUE;

		if (fRDCutsAnalysis->IsInFiducialAcceptance(d->Pt(),d->Y(4122))) {fNentries->Fill(9);}else{continue;}
		FillEffHists(kIsInFidAcc);

		Int_t ptbin=fRDCutsAnalysis->PtBin(d->Pt());
		if(ptbin==-1) {fNentries->Fill(10); continue;} //out of bounds
		FillEffHists(kPtRange);

		//idproton, pion using isSelectedPID
		//Bool_t isPID=fRDCutsAnalysis->GetIsUsePID();
		Int_t isSelectedPID=fRDCutsAnalysis->IsSelected(d,AliRDHFCuts::kPID,aod);
		if(isSelectedPID==0 ) fNentries->Fill(12);
		else if(isSelectedPID==1) fNentries->Fill(13);
		else if(isSelectedPID==2) fNentries->Fill(14);
		else fNentries->Fill(15);
		if(isSelectedPID>0) FillEffHists(kIsSelectedPID);

		Int_t selection=fRDCutsAnalysis->IsSelected(d,AliRDHFCuts::kCandidate,aod);
		if(!selection) continue;
		FillEffHists(kIsSelectedCandidate);
		fNentries->Fill(11);

		if(fFillNtuple) FillNtuple(aod,d,arrayMC,selection);
		FillEffHists(kIsSelectedNtuple);
		if(fIsLc>=1 && fIsLc <= 2) fNentries->Fill(16);
		if(unsetvtx) d->UnsetOwnPrimaryVtx();
	}
	fCounter->StoreCandidates(aod,nSelectedloose[0],kTRUE);
	fCounter->StoreCandidates(aod,nSelectedtight[0],kFALSE);


	PostData(1,fOutput); 
	PostData(3,fNentries);
	PostData(4,fCounter);

	return;
}



//________________________________________________________________________
void AliAnalysisTaskSELambdacTMVA::Terminate(Option_t */*option*/)
{
	//
	// Terminate analysis
	//
	
	if (fDebug > 1) printf("AnalysisTaskSELambdac: Terminate() \n");

	fOutput = dynamic_cast<TList*> (GetOutputData(1));
	if (!fOutput) {     
		AliError("ERROR: fOutput not available\n");
		return;
	}
	//fHistNEvents = dynamic_cast<TH1F*>(fOutput->FindObject("fHistNEvents"));

	if(fFillNtuple){
		fNtupleLambdac = dynamic_cast<TNtuple*>(GetOutputData(5));
	}


	return;
}

//________________________________________________________________________
Int_t AliAnalysisTaskSELambdacTMVA::MatchToMCLambdac(AliAODRecoDecayHF3Prong *d,TClonesArray *arrayMC) const{

	//
	// check if the candidate is a Lambdac decaying in pKpi or in the resonant channels
	//
	
	Int_t lambdacLab[3]={0,0,0};
	//  Int_t pdgs[3]={0,0,0};
	for(Int_t i=0;i<3;i++){
		AliAODTrack *daugh=(AliAODTrack*)d->GetDaughter(i);
		Int_t lab=daugh->GetLabel();
		if(lab<0) return 0;
		AliAODMCParticle *part= (AliAODMCParticle*)arrayMC->At(lab);
		if(!part) continue;
		//    pdgs[i]=part->GetPdgCode();
		Int_t partPdgcode = TMath::Abs(part->GetPdgCode());
		if(partPdgcode==211 || partPdgcode==321 || partPdgcode==2212){
			Int_t motherLabel=part->GetMother();
			if(motherLabel<0) return 0;
			AliAODMCParticle *motherPart = (AliAODMCParticle*)arrayMC->At(motherLabel);
			if(!motherPart) continue;
			Int_t motherPdg = TMath::Abs(motherPart->GetPdgCode());
			if(motherPdg==4122) {
				if(GetLambdacDaugh(motherPart,arrayMC)){lambdacLab[i]=motherLabel;continue;}
			}
			if(motherPdg==313 || motherPdg==2224 || motherPdg==3124){
				Int_t granMotherLabel=motherPart->GetMother();
				if(granMotherLabel<0) return 0;
				AliAODMCParticle *granMotherPart = (AliAODMCParticle*)arrayMC->At(granMotherLabel);
				if(!granMotherPart) continue;
				Int_t granMotherPdg  = TMath::Abs(granMotherPart->GetPdgCode());
				if(granMotherPdg ==4122) {
					if(GetLambdacDaugh(granMotherPart,arrayMC)) {lambdacLab[i]=granMotherLabel;continue;}
				}
			}
		}
	}

	if(lambdacLab[0]==lambdacLab[1] && lambdacLab[1]==lambdacLab[2]) {return lambdacLab[0];}
	return 0;

}

//-----------------------------

Int_t AliAnalysisTaskSELambdacTMVA::LambdacDaugh(AliAODMCParticle *part, TClonesArray *arrayMC, Bool_t &IsInAcc) const {
	
	// 
	// return value of Lc resonant channel, from AOD MC particle
	// Also check whether Lc daughters are in acceptance
	// 0=not Lc, 1=non resonant Lc, 2=via L(1520) + pi, 3=via K* + p, 4=via Delta++ + K  
	//

	Int_t numberOfLambdac=0;
	IsInAcc=kTRUE;
	if(TMath::Abs(part->GetPdgCode())!=4122) return 0;
	// Int_t daughTmp[2];
	// daughTmp[0]=part->GetDaughter(0);
	// daughTmp[1]=part->GetDaughter(1);
	Int_t nDaugh = (Int_t)part->GetNDaughters();
	if(nDaugh<2) return 0;
	if(nDaugh>3) return 0;
	AliAODMCParticle* pdaugh1 = (AliAODMCParticle*)arrayMC->At(part->GetDaughter(0));
	if(!pdaugh1) {return 0;}
	Int_t number1 = TMath::Abs(pdaugh1->GetPdgCode());
	AliAODMCParticle* pdaugh2 = (AliAODMCParticle*)arrayMC->At(part->GetDaughter(1));
	if(!pdaugh2) {return 0;}
	Int_t number2 = TMath::Abs(pdaugh2->GetPdgCode());

	AliDebug(2,"Is non resonant?");
	if(nDaugh==3){
		Int_t thirdDaugh=part->GetDaughter(1)-1;
		AliAODMCParticle* pdaugh3 = (AliAODMCParticle*)arrayMC->At(thirdDaugh);
		Int_t number3 = TMath::Abs(pdaugh3->GetPdgCode());
		if((number1==321 && number2==211 && number3==2212) ||
				(number1==211 && number2==321 && number3==2212) ||
				(number1==211 && number2==2212 && number3==321) ||
				(number1==321 && number2==2212 && number3==211) ||
				(number1==2212 && number2==321 && number3==211) ||
				(number1==2212 && number2==211 && number3==321)) {
			numberOfLambdac++;
			if(    TMath::Abs(pdaugh1->Eta()) > 0.9 || pdaugh1->Pt() < 0.1
					|| TMath::Abs(pdaugh2->Eta()) > 0.9 || pdaugh2->Pt() < 0.1
					|| TMath::Abs(pdaugh3->Eta()) > 0.9 || pdaugh3->Pt() < 0.1) IsInAcc=kFALSE;
			AliDebug(2,"Lc decays non-resonantly");
			return 1;
		} 
	}

	if(nDaugh==2){

		//Lambda resonant

		//Lambda -> p K*0
		//
		Int_t nfiglieK=0;

		if((number1==2212 && number2==313)){
			nfiglieK=pdaugh2->GetNDaughters();
			if(nfiglieK!=2) return 0;
			AliAODMCParticle* pdaughK1 = (AliAODMCParticle*)arrayMC->At(pdaugh2->GetDaughter(0));
			AliAODMCParticle* pdaughK2 = (AliAODMCParticle*)arrayMC->At(pdaugh2->GetDaughter(1));
			if(!pdaughK1) return 0;
			if(!pdaughK2) return 0;
			if((TMath::Abs(pdaughK1->GetPdgCode())==211 && TMath::Abs(pdaughK2->GetPdgCode())==321) || (TMath::Abs(pdaughK1->GetPdgCode())==321 && TMath::Abs(pdaughK2->GetPdgCode())==211)) {
				numberOfLambdac++;
			if(    TMath::Abs(pdaugh1->Eta()) > 0.9 || pdaugh1->Pt() < 0.1
					|| TMath::Abs(pdaughK1->Eta()) > 0.9 || pdaughK1->Pt() < 0.1
					|| TMath::Abs(pdaughK2->Eta()) > 0.9 || pdaughK2->Pt() < 0.1) IsInAcc=kFALSE;
				AliDebug(2,"Lc decays via K* p");
				return 3;
			}
		}

		if((number1==313 && number2==2212)){
			nfiglieK=pdaugh1->GetNDaughters();
			if(nfiglieK!=2) return 0;
			AliAODMCParticle* pdaughK1 = (AliAODMCParticle*)arrayMC->At(pdaugh1->GetDaughter(0));
			AliAODMCParticle* pdaughK2 = (AliAODMCParticle*)arrayMC->At(pdaugh1->GetDaughter(1));
			if(!pdaughK1) return 0;
			if(!pdaughK2) return 0;
			if((TMath::Abs(pdaughK1->GetPdgCode())==211 && TMath::Abs(pdaughK2->GetPdgCode())==321) || (TMath::Abs(pdaughK1->GetPdgCode())==321 && TMath::Abs(pdaughK2->GetPdgCode())==211)) {
				numberOfLambdac++;
			if(    TMath::Abs(pdaugh2->Eta()) > 0.9 || pdaugh2->Pt() < 0.1
					|| TMath::Abs(pdaughK1->Eta()) > 0.9 || pdaughK1->Pt() < 0.1
					|| TMath::Abs(pdaughK2->Eta()) > 0.9 || pdaughK2->Pt() < 0.1) IsInAcc=kFALSE;
				AliDebug(2,"Lc decays via K* p");
				return 3;
			}
		}

		//Lambda -> Delta++ k
		Int_t nfiglieDelta=0;
		if(number1==321 && number2==2224){
			nfiglieDelta=pdaugh2->GetNDaughters();
			if(nfiglieDelta!=2) return 0;
			AliAODMCParticle *pdaughD1=(AliAODMCParticle*)arrayMC->At(pdaugh2->GetDaughter(0));
			AliAODMCParticle *pdaughD2=(AliAODMCParticle*)arrayMC->At(pdaugh2->GetDaughter(1));
			if(!pdaughD1) return 0;
			if(!pdaughD2) return 0;
			if((TMath::Abs(pdaughD1->GetPdgCode())==211 && TMath::Abs(pdaughD2->GetPdgCode())==2212) || (TMath::Abs(pdaughD1->GetPdgCode())==2212 && TMath::Abs(pdaughD2->GetPdgCode())==211)) {
				numberOfLambdac++;
			if(    TMath::Abs(pdaugh1->Eta()) > 0.9 || pdaugh1->Pt() < 0.1
					|| TMath::Abs(pdaughD1->Eta()) > 0.9 || pdaughD1->Pt() < 0.1
					|| TMath::Abs(pdaughD2->Eta()) > 0.9 || pdaughD2->Pt() < 0.1) IsInAcc=kFALSE;
				AliDebug(2,"Lc decays via Delta++ k");
				return 4;
			}
		}
		if(number1==2224 && number2==321){
			nfiglieDelta=pdaugh1->GetNDaughters();
			if(nfiglieDelta!=2) return 0;
			AliAODMCParticle* pdaughD1 = (AliAODMCParticle*)arrayMC->At(pdaugh1->GetDaughter(0));
			AliAODMCParticle* pdaughD2 = (AliAODMCParticle*)arrayMC->At(pdaugh1->GetDaughter(1));
			if(!pdaughD1) return 0;
			if(!pdaughD2) return 0;
			if((TMath::Abs(pdaughD1->GetPdgCode())==211 && TMath::Abs(pdaughD2->GetPdgCode())==2212) || (TMath::Abs(pdaughD1->GetPdgCode())==2212 && TMath::Abs(pdaughD2->GetPdgCode())==211)) {
				numberOfLambdac++;
			if(    TMath::Abs(pdaugh2->Eta()) > 0.9 || pdaugh2->Pt() < 0.1
					|| TMath::Abs(pdaughD1->Eta()) > 0.9 || pdaughD1->Pt() < 0.1
					|| TMath::Abs(pdaughD2->Eta()) > 0.9 || pdaughD2->Pt() < 0.1) IsInAcc=kFALSE;
				AliDebug(2,"Lc decays via Delta++ k");
				return 4;
			}
		}


		//Lambdac -> Lambda(1520) pi
		Int_t nfiglieLa=0;
		if(number1==3124 && number2==211){
			nfiglieLa=pdaugh1->GetNDaughters();
			if(nfiglieLa!=2) return 0;
			AliAODMCParticle *pdaughL1=(AliAODMCParticle*)arrayMC->At(pdaugh1->GetDaughter(0));
			AliAODMCParticle *pdaughL2=(AliAODMCParticle*)arrayMC->At(pdaugh1->GetDaughter(1));
			if(!pdaughL1) return 0;
			if(!pdaughL2) return 0;
			if((TMath::Abs(pdaughL1->GetPdgCode())==321 && TMath::Abs(pdaughL2->GetPdgCode())==2212) || (TMath::Abs(pdaughL1->GetPdgCode())==2212 && TMath::Abs(pdaughL2->GetPdgCode())==321)) {
				numberOfLambdac++;
				if(    TMath::Abs(pdaugh2->Eta()) > 0.9 || pdaugh2->Pt() < 0.1
					|| TMath::Abs(pdaughL1->Eta()) > 0.9 || pdaughL1->Pt() < 0.1
					|| TMath::Abs(pdaughL2->Eta()) > 0.9 || pdaughL2->Pt() < 0.1) IsInAcc=kFALSE;
				AliDebug(2,"Lc decays via Lambda(1520) pi");
				return 2;
			}
		}
		if(number1==211 && number2==3124){
			nfiglieLa=pdaugh2->GetNDaughters();
			if(nfiglieLa!=2) return 0;
			AliAODMCParticle *pdaughL1=(AliAODMCParticle*)arrayMC->At(pdaugh2->GetDaughter(0));
			AliAODMCParticle *pdaughL2=(AliAODMCParticle*)arrayMC->At(pdaugh2->GetDaughter(1));
			if(!pdaughL1) return 0;
			if(!pdaughL2) return 0;
			if((TMath::Abs(pdaughL1->GetPdgCode())==321 && TMath::Abs(pdaughL2->GetPdgCode())==2212) || (TMath::Abs(pdaughL1->GetPdgCode())==2212 && TMath::Abs(pdaughL2->GetPdgCode())==321)) {
				numberOfLambdac++;
				if(    TMath::Abs(pdaugh1->Eta()) > 0.9 || pdaugh1->Pt() < 0.1
					|| TMath::Abs(pdaughL1->Eta()) > 0.9 || pdaughL1->Pt() < 0.1
					|| TMath::Abs(pdaughL2->Eta()) > 0.9 || pdaughL2->Pt() < 0.1) IsInAcc=kFALSE;
				AliDebug(2,"Lc decays via Lambda(1520) pi");
				return 2;
			}
		}
	}

	if(numberOfLambdac>0) {return -100; AliDebug(2,"Lc decays via one of 4 resonances!");}
	return 0;
}

//-----------------------------

void AliAnalysisTaskSELambdacTMVA::SetIsLc(AliAODRecoDecayHF3Prong *part,
		TClonesArray *arrayMC) {

	//
	// function that sets fIsLc and fIsLcResonant, from reconstructed 3 prong decay
	// fIsLc - 0 = not Lc,  1 = Lc from c, 2 = Lc from b
	// fIsLcResonant - 1= Non Resonant, 2=Lc->L1520+p, 3=Lc->K*+pi, 4=Lc->Delta+K
	//

	Int_t labDp=-1;
	fIsLc = 0;
	fIsLcResonant=0;
	if(fReadMC){ //MC, check if Lc prompt or non prompt
		AliVertexingHFUtils *util = new AliVertexingHFUtils();
		labDp = MatchToMCLambdac(part,arrayMC);
		if(labDp>0){
			fIsLc=1;
			AliAODMCParticle *partDp = (AliAODMCParticle*)arrayMC->At(labDp);
			Int_t pdgMom=util->CheckOrigin(arrayMC,partDp,kFALSE);
			if(pdgMom == 5) fIsLc=2;
		}
		delete util;


		if(fIsLc){ //sig, MC

			AliCFVertexingHF3Prong *Lcres0 = new AliCFVertexingHF3Prong(arrayMC,0,32,1); //nonresonant
			AliCFVertexingHF3Prong *Lcres1 = new AliCFVertexingHF3Prong(arrayMC,0,32,2); // Lc -> L(1520) + p
			AliCFVertexingHF3Prong *Lcres2 = new AliCFVertexingHF3Prong(arrayMC,0,32,3); //Lc -> K* + pi
			AliCFVertexingHF3Prong *Lcres3 = new AliCFVertexingHF3Prong(arrayMC,0,32,4); // Lc -> Delta + K

			Lcres0->SetRecoCandidateParam(part);
			Lcres1->SetRecoCandidateParam(part);
			Lcres2->SetRecoCandidateParam(part);
			Lcres3->SetRecoCandidateParam(part);

			Int_t mclabel0 = Lcres0->GetMCLabel();
			Int_t mclabel1 = Lcres1->GetMCLabel();
			Int_t mclabel2 = Lcres2->GetMCLabel();
			Int_t mclabel3 = Lcres3->GetMCLabel();

			Lcres0->SetMCCandidateParam(mclabel0);
			Lcres1->SetMCCandidateParam(mclabel1);
			Lcres2->SetMCCandidateParam(mclabel2);
			Lcres3->SetMCCandidateParam(mclabel3);

			if(Lcres0->CheckLc3Prong()) fIsLcResonant=1;
			else if(Lcres1->CheckLc3Prong()) fIsLcResonant=2;
			else if(Lcres2->CheckLc3Prong()) fIsLcResonant=3;
			else if(Lcres3->CheckLc3Prong()) fIsLcResonant=4;
			else fIsLcResonant=-1;

			delete Lcres0;
			delete Lcres1;
			delete Lcres2;
			delete Lcres3;
		}
	}
}

//---------------------------

void AliAnalysisTaskSELambdacTMVA::FillNtuple(AliAODEvent *aod,AliAODRecoDecayHF3Prong *part,
																							TClonesArray *arrayMC, Int_t selection)
{
	//
	// Function to fill NTuple with candidate's variables
	//

	Bool_t IsInjected   = -1;
	Bool_t IsLc		= 0;
	Bool_t IsLcfromLb	= 0;

	if(fReadMC){ 
		AliAODMCHeader *mcHeader2 = (AliAODMCHeader*)aod->GetList()->FindObject(AliAODMCHeader::StdBranchName());
		AliVertexingHFUtils *util = new AliVertexingHFUtils();
		IsInjected = util->IsCandidateInjected(part,mcHeader2,arrayMC)?1:0;
		delete util;
	}
	if(fIsLc>=1 && fIsLc<=2) IsLc=kTRUE;
	if(fIsLc==2) IsLcfromLb=kTRUE;

	Double_t invMasspKpi=-1.;
	Double_t invMasspiKp=-1.;
	//apply MC PID
	if(fReadMC && fMCPid){
		if(IspKpiMC(part,arrayMC)) invMasspKpi=part->InvMassLcpKpi();
		else if(IspiKpMC(part,arrayMC)) invMasspiKp=part->InvMassLcpiKp();
		//else return;
	}
	// apply realistic PID
	if(fRealPid){
		if(selection==1 || selection==3) invMasspKpi=part->InvMassLcpKpi();
		else if(selection>=2) invMasspiKp=part->InvMassLcpiKp();
	}
	//    Float_t centrality=fRDCutsAnalysis->GetCentrality(aod);
	//    Int_t runNumber=aod->GetRunNumber();

	//fill ntuple
	// 1 - loose pid
	// 2 - bayesian pid
	Float_t tmp[41];
	//Is Lc
	if(!IsInjected && IsLc==0) tmp[0]=0; //non-injected bkg
	else if(IsLc==1 && !IsLcfromLb) tmp[0]=1; //prompt Lc
	else if(IsLc==1 && IsLcfromLb) tmp[0]=2; //non-prompt Lc
	else if(IsInjected && IsLc==0) tmp[0]=3; //injected bkg
	else tmp[0]=-99; //should not happen

	//invariant mass
	tmp[2]=part->InvMassLcpiKp();
	tmp[1]=part->InvMassLcpKpi(); 
	//pt of decay products
	tmp[3]=part->PtProng(0);
	tmp[4]=part->PtProng(1); //pt kaon
	tmp[5]=part->PtProng(2);
	tmp[6]=part->Pt();
	tmp[7]=part->CosPointingAngle();
	tmp[8]=part->DecayLength();
	tmp[9]=part->NormalizedDecayLength();
	tmp[10]=TMath::Min(part->GetDist12toPrim(),part->GetDist23toPrim());
	tmp[11]=part->GetSigmaVert();
	Double_t dcas[3]={0};
	part->GetDCAs(dcas);
	tmp[12]=TMath::Max(dcas[0],TMath::Max(dcas[1],dcas[2]));
	tmp[13]=part->DecayLengthXY();
	tmp[14]=part->NormalizedDecayLengthXY();
	/*tmp[18]=part->Getd0Prong(0);
		tmp[19]=part->Getd0Prong(1);
		tmp[20]=part->Getd0Prong(2);
		tmp[21]=part->Getd0Prong(0)*part->Getd0Prong(0)+part->Getd0Prong(1)*part->Getd0Prong(1)+part->Getd0Prong(2)*part->Getd0Prong(2);
		tmp[22]=(Float_t)centrality;
		tmp[23]=(Float_t)runNumber;*/

	//Check resonant decays for efficiencies
	tmp[15]=fIsLcResonant; // bkg
	//PID selection
	tmp[16]=selection;

	AliVTrack *track0=dynamic_cast<AliVTrack*>(part->GetDaughter(0));
	AliVTrack *track1=dynamic_cast<AliVTrack*>(part->GetDaughter(1));
	AliVTrack *track2=dynamic_cast<AliVTrack*>(part->GetDaughter(2));
	if(fFillNtuple>=1) {
		if(fRDCutsAnalysis->GetPidHF()->GetUseCombined()) {//check this
			//bayesian probabilities
			Double_t prob0[AliPID::kSPECIES];
			Double_t prob1[AliPID::kSPECIES];
			Double_t prob2[AliPID::kSPECIES];

			if (!track0 || !track1 || !track2) {
				AliError("AliVTrack missing - wont fill Ntuple");
				return;
			}
			fRDCutsAnalysis->GetPidHF()->GetPidCombined()->ComputeProbabilities(track0,fRDCutsAnalysis->GetPidHF()->GetPidResponse(),prob0);
			fRDCutsAnalysis->GetPidHF()->GetPidCombined()->ComputeProbabilities(track1,fRDCutsAnalysis->GetPidHF()->GetPidResponse(),prob1);
			fRDCutsAnalysis->GetPidHF()->GetPidCombined()->ComputeProbabilities(track2,fRDCutsAnalysis->GetPidHF()->GetPidResponse(),prob2);
			//if(prob0[AliPID::kPion] < 0.3 && prob0[AliPID::kProton] < 0.3) return;
			//if(prob1[AliPID::kKaon] < 0.3) return;
			//if(prob2[AliPID::kPion] < 0.3 && prob2[AliPID::kProton] < 0.3) return;
			tmp[17]=prob0[AliPID::kPion];			//track 0, pion
			tmp[18]=prob0[AliPID::kKaon];     		//kaon
			tmp[19]=prob0[AliPID::kProton];			//proton
			tmp[20]=prob1[AliPID::kPion];			//track 1, pion		
			tmp[21]=prob1[AliPID::kKaon];     		//kaon
			tmp[22]=prob1[AliPID::kProton];			//proton
			tmp[23]=prob2[AliPID::kPion];			//track 2, pion
			tmp[24]=prob2[AliPID::kKaon];     		//kaon
			tmp[25]=prob2[AliPID::kProton];			//proton
		}
		else {
			//fill w 0
			for(Int_t iprob=17;iprob<=25;iprob++) {
				tmp[iprob]=-1;
			}
		}
	}
	if(fFillNtuple>=2) {
		tmp[26]=fPIDResponse->NumberOfSigmasTPC(track0,AliPID::kPion);			//track 0, pion
		tmp[27]=fPIDResponse->NumberOfSigmasTPC(track0,AliPID::kKaon);     		//kaon
		tmp[28]=fPIDResponse->NumberOfSigmasTPC(track0,AliPID::kProton);			//proton
		tmp[29]=fPIDResponse->NumberOfSigmasTPC(track1,AliPID::kPion);			//track 1, pion		
		tmp[30]=fPIDResponse->NumberOfSigmasTPC(track1,AliPID::kKaon);     		//kaon
		tmp[31]=fPIDResponse->NumberOfSigmasTPC(track1,AliPID::kProton);			//proton
		tmp[32]=fPIDResponse->NumberOfSigmasTPC(track2,AliPID::kPion);			//track 2, pion
		tmp[33]=fPIDResponse->NumberOfSigmasTPC(track2,AliPID::kKaon);     		//kaon
		tmp[34]=fPIDResponse->NumberOfSigmasTPC(track2,AliPID::kProton);	
	}
	if(fFillNtuple>=3) {
		tmp[35]=part->InvMass2Prongs(1,0,321,2212); //inv mass pK
		tmp[36]=part->InvMass2Prongs(2,1,211,321); //inv mass Kpi
		tmp[37]=part->InvMass2Prongs(2,0,211,2212);//inv mass ppi
		tmp[38]=part->InvMass2Prongs(1,2,321,2212); //inv mass Kp 
		tmp[39]=part->InvMass2Prongs(0,1,211,321); //inv mass piK
		tmp[40]=part->InvMass2Prongs(0,2,211,2212);//inv mass pip
	}

	fNtupleLambdac->Fill(tmp);
	PostData(5,fNtupleLambdac);

	return;
}

//---------------------------

void AliAnalysisTaskSELambdacTMVA::FillEffHists(Int_t kStep) {
	
	//
	// Fill histograms (pt, pt vs eta, pt vs y, pt vs phi with 
	// candidates passing each step of the analysis
	//

	//fill according to kStep - hArray[nSteps]
	if(!fIsLc) {	
		fhNBkgNI[kStep]->Fill(fCandidateVars[0]);
		fhPtEtaBkgNI[kStep]->Fill(fCandidateVars[0],fCandidateVars[1]);
		fhPtYBkgNI[kStep]->Fill(fCandidateVars[0],fCandidateVars[2]);
		fhPtPhiBkgNI[kStep]->Fill(fCandidateVars[0],fCandidateVars[3]);
	}

	else{
		fhNLc[kStep]->Fill(fCandidateVars[0]);
		fhPtEtaLc[kStep]->Fill(fCandidateVars[0],fCandidateVars[1]);
		fhPtYLc[kStep]->Fill(fCandidateVars[0],fCandidateVars[2]);
		fhPtPhiLc[kStep]->Fill(fCandidateVars[0],fCandidateVars[3]);

		if(fIsLc==1) { //prompt Lc
			fhNLcc[kStep]->Fill(fCandidateVars[0]);
			fhPtEtaLcc[kStep]->Fill(fCandidateVars[0],fCandidateVars[1]);
			fhPtYLcc[kStep]->Fill(fCandidateVars[0],fCandidateVars[2]);
			fhPtPhiLcc[kStep]->Fill(fCandidateVars[0],fCandidateVars[3]);
			if(fIsLcResonant==1){
				fhNLcNonRc[kStep]->Fill(fCandidateVars[0]);
				fhPtEtaLcNonRc[kStep]->Fill(fCandidateVars[0],fCandidateVars[1]);
				fhPtYLcNonRc[kStep]->Fill(fCandidateVars[0],fCandidateVars[2]);
				fhPtPhiLcNonRc[kStep]->Fill(fCandidateVars[0],fCandidateVars[3]);
			}
			else if(fIsLcResonant==2){
				fhNLcL1520c[kStep]->Fill(fCandidateVars[0]);
				fhPtEtaLcL1520c[kStep]->Fill(fCandidateVars[0],fCandidateVars[1]);
				fhPtYLcL1520c[kStep]->Fill(fCandidateVars[0],fCandidateVars[2]);
				fhPtPhiLcL1520c[kStep]->Fill(fCandidateVars[0],fCandidateVars[3]);
			}
			else if(fIsLcResonant==3){
				fhNLcKstarc[kStep]->Fill(fCandidateVars[0]);
				fhPtEtaLcKstarc[kStep]->Fill(fCandidateVars[0],fCandidateVars[1]);
				fhPtYLcKstarc[kStep]->Fill(fCandidateVars[0],fCandidateVars[2]);
				fhPtPhiLcKstarc[kStep]->Fill(fCandidateVars[0],fCandidateVars[3]);
			}
			else if(fIsLcResonant==4){ 
				fhNLcDeltac[kStep]->Fill(fCandidateVars[0]);
				fhPtEtaLcDeltac[kStep]->Fill(fCandidateVars[0],fCandidateVars[1]);
				fhPtYLcDeltac[kStep]->Fill(fCandidateVars[0],fCandidateVars[2]);
				fhPtPhiLcDeltac[kStep]->Fill(fCandidateVars[0],fCandidateVars[3]);
			}
		}
		else if(fIsLc==2) { //non-prompt Lc
			fhNLcb[kStep]->Fill(fCandidateVars[0]);
			fhPtEtaLcb[kStep]->Fill(fCandidateVars[0],fCandidateVars[1]);
			fhPtYLcb[kStep]->Fill(fCandidateVars[0],fCandidateVars[2]);
			fhPtPhiLcb[kStep]->Fill(fCandidateVars[0],fCandidateVars[3]);
			if(fIsLcResonant==1){
				fhNLcNonRb[kStep]->Fill(fCandidateVars[0]);
				fhPtEtaLcNonRb[kStep]->Fill(fCandidateVars[0],fCandidateVars[1]);
				fhPtYLcNonRb[kStep]->Fill(fCandidateVars[0],fCandidateVars[2]);
				fhPtPhiLcNonRb[kStep]->Fill(fCandidateVars[0],fCandidateVars[3]);
			}
			else if(fIsLcResonant==2){
				fhNLcL1520b[kStep]->Fill(fCandidateVars[0]);
				fhPtEtaLcL1520b[kStep]->Fill(fCandidateVars[0],fCandidateVars[1]);
				fhPtYLcL1520b[kStep]->Fill(fCandidateVars[0],fCandidateVars[2]);
				fhPtPhiLcL1520b[kStep]->Fill(fCandidateVars[0],fCandidateVars[3]);
			}
			else if(fIsLcResonant==3){
				fhNLcKstarb[kStep]->Fill(fCandidateVars[0]);
				fhPtEtaLcKstarb[kStep]->Fill(fCandidateVars[0],fCandidateVars[1]);
				fhPtYLcKstarb[kStep]->Fill(fCandidateVars[0],fCandidateVars[2]);
				fhPtPhiLcKstarb[kStep]->Fill(fCandidateVars[0],fCandidateVars[3]);
			}
			else if(fIsLcResonant==4){ 
				fhNLcDeltab[kStep]->Fill(fCandidateVars[0]);
				fhPtEtaLcDeltab[kStep]->Fill(fCandidateVars[0],fCandidateVars[1]);
				fhPtYLcDeltab[kStep]->Fill(fCandidateVars[0],fCandidateVars[2]);
				fhPtPhiLcDeltab[kStep]->Fill(fCandidateVars[0],fCandidateVars[3]);
			}
		}
	}
}



//-----------------------------
Bool_t AliAnalysisTaskSELambdacTMVA::IspKpiMC(AliAODRecoDecayHF3Prong *d,TClonesArray *arrayMC) const{

	//
	// Apply MC PID
	// 

	Int_t lab[3]={0,0,0},pdgs[3]={0,0,0};
	for(Int_t i=0;i<3;i++){
		AliAODTrack *daugh=(AliAODTrack*)d->GetDaughter(i);
		lab[i]=daugh->GetLabel();
		if(lab[i]<0) return kFALSE;
		AliAODMCParticle *part= (AliAODMCParticle*)arrayMC->At(lab[i]);
		if(!part) return kFALSE;
		pdgs[i]=TMath::Abs(part->GetPdgCode());
	}

	if(pdgs[0]==2212 && pdgs[1]==321 && pdgs[2]==211) return kTRUE;

	return kFALSE;
}
//-----------------------------
Bool_t AliAnalysisTaskSELambdacTMVA::IspiKpMC(AliAODRecoDecayHF3Prong *d,TClonesArray *arrayMC) const{

	//
	// Apply MC PID
	//

	Int_t lab[3]={0,0,0},pdgs[3]={0,0,0};
	for(Int_t i=0;i<3;i++){
		AliAODTrack *daugh=(AliAODTrack*)d->GetDaughter(i);
		lab[i]=daugh->GetLabel();
		if(lab[i]<0) return kFALSE;
		AliAODMCParticle *part= (AliAODMCParticle*)arrayMC->At(lab[i]);
		if(!part) return kFALSE;
		pdgs[i]=TMath::Abs(part->GetPdgCode());
	}

	if(pdgs[2]==2212 && pdgs[1]==321 && pdgs[0]==211) {return kTRUE;}

	return kFALSE;
}

