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

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskSELambdacTMVA);
/// \endcond

	//________________________________________________________________________
	AliAnalysisTaskSELambdacTMVA::AliAnalysisTaskSELambdacTMVA():
		AliAnalysisTaskSE(),
		fOutput(0), 
		fHistNEvents(0),
		fHistNEventsRejTM(0),
		fhSelectBit(0),
		fhSelectionBits(0),
		fhSelectionBitsSigc(0),
		fhSelectionBitsSigb(0),
		fhSetIsLc(0),
		fhPIDmassLcPt(0),
		fhPIDmassLcPtSig(0),
		fhPIDmassLcPtSigc(0),
		fhPIDmassLcPtSigb(0),
		fhMCmassLcPt(0),
		fhMCmassLcPtSig(0),
		fhMCmassLcPtSigc(0),
		fhMCmassLcPtSigb(0),
		fhProbmassLcPt(0),
		fhProbmassLcPtSig(0),
		fhProbmassLcPtSigc(0),
		fhProbmassLcPtSigb(0),
		fhIsLcResonantGen(0),
		fhIsLcResonantReco(0),
		fhIsLcGen(0),
		fhIsLcReco(0),
		fhRecoPDGmom(0),
		fhPtMisIdpKpi(0),
		fhPtMisIdpiKp(0),
		fhPtCorrId(0),
		fhInvMassMisIdpKpi(0),
		fhInvMassMisIdpiKp(0),
		fhPtMisIdpKpiProb(0),
		fhPtMisIdpiKpProb(0),
		fhPtCorrIdProb(0),
		fhInvMassMisIdpKpiProb(0),
		fhInvMassMisIdpiKpProb(0),
		fNtupleLambdac(0),
		fNtupleLambdacReco(0),
		fFuncWeightPythia(0),
		fFuncWeightFONLL7overLHC10f6a(0),
		fFuncWeightFONLL5overLHC13d3(0),
		fFuncWeightFONLL5overLHC10f6a(0),
		fFuncWeightFONLL5overLHC13d3Lc(0),
		fFuncWeightFONLL7overLHC11b2Lc(0),
		fFuncWeightFONLL7overLHC10f7aLc(0),
		fUseNchWeight(kFALSE),
		fHistoMCNch(0x0),
		fIsLc(0),
		fIsLcResonant(0),
		fPtLc(0.),
		fUpmasslimit(2.486),
		fLowmasslimit(2.086),
		fRDCutsAnalysis(0),
		fListCuts(0),
		fFillNtuple(0),
		fFillNtupleReco(0),
		fKeepLcNotFromQuark(kFALSE),
		fKeepBkgNt(kTRUE),
		fSyst(2),
		fReadMC(kFALSE),
		fMCPid(kFALSE),
		fRealPid(kFALSE),
		fResPid(kTRUE),
		fUseKF(kFALSE),
		fAnalysis(kFALSE),
		fVHF(0),
		fLcCut(kFALSE),
		fLcPIDCut(kFALSE),    
		fIsHijing(kFALSE),
		fNentries(0),
		fPIDResponse(0),
		fCounter(0),
		fVertUtil(0)

{
	//
	/// Default constructor
	//

	for(Int_t i=0;i<12;i++) { 
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
AliAnalysisTaskSELambdacTMVA::AliAnalysisTaskSELambdacTMVA(const char *name,Int_t fillNtuple,Int_t fillNtupleReco,AliRDHFCutsLctopKpi *lccutsana):
	AliAnalysisTaskSE(name),
	fOutput(0), 
	fHistNEvents(0),
	fHistNEventsRejTM(0),
	fhSelectBit(0),
	fhSelectionBits(0),
	fhSelectionBitsSigc(0),
	fhSelectionBitsSigb(0),
	fhSetIsLc(0),
	fhPIDmassLcPt(0),
	fhPIDmassLcPtSig(0),
	fhPIDmassLcPtSigc(0),
	fhPIDmassLcPtSigb(0),
	fhMCmassLcPt(0),
	fhMCmassLcPtSig(0),
	fhMCmassLcPtSigc(0),
	fhMCmassLcPtSigb(0),
	fhProbmassLcPt(0),
	fhProbmassLcPtSig(0),
	fhProbmassLcPtSigc(0),
	fhProbmassLcPtSigb(0),
	fhIsLcResonantGen(0),
	fhIsLcResonantReco(0),
	fhIsLcGen(0),
	fhIsLcReco(0),
	fhRecoPDGmom(0),
	fhPtMisIdpKpi(0),
	fhPtMisIdpiKp(0),
	fhPtCorrId(0),
	fhInvMassMisIdpKpi(0),
	fhInvMassMisIdpiKp(0),
	fhPtMisIdpKpiProb(0),
	fhPtMisIdpiKpProb(0),
	fhPtCorrIdProb(0),
	fhInvMassMisIdpKpiProb(0),
	fhInvMassMisIdpiKpProb(0),
	fNtupleLambdac(0),
	fNtupleLambdacReco(0),
	fFuncWeightPythia(0),
	fFuncWeightFONLL7overLHC10f6a(0),
	fFuncWeightFONLL5overLHC13d3(0),
	fFuncWeightFONLL5overLHC10f6a(0),
	fFuncWeightFONLL5overLHC13d3Lc(0),
	fFuncWeightFONLL7overLHC11b2Lc(0),
	fFuncWeightFONLL7overLHC10f7aLc(0),
	fUseNchWeight(kFALSE),
	fHistoMCNch(0x0),
	fIsLc(0),
	fIsLcResonant(0),
	fPtLc(0.),
	fUpmasslimit(2.486),
	fLowmasslimit(2.086),
	fRDCutsAnalysis(lccutsana),
	fListCuts(0),
	fFillNtuple(fillNtuple),
	fFillNtupleReco(fillNtupleReco),
	fKeepLcNotFromQuark(kFALSE),
	fKeepBkgNt(kTRUE),
	fSyst(2),
	fReadMC(kFALSE),
	fMCPid(kFALSE),
	fRealPid(kFALSE),
	fResPid(kTRUE),
	fUseKF(kFALSE),
	fAnalysis(kFALSE),
	fVHF(0),
	fLcCut(kFALSE),
	fLcPIDCut(kFALSE),    
	fIsHijing(kFALSE),
	fNentries(0),
	fPIDResponse(0),
	fCounter(0),
	fVertUtil(0)
{
	//
	/// Default constructor
	/// Output slot #1 writes into a TList container
	//
	for(Int_t i=0;i<12;i++) { 
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
	fVertUtil = new AliVertexingHFUtils();

	DefineOutput(1,TList::Class());  //My private output
	DefineOutput(2,TList::Class());
	DefineOutput(3,TH1F::Class());
	DefineOutput(4,AliNormalizationCounter::Class());
	if (fFillNtuple) {
		// Output slot #2 writes into a TNtuple container
		DefineOutput(5,TNtuple::Class());  //My private output
	}
	if(fFillNtupleReco) DefineOutput(6,TNtuple::Class());  //My private output
}


//________________________________________________________________________
AliAnalysisTaskSELambdacTMVA::~AliAnalysisTaskSELambdacTMVA()
{
	//
	/// Destructor
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
	if(fVertUtil) {
		delete fVertUtil;
		fVertUtil = 0;
	}

}  



//_________________________________________________________________
void AliAnalysisTaskSELambdacTMVA::Init()
{
	//
	/// Initialization
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
	/// Create the output container
	//
	
	if (fDebug > 1) printf("AnalysisTaskSELambdac::UserCreateOutputObjects() \n");

	// Several histograms are more conveniently managed in a TList
	fOutput = new TList();
	fOutput->SetOwner();
	fOutput->SetName("OutputHistos");

	TString hisname,histitle;

	//Lc bit QA
	fhSelectBit			= new TH1F("hSelectBit","hSelectBit",5,-0.5,4.5);
	fhSelectBit->GetXaxis()->SetBinLabel(2,"All");
	fhSelectBit->GetXaxis()->SetBinLabel(3,"SelectionMap");
	fhSelectBit->GetXaxis()->SetBinLabel(4,"!LcCut");
	fhSelectBit->GetXaxis()->SetBinLabel(5,"!LcPID");
	fhSelectBit->GetXaxis()->SetNdivisions(1,kFALSE);
	fOutput->Add(fhSelectBit);

	fHistNEvents = new TH1F("fHistNEvents", "Number of processed events; ; Events",3,-0.5,2.5);
	fHistNEvents->GetXaxis()->SetBinLabel(2,"N events");
	fHistNEvents->GetXaxis()->SetBinLabel(3,"N events (after selection)");
	fHistNEvents->GetXaxis()->SetNdivisions(1,kFALSE);
	fHistNEvents->Sumw2();
	fOutput->Add(fHistNEvents);

	fhPIDmassLcPt = new TH2F("hPIDmassLcPt","hPIDmassLcPt;3-Prong p_{T} GeV/c;Invariant mass pK#pi (GeV/c)",150,0.,15.,(fUpmasslimit-fLowmasslimit)/0.002,fLowmasslimit,fUpmasslimit);
	fhPIDmassLcPtSig = new TH2F("hPIDmassLcPtSig","hPIDmassLcPtSig;3-Prong signal p_{T} GeV/c;Invariant mass pK#pi (GeV/c)",150,0.,15.,(fUpmasslimit-fLowmasslimit)/0.002,fLowmasslimit,fUpmasslimit);
	fhPIDmassLcPtSigc = new TH2F("hPIDmassLcPtSigc","hPIDmassLcPtSigc;3-Prong signal p_{T} GeV/c;Invariant mass pK#pi (GeV/c)",150,0.,15.,(fUpmasslimit-fLowmasslimit)/0.002,fLowmasslimit,fUpmasslimit);
	fhPIDmassLcPtSigb = new TH2F("hPIDmassLcPtSigb","hPIDmassLcPtSigb;3-Prong signal p_{T} GeV/c;Invariant mass pK#pi (GeV/c)",150,0.,15.,(fUpmasslimit-fLowmasslimit)/0.002,fLowmasslimit,fUpmasslimit);
	fhMCmassLcPt = new TH2F("hMCmassLcPt","hMCmassLcPt;3-Prong p_{T} GeV/c;Invariant mass pK#pi (GeV/c)",150,0.,15.,(fUpmasslimit-fLowmasslimit)/0.002,fLowmasslimit,fUpmasslimit);
	fhMCmassLcPtSig = new TH2F("hMCmassLcPtSig","hMCmassLcPtSig;3-Prong signal p_{T} GeV/c;Invariant mass pK#pi (GeV/c)",150,0.,15.,(fUpmasslimit-fLowmasslimit)/0.002,fLowmasslimit,fUpmasslimit);
	fhMCmassLcPtSigc = new TH2F("hMCmassLcPtSigc","hMCmassLcPtSigc;3-Prong signal p_{T} GeV/c;Invariant mass pK#pi (GeV/c)",150,0.,15.,(fUpmasslimit-fLowmasslimit)/0.002,fLowmasslimit,fUpmasslimit);
	fhMCmassLcPtSigb = new TH2F("hMCmassLcPtSigb","hMCmassLcPtSigb;3-Prong signal p_{T} GeV/c;Invariant mass pK#pi (GeV/c)",150,0.,15.,(fUpmasslimit-fLowmasslimit)/0.002,fLowmasslimit,fUpmasslimit);
	fhProbmassLcPt = new TH2F("hProbmassLcPt","hProbmassLcPt;3-Prong p_{T} GeV/c;Invariant mass pK#pi (GeV/c)",150,0.,15.,(fUpmasslimit-fLowmasslimit)/0.002,fLowmasslimit,fUpmasslimit);
	fhProbmassLcPtSig = new TH2F("hProbmassLcPtSig","hProbmassLcPtSig;3-Prong signal p_{T} GeV/c;Invariant mass pK#pi (GeV/c)",150,0.,15.,(fUpmasslimit-fLowmasslimit)/0.002,fLowmasslimit,fUpmasslimit);
	fhProbmassLcPtSigc = new TH2F("hProbmassLcPtSigc","hProbmassLcPtSigc;3-Prong signal p_{T} GeV/c;Invariant mass pK#pi (GeV/c)",150,0.,15.,(fUpmasslimit-fLowmasslimit)/0.002,fLowmasslimit,fUpmasslimit);
	fhProbmassLcPtSigb = new TH2F("hProbmassLcPtSigb","hProbmassLcPtSigb;3-Prong signal p_{T} GeV/c;Invariant mass pK#pi (GeV/c)",150,0.,15.,(fUpmasslimit-fLowmasslimit)/0.002,fLowmasslimit,fUpmasslimit);
	fOutput->Add(fhPIDmassLcPt);
	fOutput->Add(fhPIDmassLcPtSig);
	fOutput->Add(fhPIDmassLcPtSigc);
	fOutput->Add(fhPIDmassLcPtSigb);
	fOutput->Add(fhMCmassLcPt);
	fOutput->Add(fhMCmassLcPtSig);
	fOutput->Add(fhMCmassLcPtSigc);
	fOutput->Add(fhMCmassLcPtSigb);
	fOutput->Add(fhProbmassLcPt);
	fOutput->Add(fhProbmassLcPtSig);
	fOutput->Add(fhProbmassLcPtSigc);
	fOutput->Add(fhProbmassLcPtSigb);

	fhIsLcResonantGen = new TH1F("hIsLcResonantGen","IsLcResonant flag gen",6,-1.5,4.5);
	fhIsLcResonantReco = new TH1F("hIsLcResonantReco","IsLcResonant flag reco",6,-1.5,4.5);
	fhIsLcGen = new TH1F("hIsLcGen","IsLc flag gen",4,-1.5,2.5);
	fhIsLcReco = new TH1F("hIsLcReco","IsLc flag reco",4,-1.5,2.5);
	fOutput->Add(fhIsLcResonantGen);
	fOutput->Add(fhIsLcResonantReco);
	fOutput->Add(fhIsLcGen);
	fOutput->Add(fhIsLcReco);

	fhRecoPDGmom = new TH1F("hRecoPDGmom","pdg of mother reco. MatchToMCLambdac",7,-0.5,6.5);
	fOutput->Add(fhRecoPDGmom);

	fhSetIsLc = new TH1F("hSetIsLc","Check candidates before/after rec. set is Lc",2,-0.5,1.5);
	fOutput->Add(fhSetIsLc);

	fhSelectionBits = new TH2F("hSelectionBits","Reconstruction + selection bit",13,-0.5,12.5,150,0,15);
	fhSelectionBits->GetXaxis()->SetBinLabel(2,"D0toKpiCuts");
	fhSelectionBits->GetXaxis()->SetBinLabel(3,"D0toKpiPID");
	fhSelectionBits->GetXaxis()->SetBinLabel(4,"D0fromDstarCuts");
	fhSelectionBits->GetXaxis()->SetBinLabel(5,"D0fromDstarPID");
	fhSelectionBits->GetXaxis()->SetBinLabel(6,"DplusCuts");
	fhSelectionBits->GetXaxis()->SetBinLabel(7,"DplusPID");
	fhSelectionBits->GetXaxis()->SetBinLabel(8,"DsCuts");
	fhSelectionBits->GetXaxis()->SetBinLabel(9,"DsPID");
	fhSelectionBits->GetXaxis()->SetBinLabel(10,"LcCuts");
	fhSelectionBits->GetXaxis()->SetBinLabel(11,"LcPID");
	fhSelectionBits->GetXaxis()->SetBinLabel(12,"DstarCuts");
	fhSelectionBits->GetXaxis()->SetBinLabel(13,"DstarPID");
	fhSelectionBits->GetYaxis()->SetTitle("p_{T} (GeV/c)");
	fhSelectionBits->GetXaxis()->SetNdivisions(1,kFALSE);
	fOutput->Add(fhSelectionBits);

	fhSelectionBitsSigc = new TH2F("hSelectionBitsSigc","Reconstruction + selection bit from c",13,-0.5,12.5,150,0,15);
	fhSelectionBitsSigc->GetXaxis()->SetBinLabel(2,"D0toKpiCuts");
	fhSelectionBitsSigc->GetXaxis()->SetBinLabel(3,"D0toKpiPID");
	fhSelectionBitsSigc->GetXaxis()->SetBinLabel(4,"D0fromDstarCuts");
	fhSelectionBitsSigc->GetXaxis()->SetBinLabel(5,"D0fromDstarPID");
	fhSelectionBitsSigc->GetXaxis()->SetBinLabel(6,"DplusCuts");
	fhSelectionBitsSigc->GetXaxis()->SetBinLabel(7,"DplusPID");
	fhSelectionBitsSigc->GetXaxis()->SetBinLabel(8,"DsCuts");
	fhSelectionBitsSigc->GetXaxis()->SetBinLabel(9,"DsPID");
	fhSelectionBitsSigc->GetXaxis()->SetBinLabel(10,"LcCuts");
	fhSelectionBitsSigc->GetXaxis()->SetBinLabel(11,"LcPID");
	fhSelectionBitsSigc->GetXaxis()->SetBinLabel(12,"DstarCuts");
	fhSelectionBitsSigc->GetXaxis()->SetBinLabel(13,"DstarPID");
	fhSelectionBitsSigc->GetYaxis()->SetTitle("p_{T} (GeV/c)");
	fhSelectionBitsSigc->GetXaxis()->SetNdivisions(1,kFALSE);
	fOutput->Add(fhSelectionBitsSigc);

	fhSelectionBitsSigb = new TH2F("hSelectionBitsSigb","Reconstruction + selection bit from b",13,-0.5,13.5,150,0,15);
	fhSelectionBitsSigb->GetXaxis()->SetBinLabel(2,"D0toKpiCuts");
	fhSelectionBitsSigb->GetXaxis()->SetBinLabel(3,"D0toKpiPID");
	fhSelectionBitsSigb->GetXaxis()->SetBinLabel(4,"D0fromDstarCuts");
	fhSelectionBitsSigb->GetXaxis()->SetBinLabel(5,"D0fromDstarPID");
	fhSelectionBitsSigb->GetXaxis()->SetBinLabel(6,"DplusCuts");
	fhSelectionBitsSigb->GetXaxis()->SetBinLabel(7,"DplusPID");
	fhSelectionBitsSigb->GetXaxis()->SetBinLabel(8,"DsCuts");
	fhSelectionBitsSigb->GetXaxis()->SetBinLabel(9,"DsPID");
	fhSelectionBitsSigb->GetXaxis()->SetBinLabel(10,"LcCuts");
	fhSelectionBitsSigb->GetXaxis()->SetBinLabel(11,"LcPID");
	fhSelectionBitsSigb->GetXaxis()->SetBinLabel(12,"DstarCuts");
	fhSelectionBitsSigb->GetXaxis()->SetBinLabel(13,"DstarPID");
	fhSelectionBitsSigb->GetYaxis()->SetTitle("p_{T} (GeV/c)");
	fhSelectionBitsSigb->GetXaxis()->SetNdivisions(1,kFALSE);
	fOutput->Add(fhSelectionBitsSigb);

// enum ESele {kD0toKpiCuts,kD0toKpiPID,kD0fromDstarCuts,kD0fromDstarPID,kDplusCuts,kDplusPID,kDsCuts,kDsPID,kLcCuts,kLcPID,kDstarCuts,kDstarPID};

	TString stepnames[12] = {"GeneratedLimAcc","GeneratedAll","Generated","GeneratedAcc","Reco3Prong","LcBit","IsSelectedTracks","IsInFidAcc","PtRange","IsSelectedCandidate","IsSelectedPID","IsSelectedNtuple"};
	for(Int_t i=0;i<12;i++) { // histograms for efficiency cross check
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


	fhPtMisIdpKpi = new TH1F("hMisIdpKpi","pKpi id'd as piKp",30,0,30);
	fhPtMisIdpiKp = new TH1F("hMisIdpiKp","piKp id'd as pKpi",30,0,30);
	fhPtCorrId = new TH1F("hCorrId","Correctly id'd pKpi",30,0,30);
	fOutput->Add(fhPtMisIdpKpi);
	fOutput->Add(fhPtMisIdpiKp);
	fOutput->Add(fhPtCorrId);
	fhInvMassMisIdpKpi = new TH2F("hInvMassIdpKpi","inv mass pKpi id'd as piKp",300,0.,30.,(fUpmasslimit-fLowmasslimit)/0.002,fLowmasslimit,fUpmasslimit);
	fhInvMassMisIdpiKp = new TH2F("hInvMassIdpiKp","inv mass piKp id'd as pKpi",300,0.,30.,(fUpmasslimit-fLowmasslimit)/0.002,fLowmasslimit,fUpmasslimit);
	fOutput->Add(fhInvMassMisIdpKpi);
	fOutput->Add(fhInvMassMisIdpiKp);

	fhPtMisIdpKpiProb = new TH1F("hMisIdpKpiProb","pKpi id'd as piKp, most prob. PID",30,0,30);
	fhPtMisIdpiKpProb = new TH1F("hMisIdpiKpProb","piKp id'd as pKpi, most prob. PID",30,0,30);
	fhPtCorrIdProb = new TH1F("hCorrIdProb","Correctly id'd pKpi, most prob. PID",30,0,30);
	fOutput->Add(fhPtMisIdpKpiProb);
	fOutput->Add(fhPtMisIdpiKpProb);
	fOutput->Add(fhPtCorrIdProb);
	fhInvMassMisIdpKpiProb = new TH2F("hInvMassIdpKpiProb","inv mass pKpi id'd as piKp most prob. PID",300,0.,30.,(fUpmasslimit-fLowmasslimit)/0.002,fLowmasslimit,fUpmasslimit);
	fhInvMassMisIdpiKpProb = new TH2F("hInvMassIdpiKpProb","inv mass piKp id'd as pKpi most prob. PID",300,0.,30.,(fUpmasslimit-fLowmasslimit)/0.002,fLowmasslimit,fUpmasslimit);
	fOutput->Add(fhInvMassMisIdpKpiProb);
	fOutput->Add(fhInvMassMisIdpiKpProb);


  // weight function from ratio of flat value (1/30) to pythia 
	// use to normalise to flat distribution (should lead to flat pT distribution
	fFuncWeightPythia=new TF1("funcWeightPythia","1./(30. *[0]*x/TMath::Power(1.+(TMath::Power((x/[1]),[3])),[2]))",0.15,30);
	fFuncWeightPythia->SetParameters(0.36609,1.94635,1.40463,2.5);

  // weight function from the ratio of the LHC10f6a MC
  // and FONLL calculations for pp data
  fFuncWeightFONLL7overLHC10f6a=new TF1("funcWeightFONLL7overLHC10f6a","([0]*x)/TMath::Power([2],(1+TMath::Power([3],x/[1])))+[4]*TMath::Exp([5]+[6]*x)+[7]*TMath::Exp([8]*x)",0.15,40.);
  fFuncWeightFONLL7overLHC10f6a->SetParameters(2.41522e+01,4.92146e+00,6.72495e+00,2.5,6.15361e-03,4.78995e+00,-4.29135e-01,3.99421e-01,-1.57220e-02);

  // weight function from the ratio of the LHC13d3 MC
  // and FONLL calculations for pp data
  fFuncWeightFONLL5overLHC13d3=new TF1("funcWeightFONLL5overLHC13d3","([0]*x)/TMath::Power([2],(1+TMath::Power([3],x/[1])))+[4]*TMath::Exp([5]+[6]*x)+[7]*TMath::Exp([8]*x)",0.15,30.);
  fFuncWeightFONLL5overLHC13d3->SetParameters(2.94999e+00,3.47032e+00,2.81278e+00,2.5,1.93370e-02,3.86865e+00,-1.54113e-01,8.86944e-02,2.56267e-02);

  // weight function from the ratio of the LHC10f6a MC
  // and FONLL calculations for pp data
  fFuncWeightFONLL5overLHC10f6a=new TF1("funcWeightFONLL5overLHC10f6a","([0]*x)/TMath::Power([2],(1+TMath::Power([3],x/[1])))+[4]*TMath::Exp([5]+[6]*x)+[7]*TMath::Exp([8]*x)",0.15,40.);
  fFuncWeightFONLL5overLHC10f6a->SetParameters(2.77730e+01,4.78942e+00,7.45378e+00,2.5,9.86255e-02,2.30120e+00,-4.16435e-01,3.43770e-01,-2.29380e-02);

  // weight function from the ratio of the LHC13d3 MC
  // and FONLL calculations for pPb data, Lc
  fFuncWeightFONLL5overLHC13d3Lc=new TF1("funcWeightFONLL5overLHC13d3Lc","([0]*x)/TMath::Power([2],(1+TMath::Power([3],x/[1])))+[4]*TMath::Exp([5]+[6]*x)+[7]*TMath::Exp([8]*x)",0.15,20.);
	fFuncWeightFONLL5overLHC13d3Lc->SetParameters(5.94428e+01,1.63585e+01,9.65555e+00,6.71944e+00,8.88338e-02,2.40477e+00,-4.88649e-02,-6.78599e-01,-2.10951e-01);

  // weight function from the ratio of the LHC11b2 MC
  // and FONLL calculations for pp data, Lc
  fFuncWeightFONLL7overLHC11b2Lc=new TF1("funcWeightFONLL7overLHC11b2Lc","([0]*x)/TMath::Power([2],(1+TMath::Power([3],x/[1])))+[4]*TMath::Exp([5]+[6]*x)+[7]*TMath::Exp([8]*x)",1.,20.);
	fFuncWeightFONLL7overLHC11b2Lc->SetParameters(2.11879e+02,3.73290e+00,2.01235e+01,1.41508e+00,1.06268e-01,1.86285e+00,-4.52956e-02,-9.90631e-01,-1.31615e+00);

  // weight function from the ratio of the LHC107a MC
  // and FONLL calculations for pp data, Lc
  fFuncWeightFONLL7overLHC10f7aLc=new TF1("funcWeightFONLL7overLHC10f7aLc","([0]*x)/TMath::Power([2],(1+TMath::Power([3],x/[1])))+[4]*TMath::Exp([5]+[6]*x)+[7]*TMath::Exp([8]*x)",01.,20.);
	fFuncWeightFONLL7overLHC10f7aLc->SetParameters(2.84268e+02,2.18850e+01,2.36298e+01,7.46144e+00,1.69747e-01,1.66993e+00,-5.54726e-02,-1.53869e+00,-1.18404e+00);

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
		TString ntName="fNtupleLambdac";
		AliAnalysisDataContainer *contnt = GetOutputSlot(5)->GetContainer();
		if(contnt)ntName=(TString)contnt->GetName();
		if(fFillNtuple==1)       fNtupleLambdac = new TNtuple(ntName.Data(), "Lc", "isLcBkg:InvMasspKpi:InvMasspiKp:Charge:PtTr0:PtTr1:PtTr2:PtLc:CosP:DecayL:DecayLSig:Dist12:SigVert:DCA:DecayLXY:DecayLXYSig:isLcResonant:selectionPID:Tr0Ppi:Tr0PK:Tr0Pp:Tr1Ppi:Tr1PK:Tr1Pp:Tr2Ppi:Tr2PK:Tr2Pp");
		//more variables
		else if(fFillNtuple==2)  fNtupleLambdac = new TNtuple(ntName.Data(), "Lc", "isLcBkg:InvMasspKpi:InvMasspiKp:Charge:PtTr0:PtTr1:PtTr2:PtLc:CosP:DecayL:DecayLSig:Dist12:SigVert:DCATr0:DecayLXY:DecayLXYSig:isLcResonant:selectionPID:Tr0Ppi:Tr0PK:Tr0Pp:Tr1Ppi:Tr1PK:Tr1Pp:Tr2Ppi:Tr2PK:Tr2Pp:d00:d01:d02:d0Squared:d00Sig:d01Sig:d02Sig:d00SigResidual:d01SigResidual:d02SigResidual:CosPXY:DCATr1:DCATr2:Dist23:RunNumber");
	//weights
		else if(fFillNtuple==3)  fNtupleLambdac = new TNtuple(ntName.Data(), "Lc", "isLcBkg:InvMasspKpi:InvMasspiKp:Charge:PtTr0:PtTr1:PtTr2:PtLc:CosP:DecayL:DecayLSig:Dist12:SigVert:DCATr0:DecayLXY:DecayLXYSig:isLcResonant:selectionPID:Tr0Ppi:Tr0PK:Tr0Pp:Tr1Ppi:Tr1PK:Tr1Pp:Tr2Ppi:Tr2PK:Tr2Pp:d00:d01:d02:d0Squared:d00Sig:d01Sig:d02Sig:d00SigResidual:d01SigResidual:d02SigResidual:CosPXY:DCATr1:DCATr2:Dist23:RunNumber:PtLcMC:weightsPythia:weights7LHC106a:weights5LHC10f6a:weights5LHC13d3:weights5LHC13d3Lc:weights7LHC11b2Lc:weights7LHC10f7aLc:multiplicity:weightsMultiplicity");
		//2 prong decay products
		else if(fFillNtuple==4)  fNtupleLambdac = new TNtuple(ntName.Data(), "Lc", "isLcBkg:InvMasspKpi:InvMasspiKp:Charge:PtTr0:PtTr1:PtTr2:PtLc:CosP:DecayL:DecayLSig:Dist12:SigVert:DCATr0:DecayLXY:DecayLXYSig:isLcResonant:selectionPID:Tr0Ppi:Tr0PK:Tr0Pp:Tr1Ppi:Tr1PK:Tr1Pp:Tr2Ppi:Tr2PK:Tr2Pp:d00:d01:d02:d0Squared:d00Sig:d01Sig:d02Sig:d00SigResidual:d01SigResidual:d02SigResidual:CosPXY:DCATr1:DCATr2:Dist23:RunNumber:PtLcMC:weightsPythia:weights7LHC106a:weights5LHC10f6a:weights5LHC13d3:weights5LHC13d3Lc:weights7LHC11b2Lc:weights7LHC10f7aLc:multiplicity:weightsMultiplicity:InvMasspK:InvMassKpi:InvMassppi:InvMassKp:InvMasspiK:InvMasspip");
		else AliFatal("Invalid fill ntuple argument");
		PostData(5,fNtupleLambdac);

	}
	if(fFillNtupleReco) {
		//Reco ntuple
		TString ntName="fNtupleLambdacReco";
		AliAnalysisDataContainer *contntrec = GetOutputSlot(6)->GetContainer();
		if(contntrec)ntName=(TString)contntrec->GetName();
		fNtupleLambdacReco = new TNtuple(ntName.Data(), "Lc Reco", "isLcBkg:PtLc:PtLcMC:PtTr0:PtTr1:PtTr2:PtTr0MC:PtTr1MC:PtTr2MC:isTOFTr0:isTOFTr1:isTOFTr2:selectionCand:selectionPID:selectionPIDprob:Charge");
		PostData(6,fNtupleLambdacReco);
	}

	return;
}

//________________________________________________________________________
void AliAnalysisTaskSELambdacTMVA::UserExec(Option_t */*option*/)
{
	//
	/// Execute analysis for current event:
	/// heavy flavor candidates association to MC truth
	//
	AliAODEvent *aod = dynamic_cast<AliAODEvent*> (InputEvent());
	//tmp
	fHistNEvents->Fill(1); // count event
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
		//if(TString(AliAODMCHeader::StdBranchName()).Contains("LHC15f2")) {
		//	AliInfo("LHC15f2 - fIsHijing==kTRUE");
		//	fIsHijing=kTRUE;
		//}
		//else {
		//	AliInfo("Not LHC15f2 - fIsHijing==kFALSE");
		//	fIsHijing=kFALSE;
		//}
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
		//return;
	}
	fCounter->StoreEvent(aod,fRDCutsAnalysis,fReadMC);
	TString trigclass=aod->GetFiredTriggerClasses();
	//if(trigclass.Contains("C0SMH-B-NOPF-ALLNOTRD") || trigclass.Contains("C0SMH-B-NOPF-ALL")) 

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
	Bool_t isInFidAcc = kFALSE; 
	Bool_t isInAcc = kFALSE;
	if(fReadMC) {
		// Event selection as done is CF task
		Double_t zPrimVertex = vtx1 ->GetZ();
		Double_t zMCVertex = mcHeader->GetVtxZ();
		if (TMath::Abs(zMCVertex) > fRDCutsAnalysis->GetMaxVtxZ()){
			return;
		}
		for (Int_t iPart=0; iPart<arrayMC->GetEntriesFast(); iPart++) { 
			fIsLc=0;
			fIsLcResonant=0;
			AliAODMCParticle* mcPart = dynamic_cast<AliAODMCParticle*>(arrayMC->At(iPart));
			if (!mcPart){
				AliError("Failed casting particle from MC array!, Skipping particle");
				continue;
			}
			//Check whether particle is Lc
			SetIsLcGen(mcPart,arrayMC);
			if(fIsLc==0) continue;
			//-- is Lc --
			isInFidAcc = fRDCutsAnalysis->IsInFiducialAcceptance(mcPart->Pt(),mcPart->Y());
			fCandidateVars[0] = mcPart->Pt();
			fCandidateVars[1] = mcPart->Eta();
			fCandidateVars[2] = mcPart->Y();
			fCandidateVars[3] = mcPart->Phi();
			Int_t imother = -1;
			if(!fIsHijing) imother=mcPart->GetMother();
			if(imother>0) { //found mother
				AliAODMCParticle* mcPartMother = dynamic_cast<AliAODMCParticle*>(arrayMC->At(imother));
				if(!mcPart){
					AliError("Failed casting mother particle, Skipping particle");
					continue;
				}
			}

			//Check daughters
			SetLambdacDaugh(mcPart,arrayMC,isInAcc);
			if(fIsLcResonant>=1 && (fIsHijing || imother>0)) { //if signal, and is LHC15f2 or has mother
				AliDebug(2,"Lc has p K pi in final state");
				FillEffHists(kGeneratedAll);
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
			else if(fIsLcResonant==0){
				AliError("no p K pi final state");
				//					TClonesArray *ares = (TClonesArray*)aod->GetList()->FindObject(AliAODMCParticle::StdBranchName());
				//					Int_t nentriesres = ares->GetEntriesFast();
				//					for(Int_t ires=0;ires<nentriesres;ires++) {
				//						AliAODMCParticle* mcPartRes = dynamic_cast<AliAODMCParticle*>(ares->At(ires));
				//						Printf("%i",ires);
				//						mcPartRes->Print();
				//					}
				//					TIter next(aod->GetList());
				//					while (TObject *obj = next())
				//						obj->Print();
			}
			else AliError(Form("Not pKpi or background - should not happen! fIsLcResonant = %i",fIsLcResonant));
			fhIsLcResonantGen->Fill(Double_t(fIsLcResonant));
			fhIsLcGen->Fill(Double_t(fIsLc));
		}
	}


	fHistNEvents->Fill(2); // count event after event selection (as CF task)

	Bool_t isEvSelAnCuts=fRDCutsAnalysis->IsEventSelected(aod);
	if(!isEvSelAnCuts){ //reject if not selected with analysis cut
		if(fRDCutsAnalysis->GetWhyRejection()==1) // rejected for pileup
			fNentries->Fill(3);
		return;
	}

	//
	//Reconstruction level
	//loop over 3 prong candidates
	//

	for (Int_t i3Prong = 0; i3Prong < n3Prong; i3Prong++) {
		AliAODRecoDecayHF3Prong *d = (AliAODRecoDecayHF3Prong*)array3Prong->UncheckedAt(i3Prong);
		//(AliAODEvent *aod,AliAODRecoDecayHF3Prong *part,TClonesArray *arrayMC)
		fhSetIsLc->Fill(0);
		SetIsLcReco(d,arrayMC);
		fhSetIsLc->Fill(1);
		fhIsLcResonantReco->Fill(Double_t(fIsLcResonant));
		fhIsLcReco->Fill(Double_t(fIsLc));
		fCandidateVars[0] = d->Pt();
		fCandidateVars[1] = d->Eta();
		fCandidateVars[2] = d->Y(4122);
		fCandidateVars[3] = d->Phi();

		FillEffHists(kReco3Prong);
		fNentries->Fill(6);




		//add histogram with filter bit
		//one at start of task, one at filter bit selection
		//check every bit, bit vs pT
		//Filter bit selection and QA:
		fhSelectBit->Fill(1);
		if(d->GetSelectionMap()){
			// enum ESele {kD0toKpiCuts,kD0toKpiPID,kD0fromDstarCuts,kD0fromDstarPID,kDplusCuts,kDplusPID,kDsCuts,kDsPID,kLcCuts,kLcPID,kDstarCuts,kDstarPID};
			if(fIsLc==0) FillSelectionBits(d,fhSelectionBits); 
			else if(fIsLc==1) FillSelectionBits(d,fhSelectionBitsSigc);
			else if(fIsLc==2) FillSelectionBits(d,fhSelectionBitsSigb);

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

		//Fill reco ntuple
		if(fReadMC && fIsLc>=1 && fFillNtupleReco) FillRecoNtuple(aod,d,arrayMC);
		//idproton, pion using isSelectedPID
		//Bool_t isPID=fRDCutsAnalysis->GetIsUsePID();
		Int_t isSelectedPID=fRDCutsAnalysis->IsSelected(d,AliRDHFCuts::kPID,aod);

		if(isSelectedPID==0 ) fNentries->Fill(12);
		else if(isSelectedPID==1) fNentries->Fill(13);
		else if(isSelectedPID==2) fNentries->Fill(14);
		else fNentries->Fill(15);
		if(isSelectedPID>0) FillEffHists(kIsSelectedPID);
		//PID selection using maximum probability configuration
		Int_t selectionProb = 0;
    if(fRDCutsAnalysis->GetPidHF()) selectionProb = GetPIDselectionMaxProb(d);

		Int_t selection=fRDCutsAnalysis->IsSelected(d,AliRDHFCuts::kCandidate,aod);
		if(!selection) continue;
		FillEffHists(kIsSelectedCandidate);
		fNentries->Fill(11);

		FillMassHists(aod,d,arrayMC,selection,selectionProb);


		if(fFillNtuple) FillNtuple(aod,d,arrayMC,selection);
		FillEffHists(kIsSelectedNtuple);
		if(fIsLc>=1 && fIsLc <= 2) fNentries->Fill(16);
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
	/// Terminate analysis
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
	if(fFillNtupleReco){
		fNtupleLambdacReco = dynamic_cast<TNtuple*>(GetOutputData(6));
	}


	return;
}

//________________________________________________________________________
Int_t AliAnalysisTaskSELambdacTMVA::MatchToMCLambdac(AliAODRecoDecayHF3Prong *d,TClonesArray *arrayMC) const{

	//
	/// check if the candidate is a Lambdac decaying in pKpi or in the resonant channels
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
	/// return value of Lc resonant channel, from AOD MC particle
	/// Also check whether Lc daughters are in acceptance
	/// 0=not Lc, 1=non resonant Lc, 2=via L(1520) + pi, 3=via K* + p, 4=via Delta++ + K
	//

	Int_t numberOfLambdac=0;
	IsInAcc=kTRUE;
	if(TMath::Abs(part->GetPdgCode())!=4122) return 0;
	// Int_t daughTmp[2];
	// daughTmp[0]=part->GetDaughterLabel(0);
	// daughTmp[1]=part->GetDaughterLabel(1);
	Int_t nDaugh = (Int_t)part->GetNDaughters();
	if(nDaugh<2) return 0;
	if(nDaugh>3) return 0;
	AliAODMCParticle* pdaugh1 = (AliAODMCParticle*)arrayMC->At(part->GetDaughterLabel(0));
	if(!pdaugh1) {return 0;}
	Int_t number1 = TMath::Abs(pdaugh1->GetPdgCode());
	AliAODMCParticle* pdaugh2 = (AliAODMCParticle*)arrayMC->At(part->GetDaughterLabel(1));
	if(!pdaugh2) {return 0;}
	Int_t number2 = TMath::Abs(pdaugh2->GetPdgCode());

	AliDebug(2,"Is non resonant?");
	if(nDaugh==3){
		Int_t thirdDaugh=part->GetDaughterLabel(1)-1;
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
			AliAODMCParticle* pdaughK1 = (AliAODMCParticle*)arrayMC->At(pdaugh2->GetDaughterLabel(0));
			AliAODMCParticle* pdaughK2 = (AliAODMCParticle*)arrayMC->At(pdaugh2->GetDaughterLabel(1));
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
			AliAODMCParticle* pdaughK1 = (AliAODMCParticle*)arrayMC->At(pdaugh1->GetDaughterLabel(0));
			AliAODMCParticle* pdaughK2 = (AliAODMCParticle*)arrayMC->At(pdaugh1->GetDaughterLabel(1));
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
			AliAODMCParticle *pdaughD1=(AliAODMCParticle*)arrayMC->At(pdaugh2->GetDaughterLabel(0));
			AliAODMCParticle *pdaughD2=(AliAODMCParticle*)arrayMC->At(pdaugh2->GetDaughterLabel(1));
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
			AliAODMCParticle* pdaughD1 = (AliAODMCParticle*)arrayMC->At(pdaugh1->GetDaughterLabel(0));
			AliAODMCParticle* pdaughD2 = (AliAODMCParticle*)arrayMC->At(pdaugh1->GetDaughterLabel(1));
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
			AliAODMCParticle *pdaughL1=(AliAODMCParticle*)arrayMC->At(pdaugh1->GetDaughterLabel(0));
			AliAODMCParticle *pdaughL2=(AliAODMCParticle*)arrayMC->At(pdaugh1->GetDaughterLabel(1));
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
			AliAODMCParticle *pdaughL1=(AliAODMCParticle*)arrayMC->At(pdaugh2->GetDaughterLabel(0));
			AliAODMCParticle *pdaughL2=(AliAODMCParticle*)arrayMC->At(pdaugh2->GetDaughterLabel(1));
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

void AliAnalysisTaskSELambdacTMVA::SetIsLcGen(AliAODMCParticle *mcPart, TClonesArray *arrayMC) {

	//
	/// Set fIsLc from AliAODMCParticle
	/// fIsLc 0 = not Lc from quark, 1 = Lc from c, 2 = Lc from b
	//

	fIsLc=0;
	if(TMath::Abs(mcPart->GetPdgCode())==4122) {
		AliDebug(2,"Found Lc! now check mother");
		fIsLc=1;
		Int_t pdgMom = 0;
		pdgMom=fVertUtil->CheckOrigin(arrayMC,mcPart,fKeepLcNotFromQuark ? kFALSE : kTRUE);
		if(pdgMom == 5){
			AliDebug(2,"Lc comes from b");
			fIsLc=2;
		}
		else if(pdgMom==4) {
			fIsLc=1;
		}
		else {
			fIsLc=0;
			fhIsLcResonantGen->Fill(-1);
		}
	}
}


//-----------------------------

void AliAnalysisTaskSELambdacTMVA::SetIsLcReco(AliAODRecoDecayHF3Prong *part,
		TClonesArray *arrayMC) {

	//
	/// function that sets fIsLc and fIsLcResonant, from reconstructed 3 prong decay
	/// fIsLc - 0 = not Lc,  1 = Lc from c, 2 = Lc from b
	/// fIsLcResonant - 1= Non Resonant, 2=Lc->L1520+p, 3=Lc->K*+pi, 4=Lc->Delta+K
	//

	Int_t labDp=-1;
	fIsLc = 0;
	fIsLcResonant=0;
	if(!fReadMC) return;
	else{ //MC, check if Lc prompt or non prompt
		Int_t pdgCand =4122;
		Int_t pdgDaughter[3]={-1,-1,-1};
		pdgDaughter[0]=2212;
		pdgDaughter[1]=321;
		pdgDaughter[2]=211;   

		labDp = part->MatchToMC(pdgCand,arrayMC,3,pdgDaughter);
		if(labDp>=0){
			AliAODMCParticle *partDp = (AliAODMCParticle*)arrayMC->At(labDp);
			Int_t pdgMom=fVertUtil->CheckOrigin(arrayMC,partDp,fKeepLcNotFromQuark ? kFALSE : kTRUE);
			fhRecoPDGmom->Fill(pdgMom);
			if(pdgMom == 4) fIsLc=1;
			else if(pdgMom == 5) fIsLc=2;
			else fIsLc=0;
			Bool_t dummy = kTRUE;
			if(fIsLc>0) SetLambdacDaugh(partDp, arrayMC, dummy);
		}

	}
}

//---------------------------

void AliAnalysisTaskSELambdacTMVA::FillMassHists(AliAODEvent *aod, AliAODRecoDecayHF3Prong *d, TClonesArray *arrayMC, Int_t selection, Int_t selectionProb) {
	/// fill mass hists
	Bool_t IsInjected = 0;
	Double_t invMassLcpKpi = d->InvMassLcpKpi();
	Double_t invMassLcpiKp = d->InvMassLcpiKp();
	Bool_t ispKpiMC =	0;
	Bool_t ispiKpMC = 0;

	if(fReadMC) {
		ispKpiMC = IspKpiMC(d,arrayMC); 
		ispiKpMC = IspiKpMC(d,arrayMC); 
		AliAODMCHeader *mcHeader2 = (AliAODMCHeader*)aod->GetList()->FindObject(AliAODMCHeader::StdBranchName());
		if(!fIsHijing) IsInjected = fVertUtil->IsCandidateInjected(d,mcHeader2,arrayMC); //for dedicated MC set to 0

		//signal
		if(fIsLc>=1) {
			//MC PID
			fhMCmassLcPtSig->Fill(d->Pt(),ispKpiMC ? invMassLcpKpi : invMassLcpiKp);
			if(fIsLc==1) fhMCmassLcPtSigc->Fill(d->Pt(),ispKpiMC ? invMassLcpKpi : invMassLcpiKp);
			else if(fIsLc==2) fhMCmassLcPtSigb->Fill(d->Pt(),ispKpiMC ? invMassLcpKpi : invMassLcpiKp);
			//Real PID
			if(selection==1){
				fhPIDmassLcPtSig->Fill(d->Pt(),invMassLcpKpi);
				if(fIsLc==1) fhPIDmassLcPtSigc->Fill(d->Pt(),invMassLcpKpi);
				else if(fIsLc==2) fhPIDmassLcPtSigb->Fill(d->Pt(),invMassLcpKpi);
			}
			else if(selection==2){
				fhPIDmassLcPtSig->Fill(d->Pt(),invMassLcpiKp);
				if(fIsLc==1) fhPIDmassLcPtSigc->Fill(d->Pt(),invMassLcpiKp);
				else if(fIsLc==2) fhPIDmassLcPtSigb->Fill(d->Pt(),invMassLcpiKp);
			}
			else if(selection==3){
				fhPIDmassLcPtSig->Fill(d->Pt(),invMassLcpiKp);
				fhPIDmassLcPtSig->Fill(d->Pt(),invMassLcpKpi);
				if(fIsLc==1) {
					fhPIDmassLcPtSigc->Fill(d->Pt(),invMassLcpiKp);
					fhPIDmassLcPtSigc->Fill(d->Pt(),invMassLcpKpi);
				}
				else if(fIsLc==2) {
					fhPIDmassLcPtSigb->Fill(d->Pt(),invMassLcpiKp);
					fhPIDmassLcPtSigb->Fill(d->Pt(),invMassLcpKpi);
				}
			}
			//max prob PID
			if(selectionProb==1) {
				fhProbmassLcPtSig->Fill(d->Pt(), invMassLcpKpi);
				if(fIsLc==1)fhProbmassLcPtSigc->Fill(d->Pt(), invMassLcpKpi);
				if(fIsLc==2)fhProbmassLcPtSigb->Fill(d->Pt(), invMassLcpKpi);
			}
			else if(selectionProb==2) {
				fhProbmassLcPtSig->Fill(d->Pt(), invMassLcpiKp);
				if(fIsLc==1)fhProbmassLcPtSigc->Fill(d->Pt(), invMassLcpiKp);
				if(fIsLc==2)fhProbmassLcPtSigb->Fill(d->Pt(), invMassLcpiKp);
			}

			// mis id'd signal candidates
			if(ispKpiMC && selection==2){ //RDHF PID
			 	fhPtMisIdpKpi->Fill(fCandidateVars[0]);
				fhInvMassMisIdpKpi->Fill(fCandidateVars[0],invMassLcpiKp);
			}
			else if(ispiKpMC && selection==1){
			 	fhPtMisIdpiKp->Fill(fCandidateVars[0]);
				fhInvMassMisIdpiKp->Fill(fCandidateVars[0],invMassLcpKpi);
			}
			else fhPtCorrId->Fill(fCandidateVars[0]);

			if(ispKpiMC && selectionProb==2){ //PID using max prob. criteria
			 	fhPtMisIdpKpiProb->Fill(fCandidateVars[0]);
				fhInvMassMisIdpKpiProb->Fill(fCandidateVars[0],invMassLcpiKp);
			}
			else if(ispiKpMC && selectionProb==1){
			 	fhPtMisIdpiKpProb->Fill(fCandidateVars[0]);
				fhInvMassMisIdpiKpProb->Fill(fCandidateVars[0],invMassLcpKpi);
			}
			else fhPtCorrIdProb->Fill(fCandidateVars[0]);

		}
	}
	//Bkg
	if(!fReadMC || (fIsLc==0 && (!IsInjected || fSyst ==0))) { //data or non-injected background or pp background
		//MC PID
		if(fReadMC) fhMCmassLcPt->Fill(d->Pt(),ispKpiMC ? invMassLcpKpi : invMassLcpiKp);
		//Real PID
		if(selection==1){
			fhPIDmassLcPt->Fill(d->Pt(),invMassLcpKpi);
		}
		else if(selection==2){
			fhPIDmassLcPt->Fill(d->Pt(),invMassLcpiKp);
		}
		else if(selection==3){
			fhPIDmassLcPt->Fill(d->Pt(),invMassLcpiKp);
			fhPIDmassLcPt->Fill(d->Pt(),invMassLcpKpi);
		}
		//max prob PID
		if(selectionProb==1) fhProbmassLcPt->Fill(d->Pt(),invMassLcpKpi);
		else if(selectionProb==2)	fhProbmassLcPt->Fill(d->Pt(),invMassLcpiKp);
	}
}

//---------------------------

void AliAnalysisTaskSELambdacTMVA::FillNtuple(AliAODEvent *aod,AliAODRecoDecayHF3Prong *part,
																							TClonesArray *arrayMC, Int_t selection)
{
	//
	/// Function to fill NTuple with candidate's variables
	//

	if(fReadMC && !fKeepBkgNt && fIsLc==0) return; //don't continue if bkg and only want signal

	Bool_t IsInjected   = 0;
	Bool_t IsLc		= 0;
	Bool_t IsLcfromLb	= 0;

	if(fReadMC){ 
		AliAODMCHeader *mcHeader3 = (AliAODMCHeader*)aod->GetList()->FindObject(AliAODMCHeader::StdBranchName());
		if(!fIsHijing) IsInjected = fVertUtil->IsCandidateInjected(part,mcHeader3,arrayMC); //for dedicated MC set to 0
	}
	if(fIsLc>=1 && fIsLc<=2) IsLc=kTRUE;
	if(fIsLc==2) IsLcfromLb=kTRUE;
	if(fReadMC && IsInjected && !IsLc && fSyst >=1 ) return; //dont fill if injected bkg, pPb or PbPb

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
	Float_t tmp[58];
	//Is Lc
	if(!IsInjected && IsLc==0) tmp[0]=0; //non-injected bkg
	else if(IsLc==1 && !IsLcfromLb) tmp[0]=1; //prompt Lc
	else if(IsLc==1 && IsLcfromLb) tmp[0]=2; //non-prompt Lc
	else if(IsInjected && IsLc==0) tmp[0]=3; //injected bkg
	else tmp[0]=-99; //should not happen

	//invariant mass
	tmp[2]=part->InvMassLcpiKp();
	tmp[1]=part->InvMassLcpKpi(); 
	tmp[3]=part->Charge();
	//pt of decay products
	tmp[4]=part->PtProng(0);
	tmp[5]=part->PtProng(1); //pt kaon
	tmp[6]=part->PtProng(2);
	Float_t ptLc=part->Pt();
	tmp[7]=ptLc;
	tmp[8]=part->CosPointingAngle();
	tmp[9]=part->DecayLength();
	tmp[10]=part->NormalizedDecayLength();
	tmp[11]=TMath::Min(part->GetDist12toPrim(),part->GetDist23toPrim());
	tmp[12]=part->GetSigmaVert();
	Double_t dcas[3]={0};
	part->GetDCAs(dcas);
	tmp[13]=TMath::Max(dcas[0],TMath::Max(dcas[1],dcas[2]));
	tmp[14]=part->DecayLengthXY();
	tmp[15]=part->NormalizedDecayLengthXY();

	//Check resonant decays for efficiencies
	tmp[16]=fIsLcResonant; // bkg
	//PID selection
	tmp[17]=selection;
	//FONLL weights

	AliVTrack *track0=dynamic_cast<AliVTrack*>(part->GetDaughter(0));
	AliVTrack *track1=dynamic_cast<AliVTrack*>(part->GetDaughter(1));
	AliVTrack *track2=dynamic_cast<AliVTrack*>(part->GetDaughter(2));
  // PID
  // fill w -1
  for(Int_t iprob=18;iprob<=26;iprob++) {
    tmp[iprob]=-1;
  }
  if(fRDCutsAnalysis->GetIsUsePID() ) {
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
      tmp[18]=prob0[AliPID::kPion];			//track 0, pion
      tmp[19]=prob0[AliPID::kKaon];     		//kaon
      tmp[20]=prob0[AliPID::kProton];			//proton
      tmp[21]=prob1[AliPID::kPion];			//track 1, pion		
      tmp[22]=prob1[AliPID::kKaon];     		//kaon
      tmp[23]=prob1[AliPID::kProton];			//proton
      tmp[24]=prob2[AliPID::kPion];			//track 2, pion
      tmp[25]=prob2[AliPID::kKaon];     		//kaon
      tmp[26]=prob2[AliPID::kProton];			//proton
    }
  }
	if(fFillNtuple>=2) { //fill with further variables
		Double_t d00 = part->Getd0Prong(0);
		Double_t d01 = part->Getd0Prong(1);
		Double_t d02 = part->Getd0Prong(2);
		Double_t d0err0 = part->Getd0errProng(0);
		Double_t d0err1 = part->Getd0errProng(1);
		Double_t d0err2 = part->Getd0errProng(2);
		tmp[27]=d00;
		tmp[28]=d01;
		tmp[29]=d02;
		tmp[30]=d00*d00 + d01*d01 + d02*d02;
		tmp[31]=d00/d0err0;
		tmp[32]=d01/d0err1;
		tmp[33]=d02/d0err2;
		Double_t dd0,edd0;
		Double_t dd1,edd1;
		Double_t dd2,edd2;
		part->Getd0MeasMinusExpProng(0,aod->GetMagneticField(),dd0,edd0);
		part->Getd0MeasMinusExpProng(1,aod->GetMagneticField(),dd1,edd1);
		part->Getd0MeasMinusExpProng(2,aod->GetMagneticField(),dd2,edd2);
		Double_t ns0=dd0/edd0;
		Double_t ns1=dd1/edd1;
		Double_t ns2=dd2/edd2;
		tmp[34]=ns0;
		tmp[35]=ns1;
		tmp[36]=ns2;
		tmp[37]=part->CosPointingAngleXY();
		tmp[13]=dcas[0];
		tmp[38]=dcas[1];
		tmp[39]=dcas[2];
		tmp[11]=part->GetDist12toPrim();
		tmp[40]=part->GetDist23toPrim();
		Int_t runNumber=aod->GetRunNumber();
		tmp[41]=(Float_t)runNumber;

		if(fFillNtuple>=3) { //fill with weights
			Double_t ptLcMC=-1;
			//D meson weights
			Double_t weightPythia=-1,weight7LHC10f6a=-1,weight5LHC10f6a=-1,weight5LHC13d3=-1;	
			//Lc weights
			Double_t weight7LHC10f7aLc=-1,weight7LHC11b2Lc=-1,weight5LHC13d3Lc=-1;	
			if(IsLc) {
				//Get MC Lc to get true pT
				Int_t labDp = -1;
				Int_t pdgCand =4122;
				Int_t pdgDaughter[3]={-1,-1,-1};
				pdgDaughter[0]=2212;
				pdgDaughter[1]=321;
				pdgDaughter[2]=211;   
				labDp = part->MatchToMC(pdgCand,arrayMC,3,pdgDaughter);
				if(labDp>=0){
					AliAODMCParticle *motherPart = (AliAODMCParticle*)arrayMC->At(labDp);
					ptLcMC = motherPart->Pt();
					weightPythia = fFuncWeightPythia->Eval(ptLcMC);
					weight7LHC10f6a = fFuncWeightFONLL7overLHC10f6a->Eval(ptLcMC);
					weight5LHC10f6a = fFuncWeightFONLL5overLHC10f6a->Eval(ptLcMC);
					weight5LHC13d3 = fFuncWeightFONLL5overLHC13d3->Eval(ptLcMC);
					weight5LHC13d3Lc = fFuncWeightFONLL5overLHC13d3Lc->Eval(ptLcMC);
					weight7LHC11b2Lc = fFuncWeightFONLL7overLHC11b2Lc->Eval(ptLcMC);
					weight7LHC10f7aLc = fFuncWeightFONLL7overLHC10f7aLc->Eval(ptLcMC);
				}
			}
			tmp[42]=ptLcMC;
			tmp[43]=weightPythia;
			tmp[44]=weight7LHC10f6a;
			tmp[45]=weight5LHC10f6a;
			tmp[46]=weight5LHC13d3;
			tmp[47]=weight5LHC13d3Lc;
			tmp[48]=weight7LHC11b2Lc;
			tmp[49]=weight7LHC10f7aLc;

			//Multiplicity weights
			Double_t multWeight = 1.;
			Double_t nTracklets = 0.;
			if(fUseNchWeight){
				//tracklets within |eta| < 1.
				nTracklets = static_cast<Int_t>(fVertUtil->GetNumberOfTrackletsInEtaRange(aod,-1.,1.));
				multWeight *= GetNchWeight(static_cast<Int_t>(nTracklets));
				AliDebug(2,Form("Using Nch weights, Mult=%f Weight=%f\n",nTracklets,multWeight));
			}
			tmp[50]=nTracklets;
			tmp[51]=multWeight;


			if(fFillNtuple>=4) { //fill with invariant mass of 2 prongs
				tmp[52]=part->InvMass2Prongs(1,0,321,2212); //inv mass pK
				tmp[53]=part->InvMass2Prongs(2,1,211,321); //inv mass Kpi
				tmp[54]=part->InvMass2Prongs(2,0,211,2212);//inv mass ppi
				tmp[55]=part->InvMass2Prongs(1,2,321,2212); //inv mass Kp 
				tmp[56]=part->InvMass2Prongs(0,1,211,321); //inv mass piK
				tmp[57]=part->InvMass2Prongs(0,2,211,2212);//inv mass pip
			}
		}
	}
	fNtupleLambdac->Fill(tmp);
	PostData(5,fNtupleLambdac);

	return;
}

//---------------------------

void AliAnalysisTaskSELambdacTMVA::FillEffHists(Int_t kStep) {

	//
	/// Fill histograms (pt, pt vs eta, pt vs y, pt vs phi with
	/// candidates passing each step of the analysis
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
	/// Apply MC PID
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
	/// Apply MC PID
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

void AliAnalysisTaskSELambdacTMVA::FillSelectionBits(AliAODRecoDecayHF3Prong *d, TH2F *hSelectionBits) {
	if(d->HasSelectionBit(AliRDHFCuts::kD0toKpiCuts)) hSelectionBits->Fill(1,fCandidateVars[0]);
	if(d->HasSelectionBit(AliRDHFCuts::kD0toKpiPID)) hSelectionBits->Fill(2,fCandidateVars[0]);
	if(d->HasSelectionBit(AliRDHFCuts::kD0fromDstarCuts)) hSelectionBits->Fill(3,fCandidateVars[0]);
	if(d->HasSelectionBit(AliRDHFCuts::kD0fromDstarPID)) hSelectionBits->Fill(4,fCandidateVars[0]);
	if(d->HasSelectionBit(AliRDHFCuts::kDplusCuts)) hSelectionBits->Fill(5,fCandidateVars[0]);
	if(d->HasSelectionBit(AliRDHFCuts::kDplusPID)) hSelectionBits->Fill(6,fCandidateVars[0]);
	if(d->HasSelectionBit(AliRDHFCuts::kDsCuts)) hSelectionBits->Fill(7,fCandidateVars[0]);
	if(d->HasSelectionBit(AliRDHFCuts::kDsPID)) hSelectionBits->Fill(8,fCandidateVars[0]);
	if(d->HasSelectionBit(AliRDHFCuts::kLcCuts)) hSelectionBits->Fill(9,fCandidateVars[0]);
	if(d->HasSelectionBit(AliRDHFCuts::kLcPID)) hSelectionBits->Fill(10,fCandidateVars[0]);
	if(d->HasSelectionBit(AliRDHFCuts::kDstarCuts)) hSelectionBits->Fill(11,fCandidateVars[0]);
	if(d->HasSelectionBit(AliRDHFCuts::kDstarPID)) hSelectionBits->Fill(12,fCandidateVars[0]);
}

//-----------------------------

Int_t AliAnalysisTaskSELambdacTMVA::GetPIDselectionMaxProb(AliAODRecoDecayHF3Prong *part) {

	Int_t selection = 0;

	if(fRDCutsAnalysis->GetPidHF()->GetUseCombined()) {// get triplet identity based on the maximum probability
		AliVTrack *track0=dynamic_cast<AliVTrack*>(part->GetDaughter(0));
		AliVTrack *track1=dynamic_cast<AliVTrack*>(part->GetDaughter(1));
		AliVTrack *track2=dynamic_cast<AliVTrack*>(part->GetDaughter(2));
		//bayesian probabilities
		Double_t prob0[AliPID::kSPECIES];
		Double_t prob1[AliPID::kSPECIES];
		Double_t prob2[AliPID::kSPECIES];

		if (!track0 || !track1 || !track2) {
			AliError("AliVTrack missing");
			return 0;
		}
		fRDCutsAnalysis->GetPidHF()->GetPidCombined()->ComputeProbabilities(track0,fRDCutsAnalysis->GetPidHF()->GetPidResponse(),prob0);
		fRDCutsAnalysis->GetPidHF()->GetPidCombined()->ComputeProbabilities(track1,fRDCutsAnalysis->GetPidHF()->GetPidResponse(),prob1);
		fRDCutsAnalysis->GetPidHF()->GetPidCombined()->ComputeProbabilities(track2,fRDCutsAnalysis->GetPidHF()->GetPidResponse(),prob2);
		if(prob0[AliPID::kProton] * prob1[AliPID::kKaon] * prob2[AliPID::kPion] > prob2[AliPID::kProton] * prob1[AliPID::kKaon] * prob0[AliPID::kPion]) selection = 1; // pKpi
		else selection = 2; // piKp
	}
	return selection;
}

Double_t AliAnalysisTaskSELambdacTMVA::GetNchWeight(Int_t nch){
  //
  //  calculates the Nch weight using the measured and generateed Nch distributions
  //
  if(nch<=0) return 0.;
  if(!fHistoMCNch) { AliError("Input histos to evaluate Nch weights missing"); return 0.; }
  Double_t pMC=fHistoMCNch->GetBinContent(fHistoMCNch->FindBin(nch));
  Double_t weight = pMC;
  return weight;
}

void AliAnalysisTaskSELambdacTMVA::FillRecoNtuple(AliAODEvent *aod,AliAODRecoDecayHF3Prong *part,
																							TClonesArray *arrayMC) 
{
	//
	/// Function to fill NTuple after reco stage with candidates variables
	/// including TOF information for each daughter track for Geant3/Fluka
	/// correction needed in pp
	//

	Bool_t IsInjected   = 0;
	Bool_t IsLc		= 0;
	Bool_t IsLcfromLb	= 0;

	if(fReadMC){ 
		AliAODMCHeader *mcHeader3 = (AliAODMCHeader*)aod->GetList()->FindObject(AliAODMCHeader::StdBranchName());
		if(!fIsHijing) IsInjected = fVertUtil->IsCandidateInjected(part,mcHeader3,arrayMC); //for dedicated MC set to 0
	}
	if(fIsLc>=1 && fIsLc<=2) IsLc=kTRUE;
	if(fIsLc==2) IsLcfromLb=kTRUE;
	//if(fReadMC && IsInjected && !IsLc && fSyst >=1 ) return; //dont fill if injected bkg, pPb or PbPb

	Float_t tmp[16];
	//Is Lc
	if(!IsInjected && IsLc==0) tmp[0]=0; //non-injected bkg
	else if(IsLc==1 && !IsLcfromLb) tmp[0]=1; //prompt Lc
	else if(IsLc==1 && IsLcfromLb) tmp[0]=2; //non-prompt Lc
	else if(IsInjected && IsLc==0) tmp[0]=3; //injected bkg
	else tmp[0]=-99; //should not happen
	tmp[1]=part->Pt();
	//Get MC Lc to get true pT
	Float_t ptLcMC=0;
	Int_t labDp = -1;
	Int_t pdgCand =4122;
	Int_t pdgDaughter[3]={-1,-1,-1};
	pdgDaughter[0]=2212;
	pdgDaughter[1]=321;
	pdgDaughter[2]=211;   
	labDp = part->MatchToMC(pdgCand,arrayMC,3,pdgDaughter);
	if(labDp>=0){
		AliAODMCParticle *motherPart = (AliAODMCParticle*)arrayMC->At(labDp);
		ptLcMC = motherPart->Pt();
	}
	tmp[2]=ptLcMC;

	tmp[3]=part->PtProng(0);
	tmp[4]=part->PtProng(1); //pt kaon
	tmp[5]=part->PtProng(2);
	// Is the track good for TOF PID (for pp proton GEANT3/FLUKA correction)
	// and MC pt of each track
	Bool_t isTOFpid[3];
	Float_t MCpt[3];
	for(Int_t i=0;i<3;i++) {
		AliAODTrack *daugh=(AliAODTrack*)part->GetDaughter(i);
		Int_t daughLab= daugh->GetLabel();
		if(daughLab<0) continue;
		else{
			AliAODMCParticle* pdaugh = dynamic_cast<AliAODMCParticle*>(arrayMC->At(daughLab));
			isTOFpid[i] = fRDCutsAnalysis->GetIsUsePID() ? fRDCutsAnalysis->GetPidHF()->CheckTOFPIDStatus(daugh) : kFALSE;
			MCpt[i] = pdaugh->Pt();
		}
	}
	tmp[6]=MCpt[0];
	tmp[7]=MCpt[1];
	tmp[8]=MCpt[2];
	tmp[9]=isTOFpid[0];
	tmp[10]=isTOFpid[1];
	tmp[11]=isTOFpid[2];

	// Selection
	Int_t selectionCand=fRDCutsAnalysis->IsSelected(part,AliRDHFCuts::kCandidate,aod);
	Int_t selectionPID=fRDCutsAnalysis->IsSelected(part,AliRDHFCuts::kPID,aod);
  Int_t selectionPIDprob = 0;
  if(fRDCutsAnalysis->GetIsUsePID()) selectionPIDprob = GetPIDselectionMaxProb(part);
	tmp[12]=selectionCand;
	tmp[13]=selectionPID;
	tmp[14]=selectionPIDprob;
	tmp[15]=part->Charge();

//	Bool_t ispKpiMC = IspKpiMC(part,arrayMC); 
//	Bool_t ispiKpMC = IspiKpMC(part,arrayMC); 

//	tmp[15] = ispKpiMC?1:ispiKpMC?2:0;

	fNtupleLambdacReco->Fill(tmp);
	PostData(6,fNtupleLambdacReco);
}
