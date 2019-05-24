/**************************************************************************
 * Copyright(c) 1998-2018, ALICE Experiment at CERN, All rights reserved. *
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

/* AliAnalysisTaskDiHadCorrelHighPt
 * The task selects candidates for K0s, Lambdas and AntiLambdas (trigger particles)
 * and calculates correlations with charged unidentified particles (associated particles) in phi and eta.
 * The charged unidentified particles are also taken as trigger particles to have a check.
 * The task works with AOD (with or without MC info) events only and containes also mixing for acceptance corrections.
 * Last update edited by Lucia Anna Husova, January 2019
 */

#include <TChain.h>
#include <TH1F.h>
#include <TH3F.h>
#include <TH1D.h>
#include <TList.h>
#include <AliAnalysisTask.h>
#include <AliAnalysisManager.h>
#include "AliAODTrack.h"
#include <AliAODEvent.h>
#include "AliAODv0.h"
#include <AliAODInputHandler.h>
#include "AliAnalysisTaskDiHadCorrelHighPt.h"
#include <TMath.h>
#include <AliMultiEventInputHandler.h>
#include <AliPIDResponse.h>
#include <AliAODPid.h>
#include <THnSparse.h>
#include <AliAODVertex.h>
#include <AliEventPoolManager.h>
#include "AliAODMCParticle.h"
#include "AliMultSelection.h"
#include "AliPPVsMultUtils.h"

class AliPIDResponse;
class AliMultSelection;
class AliAnalysisTaskDiHadCorrelHighPt;    // analysis class

using namespace std;            // std namespace: so you can do things like 'cout'

ClassImp(AliAnalysisTaskDiHadCorrelHighPt) // classimp: necessary for root

AliAnalysisTaskDiHadCorrelHighPt::AliAnalysisTaskDiHadCorrelHighPt() : AliAnalysisTaskSE(),
    fAliEventCuts(),
    fAOD(0),
    fPIDResponse(0),
    fOutputList(0),
    fHistLambdaMassPtCut(0),
    fHistK0MassPtCut(0),
    fHistAntiLambdaMassPtCut(0),
    fHistPtHard(0),
	fHistKorelacie(0),
    fHistdPhidEtaMix(0),
    fHistV0Multiplicity(0),
    fHistMultVtxz(0),
    fHistMCPtAs(0),
    fHistRCPtAs(0),
    fHistNumberOfTriggers(0),
    fHistMCMixingRec(0),
    fFillMixed(0),
    fMixingTracks(5000),
	fPoolMgr(0x0),
    fPool(0x0),
    fAnalysisMC(kTRUE),
    fOStatus(0),
    fPtTrigMin(3),
    fPtAsocMin(1),
    fMixedEvents(20),
    fV0Radius(0.5),
    fSigmaCut(3.),
    fEtaCut(0.8),
    fRapidityCut(0.5),
    fLifeTimeLam(30),
    fLifeTimeK0(20),
    fMassRejectCutK0(0.005),
    fMassRejectCutLam(0.01),
    fRejectEventPileUp(kTRUE),
    fRejectTrackPileUp(kTRUE),
    fHistKorelacieMCrec(0),
    fHistNumberOfTriggersGen(0),
    fHistNumberOfTriggersRec(0),
    fHistRecV0(0),
    fHistGenV0(0),
    fHistMCPtTrigg(0),
    fHistRCPtTrigg(0),
    fHistSelection(0),
    fHistMultipPercentile(0),
    fHistTopolCut(0),
    fHistTopolCutMC(0),
    fHistPurityCheck(0),
    fCosPointAngleK0(0.97),
    fCosPointAngleLam(0.995),
    fDCAV0Daughters(1),
    fDCAposDaughter(0.06),
    fDCAnegDaughter(0.06),
    fnumOfTPCcrossedRows(70),
    fEfficiency(kTRUE),
    fPurityCheck(kTRUE),
    fCorrelations(kTRUE),
    fHistPtResolution(0),
    fNumberOfPtBinsTrigger(12),
    fNumberOfPtBinsAssoc(14),
    fHistV0MultiplicityK0(0),
    fHistV0Lam(0),
    fHistMultiplicityALam(0),
    fHitsNTracks(0),
    fHistPhiEta(0),
    fRejectTOF(kTRUE),
    fRejectV0PileUp(kTRUE),
    fMultEstimator("V0M"),
    fCorrelationsGen(kTRUE),
    fV0hCorr(kTRUE),
    fhhCorr(kTRUE),
    fhV0Corr(kTRUE),
    fFilterBit(768),
    fRemoveLamhFromCascade(kTRUE),
    fRemoveHadrFromV0(kTRUE),
    fAacceptLambdasFromCasscade(kTRUE),
    fPurePrimHadrons(kFALSE),
    fPureV0(kFALSE),
    fHistDeltaEtaDeltaPhiLamFineBinningLowPt(0),
    fMergingCut(0.02),
    fNumberOfDeltaPhiBins(72),
    fNumberOfDeltaEtaBins(75),
    fNumberOfEtaBins(10),
    fMixing(kTRUE)
{
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
AliAnalysisTaskDiHadCorrelHighPt::AliAnalysisTaskDiHadCorrelHighPt(const char* name, Bool_t analysisMC) : AliAnalysisTaskSE(name),
    fAliEventCuts(),
    fAOD(0),
    fPIDResponse(0),
    fOutputList(0),
    fHistLambdaMassPtCut(0),
    fHistK0MassPtCut(0),
    fHistAntiLambdaMassPtCut(0),
    fHistPtHard(0),
	fHistKorelacie(0),
    fHistdPhidEtaMix(0),
    fHistV0Multiplicity(0),
    fHistMultVtxz(0),
    fHistMCPtAs(0),
    fHistRCPtAs(0),
    fHistNumberOfTriggers(0),
    fHistMCMixingRec(0),
    fFillMixed(kTRUE),
    fMixingTracks(5000),
	fPoolMgr(0x0),
    fPool(0x0),
    fAnalysisMC(analysisMC),
    fOStatus(0),
    fPtTrigMin(3),
    fPtAsocMin(1),
    fMixedEvents(20),
    fV0Radius(0.5),
    fSigmaCut(3.),
    fEtaCut(0.8),
    fRapidityCut(0.5),
    fLifeTimeLam(30),
    fLifeTimeK0(20),
    fMassRejectCutK0(0.005),
    fMassRejectCutLam(0.01),
    fRejectEventPileUp(kTRUE),
    fRejectTrackPileUp(kTRUE),
    fHistKorelacieMCrec(0),
    fHistNumberOfTriggersGen(0),
    fHistNumberOfTriggersRec(0),
    fHistRecV0(0),
    fHistGenV0(0),
    fHistMCPtTrigg(0),
    fHistRCPtTrigg(0),
    fHistSelection(0),
    fHistMultipPercentile(0),
    fHistTopolCut(0),
    fHistTopolCutMC(0),
    fHistPurityCheck(0),
    fCosPointAngleK0(0.97),
    fCosPointAngleLam(0.995),
    fDCAV0Daughters(1),
    fDCAposDaughter(0.06),
    fDCAnegDaughter(0.06),
    fnumOfTPCcrossedRows(70),
    fEfficiency(kTRUE),
    fPurityCheck(kTRUE),
    fCorrelations(kTRUE),
    fHistPtResolution(0),
    fNumberOfPtBinsTrigger(12),
    fNumberOfPtBinsAssoc(14),
    fHistV0MultiplicityK0(0),
    fHistV0Lam(0),
    fHistMultiplicityALam(0),
    fHitsNTracks(0),
    fHistPhiEta(0),
    fRejectTOF(kTRUE),
    fRejectV0PileUp(kTRUE),
    fMultEstimator("V0M"),
    fCorrelationsGen(kTRUE),
    fV0hCorr(kTRUE),
    fhhCorr(kTRUE),
    fhV0Corr(kTRUE),
    fFilterBit(768),
    fRemoveLamhFromCascade(kTRUE),
    fRemoveHadrFromV0(kTRUE),
    fAacceptLambdasFromCasscade(kTRUE),
    fPurePrimHadrons(kFALSE),
    fPureV0(kFALSE),
    fHistDeltaEtaDeltaPhiLamFineBinningLowPt(0),
    fMergingCut(0.02),
    fNumberOfDeltaPhiBins(72),
    fNumberOfDeltaEtaBins(75),
    fNumberOfEtaBins(10),
    fMixing(kTRUE)
{
    // constructor
    DefineInput(0, TChain::Class());    // define the input of the analysis: in this case we take a 'chain' of events
                                        // this chain is created by the analysis manager, so no need to worry about it, 
                                        // it does its work automatically
    DefineOutput(1, TList::Class());    // define the ouptut of the analysis: in this case it's a list of histograms 
                                        // you can add more output objects by calling DefineOutput(2, classname::Class())
                                        // if you add more output objects, make sure to call PostData for all of them, and to
                                        // make changes to your AddTask macro!
}
//_____________________________________________________________________________
AliAnalysisTaskDiHadCorrelHighPt::~AliAnalysisTaskDiHadCorrelHighPt()
{
    // destructor
    if(fOutputList) {
        delete fOutputList;     // at the end of your task, it is deleted from memory by calling this function
    }
}
//_____________________________________________________________________________
void AliAnalysisTaskDiHadCorrelHighPt::UserCreateOutputObjects()
{
	const Double_t kPi = TMath::Pi();
	Float_t kPtBins[15] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
	Int_t kMassBins = 500;
	Float_t kMassMinK = 0.4;
	Float_t kMassMaxK = 0.6;

	Float_t kMassMinLambda = 0.8;
	Float_t kMassMaxLambda = 1.2;

	Float_t kMassBinsK[kMassBins];
	Float_t kMassBinsLambda[kMassBins];

	kMassBinsK[0] = kMassMinK;
	kMassBinsLambda[0] = kMassMinLambda;
	
	for (Int_t i=0; i<kMassBins; i++){
		kMassBinsK[i+1] = kMassBinsK[i] + (kMassMaxK-kMassMinK)/kMassBins;
		kMassBinsLambda[i+1] = kMassBinsLambda[i] + (kMassMaxLambda-kMassMinLambda)/kMassBins;
	}

	Int_t kNCuts = 11;
	Float_t kCuts[kNCuts];
	kCuts[0]=0;
	for (Int_t i=0; i<kNCuts; i++){
		kCuts[i+1]=kCuts[i]+1;
	}

    Int_t bins[10]= {fNumberOfPtBinsTrigger,fNumberOfPtBinsAssoc,fNumberOfDeltaPhiBins,fNumberOfDeltaEtaBins,9,7,fNumberOfEtaBins,fNumberOfEtaBins,601,10};
    Double_t min[10] = {fPtTrigMin,fPtAsocMin, -kPi/2, -2., -10., 0.,-0.8,-0.8,0.44,0};
    Double_t max[10] = {15., 15., -kPi/2+2*kPi, 2., 10., 7.,0.8,0.8, 1.15,100};
    
	Int_t  NofCentBins  = 10;
    Double_t MBins[]={0,10,20,30,40,50,60,70,80,90,100};
         
    Int_t NofZVrtxBins  = 9;
    Double_t ZBins[10]={-10.0, -7., -5.0, -3., -1.0, 1., 3.0, 5., 7., 10.};

	Int_t bins2d[6] = {fNumberOfPtBinsTrigger,9,fNumberOfEtaBins,5,601,10};
	Double_t mis2d[6] = {fPtTrigMin,-10,-0.8,0.,0.44,0};
	Double_t maxs2d[6] = {15.,10,0.8,5.,1.15,100};
	
    // create output objects
    //
    // this function is called ONCE at the start of your analysis (RUNTIME)
    // here you ceate the histograms that you want to use 
    //
    // the histograms are in this case added to a tlist, this list is in the end saved
    // to an output file
    //
    fOutputList = new TList();          // this is a list which will contain all of your histograms
                                        // at the end of the analysis, the contents of this list are written
                                        // to the output file
    fOutputList->SetOwner(kTRUE);       // memory stuff: the list is owner of all objects it contains and will delete them
                                        // if requested (dont worry about this now)
        
    fHistK0MassPtCut = new TH3F ("fHistK0MassPtCut", "fHistK0MassPtCut", kMassBins, kMassBinsK, 14, kPtBins, kNCuts, kCuts);
    fHistLambdaMassPtCut = new TH3F("fHistLambdaMassPtCut", "fHistLambdaMassPtCut", kMassBins, kMassBinsLambda, 14, kPtBins, kNCuts, kCuts);
	fHistAntiLambdaMassPtCut = new TH3F("fHistAntiLambdaMassPtCut", "fHistAntiLambdaMassPtCut", kMassBins, kMassBinsLambda, 14, kPtBins, kNCuts, kCuts);
    fOutputList->Add(fHistK0MassPtCut); 
    fOutputList->Add(fHistLambdaMassPtCut); 
	fOutputList->Add(fHistAntiLambdaMassPtCut);
    
    fHistPtHard = new TH3F("fHistPtHard","fHistPtHard",100,0,1,10,0,100,4,0,4);
    fOutputList->Add(fHistPtHard);

	fHistKorelacie = new THnSparseF ("fHistKorelacie","fHistKorelacie", 10, bins, min, max);
    fHistKorelacie->GetAxis(0)->SetTitle("p_{T}^{trig}");
    fHistKorelacie->GetAxis(1)->SetTitle("p_{T}^{assoc}");
    fHistKorelacie->GetAxis(2)->SetTitle("#Delta#phi");
    fHistKorelacie->GetAxis(3)->SetTitle("#Delta#eta");
    fHistKorelacie->GetAxis(4)->SetTitle("p_{vz}");
    fHistKorelacie->GetAxis(5)->SetTitle("trigger");
    fHistKorelacie->GetAxis(6)->SetTitle("#eta_{trig}");
    fHistKorelacie->GetAxis(7)->SetTitle("#eta_{assoc}");
    fHistKorelacie->GetAxis(8)->SetTitle("mass");
    fHistKorelacie->GetAxis(9)->SetTitle("multiplicity percentile");
    fOutputList->Add(fHistKorelacie);
    fHistKorelacie->Sumw2();
    
    Double_t *binsMass = new Double_t [602];
    binsMass[0]=0.44;
    for(Int_t i=0;i<602;i++){
        if(i<301) binsMass[i+1] = binsMass[i]+(0.56-0.44)/300;
        if(i==301) binsMass[i+1]=1.08;
        if(i>301) binsMass[i+1] = binsMass[i]+(1.15-1.08)/300;
        
    }
    fHistKorelacie->GetAxis(8)->Set(602,binsMass);
    fHistKorelacie->GetAxis(4)->Set(9,ZBins);
    
	fHistdPhidEtaMix = new THnSparseF ("fHistdPhidEtaMix", "fHistdPhidEtaMix", 10, bins, min, max);
    fHistdPhidEtaMix->GetAxis(0)->SetTitle("p_{T}^{trig}");
    fHistdPhidEtaMix->GetAxis(1)->SetTitle("p_{T}^{assoc}");
    fHistdPhidEtaMix->GetAxis(2)->SetTitle("#Delta#phi");
    fHistdPhidEtaMix->GetAxis(3)->SetTitle("#Delta#eta");
    fHistdPhidEtaMix->GetAxis(4)->SetTitle("p_{vz}");
    fHistdPhidEtaMix->GetAxis(5)->SetTitle("trigger");
    fHistdPhidEtaMix->GetAxis(6)->SetTitle("#eta_{trig}");
    fHistdPhidEtaMix->GetAxis(7)->SetTitle("#eta_{assoc}");
    fHistdPhidEtaMix->GetAxis(8)->SetTitle("mass");
    fHistdPhidEtaMix->GetAxis(9)->SetTitle("multiplicity percentile");
    fHistdPhidEtaMix->Sumw2();
	fOutputList->Add(fHistdPhidEtaMix);
    fHistdPhidEtaMix->GetAxis(4)->Set(9,ZBins);
    fHistdPhidEtaMix->GetAxis(8)->Set(602,binsMass);
    
    fHistMCMixingRec = new THnSparseF ("fHistMCMixingRec", "fHistMCMixingRec", 10, bins, min, max);
    fHistMCMixingRec->GetAxis(0)->SetTitle("p_{T}^{trig}");
    fHistMCMixingRec->GetAxis(1)->SetTitle("p_{T}^{assoc}");
    fHistMCMixingRec->GetAxis(2)->SetTitle("#Delta#phi");
    fHistMCMixingRec->GetAxis(3)->SetTitle("#Delta#eta");
    fHistMCMixingRec->GetAxis(4)->SetTitle("p_{vz}");
    fHistMCMixingRec->GetAxis(5)->SetTitle("trigger");
    fHistMCMixingRec->GetAxis(6)->SetTitle("#eta_{trig}");
    fHistMCMixingRec->GetAxis(7)->SetTitle("#eta_{assoc}");
    fHistMCMixingRec->GetAxis(8)->SetTitle("mass");
    fHistMCMixingRec->GetAxis(9)->SetTitle("multiplicity percentile");
    fOutputList->Add(fHistMCMixingRec);
    fHistMCMixingRec->Sumw2();
    fHistMCMixingRec->GetAxis(4)->Set(9,ZBins);
    fHistMCMixingRec->GetAxis(8)->Set(602,binsMass);

    fHistMCKorelacie = new THnSparseF ("fHistMCKorelacie","fHistMCKorelacie", 10, bins, min, max);
    fHistMCKorelacie->GetAxis(0)->SetTitle("p_{T}^{trig}");
    fHistMCKorelacie->GetAxis(1)->SetTitle("p_{T}^{assoc}");
    fHistMCKorelacie->GetAxis(2)->SetTitle("#Delta#phi");
    fHistMCKorelacie->GetAxis(3)->SetTitle("#Delta#eta");
    fHistMCKorelacie->GetAxis(4)->SetTitle("p_{vz}");
    fHistMCKorelacie->GetAxis(5)->SetTitle("trigger");
    fHistMCKorelacie->GetAxis(6)->SetTitle("#eta_{trig}");
    fHistMCKorelacie->GetAxis(7)->SetTitle("#eta_{assoc}");
    fHistMCKorelacie->GetAxis(8)->SetTitle("mass");
    fHistMCKorelacie->GetAxis(9)->SetTitle("multiplicity percentile");
    fHistMCKorelacie->Sumw2();
    fOutputList->Add(fHistMCKorelacie);
    
    fHistMCKorelacie->GetAxis(8)->Set(602,binsMass);
    fHistMCKorelacie->GetAxis(4)->Set(9,ZBins);
    
    fHistKorelacieMCrec = new THnSparseF ("fHistKorelacieMCrec","fHistKorelacieMCrec", 10, bins, min, max);
    fHistKorelacieMCrec->GetAxis(0)->SetTitle("p_{T}^{trig}");
    fHistKorelacieMCrec->GetAxis(1)->SetTitle("p_{T}^{assoc}");
    fHistKorelacieMCrec->GetAxis(2)->SetTitle("#Delta#phi");
    fHistKorelacieMCrec->GetAxis(3)->SetTitle("#Delta#eta");
    fHistKorelacieMCrec->GetAxis(4)->SetTitle("p_{vz}");
    fHistKorelacieMCrec->GetAxis(5)->SetTitle("trigger");
    fHistKorelacieMCrec->GetAxis(6)->SetTitle("#eta_{trig}");
    fHistKorelacieMCrec->GetAxis(7)->SetTitle("#eta_{assoc}");
    fHistKorelacieMCrec->GetAxis(8)->SetTitle("mass");
    fHistKorelacieMCrec->GetAxis(9)->SetTitle("multiplicity percentile");
    fHistKorelacieMCrec->Sumw2();
    fOutputList->Add(fHistKorelacieMCrec);
    fHistKorelacieMCrec->GetAxis(8)->Set(602,binsMass);
    fHistKorelacieMCrec->GetAxis(4)->Set(9,ZBins);

	fHistV0Multiplicity = new TH1D ("fHistV0Multiplicity", "fHistV0Multiplicity", 100, 0, 100);
	fOutputList->Add(fHistV0Multiplicity);
    
    fHistV0MultiplicityK0 = new TH1F ("fHistV0MultiplicityK0", "fHistV0MultiplicityK0", 100, 0, 100);
    fOutputList->Add(fHistV0MultiplicityK0);
    
    fHistV0Lam = new TH1F ("fHistV0Lam", "fHistV0Lam", 100, 0, 100);
    fOutputList->Add(fHistV0Lam);
   
    fHistMultiplicityALam = new TH1F ("fHistMultiplicityALam","fHistMultiplicityALam",100,0,100);
    fOutputList->Add(fHistMultiplicityALam);

	fHistMultVtxz = new TH2D ("fHistMultVtxz","fHistMultVtxz",NofCentBins,MBins,NofZVrtxBins,ZBins);
	fOutputList->Add(fHistMultVtxz);

	fHistMCPtAs = new TH3D("fHistMCPtAs","fHistMCPtAs",fNumberOfPtBinsAssoc,fPtAsocMin,15,9,-10,10,fNumberOfEtaBins,-0.8,0.8);
	fOutputList->Add(fHistMCPtAs);
    fHistMCPtAs->Sumw2();
    fHistMCPtAs->GetYaxis()->Set(9,ZBins);
	fHistRCPtAs = new TH3D("fHistRCPtAs","fHistRCPtAs",fNumberOfPtBinsAssoc,fPtAsocMin,15,9,-10,10,fNumberOfEtaBins,-0.8,0.8);
	fOutputList->Add(fHistRCPtAs);
    fHistRCPtAs->Sumw2();
    fHistRCPtAs->GetYaxis()->Set(9,ZBins);
    
    fHistMCPtTrigg = new TH3D("fHistMCPtTrigg","fHistMCPtTrigg",fNumberOfPtBinsTrigger,fPtTrigMin,15,9,-10,10,fNumberOfEtaBins,-0.8,0.8);
    fOutputList->Add(fHistMCPtTrigg);
    fHistMCPtTrigg->Sumw2();
    fHistMCPtTrigg->GetYaxis()->Set(9,ZBins);
    fHistRCPtTrigg = new TH3D("fHistRCPtTrigg","fHistRCPtTrigg",fNumberOfPtBinsTrigger,fPtTrigMin,15,9,-10,10,fNumberOfEtaBins,-0.8,0.8);
    fOutputList->Add(fHistRCPtTrigg);
    fHistRCPtTrigg->Sumw2();
    fHistRCPtTrigg->GetYaxis()->Set(9,ZBins);
    
    Int_t binsTrig[4]={fNumberOfPtBinsAssoc,9,3,fNumberOfEtaBins};
    Double_t mintrig[4]={fPtAsocMin,-10,0,-0.8};
    Double_t maxtrig[4]={15,10,3,0.8};
    fHistGenV0 = new THnSparseF("fHistGenV0","fHistGenV0",4,binsTrig,mintrig,maxtrig);
    fOutputList->Add(fHistGenV0);
    fHistGenV0->Sumw2();
    fHistGenV0->GetAxis(1)->Set(9,ZBins);
    Int_t binsTrigRec[5]={fNumberOfPtBinsAssoc,9,3,fNumberOfEtaBins,602};
    Double_t mintrigRec[6]={fPtAsocMin,-10,0,-0.8,0.44};
    Double_t maxtrigRec[6]={15,10,3,0.8,1.15};
    fHistRecV0 = new THnSparseF("fHistRecV0","fHistRecV0",5,binsTrigRec,mintrigRec,maxtrigRec);
    fOutputList->Add(fHistRecV0);
    fHistRecV0->Sumw2();
    fHistRecV0->GetAxis(0)->SetTitle("p_{T}");
    fHistRecV0->GetAxis(1)->SetTitle("p_{vz}");
    fHistRecV0->GetAxis(2)->SetTitle("trigger");
    fHistRecV0->GetAxis(3)->SetTitle("#eta");
    fHistRecV0->GetAxis(4)->SetTitle("mass");
    fHistRecV0->GetAxis(4)->Set(602,binsMass);
    fHistRecV0->GetAxis(1)->Set(9,ZBins);

	fHistNumberOfTriggers = new THnSparseF("fHistNumberOfTriggers","fHistNumberOfTriggers",6,bins2d,mis2d, maxs2d);
    fHistNumberOfTriggers->GetAxis(0)->SetTitle("p_{T}");
    fHistNumberOfTriggers->GetAxis(1)->SetTitle("p_{vz}");
    fHistNumberOfTriggers->GetAxis(2)->SetTitle("#eta");
    fHistNumberOfTriggers->GetAxis(3)->SetTitle("trigger");
    fHistNumberOfTriggers->GetAxis(4)->SetTitle("mass");
    fHistNumberOfTriggers->GetAxis(5)->SetTitle("multiplicity percentile");
	fOutputList->Add(fHistNumberOfTriggers);
    fHistNumberOfTriggers->Sumw2();
    fHistNumberOfTriggers->GetAxis(4)->Set(602,binsMass);
    fHistNumberOfTriggers->GetAxis(1)->Set(9,ZBins);
    fHistNumberOfTriggersGen = new THnSparseF("fHistNumberOfTriggersGen","fHistNumberOfTriggersGen",6,bins2d,mis2d, maxs2d);
    fHistNumberOfTriggersGen->GetAxis(0)->SetTitle("p_{T}");
    fHistNumberOfTriggersGen->GetAxis(1)->SetTitle("p_{vz}");
    fHistNumberOfTriggersGen->GetAxis(2)->SetTitle("#eta");
    fHistNumberOfTriggersGen->GetAxis(3)->SetTitle("trigger");
    fHistNumberOfTriggersGen->GetAxis(4)->SetTitle("mass");
    fHistNumberOfTriggersGen->GetAxis(5)->SetTitle("multiplicity percentile");
    fOutputList->Add(fHistNumberOfTriggersGen);
    fHistNumberOfTriggersGen->Sumw2();
    fHistNumberOfTriggersGen->GetAxis(4)->Set(602,binsMass);
    fHistNumberOfTriggersGen->GetAxis(1)->Set(9,ZBins);
    
    fHistNumberOfTriggersRec = new THnSparseF("fHistNumberOfTriggersRec","fHistNumberOfTriggersRec",6,bins2d,mis2d,maxs2d);
    fHistNumberOfTriggersRec->GetAxis(0)->SetTitle("p_{T}");
    fHistNumberOfTriggersRec->GetAxis(1)->SetTitle("p_{vz}");
    fHistNumberOfTriggersRec->GetAxis(2)->SetTitle("#eta");
    fHistNumberOfTriggersRec->GetAxis(3)->SetTitle("trigger");
    fHistNumberOfTriggersRec->GetAxis(4)->SetTitle("mass");
    fHistNumberOfTriggersRec->GetAxis(5)->SetTitle("multiplicity percentile");
    fOutputList->Add(fHistNumberOfTriggersRec);
    fHistNumberOfTriggersRec->Sumw2();
    fHistNumberOfTriggersRec->GetAxis(4)->Set(602,binsMass);
    fHistNumberOfTriggersRec->GetAxis(1)->Set(9,ZBins);

    fHistSelection = new TH1D("fHistSelection","fHistSelection",4,0,4);
    fOutputList->Add(fHistSelection);
    
    fHistMultipPercentile = new TH1F("fHistMultipPercentile","fHistMultipPercentile",10,0,100);
    fOutputList->Add(fHistMultipPercentile);
    fHistMultipPercentile->Sumw2();
    
    Int_t binsCuts[11] = {12,601,100,100,20,100,500,200,100,3,2};
    Double_t binsMinCuts[11] = {fPtTrigMin,0.4,0.03,0.03,0,0,0.95,0.,0,0,0};
    Double_t binsMaxCuts[11] = {15,1.2,0.25,0.25,2,1.5,1,50.,0.09,3,2};
    
    fHistTopolCut = new THnSparseF("fHistTopolCut","fHistTopolCut",11,binsCuts,binsMinCuts,binsMaxCuts);
    fHistTopolCutMC = new THnSparseF("fHistTopolCutMC","fHistTopolCutMC",11,binsCuts,binsMinCuts,binsMaxCuts);
    fOutputList->Add(fHistTopolCut);
    fOutputList->Add(fHistTopolCutMC);
    fHistTopolCut->Sumw2();
    fHistTopolCutMC->Sumw2();
    
    fHistTopolCut->GetAxis(0)->SetTitle("p_{T}");
    fHistTopolCut->GetAxis(1)->SetTitle("mass");
    fHistTopolCut->GetAxis(2)->SetTitle("DCA neg");
    fHistTopolCut->GetAxis(3)->SetTitle("DCA pos");
    fHistTopolCut->GetAxis(4)->SetTitle("DCA daughters");
    fHistTopolCut->GetAxis(5)->SetTitle("V0 radius");
    fHistTopolCut->GetAxis(6)->SetTitle("cos PA");
    fHistTopolCut->GetAxis(7)->SetTitle("life time");
    fHistTopolCut->GetAxis(8)->SetTitle("mass selection");
    fHistTopolCut->GetAxis(9)->SetTitle("V0 type");
    fHistTopolCut->GetAxis(10)->SetTitle("OnFly/Offline");
    fHistTopolCut->GetAxis(1)->Set(602,binsMass);
    
    fHistTopolCutMC->GetAxis(0)->SetTitle("p_{T}");
    fHistTopolCutMC->GetAxis(1)->SetTitle("mass");
    fHistTopolCutMC->GetAxis(2)->SetTitle("DCA neg");
    fHistTopolCutMC->GetAxis(3)->SetTitle("DCA pos");
    fHistTopolCutMC->GetAxis(4)->SetTitle("DCA daughters");
    fHistTopolCutMC->GetAxis(5)->SetTitle("V0 radius");
    fHistTopolCutMC->GetAxis(6)->SetTitle("cos PA");
    fHistTopolCutMC->GetAxis(7)->SetTitle("life time");
    fHistTopolCutMC->GetAxis(8)->SetTitle("mass selection");
    fHistTopolCutMC->GetAxis(9)->SetTitle("V0 type");
    fHistTopolCutMC->GetAxis(10)->SetTitle("OnFly/Offline"); // 0.5 - "On-The-Fly", 1.5 Offline
    fHistTopolCutMC->GetAxis(1)->Set(602,binsMass);
    
    Int_t binsPur[6] = {fNumberOfPtBinsTrigger,602,4,8,8,40};
    Double_t binsPurMin[6] = {fPtTrigMin,0.44,0.,0.,0,-0.8};
    Double_t binsPurMax[6] = {15,1.15,4.,8.,8.,0.8};
    fHistPurityCheck = new THnSparseF("fHistPurityCheck","fHistPurityCheck",6,binsPur,binsPurMin,binsPurMax);
    fOutputList->Add(fHistPurityCheck);
    fHistPurityCheck->Sumw2();
    fHistPurityCheck->GetAxis(0)->SetTitle("p_{T}");
    fHistPurityCheck->GetAxis(1)->SetTitle("mass");
    fHistPurityCheck->GetAxis(2)->SetTitle("V0 type");
    fHistPurityCheck->GetAxis(3)->SetTitle("check");
    fHistPurityCheck->GetAxis(5)->SetTitle("#eta");
    fHistPurityCheck->GetAxis(1)->Set(602,binsMass);
    
    fHistPtResolution = new TH3F("fHistPtResol","fHistPtResol",144,0,18,144,0,18,4,0,4);
    fOutputList->Add(fHistPtResolution);
    
    fHitsNTracks = new TH2F("fHitsNTracks","fHitsNTracks",100,0,100,2,0,2);
    fOutputList->Add(fHitsNTracks);
    
    Int_t binsPhiEta[4]= {fNumberOfPtBinsAssoc,72,40,4};
    Double_t minsPhiEta[4] = {fPtAsocMin,0,-0.8,0};
    Double_t maxsphiEta[4] = {15,2*kPi,0.8,4};
    
    fHistPhiEta= new THnSparseF("fHistPhiEta","fHistPhiEta",4,binsPhiEta,minsPhiEta,maxsphiEta);
    fOutputList->Add(fHistPhiEta);

    fHistDeltaEtaDeltaPhiLamFineBinningLowPt = new TH2D("fHistDeltaEtaDeltaPhiLamFineBinningLowPt","fHistDeltaEtaDeltaPhiLamFineBinningLowPt",140 ,-0.7,0.7, 140,-0.7,0.7);
    fHistDeltaEtaDeltaPhiLamFineBinningLowPt->GetXaxis()->SetTitle("#Delta#phi");
    fHistDeltaEtaDeltaPhiLamFineBinningLowPt->GetYaxis()->SetTitle("#Delta#eta");
    fOutputList->Add(fHistDeltaEtaDeltaPhiLamFineBinningLowPt);
    fHistDeltaEtaDeltaPhiLamFineBinningLowPt->Sumw2();

    PostData(1, fOutputList);           // postdata will notify the analysis manager of changes / updates to the 
                                        // fOutputList object. the manager will in the end take care of writing your output to file
                                        // so it needs to know what's in the output
	// Settings for event mixing -------------------------------------
    Int_t trackDepth = fMixingTracks;
    Int_t poolSize   = 200;  // Maximum number of events, ignored in the present implemented of AliEventPoolManager
  
     fPoolMgr = new AliEventPoolManager(poolSize, trackDepth, NofCentBins, MBins, NofZVrtxBins, ZBins);
     fPoolMgr->SetTargetValues(trackDepth, 0.1, 5);
	//----------------------------------------------
}
//_____________________________________________________________________________
void AliAnalysisTaskDiHadCorrelHighPt::UserExec(Option_t *)
{
    // user exec
    // this function is called once for each event
    // the manager will take care of reading the events from file, and with the static function InputEvent() you 
    // have access to the current event. 
    // once you return from the UserExec function, the manager will retrieve the next event from the chain
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    AliAODInputHandler *inEvMain = (AliAODInputHandler *) mgr->GetInputEventHandler();

    fAOD = (AliAODEvent *) inEvMain -> GetEvent();    // get an event (called fAOD) from the input file
    if(!fAOD) return;                                   // if the pointer to the ,event is empty (getting it failed) skip this event
    
    fPIDResponse = (AliPIDResponse *) inEvMain-> GetPIDResponse();
    
    // physics selection
    fHistSelection->Fill(0.5);
	UInt_t maskIsSelected = inEvMain->IsEventSelected();

	//  data trigger selection
	Bool_t isSelected = (maskIsSelected & AliVEvent::kINT7); //pp
	if (!isSelected) return;
    fHistSelection->Fill(1.5);

	AliAODVertex *myPrimVertex = fAOD->GetPrimaryVertex();
	if (!myPrimVertex) return;
	Double_t lPVz = myPrimVertex->GetZ();

	if (TMath::Abs(lPVz)>=10) return;
    fHistSelection->Fill(2.5);
	Double_t lPVx = myPrimVertex->GetX();
	Double_t lPVy = myPrimVertex->GetY();

    if(!fAnalysisMC&&fRejectEventPileUp){
        cout << "event cuts" << endl;
        fAliEventCuts.SetupRun2pp();
        if(!fAliEventCuts.AcceptEvent(fAOD)) return;
    }
    fHistSelection->Fill(3.5);
    
    Int_t iTracks(fAOD->GetNumberOfTracks());           // see how many tracks there are in the event
    Int_t nV0(fAOD->GetNumberOfV0s());                  // see how many V0 there are in the event

	fHistV0Multiplicity->Fill(nV0);

	// Multiplicity definition
    Double_t lPercentile = 302;
    AliMultSelection *MultSelection = 0x0;
    MultSelection = (AliMultSelection * ) fAOD->FindListObject("MultSelection");
    if( !MultSelection) {
        AliWarning("AliMultSelection object not found!");
    }else{
        if(fMultEstimator=="V0M") lPercentile = MultSelection->GetMultiplicityPercentile("V0M");
        else if(fMultEstimator=="V0A") lPercentile = MultSelection->GetMultiplicityPercentile("V0A");
        else if(fMultEstimator=="SPDTracklets") lPercentile = MultSelection->GetMultiplicityPercentile("SPDTracklets");
        else if(fMultEstimator=="RefMult05") lPercentile = MultSelection->GetMultiplicityPercentile("RefMult05");
        else if(fMultEstimator=="RefMult08") lPercentile = MultSelection->GetMultiplicityPercentile("RefMult08");
    }
    if ((lPercentile<0.)||(lPercentile>100.)) return;
    
    fHistMultipPercentile->Fill(lPercentile);

    TObjArray *mcTracksSel = new TObjArray; // generated associated particles
    mcTracksSel->SetOwner(kTRUE);
    TObjArray *mcTracksTrigSel = new TObjArray; // generated trigger charged hadrons
    mcTracksTrigSel->SetOwner(kTRUE);
    TObjArray *mcTracksV0Sel = new TObjArray; // Generated V0 triggers
    mcTracksV0Sel->SetOwner(kTRUE);
    TObjArray *mcV0AssocSel = new TObjArray; // Generated V0 assoc
    mcTracksV0Sel->SetOwner(kTRUE);
    TObjArray *selectedMCassoc = new TObjArray; // all reconstructed associated particles, with reconstructed pt,phi,eta values - for raw correlation function
    selectedMCassoc->SetOwner(kTRUE);
    TObjArray *selectedMCV0assoc = new TObjArray; // all reconstructed V0 as associated particles, with reconstructed pt,phi,eta values - for raw correlation function
    selectedMCV0assoc->SetOwner(kTRUE);
    TObjArray *selectedMCtrig= new TObjArray; // all reconstructed trigger particles, with reconstructed pt,phi,eta values - for raw correlation function
    selectedMCtrig->SetOwner(kTRUE);
    TObjArray *selectedMCV0Triggersrec = new TObjArray;  // All reconstructed V0 candidates for triggers with reconstructed pt,phi,eta values - for raw correlation function
    selectedMCV0Triggersrec->SetOwner(kTRUE);
    TClonesArray *mcArray = new TClonesArray;
    mcArray->SetOwner(kTRUE);

	//=========== MC loop ===============================
    Double_t ptHard = 0.;
    
	if(fAnalysisMC){
        AliAODMCHeader *aodMCheader = (AliAODMCHeader*)fAOD->FindListObject(AliAODMCHeader::StdBranchName());
        ptHard = aodMCheader->GetPtHard();
        Float_t vzMC = aodMCheader->GetVtxZ();
        if (TMath::Abs(vzMC) >= 10.) return;
 	     //retrieve MC particles from event
        mcArray = (TClonesArray*)fAOD->FindListObject(AliAODMCParticle::StdBranchName());
        if(!mcArray){
            Printf("No MC particle branch found");
            return;
        }

		Int_t nMCAllTracks = mcArray->GetEntriesFast();
	
		for (Int_t i = 0; i < nMCAllTracks; i++){
			AliAODMCParticle *mcTrack = (AliAODMCParticle*)mcArray->At(i);
            if (!mcTrack) {
                Error("ReadEventAODMC", "Could not receive particle %d", i);
                continue;
            }
 			// track cuts for generated particles
			Double_t mcTrackEta = mcTrack->Eta();
			Double_t mcTrackPt = mcTrack->Pt();
			Bool_t TrIsPrim = mcTrack->IsPhysicalPrimary();
			Bool_t TrEtaMax = TMath::Abs(mcTrackEta)<0.8;
			Bool_t TrPtMin = mcTrackPt>fPtAsocMin; 
			Bool_t TrCharge = (mcTrack->Charge())!=0;
            
			if (TrIsPrim && TrPtMin && TrCharge && TrEtaMax) {
                
                if(fEfficiency) fHistMCPtAs->Fill(mcTrackPt,lPVz,mcTrackEta); // for recunstruction efficiency calculation
                if(fCorrelationsGen) mcTracksSel->Add(new AliV0ChParticle(mcTrack->Eta(),mcTrack->Phi(),mcTrack->Pt(),4,mcTrack->GetLabel(),mcTrack->GetLabel()));
                
                if (mcTrackPt>fPtTrigMin) {
                    if(fCorrelationsGen) mcTracksTrigSel->Add(new AliV0ChParticle(mcTrack->Eta(),mcTrack->Phi(),mcTrack->Pt(),4,mcTrack->GetLabel()));
                    if(fEfficiency) fHistMCPtTrigg->Fill(mcTrackPt,lPVz,mcTrackEta); // for recunstruction efficiency calculation
                }
            }

            //--- MC closure test - selection of V0 ----

            Int_t mcPartPdg = mcTrack->GetPdgCode();
            
            if ((mcPartPdg != 310) && (mcPartPdg != 3122) && (mcPartPdg != (-3122))) continue; // keep only Lambdas and K0S
            Bool_t IsFromCascade = kFALSE;
            
            Int_t mother  = mcTrack->GetMother();
            AliAODMCParticle *mcMotherParticle = static_cast<AliAODMCParticle*>(mcArray->At(mother));
            Int_t motherPDG = 0;
            if (mother<0) motherPDG =0;
            else motherPDG = TMath::Abs(mcMotherParticle->PdgCode());
            
           if(fAacceptLambdasFromCasscade) IsFromCascade = (((motherPDG == 3222)|| (motherPDG==3212)|| (motherPDG==3112) || (motherPDG==3224) || (motherPDG==3214) || (motherPDG==3114) || (motherPDG==3322) || (motherPDG==3312)|| (motherPDG==3324) || (motherPDG==3314) || (motherPDG==3334)) && (mcMotherParticle->IsPhysicalPrimary()));

            Bool_t IsK0 = mcPartPdg==310;
            Bool_t IsLambda = mcPartPdg==3122;
            Bool_t IsAntiLambda = mcPartPdg==-3122;
            
            Bool_t isPhysPrim = mcTrack->IsPhysicalPrimary();

            IsK0 = IsK0 && (isPhysPrim);
            IsLambda = IsLambda && (isPhysPrim||IsFromCascade);
            IsAntiLambda = IsAntiLambda && (isPhysPrim||IsFromCascade);
            
            AliAODMCParticle* daughter0 = 0x0;
            AliAODMCParticle* daughter1 = 0x0;
            Int_t dau0 = mcTrack->GetDaughterLabel(0);
            if (dau0>0) daughter0 = (AliAODMCParticle*) mcArray->At(dau0);
            Int_t dau1 = mcTrack->GetDaughterLabel(1);
            if (dau1>0) daughter1 = (AliAODMCParticle*) mcArray->At(dau1);
            
            if(!daughter0||!daughter1) continue;
            
            Int_t labelPos = -1;
            Int_t labelNeg = -1;
        
            if(daughter0->Charge()<0){
                labelPos = daughter1->GetLabel();
                labelNeg = daughter0->GetLabel();
            }
            if(daughter0->Charge()>0) {
                labelPos = daughter0->GetLabel();
                labelNeg = daughter1->GetLabel();
            }else{
                labelPos = daughter0->GetLabel();
                labelNeg = daughter1->GetLabel();
            }

            Double_t etaDau0 = daughter0->Eta();
            Double_t etaDau1 = daughter1->Eta();
            //MC V0 cuts
            Double_t V0genrapidity = mcTrack->Y();

            if (mcTrack->Pt()>fPtAsocMin&&TMath::Abs(V0genrapidity)<0.5&&TMath::Abs(etaDau0)<0.8&&TMath::Abs(etaDau1)<0.8){
                if(IsK0) {
                    if(fCorrelationsGen) mcV0AssocSel->Add(new AliV0ChParticle(mcTrack->Eta(),mcTrack->Phi(),mcTrack->Pt(),5,mcTrack->GetLabel(),labelPos,labelNeg,kFALSE,mcTrack->M()));
                    if (fEfficiency){
                        Double_t v0effic[4]={mcTrack->Pt(),lPVz,0.5,mcTrack->Eta()};
                        fHistGenV0->Fill(v0effic); // for recunstruction efficiency calculation
                    }
                }
                if(IsLambda) {
                    if(fCorrelationsGen) mcV0AssocSel->Add(new AliV0ChParticle(mcTrack->Eta(),mcTrack->Phi(),mcTrack->Pt(),6,mcTrack->GetLabel(),labelPos,labelNeg,kFALSE,mcTrack->M()));
                    if (fEfficiency){
                        Double_t v0effic[4]={mcTrack->Pt(),lPVz,1.5,mcTrack->Eta()};
                        fHistGenV0->Fill(v0effic); // for recunstruction efficiency calculation
                    }
                }
                if(IsAntiLambda) {
                    if(fCorrelationsGen) mcV0AssocSel->Add(new AliV0ChParticle(mcTrack->Eta(),mcTrack->Phi(),mcTrack->Pt(),7,mcTrack->GetLabel(),labelPos,labelNeg,kFALSE,mcTrack->M()));
                    if (fEfficiency){
                        Double_t v0effic[4]={mcTrack->Pt(),lPVz,2.5,mcTrack->Eta()};
                        fHistGenV0->Fill(v0effic); // for recunstruction efficiency calculation
                    }
                }
                
            }
            if (mcTrack->Pt()>fPtTrigMin&&fCorrelationsGen&&TMath::Abs(V0genrapidity)<0.5&&TMath::Abs(etaDau0)<0.8&&TMath::Abs(etaDau1)<0.8){
                if(IsK0) mcTracksV0Sel->Add(new AliV0ChParticle(mcTrack->Eta(),mcTrack->Phi(),mcTrack->Pt(),1,mcTrack->GetLabel(),labelPos,labelNeg,kFALSE,mcTrack->M()));
                if(IsLambda) mcTracksV0Sel->Add(new AliV0ChParticle(mcTrack->Eta(),mcTrack->Phi(),mcTrack->Pt(),2,mcTrack->GetLabel(),labelPos,labelNeg,kFALSE,mcTrack->M()));
                if(IsAntiLambda) mcTracksV0Sel->Add(new AliV0ChParticle(mcTrack->Eta(),mcTrack->Phi(),mcTrack->Pt(),3,mcTrack->GetLabel(),labelPos,labelNeg,kFALSE,mcTrack->M()));
            }

		}
        // MC closure test corellations
        if(fCorrelationsGen){
            //V0-h
            if(fV0hCorr) Corelations(mcTracksV0Sel,mcTracksSel,fHistMCKorelacie, lPVz, fHistNumberOfTriggersGen,kFALSE,kTRUE,lPercentile,fHistPtHard,ptHard,kFALSE,kTRUE);

            //h-h
            if(fhhCorr) Corelations(mcTracksTrigSel,mcTracksSel,fHistMCKorelacie, lPVz,fHistNumberOfTriggersGen, kFALSE, kFALSE,lPercentile,fHistPtHard,ptHard,kFALSE,kFALSE);
            //h-V0
            if(fhV0Corr) Corelations(mcTracksTrigSel,mcV0AssocSel,fHistMCKorelacie, lPVz,fHistNumberOfTriggersGen, kFALSE, kTRUE,lPercentile,fHistPtHard,ptHard,kTRUE,kFALSE);
        }
        //reconstructed part. 
		Int_t nTracks = fAOD->GetNumberOfTracks();

        for (Int_t i = 0; i < nTracks; i++)
        {
            AliAODTrack* tras = dynamic_cast<AliAODTrack*>(fAOD->GetTrack(i));
            if(!tras) {
                AliFatal("Not a standard AOD");
                continue;
            }

       		if ((tras->Pt())<fPtAsocMin) continue;
        	if (!(IsMyGoodPrimaryTrack(tras))) continue;
        	
            if ((tras->Charge())==0) continue;
            
            if(fRejectTrackPileUp&&fRejectTOF&&(!(tras->HasPointOnITSLayer(1) || tras->HasPointOnITSLayer(0) || tras->GetTOFBunchCrossing()==0 ))) continue; // track by track pile-up rejection
            if(fRejectTrackPileUp&&!fRejectTOF&&(!(tras->HasPointOnITSLayer(1) ||tras->HasPointOnITSLayer(0)))) continue; // track by track pile-up rejection
            
            Double_t mcPt = tras->Pt();
            Double_t mcPhi = tras->Phi();
            Double_t mcEta = tras->Eta();
            
            Double_t phiEta[4] = {mcPt,mcPhi,mcEta,3.5};
            fHistPhiEta->Fill(phiEta);
            Int_t AssocLabel = tras->GetLabel();
            
            if (AssocLabel<=0) continue;
            if(fCorrelations&&!fPurePrimHadrons) selectedMCassoc->Add(new AliV0ChParticle(mcEta,mcPhi,mcPt,4,AssocLabel,tras->GetID(),tras->Charge(),tras->Pz(),tras->E()));
            if (mcPt>fPtTrigMin) {
                if(fCorrelations&&!fPurePrimHadrons) selectedMCtrig->Add(new AliV0ChParticle(mcEta,mcPhi,mcPt,4,AssocLabel,tras->Charge(),tras->Pz(),tras->E()));
            }
            Double_t purhadr[6] = {mcPt,0,3.5,0.5,-1,mcEta};
            fHistPurityCheck->Fill(purhadr);
            
            AliAODMCParticle* mcTrack = static_cast<AliAODMCParticle*>(mcArray->At(AssocLabel));
            if(!mcTrack) continue;
            Bool_t isPhyPrim = mcTrack->IsPhysicalPrimary();
            Double_t genPt = mcTrack->Pt();
            Double_t genEta = mcTrack->Eta();
            
            if (isPhyPrim) {
                if(fCorrelations&&fPurePrimHadrons) selectedMCassoc->Add(new AliV0ChParticle(mcEta,mcPhi,mcPt,4,AssocLabel,tras->GetID(),tras->Charge(),tras->Pz(),tras->E()));
                if (mcPt>fPtTrigMin) {
                    if(fCorrelations&&fPurePrimHadrons) selectedMCtrig->Add(new AliV0ChParticle(mcEta,mcPhi,mcPt,4,AssocLabel,tras->Charge(),tras->Pz(),tras->E()));
                }
                Double_t purhadrPrim[6] = {mcPt,0,3.5,1.5,-1,mcEta};
                fHistPurityCheck->Fill(purhadrPrim);
                fHistPtResolution->Fill(genPt,mcPt,3.5);
                if(fEfficiency) fHistRCPtAs->Fill(genPt,lPVz,genEta); // for recunstruction efficiency calculation
            }
            if(mcTrack->IsSecondaryFromMaterial()){
                Double_t purhadrMater[6] = {mcPt,0,3.5,2.5,0.5,mcEta};
                fHistPurityCheck->Fill(purhadrMater);
            }
            if(mcTrack->IsSecondaryFromWeakDecay()){
                Double_t purhadrDecay[6] = {mcPt,0,3.5,2.5,1.5,mcEta};
                fHistPurityCheck->Fill(purhadrDecay);
            }
        }
    }

	//=========== end of MC loop ==========

	TObjArray * selectedTracks = new TObjArray;
	selectedTracks->SetOwner(kTRUE);

	TObjArray * selectedAssociatedTracks = new TObjArray;
	selectedAssociatedTracks->SetOwner(kTRUE);
	TObjArray * selectedTriggerTracks = new TObjArray;
	selectedTriggerTracks->SetOwner(kTRUE);
    Int_t nTrak = 0;
    Int_t nTrakBefore = 0;
    
    for(Int_t i=0; i < iTracks; i++) {                 // loop over all these tracks
        AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(i));         // get a track (type AliAODTrack) from the event
        if(!track) continue;                            // if we failed, skip this track
        
        if(!IsMyGoodPrimaryTrack(track)) continue; // hybrid track selection
        nTrakBefore+=1;
        if(fRejectTrackPileUp&&fRejectTOF&&(!(track->HasPointOnITSLayer(1) || track->HasPointOnITSLayer(0) || track->GetTOFBunchCrossing()==0 ))) continue; // track by track pile-up rejection using TOF
        if(fRejectTrackPileUp&&!fRejectTOF&&(!(track->HasPointOnITSLayer(1) ||track->HasPointOnITSLayer(0)))) continue; // track by track pile-up rejection without TOF information
        nTrak+=1;
        
		if(track->Pt()>fPtAsocMin) selectedTracks->Add(track);
        if(track->Pt()>fPtAsocMin&&!fAnalysisMC) {
            selectedAssociatedTracks-> Add(new AliV0ChParticle(track->Eta(), track->Phi(), track->Pt(), 4, 0,track->GetID(),track->Charge(),track->Pz(),track->E()));
            Double_t phiEtaData[4] = {track->Pt(),track->Phi(),track->Eta(),3.5};
            fHistPhiEta->Fill(phiEtaData);
        }
        if(track->Pt()>fPtTrigMin&&!fAnalysisMC) selectedTriggerTracks-> Add(new AliV0ChParticle(track->Eta(), track->Phi(), track->Pt(), 4,0,track->GetID(),track->Charge(),track->Pz(),track->E()));
	}
    fHitsNTracks->Fill(nTrakBefore,0.5);
    fHitsNTracks->Fill(nTrak,1.5);
    TObjArray * selectedV0 = new TObjArray;
	selectedV0->SetOwner(kTRUE); 
	TObjArray * selectedV0Triggers = new TObjArray;
	selectedV0Triggers->SetOwner(kTRUE); 
	TObjArray * selectedV0Assoc = new TObjArray;
	selectedV0Assoc->SetOwner(kTRUE);

    AliAODTrack *myTrackPos = 0x0;
    AliAODTrack *myTrackNeg = 0x0;
    Int_t nK0 =0;
    Int_t nLam =0;
    Int_t nALam =0;
	for (Int_t i=0; i<nV0; i++){
        
        AliAODv0* V0 = static_cast<AliAODv0*>(fAOD->GetV0(i));
        if(!V0) continue;
        
        Double_t massK0=V0->MassK0Short();
        Double_t massLambda=V0->MassLambda();
        Double_t massAntilambda=V0->MassAntiLambda();

        Double_t ptTrig = V0->Pt();
        if(V0->Pt()<fPtAsocMin) continue; // pt assoc cut

        Bool_t k0 = ((massK0>0.44)&&(massK0<0.56));
        Bool_t Lambda = ((massLambda>1.08)&&(massLambda<1.15));
        Bool_t Antilambda = ((massAntilambda>1.08)&&(massAntilambda<1.15));

        // PID cut--------------------------
        Float_t nSigmaPosPion   = 0.;
        Float_t nSigmaNegPion   = 0.;
        Float_t nSigmaPosProton = 0.;
        Float_t nSigmaNegProton = 0.;

        AliVTrack *trackNegTest=dynamic_cast<AliVTrack *>(V0->GetDaughter(1));
        AliVTrack *trackPosTest=dynamic_cast<AliVTrack *>(V0->GetDaughter(0));

        AliAODTrack * myTrackNegTest = (AliAODTrack*) trackNegTest;
        AliAODTrack * myTrackPosTest = (AliAODTrack*) trackPosTest;
  
        if (!myTrackPosTest || !myTrackNegTest) {
            Printf("strange analysis::UserExec:: Error:Could not retreive one of the daughter track\n");
            continue;
        }
 
        if( myTrackPosTest->Charge() ==1){
            myTrackPos = myTrackPosTest;
            myTrackNeg = myTrackNegTest;
        }
 
        if( myTrackPosTest->Charge() ==-1){
            myTrackPos = myTrackNegTest;
            myTrackNeg = myTrackPosTest;
        }
  
        const AliAODPid *pPid = myTrackPos->GetDetPid();
        const AliAODPid *nPid = myTrackNeg->GetDetPid();
        
        if (pPid){
            Double_t pdMom = pPid->GetTPCmomentum();
            if (pdMom<1.){
                nSigmaPosPion = fPIDResponse->NumberOfSigmasTPC(myTrackPos, AliPID::kPion);
                nSigmaPosProton = fPIDResponse->NumberOfSigmasTPC(myTrackPos, AliPID::kProton);
            }
        }
        
        
        if (nPid){
            Double_t ndMom = nPid->GetTPCmomentum();
            if (ndMom<1.){
                nSigmaNegPion = fPIDResponse->NumberOfSigmasTPC(myTrackNeg, AliPID::kPion);
                nSigmaNegProton = fPIDResponse->NumberOfSigmasTPC(myTrackNeg, AliPID::kProton);
            }
        }
        Bool_t bpPion = kTRUE;
        Bool_t bpProton = kTRUE;
        Bool_t bnPion = kTRUE;
        Bool_t bnProton = kTRUE;
        
        Bool_t cutK0Pid = kTRUE;
        Bool_t cutLambdaPid = kTRUE;
        Bool_t cutAntiLambdaPid = kTRUE;
        
        if(!fAnalysisMC){
            bpPion = TMath::Abs(nSigmaPosPion) <= fSigmaCut; // TPC dE/dx selection
            bpProton = TMath::Abs(nSigmaPosProton) <= fSigmaCut;
        
            bnPion = TMath::Abs(nSigmaNegPion) <= fSigmaCut;
            bnProton = TMath::Abs(nSigmaNegProton) <= fSigmaCut;
        }
        
        cutK0Pid = (bpPion && bnPion);
        cutLambdaPid = (bpProton && bnPion);
        cutAntiLambdaPid = (bpPion && bnProton);
        
        // reject bunch-off pile-up
        
        if (fRejectV0PileUp&&(!(((myTrackNeg->IsOn(AliAODTrack::kTPCrefit)&&myTrackNeg->IsOn(AliAODTrack::kITSrefit))||myTrackNeg->IsOn(AliAODTrack::kTOFout))||((myTrackPos->IsOn(AliAODTrack::kTPCrefit)&&myTrackPos->IsOn(AliAODTrack::kITSrefit))||myTrackPos->IsOn(AliAODTrack::kTOFout))))) continue;
        
        if (k0) {
            fHistK0MassPtCut->Fill(massK0,ptTrig,0.5);
            nK0+=1;
        }
        if (Lambda) {
            fHistLambdaMassPtCut->Fill(massLambda,ptTrig,0.5);
            nLam+=1;
        }
        if (Antilambda) {
            fHistAntiLambdaMassPtCut->Fill(massAntilambda,ptTrig,0.5);
            nALam+=1;
        }
        
        if (k0&&cutK0Pid) fHistK0MassPtCut->Fill(massK0,ptTrig,1.5);
        if (Lambda&&cutLambdaPid) fHistLambdaMassPtCut->Fill(massLambda,ptTrig,1.5);
        if (Antilambda&&cutAntiLambdaPid) fHistAntiLambdaMassPtCut->Fill(massAntilambda,ptTrig,1.5);
        
        Int_t oStatus = GetOStatus();
        if(!IsMyGoodV0(V0,myTrackPos,myTrackNeg,oStatus)) continue; // on fly and daughters cuts
        
        if (k0&&cutK0Pid) fHistK0MassPtCut->Fill(massK0,ptTrig,2.5);
        if (Lambda&&cutLambdaPid) fHistLambdaMassPtCut->Fill(massLambda,ptTrig,2.5);
        if (Antilambda&&cutAntiLambdaPid) fHistAntiLambdaMassPtCut->Fill(massAntilambda,ptTrig,2.5);
        
        //======= crosscheck topological cuts
        if(fCutsCrosscheck){
            Double_t DCANegToPV = V0->DcaNegToPrimVertex();
            Double_t DCAPosToPV = V0->DcaPosToPrimVertex();
            Double_t DCADaught = V0->DcaV0Daughters();
            Double_t V0radius = V0->RadiusV0();
            Double_t CosPA = V0->CosPointingAngle(myPrimVertex);
            Double_t massSellLam = TMath::Abs(massK0-0.497614);
            Double_t massSellK0 = TMath::Abs(massLambda-1.115683);
        
            Double_t *tParentVertexPosition = new Double_t[3];
            tParentVertexPosition[0]= lPVx;
            tParentVertexPosition[1]= lPVy;
            tParentVertexPosition[2]= lPVz;
        
            Double_t lenght = V0->DecayLengthV0(tParentVertexPosition);
            Double_t momentum = TMath::Sqrt(TMath::Power(V0->Px(),2)+TMath::Power(V0->Py(),2)+TMath::Power(V0->Pz(),2));
            
            delete[] tParentVertexPosition;

            Double_t lifetimeK0 = (massK0*lenght)/momentum;
            Double_t lifetimeLam = (massLambda*lenght)/momentum;
            Double_t lifetimeALam = (massAntilambda*lenght)/momentum;
        
            if(k0&&cutK0Pid&&IsMyGoodV0RapidityK0(V0)) {
                if (V0->GetOnFlyStatus()) TopologCuts(fHistTopolCut,ptTrig,massK0,DCANegToPV,DCAPosToPV,DCADaught,V0radius,CosPA,lifetimeK0,massSellK0,0.5,0.5);
                if (!V0->GetOnFlyStatus()) TopologCuts(fHistTopolCut,ptTrig,massK0,DCANegToPV,DCAPosToPV,DCADaught,V0radius,CosPA,lifetimeK0,massSellK0,0.5,1.5);
            }
        
            if (Lambda&&cutLambdaPid&&IsMyGoodV0RapidityLambda(V0)) {
                if (V0->GetOnFlyStatus()) TopologCuts(fHistTopolCut,ptTrig,massLambda,DCANegToPV,DCAPosToPV,DCADaught,V0radius,CosPA,lifetimeLam,massSellLam,1.5,0.5);
                if (!V0->GetOnFlyStatus()) TopologCuts(fHistTopolCut,ptTrig,massLambda,DCANegToPV,DCAPosToPV,DCADaught,V0radius,CosPA,lifetimeLam,massSellLam,1.5,1.5);
            }
          
            if (Antilambda&&cutAntiLambdaPid&&IsMyGoodV0RapidityLambda(V0)) {
                if (V0->GetOnFlyStatus()) TopologCuts(fHistTopolCut,ptTrig,massAntilambda,DCANegToPV,DCAPosToPV,DCADaught,V0radius,CosPA,lifetimeALam,massSellLam,2.5,0.5);
                if (!V0->GetOnFlyStatus()) TopologCuts(fHistTopolCut,ptTrig,massAntilambda,DCANegToPV,DCAPosToPV,DCADaught,V0radius,CosPA,lifetimeALam,massSellLam,2.5,1.5);
            }
        
            if(fAnalysisMC){
                AliAODTrack *daughter1= static_cast<AliAODTrack*> (V0->GetDaughter(1));
                AliAODTrack *daughter0= static_cast<AliAODTrack*> (V0->GetDaughter(0));
                if (!(daughter1->GetLabel()<0||daughter0->GetLabel()<0)){
            
                    Int_t mcmother1 = static_cast<AliAODMCParticle*>(mcArray->At(daughter1->GetLabel()))->GetMother();
                    Int_t mcmother0 = static_cast<AliAODMCParticle*>(mcArray->At(daughter0->GetLabel()))->GetMother();
            
                    if(mcmother1==mcmother0){
                        Int_t pdgV0 = static_cast<AliAODMCParticle*>(mcArray->At(mcmother0))->PdgCode();
                        Int_t pdgD1 = static_cast<AliAODMCParticle*>(mcArray->At(daughter1->GetLabel()))->PdgCode();
                        Int_t pdgD0 = static_cast<AliAODMCParticle*>(mcArray->At(daughter0->GetLabel()))->PdgCode();
            
                        Bool_t isPhyPrim = static_cast<AliAODMCParticle*>(mcArray->At(mcmother0))->IsPhysicalPrimary();
                        if(isPhyPrim) {
                            Double_t V0mcPt = static_cast<AliAODMCParticle*>(mcArray->At(mcmother0))->Pt();
                            if(V0mcPt>=fPtTrigMin) {
                                if (pdgV0==3122 && ((pdgD0==2212 && pdgD1==-211 )||(pdgD0==-211 && pdgD1==2212))&&IsMyGoodV0RapidityLambda(V0)) {
                                    if (V0->GetOnFlyStatus()) TopologCuts(fHistTopolCutMC,ptTrig,massLambda,DCANegToPV,DCAPosToPV,DCADaught,V0radius,CosPA,lifetimeLam,massSellLam,1.5,0.5);
                                    if (!V0->GetOnFlyStatus()) TopologCuts(fHistTopolCutMC,ptTrig,massLambda,DCANegToPV,DCAPosToPV,DCADaught,V0radius,CosPA,lifetimeLam,massSellLam,1.5,1.5);
                                }
            
                                if ((pdgV0==310) && (( pdgD0==211 && pdgD1==-211 )||( pdgD0==-211 && pdgD1==211 ))&&IsMyGoodV0RapidityK0(V0)){
                                    if (V0->GetOnFlyStatus()) TopologCuts(fHistTopolCutMC,ptTrig,massK0,DCANegToPV,DCAPosToPV,DCADaught,V0radius,CosPA,lifetimeK0,massSellK0,0.5,0.5);
                                    if (!V0->GetOnFlyStatus()) TopologCuts(fHistTopolCutMC,ptTrig,massK0,DCANegToPV,DCAPosToPV,DCADaught,V0radius,CosPA,lifetimeK0,massSellK0,0.5,1.5);
                                }
            
                                if ((pdgV0==-3122)&& ((pdgD0==-2212 && pdgD1==211 )||(pdgD0==211 && pdgD1==-2212 ))&&IsMyGoodV0RapidityLambda(V0)) {
                                    if (V0->GetOnFlyStatus()) TopologCuts(fHistTopolCutMC,ptTrig,massAntilambda,DCANegToPV,DCAPosToPV,DCADaught,V0radius,CosPA,lifetimeALam,massSellLam,2.5,0.5);
                                    if (!V0->GetOnFlyStatus()) TopologCuts(fHistTopolCutMC,ptTrig,massAntilambda,DCANegToPV,DCAPosToPV,DCADaught,V0radius,CosPA,lifetimeALam,massSellLam,2.5,1.5);
                                }
                            }
                        }
                    }
                }
            }
        }
        
        //------------------- V0 cuts -------------------------------
        
        if (k0&&cutK0Pid) fHistK0MassPtCut->Fill(massK0,ptTrig,3.5);
        if (Lambda&&cutLambdaPid) fHistLambdaMassPtCut->Fill(massLambda,ptTrig,3.5);
        if (Antilambda&&cutAntiLambdaPid) fHistAntiLambdaMassPtCut->Fill(massAntilambda,ptTrig,3.5);
        
        if(!IsMyGoodV0Topology(V0)) continue; //topoligical cuts

        if (k0&&cutK0Pid) fHistK0MassPtCut->Fill(massK0,ptTrig,4.5);
        if (Lambda&&cutLambdaPid) fHistLambdaMassPtCut->Fill(massLambda,ptTrig,4.5);
        if (Antilambda&&cutAntiLambdaPid) fHistAntiLambdaMassPtCut->Fill(massAntilambda,ptTrig,4.5);

        if (Lambda&&cutLambdaPid){
            
            if(IsMyGoodV0RapidityLambda(V0)){ //Rapidity
                fHistLambdaMassPtCut->Fill(massLambda,ptTrig,5.5);
            
                if(IsMyGoodLifeTimeLambda(lPVx,lPVy,lPVz,V0)) { //Proper Lifetime (mL/p)
                    fHistLambdaMassPtCut->Fill(massLambda,ptTrig,6.5);
            
                    if(IsMyGoodV0AngleLambda(V0,myPrimVertex)){ //V0 Cosine of Pointing Angle
                        fHistLambdaMassPtCut->Fill(massLambda,ptTrig,7.5);
            
                        if (TMath::Abs(massK0-0.497614)>fMassRejectCutLam){
                            fHistLambdaMassPtCut->Fill(massLambda,ptTrig,8.5);
            
                            if(fAnalysisMC){
                                FillMC(V0,mcArray,3122,2212, -211,2,massLambda,selectedMCV0Triggersrec,fHistRecV0,fHistLambdaMassPtCut,lPVz,myTrackPos,myTrackNeg,V0->GetOnFlyStatus(),fHistPurityCheck,selectedMCV0assoc,fHistPtResolution);
                                Double_t phiEtaLamMC[4] = {V0->Pt(),V0->Phi(),V0->Eta(),1.5};
                                fHistPhiEta->Fill(phiEtaLamMC);
                            }
                            if(!fAnalysisMC) {
                                if(V0->Pt()>fPtTrigMin) selectedV0Triggers-> Add(new AliV0ChParticle(V0->Eta(), V0->Phi(), V0->Pt(), 2,0,myTrackPos->GetID(),myTrackNeg->GetID(),V0->GetOnFlyStatus(),massLambda,myTrackPos->Phi(),myTrackPos->Pt(),myTrackPos->Eta(),myTrackPos->Charge(),myTrackNeg->Phi(),myTrackNeg->Pt(),myTrackNeg->Eta(),myTrackNeg->Charge()));
                                Double_t phiEtaLamData[4] = {V0->Pt(),V0->Phi(),V0->Eta(),1.5};
                                fHistPhiEta->Fill(phiEtaLamData);
                                selectedV0Assoc-> Add(new AliV0ChParticle(V0->Eta(), V0->Phi(), V0->Pt(), 6,0,myTrackPos->GetID(),myTrackNeg->GetID(),V0->GetOnFlyStatus(),massLambda,myTrackPos->Phi(),myTrackPos->Pt(),myTrackPos->Eta(),myTrackPos->Charge(),myTrackNeg->Phi(),myTrackNeg->Pt(),myTrackNeg->Eta(),myTrackNeg->Charge()));
                            }
                        }
                    }
                }
            }
        }
        
        if (k0&&cutK0Pid){
            if(IsMyGoodV0RapidityK0(V0)){ //Rapidity
                fHistK0MassPtCut->Fill(massK0,ptTrig,5.5);
            
                if(IsMyGoodLifeTimeK0(lPVx,lPVy,lPVz,V0)) { //Proper Lifetime (mL/p)
                    fHistK0MassPtCut->Fill(massK0,ptTrig,6.5);
            
                    if (IsMyGoodV0AngleK0(V0,myPrimVertex)){ //V0 Cosine of Pointing Angle
                        fHistK0MassPtCut->Fill(massK0,ptTrig,7.5);
            
                        if(TMath::Abs(massLambda-1.115683)>fMassRejectCutK0&&TMath::Abs(massAntilambda-1.115683)>fMassRejectCutK0) {
                            fHistK0MassPtCut->Fill(massK0,ptTrig,8.5);
            
                            if(fAnalysisMC){
                                FillMC(V0,mcArray,310,211, -211,1,massK0,selectedMCV0Triggersrec,fHistRecV0,fHistK0MassPtCut,lPVz,myTrackPos,myTrackNeg,V0->GetOnFlyStatus(),fHistPurityCheck,selectedMCV0assoc,fHistPtResolution);
                                Double_t phiEtaK0MC[4] = {V0->Pt(),V0->Phi(),V0->Eta(),0.5};
                                fHistPhiEta->Fill(phiEtaK0MC);
                            }
                            if(!fAnalysisMC) {
                                if(V0->Pt()>fPtTrigMin) selectedV0Triggers-> Add(new AliV0ChParticle(V0->Eta(), V0->Phi(), V0->Pt(), 1,0,myTrackPos->GetID(),myTrackNeg->GetID(),V0->GetOnFlyStatus(),massK0,myTrackPos->Phi(),myTrackPos->Pt(),myTrackPos->Eta(),myTrackPos->Charge(),myTrackNeg->Phi(),myTrackNeg->Pt(),myTrackNeg->Eta(),myTrackNeg->Charge()));
                                Double_t phiEtaK0Data[4] = {V0->Pt(),V0->Phi(),V0->Eta(),0.5};
                                fHistPhiEta->Fill(phiEtaK0Data);
                                selectedV0Assoc-> Add(new AliV0ChParticle(V0->Eta(), V0->Phi(), V0->Pt(), 5,0,myTrackPos->GetID(),myTrackNeg->GetID(),V0->GetOnFlyStatus(),massK0,myTrackPos->Phi(),myTrackPos->Pt(),myTrackPos->Eta(),myTrackPos->Charge(),myTrackNeg->Phi(),myTrackNeg->Pt(),myTrackNeg->Eta(),myTrackNeg->Charge()));
                            }
                        }
                    }
                }
            }
        }
            
        if (Antilambda&&cutAntiLambdaPid){
            if (IsMyGoodV0RapidityLambda(V0)) {
                fHistAntiLambdaMassPtCut->Fill(massAntilambda,ptTrig,5.5);

                if(IsMyGoodLifeTimeAntiLambda(lPVx,lPVy,lPVz,V0)){
                    fHistAntiLambdaMassPtCut->Fill(massAntilambda,ptTrig,6.5);

                    if (IsMyGoodV0AngleLambda(V0, myPrimVertex)){
                        fHistAntiLambdaMassPtCut->Fill(massAntilambda,ptTrig,7.5);
                
                        if(TMath::Abs(massK0-0.497614)>fMassRejectCutLam){
                            fHistAntiLambdaMassPtCut->Fill(massAntilambda,ptTrig,8.5);

            
                            if(fAnalysisMC){
                                FillMC(V0,mcArray,-3122,211, -2212,3,massAntilambda,selectedMCV0Triggersrec,fHistRecV0,fHistAntiLambdaMassPtCut,lPVz,myTrackPos,myTrackNeg,V0->GetOnFlyStatus(),fHistPurityCheck,selectedMCV0assoc,fHistPtResolution);
                                Double_t phiEtaAlamMC[4] = {V0->Pt(),V0->Phi(),V0->Eta(),2.5};
                                fHistPhiEta->Fill(phiEtaAlamMC);
                                
                            }
                            if(!fAnalysisMC) {
                                if(V0->Pt()>fPtTrigMin) selectedV0Triggers-> Add(new AliV0ChParticle(V0->Eta(), V0->Phi(), V0->Pt(), 3,0,myTrackPos->GetID(),myTrackNeg->GetID(),V0->GetOnFlyStatus(),massAntilambda,myTrackPos->Phi(),myTrackPos->Pt(),myTrackPos->Eta(),myTrackPos->Charge(),myTrackNeg->Phi(),myTrackNeg->Pt(),myTrackNeg->Eta(),myTrackNeg->Charge()));
                                Double_t phiEtaAlamData[4] = {V0->Pt(),V0->Phi(),V0->Eta(),2.5};
                                fHistPhiEta->Fill(phiEtaAlamData);
                                selectedV0Assoc-> Add(new AliV0ChParticle(V0->Eta(), V0->Phi(), V0->Pt(), 7,0,myTrackPos->GetID(),myTrackNeg->GetID(),V0->GetOnFlyStatus(),massAntilambda,myTrackPos->Phi(),myTrackPos->Pt(),myTrackPos->Eta(),myTrackPos->Charge(),myTrackNeg->Phi(),myTrackNeg->Pt(),myTrackNeg->Eta(),myTrackNeg->Charge()));
                            }
                        }
                    }
                }
            }
        }
    }

    fHistV0MultiplicityK0->Fill(nK0);
    fHistV0Lam->Fill(nLam);
    fHistMultiplicityALam->Fill(nALam);
	// Corelation ==========================================
    
     if(fAnalysisMC&&fCorrelations){
        //V0-h MC rec
        if(fV0hCorr) Corelations(selectedMCV0Triggersrec,selectedMCassoc,fHistKorelacieMCrec, lPVz, fHistNumberOfTriggersRec,kFALSE,kTRUE,lPercentile,fHistPtHard,ptHard,kFALSE,kFALSE);

        //h-h MC rec
        if(fhhCorr) Corelations(selectedMCtrig,selectedMCassoc,fHistKorelacieMCrec, lPVz, fHistNumberOfTriggersRec,kTRUE,kFALSE,lPercentile,fHistPtHard,ptHard,kFALSE,kFALSE);
         
         //MC rec h-V0
        if(fhV0Corr) Corelations(selectedMCtrig,selectedMCV0assoc,fHistKorelacieMCrec,lPVz,fHistNumberOfTriggersRec,kFALSE,kTRUE,lPercentile,fHistPtHard,ptHard,kTRUE,kFALSE);
    } else if(fCorrelations){
        //Data V0-h
        if(fV0hCorr) Corelations(selectedV0Triggers,selectedAssociatedTracks,fHistKorelacie,lPVz,fHistNumberOfTriggers,kFALSE,kTRUE,lPercentile,fHistPtHard,ptHard,kFALSE,kFALSE);

	    //Data h-h
        if(fhhCorr) Corelations(selectedTriggerTracks,selectedAssociatedTracks,fHistKorelacie,lPVz,fHistNumberOfTriggers,kFALSE,kFALSE,lPercentile,fHistPtHard,ptHard,kFALSE,kFALSE);
        
        //Data h-V0
        if(fhV0Corr) Corelations(selectedTriggerTracks,selectedV0Assoc,fHistKorelacie,lPVz,fHistNumberOfTriggers,kFALSE,kTRUE,lPercentile,fHistPtHard,ptHard,kTRUE,kFALSE);
    }


 	// Mixing ==============================================

    fHistMultVtxz->Fill(lPercentile,lPVz);

	fPool = fPoolMgr->GetEventPool(lPercentile, lPVz);
    if (!fPool) {
        AliWarning(Form("No pool found for centrality = %f, zVtx = %f", lPercentile, lPVz));
        return;
    }
 	Int_t nMix = fPool->GetCurrentNEvents();
	if (fPool->IsReady() || fPool->NTracksInPool() > fMixingTracks / 5 || nMix >= fMixedEvents)
	{
        for (Int_t jMix=0; jMix<nMix; jMix++)
 		{// loop through mixing events
 			TObjArray* bgTracks = fPool->GetEvent(jMix);
            if(fAnalysisMC&&fMixing) {
                
                if(fV0hCorr) CorelationsMixing(selectedMCV0Triggersrec,bgTracks,fHistMCMixingRec,lPVz,lPercentile);
                if(fhhCorr) CorelationsMixing(selectedMCtrig,bgTracks,fHistMCMixingRec,lPVz,lPercentile);
                if(fhV0Corr) CorelationsMixingV0h(bgTracks,selectedMCV0assoc,fHistMCMixingRec,lPVz,lPercentile);
            }else if(fMixing){
                if(fV0hCorr) CorelationsMixing(selectedV0Triggers,bgTracks,fHistdPhidEtaMix,lPVz,lPercentile);
                if(fhhCorr) CorelationsMixing(selectedTriggerTracks,bgTracks,fHistdPhidEtaMix,lPVz,lPercentile);
                if(fhV0Corr) CorelationsMixingV0h(bgTracks,selectedV0Assoc,fHistdPhidEtaMix,lPVz,lPercentile);
            }
		 }
	}
	TObjArray* cloneArray = (TObjArray *)selectedTracks->Clone();
	cloneArray->SetOwner(kTRUE);
	fPool->UpdatePool(cloneArray);
    
    // deleting TObjArrays
    mcTracksSel->Clear();
    delete mcTracksSel;
    mcTracksTrigSel->Clear();
    delete mcTracksTrigSel;
    mcTracksV0Sel->Clear();
    delete mcTracksV0Sel;
    selectedMCassoc->Clear();
    delete selectedMCassoc;
    selectedMCtrig->Clear();
    delete selectedMCtrig;
    selectedMCV0Triggersrec->Clear();
    delete selectedMCV0Triggersrec;
    mcArray->Clear("C");
    delete selectedAssociatedTracks;
    selectedTriggerTracks->Clear();
    delete selectedTriggerTracks;
    selectedV0->Clear();
    delete selectedV0;
    selectedV0Triggers->Clear();
    delete selectedV0Triggers;
    selectedV0Assoc->Clear();
    delete selectedV0Assoc;
    mcV0AssocSel->Clear();
    delete mcV0AssocSel;
    selectedMCV0assoc->Clear();
    delete selectedMCV0assoc;
    
    PostData(1, fOutputList);                           // stream the results the analysis of this event to
                                                        // the output manager which will take care of writing
                                                        // it to a file
}
//_____________________________________________________________________________
void AliAnalysisTaskDiHadCorrelHighPt::Terminate(Option_t *)
{
    // terminate
    // called at the END of the analysis (when all events are processed)
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskDiHadCorrelHighPt::IsMyGoodPrimaryTrack(const AliAODTrack *t)
 {
          // Pseudorapidity cut
          if (TMath::Abs(t->Eta())>=fEtaCut) return kFALSE;
		  if (!t->TestFilterBit(fFilterBit)) return kFALSE;
  
          return kTRUE;
 }
//____________________________________________________________________________
Bool_t AliAnalysisTaskDiHadCorrelHighPt::IsMyGoodV0AngleK0(const AliAODv0 *t, AliAODVertex *pv)
{
		//V0 Cosine of Pointing Angle
	    if(t->CosPointingAngle(pv)<=fCosPointAngleK0) return kFALSE;
    
		return kTRUE;
}
//____________________________________________________________________________
Bool_t AliAnalysisTaskDiHadCorrelHighPt::IsMyGoodV0AngleLambda(const AliAODv0 *t, AliAODVertex *pv)
{
		//V0 Cosine of Pointing Angle
		if(t->CosPointingAngle(pv)<=fCosPointAngleLam) return kFALSE;
		
		return kTRUE;
}
//____________________________________________________________________________
Bool_t AliAnalysisTaskDiHadCorrelHighPt::IsMyGoodV0RapidityLambda(const AliAODv0 *t)
{
		//Rapidity
		if(TMath::Abs(t->RapLambda())>=fRapidityCut) return kFALSE;
    //Pseudorap
   // if(TMath::Abs(t->Eta())>=0.8) return kFALSE;
		
		return kTRUE;
}
//____________________________________________________________________________
Bool_t AliAnalysisTaskDiHadCorrelHighPt::IsMyGoodV0RapidityK0(const AliAODv0 *t)
{
		//Rapidity
		if(TMath::Abs(t->RapK0Short())>=fRapidityCut) return kFALSE;
    //Pseudorap
    //if(TMath::Abs(t->Eta())>=0.8) return kFALSE;
		
		return kTRUE;
}

//____________________________________________________________________________
Bool_t AliAnalysisTaskDiHadCorrelHighPt::IsMyGoodLifeTimeK0(Double_t x,Double_t y, Double_t z, const AliAODv0 *V0)
{
		//Proper Lifetime (mL/p)
		
		Double_t *tParentVertexPosition = new Double_t[3];
		tParentVertexPosition[0]= x;
		tParentVertexPosition[1]= y;
		tParentVertexPosition[2]= z;

		Double_t lenght = V0->DecayLengthV0(tParentVertexPosition);
		Double_t momentum = TMath::Sqrt(TMath::Power(V0->Px(),2)+TMath::Power(V0->Py(),2)+TMath::Power(V0->Pz(),2));
		Double_t mass = V0->MassK0Short();
		Double_t lifetime = (mass*lenght)/momentum;
		
        delete[] tParentVertexPosition;
		if(lifetime>=fLifeTimeK0) return kFALSE;
		
		return kTRUE;
}
//____________________________________________________________________________
Bool_t AliAnalysisTaskDiHadCorrelHighPt::IsMyGoodLifeTimeLambda(Double_t x,Double_t y, Double_t z, const AliAODv0 *V0)
{
		//Proper Lifetime (mL/p)
		
		Double_t *tParentVertexPosition = new Double_t[3];
		tParentVertexPosition[0]= x;
		tParentVertexPosition[1]= y;
		tParentVertexPosition[2]= z;

		Double_t lenght = V0->DecayLengthV0(tParentVertexPosition);
		Double_t momentum = TMath::Sqrt(TMath::Power(V0->Px(),2)+TMath::Power(V0->Py(),2)+TMath::Power(V0->Pz(),2));
		Double_t mass = V0->MassLambda();
		Double_t lifetime = (mass*lenght)/momentum;
		
        delete[] tParentVertexPosition;
		if(lifetime>=fLifeTimeLam) return kFALSE;
		
		return kTRUE;
}
//____________________________________________________________________________
Bool_t AliAnalysisTaskDiHadCorrelHighPt::IsMyGoodLifeTimeAntiLambda(Double_t x,Double_t y, Double_t z, const AliAODv0 *V0)
{
    //Proper Lifetime (mL/p)
    
    Double_t *tParentVertexPosition = new Double_t[3];
    tParentVertexPosition[0]= x;
    tParentVertexPosition[1]= y;
    tParentVertexPosition[2]= z;
    
    Double_t lenght = V0->DecayLengthV0(tParentVertexPosition);
    Double_t momentum = TMath::Sqrt(TMath::Power(V0->Px(),2)+TMath::Power(V0->Py(),2)+TMath::Power(V0->Pz(),2));
    Double_t mass = V0->MassAntiLambda();
    Double_t lifetime = (mass*lenght)/momentum;
    
    delete[] tParentVertexPosition;
    if(lifetime>=fLifeTimeLam) return kFALSE;
    
    return kTRUE;
}

//____________________________________________________________________________
Bool_t AliAnalysisTaskDiHadCorrelHighPt::IsMyGoodPID(const AliAODTrack *TrackPos, const AliAODTrack *TrackNeg) {
	
	//TPC dE/dx selection
	Float_t dedxnegative = fPIDResponse->NumberOfSigmasTPC(TrackPos, AliPID::kPion);
	Float_t dedxpositive = fPIDResponse->NumberOfSigmasTPC(TrackNeg, AliPID::kPion);

	if(dedxnegative>=5) return kFALSE;
	if(dedxpositive>=5) return kFALSE;

	return kTRUE;

}
//____________________________________________________________________________
Bool_t AliAnalysisTaskDiHadCorrelHighPt::IsMyGoodDaughterTrack(const AliAODTrack *t) {
	// TPC refit
 	if (!t->IsOn(AliAODTrack::kTPCrefit)) return kFALSE;
	Float_t nCrossedRowsTPC = t->GetTPCClusterInfo(2,1);
	if (nCrossedRowsTPC < fnumOfTPCcrossedRows) return kFALSE;
	Int_t findable=t->GetTPCNclsF();
	if (findable <= 0) return kFALSE;
	if (nCrossedRowsTPC/findable < 0.8) return kFALSE;
    
    if (TMath::Abs(t->Eta())>=fEtaCut) return kFALSE;
		
	return kTRUE;
	
}
//____________________________________________________________________________
Bool_t AliAnalysisTaskDiHadCorrelHighPt::IsMyGoodV0(const AliAODv0 *v0,const AliAODTrack* myTrackPos, const AliAODTrack* myTrackNeg, Int_t oSta) {
	if (!v0) {
         AliError(Form("ERROR: Could not retrieve aodV0"));
         return kFALSE;
	}

	if (oSta==1) {if (v0->GetOnFlyStatus()) return kFALSE;} //offline
	if (oSta==3) {if (!v0->GetOnFlyStatus()) return kFALSE;} // "on fly" during the tracking
    
    if (oSta==2){         
	    if (v0->GetOnFlyStatus()){
			return kTRUE;
 	    } else {
			return kFALSE;
	    }
	}
    if (oSta==4){
	    if (!v0->GetOnFlyStatus()){
 			return kTRUE;
        } else {
 			return kFALSE;
 		}
	}
 	// Track cuts for daughter tracks
   	if ( !(IsMyGoodDaughterTrack(myTrackPos)) || !(IsMyGoodDaughterTrack(myTrackNeg)) ) return kFALSE;
    return kTRUE;
}
//______________________________________________________________________
Bool_t AliAnalysisTaskDiHadCorrelHighPt::IsMyGoodV0Topology(const AliAODv0 *v0){
    if (!v0) {
        AliError(Form("ERROR: Could not retrieve aodV0"));
        return kFALSE;
    }
	//DCA Negative Track to PV
	if(v0->DcaNegToPrimVertex()<=fDCAposDaughter) return kFALSE;
	//DCA Positive Track to PV
	if(v0->DcaPosToPrimVertex()<=fDCAnegDaughter) return kFALSE;
	//DCA V0 daughters
	if(v0->DcaV0Daughters()>=fDCAV0Daughters) return kFALSE;
	//V0 2D Decay Radius
	if(v0->RadiusV0()<=fV0Radius) return kFALSE;
	
	return kTRUE;
}
//____________________________________________________________________________
void AliAnalysisTaskDiHadCorrelHighPt::Corelations(TObjArray *triggers, TObjArray *associated, THnSparse * fHistKor, Double_t lPVz, THnSparse* fHistNumOfTrig,Bool_t hh,Bool_t V0h,Float_t perc,TH3F *fHistPtHard, Double_t ptHard,Bool_t hV0,Bool_t isMCGen){

    const Float_t kLimit = fMergingCut * 3;
    Float_t dphistarminabs = 1e5;
    const Double_t kPi = TMath::Pi();
    Int_t nAssoc = associated->GetEntriesFast();
    Int_t nTrig = triggers->GetEntriesFast();
    Double_t status = 0.;
    Double_t detaDau =0;
    Double_t assocPt =0.;
    Double_t triggPt =0.;
    Double_t asocEta =0;
    Double_t assocPhi =0;
    Double_t triggEta =0.;
    Double_t ptDaug =0.;
    Double_t assocCharge =0.;
    Double_t phiDaug =0.;
    Double_t charDaug =0.;

    for (Int_t i=0; i<nTrig; i++){
        AliV0ChParticle* trig = (AliV0ChParticle*)  triggers->At(i);
        triggPt = trig->Pt();
        triggEta = trig->Eta();
        if(triggPt<fPtTrigMin) continue;
        
        if (trig->GetRecStatus()) status=0.5;
        if (!trig->GetRecStatus()) status=1.5;
        if (ptHard!=0) fHistPtHard->Fill(triggPt/ptHard,perc,trig->WhichCandidate()-0.5);
        Double_t massTrig = 0.;
        if(trig->WhichCandidate()<4) massTrig=trig->M();
        
        if(!hV0){
            Double_t triggers[6]={triggPt,lPVz,triggEta,trig->WhichCandidate()-0.5,massTrig,perc};
            fHistNumOfTrig->Fill(triggers);
        }else{
            Double_t triggers[6]={triggPt,lPVz,triggEta,4.5,massTrig,perc};
            fHistNumOfTrig->Fill(triggers);
        }
        
        for (Int_t j=0; j<nAssoc; j++){
            AliV0ChParticle* assoc = (AliV0ChParticle*)  associated->At(j);
            if(assoc->WhichCandidate() < 4) continue;
            asocEta = assoc->Eta();
            assocPhi = assoc->Phi();
            assocCharge = assoc->Charge();

            Double_t deltaEta = triggEta - asocEta;
            Double_t deltaPhi = trig->Phi() - assocPhi;
            assocPt = assoc->Pt();
            if (deltaPhi > (1.5*kPi)) deltaPhi -= 2.0*kPi;
            if (deltaPhi < (-0.5*kPi)) deltaPhi += 2.0*kPi;

            if(triggPt<=assocPt) continue;
            
            //removing autocorrelations
            if(V0h){
                
                Int_t negID = 0;
                Int_t posID = 0;
                Int_t atrID = 0;
                
                if(!hV0){
                    negID = trig->GetIDNeg();
                    posID = trig->GetIDPos();
                    atrID = assoc->GetIDCh();
                }
                else{
                    negID = assoc->GetIDNeg();
                    posID = assoc->GetIDPos();
                    atrID = trig->GetIDCh();
                    
                }
                
                if ((TMath::Abs(negID))==(TMath::Abs(atrID))) continue;
                if ((TMath::Abs(posID))==(TMath::Abs(atrID))) continue;
                
            }
            
            Double_t massK0=10;
            Double_t massLam=10;
            Double_t massGamma=10;
            Double_t massSigmaP=10;
            Double_t massSigmaN=10;
            Double_t massXiN=10;
            Double_t massOmegaN=10;
            
            if(hh&&fRemoveHadrFromV0){
                if(trig->Charge()!=assocCharge){
                    massK0 = TMath::Sqrt(2*0.13957*0.13957+2*(trig->E()*assoc->E()-triggPt*assocPt-trig->Pz()*assoc->Pz()));
                    massLam = TMath::Sqrt(0.13957*0.13957+0.93827*0.93827+2*(trig->E()*assoc->E()-triggPt*assocPt-trig->Pz()*assoc->Pz()));
                    massGamma = TMath::Sqrt(2*0.0005109*0.0005109+2*(trig->E()*assoc->E()-triggPt*assocPt-trig->Pz()*assoc->Pz()));
                }
            }
            if(TMath::Abs( 0.497614-massK0)< 0.005) continue;
            if(TMath::Abs( 1.115683-massLam)< 0.005) continue;
            if(TMath::Abs(massGamma)<0.004) continue;
            
            if(fRemoveLamhFromCascade&&(trig->WhichCandidate()==2||trig->WhichCandidate()==3)){
                massSigmaP=TMath::Sqrt(0.13957*0.13957+1.1156*1.1156+2*(trig->E()*assoc->E()-triggPt*assocPt-trig->Pz()*assoc->Pz()));
                massSigmaN=TMath::Sqrt(0.13957*0.13957+1.1156*1.1156+2*(trig->E()*assoc->E()-triggPt*assocPt-trig->Pz()*assoc->Pz()));
                massXiN=TMath::Sqrt(0.13957*0.13957+1.1156*1.1156+2*(trig->E()*assoc->E()-triggPt*assocPt-trig->Pz()*assoc->Pz()));
                massOmegaN=TMath::Sqrt(0.4936*0.4936+1.1156*1.1156+2*(trig->E()*assoc->E()-triggPt*assocPt-trig->Pz()*assoc->Pz()));
            }
            
            if(TMath::Abs( 1.3872-massSigmaN)< 0.005||TMath::Abs( 1.3828-massSigmaP)<0.005||TMath::Abs( 1.32171-massXiN)<0.005||TMath::Abs( 1.67245-massOmegaN)<0.005) continue;

            if(!isMCGen&&trig->WhichCandidate()==2&&TMath::Abs(deltaPhi)<0.7&&TMath::Abs(deltaEta)<0.7&&triggPt<=4&&assocPt<=2) fHistDeltaEtaDeltaPhiLamFineBinningLowPt->Fill(deltaPhi,deltaEta); // check how does it look like before track merging cut
            
            if(!hV0) {
                // track merging cut (hadron with V0 daughter tracks)
                if(trig->WhichCandidate()<4){
                    if(assocCharge*trig->GetCharge1()>0) {
                        phiDaug = trig->GetPhi1();
                        ptDaug = trig->GetPt1();
                        charDaug = trig->GetCharge1();
                        detaDau = trig->GetEta1()-asocEta;
                        if(TMath::Abs(detaDau)<fMergingCut*2.5*3){ //default fMergingCut=0.02
                            Float_t dphistarLow = GetDPhiStar(phiDaug, ptDaug, charDaug, assocPhi, assocPt, assocCharge, 0.8);
                            Float_t dphistarHigh = GetDPhiStar(phiDaug, ptDaug, charDaug, assocPhi, assocPt, assocCharge, 2.5);
                            if (TMath::Abs(dphistarLow) < kLimit || TMath::Abs(dphistarHigh) < kLimit || dphistarLow * dphistarHigh < 0 ) { //searching for dphistar minimum
                                dphistarminabs = 1e5;
                                for (Double_t rad=0.8; rad<2.51; rad+=0.01) {
                                    Float_t dphistar = GetDPhiStar(phiDaug, ptDaug, charDaug, assocPhi, assocPt, assocCharge, rad);
                
                                    Float_t dphistarabs = TMath::Abs(dphistar);
        
                                    if (dphistarabs < dphistarminabs) {
                                        dphistarminabs = dphistarabs;
                                    }
                                }
          
                                if (dphistarminabs < fMergingCut && TMath::Abs(detaDau) < fMergingCut) continue;
                            } 
                        }
                    }
                    if(assocCharge*trig->GetCharge2()>0){
                        phiDaug = trig->GetPhi2();
                        ptDaug = trig->GetPt2();
                        charDaug = trig->GetCharge2();
                        detaDau = trig->GetEta2()-asocEta;
                        if(TMath::Abs(detaDau)<fMergingCut*2.5*3){ //default fMergingCut=0.02
                            Float_t dphistarLow = GetDPhiStar(phiDaug, ptDaug, charDaug, assocPhi, assocPt, assocCharge, 0.8);
                            Float_t dphistarHigh = GetDPhiStar(phiDaug, ptDaug, charDaug, assocPhi, assocPt, assocCharge, 2.5);
                            if (TMath::Abs(dphistarLow) < kLimit || TMath::Abs(dphistarHigh) < kLimit || dphistarLow * dphistarHigh < 0 ) { //searching for dphistar minimum
                                dphistarminabs = 1e5;
                                for (Double_t rad=0.8; rad<2.51; rad+=0.01) {
                                    Float_t dphistar = GetDPhiStar(phiDaug, ptDaug, charDaug, assocPhi, assocPt, assocCharge, rad);
                
                                    Float_t dphistarabs = TMath::Abs(dphistar);
        
                                    if (dphistarabs < dphistarminabs) {
                                        dphistarminabs = dphistarabs;
                                    }
                                }
          
                                if (dphistarminabs < fMergingCut && TMath::Abs(detaDau) < fMergingCut) continue;
                            } 
                        }
                    }
                }
                Double_t korel[10] = {triggPt,assocPt,deltaPhi,deltaEta, lPVz,trig->WhichCandidate()-0.5, triggEta,asocEta,massTrig,perc};
                fHistKor->Fill(korel);
            }else{
                massTrig = assoc->M();
                Double_t korel[10] = {triggPt,assocPt,deltaPhi,deltaEta, lPVz,assoc->WhichCandidate()-0.5, triggEta,asocEta,massTrig,perc};
                fHistKor->Fill(korel);
            }
            
        }
    }
    

}

//____________________________________________________________________________
void AliAnalysisTaskDiHadCorrelHighPt::CorelationsMixing(TObjArray *triggers, TObjArray *bgTracks, THnSparse * fHistKor, Double_t lPVz, Float_t perc){

    const Double_t kPi = TMath::Pi();
    Int_t nAssoc = bgTracks->GetEntriesFast();
    Int_t nTrig = triggers->GetEntriesFast();
    Double_t assocCharge =0.;
    Double_t phiDaug =0.;
    Double_t ptDaug =0.;
    Double_t detaDau =0.;
    Double_t asocEta =0.;
    Double_t assocPhi =0.;
    Double_t assocPt =0.;
    Double_t charDaug =0.;
    const Float_t kLimit = fMergingCut * 3;
    Float_t dphistarminabs = 1e5;

    for (Int_t i=0; i<nTrig; i++){
        AliV0ChParticle* trig = (AliV0ChParticle*)  triggers->At(i);
        
        Double_t massTrig = 0.;
        if(trig->WhichCandidate()<4) massTrig=trig->M();
        for (Int_t j=0; j<nAssoc; j++){

             AliAODTrack* assoc = (AliAODTrack*) bgTracks->At(j);
             assocCharge = assoc->Charge();
             asocEta = assoc->Eta();
             assocPhi = assoc -> Phi();
             assocPt = assoc->Pt();

             if (( assocPt>=trig->Pt() ) || ( assocPt<fPtAsocMin )) continue;

             if(!IsMyGoodPrimaryTrack(assoc)) continue; 

             Double_t   deltaEta = trig->Eta() - asocEta;
             Double_t   deltaPhi = trig->Phi() - assocPhi;
             

            if (deltaPhi > (1.5*kPi)) deltaPhi -= 2.0*kPi;
            if (deltaPhi < (-0.5*kPi)) deltaPhi += 2.0*kPi;

            // track merging cut (hadron with V0 daughter tracks)
                if(trig->WhichCandidate()<4){
                    if(assocCharge*trig->GetCharge1()>0) {
                        phiDaug = trig->GetPhi1();
                        ptDaug = trig->GetPt1();
                        charDaug = trig->GetCharge1();
                        detaDau = trig->GetEta1()-asocEta;
                        if(TMath::Abs(detaDau)<fMergingCut*2.5*3){ //default fMergingCut=0.02
                            Float_t dphistarLow = GetDPhiStar(phiDaug, ptDaug, charDaug, assocPhi, assocPt, assocCharge, 0.8);
                            Float_t dphistarHigh = GetDPhiStar(phiDaug, ptDaug, charDaug, assocPhi, assocPt, assocCharge, 2.5);
                            if (TMath::Abs(dphistarLow) < kLimit || TMath::Abs(dphistarHigh) < kLimit || dphistarLow * dphistarHigh < 0 ) { //searching for dphistar minimum
                                dphistarminabs = 1e5;
                                for (Double_t rad=0.8; rad<2.51; rad+=0.01) {
                                    Float_t dphistar = GetDPhiStar(phiDaug, ptDaug, charDaug, assocPhi, assocPt, assocCharge, rad);
                
                                    Float_t dphistarabs = TMath::Abs(dphistar);
        
                                    if (dphistarabs < dphistarminabs) {
                                        dphistarminabs = dphistarabs;
                                    }
                                }
          
                                if (dphistarminabs < fMergingCut && TMath::Abs(detaDau) < fMergingCut) continue;
                            } 
                        }
                    }
                    if(assocCharge*trig->GetCharge2()>0){
                        phiDaug = trig->GetPhi2();
                        ptDaug = trig->GetPt2();
                        charDaug = trig->GetCharge2();
                        detaDau = trig->GetEta2()-asocEta;
                        if(TMath::Abs(detaDau)<fMergingCut*2.5*3){ //default fMergingCut=0.02
                            Float_t dphistarLow = GetDPhiStar(phiDaug, ptDaug, charDaug, assocPhi, assocPt, assocCharge, 0.8);
                            Float_t dphistarHigh = GetDPhiStar(phiDaug, ptDaug, charDaug, assocPhi, assocPt, assocCharge, 2.5);
                            if (TMath::Abs(dphistarLow) < kLimit || TMath::Abs(dphistarHigh) < kLimit || dphistarLow * dphistarHigh < 0 ) { //searching for dphistar minimum
                                dphistarminabs = 1e5;
                                for (Double_t rad=0.8; rad<2.51; rad+=0.01) {
                                    Float_t dphistar = GetDPhiStar(phiDaug, ptDaug, charDaug, assocPhi, assocPt, assocCharge, rad);
                
                                    Float_t dphistarabs = TMath::Abs(dphistar);
        
                                    if (dphistarabs < dphistarminabs) {
                                        dphistarminabs = dphistarabs;
                                    }
                                }
          
                                if (dphistarminabs < fMergingCut && TMath::Abs(detaDau) < fMergingCut) continue;
                            } 
                        }
                    }
                }
            Double_t korel[10] = {trig->Pt(),assocPt,deltaPhi,deltaEta, lPVz,trig->WhichCandidate()-0.5,trig->Eta(),asocEta,massTrig,perc};
            fHistKor->Fill(korel);
        }
    }
    

}
//_____________________________________________________________
void AliAnalysisTaskDiHadCorrelHighPt::CorelationsMixingV0h(TObjArray *bgTracks, TObjArray *assocArray, THnSparse * fHistKor, Double_t lPVz, Float_t perc){
    
    const Double_t kPi = TMath::Pi();
    Int_t nAssoc = assocArray->GetEntriesFast();
    Int_t nTrig = bgTracks->GetEntriesFast();
    
    for (Int_t i=0; i<nTrig; i++){
        AliAODTrack* trig = (AliAODTrack*)  bgTracks->At(i);
        if(trig->Pt()<fPtTrigMin) continue;
        if(!IsMyGoodPrimaryTrack(trig)) continue;
        
        for (Int_t j=0; j<nAssoc; j++){
            AliV0ChParticle* assoc = (AliV0ChParticle*) assocArray->At(j);
            
            Double_t massAssoc = assoc->M();
            if(assoc->WhichCandidate()==5&&(massAssoc<0.486||massAssoc>0.509)) continue;
            if((assoc->WhichCandidate()==6||assoc->WhichCandidate()==7)&&(massAssoc<1.1112||massAssoc>1.12)) continue;
            
            if (( (assoc->Pt())>=trig->Pt() ) || ( (assoc->Pt())<fPtAsocMin )) continue;
            
            Double_t   deltaEta = trig->Eta() - assoc->Eta();
            Double_t   deltaPhi = trig->Phi() - assoc->Phi();
            Double_t   assocPt = assoc->Pt();
            
            if (deltaPhi > (1.5*kPi)) deltaPhi -= 2.0*kPi;
            if (deltaPhi < (-0.5*kPi)) deltaPhi += 2.0*kPi;
            
            Double_t korel[7] = {trig->Pt(),assocPt,deltaPhi,deltaEta, lPVz,assoc->WhichCandidate()-0.5,perc};
            fHistKor->Fill(korel);
        }
    }
}
//_____________________________________________________________
void AliAnalysisTaskDiHadCorrelHighPt::TopologCuts(THnSparse* fHist,Double_t pttrig,Double_t mass,Double_t dcaNeg, Double_t dcaPos,Double_t dcaDau, Double_t V0rad, Double_t cosPA,Double_t lifetime,Double_t massSell,Double_t triggType,Double_t status){
    
    Double_t topolCutsValues[11]={pttrig,mass,dcaNeg,dcaPos,dcaDau,V0rad,cosPA,lifetime,massSell,triggType,status};
    fHist->Fill(topolCutsValues);
    
}
//____________________________________________________________
void AliAnalysisTaskDiHadCorrelHighPt::FillMC(const AliAODv0 *V0,TClonesArray *mcArray,Int_t pdgV0,Int_t pdgDau1, Int_t pdgDau2,Int_t triggerType, Double_t mass, TObjArray * selectedMCV0Triggersrec,THnSparse * fHistRecV0, TH3F * fHistMassPtCut,Double_t lPVz, const AliAODTrack * myTrackPos,const AliAODTrack * myTrackNeg,Bool_t status,THnSparse * histPur, TObjArray * selectedMCV0assoc,TH3F * fHistresol){
    
    if(fPurityCheck){
        Double_t purity[6] ={V0->Pt(),mass,triggerType-0.5,0.5,-1,V0->Eta()};
        histPur->Fill(purity);
    }
    
    if(fCorrelations&&!fPureV0) { // for MC closure test - also misidentified V0 taken
       if(V0->Pt()>fPtTrigMin) selectedMCV0Triggersrec-> Add(new AliV0ChParticle(V0->Eta(), V0->Phi(), V0->Pt(), triggerType,0,myTrackPos->GetID(),myTrackNeg->GetID(),status,mass,myTrackPos->Phi(),myTrackPos->Pt(),myTrackPos->Eta(),myTrackPos->Charge(),myTrackNeg->Phi(),myTrackNeg->Pt(),myTrackNeg->Eta(),myTrackNeg->Charge())); // all reconstructed candidates for raw correlation function, with reconstructed pt
       if(V0->Pt()>fPtAsocMin) selectedMCV0assoc -> Add(new AliV0ChParticle(V0->Eta(), V0->Phi(), V0->Pt(), triggerType,0,myTrackPos->GetID(),myTrackNeg->GetID(),status,mass,myTrackPos->Phi(),myTrackPos->Pt(),myTrackPos->Eta(),myTrackPos->Charge(),myTrackNeg->Phi(),myTrackNeg->Pt(),myTrackNeg->Eta(),myTrackNeg->Charge())); // all reconstructed candidates for raw correlation function, with reconstructed pt
    }
    
    Int_t myTrackPosLabel = TMath::Abs(myTrackPos->GetLabel());
    Int_t myTrackNegLabel = TMath::Abs(myTrackNeg->GetLabel());
    
    AliAODMCParticle *mcPosTrack = (AliAODMCParticle*)mcArray->At(myTrackPosLabel);
    if (!mcPosTrack) return;
    Int_t PosTrackPdg = mcPosTrack->GetPdgCode();
    AliAODMCParticle *mcNegTrack = (AliAODMCParticle*)mcArray->At(myTrackNegLabel);
    if (!mcNegTrack) return;
    Int_t NegTrackPdg = mcNegTrack->GetPdgCode();
    
    if(fPurityCheck){
        Double_t puri[6] ={V0->Pt(),mass,triggerType-0.5,1.5,-1,V0->Eta()};
        histPur->Fill(puri);
    }
    
    Int_t myTrackPosMotherLabel = mcPosTrack->GetMother();
    Int_t myTrackNegMotherLabel = mcNegTrack->GetMother();
    
    if ((myTrackPosMotherLabel==-1)||(myTrackNegMotherLabel==-1)) return;
    
    if(fPurityCheck){
        Double_t pur[6] ={V0->Pt(),mass,triggerType-0.5,2.5,-1,V0->Eta()};
        histPur->Fill(pur);
    }
    
    if (myTrackPosMotherLabel!=myTrackNegMotherLabel) return;
    if(fPurityCheck){
        Double_t pu[6] ={V0->Pt(),mass,triggerType-0.5,3.5,-1,V0->Eta()};
        histPur->Fill(pu);
    }
    
    AliAODMCParticle *mcPosMother = (AliAODMCParticle*)mcArray->At(myTrackPosMotherLabel);
    if (!mcPosMother) return;

    Int_t MotherPdg = mcPosMother->GetPdgCode();
    Int_t MotherOfMotherLabel = mcPosMother->GetMother();
    
    Bool_t IsPhysPrim = mcPosMother->IsPhysicalPrimary();
    
    Bool_t IsFromMC = (MotherPdg==pdgV0)&&(PosTrackPdg==pdgDau1)&&(NegTrackPdg==pdgDau2);

    Bool_t IsFromCascade = kFALSE;
    
    Int_t MoMPdg = 50;
    
    if (MotherOfMotherLabel != -1)
    {
        AliAODMCParticle *mcPosMotherOfMother = (AliAODMCParticle*)mcArray->At(MotherOfMotherLabel);
        Int_t MotherOfMotherPdg = mcPosMotherOfMother->GetPdgCode();
        MoMPdg = TMath::Abs(MotherOfMotherPdg);
        if(fAacceptLambdasFromCasscade) IsFromCascade = (((MoMPdg == 3222)|| (MoMPdg==3212)|| (MoMPdg==3112) || (MoMPdg==3224) || (MoMPdg==3214) || (MoMPdg==3114) || (MoMPdg==3322) || (MoMPdg==3312)|| (MoMPdg==3324) || (MoMPdg==3314) || (MoMPdg==3334)) && (mcPosMotherOfMother->IsPhysicalPrimary()));
    }
    if(fPurityCheck){
        Double_t purit[6] ={V0->Pt(),mass,triggerType-0.5,4.5,-1,V0->Eta()};
        histPur->Fill(purit);
    }
    Bool_t isGoodID = (MotherPdg==pdgV0);
    Bool_t isFromMaterial = mcPosMother->IsSecondaryFromMaterial();
     Bool_t isFromDecay = mcPosMother->IsSecondaryFromWeakDecay();

    if(!isGoodID&&fPurityCheck) {
        Int_t ident =0;
        if(MotherPdg==211) ident =1;
        else if(MotherPdg==3122) ident=2;
        else if(MotherPdg==-3122) ident=3;
        else if(MotherPdg==-211) ident =4;
        else if(MotherPdg==310) ident =5;
        else if(MotherPdg==22) ident =6;
        else if(MotherPdg==223) ident =7;
        else ident=8;
        Double_t purity[6] ={V0->Pt(),mass,triggerType-0.5,5.5,ident-0.5,V0->Eta()};
        histPur->Fill(purity);
    }
    
    if(isGoodID&&fPurityCheck){
        Double_t purity[6] ={V0->Pt(),mass,triggerType-0.5,6.5,-1,V0->Eta()};
        histPur->Fill(purity);
    }

    if(!isFromMaterial&&isGoodID&&fPurityCheck){
        Double_t purity[6] ={V0->Pt(),mass,triggerType-0.5,7.5,-1,V0->Eta()};
        histPur->Fill(purity);
    }
    
    Bool_t IsParticleFromMC = kFALSE;
    
    if(pdgV0==310) IsParticleFromMC= (IsFromMC&&IsPhysPrim);
    else {
        IsParticleFromMC = (IsFromMC&&(IsPhysPrim||IsFromCascade));
    }
    
    Double_t V0mcPt = mcPosMother->Pt();
   
    if(IsParticleFromMC){
        if(fCorrelations&&fPureV0) { // for MC closure test - only good ID V0 taken
            if(V0->Pt()>fPtTrigMin) selectedMCV0Triggersrec-> Add(new AliV0ChParticle(V0->Eta(), V0->Phi(), V0->Pt(), triggerType,0,myTrackPos->GetID(),myTrackNeg->GetID(),status,mass,myTrackPos->Phi(),myTrackPos->Pt(),myTrackPos->Eta(),myTrackPos->Charge(),myTrackNeg->Phi(),myTrackNeg->Pt(),myTrackNeg->Eta(),myTrackNeg->Charge())); // all reconstructed candidates for raw correlation function, with reconstructed pt
            if(V0->Pt()>fPtAsocMin) selectedMCV0assoc -> Add(new AliV0ChParticle(V0->Eta(), V0->Phi(), V0->Pt(), triggerType,0,myTrackPos->GetID(),myTrackNeg->GetID(),status,mass,myTrackPos->Phi(),myTrackPos->Pt(),myTrackPos->Eta(),myTrackPos->Charge(),myTrackNeg->Phi(),myTrackNeg->Pt(),myTrackNeg->Eta(),myTrackNeg->Charge())); // all reconstructed candidates for raw correlation function, with reconstructed pt
        }
        fHistresol->Fill(V0mcPt,V0->Pt(),triggerType-0.5);
        Double_t V0mcEta = mcPosMother->Eta();
    
        if(fEfficiency){
            Double_t v0effic[6]={V0mcPt,lPVz,triggerType-0.5,V0mcEta,mass};
            fHistRecV0->Fill(v0effic);
        }

        fHistMassPtCut->Fill(mass,V0mcPt,9.5);
    }

}
//____________________________________________________________________//
Float_t AliAnalysisTaskDiHadCorrelHighPt::GetDPhiStar(Float_t phi1, Float_t pt1, Float_t charge1, Float_t phi2, Float_t pt2, Float_t charge2, Float_t radius) { 
  //
  // calculates dphistar
  //
    Float_t bSign = (fAOD->GetMagneticField() > 0) ? 1 : -1;
    Float_t dphistar = phi1 - phi2 + charge1 * bSign * TMath::ASin(0.075 * radius / pt1) - charge2 * bSign * TMath::ASin(0.075 * radius / pt2);
  
    static const Double_t kPi = TMath::Pi();
  
    if (dphistar > kPi)
        dphistar = kPi * 2 - dphistar;
    if (dphistar < -kPi)
        dphistar = -kPi * 2 - dphistar;
  
  return dphistar;
}

