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
 * Last update edited by Lucia Anna Husova, December 2018
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
    fPtTrigMin(0),
    fPtAsocMin(0),
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
    fHistKorelPurCorr(0),
    fCosPointAngleK0(0.97),
    fCosPointAngleLam(0.995),
    fDCAV0Daughters(1),
    fDCAposDaughter(0.06),
    fDCAnegDaughter(0.06),
    fnumOfTPCcrossedRows(70),
    fEfficiency(kTRUE),
    fPurityCheck(kTRUE),
    fCorrelations(kTRUE),
    fHistKorelResolCorr(0),
    fHistKorelPurCorrGen(0),
    fHistPtResolution(0)
{
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
AliAnalysisTaskDiHadCorrelHighPt::AliAnalysisTaskDiHadCorrelHighPt(const char* name, Bool_t analysisMC) : AliAnalysisTaskSE(name),
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
    fPtTrigMin(0),
    fPtAsocMin(0),
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
    fHistKorelPurCorr(0),
    fCosPointAngleK0(0.97),
    fCosPointAngleLam(0.995),
    fDCAV0Daughters(1),
    fDCAposDaughter(0.06),
    fDCAnegDaughter(0.06),
    fnumOfTPCcrossedRows(70),
    fEfficiency(kTRUE),
    fPurityCheck(kTRUE),
    fCorrelations(kTRUE),
    fHistKorelResolCorr(0),
    fHistKorelPurCorrGen(0),
    fHistPtResolution(0)
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
	Float_t kPtBins[11] = {0,1,2,3,4,5,6,7,9,11,15};
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

    Int_t bins[10]= {12,14,72,75,20,4,40,40,601,10};
    Double_t min[10] = {fPtTrigMin,fPtAsocMin, -kPi/2, -2., -10., 0.,-0.8,-0.8,0.44,0};
    Double_t max[10] = {15., 15., -kPi/2+2*kPi, 2., 10., 4.,0.8,0.8, 1.15,100};
    
    Int_t binsMix[7] = {12,14,72,75,20,4,10};
    Double_t minMix[7] ={fPtTrigMin,fPtAsocMin, -kPi/2, -2., -10., 0.,0};
    Double_t maxMix[7] = {15., 15., -kPi/2+2*kPi, 2., 10., 4.,100};
    
	Int_t  NofCentBins  = 10;
    Double_t MBins[]={0,10,20,30,40,50,60,70,80,90,100};
         
    Int_t NofZVrtxBins  = 20;
    Double_t ZBins[]={-10.0, -9., -8.0, -7., -6.0, -5., -4.0, -3., -2., -1., 0, 1., 2.0 ,3., 4.0, 5., 6.0, 7., 8.0, 9., 10.0};

	Int_t bins2d[6] = {12,20,40,4,601,10};
	Double_t mis2d[6] = {fPtTrigMin,-10,-0.8,0.,0.44,0};
	Double_t maxs2d[6] = {15.,10,0.8,4.,1.15,100};
	
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
        
    fHistK0MassPtCut = new TH3F ("fHistK0MassPtCut", "fHistK0MassPtCut", kMassBins, kMassBinsK, 10, kPtBins, kNCuts, kCuts);
    fHistLambdaMassPtCut = new TH3F("fHistLambdaMassPtCut", "fHistLambdaMassPtCut", kMassBins, kMassBinsLambda, 10, kPtBins, kNCuts, kCuts);
	fHistAntiLambdaMassPtCut = new TH3F("fHistAntiLambdaMassPtCut", "fHistAntiLambdaMassPtCut", kMassBins, kMassBinsLambda, 10, kPtBins, kNCuts, kCuts);
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
    
	fHistdPhidEtaMix = new THnSparseF ("fHistdPhidEtaMix", "fHistdPhidEtaMix", 7, binsMix, minMix, maxMix);
    fHistdPhidEtaMix->GetAxis(0)->SetTitle("p_{T}^{trig}");
    fHistdPhidEtaMix->GetAxis(1)->SetTitle("p_{T}^{assoc}");
    fHistdPhidEtaMix->GetAxis(2)->SetTitle("#Delta#phi");
    fHistdPhidEtaMix->GetAxis(3)->SetTitle("#Delta#eta");
    fHistdPhidEtaMix->GetAxis(4)->SetTitle("p_{vz}");
    fHistdPhidEtaMix->GetAxis(5)->SetTitle("trigger");
    fHistdPhidEtaMix->GetAxis(6)->SetTitle("multiplicity percentile");
    fHistdPhidEtaMix->Sumw2();
	fOutputList->Add(fHistdPhidEtaMix);
    
    fHistMCMixingRec = new THnSparseF ("fHistMCMixingRec", "fHistMCMixingRec", 7, binsMix, minMix, maxMix);
    fHistMCMixingRec->GetAxis(0)->SetTitle("p_{T}^{trig}");
    fHistMCMixingRec->GetAxis(1)->SetTitle("p_{T}^{assoc}");
    fHistMCMixingRec->GetAxis(2)->SetTitle("#Delta#phi");
    fHistMCMixingRec->GetAxis(3)->SetTitle("#Delta#eta");
    fHistMCMixingRec->GetAxis(4)->SetTitle("p_{vz}");
    fHistMCMixingRec->GetAxis(5)->SetTitle("trigger");
    fHistMCMixingRec->GetAxis(6)->SetTitle("multiplicity percentile");
    fOutputList->Add(fHistMCMixingRec);
    fHistMCMixingRec->Sumw2();

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

	fHistV0Multiplicity = new TH1D ("fHistV0Multiplicity", "fHistV0Multiplicity", 60, 0, 60);
	fOutputList->Add(fHistV0Multiplicity); 

	fHistMultVtxz = new TH2D ("fHistMultVtxz","fHistMultVtxz",NofCentBins,MBins,NofZVrtxBins,ZBins);
	fOutputList->Add(fHistMultVtxz);

	fHistMCPtAs = new TH3D("fHistMCPtAs","fHistMCPtAs",14,fPtAsocMin,15,20,-10,10,40,-0.8,0.8);
	fOutputList->Add(fHistMCPtAs);
    fHistMCPtAs->Sumw2();
	fHistRCPtAs = new TH3D("fHistRCPtAs","fHistRCPtAs",14,fPtAsocMin,15,20,-10,10,40,-0.8,0.8);
	fOutputList->Add(fHistRCPtAs);
    fHistRCPtAs->Sumw2();
    
    fHistMCPtTrigg = new TH3D("fHistMCPtTrigg","fHistMCPtTrigg",12,fPtTrigMin,15,20,-10,10,40,-0.8,0.8);
    fOutputList->Add(fHistMCPtTrigg);
    fHistMCPtTrigg->Sumw2();
    fHistRCPtTrigg = new TH3D("fHistRCPtTrigg","fHistRCPtTrigg",12,fPtTrigMin,15,20,-10,10,40,-0.8,0.8);
    fOutputList->Add(fHistRCPtTrigg);
    fHistRCPtTrigg->Sumw2();
    
    Int_t binsTrig[4]={12,20,3,40};
    Double_t mintrig[4]={fPtTrigMin,-10,0,-0.8};
    Double_t maxtrig[4]={15,10,3,0.8};
    fHistGenV0 = new THnSparseF("fHistGenV0","fHistGenV0",4,binsTrig,mintrig,maxtrig);
    fOutputList->Add(fHistGenV0);
    fHistGenV0->Sumw2();
    Int_t binsTrigRec[5]={12,20,3,40,602};
    Double_t mintrigRec[6]={fPtTrigMin,-10,0,-0.8,0.44};
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

    fHistSelection = new TH1D("fHistSelection","fHistSelection",3,0,3);
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
    
    Int_t binsPur[5] = {12,602,3,8,8};
    Double_t binsPurMin[5] = {fPtTrigMin,0.44,0.,0.,0};
    Double_t binsPurMax[5] = {15,1.15,3.,8.,8.};
    fHistPurityCheck = new THnSparseF("fHistPurityCheck","fHistPurityCheck",5,binsPur,binsPurMin,binsPurMax);
    fOutputList->Add(fHistPurityCheck);
    fHistPurityCheck->Sumw2();
    fHistPurityCheck->GetAxis(0)->SetTitle("p_{T}");
    fHistPurityCheck->GetAxis(1)->SetTitle("mass");
    fHistPurityCheck->GetAxis(2)->SetTitle("V0 type");
    fHistPurityCheck->GetAxis(3)->SetTitle("check");
    fHistPurityCheck->GetAxis(1)->Set(602,binsMass);
    
    Int_t binsPurCorr[5] = {12,72,75,3,601};
    Double_t binsPurCorrMin[5] = {fPtTrigMin,-kPi/2,-2,0,0.44};
    Double_t binsPurCorrMax[5] = {15,-kPi/2+2*kPi,2.,3,1.15};
    
    fHistKorelPurCorr = new THnSparseF("fHistKorelPurCorr","fHistKorelPurCorr",5,binsPurCorr,binsPurCorrMin,binsPurCorrMax);
    fOutputList->Add(fHistKorelPurCorr);
    fHistKorelPurCorr->Sumw2();
    fHistKorelPurCorr->GetAxis(0)->SetTitle("p_{T}");
    fHistKorelPurCorr->GetAxis(1)->SetTitle("#Delta #Phi");
    fHistKorelPurCorr->GetAxis(2)->SetTitle("#Delta #Eta");
    fHistKorelPurCorr->GetAxis(3)->SetTitle("trigger");
    fHistKorelPurCorr->GetAxis(4)->SetTitle("mass");
    fHistKorelPurCorr->GetAxis(4)->Set(602,binsMass);
    
    fHistKorelResolCorr = new THnSparseF("fHistKorelResolCorr","fHistKorelResolCorr",5,binsPurCorr,binsPurCorrMin,binsPurCorrMax);
    fOutputList->Add(fHistKorelResolCorr);
    fHistKorelResolCorr->Sumw2();
    fHistKorelResolCorr->GetAxis(0)->SetTitle("p_{T}");
    fHistKorelResolCorr->GetAxis(1)->SetTitle("#Delta #Phi");
    fHistKorelResolCorr->GetAxis(2)->SetTitle("#Delta #Eta");
    fHistKorelResolCorr->GetAxis(3)->SetTitle("trigger");
    fHistKorelResolCorr->GetAxis(4)->SetTitle("mass");
    fHistKorelResolCorr->GetAxis(4)->Set(602,binsMass);
    
    fHistKorelPurCorrGen = new THnSparseF("fHistKorelPurCorrGen","fHistKorelPurCorrGen",5,binsPurCorr,binsPurCorrMin,binsPurCorrMax);
    fOutputList->Add(fHistKorelPurCorrGen);
    fHistKorelPurCorrGen->Sumw2();
    fHistKorelPurCorrGen->GetAxis(0)->SetTitle("p_{T}");
    fHistKorelPurCorrGen->GetAxis(1)->SetTitle("#Delta #Phi");
    fHistKorelPurCorrGen->GetAxis(2)->SetTitle("#Delta #Eta");
    fHistKorelPurCorrGen->GetAxis(3)->SetTitle("trigger");
    fHistKorelPurCorrGen->GetAxis(4)->SetTitle("mass");
    fHistKorelPurCorrGen->GetAxis(4)->Set(602,binsMass);
    
    fHistPtResolution = new TH3F("fHistPtResol","fHistPtResol",96,3,15,96,3,15,4,0,4);
    fOutputList->Add(fHistPtResolution);
    
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
        lPercentile = MultSelection->GetMultiplicityPercentile("V0M");
    }
    if ((lPercentile<0.)||(lPercentile>100.)) return;
    
    fHistMultipPercentile->Fill(lPercentile);

    TObjArray *mcTracksSel = new TObjArray;
    mcTracksSel->SetOwner(kTRUE);
    TObjArray *mcTracksTrigSel = new TObjArray;
    mcTracksTrigSel->SetOwner(kTRUE);
    TObjArray *mcTracksV0Sel = new TObjArray;
    mcTracksV0Sel->SetOwner(kTRUE);
    TObjArray *selectedMCTracks = new TObjArray;
    selectedMCTracks->SetOwner(kTRUE);
    TObjArray *selectedMCassoc = new TObjArray;
    selectedMCassoc->SetOwner(kTRUE);
    TObjArray *selectedMCtrig= new TObjArray;
    selectedMCtrig->SetOwner(kTRUE);
    TObjArray *selectedMCV0Triggersrec = new TObjArray;
    selectedMCV0Triggersrec->SetOwner(kTRUE);
    TClonesArray *mcArray = new TClonesArray;
    mcArray->SetOwner(kTRUE);
    TObjArray *selectedMCV0TriggersrecGoodId = new TObjArray;
    selectedMCV0TriggersrecGoodId->SetOwner(kTRUE);
    TObjArray *selectedMCV0TriggersrecGoodIdRec = new TObjArray;
    selectedMCV0TriggersrecGoodIdRec->SetOwner(kTRUE);
    TObjArray * selectedMCV0TriggersCanGen= new TObjArray;
    selectedMCV0TriggersCanGen->SetOwner(kTRUE);

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
 		TObjArray *mcTracks = new TObjArray;
		mcTracks->SetOwner(kTRUE);
 

		for (Int_t i = 0; i < nMCAllTracks; i++){ 
 			AliAODMCParticle *mcTrack = (AliAODMCParticle*)mcArray->At(i);
 			if (!mcTrack) {
				Error("ReadEventAODMC", "Could not receive particle %d", i);
				continue;
 			}
 			mcTracks->Add(mcTrack);
 		}

		Int_t nMCTracks = mcTracks->GetEntriesFast();
	
		for (Int_t iMC = 0; iMC<nMCTracks; iMC++){
			AliAODMCParticle *mcTrack = (AliAODMCParticle*)mcTracks->At(iMC);
			if (!mcTrack) {
				Error("ReadEventAODMC", "Could not receive particle %d", iMC);
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
                mcTracksSel->Add(new AliV0ChParticle(mcTrack->Eta(),mcTrack->Phi(),mcTrack->Pt(),4,mcTrack->GetLabel(),mcTrack->GetLabel(),kFALSE));
                
                if (mcTrackPt>fPtTrigMin) {
                    mcTracksTrigSel->Add(new AliV0ChParticle(mcTrack->Eta(),mcTrack->Phi(),mcTrack->Pt(),4,mcTrack->GetLabel(),kFALSE));
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
            
            IsFromCascade = (((motherPDG == 3222)|| (motherPDG==3212)|| (motherPDG==3112) || (motherPDG==3224) || (motherPDG==3214) || (motherPDG==3114) || (motherPDG==3322) || (motherPDG==3312)|| (motherPDG==3324) || (motherPDG==3314) || (motherPDG==3334)) && (mcMotherParticle->IsPhysicalPrimary()));

            Bool_t IsK0 = mcPartPdg==310;
            Bool_t IsLambda = mcPartPdg==3122;
            Bool_t IsAntiLambda = mcPartPdg==-3122;
            
            Bool_t isPhysPrim = mcTrack->IsPhysicalPrimary();

            IsK0 = IsK0 && (isPhysPrim);
            IsLambda = IsLambda && (isPhysPrim||IsFromCascade);
            IsAntiLambda = IsAntiLambda && (isPhysPrim||IsFromCascade);
            
            AliAODMCParticle* daughter0 = 0x0;
            AliAODMCParticle* daughter1 = 0x0;
            Int_t dau0 = mcTrack->GetDaughter(0);
            if (dau0>0) daughter0 = (AliAODMCParticle*) mcTracks->At(dau0);
            Int_t dau1 = mcTrack->GetDaughter(1);
            if (dau1>0) daughter1 = (AliAODMCParticle*) mcTracks->At(dau1);
            
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

            //MC V0 cuts

            if (mcTrack->Pt()>fPtTrigMin&&TrEtaMax){
                if(IsK0) {
                    mcTracksV0Sel->Add(new AliV0ChParticle(mcTrack->Eta(),mcTrack->Phi(),mcTrack->Pt(),1,mcTrack->GetLabel(),labelPos,labelNeg,kFALSE,mcTrack->M()));
                    if (fEfficiency){
                        Double_t v0effic[4]={mcTrack->Pt(),lPVz,0.5,mcTrack->Eta()};
                        fHistGenV0->Fill(v0effic); // for recunstruction efficiency calculation
                    }
                }
                if(IsLambda) {
                    mcTracksV0Sel->Add(new AliV0ChParticle(mcTrack->Eta(),mcTrack->Phi(),mcTrack->Pt(),2,mcTrack->GetLabel(),labelPos,labelNeg,kFALSE,mcTrack->M()));
                    if (fEfficiency){
                        Double_t v0effic[4]={mcTrack->Pt(),lPVz,1.5,mcTrack->Eta()};
                        fHistGenV0->Fill(v0effic); // for recunstruction efficiency calculation
                    }
                }
                if(IsAntiLambda) {
                    mcTracksV0Sel->Add(new AliV0ChParticle(mcTrack->Eta(),mcTrack->Phi(),mcTrack->Pt(),3,mcTrack->GetLabel(),labelPos,labelNeg,kFALSE,mcTrack->M()));
                    if (fEfficiency){
                        Double_t v0effic[4]={mcTrack->Pt(),lPVz,2.5,mcTrack->Eta()};
                        fHistGenV0->Fill(v0effic); // for recunstruction efficiency calculation
                    }
                }
                
            }

		}
        // MC closure test corellations
        if(fCorrelations){
            //V0-h
            Corelations(mcTracksV0Sel,mcTracksSel,fHistMCKorelacie, lPVz, fHistNumberOfTriggersGen,kFALSE,kTRUE,lPercentile,fHistPtHard,ptHard,kFALSE);

            //h-h
            Corelations(mcTracksTrigSel,mcTracksSel,fHistMCKorelacie, lPVz,fHistNumberOfTriggersGen, kTRUE, kFALSE,lPercentile,fHistPtHard,ptHard,kFALSE);
        }
        //reconstructed part. 
		Int_t nTracks = fAOD->GetNumberOfTracks();

		for (Int_t i = 0; i < nTracks; i++)
 	    {
            AliAODTrack* tr = dynamic_cast<AliAODTrack*>(fAOD->GetTrack(i));
            if(!tr) AliFatal("Not a standard AOD");
            selectedMCTracks->Add(tr);
        }

		Int_t nRecTracks = selectedMCTracks->GetEntriesFast();

		for (Int_t i = 0; i < nRecTracks; i++){
       		AliAODTrack* tras = (AliAODTrack*)selectedMCTracks->At(i);
       		if ((tras->Pt())<fPtAsocMin) continue;
        	if (!(IsMyGoodPrimaryTrack(tras))) continue;
        	Int_t AssocLabel = tras->GetLabel();
            if ((tras->Charge())==0) continue;
            
        	if (AssocLabel<=0) continue;
            AliAODMCParticle* mcTrack = static_cast<AliAODMCParticle*>(mcArray->At(AssocLabel));
        	Bool_t isPhyPrim = mcTrack->IsPhysicalPrimary();
            Double_t mcPt = tras->Pt();
            Double_t genPt = mcTrack->Pt();
            Double_t mcPhi = tras->Phi();
            Double_t mcEta = tras->Eta();
            Double_t genEta = mcTrack->Eta();
        	if (isPhyPrim) {
                
                if(fEfficiency) fHistRCPtAs->Fill(genPt,lPVz,genEta); // for recunstruction efficiency calculation
                selectedMCassoc->Add(new AliV0ChParticle(mcEta,mcPhi,mcPt,4,AssocLabel,tras->GetID(),kFALSE));
                if (mcPt>fPtTrigMin) {
                    selectedMCtrig->Add(new AliV0ChParticle(mcEta,mcPhi,mcPt,4,AssocLabel,kFALSE));
                    if(fEfficiency) fHistRCPtTrigg->Fill(genPt,lPVz,genEta); // for recunstruction efficiency calculation
                }
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
    
    for(Int_t i=0; i < iTracks; i++) {                 // loop over all these tracks
        AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(i));         // get a track (type AliAODTrack) from the event
        if(!track) continue;                            // if we failed, skip this track
        
        if(!IsMyGoodPrimaryTrack(track)) continue; // hybrid track selection
	   
		selectedTracks->Add(track);
        if(track->Pt()>fPtAsocMin) selectedAssociatedTracks-> Add(new AliV0ChParticle(track->Eta(), track->Phi(), track->Pt(), 4, 0,track->GetID(),kFALSE));
        if(track->Pt()>fPtTrigMin) selectedTriggerTracks-> Add(new AliV0ChParticle(track->Eta(), track->Phi(), track->Pt(), 4,0,track->GetID(),kFALSE));
	}

    TObjArray * selectedV0 = new TObjArray;
	selectedV0->SetOwner(kTRUE); 
	TObjArray * selectedV0Triggers = new TObjArray;
	selectedV0Triggers->SetOwner(kTRUE); 
	TObjArray * selectedV0Assoc = new TObjArray;
	selectedV0Assoc->SetOwner(kTRUE);

    
	for (Int_t i=0; i<nV0; i++){
        
        AliAODv0* V0 = static_cast<AliAODv0*>(fAOD->GetV0(i));
        if(!V0) continue;
        
        Double_t massK0=V0->MassK0Short();
        Double_t massLambda=V0->MassLambda();
        Double_t massAntilambda=V0->MassAntiLambda();

        Double_t ptTrig = V0->Pt();
        if(V0->Pt()<fPtTrigMin) continue; // pt trigger cut

        Bool_t k0 = ((massK0>0.44)&&(massK0<0.56));
        Bool_t Lambda = ((massLambda>1.08)&&(massLambda<1.15));
        Bool_t Antilambda = ((massAntilambda>1.08)&&(massAntilambda<1.15));

        // PID cut--------------------------
        Float_t nSigmaPosPion   = 0.;
        Float_t nSigmaNegPion   = 0.;
        Float_t nSigmaPosProton = 0.;
        Float_t nSigmaNegProton = 0.;

        const AliAODTrack *myTrackPos = new AliAODTrack();
        const AliAODTrack *myTrackNeg = new AliAODTrack();
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
                nSigmaNegPion = fPIDResponse->NumberOfSigmasTPC(myTrackPos, AliPID::kPion);
                nSigmaNegProton = fPIDResponse->NumberOfSigmasTPC(myTrackPos, AliPID::kProton);
            }
        }
        Bool_t bpPion = kTRUE;
        Bool_t bpProton = kTRUE;
        Bool_t bnPion = kTRUE;
        Bool_t bnProton = kTRUE;
        
        Bool_t cutK0Pid = kTRUE;
        Bool_t cutLambdaPid = kTRUE;
        Bool_t cutAntiLambdaPid = kTRUE;
        
            bpPion = TMath::Abs(nSigmaPosPion) <= 3.; // TPC dE/dx selection
            bpProton = TMath::Abs(nSigmaPosProton) <= 3.;
        
            bnPion = TMath::Abs(nSigmaNegPion) <= 3.;
            bnProton = TMath::Abs(nSigmaNegProton) <= 3.;
        
        
        cutK0Pid = (bpPion && bnPion);
        cutLambdaPid = (bpProton && bnPion);
        cutAntiLambdaPid = (bpPion && bnProton);
        
        if (k0) fHistK0MassPtCut->Fill(massK0,ptTrig,0.5);
        if (Lambda) fHistLambdaMassPtCut->Fill(massLambda,ptTrig,0.5);
        if (Antilambda) fHistAntiLambdaMassPtCut->Fill(massAntilambda,ptTrig,0.5);
        
        if (k0&&cutK0Pid) fHistK0MassPtCut->Fill(massK0,ptTrig,1.5);
        if (Lambda&&cutLambdaPid) fHistLambdaMassPtCut->Fill(massLambda,ptTrig,1.5);
        if (Antilambda&&cutAntiLambdaPid) fHistAntiLambdaMassPtCut->Fill(massAntilambda,ptTrig,1.5);
        
        Int_t oStatus = GetOStatus();
       // cout << oStatus << endl;
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
            
                        if (TMath::Abs(massK0-0.497614)>0.01){
                            fHistLambdaMassPtCut->Fill(massLambda,ptTrig,8.5);
            
                            if(fAnalysisMC){
                                FillMC(V0,mcArray,3122,2212, -211,2,massLambda,selectedMCV0Triggersrec,fHistRecV0,fHistLambdaMassPtCut,lPVz,myTrackPos,myTrackNeg,V0->GetOnFlyStatus(),fHistPurityCheck,selectedMCV0TriggersrecGoodId,selectedMCV0TriggersrecGoodIdRec,selectedMCV0TriggersCanGen,fHistPtResolution);
                            }
                            selectedV0Triggers-> Add(new AliV0ChParticle(V0->Eta(), V0->Phi(), V0->Pt(), 2,0,myTrackPos->GetID(),myTrackNeg->GetID(),V0->GetOnFlyStatus(),massLambda));
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
            
                        if(TMath::Abs(massLambda-1.115683)>0.02) {
                            fHistK0MassPtCut->Fill(massK0,ptTrig,8.5);
            
                            if(fAnalysisMC){
                                FillMC(V0,mcArray,310,211, -211,1,massK0,selectedMCV0Triggersrec,fHistRecV0,fHistK0MassPtCut,lPVz,myTrackPos,myTrackNeg,V0->GetOnFlyStatus(),fHistPurityCheck,selectedMCV0TriggersrecGoodId,selectedMCV0TriggersrecGoodIdRec,selectedMCV0TriggersCanGen,fHistPtResolution);
                            }
                            selectedV0Triggers-> Add(new AliV0ChParticle(V0->Eta(), V0->Phi(), V0->Pt(), 1,0,myTrackPos->GetID(),myTrackNeg->GetID(),V0->GetOnFlyStatus(),massK0));
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
                
                        if(TMath::Abs(massK0-0.497614)>0.01){
                            fHistAntiLambdaMassPtCut->Fill(massAntilambda,ptTrig,8.5);

            
                            if(fAnalysisMC){
                                FillMC(V0,mcArray,-3122,211, -2212,3,massAntilambda,selectedMCV0Triggersrec,fHistRecV0,fHistAntiLambdaMassPtCut,lPVz,myTrackPos,myTrackNeg,V0->GetOnFlyStatus(),fHistPurityCheck,selectedMCV0TriggersrecGoodId,selectedMCV0TriggersrecGoodIdRec,selectedMCV0TriggersCanGen,fHistPtResolution);
                            }
                            selectedV0Triggers-> Add(new AliV0ChParticle(V0->Eta(), V0->Phi(), V0->Pt(), 3,0,myTrackPos->GetID(),myTrackNeg->GetID(),V0->GetOnFlyStatus(),massAntilambda));
                        }
                    }
                }
            }
        }
    }

	// Corelation ==========================================

    
    if(fAnalysisMC&&fEfficiency){
        //for puriry correction
        Corelations(selectedMCV0TriggersrecGoodId,selectedMCassoc,fHistKorelPurCorr, lPVz, fHistNumberOfTriggersRec,kFALSE,kTRUE,lPercentile,fHistPtHard,ptHard,kTRUE);
        //for pt resolution correction
        Corelations(selectedMCV0TriggersrecGoodIdRec,selectedMCassoc,fHistKorelResolCorr, lPVz, fHistNumberOfTriggersRec,kFALSE,kTRUE,lPercentile,fHistPtHard,ptHard,kTRUE);
        //for purity correction
        Corelations(selectedMCV0TriggersCanGen,selectedMCassoc,fHistKorelPurCorrGen, lPVz, fHistNumberOfTriggersRec,kFALSE,kTRUE,lPercentile,fHistPtHard,ptHard,kTRUE);
    }
    
     if(fAnalysisMC&&fCorrelations){
        //V0-h MC rec
        Corelations(selectedMCV0Triggersrec,selectedMCassoc,fHistKorelacieMCrec, lPVz, fHistNumberOfTriggersRec,kFALSE,kTRUE,lPercentile,fHistPtHard,ptHard,kFALSE);

        //h-h MC rec
        Corelations(selectedMCtrig,selectedMCassoc,fHistKorelacieMCrec, lPVz, fHistNumberOfTriggersRec,kTRUE,kFALSE,lPercentile,fHistPtHard,ptHard,kFALSE);
    } else if(fCorrelations){
        //Data V0-h
        Corelations(selectedV0Triggers,selectedAssociatedTracks,fHistKorelacie,lPVz,fHistNumberOfTriggers,kFALSE,kTRUE,lPercentile,fHistPtHard,ptHard,kFALSE);

	    //Data h-h
        Corelations(selectedTriggerTracks,selectedAssociatedTracks,fHistKorelacie,lPVz,fHistNumberOfTriggers,kFALSE,kFALSE,lPercentile,fHistPtHard,ptHard,kFALSE);
    }


 	// Mixing ==============================================

    fHistMultVtxz->Fill(lPercentile,lPVz);

	fPool = fPoolMgr->GetEventPool(lPercentile, lPVz);
    if (!fPool) {
        AliWarning(Form("No pool found for centrality = %f, zVtx = %f", lPercentile, lPVz));
        return;
    }
 	Int_t nMix = fPool->GetCurrentNEvents();
	if (fPool->IsReady() || fPool->NTracksInPool() > fMixingTracks / 5 || nMix >= 5)
	{
		
        for (Int_t jMix=0; jMix<nMix; jMix++)
 		{// loop through mixing events
 			TObjArray* bgTracks = fPool->GetEvent(jMix);
            if(fAnalysisMC&&fCorrelations) {
                
                CorelationsMixing(selectedMCV0Triggersrec,bgTracks,fHistMCMixingRec,lPVz,lPercentile);
                CorelationsMixing(selectedMCtrig,bgTracks,fHistMCMixingRec,lPVz,lPercentile);
            }else if(fCorrelations){
                CorelationsMixing(selectedV0Triggers,bgTracks,fHistdPhidEtaMix,lPVz,lPercentile);
                CorelationsMixing(selectedTriggerTracks,bgTracks,fHistdPhidEtaMix,lPVz,lPercentile);
            }
		 }
	}
	TObjArray* cloneArray = (TObjArray *)selectedTracks->Clone();
	cloneArray->SetOwner(kTRUE);
	fPool->UpdatePool(cloneArray);
    
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
          if (TMath::Abs(t->Eta())>=0.8) return kFALSE;
		  if (!t->TestFilterBit(768)) return kFALSE;
  
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
		//if(TMath::Abs(t->RapLambda())>=0.5) return kFALSE;
    //Pseudorap
    if(TMath::Abs(t->Eta())>=0.8) return kFALSE;
		
		return kTRUE;
}
//____________________________________________________________________________
Bool_t AliAnalysisTaskDiHadCorrelHighPt::IsMyGoodV0RapidityK0(const AliAODv0 *t)
{
		//Rapidity
		//if(TMath::Abs(t->RapK0Short())>=0.5) return kFALSE;
    //Pseudorap
    if(TMath::Abs(t->Eta())>=0.8) return kFALSE;
		
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
		
		if(lifetime>=20) return kFALSE;
		
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
		
		if(lifetime>=30) return kFALSE;
		
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
    
    if(lifetime>=30) return kFALSE;
    
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
    
    if (TMath::Abs(t->Eta())>=0.8) return kFALSE;
		
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
	if(v0->RadiusV0()<=0.5) return kFALSE;
	
	return kTRUE;
}
//____________________________________________________________________________
void AliAnalysisTaskDiHadCorrelHighPt::Corelations(TObjArray *triggers, TObjArray *associated, THnSparse * fHistKor, Double_t lPVz, THnSparse* fHistNumOfTrig,Bool_t hhMC,Bool_t V0h,Float_t perc,TH3F *fHistPtHard, Double_t ptHard,Bool_t godId){

    const Double_t kPi = TMath::Pi();
    Int_t nAssoc = associated->GetEntriesFast();
    Int_t nTrig = triggers->GetEntriesFast();
    Double_t status = 0.;

    for (Int_t i=0; i<nTrig; i++){
        AliV0ChParticle* trig = (AliV0ChParticle*)  triggers->At(i);
        if (TMath::Abs(trig->Eta())>=0.8) continue;
        
        if (trig->GetRecStatus()) status=0.5;
        if (!trig->GetRecStatus()) status=1.5;
        if (ptHard!=0&&!godId) fHistPtHard->Fill(trig->Pt()/ptHard,perc,trig->WhichCandidate()-0.5);
        Double_t massTrig = 0.;
        if(trig->WhichCandidate()<4) massTrig=trig->GetMass();
        
        if(!godId){
            Double_t triggers[6]={trig->Pt(),lPVz,trig->Eta(),trig->WhichCandidate()-0.5,massTrig,perc};
       
            fHistNumOfTrig->Fill(triggers);
        }
        for (Int_t j=0; j<nAssoc; j++){
            AliV0ChParticle* assoc = (AliV0ChParticle*)  associated->At(j);

            Double_t deltaEta = trig->Eta() - assoc->Eta();
            Double_t deltaPhi = trig->Phi() - assoc->Phi();
            if (deltaPhi > (1.5*kPi)) deltaPhi -= 2.0*kPi;
            if (deltaPhi < (-0.5*kPi)) deltaPhi += 2.0*kPi;

            if(trig->Pt()<=assoc->Pt()) continue;
            
            //removing autocorrelations
            if(V0h){
                
                Int_t negID = trig->GetIDNeg();
                Int_t posID = trig->GetIDPos();
                Int_t atrID = assoc->GetIDCh();
                
                if ((TMath::Abs(negID))==(TMath::Abs(atrID))) continue;
                if ((TMath::Abs(posID))==(TMath::Abs(atrID))) continue;
                
            }

            Int_t labelTrig = -2;
            Int_t labelAssoc =0;
            if(hhMC){
                labelTrig=trig->MyLabel();
                labelAssoc=assoc->MyLabel();
                
            }
            
            if(labelTrig==labelAssoc) continue;
            
            if(godId){
                Double_t korel[5] = {trig->Pt(),deltaPhi,deltaEta,trig->WhichCandidate()-0.5,massTrig}; //histogram for contamination correction
                fHistKor->Fill(korel);
            }
            else{
                Double_t korel[10] = {trig->Pt(),assoc->Pt(),deltaPhi,deltaEta, lPVz,trig->WhichCandidate()-0.5, trig->Eta(),assoc->Eta(),massTrig,perc};
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

    for (Int_t i=0; i<nTrig; i++){
        AliV0ChParticle* trig = (AliV0ChParticle*)  triggers->At(i);
        
        Double_t massTrig = 0.;
        if(trig->WhichCandidate()<4) massTrig=trig->GetMass();
        if(trig->WhichCandidate()==1&&(massTrig<0.486||massTrig>0.509)) continue;
        if((trig->WhichCandidate()==2||trig->WhichCandidate()==3)&&(massTrig<1.1112||massTrig>1.12)) continue;

        for (Int_t j=0; j<nAssoc; j++){
             AliAODTrack* assoc = (AliAODTrack*) bgTracks->At(j);

             if (( (assoc->Pt())>=trig->Pt() ) || ( (assoc->Pt())<fPtAsocMin )) continue;

             if(!IsMyGoodPrimaryTrack(assoc)) continue; 

             Double_t   deltaEta = trig->Eta() - assoc->Eta();
             Double_t   deltaPhi = trig->Phi() - assoc->Phi();
             Double_t   assocPt = assoc->Pt();

            if (deltaPhi > (1.5*kPi)) deltaPhi -= 2.0*kPi;
            if (deltaPhi < (-0.5*kPi)) deltaPhi += 2.0*kPi;
            
            Double_t korel[7] = {trig->Pt(),assocPt,deltaPhi,deltaEta, lPVz,trig->WhichCandidate()-0.5,perc};
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
void AliAnalysisTaskDiHadCorrelHighPt::FillMC(const AliAODv0 *V0,TClonesArray *mcArray,Int_t pdgV0,Int_t pdgDau1, Int_t pdgDau2,Int_t triggerType, Double_t mass, TObjArray * selectedMCV0Triggersrec,THnSparse * fHistRecV0, TH3F * fHistMassPtCut,Double_t lPVz, const AliAODTrack * myTrackPos,const AliAODTrack * myTrackNeg,Bool_t status,THnSparse * histPur, TObjArray * selectedMCV0TriggersrecGoodId, TObjArray * selectedMCV0TriggersrecGoodIdrec, TObjArray * selectedMCV0TriggersrecGen,TH3F * fHistresol){
    
    if(fPurityCheck){
        Double_t purity[5] ={V0->Pt(),mass,triggerType-0.5,0.5,-1};
        histPur->Fill(purity);
    }
    
    selectedMCV0Triggersrec-> Add(new AliV0ChParticle(V0->Eta(), V0->Phi(), V0->Pt(), triggerType,0,myTrackPos->GetID(),myTrackNeg->GetID(),status,mass)); // all reconstructed candidates for raw correlation function, with reconstructed pt
    
    
    Int_t myTrackPosLabel = TMath::Abs(myTrackPos->GetLabel());
    Int_t myTrackNegLabel = TMath::Abs(myTrackNeg->GetLabel());
    
    AliAODMCParticle *mcPosTrack = (AliAODMCParticle*)mcArray->At(myTrackPosLabel);
    if (!mcPosTrack) return;
    Int_t PosTrackPdg = mcPosTrack->GetPdgCode();
    AliAODMCParticle *mcNegTrack = (AliAODMCParticle*)mcArray->At(myTrackNegLabel);
    if (!mcNegTrack) return;
    Int_t NegTrackPdg = mcNegTrack->GetPdgCode();
    
    if(fPurityCheck){
        Double_t puri[5] ={V0->Pt(),mass,triggerType-0.5,1.5,-1};
        histPur->Fill(puri);
    }
    
    Int_t myTrackPosMotherLabel = mcPosTrack->GetMother();
    Int_t myTrackNegMotherLabel = mcNegTrack->GetMother();
    
    if ((myTrackPosMotherLabel==-1)||(myTrackNegMotherLabel==-1)) return;
    
    AliAODMCParticle *GenV0 = (AliAODMCParticle*)mcArray->At(myTrackPosMotherLabel);
    selectedMCV0TriggersrecGen-> Add(new AliV0ChParticle(GenV0->Eta(), GenV0->Phi(), GenV0->Pt(), triggerType,0,myTrackPos->GetID(),myTrackNeg->GetID(),status,mass)); // all reconstructed candidates with generatated pt, for contamination correction
    
    if(fPurityCheck){
        Double_t pur[5] ={V0->Pt(),mass,triggerType-0.5,2.5,-1};
        histPur->Fill(pur);
    }
    
    if (myTrackPosMotherLabel!=myTrackNegMotherLabel) return;
    if(fPurityCheck){
        Double_t pu[5] ={V0->Pt(),mass,triggerType-0.5,3.5,-1};
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
        IsFromCascade = (((MoMPdg == 3222)|| (MoMPdg==3212)|| (MoMPdg==3112) || (MoMPdg==3224) || (MoMPdg==3214) || (MoMPdg==3114) || (MoMPdg==3322) || (MoMPdg==3312)|| (MoMPdg==3324) || (MoMPdg==3314) || (MoMPdg==3334)) && (mcPosMotherOfMother->IsPhysicalPrimary()));
    }
    if(fPurityCheck){
        Double_t purit[5] ={V0->Pt(),mass,triggerType-0.5,4.5,-1};
        histPur->Fill(purit);
    }
    Bool_t isGoodID = (MotherPdg==pdgV0);
    Bool_t isFromMaterial = mcPosMother->IsSecondaryFromMaterial();

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
        Double_t purity[5] ={V0->Pt(),mass,triggerType-0.5,5.5,ident-0.5};
        histPur->Fill(purity);
    }
    
    if(isGoodID&&fPurityCheck){
        Double_t purity[5] ={V0->Pt(),mass,triggerType-0.5,6.5,-1};
        histPur->Fill(purity);
    }

    if(!isFromMaterial&&isGoodID&&fPurityCheck){
        Double_t purity[5] ={V0->Pt(),mass,triggerType-0.5,7.5,-1};
        histPur->Fill(purity);
    }
    
    Bool_t IsParticleFromMC = kFALSE;
    
    if(pdgV0==310) IsParticleFromMC= (IsFromMC&&IsPhysPrim);
    else {
        IsParticleFromMC = (IsFromMC&&(IsPhysPrim||IsFromCascade));
    }
    
    Double_t V0mcPt = mcPosMother->Pt();
    if(V0mcPt<=fPtTrigMin) return;
   
    if(IsParticleFromMC){
        fHistresol->Fill(V0mcPt,V0->Pt(),triggerType-0.5);
    Double_t V0mcEta = mcPosMother->Eta();
    Double_t V0mcPhi = mcPosMother->Phi();
    selectedMCV0TriggersrecGoodId-> Add(new AliV0ChParticle(V0mcEta, V0mcPhi, V0mcPt, triggerType,0,myTrackPos->GetID(),myTrackNeg->GetID(),status,mass)); // good identified V0 with generated pt, for contamination correction
    selectedMCV0TriggersrecGoodIdrec-> Add(new AliV0ChParticle(V0->Eta(), V0->Phi(), V0->Pt(), triggerType,0,myTrackPos->GetID(),myTrackNeg->GetID(),status,mass)); // good identified V0 with reconstructed pt, for resolution correction
    
    if(fEfficiency){
        Double_t v0effic[6]={V0mcPt,lPVz,triggerType-0.5,V0mcEta,mass};
        fHistRecV0->Fill(v0effic);
    }

    fHistMassPtCut->Fill(mass,V0mcPt,7.5);
    }

}

