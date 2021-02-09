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
 * The task selects candidates for K0s, Lambdas and AntiLambdas
 * and calculates correlations with charged unidentified particles in phi and eta.
 * The charged unidentified particles are also taken as trigger particles to have a check.
 * The task works with AOD or ESD (with or without MC info) events only and containes also mixing for acceptance corrections.
 * Last update edited by Lucia Anna Husova, January 2021
 */

#include <TChain.h>
#include <TH1F.h>
#include <TH3F.h>
#include <TH1D.h>
#include <TList.h>
#include <AliAnalysisTask.h>
#include <AliAnalysisManager.h>
#include "AliAODTrack.h"
#include "AliESDtrack.h"
#include <AliAODEvent.h>
#include <AliESDEvent.h>
#include "AliAODv0.h"
#include "AliESDv0.h"
#include "AliESDcascade.h"
#include <AliAODInputHandler.h>
#include "AliAnalysisTaskDiHadCorrelHighPt.h"
#include <TMath.h>
#include <AliMultiEventInputHandler.h>
#include <AliPIDResponse.h>
#include <AliAODPid.h>
#include <THnSparse.h>
#include <AliAODVertex.h>
#include <AliESDVertex.h>
#include <AliEventPoolManager.h>
#include "AliAODMCParticle.h"
#include "AliMultSelection.h"
#include "AliPPVsMultUtils.h"
#include <AliESDtrackCuts.h>
#include "TDatabasePDG.h"
#include "AliEventCuts.h"

class AliPIDResponse;
class AliMultSelection;
class AliAnalysisTaskDiHadCorrelHighPt;    // analysis class

using namespace std;            // std namespace: so you can do things like 'cout'

ClassImp(AliAnalysisTaskDiHadCorrelHighPt) // classimp: necessary for root

AliAnalysisTaskDiHadCorrelHighPt::AliAnalysisTaskDiHadCorrelHighPt() : AliAnalysisTaskSE(),
    fAliEventCuts(0),
    fAOD(0),
    fESD(0),
    fmcEvent(0),
    fPIDResponse(0),
    fOutputList(0),
    fHistLambdaMassPtCut(0),
    fHistK0MassPtCut(0),
    fHistAntiLambdaMassPtCut(0),
	fHistKorelacie(0),
    fHistdPhidEtaMix(0),
    fHistV0Multiplicity(0),
    fHistMultVtxz(0),
    fHistMCPtAs(0),
    fHistRCPtAs(0),
    fHistNumberOfTriggers(0),
    fHistMCKorelacie(0),
    fHistMCMixingRec(0),
    fHistMCMixingGen(0),
    fFillMixed(0),
    fMixingTracks(5000),
	fPoolMgr(0x0),
    fPool(0x0),
    fPoolMCGen(0x0),
    fAnalysisMC(kTRUE),
    fOStatus(0),
    fPtTrigMin(3),
    fPtAsocMin(1),
    fCutsCrosscheck(kFALSE),
    fMixedEvents(10),
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
    fNumberOfPtBinsTrigger(17),
    fNumberOfPtBinsAssoc(19),
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
    fFilterBit(32),
    fRemoveLamhFromCascade(kFALSE),
    fRemoveHadrFromV0(kTRUE),
    fAacceptLambdasFromCasscade(kFALSE),
    fPurePrimHadrons(kFALSE),
    fPureV0(kFALSE),
    fAnalysisAOD(kTRUE),
    fNumberOfDeltaPhiBins(72),
    fNumberOfDeltaEtaBins(75),
    fNumberOfEtaBins(40),
    fNumberPhiBins(72),
    fMixing(kTRUE),
    fMixingGen(kFALSE),
    fNumOfVzBins(9),
    fPrimaryVertexCut(10),
    fHistGenMultiplicity(0),
    fHistTPCTracksVsClusters(0),
    fHistVZeroPercentileTPCMult(0),
    fESDTrackCuts(nullptr),
    fTestPions(kFALSE),
    fHistPosNegTracks(0),
    fHistLambdaFeedDown(0),
    fAnalyseFeedDown(kFALSE),
    fPtAssocMax(20),
    fPtTrigMax(20),
    fHistXiMinusMassPtCut(0),
    fHistXiPlusMassPtCut(0),
    fHistCasMC(0),
    fMagneticField(0),
    fMaxDCAToVertexZ(2),
    fMaxDCAToVertexXY(0.3),
    fMinNCrossedRowsTPCprimtracks(70),
    fmcTracksSel(0),
    fmcGenTracksMixing(0),
    fmcTracksTrigSel(0),
    fmcTracksV0Sel(0),
    fmcV0AssocSel(0),
    fselectedMCassoc(0),
    fselectedMCV0assoc(0),
    fselectedMCtrig(0),
    fselectedMCV0Triggersrec(0),
    fselectedTracks(0),
    fselectedAssociatedTracks(0),
    fselectedTriggerTracks(0),
    fselectedV0Triggers(0),
    fselectedV0Assoc(0),
    fHistEffCorrectionHadron(0),
    fHistEffCorrectionK0(0),
    fHistEffCorrectionLam(0),
    fHistEffCorrectionAntiLam(0),
    fHistEffCorrectionNegXi(0),
    fHistEffCorrectionPosXi(0),
    fHistSecondaryCont(0),
    fEffList(0),
    fUseEff(kFALSE),
    fMixCorrect(kTRUE),
    fminBias(kTRUE),
    fhighMult(kFALSE),
    fhighMultSPD(kFALSE),
    fHistMultVZEROTracklets(0),
    fPercentileMin(0.),
    fPercetileMax(100.),
    fEventCutsQAPlots(kFALSE),
    fpp(kTRUE),
    fNMultiplicityBins(10)
{
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty

    for(Int_t i = 0; i < 3; i++ ) { fPV        [i ] = -1.; }
}
//_____________________________________________________________________________
AliAnalysisTaskDiHadCorrelHighPt::AliAnalysisTaskDiHadCorrelHighPt(const char *name, Bool_t analysisMC, Bool_t useeff): AliAnalysisTaskSE(name),
    fAliEventCuts(0),
    fAOD(0),
    fESD(0),
    fmcEvent(0),
    fPIDResponse(0),
    fOutputList(0),
    fHistLambdaMassPtCut(0),
    fHistK0MassPtCut(0),
    fHistAntiLambdaMassPtCut(0),
    fHistKorelacie(0),
    fHistdPhidEtaMix(0),
    fHistV0Multiplicity(0),
    fHistMultVtxz(0),
    fHistMCPtAs(0),
    fHistRCPtAs(0),
    fHistNumberOfTriggers(0),
    fHistMCKorelacie(0),
    fHistMCMixingRec(0),
    fHistMCMixingGen(0),
    fFillMixed(0),
    fMixingTracks(5000),
    fPoolMgr(0x0),
    fPool(0x0),
    fPoolMCGen(0x0),
    fAnalysisMC(analysisMC),
    fOStatus(0),
    fPtTrigMin(3),
    fPtAsocMin(1),
    fCutsCrosscheck(kFALSE),
    fMixedEvents(10),
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
    fNumberOfPtBinsTrigger(17),
    fNumberOfPtBinsAssoc(19),
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
    fFilterBit(32),
    fRemoveLamhFromCascade(kFALSE),
    fRemoveHadrFromV0(kTRUE),
    fAacceptLambdasFromCasscade(kFALSE),
    fPurePrimHadrons(kFALSE),
    fPureV0(kFALSE),
    fAnalysisAOD(kTRUE),
    fNumberOfDeltaPhiBins(72),
    fNumberOfDeltaEtaBins(75),
    fNumberOfEtaBins(40),
    fNumberPhiBins(72),
    fMixing(kTRUE),
    fMixingGen(kFALSE),
    fNumOfVzBins(9),
    fPrimaryVertexCut(10),
    fHistGenMultiplicity(0),
    fHistTPCTracksVsClusters(0),
    fHistVZeroPercentileTPCMult(0),
    fESDTrackCuts(nullptr),
    fTestPions(kFALSE),
    fHistPosNegTracks(0),
    fHistLambdaFeedDown(0),
    fAnalyseFeedDown(kFALSE),
    fPtAssocMax(20),
    fPtTrigMax(20),
    fHistXiMinusMassPtCut(0),
    fHistXiPlusMassPtCut(0),
    fHistCasMC(0),
    fMagneticField(0),
    fMaxDCAToVertexZ(2),
    fMaxDCAToVertexXY(0.3),
    fMinNCrossedRowsTPCprimtracks(70),
    fmcTracksSel(0),
    fmcGenTracksMixing(0),
    fmcTracksTrigSel(0),
    fmcTracksV0Sel(0),
    fmcV0AssocSel(0),
    fselectedMCassoc(0),
    fselectedMCV0assoc(0),
    fselectedMCtrig(0),
    fselectedMCV0Triggersrec(0),
    fselectedTracks(0),
    fselectedAssociatedTracks(0),
    fselectedTriggerTracks(0),
    fselectedV0Triggers(0),
    fselectedV0Assoc(0),
    fHistEffCorrectionHadron(0),
    fHistEffCorrectionK0(0),
    fHistEffCorrectionLam(0),
    fHistEffCorrectionAntiLam(0),
    fHistEffCorrectionNegXi(0),
    fHistEffCorrectionPosXi(0),
    fHistSecondaryCont(0),
    fEffList(0),
    fUseEff(useeff),
    fMixCorrect(kTRUE),
    fminBias(kTRUE),
    fhighMult(kFALSE),
    fhighMultSPD(kFALSE),
    fHistMultVZEROTracklets(0),
    fPercentileMin(0.),
    fPercetileMax(100.),
    fEventCutsQAPlots(kFALSE),
    fpp(kTRUE),
    fNMultiplicityBins(10)
{
    // constructor

    for(Int_t i = 0; i < 3; i++ ) { fPV[i] = -1.; }

    if(fUseEff) {
        DefineInput(1, TList::Class());
    }

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
    if(fESDTrackCuts){delete fESDTrackCuts; fESDTrackCuts = nullptr;}
    if(fOutputList) {
        delete fOutputList;     // at the end of your task, it is deleted from memory by calling this function
    }
    if(fmcTracksSel) delete fmcTracksSel;
    if(fmcGenTracksMixing) delete fmcGenTracksMixing;
    if(fmcTracksTrigSel) delete fmcTracksTrigSel;
    if(fmcTracksV0Sel) delete fmcTracksV0Sel;
    if(fmcV0AssocSel) delete fmcV0AssocSel;
    if(fselectedMCassoc) delete fselectedMCassoc;
    if(fselectedMCV0assoc) delete fselectedMCV0assoc;
    if(fselectedMCtrig) delete fselectedMCtrig;
    if(fselectedMCV0Triggersrec) delete fselectedMCV0Triggersrec;
    if(fselectedTracks) delete fselectedTracks;
    if(fselectedAssociatedTracks) delete fselectedAssociatedTracks;
    if(fselectedTriggerTracks) delete fselectedTriggerTracks;
    if(fselectedV0Triggers) delete fselectedV0Triggers;
    if(fselectedV0Assoc) delete fselectedV0Assoc;
}
//_____________________________________________________________________________
void AliAnalysisTaskDiHadCorrelHighPt::UserCreateOutputObjects()
{
	const Double_t kPi = TMath::Pi();

   if(fUseEff) {
        fEffList = dynamic_cast<TList*>(GetInputData(1));
    }
    if(fEffList){
      fHistEffCorrectionK0 = (THnSparseF*)fEffList->FindObject("fHistEffCorrectionK0");
      fHistEffCorrectionLam = (THnSparseF*)fEffList->FindObject("fHistEffCorrectionLam");
      fHistEffCorrectionAntiLam = (THnSparseF*)fEffList->FindObject("fHistEffCorrectionAntiLam");
      fHistEffCorrectionHadron = (THnSparseF*)fEffList->FindObject("fHistEffCorrectionHadron");
      fHistEffCorrectionNegXi = (THnSparseF*)fEffList->FindObject("fHistEffCorrectionNegXi");
      fHistEffCorrectionPosXi = (THnSparseF*)fEffList->FindObject("fHistEffCorrectionPosXi");
      fHistSecondaryCont = (TH1F*)fEffList->FindObject("fHistSecondaryCont");

        if(!fHistEffCorrectionK0 || !fHistEffCorrectionLam || !fHistEffCorrectionAntiLam || !fHistEffCorrectionHadron|| !fHistEffCorrectionNegXi || !fHistEffCorrectionPosXi || !fHistSecondaryCont){
            if(fCorrelations){
                if (!fHistEffCorrectionHadron|| !fHistSecondaryCont) {
                    cout<<"Efficiency histograms hadrons are not available!"<<endl;
                    return;
                }
                if(fV0hCorr||fhV0Corr){
                    if(!fHistEffCorrectionK0 || !fHistEffCorrectionLam || !fHistEffCorrectionAntiLam) return;
                }
                if(fAnalyseFeedDown){
                    if(!fHistEffCorrectionNegXi || !fHistEffCorrectionPosXi) return;
                }
            }
            if(fMixing&&fUseEff){
                if (!fHistEffCorrectionHadron|| !fHistSecondaryCont) {
                    cout<<"Efficiency histograms hadrons are not available for mixing!"<<endl;
                    return;
                }
            }
        }
    }

	Float_t kPtBins[20] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20};
	Int_t kMassBins = 500;
	Float_t kMassMinK = 0.4;
	Float_t kMassMaxK = 0.6;

	Float_t kMassMinLambda = 0.8;
	Float_t kMassMaxLambda = 1.2;

    Float_t kMassMinXi = 1.285;
    Float_t kMassMaxXi = 1.355;

	Float_t kMassBinsK[kMassBins];
	Float_t kMassBinsLambda[kMassBins];
    Float_t kMassBinsXi[kMassBins];

	kMassBinsK[0] = kMassMinK;
	kMassBinsLambda[0] = kMassMinLambda;
    kMassBinsXi[0] = kMassMinXi;
	
	for (Int_t i=0; i<kMassBins; i++){
		kMassBinsK[i+1] = kMassBinsK[i] + (kMassMaxK-kMassMinK)/kMassBins;
		kMassBinsLambda[i+1] = kMassBinsLambda[i] + (kMassMaxLambda-kMassMinLambda)/kMassBins;
        kMassBinsXi[i+1] = kMassBinsXi[i] + (kMassMaxXi-kMassMinXi)/kMassBins;
	}

	Int_t kNCuts = 20;
	Float_t kCuts[kNCuts+1];
	kCuts[0]=0;
	for (Int_t i=0; i<kNCuts; i++){
		kCuts[i+1]=kCuts[i]+1;
	}

    Int_t bins[9]= {fNumberOfPtBinsTrigger,fNumberOfPtBinsAssoc,fNumberOfDeltaPhiBins,fNumberOfDeltaEtaBins,fNumOfVzBins,12,902,fNMultiplicityBins,2};
    Double_t min[9] = {fPtTrigMin,fPtAsocMin, -kPi/2, -2., -10., 0.,0.44,0,-2};
    Double_t max[9] = {fPtTrigMax, fPtAssocMax, -kPi/2+2*kPi, 2., 10., 12., 1.355,100,2};
    
	Int_t  NofCentBins  = 10;
    Double_t MBins[]={0,10,20,30,40,50,60,70,80,90,100};
         
    const Int_t NofZVrtxBins  =  fNumOfVzBins;
    Double_t ZBins[NofZVrtxBins+1];
    ZBins[0]=-1*fPrimaryVertexCut;
    if(NofZVrtxBins==9) {
        ZBins[1]=-7.0;
        for(Int_t i=2;i<NofZVrtxBins;i++){
            ZBins[i]=ZBins[i-1]+2;
        }
        ZBins[9]=fPrimaryVertexCut;
    }
    else{
        Double_t binstep = Double_t(2*fPrimaryVertexCut)/NofZVrtxBins;
        for(Int_t i=1; i<NofZVrtxBins+1; i++){
            ZBins[i]=ZBins[i-1]+binstep;
        }
    }

	Int_t bins2d[5] = {fNumberOfPtBinsTrigger,fNumOfVzBins,7,902,fNMultiplicityBins};
	Double_t mis2d[5] = {fPtTrigMin,-10,0.,0.44,0};
	Double_t maxs2d[5] = {fPtTrigMax,10,7.,1.355,100};
	
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
        
    fHistK0MassPtCut = new TH3F ("fHistK0MassPtCut", "fHistK0MassPtCut", kMassBins, kMassBinsK, 19, kPtBins, kNCuts, kCuts);
    fHistLambdaMassPtCut = new TH3F("fHistLambdaMassPtCut", "fHistLambdaMassPtCut", kMassBins, kMassBinsLambda, 19, kPtBins, kNCuts, kCuts);
	fHistAntiLambdaMassPtCut = new TH3F("fHistAntiLambdaMassPtCut", "fHistAntiLambdaMassPtCut", kMassBins, kMassBinsLambda, 19, kPtBins, kNCuts, kCuts);
    fHistXiMinusMassPtCut = new TH3F ("fHistXiMinusMassPtCut", "fHistXiMinusMassPtCut", kMassBins, kMassBinsXi, 19, kPtBins, kNCuts, kCuts);
    fHistXiPlusMassPtCut = new TH3F ("fHistXiPlusMassPtCut", "fHistXiPlusMassPtCut", kMassBins, kMassBinsXi, 19, kPtBins, kNCuts, kCuts);
    fOutputList->Add(fHistK0MassPtCut); 
    fOutputList->Add(fHistLambdaMassPtCut); 
	fOutputList->Add(fHistAntiLambdaMassPtCut);
    fOutputList->Add(fHistXiPlusMassPtCut);
    fOutputList->Add(fHistXiMinusMassPtCut);

    fHistCasMC = new TH2D ("fHistCasMC","fHistCasMC",2,5,7,11,0,11);
    fOutputList->Add(fHistCasMC); 

	fHistKorelacie = new THnSparseF ("fHistKorelacie","fHistKorelacie", 9, bins, min, max);
    fHistKorelacie->GetAxis(0)->SetTitle("p_{T}^{trig}");
    fHistKorelacie->GetAxis(1)->SetTitle("p_{T}^{assoc}");
    fHistKorelacie->GetAxis(2)->SetTitle("#Delta#phi");
    fHistKorelacie->GetAxis(3)->SetTitle("#Delta#eta");
    fHistKorelacie->GetAxis(4)->SetTitle("p_{vz}");
    fHistKorelacie->GetAxis(5)->SetTitle("trigger");
    fHistKorelacie->GetAxis(6)->SetTitle("mass");
    fHistKorelacie->GetAxis(7)->SetTitle("multiplicity percentile");
    fHistKorelacie->GetAxis(8)->SetTitle("hadron charge");
    fOutputList->Add(fHistKorelacie);
    fHistKorelacie->Sumw2();
    
   Double_t *binsMass = new Double_t [903];
    binsMass[0]=0.44;
    for(Int_t i=0;i<901;i++){
        if(i<301) binsMass[i+1] = binsMass[i]+(0.56-0.44)/300;
        if(i==301) binsMass[i+1]=1.08;
        if(i>301&&i<602) binsMass[i+1] = binsMass[i]+(1.15-1.08)/300;
        if(i==602) binsMass[i+1]=1.285;
        if(i>602) binsMass[i+1] = binsMass[i]+(1.355-1.285)/300;
    }
    binsMass[902]=1.355;
    fHistKorelacie->GetAxis(6)->Set(902,binsMass);
    fHistKorelacie->GetAxis(4)->Set(NofZVrtxBins,ZBins);
    Double_t binsMult[17]={0,1,2,3,5,7,10,15,20,30,40,50,60,70,80,90,100};
    if(fNMultiplicityBins==16)fHistKorelacie->GetAxis(7)->Set(16,binsMult);
    
	fHistdPhidEtaMix = new THnSparseF ("fHistdPhidEtaMix", "fHistdPhidEtaMix", 9, bins, min, max);
    fHistdPhidEtaMix->GetAxis(0)->SetTitle("p_{T}^{trig}");
    fHistdPhidEtaMix->GetAxis(1)->SetTitle("p_{T}^{assoc}");
    fHistdPhidEtaMix->GetAxis(2)->SetTitle("#Delta#phi");
    fHistdPhidEtaMix->GetAxis(3)->SetTitle("#Delta#eta");
    fHistdPhidEtaMix->GetAxis(4)->SetTitle("p_{vz}");
    fHistdPhidEtaMix->GetAxis(5)->SetTitle("trigger");
    fHistdPhidEtaMix->GetAxis(6)->SetTitle("mass");
    fHistdPhidEtaMix->GetAxis(7)->SetTitle("multiplicity percentile");
    fHistdPhidEtaMix->GetAxis(8)->SetTitle("hadron charge");
    fHistdPhidEtaMix->Sumw2();
	fOutputList->Add(fHistdPhidEtaMix);
    fHistdPhidEtaMix->GetAxis(4)->Set(NofZVrtxBins,ZBins);
    fHistdPhidEtaMix->GetAxis(6)->Set(902,binsMass);
    
    fHistMCMixingRec = new THnSparseF ("fHistMCMixingRec", "fHistMCMixingRec", 9, bins, min, max);
    fHistMCMixingRec->GetAxis(0)->SetTitle("p_{T}^{trig}");
    fHistMCMixingRec->GetAxis(1)->SetTitle("p_{T}^{assoc}");
    fHistMCMixingRec->GetAxis(2)->SetTitle("#Delta#phi");
    fHistMCMixingRec->GetAxis(3)->SetTitle("#Delta#eta");
    fHistMCMixingRec->GetAxis(4)->SetTitle("p_{vz}");
    fHistMCMixingRec->GetAxis(5)->SetTitle("trigger");
    fHistMCMixingRec->GetAxis(6)->SetTitle("mass");
    fHistMCMixingRec->GetAxis(7)->SetTitle("multiplicity percentile");
    fHistMCMixingRec->GetAxis(8)->SetTitle("hadron charge");
    fOutputList->Add(fHistMCMixingRec);
    fHistMCMixingRec->Sumw2();
    fHistMCMixingRec->GetAxis(4)->Set(NofZVrtxBins,ZBins);
    fHistMCMixingRec->GetAxis(6)->Set(902,binsMass);
    
    fHistKorelacieMCrec = new THnSparseF ("fHistKorelacieMCrec","fHistKorelacieMCrec", 9, bins, min, max);
    fHistKorelacieMCrec->GetAxis(0)->SetTitle("p_{T}^{trig}");
    fHistKorelacieMCrec->GetAxis(1)->SetTitle("p_{T}^{assoc}");
    fHistKorelacieMCrec->GetAxis(2)->SetTitle("#Delta#phi");
    fHistKorelacieMCrec->GetAxis(3)->SetTitle("#Delta#eta");
    fHistKorelacieMCrec->GetAxis(4)->SetTitle("p_{vz}");
    fHistKorelacieMCrec->GetAxis(5)->SetTitle("trigger");
    fHistKorelacieMCrec->GetAxis(6)->SetTitle("mass");
    fHistKorelacieMCrec->GetAxis(7)->SetTitle("multiplicity percentile");
    fHistKorelacieMCrec->GetAxis(8)->SetTitle("hadron charge");
    fHistKorelacieMCrec->Sumw2();
    fOutputList->Add(fHistKorelacieMCrec);
    fHistKorelacieMCrec->GetAxis(6)->Set(902,binsMass);
    fHistKorelacieMCrec->GetAxis(4)->Set(NofZVrtxBins,ZBins);

    bins[7] = 500;
    max[7] = 500;

    fHistMCKorelacie = new THnSparseF ("fHistMCKorelacie","fHistMCKorelacie", 9, bins, min, max);
    fHistMCKorelacie->GetAxis(0)->SetTitle("p_{T}^{trig}");
    fHistMCKorelacie->GetAxis(1)->SetTitle("p_{T}^{assoc}");
    fHistMCKorelacie->GetAxis(2)->SetTitle("#Delta#phi");
    fHistMCKorelacie->GetAxis(3)->SetTitle("#Delta#eta");
    fHistMCKorelacie->GetAxis(4)->SetTitle("p_{vz}");
    fHistMCKorelacie->GetAxis(5)->SetTitle("trigger");
    fHistMCKorelacie->GetAxis(6)->SetTitle("mass");
    fHistMCKorelacie->GetAxis(7)->SetTitle("multiplicity");
    fHistMCKorelacie->GetAxis(8)->SetTitle("hadron charge");
    fHistMCKorelacie->Sumw2();
    fOutputList->Add(fHistMCKorelacie);
    
    fHistMCKorelacie->GetAxis(6)->Set(902,binsMass);
    fHistMCKorelacie->GetAxis(4)->Set(NofZVrtxBins,ZBins);

    fHistMCMixingGen = new THnSparseF ("fHistMCMixingGen", "fHistMCMixingGen", 9, bins, min, max);
    fHistMCMixingGen->GetAxis(0)->SetTitle("p_{T}^{trig}");
    fHistMCMixingGen->GetAxis(1)->SetTitle("p_{T}^{assoc}");
    fHistMCMixingGen->GetAxis(2)->SetTitle("#Delta#phi");
    fHistMCMixingGen->GetAxis(3)->SetTitle("#Delta#eta");
    fHistMCMixingGen->GetAxis(4)->SetTitle("p_{vz}");
    fHistMCMixingGen->GetAxis(5)->SetTitle("trigger");
    fHistMCMixingGen->GetAxis(6)->SetTitle("mass");
    fHistMCMixingGen->GetAxis(7)->SetTitle("multiplicity percentile");
    fHistMCMixingGen->GetAxis(8)->SetTitle("hadron charge");
    fOutputList->Add(fHistMCMixingGen);
    fHistMCMixingGen->Sumw2();
    fHistMCMixingGen->GetAxis(4)->Set(NofZVrtxBins,ZBins);
    fHistMCMixingGen->GetAxis(6)->Set(902,binsMass);

    fHistGenMultiplicity = new TH1D ("fHistGenMultiplicity","fHistGenMultiplicity",500,0,500);
    fOutputList->Add(fHistGenMultiplicity);

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

    Int_t binsHadr[5] ={fNumberOfPtBinsAssoc,9,fNumberOfEtaBins,2,fNumberPhiBins};
    Double_t binsHadrMin[5] = {fPtAsocMin,-10,-0.8,-2,0};
    Double_t binsHadrMax[5] = {fPtAssocMax,10,0.8,2,2*kPi};

	fHistMCPtAs = new THnSparseF("fHistMCPtAs","fHistMCPtAs",5,binsHadr,binsHadrMin,binsHadrMax);
	fOutputList->Add(fHistMCPtAs);
    fHistMCPtAs->Sumw2();
    fHistMCPtAs->GetAxis(1)->Set(NofZVrtxBins,ZBins);
    fHistMCPtAs->GetAxis(0)->SetTitle("p_{T}");
    fHistMCPtAs->GetAxis(1)->SetTitle("p_{vz}");
    fHistMCPtAs->GetAxis(2)->SetTitle("#eta");
    fHistMCPtAs->GetAxis(3)->SetTitle("charge");
    fHistMCPtAs->GetAxis(4)->SetTitle("#varphi");

	fHistRCPtAs = new THnSparseF("fHistRCPtAs","fHistRCPtAs",5,binsHadr,binsHadrMin,binsHadrMax);
    fOutputList->Add(fHistRCPtAs);
    fHistRCPtAs->Sumw2();
    fHistRCPtAs->GetAxis(1)->Set(NofZVrtxBins,ZBins);
    fHistRCPtAs->GetAxis(0)->SetTitle("p_{T}");
    fHistRCPtAs->GetAxis(1)->SetTitle("p_{vz}");
    fHistRCPtAs->GetAxis(2)->SetTitle("#eta");
    fHistRCPtAs->GetAxis(3)->SetTitle("charge");
    fHistRCPtAs->GetAxis(4)->SetTitle("#varphi");
    
    Int_t binsTrig[5]={fNumberOfPtBinsAssoc,9,6,fNumberOfEtaBins,fNumberPhiBins};
    Double_t mintrig[5]={fPtAsocMin,-10,0,-0.8,0};
    Double_t maxtrig[5]={fPtAssocMax,10,6,0.8,2*kPi};
    fHistGenV0 = new THnSparseF("fHistGenV0","fHistGenV0",5,binsTrig,mintrig,maxtrig);
    fOutputList->Add(fHistGenV0);
    fHistGenV0->Sumw2();
    fHistGenV0->GetAxis(0)->SetTitle("p_{T}");
    fHistGenV0->GetAxis(1)->SetTitle("p_{vz}");
    fHistGenV0->GetAxis(2)->SetTitle("trigger");
    fHistGenV0->GetAxis(3)->SetTitle("#eta");
    fHistGenV0->GetAxis(4)->SetTitle("#varphi");
    fHistGenV0->GetAxis(1)->Set(NofZVrtxBins,ZBins);
    Int_t binsTrigRec[6]={fNumberOfPtBinsAssoc,9,6,fNumberOfEtaBins,902,fNumberPhiBins};
    Double_t mintrigRec[6]={fPtAsocMin,-10,0,-0.8,0.44,0};
    Double_t maxtrigRec[6]={fPtAssocMax,10,6,0.8,1.355,2*kPi};
    fHistRecV0 = new THnSparseF("fHistRecV0","fHistRecV0",6,binsTrigRec,mintrigRec,maxtrigRec);
    fOutputList->Add(fHistRecV0);
    fHistRecV0->Sumw2();
    fHistRecV0->GetAxis(0)->SetTitle("p_{T}");
    fHistRecV0->GetAxis(1)->SetTitle("p_{vz}");
    fHistRecV0->GetAxis(2)->SetTitle("trigger");
    fHistRecV0->GetAxis(3)->SetTitle("#eta");
    fHistRecV0->GetAxis(4)->SetTitle("mass");
    fHistRecV0->GetAxis(5)->SetTitle("#varphi");
    fHistRecV0->GetAxis(4)->Set(902,binsMass);
    fHistRecV0->GetAxis(1)->Set(NofZVrtxBins,ZBins);

	fHistNumberOfTriggers = new THnSparseF("fHistNumberOfTriggers","fHistNumberOfTriggers",5,bins2d,mis2d, maxs2d);
    fHistNumberOfTriggers->GetAxis(0)->SetTitle("p_{T}");
    fHistNumberOfTriggers->GetAxis(1)->SetTitle("p_{vz}");
    fHistNumberOfTriggers->GetAxis(2)->SetTitle("trigger");
    fHistNumberOfTriggers->GetAxis(3)->SetTitle("mass");
    fHistNumberOfTriggers->GetAxis(4)->SetTitle("multiplicity percentile");
	fOutputList->Add(fHistNumberOfTriggers);
    fHistNumberOfTriggers->Sumw2();
    fHistNumberOfTriggers->GetAxis(3)->Set(902,binsMass);
    fHistNumberOfTriggers->GetAxis(1)->Set(NofZVrtxBins,ZBins);
    if(fNMultiplicityBins==16)fHistNumberOfTriggers->GetAxis(4)->Set(16,binsMult);

    fHistNumberOfTriggersRec = new THnSparseF("fHistNumberOfTriggersRec","fHistNumberOfTriggersRec",5,bins2d,mis2d,maxs2d);
    fHistNumberOfTriggersRec->GetAxis(0)->SetTitle("p_{T}");
    fHistNumberOfTriggersRec->GetAxis(1)->SetTitle("p_{vz}");
    fHistNumberOfTriggersRec->GetAxis(2)->SetTitle("trigger");
    fHistNumberOfTriggersRec->GetAxis(3)->SetTitle("mass");
    fHistNumberOfTriggersRec->GetAxis(4)->SetTitle("multiplicity percentile");
    fOutputList->Add(fHistNumberOfTriggersRec);
    fHistNumberOfTriggersRec->Sumw2();
    fHistNumberOfTriggersRec->GetAxis(3)->Set(902,binsMass);
    fHistNumberOfTriggersRec->GetAxis(1)->Set(NofZVrtxBins,ZBins);

    bins2d[4] = 500;
    maxs2d[4] = 500;

    fHistNumberOfTriggersGen = new THnSparseF("fHistNumberOfTriggersGen","fHistNumberOfTriggersGen",5,bins2d,mis2d, maxs2d);
    fHistNumberOfTriggersGen->GetAxis(0)->SetTitle("p_{T}");
    fHistNumberOfTriggersGen->GetAxis(1)->SetTitle("p_{vz}");
    fHistNumberOfTriggersGen->GetAxis(2)->SetTitle("trigger");
    fHistNumberOfTriggersGen->GetAxis(3)->SetTitle("mass");
    fHistNumberOfTriggersGen->GetAxis(4)->SetTitle("multiplicity");
    fOutputList->Add(fHistNumberOfTriggersGen);
    fHistNumberOfTriggersGen->Sumw2();
    fHistNumberOfTriggersGen->GetAxis(3)->Set(902,binsMass);
    fHistNumberOfTriggersGen->GetAxis(1)->Set(NofZVrtxBins,ZBins);
    

    fHistSelection = new TH1D("fHistSelection","fHistSelection",4,0,4);
    fOutputList->Add(fHistSelection);
    
    fHistMultipPercentile = new TH1F("fHistMultipPercentile","fHistMultipPercentile",10,0,100);
    fOutputList->Add(fHistMultipPercentile);
    fHistMultipPercentile->Sumw2();
    
    Int_t binsCuts[11] = {12,902,100,100,20,100,500,200,100,3,2};
    Double_t binsMinCuts[11] = {fPtTrigMin,0.44,0.03,0.03,0,0,0.95,0.,0,0,0};
    Double_t binsMaxCuts[11] = {15,1.355,0.25,0.25,2,1.5,1,50.,0.09,3,2};
    
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
    fHistTopolCut->GetAxis(1)->Set(902,binsMass);
    
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
    fHistTopolCutMC->GetAxis(1)->Set(902,binsMass);
    
    Int_t binsPur[6] = {fNumberOfPtBinsAssoc,903,4,8,8,40};
    Double_t binsPurMin[6] = {fPtAsocMin,0.44,0.,0.,0,-0.8};
    Double_t binsPurMax[6] = {fPtAssocMax,1.15,4.,8.,8.,0.8};
    fHistPurityCheck = new THnSparseF("fHistPurityCheck","fHistPurityCheck",6,binsPur,binsPurMin,binsPurMax);
    fOutputList->Add(fHistPurityCheck);
    fHistPurityCheck->Sumw2();
    fHistPurityCheck->GetAxis(0)->SetTitle("p_{T}");
    fHistPurityCheck->GetAxis(1)->SetTitle("mass");
    fHistPurityCheck->GetAxis(2)->SetTitle("V0 type");
    fHistPurityCheck->GetAxis(3)->SetTitle("check");
    fHistPurityCheck->GetAxis(5)->SetTitle("#eta");
    fHistPurityCheck->GetAxis(1)->Set(902,binsMass);
    
    fHistPtResolution = new TH3F("fHistPtResol","fHistPtResol",144,0,18,144,0,18,4,0,4);
    fOutputList->Add(fHistPtResolution);
    
    fHitsNTracks = new TH2F("fHitsNTracks","fHitsNTracks",100,0,100,2,0,2);
    fOutputList->Add(fHitsNTracks);
    
    Int_t binsPhiEta[4]= {fNumberOfPtBinsAssoc,72,40,4};
    Double_t minsPhiEta[4] = {fPtAsocMin,0,-0.8,0};
    Double_t maxsphiEta[4] = {fPtAssocMax,2*kPi,0.8,4};
    
    fHistPhiEta= new THnSparseF("fHistPhiEta","fHistPhiEta",4,binsPhiEta,minsPhiEta,maxsphiEta);
    fOutputList->Add(fHistPhiEta);

    fHistTPCTracksVsClusters =  new TH2D("fHistTPCTracksVsClusters","fHistTPCTracksVsClusters",400,0,400,50000,0,50000); 
    fOutputList->Add(fHistTPCTracksVsClusters);
    fHistTPCTracksVsClusters->Sumw2();
    fHistTPCTracksVsClusters->GetXaxis()->SetTitle("TPC tacks");
    fHistTPCTracksVsClusters->GetYaxis()->SetTitle("TPC Clusters");

    fHistVZeroPercentileTPCMult = new TH2D("fHistVZeroPercentileTPCMult","fHistVZeroPercentileTPCMult",200,0,100,800,0,800); 
    fOutputList->Add(fHistVZeroPercentileTPCMult);
    fHistVZeroPercentileTPCMult->Sumw2();
    fHistVZeroPercentileTPCMult->GetXaxis()->SetTitle("VZERO percentile");
    fHistVZeroPercentileTPCMult->GetYaxis()->SetTitle("TPC multiplicity");

    fHistPosNegTracks = new TH2D("fHistPosNegTracks","fHistPosNegTracks",fNumberOfPtBinsAssoc,fPtAsocMin,fPtAssocMax,3,0,3);
    fHistPosNegTracks->Sumw2();
    fOutputList->Add(fHistPosNegTracks);
    fHistPosNegTracks->GetXaxis()->SetTitle("p_{T}");
    fHistPosNegTracks->GetYaxis()->SetTitle("charge");

    fHistLambdaFeedDown = new TH3F("fHistLambdaFeedDown","fHistLambdaFeedDown",fNumberOfPtBinsAssoc,fPtAsocMin,fPtAssocMax,fNumberOfPtBinsAssoc,fPtAsocMin,fPtAssocMax,2,0,2);
    fHistLambdaFeedDown->Sumw2();
    fOutputList->Add(fHistLambdaFeedDown);
    fHistLambdaFeedDown->GetXaxis()->SetTitle("V0 p_{T}");
    fHistLambdaFeedDown->GetYaxis()->SetTitle("cascade p_{T}");
    fHistLambdaFeedDown->GetZaxis()->SetTitle("V0 type");

    fHistMultVZEROTracklets= new TH3D("fHistMultVZEROTracklets","fHistMultVZEROTracklets",1000,0,1000,100,0,100,16,0,100);
    fHistMultVZEROTracklets->Sumw2();
    fOutputList->Add(fHistMultVZEROTracklets);
    fHistMultVZEROTracklets->GetXaxis()->SetTitle("V0 Multiplicity");
    fHistMultVZEROTracklets->GetYaxis()->SetTitle("N rracklets");
    fHistMultVZEROTracklets->GetZaxis()->SetTitle("V0 percentile");
    fHistMultVZEROTracklets->GetZaxis()->Set(16,binsMult);

    if(!fAnalysisAOD){
        if(fFilterBit==32) fESDTrackCuts = AliESDtrackCuts:: GetStandardITSTPCTrackCuts2011(kTRUE);
        if(fFilterBit==16||fFilterBit==256) fESDTrackCuts = AliESDtrackCuts:: GetStandardITSTPCTrackCuts2011(kFALSE);
        if(fFilterBit==128) fESDTrackCuts = AliESDtrackCuts:: GetStandardTPCOnlyTrackCuts();
        if(fFilterBit==1){
            fESDTrackCuts = new AliESDtrackCuts;
            fESDTrackCuts->SetMaxChi2PerClusterTPC(4);
            fESDTrackCuts->SetAcceptKinkDaughters(kFALSE);
            fESDTrackCuts->SetDCAToVertex2D(kFALSE);
            
            fESDTrackCuts->SetMaxDCAToVertexZ(fMaxDCAToVertexZ);
            fESDTrackCuts->SetMaxDCAToVertexXY(fMaxDCAToVertexXY);
            fESDTrackCuts->SetRequireTPCRefit(kTRUE);
            fESDTrackCuts->SetRequireITSRefit(kTRUE);
            fESDTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSDD, AliESDtrackCuts::kAny);   
            fESDTrackCuts->SetRequireSigmaToVertex(kFALSE);

            fESDTrackCuts->SetMinNCrossedRowsTPC(fMinNCrossedRowsTPCprimtracks);
            fESDTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
        }
    } 
    fmcTracksSel = new TObjArray();
    fmcTracksSel->SetOwner(kTRUE);
    fmcGenTracksMixing = new TObjArray();
    fmcGenTracksMixing->SetOwner(kTRUE);
    fmcTracksTrigSel = new TObjArray();
    fmcTracksTrigSel->SetOwner(kTRUE);
    fmcTracksV0Sel = new TObjArray();
    fmcTracksV0Sel->SetOwner(kTRUE);
    fmcV0AssocSel = new TObjArray();
    fmcV0AssocSel->SetOwner(kTRUE);
    fselectedMCassoc = new TObjArray();
    fselectedMCassoc->SetOwner(kTRUE);
    fselectedMCV0assoc = new TObjArray();
    fselectedMCV0assoc->SetOwner(kTRUE);
    fselectedMCtrig = new TObjArray();
    fselectedMCtrig->SetOwner(kTRUE);
    fselectedMCV0Triggersrec = new TObjArray();
    fselectedMCV0Triggersrec->SetOwner(kTRUE);
    fselectedTracks = new TObjArray();
    fselectedAssociatedTracks = new TObjArray();
    fselectedAssociatedTracks->SetOwner(kTRUE);
    fselectedTriggerTracks = new TObjArray();
    fselectedTriggerTracks->SetOwner(kTRUE);
    fselectedV0Triggers = new TObjArray();
    fselectedV0Triggers->SetOwner(kTRUE);
    fselectedV0Assoc = new TObjArray();
    fselectedV0Assoc->SetOwner(kTRUE);  

    fAliEventCuts = new AliEventCuts();
    if(fEventCutsQAPlots)fAliEventCuts->AddQAplotsToList(fOutputList,kTRUE);
    if(fpp)fAliEventCuts->SetupRun2pp();
    else fAliEventCuts->SetupRun2PbPb();
    if(fhighMult) fAliEventCuts->OverrideAutomaticTriggerSelection(AliVEvent::kHighMultV0,kTRUE);
    if(fhighMultSPD) fAliEventCuts->OverrideAutomaticTriggerSelection(AliVEvent::kHighMultSPD,kTRUE);

    PostData(1, fOutputList);   // postdata will notify the analysis manager of changes / updates to the 
                                        // fOutputList object. the manager will in the end take care of writing your output to file
                                        // so it needs to know what's in the output
	// Settings for event mixing -------------------------------------
    Int_t trackDepth = fMixingTracks;
    Int_t poolSize   = 100;  // Maximum number of events, ignored in the present implemented of AliEventPoolManager
  
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
    // once you return from the UserExec function, the manager will retrieve the next event from the chain#

    Double_t lPercentile = 302;

    Int_t iTracks = 0;
    Int_t nV0 =0;
    Int_t nCascades =0;
    
    AliAODVertex *myPrimVertex = 0x0;
    AliESDVertex *myPrimVertexESD = 0x0;

    if(fCorrelations||fEfficiency||fMixing){
        AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
        AliAODInputHandler *inEvMain = (AliAODInputHandler *) mgr->GetInputEventHandler();

        if(fAnalysisAOD) {
            fAOD = (AliAODEvent *) inEvMain -> GetEvent();    // get an event (called fAOD) from the input file
            if(!fAOD) {
                cout << "ERROR: no AOD event" << endl;
                return;                                 // if the pointer to the ,event is empty (getting it failed) skip this event
            }
        }else{
            fESD = (AliESDEvent *) inEvMain -> GetEvent();
            if(!fESD) {
                cout << "ERROR: no ESD event" << endl;
                return;
            }
        }                                   
    
        fPIDResponse = (AliPIDResponse *) inEvMain-> GetPIDResponse();
    
        // physics selection
        fHistSelection->Fill(0.5);
	   UInt_t maskIsSelected = inEvMain->IsEventSelected();

	   //  data trigger selection
	   Bool_t isSelected = kFALSE;

       if(fminBias) isSelected = ((maskIsSelected & AliVEvent::kINT7) == AliVEvent::kINT7); //pp - minumum bias
       if(fhighMult) isSelected = ((maskIsSelected & AliVEvent::kHighMultV0)== AliVEvent::kHighMultV0); //pp - high multiplicity V0
       if(fhighMultSPD)isSelected = ((maskIsSelected & AliVEvent::kHighMultSPD) == AliVEvent::kHighMultSPD); //pp - high multiplicity V0
          
	   if (!isSelected) return;
        fHistSelection->Fill(1.5);

	   if(fAOD){
            myPrimVertex = (AliAODVertex *)fAOD->GetPrimaryVertex();
            if (!myPrimVertex) return;
            fPV[2] = myPrimVertex->GetZ();
            fPV[0] = myPrimVertex->GetX();
            fPV[1] = myPrimVertex->GetY();
        }
        if(fESD){
            myPrimVertexESD = (AliESDVertex *)fESD->GetPrimaryVertex();
            if (!myPrimVertexESD) return;
            fPV[2]= myPrimVertexESD->GetZ();
            fPV[0] = myPrimVertexESD->GetX();
            fPV[1] = myPrimVertexESD->GetY();
            fMagneticField = fESD->GetMagneticField( );
        }
        

        if (TMath::Abs(fPV[2])>=fPrimaryVertexCut) return;
        fHistSelection->Fill(2.5);

        if(!fAnalysisMC&&fRejectEventPileUp){
            if(fAOD) if(!fAliEventCuts->AcceptEvent(fAOD)) return;
            if(fESD) if(!fAliEventCuts->AcceptEvent(fESD)) {
                PostData(1, fOutputList); 
                return;
            }
        }
        fHistSelection->Fill(3.5);

        Int_t tpcMult = 0;
        Int_t tpcClusters = 0;
        Float_t VZEROmultiplicity =0;
        Int_t nTracklets=0;
    
        if(fAOD){
            tpcMult = fAOD->GetNumberOfTPCTracks();
            tpcClusters = fAOD->GetNumberOfTPCClusters();
            iTracks = (fAOD->GetNumberOfTracks());           // see how many tracks there are in the event
            nV0 = (fAOD->GetNumberOfV0s());                  // see how many V0 there are in the event
            nCascades = fAOD->GetNumberOfCascades();
        }
        if(fESD){
            tpcMult = fESD->GetNumberOfTPCTracks();
            iTracks = (fESD->GetNumberOfTracks());
            tpcClusters = fESD->GetNumberOfTPCClusters();
            nV0 = (fESD->GetNumberOfV0s()); 
            nCascades = fESD->GetNumberOfCascades();
            AliESDVZERO *vZERO = fESD->GetVZEROData();
            VZEROmultiplicity = vZERO->GetMTotV0A();
            VZEROmultiplicity += vZERO->GetMTotV0C();
            nTracklets=fESD->GetMultiplicity()->GetNumberOfTracklets();
        }

	   fHistV0Multiplicity->Fill(nV0);
       fHistTPCTracksVsClusters->Fill(tpcMult,tpcClusters);

	   // Multiplicity definition
        AliMultSelection *MultSelection = 0x0;
        if(fAOD) MultSelection = (AliMultSelection * ) fAOD->FindListObject("MultSelection");
        if(fESD) MultSelection = (AliMultSelection * ) fESD->FindListObject("MultSelection");
        if( !MultSelection) {
            AliWarning("AliMultSelection object not found!");
        }else{
            if(fMultEstimator=="V0M") lPercentile = MultSelection->GetMultiplicityPercentile("V0M");
            else if(fMultEstimator=="V0A") lPercentile = MultSelection->GetMultiplicityPercentile("V0A");
            else if(fMultEstimator=="SPDTracklets") lPercentile = MultSelection->GetMultiplicityPercentile("SPDTracklets");
            else if(fMultEstimator=="RefMult05") lPercentile = MultSelection->GetMultiplicityPercentile("RefMult05");
            else if(fMultEstimator=="RefMult08") lPercentile = MultSelection->GetMultiplicityPercentile("RefMult08");
        }
        if ((lPercentile<fPercentileMin)||(lPercentile>fPercetileMax)) return;
    
        fHistMultVZEROTracklets->Fill(VZEROmultiplicity,nTracklets,lPercentile);
        fHistMultipPercentile->Fill(lPercentile);
        fHistVZeroPercentileTPCMult -> Fill(lPercentile,tpcMult);
    }

	//=========== MC loop ===============================
    Int_t nAcceptedParticles =0;
    AliMCParticle *mcTrack = 0x0;

	if(fAnalysisMC){
        fmcEvent  = dynamic_cast<AliMCEvent*> (MCEvent());
        if(!fmcEvent){
            Printf("No MC particle branch found");
            return;
        }

        AliVVertex * mcVertex = (AliVVertex* ) fmcEvent->GetPrimaryVertex();
        fPV[2] = mcVertex->GetZ();
        if (TMath::Abs(fPV[2])>=fPrimaryVertexCut) return;

        Int_t nMCAllTracks = fmcEvent->GetNumberOfTracks();

        AliVTrack *genTrackMix = 0x0;

        for (Int_t i = 0; i < nMCAllTracks; i++){
            AliMCParticle *mcTrack = (AliMCParticle*)fmcEvent->GetTrack(i);
            if (!mcTrack) {
                Error("ReadEventAODMC", "Could not receive particle %d", i);
                continue;
            }

            Double_t trackPseudorap = mcTrack->Eta();

            if( mcTrack->IsPhysicalPrimary()&&mcTrack->Charge()!=0&&((trackPseudorap>-3.7&&trackPseudorap<-1.7)||(trackPseudorap>2.8&&trackPseudorap<5.1))){
                nAcceptedParticles += 1;
            }
        }

        fHistGenMultiplicity->Fill(nAcceptedParticles);

        AliMCParticle *mcMotherParticle = 0x0;
        AliMCParticle* daughter0 = 0x0;
        AliMCParticle* daughter1 = 0x0;
	
		for (Int_t i = 0; i < nMCAllTracks; i++){
			mcTrack = (AliMCParticle*)fmcEvent->GetTrack(i);
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
            Double_t cha;
            if (mcTrack->Charge()>0) cha=1.;
            else if (mcTrack->Charge()<0) cha= -1.;
            else cha =0;

			if (TrIsPrim && TrPtMin && TrCharge && TrEtaMax) {
                
                if(fEfficiency) {
                    Double_t eff[5] = {mcTrackPt,fPV[2],mcTrackEta, cha,mcTrack->Phi()};
                    fHistMCPtAs->Fill(eff); // for recunstruction efficiency calculation
                }
                if(fCorrelationsGen) fmcTracksSel->Add(new AliV0ChParticle(mcTrack->Eta(),mcTrack->Phi(),mcTrack->Pt(),4,mcTrack->GetLabel(),mcTrack->GetLabel()));
                
                if(fMixingGen){
                    genTrackMix = SetAliAODTrack(mcTrack->Theta(),mcTrack->Phi(),mcTrack->Pt(),mcTrack->Charge());
                    fmcGenTracksMixing->Add(genTrackMix);
                    lPercentile = 19;
                }
                
                if (mcTrackPt>fPtTrigMin) {
                    if(fMixingGen||fCorrelationsGen) fmcTracksTrigSel->Add(new AliV0ChParticle(mcTrack->Eta(),mcTrack->Phi(),mcTrack->Pt(),4,mcTrack->GetLabel(),mcTrack->GetLabel()));
                }
            }
            //--- MC closure test - selection of V0 ----

            Int_t mcPartPdg = mcTrack->PdgCode();
            Bool_t isPhysPrim = mcTrack->IsPhysicalPrimary();
            Double_t V0genrapidity = mcTrack->Y();

            Bool_t IsK0, IsLambda, IsAntiLambda, IsPositiveXi, IsNegativeXi;
            Double_t etaDau0 =0.;
            Double_t etaDau1 = 0.;

            Int_t labelPos = -1;
            Int_t labelNeg = -1;

            if(fTestPions){
            	
            	if ((mcPartPdg != 111) && (mcPartPdg != 211) && (mcPartPdg != (-211))) continue; // keep only pions 

            	IsK0 = mcPartPdg==111&& (isPhysPrim);
            	IsLambda = mcPartPdg==211&& (isPhysPrim);
            	IsAntiLambda = mcPartPdg==-211&& (isPhysPrim);

            }else {

            	if ((mcPartPdg != 310) && (mcPartPdg != 3122) && (mcPartPdg != (-3122)) && TMath::Abs(mcPartPdg) !=3312 ) continue; // keep only Lambdas and K0S and charged Xi (for the feedDownCorrection)
            	Bool_t IsFromCascade = kFALSE;
            	
                if(fAnalysisAOD){
                	Int_t mother  = mcTrack->GetMother();
                	mcMotherParticle = static_cast<AliMCParticle*>(fmcEvent->GetTrack(mother));
                	Int_t motherPDG = 0;
                	if (mother<0) motherPDG =0;
                	else motherPDG = TMath::Abs(mcMotherParticle->PdgCode());
                
               		if(fAacceptLambdasFromCasscade) IsFromCascade = (((motherPDG == 3222)|| (motherPDG==3212)|| (motherPDG==3112) || (motherPDG==3224) || (motherPDG==3214) || (motherPDG==3114) || (motherPDG==3322) || (motherPDG==3312)|| (motherPDG==3324) || (motherPDG==3314) || (motherPDG==3334)) && (mcMotherParticle->IsPhysicalPrimary()));

                	Int_t dau0 = mcTrack->GetDaughterLabel(0);
                	if (dau0>0) daughter0 = (AliMCParticle*) fmcEvent->GetTrack(dau0);
                	Int_t dau1 = mcTrack->GetDaughterLabel(1);
                	if (dau1>0) daughter1 = (AliMCParticle*) fmcEvent->GetTrack(dau1);
                
                	if(!daughter0||!daughter1) continue;

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

                	etaDau0 = daughter0->Eta();
                	etaDau1 = daughter1->Eta();
                }
                
                IsK0 = mcPartPdg==310&& (isPhysPrim);
                IsLambda = mcPartPdg==3122&& (isPhysPrim||IsFromCascade);
                IsAntiLambda = mcPartPdg==-3122&& (isPhysPrim||IsFromCascade);
            }
            IsPositiveXi = mcPartPdg == -3312&& (isPhysPrim);
            IsNegativeXi = mcPartPdg == 3312&& (isPhysPrim);
            
            
            if (mcTrack->Pt()>fPtAsocMin&&TMath::Abs(V0genrapidity)<0.5&&TMath::Abs(etaDau0)<0.8&&TMath::Abs(etaDau1)<0.8){
                if(IsK0) {
                    if(fMixingGen||fCorrelationsGen) fmcV0AssocSel->Add(new AliV0ChParticle(mcTrack->Eta(),mcTrack->Phi(),mcTrack->Pt(),5,mcTrack->GetLabel(),labelPos,labelNeg,kFALSE,mcTrack->M()));
                    if (fEfficiency){
                        Double_t v0effic[5]={mcTrack->Pt(),fPV[2],0.5,mcTrack->Eta(),mcTrack->Phi()};
                        fHistGenV0->Fill(v0effic); // for recunstruction efficiency calculation
                    }
                }
                if(IsLambda) {
                    if(fMixingGen||fCorrelationsGen) fmcV0AssocSel->Add(new AliV0ChParticle(mcTrack->Eta(),mcTrack->Phi(),mcTrack->Pt(),6,mcTrack->GetLabel(),labelPos,labelNeg,kFALSE,mcTrack->M()));
                    if (fEfficiency){
                        Double_t v0effic[5]={mcTrack->Pt(),fPV[2],1.5,mcTrack->Eta(),mcTrack->Phi()};
                        fHistGenV0->Fill(v0effic); // for recunstruction efficiency calculation
                    }
                }
                if(IsAntiLambda) {
                    if(fMixingGen||fCorrelationsGen) fmcV0AssocSel->Add(new AliV0ChParticle(mcTrack->Eta(),mcTrack->Phi(),mcTrack->Pt(),7,mcTrack->GetLabel(),labelPos,labelNeg,kFALSE,mcTrack->M()));
                    if (fEfficiency){
                        Double_t v0effic[5]={mcTrack->Pt(),fPV[2],2.5,mcTrack->Eta(),mcTrack->Phi()};
                        fHistGenV0->Fill(v0effic); // for recunstruction efficiency calculation
                    }
                }
                
            }
            if (mcTrack->Pt()>fPtTrigMin&&(fMixingGen||fCorrelationsGen)&&TMath::Abs(V0genrapidity)<0.5&&TMath::Abs(etaDau0)<0.8&&TMath::Abs(etaDau1)<0.8){
                if(IsK0) fmcTracksV0Sel->Add(new AliV0ChParticle(mcTrack->Eta(),mcTrack->Phi(),mcTrack->Pt(),1,mcTrack->GetLabel(),labelPos,labelNeg,kFALSE,mcTrack->M()));
                if(IsLambda) fmcTracksV0Sel->Add(new AliV0ChParticle(mcTrack->Eta(),mcTrack->Phi(),mcTrack->Pt(),2,mcTrack->GetLabel(),labelPos,labelNeg,kFALSE,mcTrack->M()));
                if(IsAntiLambda) fmcTracksV0Sel->Add(new AliV0ChParticle(mcTrack->Eta(),mcTrack->Phi(),mcTrack->Pt(),3,mcTrack->GetLabel(),labelPos,labelNeg,kFALSE,mcTrack->M()));
            }

            if (mcTrack->Pt()>fPtAsocMin&&TMath::Abs(V0genrapidity)<0.5){
                if(IsPositiveXi){
                    if (fCorrelationsGen) fmcV0AssocSel->Add(new AliV0ChParticle(mcTrack->Eta(),mcTrack->Phi(),mcTrack->Pt(),11,mcTrack->M(),1,2,3));
                    if (fEfficiency){
                        Double_t v0effic[5]={mcTrack->Pt(),fPV[2],5.5,mcTrack->Eta(),mcTrack->Phi()};
                        fHistGenV0->Fill(v0effic); // for recunstruction efficiency calculation
                    }
                }
                if(IsNegativeXi) {
                    if (fCorrelationsGen) fmcV0AssocSel->Add(new AliV0ChParticle(mcTrack->Eta(),mcTrack->Phi(),mcTrack->Pt(),10,mcTrack->M(),1,2,3));
                    if (fEfficiency){
                        Double_t v0effic[5]={mcTrack->Pt(),fPV[2],4.5,mcTrack->Eta(),mcTrack->Phi()};
                        fHistGenV0->Fill(v0effic); // for recunstruction efficiency calculation
                    }
                }   
            }
            if (mcTrack->Pt()>fPtTrigMin&&(fCorrelationsGen)&&TMath::Abs(V0genrapidity)<0.5){
                if(IsNegativeXi) fmcTracksV0Sel->Add(new AliV0ChParticle(mcTrack->Eta(),mcTrack->Phi(),mcTrack->Pt(),8,mcTrack->M(),1,2,3));
                if(IsPositiveXi) fmcTracksV0Sel->Add(new AliV0ChParticle(mcTrack->Eta(),mcTrack->Phi(),mcTrack->Pt(),9,mcTrack->M(),1,2,3));
            }

		}
        // MC closure test corellations
        if(fCorrelationsGen){
            //V0-h
            if(fV0hCorr) Corelations(fmcTracksV0Sel,fmcTracksSel,kFALSE,kTRUE,nAcceptedParticles,kFALSE);

            //h-h
            if(fhhCorr) Corelations(fmcTracksTrigSel,fmcTracksSel, kFALSE, kFALSE,nAcceptedParticles,kFALSE);
            //h-V0
            if(fhV0Corr) Corelations(fmcTracksTrigSel,fmcV0AssocSel, kFALSE, kTRUE,nAcceptedParticles,kTRUE);
            
            if(fAnalyseFeedDown){
                CorrelationsXi(fmcTracksV0Sel,fmcTracksSel,nAcceptedParticles,kTRUE);
                CorrelationsXi(fmcTracksTrigSel,fmcV0AssocSel,nAcceptedParticles,kFALSE);
            }
        }
    }

	//=========== end of MC loop ==========

    Int_t nTrak = 0;
    Int_t nTrakBefore = 0;
    Double_t ptTrack =0.;
    Double_t EtaTrack =0.;
    Double_t phiTrack =0.;

    if(fCorrelations||fEfficiency||fMixing){
        for(Int_t i=0; i < iTracks; i++) {                 // loop over all these tracks
            if(fAOD){
                AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(i));         // get a track (type AliAODTrack) from the event
                if(!track) continue;                            // if we failed, skip this track

                if(!IsMyGoodPrimaryTrack(track)) continue; // hybrid track selection
                
                ptTrack = track->Pt();
                EtaTrack = track->Eta();
                phiTrack = track->Phi();
                if(ptTrack<fPtAsocMin) continue;
        
                nTrakBefore+=1;
                if(fRejectTrackPileUp&&fRejectTOF&&(!(track->HasPointOnITSLayer(1) || track->HasPointOnITSLayer(0) || track->GetTOFBunchCrossing()==0 ))) continue; // track by track pile-up rejection using TOF
                if(fRejectTrackPileUp&&!fRejectTOF&&(!(track->HasPointOnITSLayer(1) ||track->HasPointOnITSLayer(0)))) continue; // track by track pile-up rejection without TOF information
                nTrak+=1;

                if(ptTrack>fPtAsocMin) fselectedTracks->Add(track);

                fHistPosNegTracks->Fill(ptTrack,track->Charge()+1.5);
                if(fAnalysisMC){
                    Int_t AssocLabel = track->GetLabel();
                    if (AssocLabel<=0) continue;

                    Double_t purhadr[6] = {ptTrack,0,3.5,0.5,-1,EtaTrack};
                    fHistPurityCheck->Fill(purhadr);
                
                    mcTrack = static_cast<AliMCParticle*>(fmcEvent->GetTrack(AssocLabel));
                    if(!mcTrack) continue;
                    Bool_t isPhyPrim = mcTrack->IsPhysicalPrimary();
                    Double_t genPt = mcTrack->Pt();
                    Double_t genEta = mcTrack->Eta();
                    if((fMixing||fCorrelations)&&!fPurePrimHadrons) fselectedMCassoc->Add(new AliV0ChParticle(EtaTrack,phiTrack,ptTrack,4,AssocLabel,track->GetID(),track->Charge(),track->Pz(),track->E()));
                    Double_t phiEtaData[4] = {ptTrack,phiTrack,EtaTrack,3.5};
                    fHistPhiEta->Fill(phiEtaData);
                    if (ptTrack>fPtTrigMin) {
                        if((fMixing||fCorrelations)&&!fPurePrimHadrons) fselectedMCtrig->Add(new AliV0ChParticle(EtaTrack,phiTrack,ptTrack,4,AssocLabel,track->GetID(),track->Charge(),track->Pz(),track->E()));
                    }
                    
                    if (isPhyPrim) {
                        if(fCorrelations&&fPurePrimHadrons) fselectedMCassoc->Add(new AliV0ChParticle(EtaTrack,phiTrack,ptTrack,4,AssocLabel,track->GetID(),track->Charge(),track->Pz(),track->E()));
                        if (ptTrack>fPtTrigMin) {
                            if(fCorrelations&&fPurePrimHadrons) fselectedMCtrig->Add(new AliV0ChParticle(EtaTrack,phiTrack,ptTrack,4,AssocLabel,track->Charge(),track->GetID(),track->Pz(),track->E()));
                        }
                        Double_t purhadrPrim[6] = {ptTrack,0,3.5,1.5,-1,EtaTrack};
                        fHistPurityCheck->Fill(purhadrPrim);
                        fHistPtResolution->Fill(genPt,ptTrack,3.5);
                        if(fEfficiency) { 
                            Double_t cha =0.;
                            if (track->Charge()>0) cha=1.;
                            else if (track->Charge()<0) cha= -1.;
                            Double_t eff[5] = {genPt,fPV[2],genEta,cha,mcTrack->Phi()};
                            fHistRCPtAs->Fill(eff);
                        } // for recunstruction efficiency calculation
                    }
                    if(mcTrack->IsSecondaryFromMaterial()){
                        Double_t purhadrMater[6] = {ptTrack,0,3.5,2.5,0.5,EtaTrack};
                        fHistPurityCheck->Fill(purhadrMater);
                    }
                    if(mcTrack->IsSecondaryFromWeakDecay()){
                        Double_t purhadrDecay[6] = {ptTrack,0,3.5,2.5,1.5,EtaTrack};
                        fHistPurityCheck->Fill(purhadrDecay);
                    }

                }else{
                    if(ptTrack>fPtAsocMin) {
                        fselectedAssociatedTracks-> Add(new AliV0ChParticle(EtaTrack, phiTrack, ptTrack, 4, 0,track->GetID(),track->Charge(),track->Pz(),track->E()));
                        Double_t phiEtaData[4] = {ptTrack,phiTrack,EtaTrack,3.5};
                        fHistPhiEta->Fill(phiEtaData);
                    }
                    if(ptTrack>fPtTrigMin&&!fAnalysisMC) fselectedTriggerTracks-> Add(new AliV0ChParticle(EtaTrack, phiTrack, ptTrack, 4,0,track->GetID(),track->Charge(),track->Pz(),track->E()));
                }
            }
            if(fESD){
                
                AliESDtrack * track =static_cast<AliESDtrack*>(fESD->GetTrack(i)); 
                if(!track) continue; 
                if(!IsMyGoodPrimaryTrackESD(track)) continue; 

                ptTrack = track->Pt();
                EtaTrack = track->Eta();
                phiTrack = track->Phi();

                if(ptTrack<fPtAsocMin) continue;

                nTrakBefore+=1;
                if(fRejectTrackPileUp&&fRejectTOF&&(!(track->HasPointOnITSLayer(1) || track->HasPointOnITSLayer(0) || track->GetTOFBunchCrossing()==0 ))) continue; // track by track pile-up rejection using TOF
                if(fRejectTrackPileUp&&!fRejectTOF&&(!(track->HasPointOnITSLayer(1) ||track->HasPointOnITSLayer(0)))) continue; // track by track pile-up rejection without TOF information
                nTrak+=1;

                if(ptTrack>fPtAsocMin) fselectedTracks->Add(track);

                fHistPosNegTracks->Fill(ptTrack,track->Charge()+1.5);

                if(fAnalysisMC){
                    Int_t AssocLabel = track->GetLabel();
                    if (AssocLabel<=0) continue;

                    Double_t purhadr[6] = {ptTrack,0,3.5,0.5,-1,EtaTrack};
                    fHistPurityCheck->Fill(purhadr);
                
                    mcTrack = static_cast<AliMCParticle*>(fmcEvent->GetTrack(AssocLabel));
                    if(!mcTrack) continue;
                    Bool_t isPhyPrim = mcTrack->IsPhysicalPrimary();
                    Double_t genPt = mcTrack->Pt();
                    Double_t genEta = mcTrack->Eta();

                    if((fMixing||fCorrelations)&&!fPurePrimHadrons) fselectedMCassoc->Add(new AliV0ChParticle(EtaTrack,phiTrack,ptTrack,4,AssocLabel,track->GetID(),track->Charge(),track->Pz(),track->E()));
                    Double_t phiEtaData[4] = {ptTrack,phiTrack,EtaTrack,3.5};
                    fHistPhiEta->Fill(phiEtaData);
                    if (ptTrack>fPtTrigMin) {
                        if((fMixing||fCorrelations)&&!fPurePrimHadrons) fselectedMCtrig->Add(new AliV0ChParticle(EtaTrack,phiTrack,ptTrack,4,AssocLabel,track->GetID(),track->Charge(),track->Pz(),track->E()));
                    }
                    
                    if (isPhyPrim) {
                        if((fMixing||fCorrelations)&&fPurePrimHadrons) fselectedMCassoc->Add(new AliV0ChParticle(EtaTrack,phiTrack,ptTrack,4,AssocLabel,track->GetID(),track->Charge(),track->Pz(),track->E()));
                        if (ptTrack>fPtTrigMin) {
                            if((fMixing||fCorrelations)&&fPurePrimHadrons) fselectedMCtrig->Add(new AliV0ChParticle(EtaTrack,phiTrack,ptTrack,4,AssocLabel,track->Charge(),track->GetID(),track->Pz(),track->E()));
                        }
                        Double_t purhadrPrim[6] = {ptTrack,0,3.5,1.5,-1,EtaTrack};
                        fHistPurityCheck->Fill(purhadrPrim);
                        fHistPtResolution->Fill(genPt,ptTrack,3.5);
                        if(fEfficiency) {
                            Double_t cha =0.;
                            if (track->Charge()>0) cha=1.;
                            else if (track->Charge()<0) cha= -1.;
                            Double_t eff[5] = {genPt,fPV[2],genEta,cha,mcTrack->Phi()};
                            fHistRCPtAs->Fill(eff); // for recunstruction efficiency calculation
                        }
                    }
                    if(mcTrack->IsSecondaryFromMaterial()){
                        Double_t purhadrMater[6] = {ptTrack,0,3.5,2.5,0.5,EtaTrack};
                        fHistPurityCheck->Fill(purhadrMater);
                    }
                    if(mcTrack->IsSecondaryFromWeakDecay()){
                        Double_t purhadrDecay[6] = {ptTrack,0,3.5,2.5,1.5,EtaTrack};
                        fHistPurityCheck->Fill(purhadrDecay);
                    }
                }else{
                    if(ptTrack>fPtAsocMin) {
                        fselectedAssociatedTracks-> Add(new AliV0ChParticle(EtaTrack, phiTrack, ptTrack, 4, 0,track->GetID(),track->Charge(),track->Pz(),track->E()));
                        Double_t phiEtaData[4] = {ptTrack,phiTrack,EtaTrack,3.5};
                        fHistPhiEta->Fill(phiEtaData);
                    }
                    if(ptTrack>fPtTrigMin) fselectedTriggerTracks-> Add(new AliV0ChParticle(EtaTrack, phiTrack, ptTrack, 4,0,track->GetID(),track->Charge(),track->Pz(),track->E()));
                }
            } 	           
        }

        fHitsNTracks->Fill(nTrakBefore,0.5);
        fHitsNTracks->Fill(nTrak,1.5);
    

        AliAODTrack *myTrackPosAOD = 0x0;
        AliAODTrack *myTrackNegAOD = 0x0;
        AliESDtrack *myTrackPosESD = 0x0;
        AliESDtrack *myTrackNegESD = 0x0;

        AliAODTrack * myTrackNegTest = 0x0;
        AliAODTrack * myTrackPosTest = 0x0;

        Int_t nK0 =0;
        Int_t nLam =0;
        Int_t nALam =0;

        AliAODv0* V0AOD =  0x0;
        AliESDcascade* cascadeESD =  0x0;
        AliESDv0* V0esd =  0x0;
        AliVParticle *V0 = 0x0;
        Double_t v0pt =0;
        Double_t v0phi =0;
        Double_t v0Eta =0;
        Int_t idNeg =0; 
        Int_t idPos =0;
        Int_t labelPos =0; 
        Int_t labelNeg =0;

        AliAODTrack *daughter1 =0x0;
        AliAODTrack *daughter0 =0x0;

        Double_t massK0 = 0.; 
        Double_t massLambda =0.; 
        Double_t massAntilambda =0.; 
        if((fAnalyseFeedDown&&fCorrelations)||(fAnalyseFeedDown&&fEfficiency)||(fAnalyseFeedDown&&fMixing)){
            
            Int_t labelpTrackXi, labelnTrackXi,labelbTrackXi;   
            Int_t NegTrackPdg, PosTrackPdg, BacTrackPdg, V0Pdg,cascadePDG;
            Int_t V0label, cascadelabel;
            Double_t cascadept = 0.;
            Double_t par[5] = {6,1.3,1,2,3};

            for (int iC = 0; iC < nCascades; ++iC)
            {   
                if(fAOD) {
                    cout << "Error: No AOD cascade implementation" << endl;
                    continue;
                }

                cascadeESD = (AliESDcascade*)(fESD->GetCascade(iC));
                if(!cascadeESD) {
                    cout << "Error: No ESD cascade" << endl;
                    continue;
                }
                
                cascadept = cascadeESD->Pt();

                if(cascadept < fPtAsocMin) continue;
                    
                if(!IsMyGoodXiCandidate(cascadeESD,par)) continue;
                    
                if(cascadeESD->Pt()>fPtAsocMin) {
                    fselectedV0Assoc-> Add(new AliV0ChParticle(cascadeESD->Eta(), cascadeESD->Phi(), cascadept, (Int_t) par[0]+4,par[1],(Int_t) par[2],(Int_t) par[3],(Int_t) par[4]));
                } 
                if(cascadeESD->Pt()>fPtTrigMin) {
                    fselectedV0Triggers-> Add(new AliV0ChParticle(cascadeESD->Eta(), cascadeESD->Phi(), cascadept, (Int_t) par[0]+2,par[1],(Int_t) par[2],(Int_t) par[3],(Int_t) par[4]));
                }

                fHistCasMC->Fill(par[0]-0.5,0.5);

                if(fAnalysisMC){
                    fHistCasMC->Fill(par[0]-0.5,1.5);

                    UInt_t idxPosXi = (UInt_t) TMath::Abs( cascadeESD->GetPindex() );
                    UInt_t idxNegXi = (UInt_t) TMath::Abs( cascadeESD->GetNindex() );
                    UInt_t idBachXi = (UInt_t) TMath::Abs( cascadeESD->GetBindex() );

                    AliESDtrack * pTrackXi = fESD->GetTrack( idxPosXi );
                    AliESDtrack * nTrackXi = fESD->GetTrack( idxNegXi );
                    AliESDtrack * bachTrackXi = fESD->GetTrack( idBachXi );

                    labelpTrackXi = (Int_t) TMath::Abs(pTrackXi->GetLabel());
                    labelnTrackXi = (Int_t) TMath::Abs(nTrackXi->GetLabel());
                    labelbTrackXi = (Int_t) TMath::Abs(bachTrackXi->GetLabel());

                    AliMCParticle *mcPosTrack = (AliMCParticle*)fmcEvent->GetTrack(labelpTrackXi);
                    AliMCParticle *mcNegTrack = (AliMCParticle*)fmcEvent->GetTrack(labelnTrackXi);
                    AliMCParticle *mcBacTrack = (AliMCParticle*)fmcEvent->GetTrack(labelbTrackXi);
                    if (!mcNegTrack||!mcPosTrack||!mcBacTrack) continue;

                    fHistCasMC->Fill(par[0]-0.5,2.5);

                    NegTrackPdg = mcNegTrack->PdgCode();
                    PosTrackPdg = mcPosTrack->PdgCode();
                    BacTrackPdg = mcBacTrack->PdgCode();

                    V0label = IsGoodMCV0(mcPosTrack,mcNegTrack); 

                    if(V0label < 0) continue;
        
                    fHistCasMC->Fill(par[0]-0.5,4.5);

                    AliMCParticle *mcPosMother = (AliMCParticle*)fmcEvent->GetTrack(V0label);
                    if (!mcPosMother) continue;
                    fHistCasMC->Fill(par[0]-0.5,5.5);

                    V0Pdg = mcPosMother->PdgCode();

                    cascadelabel = IsGoodMCCascade(mcPosMother,mcBacTrack);
                    if(cascadelabel<0) continue;    
                    fHistCasMC->Fill(par[0]-0.5,7.5);


                    AliMCParticle *mcBachelorMother = (AliMCParticle*)fmcEvent->GetTrack(cascadelabel);
                    if (!mcBachelorMother) continue;
                    fHistCasMC->Fill(par[0]-0.5,8.5);
                    cascadePDG = mcBachelorMother->PdgCode();
        
                    Bool_t IsPhysPrim = mcBachelorMother->IsPhysicalPrimary();
                    if(!IsPhysPrim) continue;

                    fHistCasMC->Fill(par[0]-0.5,9.5);
                        
                    Bool_t IsFromMC = ((cascadePDG==3312)&&(V0Pdg==3122)&&(PosTrackPdg==2212)&&(NegTrackPdg==-211)&&(BacTrackPdg==-211)) || ((cascadePDG==-3312)&&(V0Pdg==-3122)&&(NegTrackPdg==-2212)&&(PosTrackPdg==211)&&(BacTrackPdg==211));

                    if(!IsFromMC) continue;
                    fHistCasMC->Fill(par[0]-0.5,10.5);

                    if(fEfficiency){
                        Double_t cascadeEffic[6]={cascadept,fPV[2],par[0]-1.5,cascadeESD->Eta(),par[1],cascadeESD->Phi()};
                        fHistRecV0->Fill(cascadeEffic);
                    }
                }
            }
        }
        if(fV0hCorr||fhV0Corr){
	    for (Int_t i=0; i<nV0; i++){
        
            if(fAOD) {
                V0AOD = (AliAODv0*)(fAOD->GetV0(i));
                V0 = (AliVParticle*)(fAOD->GetV0(i));
                if(!V0AOD) {
                    cout << "Error: No AOD V0" << endl;
                    continue;
                }
            }
            if(fESD) {
                V0esd = (AliESDv0*)(fESD->GetV0(i));
                V0 = (AliVParticle*)(fESD->GetV0(i));
                if(!V0esd) continue;
            }

            if(!V0) {
                cout << "Error: No particle V0" << endl;
                continue;
            }
            
            

            if(fAOD){
                massK0 = V0AOD->MassK0Short();
                massLambda = V0AOD->MassLambda();
                massAntilambda = V0AOD->MassAntiLambda();
                v0pt = V0AOD->Pt();
                v0phi = V0AOD->Phi();
                v0Eta = V0AOD->Eta();
            }
            if(fESD){
                V0esd->ChangeMassHypothesis(310);
                massK0 = V0esd->GetEffMass();
                V0esd->ChangeMassHypothesis(3122);
                massLambda = V0esd->GetEffMass();
                V0esd->ChangeMassHypothesis(-3122);
                massAntilambda = V0esd->GetEffMass();
                v0pt = V0esd->Pt();
                v0phi = V0esd->Phi();
                v0Eta = V0esd->Eta();
            }

            if(v0pt<fPtAsocMin) continue; // pt assoc cut

            Bool_t k0 = ((massK0>0.44)&&(massK0<0.56));
            Bool_t Lambda = ((massLambda>1.08)&&(massLambda<1.15));
            Bool_t Antilambda = ((massAntilambda>1.08)&&(massAntilambda<1.15));

            // PID cut--------------------------
            Float_t nSigmaPosPion   = 0.;
            Float_t nSigmaNegPion   = 0.;
            Float_t nSigmaPosProton = 0.;
            Float_t nSigmaNegProton = 0.;

            Double_t posProp[6];
            Double_t negProp[6];

            if(fAOD){
                AliVTrack  * trackNegTest=dynamic_cast<AliVTrack *>(V0AOD->GetDaughter(1));
                AliVTrack  * trackPosTest=dynamic_cast<AliVTrack *>(V0AOD->GetDaughter(0));

                myTrackNegTest = (AliAODTrack*) trackNegTest;
                myTrackPosTest = (AliAODTrack*) trackPosTest;
      
                if (!myTrackPosTest || !myTrackNegTest) {
                    Printf("strange analysis::UserExec:: Error:Could not retreive one of the daughter track\n");
                    continue;
                }

                if(myTrackPosTest->Charge() == myTrackNegTest->Charge()) continue;

     
                if( myTrackPosTest->Charge() ==1){
                    myTrackPosAOD = myTrackPosTest;
                    myTrackNegAOD = myTrackNegTest;
                }
     
                if( myTrackPosTest->Charge() ==-1){
                    myTrackPosAOD = myTrackNegTest;
                    myTrackNegAOD = myTrackPosTest;
                }

                const AliAODPid *pPid = myTrackPosAOD->GetDetPid();
                const AliAODPid *nPid = myTrackNegAOD->GetDetPid();

                idNeg = myTrackNegAOD->GetID();
                idPos = myTrackPosAOD->GetID();

                if(fAnalysisMC){
                    labelNeg = myTrackNegAOD->GetLabel();
                    labelPos = myTrackPosAOD->GetLabel();
                }

                posProp[0] = myTrackPosAOD->Phi();
                posProp[1] = myTrackPosAOD->Pt();
                posProp[2] = myTrackPosAOD->Eta();
                posProp[3] = myTrackPosAOD->Charge();
                posProp[4] = idPos;
                posProp[5] = labelPos;
                
                negProp[0] = myTrackNegAOD->Phi();
                negProp[1] = myTrackNegAOD->Pt();
                negProp[2] = myTrackNegAOD->Eta();
                negProp[3] = myTrackNegAOD->Charge();
                negProp[4] = idNeg;
                negProp[5] = labelNeg;

                if (pPid){
                    Double_t pdMom = pPid->GetTPCmomentum();
                    if (pdMom<1.){
                        nSigmaPosPion = fPIDResponse->NumberOfSigmasTPC(myTrackPosAOD, AliPID::kPion);
                        nSigmaPosProton = fPIDResponse->NumberOfSigmasTPC(myTrackPosAOD, AliPID::kProton);
                    }
                }
            
                if (nPid){
                    Double_t ndMom = nPid->GetTPCmomentum();
                    if (ndMom<1.){
                        nSigmaNegPion = fPIDResponse->NumberOfSigmasTPC(myTrackNegAOD, AliPID::kPion);
                        nSigmaNegProton = fPIDResponse->NumberOfSigmasTPC(myTrackNegAOD, AliPID::kProton);
                    }
                }
            }
            if(fESD){

                UInt_t indexPos = (UInt_t)TMath::Abs(V0esd->GetPindex());
                UInt_t indexNeg = (UInt_t)TMath::Abs(V0esd->GetNindex());

                AliESDtrack *pTrack=(AliESDtrack*)fESD->GetTrack(indexPos);
                AliESDtrack *nTrack=(AliESDtrack*)fESD->GetTrack(indexNeg);

                if (!pTrack || !nTrack) {
                    Printf("ERROR: Could not retreive one of the daughter track");
                    continue;
                }
                if ( pTrack->GetSign() == nTrack->GetSign()) continue;

                if( pTrack->Charge() ==1){
                    myTrackPosESD = pTrack;
                    myTrackNegESD = nTrack;
                }
     
                if( nTrack->Charge() == 1){
                    myTrackPosESD = nTrack;
                    myTrackNegESD = pTrack;
                }

                idNeg = myTrackNegESD->GetID();
                idPos = myTrackPosESD->GetID();

                if(fAnalysisMC){
                    labelNeg = myTrackNegESD->GetLabel();
                    labelPos = myTrackPosESD->GetLabel();
                }

                posProp[0] = myTrackPosESD->Phi();
                posProp[1] = myTrackPosESD->Pt();
                posProp[2] = myTrackPosESD->Eta();
                posProp[3] = myTrackPosESD->Charge();
                posProp[4] = idPos;
                posProp[5] = labelPos;
                
                negProp[0] = myTrackNegESD->Phi();
                negProp[1] = myTrackNegESD->Pt();
                negProp[2] = myTrackNegESD->Eta();
                negProp[3] = myTrackNegESD->Charge();
                negProp[4] = idNeg;
                negProp[5] = labelNeg;

                if (fPIDResponse){
                    Double_t pdMom = myTrackPosESD->GetTPCmomentum();
                    if (pdMom<1.){
                        nSigmaPosPion = fPIDResponse->NumberOfSigmasTPC(myTrackPosESD, (AliPID::EParticleType)AliPID::kPion);
                        nSigmaPosProton = fPIDResponse->NumberOfSigmasTPC(myTrackPosESD, (AliPID::EParticleType)AliPID::kProton);
                    }
                }
                
                if (fPIDResponse){
                    Double_t ndMom = myTrackNegESD->GetTPCmomentum();
                    if (ndMom<1.){
                        nSigmaNegPion = fPIDResponse->NumberOfSigmasTPC(myTrackNegESD, (AliPID::EParticleType)AliPID::kPion);
                        nSigmaNegProton = fPIDResponse->NumberOfSigmasTPC(myTrackNegESD, (AliPID::EParticleType)AliPID::kProton);
                    }
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
        
            if(fAOD) {
                if (fRejectV0PileUp&&(!(((myTrackNegAOD->IsOn(AliAODTrack::kTPCrefit)&&myTrackNegAOD->IsOn(AliAODTrack::kITSrefit))||myTrackNegAOD->IsOn(AliAODTrack::kTOFout))||((myTrackPosAOD->IsOn(AliAODTrack::kTPCrefit)&&myTrackPosAOD->IsOn(AliAODTrack::kITSrefit))||myTrackPosAOD->IsOn(AliAODTrack::kTOFout))))) continue;
            }
            
            if(fESD){ 
                if (fRejectV0PileUp&&(!(((myTrackNegESD->IsOn(AliESDtrack::kTPCrefit)&&myTrackNegESD->IsOn(AliESDtrack::kITSrefit))||myTrackNegESD->IsOn(AliESDtrack::kTOFout))||((myTrackPosESD->IsOn(AliESDtrack::kTPCrefit)&&myTrackPosESD->IsOn(AliESDtrack::kITSrefit))||myTrackPosESD->IsOn(AliESDtrack::kTOFout))))) continue;
            }

            if (k0) {
                fHistK0MassPtCut->Fill(massK0,v0pt,0.5);
                nK0+=1;
            }
            if (Lambda) {
                fHistLambdaMassPtCut->Fill(massLambda,v0pt,0.5);
                nLam+=1;
            }
            if (Antilambda) {
                fHistAntiLambdaMassPtCut->Fill(massAntilambda,v0pt,0.5);
                nALam+=1;
            }
            
            if (k0&&cutK0Pid) fHistK0MassPtCut->Fill(massK0,v0pt,1.5);
            if (Lambda&&cutLambdaPid) fHistLambdaMassPtCut->Fill(massLambda,v0pt,1.5);
            if (Antilambda&&cutAntiLambdaPid) fHistAntiLambdaMassPtCut->Fill(massAntilambda,v0pt,1.5);
            
            Int_t oStatus = GetOStatus();
            if(fAOD){
                if(!IsMyGoodV0(V0AOD,myTrackPosAOD,myTrackNegAOD,oStatus)) continue; // on fly and daughters cuts
            }
            if(fESD){
                if(!IsMyGoodV0ESD(V0esd,myTrackPosESD,myTrackNegESD,oStatus)) continue;
            }

            if (k0&&cutK0Pid) fHistK0MassPtCut->Fill(massK0,v0pt,2.5);
            if (Lambda&&cutLambdaPid) fHistLambdaMassPtCut->Fill(massLambda,v0pt,2.5);
            if (Antilambda&&cutAntiLambdaPid) fHistAntiLambdaMassPtCut->Fill(massAntilambda,v0pt,2.5);
        
        //======= crosscheck topological cuts
            if(fCutsCrosscheck&&fAOD){
                Double_t DCANegToPV = V0AOD->DcaNegToPrimVertex();
                Double_t DCAPosToPV = V0AOD->DcaPosToPrimVertex();
                Double_t DCADaught = V0AOD->DcaV0Daughters();
                Double_t V0radius = V0AOD->RadiusV0();
                Double_t CosPA = V0AOD->CosPointingAngle(myPrimVertex);
                Double_t massSellLam = TMath::Abs(massK0-0.497614);
                Double_t massSellK0 = TMath::Abs(massLambda-1.115683);
            
                Double_t *tParentVertexPosition = new Double_t[3];
                tParentVertexPosition[0]= fPV[0];
                tParentVertexPosition[1]= fPV[1];
                tParentVertexPosition[2]= fPV[2];
            
                Double_t lenght = V0AOD->DecayLengthV0(tParentVertexPosition);
                Double_t momentum = TMath::Sqrt(TMath::Power(V0AOD->Px(),2)+TMath::Power(V0AOD->Py(),2)+TMath::Power(V0AOD->Pz(),2));
                
                delete[] tParentVertexPosition;

                Double_t lifetimeK0 = (massK0*lenght)/momentum;
                Double_t lifetimeLam = (massLambda*lenght)/momentum;
                Double_t lifetimeALam = (massAntilambda*lenght)/momentum;
            
                if(k0&&cutK0Pid&&IsMyGoodV0RapidityK0(V0AOD)) {
                    if (V0AOD->GetOnFlyStatus()) TopologCuts(v0pt,massK0,DCANegToPV,DCAPosToPV,DCADaught,V0radius,CosPA,lifetimeK0,massSellK0,0.5,0.5);
                    if (!V0AOD->GetOnFlyStatus()) TopologCuts(v0pt,massK0,DCANegToPV,DCAPosToPV,DCADaught,V0radius,CosPA,lifetimeK0,massSellK0,0.5,1.5);
                }
            
                if (Lambda&&cutLambdaPid&&IsMyGoodV0RapidityLambda(V0AOD)) {
                    if (V0AOD->GetOnFlyStatus()) TopologCuts(v0pt,massLambda,DCANegToPV,DCAPosToPV,DCADaught,V0radius,CosPA,lifetimeLam,massSellLam,1.5,0.5);
                    if (!V0AOD->GetOnFlyStatus()) TopologCuts(v0pt,massLambda,DCANegToPV,DCAPosToPV,DCADaught,V0radius,CosPA,lifetimeLam,massSellLam,1.5,1.5);
                }
              
                if (Antilambda&&cutAntiLambdaPid&&IsMyGoodV0RapidityLambda(V0AOD)) {
                    if (V0AOD->GetOnFlyStatus()) TopologCuts(v0pt,massAntilambda,DCANegToPV,DCAPosToPV,DCADaught,V0radius,CosPA,lifetimeALam,massSellLam,2.5,0.5);
                    if (!V0AOD->GetOnFlyStatus()) TopologCuts(v0pt,massAntilambda,DCANegToPV,DCAPosToPV,DCADaught,V0radius,CosPA,lifetimeALam,massSellLam,2.5,1.5);
                }
            
                if(fAnalysisMC){
                    daughter1= static_cast<AliAODTrack*> (V0AOD->GetDaughter(1));
                    daughter0= static_cast<AliAODTrack*> (V0AOD->GetDaughter(0));
                    if (!(daughter1->GetLabel()<0||daughter0->GetLabel()<0)){
                
                        Int_t mcmother1 = static_cast<AliAODMCParticle*>(fmcEvent->GetTrack(daughter1->GetLabel()))->GetMother();
                        Int_t mcmother0 = static_cast<AliAODMCParticle*>(fmcEvent->GetTrack(daughter0->GetLabel()))->GetMother();
                
                        if(mcmother1==mcmother0){
                            Int_t pdgV0 = static_cast<AliAODMCParticle*>(fmcEvent->GetTrack(mcmother0))->PdgCode();
                            Int_t pdgD1 = static_cast<AliAODMCParticle*>(fmcEvent->GetTrack(daughter1->GetLabel()))->PdgCode();
                            Int_t pdgD0 = static_cast<AliAODMCParticle*>(fmcEvent->GetTrack(daughter0->GetLabel()))->PdgCode();
                
                            Bool_t isPhyPrim = static_cast<AliAODMCParticle*>(fmcEvent->GetTrack(mcmother0))->IsPhysicalPrimary();
                            if(isPhyPrim) {
                                Double_t V0mcPt = static_cast<AliAODMCParticle*>(fmcEvent->GetTrack(mcmother0))->Pt();
                                if(V0mcPt>=fPtTrigMin) {
                                    if (pdgV0==3122 && ((pdgD0==2212 && pdgD1==-211 )||(pdgD0==-211 && pdgD1==2212))&&IsMyGoodV0RapidityLambda(V0AOD)) {
                                        if (V0AOD->GetOnFlyStatus()) TopologCuts(v0pt,massLambda,DCANegToPV,DCAPosToPV,DCADaught,V0radius,CosPA,lifetimeLam,massSellLam,1.5,0.5);
                                        if (!V0AOD->GetOnFlyStatus()) TopologCuts(v0pt,massLambda,DCANegToPV,DCAPosToPV,DCADaught,V0radius,CosPA,lifetimeLam,massSellLam,1.5,1.5);
                                    }
                
                                    if ((pdgV0==310) && (( pdgD0==211 && pdgD1==-211 )||( pdgD0==-211 && pdgD1==211 ))&&IsMyGoodV0RapidityK0(V0AOD)){
                                        if (V0AOD->GetOnFlyStatus()) TopologCuts(v0pt,massK0,DCANegToPV,DCAPosToPV,DCADaught,V0radius,CosPA,lifetimeK0,massSellK0,0.5,0.5);
                                        if (!V0AOD->GetOnFlyStatus()) TopologCuts(v0pt,massK0,DCANegToPV,DCAPosToPV,DCADaught,V0radius,CosPA,lifetimeK0,massSellK0,0.5,1.5);
                                    }
                
                                    if ((pdgV0==-3122)&& ((pdgD0==-2212 && pdgD1==211 )||(pdgD0==211 && pdgD1==-2212 ))&&IsMyGoodV0RapidityLambda(V0AOD)) {
                                        if (V0AOD->GetOnFlyStatus()) TopologCuts(v0pt,massAntilambda,DCANegToPV,DCAPosToPV,DCADaught,V0radius,CosPA,lifetimeALam,massSellLam,2.5,0.5);
                                        if (!V0AOD->GetOnFlyStatus()) TopologCuts(v0pt,massAntilambda,DCANegToPV,DCAPosToPV,DCADaught,V0radius,CosPA,lifetimeALam,massSellLam,2.5,1.5);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        
        //------------------- V0 cuts -------------------------------
        
            if (k0&&cutK0Pid) fHistK0MassPtCut->Fill(massK0,v0pt,3.5);
            if (Lambda&&cutLambdaPid) fHistLambdaMassPtCut->Fill(massLambda,v0pt,3.5);
            if (Antilambda&&cutAntiLambdaPid) fHistAntiLambdaMassPtCut->Fill(massAntilambda,v0pt,3.5);
        
            if(fAOD)
                if(!IsMyGoodV0Topology(V0AOD)) continue; //topoligical cuts
            if(fESD)
                if(!IsMyGoodV0TopologyESD(V0esd, myTrackPosESD, myTrackNegESD)) continue;


            if (k0&&cutK0Pid) fHistK0MassPtCut->Fill(massK0,v0pt,4.5);
            if (Lambda&&cutLambdaPid) fHistLambdaMassPtCut->Fill(massLambda,v0pt,4.5);
            if (Antilambda&&cutAntiLambdaPid) fHistAntiLambdaMassPtCut->Fill(massAntilambda,v0pt,4.5);

            Bool_t RapidityCut = kFALSE;
            Bool_t LifetimeCut = kFALSE;
            Bool_t CosPointingAngleCut = kFALSE;

            if (Lambda&&cutLambdaPid){
            
                if(fAOD) RapidityCut = IsMyGoodV0RapidityLambda(V0AOD);
                if(fESD) RapidityCut = IsMyGoodV0RapidityLambdaESD(V0esd);
                
                if(RapidityCut){ //Rapidity
                    fHistLambdaMassPtCut->Fill(massLambda,v0pt,5.5);
            
                    if(fAOD) LifetimeCut = IsMyGoodLifeTimeLambda(V0AOD);
                    if(fESD) LifetimeCut = IsMyGoodLifeTimeESD(V0esd,massLambda);
                
                    if(LifetimeCut) { //Proper Lifetime (mL/p)
                        fHistLambdaMassPtCut->Fill(massLambda,v0pt,6.5);
            
                        if(fAOD) CosPointingAngleCut = IsMyGoodV0AngleLambda(V0AOD,myPrimVertex);    
                        if(fESD) CosPointingAngleCut = IsMyGoodV0AngleLambdaESD(V0esd);

                        if(CosPointingAngleCut){ //V0 Cosine of Pointing Angle
                            fHistLambdaMassPtCut->Fill(massLambda,v0pt,7.5);
            
                            if (TMath::Abs(massK0-0.497614)>fMassRejectCutLam){
                                fHistLambdaMassPtCut->Fill(massLambda,v0pt,8.5);
            
                                if(fAnalysisMC){
                                    FillMC(V0,3122,2212, -211,2,massLambda,posProp,negProp);
                                    Double_t phiEtaLamMC[4] = {v0pt,v0phi,v0Eta,1.5};
                                    fHistPhiEta->Fill(phiEtaLamMC);
                                }
                                if(!fAnalysisMC) {
                                    if(V0->Pt()>fPtTrigMin) fselectedV0Triggers-> Add(new AliV0ChParticle(v0Eta, v0phi, v0pt, 2,0,(Int_t)posProp[4],(Int_t)negProp[4],massLambda,posProp[0],posProp[1],posProp[2],(Int_t)posProp[3],negProp[0],negProp[1],negProp[2],(Int_t)negProp[3]));
                                    Double_t phiEtaLamData[4] = {v0pt,v0phi,v0Eta,1.5};
                                    fHistPhiEta->Fill(phiEtaLamData);
                                    fselectedV0Assoc-> Add(new AliV0ChParticle(v0Eta, v0phi, v0pt, 6,0,(Int_t)posProp[4],(Int_t)negProp[4],massLambda,posProp[0],posProp[1],posProp[2],(Int_t)posProp[3],negProp[0],negProp[1],negProp[2],(Int_t)negProp[3]));
                                }
                            }
                        }
                    }
                }
            }
            RapidityCut = kFALSE;
            LifetimeCut =kFALSE;
            CosPointingAngleCut=kFALSE;
        
            if (k0&&cutK0Pid){
                if(fAOD) RapidityCut = IsMyGoodV0RapidityK0(V0AOD);
                if(fESD) RapidityCut = IsMyGoodV0RapidityK0ESD(V0esd);

                if(RapidityCut){ //Rapidity
                    fHistK0MassPtCut->Fill(massK0,v0pt,5.5);
            
                    if(fAOD) LifetimeCut = IsMyGoodLifeTimeK0(V0AOD);
                    if(fESD) LifetimeCut = IsMyGoodLifeTimeESD(V0esd,massK0);
                
                    if(LifetimeCut) { //Proper Lifetime (mL/p)
                        fHistK0MassPtCut->Fill(massK0,v0pt,6.5);
            
                        if(fAOD) CosPointingAngleCut = IsMyGoodV0AngleK0(V0AOD,myPrimVertex);    
                        if(fESD) CosPointingAngleCut = IsMyGoodV0AngleK0ESD(V0esd);
                
                        if (CosPointingAngleCut){ //V0 Cosine of Pointing Angle
                            fHistK0MassPtCut->Fill(massK0,v0pt,7.5);
            
                            if(TMath::Abs(massLambda-1.115683)>fMassRejectCutK0&&TMath::Abs(massAntilambda-1.115683)>fMassRejectCutK0) {
                                fHistK0MassPtCut->Fill(massK0,v0pt,8.5);
                
                                if(fAnalysisMC){
                                    FillMC(V0,310,211, -211,1,massK0,posProp,negProp);
                                    Double_t phiEtaK0MC[4] = {v0pt,v0phi,v0Eta,0.5};
                                    fHistPhiEta->Fill(phiEtaK0MC);
                                }
                                if(!fAnalysisMC) {
                                    if(V0->Pt()>fPtTrigMin) fselectedV0Triggers-> Add(new AliV0ChParticle(v0Eta, v0phi, v0pt,1,0,(Int_t)posProp[4],(Int_t)negProp[4],massK0,posProp[0],posProp[1],posProp[2],(Int_t)posProp[3],negProp[0],negProp[1],negProp[2],(Int_t)negProp[3]));
                                    Double_t phiEtaK0Data[4] = {v0pt,v0phi,v0Eta,0.5};
                                    fHistPhiEta->Fill(phiEtaK0Data);
                                    fselectedV0Assoc-> Add(new AliV0ChParticle(v0Eta, v0phi, v0pt,5,0,(Int_t)posProp[4],(Int_t)negProp[4],massK0,posProp[0],posProp[1],posProp[2],(Int_t)posProp[3],negProp[0],negProp[1],negProp[2],(Int_t)negProp[3]));
                                }
                            }
                        }
                    }
                }
            }
            RapidityCut = kFALSE;
            LifetimeCut =kFALSE;
            CosPointingAngleCut=kFALSE;
            
            if (Antilambda&&cutAntiLambdaPid){

                if(fAOD) RapidityCut = IsMyGoodV0RapidityLambda(V0AOD);
                if(fESD) RapidityCut = IsMyGoodV0RapidityLambdaESD(V0esd);

                if (RapidityCut) {
                    fHistAntiLambdaMassPtCut->Fill(massAntilambda,v0pt,5.5);

                    if(fAOD) LifetimeCut = IsMyGoodLifeTimeK0(V0AOD);
                    if(fESD) LifetimeCut = IsMyGoodLifeTimeESD(V0esd,massAntilambda);

                    if(LifetimeCut){
                        fHistAntiLambdaMassPtCut->Fill(massAntilambda,v0pt,6.5);

                        if(fAOD) CosPointingAngleCut = IsMyGoodV0AngleK0(V0AOD,myPrimVertex);    
                        if(fESD) CosPointingAngleCut = IsMyGoodV0AngleK0ESD(V0esd);

                        if (CosPointingAngleCut){
                            fHistAntiLambdaMassPtCut->Fill(massAntilambda,v0pt,7.5);
                
                            if(TMath::Abs(massK0-0.497614)>fMassRejectCutLam){
                                fHistAntiLambdaMassPtCut->Fill(massAntilambda,v0pt,8.5);

                
                                if(fAnalysisMC){
                                    FillMC(V0,-3122,211, -2212,3,massAntilambda,posProp,negProp);
                                    Double_t phiEtaAlamMC[4] = {v0pt,v0phi,v0Eta,2.5};
                                    fHistPhiEta->Fill(phiEtaAlamMC);
                                    
                                }
                                if(!fAnalysisMC) {
                                    if(V0->Pt()>fPtTrigMin) fselectedV0Triggers-> Add(new AliV0ChParticle(v0Eta, v0phi, v0pt, 3,0,(Int_t)posProp[4],(Int_t)negProp[4],massAntilambda,posProp[0],posProp[1],posProp[2],(Int_t)posProp[3],negProp[0],negProp[1],negProp[2],(Int_t)negProp[3]));
                                    Double_t phiEtaAlamData[4] = {v0pt,v0phi,v0Eta,2.5};
                                    fHistPhiEta->Fill(phiEtaAlamData);
                                    fselectedV0Assoc-> Add(new AliV0ChParticle(v0Eta, v0phi, v0pt, 7,0,(Int_t)posProp[4],(Int_t)negProp[4],massAntilambda,posProp[0],posProp[1],posProp[2],(Int_t)posProp[3],negProp[0],negProp[1],negProp[2],(Int_t)negProp[3]));
                                }
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
            if(fV0hCorr) Corelations(fselectedMCV0Triggersrec,fselectedMCassoc,kFALSE,kTRUE,lPercentile,kFALSE);

            //h-h MC rec
            if(fhhCorr) Corelations(fselectedMCtrig,fselectedMCassoc,kTRUE,kFALSE,lPercentile,kFALSE);
             
             //MC rec h-V0
            if(fhV0Corr) Corelations(fselectedMCtrig,fselectedMCV0assoc,kFALSE,kTRUE,lPercentile,kTRUE);
        
            if(fAnalyseFeedDown) {
                CorrelationsXi(fselectedV0Triggers,fselectedMCassoc,lPercentile,kTRUE);
                CorrelationsXi(fselectedMCtrig,fselectedV0Assoc,lPercentile,kFALSE);
            }
        } else if(fCorrelations){
            //Data V0-h
            if(fV0hCorr) Corelations(fselectedV0Triggers,fselectedAssociatedTracks,kFALSE,kTRUE,lPercentile,kFALSE);

    	    //Data h-h
            if(fhhCorr) Corelations(fselectedTriggerTracks,fselectedAssociatedTracks,kTRUE,kFALSE,lPercentile,kFALSE);
            
            //Data h-V0
            if(fhV0Corr) Corelations(fselectedTriggerTracks,fselectedV0Assoc,kFALSE,kTRUE,lPercentile,kTRUE);
        
            if(fAnalyseFeedDown) {
                CorrelationsXi(fselectedV0Triggers,fselectedAssociatedTracks,lPercentile,kTRUE);
                CorrelationsXi(fselectedTriggerTracks,fselectedV0Assoc,lPercentile,kFALSE);
            }
        }
    }    

 	// Mixing ==============================================

    fHistMultVtxz->Fill(lPercentile,fPV[2]);
    if(fMixing){
        fPool = fPoolMgr->GetEventPool(lPercentile, fPV[2]);
        if (!fPool) {
            AliWarning(Form("No pool found for centrality = %f, zVtx = %f", lPercentile, fPV[2]));
            return;
        }
        Int_t nMix = fPool->GetCurrentNEvents();
        if (fPool->IsReady() || fPool->NTracksInPool() > fMixingTracks / 10 || nMix >= fMixedEvents)
        {
            for (Int_t jMix=0; jMix<nMix; jMix++)
            {// loop through mixing events
                TObjArray* bgTracks = fPool->GetEvent(jMix);
                if(fAnalysisMC) {
                    if(fV0hCorr) CorelationsMixing(fselectedMCV0Triggersrec,bgTracks,lPercentile);
                    if(fhhCorr) CorelationsMixing(fselectedMCtrig,bgTracks,lPercentile);
                    if(fhV0Corr) CorelationsMixinghV0(bgTracks,fselectedMCV0assoc,lPercentile);
                    if(fAnalyseFeedDown){
                        CorrelationsXi(fselectedV0Triggers,bgTracks,lPercentile,kTRUE);
                        CorrelationsXi(bgTracks,fselectedV0Assoc,lPercentile,kFALSE);
                    }
                }else{
                    if(fV0hCorr) CorelationsMixing(fselectedV0Triggers,bgTracks,lPercentile);
                    if(fhhCorr) CorelationsMixing(fselectedTriggerTracks,bgTracks,lPercentile);
                    if(fhV0Corr) CorelationsMixinghV0(bgTracks,fselectedV0Assoc,lPercentile);
                    if(fAnalyseFeedDown){
                        CorrelationsXi(fselectedV0Triggers,bgTracks,lPercentile,kTRUE);
                        CorrelationsXi(bgTracks,fselectedV0Assoc,lPercentile,kFALSE);
                    }
                }
             }
        }
        TObjArray* cloneArray = (TObjArray *)fselectedTracks->Clone();
        cloneArray->SetOwner(kTRUE);
        fPool->UpdatePool(cloneArray);
    }

    if(fAnalysisMC&&fMixingGen){
        fPoolMCGen = fPoolMgr->GetEventPool(lPercentile, fPV[2]);
        if (!fPoolMCGen) {
            AliWarning(Form("No pool MC Gen found for centrality = %f, zVtx = %f", lPercentile, fPV[2]));
            return;
        }
        Int_t nMixGen = fPoolMCGen->GetCurrentNEvents();
        if (fPoolMCGen->IsReady() || fPoolMCGen->NTracksInPool() > fMixingTracks / 10 || nMixGen >= fMixedEvents)
        {   
            for (Int_t jMix=0; jMix<nMixGen; jMix++){
                TObjArray* bgTracksGen = fPoolMCGen->GetEvent(jMix);
                if(fV0hCorr) CorelationsMixing(fmcTracksV0Sel,bgTracksGen,nAcceptedParticles);
                if(fhhCorr) CorelationsMixing(fmcTracksTrigSel,bgTracksGen,nAcceptedParticles);
                if(fhV0Corr) CorelationsMixinghV0(bgTracksGen,fmcV0AssocSel,nAcceptedParticles);
            }
        }

        TObjArray* cloneArrayMCGen = (TObjArray *)fmcGenTracksMixing->Clone();
        cloneArrayMCGen->SetOwner(kTRUE);
        fPoolMCGen->UpdatePool(cloneArrayMCGen);
    }    
    // deleting TObjArrays
    fmcTracksSel->Clear();
    fmcTracksTrigSel->Clear();
    fmcTracksV0Sel->Clear();
    fselectedMCassoc->Clear();
    fselectedMCtrig->Clear();
    fselectedMCV0Triggersrec->Clear();
    fselectedAssociatedTracks->Clear();
    fselectedTriggerTracks->Clear();
    fselectedV0Triggers->Clear();
    fselectedV0Assoc->Clear();
    fmcV0AssocSel->Clear();
    fselectedMCV0assoc->Clear();
    fmcGenTracksMixing->Clear();
    fselectedTracks->Clear();
    
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
//_____________________________________________________________________________
Bool_t AliAnalysisTaskDiHadCorrelHighPt::IsMyGoodPrimaryTrackESD(const AliESDtrack *t)
 {      
          // Pseudorapidity cut
          if (TMath::Abs(t->Eta())>=fEtaCut) return kFALSE;
          const AliESDEvent* esdEvent = (AliESDEvent*) t->GetESDEvent();
          if(!esdEvent) {
            cout << "No ESD event" << endl;
            return kFALSE;
          }
          if(fFilterBit==16||fFilterBit==256){
            fESDTrackCuts->SetMaxDCAToVertexXY(2.4); 
            fESDTrackCuts->SetMaxDCAToVertexZ(3.2); 
            fESDTrackCuts->SetDCAToVertex2D(kTRUE);
            if(fFilterBit==256){
                fESDTrackCuts->SetMaxChi2TPCConstrainedGlobal(36);
                fESDTrackCuts->SetMaxFractionSharedTPCClusters(0.4);
            }
          }
          if(!fESDTrackCuts->AcceptTrack(t)) return kFALSE;
  
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
Bool_t AliAnalysisTaskDiHadCorrelHighPt::IsMyGoodV0AngleK0ESD(const AliESDv0 *t)
{
        //V0 Cosine of Pointing Angle
        if(t->GetV0CosineOfPointingAngle()<=fCosPointAngleK0) return kFALSE;
    
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
Bool_t AliAnalysisTaskDiHadCorrelHighPt::IsMyGoodV0AngleLambdaESD(const AliESDv0 *t)
{
        //V0 Cosine of Pointing Angle
        if(t->GetV0CosineOfPointingAngle()<=fCosPointAngleLam) return kFALSE;
        
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
Bool_t AliAnalysisTaskDiHadCorrelHighPt::IsMyGoodV0RapidityLambdaESD(const AliESDv0 *t)
{
        //Rapidity
        if(TMath::Abs(t->RapLambda())>=fRapidityCut) return kFALSE;
        
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
Bool_t AliAnalysisTaskDiHadCorrelHighPt::IsMyGoodV0RapidityK0ESD(const AliESDv0 *t)
{
        //Rapidity
        if(TMath::Abs(t->RapK0Short())>=fRapidityCut) return kFALSE;
        
        return kTRUE;
}
//____________________________________________________________________________
Bool_t AliAnalysisTaskDiHadCorrelHighPt::IsMyGoodLifeTimeK0(const AliAODv0 *V0)
{
		//Proper Lifetime (mL/p)

		Double_t lenght = V0->DecayLengthV0(fPV);
		Double_t momentum = TMath::Sqrt(TMath::Power(V0->Px(),2)+TMath::Power(V0->Py(),2)+TMath::Power(V0->Pz(),2));
		Double_t mass = V0->MassK0Short();
		Double_t lifetime = (mass*lenght)/momentum;

		if(lifetime>=fLifeTimeK0) return kFALSE;
		
		return kTRUE;
}
//____________________________________________________________________________
Bool_t AliAnalysisTaskDiHadCorrelHighPt::IsMyGoodLifeTimeLambda(const AliAODv0 *V0)
{
		//Proper Lifetime (mL/p)

		Double_t lenght = V0->DecayLengthV0(fPV);
		Double_t momentum = TMath::Sqrt(TMath::Power(V0->Px(),2)+TMath::Power(V0->Py(),2)+TMath::Power(V0->Pz(),2));
		Double_t mass = V0->MassLambda();
		Double_t lifetime = (mass*lenght)/momentum;
		
		if(lifetime>=fLifeTimeLam) return kFALSE;
		
		return kTRUE;
}
//____________________________________________________________________________
Bool_t AliAnalysisTaskDiHadCorrelHighPt::IsMyGoodLifeTimeAntiLambda(const AliAODv0 *V0)
{
    //Proper Lifetime (mL/p)
    
    Double_t lenght = V0->DecayLengthV0(fPV);
    Double_t momentum = TMath::Sqrt(TMath::Power(V0->Px(),2)+TMath::Power(V0->Py(),2)+TMath::Power(V0->Pz(),2));
    Double_t mass = V0->MassAntiLambda();
    Double_t lifetime = (mass*lenght)/momentum;
    
    if(lifetime>=fLifeTimeLam) return kFALSE;
    
    return kTRUE;
}
//____________________________________________________________________________
Bool_t AliAnalysisTaskDiHadCorrelHighPt::IsMyGoodLifeTimeESD(const AliESDv0 *V0, Double_t masshypotesis)
{
        //Proper Lifetime (mL/p)
        
        Double_t tDecayVertexV0[3];
        V0->GetXYZ(tDecayVertexV0[0],tDecayVertexV0[1],tDecayVertexV0[2]);

        Double_t lenght = TMath::Sqrt(
                                    TMath::Power( tDecayVertexV0[0] - fPV[0] , 2) +
                                    TMath::Power( tDecayVertexV0[1] - fPV[1] , 2) +
                                    TMath::Power( tDecayVertexV0[2] - fPV[2] , 2));
        Double_t momentum = TMath::Sqrt(TMath::Power(V0->Px(),2)+TMath::Power(V0->Py(),2)+TMath::Power(V0->Pz(),2));
        Double_t lifetime = (masshypotesis*lenght)/momentum;
        
        if(lifetime>=fLifeTimeLam) return kFALSE;
        
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
Bool_t AliAnalysisTaskDiHadCorrelHighPt::IsMyGoodDaughterTrackESD(const AliESDtrack *t) {
    // TPC refit

    if (!t->IsOn(AliESDtrack::kTPCrefit)) return kFALSE;
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
//____________________________________________________________________________
Bool_t AliAnalysisTaskDiHadCorrelHighPt::IsMyGoodV0ESD(const AliESDv0 *v0,const AliESDtrack* myTrackPos, const AliESDtrack* myTrackNeg, Int_t oSta) {
    if (!v0) {
         AliError(Form("ERROR: Could not retrieve esdV0"));
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
    if ( !(IsMyGoodDaughterTrackESD(myTrackPos)) || !(IsMyGoodDaughterTrackESD(myTrackNeg)) ) return kFALSE;
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
//______________________________________________________________________
Bool_t AliAnalysisTaskDiHadCorrelHighPt::IsMyGoodV0TopologyESD(const AliESDv0 *v0, const AliESDtrack* myTrackPos, const AliESDtrack* myTrackNeg){
    if (!v0) {
        AliError(Form("ERROR: Could not retrieve esdV0"));
        return kFALSE;
    }
    //DCA Positive Track to PV
    Double_t lDcaPosToPrimVertex = TMath::Abs(myTrackPos->GetD(fPV[0], fPV[1],fMagneticField) );
    if(lDcaPosToPrimVertex<=fDCAposDaughter) return kFALSE;
    //DCA Negative Track to PV
    Double_t lDcaNegToPrimVertex = TMath::Abs(myTrackNeg->GetD(fPV[0], fPV[1],fMagneticField) );
    if(lDcaNegToPrimVertex<=fDCAnegDaughter) return kFALSE;
    //DCA V0 daughters
    if(v0->GetDcaV0Daughters()>=fDCAV0Daughters) return kFALSE;

    //V0 2D Decay Radius

    Double_t tDecayVertexV0[3];
    v0->GetXYZ(tDecayVertexV0[0],tDecayVertexV0[1],tDecayVertexV0[2]);
    Double_t lV0Radius = TMath::Sqrt(tDecayVertexV0[0]*tDecayVertexV0[0]+tDecayVertexV0[1]*tDecayVertexV0[1]);

    if(lV0Radius<=fV0Radius) return kFALSE;
    
    return kTRUE;
}
//____________________________________________________________________________
void AliAnalysisTaskDiHadCorrelHighPt::Corelations(TObjArray *triggers, TObjArray *associated,Bool_t hh,Bool_t V0h,Float_t perc,Bool_t hV0){

    const Double_t kPi = TMath::Pi();
    Int_t nAssoc = associated->GetEntriesFast();
    Int_t nTrig = triggers->GetEntriesFast();
    Double_t assocPt =0.;
    Double_t triggPt =0.;
    Double_t asocEta =0;
    Double_t assocPhi =0;
    Double_t triggEta =0.;
    Double_t assocCharge =0.;
    Double_t weight =1;
    Double_t triggPhi =0;
    Double_t massTrig = 0.;

    for (Int_t i=0; i<nTrig; i++){
         weight =1;
        AliV0ChParticle* trig = (AliV0ChParticle*)  triggers->At(i);
        triggPt = trig->Pt();
        triggEta = trig->Eta();
        triggPhi = trig->Phi(); 
        if(triggPt<fPtTrigMin) continue;
        if(trig->WhichCandidate()>7) continue;

        if(trig->WhichCandidate()<4) massTrig=trig->M();

        Int_t idbintrigg [4];
        if(!fCorrelationsGen&&trig->WhichCandidate()<4){
            idbintrigg[0]=fHistEffCorrectionK0->GetAxis(0)->FindBin(triggPt);
            idbintrigg[1]=fHistEffCorrectionK0->GetAxis(1)->FindBin(triggEta);
            idbintrigg[2]=fHistEffCorrectionK0->GetAxis(2)->FindBin(triggPhi);
            idbintrigg[3]=fHistEffCorrectionK0->GetAxis(3)->FindBin(fPV[2]);
        }else if(!fCorrelationsGen){
            idbintrigg[0]=fHistEffCorrectionHadron->GetAxis(0)->FindBin(triggPt);
            idbintrigg[1]=fHistEffCorrectionHadron->GetAxis(1)->FindBin(triggEta);
            idbintrigg[2]=fHistEffCorrectionHadron->GetAxis(2)->FindBin(triggPhi);
            idbintrigg[3]=fHistEffCorrectionHadron->GetAxis(3)->FindBin(fPV[2]);
        } 
        
        if(!hV0){
            if(!fCorrelationsGen){
                if(trig->WhichCandidate()==1) weight = fHistEffCorrectionK0->GetBinContent(idbintrigg);
                else if (trig->WhichCandidate()==2) weight = fHistEffCorrectionLam->GetBinContent(idbintrigg);
                else if (trig->WhichCandidate()==3) weight = fHistEffCorrectionAntiLam->GetBinContent(idbintrigg);
                else weight = fHistEffCorrectionHadron->GetBinContent(idbintrigg)/(1. - fHistSecondaryCont->Interpolate(triggPt));
            }
            if(weight ==0) continue;

            Double_t triggers[5]={triggPt,fPV[2],trig->WhichCandidate()-0.5,massTrig,perc};
            if(fCorrelationsGen&&fAnalysisMC) fHistNumberOfTriggersGen->Fill(triggers); 
            else if(fAnalysisMC&&fCorrelations) fHistNumberOfTriggersRec->Fill(triggers,1./weight);
            else fHistNumberOfTriggers->Fill(triggers,1./weight);
        }else{
            if(!fCorrelationsGen) weight = fHistEffCorrectionHadron->GetBinContent(idbintrigg)/(1. - fHistSecondaryCont->Interpolate(triggPt));

            if(weight ==0) continue;
            Double_t triggers[5]={triggPt,fPV[2],4.5,massTrig,perc};
            if(fCorrelationsGen&&fAnalysisMC) fHistNumberOfTriggersGen->Fill(triggers); 
            else if(fAnalysisMC&&fCorrelations) fHistNumberOfTriggersRec->Fill(triggers,1./weight);
            else fHistNumberOfTriggers->Fill(triggers,1./weight);
        }
        
        for (Int_t j=0; j<nAssoc; j++){
             weight =1;
            AliV0ChParticle* assoc = (AliV0ChParticle*)  associated->At(j);
            if(assoc->WhichCandidate() < 4 || assoc->WhichCandidate()>7) continue;
            asocEta = assoc->Eta();
            assocPhi = assoc->Phi();
            assocCharge = assoc->Charge();

            Double_t deltaEta = triggEta - asocEta;
            Double_t deltaPhi = triggPhi - assocPhi;
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
            Double_t massRho=10;
            Double_t massPhi=10;
            Double_t massDelta =10;
            Double_t massJPsi1=10;
            Double_t massJPsi2=10;
            Double_t massKStar = 10;
            Double_t massD0 =10;

            if(hh&&fRemoveHadrFromV0&&!fCorrelationsGen){
                if(trig->Charge()!=assocCharge){
                    massK0 = TMath::Sqrt(2*0.13957*0.13957+2*(trig->E()*assoc->E()-triggPt*assocPt-trig->Pz()*assoc->Pz()));
                    massLam = TMath::Sqrt(0.13957*0.13957+0.93827*0.93827+2*(trig->E()*assoc->E()-triggPt*assocPt-trig->Pz()*assoc->Pz()));
                    massGamma = TMath::Sqrt(2*0.0005109*0.0005109+2*(trig->E()*assoc->E()-triggPt*assocPt-trig->Pz()*assoc->Pz()));
                    massRho = TMath::Sqrt(2*0.13957*0.13957+2*(trig->E()*assoc->E()-triggPt*assocPt-trig->Pz()*assoc->Pz()));
                    massPhi = TMath::Sqrt(2*0.493677*0.493677+2*(trig->E()*assoc->E()-triggPt*assocPt-trig->Pz()*assoc->Pz())); 
                    massDelta = TMath::Sqrt(0.13957*0.13957+0.93827*0.93827+2*(trig->E()*assoc->E()-triggPt*assocPt-trig->Pz()*assoc->Pz()));
                    massJPsi1 = TMath::Sqrt(2*0.000511*0.000511+2*(trig->E()*assoc->E()-triggPt*assocPt-trig->Pz()*assoc->Pz())); 
                    massJPsi2 = TMath::Sqrt(2*0.10565*0.10565+2*(trig->E()*assoc->E()-triggPt*assocPt-trig->Pz()*assoc->Pz())); 
                    massKStar = TMath::Sqrt(0.13957*0.13957+0.493677*0.493677+2*(trig->E()*assoc->E()-triggPt*assocPt-trig->Pz()*assoc->Pz()));
                    massD0 = TMath::Sqrt(0.13957*0.13957+0.493677*0.493677+2*(trig->E()*assoc->E()-triggPt*assocPt-trig->Pz()*assoc->Pz()));//*/
                }else if(trig->Charge()==assocCharge){
                    massDelta = TMath::Sqrt(0.13957*0.13957+0.93827*0.93827+2*(trig->E()*assoc->E()-triggPt*assocPt-trig->Pz()*assoc->Pz()));
                }
                if(TMath::Abs( 0.497614-massK0)< 0.005) continue;
                if(TMath::Abs( 1.115683-massLam)< 0.005) continue;
                if(TMath::Abs(massGamma)<0.004) continue;
                if(TMath::Abs(0.77549 -massRho)<0.005) continue;
                if(TMath::Abs(1.01946 -massPhi)<0.005) continue;
                if(TMath::Abs(1.232 -massDelta)<0.005) continue;
                if(TMath::Abs(3.096 -massJPsi1)<0.005) continue;
                if(TMath::Abs(3.096 -massJPsi2)<0.005) continue;
                if(TMath::Abs(0.8955 -massKStar)<0.005) continue;
                if(TMath::Abs(1.86483 -massD0)<0.005) continue;
            }
            Int_t labelTrig = -2;
            Int_t labelAssoc =0;
            if(hh&&fAnalysisMC){
                labelTrig=trig->GetLabel();
                labelAssoc=assoc->GetLabel();
                
            }else if(hh){
                labelTrig=trig->GetIDCh();
                labelAssoc=assoc->GetIDCh();
            }
            
            if(labelTrig==labelAssoc) continue;
            
            
            if(!fCorrelationsGen&&fRemoveLamhFromCascade&&(trig->WhichCandidate()==2||trig->WhichCandidate()==3)){
                massSigmaP=TMath::Sqrt(0.13957*0.13957+1.1156*1.1156+2*(trig->E()*assoc->E()-triggPt*assocPt-trig->Pz()*assoc->Pz()));
                massSigmaN=TMath::Sqrt(0.13957*0.13957+1.1156*1.1156+2*(trig->E()*assoc->E()-triggPt*assocPt-trig->Pz()*assoc->Pz()));
                massXiN=TMath::Sqrt(0.13957*0.13957+1.1156*1.1156+2*(trig->E()*assoc->E()-triggPt*assocPt-trig->Pz()*assoc->Pz()));
                massOmegaN=TMath::Sqrt(0.4936*0.4936+1.1156*1.1156+2*(trig->E()*assoc->E()-triggPt*assocPt-trig->Pz()*assoc->Pz()));
                if(TMath::Abs( 1.3872-massSigmaN)< 0.005||TMath::Abs( 1.3828-massSigmaP)<0.005||TMath::Abs( 1.32171-massXiN)<0.005||TMath::Abs( 1.67245-massOmegaN)<0.005) continue;

            }
            Int_t idbinassoc [4]; 
            if(!fCorrelationsGen&&assoc->WhichCandidate()>4){
                idbinassoc[0]=fHistEffCorrectionK0->GetAxis(0)->FindBin(assocPt);
                idbinassoc[1]=fHistEffCorrectionK0->GetAxis(1)->FindBin(asocEta);
                idbinassoc[2]=fHistEffCorrectionK0->GetAxis(2)->FindBin(assocPhi);
                idbinassoc[3]=fHistEffCorrectionK0->GetAxis(3)->FindBin(fPV[2]);
            }else if(!fCorrelationsGen&&assoc->WhichCandidate()==4){
                idbinassoc[0]=fHistEffCorrectionHadron->GetAxis(0)->FindBin(assocPt);
                idbinassoc[1]=fHistEffCorrectionHadron->GetAxis(1)->FindBin(asocEta);
                idbinassoc[2]=fHistEffCorrectionHadron->GetAxis(2)->FindBin(assocPhi);
                idbinassoc[3]=fHistEffCorrectionHadron->GetAxis(3)->FindBin(fPV[2]);
            }        
            if(!hV0) {
                if(!fCorrelationsGen){
                    if(trig->WhichCandidate()==1) weight = fHistEffCorrectionK0->GetBinContent(idbintrigg)*(fHistEffCorrectionHadron->GetBinContent(idbinassoc)/(1. - fHistSecondaryCont->Interpolate(assocPt)));
                    else if (trig->WhichCandidate()==2) weight = fHistEffCorrectionLam->GetBinContent(idbintrigg)*(fHistEffCorrectionHadron->GetBinContent(idbinassoc)/(1. - fHistSecondaryCont->Interpolate(assocPt)));
                    else if (trig->WhichCandidate()==3) weight = fHistEffCorrectionAntiLam->GetBinContent(idbintrigg)*(fHistEffCorrectionHadron->GetBinContent(idbinassoc)/(1. - fHistSecondaryCont->Interpolate(assocPt)));
                    else weight = (fHistEffCorrectionHadron->GetBinContent(idbintrigg)/(1. - fHistSecondaryCont->Interpolate(triggPt))) * (fHistEffCorrectionHadron->GetBinContent(idbinassoc)/(1. - fHistSecondaryCont->Interpolate(assocPt)));

                    if(weight ==0) continue;
                }
                
                Double_t korel[9] = {triggPt,assocPt,deltaPhi,deltaEta, fPV[2],trig->WhichCandidate()-0.5,massTrig,perc,(Double_t)assocCharge};
                if(fCorrelationsGen&&fAnalysisMC) fHistMCKorelacie->Fill(korel);
                else if(fAnalysisMC) fHistKorelacieMCrec->Fill(korel,1./weight);
                else fHistKorelacie->Fill(korel,1./weight);
            }else{
                if(!fCorrelationsGen){
                    if(assoc->WhichCandidate()==5) weight = fHistEffCorrectionK0->GetBinContent(idbinassoc)*(fHistEffCorrectionHadron->GetBinContent(idbintrigg)/(1. - fHistSecondaryCont->Interpolate(triggPt)));
                    else if (assoc->WhichCandidate()==6) weight = fHistEffCorrectionLam->GetBinContent(idbinassoc)*(fHistEffCorrectionHadron->GetBinContent(idbintrigg)/(1. - fHistSecondaryCont->Interpolate(triggPt)));
                    else if (assoc->WhichCandidate()==7) weight = fHistEffCorrectionAntiLam->GetBinContent(idbinassoc)*(fHistEffCorrectionHadron->GetBinContent(idbintrigg)/(1. - fHistSecondaryCont->Interpolate(triggPt)));
                    if(weight ==0) continue;
                }
                
                massTrig = assoc->M();
                Double_t korel[9] = {triggPt,assocPt,deltaPhi,deltaEta, fPV[2],assoc->WhichCandidate()-0.5,massTrig,perc,(Double_t)trig->Charge()};
                if(fCorrelationsGen&&fAnalysisMC) fHistMCKorelacie->Fill(korel);
                else if(fAnalysisMC) fHistKorelacieMCrec->Fill(korel,1./weight);
                else fHistKorelacie->Fill(korel,1./weight);
            }
            
        }
    }
    

}
//____________________________________________________________________________
void AliAnalysisTaskDiHadCorrelHighPt::CorrelationsXi(TObjArray *triggers,TObjArray *associated,Double_t perc,Bool_t Xih){

    const Double_t kPi = TMath::Pi();
    Int_t nAssoc =  associated->GetEntriesFast();
    Int_t nTrig = triggers->GetEntriesFast();

    Double_t triggPt,triggEta,asocEta,assocPhi,deltaEta,deltaPhi,assocPt,weight,triggPhi;
    Double_t massTrig =-1;

    AliV0ChParticle* trig = 0x0;
    AliV0ChParticle* assoc = 0x0;
    for (Int_t i=0; i<nTrig; i++){
        weight = 1.;
        trig = (AliV0ChParticle*)  triggers->At(i);
        triggPt = trig->Pt();
        triggEta = trig->Eta();
        triggPhi = trig->Phi();
        if(triggPt<fPtTrigMin) continue;
        if(trig->WhichCandidate()<4) continue;
        if(Xih) massTrig=trig->M();

        Int_t idbintrigg [4];
        if(!fCorrelationsGen&&Xih){
            idbintrigg [0] =fHistEffCorrectionNegXi->GetAxis(0)->FindBin(triggPt);
            idbintrigg [1] =fHistEffCorrectionNegXi->GetAxis(1)->FindBin(triggEta);
            idbintrigg [2] =fHistEffCorrectionNegXi->GetAxis(2)->FindBin(triggPhi);
            idbintrigg [3] =fHistEffCorrectionNegXi->GetAxis(3)->FindBin(fPV[2]);
        }else if(!fCorrelationsGen){
            idbintrigg [0] =fHistEffCorrectionHadron->GetAxis(0)->FindBin(triggPt);
            idbintrigg [1] =fHistEffCorrectionHadron->GetAxis(1)->FindBin(triggEta);
            idbintrigg [2] =fHistEffCorrectionHadron->GetAxis(2)->FindBin(triggPhi);
            idbintrigg [3] =fHistEffCorrectionHadron->GetAxis(3)->FindBin(fPV[2]);
        }
        
        if(Xih&&!fMixing){
            if(!fCorrelationsGen){
                if(trig->WhichCandidate()==8) weight = fHistEffCorrectionNegXi->GetBinContent(idbintrigg);
                else weight = fHistEffCorrectionPosXi->GetBinContent(idbintrigg);
                if(weight==0) continue;   
            }
            Double_t triggers[5]={triggPt,fPV[2],trig->WhichCandidate()-2.5,massTrig,perc};
            if(fCorrelationsGen&&fAnalysisMC) fHistNumberOfTriggersGen->Fill(triggers); 
            else if(fAnalysisMC&&fCorrelations) fHistNumberOfTriggersRec->Fill(triggers,1./weight);
            else fHistNumberOfTriggers->Fill(triggers,1./weight);
        }else if(!fMixing){
            if(!fCorrelationsGen){
                weight = fHistEffCorrectionHadron->GetBinContent(idbintrigg)/(1. - fHistSecondaryCont->Interpolate(triggPt));
                if(weight==0) continue;
            }
            Double_t triggers[5]={triggPt,fPV[2],4.5,massTrig,perc};
            if(fCorrelationsGen&&fAnalysisMC) fHistNumberOfTriggersGen->Fill(triggers); 
            else if(fAnalysisMC&&fCorrelations) fHistNumberOfTriggersRec->Fill(triggers,1./weight);
            else fHistNumberOfTriggers->Fill(triggers,1./weight);

        }
        if(fMixing&&!Xih){
            if(TMath::Abs(trig->Eta())>0.8) continue;
        }

        for (Int_t j=0; j<nAssoc; j++){
            weight = 1.;
            assoc = (AliV0ChParticle*)  associated->At(j);

            asocEta = assoc->Eta();
            assocPhi = assoc->Phi();

            deltaEta = triggEta - asocEta;
            deltaPhi = triggPhi - assocPhi;
            assocPt = assoc->Pt();
            if (deltaPhi > (1.5*kPi)) deltaPhi -= 2.0*kPi;
            if (deltaPhi < (-0.5*kPi)) deltaPhi += 2.0*kPi;

            if(triggPt<=assocPt) continue;
            
            if(fMixing&&Xih){ 
                if(TMath::Abs(assoc->Eta())>0.8) continue;
            }           
            //removing autocorrelations
                
            Int_t negID = -1;
            Int_t posID = -2;
            Int_t bachID = -3;
            Int_t atrID = -4;

            if(Xih&&!fMixing){
                negID = trig->GetIDNeg();
                posID = trig->GetIDPos();
                bachID = trig->GetIDBach();
                atrID = assoc->GetIDCh();
            }else if(!fMixing) {
                negID = assoc->GetIDNeg();
                posID = assoc->GetIDPos();
                bachID = assoc->GetIDBach();
                atrID = trig->GetIDCh();
            }

            if ((TMath::Abs(negID))==(TMath::Abs(atrID))) continue;
            if ((TMath::Abs(posID))==(TMath::Abs(atrID))) continue;
            if ((TMath::Abs(bachID))==(TMath::Abs(atrID))) continue;

            Int_t idbinassoc [4];
            
            if(Xih){
                if(!fCorrelationsGen){
                    idbinassoc[0]=fHistEffCorrectionHadron->GetAxis(0)->FindBin(assocPt);
                    idbinassoc[1]=fHistEffCorrectionHadron->GetAxis(1)->FindBin(asocEta);
                    idbinassoc[2]=fHistEffCorrectionHadron->GetAxis(2)->FindBin(assocPhi);
                    idbinassoc[3]=fHistEffCorrectionHadron->GetAxis(3)->FindBin(fPV[2]);

                    if(trig->WhichCandidate()==8) weight = fHistEffCorrectionNegXi->GetBinContent(idbintrigg)*(fHistEffCorrectionHadron->GetBinContent(idbinassoc)/(1. - fHistSecondaryCont->Interpolate(assocPt)));
                    else weight = fHistEffCorrectionPosXi->GetBinContent(idbintrigg)*(fHistEffCorrectionHadron->GetBinContent(idbinassoc)/(1. - fHistSecondaryCont->Interpolate(assocPt)));
                    if(weight==0) continue;
                }
                

                Double_t korel[9] = {triggPt,assocPt,deltaPhi,deltaEta, fPV[2],trig->WhichCandidate()-0.5,massTrig,perc,(Double_t)assoc->Charge()};
                if(fCorrelationsGen&&fAnalysisMC&&!fMixing) fHistMCKorelacie->Fill(korel);
                else if(fAnalysisMC&&!fMixing) fHistKorelacieMCrec->Fill(korel,1./weight);
                else if(!fMixing) fHistKorelacie->Fill(korel,1./weight);
                else if(fAnalysisMC&&fMixing) fHistMCMixingRec->Fill(korel);
                else if (fMixing) fHistdPhidEtaMix->Fill(korel);
                else if (fMixingGen) fHistMCMixingGen->Fill(korel);
            }else{
                if(!fCorrelationsGen){
                    idbinassoc[0]=fHistEffCorrectionNegXi->GetAxis(0)->FindBin(assocPt);
                    idbinassoc[1]=fHistEffCorrectionNegXi->GetAxis(1)->FindBin(asocEta);
                    idbinassoc[2]=fHistEffCorrectionNegXi->GetAxis(2)->FindBin(assocPhi);
                    idbinassoc[3]=fHistEffCorrectionNegXi->GetAxis(3)->FindBin(fPV[2]);

                    if(assoc->WhichCandidate()==10) weight = (fHistEffCorrectionHadron->GetBinContent(idbintrigg)/(1. - fHistSecondaryCont->Interpolate(triggPt)))*fHistEffCorrectionNegXi->GetBinContent(idbinassoc);
                    else weight = (fHistEffCorrectionHadron->GetBinContent(idbintrigg)/(1. - fHistSecondaryCont->Interpolate(triggPt)))*fHistEffCorrectionPosXi->GetBinContent(idbinassoc);
                    if(weight==0) continue;
                }
                massTrig = assoc->M();
                Double_t korel[9] = {triggPt,assocPt,deltaPhi,deltaEta, fPV[2],assoc->WhichCandidate()-0.5,massTrig,perc,(Double_t)trig->Charge()};
                if(fCorrelationsGen&&fAnalysisMC&&!fMixing) fHistMCKorelacie->Fill(korel);
                else if(fAnalysisMC&&!fMixing) fHistKorelacieMCrec->Fill(korel,1./weight);
                else if(!fMixing) fHistKorelacie->Fill(korel,1./weight);
                else if(fAnalysisMC&&fMixing) fHistMCMixingRec->Fill(korel);
                else if (fMixing) fHistdPhidEtaMix->Fill(korel);
                else if (fMixingGen) fHistMCMixingGen->Fill(korel);
            }
        }
    }

}

//____________________________________________________________________________
void AliAnalysisTaskDiHadCorrelHighPt::CorelationsMixing(TObjArray *triggers, TObjArray *bgTracks, Float_t perc){

    const Double_t kPi = TMath::Pi();
    Int_t nAssoc = bgTracks->GetEntriesFast();
    Int_t nTrig = triggers->GetEntriesFast();
    Double_t assocCharge =0.;
    Double_t asocEta =0.;
    Double_t assocPhi =0.;
    Double_t assocPt =0.;
    Double_t triggPt, triggEta, triggPhi;

    AliV0ChParticle* trig = 0x0;
    AliVTrack* assoc = 0x0;
    Int_t idbinassoc [4];
    Int_t idbintrigg [4];
    Double_t weight=1;

    for (Int_t i=0; i<nTrig; i++){
        trig = (AliV0ChParticle*)  triggers->At(i);
        
        triggPt = trig->Pt(); 
        triggEta = trig->Eta();
        triggPhi = trig->Phi();

        Double_t massTrig = 0.;

        if(fMixCorrect&&trig->WhichCandidate()<4){
            idbintrigg[0]=fHistEffCorrectionK0->GetAxis(0)->FindBin(triggPt);
            idbintrigg[1]=fHistEffCorrectionK0->GetAxis(1)->FindBin(triggEta);
            idbintrigg[2]=fHistEffCorrectionK0->GetAxis(2)->FindBin(triggPhi);
            idbintrigg[3]=fHistEffCorrectionK0->GetAxis(3)->FindBin(fPV[2]);
        }else if(fMixCorrect){
            idbintrigg[0]=fHistEffCorrectionHadron->GetAxis(0)->FindBin(triggPt);
            idbintrigg[1]=fHistEffCorrectionHadron->GetAxis(1)->FindBin(triggEta);
            idbintrigg[2]=fHistEffCorrectionHadron->GetAxis(2)->FindBin(triggPhi);
            idbintrigg[3]=fHistEffCorrectionHadron->GetAxis(3)->FindBin(fPV[2]);
        } 
        if(trig->WhichCandidate()<4) massTrig=trig->M();
        for (Int_t j=0; j<nAssoc; j++){
             assoc = (AliVTrack*) bgTracks->At(j);
             if(!assoc) continue;

            if(isnan(assoc->Eta() ))continue;
            if(TMath::Abs(assoc->Eta())>0.8) continue;

             assocCharge = assoc->Charge();
             asocEta = assoc->Eta();
             assocPhi = assoc -> Phi();
             assocPt = assoc->Pt();

             if (( assocPt>=trig->Pt() ) || ( assocPt<fPtAsocMin )) continue;

             Double_t   deltaEta = triggEta - asocEta;
             Double_t   deltaPhi = triggPhi - assocPhi;
             

            if (deltaPhi > (1.5*kPi)) deltaPhi -= 2.0*kPi;
            if (deltaPhi < (-0.5*kPi)) deltaPhi += 2.0*kPi;

            if(fMixCorrect){
            	idbinassoc[0]=fHistEffCorrectionHadron->GetAxis(0)->FindBin(assocPt);
                idbinassoc[1]=fHistEffCorrectionHadron->GetAxis(1)->FindBin(asocEta);
                idbinassoc[2]=fHistEffCorrectionHadron->GetAxis(2)->FindBin(assocPhi);
                idbinassoc[3]=fHistEffCorrectionHadron->GetAxis(3)->FindBin(fPV[2]);

                if(trig->WhichCandidate()==1) weight = fHistEffCorrectionK0->GetBinContent(idbintrigg)*(fHistEffCorrectionHadron->GetBinContent(idbinassoc)/(1. - fHistSecondaryCont->Interpolate(assocPt)));
                else if (trig->WhichCandidate()==2) weight = fHistEffCorrectionLam->GetBinContent(idbintrigg)*(fHistEffCorrectionHadron->GetBinContent(idbinassoc)/(1. - fHistSecondaryCont->Interpolate(assocPt)));
                else if (trig->WhichCandidate()==3) weight = fHistEffCorrectionAntiLam->GetBinContent(idbintrigg)*(fHistEffCorrectionHadron->GetBinContent(idbinassoc)/(1. - fHistSecondaryCont->Interpolate(assocPt)));
                else weight = (fHistEffCorrectionHadron->GetBinContent(idbintrigg)/(1. - fHistSecondaryCont->Interpolate(triggPt))) * (fHistEffCorrectionHadron->GetBinContent(idbinassoc)/(1. - fHistSecondaryCont->Interpolate(assocPt)));

               	if(weight ==0) continue;
            }
            
            
            Double_t korel[9] = {triggPt,assocPt,deltaPhi,deltaEta, fPV[2],trig->WhichCandidate()-0.5,massTrig,perc,(Double_t)assocCharge};
            if(fAnalysisMC&&fMixingGen) fHistMCMixingGen->Fill(korel);
            else if(fAnalysisMC) fHistMCMixingRec->Fill(korel,weight);
            else fHistdPhidEtaMix->Fill(korel,weight);
        }
    }


}
//_____________________________________________________________
void AliAnalysisTaskDiHadCorrelHighPt::CorelationsMixinghV0(TObjArray *bgTracks, TObjArray *assocArray, Float_t perc){
    
    const Double_t kPi = TMath::Pi();
    Int_t nAssoc = assocArray->GetEntriesFast();
    Int_t nTrig = bgTracks->GetEntriesFast();
    Double_t triggPt, triggEta, triggPhi, assocPt, asocEta, assocPhi;
    Int_t idbinassoc [4];
    Int_t idbintrigg [4];
    Double_t weight=1;

    for (Int_t i=0; i<nTrig; i++){
        AliVTrack* trig = (AliVTrack*)  bgTracks->At(i);
        if(trig->Pt()<fPtTrigMin) continue;

        triggPt = trig->Pt(); 
        triggEta = trig->Eta();
        triggPhi = trig->Phi();

        if(fMixCorrect){
            idbintrigg[0]=fHistEffCorrectionHadron->GetAxis(0)->FindBin(triggPt);
            idbintrigg[1]=fHistEffCorrectionHadron->GetAxis(1)->FindBin(triggEta);
            idbintrigg[2]=fHistEffCorrectionHadron->GetAxis(2)->FindBin(triggPhi);
            idbintrigg[3]=fHistEffCorrectionHadron->GetAxis(3)->FindBin(fPV[2]);
        }

        if(TMath::Abs(trig->Eta())>0.8) continue;
         
        for (Int_t j=0; j<nAssoc; j++){
           
            AliV0ChParticle* assoc = (AliV0ChParticle*) assocArray->At(j);
            
            Double_t massAssoc = assoc->M();

            assocPt = assoc->Pt(); 
	        asocEta = assoc->Eta();
	        assocPhi = assoc->Phi();
            
            if (( assocPt>=triggPt ) || ( assocPt<fPtAsocMin )) continue;
            
            Double_t   deltaEta = triggEta - asocEta;
            Double_t   deltaPhi = triggPhi - assocPhi;
            
            if (deltaPhi > (1.5*kPi)) deltaPhi -= 2.0*kPi;
            if (deltaPhi < (-0.5*kPi)) deltaPhi += 2.0*kPi;

            if(fMixCorrect){
            	if(assoc->WhichCandidate()>4){
		            idbinassoc[0]=fHistEffCorrectionK0->GetAxis(0)->FindBin(assocPt);
		            idbinassoc[1]=fHistEffCorrectionK0->GetAxis(1)->FindBin(asocEta);
		            idbinassoc[2]=fHistEffCorrectionK0->GetAxis(2)->FindBin(assocPhi);
		            idbinassoc[3]=fHistEffCorrectionK0->GetAxis(3)->FindBin(fPV[2]);
		        }else{
		        	idbinassoc[0]=fHistEffCorrectionHadron->GetAxis(0)->FindBin(assocPt);
                	idbinassoc[1]=fHistEffCorrectionHadron->GetAxis(1)->FindBin(asocEta);
                	idbinassoc[2]=fHistEffCorrectionHadron->GetAxis(2)->FindBin(assocPhi);
                	idbinassoc[3]=fHistEffCorrectionHadron->GetAxis(3)->FindBin(fPV[2]);
		        }

                if(assoc->WhichCandidate()==5) weight = fHistEffCorrectionK0->GetBinContent(idbinassoc)*(fHistEffCorrectionHadron->GetBinContent(idbintrigg)/(1. - fHistSecondaryCont->Interpolate(triggPt)));
                else if (assoc->WhichCandidate()==6) weight = fHistEffCorrectionLam->GetBinContent(idbinassoc)*(fHistEffCorrectionHadron->GetBinContent(idbintrigg)/(1. - fHistSecondaryCont->Interpolate(triggPt)));
                else if (assoc->WhichCandidate()==3) weight = fHistEffCorrectionAntiLam->GetBinContent(idbinassoc)*(fHistEffCorrectionHadron->GetBinContent(idbintrigg)/(1. - fHistSecondaryCont->Interpolate(triggPt)));
               
               	if(weight ==0) continue;
            }

            Double_t korel[9] = {triggPt,assocPt,deltaPhi,deltaEta, fPV[2],assoc->WhichCandidate()-0.5,massAssoc,perc,(Double_t)trig->Charge()};
            if(fAnalysisMC&&fMixingGen) fHistMCMixingGen->Fill(korel);
            else if(fAnalysisMC) fHistMCMixingRec->Fill(korel,weight);
            else fHistdPhidEtaMix->Fill(korel,weight);
        }
    }
}
//_____________________________________________________________
void AliAnalysisTaskDiHadCorrelHighPt::TopologCuts(Double_t pttrig,Double_t mass,Double_t dcaNeg, Double_t dcaPos,Double_t dcaDau, Double_t V0rad, Double_t cosPA,Double_t lifetime,Double_t massSell,Double_t triggType,Double_t status){
    
    Double_t topolCutsValues[11]={pttrig,mass,dcaNeg,dcaPos,dcaDau,V0rad,cosPA,lifetime,massSell,triggType,status};
    if(fAnalysisMC) fHistTopolCutMC->Fill(topolCutsValues);
    else fHistTopolCut->Fill(topolCutsValues);
    
}
//____________________________________________________________
void AliAnalysisTaskDiHadCorrelHighPt::FillMC(const AliVParticle *V0,Int_t pdgV0,Int_t pdgDau1, Int_t pdgDau2,Int_t triggerType, Double_t mass, Double_t posTrackProp[6],Double_t negTrackProp[6]){
    
    if(fPurityCheck){
        Double_t purity[6] ={V0->Pt(),mass,triggerType-0.5,0.5,-1,V0->Eta()};
        fHistPurityCheck->Fill(purity);
    }
    
    if((fCorrelations||fMixing)&&!fPureV0) { // for MC closure test - also misidentified V0 taken
       if(V0->Pt()>fPtTrigMin) fselectedMCV0Triggersrec-> Add(new AliV0ChParticle(V0->Eta(), V0->Phi(), V0->Pt(), triggerType,0,posTrackProp[4],negTrackProp[4],mass,posTrackProp[0],posTrackProp[1],posTrackProp[2],(Int_t)posTrackProp[3],negTrackProp[0],negTrackProp[1],negTrackProp[2],(Int_t)negTrackProp[3])); // all reconstructed candidates for raw correlation function, with reconstructed pt
       if(V0->Pt()>fPtAsocMin) fselectedMCV0assoc -> Add(new AliV0ChParticle(V0->Eta(), V0->Phi(), V0->Pt(), triggerType+4,0,posTrackProp[4],negTrackProp[4],mass,posTrackProp[0],posTrackProp[1],posTrackProp[2],(Int_t)posTrackProp[3],negTrackProp[0],negTrackProp[1],negTrackProp[2],(Int_t)negTrackProp[3])); // all reconstructed candidates for raw correlation function, with reconstructed pt
    }
    if(posTrackProp[5]!=-1) posTrackProp[5]= TMath::Abs(posTrackProp[5]);
    if(negTrackProp[5]!=-1) negTrackProp[5]= TMath::Abs(negTrackProp[5]);
    AliMCParticle *mcPosTrack = (AliMCParticle*)fmcEvent->GetTrack((Int_t)posTrackProp[5]);
    if (!mcPosTrack) return;
    Int_t PosTrackPdg = mcPosTrack->PdgCode();
    AliMCParticle *mcNegTrack = (AliMCParticle*)fmcEvent->GetTrack((Int_t)negTrackProp[5]);
    if (!mcNegTrack) return;
    Int_t NegTrackPdg = mcNegTrack->PdgCode();
    
    if(fPurityCheck){
        Double_t puri[6] ={V0->Pt(),mass,triggerType-0.5,1.5,-1,V0->Eta()};
        fHistPurityCheck->Fill(puri);
    }

    Int_t myTrackPosMotherLabel = IsGoodMCV0(mcPosTrack,mcNegTrack);
    if(myTrackPosMotherLabel < 0) return;

    if(fPurityCheck){
        Double_t pu[6] ={V0->Pt(),mass,triggerType-0.5,3.5,-1,V0->Eta()};
        fHistPurityCheck->Fill(pu);
    }
    
    AliMCParticle *mcPosMother = (AliMCParticle*)fmcEvent->GetTrack(myTrackPosMotherLabel);
    if (!mcPosMother) return;

    Int_t MotherPdg = mcPosMother->PdgCode();
    Int_t MotherOfMotherLabel = mcPosMother->GetMother();
    
    Bool_t IsPhysPrim = mcPosMother->IsPhysicalPrimary();
    
    Bool_t IsFromMC = (MotherPdg==pdgV0)&&(PosTrackPdg==pdgDau1)&&(NegTrackPdg==pdgDau2);

    Bool_t IsFromCascade = kFALSE;
    Bool_t IsFromXi = kFALSE;
    
    Int_t MoMPdg = 50;
    Double_t xiPt =0;
    if (MotherOfMotherLabel != -1)
    {
        AliMCParticle *mcPosMotherOfMother = (AliMCParticle*)fmcEvent->GetTrack(MotherOfMotherLabel);
        Int_t MotherOfMotherPdg = mcPosMotherOfMother->PdgCode();
        MoMPdg = TMath::Abs(MotherOfMotherPdg);
        if(fAacceptLambdasFromCasscade) IsFromCascade = (((MoMPdg == 3222)|| (MoMPdg==3212)|| (MoMPdg==3112) || (MoMPdg==3224) || (MoMPdg==3214) || (MoMPdg==3114) || (MoMPdg==3322) || (MoMPdg==3312)|| (MoMPdg==3324) || (MoMPdg==3314) || (MoMPdg==3334)) && (mcPosMotherOfMother->IsPhysicalPrimary()));
        IsFromXi = MoMPdg==3322 ||TMath::Abs(MoMPdg)==3312; 
        if(IsFromXi) xiPt = mcPosMotherOfMother->Pt();
    }
    if(fPurityCheck){
        Double_t purit[6] ={V0->Pt(),mass,triggerType-0.5,4.5,-1,V0->Eta()};
        fHistPurityCheck->Fill(purit);
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
        fHistPurityCheck->Fill(purity);
    }
    
    if(isGoodID&&fPurityCheck){
        Double_t purity[6] ={V0->Pt(),mass,triggerType-0.5,6.5,-1,V0->Eta()};
        fHistPurityCheck->Fill(purity);
    }

    if(!isFromDecay&&!isFromMaterial&&isGoodID&&fPurityCheck){
        Double_t purity[6] ={V0->Pt(),mass,triggerType-0.5,7.5,-1,V0->Eta()};
        fHistPurityCheck->Fill(purity);
    }
    
    Bool_t IsParticleFromMC = kFALSE;
    
    if(pdgV0==310) IsParticleFromMC= (IsFromMC&&IsPhysPrim);
    else {
        IsParticleFromMC = (IsFromMC&&(IsPhysPrim||IsFromCascade));
    }
    
    Double_t V0mcPt = mcPosMother->Pt();

    if(IsFromXi){
        fHistLambdaFeedDown->Fill(V0->Pt(),xiPt,triggerType-1.5);
    }
   
    if(IsParticleFromMC){
        if(fCorrelations&&fPureV0) { // for MC closure test - only good ID V0 taken
            if(V0->Pt()>fPtTrigMin) fselectedMCV0Triggersrec-> Add(new AliV0ChParticle(V0->Eta(), V0->Phi(), V0->Pt(), triggerType,0,(Int_t)posTrackProp[4],(Int_t)negTrackProp[4],mass,posTrackProp[0],posTrackProp[1],posTrackProp[2],(Int_t)posTrackProp[3],negTrackProp[0],negTrackProp[1],negTrackProp[2],(Int_t)negTrackProp[3])); // all reconstructed candidates for raw correlation function, with reconstructed pt
            if(V0->Pt()>fPtAsocMin) fselectedMCV0assoc -> Add(new AliV0ChParticle(V0->Eta(), V0->Phi(), V0->Pt(), triggerType+4,0,(Int_t)posTrackProp[4],(Int_t)negTrackProp[4],mass,posTrackProp[0],posTrackProp[1],posTrackProp[2],(Int_t)posTrackProp[3],negTrackProp[0],negTrackProp[1],negTrackProp[2],(Int_t)negTrackProp[3])); // all reconstructed candidates for raw correlation function, with reconstructed pt
        }
        fHistPtResolution->Fill(V0mcPt,V0->Pt(),triggerType-0.5);
        Double_t V0mcEta = mcPosMother->Eta();
    
        if(fEfficiency){
            Double_t v0effic[6]={V0mcPt,fPV[2],triggerType-0.5,V0mcEta,mass,mcPosMother->Phi()};
            fHistRecV0->Fill(v0effic);
        }

        if(triggerType==1)fHistK0MassPtCut->Fill(mass,V0mcPt,9.5);
        if(triggerType==2)fHistLambdaMassPtCut->Fill(mass,V0mcPt,9.5);
        if(triggerType==3)fHistAntiLambdaMassPtCut->Fill(mass,V0mcPt,9.5);
    }

}
//____________________________________________________________________//
AliAODTrack * AliAnalysisTaskDiHadCorrelHighPt::SetAliAODTrack(Double_t theta,Double_t phi, Double_t pt , Short_t charge){
    AliAODTrack * track = new AliAODTrack();
    track->SetPhi(phi);
    track->SetTheta(theta);
    track->SetPt(pt);
    track->SetCharge(charge);
    return track; 
}
//____________________________________________________________________//
Bool_t AliAnalysisTaskDiHadCorrelHighPt::IsNotOOBPileUp(AliESDtrack* neg, AliESDtrack* pos, AliESDtrack* bach){
    if ((neg->Pt()<=2.2) && (pos->Pt()<=2.2) && (bach->Pt()<=2.2)) return kTRUE;
    else {
        return (neg->Pt()>2.2 && neg->GetTOFsignal()*1.e-3 <100) 
        || (pos->Pt()>2.2 && pos->GetTOFsignal()*1.e-3 <100) 
        || (bach->Pt()>2.2 && bach->GetTOFsignal()*1.e-3 <100) ? kTRUE : kFALSE;
    }

    return kTRUE;

}
//____________________________________________________________________//
Bool_t AliAnalysisTaskDiHadCorrelHighPt::IsMyGoodXiCandidate(AliESDcascade * cas, Double_t * par){

    Int_t chargeXi= cas->Charge(); // xi charge

    Double_t qa =0.;
    if(chargeXi>0) {
        par[0]= 7;
        cas->ChangeMassHypothesis(qa, -3312);
    }
    qa= 0.;
    if(chargeXi<0) {
        par[0]= 6;
        cas->ChangeMassHypothesis(qa, 3312);
    }

    Double_t cascadept = cas->Pt();
    Double_t massXi = cas->GetEffMassXi();// invariant mass
                

    if (( massXi< 1.29 || massXi > 1.35 ) ) return kFALSE; // invariant mass selection

    if(chargeXi>0) fHistXiPlusMassPtCut->Fill(massXi,cascadept,0.5);
    else fHistXiMinusMassPtCut->Fill(massXi,cascadept,0.5);

    if(TMath::Abs(cas->RapXi())>0.5) return kFALSE; // rapidity selection

    if(chargeXi>0) fHistXiPlusMassPtCut->Fill(massXi,cascadept,1.5);
    else fHistXiMinusMassPtCut->Fill(massXi,cascadept,1.5);

    if( cas->GetCascadeCosineOfPointingAngle( fPV[0],fPV[1],fPV[2] )<0.97) return kFALSE; // Cascade cos PA

    if(chargeXi>0) fHistXiPlusMassPtCut->Fill(massXi,cascadept,2.5);
    else fHistXiMinusMassPtCut->Fill(massXi,cascadept,2.5);

    if(!IsMyGoodXiDecayRadius(cas)) return kFALSE;

    if(chargeXi>0) fHistXiPlusMassPtCut->Fill(massXi,cascadept,3.5);
    else fHistXiMinusMassPtCut->Fill(massXi,cascadept,3.5);

    if(!IsMyGoodDaughterV0DecayRadius(cas)) return kFALSE;

    if(chargeXi>0) fHistXiPlusMassPtCut->Fill(massXi,cascadept,4.5);
    else fHistXiMinusMassPtCut->Fill(massXi,cascadept,4.5);

    if (cas->GetV0CosineOfPointingAngle( fPV[0],fPV[1],fPV[2] ) < 0.97) return kFALSE; // V0 PA

    if(chargeXi>0) fHistXiPlusMassPtCut->Fill(massXi,cascadept,5.5);
    else fHistXiMinusMassPtCut->Fill(massXi,cascadept,5.5);

    if (cas->GetD( fPV[0],fPV[1],fPV[2] ) < 0.06) return kFALSE; // DCA V0 to PV

    if(chargeXi>0) fHistXiPlusMassPtCut->Fill(massXi,cascadept,6.5);
    else fHistXiMinusMassPtCut->Fill(massXi,cascadept,6.5);

    if (cas->GetDcaV0Daughters() >1.5) return kFALSE; // DCA V0 daughters

    if(chargeXi>0) fHistXiPlusMassPtCut->Fill(massXi,cascadept,7.5);
    else fHistXiMinusMassPtCut->Fill(massXi,cascadept,7.5);

    if(cas->GetDcaXiDaughters() > 1.3) return kFALSE; // DCA bach-V0

    if(chargeXi>0) fHistXiPlusMassPtCut->Fill(massXi,cascadept,8.5);
    else fHistXiMinusMassPtCut->Fill(massXi,cascadept,8.5);

    if(!IsMyGoodXiDaughterV0Mass(cas,chargeXi)) return kFALSE;

    if(chargeXi>0) fHistXiPlusMassPtCut->Fill(massXi,cascadept,9.5);
    else fHistXiMinusMassPtCut->Fill(massXi,cascadept,9.5);

    if(!IsMyGoodXiProperLifetime(cas)) return kFALSE;

    if(chargeXi>0) fHistXiPlusMassPtCut->Fill(massXi,cascadept,10.5);
    else fHistXiMinusMassPtCut->Fill(massXi,cascadept,10.5);

    Int_t ids[3];
    if(!IsMyGoodXiDaughterTracks(cas,chargeXi,ids,massXi)) return kFALSE;

    if(chargeXi>0) fHistXiPlusMassPtCut->Fill(massXi,cascadept,18.5);
    else fHistXiMinusMassPtCut->Fill(massXi,cascadept,18.5);

    par[1] = massXi;
    par[2] = ids[0]; 
    par[3] = ids[1];
    par[4] = ids[2];

    return kTRUE;               

}
//__________________________________________________________//
Bool_t AliAnalysisTaskDiHadCorrelHighPt::IsMyGoodXiDecayRadius(AliESDcascade * cas){
    Double_t posXi[3] = { -1000.0, -1000.0, -1000.0 };

    cas->GetXYZcascade( posXi[0],  posXi[1], posXi[2] );
    if( TMath::Sqrt( posXi[0]*posXi[0]  +  posXi[1]*posXi[1] )<0.6 ) return kFALSE; // cascade decay radius

    return kTRUE;
}
//__________________________________________________________//
Bool_t AliAnalysisTaskDiHadCorrelHighPt::IsMyGoodDaughterV0DecayRadius(AliESDcascade * cas){
    Double_t posV0Xi[3] = { -1000.0, -1000.0, -1000.0 };
    cas->GetXYZ( posV0Xi[0],  posV0Xi[1], posV0Xi[2] );
    if( TMath::Sqrt( posV0Xi[0]*posV0Xi[0]  +  posV0Xi[1]*posV0Xi[1] ) < 1.2 )return kFALSE; // V0 decay radius

    return kTRUE;
}
//_________________________________________________________//
Bool_t AliAnalysisTaskDiHadCorrelHighPt::IsMyGoodXiProperLifetime(AliESDcascade * cas){
    Double_t xiMom[3];
    Double_t posXi[3] = { -1000.0, -1000.0, -1000.0 };

    cas->GetXYZcascade( posXi[0],  posXi[1], posXi[2] );
    cas->GetPxPyPz( xiMom[0], xiMom[1], xiMom[2] );
    Double_t totalMomXi  = TMath::Sqrt( xiMom[0]*xiMom[0]   + xiMom[1]*xiMom[1]   + xiMom[2]*xiMom[2] );
                    
    Double_t ditanceToMomentum = TMath::Sqrt(TMath::Power( posXi[0] - fPV[0] , 2) +TMath::Power( posXi[1] - fPV[1] , 2) +TMath::Power( posXi[2] - fPV[2] , 2));
    ditanceToMomentum /= (totalMomXi+1e-13);

    if(ditanceToMomentum* 1.32171 > 3*4.91) return kFALSE;  // proper lifetime cut

    return kTRUE;

}
//________________________________________________________//
Bool_t AliAnalysisTaskDiHadCorrelHighPt::IsMyGoodXiDaughterTracks(AliESDcascade * cas,Int_t casCharge, Int_t ids[3], Double_t massXi){
    UInt_t idxPosXi = (UInt_t) TMath::Abs( cas->GetPindex() );
    UInt_t idxNegXi = (UInt_t) TMath::Abs( cas->GetNindex() );
    UInt_t idBachXi = (UInt_t) TMath::Abs( cas->GetBindex() );

    AliESDtrack * pTrackXi = fESD->GetTrack( idxPosXi );
    AliESDtrack * nTrackXi = fESD->GetTrack( idxNegXi );
    AliESDtrack * bachTrackXi = fESD->GetTrack( idBachXi );

    if (!pTrackXi || !nTrackXi || !bachTrackXi ) {
        AliWarning("ERROR: Could not retrieve one of the 3 ESD daughter tracks of the cascade ...");
        return kFALSE;
    }

    if(TMath::Abs(nTrackXi->Eta())>0.8 || TMath::Abs(pTrackXi->Eta())>0.8 || TMath::Abs(bachTrackXi->Eta())>0.8 )  return kFALSE; // daughter pseudorapidity cut

    Double_t cascadept = cas->Pt();

    if(casCharge>0) fHistXiPlusMassPtCut->Fill(massXi,cascadept,11.5);
    else fHistXiMinusMassPtCut->Fill(massXi,cascadept,11.5);

    if(!fAnalysisMC){
        if(!IsMyGoodXiDaughterTrackPID(nTrackXi,casCharge,-1)) return kFALSE;
        if(!IsMyGoodXiDaughterTrackPID(pTrackXi,casCharge,1)) return kFALSE; // nSigma daughter tracks
        if( fPIDResponse->NumberOfSigmasTPC( bachTrackXi, AliPID::kPion ) > 5 ) return kFALSE;
    }
    if(casCharge>0) fHistXiPlusMassPtCut->Fill(massXi,cascadept,12.5);
    else fHistXiMinusMassPtCut->Fill(massXi,cascadept,12.5);
                    
    if( pTrackXi->GetTPCNcls() < 70 || nTrackXi->GetTPCNcls() < 70 || bachTrackXi->GetTPCNcls() < 70 ) return kFALSE; // number of clusters daughter tracks

    if(casCharge>0) fHistXiPlusMassPtCut->Fill(massXi,cascadept,13.5);
    else fHistXiMinusMassPtCut->Fill(massXi,cascadept,13.5);

    ULong_t pStat = pTrackXi->GetStatus();
    ULong_t nStat = nTrackXi->GetStatus();
    ULong_t bStat = bachTrackXi->GetStatus();

    if ((pStat&AliESDtrack::kTPCrefit) == 0) return kFALSE;
    if ((nStat&AliESDtrack::kTPCrefit) == 0) return kFALSE;
    if ((bStat&AliESDtrack::kTPCrefit) == 0) return kFALSE; // TPC refit daughters tracks

    if(casCharge>0) fHistXiPlusMassPtCut->Fill(massXi,cascadept,14.5);
    else fHistXiMinusMassPtCut->Fill(massXi,cascadept,14.5);

    if(!IsMyGoodBachelorTrack(bachTrackXi)) return kFALSE;                                

    if(casCharge>0) fHistXiPlusMassPtCut->Fill(massXi,cascadept,16.5);
    else fHistXiMinusMassPtCut->Fill(massXi,cascadept,16.5);


    if(casCharge<0){
        if( TMath::Abs( pTrackXi  ->GetD( fPV[0], fPV[1],fMagneticField  ) ) < 0.03) return kFALSE;  // DCA baryon V0 track to PV
        if( TMath::Abs( nTrackXi  ->GetD( fPV[0], fPV[1],fMagneticField  ) ) < 0.04) return kFALSE; // DCA meson V0 track to PV
    }
    if(casCharge>0){
        if( TMath::Abs( pTrackXi  ->GetD( fPV[0], fPV[1],fMagneticField  ) ) < 0.04) return kFALSE;  // DCA meson V0 track to PV
        if( TMath::Abs( nTrackXi  ->GetD( fPV[0], fPV[1],fMagneticField  ) ) < 0.03) return kFALSE; // DCA baryon V0 track to PV
    }

    if(casCharge>0) fHistXiPlusMassPtCut->Fill(massXi,cascadept,17.5);
    else fHistXiMinusMassPtCut->Fill(massXi,cascadept,17.5);

    if(!IsNotOOBPileUp(nTrackXi,pTrackXi,bachTrackXi)) return kFALSE;   

    ids[0] = pTrackXi->GetID();
    ids[1] = nTrackXi->GetID();
    ids[2] = bachTrackXi->GetID();

    return kTRUE;
}
//________________________________________________________//
Bool_t AliAnalysisTaskDiHadCorrelHighPt::IsMyGoodXiDaughterTrackPID(AliESDtrack * tr, Int_t casCharge, Int_t trCharge){

    if(casCharge<0  && trCharge<0 && fPIDResponse->NumberOfSigmasTPC( tr, AliPID::kPion ) > 5 ) return kFALSE;
    if(casCharge>0 && trCharge<0 && fPIDResponse->NumberOfSigmasTPC( tr, AliPID::kProton ) > 5 ) return kFALSE;

    if(casCharge>0 && trCharge>0 && fPIDResponse->NumberOfSigmasTPC( tr, AliPID::kPion ) > 5 ) return kFALSE;
    if(casCharge<0 && trCharge>0 && fPIDResponse->NumberOfSigmasTPC( tr, AliPID::kProton ) > 5 ) return kFALSE;

    return kTRUE;
}
//________________________________________________________//
Bool_t AliAnalysisTaskDiHadCorrelHighPt::IsMyGoodBachelorTrack(AliESDtrack * bach){

    if( TMath::Abs( bach->GetD( fPV[0], fPV[1],fMagneticField  ) ) < 0.04) return kFALSE; // DCA bachelor - PV

    Float_t dca[2];
    bach->GetDZ(fPV[0], fPV[1],fPV[2],fMagneticField,dca);
    if( TMath::Abs(dca[1] ) >4) return kFALSE; // DCAz bach to PV - pile-up cut 

    return kTRUE;

}
//________________________________________________________//
Bool_t AliAnalysisTaskDiHadCorrelHighPt::IsMyGoodXiDaughterV0Mass(AliESDcascade * cas, Int_t casChar){

    Double_t lPMom[3];
    Double_t lNMom[3];

    cas->GetPPxPyPz( lPMom[0], lPMom[1], lPMom[2] );
    cas->GetNPxPyPz( lNMom[0], lNMom[1], lNMom[2] );

    Double_t m1,m2,e12,e22;

    if(casChar<0){
        //+-+ Recalculate Lambda mass from scratch
        //Under Lambda hypothesis, the positive daughter is the proton, negative pion
        m1 = TDatabasePDG::Instance()->GetParticle(kProton)->Mass();
        m2 = TDatabasePDG::Instance()->GetParticle(kPiPlus)->Mass();
        e12 = TMath::Sqrt(m1*m1+lPMom[0]*lPMom[0]+lPMom[1]*lPMom[1]+lPMom[2]*lPMom[2]);
        e22 = TMath::Sqrt(m2*m2+lNMom[0]*lNMom[0]+lNMom[1]*lNMom[1]+lNMom[2]*lNMom[2]);
    }else if(casChar>0){
        //+-+ Recalculate AntiLambda mass from scratch
        //Under AntiLambda hypothesis, the positive daughter is the pion, negative antiproton
        m1 = TDatabasePDG::Instance()->GetParticle(kPiPlus)->Mass();
        m2 = TDatabasePDG::Instance()->GetParticle(kProton)->Mass();
        e12   = TMath::Sqrt(m1*m1+lPMom[0]*lPMom[0]+lPMom[1]*lPMom[1]+lPMom[2]*lPMom[2]);
        e22   = TMath::Sqrt(m2*m2+lNMom[0]*lNMom[0]+lNMom[1]*lNMom[1]+lNMom[2]*lNMom[2]);
    }else {
        m1 = -1;
        m2 = -1;
        e12 = -1;
        e22 = -1;
    }
    Double_t invMassSquared =  m1*m1+m2*m2+2.*(e12*e22-lPMom[0]*lNMom[0]-lPMom[1]*lNMom[1]-lPMom[2]*lNMom[2]);
    if( TMath::Abs(TMath::Sqrt(TMath::Max(invMassSquared,0.)) - 1.115683) > 0.008) return kFALSE;

    return kTRUE;

}
//________________________________________________________//
Int_t AliAnalysisTaskDiHadCorrelHighPt::IsGoodMCV0(AliMCParticle *pos, AliMCParticle * neg){
    Int_t myTrackPosMotherLabel = pos->GetMother();
    Int_t myTrackNegMotherLabel = neg->GetMother();

    if ((myTrackPosMotherLabel==-1)||(myTrackNegMotherLabel==-1)) return -200;
    if (myTrackPosMotherLabel!=myTrackNegMotherLabel) return -200;

    return myTrackPosMotherLabel;
}

//________________________________________________________//
Int_t AliAnalysisTaskDiHadCorrelHighPt::IsGoodMCCascade(AliMCParticle *v0, AliMCParticle * bach){
    Int_t MotherOfV0Label = v0->GetMother();
    Int_t MotherOfBachelorLabel = (Int_t) TMath::Abs(bach->GetMother());

    if ((MotherOfV0Label==-1)||(MotherOfBachelorLabel==-1)) return -200;
    if (MotherOfV0Label!=MotherOfBachelorLabel) return -200;

    return MotherOfBachelorLabel;
}