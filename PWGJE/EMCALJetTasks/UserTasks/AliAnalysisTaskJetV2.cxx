/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

/* 
 * Jet V2 task
 *
 * this task is part of the emcal jet framework and should be run in the emcaljet train
 * the following extensions to an accepted AliVEvent are expected:
 *      - (anti-kt) jets
 *      - background estimate rho
 *      - pico tracks
 *      aod's and esd's are handled transparently
 * the task will attempt to estimate a phi-dependent background density rho 
 * by fitting vn harmonics to the dpt/dphi distribution
 *
 * author: Redmer Alexander Bertens, Utrecht Univeristy, Utrecht, Netherlands
 * rbertens@cern.ch, rbertens@nikhef.nl, r.a.bertens@uu.nl 
 */

// root includes
#include <TStyle.h>
#include <TRandom3.h>
#include <TChain.h>
#include <TMath.h>
#include <TF1.h>
#include <TF2.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TFile.h>
// aliroot includes
#include <AliAnalysisTask.h>
#include <AliAnalysisManager.h>
#include <AliCentrality.h>
#include <AliVVertex.h>
#include <AliVTrack.h>
#include <AliVVZERO.h>
#include <AliESDEvent.h>
#include <AliAODEvent.h>
#include <AliAODTrack.h>
#include <AliOADBContainer.h>
// emcal jet framework includes
#include <AliPicoTrack.h>
#include <AliEmcalJet.h>
#include <AliRhoParameter.h>
#include <AliLocalRhoParameter.h>
#include <AliAnalysisTaskJetV2.h>
#include <AliClusterContainer.h>

class AliAnalysisTaskJetV2;
using namespace std;

ClassImp(AliAnalysisTaskJetV2)

AliAnalysisTaskJetV2::AliAnalysisTaskJetV2() : AliAnalysisTaskEmcalJet("AliAnalysisTaskJetV2", kTRUE), 
    fDebug(0), fRunToyMC(kFALSE), fLocalInit(0), fAttachToEvent(kTRUE), fFillHistograms(kTRUE), fFillQAHistograms(kTRUE), fReduceBinsXByFactor(-1.), fReduceBinsYByFactor(-1.), fNoEventWeightsForQC(kTRUE), fCentralityClasses(0), fExpectedRuns(0), fExpectedSemiGoodRuns(0), fUserSuppliedV2(0), fUserSuppliedV3(0), fUserSuppliedR2(0), fUserSuppliedR3(0), fTracksCont(0), fClusterCont(0), fJetsCont(0), fLeadingJet(0), fNAcceptedTracks(0), fNAcceptedTracksQCn(0), fFitModulationType(kNoFit), fFitGoodnessTest(kChi2Poisson), fQCRecovery(kTryFit), fUsePtWeight(kTRUE), fUsePtWeightErrorPropagation(kTRUE), fDetectorType(kVZEROComb), fAnalysisType( kCharged), fFitModulationOptions("QWLI"), fRunModeType(kGrid), fDataType(kESD), fCollisionType(kPbPb), fRandom(0), fRunNumber(-1), fMappedRunNumber(0), fInCentralitySelection(-1), fFitModulation(0), fFitControl(0), fMinPvalue(0.01), fMaxPvalue(1), fNameSmallRho(""), fCachedRho(0), fSoftTrackMinPt(0.15), fSoftTrackMaxPt(5.), fSemiGoodJetMinPhi(0.), fSemiGoodJetMaxPhi(4.), fSemiGoodTrackMinPhi(0.), fSemiGoodTrackMaxPhi(4.), fAbsVertexZ(10), fHistCentrality(0), fHistVertexz(0), fHistRunnumbersPhi(0), fHistRunnumbersEta(0), fHistPvalueCDFROOT(0), fHistPvalueCDFROOTCent(0), fHistChi2ROOTCent(0), fHistPChi2Root(0),  fHistPvalueCDF(0), fHistPvalueCDFCent(0), fHistChi2Cent(0), fHistPChi2(0), fHistKolmogorovTest(0), fHistKolmogorovTestCent(0), fHistPKolmogorov(0), fHistRhoStatusCent(0), fHistUndeterminedRunQA(0), fMinDisanceRCtoLJ(0), fMaxCones(-1), fExcludeLeadingJetsFromFit(1.), fRebinSwapHistoOnTheFly(kTRUE), fPercentageOfFits(10.), fOutputList(0), fOutputListGood(0), fOutputListBad(0), fHistAnalysisSummary(0), fHistSwap(0), fProfV2(0), fProfV2Cumulant(0), fProfV3(0), fProfV3Cumulant(0), fHistPsiControl(0), fHistPsiSpread(0), fHistPsiVZEROA(0), fHistPsiVZEROC(0), fHistPsiVZERO(0), fHistPsiTPC(0), fHistPsiVZEROAV0M(0), fHistPsiVZEROCV0M(0), fHistPsiVZEROVV0M(0), fHistPsiTPCiV0M(0), fHistPsiVZEROATRK(0), fHistPsiVZEROCTRK(0), fHistPsiVZEROTRK(0), fHistPsiTPCTRK(0), fHistRhoVsMult(0), fHistRhoVsCent(0), fHistRhoAVsMult(0), fHistRhoAVsCent(0), fVZEROgainEqualization(0x0), fVZEROgainEqualizationPerRing(kFALSE), fChi2A(0x0), fChi2C(0x0), fChi3A(0x0), fChi3C(0x0), fOADB(0x0)
{
    for(Int_t i(0); i < 10; i++) {
        fProfV2Resolution[i] = 0;
        fProfV3Resolution[i] = 0;
        fHistPicoTrackPt[i] = 0;
        fHistPicoTrackMult[i] = 0;
        fHistPicoCat1[i] = 0;
        fHistPicoCat2[i] = 0;
        fHistPicoCat3[i] = 0;
        fHistClusterPt[i] = 0;
        fHistClusterEtaPhi[i] = 0;
        fHistClusterEtaPhiWeighted[i] = 0;
        fHistRhoPackage[i] = 0;
        fHistRho[i] = 0;
        fHistRCPhiEta[i] = 0;
        fHistRhoVsRCPt[i] = 0;
        fHistRCPt[i] = 0;
        fHistDeltaPtDeltaPhi2[i] = 0;
        fHistDeltaPtDeltaPhi2Rho0[i] = 0;
        fHistRCPhiEtaExLJ[i] = 0;
        fHistRhoVsRCPtExLJ[i] = 0;
        fHistRCPtExLJ[i] = 0;
        fHistDeltaPtDeltaPhi2ExLJ[i] = 0;
        fHistDeltaPtDeltaPhi2ExLJRho0[i] = 0;
        fHistJetPtRaw[i] = 0;
        fHistJetPt[i] = 0;
        fHistJetEtaPhi[i] = 0;
        fHistJetPtArea[i] = 0;
        fHistJetPtEta[i] = 0;
        fHistJetPtConstituents[i] = 0;
        fHistJetEtaRho[i] = 0;
        fHistJetPsi2Pt[i] = 0;
        fHistJetPsi2PtRho0[i] = 0;
   }
   for(Int_t i(0); i < 9; i++) {
       for(Int_t j(0); j < 2; j++) {
           for(Int_t k(0); k < 2; k++) {
               fMeanQ[i][j][k] = 0.; 
               fWidthQ[i][j][k] = 0.;  
               fMeanQv3[i][j][k] = 0.; 
               fWidthQv3[i][j][k] = 0.;
           }
       }
   }
   for(Int_t i(0); i < 4; i++) {
       fVZEROApol[i] = 0.;
       fVZEROCpol[i] = 0.;
   }
   for(Int_t i(0); i < 8; i++) fUseVZERORing[i] = kTRUE;
   // default constructor
}
//_____________________________________________________________________________
AliAnalysisTaskJetV2::AliAnalysisTaskJetV2(const char* name, runModeType type) : AliAnalysisTaskEmcalJet(name, kTRUE),
  fDebug(0), fRunToyMC(kFALSE), fLocalInit(0), fAttachToEvent(kTRUE), fFillHistograms(kTRUE), fFillQAHistograms(kTRUE), fReduceBinsXByFactor(-1.), fReduceBinsYByFactor(-1.), fNoEventWeightsForQC(kTRUE), fCentralityClasses(0), fExpectedRuns(0), fExpectedSemiGoodRuns(0), fUserSuppliedV2(0), fUserSuppliedV3(0), fUserSuppliedR2(0), fUserSuppliedR3(0), fTracksCont(0), fClusterCont(0), fJetsCont(0), fLeadingJet(0), fNAcceptedTracks(0), fNAcceptedTracksQCn(0), fFitModulationType(kNoFit), fFitGoodnessTest(kChi2Poisson), fQCRecovery(kTryFit), fUsePtWeight(kTRUE), fUsePtWeightErrorPropagation(kTRUE), fDetectorType(kVZEROComb), fAnalysisType(kCharged), fFitModulationOptions("QWLI"), fRunModeType(type), fDataType(kESD), fCollisionType(kPbPb), fRandom(0), fRunNumber(-1), fMappedRunNumber(0), fInCentralitySelection(-1), fFitModulation(0), fFitControl(0), fMinPvalue(0.01), fMaxPvalue(1), fNameSmallRho(""), fCachedRho(0), fSoftTrackMinPt(0.15), fSoftTrackMaxPt(5.), fSemiGoodJetMinPhi(0.), fSemiGoodJetMaxPhi(4.), fSemiGoodTrackMinPhi(0.), fSemiGoodTrackMaxPhi(4.), fAbsVertexZ(10), fHistCentrality(0), fHistVertexz(0), fHistRunnumbersPhi(0), fHistRunnumbersEta(0), fHistPvalueCDFROOT(0), fHistPvalueCDFROOTCent(0), fHistChi2ROOTCent(0), fHistPChi2Root(0),  fHistPvalueCDF(0), fHistPvalueCDFCent(0), fHistChi2Cent(0), fHistPChi2(0), fHistKolmogorovTest(0), fHistKolmogorovTestCent(0), fHistPKolmogorov(0), fHistRhoStatusCent(0), fHistUndeterminedRunQA(0), fMinDisanceRCtoLJ(0), fMaxCones(-1), fExcludeLeadingJetsFromFit(1.), fRebinSwapHistoOnTheFly(kTRUE), fPercentageOfFits(10.), fOutputList(0), fOutputListGood(0), fOutputListBad(0), fHistAnalysisSummary(0), fHistSwap(0), fProfV2(0), fProfV2Cumulant(0), fProfV3(0), fProfV3Cumulant(0), fHistPsiControl(0), fHistPsiSpread(0), fHistPsiVZEROA(0), fHistPsiVZEROC(0), fHistPsiVZERO(0), fHistPsiTPC(0), fHistPsiVZEROAV0M(0), fHistPsiVZEROCV0M(0), fHistPsiVZEROVV0M(0), fHistPsiTPCiV0M(0), fHistPsiVZEROATRK(0), fHistPsiVZEROCTRK(0), fHistPsiVZEROTRK(0), fHistPsiTPCTRK(0), fHistRhoVsMult(0), fHistRhoVsCent(0), fHistRhoAVsMult(0), fHistRhoAVsCent(0), fVZEROgainEqualization(0x0), fVZEROgainEqualizationPerRing(kFALSE), fChi2A(0x0), fChi2C(0x0), fChi3A(0x0), fChi3C(0x0), fOADB(0x0)
{
    for(Int_t i(0); i < 10; i++) {
        fProfV2Resolution[i] = 0;
        fProfV3Resolution[i] = 0;
        fHistPicoTrackPt[i] = 0;
        fHistPicoTrackMult[i] = 0;
        fHistPicoCat1[i] = 0;
        fHistPicoCat2[i] = 0;
        fHistPicoCat3[i] = 0;
        fHistClusterPt[i] = 0;
        fHistClusterEtaPhi[i] = 0;
        fHistClusterEtaPhiWeighted[i] = 0;
        fHistRhoPackage[i] = 0;
        fHistRho[i] = 0;
        fHistRCPhiEta[i] = 0;
        fHistRhoVsRCPt[i] = 0;
        fHistRCPt[i] = 0;
        fHistDeltaPtDeltaPhi2[i] = 0;
        fHistDeltaPtDeltaPhi2Rho0[i] = 0;
        fHistRCPhiEtaExLJ[i] = 0;
        fHistRhoVsRCPtExLJ[i] = 0;
        fHistRCPtExLJ[i] = 0;
        fHistDeltaPtDeltaPhi2ExLJ[i] = 0;
        fHistDeltaPtDeltaPhi2ExLJRho0[i] = 0;
        fHistJetPtRaw[i] = 0;
        fHistJetPt[i] = 0;
        fHistJetEtaPhi[i] = 0;
        fHistJetPtArea[i] = 0;
        fHistJetPtEta[i] = 0;
        fHistJetPtConstituents[i] = 0;
        fHistJetEtaRho[i] = 0;
        fHistJetPsi2Pt[i] = 0;
        fHistJetPsi2PtRho0[i] = 0;
   }
   for(Int_t i(0); i < 9; i++) {
       for(Int_t j(0); j < 2; j++) {
           for(Int_t k(0); k < 2; k++) {
               fMeanQ[i][j][k] = 0.; 
               fWidthQ[i][j][k] = 0.;  
               fMeanQv3[i][j][k] = 0.; 
               fWidthQv3[i][j][k] = 0.;
           }
       }
   }
   for(Int_t i(0); i < 4; i++) {
       fVZEROApol[i] = 0.;
       fVZEROCpol[i] = 0.;
   }
   for(Int_t i(0); i < 8; i++) fUseVZERORing[i] = kTRUE;

    // constructor
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());
    switch (fRunModeType) {
        case kLocal : {
            gStyle->SetOptFit(1);
            DefineOutput(2, TList::Class());
            DefineOutput(3, TList::Class());
        } break;
        default: fDebug = -1;   // suppress debug info explicitely when not running locally
    }
    switch (fCollisionType) {
        case kPythia : {
            fFitModulationType = kNoFit;
        } break;
        default : break;
    }
    if(fLocalRhoName=="") fLocalRhoName = Form("LocalRhoFrom_%s", GetName());
}
//_____________________________________________________________________________
AliAnalysisTaskJetV2::~AliAnalysisTaskJetV2()
{
    // destructor
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    if(fOutputList)             {delete fOutputList;            fOutputList = 0x0;}
    if(fOutputListGood)         {delete fOutputListGood;        fOutputListGood = 0x0;}
    if(fOutputListBad)          {delete fOutputListBad;         fOutputListBad = 0x0;}
    if(fFitModulation)          {delete fFitModulation;         fFitModulation = 0x0;}
    if(fHistSwap)               {delete fHistSwap;              fHistSwap = 0x0;}
    if(fCentralityClasses)      {delete fCentralityClasses;     fCentralityClasses = 0x0;}
    if(fExpectedRuns)           {delete fExpectedRuns;          fExpectedRuns = 0x0;}
    if(fExpectedSemiGoodRuns)   {delete fExpectedSemiGoodRuns;  fExpectedSemiGoodRuns = 0x0;}
    if(fFitControl)             {delete fFitControl;            fFitControl = 0x0;}
    if(fVZEROgainEqualization)  {delete fVZEROgainEqualization; fVZEROgainEqualization = 0x0;}
    if(fChi2A)                  {delete fChi2A;                 fChi2A = 0x0;}
    if(fChi2C)                  {delete fChi2C;                 fChi2C = 0x0;}
    if(fChi3A)                  {delete fChi3A;                 fChi3A = 0x0;}
    if(fChi3C)                  {delete fChi3C;                 fChi3C = 0x0;}
    if(fOADB && !fOADB->IsZombie()) {
        fOADB->Close();        fOADB = 0x0;
    } else if (fOADB) fOADB = 0x0;
}
//_____________________________________________________________________________
void AliAnalysisTaskJetV2::ExecOnce()
{
    // Init the analysis
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    fLocalRho = new AliLocalRhoParameter(fLocalRhoName.Data(), 0); 
    if(fAttachToEvent) {
        if(!(InputEvent()->FindListObject(fLocalRho->GetName()))) {
            InputEvent()->AddObject(fLocalRho);
        } else {
            AliFatal(Form("%s: Container with name %s already present. Aborting", GetName(), fLocalRho->GetName()));
        }
    }
    AliAnalysisTaskEmcalJet::ExecOnce();        // init the base class
    AliAnalysisTaskEmcalJet::SetVzRange(-1.*fAbsVertexZ, fAbsVertexZ);
    if(!GetJetContainer()) AliFatal(Form("%s: Couldn't find jet container. Aborting !", GetName()));
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskJetV2::Notify()
{
    // determine the run number to see if the track and jet cuts should be refreshed for semi-good TPC runs
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    if(fRunNumber != InputEvent()->GetRunNumber()) {
        fRunNumber = InputEvent()->GetRunNumber();        // set the current run number
        if(fDebug > 0) printf("__FUNC__ %s > NEW RUNNUMBER DETECTED \n ", __func__);
        // check if this is 10h or 11h data
        switch (fCollisionType) {
            case kPbPb10h : {
                if(fDebug > 0) printf(" LHC10h data, assuming full acceptance, reading VZERO calibration DB \n ");
                // for 10h data the vzero event plane calibration needs to be cached
                ReadVZEROCalibration2010h(); 
                // no need to change rho or acceptance for 10h, so we're done
                return kTRUE;
            } break;
            case kJetFlowMC : {
                return kTRUE;
            } break;
            default :  {
                if(fDebug > 0) printf(" checking runnumber to adjust acceptance on the fly \n");           
            } break;
        }
        // reset the cuts. should be a pointless operation except for the case where the run number changes
        // from semi-good back to good on one node, which is not a likely scenario (unless trains will
        // run as one masterjob)
        switch (fAnalysisType) {
            case kCharged: {
                AliAnalysisTaskEmcalJet::SetJetPhiLimits(-10., 10.);   
            } break;
            case kFull: {
                AliAnalysisTaskEmcalJet::SetJetPhiLimits(1.405 + GetJetRadius(), 3.135 - GetJetRadius());
            } break;
            default: {
                AliAnalysisTaskEmcal::SetTrackPhiLimits(-10., 10.);
            } break;
        }
        if(fCachedRho) {                // if there's a cached rho, it's the default, so switch back
            if(fDebug > 0) printf("__FUNC__ %s > replacing rho with cached rho \n ", __func__);
            fRho = fCachedRho;          // reset rho back to cached value. again, should be pointless
        }
        Bool_t flaggedAsSemiGood(kFALSE);       // not flagged as anything
        for(Int_t i(0); i < fExpectedSemiGoodRuns->GetSize(); i++) {
            if(fExpectedSemiGoodRuns->At(i) == fRunNumber) { // run is semi-good
               if(fDebug > 0) printf("__FUNC__ %s > semi-good tpc run detected, adjusting acceptance \n ", __func__);
                flaggedAsSemiGood = kTRUE;
                switch (fAnalysisType) {
                    // for full jets the jet acceptance does not have to be changed as emcal does not
                    // cover the tpc low voltage readout strips
                    case kCharged: {
                        AliAnalysisTaskEmcalJet::SetJetPhiLimits(fSemiGoodJetMinPhi, fSemiGoodJetMaxPhi);       // just an acceptance cut, jets are obtained from full azimuth, so no edge effects
                    } break;
                    default: break;
                }
                AliAnalysisTaskEmcal::SetTrackPhiLimits(fSemiGoodTrackMinPhi, fSemiGoodTrackMaxPhi);    // only affects vn extraction, NOT jet finding
                // for semi-good runs, also try to get the 'small rho' estimate, if it is available
                AliRhoParameter* tempRho(dynamic_cast<AliRhoParameter*>(InputEvent()->FindListObject(fNameSmallRho.Data())));
                if(tempRho) {
                    if(fDebug > 0) printf("__FUNC__ %s > switching to small rho, caching normal rho \n ", __func__);
                    fHistAnalysisSummary->SetBinContent(54, 1.);        // bookkeep the fact that small rho is used
                    fCachedRho = fRho;          // cache the original rho ...
                    fRho = tempRho;             // ... and use the small rho
                }
            }
        }
        if(!flaggedAsSemiGood) {
            // in case the run is not a semi-good run, check if it is recognized as another run
            // only done to catch unexpected runs
            for(Int_t i(0); i < fExpectedRuns->GetSize(); i++) {
                if(fExpectedRuns->At(i) == fRunNumber) break; // run is known, break the loop else store the number in a random bin
                fHistUndeterminedRunQA->SetBinContent(TMath::Nint(10.*gRandom->Uniform(0.,.9))+1, fRunNumber);
            }
            fHistAnalysisSummary->SetBinContent(53, 1.);                // bookkeep which rho estimate is used 
        }
    }
    return kTRUE;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskJetV2::InitializeAnalysis() 
{
    // initialize the anaysis
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    // if not set, estimate the number of cones that would fit into the selected acceptance
    if(fMaxCones <= 0) fMaxCones = TMath::CeilNint((TMath::Abs(GetJetContainer()->GetJetEtaMax()-GetJetContainer()->GetJetEtaMin())*TMath::Abs(GetJetContainer()->GetJetPhiMax()-GetJetContainer()->GetJetPhiMin()))/(TMath::Pi()*GetJetRadius()*GetJetRadius()));
    // manually 'override' the default acceptance cuts of the emcal framework (use with caution) 
    if(fMinDisanceRCtoLJ==0) fMinDisanceRCtoLJ = GetJetRadius();
    if(dynamic_cast<AliAODEvent*>(InputEvent())) fDataType = kAOD; // determine the datatype
    else if(dynamic_cast<AliESDEvent*>(InputEvent())) fDataType = kESD;
    fHistAnalysisSummary->SetBinContent(36, (int)fDataType);
    if(!fRandom) fRandom = new TRandom3(0);  // set randomizer and random seed
    switch (fFitModulationType)  {
        case kNoFit : { SetModulationFit(new TF1("fix_kNoFit", "[0]", 0, TMath::TwoPi())); } break;
        case kV2 : {
            SetModulationFit(new TF1("fit_kV2", "[0]*([1]+[2]*[3]*TMath::Cos([2]*(x-[4])))", 0, TMath::TwoPi()));
            fFitModulation->SetParameter(0, 0.);        // normalization
            fFitModulation->SetParameter(3, 0.2);       // v2
            fFitModulation->FixParameter(1, 1.);        // constant
            fFitModulation->FixParameter(2, 2.);        // constant
        } break;
        case kV3: {
            SetModulationFit(new TF1("fit_kV3", "[0]*([1]+[2]*[3]*TMath::Cos([2]*(x-[4])))", 0, TMath::TwoPi()));
            fFitModulation->SetParameter(0, 0.);        // normalization
            fFitModulation->SetParameter(3, 0.2);       // v3
            fFitModulation->FixParameter(1, 1.);        // constant
            fFitModulation->FixParameter(2, 3.);        // constant
        } break;
        default : { // for the combined fit, the 'direct fourier series' or the user supplied vn values we use v2 and v3
             SetModulationFit(new TF1("fit_kCombined", "[0]*([1]+[2]*([3]*TMath::Cos([2]*(x-[4]))+[7]*TMath::Cos([5]*(x-[6]))))", 0, TMath::TwoPi()));
             fFitModulation->SetParameter(0, 0.);       // normalization
             fFitModulation->SetParameter(3, 0.2);      // v2
             fFitModulation->FixParameter(1, 1.);       // constant
             fFitModulation->FixParameter(2, 2.);       // constant
             fFitModulation->FixParameter(5, 3.);       // constant
             fFitModulation->SetParameter(7, 0.2);      // v3
        } break;
    }
    switch (fRunModeType) {
        case kGrid : { fFitModulationOptions += "N0"; } break;
        default : break;
    }
    FillAnalysisSummaryHistogram();
    return kTRUE;
}
//_____________________________________________________________________________
TH1F* AliAnalysisTaskJetV2::BookTH1F(const char* name, const char* x, Int_t bins, Double_t min, Double_t max, Int_t c, Bool_t append)
{
    // book a TH1F and connect it to the output container
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    if(fReduceBinsXByFactor > 0 ) bins = TMath::Nint(bins/fReduceBinsXByFactor);
    if(!fOutputList) return 0x0;
    TString title(name);
    if(c!=-1) { // format centrality dependent histograms accordingly
        name = Form("%s_%i", name, c);
        title += Form("_%i-%i", (int)(fCentralityClasses->At(c)), (int)(fCentralityClasses->At((1+c))));
    }
    title += Form(";%s;[counts]", x);
    TH1F* histogram = new TH1F(name, title.Data(), bins, min, max);
    histogram->Sumw2();
    if(append) fOutputList->Add(histogram);
    return histogram;   
}
//_____________________________________________________________________________
TH2F* AliAnalysisTaskJetV2::BookTH2F(const char* name, const char* x, const char*y, Int_t binsx, Double_t minx, Double_t maxx, Int_t binsy, Double_t miny, Double_t maxy, Int_t c, Bool_t append)
{
    // book a TH2F and connect it to the output container
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    if(fReduceBinsXByFactor > 0 ) binsx = TMath::Nint(binsx/fReduceBinsXByFactor);
    if(fReduceBinsYByFactor > 0 ) binsy = TMath::Nint(binsy/fReduceBinsYByFactor);
    if(!fOutputList) return 0x0;
    TString title(name);
    if(c!=-1) { // format centrality dependent histograms accordingly
        name = Form("%s_%i", name, c);
        title += Form("_%i-%i", (int)fCentralityClasses->At(c), (int)(fCentralityClasses->At((1+c))));
    }
    title += Form(";%s;%s", x, y);
    TH2F* histogram = new TH2F(name, title.Data(), binsx, minx, maxx, binsy, miny, maxy);
    histogram->Sumw2();
    if(append) fOutputList->Add(histogram);
    return histogram;   
}
//_____________________________________________________________________________
void AliAnalysisTaskJetV2::UserCreateOutputObjects()
{
    // create output objects. also initializes some default values in case they aren't 
    // loaded via the AddTask macro
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    fOutputList = new TList();
    fOutputList->SetOwner(kTRUE);
    if(!fCentralityClasses) {   // classes must be defined at this point
        Double_t c[] = {0., 20., 40., 60., 80., 100.};
        fCentralityClasses = new TArrayD(sizeof(c)/sizeof(c[0]), c);
    }
    if(!fExpectedRuns) {        // expected runs must be defined at this point
        Int_t r[] =  {167813, 167988, 168066, 168068, 168069, 168076, 168104, 168212, 168311, 168322, 168325, 168341, 168361, 168362, 168458, 168460, 168461, 168992, 169091, 169094, 169138, 169143, 169167, 169417, 169835, 169837, 169838, 169846, 169855, 169858, 169859, 169923, 169956, 170027, 170036, 170081, /* up till here original good TPC list */169975, 169981, 170038, 170040, 170083, 170084, 170085, 170088, 170089, 170091, 170152, 170155, 170159, 170163, 170193, 170195, 170203, 170204, 170205, 170228, 170230, 170264, 170268, 170269, 170270, 170306, 170308, 170309, /* original semi-good tpc list */169415, 169411, 169035, 168988, 168984, 168826, 168777, 168512, 168511, 168467, 168464, 168342, 168310, 168115, 168108, 168107, 167987, 167915, 167903, /*new runs, good according to RCT */ 169238, 169160, 169156, 169148, 169145, 169144 /* run swith missing OROC 8 but seem ok in QA */};
        fExpectedRuns = new TArrayI(sizeof(r)/sizeof(r[0]), r);
    }
    // set default semi-good runs only for 11h data
    switch (fCollisionType) {
        case kPbPb10h : break;
        default : {
            if(!fExpectedSemiGoodRuns) {
                Int_t r[] = {169975, 169981, 170038, 170040, 170083, 170084, 170085, 170088, 170089, 170091, 170152, 170155, 170159, 170163, 170193, 170195, 170203, 170204, 170205, 170228, 170230, 170264, 170268, 170269, 170270, 170306, 170308, 170309};
                fExpectedSemiGoodRuns = new TArrayI(sizeof(r)/sizeof(r[0]), r);
            }
        }
    }

    // global QA
    fHistCentrality =           BookTH1F("fHistCentrality", "centrality", 102, -2, 100);
    fHistVertexz =              BookTH1F("fHistVertexz", "vertex z (cm)", 100, -12, 12);

    // pico track and emcal cluster kinematics
    for(Int_t i(0); i < fCentralityClasses->GetSize()-1; i++) { 
        fHistPicoTrackPt[i] =           BookTH1F("fHistPicoTrackPt", "p_{t} [GeV/c]", 100, 0, 100, i);
        fHistPicoTrackMult[i] =         BookTH1F("fHistPicoTrackMult", "multiplicity", 100, 0, 5000, i);
        if(fFillQAHistograms) {
            fHistPicoCat1[i] =          BookTH2F("fHistPicoCat1", "#eta", "#phi", 50, -1, 1, 50, 0, TMath::TwoPi(), i);
            fHistPicoCat2[i] =          BookTH2F("fHistPicoCat2", "#eta", "#phi", 50, -1, 1, 50, 0, TMath::TwoPi(), i);
            fHistPicoCat3[i] =          BookTH2F("fHistPicoCat3", "#eta", "#phi", 50, -1, 1, 50, 0, TMath::TwoPi(), i);
            if(fAnalysisType == AliAnalysisTaskJetV2::kFull) {
                fHistClusterPt[i] =     BookTH1F("fHistClusterPt", "p_{t} [GeV/c]", 100, 0, 100, i);
                fHistClusterEtaPhi[i] = BookTH2F("fHistClusterEtaPhi", "#eta", "#phi", 100, -1., 1., 100, 0, TMath::TwoPi(), i);
                fHistClusterEtaPhiWeighted[i] =    BookTH2F("fHistClusterEtaPhiWeighted", "#eta", "#phi", 100, -1., 1., 100, 0, TMath::TwoPi(), i);
            }
        }
    }

    if(fFillQAHistograms) {
        // event plane estimates and quality
        fHistPsiControl =           new TProfile("fHistPsiControl", "fHistPsiControl", 10, 0, 10);
        fHistPsiControl->Sumw2();
        fHistPsiSpread =            new TProfile("fHistPsiSpread", "fHistPsiSpread", 4, 0, 4);
        fHistPsiSpread->Sumw2();
        fHistPsiControl->GetXaxis()->SetBinLabel(1, "<#Psi_{2, VZEROA}>");
        fHistPsiControl->GetXaxis()->SetBinLabel(2, "<#Psi_{2, VZEROC}>");
        fHistPsiControl->GetXaxis()->SetBinLabel(3, "<#Psi_{2, TPC}>");
        fHistPsiControl->GetXaxis()->SetBinLabel(4, "<#Psi_{2, TPC, #eta < 0}>");
        fHistPsiControl->GetXaxis()->SetBinLabel(5, "<#Psi_{2, TPC, #eta > 0}>");
        fHistPsiControl->GetXaxis()->SetBinLabel(6, "<#Psi_{3, VZEROA}>");
        fHistPsiControl->GetXaxis()->SetBinLabel(7, "<#Psi_{3, VZEROC}>");
        fHistPsiControl->GetXaxis()->SetBinLabel(8, "<#Psi_{3, TPC}>");
        fHistPsiControl->GetXaxis()->SetBinLabel(9, "<#Psi_{3, TPC, #eta < 0}>");
        fHistPsiControl->GetXaxis()->SetBinLabel(10, "<#Psi_{3, TPC, #eta > 0}>");
        fHistPsiSpread->GetXaxis()->SetBinLabel(1, "<#Psi_{2, VZEROA} - #Psi_{2, VZEROC}>");
        fHistPsiSpread->GetXaxis()->SetBinLabel(2, "<#Psi_{2, VZEROC} - #Psi_{2, TPC}>");
        fHistPsiSpread->GetXaxis()->SetBinLabel(3, "<#Psi_{2, VZEROC} - #Psi_{2, TPC}>");
        fHistPsiSpread->GetXaxis()->SetBinLabel(4, "<#Psi_{2, TPC, #eta < 0} - #Psi_{2, TPC, #eta > 0}>");
        fOutputList->Add(fHistPsiControl);
        fOutputList->Add(fHistPsiSpread);
        fHistPsiVZEROA =            BookTH1F("fHistPsiVZEROA", "#Psi_{VZEROA}", 40, -.5*TMath::Pi(), .5*TMath::Pi());
        fHistPsiVZEROC =            BookTH1F("fHistPsiVZEROC", "#Psi_{VZEROC}", 40, -.5*TMath::Pi(), .5*TMath::Pi());
        fHistPsiVZERO =             BookTH1F("fHistPsiVZERO", "#Psi_{VZERO}", 40, -.5*TMath::Pi(), .5*TMath::Pi());
        fHistPsiTPC =               BookTH1F("fHistPsiTPC", "#Psi_{TPC}", 40, -.5*TMath::Pi(), .5*TMath::Pi());
        fHistPsiVZEROAV0M =         BookTH2F("fHistPsiVZEROAV0M", "V0M", "#Psi_{2, VZEROA}", 60, 0, 60, 40, -.5*TMath::Pi(), .5*TMath::Pi());
        fHistPsiVZEROCV0M =         BookTH2F("fHistPsiVZEROCV0M", "V0M", "#Psi_{2, VZEROC}", 60, 0, 60, 40, -.5*TMath::Pi(), .5*TMath::Pi());
        fHistPsiVZEROVV0M =         BookTH2F("fHistPsiVZEROV0M", "V0M", "#Psi_{2, VZERO}", 60, 0, 60, 40, -.5*TMath::Pi(), .5*TMath::Pi());
        fHistPsiTPCiV0M =           BookTH2F("fHistPsiTPCV0M", "V0M", "#Psi_{2, TRK}", 60, 0, 60, 40, -.5*TMath::Pi(), .5*TMath::Pi());
        fHistPsiVZEROATRK =         BookTH2F("fHistPsiVZEROATRK", "TRK", "#Psi_{2, VZEROA}", 60, 0, 60, 40, -.5*TMath::Pi(), .5*TMath::Pi());
        fHistPsiVZEROCTRK =         BookTH2F("fHistPsiVZEROCTRK", "TRK", "#Psi_{2, VZEROC}", 60, 0, 60, 40, -.5*TMath::Pi(), .5*TMath::Pi());
        fHistPsiVZEROTRK =          BookTH2F("fHistPsiVZEROTRK", "TRK", "#Psi_{2, VZERO}", 60, 0, 60, 40, -.5*TMath::Pi(), .5*TMath::Pi());
        fHistPsiTPCTRK =            BookTH2F("fHistPsiTPCTRK", "TRK", "#Psi_{2, TRK}", 60, 0, 60, 40, -.5*TMath::Pi(), .5*TMath::Pi());
    }
    // background
    for(Int_t i(0); i < fCentralityClasses->GetSize()-1; i ++) {
        fHistRhoPackage[i] =           BookTH1F("fHistRhoPackage",  "#rho [GeV/c]", 100, 0, 150, i);
        fHistRho[i] =                  BookTH1F("fHistRho", "#rho [GeV/c]", 100, 0, 150, i);
    }
    fHistRhoVsMult =            BookTH2F("fHistRhoVsMult", "multiplicity", "#rho [GeV/c]", 100, 0, 4000, 100, 0, 250);
    fHistRhoVsCent =            BookTH2F("fHistRhoVsCent", "centrality", "#rho [GeV/c]", 100, 0, 100, 100, 0, 250);
    fHistRhoAVsMult =           BookTH2F("fHistRhoAVsMult", "multiplicity", "#rho * A (jet) [GeV/c]", 100, 0, 4000, 100, 0, 50);
    fHistRhoAVsCent =           BookTH2F("fHistRhoAVsCent", "centrality", "#rho * A (jet) [GeV/c]", 100, 0, 100, 100, 0, 50);

    TString detector("");
    switch (fDetectorType) {
        case kTPC : detector+="TPC";
            break;
        case kVZEROA : detector+="VZEROA";
            break;
        case kVZEROC : detector+="VZEROC";
            break;
        case kVZEROComb : detector+="VZEROComb";
            break; 
        case kFixedEP : detector+="FixedEP";
            break;
        default: break;
    }
    // delta pt distributions
    for(Int_t i(0); i < fCentralityClasses->GetSize()-1; i ++) {
        if(fFillQAHistograms)   fHistRCPhiEta[i] = BookTH2F("fHistRCPhiEta", "#phi (RC)", "#eta (RC)", 40, 0, TMath::TwoPi(), 40, -1, 1, i);
        fHistRhoVsRCPt[i] =            BookTH2F("fHistRhoVsRCPt", "p_{t} (RC) [GeV/c]", "#rho * A (RC) [GeV/c]", 100, 0, 300, 100, 0, 350, i);
        fHistRCPt[i] =                 BookTH1F("fHistRCPt", "p_{t} (RC) [GeV/c]", 130, -20, 150, i);
        if(fFillQAHistograms)   fHistRCPhiEtaExLJ[i] = BookTH2F("fHistRCPhiEtaExLJ", "#phi (RC)", "#eta (RC)", 40, 0, TMath::TwoPi(), 40, -1, 1, i);
        fHistDeltaPtDeltaPhi2[i] =  BookTH2F("fHistDeltaPtDeltaPhi2", Form("#phi - #Psi_{2, %s}", detector.Data()), "#delta p_{t} [GeV/c]", 40, 0, TMath::Pi(), 400, -70, 130, i);
        fHistDeltaPtDeltaPhi2Rho0[i] =  BookTH2F("fHistDeltaPtDeltaPhi2Rho0", Form("#phi - #Psi_{2, %s}", detector.Data()), "#delta p_{t} [GeV/c]", 40, 0, TMath::Pi(), 400, -70, 130, i);
        fHistRhoVsRCPtExLJ[i] =        BookTH2F("fHistRhoVsRCPtExLJ", "p_{t} (RC) [GeV/c]", "#rho * A (RC) [GeV/c]", 100, 0, 300, 100, 0, 350, i);
        fHistRCPtExLJ[i] =             BookTH1F("fHistRCPtExLJ", "p_{t} (RC) [GeV/c]", 130, -20, 150, i);
        fHistDeltaPtDeltaPhi2ExLJ[i] = BookTH2F("fHistDeltaPtDeltaPhi2ExLJ", Form("#phi - #Psi_{2, %s}", detector.Data()),  "#delta p_{t} [GeV/c]", 40, 0, TMath::Pi(), 400, -70, 130, i);
        fHistDeltaPtDeltaPhi2ExLJRho0[i] = BookTH2F("fHistDeltaPtDeltaPhi2ExLJRho0", Form("#phi - #Psi_{2, %s}", detector.Data()),  "#delta p_{t} [GeV/c]", 40, 0, TMath::Pi(), 400, -70, 130, i);
        // jet histograms (after kinematic cuts)
        fHistJetPtRaw[i] =             BookTH1F("fHistJetPtRaw", "p_{t, jet} RAW [GeV/c]", 200, -50, 150, i);
        fHistJetPt[i] =                BookTH1F("fHistJetPt", "p_{t, jet} [GeV/c]", 350, -100, 250, i);
        if(fFillQAHistograms)   fHistJetEtaPhi[i] =            BookTH2F("fHistJetEtaPhi", "#eta", "#phi", 100, -1, 1, 100, 0, TMath::TwoPi(), i);
        fHistJetPtArea[i] =            BookTH2F("fHistJetPtArea", "p_{t, jet} [GeV/c]", "Area", 175, -100, 250, 30, 0, 0.9, i);
        fHistJetPtEta[i] =             BookTH2F("fHistJetPtEta", "p_{t, jet} [GeV/c]", "Eta", 175, -100, 250, 30, -0.9, 0.9, i);
        fHistJetPtConstituents[i] =    BookTH2F("fHistJetPtConstituents", "p_{t, jet} [GeV/c]", "Area", 350, -100, 250, 60, 0, 150, i);
        fHistJetEtaRho[i] =            BookTH2F("fHistJetEtaRho", "#eta", "#rho", 100, -1, 1, 100, 0, 300, i);
        // in plane and out of plane spectra
        fHistJetPsi2Pt[i] =            BookTH2F("fHistJetPsi2Pt", Form("#phi_{jet} - #Psi_{2, %s}", detector.Data()), "p_{t, jet} [GeV/c]", 40, 0., TMath::Pi(), 350, -100, 250, i);
        fHistJetPsi2PtRho0[i] =        BookTH2F("fHistJetPsi2PtRho0", Form("#phi_{jet} - #Psi_{2, %s}", detector.Data()), "p_{t, jet} [GeV/c]", 40, 0., TMath::Pi(), 350, -100, 250, i);
        // profiles for all correlator permutations which are necessary to calculate each second and third order event plane resolution
        fProfV2Resolution[i] = new TProfile(Form("fProfV2Resolution_%i", i), Form("fProfV2Resolution_%i", i), 11, -0.5, 10.5);
        fProfV2Resolution[i]->GetXaxis()->SetBinLabel(3, "<cos(2(#Psi_{VZEROA} - #Psi_{VZEROC}))>");
        fProfV2Resolution[i]->GetXaxis()->SetBinLabel(4, "<cos(2(#Psi_{VZEROC} - #Psi_{VZEROA}))>");
        fProfV2Resolution[i]->GetXaxis()->SetBinLabel(5, "<cos(2(#Psi_{VZEROA} - #Psi_{TPC}))>");
        fProfV2Resolution[i]->GetXaxis()->SetBinLabel(6, "<cos(2(#Psi_{TPC} - #Psi_{VZEROA}))>");
        fProfV2Resolution[i]->GetXaxis()->SetBinLabel(7, "<cos(2(#Psi_{VZEROC} - #Psi_{TPC}))>");
        fProfV2Resolution[i]->GetXaxis()->SetBinLabel(8, "<cos(2(#Psi_{TPC} - #Psi_{VZEROC}))>");
        fProfV2Resolution[i]->GetXaxis()->SetBinLabel(9, "<cos(2(#Psi_{VZERO} - #Psi_{TPC_A}))>");
        fProfV2Resolution[i]->GetXaxis()->SetBinLabel(10, "<cos(2(#Psi_{VZERO} - #Psi_{TPC_B}))>");
        fProfV2Resolution[i]->GetXaxis()->SetBinLabel(11, "<cos(2(#Psi_{TPC_A} - #Psi_{TPC_B}))>");
        fOutputList->Add(fProfV2Resolution[i]); 
        fProfV3Resolution[i] = new TProfile(Form("fProfV3Resolution_%i", i), Form("fProfV3Resolution_%i", i), 11, -0.5, 10.5);
        fProfV3Resolution[i]->GetXaxis()->SetBinLabel(3, "<cos(3(#Psi_{VZEROA} - #Psi_{VZEROC}))>");
        fProfV3Resolution[i]->GetXaxis()->SetBinLabel(4, "<cos(3(#Psi_{VZEROC} - #Psi_{VZEROA}))>");
        fProfV3Resolution[i]->GetXaxis()->SetBinLabel(5, "<cos(3(#Psi_{VZEROA} - #Psi_{TPC}))>");
        fProfV3Resolution[i]->GetXaxis()->SetBinLabel(6, "<cos(3(#Psi_{TPC} - #Psi_{VZEROA}))>");
        fProfV3Resolution[i]->GetXaxis()->SetBinLabel(7, "<cos(3(#Psi_{VZEROC} - #Psi_{TPC}))>");
        fProfV3Resolution[i]->GetXaxis()->SetBinLabel(8, "<cos(3(#Psi_{TPC} - #Psi_{VZEROC}))>");
        fProfV3Resolution[i]->GetXaxis()->SetBinLabel(9, "<cos(3(#Psi_{VZERO} - #Psi_{TPC_A}))>");
        fProfV3Resolution[i]->GetXaxis()->SetBinLabel(10, "<cos(3(#Psi_{VZERO} - #Psi_{TPC_B}))>");
        fProfV3Resolution[i]->GetXaxis()->SetBinLabel(11, "<cos(3(#Psi_{TPC_A} - #Psi_{TPC_B}))>");
        fOutputList->Add(fProfV3Resolution[i]); 
    }
   // vn profile
    Float_t temp[fCentralityClasses->GetSize()];
    for(Int_t i(0); i < fCentralityClasses->GetSize(); i++) temp[i] = fCentralityClasses->At(i);
    fProfV2 = new TProfile("fProfV2", "fProfV2", fCentralityClasses->GetSize()-1, temp);
    fProfV3 = new TProfile("fProfV3", "fProfV3", fCentralityClasses->GetSize()-1, temp);
    fOutputList->Add(fProfV2);
    fOutputList->Add(fProfV3);
    switch (fFitModulationType) {
        case kQC2 : {
            fProfV2Cumulant = new TProfile("fProfV2Cumulant", "fProfV2Cumulant", fCentralityClasses->GetSize()-1, temp);
            fProfV3Cumulant = new TProfile("fProfV3Cumulant", "fProfV3Cumulant", fCentralityClasses->GetSize()-1, temp);
            fOutputList->Add(fProfV2Cumulant);
            fOutputList->Add(fProfV3Cumulant);
        } break;
        case kQC4 : {
            fProfV2Cumulant = new TProfile("fProfV2Cumulant", "fProfV2Cumulant", fCentralityClasses->GetSize()-1, temp);
            fProfV3Cumulant = new TProfile("fProfV3Cumulant", "fProfV3Cumulant", fCentralityClasses->GetSize()-1, temp);
            fOutputList->Add(fProfV2Cumulant);
            fOutputList->Add(fProfV3Cumulant);
        } break;
        default : break;
    }
    // for the histograms initialized below, binning is fixed to runnumbers or flags
    fReduceBinsXByFactor = 1;
    fReduceBinsYByFactor = 1;
    if(fFillQAHistograms) {
        fHistRunnumbersEta = new TH2F("fHistRunnumbersEta", "fHistRunnumbersEta", fExpectedRuns->GetSize()+1, -.5, fExpectedRuns->GetSize()+.5, 100, -1.1, 1.1);
        fHistRunnumbersEta->Sumw2();
        fOutputList->Add(fHistRunnumbersEta);
        fHistRunnumbersPhi = new TH2F("fHistRunnumbersPhi", "fHistRunnumbersPhi", fExpectedRuns->GetSize()+1, -.5, fExpectedRuns->GetSize()+.5, 100, -0.2, TMath::TwoPi()+0.2);
        fHistRunnumbersPhi->Sumw2();
        fOutputList->Add(fHistRunnumbersPhi);
        for(Int_t i(0); i < fExpectedRuns->GetSize(); i++) { 
            fHistRunnumbersPhi->GetXaxis()->SetBinLabel(i+1, Form("%i", fExpectedRuns->At(i)));
            fHistRunnumbersEta->GetXaxis()->SetBinLabel(i+1, Form("%i", fExpectedRuns->At(i)));
        }
        fHistRunnumbersPhi->GetXaxis()->SetBinLabel(fExpectedRuns->GetSize()+1, "undetermined");
        fHistRunnumbersEta->GetXaxis()->SetBinLabel(fExpectedRuns->GetSize()+1, "undetermined");
    }
    fHistAnalysisSummary = BookTH1F("fHistAnalysisSummary", "flag", 54, -0.5, 54.5);
    fHistSwap = new TH1F("fHistSwap", "fHistSwap", 20, 0, TMath::TwoPi());
    if(fUsePtWeight) fHistSwap->Sumw2();

    if(fUserSuppliedV2) fOutputList->Add(fUserSuppliedV2);
    if(fUserSuppliedV3) fOutputList->Add(fUserSuppliedV3);
    if(fUserSuppliedR2) fOutputList->Add(fUserSuppliedR2);
    if(fUserSuppliedR3) fOutputList->Add(fUserSuppliedR3);
    // increase readability of output list
    fOutputList->Sort();
    // cdf and pdf of chisquare distribution
    fHistPvalueCDF = BookTH1F("fHistPvalueCDF", "CDF #chi^{2}", 50, 0, 1);
    fHistPvalueCDFCent = BookTH2F("fHistPvalueCDFCent", "centrality", "p-value", 40, 0, 100, 40, 0, 1);
    fHistChi2Cent = BookTH2F("fHistChi2Cent", "centrality", "#tilde{#chi^{2}}", 100, 0, 100, 100, 0, 5);
    fHistPChi2 = BookTH2F("fHistPChi2", "p-value", "#tilde{#chi^{2}}", 1000, 0, 1, 100, 0, 5);
    fHistKolmogorovTest = BookTH1F("fHistKolmogorovTest", "KolmogorovTest", 50, 0, 1);
    fHistKolmogorovTestCent = BookTH2F("fHistKolmogorovTestCent", "centrality", "Kolmogorov p", 40, 0, 100, 45, 0, 1); 
    fHistPvalueCDFROOT = BookTH1F("fHistPvalueCDFROOT", "CDF #chi^{2} ROOT", 50, 0, 1);
    fHistPvalueCDFROOTCent = BookTH2F("fHistPvalueCDFROOTCent", "centrality", "p-value ROOT", 40, 0, 100, 45, 0, 1);
    fHistChi2ROOTCent = BookTH2F("fHistChi2ROOTCent", "centrality", "#tilde{#chi^{2}}", 40, 0, 100, 45, 0, 5);
    fHistPChi2Root = BookTH2F("fHistPChi2Root", "p-value", "#tilde{#chi^{2}} ROOT", 1000, 0, 1, 100, 0, 5);
    fHistPKolmogorov = BookTH2F("fHistPKolmogorov", "p-value", "kolmogorov p",40, 0, 1, 40, 0, 1);
    fHistRhoStatusCent = BookTH2F("fHistRhoStatusCent", "centrality", "status [-1=lin was better, 0=ok, 1 = failed]", 101, -1, 100, 3, -1.5, 1.5);
    fHistUndeterminedRunQA = BookTH1F("fHistUndeterminedRunQA", "runnumber", 10, 0, 10);
 
    PostData(1, fOutputList);

    switch (fRunModeType) {
        case kLocal : {
            fOutputListGood = new TList();
            fOutputListGood->SetOwner(kTRUE);
            fOutputListBad = new TList();
            fOutputListBad->SetOwner(kTRUE);
            PostData(2, fOutputListGood);
            PostData(3, fOutputListBad);
        } break;
        default: break;
    }

    // get the containers
    fTracksCont = GetParticleContainer("Tracks");
    fClusterCont = GetClusterContainer(0);      // get the default cluster container
    fJetsCont = GetJetContainer("Jets");
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskJetV2::Run()
{
    // called for each accepted event (call made from user exec of parent class)
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    if(!fTracks||!fJets||!fRho) {
        printf(" > Failed to retrieve input objects ! < \n");
        return kFALSE;
    }
    if(!fLocalInit) fLocalInit = InitializeAnalysis();
    // reject the event if expected data is missing
    if(!PassesCuts(InputEvent())) return kFALSE;
    fLeadingJet = GetLeadingJet();      // store the leading jet
    // set the rho value 
    fLocalRho->SetVal(fRho->GetVal());
    // place holder arrays for the event planes
    //
    // [0][0] psi2a     [1,0]   psi2c
    // [0][1] psi3a     [1,1]   psi3c
    Double_t vzero[2][2];
    /* for the combined vzero event plane
     * [0] psi2         [1] psi3
     * not fully implmemented yet, use with caution ! */
    Double_t vzeroComb[2];
    // [0] psi2         [1] psi3
    Double_t tpc[2];
    // evaluate the actual event planes
    switch (fDetectorType) {
        case kFixedEP : {
            // for fixed, fix all ep's to default values
            tpc[0] = 0.;         tpc[1] = 1.;
            vzero[0][0] = 0.;    vzero[0][1] = 1.;
            vzero[1][0] = 0.;    vzero[1][1] = 1.;
            vzeroComb[0] = 0.;   vzeroComb[1] = 1.;
        } break;
        default : {
            // else grab the actual data
            CalculateEventPlaneVZERO(vzero);
            CalculateEventPlaneCombinedVZERO(vzeroComb);
            CalculateEventPlaneTPC(tpc);
        } break;
    }
    Double_t psi2(-1), psi3(-1);
    // arrays which will hold the fit parameters
    switch (fDetectorType) {    // determine the detector type for the rho fit
        case kTPC :     { psi2 = tpc[0];         psi3 = tpc[1]; }       break;
        case kVZEROA :  { psi2 = vzero[0][0];    psi3 = vzero[0][1]; }  break;  
        case kVZEROC :  { psi2 = vzero[1][0];    psi3 = vzero[1][1]; }  break;
        case kVZEROComb : { psi2 = vzeroComb[0]; psi3 = vzeroComb[1];}  break;
        case kFixedEP : { psi2 = 0.;             psi3 = 1.;}            break;
        default : break;
    }
    switch (fFitModulationType) { // do the fits
        case kNoFit : { 
             switch (fCollisionType) {
                 case kPythia : { // background is zero for pp jets
                     fFitModulation->FixParameter(0, 0);
                     fLocalRho->SetVal(0);
                 } break;
                 default :  {
                     fFitModulation->FixParameter(0, fLocalRho->GetVal()); 
                 } break;
             }
        } break;
        case kV2 : {    // only v2
            if(CorrectRho(psi2, psi3)) {
                fProfV2->Fill(fCent, fFitModulation->GetParameter(3));
                if(fUserSuppliedR2) {
                    Double_t r(fUserSuppliedR2->GetBinContent(fUserSuppliedR2->GetXaxis()->FindBin(fCent)));
                    if(r > 0) fFitModulation->SetParameter(3, fFitModulation->GetParameter(3)/r);
                }
                CalculateEventPlaneResolution(vzero, vzeroComb, tpc);
            }
        } break;
        case kV3 : {    // only v3
            if(CorrectRho(psi2, psi3)) {
                if(fUserSuppliedR3) {
                    Double_t r(fUserSuppliedR3->GetBinContent(fUserSuppliedR3->GetXaxis()->FindBin(fCent)));
                    if(r > 0) fFitModulation->SetParameter(3, fFitModulation->GetParameter(3)/r);
                }
                fProfV3->Fill(fCent, fFitModulation->GetParameter(3));
                CalculateEventPlaneResolution(vzero, vzeroComb, tpc);
            }
        } break;
        case kQC2 : {   // qc2 analysis
            if(CorrectRho(psi2, psi3)) {
                if(fUserSuppliedR2 && fUserSuppliedR3) {
                    // note for the qc method, resolution is REVERSED to go back to v2obs
                    Double_t r2(fUserSuppliedR2->GetBinContent(fUserSuppliedR2->GetXaxis()->FindBin(fCent)));
                    Double_t r3(fUserSuppliedR3->GetBinContent(fUserSuppliedR3->GetXaxis()->FindBin(fCent)));
                    if(r2 > 0) fFitModulation->SetParameter(3, fFitModulation->GetParameter(3)*r2);
                    if(r3 > 0) fFitModulation->SetParameter(7, fFitModulation->GetParameter(7)*r3);
                }
                if (fUsePtWeight) { // use weighted weights
                    Double_t dQCnM11 = (fNoEventWeightsForQC) ? 1. : QCnM11();
                    fProfV2->Fill(fCent, fFitModulation->GetParameter(3), dQCnM11);
                    fProfV3->Fill(fCent, fFitModulation->GetParameter(7), dQCnM11); 
                } else {
                    Double_t dQCnM = (fNoEventWeightsForQC) ? 2. : QCnM();
                    fProfV2->Fill(fCent, fFitModulation->GetParameter(3), dQCnM*(dQCnM-1));
                    fProfV3->Fill(fCent, fFitModulation->GetParameter(7), dQCnM*(dQCnM-1));
                }
                CalculateEventPlaneResolution(vzero, vzeroComb, tpc);
            }
        } break;
        case kQC4 : {
            if(CorrectRho(psi2, psi3)) {
                if(fUserSuppliedR2 && fUserSuppliedR3) {
                    // note for the qc method, resolution is REVERSED to go back to v2obs   
                    Double_t r2(fUserSuppliedR2->GetBinContent(fUserSuppliedR2->GetXaxis()->FindBin(fCent)));
                    Double_t r3(fUserSuppliedR3->GetBinContent(fUserSuppliedR3->GetXaxis()->FindBin(fCent)));
                    if(r2 > 0) fFitModulation->SetParameter(3, fFitModulation->GetParameter(3)*r2);
                    if(r3 > 0) fFitModulation->SetParameter(7, fFitModulation->GetParameter(7)*r3);
                }
                if (fUsePtWeight) { // use weighted weights
                    fProfV2->Fill(fCent, TMath::Power(fFitModulation->GetParameter(3),0.5)/*, QCnM1111()*/);
                    fProfV3->Fill(fCent, TMath::Power(fFitModulation->GetParameter(7),0.5)/*, QCnM1111()*/); 
                } else {
                    fProfV2->Fill(fCent, TMath::Power(fFitModulation->GetParameter(3),0.5)/*, QCnM()*(QCnM()-1)*(QCnM()-2)*(QCnM()-3)*/);
                    fProfV3->Fill(fCent, TMath::Power(fFitModulation->GetParameter(7),0.5)/*, QCnM()*(QCnM()-1)*(QCnM()-2)*(QCnM()-3)*/);
                }
            }
            CalculateEventPlaneResolution(vzero, vzeroComb, tpc);
        } break;
        default : {
            if(CorrectRho(psi2, psi3)) {
                if(fUserSuppliedR2 && fUserSuppliedR3) {
                    Double_t r2(fUserSuppliedR2->GetBinContent(fUserSuppliedR2->GetXaxis()->FindBin(fCent)));
                    Double_t r3(fUserSuppliedR3->GetBinContent(fUserSuppliedR3->GetXaxis()->FindBin(fCent)));
                    if(r2 > 0) fFitModulation->SetParameter(3, fFitModulation->GetParameter(3)/r2);
                    if(r3 > 0) fFitModulation->SetParameter(7, fFitModulation->GetParameter(7)/r3);
                }
                fProfV2->Fill(fCent, fFitModulation->GetParameter(3));
                fProfV3->Fill(fCent, fFitModulation->GetParameter(7));
                CalculateEventPlaneResolution(vzero, vzeroComb, tpc);
            }
        } break;
    }
    // if all went well, update the local rho parameter
    fLocalRho->SetLocalRho(fFitModulation);
    // fill a number of histograms. event qa needs to be filled first as it also determines the runnumber for the track qa 
    if(fFillQAHistograms)       FillQAHistograms(InputEvent());
    if(fFillHistograms)         FillHistogramsAfterSubtraction(psi2, vzero, vzeroComb, tpc);
    // send the output to the connected output container
    PostData(1, fOutputList);
    switch (fRunModeType) {
        case kLocal : {
            PostData(2, fOutputListGood);
            PostData(3, fOutputListBad);
        } break;
        default: break;
    }

    return kTRUE;
}
//_____________________________________________________________________________
void AliAnalysisTaskJetV2::Exec(Option_t* c)
{
    // for stand alone, avoid framework event setup
    switch (fCollisionType) {
        case kJetFlowMC : {
            // need to call ExecOnce as it is not loaded otherwise
            if(!fLocalRho) AliAnalysisTaskJetV2::ExecOnce();
            AliAnalysisTaskJetV2::Run();
        } break;
        default : {
            AliAnalysisTaskSE::Exec(c);
        } break;
    }
}  
//_____________________________________________________________________________
Double_t AliAnalysisTaskJetV2::CalculateEventPlaneChi(Double_t res)
{
    // return chi for given resolution to combine event plane estimates from two subevents
    // see Phys. Rev. C no. CS6346 (http://arxiv.org/abs/nucl-ex/9805001)
    Double_t chi(2.), delta(1.), con((TMath::Sqrt(TMath::Pi()))/(2.*TMath::Sqrt(2)));
    for (Int_t i(0); i < 15; i++) {
        chi = ((con*chi*TMath::Exp(-chi*chi/4.)*(TMath::BesselI0(chi*chi/4.)+TMath::BesselI1(chi*chi/4.))) < res) ? chi + delta : chi - delta;
        delta = delta / 2.;
    }
    return chi;
}
//_____________________________________________________________________________
void AliAnalysisTaskJetV2::CalculateEventPlaneVZERO(Double_t vzero[2][2]) const 
{
    // get the vzero event plane (a and c separately)
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    switch (fCollisionType) {
        case kPbPb10h : {
            // for 10h data, get the calibrated q-vector from the database
            Double_t QA2[] = {-999., -999.};
            Double_t QA3[] = {-999., -999.};
            Double_t QC2[] = {-999., -999.};
            Double_t QC3[] = {-999., -999.};
            CalculateQvectorVZERO(QA2, QA3, QC2, QC3);
            vzero[0][0] = .5*TMath::ATan2(QA2[1], QA2[0]);
            vzero[1][0] = .5*TMath::ATan2(QC2[1], QC2[0]);
            vzero[0][1] = (1./3.)*TMath::ATan2(QA3[1], QA3[0]);
            vzero[1][1] = (1./3.)*TMath::ATan2(QC3[1], QC3[0]);
        } break;
        default: {
            // by default use the ep from the event header (make sure EP selection task is enabeled!)
            Double_t a(0), b(0), c(0), d(0), e(0), f(0), g(0), h(0);
            vzero[0][0] = InputEvent()->GetEventplane()->CalculateVZEROEventPlane(InputEvent(), 8, 2, a, b);
            vzero[1][0] = InputEvent()->GetEventplane()->CalculateVZEROEventPlane(InputEvent(), 9, 2, c, d);
            vzero[0][1] = InputEvent()->GetEventplane()->CalculateVZEROEventPlane(InputEvent(), 8, 3, e, f);
            vzero[1][1] = InputEvent()->GetEventplane()->CalculateVZEROEventPlane(InputEvent(), 9, 3, g, h);
            return;
        }
    }
}
//_____________________________________________________________________________
void AliAnalysisTaskJetV2::CalculateEventPlaneCombinedVZERO(Double_t* comb) const
{
    // return the combined vzero event plane
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    switch (fCollisionType) {
        // for 10h data call calibration info
        case kPbPb10h : {
            // get the calibrated q-vectors
            Double_t Q2[] = {-999., -999.};            
            Double_t Q3[] = {-999., -999.};
            // return if something isn't ok from the calibration side
            CalculateQvectorCombinedVZERO(Q2, Q3);
            comb[0] = .5*TMath::ATan2(Q2[1], Q2[0]);
            comb[1] = (1./3.)*TMath::ATan2(Q3[1], Q3[0]);
        } break;
        default : {
            // for all other types use calibrated event plane from the event header
            Double_t a(0), b(0), c(0), d(0);
            comb[0] = InputEvent()->GetEventplane()->CalculateVZEROEventPlane(InputEvent(), 10, 2, a, b);
            comb[1] = InputEvent()->GetEventplane()->CalculateVZEROEventPlane(InputEvent(), 10, 3, c, d);
        } break;
    }
}
//_____________________________________________________________________________
void AliAnalysisTaskJetV2::CalculateEventPlaneTPC(Double_t* tpc)
{
   // grab the TPC event plane
   if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
   fNAcceptedTracks = 0;                // reset the track counter
   Double_t qx2(0), qy2(0);     // for psi2
   Double_t qx3(0), qy3(0);     // for psi3
   if(fTracksCont) {
       Float_t excludeInEta = -999;
       if(fExcludeLeadingJetsFromFit > 0 ) {    // remove the leading jet from ep estimate
           if(fLeadingJet) excludeInEta = fLeadingJet->Eta();
       }
       for(Int_t iTPC(0); iTPC < fTracksCont->GetNEntries(); iTPC++) {
           AliVParticle* track = fTracksCont->GetParticle(iTPC);
           if(!PassesCuts(track) || track->Pt() < fSoftTrackMinPt || track->Pt() > fSoftTrackMaxPt) continue;
           if(fExcludeLeadingJetsFromFit > 0 &&( (TMath::Abs(track->Eta() - excludeInEta) < GetJetContainer()->GetJetRadius()*fExcludeLeadingJetsFromFit ) || (TMath::Abs(track->Eta()) - GetJetContainer()->GetJetRadius() - GetJetContainer()->GetJetEtaMax() ) > 0 )) continue;
           fNAcceptedTracks++;
           qx2+= TMath::Cos(2.*track->Phi());
           qy2+= TMath::Sin(2.*track->Phi());
           qx3+= TMath::Cos(3.*track->Phi());
           qy3+= TMath::Sin(3.*track->Phi());
       }
   }
   tpc[0] = .5*TMath::ATan2(qy2, qx2);
   tpc[1] = (1./3.)*TMath::ATan2(qy3, qx3);
}
//_____________________________________________________________________________
void AliAnalysisTaskJetV2::CalculateEventPlaneResolution(Double_t vzero[2][2], Double_t* vzeroComb, Double_t* tpc)
{
    // fill the profiles for the resolution parameters
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    fProfV2Resolution[fInCentralitySelection]->Fill(2., TMath::Cos(2.*(vzero[0][0] - vzero[1][0])));
    fProfV2Resolution[fInCentralitySelection]->Fill(3., TMath::Cos(2.*(vzero[1][0] - vzero[0][0])));
    fProfV2Resolution[fInCentralitySelection]->Fill(4., TMath::Cos(2.*(vzero[0][0] - tpc[0])));
    fProfV2Resolution[fInCentralitySelection]->Fill(5., TMath::Cos(2.*(tpc[0] - vzero[0][0])));
    fProfV2Resolution[fInCentralitySelection]->Fill(6., TMath::Cos(2.*(vzero[1][0] - tpc[0])));
    fProfV2Resolution[fInCentralitySelection]->Fill(7., TMath::Cos(2.*(tpc[0] - vzero[1][0])));
    fProfV3Resolution[fInCentralitySelection]->Fill(2., TMath::Cos(3.*(vzero[0][0] - vzero[1][0])));
    fProfV3Resolution[fInCentralitySelection]->Fill(3., TMath::Cos(3.*(vzero[1][0] - vzero[0][0])));
    fProfV3Resolution[fInCentralitySelection]->Fill(4., TMath::Cos(3.*(vzero[0][0] - tpc[0])));
    fProfV3Resolution[fInCentralitySelection]->Fill(5., TMath::Cos(3.*(tpc[0] - vzero[0][0])));
    fProfV3Resolution[fInCentralitySelection]->Fill(6., TMath::Cos(3.*(vzero[1][0] - tpc[0])));
    fProfV3Resolution[fInCentralitySelection]->Fill(7., TMath::Cos(3.*(tpc[0] - vzero[1][0])));
    // for the resolution of the combined vzero event plane, use two tpc halves as uncorrelated subdetectors
    Double_t qx2a(0), qy2a(0);     // for psi2a, negative eta
    Double_t qx3a(0), qy3a(0);     // for psi3a, negative eta
    Double_t qx2b(0), qy2b(0);     // for psi2a, positive eta
    Double_t qx3b(0), qy3b(0);     // for psi3a, positive eta
    if(fTracks) {
       Int_t iTracks(fTracks->GetEntriesFast());
       for(Int_t iTPC(0); iTPC < iTracks; iTPC++) {
           AliVTrack* track = static_cast<AliVTrack*>(fTracks->At(iTPC));
           if(!PassesCuts(track) || track->Pt() < fSoftTrackMinPt || track->Pt() > fSoftTrackMaxPt) continue;
           if(track->Eta() < 0 ) {
               qx2a+= TMath::Cos(2.*track->Phi());
               qy2a+= TMath::Sin(2.*track->Phi());
               qx3a+= TMath::Cos(3.*track->Phi());
               qy3a+= TMath::Sin(3.*track->Phi());
           } else if (track->Eta() > 0) {
               qx2b+= TMath::Cos(2.*track->Phi());
               qy2b+= TMath::Sin(2.*track->Phi());
               qx3b+= TMath::Cos(3.*track->Phi());
               qy3b+= TMath::Sin(3.*track->Phi());
           }
       }
   }
   Double_t tpca2(.5*TMath::ATan2(qy2a, qx2a));
   Double_t tpca3((1./3.)*TMath::ATan2(qy3a, qx3a));
   Double_t tpcb2(.5*TMath::ATan2(qy2b, qx2b));
   Double_t tpcb3((1./3.)*TMath::ATan2(qy3b, qx3b));
   fProfV2Resolution[fInCentralitySelection]->Fill(8., TMath::Cos(2.*(vzeroComb[0] - tpca2)));
   fProfV2Resolution[fInCentralitySelection]->Fill(9., TMath::Cos(2.*(vzeroComb[0] - tpcb2)));
   fProfV2Resolution[fInCentralitySelection]->Fill(10., TMath::Cos(2.*(tpca2 - tpcb2))); 
   fProfV3Resolution[fInCentralitySelection]->Fill(8., TMath::Cos(3.*(vzeroComb[1] - tpca3)));
   fProfV3Resolution[fInCentralitySelection]->Fill(9., TMath::Cos(3.*(vzeroComb[1] - tpcb3)));
   fProfV3Resolution[fInCentralitySelection]->Fill(10., TMath::Cos(3.*(tpca3 - tpcb3))); 
}   
//_____________________________________________________________________________
void AliAnalysisTaskJetV2::CalculateQvectorVZERO(Double_t Qa2[2], Double_t Qc2[2], Double_t Qa3[2], Double_t Qc3[2]) const
{
    // return the calibrated 2nd and 3rd order q-vectors for vzeroa and vzeroc
    // function takes arrays as arguments, which correspond to vzero info in the following way
    // 
    // Qa2[0] = Qx2 for vzero A         Qa2[1] = Qy2 for vzero A (etc)
    
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    // placeholders for geometric information
    Double_t phi(-999.), weight(-999.); 
    // reset placeholders for Q-vector components
    Qa2[0] = 0.;    Qc2[0] = 0.;    Qa3[0] = 0.;    Qc3[0] = 0.;
    Qa2[1] = 0.;    Qc2[1] = 0.;    Qa3[1] = 0.;    Qc3[1] = 0.;
    
    for(Int_t i(0); i < 64; i++) {
        // loop over all scintillators, construct Q-vectors in the same loop
        phi     = TMath::PiOver4()*(0.5+i%8);
        weight  = 0.;
        // note that disabled rings have already been excluded in ReadVZEROCalibration2010h
        if(i<32) {    // v0c side
            if(i < 8) weight = InputEvent()->GetVZEROData()->GetMultiplicity(i)*fVZEROCpol[0]/fVZEROgainEqualization->GetBinContent(1+i);
            else if (i < 16 ) weight = InputEvent()->GetVZEROData()->GetMultiplicity(i)*fVZEROCpol[1]/fVZEROgainEqualization->GetBinContent(1+i);
            else if (i < 24 ) weight = InputEvent()->GetVZEROData()->GetMultiplicity(i)*fVZEROCpol[2]/fVZEROgainEqualization->GetBinContent(1+i);
            else if (i < 32 ) weight = InputEvent()->GetVZEROData()->GetMultiplicity(i)*fVZEROCpol[3]/fVZEROgainEqualization->GetBinContent(1+i);
            // fill Q-vectors for v0c side
            Qc2[0]+=weight*TMath::Cos(2.*phi);
            Qc3[0]+=weight*TMath::Cos(3.*phi);
            Qc2[1]+=weight*TMath::Sin(2.*phi);
            Qc3[1]+=weight*TMath::Sin(3.*phi);
        } else {       // v0a side
            if( i < 40) weight = InputEvent()->GetVZEROData()->GetMultiplicity(i)*fVZEROApol[0]/fVZEROgainEqualization->GetBinContent(1+i);
            else if ( i < 48 ) weight = InputEvent()->GetVZEROData()->GetMultiplicity(i)*fVZEROApol[1]/fVZEROgainEqualization->GetBinContent(1+i);
            else if ( i < 56 ) weight = InputEvent()->GetVZEROData()->GetMultiplicity(i)*fVZEROApol[2]/fVZEROgainEqualization->GetBinContent(1+i);
            else if ( i < 64 ) weight = InputEvent()->GetVZEROData()->GetMultiplicity(i)*fVZEROApol[3]/fVZEROgainEqualization->GetBinContent(1+i);
            // fill Q-vectors for v0a side
            Qa2[0]+=weight*TMath::Cos(2.*phi);
            Qa3[0]+=weight*TMath::Cos(3.*phi);
            Qa2[1]+=weight*TMath::Sin(2.*phi);
            Qa3[1]+=weight*TMath::Sin(3.*phi);
        }
    }
    // get the cache index and read the correction terms from the cache
    Int_t VZEROcentralityBin(GetVZEROCentralityBin());
    Double_t Qx2amean = fMeanQ[VZEROcentralityBin][1][0];
    Double_t Qx2arms  = fWidthQ[VZEROcentralityBin][1][0];
    Double_t Qy2amean = fMeanQ[VZEROcentralityBin][1][1];
    Double_t Qy2arms  = fWidthQ[VZEROcentralityBin][1][1];

    Double_t Qx2cmean = fMeanQ[VZEROcentralityBin][0][0];
    Double_t Qx2crms  = fWidthQ[VZEROcentralityBin][0][0];
    Double_t Qy2cmean = fMeanQ[VZEROcentralityBin][0][1];
    Double_t Qy2crms  = fWidthQ[VZEROcentralityBin][0][1];	

    Double_t Qx3amean = fMeanQv3[VZEROcentralityBin][1][0];
    Double_t Qx3arms  = fWidthQv3[VZEROcentralityBin][1][0];
    Double_t Qy3amean = fMeanQv3[VZEROcentralityBin][1][1];
    Double_t Qy3arms  = fWidthQv3[VZEROcentralityBin][1][1];

    Double_t Qx3cmean = fMeanQv3[VZEROcentralityBin][0][0];
    Double_t Qx3crms  = fWidthQv3[VZEROcentralityBin][0][0];
    Double_t Qy3cmean = fMeanQv3[VZEROcentralityBin][0][1];
    Double_t Qy3crms  = fWidthQv3[VZEROcentralityBin][0][1];	

    // update the weighted q-vectors with the re-centered values
    Qa2[0] = (Qa2[0] - Qx2amean)/Qx2arms;
    Qa2[1] = (Qa2[1] - Qy2amean)/Qy2arms;
    Qc2[0] = (Qc2[0] - Qx2cmean)/Qx2crms;
    Qc2[1] = (Qc2[1] - Qy2cmean)/Qy2crms;

    Qa3[0] = (Qa3[0] - Qx3amean)/Qx3arms;
    Qa3[1] = (Qa3[1] - Qy3amean)/Qy3arms;
    Qc3[0] = (Qc3[0] - Qx3cmean)/Qx3crms;
    Qc3[1] = (Qc3[0] - Qy3cmean)/Qy3crms;
}
//_____________________________________________________________________________
void AliAnalysisTaskJetV2::CalculateQvectorCombinedVZERO(Double_t Q2[2], Double_t Q3[2]) const
{
    // calculate calibrated q-vector of the combined vzeroa, vzeroc system
    // this is somewhat ugly as CalculateQvectorCombinedVZERO is called more than once per event
    // but for now it will have to do ...
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);

    // first step: retrieve the q-vectors component-wise per vzero detector
    Double_t QA2[] = {-999., -999.};
    Double_t QA3[] = {-999., -999.};
    Double_t QC2[] = {-999., -999.};
    Double_t QC3[] = {-999., -999.};
    CalculateQvectorVZERO(QA2, QA3, QC2, QC3);

    // get cache index and retrieve the chi weights for this centrality
    Int_t VZEROcentralityBin(GetVZEROCentralityBin());
    Double_t chi2A(fChi2A->At(VZEROcentralityBin));
    Double_t chi2C(fChi2C->At(VZEROcentralityBin));
    Double_t chi3A(fChi3A->At(VZEROcentralityBin));
    Double_t chi3C(fChi3C->At(VZEROcentralityBin));

    // combine the vzera and vzeroc signal
    Q2[0] = chi2A*chi2A*QA2[0]+chi2C*chi2C*QC2[0];
    Q2[1] = chi2A*chi2A*QA2[1]+chi2C*chi2C*QC2[1];
    Q3[0] = chi3A*chi3A*QA3[0]+chi3C*chi3C*QC3[0];
    Q3[1] = chi3A*chi3A*QC3[1]+chi3C*chi3C*QC3[1];
}
//_____________________________________________________________________________
void AliAnalysisTaskJetV2::CalculateRandomCone(Float_t &pt, Float_t &eta, Float_t &phi, 
        AliParticleContainer* tracksCont, AliClusterContainer* clusterCont, AliEmcalJet* jet) const
{
    // get a random cone
    if(fDebug > 1) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    pt = 0; eta = 0; phi = 0;
    Float_t etaJet(999), phiJet(999), dJet(999);        // no jet: same as jet very far away
    if(jet) { // if a leading jet is given, use its kinematic properties to exclude it
        etaJet = jet->Eta();
        phiJet = jet->Phi();
    }
    // the random cone acceptance has to equal the jet acceptance
    // this also insures safety when runnnig on the semi-good tpc runs for 11h data,
    // where jet acceptance is adjusted to reduced acceptance - hence random cone acceptance as well
    Float_t minPhi(GetJetContainer()->GetJetPhiMin()), maxPhi(GetJetContainer()->GetJetPhiMax());
    if(maxPhi > TMath::TwoPi()) maxPhi = TMath::TwoPi();
    if(minPhi < 0 ) minPhi = 0.;
    // construct a random cone and see if it's far away enough from the leading jet
    Int_t attempts(1000);
    while(kTRUE) {
        attempts--;
        eta = gRandom->Uniform(GetJetContainer()->GetJetEtaMin(), GetJetContainer()->GetJetEtaMax());
        phi = gRandom->Uniform(minPhi, maxPhi);

        dJet = TMath::Sqrt((etaJet-eta)*(etaJet-eta)+(phiJet-phi)*(phiJet-phi));
        if(dJet > fMinDisanceRCtoLJ) break;
        else if (attempts == 0) {
            printf(" > No random cone after 1000 tries, giving up ... !\n");
            return;
        }
    }
    // get the charged energy (if tracks are provided)
    if(tracksCont) {
        AliVParticle* track = tracksCont->GetNextAcceptParticle(0);
        while(track) {
            Float_t etaTrack(track->Eta()), phiTrack(track->Phi());
            // get distance from cone
            if(TMath::Abs(phiTrack-phi) > TMath::Abs(phiTrack - phi + TMath::TwoPi())) phiTrack+=TMath::TwoPi();
            if(TMath::Abs(phiTrack-phi) > TMath::Abs(phiTrack - phi - TMath::TwoPi())) phiTrack-=TMath::TwoPi();
            if(TMath::Sqrt(TMath::Abs((etaTrack-eta)*(etaTrack-eta)+(phiTrack-phi)*(phiTrack-phi))) <= GetJetRadius()) pt += track->Pt();
            track = tracksCont->GetNextAcceptParticle();
        }
    }
    // get the neutral energy (if clusters are provided)
    if(clusterCont) {
        AliVCluster* cluster = clusterCont->GetNextAcceptCluster(0);
        while(cluster) {
            TLorentzVector momentum;
            cluster->GetMomentum(momentum, const_cast<Double_t*>(fVertex));
            Float_t etaClus(momentum.Eta()), phiClus(momentum.Phi());
            // get distance from cone
            if(TMath::Abs(phiClus-phi) > TMath::Abs(phiClus - phi + TMath::TwoPi())) phiClus+=TMath::TwoPi();
            if(TMath::Abs(phiClus-phi) > TMath::Abs(phiClus - phi - TMath::TwoPi())) phiClus-=TMath::TwoPi();
            if(TMath::Sqrt(TMath::Abs((etaClus-eta)*(etaClus-eta)+(phiClus-phi)*(phiClus-phi))) <= GetJetRadius()) pt += momentum.Pt();
            cluster = clusterCont->GetNextAcceptCluster();
        }
    }
}
//_____________________________________________________________________________
Double_t AliAnalysisTaskJetV2::CalculateQC2(Int_t harm) {
    // get the second order q-cumulant, a -999 return will be caught in the qa routine of CorrectRho
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    Double_t reQ(0), imQ(0), modQ(0), M11(0), M(0);
    if(fUsePtWeight) {  // for the weighted 2-nd order q-cumulant
        QCnQnk(harm, 1, reQ, imQ);      // get the weighted 2-nd order q-vectors
        modQ = reQ*reQ+imQ*imQ;         // get abs Q-squared
        M11 = QCnM11();                 // equals S2,1 - S1,2
        return (M11 > 0) ? ((modQ - QCnS(1,2))/M11) : -999;
    } // else return the non-weighted 2-nd order q-cumulant
    QCnQnk(harm, 0, reQ, imQ);          // get the non-weighted 2-nd order q-vectors
    modQ = reQ*reQ+imQ*imQ;             // get abs Q-squared
    M = QCnM();
    return (M > 1) ? (modQ - M)/(M*(M-1)) : -999;
}
//_____________________________________________________________________________
Double_t AliAnalysisTaskJetV2::CalculateQC4(Int_t harm) {
    // get the fourth order q-cumulant, a -999 return will be caught in the qa routine of CorrectRho
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    Double_t reQn1(0), imQn1(0), reQ2n2(0), imQ2n2(0), reQn3(0), imQn3(0), M1111(0), M(0);
    Double_t a(0), b(0), c(0), d(0), e(0), f(0), g(0);  // terms of the calculation
    if(fUsePtWeight) {  // for the weighted 4-th order q-cumulant
        QCnQnk(harm, 1, reQn1, imQn1);
        QCnQnk(harm*2, 2, reQ2n2, imQ2n2);
        QCnQnk(harm, 3, reQn3, imQn3);
        // fill in the terms ...
        a = (reQn1*reQn1+imQn1*imQn1)*(reQn1*reQn1+imQn1*imQn1);
        b = reQ2n2*reQ2n2 + imQ2n2*imQ2n2;
        c = -2.*(reQ2n2*reQn1*reQn1-reQ2n2*imQn1*imQn1+2.*imQ2n2*reQn1*imQn1);
        d = 8.*(reQn3*reQn1+imQn3*imQn1);
        e = -4.*QCnS(1,2)*(reQn1*reQn1+imQn1*imQn1);
        f = -6.*QCnS(1,4);
        g = 2.*QCnS(2,2);
        M1111 = QCnM1111();
        return (M1111 > 0) ? (a+b+c+d+e+f+g)/M1111 : -999;
    }   // else return the unweighted case
    Double_t reQn(0), imQn(0), reQ2n(0), imQ2n(0);
    QCnQnk(harm, 0, reQn, imQn);
    QCnQnk(harm*2, 0, reQ2n, imQ2n);
    // fill in the terms ...
    M = QCnM();
    if(M < 4) return -999;
    a = (reQn*reQn+imQn*imQn)*(reQn*reQn+imQn*imQn);
    b = reQ2n*reQ2n + imQ2n*imQ2n;
    c = -2.*(reQ2n*reQn*reQn-reQ2n*imQn*imQn+2.*imQ2n*reQn*imQn);
    e = -4.*(M-2)*(reQn*reQn+imQn*imQn);
    f = 2.*M*(M-3);
    return (a+b+c+e+f)/(M*(M-1)*(M-2)*(M-3));
}
//_____________________________________________________________________________
void AliAnalysisTaskJetV2::QCnQnk(Int_t n, Int_t k, Double_t &reQ, Double_t &imQ) {
    // get the weighted n-th order q-vector, pass real and imaginary part as reference
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    if(!fTracks) return;
    fNAcceptedTracksQCn = 0;
    Int_t iTracks(fTracks->GetEntriesFast());
    for(Int_t iTPC(0); iTPC < iTracks; iTPC++) {
        AliVTrack* track = static_cast<AliVTrack*>(fTracks->At(iTPC));
        if(!PassesCuts(track) || track->Pt() < fSoftTrackMinPt || track->Pt() > fSoftTrackMaxPt) continue;
        fNAcceptedTracksQCn++;
        // for the unweighted case, k equals zero and the weight doesn't contribute to the equation below
        reQ += TMath::Power(track->Pt(), k) * TMath::Cos(((double)n)*track->Phi());
        imQ += TMath::Power(track->Pt(), k) * TMath::Sin(((double)n)*track->Phi());
    }
}
//_____________________________________________________________________________
void AliAnalysisTaskJetV2::QCnDiffentialFlowVectors(
        TClonesArray* pois, TArrayD* ptBins, Bool_t vpart, Double_t* repn, Double_t* impn, 
        Double_t *mp, Double_t *reqn, Double_t *imqn, Double_t* mq, Int_t n) 
{
     if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
   // get  unweighted differential flow vectors
    Int_t iPois(pois->GetEntriesFast());
    if(vpart) {
        for(Int_t i(0); i < iPois; i++) {
            for(Int_t ptBin(0); ptBin < ptBins->GetSize()-1; ptBin++) {
                AliVTrack* poi = static_cast<AliVTrack*>(pois->At(i));
                if(PassesCuts(poi)) {
                    if(poi->Pt() >= ptBins->At(ptBin) && poi->Pt() < ptBins->At(ptBin+1)) {
                            // fill the flow vectors assuming that all poi's are in the rp selection (true by design)
                            repn[ptBin]+=TMath::Cos(((double)n)*poi->Phi());
                            impn[ptBin]+=TMath::Sin(((double)n)*poi->Phi());
                            mp[ptBin]++;
                            reqn[ptBin]+=TMath::Cos(((double)n)*poi->Phi());
                            imqn[ptBin]+=TMath::Sin(((double)n)*poi->Phi());
                            mq[ptBin]++;
                    }
                }
            }
        }
    } else {
        for(Int_t i(0); i < iPois; i++) {
            for(Int_t ptBin(0); ptBin < ptBins->GetSize()-1; ptBin++) {
                AliEmcalJet* poi = static_cast<AliEmcalJet*>(pois->At(i));
                if(PassesCuts(poi)) {    
                    Double_t pt(poi->Pt()-poi->Area()*fLocalRho->GetLocalVal(poi->Phi(), GetJetContainer()->GetJetRadius(), fLocalRho->GetVal()));
                    if(pt >= ptBins->At(ptBin) && pt < ptBins->At(ptBin+1)) {    
                            repn[ptBin]+=TMath::Cos(((double)n)*poi->Phi());
                            impn[ptBin]+=TMath::Sin(((double)n)*poi->Phi());
                            mp[ptBin]++;        // qn isn't filled, no overlap between poi's and rp's
                    }
                }
            }
        }
    }
}
//_____________________________________________________________________________
Double_t AliAnalysisTaskJetV2::QCnS(Int_t i, Int_t j) {
    // get the weighted ij-th order autocorrelation correction
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    if(!fTracks || i <= 0 || j <= 0) return -999;
    Int_t iTracks(fTracks->GetEntriesFast());
    Double_t Sij(0);
    for(Int_t iTPC(0); iTPC < iTracks; iTPC++) {
        AliVTrack* track = static_cast<AliVTrack*>(fTracks->At(iTPC));
        if(!PassesCuts(track) || track->Pt() < fSoftTrackMinPt || track->Pt() > fSoftTrackMaxPt) continue;
        Sij+=TMath::Power(track->Pt(), j);
    }
    return TMath::Power(Sij, i);
}
//_____________________________________________________________________________
Double_t AliAnalysisTaskJetV2::QCnM() {
    // get multiplicity for unweighted q-cumulants. function QCnQnk should be called first
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    return (Double_t) fNAcceptedTracksQCn;
}
//_____________________________________________________________________________
Double_t AliAnalysisTaskJetV2::QCnM11() {
    // get multiplicity weights for the weighted two particle cumulant
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    return (QCnS(2,1) - QCnS(1,2));
}
//_____________________________________________________________________________
Double_t AliAnalysisTaskJetV2::QCnM1111() {
    // get multiplicity weights for the weighted four particle cumulant
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    return (QCnS(4,1)-6*QCnS(1,2)*QCnS(2,1)+8*QCnS(1,3)*QCnS(1,1)+3*QCnS(2,2)-6*QCnS(1,4));
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskJetV2::QCnRecovery(Double_t psi2, Double_t psi3) {
    // decides how to deal with the situation where c2 or c3 is negative 
    // returns kTRUE depending on whether or not a modulated rho is used for the jet background
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    if(TMath::AreEqualAbs(fFitModulation->GetParameter(3), .0, 1e-10) && TMath::AreEqualAbs(fFitModulation->GetParameter(7), .0,1e-10)) {
        fFitModulation->SetParameter(7, 0);
        fFitModulation->SetParameter(3, 0);
        fFitModulation->SetParameter(0, fLocalRho->GetVal());
        return kTRUE;   // v2 and v3 have physical null values
    }
    switch (fQCRecovery) {
        case kFixedRho : {      // roll back to the original rho
           fFitModulation->SetParameter(7, 0);
           fFitModulation->SetParameter(3, 0);
           fFitModulation->SetParameter(0, fLocalRho->GetVal());
           return kFALSE;       // rho is forced to be fixed
        }
        case kNegativeVn : {
           Double_t c2(fFitModulation->GetParameter(3));
           Double_t c3(fFitModulation->GetParameter(7));
           if( c2 < 0 ) c2 = -1.*TMath::Sqrt(-1.*c2);
           if( c3 < 0 ) c3 = -1.*TMath::Sqrt(-1.*c3);
           fFitModulation->SetParameter(3, c2);
           fFitModulation->SetParameter(7, c3);
           return kTRUE;        // is this a physical quantity ?
        }
        case kTryFit : {
           fitModulationType tempType(fFitModulationType);  // store temporarily
           fFitModulationType = kCombined;
           fFitModulation->SetParameter(7, 0);
           fFitModulation->SetParameter(3, 0);
           Bool_t pass(CorrectRho(psi2, psi3));         // do the fit and all quality checks
           fFitModulationType = tempType;               // roll back for next event
           return pass;
        }
        default : return kFALSE;
    }
    return kFALSE;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskJetV2::CorrectRho(Double_t psi2, Double_t psi3) 
{
    // get rho' -> rho(phi)
    // two routines are available, both can be used with or without pt weights
    //  [1] get vn from q-cumulants or as an integrated value from a user supplied histogram
    //      in case of cumulants, both cumulants and vn values are stored. in both cases, v2 and v3
    //      are expected. a check is performed to see if rho has no negative local minimum
    //      for full description, see Phys. Rev. C 83, 044913
    //      since the cn distribution has negative values, vn = sqrt(cn) can be imaginary sometimes
    //      in this case one can either roll back to the 'original' rixed rho, do a fit for vn or take use
    //      vn = - sqrt(|cn|) 
    //  [2] fitting a fourier expansion to the de/dphi distribution
    //      the fit can be done with either v2, v3 or a combination.
    //      in all cases, a cut can be made on the p-value of the chi-squared value of the fit
    //      and a check can be performed to see if rho has no negative local minimum
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    Int_t freeParams(2);                // free parameters of the fit (for NDF) 
    switch (fFitModulationType) {       // for approaches where no fitting is required
        case kQC2 : {
            fFitModulation->FixParameter(4, psi2); 
            fFitModulation->FixParameter(6, psi3);
            fFitModulation->FixParameter(3, CalculateQC2(2));   // set here with cn, vn = sqrt(cn)
            fFitModulation->FixParameter(7, CalculateQC2(3));
            // first fill the histos of the raw cumulant distribution
            if (fUsePtWeight) { // use weighted weights
                Double_t dQCnM11 = (fNoEventWeightsForQC) ? 1. : QCnM11();
                fProfV2Cumulant->Fill(fCent, fFitModulation->GetParameter(3), dQCnM11);
                fProfV3Cumulant->Fill(fCent, fFitModulation->GetParameter(7), dQCnM11);
            } else {
                Double_t dQCnM = (fNoEventWeightsForQC) ? 2. : QCnM();
                fProfV2Cumulant->Fill(fCent, fFitModulation->GetParameter(3), dQCnM*(dQCnM-1));
                fProfV3Cumulant->Fill(fCent, fFitModulation->GetParameter(7), dQCnM*(dQCnM-1));
            }
            // then see if one of the cn value is larger than zero and vn is readily available
            if(fFitModulation->GetParameter(3) > 0 && fFitModulation->GetParameter(7) > 0) {
                fFitModulation->FixParameter(3, TMath::Sqrt(fFitModulation->GetParameter(3)));
                fFitModulation->FixParameter(7, TMath::Sqrt(fFitModulation->GetParameter(7)));
            } else if (!QCnRecovery(psi2, psi3)) return kFALSE;  // try to recover the cumulant, this will set v2 and v3
            if(fFitModulation->GetMinimum(0, TMath::TwoPi()) < 0) {  // general check 
                fFitModulation->SetParameter(7, 0);
                fFitModulation->SetParameter(3, 0);
                fFitModulation->SetParameter(0, fLocalRho->GetVal());
                return kFALSE;
            }
            return kTRUE;
        } break;
        case kQC4 : {
            fFitModulation->FixParameter(4, psi2); 
            fFitModulation->FixParameter(6, psi3);
            fFitModulation->FixParameter(3, CalculateQC4(2));   // set here with cn, vn = sqrt(cn)
            fFitModulation->FixParameter(7, CalculateQC4(3));
            // first fill the histos of the raw cumulant distribution
            if (fUsePtWeight) { // use weighted weights
                fProfV2Cumulant->Fill(fCent, fFitModulation->GetParameter(3)/*, QCnM1111()*/);
                fProfV3Cumulant->Fill(fCent, fFitModulation->GetParameter(7)/*, QCnM1111()*/);
            } else {
                fProfV2Cumulant->Fill(fCent, fFitModulation->GetParameter(3)/*, QCnM1111()*/);
                fProfV3Cumulant->Fill(fCent, fFitModulation->GetParameter(7)/*, QCnM1111()*/);
            }
            // then see if one of the cn value is larger than zero and vn is readily available
            if(fFitModulation->GetParameter(3) > 0 && fFitModulation->GetParameter(7) > 0) {
                fFitModulation->FixParameter(3, TMath::Sqrt(fFitModulation->GetParameter(3)));
                fFitModulation->FixParameter(7, TMath::Sqrt(fFitModulation->GetParameter(7)));
            } else if (!QCnRecovery(psi2, psi3)) return kFALSE;  // try to recover the cumulant, this will set v2 and v3
            if(fFitModulation->GetMinimum(0, TMath::TwoPi()) < 0) {  // general check 
                fFitModulation->SetParameter(7, 0);
                fFitModulation->SetParameter(3, 0);
                fFitModulation->SetParameter(0, fLocalRho->GetVal());
                return kFALSE;
            }
        } break;
        case kIntegratedFlow : {
            // use v2 and v3 values from an earlier iteration over the data
            fFitModulation->FixParameter(3, fUserSuppliedV2->GetBinContent(fUserSuppliedV2->GetXaxis()->FindBin(fCent)));
            fFitModulation->FixParameter(4, psi2);
            fFitModulation->FixParameter(6, psi3);
            fFitModulation->FixParameter(7, fUserSuppliedV3->GetBinContent(fUserSuppliedV3->GetXaxis()->FindBin(fCent)));
            if(fFitModulation->GetMinimum(0, TMath::TwoPi()) < 0) { 
                fFitModulation->SetParameter(7, 0);
                fFitModulation->SetParameter(3, 0);
                fFitModulation->SetParameter(0, fLocalRho->GetVal());
                return kFALSE;
            }
            return kTRUE;
        }
        default : break;
    }
    TString detector("");
    switch (fDetectorType) {
        case kTPC : detector+="TPC";
            break;
        case kVZEROA : detector+="VZEROA";
            break;
        case kVZEROC : detector+="VZEROC";
            break;
        case kVZEROComb : detector+="VZEROComb";
            break; 
        case kFixedEP : detector+="FixedEP";
            break;
        default: break;
    }
    Int_t iTracks(fTracks->GetEntriesFast());
    Double_t excludeInEta = -999;
    Double_t excludeInPhi = -999;
    Double_t excludeInPt  = -999;
    if(iTracks <= 0 || fLocalRho->GetVal() <= 0 ) return kFALSE;   // no use fitting an empty event ...
    if(fExcludeLeadingJetsFromFit > 0 ) {
        if(fLeadingJet) {
            excludeInEta = fLeadingJet->Eta();
            excludeInPhi = fLeadingJet->Phi();
            excludeInPt = fLeadingJet->Pt();
        }
    }
    // check the acceptance of the track selection that will be used
    // if one uses e.g. semi-good tpc tracks, accepance in phi is reduced to 0 < phi < 4
    // the defaults (-10 < phi < 10) which accept all, are then overwritten
    Double_t lowBound(0.), upBound(TMath::TwoPi());     // bounds for fit
    if(GetParticleContainer()->GetParticlePhiMin() > lowBound) lowBound = GetParticleContainer()->GetParticlePhiMin();
    if(GetParticleContainer()->GetParticlePhiMax() < upBound) upBound = GetParticleContainer()->GetParticlePhiMax();
    fHistSwap->Reset(); // clear the histogram
    TH1F _tempSwap;     // on stack for quick access
    TH1F _tempSwapN;    // on stack for quick access, bookkeeping histogram
    if(fRebinSwapHistoOnTheFly) {
        if(fNAcceptedTracks < 49) fNAcceptedTracks = 49;       // avoid aliasing effects
        _tempSwap = TH1F("_tempSwap", "_tempSwap", TMath::CeilNint(TMath::Sqrt(fNAcceptedTracks)), lowBound, upBound);
        if(fUsePtWeightErrorPropagation) _tempSwapN = TH1F("_tempSwapN", "_tempSwapN", TMath::CeilNint(TMath::Sqrt(fNAcceptedTracks)), lowBound, upBound);
        if(fUsePtWeight) _tempSwap.Sumw2();
    }
    else _tempSwap = *fHistSwap;         // now _tempSwap holds the desired histo
    // non poissonian error when using pt weights
    Double_t totalpts(0.), totalptsquares(0.), totalns(0.);
    for(Int_t i(0); i < iTracks; i++) {
        AliVTrack* track = static_cast<AliVTrack*>(fTracks->At(i));
        if(fExcludeLeadingJetsFromFit > 0 &&( (TMath::Abs(track->Eta() - excludeInEta) < GetJetContainer()->GetJetRadius()*fExcludeLeadingJetsFromFit ) || (TMath::Abs(track->Eta()) - GetJetContainer()->GetJetRadius() - GetJetContainer()->GetJetEtaMax() ) > 0 )) continue;
        if(!PassesCuts(track) || track->Pt() > fSoftTrackMaxPt || track->Pt() < fSoftTrackMinPt) continue;
        if(fUsePtWeight) {
            _tempSwap.Fill(track->Phi(), track->Pt());
            if(fUsePtWeightErrorPropagation) {
                totalpts += track->Pt();
                totalptsquares += track->Pt()*track->Pt();
		totalns += 1;
                _tempSwapN.Fill(track->Phi());
            }
        }
        else _tempSwap.Fill(track->Phi());
    }
    if(fUsePtWeight && fUsePtWeightErrorPropagation) {
        // in the case of pt weights overwrite the poissonian error estimate which is assigned by root by a more sophisticated appraoch
        // the assumption here is that the bin error will be dominated by the uncertainty in the mean pt in a bin and in the uncertainty
        // of the number of tracks in a bin, the first of which will be estimated from the sample standard deviation of all tracks in the 
        // event, for the latter use a poissonian estimate. the two contrubitions are assumed to be uncorrelated
        if(totalns < 2) return kFALSE; // not one track passes the cuts > 2 avoids possible division by 0 later on
        for(Int_t l = 0; l < _tempSwap.GetNbinsX(); l++) {
            if(_tempSwapN.GetBinContent(l+1) == 0) {
                _tempSwap.SetBinContent(l+1,0);
                _tempSwap.SetBinError(l+1,0);
            }
            else {
                Double_t vartimesnsq = totalptsquares*totalns - totalpts*totalpts;
                Double_t variance = vartimesnsq/(totalns*(totalns-1.));
                Double_t SDOMSq = variance / _tempSwapN.GetBinContent(l+1);
                Double_t SDOMSqOverMeanSq = SDOMSq * _tempSwapN.GetBinContent(l+1) * _tempSwapN.GetBinContent(l+1) / (_tempSwapN.GetBinContent(l+1) * _tempSwapN.GetBinContent(l+1));
                Double_t poissonfrac = 1./_tempSwapN.GetBinContent(l+1);
                Double_t vartotalfrac = SDOMSqOverMeanSq + poissonfrac;
                Double_t vartotal = vartotalfrac * _tempSwap.GetBinContent(l+1) * _tempSwap.GetBinContent(l+1);
                if(vartotal > 0.0001) _tempSwap.SetBinError(l+1,TMath::Sqrt(vartotal));
                else {
                    _tempSwap.SetBinContent(l+1,0);
                    _tempSwap.SetBinError(l+1,0);
                }
            }
        }
    }
    fFitModulation->SetParameter(0, fLocalRho->GetVal());
    switch (fFitModulationType) {
        case kNoFit : { 
            fFitModulation->FixParameter(0, fLocalRho->GetVal() ); 
            freeParams = 0;
        } break;
        case kV2 : { 
            fFitModulation->FixParameter(4, psi2); 
            freeParams = 1;
        } break;
        case kV3 : { 
            fFitModulation->FixParameter(4, psi3); 
            freeParams = 1;
        } break;
        case kCombined : {
            fFitModulation->FixParameter(4, psi2); 
            fFitModulation->FixParameter(6, psi3);
            freeParams = 2;
        } break;
        case kFourierSeries : {
            // in this approach, an explicit calculation will be made of vn = sqrt(xn^2+yn^2)
            // where x[y] = Integrate[r(phi)cos[sin](n phi)dphi, 0, 2pi]
            Double_t cos2(0), sin2(0), cos3(0), sin3(0), sumPt(0);
            for(Int_t i(0); i < iTracks; i++) {
                AliVTrack* track = static_cast<AliVTrack*>(fTracks->At(i));
                if(!PassesCuts(track) || track->Pt() > fSoftTrackMaxPt || track->Pt() < fSoftTrackMinPt) continue;
                sumPt += track->Pt();
                cos2 += track->Pt()*TMath::Cos(2*PhaseShift(track->Phi()-psi2)); 
                sin2 += track->Pt()*TMath::Sin(2*PhaseShift(track->Phi()-psi2));
                cos3 += track->Pt()*TMath::Cos(3*PhaseShift(track->Phi()-psi3)); 
                sin3 += track->Pt()*TMath::Sin(3*PhaseShift(track->Phi()-psi3));
            }
            fFitModulation->SetParameter(3, TMath::Sqrt(cos2*cos2+sin2*sin2)/fLocalRho->GetVal());
            fFitModulation->SetParameter(4, psi2);
            fFitModulation->SetParameter(6, psi3);
            fFitModulation->SetParameter(7, TMath::Sqrt(cos3*cos3+sin3*sin3)/fLocalRho->GetVal());
        } break;
        default : break;
    }
    if(fRunToyMC) {
        // toy mc, just here to check procedure, azimuthal profile is filled from hypothesis so p-value distribution should be flat
        Int_t _bins = _tempSwap.GetXaxis()->GetNbins();
        TF1* _tempFit = new TF1("temp_fit_kCombined", "[0]*([1]+[2]*([3]*TMath::Cos([2]*(x-[4]))+[7]*TMath::Cos([5]*(x-[6]))))", 0, TMath::TwoPi());
        _tempFit->SetParameter(0, fFitModulation->GetParameter(0));       // normalization
        _tempFit->SetParameter(3, 0.1);      // v2
        _tempFit->FixParameter(1, 1.);       // constant
        _tempFit->FixParameter(2, 2.);       // constant
        _tempFit->FixParameter(5, 3.);       // constant
        _tempFit->FixParameter(4, fFitModulation->GetParameter(4));
        _tempFit->FixParameter(6, fFitModulation->GetParameter(6));
        _tempFit->SetParameter(7, 0.1);      // v3
        _tempSwap.Reset();                   // rese bin content
        for(int _binsI = 0; _binsI < _bins*_bins; _binsI++)  _tempSwap.Fill(_tempFit->GetRandom());
    }
    _tempSwap.Fit(fFitModulation, fFitModulationOptions.Data(), "", lowBound, upBound);
    // the quality of the fit is evaluated from 1 - the cdf of the chi square distribution
    // three methods are available, all with their drawbacks. all are stored, one is selected to do the cut
    Int_t NDF(_tempSwap.GetXaxis()->GetNbins()-freeParams);
    if(NDF == 0 || (float)NDF <= 0.) return kFALSE;
    Double_t CDF(1.-ChiSquareCDF(NDF, ChiSquare(_tempSwap, fFitModulation)));
    Double_t CDFROOT(1.-ChiSquareCDF(NDF, fFitModulation->GetChisquare()));
    Double_t CDFKolmogorov(KolmogorovTest(_tempSwap, fFitModulation));
    // fill the values and centrality correlation (redundant but easy on the eyes)
    fHistPvalueCDF->Fill(CDF);
    fHistPvalueCDFCent->Fill(fCent, CDF);
    fHistPvalueCDFROOT->Fill(CDFROOT);
    fHistPvalueCDFROOTCent->Fill(fCent, CDFROOT);
    fHistKolmogorovTest->Fill(CDFKolmogorov);
    fHistChi2ROOTCent->Fill(fCent, fFitModulation->GetChisquare()/((float)NDF));
    fHistChi2Cent->Fill(fCent, ChiSquare(_tempSwap, fFitModulation)/((float)NDF));
    fHistKolmogorovTestCent->Fill(fCent, CDFKolmogorov);
    fHistPChi2Root->Fill(CDFROOT, fFitModulation->GetChisquare()/((float)NDF));
    fHistPChi2->Fill(CDF, ChiSquare(_tempSwap, fFitModulation)/((float)NDF));
    fHistPKolmogorov->Fill(CDF, CDFKolmogorov);

    // variable CDF is used for making cuts, so we fill it with the selected p-value
    switch (fFitGoodnessTest) {
        case kChi2ROOT : {
            CDF = CDFROOT; 
        } break;
        case kChi2Poisson : break;      // CDF is already CDF
        case kKolmogorov : {
            CDF = CDFKolmogorov; 
        } break;
        default: break;
    }

    if(fFitControl) {
        // as an additional quality check, see if fitting a control fit has a higher significance
        _tempSwap.Fit(fFitControl, fFitModulationOptions.Data(), "", lowBound, upBound);
        Double_t CDFControl(-1.);
        switch (fFitGoodnessTest) {
            case kChi2ROOT : {
                CDFControl = 1.-ChiSquareCDF(fFitControl->GetNDF(), fFitModulation->GetChisquare());
            } break;
            case kChi2Poisson : {
                CDFControl = 1.-ChiSquareCDF(fFitControl->GetNDF(), ChiSquare(_tempSwap, fFitModulation));
            } break;
            case kKolmogorov : {
                CDFControl = KolmogorovTest(_tempSwap, fFitControl); 
            } break;
            default: break;
        }
        if(CDFControl > CDF) {
            CDF = -1.; // control fit is more significant, so throw out the 'old' fit
            fHistRhoStatusCent->Fill(fCent, -1);
        }
    }
    if(CDF >= fMinPvalue && CDF <= fMaxPvalue && ( fFitModulation->GetMinimum(0, TMath::TwoPi()) > 0)) {       
        // fit quality. not that although with limited acceptance the fit is performed on just
        // part of phase space, the requirement that energy desntiy is larger than zero is applied
        // to the FULL spectrum
        fHistRhoStatusCent->Fill(fCent, 0.);
        // for LOCAL didactic purposes, save the  best and the worst fits
        // this routine can produce a lot of output histograms (it's not memory 'safe') and will not work on GRID 
        // since the output will become unmergeable (i.e. different nodes may produce conflicting output)
        switch (fRunModeType) {
            case kLocal : {
                if(fRandom->Uniform(0, 100) > fPercentageOfFits) break;
                static Int_t didacticCounterBest(0);
                TProfile* didacticProfile = (TProfile*)_tempSwap.Clone(Form("Fit_%i_1-CDF_%.3f_cen_%i_%s", didacticCounterBest, CDF, fInCentralitySelection, detector.Data()));
                TF1* didacticFit = (TF1*)fFitModulation->Clone(Form("fit_%i_CDF_%.3f_cen_%i_%s", didacticCounterBest, CDF, fInCentralitySelection, detector.Data()));
                switch(fFitModulationType) { 
                    case kCombined : {
                        // to make a nice picture also plot the separate components (v2 and v3) of the fit
                        // only done for cobined fit where there are actually components to split ...
                        TF1* v0(new TF1("dfit_kV2", "[0]", 0, TMath::TwoPi()));
                        v0->SetParameter(0, didacticFit->GetParameter(0));        // normalization
                        v0->SetLineColor(kMagenta);
                        v0->SetLineStyle(7);
                        didacticProfile->GetListOfFunctions()->Add(v0);
                        TF1* v2(new TF1("dfit_kV2", "[0]*([1]+[2]*[3]*TMath::Cos([2]*(x-[4])))", 0, TMath::TwoPi()));
                        v2->SetParameter(0, didacticFit->GetParameter(0));        // normalization
                        v2->SetParameter(3, didacticFit->GetParameter(3));        // v2
                        v2->FixParameter(1, 1.);        // constant
                        v2->FixParameter(2, 2.);        // constant
                        v2->FixParameter(4, didacticFit->GetParameter(4));        // psi2
                        v2->SetLineColor(kGreen);
                        didacticProfile->GetListOfFunctions()->Add(v2);
                        TF1* v3(new TF1("dfit_kV3", "[0]*([1]+[2]*[3]*TMath::Cos([5]*(x-[4])))", 0, TMath::TwoPi()));
                        v3->SetParameter(0, didacticFit->GetParameter(0));        // normalization
                        v3->SetParameter(3, didacticFit->GetParameter(7));        // v3
                        v3->FixParameter(1, 1.);        // constant
                        v3->FixParameter(2, 2.);        // constant
                        v3->FixParameter(4, didacticFit->GetParameter(6));        // psi3
                        v3->FixParameter(5, 3.);        // constant
                        v3->SetLineColor(kCyan);
                        didacticProfile->GetListOfFunctions()->Add(v3);
                    }
                    default : break;
                }
                didacticProfile->GetListOfFunctions()->Add(didacticFit);
                didacticProfile->GetYaxis()->SetTitle("#frac{d #sum #it{p}_{T}}{d #varphi} [GeV/#it{c}]");
                didacticProfile->GetXaxis()->SetTitle("#varphi");
                fOutputListGood->Add(didacticProfile);
                didacticCounterBest++;
                TH2F* didacticSurface = BookTH2F(Form("surface_%s", didacticProfile->GetName()), "#phi", "#eta", 50, 0, TMath::TwoPi(), 50, -1, 1, -1, kFALSE);
                for(Int_t i(0); i < iTracks; i++) {
                    AliVTrack* track = static_cast<AliVTrack*>(fTracks->At(i));
                    if(PassesCuts(track)) {
                        if(fUsePtWeight) didacticSurface->Fill(track->Phi(), track->Eta(), track->Pt());
                        else didacticSurface->Fill(track->Phi(), track->Eta());
                    }
                }
                if(fExcludeLeadingJetsFromFit) {       // visualize the excluded region
                    TF2 *f2 = new TF2(Form("%s_LJ", didacticSurface->GetName()),"[0]*TMath::Gaus(x,[1],[2])*TMath::Gaus(y,[3],[4])", 0, TMath::TwoPi(), -1, 1);
                    f2->SetParameters(excludeInPt/3.,excludeInPhi,.1,excludeInEta,.1);
                    didacticSurface->GetListOfFunctions()->Add(f2);
                }
                fOutputListGood->Add(didacticSurface);
            } break;
            default : break;
        }
    } else {    // if the fit is of poor quality revert to the original rho estimate
        switch (fRunModeType) { // again see if we want to save the fit
            case kLocal : {
                static Int_t didacticCounterWorst(0);
                if(fRandom->Uniform(0, 100) > fPercentageOfFits) break;
                TProfile* didacticProfile = (TProfile*)_tempSwap.Clone(Form("Fit_%i_1-CDF_%.3f_cen_%i_%s", didacticCounterWorst, CDF, fInCentralitySelection, detector.Data() ));
                TF1* didacticFit = (TF1*)fFitModulation->Clone(Form("fit_%i_p_%.3f_cen_%i_%s", didacticCounterWorst, CDF, fInCentralitySelection, detector.Data()));
                didacticProfile->GetListOfFunctions()->Add(didacticFit);
                fOutputListBad->Add(didacticProfile);
                didacticCounterWorst++;
                } break;
            default : break;
        }
        switch (fFitModulationType) {
            case kNoFit : break;        // nothing to do
            case kCombined : fFitModulation->SetParameter(7, 0);        // no break
            case kFourierSeries : fFitModulation->SetParameter(7, 0);   // no break
            default : { // needs to be done if there was a poor fit
                 fFitModulation->SetParameter(3, 0);
                 fFitModulation->SetParameter(0, fLocalRho->GetVal());
            } break;
        }
        if(CDF > -.5) fHistRhoStatusCent->Fill(fCent, 1.);
        return kFALSE;  // return false if the fit is rejected
    }
    return kTRUE;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskJetV2::PassesCuts(AliVEvent* event)
{
    // event cuts
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    switch (fCollisionType) {
        case kJetFlowMC : {
            fInCentralitySelection = 0;
            return kTRUE;
    } break;
        default : break;
    }
    if(!event || !AliAnalysisTaskEmcal::IsEventSelected()) return kFALSE;
    if(TMath::Abs(InputEvent()->GetPrimaryVertex()->GetZ()) > 10.) return kFALSE;
    // aod and esd specific checks
    switch (fDataType) {
       case kESD: {
            AliESDEvent* esdEvent = static_cast<AliESDEvent*>(InputEvent());
            if( (!esdEvent) || (TMath::Abs(esdEvent->GetPrimaryVertexSPD()->GetZ() - esdEvent->GetPrimaryVertex()->GetZ()) > .5) ) return kFALSE; 
       } break;
       case kAOD: {
            AliAODEvent* aodEvent = static_cast<AliAODEvent*>(InputEvent());
            if( (!aodEvent) || (TMath::Abs(aodEvent->GetPrimaryVertexSPD()->GetZ() - aodEvent->GetPrimaryVertex()->GetZ()) > .5) ) return kFALSE; 
       } break;
       default: break;
    }
    fCent = InputEvent()->GetCentrality()->GetCentralityPercentile("V0M");
    if(fCent <= fCentralityClasses->At(0) || fCent >= fCentralityClasses->At(fCentralityClasses->GetSize()-1) || TMath::Abs(fCent-InputEvent()->GetCentrality()->GetCentralityPercentile("TRK")) > 5.) return kFALSE;
    // determine centrality class
    fInCentralitySelection = -1;
    for(Int_t i(0); i < fCentralityClasses->GetSize()-1; i++) {
        if(fCent >= fCentralityClasses->At(i) && fCent <= fCentralityClasses->At(1+i)) {
            fInCentralitySelection = i;
            break; }
    } 
    if(fInCentralitySelection<0) return kFALSE;     // should be null op
    // see if input containers are filled
    if(fTracks->GetEntries() < 1) return kFALSE;
    if(fRho->GetVal() <= 0 ) return kFALSE;
    if(fAnalysisType == AliAnalysisTaskJetV2::kFull && !fClusterCont) return kFALSE;
    return kTRUE;
}
//_____________________________________________________________________________
void AliAnalysisTaskJetV2::FillHistogramsAfterSubtraction(Double_t psi2, Double_t vzero[2][2], Double_t* vzeroComb, Double_t* tpc)
{
    // fill histograms 
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    FillTrackHistograms();
    if(fAnalysisType == AliAnalysisTaskJetV2::kFull) FillClusterHistograms();
    FillJetHistograms(psi2); 
    if(fFillQAHistograms) FillEventPlaneHistograms(vzero, vzeroComb, tpc);
    FillRhoHistograms();
    FillDeltaPtHistograms(psi2);
}
//_____________________________________________________________________________
void AliAnalysisTaskJetV2::FillTrackHistograms() const
{
    // fill track histograms
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    Int_t iTracks(fTracks->GetEntriesFast()), iAcceptedTracks(0);
    for(Int_t i(0); i < iTracks; i++) {
        AliVTrack* track = static_cast<AliVTrack*>(fTracks->At(i));
        if(!PassesCuts(track)) continue;
        iAcceptedTracks++;
        fHistPicoTrackPt[fInCentralitySelection]->Fill(track->Pt());
        if(fFillQAHistograms) FillQAHistograms(track);
    }
    fHistPicoTrackMult[fInCentralitySelection]->Fill(iAcceptedTracks);
}
//_____________________________________________________________________________
void AliAnalysisTaskJetV2::FillClusterHistograms() const
{
    // fill cluster histograms
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    if(!fClusterCont) return;
    Int_t iClusters(fClusterCont->GetNClusters());
    for(Int_t i(0); i < iClusters; i++) {
        AliVCluster* cluster = fClusterCont->GetCluster(i);
        if (!PassesCuts(cluster)) continue;
        TLorentzVector clusterLorentzVector;
        cluster->GetMomentum(clusterLorentzVector, const_cast<Double_t*>(fVertex));
        fHistClusterPt[fInCentralitySelection]->Fill(clusterLorentzVector.Pt());
        fHistClusterEtaPhi[fInCentralitySelection]->Fill(clusterLorentzVector.Eta(), clusterLorentzVector.Phi());
        fHistClusterEtaPhiWeighted[fInCentralitySelection]->Fill(clusterLorentzVector.Eta(), clusterLorentzVector.Phi(), clusterLorentzVector.Pt());
    }
    return;
}
//_____________________________________________________________________________
void AliAnalysisTaskJetV2::FillEventPlaneHistograms(Double_t vzero[2][2], Double_t* vzeroComb, Double_t* tpc) const
{
    // fill event plane histograms
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    fHistPsiControl->Fill(0.5, vzero[0][0]);    // vzero a psi2
    fHistPsiControl->Fill(1.5, vzero[1][0]);    // vzero c psi2
    fHistPsiControl->Fill(2.5, tpc[0]);         // tpc psi 2
    fHistPsiControl->Fill(5.5, vzero[0][1]);    // vzero a psi3
    fHistPsiControl->Fill(6.5, vzero[1][1]);    // vzero b psi3
    fHistPsiControl->Fill(7.5, tpc[1]);         // tpc psi 3
    fHistPsiVZEROA->Fill(vzero[0][0]);
    fHistPsiVZEROC->Fill(vzero[1][0]);
    fHistPsiVZERO->Fill(vzeroComb[0]);
    fHistPsiTPC->Fill(tpc[0]);
    fHistPsiSpread->Fill(0.5, TMath::Abs(vzero[0][0]-vzero[1][0]));
    fHistPsiSpread->Fill(1.5, TMath::Abs(vzero[0][0]-tpc[0]));
    fHistPsiSpread->Fill(2.5, TMath::Abs(vzero[1][0]-tpc[0]));
    // event plane vs centrality QA histo's to check recentering
    Double_t TRK(InputEvent()->GetCentrality()->GetCentralityPercentile("TRK"));
    Double_t V0M(InputEvent()->GetCentrality()->GetCentralityPercentile("V0M"));
    fHistPsiVZEROAV0M->Fill(V0M, vzero[0][0]);
    fHistPsiVZEROCV0M->Fill(V0M, vzero[1][0]);
    fHistPsiVZEROVV0M->Fill(V0M, vzeroComb[0]);
    fHistPsiTPCiV0M->Fill(V0M, tpc[0]);
    fHistPsiVZEROATRK->Fill(TRK, vzero[0][0]);
    fHistPsiVZEROCTRK->Fill(TRK, vzero[1][0]);
    fHistPsiVZEROTRK->Fill(TRK, vzeroComb[0]);
    fHistPsiTPCTRK->Fill(TRK, tpc[0]);
}
//_____________________________________________________________________________
void AliAnalysisTaskJetV2::FillRhoHistograms()
{
    // fill rho histograms
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    fHistRhoPackage[fInCentralitySelection]->Fill(fLocalRho->GetVal());    // save the rho estimate from the emcal jet package
    // get multiplicity FIXME inefficient
    Int_t iJets(fJets->GetEntriesFast());
    Double_t rho(fLocalRho->GetLocalVal(TMath::Pi(), TMath::Pi(), fLocalRho->GetVal()));
    fHistRho[fInCentralitySelection]->Fill(rho);
    fHistRhoVsMult->Fill(fTracks->GetEntries(), rho);
    fHistRhoVsCent->Fill(fCent, rho);
    for(Int_t i(0); i < iJets; i++) {
        AliEmcalJet* jet = static_cast<AliEmcalJet*>(fJets->At(i));
        if(!PassesCuts(jet)) continue;
        fHistRhoAVsMult->Fill(fTracks->GetEntries(), rho * jet->Area());
        fHistRhoAVsCent->Fill(fCent, rho * jet->Area());
    }
}
//_____________________________________________________________________________
void AliAnalysisTaskJetV2::FillDeltaPtHistograms(Double_t psi2) const
{
    // fill delta pt histograms
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    Int_t i(0);
    const Float_t areaRC = GetJetRadius()*GetJetRadius()*TMath::Pi();
    // we're retrieved the leading jet, now get a random cone
    for(i = 0; i < fMaxCones; i++) {
       Float_t pt(0), eta(0), phi(0);
       // get a random cone without constraints on leading jet position
       CalculateRandomCone(pt, eta, phi, fTracksCont, fClusterCont, 0x0);
       if(pt > 0) {
           if(fFillQAHistograms) fHistRCPhiEta[fInCentralitySelection]->Fill(phi, eta);
           fHistRhoVsRCPt[fInCentralitySelection]->Fill(pt, fLocalRho->GetLocalVal(phi, GetJetContainer()->GetJetRadius(), fLocalRho->GetVal())*areaRC);
           fHistRCPt[fInCentralitySelection]->Fill(pt);
           fHistDeltaPtDeltaPhi2[fInCentralitySelection]->Fill(PhaseShift(phi-psi2, 2.), pt - areaRC*fLocalRho->GetLocalVal(phi, GetJetContainer()->GetJetRadius(), fLocalRho->GetVal()));
           fHistDeltaPtDeltaPhi2Rho0[fInCentralitySelection]->Fill(PhaseShift(phi-psi2, 2.), pt - areaRC*fLocalRho->GetVal());

       }
       // get a random cone excluding leading jet area
       CalculateRandomCone(pt, eta, phi, fTracksCont, fClusterCont, fLeadingJet);
       if(pt > 0) {
           if(fFillQAHistograms) fHistRCPhiEtaExLJ[fInCentralitySelection]->Fill(phi, eta);
           fHistRhoVsRCPtExLJ[fInCentralitySelection]->Fill(pt, fLocalRho->GetLocalVal(phi, GetJetContainer()->GetJetRadius(), fLocalRho->GetVal())*areaRC);
           fHistRCPtExLJ[fInCentralitySelection]->Fill(pt);
           fHistDeltaPtDeltaPhi2ExLJ[fInCentralitySelection]->Fill(PhaseShift(phi-psi2, 2.), pt - areaRC*fLocalRho->GetLocalVal(phi, GetJetContainer()->GetJetRadius(), fLocalRho->GetVal()));
           fHistDeltaPtDeltaPhi2ExLJRho0[fInCentralitySelection]->Fill(PhaseShift(phi-psi2, 2.), pt - areaRC*fLocalRho->GetVal());
       }
    } 
}
//_____________________________________________________________________________
void AliAnalysisTaskJetV2::FillJetHistograms(Double_t psi2)
{
    // fill jet histograms
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    Int_t iJets(fJets->GetEntriesFast());
    for(Int_t i(0); i < iJets; i++) {
        AliEmcalJet* jet = static_cast<AliEmcalJet*>(fJets->At(i));
        if(PassesCuts(jet)) {
            Double_t pt(jet->Pt()), area(jet->Area()), eta(jet->Eta()), phi(jet->Phi());
            Double_t rho(fLocalRho->GetLocalVal(phi, GetJetContainer()->GetJetRadius(), fLocalRho->GetVal()));
            fHistJetPtRaw[fInCentralitySelection]->Fill(pt);
            fHistJetPt[fInCentralitySelection]->Fill(pt-area*rho);
            if(fFillQAHistograms) fHistJetEtaPhi[fInCentralitySelection]->Fill(eta, phi);
            fHistJetPtArea[fInCentralitySelection]->Fill(pt-area*rho, area);
            fHistJetPtEta[fInCentralitySelection]->Fill(pt-area*rho, eta);
            fHistJetPsi2Pt[fInCentralitySelection]->Fill(PhaseShift(phi-psi2, 2.), pt-area*rho);
            fHistJetPsi2PtRho0[fInCentralitySelection]->Fill(PhaseShift(phi-psi2, 2.), pt-area*fLocalRho->GetVal());
            fHistJetPtConstituents[fInCentralitySelection]->Fill(pt-area*rho, jet->Nch());
            fHistJetEtaRho[fInCentralitySelection]->Fill(eta, pt/area);
        } 
    }
}
//_____________________________________________________________________________
void AliAnalysisTaskJetV2::FillQAHistograms(AliVTrack* vtrack) const
{
    // fill qa histograms for pico tracks
    if(fDebug > 1) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    if(!vtrack) return;
    AliPicoTrack* track = static_cast<AliPicoTrack*>(vtrack);
    fHistRunnumbersPhi->Fill(fMappedRunNumber, track->Phi());
    fHistRunnumbersEta->Fill(fMappedRunNumber, track->Eta());
    Int_t type((int)(track->GetTrackType()));
    switch (type) {
        case 0:
           fHistPicoCat1[fInCentralitySelection]->Fill(track->Eta(), track->Phi()); 
           break;
        case 1:
           fHistPicoCat2[fInCentralitySelection]->Fill(track->Eta(), track->Phi()); 
           break;
        case 2:
           fHistPicoCat3[fInCentralitySelection]->Fill(track->Eta(), track->Phi()); 
           break;
        default: break;
    }
}
//_____________________________________________________________________________
void AliAnalysisTaskJetV2::FillQAHistograms(AliVEvent* vevent) 
{
    // fill qa histograms for events
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    if(!vevent) return;
    fHistVertexz->Fill(vevent->GetPrimaryVertex()->GetZ());
    fHistCentrality->Fill(fCent);
    Int_t runNumber(InputEvent()->GetRunNumber());
    for(fMappedRunNumber = 0; fMappedRunNumber < fExpectedRuns->GetSize(); fMappedRunNumber++) {
        if(fExpectedRuns->At(fMappedRunNumber) == runNumber) return;
    }
    if(fDebug > 0) printf("\n > TASK %s CANNOT IDENTIFY RUN - CONFIGURATION COULD BE INCORRECT < \n", GetName());
}
//_____________________________________________________________________________
void AliAnalysisTaskJetV2::FillAnalysisSummaryHistogram() const
{
    // fill the analysis summary histrogram, saves all relevant analysis settigns
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    fHistAnalysisSummary->GetXaxis()->SetBinLabel(2, "fJetRadius");
    fHistAnalysisSummary->SetBinContent(2, GetJetContainer()->GetJetRadius());
    fHistAnalysisSummary->GetXaxis()->SetBinLabel(3, "fJetEtaMin");
    fHistAnalysisSummary->SetBinContent(3, GetJetContainer()->GetJetEtaMin());
    fHistAnalysisSummary->GetXaxis()->SetBinLabel(4, "fJetEtaMax");
    fHistAnalysisSummary->SetBinContent(4, GetJetContainer()->GetJetEtaMax());
    fHistAnalysisSummary->GetXaxis()->SetBinLabel(5, "fJetPhiMin");
    fHistAnalysisSummary->SetBinContent(5, GetJetContainer()->GetJetPhiMin());
    fHistAnalysisSummary->GetXaxis()->SetBinLabel(6, "fJetPhiMax");
    fHistAnalysisSummary->SetBinContent(6, GetJetContainer()->GetJetPhiMin());
    fHistAnalysisSummary->GetXaxis()->SetBinLabel(16, "fForceBeamType");
    fHistAnalysisSummary->SetBinContent(16, fForceBeamType);
    fHistAnalysisSummary->GetXaxis()->SetBinLabel(17, "fMinCent");
    fHistAnalysisSummary->SetBinContent(17, fMinCent);
    fHistAnalysisSummary->GetXaxis()->SetBinLabel(18, "fMaxCent");
    fHistAnalysisSummary->SetBinContent(18, fMaxCent);
    fHistAnalysisSummary->GetXaxis()->SetBinLabel(19, "fMinVz");
    fHistAnalysisSummary->SetBinContent(19, fMinVz);
    fHistAnalysisSummary->GetXaxis()->SetBinLabel(20, "fMaxVz");
    fHistAnalysisSummary->SetBinContent(20, fMaxVz);
    fHistAnalysisSummary->GetXaxis()->SetBinLabel(21, "fOffTrigger");
    fHistAnalysisSummary->SetBinContent(21, fOffTrigger);
    fHistAnalysisSummary->GetXaxis()->SetBinLabel(34, "fitModulationType");
    fHistAnalysisSummary->SetBinContent(34, (int)fFitModulationType);
    fHistAnalysisSummary->GetXaxis()->SetBinLabel(35, "runModeType");
    fHistAnalysisSummary->SetBinContent(35, (int)fRunModeType);
    fHistAnalysisSummary->GetXaxis()->SetBinLabel(36, "data type");
    fHistAnalysisSummary->SetBinContent(36, (int)fDataType);
    fHistAnalysisSummary->GetXaxis()->SetBinLabel(37, "iterator");
    fHistAnalysisSummary->SetBinContent(37, 1.);
    fHistAnalysisSummary->GetXaxis()->SetBinLabel(38, "fMinPvalue");
    fHistAnalysisSummary->SetBinContent(38, fMinPvalue);
    fHistAnalysisSummary->GetXaxis()->SetBinLabel(39, "fMaxPvalue");
    fHistAnalysisSummary->SetBinContent(39, fMaxPvalue);
    fHistAnalysisSummary->GetXaxis()->SetBinLabel(40, "fExcludeLeadingJetsFromFit");
    fHistAnalysisSummary->SetBinContent(40, fExcludeLeadingJetsFromFit);
    fHistAnalysisSummary->GetXaxis()->SetBinLabel(41, "fRebinSwapHistoOnTheFly");
    fHistAnalysisSummary->SetBinContent(41, (int)fRebinSwapHistoOnTheFly);
    fHistAnalysisSummary->GetXaxis()->SetBinLabel(42, "fUsePtWeight");
    fHistAnalysisSummary->SetBinContent(42, (int)fUsePtWeight);
    fHistAnalysisSummary->GetXaxis()->SetBinLabel(44, "fSoftTrackMinPt");
    fHistAnalysisSummary->SetBinContent(44, fSoftTrackMinPt);
    fHistAnalysisSummary->GetXaxis()->SetBinLabel(45, "fSoftTrackMaxPt");
    fHistAnalysisSummary->SetBinContent(45, fSoftTrackMaxPt);
    fHistAnalysisSummary->GetXaxis()->SetBinLabel(46, "fMaxCones");
    fHistAnalysisSummary->SetBinContent(46, fMaxCones);
    fHistAnalysisSummary->GetXaxis()->SetBinLabel(47, "used rho");
    fHistAnalysisSummary->GetXaxis()->SetBinLabel(48, "used small rho");
}
//_____________________________________________________________________________
void AliAnalysisTaskJetV2::Terminate(Option_t *)
{
    // terminate
    switch (fRunModeType) {
        case kLocal : {
        if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
        AliAnalysisTaskJetV2::Dump();
        for(Int_t i(0); i < fHistAnalysisSummary->GetXaxis()->GetNbins(); i++) printf( " > flag: %s \t content %.2f \n", fHistAnalysisSummary->GetXaxis()->GetBinLabel(1+i), fHistAnalysisSummary->GetBinContent(1+i));
        } break;
        default : break;
    }
}
//_____________________________________________________________________________
void AliAnalysisTaskJetV2::SetModulationFit(TF1* fit) 
{
    // set modulation fit
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    if (fFitModulation) delete fFitModulation;
    fFitModulation = fit; 
}
//_____________________________________________________________________________
void AliAnalysisTaskJetV2::SetUseControlFit(Bool_t c)
{
    // set control fit
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    if (fFitControl) delete fFitControl;
    if (c) {
        fFitControl = new TF1("controlFit", "pol0", 0, TMath::TwoPi());
    } else fFitControl = 0x0;
}
//_____________________________________________________________________________
TH1F* AliAnalysisTaskJetV2::GetResolutionFromOuptutFile(detectorType det, Int_t h, TArrayD* cen)
{
    // INTERFACE METHOD FOR OUTPUTFILE
    // get the detector resolution, user has ownership of the returned histogram
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    if(!fOutputList) {
        printf(" > Please add fOutputList first < \n");
        return 0x0;
    }
    TH1F* r(0x0);
    (cen) ? r = new TH1F("R", "R", cen->GetSize()-1, cen->GetArray()) : r = new TH1F("R", "R", 10, 0, 10);
    if(!cen) r->GetXaxis()->SetTitle("number of centrality bin");
    r->GetYaxis()->SetTitle(Form("Resolution #Psi_{%i}", h));
    for(Int_t i(0); i < 10; i++) {
        TProfile* temp((TProfile*)fOutputList->FindObject(Form("fProfV%iResolution_%i", h, i)));
        if(!temp) break;
        Double_t a(temp->GetBinContent(3)), b(temp->GetBinContent(5)), c(temp->GetBinContent(7));
        Double_t d(temp->GetBinContent(9)), e(temp->GetBinContent(10)), f(temp->GetBinContent(11));
        Double_t _a(temp->GetBinError(3)), _b(temp->GetBinError(5)), _c(temp->GetBinError(7));
        Double_t _d(temp->GetBinError(9)), _e(temp->GetBinError(10)), _f(temp->GetBinError(11));
        if(a <= 0 || b <= 0 || c <= 0 || d <= 0 || e <= 0 || f <= 0) continue;
        switch (det) {
            case kVZEROA : {
                r->SetBinContent(1+i, TMath::Sqrt((a*b)/c));
                if(i==0) r->SetNameTitle("VZEROA resolution", "VZEROA resolution");
                r->SetBinError(1+i, TMath::Sqrt(_a*_a+_b*_b+_c*_c));
            } break;
            case kVZEROC : {
                r->SetBinContent(1+i, TMath::Sqrt((a*c)/b));
                if(i==0) r->SetNameTitle("VZEROC resolution", "VZEROC resolution");
                r->SetBinError(1+i, TMath::Sqrt(_a*_a+_b*_b+_c*_c));
            } break;
            case kTPC : {
                r->SetBinContent(1+i, TMath::Sqrt((b*c)/a));
                if(i==0) r->SetNameTitle("TPC resolution", "TPC resolution");
                r->SetBinError(1+i, TMath::Sqrt(_a*_a+_b*_b+_c*_c));
            } break;
            case kVZEROComb : {
                r->SetBinContent(1+i, TMath::Sqrt((d*e)/f));
                if(i==0) r->SetNameTitle("VZEROComb resolution", "VZEROComb resolution");
                r->SetBinError(1+i, TMath::Sqrt(_d*_d+_e*_e+_f*_f));
            } break;
            default : break;
        }
    }
    return r;
}
//_____________________________________________________________________________
TH1F* AliAnalysisTaskJetV2::CorrectForResolutionDiff(TH1F* v, detectorType det, TArrayD* cen, Int_t c, Int_t h)
{
    // INTERFACE METHOD FOR OUTPUT FILE
    // correct the supplied differential vn histogram v for detector resolution
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    TH1F* r(GetResolutionFromOuptutFile(det, h, cen));
    if(!r) {
        printf(" > Couldn't find resolution < \n");
        return 0x0;
    }
    Double_t res(1./r->GetBinContent(1+r->FindBin(c)));
    TF1* line = new TF1("line", "pol0", 0, 200);
    line->SetParameter(0, res);
    v->Multiply(line);
    return v;
}
//_____________________________________________________________________________
TH1F* AliAnalysisTaskJetV2::CorrectForResolutionInt(TH1F* v, detectorType det, TArrayD* cen, Int_t h)
{
    // INTERFACE METHOD FOR OUTPUT FILE
    // correct the supplied intetrated vn histogram v for detector resolution
    // integrated vn must have the same centrality binning as the resolotion correction
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    TH1F* r(GetResolutionFromOuptutFile(det, h, cen));
    v->Divide(v, r);
    return v;
}
//_____________________________________________________________________________
TH1F* AliAnalysisTaskJetV2::GetDifferentialQC(TProfile* refCumulants, TProfile* diffCumlants, TArrayD* ptBins, Int_t h)
{
    // get differential QC
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    Double_t r(refCumulants->GetBinContent(h-1)); // v2 reference flow
    if(r > 0) r = TMath::Sqrt(r);
    TH1F* qc = new TH1F(Form("QC2v%i", h), Form("QC2v%i", h), ptBins->GetSize()-1, ptBins->GetArray());
    Double_t a(0), b(0), c(0);  // dummy variables
    for(Int_t i(0); i < ptBins->GetSize(); i++) {
        if(r > 0) {
            a = diffCumlants->GetBinContent(1+i);
            b = diffCumlants->GetBinError(1+i);
            c = a/r;
            qc->SetBinContent(1+i, c);
            (a <= 0 || b <= 0) ? qc->SetBinError(1+i, b) : qc->SetBinError(1+i, TMath::Sqrt(c*c*b*b/(a*a)));
        }
    }
    return qc;
}

//_____________________________________________________________________________
void AliAnalysisTaskJetV2::ReadVZEROCalibration2010h()
{
    // necessary for calibration of 10h vzero event plane. code copied from flow package 
    // (duplicate, but i didn't want to introduce an ulgy dependency )
    // this function is only called when the runnumber changes 
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);

    // 1) check if the proper chi weights for merging vzero a and vzero c ep are present
    // if not, use sane defaults. centrality binning is equal to that given in the fVZEROcentralityBin snippet
    //
    // chi values can be calculated using the static helper function 
    // AliAnalysisTaskJetV2::CalculateEventPlaneChi(Double_t res) where res is the event plane
    // resolution in a given centrality bin
    //
    // the resolutions that were used for these defaults are
    // this might need a bit of updating as they were read 'by-eye' from a performance plot ..
    // Double_t R2VZEROA[] = {.35, .40, .48, .50, .48, .45, .38, .26, .16};
    // Double_t R2VZEROC[] = {.45, .60, .70, .73, .68, .60, .40, .36, .17};
    // Double_t R3VZEROA[] = {.22, .23, .22, .19, .15, .12, .08, .00, .00};
    // Double_t R3VZEROC[] = {.30, .30, .28, .25, .22, .17, .11, .00, .00};

    Double_t chiC2[] = {0.771423, 1.10236, 1.38116, 1.48077, 1.31964, 1.10236, 0.674622, 0.600403, 0.273865};
    Double_t chiA2[] = {0.582214, 0.674622, 0.832214, 0.873962, 0.832214, 0.771423, 0.637146, 0.424255, 0.257385};
    Double_t chiC3[] = {0.493347, 0.493347, 0.458557, 0.407166, 0.356628, 0.273865, 0.176208, 6.10352e-05, 6.10352e-05};
    Double_t chiA3[] = {0.356628, 0.373474, 0.356628, 0.306702, 0.24115, 0.192322, 0.127869, 6.10352e-05, 6.10352e-05};

    if(!fChi2A) fChi2A = new TArrayD(9, chiA2);
    if(!fChi2C) fChi2C = new TArrayD(9, chiC2);
    if(!fChi3A) fChi3A = new TArrayD(9, chiA3);
    if(!fChi3C) fChi3C = new TArrayD(9, chiC3);

    // 2) open database file
    fOADB = TFile::Open("$ALICE_ROOT/OADB/PWGCF/VZERO/VZEROcalibEP.root");
    if(fOADB->IsZombie()){
	printf("OADB file $ALICE_ROOT/OADB/PWGCF/VZERO/VZEROcalibEP.root cannot be opened, CALIBRATION FAILED !");
	return;
    }

    AliOADBContainer *cont = (AliOADBContainer*) fOADB->Get("hMultV0BefCorr");
    if(!cont){
        // see if database is readable
	printf("OADB object hMultV0BefCorr is not available in the file\n");
	return;	
    }
    Int_t run(fRunNumber);
    if(!(cont->GetObject(run))){
        // if the run isn't recognized fall back to a default run
	printf("OADB object hMultVZEROBefCorr is not available for run %i (used default run 137366)\n",run);
	run = 137366;
    }
    // step 3) get the proper multiplicity weights from the vzero signal
    fVZEROgainEqualization = ((TH2F*)cont->GetObject(run))->ProfileX();
    if(!fVZEROgainEqualization) {
        AliFatal(Form("%s: Fatal error, couldn't read fVZEROgainEqualization from OADB object < \n", GetName()));
        return;
    }

    TF1* fpol0 = new TF1("fpol0","pol0"); 
    if(fVZEROgainEqualizationPerRing) {
        // do the calibration per ring
        // start with the vzero c rings (segments 0 through 31)
        fVZEROgainEqualization->Fit(fpol0, "", "", 0, 8);
        (fUseVZERORing[0]) ? SetVZEROCpol(0, fpol0->GetParameter(0)) : SetVZEROCpol(0, 0.);
        fVZEROgainEqualization->Fit(fpol0, "", "", 8, 16);
        (fUseVZERORing[1]) ? SetVZEROCpol(1, fpol0->GetParameter(0)) : SetVZEROCpol(1, 0.);
        fVZEROgainEqualization->Fit(fpol0, "", "", 16, 24);
        (fUseVZERORing[2]) ? SetVZEROCpol(2, fpol0->GetParameter(0)) : SetVZEROCpol(2, 0.);
        fVZEROgainEqualization->Fit(fpol0, "", "", 24, 32);
        (fUseVZERORing[3]) ? SetVZEROCpol(3, fpol0->GetParameter(0)) : SetVZEROCpol(3, 0.);
        // same thing for vero A
        fVZEROgainEqualization->Fit(fpol0, "", "", 32, 40);
        (fUseVZERORing[4]) ? SetVZEROApol(0, fpol0->GetParameter(0)) : SetVZEROApol(0, 0.);
        fVZEROgainEqualization->Fit(fpol0, "", "", 40, 48);
        (fUseVZERORing[5]) ? SetVZEROApol(1, fpol0->GetParameter(0)) : SetVZEROApol(1, 0.);
        fVZEROgainEqualization->Fit(fpol0, "", "", 48, 56);
        (fUseVZERORing[6]) ? SetVZEROApol(2, fpol0->GetParameter(0)) : SetVZEROApol(2, 0.);
        fVZEROgainEqualization->Fit(fpol0, "", "", 56, 64);
        (fUseVZERORing[7]) ? SetVZEROApol(3, fpol0->GetParameter(0)) : SetVZEROApol(3, 0.);
    } else {
        // do the calibration in one go. the calibration will still be 
        // stored per ring, but each ring has the same weight now
        // this should be the default for the analysis as the database is tuned to this configuration
       fVZEROgainEqualization->Fit(fpol0,"","",0,31);
       for(Int_t i(0); i < 4; i++) SetVZEROCpol(i, fpol0->GetParameter(0));
       fVZEROgainEqualization->Fit(fpol0,"","",32,64);
       for(Int_t i(0); i < 4; i++) SetVZEROApol(i, fpol0->GetParameter(0));
    }

    // step 4) extract the information to re-weight the q-vectors 
    for(Int_t iside=0;iside<2;iside++){
	for(Int_t icoord=0;icoord<2;icoord++){
	    for(Int_t i=0;i  < 9;i++){
		char namecont[100];
  		if(iside==0 && icoord==0)
		  snprintf(namecont,100,"hQxc2_%i",i);
		else if(iside==1 && icoord==0)
		  snprintf(namecont,100,"hQxa2_%i",i);
		else if(iside==0 && icoord==1)
		  snprintf(namecont,100,"hQyc2_%i",i);
		else if(iside==1 && icoord==1)
		  snprintf(namecont,100,"hQya2_%i",i);

		cont = (AliOADBContainer*) fOADB->Get(namecont);
		if(!cont){
		    printf("OADB object %s is not available in the file\n",namecont);
		    return;	
		}
	
		if(!(cont->GetObject(run))){
		    printf("OADB object %s is not available for run %i (used run 137366)\n",namecont,run);
		    run = 137366;
		}

                // store info for all centralities to cache
                fMeanQ[i][iside][icoord] = ((TH1F *) cont->GetObject(run))->GetMean();
		fWidthQ[i][iside][icoord] = ((TH1F *) cont->GetObject(run))->GetRMS();

		//for v3
		if(iside==0 && icoord==0)
		  snprintf(namecont,100,"hQxc3_%i",i);
		else if(iside==1 && icoord==0)
		  snprintf(namecont,100,"hQxa3_%i",i);
		else if(iside==0 && icoord==1)
		  snprintf(namecont,100,"hQyc3_%i",i);
		else if(iside==1 && icoord==1)
		  snprintf(namecont,100,"hQya3_%i",i);

		cont = (AliOADBContainer*) fOADB->Get(namecont);
		if(!cont){
		    printf("OADB object %s is not available in the file\n",namecont);
		    return;	
		}
		
		if(!(cont->GetObject(run))){
		    printf("OADB object %s is not available for run %i (used run 137366)\n",namecont,run);
		    run = 137366;
		}
                // store info for all centralities to cache
		fMeanQv3[i][iside][icoord] = ((TH1F *) cont->GetObject(run))->GetMean();
		fWidthQv3[i][iside][icoord] = ((TH1F *) cont->GetObject(run))->GetRMS();
     	    }
	}
    }
    // cleanup. the opened file is closed in the destructor, otherwise fVZEROgainEqualization is no longer available
    delete fpol0;
}
//_____________________________________________________________________________
Int_t AliAnalysisTaskJetV2::GetVZEROCentralityBin() const
{
    // return cache index number corresponding to the event centrality
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    Float_t v0Centr(InputEvent()->GetCentrality()->GetCentralityPercentile("V0M"));
    if(v0Centr < 5) return 0;
    else if(v0Centr < 10) return 1;
    else if(v0Centr < 20) return  2;
    else if(v0Centr < 30) return  3;
    else if(v0Centr < 40) return  4;
    else if(v0Centr < 50) return  5;
    else if(v0Centr < 60) return  6;
    else if(v0Centr < 70) return  7;
    else return 8;
}
//_____________________________________________________________________________
