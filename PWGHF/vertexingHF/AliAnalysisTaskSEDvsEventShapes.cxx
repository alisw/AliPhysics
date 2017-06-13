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

/* $Id$ */

//*************************************************************************
// Class AliAnalysisTaskSEDvsEventShapes
// AliAnalysisTaskSE for the D meson vs. Event shape analysis in different mutiplicity window
// Authors: Renu Bala, Manoj Bhanudas Jadhav
//*************************************************************************

#include <TClonesArray.h>
#include <TCanvas.h>
#include <TList.h>
#include <TMath.h>
#include <TString.h>
#include <TDatabasePDG.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <THnSparse.h>
#include <TProfile.h>
#include "AliAnalysisManager.h"
#include "AliRDHFCuts.h"
#include "AliRDHFCutsDplustoKpipi.h"
#include "AliRDHFCutsDStartoKpipi.h"
#include "AliRDHFCutsD0toKpi.h"
#include "AliAODHandler.h"
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliAODRecoDecayHF.h"
#include "AliAODRecoCascadeHF.h"
#include "AliAnalysisVertexingHF.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskSEDvsEventShapes.h"
#include "AliNormalizationCounter.h"
#include "AliVertexingHFUtils.h"
#include "AliAODVZERO.h"
#include "AliESDUtils.h"
ClassImp(AliAnalysisTaskSEDvsEventShapes)

//________________________________________________________________________
AliAnalysisTaskSEDvsEventShapes::AliAnalysisTaskSEDvsEventShapes():
AliAnalysisTaskSE(),
fOutput(0),
fListCuts(0),
fOutputCounters(0),
fListProfiles(0),
fOutputEffCorr(0),
fHistNEvents(0),
fHistNtrVsZvtx(0),
fHistNtrCorrVsZvtx(0),
fHistNtrVsnTrackEvWithCand(0),
fHistNtrVsSo(0),
fHistNtrCorrVsSo(0),
fHistNtrVsSpheri(0),
fHistNtrCorrVsSpheri(0),
fHistNtrVsNchMC(0),
fHistNtrCorrVsNchMC(0),
fHistNtrVsNchMCPrimary(0),
fHistNtrCorrVsNchMCPrimary(0),
fHistNtrVsNchMCPhysicalPrimary(0),
fHistNtrCorrVsNchMCPhysicalPrimary(0),
fHistGenPrimaryParticlesInelGt0(0),
fHistNchMCVsNchMCPrimaryVsNchMCPhysicalPrimary(0),
fHistNtrCorrPSSel(0),
fHistNtrCorrEvSel(0),
fHistNtrCorrEvWithCand(0),
fHistNtrCorrEvWithD(0),
fHistnTrackvsEtavsPhi(0),
fHistnTrackvsEtavsPhiEvWithCand(0),
fHistTrueSovsMeasSo(0),
fHistTrueSovsMeasSoEvWithCand(0),
fSparseEvtShape(0),
fSparseEvtShapewithNoPid(0),
fSparseEvtShapePrompt(0),
fSparseEvtShapeFeeddown(0),
fSparseEvtShapeRecSphero(0),
fMCAccGenPrompt(0),
fMCAccGenFeeddown(0),
fMCRecoPrompt(0),
fMCRecoFeeddown(0),
fMCRecoBothPromptFD(0),
fMCAccGenPromptSpheri(0),
fMCAccGenFeeddownSpheri(0),
fMCRecoPromptSpheri(0),
fMCRecoFeeddownSpheri(0),
fMCRecoBothPromptFDSpheri(0),
fMCAccGenPromptEvSel(0),
fMCAccGenFeeddownEvSel(0),
fUpmasslimit(1.965),
fLowmasslimit(1.765),
fNMassBins(200),
fRDCutsAnalysis(0),
fCounterC(0),
fCounterU(0),
fCounterCandidates(0),
fDoImpPar(kFALSE),
fNImpParBins(400),
fLowerImpPar(-2000.),
fHigherImpPar(2000.),
fReadMC(kFALSE),
fMCOption(0),
fisPPbData(kFALSE),
fUseBit(kTRUE),
fSubtractTrackletsFromDau(kFALSE),
fCalculateSphericity(kFALSE),
fRecomputeSpherocity(kFALSE),
fRemoveD0fromDstar(kFALSE),
fUseNchWeight(0),
fHistoMCNch(0),
fHistoMeasNch(0),
fUsePtWeight(kFALSE),
fWeight(1.),
fRefMult(9.26),
fPdgMeson(411),
fMultiplicityEstimator(kNtrk10),
fMCPrimariesEstimator(kEta10),
fDoVZER0ParamVertexCorr(1),
fFillSoSparseChecks(0),
fUseQuarkTag(kTRUE),
fFillTrackHisto(kFALSE),
fEtaAccCut(0.9),
fPtAccCut(0.1),
fetaMin(-0.8),
fetaMax(0.8),
fptMin(0.15),
fptMax(10.),
fminMult(3),
ffiltbit1(256),
ffiltbit2(512),
fphiStepSizeDeg(0.1)
{
    // Default constructor
    for(Int_t i=0; i<5; i++) fHistMassPtImpPar[i]=0;
    for(Int_t i=0; i<4; i++) fMultEstimatorAvg[i]=0;
}

//________________________________________________________________________
AliAnalysisTaskSEDvsEventShapes::AliAnalysisTaskSEDvsEventShapes(const char *name, Int_t pdgMeson,AliRDHFCuts *cuts, Bool_t switchPPb):
AliAnalysisTaskSE(name),
fOutput(0),
fListCuts(0),
fOutputCounters(0),
fListProfiles(0),
fOutputEffCorr(0),
fHistNEvents(0),
fHistNtrVsZvtx(0),
fHistNtrCorrVsZvtx(0),
fHistNtrVsnTrackEvWithCand(0),
fHistNtrVsSo(0),
fHistNtrCorrVsSo(0),
fHistNtrVsSpheri(0),
fHistNtrCorrVsSpheri(0),
fHistNtrVsNchMC(0),
fHistNtrCorrVsNchMC(0),
fHistNtrVsNchMCPrimary(0),
fHistNtrCorrVsNchMCPrimary(0),
fHistNtrVsNchMCPhysicalPrimary(0),
fHistNtrCorrVsNchMCPhysicalPrimary(0),
fHistGenPrimaryParticlesInelGt0(0),
fHistNchMCVsNchMCPrimaryVsNchMCPhysicalPrimary(0),
fHistNtrCorrPSSel(0),
fHistNtrCorrEvSel(0),
fHistNtrCorrEvWithCand(0),
fHistNtrCorrEvWithD(0),
fHistnTrackvsEtavsPhi(0),
fHistnTrackvsEtavsPhiEvWithCand(0),
fHistTrueSovsMeasSo(0),
fHistTrueSovsMeasSoEvWithCand(0),
fSparseEvtShape(0),
fSparseEvtShapewithNoPid(0),
fSparseEvtShapePrompt(0),
fSparseEvtShapeFeeddown(0),
fSparseEvtShapeRecSphero(0),
fMCAccGenPrompt(0),
fMCAccGenFeeddown(0),
fMCRecoPrompt(0),
fMCRecoFeeddown(0),
fMCRecoBothPromptFD(0),
fMCAccGenPromptSpheri(0),
fMCAccGenFeeddownSpheri(0),
fMCRecoPromptSpheri(0),
fMCRecoFeeddownSpheri(0),
fMCRecoBothPromptFDSpheri(0),
fMCAccGenPromptEvSel(0),
fMCAccGenFeeddownEvSel(0),
fUpmasslimit(1.965),
fLowmasslimit(1.765),
fNMassBins(200),
fRDCutsAnalysis(cuts),
fCounterC(0),
fCounterU(0),
fCounterCandidates(0),
fDoImpPar(kFALSE),
fNImpParBins(400),
fLowerImpPar(-2000.),
fHigherImpPar(2000.),
fReadMC(kFALSE),
fMCOption(0),
fisPPbData(switchPPb),
fUseBit(kTRUE),
fSubtractTrackletsFromDau(kFALSE),
fCalculateSphericity(kFALSE),
fRecomputeSpherocity(kFALSE),
fRemoveD0fromDstar(kFALSE),
fUseNchWeight(0),
fHistoMCNch(0),
fHistoMeasNch(0),
fUsePtWeight(kFALSE),
fWeight(1.),
fRefMult(9.26),
fPdgMeson(pdgMeson),
fMultiplicityEstimator(kNtrk10),
fMCPrimariesEstimator(kEta10),
fDoVZER0ParamVertexCorr(1),
fFillSoSparseChecks(0),
fUseQuarkTag(kTRUE),
fFillTrackHisto(kFALSE),
fEtaAccCut(0.9),
fPtAccCut(0.1),
fetaMin(-0.8),
fetaMax(0.8),
fptMin(0.15),
fptMax(10.),
fminMult(3),
ffiltbit1(256),
ffiltbit2(512),
fphiStepSizeDeg(0.1)
{
    //
    // Standard constructor
    //
    
    for(Int_t i=0; i<5; i++) fHistMassPtImpPar[i]=0;
    for(Int_t i=0; i<4; i++) fMultEstimatorAvg[i]=0;
    if(fPdgMeson==413){
        fNMassBins=200;
        SetMassLimits(0.12,0.2);
    }else{
        fNMassBins=200;
        SetMassLimits(fPdgMeson,0.1);
    }
    // Default constructor
    // Otput slot #1 writes into a TList container
    DefineOutput(1,TList::Class());  //My private output
    // Output slot #2 writes cut to private output
    DefineOutput(2,TList::Class());
    // Output slot #3 writes cut to private output
    DefineOutput(3,TList::Class());
    // Output slot #4 writes cut to private output
    DefineOutput(4,TList::Class());
    // Output slot #5 writes cut to private output
    DefineOutput(5,TList::Class());
}
//________________________________________________________________________
AliAnalysisTaskSEDvsEventShapes::~AliAnalysisTaskSEDvsEventShapes()
{
    //
    // Destructor
    //
    delete fOutput;
    delete fHistNEvents;
    delete fListCuts;
    delete fListProfiles;
    delete fOutputEffCorr;
    delete fRDCutsAnalysis;
    delete fCounterC;
    delete fCounterU;
    delete fCounterCandidates;
    
    for(Int_t i=0; i<4; i++) {
        if (fMultEstimatorAvg[i]) delete fMultEstimatorAvg[i];
    }
    
    for(Int_t i=0; i<5; i++){
        delete fHistMassPtImpPar[i];
    }
    if(fHistoMCNch) delete fHistoMCNch;
    if(fHistoMeasNch) delete fHistoMeasNch;
}

//_________________________________________________________________
void  AliAnalysisTaskSEDvsEventShapes::SetMassLimits(Double_t lowlimit, Double_t uplimit){
    // set invariant mass limits
    if(uplimit>lowlimit){
        fLowmasslimit = lowlimit;
        fUpmasslimit = uplimit;
    }else{
        AliError("Wrong mass limits: upper value should be larger than lower one");
    }
}
//_________________________________________________________________
void  AliAnalysisTaskSEDvsEventShapes::SetMassLimits(Int_t pdg, Double_t range){
    // set invariant mass limits
    Double_t mass=TDatabasePDG::Instance()->GetParticle(TMath::Abs(pdg))->Mass();
    SetMassLimits(mass-range,mass+range);
}
//________________________________________________________________________
void AliAnalysisTaskSEDvsEventShapes::Init(){
    //
    // Initialization
    //
    printf("AnalysisTaskSEDvsMultiplicity_0::Init() \n");
    
    if(fUseNchWeight && !fReadMC){ AliFatal("Nch weights can only be used in MC mode"); return; }
    if(fUsePtWeight && !fReadMC){ AliFatal("pT weights can only be used in MC mode"); return; }
    if(fUsePtWeight && fUseNchWeight) { AliInfo("Beware, using at the same time pT and Nch weights, please check"); }
    if(fUseNchWeight && !fHistoMCNch){ AliFatal("Nch weights can only be used without histogram"); return; }
    if(fUseNchWeight==1 && !fHistoMeasNch) {//Nch weights
        if(fisPPbData){ AliFatal("Nch weights can only be used with MC and data histogram in pPb"); return; }
        else CreateMeasuredNchHisto();
    }
    if(fUseNchWeight==2 && !fHistoMeasNch){ AliFatal("Ntrk weights can only be used with MC and data histogram"); return; } //for pp, pPb Ntrk weights
    
    fListCuts=new TList();
    fListCuts->SetOwner();
    fListCuts->SetName("CutsList");
    
    if(fPdgMeson==411){
        AliRDHFCutsDplustoKpipi* copycut=new AliRDHFCutsDplustoKpipi(*(static_cast<AliRDHFCutsDplustoKpipi*>(fRDCutsAnalysis)));
        copycut->SetName("AnalysisCutsDplus");
        fListCuts->Add(copycut);
    }else if(fPdgMeson==421){
        AliRDHFCutsD0toKpi* copycut=new AliRDHFCutsD0toKpi(*(static_cast<AliRDHFCutsD0toKpi*>(fRDCutsAnalysis)));
        copycut->SetName("AnalysisCutsDzero");
        fListCuts->Add(copycut);
    }else if(fPdgMeson==413){
        AliRDHFCutsDStartoKpipi* copycut=new AliRDHFCutsDStartoKpipi(*(static_cast<AliRDHFCutsDStartoKpipi*>(fRDCutsAnalysis)));
        copycut->SetName("AnalysisCutsDStar");
        fListCuts->Add(copycut);
    }
    if(fHistoMeasNch) fListCuts->Add(fHistoMeasNch);
    if(fHistoMCNch) fListCuts->Add(fHistoMCNch);
    
    PostData(2,fListCuts);
    
    fListProfiles = new TList();
    fListProfiles->SetOwner();
    TString period[4];
    Int_t nProfiles=4;
    if (fisPPbData) {period[0]="LHC13b"; period[1]="LHC13c"; nProfiles = 2;}
    else {period[0]="LHC10b"; period[1]="LHC10c"; period[2]="LHC10d"; period[3]="LHC10e"; nProfiles = 4;}
    
    for(Int_t i=0; i<nProfiles; i++){
        if(fMultEstimatorAvg[i]){
            TProfile* hprof=new TProfile(*fMultEstimatorAvg[i]);
            hprof->SetName(Form("ProfileTrkVsZvtx%s\n",period[i].Data()));
            fListProfiles->Add(hprof);
        }
    }
    
    PostData(4,fListProfiles);
    
    return;
}

//________________________________________________________________________
void AliAnalysisTaskSEDvsEventShapes::UserCreateOutputObjects()
{
    // Create the output container
    //
    if(fDebug > 1) printf("AnalysisTaskSEDvsMultiplicity::UserCreateOutputObjects() \n");
    
    // Several histograms are more conveniently managed in a TList
    fOutput = new TList();
    fOutput->SetOwner();
    fOutput->SetName("OutputHistos");
    
    fOutputEffCorr = new TList();
    fOutputEffCorr->SetOwner();
    fOutputEffCorr->SetName("OutputEffCorrHistos");
    
    Int_t nMultBins = 200;
    Float_t firstMultBin = -0.5;
    Float_t lastMultBin = 199.5;
    Int_t nMultBinsNtrk = nMultBins;
    Float_t lastMultBinNtrk = lastMultBin;
    Int_t nMultBinsV0 = 400;
    Float_t lastMultBinV0 = 799.5;
    const char *estimatorName="tracklets";
    if(fisPPbData) {
        nMultBinsNtrk = 375;
        lastMultBinNtrk = 374.5;
        nMultBins = nMultBinsNtrk;
        lastMultBin = lastMultBinNtrk;
    }
    if(fMultiplicityEstimator==kVZERO || fMultiplicityEstimator==kVZEROA || fMultiplicityEstimator==kVZEROEq || fMultiplicityEstimator==kVZEROAEq) {
        nMultBins = nMultBinsV0;
        lastMultBin = lastMultBinV0;
        estimatorName = "vzero";
    }
    
    fHistNtrCorrPSSel = new TH1F("hNtrCorrPSSel",Form("Corrected %s multiplicity for PS selected events; %s ; Entries",estimatorName,estimatorName),nMultBins,firstMultBin,lastMultBin);
    fHistNtrCorrEvSel = new TH1F("hNtrCorrEvSel",Form("Corrected %s multiplicity for selected events; %s ; Entries",estimatorName,estimatorName),nMultBins,firstMultBin,lastMultBin);
    fHistNtrCorrEvWithCand = new TH1F("hNtrCorrEvWithCand", Form("%s multiplicity for events with D candidates; %s ; Entries",estimatorName,estimatorName),nMultBins,firstMultBin,lastMultBin);// Total multiplicity
    fHistNtrCorrEvWithD = new TH1F("hNtrCorrEvWithD", Form("%s multiplicity for events with D in mass region ; %s ; Entries",estimatorName,estimatorName),nMultBins,firstMultBin,lastMultBin); //
    
    fHistNtrVsZvtx = new TH2F("hNtrVsZvtx",Form("N%s vs VtxZ; VtxZ;N_{%s};",estimatorName,estimatorName),300,-15,15,nMultBins,firstMultBin,lastMultBin); //
    fHistNtrCorrVsZvtx = new TH2F("hNtrCorrVsZvtx",Form("N%s vs VtxZ; VtxZ;N_{%s};",estimatorName,estimatorName),300,-15,15,nMultBins,firstMultBin,lastMultBin); //
    fHistNtrVsnTrackEvWithCand = new TH2F("hNtrVsnTrackEvWithCand",Form("N%s vs nTracks; nTracks; N_{%s};",estimatorName,estimatorName),nMultBins,firstMultBin,lastMultBin,nMultBins,firstMultBin,lastMultBin); //
    
    TString histoNtrName;
    TString histoNtrCorrName;
    TString parNameNtr;
    
    histoNtrName = "hNtrVsSphero";
    histoNtrCorrName = "hNtrCorrVsSphero";
    parNameNtr = "Sphero";
    
    fHistNtrVsSo = new TH2F(histoNtrName.Data(),Form("N_{%s} vs %s; %s; N_{%s};",estimatorName,parNameNtr.Data(),parNameNtr.Data(),estimatorName), 20, 0., 1., nMultBins,firstMultBin,lastMultBin); //
    fHistNtrCorrVsSo = new TH2F(histoNtrCorrName.Data(),Form("N_{%s} vs %s; %s; N_{%s};",estimatorName,parNameNtr.Data(),parNameNtr.Data(),estimatorName), 20, 0., 1., nMultBins, firstMultBin,lastMultBin); //
    
    fHistnTrackvsEtavsPhi = new TH3F("hnTrackvsEtavsPhi", "Eta vs Phi vs nTracks; #eta; #phi[rad]; nTracks;", 100, -1.5, 1.5, 100, 0., 6.28,nMultBins,firstMultBin,lastMultBin);
    fHistnTrackvsEtavsPhiEvWithCand = new TH3F("hnTrackvsEtavsPhiEvWithCand", "Eta vs Phi vs nTracks; #eta; #phi[rad]; nTracks;", 100, -1.5, 1.5, 100, 0., 6.28,nMultBins,firstMultBin,lastMultBin);
    
    fHistTrueSovsMeasSo = new TH3F("hTrueSovsMeasSo", "trueSo vs measSo; S_{o} (true); S_{o} (meas); tracklets;", 100, 0., 1., 100, 0., 1., nMultBins, firstMultBin, lastMultBin);
    fHistTrueSovsMeasSoEvWithCand = new TH3F("hTrueSovsMeasSoEvWithCand", "trueSo vs measSo; S_{o} (true); S_{o} (meas); tracklets;", 100, 0., 1., 100, 0., 1., nMultBins, firstMultBin, lastMultBin);
    
    TString histoNtrSphriName;
    TString histoNtrCorrSphriName;
    TString parNameNtrSphri;
    
    histoNtrSphriName = "hNtrVsSpheri";
    histoNtrCorrSphriName = "hNtrCorrVsSpheri";
    parNameNtrSphri = "Spheri";
    
    if(fCalculateSphericity){
        fHistNtrVsSpheri = new TH2F(histoNtrSphriName.Data(),Form("N_{%s} vs %s; %s; N_{%s};",estimatorName,parNameNtrSphri.Data(),parNameNtrSphri.Data(),estimatorName), 20, 0., 1., nMultBins,firstMultBin,lastMultBin); //
        fHistNtrCorrVsSpheri = new TH2F(histoNtrCorrSphriName.Data(),Form("N_{%s} vs %s; %s; N_{%s};",estimatorName,parNameNtrSphri.Data(),parNameNtrSphri.Data(),estimatorName), 20, 0., 1., nMultBins, firstMultBin,lastMultBin); //
    }
    
    fHistNtrVsNchMC = new TH2F("hNtrVsNchMC",Form("N%s vs NchMC; Nch;N_{%s};",estimatorName,estimatorName),nMultBins,firstMultBin,lastMultBin,nMultBins,firstMultBin,lastMultBin); //
    fHistNtrCorrVsNchMC = new TH2F("hNtrCorrVsNchMC",Form("N%s vs Nch; Nch;N_{%s};",estimatorName,estimatorName),nMultBins,firstMultBin,lastMultBin,nMultBins,firstMultBin,lastMultBin); //
    
    fHistNtrVsNchMCPrimary = new TH2F("hNtrVsNchMCPrimary",Form("N%s vs Nch (Primary); Nch (Primary);N_{%s};",estimatorName,estimatorName),nMultBins,firstMultBin,lastMultBin,nMultBins,firstMultBin,lastMultBin); //
    fHistNtrCorrVsNchMCPrimary = new TH2F("hNtrCorrVsNchMCPrimary",Form("N%s vs Nch (Primary); Nch(Primary) ;N_{%s};",estimatorName,estimatorName),nMultBins,firstMultBin,lastMultBin,nMultBins,firstMultBin,lastMultBin); //
    
    fHistNtrVsNchMCPhysicalPrimary = new TH2F("hNtrVsNchMCPhysicalPrimary",Form("N%s vs Nch (Physical Primary); Nch (Physical Primary);N_{%s};",estimatorName,estimatorName),nMultBins,firstMultBin,lastMultBin,nMultBins,firstMultBin,lastMultBin); //
    fHistNtrCorrVsNchMCPhysicalPrimary = new TH2F("hNtrCorrVsMCPhysicalPrimary",Form("N%s vs Nch (Physical Primary); Nch (Physical Primary);N_{%s};",estimatorName,estimatorName),nMultBins,firstMultBin,lastMultBin,nMultBins,firstMultBin,lastMultBin); //
    
    fHistGenPrimaryParticlesInelGt0 = new TH1F("hGenPrimaryParticlesInelGt0","Multiplcity of generated charged particles ; Nparticles ; Entries",nMultBins,firstMultBin,lastMultBin);
    
    fHistNchMCVsNchMCPrimaryVsNchMCPhysicalPrimary = new TH3F("fHistNchMCVsNchMCPrimaryVsNchMCPhysicalPrimary", "MC: Nch (Physical Primary) vs Nch (Primary) vs Nch (Generated); Nch (Generated); Nch (Primary); Nch (Physical Primary)",nMultBins,firstMultBin,lastMultBin,nMultBins,firstMultBin,lastMultBin,nMultBins,firstMultBin,lastMultBin);
    
    fOutput->Add(fHistNtrCorrPSSel);
    fOutput->Add(fHistNtrCorrEvSel);
    fOutput->Add(fHistNtrCorrEvWithCand);
    fOutput->Add(fHistNtrCorrEvWithD);
    
    fOutput->Add(fHistNtrVsZvtx);
    fOutput->Add(fHistNtrCorrVsZvtx);
    fOutput->Add(fHistNtrVsnTrackEvWithCand);
    fOutput->Add(fHistnTrackvsEtavsPhi);
    fOutput->Add(fHistnTrackvsEtavsPhiEvWithCand);
    fOutput->Add(fHistTrueSovsMeasSo);
    fOutput->Add(fHistTrueSovsMeasSoEvWithCand);
    
    fOutput->Add(fHistNtrVsSo);
    fOutput->Add(fHistNtrCorrVsSo);
    if(fCalculateSphericity){
        fOutput->Add(fHistNtrVsSpheri);
        fOutput->Add(fHistNtrCorrVsSpheri);
    }
    fOutput->Add(fHistNtrVsNchMC);
    fOutput->Add(fHistNtrCorrVsNchMC);
    fOutput->Add(fHistNtrVsNchMCPrimary);
    fOutput->Add(fHistNtrCorrVsNchMCPrimary);
    fOutput->Add(fHistNtrVsNchMCPhysicalPrimary);
    fOutput->Add(fHistNtrCorrVsNchMCPhysicalPrimary);
    fOutput->Add(fHistGenPrimaryParticlesInelGt0);
    fOutput->Add(fHistNchMCVsNchMCPrimaryVsNchMCPhysicalPrimary);
    
    fHistNEvents = new TH1F("fHistNEvents", "number of events ",11,-0.5,10.5);
    fHistNEvents->GetXaxis()->SetBinLabel(1,"nEvents total");
    fHistNEvents->GetXaxis()->SetBinLabel(2,"nEvents with Z vertex");
    fHistNEvents->GetXaxis()->SetBinLabel(3,"nEvents selected");
    fHistNEvents->GetXaxis()->SetBinLabel(4,"Rejected due to trigger");
    fHistNEvents->GetXaxis()->SetBinLabel(5,"Rejected due to phys sel");
    fHistNEvents->GetXaxis()->SetBinLabel(6,"Rejected due to vertex cuts");
    fHistNEvents->GetXaxis()->SetBinLabel(7,"Rejected due to pileup");
    fHistNEvents->GetXaxis()->SetBinLabel(8,"Total no. of candidate");
    fHistNEvents->GetXaxis()->SetBinLabel(9,"no. of cand wo bitmask");
    fHistNEvents->GetXaxis()->SetBinLabel(10,"D after cuts (No PID)");
    fHistNEvents->GetXaxis()->SetBinLabel(11,"D after cuts + PID)");
    fHistNEvents->GetXaxis()->SetNdivisions(1,kFALSE);
    fHistNEvents->SetMinimum(0);
    fOutput->Add(fHistNEvents);
    
    // With flag fFillSoSparseChecks to fill THnSparse with MultUncorr and NoPid cases ( 0 = only Mult, 1 = Mult and multUncorr, 2 = NoPid and 3 is All)
    Int_t nbinsSo[4]={48, fNMassBins, 20, nMultBins};
    Double_t xminSo[4]={0., fLowmasslimit,0., firstMultBin};
    Double_t xmaxSo[4]={24., fUpmasslimit, 1., lastMultBin};
    
    Int_t nbinsSoSpheri[5]={48, fNMassBins, 20, nMultBins, 20};
    Double_t xminSoSpheri[5]={0., fLowmasslimit,0., firstMultBin, 0.};
    Double_t xmaxSoSpheri[5]={24., fUpmasslimit, 1., lastMultBin, 1.};
    
    Int_t nbinsSowithMultUncorr[5]={48, fNMassBins, 20, nMultBins, nMultBins};
    Double_t xminSowithMultUncorr[5]={0., fLowmasslimit,0., firstMultBin, firstMultBin};
    Double_t xmaxSowithMultUncorr[5]={24., fUpmasslimit, 1., lastMultBin, lastMultBin};
    
    Int_t nbinsSoSpheriwithMultUncorr[6]={48, fNMassBins, 20, nMultBins, nMultBins, 20};
    Double_t xminSoSpheriwithMultUncorr[6]={0., fLowmasslimit,0., firstMultBin, firstMultBin, 0.};
    Double_t xmaxSoSpheriwithMultUncorr[6]={24., fUpmasslimit, 1., lastMultBin, lastMultBin, 1.};
    
    TString histoName = "hSparseEvtShape";
    TString histoNameNoPid = "hSparseEvtShapewithNoPid";
    TString parNameSo = "Spherocity";
    TString parNameSpheri = "Sphericity";
    TString histoNamePrompt = "hSparseEvtShapePrompt";
    TString histoNameFeeddown = "hSparseEvtShapeFeeddown";
    TString histoNameRecSphero = "hSparseEvtShapeRecSphero";
    
    if(fFillSoSparseChecks == 1 || fFillSoSparseChecks == 3){
        if(fCalculateSphericity) fSparseEvtShape = new THnSparseD(histoName.Data(), Form("D candidates:; p_{T} [GeV/c]; InvMass [GeV/c^{2}]; %s; Multipicity; MultipicityUncorr; %s;", parNameSo.Data(), parNameSpheri.Data()), 6 , nbinsSoSpheriwithMultUncorr, xminSoSpheriwithMultUncorr, xmaxSoSpheriwithMultUncorr);
        else fSparseEvtShape = new THnSparseD(histoName.Data(), Form("D candidates:; p_{T} [GeV/c]; InvMass [GeV/c^{2}]; %s; Multipicity; MultipicityUncorr;", parNameSo.Data()), 5 , nbinsSowithMultUncorr, xminSowithMultUncorr, xmaxSowithMultUncorr);
    }
    else{
        if(fCalculateSphericity) fSparseEvtShape = new THnSparseD(histoName.Data(), Form("D candidates:; p_{T} [GeV/c]; InvMass [GeV/c^{2}]; %s; Multipicity; %s;", parNameSo.Data(), parNameSpheri.Data()), 5 , nbinsSoSpheri, xminSoSpheri, xmaxSoSpheri);
        else fSparseEvtShape = new THnSparseD(histoName.Data(), Form("D candidates:; p_{T} [GeV/c]; InvMass [GeV/c^{2}]; %s; Multipicity;", parNameSo.Data()), 4 , nbinsSo, xminSo, xmaxSo);
    }
    
    if(fFillSoSparseChecks == 2|| fFillSoSparseChecks == 3) {
        if(fCalculateSphericity) fSparseEvtShapewithNoPid = new THnSparseD(histoNameNoPid.Data(), Form("D candidates with NoPID:; p_{T} [GeV/c]; InvMass [GeV/c^{2}]; %s; Multipicity; %s;", parNameSo.Data(), parNameSpheri.Data()), 5 , nbinsSoSpheri, xminSoSpheri, xmaxSoSpheri);
        else  fSparseEvtShapewithNoPid = new THnSparseD(histoNameNoPid.Data(), Form("D candidates with NoPID:; p_{T} [GeV/c]; InvMass [GeV/c^{2}]; %s; Multipicity;", parNameSo.Data()), 4 , nbinsSo, xminSo, xmaxSo);
    }
    
    if(fRecomputeSpherocity) fSparseEvtShapeRecSphero = new THnSparseD(histoNameRecSphero.Data(), Form("D candidates:; p_{T} [GeV/c]; InvMass [GeV/c^{2}]; %s; Multipicity; RecSpherocity;", parNameSo.Data()), 5 , nbinsSoSpheri, xminSoSpheri, xmaxSoSpheri);
    
    fOutput->Add(fSparseEvtShape);
    if(fFillSoSparseChecks == 2 || fFillSoSparseChecks == 3) fOutput->Add(fSparseEvtShapewithNoPid);
    if(fRecomputeSpherocity) fOutput->Add(fSparseEvtShapeRecSphero);
    
    Int_t nbinsPrompt[4]={48, nMultBins, 20, 100};
    Int_t nbinsFeeddown[4]={48, nMultBins, 20, 100};
    Double_t xminPrompt[4] = {0.,firstMultBin, 0., -1.};
    Double_t xmaxPrompt[4] = {24.,lastMultBin, 1., 1.};
    Double_t xminFeeddown[4] = {0.,firstMultBin, 0., -1.};
    Double_t xmaxFeeddown[4] = {24.,lastMultBin, 1., 1.};
    
    Int_t nbinsRecSpheroPrompt[5]={48, nMultBins, 20, 100, 20};
    Int_t nbinsRecSpheroFeeddown[5]={48, nMultBins, 20, 100, 20};
    Double_t xminRecSpheroPrompt[5] = {0.,firstMultBin, 0., -1., 0.};
    Double_t xmaxRecSpheroPrompt[5] = {24.,lastMultBin, 1., 1., 1.};
    Double_t xminRecSpheroFeeddown[5] = {0.,firstMultBin, 0., -1., 0.};
    Double_t xmaxRecSpheroFeeddown[5] = {24.,lastMultBin, 1., 1., 1.};
    
    if(fReadMC){
        if(fRecomputeSpherocity){
            fSparseEvtShapePrompt = new THnSparseD(histoNamePrompt.Data(), Form("D candidates:; p_{T} [GeV/c]; InvMass [GeV/c^{2}]; %s; Multipicity; RecSpherocity;", parNameSo.Data()), 5 , nbinsSoSpheri, xminSoSpheri, xmaxSoSpheri);
            fSparseEvtShapeFeeddown = new THnSparseD(histoNameFeeddown.Data(), Form("D candidates:; p_{T} [GeV/c]; InvMass [GeV/c^{2}]; %s; Multipicity; RecSpherocity;", parNameSo.Data()), 5 , nbinsSoSpheri, xminSoSpheri, xmaxSoSpheri);
        }
        else{
            fSparseEvtShapePrompt = new THnSparseD(histoNamePrompt.Data(), Form("D candidates:; p_{T} [GeV/c]; InvMass [GeV/c^{2}]; %s; Multipicity;", parNameSo.Data()), 4 , nbinsSo, xminSo, xmaxSo);
            fSparseEvtShapeFeeddown = new THnSparseD(histoNameFeeddown.Data(), Form("D candidates:; p_{T} [GeV/c]; InvMass [GeV/c^{2}]; %s; Multipicity;", parNameSo.Data()), 4 , nbinsSo, xminSo, xmaxSo);
        }
        
        fOutput->Add(fSparseEvtShapePrompt);
        fOutput->Add(fSparseEvtShapeFeeddown);
        
        //Prompt
        if(fRecomputeSpherocity){
            fMCAccGenPrompt = new THnSparseD("hMCAccGenPrompt", "kStepMCAcceptance:; p_{T} [GeV/c]; Multipicity; Spherocity; y; RecSpherocity; - promptD",5,nbinsRecSpheroPrompt,xminRecSpheroPrompt,xmaxRecSpheroPrompt);
            fMCRecoPrompt = new THnSparseD("hMCRecoPrompt","kStepRecoPID:; p_{T} [GeV/c]; Multipicity; Spherocity; y; RecSpherocity; - promptD",5,nbinsRecSpheroPrompt,xminRecSpheroPrompt,xmaxRecSpheroPrompt);
        }
        else{
            fMCAccGenPrompt = new THnSparseD("hMCAccGenPrompt","kStepMCAcceptance:; p_{T} [GeV/c]; Multipicity; Spherocity; y; - promptD",4,nbinsPrompt,xminPrompt,xmaxPrompt);
            fMCRecoPrompt = new THnSparseD("hMCRecoPrompt","kStepRecoPID:; p_{T} [GeV/c]; Multipicity; Spherocity; y; - promptD",4,nbinsPrompt,xminPrompt,xmaxPrompt);
        }
        if(fCalculateSphericity){
            fMCAccGenPromptSpheri = new THnSparseD("hMCAccGenPromptSpheri","kStepMCAcceptance:; p_{T} [GeV/c]; Multipicity; Sphericity; y; - promptD",4,nbinsPrompt,xminPrompt,xmaxPrompt);
            fMCRecoPromptSpheri = new THnSparseD("hMCRecoPromptSpheri","kStepRecoPID:; p_{T} [GeV/c]; Multipicity; Sphericity; y; - promptD",4,nbinsPrompt,xminPrompt,xmaxPrompt);
        }
        fMCAccGenPromptEvSel = new THnSparseD("hMCAccGenPromptEvSel","kStepMCAcceptanceEvSel:; p_{T} [GeV/c]; Multipicity; Spherocity; y; - promptD",4,nbinsPrompt,xminPrompt,xmaxPrompt);
        
        //Feeddown
        if(fRecomputeSpherocity){
            fMCAccGenFeeddown = new THnSparseD("hMCAccGenBFeeddown","kStepMCAcceptance:; p_{T} [GeV/c]; Multipicity; Spherocity; y; RecSpherocity; - DfromB",5,nbinsRecSpheroFeeddown,xminRecSpheroFeeddown,xmaxRecSpheroFeeddown);
            fMCRecoFeeddown = new THnSparseD("hMCRecoFeeddown","kStepRecoPID:; p_{T} [GeV/c]; Multipicity; Spherocity; y; RecSpherocity; - DfromB",5,nbinsRecSpheroFeeddown,xminRecSpheroFeeddown,xmaxRecSpheroFeeddown);
        }
        else{
            fMCAccGenFeeddown = new THnSparseD("hMCAccGenBFeeddown","kStepMCAcceptance:; p_{T} [GeV/c]; Multipicity; Spherocity; y; - DfromB",4,nbinsFeeddown,xminFeeddown,xmaxFeeddown);
            fMCRecoFeeddown = new THnSparseD("hMCRecoFeeddown","kStepRecoPID:; p_{T} [GeV/c]; Multipicity; Spherocity; y; - DfromB",4,nbinsFeeddown,xminFeeddown,xmaxFeeddown);
        }
        if(fCalculateSphericity){
            fMCAccGenFeeddownSpheri = new THnSparseD("hMCAccGenBFeeddownSpheri","kStepMCAcceptance:; p_{T} [GeV/c]; Multipicity; Sphericity; y; - DfromB",4,nbinsFeeddown,xminFeeddown,xmaxFeeddown);
            fMCRecoFeeddownSpheri = new THnSparseD("hMCRecoFeeddownSpheri","kStepRecoPID:; p_{T} [GeV/c]; Multipicity; Sphericity; y; - DfromB",4,nbinsFeeddown,xminFeeddown,xmaxFeeddown);
        }
        fMCAccGenFeeddownEvSel = new THnSparseD("hMCAccGenBFeeddownEvSel","kStepMCAcceptance:; p_{T} [GeV/c]; Multipicity; Spherocity; y; - DfromB",4,nbinsFeeddown,xminFeeddown,xmaxFeeddown);
        
        //BothPromptFeeddown
        fMCRecoBothPromptFD = new THnSparseD("hMCRecoBothPromptFD","kStepRecoPID:; p_{T} [GeV/c]; Multipicity; Spherocity; y; - BothPromptFD",4,nbinsPrompt,xminPrompt,xmaxPrompt);
        
        if(fCalculateSphericity) fMCRecoBothPromptFDSpheri = new THnSparseD("hMCRecoBothPromptFDSpheri","kStepRecoPID:; p_{T} [GeV/c]; Multipicity; Sphericity; y; - BothPromptFD",4,nbinsPrompt,xminPrompt,xmaxPrompt);
        
        fOutputEffCorr->Add(fMCAccGenPrompt);
        fOutputEffCorr->Add(fMCAccGenFeeddown);
        fOutputEffCorr->Add(fMCRecoPrompt);
        fOutputEffCorr->Add(fMCRecoFeeddown);
        fOutputEffCorr->Add(fMCRecoBothPromptFD);
        if(fCalculateSphericity){
            fOutputEffCorr->Add(fMCAccGenPromptSpheri);
            fOutputEffCorr->Add(fMCAccGenFeeddownSpheri);
            fOutputEffCorr->Add(fMCRecoPromptSpheri);
            fOutputEffCorr->Add(fMCRecoFeeddownSpheri);
            fOutputEffCorr->Add(fMCRecoBothPromptFDSpheri);
        }
        fOutputEffCorr->Add(fMCAccGenPromptEvSel);
        fOutputEffCorr->Add(fMCAccGenFeeddownEvSel);
    }
    
    if(fDoImpPar) CreateImpactParameterHistos();
    
    fCounterC = new AliNormalizationCounter("NormCounterCorrMult");
    fCounterC->SetStudyMultiplicity(kTRUE,1.);
    fCounterC->SetStudySpherocity(kTRUE,20.);
    fCounterC->Init();
    
    fCounterU = new AliNormalizationCounter("NormCounterUnCorrMult");
    fCounterU->SetStudyMultiplicity(kTRUE,1.);
    fCounterU->SetStudySpherocity(kTRUE,20.);
    fCounterU->Init();
    
    fCounterCandidates = new AliNormalizationCounter("NormCounterCorrMultCandidates");
    fCounterCandidates->SetStudyMultiplicity(kTRUE,1.);
    fCounterCandidates->SetStudySpherocity(kTRUE,20.);
    fCounterCandidates->Init();
    
    fOutputCounters = new TList();
    fOutputCounters->SetOwner();
    fOutputCounters->SetName("OutputCounters");
    fOutputCounters->Add(fCounterC);
    fOutputCounters->Add(fCounterU);
    fOutputCounters->Add(fCounterCandidates);
    
    PostData(1,fOutput);
    PostData(2,fListCuts);
    PostData(3,fOutputCounters);
    PostData(4,fListProfiles);
    if(fReadMC) PostData(5,fOutputEffCorr);
    
    return;
}

//________________________________________________________________________
void AliAnalysisTaskSEDvsEventShapes::UserExec(Option_t */*option*/)
{
    // Execute analysis for current event:
    // heavy flavor candidates association to MC truth
    
    AliAODEvent *aod = dynamic_cast<AliAODEvent*> (InputEvent());
    
    TClonesArray *arrayCand = 0;
    TString arrayName="";
    UInt_t pdgDau[3];
    Int_t nDau=0;
    Int_t selbit=0;
    if(fPdgMeson==411){
        arrayName="Charm3Prong";
        pdgDau[0]=211; pdgDau[1]=321; pdgDau[2]=211;
        nDau=3;
        selbit=AliRDHFCuts::kDplusCuts;
    }else if(fPdgMeson==421){
        arrayName="D0toKpi";
        pdgDau[0]=211; pdgDau[1]=321; pdgDau[2]=0;
        nDau=2;
        selbit=AliRDHFCuts::kD0toKpiCuts;
    }else if(fPdgMeson==413){
        arrayName="Dstar";
        pdgDau[0]=321; pdgDau[1]=211; pdgDau[2]=0; // Quoting here D0 daughters (D* ones on another variable later)
        nDau=2;
        selbit=AliRDHFCuts::kDstarCuts;
    }
    
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
            arrayCand=(TClonesArray*)aodFromExt->GetList()->FindObject(arrayName.Data());
        }
    } else if(aod) {
        arrayCand=(TClonesArray*)aod->GetList()->FindObject(arrayName.Data());
    }
    
    if(!aod || !arrayCand) {
        printf("AliAnalysisTaskSEDvsEventShapes::UserExec: Charm3Prong branch not found!\n");
        return;
    }
    
    if(fisPPbData && fReadMC){
        Int_t runnumber = aod->GetRunNumber();
        if(aod->GetTriggerMask()==0 &&
           (runnumber>=195344 && runnumber<=195677)){
            AliDebug(3,"Event rejected because of null trigger mask");
            return;
        }
    }
    
    // fix for temporary bug in ESDfilter
    // the AODs with null vertex pointer didn't pass the PhysSel
    if(!aod->GetPrimaryVertex()||TMath::Abs(aod->GetMagneticField())<0.001) return;
    
    Int_t countTreta1=0, countTreta03=0, countTreta05=0, countTreta16=0;
    AliAODTracklets* tracklets=aod->GetTracklets();
    Int_t nTr=tracklets->GetNumberOfTracklets();
    for(Int_t iTr=0; iTr<nTr; iTr++){
        Double_t theta=tracklets->GetTheta(iTr);
        Double_t eta=-TMath::Log(TMath::Tan(theta/2.));
        if(eta>-0.3 && eta<0.3) countTreta03++;
        if(eta>-0.5 && eta<0.5) countTreta05++;
        if(eta>-1.0 && eta<1.0) countTreta1++;
        if(eta>-1.6 && eta<1.6) countTreta16++;
    }
    
    Int_t vzeroMult=0, vzeroMultA=0, vzeroMultC=0;
    Int_t vzeroMultEq=0, vzeroMultAEq=0, vzeroMultCEq=0;
    AliAODVZERO *vzeroAOD = (AliAODVZERO*)aod->GetVZEROData();
    if(vzeroAOD) {
        vzeroMultA = static_cast<Int_t>(vzeroAOD->GetMTotV0A());
        vzeroMultC = static_cast<Int_t>(vzeroAOD->GetMTotV0C());
        vzeroMult = vzeroMultA + vzeroMultC;
        vzeroMultAEq = static_cast<Int_t>(AliVertexingHFUtils::GetVZEROAEqualizedMultiplicity(aod));
        vzeroMultCEq = static_cast<Int_t>(AliVertexingHFUtils::GetVZEROCEqualizedMultiplicity(aod));
        vzeroMultEq = vzeroMultAEq + vzeroMultCEq;
    }
    
    Int_t countMult = countTreta1;
    if(fMultiplicityEstimator==kNtrk03) { countMult = countTreta03; }
    else if(fMultiplicityEstimator==kNtrk05) { countMult = countTreta05; }
    else if(fMultiplicityEstimator==kNtrk10to16) { countMult = countTreta16 - countTreta1; }
    else if(fMultiplicityEstimator==kVZERO) { countMult = vzeroMult; }
    else if(fMultiplicityEstimator==kVZEROA) { countMult = vzeroMultA; }
    else if(fMultiplicityEstimator==kVZEROEq) { countMult = vzeroMultEq; }
    else if(fMultiplicityEstimator==kVZEROAEq) { countMult = vzeroMultAEq; }
    
    Double_t spherocity = -0.5;
    Double_t sphericity = -0.5;
    if(fCalculateSphericity){ //When kTRUE, it calculates Sphericity and THnSparse filled for sphericity
        sphericity=AliVertexingHFUtils::GetSphericity(aod, fetaMin, fetaMax, fptMin, fptMax, ffiltbit1, ffiltbit2, fminMult);
    }
    spherocity=AliVertexingHFUtils::GetSpherocity(aod, fetaMin, fetaMax, fptMin, fptMax, ffiltbit1, ffiltbit2, fminMult, fphiStepSizeDeg);
    
    Double_t St=1;
    fCounterU->StoreEvent(aod,fRDCutsAnalysis,fReadMC,countMult,spherocity);
    fHistNEvents->Fill(0); // count event
    
    Double_t countTreta1corr=countTreta1;
    Double_t countCorr=countMult;
    AliAODVertex *vtx1 = (AliAODVertex*)aod->GetPrimaryVertex();
    // In case of VZERO multiplicity, consider the zvtx correction flag
    //  fDoVZER0ParamVertexCorr: 0= none, 1= usual d2h, 2=AliESDUtils
    Bool_t isDataDrivenZvtxCorr=kTRUE;
    Bool_t isVtxOk=kFALSE;
    Int_t vzeroMultACorr=vzeroMultA, vzeroMultCCorr=vzeroMultC, vzeroMultCorr=vzeroMult;
    Int_t vzeroMultAEqCorr=vzeroMultAEq, vzeroMultCEqCorr=vzeroMultCEq, vzeroMultEqCorr=vzeroMultEq;
    if(vtx1){
        if(vtx1->GetNContributors()>0){
            fHistNEvents->Fill(1);
            isVtxOk=kTRUE;
        }
    }
    if(isVtxOk){
        if( (fMultiplicityEstimator==kVZERO) || (fMultiplicityEstimator==kVZEROA) ||
           (fMultiplicityEstimator==kVZEROEq) || (fMultiplicityEstimator==kVZEROAEq) ){
            if(fDoVZER0ParamVertexCorr==0){
                // do not correct
                isDataDrivenZvtxCorr=kFALSE;
            } else if (fDoVZER0ParamVertexCorr==2){
                // use AliESDUtils correction
                Float_t zvtx = vtx1->GetZ();
                isDataDrivenZvtxCorr=kFALSE;
                vzeroMultACorr = static_cast<Int_t>(AliESDUtils::GetCorrV0A(vzeroMultA,zvtx));
                vzeroMultCCorr = static_cast<Int_t>(AliESDUtils::GetCorrV0C(vzeroMultC,zvtx));
                vzeroMultCorr = vzeroMultACorr + vzeroMultCCorr;
                vzeroMultAEqCorr = static_cast<Int_t>(AliESDUtils::GetCorrV0A(vzeroMultAEq,zvtx));
                vzeroMultCEqCorr =static_cast<Int_t>( AliESDUtils::GetCorrV0C(vzeroMultCEq,zvtx));
                vzeroMultEqCorr = vzeroMultAEqCorr + vzeroMultCEqCorr;
                if(fMultiplicityEstimator==kVZERO) { countCorr = vzeroMultCorr; }
                else if(fMultiplicityEstimator==kVZEROA) { countCorr = vzeroMultACorr; }
                else if(fMultiplicityEstimator==kVZEROEq) { countCorr = vzeroMultEqCorr; }
                else if(fMultiplicityEstimator==kVZEROAEq) { countCorr = vzeroMultAEqCorr; }
            }
        }
    }
    // Data driven multiplicity z-vertex correction
    if(isVtxOk && isDataDrivenZvtxCorr){
        TProfile* estimatorAvg = GetEstimatorHistogram(aod);
        if(estimatorAvg){
            countTreta1corr=static_cast<Int_t>(AliVertexingHFUtils::GetCorrectedNtracklets(estimatorAvg,countTreta1,vtx1->GetZ(),fRefMult));
            countCorr=static_cast<Int_t>(AliVertexingHFUtils::GetCorrectedNtracklets(estimatorAvg,countMult,vtx1->GetZ(),fRefMult));
        }
    }
    
    fCounterC->StoreEvent(aod,fRDCutsAnalysis,fReadMC,countCorr,spherocity);
    
    Bool_t isEvSel=fRDCutsAnalysis->IsEventSelected(aod);
    
    if(fRDCutsAnalysis->GetWhyRejection()==5) fHistNEvents->Fill(3);
    if(fRDCutsAnalysis->GetWhyRejection()==7) fHistNEvents->Fill(4);
    if(fRDCutsAnalysis->GetWhyRejection()==6) fHistNEvents->Fill(5);
    if(fRDCutsAnalysis->GetWhyRejection()==1) fHistNEvents->Fill(6);
    
    Bool_t isEvPSRejected = fRDCutsAnalysis->IsEventRejectedDuePhysicsSelection();
    
    if(!isEvPSRejected){
        fHistNtrCorrPSSel->Fill(countCorr);
    }
    
    TClonesArray *arrayMC=0;
    AliAODMCHeader *mcHeader=0;
    
    fWeight=1.;
    Double_t nchWeight=1.0;
    
    Int_t nChargedMCEta10=0, nChargedMCEta03=0, nChargedMCEta05=0, nChargedMCEta16=0, nChargedMCEtam37tm17=0, nChargedMCEta28t51=0;
    Int_t nChargedMCPrimaryEta10=0, nChargedMCPrimaryEta03=0, nChargedMCPrimaryEta05=0, nChargedMCPrimaryEta16=0, nChargedMCPrimaryEtam37tm17=0, nChargedMCPrimaryEta28t51=0;
    Int_t nChargedMCPhysicalPrimaryEta10=0, nChargedMCPhysicalPrimaryEta03=0, nChargedMCPhysicalPrimaryEta05=0, nChargedMCPhysicalPrimaryEta16=0, nChargedMCPhysicalPrimaryEtam37tm17=0, nChargedMCPhysicalPrimaryEta28t51=0;
    Int_t nChargedMC=0, nChargedMCPrimary=0, nChargedMCPhysicalPrimary=0;
    
    Double_t genspherocity = -0.5;
    // load MC particles and get weight on Nch
    if(fReadMC){
        
        arrayMC =  (TClonesArray*)aod->GetList()->FindObject(AliAODMCParticle::StdBranchName());
        if(!arrayMC) {
            printf("AliAnalysisTaskSEDvsEventShapes::UserExec: MC particles branch not found!\n");
            return;
        }
        // load MC header
        mcHeader =  (AliAODMCHeader*)aod->GetList()->FindObject(AliAODMCHeader::StdBranchName());
        if(!mcHeader) {
            printf("AliAnalysisTaskSEDvsEventShapes::UserExec: MC header branch not found!\n");
            return;
        }
        
        for(Int_t i=0; i<arrayMC->GetEntriesFast(); i++){
            AliAODMCParticle *part=(AliAODMCParticle*)arrayMC->UncheckedAt(i);
            Int_t charge = part->Charge();
            Double_t eta = part->Eta();
            Bool_t isPrim = part->IsPrimary();
            Bool_t isPhysPrim = part->IsPhysicalPrimary();
            if(charge!=0) {
                if(eta>-0.3 && eta< 0.3) {
                    nChargedMCEta03++;
                    if(isPrim) nChargedMCPrimaryEta03++;
                    if(isPhysPrim) nChargedMCPhysicalPrimaryEta03++;
                }
                if(eta>-0.5 && eta< 0.5) {
                    nChargedMCEta05++;
                    if(isPrim) nChargedMCPrimaryEta05++;
                    if(isPhysPrim) nChargedMCPhysicalPrimaryEta05++;
                }
                if(eta>-1.0 && eta< 1.0) {
                    nChargedMCEta10++;
                    if(isPrim) nChargedMCPrimaryEta10++;
                    if(isPhysPrim) nChargedMCPhysicalPrimaryEta10++;
                }
                if(eta>-1.6 && eta< 1.6) {
                    nChargedMCEta16++;
                    if(isPrim) nChargedMCPrimaryEta16++;
                    if(isPhysPrim) nChargedMCPhysicalPrimaryEta16++;
                }
                if(eta>-3.7 && eta<-1.7) {
                    nChargedMCEtam37tm17++;
                    if(isPrim) nChargedMCPrimaryEtam37tm17++;
                    if(isPhysPrim) nChargedMCPhysicalPrimaryEtam37tm17++;
                }
                if(eta> 2.8 && eta< 5.1) {
                    nChargedMCEta28t51++;
                    if(isPrim) nChargedMCPrimaryEta28t51++;
                    if(isPhysPrim) nChargedMCPhysicalPrimaryEta28t51++;
                }
            }
        }
        
        nChargedMC=nChargedMCEta10;
        nChargedMCPrimary=nChargedMCPrimaryEta10;
        nChargedMCPhysicalPrimary=nChargedMCPhysicalPrimaryEta10;
        
        // Compute the Nch weights (reference is Ntracklets within |eta|<1.0)
        if(fUseNchWeight>0){
            
            Double_t tmpweight = 1.0;
            Double_t tmpXweight=nChargedMCPhysicalPrimary; // Nch weights
            if(fUseNchWeight==2) tmpXweight=countMult;     // Ntrk weights
            
            if(tmpXweight<=0) tmpweight = 0.0;
            else{
                Double_t pMeas = fHistoMeasNch->GetBinContent(fHistoMeasNch->FindBin(tmpXweight));
                //printf(" pMeas=%2.2f  and histo MCNch %s \n",pMeas,fHistoMCNch);
                Double_t pMC = fHistoMCNch->GetBinContent(fHistoMCNch->FindBin(tmpXweight));
                tmpweight = pMC>0 ? pMeas/pMC : 0.;
            }
            nchWeight *= tmpweight;
            AliDebug(2,Form("Using Nch weights, Mult=%f Weight=%f\n",tmpXweight,nchWeight));
        }
        
        FillMCGenAccHistos(aod, arrayMC, mcHeader, countCorr, spherocity, sphericity, isEvSel, nchWeight);//Fill 2 separate THnSparses, one for prompt andf one for feeddown
        genspherocity=AliVertexingHFUtils::GetGeneratedSpherocity(arrayMC, fetaMin, fetaMax, fptMin, fptMax, fminMult, fphiStepSizeDeg);

    }
    
    if(!isEvSel)return;
    fHistNEvents->Fill(2);
    
    if(vtx1){
        fHistNtrVsZvtx->Fill(vtx1->GetZ(),countMult);
        fHistNtrCorrVsZvtx->Fill(vtx1->GetZ(),countCorr);
        fHistNtrVsSo->Fill(spherocity,countMult);
        fHistNtrCorrVsSo->Fill(spherocity,countCorr);
        if(fCalculateSphericity){
            fHistNtrVsSpheri->Fill(sphericity,countMult);
            fHistNtrCorrVsSpheri->Fill(sphericity,countCorr);
        }
    }
    
    if(fReadMC){
        // Now recompute the variables in case another MC estimator is considered
        if(fMCPrimariesEstimator==kEta10to16){
            nChargedMC = nChargedMCEta16 - nChargedMCEta10;
            nChargedMCPrimary = nChargedMCPrimaryEta16 - nChargedMCPrimaryEta10;
            nChargedMCPhysicalPrimary = nChargedMCPhysicalPrimaryEta16 - nChargedMCPhysicalPrimaryEta10;
        } else if(fMCPrimariesEstimator==kEta05){
            nChargedMC = nChargedMCEta05;
            nChargedMCPrimary = nChargedMCPrimaryEta05;
            nChargedMCPhysicalPrimary = nChargedMCPhysicalPrimaryEta05;
        } else if(fMCPrimariesEstimator==kEta03){
            nChargedMC = nChargedMCEta03;
            nChargedMCPrimary = nChargedMCPrimaryEta03;
            nChargedMCPhysicalPrimary = nChargedMCPhysicalPrimaryEta03;
        } else if(fMCPrimariesEstimator==kEtaVZERO){
            nChargedMC = nChargedMCEtam37tm17 + nChargedMCEta28t51;
            nChargedMCPrimary = nChargedMCPrimaryEtam37tm17 + nChargedMCPrimaryEta28t51;
            nChargedMCPhysicalPrimary = nChargedMCPhysicalPrimaryEtam37tm17 + nChargedMCPhysicalPrimaryEta28t51;
        } else if(fMCPrimariesEstimator==kEtaVZEROA){
            nChargedMC = nChargedMCEta28t51;
            nChargedMCPrimary = nChargedMCPrimaryEta28t51;
            nChargedMCPhysicalPrimary = nChargedMCPhysicalPrimaryEta28t51;
        }
        // Here fill the MC correlation plots
        if(nChargedMCPhysicalPrimary>0){ // INEL>0 for |eta|<1
            fHistGenPrimaryParticlesInelGt0->Fill(nChargedMCPhysicalPrimary,nchWeight);
        }
        
        fHistNtrVsNchMC->Fill(nChargedMC,countMult,nchWeight);
        fHistNtrCorrVsNchMC->Fill(nChargedMC,countCorr,nchWeight);
        
        fHistNtrVsNchMCPrimary->Fill(nChargedMCPrimary,countMult,nchWeight);
        fHistNtrCorrVsNchMCPrimary->Fill(nChargedMCPrimary,countCorr,nchWeight);
        
        fHistNtrVsNchMCPhysicalPrimary->Fill(nChargedMCPhysicalPrimary,countMult,nchWeight);
        fHistNtrCorrVsNchMCPhysicalPrimary->Fill(nChargedMCPhysicalPrimary,countCorr,nchWeight);
        
        fHistNchMCVsNchMCPrimaryVsNchMCPhysicalPrimary->Fill(nChargedMC,nChargedMCPrimary,nChargedMCPhysicalPrimary,nchWeight);
    }
    
    Int_t nCand = arrayCand->GetEntriesFast();
    Int_t nSelectedNoPID=0,nSelectedPID=0,nSelectedInMassPeak=0;
    Double_t mD0PDG = TDatabasePDG::Instance()->GetParticle(421)->Mass();
    Double_t mDplusPDG = TDatabasePDG::Instance()->GetParticle(411)->Mass();
    Double_t mDstarPDG = TDatabasePDG::Instance()->GetParticle(413)->Mass();
    
    // pdg of daughters needed for D* too
    UInt_t pdgDgDStartoD0pi[2]={421,211};
    
    Double_t aveMult=0.;
    Double_t nSelCand=0.;
    for (Int_t iCand = 0; iCand < nCand; iCand++) {
        AliAODRecoDecayHF *d = (AliAODRecoDecayHF*)arrayCand->UncheckedAt(iCand);
        AliAODRecoCascadeHF *dCascade = NULL;
        if(fPdgMeson==413) dCascade = (AliAODRecoCascadeHF*)d;
        
        fHistNEvents->Fill(7);
        if(fUseBit && !d->HasSelectionBit(selbit)){
            fHistNEvents->Fill(8);
            continue;
        }
        
        Double_t ptCand = d->Pt();
        Double_t rapid=d->Y(fPdgMeson);
        Bool_t isFidAcc=fRDCutsAnalysis->IsInFiducialAcceptance(ptCand,rapid);
        if(!isFidAcc) continue;
        
        Int_t passAllCuts=fRDCutsAnalysis->IsSelected(d,AliRDHFCuts::kAll,aod);
        Int_t passTopolCuts=fRDCutsAnalysis->GetIsSelectedCuts();
        if(passTopolCuts==0) continue;
        nSelectedNoPID++;
        fHistNEvents->Fill(9);
        if(passAllCuts){
            nSelectedPID++;
            fHistNEvents->Fill(10);
        }
        Double_t multForCand = countCorr;
        
        //Weight according to ptcand, using FONLL
        Double_t ptWeight = 1.;
        if(fUsePtWeight) ptWeight = GetPtWeight(ptCand);
        fWeight = nchWeight*ptWeight;
        
        if(fSubtractTrackletsFromDau){
            // For the D* case, subtract only the D0 daughter tracks <=== FIXME !!
            AliAODRecoDecayHF2Prong* d0fromDstar = NULL;
            if(fPdgMeson==413) d0fromDstar = (AliAODRecoDecayHF2Prong*)dCascade->Get2Prong();
            
            for(Int_t iDau=0; iDau<nDau; iDau++){
                AliAODTrack *t = NULL;
                if(fPdgMeson==413){ t = (AliAODTrack*)d0fromDstar->GetDaughter(iDau); }
                else{ t = (AliAODTrack*)d->GetDaughter(iDau); }
                if(!t) continue;
                if(t->HasPointOnITSLayer(0) && t->HasPointOnITSLayer(1)){
                    if(multForCand>0) multForCand-=1;
                }
            }
        }
        
        Bool_t isPrimary=kTRUE;
        Double_t trueImpParXY=9999.;
        Double_t impparXY=d->ImpParXY()*10000.;
        Double_t dlen=0.1; //FIXME
        Double_t mass[2];
        if(fPdgMeson==411){
            mass[0]=d->InvMass(nDau,pdgDau);
            mass[1]=-1.;
            if(TMath::Abs(mass[0]-mDplusPDG)<0.02) nSelectedInMassPeak++; //20 MeV for now... FIXME
        }else if(fPdgMeson==421){
            UInt_t pdgdaughtersD0[2]={211,321};//pi,K
            UInt_t pdgdaughtersD0bar[2]={321,211};//K,pi
            mass[0]=d->InvMass(2,pdgdaughtersD0);
            mass[1]=d->InvMass(2,pdgdaughtersD0bar);
            if(TMath::Abs(mass[0]-mD0PDG)<0.02 || TMath::Abs(mass[1]-mD0PDG)<0.02 ) nSelectedInMassPeak++; //20 MeV for now... FIXME
        }else if(fPdgMeson==413){
            // FIXME
            mass[0]=dCascade->DeltaInvMass();
            mass[1]=-1.;
            if(TMath::Abs(mass[0]-(mDstarPDG-mD0PDG))<0.0015) nSelectedInMassPeak++; //1 MeV for now... FIXME
        }
        
        Int_t labD0=-1;
        Double_t recSpherocity = -0.5;
        
        // subtract D-meson daughters from spherocity calculation !!
        if(fRecomputeSpherocity){
            Int_t totTrkToSkip = d->GetNDaughters();
            const Int_t nTrkToSkip = totTrkToSkip;
            Int_t idToSkip[nTrkToSkip];
            for(Int_t i=0; i<nTrkToSkip; i++) idToSkip[i]=-1;
            
            for(Int_t iDau=0; iDau<nTrkToSkip; iDau++){
                AliAODTrack *t = NULL;
                t = dynamic_cast<AliAODTrack*>(d->GetDaughter(iDau));
                if(!t) continue;
                idToSkip[iDau] = t->GetID();
            }
            recSpherocity=AliVertexingHFUtils::GetSpherocity(aod, fetaMin, fetaMax, fptMin, fptMax, ffiltbit1, ffiltbit2, fminMult, fphiStepSizeDeg, nTrkToSkip, idToSkip);
        }
        
        // remove D0 from Dstar at reconstruction !!
        if(fReadMC && fRemoveD0fromDstar){
            if(fPdgMeson==421){
                labD0 = d->MatchToMC(fPdgMeson,arrayMC,nDau,(Int_t*)pdgDau);
                if(labD0>=0){
                    Bool_t keep=kTRUE;
                    AliAODMCParticle* mcMoth = dynamic_cast<AliAODMCParticle*>(arrayMC->At(labD0));
                    Int_t motherD0 = mcMoth->GetMother();
                    AliAODMCParticle* mcMothD0 = dynamic_cast<AliAODMCParticle*>(arrayMC->At(motherD0));
                    if(!mcMothD0) continue;
                    if(TMath::Abs(mcMothD0->GetPdgCode())==413) keep=kFALSE;
                    if(!keep) continue;
                }
            }
        }
        
        Int_t labD=-1;
        Int_t Origin = 0;
        
        for(Int_t iHyp=0; iHyp<2; iHyp++){
            if(mass[iHyp]<0.) continue; // for D+ and D* we have 1 mass hypothesis
            Double_t invMass=mass[iHyp];
            Double_t arrayForSparse[5]={invMass,ptCand,impparXY,dlen,multForCand};
            
            if(fReadMC){
                if(fPdgMeson==413){
                    labD = dCascade->MatchToMC(fPdgMeson,421,(Int_t*)pdgDgDStartoD0pi,(Int_t*)pdgDau,arrayMC);
                } else {
                    labD = d->MatchToMC(fPdgMeson,arrayMC,nDau,(Int_t*)pdgDau);
                }
                
                Bool_t fillHisto=fDoImpPar;
                if(labD>=0){
                    AliAODMCParticle *partD = (AliAODMCParticle*)arrayMC->At(labD);
                    Int_t code=partD->GetPdgCode();
                    Origin = AliVertexingHFUtils::CheckOrigin(arrayMC,partD,fUseQuarkTag);
                    if(Origin==5) isPrimary=kFALSE;
                    if(code<0 && iHyp==0) fillHisto=kFALSE;
                    if(code>0 && iHyp==1) fillHisto=kFALSE;
                    if(!isPrimary){
                        if(fPdgMeson==411){
                            trueImpParXY=AliVertexingHFUtils::GetTrueImpactParameterDplus(mcHeader,arrayMC,partD)*10000.;
                        }else if(fPdgMeson==421){
                            trueImpParXY=AliVertexingHFUtils::GetTrueImpactParameterDzero(mcHeader,arrayMC,partD)*10000.;
                        }else if(fPdgMeson==413){
                            trueImpParXY=0.; /// FIXME
                        }
                        Double_t arrayForSparseTrue[5]={invMass,ptCand,trueImpParXY,dlen,multForCand};
                        if(fillHisto && passAllCuts){
                            fHistMassPtImpPar[2]->Fill(arrayForSparse, fWeight);
                            fHistMassPtImpPar[3]->Fill(arrayForSparseTrue, fWeight);
                        }
                    }else{
                        if(fillHisto && passAllCuts) fHistMassPtImpPar[1]->Fill(arrayForSparse, fWeight);
                    }
                }else{
                    if(fillHisto && passAllCuts)fHistMassPtImpPar[4]->Fill(arrayForSparse, fWeight);
                }
                if(TMath::Abs(labD)==fPdgMeson && fMCOption==2) continue;
                if(TMath::Abs(labD)!=fPdgMeson && fMCOption==1) continue;
            }
            if(fPdgMeson==421){
                if(iHyp==0 && !(passTopolCuts&1)) continue; // candidate not passing as D0
                if(iHyp==1 && !(passTopolCuts&2)) continue; // candidate not passing as D0bar
            }
            if(fFillSoSparseChecks == 2 || fFillSoSparseChecks == 3){   //Filling THnSparse for Spherocity without PID
                if(fCalculateSphericity){
                    Double_t arrayForSparseSoNoPid[5]={ptCand, invMass, spherocity, multForCand, sphericity};
                    fSparseEvtShapewithNoPid->Fill(arrayForSparseSoNoPid, fWeight);
                }else{
                    Double_t arrayForSparseSoNoPid[4]={ptCand, invMass, spherocity, multForCand};
                    fSparseEvtShapewithNoPid->Fill(arrayForSparseSoNoPid, fWeight);
                }
            }
            if(fPdgMeson==421){
                if(iHyp==0 && !(passAllCuts&1)) continue; // candidate not passing as D0
                if(iHyp==1 && !(passAllCuts&2)) continue; // candidate not passing as D0bar
            }
            if(passAllCuts){
                aveMult+=multForCand;
                nSelCand+=1.;
                
                if(fFillSoSparseChecks == 1 || fFillSoSparseChecks == 3){
                    if(fCalculateSphericity){
                        Double_t arrayForSparseSowithMultUnncorr[6]={ptCand, invMass, spherocity, multForCand, (Double_t)countTreta1, sphericity};
                        fSparseEvtShape->Fill(arrayForSparseSowithMultUnncorr, fWeight);
                    }else{
                        Double_t arrayForSparseSowithMultUnncorr[5]={ptCand, invMass, spherocity, multForCand, (Double_t)countTreta1};
                        fSparseEvtShape->Fill(arrayForSparseSowithMultUnncorr, fWeight);
                    }
                }
                else{
                    if(fCalculateSphericity){
                        Double_t arrayForSparseSo[5]={ptCand, invMass, spherocity, multForCand, sphericity};
                        fSparseEvtShape->Fill(arrayForSparseSo, fWeight);
                    }else{
                        Double_t arrayForSparseSo[4]={ptCand, invMass, spherocity, multForCand};
                        fSparseEvtShape->Fill(arrayForSparseSo, fWeight);
                    }
                }
                
                if(fRecomputeSpherocity){
                    Double_t arrayForSparseRecSphero[5]={ptCand, invMass, spherocity, multForCand, recSpherocity};
                    fSparseEvtShapeRecSphero->Fill(arrayForSparseRecSphero, fWeight);
                }
                
                if(fReadMC){
                    if(fRecomputeSpherocity){
                        Double_t arrayForSparseSoPromptFD[5]={ptCand, invMass, spherocity, multForCand, recSpherocity};
                        
                        if(Origin==4) fSparseEvtShapePrompt->Fill(arrayForSparseSoPromptFD, fWeight);
                        else if(Origin==5) fSparseEvtShapeFeeddown->Fill(arrayForSparseSoPromptFD, fWeight);
                    }else{
                        Double_t arrayForSparseSoPromptFD[4]={ptCand, invMass, spherocity, multForCand};
                        
                        if(Origin==4) fSparseEvtShapePrompt->Fill(arrayForSparseSoPromptFD, fWeight);
                        else if(Origin==5) fSparseEvtShapeFeeddown->Fill(arrayForSparseSoPromptFD, fWeight);
                    }
                }
                if(labD>=0){
                    Bool_t keepCase=kTRUE;
                    if(fPdgMeson==421){
                        AliAODMCParticle *partD = (AliAODMCParticle*)arrayMC->At(labD);
                        Int_t code=partD->GetPdgCode();
                        if(code<0 && iHyp==0) keepCase=kFALSE;
                        if(code>0 && iHyp==1) keepCase=kFALSE;
                    }
                    if(keepCase) FillMCMassHistos(arrayMC,labD, multForCand, spherocity, sphericity, recSpherocity, nchWeight);
                }
                if(fDoImpPar) fHistMassPtImpPar[0]->Fill(arrayForSparse, fWeight);
            }
        }
    }
    if(fSubtractTrackletsFromDau && nSelCand>0){
        aveMult/=nSelCand;
        fCounterCandidates->StoreEvent(aod,fRDCutsAnalysis,fReadMC,(Int_t)(aveMult+0.5001),spherocity);
    }else{
        fCounterCandidates->StoreEvent(aod,fRDCutsAnalysis,fReadMC,(Int_t)countCorr,spherocity);
    }
    
    fCounterCandidates->StoreCandidates(aod,nSelectedNoPID,kTRUE);
    fCounterCandidates->StoreCandidates(aod,nSelectedPID,kFALSE);
    fHistNtrCorrEvSel->Fill(countCorr,nchWeight);
    if(nSelectedPID>0)  fHistNtrCorrEvWithCand->Fill(countCorr,nchWeight);
    if(nSelectedInMassPeak>0) fHistNtrCorrEvWithD->Fill(countCorr,nchWeight);
    
    if(fFillTrackHisto && fReadMC) FillTrackControlHisto(aod, countCorr, spherocity, genspherocity, nSelectedInMassPeak);
    
    PostData(1,fOutput);
    PostData(2,fListCuts);
    PostData(3,fOutputCounters);
    if(fReadMC) PostData(5,fOutputEffCorr);
    
    return;
}

//________________________________________________________________________
void AliAnalysisTaskSEDvsEventShapes::CreateImpactParameterHistos()
{
    // Histos for impact paramter study
    // mass . pt , impact parameter , decay length , multiplicity
    
    Int_t nbins[5]={fNMassBins,200,fNImpParBins,50,100};
    Double_t xmin[5]={fLowmasslimit,0.,fLowerImpPar,0.,0.};
    Double_t xmax[5]={fUpmasslimit,20.,fHigherImpPar,1.,100.};
    
    fHistMassPtImpPar[0]=new THnSparseF("hMassPtImpParAll",
                                        "Mass vs. pt vs.imppar - All",
                                        5,nbins,xmin,xmax);
    fHistMassPtImpPar[1]=new THnSparseF("hMassPtImpParPrompt",
                                        "Mass vs. pt vs.imppar - promptD",
                                        5,nbins,xmin,xmax);
    fHistMassPtImpPar[2]=new THnSparseF("hMassPtImpParBfeed",
                                        "Mass vs. pt vs.imppar - DfromB",
                                        5,nbins,xmin,xmax);
    fHistMassPtImpPar[3]=new THnSparseF("hMassPtImpParTrueBfeed",
                                        "Mass vs. pt vs.true imppar -DfromB",
                                        5,nbins,xmin,xmax);
    fHistMassPtImpPar[4]=new THnSparseF("hMassPtImpParBkg",
                                        "Mass vs. pt vs.imppar - backgr.",
                                        5,nbins,xmin,xmax);
    for(Int_t i=0; i<5;i++){
        fOutput->Add(fHistMassPtImpPar[i]);
    }
}

//________________________________________________________________________
void AliAnalysisTaskSEDvsEventShapes::Terminate(Option_t */*option*/)
{
    // Terminate analysis
    //
    if(fDebug > 1) printf("AnalysisTaskSEDvsMultiplicity: Terminate() \n");
    
    fOutput = dynamic_cast<TList*> (GetOutputData(1));
    if (!fOutput) {
        printf("ERROR: fOutput not available\n");
        return;
    }
    
    fHistNEvents = dynamic_cast<TH1F*>(fOutput->FindObject("fHistNEvents"));
    if(!fHistNEvents){
        printf("ERROR: fHistNEvents not available\n");
        return;
    }
    printf("Number of Analyzed Events = %d\n",(Int_t)fHistNEvents->GetBinContent(3));
    
    return;
}

//____________________________________________________________________________
TProfile* AliAnalysisTaskSEDvsEventShapes::GetEstimatorHistogram(const AliVEvent* event)
{
    // Get Estimator Histogram from period event->GetRunNumber();
    //
    // If you select SPD tracklets in |eta|<1 you should use type == 1
    //
    
    Int_t runNo  = event->GetRunNumber();
    Int_t period = -1;   // pp: 0-LHC10b, 1-LHC10c, 2-LHC10d, 3-LHC10e
    // pPb: 0-LHC13b, 1-LHC13c
    if (fisPPbData) {
        if (runNo>195343 && runNo<195484) period = 0;
        if (runNo>195528 && runNo<195678) period = 1;
        if (period < 0 || period > 1) return 0;
    }
    else {
        if(runNo>114930 && runNo<117223) period = 0;
        if(runNo>119158 && runNo<120830) period = 1;
        if(runNo>122373 && runNo<126438) period = 2;
        if(runNo>127711 && runNo<130851) period = 3;
        if(period<0 || period>3) return 0;
    }
    
    return fMultEstimatorAvg[period];
}
//__________________________________________________________________________________________________

Double_t AliAnalysisTaskSEDvsEventShapes::GetPtWeight(Float_t pt)
{
    //
    // calculating the weight to fill the container
    //
    
    // FNOLL central:
    // p0 = 1.63297e-01 --> 0.322643
    // p1 = 2.96275e+00
    // p2 = 2.30301e+00
    // p3 = 2.50000e+00
    
    // PYTHIA
    // p0 = 1.85906e-01 --> 0.36609
    // p1 = 1.94635e+00
    // p2 = 1.40463e+00
    // p3 = 2.50000e+00
    
    Double_t func1[4] = {0.322643,2.96275,2.30301,2.5};
    Double_t func2[4] = {0.36609,1.94635,1.40463,2.5};
    
    Double_t dndpt_func1 = dNdptFit(pt,func1);
    //if(fUseFlatPtWeight) dndpt_func1 = 1./30.;
    Double_t dndpt_func2 = dNdptFit(pt,func2);
    AliDebug(2,Form("pt = %f, FONLL = %f, Pythia = %f, ratio = %f",pt,dndpt_func1,dndpt_func2,dndpt_func1/dndpt_func2));
    return dndpt_func1/dndpt_func2;
}

//__________________________________________________________________________________________________
Double_t AliAnalysisTaskSEDvsEventShapes::dNdptFit(Float_t pt, Double_t* par)
{
    //
    // calculating dNdpt
    //
    
    Double_t denom =  TMath::Power((pt/par[1]), par[3] );
    Double_t dNdpt = par[0]*pt/TMath::Power(1.+denom, par[2]);
    
    return dNdpt;
}


//__________________________________________________________________________________________________
void AliAnalysisTaskSEDvsEventShapes::CreateMeasuredNchHisto()
{
    // creates historgam with measured multiplcity distribution in pp 7 TeV collisions (from Eur. Phys. J. C (2010) 68: 345–354)
    //
    // for Nch  > 70 the points were obtainedwith a double NBD distribution
    // TF1 *fit1 = new TF1("fit1","[0]*(TMath::Gamma(x+[1])/(TMath::Gamma(x+1)*TMath::Gamma([1])))*(TMath::Power(([2]/[1]),x))*(TMath::Power((1+([2]/[1])),-x-[1]))"); fit1->SetParameter(0,1.);// normalization constant
    // fit1->SetParameter(1,1.63); // k parameter
    // fit1->SetParameter(2,12.8); // mean multiplicity
    Double_t nchbins[82]={0.50,1.50,2.50,3.50,4.50,5.50,6.50,7.50,8.50,9.50,
        10.50,11.50,12.50,13.50,14.50,15.50,16.50,17.50,18.50,19.50,
        20.50,21.50,22.50,23.50,24.50,25.50,26.50,27.50,28.50,29.50,
        30.50,31.50,32.50,33.50,34.50,35.50,36.50,37.50,38.50,39.50,
        40.50,41.50,42.50,43.50,44.50,45.50,46.50,47.50,48.50,49.50,
        50.50,51.50,52.50,53.50,54.50,55.50,56.50,57.50,58.50,59.50,
        60.50,62.50,64.50,66.50,68.50,70.50,72.50,74.50,76.50,78.50,
        80.50,82.50,84.50,86.50,88.50,90.50,92.50,94.50,96.50,98.50,
        100.50,102.50};
    Double_t pch[81]={0.062011,0.072943,0.070771,0.067245,0.062834,0.057383,0.051499,0.04591,0.041109,0.036954,
        0.03359,0.030729,0.028539,0.026575,0.024653,0.0229,0.021325,0.019768,0.018561,0.017187,
        0.01604,0.014836,0.013726,0.012576,0.011481,0.010393,0.009502,0.008776,0.008024,0.007452,
        0.006851,0.006428,0.00594,0.005515,0.005102,0.00469,0.004162,0.003811,0.003389,0.003071,
        0.002708,0.002422,0.002184,0.001968,0.00186,0.00165,0.001577,0.001387,0.001254,0.001118,
        0.001037,0.000942,0.000823,0.000736,0.000654,0.000579,0.000512,0.00049,0.00045,0.000355,
        0.000296,0.000265,0.000193,0.00016,0.000126,0.0000851, 0.0000676,0.0000537,0.0000426, 0.0000338,
        0.0000268,0.0000213,0.0000166,0.0000133,0.0000106,0.00000837,0.00000662, 0.00000524,0.00000414, 0.00000327,
        0.00000258};
    
    if(fHistoMeasNch) delete fHistoMeasNch;
    fHistoMeasNch=new TH1F("hMeaseNch","",81,nchbins);
    for(Int_t i=0; i<81; i++){
        fHistoMeasNch->SetBinContent(i+1,pch[i]);
        fHistoMeasNch->SetBinError(i+1,0.);
    }
}

//________________________________________________________________________
Bool_t AliAnalysisTaskSEDvsEventShapes::FillTrackControlHisto(AliAODEvent* aod, Int_t nSelTrkCorr, Double_t spherocity, Double_t genspherocity, Int_t nSelectedEvwithCand)
{
    /// Function to fill the track control histograms
    
    Int_t nTracks=aod->GetNumberOfTracks();
    Int_t nSelTracks=0;
    
    Double_t* etaArr=new Double_t[nTracks];
    Double_t* phiArr=new Double_t[nTracks];
    Double_t sumpt=0.;
    
    for(Int_t it=0; it<nTracks; it++) {
        AliAODTrack *tr=dynamic_cast<AliAODTrack*>(aod->GetTrack(it));
        if(!tr) continue;
        Float_t eta = tr->Eta();
        Float_t pt  = tr->Pt();
        Float_t phi  = tr->Phi();
        if(eta<fetaMin || eta>fetaMax) continue;
        if(pt<fptMin || pt>fptMax) continue;
        Bool_t fb1 = tr->TestFilterBit(ffiltbit1);
        Bool_t fb2 = tr->TestFilterBit(ffiltbit2);
        if( !(fb1 || fb2) ) continue;
        if(ffiltbit1==0 || ffiltbit2==0){
            Int_t status=tr->GetStatus();
            if(!(status & AliAODTrack::kTPCrefit)) continue;
        }
        etaArr[nSelTracks]=eta;
        phiArr[nSelTracks]=phi;
        nSelTracks++;
        sumpt+=pt;
    }
    
    if(nSelTracks<fminMult) return kFALSE;
    
    if(nSelectedEvwithCand>0) fHistNtrVsnTrackEvWithCand->Fill(nSelTracks, nSelTrkCorr);
    
    for(Int_t iselt=0; iselt<nSelTracks; iselt++) {
        fHistnTrackvsEtavsPhi->Fill(etaArr[iselt], phiArr[iselt], nSelTracks);
        fHistTrueSovsMeasSo->Fill(genspherocity, spherocity, nSelTrkCorr);
        if(nSelectedEvwithCand>0){
            fHistnTrackvsEtavsPhiEvWithCand->Fill(etaArr[iselt], phiArr[iselt], nSelTracks);
            fHistTrueSovsMeasSoEvWithCand->Fill(genspherocity, spherocity, nSelTrkCorr);
        }
    }
    
    return kTRUE;
    
}

//__________________________________________________________________________________________________
void AliAnalysisTaskSEDvsEventShapes::FillMCMassHistos(TClonesArray *arrayMC, Int_t labD, Double_t countMult, Double_t spherocity, Double_t sphericity, Double_t recSpherocity, Double_t nchWeight)
{
    // Function to fill the true MC signal
    
    AliAODMCParticle *partD = (AliAODMCParticle*)arrayMC->At(labD);
    
    Double_t mass = partD->M();
    Double_t pt = partD->Pt();
    Double_t rapid = partD->Y();
    
    //Weight according to pt, using FONLL
    Double_t ptWeight = 1.;
    if(fUsePtWeight) ptWeight = GetPtWeight(pt);
    fWeight = nchWeight*ptWeight;
    
    Int_t orig=AliVertexingHFUtils::CheckOrigin(arrayMC,partD,fUseQuarkTag); // Prompt = 4, FeedDown = 5
    
    //for prompt
    if(orig == 4){
        //fill histo for prompt
        Double_t arrayMCRecoRecSpheroPrompt[5] = {pt, countMult, spherocity, rapid, recSpherocity};
        Double_t arrayMCRecoPrompt[4] = {pt, countMult, spherocity, rapid};
        if(fRecomputeSpherocity) fMCRecoPrompt->Fill(arrayMCRecoRecSpheroPrompt, fWeight);
        else fMCRecoPrompt->Fill(arrayMCRecoPrompt, fWeight);
        if(fCalculateSphericity) {
            Double_t arrayMCRecoPromptSpheri[4] = {pt, countMult, sphericity, rapid};
            fMCRecoPromptSpheri->Fill(arrayMCRecoPromptSpheri, fWeight);
        }
    }
    //for FD
    else if(orig == 5){
        //fill histo for FD
        Double_t arrayMCRecoRecSpheroFeeddown[5] = {pt, countMult, spherocity, rapid, recSpherocity};
        Double_t arrayMCRecoFeeddown[4] = {pt, countMult, spherocity, rapid};
        if(fRecomputeSpherocity) fMCRecoFeeddown->Fill(arrayMCRecoRecSpheroFeeddown, fWeight);
        else fMCRecoFeeddown->Fill(arrayMCRecoFeeddown, fWeight);
        if(fCalculateSphericity){
            Double_t arrayMCRecoFeeddownSpheri[4] = {pt, countMult, sphericity, rapid};
            fMCRecoFeeddownSpheri->Fill(arrayMCRecoFeeddownSpheri, fWeight);
        }
    }
    
    Double_t arrayMCReco[4] = {pt, countMult, spherocity, rapid};
    fMCRecoBothPromptFD->Fill(arrayMCReco, fWeight);
    if(fCalculateSphericity){
        Double_t arrayMCRecoSpheri[4] = {pt, countMult, sphericity, rapid};
        fMCRecoBothPromptFDSpheri->Fill(arrayMCRecoSpheri, fWeight);
    }
    
}

//__________________________________________________________________________________________________
void AliAnalysisTaskSEDvsEventShapes::FillMCGenAccHistos(AliAODEvent* aod, TClonesArray *arrayMC, AliAODMCHeader *mcHeader, Double_t countMult, Double_t spherocity, Double_t sphericity, Bool_t isEvSel, Double_t nchWeight)
{
    
    /// Fill MC acceptance histos at generator level
    
    Int_t nProng=2;
    Int_t totPart = arrayMC->GetEntriesFast(); //number of particles
    Int_t totTracks = aod->GetNumberOfTracks(); // number of tracks
    Double_t recSpherocity = -0.5;
    
    const Int_t nPart = totPart;
    Int_t trkToSkip[nPart];
    for(Int_t i=0; i<nPart; i++) trkToSkip[i]=-1; //stores ID of tracks at i'th label of particle
    
    if(fRecomputeSpherocity){
        for(Int_t it=0; it<totTracks; it++){
            AliAODTrack *t=dynamic_cast<AliAODTrack*>(aod->GetTrack(it));
            if(!t) continue;
            if(!(t->TestFilterMask(BIT(4)))) continue;
            //if(!(t->HasPointOnITSLayer(0)) && !(t->HasPointOnITSLayer(1))) continue;
            
            Int_t lab=TMath::Abs(t->GetLabel());
            Int_t id=t->GetID();
            Double_t pt = t->Pt();
            if(id<0) continue;
            if(pt<0.3) continue;
            AliAODMCParticle* genDup = dynamic_cast<AliAODMCParticle*>(arrayMC->At(lab));
            Double_t xver = genDup->Xv();
            Double_t yver = genDup->Yv();
            Double_t rver = TMath::Sqrt(xver*xver + yver*yver);
            if(rver>3) continue;
            trkToSkip[lab] = id;
        }
    }
    
    if(fPdgMeson==421){
        nProng=2;
    }else if(fPdgMeson==411 || fPdgMeson==431){
        nProng=3;
    }
    
    Int_t nTrkToSkip = nProng;
    Int_t idToSkip[nTrkToSkip];
    for(Int_t i=0; i<nTrkToSkip; i++) idToSkip[i]=-1;
    
    Double_t zMCVertex = mcHeader->GetVtxZ(); //vertex MC
    
    for(Int_t iPart=0; iPart<totPart; iPart++){
        AliAODMCParticle* mcGenPart = dynamic_cast<AliAODMCParticle*>(arrayMC->At(iPart));
        
        if (TMath::Abs(mcGenPart->GetPdgCode()) == fPdgMeson){
            Int_t orig=AliVertexingHFUtils::CheckOrigin(arrayMC,mcGenPart,fUseQuarkTag);//Prompt = 4, FeedDown = 5
            
            Int_t deca = 0;
            Bool_t isGoodDecay=kFALSE;
            Int_t labDau[4]={-1,-1,-1,-1};
            Bool_t isInAcc = kFALSE;
            Bool_t isFidAcc = kFALSE;
            
            if(fPdgMeson==421){
                deca=AliVertexingHFUtils::CheckD0Decay(arrayMC,mcGenPart,labDau);
                if(mcGenPart->GetNDaughters()!=2) continue;
                if(deca==1) isGoodDecay=kTRUE;
            }else if(fPdgMeson==411){
                deca=AliVertexingHFUtils::CheckDplusDecay(arrayMC,mcGenPart,labDau);
                if(deca>0) isGoodDecay=kTRUE;
            }else if(fPdgMeson==431){
                deca=AliVertexingHFUtils::CheckDsDecay(arrayMC,mcGenPart,labDau);
                if(deca==1) isGoodDecay=kTRUE;
            }else if(fPdgMeson==413){
                deca=AliVertexingHFUtils::CheckDstarDecay(arrayMC,mcGenPart,labDau);
                if(deca==1) isGoodDecay=kTRUE;
            }
            
            if(labDau[0]==-1){
                continue; //protection against unfilled array of labels
            }
            
            if(fRecomputeSpherocity && isGoodDecay){
                for(Int_t iDau=0; iDau<nTrkToSkip; iDau++){
                    Int_t indexDau = TMath::Abs(mcGenPart->GetDaughter(iDau));  //index of daughter i.e. label
                    idToSkip[iDau] = trkToSkip[indexDau];
                }
                recSpherocity=AliVertexingHFUtils::GetSpherocity(aod, fetaMin, fetaMax, fptMin, fptMax, ffiltbit1, ffiltbit2, fminMult, fphiStepSizeDeg, nTrkToSkip, idToSkip);
            }
            if(fPdgMeson==421){
                //Removal of D0 from Dstar at Generation !!
                if(fRemoveD0fromDstar){
                    Int_t mother = mcGenPart->GetMother();
                    AliAODMCParticle* mcMoth = dynamic_cast<AliAODMCParticle*>(arrayMC->At(mother));
                    if(!mcMoth) continue;
                    if(TMath::Abs(mcMoth->GetPdgCode())==413) continue;
                }
            }
            
            Double_t pt = mcGenPart->Pt();
            Double_t rapid = mcGenPart->Y();
            
            //Weight according to pt, using FONLL
            Double_t ptWeight = 1.;
            if(fUsePtWeight) ptWeight = GetPtWeight(pt);
            fWeight = nchWeight*ptWeight;
            
            isFidAcc=fRDCutsAnalysis->IsInFiducialAcceptance(pt,rapid);
            isInAcc=CheckGenAcc(arrayMC,nProng,labDau);
            
            if(isGoodDecay && TMath::Abs(zMCVertex) < fRDCutsAnalysis->GetMaxVtxZ() && isFidAcc && isInAcc) {
                //for prompt
                if(orig == 4){
                    //fill histo for prompt
                    Double_t arrayMCGenRecSpheroPrompt[5] = {pt, countMult, spherocity, rapid, recSpherocity};
                    Double_t arrayMCGenPrompt[4] = {pt, countMult, spherocity, rapid};
                    if(fRecomputeSpherocity) fMCAccGenPrompt->Fill(arrayMCGenRecSpheroPrompt, fWeight);
                    else fMCAccGenPrompt->Fill(arrayMCGenPrompt, fWeight);
                    if(isEvSel) fMCAccGenPromptEvSel->Fill(arrayMCGenPrompt, fWeight);
                    if(fCalculateSphericity){
                        Double_t arrayMCGenPromptSpheri[4] = {pt, countMult, sphericity, rapid};
                        fMCAccGenPromptSpheri->Fill(arrayMCGenPromptSpheri, fWeight);
                    }
                }
                //for FD
                else if(orig == 5){
                    //fill histo for FD
                    Double_t arrayMCGenRecSpheroFeeddown[5] = {pt, countMult, spherocity, rapid, recSpherocity};
                    Double_t arrayMCGenFeeddown[4] = {pt, countMult, spherocity, rapid};
                    if(fRecomputeSpherocity) fMCAccGenFeeddown->Fill(arrayMCGenRecSpheroFeeddown, fWeight);
                    else fMCAccGenFeeddown->Fill(arrayMCGenFeeddown, fWeight);
                    if(isEvSel) fMCAccGenFeeddownEvSel->Fill(arrayMCGenFeeddown, fWeight);
                    if(fCalculateSphericity){
                        Double_t arrayMCGenFeeddownSpheri[4] = {pt, countMult, sphericity, rapid};
                        fMCAccGenFeeddownSpheri->Fill(arrayMCGenFeeddownSpheri, fWeight);
                    }
                }
                else
                    continue;
            }
        }
    }
}

//_________________________________________________________________
Bool_t AliAnalysisTaskSEDvsEventShapes::CheckGenAcc(TClonesArray* arrayMC, Int_t nProng, Int_t *labDau)
{
    /// check if the decay products are in the good eta and pt range
    for (Int_t iProng = 0; iProng<nProng; iProng++){
        AliAODMCParticle* mcPartDaughter=dynamic_cast<AliAODMCParticle*>(arrayMC->At(labDau[iProng]));
        if(!mcPartDaughter) return kFALSE;
        Double_t eta = mcPartDaughter->Eta();
        Double_t pt = mcPartDaughter->Pt();
        if (TMath::Abs(eta) > fEtaAccCut || pt < fPtAccCut) return kFALSE;
    }
    return kTRUE;
}
