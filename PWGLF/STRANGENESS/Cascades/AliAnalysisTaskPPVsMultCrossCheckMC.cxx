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

// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//
// This task is meant to explore the normalization issues in identified
// particle spectra by looking at monte carlo and checking if spectra
// measured in the INEL>0 (true) and INEL>0 (SPD) classes are radically
// different.
//
// Main Initial objective is to compare:
//
//    (A) Spectra in the INEL>0 (tracklets) event class as a function of
//        V0M percentiles (cuts applied: all analysis cuts)
//
//    (B) Spectra in the true INEL>0 class as a function of V0M percentiles
//
// Further functionality being added on demand.
//
// --- david.dobrigkeit.chinellato@cern.ch
//
// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class TTree;
class TParticle;
class TVector3;

//class AliMCEventHandler;
//class AliMCEvent;
//class AliStack;

class AliESDVertex;
class AliAODVertex;
class AliESDv0;
class AliAODv0;

#include <Riostream.h>
#include "TList.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TFile.h"
#include "THnSparse.h"
#include "TVector3.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TLegend.h"
//#include "AliLog.h"

#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliV0vertexer.h"
#include "AliCascadeVertexer.h"
#include "AliESDpid.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliInputEventHandler.h"
#include "AliAnalysisManager.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliCentrality.h"
#include "AliPPVsMultUtils.h"
#include "AliPWG0Helper.h"
#include "AliCFContainer.h"
#include "AliMultiplicity.h"
#include "AliAODMCParticle.h"
#include "AliESDcascade.h"
#include "AliAODcascade.h"
#include "AliESDUtils.h"
#include "AliGenEventHeader.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisUtils.h"
#include "AliAnalysisTaskPPVsMultCrossCheckMC.h"
////
#include "AliMultEstimator.h"
#include "AliMultVariable.h"
#include "AliMultInput.h"
#include "AliMultSelection.h"


using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskPPVsMultCrossCheckMC)

AliAnalysisTaskPPVsMultCrossCheckMC::AliAnalysisTaskPPVsMultCrossCheckMC()
: AliAnalysisTaskSE(), fListHist(0), fPIDResponse(0), fESDtrackCuts(0), fPPVsMultUtils(0), fUtils(0),
fHistEventCounter(0),lPureMonteCarlo(kFALSE), fCheckVtxZMC(kTRUE), fAlternateMCSelection(kFALSE), fSkipPS(kFALSE), fUseRecoVtxZ(kFALSE),
fHistV0M_DataSelection(0), fHistV0M_MCSelection(0), fHistV0MAmplitude_DataSelection(0), fHistV0MAmplitude_MCSelection(0),
fHistV0MTrue_DataSelection(0), fHistV0MTrue_MCSelection(0),
fHistV0MVsMidRapidityTrue_DataSelection(0), fHistV0MAmplitudeVsMidRapidityTrue_DataSelection(0), fHistV0MTrueVsMidRapidityTrue_DataSelection(0),
fHistV0MVsMidRapidityTrue_MCSelection(0), fHistV0MAmplitudeVsMidRapidityTrue_MCSelection(0), fHistV0MTrueVsMidRapidityTrue_MCSelection(0),
/////
fHistTracklets08Cent_DataSelection(0), fHistTracklets08Cent_MCSelection(0), fHistTracklets08_DataSelection(0), fHistTracklets08_MCSelection(0),
fHistTracklets08True_DataSelection(0), fHistTracklets08True_MCSelection(0),
fHistTracklets08CentVsMidRapidityTrue_DataSelection(0), fHistTracklets08VsMidRapidityTrue_DataSelection(0), fHistTracklets08TrueVsMidRapidityTrue_DataSelection(0),
fHistTracklets08CentVsMidRapidityTrue_MCSelection(0), fHistTracklets08VsMidRapidityTrue_MCSelection(0), fHistTracklets08TrueVsMidRapidityTrue_MCSelection(0),
/////
fHistTracklets0815Cent_DataSelection(0), fHistTracklets0815Cent_MCSelection(0), fHistTracklets0815_DataSelection(0), fHistTracklets0815_MCSelection(0),
fHistTracklets0815True_DataSelection(0), fHistTracklets0815True_MCSelection(0),
fHistTracklets0815CentVsMidRapidityTrue_DataSelection(0), fHistTracklets0815VsMidRapidityTrue_DataSelection(0), fHistTracklets0815TrueVsMidRapidityTrue_DataSelection(0),
fHistTracklets0815CentVsMidRapidityTrue_MCSelection(0), fHistTracklets0815VsMidRapidityTrue_MCSelection(0), fHistTracklets0815TrueVsMidRapidityTrue_MCSelection(0),
////
fkMultSelection ( kFALSE ), fTrigType(AliVEvent::kMB), fTrigName(""), fkSelectTriggerByName ( kFALSE )
{
    //Empty constructor (not to be used, always pass name...)
    for(Int_t ih=0; ih<9; ih++){
        fHistPt_Generated[ih] = 0x0;
        fHistPt_DataSelection[ih] = 0x0;
        fHistPt_MCSelection[ih] = 0x0;
        fHistPtVsV0M_Generated[ih] = 0x0;
        fHistPtVsV0M_DataSelection[ih] = 0x0;
        fHistPtVsV0M_MCSelection[ih] = 0x0;
        fHistPtVsV0MAmplitude_Generated[ih] = 0x0;
        fHistPtVsV0MAmplitude_DataSelection[ih] = 0x0;
        fHistPtVsV0MAmplitude_MCSelection[ih] = 0x0;
        fHistPtVsV0MTrue_Generated[ih] = 0x0;
        fHistPtVsV0MTrue_DataSelection[ih] = 0x0;
        fHistPtVsV0MTrue_MCSelection[ih] = 0x0;
        /////
        fHistPtVsTracklets08Cent_Generated[ih] = 0x0;
        fHistPtVsTracklets08Cent_DataSelection[ih] = 0x0;
        fHistPtVsTracklets08Cent_MCSelection[ih] = 0x0;
        fHistPtVsTracklets08_Generated[ih] = 0x0;
        fHistPtVsTracklets08_DataSelection[ih] = 0x0;
        fHistPtVsTracklets08_MCSelection[ih] = 0x0;
        fHistPtVsTracklets08True_Generated[ih] = 0x0;
        fHistPtVsTracklets08True_DataSelection[ih] = 0x0;
        fHistPtVsTracklets08True_MCSelection[ih] = 0x0; 
        /////
        fHistPtVsTracklets0815Cent_Generated[ih] = 0x0;
        fHistPtVsTracklets0815Cent_DataSelection[ih] = 0x0;
        fHistPtVsTracklets0815Cent_MCSelection[ih] = 0x0;
        fHistPtVsTracklets0815_Generated[ih] = 0x0;
        fHistPtVsTracklets0815_DataSelection[ih] = 0x0;
        fHistPtVsTracklets0815_MCSelection[ih] = 0x0;
        fHistPtVsTracklets0815True_Generated[ih] = 0x0;
        fHistPtVsTracklets0815True_DataSelection[ih] = 0x0;
        fHistPtVsTracklets0815True_MCSelection[ih] = 0x0;
   }
}

AliAnalysisTaskPPVsMultCrossCheckMC::AliAnalysisTaskPPVsMultCrossCheckMC(const char *name)
: AliAnalysisTaskSE(name), fListHist(0), fPIDResponse(0), fESDtrackCuts(0), fPPVsMultUtils(0), fUtils(0),
fHistEventCounter(0),lPureMonteCarlo(kFALSE), fCheckVtxZMC(kTRUE), fAlternateMCSelection(kFALSE), fSkipPS(kFALSE), fUseRecoVtxZ(kFALSE),
fHistV0M_DataSelection(0), fHistV0M_MCSelection(0), fHistV0MAmplitude_DataSelection(0), fHistV0MAmplitude_MCSelection(0),
fHistV0MTrue_DataSelection(0), fHistV0MTrue_MCSelection(0),
fHistV0MVsMidRapidityTrue_DataSelection(0), fHistV0MAmplitudeVsMidRapidityTrue_DataSelection(0), fHistV0MTrueVsMidRapidityTrue_DataSelection(0),
fHistV0MVsMidRapidityTrue_MCSelection(0), fHistV0MAmplitudeVsMidRapidityTrue_MCSelection(0), fHistV0MTrueVsMidRapidityTrue_MCSelection(0),
/////
fHistTracklets08Cent_DataSelection(0), fHistTracklets08Cent_MCSelection(0), fHistTracklets08_DataSelection(0), fHistTracklets08_MCSelection(0),
fHistTracklets08True_DataSelection(0), fHistTracklets08True_MCSelection(0),
fHistTracklets08CentVsMidRapidityTrue_DataSelection(0), fHistTracklets08VsMidRapidityTrue_DataSelection(0), fHistTracklets08TrueVsMidRapidityTrue_DataSelection(0),
fHistTracklets08CentVsMidRapidityTrue_MCSelection(0), fHistTracklets08VsMidRapidityTrue_MCSelection(0), fHistTracklets08TrueVsMidRapidityTrue_MCSelection(0),
/////
fHistTracklets0815Cent_DataSelection(0), fHistTracklets0815Cent_MCSelection(0), fHistTracklets0815_DataSelection(0), fHistTracklets0815_MCSelection(0),
fHistTracklets0815True_DataSelection(0), fHistTracklets0815True_MCSelection(0),
fHistTracklets0815CentVsMidRapidityTrue_DataSelection(0), fHistTracklets0815VsMidRapidityTrue_DataSelection(0), fHistTracklets0815TrueVsMidRapidityTrue_DataSelection(0),
fHistTracklets0815CentVsMidRapidityTrue_MCSelection(0), fHistTracklets0815VsMidRapidityTrue_MCSelection(0), fHistTracklets0815TrueVsMidRapidityTrue_MCSelection(0),
////
fkMultSelection ( kFALSE ), fTrigType(AliVEvent::kMB), fTrigName(""), fkSelectTriggerByName ( kFALSE )
{
    for(Int_t ih=0; ih<9; ih++){
        fHistPt_Generated[ih] = 0x0;
        fHistPt_DataSelection[ih] = 0x0;
        fHistPt_MCSelection[ih] = 0x0;
        fHistPtVsV0M_Generated[ih] = 0x0;
        fHistPtVsV0M_DataSelection[ih] = 0x0;
        fHistPtVsV0M_MCSelection[ih] = 0x0;
        fHistPtVsV0MAmplitude_Generated[ih] = 0x0;
        fHistPtVsV0MAmplitude_DataSelection[ih] = 0x0;
        fHistPtVsV0MAmplitude_MCSelection[ih] = 0x0;
        fHistPtVsV0MTrue_Generated[ih] = 0x0;
        fHistPtVsV0MTrue_DataSelection[ih] = 0x0;
        fHistPtVsV0MTrue_MCSelection[ih] = 0x0;
        /////
        fHistPtVsTracklets08Cent_Generated[ih] = 0x0;
        fHistPtVsTracklets08Cent_DataSelection[ih] = 0x0;
        fHistPtVsTracklets08Cent_MCSelection[ih] = 0x0;
        fHistPtVsTracklets08_Generated[ih] = 0x0;
        fHistPtVsTracklets08_DataSelection[ih] = 0x0;
        fHistPtVsTracklets08_MCSelection[ih] = 0x0;
        fHistPtVsTracklets08True_Generated[ih] = 0x0;
        fHistPtVsTracklets08True_DataSelection[ih] = 0x0;
        fHistPtVsTracklets08True_MCSelection[ih] = 0x0;
        /////
        fHistPtVsTracklets0815Cent_Generated[ih] = 0x0;
        fHistPtVsTracklets0815Cent_DataSelection[ih] = 0x0;
        fHistPtVsTracklets0815Cent_MCSelection[ih] = 0x0;
        fHistPtVsTracklets0815_Generated[ih] = 0x0;
        fHistPtVsTracklets0815_DataSelection[ih] = 0x0;
        fHistPtVsTracklets0815_MCSelection[ih] = 0x0;
        fHistPtVsTracklets0815True_Generated[ih] = 0x0;
        fHistPtVsTracklets0815True_DataSelection[ih] = 0x0;
        fHistPtVsTracklets0815True_MCSelection[ih] = 0x0;
    }
    DefineOutput(1, TList::Class()); // Event Counter Histo
}


AliAnalysisTaskPPVsMultCrossCheckMC::~AliAnalysisTaskPPVsMultCrossCheckMC()
{
    //------------------------------------------------
    // DESTRUCTOR
    //------------------------------------------------
    
    if (fListHist) {
        delete fListHist;
        fListHist = 0x0;
    }
    if (fPPVsMultUtils) {
        delete fPPVsMultUtils;
        fPPVsMultUtils = 0x0;
    }
    if (fUtils) {
        delete fUtils;
        fUtils = 0x0;
    }
}

//________________________________________________________________________
void AliAnalysisTaskPPVsMultCrossCheckMC::UserCreateOutputObjects()
{
    //------------------------------------------------
    // Particle Identification Setup
    //------------------------------------------------
    
    if ( !lPureMonteCarlo ){
        AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
        AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
        fPIDResponse = inputHandler->GetPIDResponse();
        inputHandler->SetNeedField();
    }
    
    // Multiplicity
    if(! fESDtrackCuts ) {
        fESDtrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE,kFALSE);
        fESDtrackCuts->SetPtRange(0.15);  // adding pt cut
        fESDtrackCuts->SetEtaRange(-1.0, 1.0);
    }
    //Helper
    if(! fPPVsMultUtils ) {
        fPPVsMultUtils = new AliPPVsMultUtils();
    }
    //Analysis Utils
    if(! fUtils ) {
        fUtils = new AliAnalysisUtils();
    }
    
    //------------------------------------------------
    // Histograms: Basic Analysis Output
    //------------------------------------------------
    
    // Create histograms
    fListHist = new TList();
    fListHist->SetOwner();  // See http://root.cern.ch/root/html/TCollection.html#TCollection:SetOwner
    
    if(! fHistEventCounter ) {
        //Histogram Output: Event-by-Event
        fHistEventCounter = new TH1F( "fHistEventCounter", ";Evt. Sel. Step;Count",3,0,3);
        //Keeps track of some basics
        fHistEventCounter->GetXaxis()->SetBinLabel(1, "Processed");
        fHistEventCounter->GetXaxis()->SetBinLabel(2, "Data Selection");
        fHistEventCounter->GetXaxis()->SetBinLabel(3, "MC Selection");
        fListHist->Add(fHistEventCounter);
    }
    
    //Identified Particles
    Int_t lPDGCodes[9] = {211, 321, 2212, 310, 3122, 3312, 3334, 333, 313};
    TString lPartNames[9] = {
        "Pion", "Kaon", "Proton", "K0Short", "Lambda", "Xi", "Omega", "Phi", "KStar"
    };
    
    //Settings for transverse momentum
    Int_t lNPtBins = 300; //50MeV/c precision
    Double_t lMaxPt = 15.0;
    
    //Settings for charged particle counters (integers!)
    Int_t lNNchBins = 400;
    Double_t lLowNchBound  = -0.5;
    Double_t lHighNchBound = -0.5 + ((double)(lNNchBins));
    
    //Settings for V0M amplitudes (floating points)
    Int_t lNAmplitudeBins = 2000;
    Double_t lUpperAmplitude = 2000; //may require tuning
 
    //Setting for tracklets
    Int_t nTracklBins = 200;
    Double_t lUpperLimit = 200;
   
    //Main Output: Histograms
    // V0M
    if(! fHistV0M_DataSelection ) { 
        fHistV0M_DataSelection = new TH1F("fHistV0M_DataSelection","",100,0,100);
        fListHist->Add(fHistV0M_DataSelection);
    }
    if(! fHistV0M_MCSelection ) { 
        fHistV0M_MCSelection = new TH1F("fHistV0M_MCSelection","",100,0,100);
        fListHist->Add(fHistV0M_MCSelection);
    }
    if(! fHistV0MAmplitude_DataSelection ) { 
        fHistV0MAmplitude_DataSelection = new TH1F("fHistV0MAmplitude_DataSelection","",lNAmplitudeBins,0,lUpperAmplitude);
        fListHist->Add(fHistV0MAmplitude_DataSelection);
    }
    if(! fHistV0MAmplitude_MCSelection ) { 
        fHistV0MAmplitude_MCSelection = new TH1F("fHistV0MAmplitude_MCSelection","",lNAmplitudeBins,0,lUpperAmplitude);
        fListHist->Add(fHistV0MAmplitude_MCSelection);
    }
    if(! fHistV0MTrue_DataSelection ) { 
        fHistV0MTrue_DataSelection = new TH1F("fHistV0MTrue_DataSelection","",lNNchBins,lLowNchBound,lHighNchBound);
        fListHist->Add(fHistV0MTrue_DataSelection);
    }
    if(! fHistV0MTrue_MCSelection ) { 
        fHistV0MTrue_MCSelection = new TH1F("fHistV0MTrue_MCSelection","",lNNchBins,lLowNchBound,lHighNchBound);
        fListHist->Add(fHistV0MTrue_MCSelection);
    }
    // |Ntrckl| < 0.8
    if(! fHistTracklets08Cent_DataSelection ) { 
        fHistTracklets08Cent_DataSelection = new TH1F("fHistTracklets08Cent_DataSelection","",100,0,100);
        fListHist->Add(fHistTracklets08Cent_DataSelection);
    }
    if(! fHistTracklets08Cent_MCSelection ) {
        fHistTracklets08Cent_MCSelection = new TH1F("fHistTracklets08Cent_MCSelection","",100,0,100);
        fListHist->Add(fHistTracklets08Cent_MCSelection);
    }
    if(! fHistTracklets08_DataSelection ) { 
        fHistTracklets08_DataSelection = new TH1F("fHistTracklets08_DataSelection","",nTracklBins,0,lUpperLimit);
        fListHist->Add(fHistTracklets08_DataSelection);
    }
    if(! fHistTracklets08_MCSelection ) { 
        fHistTracklets08_MCSelection = new TH1F("fHistTracklets08_MCSelection","",nTracklBins,0,lUpperLimit);
        fListHist->Add(fHistTracklets08_MCSelection);
    }
    if(! fHistTracklets08True_DataSelection ) { 
        fHistTracklets08True_DataSelection = new TH1F("fHistTracklets08True_DataSelection","",lNNchBins,lLowNchBound,lHighNchBound);
        fListHist->Add(fHistTracklets08True_DataSelection);
    }
    if(! fHistTracklets08True_MCSelection ) { 
        fHistTracklets08True_MCSelection = new TH1F("fHistTracklets08True_MCSelection","",lNNchBins,lLowNchBound,lHighNchBound);
        fListHist->Add(fHistTracklets08True_MCSelection);
    }
    // 0.8 < |Ntrckl| < 1.5
    if(! fHistTracklets0815Cent_DataSelection ) { 
        fHistTracklets0815Cent_DataSelection = new TH1F("fHistTracklets0815Cent_DataSelection","",100,0,100);
        fListHist->Add(fHistTracklets0815Cent_DataSelection);
    }
    if(! fHistTracklets0815Cent_MCSelection ) { 
        fHistTracklets0815Cent_MCSelection = new TH1F("fHistTracklets0815Cent_MCSelection","",100,0,100);
        fListHist->Add(fHistTracklets0815Cent_MCSelection);
    }
    if(! fHistTracklets0815_DataSelection ) { 
        fHistTracklets0815_DataSelection = new TH1F("fHistTracklets0815_DataSelection","",nTracklBins,0,lUpperLimit);
        fListHist->Add(fHistTracklets0815_DataSelection);
    }
    if(! fHistTracklets0815_MCSelection ) { 
        fHistTracklets0815_MCSelection = new TH1F("fHistTracklets0815_MCSelection","",nTracklBins,0,lUpperLimit);
        fListHist->Add(fHistTracklets0815_MCSelection);
    }
    if(! fHistTracklets0815True_DataSelection ) { 
        fHistTracklets0815True_DataSelection = new TH1F("fHistTracklets0815True_DataSelection","",lNNchBins,lLowNchBound,lHighNchBound);
        fListHist->Add(fHistTracklets0815True_DataSelection);
    }
    if(! fHistTracklets0815True_MCSelection ) { 
        fHistTracklets0815True_MCSelection = new TH1F("fHistTracklets0815True_MCSelection","",lNNchBins,lLowNchBound,lHighNchBound);
        fListHist->Add(fHistTracklets0815True_MCSelection);
    }

    // Correlations V0M 
    //Correlation between mid-rapidity and forward V0M percentile
    if(! fHistV0MVsMidRapidityTrue_DataSelection ) {
        fHistV0MVsMidRapidityTrue_DataSelection = new TH2F("fHistV0MVsMidRapidityTrue_DataSelection","",100,0,100,lNNchBins,lLowNchBound,lHighNchBound);
        fListHist->Add(fHistV0MVsMidRapidityTrue_DataSelection);
    }
    //Correlation between mid-rapidity and forward V0M percentile
    if(! fHistV0MAmplitudeVsMidRapidityTrue_DataSelection ) {
        fHistV0MAmplitudeVsMidRapidityTrue_DataSelection = new TH2F("fHistV0MAmplitudeVsMidRapidityTrue_DataSelection","",lNAmplitudeBins,0,lUpperAmplitude,lNNchBins,lLowNchBound,lHighNchBound);
        fListHist->Add(fHistV0MAmplitudeVsMidRapidityTrue_DataSelection);
    }
    //Correlation between mid-rapidity and forward true multiplicity
    if(! fHistV0MTrueVsMidRapidityTrue_DataSelection ) {
        fHistV0MTrueVsMidRapidityTrue_DataSelection = new TH2F("fHistV0MTrueVsMidRapidityTrue_DataSelection","",lNNchBins,lLowNchBound,lHighNchBound,lNNchBins,lLowNchBound,lHighNchBound);
        fListHist->Add(fHistV0MTrueVsMidRapidityTrue_DataSelection);
    }
    
    //Correlation between mid-rapidity and forward V0M percentile
    if(! fHistV0MVsMidRapidityTrue_MCSelection ) {
        fHistV0MVsMidRapidityTrue_MCSelection = new TH2F("fHistV0MVsMidRapidityTrue_MCSelection","",100,0,100,lNNchBins,lLowNchBound,lHighNchBound);
        fListHist->Add(fHistV0MVsMidRapidityTrue_MCSelection);
    }
    //Correlation between mid-rapidity and forward V0M percentile
    if(! fHistV0MAmplitudeVsMidRapidityTrue_MCSelection ) {
        fHistV0MAmplitudeVsMidRapidityTrue_MCSelection = new TH2F("fHistV0MAmplitudeVsMidRapidityTrue_MCSelection","",lNAmplitudeBins,0,lUpperAmplitude,lNNchBins,lLowNchBound,lHighNchBound);
        fListHist->Add(fHistV0MAmplitudeVsMidRapidityTrue_MCSelection);
    }
    //Correlation between mid-rapidity and forward true multiplicity
    if(! fHistV0MTrueVsMidRapidityTrue_MCSelection ) {
        fHistV0MTrueVsMidRapidityTrue_MCSelection = new TH2F("fHistV0MTrueVsMidRapidityTrue_MCSelection","",lNNchBins,lLowNchBound,lHighNchBound,lNNchBins,lLowNchBound,lHighNchBound);
        fListHist->Add(fHistV0MTrueVsMidRapidityTrue_MCSelection);
    }
    
    // Correlations |Ntrackl| < 0.8
    //Correlation between mid-rapidity and central Ntrkl08 percentile
    if(! fHistTracklets08CentVsMidRapidityTrue_DataSelection ) {
        fHistTracklets08CentVsMidRapidityTrue_DataSelection = new TH2F("fHistTracklets08CentVsMidRapidityTrue_DataSelection","",100,0,100,lNNchBins,lLowNchBound,lHighNchBound);
        fListHist->Add(fHistTracklets08CentVsMidRapidityTrue_DataSelection);
    }
    //Correlation between mid-rapidity and central Ntrkl08 percentile
    if(! fHistTracklets08VsMidRapidityTrue_DataSelection ) {
        fHistTracklets08VsMidRapidityTrue_DataSelection = new TH2F("fHistTracklets08VsMidRapidityTrue_DataSelection","",nTracklBins,0,lUpperLimit,lNNchBins,lLowNchBound,lHighNchBound);
        fListHist->Add(fHistTracklets08VsMidRapidityTrue_DataSelection);
    }
    //Correlation between mid-rapidity and central true (|eta|<0.8) multiplicity
    if(! fHistTracklets08TrueVsMidRapidityTrue_DataSelection ) {
        fHistTracklets08TrueVsMidRapidityTrue_DataSelection = new TH2F("fHistTracklets08TrueVsMidRapidityTrue_DataSelection","",lNNchBins,lLowNchBound,lHighNchBound,lNNchBins,lLowNchBound,lHighNchBound);
        fListHist->Add(fHistTracklets08TrueVsMidRapidityTrue_DataSelection);
    }

    //Correlation between mid-rapidity and central Ntrkl08 percentile
    if(! fHistTracklets08CentVsMidRapidityTrue_MCSelection ) {
        fHistTracklets08CentVsMidRapidityTrue_MCSelection = new TH2F("fHistTracklets08CentVsMidRapidityTrue_MCSelection","",100,0,100,lNNchBins,lLowNchBound,lHighNchBound);
        fListHist->Add(fHistTracklets08CentVsMidRapidityTrue_MCSelection);
    }
    //Correlation between mid-rapidity and central Ntrkl08 percentile
    if(! fHistTracklets08VsMidRapidityTrue_MCSelection ) {
        fHistTracklets08VsMidRapidityTrue_MCSelection = new TH2F("fHistTracklets08VsMidRapidityTrue_MCSelection","",nTracklBins,0,lUpperLimit,lNNchBins,lLowNchBound,lHighNchBound);
        fListHist->Add(fHistTracklets08VsMidRapidityTrue_MCSelection);
    }
    //Correlation between mid-rapidity and central true (|eta|<0.8) multiplicity
    if(! fHistTracklets08TrueVsMidRapidityTrue_MCSelection ) {
        fHistTracklets08TrueVsMidRapidityTrue_MCSelection = new TH2F("fHistTracklets08TrueVsMidRapidityTrue_MCSelection","",lNNchBins,lLowNchBound,lHighNchBound,lNNchBins,lLowNchBound,lHighNchBound);
        fListHist->Add(fHistTracklets08TrueVsMidRapidityTrue_MCSelection);
    }
    
    // Correlations 0.8 < |Ntrackl| < 1.5
    //Correlation between mid-rapidity and central Ntrkl0815 percentile
    if(! fHistTracklets0815CentVsMidRapidityTrue_DataSelection ) {
        fHistTracklets0815CentVsMidRapidityTrue_DataSelection = new TH2F("fHistTracklets0815CentVsMidRapidityTrue_DataSelection","",100,0,100,lNNchBins,lLowNchBound,lHighNchBound);
        fListHist->Add(fHistTracklets0815CentVsMidRapidityTrue_DataSelection);
    }
    //Correlation between mid-rapidity and central Ntrkl0815 percentile
    if(! fHistTracklets0815VsMidRapidityTrue_DataSelection ) {
        fHistTracklets0815VsMidRapidityTrue_DataSelection = new TH2F("fHistTracklets0815VsMidRapidityTrue_DataSelection","",nTracklBins,0,lUpperLimit,lNNchBins,lLowNchBound,lHighNchBound);
        fListHist->Add(fHistTracklets0815VsMidRapidityTrue_DataSelection);
    }
    //Correlation between mid-rapidity and central true (|eta|<0.8) multiplicity
    if(! fHistTracklets0815TrueVsMidRapidityTrue_DataSelection ) {
        fHistTracklets0815TrueVsMidRapidityTrue_DataSelection = new TH2F("fHistTracklets0815TrueVsMidRapidityTrue_DataSelection","",lNNchBins,lLowNchBound,lHighNchBound,lNNchBins,lLowNchBound,lHighNchBound);
        fListHist->Add(fHistTracklets0815TrueVsMidRapidityTrue_DataSelection);
    }

    //Correlation between mid-rapidity and central Ntrkl0815 percentile
    if(! fHistTracklets0815CentVsMidRapidityTrue_MCSelection ) {
        fHistTracklets0815CentVsMidRapidityTrue_MCSelection = new TH2F("fHistTracklets0815CentVsMidRapidityTrue_MCSelection","",100,0,100,lNNchBins,lLowNchBound,lHighNchBound);
        fListHist->Add(fHistTracklets0815CentVsMidRapidityTrue_MCSelection);
    }
    //Correlation between mid-rapidity and central Ntrkl0815 percentile
    if(! fHistTracklets0815VsMidRapidityTrue_MCSelection ) {
        fHistTracklets0815VsMidRapidityTrue_MCSelection = new TH2F("fHistTracklets0815VsMidRapidityTrue_MCSelection","",nTracklBins,0,lUpperLimit,lNNchBins,lLowNchBound,lHighNchBound);
        fListHist->Add(fHistTracklets0815VsMidRapidityTrue_MCSelection);
    }
    //Correlation between mid-rapidity and central true (|eta|<0.8) multiplicity
    if(! fHistTracklets0815TrueVsMidRapidityTrue_MCSelection ) {
        fHistTracklets0815TrueVsMidRapidityTrue_MCSelection = new TH2F("fHistTracklets0815TrueVsMidRapidityTrue_MCSelection","",lNNchBins,lLowNchBound,lHighNchBound,lNNchBins,lLowNchBound,lHighNchBound);
        fListHist->Add(fHistTracklets0815TrueVsMidRapidityTrue_MCSelection);
    }
 
    //Main Output: Histograms
    for(Int_t ih=0; ih<9; ih++){
        if(! fHistPt_Generated[ih] ) {
            fHistPt_Generated[ih] = new TH1F(Form("fHistPt_Generated_%s",lPartNames[ih].Data()),    "Generated;p_{T} (GeV/c)",lNPtBins,0,lMaxPt);
            fListHist->Add(fHistPt_Generated[ih]);
        }
    }
    for(Int_t ih=0; ih<9; ih++){
        if(! fHistPt_DataSelection[ih] ) {
            fHistPt_DataSelection[ih] = new TH1F(Form("fHistPt_DataSelection_%s",lPartNames[ih].Data()),    "Generated;p_{T} (GeV/c)",lNPtBins,0,lMaxPt);
            fListHist->Add(fHistPt_DataSelection[ih]);
        }
    }
    for(Int_t ih=0; ih<9; ih++){
        if(! fHistPt_MCSelection[ih] ) {
            fHistPt_MCSelection[ih] = new TH1F(Form("fHistPt_MCSelection_%s",lPartNames[ih].Data()),    "Generated;p_{T} (GeV/c)",lNPtBins,0,lMaxPt);
            fListHist->Add(fHistPt_MCSelection[ih]);
        }
    }
    
    //2-Dimensional Histograms
    for(Int_t ih=0; ih<9; ih++){
        if(! fHistPtVsV0M_Generated[ih] ) {
            fHistPtVsV0M_Generated[ih] = new TH2F(Form("fHistPtVsV0M_Generated_%s",lPartNames[ih].Data()),    "Generated;p_{T} (GeV/c)",lNPtBins,0,lMaxPt,100,0,100);
            fListHist->Add(fHistPtVsV0M_Generated[ih]);
        }
    }
    for(Int_t ih=0; ih<9; ih++){
        if(! fHistPtVsV0M_DataSelection[ih] ) {
            fHistPtVsV0M_DataSelection[ih] = new TH2F(Form("fHistPtVsV0M_DataSelection_%s",lPartNames[ih].Data()),    "Generated;p_{T} (GeV/c)",lNPtBins,0,lMaxPt,100,0,100);
            fListHist->Add(fHistPtVsV0M_DataSelection[ih]);
        }
    }
    for(Int_t ih=0; ih<9; ih++){
        if(! fHistPtVsV0M_MCSelection[ih] ) {
            fHistPtVsV0M_MCSelection[ih] = new TH2F(Form("fHistPtVsV0M_MCSelection_%s",lPartNames[ih].Data()),    "Generated;p_{T} (GeV/c)",lNPtBins,0,lMaxPt,100,0,100);
            fListHist->Add(fHistPtVsV0M_MCSelection[ih]);
        }
    }
   
 
    //2-Dimensional Histograms with V0M Amplitudes
    for(Int_t ih=0; ih<9; ih++){
        if(! fHistPtVsV0MAmplitude_Generated[ih] ) {
            fHistPtVsV0MAmplitude_Generated[ih] = new TH2F(Form("fHistPtVsV0MAmplitude_Generated_%s",lPartNames[ih].Data()),    "Generated;p_{T} (GeV/c)",lNPtBins,0,lMaxPt,lNAmplitudeBins,0,lUpperAmplitude);
            fListHist->Add(fHistPtVsV0MAmplitude_Generated[ih]);
        }
    }
    for(Int_t ih=0; ih<9; ih++){
        if(! fHistPtVsV0MAmplitude_DataSelection[ih] ) {
            fHistPtVsV0MAmplitude_DataSelection[ih] = new TH2F(Form("fHistPtVsV0MAmplitude_DataSelection_%s",lPartNames[ih].Data()),    "Generated;p_{T} (GeV/c)",lNPtBins,0,lMaxPt,lNAmplitudeBins,0,lUpperAmplitude);
            fListHist->Add(fHistPtVsV0MAmplitude_DataSelection[ih]);
        }
    }
    for(Int_t ih=0; ih<9; ih++){
        if(! fHistPtVsV0MAmplitude_MCSelection[ih] ) {
            fHistPtVsV0MAmplitude_MCSelection[ih] = new TH2F(Form("fHistPtVsV0MAmplitude_MCSelection_%s",lPartNames[ih].Data()),    "Generated;p_{T} (GeV/c)",lNPtBins,0,lMaxPt,lNAmplitudeBins,0,lUpperAmplitude);
            fListHist->Add(fHistPtVsV0MAmplitude_MCSelection[ih]);
        }
    }
    
    //2-Dimensional Histograms with True V0M Multiplicity
    for(Int_t ih=0; ih<9; ih++){
        if(! fHistPtVsV0MTrue_Generated[ih] ) {
            fHistPtVsV0MTrue_Generated[ih] = new TH2F(Form("fHistPtVsV0MTrue_Generated_%s",lPartNames[ih].Data()),    "Generated;p_{T} (GeV/c)",lNPtBins,0,lMaxPt,lNNchBins,lLowNchBound,lHighNchBound);
            fListHist->Add(fHistPtVsV0MTrue_Generated[ih]);
        }
    }
    for(Int_t ih=0; ih<9; ih++){
        if(! fHistPtVsV0MTrue_DataSelection[ih] ) {
            fHistPtVsV0MTrue_DataSelection[ih] = new TH2F(Form("fHistPtVsV0MTrue_DataSelection_%s",lPartNames[ih].Data()),    "Generated;p_{T} (GeV/c)",lNPtBins,0,lMaxPt,lNNchBins,lLowNchBound,lHighNchBound);
            fListHist->Add(fHistPtVsV0MTrue_DataSelection[ih]);
        }
    }
    for(Int_t ih=0; ih<9; ih++){
        if(! fHistPtVsV0MTrue_MCSelection[ih] ) {
            fHistPtVsV0MTrue_MCSelection[ih] = new TH2F(Form("fHistPtVsV0MTrue_MCSelection_%s",lPartNames[ih].Data()),    "Generated;p_{T} (GeV/c)",lNPtBins,0,lMaxPt,lNNchBins,lLowNchBound,lHighNchBound);
            fListHist->Add(fHistPtVsV0MTrue_MCSelection[ih]);
        }
    }
    
    //// Ntrackl08
    //2-Dimensional Histograms
    for(Int_t ih=0; ih<9; ih++){
        if(! fHistPtVsTracklets08Cent_Generated[ih] ) {
            fHistPtVsTracklets08Cent_Generated[ih] = new TH2F(Form("fHistPtVsTracklets08Cent_Generated_%s",lPartNames[ih].Data()),    "Generated;p_{T} (GeV/c)",lNPtBins,0,lMaxPt,100,0,100);
            fListHist->Add(fHistPtVsTracklets08Cent_Generated[ih]);
        } 
        if(! fHistPtVsTracklets08Cent_DataSelection[ih] ) {
            fHistPtVsTracklets08Cent_DataSelection[ih] = new TH2F(Form("fHistPtVsTracklets08Cent_DataSelection_%s",lPartNames[ih].Data()),    "Generated;p_{T} (GeV/c)",lNPtBins,0,lMaxPt,100,0,100);
            fListHist->Add(fHistPtVsTracklets08Cent_DataSelection[ih]);
        }
        if(! fHistPtVsTracklets08Cent_MCSelection[ih] ) {
            fHistPtVsTracklets08Cent_MCSelection[ih] = new TH2F(Form("fHistPtVsTracklets08Cent_MCSelection_%s",lPartNames[ih].Data()),    "Generated;p_{T} (GeV/c)",lNPtBins,0,lMaxPt,100,0,100);
            fListHist->Add(fHistPtVsTracklets08Cent_MCSelection[ih]);
        }
       if(! fHistPtVsTracklets08_Generated[ih] ) {
            fHistPtVsTracklets08_Generated[ih] = new TH2F(Form("fHistPtVsTracklets08_Generated_%s",lPartNames[ih].Data()),    "Generated;p_{T} (GeV/c)",lNPtBins,0,lMaxPt,nTracklBins,0,lUpperLimit);
            fListHist->Add(fHistPtVsTracklets08_Generated[ih]);
        }
        if(! fHistPtVsTracklets08_DataSelection[ih] ) {
            fHistPtVsTracklets08_DataSelection[ih] = new TH2F(Form("fHistPtVsTracklets08_DataSelection_%s",lPartNames[ih].Data()),    "Generated;p_{T} (GeV/c)",lNPtBins,0,lMaxPt,nTracklBins,0,lUpperLimit);
            fListHist->Add(fHistPtVsTracklets08_DataSelection[ih]);
        }
        if(! fHistPtVsTracklets08_MCSelection[ih] ) {
            fHistPtVsTracklets08_MCSelection[ih] = new TH2F(Form("fHistPtVsTracklets08_MCSelection_%s",lPartNames[ih].Data()),    "Generated;p_{T} (GeV/c)",lNPtBins,0,lMaxPt,nTracklBins,0,lUpperLimit);
            fListHist->Add(fHistPtVsTracklets08_MCSelection[ih]);
        }
        if(! fHistPtVsTracklets08True_Generated[ih] ) {
            fHistPtVsTracklets08True_Generated[ih] = new TH2F(Form("fHistPtVsTracklets08True_Generated_%s",lPartNames[ih].Data()),    "Generated;p_{T} (GeV/c)",lNPtBins,0,lMaxPt,lNNchBins,lLowNchBound,lHighNchBound);
            fListHist->Add(fHistPtVsTracklets08True_Generated[ih]);
        }
        if(! fHistPtVsTracklets08True_DataSelection[ih] ) {
            fHistPtVsTracklets08True_DataSelection[ih] = new TH2F(Form("fHistPtVsTracklets08True_DataSelection_%s",lPartNames[ih].Data()),    "Generated;p_{T} (GeV/c)",lNPtBins,0,lMaxPt,lNNchBins,lLowNchBound,lHighNchBound);
            fListHist->Add(fHistPtVsTracklets08True_DataSelection[ih]);
        }
        if(! fHistPtVsTracklets08True_MCSelection[ih] ) {
            fHistPtVsTracklets08True_MCSelection[ih] = new TH2F(Form("fHistPtVsTracklets08True_MCSelection_%s",lPartNames[ih].Data()),    "Generated;p_{T} (GeV/c)",lNPtBins,0,lMaxPt,lNNchBins,lLowNchBound,lHighNchBound);
            fListHist->Add(fHistPtVsTracklets08True_MCSelection[ih]);
        }
    }
     //// Ntrackl0815
    //2-Dimensional Histograms
    for(Int_t ih=0; ih<9; ih++){
        if(! fHistPtVsTracklets0815Cent_Generated[ih] ) {
            fHistPtVsTracklets0815Cent_Generated[ih] = new TH2F(Form("fHistPtVsTracklets0815Cent_Generated_%s",lPartNames[ih].Data()),    "Generated;p_{T} (GeV/c)",lNPtBins,0,lMaxPt,100,0,100);
            fListHist->Add(fHistPtVsTracklets0815Cent_Generated[ih]);
        }
        if(! fHistPtVsTracklets0815Cent_DataSelection[ih] ) {
            fHistPtVsTracklets0815Cent_DataSelection[ih] = new TH2F(Form("fHistPtVsTracklets0815Cent_DataSelection_%s",lPartNames[ih].Data()),    "Generated;p_{T} (GeV/c)",lNPtBins,0,lMaxPt,100,0,100);
            fListHist->Add(fHistPtVsTracklets0815Cent_DataSelection[ih]);
        }
        if(! fHistPtVsTracklets0815Cent_MCSelection[ih] ) {
            fHistPtVsTracklets0815Cent_MCSelection[ih] = new TH2F(Form("fHistPtVsTracklets0815Cent_MCSelection_%s",lPartNames[ih].Data()),    "Generated;p_{T} (GeV/c)",lNPtBins,0,lMaxPt,100,0,100);
            fListHist->Add(fHistPtVsTracklets0815Cent_MCSelection[ih]);
        }
       if(! fHistPtVsTracklets0815_Generated[ih] ) {
            fHistPtVsTracklets0815_Generated[ih] = new TH2F(Form("fHistPtVsTracklets0815_Generated_%s",lPartNames[ih].Data()),    "Generated;p_{T} (GeV/c)",lNPtBins,0,lMaxPt,nTracklBins,0,lUpperLimit);
            fListHist->Add(fHistPtVsTracklets0815_Generated[ih]);
        }
        if(! fHistPtVsTracklets0815_DataSelection[ih] ) {
            fHistPtVsTracklets0815_DataSelection[ih] = new TH2F(Form("fHistPtVsTracklets0815_DataSelection_%s",lPartNames[ih].Data()),    "Generated;p_{T} (GeV/c)",lNPtBins,0,lMaxPt,nTracklBins,0,lUpperLimit);
            fListHist->Add(fHistPtVsTracklets0815_DataSelection[ih]);
        }
        if(! fHistPtVsTracklets0815_MCSelection[ih] ) {
            fHistPtVsTracklets0815_MCSelection[ih] = new TH2F(Form("fHistPtVsTracklets0815_MCSelection_%s",lPartNames[ih].Data()),    "Generated;p_{T} (GeV/c)",lNPtBins,0,lMaxPt,nTracklBins,0,lUpperLimit);
            fListHist->Add(fHistPtVsTracklets0815_MCSelection[ih]);
        }
        if(! fHistPtVsTracklets0815True_Generated[ih] ) {
            fHistPtVsTracklets0815True_Generated[ih] = new TH2F(Form("fHistPtVsTracklets0815True_Generated_%s",lPartNames[ih].Data()),    "Generated;p_{T} (GeV/c)",lNPtBins,0,lMaxPt,lNNchBins,lLowNchBound,lHighNchBound);
            fListHist->Add(fHistPtVsTracklets0815True_Generated[ih]);
        }
        if(! fHistPtVsTracklets0815True_DataSelection[ih] ) {
            fHistPtVsTracklets0815True_DataSelection[ih] = new TH2F(Form("fHistPtVsTracklets0815True_DataSelection_%s",lPartNames[ih].Data()),    "Generated;p_{T} (GeV/c)",lNPtBins,0,lMaxPt,lNNchBins,lLowNchBound,lHighNchBound);
            fListHist->Add(fHistPtVsTracklets0815True_DataSelection[ih]);
        }
        if(! fHistPtVsTracklets0815True_MCSelection[ih] ) {
            fHistPtVsTracklets0815True_MCSelection[ih] = new TH2F(Form("fHistPtVsTracklets0815True_MCSelection_%s",lPartNames[ih].Data()),    "Generated;p_{T} (GeV/c)",lNPtBins,0,lMaxPt,lNNchBins,lLowNchBound,lHighNchBound);
            fListHist->Add(fHistPtVsTracklets0815True_MCSelection[ih]);
        }
    }


    //List of Histograms: Normal
    PostData(1, fListHist);
    
}// end UserCreateOutputObjects


//________________________________________________________________________
void AliAnalysisTaskPPVsMultCrossCheckMC::UserExec(Option_t *)
{
    // Main loop
    // Called for each event
    
    AliESDEvent *lESDevent = 0x0;
    AliMCEvent  *lMCevent  = 0x0;
    AliStack    *lMCstack  = 0x0;
    
    
    // Connect to the InputEvent
    // After these lines, we should have an ESD/AOD event + the number of V0s in it.
    
    // Appropriate for ESD analysis!
    
    if ( !lPureMonteCarlo ){
        lESDevent = dynamic_cast<AliESDEvent*>( InputEvent() );
        if (!lESDevent) {
            AliWarning("AliESDevent not available, this should be in pure Monte Carlo mode... \n");
            return;
        }
    }
    
    lMCevent = MCEvent();
    if (!lMCevent) {
        Printf("ERROR: Could not retrieve MC event \n");
        cout << "Name of the file with pb :" <<  fInputHandler->GetTree()->GetCurrentFile()->GetName() << endl;
        return;
    }
    
    //-------------------------------------------------------------------
    //Code for the acquisition of the 'perfect' primary vertex position
    Float_t lVertexZMC = -100;
    
    if( fCheckVtxZMC ){
        TArrayF mcPrimaryVtx;
        AliGenEventHeader* mcHeader=lMCevent->GenEventHeader();
        if(!mcHeader){
            Printf("ERROR: Could not retrieve MC header \n");
            return;
        }
        mcHeader->PrimaryVertex(mcPrimaryVtx);
        lVertexZMC = mcPrimaryVtx.At(2);
    }
    if( !fCheckVtxZMC ){lVertexZMC = 0.0; } //Put it at the center if not asked to check
    //-------------------------------------------------------------------
    
    lMCstack = lMCevent->Stack();
    if (!lMCstack) {
        Printf("ERROR: Could not retrieve MC stack \n");
        cout << "Name of the file with pb :" <<  fInputHandler->GetTree()->GetCurrentFile()->GetName() << endl;
        return;
    }
    Double_t lMagneticField = -10;
    if (!lPureMonteCarlo){
        lMagneticField = lESDevent->GetMagneticField( );
    }
    //------------------------------------------------
    // Event Selection
    //------------------------------------------------
    
    Bool_t isSelected = 0;
    Bool_t lEvSel_Triggered = kFALSE;
    
    if( !lPureMonteCarlo ) {
        if(fkSelectTriggerByName) isSelected = AliPPVsMultUtils::IsSelectedTrigger (lESDevent, fTrigName);
        else isSelected = AliPPVsMultUtils::IsSelectedTrigger (lESDevent, fTrigType);
    }

    //Physics Selection
    lEvSel_Triggered = isSelected;
    
    //Override physics selection requirement if running on pure monte carlo
    if( lPureMonteCarlo ) lEvSel_Triggered = kTRUE;
    
    //------------------------------------------------
    // Primary Vertex Requirements Section:
    //  ---> pp: has vertex, |z|<10cm
    //------------------------------------------------
    
    //classical Proton-proton like selection
    Float_t lEvSel_VtxZ = -100;
    Bool_t lEvSel_IsNotPileup = kFALSE;
    Bool_t lEvSel_IsNotPileupInMultBins = kFALSE;

    if( !lPureMonteCarlo ){
        const AliESDVertex *lPrimaryBestESDVtx     = lESDevent->GetPrimaryVertex();
        const AliESDVertex *lPrimaryTrackingESDVtx = lESDevent->GetPrimaryVertexTracks();
        const AliESDVertex *lPrimarySPDVtx         = lESDevent->GetPrimaryVertexSPD();
    
        Double_t lBestPrimaryVtxPos[3]          = {-100.0, -100.0, -100.0};
        lPrimaryBestESDVtx->GetXYZ( lBestPrimaryVtxPos );
        lEvSel_VtxZ = lBestPrimaryVtxPos[2];
    
        if( !lESDevent->IsPileupFromSPD()           ) lEvSel_IsNotPileup           = kTRUE;
        if( !lESDevent->IsPileupFromSPDInMultBins() ) lEvSel_IsNotPileupInMultBins = kTRUE;
    }
    
    //------------------------------------------------
    // Multiplicity Information Acquistion
    //------------------------------------------------
    
    //Monte Carlo Level information !
    //--------- GENERATED NUMBER OF CHARGED PARTICLES
    // ---> Variable Definition
    
    Long_t lNchEta5   = 0;
    Long_t lNchEta8   = 0;
    Long_t lNchEta8to15   = 0;
    Long_t lNchEta10  = 0;
    Long_t lNchVZEROA = 0;
    Long_t lNchVZEROC = 0;
    Bool_t lEvSel_INELgtZEROStackPrimaries=kFALSE;
    
    //----- Loop on Stack ----------------------------------------------------------------
    for (Int_t iCurrentLabelStack = 0;  iCurrentLabelStack < (lMCstack->GetNtrack()); iCurrentLabelStack++)
    {   // This is the begining of the loop on tracks
        TParticle* particleOne = lMCstack->Particle(iCurrentLabelStack);
        if(!particleOne) continue;
        if(!particleOne->GetPDG()) continue;
        Double_t lThisCharge = particleOne->GetPDG()->Charge()/3.;
        if(TMath::Abs(lThisCharge)<0.001) continue;
        if(! (lMCstack->IsPhysicalPrimary(iCurrentLabelStack)) ) continue;
        
        //Double_t gpt = particleOne -> Pt();
        Double_t geta = particleOne -> Eta();
        
        if( TMath::Abs(geta) < 0.5 ) lNchEta5++;
        if( TMath::Abs(geta) < 0.8 ) lNchEta8++;
        if( (TMath::Abs(geta) > 0.8) && (TMath::Abs(geta) < 1.5) ) lNchEta8to15++;
        if( TMath::Abs(geta) < 1.0 ) lNchEta10++;
        if( TMath::Abs(geta) < 1.0 ) lEvSel_INELgtZEROStackPrimaries = kTRUE;
        if( 2.8 < geta && geta < 5.1 ) lNchVZEROA++;
        if(-3.7 < geta && geta <-1.7 ) lNchVZEROC++;
    }//End of loop on tracks
    //----- End Loop on Stack ------------------------------------------------------------
    
    //Standard GetReferenceMultiplicity Estimator (0.5 and 0.8)
    
    Float_t fCentrality_V0M   = -100;
    Float_t fCentrality_V0MUnselected = -100;
    Float_t fCentrality_Tracklets08   = -100;
    Float_t fCentrality_Tracklets08Unselected = -100;
    Float_t fCentrality_Tracklets0815   = -100;
    Float_t fCentrality_Tracklets0815Unselected = -100;
 
    if ( ! lPureMonteCarlo ) {
     if(!fkMultSelection){
     fCentrality_V0M   = fPPVsMultUtils -> GetMultiplicityPercentile(lESDevent, "V0M"   );
     fCentrality_V0MUnselected = fPPVsMultUtils -> GetMultiplicityPercentile(lESDevent, "V0M" , kFALSE );
     }else{
       AliMultSelection *MultSelection = (AliMultSelection*) lESDevent -> FindListObject("MultSelection");
        if( MultSelection ){
            fCentrality_V0M = MultSelection->GetMultiplicityPercentile("V0M",kTRUE);
            fCentrality_V0MUnselected = MultSelection->GetMultiplicityPercentile("V0M",kFALSE);
            fCentrality_Tracklets08 = MultSelection->GetMultiplicityPercentile("SPDTracklets08Corr",kTRUE);
            fCentrality_Tracklets08Unselected = MultSelection->GetMultiplicityPercentile("SPDTracklets08Corr",kFALSE);
            fCentrality_Tracklets0815 = MultSelection->GetMultiplicityPercentile("SPDTracklets08to15Corr",kTRUE);
            fCentrality_Tracklets0815Unselected = MultSelection->GetMultiplicityPercentile("SPDTracklets08to15Corr",kFALSE);
            }
          }
     } 

    //Amplitude if not pure
    Float_t fV0MAmplitude = -1;
    Int_t nSPDtrackl08 = -1;
    Int_t nSPDtrackl0815 = -1;
    if ( !lPureMonteCarlo ) 
    {
     fV0MAmplitude = GetV0MAmplitude( lESDevent );
     if(fkMultSelection){
       AliMultSelection *MultSelection = (AliMultSelection*) lESDevent -> FindListObject("MultSelection");
       nSPDtrackl08 = MultSelection->GetEstimator("SPDTracklets08")->GetValue(); // raw value not corrected for z-vtx
       nSPDtrackl0815 = MultSelection->GetEstimator("SPDTracklets08to15")->GetValue(); // raw value not corrected for z-vtx
     }
     
    }
    //------------------------------------------------
    // Get All Conditionals
    //------------------------------------------------
    Bool_t lIsINELgtZEROtracklets    = kFALSE;
    Bool_t lIsAcceptedVertexPosition = kFALSE;
    Bool_t lIsNotPileupInMultBins    = kFALSE;
    Bool_t lConsistentVertices       = kFALSE;
    Bool_t lIsNotSPDClusterVsTrackletBG = kFALSE; // (used only for 13 TeV analysis) set it as optional cut ? 


    if( !lPureMonteCarlo ){
        //1) Physics Selection
        //   (already exists, it's called "lEvSel_Triggered" here)
        //2) INEL>0 (data)
        lIsINELgtZEROtracklets    = AliPPVsMultUtils::IsINELgtZERO( lESDevent );
        //3) Is accepted vertex position
        lIsAcceptedVertexPosition = AliPPVsMultUtils::IsAcceptedVertexPosition( lESDevent );
        //4) IsNotPileupInMultBins
        lIsNotPileupInMultBins = AliPPVsMultUtils::IsNotPileupSPDInMultBins( lESDevent );
        //5) Consistent track / SPD vertices
        lConsistentVertices    = AliPPVsMultUtils::HasNoInconsistentSPDandTrackVertices( lESDevent );
        //6) Tracklets vs Cluster SPD
        lIsNotSPDClusterVsTrackletBG = ! ( fUtils->IsSPDClusterVsTrackletBG( lESDevent) );
    }
    
    //Monte-Carlo Specifics
    //Alpha: INELgtZERO in the MC stack
    // (already exists, called "fEvSel_INELgtZEROStackPrimaries")
    //Beta: Vertex position
    Bool_t lIsAcceptedVertexPositionMC = (TMath::Abs(lVertexZMC)<10.0); //true if within desired range

    if( fUseRecoVtxZ ) lIsAcceptedVertexPositionMC = lIsAcceptedVertexPosition; //override!
    
    //Merge all conditionals
    Bool_t lDataSelection = ( lEvSel_Triggered && lIsINELgtZEROtracklets && lIsAcceptedVertexPosition && lIsNotPileupInMultBins && lConsistentVertices && lIsNotSPDClusterVsTrackletBG );
    Bool_t lMCSelection   = ( ( fSkipPS ||lEvSel_Triggered ) && lEvSel_INELgtZEROStackPrimaries && lIsAcceptedVertexPositionMC );
    
    //Alternate Selection: Factor out only INEL>0 selection (extra cross-check) 
    if ( fAlternateMCSelection ){
        lMCSelection   = ( ( fSkipPS ||lEvSel_Triggered ) && lEvSel_INELgtZEROStackPrimaries && lIsAcceptedVertexPosition && lIsNotPileupInMultBins && lConsistentVertices );
    }
    
    //------------------------------------------------
    // Fill Event Counters
    //------------------------------------------------

    //Basics: All Processed
    fHistEventCounter->Fill(0.5);
    if( lDataSelection ) fHistEventCounter -> Fill(1.5);
    if( lMCSelection   ) fHistEventCounter -> Fill(2.5);
    
    if( lDataSelection ) {
        fHistV0M_DataSelection -> Fill( fCentrality_V0MUnselected );
        fHistV0MAmplitude_DataSelection -> Fill( fV0MAmplitude );
        fHistV0MTrue_DataSelection -> Fill( lNchVZEROA+lNchVZEROC ); 
        ////
        fHistTracklets08Cent_DataSelection -> Fill( fCentrality_Tracklets08Unselected );
        fHistTracklets08_DataSelection -> Fill( nSPDtrackl08);
        fHistTracklets08True_DataSelection -> Fill( lNchEta8 ); 
        ////
        fHistTracklets0815Cent_DataSelection -> Fill( fCentrality_Tracklets0815Unselected );
        fHistTracklets0815_DataSelection -> Fill( nSPDtrackl0815);
        fHistTracklets0815True_DataSelection -> Fill( lNchEta8to15 ); 
    }
    if( lMCSelection   ){
        fHistV0M_MCSelection   -> Fill( fCentrality_V0MUnselected );
        fHistV0MAmplitude_MCSelection -> Fill( fV0MAmplitude );
        fHistV0MTrue_MCSelection -> Fill( lNchVZEROA+lNchVZEROC );
        ////
        fHistTracklets08Cent_MCSelection -> Fill( fCentrality_Tracklets08Unselected );
        fHistTracklets08_MCSelection -> Fill( nSPDtrackl08);
        fHistTracklets08True_MCSelection -> Fill( lNchEta8 );
        ////
        fHistTracklets0815Cent_MCSelection -> Fill( fCentrality_Tracklets0815Unselected );
        fHistTracklets0815_MCSelection -> Fill( nSPDtrackl0815);
        fHistTracklets0815True_MCSelection -> Fill( lNchEta8to15 );
    }

    if( lDataSelection ) {
        fHistV0MVsMidRapidityTrue_DataSelection->Fill( fCentrality_V0MUnselected, lNchEta5 );
        fHistV0MTrueVsMidRapidityTrue_DataSelection->Fill( lNchVZEROA+lNchVZEROC, lNchEta5 );
        fHistV0MAmplitudeVsMidRapidityTrue_DataSelection->Fill( fV0MAmplitude, lNchEta5 );
        ////
        fHistTracklets08CentVsMidRapidityTrue_DataSelection->Fill( fCentrality_Tracklets08Unselected, lNchEta5 );
        fHistTracklets08TrueVsMidRapidityTrue_DataSelection->Fill( lNchEta8, lNchEta5 );
        fHistTracklets08VsMidRapidityTrue_DataSelection->Fill( nSPDtrackl08, lNchEta5 );
        ////
        fHistTracklets0815CentVsMidRapidityTrue_DataSelection->Fill( fCentrality_Tracklets0815Unselected, lNchEta5 );
        fHistTracklets0815TrueVsMidRapidityTrue_DataSelection->Fill( lNchEta8to15, lNchEta5 );
        fHistTracklets0815VsMidRapidityTrue_DataSelection->Fill( nSPDtrackl0815, lNchEta5 );
    }
    
    if( lMCSelection ) {
        fHistV0MVsMidRapidityTrue_MCSelection->Fill( fCentrality_V0MUnselected, lNchEta5 );
        fHistV0MTrueVsMidRapidityTrue_MCSelection->Fill( lNchVZEROA+lNchVZEROC, lNchEta5 );
        fHistV0MAmplitudeVsMidRapidityTrue_MCSelection->Fill( fV0MAmplitude, lNchEta5 );
        ////
        fHistTracklets08CentVsMidRapidityTrue_MCSelection->Fill( fCentrality_Tracklets08Unselected, lNchEta5 );
        fHistTracklets08TrueVsMidRapidityTrue_MCSelection->Fill( lNchEta8, lNchEta5 );
        fHistTracklets08VsMidRapidityTrue_MCSelection->Fill( nSPDtrackl08, lNchEta5 );
        ////
        fHistTracklets0815CentVsMidRapidityTrue_MCSelection->Fill( fCentrality_Tracklets0815Unselected, lNchEta5 );
        fHistTracklets0815TrueVsMidRapidityTrue_MCSelection->Fill( lNchEta8to15, lNchEta5 );
        fHistTracklets0815VsMidRapidityTrue_MCSelection->Fill( nSPDtrackl0815, lNchEta5 );
    }
    
    //------------------------------------------------
    // Fill Spectra as Needed
    //------------------------------------------------
    
    //~All relevant PWG-LF Identified Particle Information (for looping)
    Int_t lPDGCodes[9] = {211, 321, 2212, 310, 3122, 3312, 3334, 333, 313};
    TString lPartNames[9] = {
        "Pion", "Kaon", "Proton", "K0Short", "Lambda", "Xi", "Omega", "Phi", "KStar"
    };
    Bool_t lCheckIsPhysicalPrimary[9] = { kTRUE, kTRUE, kTRUE, kTRUE, kTRUE, kTRUE, kTRUE, kFALSE, kFALSE };
    
    Int_t lThisPDG  = 0;
    Double_t lThisRap  = 0;
    Double_t lThisPt   = 0;
    Bool_t lIsPhysicalPrimary = kFALSE;
    
    //----- Loop on Stack Starts Here ---------------
    for (Int_t ilab = 0;  ilab < (lMCstack->GetNtrack()); ilab++)
    {   // This is the begining of the loop on tracks
        
        TParticle* lPart = 0x0;
        lPart = lMCstack->Particle( ilab );
        if(!lPart) {
            Printf("Generated loop %d - MC TParticle pointer to current stack particle = 0x0 ! Skip ...\n", ilab );
            continue;
        }
        
        lThisPDG = lPart->GetPdgCode();
        
        //Continue if this is not a particle of the right PDG Code (avoids y-calculation problems)
        Bool_t lContinue = kTRUE;
        for(Int_t ih=0; ih<9; ih++) if( TMath::Abs(lThisPDG) == lPDGCodes[ih] ) lContinue = kFALSE;
        if ( lContinue ) continue;
            
        lThisRap   = MyRapidity(lPart->Energy(),lPart->Pz());
        lThisPt    = lPart->Pt();
        
        //Use Physical Primaries only for filling These Histos
        //if ( lMCstack->IsPhysicalPrimary(ilab)!=kTRUE ) continue;
        lIsPhysicalPrimary = lMCstack->IsPhysicalPrimary(ilab);
        
        for(Int_t ih=0; ih<9; ih++){
            if( TMath::Abs(lThisPDG) == lPDGCodes[ih] && TMath::Abs(lThisRap) < 0.5 ) {
                //Check if primary (if needed) and if not don't use this particle
                if( lCheckIsPhysicalPrimary[ih] == kTRUE && lIsPhysicalPrimary == kFALSE ) continue;
                //Fill Histograms
                fHistPt_Generated     [ih] -> Fill(lThisPt);
                fHistPtVsV0M_Generated[ih] -> Fill(lThisPt,fCentrality_V0MUnselected);
                fHistPtVsV0MAmplitude_Generated[ih] -> Fill(lThisPt,fV0MAmplitude);
                fHistPtVsV0MTrue_Generated     [ih] -> Fill(lThisPt,lNchVZEROA+lNchVZEROC);
                ///// 
                fHistPtVsTracklets08Cent_Generated[ih] -> Fill(lThisPt,fCentrality_Tracklets08Unselected);
                fHistPtVsTracklets08_Generated[ih] -> Fill(lThisPt, nSPDtrackl08);
                fHistPtVsTracklets08True_Generated     [ih] -> Fill(lThisPt, lNchEta8);
                ///// 
                fHistPtVsTracklets0815Cent_Generated[ih] -> Fill(lThisPt,fCentrality_Tracklets0815Unselected);
                fHistPtVsTracklets0815_Generated[ih] -> Fill(lThisPt, nSPDtrackl0815);
                fHistPtVsTracklets0815True_Generated     [ih] -> Fill(lThisPt, lNchEta8to15);
                if( lDataSelection ){
                    fHistPt_DataSelection     [ih] -> Fill(lThisPt);
                    fHistPtVsV0M_DataSelection[ih] -> Fill(lThisPt,fCentrality_V0MUnselected);
                    fHistPtVsV0MAmplitude_DataSelection[ih] -> Fill(lThisPt,fV0MAmplitude);
                    fHistPtVsV0MTrue_DataSelection     [ih] -> Fill(lThisPt,lNchVZEROA+lNchVZEROC);
                    ////
                    fHistPtVsTracklets08Cent_DataSelection[ih] -> Fill(lThisPt,fCentrality_Tracklets08Unselected);
                    fHistPtVsTracklets08_DataSelection[ih] -> Fill(lThisPt,nSPDtrackl08);
                    fHistPtVsTracklets08True_DataSelection     [ih] -> Fill(lThisPt,lNchEta8);
                    ////
                    fHistPtVsTracklets0815Cent_DataSelection[ih] -> Fill(lThisPt,fCentrality_Tracklets0815Unselected);
                    fHistPtVsTracklets0815_DataSelection[ih] -> Fill(lThisPt,nSPDtrackl0815);
                    fHistPtVsTracklets0815True_DataSelection     [ih] -> Fill(lThisPt,lNchEta8to15);
                }
                if( lMCSelection   ){
                    fHistPt_MCSelection       [ih] -> Fill(lThisPt);
                    fHistPtVsV0M_MCSelection  [ih] -> Fill(lThisPt,fCentrality_V0MUnselected);
                    fHistPtVsV0MAmplitude_MCSelection  [ih] -> Fill(lThisPt,fV0MAmplitude);
                    fHistPtVsV0MTrue_MCSelection       [ih] -> Fill(lThisPt,lNchVZEROA+lNchVZEROC);
                    ////
                    fHistPtVsTracklets08Cent_MCSelection[ih] -> Fill(lThisPt,fCentrality_Tracklets08Unselected);
                    fHistPtVsTracklets08_MCSelection[ih] -> Fill(lThisPt,nSPDtrackl08);
                    fHistPtVsTracklets08True_MCSelection     [ih] -> Fill(lThisPt,lNchEta8);
                    ////
                    fHistPtVsTracklets0815Cent_MCSelection[ih] -> Fill(lThisPt,fCentrality_Tracklets0815Unselected);
                    fHistPtVsTracklets0815_MCSelection[ih] -> Fill(lThisPt,nSPDtrackl0815);
                    fHistPtVsTracklets0815True_MCSelection     [ih] -> Fill(lThisPt,lNchEta8to15);
                }
            }
        }
    }//End of loop on tracks
    //----- End Loop on Stack ----------------------
    
    // Post output data.
    PostData(1, fListHist);
}

//________________________________________________________________________
void AliAnalysisTaskPPVsMultCrossCheckMC::Terminate(Option_t *)
{
    // Draw result to the screen
    // Called once at the end of the query
    
    TList *cRetrievedList = 0x0;
    cRetrievedList = (TList*)GetOutputData(1);
    if(!cRetrievedList) {
        Printf("ERROR - AliAnalysisTaskPPVsMultCrossCheckMC : ouput data container list not available\n");
        return;
    }
    
    fHistEventCounter = dynamic_cast<TH1F*> (  cRetrievedList->FindObject("fHistEventCounter")  );
    if (!fHistEventCounter) {
        Printf("ERROR - AliAnalysisTaskPPVsMultCrossCheckMC : fHistEventCounter not available");
        return;
    }
    
    TCanvas *canCheck = new TCanvas("AliAnalysisTaskPPVsMultCrossCheckMC","Event Multiplicity",10,10,510,510);
    canCheck->cd(1)->SetLogy();
    
    fHistEventCounter->SetMarkerStyle(22);
    fHistEventCounter->DrawCopy("E");
}

//----------------------------------------------------------------------------
Double_t AliAnalysisTaskPPVsMultCrossCheckMC::MyRapidity(Double_t rE, Double_t rPz) const
{
    // Local calculation for rapidity
    Double_t ReturnValue = -100;
    if( (rE-rPz+1.e-13) != 0 && (rE+rPz) != 0 ) {
        ReturnValue =  0.5*TMath::Log((rE+rPz)/(rE-rPz+1.e-13));
    }
    return ReturnValue;
}

//----------------------------------------------------------------------------
Double_t AliAnalysisTaskPPVsMultCrossCheckMC::GetV0MAmplitude( AliESDEvent *lInputEvent ) const
{
    // Get Amplitude from ESD
    //Get VZERO Information for multiplicity later
    AliVVZERO* esdV0 = lInputEvent->GetVZEROData();
    if (!esdV0) {
        AliError("AliVVZERO not available");
        return -1;
    }
    Double_t lReturnValue = esdV0->GetMTotV0A() + esdV0->GetMTotV0C();
    return lReturnValue;
}
