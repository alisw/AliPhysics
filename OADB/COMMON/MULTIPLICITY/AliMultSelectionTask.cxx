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
// This Task is meant to be used as a Multiplicity Determination task.
// It can also be used in calibration mode to acquire a TTree object with
// simple event characteristics, which can later be used for debugging,
// calibration, QA and calibration tests.
//
//  --- david.dobrigkeit.chinellato@cern.ch
//  --- alberica.toia@cern.ch
//  --- tatiana.drozhzhova@cern.ch
//
// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class TTree;
class TParticle;
class TVector3;

class AliESDVertex;
class AliAODVertex;
class AliESDAD; //AD


#include <AliVAD.h> //AD
#include <Riostream.h>
#include "TList.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"
#include "TFile.h"
#include "THnSparse.h"
#include "TVector3.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TLegend.h"
#include "TObjectTable.h"
#include "TSystem.h"
//#include "AliLog.h"

#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliAODHandler.h"
#include "AliESDpid.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliInputEventHandler.h"
#include "AliAnalysisManager.h"
#include "AliGenEventHeader.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliGenHijingEventHeader.h"
#include "AliGenDPMjetEventHeader.h"
#include "AliGenCocktailEventHeader.h"
#include "AliMCParticle.h"

#include "AliESDAD.h" //AD
#include "AliVZDC.h" //AD

//#include "AliCFContainer.h"
#include "AliMultiplicity.h"
#include "AliESDUtils.h"
#include "AliAnalysisUtils.h"
#include "AliAnalysisTaskSE.h"

//For MultSelection Framework
#include "AliOADBContainer.h"
#include "AliOADBMultSelection.h"
#include "AliMultEstimator.h"
#include "AliMultVariable.h"
#include "AliMultInput.h"
#include "AliMultSelection.h"
#include "AliMultSelectionCuts.h"

//task header
#include "AliMultSelectionTask.h"

using std::cout;
using std::endl;

ClassImp(AliMultSelectionTask)

AliMultSelectionTask::AliMultSelectionTask()
    : AliAnalysisTaskSE(), fListHist(0), fTreeEvent(0),fESDtrackCuts(0), fTrackCuts(0), fUtils(0),
      fkCalibration ( kFALSE ), fkAddInfo(kTRUE), fkFilterMB(kTRUE), fkAttached(0), fkDebug(kTRUE),
      fkDebugAliCentrality ( kFALSE ), fkDebugAliPPVsMultUtils( kFALSE ), fkDebugIsMC( kFALSE ),
      fkUseDefaultCalib (kFALSE), fkUseDefaultMCCalib (kFALSE),
      fkTrigger(AliVEvent::kINT7), fAlternateOADBForEstimators(""),
      fAlternateOADBFullManualBypass(""),fAlternateOADBFullManualBypassMC(""),
      fZncEnergy(0),
      fZpcEnergy(0),
      fZnaEnergy(0),
      fZpaEnergy(0),
      fZem1Energy(0),
      fZem2Energy(0),
      fZnaTower(0),
      fZncTower(0),
      fZpaTower(0),
      fZpcTower(0),
      fZnaFired(0),
      fZncFired(0),
      fZpaFired(0),
      fZpcFired(0),
      fNTracks(0),
      fCurrentRun(-1),
      fMultiplicity_ADA (0),
      fMultiplicity_ADC (0),
      fAmplitude_V0A   (0),
      fAmplitude_V0C   (0),
      fAmplitude_V0Apartial   (0),
      fAmplitude_V0Cpartial   (0),
      fAmplitude_V0AEq (0),
      fAmplitude_V0CEq (0),
      fAmplitude_OnlineV0A(0),
      fAmplitude_OnlineV0C(0),
      fAmplitude_V0A1(0),
      fAmplitude_V0A2(0),
      fAmplitude_V0A3(0),
      fAmplitude_V0A4(0),
      fAmplitude_V0C1(0),
      fAmplitude_V0C2(0),
      fAmplitude_V0C3(0),
      fAmplitude_V0C4(0),
      fAmplitude_V0AADC   (0),
      fAmplitude_V0CADC   (0),
      fnSPDClusters(0),
      fnTracklets(0),
      fnSPDClusters0(0),
      fnSPDClusters1(0),
      fnContributors(0),
      fRefMultEta5(0),
      fRefMultEta8(0),
      fRunNumber(0),
      fEvSel_VtxZCut(0),
      fEvSel_IsNotPileup(0),
      fEvSel_IsNotPileupMV(0),
      fEvSel_IsNotPileupInMultBins(0),
      fEvSel_Triggered(0),
      fEvSel_INELgtZERO(0),
      fEvSel_HasNoInconsistentVertices(0),
      fEvSel_PassesTrackletVsCluster(0),
      fEvSel_VtxZ(0),
      fEvSelCode(0),
      fNDebug(1),
      fAliCentralityV0M(0),
      fPPVsMultUtilsV0M(0),
      fMC_NColl(-1),
      fMC_NPart(-1),
      fMC_NchV0A(-1),
      fMC_NchV0C(-1),
      fMC_NchEta05(-1),
      fMC_NchEta08(-1),
      fMC_NchEta10(-1),
      
      //Histos
      fHistEventCounter(0),
      fHistQA_V0M(0),
      fHistQA_CL0(0),
      fHistQA_CL1(0),
      fHistQA_TrackletsVsV0M(0),
      fHistQA_TrackletsVsCL0(0),
      fHistQA_TrackletsVsCL1(0),
      fHistQASelected_V0M(0),
      fHistQASelected_CL0(0),
      fHistQASelected_CL1(0),
      fHistQASelected_TrackletsVsV0M(0),
      fHistQASelected_TrackletsVsCL0(0),
      fHistQASelected_TrackletsVsCL1(0),

      //Objects
      fOadbMultSelection(0),
      fInput(0)
//------------------------------------------------
// Tree Variables
{

}

AliMultSelectionTask::AliMultSelectionTask(const char *name, TString lExtraOptions, Bool_t lCalib, Int_t lNDebugEstimators)
    : AliAnalysisTaskSE(name), fListHist(0), fTreeEvent(0), fESDtrackCuts(0), fTrackCuts(0), fUtils(0),
      fkCalibration ( lCalib ), fkAddInfo(kTRUE), fkFilterMB(kTRUE), fkAttached(0), fkDebug(kTRUE),
      fkDebugAliCentrality ( kFALSE ), fkDebugAliPPVsMultUtils( kFALSE ), fkDebugIsMC ( kFALSE ),
      fkUseDefaultCalib (kFALSE), fkUseDefaultMCCalib (kFALSE),
      fkTrigger(AliVEvent::kINT7), fAlternateOADBForEstimators(""),
      fAlternateOADBFullManualBypass(""),fAlternateOADBFullManualBypassMC(""),
      fZncEnergy(0),
      fZpcEnergy(0),
      fZnaEnergy(0),
      fZpaEnergy(0),
      fZem1Energy(0),
      fZem2Energy(0),
      fZnaTower(0),
      fZncTower(0),
      fZpaTower(0),
      fZpcTower(0),
      fZnaFired(0),
      fZncFired(0),
      fZpaFired(0),
      fZpcFired(0),
      fNTracks(0),
      fCurrentRun(-1),
      fMultiplicity_ADA (0),
      fMultiplicity_ADC (0),
      fAmplitude_V0A   (0),
      fAmplitude_V0C   (0),
      fAmplitude_V0Apartial   (0),
      fAmplitude_V0Cpartial   (0),
      fAmplitude_V0AEq (0),
      fAmplitude_V0CEq (0),
      fAmplitude_OnlineV0A(0),
      fAmplitude_OnlineV0C(0),
      fAmplitude_V0A1(0),
      fAmplitude_V0A2(0),
      fAmplitude_V0A3(0),
      fAmplitude_V0A4(0),
      fAmplitude_V0C1(0),
      fAmplitude_V0C2(0),
      fAmplitude_V0C3(0),
      fAmplitude_V0C4(0),
      fAmplitude_V0AADC   (0),
      fAmplitude_V0CADC   (0),
      fnSPDClusters(0),
      fnTracklets(0),
      fnSPDClusters0(0),
      fnSPDClusters1(0),
      fnContributors(0),
      fRefMultEta5(0),
      fRefMultEta8(0),
      fRunNumber(0),
      fEvSel_VtxZCut(0),
      fEvSel_IsNotPileup(0),
      fEvSel_IsNotPileupMV(0),
      fEvSel_IsNotPileupInMultBins(0),
      fEvSel_Triggered(0),
      fEvSel_INELgtZERO(0),
      fEvSel_HasNoInconsistentVertices(0),
      fEvSel_PassesTrackletVsCluster(0),
      fEvSel_VtxZ(0),
      fEvSelCode(0),
      fNDebug(1),
      fAliCentralityV0M(0),
      fPPVsMultUtilsV0M(0),
      fMC_NColl(-1),
      fMC_NPart(-1),
      fMC_NchV0A(-1),
      fMC_NchV0C(-1),
      fMC_NchEta05(-1),
      fMC_NchEta08(-1),
      fMC_NchEta10(-1),
      
      //Histos
      fHistEventCounter(0),
      fHistQA_V0M(0),
      fHistQA_CL0(0),
      fHistQA_CL1(0),
      fHistQA_TrackletsVsV0M(0),
      fHistQA_TrackletsVsCL0(0),
      fHistQA_TrackletsVsCL1(0),
      fHistQASelected_V0M(0),
      fHistQASelected_CL0(0),
      fHistQASelected_CL1(0),
      fHistQASelected_TrackletsVsV0M(0),
      fHistQASelected_TrackletsVsCL0(0),
      fHistQASelected_TrackletsVsCL1(0),
      
      //Objects
      fOadbMultSelection(0),
      fInput(0)
{

    for( Int_t iq=0; iq<100; iq++ ) fQuantiles[iq] = -1 ;

    DefineOutput(1, TList::Class()); // Event Counter Histo
    if (fkCalibration) DefineOutput(2, TTree::Class()); // Event Tree

    //Create output slots for debugging estimators
    //Default: Save only first (should typically be V0M, depends on OADB!)
    fNDebug = lNDebugEstimators;

    //Special Debug Options (more to be added as needed)
    // A - Debug AliCentrality
    // B - Debug AliPPVsMultUtils
    // M - Extra MC variables

    if ( lExtraOptions.Contains("A") ) fkDebugAliCentrality = kTRUE;
    if ( lExtraOptions.Contains("B") ) fkDebugAliPPVsMultUtils = kTRUE;
    if ( lExtraOptions.Contains("M") ) fkDebugIsMC = kTRUE;
}


AliMultSelectionTask::~AliMultSelectionTask()
{
    //------------------------------------------------
    // DESTRUCTOR
    //------------------------------------------------//

    //if (fTreeEvent) {
    //    delete fTreeEvent;
    //    fTreeEvent = 0x0;
    //}
    if (fESDtrackCuts) {
        delete fESDtrackCuts;
        fESDtrackCuts = 0x0;
    }
    if(fTrackCuts) delete fTrackCuts;
    if (fUtils) {
        delete fUtils;
        fUtils = 0x0;
    }
}



//________________________________________________________________________
void AliMultSelectionTask::UserCreateOutputObjects()
{
    //------------------------------------------------

    //Create Input Information
    fInput = new AliMultInput();

    //Create input variables in AliMultInput Class
    //V0 related
    fAmplitude_V0A        = new AliMultVariable("fAmplitude_V0A");
    fAmplitude_V0A1       = new AliMultVariable("fAmplitude_V0A1");
    fAmplitude_V0A2       = new AliMultVariable("fAmplitude_V0A2");
    fAmplitude_V0A3       = new AliMultVariable("fAmplitude_V0A3");
    fAmplitude_V0A4       = new AliMultVariable("fAmplitude_V0A4");
    fAmplitude_V0C        = new AliMultVariable("fAmplitude_V0C");
    fAmplitude_V0C1       = new AliMultVariable("fAmplitude_V0C1");
    fAmplitude_V0C2       = new AliMultVariable("fAmplitude_V0C2");
    fAmplitude_V0C3       = new AliMultVariable("fAmplitude_V0C3");
    fAmplitude_V0C4       = new AliMultVariable("fAmplitude_V0C4");
    fAmplitude_V0Apartial = new AliMultVariable("fAmplitude_V0Apartial");
    fAmplitude_V0Cpartial = new AliMultVariable("fAmplitude_V0Cpartial");
    fAmplitude_V0AEq      = new AliMultVariable("fAmplitude_V0AEq");
    fAmplitude_V0CEq      = new AliMultVariable("fAmplitude_V0CEq");
    fAmplitude_OnlineV0A  = new AliMultVariable("fAmplitude_OnlineV0A");
    fAmplitude_OnlineV0C  = new AliMultVariable("fAmplitude_OnlineV0C");
    fAmplitude_V0AADC        = new AliMultVariable("fAmplitude_V0AADC");
    fAmplitude_V0CADC        = new AliMultVariable("fAmplitude_V0CADC");
    //SPD Related
    fnSPDClusters         = new AliMultVariable("fnSPDClusters");
    fnSPDClusters->SetIsInteger( kTRUE );
    fnSPDClusters0         = new AliMultVariable("fnSPDClusters0");
    fnSPDClusters0->SetIsInteger( kTRUE );
    fnSPDClusters1         = new AliMultVariable("fnSPDClusters1");
    fnSPDClusters1->SetIsInteger( kTRUE );

    //AD Related
    fMultiplicity_ADA     = new AliMultVariable("fMultiplicity_ADA");
    fMultiplicity_ADC     = new AliMultVariable("fMultiplicity_ADC");

    fnTracklets = new AliMultVariable("fnTracklets");
    fnTracklets ->SetIsInteger( kTRUE );
    fRefMultEta5 = new AliMultVariable("fRefMultEta5");
    fRefMultEta5 ->SetIsInteger( kTRUE );
    fRefMultEta8 = new AliMultVariable("fRefMultEta8");
    fRefMultEta8 ->SetIsInteger( kTRUE );

    //ZDC Related
    fZncEnergy = new AliMultVariable("fZncEnergy");
    fZpcEnergy = new AliMultVariable("fZpcEnergy");
    fZnaEnergy = new AliMultVariable("fZnaEnergy");
    fZpaEnergy = new AliMultVariable("fZpaEnergy");
    fZem1Energy = new AliMultVariable("fZem1Energy");
    fZem2Energy = new AliMultVariable("fZem2Energy");

    fZnaTower = new AliMultVariable("fZnaTower");
    fZncTower = new AliMultVariable("fZncTower");
    fZpaTower = new AliMultVariable("fZpaTower");
    fZpcTower = new AliMultVariable("fZpcTower");

    //Fired or not booleans (stored as integer for compatibility)
    fZnaFired = new AliMultVariable("fZnaFired"); 
    fZnaFired->SetIsInteger(kTRUE); 
    fZncFired = new AliMultVariable("fZncFired");
    fZncFired->SetIsInteger(kTRUE); 
    fZpaFired = new AliMultVariable("fZpaFired");
    fZpaFired->SetIsInteger(kTRUE);  
    fZpcFired = new AliMultVariable("fZpcFired"); 
    fZpcFired->SetIsInteger(kTRUE); 

    fEvSel_VtxZ = new AliMultVariable("fEvSel_VtxZ");

    //Add to AliMultInput Object, will later bind to TTree object in a loop
    fInput->AddVariable( fAmplitude_V0A );
    fInput->AddVariable( fAmplitude_V0A1 );
    fInput->AddVariable( fAmplitude_V0A2 );
    fInput->AddVariable( fAmplitude_V0A3 );
    fInput->AddVariable( fAmplitude_V0A4 );
    fInput->AddVariable( fAmplitude_V0C );
    fInput->AddVariable( fAmplitude_V0C1 );
    fInput->AddVariable( fAmplitude_V0C2 );
    fInput->AddVariable( fAmplitude_V0C3 );
    fInput->AddVariable( fAmplitude_V0C4 );
    fInput->AddVariable( fAmplitude_V0Apartial );
    fInput->AddVariable( fAmplitude_V0Cpartial );
    fInput->AddVariable( fAmplitude_V0AEq );
    fInput->AddVariable( fAmplitude_V0CEq );
    fInput->AddVariable( fAmplitude_OnlineV0A );
    fInput->AddVariable( fAmplitude_OnlineV0C );
    fInput->AddVariable( fAmplitude_V0AADC );
    fInput->AddVariable( fAmplitude_V0CADC );
    fInput->AddVariable( fnSPDClusters );
    fInput->AddVariable( fnSPDClusters0 );
    fInput->AddVariable( fnSPDClusters1 );
    fInput->AddVariable( fnTracklets );
    fInput->AddVariable( fRefMultEta5 );
    fInput->AddVariable( fRefMultEta8 );
    fInput->AddVariable( fMultiplicity_ADA );
    fInput->AddVariable( fMultiplicity_ADC );
    fInput->AddVariable( fZncEnergy );
    fInput->AddVariable( fZpcEnergy );
    fInput->AddVariable( fZnaEnergy );
    fInput->AddVariable( fZpaEnergy );
    fInput->AddVariable( fZem1Energy );
    fInput->AddVariable( fZem2Energy );
    fInput->AddVariable( fZnaTower );
    fInput->AddVariable( fZncTower );
    fInput->AddVariable( fZpaTower );
    fInput->AddVariable( fZpcTower );
    fInput->AddVariable( fZnaFired ); 
    fInput->AddVariable( fZncFired ); 
    fInput->AddVariable( fZpaFired ); 
    fInput->AddVariable( fZpcFired ); 
    fInput->AddVariable( fEvSel_VtxZ );

    if( fkCalibration ) {
        fTreeEvent = new TTree("fTreeEvent","Event");

        //------------------------------------------------
        // fTree Branch definitions - Event by Event info
        //------------------------------------------------

        //-----------BASIC-INFO---------------------------
        //Run Number
        fTreeEvent->Branch("fRunNumber", &fRunNumber, "fRunNumber/I");

        //Booleans for Event Selection
        fTreeEvent->Branch("fEvSel_VtxZCut", &fEvSel_VtxZCut, "fEvSel_VtxZCut/O");
        fTreeEvent->Branch("fEvSel_IsNotPileup", &fEvSel_IsNotPileup, "fEvSel_IsNotPileup/O");
        fTreeEvent->Branch("fEvSel_IsNotPileupMV", &fEvSel_IsNotPileupMV, "fEvSel_IsNotPileupMV/O");
        fTreeEvent->Branch("fEvSel_IsNotPileupInMultBins", &fEvSel_IsNotPileupInMultBins, "fEvSel_IsNotPileupInMultBins/O");
        fTreeEvent->Branch("fEvSel_Triggered", &fEvSel_Triggered, "fEvSel_Triggered/O");
        fTreeEvent->Branch("fEvSel_INELgtZERO", &fEvSel_INELgtZERO, "fEvSel_INELgtZERO/O");
        fTreeEvent->Branch("fEvSel_HasNoInconsistentVertices", &fEvSel_HasNoInconsistentVertices, "fEvSel_HasNoInconsistentVertices/O");
        fTreeEvent->Branch("fEvSel_PassesTrackletVsCluster", &fEvSel_PassesTrackletVsCluster, "fEvSel_PassesTrackletVsCluster/O");

        //A.T. FIXME change into AliMultVariable
        fTreeEvent->Branch("fnContributors", &fnContributors, "fnContributors/I");

        fTreeEvent->Branch("fNTracks",      &fNTracks, "fNTracks/I");

        if( fkDebugIsMC ) {
            fTreeEvent->Branch("fMC_NPart",      &fMC_NPart, "fMC_NPart/I");
            fTreeEvent->Branch("fMC_NColl",      &fMC_NColl, "fMC_NColl/I");
            fTreeEvent->Branch("fMC_NchV0A",      &fMC_NchV0A, "fMC_NchV0A/I");
            fTreeEvent->Branch("fMC_NchV0C",      &fMC_NchV0C, "fMC_NchV0C/I");
            fTreeEvent->Branch("fMC_NchEta05",      &fMC_NchEta05, "fMC_NchEta05/I");
            fTreeEvent->Branch("fMC_NchEta08",      &fMC_NchEta08, "fMC_NchEta08/I");
            fTreeEvent->Branch("fMC_NchEta10",      &fMC_NchEta10, "fMC_NchEta10/I");
        }

        //Automatic Loop for linking directly to AliMultInput
        for( Long_t iVar=0; iVar<fInput->GetNVariables(); iVar++) {
            if( !fInput->GetVariable(iVar)->IsInteger()  ) {
                fTreeEvent->Branch(fInput->GetVariable(iVar)->GetName(), &fInput->GetVariable(iVar)->GetRValue(), Form("%s/F",fInput->GetVariable(iVar)->GetName()) );
            } else {
                fTreeEvent->Branch(fInput->GetVariable(iVar)->GetName(), &fInput->GetVariable(iVar)->GetRValueInteger(), Form("%s/I",fInput->GetVariable(iVar)->GetName()) );
            }
        }

        if( fkDebug ) {
            fTreeEvent->Branch("fEvSelCode",      &fEvSelCode, "fEvSelCode/I");
            //Fixme: Save first 5 quantiles, should be enough for debugging
            for ( Int_t iq=0; iq<fNDebug; iq++) {
                fTreeEvent->Branch(Form("fDebug_Percentile_%i",iq), &fQuantiles[iq], Form("fDebug_Percentile_%i/F",iq));
            }
        }
        //Debug functionality
        if ( fkDebugAliCentrality ) {
            fTreeEvent->Branch("fAliCentralityV0M",&fAliCentralityV0M,"fAliCentralityV0M/F");
        }
        if ( fkDebugAliPPVsMultUtils ) {
            fTreeEvent->Branch("fPPVsMultUtilsV0M",&fPPVsMultUtilsV0M,"fPPVsMultUtilsV0M/F");
        }
    }
    //------------------------------------------------
    // Set up objects
    //------------------------------------------------

    // Multiplicity
    if(! fESDtrackCuts ) {
        fESDtrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE,kFALSE);
        fESDtrackCuts->SetPtRange(0.15);  // adding pt cut
        fESDtrackCuts->SetEtaRange(-1.0, 1.0);
    }
    
    //Create TPC only track cuts
    if(!fTrackCuts) fTrackCuts = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
    
    if(! fUtils ) {
        fUtils = new AliAnalysisUtils();
    }

    AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
    inputHandler->SetNeedField();

    //------------------------------------------------
    // Histograms
    //------------------------------------------------
    // Create histograms
    OpenFile(1);
    fListHist = new TList();
    fListHist->SetOwner();  // See http://root.cern.ch/root/html/TCollection.html#TCollection:SetOwner

    if(! fHistEventCounter ) {
        //Histogram Output: Event-by-Event
        fHistEventCounter = new TH1D( "fHistEventCounter", ";Evt. Sel. Step;Count",5,0,5);
        fHistEventCounter->GetXaxis()->SetBinLabel(1, "Processed");
        fHistEventCounter->GetXaxis()->SetBinLabel(2, "Phys-Sel");
        fHistEventCounter->GetXaxis()->SetBinLabel(3, "Has Vtx");
        fHistEventCounter->GetXaxis()->SetBinLabel(4, "Vtx |z|<10cm");
        fHistEventCounter->GetXaxis()->SetBinLabel(5, "Isn't Pileup (in mult bins)");
        fListHist->Add(fHistEventCounter);
    }

    //QA Histograms - User-side output for cross-checking
    if ( !fHistQA_V0M ) {
        fHistQA_V0M = new TH1D("fHistQA_V0M", ";V0M Percentile;Count", 105,0,105);
        fListHist->Add(fHistQA_V0M);
    }
    if ( !fHistQA_CL0 ) {
        fHistQA_CL0 = new TH1D("fHistQA_CL0", ";CL0 Percentile;Count", 105,0,105);
        fListHist->Add(fHistQA_CL0);
    }
    if ( !fHistQA_CL1 ) {
        fHistQA_CL1 = new TH1D("fHistQA_CL1", ";CL1 Percentile;Count", 105,0,105);
        fListHist->Add(fHistQA_CL1);
    }
    //To compare SPD tracklets in data and MC
    if ( !fHistQA_TrackletsVsV0M ) {
        fHistQA_TrackletsVsV0M = new TProfile("fHistQA_TrackletsVsV0M", ";V0M Percentile;#LTSPD Tracklets#GT", 105,0,105);
        fListHist->Add(fHistQA_TrackletsVsV0M);
    }
    if ( !fHistQA_TrackletsVsCL0 ) {
        fHistQA_TrackletsVsCL0 = new TProfile("fHistQA_TrackletsVsCL0", ";CL0 Percentile;#LTSPD Tracklets#GT", 105,0,105);
        fListHist->Add(fHistQA_TrackletsVsCL0);
    }
    if ( !fHistQA_TrackletsVsCL1 ) {
        fHistQA_TrackletsVsCL1 = new TProfile("fHistQA_TrackletsVsCL1", ";CL1 Percentile;#LTSPD Tracklets#GT", 105,0,105);
        fListHist->Add(fHistQA_TrackletsVsCL1);
    }
    
    //Histograms filled with event selection embedded 
    if ( !fHistQASelected_V0M ) {
        fHistQASelected_V0M = new TH1D("fHistQASelected_V0M", ";V0M Percentile;Count", 105,0,105);
        fListHist->Add(fHistQASelected_V0M);
    }
    if ( !fHistQASelected_CL0 ) {
        fHistQASelected_CL0 = new TH1D("fHistQASelected_CL0", ";CL0 Percentile;Count", 105,0,105);
        fListHist->Add(fHistQASelected_CL0);
    }
    if ( !fHistQASelected_CL1 ) {
        fHistQASelected_CL1 = new TH1D("fHistQASelected_CL1", ";CL1 Percentile;Count", 105,0,105);
        fListHist->Add(fHistQASelected_CL1);
    }
    //To compare SPD tracklets in data and MC
    if ( !fHistQASelected_TrackletsVsV0M ) {
        fHistQASelected_TrackletsVsV0M = new TProfile("fHistQASelected_TrackletsVsV0M", ";V0M Percentile;#LTSPD Tracklets#GT", 105,0,105);
        fListHist->Add(fHistQASelected_TrackletsVsV0M);
    }
    if ( !fHistQASelected_TrackletsVsCL0 ) {
        fHistQASelected_TrackletsVsCL0 = new TProfile("fHistQASelected_TrackletsVsCL0", ";CL0 Percentile;#LTSPD Tracklets#GT", 105,0,105);
        fListHist->Add(fHistQASelected_TrackletsVsCL0);
    }
    if ( !fHistQASelected_TrackletsVsCL1 ) {
        fHistQASelected_TrackletsVsCL1 = new TProfile("fHistQASelected_TrackletsVsCL1", ";CL1 Percentile;#LTSPD Tracklets#GT", 105,0,105);
        fListHist->Add(fHistQASelected_TrackletsVsCL1);
    }

    //List of Histograms: Normal
    PostData(1, fListHist);

    //TTree Object: Saved to base directory. Should cache to disk while saving.
    //(Important to avoid excessive memory usage, particularly when merging)
    if ( fkCalibration ) PostData(2, fTreeEvent);
}// end UserCreateOutputObjects


//________________________________________________________________________
void AliMultSelectionTask::UserExec(Option_t *)
{
    // Main loop
    // Called for each event

    Bool_t lVerbose = kFALSE ;

    //Debugging / Memory usage tests
    //gObjectTable->Print();

    //Modifications for running on AODs or ESDs
    AliVEvent *lVevent = 0x0;

    //Zero all booleans, etc: safe initialization per event
    fEvSel_VtxZCut                = kFALSE;
    fEvSel_IsNotPileup            = kFALSE;
    fEvSel_IsNotPileupMV          = kFALSE;
    fEvSel_IsNotPileupInMultBins  = kFALSE;
    fEvSel_Triggered              = kFALSE;
    fEvSel_PassesTrackletVsCluster   = kFALSE;
    fEvSel_HasNoInconsistentVertices = kFALSE;
    fEvSel_INELgtZERO             = kFALSE;
    //fnSPDClusters = -1;
    fnSPDClusters -> SetValueInteger(-1);
    fnSPDClusters0 -> SetValueInteger( -1) ;
    fnSPDClusters1 -> SetValueInteger( -1) ;
    fEvSel_VtxZ ->SetValue( -100 );

    Float_t multADA =0;
    Float_t multADC =0;
    Float_t multAD =0;

    fMC_NchV0A = -1;
    fMC_NchV0C = -1;
    fMC_NchEta05 = -1;
    fMC_NchEta08 = -1;
    fMC_NchEta10 = -1;

    // Connect to the InputEvent
    // Appropriate for ESD analysis ..

    if(lVerbose) Printf("Casting AliVEvent...");

    lVevent = dynamic_cast<AliVEvent*>( InputEvent() );
    if (!lVevent) {
        AliWarning("ERROR: ESD / AOD event not available \n");
        return;
    }
    if(lVerbose) Printf("Casting AliVVZERO...");

    //Get VZERO Information for multiplicity later
    AliVVZERO* lVV0 = lVevent->GetVZEROData();
    if (!lVV0) {
        AliError("AliVVZERO not available");
        return;
    }
    if(lVerbose) Printf("Casting AliVAD...");
    //Get AD Multiplicity Information
    AliVAD *lVAD = lVevent->GetADData();
    if(!lVAD) {
        //commented out for smaller logs!
        //AliWarning("ERROR:lVAD not available\n");
    }

    //Not Acquired!
    fAliCentralityV0M = -1;
    //if requested, grab AliCentrality value for this event
    if( fkDebugAliCentrality ) {
        AliCentrality* centrality;
        centrality = lVevent->GetCentrality();
        if ( centrality ) {
            fAliCentralityV0M = centrality->GetCentralityPercentile( "V0M" );
        }
    }

    //Not Acquired!
    fPPVsMultUtilsV0M = -1;
    //if requested, grab AliCentrality value for this event
    if( fkDebugAliPPVsMultUtils ) {
        fPPVsMultUtilsV0M = fUtils->GetMultiplicityPercentile( lVevent, "V0M" );
    }

    //------------------------------------------------
    //Information from MC (thanks to Alberica)
    //Don't forget to set: some of the "ifs" may not be there
    fMC_NColl = -1;
    fMC_NPart = -1;

    if ( fkDebugIsMC ) {
        AliAnalysisManager* anMan = AliAnalysisManager::GetAnalysisManager();
        AliMCEventHandler* eventHandler = (AliMCEventHandler*)anMan->GetMCtruthEventHandler();
        AliStack*    stack=0;
        AliMCEvent*  mcEvent=0;


        if (eventHandler && (mcEvent=eventHandler->MCEvent()) && (stack=mcEvent->Stack())) {

            //Npart and Ncoll information
            AliGenHijingEventHeader* hHijing=0;
            AliGenDPMjetEventHeader* dpmHeader=0;
            AliGenEventHeader* mcGenH = mcEvent->GenEventHeader();
            if (mcGenH->InheritsFrom(AliGenHijingEventHeader::Class()))
                hHijing = (AliGenHijingEventHeader*)mcGenH;
            else if (mcGenH->InheritsFrom(AliGenCocktailEventHeader::Class())) {
                TList* headers = ((AliGenCocktailEventHeader*)mcGenH)->GetHeaders();
                hHijing = dynamic_cast<AliGenHijingEventHeader*>(headers->FindObject("Hijing"));
                if (!hHijing) hHijing = dynamic_cast<AliGenHijingEventHeader*>(headers->FindObject("Hijing pPb_0"));
                if (!hHijing) hHijing = dynamic_cast<AliGenHijingEventHeader*>(headers->FindObject("Hijing_0"));
            }
            else if (mcGenH->InheritsFrom(AliGenDPMjetEventHeader::Class())) {
                dpmHeader = (AliGenDPMjetEventHeader*)mcGenH;
            }
            if(hHijing)   {
                fMC_NPart = hHijing->ProjectileParticipants()+hHijing->TargetParticipants();
                fMC_NColl = hHijing->NN()+hHijing->NNw()+hHijing->NwN()+hHijing->NwNw();
            }
            if(dpmHeader) {
                fMC_NPart =dpmHeader->ProjectileParticipants()+dpmHeader->TargetParticipants();
                fMC_NColl =dpmHeader->NN()+dpmHeader->NNw()+dpmHeader->NwN()+dpmHeader->NwNw();
            }

            //Nch information in V0A and V0C acceptance
            //Initialize counters to valid!
            fMC_NchV0A = 0;
            fMC_NchV0C = 0;
            fMC_NchEta05 = 0;
            fMC_NchEta08 = 0;
            fMC_NchEta10 = 0;
            //----- Loop on Stack ----------------------------------------------------------------
            for (Int_t iCurrentLabelStack = 0;  iCurrentLabelStack < (stack->GetNtrack()); iCurrentLabelStack++)
            {   // This is the begining of the loop on tracks
                TParticle* particleOne = stack->Particle(iCurrentLabelStack);
                if(!particleOne) continue;
                if(!particleOne->GetPDG()) continue;
                Double_t lThisCharge = particleOne->GetPDG()->Charge()/3.;
                if(TMath::Abs(lThisCharge)<0.001) continue;
                if(! (stack->IsPhysicalPrimary(iCurrentLabelStack)) ) continue;

                //Double_t gpt = particleOne -> Pt();
                Double_t geta = particleOne -> Eta();

                if( 2.8 < geta && geta < 5.1 ) fMC_NchV0A++;
                if(-3.7 < geta && geta <-1.7 ) fMC_NchV0C++;
                if(TMath::Abs( geta ) < 1.0 ) fMC_NchEta10++;
                if(TMath::Abs( geta ) < 0.8 ) fMC_NchEta08++;
                if(TMath::Abs( geta ) < 0.5 ) fMC_NchEta05++;
            }//End of loop on tracks
            //----- End Loop on Stack ------------------------------------------------------------
        }
    }
    //------------------------------------------------

    if(lVerbose) Printf("Starting...");
    fRunNumber = lVevent->GetRunNumber();
    Double_t lMagneticField = -10;
    lMagneticField = lVevent->GetMagneticField( );

    //------------------------------------------------
    // Physics Selection
    //------------------------------------------------

    fHistEventCounter->Fill(0.5);

    //===============================================
    // Event Selection Variables (fEvSel_xxx)
    //===============================================

    //------------------------------------------------
    // Done Via one-line functions
    // (static if possible)
    //------------------------------------------------
    if(lVerbose) Printf("Doing Event Selections...");
    fEvSel_Triggered                 = IsSelectedTrigger                   (lVevent, fkTrigger);
    fEvSel_IsNotPileup               = IsNotPileupSPD                      (lVevent);
    fEvSel_IsNotPileupInMultBins     = IsNotPileupSPDInMultBins            (lVevent);
    fEvSel_IsNotPileupMV             = IsNotPileupMV                       (lVevent);
    fEvSel_PassesTrackletVsCluster   = PassesTrackletVsCluster             (lVevent);
    fEvSel_HasNoInconsistentVertices = HasNoInconsistentSPDandTrackVertices(lVevent);
    fEvSel_INELgtZERO                = IsINELgtZERO                        (lVevent);

    //classical Proton-proton like selection
    const AliVVertex *lPrimaryBestESDVtx     = lVevent->GetPrimaryVertex();
    const AliVVertex *lPrimarySPDVtx         = lVevent->GetPrimaryVertexSPD();

    Double_t lBestPrimaryVtxPos[3]          = {-100.0, -100.0, -100.0};
    lPrimaryBestESDVtx->GetXYZ( lBestPrimaryVtxPos );

    //Number of contributors from best primary vertex
    fnContributors = -1;
    fnContributors = lPrimaryBestESDVtx -> GetNContributors();

    if(TMath::Abs(lBestPrimaryVtxPos[2]) <= 10.0 ) {
        //FIXME Passed default 10.0cm selection!
        fEvSel_VtxZCut = kTRUE;
    }
    fEvSel_VtxZ -> SetValue( lBestPrimaryVtxPos[2] ); //Set for later use

    //===============================================
    // End Event Selection Variables Section
    //===============================================
    if(lVerbose) Printf("Doing Multiplicity Calculations...");
    //------------------------------------------------
    // Multiplicity Information from AD
    //------------------------------------------------
    Float_t fMultiplicityAD[16];
    Float_t fMultiplicityADA[8];
    Float_t fMultiplicityADC[8];
    multADA =0;
    multADC =0;
    multAD =0;

    if (lVAD) {
        //Get Multiplicity info per AD 16 channel: C-side : 0-7, A-side 8-15
        for (Int_t i=0; i<8; i++)
        {
            fMultiplicityAD[i]= lVAD->GetMultiplicity(i);
            fMultiplicityADA[i]= lVAD->GetMultiplicityADA(i);
            multADA += fMultiplicityADA[i];
            multAD += fMultiplicityAD[i];
        }
        for (Int_t i=8; i<16; i++)
        {   fMultiplicityAD[i]= lVAD->GetMultiplicity(i);
            fMultiplicityADC[i]=lVAD->GetMultiplicityADC(i-8);
            multADC += fMultiplicityADC[i];
            multAD+=fMultiplicityAD[i];
        }
    }
    //------------------------------------------------
    // Multiplicity Information Acquistion
    //------------------------------------------------

    // VZERO PART
    Float_t  multV0A  = 0;            //  multiplicity from V0 reco side A
    Float_t  multV0C  = 0;            //  multiplicity from V0 reco side C
    Float_t  multV0AEq  = 0;          //  multiplicity from V0 reco side A
    Float_t  multV0CEq  = 0;          //  multiplicity from V0 reco side C
    Float_t  multV0ACorr  = 0;            //  multiplicity from V0 reco side A
    Float_t  multV0CCorr  = 0;            //  multiplicity from V0 reco side C
    Float_t multonlineV0A = 0; //charge
    Float_t multonlineV0C = 0; //charge
    //Special Ring selection (for eta symmetry)

    //Selection we want to perform here:
    // --- V0A Rings 3+4
    // --- V0C Rings 1+2

    Float_t multV0Apartial = 0;
    Float_t multV0Cpartial = 0;
    for(Int_t iCh = 48; iCh < 64; iCh++) {
        Double_t mult = lVV0->GetMultiplicity(iCh);
        multV0Apartial += mult;
    }
    for(Int_t iCh = 0; iCh < 16; iCh++) {
        Double_t mult = lVV0->GetMultiplicity(iCh);
        multV0Cpartial += mult;
    }

    Float_t multV0A1 = 0;
    Float_t multV0A2 = 0;
    Float_t multV0A3 = 0;
    Float_t multV0A4 = 0;
    Float_t multV0C1 = 0;
    Float_t multV0C2 = 0;
    Float_t multV0C3 = 0;
    Float_t multV0C4 = 0;
    for(Int_t iCh = 0; iCh < 8; iCh++) {
        Double_t mult = lVV0->GetMultiplicity(iCh);
        multV0C1 += mult;
    }
    for(Int_t iCh = 8; iCh < 16; iCh++) {
        Double_t mult = lVV0->GetMultiplicity(iCh);
        multV0C2 += mult;
    }
    for(Int_t iCh = 16; iCh < 24; iCh++) {
        Double_t mult = lVV0->GetMultiplicity(iCh);
        multV0C3 += mult;
    }
    for(Int_t iCh = 24; iCh < 32; iCh++) {
        Double_t mult = lVV0->GetMultiplicity(iCh);
        multV0C4 += mult;
    }

    for(Int_t iCh = 32; iCh < 40; iCh++) {
        Double_t mult = lVV0->GetMultiplicity(iCh);
        multV0A1 += mult;
    }
    for(Int_t iCh = 40; iCh < 48; iCh++) {
        Double_t mult = lVV0->GetMultiplicity(iCh);
        multV0A2 += mult;
    }
    for(Int_t iCh = 48; iCh < 56; iCh++) {
        Double_t mult = lVV0->GetMultiplicity(iCh);
        multV0A3 += mult;
    }
    for(Int_t iCh = 56; iCh < 64; iCh++) {
        Double_t mult = lVV0->GetMultiplicity(iCh);
        multV0A4 += mult;
    }

    //Non-Equalized Signal: copy of multV0ACorr and multV0CCorr from AliCentralitySelectionTask
    //Getters for uncorrected multiplicity
    multV0A=lVV0->GetMTotV0A();
    multV0C=lVV0->GetMTotV0C();
    //charge V0
    multonlineV0A = lVV0->GetTriggerChargeA(); //charge
    multonlineV0C = lVV0->GetTriggerChargeC(); //charge

    //Get Z vertex position of SPD vertex (why not Tracking if available?)
    Float_t zvtx = lPrimarySPDVtx->GetZ();

    //Acquire Corrected multV0A
    multV0ACorr = AliESDUtils::GetCorrV0A(multV0A,zvtx);
    multV0CCorr = AliESDUtils::GetCorrV0C(multV0C,zvtx);

    //Set Desired Variables

    fAmplitude_V0A->SetValue(multV0A);
    fAmplitude_V0C->SetValue(multV0C);

    //Implementation of V0 ADC information
    // FIXME: THIS ONLY WORKS IN ESDS FOR NOW
    Float_t  multV0AADC  = 0;            //  multiplicity from V0 reco side A from ADC
    Float_t  multV0CADC  = 0;            //  multiplicity from V0 reco side C from ADC
    fAmplitude_V0AADC->SetValue(0);
    fAmplitude_V0CADC->SetValue(0);

    /* FIXME: THIS DOES NOT WORK !!
     for(Int_t iCh = 0; iCh < 32; iCh++){
     Double_t mult = lVV0->GetADC(iCh);
     multV0CADC += mult;
     }
     for(Int_t iCh = 32; iCh < 64; iCh++){
     Double_t mult = lVV0->GetMultiplicity(iCh);
     multV0AADC += mult;
     }
     */

    fAmplitude_V0AADC->SetValue(multV0AADC);
    fAmplitude_V0CADC->SetValue(multV0CADC);

    if ( lVerbose ) {
        Printf(" V0A Amplitude: %.5f", multV0A );
        Printf(" V0C Amplitude: %.5f", multV0C );
    }

    fAmplitude_OnlineV0A->SetValue(multonlineV0A);
    fAmplitude_OnlineV0C->SetValue(multonlineV0C);

    fAmplitude_V0Apartial->SetValue(multV0Apartial);
    fAmplitude_V0Cpartial->SetValue(multV0Cpartial);

    //A.T. (vertex correction for all rings?!?)
    fAmplitude_V0A1 -> SetValue( multV0A1 );
    fAmplitude_V0A2 -> SetValue( multV0A2 );
    fAmplitude_V0A3 -> SetValue( multV0A3 );
    fAmplitude_V0A4 -> SetValue( multV0A4 );
    fAmplitude_V0C1 -> SetValue( multV0C1 );
    fAmplitude_V0C2 -> SetValue( multV0C2 );
    fAmplitude_V0C3 -> SetValue( multV0C3 );
    fAmplitude_V0C4 -> SetValue( multV0C4 );

    //AD scintillator Data added to Event Tree
    if (lVAD) {
        fMultiplicity_ADA->SetValue(multADA);
        fMultiplicity_ADC->SetValue(multADC);
    }

    // Equalized signals // From AliCentralitySelectionTask // Updated
    for(Int_t iCh = 32; iCh < 64; ++iCh) {
        Double_t mult = lVevent->GetVZEROEqMultiplicity(iCh);
        multV0AEq += mult;
    }
    for(Int_t iCh = 0; iCh < 32; ++iCh) {
        Double_t mult = lVevent->GetVZEROEqMultiplicity(iCh);
        multV0CEq += mult;
    }
    fAmplitude_V0AEq->SetValue(multV0AEq);
    fAmplitude_V0CEq->SetValue(multV0CEq);

    //Integer Estimators
    fnTracklets->SetValueInteger(lVevent->GetMultiplicity()->GetNumberOfTracklets());
    fnSPDClusters->SetValueInteger(lVevent->GetNumberOfITSClusters(0) + lVevent->GetNumberOfITSClusters(1));
    fnSPDClusters0 -> SetValueInteger(lVevent->GetNumberOfITSClusters(0));
    fnSPDClusters1 -> SetValueInteger(lVevent->GetNumberOfITSClusters(1));
    //===============================================
    //This part requires separation of AOD and ESD
    //===============================================

    //Setting variables to non-sense values
    fRefMultEta5 -> SetValueInteger ( -5 ); //not acquired
    fRefMultEta8 -> SetValueInteger ( -5 ); //not acquired
    fNTracks         = -10;

    //Set ZDC variables to defaults
    fZncEnergy->SetValue(-1e6);
    fZpcEnergy->SetValue(-1e6);
    fZnaEnergy->SetValue(-1e6);
    fZpaEnergy->SetValue(-1e6);
    fZem1Energy->SetValue(-1e6);
    fZem2Energy->SetValue(-1e6);
    fZnaTower->SetValue(-1e6);
    fZncTower->SetValue(-1e6);
    fZpaTower->SetValue(-1e6);
    fZpcTower->SetValue(-1e6);
    fZnaFired->SetValueInteger( 0 );
    fZncFired->SetValueInteger( 0 );
    fZpaFired->SetValueInteger( 0 );
    fZpcFired->SetValueInteger( 0 );
    if(lVerbose) Printf("Doing ESD/AOD part...");
    if (lVevent->InheritsFrom("AliESDEvent")) {
        AliESDEvent *esdevent = dynamic_cast<AliESDEvent *>(lVevent);

        //Standard GetReferenceMultiplicity Estimator (0.5 and 0.8)
        fRefMultEta5 -> SetValueInteger ( fESDtrackCuts->GetReferenceMultiplicity(esdevent, AliESDtrackCuts::kTrackletsITSTPC,0.5) );
        fRefMultEta8 -> SetValueInteger ( fESDtrackCuts->GetReferenceMultiplicity(esdevent, AliESDtrackCuts::kTrackletsITSTPC,0.8) );

        //Use fallback in case of return value of -3 or -4
        //This is what will happen in AODs: -3 will use fallback
        //HOWEVER: -4 will not. Inconsistency requires use of "HasNoInconsistentSPDandTrackVertices"!
        if ( fRefMultEta5 -> GetValueInteger() < -2 ) {
            fRefMultEta5 -> SetValueInteger ( fESDtrackCuts->GetReferenceMultiplicity(esdevent, AliESDtrackCuts::kTracklets,0.5) );
        }
        if ( fRefMultEta8 -> GetValueInteger() < -2 ) {
            fRefMultEta8 -> SetValueInteger ( fESDtrackCuts->GetReferenceMultiplicity(esdevent, AliESDtrackCuts::kTracklets,0.8) );
        }

        //A.T.
        fNTracks    = fTrackCuts ? (Short_t)fTrackCuts->GetReferenceMultiplicity(esdevent,kTRUE):-1;

        // ***** ZDC info
        AliESDZDC *lESDZDC = esdevent->GetESDZDC();
        Float_t CalF=0;
        if (lESDZDC->AliESDZDC::TestBit(AliESDZDC::kEnergyCalibratedSignal))  CalF=1.0; //! if zdc is calibrated (in pass2)
        else CalF=8.0;

        fZncEnergy -> SetValue ( (Float_t) (lESDZDC->GetZDCN1Energy())/CalF );
        fZpcEnergy -> SetValue ( (Float_t) (lESDZDC->GetZDCP1Energy())/CalF );
        fZnaEnergy -> SetValue ( (Float_t) (lESDZDC->GetZDCN2Energy())/CalF );
        fZpaEnergy -> SetValue ( (Float_t) (lESDZDC->GetZDCP2Energy())/CalF );

        fZem1Energy -> SetValue ( (Float_t) (lESDZDC->GetZDCEMEnergy(0))/CalF );
        fZem2Energy -> SetValue ( (Float_t) (lESDZDC->GetZDCEMEnergy(1))/CalF );

        Int_t detCh_ZNA = lESDZDC->GetZNATDCChannel();
        Int_t detCh_ZNC = lESDZDC->GetZNCTDCChannel();
        Int_t detCh_ZPA = lESDZDC->GetZPATDCChannel();
        Int_t detCh_ZPC = lESDZDC->GetZPCTDCChannel();

        for (Int_t j = 0; j < 4; ++j) {
	   if (lESDZDC->GetZDCTDCData(detCh_ZNA,j) != 0)      fZnaFired -> SetValueInteger(1);
	   if (lESDZDC->GetZDCTDCData(detCh_ZNC,j) != 0)      fZncFired -> SetValueInteger(1);
	   if (lESDZDC->GetZDCTDCData(detCh_ZPA,j) != 0)      fZpaFired -> SetValueInteger(1);
	   if (lESDZDC->GetZDCTDCData(detCh_ZPC,j) != 0)      fZpcFired -> SetValueInteger(1);
        }

        const Double_t *ZNAtower = lESDZDC->GetZNATowerEnergy();
        const Double_t *ZNCtower = lESDZDC->GetZNCTowerEnergy();
        const Double_t *ZPAtower = lESDZDC->GetZPATowerEnergy();
        const Double_t *ZPCtower = lESDZDC->GetZPCTowerEnergy();
        if (fZnaFired->GetValueInteger() > 0) fZnaTower -> SetValue ( (Float_t) ZNAtower[0] );
        if (fZncFired->GetValueInteger() > 0) fZncTower -> SetValue ( (Float_t) ZNCtower[0] );
        if (fZpaFired->GetValueInteger() > 0) fZpaTower -> SetValue ( (Float_t) ZPAtower[0] );
        if (fZpcFired->GetValueInteger() > 0) fZpcTower -> SetValue ( (Float_t) ZPCtower[0] );

    } else if (lVevent->InheritsFrom("AliAODEvent")) {
        AliAODEvent *aodevent = dynamic_cast<AliAODEvent *>(lVevent);
        AliAODHeader * header = dynamic_cast<AliAODHeader*>(aodevent->GetHeader());
        fRefMultEta5 -> SetValueInteger ( header->GetRefMultiplicityComb05() );
        fRefMultEta8 -> SetValueInteger ( header->GetRefMultiplicityComb08() );

        //FIXME: get ZDC information in AOD in a fully consistent way
        AliAODZDC *lAODZDC = aodevent->GetZDCData();

        if (TMath::Abs(lAODZDC->GetZNATime()) < 990.) fZnaFired -> SetValueInteger(1);
        if (TMath::Abs(lAODZDC->GetZNCTime()) < 990.) fZncFired -> SetValueInteger(1);
        if (TMath::Abs(lAODZDC->GetZPATime()) < 990.) fZpaFired -> SetValueInteger(1);
        if (TMath::Abs(lAODZDC->GetZPCTime()) < 990.) fZpcFired -> SetValueInteger(1);

        const Double_t *ZNAtower = lAODZDC->GetZNATowerEnergy();
        const Double_t *ZNCtower = lAODZDC->GetZNCTowerEnergy();
        const Double_t *ZPAtower = lAODZDC->GetZPATowerEnergy();
        const Double_t *ZPCtower = lAODZDC->GetZPCTowerEnergy();
        if (fZnaFired->GetValueInteger() > 0) fZnaTower -> SetValue ( (Float_t) ZNAtower[0] );
        if (fZncFired->GetValueInteger() > 0) fZncTower -> SetValue ( (Float_t) ZNCtower[0] );
        if (fZpaFired->GetValueInteger() > 0) fZpaTower -> SetValue ( (Float_t) ZPAtower[0] );
        if (fZpcFired->GetValueInteger() > 0) fZpcTower -> SetValue ( (Float_t) ZPCtower[0] );
    }

    //===============================================
    // End part which requires AOD/ESD separation
    //===============================================
    if(lVerbose) Printf("Add info if asked...");
    if ( fkAddInfo ) { //Master switch for users
        //===============================================
        // Compute Percentiles
        //===============================================
        //Make sure OADB is loaded
        SetupRun( lVevent );

        //===============================================
        // I/O: Create object for storing, add
        //===============================================

        if(lVerbose) Printf( "--- Evaluate -1-");
        //Evaluate Estimators from Variables
        AliMultSelection*     lSelection = fOadbMultSelection->GetMultSelection();
        AliMultSelectionCuts* lMultCuts  = fOadbMultSelection->GetEventCuts();
        if(lVerbose) Printf( "--- Evaluate -2-");
        lSelection -> Evaluate (fInput);
        if(lVerbose) Printf( "--- INPUT --- ");
        if(lVerbose) fInput -> Print("V") ;
        if(lVerbose) Printf( "--- OUTPUT --- ");
        if(lVerbose) lSelection -> PrintInfo();
        if(lVerbose) Printf( "--- Evaluate -3-");

        //Event Selection Code: No need to do this for all estimators ...
        lSelection->SetEvSelCode(0); //No Problem!

        if( lMultCuts->GetTriggerCut()    && ! fEvSel_Triggered           )
            lSelection->SetEvSelCode(AliMultSelectionCuts::kRejTrigger);

        if( lMultCuts->GetINELgtZEROCut() && ! fEvSel_INELgtZERO          )
            lSelection->SetEvSelCode(AliMultSelectionCuts::kRejINELgtZERO);

        if( TMath::Abs(fEvSel_VtxZ->GetValue() ) > lMultCuts->GetVzCut()      )
            lSelection->SetEvSelCode(AliMultSelectionCuts::kRejVzCut);

        if( lMultCuts->GetRejectPileupInMultBinsCut() && ! fEvSel_IsNotPileupInMultBins      )
            lSelection->SetEvSelCode(AliMultSelectionCuts::kRejPileupInMultBins);

        if( lMultCuts->GetVertexConsistencyCut()      && ! fEvSel_HasNoInconsistentVertices  )
            lSelection->SetEvSelCode(AliMultSelectionCuts::kRejConsistencySPDandTrackVertices);

        if( lMultCuts->GetTrackletsVsClustersCut()    && ! fEvSel_PassesTrackletVsCluster    )
            lSelection->SetEvSelCode(AliMultSelectionCuts::kRejTrackletsVsClusters);

        if( lMultCuts->GetNonZeroNContribs()    && fnContributors < 1    )
            lSelection->SetEvSelCode(AliMultSelectionCuts::kRejNonZeroNContribs);

        //Just in case you want to store it for debugging
        fEvSelCode = lSelection->GetEvSelCode();

        //Determine Quantiles from calibration histogram
        TH1F *lThisCalibHisto = 0x0;
        TString lThisCalibHistoName;
        Float_t lThisQuantile = -1;
        for(Long_t iEst=0; iEst<lSelection->GetNEstimators(); iEst++) {
            //Changed: no need for run number, object already matches required one
            lThisCalibHistoName = Form("hCalib_%s",lSelection->GetEstimator(iEst)->GetName());
            lThisCalibHisto = 0x0;
            lThisCalibHisto = fOadbMultSelection->GetCalibHisto( lThisCalibHistoName );
            if ( ! lThisCalibHisto ) {
                lThisQuantile = AliMultSelectionCuts::kNoCalib;
                if( iEst < fNDebug ) fQuantiles[iEst] = lThisQuantile;
                lSelection->GetEstimator(iEst)->SetPercentile(lThisQuantile);
            } else {
                lThisQuantile = lThisCalibHisto->GetBinContent( lThisCalibHisto->FindBin( lSelection->GetEstimator(iEst)->GetValue() ));
                if( iEst < fNDebug ) {
                    fQuantiles[iEst] = lThisQuantile; //Debug, please
                }
                lSelection->GetEstimator(iEst)->SetPercentile(lThisQuantile);
            }
        }

        //=============================================================================
        // Fill in quick debug information (available in any execution)

        Float_t lV0M = lSelection->GetMultiplicityPercentile("V0M");
        Float_t lCL0 = lSelection->GetMultiplicityPercentile("CL0");
        Float_t lCL1 = lSelection->GetMultiplicityPercentile("CL1");
        Int_t ltracklets = fnTracklets->GetValueInteger();

        fHistQA_V0M -> Fill( lV0M );
        fHistQA_CL0 -> Fill( lCL0 );
        fHistQA_CL1 -> Fill( lCL1 );

        fHistQA_TrackletsVsV0M -> Fill( lV0M, ltracklets );
        fHistQA_TrackletsVsCL0 -> Fill( lCL0, ltracklets );
        fHistQA_TrackletsVsCL1 -> Fill( lCL1, ltracklets );

        lV0M = lSelection->GetMultiplicityPercentile("V0M",kTRUE);
        lCL0 = lSelection->GetMultiplicityPercentile("CL0",kTRUE);
        lCL1 = lSelection->GetMultiplicityPercentile("CL1",kTRUE);

        fHistQASelected_V0M -> Fill( lV0M );
        fHistQASelected_CL0 -> Fill( lCL0 );
        fHistQASelected_CL1 -> Fill( lCL1 );

        fHistQASelected_TrackletsVsV0M -> Fill( lV0M, ltracklets );
        fHistQASelected_TrackletsVsCL0 -> Fill( lCL0, ltracklets );
        fHistQASelected_TrackletsVsCL1 -> Fill( lCL1, ltracklets );
        //=============================================================================

        //Add to AliVEvent
        //if( (!(InputEvent()->FindListObject("MultSelection")) ) && !fkAttached ) {
        //    InputEvent()->AddObject(fSelection);
        //    fkAttached = kTRUE;
        //}
        // Here, we need to get the object directly from the event.
        // The object is deleted on change of input file by the event
        // reset procedure, so in case the object has disappeared, we
        // need to re-add it.  This is particulary needed on Proof(Lite)
        AliVEvent*        input = InputEvent();

        //=============================================================================
        //if on-the-fly AOD generation, save AliMultSelection also there, please!
        AliVEventHandler *lhandler = 0x0;
        lhandler = AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler();
        AliAODHandler* lAodOutputHandler = dynamic_cast<AliAODHandler*> (lhandler);
        if ( lAodOutputHandler ) {
            //Add to output as well as input
            AliVEvent *lOutputEv = lAodOutputHandler->GetAOD();
            TObject*          outaodO  = lOutputEv->FindListObject("MultSelection");
            AliMultSelection* outaodS  = 0;
            if (!outaodO) {
                outaodS = new AliMultSelection(*lSelection);
                outaodS->SetName("MultSelection");
                lAodOutputHandler->AddBranch("AliMultSelection",&outaodS);
            }
            else {
                outaodS = static_cast<AliMultSelection*>(outaodO);
                outaodS->Set(lSelection);
            }
        }
        //=============================================================================

        TObject*          outO  = input->FindListObject("MultSelection");
        AliMultSelection* outS  = 0;
        if (!outO) {
            outS = new AliMultSelection(*lSelection);
            outS->SetName("MultSelection");
            input->AddObject(outS);
        }
        else {
            outS = static_cast<AliMultSelection*>(outO);
            outS->Set(lSelection);
        }
    }
    if(lVerbose) Printf( "--- INTERMEDIATE --- ");
    if(lVerbose) fInput -> Print("V") ;
    //Event-level fill
    if ( fkCalibration ) {
        //Pre-filter on triggered (kMB) events for saving info
        if( !fkFilterMB || (fkFilterMB && fEvSel_Triggered) ) {
            if(lVerbose) Printf( "--- FILLTREE --- ");
            if(lVerbose) fInput -> Print("V") ;
            fTreeEvent->Fill() ;
        }
    }

    // Post output data.
    PostData(1, fListHist);
    if ( fkCalibration ) PostData(2, fTreeEvent);
}

//________________________________________________________________________
void AliMultSelectionTask::Terminate(Option_t *)
{
    // Draw result to the screen
    // Called once at the end of the query

    TList *cRetrievedList = 0x0;
    cRetrievedList = (TList*)GetOutputData(1);
    if(!cRetrievedList) {
        Printf("ERROR - AliMultSelectionTask : ouput data container list not available\n");
        return;
    }

    fHistEventCounter = dynamic_cast<TH1D*> (  cRetrievedList->FindObject("fHistEventCounter")  );
    if (!fHistEventCounter) {
        Printf("ERROR - AliMultSelectionTask : fHistEventCounter not available");
        return;
    }

    TCanvas *canCheck = new TCanvas("AliMultSelectionTaskQA","Event Counters",10,10,510,510);
    canCheck->cd(1)->SetLogy();

    fHistEventCounter->SetMarkerStyle(22);
    fHistEventCounter->DrawCopy("E");
}

//________________________________________________________________________
Int_t AliMultSelectionTask::SetupRun(const AliVEvent* const esd)
{
    // Setup files for run

    if (!esd)
        return -1;

    // check if something to be done
    if (fCurrentRun == esd->GetRunNumber())
        return 0;
    else
        fCurrentRun = esd->GetRunNumber();
    AliInfoF("Detected run number: %i",fCurrentRun);

    TString lPathInput = CurrentFileName();

    //Autodetect Period Name
    TString lPeriodName     = GetPeriodNameByRunNumber();
    AliWarning("Autodetecting production name via LPM tag...");
    TString lProductionName = GetPeriodNameByLPM();
    if( lProductionName.EqualTo("") ){
	AliWarning("LPM Tag didn't work, autodetecting via path..."); 
	lProductionName = GetPeriodNameByPath( lPathInput );
	if (!lProductionName.Contains("LHC")){
	   AliWarning("Autodetect via path seems to not have worked?"); 
	}
    }

    AliWarning("==================================================");
    AliWarning(Form(" Period Name (by run number)....: %s", lPeriodName.Data()));
    AliWarning(Form(" Production Name (by path)......: %s", lProductionName.Data()));
    AliWarning("==================================================");
    if ( lPeriodName.EqualTo(lProductionName.Data()) == kTRUE ) {
        AliWarning( Form(" Assumed to be DATA ANALYSIS on period %s",lPeriodName.Data() ));
        AliWarning("==================================================");
    }
    if ( fAlternateOADBForEstimators.EqualTo("")==kTRUE && lPeriodName.EqualTo(lProductionName.Data()) == kFALSE ) {
        AliWarning(" Auto-detected that this is MC, but you didn't provide a production name!");
        AliWarning(Form(" Auto-detection production name is %s, checking if there...",lProductionName.Data()));
	Bool_t lItsThere = CheckOADB( lProductionName ); 
	if( !lItsThere ){ 
		//Override options: go to default OADBs depending on generator
		AliWarning(" OADB for this MC does not exist! Checking generator type..."); 
		Bool_t lItsHijing = IsHijing(); 
		Bool_t lItsDPMJet = IsDPMJet(); 
		if ( lItsHijing ){ 
			lProductionName = Form("%s-DefaultMC-Hijing",lPeriodName.Data());
			AliWarning(Form(" This is HIJING! Will use OADB named %s",lProductionName.Data()));
		}
		if ( lItsDPMJet ){ 
			lProductionName = Form("%s-DefaultMC-DPMJet",lPeriodName.Data());
			AliWarning(Form(" This is DPMJet! Will use OADB named %s",lProductionName.Data())); 
		}
		if ( (!lItsHijing) && (!lItsDPMJet) ){ 
			AliWarning(" Unable to detect generator type from header. Sorry."); 
		}	
	}else{ 
		AliWarning(" OADB for this period exists. Proceeding as usual."); 
	}	
        fAlternateOADBForEstimators = lProductionName;
        AliInfo("==================================================");
    }

    if ( lPeriodName.EqualTo("Empty") ){ 
	//This is uncalibrated data, skip all
	AliWarning("This is uncalibrated data, will generate empty OADB!");
	AliWarning("This means kNoCalib (=199) will be returned everywhere.");  
	CreateEmptyOADB(); 
	return -1;  
    }

    //Determine location of file to open: default OADB
    TString fileName =(Form("%s/COMMON/MULTIPLICITY/data/OADB-%s.root", AliAnalysisManager::GetOADBPath(), lPeriodName.Data() ));
    AliInfo(Form("Setup Multiplicity Selection for run %d with file %s, period: %s\n",fCurrentRun,fileName.Data(),lPeriodName.Data()));

    //Full Manual Bypass Mode (DEBUG ONLY)
    if ( fAlternateOADBFullManualBypass.EqualTo("")==kFALSE ) {
        AliWarning(" Extra option detected: FULL MANUAL BYPASS of DATA OADB Location ");
        AliWarning(" --- Warning: Use with care ---");
        AliWarning(Form(" New complete path: %s", fAlternateOADBFullManualBypass.Data() ));
        fileName = Form("%s", fAlternateOADBFullManualBypass.Data() );
    }

    AliOADBContainer *con = new AliOADBContainer("OADB");
    Int_t lFoundFile = con->InitFromFile(fileName,"MultSel");

    //FIXME: Here one should open an empty file instead...
    TString fileNameDef =(Form("%s/COMMON/MULTIPLICITY/data/OADB-LHC15f.root", AliAnalysisManager::GetOADBPath() ));
    if ( lFoundFile == 1 ) {
        lFoundFile = con->InitFromFile(fileNameDef,"MultSel");
        if (lFoundFile == 1)
            AliFatal("Check AliPhysics Installation: nothing found in LHC15f calibration!");
    }

    //Get Object for this run!
    TObject *lObjAcquired = 0x0;

    lObjAcquired = con->GetObject(fCurrentRun, "Default");

    if (!lObjAcquired) {
        if ( fkUseDefaultCalib ) {
            AliWarning("======================================================================");
            AliWarning(Form(" Multiplicity OADB does not exist for run %d, using Default \n",fCurrentRun ));
            AliWarning(" This is only a 'good guess'! Use with Care! ");
            AliWarning(" To Switch off this good guess, use SetUseDefaultCalib(kFALSE)");
            AliWarning("======================================================================");
            lObjAcquired  = con->GetDefaultObject("oadbDefault");
        } else {
            AliWarning("======================================================================");
            AliWarning(Form(" Multiplicity OADB does not exist for run %d, will return kNoCalib!",fCurrentRun ));
            AliWarning("======================================================================");
            //Create an empty OADB for us, please !
            CreateEmptyOADB();
            return -1;
        }
    }
    if (!lObjAcquired) {
        // This should not happen...
        AliFatal("Really cannot find any OADB object - giving up!");
    }

    AliOADBMultSelection *lObjTypecast = (AliOADBMultSelection*) lObjAcquired;

    fOadbMultSelection = new AliOADBMultSelection(*lObjTypecast);
    // De-couple histograms from the underlying file
    fOadbMultSelection->Dissociate();
    // Update cache map from estimator to histogram for fast look-up
    fOadbMultSelection->Setup();

    //Make sure naming convention is followed!
    AliMultSelection* sel = fOadbMultSelection->GetMultSelection();
    if (sel) {
        sel->SetName("MultSelection");
    }

    //=====================================================================
    //Option to override estimators from alternate oadb file
    if ( fAlternateOADBForEstimators.EqualTo("")==kFALSE ) {
        AliWarning("Extra option detected: Load estimators from OADB file called: ");
        AliWarning(Form(" path: %s", fAlternateOADBForEstimators.Data() ));

        TString fileNameAlter =(Form("%s/COMMON/MULTIPLICITY/data/OADB-%s.root", AliAnalysisManager::GetOADBPath(), fAlternateOADBForEstimators.Data() ));

        //Full Manual Bypass Mode (DEBUG ONLY)
        if ( fAlternateOADBFullManualBypassMC.EqualTo("")==kFALSE ) {
            AliWarning(" Extra option detected: FULL MANUAL BYPASS of MONTE CARLO OADB Location ");
            AliWarning(" --- Warning: Use with care ---");
            AliWarning(Form(" New complete path: %s", fAlternateOADBFullManualBypassMC.Data() ));
            fileNameAlter = Form("%s", fAlternateOADBFullManualBypassMC.Data() );
        }

        AliOADBContainer *conAlter = new AliOADBContainer("OADB-Alternate");
        Int_t lFoundFileAlter = conAlter->InitFromFile(fileNameAlter,"MultSel");
        if ( lFoundFileAlter == 1 ) {
            AliFatal("Couldn't find requested alternate calibration! Quitting!");
        }

        //Get Object for this run
        TObject *lObjAcquiredAlter = 0x0;
        lObjAcquiredAlter = conAlter->GetObject(fCurrentRun, "Default");
        if (!lObjAcquiredAlter) {
            if ( fkUseDefaultMCCalib ) {
                AliWarning("======================================================================");
                AliWarning(Form(" MC Multiplicity OADB does not exist for run %d, using Default \n",fCurrentRun ));
                AliWarning(" This is usually only approximately OK! Use with Care! ");
                AliWarning(" To Switch off this good guess, use SetUseDefaultMCCalib(kFALSE)");
                AliWarning("======================================================================");
                lObjAcquiredAlter  = conAlter->GetDefaultObject("oadbDefault");
            } else {
                AliWarning("======================================================================");
                AliWarning(Form(" MC Multiplicity OADB does not exist for run %d, will return kNoCalib!",fCurrentRun ));
                AliWarning("======================================================================");
                //Create an empty OADB for us, please !
                CreateEmptyOADB();
                return -1;
            }
        }
        if (!lObjAcquiredAlter) {
            // This should not happen...
            AliFatal("Really cannot find any OADB object - giving up!");
        }

        //Actually, it's not required that we keep a copy of this object in memory. We only need to grab
        //the definitions... This can be much optimized!
        AliOADBMultSelection *fOadbMultSelectionAlter = (AliOADBMultSelection*) lObjAcquiredAlter;
        AliMultSelection* selAlter = fOadbMultSelectionAlter->GetMultSelection();

        //Sweep all estimators from standard OADB and replace their definitions...
        TString lTempStr;
        for(Int_t iEst=0; iEst<sel->GetNEstimators(); iEst++) {
            lTempStr = sel->GetEstimator(iEst)->GetName();
            AliMultEstimator *lEstim = selAlter->GetEstimator( lTempStr.Data() );
            if ( !lEstim ) {
                AliWarning(Form("Estimator named %s: no scaling factor applied!",sel->GetEstimator(iEst)->GetName()));
            } else {
                sel->GetEstimator(iEst)->SetDefinition ( lEstim->GetDefinition().Data() );
            }
        }
        //That should be it...
    }
    //=====================================================================

    if (sel) {
        sel->SetName("MultSelection");
        //Optimize evaluation
        sel->Setup(fInput);
    }

    AliInfo("---> Successfully set up! Inspect MultSelection:");
    if (sel) {
        sel->PrintInfo();
    }
    else {
        AliWarning("Weird! No AliMultSelection found...");
    }
    AliInfo("---> Inspect Event Selection Criteria:");
    AliMultSelectionCuts* selcuts = fOadbMultSelection->GetEventCuts();
    if (selcuts) {
        selcuts->Print();
    }
    else {
        AliWarning("Weird! No AliMultSelectionCuts found...");
    }

    return 0;
}


//______________________________________________________________________
Bool_t AliMultSelectionTask::IsSelectedTrigger(AliVEvent* event, AliVEvent::EOfflineTriggerTypes lCheckedTrig)
// Function to check for a specific trigger class available in AliVEvent (default AliVEvent::kMB)
{
    //Code to reject events that aren't trigType
    UInt_t maskIsSelected = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
    Bool_t isSelected = 0;
    isSelected = maskIsSelected & lCheckedTrig;
    return isSelected;
}

//______________________________________________________________________
Bool_t AliMultSelectionTask::IsINELgtZERO(AliVEvent *event)
// Function to check for INEL > 0 condition
// Makes use of tracklets and requires at least and SPD vertex
{
    Bool_t lReturnValue = kFALSE;
    //Use Ref.Mult. code...
    if (event->InheritsFrom("AliESDEvent")) {
        AliESDEvent *esdevent = dynamic_cast<AliESDEvent *>(event);
        if (!esdevent) return kFALSE;
        if ( AliESDtrackCuts::GetReferenceMultiplicity(esdevent, AliESDtrackCuts::kTracklets, 1.0) >= 1 ) lReturnValue = kTRUE;
    }
    //Redo equivalent test
    else if (event->InheritsFrom("AliAODEvent")) {
        AliAODEvent *aodevent = dynamic_cast<AliAODEvent *>(event);
        if (!aodevent) return kFALSE;

        //FIXME --- Actually, here we can come up with a workaround.
        // We can check for the reference multiplicity stored and look for error codes!

        AliAODHeader * header = dynamic_cast<AliAODHeader*>(aodevent->GetHeader());
        Int_t lStoredRefMult = header->GetRefMultiplicityComb08();

        //Get Multiplicity object
        AliAODTracklets *spdmult = aodevent->GetMultiplicity();
        for (Int_t i=0; i<spdmult->GetNumberOfTracklets(); ++i)
        {
            if ( lStoredRefMult != -1 && lStoredRefMult != -2 && TMath::Abs(spdmult->GetEta(i)) < 1.0 ) lReturnValue = kTRUE;
        }
    }
    return lReturnValue;
}

//______________________________________________________________________
Bool_t AliMultSelectionTask::IsAcceptedVertexPosition(AliVEvent *event)
// Simple check for the best primary vertex Z position:
// Will accept events only if |z| < 10cm
{
    Bool_t lReturnValue = kFALSE;
    //Getting around to the best vertex -> typecast to ESD/AOD
    const AliVVertex *lPrimaryVtx = NULL;
    /* get ESD vertex */
    if (event->InheritsFrom("AliESDEvent")) {
        AliESDEvent *esdevent = dynamic_cast<AliESDEvent *>(event);
        if (!esdevent) return kFALSE;
        lPrimaryVtx = esdevent->GetPrimaryVertex();
    }
    /* get AOD vertex */
    else if (event->InheritsFrom("AliAODEvent")) {
        AliAODEvent *aodevent = dynamic_cast<AliAODEvent *>(event);
        if (!aodevent) return kFALSE;
        lPrimaryVtx = aodevent->GetPrimaryVertex();
    }
    if ( TMath::Abs( lPrimaryVtx->GetZ() ) <= 10.0 ) lReturnValue = kTRUE;
    return lReturnValue;
}

//______________________________________________________________________
Bool_t AliMultSelectionTask::HasNoInconsistentSPDandTrackVertices(AliVEvent *event)
// This function checks if track and SPD vertices are consistent.
// N.B.: It is rigorously a "Not Inconsistent" function which will
// let events with only SPD vertex go through without troubles.
{
    //It's consistent until proven otherwise...
    Bool_t lReturnValue = kTRUE;

    //Getting around to the best vertex -> typecast to ESD/AOD


    /* get ESD vertex */
    if (event->InheritsFrom("AliESDEvent")) {
        AliESDEvent *esdevent = dynamic_cast<AliESDEvent *>(event);
        if (!esdevent) return kFALSE;
        const AliESDVertex *lPrimaryVtxSPD    = NULL;
        const AliESDVertex *lPrimaryVtxTracks = NULL;

        lPrimaryVtxSPD    = esdevent->GetPrimaryVertexSPD   ();
        lPrimaryVtxTracks = esdevent->GetPrimaryVertexTracks();

        //Only continue if track vertex defined
        if( lPrimaryVtxTracks->GetStatus() && lPrimaryVtxSPD->GetStatus() ) {
            //Copy-paste from refmult estimator
            // TODO value of displacement to be studied
            const Float_t maxDisplacement = 0.5;
            //check for displaced vertices
            Double_t displacement = TMath::Abs(lPrimaryVtxSPD->GetZ() - lPrimaryVtxTracks->GetZ());
            if (displacement > maxDisplacement) lReturnValue = kFALSE;
        }
    }
    /* get AOD vertex */
    else if (event->InheritsFrom("AliAODEvent")) {
        AliAODEvent *aodevent = dynamic_cast<AliAODEvent *>(event);
        if (!aodevent) return kFALSE;

        //FIXME - Hack to deal with the fact that no
        //        AliAODEvent::GetPrimaryVertexTracks() exists...
        AliAODHeader * header = dynamic_cast<AliAODHeader*>(aodevent->GetHeader());
        Int_t lStoredRefMult = header->GetRefMultiplicityComb08();
        if( lStoredRefMult == -4 ) lReturnValue = kFALSE;
    }
    return lReturnValue;
}

//______________________________________________________________________
Bool_t AliMultSelectionTask::IsNotPileupSPDInMultBins(AliVEvent *event)
// Checks if not pileup from SPD (via IsPileupFromSPDInMultBins)
{
    Bool_t lReturnValue = kTRUE;
    //Getting around to the SPD vertex -> typecast to ESD/AOD
    if (event->InheritsFrom("AliESDEvent")) {
        AliESDEvent *esdevent = dynamic_cast<AliESDEvent *>(event);
        if (!esdevent) return kFALSE;

        //return kTRUE in case no AliMultiplicity object is present!
        //Warning: this essentially means pileup went unchecked!
        if ( !(esdevent->GetMultiplicity()) ) return kTRUE;

        if ( esdevent->IsPileupFromSPDInMultBins() == kTRUE ) lReturnValue = kFALSE;
    }
    else if (event->InheritsFrom("AliAODEvent")) {
        AliAODEvent *aodevent = dynamic_cast<AliAODEvent *>(event);
        if (!aodevent) return kFALSE;

        //return kTRUE in case no AliMultiplicity object is present!
        //Warning: this essentially means pileup went unchecked!
        if ( !(aodevent->GetTracklets()) ) return kTRUE;

        if ( aodevent->IsPileupFromSPDInMultBins() == kTRUE ) lReturnValue = kFALSE;
    }
    return lReturnValue;
}

//______________________________________________________________________
Bool_t AliMultSelectionTask::IsNotPileupSPD(AliVEvent *event)
// Checks if not pileup from SPD (via IsNotPileupSPD)
{
    Bool_t lReturnValue = kTRUE;
    //Getting around to the SPD vertex -> typecast to ESD/AOD
    if (event->InheritsFrom("AliESDEvent")) {
        AliESDEvent *esdevent = dynamic_cast<AliESDEvent *>(event);
        if (!esdevent) return kFALSE;
        if ( esdevent->IsPileupFromSPD() == kTRUE ) lReturnValue = kFALSE;
    }
    else if (event->InheritsFrom("AliAODEvent")) {
        AliAODEvent *aodevent = dynamic_cast<AliAODEvent *>(event);
        if (!aodevent) return kFALSE;
        if ( aodevent->IsPileupFromSPD() == kTRUE ) lReturnValue = kFALSE;
    }
    return lReturnValue;
}

//______________________________________________________________________
Bool_t AliMultSelectionTask::IsNotPileupMV(AliVEvent *event)
// Checks if not pileup from SPD (via IsNotPileupSPD)
{
    //Negation: this is kTRUE if this is NOT pileup
    Bool_t lReturnValue = !(fUtils->IsPileUpMV( event ) );
    return lReturnValue;
}

//______________________________________________________________________
Bool_t AliMultSelectionTask::PassesTrackletVsCluster(AliVEvent* event)
{
    //Negation: this is kTRUE if this is NOT background
    Bool_t lReturnValue = !(fUtils->IsSPDClusterVsTrackletBG ( event ) );
    return lReturnValue;
}

//______________________________________________________________________
TString AliMultSelectionTask::GetPeriodNameByLPM() 
{
    //==================================
    // Setup initial Info
    Bool_t lLocated = kFALSE;
    TString lTag = "LPMProductionTag";
    TString lProductionName = "";

    //==================================
    // Get alirootVersion object title
    AliInputEventHandler* handler = dynamic_cast<AliInputEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    if (!handler) return lProductionName; //failed!
    TObject* prodInfoData = handler->GetUserInfo()->FindObject("alirootVersion");
    if (!prodInfoData) return lProductionName; //failed!
    TString lAlirootVersion(prodInfoData->GetTitle());
    
    //==================================
    // Get Production name
    TObjArray* lArrStr = lAlirootVersion.Tokenize(";");
    if(lArrStr->GetEntriesFast()) {
        TIter iString(lArrStr);
        TObjString* os=0;
        Int_t j=0;
        while ((os=(TObjString*)iString())) {
            if( os->GetString().Contains(lTag.Data()) ){
                lLocated = kTRUE;
                lProductionName = os->GetString().Data();
                //Remove Label
                lProductionName.ReplaceAll(lTag.Data(),"");
                //Remove any remaining whitespace (just in case)
                lProductionName.ReplaceAll("=","");
                lProductionName.ReplaceAll(" ","");
            }
            j++;
        }
    }
    //Memory cleanup
    delete lArrStr;
    //Return production name
    return lProductionName;
}

//______________________________________________________________________
TString AliMultSelectionTask::GetPeriodNameByPath(const TString lPath) const
{
    //============================================================
    //This function is meant to get the period name.
    //IMPORTANT: this cannot solely depend on run number.
    //This has to load different settings in case the period name
    //refers to a Monte Carlo production.
    //============================================================

    //==================================
    // Setup initial Info
    TString lProductionName = lPath.Data();
    //==================================
    AliInfoF(" Autodetecting production name from filename: %s", lProductionName.Data() );
    //==================================
    // Get Production name
    Long_t iOcurrence = 0;

    // 1) Check first occurence of LHC
    iOcurrence = lProductionName.Index("LHC");

    // 2) Remove Anything prior to this, please
    lProductionName.Remove(0,iOcurrence);

    // 3) Check first occurrence of "/" and strip
    iOcurrence = lProductionName.Index("/");
    lProductionName.Remove(iOcurrence, lProductionName.Length() );

    // 4) Check first occurrence of "__" (LEGO Train Executions)
    iOcurrence = lProductionName.Index("__");
    lProductionName.Remove(iOcurrence, lProductionName.Length() );

    //======================================================================
    //
    // --> More checks to be added for generality ?
    //
    // --> Note: if all else fails, this can always be bypassed by explicitly
    //     defining production name calling
    //
    //     --- SetAlternateOADBforEstimators ( production name )
    //
    //     where production name can be real data (e.g. LHC15o)
    //     or MC (e.g. LHC15k1_plus) and all will be fine...
    //
    //======================================================================

    return lProductionName;
}
//______________________________________________________________________
TString AliMultSelectionTask::GetPeriodNameByRunNumber() const
{
    //============================================================
    // This function is meant to get the period name.
    // Note: This corresponds to the period to which this data
    // was anchored, will read from run number!
    //
    // N.B. This will require some bookkeeping of all enabled productions
    //
    //============================================================

    //Will make anything return AliMultSelectionCuts::kNoCalib
    TString lProductionName = "Empty";

    //Registered Productions : Run 1 pp
    if ( fCurrentRun >= 114751 && fCurrentRun <= 117222 ) lProductionName = "LHC10b";
    if ( fCurrentRun >= 118903 && fCurrentRun <= 120829 ) lProductionName = "LHC10c";
    if ( fCurrentRun >= 122374 && fCurrentRun <= 126437 ) lProductionName = "LHC10d";
    if ( fCurrentRun >= 127712 && fCurrentRun <= 130840 ) lProductionName = "LHC10e";

    //Registered Productions : Run 1 Pb-Pb
    if ( fCurrentRun >= 136851 && fCurrentRun <= 139517 ) lProductionName = "LHC10h";

    //Registered Productions : Run 2 pp
    if ( fCurrentRun >= 225000 && fCurrentRun <= 226606 ) lProductionName = "LHC15f";
    if ( fCurrentRun >= 235196 && fCurrentRun <= 236866 ) lProductionName = "LHC15i";
    if ( fCurrentRun >= 237003 && fCurrentRun <= 238622 ) lProductionName = "LHC15j";
    if ( fCurrentRun >= 239319 && fCurrentRun <= 241541 ) lProductionName = "LHC15l";

    //Registered Productions : Run 2 Pb-Pb
    if ( fCurrentRun >= 243395 && fCurrentRun <= 243984 ) lProductionName = "LHC15m";
    if ( fCurrentRun >= 244917 && fCurrentRun <= 300000 ) lProductionName = "LHC15o";

    return lProductionName;
}
//______________________________________________________________________
Bool_t AliMultSelectionTask::CheckOADB(TString lProdName) const { 
    //This helper function checks if an OADB exists for the production named lProdName
    //Determine file name 
    TString fileName = Form("%s/COMMON/MULTIPLICITY/data/OADB-%s.root", AliAnalysisManager::GetOADBPath(), lProdName.Data() );

    Bool_t lInverseThis = !gSystem->AccessPathName(fileName.Data());
    return lInverseThis;
}

//______________________________________________________________________
Bool_t AliMultSelectionTask::IsHijing() const { 
	//Function to check if this is Hijing MC
	Bool_t lReturnValue = kFALSE; 
        AliAnalysisManager* anMan = AliAnalysisManager::GetAnalysisManager();
        AliMCEventHandler* eventHandler = (AliMCEventHandler*)anMan->GetMCtruthEventHandler();
        AliStack*    stack=0;
        AliMCEvent*  mcEvent=0;
        AliGenHijingEventHeader* hHijing=0;
        if (eventHandler && (mcEvent=eventHandler->MCEvent()) && (stack=mcEvent->Stack())) {
            AliGenEventHeader* mcGenH = mcEvent->GenEventHeader();
            if (mcGenH->InheritsFrom(AliGenHijingEventHeader::Class())){
		//Option 1: Just Hijing
		lReturnValue = kTRUE;
	    } else if (mcGenH->InheritsFrom(AliGenCocktailEventHeader::Class())) {
		//Option 2: cocktail involving Hijing
                TList* headers = ((AliGenCocktailEventHeader*)mcGenH)->GetHeaders();
                hHijing = dynamic_cast<AliGenHijingEventHeader*>(headers->FindObject("Hijing"));
		if ( hHijing ) lReturnValue = kTRUE;  
            }		

	}	
	return lReturnValue;
}

//______________________________________________________________________
Bool_t AliMultSelectionTask::IsDPMJet() const { 
    //Function to check if this is DPMJet
        //Function to check if this is Hijing MC
        Bool_t lReturnValue = kFALSE;
        AliAnalysisManager* anMan = AliAnalysisManager::GetAnalysisManager();
        AliMCEventHandler* eventHandler = (AliMCEventHandler*)anMan->GetMCtruthEventHandler();
        AliStack*    stack=0;
        AliMCEvent*  mcEvent=0;
        if (eventHandler && (mcEvent=eventHandler->MCEvent()) && (stack=mcEvent->Stack())) {
            AliGenEventHeader* mcGenH = mcEvent->GenEventHeader();
            if (mcGenH->InheritsFrom(AliGenDPMjetEventHeader::Class())) {
		//DPMJet Header is there!
                lReturnValue = kTRUE;
            }
        }
        return lReturnValue;
}

//______________________________________________________________________
void AliMultSelectionTask::CreateEmptyOADB()
{
    //This will completely reset the OADB, such that any attempt to use
    //the framework will return kNoCalib everywhere: fully safe mode of operation!
    if( fOadbMultSelection ) {
        delete fOadbMultSelection;
    }
    fOadbMultSelection = new AliOADBMultSelection();
    AliMultSelectionCuts *cuts              = new AliMultSelectionCuts(); //irrelevant
    AliMultSelection *fsels                 = new AliMultSelection    (); //is empty, will return kNoCalib always

    fOadbMultSelection->SetEventCuts        ( cuts  );
    fOadbMultSelection->SetMultSelection    ( fsels );
}
