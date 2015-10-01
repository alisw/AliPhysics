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
#include "TFile.h"
#include "THnSparse.h"
#include "TVector3.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TLegend.h"
#include "TObjectTable.h"
//#include "AliLog.h"

#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliESDpid.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliInputEventHandler.h"
#include "AliAnalysisManager.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliESDAD.h" //AD

//#include "AliCFContainer.h"
#include "AliMultiplicity.h"
#include "AliESDUtils.h"
#include "AliGenEventHeader.h"
#include "AliAnalysisTaskSE.h"

//For MultSelection Framework
#include "AliOADBContainer.h"
#include "AliOADBMultSelection.h"
#include "AliMultEstimator.h"
#include "AliMultVariable.h"
#include "AliMultInput.h"
#include "AliMultSelection.h"

//task header
#include "AliMultSelectionTask.h"

using std::cout;
using std::endl;

ClassImp(AliMultSelectionTask)

AliMultSelectionTask::AliMultSelectionTask()
: AliAnalysisTaskSE(), fListHist(0), fTreeEvent(0),fESDtrackCuts(0), fTrackCuts(0),
fkCalibration ( kTRUE ), fkAddInfo(kTRUE), fkFilterMB(kTRUE), fkAttached(0),
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
fnSPDClusters(0),
fnSPDClusters0(0),
fnSPDClusters1(0),
fRefMultEta5(0),
fRefMultEta8(0),
fRunNumber(0),
fEvSel_HasAtLeastSPDVertex(0),
fEvSel_VtxZCut(0),
fEvSel_IsNotPileup(0),
fEvSel_IsNotPileupMV(0),
fEvSel_IsNotPileupInMultBins(0),
fEvSel_Triggered(0),
fEvSel_INELgtZERO(0),
fEvSel_nTracklets(0),
fEvSel_nTrackletsEta10(0),
fEvSel_VtxZ(0),
fEvSel_nSPDPrimVertices(0),
fEvSel_distZ(0),
fEvSel_nContributors(0),
//Histos
fHistEventCounter(0),
oadbMultSelection(0),
fMultCuts(0),
fSelection(0),
fInput(0)
//------------------------------------------------
// Tree Variables
{

}

AliMultSelectionTask::AliMultSelectionTask(const char *name)
    : AliAnalysisTaskSE(name), fListHist(0), fTreeEvent(0), fESDtrackCuts(0), fTrackCuts(0),
fkCalibration ( kTRUE ), fkAddInfo(kTRUE), fkFilterMB(kTRUE), fkAttached(0),
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
fnSPDClusters(0),
fnSPDClusters0(0),
fnSPDClusters1(0),
fRefMultEta5(0),
fRefMultEta8(0),
fRunNumber(0),
fEvSel_HasAtLeastSPDVertex(0),
fEvSel_VtxZCut(0),
fEvSel_IsNotPileup(0),
fEvSel_IsNotPileupMV(0),
fEvSel_IsNotPileupInMultBins(0),
fEvSel_Triggered(0),
fEvSel_INELgtZERO(0),
fEvSel_nTracklets(0),
fEvSel_nTrackletsEta10(0),
fEvSel_VtxZ(0),
fEvSel_nSPDPrimVertices(0),
fEvSel_distZ(0),
fEvSel_nContributors(0),
//Histos
fHistEventCounter(0),
oadbMultSelection(0),
fMultCuts(0),
fSelection(0),
fInput(0)
{
    DefineOutput(1, TList::Class()); // Event Counter Histo
    if (fkCalibration) DefineOutput(2, TTree::Class()); // Event Tree
}


AliMultSelectionTask::~AliMultSelectionTask()
{
    //------------------------------------------------
    // DESTRUCTOR
    //------------------------------------------------

    if (fListHist) {
        delete fListHist;
        fListHist = 0x0;
    }
    if (fTreeEvent) {
        delete fTreeEvent;
        fTreeEvent = 0x0;
    }
}

 

//________________________________________________________________________
void AliMultSelectionTask::UserCreateOutputObjects()
{

    OpenFile(2);
    // Called once

    //------------------------------------------------

    
    //Create Input Information
    fInput = new AliMultInput();
    
    //Create input variables in AliMultInput Class
    //V0 related
    fAmplitude_V0A        = new AliMultVariable("fAmplitude_V0A");
    fAmplitude_V0C        = new AliMultVariable("fAmplitude_V0C");
    fAmplitude_V0Apartial = new AliMultVariable("fAmplitude_V0Apartial");
    fAmplitude_V0Cpartial = new AliMultVariable("fAmplitude_V0Cpartial");
    fAmplitude_V0AEq      = new AliMultVariable("fAmplitude_V0AEq");
    fAmplitude_V0CEq      = new AliMultVariable("fAmplitude_V0CEq");
    fAmplitude_OnlineV0A  = new AliMultVariable("fAmplitude_OnlineV0A");
    fAmplitude_OnlineV0C  = new AliMultVariable("fAmplitude_OnlineV0C");
    //SPD Related
    fnSPDClusters         = new AliMultVariable("fnSPDClusters");
    fnSPDClusters->SetIsInteger( kTRUE );
    //AD Related
    fMultiplicity_ADA     = new AliMultVariable("fMultiplicity_ADA");
    fMultiplicity_ADC     = new AliMultVariable("fMultiplicity_ADC");
    
    //Add to AliMultInput Object, will later bind to TTree object
    fInput->AddVariable( fAmplitude_V0A );
    fInput->AddVariable( fAmplitude_V0C );
    fInput->AddVariable( fAmplitude_V0Apartial );
    fInput->AddVariable( fAmplitude_V0Cpartial );
    fInput->AddVariable( fAmplitude_V0AEq );
    fInput->AddVariable( fAmplitude_V0CEq );
    fInput->AddVariable( fAmplitude_OnlineV0A );
    fInput->AddVariable( fAmplitude_OnlineV0C );
    fInput->AddVariable( fnSPDClusters );
    fInput->AddVariable( fMultiplicity_ADA );
    fInput->AddVariable( fMultiplicity_ADC );
    
    fTreeEvent = new TTree("fTreeEvent","Event");

    //------------------------------------------------
    // fTree Branch definitions - Event by Event info
    //------------------------------------------------

    //-----------BASIC-INFO---------------------------
    
    // A.T. (my suggestion for V0: 1..4 replacing partial)
    fTreeEvent->Branch("fAmplitude_V0A1",&fAmplitude_V0A1,"fAmplitude_V0A1/F");
    fTreeEvent->Branch("fAmplitude_V0A2",&fAmplitude_V0A2,"fAmplitude_V0A2/F");
    fTreeEvent->Branch("fAmplitude_V0A3",&fAmplitude_V0A3,"fAmplitude_V0A3/F");
    fTreeEvent->Branch("fAmplitude_V0A4",&fAmplitude_V0A4,"fAmplitude_V0A4/F");
    fTreeEvent->Branch("fAmplitude_V0C1",&fAmplitude_V0C1,"fAmplitude_V0C1/F");
    fTreeEvent->Branch("fAmplitude_V0C2",&fAmplitude_V0C2,"fAmplitude_V0C2/F");
    fTreeEvent->Branch("fAmplitude_V0C3",&fAmplitude_V0C3,"fAmplitude_V0C3/F");
    fTreeEvent->Branch("fAmplitude_V0C4",&fAmplitude_V0C4,"fAmplitude_V0C4/F");

    //Official GetReferenceMultiplicity
    fTreeEvent->Branch("fRefMultEta5",&fRefMultEta5,"fRefMultEta5/I");
    fTreeEvent->Branch("fRefMultEta8",&fRefMultEta8,"fRefMultEta8/I");

    //Run Number
    fTreeEvent->Branch("fRunNumber", &fRunNumber, "fRunNumber/I");

    //Booleans for Event Selection
    fTreeEvent->Branch("fEvSel_HasAtLeastSPDVertex", &fEvSel_HasAtLeastSPDVertex, "fEvSel_HasAtLeastSPDVertex/O");
    fTreeEvent->Branch("fEvSel_VtxZCut", &fEvSel_VtxZCut, "fEvSel_VtxZCut/O");
    fTreeEvent->Branch("fEvSel_IsNotPileup", &fEvSel_IsNotPileup, "fEvSel_IsNotPileup/O");
    fTreeEvent->Branch("fEvSel_IsNotPileupMV", &fEvSel_IsNotPileupMV, "fEvSel_IsNotPileupMV/O");
    fTreeEvent->Branch("fEvSel_IsNotPileupInMultBins", &fEvSel_IsNotPileupInMultBins, "fEvSel_IsNotPileupInMultBins/O");
    fTreeEvent->Branch("fEvSel_Triggered", &fEvSel_Triggered, "fEvSel_Triggered/O");
    fTreeEvent->Branch("fEvSel_INELgtZERO", &fEvSel_INELgtZERO, "fEvSel_INELgtZERO/O");

    //Tracklets vs clusters test
    fTreeEvent->Branch("fEvSel_nTracklets",      &fEvSel_nTracklets, "fEvSel_nTracklets/I");
    fTreeEvent->Branch("fEvSel_nTrackletsEta10", &fEvSel_nTrackletsEta10, "fEvSel_nTrackletsEta10/I");
    //A.T.
    //fTreeEvent->Branch("fnSPDClusters0", &nSPDClusters0, "fnSPDClusters0/I");
    //fTreeEvent->Branch("fnSPDClusters1", &nSPDClusters1, "fnSPDClusters1/I");
    fTreeEvent->Branch("fNTracks",      &fNTracks, "fNTracks/I");

    fTreeEvent->Branch("fEvSel_VtxZ", &fEvSel_VtxZ, "fEvSel_VtxZ/F");
    fTreeEvent->Branch("fEvSel_nSPDPrimVertices", &fEvSel_nSPDPrimVertices, "fEvSel_nSPDPrimVertices/I");
    fTreeEvent->Branch("fEvSel_distZ", &fEvSel_distZ, "fEvSel_distZ/F");
    fTreeEvent->Branch("fEvSel_nContributors", &fEvSel_nContributors, "fEvSel_nContributors/I");
    
    //A.T.
    //ZDC info
    fTreeEvent->Branch("fZncEnergy", &fZncEnergy, "fZncEnergy/F");
    fTreeEvent->Branch("fZpcEnergy", &fZpcEnergy, "fZpcEnergy/F");
    fTreeEvent->Branch("fZnaEnergy", &fZnaEnergy, "fZnaEnergy/F");
    fTreeEvent->Branch("fZpaEnergy", &fZpaEnergy, "fZpaEnergy/F");
    fTreeEvent->Branch("fZem1Energy", &fZem1Energy, "fZem1Energy/F");
    fTreeEvent->Branch("fZem2Energy", &fZem2Energy, "fZem2Energy/F");
    fTreeEvent->Branch("fZnaTower", &fZnaTower, "fZnaTower/F");
    fTreeEvent->Branch("fZncTower", &fZncTower, "fZncTower/F");
    fTreeEvent->Branch("fZpaTower", &fZpaTower, "fZpaTower/F");
    fTreeEvent->Branch("fZpcTower", &fZpcTower, "fZpcTower/F");
    fTreeEvent->Branch("fZnaFired", &fZnaFired, "fZnaFired/O");    //Booleans for Event Selection
    fTreeEvent->Branch("fZncFired", &fZncFired, "fZncFired/O");    //Booleans for Event Selection
    fTreeEvent->Branch("fZpaFired", &fZpaFired, "fZpaFired/O");    //Booleans for Event Selection
    fTreeEvent->Branch("fZpcFired", &fZpcFired, "fZpcFired/O");    //Booleans for Event Selection
    
    //Automatic Loop for linking directly to AliMultInput
    for( Long_t iVar=0; iVar<fInput->GetNVariables(); iVar++){
        if( !fInput->GetVariable(iVar)->IsInteger()  ){
            fTreeEvent->Branch(fInput->GetVariable(iVar)->GetName(), &fInput->GetVariable(iVar)->GetRValue(), Form("%s/F",fInput->GetVariable(iVar)->GetName()) );
        }else{
            fTreeEvent->Branch(fInput->GetVariable(iVar)->GetName(), &fInput->GetVariable(iVar)->GetRValueInteger(), Form("%s/I",fInput->GetVariable(iVar)->GetName()) );
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

    //Debugging / Memory usage tests
    //gObjectTable->Print();
    
    AliESDEvent *lESDevent = 0x0;

    //Zero all booleans, etc: safe initialization per event
    fEvSel_HasAtLeastSPDVertex    = kFALSE;
    fEvSel_VtxZCut                = kFALSE;
    fEvSel_IsNotPileup            = kFALSE;
    fEvSel_IsNotPileupMV          = kFALSE;
    fEvSel_IsNotPileupInMultBins  = kFALSE;
    fEvSel_Triggered              = kFALSE;
	fEvSel_nTracklets   = -1;
    //fnSPDClusters = -1;
    fnSPDClusters0 = -1;
    fnSPDClusters1 = -1;
    fEvSel_nContributors = -1;
    fEvSel_nSPDPrimVertices = -1;
    
    fEvSel_distZ = -100;
    fEvSel_VtxZ = -100;
    
    Float_t multADA =0;
	Float_t multADC =0;
	Float_t multAD =0;
	
    // Connect to the InputEvent
    // Appropriate for ESD analysis ..

    lESDevent = dynamic_cast<AliESDEvent*>( InputEvent() );
    if (!lESDevent) {
        AliWarning("ERROR: lESDevent not available \n");
        return;
    }
    
    //Get VZERO Information for multiplicity later
    AliVVZERO* esdV0 = lESDevent->GetVZEROData();
    if (!esdV0) {
        AliError("AliVVZERO not available");
        return;
    }
    //Get AD Multiplicity Information
    AliESDAD *fesdAD = lESDevent->GetADData();
    if(!fesdAD){
        AliWarning("ERROR:fesdAD not available\n");
        
    }
    
    
    fRunNumber = lESDevent->GetRunNumber();
    Double_t lMagneticField = -10;
    lMagneticField = lESDevent->GetMagneticField( );

    //------------------------------------------------
    // Physics Selection
    //------------------------------------------------

    fHistEventCounter->Fill(0.5);

    UInt_t maskIsSelected = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
    Bool_t isSelected = 0;
    isSelected = (maskIsSelected & AliVEvent::kMB) == AliVEvent::kMB;
    fEvSel_Triggered = isSelected;

    // FIXME Tracklets versus clusters cut
    /*
    //Tracklets vs Clusters cut via AliAnalysisUtils
    if ( fkApplyTrackletsVsClustersCut && (! fkSkipEventSelection ) ) {
        if( fUtils->IsSPDClusterVsTrackletBG( lESDevent ) ) {
            PostData(1, fListHist);
            PostData(2, fTreeEvent);
            PostData(3, fTreeV0);
            PostData(4, fTreeCascade);
            return;
        }
    }
    */

    if (  fEvSel_Triggered  ) fHistEventCounter->Fill(1.5);

    //------------------------------------------------
    // Primary Vertex Requirements Section:
    //  ---> pp: has vertex, |z|<10cm
    //------------------------------------------------

    //classical Proton-proton like selection
    const AliESDVertex *lPrimaryBestESDVtx     = lESDevent->GetPrimaryVertex();
    const AliESDVertex *lPrimaryTrackingESDVtx = lESDevent->GetPrimaryVertexTracks();
    const AliESDVertex *lPrimarySPDVtx         = lESDevent->GetPrimaryVertexSPD();

    Double_t lBestPrimaryVtxPos[3]          = {-100.0, -100.0, -100.0};
    lPrimaryBestESDVtx->GetXYZ( lBestPrimaryVtxPos );

    if(! (!lPrimarySPDVtx->GetStatus() && !lPrimaryTrackingESDVtx->GetStatus()) ) {
        //Passed selection!
        fEvSel_HasAtLeastSPDVertex = kTRUE;
    }

    //Has SPD or Tracking Vertex
    if ( fEvSel_Triggered && fEvSel_HasAtLeastSPDVertex ) fHistEventCounter -> Fill(2.5);

    if(TMath::Abs(lBestPrimaryVtxPos[2]) <= 10.0 ) {
        //Passed selection!
        fEvSel_VtxZCut = kTRUE;
    }
    fEvSel_VtxZ = lBestPrimaryVtxPos[2] ; //Set

    //Fill Event selected counter
    if ( fEvSel_Triggered && fEvSel_HasAtLeastSPDVertex && fEvSel_VtxZCut ) fHistEventCounter -> Fill(3.5);

    //------------------------------------------------
    // Check if this isn't pileup
    //------------------------------------------------

    if( !lESDevent->IsPileupFromSPD()           ) fEvSel_IsNotPileup           = kTRUE;
    if( !lESDevent->IsPileupFromSPDInMultBins() ) fEvSel_IsNotPileupInMultBins = kTRUE;

    //Acquire information to compute residual pileup
    fEvSel_nSPDPrimVertices = lESDevent->GetNumberOfPileupVerticesSPD();

    //Long_t lNcontributorsSPDvtx = lPrimarySPDVtx -> GetNContributors();
    Long_t lNcontributorsSecondLargest = -1;
    Long_t lIndexSecondLargest = -1;
    //Look for the two events with the largest numbers of contributors...
    for(Int_t i=0; i<fEvSel_nSPDPrimVertices; i++) {
        const AliESDVertex* pv=lESDevent -> GetPileupVertexSPD(i);
        if( pv->GetNContributors() > lNcontributorsSecondLargest ) {
            lNcontributorsSecondLargest = pv->GetNContributors();
            lIndexSecondLargest = i;
        }
    }
    fEvSel_nContributors = lPrimaryBestESDVtx -> GetNContributors();
    if( fEvSel_nSPDPrimVertices > 0 && lIndexSecondLargest > -1) {
        const AliESDVertex* largestpv=lESDevent ->GetPileupVertexSPD(lIndexSecondLargest);
        fEvSel_distZ = lPrimarySPDVtx->GetZ() - largestpv->GetZ();
    }

    //First implementation of pileup from multi-vertexer (simple use of analysis utils)
    //if ( !fUtils->IsPileUpMV( lESDevent ) ) fEvSel_IsNotPileupMV = kTRUE;

    //Fill Event isn't pileup counter
    if ( fEvSel_Triggered && fEvSel_HasAtLeastSPDVertex && fEvSel_VtxZCut && fEvSel_IsNotPileupInMultBins ) fHistEventCounter -> Fill(4.5);
	
	
	//------------------------------------------------
    // Multiplicity Information from AD
    //------------------------------------------------
    Float_t fMultiplicityAD[16];
    Float_t fMultiplicityADA[8];
    Float_t fMultiplicityADC[8];
    multADA =0;
    multADC =0;
    multAD =0;
    
    if (fesdAD) {
        //Get Multiplicity info per AD 16 channel: C-side : 0-7, A-side 8-15
        for (Int_t i=0;i<8; i++)
        {
            fMultiplicityAD[i]= fesdAD->GetMultiplicity(i);
            fMultiplicityADA[i]= fesdAD->GetMultiplicityADA(i);
            multADA += fMultiplicityADA[i];
            multAD += fMultiplicityAD[i];
        }
        for (Int_t i=8; i<16; i++)
        {	fMultiplicityAD[i]= fesdAD->GetMultiplicity(i);
            fMultiplicityADC[i]=fesdAD->GetMultiplicityADC(i-8);
            multADC += fMultiplicityADC[i];
            multAD+=fMultiplicityAD[i];
        }
    }
    //------------------------------------------------
    // Multiplicity Information Acquistion
    //------------------------------------------------

    //Standard GetReferenceMultiplicity Estimator (0.5 and 0.8)
    fRefMultEta5 = fESDtrackCuts->GetReferenceMultiplicity(lESDevent, AliESDtrackCuts::kTrackletsITSTPC,0.5);
    fRefMultEta8 = fESDtrackCuts->GetReferenceMultiplicity(lESDevent, AliESDtrackCuts::kTrackletsITSTPC,0.8);
    
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
    for(Int_t iCh = 48; iCh < 64; iCh++){
      Double_t mult = esdV0->GetMultiplicity(iCh); 
      multV0Apartial += mult; 
    }
    for(Int_t iCh = 0; iCh < 16; iCh++){ 
      Double_t mult = esdV0->GetMultiplicity(iCh); 
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
    for(Int_t iCh = 0; iCh < 8; iCh++){ 
      Double_t mult = esdV0->GetMultiplicity(iCh); 
      multV0C1 += mult; 
    }    
    for(Int_t iCh = 8; iCh < 16; iCh++){ 
      Double_t mult = esdV0->GetMultiplicity(iCh); 
      multV0C2 += mult; 
    }    
    for(Int_t iCh = 16; iCh < 24; iCh++){ 
      Double_t mult = esdV0->GetMultiplicity(iCh); 
      multV0C3 += mult; 
    }    
    for(Int_t iCh = 24; iCh < 32; iCh++){ 
      Double_t mult = esdV0->GetMultiplicity(iCh); 
      multV0C4 += mult; 
    }    

    for(Int_t iCh = 32; iCh < 40; iCh++){ 
      Double_t mult = esdV0->GetMultiplicity(iCh); 
      multV0A1 += mult; 
    }    
    for(Int_t iCh = 40; iCh < 48; iCh++){ 
      Double_t mult = esdV0->GetMultiplicity(iCh); 
      multV0A2 += mult; 
    }    
    for(Int_t iCh = 48; iCh < 56; iCh++){ 
      Double_t mult = esdV0->GetMultiplicity(iCh); 
      multV0A3 += mult;
    }    
    for(Int_t iCh = 56; iCh < 64; iCh++){ 
      Double_t mult = esdV0->GetMultiplicity(iCh); 
      multV0A4 += mult; 
    }
    
    //Non-Equalized Signal: copy of multV0ACorr and multV0CCorr from AliCentralitySelectionTask
    //Getters for uncorrected multiplicity
    multV0A=esdV0->GetMTotV0A();
    multV0C=esdV0->GetMTotV0C();
	//charge V0
	multonlineV0A = esdV0->GetTriggerChargeA(); //charge
	multonlineV0C = esdV0->GetTriggerChargeC(); //charge
	
    //Get Z vertex position of SPD vertex (why not Tracking if available?)
    Float_t zvtx = lPrimarySPDVtx->GetZ();

    //Acquire Corrected multV0A
    multV0ACorr = AliESDUtils::GetCorrV0A(multV0A,zvtx);
    multV0CCorr = AliESDUtils::GetCorrV0C(multV0C,zvtx);

    //Copy to Event Tree for extra information
    //FIXME: no z-vertex dependence yet, these are raw amplitudes...
    //FIXME: Would anyhow require OCDB entry
    //fAmplitude_V0A = multV0A;
    //fAmplitude_V0C = multV0C;
    
    fAmplitude_V0A->SetValue(multV0A);
    fAmplitude_V0C->SetValue(multV0C);
    
    fAmplitude_OnlineV0A->SetValue(multonlineV0A); //charge
    fAmplitude_OnlineV0C->SetValue(multonlineV0C); //charge
    //fTriggerChargeV0M = triggerChargeA+triggerChargeC; //charge
  
    fAmplitude_V0Apartial->SetValue(multV0Apartial);
    fAmplitude_V0Cpartial->SetValue(multV0Cpartial);
    
    //A.T. (vertex correction for all rings?!?)
    fAmplitude_V0A1 = multV0A1;
    fAmplitude_V0A2 = multV0A2;
    fAmplitude_V0A3 = multV0A3;
    fAmplitude_V0A4 = multV0A4;
    fAmplitude_V0C1 = multV0C1;
    fAmplitude_V0C2 = multV0C2;
    fAmplitude_V0C3 = multV0C3;
    fAmplitude_V0C4 = multV0C4;
    
    //AD scintillator Data added to Event Tree
    if (fesdAD) {
		//fMultiplicity_AD =multAD;
		fMultiplicity_ADA->SetValue(multADA);
		fMultiplicity_ADC->SetValue(multADC);
	}

    // Equalized signals // From AliCentralitySelectionTask // Updated
    for(Int_t iCh = 32; iCh < 64; ++iCh) {
        Double_t mult = lESDevent->GetVZEROEqMultiplicity(iCh);
        multV0AEq += mult;
    }
    for(Int_t iCh = 0; iCh < 32; ++iCh) {
        Double_t mult = lESDevent->GetVZEROEqMultiplicity(iCh);
        multV0CEq += mult;
    }
    fAmplitude_V0AEq->SetValue(multV0AEq);
    fAmplitude_V0CEq->SetValue(multV0CEq);


    //A.T.
    fTrackCuts = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
    fNTracks    = fTrackCuts ? (Short_t)fTrackCuts->GetReferenceMultiplicity(lESDevent,kTRUE):-1;

    // ***** ZDC info
    AliESDZDC *esdZDC = lESDevent->GetESDZDC();
    Float_t CalF=0;
    if (esdZDC->AliESDZDC::TestBit(AliESDZDC::kEnergyCalibratedSignal))  CalF=1.0; //! if zdc is calibrated (in pass2)
    else CalF=8.0;

    fZncEnergy = (Float_t) (esdZDC->GetZDCN1Energy())/CalF;
    fZpcEnergy = (Float_t) (esdZDC->GetZDCP1Energy())/CalF;
    fZnaEnergy = (Float_t) (esdZDC->GetZDCN2Energy())/CalF;
    fZpaEnergy = (Float_t) (esdZDC->GetZDCP2Energy())/CalF;

    fZem1Energy = (Float_t) (esdZDC->GetZDCEMEnergy(0))/CalF;
    fZem2Energy = (Float_t) (esdZDC->GetZDCEMEnergy(1))/CalF;

    Int_t detCh_ZNA = 0; //esdZDC->GetZNATDCChannel();
    Int_t detCh_ZNC = 0; //esdZDC->GetZNCTDCChannel();
    Int_t detCh_ZPA = 0; //esdZDC->GetZPATDCChannel();
    Int_t detCh_ZPC = 0; //esdZDC->GetZPCTDCChannel();

    for (Int_t j = 0; j < 4; ++j) {
      if (esdZDC->GetZDCTDCData(detCh_ZNA,j) != 0) 	fZnaFired = kTRUE;
      if (esdZDC->GetZDCTDCData(detCh_ZNC,j) != 0) 	fZncFired = kTRUE;
      if (esdZDC->GetZDCTDCData(detCh_ZPA,j) != 0) 	fZpaFired = kTRUE;
      if (esdZDC->GetZDCTDCData(detCh_ZPC,j) != 0) 	fZpcFired = kTRUE;
    }

    const Double_t *ZNAtower = esdZDC->GetZN2TowerEnergy(); 
    const Double_t *ZNCtower = esdZDC->GetZN1TowerEnergy();
    const Double_t *ZPAtower = esdZDC->GetZP2TowerEnergy(); 
    const Double_t *ZPCtower = esdZDC->GetZP1TowerEnergy();
    if (fZnaFired) fZnaTower = ZNAtower[0];
    if (fZncFired) fZncTower = ZNCtower[0];
    if (fZpaFired) fZpaTower = ZPAtower[0];
    if (fZpcFired) fZpcTower = ZPCtower[0];

    /////////////////

    //Tracklets vs Clusters Exploratory data
    fEvSel_nTracklets     = lESDevent->GetMultiplicity()->GetNumberOfTracklets();
    fnSPDClusters->SetValueInteger(lESDevent->GetNumberOfITSClusters(0) + lESDevent->GetNumberOfITSClusters(1));

    fEvSel_nTrackletsEta10 = fESDtrackCuts->GetReferenceMultiplicity(lESDevent, AliESDtrackCuts::kTracklets, 1.0); 
    if ( fEvSel_nTrackletsEta10 >= 1 ) fEvSel_INELgtZERO = kTRUE;

    //Event-level fill
    if ( fkCalibration ){
        //Pre-filter on triggered (kMB) events for saving info
        if( !fkFilterMB || (fkFilterMB && fEvSel_Triggered) ){
            fTreeEvent->Fill() ;
        }
    }

    if ( fkAddInfo ){//Master switch for users
        //===============================================
        // Compute Percentiles
        //===============================================
        //Make sure OADB is loaded
        SetupRun( lESDevent );
        
        //===============================================
        // I/O: Create object for storing, add
        //===============================================
        
        //Evaluate Estimators from Variables
        fSelection -> Evaluate (fInput);
        
        //Determine Quantiles from calibration histogram
        TH1F *lThisCalibHisto;
        TString lThisCalibHistoName;
        Float_t lThisQuantile = -1;
        for(Long_t iEst=0; iEst<fSelection->GetNEstimators(); iEst++){
            lThisCalibHistoName = Form("hCalib_%i_%s",fRunNumber,fSelection->GetEstimator(iEst)->GetName());
            lThisCalibHisto = oadbMultSelection->GetCalibHisto( lThisCalibHistoName );
            lThisQuantile = lThisCalibHisto->GetBinContent( lThisCalibHisto->FindBin( fSelection->GetEstimator(iEst)->GetValue() ));
            
            //cleanup: discard events
            if(!IsEventSelected(lESDevent))lThisQuantile = 200;
            
            fSelection->GetEstimator(iEst)->SetPercentile(lThisQuantile);
        }
        
        //Add to AliVEvent
        if( (!(InputEvent()->FindListObject("MultSelection")) ) && !fkAttached ){
            InputEvent()->AddObject(fSelection);
            fkAttached = kTRUE;
        }else{
            //AliInfo("Already there!"); //do nothing
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
Bool_t AliMultSelectionTask::IsEventSelected( AliESDEvent *lEvent )
{
    return fMultCuts->IsEventSelected( lEvent );
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
    
    TString fileName =(Form("%s/COMMON/MULTIPLICITY/data/OADB-01.root", AliAnalysisManager::GetOADBPath()));
    AliInfo(Form("Setup Multiplicity Selection for run %d with file %s\n",fCurrentRun,fileName.Data()));
    
    AliOADBContainer *con = new AliOADBContainer("OADB");
    con->InitFromFile(fileName,"MultSel");
    
    //Get Object for this run!
    oadbMultSelection = (AliOADBMultSelection* )con->GetObject(fCurrentRun, "Default");
    
    if (!oadbMultSelection) {
        AliWarning(Form("Multiplicity OADB does not exist for run %d, using Default \n",fCurrentRun ));
        oadbMultSelection  = (AliOADBMultSelection*)(con->GetDefaultObject("oadbDefault"));
    }
    
    //Ensure Selection is possible using simple IsEventSelected rationale
    fMultCuts  = oadbMultSelection->GetEventCuts();
    fSelection = new AliMultSelection( oadbMultSelection->GetMultSelection() ); //allocate new

    //Make sure naming convention is followed!
    fSelection->SetName("MultSelection");
    return 0;
}
