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
#include "AliVZDC.h" //AD

//#include "AliCFContainer.h"
#include "AliMultiplicity.h"
#include "AliESDUtils.h"
#include "AliAnalysisUtils.h"
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
: AliAnalysisTaskSE(), fListHist(0), fTreeEvent(0),fESDtrackCuts(0), fTrackCuts(0), fUtils(0), 
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
fnTracklets(0), 
fnSPDClusters0(0),
fnSPDClusters1(0),
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
    : AliAnalysisTaskSE(name), fListHist(0), fTreeEvent(0), fESDtrackCuts(0), fTrackCuts(0), fUtils(0), 
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
fnTracklets(0), 
fnSPDClusters0(0),
fnSPDClusters1(0),
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
    if (fESDtrackCuts) {
      delete fESDtrackCuts;
      fESDtrackCuts = 0x0;
    }
    if (fUtils) {
      delete fUtils;
      fUtils = 0x0;
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
    
    fnTracklets = new AliMultVariable("fnTracklets");
    fnTracklets ->SetIsInteger( kTRUE );
    fRefMultEta5 = new AliMultVariable("fRefMultEta5");
    fRefMultEta5 ->SetIsInteger( kTRUE );
    fRefMultEta8 = new AliMultVariable("fRefMultEta8");
    fRefMultEta8 ->SetIsInteger( kTRUE );
    
    //Add to AliMultInput Object, will later bind to TTree object in a loop
    fInput->AddVariable( fAmplitude_V0A );
    fInput->AddVariable( fAmplitude_V0C );
    fInput->AddVariable( fAmplitude_V0Apartial );
    fInput->AddVariable( fAmplitude_V0Cpartial );
    fInput->AddVariable( fAmplitude_V0AEq );
    fInput->AddVariable( fAmplitude_V0CEq );
    fInput->AddVariable( fAmplitude_OnlineV0A );
    fInput->AddVariable( fAmplitude_OnlineV0C );
    fInput->AddVariable( fnSPDClusters );
    fInput->AddVariable( fnTracklets );
    fInput->AddVariable( fRefMultEta5 );
    fInput->AddVariable( fRefMultEta8 );
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
    fTreeEvent->Branch("fEvSel_VtxZ", &fEvSel_VtxZ, "fEvSel_VtxZ/F");
    
    //A.T. FIXME change into AliMultVariable
    //fTreeEvent->Branch("fnSPDClusters0", &nSPDClusters0, "fnSPDClusters0/I");
    //fTreeEvent->Branch("fnSPDClusters1", &nSPDClusters1, "fnSPDClusters1/I");
    fTreeEvent->Branch("fNTracks",      &fNTracks, "fNTracks/I");
    
    //A.T. FIXME change into AliMultVariable
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
    fnSPDClusters0 = -1;
    fnSPDClusters1 = -1;
    fEvSel_VtxZ = -100;

    Float_t multADA =0;
    Float_t multADC =0;
    Float_t multAD =0;

    // Connect to the InputEvent
    // Appropriate for ESD analysis ..

    lVevent = dynamic_cast<AliVEvent*>( InputEvent() );
    if (!lVevent) {
        AliWarning("ERROR: ESD / AOD event not available \n");
        return;
    }

    //Get VZERO Information for multiplicity later
    AliVVZERO* lVV0 = lVevent->GetVZEROData();
    if (!lVV0) {
        AliError("AliVVZERO not available");
        return;
    }
    //Get AD Multiplicity Information
    AliVAD *lVAD = lVevent->GetADData();
    if(!lVAD) {
        AliWarning("ERROR:lVAD not available\n");

    }

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
    
    fEvSel_Triggered                 = IsSelectedTrigger                   (lVevent);
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

    if(TMath::Abs(lBestPrimaryVtxPos[2]) <= 10.0 ) {
        //FIXME Passed default 10.0cm selection!
        fEvSel_VtxZCut = kTRUE;
    }
    fEvSel_VtxZ = lBestPrimaryVtxPos[2] ; //Set for later use 
    
    //===============================================
    // End Event Selection Variables Section
    //===============================================
    
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
        for (Int_t i=0;i<8; i++)
        {
            fMultiplicityAD[i]= lVAD->GetMultiplicity(i);
            fMultiplicityADA[i]= lVAD->GetMultiplicityADA(i);
            multADA += fMultiplicityADA[i];
            multAD += fMultiplicityAD[i];
        }
        for (Int_t i=8; i<16; i++)
        {	fMultiplicityAD[i]= lVAD->GetMultiplicity(i);
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
    for(Int_t iCh = 48; iCh < 64; iCh++){
      Double_t mult = lVV0->GetMultiplicity(iCh); 
      multV0Apartial += mult; 
    }
    for(Int_t iCh = 0; iCh < 16; iCh++){ 
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
    for(Int_t iCh = 0; iCh < 8; iCh++){ 
      Double_t mult = lVV0->GetMultiplicity(iCh); 
      multV0C1 += mult; 
    }    
    for(Int_t iCh = 8; iCh < 16; iCh++){ 
      Double_t mult = lVV0->GetMultiplicity(iCh); 
      multV0C2 += mult; 
    }    
    for(Int_t iCh = 16; iCh < 24; iCh++){ 
      Double_t mult = lVV0->GetMultiplicity(iCh); 
      multV0C3 += mult; 
    }    
    for(Int_t iCh = 24; iCh < 32; iCh++){ 
      Double_t mult = lVV0->GetMultiplicity(iCh); 
      multV0C4 += mult; 
    }    

    for(Int_t iCh = 32; iCh < 40; iCh++){ 
      Double_t mult = lVV0->GetMultiplicity(iCh); 
      multV0A1 += mult; 
    }    
    for(Int_t iCh = 40; iCh < 48; iCh++){ 
      Double_t mult = lVV0->GetMultiplicity(iCh); 
      multV0A2 += mult; 
    }    
    for(Int_t iCh = 48; iCh < 56; iCh++){ 
      Double_t mult = lVV0->GetMultiplicity(iCh); 
      multV0A3 += mult;
    }    
    for(Int_t iCh = 56; iCh < 64; iCh++){ 
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
    
    fAmplitude_OnlineV0A->SetValue(multonlineV0A);
    fAmplitude_OnlineV0C->SetValue(multonlineV0C);

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

    //===============================================
    //This part requires separation of AOD and ESD
    //===============================================
    
    //Setting variables to non-sense values
    fRefMultEta5 -> SetValueInteger ( -5 ); //not acquired 
    fRefMultEta8 -> SetValueInteger ( -5 ); //not acquired 
    fNTracks         = -10;

    if (lVevent->InheritsFrom("AliESDEvent")) {
        AliESDEvent *esdevent = dynamic_cast<AliESDEvent *>(lVevent);
        //Standard GetReferenceMultiplicity Estimator (0.5 and 0.8)
        fRefMultEta5 -> SetValueInteger ( fESDtrackCuts->GetReferenceMultiplicity(esdevent, AliESDtrackCuts::kTrackletsITSTPC,0.5) );
        fRefMultEta8 -> SetValueInteger ( fESDtrackCuts->GetReferenceMultiplicity(esdevent, AliESDtrackCuts::kTrackletsITSTPC,0.8) );
        
        //Use fallback in case of return value of -3 or -4
        //This is what will happen in AODs: -3 will use fallback
        //HOWEVER: -4 will not. Inconsistency requires use of "HasNoInconsistentSPDandTrackVertices"!
        if ( fRefMultEta5 -> GetValueInteger() < -2 ){
            fRefMultEta5 -> SetValueInteger ( fESDtrackCuts->GetReferenceMultiplicity(esdevent, AliESDtrackCuts::kTracklets,0.5) );
        }
        if ( fRefMultEta8 -> GetValueInteger() < -2 ){
            fRefMultEta8 -> SetValueInteger ( fESDtrackCuts->GetReferenceMultiplicity(esdevent, AliESDtrackCuts::kTracklets,0.8) );
        }

        //A.T.
        fTrackCuts = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
        fNTracks    = fTrackCuts ? (Short_t)fTrackCuts->GetReferenceMultiplicity(esdevent,kTRUE):-1;

        // ***** ZDC info
        AliESDZDC *lESDZDC = esdevent->GetESDZDC();
        Float_t CalF=0;
        if (lESDZDC->AliESDZDC::TestBit(AliESDZDC::kEnergyCalibratedSignal))  CalF=1.0; //! if zdc is calibrated (in pass2)
        else CalF=8.0;

        fZncEnergy = (Float_t) (lESDZDC->GetZDCN1Energy())/CalF;
        fZpcEnergy = (Float_t) (lESDZDC->GetZDCP1Energy())/CalF;
        fZnaEnergy = (Float_t) (lESDZDC->GetZDCN2Energy())/CalF;
        fZpaEnergy = (Float_t) (lESDZDC->GetZDCP2Energy())/CalF;

        fZem1Energy = (Float_t) (lESDZDC->GetZDCEMEnergy(0))/CalF;
        fZem2Energy = (Float_t) (lESDZDC->GetZDCEMEnergy(1))/CalF;

        Int_t detCh_ZNA = 0; //lESDZDC->GetZNATDCChannel();
        Int_t detCh_ZNC = 0; //lESDZDC->GetZNCTDCChannel();
        Int_t detCh_ZPA = 0; //lESDZDC->GetZPATDCChannel();
        Int_t detCh_ZPC = 0; //lESDZDC->GetZPCTDCChannel();

        for (Int_t j = 0; j < 4; ++j) {
            if (lESDZDC->GetZDCTDCData(detCh_ZNA,j) != 0) 	fZnaFired = kTRUE;
            if (lESDZDC->GetZDCTDCData(detCh_ZNC,j) != 0) 	fZncFired = kTRUE;
            if (lESDZDC->GetZDCTDCData(detCh_ZPA,j) != 0) 	fZpaFired = kTRUE;
            if (lESDZDC->GetZDCTDCData(detCh_ZPC,j) != 0) 	fZpcFired = kTRUE;
        }

        const Double_t *ZNAtower = lESDZDC->GetZN2TowerEnergy();
        const Double_t *ZNCtower = lESDZDC->GetZN1TowerEnergy();
        const Double_t *ZPAtower = lESDZDC->GetZP2TowerEnergy();
        const Double_t *ZPCtower = lESDZDC->GetZP1TowerEnergy();
        if (fZnaFired) fZnaTower = ZNAtower[0];
        if (fZncFired) fZncTower = ZNCtower[0];
        if (fZpaFired) fZpaTower = ZPAtower[0];
        if (fZpcFired) fZpcTower = ZPCtower[0];

    } else if (lVevent->InheritsFrom("AliAODEvent")) {
        AliAODEvent *aodevent = dynamic_cast<AliAODEvent *>(lVevent);
        AliAODHeader * header = dynamic_cast<AliAODHeader*>(aodevent->GetHeader());
        fRefMultEta5 -> SetValueInteger ( header->GetRefMultiplicityComb05() );
        fRefMultEta8 -> SetValueInteger ( header->GetRefMultiplicityComb08() );

        //FIXME: get ZDC information in AOD in a fully consistent way
    }

    //===============================================
    // End part which requires AOD/ESD separation
    //===============================================

    //Event-level fill
    if ( fkCalibration ) {
        //Pre-filter on triggered (kMB) events for saving info
        if( !fkFilterMB || (fkFilterMB && fEvSel_Triggered) ) {
            fTreeEvent->Fill() ;
        }
    }

    if ( fkAddInfo ) { //Master switch for users
        //===============================================
        // Compute Percentiles
        //===============================================
        //Make sure OADB is loaded
        SetupRun( lVevent );

        //===============================================
        // I/O: Create object for storing, add
        //===============================================
        
        //Evaluate Estimators from Variables
        fSelection -> Evaluate (fInput);

        //Determine Quantiles from calibration histogram
        TH1F *lThisCalibHisto;
        TString lThisCalibHistoName;
        Float_t lThisQuantile = -1;
        for(Long_t iEst=0; iEst<fSelection->GetNEstimators(); iEst++) {
            lThisCalibHistoName = Form("hCalib_%i_%s",fRunNumber,fSelection->GetEstimator(iEst)->GetName());
            lThisCalibHisto = oadbMultSelection->GetCalibHisto( lThisCalibHistoName );
            lThisQuantile = lThisCalibHisto->GetBinContent( lThisCalibHisto->FindBin( fSelection->GetEstimator(iEst)->GetValue() ));

            //cleanup: discard events according to criteria stored in OADB object
            //Check Selections as they are in the fMultSelectionCuts Object
            if( fMultCuts->GetTriggerCut()    && ! fEvSel_Triggered           ) lThisQuantile = 200;
            if( fMultCuts->GetINELgtZEROCut() && ! fEvSel_INELgtZERO          ) lThisQuantile = 201;
            if( TMath::Abs(fEvSel_VtxZ)          > fMultCuts->GetVzCut()      ) lThisQuantile = 202;
            if( fMultCuts->GetRejectPileupInMultBinsCut() && ! fEvSel_IsNotPileupInMultBins      ) lThisQuantile = 203;
            if( fMultCuts->GetVertexConsistencyCut()      && ! fEvSel_HasNoInconsistentVertices  ) lThisQuantile = 204;
            if( fMultCuts->GetTrackletsVsClustersCut()    && ! fEvSel_PassesTrackletVsCluster    ) lThisQuantile = 205;

            fSelection->GetEstimator(iEst)->SetPercentile(lThisQuantile);
        }

        //Add to AliVEvent
        if( (!(InputEvent()->FindListObject("MultSelection")) ) && !fkAttached ) {
            InputEvent()->AddObject(fSelection);
            fkAttached = kTRUE;
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


//______________________________________________________________________
Bool_t AliMultSelectionTask::IsSelectedTrigger(AliVEvent* event, AliVEvent::EOfflineTriggerTypes trigType)
// Function to check for a specific trigger class available in AliVEvent (default AliVEvent::kMB)
{
    //Code to reject events that aren't trigType
    UInt_t maskIsSelected = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
    Bool_t isSelected = 0;
    isSelected = (maskIsSelected & trigType) == trigType;
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
        if( lPrimaryVtxTracks->GetStatus() && lPrimaryVtxSPD->GetStatus() ){
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
        if ( esdevent->IsPileupFromSPDInMultBins() == kTRUE ) lReturnValue = kFALSE;
    }
    else if (event->InheritsFrom("AliAODEvent")) {
        AliAODEvent *aodevent = dynamic_cast<AliAODEvent *>(event);
        if (!aodevent) return kFALSE;
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

