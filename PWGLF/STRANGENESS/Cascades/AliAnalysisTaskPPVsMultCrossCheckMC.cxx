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


using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskPPVsMultCrossCheckMC)

AliAnalysisTaskPPVsMultCrossCheckMC::AliAnalysisTaskPPVsMultCrossCheckMC()
: AliAnalysisTaskSE(), fListHist(0), fPIDResponse(0), fESDtrackCuts(0), fPPVsMultUtils(0), fUtils(0),
fHistEventCounter(0),fHistV0M_DataSelection(0), fHistV0M_MCSelection(0), fHistV0MVsMidRapidityTrue(0),
fHistV0MTrueVsMidRapidityTrue(0)
{
    //Empty constructor (not to be used, always pass name...)
    for(Int_t ih=0; ih<9; ih++){
        fHistPt_Generated[ih] = 0x0;
        fHistPt_DataSelection[ih] = 0x0;
        fHistPt_MCSelection[ih] = 0x0;
        fHistPtVsV0M_Generated[ih] = 0x0;
        fHistPtVsV0M_DataSelection[ih] = 0x0;
        fHistPtVsV0M_MCSelection[ih] = 0x0;
        fHistPtVsV0MTrue_Generated[ih] = 0x0;
        fHistPtVsV0MTrue_DataSelection[ih] = 0x0;
        fHistPtVsV0MTrue_MCSelection[ih] = 0x0;
    }
}

AliAnalysisTaskPPVsMultCrossCheckMC::AliAnalysisTaskPPVsMultCrossCheckMC(const char *name)
: AliAnalysisTaskSE(name), fListHist(0), fPIDResponse(0), fESDtrackCuts(0), fPPVsMultUtils(0), fUtils(0),
fHistEventCounter(0),fHistV0M_DataSelection(0), fHistV0M_MCSelection(0), fHistV0MVsMidRapidityTrue(0),
fHistV0MTrueVsMidRapidityTrue(0)
{
    for(Int_t ih=0; ih<9; ih++){
        fHistPt_Generated[ih] = 0x0;
        fHistPt_DataSelection[ih] = 0x0;
        fHistPt_MCSelection[ih] = 0x0;
        fHistPtVsV0M_Generated[ih] = 0x0;
        fHistPtVsV0M_DataSelection[ih] = 0x0;
        fHistPtVsV0M_MCSelection[ih] = 0x0;
        fHistPtVsV0MTrue_Generated[ih] = 0x0;
        fHistPtVsV0MTrue_DataSelection[ih] = 0x0;
        fHistPtVsV0MTrue_MCSelection[ih] = 0x0;
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
    
    AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
    fPIDResponse = inputHandler->GetPIDResponse();
    inputHandler->SetNeedField();
    
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
    
    //Main Output: Histograms
    if(! fHistV0M_DataSelection ) {
        fHistV0M_DataSelection = new TH1F("fHistV0M_DataSelection","",100,0,100);
        fListHist->Add(fHistV0M_DataSelection);
    }
    if(! fHistV0M_MCSelection ) {
        fHistV0M_MCSelection = new TH1F("fHistV0M_MCSelection","",100,0,100);
        fListHist->Add(fHistV0M_MCSelection);
    }

    //Correlation between mid-rapidity and forward
    if(! fHistV0MVsMidRapidityTrue ) {
        fHistV0MVsMidRapidityTrue = new TH2F("fHistV0MVsMidRapidityTrue","",100,0,100,1000,0,1000);
        fListHist->Add(fHistV0MVsMidRapidityTrue);
    }
    //Correlation between mid-rapidity and forward
    if(! fHistV0MTrueVsMidRapidityTrue ) {
        fHistV0MTrueVsMidRapidityTrue = new TH2F("fHistV0MTrueVsMidRapidityTrue","",1000,0,1000,1000,0,1000);
        fListHist->Add(fHistV0MTrueVsMidRapidityTrue);
    }
    
    //Main Output: Histograms
    for(Int_t ih=0; ih<9; ih++){
        if(! fHistPt_Generated[ih] ) {
            fHistPt_Generated[ih] = new TH1F(Form("fHistPt_Generated_%s",lPartNames[ih].Data()),    "Generated;p_{T} (GeV/c)",200,0,20);
            fListHist->Add(fHistPt_Generated[ih]);
        }
    }
    for(Int_t ih=0; ih<9; ih++){
        if(! fHistPt_DataSelection[ih] ) {
            fHistPt_DataSelection[ih] = new TH1F(Form("fHistPt_DataSelection_%s",lPartNames[ih].Data()),    "Generated;p_{T} (GeV/c)",200,0,20);
            fListHist->Add(fHistPt_DataSelection[ih]);
        }
    }
    for(Int_t ih=0; ih<9; ih++){
        if(! fHistPt_MCSelection[ih] ) {
            fHistPt_MCSelection[ih] = new TH1F(Form("fHistPt_MCSelection_%s",lPartNames[ih].Data()),    "Generated;p_{T} (GeV/c)",200,0,20);
            fListHist->Add(fHistPt_MCSelection[ih]);
        }
    }
    
    //2-Dimensional Histograms
    for(Int_t ih=0; ih<9; ih++){
        if(! fHistPtVsV0M_Generated[ih] ) {
            fHistPtVsV0M_Generated[ih] = new TH2F(Form("fHistPtVsV0M_Generated_%s",lPartNames[ih].Data()),    "Generated;p_{T} (GeV/c)",200,0,20,100,0,100);
            fListHist->Add(fHistPtVsV0M_Generated[ih]);
        }
    }
    for(Int_t ih=0; ih<9; ih++){
        if(! fHistPtVsV0M_DataSelection[ih] ) {
            fHistPtVsV0M_DataSelection[ih] = new TH2F(Form("fHistPtVsV0M_DataSelection_%s",lPartNames[ih].Data()),    "Generated;p_{T} (GeV/c)",200,0,20,100,0,100);
            fListHist->Add(fHistPtVsV0M_DataSelection[ih]);
        }
    }
    for(Int_t ih=0; ih<9; ih++){
        if(! fHistPtVsV0M_MCSelection[ih] ) {
            fHistPtVsV0M_MCSelection[ih] = new TH2F(Form("fHistPtVsV0M_MCSelection_%s",lPartNames[ih].Data()),    "Generated;p_{T} (GeV/c)",200,0,20,100,0,100);
            fListHist->Add(fHistPtVsV0M_MCSelection[ih]);
        }
    }
    
    //2-Dimensional Histograms with True V0M Multiplicity
    for(Int_t ih=0; ih<9; ih++){
        if(! fHistPtVsV0MTrue_Generated[ih] ) {
            fHistPtVsV0MTrue_Generated[ih] = new TH2F(Form("fHistPtVsV0MTrue_Generated_%s",lPartNames[ih].Data()),    "Generated;p_{T} (GeV/c)",200,0,20,1000,0,1000);
            fListHist->Add(fHistPtVsV0MTrue_Generated[ih]);
        }
    }
    for(Int_t ih=0; ih<9; ih++){
        if(! fHistPtVsV0MTrue_DataSelection[ih] ) {
            fHistPtVsV0MTrue_DataSelection[ih] = new TH2F(Form("fHistPtVsV0MTrue_DataSelection_%s",lPartNames[ih].Data()),    "Generated;p_{T} (GeV/c)",200,0,20,1000,0,1000);
            fListHist->Add(fHistPtVsV0MTrue_DataSelection[ih]);
        }
    }
    for(Int_t ih=0; ih<9; ih++){
        if(! fHistPtVsV0MTrue_MCSelection[ih] ) {
            fHistPtVsV0MTrue_MCSelection[ih] = new TH2F(Form("fHistPtVsV0MTrue_MCSelection_%s",lPartNames[ih].Data()),    "Generated;p_{T} (GeV/c)",200,0,20,1000,0,1000);
            fListHist->Add(fHistPtVsV0MTrue_MCSelection[ih]);
        }
    }
    
    //2D Correlation plot between multiplicities
    
    
    
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
    
    Bool_t lPureMonteCarlo = kFALSE;
    
    // Connect to the InputEvent
    // After these lines, we should have an ESD/AOD event + the number of V0s in it.
    
    // Appropriate for ESD analysis!
    
    lESDevent = dynamic_cast<AliESDEvent*>( InputEvent() );
    if (!lESDevent) {
        AliWarning("AliESDevent not available, going to Pure Monte Carlo Mode! \n");
        lPureMonteCarlo = kTRUE;
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
    TArrayF mcPrimaryVtx;
    AliGenEventHeader* mcHeader=lMCevent->GenEventHeader();
    if(!mcHeader) return;
    mcHeader->PrimaryVertex(mcPrimaryVtx);
    lVertexZMC = mcPrimaryVtx.At(2);
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
        UInt_t maskIsSelected = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
        isSelected = (maskIsSelected & AliVEvent::kMB) == AliVEvent::kMB;
    }
    //Physics Selection
    lEvSel_Triggered = isSelected;
    
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
        if( TMath::Abs(geta) < 1.0 ) lNchEta10++;
        if( TMath::Abs(geta) < 1.0 ) lEvSel_INELgtZEROStackPrimaries = kTRUE;
        if( 2.8 < geta && geta < 5.1 ) lNchVZEROA++;
        if(-3.7 < geta && geta <-1.7 ) lNchVZEROC++;
    }//End of loop on tracks
    //----- End Loop on Stack ------------------------------------------------------------
    
    //Standard GetReferenceMultiplicity Estimator (0.5 and 0.8)
    
    Float_t fCentrality_V0M   = -100;
    if ( ! lPureMonteCarlo ) fCentrality_V0M   = fPPVsMultUtils -> GetMultiplicityPercentile(lESDevent, "V0M"   );
    Float_t fCentrality_V0MUnselected = -100;
    if ( ! lPureMonteCarlo ) fCentrality_V0MUnselected = fPPVsMultUtils -> GetMultiplicityPercentile(lESDevent, "V0M" , kFALSE );
    
    //------------------------------------------------
    // Get All Conditionals
    //------------------------------------------------
    Bool_t lIsINELgtZEROtracklets    = kFALSE;
    Bool_t lIsAcceptedVertexPosition = kFALSE;
    Bool_t lIsNotPileupInMultBins    = kFALSE;
    Bool_t lConsistentVertices       = kFALSE;

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
    }
    
    //Monte-Carlo Specifics
    //Alpha: INELgtZERO in the MC stack
    // (already exists, called "fEvSel_INELgtZEROStackPrimaries")
    //Beta: Vertex position
    Bool_t lIsAcceptedVertexPositionMC = (TMath::Abs(lVertexZMC)<10.0); //true if within desired range
    
    //Merge all conditionals
    
    Bool_t lDataSelection = ( lEvSel_Triggered && lIsINELgtZEROtracklets && lIsAcceptedVertexPosition && lIsNotPileupInMultBins && lConsistentVertices );
    
    Bool_t lMCSelection   = ( lEvSel_Triggered && lEvSel_INELgtZEROStackPrimaries && lIsAcceptedVertexPositionMC );
    
    //------------------------------------------------
    // Fill Event Counters
    //------------------------------------------------

    //Basics: All Processed
    fHistEventCounter->Fill(0.5);
    if( lDataSelection ) fHistEventCounter -> Fill(1.5);
    if( lMCSelection   ) fHistEventCounter -> Fill(2.5);
    
    if( lDataSelection ) fHistV0M_DataSelection -> Fill( fCentrality_V0MUnselected );
    if( lMCSelection   ) fHistV0M_MCSelection   -> Fill( fCentrality_V0MUnselected );

    if( lMCSelection ) {
        fHistV0MVsMidRapidityTrue->Fill( fCentrality_V0MUnselected, lNchEta5 );
        fHistV0MTrueVsMidRapidityTrue->Fill( lNchVZEROA+lNchVZEROC, lNchEta5 );
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
                fHistPtVsV0MTrue_Generated[ih] -> Fill(lThisPt,lNchVZEROA+lNchVZEROC);
                if( lDataSelection ){
                    fHistPt_DataSelection     [ih] -> Fill(lThisPt);
                    fHistPtVsV0M_DataSelection[ih] -> Fill(lThisPt,fCentrality_V0MUnselected);
                    fHistPtVsV0MTrue_DataSelection[ih] -> Fill(lThisPt,lNchVZEROA+lNchVZEROC);
                }
                if( lMCSelection   ){
                    fHistPt_MCSelection       [ih] -> Fill(lThisPt);
                    fHistPtVsV0M_MCSelection  [ih] -> Fill(lThisPt,fCentrality_V0MUnselected);
                    fHistPtVsV0MTrue_MCSelection  [ih] -> Fill(lThisPt,lNchVZEROA+lNchVZEROC);
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