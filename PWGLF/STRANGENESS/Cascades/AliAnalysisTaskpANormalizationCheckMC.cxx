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
//        V0A percentiles (cuts applied: all analysis cuts)
//
//    (B) Spectra in the true INEL>0 class as a function of V0A percentiles
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
#include "AliHeader.h"
#include "AliGenDPMjetEventHeader.h"
#include "AliAnalysisTaskpANormalizationCheckMC.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskpANormalizationCheckMC)

AliAnalysisTaskpANormalizationCheckMC::AliAnalysisTaskpANormalizationCheckMC()
: AliAnalysisTaskSE(), fListHist(0),fUtils(0),
fHistEventCounter(0),
fHistV0A_DataSelection(0), fHistV0A_MCSelection(0)
{
    //Empty constructor (not to be used, always pass name...)
    for(Int_t ih=0; ih<9; ih++){
        fHistPt_Generated[ih] = 0x0;
        fHistPt_DataSelection[ih] = 0x0;
        fHistPt_MCSelection[ih] = 0x0;
        fHistPtVsV0A_Generated[ih] = 0x0;
        fHistPtVsV0A_DataSelection[ih] = 0x0;
        fHistPtVsV0A_MCSelection[ih] = 0x0;
    }
}

AliAnalysisTaskpANormalizationCheckMC::AliAnalysisTaskpANormalizationCheckMC(const char *name)
: AliAnalysisTaskSE(name), fListHist(0),fUtils(0),
fHistEventCounter(0),
fHistV0A_DataSelection(0), fHistV0A_MCSelection(0)
{
    for(Int_t ih=0; ih<9; ih++){
        fHistPt_Generated[ih] = 0x0;
        fHistPt_DataSelection[ih] = 0x0;
        fHistPt_MCSelection[ih] = 0x0;
        fHistPtVsV0A_Generated[ih] = 0x0;
        fHistPtVsV0A_DataSelection[ih] = 0x0;
        fHistPtVsV0A_MCSelection[ih] = 0x0;
    }
    DefineOutput(1, TList::Class()); // Event Counter Histo
}


AliAnalysisTaskpANormalizationCheckMC::~AliAnalysisTaskpANormalizationCheckMC()
{
    //------------------------------------------------
    // DESTRUCTOR
    //------------------------------------------------
    
    if (fListHist) {
        delete fListHist;
        fListHist = 0x0;
    }
    if (fUtils) {
        delete fUtils;
        fUtils = 0x0;
    }
}

//________________________________________________________________________
void AliAnalysisTaskpANormalizationCheckMC::UserCreateOutputObjects()
{
    //------------------------------------------------
    // Particle Identification Setup
    //------------------------------------------------
    
    AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
    inputHandler->SetNeedField();
    
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
    
    //Main Output: Histograms
    if(! fHistV0A_DataSelection ) {
        fHistV0A_DataSelection = new TH1F("fHistV0A_DataSelection","",100,0,100);
        fListHist->Add(fHistV0A_DataSelection);
    }
    if(! fHistV0A_MCSelection ) {
        fHistV0A_MCSelection = new TH1F("fHistV0A_MCSelection","",100,0,100);
        fListHist->Add(fHistV0A_MCSelection);
    }
    if(! fHistV0AVsNch_DataSelection ) {
        fHistV0AVsNch_DataSelection = new TH2F("fHistV0AVsNch_DataSelection","",100,0,100,200,-0.5,199.5);
        fListHist->Add(fHistV0AVsNch_DataSelection);
    }
    if(! fHistV0AVsNch_MCSelection ) {
        fHistV0AVsNch_MCSelection = new TH2F("fHistV0AVsNch_MCSelection","",100,0,100,200,-0.5,199.5);
        fListHist->Add(fHistV0AVsNch_MCSelection);
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
        if(! fHistPtVsV0A_Generated[ih] ) {
            fHistPtVsV0A_Generated[ih] = new TH2F(Form("fHistPtVsV0A_Generated_%s",lPartNames[ih].Data()),    "Generated;p_{T} (GeV/c)",lNPtBins,0,lMaxPt,100,0,100);
            fListHist->Add(fHistPtVsV0A_Generated[ih]);
        }
    }
    for(Int_t ih=0; ih<9; ih++){
        if(! fHistPtVsV0A_DataSelection[ih] ) {
            fHistPtVsV0A_DataSelection[ih] = new TH2F(Form("fHistPtVsV0A_DataSelection_%s",lPartNames[ih].Data()),    "Generated;p_{T} (GeV/c)",lNPtBins,0,lMaxPt,100,0,100);
            fListHist->Add(fHistPtVsV0A_DataSelection[ih]);
        }
    }
    for(Int_t ih=0; ih<9; ih++){
        if(! fHistPtVsV0A_MCSelection[ih] ) {
            fHistPtVsV0A_MCSelection[ih] = new TH2F(Form("fHistPtVsV0A_MCSelection_%s",lPartNames[ih].Data()),    "Generated;p_{T} (GeV/c)",lNPtBins,0,lMaxPt,100,0,100);
            fListHist->Add(fHistPtVsV0A_MCSelection[ih]);
        }
    }
    
    //List of Histograms: Normal
    PostData(1, fListHist);
    
}// end UserCreateOutputObjects


//________________________________________________________________________
void AliAnalysisTaskpANormalizationCheckMC::UserExec(Option_t *)
{
    // Main loop
    // Called for each event
    
    AliESDEvent *lESDevent = 0x0;
    AliMCEvent  *lMCevent  = 0x0;
    AliStack    *lMCstack  = 0x0;
    
    // Connect to the InputEvent
    // After these lines, we should have an ESD/AOD event + the number of V0s in it.
    
    // Appropriate for ESD analysis!
    
    lESDevent = dynamic_cast<AliESDEvent*>( InputEvent() );
    if (!lESDevent) {
        AliWarning("AliESDevent not available, this should be in pure Monte Carlo mode... \n");
        return;
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
    if(!mcHeader){
        Printf("ERROR: Could not retrieve MC header \n");
        return;
    }
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
    lMagneticField = lESDevent->GetMagneticField( );

    //------------------------------------------------
    // Event Selection
    //------------------------------------------------
    // Copy of pA spectra analysis snippet
    //------------------------------------------------

    //Booleans to store selections
    Bool_t lbIsPS    = kTRUE;
    Bool_t lbHasPV   = kTRUE;
    Bool_t lbFstCk   = kTRUE;
    Bool_t lbVtxZ    = kTRUE;
    Bool_t lbCent    = kTRUE;
    
    UInt_t maskIsSelected = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
    Bool_t isSelected = 0;
    isSelected = (maskIsSelected & AliVEvent::kINT7) == AliVEvent::kINT7;

    //pA Selection
    if ( ! isSelected ) {
        lbIsPS = kFALSE;
    }

    //Roberto's PV selection criteria, implemented 17th April 2013
    
    /* vertex selection */
    Bool_t fHasVertex = kFALSE;
    Bool_t lWouldHaveBeenRemoved = kFALSE;
    
    const AliESDVertex *vertex = lESDevent->GetPrimaryVertexTracks();
    if (vertex->GetNContributors() < 1) {
        vertex = lESDevent->GetPrimaryVertexSPD();
        if (vertex->GetNContributors() < 1) fHasVertex = kFALSE;
        else fHasVertex = kTRUE;
        Double_t cov[6]= {0};
        vertex->GetCovarianceMatrix(cov);
        Double_t zRes = TMath::Sqrt(cov[5]);
        if (vertex->IsFromVertexerZ() && (zRes>0.25)) {
            if( fHasVertex ) lWouldHaveBeenRemoved = kTRUE;
            fHasVertex = kFALSE;
        }
    }
    else fHasVertex = kTRUE;
    
    //Is First event in chunk rejection: Still present!
    if(fHasVertex == kFALSE) {
        lbHasPV = kFALSE;
    }
    //Is First event in chunk rejection: Still present!
    if(fUtils->IsFirstEventInChunk(lESDevent)) {
        lbFstCk = kFALSE;
    }
    
    Double_t lBestPrimaryVtxPos[3]          = {-100.0, -100.0, -100.0};
    const AliESDVertex *lPrimaryBestESDVtx = lESDevent->GetPrimaryVertex();
    // get the best primary vertex available for the event
    // As done in AliCascadeVertexer, we keep the one which is the best one available.
    // between : Tracking vertex > SPD vertex > TPC vertex > default SPD vertex
    // This one will be used for next calculations (DCA essentially)
    lPrimaryBestESDVtx->GetXYZ( lBestPrimaryVtxPos );
    
    //17 April Fix: Always do primary vertex Z selection, after pA vertex selection from Roberto
    if(TMath::Abs(lBestPrimaryVtxPos[2]) > 10.0) {
        lbVtxZ = kFALSE;
    }
    
    //------------------------------------------------
    // Multiplicity Information Acquistion
    //------------------------------------------------
    
    //Only V0A Explored in this case
    //Requires AliCentrality Object to be defined
    AliCentrality* centrality;
    centrality = lESDevent->GetCentrality();
    Float_t fCentrality_V0A = ( ( Int_t ) ( centrality->GetCentralityPercentile( "V0A" ) ) );
    if (centrality->GetQuality()>1) {
        lbCent = kFALSE;
    }
    
    //MC Information: NSD + Vertex generated within |z|<10cm
    Bool_t lbMCIsNSD = IsNSD( lMCevent );
    Bool_t lbMCVtxZ = kTRUE;
    if(TMath::Abs(lVertexZMC) > 10.0) {
        lbMCVtxZ = kFALSE;
    }
    
    //pA Spectra Analysis Conditionals
    Bool_t lDataSelection = ( lbIsPS && lbHasPV && lbVtxZ && lbFstCk && lbCent );
    
    //NSD-generated Conditionals
    Bool_t lMCSelection   = ( lbIsPS && lbFstCk && lbMCVtxZ && lbCent );
    
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
    
    //------------------------------------------------
    // Fill Event Counters
    //------------------------------------------------

    //Basics: All Processed
    fHistEventCounter->Fill(0.5);
    if( lDataSelection ) fHistEventCounter -> Fill(1.5);
    if( lMCSelection   ) fHistEventCounter -> Fill(2.5);
    
    if( lDataSelection ) {
        fHistV0A_DataSelection -> Fill( fCentrality_V0A );
        fHistV0AVsNch_DataSelection -> Fill( fCentrality_V0A, lNchEta5 );
    }
    if( lMCSelection   ){
        fHistV0A_MCSelection   -> Fill( fCentrality_V0A );
        fHistV0AVsNch_MCSelection   -> Fill( fCentrality_V0A, lNchEta5 );
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
                fHistPtVsV0A_Generated[ih] -> Fill(lThisPt,fCentrality_V0A);
                if( lDataSelection ){
                    fHistPt_DataSelection     [ih] -> Fill(lThisPt);
                    fHistPtVsV0A_DataSelection[ih] -> Fill(lThisPt,fCentrality_V0A);
                }
                if( lMCSelection   ){
                    fHistPt_MCSelection       [ih] -> Fill(lThisPt);
                    fHistPtVsV0A_MCSelection  [ih] -> Fill(lThisPt,fCentrality_V0A);
                }
            }
        }
    }//End of loop on tracks
    //----- End Loop on Stack ----------------------
    
    // Post output data.
    PostData(1, fListHist);
}

//________________________________________________________________________
void AliAnalysisTaskpANormalizationCheckMC::Terminate(Option_t *)
{
    // Draw result to the screen
    // Called once at the end of the query
    
    TList *cRetrievedList = 0x0;
    cRetrievedList = (TList*)GetOutputData(1);
    if(!cRetrievedList) {
        Printf("ERROR - AliAnalysisTaskpANormalizationCheckMC : ouput data container list not available\n");
        return;
    }
    
    fHistEventCounter = dynamic_cast<TH1F*> (  cRetrievedList->FindObject("fHistEventCounter")  );
    if (!fHistEventCounter) {
        Printf("ERROR - AliAnalysisTaskpANormalizationCheckMC : fHistEventCounter not available");
        return;
    }
    
    TCanvas *canCheck = new TCanvas("AliAnalysisTaskpANormalizationCheckMC","Event Multiplicity",10,10,510,510);
    canCheck->cd(1)->SetLogy();
    
    fHistEventCounter->SetMarkerStyle(22);
    fHistEventCounter->DrawCopy("E");
}

//----------------------------------------------------------------------------
Double_t AliAnalysisTaskpANormalizationCheckMC::MyRapidity(Double_t rE, Double_t rPz) const
{
    // Local calculation for rapidity
    Double_t ReturnValue = -100;
    if( (rE-rPz+1.e-13) != 0 && (rE+rPz) != 0 ) {
        ReturnValue =  0.5*TMath::Log((rE+rPz)/(rE-rPz+1.e-13));
    }
    return ReturnValue;
}

//----------------------------------------------------------------------------
Bool_t AliAnalysisTaskpANormalizationCheckMC::IsNSD( AliMCEvent *lMCevent ) const
{
    //Adaptation of snippet received from Alexander +
    // this: http://svnweb.cern.ch/world/wsvn/AliRoot/trunk/PWGLF/SPECTRA/ChargedHadrons/dNdPt/AlidNdPtHelper.cxx?rev=61295

    //Acquire Header
    AliHeader * header = lMCevent->Header();
    AliGenDPMjetEventHeader* dpmHeader = dynamic_cast<AliGenDPMjetEventHeader*>(header->GenEventHeader());
    if (!dpmHeader) return kFALSE;

    Int_t nsd1 = 0;
    Int_t nsd2 = 0;
    Int_t ndd  = 0;
    if (dpmHeader) dpmHeader->GetNDiffractive(nsd1, nsd2, ndd);
    
    if((dpmHeader->ProjectileParticipants()==nsd1) && (ndd==0)) { return kFALSE; }
    else if ((dpmHeader->ProjectileParticipants()==nsd2) && (ndd==0)) { return kFALSE; }
    else { return kTRUE; }
}