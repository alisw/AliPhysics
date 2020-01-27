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
//last edition:27.01.2020

#include <stdio.h>
#include <Riostream.h>

#include "TLorentzVector.h"
#include "TTree.h"
#include "TMath.h"
#include "TVector3.h"
#include "TRandom3.h"
#include "TClonesArray.h"
#include "TF1.h"
#include "TFile.h"
#include "TObject.h"

//AliPhysics/ROOT
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisUtils.h"
#include "AliCascadeVertexer.h"
#include "AliCentrality.h"
#include "AliESDEvent.h"
#include "AliESDVertex.h"
#include "AliESDcascade.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliExternalTrackParam.h"
#include "AliInputEventHandler.h"
#include "AliMultSelection.h"
#include "AliLightCascadeVertexer.h"
#include "AliV0vertexer.h"
#include "AliPIDResponse.h"
//#include "AliV0HypSel.h"

//my task
#include "AliAnalysisTaskStrangeCascadesDiscrete.h"



ClassImp(AliRunningCascadeCandidate)
ClassImp(AliRunningCascadeEvent)


ClassImp(AliAnalysisTaskStrangeCascadesDiscrete) //not sure, if needed at all

//default constructor
AliAnalysisTaskStrangeCascadesDiscrete::AliAnalysisTaskStrangeCascadesDiscrete()
: AliAnalysisTaskSE(),
fguard_CheckTrackQuality(kTRUE),
fguard_CheckCascadeQuality(kTRUE),
fguard_CheckTPCPID(kTRUE),
fkRunV0Vertexers(kTRUE),
fkRunVertexers(kTRUE),
fkUseLightVertexer(kFALSE),
fkUseOnTheFlyV0Cascading(kFALSE),
fPIDResponse(0),
fESDtrackCuts(0),
fESDtrackCutsITSsa2010(0),
fESDtrackCutsGlobal2015(0),
fUtils(0),
fRand(0),
fkDebugOOBPileup(kFALSE),
man(0x0),
inputHandler(0x0),
lESDevent(0x0),
MultSelection(0x0),
centrality(0x0),
lPrimaryBestESDVtx(0x0),
esdtrackcuts(0x0),
xi(0x0),
pTrackXi(0x0),
nTrackXi(0x0),
bachTrackXi(0x0),
hCascTraj(0x0),
lBaryonTrack(0x0),
Cascade_Track(0x0),
Cascade_Event(0x0),
fTreeCascadeAsEvent(0x0),
fMagneticField(-100.),
fPV_X(0), fPV_Y(0), fPV_Z(0), sigmamaxrunning(-100.),fkOmegaCleanMassWindow(0.1)
{
    
}

//constructor with task configurations
AliAnalysisTaskStrangeCascadesDiscrete::AliAnalysisTaskStrangeCascadesDiscrete(
                                                                               Bool_t lRunV0Vertexers,
                                                                               Bool_t lRunVertexers,
                                                                               Bool_t lUseLightVertexer,
                                                                               Bool_t lUseOnTheFlyV0Cascading,
                                                                               Bool_t lguard_CheckTrackQuality,
                                                                               Bool_t lguard_CheckCascadeQuality,
                                                                               Bool_t lguard_CheckTPCPID,
                                                                               
                                                                               Double_t lV0MaxChi2,
                                                                               Double_t lV0minDCAfirst,
                                                                               Double_t lV0minDCAsecond,
                                                                               Double_t lV0maxDCAdaughters,
                                                                               
                                                                               Double_t lV0minCosAngle,
                                                                               Double_t lV0minRadius,
                                                                               Double_t lV0maxRadius,
                                                                               
                                                                               Double_t lCascaderMaxChi2,
                                                                               Double_t lCascaderV0MinImpactParam,
                                                                               Double_t lCascaderV0MassWindow,
                                                                               Double_t lCascaderBachMinImpactParam,
                                                                               Double_t lCascaderMaxDCAV0andBach,
                                                                               Double_t lCascaderMinCosAngle,
                                                                               Double_t lCascaderMinRadius,
                                                                               Double_t lCascaderMaxRadius,
                                                                               Float_t sigmaRangeTPC,
                                                                               Float_t lOmegaCleanMassWindow,
                                                                               const char *name
                                                                               )
:AliAnalysisTaskSE(name),
fguard_CheckTrackQuality(kTRUE),
fguard_CheckCascadeQuality(kTRUE),
fguard_CheckTPCPID(kTRUE),
fkRunV0Vertexers(kTRUE),
fkRunVertexers(kTRUE),
fkUseLightVertexer(kFALSE),
fkUseOnTheFlyV0Cascading(kFALSE), //assumed to be false! do not play with it
fPIDResponse(0),
fESDtrackCuts(0),
fESDtrackCutsITSsa2010(0),
fESDtrackCutsGlobal2015(0),
fUtils(0),
fRand(0),
fkDebugOOBPileup(kFALSE),
man(0x0),
inputHandler(0x0),
lESDevent(0x0),
MultSelection(0x0),
centrality(0x0),
lPrimaryBestESDVtx(0x0),
esdtrackcuts(0x0),
xi(0x0),
pTrackXi(0x0),
nTrackXi(0x0),
bachTrackXi(0x0),
hCascTraj(0x0),
lBaryonTrack(0x0),
Cascade_Track(0x0),
Cascade_Event(0x0),
fTreeCascadeAsEvent(0x0),
fMagneticField(-100.),
fPV_X(0), fPV_Y(0), fPV_Z(0), sigmamaxrunning(-100.),fkOmegaCleanMassWindow(0.1)
{
    
    
    fkUseLightVertexer = lUseLightVertexer;
    fkRunV0Vertexers = lRunV0Vertexers;
    fkRunVertexers = lRunVertexers;
    fkUseOnTheFlyV0Cascading = lUseOnTheFlyV0Cascading;
    fguard_CheckTrackQuality = lguard_CheckTrackQuality;
    fguard_CheckCascadeQuality = lguard_CheckTrackQuality;
    fguard_CheckTPCPID = lguard_CheckTPCPID;
    
    //set the variables to rerun V0 vertexer! (Keep in mind, that these V0 are daughters and do not have to point back to the primary vertex!
    /* fV0VertexerSels[0] =  33.  ;  // max allowed chi2
     fV0VertexerSels[1] =   0.02;  // min allowed impact parameter for the 1st daughter (LHC09a4 : 0.05)
     fV0VertexerSels[2] =   0.02;  // min allowed impact parameter for the 2nd daughter (LHC09a4 : 0.05)
     fV0VertexerSels[3] =   2.0 ;  // max allowed DCA between the daughter tracks       (LHC09a4 : 0.5)
     fV0VertexerSels[4] =   0.8;  // min allowed cosine of V0's pointing angle         (LHC09a4 : 0.99)
     fV0VertexerSels[5] =   1.0 ;  // min radius of the fiducial volume                 (LHC09a4 : 0.2)
     fV0VertexerSels[6] = 200.  ;  // max radius of the fiducial volume                 (LHC09a4 : 100.0)*/
    
    fV0VertexerSels[0] = lV0MaxChi2;
    fV0VertexerSels[1] = lV0minDCAfirst;
    fV0VertexerSels[2] = lV0minDCAsecond;
    fV0VertexerSels[3] = lV0maxDCAdaughters;
    fV0VertexerSels[4] = lV0minCosAngle;
    fV0VertexerSels[5] = lV0minRadius;
    fV0VertexerSels[6] = lV0maxRadius;
    
    
    //TPC sigma range, usually just set to 3. Looser cuts? -> set to 4-5
    sigmamaxrunning = sigmaRangeTPC;
    
    //Mass window for Omega baryons, when extra clean up needed
    fkOmegaCleanMassWindow = lOmegaCleanMassWindow;
    
    //set the variables for the rerun of cascades
    fCascadeVertexerSels[0] = lCascaderMaxChi2  ;  // max allowed chi2 (same as PDC07) (33)
    fCascadeVertexerSels[1] = lCascaderV0MinImpactParam ;  // min allowed V0 impact parameter (0.05) (PDC07: 0.05/ LHC09a4: 0.025 )
    fCascadeVertexerSels[2] = lCascaderV0MassWindow;  // "window" around the Lambda mass  (0.010) (PDC07: 0.008  / LHC09a4 : 0.010 )
    fCascadeVertexerSels[3] = lCascaderBachMinImpactParam;  // min allowed bachelor's impact parameter (0.03) (PDC07: 0.035  / LHC09a4 : 0.025 )
    fCascadeVertexerSels[4] = lCascaderMaxDCAV0andBach;  // max allowed DCA between the V0 and the bachelor (2.0)  (PDC07: 0.1    / LHC09a4 : 0.2   )
    fCascadeVertexerSels[5] = lCascaderMinCosAngle;  // min allowed cosine of the cascade pointing angle (0.95)  (PDC07 : 0.9985 / LHC09a4 : 0.998 )
    fCascadeVertexerSels[6] = lCascaderMinRadius;  // min radius of the fiducial volume (0.4)   (PDC07 : 0.9    / LHC09a4 : 0.2   )
    fCascadeVertexerSels[7] = lCascaderMaxRadius;  // max radius of the fiducial volume (100.)  (PDC07 : 100    / LHC09a4 : 100   )
    
    DefineOutput(1, TTree::Class());
    
}

//destructor
AliAnalysisTaskStrangeCascadesDiscrete::~AliAnalysisTaskStrangeCascadesDiscrete()
{
    if (fPIDResponse) {
        delete fPIDResponse;
        fPIDResponse = nullptr;
    }
    if (fESDtrackCuts) {
        delete fESDtrackCuts;
        fESDtrackCuts = nullptr;
    }
    
    if (fESDtrackCutsITSsa2010) {
        delete fESDtrackCutsITSsa2010;
        fESDtrackCutsITSsa2010 = nullptr;
    }
    
    if (fESDtrackCutsGlobal2015) {
        delete fESDtrackCutsGlobal2015;
        fESDtrackCutsGlobal2015 = nullptr;
    }
    
    if (fUtils) {
        delete fUtils;
        fUtils = nullptr;
    }
    
    if (fRand) {
        delete fRand;
        fRand = nullptr;
    }
    
    if (man) {
        delete man;
        man = nullptr;
    }
    if (inputHandler) {
        delete inputHandler;
        inputHandler = nullptr;
    }
    if (lESDevent) {
        delete lESDevent;
        lESDevent = nullptr;
    }
    if (MultSelection) {
        delete MultSelection;
        MultSelection = nullptr;
    }
    
    if (centrality) {
        delete centrality;
        centrality = nullptr;
    }
    if (lPrimaryBestESDVtx) {
        delete lPrimaryBestESDVtx;
        lPrimaryBestESDVtx = nullptr;
    }
    if (esdtrackcuts) {
        delete esdtrackcuts;
        esdtrackcuts = nullptr;
    }
    if (xi) {
        delete xi;
        xi = nullptr;
    }
    
    if (pTrackXi) {
        delete pTrackXi;
        pTrackXi = nullptr;
    }
    
    if (nTrackXi) {
        delete nTrackXi;
        nTrackXi = nullptr;
    }
    
    if (bachTrackXi) {
        delete bachTrackXi;
        bachTrackXi = nullptr;
    }
    
    if (hCascTraj) {
        delete hCascTraj;
        hCascTraj = nullptr;
    }
    if (lBaryonTrack) {
        delete lBaryonTrack;
        lBaryonTrack = nullptr;
    }
    
}

//----------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------


void AliAnalysisTaskStrangeCascadesDiscrete::UserCreateOutputObjects()
{
    //------------------------------------------------------------------------------------------------------------------------------
    man=AliAnalysisManager::GetAnalysisManager();
    inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
    if (!inputHandler) {
        AliWarning("--------------------------------------- NO INPUT HANDLER FOUND---------------------------------------- \n");
        return;
    }
    else AliWarning("--------------------------------------- INPUT HANDLER FOUND indeed---------------------------------------- \n");
    
    fPIDResponse = inputHandler->GetPIDResponse();
    if (!fPIDResponse) {
        AliWarning("-------------------------------------------- NO FPIDRESPONSE FOUND!----------------------------------------- \n");
        return;
    }
    else AliWarning("-------------------------------------------- FPIDRESPONSE FOUND indeed!----------------------------------------- \n");
    inputHandler->SetNeedField();
    
    // Multiplicity
    if(! fESDtrackCuts ) {
        fESDtrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE,kFALSE);
        fESDtrackCuts->SetPtRange(0.15);  // adding pt cut
        fESDtrackCuts->SetEtaRange(-1.0, 1.0);
    }
    //Analysis Utils
    if(! fUtils ) {
        fUtils = new AliAnalysisUtils();
    }
    if(! fRand ){
        fRand = new TRandom3();
        fRand->SetSeed(0);
    }
    
    // OOB Pileup in pp 2016
    if( !fESDtrackCutsGlobal2015 && fkDebugOOBPileup ) {
        fESDtrackCutsGlobal2015 = AliESDtrackCuts::GetStandardITSTPCTrackCuts2015PbPb(kTRUE,kFALSE);
        //Initial set of cuts - to be adjusted
        fESDtrackCutsGlobal2015->SetPtRange(0.15);
        fESDtrackCutsGlobal2015->SetEtaRange(-1.0, 1.0);
    }
    if( !fESDtrackCutsITSsa2010 && fkDebugOOBPileup ) {
        fESDtrackCutsITSsa2010 = AliESDtrackCuts::GetStandardITSSATrackCuts2010();
    }
    //----------------------------------------------------------------------------------------
    
    
    Int_t iii = 0;
    Cascade_Event = new AliRunningCascadeEvent(iii);
    Cascade_Track = new AliRunningCascadeCandidate(iii);
    //  fTreeCascadeAsEvent = NULL;
    fTreeCascadeAsEvent = new TTree("fTreeCascadeAsEvent", "cascade tree as event");
    fTreeCascadeAsEvent->Branch("fTreeCascadeAsEvent_branch", "Cascade_Event" , Cascade_Event);
    
    PostData(1, fTreeCascadeAsEvent);
}



void AliAnalysisTaskStrangeCascadesDiscrete::UserExec(Option_t *option)
{
    //----------------------EVENT SELECTION--------------------------------------------------------
    //First of all, call the events and deal with events, only "good events" should be selected!
    lESDevent = dynamic_cast<AliESDEvent*>(InputEvent());
    if (!lESDevent) {
        AliWarning("ERROR: lESDevent not available \n");
        return;
    }
    
    //call all the pointers needed for the check of the event here!
    lPrimaryBestESDVtx = lESDevent -> GetPrimaryVertex();
    esdtrackcuts = new AliESDtrackCuts("esdtrackcuts", "esdtrackcuts");
    
    //check if it is a good event
    Bool_t goodevent = kTRUE;
    goodevent = GoodESDEvent(lESDevent, lPrimaryBestESDVtx, esdtrackcuts);
    if (!goodevent) return;
    //========================================================================================
    //------------------------------------------------
    // Multiplicity Information Acquistion
    //------------------------------------------------
    //Comment: fCentrality, lEvSelCode, fkUseOldCentrality should be declared in the header, so that the variables can be stored
    // we do not use any information about centralities.
    Float_t lPercentile = 500;
    Float_t fCentrality = -10;
    // Int_t lEvSelCode = 100;
    MultSelection = (AliMultSelection*) lESDevent -> FindListObject("MultSelection");
    if( !MultSelection) {
        //If you get this warning (and lPercentiles 300) please check that the AliMultSelectionTask actually ran (before your task)
        AliWarning("AliMultSelection object not found!");
    } else {
        //V0M Multiplicity Percentile
        lPercentile = MultSelection->GetMultiplicityPercentile("V0M");
        // lEvSelCode = MultSelection->GetEvSelCode();
    }
    
    //just ask AliMultSelection. It will know.
    Bool_t fMVPileupFlag = kFALSE;
    fMVPileupFlag = MultSelection->GetThisEventIsNotPileupMV();
    
    fCentrality = lPercentile;
    
    //===================================================================
    //Override centrality with equivalent run 1 info if requested, please
    Bool_t fkUseOldCentrality = kFALSE; //basically should be declared in the header, but it is not used anyway.
    if (fkUseOldCentrality) {
        centrality = lESDevent->GetCentrality();
        if ( centrality ) {
            fCentrality = centrality->GetCentralityPercentile( "V0M" );
        }
    }
    //=============================================================================================
    
    
    //Access the variables needed from the event:
    //1)Primary vertex; we deal with pp collisions, so vertex exists and |z| < 10 cm
    Double_t lBestPrimaryVtxPos[3] = {-100.,-100.,-100.};
    lPrimaryBestESDVtx->GetXYZ(lBestPrimaryVtxPos);
    fPV_X = lBestPrimaryVtxPos[0]; //originally: fTreeVariablePrimVertexX
    fPV_Y = lBestPrimaryVtxPos[1];
    fPV_Z = lBestPrimaryVtxPos[2];
    
    //2) Magnetic field
    fMagneticField = lESDevent->GetMagneticField();
    //3) number of tracks, id of the event, pileup flag, centrality(above), multiplicty(above), trigger word
    Int_t fNtracks  = lESDevent ->GetNumberOfTracks();
    Int_t eventid = lESDevent -> GetRunNumber();
    UInt_t fperiodnumber = lESDevent -> GetPeriodNumber();
    Int_t fkESDTrackMultiplicity = -1;
    fkESDTrackMultiplicity = esdtrackcuts->GetReferenceMultiplicity(lESDevent, AliESDtrackCuts::kTracklets, 0.8, 0.);
    Long64_t fkESDEventTriggerWord =lESDevent->GetTriggerMask();
    
    
    //set all the event variables into the cascade event, which will be stored.
    Int_t sizeeventbefore(0);
    sizeeventbefore = Cascade_Event -> GetSizeEvent();
    //   std::cout << "Size of the Cascade_Event array BEFORE clearing= " << sizeeventbefore << std::endl;
    
    Cascade_Event -> ClearTrackList();
    Int_t sizeeventafter(0);
    sizeeventafter = Cascade_Event -> GetSizeEvent();
    //   std::cout << "Size of the Cascade_Event array AFTER clearing = " << sizeeventafter << std::endl;
    
    Cascade_Event -> setx(fPV_X);
    Cascade_Event -> sety(fPV_Y);
    Cascade_Event -> setz(fPV_Z);
    Cascade_Event -> setN_tracks(fNtracks);
    Cascade_Event -> setid(eventid);
    Cascade_Event -> set_periodnumber(fperiodnumber);
    Cascade_Event -> setcentrality(fCentrality);
    Cascade_Event -> setMVPPileUpFlag(fMVPileupFlag);
    Cascade_Event -> setmultiplicity(fkESDTrackMultiplicity);
    Cascade_Event -> settrigger_word(fkESDEventTriggerWord);
    Cascade_Event -> setmagfield(fMagneticField);
    
    
    
    //=======================================================================================================
    //======================================EVENT DATA AQUISITION IS OVER HERE================================
    //=========================SELECTION OF CASCADES BEGINS!==================================================
    
    
    
    //HINT: light vertexers might be still in the development mode, so use them at your own risk
    
    //-----------------------------------------------------
    //-------- 1)Rerun the V0 vertexer! -------------------
    //-----------------------------------------------------
    if( fkRunV0Vertexers ) {
        //Remove existing cascades
        lESDevent->ResetV0s();
        AliV0vertexer lV0Vtxer;
        lV0Vtxer.SetDefaultCuts(fV0VertexerSels);
        lV0Vtxer.SetCuts(fV0VertexerSels);
        lV0Vtxer.Tracks2V0vertices(lESDevent);
    }
    
    //-----------------------------------------------------
    //-------- 2)Rerun the cascade vertexer! --------------
    //-----------------------------------------------------
    if( fkRunVertexers ) {
        //Remove existing cascades
        lESDevent->ResetCascades();
        //Decide between regular and light vertexer (default: light)
        if ( !fkUseLightVertexer ){
            //Instantiate vertexer object
            AliCascadeVertexer lCascVtxer;
            lCascVtxer.SetDefaultCuts(fCascadeVertexerSels);
            lCascVtxer.SetCuts(fCascadeVertexerSels);
            lCascVtxer.V0sTracks2CascadeVertices(lESDevent);
        } else {
            AliLightCascadeVertexer lCascVtxer;
            lCascVtxer.SetDefaultCuts(fCascadeVertexerSels);
            lCascVtxer.SetCuts(fCascadeVertexerSels);
            if( fkUseOnTheFlyV0Cascading ) lCascVtxer.SetUseOnTheFlyV0(kTRUE);
            lCascVtxer.V0sTracks2CascadeVertices(lESDevent);
        }
    }
    
    //-----------------------------------------------------------------
    //---------- 2) Start the loop over the cascades ------------------
    //-----------------------------------------------------------------
    
    Int_t ncascades = lESDevent -> GetNumberOfCascades();
    if(ncascades <= 0) return;
    Cascade_Event -> setNumCascadeCandidates((UShort_t)ncascades);
    
    Short_t lChargeXi = 3.;
    UShort_t iXi_passed(0);

    //--------VARIABLES USED FOR CLEANUP--------------
    Bool_t goodPos = kFALSE; Bool_t goodNeg = kFALSE; Bool_t goodBach = kFALSE;
    Double_t rapidity_range(0.5), eta_range(0.8);
    Float_t lChargePos(0.), lChargeNeg(0.);
    Bool_t extracleanupOmega = kFALSE;
    Bool_t extracleanupTPCPID = kFALSE;
   
    
    //preliminary loop is over
    for (Int_t iXi=0; iXi<ncascades; iXi++) {
        
        xi = lESDevent -> GetCascade(iXi);
        if(!xi) continue;
        lChargeXi = xi -> Charge();
        
        if (TMath::Abs(lChargeXi) != 1) continue;
        
        UInt_t lIdxPosXi     = (UInt_t) TMath::Abs( xi->GetPindex() );
        UInt_t lIdxNegXi     = (UInt_t) TMath::Abs( xi->GetNindex() );
        UInt_t lBachIdx     = (UInt_t) TMath::Abs( xi->GetBindex() );
        
        //cross check, should be fine.
        if(lBachIdx == lIdxNegXi) {
            AliWarning("Pb / Idx(Bach. track) = Idx(Neg. track) ... continue!");
            continue;
        }
        if(lBachIdx == lIdxPosXi) {
            AliWarning("Pb / Idx(Bach. track) = Idx(Pos. track) ... continue!");
            continue;
        }
        
        pTrackXi        = lESDevent->GetTrack( lIdxPosXi );
        nTrackXi        = lESDevent->GetTrack( lIdxNegXi );
        bachTrackXi    = lESDevent->GetTrack( lBachIdx );
        
        //------------CLEANUP guards----------------------------------------------------------------------------------------
        
        
        //check if tracks are good!
        goodPos = kFALSE; goodNeg = kFALSE; goodBach = kFALSE;
        rapidity_range = 0.5;
        eta_range = 0.8;
        goodPos = GoodESDTrack(pTrackXi,rapidity_range, eta_range);
        goodNeg = GoodESDTrack(nTrackXi,rapidity_range, eta_range);
        goodBach = GoodESDTrack(bachTrackXi,rapidity_range, eta_range);
        if((goodPos != kTRUE) && fguard_CheckTrackQuality) continue;
        if((goodNeg != kTRUE) && fguard_CheckTrackQuality) continue;
        if((goodBach != kTRUE) && fguard_CheckTrackQuality) continue;
        
        //cross-check for the charge consistency
        if(lChargeXi != bachTrackXi->Charge()) continue;
        //  Float_t lChargePos(0.), lChargeNeg(0.);
        lChargePos = pTrackXi->Charge(); lChargeNeg = nTrackXi->Charge();
        if((lChargePos*lChargeNeg) != -1) continue; //cross-check
        if(lChargePos < lChargeNeg) continue; //cross-check
        
        //extra cleanup for the omega!
        extracleanupOmega = kFALSE;
        //  Double_t Omegamasswindow = 0.1;
        extracleanupOmega = ExtraCleanupCascade(xi, pTrackXi, nTrackXi, bachTrackXi, fkOmegaCleanMassWindow);
        if((extracleanupOmega == kFALSE) && fguard_CheckCascadeQuality) continue;
        
        //TPC PID (might be the strongest one)
        extracleanupTPCPID = kFALSE;
        extracleanupTPCPID = GoodCandidatesTPCPID(pTrackXi, nTrackXi, bachTrackXi, sigmamaxrunning);
        if((extracleanupTPCPID == kFALSE) && fguard_CheckTPCPID) continue;
        
        
        
        //----------cleanup guards over--------------------------------------------------------------------------------------
        
        //Some of the Cascade properties. We need the hCascTraj object for DCA
        Double_t xyzCascade[3], pxpypzCascade[3], cvCascade[21];
        for(Int_t ii=0;ii<21;ii++) cvCascade[ii]=0.0; //something small
        
        xi->GetXYZcascade( xyzCascade[0],  xyzCascade[1], xyzCascade[2] );
        xi->GetPxPyPz( pxpypzCascade[0], pxpypzCascade[1], pxpypzCascade[2] );
        AliExternalTrackParam lCascTrajObject(xyzCascade,pxpypzCascade,cvCascade,lChargeXi);
        hCascTraj = &lCascTrajObject;
        
        //Get all the information for the tracks needed for the tree here:
        //First, declare all of them
        //1) Momentum of the daughters
        Double_t lBMom[3], lNMom[3], lPMom[3];
        //    Float_t clBMom[3], clNMom[3], clPMom[3];
        //2) dcas
        Float_t RunningDCAxy_z[2] = {-1.,-1.}; //running one, the first component is the transverse DCA and the second - z component
        Float_t fTreeCascVarDCAPosToPrimVtx(-1), //transverse dca of positive track to prim. vertex.
        fTreeCascVarDCAPosToPrimVtxZ(-1), //z component of the dca of positive track to prim. vertex.
        fTreeCascVarDCANegToPrimVtx(-1), //transverse dca of negative track to prim. vertex.
        fTreeCascVarDCANegToPrimVtxZ(-1), //z component of the dca of negative track to prim. vertex.
        fTreeCascVarDCABachToPrimVtx(-1), //transverse dca of bachelor track to prim. vertex.
        fTreeCascVarDCABachToPrimVtxZ(-1), //z component of the dca of bachelor track to prim. vertex.
        fTreeCascVarDCAV0ToPrimVtx(-1), //transverse DCA of the V0 to prim vertex
        fTreeCascVarCascDCAtoPVxy(-1),  //transverse DCA of the Cascade to PV
        // fTreeCascVarCascDCAtoPVz(-1),   //z component of the DCA of Cascade to PV
        fTreeCascVarDCAV0Daughters(-1),  //DCA between positive and negative tracks
        fTreeCascVarDCACascDaughters(-1), //DCA between Lambda and bach
        fTreeCascVarDCABachToBaryon(-1); //DCA of the bachelor track to the proton (baryon in general from v0)
        //Info for ITS
        Float_t fTreeCascVarPosITSNSigmaProton(-10.),
        fTreeCascVarPosITSNSigmaPion(-10.),
        fTreeCascVarNegITSNSigmaPion(-10.),
        fTreeCascVarNegITSNSigmaProton(-10.),
        fTreeCascVarBachITSNSigmaKaon(-10.),
        fTreeCascVarBachITSNSigmaPion(-10.);
        //TPC
        Float_t fTreeCascVarPosNSigmaProton(-10.),
        fTreeCascVarPosNSigmaPion(-10.),
        fTreeCascVarNegNSigmaPion(-10.),
        fTreeCascVarNegNSigmaProton(-10.),
        fTreeCascVarBachNSigmaKaon(-10.),
        fTreeCascVarBachNSigmaPion(-10.);
        //TOF
        Float_t fTreeCascVarPosTOFNSigmaProton(-10.),
        fTreeCascVarPosTOFNSigmaPion(-10.),
        fTreeCascVarNegTOFNSigmaPion(-10.),
        fTreeCascVarNegTOFNSigmaProton(-10.),
        fTreeCascVarBachTOFNSigmaKaon(-10.),
        fTreeCascVarBachTOFNSigmaPion(-10.);
        //positions of the cascade and V0;
        Float_t fTreeCascVarCascCosPointingAngle(-10.); //NOTE: here we also save the charge of the BACHELOR track!
        Float_t fTreeCascVarCascadeDecayX(1000.),
        fTreeCascVarCascadeDecayY(1000.),
        fTreeCascVarCascadeDecayZ(1000.),
        fTreeCascVarV0DecayX(1000.),
        fTreeCascVarV0DecayY(1000.),
        fTreeCascVarV0DecayZ(1000.);
        
        //  AliPIDResponse::EDetPidStatus pidstatus; //in order to access track probabilities in detectors; don't need this!
        
        
        //and now find them!
        
        
        //momenta
        xi->GetBPxPyPz( lBMom[0], lBMom[1], lBMom[2] ); //has to be Double_t
        xi->GetPPxPyPz( lPMom[0], lPMom[1], lPMom[2] );
        xi->GetNPxPyPz( lNMom[0], lNMom[1], lNMom[2] );
        //dca pos
        pTrackXi->GetDZ(fPV_X,fPV_Y,fPV_Z,fMagneticField, RunningDCAxy_z);
        fTreeCascVarDCAPosToPrimVtx = RunningDCAxy_z[0];
        fTreeCascVarDCAPosToPrimVtxZ = RunningDCAxy_z[1];
        RunningDCAxy_z[0] = -1.; RunningDCAxy_z[1] = -1.;
        //dca neg
        nTrackXi->GetDZ(fPV_X,fPV_Y,fPV_Z,fMagneticField, RunningDCAxy_z);
        fTreeCascVarDCANegToPrimVtx = RunningDCAxy_z[0];
        fTreeCascVarDCANegToPrimVtxZ = RunningDCAxy_z[1];
        RunningDCAxy_z[0] = -1.; RunningDCAxy_z[1] = -1.;
        //dca bach
        bachTrackXi->GetDZ(fPV_X,fPV_Y,fPV_Z,fMagneticField, RunningDCAxy_z);
        fTreeCascVarDCABachToPrimVtx = RunningDCAxy_z[0] * bachTrackXi->Charge(); //Also the CHARGE OF THE BACHELOR TRACK IS SAVED HERE!
        fTreeCascVarDCABachToPrimVtxZ = RunningDCAxy_z[1];
        RunningDCAxy_z[0] = -1.; RunningDCAxy_z[1] = -1.;
        //dca V0 prim
        fTreeCascVarDCAV0ToPrimVtx = xi->GetD(lBestPrimaryVtxPos[0], lBestPrimaryVtxPos[1], lBestPrimaryVtxPos[2]);
        
        //dca Casc prim
        hCascTraj->GetDZ(lBestPrimaryVtxPos[0],lBestPrimaryVtxPos[1],lBestPrimaryVtxPos[2], fMagneticField, RunningDCAxy_z);
        fTreeCascVarCascDCAtoPVxy = RunningDCAxy_z[0];
        //  fTreeCascVarCascDCAtoPVz = RunningDCAxy_z[1];
        RunningDCAxy_z[0] = -1.;
        //RunningDCAxy_z[1] = -1.;
        
        //pos and neg
        fTreeCascVarDCAV0Daughters = xi->GetDcaV0Daughters();
        //V0 and bach
        fTreeCascVarDCACascDaughters = xi->GetDcaXiDaughters();
        
        //baryon to bachelor dca
        if ( lChargeXi == -1 ){
            lBaryonTrack = pTrackXi;
        }
        if ( lChargeXi == +1 ){
            lBaryonTrack = nTrackXi;
        }
        Double_t xn, xp; //reference planes at the DCA
        fTreeCascVarDCABachToBaryon = lBaryonTrack->GetDCA(bachTrackXi, fMagneticField, xn, xp);
        
        
        //positions
        fTreeCascVarCascCosPointingAngle = xi->GetCascadeCosineOfPointingAngle(lBestPrimaryVtxPos[0], lBestPrimaryVtxPos[1], lBestPrimaryVtxPos[2]);
        Double_t PosCasc[3], PosV0[3];
        xi->GetXYZcascade( PosCasc[0], PosCasc[1], PosCasc[2]);
        xi->GetXYZ( PosV0[0], PosV0[1], PosV0[2]);
        fTreeCascVarCascadeDecayX = PosCasc[0];
        fTreeCascVarCascadeDecayY = PosCasc[1];
        fTreeCascVarCascadeDecayZ = PosCasc[2];
        fTreeCascVarV0DecayX = PosV0[0];
        fTreeCascVarV0DecayY = PosV0[1];
        fTreeCascVarV0DecayZ = PosV0[2];
        
        Float_t TrackLengthInActiveZone[3];
        // mode = 1  - Track parameters estimated at the inner wall of TPC at PropagateBack stage
        // 2.0 - deadzone, dy (user defined)
        // dz = 220
        
        TrackLengthInActiveZone[0]= pTrackXi->GetLengthInActiveZone(1, 2.0, 220.0, lESDevent->GetMagneticField());
        TrackLengthInActiveZone[1]= nTrackXi->GetLengthInActiveZone(1, 2.0, 220.0, lESDevent->GetMagneticField());
        TrackLengthInActiveZone[2]= bachTrackXi->GetLengthInActiveZone(1, 2.0, 220.0, lESDevent->GetMagneticField());
        
        
        //----------------------------ACCESS DATA: ITS------------------------------------------
        //access status on all ITS layers
        Int_t pITS[6], nITS[6], bITS[6];
        Int_t iscpos[6], iscneg[6], iscbach[6]; //number of ITS shared clusters
        Bool_t hasp(kFALSE), hasn(kFALSE), hasbach(kFALSE); //booleans for: shared clusters in its.
        for (int k=0; k<6; k++) {
            pITS[k] = GetITSstatus(pTrackXi, k);
            nITS[k] = GetITSstatus(nTrackXi, k);
            bITS[k] = GetITSstatus(bachTrackXi, k);
            
            hasp = pTrackXi -> HasSharedPointOnITSLayer(k);
            hasn = nTrackXi -> HasSharedPointOnITSLayer(k);
            hasbach = bachTrackXi -> HasSharedPointOnITSLayer(k);
            
            if(hasp == kTRUE)iscpos[k] = 1;
            else iscpos[k] =0;
            if (hasn == kTRUE)iscneg[k] =1;
            else iscpos[k] = 0;
            if(hasbach == kTRUE)iscbach[k] = 1;
            else iscbach[k] = 0;
            
        }
        //if the track was refitted in the ITS
        ULong_t status[3]; Int_t ITSrefitflag[4];
        status[0] = pTrackXi ->GetStatus(); status[1] = nTrackXi ->GetStatus(); status[2] = bachTrackXi ->GetStatus();
        //refitflag
        if((status[0] & AliESDtrack::kITSrefit) == 0) ITSrefitflag[0] = 0;
        else ITSrefitflag[0]=1;
        if((status[1] & AliESDtrack::kITSrefit) == 0) ITSrefitflag[1] = 0;
        else ITSrefitflag[1]=1;
        if((status[2] & AliESDtrack::kITSrefit) == 0) ITSrefitflag[2] = 0;
        else ITSrefitflag[2]=1;
        
        
        //its chi2
        Double_t itschi2[3];
        itschi2[0] = pTrackXi -> GetITSchi2();
        itschi2[1] = nTrackXi -> GetITSchi2();
        itschi2[2] = bachTrackXi -> GetITSchi2();
        
        //its signal
        Double_t itssignal[3];
        itssignal[0] = pTrackXi -> GetITSsignal();
        itssignal[1] = nTrackXi -> GetITSsignal();
        itssignal[2] = bachTrackXi -> GetITSsignal();
        
        //number of sigmas
        fTreeCascVarPosITSNSigmaPion   = fPIDResponse->NumberOfSigmasITS( pTrackXi, AliPID::kPion );
        fTreeCascVarPosITSNSigmaProton = fPIDResponse->NumberOfSigmasITS( pTrackXi, AliPID::kProton );
        fTreeCascVarNegITSNSigmaPion   = fPIDResponse->NumberOfSigmasITS( nTrackXi, AliPID::kPion );
        fTreeCascVarNegITSNSigmaProton = fPIDResponse->NumberOfSigmasITS( nTrackXi, AliPID::kProton );
        fTreeCascVarBachITSNSigmaPion  = fPIDResponse->NumberOfSigmasITS( bachTrackXi, AliPID::kPion );
        fTreeCascVarBachITSNSigmaKaon  = fPIDResponse->NumberOfSigmasITS( bachTrackXi, AliPID::kKaon );
        
        
        //----------------------------------ACCESS DATA: TPC---------------------------------------------------------
        //PID
        fTreeCascVarPosNSigmaPion   = fPIDResponse->NumberOfSigmasTPC( pTrackXi, AliPID::kPion );
        fTreeCascVarPosNSigmaProton = fPIDResponse->NumberOfSigmasTPC( pTrackXi, AliPID::kProton );
        fTreeCascVarNegNSigmaPion   = fPIDResponse->NumberOfSigmasTPC( nTrackXi, AliPID::kPion );
        fTreeCascVarNegNSigmaProton = fPIDResponse->NumberOfSigmasTPC( nTrackXi, AliPID::kProton );
        fTreeCascVarBachNSigmaPion  = fPIDResponse->NumberOfSigmasTPC( bachTrackXi, AliPID::kPion );
        fTreeCascVarBachNSigmaKaon  = fPIDResponse->NumberOfSigmasTPC( bachTrackXi, AliPID::kKaon );
        
        
        
        //-------------continue with TPC signal, etc.-------------
        Double_t tpcsignal[3];
        tpcsignal[0] = pTrackXi -> GetTPCsignal();
        tpcsignal[1] = nTrackXi -> GetTPCsignal();
        tpcsignal[2] = bachTrackXi -> GetTPCsignal();
        
        Double_t tpcchi2[3];
        tpcchi2[0] = pTrackXi->GetTPCchi2();
        tpcchi2[1] = nTrackXi->GetTPCchi2();
        tpcchi2[2] = bachTrackXi->GetTPCchi2();
        
        Int_t tpccls[3];
        tpccls[0] = pTrackXi -> GetTPCNcls();
        tpccls[1] = nTrackXi -> GetTPCNcls();
        tpccls[2] = bachTrackXi -> GetTPCNcls();
        
        Int_t tpcclsF[3]; //findable clusters, must be bigger or equal to number of clusters
        tpcclsF[0] = pTrackXi -> GetTPCNclsF();
        tpcclsF[1] = nTrackXi -> GetTPCNclsF();
        tpcclsF[2] = bachTrackXi -> GetTPCNclsF();
        
        //------data from TOF------------------------
        ULong_t statustof[3]; Int_t TOFrefitflag[3];
        statustof[0] = pTrackXi ->GetStatus(); statustof[1] = nTrackXi ->GetStatus(); statustof[2] = bachTrackXi ->GetStatus();
        //refitflag
        if((statustof[0] & AliESDtrack::kTOFrefit) == 0) TOFrefitflag[0] = 0;
        else TOFrefitflag[0]=1;
        if((statustof[1] & AliESDtrack::kTOFrefit) == 0) TOFrefitflag[1] = 0;
        else TOFrefitflag[1]=1;
        if((statustof[2] & AliESDtrack::kTOFrefit) == 0) TOFrefitflag[2] = 0;
        else TOFrefitflag[2]=1;
        
        Double_t tofsignal[3];
        tofsignal[0] = pTrackXi -> GetTOFsignal();
        tofsignal[1] = nTrackXi -> GetTOFsignal();
        tofsignal[2] = bachTrackXi -> GetTOFsignal();
        
        
        
        fTreeCascVarNegTOFNSigmaPion   = fPIDResponse->NumberOfSigmasTOF( nTrackXi, AliPID::kPion );
        fTreeCascVarNegTOFNSigmaProton = fPIDResponse->NumberOfSigmasTOF( nTrackXi, AliPID::kProton );
        fTreeCascVarPosTOFNSigmaPion   = fPIDResponse->NumberOfSigmasTOF( pTrackXi, AliPID::kPion );
        fTreeCascVarPosTOFNSigmaProton = fPIDResponse->NumberOfSigmasTOF( pTrackXi, AliPID::kProton );
        fTreeCascVarBachTOFNSigmaPion  = fPIDResponse->NumberOfSigmasTOF( bachTrackXi, AliPID::kPion );
        fTreeCascVarBachTOFNSigmaKaon  = fPIDResponse->NumberOfSigmasTOF( bachTrackXi, AliPID::kKaon );
        
        
        
        //-------------------set the variables to the track object and fill the tree!--------------------------------
        
        
        if(iXi_passed > 10000) continue; //see the AliRunningCascadeEvent class and the boundaries of Track array!
        
        Cascade_Track = Cascade_Event -> AddCandidate(iXi_passed);
        
        Cascade_Track -> set_PMom(lPMom[0],lPMom[1],lPMom[2]);
        Cascade_Track -> set_NMom(lNMom[0],lNMom[1],lNMom[2]);
        Cascade_Track -> set_BMom(lBMom[0],lBMom[1],lBMom[2]);
        //DCA
        Cascade_Track -> set_dca_pos_to_prim(fTreeCascVarDCAPosToPrimVtx, fTreeCascVarDCAPosToPrimVtxZ);
        Cascade_Track -> set_dca_neg_to_prim(fTreeCascVarDCANegToPrimVtx, fTreeCascVarDCANegToPrimVtxZ);
        Cascade_Track -> set_dca_bach_to_prim(fTreeCascVarDCABachToPrimVtx, fTreeCascVarDCABachToPrimVtxZ);
        Cascade_Track -> set_dca_V0_to_prim(fTreeCascVarDCAV0ToPrimVtx);
        Cascade_Track -> set_dca_Omega_to_prim(fTreeCascVarCascDCAtoPVxy);
        Cascade_Track -> set_dca_pos_to_neg(fTreeCascVarDCAV0Daughters);
        Cascade_Track -> set_dca_bach_to_Lambda(fTreeCascVarDCACascDaughters);
        Cascade_Track -> set_dca_bach_to_baryon(fTreeCascVarDCABachToBaryon);
        
        
        
        //fill with the positions of the V0 and Xi.
        Cascade_Track -> set_CosPointingAngle(fTreeCascVarCascCosPointingAngle);
        Cascade_Track -> set_CascadeDecayPos(fTreeCascVarCascadeDecayX, fTreeCascVarCascadeDecayY, fTreeCascVarCascadeDecayZ);
        Cascade_Track -> set_V0fromCascadePos(fTreeCascVarV0DecayX, fTreeCascVarV0DecayY, fTreeCascVarV0DecayZ);
        Cascade_Track -> set_TrackLengthTPC(TrackLengthInActiveZone[0], TrackLengthInActiveZone[1], TrackLengthInActiveZone[2]);
        
        //Information from ITS
        //ITS status for the daughters
        Cascade_Track -> set_ITSstatusPosTrack(pITS[0],pITS[1],pITS[2],pITS[3],pITS[4],pITS[5]);
        Cascade_Track -> set_ITSstatusNegTrack(nITS[0],nITS[1],nITS[2],nITS[3],nITS[4],nITS[5]);
        Cascade_Track -> set_ITSstatusBachTrack(bITS[0],bITS[1],bITS[2],bITS[3],bITS[4],bITS[5]);
        Cascade_Track -> set_ITSrefitFlag(ITSrefitflag[0], ITSrefitflag[1], ITSrefitflag[2]);
        Cascade_Track -> set_ITSchi2(itschi2[0], itschi2[1], itschi2[2]);
        Cascade_Track -> set_ITSsignal(itssignal[0], itssignal[1], itssignal[2]);
        
        Cascade_Track -> set_nSigma_ITS_pos(fTreeCascVarPosITSNSigmaProton, fTreeCascVarPosITSNSigmaPion);
        Cascade_Track -> set_nSigma_ITS_neg(fTreeCascVarNegITSNSigmaPion, fTreeCascVarNegITSNSigmaProton);
        Cascade_Track -> set_nSigma_ITS_bach(fTreeCascVarBachITSNSigmaKaon, fTreeCascVarBachITSNSigmaPion);
        
        Cascade_Track -> set_ITSPosSharedPoints(iscpos[0],iscpos[1],iscpos[2],iscpos[3],iscpos[4],iscpos[5]);
        Cascade_Track -> set_ITSNegSharedPoints(iscneg[0],iscneg[1],iscneg[2],iscneg[3],iscneg[4],iscneg[5]);
        Cascade_Track -> set_ITSBachSharedPoints(iscbach[0],iscbach[1],iscbach[2],iscbach[3],iscbach[4],iscbach[5]);
        
        //Information from TPC
        Cascade_Track -> set_TPCsignal(tpcsignal[0], tpcsignal[1], tpcsignal[2]);
        Cascade_Track -> set_TPCcls(tpccls[0], tpccls[1], tpccls[2]);
        Cascade_Track -> set_TPCclsF(tpcclsF[0], tpcclsF[1], tpcclsF[2]);
        Cascade_Track -> set_TPCchi2(tpcchi2[0], tpcchi2[1], tpcchi2[2]);
        Cascade_Track -> set_nSigma_dEdx_pos(fTreeCascVarPosNSigmaProton, fTreeCascVarPosNSigmaPion); //[0] Casc charge = -1, [1] Casc +
        Cascade_Track -> set_nSigma_dEdx_neg(fTreeCascVarNegNSigmaPion,  fTreeCascVarNegNSigmaProton);
        Cascade_Track -> set_nSigma_dEdx_bach(fTreeCascVarBachNSigmaKaon, fTreeCascVarBachNSigmaPion);
        
        
        //TOF
        Cascade_Track -> set_TOFsignal(tofsignal[0], tofsignal[1], tofsignal[2]);
        Cascade_Track -> set_TOFrefitFlag(TOFrefitflag[0],TOFrefitflag[1],TOFrefitflag[2]);
        Cascade_Track -> set_nSigma_TOF_pos(fTreeCascVarPosTOFNSigmaProton, fTreeCascVarPosTOFNSigmaPion);
        Cascade_Track -> set_nSigma_TOF_neg(fTreeCascVarNegTOFNSigmaPion, fTreeCascVarNegTOFNSigmaProton);
        Cascade_Track -> set_nSigma_TOF_bach(fTreeCascVarBachTOFNSigmaKaon, fTreeCascVarBachTOFNSigmaPion);
        
        iXi_passed++;

    }//end of the cascade loop.
      
    Cascade_Event -> setNumSelectedCascades(iXi_passed);
    if(iXi_passed < 1) return; //so we do not need to loop again.
    fTreeCascadeAsEvent -> Fill();
    
    PostData(1, fTreeCascadeAsEvent);
    
    
}

void AliAnalysisTaskStrangeCascadesDiscrete::Terminate(Option_t *){}

//----------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------
//----------------------------------------------extra functions---------------------------------------------
//----------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------

Bool_t AliAnalysisTaskStrangeCascadesDiscrete::GoodESDEvent(AliESDEvent* lESDevent, const AliESDVertex *lPrimaryBestESDVtx,AliESDtrackCuts* esdtrackcuts1){
    
    
    //general var to return
    Bool_t goodevent = kTRUE;
    
    //check for MAGNETIC FIELD
    Float_t lmagneticfield = -100.;
    lmagneticfield = lESDevent->GetMagneticField();
    if(!lmagneticfield) goodevent = kFALSE;
    
    //PRIMARY VERTEX can be got, should be changed for PbPb collisions
    // const AliESDVertex *lPrimaryBestESDVtx  = lESDevent->GetPrimaryVertex();
    if (!lPrimaryBestESDVtx) goodevent = kFALSE;
    Double_t lBestPrimaryVtxPos[3] = {-100.0, -100.0, -100.0};
    lPrimaryBestESDVtx->GetXYZ( lBestPrimaryVtxPos );
    //check |z| < 10
    if (TMath::Abs(lBestPrimaryVtxPos[2]) > 10.) goodevent = kFALSE;
    
    //TRACK MULTIPLICITY should be logically greater than 0.
    //gets multiplicity estimate based on TPC/ITS tracks and tracklets.
    // The negative (int) value is returned if: -1 SPD vertex missing, -2 SPD vertex z dispersion too large
    // -3 Track vertex missing (only for tracks, not checked for tracklets), -4 distance between SPD and track vertices
    // too large, not checked for tracklets.
    // AliESDtrackCuts::kTracklets -> MultEstTrackType trackType, eta range = |eta_max| = 0.8, eta central = 0.
    //MultEstTrackType trackType: kTracklets, kTrackletsITSTPC, kTrackletsITSSA (ITS SA pure tracklets).
    // AliESDtrackCuts* esdtrackcuts1 = new AliESDtrackCuts("esdtrackcuts1", "esdtrackcuts1");
    Int_t fkESDTrackMult = -1;
    fkESDTrackMult = esdtrackcuts1->GetReferenceMultiplicity(lESDevent, AliESDtrackCuts::kTracklets, 0.8, 0.);
    if (fkESDTrackMult < 0.) goodevent = kFALSE;
    
    //NUMBER OF TRACKS, check anyway
    Int_t fNtracks  = lESDevent ->GetNumberOfTracks();
    if(fNtracks < 0) goodevent = kFALSE;
    
    //Number of cascades
    Int_t ncascades = lESDevent -> GetNumberOfCascades();
    if(ncascades < 1) goodevent = kFALSE;
    
    return goodevent;
}


Bool_t AliAnalysisTaskStrangeCascadesDiscrete::GoodESDTrack(AliESDtrack* trackESD, Double_t raprange, Double_t etarange)
{
    //raprange = y_max, here set it to 0.5
    //etarange = eta_max, here set it to 0.8
    //it is assumed that the centre_rap = centre_eta = 0.
    Bool_t goodtrack = kTRUE;
    
    //is track stored?
    if(!trackESD) goodtrack = kFALSE;
    
    //check the electric charge
    Float_t charge = trackESD->Charge();
    if(TMath::Abs(charge)!=1) goodtrack = kFALSE;
    
    //check the TPC refit. (optionally we could also check for ITS refit)
    ULong_t status = trackESD->GetStatus();
    if((status & AliESDtrack::kTPCrefit) == 0) goodtrack = kFALSE;
    
    //TPC clusters
    Int_t tpcnclus=trackESD->GetTPCNcls();
    Int_t tpcnclusF = trackESD -> GetTPCNclsF();
    if(tpcnclusF<70) goodtrack = kFALSE;
    if(tpcnclus<70) goodtrack = kFALSE; //looser because of V0
    
    //filter out the noise, cross check (done automatically?) AliPerformanceTPC
    if(trackESD -> GetTPCsignal() < 5) goodtrack = kFALSE;
    
    //maximal chi2 per TPC cluster should be smaller than 4.
    Float_t tpcchi2 = 1000;
    if(tpcnclus!=0) tpcchi2= trackESD->GetTPCchi2()/tpcnclus; //was with (Float_t)trackESD->...
    if(tpcchi2 > 4)goodtrack = kFALSE;
    
    //pseudorapity
    Double_t eta = trackESD->Eta();
    if (!eta || TMath::Abs(eta) > etarange) goodtrack = kFALSE;
    
    //kink index
    if(trackESD->GetKinkIndex(0)>0 || trackESD -> GetKinkIndex(1) > 0) goodtrack = kFALSE;
    
    return goodtrack;
    
}

Int_t AliAnalysisTaskStrangeCascadesDiscrete::GetITSstatus(AliESDtrack* esdtrack, Int_t layer) const {
    //
    // Check ITS layer status
    //
    Int_t status = -2;
    Int_t fidet(0);
    Float_t fxloc(-1.), fzloc(-1.);
    if(!esdtrack) return status;
    Bool_t GetInfo = kFALSE;
    GetInfo = esdtrack -> GetITSModuleIndexInfo(layer, fidet, status, fxloc, fzloc);
    if(GetInfo == kFALSE) return -1;
    return status;
}

Bool_t AliAnalysisTaskStrangeCascadesDiscrete::ExtraCleanupCascade(AliESDcascade* Xi,
                                                                   AliESDtrack* gPtrack,
                                                                   AliESDtrack* gNtrack,
                                                                   AliESDtrack* gBtrack,
                                                                   Double_t OmegaMassWindow)
{
    //extra cleanup is needed to reduce the weight of the output file.
    //in this analysis we want to study the central prod. omega baryons, d.h. rapidity of omega |y| < 0.5
    //this cut won't be made for the Xi (using the same trio of tracks as for the omega under change of mass hypothesis!)
    //note: input tracks should have gone through the good track selection (Bool_t GoodESDTrack)
    //note: since magnetic field is just a parameter, the charges of the tracks are REAL PHYSICAL charges, so do not check for the magnetic field polarity
    Bool_t cleaned = kTRUE;
    
    if(!Xi){AliWarning("No cascade found, no clean up"); cleaned = kFALSE; }
    if(!gPtrack){AliWarning("No positive track found, no clean up"); cleaned = kFALSE; }
    if(!gNtrack){AliWarning("No negative track found, no clean up"); cleaned = kFALSE; }
    if(!gBtrack){AliWarning("No bachelor track found, no clean up"); cleaned = kFALSE; }
    
    
    TLorentzVector vXi;
    Double_t y(1.);
    Double_t mOmega(-1.);
    
    Double_t massproton = 0.9382;
    Double_t masspion = 0.1395;
    Double_t masskaon = 0.4937;
    
    Double_t Epos(0), Eneg(0),Ebach(0), EOmega(0);
    Double_t Ppos(0), Pneg(0), Pbach(0);
    //momentum of omega
    Double_t POmega[4];
    Xi -> GetPxPyPz(POmega[0], POmega[1], POmega[2]);
    POmega[3] = TMath::Sqrt(POmega[0]*POmega[0] + POmega[1]*POmega[1] + POmega[2]*POmega[2]);
    //momenta of tracks
    Ppos = gPtrack -> P();
    Pneg = gNtrack -> P();
    Pbach = gBtrack -> P();
    
    if(gBtrack -> Charge() == -1){ //omega:proton-pos, pion-neg, kaon-bach
        Epos = TMath::Sqrt(massproton*massproton + Ppos*Ppos);
        Eneg = TMath::Sqrt(masspion*masspion + Pneg*Pneg);
        Ebach = TMath::Sqrt(masskaon*masskaon + Pbach*Pbach);
        EOmega = Epos + Eneg + Ebach;
        if(EOmega <= POmega[3]) cleaned = kFALSE;
        mOmega = TMath::Sqrt(TMath::Abs(EOmega*EOmega - POmega[3]*POmega[3]));
        vXi.SetPxPyPzE(POmega[0], POmega[1], POmega[2], EOmega);
    }
    
    if(gBtrack -> Charge() == +1){ //omega:proton-pos, pion-neg, kaon-bach
        Epos = TMath::Sqrt(masspion*masspion + Ppos*Ppos);
        Eneg = TMath::Sqrt(massproton*massproton + Pneg*Pneg);
        Ebach = TMath::Sqrt(masskaon*masskaon + Pbach*Pbach);
        EOmega = Epos + Eneg + Ebach;
        if(EOmega <= POmega[3]) cleaned = kFALSE;
        mOmega = TMath::Sqrt(TMath::Abs(EOmega*EOmega - POmega[3]*POmega[3]));
        vXi.SetPxPyPzE(POmega[0], POmega[1], POmega[2], EOmega);
    }
    
    //RAPIDITY check
    y = vXi.Rapidity();
    if (TMath::Abs(y)>0.5) cleaned = kFALSE;
    
    //MASS WINDOW for OMEGA
    mOmega = TMath::Abs(vXi.M()); //M() already returns positive values only, but just in case
    if(TMath::Abs(mOmega - 1.672) > TMath::Abs(OmegaMassWindow)) cleaned = kFALSE;
    
    
    return cleaned;
}

Bool_t AliAnalysisTaskStrangeCascadesDiscrete::GoodCandidatesTPCPID(AliESDtrack* pTrackXi, AliESDtrack* nTrackXi, AliESDtrack* bachTrackXi, Float_t sigmamax)
{
    //here the TPC Pid of the Tracks is checked. There is not restriction for TOF PID, since track is not
    //necesseraly in the TOF volume (low energy tracks are prefered)
    //ITS PID cuts won't be applied either, since cascades usually fly some distance; not really seen in the ITS?
    
    //note that the daughter tracks in question fall into OMEGA mass window (pm 0.05 MeV)! Not certainly in Xi!
    
    Bool_t goodTPCtrio = kTRUE;
    if (!pTrackXi)goodTPCtrio = kFALSE;
    if (!nTrackXi)goodTPCtrio = kFALSE;
    if (!bachTrackXi)goodTPCtrio = kFALSE;
    
    Short_t chargeCascade = bachTrackXi -> Charge();
    if(TMath::Abs(chargeCascade)!=1)goodTPCtrio = kFALSE; //cross-check
    Float_t PID_Pos_Pion  = fPIDResponse->NumberOfSigmasTPC( pTrackXi, AliPID::kPion );
    Float_t PID_Pos_Proton = fPIDResponse->NumberOfSigmasTPC( pTrackXi, AliPID::kProton );
    Float_t PID_Neg_Pion   = fPIDResponse->NumberOfSigmasTPC( nTrackXi, AliPID::kPion );
    Float_t PID_Neg_Proton = fPIDResponse->NumberOfSigmasTPC( nTrackXi, AliPID::kProton );
    // Float_t PID_Bach_Pion  = fPIDResponse->NumberOfSigmasTPC( bachTrackXi, AliPID::kPion );
    Float_t PID_Bach_Kaon  = fPIDResponse->NumberOfSigmasTPC( bachTrackXi, AliPID::kKaon );
    
    if (chargeCascade == -1)
    { //positive track = proton, negative track = pion, bach track = kaon
        
        if(TMath::Abs(PID_Pos_Proton) > sigmamax) goodTPCtrio = kFALSE;
        if(TMath::Abs(PID_Neg_Pion) > sigmamax) goodTPCtrio = kFALSE;
        if(TMath::Abs(PID_Bach_Kaon) > sigmamax) goodTPCtrio = kFALSE;
        
    }
    
    if (chargeCascade == +1)
    { //positive = pion, negative = proton, bach = kaon
        
        if(TMath::Abs(PID_Pos_Pion) > sigmamax) goodTPCtrio = kFALSE;
        if(TMath::Abs(PID_Neg_Proton) > sigmamax) goodTPCtrio = kFALSE;
        if(TMath::Abs(PID_Bach_Kaon) > sigmamax) goodTPCtrio = kFALSE;
        
    }
    
    return goodTPCtrio;
}

