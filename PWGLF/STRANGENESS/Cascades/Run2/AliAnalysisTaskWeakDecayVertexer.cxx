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
// This task is meant to encapsulate weak decay revertexing in a fully
// configurable way, so that performance-enhancing additions can be added
// already at the Tracks2Vertices level inside the V0 and cascade finding
// loops.
//
// Disclaimer: This is work in progress! Use at your own discretion!
//
// For questions, comments, etc, please write to:
//      david.dobrigkeit.chinellato@cern.ch
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
#include "TRandom3.h"
#include "TLorentzVector.h"
//#include "AliLog.h"

#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliV0vertexer.h"
#include "AliCascadeVertexer.h"
#include "AliLightV0vertexer.h"
#include "AliLightCascadeVertexer.h"
#include "AliESDpid.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliInputEventHandler.h"
#include "AliAnalysisManager.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliCentrality.h"
#include "AliOADBContainer.h"
#include "AliOADBMultSelection.h"
#include "AliMultEstimator.h"
#include "AliMultVariable.h"
#include "AliMultInput.h"
#include "AliMultSelection.h"

#include "AliCFContainer.h"
#include "AliMultiplicity.h"
#include "AliAODMCParticle.h"
#include "AliESDcascade.h"
#include "AliAODcascade.h"
#include "AliESDUtils.h"
#include "AliGenEventHeader.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisUtils.h"
#include "AliEventCuts.h"
#include "AliV0Result.h"
#include "AliCascadeResult.h"
#include "AliAnalysisTaskWeakDecayVertexer.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskWeakDecayVertexer)

AliAnalysisTaskWeakDecayVertexer::AliAnalysisTaskWeakDecayVertexer()
    : AliAnalysisTaskSE(), fListHist(0), fPIDResponse(0),
//________________________________________________
//Options for general task operation
fTrigType(AliVEvent::kMB),
fkDoExtraEvSels(kTRUE),
fMinCentrality(0.0),
fMaxCentrality(90.0),
//________________________________________________
//Flags for both V0+cascade vertexer
fkPreselectDedx ( kTRUE ),
fkPreselectDedxLambda ( kTRUE ),
fkExtraCleanup    ( kTRUE ), //extra cleanup: eta, etc
//________________________________________________
//Flags for V0 vertexer
fkRunV0Vertexer (kFALSE),
fkDoV0Refit       ( kTRUE ),
//________________________________________________
//Flags for cascade vertexer
fkRunCascadeVertexer    ( kFALSE ),
fkUseOnTheFlyV0Cascading( kFALSE ),
fMinPtCascade(   0.3 ),
fMaxPtCascade( 100.00 ),
fMassWindowAroundCascade(0.060),
//________________________________________________
//Histos
fHistEventCounter(0),
fHistCentrality(0),
fHistNumberOfCandidates(0) //bookkeep total number of candidates analysed
//________________________________________________
{

}

AliAnalysisTaskWeakDecayVertexer::AliAnalysisTaskWeakDecayVertexer(const char *name, TString lExtraOptions)
: AliAnalysisTaskSE(name), fListHist(0), fPIDResponse(0),
//________________________________________________
//Options for general task operation
fTrigType(AliVEvent::kMB),
fkDoExtraEvSels(kTRUE),
fMinCentrality(0.0),
fMaxCentrality(90.0),
//________________________________________________
//Flags for both V0+cascade vertexer
fkPreselectDedx ( kTRUE ),
fkPreselectDedxLambda ( kTRUE ),
fkExtraCleanup    ( kTRUE ), //extra cleanup: eta, etc
//________________________________________________
//Flags for V0 vertexer
fkRunV0Vertexer (kFALSE),
fkDoV0Refit       ( kTRUE ),
//________________________________________________
//Flags for cascade vertexer
fkRunCascadeVertexer    ( kFALSE ),
fkUseOnTheFlyV0Cascading( kFALSE ),
fMinPtCascade(   0.3 ), //pre-selection
fMaxPtCascade( 100.00 ),
fMassWindowAroundCascade(0.060),
//________________________________________________
//Histos
fHistEventCounter(0),
fHistCentrality(0),
fHistNumberOfCandidates(0) //bookkeep total number of candidates analysed
//________________________________________________
{

    //Re-vertex: Will only apply for cascade candidates

    fV0VertexerSels[0] =  33.  ;  // max allowed chi2
    fV0VertexerSels[1] =   0.02;  // min allowed impact parameter for the 1st daughter (LHC09a4 : 0.05)
    fV0VertexerSels[2] =   0.02;  // min allowed impact parameter for the 2nd daughter (LHC09a4 : 0.05)
    fV0VertexerSels[3] =   2.0 ;  // max allowed DCA between the daughter tracks       (LHC09a4 : 0.5)
    fV0VertexerSels[4] =   0.95;  // min allowed cosine of V0's pointing angle         (LHC09a4 : 0.99)
    fV0VertexerSels[5] =   1.0 ;  // min radius of the fiducial volume                 (LHC09a4 : 0.2)
    fV0VertexerSels[6] = 200.  ;  // max radius of the fiducial volume                 (LHC09a4 : 100.0)

    fCascadeVertexerSels[0] =  33.   ;  // max allowed chi2 (same as PDC07)
    fCascadeVertexerSels[1] =   0.05 ;  // min allowed V0 impact parameter                    (PDC07 : 0.05   / LHC09a4 : 0.025 )
    fCascadeVertexerSels[2] =   0.010;  // "window" around the Lambda mass                    (PDC07 : 0.008  / LHC09a4 : 0.010 )
    fCascadeVertexerSels[3] =   0.03 ;  // min allowed bachelor's impact parameter            (PDC07 : 0.035  / LHC09a4 : 0.025 )
    fCascadeVertexerSels[4] =   2.0  ;  // max allowed DCA between the V0 and the bachelor    (PDC07 : 0.1    / LHC09a4 : 0.2   )
    fCascadeVertexerSels[5] =   0.95 ;  // min allowed cosine of the cascade pointing angle   (PDC07 : 0.9985 / LHC09a4 : 0.998 )
    fCascadeVertexerSels[6] =   0.4  ;  // min radius of the fiducial volume                  (PDC07 : 0.9    / LHC09a4 : 0.2   )
    fCascadeVertexerSels[7] = 200.   ;  // max radius of the fiducial volume                  (PDC07 : 100    / LHC09a4 : 100   )

    DefineOutput(1, TList::Class()); // Basic Histograms
}


AliAnalysisTaskWeakDecayVertexer::~AliAnalysisTaskWeakDecayVertexer()
{
    //------------------------------------------------
    // DESTRUCTOR
    //------------------------------------------------

    //Destroy output objects if present
    if (fListHist) {
        delete fListHist;
        fListHist = 0x0;
    }
}

//________________________________________________________________________
void AliAnalysisTaskWeakDecayVertexer::UserCreateOutputObjects()
{
    //------------------------------------------------
    // Particle Identification Setup
    //------------------------------------------------

    AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
    fPIDResponse = inputHandler->GetPIDResponse();
    inputHandler->SetNeedField();

    //------------------------------------------------
    // V0 Multiplicity Histograms
    //------------------------------------------------

    // Create histograms
    fListHist = new TList();
    fListHist->SetOwner();  // See http://root.cern.ch/root/html/TCollection.html#TCollection:SetOwner

    fEventCuts.AddQAplotsToList(fListHist);

    if(! fHistEventCounter ) {
        //Histogram Output: Event-by-Event
        fHistEventCounter = new TH1D( "fHistEventCounter", ";Evt. Sel. Step;Count",2,0,2);
        fHistEventCounter->GetXaxis()->SetBinLabel(1, "Processed");
        fHistEventCounter->GetXaxis()->SetBinLabel(2, "Selected");
        fListHist->Add(fHistEventCounter);
    }

    if(! fHistCentrality ) {
        //Histogram Output: Event-by-Event
        fHistCentrality = new TH1D( "fHistCentrality", "WARNING: no pileup rejection applied!;Centrality;Event Count",100,0,100);
        fListHist->Add(fHistCentrality);
    }
    
    if(! fHistNumberOfCandidates ) {
        //Histogram Output: Event-by-Event
        fHistNumberOfCandidates = new TH1D( "fHistNumberOfCandidates", "Candidate count;Centrality;Event Count",4,0,4);
        fHistNumberOfCandidates->GetXaxis()->SetBinLabel(1, "V0s: original");
        fHistNumberOfCandidates->GetXaxis()->SetBinLabel(2, "V0s: re-vertexed");
        fHistNumberOfCandidates->GetXaxis()->SetBinLabel(3, "Cascades: original");
        fHistNumberOfCandidates->GetXaxis()->SetBinLabel(4, "Cascades: re-vertexed");
        fListHist->Add(fHistNumberOfCandidates);
    }

    PostData(1, fListHist    );
}// end UserCreateOutputObjects


//________________________________________________________________________
void AliAnalysisTaskWeakDecayVertexer::UserExec(Option_t *)
{
    // Main loop
    // Called for each event

    AliESDEvent *lESDevent = 0x0;

    // Connect to the InputEvent
    // After these lines, we should have an ESD/AOD event + the number of V0s in it.

    // Appropriate for ESD analysis!

    lESDevent = dynamic_cast<AliESDEvent*>( InputEvent() );
    if (!lESDevent) {
        AliWarning("ERROR: lESDevent not available \n");
        return;
    }

    Double_t lMagneticField = -10;
    lMagneticField = lESDevent->GetMagneticField( );

    fHistEventCounter->Fill(0.5);

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

    //------------------------------------------------
    // Multiplicity Information Acquistion, basic event selection
    //------------------------------------------------

    Float_t lPercentile = 500;
    Int_t lEvSelCode = 100;
    AliMultSelection *MultSelection = (AliMultSelection*) lESDevent -> FindListObject("MultSelection");
    if( !MultSelection) {
        //If you get this warning (and lPercentiles 300) please check that the AliMultSelectionTask actually ran (before your task)
        AliWarning("AliMultSelection object not found!");
    } else {
        //V0M Multiplicity Percentile
        lPercentile = MultSelection->GetMultiplicityPercentile("V0M");
        //Event Selection Code
        lEvSelCode = MultSelection->GetEvSelCode();
    }

    if( lEvSelCode != 0 ) {
        PostData(1, fListHist    );
        return;
    }

    AliVEvent *ev = InputEvent();
    if( fkDoExtraEvSels ) {
        if( !fEventCuts.AcceptEvent(ev) ) {
            PostData(1, fListHist    );
            return;
        }
    }

    fHistEventCounter->Fill(1.5);

    //Fill centrality histogram
    fHistCentrality->Fill(lPercentile);

    
    //============================================================
    // V0 Revertexing part
    //============================================================
    
    //bookkeep original number of v0s
    Int_t nv0s = 0;
    nv0s = lESDevent->GetNumberOfV0s();
    fHistNumberOfCandidates->Fill(0.5, nv0s);
    
    if( fkRunV0Vertexer ){
        lESDevent->ResetV0s();
        //Only regenerate candidates if within interesting interval
        if( lPercentile>fMinCentrality && lPercentile<fMaxCentrality ){
            Tracks2V0vertices(lESDevent);
        }
    }
    
    nv0s = lESDevent->GetNumberOfV0s();
    fHistNumberOfCandidates->Fill(1.5, nv0s);

    //============================================================
    // Cascade revertexing part
    //============================================================
    
    //bookkeep original number of cascades
    Long_t ncascades = 0;
    ncascades = lESDevent->GetNumberOfCascades();
    fHistNumberOfCandidates->Fill(2.5, ncascades);

    if( fkRunCascadeVertexer ){
        lESDevent->ResetCascades();
        //Only regenerate candidates if within interesting interval
        if( lPercentile>fMinCentrality && lPercentile<fMaxCentrality ){
            V0sTracks2CascadeVertices(lESDevent);
        }
    }
    
    ncascades = lESDevent->GetNumberOfCascades();
    fHistNumberOfCandidates->Fill(3.5, ncascades);

    // Post output data: end of times
    PostData(1, fListHist    );
}

//________________________________________________________________________
void AliAnalysisTaskWeakDecayVertexer::Terminate(Option_t *)
{
    // Draw result to the screen
    // Called once at the end of the query

    TList *cRetrievedList = 0x0;
    cRetrievedList = (TList*)GetOutputData(1);
    if(!cRetrievedList) {
        Printf("ERROR - AliAnalysisTaskWeakDecayVertexer : ouput data container list not available\n");
        return;
    }

    fHistEventCounter = dynamic_cast<TH1D*> (  cRetrievedList->FindObject("fHistEventCounter")  );
    if (!fHistEventCounter) {
        Printf("ERROR - AliAnalysisTaskWeakDecayVertexer : fHistEventCounter not available");
        return;
    }

    TCanvas *canCheck = new TCanvas("AliAnalysisTaskWeakDecayVertexer","V0 Multiplicity",10,10,510,510);
    canCheck->cd(1)->SetLogy();

    fHistEventCounter->SetMarkerStyle(22);
    fHistEventCounter->DrawCopy("E");
}

//________________________________________________________________________
void AliAnalysisTaskWeakDecayVertexer::SetupStandardVertexing()
//Meant to store standard re-vertexing configuration
{
    //Tell the task to re-run vertexers
    SetRunV0Vertexer(kTRUE);
    SetRunCascadeVertexer(kTRUE);
    SetDoV0Refit(kTRUE);

    //V0-Related topological selections
    SetV0VertexerDCAFirstToPV(0.05);
    SetV0VertexerDCASecondtoPV(0.05);
    SetV0VertexerDCAV0Daughters(1.20);
    SetV0VertexerCosinePA(0.98);
    SetV0VertexerMinRadius(0.9);
    SetV0VertexerMaxRadius(200);

    //Cascade-Related topological selections
    SetCascVertexerMinV0ImpactParameter(0.05);
    SetCascVertexerV0MassWindow(0.006);
    SetCascVertexerDCABachToPV(0.02);
    SetCascVertexerDCACascadeDaughters(1.2);
    SetCascVertexerCascadeMinRadius(.8);
    SetCascVertexerCascadeCosinePA(.98);
}

//________________________________________________________________________
void AliAnalysisTaskWeakDecayVertexer::SetupLooseVertexing()
//Meant to store standard re-vertexing configuration
{
    //Tell the task to re-run vertexers
    SetRunV0Vertexer(kTRUE);
    SetRunCascadeVertexer(kTRUE);
    SetDoV0Refit(kTRUE);
    
    //V0-Related topological selections
    SetV0VertexerDCAFirstToPV(0.1);
    SetV0VertexerDCASecondtoPV(0.1);
    SetV0VertexerDCAV0Daughters(1.40);
    SetV0VertexerCosinePA(0.95);
    SetV0VertexerMinRadius(0.9);
    SetV0VertexerMaxRadius(200);
    
    //Cascade-Related topological selections
    SetCascVertexerMinV0ImpactParameter(0.05);
    SetCascVertexerV0MassWindow(0.006);
    SetCascVertexerDCABachToPV(0.02);
    SetCascVertexerDCACascadeDaughters(1.4);
    SetCascVertexerCascadeMinRadius(.5);
    SetCascVertexerCascadeCosinePA(.95);
}

//________________________________________________________________________
Long_t AliAnalysisTaskWeakDecayVertexer::Tracks2V0vertices(AliESDEvent *event) {
    //--------------------------------------------------------------------
    //This function reconstructs V0 vertices
    //--------------------------------------------------------------------
    
    const AliESDVertex *vtxT3D=event->GetPrimaryVertex();
    
    Double_t xPrimaryVertex=vtxT3D->GetX();
    Double_t yPrimaryVertex=vtxT3D->GetY();
    Double_t zPrimaryVertex=vtxT3D->GetZ();
    
    Long_t nentr=event->GetNumberOfTracks();
    Double_t b=event->GetMagneticField();
    
    if (nentr<2) return 0;
    
    TArrayI neg(nentr);
    TArrayI pos(nentr);
    
    Long_t nneg=0, npos=0, nvtx=0;
    
    Long_t i;
    for (i=0; i<nentr; i++) {
        AliESDtrack *esdTrack=event->GetTrack(i);
        ULong_t status=esdTrack->GetStatus();
        
        //if ((status&AliESDtrack::kITSrefit)==0)//not to accept the ITS SA tracks
        if ((status&AliESDtrack::kTPCrefit)==0) continue;
        
        //Track pre-selection: clusters
        if (esdTrack->GetTPCNcls() < 70 && fkExtraCleanup ) continue;

        Double_t d=esdTrack->GetD(xPrimaryVertex,yPrimaryVertex,b);
        if (TMath::Abs(d)<fV0VertexerSels[2]) continue;
        if (TMath::Abs(d)>fV0VertexerSels[6]) continue;
        
        if (esdTrack->GetSign() < 0.) neg[nneg++]=i;
        else pos[npos++]=i;
    }
    
    
    for (i=0; i<nneg; i++) {
        Long_t nidx=neg[i];
        AliESDtrack *ntrk=event->GetTrack(nidx);
        
        for (Int_t k=0; k<npos; k++) {
            Int_t pidx=pos[k];
            AliESDtrack *ptrk=event->GetTrack(pidx);
            
            //Pre-select dE/dx: only proceed if at least one of these tracks looks like a proton
            if(fkPreselectDedxLambda){
                Double_t lNSigPproton = TMath::Abs(fPIDResponse->NumberOfSigmasTPC( ptrk, AliPID::kProton ));
                Double_t lNSigNproton = TMath::Abs(fPIDResponse->NumberOfSigmasTPC( ntrk, AliPID::kProton ));
                if( lNSigPproton>5.0 && lNSigNproton>5.0 ) continue; 
            }
            
            //Track pre-selection: clusters
            if (ptrk->GetTPCNcls() < 70 &&fkExtraCleanup ) continue;
            
            if (TMath::Abs(ntrk->GetD(xPrimaryVertex,yPrimaryVertex,b))<fV0VertexerSels[1])
                if (TMath::Abs(ptrk->GetD(xPrimaryVertex,yPrimaryVertex,b))<fV0VertexerSels[2]) continue;
            
            Double_t xn, xp, dca=ntrk->GetDCA(ptrk,b,xn,xp);
            if (dca > fV0VertexerSels[3]) continue;
            if ((xn+xp) > 2*fV0VertexerSels[6]) continue;
            if ((xn+xp) < 2*fV0VertexerSels[5]) continue;
            
            AliExternalTrackParam nt(*ntrk), pt(*ptrk);
            Bool_t corrected=kFALSE;
            if ((nt.GetX() > 3.) && (xn < 3.)) {
                //correct for the beam pipe material
                corrected=kTRUE;
            }
            if ((pt.GetX() > 3.) && (xp < 3.)) {
                //correct for the beam pipe material
                corrected=kTRUE;
            }
            if (corrected) {
                dca=nt.GetDCA(&pt,b,xn,xp);
                if (dca > fV0VertexerSels[3]) continue;
                if ((xn+xp) > 2*fV0VertexerSels[6]) continue;
                if ((xn+xp) < 2*fV0VertexerSels[5]) continue;
            }
            
            nt.PropagateTo(xn,b); pt.PropagateTo(xp,b);
            
            //select maximum eta range (after propagation)
            if (TMath::Abs(nt.Eta())>0.8&&fkExtraCleanup) continue;
            if (TMath::Abs(pt.Eta())>0.8&&fkExtraCleanup) continue;
            
            AliESDv0 vertex(nt,nidx,pt,pidx);
            
            //Experimental: refit V0 if asked to do so
            if( fkDoV0Refit ) vertex.Refit();
            
            //No selection: it was not previously applied, don't  apply now.
            //if (vertex.GetChi2V0() > fChi2max) continue;
            
            Double_t x=vertex.Xv(), y=vertex.Yv();
            Double_t r2=x*x + y*y;
            if (r2 < fV0VertexerSels[5]*fV0VertexerSels[5]) continue;
            if (r2 > fV0VertexerSels[6]*fV0VertexerSels[6]) continue;
            
            Float_t cpa=vertex.GetV0CosineOfPointingAngle(xPrimaryVertex,yPrimaryVertex,zPrimaryVertex);
            
            //Simple cosine cut (no pt dependence for now)
            if (cpa < fV0VertexerSels[4]) continue;
            
            vertex.SetDcaV0Daughters(dca);
            vertex.SetV0CosineOfPointingAngle(cpa);
            vertex.ChangeMassHypothesis(kK0Short);
            
            event->AddV0(&vertex);
            
            nvtx++;
        }
    }
    Info("Tracks2V0vertices","Number of reconstructed V0 vertices: %ld",nvtx);
    return nvtx;
}

//________________________________________________________________________
Long_t AliAnalysisTaskWeakDecayVertexer::V0sTracks2CascadeVertices(AliESDEvent *event) {
    //--------------------------------------------------------------------
    // This function reconstructs cascade vertices
    //      Adapted to the ESD by I.Belikov (Jouri.Belikov@cern.ch)
    //--------------------------------------------------------------------
    const AliESDVertex *vtxT3D=event->GetPrimaryVertex();
    
    Double_t xPrimaryVertex=vtxT3D->GetX();
    Double_t yPrimaryVertex=vtxT3D->GetY();
    Double_t zPrimaryVertex=vtxT3D->GetZ();
    
    Double_t b=event->GetMagneticField();
    Int_t nV0=(Int_t)event->GetNumberOfV0s();
    
    //stores relevant V0s in an array
    TObjArray vtcs(nV0);
    Long_t i;
    for (i=0; i<nV0; i++) {
        AliESDv0 *v=event->GetV0(i);
        if ( v->GetOnFlyStatus() && !fkUseOnTheFlyV0Cascading) continue;
        if (!v->GetOnFlyStatus() &&  fkUseOnTheFlyV0Cascading) continue;
        
        //Fix incorrect storing of charges in on-the-fly V0s
        if( fkUseOnTheFlyV0Cascading ){
            //Fix charge ordering
            CheckChargeV0( v );
            //Remove like-sign
            if( v->GetParamN()->Charge() > 0 && v->GetParamP()->Charge() > 0 ){
                continue;
            }
            if( v->GetParamN()->Charge() < 0 && v->GetParamP()->Charge() < 0 ){
                continue;
            }
            //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
            //Pre-filter on-the-fly candidates for CPU time reduction
            //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
            
            //1) DCA V0 Dau
            if( v->GetDcaV0Daughters() > fV0VertexerSels[3] ) continue;
            
            //2) CosPA
            if( v->GetV0CosineOfPointingAngle(xPrimaryVertex,yPrimaryVertex,zPrimaryVertex)< fV0VertexerSels[4] ) continue;
            
            //3) DCA neg/pos to PV
            UInt_t lKeyPos = (UInt_t)TMath::Abs(v->GetPindex());
            UInt_t lKeyNeg = (UInt_t)TMath::Abs(v->GetNindex());
            AliESDtrack *pTrack=((AliESDEvent*)event)->GetTrack(lKeyPos);
            AliESDtrack *nTrack=((AliESDEvent*)event)->GetTrack(lKeyNeg);
            
            Double_t lDcaPosToPrimVertex = TMath::Abs(pTrack->GetD(xPrimaryVertex,
                                                                   yPrimaryVertex,
                                                                   b) );
            
            Double_t lDcaNegToPrimVertex = TMath::Abs(nTrack->GetD(xPrimaryVertex,
                                                                   yPrimaryVertex,
                                                                   b) );
            if(lDcaNegToPrimVertex<fV0VertexerSels[1]||lDcaPosToPrimVertex<fV0VertexerSels[2]) continue; //ignore if too close
            
            //4) V0 Decay Radius
            Double_t tDecayVertexV0[3];
            v->GetXYZ(tDecayVertexV0[0],tDecayVertexV0[1],tDecayVertexV0[2]);
            Double_t lV0Radius = TMath::Sqrt(tDecayVertexV0[0]*tDecayVertexV0[0]+tDecayVertexV0[1]*tDecayVertexV0[1]);
            if(lV0Radius<fV0VertexerSels[5]) continue;
            
            //5) kTPC refit check
            // TPC refit condition (done during reconstruction for Offline but not for On-the-fly)
            if( !(pTrack->GetStatus() & AliESDtrack::kTPCrefit)) continue;
            if( !(nTrack->GetStatus() & AliESDtrack::kTPCrefit)) continue;
            
            //6) 70 clusters
            if ( ( ( pTrack->GetTPCClusterInfo(2,1) ) < 70 ) || ( ( nTrack->GetTPCClusterInfo(2,1) ) < 70 ) ) continue;
            
            //7) Daughter eta
            Double_t lNegEta = nTrack->Eta();
            Double_t lPosEta = pTrack->Eta();
            if( TMath::Abs(lNegEta)>0.8 || TMath::Abs(lPosEta)>0.8 ) continue;
            
            //8) dE/dx
            //Pre-select dE/dx: only proceed if at least one of these tracks looks like a proton
            if(fkPreselectDedxLambda){
                Double_t lNSigPproton = TMath::Abs(fPIDResponse->NumberOfSigmasTPC( pTrack, AliPID::kProton ));
                Double_t lNSigNproton = TMath::Abs(fPIDResponse->NumberOfSigmasTPC( nTrack, AliPID::kProton ));
                if( lNSigPproton>5.0 && lNSigNproton>5.0 ) continue;
            }
            //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
        }
        
        if (v->GetD(xPrimaryVertex,yPrimaryVertex,zPrimaryVertex)<fCascadeVertexerSels[1]) continue;
        vtcs.AddLast(v);
    }
    nV0=vtcs.GetEntriesFast();
    
    // stores relevant tracks in another array
    Long_t nentr=(Int_t)event->GetNumberOfTracks();
    TArrayI trk(nentr); Long_t ntr=0;
    for (i=0; i<nentr; i++) {
        AliESDtrack *esdtr=event->GetTrack(i);
        ULong_t status=esdtr->GetStatus();
        
        if ((status&AliESDtrack::kITSrefit)==0)
            if ((status&AliESDtrack::kTPCrefit)==0) continue;
        
        //Track pre-selection: clusters
        if (esdtr->GetTPCNcls() < 70 ) continue;
        
        if (TMath::Abs(esdtr->GetD(xPrimaryVertex,yPrimaryVertex,b))<fCascadeVertexerSels[3]) continue;
        trk[ntr++]=i;
    }
    
    Double_t massLambda=1.11568;
    Long_t ncasc=0;
    
    // Looking for the cascades...
    for (i=0; i<nV0; i++) { //loop on V0s
        AliESDv0 *v=(AliESDv0*)vtcs.UncheckedAt(i);
        AliESDv0 v0(*v);
        v0.ChangeMassHypothesis(kLambda0); // the v0 must be Lambda
        if (TMath::Abs(v0.GetEffMass()-massLambda)>fCascadeVertexerSels[2]) continue;
        for (Int_t j=0; j<ntr; j++) {//loop on tracks
            Int_t bidx=trk[j];
            //Bo:   if (bidx==v->GetNindex()) continue; //bachelor and v0's negative tracks must be different
            if (bidx==v0.GetIndex(0)) continue; //Bo:  consistency 0 for neg
            
            AliESDtrack *btrk=event->GetTrack(bidx);
            
            if (btrk->GetSign()>0) continue;  // bachelor's charge
            
            AliESDv0 *pv0=&v0;
            AliExternalTrackParam bt(*btrk), *pbt=&bt;
            
            Double_t dca=PropagateToDCA(pv0,pbt,b);
            if (dca > fCascadeVertexerSels[4]) continue;
            
            //eta cut - test
            if (TMath::Abs(pbt->Eta())>0.8) continue;
            
            AliESDcascade cascade(*pv0,*pbt,bidx);//constucts a cascade candidate
            //PH        if (cascade.GetChi2Xi() > fChi2max) continue;
            
            Double_t x,y,z; cascade.GetXYZcascade(x,y,z); // Bo: bug correction
            Double_t r2=x*x + y*y;
            if (r2 > fCascadeVertexerSels[7]*fCascadeVertexerSels[7]) continue;   // condition on fiducial zone
            if (r2 < fCascadeVertexerSels[6]*fCascadeVertexerSels[6]) continue;
            
            Double_t pxV0,pyV0,pzV0;
            pv0->GetPxPyPz(pxV0,pyV0,pzV0);
            if (x*pxV0+y*pyV0+z*pzV0 < 0) continue; //causality
            
            Double_t x1,y1,z1; pv0->GetXYZ(x1,y1,z1);
            if (r2 > (x1*x1+y1*y1)) continue;
            
            if (cascade.GetCascadeCosineOfPointingAngle(xPrimaryVertex,yPrimaryVertex,zPrimaryVertex) <fCascadeVertexerSels[5]) continue; //condition on the cascade pointing angle
            
            //Filter masses: anti-cascade hypotheses
            Double_t lV0quality = 0.;
            cascade.ChangeMassHypothesis(lV0quality , 3312); // pdg code -3312 = Xi+
            Double_t lInvMassXi = cascade.GetEffMassXi();
            cascade.ChangeMassHypothesis(lV0quality , 3334); // pdg code -3312 = Xi+
            Double_t lInvMassOmega = cascade.GetEffMassXi();
            
            //Remove if outside window of interest
            if(TMath::Abs(lInvMassXi   -1.322)>fMassWindowAroundCascade &&
               TMath::Abs(lInvMassOmega-1.672)>fMassWindowAroundCascade ) continue;
            
            cascade.SetDcaXiDaughters(dca);
            event->AddCascade(&cascade);
            ncasc++;
        } // end loop tracks
    } // end loop V0s
    
    // Looking for the anti-cascades...
    for (i=0; i<nV0; i++) { //loop on V0s
        AliESDv0 *v=(AliESDv0*)vtcs.UncheckedAt(i);
        AliESDv0 v0(*v);
        v0.ChangeMassHypothesis(kLambda0Bar); //the v0 must be anti-Lambda
        if (TMath::Abs(v0.GetEffMass()-massLambda)>fCascadeVertexerSels[2]) continue;
        
        for (Int_t j=0; j<ntr; j++) {//loop on tracks
            Int_t bidx=trk[j];
            if (bidx==v0.GetIndex(1)) continue; //Bo:  consistency 1 for pos
            
            AliESDtrack *btrk=event->GetTrack(bidx);
            
            if (btrk->GetSign()<0) continue;  // bachelor's charge
            
            AliESDv0 *pv0=&v0;
            AliExternalTrackParam bt(*btrk), *pbt=&bt;
            
            Double_t dca=PropagateToDCA(pv0,pbt,b);
            if (dca > fCascadeVertexerSels[4]) continue;
            
            //eta cut - test
            if (TMath::Abs(pbt->Eta())>0.8) continue;
            
            AliESDcascade cascade(*pv0,*pbt,bidx); //constucts a cascade candidate
            //PH         if (cascade.GetChi2Xi() > fChi2max) continue;
            
            Double_t x,y,z; cascade.GetXYZcascade(x,y,z); // Bo: bug correction
            Double_t r2=x*x + y*y;
            if (r2 > fCascadeVertexerSels[7]*fCascadeVertexerSels[7]) continue;   // condition on fiducial zone
            if (r2 < fCascadeVertexerSels[6]*fCascadeVertexerSels[6]) continue;
            
            Double_t pxV0,pyV0,pzV0;
            pv0->GetPxPyPz(pxV0,pyV0,pzV0);
            if (x*pxV0+y*pyV0+z*pzV0 < 0) continue; //causality
            
            Double_t x1,y1,z1; pv0->GetXYZ(x1,y1,z1);
            if (r2 > (x1*x1+y1*y1)) continue;
            
            if (cascade.GetCascadeCosineOfPointingAngle(xPrimaryVertex,yPrimaryVertex,zPrimaryVertex) < fCascadeVertexerSels[5]) continue; //condition on the cascade pointing angle
            
            //pre-select on pT
            Double_t lXiMomX       = 0. , lXiMomY = 0., lXiMomZ = 0.;
            Double_t lXiTransvMom  = 0. ;
            cascade.GetPxPyPz( lXiMomX, lXiMomY, lXiMomZ );
            lXiTransvMom  	= TMath::Sqrt( lXiMomX*lXiMomX   + lXiMomY*lXiMomY );
            if(lXiTransvMom<fMinPtCascade) continue;
            if(lXiTransvMom>fMaxPtCascade) continue;
            
            //Filter masses: anti-cascade hypotheses
            Double_t lV0quality = 0.;
            cascade.ChangeMassHypothesis(lV0quality , -3312); // pdg code -3312 = Xi+
            Double_t lInvMassXi = cascade.GetEffMassXi();
            cascade.ChangeMassHypothesis(lV0quality , -3334); // pdg code -3312 = Xi+
            Double_t lInvMassOmega = cascade.GetEffMassXi();
            
            //Remove if outside window of interest
            if(TMath::Abs(lInvMassXi   -1.322)>fMassWindowAroundCascade &&
               TMath::Abs(lInvMassOmega-1.672)>fMassWindowAroundCascade ) continue;
            
            cascade.SetDcaXiDaughters(dca);
            event->AddCascade(&cascade);
            ncasc++;
            
        } // end loop tracks
    } // end loop V0s
    
    Info("V0sTracks2CascadeVertices","Number of reconstructed cascades: %ld",ncasc);
    
    return ncasc;
}

//________________________________________________________________________
Double_t AliAnalysisTaskWeakDecayVertexer::Det(Double_t a00, Double_t a01, Double_t a10, Double_t a11) const {
    //--------------------------------------------------------------------
    // This function calculates locally a 2x2 determinant
    //--------------------------------------------------------------------
    return a00*a11 - a01*a10;
}

//________________________________________________________________________
Double_t AliAnalysisTaskWeakDecayVertexer::Det(Double_t a00,Double_t a01,Double_t a02,
                                      Double_t a10,Double_t a11,Double_t a12,
                                      Double_t a20,Double_t a21,Double_t a22) const {
    //--------------------------------------------------------------------
    // This function calculates locally a 3x3 determinant
    //--------------------------------------------------------------------
    return  a00*Det(a11,a12,a21,a22)-a01*Det(a10,a12,a20,a22)+a02*Det(a10,a11,a20,a21);
}

//________________________________________________________________________
Double_t AliAnalysisTaskWeakDecayVertexer::PropagateToDCA(AliESDv0 *v, AliExternalTrackParam *t, Double_t b) {
    //--------------------------------------------------------------------
    // This function returns the DCA between the V0 and the track
    //--------------------------------------------------------------------
    Double_t alpha=t->GetAlpha(), cs1=TMath::Cos(alpha), sn1=TMath::Sin(alpha);
    Double_t r[3]; t->GetXYZ(r);
    Double_t x1=r[0], y1=r[1], z1=r[2];
    Double_t p[3]; t->GetPxPyPz(p);
    Double_t px1=p[0], py1=p[1], pz1=p[2];
    
    Double_t x2,y2,z2;     // position and momentum of V0
    Double_t px2,py2,pz2;
    
    v->GetXYZ(x2,y2,z2);
    v->GetPxPyPz(px2,py2,pz2);
    
    // calculation dca
    
    Double_t dd= Det(x2-x1,y2-y1,z2-z1,px1,py1,pz1,px2,py2,pz2);
    Double_t ax= Det(py1,pz1,py2,pz2);
    Double_t ay=-Det(px1,pz1,px2,pz2);
    Double_t az= Det(px1,py1,px2,py2);
    
    Double_t dca=TMath::Abs(dd)/TMath::Sqrt(ax*ax + ay*ay + az*az);
    
    //points of the DCA
    Double_t t1 = Det(x2-x1,y2-y1,z2-z1,px2,py2,pz2,ax,ay,az)/
    Det(px1,py1,pz1,px2,py2,pz2,ax,ay,az);
    
    x1 += px1*t1; y1 += py1*t1; //z1 += pz1*t1;
    
    
    //propagate track to the points of DCA
    
    x1=x1*cs1 + y1*sn1;
    if (!t->PropagateTo(x1,b)) {
        Error("PropagateToDCA","Propagation failed !");
        return 1.e+33;
    }
    
    return dca;
}

//________________________________________________________________________
void AliAnalysisTaskWeakDecayVertexer::CheckChargeV0(AliESDv0 *v0)
{
    // This function checks charge of negative and positive daughter tracks.
    // If incorrectly defined (onfly vertexer), swaps out.
    if( v0->GetParamN()->Charge() > 0 && v0->GetParamP()->Charge() < 0 ) {
        //V0 daughter track swapping is required! Note: everything is swapped here... P->N, N->P
        Long_t lCorrectNidx = v0->GetPindex();
        Long_t lCorrectPidx = v0->GetNindex();
        Double32_t	lCorrectNmom[3];
        Double32_t	lCorrectPmom[3];
        v0->GetPPxPyPz( lCorrectNmom[0], lCorrectNmom[1], lCorrectNmom[2] );
        v0->GetNPxPyPz( lCorrectPmom[0], lCorrectPmom[1], lCorrectPmom[2] );
        
        AliExternalTrackParam	lCorrectParamN(
                                               v0->GetParamP()->GetX() ,
                                               v0->GetParamP()->GetAlpha() ,
                                               v0->GetParamP()->GetParameter() ,
                                               v0->GetParamP()->GetCovariance()
                                               );
        AliExternalTrackParam	lCorrectParamP(
                                               v0->GetParamN()->GetX() ,
                                               v0->GetParamN()->GetAlpha() ,
                                               v0->GetParamN()->GetParameter() ,
                                               v0->GetParamN()->GetCovariance()
                                               );
        lCorrectParamN.SetMostProbablePt( v0->GetParamP()->GetMostProbablePt() );
        lCorrectParamP.SetMostProbablePt( v0->GetParamN()->GetMostProbablePt() );
        
        //Get Variables___________________________________________________
        Double_t lDcaV0Daughters = v0 -> GetDcaV0Daughters();
        Double_t lCosPALocal     = v0 -> GetV0CosineOfPointingAngle();
        Bool_t lOnFlyStatusLocal = v0 -> GetOnFlyStatus();
        
        //Create Replacement Object_______________________________________
        AliESDv0 *v0correct = new AliESDv0(lCorrectParamN,lCorrectNidx,lCorrectParamP,lCorrectPidx);
        v0correct->SetDcaV0Daughters          ( lDcaV0Daughters   );
        v0correct->SetV0CosineOfPointingAngle ( lCosPALocal       );
        v0correct->ChangeMassHypothesis       ( kK0Short          );
        v0correct->SetOnFlyStatus             ( lOnFlyStatusLocal );
        
        //Reverse Cluster info..._________________________________________
        v0correct->SetClusters( v0->GetClusters( 1 ), v0->GetClusters ( 0 ) );
        
        *v0 = *v0correct;
        //Proper cleanup..._______________________________________________
        v0correct->Delete();
        v0correct = 0x0;
        
        //Just another cross-check and output_____________________________
        if( v0->GetParamN()->Charge() > 0 && v0->GetParamP()->Charge() < 0 ) {
            AliWarning("Found Swapped Charges, tried to correct but something FAILED!");
        } else {
            //AliWarning("Found Swapped Charges and fixed.");
        }
        //________________________________________________________________
    } else {
        //Don't touch it! ---
        //Printf("Ah, nice. Charges are already ordered...");
    }
    return;
}


