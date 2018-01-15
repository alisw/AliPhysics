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
// This task is meant to acquire MC-level predictions in general
// for several LF-related particle species. First deployed to deal with
// the Pb-Pb 5 TeV strangeness analysis.
//
// Please report any bugs, complaints, suggestions to:
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
#include "AliAnalysisTaskMCPredictions.h"
////
#include "AliMultEstimator.h"
#include "AliMultVariable.h"
#include "AliMultInput.h"
#include "AliMultSelection.h"


using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskMCPredictions)

AliAnalysisTaskMCPredictions::AliAnalysisTaskMCPredictions()
: AliAnalysisTaskSE(),
fListHist(0),
fHistEventCounter(0),
fHistV0MMult(0),
fHistNchVsV0MMult(0)
{

}

AliAnalysisTaskMCPredictions::AliAnalysisTaskMCPredictions(const char *name)
: AliAnalysisTaskSE(name),
fListHist(0),
fHistEventCounter(0),
fHistV0MMult(0),
fHistNchVsV0MMult(0)
{
    for(Int_t ih=0; ih<9; ih++){
        fHistPt[ih]          = 0x0;
        fHistPtVsV0MMult[ih] = 0x0;
    }
    DefineOutput(1, TList::Class()); // Event Counter Histo
}


AliAnalysisTaskMCPredictions::~AliAnalysisTaskMCPredictions()
{
    //------------------------------------------------
    // DESTRUCTOR
    //------------------------------------------------
    
    if (fListHist) {
        delete fListHist;
        fListHist = 0x0;
    }
}

//________________________________________________________________________
void AliAnalysisTaskMCPredictions::UserCreateOutputObjects()
{
    //------------------------------------------------
    // Histograms: Basic Analysis Output
    //------------------------------------------------
    
    // Create histograms
    fListHist = new TList();
    fListHist->SetOwner();  // See http://root.cern.ch/root/html/TCollection.html#TCollection:SetOwner
    
    if(! fHistEventCounter ) {
        //Histogram Output: Event-by-Event
        fHistEventCounter = new TH1D( "fHistEventCounter", ";Evt. Sel. Step;Count",1,0,1);
        //Keeps track of some basics
        fHistEventCounter->GetXaxis()->SetBinLabel(1, "Processed");
        fListHist->Add(fHistEventCounter);
    }
    
    //Identified Particles
    Int_t lPDGCodes[9] = {211, 321, 2212, 310, 3122, 3312, 3334, 333, 313};
    TString lPartNames[9] = {
        "Pion", "Kaon", "Proton", "K0Short", "Lambda", "Xi", "Omega", "Phi", "KStar"
    };
    
    //Settings for transverse momentum
    Int_t lNPtBins = 200;
    Double_t lMaxPt = 20.0;
    
    //Settings for charged particle counters (integers!)
    Int_t lNNchBins = 3000;
    Double_t lLowNchBound  = -0.5;
    Double_t lHighNchBound = -0.5 + ((double)(lNNchBins));
    
    //Main Output: Histograms
    
    //Event counter histogram: Multiplicity, Npart, b (if available)
    for(Long_t ih=0;ih<9;ih++){
        if(! fHistPt[ih] ) {
            fHistPt[ih] = new TH1D(Form("fHistPt_%s",lPartNames[ih].Data()),    "Generated;p_{T} (GeV/c)",lNPtBins,0,lMaxPt);
            fListHist->Add(fHistPt[ih]);
        }
        if(! fHistPtVsV0MMult[ih] ) {
            fHistPtVsV0MMult[ih] = new TH2D(Form("fHistPtVsV0MMult_%s",lPartNames[ih].Data()),    "Generated;p_{T} (GeV/c)",lNNchBins,lLowNchBound,lHighNchBound,lNPtBins,0,lMaxPt);
            fListHist->Add(fHistPtVsV0MMult[ih]);
        }
    }

    //List of Histograms: Normal
    PostData(1, fListHist);
    
}// end UserCreateOutputObjects


//________________________________________________________________________
void AliAnalysisTaskMCPredictions::UserExec(Option_t *)
{
    // Main loop
    // Called for each event
    
    AliMCEvent  *lMCevent  = 0x0;
    AliStack    *lMCstack  = 0x0;
    
    
    // Connect to the InputEvent
    // After these lines, we should have an ESD/AOD event + the number of V0s in it.
    
    // Appropriate for ESD analysis!
    
    lMCevent = MCEvent();
    if (!lMCevent) {
        Printf("ERROR: Could not retrieve MC event \n");
        cout << "Name of the file with pb :" <<  fInputHandler->GetTree()->GetCurrentFile()->GetName() << endl;
        return;
    }
    
    lMCstack = lMCevent->Stack();
    if (!lMCstack) {
        Printf("ERROR: Could not retrieve MC stack \n");
        cout << "Name of the file with pb :" <<  fInputHandler->GetTree()->GetCurrentFile()->GetName() << endl;
        return;
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
    
    //------------------------------------------------
    // Fill Event Counters
    //------------------------------------------------

    //Basics: All Processed
    fHistEventCounter->Fill(0.5);
    
    fHistV0MMult -> Fill ( lNchVZEROA+lNchVZEROC );
    fHistNchVsV0MMult -> Fill ( lNchVZEROA+lNchVZEROC, lNchEta5  );
    //fHistNpart
    //fHistNchVsNpart
    //fHistB
    //fHistNchVsB
    
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
                fHistPt[ih]->Fill(lThisPt);
                fHistPtVsV0MMult[ih]->Fill(lNchVZEROA+lNchVZEROC,lThisPt);
            }
        }
    }//End of loop on tracks
    //----- End Loop on Stack ----------------------
    
    // Post output data.
    PostData(1, fListHist);
}

//________________________________________________________________________
void AliAnalysisTaskMCPredictions::Terminate(Option_t *)
{
    // Draw result to the screen
    // Called once at the end of the query
    
    TList *cRetrievedList = 0x0;
    cRetrievedList = (TList*)GetOutputData(1);
    if(!cRetrievedList) {
        Printf("ERROR - AliAnalysisTaskMCPredictions : ouput data container list not available\n");
        return;
    }
    
    fHistEventCounter = dynamic_cast<TH1D*> (  cRetrievedList->FindObject("fHistEventCounter")  );
    if (!fHistEventCounter) {
        Printf("ERROR - AliAnalysisTaskMCPredictions : fHistEventCounter not available");
        return;
    }
    
    TCanvas *canCheck = new TCanvas("AliAnalysisTaskMCPredictions","Event Multiplicity",10,10,510,510);
    canCheck->cd(1)->SetLogy();
    
    fHistEventCounter->SetMarkerStyle(22);
    fHistEventCounter->DrawCopy("E");
}

//----------------------------------------------------------------------------
Double_t AliAnalysisTaskMCPredictions::MyRapidity(Double_t rE, Double_t rPz) const
{
    // Local calculation for rapidity
    Double_t ReturnValue = -100;
    if( (rE-rPz+1.e-13) != 0 && (rE+rPz) != 0 ) {
        ReturnValue =  0.5*TMath::Log((rE+rPz)/(rE-rPz+1.e-13));
    }
    return ReturnValue;
}
