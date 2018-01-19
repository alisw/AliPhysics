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

#include "AliGenHijingEventHeader.h"
#include "AliGenDPMjetEventHeader.h"
#include "AliGenCocktailEventHeader.h"
#include "AliGenHepMCEventHeader.h"


using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskMCPredictions)

AliAnalysisTaskMCPredictions::AliAnalysisTaskMCPredictions()
: AliAnalysisTaskSE(),
fListHist(0),
fHistEventCounter(0),
fHistV0MMult(0),
fHistNchVsV0MMult(0),
fHistNpart(0),
fHistNchVsNpart(0),
fHistB(0),
fHistNchVsB(0)
{

}

AliAnalysisTaskMCPredictions::AliAnalysisTaskMCPredictions(const char *name)
: AliAnalysisTaskSE(name),
fListHist(0),
fHistEventCounter(0),
fHistV0MMult(0),
fHistNchVsV0MMult(0),
fHistNpart(0),
fHistNchVsNpart(0),
fHistB(0),
fHistNchVsB(0)
{
    for(Int_t ih=0; ih<9; ih++){
        fHistPt[ih]          = 0x0;
        fHistPtVsV0MMult[ih] = 0x0;
        fHistPtVsNpart[ih]   = 0x0;
        fHistPtVsB[ih]       = 0x0;
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
    
    //Settings for transverse momentum
    Int_t lNPtBins = 200;
    Double_t lMaxPt = 20.0;
    
    //Settings for charged particle counters (integers!)
    Int_t lNNchBins = 3000;
    Double_t lLowNchBound  = -0.5;
    Double_t lHighNchBound = -0.5 + ((double)(lNNchBins));
    
    if(! fHistEventCounter ) {
        //Histogram Output: Event-by-Event
        fHistEventCounter = new TH1D( "fHistEventCounter", ";Evt. Sel. Step;Count",1,0,1);
        //Keeps track of some basics
        fHistEventCounter->GetXaxis()->SetBinLabel(1, "Processed");
        fListHist->Add(fHistEventCounter);
    }
    //___________________________________________________
    if(! fHistV0MMult ) {
        //Histogram Output: Event-by-Event
        fHistV0MMult = new TH1D( "fHistV0MMult", ";V0M Mult;Count",lNNchBins,lLowNchBound,lHighNchBound);
        //Keeps track of some basics
        fListHist->Add(fHistV0MMult);
    }
    if(! fHistNchVsV0MMult ) {
        //Histogram Output: Event-by-Event
        fHistNchVsV0MMult = new TH2D( "fHistNchVsV0MMult", ";V0M Mult;Count",
                                     lNNchBins,lLowNchBound,lHighNchBound,
                                     lNNchBins,lLowNchBound,lHighNchBound);
        //Keeps track of some basics
        fListHist->Add(fHistNchVsV0MMult);
    }
    //___________________________________________________
    if(! fHistNpart ) {
        //Histogram Output: Event-by-Event
        fHistNpart = new TH1D( "fHistNpart", ";N_{part};Count",500,-0.5,499.5);
        //Keeps track of some basics
        fListHist->Add(fHistNpart);
    }
    if(! fHistNchVsNpart ) {
        //Histogram Output: Event-by-Event
        fHistNchVsNpart = new TH2D( "fHistNchVsNpart", ";N_{part};Count",500,-0.5,499.5,lNNchBins,lLowNchBound,lHighNchBound);
        //Keeps track of some basics
        fListHist->Add(fHistNchVsNpart);
    }
    //___________________________________________________
    if(! fHistB ) {
        //Histogram Output: Event-by-Event
        fHistB = new TH1D( "fHistB", ";b;Count",400,0,20);
        //Keeps track of some basics
        fListHist->Add(fHistB);
    }
    if(! fHistNchVsB ) {
        //Histogram Output: Event-by-Event
        fHistNchVsB = new TH2D( "fHistNchVsB", ";b;Count",400,0,20,lNNchBins,lLowNchBound,lHighNchBound);
        //Keeps track of some basics
        fListHist->Add(fHistNchVsB);
    }
    
    //Identified Particles
    Int_t lPDGCodes[9] = {211, 321, 2212, 310, 3122, 3312, 3334, 333, 313};
    TString lPartNames[9] = {
        "Pion", "Kaon", "Proton", "K0Short", "Lambda", "Xi", "Omega", "Phi", "KStar"
    };

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
        if(! fHistPtVsNpart[ih] ) {
            fHistPtVsNpart[ih] = new TH2D(Form("fHistPtVsNpart_%s",lPartNames[ih].Data()),    "Generated;p_{T} (GeV/c)",500,-0.5,499.5,lNPtBins,0,lMaxPt);
            fListHist->Add(fHistPtVsNpart[ih]);
        }
        if(! fHistPtVsB[ih] ) {
            fHistPtVsB[ih] = new TH2D(Form("fHistPtVsB_%s",lPartNames[ih].Data()),    "Generated;p_{T} (GeV/c)",400,0,20,lNPtBins,0,lMaxPt);
            fListHist->Add(fHistPtVsB[ih]);
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
    // Acquire information on Npart, Ncoll, b
    //------------------------------------------------
    
    //Npart and Ncoll information
    AliGenHijingEventHeader* hHijing=0;
    AliGenDPMjetEventHeader* dpmHeader=0;
    AliGenEventHeader* mcGenH = lMCevent->GenEventHeader();
    
    Int_t fMC_NPart = -1;
    Int_t fMC_NColl = -1;
    Float_t fMC_b = -1;
    
    //DPMJet/HIJING info if available
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
    
    //check EPOS info, if available
    if ( IsEPOSLHC() ){
        AliGenHepMCEventHeader *lHepMCHeader = 0x0;
        if (mcGenH->InheritsFrom(AliGenHepMCEventHeader::Class()))
        lHepMCHeader = (AliGenHepMCEventHeader*)mcGenH;
        
        if (lHepMCHeader ){
            fMC_NPart = lHepMCHeader->Npart_proj()+lHepMCHeader->Npart_targ();
            fMC_NColl = lHepMCHeader->N_Nwounded_collisions() +
            lHepMCHeader->Nwounded_N_collisions() +
            lHepMCHeader->Nwounded_Nwounded_collisions();
            
            fMC_b = lHepMCHeader->impact_parameter();
        }
    }
    
    //------------------------------------------------
    // Fill Event Counters
    //------------------------------------------------

    //Basics: All Processed
    fHistEventCounter->Fill(0.5);
    
    fHistV0MMult        -> Fill ( lNchVZEROA+lNchVZEROC );
    fHistNchVsV0MMult   -> Fill ( lNchVZEROA+lNchVZEROC, lNchEta5  );
    fHistNpart          -> Fill ( fMC_NPart );
    fHistNchVsNpart     -> Fill ( fMC_NPart, lNchEta5  );
    fHistB              -> Fill ( fMC_b );
    fHistNchVsB         -> Fill ( fMC_b, lNchEta5  );
    
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
                fHistPtVsNpart[ih]->Fill(fMC_NPart,lThisPt);
                fHistPtVsB[ih]->Fill(fMC_b,lThisPt);
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

//______________________________________________________________________
Double_t AliAnalysisTaskMCPredictions::MyRapidity(Double_t rE, Double_t rPz) const
{
    // Local calculation for rapidity
    Double_t ReturnValue = -100;
    if( (rE-rPz+1.e-13) != 0 && (rE+rPz) != 0 ) {
        ReturnValue =  0.5*TMath::Log((rE+rPz)/(rE-rPz+1.e-13));
    }
    return ReturnValue;
}


//______________________________________________________________________
Bool_t AliAnalysisTaskMCPredictions::IsHijing() const {
    //Function to check if this is Hijing MC
    Bool_t lReturnValue = kFALSE;
    AliMCEvent*  mcEvent = MCEvent();
    if (mcEvent) {
        AliGenEventHeader* mcGenH = mcEvent->GenEventHeader();
        if (mcGenH->InheritsFrom(AliGenHijingEventHeader::Class())){
            //Option 1: Just Hijing
            lReturnValue = kTRUE;
        } else if (mcGenH->InheritsFrom(AliGenCocktailEventHeader::Class())) {
            //Option 2: cocktail involving Hijing
            TList* headers = ((AliGenCocktailEventHeader*)mcGenH)->GetHeaders();
            TIter next(headers);
            while (const TObject *obj=next()){
                //Look for an object inheriting from the hijing header class
                if ( obj->InheritsFrom(AliGenHijingEventHeader::Class()) ){ lReturnValue = kTRUE; }
            }
        }
    }
    return lReturnValue;
}

//______________________________________________________________________
Bool_t AliAnalysisTaskMCPredictions::IsDPMJet() const {
    //Function to check if this is DPMJet
    Bool_t lReturnValue = kFALSE;
    AliMCEvent*  mcEvent = MCEvent();
    if (mcEvent) {
        AliGenEventHeader* mcGenH = mcEvent->GenEventHeader();
        if (mcGenH->InheritsFrom(AliGenDPMjetEventHeader::Class())) {
            //DPMJet Header is there!
            lReturnValue = kTRUE;
        }
    }
    return lReturnValue;
}

//______________________________________________________________________
Bool_t AliAnalysisTaskMCPredictions::IsEPOSLHC() const {
    //Function to check if this is DPMJet
    Bool_t lReturnValue = kFALSE;
    AliMCEvent*  mcEvent = MCEvent();
    if (mcEvent) {
        AliGenEventHeader* mcGenH = mcEvent->GenEventHeader();
        //A bit uncivilized, but hey, if it works...
        TString lHeaderTitle = mcGenH->GetName();
        if (lHeaderTitle.Contains("EPOSLHC")) {
            //This header has "EPOS" in its title!
            lReturnValue = kTRUE;
        }
    }
    return lReturnValue;
}

