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
#include "AliGenPythiaEventHeader.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskMCPredictions)

AliAnalysisTaskMCPredictions::AliAnalysisTaskMCPredictions()
: AliAnalysisTaskSE(),
fListHist(0),
fHistEventCounter(0),
fHistChargedEta(0),
fSmallMultRange(1000),
fLargeMultRange(2000),
fRebinFactor(1),
fkSelectINELgtZERO(kTRUE),
fHistV0MMult(0),
fHistSPDMult(0),
fHistNchVsV0MMult(0),
fHistNchVsSPDMult(0),
fHistNpart(0),
fHistNchVsNpart(0),
fHistB(0),
fHistNchVsB(0),
fHistNMPI(0),
fHistNchVsNMPI(0),
fkDo2pc(kTRUE),
fMinPtTriggerCharged(0.15),
fMinPtTriggerXi(0.8),
fMinPtTriggerPhi(0.4),
fEtaTriggerCharged(0),
fEtaTriggerXi(0),
fEtaTriggerPhi(0)
{
    for(Int_t ih=0; ih<27; ih++){
        fHistPt[ih]          = 0x0;
        fHistEta[ih]         = 0x0;
        fHistPtVsV0MMult[ih] = 0x0;
        fHistPtVsSPDMult[ih] = 0x0;
        fHistPtVsNpart[ih]   = 0x0;
        fHistPtVsB[ih]       = 0x0;
        fHistPtVsNMPI[ih]   = 0x0;
        fHist3d2pcSE[ih]     = 0x0;
        fHist3d2pcXiSE[ih]   = 0x0;
        fHist3d2pcPhiSE[ih]  = 0x0;
    }
}

AliAnalysisTaskMCPredictions::AliAnalysisTaskMCPredictions(const char *name, Int_t lNSmallBinning, Int_t lNLargeBinning, Int_t lRebinFactor)
: AliAnalysisTaskSE(name),
fListHist(0),
fHistEventCounter(0),
fHistChargedEta(0),
fSmallMultRange(lNSmallBinning),
fLargeMultRange(lNLargeBinning),
fRebinFactor(lRebinFactor),
fkSelectINELgtZERO(kTRUE),
fHistV0MMult(0),
fHistSPDMult(0),
fHistNchVsV0MMult(0),
fHistNchVsSPDMult(0),
fHistNpart(0),
fHistNchVsNpart(0),
fHistB(0),
fHistNchVsB(0),
fHistNMPI(0),
fHistNchVsNMPI(0),
fkDo2pc(kTRUE),
fMinPtTriggerCharged(0.15),
fMinPtTriggerXi(0.8),
fMinPtTriggerPhi(0.4),
fEtaTriggerCharged(0),
fEtaTriggerXi(0),
fEtaTriggerPhi(0)
{
    for(Int_t ih=0; ih<27; ih++){
        fHistPt[ih]          = 0x0;
        fHistEta[ih]         = 0x0;
        fHistPtVsV0MMult[ih] = 0x0;
        fHistPtVsSPDMult[ih] = 0x0;
        fHistPtVsNpart[ih]   = 0x0;
        fHistPtVsB[ih]       = 0x0;
        fHistPtVsNMPI[ih]   = 0x0;
        fHist3d2pcSE[ih]     = 0x0;
        fHist3d2pcXiSE[ih]   = 0x0;
        fHist3d2pcPhiSE[ih]  = 0x0;
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
    
    Int_t lNEtaBins = 400;
    Double_t lMaxAbsEta = 2;
    
    //Settings for charged particle counters (integers!)
    Int_t lNNchBins = fSmallMultRange/fRebinFactor;
    Double_t lLowNchBound  = -0.5;
    Double_t lHighNchBound = -0.5 + ((double)(lNNchBins));
    
    Int_t lNNchBinsV0M = fLargeMultRange/fRebinFactor;
    Double_t lLowNchBoundV0M  = -0.5;
    Double_t lHighNchBoundV0M = -0.5 + ((double)(lNNchBinsV0M));
    
    if(! fHistEventCounter ) {
        //Histogram Output: Event-by-Event
        fHistEventCounter = new TH1D( "fHistEventCounter", ";Evt. Sel. Step;Count",1,0,1);
        //Keeps track of some basics
        fHistEventCounter->GetXaxis()->SetBinLabel(1, "Processed");
        fListHist->Add(fHistEventCounter);
    }
    if(! fHistChargedEta ) {
        //Histogram Output: Event-by-Event
        fHistChargedEta = new TH1D( "fHistChargedEta", ";#eta;Count",lNEtaBins,-lMaxAbsEta,+lMaxAbsEta);
        fListHist->Add(fHistChargedEta);
    }
    //___________________________________________________
    if(! fHistV0MMult ) {
        //Histogram Output: Event-by-Event
        fHistV0MMult = new TH1D( "fHistV0MMult", ";V0M Mult;Count",lNNchBinsV0M,lLowNchBoundV0M,lHighNchBoundV0M);
        //Keeps track of some basics
        fListHist->Add(fHistV0MMult);
    }
    if(! fHistSPDMult ) {
        //Histogram Output: Event-by-Event
        fHistSPDMult = new TH1D( "fHistSPDMult", ";SPD Mult;Count",lNNchBinsV0M,lLowNchBoundV0M,lHighNchBoundV0M);
        //Keeps track of some basics
        fListHist->Add(fHistSPDMult);
    }
    if(! fHistNchVsV0MMult ) {
        //Histogram Output: Event-by-Event
        fHistNchVsV0MMult = new TH2D( "fHistNchVsV0MMult", ";V0M Mult;Count",
                                     lNNchBinsV0M,lLowNchBoundV0M,lHighNchBoundV0M,
                                     lNNchBins,lLowNchBound,lHighNchBound);
        //Keeps track of some basics
        fListHist->Add(fHistNchVsV0MMult);
    }
    if(! fHistNchVsSPDMult ) {
        //Histogram Output: Event-by-Event
        fHistNchVsSPDMult = new TH2D( "fHistNchVsSPDMult", ";SPD Mult;Count",
                                     lNNchBinsV0M,lLowNchBoundV0M,lHighNchBoundV0M,
                                     lNNchBins,lLowNchBound,lHighNchBound);
        //Keeps track of some basics
        fListHist->Add(fHistNchVsSPDMult);
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
    //___________________________________________________
    if(! fHistNMPI ) {
        //Histogram Output: Event-by-Event
        fHistNMPI = new TH1D( "fHistNMPI", ";N_{MPI};Count",50,-0.5,49.5);
        //Keeps track of some basics
        fListHist->Add(fHistNMPI);
    }
    if(! fHistNchVsNMPI ) {
        //Histogram Output: Event-by-Event
        fHistNchVsNMPI = new TH2D( "fHistNchVsNMPI", ";N_{part};Count",50,-0.5,49.5,lNNchBins,lLowNchBound,lHighNchBound);
        //Keeps track of some basics
        fListHist->Add(fHistNchVsNMPI);
    }
    //___________________________________________________
    
    //Identified Particles
    TString lPartNames[27] = {
        "PiPlus", "PiMinus", "KaPlus", "KaMinus", "Proton", "AntiProton",
        "K0Short", "Lambda", "AntiLambda",
        "XiMinus", "XiPlus", "OmegaMinus", "OmegaPlus",
        "Phi", "KStar", "AntiKStar",
        "D0", "AntiD0", "D0s", "AntiD0s",
        "Lambdac", "AntiLambdac", "JPsi",
        "Pi0", "AntiPi0", "Eta", "AntiEta"
    };

    //Main Output: Histograms
    
    //Event counter histogram: Multiplicity, Npart, b (if available)
    for(Int_t ih=0; ih<27; ih++){
        if(! fHistPt[ih] ) {
            fHistPt[ih] = new TH1D(Form("fHistPt_%s",lPartNames[ih].Data()),    "Generated;p_{T} (GeV/c)",lNPtBins,0,lMaxPt);
            fListHist->Add(fHistPt[ih]);
        }
    }
    for(Int_t ih=0; ih<27; ih++){
        if(! fHistEta[ih] ) {
            fHistEta[ih] = new TH1D(Form("fHistEta_%s",lPartNames[ih].Data()),    "Generated;#eta",lNEtaBins,-lMaxAbsEta,+lMaxAbsEta);
            fListHist->Add(fHistEta[ih]);
        }
    }
    for(Int_t ih=0; ih<27; ih++){
        if(! fHistPtVsV0MMult[ih] ) {
            fHistPtVsV0MMult[ih] = new TH2D(Form("fHistPtVsV0MMult_%s",lPartNames[ih].Data()),    "Generated;p_{T} (GeV/c)",lNNchBinsV0M,lLowNchBoundV0M,lHighNchBoundV0M,lNPtBins,0,lMaxPt);
            fListHist->Add(fHistPtVsV0MMult[ih]);
        }
    }
    for(Int_t ih=0; ih<27; ih++){
        if(! fHistPtVsSPDMult[ih] ) {
            fHistPtVsSPDMult[ih] = new TH2D(Form("fHistPtVsSPDMult_%s",lPartNames[ih].Data()),    "Generated;p_{T} (GeV/c)",lNNchBinsV0M,lLowNchBoundV0M,lHighNchBoundV0M,lNPtBins,0,lMaxPt);
            fListHist->Add(fHistPtVsSPDMult[ih]);
        }
    }
    for(Int_t ih=0; ih<27; ih++){
        if(! fHistPtVsNpart[ih] ) {
            fHistPtVsNpart[ih] = new TH2D(Form("fHistPtVsNpart_%s",lPartNames[ih].Data()),    "Generated;p_{T} (GeV/c)",500,-0.5,499.5,lNPtBins,0,lMaxPt);
            fListHist->Add(fHistPtVsNpart[ih]);
        }
    }
    for(Int_t ih=0; ih<27; ih++){
        if(! fHistPtVsB[ih] ) {
            fHistPtVsB[ih] = new TH2D(Form("fHistPtVsB_%s",lPartNames[ih].Data()),    "Generated;p_{T} (GeV/c)",400,0,20,lNPtBins,0,lMaxPt);
            fListHist->Add(fHistPtVsB[ih]);
        }
    }
    for(Int_t ih=0; ih<27; ih++){
        if(! fHistPtVsNMPI[ih] ) {
            fHistPtVsNMPI[ih] = new TH2D(Form("fHistPtVsNMPI_%s",lPartNames[ih].Data()),    "Generated;p_{T} (GeV/c)",50,-0.5,49.5,lNPtBins,0,lMaxPt);
            fListHist->Add(fHistPtVsNMPI[ih]);
        }
    }
    
    //2pc histograms
    for(Int_t ih=0; ih<27; ih++){
        if(! fHist3d2pcSE[ih] ) {
            fHist3d2pcSE[ih] = new TH3D(Form("fHist3d2pcSE_%s",lPartNames[ih].Data()),"",64,-1.6,1.6,80,-0.5*TMath::Pi(), 1.5*TMath::Pi(),lNPtBins,0,lMaxPt);
            fListHist->Add(fHist3d2pcSE[ih]);
        }
    }
    //2pc histograms
    for(Int_t ih=0; ih<27; ih++){
        if(! fHist3d2pcPhiSE[ih] ) {
            fHist3d2pcPhiSE[ih] = new TH3D(Form("fHist3d2pcPhiSE_%s",lPartNames[ih].Data()),"",64,-1.6,1.6,80,-0.5*TMath::Pi(), 1.5*TMath::Pi(),lNPtBins,0,lMaxPt);
            fListHist->Add(fHist3d2pcPhiSE[ih]);
        }
    }
    //2pc histograms
    for(Int_t ih=0; ih<27; ih++){
        if(! fHist3d2pcXiSE[ih] ) {
            fHist3d2pcXiSE[ih] = new TH3D(Form("fHist3d2pcXiSE_%s",lPartNames[ih].Data()),"",64,-1.6,1.6,80,-0.5*TMath::Pi(), 1.5*TMath::Pi(),lNPtBins,0,lMaxPt);
            fListHist->Add(fHist3d2pcXiSE[ih]);
        }
    }
    if(! fEtaTriggerCharged ) {
        //Histogram Output: Event-by-Event
        fEtaTriggerCharged = new TH1D( "fEtaTriggerCharged", ";#eta;Count",lNEtaBins,-lMaxAbsEta,+lMaxAbsEta);
        fListHist->Add(fEtaTriggerCharged);
    }
    if(! fEtaTriggerXi ) {
        //Histogram Output: Event-by-Event
        fEtaTriggerXi = new TH1D( "fEtaTriggerXi", ";#eta;Count",lNEtaBins,-lMaxAbsEta,+lMaxAbsEta);
        fListHist->Add(fEtaTriggerXi);
    }
    if(! fEtaTriggerPhi ) {
        //Histogram Output: Event-by-Event
        fEtaTriggerPhi = new TH1D( "fEtaTriggerPhi", ";#eta;Count",lNEtaBins,-lMaxAbsEta,+lMaxAbsEta);
        fListHist->Add(fEtaTriggerPhi);
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
    Long_t lNchEta14  = 0;
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
        
        Double_t gpt = particleOne -> Pt();
        Double_t geta = particleOne -> Eta();
        
        //keep track of base eta distribution
        if ( gpt > fMinPtTriggerCharged ) fHistChargedEta->Fill( geta );
        
        if( TMath::Abs(geta) < 0.5 ) lNchEta5++;
        if( TMath::Abs(geta) < 0.8 ) lNchEta8++;
        if( (TMath::Abs(geta) > 0.8) && (TMath::Abs(geta) < 1.5) ) lNchEta8to15++;
        if( TMath::Abs(geta) < 1.0 ) lNchEta10++;
        if( TMath::Abs(geta) < 1.4 ) lNchEta14++;
        if( TMath::Abs(geta) < 1.0 ) lEvSel_INELgtZEROStackPrimaries = kTRUE;
        if( 2.8 < geta && geta < 5.1 ) lNchVZEROA++;
        if(-3.7 < geta && geta <-1.7 ) lNchVZEROC++;
    }//End of loop on tracks
    //----- End Loop on Stack ------------------------------------------------------------
    
    //Reject non-INEL>0 if requested
    if( !lEvSel_INELgtZEROStackPrimaries && fkSelectINELgtZERO ) return;
    
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
    Int_t fMC_NMPI = -1;
    
    if (mcGenH->InheritsFrom(AliGenPythiaEventHeader::Class())){
        AliGenPythiaEventHeader *fMcPythiaHeader = dynamic_cast <AliGenPythiaEventHeader*> (mcGenH);
        if(fMcPythiaHeader){
            fMC_NMPI = fMcPythiaHeader->GetNMPI();
        }
    }
    
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
    fHistSPDMult        -> Fill ( lNchEta14 );
    fHistNchVsV0MMult   -> Fill ( lNchVZEROA+lNchVZEROC, lNchEta5  );
    fHistNchVsSPDMult   -> Fill ( lNchEta14, lNchEta5  );
    fHistNpart          -> Fill ( fMC_NPart );
    fHistNchVsNpart     -> Fill ( fMC_NPart, lNchEta5  );
    fHistB              -> Fill ( fMC_b );
    fHistNchVsB         -> Fill ( fMC_b, lNchEta5  );
    fHistNMPI           -> Fill ( fMC_NMPI );
    fHistNchVsNMPI      -> Fill ( fMC_NMPI, lNchEta5  );
    
    //------------------------------------------------
    // Fill Spectra as Needed
    //------------------------------------------------
    
    //~All relevant PWG-LF Identified Particle Information (for looping)
    Int_t lPDGCodes[27] = {
        211, -211, 321, -321, 2212, -2212,
        310, 3122, -3122,
        3312, -3312, 3334, -3334,
        333, 313, -313,
        421, -421, 431, -431,
        4122, -4122, 443,
        111,-111,221,-221
    };
    TString lPartNames[27] = {
        "PiPlus", "PiMinus", "KaPlus", "KaMinus", "Proton", "AntiProton",
        "K0Short", "Lambda", "AntiLambda",
        "XiMinus", "XiPlus", "OmegaMinus", "OmegaPlus",
        "Phi", "KStar", "AntiKStar",
        "D0", "AntiD0", "D0s", "AntiD0s",
        "Lambdac", "AntiLambdac", "JPsi",
        "Pi0", "AntiPi0", "Eta", "AntiEta"
    };
    Bool_t lCheckIsPhysicalPrimary[27] = {
        kTRUE, kTRUE, kTRUE, kTRUE, kTRUE, kTRUE,
        kTRUE, kTRUE, kTRUE,
        kTRUE, kTRUE, kTRUE, kTRUE,
        kFALSE, kFALSE, kFALSE,
        kFALSE, kFALSE, kFALSE, kFALSE,
        kFALSE, kFALSE, kFALSE,
        kTRUE, kTRUE, kTRUE, kTRUE
    };
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
        for(Int_t ih=0; ih<27; ih++) if( lThisPDG == lPDGCodes[ih] ) lContinue = kFALSE;
        if ( lContinue ) continue;
            
        lThisRap   = MyRapidity(lPart->Energy(),lPart->Pz());
        //lThisRap   = lPart->Y();
        lThisPt    = lPart->Pt();
        
        //Use Physical Primaries only for filling These Histos
        //if ( lMCstack->IsPhysicalPrimary(ilab)!=kTRUE ) continue;
        lIsPhysicalPrimary = lMCstack->IsPhysicalPrimary(ilab);
        
        for(Int_t ih=0; ih<27; ih++){
            if( lThisPDG == lPDGCodes[ih] ) {
                //Check if primary (if needed) and if not don't use this particle
                if( lCheckIsPhysicalPrimary[ih] == kTRUE && lIsPhysicalPrimary == kFALSE ) continue;
                //Fill Histograms
                fHistEta[ih] -> Fill ( lPart -> Eta() );
                if( TMath::Abs(lThisRap) < 0.5 ) {
                    fHistPt[ih]->Fill(lThisPt);
                    fHistPtVsV0MMult[ih]->Fill(lNchVZEROA+lNchVZEROC,lThisPt);
                    fHistPtVsSPDMult[ih]->Fill(lNchEta14,lThisPt);
                    fHistPtVsNpart[ih]->Fill(fMC_NPart,lThisPt);
                    fHistPtVsB[ih]->Fill(fMC_b,lThisPt);
                    fHistPtVsNMPI[ih]->Fill(fMC_NMPI,lThisPt);
                }
            }
        }
    }//End of loop on tracks
    //----- End Loop on Stack ----------------------
    
    //===== Start 2pc nested loops =================
    if( fkDo2pc ) {
        //Apply the eta cut first or go home
        Long_t lNValidParticles = 0;
        TArrayI lValidParticles(lMCstack->GetNtrack());
        Long_t lNValidPhi = 0;
        TArrayI lValidPhi(lMCstack->GetNtrack());
        Long_t lNValidXi = 0;
        TArrayI lValidXi(lMCstack->GetNtrack());
        //----- Determine valid triggers ----------------------------------------------------------------
        for (Int_t iCurrentLabelStack = 0;  iCurrentLabelStack < (lMCstack->GetNtrack()); iCurrentLabelStack++)
        {
            // Determine if within acceptance, otherwise fully reject from list
            // done such that this check is done O(N) and not O(N^2)
            TParticle* lThisParticle = lMCstack->Particle(iCurrentLabelStack);
            if(!lThisParticle) continue;
            Double_t geta = lThisParticle -> Eta();
            if( TMath::Abs(geta)<0.8 ) lValidParticles[lNValidParticles++]=iCurrentLabelStack;
        }
        //----- Loop on Stack ----------------------------------------------------------------
        for (Int_t iCurrentLabelStack = 0;  iCurrentLabelStack < lNValidParticles; iCurrentLabelStack++)
        {   // This is the begining of the loop on tracks
            TParticle* lTriggerParticle = lMCstack->Particle(lValidParticles[iCurrentLabelStack]);
            if(!lTriggerParticle) continue;
            if(!lTriggerParticle->GetPDG()) continue;
            Double_t lThisCharge = lTriggerParticle->GetPDG()->Charge()/3.;
            //if(TMath::Abs(lThisCharge)<0.001) continue;
            //if(! (lMCstack->IsPhysicalPrimary(lValidParticles[iCurrentLabelStack])) ) continue;
            
            Bool_t lTrigIsCharged = kTRUE;
            if( TMath::Abs(lThisCharge)<0.001 ) lTrigIsCharged = kFALSE;
            Bool_t lTrigIsPrimary = kTRUE;
            if ( !lMCstack->IsPhysicalPrimary(lValidParticles[iCurrentLabelStack]) ) lTrigIsPrimary = kFALSE;
            Bool_t lTrigIsPhi = kTRUE;
            if (lTriggerParticle->GetPdgCode()!=333) lTrigIsPhi = kFALSE;
            
            if( ((!lTrigIsCharged)||(!lTrigIsPrimary)) && !lTrigIsPhi ) continue;
            
            Double_t geta = lTriggerParticle -> Eta();
            Double_t gphi = lTriggerParticle -> Phi();
            
            if( lTriggerParticle -> Pt() > fMinPtTriggerCharged && lTrigIsCharged && lTrigIsPrimary )
                fEtaTriggerCharged -> Fill( geta );
            if( lTriggerParticle -> Pt() > fMinPtTriggerXi && lTrigIsPrimary && TMath::Abs(lTriggerParticle->GetPdgCode())==3312 )
                fEtaTriggerXi      -> Fill( geta );
            if( lTriggerParticle -> Pt() > fMinPtTriggerPhi && TMath::Abs(lTriggerParticle->GetPdgCode())==333 )
                fEtaTriggerPhi     -> Fill( geta );
            
            for (Int_t ilab = 0;  ilab < lNValidParticles; ilab++)
            {   // This is the begining of the loop on tracks
                
                if(ilab == iCurrentLabelStack) continue; //remove auto-correlations
                TParticle* lAssociatedParticle = 0x0;
                lAssociatedParticle = lMCstack->Particle( lValidParticles[ilab] );
                if(!lAssociatedParticle) {
                    Printf("Generated loop %d - MC TParticle pointer to current stack particle = 0x0 ! Skip ...\n", ilab );
                    continue;
                }
                
                lThisPDG = lAssociatedParticle->GetPdgCode();
                
                //Continue if this is not a particle of the right PDG Code (avoids y-calculation problems)
                Bool_t lContinue = kTRUE;
                for(Int_t ih=0; ih<27; ih++) if( lThisPDG == lPDGCodes[ih] ) lContinue = kFALSE;
                if ( lContinue ) continue;
                
                Double_t geta2 = lAssociatedParticle -> Eta();
                Double_t gphi2 = lAssociatedParticle -> Phi();

                lThisPt    = lAssociatedParticle->Pt();

                lIsPhysicalPrimary = lMCstack->IsPhysicalPrimary(lValidParticles[ilab]);
                
                if( lTrigIsCharged && lTrigIsPrimary ){
                    for(Int_t ih=0; ih<27; ih++){
                        if( lThisPDG == lPDGCodes[ih] && TMath::Abs(geta2) < 0.8 ) {
                            //Check if primary (if needed) and if not don't use this particle
                            if( lCheckIsPhysicalPrimary[ih] == kTRUE && lIsPhysicalPrimary == kFALSE ) continue;
                            //Fill 2pc same-event histograms, please
                            fHist3d2pcSE[ih]->Fill(geta2-geta, ComputeDeltaPhi(gphi,gphi2), lThisPt) ;
                        }
                    }
                }
                if( lTriggerParticle->GetPdgCode() == 3312 ){
                    for(Int_t ih=0; ih<27; ih++){
                        if( lThisPDG == lPDGCodes[ih] && TMath::Abs(geta2) < 0.8 ) {
                            //Check if primary (if needed) and if not don't use this particle
                            if( lCheckIsPhysicalPrimary[ih] == kTRUE && lIsPhysicalPrimary == kFALSE ) continue;
                            //Fill 2pc same-event histograms, please
                            fHist3d2pcXiSE[ih]->Fill(geta2-geta, ComputeDeltaPhi(gphi,gphi2), lThisPt) ;
                        }
                    }
                }
                if( TMath::Abs( lTriggerParticle->GetPdgCode() ) == 333 ){
                    for(Int_t ih=0; ih<27; ih++){
                        if( lThisPDG == lPDGCodes[ih] && TMath::Abs(geta2) < 0.8 ) {
                            //Check if primary (if needed) and if not don't use this particle
                            if( lCheckIsPhysicalPrimary[ih] == kTRUE && lIsPhysicalPrimary == kFALSE ) continue;
                            //Fill 2pc same-event histograms, please
                            fHist3d2pcPhiSE[ih]->Fill(geta2-geta, ComputeDeltaPhi(gphi,gphi2), lThisPt) ;
                        }
                    }
                }
            }//End of loop on tracks
        }//End of loop on tracks
        //----- End Loop on Stack ------------------------------------------------------------
    }
    //===== End 2pc nested loops ===================
    
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

//______________________________________________________________________
Double_t AliAnalysisTaskMCPredictions::ComputeDeltaPhi( Double_t phi1, Double_t phi2) const {
    //To be completely sure, use inner products
    Double_t x1, y1, x2, y2;
    x1 = TMath::Cos( phi1 );
    y1 = TMath::Sin( phi1 );
    x2 = TMath::Cos( phi2 );
    y2 = TMath::Sin( phi2 );
    Double_t lInnerProd  = x1*x2 + y1*y2;
    Double_t lVectorProd = x1*y2 - x2*y1;
    
    Double_t lReturnVal = 0;
    if( lVectorProd > 1e-8 ){
        lReturnVal = TMath::ACos(lInnerProd);
    }
    if( lVectorProd < -1e-8 ){
        lReturnVal = -TMath::ACos(lInnerProd);
    }
    if( lReturnVal < -TMath::Pi()/2 ) lReturnVal += 2*TMath::Pi();
    return lReturnVal;
}
