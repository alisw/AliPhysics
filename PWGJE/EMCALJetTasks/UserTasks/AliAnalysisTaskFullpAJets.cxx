#include "AliAnalysisTaskFullpAJets.h"

#include <Riostream.h>
#include <ctime>
#include <TString.h>
#include <TChain.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TCanvas.h>
#include <TList.h>
#include <TLorentzVector.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TProfile3D.h>
#include <TRandom.h>
#include <TRandom3.h>

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliStack.h"
#include "AliESDtrackCuts.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliAODEvent.h"
#include "AliMCEvent.h"
#include "AliEmcalJet.h"
#include "AliEMCALGeometry.h"

ClassImp(AliAnalysisTaskFullpAJets)

//________________________________________________________________________
AliAnalysisTaskFullpAJets::AliAnalysisTaskFullpAJets() : 
    AliAnalysisTaskSE(),

    fOutput(0),
    fTrackCuts(0),
    fhTrackPt(0),
    fhTrackEta(0),
    fhTrackPhi(0),
    fhClusterPt(0),
    fhClusterEta(0),
    fhClusterPhi(0),
    fhCentrality(0),
    fhBckgMult(0),
    fhBckgFluc(0),
    fhChargedJetPt(0),
    fhChargedJetPtAreaCut(0),
    fhJetPtEMCal(0),
    fhJetPtEMCalAreaCut(0),
    fhJetPtEMCalAreaCutSignal(0),
    fhJetPtTPC(0),
    fhJetPtTPCAreaCut(0),
    fhJetTPtRhoTotal(0),
    fhJetTPtRhoTotalSignal(0),
    fhJetTPtRhoNoLeading(0),
    fhJetTPtRhoNoLeadingSignal(0),
    fhJetTPt1B(0),
    fhJetTPt1BSignal(0),
    fhEMCalBckg1B(0),
    fhJetTPt1C(0),
    fhEMCalBckg1C(0),
    fhEMCalJet2A(0),
    fhJetTPt2B(0),
    fhEMCalBckg2B(0),
    fhJetTPt3(0),
    fhDeltaPtTotal(0),
    fhDeltaPtNoLeading(0),
    fhDeltaPt1B(0),
    fhDeltaRho01(0),
    fhEMCalCellCounts(0),
    fh020RhoTotal(0),
    fh020RhoNoLeading(0),
    fh020Rho1B(0),
    fh020Rho2B(0),
    fh020Rho3(0),
    fh020JetPtEMCal(0),
    fh020JetPtEMCalAreaCut(0),
    fh020JetPtEMCalAreaCutSignal(0),
    fh020JetTPtRhoTotal(0),
    fh020JetTPtRhoTotalSignal(0),
    fh020JetTPtRhoNoLeading(0),
    fh020JetTPtRhoNoLeadingSignal(0),
    fh020JetTPt1B(0),
    fh020JetTPt1BSignal(0),
    fh020JetTPt1C(0),
    fh020JetTPt2B(0),
    fh020JetTPt3(0),
    fhDeltaPt2B(0),
    fhDeltaPtkT(0),
    fh020DiJetAsy(0),
    fh020RhokT(0),
    fh020EMCalkTClusters(0),
    fh020EMCalAkTJets(0),
    fh020DiJetDeltaPhi(0),
    fhDiJetEMCalLeadingPt(0),
    fhDiJetEMCalLeadingDeltaPhi(0),
    fh020EMCalJet2A(0),
    fhDeltaRho0DiJet(0),
    fh020Rho2BCore(0),
    fh020Rho3NoJets(0),
    fh020Rho3DiJets(0),
    fh020Rho3Perp(0),

    fhTrackEtaPhi(0),
    fhClusterEtaPhi(0),
    fhJetPtArea(0),
    fhRhoTotal(0),
    fhRhoNoLeading(0),
    fhRho1B(0),
    fhRho1C(0),
    fhRho2B(0),
    fhRho3(0),
    fhJetConstituentPt(0),
    fhJetPtCenEMCal(0),
    fhJetPtCenEMCalAreaCut(0),
    fhJetPtCenEMCalAreaCutSignal(0),
    fhJetTPtCenRhoTotal(0),
    fhJetTPtCenRhoTotalSignal(0),
    fhJetTPtCenRhoNoLeading(0),
    fhJetTPtCenRhoNoLeadingSignal(0),
    fhJetTPtCen1B(0),
    fhJetTPtCen1BSignal(0),
    fhJetTPtCen1C(0),
    fhJetTPtCen2B(0),
    fhJetTPtCen3(0),
    fhDiJetCenAsy(0),
    fhDiJetCenDeltaPhi(0),
    fhEMCalCenJet2A(0),
    fhRho2BCore(0),
    fhRho3NoJets(0),
    fhRho3DiJets(0),
    fhRho3Perp(0),

    fhJetTrigR1A(0),

    fpEventMult(0),
    fpRhoTotal(0),
    fpRhoNoLeading(0),
    fpRho1B(0),
    fpRho2B(0),
    fpRho3(0),
    fpRhoScale(0),
    fpRhokT(0),
    fpJetPtRhoTotal(0),
    fpJetPtRhoNoLeading(0),
    fpJetPtRhokT(0),
    fpRhoChargedkT(0),
    fpRhoScalekT(0),
    fpRho2BCore(0),
    fpRho3NoJets(0),
    fpRho3DiJets(0),
    fpRho3Perp(0),

    fpTrackPtProfile(0),
    fpClusterPtProfile(0),

    fIsInitialized(0),
    fRJET(0),
    fnEvents(0),
    fnEventsCharged(0),
    fnDiJetEvents(0),
    fEMCalPhiMin(0),
    fEMCalPhiMax(0),
    fEMCalPhiTotal(0),
    fEMCalEtaMin(0),
    fEMCalEtaMax(0),
    fEMCalEtaTotal(0),
    fEMCalArea(0),
    fTPCPhiMin(0),
    fTPCPhiMax(0),
    fTPCPhiTotal(0),
    fTPCEtaMin(0),
    fTPCEtaMax(0),
    fTPCEtaTotal(0),
    fTPCArea(0),
    fJetR(0),
    fJetAreaCutFrac(0),
    fJetAreaThreshold(0),
    fDeltaRho01(0),
    fnEMCalCells(0),
    fCentralityBins(0),
    fCentralityLow(0),
    fCentralityUp(0),
    fEventCentrality(0),
    fRhoTotal(0),
    fRhoCharged(0),
    fRhokTTotal(0),
    fRhokTCharged(0),
    fRhoAkTTotal(0),
    fEtaProfileBins(0),
    fEtaProfileLow(0),
    fEtaProfileUp(0),
    fEDProfileRBins(0),
    fEDProfileRLow(0),
    fEDProfileRUp(0),
    fEDProfilePtBins(0),
    fEDProfilePtLow(0),
    fEDProfilePtUp(0),
    fEDProfileEtaBins(0),
    fEDProfileEtaLow(0),
    fEDProfileEtaUp(0),
    fnTracks(0),
    fnClusters(0),
    fnAKTFullJets(0),
    fnAKTChargedJets(0),
    fnKTFullJets(0),
    fnKTChargedJets(0),
    fnBckgClusters(0),
    fTPCJetThreshold(0),
    fEMCalJetThreshold(0),
    fVertexWindow(0),
    fVertexMaxR(0),
    fnJetsPtCut(0),
    fnJetsPtTPCCut(0),
    fnJetsPtTotalCut(0),
    fnJetsChargedPtCut(0),
    fnJetskTEMCalFull(0),
    fnJetskTTPCFull(0),
    fPtMaxID(0),
    fPtFullMaxID(0),
    fPtTPCMaxID(0),
    fPtFullTPCMaxID(0),
    fPtTotalMaxID(0),
    fPtChargedMaxID(0),
    fPtMax(0),
    fPtFullMax(0),
    fPtTPCMax(0),
    fPtFullTPCMax(0),
    fPtTotalMax(0),
    fPtChargedMax(0),
    fChargedBackJetID(0),
    fmyTracks(0),
    fmyClusters(0),
    fmyAKTFullJets(0),
    fmyAKTChargedJets(0),
    fmyKTFullJets(0),
    fmyKTChargedJets(0),
    fJetPtCutID(0),
    fJetPtTPCCutID(0),
    fJetPtTotalCutID(0),
    fJetPtChargedCutID(0),
    fJetkTEMCalFullID(0),
    fJetkTTPCFullID(0),
    fInEMCal(0),
    fInEMCalFull(0),
    fInTPCFull(0),
    fInTPCChargedFull(0),
    fRCBckgFluc(0)
{
    // Dummy constructor ALWAYS needed for I/O.
    fpJetEtaProfile = new TProfile *[14];
    fpJetAbsEtaProfile = new TProfile *[14];
    fpChargedJetRProfile = new TProfile *[8];
    fpJetRProfile= new TProfile *[4];
    fpChargedJetEDProfile= new TProfile3D *[10];
    fpJetEDProfile= new TProfile3D *[10];
    fvertex[0]=0.0,fvertex[1]=0.0,fvertex[2]=0.0;
}

//________________________________________________________________________
AliAnalysisTaskFullpAJets::AliAnalysisTaskFullpAJets(const char *name) :
    AliAnalysisTaskSE(name),

    fOutput(0),
    fTrackCuts(0),
    fhTrackPt(0),
    fhTrackEta(0),
    fhTrackPhi(0),
    fhClusterPt(0),
    fhClusterEta(0),
    fhClusterPhi(0),
    fhCentrality(0),
    fhBckgMult(0),
    fhBckgFluc(0),
    fhChargedJetPt(0),
    fhChargedJetPtAreaCut(0),
    fhJetPtEMCal(0),
    fhJetPtEMCalAreaCut(0),
    fhJetPtEMCalAreaCutSignal(0),
    fhJetPtTPC(0),
    fhJetPtTPCAreaCut(0),
    fhJetTPtRhoTotal(0),
    fhJetTPtRhoTotalSignal(0),
    fhJetTPtRhoNoLeading(0),
    fhJetTPtRhoNoLeadingSignal(0),
    fhJetTPt1B(0),
    fhJetTPt1BSignal(0),
    fhEMCalBckg1B(0),
    fhJetTPt1C(0),
    fhEMCalBckg1C(0),
    fhEMCalJet2A(0),
    fhJetTPt2B(0),
    fhEMCalBckg2B(0),
    fhJetTPt3(0),
    fhDeltaPtTotal(0),
    fhDeltaPtNoLeading(0),
    fhDeltaPt1B(0),
    fhDeltaRho01(0),
    fhEMCalCellCounts(0),
    fh020RhoTotal(0),
    fh020RhoNoLeading(0),
    fh020Rho1B(0),
    fh020Rho2B(0),
    fh020Rho3(0),
    fh020JetPtEMCal(0),
    fh020JetPtEMCalAreaCut(0),
    fh020JetPtEMCalAreaCutSignal(0),
    fh020JetTPtRhoTotal(0),
    fh020JetTPtRhoTotalSignal(0),
    fh020JetTPtRhoNoLeading(0),
    fh020JetTPtRhoNoLeadingSignal(0),
    fh020JetTPt1B(0),
    fh020JetTPt1BSignal(0),
    fh020JetTPt1C(0),
    fh020JetTPt2B(0),
    fh020JetTPt3(0),
    fhDeltaPt2B(0),
    fhDeltaPtkT(0),
    fh020DiJetAsy(0),
    fh020RhokT(0),
    fh020EMCalkTClusters(0),
    fh020EMCalAkTJets(0),
    fh020DiJetDeltaPhi(0),
    fhDiJetEMCalLeadingPt(0),
    fhDiJetEMCalLeadingDeltaPhi(0),
    fh020EMCalJet2A(0),
    fhDeltaRho0DiJet(0),
    fh020Rho2BCore(0),
    fh020Rho3NoJets(0),
    fh020Rho3DiJets(0),
    fh020Rho3Perp(0),

    fhTrackEtaPhi(0),
    fhClusterEtaPhi(0),
    fhJetPtArea(0),
    fhRhoTotal(0),
    fhRhoNoLeading(0),
    fhRho1B(0),
    fhRho1C(0),
    fhRho2B(0),
    fhRho3(0),
    fhJetConstituentPt(0),
    fhJetPtCenEMCal(0),
    fhJetPtCenEMCalAreaCut(0),
    fhJetPtCenEMCalAreaCutSignal(0),
    fhJetTPtCenRhoTotal(0),
    fhJetTPtCenRhoTotalSignal(0),
    fhJetTPtCenRhoNoLeading(0),
    fhJetTPtCenRhoNoLeadingSignal(0),
    fhJetTPtCen1B(0),
    fhJetTPtCen1BSignal(0),
    fhJetTPtCen1C(0),
    fhJetTPtCen2B(0),
    fhJetTPtCen3(0),
    fhDiJetCenAsy(0),
    fhDiJetCenDeltaPhi(0),
    fhEMCalCenJet2A(0),
    fhRho2BCore(0),
    fhRho3NoJets(0),
    fhRho3DiJets(0),
    fhRho3Perp(0),

    fhJetTrigR1A(0),

    fpEventMult(0),
    fpRhoTotal(0),
    fpRhoNoLeading(0),
    fpRho1B(0),
    fpRho2B(0),
    fpRho3(0),
    fpRhoScale(0),
    fpRhokT(0),
    fpJetPtRhoTotal(0),
    fpJetPtRhoNoLeading(0),
    fpJetPtRhokT(0),
    fpRhoChargedkT(0),
    fpRhoScalekT(0),
    fpRho2BCore(0),
    fpRho3NoJets(0),
    fpRho3DiJets(0),
    fpRho3Perp(0),

    fpTrackPtProfile(0),
    fpClusterPtProfile(0),

    fIsInitialized(0),
    fRJET(0),
    fnEvents(0),
    fnEventsCharged(0),
    fnDiJetEvents(0),
    fEMCalPhiMin(0),
    fEMCalPhiMax(0),
    fEMCalPhiTotal(0),
    fEMCalEtaMin(0),
    fEMCalEtaMax(0),
    fEMCalEtaTotal(0),
    fEMCalArea(0),
    fTPCPhiMin(0),
    fTPCPhiMax(0),
    fTPCPhiTotal(0),
    fTPCEtaMin(0),
    fTPCEtaMax(0),
    fTPCEtaTotal(0),
    fTPCArea(0),
    fJetR(0),
    fJetAreaCutFrac(0),
    fJetAreaThreshold(0),
    fDeltaRho01(0),
    fnEMCalCells(0),
    fCentralityBins(0),
    fCentralityLow(0),
    fCentralityUp(0),
    fEventCentrality(0),
    fRhoTotal(0),
    fRhoCharged(0),
    fRhokTTotal(0),
    fRhokTCharged(0),
    fRhoAkTTotal(0),
    fEtaProfileBins(0),
    fEtaProfileLow(0),
    fEtaProfileUp(0),
    fEDProfileRBins(0),
    fEDProfileRLow(0),
    fEDProfileRUp(0),
    fEDProfilePtBins(0),
    fEDProfilePtLow(0),
    fEDProfilePtUp(0),
    fEDProfileEtaBins(0),
    fEDProfileEtaLow(0),
    fEDProfileEtaUp(0),
    fnTracks(0),
    fnClusters(0),
    fnAKTFullJets(0),
    fnAKTChargedJets(0),
    fnKTFullJets(0),
    fnKTChargedJets(0),
    fnBckgClusters(0),
    fTPCJetThreshold(0),
    fEMCalJetThreshold(0),
    fVertexWindow(0),
    fVertexMaxR(0),
    fnJetsPtCut(0),
    fnJetsPtTPCCut(0),
    fnJetsPtTotalCut(0),
    fnJetsChargedPtCut(0),
    fnJetskTEMCalFull(0),
    fnJetskTTPCFull(0),
    fPtMaxID(0),
    fPtFullMaxID(0),
    fPtTPCMaxID(0),
    fPtFullTPCMaxID(0),
    fPtTotalMaxID(0),
    fPtChargedMaxID(0),
    fPtMax(0),
    fPtFullMax(0),
    fPtTPCMax(0),
    fPtFullTPCMax(0),
    fPtTotalMax(0),
    fPtChargedMax(0),
    fChargedBackJetID(0),
    fmyTracks(0),
    fmyClusters(0),
    fmyAKTFullJets(0),
    fmyAKTChargedJets(0),
    fmyKTFullJets(0),
    fmyKTChargedJets(0),
    fJetPtCutID(0),
    fJetPtTPCCutID(0),
    fJetPtTotalCutID(0),
    fJetPtChargedCutID(0),
    fJetkTEMCalFullID(0),
    fJetkTTPCFullID(0),
    fInEMCal(0),
    fInEMCalFull(0),
    fInTPCFull(0),
    fInTPCChargedFull(0),
    fRCBckgFluc(0)
{
    // Constructor
    // Define input and output slots here (never in the dummy constructor)
    // Input slot #0 works with a TChain - it is connected to the default input container
    // Output slot #1 writes into a TH1 container
    fpJetEtaProfile = new TProfile *[14];
    fpJetAbsEtaProfile = new TProfile *[14];
    fpChargedJetRProfile = new TProfile *[8];
    fpJetRProfile = new TProfile *[4];
    fpChargedJetEDProfile= new TProfile3D *[10];
    fpJetEDProfile= new TProfile3D *[10];
    fvertex[0]=0.0,fvertex[1]=0.0,fvertex[2]=0.0;

    DefineOutput(1,TList::Class());  // for output list
}

//________________________________________________________________________
AliAnalysisTaskFullpAJets::~AliAnalysisTaskFullpAJets()
{
  // Destructor. Clean-up the output list, but not the histograms that are put inside
  // (the list is owner and will clean-up these histograms). Protect in PROOF case.
    if (fOutput && !AliAnalysisManager::GetAnalysisManager()->IsProofMode())
    {
        delete fOutput;
    }
    delete fTrackCuts;
}

//________________________________________________________________________
void AliAnalysisTaskFullpAJets::UserCreateOutputObjects()
{
    // Create histograms
    // Called once (on the worker node)
    fIsInitialized=kFALSE;
    fOutput = new TList();
    fOutput->SetOwner();  // IMPORTANT!
   
    fTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE);

    // Initialize Global Variables
    fnEvents=0;
    fnEventsCharged=0;
    fnDiJetEvents=0;
    
    // fRJET=4 -> fJetR=0.4 && fRJET=25 -> fJetR=0.25, but for writing files, should be 4 and 25 respectively
    if (fRJET>10)
    {
        fJetR=(Double_t)fRJET/100.0;
    }
    else
    {
        fJetR=(Double_t)fRJET/10.0;
    }

    fEMCalPhiMin=(80/(double)360)*2*TMath::Pi();
    fEMCalPhiMax=(187/(double)360)*2*TMath::Pi();
    fEMCalPhiTotal= fEMCalPhiMax-fEMCalPhiMin;
    fEMCalEtaMin=-0.7;
    fEMCalEtaMax=0.7;
    fEMCalEtaTotal=fEMCalEtaMax-fEMCalEtaMin;
    fEMCalArea=fEMCalPhiTotal*fEMCalEtaTotal;

    fTPCPhiMin=(0/(double)360)*2*TMath::Pi();
    fTPCPhiMax=(360/(double)360)*2*TMath::Pi();
    fTPCPhiTotal= fTPCPhiMax-fTPCPhiMin;
    fTPCEtaMin=-0.9;
    fTPCEtaMax=0.9;
    fTPCEtaTotal=fTPCEtaMax-fTPCEtaMin;
    fTPCArea=fTPCPhiTotal*fTPCEtaTotal;

    fCentralityBins=10;
    fCentralityLow=0.0;
    fCentralityUp=100.0;
    Int_t CentralityBinMult=10;
    
    fJetAreaCutFrac =0.6;  // Fudge factor for selecting on jets with threshold Area or higher
    fJetAreaThreshold=fJetAreaCutFrac*TMath::Pi()*fJetR*fJetR;
    fTPCJetThreshold=5.0;  //  Threshold required for an Anti-kt jet to be considered a "true" jet
    fEMCalJetThreshold=5.0;
    fVertexWindow=10.0;
    fVertexMaxR=1.0;
    
    fnBckgClusters=TMath::FloorNint(fEMCalArea/(TMath::Pi()*fJetR*fJetR));  // Select the number of RC that is as close as possible to the area of the EMCal.
    fRCBckgFluc = new Double_t[fnBckgClusters];
    for (Int_t i=0;i<fnBckgClusters;i++)
    {
        fRCBckgFluc[i]=0.0;
    }

    fDeltaRho01=0.0;
    fnEMCalCells=12288;  // sMods 1-10 have 24x48 cells, sMods 11&12 have 8x48 cells...
    
    // Histograms
    Int_t JetPtBins = 200;
    Double_t JetPtLow = 0.0;
    Double_t JetPtUp = 200.0;

    Int_t TCBins=100;
    
    // QA Plots
    fhTrackPt = new TH1D("fhTrackPt","p_{T} distribution of tracks in event",10*JetPtBins,JetPtLow,JetPtUp);
    fhTrackPt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fhTrackPt->GetYaxis()->SetTitle("1/N_{Events} dN/dp_{T}");
    fhTrackPt->Sumw2();

    fhTrackPhi = new TH1D("fhTrackPhi","#phi distribution of tracks in event",TCBins,fTPCPhiMin,fTPCPhiMax);
    fhTrackPhi->GetXaxis()->SetTitle("#phi");
    fhTrackPhi->GetYaxis()->SetTitle("1/N_{Events} dN/d#phi");
    fhTrackPhi->Sumw2();

    fhTrackEta = new TH1D("fhTrackEta","#eta distribution of tracks in event",TCBins,fTPCEtaMin,fTPCEtaMax);
    fhTrackEta->GetXaxis()->SetTitle("#eta");
    fhTrackEta->GetYaxis()->SetTitle("1/N_{Events} dN/d#eta");
    fhTrackEta->Sumw2();

    fhTrackEtaPhi = new TH2D("fhTrackEtaPhi","#eta-#phi distribution of tracks in event",TCBins,fTPCEtaMin,fTPCEtaMax,TCBins,fTPCPhiMin,fTPCPhiMax);
    fhTrackEtaPhi->GetXaxis()->SetTitle("#eta");
    fhTrackEtaPhi->GetYaxis()->SetTitle("#phi");
    fhTrackEtaPhi->GetZaxis()->SetTitle("1/N_{Events} dN/d#etad#phi");
    fhTrackEtaPhi->Sumw2();

    fhClusterPt = new TH1D("fhClusterPt","p_{T} distribution of clusters in event",10*JetPtBins,JetPtLow,JetPtUp);
    fhClusterPt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fhClusterPt->GetYaxis()->SetTitle("1/N_{Events} dN/dp_{T}");
    fhClusterPt->Sumw2();
    
    fhClusterPhi = new TH1D("fhClusterPhi","#phi distribution of clusters in event",TCBins,fTPCPhiMin,fTPCPhiMax);
    fhClusterPhi->GetXaxis()->SetTitle("#phi");
    fhClusterPhi->GetYaxis()->SetTitle("1/N_{Events} dN/d#phi");
    fhClusterPhi->Sumw2();
    
    fhClusterEta = new TH1D("fhClusterEta","#eta distribution of clusters in event",TCBins,fTPCEtaMin,fTPCEtaMax);
    fhClusterEta->GetXaxis()->SetTitle("#eta");
    fhClusterEta->GetYaxis()->SetTitle("1/N_{Events} dN/d#eta");
    fhClusterEta->Sumw2();
    
    fhClusterEtaPhi = new TH2D("fhClusterEtaPhi","#eta-#phi distribution of clusters in event",TCBins,fEMCalEtaMin,fEMCalEtaMax,TCBins,fEMCalPhiMin,fEMCalPhiMax);
    fhClusterEtaPhi->GetXaxis()->SetTitle("#eta");
    fhClusterEtaPhi->GetYaxis()->SetTitle("#phi");
    fhClusterEtaPhi->GetZaxis()->SetTitle("1/N_{Events} dN/d#etad#phi");
    fhClusterEtaPhi->Sumw2();

    fhCentrality = new TH1D("fhCentrality","Event Centrality Distribution",fCentralityBins*CentralityBinMult,fCentralityLow,fCentralityUp);
    fhCentrality->GetXaxis()->SetTitle("Centrality (V0M)");
    fhCentrality->GetYaxis()->SetTitle("1/N_{Events}");
    fhCentrality->Sumw2();
    
    fhBckgFluc = new TH1D("fhBckgFluc",Form("p_{T} distribution of Background Clusters in near central events at center of EMCal with R=%g",fJetR),JetPtBins,JetPtLow,JetPtUp);
    fhBckgFluc->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fhBckgFluc->GetYaxis()->SetTitle("1/N_{Events} dN/dp_{T}");
    fhBckgFluc->Sumw2();

    fhBckgMult = new TH1D("fhBckgMult",Form("Multiplicity distribution of Background Clusters in near central events at center of EMCal with R=%g",fJetR),JetPtBins,JetPtLow,JetPtUp);
    fhBckgMult->GetXaxis()->SetTitle("Multiplicity");
    fhBckgMult->GetYaxis()->SetTitle("1/N_{Events}");
    fhBckgMult->Sumw2();
    
    fhJetConstituentPt= new TH2D("fhJetConstituentPt","Jet constituents p_{T} distribution",JetPtBins, JetPtLow, JetPtUp,10*JetPtBins, JetPtLow, JetPtUp);
    fhJetConstituentPt->GetXaxis()->SetTitle("Jet p_{T} (GeV/c)");
    fhJetConstituentPt->GetYaxis()->SetTitle("Constituent p_{T} (GeV/c)");
    fhJetConstituentPt->Sumw2();
    
    fhEMCalCellCounts = new TH1D("fhEMCalCellCounts","Distribtuion of cluster counts across the EMCal",fnEMCalCells,1,fnEMCalCells);
    fhEMCalCellCounts->GetXaxis()->SetTitle("Absoulute Cell Id");
    fhEMCalCellCounts->GetYaxis()->SetTitle("Counts per Event");
    fhEMCalCellCounts->Sumw2();

    // Raw Jet Spectra
    fhChargedJetPt = new TH1D("fhChargedJetPt","Charged Jet p_{T} distribution for reconstructed Jets",JetPtBins, JetPtLow, JetPtUp);
    fhChargedJetPt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fhChargedJetPt->GetYaxis()->SetTitle("1/N_{Events} dN/dp_{T}");
    fhChargedJetPt->Sumw2();
    
    fhChargedJetPtAreaCut = new TH1D("fhChargedJetPtAreaCut"," Charged Jet p_{T} distribution for reconstructed Jets with Standard Area Cut",JetPtBins, JetPtLow, JetPtUp);
    fhChargedJetPtAreaCut->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fhChargedJetPtAreaCut->GetYaxis()->SetTitle("1/N_{Events} dN/dp_{T}");
    fhChargedJetPtAreaCut->Sumw2();

    fhJetPtTPC = new TH1D("fhJetPtTPC","Jet p_{T} distribution for reconstructed Jets",JetPtBins, JetPtLow, JetPtUp);
    fhJetPtTPC->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fhJetPtTPC->GetYaxis()->SetTitle("1/N_{Events} dN/dp_{T}");
    fhJetPtTPC->Sumw2();
    
    fhJetPtTPCAreaCut = new TH1D("fhJetPtTPCAreaCut","Jet p_{T} distribution for reconstructed Jets",JetPtBins, JetPtLow, JetPtUp);
    fhJetPtTPCAreaCut->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fhJetPtTPCAreaCut->GetYaxis()->SetTitle("1/N_{Events} dN/dp_{T}");
    fhJetPtTPCAreaCut->Sumw2();
    
    fhJetPtEMCal = new TH1D("fhJetPtEMCal","Jet p_{T} distribution for reconstructed Jets within the EMCal",JetPtBins, JetPtLow, JetPtUp);
    fhJetPtEMCal->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fhJetPtEMCal->GetYaxis()->SetTitle("1/N_{Events} dN/dp_{T}");
    fhJetPtEMCal->Sumw2();

    fhJetPtEMCalAreaCut = new TH1D("fhJetPtEMCalAreaCut","Jet p_{T} distribution for reconstructed Jets within the EMCal with Standard Area Cut",JetPtBins, JetPtLow, JetPtUp);
    fhJetPtEMCalAreaCut->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fhJetPtEMCalAreaCut->GetYaxis()->SetTitle("1/N_{Events} dN/dp_{T}");
    fhJetPtEMCalAreaCut->Sumw2();

    fhJetPtEMCalAreaCutSignal = new TH1D("fhJetPtEMCalAreaCutSignal","Jet p_{T} distribution for reconstructed Jets within the EMCal with Standard Area Cut and Signal Cut",JetPtBins, JetPtLow, JetPtUp);
    fhJetPtEMCalAreaCutSignal->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fhJetPtEMCalAreaCutSignal->GetYaxis()->SetTitle("1/N_{Events} dN/dp_{T}");
    fhJetPtEMCalAreaCutSignal->Sumw2();

    fh020JetPtEMCal = new TH1D("fh020JetPtEMCal","0-20% Jet p_{T} distribution for reconstructed Jets within the EMCal",JetPtBins, JetPtLow, JetPtUp);
    fh020JetPtEMCal->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fh020JetPtEMCal->GetYaxis()->SetTitle("1/N_{Events} dN/dp_{T}");
    fh020JetPtEMCal->Sumw2();
    
    fh020JetPtEMCalAreaCut = new TH1D("fh020JetPtEMCalAreaCut","0-20% Jet p_{T} distribution for reconstructed Jets within the EMCal with Standard Area Cut",JetPtBins, JetPtLow, JetPtUp);
    fh020JetPtEMCalAreaCut->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fh020JetPtEMCalAreaCut->GetYaxis()->SetTitle("1/N_{Events} dN/dp_{T}");
    fh020JetPtEMCalAreaCut->Sumw2();
    
    fh020JetPtEMCalAreaCutSignal = new TH1D("fh020JetPtEMCalAreaCutSignal","0-20% Jet p_{T} distribution for reconstructed Jets within the EMCal with Standard Area Cut and Signal Cut",JetPtBins, JetPtLow, JetPtUp);
    fh020JetPtEMCalAreaCutSignal->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fh020JetPtEMCalAreaCutSignal->GetYaxis()->SetTitle("1/N_{Events} dN/dp_{T}");
    fh020JetPtEMCalAreaCutSignal->Sumw2();

    
    fhJetPtCenEMCal = new TH2D("fhJetPtCenEMCal","Jet p_{T} distribution for reconstructed Jets within the EMCal vs Centrality",JetPtBins, JetPtLow, JetPtUp, fCentralityBins*CentralityBinMult,fCentralityLow,fCentralityUp);
    fhJetPtCenEMCal->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fhJetPtCenEMCal->GetYaxis()->SetTitle("Centrality (V0M)");
    fhJetPtCenEMCal->GetZaxis()->SetTitle("1/N_{Events} dN/dp_{T}");
    fhJetPtCenEMCal->Sumw2();
    
    fhJetPtCenEMCalAreaCut = new TH2D("fhJetPtCenEMCalAreaCut","Jet p_{T} distribution for reconstructed Jets within the EMCal with Area Cut vs Centrality",JetPtBins, JetPtLow, JetPtUp, fCentralityBins*CentralityBinMult,fCentralityLow,fCentralityUp);
    fhJetPtCenEMCalAreaCut->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fhJetPtCenEMCalAreaCut->GetYaxis()->SetTitle("Centrality (V0M)");
    fhJetPtCenEMCalAreaCut->GetZaxis()->SetTitle("1/N_{Events} dN/dp_{T}");
    fhJetPtCenEMCalAreaCut->Sumw2();

    fhJetPtCenEMCalAreaCutSignal = new TH2D("fhJetPtCenEMCalAreaCutSignal","Jet p_{T} distribution for reconstructed Jets within the EMCal with Area and Signal Cut vs Centrality",JetPtBins, JetPtLow, JetPtUp, fCentralityBins*CentralityBinMult,fCentralityLow,fCentralityUp);
    fhJetPtCenEMCalAreaCutSignal->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fhJetPtCenEMCalAreaCutSignal->GetYaxis()->SetTitle("Centrality (V0M)");
    fhJetPtCenEMCalAreaCutSignal->GetZaxis()->SetTitle("1/N_{Events} dN/dp_{T}");
    fhJetPtCenEMCalAreaCutSignal->Sumw2();
    
    // Jet Area vs pT Distribution
    Int_t JetPtAreaBins=200;
    Double_t JetPtAreaLow=0.0;
    Double_t JetPtAreaUp=2.0;
    
    fhJetPtArea = new TH2D("fhJetPtArea","Jet Area Distribution",JetPtBins, JetPtLow,JetPtUp,JetPtAreaBins,JetPtAreaLow,JetPtAreaUp);
    fhJetPtArea->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fhJetPtArea->GetYaxis()->SetTitle("A_{jet}");
    fhJetPtArea->GetZaxis()->SetTitle("1/N_{Events} dN/dA_{jet}dp_{T}");
    fhJetPtArea->Sumw2();

    Int_t A1_PtBins = 100;
    Double_t A1_PtLow=0.0;
    Double_t A1_PtUp=100.0;
    Int_t A1_TrigBins = 100;
    Double_t A1_TrigLow=0.0;
    Double_t A1_TrigUp=100.0;
    Int_t A1_RBins = 16;
    Double_t A1_RLow=0.4;
    Double_t A1_RUp=2.0;

    fhJetTrigR1A = new TH3D("fhJetTrigR1A","Jet p_{T} distribution for reconstructed Jets vs Trigger Jet p_{T} vs #DeltaR",A1_PtBins, A1_PtLow, A1_PtUp,A1_TrigBins, A1_TrigLow, A1_TrigUp,A1_RBins, A1_RLow, A1_RUp);
    fhJetTrigR1A->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fhJetTrigR1A->GetYaxis()->SetTitle("Trigger Jet p_{T} (GeV/c)");
    fhJetTrigR1A->GetZaxis()->SetTitle("R");
    fhJetTrigR1A->Sumw2();

    // Corrected Jet Spectra
    Int_t JetTPtBins = 250;
    Double_t JetTPtLow=-50.0;
    Double_t JetTPtUp=200.0;
    
    fhJetTPtRhoTotal = new TH1D("fhJetTPtRhoTotal","True Jet p_{T} distribution for reconstructed Jets in the EMCal",JetTPtBins, JetTPtLow, JetTPtUp);
    fhJetTPtRhoTotal->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fhJetTPtRhoTotal->GetYaxis()->SetTitle("1/N_{Events} dN/dp_{T}");
    fhJetTPtRhoTotal->Sumw2();

    fhJetTPtRhoTotalSignal = new TH1D("fhJetTPtRhoTotalSignal","True Jet p_{T} distribution for reconstructed Jets in the EMCal with Signal cut",JetTPtBins, JetTPtLow, JetTPtUp);
    fhJetTPtRhoTotalSignal->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fhJetTPtRhoTotalSignal->GetYaxis()->SetTitle("1/N_{Events} dN/dp_{T}");
    fhJetTPtRhoTotalSignal->Sumw2();

    fhJetTPtRhoNoLeading = new TH1D("fhJetTPtRhoNoLeading","True Jet p_{T} distribution for reconstructed Jets in the EMCal",JetTPtBins, JetTPtLow, JetTPtUp);
    fhJetTPtRhoNoLeading->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fhJetTPtRhoNoLeading->GetYaxis()->SetTitle("1/N_{Events} dN/dp_{T}");
    fhJetTPtRhoNoLeading->Sumw2();

    fhJetTPtRhoNoLeadingSignal = new TH1D("fhJetTPtRhoNoLeadingSignal","True Jet p_{T} distribution for reconstructed Jets in the EMCal with Signal cut",JetTPtBins, JetTPtLow, JetTPtUp);
    fhJetTPtRhoNoLeadingSignal->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fhJetTPtRhoNoLeadingSignal->GetYaxis()->SetTitle("1/N_{Events} dN/dp_{T}");
    fhJetTPtRhoNoLeadingSignal->Sumw2();

    fhJetTPt1B = new TH1D("fhJetTPt1B","True Jet p_{T} distribution for reconstructed Jets in the EMCal",JetTPtBins, JetTPtLow, JetTPtUp);
    fhJetTPt1B->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fhJetTPt1B->GetYaxis()->SetTitle("1/N_{Events} dN/dp_{T}");
    fhJetTPt1B->Sumw2();

    fhJetTPt1BSignal = new TH1D("fhJetTPt1BSignal","True Jet p_{T} distribution for reconstructed Jets in the EMCal",JetTPtBins, JetTPtLow, JetTPtUp);
    fhJetTPt1BSignal->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fhJetTPt1BSignal->GetYaxis()->SetTitle("1/N_{Events} dN/dp_{T}");
    fhJetTPt1BSignal->Sumw2();

    fhJetTPt1C = new TH1D("fhJetTPt1C","True Jet p_{T} distribution for reconstructed Jets in the EMCal",JetPtBins,JetPtLow,JetPtUp);
    fhJetTPt1C->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fhJetTPt1C->GetYaxis()->SetTitle("1/N_{Events} dN/dp_{T}");
    fhJetTPt1C->Sumw2();

    fhJetTPt2B = new TH1D("fhJetTPt2B","True Jet p_{T} distribution for reconstructed Jets in the EMCal with dijet Trigger in TPC",JetTPtBins, JetTPtLow, JetTPtUp);
    fhJetTPt2B->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fhJetTPt2B->GetYaxis()->SetTitle("1/N_{Events} dN/dp_{T}");
    fhJetTPt2B->Sumw2();
    
    fhJetTPt3 = new TH1D("fhJetTPt3","True Charged jet p_{T} distribution for reconstructed Jets in the TPC",JetTPtBins, JetTPtLow, JetTPtUp);
    fhJetTPt3->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fhJetTPt3->GetYaxis()->SetTitle("1/N_{Events} dN/dp_{T}");
    fhJetTPt3->Sumw2();
    
    // 0-20% Centrality Corrected Spectra (fh020JetTPt...)
    fh020JetTPtRhoTotal = new TH1D("fh020JetTPtRhoTotal","0-20% True Jet p_{T} distribution for reconstructed Jets in the EMCal",JetTPtBins, JetTPtLow, JetTPtUp);
    fh020JetTPtRhoTotal->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fh020JetTPtRhoTotal->GetYaxis()->SetTitle("1/N_{Events} dN/dp_{T}");
    fh020JetTPtRhoTotal->Sumw2();
    
    fh020JetTPtRhoTotalSignal = new TH1D("fh020JetTPtRhoTotalSignal","0-20% True Jet p_{T} distribution for reconstructed Jets in the EMCal with Signal cut",JetTPtBins, JetTPtLow, JetTPtUp);
    fh020JetTPtRhoTotalSignal->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fh020JetTPtRhoTotalSignal->GetYaxis()->SetTitle("1/N_{Events} dN/dp_{T}");
    fh020JetTPtRhoTotalSignal->Sumw2();
    
    fh020JetTPtRhoNoLeading = new TH1D("fh020JetTPtRhoNoLeading","0-20% True Jet p_{T} distribution for reconstructed Jets in the EMCal",JetTPtBins, JetTPtLow, JetTPtUp);
    fh020JetTPtRhoNoLeading->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fh020JetTPtRhoNoLeading->GetYaxis()->SetTitle("1/N_{Events} dN/dp_{T}");
    fh020JetTPtRhoNoLeading->Sumw2();
    
    fh020JetTPtRhoNoLeadingSignal = new TH1D("fh020JetTPtRhoNoLeadingSignal","0-20% True Jet p_{T} distribution for reconstructed Jets in the EMCal with Signal cut",JetTPtBins, JetTPtLow, JetTPtUp);
    fh020JetTPtRhoNoLeadingSignal->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fh020JetTPtRhoNoLeadingSignal->GetYaxis()->SetTitle("1/N_{Events} dN/dp_{T}");
    fh020JetTPtRhoNoLeadingSignal->Sumw2();
    
    fh020JetTPt1B = new TH1D("fh020JetTPt1B","0-20% True Jet p_{T} distribution for reconstructed Jets in the EMCal",JetTPtBins, JetTPtLow, JetTPtUp);
    fh020JetTPt1B->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fh020JetTPt1B->GetYaxis()->SetTitle("1/N_{Events} dN/dp_{T}");
    fh020JetTPt1B->Sumw2();
    
    fh020JetTPt1BSignal = new TH1D("fh020JetTPt1BSignal","0-20% True Jet p_{T} distribution for reconstructed Jets in the EMCal",JetTPtBins, JetTPtLow, JetTPtUp);
    fh020JetTPt1BSignal->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fh020JetTPt1BSignal->GetYaxis()->SetTitle("1/N_{Events} dN/dp_{T}");
    fh020JetTPt1BSignal->Sumw2();
    
    fh020JetTPt1C = new TH1D("fh020JetTPt1C","0-20% True Jet p_{T} distribution for reconstructed Jets in the EMCal",JetPtBins,JetPtLow,JetPtUp);
    fh020JetTPt1C->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fh020JetTPt1C->GetYaxis()->SetTitle("1/N_{Events} dN/dp_{T}");
    fh020JetTPt1C->Sumw2();
    
    fh020JetTPt2B = new TH1D("fh020JetTPt2B","0-20% True Jet p_{T} distribution for reconstructed Jets in the EMCal with dijet Trigger in TPC",JetTPtBins, JetTPtLow, JetTPtUp);
    fh020JetTPt2B->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fh020JetTPt2B->GetYaxis()->SetTitle("1/N_{Events} dN/dp_{T}");
    fh020JetTPt2B->Sumw2();
    
    fh020JetTPt3 = new TH1D("fh020JetTPt3","0-20% True Charged jet p_{T} distribution for reconstructed Jets in the TPC",JetTPtBins, JetTPtLow, JetTPtUp);
    fh020JetTPt3->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fh020JetTPt3->GetYaxis()->SetTitle("1/N_{Events} dN/dp_{T}");
    fh020JetTPt3->Sumw2();
    
    // 2D Corrected Spectra (fhJetTPtCen...)
    fhJetTPtCenRhoTotal = new TH2D("fhJetTPtCenRhoTotal","True Jet p_{T} distribution for reconstructed Jets in the EMCal vs Centrality",JetTPtBins, JetTPtLow, JetTPtUp,fCentralityBins*CentralityBinMult,fCentralityLow,fCentralityUp);
    fhJetTPtCenRhoTotal->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fhJetTPtCenRhoTotal->GetYaxis()->SetTitle("Centrality (V0M)");
    fhJetTPtCenRhoTotal->GetZaxis()->SetTitle("1/N_{Events} dN/dp_{T}");
    fhJetTPtCenRhoTotal->Sumw2();
    
    fhJetTPtCenRhoTotalSignal = new TH2D("fhJetTPtCenRhoTotalSignal","True Jet p_{T} distribution for reconstructed Jets in the EMCal with Signal cut vs Centrality",JetTPtBins, JetTPtLow, JetTPtUp,fCentralityBins*CentralityBinMult,fCentralityLow,fCentralityUp);
    fhJetTPtCenRhoTotalSignal->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fhJetTPtCenRhoTotalSignal->GetYaxis()->SetTitle("Centrality (V0M)");
    fhJetTPtCenRhoTotalSignal->GetZaxis()->SetTitle("1/N_{Events} dN/dp_{T}");
    fhJetTPtCenRhoTotalSignal->Sumw2();
    
    fhJetTPtCenRhoNoLeading = new TH2D("fhJetTPtCenRhoNoLeading","True Jet p_{T} distribution for reconstructed Jets in the EMCal vs Centrality",JetTPtBins, JetTPtLow, JetTPtUp,fCentralityBins*CentralityBinMult,fCentralityLow,fCentralityUp);
    fhJetTPtCenRhoNoLeading->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fhJetTPtCenRhoNoLeading->GetYaxis()->SetTitle("Centrality (V0M)");
    fhJetTPtCenRhoNoLeading->GetZaxis()->SetTitle("1/N_{Events} dN/dp_{T}");
    fhJetTPtCenRhoNoLeading->Sumw2();
    
    fhJetTPtCenRhoNoLeadingSignal = new TH2D("fhJetTPtCenRhoNoLeadingSignal","True Jet p_{T} distribution for reconstructed Jets in the EMCal with Signal cut vs Centrality",JetTPtBins, JetTPtLow, JetTPtUp,fCentralityBins*CentralityBinMult,fCentralityLow,fCentralityUp);
    fhJetTPtCenRhoNoLeadingSignal->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fhJetTPtCenRhoNoLeadingSignal->GetYaxis()->SetTitle("Centrality (V0M)");
    fhJetTPtCenRhoNoLeadingSignal->GetZaxis()->SetTitle("1/N_{Events} dN/dp_{T}");
    fhJetTPtCenRhoNoLeadingSignal->Sumw2();
    
    fhJetTPtCen1B = new TH2D("fhJetTPtCen1B","True Jet p_{T} distribution for reconstructed Jets in the EMCal vs Centrality",JetTPtBins, JetTPtLow, JetTPtUp,fCentralityBins*CentralityBinMult,fCentralityLow,fCentralityUp);
    fhJetTPtCen1B->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fhJetTPtCen1B->GetYaxis()->SetTitle("Centrality (V0M)");
    fhJetTPtCen1B->GetZaxis()->SetTitle("1/N_{Events} dN/dp_{T}");
    fhJetTPtCen1B->Sumw2();
    
    fhJetTPtCen1BSignal = new TH2D("fhJetTPtCen1BSignal","True Jet p_{T} distribution for reconstructed Jets in the EMCal vs Centrality",JetTPtBins, JetTPtLow, JetTPtUp,fCentralityBins*CentralityBinMult,fCentralityLow,fCentralityUp);
    fhJetTPtCen1BSignal->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fhJetTPtCen1BSignal->GetYaxis()->SetTitle("Centrality (V0M)");
    fhJetTPtCen1BSignal->GetZaxis()->SetTitle("1/N_{Events} dN/dp_{T}");
    fhJetTPtCen1BSignal->Sumw2();
    
    fhJetTPtCen1C = new TH2D("fhJetTPtCen1C","True Jet p_{T} distribution for reconstructed Jets in the EMCal vs Centrality",JetPtBins,JetPtLow,JetPtUp,fCentralityBins*CentralityBinMult,fCentralityLow,fCentralityUp);
    fhJetTPtCen1C->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fhJetTPtCen1C->GetYaxis()->SetTitle("Centrality (V0M)");
    fhJetTPtCen1C->GetZaxis()->SetTitle("1/N_{Events} dN/dp_{T}");
    fhJetTPtCen1C->Sumw2();
    
    fhJetTPtCen2B = new TH2D("fhJetTPtCen2B","True Jet p_{T} distribution for reconstructed Jets in the EMCal with dijet Trigger in TPC vs Centrality",JetTPtBins, JetTPtLow, JetTPtUp,fCentralityBins*CentralityBinMult,fCentralityLow,fCentralityUp);
    fhJetTPtCen2B->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fhJetTPtCen2B->GetYaxis()->SetTitle("Centrality (V0M)");
    fhJetTPtCen2B->GetZaxis()->SetTitle("1/N_{Events} dN/dp_{T}");
    fhJetTPtCen2B->Sumw2();
    
    fhJetTPtCen3 = new TH2D("fhJetTPtCen3","True Charged jet p_{T} distribution for reconstructed Jets in the TPC vs Centrality",JetTPtBins, JetTPtLow, JetTPtUp,fCentralityBins*CentralityBinMult,fCentralityLow,fCentralityUp);
    fhJetTPtCen3->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fhJetTPtCen3->GetYaxis()->SetTitle("Centrality (V0M)");
    fhJetTPtCen3->GetZaxis()->SetTitle("1/N_{Events} dN/dp_{T}");
    fhJetTPtCen3->Sumw2();
    
    // Cluster Plots
    fhEMCalBckg1B = new TH1D("fhEMCalBckg1B","Cluster p_{T} distribution for reconstructed Tracks and Calocluster at least 0.5 away from all jets in the EMCal",JetPtBins, JetPtLow,JetPtUp);
    fhEMCalBckg1B->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fhEMCalBckg1B->GetYaxis()->SetTitle("1/N_{Events} dN/dp_{T}");
    fhEMCalBckg1B->Sumw2();
    
    fhEMCalBckg1C = new TH1D("fhEMCalBckg1C","Cluster p_{T} distribution for reconstructed Tracks and Calocluster in Transverse area with R=0.4",JetPtBins, JetPtLow, JetPtUp);
    fhEMCalBckg1C->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fhEMCalBckg1C->GetYaxis()->SetTitle("1/N_{Events} dN/dp_{T}");
    fhEMCalBckg1C->Sumw2();

    fhEMCalJet2A = new TH1D("fhEMCalJet2A","Cluster p_{T} distribution for jets within EMCal from di-jets in TPC",A1_PtBins,A1_PtLow,A1_PtUp);
    fhEMCalJet2A->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fhEMCalJet2A->GetYaxis()->SetTitle("1/N_{Events} dN/dp_{T}");
    fhEMCalJet2A->Sumw2();

    fh020EMCalJet2A = new TH1D("fh020EMCalJet2A","0-20% Centrality, Cluster p_{T} distribution for jets within EMCal from di-jets in TPC",A1_PtBins,A1_PtLow,A1_PtUp);
    fh020EMCalJet2A->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fh020EMCalJet2A->GetYaxis()->SetTitle("1/N_{Events} dN/dp_{T}");
    fh020EMCalJet2A->Sumw2();

    fhEMCalCenJet2A = new TH2D("fhEMCalCenJet2A","Cluster p_{T} distribution for jets within EMCal from di-jets in TPC",A1_PtBins,A1_PtLow,A1_PtUp,fCentralityBins*CentralityBinMult,fCentralityLow,fCentralityUp);
    fhEMCalCenJet2A->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fhEMCalCenJet2A->GetYaxis()->SetTitle("Centrality (V0M)");
    fhEMCalCenJet2A->GetZaxis()->SetTitle("1/N_{Events} dN/dp_{T}");
    fhEMCalCenJet2A->Sumw2();

    fhEMCalBckg2B = new TH1D("fhEMCalBckg2B","Cluster p_{T} distribution for reconstructed Tracks and Calocluster with dijet Trigger in TPC with R=0.4",JetPtBins, JetPtLow, JetPtUp);
    fhEMCalBckg2B->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fhEMCalBckg2B->GetYaxis()->SetTitle("1/N_{Events} dN/dp_{T}");
    fhEMCalBckg2B->Sumw2();
    
    
    Int_t EMCalClusterBins=100;
    Int_t EMCalClusterLow=0;
    Int_t EMCalClusterUp=100;
    
    fh020EMCalkTClusters = new TH1D("fh020EMCalkTClusters","0-20 % Centrality, Number of k_{T} clusters per event",EMCalClusterBins,EMCalClusterLow,EMCalClusterUp);
    fh020EMCalkTClusters->GetXaxis()->SetTitle("# of Clusters");
    fh020EMCalkTClusters->GetYaxis()->SetTitle("1/N_{Events}");
    fh020EMCalkTClusters->Sumw2();

    fh020EMCalAkTJets = new TH1D("fh020EMCalAkTJets","0-20 % Centrality, Number of anti-k_{T} jets per event",EMCalClusterBins,EMCalClusterLow,EMCalClusterUp);
    fh020EMCalAkTJets->GetXaxis()->SetTitle("# of jets");
    fh020EMCalAkTJets->GetYaxis()->SetTitle("1/N_{Events}");
    fh020EMCalAkTJets->Sumw2();

    // Background Density Plots
    // 2D Rho plots (fhRho...)
    Int_t RhoPtBins = 500;
    Double_t RhoPtLow=0.0;
    Double_t RhoPtUp=50.0;
    
    Int_t DeltaRhoPtBins=100;
    Double_t DeltaRhoPtLow=-5.0;
    Double_t DeltaRhoPtUp=5.0;
    
    fhRhoTotal= new TH2D("fhRhoTotal","Background Density #rho_{0}",RhoPtBins,RhoPtLow,RhoPtUp,fCentralityBins*CentralityBinMult,fCentralityLow,fCentralityUp);
    fhRhoTotal->GetXaxis()->SetTitle("p_{T}/Area (GeV/c)");
    fhRhoTotal->GetYaxis()->SetTitle("Centrality (V0M)");
    fhRhoTotal->GetZaxis()->SetTitle("1/N_{Events} dN/d#rho");
    fhRhoTotal->Sumw2();

    fhRhoNoLeading= new TH2D("fhRhoNoLeading","Background Density #rho_{1}",RhoPtBins,RhoPtLow,RhoPtUp,fCentralityBins*CentralityBinMult,fCentralityLow,fCentralityUp);
    fhRhoNoLeading->GetXaxis()->SetTitle("p_{T}/Area (GeV/c)");
    fhRhoNoLeading->GetYaxis()->SetTitle("Centrality (V0M)");
    fhRhoNoLeading->GetZaxis()->SetTitle("1/N_{Events} dN/d#rho");
    fhRhoNoLeading->Sumw2();

    fhRho1B = new TH2D("fhRho1B","Background Density #rho_{n} ",RhoPtBins, RhoPtLow, RhoPtUp,fCentralityBins*CentralityBinMult,fCentralityLow,fCentralityUp);
    fhRho1B->GetXaxis()->SetTitle("p_{T}/Area (GeV/c)");
    fhRho1B->GetYaxis()->SetTitle("Centrality (V0M)");
    fhRho1B->GetZaxis()->SetTitle("1/N_{Events} dN/d#rho");
    fhRho1B->Sumw2();

    fhRho1C = new TH2D("fhRho1C","Background Density #rho (Method 1C)",RhoPtBins, RhoPtLow, RhoPtUp,fCentralityBins*CentralityBinMult,fCentralityLow,fCentralityUp);
    fhRho1C->GetXaxis()->SetTitle("p_{T}/Area (GeV/c)");
    fhRho1C->GetYaxis()->SetTitle("Centrality (V0M)");
    fhRho1C->GetZaxis()->SetTitle("1/N_{Events} dN/d#rho");
    fhRho1C->Sumw2();

    fhRho2B = new TH2D("fhRho2B","Background Density #rho_{dijet}",RhoPtBins, RhoPtLow, RhoPtUp,fCentralityBins*CentralityBinMult,fCentralityLow,fCentralityUp);
    fhRho2B->GetXaxis()->SetTitle("p_{T}/Area (GeV/c)");
    fhRho2B->GetYaxis()->SetTitle("Centrality (V0M)");
    fhRho2B->GetZaxis()->SetTitle("1/N_{Events} dN/d#rho");
    fhRho2B->Sumw2();

    fhRho2BCore = new TH2D("fhRho2BCore","Background Density #rho_{dijet}",RhoPtBins, RhoPtLow, RhoPtUp,fCentralityBins*CentralityBinMult,fCentralityLow,fCentralityUp);
    fhRho2BCore->GetXaxis()->SetTitle("p_{T}/Area (GeV/c)");
    fhRho2BCore->GetYaxis()->SetTitle("Centrality (V0M)");
    fhRho2BCore->GetZaxis()->SetTitle("1/N_{Events} dN/d#rho");
    fhRho2BCore->Sumw2();

    fhRho3 = new TH2D("fhRho3","Charged Background Density #rho_{char}",RhoPtBins, RhoPtLow, RhoPtUp,fCentralityBins*CentralityBinMult,fCentralityLow,fCentralityUp);
    fhRho3->GetXaxis()->SetTitle("p_{T}/Area (GeV/c)");
    fhRho3->GetYaxis()->SetTitle("Centrality (V0M)");
    fhRho3->GetZaxis()->SetTitle("1/N_{Events} dN/d#rho");
    fhRho3->Sumw2();

    fhRho3NoJets = new TH2D("fhRho3NoJets","Charged Background Density #rho_{char} for Events with No Signal Jets",RhoPtBins, RhoPtLow, RhoPtUp,fCentralityBins*CentralityBinMult,fCentralityLow,fCentralityUp);
    fhRho3NoJets->GetXaxis()->SetTitle("p_{T}/Area (GeV/c)");
    fhRho3NoJets->GetYaxis()->SetTitle("Centrality (V0M)");
    fhRho3NoJets->GetZaxis()->SetTitle("1/N_{Events} dN/d#rho");
    fhRho3NoJets->Sumw2();

    fhRho3DiJets = new TH2D("fhRho3DiJets","Charged Background Density #rho_{char} for Events with a DiJet",RhoPtBins, RhoPtLow, RhoPtUp,fCentralityBins*CentralityBinMult,fCentralityLow,fCentralityUp);
    fhRho3DiJets->GetXaxis()->SetTitle("p_{T}/Area (GeV/c)");
    fhRho3DiJets->GetYaxis()->SetTitle("Centrality (V0M)");
    fhRho3DiJets->GetZaxis()->SetTitle("1/N_{Events} dN/d#rho");
    fhRho3DiJets->Sumw2();

    fhRho3Perp = new TH2D("fhRho3Perp","Charged Background Density #rho_{char} for Events with a DiJet at Median Angle Between Dijets",RhoPtBins, RhoPtLow, RhoPtUp,fCentralityBins*CentralityBinMult,fCentralityLow,fCentralityUp);
    fhRho3Perp->GetXaxis()->SetTitle("p_{T}/Area (GeV/c)");
    fhRho3Perp->GetYaxis()->SetTitle("Centrality (V0M)");
    fhRho3Perp->GetZaxis()->SetTitle("1/N_{Events} dN/d#rho");
    fhRho3Perp->Sumw2();

    // 0-20% Centrality Plots (rh020Rho...)
    fh020RhoTotal = new TH1D("fh020RhoTotal","0-20% Background Density #rho_{0}",RhoPtBins,RhoPtLow,RhoPtUp);
    fh020RhoTotal->GetXaxis()->SetTitle("p_{T}/Area (GeV/c)");
    fh020RhoTotal->GetYaxis()->SetTitle("1/N_{Events} dN/d#rho");
    fh020RhoTotal->Sumw2();

    fh020RhoNoLeading = new TH1D("fh020RhoNoLeading","0-20% Background Density #rho_{1}",RhoPtBins,RhoPtLow,RhoPtUp);
    fh020RhoNoLeading->GetXaxis()->SetTitle("p_{T}/Area (GeV/c)");
    fh020RhoNoLeading->GetYaxis()->SetTitle("1/N_{Events} dN/d#rho");
    fh020RhoNoLeading->Sumw2();
    
    fh020Rho1B = new TH1D("fh020Rho1B","0-20% Background Density #rho_{n}",RhoPtBins,RhoPtLow,RhoPtUp);
    fh020Rho1B->GetXaxis()->SetTitle("p_{T}/Area (GeV/c)");
    fh020Rho1B->GetYaxis()->SetTitle("1/N_{Events} dN/d#rho");
    fh020Rho1B->Sumw2();
    
    fh020Rho2B = new TH1D("fh020Rho2B","0-20% Background Density #rho_{dijet}",RhoPtBins,RhoPtLow,RhoPtUp);
    fh020Rho2B->GetXaxis()->SetTitle("p_{T}/Area (GeV/c)");
    fh020Rho2B->GetYaxis()->SetTitle("1/N_{Events} dN/d#rho");
    fh020Rho2B->Sumw2();

    fh020Rho2BCore = new TH1D("fh020Rho2BCore","0-20% Background Density #rho_{dijet}",RhoPtBins,RhoPtLow,RhoPtUp);
    fh020Rho2BCore->GetXaxis()->SetTitle("p_{T}/Area (GeV/c)");
    fh020Rho2BCore->GetYaxis()->SetTitle("1/N_{Events} dN/d#rho");
    fh020Rho2BCore->Sumw2();

    fh020Rho3 = new TH1D("fh020Rho3","0-20% Charged Background Density #rho_{char}",RhoPtBins,RhoPtLow,RhoPtUp);
    fh020Rho3->GetXaxis()->SetTitle("p_{T}/Area (GeV/c)");
    fh020Rho3->GetYaxis()->SetTitle("1/N_{Events} dN/d#rho");
    fh020Rho3->Sumw2();

    fh020Rho3NoJets = new TH1D("fh020Rho3NoJets","0-20% Charged Background Density #rho_{char} for Events with No Signal Jets",RhoPtBins,RhoPtLow,RhoPtUp);
    fh020Rho3NoJets->GetXaxis()->SetTitle("p_{T}/Area (GeV/c)");
    fh020Rho3NoJets->GetYaxis()->SetTitle("1/N_{Events} dN/d#rho");
    fh020Rho3NoJets->Sumw2();

    fh020Rho3DiJets = new TH1D("fh020Rho3DiJets","0-20% Charged Background Density #rho_{char} for Events with a DiJet",RhoPtBins,RhoPtLow,RhoPtUp);
    fh020Rho3DiJets->GetXaxis()->SetTitle("p_{T}/Area (GeV/c)");
    fh020Rho3DiJets->GetYaxis()->SetTitle("1/N_{Events} dN/d#rho");
    fh020Rho3DiJets->Sumw2();

    fh020Rho3Perp = new TH1D("fh020Rho3Perp","0-20% Charged Background Density #rho_{char} for Events with a DiJet at Median Angle Between Dijets",RhoPtBins,RhoPtLow,RhoPtUp);
    fh020Rho3Perp->GetXaxis()->SetTitle("p_{T}/Area (GeV/c)");
    fh020Rho3Perp->GetYaxis()->SetTitle("1/N_{Events} dN/d#rho");
    fh020Rho3Perp->Sumw2();

    fh020RhokT = new TH1D("fh020RhokT","0-20% Background Density #rho_{k_{T}}",RhoPtBins,RhoPtLow,RhoPtUp);
    fh020RhokT->GetXaxis()->SetTitle("p_{T}/Area (GeV/c)");
    fh020RhokT->GetYaxis()->SetTitle("1/N_{Events} dN/d#rho");
    fh020RhokT->Sumw2();

    // Delta pT Plots
    Int_t DeltaPtBins=150;
    Double_t DeltaPtLow=-50.0;
    Double_t DeltaPtUp=100.0;

    fhDeltaPtTotal = new TH1D("fhDeltaPtTotal","#deltap_{T} distribution for Random Cones using total #rho",DeltaPtBins,DeltaPtLow,DeltaPtUp);
    fhDeltaPtTotal->GetXaxis()->SetTitle("#deltap_{T} (GeV/c)");
    fhDeltaPtTotal->GetYaxis()->SetTitle("1/N_{Events} dN/dp_{T}");
    fhDeltaPtTotal->Sumw2();

    fhDeltaPtNoLeading = new TH1D("fhDeltaPtNoLeading","#deltap_{T} distribution for Random Cones using leading jet #rho",DeltaPtBins,DeltaPtLow,DeltaPtUp);
    fhDeltaPtNoLeading->GetXaxis()->SetTitle("#deltap_{T} (GeV/c)");
    fhDeltaPtNoLeading->GetYaxis()->SetTitle("1/N_{Events} dN/dp_{T}");
    fhDeltaPtNoLeading->Sumw2();

    fhDeltaPt1B = new TH1D("fhDeltaPt1B","#deltap_{T} distribution for Random Cones using method 1B #rho",DeltaPtBins,DeltaPtLow,DeltaPtUp);
    fhDeltaPt1B->GetXaxis()->SetTitle("#deltap_{T} (GeV/c)");
    fhDeltaPt1B->GetYaxis()->SetTitle("1/N_{Events} dN/dp_{T}");
    fhDeltaPt1B->Sumw2();

    fhDeltaPt2B = new TH1D("fhDeltaPt2B","#deltap_{T} distribution for Random Cones using method 2B #rho",DeltaPtBins,DeltaPtLow,DeltaPtUp);
    fhDeltaPt2B->GetXaxis()->SetTitle("#deltap_{T} (GeV/c)");
    fhDeltaPt2B->GetYaxis()->SetTitle("1/N_{Events} dN/dp_{T}");
    fhDeltaPt2B->Sumw2();

    fhDeltaPtkT = new TH1D("fhDeltaPtkT","#deltap_{T} distribution for Random Cones using kT #rho",DeltaPtBins,DeltaPtLow,DeltaPtUp);
    fhDeltaPtkT->GetXaxis()->SetTitle("#deltap_{T} (GeV/c)");
    fhDeltaPtkT->GetYaxis()->SetTitle("1/N_{Events} dN/dp_{T}");
    fhDeltaPtkT->Sumw2();

    fhDeltaRho01 = new TH1D("fhDeltaRho01","Event by Event Differential between #rho_{0} and #rho_{1}",RhoPtBins,RhoPtLow,RhoPtUp);
    fhDeltaRho01->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fhDeltaRho01->GetYaxis()->SetTitle("1/N_{Events} dN/dp_{T}");
    fhDeltaRho01->Sumw2();

    fhDeltaRho0DiJet = new TH1D("fhDeltaRho0DiJet","Event by Event Differential between #rho_{0} and #rho_{dijet}",DeltaRhoPtBins,DeltaRhoPtLow,DeltaRhoPtUp);
    fhDeltaRho0DiJet->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fhDeltaRho0DiJet->GetYaxis()->SetTitle("1/N_{Events} dN/dp_{T}");
    fhDeltaRho0DiJet->Sumw2();

    // DiJet Plots
    Int_t DiJetBins=20;
    Double_t DiJetLow=0.0;
    Double_t DiJetUp=1.0;
    
    fh020DiJetAsy = new TH1D("fh020DiJetAsy","0-20% Centrality, Di-Jet Asymmetry",DiJetBins,DiJetLow,DiJetUp);
    fh020DiJetAsy->GetXaxis()->SetTitle("A_{j}");
    fh020DiJetAsy->GetYaxis()->SetTitle("1/N_{Events} dN/dA_{j}");
    fh020DiJetAsy->Sumw2();

    fhDiJetCenAsy = new TH2D("fhDiJetCenAsy","Di-Jet Asymmetry",DiJetBins,DiJetLow,DiJetUp,fCentralityBins*CentralityBinMult,fCentralityLow,fCentralityUp);
    fhDiJetCenAsy->GetXaxis()->SetTitle("A_{j}");
    fhDiJetCenAsy->GetYaxis()->SetTitle("Centrality (V0M)");
    fhDiJetCenAsy->GetZaxis()->SetTitle("1/N_{Events} dN/dA_{j}");
    fhDiJetCenAsy->Sumw2();

    fh020DiJetDeltaPhi = new TH1D("fh020DiJetDeltaPhi","0-20% Centrality, Di-Jet #Delta#phi",TCBins,fTPCPhiMin,0.5*fTPCPhiMax);
    fh020DiJetDeltaPhi->GetXaxis()->SetTitle("#Delta#phi");
    fh020DiJetDeltaPhi->GetYaxis()->SetTitle("1/N_{Events} dN/d#Delta#phi");
    fh020DiJetDeltaPhi->Sumw2();

    fhDiJetCenDeltaPhi = new TH2D("fhDiJetCenDeltaPhi","Di-Jet #Delta#phi",TCBins,fTPCPhiMin,0.5*fTPCPhiMax,fCentralityBins*CentralityBinMult,fCentralityLow,fCentralityUp);
    fhDiJetCenDeltaPhi->GetXaxis()->SetTitle("#Delta#phi");
    fhDiJetCenDeltaPhi->GetYaxis()->SetTitle("Centrality (V0M)");
    fhDiJetCenDeltaPhi->GetZaxis()->SetTitle("1/N_{Events} dN/d#Delta#phi");
    fhDiJetCenDeltaPhi->Sumw2();

    fhDiJetEMCalLeadingPt = new TH1D("fhDiJetEMCalLeadingPt","Leading EMCal anti-k_{T} cluster p_{T} in a dijet event",EMCalClusterBins,EMCalClusterLow,EMCalClusterUp);
    fhDiJetEMCalLeadingPt->GetXaxis()->SetTitle("p_{T}");
    fhDiJetEMCalLeadingPt->GetYaxis()->SetTitle("1/N_{Events} dN/dp_{T}");
    fhDiJetEMCalLeadingPt->Sumw2();

    fhDiJetEMCalLeadingDeltaPhi = new TH1D("fhDiJetEMCalLeadingDeltaPhi","Leading EMCal anti-k_{T} cluster #Delta#phi in a dijet event",TCBins,fTPCPhiMin,0.5*fTPCPhiMax);
    fhDiJetEMCalLeadingDeltaPhi->GetXaxis()->SetTitle("#Delta#phi");
    fhDiJetEMCalLeadingDeltaPhi->GetYaxis()->SetTitle("1/N_{Events} dN/d#Delta#phi");
    fhDiJetEMCalLeadingDeltaPhi->Sumw2();
    
    // Profiles
    // Background Density vs Centrality (fpRho...)
    fpRhoTotal = new TProfile("fpRhoTotal","Background Profile vs Centrality",CentralityBinMult*fCentralityBins,fCentralityLow,fCentralityUp);
    fpRhoTotal->GetXaxis()->SetTitle("Centrality (V0M)");
    fpRhoTotal->GetYaxis()->SetTitle("p_{T}/Area (GeV/c)");

    fpRhoNoLeading = new TProfile("fpRhoNoLeading","Background Profile vs Centrality",CentralityBinMult*fCentralityBins,fCentralityLow,fCentralityUp);
    fpRhoNoLeading->GetXaxis()->SetTitle("Centrality (V0M)");
    fpRhoNoLeading->GetYaxis()->SetTitle("p_{T}/Area (GeV/c)");

    fpRho1B = new TProfile("fpRho1B","Background Profile vs Centrality",CentralityBinMult*fCentralityBins,fCentralityLow,fCentralityUp);
    fpRho1B->GetXaxis()->SetTitle("Centrality (V0M)");
    fpRho1B->GetYaxis()->SetTitle("p_{T}/Area (GeV/c)");

    fpRho2B = new TProfile("fpRho2B","Background Profile vs Centrality",CentralityBinMult*fCentralityBins,fCentralityLow,fCentralityUp);
    fpRho2B->GetXaxis()->SetTitle("Centrality (V0M)");
    fpRho2B->GetYaxis()->SetTitle("p_{T}/Area (GeV/c)");

    fpRho2BCore = new TProfile("fpRho2BCore","Background Profile vs Centrality",CentralityBinMult*fCentralityBins,fCentralityLow,fCentralityUp);
    fpRho2BCore->GetXaxis()->SetTitle("Centrality (V0M)");
    fpRho2BCore->GetYaxis()->SetTitle("p_{T}/Area (GeV/c)");

    fpRho3 = new TProfile("fpRho3","Background Profile vs Centrality",CentralityBinMult*fCentralityBins,fCentralityLow,fCentralityUp);
    fpRho3->GetXaxis()->SetTitle("Centrality (V0M)");
    fpRho3->GetYaxis()->SetTitle("p_{T}/Area (GeV/c)");

    fpRho3NoJets = new TProfile("fpRho3NoJets","Background Profile vs Centrality",CentralityBinMult*fCentralityBins,fCentralityLow,fCentralityUp);
    fpRho3NoJets->GetXaxis()->SetTitle("Centrality (V0M)");
    fpRho3NoJets->GetYaxis()->SetTitle("p_{T}/Area (GeV/c)");

    fpRho3DiJets = new TProfile("fpRho3DiJets","Background Profile vs Centrality",CentralityBinMult*fCentralityBins,fCentralityLow,fCentralityUp);
    fpRho3DiJets->GetXaxis()->SetTitle("Centrality (V0M)");
    fpRho3DiJets->GetYaxis()->SetTitle("p_{T}/Area (GeV/c)");

    fpRho3Perp = new TProfile("fpRho3Perp","Background Profile vs Centrality",CentralityBinMult*fCentralityBins,fCentralityLow,fCentralityUp);
    fpRho3Perp->GetXaxis()->SetTitle("Centrality (V0M)");
    fpRho3Perp->GetYaxis()->SetTitle("p_{T}/Area (GeV/c)");

    fpRhoScale = new TProfile("fpRhoScale","Scaling Factor Profile vs Centrality",CentralityBinMult*fCentralityBins,fCentralityLow,fCentralityUp);
    fpRhoScale->GetXaxis()->SetTitle("Centrality (V0M)");
    fpRhoScale->GetYaxis()->SetTitle("Scale Factor");
    
    fpRhokT = new TProfile("fpRhokT","Background Profile vs Centrality",CentralityBinMult*fCentralityBins,fCentralityLow,fCentralityUp);
    fpRhokT->GetXaxis()->SetTitle("Centrality (V0M)");
    fpRhokT->GetYaxis()->SetTitle("p_{T}/Area (GeV/c)");

    fpRhoChargedkT = new TProfile("fpRhoChargedkT","Background Profile vs Centrality",CentralityBinMult*fCentralityBins,fCentralityLow,fCentralityUp);
    fpRhoChargedkT->GetXaxis()->SetTitle("Centrality (V0M)");
    fpRhoChargedkT->GetYaxis()->SetTitle("p_{T}/Area (GeV/c)");

    fpRhoScalekT = new TProfile("fpRhoScalekT","Scaling Factor Profile vs Centrality",CentralityBinMult*fCentralityBins,fCentralityLow,fCentralityUp);
    fpRhoScalekT->GetXaxis()->SetTitle("Centrality (V0M)");
    fpRhoScalekT->GetYaxis()->SetTitle("Scale Factor");

    // Others
    fpEventMult = new TProfile("fpEventMult","Event Multiplcity vs Centrality",CentralityBinMult*fCentralityBins,fCentralityLow,fCentralityUp);
    fpEventMult->GetXaxis()->SetTitle("Centrality (V0M)");
    fpEventMult->GetYaxis()->SetTitle("Multiplicity");

    fpJetPtRhoTotal = new TProfile("fpJetPtRhoTotal","#rho_{0} Profile vs Leading jet p_{T}",JetPtBins,JetPtLow,JetPtUp);
    fpJetPtRhoTotal->GetXaxis()->SetTitle("Leading jet p_{T} (GeV/c)");
    fpJetPtRhoTotal->GetYaxis()->SetTitle("p_{T}/Area (GeV/c)");

    fpJetPtRhoNoLeading = new TProfile("fpJetPtRhoNoLeading","#rho_{1} Profile vs Leading jet p_{T}",JetPtBins,JetPtLow,JetPtUp);
    fpJetPtRhoNoLeading->GetXaxis()->SetTitle("Leading jet p_{T} (GeV/c)");
    fpJetPtRhoNoLeading->GetYaxis()->SetTitle("p_{T}/Area (GeV/c)");
    
    fpJetPtRhokT = new TProfile("fpJetPtRhokT","#rho_{k_{T}} Profile vs Leading jet p_{T}",JetPtBins,JetPtLow,JetPtUp);
    fpJetPtRhokT->GetXaxis()->SetTitle("Leading jet p_{T} (GeV/c)");
    fpJetPtRhokT->GetYaxis()->SetTitle("p_{T}/Area (GeV/c)");

    // QA::2D Energy Density Profiles for Tracks and CLusters
    fpTrackPtProfile = new TProfile2D("fpTrackPtProfile","2D Profile of track pT density throughout the TPC",TCBins,fTPCEtaMin,fTPCEtaMax,TCBins,fTPCPhiMin,fTPCPhiMax);
    fpTrackPtProfile->GetXaxis()->SetTitle("#eta");
    fpTrackPtProfile->GetYaxis()->SetTitle("#phi");
    fpTrackPtProfile->GetZaxis()->SetTitle("p_{T} density (GeV/Area)");
    
    fpClusterPtProfile = new TProfile2D("fpClusterPtProfile","2D Profile of cluster pT density throughout the EMCal",TCBins,fEMCalEtaMin,fEMCalEtaMax,TCBins,fEMCalPhiMin,fEMCalPhiMax);
    fpClusterPtProfile->GetXaxis()->SetTitle("#eta");
    fpClusterPtProfile->GetYaxis()->SetTitle("#phi");
    fpClusterPtProfile->GetZaxis()->SetTitle("p_{T} density (GeV/Area)");

    TString temp_name="";
    TString title_name="";
    
    fEtaProfileBins=14;
    fEtaProfileLow=-0.7;
    fEtaProfileUp=0.7;
    
    for (Int_t i=0;i<14;i++)
    {
        temp_name=Form("fpJetEtaProfile%d",i);
        title_name=Form("Jet Energy Density #eta Profile for ALL p_{T}, 0-20 Centrality, and eta=%g to %g",(i-7)/10.,(i-6)/10.);
        
        fpJetEtaProfile[i]= new TProfile(temp_name,title_name,fEtaProfileBins,fEtaProfileLow,fEtaProfileUp);
        fpJetEtaProfile[i]->GetXaxis()->SetTitle("#eta");
        fpJetEtaProfile[i]->GetYaxis()->SetTitle("p_{T}/Area (GeV/c)");
        fOutput->Add(fpJetEtaProfile[i]);
        temp_name="";
        title_name="";

        temp_name=Form("fpJetAbsEtaProfile%d",i);
        title_name=Form("Jet Energy Density #Delta #eta Profile for ALL p_{T}, 0-20 Centrality, and eta=%g to %g",(i-7)/10.,(i-6)/10.);
        
        fpJetAbsEtaProfile[i]= new TProfile(temp_name,title_name,fEtaProfileBins,0,2*fEtaProfileUp);
        fpJetAbsEtaProfile[i]->GetXaxis()->SetTitle("#Delta#eta");
        fpJetAbsEtaProfile[i]->GetYaxis()->SetTitle("p_{T}/Area (GeV/c)");
        fOutput->Add(fpJetAbsEtaProfile[i]);
        temp_name="";
        title_name="";
}
    
    fEDProfileRBins=50;
    fEDProfileRLow=0.0;
    fEDProfileRUp=0.5;
    fEDProfilePtBins=100;
    fEDProfilePtLow=0.0;
    fEDProfilePtUp=100.0;
    fEDProfileEtaBins=4;
    fEDProfileEtaLow=-0.2;
    fEDProfileEtaUp=0.2;
    
    for (Int_t i=0;i<8;i++)
    {
        temp_name=Form("fpChargedJetRProfile%d",i);
        title_name=Form("Charged Jet Energy Density Radius Profile for ALL p_{T}, 0-20 Centrality, and eta=%g to %g",(i-4)/10.,(i-3)/10.);
        
        fpChargedJetRProfile[i]= new TProfile(temp_name,title_name,fEDProfileRBins,fEDProfileRLow,fEDProfileRUp);
        fpChargedJetRProfile[i]->GetXaxis()->SetTitle("Radius");
        fpChargedJetRProfile[i]->GetYaxis()->SetTitle("p_{T}/Area (GeV/c)");
        fOutput->Add(fpChargedJetRProfile[i]);
        temp_name="";
        title_name="";
    }
    
    for (Int_t i=0;i<4;i++)
    {
        temp_name=Form("fpJetRProfile%d",i);
        title_name=Form("Jet Energy Density Radius Profile for ALL p_{T}, 0-20 Centrality, and eta=%g to %g",(i-2)/10.,(i-1)/10.);
        
        fpJetRProfile[i]= new TProfile(temp_name,title_name,fEDProfileRBins,fEDProfileRLow,fEDProfileRUp);
        fpJetRProfile[i]->GetXaxis()->SetTitle("Radius");
        fpJetRProfile[i]->GetYaxis()->SetTitle("p_{T}/Area (GeV/c)");
        fOutput->Add(fpJetRProfile[i]);
        temp_name="";
        title_name="";
    }

    for (Int_t i=0;i<10;i++)
    {
        temp_name=Form("fpChargedJetEDProfile%d0",i);
        title_name=Form("Charged Jet Energy Density Profile for %d0-%d0 Centrality",i,i+1);
        
        fpChargedJetEDProfile[i]= new TProfile3D(temp_name,title_name,fEDProfilePtBins,fEDProfilePtLow,fEDProfilePtUp,fEDProfileEtaBins+4,fEDProfileEtaLow-0.2,fEDProfileEtaUp+0.2,fEDProfileRBins,fEDProfileRLow,fEDProfileRUp);
        fpChargedJetEDProfile[i]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
        fpChargedJetEDProfile[i]->GetYaxis()->SetTitle("Pseudorapidity");
        fpChargedJetEDProfile[i]->GetZaxis()->SetTitle("Radius");
        fOutput->Add(fpChargedJetEDProfile[i]);
        temp_name="";
        title_name="";

        temp_name=Form("fpJetEDProfile%d0",i);
        title_name=Form("Jet Energy Density Profile for %d0-%d0 Centrality",i,i+1);
        
        fpJetEDProfile[i]= new TProfile3D(temp_name,title_name,fEDProfilePtBins,fEDProfilePtLow,fEDProfilePtUp,fEDProfileEtaBins,fEDProfileEtaLow,fEDProfileEtaUp,fEDProfileRBins,fEDProfileRLow,fEDProfileRUp);
        fpJetEDProfile[i]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
        fpJetEDProfile[i]->GetYaxis()->SetTitle("Pseudorapidity");
        fpJetEDProfile[i]->GetZaxis()->SetTitle("Radius");
        fOutput->Add(fpJetEDProfile[i]);
        temp_name="";
        title_name="";
    }
    
    fOutput->Add(fhTrackPt);
    fOutput->Add(fhTrackEta);
    fOutput->Add(fhTrackPhi);
    fOutput->Add(fhTrackEtaPhi);
    fOutput->Add(fhClusterPt);
    fOutput->Add(fhClusterEta);
    fOutput->Add(fhClusterPhi);
    fOutput->Add(fhClusterEtaPhi);
    fOutput->Add(fhCentrality);
    fOutput->Add(fhBckgMult);
    fOutput->Add(fhBckgFluc);
    fOutput->Add(fhChargedJetPt);
    fOutput->Add(fhChargedJetPtAreaCut);
    fOutput->Add(fhJetPtEMCal);
    fOutput->Add(fhJetPtEMCalAreaCut);
    fOutput->Add(fhJetPtEMCalAreaCutSignal);
    fOutput->Add(fhJetPtTPC);
    fOutput->Add(fhJetPtTPCAreaCut);
    fOutput->Add(fhJetPtArea);
    fOutput->Add(fhRhoTotal);
    fOutput->Add(fhRhoNoLeading);
    fOutput->Add(fhJetTPtRhoTotal);
    fOutput->Add(fhJetTPtRhoTotalSignal);
    fOutput->Add(fhJetTPtRhoNoLeading);
    fOutput->Add(fhJetTPtRhoNoLeadingSignal);
    fOutput->Add(fhJetTrigR1A);
    fOutput->Add(fhJetTPt1B);
    fOutput->Add(fhJetTPt1BSignal);
    fOutput->Add(fhEMCalBckg1B);
    fOutput->Add(fhRho1B);
    fOutput->Add(fhJetTPt1C);
    fOutput->Add(fhRho1C);
    fOutput->Add(fhEMCalBckg1C);
    fOutput->Add(fhEMCalJet2A);
    fOutput->Add(fh020EMCalJet2A);
    fOutput->Add(fhJetTPt2B);
    fOutput->Add(fhEMCalBckg2B);
    fOutput->Add(fhRho2B);
    fOutput->Add(fhRho3);
    fOutput->Add(fhJetTPt3);
    fOutput->Add(fhDeltaPtTotal);
    fOutput->Add(fhDeltaPtNoLeading);
    fOutput->Add(fhDeltaPt1B);
    fOutput->Add(fhJetConstituentPt);
    fOutput->Add(fhDeltaRho01);
    fOutput->Add(fhEMCalCellCounts);
    fOutput->Add(fh020RhoTotal);
    fOutput->Add(fh020RhoNoLeading);
    fOutput->Add(fh020Rho1B);
    fOutput->Add(fh020Rho2B);
    fOutput->Add(fh020Rho3);
    fOutput->Add(fh020JetPtEMCal);
    fOutput->Add(fh020JetPtEMCalAreaCut);
    fOutput->Add(fh020JetPtEMCalAreaCutSignal);
    fOutput->Add(fh020JetTPtRhoTotal);
    fOutput->Add(fh020JetTPtRhoTotalSignal);
    fOutput->Add(fh020JetTPtRhoNoLeading);
    fOutput->Add(fh020JetTPtRhoNoLeadingSignal);
    fOutput->Add(fh020JetTPt1B);
    fOutput->Add(fh020JetTPt1BSignal);
    fOutput->Add(fh020JetTPt1C);
    fOutput->Add(fh020JetTPt2B);
    fOutput->Add(fh020JetTPt3);
    fOutput->Add(fhDeltaPt2B);
    fOutput->Add(fhDeltaPtkT);
    fOutput->Add(fh020DiJetAsy);
    fOutput->Add(fh020RhokT);
    fOutput->Add(fh020EMCalkTClusters);
    fOutput->Add(fh020EMCalAkTJets);
    fOutput->Add(fh020DiJetDeltaPhi);
    fOutput->Add(fhDiJetEMCalLeadingPt);
    fOutput->Add(fhDiJetEMCalLeadingDeltaPhi);
    fOutput->Add(fhJetPtCenEMCal);
    fOutput->Add(fhJetPtCenEMCalAreaCut);
    fOutput->Add(fhJetPtCenEMCalAreaCutSignal);
    fOutput->Add(fhJetTPtCenRhoTotal);
    fOutput->Add(fhJetTPtCenRhoTotalSignal);
    fOutput->Add(fhJetTPtCenRhoNoLeading);
    fOutput->Add(fhJetTPtCenRhoNoLeadingSignal);
    fOutput->Add(fhJetTPtCen1B);
    fOutput->Add(fhJetTPtCen1BSignal);
    fOutput->Add(fhJetTPtCen1C);
    fOutput->Add(fhJetTPtCen2B);
    fOutput->Add(fhJetTPtCen3);
    fOutput->Add(fhDiJetCenAsy);
    fOutput->Add(fhDiJetCenDeltaPhi);
    fOutput->Add(fhEMCalCenJet2A);
    fOutput->Add(fhDeltaRho0DiJet);
    fOutput->Add(fh020Rho2BCore);
    fOutput->Add(fhRho2BCore);
    fOutput->Add(fh020Rho3NoJets);
    fOutput->Add(fh020Rho3DiJets);
    fOutput->Add(fhRho3NoJets);
    fOutput->Add(fhRho3DiJets);
    fOutput->Add(fh020Rho3Perp);
    fOutput->Add(fhRho3Perp);
    
    fOutput->Add(fpEventMult);
    fOutput->Add(fpRhoTotal);
    fOutput->Add(fpRhoNoLeading);
    fOutput->Add(fpRho1B);
    fOutput->Add(fpRho2B);
    fOutput->Add(fpRho2BCore);
    fOutput->Add(fpRho3);
    fOutput->Add(fpRho3NoJets);
    fOutput->Add(fpRho3DiJets);
    fOutput->Add(fpRho3Perp);
    fOutput->Add(fpRhoScale);
    fOutput->Add(fpRhokT);
    fOutput->Add(fpRhoChargedkT);
    fOutput->Add(fpRhoScalekT);
    fOutput->Add(fpJetPtRhoTotal);
    fOutput->Add(fpJetPtRhoNoLeading);
    fOutput->Add(fpJetPtRhokT);
    
    fOutput->Add(fpTrackPtProfile);
    fOutput->Add(fpClusterPtProfile);
    
    // Post data for ALL output slots >0 here,
    // To get at least an empty histogram
    // 1 is the outputnumber of a certain weg of task 1
    PostData(1, fOutput);
}

void AliAnalysisTaskFullpAJets::UserExecOnce()
{
    // Get the event tracks from PicoTracks
    TString track_name="PicoTracks";
    fmyTracks =dynamic_cast <TClonesArray*>(InputEvent()->FindListObject(track_name));
    
    // Get the event caloclusters from CaloClustersCorr
    TString cluster_name="CaloClustersCorr";
    fmyClusters =dynamic_cast <TClonesArray*>(InputEvent()->FindListObject(cluster_name));
    
    // Get charged jets
    TString jet_algorithm=Form("Jet_AKTChargedR0%d0_PicoTracks_pT0150",fRJET);
    fmyAKTChargedJets =dynamic_cast <TClonesArray*>(InputEvent()->FindListObject(jet_algorithm));

    jet_algorithm=Form("Jet_KTChargedR0%d0_PicoTracks_pT0150",fRJET);
    fmyKTChargedJets =dynamic_cast <TClonesArray*>(InputEvent()->FindListObject(jet_algorithm));

    // Get the full jets
    jet_algorithm=Form("Jet_AKTFullR0%d0_PicoTracks_pT0150_CaloClustersCorr_ET0300",fRJET);
    fmyAKTFullJets =dynamic_cast <TClonesArray*>(InputEvent()->FindListObject(jet_algorithm));
    
    jet_algorithm=Form("Jet_KTFullR0%d0_PicoTracks_pT0150_CaloClustersCorr_ET0300",fRJET);
    fmyKTFullJets =dynamic_cast <TClonesArray*>(InputEvent()->FindListObject(jet_algorithm));
    
    jet_algorithm="";
    fIsInitialized=kTRUE;
}
//________________________________________________________________________
void AliAnalysisTaskFullpAJets::UserExec(Option_t *) 
{
    if (fIsInitialized==kFALSE)
    {
        UserExecOnce();
    }

    // Get pointer to reconstructed event
    AliVEvent *event = InputEvent();
    if (!event)
    {
        AliError("Pointer == 0, this can not happen!");
        return;
    }

    AliESDEvent* esd = dynamic_cast<AliESDEvent*>(event);
    AliAODEvent* aod = dynamic_cast<AliAODEvent*>(event);
    
    if (esd)
    {
        fEventCentrality=esd->GetCentrality()->GetCentralityPercentile("V0M");
        
        if (esd->GetPrimaryVertex()->GetNContributors()<1 || (TMath::Abs(esd->GetPrimaryVertex()->GetZ())>fVertexWindow))
        {
            return;
        }
        if (TMath::Sqrt(TMath::Power(esd->GetPrimaryVertex()->GetX(),2)+TMath::Power(esd->GetPrimaryVertex()->GetY(),2))>fVertexMaxR)
        {
            return;
        }

        esd->GetPrimaryVertex()->GetXYZ(fvertex);
    }
    else if (aod)
    {
        //cout<<"Hello AOD"<<endl;
        
        fEventCentrality=aod->GetCentrality()->GetCentralityPercentile("V0M");
        
        if (aod->GetPrimaryVertex()->GetNContributors()<1 || (TMath::Abs(aod->GetPrimaryVertex()->GetZ())>fVertexWindow))
        {
            return;
        }
        if (TMath::Sqrt(TMath::Power(aod->GetPrimaryVertex()->GetX(),2)+TMath::Power(aod->GetPrimaryVertex()->GetY(),2))>fVertexMaxR)
        {
            return;
        }

        aod->GetPrimaryVertex()->GetXYZ(fvertex);
    }
    else
    {
        AliError("Cannot get AOD/ESD event. Rejecting Event");
        return;
    }

    // Make sure Centrality isn't exactly 100% (to avoid bin filling errors in profile plots. Make it 99.99
    if (fEventCentrality>99.99)
    {
        fEventCentrality=99.99;
    }
    fhCentrality->Fill(fEventCentrality);
    
    fnTracks = fmyTracks->GetEntries();
    //cout<<"n Tracks:"<<fnTracks<<endl;
    // Reject any event that doesn't have any tracks, i.e. TPC is off
    if (fnTracks<1)
    {
        AliWarning("No PicoTracks, Rejecting Event");
        return;
    }

    fnClusters = fmyClusters->GetEntries();
    //cout<<"n Cluster:"<<fnClusters<<endl;
    if (fnClusters<1)
    {
        AliInfo("No Corrected CaloClusters, using only charged jets");
       
        TrackHisto();
        InitChargedJets();
        Method3(kFALSE);
        JetPtChargedProfile();
        DeleteArrays(kFALSE);
        
        fnEventsCharged++;
        return;
    }
    
    TrackHisto();
    ClusterHisto();
    
    // Prep the jets
    InitChargedJets();
    InitFullJets();
    EventHistos();
    
    EstimateTotalBackground();
    EstimateBackgoundMinusLJet();

    Method1A();
    Method1B();
    Method1C();
    
    // Method 2
    if (IsDiJetEvent()==kTRUE)
    {
        Method2A();
        Method2B();
        Method3DiJet();
        Method3Perp();
    }
    
    Method3(kTRUE);
    
    // Compute Jet Energy Density Profile
    JetPtChargedProfile();
    JetPtFullProfile();
    JetPtEtaProfile();
    
    // Delete Dynamic Arrays
    DeleteArrays(kTRUE);
    fnEvents++;
    
    PostData(1, fOutput);
}

//________________________________________________________________________
void AliAnalysisTaskFullpAJets::Terminate(Option_t *) //specify what you want to have done
{
    // Called once at the end of the query. Done nothing here.
}

/////////////////////////////////////////////////////////////////////////////////////////
/////////////////     User Defined Sub_Routines   ///////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////

void AliAnalysisTaskFullpAJets::TrackHisto()
{
    // Fill track histograms: Phi,Eta,Pt
    Int_t i,j;
    Int_t TCBins=100;
    TH2D *hdummypT= new TH2D("hdummypT","",TCBins,fTPCEtaMin,fTPCEtaMax,TCBins,fTPCPhiMin,fTPCPhiMax);  //!

    for (i=0;i<fnTracks;i++)
    {
        AliVTrack* vtrack =(AliVTrack*) fmyTracks->At(i);
        fhTrackPt->Fill(vtrack->Pt());
        fhTrackEta->Fill(vtrack->Eta());
        fhTrackPhi->Fill(vtrack->Phi());
        fhTrackEtaPhi->Fill(vtrack->Eta(),vtrack->Phi());
        hdummypT->Fill(vtrack->Eta(),vtrack->Phi(),vtrack->Pt());
    }
    for (i=1;i<=TCBins;i++)
    {
        for (j=1;j<=TCBins;j++)
        {
            fpTrackPtProfile->Fill(hdummypT->GetXaxis()->GetBinCenter(i),hdummypT->GetYaxis()->GetBinCenter(j),fTPCArea*TMath::Power(TCBins,-2)*hdummypT->GetBinContent(i,j));
        }
    }
    delete hdummypT;
}

void AliAnalysisTaskFullpAJets::ClusterHisto()
{
    // Fill cluster histograms: Phi,Eta,Pt
    Int_t i,j;
    Int_t TCBins=100;
    TH2D *hdummypT= new TH2D("hdummypT","",TCBins,fEMCalEtaMin,fEMCalEtaMax,TCBins,fEMCalPhiMin,fEMCalPhiMax);  //!
    AliEMCALGeometry *myAliEMCGeo = AliEMCALGeometry::GetInstance();
    Int_t myCellID=-2;
    
    for (i=0;i<fnClusters;i++)
    {
        AliVCluster* vcluster = (AliVCluster*) fmyClusters->At(i);
        TLorentzVector *cluster_vec = new TLorentzVector;
        vcluster->GetMomentum(*cluster_vec,fvertex);
        
        fhClusterPt->Fill(cluster_vec->Pt());
        fhClusterEta->Fill(cluster_vec->Eta());
        fhClusterPhi->Fill(cluster_vec->Phi());
        fhClusterEtaPhi->Fill(cluster_vec->Eta(),cluster_vec->Phi());
        hdummypT->Fill(cluster_vec->Eta(),cluster_vec->Phi(),cluster_vec->Pt());
        myAliEMCGeo->GetAbsCellIdFromEtaPhi(cluster_vec->Eta(),cluster_vec->Phi(),myCellID);
        fhEMCalCellCounts->Fill(myCellID);
        //cout<<"Cluster ID:"<<i<<"  (Eta,Phi): ("<<cluster_vec->Eta()<<","<<cluster_vec->Phi()<<")  Cell ID:"<<myCellID<<endl;
        delete cluster_vec;
    }
    for (i=1;i<=TCBins;i++)
    {
        for (j=1;j<=TCBins;j++)
        {
            fpClusterPtProfile->Fill(hdummypT->GetXaxis()->GetBinCenter(i),hdummypT->GetYaxis()->GetBinCenter(j),fEMCalArea*TMath::Power(TCBins,-2)*hdummypT->GetBinContent(i,j));
        }
    }
    //myAliEMCGeo->GetAbsCellIdFromEtaPhi(0.38,2.88,myCellID);
    //cout<<"Cell ID Test:"<<myCellID<<endl;
    delete hdummypT;
}

void AliAnalysisTaskFullpAJets::EventHistos()
{
    Int_t i,j;
    Double_t E_tracks_total=0.;
    Double_t E_caloclusters_total=0.;
    TRandom3 u(time(NULL));
    
    Double_t Eta_Center=0.5*(fEMCalEtaMin+fEMCalEtaMax);
    Double_t Phi_Center=0.5*(fEMCalPhiMin+fEMCalPhiMax);
    Int_t event_mult=0;
    Int_t clus_mult=0;
    
    TLorentzVector *dummy= new TLorentzVector;
    
    for (j=0;j<fnBckgClusters;j++)
    {
        event_mult=0;
        clus_mult=0;
        E_tracks_total=0.;
        E_caloclusters_total=0.;
        dummy->SetPtEtaPhiE(1,u.Uniform(Eta_Center-fJetR,Eta_Center+fJetR),u.Uniform(Phi_Center-fJetR,Phi_Center+fJetR),0);
        
        //  Loop over all tracks
        for (i=0;i<fnTracks;i++)
        {
            AliVTrack* vtrack =(AliVTrack*) fmyTracks->At(i);
            if (IsInEMCal(vtrack->Phi(),vtrack->Eta())==kTRUE)
            {
                event_mult++;
                TLorentzVector *track_vec = new TLorentzVector;
                track_vec->SetPtEtaPhiE(vtrack->Pt(),vtrack->Eta(),vtrack->Phi(),vtrack->E());
                if (dummy->DeltaR(*track_vec)<fJetR)
                {
                    clus_mult++;
                    E_tracks_total+=vtrack->Pt();
                }
                delete track_vec;
            }
        }
        
        //  Loop over all caloclusters
        for (i=0;i<fnClusters;i++)
        {
            AliVCluster* vcluster = (AliVCluster*) fmyClusters->At(i);
            TLorentzVector *cluster_vec = new TLorentzVector;
            vcluster->GetMomentum(*cluster_vec,fvertex);
            event_mult++;
            if (dummy->DeltaR(*cluster_vec)<fJetR)
            {
                clus_mult++;
                E_caloclusters_total+=vcluster->E();
            }
            delete cluster_vec;
        }
        //  Fill Histograms
        if (fEventCentrality<=20)
        {
            fhBckgMult->Fill(clus_mult);
            fhBckgFluc->Fill(E_tracks_total+E_caloclusters_total);
            fRCBckgFluc[j]=E_tracks_total+E_caloclusters_total;
        }
    }
    
    fpEventMult->Fill(fEventCentrality,event_mult);
    delete dummy;
}

void AliAnalysisTaskFullpAJets::InitChargedJets()
{
    // Preliminary Jet Placement and Selection Cuts
    Int_t i,j;
    Double_t kTRho=0.0;
    Double_t delta_phi=0.0;
    
    fnAKTChargedJets = fmyAKTChargedJets->GetEntries();
    fnKTChargedJets = fmyKTChargedJets->GetEntries();
    fJetPtChargedCutID = new Int_t[fnAKTChargedJets+1];
    fJetkTTPCFullID = new Int_t[fnKTChargedJets+1];
    fnJetsChargedPtCut=0;
    fnJetskTTPCFull=0;
    fPtChargedMaxID=-1;  // Initialize to not have any jet(s) fully contained within
    fPtChargedMax=0.0;
    
    fInTPCChargedFull = new Bool_t[fnAKTChargedJets+1];
    
    for (i=0;i<fnAKTChargedJets;i++)
    {
        AliEmcalJet *myJet =(AliEmcalJet*) fmyAKTChargedJets->At(i);
        
        // Determine where the jet is
        fInTPCChargedFull[i]=IsInTPC(fJetR,myJet->Phi(),myJet->Eta(),kTRUE);
        
        // Fill Histograms with appropriate Jet Kinematics
        if (fInTPCChargedFull[i]==kTRUE)
        {
            fhChargedJetPt->Fill(myJet->Pt());
            
            if (myJet->Pt()>=fPtChargedMax)
            {
                fPtChargedMax=myJet->Pt();
                fPtChargedMaxID=i;
            }
            //  Now determine if the jet is above the EMCal Jet Pt Threshold
            if (myJet->Area()>=fJetAreaThreshold)
            {
                fhChargedJetPtAreaCut->Fill(myJet->Pt());
            }
            if (myJet->Pt()>=fTPCJetThreshold && myJet->Area()>=fJetAreaThreshold)
            {
                fJetPtChargedCutID[fnJetsChargedPtCut]=i;
                fnJetsChargedPtCut++;
            }
        }
    }
    
    // Fill dijet delta phi plots
    if (fnJetsChargedPtCut>1)
    {
        AliEmcalJet *myhJet =(AliEmcalJet*) fmyAKTChargedJets->At(fPtChargedMaxID);
        j=0;
        while (j<fnJetsChargedPtCut)
        {
            if (fJetPtChargedCutID[j]==fPtChargedMaxID)
            {
                j++;
            }
            else
            {
                AliEmcalJet *myJet =(AliEmcalJet*) fmyAKTChargedJets->At(fJetPtChargedCutID[j]);
                
                delta_phi=TMath::Min(TMath::Abs(myhJet->Phi()-myJet->Phi()),2*TMath::Pi()-TMath::Abs(myhJet->Phi()-myJet->Phi()));
                fhDiJetCenDeltaPhi->Fill(delta_phi,fEventCentrality);
                if (fEventCentrality<=20)
                {
                    fh020DiJetDeltaPhi->Fill(delta_phi);
                }
                j++;
            }
        }
    }
    
    fRhokTCharged=0.0;
    // kT jets. Used for calculating rho
    Int_t nkT_mid=-1;
    for (i=0;i<fnKTChargedJets;i++)
    {
        AliEmcalJet *myJet =(AliEmcalJet*) fmyKTChargedJets->At(i);
        
        if (IsInTPC(fJetR,myJet->Phi(),myJet->Eta(),kTRUE)==kTRUE)
        {
            fJetkTTPCFullID[fnJetskTTPCFull]=i;
            fnJetskTTPCFull++;
        }
    }
   
    if (fnJetskTTPCFull>=2)
    {
        nkT_mid=fnJetskTTPCFull/2;
    }
    else if (fnJetskTTPCFull==1)
    {
        nkT_mid=0;
    }
    
    if (nkT_mid>=0)
    {
        AliEmcalJet *myJet =(AliEmcalJet*) fmyKTChargedJets->At(fJetkTTPCFullID[nkT_mid]);
        kTRho=myJet->Pt()/myJet->Area();
        fRhokTCharged=kTRho;
        fpRhoChargedkT->Fill(fEventCentrality,kTRho);
    }
}

void AliAnalysisTaskFullpAJets::InitFullJets()
{
    // Preliminary Jet Placement and Selection Cuts
    Int_t i;
    Int_t EMCalFullCount=0;
    Double_t kTRho=0.0;
    
    fnAKTFullJets = fmyAKTFullJets->GetEntries();
    fnKTFullJets = fmyKTFullJets->GetEntries();
    
    fJetPtCutID = new Int_t[fnAKTFullJets+1];
    fJetPtTPCCutID= new Int_t[fnAKTFullJets+1];
    fJetPtTotalCutID= new Int_t[fnAKTFullJets+1];
    fJetkTEMCalFullID= new Int_t[fnKTFullJets+1];
    
    fnJetsPtCut=0;
    fnJetsPtTPCCut=0;
    fnJetsPtTotalCut=0;
    fnJetskTEMCalFull=0;
    
    fPtMaxID=-1;  // Initialize to not have any jet(s) in EMCal
    fPtFullMaxID=-1;  // Initialize to not have any jet(s) fully contained within EMCal
    fPtTPCMaxID=-1;  // Initialize to not have any jet(s) in TPC
    fPtFullTPCMaxID=-1;  // Initialize to not have any jet(s) fully contained within TPC
    fPtTotalMaxID=-1;  // Initialize to not have any jet(s) in Total Acceptance
    
    fPtMax=0.;
    fPtFullMax=0.;
    fPtTPCMax=0.;
    fPtFullTPCMax=0.;
    fPtTotalMax=0.;
    
    fInEMCal = new Bool_t[fnAKTFullJets+1];
    fInEMCalFull = new Bool_t[fnAKTFullJets+1];
    fInTPCFull = new Bool_t[fnAKTFullJets+1];

    for (i=0;i<fnAKTFullJets;i++)
    {
        AliEmcalJet *myJet =(AliEmcalJet*) fmyAKTFullJets->At(i);
        
        // Area distribution of the jet
        fhJetPtArea->Fill(myJet->Pt(),myJet->Area());
        
        // Determine where the jet is
        fInEMCal[i]=IsInEMCalPart(fJetR,myJet->Phi(),myJet->Eta());
        fInEMCalFull[i]=IsInEMCalFull(fJetR,myJet->Phi(),myJet->Eta());
        fInTPCFull[i]=IsInTPCFull(fJetR,myJet->Phi(),myJet->Eta());
        
        // Fill Histograms with appropriate Jet Kinematics
        if (fInEMCal[i]==kTRUE)
        {
            // Method 1A
            if (myJet->Pt()>fPtMax)
            {
                fPtMax=myJet->Pt();
                fPtMaxID=i;
            }
            
            // Method 1B
            if (fInEMCalFull[i]==kTRUE)
            {
                // Fill Jet Pt Distribution
                fhJetPtEMCal->Fill(myJet->Pt());
                fhJetPtCenEMCal->Fill(myJet->Pt(),fEventCentrality);
                if (fEventCentrality<=20)
                {
                    fh020JetPtEMCal->Fill(myJet->Pt());
                }
                if (myJet->Area()>=fJetAreaThreshold)
                {
                    fhJetPtEMCalAreaCut->Fill(myJet->Pt());
                    fhJetPtCenEMCalAreaCut->Fill(myJet->Pt(),fEventCentrality);
                    if (fEventCentrality<=20)
                    {
                        fh020JetPtEMCalAreaCut->Fill(myJet->Pt());
                        EMCalFullCount++;
                    }
                    if (myJet->Pt()>=fEMCalJetThreshold)
                    {
                        fhJetPtEMCalAreaCutSignal->Fill(myJet->Pt());
                        fhJetPtCenEMCalAreaCutSignal->Fill(myJet->Pt(),fEventCentrality);
                        if (fEventCentrality<=20)
                        {
                            fh020JetPtEMCalAreaCutSignal->Fill(myJet->Pt());
                        }
                    }
                }

                if (myJet->Pt()>=fPtFullMax)
                {
                    fPtFullMax=myJet->Pt();
                    fPtFullMaxID=i;
                }

                //  Now determine if the jet is above the EMCal Jet Pt Threshold
                if (myJet->Pt()>=fEMCalJetThreshold)
                {
                    fJetPtCutID[fnJetsPtCut]=i;
                    fnJetsPtCut++;
                }
            }
        }
        else
        {
            // Method 2A
            if (myJet->Pt()>fPtTPCMax)
            {
                fPtTPCMax=myJet->Pt();
                fPtTPCMaxID=i;
            }
            if (fInTPCFull[i]==kTRUE)
            {
                // Jet Pt distribution
                fhJetPtTPC->Fill(myJet->Pt());
                if (myJet->Area()>=fJetAreaThreshold)
                {
                    fhJetPtTPCAreaCut->Fill(myJet->Pt());
                }
                
                if (myJet->Pt()>fPtFullTPCMax)
                {
                    fPtFullTPCMax=myJet->Pt();
                    fPtFullTPCMaxID=i;
                }
            }
            
            //  Now determine if the jet is above the TPC Jet Pt Threshold
            if (myJet->Pt()>=fTPCJetThreshold)
            {
                fJetPtTPCCutID[fnJetsPtTPCCut]=i;
                fnJetsPtTPCCut++;
            }
        }
        // Find all jet(s) above the threshold within the Detector (TPC+EMCal; No Fudicial cut) for Plan2
        if (myJet->Pt()>fPtTotalMax)
        {
            fPtTotalMax=myJet->Pt();
            fPtTotalMaxID=i;
        }
        // And if they are above the threshold?
        if (myJet->Pt()>=fTPCJetThreshold)
        {
            fJetPtTotalCutID[fnJetsPtTotalCut]=i;
            fnJetsPtTotalCut++;
        }
    }
    fh020EMCalAkTJets->Fill(EMCalFullCount);
    EMCalFullCount=0;
    
    fRhokTTotal=0.0;
    // kT jets. Used for calculating rho
    Int_t nkT_mid=-1;
    for (i=0;i<fnKTFullJets;i++)
    {
        AliEmcalJet *myJet =(AliEmcalJet*) fmyKTFullJets->At(i);
        
        if (IsInEMCalFull(fJetR,myJet->Phi(),myJet->Eta())==kTRUE)
        {
            fJetkTEMCalFullID[fnJetskTEMCalFull]=i;
            fnJetskTEMCalFull++;
        }
    }
    if (fEventCentrality<=20)
    {
        fh020EMCalkTClusters->Fill(fnJetskTEMCalFull);
    }
    
    if (fnJetskTEMCalFull>=2)
    {
        nkT_mid=fnJetskTEMCalFull/2;
    }
    else if (fnJetskTEMCalFull==1)
    {
        nkT_mid=0;
    }

    if (nkT_mid>=0)
    {
        AliEmcalJet *myJet =(AliEmcalJet*) fmyKTFullJets->At(fJetkTEMCalFullID[nkT_mid]);
        kTRho=myJet->Pt()/myJet->Area();
        fRhokTTotal=kTRho;
        fpRhokT->Fill(fEventCentrality,kTRho);
        fpJetPtRhokT->Fill(fPtFullMax,kTRho);
        if (fRhokTCharged!=0.0)
        {
            fpRhoScalekT->Fill(fEventCentrality,fRhokTTotal/fRhokTCharged);
        }
        if (fEventCentrality<=20)
        {
            FillBckgFlucDeltaPt(fhDeltaPtkT,kTRho);
            fh020RhokT->Fill(kTRho);
        }
    }
}

void AliAnalysisTaskFullpAJets::EstimateTotalBackground()
{
    Int_t i;
    Double_t E_tracks_total=0.;
    Double_t E_caloclusters_total=0.;
    Double_t EMCal_rho=0.;
    fDeltaRho01=0.0;
    
    //  Loop over all tracks
    for (i=0;i<fnTracks;i++)
    {
        AliVTrack* vtrack =(AliVTrack*) fmyTracks->At(i);
        if (IsInEMCal(vtrack->Phi(),vtrack->Eta())==kTRUE)
        {
            E_tracks_total+=vtrack->Pt();
        }
    }
    
    //  Loop over all caloclusters
    for (i=0;i<fnClusters;i++)
    {
        AliVCluster* vcluster = (AliVCluster*) fmyClusters->At(i);
        E_caloclusters_total+=vcluster->E();
    }
    
    //  Calculate the mean Background density
    EMCal_rho=(E_tracks_total+E_caloclusters_total)/fEMCalArea;
    fRhoAkTTotal=EMCal_rho;
    
    //  Fill histograms
    fhRhoTotal->Fill(EMCal_rho,fEventCentrality);
    fpRhoTotal->Fill(fEventCentrality,EMCal_rho);
    fpJetPtRhoTotal->Fill(fPtFullMax,EMCal_rho);
    FillFullCorrJetPt(fhJetTPtRhoTotal,EMCal_rho,kFALSE);
    FillFullCorrJetPt(fhJetTPtRhoTotalSignal,EMCal_rho,kTRUE);
    fDeltaRho01=EMCal_rho;
    
    FillFullCorrJetPt(fhJetTPtCenRhoTotal,EMCal_rho,kFALSE);
    FillFullCorrJetPt(fhJetTPtCenRhoTotalSignal,EMCal_rho,kTRUE);
    // Fill Background fluctuation spectrum F(A)
    if (fEventCentrality<=20)
    {
        FillFullCorrJetPt(fh020JetTPtRhoTotal,EMCal_rho,kFALSE);
        FillFullCorrJetPt(fh020JetTPtRhoTotalSignal,EMCal_rho,kTRUE);
        FillBckgFlucDeltaPt(fhDeltaPtTotal,EMCal_rho);
        fh020RhoTotal->Fill(EMCal_rho);
    }
}

void AliAnalysisTaskFullpAJets::EstimateBackgoundMinusLJet()
{
    Int_t i;
    Double_t E_tracks_total=0.;
    Double_t E_caloclusters_total=0.;
    Double_t EMCal_rho=0.;
    
    if (fPtFullMax>=fEMCalJetThreshold)
    {
        AliEmcalJet *myJet =(AliEmcalJet*) fmyAKTFullJets->At(fPtMaxID);
        TLorentzVector *temp_jet= new TLorentzVector;
        myJet->GetMom(*temp_jet);
        
        //  Loop over all tracks
        for (i=0;i<fnTracks;i++)
        {
            AliVTrack* vtrack =(AliVTrack*) fmyTracks->At(i);
            if (IsInEMCal(vtrack->Phi(),vtrack->Eta())==kTRUE)
            {
                TLorentzVector *track_vec = new TLorentzVector;
                track_vec->SetPtEtaPhiE(vtrack->Pt(),vtrack->Eta(),vtrack->Phi(),vtrack->E());
                if (temp_jet->DeltaR(*track_vec)>fJetR)
                {
                    E_tracks_total+=vtrack->Pt();
                }
                delete track_vec;
            }
        }
        
        //  Loop over all caloclusters
        for (i=0;i<fnClusters;i++)
        {
            AliVCluster* vcluster = (AliVCluster*) fmyClusters->At(i);
            TLorentzVector *cluster_vec = new TLorentzVector;
            vcluster->GetMomentum(*cluster_vec,fvertex);
            if (temp_jet->DeltaR(*cluster_vec)>fJetR)
            {
                E_caloclusters_total+=vcluster->E();
            }
            delete cluster_vec;
        }
        delete temp_jet;
        //  Calculate the mean Background density
        EMCal_rho=(E_tracks_total+E_caloclusters_total)/(fEMCalArea-TMath::Pi()*TMath::Power(fJetR,2));
    }
    else  // i.e. No signal jets-> same as total background density
    {
        //  Loop over all tracks
        for (i=0;i<fnTracks;i++)
        {
            AliVTrack* vtrack =(AliVTrack*) fmyTracks->At(i);
            if (IsInEMCal(vtrack->Phi(),vtrack->Eta())==kTRUE)
            {
                E_tracks_total+=vtrack->Pt();
            }
        }
        
        //  Loop over all caloclusters
        for (i=0;i<fnClusters;i++)
        {
            AliVCluster* vcluster = (AliVCluster*) fmyClusters->At(i);
            E_caloclusters_total+=vcluster->E();
        }
        //  Calculate the mean Background density
        EMCal_rho=(E_tracks_total+E_caloclusters_total)/fEMCalArea;
    }
    
    //  Fill histograms
    fhRhoNoLeading->Fill(EMCal_rho,fEventCentrality);
    fpRhoNoLeading->Fill(fEventCentrality,EMCal_rho);
    fpJetPtRhoNoLeading->Fill(fPtFullMax,EMCal_rho);
    FillFullCorrJetPt(fhJetTPtRhoNoLeading,EMCal_rho,kFALSE);
    FillFullCorrJetPt(fhJetTPtRhoNoLeadingSignal,EMCal_rho,kTRUE);
    fDeltaRho01-=EMCal_rho;
    FillFullDeltaRho(fhDeltaRho01,fDeltaRho01,kTRUE);
    fDeltaRho01=0.0;
    
    FillFullCorrJetPt(fhJetTPtCenRhoNoLeading,EMCal_rho,kFALSE);
    FillFullCorrJetPt(fhJetTPtCenRhoNoLeadingSignal,EMCal_rho,kTRUE);
    // Fill Background fluctuation spectrum F(A)
    if (fEventCentrality<=20)
    {
        FillFullCorrJetPt(fh020JetTPtRhoNoLeading,EMCal_rho,kFALSE);
        FillFullCorrJetPt(fh020JetTPtRhoNoLeadingSignal,EMCal_rho,kTRUE);
        FillBckgFlucDeltaPt(fhDeltaPtNoLeading,EMCal_rho);
        fh020RhoNoLeading->Fill(EMCal_rho);
    }
}

void AliAnalysisTaskFullpAJets::Method1A()
{
    Int_t i;
    Double_t delta_R;

    if (fPtMax>=fEMCalJetThreshold && fnAKTFullJets>1)
    {
        TLorentzVector *high_jet= new TLorentzVector;
        TLorentzVector *temp_jet= new TLorentzVector;
        
        AliEmcalJet *myJet =(AliEmcalJet*) fmyAKTFullJets->At(fPtMaxID);
        myJet->GetMom(*high_jet);
        //cout<<"HJ Phi="<<high_jet->Phi()<<"   HJ Eta="<<high_jet->Eta()<<endl;
        for(i=0;i<fnAKTFullJets;i++)
        {
            if (i!=fPtMaxID && fInEMCalFull[i]==kTRUE)
            {
                AliEmcalJet *myBckg =(AliEmcalJet*) fmyAKTFullJets->At(i);
                myBckg->GetMom(*temp_jet);
                delta_R=temp_jet->DeltaR(*high_jet);
                //cout<<"TJ Phi="<<temp_jet->Phi()<<"   TJ Eta="<<temp_jet->Eta()<<endl;
                //cout<<"Delta R="<<delta_R<<endl;
                if (delta_R>=(2*fJetR))
                {
                    fhJetTrigR1A->Fill(myBckg->Pt(),fPtMax,delta_R);
                }
            }
        }
        delete high_jet;
        delete temp_jet;
    }
}

void AliAnalysisTaskFullpAJets::Method1B()
{
    Int_t i,j;
    Bool_t track_away_from_jet;
    Bool_t cluster_away_from_jet;
    Double_t E_tracks_total=0.;
    Double_t E_caloclusters_total=0.;
    Double_t EMCal_rho=0.;
    Double_t jet_area_total=0.;
    
    // First, sum all tracks within the EMCal that are away from jet(s) above Pt Threshold
    fRhoTotal=0;
    for (i=0;i<fnTracks;i++)
    {
        // First, check if track is in the EMCal!!
        AliVTrack* vtrack = (AliVTrack*) fmyTracks->At(i); // pointer to reconstructed to track
        if (IsInEMCal(vtrack->Phi(),vtrack->Eta())==kTRUE)
        {
            if (fnJetsPtCut<1)
            {
                E_tracks_total+=vtrack->Pt();
            }
            else
            {
                track_away_from_jet=kTRUE;
                j=0;
                TLorentzVector *track_vec = new TLorentzVector;
                track_vec->SetPtEtaPhiE(vtrack->Pt(),vtrack->Eta(),vtrack->Phi(),vtrack->E());
                //cout<<endl<<endl<<endl<<"Event # :"<<fnEvents+1<<"  njets="<<fnAKTFullJets<<"  Threshold EMCal jets="<<fnJetsPtCut<<"  tracks # :"<<i<<","<<fnTracks<<endl;
                while (track_away_from_jet==kTRUE && j<fnJetsPtCut)
                {
                    //cout<<"j="<<j<<endl<<endl<<endl;
                    
                    TLorentzVector *jet_vec= new TLorentzVector;
                    AliEmcalJet *myJet =(AliEmcalJet*) fmyAKTFullJets->At(fJetPtCutID[j]);
                    myJet->GetMom(*jet_vec);
                    if (track_vec->DeltaR(*jet_vec)<=(fJetR))
                    {
                        track_away_from_jet=kFALSE;
                    }
                    delete jet_vec;
                    j++;
                }
                if (track_away_from_jet==kTRUE)
                {
                    E_tracks_total+=vtrack->Pt();
                }
                delete track_vec;
            }
        }
    }
    
    // Next, sum all CaloClusters within the EMCal (obviously all clusters must be within EMCal!!) that are away from jet(s) above Pt Threshold
    
    for (i=0;i<fnClusters;i++)
    {
        AliVCluster* vcluster = (AliVCluster*) fmyClusters->At(i); // pointer to reconstructed to cluster
        if (fnJetsPtCut<1)
        {
            E_caloclusters_total+=vcluster->E();
        }
        else
        {
            cluster_away_from_jet=kTRUE;
            j=0;
            
            TLorentzVector *cluster_vec = new TLorentzVector;
            vcluster->GetMomentum(*cluster_vec,fvertex);
            //cout<<endl<<endl<<endl<<"Event # :"<<fnEvents+1<<"  njets="<<fnAKTFullJets<<"  Threshold EMCal jets="<<fnJetsPtCut<<"  cluster # :"<<i<<","<<fnClusters<<endl;
            
            while (cluster_away_from_jet==kTRUE && j<fnJetsPtCut)
            {
                TLorentzVector *jet_vec= new TLorentzVector;
                AliEmcalJet *myJet =(AliEmcalJet*) fmyAKTFullJets->At(fJetPtCutID[j]);
                myJet->GetMom(*jet_vec);
                if (cluster_vec->DeltaR(*jet_vec)<=(fJetR))
                {
                    cluster_away_from_jet=kFALSE;
                }
                delete jet_vec;
                j++;
            }
            if (cluster_away_from_jet==kTRUE)
            {
                E_caloclusters_total+=vcluster->E();
            }
            delete cluster_vec;
        }
    }
    
    // Determine area of all Jets that are within the EMCal
    if (fnJetsPtCut==0)
    {
        jet_area_total=0.;
    }
    else
    {
        for (i=0;i<fnJetsPtCut;i++)
        {
            AliEmcalJet *myJet = (AliEmcalJet*) fmyAKTFullJets->At(fJetPtCutID[i]);
            jet_area_total+=AreaWithinEMCal(fJetR,myJet->Phi(),myJet->Eta());
        }
    }
    
    // Calculate Rho
    EMCal_rho=(E_tracks_total+E_caloclusters_total)/(fEMCalArea-jet_area_total);
    fRhoTotal=EMCal_rho;
    
    // Fill Background Histogram
    fhEMCalBckg1B->Fill(EMCal_rho*TMath::Pi()*TMath::Power(fJetR,2));
    fhRho1B->Fill(EMCal_rho,fEventCentrality);
    fpRho1B->Fill(fEventCentrality,EMCal_rho);
    FillFullCorrJetPt(fhJetTPt1B,EMCal_rho,kFALSE);
    FillFullCorrJetPt(fhJetTPt1BSignal,EMCal_rho,kTRUE);
    
    FillFullCorrJetPt(fhJetTPtCen1B,EMCal_rho,kFALSE);
    FillFullCorrJetPt(fhJetTPtCen1BSignal,EMCal_rho,kTRUE);
    // Fill delta pT
    if (fEventCentrality<=20)
    {
        FillFullCorrJetPt(fh020JetTPt1B,EMCal_rho,kFALSE);
        FillFullCorrJetPt(fh020JetTPt1BSignal,EMCal_rho,kTRUE);
        FillBckgFlucDeltaPt(fhDeltaPt1B,EMCal_rho);
        fh020Rho1B->Fill(EMCal_rho);
    }
}

void AliAnalysisTaskFullpAJets::Method1C()
{
    const Double_t psi_ref=0.5*(45/360.)*2*TMath::Pi(); //Input the total acceptance within the paraenthesis to be +/- psi_ref
    Int_t i;
    Double_t E_tracks_total=0.;
    Double_t E_caloclusters_total=0.;
    Double_t EMCal_rho=0.;
    Double_t psi;
    
    if (fPtFullMaxID !=-1)
    {
        AliEmcalJet *myJet = (AliEmcalJet*) fmyAKTFullJets->At(fPtFullMaxID);
        
        if (myJet->Area()>=fJetAreaThreshold)
        {
            //  Loop over all tracks
            for (i=0;i<fnTracks;i++)
            {
                AliVTrack* vtrack =(AliVTrack*) fmyTracks->At(i);
                if (IsInEMCal(vtrack->Phi(),vtrack->Eta())==kTRUE)
                {
                    if ((vtrack->Eta()>=(myJet->Eta()+fJetR)) || (vtrack->Eta()<=(myJet->Eta()-fJetR)))
                    {
                        psi=TMath::ATan((vtrack->Phi()-myJet->Phi())/(vtrack->Eta()-myJet->Eta()));
                        if ((psi>=(-1*psi_ref)) && (psi<=psi_ref))
                        {
                            //  The Track made the Cut!!
                            E_tracks_total+=vtrack->Pt();
                        }
                    }
                }
            }
            
            //  Loop over all caloclusters
            for (i=0;i<fnClusters;i++)
            {
                AliVCluster* vcluster = (AliVCluster*) fmyClusters->At(i);
                TLorentzVector *cluster_vec = new TLorentzVector;
                vcluster->GetMomentum(*cluster_vec,fvertex);
                
                if ((cluster_vec->Eta()>=(myJet->Eta()+fJetR)) || (cluster_vec->Eta()<=(myJet->Eta()-fJetR)))
                {
                    psi=TMath::ATan((cluster_vec->Phi()-myJet->Phi())/(cluster_vec->Eta()-myJet->Eta()));
                    if ((psi>=(-1*psi_ref)) && (psi<=psi_ref))
                    {
                        //  The CaloCluster made the Cut!!
                        E_caloclusters_total+=vcluster->E();
                    }
                }
            }
            
            //  Calculate the mean Background density
            EMCal_rho=(E_tracks_total+E_caloclusters_total)/TransverseArea(fJetR,psi_ref,myJet->Phi(),myJet->Eta());
            
            //  Fill histograms
            fhEMCalBckg1C->Fill(EMCal_rho*TMath::Pi()*fJetR*fJetR);
            fhRho1C->Fill(EMCal_rho,fEventCentrality);
            fhJetTPt1C->Fill(myJet->Pt()-EMCal_rho*myJet->Area());
            FillFullCorrJetPt(fhJetTPt1C,EMCal_rho,kFALSE);
            FillFullCorrJetPt(fhJetTPtCen1C,EMCal_rho,kFALSE);
            
            
            if (fEventCentrality<=20)
            {
                FillFullCorrJetPt(fh020JetTPt1C,EMCal_rho,kFALSE);
            }
        }
    }
}

void AliAnalysisTaskFullpAJets::Method2A()
{
    Int_t i;
    
    for (i=0;i<fnAKTFullJets;i++)
    {
        if (fInEMCalFull[i]==kTRUE)
        {
            AliEmcalJet *myJet =(AliEmcalJet*) fmyAKTFullJets->At(i);
            if (myJet->Area()>=fJetAreaThreshold)
            {
                fhEMCalCenJet2A->Fill(myJet->Pt(),fEventCentrality);
                fhEMCalJet2A->Fill(myJet->Pt());
                if (fEventCentrality<=20)
                {
                    fh020EMCalJet2A->Fill(myJet->Pt());
                }
            }
        }
    }
}

void AliAnalysisTaskFullpAJets::Method2B()
{
    Int_t i;
    
    Double_t E_tracks_total=0.0;
    Double_t E_caloclusters_total=0.0;
    Double_t EMCal_rho=0.0;
    Double_t E_tracks_core_total=0.0;
    Double_t E_caloclusters_core_total=0.0;
    Double_t EMCal_core_rho=0.0;
    Double_t RCore=0.4;
    
    //  Loop over all tracks
    for (i=0;i<fnTracks;i++)
    {
        AliVTrack* vtrack =(AliVTrack*) fmyTracks->At(i);
        if (IsInEMCal(vtrack->Phi(),vtrack->Eta())==kTRUE)
        {
            E_tracks_total+=vtrack->Pt();
            if (InsideRect(vtrack->Phi(),fEMCalPhiMin+RCore,fEMCalPhiMax-RCore,vtrack->Eta(),fEMCalEtaMin+RCore,fEMCalEtaMax-RCore)==kTRUE)
            {
                E_tracks_core_total+=vtrack->Pt();
            }
        }
    }
    
    //  Loop over all caloclusters
    for (i=0;i<fnClusters;i++)
    {
        AliVCluster* vcluster = (AliVCluster*) fmyClusters->At(i);
        E_caloclusters_total+=vcluster->E();
        TLorentzVector *cluster_vec = new TLorentzVector;
        vcluster->GetMomentum(*cluster_vec,fvertex);
        if (InsideRect(cluster_vec->Phi(),fEMCalPhiMin+RCore,fEMCalPhiMax-RCore,cluster_vec->Eta(),fEMCalEtaMin+RCore,fEMCalEtaMax-RCore)==kTRUE)
        {
            E_caloclusters_core_total+=vcluster->E();
        }
        delete cluster_vec;
    }
    
    //  Calculate the mean Background density
    EMCal_rho=(E_tracks_total+E_caloclusters_total)/fEMCalArea;
    EMCal_core_rho=(E_tracks_core_total+E_caloclusters_core_total)/((fEMCalPhiTotal-2*RCore)*(fEMCalEtaTotal-2*RCore));
    
    //Fill Background Cluster Histogram
    fhEMCalBckg2B->Fill(EMCal_rho*TMath::Pi()*TMath::Power(fJetR,2));
    fhRho2B->Fill(EMCal_rho,fEventCentrality);
    fpRho2B->Fill(fEventCentrality,EMCal_rho);
    FillFullCorrJetPt(fhJetTPt2B,EMCal_rho,kFALSE);
    FillFullCorrJetPt(fhJetTPtCen2B,EMCal_rho,kFALSE);
    fhDeltaRho0DiJet->Fill(fRhoAkTTotal-EMCal_rho);
    fhRho2BCore->Fill(EMCal_core_rho,fEventCentrality);
    fpRho2BCore->Fill(fEventCentrality,EMCal_core_rho);
    
    if (fEventCentrality<=20)
    {
        FillFullCorrJetPt(fh020JetTPt2B,EMCal_rho,kFALSE);
        FillBckgFlucDeltaPt(fhDeltaPt2B,EMCal_rho);
        fh020Rho2B->Fill(EMCal_rho);
        fh020Rho2BCore->Fill(EMCal_core_rho);
    }
}

void AliAnalysisTaskFullpAJets::Method3(Bool_t EMCalOn)
{
    Int_t i,j;
    Bool_t track_away_from_jet;
    Double_t E_tracks_total=0.0;
    Double_t TPC_rho=0.0;
    Double_t jet_area_total=0.0;
    
    // Calculate charged background density in events with no signal jets
    if (fnJetsChargedPtCut==0)
    {
        //  Loop over all tracks
        for (i=0;i<fnTracks;i++)
        {
            AliVTrack* vtrack =(AliVTrack*) fmyTracks->At(i);
            E_tracks_total+=vtrack->Pt();
        }
        TPC_rho=E_tracks_total/fTPCArea;
        
        fhRho3NoJets->Fill(TPC_rho,fEventCentrality);
        fpRho3NoJets->Fill(fEventCentrality,TPC_rho);
        if (fEventCentrality<=20)
        {
            fh020Rho3NoJets->Fill(TPC_rho);
        }
        E_tracks_total=0.0;
        TPC_rho=0.0;
    }
    
    // First, sum all tracks within the TPC that are away from jet(s) above Pt Threshold
    fRhoCharged=0;
    for (i=0;i<fnTracks;i++)
    {
        AliVTrack* vtrack = (AliVTrack*) fmyTracks->At(i);
        if (fnJetsChargedPtCut<1)
        {
            E_tracks_total+=vtrack->Pt();
        }
        else
        {
            track_away_from_jet=kTRUE;
            j=0;
            TLorentzVector *track_vec = new TLorentzVector;
            track_vec->SetPtEtaPhiE(vtrack->Pt(),vtrack->Eta(),vtrack->Phi(),vtrack->E());
            while (track_away_from_jet==kTRUE && j<fnJetsChargedPtCut)
            {
                TLorentzVector *jet_vec= new TLorentzVector;
                AliEmcalJet *myJet =(AliEmcalJet*) fmyAKTChargedJets->At(fJetPtChargedCutID[j]);
                myJet->GetMom(*jet_vec);
                if (track_vec->DeltaR(*jet_vec)<=fJetR)
                {
                    track_away_from_jet=kFALSE;
                }
                delete jet_vec;
                j++;
            }
            if (track_away_from_jet==kTRUE)
            {
                E_tracks_total+=vtrack->Pt();
            }
            delete track_vec;
        }
    }
    
    // Determine area of all Jets that are within the EMCal
    if (fnJetsChargedPtCut==0)
    {
        jet_area_total=0.;
    }
    else
    {
        for (i=0;i<fnJetsChargedPtCut;i++)
        {
            AliEmcalJet *myJet = (AliEmcalJet*) fmyAKTChargedJets->At(fJetPtChargedCutID[i]);
            jet_area_total+=AreaWithinTPC(fJetR,myJet->Eta());
        }
    }

    //Calculate Rho
    TPC_rho=(E_tracks_total)/(fTPCArea-jet_area_total);
    fRhoCharged=TPC_rho;
    
    //Fill Background Histogram
    fhRho3->Fill(TPC_rho,fEventCentrality);
    fpRho3->Fill(fEventCentrality,TPC_rho);
    
    //Fill "True" Jet Pt Spectrum
    for (i=0;i<fnAKTChargedJets;i++)
    {
        if (fInTPCChargedFull[i]==kTRUE)
        {
            AliEmcalJet *myJet = (AliEmcalJet*) fmyAKTChargedJets->At(i);
            if (myJet->Area()>=fJetAreaThreshold && myJet->Pt()>=fTPCJetThreshold)
            {
                fhJetTPt3->Fill(myJet->Pt()-(TPC_rho*myJet->Area()));
                fhJetTPtCen3->Fill(myJet->Pt()-(TPC_rho*myJet->Area()),fEventCentrality);
                if (fEventCentrality<=20)
                {
                    fh020JetTPt3->Fill(myJet->Pt()-(TPC_rho*myJet->Area()));
                }
            }
        }
    }
    
    // Fill Background Scaling factor profile.
    if (EMCalOn==kTRUE && fRhoCharged!=0)
    {
        fpRhoScale->Fill(fEventCentrality,(fRhoTotal/fRhoCharged));
    }
    
    if (fEventCentrality<=20)
    {
        fh020Rho3->Fill(TPC_rho);
    }
}

void AliAnalysisTaskFullpAJets::Method3DiJet()
{
    Int_t i;
    Double_t E_tracks_total=0.0;
    Double_t TPC_rho=0.0;
    Double_t jet_area_total=0.0;
    Double_t jet_pT_total=0.0;
    
    //  Loop over all tracks
    for (i=0;i<fnTracks;i++)
    {
        AliVTrack* vtrack =(AliVTrack*) fmyTracks->At(i);
        E_tracks_total+=vtrack->Pt();
    }
    
    AliEmcalJet *myhJet =(AliEmcalJet*) fmyAKTChargedJets->At(fPtChargedMaxID);
    AliEmcalJet *mybJet =(AliEmcalJet*) fmyAKTChargedJets->At(fChargedBackJetID);

    jet_area_total=myhJet->Area()+mybJet->Area();
    jet_pT_total=myhJet->Pt()+mybJet->Pt();
    TPC_rho=(E_tracks_total-jet_pT_total)/(fTPCArea-jet_area_total);
    
    fhRho3DiJets->Fill(TPC_rho,fEventCentrality);
    fpRho3DiJets->Fill(fEventCentrality,TPC_rho);
    if (fEventCentrality<=20)
    {
        fh020Rho3DiJets->Fill(TPC_rho);
    }
}

void AliAnalysisTaskFullpAJets::Method3Perp()
{
    Int_t i;
    Double_t E_tracks_total=0.0;
    Double_t TPC_rho=0.0;
    Double_t phi_perp=0.0;  // By construction, this angle must be between 90 and 270
    const Double_t delta_phi_prep=(30/360.0)*2*TMath::Pi();  // Half angle...
    
    AliEmcalJet *myhJet =(AliEmcalJet*) fmyAKTChargedJets->At(fPtChargedMaxID);
    AliEmcalJet *mybJet =(AliEmcalJet*) fmyAKTChargedJets->At(fChargedBackJetID);
    
    phi_perp=(0.5*(myhJet->Phi()+mybJet->Phi())/360.0)*2*TMath::Pi();
    
    //  Loop over all tracks
    for (i=0;i<fnTracks;i++)
    {
        AliVTrack* vtrack =(AliVTrack*) fmyTracks->At(i);
        // Select only those tracks that are within the "strip"
        if (vtrack->Phi() >= (phi_perp-delta_phi_prep) && vtrack->Phi() <= (phi_perp+delta_phi_prep))
        {
            E_tracks_total+=vtrack->Pt();
        }
    }
    
    TPC_rho=(E_tracks_total)/(fTPCEtaMax*2*delta_phi_prep);
    
    fhRho3Perp->Fill(TPC_rho,fEventCentrality);
    fpRho3Perp->Fill(fEventCentrality,TPC_rho);
    if (fEventCentrality<=20)
    {
        fh020Rho3Perp->Fill(TPC_rho);
    }
}

void AliAnalysisTaskFullpAJets::JetPtChargedProfile()
{
    Int_t i,j;
    Double_t delta_R;
    Double_t ED_pT[fEDProfileRBins];
    
    for (i=0;i<fnJetsChargedPtCut;i++)
    {
        AliEmcalJet *myJet = (AliEmcalJet*) fmyAKTChargedJets->At(fJetPtChargedCutID[i]);
        if (InsideRect(myJet->Phi(),fTPCPhiMin,fTPCPhiMax,myJet->Eta(),fTPCEtaMin+fEDProfileRUp,fTPCEtaMax-fEDProfileRUp)==kTRUE)
        {
            for (j=0;j<fEDProfileRBins;j++)
            {
                ED_pT[j]=0;
            }
            TLorentzVector *jet_vec= new TLorentzVector;
            myJet->GetMom(*jet_vec);
            // Sum all tracks in concentric rings around jet vertex
            for (j=0;j<fnTracks;j++)
            {
                AliVTrack* vtrack = (AliVTrack*) fmyTracks->At(j);
                TLorentzVector *track_vec = new TLorentzVector;
                track_vec->SetPtEtaPhiE(vtrack->Pt(),vtrack->Eta(),vtrack->Phi(),vtrack->E());
                delta_R=jet_vec->DeltaR(*track_vec);
                if (delta_R<=fEDProfileRUp)
                {
                    ED_pT[TMath::FloorNint((fEDProfileRBins/fEDProfileRUp)*delta_R)]+=vtrack->Pt();
                }
                delete track_vec;
            }
            
            //cout<<"Event Centrality:"<<fEventCentrality<<endl;
            //cout<<endl<<endl<<"Event Centrality bin:"<<TMath::FloorNint(fEventCentrality/10.)<<endl;
            for (j=0;j<fEDProfileRBins;j++)
            {
                ED_pT[j]/=TMath::Pi()*TMath::Power((fEDProfileRUp/fEDProfileRBins),2)*(2*j+1);
                //cout<<"Strip "<<j<<"  ED="<<ED_pT[j]<<endl;
                fpChargedJetEDProfile[TMath::FloorNint(fEventCentrality/10.)]->Fill(myJet->Pt(),myJet->Eta(),(fEDProfileRUp/fEDProfileRBins)*j,ED_pT[j]);
                if (fEventCentrality<=20)
                {
                    fpChargedJetRProfile[4+TMath::FloorNint(myJet->Eta()*10.)]->Fill((fEDProfileRUp/fEDProfileRBins)*j,ED_pT[j]);
                }
            }
            delete jet_vec;
        }
    }
}

void AliAnalysisTaskFullpAJets::JetPtFullProfile()
{
    Int_t i,j;
    Double_t delta_R;
    Double_t ED_pT[fEDProfileRBins];
    
    for (i=0;i<fnJetsPtCut;i++)
    {
        AliEmcalJet *myJet = (AliEmcalJet*) fmyAKTFullJets->At(fJetPtCutID[i]);
        if (InsideRect(myJet->Phi(),fEMCalPhiMin+fEDProfileRUp,fEMCalPhiMax-fEDProfileRUp,myJet->Eta(),fEMCalEtaMin+fEDProfileRUp,fEMCalEtaMax-fEDProfileRUp)==kTRUE)
        {
            for (j=0;j<fEDProfileRBins;j++)
            {
                ED_pT[j]=0;
            }
            TLorentzVector *jet_vec= new TLorentzVector;
            myJet->GetMom(*jet_vec);
            // Sum all tracks in concentric rings around jet vertex
            for (j=0;j<fnTracks;j++)
            {
                AliVTrack* vtrack = (AliVTrack*) fmyTracks->At(j);
                TLorentzVector *track_vec = new TLorentzVector;
                track_vec->SetPtEtaPhiE(vtrack->Pt(),vtrack->Eta(),vtrack->Phi(),vtrack->E());
                delta_R=jet_vec->DeltaR(*track_vec);
                if (delta_R<=fEDProfileRUp)
                {
                    ED_pT[TMath::FloorNint((fEDProfileRBins/fEDProfileRUp)*delta_R)]+=vtrack->Pt();
                }
                delete track_vec;
            }
            
            // Sum all clusters in concentric rings around jet vertex
            for (j=0;j<fnClusters;j++)
            {
                AliVCluster* vcluster = (AliVCluster*) fmyClusters->At(j);
                TLorentzVector *cluster_vec = new TLorentzVector;
                vcluster->GetMomentum(*cluster_vec,fvertex);
                delta_R=jet_vec->DeltaR(*cluster_vec);
                if (delta_R<=fEDProfileRUp)
                {
                    ED_pT[TMath::FloorNint((fEDProfileRBins/fEDProfileRUp)*delta_R)]+=vcluster->E();
                }
                delete cluster_vec;
            }
            
            for (j=0;j<fEDProfileRBins;j++)
            {
                ED_pT[j]/=TMath::Pi()*TMath::Power((fEDProfileRUp/fEDProfileRBins),2)*(2*j+1);
                fpJetEDProfile[TMath::FloorNint(fEventCentrality/10.)]->Fill(myJet->Pt(),myJet->Eta(),(fEDProfileRUp/fEDProfileRBins)*j,ED_pT[j]);
                // Fill profile if a "most" central event (0-20%)
                if (fEventCentrality<=20)
                {
                    fpJetRProfile[2+TMath::FloorNint(myJet->Eta()*10.)]->Fill((fEDProfileRUp/fEDProfileRBins)*j,ED_pT[j]);
                }
            }
            delete jet_vec;
            
            // Fill constituent histogram
            for (j=0;j<myJet->GetNumberOfTracks();j++)
            {
                AliVParticle* vparticle = (AliVParticle*) myJet->TrackAt(j,fmyTracks);
                fhJetConstituentPt->Fill(myJet->Pt(),vparticle->Pt());
            }
            
            for (j=0;j<myJet->GetNumberOfClusters();j++)
            {
                AliVCluster* vcluster = (AliVCluster*) fmyClusters->At(myJet->ClusterAt(j));
                fhJetConstituentPt->Fill(myJet->Pt(),vcluster->E());
            }
        }
    }
}

void AliAnalysisTaskFullpAJets::JetPtEtaProfile()
{
    Int_t i,j;
    Double_t eta;
    Double_t delta_eta;
    Double_t Eta_pT[fEtaProfileBins];
    Double_t Eta_abs_pT[Int_t(0.5*fEtaProfileBins)];
    
    for (i=0;i<fnJetsPtCut;i++)
    {
        AliEmcalJet *myJet = (AliEmcalJet*) fmyAKTFullJets->At(fJetPtCutID[i]);
        if (IsInEMCal(myJet->Phi(),myJet->Eta())==kTRUE)
        {
            for (j=0;j<fEtaProfileBins;j++)
            {
                Eta_pT[j]=0;
                Eta_abs_pT[j]=0;
            }

            // Sum all tracks in strips of eta away from the jet vertex
            for (j=0;j<fnTracks;j++)
            {
                AliVTrack* vtrack = (AliVTrack*) fmyTracks->At(j);
                eta=vtrack->Eta();
                delta_eta=TMath::Abs(vtrack->Eta()-myJet->Eta());
                if (IsInEMCal(vtrack->Phi(),vtrack->Eta())==kTRUE)
                {
                    Eta_pT[Int_t(0.5*fEtaProfileBins)+TMath::FloorNint(10*eta)]+=vtrack->Pt();
                    Eta_abs_pT[TMath::FloorNint(10*delta_eta)]+=vtrack->Pt();
                }
            }
            
            // Sum all clusters in strips of eta away from the jet vertex
            for (j=0;j<fnClusters;j++)
            {
                AliVCluster* vcluster = (AliVCluster*) fmyClusters->At(j);
                TLorentzVector *cluster_vec = new TLorentzVector;
                vcluster->GetMomentum(*cluster_vec,fvertex);
                eta=cluster_vec->Eta();
                delta_eta=TMath::Abs(cluster_vec->Eta()-myJet->Eta());
                Eta_pT[Int_t(0.5*fEtaProfileBins)+TMath::FloorNint(10*eta)]+=vcluster->E();
                Eta_abs_pT[TMath::FloorNint(10*delta_eta)]+=vcluster->E();
                delete cluster_vec;
            }
            
            for (j=0;j<fEtaProfileBins;j++)
            {
                Eta_pT[j]/=0.1*fEMCalPhiTotal;
                // Fill profile if a "most" central event (0-20%)
                if (j<(10*(fEMCalEtaMax-TMath::Abs(myJet->Eta()))))
                {
                    Eta_abs_pT[j]/=0.2*fEMCalPhiTotal;
                }
                else
                {
                    Eta_abs_pT[j]/=0.1*fEMCalPhiTotal;
                }
                // Fill profile if a "most" central event (0-20%)
                if (fEventCentrality<=20)
                {
                    fpJetAbsEtaProfile[7+TMath::FloorNint(myJet->Eta()*10.)]->Fill(0.1*j,Eta_abs_pT[j]);
                    fpJetEtaProfile[7+TMath::FloorNint(myJet->Eta()*10.)]->Fill(0.1*(j-7),Eta_pT[j]);
                }
            }
        }
    }
}

void AliAnalysisTaskFullpAJets::FillFullCorrJetPt(TH1D *myHisto,Double_t rho, Bool_t signal_cut)
{
    Int_t i;
    // Fill "True" Jet Pt Spectrum. Always demand that jet area is greater then area threshold.
    for (i=0;i<fnAKTFullJets;i++)
    {
        if (fInEMCalFull[i]==kTRUE)
        {
            AliEmcalJet *myJet = (AliEmcalJet*) fmyAKTFullJets->At(i);
            if (myJet->Area()>=fJetAreaThreshold)
            {
                if (signal_cut==kTRUE && myJet->Pt()>=fEMCalJetThreshold)
                {
                    myHisto->Fill(myJet->Pt()-(rho*myJet->Area()));
                }
                else if (signal_cut==kFALSE)
                {
                    myHisto->Fill(myJet->Pt()-(rho*myJet->Area()));
                }
            }
        }
    }
}

void AliAnalysisTaskFullpAJets::FillFullCorrJetPt(TH2D *myHisto,Double_t rho, Bool_t signal_cut)
{
    Int_t i;
    // Fill "True" Jet Pt Spectrum. Always demand that jet area is greater then area threshold.
    for (i=0;i<fnAKTFullJets;i++)
    {
        if (fInEMCalFull[i]==kTRUE)
        {
            AliEmcalJet *myJet = (AliEmcalJet*) fmyAKTFullJets->At(i);
            if (myJet->Area()>=fJetAreaThreshold)
            {
                if (signal_cut==kTRUE && myJet->Pt()>=fEMCalJetThreshold)
                {
                    myHisto->Fill(myJet->Pt()-(rho*myJet->Area()),fEventCentrality);
                }
                else if (signal_cut==kFALSE)
                {
                    myHisto->Fill(myJet->Pt()-(rho*myJet->Area()),fEventCentrality);
                }
            }
        }
    }
}

void AliAnalysisTaskFullpAJets::FillFullDeltaRho(TH1D *myHisto,Double_t delta_rho,Bool_t signal_cut)
{
    Int_t i;
    // Fill "True" Jet Pt Spectrum. Always demand that jet area is greater then area threshold.
    for (i=0;i<fnAKTFullJets;i++)
    {
        if (fInEMCalFull[i]==kTRUE)
        {
            AliEmcalJet *myJet = (AliEmcalJet*) fmyAKTFullJets->At(i);
            if (myJet->Area()>=fJetAreaThreshold)
            {
                if (signal_cut==kTRUE && myJet->Pt()>=fEMCalJetThreshold)
                {
                    myHisto->Fill(delta_rho);
                }
                else if (signal_cut==kFALSE)
                {
                    myHisto->Fill(delta_rho);
                }
            }
        }
    }
}

void AliAnalysisTaskFullpAJets::FillBckgFlucDeltaPt(TH1D *myHisto, Double_t rho)
{
    Int_t i;
    
    for (i=0;i<fnBckgClusters;i++)
    {
        myHisto->Fill(fRCBckgFluc[i]-rho*TMath::Pi()*fJetR*fJetR);
    }
}


void AliAnalysisTaskFullpAJets::DeleteArrays(Bool_t EMCalOn)
{
    delete [] fJetPtChargedCutID;
    delete [] fInTPCChargedFull;
    delete [] fJetkTTPCFullID;
    if (EMCalOn==kTRUE)
    {
        delete [] fJetPtCutID;
        delete [] fJetPtTPCCutID;
        delete [] fJetPtTotalCutID;
        delete [] fJetkTEMCalFullID;
        delete [] fInEMCal;
        delete [] fInEMCalFull;
        delete [] fInTPCFull;
    }
}

/////////////////////////////////////////////////////////////////////////////////////////
/////////////////     User Defined Functions      ///////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////

Bool_t AliAnalysisTaskFullpAJets::IsDiJetEvent()
{
    // Determine if event contains a di-jet within the detector. Uses charged jets.
    // Requires the delta phi of the jets to be 180 +/- 15 degrees.
    // Requires both jets to be outside of the EMCal
    // Requires both jets to be signal jets

    Int_t i;
    const Double_t dijet_delta_phi=(180/360.)*2*TMath::Pi();
    const Double_t dijet_phi_acceptance=0.5*(30/360.)*2*TMath::Pi(); //Input the total acceptance within the paraenthesis to be +/- dijet_phi_acceptance
    Double_t dummy_phi=0.0;
    Double_t dijet_asymmetry=0.0;
    Double_t delta_phi=0.0;
    fChargedBackJetID=-1;
    
    if (fnJetsChargedPtCut>1)
    {
        AliEmcalJet *myhJet =(AliEmcalJet*) fmyAKTChargedJets->At(fPtChargedMaxID);
        i=0;
        if (IsInTPCFull(fJetR,myhJet->Phi(),myhJet->Eta())==kFALSE)
        {
            return kFALSE;
        }
        while (i<fnJetsChargedPtCut)
        {
            AliEmcalJet *myJet =(AliEmcalJet*) fmyAKTChargedJets->At(fJetPtChargedCutID[i]);
            dummy_phi=TMath::Min(TMath::Abs(myhJet->Phi()-myJet->Phi()),2*TMath::Pi()-TMath::Abs(myhJet->Phi()-myJet->Phi()));
            if ((dummy_phi>(dijet_delta_phi-dijet_phi_acceptance)) && (IsInTPCFull(fJetR,myJet->Phi(),myJet->Eta())==kTRUE))
            {
                fChargedBackJetID=fJetPtChargedCutID[i];
                fnDiJetEvents++;
                dijet_asymmetry=(myhJet->Pt()-myJet->Pt())/(myhJet->Pt()+myJet->Pt());
                fhDiJetCenAsy->Fill(dijet_asymmetry,fEventCentrality);
                // Add Anti-kT Plots here...
                if (fPtFullMaxID!=-1)
                {
                    AliEmcalJet *myTestJet = (AliEmcalJet*) fmyAKTFullJets->At(fPtFullMaxID);
                    
                    fhDiJetEMCalLeadingPt->Fill(myTestJet->Pt());
                    delta_phi=TMath::Min(TMath::Abs(myhJet->Phi()-myTestJet->Phi()),2*TMath::Pi()-TMath::Abs(myhJet->Phi()-myTestJet->Phi()));
                    fhDiJetEMCalLeadingDeltaPhi->Fill(delta_phi);
                }
                if (fEventCentrality<=20)
                {
                    fh020DiJetAsy->Fill(dijet_asymmetry);
                }
                return kTRUE;
            }
            i++;
        }
    }
    return kFALSE;
}

Bool_t AliAnalysisTaskFullpAJets::InsideRect(Double_t phi,Double_t phi_min,Double_t phi_max,Double_t eta,Double_t eta_min,Double_t eta_max)
{
    if (phi>phi_min && phi<phi_max)
    {
        if (eta>eta_min && eta<eta_max)
        {
            return kTRUE;
        }
    }
    return kFALSE;
}

Bool_t AliAnalysisTaskFullpAJets::IsInEMCal(Double_t phi,Double_t eta)
{
    return InsideRect(phi,fEMCalPhiMin,fEMCalPhiMax,eta,fEMCalEtaMin,fEMCalEtaMax);
}

Bool_t AliAnalysisTaskFullpAJets::IsInEMCalFull(Double_t r,Double_t phi,Double_t eta)
{
    return InsideRect(phi,fEMCalPhiMin+r,fEMCalPhiMax-r,eta,fEMCalEtaMin+r,fEMCalEtaMax-r);
}

Bool_t AliAnalysisTaskFullpAJets::IsInEMCalPart(Double_t r,Double_t phi,Double_t eta)
{
    return InsideRect(phi,fEMCalPhiMin-r,fEMCalPhiMax+r,eta,fEMCalEtaMin-r,fEMCalEtaMax+r);
}

Bool_t AliAnalysisTaskFullpAJets::IsInTPCFull(Double_t r,Double_t phi,Double_t eta)
{
    Bool_t in_EMCal= InsideRect(phi,fEMCalPhiMin-r,fEMCalPhiMax+r,eta,fEMCalEtaMin-r,fEMCalEtaMax+r);
    Bool_t in_TPC= InsideRect(phi,fTPCPhiMin,fTPCPhiMax,eta,fTPCEtaMin+r,fTPCEtaMax-r);

    if (in_EMCal==kFALSE && in_TPC==kTRUE)
    {
        return kTRUE;
    }
    return kFALSE;
}

Bool_t AliAnalysisTaskFullpAJets::IsInTPC(Double_t r,Double_t phi,Double_t eta,Bool_t Complete)
{
    if (Complete==kTRUE)
    {
        return InsideRect(phi,fTPCPhiMin,fTPCPhiMax,eta,fTPCEtaMin+r,fTPCEtaMax-r);
    }
    return InsideRect(phi,fTPCPhiMin,fTPCPhiMax,eta,fTPCEtaMin,fTPCEtaMax);
}

Double_t AliAnalysisTaskFullpAJets::AreaWithinTPC(Double_t r,Double_t eta)
{
    Double_t z;
    if (eta<(fTPCEtaMin+r))
    {
        z=eta-fTPCEtaMin;
    }
    else if(eta>(fTPCEtaMax-r))
    {
        z=fTPCEtaMax-eta;
    }
    else
    {
        z=r;
    }
    return r*r*TMath::Pi()-AreaEdge(r,z);
}

Double_t AliAnalysisTaskFullpAJets::AreaWithinEMCal(Double_t r,Double_t phi,Double_t eta)
{
    Double_t x,y;
    
    if (phi<(fEMCalPhiMin-r) || phi>(fEMCalPhiMax+r))
    {
        x=-r;
    }
    else if (phi<(fEMCalPhiMin+r))
    {
        x=phi-fEMCalPhiMin;
    }
    else if (phi>(fEMCalPhiMin+r) && phi<(fEMCalPhiMax-r))
    {
        x=r;
    }
    else
    {
        x=fEMCalPhiMax-phi;
    }
    
    if (eta<(fEMCalEtaMin-r) || eta>(fEMCalEtaMax+r))
    {
        y=-r;
    }
    else if (eta<(fEMCalEtaMin+r))
    {
        y=eta-fEMCalEtaMin;
    }
    else if (eta>(fEMCalEtaMin+r) && eta<(fEMCalEtaMax-r))
    {
        y=r;
    }
    else
    {
        y=fEMCalEtaMax-eta;
    }

    if (x>=0 && y>=0)
    {
        if (TMath::Sqrt(x*x+y*y)>=r)
        {
            return r*r*TMath::Pi()-AreaEdge(r,x)-AreaEdge(r,y);
        }
        return r*r*TMath::Pi()-AreaEdge(r,x)-AreaEdge(r,y)+AreaOverlap(r,x,y);
    }
    else if ((x>=r && y<0) || (y>=r && x<0))
    {
        return r*r*TMath::Pi()-AreaEdge(r,x)-AreaEdge(r,y);
    }
    else if (x>0 && x<r && y<0)
    {
        Double_t a=TMath::Sqrt(r*r-x*x);
        Double_t b=TMath::Sqrt(r*r-y*y);
        if ((x-b)>0)
        {
            return r*r*TMath::ASin(b/r)+y*b;
        }
        else
        {
            return 0.5*x*a+0.5*r*r*TMath::ASin(x/r)+0.5*y*b+x*y+0.5*r*r*TMath::ASin(b/r);
        }
    }
    else if (y>0 && y<r && x<0)
    {
        Double_t a=TMath::Sqrt(r*r-x*x);
        Double_t b=TMath::Sqrt(r*r-y*y);
        if ((y-a)>0)
        {
            return r*r*TMath::ASin(a/r)+x*a;
        }
        else
        {
            return 0.5*y*b+0.5*r*r*TMath::ASin(y/r)+0.5*x*a+x*y+0.5*r*r*TMath::ASin(a/r);
        }
    }
    else
    {
        Double_t a=TMath::Sqrt(r*r-x*x);
        Double_t b=TMath::Sqrt(r*r-y*y);
        if ((x+b)<0)
        {
            return 0;
        }
        else
        {
            return 0.5*x*a+0.5*r*r*TMath::ASin(x/r)+0.5*y*b+x*y+0.5*r*r*TMath::ASin(b/r);
        }
    }
}

Double_t AliAnalysisTaskFullpAJets::AreaEdge(Double_t r,Double_t z)
{
    Double_t a=TMath::Sqrt(r*r-z*z);
    return r*r*TMath::ASin(a/r)-a*z;
}

Double_t AliAnalysisTaskFullpAJets::AreaOverlap(Double_t r,Double_t x,Double_t y)
{
    Double_t a=TMath::Sqrt(r*r-x*x);
    Double_t b=TMath::Sqrt(r*r-y*y);
    return x*y-0.5*(x*a+y*b)+0.5*r*r*(TMath::ASin(b/r)-TMath::ASin(x/r));
}

Double_t AliAnalysisTaskFullpAJets::TransverseArea(Double_t r,Double_t psi0,Double_t phi,Double_t eta)
{
    Double_t area_left,area_right;
    Double_t eta_a,eta_b,eta_up,eta_down;
    
    Double_t u=eta-fEMCalEtaMin;
    Double_t v=fEMCalEtaMax-eta;
    
    Double_t phi1=phi+u*TMath::Tan(psi0);
    Double_t phi2=phi-u*TMath::Tan(psi0);
    Double_t phi3=phi+v*TMath::Tan(psi0);
    Double_t phi4=phi-v*TMath::Tan(psi0);
    
    //Calculate the Left side area
    if (phi1>=fEMCalPhiMax)
    {
        eta_a=eta-u*((fEMCalPhiMax-phi)/(phi1-phi));
    }
    if (phi2<=fEMCalPhiMin)
    {
        eta_b=eta-u*((phi-fEMCalPhiMin)/(phi-phi2));
    }
    
    if ((phi1>=fEMCalPhiMax) && (phi2<=fEMCalPhiMin))
    {
        eta_up=TMath::Max(eta_a,eta_b);
        eta_down=TMath::Min(eta_a,eta_b);
        
        area_left=(eta_down-fEMCalEtaMin)*fEMCalPhiTotal + 0.5*(fEMCalPhiTotal+2*(eta-eta_up)*TMath::Tan(psi0))*(eta_up-eta_down) + (eta-eta_up+r)*TMath::Tan(psi0)*(eta-eta_up-r);
    }
    else if (phi1>=fEMCalPhiMax)
    {
        area_left=0.5*(fEMCalPhiMax-phi2+2*(eta-eta_a)*TMath::Tan(psi0))*(eta_a-fEMCalEtaMin) + (eta-eta_a+r)*TMath::Tan(psi0)*(eta-eta_a-r);
    }
    else if (phi2<=fEMCalPhiMin)
    {
        area_left=0.5*(phi1-fEMCalPhiMin+2*(eta-eta_b)*TMath::Tan(psi0))*(eta_b-fEMCalEtaMin) + (eta-eta_b+r)*TMath::Tan(psi0)*(eta-eta_b-r);
    }
    else
    {
        area_left=0.5*(phi1-phi2+2*r*TMath::Tan(psi0))*(u-r);
    }

    //Calculate the Right side area
    if (phi3>=fEMCalPhiMax)
    {
        eta_a=eta+v*((fEMCalPhiMax-phi)/(phi3-phi));
    }
    if (phi4<=fEMCalPhiMin)
    {
        eta_b=eta+v*((phi-fEMCalPhiMin)/(phi-phi4));
    }
    
    if ((phi3>=fEMCalPhiMax) && (phi4<=fEMCalPhiMin))
    {
        eta_up=TMath::Max(eta_a,eta_b);
        eta_down=TMath::Min(eta_a,eta_b);
        
        area_right=(fEMCalEtaMax-eta_up)*fEMCalPhiTotal + 0.5*(fEMCalPhiTotal+2*(eta_down-eta)*TMath::Tan(psi0))*(eta_down-eta_up) + (eta_down-eta+r)*TMath::Tan(psi0)*(eta_up-eta-r);
    }
    else if (phi3>=fEMCalPhiMax)
    {
        area_right=0.5*(fEMCalPhiMax-phi4+2*(eta_a-eta)*TMath::Tan(psi0))*(fEMCalEtaMax-eta_a) + (eta_a-eta+r)*TMath::Tan(psi0)*(eta_a-eta-r);
    }
    else if (phi4<=fEMCalPhiMin)
    {
        area_right=0.5*(phi3-fEMCalPhiMin+2*(eta_b-eta)*TMath::Tan(psi0))*(fEMCalEtaMax-eta_b) + (eta_b-eta+r)*TMath::Tan(psi0)*(eta_b-eta-r);
    }
    else
    {
        area_right=0.5*(phi3-phi4+2*r*TMath::Tan(psi0))*(v-r);
    }
    return area_left+area_right;
}
