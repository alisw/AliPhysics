#include "AliAnalysisTaskFullpAJets.h"

#include <Riostream.h>
#include <ctime>
#include <TString.h>
#include <TChain.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <THnSparse.h>
#include <TCanvas.h>
#include <TList.h>
#include <TLorentzVector.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TProfile3D.h>
#include <TRandom.h>
#include <TRandom3.h>
#include <TClonesArray.h>
#include <TObjArray.h>

#include "AliAnalysisTaskEmcal.h"
#include "AliAnalysisTaskEmcalJet.h"
#include "AliAnalysisManager.h"
#include "AliStack.h"
#include "AliESDtrackCuts.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliAODEvent.h"
#include "AliVEvent.h"
#include "AliMCEvent.h"
#include "AliVTrack.h"
#include "AliVCluster.h"
#include "AliEmcalJet.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALRecoUtils.h"
#include "AliVCaloCells.h"
#include "AliPicoTrack.h"
#include "AliAnalysisTaskEmcal.h"
#include "AliAnalysisTaskEmcalJet.h"
#include "Rtypes.h"

ClassImp(AliAnalysisTaskFullpAJets)

//________________________________________________________________________
AliAnalysisTaskFullpAJets::AliAnalysisTaskFullpAJets() : 
    AliAnalysisTaskEmcalJet(),

    fOutput(0),
    flTrack(0),
    flCluster(0),
    fhTrackPt(0),
    fhTrackEta(0),
    fhTrackPhi(0),
    fhGlobalTrackPt(0),
    fhGlobalTrackEta(0),
    fhGlobalTrackPhi(0),
    fhComplementaryTrackPt(0),
    fhComplementaryTrackEta(0),
    fhComplementaryTrackPhi(0),
    fhClusterPt(0),
    fhClusterEta(0),
    fhClusterPhi(0),
    fhCentrality(0),
    fhEMCalCellCounts(0),

    fhChargeAndNeutralEvents(0),
    fhChargeOnlyEvents(0),
    fhNeutralOnlyEvents(0),
    fhNothingEvents(0),
    fhEMCalChargeAndNeutralEvents(0),
    fhEMCalChargeOnlyEvents(0),
    fhEMCalNeutralOnlyEvents(0),
    fhEMCalNothingEvents(0),

    fhTrackEtaPhi(0),
    fhTrackPhiPt(0),
    fhTrackEtaPt(0),
    fhGlobalTrackEtaPhi(0),
    fhGlobalTrackPhiPt(0),
    fhGlobalTrackEtaPt(0),
    fhComplementaryTrackEtaPhi(0),
    fhComplementaryTrackPhiPt(0),
    fhComplementaryTrackEtaPt(0),
    fhClusterEtaPhi(0),
    fhClusterPhiPt(0),
    fhClusterEtaPt(0),
    fhEMCalEventMult(0),
    fhTPCEventMult(0),
    fhEMCalTrackEventMult(0),

    fpEMCalEventMult(0),
    fpTPCEventMult(0),

    fpTrackPtProfile(0),
    fpClusterPtProfile(0),

    fpFullJetEDProfile(0),
    fpChargedJetEDProfile(0),
    fpChargedJetEDProfileScaled(0),

    fTPCRawJets(0),
    fEMCalRawJets(0),
    fRhoChargedCMSScale(0),
    fRhoChargedScale(0),
    fRhoFull0(0),
    fRhoFull1(0),
    fRhoFull2(0),
    fRhoFullN(0),
    fRhoFullDijet(0),
    fRhoFullkT(0),
    fRhoFullCMS(0),
    fRhoCharged0(0),
    fRhoCharged1(0),
    fRhoCharged2(0),
    fRhoChargedN(0),
    fRhoChargedkT(0),
    fRhoChargedkTScale(0),
    fRhoChargedCMS(0),

    fTPCJet(0),
    fTPCFullJet(0),
    fTPCOnlyJet(0),
    fTPCJetUnbiased(0),
    fTPCkTFullJet(0),
    fEMCalJet(0),
    fEMCalFullJet(0),
    fEMCalPartJet(0),
    fEMCalPartJetUnbiased(0),
    fEMCalkTFullJet(0),

    fIsInitialized(0),
    fRJET(4),
    fnEvents(0),
    fnEventsCharged(0),
    fnDiJetEvents(0),
    fEvent(0),
    fRecoUtil(0),
    fEMCALGeometry(0),
    fCells(0),
    fDoNEF(0),
    fDoNEFSignalOnly(1),
    fSignalTrackBias(0),
    fTrackQA(0),
    fClusterQA(0),
    fCalculateRhoJet(0),
    fDoVertexRCut(0),
    fMCPartLevel(0),
    fDoTHnSparse(0),
    fDoJetRhoDensity(0),
    fDo3DHistos(0),
    fEMCalPhiMin(1.39626),
    fEMCalPhiMax(3.26377),
    fEMCalPhiTotal(1.86750),
    fEMCalEtaMin(-0.7),
    fEMCalEtaMax(0.7),
    fEMCalEtaTotal(1.4),
    fEMCalArea(2.61450),
    fTPCPhiMin(0),
    fTPCPhiMax(6.28319),
    fTPCPhiTotal(6.28319),
    fTPCEtaMin(-0.9),
    fTPCEtaMax(0.9),
    fTPCEtaTotal(1.8),
    fTPCArea(11.30973),
    fParticlePtLow(0.0),
    fParticlePtUp(200.0),
    fParticlePtBins(200),
    fJetR(0.4),
    fJetRAccept(0.4),
    fFullEDJetR(0.6),
    fChargedEDJetR(0.8),
    fJetRForRho(0.5),
    fJetAreaCutFrac(0.6),
    fJetAreaThreshold(0.30159),
    fnEMCalCells(12288),
    fScaleFactor(1.28),
    fNColl(7),
    fTrackMinPt(0.15),
    fClusterMinPt(0.3),
    fNEFSignalJetCut(0.9),
    fCentralityTag("V0A"),
    fCentralityBins(10),
    fCentralityLow(0),
    fCentralityUp(100),
    fEventCentrality(0),
    fRhoFull(0),
    fRhoCharged(0),
    fnTracks(0),
    fnEMCalTracks(0),
    fnClusters(0),
    fnCaloClusters(0),
    fnAKTFullJets(0),
    fnAKTChargedJets(0),
    fnKTFullJets(0),
    fnKTChargedJets(0),
    fnBckgClusters(0),
    fTPCJetThreshold(5),
    fEMCalJetThreshold(5),
    fVertexWindow(10),
    fVertexMaxR(1),
    fTrackName(0),
    fClusName(0),
    fkTChargedName(0),
    fAkTChargedName(0),
    fkTFullName(0),
    fAkTFullName(0),
    fOrgTracks(0),
    fOrgClusters(0),
    fmyAKTFullJets(0),
    fmyAKTChargedJets(0),
    fmyKTFullJets(0),
    fmyKTChargedJets(0),
    fmyTracks(0),
    fmyClusters(0),
    fEMCalRCBckgFluc(0),
    fTPCRCBckgFluc(0),
    fEMCalRCBckgFlucSignal(0),
    fTPCRCBckgFlucSignal(0),
    fEMCalRCBckgFlucNColl(0),
    fTPCRCBckgFlucNColl(0)
{
    // Dummy constructor ALWAYS needed for I/O.
    fVertex[0]=0.0,fVertex[1]=0.0,fVertex[2]=0.0;
}

//________________________________________________________________________
AliAnalysisTaskFullpAJets::AliAnalysisTaskFullpAJets(const char *name) :
    AliAnalysisTaskEmcalJet(name),

    fOutput(0),
    flTrack(0),
    flCluster(0),
    fhTrackPt(0),
    fhTrackEta(0),
    fhTrackPhi(0),
    fhGlobalTrackPt(0),
    fhGlobalTrackEta(0),
    fhGlobalTrackPhi(0),
    fhComplementaryTrackPt(0),
    fhComplementaryTrackEta(0),
    fhComplementaryTrackPhi(0),
    fhClusterPt(0),
    fhClusterEta(0),
    fhClusterPhi(0),
    fhCentrality(0),
    fhEMCalCellCounts(0),

    fhChargeAndNeutralEvents(0),
    fhChargeOnlyEvents(0),
    fhNeutralOnlyEvents(0),
    fhNothingEvents(0),
    fhEMCalChargeAndNeutralEvents(0),
    fhEMCalChargeOnlyEvents(0),
    fhEMCalNeutralOnlyEvents(0),
    fhEMCalNothingEvents(0),

    fhTrackEtaPhi(0),
    fhTrackPhiPt(0),
    fhTrackEtaPt(0),
    fhGlobalTrackEtaPhi(0),
    fhGlobalTrackPhiPt(0),
    fhGlobalTrackEtaPt(0),
    fhComplementaryTrackEtaPhi(0),
    fhComplementaryTrackPhiPt(0),
    fhComplementaryTrackEtaPt(0),
    fhClusterEtaPhi(0),
    fhClusterPhiPt(0),
    fhClusterEtaPt(0),
    fhEMCalEventMult(0),
    fhTPCEventMult(0),
    fhEMCalTrackEventMult(0),

    fpEMCalEventMult(0),
    fpTPCEventMult(0),

    fpTrackPtProfile(0),
    fpClusterPtProfile(0),

    fpFullJetEDProfile(0),
    fpChargedJetEDProfile(0),
    fpChargedJetEDProfileScaled(0),

    fTPCRawJets(0),
    fEMCalRawJets(0),
    fRhoChargedCMSScale(0),
    fRhoChargedScale(0),
    fRhoFull0(0),
    fRhoFull1(0),
    fRhoFull2(0),
    fRhoFullN(0),
    fRhoFullDijet(0),
    fRhoFullkT(0),
    fRhoFullCMS(0),
    fRhoCharged0(0),
    fRhoCharged1(0),
    fRhoCharged2(0),
    fRhoChargedN(0),
    fRhoChargedkT(0),
    fRhoChargedkTScale(0),
    fRhoChargedCMS(0),

    fTPCJet(0),
    fTPCFullJet(0),
    fTPCOnlyJet(0),
    fTPCJetUnbiased(0),
    fTPCkTFullJet(0),
    fEMCalJet(0),
    fEMCalFullJet(0),
    fEMCalPartJet(0),
    fEMCalPartJetUnbiased(0),
    fEMCalkTFullJet(0),

    fIsInitialized(0),
    fRJET(4),
    fnEvents(0),
    fnEventsCharged(0),
    fnDiJetEvents(0),
    fEvent(0),
    fRecoUtil(0),
    fEMCALGeometry(0),
    fCells(0),
    fDoNEF(0),
    fDoNEFSignalOnly(1),
    fSignalTrackBias(0),
    fTrackQA(0),
    fClusterQA(0),
    fCalculateRhoJet(0),
    fDoVertexRCut(0),
    fMCPartLevel(0),
    fDoTHnSparse(0),
    fDoJetRhoDensity(0),
    fDo3DHistos(0),
    fEMCalPhiMin(1.39626),
    fEMCalPhiMax(3.26377),
    fEMCalPhiTotal(1.86750),
    fEMCalEtaMin(-0.7),
    fEMCalEtaMax(0.7),
    fEMCalEtaTotal(1.4),
    fEMCalArea(2.61450),
    fTPCPhiMin(0),
    fTPCPhiMax(6.28319),
    fTPCPhiTotal(6.28319),
    fTPCEtaMin(-0.9),
    fTPCEtaMax(0.9),
    fTPCEtaTotal(1.8),
    fTPCArea(11.30973),
    fParticlePtLow(0.0),
    fParticlePtUp(200.0),
    fParticlePtBins(2000),
    fJetR(0.4),
    fJetRAccept(0.4),
    fFullEDJetR(0.6),
    fChargedEDJetR(0.8),
    fJetRForRho(0.5),
    fJetAreaCutFrac(0.6),
    fJetAreaThreshold(0.30159),
    fnEMCalCells(12288),
    fScaleFactor(1.28),
    fNColl(7),
    fTrackMinPt(0.15),
    fClusterMinPt(0.3),
    fNEFSignalJetCut(0.9),
    fCentralityTag("V0A"),
    fCentralityBins(10),
    fCentralityLow(0),
    fCentralityUp(100),
    fEventCentrality(0),
    fRhoFull(0),
    fRhoCharged(0),
    fnTracks(0),
    fnEMCalTracks(0),
    fnClusters(0),
    fnCaloClusters(0),
    fnAKTFullJets(0),
    fnAKTChargedJets(0),
    fnKTFullJets(0),
    fnKTChargedJets(0),
    fnBckgClusters(0),
    fTPCJetThreshold(5),
    fEMCalJetThreshold(5),
    fVertexWindow(10),
    fVertexMaxR(1),
    fTrackName(0),
    fClusName(0),
    fkTChargedName(0),
    fAkTChargedName(0),
    fkTFullName(0),
    fAkTFullName(0),
    fOrgTracks(0),
    fOrgClusters(0),
    fmyAKTFullJets(0),
    fmyAKTChargedJets(0),
    fmyKTFullJets(0),
    fmyKTChargedJets(0),
    fmyTracks(0),
    fmyClusters(0),
    fEMCalRCBckgFluc(0),
    fTPCRCBckgFluc(0),
    fEMCalRCBckgFlucSignal(0),
    fTPCRCBckgFlucSignal(0),
    fEMCalRCBckgFlucNColl(0),
    fTPCRCBckgFlucNColl(0)
{
    // Constructor
    // Define input and output slots here (never in the dummy constructor)
    // Input slot #0 works with a TChain - it is connected to the default input container
    // Output slot #1 writes into a TH1 container
    fVertex[0]=0.0,fVertex[1]=0.0,fVertex[2]=0.0;

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
}

//________________________________________________________________________
void AliAnalysisTaskFullpAJets::UserCreateOutputObjects()
{
    // Create histograms
    // Called once (on the worker node)
    fIsInitialized=kFALSE;
    fOutput = new TList();
    fOutput->SetOwner();  // IMPORTANT!
   
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
    fJetRForRho=0.5;
    
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
    
    fParticlePtLow=0.0;
    fParticlePtUp=200.0;
    fParticlePtBins=Int_t(fParticlePtUp-fParticlePtLow);
    
    fCentralityBins=10;
    fCentralityLow=0.0;
    fCentralityUp=100.0;
    Int_t CentralityBinMult=10;
    
    fJetAreaCutFrac =0.6; // Fudge factor for selecting on jets with threshold Area or higher
    fJetAreaThreshold=fJetAreaCutFrac*TMath::Pi()*fJetR*fJetR;
    fTPCJetThreshold=5.0; // Threshold required for an Anti-kT Charged jet to be considered a "true" jet in TPC
    fEMCalJetThreshold=5.0; // Threshold required for an Anti-kT Charged+Neutral jet to be considered a "true" jet in EMCal
    fVertexWindow=10.0;
    fVertexMaxR=1.0;
    
    fnBckgClusters=1;
    fEMCalRCBckgFluc = new Double_t[fnBckgClusters];
    fTPCRCBckgFluc = new Double_t[fnBckgClusters];
    fEMCalRCBckgFlucSignal = new Double_t[fnBckgClusters];
    fTPCRCBckgFlucSignal = new Double_t[fnBckgClusters];
    fEMCalRCBckgFlucNColl = new Double_t[fnBckgClusters];
    fTPCRCBckgFlucNColl = new Double_t[fnBckgClusters];
    for (Int_t i=0;i<fnBckgClusters;i++)
    {
        fEMCalRCBckgFluc[i]=0.0;
        fTPCRCBckgFluc[i]=0.0;
        fEMCalRCBckgFlucSignal[i]=0.0;
        fTPCRCBckgFlucSignal[i]=0.0;
        fEMCalRCBckgFlucNColl[i]=0.0;
        fTPCRCBckgFlucNColl[i]=0.0;
    }

    fnEMCalCells=12288;  // sMods 1-10 have 24x48 cells, sMods 11&12 have 8x48 cells...
    
    Int_t TCBins=100;
    Int_t multBins=200;
    Double_t multLow=0;
    Double_t multUp=200;

    // Track QA Plots
    if (fTrackQA==kTRUE)
    {
        flTrack = new TList();
        flTrack->SetName("TrackQA");
        
        // Hybrid Tracks
        fhTrackPt = new TH1F("fhTrackPt","p_{T} distribution of tracks in event",10*fParticlePtBins,fParticlePtLow,fParticlePtUp);
        fhTrackPt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
        fhTrackPt->GetYaxis()->SetTitle("1/N_{Events} dN/dp_{T}");
        fhTrackPt->Sumw2();
        
        fhTrackPhi = new TH1F("fhTrackPhi","#varphi distribution of tracks in event",TCBins,fTPCPhiMin,fTPCPhiMax);
        fhTrackPhi->GetXaxis()->SetTitle("#varphi");
        fhTrackPhi->GetYaxis()->SetTitle("1/N_{Events} dN/d#varphi");
        fhTrackPhi->Sumw2();
        
        fhTrackEta = new TH1F("fhTrackEta","#eta distribution of tracks in event",TCBins,fTPCEtaMin,fTPCEtaMax);
        fhTrackEta->GetXaxis()->SetTitle("#eta");
        fhTrackEta->GetYaxis()->SetTitle("1/N_{Events} dN/d#eta");
        fhTrackEta->Sumw2();
        
        fhTrackEtaPhi = new TH2F("fhTrackEtaPhi","#eta-#varphi distribution of tracks in event",TCBins,fTPCEtaMin,fTPCEtaMax,TCBins,fTPCPhiMin,fTPCPhiMax);
        fhTrackEtaPhi->GetXaxis()->SetTitle("#eta");
        fhTrackEtaPhi->GetYaxis()->SetTitle("#varphi");
        fhTrackEtaPhi->GetZaxis()->SetTitle("1/N_{Events} dN/d#etad#varphi");
        fhTrackEtaPhi->Sumw2();
        
        fhTrackPhiPt = new TH2F("fhTrackPhiPt","#varphi-p_{T} distribution of tracks in event",TCBins,fTPCPhiMin,fTPCPhiMax,fParticlePtBins,fParticlePtLow,fParticlePtUp);
        fhTrackPhiPt->GetXaxis()->SetTitle("#varphi");
        fhTrackPhiPt->GetYaxis()->SetTitle("p_{T} (GeV/c)");
        fhTrackPhiPt->GetZaxis()->SetTitle("1/N_{Events} dN/d#varphidp_{T}");
        fhTrackPhiPt->Sumw2();
        
        fhTrackEtaPt = new TH2F("fhTrackEtaPt","#eta-p_{T} distribution of tracks in event",TCBins,fTPCEtaMin,fTPCEtaMax,fParticlePtBins,fParticlePtLow,fParticlePtUp);
        fhTrackEtaPt->GetXaxis()->SetTitle("#varphi");
        fhTrackEtaPt->GetYaxis()->SetTitle("p_{T} (GeV/c)");
        fhTrackEtaPt->GetZaxis()->SetTitle("1/N_{Events} dN/d#etadp_{T}");
        fhTrackEtaPt->Sumw2();
        
        // Global Tracks
        fhGlobalTrackPt = new TH1F("fhGlobalTrackPt","Global p_{T} distribution of tracks in event",10*fParticlePtBins,fParticlePtLow,fParticlePtUp);
        fhGlobalTrackPt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
        fhGlobalTrackPt->GetYaxis()->SetTitle("1/N_{Events} dN/dp_{T}");
        fhGlobalTrackPt->Sumw2();
        
        fhGlobalTrackPhi = new TH1F("fhGlobalTrackPhi","Global #varphi distribution of tracks in event",TCBins,fTPCPhiMin,fTPCPhiMax);
        fhGlobalTrackPhi->GetXaxis()->SetTitle("#varphi");
        fhGlobalTrackPhi->GetYaxis()->SetTitle("1/N_{Events} dN/d#varphi");
        fhGlobalTrackPhi->Sumw2();
        
        fhGlobalTrackEta = new TH1F("fhGlobalTrackEta","Global #eta distribution of tracks in event",TCBins,fTPCEtaMin,fTPCEtaMax);
        fhGlobalTrackEta->GetXaxis()->SetTitle("#eta");
        fhGlobalTrackEta->GetYaxis()->SetTitle("1/N_{Events} dN/d#eta");
        fhGlobalTrackEta->Sumw2();
        
        fhGlobalTrackEtaPhi = new TH2F("fhGlobalTrackEtaPhi","Global #eta-#varphi distribution of tracks in event",TCBins,fTPCEtaMin,fTPCEtaMax,TCBins,fTPCPhiMin,fTPCPhiMax);
        fhGlobalTrackEtaPhi->GetXaxis()->SetTitle("#eta");
        fhGlobalTrackEtaPhi->GetYaxis()->SetTitle("#varphi");
        fhGlobalTrackEtaPhi->GetZaxis()->SetTitle("1/N_{Events} dN/d#etad#varphi");
        fhGlobalTrackEtaPhi->Sumw2();
        
        fhGlobalTrackPhiPt = new TH2F("fhGlobalTrackPhiPt","Global #varphi-p_{T} distribution of tracks in event",TCBins,fTPCPhiMin,fTPCPhiMax,fParticlePtBins,fParticlePtLow,fParticlePtUp);
        fhGlobalTrackPhiPt->GetXaxis()->SetTitle("#varphi");
        fhGlobalTrackPhiPt->GetYaxis()->SetTitle("p_{T} (GeV/c)");
        fhGlobalTrackPhiPt->GetZaxis()->SetTitle("1/N_{Events} dN/d#varphidp_{T}");
        fhGlobalTrackPhiPt->Sumw2();
        
        fhGlobalTrackEtaPt = new TH2F("fhGlobalTrackEtaPt","Global #eta-p_{T} distribution of tracks in event",TCBins,fTPCEtaMin,fTPCEtaMax,fParticlePtBins,fParticlePtLow,fParticlePtUp);
        fhGlobalTrackEtaPt->GetXaxis()->SetTitle("#varphi");
        fhGlobalTrackEtaPt->GetYaxis()->SetTitle("p_{T} (GeV/c)");
        fhGlobalTrackEtaPt->GetZaxis()->SetTitle("1/N_{Events} dN/d#etadp_{T}");
        fhGlobalTrackEtaPt->Sumw2();
        
        // Complementary Tracks
        fhComplementaryTrackPt = new TH1F("fhComplementaryTrackPt","Complementary p_{T} distribution of tracks in event",10*fParticlePtBins,fParticlePtLow,fParticlePtUp);
        fhComplementaryTrackPt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
        fhComplementaryTrackPt->GetYaxis()->SetTitle("1/N_{Events} dN/dp_{T}");
        fhComplementaryTrackPt->Sumw2();
        
        fhComplementaryTrackPhi = new TH1F("fhComplementaryTrackPhi","Complementary #varphi distribution of tracks in event",TCBins,fTPCPhiMin,fTPCPhiMax);
        fhComplementaryTrackPhi->GetXaxis()->SetTitle("#varphi");
        fhComplementaryTrackPhi->GetYaxis()->SetTitle("1/N_{Events} dN/d#varphi");
        fhComplementaryTrackPhi->Sumw2();
        
        fhComplementaryTrackEta = new TH1F("fhComplementaryTrackEta","Complementary #eta distribution of tracks in event",TCBins,fTPCEtaMin,fTPCEtaMax);
        fhComplementaryTrackEta->GetXaxis()->SetTitle("#eta");
        fhComplementaryTrackEta->GetYaxis()->SetTitle("1/N_{Events} dN/d#eta");
        fhComplementaryTrackEta->Sumw2();
        
        fhComplementaryTrackEtaPhi = new TH2F("fhComplementaryTrackEtaPhi","Complementary #eta-#varphi distribution of tracks in event",TCBins,fTPCEtaMin,fTPCEtaMax,TCBins,fTPCPhiMin,fTPCPhiMax);
        fhComplementaryTrackEtaPhi->GetXaxis()->SetTitle("#eta");
        fhComplementaryTrackEtaPhi->GetYaxis()->SetTitle("#varphi");
        fhComplementaryTrackEtaPhi->GetZaxis()->SetTitle("1/N_{Events} dN/d#etad#varphi");
        fhComplementaryTrackEtaPhi->Sumw2();
        
        fhComplementaryTrackPhiPt = new TH2F("fhComplementaryTrackPhiPt","Complementary #varphi-p_{T} distribution of tracks in event",TCBins,fTPCPhiMin,fTPCPhiMax,fParticlePtBins,fParticlePtLow,fParticlePtUp);
        fhComplementaryTrackPhiPt->GetXaxis()->SetTitle("#varphi");
        fhComplementaryTrackPhiPt->GetYaxis()->SetTitle("p_{T} (GeV/c)");
        fhComplementaryTrackPhiPt->GetZaxis()->SetTitle("1/N_{Events} dN/d#varphidp_{T}");
        fhComplementaryTrackPhiPt->Sumw2();
        
        fhComplementaryTrackEtaPt = new TH2F("fhComplementaryTrackEtaPt","Complementary #eta-p_{T} distribution of tracks in event",TCBins,fTPCEtaMin,fTPCEtaMax,fParticlePtBins,fParticlePtLow,fParticlePtUp);
        fhComplementaryTrackEtaPt->GetXaxis()->SetTitle("#varphi");
        fhComplementaryTrackEtaPt->GetYaxis()->SetTitle("p_{T} (GeV/c)");
        fhComplementaryTrackEtaPt->GetZaxis()->SetTitle("1/N_{Events} dN/d#etadp_{T}");
        fhComplementaryTrackEtaPt->Sumw2();
        
        fhTPCEventMult = new TH2F("fhTPCEventMult","TPC Event Multiplcity vs Centrality",CentralityBinMult*fCentralityBins,fCentralityLow,fCentralityUp,multBins,multLow,multUp);
        fhTPCEventMult->GetXaxis()->SetTitle(fCentralityTag.Data());
        fhTPCEventMult->GetYaxis()->SetTitle("Multiplicity");
        fhTPCEventMult->GetZaxis()->SetTitle("1/N_{Events} dN/dCentdN_{Charged}");
        fhTPCEventMult->Sumw2();
        
        fhEMCalTrackEventMult = new TH2F("fhEMCalTrackEventMult","EMCal Track Event Multiplcity vs Centrality",CentralityBinMult*fCentralityBins,fCentralityLow,fCentralityUp,multBins,multLow,multUp);
        fhEMCalTrackEventMult->GetXaxis()->SetTitle(fCentralityTag.Data());
        fhEMCalTrackEventMult->GetYaxis()->SetTitle("Multiplicity");
        fhEMCalTrackEventMult->GetZaxis()->SetTitle("1/N_{Events} dN/dCentdN_{Neutral}");
        fhEMCalTrackEventMult->Sumw2();

        fpTPCEventMult = new TProfile("fpTPCEventMult","TPC Event Multiplcity vs Centrality",CentralityBinMult*fCentralityBins,fCentralityLow,fCentralityUp);
        fpTPCEventMult->GetXaxis()->SetTitle(fCentralityTag.Data());
        fpTPCEventMult->GetYaxis()->SetTitle("Multiplicity");

        // QA::2D Energy Density Profiles for Tracks and Clusters
        fpTrackPtProfile = new TProfile2D("fpTrackPtProfile","2D Profile of track pT density throughout the TPC",TCBins,fTPCEtaMin,fTPCEtaMax,TCBins,fTPCPhiMin,fTPCPhiMax);
        fpTrackPtProfile->GetXaxis()->SetTitle("#eta");
        fpTrackPtProfile->GetYaxis()->SetTitle("#varphi");
        fpTrackPtProfile->GetZaxis()->SetTitle("p_{T} density (GeV/Area)");

        flTrack->Add(fhTrackPt);
        flTrack->Add(fhTrackEta);
        flTrack->Add(fhTrackPhi);
        flTrack->Add(fhTrackEtaPhi);
        flTrack->Add(fhTrackPhiPt);
        flTrack->Add(fhTrackEtaPt);
        flTrack->Add(fhGlobalTrackPt);
        flTrack->Add(fhGlobalTrackEta);
        flTrack->Add(fhGlobalTrackPhi);
        flTrack->Add(fhGlobalTrackEtaPhi);
        flTrack->Add(fhGlobalTrackPhiPt);
        flTrack->Add(fhGlobalTrackEtaPt);
        flTrack->Add(fhComplementaryTrackPt);
        flTrack->Add(fhComplementaryTrackEta);
        flTrack->Add(fhComplementaryTrackPhi);
        flTrack->Add(fhComplementaryTrackEtaPhi);
        flTrack->Add(fhComplementaryTrackPhiPt);
        flTrack->Add(fhComplementaryTrackEtaPt);
        flTrack->Add(fhTPCEventMult);
        flTrack->Add(fhEMCalTrackEventMult);
        flTrack->Add(fpTPCEventMult);
        flTrack->Add(fpTrackPtProfile);
        fOutput->Add(flTrack);
    }

    if (fClusterQA==kTRUE)
    {
        flCluster = new TList();
        flCluster->SetName("ClusterQA");

        fhClusterPt = new TH1F("fhClusterPt","p_{T} distribution of clusters in event",10*fParticlePtBins,fParticlePtLow,fParticlePtUp);
        fhClusterPt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
        fhClusterPt->GetYaxis()->SetTitle("1/N_{Events} dN/dp_{T}");
        fhClusterPt->Sumw2();
        
        fhClusterPhi = new TH1F("fhClusterPhi","#varphi distribution of clusters in event",TCBins,fTPCPhiMin,fTPCPhiMax);
        fhClusterPhi->GetXaxis()->SetTitle("#varphi");
        fhClusterPhi->GetYaxis()->SetTitle("1/N_{Events} dN/d#varphi");
        fhClusterPhi->Sumw2();
        
        fhClusterEta = new TH1F("fhClusterEta","#eta distribution of clusters in event",TCBins,fTPCEtaMin,fTPCEtaMax);
        fhClusterEta->GetXaxis()->SetTitle("#eta");
        fhClusterEta->GetYaxis()->SetTitle("1/N_{Events} dN/d#eta");
        fhClusterEta->Sumw2();
        
        fhClusterEtaPhi = new TH2F("fhClusterEtaPhi","#eta-#varphi distribution of clusters in event",TCBins,fEMCalEtaMin,fEMCalEtaMax,TCBins,fEMCalPhiMin,fEMCalPhiMax);
        fhClusterEtaPhi->GetXaxis()->SetTitle("#eta");
        fhClusterEtaPhi->GetYaxis()->SetTitle("#varphi");
        fhClusterEtaPhi->GetZaxis()->SetTitle("1/N_{Events} dN/d#etad#varphi");
        fhClusterEtaPhi->Sumw2();
        
        fhClusterPhiPt = new TH2F("fhClusterPhiPt","#varphi-p_{T} distribution of clusters in event",TCBins,fEMCalPhiMin,fEMCalPhiMax,fParticlePtBins,fParticlePtLow,fParticlePtUp);
        fhClusterPhiPt->GetXaxis()->SetTitle("#varphi");
        fhClusterPhiPt->GetYaxis()->SetTitle("p_{T} (GeV/c)");
        fhClusterPhiPt->GetZaxis()->SetTitle("1/N_{Events} dN/d#varphidp_{T}");
        fhClusterPhiPt->Sumw2();
        
        fhClusterEtaPt = new TH2F("fhClusterEtaPt","#eta-p_{T} distribution of clusters in event",TCBins,fEMCalEtaMin,fEMCalEtaMax,fParticlePtBins,fParticlePtLow,fParticlePtUp);
        fhClusterEtaPt->GetXaxis()->SetTitle("#varphi");
        fhClusterEtaPt->GetYaxis()->SetTitle("p_{T} (GeV/c)");
        fhClusterEtaPt->GetZaxis()->SetTitle("1/N_{Events} dN/d#etadp_{T}");
        fhClusterEtaPt->Sumw2();
        
        fhEMCalEventMult = new TH2F("fhEMCalEventMult","EMCal Event Multiplcity vs Centrality",CentralityBinMult*fCentralityBins,fCentralityLow,fCentralityUp,multBins,multLow,multUp);
        fhEMCalEventMult->GetXaxis()->SetTitle(fCentralityTag.Data());
        fhEMCalEventMult->GetYaxis()->SetTitle("Multiplicity");
        fhEMCalEventMult->GetZaxis()->SetTitle("1/N_{Events} dN/dCentdN_{Neutral}");
        fhEMCalEventMult->Sumw2();
        
        fpEMCalEventMult = new TProfile("fpEMCalEventMult","EMCal Event Multiplcity vs Centrality",CentralityBinMult*fCentralityBins,fCentralityLow,fCentralityUp);
        fpEMCalEventMult->GetXaxis()->SetTitle(fCentralityTag.Data());
        fpEMCalEventMult->GetYaxis()->SetTitle("Multiplicity");
        
        fpClusterPtProfile = new TProfile2D("fpClusterPtProfile","2D Profile of cluster pT density throughout the EMCal",TCBins,fEMCalEtaMin,fEMCalEtaMax,TCBins,fEMCalPhiMin,fEMCalPhiMax);
        fpClusterPtProfile->GetXaxis()->SetTitle("#eta");
        fpClusterPtProfile->GetYaxis()->SetTitle("#varphi");
        fpClusterPtProfile->GetZaxis()->SetTitle("p_{T} density (GeV/Area)");

        fhEMCalCellCounts = new TH1F("fhEMCalCellCounts","Distribtuion of cluster counts across the EMCal",fnEMCalCells,1,fnEMCalCells);
        fhEMCalCellCounts->GetXaxis()->SetTitle("Absoulute Cell Id");
        fhEMCalCellCounts->GetYaxis()->SetTitle("Counts per Event");
        fhEMCalCellCounts->Sumw2();

        flCluster->Add(fhClusterPt);
        flCluster->Add(fhClusterEta);
        flCluster->Add(fhClusterPhi);
        flCluster->Add(fhClusterEtaPhi);
        flCluster->Add(fhClusterPhiPt);
        flCluster->Add(fhClusterEtaPt);
        flCluster->Add(fhEMCalEventMult);
        flCluster->Add(fpEMCalEventMult);
        flCluster->Add(fpClusterPtProfile);
        flCluster->Add(fhEMCalCellCounts);
        fOutput->Add(flCluster);
    }
    
    if (fCalculateRhoJet>=0 && fMCPartLevel==kFALSE) // Default Rho & Raw Jet Spectra
    {
        fEMCalRawJets = new AlipAJetHistos("EMCalRawJets",fCentralityTag.Data(),fDoNEF,fDoNEFSignalOnly,fDoTHnSparse,fDo3DHistos);
        
        fRhoChargedCMSScale = new AlipAJetHistos("RhoChargedCMSScale",fCentralityTag.Data(),fDoNEF,fDoNEFSignalOnly,fDoTHnSparse,fDo3DHistos);
        fRhoChargedCMSScale->SetSignalTrackPtBias(fSignalTrackBias);
        fOutput->Add(fEMCalRawJets->GetOutputHistos());
        fOutput->Add(fRhoChargedCMSScale->GetOutputHistos());

    }
    if (fCalculateRhoJet>=1) // Basic Rho & Raw Jet Spectra
    {
        fRhoChargedScale = new AlipAJetHistos("RhoChargedScale",fCentralityTag.Data(),fDoNEF,fDoNEFSignalOnly,fDoTHnSparse,fDo3DHistos);
        fRhoChargedScale->SetSignalTrackPtBias(fSignalTrackBias);
        
        if (fMCPartLevel==kTRUE)
        {
            fTPCRawJets = new AlipAJetHistos("TPCRawJets",fCentralityTag.Data(),fDoNEF,fDoNEFSignalOnly,fDoTHnSparse,fDo3DHistos);
            fRhoChargedCMS = new AlipAJetHistos("RhoChargedCMS",fCentralityTag.Data(),fDoNEF,fDoNEFSignalOnly,fDoTHnSparse,fDo3DHistos);
            fRhoChargedCMS->SetSignalTrackPtBias(fSignalTrackBias);
            
            fOutput->Add(fTPCRawJets->GetOutputHistos());
            fOutput->Add(fRhoChargedCMS->GetOutputHistos());
        }

        if (fDoJetRhoDensity == kTRUE && fDo3DHistos == kTRUE)
        {
            Double_t fullEDR = fFullEDJetR *100;
            Double_t chargedEDR = fChargedEDJetR *100;
            const Int_t fullEDBins = (Int_t) fullEDR;
            const Int_t chargedEDBins = (Int_t) chargedEDR;
            
            fpFullJetEDProfile = new TProfile3D("fpFullJetEDProfile","Jet ED Profile for #varphi_{jet}#in(#varphi_{EMCal min} + R,#varphi_{EMCal max} - R) & #eta_{jet}#in(#eta_{EMCal min} + R,#eta_{EMCal max} - R)",fParticlePtBins,fParticlePtLow,fParticlePtUp,fCentralityBins,fCentralityLow,fCentralityUp,fullEDBins,0,fFullEDJetR);
            fpFullJetEDProfile->GetXaxis()->SetTitle("p_{T,jet}^{ch+em} (GeV)");
            fpFullJetEDProfile->GetYaxis()->SetTitle(fCentralityTag.Data());
            fpFullJetEDProfile->GetZaxis()->SetTitle("R");
            
            fpChargedJetEDProfile = new TProfile3D("fpChargedJetEDProfile","Charged Jet ED Profile for #eta_{jet}#in(#eta_{TPC min} + R,#eta_{TPC max} - R)",fParticlePtBins,fParticlePtLow,fParticlePtUp,fCentralityBins,fCentralityLow,fCentralityUp,chargedEDBins,0,fChargedEDJetR);
            fpChargedJetEDProfile->GetXaxis()->SetTitle("p_{T,jet}^{ch} (GeV)");
            fpChargedJetEDProfile->GetYaxis()->SetTitle(fCentralityTag.Data());
            fpChargedJetEDProfile->GetZaxis()->SetTitle("R");
            
            fpChargedJetEDProfileScaled = new TProfile3D("fpChargedJetEDProfileScaled","Charged Jet ED Profile x SF for #eta_{jet}#in(#eta_{TPC min} + R,#eta_{TPC max} - R)",fParticlePtBins,fParticlePtLow,fParticlePtUp,fCentralityBins,fCentralityLow,fCentralityUp,chargedEDBins,0,fChargedEDJetR);
            fpChargedJetEDProfileScaled->GetXaxis()->SetTitle("p_{T,jet}^{ch} (GeV)");
            fpChargedJetEDProfileScaled->GetYaxis()->SetTitle(fCentralityTag.Data());
            fpChargedJetEDProfileScaled->GetZaxis()->SetTitle("R");
            
            fOutput->Add(fpFullJetEDProfile);
            fOutput->Add(fpChargedJetEDProfile);
            fOutput->Add(fpChargedJetEDProfileScaled);
        }
        fOutput->Add(fRhoChargedScale->GetOutputHistos());
    }
    if (fCalculateRhoJet>=2 && fMCPartLevel==kFALSE) // Charged Rho & Raw Jet Spectra
    {
        fTPCRawJets = new AlipAJetHistos("TPCRawJets",fCentralityTag.Data(),fDoNEF,fDoNEFSignalOnly,fDoTHnSparse,fDo3DHistos);
        fRhoChargedCMS = new AlipAJetHistos("RhoChargedCMS",fCentralityTag.Data(),fDoNEF,fDoNEFSignalOnly,fDoTHnSparse,fDo3DHistos);
        fRhoChargedCMS->SetSignalTrackPtBias(fSignalTrackBias);

        fOutput->Add(fTPCRawJets->GetOutputHistos());
        fOutput->Add(fRhoChargedCMS->GetOutputHistos());
    }
    if (fCalculateRhoJet>=3 && fMCPartLevel==kFALSE) // Basic Rho & Raw Jet Spectra
    {
        fRhoFull0 = new AlipAJetHistos("RhoFull0",fCentralityTag.Data());
        fRhoFull0->SetSignalTrackPtBias(fSignalTrackBias);
        fRhoFull1 = new AlipAJetHistos("RhoFull1",fCentralityTag.Data());
        fRhoFull1->SetSignalTrackPtBias(fSignalTrackBias);
        fRhoFull2 = new AlipAJetHistos("RhoFull2",fCentralityTag.Data());
        fRhoFull2->SetSignalTrackPtBias(fSignalTrackBias);
        fRhoFullN = new AlipAJetHistos("RhoFullN",fCentralityTag.Data());
        fRhoFullN->SetSignalTrackPtBias(fSignalTrackBias);
        fRhoFullDijet = new AlipAJetHistos("RhoFullDijet",fCentralityTag.Data());
        fRhoFullDijet->SetSignalTrackPtBias(fSignalTrackBias);
        fRhoFullkT = new AlipAJetHistos("RhoFullkT",fCentralityTag.Data());
        fRhoFullkT->SetSignalTrackPtBias(fSignalTrackBias);
        fRhoFullCMS = new AlipAJetHistos("RhoFullCMS",fCentralityTag.Data());
        fRhoFullCMS->SetSignalTrackPtBias(fSignalTrackBias);
        fRhoCharged0 = new AlipAJetHistos("RhoCharged0",fCentralityTag.Data());
        fRhoCharged0->SetSignalTrackPtBias(fSignalTrackBias);
        fRhoCharged1 = new AlipAJetHistos("RhoCharged1",fCentralityTag.Data());
        fRhoCharged1->SetSignalTrackPtBias(fSignalTrackBias);
        fRhoCharged2 = new AlipAJetHistos("RhoCharged2",fCentralityTag.Data());
        fRhoCharged2->SetSignalTrackPtBias(fSignalTrackBias);
        fRhoChargedN = new AlipAJetHistos("RhoChargedN",fCentralityTag.Data());
        fRhoChargedN->SetSignalTrackPtBias(fSignalTrackBias);
        fRhoChargedkT = new AlipAJetHistos("RhoChargedkT",fCentralityTag.Data());
        fRhoChargedkT->SetSignalTrackPtBias(fSignalTrackBias);
        fRhoChargedkTScale = new AlipAJetHistos("RhoChargedkTScale",fCentralityTag.Data());
        fRhoChargedkTScale->SetSignalTrackPtBias(fSignalTrackBias);

        fOutput->Add(fRhoFull0->GetOutputHistos());
        fOutput->Add(fRhoFull1->GetOutputHistos());
        fOutput->Add(fRhoFull2->GetOutputHistos());
        fOutput->Add(fRhoFullN->GetOutputHistos());
        fOutput->Add(fRhoFullDijet->GetOutputHistos());
        fOutput->Add(fRhoFullkT->GetOutputHistos());
        fOutput->Add(fRhoFullCMS->GetOutputHistos());
        fOutput->Add(fRhoCharged0->GetOutputHistos());
        fOutput->Add(fRhoCharged1->GetOutputHistos());
        fOutput->Add(fRhoCharged2->GetOutputHistos());
        fOutput->Add(fRhoChargedN->GetOutputHistos());
        fOutput->Add(fRhoChargedkT->GetOutputHistos());
        fOutput->Add(fRhoChargedkTScale->GetOutputHistos());
    }
    
    fhCentrality = new TH1F("fhCentrality","Event Centrality Distribution",fCentralityBins*CentralityBinMult,fCentralityLow,fCentralityUp);
    fhCentrality->GetXaxis()->SetTitle(fCentralityTag.Data());
    fhCentrality->GetYaxis()->SetTitle("1/N_{Events}");
    fhCentrality->Sumw2();
    
    fOutput->Add(fhCentrality);

    if (fMCPartLevel == kFALSE)
    {
        fhChargeAndNeutralEvents = new TH1F("fhChargeAndNeutralEvents","Events with n_{tracks}>0 & n_{clusters}>0 vs Centrality",fCentralityBins*CentralityBinMult,fCentralityLow,fCentralityUp);
        fhChargeAndNeutralEvents->GetXaxis()->SetTitle(fCentralityTag.Data());
        fhChargeAndNeutralEvents->GetYaxis()->SetTitle("1/N_{Events}");
        fhChargeAndNeutralEvents->Sumw2();
        
        fhChargeOnlyEvents = new TH1F("fhChargeOnlyEvents","Events with n_{tracks}>0 & n_{clusters}=0 vs Centrality",fCentralityBins*CentralityBinMult,fCentralityLow,fCentralityUp);
        fhChargeOnlyEvents->GetXaxis()->SetTitle(fCentralityTag.Data());
        fhChargeOnlyEvents->GetYaxis()->SetTitle("1/N_{Events}");
        fhChargeOnlyEvents->Sumw2();
        
        fhNeutralOnlyEvents = new TH1F("fhNeutralOnlyEvents","Events with n_{tracks}=0 & n_{clusters}>0 vs Centrality",fCentralityBins*CentralityBinMult,fCentralityLow,fCentralityUp);
        fhNeutralOnlyEvents->GetXaxis()->SetTitle(fCentralityTag.Data());
        fhNeutralOnlyEvents->GetYaxis()->SetTitle("1/N_{Events}");
        fhNeutralOnlyEvents->Sumw2();
        
        fhNothingEvents = new TH1F("fhNothingEvents","Events with n_{tracks}=n_{clusters}=0 vs Centrality",fCentralityBins*CentralityBinMult,fCentralityLow,fCentralityUp);
        fhNothingEvents->GetXaxis()->SetTitle(fCentralityTag.Data());
        fhNothingEvents->GetYaxis()->SetTitle("1/N_{Events}");
        fhNothingEvents->Sumw2();
        
        fhEMCalChargeAndNeutralEvents = new TH1F("fhEMCalChargeAndNeutralEvents","Events with n_{emcal tracks}>0 & n_{clusters}>0 vs Centrality",fCentralityBins*CentralityBinMult,fCentralityLow,fCentralityUp);
        fhEMCalChargeAndNeutralEvents->GetXaxis()->SetTitle(fCentralityTag.Data());
        fhEMCalChargeAndNeutralEvents->GetYaxis()->SetTitle("1/N_{Events}");
        fhEMCalChargeAndNeutralEvents->Sumw2();
        
        fhEMCalChargeOnlyEvents = new TH1F("fhEMCalChargeOnlyEvents","Events with n_{emcal tracks}>0 & n_{clusters}=0 vs Centrality",fCentralityBins*CentralityBinMult,fCentralityLow,fCentralityUp);
        fhEMCalChargeOnlyEvents->GetXaxis()->SetTitle(fCentralityTag.Data());
        fhEMCalChargeOnlyEvents->GetYaxis()->SetTitle("1/N_{Events}");
        fhEMCalChargeOnlyEvents->Sumw2();
        
        fhEMCalNeutralOnlyEvents = new TH1F("fhEMCalNeutralOnlyEvents","Events with n_{tracks}>0 & n_{emcal tracks}=0 & n_{clusters}>0 vs Centrality",fCentralityBins*CentralityBinMult,fCentralityLow,fCentralityUp);
        fhEMCalNeutralOnlyEvents->GetXaxis()->SetTitle(fCentralityTag.Data());
        fhEMCalNeutralOnlyEvents->GetYaxis()->SetTitle("1/N_{Events}");
        fhEMCalNeutralOnlyEvents->Sumw2();
        
        fhEMCalNothingEvents = new TH1F("fhEMCalNothingEvents","Events with n_{emcal tracks}=n_{clusters}=0 vs Centrality",fCentralityBins*CentralityBinMult,fCentralityLow,fCentralityUp);
        fhEMCalNothingEvents->GetXaxis()->SetTitle(fCentralityTag.Data());
        fhEMCalNothingEvents->GetYaxis()->SetTitle("1/N_{Events}");
        fhEMCalNothingEvents->Sumw2();
        
        fOutput->Add(fhChargeAndNeutralEvents);
        fOutput->Add(fhChargeOnlyEvents);
        fOutput->Add(fhNeutralOnlyEvents);
        fOutput->Add(fhNothingEvents);
        fOutput->Add(fhEMCalChargeAndNeutralEvents);
        fOutput->Add(fhEMCalChargeOnlyEvents);
        fOutput->Add(fhEMCalNeutralOnlyEvents);
        fOutput->Add(fhEMCalNothingEvents);
    }
    
    // Post data for ALL output slots >0 here,
    // To get at least an empty histogram
    // 1 is the outputnumber of a certain weg of task 1
    PostData(1, fOutput);
}

void AliAnalysisTaskFullpAJets::UserExecOnce()
{
    // Get the event tracks
    fOrgTracks = dynamic_cast <TClonesArray*>(InputEvent()->FindListObject(fTrackName));
    
    // Get the event caloclusters
    fOrgClusters = dynamic_cast <TClonesArray*>(InputEvent()->FindListObject(fClusName));
    
    // Get charged jets
    fmyKTChargedJets = dynamic_cast <TClonesArray*>(InputEvent()->FindListObject(fkTChargedName));
    fmyAKTChargedJets = dynamic_cast <TClonesArray*>(InputEvent()->FindListObject(fAkTChargedName));

    // Get the full jets
    fmyKTFullJets = dynamic_cast <TClonesArray*>(InputEvent()->FindListObject(fkTFullName));
    fmyAKTFullJets = dynamic_cast <TClonesArray*>(InputEvent()->FindListObject(fAkTFullName));

    fIsInitialized=kTRUE;
}
//________________________________________________________________________
Bool_t AliAnalysisTaskFullpAJets::Run()
{
    if (fIsInitialized==kFALSE)
    {
        UserExecOnce();
    }

    // Get pointer to reconstructed event
    fEvent = InputEvent();
    if (!fEvent)
    {
        AliError("Pointer == 0, this can not happen!");
        return 0;
    }

    AliESDEvent* esd = dynamic_cast<AliESDEvent*>(fEvent);
    AliAODEvent* aod = dynamic_cast<AliAODEvent*>(fEvent);
    
    fRecoUtil = new AliEMCALRecoUtils();
    fEMCALGeometry = AliEMCALGeometry::GetInstance();
    
    if (esd)
    {
        fEventCentrality=esd->GetCentrality()->GetCentralityPercentile(fCentralityTag.Data());
        
        if (fDoVertexRCut == kTRUE)
        {
            // Vertex R cut not done in AliAnalysisTaskEmcal
            if (TMath::Sqrt(TMath::Power(esd->GetPrimaryVertex()->GetX(),2)+TMath::Power(esd->GetPrimaryVertex()->GetY(),2))>fVertexMaxR)
            {
                DeleteJetData(0);
                return 0;
            }
        }

        esd->GetPrimaryVertex()->GetXYZ(fVertex);
        fCells = (AliVCaloCells*) esd->GetEMCALCells();
        fnCaloClusters = esd->GetNumberOfCaloClusters();
    }
    else if (aod)
    {
        fEventCentrality=aod->GetCentrality()->GetCentralityPercentile(fCentralityTag.Data());
        
        if (fDoVertexRCut == kTRUE)
        {
            // Vertex R cut not done in AliAnalysisTaskEmcal
            if (TMath::Sqrt(TMath::Power(aod->GetPrimaryVertex()->GetX(),2)+TMath::Power(aod->GetPrimaryVertex()->GetY(),2))>fVertexMaxR)
            {
                DeleteJetData(0);
                return 0;
            }
        }
        
        aod->GetPrimaryVertex()->GetXYZ(fVertex);
        fCells = (AliVCaloCells*) aod->GetEMCALCells();
        fnCaloClusters = aod->GetNumberOfCaloClusters();
    }
    else
    {
        AliError("Cannot get AOD/ESD event. Rejecting Event");
        DeleteJetData(0);
        return 0;
    }
    
    // Make sure Centrality isn't exactly 100% (to avoid bin filling errors in profile plots. Make it 99.99
    if (fEventCentrality>99.99)
    {
        fEventCentrality=99.99;
    }
    fhCentrality->Fill(fEventCentrality); // Counts total events

    
    // Count Events
    EventCounts();
    
    // Reject any event that doesn't have any tracks, i.e. TPC is off
    if (fnTracks<1)
    {
        if (fTrackQA==kTRUE)
        {
            fhTPCEventMult->Fill(fEventCentrality,0.0);
            fpTPCEventMult->Fill(fEventCentrality,0.0);
            fhEMCalTrackEventMult->Fill(fEventCentrality,0.0);
        }
        AliWarning("No PicoTracks, Rejecting Event");
        DeleteJetData(1);
        PostData(1,fOutput);
        return 0;
    }
    
    if (fnClusters<1)
    {
        //AliInfo("No Corrected CaloClusters, using only charged jets");
       
        if (fTrackQA==kTRUE)
        {
            TrackHisto();
        }
        InitChargedJets();
        GenerateTPCRandomConesPt();
        
        if (fClusterQA==kTRUE)
        {
            fhEMCalEventMult->Fill(fEventCentrality,0.0);
            fpEMCalEventMult->Fill(fEventCentrality,0.0);
        }

        // Rho's
        if (fCalculateRhoJet>=1)
        {
            if (fDoJetRhoDensity == kTRUE && fDo3DHistos == kTRUE)
            {
                ChargedJetEnergyDensityProfile();
            }
            if (fMCPartLevel==kTRUE)
            {
                EstimateChargedRhoCMS();
            }
        }
        if (fCalculateRhoJet>=2 && fMCPartLevel==kFALSE)
        {
            EstimateChargedRhoCMS();
        }
        if (fCalculateRhoJet>=3 && fMCPartLevel==kFALSE)
        {
            EstimateChargedRho0();
            EstimateChargedRho1();
            EstimateChargedRho2();
            EstimateChargedRhoN();
            EstimateChargedRhokT();
        }
        
        DeleteJetData(2);
        fnEventsCharged++;
        PostData(1,fOutput);
        return 0;
    }

    if (fTrackQA==kTRUE)
    {
        TrackHisto();
    }
    if (fClusterQA==kTRUE)
    {
        ClusterHisto();
    }
    
    // Prep the jets
    InitChargedJets();
    GenerateTPCRandomConesPt();

    InitFullJets();
    GenerateEMCalRandomConesPt();
    
    if (fCalculateRhoJet>=0)
    {
        EstimateChargedRhoCMSScale();
    }
    if (fCalculateRhoJet>=1)
    {
        EstimateChargedRhoScale();
        if (fDoJetRhoDensity == kTRUE && fDo3DHistos == kTRUE)
        {
            ChargedJetEnergyDensityProfile();
            FullJetEnergyDensityProfile();
        }
    }
    if (fCalculateRhoJet>=2)
    {
        EstimateChargedRhoCMS();
    }
    if (fCalculateRhoJet>=3)
    {
        EstimateChargedRho0();
        EstimateChargedRho1();
        EstimateChargedRho2();
        EstimateChargedRhoN();
        EstimateChargedRhokT();
        
        EstimateFullRho0();
        EstimateFullRho1();
        EstimateFullRho2();
        EstimateFullRhoN();
        EstimateFullRhokT();
        EstimateFullRhoCMS();
        
        EstimateChargedRhokTScale();

        // Dijet
        if (IsDiJetEvent()==kTRUE)
        {
            EstimateFullRhoDijet();
        }
    }

    // Delete Dynamic Arrays
    DeleteJetData(3);
    fnEvents++;
    PostData(1,fOutput);
    return 1;
}

//________________________________________________________________________
void AliAnalysisTaskFullpAJets::Terminate(Option_t *) //specify what you want to have done
{
    // Called once at the end of the query. Done nothing here.
}

/////////////////////////////////////////////////////////////////////////////////////////
/////////////////     User Defined Sub_Routines   ///////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////

void AliAnalysisTaskFullpAJets::TrackCuts()
{
    // Fill a TObjArray with the tracks from a TClonesArray which grabs the picotracks.
    Int_t i;
    Int_t j=0;
    
    fmyTracks = new TObjArray();
    for (i=0;i<fOrgTracks->GetEntries();i++)
    {
        AliVTrack* vtrack = (AliVTrack*) fOrgTracks->At(i);
        if (vtrack->Pt()>=fTrackMinPt)
        {
            fmyTracks->Add(vtrack);
            if (IsInEMCal(vtrack->Phi(),vtrack->Eta()) == kTRUE)
            {
                j++;
            }
        }
    }
    fnTracks = fmyTracks->GetEntries();
    fnEMCalTracks = j;
}

void AliAnalysisTaskFullpAJets::ClusterCuts()
{
    // Fill a TObjArray with the clusters from a TClonesArray which grabs the caloclusterscorr.
    Int_t i;
    
    fmyClusters = new TObjArray();
    TLorentzVector *cluster_vec = new TLorentzVector;
    for (i=0;i<fOrgClusters->GetEntries();i++)
    {
        AliVCluster* vcluster = (AliVCluster*) fOrgClusters->At(i);
        vcluster->GetMomentum(*cluster_vec,fVertex);
        
        if (cluster_vec->Pt()>=fClusterMinPt && vcluster->IsEMCAL()==kTRUE)
        {
            fmyClusters->Add(vcluster);
        }
    }
    fnClusters = fmyClusters->GetEntries();
    delete cluster_vec;
}

void AliAnalysisTaskFullpAJets::EventCounts()
{
    TrackCuts();
    if (fMCPartLevel == kTRUE)
    {
        return;
    }
    
    ClusterCuts();
    if (fnTracks>0)
    {
        if (fnClusters>0)
        {
            fhChargeAndNeutralEvents->Fill(fEventCentrality);
        }
        else
        {
            fhChargeOnlyEvents->Fill(fEventCentrality);
        }
        if (fnEMCalTracks>0)
        {
            if (fnClusters>0)
            {
                fhEMCalChargeAndNeutralEvents->Fill(fEventCentrality);
            }
            else
            {
                fhEMCalChargeOnlyEvents->Fill(fEventCentrality);
            }
        }
        else
        {
            if (fnClusters>0)
            {
                fhEMCalNeutralOnlyEvents->Fill(fEventCentrality);
            }
            else
            {
                fhEMCalNothingEvents->Fill(fEventCentrality);
            }
        }
    }
    else
    {
        if (fnClusters>0)
        {
            fhNeutralOnlyEvents->Fill(fEventCentrality);
        }
        else
        {
            fhNothingEvents->Fill(fEventCentrality);
        }
    }
}

void AliAnalysisTaskFullpAJets::TrackHisto()
{
    // Fill track histograms: Phi,Eta,Pt
    Int_t i,j;
    Int_t TCBins=100;
    TH2F *hdummypT= new TH2F("hdummypT","",TCBins,fTPCEtaMin,fTPCEtaMax,TCBins,fTPCPhiMin,fTPCPhiMax);  //!

    for (i=0;i<fnTracks;i++)
    {
        AliPicoTrack* vtrack = (AliPicoTrack*) fmyTracks->At(i);
        fhTrackPt->Fill(vtrack->Pt());
        fhTrackEta->Fill(vtrack->Eta());
        fhTrackPhi->Fill(vtrack->Phi());
        fhTrackEtaPhi->Fill(vtrack->Eta(),vtrack->Phi());
        fhTrackPhiPt->Fill(vtrack->Phi(),vtrack->Pt());
        fhTrackEtaPt->Fill(vtrack->Eta(),vtrack->Pt());
        
        // Fill Associated Track Distributions for AOD QA Productions
        // Global Tracks
        if (vtrack->GetTrackType()==0)
        {
            fhGlobalTrackPt->Fill(vtrack->Pt());
            fhGlobalTrackEta->Fill(vtrack->Eta());
            fhGlobalTrackPhi->Fill(vtrack->Phi());
            fhGlobalTrackEtaPhi->Fill(vtrack->Eta(),vtrack->Phi());
            fhGlobalTrackPhiPt->Fill(vtrack->Phi(),vtrack->Pt());
            fhGlobalTrackEtaPt->Fill(vtrack->Eta(),vtrack->Pt());
        }
        // Complementary Tracks
        else if (vtrack->GetTrackType()==1)
        {
            fhComplementaryTrackPt->Fill(vtrack->Pt());
            fhComplementaryTrackEta->Fill(vtrack->Eta());
            fhComplementaryTrackPhi->Fill(vtrack->Phi());
            fhComplementaryTrackEtaPhi->Fill(vtrack->Eta(),vtrack->Phi());
            fhComplementaryTrackPhiPt->Fill(vtrack->Phi(),vtrack->Pt());
            fhComplementaryTrackEtaPt->Fill(vtrack->Eta(),vtrack->Pt());
        }
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
    TH2F *hdummypT= new TH2F("hdummypT","",TCBins,fEMCalEtaMin,fEMCalEtaMax,TCBins,fEMCalPhiMin,fEMCalPhiMax);  //!
    Int_t myCellID=-2;
    TLorentzVector *cluster_vec = new TLorentzVector;
    
    for (i=0;i<fnClusters;i++)
    {
        AliVCluster* vcluster = (AliVCluster*) fmyClusters->At(i);
        vcluster->GetMomentum(*cluster_vec,fVertex);
        
        fhClusterPt->Fill(cluster_vec->Pt());
        fhClusterEta->Fill(cluster_vec->Eta());
        fhClusterPhi->Fill(cluster_vec->Phi());
        fhClusterEtaPhi->Fill(cluster_vec->Eta(),cluster_vec->Phi());
        fhClusterPhiPt->Fill(cluster_vec->Phi(),cluster_vec->Pt());
        fhClusterEtaPt->Fill(cluster_vec->Eta(),cluster_vec->Pt());
        hdummypT->Fill(cluster_vec->Eta(),cluster_vec->Phi(),cluster_vec->Pt());
        fEMCALGeometry->GetAbsCellIdFromEtaPhi(cluster_vec->Eta(),cluster_vec->Phi(),myCellID);
        fhEMCalCellCounts->Fill(myCellID);
    }
    for (i=1;i<=TCBins;i++)
    {
        for (j=1;j<=TCBins;j++)
        {
            fpClusterPtProfile->Fill(hdummypT->GetXaxis()->GetBinCenter(i),hdummypT->GetYaxis()->GetBinCenter(j),fEMCalArea*TMath::Power(TCBins,-2)*hdummypT->GetBinContent(i,j));
        }
    }
    delete hdummypT;
    delete cluster_vec;
}

void AliAnalysisTaskFullpAJets::InitChargedJets()
{
    // Preliminary Jet Placement and Selection Cuts
    Int_t i;
    
    fnAKTChargedJets = fmyAKTChargedJets->GetEntries();
    fnKTChargedJets = fmyKTChargedJets->GetEntries();
    
    fTPCJet = new AlipAJetData("fTPCJet",kFALSE,fnAKTChargedJets);
    fTPCFullJet = new AlipAJetData("fTPCFullJet",kFALSE,fnAKTChargedJets);
    fTPCOnlyJet = new AlipAJetData("fTPCOnlyJet",kFALSE,fnAKTChargedJets);
    fTPCJetUnbiased = new AlipAJetData("fTPCJetUnbiased",kFALSE,fnAKTChargedJets);

    fTPCJet->SetSignalCut(fTPCJetThreshold);
    fTPCJet->SetAreaCutFraction(fJetAreaCutFrac);
    fTPCJet->SetJetR(fJetR);
    fTPCJet->SetSignalTrackPtBias(fSignalTrackBias);
    fTPCFullJet->SetSignalCut(fTPCJetThreshold);
    fTPCFullJet->SetAreaCutFraction(fJetAreaCutFrac);
    fTPCFullJet->SetJetR(fJetR);
    fTPCFullJet->SetSignalTrackPtBias(fSignalTrackBias);
    fTPCOnlyJet->SetSignalCut(fTPCJetThreshold);
    fTPCOnlyJet->SetAreaCutFraction(fJetAreaCutFrac);
    fTPCOnlyJet->SetJetR(fJetR);
    fTPCOnlyJet->SetSignalTrackPtBias(fSignalTrackBias);
    fTPCJetUnbiased->SetSignalCut(fTPCJetThreshold);
    fTPCJetUnbiased->SetAreaCutFraction(fJetAreaCutFrac);
    fTPCJetUnbiased->SetJetR(fJetR);
    fTPCJetUnbiased->SetSignalTrackPtBias(kFALSE);
    
    // Initialize Jet Data
    for (i=0;i<fnAKTChargedJets;i++)
    {
        AliEmcalJet *myJet =(AliEmcalJet*) fmyAKTChargedJets->At(i);
        
        fTPCJet->SetIsJetInArray(IsInTPC(fJetR,myJet->Phi(),myJet->Eta(),kFALSE),i);
        fTPCFullJet->SetIsJetInArray(IsInTPC(fJetRAccept,myJet->Phi(),myJet->Eta(),kTRUE),i);
        fTPCOnlyJet->SetIsJetInArray(IsInTPCFull(fJetR,myJet->Phi(),myJet->Eta()),i);
        fTPCJetUnbiased->SetIsJetInArray(IsInTPC(fJetR,myJet->Phi(),myJet->Eta(),kFALSE),i);
    }
    
    fTPCJet->InitializeJetData(fmyAKTChargedJets,fnAKTChargedJets);
    fTPCFullJet->InitializeJetData(fmyAKTChargedJets,fnAKTChargedJets);
    fTPCOnlyJet->InitializeJetData(fmyAKTChargedJets,fnAKTChargedJets);
    fTPCJetUnbiased->InitializeJetData(fmyAKTChargedJets,fnAKTChargedJets);
    
    // kT Jets
    fTPCkTFullJet = new AlipAJetData("fTPCkTFullJet",kFALSE,fnKTChargedJets);
    fTPCkTFullJet->SetSignalCut(fTPCJetThreshold);
    fTPCkTFullJet->SetAreaCutFraction(0.25*fJetAreaCutFrac);
    fTPCkTFullJet->SetJetR(fJetR);

    for (i=0;i<fnKTChargedJets;i++)
    {
        AliEmcalJet *myJet =(AliEmcalJet*) fmyKTChargedJets->At(i);
        fTPCkTFullJet->SetIsJetInArray(IsInTPC(fJetR,myJet->Phi(),myJet->Eta(),kTRUE),i);
    }
    fTPCkTFullJet->InitializeJetData(fmyKTChargedJets,fnKTChargedJets);

    // Raw Charged Jet Spectra
    if (fCalculateRhoJet>=2 || (fCalculateRhoJet>=1 && fMCPartLevel==kTRUE))
    {
        fTPCRawJets->FillBSJS(fEventCentrality,0.0,fTPCJetThreshold,fmyAKTChargedJets,fTPCFullJet->GetJets(),fTPCFullJet->GetTotalJets());
    }
}

void AliAnalysisTaskFullpAJets::InitFullJets()
{
    // Preliminary Jet Placement and Selection Cuts
    Int_t i;
    
    fnAKTFullJets = fmyAKTFullJets->GetEntries();
    fnKTFullJets = fmyKTFullJets->GetEntries();
    
    fEMCalJet = new AlipAJetData("fEMCalJet",kTRUE,fnAKTFullJets);
    fEMCalFullJet = new AlipAJetData("fEMCalFullJet",kTRUE,fnAKTFullJets);
    fEMCalPartJet = new AlipAJetData("fEMCalPartJet",kTRUE,fnAKTFullJets);
    fEMCalPartJetUnbiased = new AlipAJetData("fEMCalPartJetUnbiased",kTRUE,fnAKTFullJets);
    
    fEMCalJet->SetSignalCut(fEMCalJetThreshold);
    fEMCalJet->SetAreaCutFraction(fJetAreaCutFrac);
    fEMCalJet->SetJetR(fJetR);
    fEMCalJet->SetNEF(fNEFSignalJetCut);
    fEMCalJet->SetSignalTrackPtBias(fSignalTrackBias);
    fEMCalFullJet->SetSignalCut(fEMCalJetThreshold);
    fEMCalFullJet->SetAreaCutFraction(fJetAreaCutFrac);
    fEMCalFullJet->SetJetR(fJetR);
    fEMCalFullJet->SetNEF(fNEFSignalJetCut);
    fEMCalFullJet->SetSignalTrackPtBias(fSignalTrackBias);
    fEMCalPartJet->SetSignalCut(fEMCalJetThreshold);
    fEMCalPartJet->SetAreaCutFraction(fJetAreaCutFrac);
    fEMCalPartJet->SetJetR(fJetR);
    fEMCalPartJet->SetNEF(fNEFSignalJetCut);
    fEMCalPartJet->SetSignalTrackPtBias(fSignalTrackBias);
    fEMCalPartJetUnbiased->SetSignalCut(fEMCalJetThreshold);
    fEMCalPartJetUnbiased->SetAreaCutFraction(fJetAreaCutFrac);
    fEMCalPartJetUnbiased->SetJetR(fJetR);
    fEMCalPartJetUnbiased->SetNEF(1.0);
    fEMCalPartJetUnbiased->SetSignalTrackPtBias(kFALSE);
    
    // Initialize Jet Data
    for (i=0;i<fnAKTFullJets;i++)
    {
        AliEmcalJet *myJet =(AliEmcalJet*) fmyAKTFullJets->At(i);
        
        fEMCalJet->SetIsJetInArray(IsInEMCal(myJet->Phi(),myJet->Eta()),i);
        fEMCalFullJet->SetIsJetInArray(IsInEMCalFull(fJetRAccept,myJet->Phi(),myJet->Eta()),i);
        fEMCalPartJet->SetIsJetInArray(IsInEMCalPart(fJetR,myJet->Phi(),myJet->Eta()),i);
        fEMCalPartJetUnbiased->SetIsJetInArray(IsInEMCalPart(fJetR,myJet->Phi(),myJet->Eta()),i);
    }
    fEMCalJet->InitializeJetData(fmyAKTFullJets,fnAKTFullJets);
    fEMCalFullJet->InitializeJetData(fmyAKTFullJets,fnAKTFullJets);
    fEMCalPartJet->InitializeJetData(fmyAKTFullJets,fnAKTFullJets);
    fEMCalPartJetUnbiased->InitializeJetData(fmyAKTFullJets,fnAKTFullJets);

    // kT Jets
    fEMCalkTFullJet = new AlipAJetData("fEMCalkTFullJet",kTRUE,fnKTFullJets);
    fEMCalkTFullJet->SetSignalCut(fEMCalJetThreshold);
    fEMCalkTFullJet->SetAreaCutFraction(0.25*fJetAreaCutFrac);
    fEMCalkTFullJet->SetJetR(fJetR);
    fEMCalkTFullJet->SetNEF(fNEFSignalJetCut);
    fEMCalkTFullJet->SetSignalTrackPtBias(fSignalTrackBias);
    
    for (i=0;i<fnKTFullJets;i++)
    {
        AliEmcalJet *myJet =(AliEmcalJet*) fmyKTFullJets->At(i);
        fEMCalkTFullJet->SetIsJetInArray(IsInEMCalFull(fJetR,myJet->Phi(),myJet->Eta()),i);
    }
    fEMCalkTFullJet->InitializeJetData(fmyKTFullJets,fnKTFullJets);

    // Raw Full Jet Spectra
    if (fMCPartLevel==kFALSE)
    {
        fEMCalRawJets->FillBSJS(fEventCentrality,0.0,fEMCalJetThreshold,fmyAKTFullJets,fEMCalFullJet->GetJets(),fEMCalFullJet->GetTotalJets());
    }
}

void AliAnalysisTaskFullpAJets::GenerateTPCRandomConesPt()
{
    Int_t i,j;
    Double_t E_tracks_total=0.;
    TRandom3 u(time(NULL));
    
    Double_t Eta_Center=0.5*(fTPCEtaMin+fTPCEtaMax);
    Double_t Phi_Center=0.5*(fTPCPhiMin+fTPCPhiMax);
    Int_t event_mult=0;
    Int_t event_track_mult=0;
    Int_t clus_mult=0;
    
    for (i=0;i<fnBckgClusters;i++)
    {
        fTPCRCBckgFluc[i]=0.0;
        fTPCRCBckgFlucSignal[i]=0.0;
        fTPCRCBckgFlucNColl[i]=0.0;
    }
    
    TLorentzVector *dummy= new TLorentzVector;
    TLorentzVector *temp_jet= new TLorentzVector;
    TLorentzVector *track_vec = new TLorentzVector;

    // First, consider the RC with no spatial restrictions
    for (j=0;j<fnBckgClusters;j++)
    {
        E_tracks_total=0.;
        
        dummy->SetPtEtaPhiE(1,u.Uniform(Eta_Center-fJetR,Eta_Center+fJetR),u.Uniform(Phi_Center-fJetR,Phi_Center+fJetR),0);
        //  Loop over all tracks
        for (i=0;i<fnTracks;i++)
        {
            AliVTrack* vtrack =(AliVTrack*) fmyTracks->At(i);
            if (IsInTPC(fJetR,vtrack->Phi(),vtrack->Eta(),kFALSE)==kTRUE)
            {
                track_vec->SetPtEtaPhiE(vtrack->Pt(),vtrack->Eta(),vtrack->Phi(),vtrack->E());
                if (dummy->DeltaR(*track_vec)<fJetR)
                {
                    E_tracks_total+=vtrack->Pt();
                }
            }
        }
        fTPCRCBckgFlucSignal[j]=E_tracks_total;
    }
    
    // Now, consider the RC where the vertex of RC is at least 2R away from the leading signal
    E_tracks_total=0.0;
    if (fTPCJetUnbiased->GetLeadingPt()<0.0)
    {
        temp_jet->SetPtEtaPhiE(1,u.Uniform(Eta_Center-fJetR,Eta_Center+fJetR),u.Uniform(Phi_Center-fJetR,Phi_Center+fJetR),0);
    }
    else
    {
        AliEmcalJet *myJet = (AliEmcalJet*) fmyAKTChargedJets->At(fTPCJetUnbiased->GetLeadingIndex());
        myJet->GetMom(*temp_jet);
    }
    
    for (j=0;j<fnBckgClusters;j++)
    {
        event_mult=0;
        event_track_mult=0;
        clus_mult=0;
        E_tracks_total=0.;
        
        dummy->SetPtEtaPhiE(1,u.Uniform(Eta_Center-fJetR,Eta_Center+fJetR),u.Uniform(Phi_Center-fJetR,Phi_Center+fJetR),0);
        while (temp_jet->DeltaR(*dummy)<fJetR)
        {
            dummy->SetPtEtaPhiE(1,u.Uniform(Eta_Center-fJetR,Eta_Center+fJetR),u.Uniform(Phi_Center-fJetR,Phi_Center+fJetR),0);
        }
        //  Loop over all tracks
        for (i=0;i<fnTracks;i++)
        {
            AliVTrack* vtrack =(AliVTrack*) fmyTracks->At(i);
            if (IsInTPC(fJetR,vtrack->Phi(),vtrack->Eta(),kFALSE) == kTRUE)
            {
                event_mult++;
                if (IsInEMCal(vtrack->Phi(),vtrack->Eta()) == kTRUE)
                {
                    event_track_mult++;
                }
                track_vec->SetPtEtaPhiE(vtrack->Pt(),vtrack->Eta(),vtrack->Phi(),vtrack->E());
                if (dummy->DeltaR(*track_vec)<fJetR)
                {
                    clus_mult++;
                    E_tracks_total+=vtrack->Pt();
                }
            }
        }
        fTPCRCBckgFluc[j]=E_tracks_total;
    }
    if (fTrackQA==kTRUE)
    {
        fhTPCEventMult->Fill(fEventCentrality,event_mult);
        fpTPCEventMult->Fill(fEventCentrality,event_mult);
        fhEMCalTrackEventMult->Fill(fEventCentrality,event_track_mult);
    }
    if (fCalculateRhoJet>=2 || (fCalculateRhoJet>=1 && fMCPartLevel==kTRUE))
    {
        fTPCRawJets->FillDeltaPt(fEventCentrality,0.0,fJetR,fTPCRCBckgFluc,1);
    }
    
    // For the case of partial exclusion, merely allow a superposition of full and no exclusion with probability p=1/Ncoll
    Double_t exclusion_prob;
    for (j=0;j<fnBckgClusters;j++)
    {
        exclusion_prob = u.Uniform(0,1);
        if (exclusion_prob<(1/fNColl))
        {
            fTPCRCBckgFlucNColl[j]=fTPCRCBckgFlucSignal[j];
        }
        else
        {
            fTPCRCBckgFlucNColl[j]=fTPCRCBckgFluc[j];
        }
    }
    
    delete dummy;
    delete temp_jet;
    delete track_vec;
}

void AliAnalysisTaskFullpAJets::GenerateEMCalRandomConesPt()
{
    Int_t i,j;
    Double_t E_tracks_total=0.;
    Double_t E_caloclusters_total=0.;
    TRandom3 u(time(NULL));
    
    Double_t Eta_Center=0.5*(fEMCalEtaMin+fEMCalEtaMax);
    Double_t Phi_Center=0.5*(fEMCalPhiMin+fEMCalPhiMax);
    Int_t event_mult=0;
    Int_t clus_mult=0;
    
    for (i=0;i<fnBckgClusters;i++)
    {
        fEMCalRCBckgFluc[i]=0.0;
        fEMCalRCBckgFlucSignal[i]=0.0;
        fEMCalRCBckgFlucNColl[i]=0.0;
    }
    
    TLorentzVector *dummy= new TLorentzVector;
    TLorentzVector *temp_jet= new TLorentzVector;
    TLorentzVector *track_vec = new TLorentzVector;
    TLorentzVector *cluster_vec = new TLorentzVector;

    // First, consider the RC with no spatial restrictions
    for (j=0;j<fnBckgClusters;j++)
    {
        E_tracks_total=0.;
        E_caloclusters_total=0.;
        
        dummy->SetPtEtaPhiE(1,u.Uniform(Eta_Center-fJetR,Eta_Center+fJetR),u.Uniform(Phi_Center-fJetR,Phi_Center+fJetR),0);
        //  Loop over all tracks
        for (i=0;i<fnTracks;i++)
        {
            AliVTrack* vtrack =(AliVTrack*) fmyTracks->At(i);
            if (IsInEMCal(vtrack->Phi(),vtrack->Eta())==kTRUE)
            {
                track_vec->SetPtEtaPhiE(vtrack->Pt(),vtrack->Eta(),vtrack->Phi(),vtrack->E());
                if (dummy->DeltaR(*track_vec)<fJetR)
                {
                    E_tracks_total+=vtrack->Pt();
                }
            }
        }
        
        //  Loop over all caloclusters
        for (i=0;i<fnClusters;i++)
        {
            AliVCluster* vcluster = (AliVCluster*) fmyClusters->At(i);
            vcluster->GetMomentum(*cluster_vec,fVertex);
            if (dummy->DeltaR(*cluster_vec)<fJetR)
            {
                clus_mult++;
                E_caloclusters_total+=cluster_vec->Pt();
            }
        }
        fEMCalRCBckgFlucSignal[j]=E_tracks_total+E_caloclusters_total;
    }

    // Now, consider the RC where the vertex of RC is at least 2R away from the leading signal
    E_tracks_total=0.;
    E_caloclusters_total=0.;
    if (fEMCalPartJetUnbiased->GetLeadingPt()<0.0)
    {
        temp_jet->SetPtEtaPhiE(1,u.Uniform(Eta_Center-fJetR,Eta_Center+fJetR),u.Uniform(Phi_Center-fJetR,Phi_Center+fJetR),0);
    }
    else
    {
        AliEmcalJet *myJet = (AliEmcalJet*) fmyAKTFullJets->At(fEMCalPartJetUnbiased->GetLeadingIndex());
        myJet->GetMom(*temp_jet);
    }
    
    for (j=0;j<fnBckgClusters;j++)
    {
        event_mult=0;
        clus_mult=0;
        E_tracks_total=0.;
        E_caloclusters_total=0.;
        
        dummy->SetPtEtaPhiE(1,u.Uniform(Eta_Center-fJetR,Eta_Center+fJetR),u.Uniform(Phi_Center-fJetR,Phi_Center+fJetR),0);
        while (temp_jet->DeltaR(*dummy)<fJetR)
        {
            dummy->SetPtEtaPhiE(1,u.Uniform(Eta_Center-fJetR,Eta_Center+fJetR),u.Uniform(Phi_Center-fJetR,Phi_Center+fJetR),0);
        }
        //  Loop over all tracks
        for (i=0;i<fnTracks;i++)
        {
            AliVTrack* vtrack =(AliVTrack*) fmyTracks->At(i);
            if (IsInEMCal(vtrack->Phi(),vtrack->Eta())==kTRUE)
            {
                event_mult++;
                track_vec->SetPtEtaPhiE(vtrack->Pt(),vtrack->Eta(),vtrack->Phi(),vtrack->E());
                if (dummy->DeltaR(*track_vec)<fJetR)
                {
                    clus_mult++;
                    E_tracks_total+=vtrack->Pt();
                }
            }
        }
        
        //  Loop over all caloclusters
        for (i=0;i<fnClusters;i++)
        {
            AliVCluster* vcluster = (AliVCluster*) fmyClusters->At(i);
            vcluster->GetMomentum(*cluster_vec,fVertex);
            event_mult++;
            if (dummy->DeltaR(*cluster_vec)<fJetR)
            {
                clus_mult++;
                E_caloclusters_total+=cluster_vec->Pt();
            }
        }
        fEMCalRCBckgFluc[j]=E_tracks_total+E_caloclusters_total;
    }
    if (fClusterQA==kTRUE)
    {
        fhEMCalEventMult->Fill(fEventCentrality,event_mult);
        fpEMCalEventMult->Fill(fEventCentrality,event_mult);
    }
    fEMCalRawJets->FillDeltaPt(fEventCentrality,0.0,fJetR,fEMCalRCBckgFluc,1);
    
    // For the case of partial exclusion, merely allow a superposition of full and no exclusion with probability p=1/Ncoll
    Double_t exclusion_prob;
    for (j=0;j<fnBckgClusters;j++)
    {
        exclusion_prob = u.Uniform(0,1);
        if (exclusion_prob<(1/fNColl))
        {
            fEMCalRCBckgFlucNColl[j]=fEMCalRCBckgFlucSignal[j];
        }
        else
        {
            fEMCalRCBckgFlucNColl[j]=fEMCalRCBckgFluc[j];
        }
    }

    delete dummy;
    delete temp_jet;
    delete track_vec;
    delete cluster_vec;
}

// Charged Rho's
void AliAnalysisTaskFullpAJets::EstimateChargedRho0()
{
    Int_t i;
    Double_t E_tracks_total=0.0;
    Double_t TPC_rho=0.0;
    
    //  Loop over all tracks
    for (i=0;i<fnTracks;i++)
    {
        AliVTrack* vtrack =(AliVTrack*) fmyTracks->At(i);
        if (IsInTPC(fJetR,vtrack->Phi(),vtrack->Eta(),kFALSE)==kTRUE)
        {
            E_tracks_total+=vtrack->Pt();
        }
    }
    
    //  Calculate the mean Background density
    TPC_rho=E_tracks_total/fTPCArea;
    
    fRhoCharged = TPC_rho;
    
    // Fill Histograms
    fRhoCharged0->FillRho(fEventCentrality,TPC_rho);
    fRhoCharged0->FillBSJS(fEventCentrality,TPC_rho,fTPCJetThreshold,fmyAKTChargedJets,fTPCJet->GetJets(),fTPCJet->GetTotalJets());
    fRhoCharged0->FillDeltaPt(fEventCentrality,TPC_rho,fJetR,fTPCRCBckgFluc,1);
    fRhoCharged0->FillDeltaPtSignal(fEventCentrality,TPC_rho,fJetR,fTPCRCBckgFlucSignal,1);
    fRhoCharged0->FillDeltaPtNColl(fEventCentrality,TPC_rho,fJetR,fTPCRCBckgFlucNColl,1);
    fRhoCharged0->FillBackgroundFluctuations(fEventCentrality,TPC_rho,fJetR);
    fRhoCharged0->FillLeadingJetPtRho(fTPCJet->GetLeadingPt(),TPC_rho);
    
}

void AliAnalysisTaskFullpAJets::EstimateChargedRho1()
{
    Int_t i;
    Double_t E_tracks_total=0.0;
    Double_t TPC_rho=0.;
    
    TLorentzVector *temp_jet= new TLorentzVector;
    TLorentzVector *track_vec = new TLorentzVector;

    if (fTPCJetUnbiased->GetLeadingPt()>=fTPCJetThreshold)
    {
        AliEmcalJet *myJet =(AliEmcalJet*) fmyAKTChargedJets->At(fTPCJet->GetLeadingIndex());
        myJet->GetMom(*temp_jet);
        
        //  Loop over all tracks
        for (i=0;i<fnTracks;i++)
        {
            AliVTrack* vtrack =(AliVTrack*) fmyTracks->At(i);
            if (IsInTPC(fJetR,vtrack->Phi(),vtrack->Eta(),kFALSE)==kTRUE)
            {
                track_vec->SetPtEtaPhiE(vtrack->Pt(),vtrack->Eta(),vtrack->Phi(),vtrack->E());
                if (temp_jet->DeltaR(*track_vec)>fJetRForRho)
                {
                    E_tracks_total+=vtrack->Pt();
                }
            }
        }
        
        //  Calculate the mean Background density
        TPC_rho=E_tracks_total/(fTPCArea-AreaWithinTPC(fJetR,myJet->Eta()));
    }
    else  // i.e. No signal jets -> same as total background density
    {
        //  Loop over all tracks
        for (i=0;i<fnTracks;i++)
        {
            AliVTrack* vtrack =(AliVTrack*) fmyTracks->At(i);
            if (IsInTPC(fJetR,vtrack->Phi(),vtrack->Eta(),kFALSE)==kTRUE)
            {
                E_tracks_total+=vtrack->Pt();
            }
        }
        //  Calculate the mean Background density
        TPC_rho=E_tracks_total/fTPCArea;
    }
    delete track_vec;
    delete temp_jet;

    // Fill histograms
    fRhoCharged1->FillRho(fEventCentrality,TPC_rho);
    fRhoCharged1->FillBSJS(fEventCentrality,TPC_rho,fTPCJetThreshold,fmyAKTChargedJets,fTPCFullJet->GetJets(),fTPCFullJet->GetTotalJets());
    fRhoCharged1->FillDeltaPt(fEventCentrality,TPC_rho,fJetR,fTPCRCBckgFluc,1);
    fRhoCharged1->FillDeltaPtSignal(fEventCentrality,TPC_rho,fJetR,fTPCRCBckgFlucSignal,1);
    fRhoCharged1->FillDeltaPtNColl(fEventCentrality,TPC_rho,fJetR,fTPCRCBckgFlucNColl,1);
    fRhoCharged1->FillBackgroundFluctuations(fEventCentrality,TPC_rho,fJetR);
    fRhoCharged1->FillLeadingJetPtRho(fTPCFullJet->GetLeadingPt(),TPC_rho);
}

void AliAnalysisTaskFullpAJets::EstimateChargedRho2()
{
    Int_t i;
    Double_t E_tracks_total=0.0;
    Double_t TPC_rho=0.;

    TLorentzVector *temp_jet1= new TLorentzVector;
    TLorentzVector *temp_jet2= new TLorentzVector;
    TLorentzVector *track_vec = new TLorentzVector;

    if ((fTPCJetUnbiased->GetLeadingPt()>=fTPCJetThreshold) && (fTPCJetUnbiased->GetSubLeadingPt()>=fTPCJetThreshold))
    {
        AliEmcalJet *myhJet =(AliEmcalJet*) fmyAKTChargedJets->At(fTPCJet->GetLeadingIndex());
        myhJet->GetMom(*temp_jet1);

        AliEmcalJet *myJet =(AliEmcalJet*) fmyAKTChargedJets->At(fTPCJet->GetSubLeadingIndex());
        myJet->GetMom(*temp_jet2);

        //  Loop over all tracks
        for (i=0;i<fnTracks;i++)
        {
            AliVTrack* vtrack =(AliVTrack*) fmyTracks->At(i);
            if (IsInTPC(fJetR,vtrack->Phi(),vtrack->Eta(),kFALSE)==kTRUE)
            {
                track_vec->SetPtEtaPhiE(vtrack->Pt(),vtrack->Eta(),vtrack->Phi(),vtrack->E());
                if ((temp_jet1->DeltaR(*track_vec)>fJetRForRho) && (temp_jet2->DeltaR(*track_vec)>fJetRForRho))
                {
                    E_tracks_total+=vtrack->Pt();
                }
            }
        }
        
        //  Calculate the mean Background density
        TPC_rho=E_tracks_total/(fTPCArea-AreaWithinTPC(fJetR,myhJet->Eta())-AreaWithinTPC(fJetR,myJet->Eta()));
    }
    else if (fTPCJetUnbiased->GetLeadingPt()>=fTPCJetThreshold)
    {
        AliEmcalJet *myJet =(AliEmcalJet*) fmyAKTChargedJets->At(fTPCJet->GetLeadingIndex());
        myJet->GetMom(*temp_jet1);
        
        //  Loop over all tracks
        for (i=0;i<fnTracks;i++)
        {
            AliVTrack* vtrack =(AliVTrack*) fmyTracks->At(i);
            if (IsInTPC(fJetR,vtrack->Phi(),vtrack->Eta(),kFALSE)==kTRUE)
            {
                track_vec->SetPtEtaPhiE(vtrack->Pt(),vtrack->Eta(),vtrack->Phi(),vtrack->E());
                if (temp_jet1->DeltaR(*track_vec)>fJetRForRho)
                {
                    E_tracks_total+=vtrack->Pt();
                }
            }
        }
        
        //  Calculate the mean Background density
        TPC_rho=E_tracks_total/(fTPCArea-AreaWithinTPC(fJetR,myJet->Eta()));
    }
    else  // i.e. No signal jets -> same as total background density
    {
        //  Loop over all tracks
        for (i=0;i<fnTracks;i++)
        {
            AliVTrack* vtrack =(AliVTrack*) fmyTracks->At(i);
            if (IsInTPC(fJetR,vtrack->Phi(),vtrack->Eta(),kFALSE)==kTRUE)
            {
                E_tracks_total+=vtrack->Pt();
            }
        }
        
        //  Calculate the mean Background density
        TPC_rho=E_tracks_total/fTPCArea;
    }
    delete temp_jet1;
    delete temp_jet2;
    delete track_vec;

    // Fill histograms
    fRhoCharged2->FillRho(fEventCentrality,TPC_rho);
    fRhoCharged2->FillBSJS(fEventCentrality,TPC_rho,fTPCJetThreshold,fmyAKTChargedJets,fTPCFullJet->GetJets(),fTPCFullJet->GetTotalJets());
    fRhoCharged2->FillDeltaPt(fEventCentrality,TPC_rho,fJetR,fTPCRCBckgFluc,1);
    fRhoCharged2->FillDeltaPtSignal(fEventCentrality,TPC_rho,fJetR,fTPCRCBckgFlucSignal,1);
    fRhoCharged2->FillDeltaPtNColl(fEventCentrality,TPC_rho,fJetR,fTPCRCBckgFlucNColl,1);
    fRhoCharged2->FillBackgroundFluctuations(fEventCentrality,TPC_rho,fJetR);
    fRhoCharged2->FillLeadingJetPtRho(fTPCFullJet->GetLeadingPt(),TPC_rho);
}

void AliAnalysisTaskFullpAJets::EstimateChargedRhoN()
{
    Int_t i,j;
    Bool_t track_away_from_jet;
    Double_t E_tracks_total=0.0;
    Double_t TPC_rho=0.0;
    Double_t jet_area_total=0.0;
    
    TLorentzVector *jet_vec= new TLorentzVector;
    TLorentzVector *track_vec = new TLorentzVector;

    // First, sum all tracks within the EMCal that are away from jet(s) above Pt Threshold
    for (i=0;i<fnTracks;i++)
    {
        // First, check if track is in the EMCal!!
        AliVTrack* vtrack = (AliVTrack*) fmyTracks->At(i);
        if (IsInTPC(fJetR,vtrack->Phi(),vtrack->Eta(),kFALSE)==kTRUE)
        {
            if (fTPCJetUnbiased->GetTotalSignalJets()<1)
            {
                E_tracks_total+=vtrack->Pt();
            }
            else
            {
                track_away_from_jet=kTRUE;
                j=0;
                track_vec->SetPtEtaPhiE(vtrack->Pt(),vtrack->Eta(),vtrack->Phi(),vtrack->E());
                while (track_away_from_jet==kTRUE && j<fTPCJetUnbiased->GetTotalSignalJets())
                {
                    AliEmcalJet *myJet =(AliEmcalJet*) fmyAKTChargedJets->At(fTPCJetUnbiased->GetSignalJetIndex(j));
                    myJet->GetMom(*jet_vec);
                    if (track_vec->DeltaR(*jet_vec)<=fJetRForRho)
                    {
                        track_away_from_jet=kFALSE;
                    }
                    j++;
                }
                if (track_away_from_jet==kTRUE)
                {
                    E_tracks_total+=vtrack->Pt();
                }
            }
        }
    }
    
    // Determine area of all Jets that are within the EMCal
    if (fTPCJetUnbiased->GetTotalSignalJets()==0)
    {
        jet_area_total=0.0;
    }
    else
    {
        for (i=0;i<fTPCJetUnbiased->GetTotalSignalJets();i++)
        {
            AliEmcalJet *myJet = (AliEmcalJet*) fmyAKTChargedJets->At(fTPCJetUnbiased->GetSignalJetIndex(i));
            jet_area_total+=AreaWithinTPC(fJetR,myJet->Eta());
        }
    }
    delete jet_vec;
    delete track_vec;

    // Calculate Rho
    TPC_rho = E_tracks_total/(fTPCArea-jet_area_total);
    
    // Fill Histogram
    fRhoChargedN->FillRho(fEventCentrality,TPC_rho);
    fRhoChargedN->FillBSJS(fEventCentrality,TPC_rho,fTPCJetThreshold,fmyAKTChargedJets,fTPCFullJet->GetJets(),fTPCFullJet->GetTotalJets());
    fRhoChargedN->FillDeltaPt(fEventCentrality,TPC_rho,fJetR,fTPCRCBckgFluc,1);
    fRhoChargedN->FillDeltaPtSignal(fEventCentrality,TPC_rho,fJetR,fTPCRCBckgFlucSignal,1);
    fRhoChargedN->FillDeltaPtNColl(fEventCentrality,TPC_rho,fJetR,fTPCRCBckgFlucNColl,1);
    fRhoChargedN->FillBackgroundFluctuations(fEventCentrality,TPC_rho,fJetR);
    fRhoChargedN->FillLeadingJetPtRho(fTPCFullJet->GetLeadingPt(),TPC_rho);

}

void AliAnalysisTaskFullpAJets::EstimateChargedRhoScale()
{
    Int_t i,j;
    Bool_t track_away_from_jet;
    Double_t E_tracks_total=0.0;
    Double_t TPC_rho=0.0;
    Double_t jet_area_total=0.0;
    
    TLorentzVector *jet_vec= new TLorentzVector;
    TLorentzVector *track_vec = new TLorentzVector;

    // First, sum all tracks within the EMCal that are away from jet(s) above Pt Threshold
    for (i=0;i<fnTracks;i++)
    {
        // First, check if track is in the EMCal!!
        AliVTrack* vtrack = (AliVTrack*) fmyTracks->At(i);
        if (IsInTPC(fJetR,vtrack->Phi(),vtrack->Eta(),kFALSE)==kTRUE)
        {
            if (fTPCJetUnbiased->GetTotalSignalJets()<1)
            {
                E_tracks_total+=vtrack->Pt();
            }
            else
            {
                track_away_from_jet=kTRUE;
                j=0;
                track_vec->SetPtEtaPhiE(vtrack->Pt(),vtrack->Eta(),vtrack->Phi(),vtrack->E());
                while (track_away_from_jet==kTRUE && j<fTPCJetUnbiased->GetTotalSignalJets())
                {
                    AliEmcalJet *myJet =(AliEmcalJet*) fmyAKTChargedJets->At(fTPCJetUnbiased->GetSignalJetIndex(j));
                    myJet->GetMom(*jet_vec);
                    if (track_vec->DeltaR(*jet_vec)<=fJetRForRho)
                    {
                        track_away_from_jet=kFALSE;
                    }
                    j++;
                }
                if (track_away_from_jet==kTRUE)
                {
                    E_tracks_total+=vtrack->Pt();
                }
            }
        }
    }
    
    // Determine area of all Jets that are within the TPC
    if (fTPCJetUnbiased->GetTotalSignalJets()==0)
    {
        jet_area_total=0.0;
    }
    else
    {
        for (i=0;i<fTPCJetUnbiased->GetTotalSignalJets();i++)
        {
            AliEmcalJet *myJet = (AliEmcalJet*) fmyAKTChargedJets->At(fTPCJetUnbiased->GetSignalJetIndex(i));
            jet_area_total+=AreaWithinTPC(fJetR,myJet->Eta());
        }
    }
    delete jet_vec;
    delete track_vec;

    // Calculate Rho
    TPC_rho = E_tracks_total/(fTPCArea-jet_area_total);
    TPC_rho*=fScaleFactor;
    
    // Fill Histogram
    fRhoChargedScale->FillRho(fEventCentrality,TPC_rho);
    fRhoChargedScale->FillBSJS(fEventCentrality,TPC_rho,fEMCalJetThreshold,fmyAKTFullJets,fEMCalFullJet->GetJets(),fEMCalFullJet->GetTotalJets());
    fRhoChargedScale->FillDeltaPt(fEventCentrality,TPC_rho,fJetR,fEMCalRCBckgFluc,1);
    fRhoChargedScale->FillDeltaPtSignal(fEventCentrality,TPC_rho,fJetR,fEMCalRCBckgFlucSignal,1);
    fRhoChargedScale->FillDeltaPtNColl(fEventCentrality,TPC_rho,fJetR,fEMCalRCBckgFlucNColl,1);
    fRhoChargedScale->FillBackgroundFluctuations(fEventCentrality,TPC_rho,fJetR);
    fRhoChargedScale->FillLeadingJetPtRho(fEMCalFullJet->GetLeadingPt(),TPC_rho);
    fRhoChargedScale->FillMiscJetStats(fmyAKTFullJets,fEMCalFullJet->GetJets(),fEMCalFullJet->GetTotalJets(),fOrgTracks,fOrgClusters,fVertex);
}

void AliAnalysisTaskFullpAJets::EstimateChargedRhokT()
{
    Int_t i;
    Double_t kTRho = 0.0;
    Double_t *pTArray = new Double_t[fTPCkTFullJet->GetTotalJets()];
    Double_t *RhoArray = new Double_t[fTPCkTFullJet->GetTotalJets()];
    
    for (i=0;i<fTPCkTFullJet->GetTotalJets();i++)
    {
        AliEmcalJet *myJet =(AliEmcalJet*) fmyKTChargedJets->At(fTPCkTFullJet->GetJetIndex(i));
        pTArray[i]=myJet->Pt();
        RhoArray[i]=myJet->Pt()/myJet->Area();
    }
    
    if (fTPCkTFullJet->GetTotalJets()>=2)
    {
        kTRho=MedianRhokT(pTArray,RhoArray,fTPCkTFullJet->GetTotalJets());
        
        fRhoChargedkT->FillRho(fEventCentrality,kTRho);
        fRhoChargedkT->FillBSJS(fEventCentrality,kTRho,fTPCJetThreshold,fmyAKTChargedJets,fTPCFullJet->GetJets(),fTPCFullJet->GetTotalJets());
        fRhoChargedkT->FillDeltaPt(fEventCentrality,kTRho,fJetR,fTPCRCBckgFluc,1);
        fRhoChargedkT->FillDeltaPtSignal(fEventCentrality,kTRho,fJetR,fTPCRCBckgFlucSignal,1);
        fRhoChargedkT->FillDeltaPtNColl(fEventCentrality,kTRho,fJetR,fTPCRCBckgFlucNColl,1);
        fRhoChargedkT->FillBackgroundFluctuations(fEventCentrality,kTRho,fJetR);
        fRhoChargedkT->FillLeadingJetPtRho(fTPCFullJet->GetLeadingPt(),kTRho);
    }
    delete [] RhoArray;
    delete [] pTArray;
}

void AliAnalysisTaskFullpAJets::EstimateChargedRhokTScale()
{
    Int_t i;
    Double_t kTRho = 0.0;
    Double_t *pTArray = new Double_t[fTPCkTFullJet->GetTotalJets()];
    Double_t *RhoArray = new Double_t[fTPCkTFullJet->GetTotalJets()];
    
    for (i=0;i<fTPCkTFullJet->GetTotalJets();i++)
    {
        AliEmcalJet *myJet =(AliEmcalJet*) fmyKTChargedJets->At(fTPCkTFullJet->GetJetIndex(i));
        pTArray[i]=myJet->Pt();
        RhoArray[i]=myJet->Pt()/myJet->Area();
    }
    
    if (fTPCkTFullJet->GetTotalJets()>=2)
    {
        kTRho=MedianRhokT(pTArray,RhoArray,fTPCkTFullJet->GetTotalJets());
        kTRho*=fScaleFactor;
        
        fRhoChargedkTScale->FillRho(fEventCentrality,kTRho);
        fRhoChargedkTScale->FillBSJS(fEventCentrality,kTRho,fEMCalJetThreshold,fmyAKTFullJets,fEMCalFullJet->GetJets(),fEMCalFullJet->GetTotalJets());
        fRhoChargedkTScale->FillDeltaPt(fEventCentrality,kTRho,fJetR,fEMCalRCBckgFluc,1);
        fRhoChargedkTScale->FillDeltaPtSignal(fEventCentrality,kTRho,fJetR,fEMCalRCBckgFlucSignal,1);
        fRhoChargedkTScale->FillDeltaPtNColl(fEventCentrality,kTRho,fJetR,fEMCalRCBckgFlucNColl,1);
        fRhoChargedkTScale->FillBackgroundFluctuations(fEventCentrality,kTRho,fJetR);
        fRhoChargedkTScale->FillLeadingJetPtRho(fEMCalFullJet->GetLeadingPt(),kTRho);
    }
    delete [] RhoArray;
    delete [] pTArray;
}

void AliAnalysisTaskFullpAJets::EstimateChargedRhoCMS()
{
    Int_t i,k;
    Double_t kTRho = 0.0;
    Double_t CMSTotalkTArea = 0.0;
    Double_t CMSTrackArea = 0.0;
    Double_t CMSCorrectionFactor = 1.0;
    Double_t *pTArray = new Double_t[fTPCkTFullJet->GetTotalJets()];
    Double_t *RhoArray = new Double_t[fTPCkTFullJet->GetTotalJets()];

    k=0;
    if ((fTPCJet->GetLeadingPt()>=fTPCJetThreshold) && (fTPCJet->GetSubLeadingPt()>=fTPCJetThreshold))
    {
        AliEmcalJet *myJet1 =(AliEmcalJet*) fmyAKTChargedJets->At(fTPCJet->GetLeadingIndex());
        AliEmcalJet *myJet2 =(AliEmcalJet*) fmyAKTChargedJets->At(fTPCJet->GetSubLeadingIndex());
        
        for (i=0;i<fTPCkTFullJet->GetTotalJets();i++)
        {
            AliEmcalJet *myJet =(AliEmcalJet*) fmyKTChargedJets->At(fTPCkTFullJet->GetJetIndex(i));
            
            CMSTotalkTArea+=myJet->Area();
            if (myJet->GetNumberOfTracks()>0)
            {
                CMSTrackArea+=myJet->Area();
            }
            if (IsJetOverlap(myJet,myJet1,kFALSE)==kFALSE && IsJetOverlap(myJet,myJet2,kFALSE)==kFALSE)
            {
                pTArray[k]=myJet->Pt();
                RhoArray[k]=myJet->Pt()/myJet->Area();
                k++;
            }
        }
        if (k>0)
        {
            kTRho=MedianRhokT(pTArray,RhoArray,k);
        }
        else
        {
            kTRho=0.0;
        }
    }
    else if (fTPCJet->GetLeadingPt()>=fTPCJetThreshold)
    {
        AliEmcalJet *myJet1 =(AliEmcalJet*) fmyAKTChargedJets->At(fTPCJet->GetLeadingIndex());
        
        for (i=0;i<fTPCkTFullJet->GetTotalJets();i++)
        {
            AliEmcalJet *myJet =(AliEmcalJet*) fmyKTChargedJets->At(fTPCkTFullJet->GetJetIndex(i));
            
            CMSTotalkTArea+=myJet->Area();
            if (myJet->GetNumberOfTracks()>0)
            {
                CMSTrackArea+=myJet->Area();
            }
            if (IsJetOverlap(myJet,myJet1,kFALSE)==kFALSE)
            {
                pTArray[k]=myJet->Pt();
                RhoArray[k]=myJet->Pt()/myJet->Area();
                k++;
            }
        }
        if (k>0)
        {
            kTRho=MedianRhokT(pTArray,RhoArray,k);
        }
        else
        {
            kTRho=0.0;
        }
    }
    else
    {
        for (i=0;i<fTPCkTFullJet->GetTotalJets();i++)
        {
            AliEmcalJet *myJet =(AliEmcalJet*) fmyKTChargedJets->At(fTPCkTFullJet->GetJetIndex(i));
            
            CMSTotalkTArea+=myJet->Area();
            if (myJet->GetNumberOfTracks()>0)
            {
                CMSTrackArea+=myJet->Area();
            }
            pTArray[k]=myJet->Pt();
            RhoArray[k]=myJet->Pt()/myJet->Area();
            k++;
        }
        if (k>0)
        {
            kTRho=MedianRhokT(pTArray,RhoArray,k);
        }
        else
        {
            kTRho=0.0;
        }
    }
    // Scale CMS Rho by Correction factor
    if (CMSTotalkTArea==0.0)
    {
        CMSCorrectionFactor = 1.0;
    }
    else
    {
        //CMSCorrectionFactor = CMSTrackArea/CMSTotalkTArea;
        CMSCorrectionFactor = CMSTrackArea/(fTPCPhiTotal*(fTPCEtaTotal-2*fJetR));  // The total physical area should be reduced by the eta cut due to looping over only fully contained kT jets within the TPC
    }
    kTRho*=CMSCorrectionFactor;
    fRhoChargedCMS->FillRho(fEventCentrality,kTRho);
    fRhoChargedCMS->FillBSJS(fEventCentrality,kTRho,fTPCJetThreshold,fmyAKTChargedJets,fTPCFullJet->GetJets(),fTPCFullJet->GetTotalJets());
    fRhoChargedCMS->FillDeltaPt(fEventCentrality,kTRho,fJetR,fTPCRCBckgFluc,1);
    fRhoChargedCMS->FillDeltaPtSignal(fEventCentrality,kTRho,fJetR,fTPCRCBckgFlucSignal,1);
    fRhoChargedCMS->FillDeltaPtNColl(fEventCentrality,kTRho,fJetR,fTPCRCBckgFlucNColl,1);
    fRhoChargedCMS->FillBackgroundFluctuations(fEventCentrality,kTRho,fJetR);
    fRhoChargedCMS->FillLeadingJetPtRho(fTPCFullJet->GetLeadingPt(),kTRho);
    delete [] RhoArray;
    delete [] pTArray;
}

void AliAnalysisTaskFullpAJets::EstimateChargedRhoCMSScale()
{
    Int_t i,k;
    Double_t kTRho = 0.0;
    Double_t CMSTotalkTArea = 0.0;
    Double_t CMSTrackArea = 0.0;
    Double_t CMSCorrectionFactor = 1.0;
    Double_t *pTArray = new Double_t[fTPCkTFullJet->GetTotalJets()];
    Double_t *RhoArray = new Double_t[fTPCkTFullJet->GetTotalJets()];
    
    k=0;
    if ((fTPCJet->GetLeadingPt()>=fTPCJetThreshold) && (fTPCJet->GetSubLeadingPt()>=fTPCJetThreshold))
    {
        AliEmcalJet *myJet1 =(AliEmcalJet*) fmyAKTChargedJets->At(fTPCJet->GetLeadingIndex());
        AliEmcalJet *myJet2 =(AliEmcalJet*) fmyAKTChargedJets->At(fTPCJet->GetSubLeadingIndex());
        
        for (i=0;i<fTPCkTFullJet->GetTotalJets();i++)
        {
            AliEmcalJet *myJet =(AliEmcalJet*) fmyKTChargedJets->At(fTPCkTFullJet->GetJetIndex(i));
            
            CMSTotalkTArea+=myJet->Area();
            if (myJet->GetNumberOfTracks()>0)
            {
                CMSTrackArea+=myJet->Area();
            }
            if (IsJetOverlap(myJet,myJet1,kFALSE)==kFALSE && IsJetOverlap(myJet,myJet2,kFALSE)==kFALSE)
            {
                pTArray[k]=myJet->Pt();
                RhoArray[k]=myJet->Pt()/myJet->Area();
                k++;
            }
        }
        if (k>0)
        {
            kTRho=MedianRhokT(pTArray,RhoArray,k);
        }
        else
        {
            kTRho=0.0;
        }
    }
    else if (fTPCJet->GetLeadingPt()>=fTPCJetThreshold)
    {
        AliEmcalJet *myJet1 =(AliEmcalJet*) fmyAKTChargedJets->At(fTPCJet->GetLeadingIndex());
        
        for (i=0;i<fTPCkTFullJet->GetTotalJets();i++)
        {
            AliEmcalJet *myJet =(AliEmcalJet*) fmyKTChargedJets->At(fTPCkTFullJet->GetJetIndex(i));
            
            CMSTotalkTArea+=myJet->Area();
            if (myJet->GetNumberOfTracks()>0)
            {
                CMSTrackArea+=myJet->Area();
            }
            if (IsJetOverlap(myJet,myJet1,kFALSE)==kFALSE)
            {
                pTArray[k]=myJet->Pt();
                RhoArray[k]=myJet->Pt()/myJet->Area();
                k++;
            }
        }
        if (k>0)
        {
            kTRho=MedianRhokT(pTArray,RhoArray,k);
        }
        else
        {
            kTRho=0.0;
        }
    }
    else
    {
        for (i=0;i<fTPCkTFullJet->GetTotalJets();i++)
        {
            AliEmcalJet *myJet =(AliEmcalJet*) fmyKTChargedJets->At(fTPCkTFullJet->GetJetIndex(i));
            
            CMSTotalkTArea+=myJet->Area();
            if (myJet->GetNumberOfTracks()>0)
            {
                CMSTrackArea+=myJet->Area();
            }
            pTArray[k]=myJet->Pt();
            RhoArray[k]=myJet->Pt()/myJet->Area();
            k++;
        }
        if (k>0)
        {
            kTRho=MedianRhokT(pTArray,RhoArray,k);
        }
        else
        {
            kTRho=0.0;
        }
    }
    kTRho*=fScaleFactor;
    // Scale CMS Rho by Correction factor
    if (CMSTotalkTArea==0.0)
    {
        CMSCorrectionFactor = 1.0;
    }
    else
    {
        CMSCorrectionFactor = CMSTrackArea/(fTPCPhiTotal*(fTPCEtaTotal-2*fJetR));  // The total physical area should be reduced by the eta cut due to looping over only fully contained kT jets within the TPC
    }
    kTRho*=CMSCorrectionFactor;
    
    fRhoChargedCMSScale->FillRho(fEventCentrality,kTRho);
    fRhoChargedCMSScale->FillBSJS(fEventCentrality,kTRho,fEMCalJetThreshold,fmyAKTFullJets,fEMCalFullJet->GetJets(),fEMCalFullJet->GetTotalJets());
    fRhoChargedCMSScale->FillDeltaPt(fEventCentrality,kTRho,fJetR,fEMCalRCBckgFluc,1);
    fRhoChargedCMSScale->FillDeltaPtSignal(fEventCentrality,kTRho,fJetR,fEMCalRCBckgFlucSignal,1);
    fRhoChargedCMSScale->FillDeltaPtNColl(fEventCentrality,kTRho,fJetR,fEMCalRCBckgFlucNColl,1);
    fRhoChargedCMSScale->FillBackgroundFluctuations(fEventCentrality,kTRho,fJetR);
    fRhoChargedCMSScale->FillLeadingJetPtRho(fEMCalFullJet->GetLeadingPt(),kTRho);
    fRhoChargedCMSScale->DoNEFAnalysis(fNEFSignalJetCut,fEMCalJetThreshold,fmyAKTFullJets,fEMCalFullJet->GetJets(),fEMCalFullJet->GetTotalJets(),fmyClusters,fOrgClusters,fEvent,fEMCALGeometry,fRecoUtil,fCells);
    fRhoChargedCMSScale->FillMiscJetStats(fmyAKTFullJets,fEMCalFullJet->GetJets(),fEMCalFullJet->GetTotalJets(),fOrgTracks,fOrgClusters,fVertex);
    fRhoChargedCMSScale->FillJetEventCentrality(fEMCalFullJet->GetLeadingPt(),fEvent);
    
    delete [] RhoArray;
    delete [] pTArray;
}

// Full Rho's
void AliAnalysisTaskFullpAJets::EstimateFullRho0()
{
    Int_t i;
    Double_t E_tracks_total=0.0;
    Double_t E_caloclusters_total=0.0;
    Double_t EMCal_rho=0.0;
    
    TLorentzVector *cluster_vec = new TLorentzVector;

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
        vcluster->GetMomentum(*cluster_vec,fVertex);
        E_caloclusters_total+=cluster_vec->Pt();
    }
    delete cluster_vec;
    
    //  Calculate the mean Background density
    EMCal_rho=(E_tracks_total+E_caloclusters_total)/fEMCalArea;
    fRhoFull=EMCal_rho;
    
    // Fill Histograms
    fRhoFull0->FillRho(fEventCentrality,EMCal_rho);
    fRhoFull0->FillBSJS(fEventCentrality,EMCal_rho,fEMCalJetThreshold,fmyAKTFullJets,fEMCalFullJet->GetJets(),fEMCalFullJet->GetTotalJets());
    fRhoFull0->FillDeltaPt(fEventCentrality,EMCal_rho,fJetR,fEMCalRCBckgFluc,1);
    fRhoFull0->FillDeltaPtSignal(fEventCentrality,EMCal_rho,fJetR,fEMCalRCBckgFlucSignal,1);
    fRhoFull0->FillDeltaPtNColl(fEventCentrality,EMCal_rho,fJetR,fEMCalRCBckgFlucNColl,1);
    fRhoFull0->FillBackgroundFluctuations(fEventCentrality,EMCal_rho,fJetR);
    fRhoFull0->FillLeadingJetPtRho(fEMCalFullJet->GetLeadingPt(),EMCal_rho);
}

void AliAnalysisTaskFullpAJets::EstimateFullRho1()
{
    Int_t i;
    Double_t E_tracks_total=0.0;
    Double_t E_caloclusters_total=0.0;
    Double_t EMCal_rho=0.0;
    
    TLorentzVector *temp_jet= new TLorentzVector;
    TLorentzVector *track_vec = new TLorentzVector;
    TLorentzVector *cluster_vec = new TLorentzVector;

    if (fEMCalPartJetUnbiased->GetLeadingPt()>=fEMCalJetThreshold)
    {
        AliEmcalJet *myJet =(AliEmcalJet*) fmyAKTFullJets->At(fEMCalPartJetUnbiased->GetLeadingIndex());
        myJet->GetMom(*temp_jet);
        
        //  Loop over all tracks
        for (i=0;i<fnTracks;i++)
        {
            AliVTrack* vtrack =(AliVTrack*) fmyTracks->At(i);
            if (IsInEMCal(vtrack->Phi(),vtrack->Eta())==kTRUE)
            {
                track_vec->SetPtEtaPhiE(vtrack->Pt(),vtrack->Eta(),vtrack->Phi(),vtrack->E());
                if (temp_jet->DeltaR(*track_vec)>fJetRForRho)
                {
                    E_tracks_total+=vtrack->Pt();
                }
            }
        }
        
        //  Loop over all caloclusters
        for (i=0;i<fnClusters;i++)
        {
            AliVCluster* vcluster = (AliVCluster*) fmyClusters->At(i);
            vcluster->GetMomentum(*cluster_vec,fVertex);
            if (temp_jet->DeltaR(*cluster_vec)>fJetRForRho)
            {
                E_caloclusters_total+=cluster_vec->Pt();
            }
        }
        //  Calculate the mean Background density
        EMCal_rho=(E_tracks_total+E_caloclusters_total)/(fEMCalArea-AreaWithinEMCal(fJetR,myJet->Phi(),myJet->Eta()));
    }
    else  // i.e. No signal jets -> same as total background density
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
            vcluster->GetMomentum(*cluster_vec,fVertex);
            E_caloclusters_total+=cluster_vec->Pt();
        }
        //  Calculate the mean Background density
        EMCal_rho=(E_tracks_total+E_caloclusters_total)/fEMCalArea;
    }
    delete temp_jet;
    delete track_vec;
    delete cluster_vec;

    // Fill histograms
    fRhoFull1->FillRho(fEventCentrality,EMCal_rho);
    fRhoFull1->FillBSJS(fEventCentrality,EMCal_rho,fEMCalJetThreshold,fmyAKTFullJets,fEMCalFullJet->GetJets(),fEMCalFullJet->GetTotalJets());
    fRhoFull1->FillDeltaPt(fEventCentrality,EMCal_rho,fJetR,fEMCalRCBckgFluc,1);
    fRhoFull1->FillDeltaPtSignal(fEventCentrality,EMCal_rho,fJetR,fEMCalRCBckgFlucSignal,1);
    fRhoFull1->FillDeltaPtNColl(fEventCentrality,EMCal_rho,fJetR,fEMCalRCBckgFlucNColl,1);
    fRhoFull1->FillBackgroundFluctuations(fEventCentrality,EMCal_rho,fJetR);
    fRhoFull1->FillLeadingJetPtRho(fEMCalFullJet->GetLeadingPt(),EMCal_rho);
}

void AliAnalysisTaskFullpAJets::EstimateFullRho2()
{
    Int_t i;
    Double_t E_tracks_total=0.0;
    Double_t E_caloclusters_total=0.0;
    Double_t EMCal_rho=0.0;
    
    TLorentzVector *temp_jet1 = new TLorentzVector;
    TLorentzVector *temp_jet2 = new TLorentzVector;
    TLorentzVector *track_vec = new TLorentzVector;
    TLorentzVector *cluster_vec = new TLorentzVector;

    if ((fEMCalPartJetUnbiased->GetLeadingPt()>=fEMCalJetThreshold) && (fEMCalPartJetUnbiased->GetSubLeadingPt()>=fEMCalJetThreshold))
    {
        AliEmcalJet *myhJet =(AliEmcalJet*) fmyAKTFullJets->At(fEMCalPartJetUnbiased->GetLeadingIndex());
        myhJet->GetMom(*temp_jet1);
        
        AliEmcalJet *myJet =(AliEmcalJet*) fmyAKTFullJets->At(fEMCalPartJetUnbiased->GetSubLeadingIndex());
        myJet->GetMom(*temp_jet2);
     
        //  Loop over all tracks
        for (i=0;i<fnTracks;i++)
        {
            AliVTrack* vtrack =(AliVTrack*) fmyTracks->At(i);
            if (IsInEMCal(vtrack->Phi(),vtrack->Eta())==kTRUE)
            {
                track_vec->SetPtEtaPhiE(vtrack->Pt(),vtrack->Eta(),vtrack->Phi(),vtrack->E());
                if ((temp_jet1->DeltaR(*track_vec)>fJetRForRho) && (temp_jet2->DeltaR(*track_vec)>fJetRForRho))
                {
                    E_tracks_total+=vtrack->Pt();
                }
            }
        }
        
        //  Loop over all caloclusters
        for (i=0;i<fnClusters;i++)
        {
            AliVCluster* vcluster = (AliVCluster*) fmyClusters->At(i);
            vcluster->GetMomentum(*cluster_vec,fVertex);
            if ((temp_jet1->DeltaR(*cluster_vec)>fJetRForRho) && (temp_jet2->DeltaR(*cluster_vec)>fJetRForRho))
            {
                E_caloclusters_total+=cluster_vec->Pt();
            }
        }

        //  Calculate the mean Background density
        EMCal_rho=(E_tracks_total+E_caloclusters_total)/(fEMCalArea-AreaWithinEMCal(fJetR,myhJet->Phi(),myhJet->Eta())-AreaWithinEMCal(fJetR,myJet->Phi(),myJet->Eta()));
    }
    else if (fEMCalPartJetUnbiased->GetLeadingPt()>=fEMCalJetThreshold)
    {
        AliEmcalJet *myJet =(AliEmcalJet*) fmyAKTFullJets->At(fEMCalPartJetUnbiased->GetLeadingIndex());
        myJet->GetMom(*temp_jet1);
        
        //  Loop over all tracks
        for (i=0;i<fnTracks;i++)
        {
            AliVTrack* vtrack =(AliVTrack*) fmyTracks->At(i);
            if (IsInEMCal(vtrack->Phi(),vtrack->Eta())==kTRUE)
            {
                track_vec->SetPtEtaPhiE(vtrack->Pt(),vtrack->Eta(),vtrack->Phi(),vtrack->E());
                if (temp_jet1->DeltaR(*track_vec)>fJetRForRho)
                {
                    E_tracks_total+=vtrack->Pt();
                }
            }
        }
        
        //  Loop over all caloclusters
        for (i=0;i<fnClusters;i++)
        {
            AliVCluster* vcluster = (AliVCluster*) fmyClusters->At(i);
            vcluster->GetMomentum(*cluster_vec,fVertex);
            if (temp_jet1->DeltaR(*cluster_vec)>fJetRForRho)
            {
                E_caloclusters_total+=cluster_vec->Pt();
            }
        }
        //  Calculate the mean Background density
        EMCal_rho=(E_tracks_total+E_caloclusters_total)/(fEMCalArea-AreaWithinEMCal(fJetR,myJet->Phi(),myJet->Eta()));
    }
    else  // i.e. No signal jets -> same as total background density
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
            vcluster->GetMomentum(*cluster_vec,fVertex);
            E_caloclusters_total+=cluster_vec->Pt();
        }
        //  Calculate the mean Background density
        EMCal_rho=(E_tracks_total+E_caloclusters_total)/fEMCalArea;
    }
    delete temp_jet1;
    delete temp_jet2;
    delete track_vec;
    delete cluster_vec;

    // Fill histograms
    fRhoFull2->FillRho(fEventCentrality,EMCal_rho);
    fRhoFull2->FillBSJS(fEventCentrality,EMCal_rho,fEMCalJetThreshold,fmyAKTFullJets,fEMCalFullJet->GetJets(),fEMCalFullJet->GetTotalJets());
    fRhoFull2->FillDeltaPt(fEventCentrality,EMCal_rho,fJetR,fEMCalRCBckgFluc,1);
    fRhoFull2->FillDeltaPtSignal(fEventCentrality,EMCal_rho,fJetR,fEMCalRCBckgFlucSignal,1);
    fRhoFull2->FillDeltaPtNColl(fEventCentrality,EMCal_rho,fJetR,fEMCalRCBckgFlucNColl,1);
    fRhoFull2->FillBackgroundFluctuations(fEventCentrality,EMCal_rho,fJetR);
    fRhoFull2->FillLeadingJetPtRho(fEMCalFullJet->GetLeadingPt(),EMCal_rho);
}

void AliAnalysisTaskFullpAJets::EstimateFullRhoN()
{
    Int_t i,j;
    Bool_t track_away_from_jet;
    Bool_t cluster_away_from_jet;
    Double_t E_tracks_total=0.0;
    Double_t E_caloclusters_total=0.0;
    Double_t EMCal_rho=0.0;
    Double_t jet_area_total=0.0;
    
    TLorentzVector *jet_vec= new TLorentzVector;
    TLorentzVector *track_vec = new TLorentzVector;
    TLorentzVector *cluster_vec = new TLorentzVector;

    // First, sum all tracks within the EMCal that are away from jet(s) above Pt Threshold
    for (i=0;i<fnTracks;i++)
    {
        // First, check if track is in the EMCal!!
        AliVTrack* vtrack = (AliVTrack*) fmyTracks->At(i);
        if (IsInEMCal(vtrack->Phi(),vtrack->Eta())==kTRUE)
        {
            if (fEMCalPartJetUnbiased->GetTotalSignalJets()<1)
            {
                E_tracks_total+=vtrack->Pt();
            }
            else
            {
                track_away_from_jet=kTRUE;
                j=0;
                track_vec->SetPtEtaPhiE(vtrack->Pt(),vtrack->Eta(),vtrack->Phi(),vtrack->E());
                while (track_away_from_jet==kTRUE && j<fEMCalPartJetUnbiased->GetTotalSignalJets())
                {
                    AliEmcalJet *myJet =(AliEmcalJet*) fmyAKTFullJets->At(fEMCalPartJetUnbiased->GetSignalJetIndex(j));
                    myJet->GetMom(*jet_vec);
                    if (track_vec->DeltaR(*jet_vec)<=fJetRForRho)
                    {
                        track_away_from_jet=kFALSE;
                    }
                    j++;
                }
                if (track_away_from_jet==kTRUE)
                {
                    E_tracks_total+=vtrack->Pt();
                }
            }
        }
    }
    
    // Next, sum all CaloClusters within the EMCal (obviously all clusters must be within EMCal!!) that are away from jet(s) above Pt Threshold
    for (i=0;i<fnClusters;i++)
    {
        AliVCluster* vcluster = (AliVCluster*) fmyClusters->At(i);
        vcluster->GetMomentum(*cluster_vec,fVertex);
        if (fEMCalPartJet->GetTotalSignalJets()<1)
        {
            E_caloclusters_total+=cluster_vec->Pt();
        }
        else
        {
            cluster_away_from_jet=kTRUE;
            j=0;
            
            while (cluster_away_from_jet==kTRUE && j<fEMCalPartJetUnbiased->GetTotalSignalJets())
            {
                AliEmcalJet *myJet =(AliEmcalJet*) fmyAKTFullJets->At(fEMCalPartJetUnbiased->GetSignalJetIndex(j));
                myJet->GetMom(*jet_vec);
                if (cluster_vec->DeltaR(*jet_vec)<=fJetRForRho)
                {
                    cluster_away_from_jet=kFALSE;
                }
                j++;
            }
            if (cluster_away_from_jet==kTRUE)
            {
                E_caloclusters_total+=cluster_vec->Pt();
            }
        }
    }
    
    // Determine area of all Jets that are within the EMCal
    if (fEMCalPartJet->GetTotalSignalJets()==0)
    {
        jet_area_total=0.0;
    }
    else
    {
        for (i=0;i<fEMCalPartJetUnbiased->GetTotalSignalJets();i++)
        {
            AliEmcalJet *myJet = (AliEmcalJet*) fmyAKTFullJets->At(fEMCalPartJetUnbiased->GetSignalJetIndex(i));
            jet_area_total+=AreaWithinEMCal(fJetR,myJet->Phi(),myJet->Eta());
        }
    }
    delete jet_vec;
    delete track_vec;
    delete cluster_vec;

    // Calculate Rho
    EMCal_rho=(E_tracks_total+E_caloclusters_total)/(fEMCalArea-jet_area_total);
    
    // Fill Histogram
    fRhoFullN->FillRho(fEventCentrality,EMCal_rho);
    fRhoFullN->FillBSJS(fEventCentrality,EMCal_rho,fEMCalJetThreshold,fmyAKTFullJets,fEMCalFullJet->GetJets(),fEMCalFullJet->GetTotalJets());
    fRhoFullN->FillDeltaPt(fEventCentrality,EMCal_rho,fJetR,fEMCalRCBckgFluc,1);
    fRhoFullN->FillDeltaPtSignal(fEventCentrality,EMCal_rho,fJetR,fEMCalRCBckgFlucSignal,1);
    fRhoFullN->FillDeltaPtNColl(fEventCentrality,EMCal_rho,fJetR,fEMCalRCBckgFlucNColl,1);
    fRhoFullN->FillBackgroundFluctuations(fEventCentrality,EMCal_rho,fJetR);
    fRhoFullN->FillLeadingJetPtRho(fEMCalFullJet->GetLeadingPt(),EMCal_rho);
}

void AliAnalysisTaskFullpAJets::EstimateFullRhoDijet()
{
    Int_t i;
    Double_t E_tracks_total=0.0;
    Double_t E_caloclusters_total=0.0;
    Double_t EMCal_rho=0.0;
    
    TLorentzVector *cluster_vec = new TLorentzVector;

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
        vcluster->GetMomentum(*cluster_vec,fVertex);
        E_caloclusters_total+=cluster_vec->Pt();
    }
    
    delete cluster_vec;

    //  Calculate the mean Background density
    EMCal_rho=(E_tracks_total+E_caloclusters_total)/fEMCalArea;
    
    // Fill Histograms
    fRhoFullDijet->FillRho(fEventCentrality,EMCal_rho);
    fRhoFullDijet->FillBSJS(fEventCentrality,EMCal_rho,fEMCalJetThreshold,fmyAKTFullJets,fEMCalFullJet->GetJets(),fEMCalFullJet->GetTotalJets());
    fRhoFullDijet->FillDeltaPt(fEventCentrality,EMCal_rho,fJetR,fEMCalRCBckgFluc,1);
    fRhoFullDijet->FillDeltaPtSignal(fEventCentrality,EMCal_rho,fJetR,fEMCalRCBckgFlucSignal,1);
    fRhoFullDijet->FillDeltaPtNColl(fEventCentrality,EMCal_rho,fJetR,fEMCalRCBckgFlucNColl,1);
    fRhoFullDijet->FillBackgroundFluctuations(fEventCentrality,EMCal_rho,fJetR);
    fRhoFullDijet->FillLeadingJetPtRho(fEMCalFullJet->GetLeadingPt(),EMCal_rho);
}

void AliAnalysisTaskFullpAJets::EstimateFullRhokT()
{
    Int_t i;
    Double_t kTRho = 0.0;
    Double_t *pTArray = new Double_t[fEMCalkTFullJet->GetTotalJets()];
    Double_t *RhoArray = new Double_t[fEMCalkTFullJet->GetTotalJets()];
    
    for (i=0;i<fEMCalkTFullJet->GetTotalJets();i++)
    {
        AliEmcalJet *myJet =(AliEmcalJet*) fmyKTFullJets->At(fEMCalkTFullJet->GetJetIndex(i));
        pTArray[i]=myJet->Pt();
        RhoArray[i]=myJet->Pt()/myJet->Area();
    }
    
    if (fEMCalkTFullJet->GetTotalJets()>0)
    {
        kTRho=MedianRhokT(pTArray,RhoArray,fEMCalkTFullJet->GetTotalJets());
    }
    else
    {
        kTRho=0.0;
    }
    fRhoFullkT->FillRho(fEventCentrality,kTRho);
    fRhoFullkT->FillBSJS(fEventCentrality,kTRho,fEMCalJetThreshold,fmyAKTFullJets,fEMCalFullJet->GetJets(),fEMCalFullJet->GetTotalJets());
    fRhoFullkT->FillDeltaPt(fEventCentrality,kTRho,fJetR,fEMCalRCBckgFluc,1);
    fRhoFullkT->FillDeltaPtSignal(fEventCentrality,kTRho,fJetR,fEMCalRCBckgFlucSignal,1);
    fRhoFullkT->FillDeltaPtNColl(fEventCentrality,kTRho,fJetR,fEMCalRCBckgFlucNColl,1);
    fRhoFullkT->FillBackgroundFluctuations(fEventCentrality,kTRho,fJetR);
    fRhoFullkT->FillLeadingJetPtRho(fEMCalFullJet->GetLeadingPt(),kTRho);
    delete [] RhoArray;
    delete [] pTArray;
}

void AliAnalysisTaskFullpAJets::EstimateFullRhoCMS()
{
    Int_t i,k;
    Double_t kTRho = 0.0;
    Double_t CMSTotalkTArea = 0.0;
    Double_t CMSParticleArea = 0.0;
    Double_t CMSCorrectionFactor = 1.0;
    Double_t *pTArray = new Double_t[fEMCalkTFullJet->GetTotalJets()];
    Double_t *RhoArray = new Double_t[fEMCalkTFullJet->GetTotalJets()];
    
    k=0;
    if ((fEMCalPartJet->GetLeadingPt()>=fEMCalJetThreshold) && (fEMCalPartJet->GetSubLeadingPt()>=fEMCalJetThreshold))
    {
        AliEmcalJet *myJet1 =(AliEmcalJet*) fmyAKTFullJets->At(fEMCalPartJet->GetLeadingIndex());
        AliEmcalJet *myJet2 =(AliEmcalJet*) fmyAKTFullJets->At(fEMCalPartJet->GetSubLeadingIndex());
        
        for (i=0;i<fEMCalkTFullJet->GetTotalJets();i++)
        {
            AliEmcalJet *myJet =(AliEmcalJet*) fmyKTFullJets->At(fEMCalkTFullJet->GetJetIndex(i));
            
            CMSTotalkTArea+=myJet->Area();
            if (myJet->GetNumberOfTracks()>0 || myJet->GetNumberOfClusters()>0)
            {
                CMSParticleArea+=myJet->Area();
            }
            if (IsJetOverlap(myJet,myJet1,kTRUE)==kFALSE && IsJetOverlap(myJet,myJet2,kFALSE)==kTRUE)
            {
                pTArray[k]=myJet->Pt();
                RhoArray[k]=myJet->Pt()/myJet->Area();
                k++;
            }
        }
        if (k>0)
        {
            kTRho=MedianRhokT(pTArray,RhoArray,k);
        }
        else
        {
            kTRho=0.0;
        }
    }
    else if (fEMCalJet->GetLeadingPt()>=fEMCalJetThreshold)
    {
        AliEmcalJet *myJet1 =(AliEmcalJet*) fmyAKTFullJets->At(fEMCalJet->GetLeadingIndex());
        
        for (i=0;i<fEMCalkTFullJet->GetTotalJets();i++)
        {
            AliEmcalJet *myJet =(AliEmcalJet*) fmyKTFullJets->At(fEMCalkTFullJet->GetJetIndex(i));
            
            CMSTotalkTArea+=myJet->Area();
            if (myJet->GetNumberOfTracks()>0 || myJet->GetNumberOfClusters()>0)
            {
                CMSParticleArea+=myJet->Area();
            }
            if (IsJetOverlap(myJet,myJet1,kTRUE)==kFALSE)
            {
                pTArray[k]=myJet->Pt();
                RhoArray[k]=myJet->Pt()/myJet->Area();
                k++;
            }
        }
        if (k>0)
        {
            kTRho=MedianRhokT(pTArray,RhoArray,k);
        }
        else
        {
            kTRho=0.0;
        }
    }
    else
    {
        for (i=0;i<fEMCalkTFullJet->GetTotalJets();i++)
        {
            AliEmcalJet *myJet =(AliEmcalJet*) fmyKTFullJets->At(fEMCalkTFullJet->GetJetIndex(i));
            
            CMSTotalkTArea+=myJet->Area();
            if (myJet->GetNumberOfTracks()>0 || myJet->GetNumberOfClusters()>0)
            {
                CMSParticleArea+=myJet->Area();
            }
            pTArray[k]=myJet->Pt();
            RhoArray[k]=myJet->Pt()/myJet->Area();
            k++;
        }
        if (k>0)
        {
            kTRho=MedianRhokT(pTArray,RhoArray,k);
        }
        else
        {
            kTRho=0.0;
        }
    }
    // Scale CMS Rho by Correction factor
    if (CMSTotalkTArea==0.0)
    {
        CMSCorrectionFactor = 1.0;
    }
    else
    {
        //CMSCorrectionFactor = CMSTrackArea/CMSTotalkTArea;
        CMSCorrectionFactor = CMSParticleArea/((fEMCalPhiTotal-2*fJetR)*(fEMCalEtaTotal-2*fJetR));  // The total physical area should be reduced by the eta & phi cuts due to looping over only fully contained kT jets within the EMCal
    }
    kTRho*=CMSCorrectionFactor;

    fRhoFullCMS->FillRho(fEventCentrality,kTRho);
    fRhoFullCMS->FillBSJS(fEventCentrality,kTRho,fEMCalJetThreshold,fmyAKTFullJets,fEMCalFullJet->GetJets(),fEMCalFullJet->GetTotalJets());
    fRhoFullCMS->FillDeltaPt(fEventCentrality,kTRho,fJetR,fEMCalRCBckgFluc,1);
    fRhoFullCMS->FillDeltaPtSignal(fEventCentrality,kTRho,fJetR,fEMCalRCBckgFlucSignal,1);
    fRhoFullCMS->FillDeltaPtNColl(fEventCentrality,kTRho,fJetR,fEMCalRCBckgFlucNColl,1);
    fRhoFullCMS->FillBackgroundFluctuations(fEventCentrality,kTRho,fJetR);
    fRhoFullCMS->FillLeadingJetPtRho(fEMCalFullJet->GetLeadingPt(),kTRho);
    delete [] RhoArray;
    delete [] pTArray;
}

void AliAnalysisTaskFullpAJets::FullJetEnergyDensityProfile()
{
    Int_t i,j;
    Double_t delta_R;
    Double_t dR=0.0;
    Double_t fullEDR = fFullEDJetR*100;
    const Int_t fullEDBins = (Int_t) fullEDR;
    Double_t ED_pT[fullEDBins];
    dR = fFullEDJetR/Double_t(fullEDBins);  // Should be 0.01 be default
    
    TLorentzVector *jet_vec= new TLorentzVector;
    TLorentzVector *track_vec = new TLorentzVector;
    TLorentzVector *cluster_vec = new TLorentzVector;

    for (i=0;i<fEMCalFullJet->GetTotalSignalJets();i++)
    {
        AliEmcalJet *myJet = (AliEmcalJet*) fmyAKTFullJets->At(fEMCalFullJet->GetSignalJetIndex(i));

        if (IsInEMCalFull(fFullEDJetR,myJet->Phi(),myJet->Eta())==kTRUE)
        {
            for (j=0;j<fullEDBins;j++)
            {
                ED_pT[j]=0;
            }
            myJet->GetMom(*jet_vec);
            
            // Sum all tracks in concentric rings around jet vertex
            for (j=0;j<fnTracks;j++)
            {
                AliVTrack* vtrack = (AliVTrack*) fmyTracks->At(j);
                track_vec->SetPtEtaPhiE(vtrack->Pt(),vtrack->Eta(),vtrack->Phi(),vtrack->E());
                delta_R=jet_vec->DeltaR(*track_vec);
                if (delta_R<fFullEDJetR)
                {
                    ED_pT[TMath::FloorNint((delta_R/fFullEDJetR)*fullEDBins)]+=vtrack->Pt();
                }
            }
            
            // Sum all clusters in concentric rings around jet vertex
            for (j=0;j<fnClusters;j++)
            {
                AliVCluster* vcluster = (AliVCluster*) fmyClusters->At(j);
                vcluster->GetMomentum(*cluster_vec,fVertex);
                delta_R=jet_vec->DeltaR(*cluster_vec);
                if (delta_R<fFullEDJetR)
                {
                    ED_pT[TMath::FloorNint((delta_R/fFullEDJetR)*fullEDBins)]+=cluster_vec->Pt();
                }
            }
            
            for (j=0;j<fullEDBins;j++)
            {
                ED_pT[j] /= TMath::Pi()*dR*dR*(2*j+1);
                fpFullJetEDProfile->Fill(myJet->Pt(),fEventCentrality,j*dR,ED_pT[j]);
            }
        }
    }
    delete cluster_vec;
    delete track_vec;
    delete jet_vec;
}

// Although this calculates the charged pT density in concentric rings around a charged jet contained fudically within the TPC (R=0.8), we actually scale the density by the scale factor (fScaleFactor) in order to compare to 'Full Charge Density' and 'Rho Median Occupancy Approach' (RhoChargedCMSScaled). We still keep the charged value for completeness.
void AliAnalysisTaskFullpAJets::ChargedJetEnergyDensityProfile()
{
    Int_t i,j;
    Double_t delta_R;
    Double_t dR=0.0;
    Double_t chargedEDR = fChargedEDJetR*100;
    const Int_t chargedEDBins = (Int_t) chargedEDR;
    Double_t ED_pT[chargedEDBins];
    dR = fChargedEDJetR/Double_t(chargedEDBins);  // Should be 0.01 be default

    TLorentzVector *jet_vec= new TLorentzVector;
    TLorentzVector *track_vec = new TLorentzVector;

    for (i=0;i<fTPCFullJet->GetTotalSignalJets();i++)
    {
        AliEmcalJet *myJet = (AliEmcalJet*) fmyAKTChargedJets->At(fTPCFullJet->GetSignalJetIndex(i));
        if (IsInTPC(fChargedEDJetR,myJet->Phi(),myJet->Eta(),kTRUE)==kTRUE)
        {
            for (j=0;j<chargedEDBins;j++)
            {
                ED_pT[j]=0;
            }
            myJet->GetMom(*jet_vec);
            
            // Sum all tracks in concentric rings around jet vertex
            for (j=0;j<fnTracks;j++)
            {
                AliVTrack* vtrack = (AliVTrack*) fmyTracks->At(j);
                track_vec->SetPtEtaPhiE(vtrack->Pt(),vtrack->Eta(),vtrack->Phi(),vtrack->E());
                delta_R=jet_vec->DeltaR(*track_vec);
                if (delta_R<fChargedEDJetR)
                {
                    ED_pT[TMath::FloorNint((delta_R/fChargedEDJetR)*chargedEDBins)]+=vtrack->Pt();
                }
            }
            for (j=0;j<chargedEDBins;j++)
            {
                ED_pT[j] /= TMath::Pi()*dR*dR*(2*j+1);
                fpChargedJetEDProfile->Fill(myJet->Pt(),fEventCentrality,j*dR,ED_pT[j]);
                fpChargedJetEDProfileScaled->Fill(myJet->Pt(),fEventCentrality,j*dR,fScaleFactor*ED_pT[j]);
            }
        }
    }
    delete track_vec;
    delete jet_vec;
}

void AliAnalysisTaskFullpAJets::DeleteJetData(Int_t delOption)
{
    if (delOption ==0)
    {
        delete fRecoUtil;
        
        fRecoUtil = NULL;
    }
    else if (delOption==1)
    {
        delete fRecoUtil;
        delete fmyTracks;
        delete fmyClusters;
        
        fRecoUtil = NULL;
        fmyTracks = NULL;
        fmyClusters = NULL;
    }
    else if (delOption==2)
    {
        delete fRecoUtil;
        delete fmyTracks;
        delete fmyClusters;
        delete fTPCJet;
        delete fTPCFullJet;
        delete fTPCOnlyJet;
        delete fTPCJetUnbiased;
        delete fTPCkTFullJet;
        
        fRecoUtil = NULL;
        fmyTracks = NULL;
        fmyClusters = NULL;
        fTPCJet = NULL;
        fTPCFullJet = NULL;
        fTPCOnlyJet = NULL;
        fTPCJetUnbiased = NULL;
        fTPCkTFullJet = NULL;
    }
    if (delOption==3)
    {
        delete fRecoUtil;
        delete fmyTracks;
        delete fmyClusters;
        delete fTPCJet;
        delete fTPCFullJet;
        delete fTPCOnlyJet;
        delete fTPCJetUnbiased;
        delete fTPCkTFullJet;
        delete fEMCalJet;
        delete fEMCalFullJet;
        delete fEMCalPartJet;
        delete fEMCalPartJetUnbiased;
        delete fEMCalkTFullJet;
        
        fRecoUtil = NULL;
        fmyTracks = NULL;
        fmyClusters = NULL;
        fTPCJet = NULL;
        fTPCFullJet = NULL;
        fTPCOnlyJet = NULL;
        fTPCJetUnbiased = NULL;
        fEMCalJet = NULL;
        fEMCalFullJet = NULL;
        fEMCalPartJet = NULL;
        fEMCalPartJetUnbiased = NULL;
        fEMCalkTFullJet = NULL;
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

    const Double_t dijet_delta_phi=(180/360.)*2*TMath::Pi();
    const Double_t dijet_phi_acceptance=0.5*(30/360.)*2*TMath::Pi(); //Input the total acceptance within the paraenthesis to be +/- dijet_phi_acceptance
    Double_t dummy_phi=0.0;
    
    if (fTPCOnlyJet->GetTotalSignalJets()>1)
    {
        AliEmcalJet *myhJet =(AliEmcalJet*) fmyAKTChargedJets->At(fTPCOnlyJet->GetLeadingIndex());
        AliEmcalJet *myJet =(AliEmcalJet*) fmyAKTChargedJets->At(fTPCOnlyJet->GetSubLeadingIndex());
        dummy_phi=TMath::Min(TMath::Abs(myhJet->Phi()-myJet->Phi()),2*TMath::Pi()-TMath::Abs(myhJet->Phi()-myJet->Phi()));
        if (dummy_phi>(dijet_delta_phi-dijet_phi_acceptance))
        {
            fnDiJetEvents++;
            return kTRUE;
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

Bool_t AliAnalysisTaskFullpAJets::IsJetOverlap(AliEmcalJet *jet1,AliEmcalJet *jet2,Bool_t EMCalOn)
{
    Int_t i,j;
    Int_t jetTrack1=0;
    Int_t jetTrack2=0;
    Int_t jetCluster1=0;
    Int_t jetCluster2=0;
    
    for (i=0;i<jet1->GetNumberOfTracks();i++)
    {
        jetTrack1=jet1->TrackAt(i);
        for (j=0;j<jet2->GetNumberOfTracks();j++)
        {
            jetTrack2=jet2->TrackAt(j);
            if (jetTrack1 == jetTrack2)
            {
                return kTRUE;
            }
        }
    }
    if (EMCalOn == kTRUE)
    {
        for (i=0;i<jet1->GetNumberOfClusters();i++)
        {
            jetCluster1=jet1->ClusterAt(i);
            for (j=0;j<jet2->GetNumberOfClusters();j++)
            {
                jetCluster2=jet2->ClusterAt(j);
                if (jetCluster1 == jetCluster2)
                {
                    return kTRUE;
                }
            }
        }
    }
    return kFALSE;
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
    Double_t area_left=0;
    Double_t area_right=0;
    Double_t eta_a=0;
    Double_t eta_b=0;
    Double_t eta_up=0;
    Double_t eta_down=0;
    
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

    // Calculate the Right side area
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

Double_t AliAnalysisTaskFullpAJets::MedianRhokT(Double_t *pTkTEntries, Double_t *RhokTEntries, Int_t nEntries)
{
    // This function is used to calculate the median Rho kT value. The procedure is:
    // - Order the kT cluster array from highest rho value to lowest
    // - Exclude highest rho kT cluster
    // - Return the median rho value of the remaining subset
    
    // Sort Array
    const Double_t rho_min=-9.9999E+99;
    Int_t j,k;
    Double_t w[nEntries];  // Used for sorting
    Double_t smax=rho_min;
    Int_t sindex=-1;
    
    Double_t pTmax=0.0;
    Int_t pTmaxID=-1;
    
    for (j=0;j<nEntries;j++)
    {
        w[j]=0.0;
    }

    for (j=0;j<nEntries;j++)
    {
        if (pTkTEntries[j]>pTmax)
        {
            pTmax=pTkTEntries[j];
            pTmaxID=j;
        }
    }
    
    for (j=0;j<nEntries;j++)
    {
        for (k=0;k<nEntries;k++)
        {
            if (RhokTEntries[k]>smax)
            {
                smax=RhokTEntries[k];
                sindex=k;
            }
        }
        w[j]=smax;
        RhokTEntries[sindex]=rho_min;
        smax=rho_min;
        sindex=-1;
    }
    return w[nEntries/2];
}


// AlipAJetData Class Member Defs
// Constructors
AliAnalysisTaskFullpAJets::AlipAJetData::AlipAJetData() :

    fName(0),
    fIsJetsFull(0),
    fnTotal(0),
    fnJets(0),
    fnJetsSC(0),
    fJetR(0),
    fSignalPt(0),
    fAreaCutFrac(0.6),
    fNEF(1.0),
    fSignalTrackBias(0),
    fPtMaxIndex(0),
    fPtMax(0),
    fPtSubLeadingIndex(0),
    fPtSubLeading(0),
    fJetsIndex(0),
    fJetsSCIndex(0),
    fIsJetInArray(0),
    fJetMaxChargedPt(0)
{
    fnTotal=0;
    // Dummy constructor ALWAYS needed for I/O.
}

AliAnalysisTaskFullpAJets::AlipAJetData::AlipAJetData(const char *name, Bool_t isFull, Int_t nEntries) :

    fName(0),
    fIsJetsFull(0),
    fnTotal(0),
    fnJets(0),
    fnJetsSC(0),
    fJetR(0),
    fSignalPt(0),
    fAreaCutFrac(0.6),
    fNEF(1.0),
    fSignalTrackBias(0),
    fPtMaxIndex(0),
    fPtMax(0),
    fPtSubLeadingIndex(0),
    fPtSubLeading(0),
    fJetsIndex(0),
    fJetsSCIndex(0),
    fIsJetInArray(0),
    fJetMaxChargedPt(0)
{
    SetName(name);
    SetIsJetsFull(isFull);
    SetTotalEntries(nEntries);
    SetLeading(0,-9.99E+099);
    SetSubLeading(0,-9.99E+099);
    SetSignalCut(0);
    SetAreaCutFraction(0.6);
    SetJetR(fJetR);
    SetSignalTrackPtBias(0);
}

// Destructor
AliAnalysisTaskFullpAJets::AlipAJetData::~AlipAJetData()
{
    SetName("");
    SetIsJetsFull(kFALSE);
    SetTotalJets(0);
    SetTotalSignalJets(0);
    SetLeading(0,0);
    SetSubLeading(0,0);
    SetSignalCut(0);
    SetAreaCutFraction(0);
    SetJetR(0);
    SetNEF(0);
    SetSignalTrackPtBias(kFALSE);
    
    delete [] fJetsIndex;
    delete [] fJetsSCIndex;
    delete [] fIsJetInArray;
    delete [] fJetMaxChargedPt;
    
    fnTotal = 0;
}

// User Defined Sub-Routines
void AliAnalysisTaskFullpAJets::AlipAJetData::InitializeJetData(TClonesArray *jetList, Int_t nEntries)
{
    Int_t i=0;
    Int_t k=0;
    Int_t l=0;
    Double_t AreaThreshold = fAreaCutFrac*TMath::Pi()*TMath::Power(fJetR,2);
    
    // Initialize Jet Data
    for (i=0;i<nEntries;i++)
    {
        AliEmcalJet *myJet =(AliEmcalJet*) jetList->At(i);
        
        if (fIsJetInArray[i]==kTRUE && myJet->Area()>AreaThreshold)
        {
            SetJetIndex(i,k);
            if (myJet->Pt()>fPtMax)
            {
                SetSubLeading(fPtMaxIndex,fPtMax);
                SetLeading(i,myJet->Pt());
            }
            else if (myJet->Pt()>fPtSubLeading)
            {
                SetSubLeading(i,myJet->Pt()); 
            }
            // require leading charged constituent to have a pT greater then the signal threshold & Jet NEF to be less then the Signal Jet NEF cut
            fJetMaxChargedPt[i] = myJet->MaxTrackPt();
            if (fSignalTrackBias==kTRUE)
            {
                if (fJetMaxChargedPt[i]>=fSignalPt && myJet->NEF()<=fNEF)
                {
                    SetSignalJetIndex(i,l);
                    l++;
                }
            }
            else
            {
                if (myJet->Pt()>=fSignalPt && myJet->NEF()<=fNEF)
                {
                    SetSignalJetIndex(i,l);
                    l++;
                }
            }
            k++;
        }
    }
    SetTotalJets(k);
    SetTotalSignalJets(l);
}

// Setters
void AliAnalysisTaskFullpAJets::AlipAJetData::SetName(const char *name)
{
    fName = name;
}

void AliAnalysisTaskFullpAJets::AlipAJetData::SetIsJetsFull(Bool_t isFull)
{
    fIsJetsFull = isFull;
}

void AliAnalysisTaskFullpAJets::AlipAJetData::SetTotalEntries(Int_t nEntries)
{
    fnTotal = nEntries;
    fJetsIndex = new Int_t[fnTotal];
    fJetsSCIndex = new Int_t[fnTotal];
    fIsJetInArray = new Bool_t[fnTotal];
    fJetMaxChargedPt = new Double_t[fnTotal];
}

void AliAnalysisTaskFullpAJets::AlipAJetData::SetTotalJets(Int_t nJets)
{
    fnJets = nJets;
}

void AliAnalysisTaskFullpAJets::AlipAJetData::SetTotalSignalJets(Int_t nSignalJets)
{
    fnJetsSC = nSignalJets;
}

void AliAnalysisTaskFullpAJets::AlipAJetData::SetSignalCut(Double_t Pt)
{
    fSignalPt = Pt;
}

void AliAnalysisTaskFullpAJets::AlipAJetData::SetLeading(Int_t index, Double_t Pt)
{
    fPtMaxIndex = index;
    fPtMax = Pt;
}

void AliAnalysisTaskFullpAJets::AlipAJetData::SetSubLeading(Int_t index, Double_t Pt)
{
    fPtSubLeadingIndex = index;
    fPtSubLeading = Pt;
}

void AliAnalysisTaskFullpAJets::AlipAJetData::SetJetIndex(Int_t index, Int_t At)
{
    fJetsIndex[At] = index;
}

void AliAnalysisTaskFullpAJets::AlipAJetData::SetSignalJetIndex(Int_t index, Int_t At)
{
    fJetsSCIndex[At] = index;
}

void AliAnalysisTaskFullpAJets::AlipAJetData::SetIsJetInArray(Bool_t isInArray, Int_t At)
{
    fIsJetInArray[At] = isInArray;
}

void AliAnalysisTaskFullpAJets::AlipAJetData::SetAreaCutFraction(Double_t areaFraction)
{
    fAreaCutFrac = areaFraction;
}

void AliAnalysisTaskFullpAJets::AlipAJetData::SetJetR(Double_t jetR)
{
    fJetR = jetR;
}

void AliAnalysisTaskFullpAJets::AlipAJetData::SetNEF(Double_t nef)
{
    fNEF = nef;
}

void AliAnalysisTaskFullpAJets::AlipAJetData::SetSignalTrackPtBias(Bool_t chargedBias)
{
    fSignalTrackBias = chargedBias;
}

// Getters
Int_t AliAnalysisTaskFullpAJets::AlipAJetData::GetTotalEntries()
{
    return fnTotal;
}

Int_t AliAnalysisTaskFullpAJets::AlipAJetData::GetTotalJets()
{
    return fnJets;
}

Int_t AliAnalysisTaskFullpAJets::AlipAJetData::GetTotalSignalJets()
{
    return fnJetsSC;
}

Double_t AliAnalysisTaskFullpAJets::AlipAJetData::GetSignalCut()
{
    return fSignalPt;
}

Int_t AliAnalysisTaskFullpAJets::AlipAJetData::GetLeadingIndex()
{
    return fPtMaxIndex;
}

Double_t AliAnalysisTaskFullpAJets::AlipAJetData::GetLeadingPt()
{
    return fPtMax;
}

Int_t AliAnalysisTaskFullpAJets::AlipAJetData::GetSubLeadingIndex()
{
    return fPtSubLeadingIndex;
}

Double_t AliAnalysisTaskFullpAJets::AlipAJetData::GetSubLeadingPt()
{
    return fPtSubLeading;
}

Int_t AliAnalysisTaskFullpAJets::AlipAJetData::GetJetIndex(Int_t At)
{
    return fJetsIndex[At];
}

Int_t AliAnalysisTaskFullpAJets::AlipAJetData::GetSignalJetIndex(Int_t At)
{
    return fJetsSCIndex[At];
}

Bool_t AliAnalysisTaskFullpAJets::AlipAJetData::GetIsJetInArray(Int_t At)
{
    return fIsJetInArray[At];
}

Double_t AliAnalysisTaskFullpAJets::AlipAJetData::GetJetMaxChargedPt(Int_t At)
{
    return fJetMaxChargedPt[At];
}

Double_t AliAnalysisTaskFullpAJets::AlipAJetData::GetNEF()
{
    return fNEF;
}

// AlipAJetHistos Class Member Defs
// Constructors
AliAnalysisTaskFullpAJets::AlipAJetHistos::AlipAJetHistos() :

    fOutput(0),

    fh020Rho(0),
    fh80100Rho(0),
    fhRho(0),
    fhRhoCen(0),
    fh020BSPt(0),
    fh80100BSPt(0),
    fhBSPt(0),
    fhBSPtCen(0),
    fh020BSPtSignal(0),
    fh80100BSPtSignal(0),
    fhBSPtSignal(0),
    fhBSPtCenSignal(0),
    fh020DeltaPt(0),
    fh80100DeltaPt(0),
    fhDeltaPt(0),
    fhDeltaPtCen(0),
    fh020DeltaPtSignal(0),
    fh80100DeltaPtSignal(0),
    fhDeltaPtSignal(0),
    fhDeltaPtCenSignal(0),
    fh020DeltaPtNColl(0),
    fh80100DeltaPtNColl(0),
    fhDeltaPtNColl(0),
    fhDeltaPtCenNColl(0),
    fh020BckgFlucPt(0),
    fh80100BckgFlucPt(0),
    fhBckgFlucPt(0),
    fhBckgFlucPtCen(0),

    fpRho(0),
    fpLJetRho(0),

    fhJetPtEtaPhi(0),
    fhJetPtArea(0),
    fhJetConstituentPt(0),
    fhJetTracksPt(0),
    fhJetClustersPt(0),
    fhJetConstituentCounts(0),
    fhJetTracksCounts(0),
    fhJetClustersCounts(0),
    fhJetPtZConstituent(0),
    fhJetPtZTrack(0),
    fhJetPtZCluster(0),
    fhJetPtZLeadingConstituent(0),
    fhJetPtZLeadingTrack(0),
    fhJetPtZLeadingCluster(0),

    fhEventCentralityVsZNA(0),
    fhEventCentralityVsZNAPt(0),

    fNEFOutput(0),
    fhJetPtNEF(0),
    fhJetNEFInfo(0),
    fhJetNEFSignalInfo(0),
    fhClusterNEFInfo(0),
    fhClusterNEFSignalInfo(0),
    fhClusterShapeAll(0),
    fhClusterPtCellAll(0),

    fName(0),
    fCentralityTag(0),
    fCentralityBins(0),
    fCentralityLow(0),
    fCentralityUp(0),
    fPtBins(0),
    fPtLow(0),
    fPtUp(0),
    fRhoPtBins(0),
    fRhoPtLow(0),
    fRhoPtUp(0),
    fDeltaPtBins(0),
    fDeltaPtLow(0),
    fDeltaPtUp(0),
    fBckgFlucPtBins(0),
    fBckgFlucPtLow(0),
    fBckgFlucPtUp(0),
    fLJetPtBins(0),
    fLJetPtLow(0),
    fLJetPtUp(0),
    fRhoValue(0),
    fLChargedTrackPtBins(0),
    fLChargedTrackPtLow(0),
    fLChargedTrackPtUp(0),
    fDoNEFQAPlots(0),
    fDoNEFSignalOnly(1),
    fSignalTrackBias(0),
    fDoTHnSparse(0),
    fDo3DHistos(0),
    fNEFBins(0),
    fNEFLow(0),
    fNEFUp(0),
    fnDimJet(0),
    fnDimCluster(0),
    fEMCalPhiMin(1.39626),
    fEMCalPhiMax(3.26377),
    fEMCalEtaMin(-0.7),
    fEMCalEtaMax(0.7)

{
    // Dummy constructor ALWAYS needed for I/O.
}

AliAnalysisTaskFullpAJets::AlipAJetHistos::AlipAJetHistos(const char *name) :

    fOutput(0),

    fh020Rho(0),
    fh80100Rho(0),
    fhRho(0),
    fhRhoCen(0),
    fh020BSPt(0),
    fh80100BSPt(0),
    fhBSPt(0),
    fhBSPtCen(0),
    fh020BSPtSignal(0),
    fh80100BSPtSignal(0),
    fhBSPtSignal(0),
    fhBSPtCenSignal(0),
    fh020DeltaPt(0),
    fh80100DeltaPt(0),
    fhDeltaPt(0),
    fhDeltaPtCen(0),
    fh020DeltaPtSignal(0),
    fh80100DeltaPtSignal(0),
    fhDeltaPtSignal(0),
    fhDeltaPtCenSignal(0),
    fh020DeltaPtNColl(0),
    fh80100DeltaPtNColl(0),
    fhDeltaPtNColl(0),
    fhDeltaPtCenNColl(0),
    fh020BckgFlucPt(0),
    fh80100BckgFlucPt(0),
    fhBckgFlucPt(0),
    fhBckgFlucPtCen(0),

    fpRho(0),
    fpLJetRho(0),

    fhJetPtEtaPhi(0),
    fhJetPtArea(0),
    fhJetConstituentPt(0),
    fhJetTracksPt(0),
    fhJetClustersPt(0),
    fhJetConstituentCounts(0),
    fhJetTracksCounts(0),
    fhJetClustersCounts(0),
    fhJetPtZConstituent(0),
    fhJetPtZTrack(0),
    fhJetPtZCluster(0),
    fhJetPtZLeadingConstituent(0),
    fhJetPtZLeadingTrack(0),
    fhJetPtZLeadingCluster(0),

    fhEventCentralityVsZNA(0),
    fhEventCentralityVsZNAPt(0),

    fNEFOutput(0),
    fhJetPtNEF(0),
    fhJetNEFInfo(0),
    fhJetNEFSignalInfo(0),
    fhClusterNEFInfo(0),
    fhClusterNEFSignalInfo(0),
    fhClusterShapeAll(0),
    fhClusterPtCellAll(0),

    fName(0),
    fCentralityTag(0),
    fCentralityBins(0),
    fCentralityLow(0),
    fCentralityUp(0),
    fPtBins(0),
    fPtLow(0),
    fPtUp(0),
    fRhoPtBins(0),
    fRhoPtLow(0),
    fRhoPtUp(0),
    fDeltaPtBins(0),
    fDeltaPtLow(0),
    fDeltaPtUp(0),
    fBckgFlucPtBins(0),
    fBckgFlucPtLow(0),
    fBckgFlucPtUp(0),
    fLJetPtBins(0),
    fLJetPtLow(0),
    fLJetPtUp(0),
    fRhoValue(0),
    fLChargedTrackPtBins(0),
    fLChargedTrackPtLow(0),
    fLChargedTrackPtUp(0),
    fDoNEFQAPlots(0),
    fDoNEFSignalOnly(1),
    fSignalTrackBias(0),
    fDoTHnSparse(0),
    fDo3DHistos(0),
    fNEFBins(0),
    fNEFLow(0),
    fNEFUp(0),
    fnDimJet(0),
    fnDimCluster(0),
    fEMCalPhiMin(1.39626),
    fEMCalPhiMax(3.26377),
    fEMCalEtaMin(-0.7),
    fEMCalEtaMax(0.7)

{
    SetName(name);
    SetCentralityTag("ZNA");
    SetCentralityRange(100,0,100);
    SetPtRange(250,-50,200);
    SetRhoPtRange(500,0,50);
    SetDeltaPtRange(200,-100,100);
    SetBackgroundFluctuationsPtRange(100,0,100);
    SetLeadingJetPtRange(200,0,200);
    SetLeadingChargedTrackPtRange(100,0,100);
    SetNEFRange(100,0,1);
    DoNEFQAPlots(kFALSE);
    
    Init();
}

AliAnalysisTaskFullpAJets::AlipAJetHistos::AlipAJetHistos(const char *name, TString centag, Bool_t doNEF) :

    fOutput(0),

    fh020Rho(0),
    fh80100Rho(0),
    fhRho(0),
    fhRhoCen(0),
    fh020BSPt(0),
    fh80100BSPt(0),
    fhBSPt(0),
    fhBSPtCen(0),
    fh020BSPtSignal(0),
    fh80100BSPtSignal(0),
    fhBSPtSignal(0),
    fhBSPtCenSignal(0),
    fh020DeltaPt(0),
    fh80100DeltaPt(0),
    fhDeltaPt(0),
    fhDeltaPtCen(0),
    fh020DeltaPtSignal(0),
    fh80100DeltaPtSignal(0),
    fhDeltaPtSignal(0),
    fhDeltaPtCenSignal(0),
    fh020DeltaPtNColl(0),
    fh80100DeltaPtNColl(0),
    fhDeltaPtNColl(0),
    fhDeltaPtCenNColl(0),
    fh020BckgFlucPt(0),
    fh80100BckgFlucPt(0),
    fhBckgFlucPt(0),
    fhBckgFlucPtCen(0),

    fpRho(0),
    fpLJetRho(0),

    fhJetPtEtaPhi(0),
    fhJetPtArea(0),
    fhJetConstituentPt(0),
    fhJetTracksPt(0),
    fhJetClustersPt(0),
    fhJetConstituentCounts(0),
    fhJetTracksCounts(0),
    fhJetClustersCounts(0),
    fhJetPtZConstituent(0),
    fhJetPtZTrack(0),
    fhJetPtZCluster(0),
    fhJetPtZLeadingConstituent(0),
    fhJetPtZLeadingTrack(0),
    fhJetPtZLeadingCluster(0),

    fhEventCentralityVsZNA(0),
    fhEventCentralityVsZNAPt(0),

    fNEFOutput(0),
    fhJetPtNEF(0),
    fhJetNEFInfo(0),
    fhJetNEFSignalInfo(0),
    fhClusterNEFInfo(0),
    fhClusterNEFSignalInfo(0),
    fhClusterShapeAll(0),
    fhClusterPtCellAll(0),

    fName(0),
    fCentralityTag(0),
    fCentralityBins(0),
    fCentralityLow(0),
    fCentralityUp(0),
    fPtBins(0),
    fPtLow(0),
    fPtUp(0),
    fRhoPtBins(0),
    fRhoPtLow(0),
    fRhoPtUp(0),
    fDeltaPtBins(0),
    fDeltaPtLow(0),
    fDeltaPtUp(0),
    fBckgFlucPtBins(0),
    fBckgFlucPtLow(0),
    fBckgFlucPtUp(0),
    fLJetPtBins(0),
    fLJetPtLow(0),
    fLJetPtUp(0),
    fRhoValue(0),
    fLChargedTrackPtBins(0),
    fLChargedTrackPtLow(0),
    fLChargedTrackPtUp(0),
    fDoNEFQAPlots(0),
    fDoNEFSignalOnly(1),
    fSignalTrackBias(0),
    fDoTHnSparse(0),
    fDo3DHistos(0),
    fNEFBins(0),
    fNEFLow(0),
    fNEFUp(0),
    fnDimJet(0),
    fnDimCluster(0),
    fEMCalPhiMin(1.39626),
    fEMCalPhiMax(3.26377),
    fEMCalEtaMin(-0.7),
    fEMCalEtaMax(0.7)

{
    SetName(name);
    SetCentralityTag(centag.Data());
    SetCentralityRange(100,0,100);
    SetPtRange(250,-50,200);
    SetRhoPtRange(500,0,50);
    SetDeltaPtRange(200,-100,100);
    SetBackgroundFluctuationsPtRange(100,0,100);
    SetLeadingJetPtRange(200,0,200);
    SetLeadingChargedTrackPtRange(100,0,100);
    SetNEFRange(100,0,1);
    DoNEFQAPlots(doNEF);

    Init();
}

AliAnalysisTaskFullpAJets::AlipAJetHistos::AlipAJetHistos(const char *name, TString centag, Bool_t doNEF, Bool_t doNEFSignalOnly, Bool_t doTHnSparse, Bool_t do3DPlotting) :

    fOutput(0),

    fh020Rho(0),
    fh80100Rho(0),
    fhRho(0),
    fhRhoCen(0),
    fh020BSPt(0),
    fh80100BSPt(0),
    fhBSPt(0),
    fhBSPtCen(0),
    fh020BSPtSignal(0),
    fh80100BSPtSignal(0),
    fhBSPtSignal(0),
    fhBSPtCenSignal(0),
    fh020DeltaPt(0),
    fh80100DeltaPt(0),
    fhDeltaPt(0),
    fhDeltaPtCen(0),
    fh020DeltaPtSignal(0),
    fh80100DeltaPtSignal(0),
    fhDeltaPtSignal(0),
    fhDeltaPtCenSignal(0),
    fh020DeltaPtNColl(0),
    fh80100DeltaPtNColl(0),
    fhDeltaPtNColl(0),
    fhDeltaPtCenNColl(0),
    fh020BckgFlucPt(0),
    fh80100BckgFlucPt(0),
    fhBckgFlucPt(0),
    fhBckgFlucPtCen(0),

    fpRho(0),
    fpLJetRho(0),

    fhJetPtEtaPhi(0),
    fhJetPtArea(0),
    fhJetConstituentPt(0),
    fhJetTracksPt(0),
    fhJetClustersPt(0),
    fhJetConstituentCounts(0),
    fhJetTracksCounts(0),
    fhJetClustersCounts(0),
    fhJetPtZConstituent(0),
    fhJetPtZTrack(0),
    fhJetPtZCluster(0),
    fhJetPtZLeadingConstituent(0),
    fhJetPtZLeadingTrack(0),
    fhJetPtZLeadingCluster(0),

    fhEventCentralityVsZNA(0),
    fhEventCentralityVsZNAPt(0),

    fNEFOutput(0),
    fhJetPtNEF(0),
    fhJetNEFInfo(0),
    fhJetNEFSignalInfo(0),
    fhClusterNEFInfo(0),
    fhClusterNEFSignalInfo(0),
    fhClusterShapeAll(0),
    fhClusterPtCellAll(0),

    fName(0),
    fCentralityTag(0),
    fCentralityBins(0),
    fCentralityLow(0),
    fCentralityUp(0),
    fPtBins(0),
    fPtLow(0),
    fPtUp(0),
    fRhoPtBins(0),
    fRhoPtLow(0),
    fRhoPtUp(0),
    fDeltaPtBins(0),
    fDeltaPtLow(0),
    fDeltaPtUp(0),
    fBckgFlucPtBins(0),
    fBckgFlucPtLow(0),
    fBckgFlucPtUp(0),
    fLJetPtBins(0),
    fLJetPtLow(0),
    fLJetPtUp(0),
    fRhoValue(0),
    fLChargedTrackPtBins(0),
    fLChargedTrackPtLow(0),
    fLChargedTrackPtUp(0),
    fDoNEFQAPlots(0),
    fDoNEFSignalOnly(1),
    fSignalTrackBias(0),
    fDoTHnSparse(0),
    fDo3DHistos(0),
    fNEFBins(0),
    fNEFLow(0),
    fNEFUp(0),
    fnDimJet(0),
    fnDimCluster(0),
    fEMCalPhiMin(1.39626),
    fEMCalPhiMax(3.26377),
    fEMCalEtaMin(-0.7),
    fEMCalEtaMax(0.7)

{
    SetName(name);
    SetCentralityTag(centag.Data());
    SetCentralityRange(100,0,100);
    SetPtRange(250,-50,200);
    SetRhoPtRange(500,0,50);
    SetDeltaPtRange(200,-100,100);
    SetBackgroundFluctuationsPtRange(100,0,100);
    SetLeadingJetPtRange(200,0,200);
    SetLeadingChargedTrackPtRange(100,0,100);
    SetNEFRange(100,0,1);
    DoNEFQAPlots(doNEF);
    DoNEFSignalOnly(doNEFSignalOnly);
    DoTHnSparse(doTHnSparse);
    Do3DPlotting(do3DPlotting);
    
    Init();
}

// Destructor
AliAnalysisTaskFullpAJets::AlipAJetHistos::~AlipAJetHistos()
{
    if (fOutput)
    {
        delete fOutput;
    }
}

void AliAnalysisTaskFullpAJets::AlipAJetHistos::Init()
{
    Int_t i;
    
    // Initialize Private Variables
    Int_t TCBins = 100;
    fEMCalPhiMin=(80/(double)360)*2*TMath::Pi();
    fEMCalPhiMax=(187/(double)360)*2*TMath::Pi();
    fEMCalEtaMin=-0.7;
    fEMCalEtaMax=0.7;
    
    fOutput = new TList();
    fOutput->SetOwner();
    fOutput->SetName(fName);
    
    TString RhoString="";
    TString PtString="";
    TString DeltaPtString="";
    TString BckgFlucPtString="";
    TString CentralityString;
    TString EventCentralityString;
    CentralityString = Form("Centrality (%s) ",fCentralityTag.Data());
    EventCentralityString = Form("%s vs ZNA Centrality ",fCentralityTag.Data());
    
    // Rho Spectral Plots
    RhoString = Form("%d-%d Centrality, Rho Spectrum",0,20);
    fh020Rho = new TH1F("fh020Rho",RhoString,fRhoPtBins,fRhoPtLow,fRhoPtUp);
    fh020Rho->GetXaxis()->SetTitle("p_{T}/Area (GeV/c)");
    fh020Rho->GetYaxis()->SetTitle("1/N_{Events} dN/d#rho");
    fh020Rho->Sumw2();
    
    RhoString = Form("%d-%d Centrality, Rho Spectrum",80,100);
    fh80100Rho = new TH1F("fh80100Rho",RhoString,fRhoPtBins,fRhoPtLow,fRhoPtUp);
    fh80100Rho->GetXaxis()->SetTitle("p_{T}/Area (GeV/c)");
    fh80100Rho->GetYaxis()->SetTitle("1/N_{Events} dN/d#rho");
    fh80100Rho->Sumw2();
    
    RhoString = Form("%d-%d Centrality, Rho Spectrum",0,100);
    fhRho = new TH1F("fhRho",RhoString,fRhoPtBins,fRhoPtLow,fRhoPtUp);
    fhRho->GetXaxis()->SetTitle("p_{T}/Area (GeV/c)");
    fhRho->GetYaxis()->SetTitle("1/N_{Events} dN/d#rho");
    fhRho->Sumw2();
    
    RhoString = "Rho Spectrum vs Centrality";
    fhRhoCen = new TH2F("fhRhoCen",RhoString,fRhoPtBins,fRhoPtLow,fRhoPtUp,fCentralityBins,fCentralityLow,fCentralityUp);
    fhRhoCen->GetXaxis()->SetTitle("p_{T}/Area (GeV/c)");
    fhRhoCen->GetYaxis()->SetTitle(Form("%s",CentralityString.Data()));
    fhRhoCen->GetZaxis()->SetTitle("1/N_{Events} dN/d#rho");
    fhRhoCen->Sumw2();
    
    // Background Subtracted Plots
    PtString = Form("%d-%d Centrality, Background Subtracted Jet Spectrum",0,20);
    fh020BSPt = new TH1F("fh020BSPt",PtString,fPtBins,fPtLow,fPtUp);
    fh020BSPt->GetXaxis()->SetTitle("p_{T} - #rhoA (GeV/c)");
    fh020BSPt->GetYaxis()->SetTitle("1/N_{Events} dN/dp_{T}d#etad#varphi");
    fh020BSPt->Sumw2();
    
    PtString = Form("%d-%d Centrality, Background Subtracted Jet Spectrum",80,100);
    fh80100BSPt = new TH1F("fh80100BSPt",PtString,fPtBins,fPtLow,fPtUp);
    fh80100BSPt->GetXaxis()->SetTitle("p_{T} - #rhoA (GeV/c)");
    fh80100BSPt->GetYaxis()->SetTitle("1/N_{Events} dN/dp_{T}d#etad#varphi");
    fh80100BSPt->Sumw2();
    
    PtString = Form("%d-%d Centrality, Background Subtracted Jet Spectrum",0,100);
    fhBSPt = new TH1F("fhBSPt",PtString,fPtBins,fPtLow,fPtUp);
    fhBSPt->GetXaxis()->SetTitle("p_{T} - #rhoA (GeV/c)");
    fhBSPt->GetYaxis()->SetTitle("1/N_{Events} dN/dp_{T}d#etad#varphi");
    fhBSPt->Sumw2();
    
    PtString = "Background Subtracted Jet Spectrum vs Centrality";
    fhBSPtCen = new TH2F("fhBSPtCen",PtString,fPtBins,fPtLow,fPtUp,fCentralityBins,fCentralityLow,fCentralityUp);
    fhBSPtCen->GetXaxis()->SetTitle("p_{T} - #rhoA (GeV/c)");
    fhBSPtCen->GetYaxis()->SetTitle(Form("%s",CentralityString.Data()));
    fhBSPtCen->GetZaxis()->SetTitle("1/N_{Events} dN/dp_{T}d#etad#varphi");
    fhBSPtCen->Sumw2();

    PtString = Form("%d-%d Centrality, Background Subtracted Signal Jet Spectrum",0,20);
    fh020BSPtSignal = new TH1F("fh020BSPtSignal",PtString,fPtBins,fPtLow,fPtUp);
    fh020BSPtSignal->GetXaxis()->SetTitle("p_{T} - #rhoA (GeV/c)");
    fh020BSPtSignal->GetYaxis()->SetTitle("1/N_{Events} dN/dp_{T}d#etad#varphi");
    fh020BSPtSignal->Sumw2();
    
    PtString = Form("%d-%d Centrality, Background Subtracted Signal Jet Spectrum",80,100);
    fh80100BSPtSignal = new TH1F("fh80100BSPtSignal",PtString,fPtBins,fPtLow,fPtUp);
    fh80100BSPtSignal->GetXaxis()->SetTitle("p_{T} - #rhoA (GeV/c)");
    fh80100BSPtSignal->GetYaxis()->SetTitle("1/N_{Events} dN/dp_{T}d#etad#varphi");
    fh80100BSPtSignal->Sumw2();
    
    PtString = Form("%d-%d Centrality, Background Subtracted Signal Jet Spectrum",0,100);
    fhBSPtSignal = new TH1F("fhBSPtSignal",PtString,fPtBins,fPtLow,fPtUp);
    fhBSPtSignal->GetXaxis()->SetTitle("p_{T} - #rhoA (GeV/c)");
    fhBSPtSignal->GetYaxis()->SetTitle("1/N_{Events} dN/dp_{T}d#etad#varphi");
    fhBSPtSignal->Sumw2();
    
    PtString = "Background Subtracted Signal Jet Spectrum vs Centrality";
    fhBSPtCenSignal = new TH2F("fhBSPtCenSignal",PtString,fPtBins,fPtLow,fPtUp,fCentralityBins,fCentralityLow,fCentralityUp);
    fhBSPtCenSignal->GetXaxis()->SetTitle("p_{T} - #rhoA (GeV/c)");
    fhBSPtCenSignal->GetYaxis()->SetTitle(Form("%s",CentralityString.Data()));
    fhBSPtCenSignal->GetZaxis()->SetTitle("1/N_{Events} dN/dp_{T}d#etad#varphi");
    fhBSPtCenSignal->Sumw2();
    
    // Delta Pt Plots with RC at least 2R away from Leading Signal
    DeltaPtString = Form("%d-%d Centrality, #deltap_{T} Spectrum",0,20);
    fh020DeltaPt = new TH1F("fh020DeltaPt",DeltaPtString,fDeltaPtBins,fDeltaPtLow,fDeltaPtUp);
    fh020DeltaPt->GetXaxis()->SetTitle("#deltap_{T} (GeV/c)");
    fh020DeltaPt->GetYaxis()->SetTitle("Probability Density");
    fh020DeltaPt->Sumw2();
    
    DeltaPtString = Form("%d-%d Centrality, #deltap_{T} Spectrum",80,100);
    fh80100DeltaPt = new TH1F("fh80100DeltaPt",DeltaPtString,fDeltaPtBins,fDeltaPtLow,fDeltaPtUp);
    fh80100DeltaPt->GetXaxis()->SetTitle("#deltap_{T} (GeV/c)");
    fh80100DeltaPt->GetYaxis()->SetTitle("Probability Density");
    fh80100DeltaPt->Sumw2();
    
    DeltaPtString = Form("%d-%d Centrality, #deltap_{T} Spectrum",0,100);
    fhDeltaPt = new TH1F("fhDeltaPt",DeltaPtString,fDeltaPtBins,fDeltaPtLow,fDeltaPtUp);
    fhDeltaPt->GetXaxis()->SetTitle("#deltap_{T} (GeV/c)");
    fhDeltaPt->GetYaxis()->SetTitle("Probability Density");
    fhDeltaPt->Sumw2();
    
    DeltaPtString = "#deltap_{T} Spectrum vs Centrality";
    fhDeltaPtCen = new TH2F("fhDeltaPtCen",DeltaPtString,fDeltaPtBins,fDeltaPtLow,fDeltaPtUp,fCentralityBins,fCentralityLow,fCentralityUp);
    fhDeltaPtCen->GetXaxis()->SetTitle("#deltap_{T} (GeV/c)");
    fhDeltaPtCen->GetYaxis()->SetTitle(Form("%s",CentralityString.Data()));
    fhDeltaPtCen->GetZaxis()->SetTitle("Probability Density");
    fhDeltaPtCen->Sumw2();
    
    // Delta Pt Plots with no spatial restrictions on RC
    DeltaPtString = Form("%d-%d Centrality, #deltap_{T} Spectrum",0,20);
    fh020DeltaPtSignal = new TH1F("fh020DeltaPtSignal",DeltaPtString,fDeltaPtBins,fDeltaPtLow,fDeltaPtUp);
    fh020DeltaPtSignal->GetXaxis()->SetTitle("#deltap_{T} (GeV/c)");
    fh020DeltaPtSignal->GetYaxis()->SetTitle("Probability Density");
    fh020DeltaPtSignal->Sumw2();
    
    DeltaPtString = Form("%d-%d Centrality, #deltap_{T} Spectrum",80,100);
    fh80100DeltaPtSignal = new TH1F("fh80100DeltaPtSignal",DeltaPtString,fDeltaPtBins,fDeltaPtLow,fDeltaPtUp);
    fh80100DeltaPtSignal->GetXaxis()->SetTitle("#deltap_{T} (GeV/c)");
    fh80100DeltaPtSignal->GetYaxis()->SetTitle("Probability Density");
    fh80100DeltaPtSignal->Sumw2();
    
    DeltaPtString = Form("%d-%d Centrality, #deltap_{T} Spectrum",0,100);
    fhDeltaPtSignal = new TH1F("fhDeltaPtSignal",DeltaPtString,fDeltaPtBins,fDeltaPtLow,fDeltaPtUp);
    fhDeltaPtSignal->GetXaxis()->SetTitle("#deltap_{T} (GeV/c)");
    fhDeltaPtSignal->GetYaxis()->SetTitle("Probability Density");
    fhDeltaPtSignal->Sumw2();
    
    DeltaPtString = "#deltap_{T} Spectrum vs Centrality";
    fhDeltaPtCenSignal = new TH2F("fhDeltaPtCenSignal",DeltaPtString,fDeltaPtBins,fDeltaPtLow,fDeltaPtUp,fCentralityBins,fCentralityLow,fCentralityUp);
    fhDeltaPtCenSignal->GetXaxis()->SetTitle("#deltap_{T} (GeV/c)");
    fhDeltaPtCenSignal->GetYaxis()->SetTitle(Form("%s",CentralityString.Data()));
    fhDeltaPtCenSignal->GetZaxis()->SetTitle("Probability Density");
    fhDeltaPtCenSignal->Sumw2();

    // Delta Pt Plots with NColl restrictions on RC
    DeltaPtString = Form("%d-%d Centrality, #deltap_{T} Spectrum",0,20);
    fh020DeltaPtNColl = new TH1F("fh020DeltaPtNColl",DeltaPtString,fDeltaPtBins,fDeltaPtLow,fDeltaPtUp);
    fh020DeltaPtNColl->GetXaxis()->SetTitle("#deltap_{T} (GeV/c)");
    fh020DeltaPtNColl->GetYaxis()->SetTitle("Probability Density");
    fh020DeltaPtNColl->Sumw2();
    
    DeltaPtString = Form("%d-%d Centrality, #deltap_{T} Spectrum",80,100);
    fh80100DeltaPtNColl = new TH1F("fh80100DeltaPtNColl",DeltaPtString,fDeltaPtBins,fDeltaPtLow,fDeltaPtUp);
    fh80100DeltaPtNColl->GetXaxis()->SetTitle("#deltap_{T} (GeV/c)");
    fh80100DeltaPtNColl->GetYaxis()->SetTitle("Probability Density");
    fh80100DeltaPtNColl->Sumw2();
    
    DeltaPtString = Form("%d-%d Centrality, #deltap_{T} Spectrum",0,100);
    fhDeltaPtNColl = new TH1F("fhDeltaPtNColl",DeltaPtString,fDeltaPtBins,fDeltaPtLow,fDeltaPtUp);
    fhDeltaPtNColl->GetXaxis()->SetTitle("#deltap_{T} (GeV/c)");
    fhDeltaPtNColl->GetYaxis()->SetTitle("Probability Density");
    fhDeltaPtNColl->Sumw2();
    
    DeltaPtString = "#deltap_{T} Spectrum vs Centrality";
    fhDeltaPtCenNColl = new TH2F("fhDeltaPtCenNColl",DeltaPtString,fDeltaPtBins,fDeltaPtLow,fDeltaPtUp,fCentralityBins,fCentralityLow,fCentralityUp);
    fhDeltaPtCenNColl->GetXaxis()->SetTitle("#deltap_{T} (GeV/c)");
    fhDeltaPtCenNColl->GetYaxis()->SetTitle(Form("%s",CentralityString.Data()));
    fhDeltaPtCenNColl->GetZaxis()->SetTitle("Probability Density");
    fhDeltaPtCenNColl->Sumw2();

    // Background Fluctuations Pt Plots
    BckgFlucPtString = Form("%d-%d Centrality, Background Fluctuation p_{T} Spectrum",0,20);
    fh020BckgFlucPt = new TH1F("fh020BckgFlucPt",PtString,fPtBins,fPtLow,fPtUp);
    fh020BckgFlucPt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fh020BckgFlucPt->GetYaxis()->SetTitle("1/N_{Events} dN/dp_{T}d#etad#varphi");
    fh020BckgFlucPt->Sumw2();
    
    BckgFlucPtString = Form("%d-%d Centrality, Background Fluctuation p_{T} Spectrum",80,100);
    fh80100BckgFlucPt = new TH1F("fh80100BckgFlucPt",BckgFlucPtString,fBckgFlucPtBins,fBckgFlucPtLow,fBckgFlucPtUp);
    fh80100BckgFlucPt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fh80100BckgFlucPt->GetYaxis()->SetTitle("1/N_{Events} dN/dp_{T}d#etad#varphi");
    fh80100BckgFlucPt->Sumw2();
    
    BckgFlucPtString = Form("%d-%d Centrality, Background Fluctuation p_{T} Spectrum",0,100);
    fhBckgFlucPt = new TH1F("fhBckgFlucPt",BckgFlucPtString,fBckgFlucPtBins,fBckgFlucPtLow,fBckgFlucPtUp);
    fhBckgFlucPt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fhBckgFlucPt->GetYaxis()->SetTitle("1/N_{Events} dN/dp_{T}d#etad#varphi");
    fhBckgFlucPt->Sumw2();
    
    BckgFlucPtString = "Background Fluctuation p_{T} Spectrum vs Centrality";
    fhBckgFlucPtCen = new TH2F("fhBckgFlucPtCen",BckgFlucPtString,fBckgFlucPtBins,fBckgFlucPtLow,fBckgFlucPtUp,fCentralityBins,fCentralityLow,fCentralityUp);
    fhBckgFlucPtCen->GetXaxis()->SetTitle("#p_{T} (GeV/c)");
    fhBckgFlucPtCen->GetYaxis()->SetTitle(Form("%s",CentralityString.Data()));
    fhBckgFlucPtCen->GetZaxis()->SetTitle("1/N_{Events} dN/dp_{T}d#etad#varphi");
    fhBckgFlucPtCen->Sumw2();
    
    // Background Density vs Centrality Profile
    RhoString = "Background Density vs Centrality";
    fpRho = new TProfile("fpRho",RhoString,fCentralityBins,fCentralityLow,fCentralityUp);
    fpRho->GetXaxis()->SetTitle(Form("%s",CentralityString.Data()));
    fpRho->GetYaxis()->SetTitle("p_{T}/Area (GeV/c)");
    
    // Background Density vs Leading Jet Profile
    fpLJetRho = new TProfile("fpLJetRho","#rho vs Leading Jet p_{T}",fLJetPtBins,fLJetPtLow,fLJetPtUp);
    fpLJetRho->GetXaxis()->SetTitle("Leading Jet p_{T}");
    fpLJetRho->GetYaxis()->SetTitle("p_{T}/Area (GeV/c)");

    // Jet pT vs Area
    Int_t JetPtAreaBins=200;
    Double_t JetPtAreaLow=0.0;
    Double_t JetPtAreaUp=2.0;

    fhJetPtArea = new TH2F("fhJetPtArea","Jet Area Distribution",fPtBins,fPtLow,fPtUp,JetPtAreaBins,JetPtAreaLow,JetPtAreaUp);
    fhJetPtArea->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fhJetPtArea->GetYaxis()->SetTitle("A_{jet}");
    fhJetPtArea->GetZaxis()->SetTitle("1/N_{Events} dN/dA_{jet}dp_{T}");
    fhJetPtArea->Sumw2();

    // Jet pT vs Constituent pT
    fhJetConstituentPt = new TH2F("fhJetConstituentPt","Jet constituents p_{T} distribution",fPtBins,fPtLow,fPtUp,10*fPtBins,fPtLow,fPtUp);
    fhJetConstituentPt->GetXaxis()->SetTitle("Jet p_{T} (GeV/c)");
    fhJetConstituentPt->GetYaxis()->SetTitle("Constituent p_{T} (GeV/c)");
    fhJetConstituentPt->GetZaxis()->SetTitle("1/N_{Events} dN/dp_{T,jet}dp_{T,track}");
    fhJetConstituentPt->Sumw2();

    fhJetTracksPt = new TH2F("fhJetTracksPt","Jet constituents Tracks p_{T} distribution",fPtBins,fPtLow,fPtUp,10*fPtBins,fPtLow,fPtUp);
    fhJetTracksPt->GetXaxis()->SetTitle("Jet p_{T} (GeV/c)");
    fhJetTracksPt->GetYaxis()->SetTitle("Constituent Track p_{T} (GeV/c)");
    fhJetTracksPt->GetZaxis()->SetTitle("1/N_{Events} dN/dp_{T,jet}dp_{T,track}");
    fhJetTracksPt->Sumw2();

    fhJetClustersPt = new TH2F("fhJetClustersPt","Jet constituents Clusters p_{T} distribution",fPtBins,fPtLow,fPtUp,10*fPtBins,fPtLow,fPtUp);
    fhJetClustersPt->GetXaxis()->SetTitle("Jet p_{T} (GeV/c)");
    fhJetClustersPt->GetYaxis()->SetTitle("Constituent Cluster p_{T} (GeV/c)");
    fhJetClustersPt->GetZaxis()->SetTitle("1/N_{Events} dN/dp_{T,jet}dp_{T,cluster}");
    fhJetClustersPt->Sumw2();

    // Jet pT vs Constituent Counts
    fhJetConstituentCounts = new TH2F("fhJetConstituentCounts","Jet constituents distribution",fPtBins,fPtLow,fPtUp,TCBins,0,(Double_t)TCBins);
    fhJetConstituentCounts->GetXaxis()->SetTitle("Jet p_{T} (GeV/c)");
    fhJetConstituentCounts->GetYaxis()->SetTitle("Constituent Count");
    fhJetConstituentCounts->GetZaxis()->SetTitle("1/N_{Events} dN/dp_{T,jet}dN_{constituent}");
    fhJetConstituentCounts->Sumw2();
    
    fhJetTracksCounts = new TH2F("fhJetTracksCounts","Jet constituents Tracks distribution",fPtBins,fPtLow,fPtUp,TCBins,0,(Double_t)TCBins);
    fhJetTracksCounts->GetXaxis()->SetTitle("Jet p_{T} (GeV/c)");
    fhJetTracksCounts->GetYaxis()->SetTitle("Constituent Track Count");
    fhJetTracksCounts->GetZaxis()->SetTitle("1/N_{Events} dN/dp_{T,jet}dN_{track}");
    fhJetTracksCounts->Sumw2();
    
    fhJetClustersCounts = new TH2F("fhJetClustersCounts","Jet constituents Clusters distribution",fPtBins,fPtLow,fPtUp,TCBins,0,(Double_t)TCBins);
    fhJetClustersCounts->GetXaxis()->SetTitle("Jet p_{T} (GeV/c)");
    fhJetClustersCounts->GetYaxis()->SetTitle("Constituent Cluster Count");
    fhJetClustersCounts->GetZaxis()->SetTitle("1/N_{Events} dN/dp_{T,jet}dN_{cluster}");
    fhJetClustersCounts->Sumw2();

    // Jet pT vs z_{constituent/track/cluster}
    fhJetPtZConstituent = new TH2F("fhJetPtZConstituent","Jet z_{constituent} distribution",fPtBins,fPtLow,fPtUp,TCBins,0,1.0);
    fhJetPtZConstituent->GetXaxis()->SetTitle("Jet p_{T} (GeV/c)");
    fhJetPtZConstituent->GetYaxis()->SetTitle("z_{constituent}");
    fhJetPtZConstituent->GetZaxis()->SetTitle("1/N_{Events} dN_{jet}/dp_{T,jet}dz_{constituent}");
    fhJetPtZConstituent->Sumw2();

    fhJetPtZTrack = new TH2F("fhJetPtZTrack","Jet z_{track} distribution",fPtBins,fPtLow,fPtUp,TCBins,0,1.0);
    fhJetPtZTrack->GetXaxis()->SetTitle("Jet p_{T} (GeV/c)");
    fhJetPtZTrack->GetYaxis()->SetTitle("z_{track}");
    fhJetPtZTrack->GetZaxis()->SetTitle("1/N_{Events} dN_{jet}/dp_{T,jet}dz_{track}");
    fhJetPtZTrack->Sumw2();
    
    fhJetPtZCluster = new TH2F("fhJetPtZCluster","Jet z_{cluster} distribution",fPtBins,fPtLow,fPtUp,TCBins,0,1.0);
    fhJetPtZCluster->GetXaxis()->SetTitle("Jet p_{T} (GeV/c)");
    fhJetPtZCluster->GetYaxis()->SetTitle("z_{cluster}");
    fhJetPtZCluster->GetZaxis()->SetTitle("1/N_{Events} dN_{jet}/dp_{T,jet}dz_{cluster}");
    fhJetPtZCluster->Sumw2();

    // Jet pT vs z_Leading{constituent/track/cluster}
    fhJetPtZLeadingConstituent = new TH2F("fhJetPtZLeadingConstituent","Jet z_{Leading,constituent} distribution",fPtBins,fPtLow,fPtUp,TCBins,0,1.0);
    fhJetPtZLeadingConstituent->GetXaxis()->SetTitle("Jet p_{T} (GeV/c)");
    fhJetPtZLeadingConstituent->GetYaxis()->SetTitle("z_{Leading,constituent}");
    fhJetPtZLeadingConstituent->GetZaxis()->SetTitle("1/N_{Events} dN_{jet}/dp_{T,jet}dz_{constituent}");
    fhJetPtZLeadingConstituent->Sumw2();
    
    fhJetPtZLeadingTrack = new TH2F("fhJetPtZLeadingTrack","Jet z_{Leading,track} distribution",fPtBins,fPtLow,fPtUp,TCBins,0,1.0);
    fhJetPtZLeadingTrack->GetXaxis()->SetTitle("Jet p_{T} (GeV/c)");
    fhJetPtZLeadingTrack->GetYaxis()->SetTitle("z_{Leading,track}");
    fhJetPtZLeadingTrack->GetZaxis()->SetTitle("1/N_{Events} dN_{jet}/dp_{T,jet}dz_{track}");
    fhJetPtZLeadingTrack->Sumw2();
    
    fhJetPtZLeadingCluster = new TH2F("fhJetPtZLeadingCluster","Jet z_{Leading,cluster} distribution",fPtBins,fPtLow,fPtUp,TCBins,0,1.0);
    fhJetPtZLeadingCluster->GetXaxis()->SetTitle("Jet p_{T} (GeV/c)");
    fhJetPtZLeadingCluster->GetYaxis()->SetTitle("z_{Leading,cluster}");
    fhJetPtZLeadingCluster->GetZaxis()->SetTitle("1/N_{Events} dN_{jet}/dp_{T,jet}dz_{cluster}");
    fhJetPtZLeadingCluster->Sumw2();
    
    // Event Centralities vs Leading Jet Pt
    fhEventCentralityVsZNA = new TH2F("fhEventCentralityVsZNA",EventCentralityString,fCentralityBins,fCentralityLow,fCentralityUp,fCentralityBins,fCentralityLow,fCentralityUp);
    fhEventCentralityVsZNA->GetXaxis()->SetTitle(Form("%s",CentralityString.Data()));
    fhEventCentralityVsZNA->GetYaxis()->SetTitle("Centrality (ZNA)");
    fhEventCentralityVsZNA->GetZaxis()->SetTitle("Probability Density");
    fhEventCentralityVsZNA->Sumw2();

    // 3D Histos
    if (fDo3DHistos == kTRUE)
    {
        // Jet pT, Eta, Phi
        fhJetPtEtaPhi = new TH3F("fhJetPtEtaPhi","Jet p_{T} vs #eta-#varphi",fPtBins,fPtLow,fPtUp,TCBins,fEMCalEtaMin,fEMCalEtaMax,TCBins,fEMCalPhiMin,fEMCalPhiMax);
        fhJetPtEtaPhi->GetXaxis()->SetTitle("p_{T} (GeV/c)");
        fhJetPtEtaPhi->GetYaxis()->SetTitle("#eta");
        fhJetPtEtaPhi->GetZaxis()->SetTitle("#varphi");
        fhJetPtEtaPhi->Sumw2();

        EventCentralityString = Form("%s vs ZNA Centrality vs Leading Jet p_{T} ",fCentralityTag.Data());
        fhEventCentralityVsZNAPt = new TH3F("fhEventCentralityVsZNAPt",EventCentralityString,fCentralityBins,fCentralityLow,fCentralityUp,fCentralityBins,fCentralityLow,fCentralityUp,fPtBins,fPtLow,fPtUp);
        fhEventCentralityVsZNAPt->GetXaxis()->SetTitle(Form("%s",CentralityString.Data()));
        fhEventCentralityVsZNAPt->GetYaxis()->SetTitle("Centrality (ZNA)");
        fhEventCentralityVsZNAPt->GetZaxis()->SetTitle("Leading Jet p_{T} (GeV/c)");
        fhEventCentralityVsZNAPt->Sumw2();
        
        fOutput->Add(fhJetPtEtaPhi);
        fOutput->Add(fhEventCentralityVsZNAPt);
    }

    // Neutral Energy Fraction Histograms & QA
    if (fDoNEFQAPlots==kTRUE)
    {
        fNEFOutput = new TList();
        fNEFOutput->SetOwner();
        fNEFOutput->SetName("ListNEFQAPlots");
        
        fhJetPtNEF = new TH2F("fhJetPtNEF","Jet p_{T} vs NEF",fPtBins,fPtLow,fPtUp,fNEFBins,fNEFLow,fNEFUp);
        fhJetPtNEF->GetXaxis()->SetTitle("Jet p_{T} (GeV/c)");
        fhJetPtNEF->GetYaxis()->SetTitle("Neutral Energy Fraction");
        fhJetPtNEF->GetZaxis()->SetTitle("1/N_{Events} dN_{jet}/dp_{T}dNEF");
        fhJetPtNEF->Sumw2();
        
        SetNEFJetDimensions(10); // Order: nef,Jet Pt,Eta,Phi,Centrality,Constituent mult,Charged mult, Neutral mult, z_leading
        SetNEFClusterDimensions(11); // Order: nef,Jet Pt,Eta,Phi,Centrality,F_Cross,z_leading,t_cell,deltat_cell,Cell Count, E_Cluster

        Int_t dimJetBins[fnDimJet];
        Double_t dimJetLow[fnDimJet];
        Double_t dimJetUp[fnDimJet];

        Int_t dimClusterBins[fnDimCluster];
        Double_t dimClusterLow[fnDimCluster];
        Double_t dimClusterUp[fnDimCluster];

        // Establish dimensinality bin counts and bounds
        // NEF
        dimJetBins[0] = dimClusterBins[0] = fNEFBins;
        dimJetLow[0] = dimClusterLow[0] = fNEFLow;
        dimJetUp[0] = dimClusterUp[0] = fNEFUp;

        // Jet Pt
        dimJetBins[1] = dimClusterBins[1] = fPtBins;
        dimJetLow[1] = dimClusterLow[1] = fPtLow;
        dimJetUp[1] = dimClusterUp[1] = fPtUp;

        // Eta-Phi
        dimJetBins[2] = dimJetBins[3] = dimClusterBins[2] = dimClusterBins[3] = TCBins;
        dimJetLow[2] = dimClusterLow[2] = fEMCalEtaMin;
        dimJetUp[2] = dimClusterUp[2] = fEMCalEtaMax;
        dimJetLow[3] = dimClusterLow[3] = fEMCalPhiMin;
        dimJetUp[3] = dimClusterUp[3] = fEMCalPhiMax;
        
        // Centrality
        dimJetBins[4] = dimClusterBins[4] = fCentralityBins;
        dimJetLow[4] = dimClusterLow[4] = fCentralityLow;
        dimJetUp[4] = dimClusterUp[4] = fCentralityUp;
        
        // z_leading
        dimJetBins[5] = dimClusterBins[5] = TCBins;
        dimJetLow[5] = dimClusterLow[5] = 0.0;
        dimJetUp[5] = dimClusterUp[5] = 1.0;
        
        // Jets Constituent Multiplicity Info {Total,Charged,Neutral}
        for (i=6;i<9;i++)
        {
            dimJetBins[i] = TCBins;
            dimJetLow[i] = 0.0;
            dimJetUp[i] = (Double_t)TCBins;
        }
        
        // z_leading^track
        dimJetBins[9] = TCBins;
        dimJetLow[9] = 0.0;
        dimJetUp[9] = 1.0;

        // Cluster E
        dimClusterBins[6] = fPtBins;
        dimClusterLow[6] = fPtLow;
        dimClusterUp[6] = fPtUp;
        
        // Cluster F_Cross
        dimClusterBins[7] = TCBins;
        dimClusterLow[7] = 0.0;
        dimClusterUp[7] = 1.0;
        
        // Cluster t_cell
        dimClusterBins[8] = 400;
        dimClusterLow[8] = -2e-07;
        dimClusterUp[8] = 2e-07;

        // Cluster delta t_cell
        dimClusterBins[9] = 100;
        dimClusterLow[9] = 0.0;
        dimClusterUp[9] = 1e-07;
        
        // Cluster Cell Count
        dimClusterBins[10] = TCBins;
        dimClusterLow[10] = 0.0;
        dimClusterUp[10] = 100.0;

        if (fDoTHnSparse == kTRUE)
        {
            fhJetNEFSignalInfo = new THnSparseF("fhJetNEFSignalInfo","Signal Jet NEF Information Histogram",fnDimJet,dimJetBins,dimJetLow,dimJetUp);
            fhJetNEFSignalInfo->Sumw2();
            
            fhClusterNEFSignalInfo = new THnSparseF("fhClusterNEFSignalInfo","Signal Jet NEF Cluster Information Histogram",fnDimCluster,dimClusterBins,dimClusterLow,dimClusterUp);
            fhClusterNEFSignalInfo->Sumw2();
            
            // Cluster Shape QA
            fhClusterShapeAll = new TH1F("fhClusterShapeAll","Cluster Shape of all CaloClustersCorr",10*TCBins,0,10*TCBins);
            fhClusterShapeAll->GetXaxis()->SetTitle("Cells");
            fhClusterShapeAll->GetYaxis()->SetTitle("1/N_{Events} dN/dCells");
            fhClusterShapeAll->Sumw2();
            
            fhClusterPtCellAll = new TH2F("fhClusterPtCellAll","Cluster p_{T} vs Cluster Shape of all CaloClustersCorr",fPtBins,fPtLow,fPtUp,10*TCBins,0,10*TCBins);
            fhClusterPtCellAll->GetXaxis()->SetTitle("p_{T} (GeV/c)");
            fhClusterPtCellAll->GetYaxis()->SetTitle("Cells");
            fhClusterPtCellAll->GetZaxis()->SetTitle("1/N_{Events} dN/dp_{T}dCells");
            fhClusterPtCellAll->Sumw2();
            
            if (fDoNEFSignalOnly == kFALSE)
            {
                fhJetNEFInfo = new THnSparseF("fhJetNEFInfo","Jet NEF Information Histogram",fnDimJet,dimJetBins,dimJetLow,dimJetUp);
                fhJetNEFInfo->Sumw2();
                
                fhClusterNEFInfo = new THnSparseF("fhClusterNEFInfo","Jet NEF Cluster Information Histogram",fnDimCluster,dimClusterBins,dimClusterLow,dimClusterUp);
                fhClusterNEFInfo->Sumw2();

                fNEFOutput->Add(fhJetNEFInfo);
                fNEFOutput->Add(fhClusterNEFInfo);
            }
            fNEFOutput->Add(fhJetNEFSignalInfo);
            fNEFOutput->Add(fhClusterNEFSignalInfo);
            fNEFOutput->Add(fhClusterShapeAll);
            fNEFOutput->Add(fhClusterPtCellAll);
        }
        fNEFOutput->Add(fhJetPtNEF);
        fOutput->Add(fNEFOutput);
    }
    
    // Add Histos & Profiles to List
    fOutput->Add(fh020Rho);
    fOutput->Add(fh80100Rho);
    fOutput->Add(fhRho);
    fOutput->Add(fhRhoCen);
    fOutput->Add(fh020BSPt);
    fOutput->Add(fh80100BSPt);
    fOutput->Add(fhBSPt);
    fOutput->Add(fhBSPtCen);
    fOutput->Add(fh020BSPtSignal);
    fOutput->Add(fh80100BSPtSignal);
    fOutput->Add(fhBSPtSignal);
    fOutput->Add(fhBSPtCenSignal);
    fOutput->Add(fh020DeltaPt);
    fOutput->Add(fh80100DeltaPt);
    fOutput->Add(fhDeltaPt);
    fOutput->Add(fhDeltaPtCen);
    fOutput->Add(fh020DeltaPtSignal);
    fOutput->Add(fh80100DeltaPtSignal);
    fOutput->Add(fhDeltaPtSignal);
    fOutput->Add(fhDeltaPtCenSignal);
    fOutput->Add(fh020DeltaPtNColl);
    fOutput->Add(fh80100DeltaPtNColl);
    fOutput->Add(fhDeltaPtNColl);
    fOutput->Add(fhDeltaPtCenNColl);
    fOutput->Add(fh020BckgFlucPt);
    fOutput->Add(fh80100BckgFlucPt);
    fOutput->Add(fhBckgFlucPt);
    fOutput->Add(fhBckgFlucPtCen);
    fOutput->Add(fpRho);
    fOutput->Add(fpLJetRho);
    fOutput->Add(fhJetPtArea);
    fOutput->Add(fhJetConstituentPt);
    fOutput->Add(fhJetTracksPt);
    fOutput->Add(fhJetClustersPt);
    fOutput->Add(fhJetConstituentCounts);
    fOutput->Add(fhJetTracksCounts);
    fOutput->Add(fhJetClustersCounts);
    fOutput->Add(fhJetPtZConstituent);
    fOutput->Add(fhJetPtZTrack);
    fOutput->Add(fhJetPtZCluster);
    fOutput->Add(fhJetPtZLeadingConstituent);
    fOutput->Add(fhJetPtZLeadingTrack);
    fOutput->Add(fhJetPtZLeadingCluster);
    fOutput->Add(fhEventCentralityVsZNA);
}

void AliAnalysisTaskFullpAJets::AlipAJetHistos::SetName(const char *name)
{
    fName = name;
}

void AliAnalysisTaskFullpAJets::AlipAJetHistos::SetCentralityTag(TString name)
{
    fCentralityTag = name.Data();
}

void AliAnalysisTaskFullpAJets::AlipAJetHistos::SetCentralityRange(Int_t bins, Double_t low, Double_t up)
{
    fCentralityBins=bins;
    fCentralityLow=low;
    fCentralityUp=up;
}

void AliAnalysisTaskFullpAJets::AlipAJetHistos::SetPtRange(Int_t bins, Double_t low, Double_t up)
{
    fPtBins=bins;
    fPtLow=low;
    fPtUp=up;
}

void AliAnalysisTaskFullpAJets::AlipAJetHistos::SetRhoPtRange(Int_t bins, Double_t low, Double_t up)
{
    fRhoPtBins=bins;
    fRhoPtLow=low;
    fRhoPtUp=up;
}

void AliAnalysisTaskFullpAJets::AlipAJetHistos::SetDeltaPtRange(Int_t bins, Double_t low, Double_t up)
{
    fDeltaPtBins=bins;
    fDeltaPtLow=low;
    fDeltaPtUp=up;
}

void AliAnalysisTaskFullpAJets::AlipAJetHistos::SetBackgroundFluctuationsPtRange(Int_t bins, Double_t low, Double_t up)
{
    fBckgFlucPtBins=bins;
    fBckgFlucPtLow=low;
    fBckgFlucPtUp=up;
}

void AliAnalysisTaskFullpAJets::AlipAJetHistos::SetLeadingJetPtRange(Int_t bins, Double_t low, Double_t up)
{
    fLJetPtBins=bins;
    fLJetPtLow=low;
    fLJetPtUp=up;
}

void AliAnalysisTaskFullpAJets::AlipAJetHistos::SetLeadingChargedTrackPtRange(Int_t bins, Double_t low, Double_t up)
{
    fLChargedTrackPtBins=bins;
    fLChargedTrackPtLow=low;
    fLChargedTrackPtUp=up;
}

void AliAnalysisTaskFullpAJets::AlipAJetHistos::SetNEFRange(Int_t bins, Double_t low, Double_t up)
{
    fNEFBins=bins;
    fNEFLow=low;
    fNEFUp=up;
}

void AliAnalysisTaskFullpAJets::AlipAJetHistos::SetSignalTrackPtBias(Bool_t chargedBias)
{
    fSignalTrackBias = chargedBias;
}

void AliAnalysisTaskFullpAJets::AlipAJetHistos::SetNEFJetDimensions(Int_t n)
{
    fnDimJet = n;
}

void AliAnalysisTaskFullpAJets::AlipAJetHistos::SetNEFClusterDimensions(Int_t n)
{
    fnDimCluster = n;
}

void AliAnalysisTaskFullpAJets::AlipAJetHistos::SetRhoValue(Double_t value)
{
    fRhoValue = value;
}

void AliAnalysisTaskFullpAJets::AlipAJetHistos::DoTHnSparse(Bool_t doTHnSparse)
{
    fDoTHnSparse = doTHnSparse;
}

void AliAnalysisTaskFullpAJets::AlipAJetHistos::Do3DPlotting(Bool_t do3DPlotting)
{
    fDo3DHistos = do3DPlotting;
}

void AliAnalysisTaskFullpAJets::AlipAJetHistos::FillRho(Double_t eventCentrality, Double_t rho)
{
    SetRhoValue(rho);
    fhRho->Fill(rho);
    fhRhoCen->Fill(rho,eventCentrality);
    fpRho->Fill(eventCentrality,rho);
    
    if (eventCentrality<=20)
    {
        fh020Rho->Fill(rho);
    }
    else if (eventCentrality>=80)
    {
        fh80100Rho->Fill(rho);
    }
}

void AliAnalysisTaskFullpAJets::AlipAJetHistos::FillBSJS(Double_t eventCentrality, Double_t rho, Double_t signalCut, TClonesArray *jetList, Int_t *indexJetList, Int_t nIndexJetList)
{
    Int_t i;
    Double_t tempPt=0.0;
    Double_t tempChargedHighPt=0.0;
    
    for (i=0;i<nIndexJetList;i++)
    {
        AliEmcalJet *myJet = (AliEmcalJet*) jetList->At(indexJetList[i]);
        tempPt=myJet->Pt()-rho*myJet->Area();
        tempChargedHighPt = myJet->MaxTrackPt();
        
        fhBSPt->Fill(tempPt);
        fhBSPtCen->Fill(tempPt,eventCentrality);
        if (eventCentrality<=20)
        {
            fh020BSPt->Fill(tempPt);
        }
        else if (eventCentrality>=80)
        {
            fh80100BSPt->Fill(tempPt);
        }
        if (fSignalTrackBias==kTRUE)
        {
            if (tempChargedHighPt>=signalCut)
            {
                fhBSPtSignal->Fill(tempPt);
                fhBSPtCenSignal->Fill(tempPt,eventCentrality);
                if (eventCentrality<=20)
                {
                    fh020BSPtSignal->Fill(tempPt);
                }
                else if (eventCentrality>=80)
                {
                    fh80100BSPtSignal->Fill(tempPt);
                }
            }
        }
        else
        {
            if (tempPt>=signalCut)
            {
                fhBSPtSignal->Fill(tempPt);
                fhBSPtCenSignal->Fill(tempPt,eventCentrality);
                if (eventCentrality<=20)
                {
                    fh020BSPtSignal->Fill(tempPt);
                }
                else if (eventCentrality>=80)
                {
                    fh80100BSPtSignal->Fill(tempPt);
                }
            }
        }
        tempPt=0.0;
        tempChargedHighPt=0.0;
    }
}

void AliAnalysisTaskFullpAJets::AlipAJetHistos::FillDeltaPt(Double_t eventCentrality, Double_t rho, Double_t jetRadius, Double_t *RCArray, Int_t nRC)
{
    Int_t i;
    Double_t tempPt=0.0;
    
    for (i=0;i<nRC;i++)
    {
        tempPt=RCArray[i]-rho*TMath::Power(jetRadius,2);
        fhDeltaPt->Fill(tempPt);
        fhDeltaPtCen->Fill(tempPt,eventCentrality);
        if (eventCentrality<=20)
        {
            fh020DeltaPt->Fill(tempPt);
        }
        else if (eventCentrality>=80)
        {
            fh80100DeltaPt->Fill(tempPt);
        }
        tempPt=0.0;
    }
}

void AliAnalysisTaskFullpAJets::AlipAJetHistos::FillDeltaPtSignal(Double_t eventCentrality, Double_t rho, Double_t jetRadius, Double_t *RCArray, Int_t nRC)
{
    Int_t i;
    Double_t tempPt=0.0;
    
    for (i=0;i<nRC;i++)
    {
        tempPt=RCArray[i]-rho*TMath::Power(jetRadius,2);
        fhDeltaPtSignal->Fill(tempPt);
        fhDeltaPtCenSignal->Fill(tempPt,eventCentrality);
        if (eventCentrality<=20)
        {
            fh020DeltaPtSignal->Fill(tempPt);
        }
        else if (eventCentrality>=80)
        {
            fh80100DeltaPtSignal->Fill(tempPt);
        }
        tempPt=0.0;
    }
}

void AliAnalysisTaskFullpAJets::AlipAJetHistos::FillDeltaPtNColl(Double_t eventCentrality, Double_t rho, Double_t jetRadius, Double_t *RCArray, Int_t nRC)
{
    Int_t i;
    Double_t tempPt=0.0;
    
    for (i=0;i<nRC;i++)
    {
        tempPt=RCArray[i]-rho*TMath::Power(jetRadius,2);
        fhDeltaPtNColl->Fill(tempPt);
        fhDeltaPtCenNColl->Fill(tempPt,eventCentrality);
        if (eventCentrality<=20)
        {
            fh020DeltaPtNColl->Fill(tempPt);
        }
        else if (eventCentrality>=80)
        {
            fh80100DeltaPtNColl->Fill(tempPt);
        }
        tempPt=0.0;
    }
}

void AliAnalysisTaskFullpAJets::AlipAJetHistos::FillBackgroundFluctuations(Double_t eventCentrality, Double_t rho, Double_t jetRadius)
{
    Double_t tempPt=0.0;
    
    tempPt=rho*TMath::Power(jetRadius,2);
    fhBckgFlucPt->Fill(tempPt);
    fhBckgFlucPtCen->Fill(tempPt,eventCentrality);
    if (eventCentrality<=20)
    {
        fh020BckgFlucPt->Fill(tempPt);
    }
    else if (eventCentrality>=80)
    {
        fh80100BckgFlucPt->Fill(tempPt);
    }
}

void AliAnalysisTaskFullpAJets::AlipAJetHistos::FillLeadingJetPtRho(Double_t jetPt, Double_t rho)
{
    fpLJetRho->Fill(jetPt,rho);
}

void AliAnalysisTaskFullpAJets::AlipAJetHistos::DoNEFQAPlots(Bool_t doNEFAna)
{
    fDoNEFQAPlots = doNEFAna;
}

void AliAnalysisTaskFullpAJets::AlipAJetHistos::DoNEFSignalOnly(Bool_t doNEFSignalOnly)
{
    fDoNEFSignalOnly = doNEFSignalOnly;
}

void AliAnalysisTaskFullpAJets::AlipAJetHistos::DoNEFAnalysis(Double_t nefCut, Double_t signalCut, TClonesArray *jetList, Int_t *indexJetList, Int_t nIndexJetList, TObjArray *clusterList, TClonesArray *orgClusterList, AliVEvent *event, AliEMCALGeometry *geometry, AliEMCALRecoUtils *recoUtils, AliVCaloCells *cells)
{
    if (fDoNEFQAPlots==kFALSE)
    {
        return;
    }
    
    if (fDoTHnSparse == kFALSE)
    {
        Int_t i;
        
        Double_t nef=0.0;
        Double_t jetPt=0.0;
        Double_t tempChargedHighPt=0.0;

        for (i=0;i<nIndexJetList;i++)
        {
            AliEmcalJet *myJet = (AliEmcalJet*) jetList->At(indexJetList[i]);
            tempChargedHighPt = myJet->MaxTrackPt();
            nef=myJet->NEF();
            jetPt=myJet->Pt();
            
            if (fSignalTrackBias==kTRUE)
            {
                if (tempChargedHighPt>=signalCut && nef<=nefCut)
                {
                    fhJetPtNEF->Fill(jetPt,nef);
                }
            }
            else
            {
                if (jetPt>=signalCut && nef<=nefCut)
                {
                    fhJetPtNEF->Fill(jetPt,nef);
                }
            }
            
            nef=0.0;
            jetPt=0.0;
            tempChargedHighPt=0.0;
        }
        return;
    }

    Int_t i,j,k;
    
    Double_t valJet[fnDimJet];
    Double_t valCluster[fnDimJet];
    for (i=0;i<fnDimJet;i++)
    {
        valJet[i]=0.0;
        valCluster[i]=0.0;
    }
    
    Double_t nef=0.0;
    Double_t jetPt=0.0;
    Double_t tempChargedHighPt=0.0;
    Double_t eta=0.0;
    Double_t phi=0.0;
    Int_t totalMult=0;
    Int_t chargedMult=0;
    Int_t neutralMult=0;
    Int_t iSupMod = -1, absId = -1, ieta = -1, iphi = -1;
    Bool_t shared = kFALSE;
    
    Double_t zLeading=0.0;
    Double_t zLeadingTrack=0.0;
    Double_t ECluster=0.0;
    Double_t eSeed=0.0;
    Double_t tCell=0.0;
    Double_t eCross=0.0;
    Double_t FCross=0.0;
    
    Double_t lowTime=9.99e99;
    Double_t upTime=-9.99e99;
    Int_t tempCellID=0;
    Double_t tempCellTime=0.0;
    
    Double_t event_centrality = event->GetCentrality()->GetCentralityPercentile(fCentralityTag);
    valJet[4] = valCluster[4] = event_centrality;

    // First, do Jet QA
    for (i=0;i<nIndexJetList;i++)
    {
        AliEmcalJet *myJet = (AliEmcalJet*) jetList->At(indexJetList[i]);
        tempChargedHighPt = myJet->MaxTrackPt();
        nef=myJet->NEF();
        jetPt=myJet->Pt();
        eta=myJet->Eta();
        phi=myJet->Phi();
        totalMult=myJet->GetNumberOfConstituents();
        chargedMult = myJet->GetNumberOfTracks();
        neutralMult=myJet->GetNumberOfClusters();
        zLeading=myJet->MaxPartPt()/myJet->Pt();
        zLeadingTrack=myJet->MaxTrackPt()/myJet->Pt();
        
        valJet[0] = valCluster[0] = nef;
        valJet[1] = valCluster[1] = jetPt;
        valJet[2] = valCluster[2] = eta;
        valJet[3] = valCluster[3] = phi;
        valJet[5] = valCluster[5] = zLeading;
        valJet[9] = zLeadingTrack;
        
        valJet[6] = totalMult;
        valJet[7] = chargedMult;
        valJet[8] = neutralMult;

        // Supress filling of this histogram due to memory size of THnSparse when running over large datasets
        if (fDoNEFSignalOnly == kFALSE)
        {
            fhJetNEFInfo->Fill(valJet);
        }
        
        if (fSignalTrackBias==kTRUE)
        {
            if (tempChargedHighPt>=signalCut && nef<=nefCut)
            {
                fhJetNEFSignalInfo->Fill(valJet);
                fhJetPtNEF->Fill(jetPt,nef);
            }
        }
        else
        {
            if (jetPt>=signalCut && nef<=nefCut)
            {
                fhJetNEFSignalInfo->Fill(valJet);
                fhJetPtNEF->Fill(jetPt,nef);
            }
        }

        for (j=0;j<neutralMult;j++)
        {
            AliVCluster* vcluster = (AliVCluster*) orgClusterList->At(myJet->ClusterAt(j));
            ECluster = vcluster->E();
            recoUtils->GetMaxEnergyCell(geometry,cells,vcluster,absId,iSupMod,ieta,iphi,shared);
            eSeed = cells->GetCellAmplitude(absId);
            tCell = cells->GetCellTime(absId);
            eCross = recoUtils->GetECross(absId,tCell,cells,event->GetBunchCrossNumber());
            FCross = 1 - eCross/eSeed;
            
            // Obtain Delta t of Cluster
            lowTime=9.99e99;
            upTime=-9.99e99;
            for (k=0;k<vcluster->GetNCells();k++)
            {
                tempCellID=vcluster->GetCellAbsId(k);
                tempCellTime=cells->GetCellTime(tempCellID);
                if (tempCellTime>upTime)
                {
                    upTime=tempCellTime;
                }
                if (tempCellTime<lowTime)
                {
                    lowTime=tempCellTime;
                }
            }
            valCluster[6] = ECluster;
            valCluster[7] = FCross;
            valCluster[8] = tCell;
            valCluster[9] = upTime-lowTime;
            valCluster[10] = vcluster->GetNCells();
            
            if (fDoNEFSignalOnly == kFALSE)
            {
                fhClusterNEFInfo->Fill(valCluster);
            }

            if (fSignalTrackBias==kTRUE)
            {
                if (tempChargedHighPt>=signalCut && nef<=nefCut)
                {
                    fhClusterNEFSignalInfo->Fill(valCluster);
                }
            }
            else
            {
                if (myJet->Pt()>=signalCut && nef<=nefCut)
                {
                    fhClusterNEFSignalInfo->Fill(valCluster);
                }
            }
            tempCellID=0;
            tempCellTime=0.0;
            eSeed=0.0;
            tCell=0.0;
            eCross=0.0;
            FCross=0.0;
            iSupMod=-1,absId=-1,ieta=-1,iphi=-1;
        }

        nef=0.0;
        jetPt=0.0;
        tempChargedHighPt=0.0;
        eta=0.0;
        phi=0.0;
        totalMult=0;
        chargedMult=0;
        neutralMult=0;
        zLeading=0.0;
        zLeadingTrack=0.0;
        ECluster=0.0;
    }
    
    // Now do Cluster QA
    for (i=0;i<clusterList->GetEntries();i++)
    {
        AliVCluster *vcluster = (AliVCluster*) clusterList->At(i);
        fhClusterShapeAll->Fill(vcluster->GetNCells());
        fhClusterPtCellAll->Fill(vcluster->E(),vcluster->GetNCells());
    }
}

void AliAnalysisTaskFullpAJets::AlipAJetHistos::FillMiscJetStats(TClonesArray *jetList, Int_t *indexJetList, Int_t nIndexJetList, TClonesArray *trackList, TClonesArray *clusterList, Double_t *vertex)
{
    Int_t i,j;

    TLorentzVector *cluster_vec = new TLorentzVector;
    for (i=0;i<nIndexJetList;i++)
    {
        AliEmcalJet *myJet = (AliEmcalJet*) jetList->At(indexJetList[i]);

        if (fDo3DHistos == kTRUE)
        {
            fhJetPtEtaPhi->Fill(myJet->Pt(),myJet->Eta(),myJet->Phi());
        }
        fhJetPtArea->Fill(myJet->Pt(),myJet->Area());
        fhJetConstituentCounts->Fill(myJet->Pt(),myJet->GetNumberOfConstituents());
        fhJetTracksCounts->Fill(myJet->Pt(),myJet->GetNumberOfTracks());
        fhJetClustersCounts->Fill(myJet->Pt(),myJet->GetNumberOfClusters());
        fhJetPtZLeadingConstituent->Fill(myJet->Pt(),myJet->MaxPartPt()/myJet->Pt());
        fhJetPtZLeadingTrack->Fill(myJet->Pt(),myJet->MaxTrackPt()/myJet->Pt());
        fhJetPtZLeadingCluster->Fill(myJet->Pt(),myJet->MaxClusterPt()/myJet->Pt());
        for (j=0;j<myJet->GetNumberOfTracks();j++)
        {
            AliVTrack *vtrack = (AliVTrack*) myJet->TrackAt(j,trackList);
            fhJetConstituentPt->Fill(myJet->Pt(),vtrack->Pt());
            fhJetTracksPt->Fill(myJet->Pt(),vtrack->Pt());
            fhJetPtZTrack->Fill(myJet->Pt(),vtrack->Pt()/myJet->Pt());
            fhJetPtZConstituent->Fill(myJet->Pt(),vtrack->Pt()/myJet->Pt());
        }
        for (j=0;j<myJet->GetNumberOfClusters();j++)
        {
            AliVCluster *vcluster = (AliVCluster*) myJet->ClusterAt(j,clusterList);
            vcluster->GetMomentum(*cluster_vec,vertex);
            fhJetConstituentPt->Fill(myJet->Pt(),cluster_vec->Pt());
            fhJetClustersPt->Fill(myJet->Pt(),vcluster->E());
            fhJetPtZCluster->Fill(myJet->Pt(),cluster_vec->Pt()/myJet->Pt());
            fhJetPtZConstituent->Fill(myJet->Pt(),cluster_vec->Pt()/myJet->Pt());
        }
    }
    delete cluster_vec;
}

void AliAnalysisTaskFullpAJets::AlipAJetHistos::FillJetEventCentrality(Double_t leadingJetPt, AliVEvent *event)
{
    Double_t event_centrality = event->GetCentrality()->GetCentralityPercentile(fCentralityTag);
    Double_t event_centrality_ZNA = event->GetCentrality()->GetCentralityPercentile("ZNA");
    
    fhEventCentralityVsZNA->Fill(event_centrality,event_centrality_ZNA);
    if (fDo3DHistos == kTRUE)
    {
        fhEventCentralityVsZNAPt->Fill(event_centrality,event_centrality_ZNA,leadingJetPt);
    }
}

TList* AliAnalysisTaskFullpAJets::AlipAJetHistos::GetOutputHistos()
{
    return fOutput;
}

Double_t AliAnalysisTaskFullpAJets::AlipAJetHistos::GetRho()
{
    return fRhoValue;
}

