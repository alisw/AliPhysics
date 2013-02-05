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
    fhTrackEtaPhi(0),
    fhClusterEtaPhi(0),
    fhJetPtArea(0),
    fhRhoCenTotal(0),
    fhRhoCenNoLeading(0),
    fhRho1B(0),
    fhRho1C(0),
    fhRho2B(0),
    fhRho3(0),
    fhJetConstituentPt(0),
    fhJetTrigR1A(0),
    fpEventMult(0),
    fpRhoTotal(0),
    fpRhoNoLeading(0),
    fpRho1B(0),
    fpRho3(0),
    fpRhoScale(0),
    fpRhokT(0),
    fpJetPtRhoTotal(0),
    fpJetPtRhoNoLeading(0),
    fpTrackPtProfile(0),
    fpClusterPtProfile(0)
{
    // Dummy constructor ALWAYS needed for I/O.
    fpJetEtaProfile = new TProfile *[14];
    fpJetAbsEtaProfile = new TProfile *[14];
    fpChargedJetRProfile = new TProfile *[8];
    fpJetRProfile= new TProfile *[4];
    fpChargedJetEDProfile= new TProfile3D *[10];
    fpJetEDProfile= new TProfile3D *[10];
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
    fhTrackEtaPhi(0),
    fhClusterEtaPhi(0),
    fhJetPtArea(0),
    fhRhoCenTotal(0),
    fhRhoCenNoLeading(0),
    fhRho1B(0),
    fhRho1C(0),
    fhRho2B(0),
    fhRho3(0),
    fhJetConstituentPt(0),
    fhJetTrigR1A(0),
    fpEventMult(0),
    fpRhoTotal(0),
    fpRhoNoLeading(0),
    fpRho1B(0),
    fpRho3(0),
    fpRhoScale(0),
    fpRhokT(0),
    fpJetPtRhoTotal(0),
    fpJetPtRhoNoLeading(0),
    fpTrackPtProfile(0),
    fpClusterPtProfile(0)
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

    DefineOutput(1,TList::Class());                                            // for output list
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
    fJetR=(Double_t)fRJET/10.;
    
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
    
    // Create histograms
    Double_t A1_Ptlow=0.0;
    Double_t A1_Ptup=100.0;
    Double_t A1_Triglow=0.0;
    Double_t A1_Trigup=100.0;
    Double_t A1_Rlow=0.4;
    Double_t A1_Rup=2.0;
    
    Int_t jetPt_Emcalbins = 200;
    Double_t jetPt_Emcallow = 0.0;
    Double_t jetPt_Emcalup = 200.0;
    
    Int_t TCBins=100;
    
    fhTrackPt = new TH1D("fhTrackPt","p_{T} distribution of tracks in event",10*jetPt_Emcalbins,jetPt_Emcallow,jetPt_Emcalup);
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

    fhClusterPt = new TH1D("fhClusterPt","p_{T} distribution of clusters in event",10*jetPt_Emcalbins,jetPt_Emcallow,jetPt_Emcalup);
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

    fhBckgFluc = new TH1D("fhBckgFluc",Form("p_{T} distribution of Background Clusters in near central events at center of EMCal with R=%g",fJetR),jetPt_Emcalbins,jetPt_Emcallow,jetPt_Emcalup);
    fhBckgFluc->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fhBckgFluc->GetYaxis()->SetTitle("1/N_{Events} dN/dp_{T}");
    fhBckgFluc->Sumw2();

    fhBckgMult = new TH1D("fhBckgMult",Form("Multiplicity distribution of Background Clusters in near central events at center of EMCal with R=%g",fJetR),jetPt_Emcalbins,jetPt_Emcallow,jetPt_Emcalup);
    fhBckgMult->GetXaxis()->SetTitle("Multiplicity");
    fhBckgMult->GetYaxis()->SetTitle("1/N_{Events}");
    fhBckgMult->Sumw2();
    
    fhChargedJetPt = new TH1D("fhChargedJetPt","Charged Jet p_{T} distribution for reconstructed Jets",jetPt_Emcalbins, jetPt_Emcallow, jetPt_Emcalup);
    fhChargedJetPt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fhChargedJetPt->GetYaxis()->SetTitle("counts");
    fhChargedJetPt->Sumw2();
    
    fhChargedJetPtAreaCut = new TH1D("fhChargedJetPtAreaCut"," Charged Jet p_{T} distribution for reconstructed Jets with Standard Area Cut",jetPt_Emcalbins, jetPt_Emcallow, jetPt_Emcalup);
    fhChargedJetPtAreaCut->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fhChargedJetPtAreaCut->GetYaxis()->SetTitle("counts");
    fhChargedJetPtAreaCut->Sumw2();

    fhJetPtEMCal = new TH1D("fhJetPtEMCal","Jet p_{T} distribution for reconstructed Jets within the EMCal",jetPt_Emcalbins, jetPt_Emcallow, jetPt_Emcalup);
    fhJetPtEMCal->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fhJetPtEMCal->GetYaxis()->SetTitle("counts");
    fhJetPtEMCal->Sumw2();

    fhJetPtEMCalAreaCut = new TH1D("fhJetPtEMCalAreaCut","Jet p_{T} distribution for reconstructed Jets within the EMCal with Standard Area Cut",jetPt_Emcalbins, jetPt_Emcallow, jetPt_Emcalup);
    fhJetPtEMCalAreaCut->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fhJetPtEMCalAreaCut->GetYaxis()->SetTitle("counts");
    fhJetPtEMCalAreaCut->Sumw2();

    fhJetPtEMCalAreaCutSignal = new TH1D("fhJetPtEMCalAreaCutSignal","Jet p_{T} distribution for reconstructed Jets within the EMCal with Standard Area Cut and Signal Cut",jetPt_Emcalbins, jetPt_Emcallow, jetPt_Emcalup);
    fhJetPtEMCalAreaCutSignal->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fhJetPtEMCalAreaCutSignal->GetYaxis()->SetTitle("counts");
    fhJetPtEMCalAreaCutSignal->Sumw2();

    Int_t jetPtbins = 200;
    Double_t jetPtlow = 0.0;
    Double_t jetPtup = 200.0;
    
    fhJetPtTPC = new TH1D("fhJetPtTPC","Jet p_{T} distribution for reconstructed Jets",jetPtbins, jetPtlow, jetPtup);
    fhJetPtTPC->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fhJetPtTPC->GetYaxis()->SetTitle("counts");
    fhJetPtTPC->Sumw2();

    fhJetPtTPCAreaCut = new TH1D("fhJetPtTPCAreaCut","Jet p_{T} distribution for reconstructed Jets",jetPtbins, jetPtlow, jetPtup);
    fhJetPtTPCAreaCut->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fhJetPtTPCAreaCut->GetYaxis()->SetTitle("counts");
    fhJetPtTPCAreaCut->Sumw2();

    Int_t jetPtAreabins=200;
    Double_t jetPtArealow=0.0;
    Double_t jetPtAreaup=2.0;
    
    fhJetPtArea = new TH2D("fhJetPtArea","Jet Area distribution vs Jet_{p_{T}}",jetPtbins, jetPtlow,jetPtup,jetPtAreabins,jetPtArealow,jetPtAreaup);
    fhJetPtArea->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fhJetPtArea->GetYaxis()->SetTitle("Jet Area");
    fhJetPtArea->Sumw2();

    Int_t A1_Ptbins = 100;  // Allow Hi-Res. Should be 95 for bin width of 1 GeV
    Int_t A1_Trigbins = 100; // Trigger jet from 5-100 GeV
    Int_t A1_Rbins = 16;     // 0.8 to 2 in 0.1 bin width
  
    fhJetTrigR1A = new TH3D("fhJetTrigR1A","Jet p_{T} distribution for reconstructed Jets vs Trigger Jet p_{T} vs #DeltaR",A1_Ptbins, A1_Ptlow, A1_Ptup,A1_Trigbins, A1_Triglow, A1_Trigup,A1_Rbins, A1_Rlow, A1_Rup);
    fhJetTrigR1A->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fhJetTrigR1A->Sumw2();

    Int_t B1_Ptbins = 250;
    Double_t B1_Ptlow=-50.0;
    Double_t B1_Ptup=200.0;
    
    fhJetTPtRhoTotal = new TH1D("fhJetTPtRhoTotal","True Jet p_{T} distribution for reconstructed Jets in the EMCal",B1_Ptbins, B1_Ptlow, B1_Ptup);
    fhJetTPtRhoTotal->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fhJetTPtRhoTotal->Sumw2();

    fhJetTPtRhoTotalSignal = new TH1D("fhJetTPtRhoTotalSignal","True Jet p_{T} distribution for reconstructed Jets in the EMCal with Signal cut",B1_Ptbins, B1_Ptlow, B1_Ptup);
    fhJetTPtRhoTotalSignal->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fhJetTPtRhoTotalSignal->Sumw2();

    fhJetTPtRhoNoLeading = new TH1D("fhJetTPtRhoNoLeading","True Jet p_{T} distribution for reconstructed Jets in the EMCal",B1_Ptbins, B1_Ptlow, B1_Ptup);
    fhJetTPtRhoNoLeading->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fhJetTPtRhoNoLeading->Sumw2();

    fhJetTPtRhoNoLeadingSignal = new TH1D("fhJetTPtRhoNoLeadingSignal","True Jet p_{T} distribution for reconstructed Jets in the EMCal with Signal cut",B1_Ptbins, B1_Ptlow, B1_Ptup);
    fhJetTPtRhoNoLeadingSignal->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fhJetTPtRhoNoLeadingSignal->Sumw2();

    fhJetTPt1B = new TH1D("fhJetTPt1B","True Jet p_{T} distribution for reconstructed Jets in the EMCal",B1_Ptbins, B1_Ptlow, B1_Ptup);
    fhJetTPt1B->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fhJetTPt1B->Sumw2();

    fhJetTPt1BSignal = new TH1D("fhJetTPt1BSignal","True Jet p_{T} distribution for reconstructed Jets in the EMCal",B1_Ptbins, B1_Ptlow, B1_Ptup);
    fhJetTPt1BSignal->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fhJetTPt1BSignal->Sumw2();

    Int_t EMCal_bckg_Ptbins = 200;
    Double_t EMCal_bckg_Ptlow=0.0;
    Double_t EMCal_bckg_Ptup=200.0;
    
    fhEMCalBckg1B = new TH1D("fhEMCalBckg1B","Cluster p_{T} distribution for reconstructed Tracks and Calocluster at least 0.5 away from all jets in the EMCal with R=0.4",EMCal_bckg_Ptbins, EMCal_bckg_Ptlow, EMCal_bckg_Ptup);
    fhEMCalBckg1B->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fhEMCalBckg1B->Sumw2();
    
    Int_t C1_Ptbins = 200;
    Double_t C1_Ptlow=0.0;
    Double_t C1_Ptup=200.0;
    
    fhJetTPt1C = new TH1D("fhJetTPt1C","True Jet p_{T} distribution for reconstructed Jets in the EMCal",C1_Ptbins, C1_Ptlow, C1_Ptup);
    fhJetTPt1C->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fhJetTPt1C->Sumw2();

    Int_t C1_rhoPtbins = 500;
    Double_t C1_rhoPtlow=0.0;
    Double_t C1_rhoPtup=50.0;
    
    fhRho1B = new TH2D("fhRho1B","Background p_{T}/A_{EMCal} distribution for events in EMCal vs Centrality",C1_rhoPtbins, C1_rhoPtlow, C1_rhoPtup,fCentralityBins,fCentralityLow,fCentralityUp);
    fhRho1B->GetXaxis()->SetTitle("p_{T}/Area (GeV/c)");
    fhRho1B->Sumw2();

    fhRho1C = new TH2D("fhRho1C","Background p_{T}/A_{EMCal} distribution for events in EMCal vs Centrality",C1_rhoPtbins, C1_rhoPtlow, C1_rhoPtup,fCentralityBins,fCentralityLow,fCentralityUp);
    fhRho1C->GetXaxis()->SetTitle("p_{T}/Area (GeV/c)");
    fhRho1C->Sumw2();

    fhRhoCenTotal= new TH2D("fhRhoCenTotal","Background Density #rho",C1_rhoPtbins,C1_rhoPtlow,C1_rhoPtup,fCentralityBins,fCentralityLow,fCentralityUp);
    fhRhoCenTotal->GetXaxis()->SetTitle("p_{T}/Area (GeV/c)");
    fhRhoCenTotal->Sumw2();

    fhRhoCenNoLeading= new TH2D("fhRhoCenNoLeading","Background Density #rho without leading jet",C1_rhoPtbins,C1_rhoPtlow,C1_rhoPtup,fCentralityBins,fCentralityLow,fCentralityUp);
    fhRhoCenNoLeading->GetXaxis()->SetTitle("p_{T}/Area (GeV/c)");
    fhRhoCenNoLeading->Sumw2();

    fhEMCalBckg1C = new TH1D("fhEMCalBckg1C","Cluster p_{T} distribution for reconstructed Tracks and Calocluster in Transverse area with R=0.4",EMCal_bckg_Ptbins, EMCal_bckg_Ptlow, EMCal_bckg_Ptup);
    fhEMCalBckg1C->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fhEMCalBckg1C->Sumw2();

    Int_t Bckg_Ptbins=150;
    Double_t Bckg_Ptlow=-50.0;
    Double_t Bckg_Ptup=100.0;

    fhDeltaPtTotal = new TH1D("fhDeltaPtTotal","F(A) p_{T} distribution for Random Cones using total #rho",Bckg_Ptbins,Bckg_Ptlow,Bckg_Ptup);
    fhDeltaPtTotal->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fhDeltaPtTotal->Sumw2();

    fhDeltaPtNoLeading = new TH1D("fhDeltaPtNoLeading","F(A) p_{T} distribution for Random Cones using leading jet #rho",Bckg_Ptbins,Bckg_Ptlow,Bckg_Ptup);
    fhDeltaPtNoLeading->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fhDeltaPtNoLeading->Sumw2();

    fhDeltaPt1B = new TH1D("fhDeltaPt1B","F(A) p_{T} distribution for Random Cones using method 1B #rho",Bckg_Ptbins,Bckg_Ptlow,Bckg_Ptup);
    fhDeltaPt1B->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fhDeltaPt1B->Sumw2();
    
    Int_t DeltaRho_Ptbins=500;
    Double_t DeltaRho_Ptlow=0.0;
    Double_t DeltaRho_Ptup=50.0;
    
    fhDeltaRho01 = new TH1D("fhDeltaRho01","Event by Event Differential between #rho_{0} and #rho_{1} of signal jets",DeltaRho_Ptbins,DeltaRho_Ptlow,DeltaRho_Ptup);
    fhDeltaRho01->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fhDeltaRho01->Sumw2();
    
    fhEMCalJet2A = new TH1D("fhEMCalJet2A","Cluster p_{T} distribution for jets within EMCal from di-jets in TPC",A1_Ptbins,A1_Ptlow,A1_Ptup);
    fhEMCalJet2A->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fhEMCalJet2A->Sumw2();
    
    fhJetTPt2B = new TH1D("fhJetTPt2B","True Jet p_{T} distribution for reconstructed Jets in the EMCal with dijet Trigger in TPC",B1_Ptbins, B1_Ptlow, B1_Ptup);
    fhJetTPt2B->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fhJetTPt2B->Sumw2();
    
    fhEMCalBckg2B = new TH1D("fhEMCalBckg2B","Cluster p_{T} distribution for reconstructed Tracks and Calocluster with dijet Trigger in TPC with R=0.4",EMCal_bckg_Ptbins, EMCal_bckg_Ptlow, EMCal_bckg_Ptup);
    fhEMCalBckg2B->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fhEMCalBckg2B->Sumw2();

    fhRho2B = new TH2D("fhRho2B","Background p_{T}/A_{EMCal} distribution for events in EMCal vs Centrality",C1_rhoPtbins, C1_rhoPtlow, C1_rhoPtup,fCentralityBins,fCentralityLow,fCentralityUp);
    fhRho2B->GetXaxis()->SetTitle("p_{T}/Area (GeV/c)");
    fhRho2B->Sumw2();

    fhRho3 = new TH2D("fhRho3","Charged background p_{T}/A_{TPC} distribution for events in TPC vs Centrality",C1_rhoPtbins, C1_rhoPtlow, C1_rhoPtup,fCentralityBins,fCentralityLow,fCentralityUp);
    fhRho3->GetXaxis()->SetTitle("p_{T}/Area (GeV/c)");
    fhRho3->Sumw2();
    
    fhJetTPt3 = new TH1D("fhJetTPt3","True Charged jet p_{T} distribution for reconstructed Jets in the TPC",B1_Ptbins, B1_Ptlow, B1_Ptup);
    fhJetTPt3->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fhJetTPt3->Sumw2();

    fhJetConstituentPt= new TH2D("fhJetConstituentPt","Jet constituents p_{T} distribution",jetPtbins, jetPtlow, jetPtup,10*jetPtbins, jetPtlow, jetPtup);
    fhJetConstituentPt->GetXaxis()->SetTitle("Jet p_{T} (GeV/c)");
    fhJetConstituentPt->GetYaxis()->SetTitle("Constituent p_{T} (GeV/c)");
    fhJetConstituentPt->Sumw2();
    
    fhEMCalCellCounts = new TH1D("fhEMCalCellCounts","Distribtuion of cluster counts across the EMCal",fnEMCalCells,1,fnEMCalCells);
    fhEMCalCellCounts->GetXaxis()->SetTitle("Absoulute Cell Id");
    fhEMCalCellCounts->GetYaxis()->SetTitle("Coutns");
    fhEMCalCellCounts->Sumw2();
    
    fpEventMult = new TProfile("fpEventMult","Event Multiplcity vs Centrality",10*fCentralityBins,fCentralityLow,fCentralityUp);
    fpEventMult->GetXaxis()->SetTitle("Centrality (V0M)");
    fpEventMult->GetYaxis()->SetTitle("Multiplicity");

    fpRhoTotal = new TProfile("fpRhoTotal","Background Profile vs Centrality",10*fCentralityBins,fCentralityLow,fCentralityUp);
    fpRhoTotal->GetXaxis()->SetTitle("Centrality (V0M)");
    fpRhoTotal->GetYaxis()->SetTitle("#p_{T}/Area (GeV/c)");

    fpRhoNoLeading = new TProfile("fpRhoNoLeading","Background Profile vs Centrality",10*fCentralityBins,fCentralityLow,fCentralityUp);
    fpRhoNoLeading->GetXaxis()->SetTitle("Centrality (V0M)");
    fpRhoNoLeading->GetYaxis()->SetTitle("p_{T}/Area (GeV/c)");

    fpRho1B = new TProfile("fpRho1B","Background Profile vs Centrality",10*fCentralityBins,fCentralityLow,fCentralityUp);
    fpRho1B->GetXaxis()->SetTitle("Centrality (V0M)");
    fpRho1B->GetYaxis()->SetTitle("p_{T}/Area (GeV/c)");

    fpRho3 = new TProfile("fpRho3","Background Profile vs Centrality",10*fCentralityBins,fCentralityLow,fCentralityUp);
    fpRho3->GetXaxis()->SetTitle("Centrality (V0M)");
    fpRho3->GetYaxis()->SetTitle("p_{T}/Area (GeV/c)");

    fpRhoScale = new TProfile("fpRhoScale","Scaling Factor Profile vs Centrality",10*fCentralityBins,fCentralityLow,fCentralityUp);
    fpRhoScale->GetXaxis()->SetTitle("Centrality (V0M)");
    fpRhoScale->GetYaxis()->SetTitle("Scale Factor");
    
    fpJetPtRhoTotal = new TProfile("fpJetPtRhoTotal","#rho_{0} Profile vs Leading jet p_{T}",jetPtbins,jetPtlow,jetPtup);
    fpJetPtRhoTotal->GetXaxis()->SetTitle("Leading jet p_{T} (GeV/c)");
    fpJetPtRhoTotal->GetYaxis()->SetTitle("p_{T}/Area (GeV/c)");

    fpJetPtRhoNoLeading = new TProfile("fpJetPtRhoNoLeading","#rho_{1} Profile vs Leading jet p_{T}",jetPtbins,jetPtlow,jetPtup);
    fpJetPtRhoNoLeading->GetXaxis()->SetTitle("Leading jet p_{T} (GeV/c)");
    fpJetPtRhoNoLeading->GetYaxis()->SetTitle("p_{T}/Area (GeV/c)");

    fpRhokT = new TProfile("fpRhokT","Background Profile vs Centrality",10*fCentralityBins,fCentralityLow,fCentralityUp);
    fpRhokT->GetXaxis()->SetTitle("Centrality (V0M)");
    fpRhokT->GetYaxis()->SetTitle("p_{T}/Area (GeV/c)");

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
    
    fpTrackPtProfile = new TProfile2D("fpTrackPtProfile","2D Profile of track pT density throughout the TPC",TCBins,fTPCEtaMin,fTPCEtaMax,TCBins,fTPCPhiMin,fTPCPhiMax);
    fpTrackPtProfile->GetXaxis()->SetTitle("#eta");
    fpTrackPtProfile->GetYaxis()->SetTitle("#phi");
    fpTrackPtProfile->GetZaxis()->SetTitle("p_{T} density (GeV/Area)");

    fpClusterPtProfile = new TProfile2D("fpClusterPtProfile","2D Profile of cluster pT density throughout the EMCal",TCBins,fEMCalEtaMin,fEMCalEtaMax,TCBins,fEMCalPhiMin,fEMCalPhiMax);
    fpClusterPtProfile->GetXaxis()->SetTitle("#eta");
    fpClusterPtProfile->GetYaxis()->SetTitle("#phi");
    fpClusterPtProfile->GetZaxis()->SetTitle("p_{T} density (GeV/Area)");

    fOutput->Add(fhTrackPt);
    fOutput->Add(fhTrackEta);
    fOutput->Add(fhTrackPhi);
    fOutput->Add(fhTrackEtaPhi);
    fOutput->Add(fhClusterPt);
    fOutput->Add(fhClusterEta);
    fOutput->Add(fhClusterPhi);
    fOutput->Add(fhClusterEtaPhi);
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
    fOutput->Add(fhRhoCenTotal);
    fOutput->Add(fhRhoCenNoLeading);
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
    fOutput->Add(fpEventMult);
    fOutput->Add(fpRhoTotal);
    fOutput->Add(fpRhoNoLeading);
    fOutput->Add(fpRho1B);
    fOutput->Add(fpRho3);
    fOutput->Add(fpRhoScale);
    fOutput->Add(fpJetPtRhoTotal);
    fOutput->Add(fpJetPtRhoNoLeading);
    fOutput->Add(fpRhokT);
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
    
    /*
     // First, make sure there are no signal jets centered within 2R from the center of the EMCal
     for (i=0;i<fnJetsPtTotalCut;i++)
     {
     AliEmcalJet *myJet =(AliEmcalJet*) fmyAKTFullJets->At(fJetPtTotalCutID[i]);
     TLorentzVector *jet_vec= new TLorentzVector;
     myJet->GetMom(*jet_vec);
     if (dummy->DeltaR(*jet_vec)<=(2*fJetR))
     {
     return;
     }
     delete jet_vec;
     }
    */
    
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
    Int_t i;
    
    fnAKTChargedJets = fmyAKTChargedJets->GetEntries();
    fJetPtChargedCutID = new Int_t[fnAKTChargedJets+1];
    fnJetsChargedPtCut=0;
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
            if (myJet->Pt()>=fTPCJetThreshold)
            {
                fJetPtChargedCutID[fnJetsChargedPtCut]=i;
                fnJetsChargedPtCut++;
            }
        }
    }
}

void AliAnalysisTaskFullpAJets::InitFullJets()
{
    // Preliminary Jet Placement and Selection Cuts
    Int_t i;
    
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
            
            // Method 1C
            if (fInEMCalFull[i]==kTRUE)
            {
                // Fill Jet Pt Distribution
                fhJetPtEMCal->Fill(myJet->Pt());
                if (myJet->Area()>=fJetAreaThreshold)
                {
                    fhJetPtEMCalAreaCut->Fill(myJet->Pt());
                    if (myJet->Pt()>=fEMCalJetThreshold)
                    {
                        fhJetPtEMCalAreaCutSignal->Fill(myJet->Pt());
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
        fpRhokT->Fill(fEventCentrality,myJet->Pt()/myJet->Area());
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
    
    //  Fill histograms
    fhRhoCenTotal->Fill(EMCal_rho,fEventCentrality);
    fpRhoTotal->Fill(fEventCentrality,EMCal_rho);
    fpJetPtRhoTotal->Fill(fPtFullMax,EMCal_rho);
    FillFullCorrJetPt(fhJetTPtRhoTotal,EMCal_rho,kFALSE);
    FillFullCorrJetPt(fhJetTPtRhoTotalSignal,EMCal_rho,kTRUE);
    fDeltaRho01=EMCal_rho;
    
    // Fill Background fluctuation spectrum F(A)
    if (fEventCentrality<=20)
    {
        FillBckgFlucDeltaPt(fhDeltaPtTotal,EMCal_rho);
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
    fhRhoCenNoLeading->Fill(EMCal_rho,fEventCentrality);
    fpRhoNoLeading->Fill(fEventCentrality,EMCal_rho);
    fpJetPtRhoNoLeading->Fill(fPtFullMax,EMCal_rho);
    FillFullCorrJetPt(fhJetTPtRhoNoLeading,EMCal_rho,kFALSE);
    FillFullCorrJetPt(fhJetTPtRhoNoLeadingSignal,EMCal_rho,kTRUE);
    fDeltaRho01-=EMCal_rho;
    FillFullDeltaRho(fhDeltaRho01,fDeltaRho01,kTRUE);
    fDeltaRho01=0.0;
    
    // Fill Background fluctuation spectrum F(A)
    if (fEventCentrality<=20)
    {
        FillBckgFlucDeltaPt(fhDeltaPtNoLeading,EMCal_rho);
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
    
    // Fill Background fluctuation spectrum F(A)
    if (fEventCentrality<=20)
    {
        FillBckgFlucDeltaPt(fhDeltaPt1B,EMCal_rho);
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
        }
    }
}

void AliAnalysisTaskFullpAJets::Method2A()
{
    Int_t i;
    
    if (fChargedFullMatch==kTRUE)
    {
        for (i=0;i<fnAKTFullJets;i++)
        {
            if (fInEMCalFull[i]==kTRUE && i!= fLeadingJetID && i!=fBackJetID)
            {
                AliEmcalJet *myJet =(AliEmcalJet*) fmyAKTFullJets->At(i);
                fhEMCalJet2A->Fill(myJet->Pt());
            }
        }
    }
}

void AliAnalysisTaskFullpAJets::Method2B()
{
    Int_t i,j;
    Bool_t track_away_from_jet;
    Bool_t cluster_away_from_jet;
    Double_t E_tracks_total=0.;
    Double_t E_caloclusters_total=0.;
    Double_t EMCal_rho=0.;
    Double_t jet_area_total=0.;
    
    for (i=0;i<fnTracks;i++)
    {
        // First, check if track is in the EMCal!!
        AliVTrack* vtrack = (AliVTrack*) fmyTracks->At(i); // pointer to reconstructed to track
        if (IsInEMCal(vtrack->Phi(),vtrack->Eta())==kTRUE)
        {
            track_away_from_jet=kTRUE;
            j=0;
            TLorentzVector *track_vec = new TLorentzVector;
            track_vec->SetPtEtaPhiE(vtrack->Pt(),vtrack->Eta(),vtrack->Phi(),vtrack->E());
            while (track_away_from_jet==kTRUE && j<fnJetsPtTotalCut)
            {
                //cout<<"j="<<j<<endl<<endl<<endl;
                
                TLorentzVector *jet_vec= new TLorentzVector;
                AliEmcalJet *myJet =(AliEmcalJet*) fmyAKTFullJets->At(fJetPtTotalCutID[j]);
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
    
    // Next, sum all CaloClusters within the EMCal (obviously all clusters must be within EMCal!!) that are away from jet(s) above Pt Threshold
    
    for (i=0;i<fnClusters;i++)
    {
        AliVCluster* vcluster = (AliVCluster*) fmyClusters->At(i); // pointer to reconstructed to cluster
        
        cluster_away_from_jet=kTRUE;
        j=0;
        
        TLorentzVector *cluster_vec = new TLorentzVector;
        vcluster->GetMomentum(*cluster_vec,fvertex);
        
        while (cluster_away_from_jet==kTRUE && j<fnJetsPtTotalCut)
        {
            TLorentzVector *jet_vec= new TLorentzVector;
            AliEmcalJet *myJet =(AliEmcalJet*) fmyAKTFullJets->At(fJetPtTotalCutID[j]);
            myJet->GetMom(*jet_vec);
            if (cluster_vec->DeltaR(*jet_vec)<=fJetR)
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
    
    // Determine area of all Jets that are within the EMCal
    for (i=0;i<fnJetsPtTotalCut;i++)
    {
        AliEmcalJet *myJet = (AliEmcalJet*) fmyAKTFullJets->At(fJetPtTotalCutID[i]);
        jet_area_total+=AreaWithinEMCal(fJetR,myJet->Phi(),myJet->Eta());
    }
    
    //Calculate Rho
    EMCal_rho=(E_tracks_total+E_caloclusters_total)/(fEMCalArea-jet_area_total);
    
    //Fill Background Cluster Histogram
    fhEMCalBckg2B->Fill(EMCal_rho*TMath::Pi()*TMath::Power(fJetR,2));
    fhRho2B->Fill(EMCal_rho,fEventCentrality);
    
    //Fill "True" Jet Pt Spectrum
    for (i=0;i<fnAKTFullJets;i++)
    {
        if (fInEMCalFull[i]==kTRUE)
        {
            AliEmcalJet *myJet = (AliEmcalJet*) fmyAKTFullJets->At(i);
            if (myJet->Area()>=fJetAreaThreshold && myJet->Pt()>=fEMCalJetThreshold)
            {
                fhJetTPt2B->Fill(myJet->Pt()-(EMCal_rho*myJet->Area()));
            }
        }
    }
}

void AliAnalysisTaskFullpAJets::Method3(Bool_t EMCalOn)
{
    Int_t i,j;
    Bool_t track_away_from_jet;
    Double_t E_tracks_total=0.;
    Double_t TPC_rho=0.;
    Double_t jet_area_total=0.;
    
    // First, sum all tracks within the EMCal that are away from jet(s) above Pt Threshold
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
            }
        }
    }
    
    // Fill Background Scaling factor profile.
    if (EMCalOn==kTRUE && fRhoCharged!=0)
    {
        fpRhoScale->Fill(fEventCentrality,(fRhoTotal/fRhoCharged));
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
    // Determine if event contains a di-jet within the detector. One of the jets could be partially contained within the EMCal. Furthermore, not required for the entire jet to be within the detector (No Fiducial Cut)
    Int_t i,j;
    const Double_t dijet_delta_phi=(180/360.)*2*TMath::Pi();
    const Double_t dijet_phi_acceptance=0.5*(30/360.)*2*TMath::Pi(); //Input the total acceptance within the paraenthesis to be +/- dijet_phi_acceptance
    Double_t dummy_phi=0.;
    Double_t dummyR=0.;
    
    fLeadingJetID=-1;
    fBackJetID=-1;
    fChargedBackJetID=-1;
    fChargedFullMatch=kFALSE;
    
    if (fnJetsChargedPtCut>1)
    {
        AliEmcalJet *myhJet =(AliEmcalJet*) fmyAKTChargedJets->At(fPtChargedMaxID);
        j=0;
        while (j<fnJetsChargedPtCut)
        {
            AliEmcalJet *myJet =(AliEmcalJet*) fmyAKTChargedJets->At(fJetPtChargedCutID[j]);
            dummy_phi=TMath::Abs(myhJet->Phi() - myJet->Phi());
            if ((dummy_phi>(dijet_delta_phi-dijet_phi_acceptance)) && (dummy_phi<(dijet_delta_phi+dijet_phi_acceptance)))
            {
                fChargedBackJetID=fJetPtTotalCutID[j];
                fnDiJetEvents++;
                
                // Now perform mathing between charged dijets to full dijets
                // Matched full jet must contain at least 70% of the Pt of the charged jet and must be no futher than R_Jet away
                for (i=0;i<fnAKTFullJets;i++)
                {
                    AliEmcalJet *myTestJet =(AliEmcalJet*) fmyAKTFullJets->At(i);
                    dummyR=TMath::Sqrt(TMath::Power(myhJet->Phi()-myTestJet->Phi(),2)+TMath::Power(myhJet->Eta()-myTestJet->Eta(),2));
                    if (dummyR<fJetR && myTestJet->Pt()>(0.7*myhJet->Pt()))
                    {
                        fLeadingJetID=i;
                    }
                    dummyR=TMath::Sqrt(TMath::Power(myJet->Phi()-myTestJet->Phi(),2)+TMath::Power(myJet->Eta()-myTestJet->Eta(),2));
                    if (dummyR<fJetR && myTestJet->Pt()>(0.7*myJet->Pt()))
                    {
                        fBackJetID=i;
                    }
                }
                if (fLeadingJetID !=-1 && fBackJetID !=-1)
                {
                    fChargedFullMatch=kTRUE;
                }
                return kTRUE;
            }
            j++;
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
