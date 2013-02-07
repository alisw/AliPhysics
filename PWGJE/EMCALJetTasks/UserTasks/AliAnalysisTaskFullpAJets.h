#ifndef AliAnalysisTaskFullpAJets_H
#define AliAnalysisTaskFullpAJets_H

class TH1D;
class TH2D;
class TH3D;
class TList;
class TProfile;
class TProfile2D;
class TProfile3D;
class TClonesArray;
class AliESDtrackCuts;

#ifndef ALIANALYSISTASKSE_H
#include "AliAnalysisTaskSE.h"
#endif

class AliAnalysisTaskFullpAJets : public AliAnalysisTaskSE
{
    public:
    AliAnalysisTaskFullpAJets();
    AliAnalysisTaskFullpAJets(const char *name);
    virtual ~AliAnalysisTaskFullpAJets();
    
    virtual void UserCreateOutputObjects();
    virtual void UserExecOnce();
    virtual void UserExec(Option_t *option);
    virtual void Terminate(Option_t *);

    // User Defined Sub-Routines
    void TrackHisto();
    void ClusterHisto();
    void EventHistos();
    void InitChargedJets();
    void InitFullJets();
    void EstimateTotalBackground();
    void EstimateBackgoundMinusLJet();
    void Method1A();
    void Method1B();
    void Method1C();
    void Method2A();
    void Method2B();
    void Method3(Bool_t EMCalOn);
    void JetPtFullProfile();
    void JetPtChargedProfile();
    void JetPtEtaProfile();
    void FillFullCorrJetPt(TH1D *myHisto,Double_t rho,Bool_t signal_cut);
    void FillFullCorrJetPt(TH2D *myHisto,Double_t rho, Bool_t signal_cut);
    void FillFullDeltaRho(TH1D *myHisto,Double_t delta_rho,Bool_t signal_cut);
    void FillBckgFlucDeltaPt(TH1D *myHisto, Double_t rho);
    void DeleteArrays(Bool_t EMCalOn);
    
    // User Defined Functions
    Bool_t IsDiJetEvent();
    Bool_t InsideRect(Double_t phi,Double_t phi_min,Double_t phi_max,Double_t eta,Double_t eta_min,Double_t eta_max);
    Bool_t IsInEMCal(Double_t phi,Double_t eta);
    Bool_t IsInEMCalFull(Double_t r,Double_t phi,Double_t eta);
    Bool_t IsInEMCalPart(Double_t r,Double_t phi,Double_t eta);
    Bool_t IsInTPCFull(Double_t r,Double_t phi,Double_t eta);
    Bool_t IsInTPC(Double_t r,Double_t phi,Double_t eta,Bool_t Complete);

    Double_t AreaWithinTPC(Double_t r,Double_t eta);
    Double_t AreaWithinEMCal(Double_t r,Double_t phi,Double_t eta);
    Double_t AreaEdge(Double_t r,Double_t z);
    Double_t AreaOverlap(Double_t r,Double_t x,Double_t y);
    Double_t TransverseArea(Double_t r,Double_t psi0,Double_t phi,Double_t eta);
    
    // Used to set the R for the jet finders
    inline void SetRjet(Int_t r)
    {
        fRJET = r;
    };

    private:
    TList *fOutput; // Output list
    AliESDtrackCuts *fTrackCuts; // Track cuts
    
    TH1D *fhTrackPt;  //!
    TH1D *fhTrackEta;  //!
    TH1D *fhTrackPhi;  //!
    TH1D *fhClusterPt;  //!
    TH1D *fhClusterEta;  //!
    TH1D *fhClusterPhi;  //!
    TH1D *fhCentrality; //!
    TH1D *fhBckgMult;  //!
    TH1D *fhBckgFluc;  //!
    TH1D *fhChargedJetPt; //! Charged Jet Pt distribution
    TH1D *fhChargedJetPtAreaCut; //! Charged Jet Pt distribution with standard Area cut applied
    TH1D *fhJetPtEMCal; //! Jet Pt distribution of Jets inside the EMCal
    TH1D *fhJetPtEMCalAreaCut; //! Jet Pt distribution of Jets inside the EMCal with standard Area cut applied
    TH1D *fhJetPtEMCalAreaCutSignal; //! Jet Pt distribution of Jets inside the EMCal with standard Area and Signal Threshold cut applied
    TH1D *fhJetPtTPC; //! Jet Pt distribution of Jet outside EMCal
    TH1D *fhJetPtTPCAreaCut; //! Jet Pt distribution of Jet outside EMCal with standard Area cut applied
    TH1D *fhJetTPtRhoTotal; //! Jet Pt distribution corrected for full background
    TH1D *fhJetTPtRhoTotalSignal; //! Jet Pt distribution corrected for full background with signal cut
    TH1D *fhJetTPtRhoNoLeading; //! Jet Pt distribution corrected for full background minus leading jet
    TH1D *fhJetTPtRhoNoLeadingSignal; //! Jet Pt distribution corrected for full background minus leading jet with signal cut
    TH1D *fhJetTPt1B; //! Jet Pt distribution corrected for background using Method 1B
    TH1D *fhJetTPt1BSignal; //! Jet Pt distribution corrected for background using Method 1B with signal cut
    TH1D *fhEMCalBckg1B; //! Cluster Pt distribution of Tracks+CaloCluster with R=0.4 Method 1B
    TH1D *fhJetTPt1C; //! Jet Pt distribution corrected for background using Method 1C
    TH1D *fhEMCalBckg1C; //! Cluster Pt distribution of Tracks+CaloCluster with R=0.4 Method 1C
    TH1D *fhEMCalJet2A; //! Clusters within the EMCal from di-jets within the TPC Satisfying certian criteria
    TH1D *fhJetTPt2B; //! Jet Pt distribution corrected for background using Method 2B
    TH1D *fhEMCalBckg2B; //! Cluster Pt distribution of Tracks+CaloCluster with R=0.4 Method 2B
    TH1D *fhJetTPt3; //! Charged jet Pt distribution corrected for background using Method 3
    TH1D *fhDeltaPtTotal;  //! Delta pT spectrum with total rho used
    TH1D *fhDeltaPtNoLeading;  //! Delta pT spectrum with total - leading jet rho used
    TH1D *fhDeltaPt1B;  //! Delta pT spectrum with all signal jets subtracted rho used
    TH1D *fhDeltaRho01;  //! Differential between rho_0 to rho_1 event by event
    TH1D *fhEMCalCellCounts;  //! Plots the distribution of cluster counts in the EMCal. Used to determine which cells are hot (if any...)
    TH1D *fh020RhoTotal; //! 0-20% Centrality rho plot for rho_0
    TH1D *fh020RhoNoLeading; //! 0-20% Centrality rho plot for rho_1
    TH1D *fh020Rho1B; //! 0-20% Centrality rho plot for rho_n
    TH1D *fh020Rho2B; //! 0-20% Centrality rho plot for di-jet rho
    TH1D *fh020Rho3; //! 0-20% Centrality rho plot for charged rho
    TH1D *fh020JetPtEMCal; //!
    TH1D *fh020JetPtEMCalAreaCut; //!
    TH1D *fh020JetPtEMCalAreaCutSignal; //!
    TH1D *fh020JetTPtRhoTotal;  //!
    TH1D *fh020JetTPtRhoTotalSignal;  //!
    TH1D *fh020JetTPtRhoNoLeading;  //!
    TH1D *fh020JetTPtRhoNoLeadingSignal;  //!
    TH1D *fh020JetTPt1B;  //!
    TH1D *fh020JetTPt1BSignal;  //!
    TH1D *fh020JetTPt1C;  //!
    TH1D *fh020JetTPt2B;  //!
    TH1D *fh020JetTPt3;  //!
    TH1D *fhDeltaPt2B;  //! Delta pT spectrum with Method 2B used for rho
    TH1D *fhDeltaPtkT;  //! Delta pT spectrum with kT jets used to calculate rho

    TH2D *fhTrackEtaPhi;  //!
    TH2D *fhClusterEtaPhi; //!
    TH2D *fhJetPtArea; //! Jet Area distribution vs Pt
    TH2D *fhRhoTotal;  //! Energy density of the EMCal (No jet exclusion)
    TH2D *fhRhoNoLeading;  //! Energy density of EMCal - leading fiducial jet
    TH2D *fhRho1B; //! Background estimate vs Centrality
    TH2D *fhRho1C; //! Background estimate vs Centrality
    TH2D *fhRho2B; //! Background estimate vs Centrality
    TH2D *fhRho3; //! Background estimate vs Centrality
    TH2D *fhJetConstituentPt; //! Pt distribution of jet constituents
    TH2D *fhJetPtCenEMCal;  //!
    TH2D *fhJetPtCenEMCalAreaCut;  //!
    TH2D *fhJetPtCenEMCalAreaCutSignal;  //!
    TH2D *fhJetTPtCenRhoTotal;  //!
    TH2D *fhJetTPtCenRhoTotalSignal;  //!
    TH2D *fhJetTPtCenRhoNoLeading;  //!
    TH2D *fhJetTPtCenRhoNoLeadingSignal;  //!
    TH2D *fhJetTPtCen1B;  //!
    TH2D *fhJetTPtCen1BSignal;  //!
    TH2D *fhJetTPtCen1C;  //!
    TH2D *fhJetTPtCen2B;  //!
    TH2D *fhJetTPtCen3;  //!
    
    TH3D *fhJetTrigR1A; //! Clusters from events with high Pt trigger as a funtion of trigger Pt and delta_R

    TProfile *fpEventMult;  //!
    TProfile *fpRhoTotal;  //!
    TProfile *fpRhoNoLeading;  //!
    TProfile *fpRho1B;  //!
    TProfile *fpRho2B;  //!
    TProfile *fpRho3;  //!
    TProfile *fpRhoScale; //! Scale of rho_total/rho_charged event/event vs centrality
    TProfile *fpRhokT;  //! Rho profile using rho from median kT jet
    TProfile *fpJetPtRhoTotal;  //!
    TProfile *fpJetPtRhoNoLeading;  //!
    
    TProfile **fpJetEtaProfile; //!
    TProfile **fpJetAbsEtaProfile; //!
    TProfile **fpChargedJetRProfile; //!
    TProfile **fpJetRProfile; //!

    TProfile2D *fpTrackPtProfile;  //!
    TProfile2D *fpClusterPtProfile;  //!
    
    TProfile3D **fpChargedJetEDProfile;  //! Profile of Charged Jet Energy Density as a function of Jet pT, jet Eta, and distance from jet center in bins of 10% centrality cuts. Jet profile must be fully contained within TPC
    TProfile3D **fpJetEDProfile;  //! Profile of Jet Energy Density as a function of Jet pT, jet Eta, and distance from jet center in bins of 10% centrality cuts. Jet profile must be fully contained within EMCal
    
    Bool_t fIsInitialized;
    Int_t fRJET;  // Used to set Anti_kt R. Called from AddTask Macro
    Long_t fnEvents;  // Counter for the number of events that made the physics selection with TPC+EMCal
    Long_t fnEventsCharged;  // Counter for the number of events that made the physics selection with TPC only
    Long_t fnDiJetEvents;  // Counter for the number of dijet events
    
    // Protected Global Variables
    Double_t fEMCalPhiMin;
    Double_t fEMCalPhiMax;
    Double_t fEMCalPhiTotal;
    Double_t fEMCalEtaMin;
    Double_t fEMCalEtaMax;
    Double_t fEMCalEtaTotal;
    Double_t fEMCalArea;
    
    Double_t fTPCPhiMin;
    Double_t fTPCPhiMax;
    Double_t fTPCPhiTotal;
    Double_t fTPCEtaMin;
    Double_t fTPCEtaMax;
    Double_t fTPCEtaTotal;
    Double_t fTPCArea;
    
    Double_t fJetR;
    Double_t fJetAreaCutFrac;  // Fudge factor for selecting on jets with threshold Pt or higher
    Double_t fJetAreaThreshold;
    Double_t fDeltaRho01;
    Int_t fnEMCalCells;  // Total number of cells in the EMCal
    
    Int_t fCentralityBins;
    Double_t fCentralityLow;
    Double_t fCentralityUp;
    Double_t fEventCentrality;
    
    Double_t fRhoTotal;
    Double_t fRhoCharged;
    
    // Jet profile variables
    Int_t fEtaProfileBins;
    Double_t fEtaProfileLow;
    Double_t fEtaProfileUp;

    Int_t fEDProfileRBins;
    Double_t fEDProfileRLow;
    Double_t fEDProfileRUp;
    
    Int_t fEDProfilePtBins;
    Double_t fEDProfilePtLow;
    Double_t fEDProfilePtUp;
    
    Int_t fEDProfileEtaBins;
    Double_t fEDProfileEtaLow;
    Double_t fEDProfileEtaUp;

    // General Global variables
    Int_t fnTracks;
    Int_t fnClusters;
    Int_t fnAKTFullJets;
    Int_t fnAKTChargedJets;
    Int_t fnKTFullJets;
    Int_t fnBckgClusters;
    
    Double_t fTPCJetThreshold;
    Double_t fEMCalJetThreshold;
    
    Double_t fvertex[3];
    Double_t fVertexWindow;
    Double_t fVertexMaxR;
    
    Int_t fnJetsPtCut;
    Int_t fnJetsPtTPCCut;
    Int_t fnJetsPtTotalCut;
    Int_t fnJetsChargedPtCut;
    Int_t fnJetskTEMCalFull;
    
    Int_t fPtMaxID;
    Int_t fPtFullMaxID;
    Int_t fPtTPCMaxID;
    Int_t fPtFullTPCMaxID;
    Int_t fPtTotalMaxID;
    Int_t fPtChargedMaxID;
    
    Double_t fPtMax;
    Double_t fPtFullMax;
    Double_t fPtTPCMax;
    Double_t fPtFullTPCMax;
    Double_t fPtTotalMax;
    Double_t fPtChargedMax;
    
    Int_t fChargedBackJetID;
    Bool_t fChargedFullMatch; // True if a match is found
    // These two variables are to mathes the IDs between the dijets from charged jet array to the corresponding Full jets array
    Int_t fLeadingJetID;
    Int_t fBackJetID;

    // Dynamic Array variables
    TClonesArray *fmyTracks; //!
    TClonesArray *fmyClusters; //!
    TClonesArray *fmyAKTFullJets; //!
    TClonesArray *fmyAKTChargedJets; //!
    TClonesArray *fmyKTFullJets; //!

    Int_t *fJetPtCutID; //!  Stores the jets(ID) above a Threshold Pt for EMCal
    Int_t *fJetPtTPCCutID;  //!  Stores the jets above a Threshold Pt for TPC
    Int_t *fJetPtTotalCutID;  //!  Stores the jets(ID) above a Threshold Pt
    Int_t *fJetPtChargedCutID; //!  Stores the jets(ID) above a Threshold Pt for TPC in events without the EMCal on 
    Int_t *fJetkTEMCalFullID;  //!
    
    Bool_t *fInEMCal; //!
    Bool_t *fInEMCalFull; //!
    Bool_t *fInTPCFull; //!
    Bool_t *fInTPCChargedFull; //!
    
    Double_t *fRCBckgFluc; //! Stores the pT of RC Background clusters in EMCal

    AliAnalysisTaskFullpAJets(const AliAnalysisTaskFullpAJets&); // not implemented
    AliAnalysisTaskFullpAJets& operator=(const AliAnalysisTaskFullpAJets&); // not implemented
    
    ClassDef(AliAnalysisTaskFullpAJets, 1); // example of analysis
};
#endif
