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
class TObjArray;
class TLorentzVector;
class AliESDtrackCuts;
class AliEmcalJet;
class AliEMCALGeometry;

#ifndef ALIANALYSISTASKSE_H
#include "AliAnalysisTaskSE.h"
#endif

class AliAnalysisTaskFullpAJets : public AliAnalysisTaskSE
{
    // AlipAJetData Helper Class
    class AlipAJetData
    {
    public:
        AlipAJetData();
        AlipAJetData(const char *name, Bool_t isFull, Int_t nEntries);
        virtual ~AlipAJetData();
        
        // User Defined Sub-Routines
        void InitializeJetData(TClonesArray *jetList, Int_t nEntries);
        
        // Setters
        void SetName(const char *name);
        void SetIsJetsFull(Bool_t isFull);
        void SetTotalEntries(Int_t nEntries);
        void SetTotalJets(Int_t nJets);
        void SetTotalSignalJets(Int_t nSignalJets);
        void SetSignalCut(Double_t Pt);
        void SetLeading(Int_t index, Double_t Pt);
        void SetSubLeading(Int_t index, Double_t Pt);
        void SetJetIndex(Int_t index, Int_t At);
        void SetSignalJetIndex(Int_t index, Int_t At);
        void SetIsJetInArray(Bool_t isInArray, Int_t At);
        void SetAreaCutFraction(Double_t areaFraction);
        void SetJetR(Double_t jetR);
        
        // Getters
        Int_t GetTotalEntries();
        Int_t GetTotalJets();
        Int_t GetTotalSignalJets();
        Double_t GetSignalCut();
        Int_t GetLeadingIndex();
        Double_t GetLeadingPt();
        Int_t GetSubLeadingIndex();
        Double_t GetSubLeadingPt();
        Int_t GetJetIndex(Int_t At);
        Int_t GetSignalJetIndex(Int_t At);
        Bool_t GetIsJetInArray(Int_t At);
        
        Int_t* GetJets() const {return fJetsIndex;}
        Int_t* GetSignalJets() const {return fJetsSCIndex;}
        
    private:
        
        // Variables
        const char *fName;
        Bool_t fIsJetsFull;
        Int_t fnTotal;
        Int_t fnJets;
        Int_t fnJetsSC;
        
        Double_t fJetR;
        Double_t fSignalPt;
        Double_t fAreaCutFrac;
        
        Int_t fPtMaxIndex;
        Double_t fPtMax;
        Int_t fPtSubLeadingIndex;
        Double_t fPtSubLeading;
        
        Int_t *fJetsIndex;
        Int_t *fJetsSCIndex;
        Bool_t *fIsJetInArray;
    };
    
    // AlipAJetHistos Helper Class
    class AlipAJetHistos
    {
    public:
        AlipAJetHistos();
        AlipAJetHistos(const char *name);
        AlipAJetHistos(const char *name, const char *centag);
        virtual ~AlipAJetHistos();
        
        // User Defined Sub-Routines
        void Init();
        void FillRho(Double_t eventCentrality, Double_t rho);
        void FillBSJS(Double_t eventCentrality, Double_t rho, Double_t signalCut, TClonesArray *jetList, Int_t *indexJetList, Int_t nIndexJetList);
        void FillDeltaPt(Double_t eventCentrality, Double_t rho, Double_t jetRadius, Double_t *RCArray, Int_t nRC);
        void FillDeltaPtSignal(Double_t eventCentrality, Double_t rho, Double_t jetRadius, Double_t *RCArray, Int_t nRC);
        void FillBackgroundFluctuations(Double_t eventCentrality, Double_t rho, Double_t jetRadius);
        void FillLeadingJetPtRho(Double_t jetPt, Double_t rho);
        
        // Setters
        void SetName(const char *name);
        void SetCentralityTag(const char *name);
        void SetCentralityRange(Int_t bins, Double_t low, Double_t up);
        void SetPtRange(Int_t bins, Double_t low, Double_t up);
        void SetRhoPtRange(Int_t bins, Double_t low, Double_t up);
        void SetDeltaPtRange(Int_t bins, Double_t low, Double_t up);
        void SetBackgroundFluctuationsPtRange(Int_t bins, Double_t low, Double_t up);
        void SetLeadingJetPtRange(Int_t bins, Double_t low, Double_t up);
        
        // User Defined Functions
        TList* GetOutputHistos();
        
    private:
        TList *fOutput; // Output list
        
        // Histograms
        // This set of Histograms is for filling the Rho Spectral distributions
        TH1D *fh020Rho; //!
        TH1D *fh80100Rho; //!
        TH1D *fhRho; //!
        TH2D *fhRhoCen; //!
        
        // This set of Histograms is for filling the Background Subtracted Jet Spectra
        TH1D *fh020BSPt; //!
        TH1D *fh80100BSPt; //!
        TH1D *fhBSPt; //!
        TH2D *fhBSPtCen; //!
        
        // This set of Histograms is for filling the Background Subtracted Signal Jet Spectra
        TH1D *fh020BSPtSignal; //!
        TH1D *fh80100BSPtSignal; //!
        TH1D *fhBSPtSignal; //!
        TH2D *fhBSPtCenSignal; //!
        
        // This set of Histograms is for filling Delta Pt where the RC are at least 2R away from the leading Signal
        TH1D *fh020DeltaPt; //!
        TH1D *fh80100DeltaPt; //!
        TH1D *fhDeltaPt; //!
        TH2D *fhDeltaPtCen; //!
        
        // This set of Histograms is for filling Delta Pt where the RC have to spatial restrictions
        TH1D *fh020DeltaPtSignal; //!
        TH1D *fh80100DeltaPtSignal; //!
        TH1D *fhDeltaPtSignal; //!
        TH2D *fhDeltaPtCenSignal; //!
        
        // This set of Histograms is for filling Background Fluctuations Spectra
        TH1D *fh020BckgFlucPt; //!
        TH1D *fh80100BckgFlucPt; //!
        TH1D *fhBckgFlucPt; //!
        TH2D *fhBckgFlucPtCen; //!
        
        // Profiles
        TProfile *fpRho; //!
        TProfile *fpLJetRho; //!
        
        // Variables
        const char *fName;
        const char *fCentralityTag;
        
        Int_t fCentralityBins;
        Double_t fCentralityLow;
        Double_t fCentralityUp;
        
        Int_t fPtBins;
        Double_t fPtLow;
        Double_t fPtUp;
        
        Int_t fRhoPtBins;
        Double_t fRhoPtLow;
        Double_t fRhoPtUp;
        
        Int_t fDeltaPtBins;
        Double_t fDeltaPtLow;
        Double_t fDeltaPtUp;
        
        Int_t fBckgFlucPtBins;
        Double_t fBckgFlucPtLow;
        Double_t fBckgFlucPtUp;
        
        Int_t fLJetPtBins;
        Double_t fLJetPtLow;
        Double_t fLJetPtUp;
    };

    // AliAnalysisTaskFullpAJets
    public:
    AliAnalysisTaskFullpAJets();
    AliAnalysisTaskFullpAJets(const char *name);
    virtual ~AliAnalysisTaskFullpAJets();
    
    virtual void UserCreateOutputObjects();
    virtual void UserExecOnce();
    virtual void UserExec(Option_t *option);
    virtual void Terminate(Option_t *);

    // User Defined Sub-Routines
    void TrackCuts();
    void ClusterCuts();
    void TrackHisto();
    void ClusterHisto();
    void InitChargedJets();
    void InitFullJets();
    void JetPtArea();
    void GenerateTPCRandomConesPt();
    void GenerateEMCalRandomConesPt();
    
    void EstimateChargedRho0();
    void EstimateChargedRho1();
    void EstimateChargedRho2();
    void EstimateChargedRhoN();
    void EstimateChargedRhoScale();
    void EstimateChargedRhokT();
    void EstimateChargedRhokTScale();
    void EstimateChargedRhoCMS();
    void EstimateChargedRhoCMSScale();
    
    void EstimateFullRho0();
    void EstimateFullRho1();
    void EstimateFullRho2();
    void EstimateFullRhoN();
    void EstimateFullRhoDijet();
    void EstimateFullRhokT();
    void EstimateFullRhoCMS();
    
    void JetPtFullProfile();
    void JetPtChargedProfile();
    void JetPtEtaProfile();
    void DeleteJetData(Bool_t EMCalOn);
    
    // User Defined Functions
    Bool_t IsDiJetEvent();
    Bool_t InsideRect(Double_t phi,Double_t phi_min,Double_t phi_max,Double_t eta,Double_t eta_min,Double_t eta_max);
    Bool_t IsInEMCal(Double_t phi,Double_t eta);
    Bool_t IsInEMCalFull(Double_t r,Double_t phi,Double_t eta);
    Bool_t IsInEMCalPart(Double_t r,Double_t phi,Double_t eta);
    Bool_t IsInTPCFull(Double_t r,Double_t phi,Double_t eta);
    Bool_t IsInTPC(Double_t r,Double_t phi,Double_t eta,Bool_t Complete);
    Bool_t IsJetOverlap(AliEmcalJet *jet1,AliEmcalJet *jet2,Bool_t EMCalOn);
    
    Double_t AreaWithinTPC(Double_t r,Double_t eta);
    Double_t AreaWithinEMCal(Double_t r,Double_t phi,Double_t eta);
    Double_t AreaEdge(Double_t r,Double_t z);
    Double_t AreaOverlap(Double_t r,Double_t x,Double_t y);
    Double_t TransverseArea(Double_t r,Double_t psi0,Double_t phi,Double_t eta);
    Double_t MedianRhokT(Double_t *pTkTEntries, Double_t *RhokTEntries, Int_t nEntries);
    
    // Used to set the R for the jet finders
    inline void SetRjet(Int_t r)
    {
        fRJET = r;
    };
    
    // Used to set the Centrality Tag
    inline void SetCentralityTag(const char *centag)
    {
        fCentralityTag = centag;
    };
    
    // Used to set apriori Scaling Factor
    inline void SetScaleFactor(Double_t factor)
    {
        fScaleFactor = factor;
    };
    
    // Used to set the minimum pT required to analyize picotracks
    inline void SetTrackPtCut(Double_t pT)
    {
        fTrackMinPt = pT;
    }

    // Used to set the minimum pT required to analyize picotracks
    inline void SetClusterPtCut(Double_t pT)
    {
        fClusterMinPt = pT;
    }
    /*
    inline Double_t GetScaleFactor()
    {
        if (fRhoCharged>0)
        {
            return fRhoFull/fRhoCharged;
        }
        return -1;
    };
    */
    private:
    TList *fOutput; // Output list
    
    TH1D *fhTrackPt;  //!
    TH1D *fhTrackEta;  //!
    TH1D *fhTrackPhi;  //!
    TH1D *fhClusterPt;  //!
    TH1D *fhClusterEta;  //!
    TH1D *fhClusterPhi;  //!
    TH1D *fhCentrality; //!
    TH1D *fhEMCalCellCounts;  //! Plots the distribution of cluster counts in the EMCal. Used to determine which cells are hot (if any...)

    TH2D *fhTrackEtaPhi;  //!
    TH2D *fhClusterEtaPhi; //!
    TH2D *fhJetPtArea; //! Jet Area distribution vs Pt
    TH2D *fhJetConstituentPt; //! Pt distribution of jet constituents
    TH2D *fhRhoScale;  //!
    
    TProfile *fpEMCalEventMult;  //!
    TProfile *fpTPCEventMult;  //!
    TProfile *fpRhoScale; //! Scale of rho_total/rho_charged event/event vs centrality
    
    TProfile **fpJetEtaProfile; //!
    TProfile **fpJetAbsEtaProfile; //!
    TProfile **fpChargedJetRProfile; //!
    TProfile **fpJetRProfile; //!

    TProfile2D *fpTrackPtProfile;  //!
    TProfile2D *fpClusterPtProfile;  //!
    
    TProfile3D **fpChargedJetEDProfile;  //! Profile of Charged Jet Energy Density as a function of Jet pT, jet Eta, and distance from jet center in bins of 10% centrality cuts. Jet profile must be fully contained within TPC
    TProfile3D **fpJetEDProfile;  //! Profile of Jet Energy Density as a function of Jet pT, jet Eta, and distance from jet center in bins of 10% centrality cuts. Jet profile must be fully contained within EMCal
    
    AlipAJetHistos *fTPCRawJets;  //!
    AlipAJetHistos *fEMCalRawJets;  //!
    
    AlipAJetHistos *fRhoFull0;  //!
    AlipAJetHistos *fRhoFull1;  //!
    AlipAJetHistos *fRhoFull2;  //!
    AlipAJetHistos *fRhoFullN;  //!
    AlipAJetHistos *fRhoFullDijet;  //!
    AlipAJetHistos *fRhoFullkT;  //!
    AlipAJetHistos *fRhoFullCMS;  //!

    AlipAJetHistos *fRhoCharged0;  //!
    AlipAJetHistos *fRhoCharged1;  //!
    AlipAJetHistos *fRhoCharged2;  //!
    AlipAJetHistos *fRhoChargedN;  //!
    AlipAJetHistos *fRhoChargedScale;  //!
    AlipAJetHistos *fRhoChargedkT;  //!
    AlipAJetHistos *fRhoChargedkTScale;  //!
    AlipAJetHistos *fRhoChargedCMS;  //!
    AlipAJetHistos *fRhoChargedCMSScale;  //!

    AlipAJetData *fTPCJet;  //!
    AlipAJetData *fTPCFullJet;  //!
    AlipAJetData *fTPCOnlyJet;  //!
    AlipAJetData *fTPCkTFullJet;  //!
    AlipAJetData *fEMCalJet;  //!
    AlipAJetData *fEMCalFullJet;  //!
    AlipAJetData *fEMCalPartJet;  //!
    AlipAJetData *fEMCalkTFullJet;  //!

    // Variables
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
    Double_t fJetRForRho;  // Required distance a track/cluster must be away from a jet for rho calculation
    Double_t fJetAreaCutFrac;  // Fudge factor for selecting on jets with threshold Pt or higher
    Double_t fJetAreaThreshold;
    Int_t fnEMCalCells;  // Total number of cells in the EMCal
    Double_t fScaleFactor;  // Scale Factor obtained from Megan/Rosi
    Double_t fTrackMinPt;
    Double_t fClusterMinPt;
    
    const char *fCentralityTag;
    Int_t fCentralityBins;
    Double_t fCentralityLow;
    Double_t fCentralityUp;
    Double_t fEventCentrality;
    
    Double_t fRhoFull;  // From Full Rho 0
    Double_t fRhoCharged;  // From Charged Rho 0
    
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
    Int_t fnKTChargedJets;
    Int_t fnBckgClusters;
    
    Double_t fTPCJetThreshold;
    Double_t fEMCalJetThreshold;
    
    Double_t fvertex[3];
    Double_t fVertexWindow;
    Double_t fVertexMaxR;
    
    // Dynamic Array variables
    TClonesArray *fOrgTracks; //!
    TClonesArray *fOrgClusters; //!
    TClonesArray *fmyAKTFullJets; //!
    TClonesArray *fmyAKTChargedJets; //!
    TClonesArray *fmyKTFullJets; //!
    TClonesArray *fmyKTChargedJets; //!
    
    TObjArray *fmyTracks;  //!
    TObjArray *fmyClusters; //!
    
    Double_t *fEMCalRCBckgFluc; //! Stores the pT of RC Background clusters in EMCal at least 2R away from Leading Signal
    Double_t *fTPCRCBckgFluc; //! Stores the pT of RC Background clusters in TPC at least 2R away from Leading Signal
    Double_t *fEMCalRCBckgFlucSignal; //! Stores the pT of RC Background clusters in EMCal with no spatial restrictions
    Double_t *fTPCRCBckgFlucSignal; //! Stores the pT of RC Background clusters in TPC with no spatial restrictionsl

    AliAnalysisTaskFullpAJets(const AliAnalysisTaskFullpAJets&); // not implemented
    AliAnalysisTaskFullpAJets& operator=(const AliAnalysisTaskFullpAJets&); // not implemented
    
    ClassDef(AliAnalysisTaskFullpAJets, 1); // example of analysis
};
#endif
