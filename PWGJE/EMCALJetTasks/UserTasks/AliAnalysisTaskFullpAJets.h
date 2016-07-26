#ifndef AliAnalysisTaskFullpAJets_H
#define AliAnalysisTaskFullpAJets_H

class TString;
class TH1F;
class TH2F;
class TH3F;
class THnSparse;
class TList;
class TProfile;
class TProfile2D;
class TProfile3D;
class TClonesArray;
class TObjArray;
class TLorentzVector;
class AliESDtrackCuts;
class AliEmcalJet;
class AliVEvent;
class AliEMCALGeometry;
class AliEMCALRecoUtils;
class AliVCaloCells;
class AliPicoTrack;

#include "AliAnalysisTaskEmcalJet.h"
#include "AliAnalysisTaskEmcal.h"

class AliAnalysisTaskFullpAJets : public AliAnalysisTaskEmcalJet
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
        void SetNEF(Double_t nef);
        void SetSignalTrackPtBias(Bool_t chargedBias);
        
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
        Double_t GetJetMaxChargedPt(Int_t At);
        Double_t GetNEF();
        
        Int_t* GetJets() const {return fJetsIndex;}  //!
        Int_t* GetSignalJets() const {return fJetsSCIndex;}  //!
        
    private:
        
        // Variables
        const char *fName;  //!
        Bool_t fIsJetsFull;
        Int_t fnTotal;
        Int_t fnJets;
        Int_t fnJetsSC;
        
        Double_t fJetR;
        Double_t fSignalPt;
        Double_t fAreaCutFrac;
        Double_t fNEF;
        Bool_t fSignalTrackBias;
        
        Int_t fPtMaxIndex;
        Double_t fPtMax;
        Int_t fPtSubLeadingIndex;
        Double_t fPtSubLeading;
        
        Int_t *fJetsIndex;  //!
        Int_t *fJetsSCIndex;  //!
        Bool_t *fIsJetInArray;  //!
        Double_t *fJetMaxChargedPt; //!
    };
    
    // AlipAJetHistos Helper Class
    class AlipAJetHistos
    {
    public:
        AlipAJetHistos();
        AlipAJetHistos(const char *name);
        AlipAJetHistos(const char *name, TString centag, Bool_t doNEF = kFALSE);
        AlipAJetHistos(const char *name, TString centag, Bool_t doNEF, Bool_t doNEFSignalOnly, Bool_t doTHnSparse, Bool_t do3DPlotting);

        virtual ~AlipAJetHistos();
        
        // User Defined Sub-Routines
        void Init();
        void FillRho(Double_t eventCentrality, Double_t rho);
        void FillBSJS(Double_t eventCentrality, Double_t rho, Double_t signalCut, TClonesArray *jetList, Int_t *indexJetList, Int_t nIndexJetList);
        void FillDeltaPt(Double_t eventCentrality, Double_t rho, Double_t jetRadius, Double_t *RCArray, Int_t nRC);
        void FillDeltaPtSignal(Double_t eventCentrality, Double_t rho, Double_t jetRadius, Double_t *RCArray, Int_t nRC);
        void FillDeltaPtNColl(Double_t eventCentrality, Double_t rho, Double_t jetRadius, Double_t *RCArray, Int_t nRC);
        void FillBackgroundFluctuations(Double_t eventCentrality, Double_t rho, Double_t jetRadius);
        void FillLeadingJetPtRho(Double_t jetPt, Double_t rho);
        void DoNEFQAPlots(Bool_t doNEFAna);
        void DoNEFSignalOnly(Bool_t doNEFSignalOnly);
        void DoNEFAnalysis(Double_t nefCut, Double_t signalCut, TClonesArray *jetList, Int_t *indexJetList, Int_t nIndexJetList, TObjArray *clusterList, TClonesArray *orgClusterList, AliVEvent *event, AliEMCALGeometry *geometry, AliEMCALRecoUtils *recoUtils, AliVCaloCells *cells);
        void FillMiscJetStats(TClonesArray *jetList, Int_t *indexJetList, Int_t nIndexJetList, TClonesArray *trackList, TClonesArray *clusterList, Double_t *vertex);
        void FillJetEventCentrality(Double_t leadingJetPt, AliVEvent *event);
        
        // Setters
        void SetName(const char *name);
        void SetCentralityTag(TString name);
        void SetCentralityRange(Int_t bins, Double_t low, Double_t up);
        void SetPtRange(Int_t bins, Double_t low, Double_t up);
        void SetRhoPtRange(Int_t bins, Double_t low, Double_t up);
        void SetDeltaPtRange(Int_t bins, Double_t low, Double_t up);
        void SetBackgroundFluctuationsPtRange(Int_t bins, Double_t low, Double_t up);
        void SetLeadingJetPtRange(Int_t bins, Double_t low, Double_t up);
        void SetLeadingChargedTrackPtRange(Int_t bins, Double_t low, Double_t up);
        void SetNEFRange(Int_t bins, Double_t low, Double_t up);
        void SetSignalTrackPtBias(Bool_t chargedBias);
        void SetNEFJetDimensions(Int_t n);
        void SetNEFClusterDimensions(Int_t n);
        void SetRhoValue(Double_t value);
        void DoTHnSparse(Bool_t doTHnSparse);
        void Do3DPlotting(Bool_t do3DPlotting);
        
        // User Defined Functions
        TList* GetOutputHistos();  //!
        Double_t GetRho();
        
    private:
        TList *fOutput; //! Output list
        
        // Histograms
        // This set of Histograms is for filling the Rho Spectral distributions
        TH1F *fh020Rho; //!
        TH1F *fh80100Rho; //!
        TH1F *fhRho; //!
        TH2F *fhRhoCen; //!
        
        // This set of Histograms is for filling the Background Subtracted Jet Spectra
        TH1F *fh020BSPt; //!
        TH1F *fh80100BSPt; //!
        TH1F *fhBSPt; //!
        TH2F *fhBSPtCen; //!
        
        // This set of Histograms is for filling the Background Subtracted Signal Jet Spectra
        TH1F *fh020BSPtSignal; //!
        TH1F *fh80100BSPtSignal; //!
        TH1F *fhBSPtSignal; //!
        TH2F *fhBSPtCenSignal; //!
        
        // This set of Histograms is for filling Delta Pt where the RC are at least 2R away from the leading Signal
        TH1F *fh020DeltaPt; //!
        TH1F *fh80100DeltaPt; //!
        TH1F *fhDeltaPt; //!
        TH2F *fhDeltaPtCen; //!
        
        // This set of Histograms is for filling Delta Pt where the RC have to spatial restrictions
        TH1F *fh020DeltaPtSignal; //!
        TH1F *fh80100DeltaPtSignal; //!
        TH1F *fhDeltaPtSignal; //!
        TH2F *fhDeltaPtCenSignal; //!

        // This set of Histograms is for filling Delta Pt with NColl
        TH1F *fh020DeltaPtNColl; //!
        TH1F *fh80100DeltaPtNColl; //!
        TH1F *fhDeltaPtNColl; //!
        TH2F *fhDeltaPtCenNColl; //!
        
        // This set of Histograms is for filling Background Fluctuations Spectra
        TH1F *fh020BckgFlucPt; //!
        TH1F *fh80100BckgFlucPt; //!
        TH1F *fhBckgFlucPt; //!
        TH2F *fhBckgFlucPtCen; //!
        
        // Profiles
        TProfile *fpRho; //!
        TProfile *fpLJetRho; //!
        
        // Jet Profile
        TH3F *fhJetPtEtaPhi; //!
        TH2F *fhJetPtArea; //!
        TH2F *fhJetConstituentPt; //!
        TH2F *fhJetTracksPt; //!
        TH2F *fhJetClustersPt; //!
        TH2F *fhJetConstituentCounts; //!
        TH2F *fhJetTracksCounts; //!
        TH2F *fhJetClustersCounts; //!
        TH2F *fhJetPtZConstituent; //!
        TH2F *fhJetPtZTrack; //!
        TH2F *fhJetPtZCluster; //!
        TH2F *fhJetPtZLeadingConstituent; //!
        TH2F *fhJetPtZLeadingTrack; //!
        TH2F *fhJetPtZLeadingCluster; //!
        
        // Event Centrality vs Leading Jet Pt
        TH2F *fhEventCentralityVsZNA; //!
        TH3F *fhEventCentralityVsZNAPt; //!
        
        // Histograms for Neutral Energy Fraction
        TList *fNEFOutput; //! NEF QA Plots
        
        TH2F *fhJetPtNEF; //!
        
        THnSparse *fhJetNEFInfo; //! Jet NEF Information Histogram
        THnSparse *fhJetNEFSignalInfo; //! Signal Jet NEF Information Histogram
        THnSparse *fhClusterNEFInfo; //! Cluster Jet NEF Information Histogram
        THnSparse *fhClusterNEFSignalInfo; //! Cluster Signal Jet NEF Information Histogram

        TH1F *fhClusterShapeAll; //!
        TH2F *fhClusterPtCellAll; //!

        // Variables
        const char *fName;  //!
        TString fCentralityTag; 
        
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
        
        Double_t fRhoValue;
        
        Int_t fLChargedTrackPtBins;
        Double_t fLChargedTrackPtLow;
        Double_t fLChargedTrackPtUp;
        
        Bool_t fDoNEFQAPlots;
        Bool_t fDoNEFSignalOnly;
        Bool_t fSignalTrackBias;
        Bool_t fDoTHnSparse;
        Bool_t fDo3DHistos;
        
        Int_t fNEFBins;
        Double_t fNEFLow;
        Double_t fNEFUp;
        
        Int_t fnDimJet;
        Int_t fnDimCluster;
        
        // These members are 'sourced' from the base class and are initalized in the constructor
        Double_t fEMCalPhiMin;
        Double_t fEMCalPhiMax;
        Double_t fEMCalEtaMin;
        Double_t fEMCalEtaMax;
    };

    // AliAnalysisTaskFullpAJets
    public:
    AliAnalysisTaskFullpAJets();
    AliAnalysisTaskFullpAJets(const char *name);
    virtual ~AliAnalysisTaskFullpAJets();
    
    virtual void UserCreateOutputObjects();
    virtual void UserExecOnce();

    virtual Bool_t Run();
    virtual void Terminate(Option_t *);

    // User Defined Sub-Routines
    void TrackCuts();
    void ClusterCuts();
    void EventCounts();
    void TrackHisto();
    void ClusterHisto();
    void InitChargedJets();
    void InitFullJets();
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
    
    void FullJetEnergyDensityProfile();
    void ChargedJetEnergyDensityProfile();
    
    void DeleteJetData(Int_t delOption);
    
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
    inline void SetCentralityTag(TString centag)
    {
        fCentralityTag = centag.Data();
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
    
    // Used to set Tracks collection name
    inline void SetTrackName(TString name)
    {
        fTrackName = name;
    };
    
    // Used to set Clusters collection name
    inline void SetClusterName(TString name)
    {
        fClusName = name;
    };

    // Used to set kT Charged Jet collection name
    inline void SetkTChargedJetName(TString name)
    {
        fkTChargedName = name;
    };

    // Used to set anti-kT Charged Jet collection name
    inline void SetAkTChargedJetName(TString name)
    {
        fAkTChargedName = name;
    };

    // Used to set kT Full Jet collection name
    inline void SetkTFullJetName(TString name)
    {
        fkTFullName = name;
    };
    
    // Used to set anti-kT Full Jet collection name
    inline void SetAkTFullJetName(TString name)
    {
        fAkTFullName = name;
    };
    
    // Used to set Ncoll for Delta Pt
    inline void SetNColl(Double_t ncoll)
    {
        fNColl = ncoll;
    };

    inline void SetNEFSignalJetCut(Double_t nef)
    {
        fNEFSignalJetCut = nef;
    };
    
    inline void DoNEFCalibration(Bool_t doNEF)
    {
        fDoNEF = doNEF;
    };
    
    inline void SetJetChargeBias(Bool_t trackBias)
    {
        fSignalTrackBias = trackBias;
    };
    
    inline void DoTrackQA(Bool_t doQA)
    {
        fTrackQA = doQA;
    };

    inline void DoClusterQA(Bool_t doQA)
    {
        fClusterQA = doQA;
    };
    
    inline void CalculateRhoJet(Int_t doRhoJet)
    {
        fCalculateRhoJet = doRhoJet;
    };
    
    inline void DoNEFSignalOnly(Bool_t doNEF)
    {
        fDoNEFSignalOnly = doNEF;
    };

    inline void DoVertexRCut(Bool_t doCut)
    {
        fDoVertexRCut = doCut;
    };
    
    inline void SetMCParticleLevel(Bool_t mcPartLevel)
    {
        fMCPartLevel = mcPartLevel;
    };
    
    inline void SetJetRAcceptance(Double_t r)
    {
        fJetRAccept = r;
    };
    
    inline void DoTHnSparse(Bool_t doTHnSparse)
    {
        fDoTHnSparse = doTHnSparse;
    };
    
    inline void DoJetRhoDensity(Bool_t doJetRhoDensity)
    {
        fDoJetRhoDensity = doJetRhoDensity;
    };
    
    inline void Do3DPlotting(Bool_t do3DPlotting)
    {
        fDo3DHistos = do3DPlotting;
    };
    private:
    TList *fOutput; //! Output list
    TList *flTrack; //! Track QA List
    TList *flCluster; //! Cluster QA List
    
    TH1F *fhTrackPt;  //!
    TH1F *fhTrackEta;  //!
    TH1F *fhTrackPhi;  //!
    TH1F *fhGlobalTrackPt;  //!
    TH1F *fhGlobalTrackEta;  //!
    TH1F *fhGlobalTrackPhi;  //!
    TH1F *fhComplementaryTrackPt;  //!
    TH1F *fhComplementaryTrackEta;  //!
    TH1F *fhComplementaryTrackPhi;  //!
    TH1F *fhClusterPt;  //!
    TH1F *fhClusterEta;  //!
    TH1F *fhClusterPhi;  //!
    TH1F *fhCentrality; //!
    TH1F *fhEMCalCellCounts;  //! Plots the distribution of cluster counts in the EMCal. Used to determine which cells are hot (if any...)
    
    TH1F *fhChargeAndNeutralEvents; //!
    TH1F *fhChargeOnlyEvents; //!
    TH1F *fhNeutralOnlyEvents; //!
    TH1F *fhNothingEvents; //!
    TH1F *fhEMCalChargeAndNeutralEvents; //!
    TH1F *fhEMCalChargeOnlyEvents; //!
    TH1F *fhEMCalNeutralOnlyEvents; //!
    TH1F *fhEMCalNothingEvents; //!

    TH2F *fhTrackEtaPhi;  //!
    TH2F *fhTrackPhiPt;  //!
    TH2F *fhTrackEtaPt;  //!
    TH2F *fhGlobalTrackEtaPhi;  //!
    TH2F *fhGlobalTrackPhiPt;  //!
    TH2F *fhGlobalTrackEtaPt;  //!
    TH2F *fhComplementaryTrackEtaPhi;  //!
    TH2F *fhComplementaryTrackPhiPt;  //!
    TH2F *fhComplementaryTrackEtaPt;  //!

    TH2F *fhClusterEtaPhi; //!
    TH2F *fhClusterPhiPt;  //!
    TH2F *fhClusterEtaPt;  //!
    
    TH2F *fhEMCalEventMult; //!
    TH2F *fhTPCEventMult; //!
    TH2F *fhEMCalTrackEventMult; //!
    
    TProfile *fpEMCalEventMult;  //!
    TProfile *fpTPCEventMult;  //!
    
    TProfile2D *fpTrackPtProfile;  //!
    TProfile2D *fpClusterPtProfile;  //!
    
    TProfile3D *fpFullJetEDProfile;  //!
    TProfile3D *fpChargedJetEDProfile;  //!
    TProfile3D *fpChargedJetEDProfileScaled;  //!
    
    AlipAJetHistos *fTPCRawJets;  //!
    AlipAJetHistos *fEMCalRawJets;  //!

    AlipAJetHistos *fRhoChargedCMSScale;  //!
    AlipAJetHistos *fRhoChargedScale;  //!
    
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
    AlipAJetHistos *fRhoChargedkT;  //!
    AlipAJetHistos *fRhoChargedkTScale;  //!
    AlipAJetHistos *fRhoChargedCMS;  //!

    AlipAJetData *fTPCJet;  //!
    AlipAJetData *fTPCFullJet;  //!
    AlipAJetData *fTPCOnlyJet;  //!
    AlipAJetData *fTPCJetUnbiased;  //!
    AlipAJetData *fTPCkTFullJet;  //!
    AlipAJetData *fEMCalJet;  //!
    AlipAJetData *fEMCalFullJet;  //!
    AlipAJetData *fEMCalPartJet;  //!
    AlipAJetData *fEMCalPartJetUnbiased;  //!
    AlipAJetData *fEMCalkTFullJet;  //!

    // Variables
    Bool_t fIsInitialized;
    Int_t fRJET;  // Used to set Anti_kt R. Called from AddTask Macro
    Long_t fnEvents;  // Counter for the number of events that made the physics selection with TPC+EMCal
    Long_t fnEventsCharged;  // Counter for the number of events that made the physics selection with TPC only
    Long_t fnDiJetEvents;  // Counter for the number of dijet events
    AliVEvent *fEvent;  //!
    AliEMCALRecoUtils *fRecoUtil;  //!
    AliEMCALGeometry *fEMCALGeometry;  //!
    AliVCaloCells *fCells;  //!
    Bool_t fDoNEF;
    Bool_t fDoNEFSignalOnly;
    Bool_t fSignalTrackBias;
    Bool_t fTrackQA;
    Bool_t fClusterQA;
    Int_t fCalculateRhoJet;
    Bool_t fDoVertexRCut;
    Bool_t fMCPartLevel;
    Bool_t fDoTHnSparse;
    Bool_t fDoJetRhoDensity;
    Bool_t fDo3DHistos;
    
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
    
    Double_t fParticlePtLow;
    Double_t fParticlePtUp;
    Int_t fParticlePtBins;

    Double_t fJetR;  // Jet resolution parameter
    Double_t fJetRAccept;  // Jet radius used for additional acceptance constraint. Used to restrict eta acceptance in full ana
    Double_t fFullEDJetR;  // Radius used to calculate boundaries for jet within EMCal
    Double_t fChargedEDJetR;  // Radius used to calculate boundaries (in eta) for jet within TPC
    Double_t fJetRForRho;  // Required distance a track/cluster must be away from a jet for rho calculation
    Double_t fJetAreaCutFrac;  // Fudge factor for selecting on jets with threshold Pt or higher
    Double_t fJetAreaThreshold;
    Int_t fnEMCalCells;  // Total number of cells in the EMCal
    Double_t fScaleFactor;  // Scale Factor obtained from Megan/Rosi
    Double_t fNColl;  // Used for partial rejection of signal from RC. Obtained via Glauber Calculations
    Double_t fTrackMinPt;
    Double_t fClusterMinPt;
    Double_t fNEFSignalJetCut;
    
    TString fCentralityTag; 
    Int_t fCentralityBins;
    Double_t fCentralityLow;
    Double_t fCentralityUp;
    Double_t fEventCentrality;
    
    Double_t fRhoFull;  // From Full Rho 0
    Double_t fRhoCharged;  // From Charged Rho 0
    
    // General Global variables
    Int_t fnTracks;
    Int_t fnEMCalTracks;
    Int_t fnClusters;
    Int_t fnCaloClusters;
    Int_t fnAKTFullJets;
    Int_t fnAKTChargedJets;
    Int_t fnKTFullJets;
    Int_t fnKTChargedJets;
    Int_t fnBckgClusters;
    
    Double_t fTPCJetThreshold;
    Double_t fEMCalJetThreshold;
    
    Double_t fVertex[3];
    Double_t fVertexWindow;
    Double_t fVertexMaxR;
    
    // Tracks and Clusters string names
    TString fTrackName;
    TString fClusName;
    
    // Jet Finder string names
    TString fkTChargedName;
    TString fAkTChargedName;
    TString fkTFullName;
    TString fAkTFullName;
    
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
    Double_t *fTPCRCBckgFlucSignal; //! Stores the pT of RC Background clusters in TPC with no spatial restrictions
    Double_t *fEMCalRCBckgFlucNColl; //! Stores the pT of RC Background clusters in EMCal with no spatial restrictions
    Double_t *fTPCRCBckgFlucNColl; //! Stores the pT of RC Background clusters in TPC with no spatial restrictions

    AliAnalysisTaskFullpAJets(const AliAnalysisTaskFullpAJets&); // not implemented
    AliAnalysisTaskFullpAJets& operator=(const AliAnalysisTaskFullpAJets&); // not implemented
    
    ClassDef(AliAnalysisTaskFullpAJets, 1); // example of analysis
};
#endif
