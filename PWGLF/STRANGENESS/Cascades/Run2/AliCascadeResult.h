#ifndef AliCascadeResult_H
#define AliCascadeResult_H
#include <TNamed.h>
#include <TList.h>
#include <TH3F.h>
#include <TProfile.h>
#include "AliVWeakResult.h"

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// TObject to hold Cascade configuration + results histogram
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class AliCascadeResult : public AliVWeakResult {
    
public:
    enum EMassHypo {
        kXiMinus    = 0,
        kXiPlus     = 1,
        kOmegaMinus = 2,
        kOmegaPlus  = 3
    };
    
    //Dummy Constructor
    AliCascadeResult();
    
    //Standard Constructor
    AliCascadeResult(const char * name, AliCascadeResult::EMassHypo lMassHypo = AliCascadeResult::kXiMinus, const char * title = "Cascade Result");
    
    //Variable-Binning Constructor:
    // Binning in ( centrality , momentum ) can be chosen and invariant mass is fixed at defaults 
    AliCascadeResult(const char * name, AliCascadeResult::EMassHypo lMassHypo, const char * title, Long_t lNCentBins, Double_t *lCentBins, Long_t lNPtBins, Double_t *lPtBins);
    
    //Settable centrality / momentum / invmass binning
    AliCascadeResult(const char * name, AliCascadeResult::EMassHypo lMassHypo, const char * title, Long_t lNCentBins, Double_t *lCentBins, Long_t lNPtBins, Double_t *lPtBins, Long_t lNMassBins, Double_t lMinMass, Double_t lMaxMass);
    
    //Specific uses
    AliCascadeResult(AliCascadeResult *lCopyMe, TString lNewName);
    AliCascadeResult(const AliCascadeResult& lCopyMe, TString lNewName);
    ~AliCascadeResult();
    
    void Clear(Option_t* = "") {}; //dummy
    
    AliCascadeResult& operator=(const AliCascadeResult& lCopyMe);
    
    Long64_t Merge(TCollection *hlist);
    
    //Acceptance
    void SetCutMinRapidity      ( Double_t lCut ) { fCutMinRapidity       = lCut; }
    void SetCutMaxRapidity      ( Double_t lCut ) { fCutMaxRapidity       = lCut; }
    
    //Setters for V0 Cuts
    void SetCutDCANegToPV     ( Double_t lCut ) { fCutDCANegToPV       = lCut; }
    void SetCutDCAPosToPV     ( Double_t lCut ) { fCutDCAPosToPV       = lCut; }
    void SetCutDCAV0Daughters ( Double_t lCut ) { fCutDCAV0Daughters   = lCut; }
    void SetCutV0CosPA        ( Double_t lCut ) { fCutV0CosPA          = lCut; }
    void SetCutV0Radius       ( Double_t lCut ) { fCutV0Radius         = lCut; }
    
    //Setters for Cascade Cuts
    void SetCutDCAV0ToPV        ( Double_t lCut ) { fCutDCAV0ToPV         = lCut; }
    void SetCutV0Mass           ( Double_t lCut ) { fCutV0Mass            = lCut; }
    void SetCutV0MassSigma      ( Double_t lCut ) { fCutV0MassSigma       = lCut; }
    void SetCutDCABachToPV      ( Double_t lCut ) { fCutDCABachToPV       = lCut; }
    void SetCutDCACascDaughters ( Double_t lCut ) { fCutDCACascDaughters  = lCut; }
    void SetCutCascCosPA        ( Double_t lCut ) { fCutCascCosPA         = lCut; }
    void SetCutCascRadius       ( Double_t lCut ) { fCutCascRadius        = lCut; }
    void SetCutDCABachToBaryon  ( Double_t lCut ) { fCutDCABachToBaryon   = lCut; }
    void SetCutBachBaryonCosPA  ( Double_t lCut ) { fCutBachBaryonCosPA   = lCut; }
    void SetCutMinV0Lifetime    ( Double_t lCut ) { fCutMinV0Lifetime     = lCut; }
    void SetCutMaxV0Lifetime    ( Double_t lCut ) { fCutMaxV0Lifetime     = lCut; }
    
    //Miscellaneous
    void SetCutProperLifetime        ( Double_t lCut ) { fCutProperLifetime        = lCut; }
    void SetCutTPCdEdx               ( Double_t lCut ) { fCutTPCdEdx               = lCut; }
    void SetCutXiRejection           ( Double_t lCut ) { fCutXiRejection           = lCut; }
    
    //MC specific
    void SetCutMCPhysicalPrimary    ( Bool_t lCut ) { fCutMCPhysicalPrimary    = lCut; }
    void SetCutMCPDGCodeAssociation ( Bool_t lCut ) { fCutMCPDGCodeAssociation = lCut; }
    void SetCutMCUseMCProperties    ( Bool_t lCut ) { fCutMCUseMCProperties    = lCut; }
    void SetCutMCSelectBump         ( Bool_t lCut ) { fCutMCSelectBump         = lCut; }
    
    //Rsn-like bg subtraction (experimental
    void SetSwapBachelorCharge         ( Bool_t lCut ) { fSwapBachCharge = lCut; }
    void SetSwapBaryon                 ( Bool_t lCut ) { fSwapBaryon        = lCut; }
    void SetSwapV0MesonCharge          ( Bool_t lCut ) { fSwapV0MesonCharge = lCut; }
    void SetSwapV0BaryonCharge         ( Bool_t lCut ) { fSwapV0BaryonCharge = lCut; }
    
    //Track Quality
    void SetCutUseITSRefitTracks    ( Bool_t lCut ) { fCutUseITSRefitTracks    = lCut; }
    void SetCutLeastNumberOfClusters ( Double_t lCut ) { fCutLeastNumberOfClusters = lCut; }
    void SetCutMinEtaTracks ( Double_t lCut ) { fCutMinEtaTracks = lCut; }
    void SetCutMaxEtaTracks ( Double_t lCut ) { fCutMaxEtaTracks = lCut; }
    void SetCutMaxChi2PerCluster ( Double_t lCut ) { fCutMaxChi2PerCluster = lCut; }
    void SetCutMinTrackLength    ( Double_t lCut ) { fCutMinTrackLength    = lCut; }
    
    //Variable CascCosPA
    void SetCutUseVarCascCosPA      ( Bool_t lCut )   { fCutUseVariableCascCosPA     = lCut; }
    void SetCutVarCascCosPAExp0Const( Double_t lCut ) { fCutVarCascCosPA_Exp0Const = lCut; }
    void SetCutVarCascCosPAExp0Slope( Double_t lCut ) { fCutVarCascCosPA_Exp0Slope = lCut; }
    void SetCutVarCascCosPAExp1Const( Double_t lCut ) { fCutVarCascCosPA_Exp1Const = lCut; }
    void SetCutVarCascCosPAExp1Slope( Double_t lCut ) { fCutVarCascCosPA_Exp1Slope = lCut; }
    void SetCutVarCascCosPAConst    ( Double_t lCut ) { fCutVarCascCosPA_Const     = lCut; }
    void SetCutVarCascCosPA ( Double_t l1, Double_t l2, Double_t l3, Double_t l4, Double_t l5 ){
        fCutUseVariableCascCosPA = kTRUE; //Automatically switch on!
        fCutVarCascCosPA_Exp0Const = l1;
        fCutVarCascCosPA_Exp0Slope = l2;
        fCutVarCascCosPA_Exp1Const = l3;
        fCutVarCascCosPA_Exp1Slope = l4;
        fCutVarCascCosPA_Const     = l5;
    }
    
    //Variable V0CosPA
    void SetCutUseVarV0CosPA      ( Bool_t lCut )   { fCutUseVariableV0CosPA     = lCut; }
    void SetCutVarV0CosPAExp0Const( Double_t lCut ) { fCutVarV0CosPA_Exp0Const = lCut; }
    void SetCutVarV0CosPAExp0Slope( Double_t lCut ) { fCutVarV0CosPA_Exp0Slope = lCut; }
    void SetCutVarV0CosPAExp1Const( Double_t lCut ) { fCutVarV0CosPA_Exp1Const = lCut; }
    void SetCutVarV0CosPAExp1Slope( Double_t lCut ) { fCutVarV0CosPA_Exp1Slope = lCut; }
    void SetCutVarV0CosPAConst    ( Double_t lCut ) { fCutVarV0CosPA_Const     = lCut; }
    void SetCutVarV0CosPA ( Double_t l1, Double_t l2, Double_t l3, Double_t l4, Double_t l5 ){
        fCutUseVariableV0CosPA = kTRUE; //Automatically switch on!
        fCutVarV0CosPA_Exp0Const = l1;
        fCutVarV0CosPA_Exp0Slope = l2;
        fCutVarV0CosPA_Exp1Const = l3;
        fCutVarV0CosPA_Exp1Slope = l4;
        fCutVarV0CosPA_Const     = l5;
    }
    
    //Variable BBCosPA
    void SetCutUseVarBBCosPA      ( Bool_t lCut )   { fCutUseVariableBBCosPA     = lCut; }
    void SetCutVarBBCosPAExp0Const( Double_t lCut ) { fCutVarBBCosPA_Exp0Const = lCut; }
    void SetCutVarBBCosPAExp0Slope( Double_t lCut ) { fCutVarBBCosPA_Exp0Slope = lCut; }
    void SetCutVarBBCosPAExp1Const( Double_t lCut ) { fCutVarBBCosPA_Exp1Const = lCut; }
    void SetCutVarBBCosPAExp1Slope( Double_t lCut ) { fCutVarBBCosPA_Exp1Slope = lCut; }
    void SetCutVarBBCosPAConst    ( Double_t lCut ) { fCutVarBBCosPA_Const     = lCut; }
    void SetCutVarBBCosPA ( Double_t l1, Double_t l2, Double_t l3, Double_t l4, Double_t l5 ){
        fCutUseVariableBBCosPA = kTRUE; //Automatically switch on!
        fCutVarBBCosPA_Exp0Const = l1;
        fCutVarBBCosPA_Exp0Slope = l2;
        fCutVarBBCosPA_Exp1Const = l3;
        fCutVarBBCosPA_Exp1Slope = l4;
        fCutVarBBCosPA_Const     = l5;
    }
    
    AliCascadeResult::EMassHypo GetMassHypothesis () const { return fMassHypo; }
    Double_t GetMass() const;
    TString GetParticleName() const; 

    //Getters for V0 Cuts
    Double_t GetCutMinRapidity     () const { return fCutMinRapidity; }
    Double_t GetCutMaxRapidity     () const { return fCutMaxRapidity; }
    
    //Getters for V0 Cuts
    Double_t GetCutDCANegToPV     () const { return fCutDCANegToPV; }
    Double_t GetCutDCAPosToPV     () const { return fCutDCAPosToPV; }
    Double_t GetCutDCAV0Daughters () const { return fCutDCAV0Daughters; }
    Double_t GetCutV0CosPA        () const { return fCutV0CosPA; }
    Double_t GetCutV0Radius       () const { return fCutV0Radius; }
    
    //Getters for Cascade Cuts
    Double_t GetCutDCAV0ToPV        () const { return fCutDCAV0ToPV; }
    Double_t GetCutV0Mass           () const { return fCutV0Mass; }
    Double_t GetCutV0MassSigma      () const { return fCutV0MassSigma; }
    Double_t GetCutDCABachToPV      () const { return fCutDCABachToPV; }
    Double_t GetCutDCACascDaughters () const { return fCutDCACascDaughters; }
    Double_t GetCutCascCosPA        () const { return fCutCascCosPA; }
    Double_t GetCutCascRadius       () const { return fCutCascRadius; }
    Double_t GetCutDCABachToBaryon  () const { return fCutDCABachToBaryon; }
    Double_t GetCutBachBaryonCosPA  () const { return fCutBachBaryonCosPA; }
    Double_t GetCutMinV0Lifetime    () const { return fCutMinV0Lifetime; }
    Double_t GetCutMaxV0Lifetime    () const { return fCutMaxV0Lifetime; }
    
    //Miscellaneous
    Double_t GetCutProperLifetime () const { return fCutProperLifetime; }
    Double_t GetCutTPCdEdx () const { return fCutTPCdEdx; }
    Double_t GetCutXiRejection () const { return fCutXiRejection; }
    
    Bool_t GetCutMCPhysicalPrimary    () const { return fCutMCPhysicalPrimary; }
    Bool_t GetCutMCPDGCodeAssociation () const { return fCutMCPDGCodeAssociation; }
    Bool_t GetCutMCUseMCProperties    () const { return fCutMCUseMCProperties; }
    Bool_t GetCutMCSelectBump         () const { return fCutMCSelectBump; }
    
    //Rsn-like bg subtraction (experimental
    Bool_t GetSwapBachelorCharge         () const { return fSwapBachCharge; }
    Bool_t GetSwapBaryon                 () const { return fSwapBaryon;     }
    Bool_t GetSwapV0MesonCharge          () const { return fSwapV0MesonCharge;     }
    Bool_t GetSwapV0BaryonCharge         () const { return fSwapV0BaryonCharge;     }
    
    //Track Quality
    Bool_t GetCutUseITSRefitTracks    () const { return fCutUseITSRefitTracks; }
    Double_t GetCutLeastNumberOfClusters () const { return fCutLeastNumberOfClusters; }
    Double_t GetCutMinEtaTracks () const { return fCutMinEtaTracks; }
    Double_t GetCutMaxEtaTracks () const { return fCutMaxEtaTracks; }
    Double_t GetCutMaxChi2PerCluster () const { return fCutMaxChi2PerCluster; }
    Double_t GetCutMinTrackLength    () const { return fCutMinTrackLength; }
    
    //Variable CascCosPA
    Bool_t GetCutUseVarCascCosPA        () const { return fCutUseVariableCascCosPA;   }
    Double_t GetCutVarCascCosPAExp0Const() const { return fCutVarCascCosPA_Exp0Const; }
    Double_t GetCutVarCascCosPAExp0Slope() const { return fCutVarCascCosPA_Exp0Slope; }
    Double_t GetCutVarCascCosPAExp1Const() const { return fCutVarCascCosPA_Exp1Const; }
    Double_t GetCutVarCascCosPAExp1Slope() const { return fCutVarCascCosPA_Exp1Slope; }
    Double_t GetCutVarCascCosPAConst    () const { return fCutVarCascCosPA_Const;     }

    //Variable V0CosPA
    Bool_t GetCutUseVarV0CosPA        () const { return fCutUseVariableV0CosPA;   }
    Double_t GetCutVarV0CosPAExp0Const() const { return fCutVarV0CosPA_Exp0Const; }
    Double_t GetCutVarV0CosPAExp0Slope() const { return fCutVarV0CosPA_Exp0Slope; }
    Double_t GetCutVarV0CosPAExp1Const() const { return fCutVarV0CosPA_Exp1Const; }
    Double_t GetCutVarV0CosPAExp1Slope() const { return fCutVarV0CosPA_Exp1Slope; }
    Double_t GetCutVarV0CosPAConst    () const { return fCutVarV0CosPA_Const;     }
    
    //Variable BBCosPA
    Bool_t GetCutUseVarBBCosPA        () const { return fCutUseVariableBBCosPA;   }
    Double_t GetCutVarBBCosPAExp0Const() const { return fCutVarBBCosPA_Exp0Const; }
    Double_t GetCutVarBBCosPAExp0Slope() const { return fCutVarBBCosPA_Exp0Slope; }
    Double_t GetCutVarBBCosPAExp1Const() const { return fCutVarBBCosPA_Exp1Const; }
    Double_t GetCutVarBBCosPAExp1Slope() const { return fCutVarBBCosPA_Exp1Slope; }
    Double_t GetCutVarBBCosPAConst    () const { return fCutVarBBCosPA_Const;     }
    
    TH3F* GetHistogram       ()       { return fHisto; }
    TH3F* GetHistogramToCopy () const { return fHisto; }
    
    TProfile *GetProtonProfile       ()       { return fProtonProfile; }
    TProfile *GetProtonProfileToCopy () const { return fProtonProfile; }
    void InitializeProtonProfile(Long_t lNPtBins, Double_t *lPtBins); //Initialize profile, otherwise not stored
    
    //No such thing as feeddown corrections for this result
    //Kept to satisfy AliVWeakResult published interface, don't use in this case
    TH3F* GetHistogramFeeddown       ()       { return 0x0; }
    TH3F* GetHistogramFeeddownToCopy () const { return 0x0; }
    
    Bool_t HasSameCuts( AliVWeakResult *lCompare, Bool_t lCheckdEdx = kTRUE );
    void Print();
    
    
private:

    AliCascadeResult::EMassHypo fMassHypo; //For determining invariant mass

    //Basic acceptance criteria
    Double_t fCutMinRapidity; //min rapidity
    Double_t fCutMaxRapidity; //max rapidity
    
    //V0 Selection Criteria
    Double_t fCutDCANegToPV;    //v0 vertexer 1
    Double_t fCutDCAPosToPV;    //v0 vertexer 2
    Double_t fCutDCAV0Daughters;//v0 vertexer 3
    Double_t fCutV0CosPA;       //v0 vertexer 4
    Double_t fCutV0Radius;      //v0 vertexer 5
    
    //Cascade Selection Criteria
    Double_t fCutDCAV0ToPV;        //ca vertexer 1
    Double_t fCutV0Mass;           //ca vertexer 2
    Double_t fCutV0MassSigma;           //ca vertexer 2bis
    Double_t fCutDCABachToPV;      //ca vertexer 3
    Double_t fCutDCACascDaughters; //ca vertexer 4
    Double_t fCutCascCosPA;        //ca vertexer 5
    Double_t fCutCascRadius;       //ca vertexer 6
    Double_t fCutDCABachToBaryon;  //extra selection on dca bach-pos (experimental)
    Double_t fCutBachBaryonCosPA;  //extra selection on bach-baryon CosPA (experimental)
    Double_t fCutMinV0Lifetime; //min V0 lifetime (cm/c)
    Double_t fCutMaxV0Lifetime; //max V0 lifetime (cm/c)
    
    Double_t fCutProperLifetime;
    Double_t fCutTPCdEdx;
    Double_t fCutXiRejection; //Xi rejection (for omega analysis only!)
    
    Bool_t fCutMCPhysicalPrimary; //IsPhysicalPrimary requirement
    Bool_t fCutMCPDGCodeAssociation; //Associate with correct PDG code
    Bool_t fCutMCUseMCProperties; //Use true MC pT, y
    Bool_t fCutMCSelectBump; //select bachelor and baryon from a single lambda decay
    
    Bool_t fSwapBachCharge; //select bachelor with improper signal for desired particle (bg)
    Bool_t fSwapBaryon; //select lambda/antilambda improperly for desired particle (bg)
    Bool_t fSwapV0MesonCharge; //swap V0 meson daughter charge
    Bool_t fSwapV0BaryonCharge; //swap V0 baryon daughter charge
    
    //Track selections
    Bool_t fCutUseITSRefitTracks; //Use ITS refit tracks (will kill efficiency at high pT!)
    Double_t fCutLeastNumberOfClusters; //min number of TPC clusters
    Double_t fCutMinEtaTracks; //Minimum eta value for daughter tracks (usually -0.8)
    Double_t fCutMaxEtaTracks; //Maximum eta value for daughter tracks (usually +0.8)
    Double_t fCutMaxChi2PerCluster; //Max chi2/clusters
    Double_t fCutMinTrackLength; //Minimum track length in the active TPC zone
    
    //Experimental: pt-variable cascade cosPA
    //Warning: if this cut is tighter than fCutCascCosPA, this gets used instead!
    Bool_t fCutUseVariableCascCosPA;
    Double_t fCutVarCascCosPA_Exp0Const;
    Double_t fCutVarCascCosPA_Exp0Slope;
    Double_t fCutVarCascCosPA_Exp1Const;
    Double_t fCutVarCascCosPA_Exp1Slope;
    Double_t fCutVarCascCosPA_Const;
    
    //Experimental: pt-variable V0 cosPA
    //Warning: if this cut is tighter than fCutV0CosPA, this gets used instead!
    Bool_t fCutUseVariableV0CosPA;
    Double_t fCutVarV0CosPA_Exp0Const;
    Double_t fCutVarV0CosPA_Exp0Slope;
    Double_t fCutVarV0CosPA_Exp1Const;
    Double_t fCutVarV0CosPA_Exp1Slope;
    Double_t fCutVarV0CosPA_Const;
    
    //Experimental: pt-variable BB cosPA
    //Warning: if this cut is tighter than fCutBachBaryonCosPA, this gets used instead!
    Bool_t fCutUseVariableBBCosPA;
    Double_t fCutVarBBCosPA_Exp0Const;
    Double_t fCutVarBBCosPA_Exp0Slope;
    Double_t fCutVarBBCosPA_Exp1Const;
    Double_t fCutVarBBCosPA_Exp1Slope;
    Double_t fCutVarBBCosPA_Const;
    
    TH3F *fHisto; //Histogram for storing output with these configurations
    
    TProfile *fProtonProfile; //Histogram for bookkeeping proton momenta
    
    ClassDef(AliCascadeResult, 26)
    // 1 - original implementation
    // 2 - MC association implementation (disabled in real data analysis)
    // 3 - Variable binning constructor + re-order variables in main output for convenience
    // 4 - Xi rejection added
    // 5 - fixes to constructor, destructor, tuning
    // 6 - addition of UseMCProperties flag
    // 7 - Adjustments, tuning, constructor improvements
    // 8 - Experimental lambda cut added
    // 9 - Cleanup + experimental bach-pos DCA cut added
    // 10 - Variable CascCosPA parametrization added (exp0 + exp1 + pol0)
    // 11 - Variable V0CosPA parametrization added (exp0 + exp1 + pol0)
    // 12 - Added explicit bump selection flag
    // 13 - Added Bach Baryon CosPA (experimental)
    // 14 - Addition of GetMass + inherit from AliVWeakResult
    // 15 - Adjustments for gen-purpose functionality
    // 16 - Addition of ITSrefit track requirement for cross-checks
    // 17 - Addition of eta window selection
    // 18 - Max chi2/clusters, min track length added for cross-checking
    // 19 - Settable invariant mass binning constructor
    // 20 - Configuration flags for rsn-like bg estimation (experimental)
    // 21 - swap v0 meson charge addition
    // 22 - swap v0 baryon charge addition
    // 23 - Parametric Bach Baryon CosPA
    // 24 - addition of NSigma cut for Lambda mass (requires pre-configured task!)
    // 25 - addition of rapidity selection (to enable 2.76 TeV re-analysis corrections)
    // 26 - addition of proton momenta histogram (for G3/F correction, 2.76 TeV re-analysis corrections)
};
#endif
