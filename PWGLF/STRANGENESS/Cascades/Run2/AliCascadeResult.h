#ifndef AliCascadeResult_H
#define AliCascadeResult_H
#include <TNamed.h>
#include <TList.h>
#include <TH3F.h>

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// TObject to hold Cascade configuration + results histogram
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class AliCascadeResult : public TNamed {
    
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
    
    //Specific uses
    AliCascadeResult(AliCascadeResult *lCopyMe, TString lNewName);
    AliCascadeResult(const AliCascadeResult& lCopyMe, TString lNewName);
    ~AliCascadeResult();
    
    void Clear(Option_t* = "") {}; //dummy
    
    AliCascadeResult& operator=(const AliCascadeResult& lCopyMe);
    
    Long64_t Merge(TCollection *hlist);
    
    //Setters for V0 Cuts
    void SetCutDCANegToPV     ( Double_t lCut ) { fCutDCANegToPV       = lCut; }
    void SetCutDCAPosToPV     ( Double_t lCut ) { fCutDCAPosToPV       = lCut; }
    void SetCutDCAV0Daughters ( Double_t lCut ) { fCutDCAV0Daughters   = lCut; }
    void SetCutV0CosPA        ( Double_t lCut ) { fCutV0CosPA          = lCut; }
    void SetCutV0Radius       ( Double_t lCut ) { fCutV0Radius         = lCut; }
    
    //Setters for Cascade Cuts
    void SetCutDCAV0ToPV        ( Double_t lCut ) { fCutDCAV0ToPV         = lCut; }
    void SetCutV0Mass           ( Double_t lCut ) { fCutV0Mass            = lCut; }
    void SetCutDCABachToPV      ( Double_t lCut ) { fCutDCABachToPV       = lCut; }
    void SetCutDCACascDaughters ( Double_t lCut ) { fCutDCACascDaughters  = lCut; }
    void SetCutCascCosPA        ( Double_t lCut ) { fCutCascCosPA         = lCut; }
    void SetCutCascRadius       ( Double_t lCut ) { fCutCascRadius        = lCut; }
    void SetCutDCABachToBaryon  ( Double_t lCut ) { fCutDCABachToBaryon   = lCut; }
    
    //Miscellaneous
    void SetCutProperLifetime        ( Double_t lCut ) { fCutProperLifetime        = lCut; }
    void SetCutLeastNumberOfClusters ( Double_t lCut ) { fCutLeastNumberOfClusters = lCut; }
    void SetCutTPCdEdx               ( Double_t lCut ) { fCutTPCdEdx               = lCut; }
    void SetCutXiRejection           ( Double_t lCut ) { fCutXiRejection           = lCut; }
    
    //MC specific
    void SetCutMCPhysicalPrimary    ( Bool_t lCut ) { fCutMCPhysicalPrimary    = lCut; }
    void SetCutMCPDGCodeAssociation ( Bool_t lCut ) { fCutMCPDGCodeAssociation = lCut; }
    void SetCutMCUseMCProperties    ( Bool_t lCut ) { fCutMCUseMCProperties    = lCut; }
    
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
    
    AliCascadeResult::EMassHypo GetMassHypothesis () const { return fMassHypo; }
    
    //Getters for V0 Cuts
    Double_t GetCutDCANegToPV     () const { return fCutDCANegToPV; }
    Double_t GetCutDCAPosToPV     () const { return fCutDCAPosToPV; }
    Double_t GetCutDCAV0Daughters () const { return fCutDCAV0Daughters; }
    Double_t GetCutV0CosPA        () const { return fCutV0CosPA; }
    Double_t GetCutV0Radius       () const { return fCutV0Radius; }
    
    //Getters for Cascade Cuts
    Double_t GetCutDCAV0ToPV        () const { return fCutDCAV0ToPV; }
    Double_t GetCutV0Mass           () const { return fCutV0Mass; }
    Double_t GetCutDCABachToPV      () const { return fCutDCABachToPV; }
    Double_t GetCutDCACascDaughters () const { return fCutDCACascDaughters; }
    Double_t GetCutCascCosPA        () const { return fCutCascCosPA; }
    Double_t GetCutCascRadius       () const { return fCutCascRadius; }
    Double_t GetCutDCABachToBaryon  () const { return fCutDCABachToBaryon; }
    
    //Miscellaneous
    Double_t GetCutProperLifetime () const { return fCutProperLifetime; }
    Double_t GetCutLeastNumberOfClusters () const { return fCutLeastNumberOfClusters; }
    Double_t GetCutTPCdEdx () const { return fCutTPCdEdx; }
    Double_t GetCutXiRejection () const { return fCutXiRejection; }
    
    Bool_t GetCutMCPhysicalPrimary    () const { return fCutMCPhysicalPrimary; }
    Bool_t GetCutMCPDGCodeAssociation () const { return fCutMCPDGCodeAssociation; }
    Bool_t GetCutMCUseMCProperties    () const { return fCutMCUseMCProperties; }
    
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
    
    TH3F* GetHistogram       ()       { return fHisto; }
    TH3F* GetHistogramToCopy () const { return fHisto; }
    
    Bool_t HasSameCuts( AliCascadeResult *lCompare );
    void Print();
    
    
private:

    AliCascadeResult::EMassHypo fMassHypo; //For determining invariant mass

    //V0 Selection Criteria
    Double_t fCutDCANegToPV;    //v0 vertexer 1
    Double_t fCutDCAPosToPV;    //v0 vertexer 2
    Double_t fCutDCAV0Daughters;//v0 vertexer 3
    Double_t fCutV0CosPA;       //v0 vertexer 4
    Double_t fCutV0Radius;      //v0 vertexer 5
    
    //Cascade Selection Criteria
    Double_t fCutDCAV0ToPV;        //ca vertexer 1
    Double_t fCutV0Mass;           //ca vertexer 2
    Double_t fCutDCABachToPV;         //ca vertexer 3
    Double_t fCutDCACascDaughters; //ca vertexer 4
    Double_t fCutCascCosPA;        //ca vertexer 5
    Double_t fCutCascRadius;       //ca vertexer 6
    Double_t fCutDCABachToBaryon;  //extra selection on dca bach-pos (experimental)
    
    Double_t fCutProperLifetime;
    Double_t fCutLeastNumberOfClusters;
    Double_t fCutTPCdEdx;
    Double_t fCutXiRejection; //Xi rejection (for omega analysis only!)
    
    Bool_t fCutMCPhysicalPrimary; //IsPhysicalPrimary requirement
    Bool_t fCutMCPDGCodeAssociation; //Associate with correct PDG code
    Bool_t fCutMCUseMCProperties; //Use true MC pT, y
    
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
    
    TH3F *fHisto; //Histogram for storing output with these configurations
    
    ClassDef(AliCascadeResult, 11)
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
};
#endif
