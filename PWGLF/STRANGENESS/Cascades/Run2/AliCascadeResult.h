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
    AliCascadeResult(AliCascadeResult *lCopyMe);
    AliCascadeResult(const AliCascadeResult& lCopyMe);
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
    
    //Miscellaneous
    void SetCutProperLifetime        ( Double_t lCut ) { fCutProperLifetime        = lCut; }
    void SetCutLeastNumberOfClusters ( Double_t lCut ) { fCutLeastNumberOfClusters = lCut; }
    void SetCutTPCdEdx               ( Double_t lCut ) { fCutTPCdEdx               = lCut; }
    
    //MC specific
    void SetCutMCPhysicalPrimary    ( Bool_t lCut ) { fCutMCPhysicalPrimary    = lCut; }
    void SetCutMCPDGCodeAssociation ( Bool_t lCut ) { fCutMCPDGCodeAssociation = lCut; }

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
    Double_t GetCutDCABachToPV         () const { return fCutDCABachToPV; }
    Double_t GetCutDCACascDaughters () const { return fCutDCACascDaughters; }
    Double_t GetCutCascCosPA        () const { return fCutCascCosPA; }
    Double_t GetCutCascRadius       () const { return fCutCascRadius; }
    
    //Miscellaneous
    Double_t GetCutProperLifetime () const { return fCutProperLifetime; }
    Double_t GetCutLeastNumberOfClusters () const { return fCutLeastNumberOfClusters; }
    Double_t GetCutTPCdEdx () const { return fCutTPCdEdx; }

    Bool_t GetCutMCPhysicalPrimary    () const { return fCutMCPhysicalPrimary; }
    Bool_t GetCutMCPDGCodeAssociation () const { return fCutMCPDGCodeAssociation; }

    TH3F* GetHistogram () { return fHisto; }
    
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
    
    Double_t fCutProperLifetime;
    Double_t fCutLeastNumberOfClusters;
    Double_t fCutTPCdEdx;
    
    Bool_t fCutMCPhysicalPrimary; //IsPhysicalPrimary requirement
    Bool_t fCutMCPDGCodeAssociation; //Associate with correct PDG code
    
    TH3F *fHisto; //Histogram for storing output with these configurations
    
    ClassDef(AliCascadeResult, 3)
    // 1 - original implementation
    // 2 - MC association implementation (disabled in real data analysis)
    // 3 - Variable binning constructor + re-order variables in main output for convenience
};
#endif
