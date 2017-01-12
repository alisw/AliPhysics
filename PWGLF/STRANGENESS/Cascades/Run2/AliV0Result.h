#ifndef AliV0Result_H
#define AliV0Result_H
#include <TNamed.h>
#include <TList.h>
#include <TH3F.h>
#include "AliVWeakResult.h"

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// TObject to hold V0 configuration + results histogram 
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class AliV0Result : public AliVWeakResult {
    
public:
    enum EMassHypo {
        kK0Short   = 0,
        kLambda    = 1,
        kAntiLambda= 2
    };
    
    //Dummy Constructor
    AliV0Result();
    
    //Standard Constructor
    AliV0Result(const char * name, AliV0Result::EMassHypo lMassHypo = AliV0Result::kK0Short, const char * title = "V0 Result");
    
    //Variable-Binning Constructor:
    // Binning in ( centrality , momentum ) can be chosen and invariant mass is fixed at defaults
    AliV0Result(const char * name, AliV0Result::EMassHypo lMassHypo, const char * title, Long_t lNCentBins, Double_t *lCentBins, Long_t lNPtBins, Double_t *lPtBins);
    
    //Specific uses
    AliV0Result(AliV0Result *lCopyMe, TString lNewName );
    AliV0Result(const AliV0Result& lCopyMe, TString lNewName );
    ~AliV0Result();
    
    void Clear(Option_t* = "") {}; //dummy
    
    AliV0Result& operator=(const AliV0Result& lCopyMe);
    
    Long64_t Merge(TCollection *hlist);
    
    void SetCutV0Radius       ( Double_t lCut ) { fCutV0Radius         = lCut; }
    void SetCutDCANegToPV     ( Double_t lCut ) { fCutDCANegToPV       = lCut; }
    void SetCutDCAPosToPV     ( Double_t lCut ) { fCutDCAPosToPV       = lCut; }
    void SetCutDCAV0Daughters ( Double_t lCut ) { fCutDCAV0Daughters   = lCut; }
    void SetCutV0CosPA        ( Double_t lCut ) { fCutV0CosPA          = lCut; }
    void SetCutProperLifetime    ( Double_t lCut ) { fCutProperLifetime   = lCut; }
    void SetCutLeastNumberOfCrossedRows             ( Double_t lCut ) { fCutLeastNumberOfCrossedRows = lCut; }
    void SetCutLeastNumberOfCrossedRowsOverFindable ( Double_t lCut ) { fCutLeastNumberOfCrossedRowsOverFindable = lCut; }
    void SetCutCompetingV0Rejection ( Double_t lCut ) { fCutCompetingV0Rejection   = lCut; }
    void SetCutArmenteros           ( Bool_t lCut   ) { fCutArmenteros        = lCut; }
    void SetCutTPCdEdx              ( Double_t lCut ) { fCutTPCdEdx           = lCut; }
    void SetCutMinBaryonMomentum    ( Double_t lCut ) { fCutMinBaryonMomentum = lCut; }
    
    //MC specific
    void SetCutMCPhysicalPrimary    ( Bool_t lCut ) { fCutMCPhysicalPrimary    = lCut; }
    void SetCutMCLambdaFromPrimaryXi( Bool_t lCut ) { fCutMCLambdaFromPrimaryXi= lCut; }
    void SetCutMCPDGCodeAssociation ( Bool_t lCut ) { fCutMCPDGCodeAssociation = lCut; }
    void SetCutMCUseMCProperties    ( Bool_t lCut ) { fCutMCUseMCProperties    = lCut; }
    
    //Track Quality
    void SetCutUseITSRefitTracks    ( Bool_t lCut ) { fCutUseITSRefitTracks    = lCut; }
    
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
    
    //Feeddown matrix initializer
    void InitializeFeeddownMatrix(Long_t lNLambdaPtBins, Double_t *lLambdaPtBins,
                                  Long_t lNXiPtPins, Double_t *lXiPtPins,
                                  Long_t lNCentBins, Double_t *lCentBins );
    
    AliV0Result::EMassHypo GetMassHypothesis () const { return fMassHypo; }
    Double_t GetMass() const;
    TString GetParticleName() const; 
    
    Double_t GetCutV0Radius       () const { return fCutV0Radius; }
    Double_t GetCutDCANegToPV     () const { return fCutDCANegToPV; }
    Double_t GetCutDCAPosToPV     () const { return fCutDCAPosToPV; }
    Double_t GetCutDCAV0Daughters () const { return fCutDCAV0Daughters; }
    Double_t GetCutV0CosPA        () const { return fCutV0CosPA; }
    Double_t GetCutProperLifetime () const { return fCutProperLifetime; }
    Double_t GetCutLeastNumberOfCrossedRows             () const { return fCutLeastNumberOfCrossedRows; }
    Double_t GetCutLeastNumberOfCrossedRowsOverFindable () const { return fCutLeastNumberOfCrossedRowsOverFindable; }
    Double_t GetCutCompetingV0Rejection () const { return fCutCompetingV0Rejection; }
    Bool_t   GetCutArmenteros           () const { return fCutArmenteros; }
    Double_t GetCutTPCdEdx              () const { return fCutTPCdEdx; }
    Double_t GetCutMinBaryonMomentum    () const { return fCutMinBaryonMomentum; }
    
    Bool_t GetCutMCPhysicalPrimary    () const { return fCutMCPhysicalPrimary; }
    Bool_t GetCutMCLambdaFromPrimaryXi() const { return fCutMCLambdaFromPrimaryXi; }
    Bool_t GetCutMCPDGCodeAssociation () const { return fCutMCPDGCodeAssociation; }
    Bool_t GetCutMCUseMCProperties    () const { return fCutMCUseMCProperties; }
    
    //Track Quality
    Bool_t GetCutUseITSRefitTracks    () const { return fCutUseITSRefitTracks; }
    
    //Variable V0CosPA
    Bool_t GetCutUseVarV0CosPA        () const { return fCutUseVariableV0CosPA;   }
    Double_t GetCutVarV0CosPAExp0Const() const { return fCutVarV0CosPA_Exp0Const; }
    Double_t GetCutVarV0CosPAExp0Slope() const { return fCutVarV0CosPA_Exp0Slope; }
    Double_t GetCutVarV0CosPAExp1Const() const { return fCutVarV0CosPA_Exp1Const; }
    Double_t GetCutVarV0CosPAExp1Slope() const { return fCutVarV0CosPA_Exp1Slope; }
    Double_t GetCutVarV0CosPAConst    () const { return fCutVarV0CosPA_Const;     }
    
    TH3F* GetHistogram       ()       { return fHisto; }
    TH3F* GetHistogramToCopy () const { return fHisto; }

    TH3F* GetHistogramFeeddown       ()       { return fHistoFeeddown; }
    TH3F* GetHistogramFeeddownToCopy () const { return fHistoFeeddown; }
    
    Bool_t HasSameCuts( AliVWeakResult *lCompare, Bool_t lCheckdEdx = kTRUE );
    void Print();
    
private:
    //V0 Selection Criteria
    AliV0Result::EMassHypo fMassHypo; //For determining invariant mass

    Double_t fCutV0Radius;
    Double_t fCutDCANegToPV;
    Double_t fCutDCAPosToPV;
    Double_t fCutDCAV0Daughters;
    Double_t fCutV0CosPA;
    Double_t fCutProperLifetime;
    Double_t fCutLeastNumberOfCrossedRows;
    Double_t fCutLeastNumberOfCrossedRowsOverFindable;
    Double_t fCutCompetingV0Rejection;
    Bool_t fCutArmenteros;
    Double_t fCutTPCdEdx;
    Double_t fCutMinBaryonMomentum;
    
    Bool_t fCutMCPhysicalPrimary; //IsPhysicalPrimary requirement
    Bool_t fCutMCLambdaFromPrimaryXi; //Checking for feeddown contributions
    Bool_t fCutMCPDGCodeAssociation; //Associate with correct PDG code
    Bool_t fCutMCUseMCProperties; //Use true MC pT, y

    Bool_t fCutUseITSRefitTracks; //Use ITS refit tracks (will kill efficiency at high pT!)
    
    //Experimental: pt-variable V0 cosPA
    //Warning: if this cut is tighter than fCutV0CosPA, this gets used instead!
    Bool_t fCutUseVariableV0CosPA;
    Double_t fCutVarV0CosPA_Exp0Const;
    Double_t fCutVarV0CosPA_Exp0Slope;
    Double_t fCutVarV0CosPA_Exp1Const;
    Double_t fCutVarV0CosPA_Exp1Slope;
    Double_t fCutVarV0CosPA_Const;
    
    TH3F *fHisto; //Histogram for storing output with these configurations
    TH3F *fHistoFeeddown; //Feeddown matrix (optional)
    
    ClassDef(AliV0Result, 11)
    // 1 - original implementation
    // 2 - first implementation of MC association (to be adjusted)
    // 3 - Variable binning constructor + re-order variables in main output for convenience
    // 4 - fixes to constructor, destructor, tuning
    // 5 - Use MC true pT, y flag added
    // 6 - Adjustments, tuning, constructor improvements
    // 7 - First implementation of feeddown matrix as standard output
    // 8 - Addition of minimum momentum cut for baryon daughters
    // 9 - Addition of GetMass + inherit from AliVWeakResult
    //10 - Adjustments for gen-purpose functionality
    //11 - Addition of variable CosPA, ITSrefit requirement
};
#endif
