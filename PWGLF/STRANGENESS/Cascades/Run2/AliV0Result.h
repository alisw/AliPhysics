#ifndef AliV0Result_H
#define AliV0Result_H
#include <TNamed.h>
#include <TList.h>
#include <TH3F.h>
#include <TProfile.h>
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

    //Settable centrality / momentum / invmass binning
    AliV0Result(const char * name, AliV0Result::EMassHypo lMassHypo, const char * title, Long_t lNCentBins, Double_t *lCentBins, Long_t lNPtBins, Double_t *lPtBins, Long_t lNMassBins, Double_t lMinMass, Double_t lMaxMass);
    
    //Specific uses
    AliV0Result(AliV0Result *lCopyMe, TString lNewName );
    AliV0Result(const AliV0Result& lCopyMe, TString lNewName );
    ~AliV0Result();
    
    void Clear(Option_t* = "") {}; //dummy
    
    AliV0Result& operator=(const AliV0Result& lCopyMe);
    
    Long64_t Merge(TCollection *hlist);
    
    //Acceptance
    void SetCutMinRapidity      ( Double_t lCut ) { fCutMinRapidity       = lCut; }
    void SetCutMaxRapidity      ( Double_t lCut ) { fCutMaxRapidity       = lCut; }
    
    void SetCutV0Radius       ( Double_t lCut ) { fCutV0Radius         = lCut; }
    void SetCutDCANegToPV     ( Double_t lCut ) { fCutDCANegToPV       = lCut; }
    void SetCutDCAPosToPV     ( Double_t lCut ) { fCutDCAPosToPV       = lCut; }
    void SetCutDCAV0Daughters ( Double_t lCut ) { fCutDCAV0Daughters   = lCut; }
    void SetCutV0CosPA        ( Double_t lCut ) { fCutV0CosPA          = lCut; }
    void SetCutProperLifetime    ( Double_t lCut ) { fCutProperLifetime   = lCut; }
    
    void SetCutCompetingV0Rejection ( Double_t lCut ) { fCutCompetingV0Rejection   = lCut; }
    void SetCutArmenteros           ( Bool_t lCut   ) { fCutArmenteros        = lCut; }
    void SetCutArmenterosParameter  ( Double_t lCut ) { fCutArmenterosParameter        = lCut; }
    void SetCutTPCdEdx              ( Double_t lCut ) { fCutTPCdEdx           = lCut; }
    void SetCutMinBaryonMomentum    ( Double_t lCut ) { fCutMinBaryonMomentum = lCut; }
    
    //MC specific
    void SetCutMCPhysicalPrimary    ( Bool_t lCut ) { fCutMCPhysicalPrimary    = lCut; }
    void SetCutMCLambdaFromPrimaryXi( Bool_t lCut ) { fCutMCLambdaFromPrimaryXi= lCut; }
    void SetCutMCPDGCodeAssociation ( Bool_t lCut ) { fCutMCPDGCodeAssociation = lCut; }
    void SetCutMCUseMCProperties    ( Bool_t lCut ) { fCutMCUseMCProperties    = lCut; }
    
    //Track Quality
    void SetCutUseITSRefitTracks    ( Bool_t lCut ) { fCutUseITSRefitTracks    = lCut; }
    void SetCutLeastNumberOfCrossedRows             ( Double_t lCut ) { fCutLeastNumberOfCrossedRows = lCut; }
    void SetCutLeastNumberOfCrossedRowsOverFindable ( Double_t lCut ) { fCutLeastNumberOfCrossedRowsOverFindable = lCut; }
    void SetCutMinEtaTracks      ( Double_t lCut ) { fCutMinEtaTracks      = lCut; }
    void SetCutMaxEtaTracks      ( Double_t lCut ) { fCutMaxEtaTracks      = lCut; }
    void SetCutMaxChi2PerCluster ( Double_t lCut ) { fCutMaxChi2PerCluster = lCut; }
    void SetCutMinTrackLength    ( Double_t lCut ) { fCutMinTrackLength    = lCut; }
    
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
    
    //Use OTF V0s
    void SetUseOnTheFly ( Bool_t lCut ) { fUseOnTheFly = lCut; } 
    
    //Feeddown matrix initializer
    void InitializeFeeddownMatrix(Long_t lNLambdaPtBins, Double_t *lLambdaPtBins,
                                  Long_t lNXiPtPins, Double_t *lXiPtPins,
                                  Long_t lNCentBins, Double_t *lCentBins );
    
    AliV0Result::EMassHypo GetMassHypothesis () const { return fMassHypo; }
    Double_t GetMass() const;
    TString GetParticleName() const;
    
    //Getters for V0 Cuts
    Double_t GetCutMinRapidity     () const { return fCutMinRapidity; }
    Double_t GetCutMaxRapidity     () const { return fCutMaxRapidity; }
    
    Double_t GetCutV0Radius       () const { return fCutV0Radius; }
    Double_t GetCutDCANegToPV     () const { return fCutDCANegToPV; }
    Double_t GetCutDCAPosToPV     () const { return fCutDCAPosToPV; }
    Double_t GetCutDCAV0Daughters () const { return fCutDCAV0Daughters; }
    Double_t GetCutV0CosPA        () const { return fCutV0CosPA; }
    Double_t GetCutProperLifetime () const { return fCutProperLifetime; }

    Double_t GetCutCompetingV0Rejection () const { return fCutCompetingV0Rejection; }
    Bool_t   GetCutArmenteros           () const { return fCutArmenteros; }
    Double_t GetCutArmenterosParameter  () const { return fCutArmenterosParameter; }
    Double_t GetCutTPCdEdx              () const { return fCutTPCdEdx; }
    Double_t GetCutMinBaryonMomentum    () const { return fCutMinBaryonMomentum; }
    
    Bool_t GetCutMCPhysicalPrimary    () const { return fCutMCPhysicalPrimary; }
    Bool_t GetCutMCLambdaFromPrimaryXi() const { return fCutMCLambdaFromPrimaryXi; }
    Bool_t GetCutMCPDGCodeAssociation () const { return fCutMCPDGCodeAssociation; }
    Bool_t GetCutMCUseMCProperties    () const { return fCutMCUseMCProperties; }
    
    //Track Quality
    Bool_t GetCutUseITSRefitTracks    () const { return fCutUseITSRefitTracks; }
    Double_t GetCutLeastNumberOfCrossedRows             () const { return fCutLeastNumberOfCrossedRows; }
    Double_t GetCutLeastNumberOfCrossedRowsOverFindable () const { return fCutLeastNumberOfCrossedRowsOverFindable; }
    Double_t GetCutMinEtaTracks      () const { return fCutMinEtaTracks; }
    Double_t GetCutMaxEtaTracks      () const { return fCutMaxEtaTracks; }
    Double_t GetCutMaxChi2PerCluster () const { return fCutMaxChi2PerCluster; }
    Double_t GetCutMinTrackLength    () const { return fCutMinTrackLength; }
    
    //Variable V0CosPA
    Bool_t GetCutUseVarV0CosPA        () const { return fCutUseVariableV0CosPA;   }
    Double_t GetCutVarV0CosPAExp0Const() const { return fCutVarV0CosPA_Exp0Const; }
    Double_t GetCutVarV0CosPAExp0Slope() const { return fCutVarV0CosPA_Exp0Slope; }
    Double_t GetCutVarV0CosPAExp1Const() const { return fCutVarV0CosPA_Exp1Const; }
    Double_t GetCutVarV0CosPAExp1Slope() const { return fCutVarV0CosPA_Exp1Slope; }
    Double_t GetCutVarV0CosPAConst    () const { return fCutVarV0CosPA_Const;     }
    
    //Use OTF V0s
    Bool_t GetUseOnTheFly() const { return fUseOnTheFly; }
    
    TH3F* GetHistogram       ()       { return fHisto; }
    TH3F* GetHistogramToCopy () const { return fHisto; }
    
    //Proton Profile - not implemented for V0s so far
    TProfile* GetProtonProfile       ()       { return 0x0; }
    TProfile* GetProtonProfileToCopy () const { return 0x0; }

    TH3F* GetHistogramFeeddown       ()       { return fHistoFeeddown; }
    TH3F* GetHistogramFeeddownToCopy () const { return fHistoFeeddown; }
    
    Bool_t HasSameCuts( AliVWeakResult *lCompare, Bool_t lCheckdEdx = kTRUE );
    void Print();
    
private:
    //V0 Selection Criteria
    AliV0Result::EMassHypo fMassHypo; //For determining invariant mass

    //Basic acceptance criteria
    Double_t fCutMinRapidity; //min rapidity
    Double_t fCutMaxRapidity; //max rapidity
    
    Double_t fCutV0Radius;
    Double_t fCutDCANegToPV;
    Double_t fCutDCAPosToPV;
    Double_t fCutDCAV0Daughters;
    Double_t fCutV0CosPA;
    Double_t fCutProperLifetime;
    Double_t fCutCompetingV0Rejection;
    Bool_t fCutArmenteros;
    Double_t fCutArmenterosParameter;
    Double_t fCutTPCdEdx;
    Double_t fCutMinBaryonMomentum;
    
    Bool_t fCutMCPhysicalPrimary; //IsPhysicalPrimary requirement
    Bool_t fCutMCLambdaFromPrimaryXi; //Checking for feeddown contributions
    Bool_t fCutMCPDGCodeAssociation; //Associate with correct PDG code
    Bool_t fCutMCUseMCProperties; //Use true MC pT, y

    //Track selections
    Double_t fCutLeastNumberOfCrossedRows;
    Double_t fCutLeastNumberOfCrossedRowsOverFindable;
    Bool_t fCutUseITSRefitTracks; //Use ITS refit tracks (will kill efficiency at high pT!)
    Double_t fCutMinEtaTracks; //Minimum eta value for daughter tracks (usually -0.8)
    Double_t fCutMaxEtaTracks; //Maximum eta value for daughter tracks (usually +0.8)
    Double_t fCutMaxChi2PerCluster; //Max chi2/clusters
    Double_t fCutMinTrackLength; //Minimum track length in the active TPC zone
    
    //Experimental: pt-variable V0 cosPA
    //Warning: if this cut is tighter than fCutV0CosPA, this gets used instead!
    Bool_t fCutUseVariableV0CosPA;
    Double_t fCutVarV0CosPA_Exp0Const;
    Double_t fCutVarV0CosPA_Exp0Slope;
    Double_t fCutVarV0CosPA_Exp1Const;
    Double_t fCutVarV0CosPA_Exp1Slope;
    Double_t fCutVarV0CosPA_Const;
    
    //Master switch to use on-the-fly candidates
    Bool_t fUseOnTheFly; //if zero -> offline, if kTRUE -> go on-the-fly
    
    TH3F *fHisto; //Histogram for storing output with these configurations
    TH3F *fHistoFeeddown; //Feeddown matrix (optional)
    
    ClassDef(AliV0Result, 16)
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
    //12 - Addition of eta window selection
    //13 - Max chi2/clusters, min track length for checking
    //14 - added possibility to select on-the-fly V0 candidates
    //15 - added proton profile (dummy as of now)
    //16 - added configurable AP cut
};
#endif
