#ifndef AliV0Result_H
#define AliV0Result_H
#include <TNamed.h>
#include <TList.h>
#include <TH3F.h>

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// TObject to hold V0 configuration + results histogram
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class AliV0Result : public TNamed {
    
public:
    enum EMassHypo {
        kK0Short   = 0,
        kLambda    = 1,
        kAntiLambda= 2
    };
    
    AliV0Result();
    AliV0Result(const char * name, AliV0Result::EMassHypo lMassHypo = AliV0Result::kK0Short, const char * title = "V0 Result");
    AliV0Result(AliV0Result *lCopyMe);
    AliV0Result(const AliV0Result& lCopyMe);
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
    void SetCutArmenteros           ( Bool_t lCut   ) { fCutArmenteros = lCut; }
    void SetCutTPCdEdx              ( Double_t lCut ) { fCutTPCdEdx    = lCut; }
    
    AliV0Result::EMassHypo GetMassHypothesis () const { return fMassHypo; }
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
    
    TH3F* GetHistogram () { return fHisto; } 
    
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
    
    TH3F *fHisto; //Histogram for storing output with these configurations
    
    ClassDef(AliV0Result, 1)
    // 1 - original implementation
};
#endif
