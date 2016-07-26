#ifndef AliMultSelection_H
#define AliMultSelection_H
#include <TNamed.h>
#include <TList.h>
#include "AliMultEstimator.h"

class AliMultInput;

class AliMultSelection : public TNamed {
    
public:
    AliMultSelection();
    AliMultSelection(const char * name, const char * title = "Mult Estimator");
    AliMultSelection(AliMultSelection *lCopyMe);
    AliMultSelection(const AliMultSelection& lCopyMe);
    void Set(AliMultSelection* s);
    ~AliMultSelection();
    
    void Clear(Option_t* = "") {}; //dummy
    void CleanUp();
    AliMultSelection& operator=(const AliMultSelection& lCopyMe);
    void     PrintInfo();

    //Estimator Handling
    void     AddEstimator ( AliMultEstimator *lEst );
    AliMultEstimator* GetEstimator (const TString& lName) const;
    AliMultEstimator* GetEstimator (Long_t lEstIdx) const;
    Long_t GetNEstimators () { return fNEsts; }
    
    //User Functions to get percentiles
    Float_t GetMultiplicityPercentile(TString lName, Bool_t lEmbedEvSel = kFALSE);
    Float_t GetZ(TString lName) { return GetEstimator(lName.Data())->GetZ(); }
    
    //Setter and Getter for Event Selection code
    Int_t GetEvSelCode() const { return fEvSelCode; }
    void SetEvSelCode(Int_t lEvSelCodeProv) { fEvSelCode = lEvSelCodeProv; }
    
    //Getter for IsEventSelected 
    Bool_t IsEventSelected();  
    
    //Master "Evaluate"
    void Evaluate ( AliMultInput *lInput );
    
    //Get ready: prepare/optimize TFormulas
    void Setup(const AliMultInput *lInput);
    
    TList *GetEstimatorList() { return fEstimatorList; } 
    
    //Setters/Getters for Event Selection Codes
    //Could have been done better, but this will keep backwards compatibility
    void         SetThisEventVtxZCut ( Bool_t lBoo ) { fThisEvent_VtxZCut = lBoo; }
    Bool_t GetThisEventVtxZCut () { return fThisEvent_VtxZCut; }
    void         SetThisEventIsNotPileup ( Bool_t lBoo ) { fThisEvent_IsNotPileup = lBoo; }
    Bool_t GetThisEventIsNotPileup () { return fThisEvent_IsNotPileup; }
    void         SetThisEventIsNotPileupMV ( Bool_t lBoo ) { fThisEvent_IsNotPileupMV = lBoo; }
    Bool_t GetThisEventIsNotPileupMV () { return fThisEvent_IsNotPileupMV; }
    void         SetThisEventIsNotPileupInMultBins ( Bool_t lBoo ) { fThisEvent_IsNotPileupInMultBins = lBoo; } 
    Bool_t GetThisEventIsNotPileupInMultBins () { return fThisEvent_IsNotPileupInMultBins; }
    void         SetThisEventTriggered ( Bool_t lBoo ) { fThisEvent_Triggered = lBoo; } 
    Bool_t GetThisEventTriggered () { return fThisEvent_Triggered; }
    void         SetThisEventINELgtZERO ( Bool_t lBoo ) { fThisEvent_INELgtZERO = lBoo; } 
    Bool_t GetThisEventINELgtZERO () { return fThisEvent_INELgtZERO; }
    void         SetThisEventHasNoInconsistentVertices ( Bool_t lBoo ) { fThisEvent_HasNoInconsistentVertices = lBoo; } 
    Bool_t GetThisEventHasNoInconsistentVertices () { return fThisEvent_HasNoInconsistentVertices; }
    void         SetThisEventPassesTrackletVsCluster ( Bool_t lBoo ) { fThisEvent_PassesTrackletVsCluster = lBoo; } 
    Bool_t GetThisEventPassesTrackletVsCluster () { return fThisEvent_PassesTrackletVsCluster; }
    void         SetThisEventIsNotAsymmetricInVZERO ( Bool_t lBoo ) { fThisEvent_IsNotAsymmetricInVZERO = lBoo; } 
    Bool_t GetThisEventIsNotAsymmetricInVZERO () { return fThisEvent_IsNotAsymmetricInVZERO; }
    void         SetThisEventIsNotIncompleteDAQ ( Bool_t lBoo ) { fThisEvent_IsNotIncompleteDAQ = lBoo; } 
    Bool_t GetThisEventIsNotIncompleteDAQ () { return fThisEvent_IsNotIncompleteDAQ; }
    
private:
    Long_t fNEsts;    //Number of estimators
    Int_t fEvSelCode; //Event Selection code
    TList *fEstimatorList; //List containing all AliMultEstimators
    
    //Event Characterization Variables - optional
    Bool_t fThisEvent_VtxZCut;                  //!
    Bool_t fThisEvent_IsNotPileup;              //!
    Bool_t fThisEvent_IsNotPileupMV;            //!
    Bool_t fThisEvent_IsNotPileupInMultBins;    //!
    Bool_t fThisEvent_Triggered;                //!
    Bool_t fThisEvent_INELgtZERO;               //! //done with SPD tracklets
    Bool_t fThisEvent_HasNoInconsistentVertices;//!
    Bool_t fThisEvent_PassesTrackletVsCluster;  //!
    Bool_t fThisEvent_IsNotAsymmetricInVZERO;   //!
    Bool_t fThisEvent_IsNotIncompleteDAQ;   //!
    
    ClassDef(AliMultSelection, 4)
    // 1 - original implementation
    // 2 - added fEvSelCode for EvSel bypass + getter changed
    // 3 - added booleans to classify which event criteria are satisfied
    // 4 - added IsEventSelected 
};
#endif
