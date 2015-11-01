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
    Float_t GetMultiplicityPercentile(TString lName);
    Float_t GetZ(TString lName) { return GetEstimator(lName.Data())->GetZ(); }
    
    //Master "Evaluate"
    void Evaluate ( AliMultInput *lInput );
    
    //Get ready: prepare/optimize TFormulas
    void Setup(const AliMultInput *lInput);
    
    TList *GetEstimatorList() { return fEstimatorList; } 
    
private:
    Long_t fNEsts; //Number of estimators
    TList *fEstimatorList; //List containing all AliMultEstimators
    
    ClassDef(AliMultSelection, 1)
};
#endif
