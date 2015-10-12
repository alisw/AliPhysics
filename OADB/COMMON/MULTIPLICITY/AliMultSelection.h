#ifndef AliMultSelection_H
#define AliMultSelection_H

#include <iostream>
#include "TNamed.h"

using namespace std;

class AliMultSelection : public TNamed {
    
public:
    AliMultSelection();
    AliMultSelection(const char * name, const char * title = "Mult Estimator");
    AliMultSelection(AliMultSelection *lCopyMe);
    ~AliMultSelection();
    void Clear(Option_t* = "") {}; //dummy

    void     PrintInfo();
    
    void     AddEstimator ( AliMultEstimator *lEst ) { fEstimatorList->Add(lEst); fNEsts++; }
    AliMultEstimator* GetEstimator (TString lName) { return ((AliMultEstimator*)fEstimatorList->FindObject(lName.Data()) ); }
    AliMultEstimator* GetEstimator (Long_t lEstIdx) { return ((AliMultEstimator*)fEstimatorList->At( lEstIdx ) ); }
    Long_t GetNEstimators () { return fNEsts; }
    
    //User Functions to get percentiles
    Float_t GetMultiplicityPercentile(TString lName);
    Float_t GetZ(TString lName) { return GetEstimator(lName.Data())->GetZ(); }
    
    //Master "Evaluate"
    void Evaluate ( AliMultInput *lInput );
    
    TList *GetEstimatorList() { return fEstimatorList; } 
    
private:
    Long_t fNEsts; //Number of estimators
    TList *fEstimatorList; //List containing all AliMultEstimators
    
    ClassDef(AliMultSelection, 1)
};
#endif
