#ifndef AliPPVsMultUtils_H
#define AliPPVsMultUtils_H

#include "TObject.h"

class AliVEvent;
class AliVVertex;

class AliPPVsMultUtils : public TObject {
    
public:
    
    AliPPVsMultUtils();
    virtual ~AliPPVsMultUtils(){};
    
    //Extra const
    AliPPVsMultUtils(const AliPPVsMultUtils& pd);
    AliPPVsMultUtils &operator=(const AliPPVsMultUtils &c);

    //Utility functions
    //for the base virtual event class: all methods are common
    Float_t GetMultiplicityPercentile(AliESDEvent *event, TString lMethod);
    Float_t GetMultiplicityPercentile(AliAODEvent *event, TString lMethod);
    Bool_t LoadCalibration(Int_t lLoadThisCalibration);
    
private:
    
    Int_t fRunNumber; // minimum vertex contributors
    
    TH1F *fBoundaryHisto_V0M;
    TH1F *fBoundaryHisto_V0A;
    TH1F *fBoundaryHisto_V0C;
    TH1F *fBoundaryHisto_V0MEq;
    TH1F *fBoundaryHisto_V0AEq;
    TH1F *fBoundaryHisto_V0CEq;
    
    ClassDef(AliPPVsMultUtils,1) // base helper class
};
#endif

