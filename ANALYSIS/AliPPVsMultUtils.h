#ifndef AliPPVsMultUtils_H
#define AliPPVsMultUtils_H

#include "TObject.h"

class AliVEvent;
class AliVVertex;
class AliESDEvent; 
class AliAODEvent;

class AliPPVsMultUtils : public TObject {
    
public:
    
    AliPPVsMultUtils();
    virtual ~AliPPVsMultUtils(){};
    
    //Extra const
    AliPPVsMultUtils(const AliPPVsMultUtils& pd);
    AliPPVsMultUtils &operator=(const AliPPVsMultUtils &c);

    //Utility functions
    //for the base virtual event class: all methods are common
    Float_t GetMultiplicityPercentile(AliVEvent *event, TString lMethod = "V0M");
    Bool_t LoadCalibration(Int_t lLoadThisCalibration);
    
private:
    
    Int_t fRunNumber; // for control of run changes
    Bool_t fCalibrationLoaded; // control flag
    
    TH1F *fBoundaryHisto_V0M;
    TH1F *fBoundaryHisto_V0A;
    TH1F *fBoundaryHisto_V0C;
    TH1F *fBoundaryHisto_V0MEq;
    TH1F *fBoundaryHisto_V0AEq;
    TH1F *fBoundaryHisto_V0CEq;
    
    ClassDef(AliPPVsMultUtils,1) // base helper class
};
#endif

