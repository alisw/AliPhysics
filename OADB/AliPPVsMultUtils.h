#ifndef AliPPVsMultUtils_H
#define AliPPVsMultUtils_H

#include "TObject.h"
#include "AliVEvent.h"

class AliVEvent;
class AliVVertex;
class AliESDEvent;
class AliAODEvent;

class AliPPVsMultUtils : public TObject {

public:

    AliPPVsMultUtils();
    virtual ~AliPPVsMultUtils() {};

    //Extra const
    AliPPVsMultUtils(const AliPPVsMultUtils& pd);
    AliPPVsMultUtils &operator=(const AliPPVsMultUtils &c);

    //Utility functions
    //for the base virtual event class: all methods are common
    Float_t GetMultiplicityPercentile(AliVEvent *event, TString lMethod = "V0M", Bool_t lEmbedEventSelection = kTRUE);

    //Called internally (automatically)
    Bool_t LoadCalibration(Int_t lLoadThisCalibration);

    //static EvSel Snippets
    static Bool_t IsMinimumBias(AliVEvent* event);
    static Bool_t IsSelectedTrigger                        (AliVEvent* event, AliVEvent::EOfflineTriggerTypes trigType = AliVEvent::kMB);
    static Bool_t IsSelectedTrigger                        (AliVEvent* event, TString trigName);
    static Bool_t IsINELgtZERO                         (AliVEvent *event);
    static Bool_t IsAcceptedVertexPosition             (AliVEvent *event);
    static Bool_t IsNotPileupSPDInMultBins             (AliVEvent *event);
    static Bool_t HasNoInconsistentSPDandTrackVertices (AliVEvent *event);
    static Bool_t IsEventSelected(AliVEvent *event,AliVEvent::EOfflineTriggerTypes trigType = AliVEvent::kMB);
    static Bool_t IsEventSelected(AliVEvent *event, TString trigName);
    
    //Wrapper with fallback to tracklets
    static Int_t GetStandardReferenceMultiplicity (AliVEvent *event, Bool_t lEmbedEventSelection = kTRUE);
    
    Float_t MinVal( Float_t A, Float_t B ); 

private:

    Int_t fRunNumber; // for control of run changes
    Bool_t fCalibrationLoaded; // control flag

    //To store calibration boundaries
    TH1F *fBoundaryHisto_V0M;
    TH1F *fBoundaryHisto_V0A;
    TH1F *fBoundaryHisto_V0C;
    TH1F *fBoundaryHisto_V0MEq;
    TH1F *fBoundaryHisto_V0AEq;
    TH1F *fBoundaryHisto_V0CEq;
    TH1F *fBoundaryHisto_V0B;
    TH1F *fBoundaryHisto_V0Apartial;
    TH1F *fBoundaryHisto_V0Cpartial;
    TH1F *fBoundaryHisto_V0S;
    TH1F *fBoundaryHisto_V0SB;

    //To Store <V0A>, <V0C>, <V0Apartial> and <V0Cpartial> on a run-per-run basis
    TH1D *fAverageAmplitudes; 
    
    ClassDef(AliPPVsMultUtils,3) // base helper class
};
#endif

