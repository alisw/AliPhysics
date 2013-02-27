#ifndef __ALICONVERSIONAODBGHANDLERRP_H__
#define __ALICONVERSIONAODBGHANDLERRP_H__

#include "AliLog.h"
#include "TObject.h"
#include "AliAODConversionPhoton.h"
#include "TObjArray.h"
#include "TList.h"
using namespace std;

typedef vector<AliAODConversionPhoton*> AliGammaConversionPhotonVector;   // Vector containing photons
typedef vector<AliGammaConversionPhotonVector> AliGammaConversionBGEventVector;       // Event contains vector of gammas (AliConversionPhotons)
typedef vector<AliGammaConversionBGEventVector> AliGammaConversionMultiplicityVector;  // Multiplicity classes containing event vectors
typedef vector<AliGammaConversionMultiplicityVector> AliGammaConversionBGVector;       // z vertex position ...


class AliConversionAODBGHandlerRP: public TObject{

public:

    AliConversionAODBGHandlerRP(Bool_t IsHeavyIon=kFALSE,Bool_t UseChargedTrackMult=kTRUE,Int_t NEvents=10);
    
    virtual ~AliConversionAODBGHandlerRP();

    Int_t GetZBinIndex(Double_t z) const;
    Int_t GetMultiplicityBinIndex(Int_t mult) const;
    void Initialize();
    Bool_t FindBins(TObjArray * const eventGammas,AliVEvent *fInputEvent,Int_t &zbin,Int_t &mbin);
    Bool_t FindBins(TList * const eventGammas,AliVEvent *fInputEvent,Int_t &zbin,Int_t &mbin);

 
    AliGammaConversionPhotonVector* GetBGGoodGammas(TObjArray * const eventGammas,AliVEvent *fInputEvent,Int_t event);
    AliGammaConversionPhotonVector* GetBGGoodGammas(TList * const eventGammas,AliVEvent *fInputEvent,Int_t event);
    void AddEvent(TObjArray * const eventGammas,AliVEvent *fInputEvent);
    void AddEvent(TList * const eventGammas,AliVEvent *fInputEvent);
    Int_t GetNBGEvents()const {return fNEvents;} // Size of the Pool (20)
    Int_t GetNBGEvents(TObjArray * const eventGammas,AliVEvent *fInputEvent);
    Int_t GetNBGEvents(TList * const eventGammas,AliVEvent *fInputEvent);
    Int_t GetNZBins()const{return fNBinsZ;};
    Int_t GetNMultiplicityBins()const{return fNBinsMultiplicity;};

private:
    Bool_t fIsHeavyIon;
    Bool_t fUseChargedTrackMult;
    Int_t fNEvents;
    Int_t **fBGEventCounter; //! bg counter
    Int_t **fNBGEvents;
    Int_t fNBinsZ; //n z bins
    Int_t fNBinsMultiplicity; //n bins multiplicity
    Double_t *fBinLimitsArrayZ;//! bin limits z array
    Double_t *fBinLimitsArrayMultiplicity;//! bin limit multiplicity array
    AliGammaConversionBGVector fBGPool; //background events

    AliConversionAODBGHandlerRP(AliConversionAODBGHandlerRP &original);
    AliConversionAODBGHandlerRP &operator=(const AliConversionAODBGHandlerRP &ref);

 ClassDef(AliConversionAODBGHandlerRP,0);

};
#endif
