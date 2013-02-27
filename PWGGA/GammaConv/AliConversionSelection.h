#ifndef AliConversionSelection_cxx
#define AliConversionSelection_cxx

#include "AliAODConversionMother.h"
#include "AliConversionPhotonBase.h"
#include "AliAODConversionPhoton.h"
#include "AliConversionAODBGHandlerRP.h"
#include "AliConversionCuts.h"
#include "AliConversionMesonCuts.h"
#include "TRandom3.h"
#include "AliVEvent.h"
#include "AliMCEvent.h"
#include "TClonesArray.h"
#include "AliESDtrackCuts.h"

class AliConversionSelection : public TObject{

public:

    AliConversionSelection(AliConversionCuts *convCut, AliConversionMesonCuts *mesonCut);
    AliConversionSelection(TString convCut, TString mesonCut);
    virtual ~AliConversionSelection();

    // Main Functions
    Bool_t ProcessEvent(TClonesArray *photons,AliVEvent *inputEvent,AliMCEvent *mcEvent);
   
    // public getter functions

    Int_t GetNumberOfPi0s(){return fPi0Candidates->GetEntriesFast();}
    Int_t GetNumberOfBGs(){return fBGPi0s->GetEntriesFast();}
    Int_t GetNumberOfPhotons(){return fGoodGammas->GetEntriesFast();}

    Double_t GetMultiplicity(AliVEvent *inputEvent);

    AliAODConversionMother* GetPi0(Int_t index);
    AliAODConversionMother* GetBG(Int_t index);
    AliAODConversionPhoton* GetPhoton(Int_t index);

    TClonesArray *GetPi0Candidates(){return fPi0Candidates;}
    TClonesArray *GetBGPi0s(){return fBGPi0s;}

    // public Set Functions
    void SetInvMassRange(Double_t low,Double_t up){fInvMassRange[0]=low;fInvMassRange[1]=up;};
    void SetInvMassRange(Double_t range[2]){SetInvMassRange(range[0],range[1]);};

    Double_t* GetInvMassRange(){return fInvMassRange;}

    TObjArray *GetGoodGammas(){return fGoodGammas;}

    Int_t GetNumberOfChargedTracks(AliVEvent *inputEvent);
    Double_t GetSPDMult(AliVEvent *inputEvent);
    Double_t GetVZEROMult(AliVEvent *inputEvent);

    Int_t GetEventNumber(AliVEvent *inputEvent);

    TString GetCutString();

protected:
   
    void InitializeBGHandler();
    void CalculatePi0Candidates();
    void CalculateBackground();

    void RotateParticle(AliAODConversionPhoton *gamma,Int_t nDegreesPMBackground);

    Bool_t MesonInMassWindow(AliAODConversionMother *pi0cand);

    AliVEvent *fInputEvent;
    AliMCEvent *fMCEvent;
    AliConversionCuts *fConversionCut;
    AliConversionMesonCuts *fMesonCut;
    AliESDtrackCuts *fESDTrackCuts;
    TObjArray *fGoodGammas; // Pointer to selected photons
    TClonesArray *fPi0Candidates;
    TClonesArray *fBGPi0s;
    TRandom3 *fRandomizer; // Randomizer for Rotation
    AliConversionAODBGHandlerRP *fBGHandler;
    Double_t *fInvMassRange;
    Double_t *fUnsmearedPx;
    Double_t *fUnsmearedPy;
    Double_t *fUnsmearedPz;
    Double_t *fUnsmearedE;
    Int_t fCurrentEventNumber; // Current Event Number
    Bool_t fIsOwner; // Cuts will be deleted when the destructor is called

    AliConversionSelection(const AliConversionSelection&); // not implemented
    AliConversionSelection& operator=(const AliConversionSelection&); // not implemented
  
    ClassDef(AliConversionSelection, 1); // example of analysis
};

#endif

