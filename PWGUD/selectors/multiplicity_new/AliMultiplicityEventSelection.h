#ifndef ALIMULTEVENTSELECTION_H
#define ALIMULTEVENTSELECTION_H

#include <TNamed.h>
#include <TParticle.h>
#include <TH1.h>
#include <AliESDEvent.h>
#include <AliMCEvent.h>
#include <AliMultiplicity.h>
#include <TList.h>
#include <TFile.h>
#include <TCanvas.h>
#include "AliTriggerAnalysis.h"
#include "AliInputEventHandler.h"
#include "AliAnalysisManager.h"

class AliMultiplicityEventSelection : public TNamed {

public:
    enum mEstimators {
        kITSTPC = 1,
        kITSSA = 2,
        kSPD = 3
    };

    AliMultiplicityEventSelection ( const char* name = "DefaultSelection", const char* title = "Default Selection" );

    virtual Bool_t AcceptESDEvent ( const AliESDEvent*, UInt_t ) {
        fCurrentESDEventAccepted = kTRUE;
        return fCurrentESDEventAccepted;
    };
    virtual Bool_t AcceptMCEvent ( const AliMCEvent*, UInt_t ) {
        fCurrentMCEventAccepted = kTRUE;
        return fCurrentMCEventAccepted;
    };
    virtual Bool_t AcceptESDandMCEvent ( const AliESDEvent*, AliMCEvent*, Bool_t, UInt_t ) {
        fCurrentESDEventAccepted = fCurrentMCEventAccepted = kTRUE;
        return kTRUE;
    };

    virtual Long64_t Merge ( const TCollection* mergelist ) = 0;
    virtual Long64_t MergeStats ( const TCollection* mergelist );

    virtual void CreateHistograms();
    virtual void CreateSelectionHistograms() = 0;
    virtual void CreateStats();

    virtual void Result() = 0;

    virtual void SaveHistograms();
    virtual void SaveSelectionHistograms() = 0;

    virtual TCanvas* GetCanvas();

    void SetNeedsMC ( Bool_t needs );
    Bool_t GetNeedsMC();

    void SetSaveHistograms ( Bool_t save );
    Bool_t GetSaveHistograms();

    void SetBatchMode ( Bool_t batch );
    Bool_t GetBatchMode();

    Bool_t IsCurrentESDEventAccepted();
    Bool_t IsCurrentMCEventAccepted();
    virtual void ResetAccepted();
    void UpdateStats();

    virtual void SetCorrelate ( Bool_t correlate );
    virtual Bool_t GetCorrelate();
    virtual ~AliMultiplicityEventSelection() {};

    Bool_t CollectEvent();
    void SetCollectEvent ( Bool_t collect );

    virtual void SetSampleDef ( UInt_t def ) {
        sampleDef = def;
    };

    virtual Bool_t CollectCondition( ) {
        return kFALSE;
    }

    void SetFollowESDselection ( Bool_t follow = kTRUE ) {
        fFollowESDselection = follow;
    };

    Bool_t GetFollowESDselection ( ) {
        return fFollowESDselection;
    }

    void SetFollowMCselection ( Bool_t follow = kTRUE ) {
        fFollowMCselection = follow;
    };

    Bool_t GetFollowMCselection ( ) {
        return fFollowMCselection;
    }

protected:
    Bool_t fNeedsMC;
    Bool_t fCurrentESDEventAccepted;
    Bool_t fCurrentMCEventAccepted;
    Bool_t fSaveHistograms;
    Bool_t fBatchMode;

    Bool_t fCorrelate;

    Bool_t fCollectEvent;

    TH1D* hStats;						//accepted statistics
    UInt_t sampleDef;
    UInt_t currentMask;

    Bool_t fFollowESDselection;
    Bool_t fFollowMCselection;

    Double_t fZcache[3];	//!
    Double_t fZcacheGen;	//!
private:

    ClassDef ( AliMultiplicityEventSelection, 8 );
};

#endif // ALIMULTEVENTSELECTION_H
