//Create by Christine Nattrass, Rebecca Scott, Irakli Martashvili
//University of Tennessee at Knoxville
#ifndef AliAnalysisTaskHadEt_cxx
#define AliAnalysisTaskHadEt_cxx

class AliAnalysisHadEt;
class TTree;
class AliVParticle;
class TH1F;
class TH2F;
class TNtuple;
class TObjArray;
class AliESDEvent;
class AliMCParticle;
class TDatabasePDG;

#include "AliAnalysisTaskSE.h"
#include "TObject.h"
#include "AliAnalysisHadEtReconstructed.h"
#include "AliAnalysisHadEtMonteCarlo.h"

/* class ParticleVars : public TObject        // Inherit from TObject to put in TClonesArray */
/*     { */
/*        public: */
	  
/* 	  ParticleVars() : TObject(){} */
/*        Int_t fPdgCode; // from MC */
/*        Int_t fPid; //from ESDs */
/*        Int_t fMass; */
/*        Int_t fCharge; */
/*        Double_t fEt; */
/*        Double_t fPhi; */
/*        Double_t fEta; */

/*        ClassDef(ParticleVars, 1); */
       
/*     }; */

class AliAnalysisTaskHadEt : public AliAnalysisTaskSE {
public:
    AliAnalysisTaskHadEt(const char *name = "AliAnalysisTaskHadEt");
    virtual ~AliAnalysisTaskHadEt() {}

    //  virtual void   ConnectInputData(Option_t *);
    virtual void   UserCreateOutputObjects();
    virtual void   UserExec(Option_t *option);
    virtual void   Terminate(Option_t *);

    virtual void SetTriggerSelection(Bool_t v) {
        fTriggerSelection = v;
    }

    AliESDtrackCuts* GetTPCITSTrackCuts(){return (AliESDtrackCuts*) fOutputList->FindObject("fEsdTrackCuts");}
    AliESDtrackCuts* GetTPCOnlyTrackCuts(){return (AliESDtrackCuts*) fOutputList->FindObject("fEsdTrackCutsTPCOnly");}
    AliESDtrackCuts* GetITSTrackCuts(){return (AliESDtrackCuts*) fOutputList->FindObject("fEsdTrackCutsITS");}

private:

  //Declare it private to avoid compilation warning
    AliAnalysisTaskHadEt & operator = (const AliAnalysisTaskHadEt & g) ;//cpy assignment
    AliAnalysisTaskHadEt(const AliAnalysisTaskHadEt & g) ; // cpy ctor

    bool CheckGoodVertex(AliVParticle *track);
    bool TrackHits(AliVParticle *track, Double_t magField);


    AliESDEvent *fESD;    //ESD object

    TList *fOutputList;

    AliAnalysisHadEtReconstructed *fRecAnalysis;
    AliAnalysisHadEtMonteCarlo *fMCAnalysis;

    TH2F *fHistEtRecvsEtMC;
    
    Bool_t fTriggerSelection;

    Int_t fCount;

    const int fkPhotonPdg;

    const Float_t fkProtonMass;

    TDatabasePDG *fPdgDB;

    class EventVars
    {
       public:
        Double_t fTotEt;
	Double_t fTotEtAcc;
        Double_t fTotEnergy;

        Double_t fTotNeutralEt;
	Double_t fTotNeutralEtAcc;

        Double_t fTotChargedEt;
        Double_t fTotChargedEtAcc;

        Int_t fChargedMultiplicity;
        Int_t fNeutralMultiplicity;

    };
    
    EventVars *fRecEventVars;
    EventVars *fSimEventVars;
    AliESDtrackCuts* esdtrackCutsITSTPC;
    AliESDtrackCuts* esdtrackCutsTPC;
    AliESDtrackCuts* esdtrackCutsITS;
    
    ClassDef(AliAnalysisTaskHadEt, 1); // example of analysis
};

#endif
