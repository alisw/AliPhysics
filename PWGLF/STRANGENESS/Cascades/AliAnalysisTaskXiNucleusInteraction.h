#ifndef AliAnalysisTaskXiNucleusInteraction_cxx
#define AliAnalysisTaskXiNucleusInteraction_cxx

#include "AliAnalysisTaskSE.h"
#include "AliMultSelection.h"
#include "AliAnalysisUtils.h"
#include "AliESDtrackCuts.h"
#include "AliAnalysisTask.h"
#include "AliPIDResponse.h"
#include "AliESDcascade.h"
#include "AliESDVertex.h"
#include "AliEventCuts.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "TVector2.h"
#include "TVector3.h"
#include "AliESDv0.h"
#include "TList.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"

//________________________________________________________________________________________________________________________________________________
class AliAnalysisTaskXiNucleusInteraction : public AliAnalysisTaskSE {
    
public:
    AliAnalysisTaskXiNucleusInteraction();
    AliAnalysisTaskXiNucleusInteraction(const char *name);
    virtual ~AliAnalysisTaskXiNucleusInteraction();
    
    //General Functions
    virtual void UserCreateOutputObjects();
    virtual void UserExec  (Option_t *option);
    virtual void Terminate (Option_t *);
    
    //Set Multiplicity
    void SetMultiplicity (UInt_t trigger, Double_t multLow, Double_t multHigh) {
        
        fTrigger  = trigger;
        fMultLow  = multLow;
        fMultHigh = multHigh;
    }
    
    //Track Selection
    void SetTrackSelectionCriteria();

    //User Functions
    Bool_t GetEvent ();
    Bool_t PassedTrackQualityCuts     (AliESDtrack *track);
    Bool_t PassedCascadeSelectionCuts (AliESDcascade *cascade);
    Bool_t PassedPIDSelection         (AliESDtrack *track, AliPID::EParticleType type);
    Bool_t IsXiCandidate              (AliESDcascade *cascade, AliESDtrack *pos, AliESDtrack *neg, AliESDtrack *bac, Double_t &m, TVector3 &momentum);
    Bool_t IsAntiXiCandidate          (AliESDcascade *cascade, AliESDtrack *pos, AliESDtrack *neg, AliESDtrack *bac, Double_t &m, TVector3 &momentum);
    TLorentzVector Boost              (TLorentzVector R, TVector3 beta_vect);
    
    //Standard Pile-up Rejection & Event Cuts
    AliEventCuts fESDeventCuts;//

private:
    //Event Pointers
    AliESDEvent       *fESDevent;//!
    AliPIDResponse    *fPIDResponse;//!
    AliAnalysisUtils  *fUtils;//!
    TList             *fOutputList;//!
    TList             *fQAList;//!
    
    //Track Quality Cuts
    AliESDtrackCuts *fESDtrackCuts;//!
    
    //Trigger
    UInt_t   fTrigger;//
    Double_t fMultLow;//
    Double_t fMultHigh;//

    //Histograms
    TH1F *hNumberOfEvents;//!
    TH1F *hNumberOfXi;//!
    TH1F *hNumberOfAntiXi;//!
    
    //Invariant-Mass Plots
    TH2F *hMassXi_vs_P[7];//!
    TH2F *hMassAntiXi_vs_P[7];//!
    
    //Pointing Angle
    TH2F *hXiPointingAngle_vs_P[7];//!
    TH2F *hAntiXiPointingAngle_vs_P[7];//!
    
    //Scattering Angle
    TH2F *hXiScatteringAngle_vs_P[7];//!
    TH2F *hAntiXiScatteringAngle_vs_P[7];//!

    
    AliAnalysisTaskXiNucleusInteraction(const AliAnalysisTaskXiNucleusInteraction&);
    AliAnalysisTaskXiNucleusInteraction& operator=(const AliAnalysisTaskXiNucleusInteraction&);
    
    ClassDef(AliAnalysisTaskXiNucleusInteraction, 1);
};
//________________________________________________________________________________________________________________________________________________

#endif
