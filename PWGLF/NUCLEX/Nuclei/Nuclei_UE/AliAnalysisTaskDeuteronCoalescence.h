#ifndef AliAnalysisTaskDeuteronCoalescence_cxx
#define AliAnalysisTaskDeuteronCoalescence_cxx

#include "AliMCEventHandler.h"
#include "AliAnalysisTaskSE.h"
#include "AliMultSelection.h"
#include "AliAnalysisUtils.h"
#include "AliAnalysisTask.h"
#include "AliMCParticle.h"
#include "AliEventCuts.h"
#include "AliESDEvent.h"
#include "AliMCEvent.h"
#include "TVector2.h"
#include "TVector3.h"
#include "AliStack.h"
#include "TList.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"

//_____________________________________________________________________________________________________________________________
class AliAnalysisTaskDeuteronCoalescence : public AliAnalysisTaskSE {
        
public:
    AliAnalysisTaskDeuteronCoalescence();
    AliAnalysisTaskDeuteronCoalescence(const char *name);
    virtual ~AliAnalysisTaskDeuteronCoalescence();
        
    //General Functions
    virtual void UserCreateOutputObjects();
    virtual void UserExec  (Option_t *option);
    virtual void Terminate (Option_t *);
    
    //Set Average Charged-Particle Multiplicity & Weights
    void SetAverageTransverseMultiplicity (Double_t Nch_Transv) { fAverage_Nch_Transv = Nch_Transv; }
    void SetReshapingProtons (TH1D *hWeight) { hProtonWeights = hWeight; }

    //User Functions
    Bool_t    GetEvent ();
    Bool_t    IsINELgtZERO ();
    Int_t     GetLeadingParticle ();
    Int_t     GetTransverseMultiplicity    (Int_t leading_particle_ID);
    Bool_t    IsParticleInTowardRegion     (Double_t phi, Double_t phi_leading);
    Bool_t    IsParticleInTransverseRegion (Double_t phi, Double_t phi_leading);
    Bool_t    IsParticleInAwayRegion       (Double_t phi, Double_t phi_leading);
    Double_t  GetProtonWeight              (Double_t pt);
    Double_t  GetDeuteronWeight            (Double_t pt_prot, Double_t pt_neut);
    Bool_t    IsInjectedParticle           (AliMCParticle *particle);
    Double_t  GetRapidity                  (TVector3 momentum);

    AliEventCuts  fESDeventCuts;//

private:
    AliESDEvent       *fESDevent;//!
    AliMCEvent        *fMCEvent;//!
    AliStack          *fMCstack;//!
    AliMCEventHandler *fMCEventHandler;//!
    TList             *fOutputList;//!
    TList             *fQAList;//!
    
    //Average Multiplicity
    Double_t fAverage_Nch_Transv;//
    
    //Re-shaping Protons 
    TH1D *hProtonWeights;//
    
    //Event Counter
    TH1D *hNumberOfEvents;//!
   
    //General Histograms
    TH1D *hTransverseMult;//!
    TH1D *hRtDistribution;//!

    //p_{T} Spectra: Protons
    TH1D *hProtonsINELgtZERO;//!
    TH1D *hProtonsINELgtZERO_reshaped;//!
    TH1D *hProtons_Toward;//!
    TH1D *hProtons_Transv;//!
    TH1D *hProtons_Away;//!
    
    //p_{T} Spectra: Neutrons
    TH1D *hNeutronsINELgtZERO;//!
    TH1D *hNeutronsINELgtZERO_reshaped;//!
    TH1D *hNeutrons_Toward;//!
    TH1D *hNeutrons_Transv;//!
    TH1D *hNeutrons_Away;//!
    
    //p_{T} Spectra: Deuterons
    TH1D *hDeuteronsINELgtZERO[10];//!
    TH1D *hDeuterons_Toward[10];//!
    TH1D *hDeuterons_Transv[10];//!
    TH1D *hDeuterons_Away[10];//!
    
    //QA Histograms & Debug
    TH1D *hRparticles;//!
    
    //Rapidity Distributions
    TH1D *hRapidityProtons;//!
    TH1D *hRapidityNeutrons;//!

    
    AliAnalysisTaskDeuteronCoalescence(const AliAnalysisTaskDeuteronCoalescence&);
    AliAnalysisTaskDeuteronCoalescence& operator=(const AliAnalysisTaskDeuteronCoalescence&);
        
    ClassDef(AliAnalysisTaskDeuteronCoalescence, 1);
    
};
//_____________________________________________________________________________________________________________________________

#endif
