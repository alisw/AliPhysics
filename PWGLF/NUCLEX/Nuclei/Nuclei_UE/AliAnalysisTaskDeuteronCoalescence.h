#ifndef AliAnalysisTaskDeuteronCoalescence_cxx
#define AliAnalysisTaskDeuteronCoalescence_cxx

#include "AliMCEventHandler.h"
#include "AliAnalysisTaskSE.h"
#include "AliMultSelection.h"
#include "AliAnalysisUtils.h"
#include "AliAnalysisTask.h"
#include "TLorentzVector.h"
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
#include "TH1F.h"
#include "TH2D.h"
#include "TF1.h"

//___________________________________________________________________________________________________________________________________________________
class AliAnalysisTaskDeuteronCoalescence : public AliAnalysisTaskSE {
        
public:
    AliAnalysisTaskDeuteronCoalescence();
    AliAnalysisTaskDeuteronCoalescence(const char *name);
    virtual ~AliAnalysisTaskDeuteronCoalescence();
        
    //General Functions
    virtual void UserCreateOutputObjects();
    virtual void UserExec  (Option_t *option);
    virtual void Terminate (Option_t *);
    
    //Set Average Charged-Particle Multiplicity, Weights & Deuteron Wave Function
    void SetAverageTransverseMultiplicity (Double_t Nch_Transv) { fAverage_Nch_Transv = Nch_Transv; }
    void SetReshapingProtons (TH1D *hWeight)   { hProtonWeights = hWeight; }
    void SetDeuteronWaveFunc (TF1 *func)       { fDeuteronWF = func; }
    void SetSourceSizeRadius (TH1F *hSourceR0) { hSourceSize = hSourceR0; }

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
    Double_t  GetSpatialDistance           (TLorentzVector p_proton, TLorentzVector p_neutron, TVector3 beta_vect);
    Bool_t    DoCoalescence                (Double_t deltaX, Double_t deltaP, Double_t sigma_p, const char *func);
    TLorentzVector LorentzTransform        (TLorentzVector R, TVector3 beta_vect);

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
    
    //Deuteron Wave Function
    TF1 *fDeuteronWF;//
    
    //Source Radius
    TH1F *hSourceSize;//
    
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
    
    //p_{T} Spectra: Deuterons (Simple Coalescence)
    TH1D *hDeuteronsINELgtZERO_simpleCoal[50];//!
    TH1D *hDeuterons_Toward_simpleCoal[50];//!
    TH1D *hDeuterons_Transv_simpleCoal[50];//!
    
    //p_{T} Spectra: Deuterons (Wigner Gaussian)
    TH1D *hDeuteronsINELgtZERO_wignerGaus[50];//!
    TH1D *hDeuterons_Toward_wignerGaus[50];//!
    TH1D *hDeuterons_Transv_wignerGaus[50];//!
    
    //p_{T} Spectra: Deuterons (Wigner Argonne)
    TH1D *hDeuteronsINELgtZERO_wignerArg[50];//!
    TH1D *hDeuterons_Toward_wignerArg[50];//!
    TH1D *hDeuterons_Transv_wignerArg[50];//!

    
    //QA Histograms & Debug
    TH1D *hRparticles;//!
    
    //Rapidity Distributions
    TH1D *hRapidityProtons;//!
    TH1D *hRapidityNeutrons;//!
    
    //DeltaP Distribution
    TH1D *hDeltaP;//!
    
    //Source Radii
    TH1D *hSourceRadius_Prot;//!
    TH1D *hSourceRadius_Neut;//!
    
    //Control Histograms
    TH1D *hDistanceLab;//!
    TH1D *hDistanceDeut;//!
    TH1D *hDistanceDiff;//!


    
    AliAnalysisTaskDeuteronCoalescence(const AliAnalysisTaskDeuteronCoalescence&);
    AliAnalysisTaskDeuteronCoalescence& operator=(const AliAnalysisTaskDeuteronCoalescence&);
        
    ClassDef(AliAnalysisTaskDeuteronCoalescence, 1);
    
};
//___________________________________________________________________________________________________________________________________________________

#endif
