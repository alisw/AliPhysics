#ifndef AliAnalysisTaskPythiaCoalescence_cxx
#define AliAnalysisTaskPythiaCoalescence_cxx

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTask.h"
#include "TLorentzVector.h"
#include "AliMCParticle.h"
#include "AliAODEvent.h"
#include "AliMCEvent.h"
#include "TVector3.h"
#include "TList.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TF1.h"

//____________________________________________________________________________________________________________________________________
class AliAnalysisTaskPythiaCoalescence : public AliAnalysisTaskSE {
        
public:
    AliAnalysisTaskPythiaCoalescence();
    AliAnalysisTaskPythiaCoalescence(const char *name);
    virtual ~AliAnalysisTaskPythiaCoalescence();
        
    //General Functions
    virtual void UserCreateOutputObjects();
    virtual void UserExec  (Option_t *option);
    virtual void Terminate (Option_t *);
    
    //Set Average Charged-Particle Multiplicity, Weights & Deuteron Wave Function
    void SetAverageTransverseMultiplicity (Double_t Nch_Transv) { fAverage_Nch_Transv = Nch_Transv; }
    void SetReshapingProtons (TH1D *hWeight, TF1 *fWeight)   { hProtonWeights = hWeight; fProtonWeights = fWeight; }
    void SetReshapingProtonsDiff (TH1D *hToward, TH1D *hTransv, TH1D *hAway) {
        
        hWeightToward = hToward;
        hWeightTransv = hTransv;
        hWeightAway   = hAway;
    }
    
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
    Double_t  GetProtonWeight              (Double_t pt, Int_t region);
    Double_t  GetDeuteronWeight            (Double_t pt_prot, Double_t pt_neut);
    Double_t  GetDeuteronWeight            (Double_t pt_prot, Double_t pt_neut, Int_t prot_Reg, Int_t neut_Reg);
    Bool_t    IsInjectedParticle           (AliMCParticle *particle);
    Double_t  GetSpatialDistance           (TLorentzVector p_proton, TLorentzVector p_neutron, TVector3 beta_vect);
    Bool_t    DoCoalescence                (Double_t deltaX, Double_t deltaP, Double_t sigma_p, const char *func);
    TLorentzVector LorentzTransform        (TLorentzVector R, TVector3 beta_vect);

private:
    AliAODEvent *fAODevent;//!
    AliMCEvent  *fMCEvent;//!
    TList       *fOutputList;//!
    TList       *fQAList;//!
    
    //Average Multiplicity
    Double_t fAverage_Nch_Transv;//
    
    //Re-shaping Protons
    TH1D *hProtonWeights;//
    TF1  *fProtonWeights;//
    TH1D *hWeightToward;//
    TH1D *hWeightTransv;//
    TH1D *hWeightAway;//

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
    TH1D *hProtons_Toward_reshaped;//!
    TH1D *hProtons_Transv_reshaped;//!
    TH1D *hProtons_Away_reshaped;//!

    
    //p_{T} Spectra: Neutrons
    TH1D *hNeutronsINELgtZERO;//!
    TH1D *hNeutronsINELgtZERO_reshaped;//!
    TH1D *hNeutrons_Toward;//!
    TH1D *hNeutrons_Transv;//!
    TH1D *hNeutrons_Away;//!
    TH1D *hNeutrons_Toward_reshaped;//!
    TH1D *hNeutrons_Transv_reshaped;//!
    TH1D *hNeutrons_Away_reshaped;//!

    
    //p_{T} Spectra: Deuterons (Simple Coalescence)
    TH1D *hDeuteronsINELgtZERO_simpleCoal[50];//!
    TH1D *hDeuterons_Toward_simpleCoal[50];//!
    TH1D *hDeuterons_Transv_simpleCoal[50];//!
    
    //p_{T} Spectra: Deuterons (Wigner Gaussian)
    TH1D *hDeuteronsINELgtZERO_wignerGaus[50];//!
    TH1D *hDeuterons_Toward_wignerGaus[50];//!
    TH1D *hDeuterons_Transv_wignerGaus[50];//!

    //p_{T} Spectra: Deuterons (Wigner Gaussian) true
    TH1D *hDeuteronsINELgtZERO_TruewignerGaus;//!
    TH1D *hDeuterons_Toward_TruewignerGaus[50];//!
    TH1D *hDeuterons_Transv_TruewignerGaus[50];//!

    //p_{T} Spectra: Deuterons (Wigner Gaussian)
    TH1D *hDeuteronsINELgtZERO_TruewignerDoubleGaus;//!
    TH1D *hDeuterons_Toward_TruewignerDoubleGaus[50];//!
    TH1D *hDeuterons_Transv_TruewignerDoubleGaus[50];//!

    //p_{T} Spectra: Deuterons (Wigner Argonne)
    TH1D *hDeuteronsINELgtZERO_wignerArg[50];//!
    TH1D *hDeuterons_Toward_wignerArg[50];//!
    TH1D *hDeuterons_Transv_wignerArg[50];//!
    
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
    TH1D *hPtProtonsFirstBinDeut;//!
    
    //Angular Distributions
    TH2D *hDeltaPhi_Toward;//!
    TH2D *hDeltaPhi_Transv;//!
    TH2D *hDeltaPhi_Away;//!
    TH2D *hDeltaPhi_INELgtZERO;//!


    
    AliAnalysisTaskPythiaCoalescence(const AliAnalysisTaskPythiaCoalescence&);
    AliAnalysisTaskPythiaCoalescence& operator=(const AliAnalysisTaskPythiaCoalescence&);
        
    ClassDef(AliAnalysisTaskPythiaCoalescence, 1);
    
};
//____________________________________________________________________________________________________________________________________

#endif
