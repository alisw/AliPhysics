#ifndef ALIANALYSISET_H
#define ALIANALYSISET_H
//_________________________________________________________________________
//  Utility Class for transverse energy studies
//  Base class for ESD & MC analysis
//  - reconstruction and MonteCarlo output
//
//*-- Authors: Oystein Djuvsland (Bergen), David Silvermyr (ORNL)
//_________________________________________________________________________

#include "TString.h"

class TTree;
class TH2F;
class TH1F;
class AliVEvent;
class TList;
class AliESDtrackCuts;
class Rtypes;
class TDatabasePDG;
class AliAnalysisEtCuts;

class AliAnalysisEt
{
public:
   
    AliAnalysisEt();
    virtual ~AliAnalysisEt();
  
public:
  
    /** Analyse the event! */
    virtual Int_t AnalyseEvent(AliVEvent *event);

    /** Fill the objects you want to output, classes which add new histograms should overload this. */
    virtual void FillOutputList(TList* list);

    /** Initialise the analysis, must be overloaded. */
    virtual void Init();

    /** 
    * Creates the histograms, must be overloaded if you want to add your own. 
    * Uses the fHistogramNameSuffix to create proper histogram names
    */
    virtual void CreateHistograms();
    virtual void CreateTrees();
    
    /** Fills the histograms, must be overloaded if you want to add your own */
    virtual void FillHistograms();

    /** Reset event specific values (Et etc.) */
    virtual void ResetEventValues();

    /** Set Particle codes/mass */
    virtual void SetParticleCodes();
    
    /** Cuts info */
    AliAnalysisEtCuts * GetCuts() const { return fCuts; } 
    virtual void SetCuts(const AliAnalysisEtCuts *cuts) 
    { fCuts = (AliAnalysisEtCuts *) cuts; } 

    /** Total Et in the event (without acceptance cuts) */
    Double_t GetTotEt() const { return fTotEt; }

    /** Total Et in the event within the acceptance cuts */
    Double_t GetTotEtAcc() const { return fTotEtAcc; }

   /** Total neutral Et in the event (without acceptance cuts) */
    Double_t GetTotNeutralEt() const { return fTotNeutralEt; }

    /** Total neutral Et in the event within the acceptance cuts */
    Double_t GetTotNeutralEtAcc() const { return fTotNeutralEtAcc; }
    
    /** Total charged Et in the event (without acceptance cuts) */
    Double_t GetTotChargedEt() const { return fTotChargedEt; }

    /** Total charged Et in the event within the acceptance cuts */
    Double_t GetTotChargedEtAcc() const { return fTotChargedEtAcc; }

    void SetTPCOnlyTrackCuts(const AliESDtrackCuts *cuts){ fEsdtrackCutsTPC = (AliESDtrackCuts *) cuts;}

protected:
       
    TString fHistogramNameSuffix; /** The suffix for the histogram names */

    AliAnalysisEtCuts *fCuts; // keeper of basic cuts

    /** PDG Database */
    TDatabasePDG *fPdgDB;//data base used for looking up pdg codes
    //these codes are stored as variables because otherwise there were issues using this with the plugin
    Int_t fPiPlusCode;//pdg pi plus code
    Int_t fPiMinusCode;//pdg pi minus code
    Int_t fKPlusCode;// pdg k plus code
    Int_t fKMinusCode;//pdg k minus code
    Int_t fProtonCode;//pdg proton code
    Int_t fAntiProtonCode;//pdg antiproton code
    Int_t fLambdaCode;// pdg lambda code
    Int_t fAntiLambdaCode;//pdg antilambda code
    Int_t fK0SCode;//pdg k0 short code
    Int_t fOmegaCode;//pdg omega code
    Int_t fAntiOmegaCode;//pdg anti-omega code
    Int_t fXi0Code;//pdg xi-0 code
    Int_t fAntiXi0Code;//pdg anti-xi0 code
    Int_t fXiCode;//pdg xi code
    Int_t fAntiXiCode;//pdg anti-xi code
    Int_t fSigmaCode;//pdg sigma code
    Int_t fAntiSigmaCode;//pdg anti-sigma code
    Int_t fK0LCode;//pdg k0 long code
    Int_t fNeutronCode;//pdg neutron code
    Int_t fAntiNeutronCode;//pdg anti-neutron code
    Int_t fEPlusCode;//pdg positron code
    Int_t fEMinusCode;//pdg electron code
    Int_t fMuPlusCode; // pdg muon + code
    Int_t fMuMinusCode; // pdg muon - code
    Int_t fGammaCode; // pdg gamma code
    Float_t fPionMass;//pdg pion mass

    Double_t fTotEt;/** Total Et in the event (without acceptance cuts) */    
    Double_t fTotEtAcc;/** Total Et in the event within the acceptance cuts */
    
    Double_t fTotNeutralEt;/** Total neutral Et in the event */    
    Double_t fTotNeutralEtAcc;/** Total neutral Et in the event within the acceptance cuts */    
    Double_t fTotChargedEt;/** Total charged Et in the event */    
    Double_t fTotChargedEtAcc;/** Total charged Et in the event within the acceptance cuts */

    Int_t fMultiplicity;/** Multiplicity of particles in the event */    
    Int_t fChargedMultiplicity;/** Multiplicity of charged particles in the event */    
    Int_t fNeutralMultiplicity; /** Multiplicity of neutral particles in the event */
    
    Double_t fBaryonEt;     /** Et of identified baryons; calo based (Rec only for now) */    
    Double_t fAntiBaryonEt; /** Et of identified anti-baryons; calo based (Rec only for now) */
    Double_t fMesonEt;     /** Et of identified mesons; calo based (Rec only for now) */

    Double_t fProtonEt; /** Et of identified protons */
    Double_t fPionEt; /** Et of identified pions */
    Double_t fChargedKaonEt; /** Et of identified charged kaons */
    Double_t fMuonEt; /** Et of identified muons */
    Double_t fElectronEt; /** Et of identified electrons */
    Double_t fNeutronEt; /** Et of neutrons (MC only for now) */
    Double_t fAntiNeutronEt; /** Et of anti-neutrons (MC only for now) */
    Double_t fGammaEt; /** Et of identified electrons (MC only for now) */
   
    Double_t fProtonEtAcc; /** Et of identified protons in calorimeter acceptance */    
    Double_t fPionEtAcc; /** Et of identified pions in calorimeter acceptance */    
    Double_t fChargedKaonEtAcc; /** Et of identified charged kaons in calorimeter acceptance */    
    Double_t fMuonEtAcc; /** Et of identified muons in calorimeter acceptance */
    Double_t fElectronEtAcc; /** Et of identified electrons in calorimeter acceptance */
    
    Float_t fEnergyDeposited; /** Energy deposited in calorimeter */
    Float_t fEnergyTPC; /** Energy measured in TPC */
    Short_t fCharge; /** Charge of the particle */
    Short_t fParticlePid; /** Particle PID */
    Float_t fPidProb; /** Probability of PID */
    Bool_t fTrackPassedCut; /** The track is accepted by ESDTrackCuts */
   
        
    Double_t fEtaCut;/** Cut in eta (standard |eta| < 0.5 )*/

    /** Eta cut for our acceptance */
    Double_t fEtaCutAcc; // Eta cut for our acceptance
    
    /** Min phi cut for our acceptance in radians */    
    Double_t fPhiCutAccMin; // Min phi cut for our acceptance in radians     
    
    /** Max phi cut for our acceptance in radians */
    Double_t fPhiCutAccMax; // Max phi cut for our acceptance in radians 
    
    /** Detector radius */
    Double_t fDetectorRadius; // Detector radius 
    
    /** Cut on the cluster energy */    
    Double_t fClusterEnergyCut; // Cut on the cluster energy 
    
    /** Minimum energy to cut on single cell cluster */
    Double_t fSingleCellEnergyCut;  // Minimum energy to cut on single cell cluster

    // Declare the histograms

    /** The full Et spectrum measured */
    TH1F *fHistEt; //Et spectrum

    /** The full charged Et spectrum measured */
    TH1F *fHistChargedEt; //Charged Et spectrum 

    /** The full neutral Et spectrum measured */
    TH1F *fHistNeutralEt; //Neutral Et spectrum

    /** The Et spectrum within the calorimeter acceptance */
    TH1F *fHistEtAcc; //Et in acceptance

    /** The charged Et spectrum within the calorimeter acceptance */
    TH1F *fHistChargedEtAcc; //Charged Et in acceptance

    /** The neutral Et spectrum within the calorimeter acceptance */
    TH1F *fHistNeutralEtAcc; //Et in acceptance

    /** Multiplicity of particles in the events */
    TH1F *fHistMult; //Multiplicity

    /** Charged multiplicity of particles in the events */
    TH1F *fHistChargedMult; //Charged multiplicity

    /** Neutral multiplicity of particles in the events */
    TH1F *fHistNeutralMult; //Neutral multiplicity

    /* Acceptance plots */
    TH2F *fHistPhivsPtPos; //phi vs pT plot for positive tracks
    TH2F *fHistPhivsPtNeg; //phi vs pT plot for negative tracks

    /* PID plots */
    TH1F *fHistBaryonEt; /** Et of identified baryons */    
    TH1F *fHistAntiBaryonEt; /** Et of identified anti-baryons */
    TH1F *fHistMesonEt; /** Et of identified mesons */

    TH1F *fHistProtonEt; /** Et of identified protons */
    TH1F *fHistPionEt; /** Et of identified protons */
    TH1F *fHistChargedKaonEt; /** Et of identified charged kaons */
    TH1F *fHistMuonEt; /** Et of identified muons */
    TH1F *fHistElectronEt; /** Et of identified electrons */
    TH1F *fHistNeutronEt; /** Et of neutrons (MC only for now) */
    TH1F *fHistAntiNeutronEt; /** Et of anti-neutrons (MC only for now) */
    TH1F *fHistGammaEt; /** Et of gammas (MC only for now) */
    
    TH1F *fHistProtonEtAcc; /** Et of identified protons in calorimeter acceptance */    
    TH1F *fHistPionEtAcc; /** Et of identified protons in calorimeter acceptance */    
    TH1F *fHistChargedKaonEtAcc; /** Et of identified charged kaons in calorimeter acceptance */    
    TH1F *fHistMuonEtAcc; /** Et of identified muons in calorimeter acceptance */
    TH1F *fHistElectronEtAcc; /** Et of identified electrons in calorimeter acceptance */
    
    /* Correction plots */
    TH1F *fHistTMDeltaR; /* Track matching plots; Rec only for now */

    TTree *fTree; // optional TTree
    TTree *fTreeDeposit; // optional TTree for energy deposit measurements

    AliESDtrackCuts* fEsdtrackCutsTPC;//esd track cuts for TPC tracks (which may also contain ITS hits)

private:
    //Declare private to avoid compilation warning
    AliAnalysisEt & operator = (const AliAnalysisEt & g) ;//cpy assignment
    AliAnalysisEt(const AliAnalysisEt & g) ; // cpy ctor

    ClassDef(AliAnalysisEt, 1);
};

#endif // ALIANALYSISET_H
