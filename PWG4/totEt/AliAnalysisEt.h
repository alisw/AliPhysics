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

class TH2F;
class TH1F;
class AliVEvent;
class TList;
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

    /** Sum of the total Et for all events */
    Double_t GetSumEt() const { return fSumEt; }

    /** Sum of the total Et within our acceptance for all events */
    Double_t GetSumEtAcc() const { return fSumEtAcc; }

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
    Float_t fPionMass;//pdg pion mass

    Double_t fSumEt;/** Sum of the total Et for all events */
    Double_t fSumEtAcc;/** Sum of the total Et within our acceptance for all events */    
    Double_t fTotEt;/** Total Et in the event (without acceptance cuts) */    
    Double_t fTotEtAcc;/** Total Et in the event within the acceptance cuts */
    
    Double_t fTotNeutralEt;/** Total neutral Et in the event */    
    Double_t fTotNeutralEtAcc;/** Total neutral Et in the event within the acceptance cuts */    
    Double_t fTotChargedEt;/** Total charged Et in the event */    
    Double_t fTotChargedEtAcc;/** Total charged Et in the event within the acceptance cuts */

    Int_t fMultiplicity;/** Multiplicity of particles in the event */    
    Int_t fChargedMultiplicity;/** Multiplicity of charged particles in the event */    
    Int_t fNeutralMultiplicity; /** Multiplicity of neutral particles in the event */

    Double_t fEtaCutAcc;/** Eta cut for our acceptance */
    Double_t fPhiCutAccMin; /** Min phi cut for our acceptance in radians */    
    Double_t fPhiCutAccMax; /** Max phi cut for our acceptance in radians */
    Double_t fDetectorRadius; /** Detector radius */
    
    Double_t fClusterEnergyCut; /** Cut on the cluster energy */    
    Double_t fSingleCellEnergyCut;  /** Minimum energy to cut on single cell cluster */

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

    /* PID plots, Et */
    TH1F *fHistBaryonEt; // baryon
    TH1F *fHistAntiBaryonEt; // anti-baryon
    TH1F *fHistMesonEt; // meson

    TH1F *fHistBaryonEtAcc; // baryon, acc
    TH1F *fHistAntiBaryonEtAcc; // anti-baryon, acc
    TH1F *fHistMesonEtAcc; // meson, acc

    TH1F *fHistTMDeltaR;  /* Track matching plots */

private:
    //Declare private to avoid compilation warning
    AliAnalysisEt & operator = (const AliAnalysisEt & g) ;//cpy assignment
    AliAnalysisEt(const AliAnalysisEt & g) ; // cpy ctor

    ClassDef(AliAnalysisEt, 0);
};

#endif // ALIANALYSISET_H
