#ifndef ALIANALYSISET_H
#define ALIANALYSISET_H

#include "TString.h"
#include "TDatabasePDG.h"
#include "Rtypes.h"

class TH2F;
class TH1F;
class AliVEvent;
class TList;

class AliAnalysisEt
{
public:
   
    AliAnalysisEt();
    virtual ~AliAnalysisEt();

private:
  //Declare it private to avoid compilation warning
  AliAnalysisEt & operator = (const AliAnalysisEt & g) ;//cpy assignment
  AliAnalysisEt(const AliAnalysisEt & g) ; // cpy ctor
  
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
    
    /** Sum of the total Et for all events */
    Double_t GetSumEt() { return fSumEt; }

    /** Sum of the total Et within our acceptance for all events */
    Double_t GetSumEtAcc() { return fSumEtAcc; }

    /** Total Et in the event (without acceptance cuts) */
    Double_t GetTotEt() { return fTotEt; }

    /** Total Et in the event within the acceptance cuts */
    Double_t GetTotEtAcc() { return fTotEtAcc; }

   /** Total neutral Et in the event (without acceptance cuts) */
    Double_t GetTotNeutralEt() { return fTotNeutralEt; }

    /** Total neutral Et in the event within the acceptance cuts */
    Double_t GetTotNeutralEtAcc() { return fTotNeutralEtAcc; }
    
    /** Total charged Et in the event (without acceptance cuts) */
    Double_t GetTotChargedEt() { return fTotChargedEt; }

    /** Total charged Et in the event within the acceptance cuts */
    Double_t GetTotChargedEtAcc() { return fTotChargedEtAcc; }


protected:
   
    /** The suffix for the histogram names */
    TString fHistogramNameSuffix;

    /** PDG Database */
    TDatabasePDG *fPdgDB;

    Int_t PiPlusCode;
    Int_t PiMinusCode;
    Int_t KPlusCode;
    Int_t KMinusCode;
    Int_t ProtonCode;
    Int_t AntiProtonCode;
    Int_t LambdaCode;
    Int_t AntiLambdaCode;
    Int_t K0SCode;
    Int_t OmegaCode;
    Int_t AntiOmegaCode;
    Int_t Xi0Code;
    Int_t AntiXi0Code;
    Int_t XiCode;
    Int_t AntiXiCode;
    Int_t SigmaCode;
    Int_t AntiSigmaCode;
    Int_t K0LCode;
    Int_t NeutronCode;
    Int_t AntiNeutronCode;
    Int_t EPlusCode;
    Int_t EMinusCode;
    Float_t PionMass;

    /** Sum of the total Et for all events */
    Double_t fSumEt;

    /** Sum of the total Et within our acceptance for all events */
    Double_t fSumEtAcc;

    /** Total Et in the event (without acceptance cuts) */
    Double_t fTotEt;

    /** Total Et in the event within the acceptance cuts */
    Double_t fTotEtAcc;

    /** Total neutral Et in the event */
    Double_t fTotNeutralEt;

    /** Total neutral Et in the event within the acceptance cuts */
    Double_t fTotNeutralEtAcc;

    /** Total charged Et in the event */
    Double_t fTotChargedEt;

    /** Total charged Et in the event within the acceptance cuts */
    Double_t fTotChargedEtAcc;

    /** Multiplicity of particles in the event */
    Int_t fMultiplicity;
    
    /** Multiplicity of charged particles in the event */
    Int_t fChargedMultiplicity;
    
    /** Multiplicity of neutral particles in the event */
    Int_t fNeutralMultiplicity; 
    
    /** Cut in eta ( normally |eta| < 0.5 */
    Double_t fEtaCut;

    /** Eta cut for our acceptance */
    Double_t fEtaCutAcc;

    /** Min phi cut for our acceptance in radians */
    Double_t fPhiCutAccMin;

    /** Max phi cut for our acceptance in radians */
    Double_t fPhiCutAccMax;

    /** Detector radius */
    Double_t fDetectorRadius;

    /** Vertex cuts */
    Double_t fVertexXCut;
    Double_t fVertexYCut;
    Double_t fVertexZCut;

    /** Impact parameter cuts */
    Double_t fIPxyCut;
    Double_t fIPzCut;

    /** Cut on the cluster energy */
    Double_t fClusterEnergyCut;

    /** Cut on track pt */
    Double_t fTrackPtCut;

    /** Minimum energy to cut on single cell cluster */
    Double_t fSingleCellEnergyCut;
    
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
    TH1F *fHistBaryonEt;
    TH1F *fHistAntiBaryonEt;
    TH1F *fHistMesonEt;

    TH1F *fHistBaryonEtAcc;
    TH1F *fHistAntiBaryonEtAcc;
    TH1F *fHistMesonEtAcc;

    /* Correction plots */
    TH2F *fHistEtRecvsEtMC; //Reconstructed Et versus MC Et

    /* Track matching plots */
    TH1F *fHistTMDeltaR;

    ClassDef(AliAnalysisEt, 0);
};

#endif // ALIANALYSISET_H
