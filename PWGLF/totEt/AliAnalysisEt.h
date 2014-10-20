#ifndef ALIANALYSISET_H
#define ALIANALYSISET_H
//_________________________________________________________________________
//  Utility Class for transverse energy studies
//  Base class for ESD & MC analysis
//  - reconstruction and MonteCarlo output
//
//*-- Authors: Oystein Djuvsland (Bergen), David Silvermyr (ORNL)
//_________________________________________________________________________

#include "AliAnalysisEtCommon.h"
#include "THnSparse.h"
#include "AliESDCaloCluster.h"
#include "AliAnalysisEtCuts.h"
#include "AliAnalysisEtTrackMatchCorrections.h"
#include <vector>
#include "Rtypes.h"
#include "AliAnalysisEtSelector.h"
#include "AliAnalysisEtSelectorEmcal.h"
#include "AliAnalysisEtSelectorPhos.h"

class AliAnalysisEtRecEffCorrection;
class AliAnalysisEtTrackMatchCorrections;
class AliAnalysisEtSelector;
class AliCentrality;
class TString;
class TTree;
class TH2F;
class TH1F;
class TH1I;
class AliVEvent;
class TList;
class TString;
class AliESDtrackCuts;
class AliAnalysisEtCuts;
class AliESDCaloCluster;
//class THnSparseD;
class AliPIDResponse;

class AliAnalysisEt : public AliAnalysisEtCommon
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
    TH2F* CreateEtaEHisto2D(TString name, TString title, TString ztitle);
    TH2F* CreateEtaPtHisto2D(TString name, TString title, TString ztitle);
    TH2F* CreateEtaEtHisto2D(TString name, TString title, TString ztitle);
    TH2F* CreateResEHisto2D(TString name, TString title, TString ztitle);
    TH2F* CreateResPtHisto2D(TString name, TString title, TString ztitle);
    THnSparseF* CreateClusterHistoSparse(TString name, TString title);
    THnSparseF* CreateNeutralPartHistoSparse(TString name, TString title);
    THnSparseF* CreateChargedPartHistoSparse(TString name, TString title);

    /** Fills the histograms, must be overloaded if you want to add your own */
    virtual void FillHistograms();

    /** Reset event specific values (Et etc.) */
    virtual void ResetEventValues();

    /** Total Et in the event (without acceptance cuts) */
    Double_t GetTotEt() const {
        return fTotEt;
    }

    /** Total neutral Et in the event (without acceptance cuts) */
    Double_t GetTotNeutralEt() const {
        return fTotNeutralEt;
    }

    /** Total charged Et in the event (without acceptance cuts) */
    Double_t GetTotChargedEt() const {
        return fTotChargedEt;
    }

    void SetTPCOnlyTrackCuts(const AliESDtrackCuts *cuts) {
        fEsdtrackCutsTPC = (AliESDtrackCuts *) cuts;
    }

    /** Set the centrality object */
    void SetCentralityObject(AliCentrality *cent) {
        fCentrality = cent;
    }

    /** Get min ET correction */
    Double_t GetMinEtCorrection(Int_t clusterMultiplicity) {
        return fTmCorrections->GetMinEtCorrection(clusterMultiplicity);
    }
    /** Get contribution from non-removed neutrons */
    Double_t GetNeutronContribution(Int_t clusterMultiplicity) {
        return fTmCorrections->GetNeutronCorrection(clusterMultiplicity);
    }
    /** Get contribution from non-removed hadrons */
    Double_t GetHadronContribution(Int_t clusterMultiplicity) {
        return fTmCorrections->GetHadronCorrection(clusterMultiplicity);
    }
    /** Get contribution from non-removed kaons */
    Double_t GetKaonContribution(Int_t clusterMultiplicity) {
        return fTmCorrections->GetKaonCorrection(clusterMultiplicity);
    }//hadron
    /** Get contribution from non-removed secondarys */
    Double_t GetSecondaryContribution(Int_t clusterMultiplicity) {
        return fTmCorrections->GetSecondaryCorrection(clusterMultiplicity);
    }//hadron


    /** Get contribution from non-removed charged particles */
    Double_t GetChargedContribution(Int_t clusterMultiplicity) {
        return fTmCorrections->ChargedContr(clusterMultiplicity);
    }

    /** Get contribution from non-removed neutral particles */
    Double_t GetNeutralContribution(Int_t clusterMultiplicity) {
        return fTmCorrections->NeutralContr(clusterMultiplicity);
    }

    /** Get contribution from removed gammas */
    Double_t GetGammaContribution(Int_t clusterMultiplicity) {
        return fTmCorrections->GammaContr(clusterMultiplicity);
    }

    void MakeSparseHistograms() {
        fMakeSparse=kTRUE;
    }
    
    AliAnalysisEtCuts * GetCuts() const { return fCuts; }
    AliAnalysisEtSelector *GetSelector() const {return fSelector;}
    

    // Read in corrections
    Int_t ReadCorrections(TString filename);  // Read in corrections
    
    void SetFsub(Float_t val){fsub=val;};//function for setting fsub for EMCal/PHOS hadronic corrections
    void SetFsubForMeanHadE(Float_t val){fsubmeanhade=val;};//function for setting fsub for EMCal/PHOS hadronic corrections

protected:

    //AliAnalysisEtCuts *fCuts; // keeper of basic cuts
    
    // Return corrected cluster E_T
    Double_t CorrectForReconstructionEfficiency(const AliESDCaloCluster &cluster,Int_t cent = 0);
    Double_t CorrectForReconstructionEfficiency(const AliESDCaloCluster &cluster, Float_t eReco,Int_t cent = 0);
    
    // Track matching (hadrdonic contamination) corrections
    AliAnalysisEtTrackMatchCorrections *fTmCorrections;
    
    // Reconstruction efficiency corrections
    AliAnalysisEtRecEffCorrection *fReCorrections;
    
    TTree *fEventSummaryTree; //! Contains event level information

    TTree *fAcceptedTree; //! Tree for information about accepted particles
    
    TTree *fDepositTree; //! optional TTree for energy deposit measurements
    
    Double_t fTotEt;/** Total Et in the event (without acceptance cuts) */

    Double_t fTotEtAcc;/** Total Et in the event (without acceptance cuts) */

    Double_t fTotNeutralEt;/** Total neutral Et in the event */

    Double_t fTotNeutralEtAcc;/** Total neutral Et in the event */

    Double_t fTotChargedEt;/** Total charged Et in the event */

    Double_t fTotChargedEtAcc;/** Total charged Et in the event */

    Int_t fMultiplicity;/** Multiplicity of particles in the event */
    Int_t fChargedMultiplicity;/** Multiplicity of charged particles in the event */
    Int_t fNeutralMultiplicity; /** Multiplicity of neutral particles in the event */

    Double_t fProtonEt; /** Et of identified protons */
    Double_t fAntiProtonEt; /** Et of identified protons */

    Double_t fNeutronEt; /** Et of neutrons (MC only for now) */
    Double_t fAntiNeutronEt; /** Et of anti-neutrons (MC only for now) */
    
    Double_t fPi0Et; // Et of identified pi0
    Double_t fPiPlusEt; // Et of identified pi+
    Double_t fPiMinusEt; // Et of identified pi-
    
    Double_t fKPlusEt; // Et of identified K+ 
    Double_t fKMinusEt; // Et of identified K- 
    Double_t fK0sEt; // Et of identified K0 short
    Double_t fK0lEt; // Et of identified K0 long
    
    Double_t fMuMinusEt; // Et of identified mu- 
    Double_t fMuPlusEt; // Et of identified mu+ 

    Double_t fEMinusEt; // Et of identified e-
    Double_t fEPlusEt; // Et of identified e+
    
    Double_t fGammaEt; /** Et of identified electrons (MC only for now) */

    Double_t fProtonRemovedEt; /** Et of identified protons */
    Double_t fAntiProtonRemovedEt; /** Et of identified protons */

    Double_t fNeutronRemovedEt; /** Et of neutrons (MC only for now) */
    Double_t fAntiNeutronRemovedEt; /** Et of anti-neutrons (MC only for now) */
    
    Double_t fPi0RemovedEt; // Removed Et of identified pi0
    Double_t fPiPlusRemovedEt; // Removed Et of identified pi+
    Double_t fPiMinusRemovedEt; // Removed Et of identified pi-

    Double_t fKPlusRemovedEt; // Removed Et of identified K+ 
    Double_t fKMinusRemovedEt; // Removed Et of identified K- 
    Double_t fK0sRemovedEt; // Removed Et of identified K0 short
    Double_t fK0lRemovedEt; // Removed Et of identified K0 long
    
    Double_t fMuMinusRemovedEt; // Removed Et of identified mu- 
    Double_t fMuPlusRemovedEt; // Removed Et of identified mu+ 
    
    Double_t fEMinusRemovedEt; // Removed Et of identified e-
    Double_t fEPlusRemovedEt; // Removed Et of identified e+
    
    Double_t fGammaRemovedEt; /** Removed Et of identified electrons (MC only for now) */

    Double_t fProtonMult; /** Mult of identified protons */
    Double_t fAntiProtonMult; /** Mult of identified protons */

    Double_t fNeutronMult; /** Mult of neutrons (MC only for now) */
    Double_t fAntiNeutronMult; /** Mult of anti-neutrons (MC only for now) */
    
    Double_t fPi0Mult; // Mult of identified pi0
    Double_t fPiPlusMult; // Mult of identified pi+
    Double_t fPiMinusMult; // Mult of identified pi-
    
    Double_t fKPlusMult; // Mult of identified K+ 
    Double_t fKMinusMult; // Mult of identified K- 
    Double_t fK0sMult; // Mult of identified K0 short
    Double_t fK0lMult; // Mult of identified K0 long
    
    Double_t fMuMinusMult; // Mult of identified mu- 
    Double_t fMuPlusMult; // Mult of identified mu+ 

    Double_t fEMinusMult; // Mult of identified e-
    Double_t fEPlusMult; // Mult of identified e+
    
    Double_t fGammaMult; /** Mult of identified electrons (MC only for now) */

    Double_t fProtonRemovedMult; /** Mult of identified protons */
    Double_t fAntiProtonRemovedMult; /** Mult of identified protons */

    Double_t fNeutronRemovedMult; /** Mult of neutrons (MC only for now) */
    Double_t fAntiNeutronRemovedMult; /** Mult of anti-neutrons (MC only for now) */
    
    Double_t fPi0RemovedMult; // Removed Mult of identified pi0
    Double_t fPiPlusRemovedMult; // Removed Mult of identified pi+
    Double_t fPiMinusRemovedMult; // Removed Mult of identified pi-

    Double_t fKPlusRemovedMult; // Removed Mult of identified K+ 
    Double_t fKMinusRemovedMult; // Removed Mult of identified K- 
    Double_t fK0sRemovedMult; // Removed Mult of identified K0 short
    Double_t fK0lRemovedMult; // Removed Mult of identified K0 long
    
    Double_t fMuMinusRemovedMult; // Removed Mult of identified mu- 
    Double_t fMuPlusRemovedMult; // Removed Mult of identified mu+ 
    
    Double_t fEMinusRemovedMult; // Removed Mult of identified e-
    Double_t fEPlusRemovedMult; // Removed Mult of identified e+
    
    Double_t fGammaRemovedMult; /** Removed Mult of identified electrons (MC only for now) */

    Float_t fEnergyDeposited; /** Energy deposited in calorimeter */
    Float_t fMomentumTPC; /** Momentum measured in TPC */
    Short_t fCharge; /** Charge of the particle */
    Short_t fParticlePid; /** Particle PID */
    Float_t fPidProb; /** Probability of PID */
    Bool_t fTrackPassedCut; /** The track is accepted by ESDTrackCuts */

    Int_t fCentClass; // centrality class

    /** Detector radius */
    Double_t fDetectorRadius; // Detector radius

    /** Minimum energy to cut on single cell cluster */
    Double_t fSingleCellEnergyCut;  // Minimum energy to cut on single cell cluster

    Double_t fChargedEnergyRemoved; // Charged energy removed
    Double_t fNeutralEnergyRemoved; // Neutral energy removed
    Double_t fGammaEnergyAdded; // gamma energy added

    // Declare the histograms

    /** The EM Et spectrum measured */
    TH1F *fHistEt; //!Et spectrum

    /** Multiplicity of neutral particles in the events */
    TH1F *fHistNeutralMult; //!Multiplicity

    // Acceptance plots 
    TH2F *fHistPhivsPtPos; //!phi vs pT plot for positive tracks
    TH2F *fHistPhivsPtNeg; //!phi vs pT plot for negative tracks

    /* Auxiliary Histogram variables */
    static Float_t fgEtaAxis[17];//bins for eta axis of histograms
    static Int_t fgnumOfEtaBins;//number of eta bins
    static Float_t fgPtAxis[117];//bins for pt axis of histograms
    static Int_t fgNumOfPtBins;//number of pt bins
    static Float_t fgEAxis[79];//bins for pt axis of histograms
    static Int_t fgNumOfEBins;//number of pt bins
    static Float_t fgRAxis[48];//bins for R axis
    static Int_t fgNumOfRBins;//number of R bins

    /** Centrality object */
    AliCentrality *fCentrality; //!Centrality object

    Bool_t fMakeSparse;//Boolean for whether or not to make sparse histograms

    TH1I *fCutFlow; //! Cut flow
    
    AliAnalysisEtSelector *fSelector; // Selector class

    AliPIDResponse *fPIDResponse;//

    Float_t fsub;
    Float_t fsubmeanhade;
    


private:
   
  
    
    //Declare private to avoid compilation warning
    AliAnalysisEt & operator = (const AliAnalysisEt & g) ;//cpy assignment
    AliAnalysisEt(const AliAnalysisEt & g) ; // cpy ctor

    ClassDef(AliAnalysisEt, 3);
};

#endif // ALIANALYSISET_H
