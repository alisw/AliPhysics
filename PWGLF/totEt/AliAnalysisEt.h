//_________________________________________________________________________
//  Utility Class for transverse energy studies
//  Base class for ESD & MC analysis
//  - reconstruction and MonteCarlo output
//
//*-- Authors: Oystein Djuvsland (Bergen), David Silvermyr (ORNL)
//_________________________________________________________________________

#ifndef ALIANALYSISET_H
#define ALIANALYSISET_H

class AliCentrality;
#include "AliAnalysisEtCommon.h"
#include "THnSparse.h"

class TString;
class TTree;
class TH2F;
class TH1F;
class AliVEvent;
class TList;
class TString;
class AliESDtrackCuts;
class AliAnalysisEtCuts;
class AliESDCaloCluster;
//class THnSparseD;

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
    
    /** Set the centrality object */
    void SetCentralityObject(AliCentrality *cent) { fCentrality = cent; }
    
    /** Get contribution from non-removed charged particles */
    virtual Double_t GetChargedContribution(Int_t /*clusterMultiplicity*/) {return 0;}

    /** Get contribution from non-removed neutral particles */
    virtual Double_t GetNeutralContribution(Int_t /*clusterMultiplicity*/) {return 0;}
    
    /** Get contribution from removed gammas */
    virtual Double_t GetGammaContribution(Int_t /*clusterMultiplicity*/) {return 0;}

    void MakeSparseHistograms(){fMakeSparse=kTRUE;}

protected:

    //AliAnalysisEtCuts *fCuts; // keeper of basic cuts
    Double_t CalculateTransverseEnergy(AliESDCaloCluster *cluster);

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

    Int_t fCentClass; // centrality class
        
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
    
    Double_t fTrackDistanceCut; // cut on track distance    
    
    Double_t fTrackDxCut; // cut on track distance in x
    
    Double_t fTrackDzCut; // cut on track distance in z
    

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
    TH2F *fHistTMDxDz; /* Track matching plots; Rec only for now */
  
    /* Auxiliary Histogram variables */
    static Float_t fgEtaAxis[17];//bins for eta axis of histograms
    static Int_t fgnumOfEtaBins;//number of eta bins
    static Float_t fgPtAxis[117];//bins for pt axis of histograms
    static Int_t fgNumOfPtBins;//number of pt bins
    static Float_t fgEAxis[79];//bins for pt axis of histograms
    static Int_t fgNumOfEBins;//number of pt bins
    static Float_t fgRAxis[48];//bins for R axis
    static Int_t fgNumOfRBins;//number of R bins
	
	
    TTree *fTree; // optional TTree
    TTree *fTreeDeposit; // optional TTree for energy deposit measurements

   /** Centrality object */
    AliCentrality *fCentrality; //Centrality object
    
    Short_t fDetector;     /** Which detector? (-1 -> PHOS, 1 -> EMCAL)*/

    Bool_t fMakeSparse;//Boolean for whether or not to make sparse histograms
    
    THnSparseF *fSparseHistTracks;     /** THnSparse histograms */
    
    THnSparseF *fSparseHistClusters;     /** THnSparse histograms */
    
    /** ET sparse valuses */
    THnSparseF *fSparseHistEt; //!
       
    /** Values for sparse hists */
    Double_t *fSparseTracks; //!
    
    /** Values for sparse hists */
    Double_t *fSparseClusters; //!
    
    /** ET sparse valuses */
    Double_t *fSparseEt; //!
    
    
    

private:
    //Declare private to avoid compilation warning
    AliAnalysisEt & operator = (const AliAnalysisEt & g) ;//cpy assignment
    AliAnalysisEt(const AliAnalysisEt & g) ; // cpy ctor

    ClassDef(AliAnalysisEt, 1);
};

#endif // ALIANALYSISET_H
