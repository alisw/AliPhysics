#ifndef AliAnalysisTaskPhiCorrelations_H
#define AliAnalysisTaskPhiCorrelations_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////////////////////////////
//
// Analysis class for Underlying Event studies w.r.t. leading track
//
// Look for correlations on the tranverse regions w.r.t
// the leading track in the event
//
// This class needs input AODs.
// The output is a list of analysis-specific containers.
//
// The AOD can be either connected to the InputEventHandler  
// for a chain of AOD files 
// or 
// to the OutputEventHandler
// for a chain of ESD files,
// in this case the class should be in the train after the jet-finder
//
//    Authors:
//    Jan Fiete Grosse-Oetringhaus
// 
////////////////////////////////////////////////////////////////////////

#include "AliAnalysisTask.h"
#include "AliUEHist.h"
#include "TString.h"
#include "AliVParticle.h"
#include "AliLog.h"
#include "THn.h" // in cxx file causes .../THn.h:257: error: conflicting declaration ‘typedef class THnT<float> THnF’

class AliAODEvent;
class AliAnalyseLeadingTrackUE;
class AliInputEventHandler;
class AliMCEvent;
class AliMCEventHandler;
class AliUEHistograms;
class AliVParticle;
class TH1D;
class TObjArray;
class AliEventPoolManager;
class AliESDEvent;
class AliHelperPID;


class  AliAnalysisTaskPhiCorrelations : public AliAnalysisTask
  {
  public:
    AliAnalysisTaskPhiCorrelations(const char* name="AliAnalysisTaskPhiCorrelations");
    virtual           ~AliAnalysisTaskPhiCorrelations();
       
      
    // Implementation of interace methods
    virtual     void   ConnectInputData(Option_t *);
    virtual     void   CreateOutputObjects();
    virtual     void   Exec(Option_t *option);

    // Setters/Getters
    // general configuration
    virtual     void    SetDebugLevel( Int_t level )  { fDebug = level; }
    virtual     void    SetMode(Int_t mode)           { fMode  = mode;  }
    virtual     void    SetReduceMemoryFootprint(Bool_t flag) { fReduceMemoryFootprint = flag; }
    virtual	void	SetEventMixing(Bool_t flag) { fFillMixed = flag; }
    virtual	void    SetMixingTracks(Int_t tracks) { fMixingTracks = tracks; }
    virtual	void	SetCompareCentralities(Bool_t flag) { fCompareCentralities = flag; }
    virtual	void	SetTwoTrackEfficiencyStudy(Bool_t flag) { fTwoTrackEfficiencyStudy = flag; }
    virtual	void	SetTwoTrackEfficiencyCut(Float_t value = 0.02) { fTwoTrackEfficiencyCut = value; }
    virtual	void	SetUseVtxAxis(Int_t flag) { fUseVtxAxis = flag; }
    virtual	void	SetCourseCentralityBinning(Bool_t flag) { fCourseCentralityBinning = flag; }
    virtual     void    SetSkipTrigger(Bool_t flag) { fSkipTrigger = flag; }
    virtual     void    SetInjectedSignals(Bool_t flag) { fInjectedSignals = flag; }
    
    // histogram settings
    void SetEfficiencyCorrection(THnF* hist, Bool_t correctTriggers) { fEfficiencyCorrection = hist; fCorrectTriggers = correctTriggers; }

    // for event QA
    void   SetTracksInVertex( Int_t val ){ fnTracksVertex = val; }
    void   SetZVertex( Double_t val )    { fZVertex = val; }
    
    // track cuts
    void   SetTrackEtaCut( Double_t val )    { fTrackEtaCut = val; }
    void   SetOnlyOneEtaSide(Int_t flag)    { fOnlyOneEtaSide = flag; }
    void   SetPtMin(Double_t val)            { fPtMin = val; }
    void   SetFilterBit( UInt_t val )        { fFilterBit = val;  }
    
    void   SetEventSelectionBit( UInt_t val )        { fSelectBit = val;  }
    void   SetUseChargeHadrons( Bool_t val ) { fUseChargeHadrons = val; }
    void   SetSelectParticleSpecies( Int_t trigger, Int_t associated ) { fParticleSpeciesTrigger = trigger; fParticleSpeciesAssociated = associated; }
    void   SetSelectCharge(Int_t selectCharge) { fSelectCharge = selectCharge; }
    void   SetSelectTriggerCharge(Int_t selectCharge) { fTriggerSelectCharge = selectCharge; }
    void   SetSelectAssociatedCharge(Int_t selectCharge) { fAssociatedSelectCharge = selectCharge; }
    void   SetTriggerRestrictEta(Float_t eta) { fTriggerRestrictEta = eta; }
    void   SetEtaOrdering(Bool_t flag) { fEtaOrdering = flag; }
    void   SetPairCuts(Bool_t conversions, Bool_t resonances) { fCutConversions = conversions; fCutResonances = resonances; }
    void   SetCentralityMethod(const char* method) { fCentralityMethod = method; }
    void   SetFillpT(Bool_t flag) { fFillpT = flag; }
    void   SetStepsFillSkip(Bool_t step0, Bool_t step6) { fFillOnlyStep0 = step0; fSkipStep6 = step6; }
    void   SetRejectCentralityOutliers(Bool_t flag = kTRUE) { fRejectCentralityOutliers = flag; }
    void   SetRemoveWeakDecays(Bool_t flag = kTRUE) { fRemoveWeakDecays = flag; }
    void   SetRemoveDuplicates(Bool_t flag = kTRUE) { fRemoveDuplicates = flag; }
    void   SetSkipFastCluster(Bool_t flag = kTRUE)  { fSkipFastCluster = flag; }
    void   SetWeightPerEvent(Bool_t flag = kTRUE)   { fWeightPerEvent = flag; }
    void   SetCustomBinning(const char* binningStr) { fCustomBinning = binningStr; }
    
    AliHelperPID* GetHelperPID() { return fHelperPID; }
    void   SetHelperPID(AliHelperPID* pid){ fHelperPID = pid; }
    
  private:
    AliAnalysisTaskPhiCorrelations(const  AliAnalysisTaskPhiCorrelations &det);
    AliAnalysisTaskPhiCorrelations&   operator=(const  AliAnalysisTaskPhiCorrelations &det);
    void            AddSettingsTree();                                  // add list of settings to output list
    // Analysis methods
    void            AnalyseCorrectionMode();                            // main algorithm to get correction maps
    void            AnalyseDataMode();                                  // main algorithm to get raw distributions
    void            Initialize(); 			                // initialize some common pointer
    TObjArray* CloneAndReduceTrackList(TObjArray* tracks);
    void RemoveDuplicates(TObjArray* tracks);
    void CleanUp(TObjArray* tracks, TObject* mcObj, Int_t maxLabel);

    // General configuration
    Int_t               fDebug;           //  Debug flag
    Int_t 	        fMode;            //  fMode = 0: data-like analysis 
    				          //  fMode = 1: corrections analysis	
    Bool_t              fReduceMemoryFootprint; // reduce memory consumption by writing less debug histograms
    Bool_t		fFillMixed;		// enable event mixing (default: ON)
    Int_t  		fMixingTracks;		// size of track buffer for event mixing
    Bool_t		fCompareCentralities;	// use the z vtx axis for a centrality comparison
    Bool_t		fTwoTrackEfficiencyStudy; // two-track efficiency study on
    Float_t		fTwoTrackEfficiencyCut;   // enable two-track efficiency cut
    Int_t		fUseVtxAxis;              // use z vtx as axis (needs 7-10 times more memory!)
    Bool_t		fCourseCentralityBinning; // less centrality bins
    Bool_t		fSkipTrigger;		  // skip trigger selection
    Bool_t		fInjectedSignals;	  // check header to skip injected signals in MC
    
    // Pointers to external UE classes
    AliHelperPID*     fHelperPID;      // points to class containing common analysis algorithms
    AliAnalyseLeadingTrackUE*     fAnalyseUE;      //! points to class containing common analysis algorithms
    AliUEHistograms*  fHistos;       //! points to class to handle histograms/containers  
    AliUEHistograms*  fHistosMixed;       //! points to class to handle mixed histograms/containers  
    
    THnF* fEfficiencyCorrection;   // if non-0 this efficiency correction is applied on the fly to the filling for associated particles. The factor is multiplicative, i.e. should contain 1/efficiency. Axes: eta, pT, centrality, z-vtx
    Bool_t fCorrectTriggers;	// if true correct also trigger particles

    // Handlers and events
    AliAODEvent*             fAOD;             //! AOD Event 
    AliESDEvent*             fESD;             //! ESD Event 
    TClonesArray*            fArrayMC;         //! Array of MC particles 
    AliInputEventHandler*    fInputHandler;    //! Generic InputEventHandler 
    AliMCEvent*              fMcEvent;         //! MC event
    AliMCEventHandler*       fMcHandler;       //! MCEventHandler 
    AliEventPoolManager*     fPoolMgr;         //! event pool manager
    
    // Histogram settings
    TList*              fListOfHistos;    //  Output list of containers 
    
    // Event QA cuts
    Int_t          	fnTracksVertex;        // QA tracks pointing to principal vertex
    Double_t       	fZVertex;              // Position of Vertex in Z direction
    TString             fCentralityMethod;     // Method to determine centrality
    
    // Track cuts
    Double_t      	fTrackEtaCut;          // Eta cut on particles
    Int_t 		fOnlyOneEtaSide;       // decides that only trigger particle from one eta side are considered (0 = all; -1 = negative, 1 = positive)
    Double_t            fPtMin;                // Min pT to start correlations
    UInt_t         	fFilterBit;            // Select tracks from an specific track cut 
    UInt_t         	fSelectBit;            // Select events according to AliAnalysisTaskJetServices bit maps 
    Bool_t         	fUseChargeHadrons;     // Only use charge hadrons
    Int_t               fParticleSpeciesTrigger; // Select which particle to use for the trigger [ -1 (all, default) 0 (pions) 1 (kaons) 2 (protons) 3 (others) particles ]
    Int_t               fParticleSpeciesAssociated; // Select which particle to use for the associated [ -1 (all, default) 0 (pions) 1 (kaons) 2 (protons) 3 (others) particles ]

    Int_t fSelectCharge;           // (un)like sign selection when building correlations: 0: no selection; 1: unlike sign; 2: like sign
    Int_t fTriggerSelectCharge;    // select charge of trigger particle: 1: positive; -1 negative
    Int_t fAssociatedSelectCharge; // select charge of associated particle: 1: positive; -1 negative
    Float_t fTriggerRestrictEta;   // restrict eta range for trigger particle (default: -1 [off])
    Bool_t fEtaOrdering;           // eta ordering, see AliUEHistograms.h for documentation
    Bool_t fCutConversions;        // cut on conversions (inv mass)
    Bool_t fCutResonances;         // cut on resonances (inv mass)
    Bool_t fFillOnlyStep0; 	   // fill only step 0
    Bool_t fSkipStep6;		   // skip step 6 when filling
    Bool_t fRejectCentralityOutliers;  // enable rejection of outliers in centrality vs no track correlation
    Bool_t fRemoveWeakDecays;	   // remove secondaries from weak decays from tracks and particles
    Bool_t fRemoveDuplicates;      // remove particles with the same label (double reconstruction)
    Bool_t fSkipFastCluster;	   // skip kFastOnly flagged events (only for data)
    Bool_t fWeightPerEvent;	   // weight with the number of trigger particles per event
    TString fCustomBinning;	   // supersedes default binning if set, see AliUEHist::GetBinning or AliUEHistograms::AliUEHistograms for syntax and examples
    
    Bool_t fFillpT;                // fill sum pT instead of number density
    
    ClassDef( AliAnalysisTaskPhiCorrelations, 27); // Analysis task for delta phi correlations
  };

class AliDPhiBasicParticle : public AliVParticle
{
  public:
    AliDPhiBasicParticle(Float_t eta, Float_t phi, Float_t pt, Short_t charge)
      : fEta(eta), fPhi(phi), fpT(pt), fCharge(charge)
    {
    }
    ~AliDPhiBasicParticle() {}
    
    // kinematics
    virtual Double_t Px() const { AliFatal("Not implemented"); return 0; }
    virtual Double_t Py() const { AliFatal("Not implemented"); return 0; }
    virtual Double_t Pz() const { AliFatal("Not implemented"); return 0; }
    virtual Double_t Pt() const { return fpT; }
    virtual Double_t P() const { AliFatal("Not implemented"); return 0; }
    virtual Bool_t   PxPyPz(Double_t[3]) const { AliFatal("Not implemented"); return 0; }

    virtual Double_t Xv() const { AliFatal("Not implemented"); return 0; }
    virtual Double_t Yv() const { AliFatal("Not implemented"); return 0; }
    virtual Double_t Zv() const { AliFatal("Not implemented"); return 0; }
    virtual Bool_t   XvYvZv(Double_t[3]) const { AliFatal("Not implemented"); return 0; }

    virtual Double_t OneOverPt()  const { AliFatal("Not implemented"); return 0; }
    virtual Double_t Phi()        const { return fPhi; }
    virtual Double_t Theta()      const { AliFatal("Not implemented"); return 0; }


    virtual Double_t E()          const { AliFatal("Not implemented"); return 0; }
    virtual Double_t M()          const { AliFatal("Not implemented"); return 0; }
    
    virtual Double_t Eta()        const { return fEta; }
    virtual Double_t Y()          const { AliFatal("Not implemented"); return 0; }
    
    virtual Short_t Charge()      const { return fCharge; }
    virtual Int_t   GetLabel()    const { AliFatal("Not implemented"); return 0; }
    // PID
    virtual Int_t   PdgCode()     const { AliFatal("Not implemented"); return 0; }      
    virtual const Double_t *PID() const { AliFatal("Not implemented"); return 0; }
    
  private:
    Float_t fEta;      // eta
    Float_t fPhi;      // phi
    Float_t fpT;       // pT
    Short_t fCharge;   // charge
    
    ClassDef( AliDPhiBasicParticle, 1); // class which contains only quantities requires for this analysis to reduce memory consumption for event mixing
};
  
#endif

    
