/**************************************************************
 *                                                            *
 * class used for the extraction of J/psi-hadron correlations *
 *                                                            *
 * authors: Lucas Altenkamper (lucas.altenkamper@cern.ch)     *
 *          Antoine Lardeux   (antoine.lardeux@cern.ch)       *
 *                                                            *
 * 12/08/2018                                                 *
 *                                                            *
 **************************************************************/

#ifndef ALICORRELATIONEXTRACTION_H
#define ALICORRELATIONEXTRACTION_H

#include <TObject.h>
#include <THn.h>
#include <THnSparse.h>
#include <TF1.h>

#include "AliResonanceFits.h"

class TMinuit;
class TH1;

//_______________________________________________________________________________
class AliCorrelationExtraction : public TObject {
  
  public:

    AliCorrelationExtraction();
    virtual ~AliCorrelationExtraction();
  
    enum BackgroundMethods {
      kBkgNone=-1,
      kBkgFitting,
      kBkgSideband,
      kBkgLikeSign,
      kBkgInterpolation,
      kBkgSuperposition,
      kBkgSuperpositionTwoComponent,
      kNBackgroundMethods=6
    };
    enum TriggerValues {
      kSig=0,
      kSigErr,
      kBkg,
      kBkgErr,
      kSplusB,
      kSplusBErr,
      kNTrigSEPP,         // no. trigger SE-PP (N_++)
      kNTrigSEPPErr,
      kNTrigSEMM,         // no. trigger SE-MM (N_--)
      kNTrigSEMMErr,
      kNTrigMEOS,         // no. trigger ME-OS (M_+-)
      kNTrigMEOSErr,
      kNTrigMEPP,         // no. trigger ME-PP (M_++)
      kNTrigMEPPErr,
      kNTrigMEMM,         // no. trigger ME-MM (M_--)
      kNTrigMEMMErr,
      kR,                 // M_+-/2*sqrt(M_++*M_--)
      kRErr,
      kBkgComb,           // 2*R*sqrt(N_++*N_--)
      kBkgCombErr,
      kSigFrac,           // S/S+B
      kSigFracErr,
      kBkgFrac,           // B/S+B
      kBkgFracErr,
      kBkgCombFrac,       // B_comb/S+B
      kBkgCombFracErr,
      kNTriggerValues=26
    };

    // maximum array dimensions
    static const Int_t kNMaxVariables             = 10;   // variables (i.e. THnF dimensions)
    static const Int_t kNMaxMixingVariables       = 10;   // mixing variables (i.e. THnF dimensions)
    static const Int_t kNMaxBackgroundMassRanges  = 10;   // background mass windows
    static const Int_t kNMaxMixingVarBins         = 100;  // mixing variable bins in THnF
    static const Int_t kNMaxDeltaPhiBins          = 100;  // delta phi bins in THnF
    static const Int_t kNMaxDeltaEtaBins          = 100;  // delta eta bins in THnF
  
    // setters
    void SetHistograms(THnF* seos, THnF* sepp=0x0, THnF* semm=0x0, THnF* meos=0x0, THnF* mepp=0x0, THnF* memm=0x0);
    void SetHistograms(THnSparseF* seos, THnSparseF* sepp=0x0, THnSparseF* semm=0x0, THnSparseF* meos=0x0, THnSparseF* mepp=0x0, THnSparseF* memm=0x0);
    void SetPairHistograms(THnF* sepp=0x0, THnF* semm=0x0, THnF* meos=0x0, THnF* mepp=0x0, THnF* memm=0x0);
    void SetSEOSHistogram(THnF* h) {fSEOS = h; fProcessDone = kFALSE;}
    void SetSEOSHistogram(THnSparseF* h) {fSEOSSparse = h; fProcessDone = kFALSE;}
    void SetMEOSHistogram(THnF* h) {fMEOS = h; fProcessDone = kFALSE;}
    void SetMEOSHistogram(THnSparseF* h) {fMEOSSparse = h; fProcessDone = kFALSE;}
    void SetMEOS2Histogram(THnF* h) {fMEOS2 = h; fProcessDone = kFALSE;}
    void SetMEOS2Histogram(THnSparseF* h) {fMEOS2Sparse = h; fProcessDone = kFALSE;}
    void SetSELSHistogram(THnF* hpp, THnF* hmm) {fSEPP = hpp; fSEMM = hmm; fProcessDone = kFALSE;}
    void SetSELSHistogram(THnSparseF* hpp, THnSparseF* hmm) {fSEPPSparse = hpp; fSEMMSparse = hmm; fProcessDone = kFALSE;}
    void SetMELSHistogram(THnF* hpp, THnF* hmm) {fMEPP = hpp; fMEMM = hmm; fProcessDone = kFALSE;}
    void SetMELSHistogram(THnSparseF* hpp, THnSparseF* hmm) {fMEPPSparse = hpp; fMEMMSparse = hmm; fProcessDone = kFALSE;}
    void SetMEOSPairHistogram(THnF* h) {fMEOSPair = h; fProcessDone = kFALSE;}
    void SetSELSPairHistogram(THnF* hpp, THnF* hmm) {fSEPPPair = hpp; fSEMMPair = hmm; fProcessDone = kFALSE;}
    void SetMELSPairHistogram(THnF* hpp, THnF* hmm) {fMEPPPair = hpp; fMEMMPair = hmm; fProcessDone = kFALSE;}
    void SetAliResonanceFitsObject(AliResonanceFits* resonanceFits) {fResonanceFits = (AliResonanceFits*)resonanceFits->Clone("ResonanceFits"); fProcessDone = kFALSE;}
    void SetBackgroundMethod(Int_t method) {fOptionBkgMethod = method; fProcessDone = kFALSE;}
    void SetBackgroundFitFunction(TF1* fitFunc) {fBkgFitFunction = (TF1*)fitFunc->Clone("BkgFitFunction"); fProcessDone = kFALSE;}
    void SetSignalMassWindow(Double_t min, Double_t max) {fMassSignalRange[0] = min; fMassSignalRange[1] = max; fProcessDone = kFALSE;}
    void SetBackgroundMassWindows(Int_t n, Double_t* min, Double_t* max);
    void SetVerbose() {fVerboseFlag = kTRUE;}
  
    // add variables and set ranges on the THnF
    void AddVariables(Int_t nVars, Int_t* vars, Int_t* indices);
    void AddVariable(Int_t var, Int_t index);
    void SetVarRange(Int_t var, Double_t* lims);
    void SetVarRange(Int_t var, Double_t min, Double_t max);
    void AddMixingVariables(Int_t nVars, Int_t* vars, Int_t* indices);
    void AddMixingVariable(Int_t var, Int_t index);

    // indicate mass, delta phi and delta eta variables
    void SetMassVariable(Int_t var) {fMassVariable = var; fProcessDone = kFALSE;}
    void SetDeltaPhiVariable(Int_t var) {fDeltaPhiVariable = var; fProcessDone = kFALSE;}
    void SetDeltaEtaVariable(Int_t var) {fDeltaEtaVariable = var; fProcessDone = kFALSE;}
    void SetMassVariablePairIndex(Int_t index) {fMassVariableIndexPair = index; fProcessDone = kFALSE;}

    // getters
    TH1D*             GetInclusiveCF1D() const {return (fProcessDone ? fInclusiveCF1D : 0x0);}
    TH1D*             GetInclusiveCF1D(Int_t massWindow) const {return (fProcessDone ? fInclusiveCF1DBackgroundMassWindow[massWindow] : 0x0);}
    TH1D*             GetInclusiveCF1D(Int_t phiBin, Int_t etaBin) const {return (fProcessDone ? fInclusiveCF1DInvMass[phiBin][etaBin] : 0x0);}
    TH1D*             GetInclusiveCF1DBackground(Int_t phiBin, Int_t etaBin) const {return (fProcessDone ? fInclusiveCF1DInvMassBackground[phiBin][etaBin] : 0x0);}
    TH2D*             GetInclusiveCF2D() const {return (fProcessDone ? fInclusiveCF2D : 0x0);}
    TH2D*             GetInclusiveCF2D(Int_t massWindow) const {return (fProcessDone ? fInclusiveCF2DBackgroundMassWindow[massWindow] : 0x0);}
    TH3D*             GetInclusiveCF3D() const {return (fProcessDone ? fInclusiveCF3D : 0x0);}
    TH1D*             GetBackgroundCF1D() const {return (fProcessDone ? fBackgroundCF1D : 0x0);}
    TH1D*             GetBackgroundCF1D(Int_t phiBin, Int_t etaBin) const {return (fProcessDone ? fBackgroundCF1DInvMass[phiBin][etaBin] : 0x0);}
    TH2D*             GetBackgroundCF2D() const {return (fProcessDone ? fBackgroundCF2D : 0x0);}
    TH1D*             GetCombinatorialBackgroundCF1D(Int_t massWindow) const {return (fProcessDone ? fCombinatorialBackgroundCF1D[massWindow] : 0x0);}
    TH2D*             GetCombinatorialBackgroundCF2D(Int_t massWindow) const {return (fProcessDone ? fCombinatorialBackgroundCF2D[massWindow] : 0x0);}
    TH1D*             GetSignalCF1D() const {return (fProcessDone ? fSignalCF1D : 0x0);}
    TH1D*             GetSignalCF1D(Int_t phiBin, Int_t etaBin) const {return (fProcessDone ? fSignalCF1DInvMass[phiBin][etaBin] : 0x0);}
    TH2D*             GetSignalCF2D() const {return (fProcessDone ? fSignalCF2D : 0x0);}
    AliResonanceFits* GetAliResonanceFitsObject() const {return (fProcessDone ? fResonanceFits : 0x0);}
    Int_t             GetBackgroundMethod() const {return (fProcessDone ? fOptionBkgMethod : 0x0);}
    TF1*              GetBackgroundFitFunction() const {return (fProcessDone ? fBkgFitFunction : 0x0);}
    TF1*              GetBackgroundFitFunction(Int_t phiBin, Int_t etaBin) const {return (fProcessDone ? fBackgroundCF1DInvMassFit[phiBin][etaBin] : 0x0);}
    TH1D*             GetPairInvMassSEPP() const {return (fProcessDone ? fSEPPPairInvMass : 0x0);}
    TH1D*             GetPairInvMassSEMM() const {return (fProcessDone ? fSEMMPairInvMass : 0x0);}
    TH1D*             GetPairInvMassMEOS() const {return (fProcessDone ? fMEOSPairInvMass : 0x0);}
    TH1D*             GetPairInvMassMEPP() const {return (fProcessDone ? fMEPPPairInvMass : 0x0);}
    TH1D*             GetPairInvMassMEMM() const {return (fProcessDone ? fMEMMPairInvMass : 0x0);}
    TH2D*             GetSEOS() const {return (fProcessDone ? fSEOSNorm : 0x0);}
    TH2D*             GetSEOS(Int_t massWindow) const {return (fProcessDone ? fSEOSNormBackgroundMassWindow[massWindow] : 0x0);}
    TH2D*             GetSEPP(Int_t massWindow) const {return (fProcessDone ? fSEPPNormBackgroundMassWindow[massWindow] : 0x0);}
    TH2D*             GetSEMM(Int_t massWindow) const {return (fProcessDone ? fSEMMNormBackgroundMassWindow[massWindow] : 0x0);}
    TH2D*             GetMEOS() const {return (fProcessDone ? fMEOSNorm : 0x0);}
    TH2D*             GetMEOS(Int_t massWindow) const {return (fProcessDone ? fMEOSNormBackgroundMassWindow[massWindow] : 0x0);}
    TH2D*             GetMEOS2(Int_t massWindow) const {return (fProcessDone ? fMEOS2NormBackgroundMassWindow[massWindow] : 0x0);}
    TH2D*             GetMEPP(Int_t massWindow) const {return (fProcessDone ? fMEPPNormBackgroundMassWindow[massWindow] : 0x0);}
    TH2D*             GetMEMM(Int_t massWindow) const {return (fProcessDone ? fMEMMNormBackgroundMassWindow[massWindow] : 0x0);}
    const Double_t*   GetTriggerValuesSignal() const {return (fProcessDone ? fTrigValSig : 0x0);}
    const Double_t*   GetTriggerValuesBackground(Int_t massWindow) const {return (fProcessDone ? fTrigValBkg[massWindow] : 0x0);}

    // processing
    Bool_t  Process();
    void    PrintUserOptions();  // print summary of user options
  
  private:
  
    // input histograms
    THnF*       fSEOS;
    THnSparseF* fSEOSSparse;
    THnF*       fSEPP;
    THnSparseF* fSEPPSparse;
    THnF*       fSEMM;
    THnSparseF* fSEMMSparse;
    THnF*       fMEOS;
    THnSparseF* fMEOSSparse;
    THnF*       fMEOS2;       // mixed event histogram for 2-component superposition method (NOTE: histogram name not ideal)
                              // NOTE: needs 1st dielectron leg and hadron from one event, 2nd dielectron leg from different event
    THnSparseF* fMEOS2Sparse;
    THnF*       fMEPP;
    THnSparseF* fMEPPSparse;
    THnF*       fMEMM;
    THnSparseF* fMEMMSparse;
    THnF*       fSEPPPair;
    THnF*       fSEMMPair;
    THnF*       fMEOSPair;
    THnF*       fMEPPPair;
    THnF*       fMEMMPair;
  
    // output histograms
    TH2D* fSEOSNorm;
    TH2D* fSEOSNormBackgroundMassWindow[kNMaxBackgroundMassRanges];
    TH2D* fSEPPNormBackgroundMassWindow[kNMaxBackgroundMassRanges];
    TH2D* fSEMMNormBackgroundMassWindow[kNMaxBackgroundMassRanges];
    TH1D* fSEPPPairInvMass;
    TH1D* fSEMMPairInvMass;
    TH1D* fMEOSPairInvMass;
    TH1D* fMEPPPairInvMass;
    TH1D* fMEMMPairInvMass;
    TH2D* fMEOSNorm;
    TH2D* fMEOSNormBackgroundMassWindow[kNMaxBackgroundMassRanges];
    TH2D* fMEOS2NormBackgroundMassWindow[kNMaxBackgroundMassRanges];
    TH2D* fMEPPNormBackgroundMassWindow[kNMaxBackgroundMassRanges];
    TH2D* fMEMMNormBackgroundMassWindow[kNMaxBackgroundMassRanges];
    TH1D* fInclusiveCF1D;
    TH1D* fInclusiveCF1DInvMass[kNMaxDeltaPhiBins][kNMaxDeltaEtaBins];
    TH1D* fInclusiveCF1DInvMassBackground[kNMaxDeltaPhiBins][kNMaxDeltaEtaBins];
    TH1D* fInclusiveCF1DBackgroundMassWindow[kNMaxBackgroundMassRanges];
    TH2D* fInclusiveCF2D;
    TH2D* fInclusiveCF2DBackgroundMassWindow[kNMaxBackgroundMassRanges];
    TH3D* fInclusiveCF3D;
    TH1D* fBackgroundCF1D;
    TF1*  fBackgroundCF1DInvMassFit[kNMaxDeltaPhiBins][kNMaxDeltaEtaBins];
    TH1D* fBackgroundCF1DInvMass[kNMaxDeltaPhiBins][kNMaxDeltaEtaBins];
    TH2D* fBackgroundCF2D;
    TH1D* fCombinatorialBackgroundCF1D[kNMaxBackgroundMassRanges];
    TH2D* fCombinatorialBackgroundCF2D[kNMaxBackgroundMassRanges];
    TH1D* fSignalCF1D;
    TH1D* fSignalCF1DInvMass[kNMaxDeltaPhiBins][kNMaxDeltaEtaBins];
    TH2D* fSignalCF2D;
  
    // variables
    Int_t     fNVariables;                                    // number of variables to be handled
    Int_t     fVariables[kNMaxVariables];                     // list of variables
    Double_t  fVarLimits[kNMaxVariables][2];                  // variable limits (i.e. "cuts" or selection ranges)
    Int_t     fVarIndices[kNMaxVariables];                    // indices of variables in THnF
    Int_t     fNMixingVariables;                              // number of variables used in event mixing
    Int_t     fMixingVariables[kNMaxMixingVariables];         // list of mixing variables
    Int_t     fMixingVarIndices[kNMaxMixingVariables];        // indices of mixing variables in THnF
    Int_t     fNMixingVarBins[kNMaxMixingVariables];          // number of bins per mixing variable
    Double_t  fMixingVarBinLimits[kNMaxMixingVariables][kNMaxMixingVarBins][2]; // mixing variable bin limits
    Int_t     fMassVariable;
    Int_t     fMassVariableIndex;
    Int_t     fMassVariableIndexPair;                         // mass variable index in LS pair histogram: fSEPPPair and fSEMMPair
    Int_t     fDeltaPhiVariable;
    Int_t     fDeltaPhiVariableIndex;
    Int_t     fDeltaEtaVariable;
    Int_t     fDeltaEtaVariableIndex;

    // user options
    Bool_t            fVerboseFlag;
    Bool_t            fUseMixingVars;
    AliResonanceFits* fResonanceFits;
    Int_t             fOptionBkgMethod;
    TF1*              fBkgFitFunction;
  
    // mass ranges
    Int_t     fNBackgroundMassRanges;                               // number of background mass windows
    Double_t  fBackgroundMassRanges[kNMaxBackgroundMassRanges][2];  // background mass windows
    Double_t  fMassSignalRange[2];                                  // signal mass window
  
    // values
    Double_t fTrigValSig[kNTriggerValues];                            // trigger values signal mass range
    Double_t fTrigValBkg[kNMaxBackgroundMassRanges][kNTriggerValues]; // trigger values background mass ranges
  
    // process flag
    Bool_t fProcessDone;
  
    // member functions
    void    ApplyUserRanges(THnBase* h);
    Bool_t  Initialize();
    Bool_t  NormalizeToNearSidePeak(TH2D* h);
    Bool_t  InBackgroundRange(Double_t min, Double_t max, Int_t& index);
    TH1D*   ProjectToDeltaPhi(TH2D* hIn, TString name);
    Bool_t  CalculateInclusiveCorrelationInMixingBins(Int_t currentVar, Int_t& nCalls,
                                                      THnBase* seos, THnBase* meos, TH2D* (&inclCF));
    Bool_t  CalculateInclusiveCorrelation(Double_t minMass, Double_t maxMass,
                                          TH2D* (&seos), TH2D* (&meos),
                                          TH2D* (&incl2D), TH1D* (&incl1D));
    Bool_t  CalculateBackgroundCorrelationFitting();
    Bool_t  CalculateBackgroundCorrelationSideband();
    Bool_t  CalculateBackgroundCorrelationLikeSign();
    Bool_t  CalculateBackgroundCorrelationInterpolation();
    Bool_t  CalculateBackgroundCorrelationSuperposition();
    Bool_t  CalculateBackgroundCorrelationSuperpositionTwoComponent();
    Bool_t  CalculateSignalCorrelation();
  
  ClassDef(AliCorrelationExtraction, 0);
};

#endif
