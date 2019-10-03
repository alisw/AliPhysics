#ifndef ALIANALYSISTASKSEDSTARCHARMFRACTION_H
#define ALIANALYSISTASKSEDSTARCHARMFRACTION_H

/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//*************************************************************************
// Class AliAnalysisTaskSEDStarCharmFraction
// AliAnalysisTask for the extraction of the fraction of prompt charm for D*
// using the charm hadron impact parameter to the primary vertex
//
// Author: Jasper van der Maarel <J.vanderMaarel@uu.nl>
//*************************************************************************

class TH1D;
class AliRDHFCutsDStartoKpipi;
class AliNormalizationCounter;
class TList;
class AliAODRecoCascadeHF;
class TClonesArray;
class AliAODMCParticle;
class AliAODEvent;
class AliAODVertex;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskSEDStarCharmFraction : public AliAnalysisTaskSE {
  public:
    AliAnalysisTaskSEDStarCharmFraction();
    AliAnalysisTaskSEDStarCharmFraction(const char *name, AliRDHFCutsDStartoKpipi *cuts);

    virtual ~AliAnalysisTaskSEDStarCharmFraction();

    // Implementation of interface methods
    virtual void UserCreateOutputObjects();
    virtual void Init();
    virtual void LocalInit() {Init();}
    virtual void UserExec(Option_t *option);
    virtual void Terminate(Option_t *option);

    void SetReadMC(Bool_t readMC = kTRUE) { fReadMC = readMC; }
    Bool_t GetReadMC() { return fReadMC; }
    void SetSkipHijing(Bool_t skipHijing = kTRUE) { fSkipHijing = skipHijing; }
    Bool_t GetSkipHijing() { return fSkipHijing; }
    void SetSingleSideband(Bool_t singleSideband = kTRUE) { fSingleSideband = singleSideband; }
    Bool_t GetSingleSideband() { return fSingleSideband; }
    void SetPeakCut(Double_t *peakCut) { for (Int_t i=0;i<fNPtBins;i++) { fPeakCut[i] = peakCut[i]; } }
    Double_t *GetPeakCut() { return fPeakCut; }
    void SetSidebandCut(Double_t *sidebandCut) { for (Int_t i=0;i<fNPtBins;i++) { fSidebandCut[i] = sidebandCut[i]; } }
    Double_t *GetSidebandCut() { return fSidebandCut; }
    void SetSidebandWindow(Double_t *sidebandWindow) { for (Int_t i=0;i<fNPtBins;i++) { fSidebandWindow[i] = sidebandWindow[i]; } }
    Double_t *GetSidebandWindow() { return fSidebandWindow; }
    void SetCuts(AliRDHFCutsDStartoKpipi *cuts) { fCuts = new AliRDHFCutsDStartoKpipi(*cuts); }
    AliRDHFCutsDStartoKpipi *GetCuts() { return fCuts; }
    void SetImpParCut(Double_t impParCut) { fImpParCut = TMath::Abs(impParCut); }
    Double_t GetImpParCut() { return fImpParCut; }

  private:
    AliRDHFCutsDStartoKpipi *fCuts;                   // Cut object
    AliNormalizationCounter *fCounter;                // Counter for normalization
    Bool_t fReadMC;                                   // Flag to switch on/off the access of MC
    Bool_t fSkipHijing;                               // Flag to switch on/off the skipping of D from Hijing
    Bool_t fSingleSideband;                           // Flag to switch on/off single sideband instead of double sideband
    Double_t fImpParCut;                              // Lower value for impact parameter cut
    Double_t fPDGMDStarD0;                            // Difference between masses of D* and D0 from PDG
    Int_t fNPtBins;                                   // Number of pt bins
    TH1D *fNEvents;                                   // Event count histogram
    TTree *fTreeCandidate;                            // Tree for D* candidates
    TList *fListCandidate;                            // List for D* candidates
    TList *fListSignal;                               // List for D* signal (MC)
    TList *fListSignalPrompt;                         // List for prompt D* (MC)
    TList *fListSignalFromB;                          // List for D* from B (MC)
    TList *fListBackground;                           // List for D* background (MC)

    Double_t fPeakCut[30];                            // Max distance from delta inv mass for peak region
    Double_t fSidebandCut[30];                        // Sideband region starts at this distance from peak
    Double_t fSidebandWindow[30];                     // Size of sideband region

    Bool_t fIsSideband;                               // Candidate is in sideband
    Bool_t fIsPeak;                                   // Candidate is in peak region
    Bool_t fIsSignal;                                 // Candidate is actual D* (from MC)
    Bool_t fIsSignalPrompt;                           // Candidate is actual prompt D* (from MC)
    Bool_t fIsSignalFromB;                            // Candidate is actual D* from B (from MC)
    Bool_t fIsBackground;                             // Candidate is fake (from MC)

    AliAODVertex *fNewPrimVtx;                        // Newly determined primary vertex without D* daughters
    AliAODVertex *fDStarVtx;                          // D* vertex
    Double_t fMagneticField;                          // Magnetic field strength in current event
    Double_t fPtForTree;                              // Use this to fill pT branch
    Double_t fInvMassForTree;                         // Use this to fill invariant mass branch
    Double_t fImpParForTree;                          // Use this to fill impact parameter branch
    Double_t fTrueImpParForTree;                      // Use this to fill true impact parameter branch
    Short_t fTypeForTree;                             // Use this to fill the type branch (0=other, 1=peak, 2=SB left, 3=SB right)
    Short_t fSignalTypeForTree;                       // Use this to fill the signal type branch (data: -1, MC: 0=background, 1=signal prompt, 2=signal non-prompt)

    void SetUpList(TList *list); // Fill a TList with histograms
    void CheckInvMassDStar(AliAODRecoCascadeHF *cand); // Check if D* candidate falls within peak or sideband region
    Bool_t IsFromB(TClonesArray *arrayMC, const AliAODMCParticle *mcPartCandidate); // Check if the MC particle comes from a B
    Bool_t IsFromHijing(TClonesArray *arrayMC, const AliAODMCParticle *mcPartCandidate); // Check if the MC particle is from Hijing
    Bool_t CalculateImpactParameter(AliAODTrack *track, Double_t &d0, Double_t &d0Err); // Calculate impact parameter for a track
    Double_t CalculateImpactParameterDStar(AliAODRecoCascadeHF *cand); // Calculate impact parameter of the D*
    Double_t CalculateTrueImpactParameterDStar(AliAODMCHeader *headerMC, TClonesArray *arrayMC, AliAODRecoCascadeHF* cand); // Calculate true impact parameter of the D*
    void FillHistograms(AliAODRecoCascadeHF *cand); // Fill histograms for a D* candidate
    void FillHistogram(const char *name, Double_t value); // Fill a specific histogram in multiple lists
    void FillRegionHistogram(const char *name, Double_t value); // Fill a specific histogram in multiple lists, in the approprate regions (all, peak region, sideband region)
    void FillTrueImpactParameter(AliAODRecoCascadeHF* cand); // Fill histogram with true impact parameter distribution for D from B
    AliAODVertex *RemoveDaughtersFromPrimaryVtx(AliAODEvent *aod, AliAODRecoCascadeHF *cand); // Determine primary vertex without D* daughters
    AliAODVertex *ReconstructDStarVtx(AliAODRecoCascadeHF *cand); //Determine the D* vertex

    AliAnalysisTaskSEDStarCharmFraction(const AliAnalysisTaskSEDStarCharmFraction&); // Not implemented
    AliAnalysisTaskSEDStarCharmFraction& operator=(const AliAnalysisTaskSEDStarCharmFraction&); // Not implemented
  
    ClassDef(AliAnalysisTaskSEDStarCharmFraction, 2); // Analysis task for D* prompt charm fraction
};

#endif
