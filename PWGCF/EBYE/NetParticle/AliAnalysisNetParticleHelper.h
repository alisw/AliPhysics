//-*- Mode: C++ -*-

#ifndef ALIANALYSISNETPARTICLEHELPER_H
#define ALIANALYSISNETPARTICLEHELPER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
// Helper Class for for NetParticle Distributions
// Authors: Jochen Thaeder <jochen@thaeder.de>

#include "THnSparse.h"
#include "TParticle.h"
#include "TH1F.h"
#include "TF1.h"

class AliESDtrack;
class AliMCEvent;
class AliStack;
class AliPIDResponse;
class AliESDInputHandler;
class AliAODInputHandler;
class AliAODEvent;
class AliAODTrack;
class AliAODMCParticle;

class AliAnalysisNetParticleHelper : public TNamed {

 public:

  /*
   * ---------------------------------------------------------------------------------
   *                            Constructor / Destructor
   * ---------------------------------------------------------------------------------
   */

  AliAnalysisNetParticleHelper();
  virtual ~AliAnalysisNetParticleHelper();

  /*
   * ---------------------------------------------------------------------------------
   *                                    Setter
   * ---------------------------------------------------------------------------------
   */

  void SetCentralityBinMax(Int_t d)                  {fCentralityBinMax    = d;}
  void SetVertexZMax(Float_t f)                      {fVertexZMax          = f;}
  void SetRapidityMax(Float_t f)                     {fRapidityMax         = f;}
  void SetMinTrackLengthMC(Float_t f)                {fMinTrackLengthMC    = f;}
  void SetNSigmaMaxCdd(Float_t f)                    {fNSigmaMaxCdd        = f;}
  void SetNSigmaMaxCzz(Float_t f)                    {fNSigmaMaxCzz        = f;}

  void SetParticleSpecies(AliPID::EParticleType pid) {fParticleSpecies     = pid;}
  void SetControlParticleSpecies(Int_t pdgCode, Bool_t isNeutral, TString name) {
    fControlParticleCode = pdgCode;
    fControlParticleIsNeutral = isNeutral;
    fControlParticleName = name;
  }

  void SetNSigmaMaxTPC(Float_t f)                    {fNSigmaMaxTPC        = f;}
  void SetNSigmaMaxTOF(Float_t f)                    {fNSigmaMaxTOF        = f;}
  void SetMinPtForTOFRequired(Float_t f)             {fMinPtForTOFRequired = f;}

  /*
   * ---------------------------------------------------------------------------------
   *                                    Getter
   * ---------------------------------------------------------------------------------
   */
  
  AliPID::EParticleType GetParticleSpecies(){return fParticleSpecies;}

  TH1F*    GetHEventStat0()                  {return fHEventStat0;}
  TH1F*    GetHEventStat1()                  {return fHEventStat1;}
  TH1F*    GetHTriggerStat()                 {return fHTriggerStat;}
  TH1F*    GetHCentralityStat()              {return fHCentralityStat;}

  Int_t    GetCentralityBin()                {return fCentralityBin;}
  Float_t  GetCentralityPercentile()         {return fCentralityPercentile;}

  Int_t    GetControlParticleCode()          {return fControlParticleCode;}
  Bool_t   IsControlParticleNeutral()        {return fControlParticleIsNeutral;}
  TString& GetControlParticleName()          {return fControlParticleName;}

  /*
   * ---------------------------------------------------------------------------------
   *                                 Public Methods
   * ---------------------------------------------------------------------------------
   */

  /** Initialize Helper */
  Int_t Initialize(Bool_t isMC);

  /** Setup Event */
  Int_t SetupEvent(AliESDInputHandler *esdHandler, AliAODInputHandler *aodHandler, AliMCEvent *mcEvent);

  /*
   * ---------------------------------------------------------------------------------
   *                         Event / Trigger Statistics
   * ---------------------------------------------------------------------------------
   */

  /** Check if event is triggred */
  Bool_t IsEventTriggered();

  /** Fill event cut statistics */
  Bool_t IsEventRejected();

  /*
   * ---------------------------------------------------------------------------------
   *                         Accept Particle Methods - private
   * ---------------------------------------------------------------------------------
   */
  
  /** Check if charged MC particle is accepted for basic parameters */
  /** NOT possible for AODs (AliAODMCParticle NOT from TParticle)*/
  Bool_t IsParticleAcceptedBasicCharged(TParticle *particle, Int_t idxMC);
  Bool_t IsParticleAcceptedBasicCharged(AliAODMCParticle *particle);

  /** Check if neutral MC particle is accepted for basic parameters */
  /** NOT possible for AODs (AliAODMCParticle NOT from TParticle)*/
  Bool_t IsParticleAcceptedBasicNeutral(TParticle *particle, Int_t idxMC);
  Bool_t IsParticleAcceptedBasicNeutral(AliAODMCParticle *particle);
 
  /** Check if MC particle is accepted for Rapidity */
  Bool_t IsParticleAcceptedRapidity(TParticle *particle, Double_t &yP);

  /** Check if MC particle is findable tracks */
  Bool_t IsParticleFindable(Int_t label);
    
  /*
   * ---------------------------------------------------------------------------------
   *                            Accept Track Methods - public
   * ---------------------------------------------------------------------------------
   */
  
  /** Check if track is accepted for basic parameters */
  /** NOT possible with AliVTrack (GetInnerParam returns NULL) */
  Bool_t IsTrackAcceptedBasicCharged(AliESDtrack *track);
  Bool_t IsTrackAcceptedBasicCharged(AliAODTrack *track);
  
  /** Check if track is accepted for Rapidity */
  Bool_t IsTrackAcceptedRapidity(AliVTrack *track, Double_t &yP);

  /** Check if track is accepted for DCA */
  Bool_t IsTrackAcceptedDCA(AliESDtrack *track);

  /** Check if track is accepted for PID */
  Bool_t IsTrackAcceptedPID(AliVTrack *track, Double_t *pid);

  /*
   * ---------------------------------------------------------------------------------
   *                         Helper Methods
   * ---------------------------------------------------------------------------------
   */

  /** Update eta corrected TPC pid */
  void  UpdateEtaCorrectedTPCPid();

  /** Get efficiency correctionf of particle dependent on (eta, phi, pt, centrality) */
  Double_t GetTrackbyTrackCorrectionFactor(Double_t *aTrack,  Int_t flag);
  
  /** Method for the correct logarithmic binning of histograms */
  void BinLogAxis(const THnSparseF *h, Int_t axisNumber);

  ///////////////////////////////////////////////////////////////////////////////////

 private:

  AliAnalysisNetParticleHelper(const AliAnalysisNetParticleHelper&); // not implemented
  AliAnalysisNetParticleHelper& operator=(const AliAnalysisNetParticleHelper&); // not implemented

  /*
   * ---------------------------------------------------------------------------------
   *                           Initialize - Private
   * ---------------------------------------------------------------------------------
   */

  /**  Initialize event cut statistics */
  void InitializeEventStats();

  /**  Initialize trigger statistics */
  void InitializeTriggerStats();

  /**  Initialize centrality statistics */
  void InitializeCentralityStats();

  /** Initialize eta correction maps for TPC pid */
  Int_t InitializeEtaCorrection(Bool_t isMC);

  /** Initialize track by track correction matrices */
  Int_t InitializeTrackbyTrackCorrection();

  /*
   * ---------------------------------------------------------------------------------
   *                         Event / Trigger Statistics - private
   * ---------------------------------------------------------------------------------
   */
  
  /** Fill event cut statistics */
  Bool_t FillEventStats(Int_t *aEventCuts);

  /*
   * ---------------------------------------------------------------------------------
   *                             Members - private
   * ---------------------------------------------------------------------------------
   */

  AliESDInputHandler   *fESDHandler;               //! Ptr to ESD handler 
  AliPIDResponse       *fPIDResponse;              //! Ptr to PID response Object
  AliESDEvent          *fESD;                      //! Ptr to ESD event
  AliAODInputHandler   *fAODHandler;               //! Ptr to AOD handler 
  AliAODEvent          *fAOD;                      //! Ptr to AOD event
  AliMCEvent           *fMCEvent;                  //! Ptr to MC event
  AliStack             *fStack;                    //! Ptr to stack

  // =======================================================================

  Int_t                 fCentralityBin;            //  Centrality bin of current event within max centrality bin
  Float_t               fCentralityPercentile;     //  Centrality percentile of current event
  // ----------------------------------------------------------------------
  Int_t                 fCentralityBinMax;         //  Max centrality bin to be used
  Float_t               fVertexZMax;               //  VertexZ cut
  Float_t               fRapidityMax;              //  Rapidity cut
  Float_t               fMinTrackLengthMC;         //  Min track length for MC tracks
  Float_t               fNSigmaMaxCdd;             //  N Sigma for dcar / sqrt(cdd) - turn off with 0.
  Float_t               fNSigmaMaxCzz;             //  N Sigma for dcaz / sqrt(czz) - turn off with 0.
  // -----------------------------------------------------------------------
  AliPID::EParticleType fParticleSpecies;          //  Particle species on basis of AliPID
  Int_t                 fControlParticleCode;      //  PDG code control particle
  Bool_t                fControlParticleIsNeutral; //  Is control particle neutral
  TString               fControlParticleName;      //  Name of control particle
  // -----------------------------------------------------------------------
  Float_t               fNSigmaMaxTPC;             //  N Sigma for TPC PID
  Float_t               fNSigmaMaxTOF;             //  N Sigma for TOF PID
  Float_t               fMinPtForTOFRequired;      //  Min pt from where TOF is required

  // =======================================================================

  TH1F                 *fHEventStat0;              //  Event cut statistics
  TH1F                 *fHEventStat1;              //  Event cut statistics - incremental
  Int_t                 fHEventStatMax;            //  Max N cuts to be included in HEventStat
  // -----------------------------------------------------------------------
  TH1F                 *fHTriggerStat;             //  Trigger statistics
  Int_t                 fNTriggers;                //  N triggers used
  // -----------------------------------------------------------------------
  TH1F                 *fHCentralityStat;          //  Centrality statistics
  Int_t                 fNCentralityBins;          //  N centrality bins used

  // =======================================================================

  TF1                  *fEtaCorrFunc;              //! Eta correction function for TPC dE/dx  
  THnSparseF         ***fCorr0;                    // Correction matrices for particle / anti-particle
  THnSparseF         ***fCorr1;                    // Correction matrices [cross section corrected] matrices for particle/ anti-particle
  // -----------------------------------------------------------------------

  ClassDef(AliAnalysisNetParticleHelper, 1);
};

#endif
