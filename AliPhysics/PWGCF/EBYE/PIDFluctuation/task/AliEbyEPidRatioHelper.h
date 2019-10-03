#ifndef ALIEBYEPIDRATIOHELPER_H
#define ALIEBYEPIDRATIOHELPER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
//=========================================================================//
//             AliEbyE Analysis for Particle Ratio Fluctuation             //
//                   Deepika Rathee  | Satyajit Jena                       //
//                   drathee@cern.ch | sjena@cern.ch                       //
//                  Date: Wed Jul  9 18:38:30 CEST 2014                    // 
//          New approch to find particle ratio to reduce memory            //
//                             (Test Only)                                 //
//        Copied from NetParticle Classes
//        Origin: Authors: Jochen Thaeder <jochen@thaeder.de>
//                         Michael Weber <m.weber@cern.ch>
//=========================================================================//

#include "THnBase.h"
#include "THn.h"
#include "TH1F.h"
#include "TF1.h"
#include "TProfile2D.h"
#include "TRandom3.h"

class AliESDtrack;
class AliMCEvent;
class AliStack;
class AliPIDResponse;
class AliESDtrackCuts;
class AliInputEventHandler;
class AliESDInputHandler;
class AliAODInputHandler;
class AliAODEvent;
class AliAODTrack;
class AliAODMCParticle;
class AliMCParticle;

class AliEbyEPidRatioHelper : public TNamed {

 public:

  AliEbyEPidRatioHelper();
  virtual ~AliEbyEPidRatioHelper();

  void SetCentralityBinMax(Int_t d)            {fCentralityBinMax    = d;}
  void SetVertexZMax(Float_t f)                {fVertexZMax          = f;}
  void SetRapidityMax(Float_t f)               {fRapidityMax         = f;}
  void SetMinTrackLengthMC(Float_t f)          {fMinTrackLengthMC    = f;}
  void SetNSigmaMaxCdd(Float_t f)              {fNSigmaMaxCdd        = f;}
  void SetNSigmaMaxCzz(Float_t f)              {fNSigmaMaxCzz        = f;}
  void SetPIDStrategy(Int_t i)                 {fPIDStrategy         = i;}
  void SetNSigmaMaxITS(Float_t f)              {fNSigmaMaxITS        = f;}
  void SetNSigmaMaxTPC(Float_t f)              {fNSigmaMaxTPC        = f;}
  void SetNSigmaMaxTPClow(Float_t f)           {fNSigmaMaxTPClow     = f;}
  void SetNSigmaMaxTOF(Float_t f)              {fNSigmaMaxTOF        = f;}
  void SetMinPtForTOFRequired(Float_t f)       {fMinPtForTOFRequired = f;}
  void SetMaxPtForTPClow(Float_t f)            {fMaxPtForTPClow      = f;}
  void SetNSubSamples(Int_t i)                 {fSubSamples         = i;}
  
  TH1F*    GetHEventStat0()                    {return fHEventStat0;     }
  TH1F*    GetHEventStat1()                    {return fHEventStat1;     }
  TH1F*    GetHTriggerStat()                   {return fHTriggerStat;    }
  TH1F*    GetHCentralityStat()                {return fHCentralityStat; }
  TH1F*    GetHCentralityPercentile()          {return fHCentralityPer; }
  TH1F*    GetHCentralityPercentileAll()       {return fHCentralityPerAll; }
  Int_t    GetCentralityBin()                  {return fCentralityBin;   }
  Float_t  GetMaxPtForTPClow()                 {return fMaxPtForTPClow;  }
  Float_t  GetRapidityMax()                    {return fRapidityMax;     }
  Float_t  GetPhiMin()                         {return fPhiMin;          }
  Float_t  GetPhiMax()                         {return fPhiMax;          }
  AliESDtrackCuts* GetESDTrackCuts()           {return fESDTrackCuts;    }
  Bool_t           GetIsMC()                   {return fIsMC;            }
  Bool_t           GetIsRatio()                {return fIsRatio;         }
  Bool_t           GetIsPtBin()                {return fIsPtBin;         }
  Bool_t           GetDetWise()                {return fIsDetectorWise ; }


  Bool_t   GetUsePID(Int_t i)                         { if (i == 0) return kFALSE; else return kTRUE;}

  AliPID::EParticleType GetParticleSpecies(Int_t i) {
    if (i == 0) return AliPID::kUnknown;
    else if (i == 1) return  AliPID::kPion;  
    else if (i == 2) return  AliPID::kKaon;  
    else if (i == 3) return  AliPID::kProton;  
    else return AliPID::kUnknown;
  }
  
  Int_t GetPdg(Int_t i) {
    if (i == 0) return 0;
    else if (i == 1) return  211;  
    else if (i == 2) return  321;  
    else if (i == 3) return  2212;  
    else return 0;
  }



  Int_t            GetAODtrackCutBit()         {return fAODtrackCutBit;  }
  AliMCEvent*           GetMCEvent()           {return fMCEvent;         }
  AliInputEventHandler* GetInputEventHandler() {return fInputEventHandler;}
  void SetPhiRange(Float_t f1, Float_t f2);
  Float_t  GetCentralityPercentile()           {return fCentralityPercentile;}
  Float_t  GetMinPtForTOFRequired()            {return fMinPtForTOFRequired;}
  Int_t    GetSubSampleIdx()                   {return fSubSampleIdx;}
  Int_t    GetNSubSamples()                    {return fSubSamples;}
  /** Initialize Helper */
  Int_t Initialize(AliESDtrackCuts *cuts, Bool_t isMC, Bool_t isRatio, Bool_t isPtBin, Bool_t isDetWise, Int_t trackCutBit, Int_t modeDistCreation);

  /** Setup Event */
  Int_t SetupEvent(AliESDInputHandler *esdHandler, AliAODInputHandler *aodHandler, AliMCEvent *mcEvent);

  /** Check if event is triggred */
  Bool_t IsEventTriggered();

  /** Fill event cut statistics */
  Bool_t IsEventRejected();

  /** Check if charged MC particle is accepted for basic parameters */
  Bool_t IsParticleAcceptedBasicCharged(AliVParticle *particle, Int_t idxMC);

  /** Check if neutral MC particle is accepted for basic parameters */
  Bool_t IsParticleAcceptedBasicNeutral(AliVParticle *particle, Int_t idxMC);
 
  /** Check if MC particle is accepted for Rapidity */
  Bool_t IsParticleAcceptedRapidity(AliVParticle *particle, Double_t &yP, Int_t gCurPid);

  /** Check if MC particle is accepted for Phi */
  Bool_t IsParticleAcceptedPhi(AliVParticle *particle);

  /** Check if MC particle is findable tracks */
  Bool_t IsParticleFindable(Int_t label);
  
  /** Check if track is accepted for basic parameters */
  Bool_t IsTrackAcceptedBasicCharged(AliVTrack *track);
  
  /** Check if track is accepted for Rapidity */
  Bool_t IsTrackAcceptedRapidity(AliVTrack *track, Double_t &yP, Int_t gCurPid);

  /** Check if track is accepted for DCA */
  Bool_t IsTrackAcceptedDCA(AliVTrack *track);

  /** Check if track is accepted for PID */
  Bool_t IsTrackAcceptedPID(AliVTrack *track, Double_t *pid, AliPID::EParticleType gCurPid);

  /** Check if trackis  accepted for Phi */
  Bool_t IsTrackAcceptedPhi(AliVTrack *track);

  /** Method for the correct logarithmic binning of histograms 
   *  and Update MinPtForTOFRequired, using the pT log-scale 
   */
  void BinLogAxis(const THnBase *h, Int_t axisNumber, AliESDtrackCuts* cuts = NULL);

  // void SetIsRatio(Bool_t b)                    {fIsRatio             = b;}  
  // void SetIsPtBin(Bool_t   b)                  {fIsPtBin             = b;}  
  
  static const Float_t fgkfHistBinWitdthRap;   // Histogram std bin width for rapidity/eta
  static const Float_t fgkfHistBinWitdthPt;    // Histogram std bin width for pt

  static const Float_t fgkfHistRangeCent[];    // Histogram range for centrality
  static const Int_t   fgkfHistNBinsCent;      // Histogram N bins for centrality
  static const Float_t fgkfHistRangeEta[];     // Histogram range for eta
  static const Int_t   fgkfHistNBinsEta;       // Histogram N bins for eta
  static const Float_t fgkfHistRangeRap[];     // Histogram range for rapidity
  static const Int_t   fgkfHistNBinsRap;       // Histogram N bins for rapidity
  static const Float_t fgkfHistRangePhi[];     // Histogram range for phi
  static const Int_t   fgkfHistNBinsPhi;       // Histogram N bins for phi
  static const Float_t fgkfHistRangePt[];      // Histogram range for pt
  static const Int_t   fgkfHistNBinsPt;        // Histogram N bins for pt
  static const Float_t fgkfHistRangeSign[];    // Histogram range for sign
  static const Int_t   fgkfHistNBinsSign;      // Histogram N bins for sign

  static const Char_t* fgkEventNames[];         // Event names 
  static const Char_t* fgkCentralityMaxNames[]; // Centrality names 
  static const Char_t* fgkTriggerNames[];       // Trigger names 
  static const Char_t* fgkCentralityNames[];    // Centrality names  
  static const Char_t* fgkPidName[4]; //
  static const Char_t* fgkPidLatex[4][2];//
  static const Char_t* fgkPidShLatex[4];//
  static const Char_t* fgkPidTitles[4][2];//
 

  static const Char_t *fgkNetHistName[4]; //
  static const Char_t *fgkNetHistLatex[4]; //
  static const Int_t fgkfNetHistBin[4][4]; //




 private:

  AliEbyEPidRatioHelper(const AliEbyEPidRatioHelper&); // not implemented
  AliEbyEPidRatioHelper& operator=(const AliEbyEPidRatioHelper&); // not implemented

  void InitializeEventStats();
  void InitializeTriggerStats();
  void InitializeCentralityStats();
  Bool_t FillEventStats(Int_t *aEventCuts);

  Int_t                 fModeDistCreation;         //  Dist creation mode       : 1 = on | 0 = off 
  AliInputEventHandler *fInputEventHandler;        //! Ptr to input event handler (ESD or AOD)
  AliPIDResponse       *fPIDResponse;              //! Ptr to PID response Object
  AliESDEvent          *fESD;                      //! Ptr to ESD event
  AliESDtrackCuts      *fESDTrackCuts;             //! Ptr to ESD cuts  
  AliAODEvent          *fAOD;                      //! Ptr to AOD event
  Int_t                 fAODtrackCutBit;           //  Track filter bit for AOD tracks
  Bool_t                fIsMC;                     //  Is MC event
  AliMCEvent           *fMCEvent;                  //! Ptr to MC event
  AliStack             *fStack;                    //! Ptr to stack
  Int_t                 fCentralityBin;            //  Centrality bin of current event within max centrality bin
  Float_t               fCentralityPercentile;     //  Centrality percentile of current event
  Int_t                 fCentralityBinMax;         //  Max centrality bin to be used
  Float_t               fVertexZMax;               //  VertexZ cut
  Float_t               fRapidityMax;              //  Rapidity cut
  Float_t               fPhiMin;                   //  Phi min cut
  Float_t               fPhiMax;                   //  Phi max cut
  Float_t               fMinTrackLengthMC;         //  Min track length for MC tracks
  Float_t               fNSigmaMaxCdd;             //  N Sigma for dcar / sqrt(cdd) - turn off with 0.
  Float_t               fNSigmaMaxCzz;             //  N Sigma for dcaz / sqrt(czz) - turn off with 0.
 
  Int_t                 fPIDStrategy;              //  PID Strategy to be used
  Float_t               fNSigmaMaxITS;             //  N Sigma for ITS PID
  Float_t               fNSigmaMaxTPC;             //  N Sigma for TPC PID
  Float_t               fNSigmaMaxTPClow;          //  N Sigma for TPC PID lower part
  Float_t               fNSigmaMaxTOF;             //  N Sigma for TOF PID
  Float_t               fMinPtForTOFRequired;      //  Min pt from where TOF is required
  Float_t               fMaxPtForTPClow;           //  Max pt until TPClow is used
  TH1F                 *fHEventStat0;              //  Event cut statistics
  TH1F                 *fHEventStat1;              //  Event cut statistics - incremental
  Int_t                 fHEventStatMax;            //  Max N cuts to be included in HEventStat
  TH1F                 *fHTriggerStat;             //  Trigger statistics
  Int_t                 fNTriggers;                //  N triggers used
  TH1F                 *fHCentralityStat;          //  Centrality statistics
  TH1F                 *fHCentralityPer;           //  Centrality Percentile
  TH1F                 *fHCentralityPerAll;        //  Centrality Percentile
  Int_t                 fNCentralityBins;          //  N centrality bins used

  Int_t                 fSubSamples;              //  N subsamples
  TRandom3             *fRandom;                   //  Random generator
  Int_t                 fSubSampleIdx;             //  Subsample idx for current event
  
  Bool_t                fIsRatio;                  //  Is ratio
  Bool_t                fIsPtBin;                  //  Is Pt Bin
  Bool_t                fIsDetectorWise;           // is detector wise

  ClassDef(AliEbyEPidRatioHelper, 1);
};

#endif
