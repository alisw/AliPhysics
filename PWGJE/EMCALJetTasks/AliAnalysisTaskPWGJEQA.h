/**
 * \file AliAnalysisTaskPWGJEQA.h
 * \brief Declaration of class AliAnalysisTaskPWGJEQA
 *
 * In this header file the class AliAnalysisTaskPWGJEQA is declared.
 *
 * \author James Mulligan <james.mulligan@yale.edu>, Yale University
 * \date Dec 1, 2016
 */

#ifndef ALIANALYSISTASKPWGJEQA_H
#define ALIANALYSISTASKPWGJEQA_H

/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

class TH1;
class TH2;
class TH3;
class THnSparse;
class AliPHOSGeometry;

#include "AliEventCuts.h"
#include "THistManager.h"
#include "AliTLorentzVector.h"
#include "AliAnalysisTaskEmcalJet.h"

/**
 * \class AliAnalysisTaskPWGJEQA
 * \brief This is a task used to do basic PWGJE QA on tracks, clusters, and jets.
 *
 * Set the names of the tracks/clusters/cells in the AddTask, as well as "mcparticles" if MC production
 * (or "" if not). Set also flags for whether to perform track/calo/jet/event QA.
 *
 * For Pt-hard productions, you should further call the function SetIsPtHard(kTRUE), and you can
 * reject outliers with SetRejectOutlierEvents(kTRUE).
 *
 * There exist post-processing scripts to efficiently plot the QA: 
 * see http://alidoc.cern.ch/AliPhysics/master/_r_e_a_d_m_ejetfw.html
 *
 * This task is based on code from Salvatore Aiola: See the tasks AliAnalysisTaskEmcalJetQA (clusters),
 * AliAnalysisTaskEmcalJetSpectraQA (jets), and AliEmcalTrackingQATask (tracks) for more detailed histograms.
 */
class AliAnalysisTaskPWGJEQA : public AliAnalysisTaskEmcalJet {
  
public:
  
  struct EventQA_t {
    EventQA_t() : fCent(0), fNTracks(0), fMaxTrack() { fNClusters[0] = 0; fNClusters[1] = 0; fNClusters[2] = 0;}
    Float_t fCent;
    Int_t fNTracks;
    Int_t fNClusters[3];
    AliTLorentzVector fMaxTrack;
    AliTLorentzVector fMaxCluster[3];
  };
  
  enum ClusterType {
    kNA       = -1,//!< Undefined
    kEMCal    = 0, //!< EMCal
    kDCal     = 1, //!< DCal
    kPHOS     = 2  //!< PHOS
  };

  AliAnalysisTaskPWGJEQA();
  AliAnalysisTaskPWGJEQA(const char *name);
  virtual ~AliAnalysisTaskPWGJEQA();

  static AliAnalysisTaskPWGJEQA* AddTaskPWGJEQA(
                                         const char* ntracks            = "usedefault",
                                         const char* nclusters          = "usedefault",
                                         const char* ncells             = "usedefault",
                                         const char *nGenLev            = "mcparticles",
                                         Bool_t      doTrackQA          = kTRUE,
                                         Bool_t      doCaloQA           = kTRUE,
                                         Bool_t      doJetQA            = kTRUE,
                                         Bool_t      doEventQA          = kTRUE,
                                         Double_t    trackPtCut         = 0.15,
                                         Double_t    clusECut           = 0.30,
                                         const char* suffix             = ""
   	   	                                       );
  
  void                        UserCreateOutputObjects();

  void                        SetUseAliEventCuts(Bool_t b)                         { fUseAliEventCuts    = b          ; }
  void                        SetUseManualEvtCuts(Bool_t input)                    { fUseManualEventCuts = input      ; }
  void                        SetCellEnergyCut(Float_t cut)                        { fCellEnergyCut      = cut        ; }
  void                        SetGeneratorLevelName(const char* name)              { fGeneratorLevelName = name       ; }
  void                        SetDetectorLevelName(const char* name)               { fDetectorLevelName  = name       ; }
  
  void                        SetDoTrackQA(Bool_t b) { fDoTrackQA = b; }
  void                        SetDoCaloQA(Bool_t b)  { fDoCaloQA  = b; }
  void                        SetDoJetQA(Bool_t b)   { fDoJetQA   = b; }
  void                        SetDoEventQA(Bool_t b) { fDoEventQA = b; }
  void                        SetRejectOutlierEvents(Bool_t b) {fRejectOutlierEvents = b; }
  void                        SetIsPtHard(Bool_t b)            {fIsPtHard = b; }
  void                        SetMaxPtBin(Bool_t b)  { fMaxPt = b; }
  
protected:
  
  void                        ExecOnce()                                                    ;
  Bool_t                      IsEventSelected()                                             ;
  Bool_t                      FillHistograms()                                              ;
  Bool_t                      RetrieveEventObjects()                                        ;
  Bool_t                      UserNotify()                                                  ;
  
  void                        AllocateTrackHistograms()                                     ;
  void                        AllocateCellHistograms()                                      ;
  void                        AllocateClusterHistograms()                                   ;
  void                        AllocateJetHistograms()                                       ;
  void                        AllocateEventQAHistograms()                                   ;
  
  void                        FillTrackHistograms()                                         ;
  void                        FillCellHistograms()                                          ;
  void                        FillClusterHistograms()                                       ;
  void                        FillJetHistograms()                                           ;
  void                        FillEventQAHistograms()                                       ;
  
  void                        GenerateHistoBins()                                           ;
  void                        AllocateDetectorLevelTHnSparse()                              ;
  void                        AllocateGeneratorLevelTHnSparse()                             ;
  void                        AllocateMatchedParticlesTHnSparse()                           ;
  void                        FillDetectorLevelTHnSparse(Double_t cent, Double_t trackEta, Double_t trackPhi, Double_t trackPt,
                                                         Double_t sigma1OverPt, Byte_t trackType);
  void                        FillGeneratorLevelTHnSparse(Double_t cent, Double_t partEta, Double_t partPhi, Double_t partPt, Byte_t findable);
  void                        FillMatchedParticlesTHnSparse(Double_t cent, Double_t partEta, Double_t partPhi, Double_t partPt,
                                                            Double_t trackEta, Double_t trackPhi, Double_t trackPt, Byte_t trackType);
  // Event selection
  Bool_t                      fUseAliEventCuts;          ///< Flag to use AliEventCuts (otherwise AliAnalysisTaskEmcal will be used)
  AliEventCuts                fEventCuts;                ///< event selection utility
  TList                      *fEventCutList;             //!<! Output list for event cut histograms
  Bool_t                      fUseManualEventCuts;       ///< Flag to use manual event cuts

  Float_t                     fCellEnergyCut;            ///< Energy cell cut
  Float_t                     fMaxPt;                    ///< Histogram pt limit
  Int_t                       fNTotClusters[3];          //!<!Total number of accepted clusters in current event (DCal/EMCal)
  AliTLorentzVector           fLeadingCluster[3];        //!<!Leading cluster in current event (EMCal/DCal)
  Int_t                       fNTotTracks;               //!<!Total number of accepted tracks in current event
  AliTLorentzVector           fLeadingTrack;             //!<!Leading track in current event
  Bool_t                      fDoTrackQA;                ///< Set whether to enable track QA
  Bool_t                      fDoCaloQA;                 ///< Set whether to enable cell/cluster QA
  Bool_t                      fDoJetQA;                  ///< Set whether to enable jet QA
  Bool_t                      fDoEventQA;                ///< Set whether to enable event QA
  TString                     fGeneratorLevelName;       ///< generator level container name
  TString                     fDetectorLevelName;        ///< detector level container name
  Bool_t                      fRejectOutlierEvents;      ///< flag to reject pythia pt-hard jet outlier events
  Bool_t                      fIsPtHard;                 ///< flag to enable pt-hard histos and make available outlier cut
  
  // Service fields (non-streamed)
  AliMCParticleContainer* fGeneratorLevel      ; //!<! generator level container
  AliTrackContainer*    fDetectorLevel         ; //!<! detector level container
  Int_t                 fNPtHistBins           ; //!<! number of pt bins
  Double_t*             fPtHistBins            ; //!<! pt bins
  Int_t                 fNEtaHistBins          ; //!<! number of eta bins
  Double_t*             fEtaHistBins           ; //!<! eta bins
  Int_t                 fNPhiHistBins          ; //!<! number of phi bins
  Double_t*             fPhiHistBins           ; //!<! phi bins
  Int_t                 fNCentHistBins         ; //!<! number of cent bins
  Double_t*             fCentHistBins          ; //!<! cent bins
  Int_t                 fNPtRelDiffHistBins    ; //!<! number of pt relative difference bins
  Double_t*             fPtRelDiffHistBins     ; //!<! pt relative difference bins
  Int_t                 fNPtResHistBins        ; //!<! number of pt res bins
  Double_t*             fPtResHistBins         ; //!<! pt res bins
  Int_t                 fNIntegerHistBins      ; //!<! number of integer bins
  Double_t*             fIntegerHistBins       ; //!<! integer bins
  AliPHOSGeometry*      fPHOSGeo               ; //!<! phos geometry
  
  THistManager          fHistManager           ; //!< Histogram manager
  
private:
  AliAnalysisTaskPWGJEQA(const AliAnalysisTaskPWGJEQA&);            // not implemented
  AliAnalysisTaskPWGJEQA &operator=(const AliAnalysisTaskPWGJEQA&); // not implemented
  
  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskPWGJEQA, 5);
  /// \endcond
};
#endif
