/**
 * \file AliAnalysisTaskEmcalJetQA.h
 * \brief Declaration of class AliAnalysisTaskEmcalJetQA
 *
 * In this header file the class AliAnalysisTaskEmcalJetQA is declared.
 * This is a task used to do basic QA on track and clusters for EMCal jet analysis.
 *
 * \author Salvatore Aiola <salvatore.aiola@cern.ch>, Yale University
 * \date Feb 1, 2016
 */

#ifndef ALIANALYSISTASKEMCALJETQA_H
#define ALIANALYSISTASKEMCALJETQA_H

/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

class TH1;
class TH2;
class TH3;
class THnSparse;
class AliVVZERO;

#include "THistManager.h"
#include "AliTLorentzVector.h"
#include "AliAnalysisTaskEmcalLight.h"

/**
 * \class AliAnalysisTaskEmcalJetQA
 * \brief Implementation of a task to perform basic QA on tracks and clusters
 *
 * Implementation of a task that performs basic QA on tracks and clusters
 * for EMCal jet analysis.
 */
class AliAnalysisTaskEmcalJetQA : public AliAnalysisTaskEmcalLight {

public:

  struct EventQA_t {
    EventQA_t() : fCent(0), fCent2(0), fCent3(0), fV0A(0), fV0C(0), fEP(0), fNTracks(0), fNCells(0), fMaxTrack() { fNClusters[0] = 0; fNClusters[1] = 0; }

    Float_t fCent;
    Float_t fCent2;
    Float_t fCent3;
    Float_t fV0A;
    Float_t fV0C;
    Float_t fEP;

    Int_t fNTracks;
    Int_t fNClusters[2];
    Int_t fNCells;

    AliTLorentzVector fMaxTrack;
    AliTLorentzVector fMaxCluster[2];
  };

  AliAnalysisTaskEmcalJetQA();
  AliAnalysisTaskEmcalJetQA(const char *name);
  virtual ~AliAnalysisTaskEmcalJetQA();

  void                        UserCreateOutputObjects();

  void                        SetCellEnergyCut(Float_t cut)                        { fCellEnergyCut      = cut        ; }
  void                        SetParticleLevel(Bool_t s)                           { fParticleLevel      = s          ; }
  void                        SetMC(Bool_t m)                                      { fIsMC               = m          ; }
  void                        SetAdditionalCentEst(const char* meth2, const char* meth3="") { fCentMethod2 = meth2; fCentMethod3 = meth3; }
  void                        SetDoV0QA(Int_t b)                                   { fDoV0QA             = b          ; }
  void                        SetDoEPQA(Int_t b)                                   { fDoEPQA             = b          ; }
  void                        SetMaxCellsInCluster(Int_t b)                        { fMaxCellsInCluster  = b          ; }
  void                        SetDoLeadingObjectPosition(Int_t b)                  { fDoLeadingObjectPosition = b     ; }
  void                        SetSeparateEMCalDCal(Bool_t b)                       { fSeparateEMCalDCal = b           ; }
  void                        SetPtBin(Float_t w, Float_t max)                     { fPtBinWidth        = w; fMaxPt = max ; }
  void                        SetIsEmbedded(Bool_t i)                              { fIsEmbedded        = i           ; }

protected:

  void                        ExecOnce()                                                    ;
  Bool_t                      FillHistograms()                                              ;
  void                        FillEventQAHisto(const EventQA_t& eventQA);
  Bool_t                      RetrieveEventObjects()                                        ;
  Int_t                       DoCellLoop()                                      ;
  void                        DoTrackLoop()                                                 ;
  void                        DoClusterLoop()                                               ;
  Double_t                    GetFcross(AliVCluster *cluster, AliVCaloCells *cells)         ;

  Float_t                     fCellEnergyCut;            ///< Energy cell cut
  Bool_t                      fParticleLevel;            ///< Set particle level analysis
  Bool_t                      fIsMC;                     ///< MC analysis
  TString                     fCentMethod2;              ///< Centrality method 2
  TString                     fCentMethod3;              ///< Centrality method 3
  Int_t                       fDoV0QA;                   ///< Add V0 QA histograms
  Int_t                       fDoEPQA;                   ///< Add event plane QA histograms
  Int_t                       fDoLeadingObjectPosition;  ///< Add axis for leading object position (eta-phi)
  Int_t                       fMaxCellsInCluster;        ///< Maximum number (approx) of cells in a cluster
  Float_t                     fPtBinWidth;               ///< Histogram pt bin width
  Float_t                     fMaxPt;                    ///< Histogram pt limit
  Bool_t                      fSeparateEMCalDCal;        ///< Separate EMCal from DCal in QA plots
  Bool_t                      fIsEmbedded;               ///< Embedded data present
  Double_t                    fCent2;                    //!<!Event centrality with method 2
  Double_t                    fCent3;                    //!<!Event centrality with method 3
  AliVVZERO                  *fVZERO;                    //!<!Event V0 object
  Double_t                    fV0ATotMult;               //!<!Event V0A total multiplicity
  Double_t                    fV0CTotMult;               //!<!Event V0C total multiplicity
  Int_t                       fNTotClusters[2];          //!<!Total number of accepted clusters in current event (DCal/EMCal)
  AliTLorentzVector           fLeadingCluster[2];        //!<!Leading cluster in current event (EMCal/DCal)
  Int_t                       fNTotTracks;               //!<!Total number of accepted tracks in current event
  AliTLorentzVector           fLeadingTrack;             //!<!Leading track in current event

  THistManager                fHistManager;              //!< Histogram manager

private:
  AliAnalysisTaskEmcalJetQA(const AliAnalysisTaskEmcalJetQA&);            // not implemented
  AliAnalysisTaskEmcalJetQA &operator=(const AliAnalysisTaskEmcalJetQA&); // not implemented

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskEmcalJetQA, 2)
  /// \endcond
};
#endif
