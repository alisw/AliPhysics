#ifndef ALIANALYSISTASKEMCALJETDIJETMASS_H
#define ALIANALYSISTASKEMCALJETDIJETMASS_H
/**
 * \file AliAnalysisTaskEmcalJetDijetMass.h
 * \brief Declaration of class AliAnalysisTaskEmcalJetDijetMass
 *
 * In this header file the class AliAnalysisTaskEmcalJetDijetMass is declared.
 *
 * \author Oskari Saarimäki <oskari.antti.matti.saarimaki@cern.ch>, Jyväskylä University
 * \date Jul 10, 2020
 */

/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <iomanip>
#include <TDirectory.h>
#include <TComplex.h>
#include <AliLog.h>
#include <AliAnalysisTaskSE.h>
#include "AliAnalysisTaskEmcalJet.h"
#include "THistManager.h"
//#include "AliAnalysisTaskEmcalJetDijetMassHisto.h"
#include "AliJCDijetHistos.h"
#include "AliJEmcalDijetAna.h"
//#include "AliAnalysisTaskEmcalJetDijetMassHistos.h"

/**
 * \class AliAnalysisTaskEmcalJetDijetMass
 * \brief Implementation of a dijet invariant mass class
 *
 * Dijet mass analysis.
 */
class AliAnalysisTaskEmcalJetDijetMass : public AliAnalysisTaskEmcalJet {
 public:

  AliAnalysisTaskEmcalJetDijetMass()                                               ;
  AliAnalysisTaskEmcalJetDijetMass(const char *name)                               ;
  virtual ~AliAnalysisTaskEmcalJetDijetMass()                                      ;

  void                        UserCreateOutputObjects()                         ;
  void                        Terminate(Option_t *option)                       ;
  //void    SetCuts(double particleEta,
  //        double particlePt,
  //        double leadingJet,
  //        double subleadingJet,
  //        double constituent,
  //        double deltaPhi,
  //        double matchingR,
  //        double trackingIneff,
  //        double minJetPt) {
  //    fparticleEtaCut=particleEta;
  //    fparticlePtCut=particlePt;
  //    fleadingJetCut=leadingJet;
  //    fsubleadingJetCut=subleadingJet;
  //    fconstituentCut=constituent;
  //    fdeltaPhiCut=deltaPhi;
  //    fmatchingR = matchingR;
  //    ftrackingIneff = trackingIneff;
  //    fMinJetPt = minJetPt;
 // }

static AliAnalysisTaskEmcalJetDijetMass *AddTaskEmcalJetDijetMass(
                                    const char *ntracks,
                                    const char *nclusters,
                                    const char* ncells,
                                    const char *suffix,
                                    TString taskName,
                                    Bool_t isMC,
                                    TString sJCatalyst,
                                    TString sJCatalystDetMC,
                                    UInt_t flags,
                                    TString centBins,
                                    double jetCone,
                                    double ktjetCone,
                                    int ktScheme,
                                    int antiktScheme,
                                    Bool_t usePionMass,
                                    Bool_t useDeltaPhiBGSubtr,
                                    double particleEtaCut,
                                    double particlePtCut,
                                    double leadingJetCut,
                                    double subleadingJetCut,
                                    double minJetPt,
                                    double constituentCut,
                                    double deltaPhiCut,
                                    double matchingR,
                                    double trackingIneff);

 protected:
  void                        ExecOnce()                                        ;
  Bool_t                      FillHistograms()                                  ;
  Bool_t                      Run()                                             ;

  void                        AllocateJetHistograms()                           ;
  void                        AllocateTrackHistograms()                         ;
//void                        AllocateClusterHistograms()                       ;
//void                        AllocateCellHistograms()                          ;

  void                        DoJetLoop()                                       ;
  void                        DoTrackLoop()                                     ;
//void                        DoClusterLoop()                                   ;
//void                        DoCellLoop()                                      ;

  THistManager                fHistManager                                      ;///< Histogram manager
//AliAnalysisTaskEmcalJetDijetMassHisto *fhistos;
  AliJEmcalDijetAna *fana;
  //AliJCDijetAna *fanaMC;
  //AliJCDijetHistos *fhistos;  ///<
  //AliAnalysisTaskEmcalJetDijetMassHistos *fEmcalHistos; ///<
  TH1D* hTest; ///<
  //AliJCDijetHistos *fhistosDetMC;

 private:
  AliAnalysisTaskEmcalJetDijetMass(const AliAnalysisTaskEmcalJetDijetMass&)           ; // not implemented
  AliAnalysisTaskEmcalJetDijetMass &operator=(const AliAnalysisTaskEmcalJetDijetMass&); // not implemented
  TDirectory *fJOutput; // Output directory
  double fjetCone;
  double fktJetCone;
  int  fktScheme;
  int  fantiktScheme;
  Bool_t fusePionMass;
  Bool_t fuseDeltaPhiBGSubtr;
  Bool_t fIsMC;       // MC data or real data
  double fparticleEtaCut;
  double fparticlePtCut;
  double fleadingJetCut;
  double fsubleadingJetCut;
  double fMinJetPt;
  double fconstituentCut;
  double fdeltaPhiCut;
  double fmatchingR;
  double ftrackingIneff;
  //TH1F *hBinning;
  //TH1F *hBinningALICE;


  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskEmcalJetDijetMass, 7);
  /// \endcond
};
#endif
