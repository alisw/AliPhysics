#ifndef ALIANALYSISTASKEMCALJETSPECTRA8TEVTRIGGERQA_H
#define ALIANALYSISTASKEMCALJETSPECTRA8TEVTRIGGERQA_H

/**********************************************************************************
 * Copyright (C) 2018, Copyright Holders of the ALICE Collaboration                *
 * All rights reserved.                                                            *
 *                                                                                 *
 * Redistribution and use in source and binary forms, with or without              *
 * modification, are permitted provided that the following conditions are met:     *
 *   * Redistributions of source code must retain the above copyright              *
 *     notice, this list of conditions and the following disclaimer.               *
 *   * Redistributions in binary form must reproduce the above copyright           *
 *     notice, this list of conditions and the following disclaimer in the         *
 *     documentation and/or other materials provided with the distribution.        *
 *   * Neither the name of the <organization> nor the                              *
 *     names of its contributors may be used to endorse or promote products        *
 *     derived from this software without specific prior written permission.       *
 *                                                                                 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND *
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED   *
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE          *
 * DISCLAIMED. IN NO EVENT SHALL ALICE COLLABORATION BE LIABLE FOR ANY             *
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES      *
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;    *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND     *
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT      *
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS   *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                    *
 * *********************************************************************************/

/**
 * \file AliAnalysisTaskEmcalJetSpectra8TeVTriggerQA.h
 * \brief Declaration of class AliAnalysisTaskEmcalJetSpectra8TeVTriggerQA
 *
 * In this header file the class AliAnalysisTaskEmcalJetSpectra8TeVTriggerQA is declared.
 * This is a jet task that calculates the EMCal triggered spectra and does
 * additional QA for the EMCal trigger for the 8 TeV jet spectra analysis
 *
 * \author Andrew Castro <andrew.john.castro@cern.ch>, University of Tennessee
 * \date April 14, 2018
 */

/* Copyright(c) 1998-2018, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */



#include "AliAnalysisTaskEmcalJet.h"
#include "THistManager.h"
#include <string>
#include <vector>
#include "AliAnalysisTaskEmcalTriggerBase.h"
#include "AliCutValueRange.h"
#include "AliEventCuts.h"
#include <TCustomBinning.h>
#include <TString.h>

class AliOADBContainer;
class AliEMCALTriggerPatchInfo;
class THistManager;
class TObjArray;

/**
 * \class AliAnalysisTaskEmcalJetSpectra8TeVTriggerQA
 * \brief Implementation of a EMCal spectra task and QA for EMCal triggers
 *
 * This class in an implementation of a trigger QA for the EMCal
 * It derives from AliAnalysisTaskEmcalJet.
 * It also performs a QA of the cluster-track matching.
 * Note: if jets are not used this class can be simplified by deriving
 * from AliAnalysisTaskEmcal and removing the functions DoJetLoop()
 * and AllocateJetHistograms().
 */
class AliAnalysisTaskEmcalJetSpectra8TeVTriggerQA : public AliAnalysisTaskEmcalJet {
public:

    AliAnalysisTaskEmcalJetSpectra8TeVTriggerQA()                                               ;
    AliAnalysisTaskEmcalJetSpectra8TeVTriggerQA(const char *name)                               ;
    virtual ~AliAnalysisTaskEmcalJetSpectra8TeVTriggerQA()                                      ;

    void                        UserCreateOutputObjects()                         ;
    void                        Terminate(Option_t *option)                       ;
    void                        ExtractMainPatch()                                ;
    
    //setters
    void                        SetUseSumw2(Bool_t b) { fUseSumw2          = b   ; }

protected:
    void                        ExecOnce()                                        ;
    Bool_t                      FillHistograms()                                  ;
    Bool_t                      Run()                                             ;

    void                        AllocateJetHistograms()                           ;///< Jet Histograms
    void                        AllocateTrackHistograms()                         ;///< ITS-TPC Track Histograms
    void                        AllocateClusterHistograms()                       ;///< EMCal Cluster Histograms
    void                        AllocateCellHistograms()                          ;///< EMCal Tower Histograms
    void                        AllocateParticleHistograms()                      ;///< Generator Level MC Histograms
    //void                        AllocateMCJetHistograms()                         ;///< MC truth Jet Histograms
    
    

    void                        DoJetLoop()                                       ;
    void                        DoTrackLoop()                                     ;
    void                        DoClusterLoop()                                   ;
    void                        DoCellLoop()                                      ;
    //void                        DoParticleLoop()                                  ;
    //void                        DoMCJetLoop()                                     ;
    //Bool_t                      IsLEDEvent() const                                ;
    
    
    Bool_t                      fUseRecalcPatches                                 ;///<                  Switch between offline (FEE) and recalc (L1) patches
    Bool_t                      SelectSingleShowerPatch(const AliEMCALTriggerPatchInfo *patch) const;
    Bool_t                      SelectJetPatch(const AliEMCALTriggerPatchInfo  *patch) const;
    THistManager                fHistManager                                      ;///<                 Histogram manager
    
    // Sparse Definition
    virtual THnSparse*     NewTHnSparseF(const char* name, UInt_t entries);
    virtual void           GetDimParams(Int_t iEntry,TString &label, Int_t &nbins, Double_t &xmin, Double_t &xmax);
    // Binning helper functions
    
    //  Utilities


private:
    

    Bool_t                        fUseSumw2                  ;//!          activate sumw2 for output histograms

    TH1F                          *fHistNumbJets             ;//!          Numb Jets Per Event
    TH1F                          *fHistJetPt                ;//!          Jet Pt Dist
    TH1F                          *fHistJetJetPatchE         ;//!          Jet - Jet Trigger Patch E
    TH1F                          *fHistJetGammaPatchE       ;//!          Jet - Gamma Trigger Patch E
    TH1F                          *fHistJetJetPatchPt        ;//!          Jet - Jet Trigger Patch Pt
    TH1F                          *fHistJetGammaPatchPt      ;//!          Jet - Gamma Trigger Patch Pt
    TH1F                          *fHistTriggerPatchE        ;//!          EMCal Trigger Patch E


    AliAnalysisTaskEmcalJetSpectra8TeVTriggerQA(const AliAnalysisTaskEmcalJetSpectra8TeVTriggerQA&)           ; // not implemented
    AliAnalysisTaskEmcalJetSpectra8TeVTriggerQA &operator=(const AliAnalysisTaskEmcalJetSpectra8TeVTriggerQA&); // not implemented

    /// \cond CLASSIMP
    ClassDef(AliAnalysisTaskEmcalJetSpectra8TeVTriggerQA, 7);
    /// \endcond
};
#endif
