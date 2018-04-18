#ifndef ALIANALYSISTASKEMCALJETSPECTRA8TEVTRIGGERQA_H
#define ALIANALYSISTASKEMCALJETSPECTRA8TEVTRIGGERQA_H
/**
 * \file AliAnalysisTaskEmcalJetSpectra8TeVTriggerQA.h
 * \brief Declaration of class AliAnalysisTaskEmcalJetSpectra8TeVTriggerQA
 *
 * In this header file the class AliAnalysisTaskEmcalJetSpectra8TeVTriggerQA is declared.
 * This is a jet task that calculates the EMCal trigger spectra and does
 * additional QA for the EMCal trigger for the 8 TeV jet spectra analysis
 *
 * \author Andrew Castro <andrew.john.castro@cern.ch>, University of Tennessee
 * \date April 18, 2018
 */

/* Copyright(c) 1998-2018, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */



#include "AliAnalysisTaskEmcalJet.h"
#include "THistManager.h"
#include <string>
#include <vector>
#include "Tracks/AliAnalysisTaskEmcalTriggerBase.h"
#include "Tracks/AliCutValueRange.h"
#include <TCustomBinning.h>
#include <TString.h>

class AliOADBContainer;
class AliEMCALTriggerPatchInfo;
class THistManager;
class TObjArray;
class TFormula;
class TH2F;
class TH1F;
class TF1;
class THnSparse;
class TRandom3;
class TObjArray;
class TClonesArray;
class TObject;
class TString;
class TProfile2D;
class AliAODEvent;
class AliESDEvent;
class AliMCEvent;
class AliVEvent;
class AliStack;
class AliEMCALGeometry;
//class AliEMCalJet;
class AliEMCALRecoUtils;
class AliESDCaloCluster;
class AliVTrack;

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

    void                        AllocateJetHistograms()                           ;
    void                        AllocateTrackHistograms()                         ;
    void                        AllocateClusterHistograms()                       ;
    void                        AllocateCellHistograms()                          ;

    void                        DoJetLoop()                                       ;
    void                        DoTrackLoop()                                     ;
    void                        DoClusterLoop()                                   ;
    void                        DoCellLoop()                                      ;
    Bool_t                      IsLEDEvent() const                                ;
    Bool_t                      fUseRecalcPatches                                 ;///<                  Switch between offline (FEE) and recalc (L1) patches
    Bool_t                      SelectSingleShowerPatch(const AliEMCALTriggerPatchInfo *patch) const;
    Bool_t                      SelectJetPatch(const AliEMCALTriggerPatchInfo  *patch) const;
    THistManager                fHistManager                                      ;//!<!                 Histogram manager
    
    //  Utilities
    Double_t                    GetZ(const Double_t trkPx, const Double_t trkPy, const Double_t trkPz, const Double_t jetPx, const Double_t jetPy, const Double_t jetPz) const;
    Double_t                    GetSmearedTrackPt(AliVTrack *track)                ; //!<!               smeared Pt resolution with gaussian
    Double_t                    GetFcross(const AliVCluster *cluster, AliVCaloCells *cells);
    
    AliVEvent                   *fRecevent                                         ;//!<!                Reconstructed event
    AliMCEvent                  *fMCevent                                          ;//!<!                Monte-Carlo event

private:
    
    //AliEMCALGeometry              *fGeom                     ;//!<!          EMCal goemetry utility
    AliEMCALRecoUtils             *fRecoUtil                 ;//!<!          Reco utility
    TF1                           *fClusterEResolution       ;//!<!          Parameterization of cluster energy resolution from 2010 test beam results a = 4.35 b = 9.07 c = 1.63
    Double_t                      fVaryTrkPtRes              ;//!<!         Variation of tracking momentum resolution
    Bool_t                        fUseSumw2                  ;//!<!         activate sumw2 for output histograms

    TH1F                          *fHistNumbJets;//!<!                     Numb Jets Per Event
    TH1F                          *fHistJetPt;//!<!                        Jet Pt Dist
    TH1F                          *fHistJetJetPatchE;//!<!                 Jet - Jet Trigger Patch E
    TH1F                          *fHistJetGammaPatchE;//!<!               Jet-Gamma Trigger Patch E
    TH1F                          *fHistJetJetPatchPt;//!<!                Jet - Jet Trigger Patch Pt
    TH1F                          *fHistJetGammaPatchPt;//!<!              Jet-Gamma Trigger Patch Pt
    TH1F                          *fHistTriggerPatchE;//!<!                EMCal Trigger Patch E
    TH1F                          *fHistEMCalTowerMult[9];//!<!            EMCal Tower Multiplicity by SM

    AliAnalysisTaskEmcalJetSpectra8TeVTriggerQA(const AliAnalysisTaskEmcalJetSpectra8TeVTriggerQA&)           ; // not implemented
    AliAnalysisTaskEmcalJetSpectra8TeVTriggerQA &operator=(const AliAnalysisTaskEmcalJetSpectra8TeVTriggerQA&); // not implemented

    /// \cond CLASSIMP
    ClassDef(AliAnalysisTaskEmcalJetSpectra8TeVTriggerQA, 7);
    /// \endcond
};
#endif
