/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/
//
// Task for Heavy Flavour Electron Analysis
// Fills a single-inclusive electron pt-spectrum
// For further information see implementation file
//
#ifndef ALIANALYSISTASKHFE_H
#define ALIANALYSISTASKHFE_H

#ifndef ALIANALYSISTASKSE_H
#include "AliAnalysisTaskSE.h"
#endif

#ifndef ROOT_THnSparse
#include <THnSparse.h>
#endif

class AliHFEpid;
class AliHFEcuts;
class AliHFEmcQA;
class AliHFEsecVtx;
class AliHFEelecbackground;
class AliHFEcollection;
class AliCFManager;
class AliVEvent;
class AliMCEvent;
class AliVParticle;
class AliTriggerAnalysis;
class TH1I; 
class TList;
class TH3D;
class TF1;

class AliAnalysisTaskHFE : public AliAnalysisTaskSE{
  public:
    enum{
      kPIDqa = 0,
      kMCqa =1 
    };
    enum{
      kPriVtx = 0,
      kSecVtx = 1,
      kIsElecBackGround = 2,
      kPostProcess = 3
    };
    enum CreationProcess_t{
      kSignalCharm = 0,
      kSignalBeauty = 1,
      kGammaConv = 2,
      kOther = 3
    };
    AliAnalysisTaskHFE();
    AliAnalysisTaskHFE(const char * name);
    AliAnalysisTaskHFE(const AliAnalysisTaskHFE &ref);
    AliAnalysisTaskHFE& operator=(const AliAnalysisTaskHFE &ref);
    virtual void Copy(TObject &o) const;
    virtual ~AliAnalysisTaskHFE();

    virtual void UserCreateOutputObjects();
    virtual void UserExec(Option_t *);
    virtual void Terminate(Option_t *);

    virtual Bool_t IsEventInBinZero();

    Bool_t IsQAOn(Int_t qaLevel) const { return TESTBIT(fQAlevel, qaLevel); };
    Bool_t IsAODanalysis() const { return TestBit(kAODanalysis); };
    Bool_t IsESDanalysis() const { return !TestBit(kAODanalysis); };
    Bool_t HasMCData() const { return TestBit(kHasMCdata); }
    Bool_t GetPlugin(Int_t plug) const { return TESTBIT(fPlugins, plug); };
    void SetHFECuts(AliHFEcuts * const cuts) { fCuts = cuts; };
    void SetHFECutsPreselect(AliHFEcuts * const cuts) { fCutspreselect = cuts; };
    void SetHFEElecBackGround(AliHFEelecbackground * const elecBackGround) { fElecBackGround = elecBackGround; };
    void SetQAOn(Int_t qaLevel) { SETBIT(fQAlevel, qaLevel); };
    void SwitchOnPlugin(Int_t plug);
    void SetHasMCData(Bool_t hasMC = kTRUE) { SetBit(kHasMCdata, hasMC); };
    void SetPIDdetectors(Char_t * const detectors){ fPIDdetectors = detectors; }
    void SetPIDStrategy(UInt_t strategy) { fPIDstrategy = strategy; }
    void SetPIDPreselect(AliHFEpid * const cuts) { fPIDpreselect = cuts; };
    void AddPIDdetector(TString detector);
    void SetAODAnalysis() { SetBit(kAODanalysis, kTRUE); };
    void SetESDAnalysis() { SetBit(kAODanalysis, kFALSE); };
    void SetWeightFactors(TH3D * const weightFactors);
    void SetWeightFactorsFunction(TF1 * const weightFactorsFunction);
    void SetBackGroundFactorsFunction(TF1 * const backGroundFactorsFunction) { fBackGroundFactorsFunction = backGroundFactorsFunction; };
    void PrintStatus() const;
 
  private:
    enum{
      kHasMCdata = BIT(19),
      kAODanalysis = BIT(20)
    };
    class LabelContainer{
      public:
        LabelContainer(Int_t capacity);
        ~LabelContainer() {delete[] fContainer; };

        Bool_t Append(Int_t label);
        Bool_t Find(Int_t Label) const;
        Int_t Next();
        void ResetIterator(){ fCurrent = fBegin; }

      private:
        LabelContainer(const LabelContainer &);
        LabelContainer &operator=(const LabelContainer &);
        Int_t *fContainer;    // the Container for the labels
        Int_t *fBegin;        // Pointer to the first entry
        Int_t *fEnd;          // Pointer to the end of the container
        Int_t *fLast;         // Pointer to the last entry
        Int_t *fCurrent;      // Current entry to mimic an iterator
    };

    Bool_t IsGammaElectron(const AliVParticle * const track) const;
    Int_t  IsSignalElectron(const AliVParticle * const track) const;
    Bool_t FillProductionVertex(const AliVParticle * const track) const;
    Double_t FindWeight(Double_t pt, Double_t eta, Double_t phi) const;
    void MakeParticleContainer();
    void MakeEventContainer();
    void ProcessMC();
    void ProcessESD();
    void ProcessAOD();
    void FilterTaggedTrack(AliESDtrack *track, Int_t species);
    Bool_t PreSelectTrack(AliESDtrack *track) const;
    Bool_t ProcessMCtrack(AliVParticle *track);
    Bool_t ProcessCutStep(Int_t cutStep, AliVParticle *track, Double_t *container, Bool_t signal, Bool_t alreadyseen,Double_t weight);
    
    ULong_t fQAlevel;                     // QA level
    TString fPIDdetectors;                // Detectors for Particle Identification
    UInt_t fPIDstrategy;                  // PID Strategy
    Double_t fTPCBetheBlochParameters[5]; // TPC Bethe-Bloch Parameters
    UShort_t fPlugins;                    // Enabled Plugins
    Bool_t  fWeighting;                   // Weighting or not for the efficiency maps
    TH3D *fWeightFactors;                 // Weight factors
    TF1  *fWeightFactorsFunction;         // Weight factors
    TF1  *fBackGroundFactorsFunction;     // BackGround factors
    AliCFManager *fCFM;                   //! Correction Framework Manager
    AliCFManager *fV0CF;                  //! Correction Framework Manager for V0-tagged tracks
    AliTriggerAnalysis *fTriggerAnalysis; //! Trigger Analysis for Normalisation
    AliCFContainer * fHadronicBackground; //! Container for hadronic Background
    TList *fCorrelation;                  //! response matrix for unfolding  
    THnSparseF *fPIDperformance;          //! info on contamination and yield of electron spectra
    THnSparseF *fSignalToBackgroundMC;    //! Signal To Background Studies on pure MC information
    AliHFEpid *fPID;                      //! PID
    AliHFEpid *fPIDtagged;                // PID oject for tagged tracks (identical clone without QA)
    AliHFEpid *fPIDpreselect;             // PID oject for pre-selected tracks (without QA)
    AliHFEcuts *fCuts;                    // Cut Collection
    AliHFEcuts *fCutsTagged;              // Cut Collection for tagged tracks
    AliHFEcuts *fCutspreselect;           // Cut Collection for pre-selected tracks
    AliHFEsecVtx *fSecVtx;                //! Secondary Vertex Analysis
    AliHFEelecbackground *fElecBackGround;//! Background analysis
    AliHFEmcQA *fMCQA;                    //! MC QA
    TH1I *fNEvents;                       //! counter for the number of Events
    TH1I *fNElectronTracksEvent;          //! Number of Electron candidates after PID decision per Event
    TList *fQA;                           //! QA histos for the cuts
    TList *fOutput;                       //! Container for Task Output
    TList *fHistMCQA;                     //! Output container for MC QA histograms 
    TList *fHistSECVTX;                   //! Output container for sec. vertexing results
    TList *fHistELECBACKGROUND;           //! Output container for electron background analysis
    AliHFEcollection *fQACollection;      //! Tasks own QA collection

    ClassDef(AliAnalysisTaskHFE, 2)       // The electron Analysis Task
};
#endif

