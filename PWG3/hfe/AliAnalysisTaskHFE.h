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

#ifndef ALIANALYSISTASK_H
#include "AliAnalysisTask.h"
#endif

#ifndef ROOT_THnSparse
#include <THnSparse.h>
#endif

class AliHFEpid;
class AliHFEcuts;
class AliHFEmcQA;
class AliHFEsecVtx;
class AliCFManager;
class AliESDEvent;
class AliESDtrackCuts;
class AliMCEvent;
class AliVParticle;
class TH1I; 
class TList;

class AliAnalysisTaskHFE : public AliAnalysisTask{
  public:
  enum{
    kPIDqa = 0,
    kMCqa =1 
  };
    AliAnalysisTaskHFE();
    AliAnalysisTaskHFE(const char * name);
    AliAnalysisTaskHFE(const AliAnalysisTaskHFE &ref);
    AliAnalysisTaskHFE& operator=(const AliAnalysisTaskHFE &ref);
    virtual ~AliAnalysisTaskHFE();

    virtual void ConnectInputData(Option_t *);
    virtual void CreateOutputObjects();
    virtual void Exec(Option_t *);
    virtual void Terminate(Option_t *);

    Bool_t IsQAOn(Int_t qaLevel) const { return TESTBIT(fQAlevel, qaLevel); };
    Bool_t IsSecVtxOn() const { return TestBit(kIsSecVtxOn); };
    Bool_t IsPriVtxOn() const { return TestBit(kIsPriVtxOn); };
    Bool_t IsRunningPostProcess() const { return TestBit(kIsRunningPostProcess); };
    Bool_t HasMCData() const { return TestBit(kHasMCdata); }
    Int_t IsSignalElectron(AliVParticle *fTrack) const;
    void Load(TString filename = "HFEtask.root");
    void PostProcess();
    void SetHFECuts(AliHFEcuts * const cuts) { fCuts = cuts; };
    void SetQAOn(Int_t qaLevel) { SETBIT(fQAlevel, qaLevel); };
    void SetHasMCData(Bool_t hasMC = kTRUE) { SetBit(kHasMCdata, hasMC); };
    void SetPriVtxOn(Bool_t option = kTRUE)        { SetBit(kIsPriVtxOn, option); };
    void SetSecVtxOn(Bool_t option = kTRUE)        { SetBit(kIsSecVtxOn, option); };
    void SetRunPostProcess(Bool_t option = kTRUE)  { SetBit(kIsRunningPostProcess, option); };
    void SetPIDdetectors(Char_t * const detectors){ fPIDdetectors = detectors; }
    void SetPIDStrategy(UInt_t strategy) { fPIDstrategy = strategy; }
    void AddPIDdetector(Char_t *detector);
    void PrintStatus() const;
    Float_t GetRapidity(TParticle *part) const;
 
  private:
    enum{
      kIsSecVtxOn = BIT(19),
      kIsPriVtxOn = BIT(20),
      kIsRunningPostProcess = BIT(21),
      kHasMCdata = BIT(22)
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
    void MakeParticleContainer();
    
    ULong_t fQAlevel;                     // QA level
    TString fPIDdetectors;                // Detectors for Particle Identification
    UInt_t fPIDstrategy;                  // PID Strategy
    AliESDEvent *fESD;                    //! The ESD Event
    AliMCEvent *fMC;                      //! The MC Event
    AliCFManager *fCFM;                   //! Correction Framework Manager
    TList *fCorrelation;                  //! response matrix for unfolding  
    THnSparseF *fPIDperformance;          //! info on contamination and yield of electron spectra
    AliHFEpid *fPID;                      //! PID
    AliHFEcuts *fCuts;                    // Cut Collection
    AliHFEsecVtx *fSecVtx;                //! Secondary Vertex Analysis
    AliHFEmcQA *fMCQA;                    //! MC QA
    TH1I *fNEvents;                       //! counter for the number of Events
    TH1I *fNElectronTracksEvent;          //! Number of Electron candidates after PID decision per Event
    TList *fQA;                           //! QA histos for the cuts
    TList *fOutput;                       //! Container for Task Output
    TList *fHistMCQA;                     //! Output container for MC QA histograms 
    TList *fHistSECVTX;                   //! Output container for sec. vertexing results

    ClassDef(AliAnalysisTaskHFE, 1)       // The electron Analysis Task
};
#endif

