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
class TH1I; 
class TList;

class AliAnalysisTaskHFE : public AliAnalysisTask{
  enum{
    kIsSecVtxOn = BIT(19),
    kIsPriVtxOn = BIT(20)
  };
  public:
  enum{
    kPIDqa = 0,
    kCUTqa = 1,
    kMCqa = 2
  };
    AliAnalysisTaskHFE();
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
    Int_t IsSignalElectron(AliESDtrack *) const;
    void SetQAOn(Int_t qaLevel) { SETBIT(fQAlevel, qaLevel); };
    void SetPriVtxOn()        { SetBit(kIsPriVtxOn, kTRUE); };
    void SetSecVtxOn()        { SetBit(kIsSecVtxOn, kTRUE); };
    void SetPIDdetectors(Char_t *detectors){ fPIDdetectors = detectors; }
    void AddPIDdetector(Char_t *detector);
    void PrintStatus();
 
  private:
    class LabelContainer{
      public:
        LabelContainer(Int_t capacity);
        ~LabelContainer() {delete[] fContainer; };

        Bool_t Append(Int_t label);
        Bool_t Find(Int_t Label);
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
    AliESDEvent *fESD;                    //! The ESD Event
    AliMCEvent *fMC;                      //! The MC Event
    AliCFManager *fCFM;                   //! Correction Framework Manager
    TList *fCorrelation;                  //! response matrix for unfolding  
    THnSparseF *fPIDperformance;          //! info on contamination and yield of electron spectra
    AliHFEpid *fPID;                      //! PID
    AliHFEcuts *fCuts;                    //! Cut Collection
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

