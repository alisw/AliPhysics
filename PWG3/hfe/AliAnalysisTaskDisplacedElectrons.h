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
// Analysis task: 
// study displaced electrons from beauty and charm 
// with cut on impact parameters in various pT bins
// 


#ifndef ALIANALYSISTASKDISPLACEDELECTRONS_H
#define ALIANALYSISTASKDISPLACEDELECTRONS_H

#ifndef ALIANALYSISTASKSE_H
#include "AliAnalysisTaskSE.h"
#endif

#ifndef ROOT_THnSparse
#include <THnSparse.h>
#endif

class TH1I; 
class TH1F;
class TList;
class AliLog;

class AliCFManager;
class AliESDEvent;
class AliESDtrackCuts;
class AliMCEvent;
class AliVParticle;

class AliStack;

class AliHFEpid;
class AliHFEcuts;
class AliHFEdisplacedElectrons;

class AliAnalysisTaskDisplacedElectrons : public AliAnalysisTaskSE{
 public:

  typedef enum{
    kPostProcess = 0,
    kDisplacedElectrons = 1, 
    kCorrection = 2
  }Switches_t;

  enum{
    kHasMCdata = BIT(19),
    kAODanalysis = BIT(20)
  };
  
  AliAnalysisTaskDisplacedElectrons();
  AliAnalysisTaskDisplacedElectrons(const char * name);
  AliAnalysisTaskDisplacedElectrons(const AliAnalysisTaskDisplacedElectrons &ref);
  AliAnalysisTaskDisplacedElectrons& operator=(const AliAnalysisTaskDisplacedElectrons &ref);
  virtual ~AliAnalysisTaskDisplacedElectrons();
  
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *);
  //  virtual void ConnectInputData(Option_t *);
  virtual void Terminate(Option_t *);

  void PrintStatus() const;
  
  Bool_t GetPlugin(Int_t plug) const { return TESTBIT(fPlugins, plug); };
  void SwitchOnPlugin(Int_t plug);

  void SetHFECuts(AliHFEcuts * const cuts) { fCuts = cuts; };
  
  void SetPIDdetectors(Char_t * const detectors){ fPIDdetectors = detectors; };
  void SetPIDStrategy(UInt_t strategy) { fPIDstrategy = strategy; };
  void SetDBLevel(UInt_t debugLevel) { fDebugLevel = debugLevel; };
  void AddPIDdetector(TString detector); 
 
  Bool_t HasMCData() const { return TestBit(kHasMCdata); };
  void SetHasMCData(Bool_t hasMC = kTRUE) { SetBit(kHasMCdata, hasMC); };
  
  
 private:
  
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
  void MakeEventContainer();

  UInt_t fDebugLevel;                   // debug level
  
  TString fPIDdetectors;                // Detectors for Particle Identification
  UInt_t fPIDstrategy;                  // PID Strategy
  

  UShort_t fPlugins;                    // Enabled Plugins                                    
  AliESDEvent *fESD;                    //! The ESD Event
  AliMCEvent *fMC;                      //! The MC Event

  AliHFEcuts *fCuts;                    // Cut Collection
  AliHFEpid *fPID;                      // PID method
  AliCFManager *fCFM;                   //! Correction Framework Manager
                                               
  TH1I *fNEvents;                       //! counter for the number of Events
  TH1F *fElectronsPt;                       //! pt distribution of electrons (hfepid)
  TList *fOutput;                       //! Container for this Task Output
  TList *fCorrection;                       //! Container for correction  Output
  AliHFEdisplacedElectrons *fDisplacedElectrons;        //! HFE displaced Electrons pointer 
  TList *fHistDisplacedElectrons;                      //! list of outputs
 
  ClassDef(AliAnalysisTaskDisplacedElectrons, 1);      // The DisplacedElectrons Analysis Task
};
#endif

