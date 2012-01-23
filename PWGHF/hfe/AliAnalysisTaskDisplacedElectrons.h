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
//
// Authors:
//  Hongyan Yang <hongyan@physi.uni-heidelberg.de>
//  Carlo Bombonati <Carlo.Bombonati@cern.ch>
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
    kCorrection = 2,
    kDePidQA = 3
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
  
  Bool_t GetPlugin(Int_t plug) const { return TESTBIT(fDePlugins, plug); };
  void SwitchOnPlugin(Int_t plug);

  void SetHFECuts(AliHFEcuts * const cuts) { fDeCuts = cuts; };
  void SetNclustersITS(Int_t nITSclusters){fNminITSCluster = nITSclusters;};
  void SetMinPrimVtxContrib(Int_t nPrimVtxContrib){fNminPrimVtxContrib = nPrimVtxContrib;};
  void SetPIDdetectors(Char_t * const detectors){ fDePIDdetectors = detectors; };
  void SetPIDStrategy(UInt_t strategy) { fDePIDstrategy = strategy; };
  void SetDBLevel(UInt_t debugLevel) { fDeDebugLevel = debugLevel; };
  void AddPIDdetector(TString detector); 
 
  Bool_t IsAODanalysis() const { return TestBit(kAODanalysis); };
  Bool_t IsESDanalysis() const { return !TestBit(kAODanalysis); };
  Bool_t HasMCData() const { return TestBit(kHasMCdata); }
  void SetHasMCData(Bool_t hasMC = kTRUE) { SetBit(kHasMCdata, hasMC); };
  
  void SetAODAnalysis() { SetBit(kAODanalysis, kTRUE); };
  void SetESDAnalysis() { SetBit(kAODanalysis, kFALSE); };

  void ProcessMC();
  void ProcessESD();
  void ProcessData();
  

  
 private:
  
  class AliLabelContainer{
  public:
    AliLabelContainer(Int_t capacity);
    ~AliLabelContainer() {delete[] fContainer; };
    
    Bool_t Append(Int_t label);
    Bool_t Find(Int_t Label) const;
    Int_t Next();
    void ResetIterator(){ fCurrent = fBegin; }
    
  private:
    AliLabelContainer(const AliLabelContainer &);
    AliLabelContainer &operator=(const AliLabelContainer &);
    Int_t *fContainer;    // the Container for the labels
    Int_t *fBegin;        // Pointer to the first entry
    Int_t *fEnd;          // Pointer to the end of the container
    Int_t *fLast;         // Pointer to the last entry
    Int_t *fCurrent;      // Current entry to mimic an iterator
  };
    
  void MakeParticleContainer();
  void MakeEventContainer();

  UInt_t fDeDebugLevel;                  // debug level
  Int_t fNminITSCluster;                 // number of clusters in ITS
  Int_t fNminPrimVtxContrib;             // number of ncontributor in ITS for prim vtx
  TString fDePIDdetectors;                // Detectors for Particle Identification
  UInt_t fDePIDstrategy;                  // PID Strategy

  UShort_t fDePlugins;                    // Enabled Plugins    

  AliHFEcuts *fDeCuts;                    // Cut Collection
  AliHFEpid *fDePID;                      //! PID method
  AliCFManager *fDeCFM;                   //! Correction Framework Manager
  AliHFEdisplacedElectrons *fDisplacedElectrons;        //! HFE displaced Electrons pointer 
                                
  TH1I *fDeNEvents;                       //! counter for the number of Events
  TH1F *fElectronsMcPt;                   //! pt distribution of MC electrons (mcpid)
  TH1F *fElectronsEsdPt;                  //! pt distribution of ESD electrons (hfepid)
  TH1F *fElectronsDataPt;                 //! pt distribution of DATA electrons (hfepid)
  TList *fDeCorrection;                   //! Container for correction  Outpu  
  TList *fDeQA;                          //! container for the PID qa 
  TList *fHistDisplacedElectrons;                      //! list of outputs
 
  ClassDef(AliAnalysisTaskDisplacedElectrons, 1);      // The DisplacedElectrons Analysis Task
};
#endif

