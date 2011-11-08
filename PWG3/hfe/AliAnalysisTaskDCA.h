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
// Task for impact parameter (DCA) analysis
// study DCA in rphi (xy) and z: resolution and pull 
// For further information see implementation file


#ifndef ALIANALYSISTASKDCA_H
#define ALIANALYSISTASKDCA_H
#ifndef ALIANALYSISTASKSE_H
#include "AliAnalysisTaskSE.h"
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

class AliVEvent;
class AliVertexerTracks;
class AliHFEpid;
class AliHFEcuts;
class AliHFEextraCuts;

class AliAnalysisTaskDCA : public AliAnalysisTaskSE{
 public:

  typedef enum{
    kPostProcess = 0,
    kImpactPar = 1,
    kPrimVtx = 2,
    kCombinedPid = 3,
    kHFEpid = 4, 
    kKFdca = 5
  }Switches_t;

  enum{
    kHasMCdata = BIT(19),
    kAODanalysis = BIT(20)
  };
  
  AliAnalysisTaskDCA();
  AliAnalysisTaskDCA(const char * name);
  AliAnalysisTaskDCA(const AliAnalysisTaskDCA &ref);
  AliAnalysisTaskDCA& operator=(const AliAnalysisTaskDCA &ref);
  virtual ~AliAnalysisTaskDCA();
  
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *);
  virtual void Terminate(Option_t *);


  void PrintStatus() const;
  void Load(TString filename = "impactPar.root");
  void PostProcess();
  void SetHFECuts(AliHFEcuts * const cuts) { fCuts = cuts; };

  void SetNclustersITS(Int_t nITSclusters){ fNclustersITS = nITSclusters;};
  void SetMinPrimVtxContrib( Int_t nPrimVtxContrib){ fMinNprimVtxContrbutor = nPrimVtxContrib; };  

  Bool_t GetPlugin(Int_t plug) const { return TESTBIT(fPlugins, plug); };
  void SwitchOnPlugin(Int_t plug);

  Bool_t IsAODanalysis() const { return TestBit(kAODanalysis); };
  Bool_t IsESDanalysis() const { return !TestBit(kAODanalysis); };
  Bool_t HasMCData() const { return TestBit(kHasMCdata); }
  void SetHasMCData(Bool_t hasMC = kTRUE) { SetBit(kHasMCdata, hasMC); };

  void SetPIDdetectors(Char_t * const detectors){ fPIDdetectors = detectors; }
  void SetPIDStrategy(UInt_t strategy) { fPIDstrategy = strategy; }
  void AddPIDdetector(TString detector);

  void SetAODAnalysis() { SetBit(kAODanalysis, kTRUE); };
  void SetESDAnalysis() { SetBit(kAODanalysis, kFALSE); };


 private:

  void MakeParticleContainer();
  void ProcessDcaAnalysis();
    
  UShort_t fPlugins;                   // Enabled Plugins 
  AliHFEcuts *fCuts;                   // Cut Collection

  AliHFEpid *fHFEpid;                  //! PID
  TString fPIDdetectors;               // Detectors for Particle Identification
  UInt_t fPIDstrategy;                 // PID Strategy

  AliCFManager *fCFM;                  //! Correction Framework Manager
  AliHFEdca *fDCA;                     // fDCA 

  Int_t fNclustersITS;                 // ITS clusters
  Int_t fMinNprimVtxContrbutor;        // minimum number of primary contributors
                                               
  TH1I *fNEvents;                      //! counter for the number of Events
  TList *fResidualList;                //! histograms for the residuals 
  TList *fPullList;                    //! histograms for the pull
  TList *fDcaList;                     //! histograms for the dca
  TList *fKfDcaList;                     //! histograms for the kf dca
  TList *fMcVertexList;                //! histograms for the MC vertex  
  TList *fDataDcaList;                 //! histograms for the data dca
  TList *fDataVertexList;              //! histograms for the data vertex
  TList *fDataPullList;                //! histograms for the data pull
  TList *fMcPidList;                   //! pid - MC: ESD combined pid
  TList *fDataPidList;                 //! pid -Data: ESD combined pid

  TList *fHfeDcaList;                  //! hfe pid: mc dca 
  TList *fHfeDataDcaList;              //! hfe pid: data dca 

  TList *fOutput;                      //! Container for Task Output

  ClassDef(AliAnalysisTaskDCA, 1);     // The DCA Analysis Task
};
#endif

