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

#ifndef ALIANALYSISTASK_H
#include "AliAnalysisTask.h"
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

class AliHFEpid;
class AliHFEcuts;
class AliHFEextraCuts;

class AliAnalysisTaskDCA : public AliAnalysisTask{
 public:

  typedef enum{
    kPostProcess = 0,
    kImpactPar = 1
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
  
  virtual void ConnectInputData(Option_t *);
  virtual void CreateOutputObjects();
  virtual void Exec(Option_t *);
  virtual void Terminate(Option_t *);


  void PrintStatus() const;
  void Load(TString filename = "impactPar.root");
  void PostProcess();
  void SetHFECuts(AliHFEcuts * const cuts) { fCuts = cuts; };

  void SetPixelStatus(Int_t pixelStatus){ fPixelStatus = pixelStatus;};
  void SetNclustersITS(Int_t nITSclusters){ fNclustersITS = nITSclusters;};
  
  Bool_t GetPlugin(Int_t plug) const { return TESTBIT(fPlugins, plug); };
  void SwitchOnPlugin(Int_t plug);
 
  Bool_t HasMCData() const { return TestBit(kHasMCdata); }
  void SetHasMCData(Bool_t hasMC = kTRUE) { SetBit(kHasMCdata, hasMC); };
  void SetAODAnalysis() { SetBit(kAODanalysis, kTRUE); };
  void SetESDAnalysis() { SetBit(kAODanalysis, kFALSE); };


 private:

  void MakeParticleContainer();
  UShort_t fPlugins;                    // Enabled Plugins                                    
  AliESDEvent *fESD;                    //! The ESD Event
  AliMCEvent *fMC;                      //! The MC Event

  AliHFEcuts *fCuts;                    // Cut Collection
  AliCFManager *fCFM;                   //! Correction Framework Manager
  AliHFEdca *fDCA;                      // fDCA 
  
  
  Int_t fPixelStatus;                      // pixel layer
  Int_t fNclustersITS;                    // ITS clusters
                                               
  TH1I *fNEvents;                       //! counter for the number of Events
  TList *fResidualList;                 //! histograms for the residuals 
  TList *fPullList;                     //! histograms for the pull
  TList *fOutput;                       //! Container for Task Output

  ClassDef(AliAnalysisTaskDCA, 1);      // The DCA Analysis Task
};
#endif

