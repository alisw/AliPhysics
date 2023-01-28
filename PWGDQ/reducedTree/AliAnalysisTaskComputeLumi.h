// Analysis task for creating a reduced tree containing event, track and resonance candidate information 
// Creation date: 2012/06/21 
// Authors: Ionut-Cristian Arsene (i.c.arsene@gsi.de,i.c.arsene@cern.ch)                                 

#ifndef ALIANALYSISTASKCOMPUTELUMI_H
#define ALIANALYSISTASKCOMPUTELUMI_H 

#include <TList.h>
#include <AliAnalysisTaskSE.h>
#include <AliTriggerAnalysis.h>
#include <TH2I.h>

class AliAnalysisCuts;
class TTree;
class TFile;
class TBits;
class AliESDtrack;
class AliAODTrack;
class AliAnalysisUtils;

//_________________________________________________________________________
class AliAnalysisTaskComputeLumi : public AliAnalysisTaskSE {
  
public:
  AliAnalysisTaskComputeLumi();
  AliAnalysisTaskComputeLumi(const char *name);
  virtual ~AliAnalysisTaskComputeLumi(){  }

  virtual void UserExec(Option_t *option);
  virtual void UserCreateOutputObjects();
  virtual void FinishTaskOutput();
  
  // Cuts for selection of event to be written to tree
  void SetEventFilter(AliAnalysisCuts * const filter) {fEventFilter=filter;}
    
 private:

  AliAnalysisCuts *fEventFilter;      // event filter
  TList*           fHistogramList;    // histogram list
  
  AliTriggerAnalysis* fTriggerAnalysis;  // trigger analysis object
  TH2I* fTrigAliasVsCentV0M_before;            // 
  TH2I* fL0TriggerInputsVsCent_before;         //
  TH2I* fL0TriggerInputsVsTrigAlias_before;      //
  TH2I* fVtxVsCentV0_before;                    //
  TH2I* fTrigAliasVsCentV0M_acc;            // 
  TH2I* fL0TriggerInputsVsCent_acc;         //
  TH2I* fL0TriggerInputsVsTrigAlias_acc;      //
  TH2I* fVtxVsCentV0_acc;                    //
  
  AliAnalysisTaskComputeLumi(const AliAnalysisTaskComputeLumi &c);
  AliAnalysisTaskComputeLumi& operator= (const AliAnalysisTaskComputeLumi &c);

  ClassDef(AliAnalysisTaskComputeLumi, 2); //Analysis Task for computing integrated lumi
};
#endif
